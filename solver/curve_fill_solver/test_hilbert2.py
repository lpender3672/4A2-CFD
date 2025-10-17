#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import numba

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import matplotlib.cm as cm

from lod import calc_lod
# ======================================================
# Geometry helpers
# ======================================================


def rotate(points, angle_rad):
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    R = np.array([[c, -s], [s, c]])
    return points @ R.T

def naca4_airfoil(code="2412", n=201, closed_te=True):
    """
    Return a closed polyline (Nx2) for a NACA 4-digit airfoil on chord [0,1], origin at LE.
    Cosine spacing for x. Upper surface 0->1, lower 1->0.
    """
    assert len(code) == 4 and code.isdigit()
    m = int(code[0]) / 100.0            # max camber
    p = int(code[1]) / 10.0             # location of max camber
    t = int(code[2:]) / 100.0           # thickness

    # Cosine-spaced x from 0..1 (cluster near leading edge)
    beta = np.linspace(0.0, np.pi, n)
    x = 0.5 * (1.0 - np.cos(beta))

    # Thickness distribution (closed or open TE)
    k5 = -0.1015 if closed_te else -0.1036
    yt = 5 * t * (0.2969*np.sqrt(np.clip(x, 0, 1)) - 0.1260*x
                  - 0.3516*x**2 + 0.2843*x**3 + k5*x**4)

    # Camber line and slope
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    if m > 0 and p > 0:
        idx1 = x < p
        idx2 = ~idx1
        yc[idx1] = m/(p**2) * (2*p*x[idx1] - x[idx1]**2)
        yc[idx2] = m/((1-p)**2) * ((1 - 2*p) + 2*p*x[idx2] - x[idx2]**2)
        dyc_dx[idx1] = 2*m/(p**2) * (p - x[idx1])
        dyc_dx[idx2] = 2*m/((1-p)**2) * (p - x[idx2])

    theta = np.arctan(dyc_dx)

    xu = x - yt*np.sin(theta)
    yu = yc + yt*np.cos(theta)
    xl = x + yt*np.sin(theta)
    yl = yc - yt*np.cos(theta)

    upper = np.stack([xu, yu], axis=1)
    lower = np.stack([xl[::-1], yl[::-1]], axis=1)  # go back along lower
    poly = np.vstack([upper, lower])
    return poly  # closed path (first != last, but segment wraps naturally)


def transform_airfoil(poly01, chord=1.0, alpha_deg=0.0, origin=(0.0, 0.0)):
    """
    Scale to chord, rotate by alpha (deg), and translate to origin=(x_le, y_le).
    """
    pts = poly01 * chord
    pts = rotate(pts, np.deg2rad(alpha_deg))
    pts[:, 0] += origin[0]
    pts[:, 1] += origin[1]
    return pts

@numba.njit
def distance_point_to_polyline(px, py, poly):
    """
    Compute unsigned Euclidean distance from point (px,py) to a polyline (Nx2).
    Vectorized over segments. Returns scalar distance.
    """
    P = poly[:-1]      # (N-1, 2)
    Q = poly[1:]       # (N-1, 2)
    vx = Q[:, 0] - P[:, 0]
    vy = Q[:, 1] - P[:, 1]

    # Project point onto each segment, clamp to [0,1]
    wx = px - P[:, 0]
    wy = py - P[:, 1]
    vv = vx*vx + vy*vy
    vv = np.where(vv == 0.0, 1e-16, vv)
    t = (wx*vx + wy*vy) / vv
    t = np.clip(t, 0.0, 1.0)

    projx = P[:, 0] + t*vx
    projy = P[:, 1] + t*vy
    dx = px - projx
    dy = py - projy
    d2 = dx*dx + dy*dy
    return float(np.sqrt(np.min(d2)))

def curvature_at_nearest_point(px, py, poly):
    """
    Return (distance, curvature) of the nearest point on a polyline to (px, py).

    poly : (N,2) array of coordinates.
    curvature computed via centered finite differences.

    Curvature κ = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
    """
    poly = np.asarray(poly)
    P = poly[:-1]
    Q = poly[1:]

    vx = Q[:, 0] - P[:, 0]
    vy = Q[:, 1] - P[:, 1]
    wx = px - P[:, 0]
    wy = py - P[:, 1]
    vv = vx * vx + vy * vy
    vv = np.where(vv == 0.0, 1e-16, vv)

    # projection parameter along each segment
    t = (wx * vx + wy * vy) / vv
    t = np.clip(t, 0.0, 1.0)

    projx = P[:, 0] + t * vx
    projy = P[:, 1] + t * vy
    dx = px - projx
    dy = py - projy
    d2 = dx * dx + dy * dy
    i_min = np.argmin(d2)
    dist = np.sqrt(d2[i_min])

    # index of nearest vertex (for curvature estimation)
    # choose the closer endpoint of that segment
    if t[i_min] < 0.5:
        idx = i_min
    else:
        idx = i_min + 1
    idx = np.clip(idx, 1, len(poly) - 2)

    # local curvature (central finite difference)
    x_prev, x_curr, x_next = poly[idx - 1: idx + 2, 0]
    y_prev, y_curr, y_next = poly[idx - 1: idx + 2, 1]

    dx1 = x_next - x_prev
    dy1 = y_next - y_prev
    ddx = x_next - 2 * x_curr + x_prev
    ddy = y_next - 2 * y_curr + y_prev
    denom = (dx1 * dx1 + dy1 * dy1) ** 1.5
    denom = 1e-16 if denom == 0.0 else denom
    curvature = abs(dx1 * ddy - dy1 * ddx) / denom

    return dist, curvature

# ======================================================
# Point-in-polygon (with on-edge = inside)
# ======================================================

@numba.njit
def _poly_close(poly, eps=1e-12):
    # poly: (N, D) float array (e.g., D=2)
    n, d = poly.shape

    # Is first point already (approximately) the last?
    same = True
    for j in range(d):
        if abs(poly[0, j] - poly[n - 1, j]) > eps:
            same = False
            break

    if same:
        return poly  # safe to return as-is

    # Append first row to the end
    out = np.empty((n + 1, d), dtype=poly.dtype)
    for i in range(n):
        for j in range(d):
            out[i, j] = poly[i, j]
    for j in range(d):
        out[n, j] = poly[0, j]
    return out

@numba.njit
def point_in_polygon(x, y, poly, eps=1e-12):
    """
    Ray-casting with on-edge detection. poly is (K,2) closed or open.
    Returns True if inside OR on boundary.
    """
    P = _poly_close(poly)

    # On-edge quick check
    for i in range(len(P) - 1):
        x1, y1 = P[i]
        x2, y2 = P[i+1]
        vx, vy = x2 - x1, y2 - y1
        wx, wy = x - x1, y - y1
        cross = vx*wy - vy*wx
        if abs(cross) <= eps:
            dot = vx*wx + vy*wy
            if -eps <= dot <= (vx*vx + vy*vy) + eps:
                return True

    # Ray cast to +x
    inside = False
    for i in range(len(P) - 1):
        x1, y1 = P[i]
        x2, y2 = P[i+1]
        if ((y1 > y) != (y2 > y)):
            x_int = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x_int >= x - 1e-12:
                inside = not inside
    return inside

# ======================================================
# Polygon/rectangle clipping + area
# ======================================================


def polygon_area(poly):
    """Signed area; positive for CCW. Accepts open or closed poly."""
    P = _poly_close(poly)
    x = P[:,0]; y = P[:,1]
    return 0.5 * np.sum(x[:-1]*y[1:] - x[1:]*y[:-1])


def rect_clip_polygon(poly, x0, y0, x1, y1):
    """
    Sutherland-Hodgman clip of 'poly' against axis-aligned rectangle [x0,x1]x[y0,y1].
    Returns Nx2 array (possibly empty). Accepts open or closed input poly.
    """
    P = _poly_close(poly)

    
    def clip_edge(PA, inside_fn, intersect_fn):
        if len(PA) == 0:
            return PA
        out = []
        for i in range(len(PA)-1):
            A = PA[i]; B = PA[i+1]
            Ain = inside_fn(A); Bin = inside_fn(B)
            if Ain and Bin:
                out.append(B)
            elif Ain and (not Bin):
                out.append(intersect_fn(A, B))
            elif (not Ain) and Bin:
                out.append(intersect_fn(A, B))
                out.append(B)
        if out:
            if not np.allclose(out[0], out[-1]):
                out = [out[0], *out, out[0]]
            return np.array(out)
        return np.empty((0,2))

    # Left: x >= x0
    
    def inside_L(Pt):
        return Pt[0] >= x0
    
    
    def intersect_L(A, B):
        t = (x0 - A[0]) / (B[0] - A[0])
        return np.array([x0, A[1] + t*(B[1] - A[1])])

    # Right: x <= x1
    
    def inside_R(Pt):
        return Pt[0] <= x1
    
    
    def intersect_R(A, B):
        t = (x1 - A[0]) / (B[0] - A[0])
        return np.array([x1, A[1] + t*(B[1] - A[1])])

    # Bottom: y >= y0
    
    def inside_B(Pt):
        return Pt[1] >= y0
    
    
    def intersect_B(A, B):
        t = (y0 - A[1]) / (B[1] - A[1])
        return np.array([A[0] + t*(B[0] - A[0]), y0])

    # Top: y <= y1
    
    def inside_T(Pt):
        return Pt[1] <= y1
    
    
    def intersect_T(A, B):
        t = (y1 - A[1]) / (B[1] - A[1])
        return np.array([A[0] + t*(B[0] - A[0]), y1])

    P = clip_edge(P, inside_L, intersect_L)
    P = clip_edge(P, inside_R, intersect_R)
    P = clip_edge(P, inside_B, intersect_B)
    P = clip_edge(P, inside_T, intersect_T)
    return P

def rect_solid_fraction(poly, x0, y0, x1, y1, poly_bbox=None):
    """
    Fraction of rectangle covered by 'poly' (airfoil).
    Fast checks:
      - bbox reject if provided
      - all 4 rect corners inside -> 1 (airfoil has no holes)
    Otherwise: S-H clip + area / rect area.
    """
    if poly_bbox is not None:
        (bx0, by0, bx1, by1) = poly_bbox
        if (x1 <= bx0) or (bx1 <= x0) or (y1 <= by0) or (by1 <= y0):
            return 0.0

    corners = [(x0,y0),(x1,y0),(x1,y1),(x0,y1)]
    if all(point_in_polygon(cx, cy, poly) for (cx, cy) in corners):
        return 1.0

    clipped = rect_clip_polygon(poly, x0, y0, x1, y1)
    if clipped.size == 0:
        return 0.0
    A_int = abs(polygon_area(clipped))
    A_rect = (x1 - x0) * (y1 - y0)
    if A_rect <= 0:
        return 0.0
    return float(np.clip(A_int / A_rect, 0.0, 1.0))


def build_uniform_bins(poly: np.ndarray, bin_w: float, bin_h: float):
    """
    Bin polyline segments into a uniform grid for fast candidate lookup.
    Returns: (bins, nx, ny, x0, y0, bin_w, bin_h)
      - bins: list of lists; length = nx*ny; bins[k] holds segment indices touching that bin
    """
    P = _poly_close(poly)
    segs = np.c_[P[:-1], P[1:]]  # [x1,y1,x2,y2] per segment
    x0, y0 = np.min(P[:,0]), np.min(P[:,1])
    x1, y1 = np.max(P[:,0]), np.max(P[:,1])
    nx = max(1, int(np.ceil((x1 - x0) / bin_w)))
    ny = max(1, int(np.ceil((y1 - y0) / bin_h)))
    bins = [list() for _ in range(nx*ny)]

    def bin_index(ix, iy): return iy*nx + ix

    for si, (x1s,y1s,x2s,y2s) in enumerate(segs):
        bx0 = min(x1s, x2s); bx1 = max(x1s, x2s)
        by0 = min(y1s, y2s); by1 = max(y1s, y2s)
        ix0 = int(np.clip(np.floor((bx0 - x0)/bin_w), 0, nx-1))
        ix1 = int(np.clip(np.floor((bx1 - x0)/bin_w), 0, nx-1))
        iy0 = int(np.clip(np.floor((by0 - y0)/bin_h), 0, ny-1))
        iy1 = int(np.clip(np.floor((by1 - y0)/bin_h), 0, ny-1))
        for ix in range(ix0, ix1+1):
            for iy in range(iy0, iy1+1):
                bins[bin_index(ix,iy)].append(si)
    return bins, nx, ny, x0, y0, bin_w, bin_h, segs

def rect_solid_fraction_with_candidates(poly, x0, y0, x1, y1, candidates, segs):
    """
    Same as rect_solid_fraction but only considers 'candidates' segments.
    We still run Sutherland–Hodgman (needs the full poly order),
    yet the quick tests use candidates for the early-exit decisions.
    """
    # Quick bbox reject using candidate list
    if len(candidates) == 0:
        return 0.0

    # Corner-in-poly check (cheap and decisive for φ=1)
    corners = [(x0,y0),(x1,y0),(x1,y1),(x0,y1)]
    if all(point_in_polygon(cx, cy, poly) for (cx, cy) in corners):
        return 1.0

    # Clip and area
    clipped = rect_clip_polygon(poly, x0, y0, x1, y1)
    if clipped.size == 0:
        return 0.0
    A_int = abs(polygon_area(clipped))
    A_rect = (x1 - x0) * (y1 - y0)
    return float(np.clip(A_int / A_rect, 0.0, 1.0)) if A_rect > 0 else 0.0


def compute_phi_batched(cells, poly, bin_w=None, bin_h=None):
    """
    Batched φ for all cells using a uniform bin accelerator.
    cells: list of (py, px, side[, ...]) in any order (Hilbert order preserved)
    poly : (K,2) airfoil polyline (closed or open)
    bin_w/bin_h: choose ~ typical cell size near interface (or auto)
    Returns: phis (np.ndarray) aligned with 'cells'
    """
    # Choose bin sizes: default to median cell side for decent balance
    sides = np.array([c[2] for c in cells], dtype=float)
    s_med = np.median(sides) if sides.size else 1.0
    if bin_w is None: bin_w = s_med
    if bin_h is None: bin_h = s_med

    bins, nx, ny, gx0, gy0, bw, bh, segs = build_uniform_bins(poly, bin_w, bin_h)

    def bins_for_rect(x0, y0, x1, y1):
        ix0 = int(np.clip(np.floor((x0 - gx0)/bw), 0, nx-1))
        ix1 = int(np.clip(np.floor((x1 - gx0)/bw), 0, nx-1))
        iy0 = int(np.clip(np.floor((y0 - gy0)/bh), 0, ny-1))
        iy1 = int(np.clip(np.floor((y1 - gy0)/bh), 0, ny-1))
        return ix0, ix1, iy0, iy1

    phis = np.empty(len(cells), dtype=float)
    for i, (py, px, side, *_) in enumerate(cells):
        x0 = px - 0.5*side; x1 = px + 0.5*side
        y0 = py - 0.5*side; y1 = py + 0.5*side
        ix0, ix1, iy0, iy1 = bins_for_rect(x0, y0, x1, y1)

        # Union of candidate segments from overlapping bins
        cand = set()
        for ix in range(ix0, ix1+1):
            for iy in range(iy0, iy1+1):
                cand.update(bins[iy*nx + ix])

        # Early reject: if no nearby segments AND any corner outside, φ likely 0.
        # We still allow corner-in-poly to detect φ=1.
        if not cand:
            # quick φ=1 check then φ=0 fallback
            if all(point_in_polygon(cx, cy, poly) for (cx, cy) in [(x0,y0),(x1,y0),(x1,y1),(x0,y1)]):
                phis[i] = 1.0
            else:
                phis[i] = 0.0
            continue

        phis[i] = rect_solid_fraction_with_candidates(poly, x0, y0, x1, y1, cand, segs)

    return phis

# ---------- Segment geometry ----------

def _seg_nearest_t(p, a, b):
    ab = b - a
    den = np.sum(ab*ab, axis=-1, keepdims=True)
    den = np.where(den == 0.0, 1.0, den)     # avoid 0/0
    t = np.sum((p - a)*ab, axis=-1, keepdims=True) / den
    return np.clip(t, 0.0, 1.0)


def _seg_dist2_to_point(p, a, b):
    """Squared distance from point p to segment a->b."""
    t = _seg_nearest_t(p, a, b)
    q = a + t*(b - a)
    d = p - q
    return d.dot(d), t, q


def _signed_distance_to_line(p, a, b):
    """
    Signed distance from point p to the infinite line through a->b.
    Positive if p is to the 'left' of the directed edge a->b.
    """
    v = b - a
    n = np.array([-v[1], v[0]])  # left normal
    ln = np.linalg.norm(v)
    if ln == 0.0:
        return 0.0
    return np.dot(p - a, n) / ln


def _signed_distance_to_line_pts(P, A, B):
    """
    Vectorized signed distance of points P (...,2) to line A->B (...,2).
    Positive on the 'left' of A->B.
    """
    V = B - A
    N = np.stack([-V[...,1], V[...,0]], axis=-1)
    ln = np.linalg.norm(V, axis=-1, keepdims=True)
    ln = np.where(ln == 0.0, 1.0, ln)
    return np.sum((P - A) * N, axis=-1, keepdims=True) / ln

# ---------- Half-plane clip (rectangle ∩ {sd*s_keep >= 0}) ----------
def _rect_halfplane_fraction(x0, y0, x1, y1, a, b, keep_sign):
    """
    Single-cell: area fraction of [x0,x1]×[y0,y1] in half-plane of line a->b.
    keep_sign ∈ {+1,-1}: keep points with keep_sign * signed_distance >= 0.
    """
    rect = np.array([[x0,y0],[x1,y0],[x1,y1],[x0,y1]], dtype=float)
    # signed distances (column vector)
    sd = _signed_distance_to_line_pts(rect[None, ...], a[None, ...], b[None, ...]).ravel()
    sd *= keep_sign
    inside = sd >= 0.0
    if inside.all():   return 1.0
    if (~inside).all(): return 0.0

    # Sutherland–Hodgman vs single half-plane (4 vertices → up to 6)
    def interp(p1, p2, s1, s2):
        t = s1 / (s1 - s2)
        return p1 + t*(p2 - p1)

    poly = []
    for i in range(4):
        j = (i + 1) % 4
        p_i, p_j = rect[i], rect[j]
        s_i, s_j = sd[i], sd[j]
        in_i, in_j = inside[i], inside[j]

        if in_i and in_j:
            poly.append(p_j)
        elif in_i and not in_j:
            poly.append(interp(p_i, p_j, s_i, s_j))
        elif (not in_i) and in_j:
            poly.append(interp(p_i, p_j, s_i, s_j))
            poly.append(p_j)

    if len(poly) < 3:
        return 0.0
    P = np.array(poly)
    x, y = P[:,0], P[:,1]
    area_clip = 0.5*np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    area_rect = (x1 - x0) * (y1 - y0)
    return float(np.clip(area_clip / area_rect, 0.0, 1.0))

# ---------- Fast φ using linear wall per cell ----------

def compute_phi_batched_linear_approx(
    cells, poly, bins_tuple=None, bin_w=None, bin_h=None,
    far_factor=0.75
):
    """
    Compute φ for all cells using a linearized wall (nearest segment) per cell.
    - If 'bins_tuple' is None, builds uniform bins (auto bin size ~ median side).
    - far_factor: if dist_center >= far_factor * (cell_diag/2), snap φ to 0/1 via center-in-poly.

    Returns: np.ndarray φ aligned with 'cells'.
    """
    # Prepare bins if needed
    if bins_tuple is None:
        bins_tuple = build_uniform_bins(poly, 
                                        np.median([c[2] for c in cells]) if cells else 1.0,
                                        np.median([c[2] for c in cells]) if cells else 1.0)
    bins, nx, ny, gx0, gy0, bw, bh, segs = bins_tuple
    P = _poly_close(poly)

    def bins_for_rect(x0, y0, x1, y1):
        ix0 = int(np.clip(np.floor((x0 - gx0)/bw), 0, nx-1))
        ix1 = int(np.clip(np.floor((x1 - gx0)/bw), 0, nx-1))
        iy0 = int(np.clip(np.floor((y0 - gy0)/bh), 0, ny-1))
        iy1 = int(np.clip(np.floor((y1 - gy0)/bh), 0, ny-1))
        return ix0, ix1, iy0, iy1

    phis = np.empty(len(cells), dtype=float)
    far2_cache = []  # optional, for analysis

    for i, (py, px, side, *_) in enumerate(cells):
        x0 = px - 0.5*side; x1 = px + 0.5*side
        y0 = py - 0.5*side; y1 = py + 0.5*side
        c = np.array([px, py])
        half_diag = 0.5*np.sqrt(2.0)*side  # radius of circumscribed circle

        # Candidate segments from overlapping bins
        ix0, ix1, iy0, iy1 = bins_for_rect(x0, y0, x1, y1)
        cand = set()
        for ix in range(ix0, ix1+1):
            for iy in range(iy0, iy1+1):
                cand.update(bins[iy*nx + ix])

        if not cand:
            # No nearby segments → snap to 0/1 via center-in-poly
            phis[i] = 1.0 if point_in_polygon(px, py, poly) else 0.0
            continue

        # Find nearest segment (by center distance)
        best_d2 = np.inf
        best_a = None
        best_b = None
        for si in cand:
            x1s,y1s,x2s,y2s = segs[si]
            a = np.array([x1s, y1s]); b = np.array([x2s, y2s])
            d2, t, q = _seg_dist2_to_point(c, a, b)
            if d2 < best_d2:
                best_d2 = d2; best_a = a; best_b = b

        # Far test: if the line is far relative to the cell size, snap to 0/1
        if best_d2 >= (far_factor*half_diag)**2:
            phis[i] = 1.0 if point_in_polygon(px, py, poly) else 0.0
            continue

        # Decide which side of the line is 'solid' using the center point
        sd_center = _signed_distance_to_line(c, best_a, best_b)
        keep_sign = +1 if point_in_polygon(px, py, poly) else -1
        # If center is inside, we keep the same sign as center; else keep the opposite
        # (Equivalently: keep points with sign(sd)*keep_sign >= 0)
        phi = _rect_halfplane_fraction(x0, y0, x1, y1, best_a, best_b, keep_sign=np.sign(sd_center) if keep_sign>0 else -np.sign(sd_center))
        # Numerical jitter guard
        phis[i] = 0.0 if phi < 1e-12 else 1.0 if phi > 1.0 - 1e-12 else phi

    return phis

def compute_phi_batched_linear_numpy(
    cells,
    poly,
    bin_w=None,
    bin_h=None,
    neighborhood=1,     # 0=center bin only, 1 = 3×3 around center (recommended)
    kmax=128,           # cap on candidate segments per cell (after 3×3 union)
    chunk=4096,         # process cells in chunks to keep memory in check
    far_factor=0.75     # snap φ to 0/1 if line is far vs cell half-diagonal
):
    """
    Vectorization strategy:
      - Vectorize nearest-segment search across each chunk via broadcasting.
      - Use uniform bins to build small candidate sets per cell (Python-level gather,
        but cheap), then compute distances in one big NumPy op per chunk.
      - Final half-plane area per cell is a tiny O(1) function.

    Returns: np.ndarray of φ aligned with 'cells'.
    """
    P = _poly_close(poly)
    segs = np.c_[P[:-1], P[1:]]  # (Ns, 4) → [x1,y1,x2,y2]
    Ns = segs.shape[0]

    # Choose bin sizes from median cell side if not given
    sides = np.array([c[2] for c in cells], dtype=float)
    s_med = np.median(sides) if sides.size else 1.0
    if bin_w is None: bin_w = s_med
    if bin_h is None: bin_h = s_med

    # Build bins (segment -> bins mapping). Ns ~ O(1e3), fine to do in Python.
    bins, nx, ny, gx0, gy0, bw, bh, segs = build_uniform_bins(poly, bin_w, bin_h)

    # Precompute per-bin CSR structure for fast gather (optional)
    # Flat list per bin is fine; we’ll gather with Python but keep arrays small.

    # Cells arrays
    C = np.array([(c[1], c[0]) for c in cells], dtype=float)   # centers as (x,y)
    S = sides.copy()
    n_cells = len(cells)

    # Bin ids for centers
    ix_c = np.clip(np.floor((C[:,0] - gx0)/bw).astype(int), 0, nx-1)
    iy_c = np.clip(np.floor((C[:,1] - gy0)/bh).astype(int), 0, ny-1)

    # Neighbor offsets (vectorized 3×3 or (2n+1)^2 neighborhood)
    offs = np.arange(-neighborhood, neighborhood+1, dtype=int)
    dxy = np.stack(np.meshgrid(offs, offs, indexing='ij'), axis=-1).reshape(-1,2)  # (Nb,2)
    Nb = dxy.shape[0]

    phis = np.empty(n_cells, dtype=float)

    # Process in chunks
    for start in range(0, n_cells, chunk):
        end = min(start + chunk, n_cells)
        idx = slice(start, end)

        Cb = C[idx]         # (B,2)
        Sb = S[idx]         # (B,)
        ixb = ix_c[idx]; iyb = iy_c[idx]

        # Build candidate segment sets for each cell (gather from neighborhood bins)
        cand_lists = []
        for i in range(Cb.shape[0]):
            seg_set = set()
            ix0 = ixb[i];  iy0 = iyb[i]
            for dx, dy in dxy:
                ix = int(np.clip(ix0 + dx, 0, nx-1))
                iy = int(np.clip(iy0 + dy, 0, ny-1))
                seg_set.update(bins[iy*nx + ix])
                if len(seg_set) >= kmax:
                    break
            # cap and store
            cl = np.fromiter(seg_set, dtype=int, count=min(kmax, len(seg_set)))
            cand_lists.append(cl)

        # Pad candidate arrays to (B, K) with -1; keep lengths
        K = max((len(cl) for cl in cand_lists), default=0)
        K = min(K, kmax)
        if K == 0:
            # No candidates for all cells in this chunk → snap 0/1 by center
            inside_mask = np.array([point_in_polygon(Cb[i,0], Cb[i,1], poly)
                                    for i in range(Cb.shape[0])], dtype=bool)
            phis[start:end] = np.where(inside_mask, 1.0, 0.0)
            continue

        cand_idx = -np.ones((Cb.shape[0], K), dtype=int)
        for i, cl in enumerate(cand_lists):
            L = min(len(cl), K)
            if L: cand_idx[i, :L] = cl[:L]

        valid_mask = cand_idx >= 0

        # Gather segment endpoints for all candidates in the chunk
        # segs: (Ns,4) → A=(x1,y1), B=(x2,y2)
        A = np.zeros((Cb.shape[0], K, 2), dtype=float)
        B = np.zeros((Cb.shape[0], K, 2), dtype=float)
        sel = cand_idx[valid_mask]
        if sel.size:
            AB = segs[sel]   # (Ncand, 4)
            A[valid_mask] = AB[:, :2]
            B[valid_mask] = AB[:, 2:]

        # Compute squared distance from each center to each candidate segment (vectorized)
        # Project with clamped t
        Pcent = Cb[:, None, :]                # (B, 1, 2)
        AB = B - A                            # (B, K, 2)
        den = np.sum(AB*AB, axis=-1, keepdims=True)
        den = np.where(den == 0.0, 1.0, den)
        t = np.sum((Pcent - A) * AB, axis=-1, keepdims=True) / den
        t = np.clip(t, 0.0, 1.0)
        Q = A + t * AB                        # nearest points (B, K, 2)
        D2 = np.sum((Pcent - Q)**2, axis=-1)  # (B, K)

        # Mask out invalid candidates by setting distance large
        D2 = np.where(valid_mask, D2, np.inf)

        # Pick nearest segment per cell
        jmin = np.argmin(D2, axis=1)          # (B,)
        a_near = A[np.arange(Cb.shape[0]), jmin]  # (B,2)
        b_near = B[np.arange(Cb.shape[0]), jmin]  # (B,2)

        # Far snap: if nearest line far relative to cell size → 0/1 by center only
        half_diag = 0.5*np.sqrt(2.0) * Sb
        dmin = np.sqrt(np.min(D2, axis=1))
        snap = dmin >= (far_factor * half_diag)

        # Compute φ for each cell in this chunk
        out = np.empty(Cb.shape[0], dtype=float)
        # First: snap 0/1 quickly
        if np.any(snap):
            inside_snap = np.array([point_in_polygon(Cb[i,0], Cb[i,1], poly)
                                    for i in np.where(snap)[0]], dtype=bool)
            out[snap] = np.where(inside_snap, 1.0, 0.0)

        # Remaining: half-plane area with nearest segment
        idx_hp = np.where(~snap)[0]
        for ii in idx_hp:
            px, py = Cb[ii]
            side = Sb[ii]
            x0, x1 = px - 0.5*side, px + 0.5*side
            y0, y1 = py - 0.5*side, py + 0.5*side
            a = a_near[ii]; b = b_near[ii]
            # decide which side to keep: inside half-plane should include the center if center is inside
            center_inside = point_in_polygon(px, py, poly)
            sd_center = float(_signed_distance_to_line_pts(np.array([[px,py]]), a[None,:], b[None,:]))
            keep_sign = +1 if (center_inside and sd_center >= 0) or ((not center_inside) and sd_center < 0) else -1
            out[ii] = _rect_halfplane_fraction(x0, y0, x1, y1, a, b, keep_sign)

        phis[start:end] = np.clip(out, 0.0, 1.0)

    return phis

def polygon_to_edges(polygon):
    """Return array of edges [(x0, y0, x1, y1)] from polygon vertices."""
    poly = np.asarray(polygon)
    if not np.allclose(poly[0], poly[-1]):
        poly = np.vstack([poly, poly[0]])
    return np.hstack([poly[:-1], poly[1:]])

def fast_shading_fraction(cells, polygon):
    edges = polygon_to_edges(polygon)
    results = []

    for py, px, side in cells:
        # Cell corners (clockwise)
        half = side / 2.0
        corners = np.array([
            [px - half, py - half],
            [px + half, py - half],
            [px + half, py + half],
            [px - half, py + half]
        ])

        # Signed distance of corners to each edge (positive outside)
        inside_mask = np.zeros(len(corners), dtype=bool)
        for x0, y0, x1, y1 in edges:
            # Edge normal pointing outward (assuming CCW polygon)
            ex, ey = x1 - x0, y1 - y0
            nx, ny = ey, -ex
            d = (corners[:, 0] - x0) * nx + (corners[:, 1] - y0) * ny
            if np.all(d <= 0):  # All corners inside this edge
                continue
            inside_mask |= (d <= 0)

        inside_frac = inside_mask.sum() / 4.0

        # --- Wall normal estimation (only for partially shaded cells)
        nx = ny = 0.0
        if 0 < inside_frac < 1:
            intersections = []
            for (x0, y0, x1, y1) in edges:
                for i in range(4):
                    xA, yA = corners[i]
                    xB, yB = corners[(i + 1) % 4]
                    denom = (x1 - x0) * (yB - yA) - (y1 - y0) * (xB - xA)
                    if abs(denom) < 1e-12:
                        continue
                    t = ((xA - x0) * (yB - yA) - (yA - y0) * (xB - xA)) / denom
                    u = ((xA - x0) * (y1 - y0) - (yA - y0) * (x1 - x0)) / denom
                    if 0 <= t <= 1 and 0 <= u <= 1:
                        ix = x0 + t * (x1 - x0)
                        iy = y0 + t * (y1 - y0)
                        intersections.append((ix, iy))
            if len(intersections) >= 2:
                (xA, yA), (xB, yB) = intersections[:2]
                nx, ny = yB - yA, -(xB - xA)
                nlen = np.hypot(nx, ny)
                if nlen > 0:
                    nx /= nlen
                    ny /= nlen

        results.append((inside_frac, nx, ny))
    return np.array(results)

# ======================================================
# Hilbert traversal + LOD
# ======================================================

@numba.njit
def rect_intersects_rect(ax0, ay0, ax1, ay1, bx0, by0, bx1, by1):
    return not (ax1 <= bx0 or bx1 <= ax0 or ay1 <= by0 or by1 <= ay0)


def generate_hilbert_airfoil_lod(n, m, airfoil_poly, band_spec, extra_global_levels=1):
    """
    n, m: domain height/width (grid units)
    airfoil_poly: (K,2) array in same coordinate system as domain
    band_spec: list[(radius, stop_level)] in domain units, sorted by radius ascending.
               Smaller stop_level → finer detail.
    extra_global_levels: how many extra Hilbert levels above the base power-of-two pad.
    Returns:
      cells: list of (py, px, side, phi) in Hilbert order
      fine_bits: int
    """
    base_bits = int(np.ceil(np.log2(max(n, m))))
    fine_bits = base_bits + int(extra_global_levels)
    N = float(1 << fine_bits)

    # Precompute polygon bbox for quick rejects
    af = _poly_close(airfoil_poly)
    bx0, by0 = np.min(af[:,0]), np.min(af[:,1])
    bx1, by1 = np.max(af[:,0]), np.max(af[:,1])
    poly_bbox = (bx0, by0, bx1, by1)

    max_level = 5

    out_cells = []

    # calc lod contours for interpolation
    X, Y, LOD = calc_lod(n, m, airfoil_poly, max_level)

    def lod_from_xy(x, y):
        """Simple integer lookup of LOD at (x,y) in domain units."""

        ix = int(np.clip(np.floor(x * (X.shape[1]-1) / m), 0, X.shape[1]-1))
        iy = int(np.clip(np.floor(y * (Y.shape[0]-1) / n), 0, Y.shape[0]-1))
        return int(LOD[iy, ix])

    @numba.njit
    def stop_level_from_distance(dist, bands):
        """
        bands: list of (radius, stop_level), sorted by radius.
        Returns the stop_level for the first band with dist ≤ radius.
        """
        for r, lvl in bands:
            if dist <= r:
                return int(lvl)
        return int(bands[-1][1])
    
    def stop_level_from_curvature(kappa, dist, bands):
        """
        bands: list of (threshold, level)
        Higher curvature → deeper refinement
        """
        alpha = 200
        beta = 0.7

        M = (kappa ** beta) / (1.0 + alpha * dist) * 1e7

        # Assign level based on descending thresholds
        for thresh, lvl in bands:
            if M <= thresh:
                return int(lvl)
        return int(bands[-1][1])
    
    domain_w=float(m)
    domain_h=float(n)

    def hilbert_adaptive_airfoil(x0, y0, xi, xj, yi, yj, level):
        """
        Recursive Hilbert traversal that emits leaf cells as (py, px, side, phi),
        where phi is the solid (airfoil) area fraction within the cell.
        """
        # Region AABB (fine-grid coords)
        rx0 = min(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        rx1 = max(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        ry0 = min(y0, y0 + xj, y0 + yj, y0 + xj + yj)
        ry1 = max(y0, y0 + xj, y0 + yj, y0 + xj + yj)

        # Completely outside domain? skip
        if not rect_intersects_rect(rx0, ry0, rx1, ry1, 0, 0, domain_w, domain_h):
            return

        # Region center
        px = x0 + 0.5*(xi + yi)
        py = y0 + 0.5*(xj + yj)

        stop_level = lod_from_xy(px, py)

        if level <= stop_level:
            # Leaf cell → compute solid fraction over [rx0,rx1]×[ry0,ry1]

            side = max(rx1 - rx0, ry1 - ry0)  # == 2**stop_level in fine-grid units
            if 0 <= px < domain_w and 0 <= py < domain_h:
                out_cells.append((py, px, side))
            return

        # Recurse in Hilbert order
        hilbert_adaptive_airfoil(x0,                     y0,                      yi/2,  yj/2,  xi/2,  xj/2,  level-1)
        hilbert_adaptive_airfoil(x0 + xi/2,             y0 + xj/2,               xi/2,  xj/2,  yi/2,  yj/2,  level-1)
        hilbert_adaptive_airfoil(x0 + xi/2 + yi/2,      y0 + xj/2 + yj/2,        xi/2,  xj/2,  yi/2,  yj/2,  level-1)
        hilbert_adaptive_airfoil(x0 + xi/2 + yi,        y0 + xj/2 + yj,         -yi/2, -yj/2, -xi/2, -xj/2, level-1)

    hilbert_adaptive_airfoil(
        x0=0.0, y0=0.0,
        xi=N, xj=0.0,
        yi=0.0, yj=N,
        level=fine_bits,
    )
    return out_cells, fine_bits

def bands_from_chord_fractions(chord, pairs):
    """pairs: [(frac, stop_level), ...] sorted asc."""
    return [(chord*frac, lvl) for (frac, lvl) in pairs]

# ======================================================
# Postprocessing helpers
# ======================================================

def cut_segments_by_phi(cells, keep_pred):
    """
    Split Hilbert polyline into contiguous segments based on a predicate on phi.
    keep_pred(phi) -> bool
    Returns (kept_segments, removed_segments), each a list of arrays of shape (L,2)
    """
    kept, removed = [], []
    curr_k, curr_r = [], []
    keeping = None
    for (py, px, side, phi) in cells:
        pt = (px, py)
        km = bool(keep_pred(phi))
        if keeping is None or km != keeping:
            if keeping is True and curr_k: kept.append(np.array(curr_k)); curr_k=[]
            if keeping is False and curr_r: removed.append(np.array(curr_r)); curr_r=[]
            keeping = km
        if km: curr_k.append(pt)
        else:  curr_r.append(pt)
    if keeping is True and curr_k: kept.append(np.array(curr_k))
    if keeping is False and curr_r: removed.append(np.array(curr_r))
    return kept, removed

# ======================================================
# Plotting demo
# ======================================================

def plot_hilbert_with_airfoil_lod(n=256, m=128,
                                  naca_code="2412",
                                  chord=90.0,
                                  alpha_deg=5.0,
                                  origin=(20.0, 64.0),
                                  band_fracs=((0.01, 1), (0.03, 2), (0.08, 3), (1e9, 4)),
                                  airfoil_pts=401,
                                  extra_global_levels=1):
    """
    band_fracs: list[(fraction_of_chord, stop_level)] sorted asc;
                last one may use a very large number for "everything else".
    """
    rotated = False
    if m > n:
        rotated = True
        n, m = m, n
        alpha_deg += 90
        origin = (origin[1], origin[0])

    # Build airfoil in domain coords (x right, y up — plotting later inverts y)
    poly01 = naca4_airfoil(naca_code, n=airfoil_pts, closed_te=True)
    foil = transform_airfoil(poly01, chord=chord, alpha_deg=alpha_deg, origin=origin)
    if not np.allclose(foil[0], foil[-1]):
        foil = np.vstack([foil, foil[0]])

    # Bands from chord fractions → absolute units
    bands = bands_from_chord_fractions(chord, band_fracs)

    # Generate Hilbert cells with solid fraction
    cells, fine_bits = generate_hilbert_airfoil_lod(n, m, foil, bands, extra_global_levels=extra_global_levels)

    print("hilbert traversed with", len(cells), "cells at fine_bits =", fine_bits)

    bins_tuple = build_uniform_bins(foil, 
                                bin_w=np.median([c[2] for c in cells]),
                                bin_h=np.median([c[2] for c in cells]))

    phis = compute_phi_batched_linear_numpy(cells, foil,
                                                bin_w=np.median([c[2] for c in cells]),
                                                bin_h=np.median([c[2] for c in cells]))

    #ans = fast_shading_fraction(cells, foil)
    #phis = ans[:,0]
        
    print("phis computed")
    
    cells_with_phi = [(py, px, side, phi) for (phi, (py, px, side, *rest)) in zip(phis, cells)]

    # Option A: cut out fully-solid segments (phi==1.0)
    kept, removed = cut_segments_by_phi(cells_with_phi, keep_pred=lambda phi: phi < 1.0)


    # Prepare figure
    fig, ax = plt.subplots(figsize=(9, 6))

    def maybe_swap_seg(seg):
        if not rotated: return seg
        if seg.size == 0: return seg
        return np.stack([seg[:,1], seg[:,0]], axis=1)

    # Draw kept (outside or boundary-cut) segments
    combined = np.vstack(kept) if kept else np.empty((0,2))
    combined = maybe_swap_seg(combined)

    for seg in kept:
        segp = maybe_swap_seg(seg)
        #ax.plot(segp[:, 0], segp[:, 1], '-', lw=1)

    ax.plot(combined[:, 0], combined[:, 1], '-', lw=1, color='C0', label='exterior + boundary segments')

    # Optional: show removed interior segments faintly for QA
    # for seg in removed:
    #     segp = maybe_swap_seg(seg)
    #     ax.plot(segp[:, 0], segp[:, 1], lw=1, alpha=0.15)

    # Draw airfoil
    foil_plot = foil[:, ::-1] if rotated else foil
    ax.plot(foil_plot[:, 0], foil_plot[:, 1], lw=2)

    # Limits/aspect
    if rotated:
        ax.set_xlim(0, n); ax.set_ylim(0, m)
    else:
        ax.set_xlim(0, m); ax.set_ylim(0, n)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(f"Hilbert LOD around NACA {naca_code} — interior segments removed")

    # Annotate band spec
    txt = "\n".join([f"≤ {r:.2f} → stop={lvl}" if r < 1e8 else f"> last → stop={lvl}"
                     for (r, lvl) in bands])
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va='top', ha='left',
            bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="0.5", alpha=0.9), fontsize=9)
    
    
    shade_cells_by_phi(
        cells_with_phi,
        n=n, m=m,
        ax=ax,
        rotated=rotated,
        filter_pred=None,          # or: lambda phi: 0.0 < phi < 1.0  (just boundary cells)
        cmap_name="Oranges",
        vmin=0.0, vmax=1.0,
        alpha=0.8,
        draw_edges=True
    )

    plt.show()

# helper plotting func

def shade_cells_by_phi(
    cells,
    n, m,
    ax=None,
    rotated=False,
    filter_pred=None,        # e.g., lambda phi: 0.0 < phi < 1.0
    cmap_name="viridis",     # choose any Matplotlib colormap
    vmin=0.0, vmax=1.0,
    alpha=1.0,
    draw_edges=False,
    edgecolor="black",
    linewidth=0.5,
    add_colorbar=True,
    colorbar_label="solid fraction φ"
):
    """
    Render leaf cells (py, px, side, phi) as shaded squares colored by phi.

    cells: list[(py, px, side, phi)] in Hilbert order (output of generate_hilbert_airfoil_lod)
    n, m: domain (height, width) in the same coordinate system as cells
    rotated: if True, swap x<->y for display (to match your rotated plotting path)
    filter_pred: optional predicate on phi; if provided, only draw matching cells
    cmap_name: Matplotlib colormap name
    vmin, vmax: normalization range for phi
    alpha: face transparency
    draw_edges: whether to draw rectangle edges
    edgecolor, linewidth: edge styling (ignored if draw_edges=False)
    add_colorbar: attach a colorbar to the same axes
    colorbar_label: label for the colorbar

    Returns (pcoll, cbar) where cbar may be None if add_colorbar=False.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # Build rectangles + phi vector (optionally filtered)
    patches = []
    phis = []

    for (py, px, side, phi) in cells:
        if (filter_pred is not None) and (not filter_pred(phi)):
            continue

        # bottom-left of the cell in "fine-grid" display units
        x0 = px - 0.5 * side
        y0 = py - 0.5 * side
        w = side
        h = side

        if rotated:
            # swap axes and dims for display
            rect = Rectangle((y0, x0), h, w)  # note swapped origin and size
        else:
            rect = Rectangle((x0, y0), w, h)

        patches.append(rect)
        phis.append(phi)

    if not patches:
        # Nothing to draw
        return None, None

    # Normalize + colormap
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    cmap = cm.get_cmap(cmap_name)

    # Patch collection for efficient draw
    pcoll = PatchCollection(
        patches,
        cmap=cmap,
        norm=norm,
        match_original=False
    )
    pcoll.set_array(np.asarray(phis, dtype=float))
    pcoll.set_alpha(alpha)
    if draw_edges:
        pcoll.set_edgecolor(edgecolor)
        pcoll.set_linewidth(linewidth)
        pcoll.set_linestyle('solid')
    else:
        pcoll.set_edgecolor("none")

    ax.add_collection(pcoll)

    # Axes limits & aspect to match your convention
    if rotated:
        ax.set_xlim(0, n)
        ax.set_ylim(0, m)
    else:
        ax.set_xlim(0, m)
        ax.set_ylim(0, n)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    cbar = None
    if add_colorbar:
        cbar = plt.colorbar(pcoll, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(colorbar_label)

    return pcoll, cbar

if __name__ == "__main__":
    # Domain (height n, width m)
    n, m = 512, 1024

    # Airfoil placement: chord 45% of width, centered in height, slight incidence
    chord = 0.45 * m
    origin = (m*0.15, n*0.55)  # leading edge position
    alpha_deg = 6.0

    # LOD bands as fractions of chord (closer gets smaller stop_level = finer)
    band_fracs = [
        (0.06, 4),    # ultra-fine skin
        (0.15, 3),
        (0.30, 2),
        (1e9,  1),    # everything else
    ]

    plot_hilbert_with_airfoil_lod(
        n=n, m=m,
        naca_code="2412",
        chord=chord,
        alpha_deg=alpha_deg,
        origin=origin,
        band_fracs=band_fracs,
        airfoil_pts=801,
        extra_global_levels=1,
    )
