import numpy as np
import matplotlib.pyplot as plt


# ---------------------------
# Geometry helpers
# ---------------------------

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
    # Closed TE constant: -0.1015; open TE: -0.1036
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
    # Avoid divide-by-zero for degenerate segments
    vv = np.where(vv == 0.0, 1e-16, vv)
    t = (wx*vx + wy*vy) / vv
    t = np.clip(t, 0.0, 1.0)

    projx = P[:, 0] + t*vx
    projy = P[:, 1] + t*vy
    dx = px - projx
    dy = py - projy
    d2 = dx*dx + dy*dy
    return float(np.sqrt(np.min(d2)))


# ---------------------------
# Hilbert traversal (same structure as before)
# ---------------------------

def rect_intersects_rect(ax0, ay0, ax1, ay1, bx0, by0, bx1, by1):
    return not (ax1 <= bx0 or bx1 <= ax0 or ay1 <= by0 or by1 <= ay0)

def stop_level_from_distance(dist, bands):
    """
    bands: list of (radius, stop_level), sorted by radius.
    Returns the stop_level for the first band with dist ≤ radius.
    """
    for r, lvl in bands:
        if dist <= r:
            return int(lvl)
    return int(bands[-1][1])

def hilbert_adaptive_airfoil(x0, y0, xi, xj, yi, yj, level,
                             domain_w, domain_h,
                             bands, airfoil_poly,
                             points):
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

    # Distance to airfoil polyline in domain units
    dist = distance_point_to_polyline(px, py, airfoil_poly)

    # Map distance to stop level
    stop_level = stop_level_from_distance(dist, bands)

    if level <= stop_level:
        if 0 <= px < domain_w and 0 <= py < domain_h:
            points.append((py, px))
        return

    # Recurse in Hilbert order
    hilbert_adaptive_airfoil(x0,                     y0,                      yi/2,  yj/2,  xi/2,  xj/2,  level-1,
                             domain_w, domain_h, bands, airfoil_poly, points)
    hilbert_adaptive_airfoil(x0 + xi/2,             y0 + xj/2,               xi/2,  xj/2,  yi/2,  yj/2,  level-1,
                             domain_w, domain_h, bands, airfoil_poly, points)
    hilbert_adaptive_airfoil(x0 + xi/2 + yi/2,      y0 + xj/2 + yj/2,        xi/2,  xj/2,  yi/2,  yj/2,  level-1,
                             domain_w, domain_h, bands, airfoil_poly, points)
    hilbert_adaptive_airfoil(x0 + xi/2 + yi,        y0 + xj/2 + yj,         -yi/2, -yj/2, -xi/2, -xj/2, level-1,
                             domain_w, domain_h, bands, airfoil_poly, points)

def generate_hilbert_airfoil_lod(n, m, airfoil_poly, band_spec, extra_global_levels=1):
    """
    n, m: domain height/width (grid units)
    airfoil_poly: (K,2) array in same coordinate system as domain (x right, y down or up; we only use distances)
    band_spec: list[(radius, stop_level)] in domain units, sorted by radius ascending.
               Smaller stop_level → finer detail.
    extra_global_levels: how many extra Hilbert levels above the base power-of-two pad to enable fine cells.
    """
    base_bits = int(np.ceil(np.log2(max(n, m))))
    fine_bits = base_bits + int(extra_global_levels)
    N = float(1 << fine_bits)

    points = []
    hilbert_adaptive_airfoil(
        x0=0.0, y0=0.0,
        xi=N, xj=0.0,
        yi=0.0, yj=N,
        level=fine_bits,
        domain_w=float(m), domain_h=float(n),
        bands=band_spec,
        airfoil_poly=airfoil_poly,
        points=points,
    )
    return points

def bands_from_chord_fractions(chord, pairs):
    """
    Convenience: specify bands as fractions of chord length.
    pairs: [(frac, stop_level), ...] sorted by frac ascending.
    """
    return [(chord*frac, lvl) for (frac, lvl) in pairs]


# ---------------------------
# Plotting demo
# ---------------------------

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
    # Build airfoil in domain coords (x right, y up — but plotting later inverts y)
    poly01 = naca4_airfoil(naca_code, n=airfoil_pts, closed_te=True)
    foil = transform_airfoil(poly01, chord=chord, alpha_deg=alpha_deg, origin=origin)

    # Make it explicitly closed for distance (wrap)
    if not np.allclose(foil[0], foil[-1]):
        foil = np.vstack([foil, foil[0]])

    # Bands from chord fractions → absolute units
    bands = bands_from_chord_fractions(chord, band_fracs)

    pts = generate_hilbert_airfoil_lod(n, m, foil, bands, extra_global_levels=extra_global_levels)
    xs = [p[1] for p in pts]
    ys = [p[0] for p in pts]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(xs, ys, '-', lw=1)
    ax.plot(foil[:, 0], foil[:, 1], lw=2)  # airfoil outline

    # Draw band radii as offsets from airfoil? (for reference we’ll just annotate text)
    ax.set_title(f"Hilbert with Quantized LOD around NACA {naca_code}")
    ax.set_xlim(0, m)
    ax.set_ylim(0, n)
    ax.set_aspect('equal')
    ax.invert_yaxis()  # make y downwards like raster grids

    # Pretty annotation
    txt = "\n".join([f"≤ {r:.2f} → stop={lvl}" if r < 1e8 else f"> last → stop={lvl}"
                     for (r, lvl) in bands])
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va='top', ha='left',
            bbox=dict(boxstyle="round,pad=0.3", fc="w", ec="0.5", alpha=0.9), fontsize=9)

    plt.show()


# ---------------------------
# Example
# ---------------------------
if __name__ == "__main__":
    # Domain (height n, width m)
    n, m = 256, 512

    # Airfoil placement: chord 45% of width, centered in height, slight incidence
    chord = 0.45 * m
    origin = (m*0.15, n*0.55)  # leading edge position
    alpha_deg = 6.0

    # LOD bands as fractions of chord (closer gets smaller stop_level = finer)
    band_fracs = [
        (0.06, 1),    # ultra-fine skin
        (0.15, 2),
        (0.3, 3),
        (1e9, 4),      # everything else
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