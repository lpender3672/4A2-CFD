import numpy as np
import matplotlib.pyplot as plt


def rect_intersects_rect(ax0, ay0, ax1, ay1, bx0, by0, bx1, by1):
    return not (ax1 <= bx0 or bx1 <= ax0 or ay1 <= by0 or by1 <= ay0)

def rect_intersects_circle(rx0, ry0, rx1, ry1, cx, cy, r):
    # clamp circle center to rect
    nx = min(max(cx, rx0), rx1)
    ny = min(max(cy, ry0), ry1)
    dx, dy = nx - cx, ny - cy
    return (dx*dx + dy*dy) <= r*r

def hilbert_adaptive(x0, y0, xi, xj, yi, yj, level,
                     domain_w, domain_h, cx, cy, lod_radius,
                     points, coarse_stop_level=2):

    # Region AABB (in fine-grid coordinates)
    rx0 = min(x0, x0 + xi, x0 + yi, x0 + xi + yi)
    rx1 = max(x0, x0 + xi, x0 + yi, x0 + xi + yi)
    ry0 = min(y0, y0 + xj, y0 + yj, y0 + xj + yj)
    ry1 = max(y0, y0 + xj, y0 + yj, y0 + xj + yj)

    # If the region is completely outside the true n×m domain, skip
    if not rect_intersects_rect(rx0, ry0, rx1, ry1, 0, 0, domain_w, domain_h):
        return

    # Decide whether this region needs fine detail (intersects the LOD circle)
    #needs_fine = rect_intersects_circle(rx0, ry0, rx1, ry1, cx, cy, lod_radius)

    px = x0 + (xi + yi) * 0.5
    py = y0 + (xj + yj) * 0.5
    dist = np.hypot(px - cx, py - cy)

    needs_fine = (dist <= lod_radius)

    # Pick the stop level for this region
    stop_level = 1 if needs_fine else coarse_stop_level

    if level <= stop_level:
        # Emit the center of this region
        px = x0 + (xi + yi) * 0.5
        py = y0 + (xj + yj) * 0.5
        # Keep only if within domain
        if 0 <= px < domain_w and 0 <= py < domain_h:
            points.append((py, px))
        return

    # Recurse in Hilbert order (classic vector form)
    hilbert_adaptive(x0,                     y0,                      yi/2,  yj/2,  xi/2,  xj/2,  level-1,
                     domain_w, domain_h, cx, cy, lod_radius, points, coarse_stop_level)
    hilbert_adaptive(x0 + xi/2,             y0 + xj/2,               xi/2,  xj/2,  yi/2,  yj/2,  level-1,
                     domain_w, domain_h, cx, cy, lod_radius, points, coarse_stop_level)
    hilbert_adaptive(x0 + xi/2 + yi/2,      y0 + xj/2 + yj/2,        xi/2,  xj/2,  yi/2,  yj/2,  level-1,
                     domain_w, domain_h, cx, cy, lod_radius, points, coarse_stop_level)
    hilbert_adaptive(x0 + xi/2 + yi,        y0 + xj/2 + yj,         -yi/2, -yj/2, -xi/2, -xj/2, level-1,
                     domain_w, domain_h, cx, cy, lod_radius, points, coarse_stop_level)

def generate_hilbert_with_lod_recursive(n, m, lod_radius=2.0):
    # base bits for the rectangular grid (square pad to power of two)
    base_bits = int(np.ceil(np.log2(max(n, m))))
    # one extra level globally so fine cells exist
    fine_bits = base_bits + 1
    N = 1 << fine_bits  # fine-grid side length (covers the padded square)

    # Domain center in the same fine-grid units (no scaling needed)
    cx, cy = m / 2.0, n / 2.0

    points = []
    # Start with an axis-aligned square: columns along +x, rows along +y
    hilbert_adaptive(
        x0=0.0, y0=0.0,
        xi=float(N), xj=0.0,    # x-axis vector
        yi=0.0,     yj=float(N),# y-axis vector
        level=fine_bits,
        domain_w=float(m), domain_h=float(n),
        cx=cx, cy=cy, lod_radius=lod_radius,
        points=points,
        coarse_stop_level=2,     # => outside LOD stops at 2 (2×2 block → 1 coarse cell)
    )
    return points  # list of (y, x) in visiting order

def plot_hilbert_with_lod_recursive(n, m, lod_radius=2.0):
    pts = generate_hilbert_with_lod_recursive(n, m, lod_radius)
    xs = [p[1] for p in pts]
    ys = [p[0] for p in pts]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(xs, ys, '-o', lw=1, ms=3)
    ax.set_title(f"Recursive Hilbert with 2-level LOD (radius={lod_radius})")
    ax.set_xlim(0, m)
    ax.set_ylim(0, n)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    circle = plt.Circle((m/2, n/2), lod_radius, fill=False, lw=1.5, linestyle='--')
    ax.add_patch(circle)
    plt.show()

plot_hilbert_with_lod_recursive(64, 32, lod_radius=8)
