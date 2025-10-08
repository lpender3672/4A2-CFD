import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------
# Morton / Z-order utilities
# -------------------------------------------------------
def expand_bits(v: int) -> int:
    """Expand 16-bit integer into 32 bits by inserting zeros between bits."""
    v = (v | (v << 8)) & 0x00FF00FF
    v = (v | (v << 4)) & 0x0F0F0F0F
    v = (v | (v << 2)) & 0x33333333
    v = (v | (v << 1)) & 0x55555555
    return v

def morton2D(i: int, j: int) -> int:
    """Compute 2D Morton code (Z-order curve index) for integer coordinates (i, j)."""
    return (expand_bits(i) << 1) | expand_bits(j)

# -------------------------------------------------------
# Generate base + LOD grid traversal
# -------------------------------------------------------
def generate_zorder_with_lod(n: int, m: int, lod_radius: float = 2.0):
    cx, cy = n / 2, m / 2
    coords = []

    # Base Morton traversal of all coarse cells
    base_cells = [(i, j, morton2D(i, j)) for i in range(n) for j in range(m)]
    base_cells.sort(key=lambda x: x[2])

    for i, j, code in base_cells:
        # Cell center
        x_center = j + 0.5
        y_center = i + 0.5
        dist = np.hypot(x_center - cx, y_center - cy)

        if dist <= lod_radius:
            # Subdivide this cell into 2×2 subcells in world coordinates
            for fi in range(2):
                for fj in range(2):
                    sub_x_center = j + (fj + 0.5) / 2
                    sub_y_center = i + (fi + 0.5) / 2
                    local_code = morton2D(fi, fj)
                    coords.append((sub_y_center, sub_x_center, code * 4 + local_code))
        else:
            # Keep coarse cell as-is
            coords.append((y_center, x_center, code * 4))

    coords.sort(key=lambda x: x[2])
    return coords

# -------------------------------------------------------
# Plot function
# -------------------------------------------------------
def plot_zorder_with_lod(n: int, m: int, lod_radius: float = 2.0):
    coords = generate_zorder_with_lod(n, m, lod_radius)
    xs = [c[1] for c in coords]
    ys = [c[0] for c in coords]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(xs, ys, '-o', lw=1, ms=3)
    ax.set_title(f"Z-Order Curve with Local LOD (radius={lod_radius}) — no overlap")
    ax.set_xlim(0, m)
    ax.set_ylim(0, n)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # draw LOD region
    circle = plt.Circle((m / 2, n / 2), lod_radius, color='r', fill=False, lw=1.5, linestyle='--')
    ax.add_patch(circle)

def rot(n, x, y, rx, ry):
    """Rotate/flip a quadrant appropriately."""
    if ry == 0:
        if rx == 1:
            x = n - 1 - x
            y = n - 1 - y
        x, y = y, x
    return x, y


def xy_to_hilbert(x, y, bits):
    """Return index and orientation (rotation id 0–3)."""
    d = 0
    n = 1 << bits
    s = n >> 1
    rot_id = 0  # 0–3 (orientation)
    while s > 0:
        rx = 1 if (x & s) > 0 else 0
        ry = 1 if (y & s) > 0 else 0
        d += s * s * ((3 * rx) ^ ry)
        x, y = rot(s, x, y, rx, ry)
        # track rotation
        if ry == 0:
            if rx == 1:
                rot_id ^= 1
            rot_id ^= 3  # quarter-turn swap
        s >>= 1
    return d, rot_id % 4


def oriented_hilbert_subcells(i, j, orientation):
    """Return 4 subcells (2×2) in Hilbert order for a given orientation."""
    # canonical Hilbert order (2×2)
    base = [(0,0), (0,1), (1,1), (1,0)]
    # rotate pattern according to parent orientation
    for _ in range(orientation):
        base = [(y, 1 - x) for x, y in base]
    return [(i + fi / 2 + 0.25, j + fj / 2 + 0.25) for fi, fj in base]


# -------------------------------------------------------
# Generate Hilbert traversal with local LOD
# -------------------------------------------------------
def generate_hilbert_with_lod(n, m, lod_radius=2.0):
    cx, cy = n / 2, m / 2
    coords = []

    bits = int(np.ceil(np.log2(max(n, m))))
    base_cells = []
    for i in range(n):
        for j in range(m):
            code, rot_id = xy_to_hilbert(i, j, bits)
            base_cells.append((i, j, code, rot_id))

    base_cells.sort(key=lambda x: x[2])

    for i, j, code, orientation in base_cells:
        x_center = j + 0.5
        y_center = i + 0.5
        dist = np.hypot(x_center - cx, y_center - cy)

        if dist <= lod_radius:
            # Subdivide with proper orientation
            subcells = oriented_hilbert_subcells(i, j, orientation)
            for k, (y_sub, x_sub) in enumerate(subcells):
                coords.append((y_sub, x_sub, code * 4 + k))
        else:
            coords.append((y_center, x_center, code * 4))

    coords.sort(key=lambda x: x[2])
    return coords


# -------------------------------------------------------
# Plot function
# -------------------------------------------------------
def plot_hilbert_with_lod(n, m, lod_radius=2.0):
    coords = generate_hilbert_with_lod(n, m, lod_radius)
    xs = [c[1] for c in coords]
    ys = [c[0] for c in coords]

    fig, ax = plt.subplots(figsize=(7,7))
    ax.plot(xs, ys, '-o', lw=1, ms=3)
    ax.set_title(f"Hilbert Curve with Local LOD (radius={lod_radius})")
    ax.set_xlim(0, m)
    ax.set_ylim(0, n)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.add_patch(plt.Circle((m/2, n/2), lod_radius, color='r', fill=False, lw=1.5, linestyle='--'))
    plt.show()

plot_zorder_with_lod(20, 20, lod_radius=3.0)

plot_hilbert_with_lod(5, 5, lod_radius=1.5)

plt.show()


