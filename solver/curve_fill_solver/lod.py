# generate level of detail based on curvature smoothing
import numpy as np



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

def gaussian_kernel_2d(size: int, sigma: float) -> np.ndarray:
    """Generate a 2D Gaussian kernel."""
    ax = np.arange(-size // 2 + 1, size // 2 + 1)
    xx, yy = np.meshgrid(ax, ax)
    kernel = np.exp(-(xx**2 + yy**2) / (2.0 * sigma**2))
    kernel /= kernel.sum()
    return kernel

def smooth2d(data: np.ndarray, size: int = 5, sigma: float = 1.0) -> np.ndarray:
    """Apply Gaussian smoothing using pure NumPy (manual convolution)."""
    kernel = gaussian_kernel_2d(size, sigma)
    pad = size // 2
    padded = np.pad(data, pad, mode='reflect')

    out = np.zeros_like(data)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            region = padded[i:i+size, j:j+size]
            out[i, j] = np.sum(region * kernel)
    return out

def calc_lod(n, m, foil, max_level):

    x = np.linspace(0, m, m)
    y = np.linspace(0, n, n)

    X,Y = np.meshgrid(x, y)

    # Generate LOD grid
    KAPPA = np.zeros_like(X)
    DIST = np.zeros_like(X)
    LOD = np.zeros_like(X, dtype=np.int32)
    for i in range(n):
        for j in range(m):
            px = X[i, j]
            py = Y[i, j]
            dist, kappa = curvature_at_nearest_point(px, py, foil)
            DIST[i, j] = dist
            KAPPA[i, j] = kappa

    LKAPPA = np.log(KAPPA / np.min(KAPPA[KAPPA>0]))  # avoid log(0)
    IDIST = 1 - DIST / np.max(DIST)

    LOD = IDIST * LKAPPA
    
    NLOD = smooth2d(LOD, size=101, sigma=10.0)

    NLOD = NLOD - np.min(NLOD)
    NLOD = (NLOD / np.max(NLOD) * max_level)
    ILOD = max_level - np.clip(np.round(NLOD), 0, max_level)

    return X, Y, ILOD

def demo():

    import matplotlib.pyplot as plt

    n = 400
    m = 1000
    airfoil_poly = transform_airfoil(
        naca4_airfoil("2412", n=201, closed_te=True),
        chord=500.0,
        alpha_deg=2.0,
        origin=(250.0, 200.0),
    )

    fine_bits = 4

    X, Y, LOD = calc_lod(n, m, airfoil_poly, fine_bits)

    plt.figure(figsize=(10,6))
    plt.contourf(X, Y, LOD, levels=fine_bits+1, cmap='viridis_r')
    plt.colorbar(label='Level of Detail')
    plt.plot(airfoil_poly[:,0], airfoil_poly[:,1], 'k-', linewidth=2)
    plt.title('Level of Detail based on Airfoil Curvature')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.show()
    

if __name__ == "__main__":
    demo()
    