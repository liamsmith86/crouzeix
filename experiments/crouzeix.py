"""Core library for Crouzeix ratio computation.

Conventions:
- Polynomials: coeffs c[0] + c[1] z + ... + c[d] z^d (lowest first), complex.
- R(A,p) = ||p(A)||_2 / max_{z in W(A)} |p(z)|.

Approximation directions (critical for rigor):
- INNER: boundary points z_k = x_k* A x_k from eigenvectors of Re(e^{-i th}A)
  lie in W(A). max_k |p(z_k)| <= true max  =>  R_inner >= R_true.
  R_inner is an OVERestimate of R: good for steering searches, never for claims.
- OUTER: W(A) is contained in the polygon cut by support half-planes
  Re(e^{-i th_k} z) <= lambda_max(Re(e^{-i th_k} A)). max over that polygon
  >= true max  =>  R_outer <= R_true. R_outer > 2 is a genuine violation
  (up to floating point, to then be certified independently).
"""
import numpy as np

def poly_A(c, A):
    """p(A) by Horner, coeffs lowest-first."""
    n = A.shape[0]
    P = np.eye(n, dtype=complex) * c[-1]
    for k in range(len(c) - 2, -1, -1):
        P = P @ A + c[k] * np.eye(n)
    return P

def poly_z(c, z):
    z = np.asarray(z, dtype=complex)
    out = np.full(z.shape, c[-1], dtype=complex)
    for k in range(len(c) - 2, -1, -1):
        out = out * z + c[k]
    return out

def opnorm(M):
    return np.linalg.norm(M, 2)

def nr_support(A, thetas):
    """Support data of W(A): for each theta return (h(theta), z(theta)) where
    h = lambda_max(Re(e^{-i th}A)) is the support function and z = x*Ax is an
    attaining boundary point (inner)."""
    hs = np.empty(len(thetas))
    zs = np.empty(len(thetas), dtype=complex)
    for i, th in enumerate(thetas):
        H = (np.exp(-1j * th) * A + np.exp(1j * th) * A.conj().T) / 2
        w, V = np.linalg.eigh(H)
        hs[i] = w[-1]
        x = V[:, -1]
        zs[i] = x.conj() @ A @ x
    return hs, zs

def inner_boundary_points(A, ntheta=256):
    thetas = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    _, zs = nr_support(A, thetas)
    return zs

def outer_polygon(A, ntheta=256):
    """Vertices of the circumscribing polygon from ntheta support lines.
    Vertex v_k solves Re(e^{-i th_k} v)=h_k, Re(e^{-i th_{k+1}} v)=h_{k+1}."""
    thetas = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    hs, _ = nr_support(A, thetas)
    verts = []
    for k in range(ntheta):
        t1, t2 = thetas[k], thetas[(k + 1) % ntheta]
        h1, h2 = hs[k], hs[(k + 1) % ntheta]
        # Solve [cos t1, sin t1; cos t2, sin t2] [x;y] = [h1;h2]
        M = np.array([[np.cos(t1), np.sin(t1)], [np.cos(t2), np.sin(t2)]])
        det = np.linalg.det(M)
        if abs(det) < 1e-14:
            continue
        xy = np.linalg.solve(M, np.array([h1, h2]))
        verts.append(xy[0] + 1j * xy[1])
    return np.array(verts)

def max_abs_on_segments(c, verts, per_edge=64):
    """Max |p| over closed polygon boundary via per-edge sampling + one
    refinement pass. (Sampling-based; for certification use interval version.)"""
    best = 0.0
    m = len(verts)
    for k in range(m):
        a, b = verts[k], verts[(k + 1) % m]
        ts = np.linspace(0, 1, per_edge)
        pts = a + (b - a) * ts
        vals = np.abs(poly_z(c, pts))
        j = int(np.argmax(vals))
        # local refine around best t
        lo, hi = max(0, j - 1), min(per_edge - 1, j + 1)
        ts2 = np.linspace(ts[lo], ts[hi], 64)
        pts2 = a + (b - a) * ts2
        best = max(best, float(np.max(np.abs(poly_z(c, pts2)))), float(vals[j]))
    return best

def ratio_inner(A, c, ntheta=256):
    """Overestimate of R (search guidance)."""
    zs = inner_boundary_points(A, ntheta)
    M = float(np.max(np.abs(poly_z(c, zs))))
    if M == 0:
        return np.inf
    return opnorm(poly_A(c, A)) / M

def ratio_outer(A, c, ntheta=512, per_edge=64):
    """Underestimate of R (safe direction for claiming R > 2)."""
    verts = outer_polygon(A, ntheta)
    M = max_abs_on_segments(c, verts, per_edge)
    if M == 0:
        return np.inf
    return opnorm(poly_A(c, A)) / M

def crabb_matrix(k):
    """(k+1)x(k+1) nilpotent with superdiagonal (sqrt2, 1, ..., 1, sqrt2):
    reported to satisfy W(A)=unit disk, ||A^k||=2 (Crabb / Greenbaum-Overton)."""
    n = k + 1
    A = np.zeros((n, n), dtype=complex)
    w = [np.sqrt(2)] + [1.0] * (k - 2) + [np.sqrt(2)] if k >= 2 else [2.0]
    for i, wi in enumerate(w):
        A[i, i + 1] = wi
    return A
