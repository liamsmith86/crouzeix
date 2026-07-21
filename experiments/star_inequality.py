"""Test inequality (*): c <= 2L at near-extremal configurations.

  c = |<(Phi(p0) p0)(M) x0, x0>|,  L = contour integral of lambda_min(P(sigma)) ds,
where Omega = outward eps-offset of W(M), p0 = (near-)extremal polynomial for
K(M,Omega) = sup{||p(M)|| : ||p||_Omega <= 1}, x0 = top right singular vector.

Theory (Schwenninger-de Vries Thm 5 + Caldwell-Greenbaum-Li Prop 8):
  K <= (1 - L/2) + sqrt((1 - L/2)^2 + c),  and  c <= 2L  ==>  K <= 2.
If (*) holds at all extremal configurations, Crouzeix's conjecture follows.

Numerics:
- boundary of W(M): z_W(theta) = x_theta* M x_theta (strictly convex generic case);
- Omega boundary: z(theta) = z_W(theta) + eps*e^{i theta}; z' by FFT differentiation;
- P(theta)*|z'| = (1/2pi i) z'(theta) (z(theta)-M)^{-1} + h.c.;
- Phi(p0) evaluated on an INNER contour (eps/2-offset) by plain quadrature over
  the outer contour (no principal values needed);
- (Phi(p0)*p0)(M) via Cauchy integral over the inner contour (spectrum of M is
  inside W(M), hence inside the inner contour).
"""
import numpy as np

from crouzeix import nr_support, poly_z, poly_A, opnorm

def boundary_W(M, N):
    thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
    _, zs = nr_support(M, thetas)
    return thetas, zs

def offset_curve(M, eps, N):
    """z(theta) = boundary of W(M) offset by eps along e^{i theta}; z' via FFT."""
    thetas, zw = boundary_W(M, N)
    z = zw + eps * np.exp(1j * thetas)
    k = np.fft.fftfreq(N, d=1.0 / N)  # integer frequencies
    zp = np.fft.ifft(1j * k * np.fft.fft(z))
    return thetas, z, zp

def dlp_min_eig_integral(M, z, zp):
    """L = sum over theta grid of lambda_min(Q(theta)) * dtheta, Q = P*|z'|."""
    N = len(z)
    n = M.shape[0]
    I = np.eye(n, dtype=complex)
    L = 0.0
    lmins = np.empty(N)
    for j in range(N):
        Rj = np.linalg.inv(z[j] * I - M)
        B = (zp[j] / (2j * np.pi)) * Rj
        Q = B + B.conj().T
        lmins[j] = np.linalg.eigvalsh(Q)[0]
        L += lmins[j]
    L *= 2 * np.pi / N
    # sanity: total mass integral should be 2I
    S = np.zeros((n, n), dtype=complex)
    for j in range(N):
        Rj = np.linalg.inv(z[j] * I - M)
        B = (zp[j] / (2j * np.pi)) * Rj
        S += B + B.conj().T
    S *= 2 * np.pi / N
    mass_err = opnorm(S - 2 * np.eye(n))
    return L, lmins, mass_err

def cauchy_matrix(fvals, z, zp, M):
    """(1/2pi i) sum f(z_j) (z_j - M)^{-1} z'_j dtheta  (f holomorphic inside)."""
    N = len(z)
    n = M.shape[0]
    I = np.eye(n, dtype=complex)
    S = np.zeros((n, n), dtype=complex)
    for j in range(N):
        S += fvals[j] * zp[j] * np.linalg.inv(z[j] * I - M)
    return S * (2 * np.pi / N) / (2j * np.pi)

def phi_on_inner(pvals_outer, z_out, zp_out, z_in):
    """Phi(p)(w) = (1/2pi i) oint conj(p(sigma)) (sigma - w)^{-1} dsigma, w on inner contour."""
    N = len(z_out)
    g = np.conj(pvals_outer)
    out = np.empty(len(z_in), dtype=complex)
    for i, w in enumerate(z_in):
        out[i] = np.sum(g * zp_out / (z_out - w)) * (2 * np.pi / N) / (2j * np.pi)
    return out

def analyze(M, c_poly, eps=0.05, N=1024, label=""):
    """Given near-extremal polynomial coeffs (lowest first) for domain Omega_eps,
    compute K, c, L and check (*)."""
    th, z_out, zp_out = offset_curve(M, eps, N)
    _, z_in, zp_in = offset_curve(M, eps / 2, N)
    # normalize p on outer boundary
    pv_out = poly_z(c_poly, z_out)
    m = np.max(np.abs(pv_out))
    c_poly = [ci / m for ci in c_poly]
    pv_out = pv_out / m
    pM = poly_A(c_poly, M)
    K = opnorm(pM)
    # x0 = top right singular vector
    U, s, Vh = np.linalg.svd(pM)
    x0 = Vh[0].conj()
    diag_orth = abs(np.vdot(x0, pM @ x0))  # |<p(M)x0, x0>|, ~0 at true extremal
    # Phi(p) on inner contour, then (Phi(p)*p)(M)
    phi_in = phi_on_inner(pv_out, z_out, zp_out, z_in)
    pv_in = poly_z(c_poly, z_in)
    prodM = cauchy_matrix(phi_in * pv_in, z_in, zp_in, M)
    c_val = abs(np.vdot(x0, prodM @ x0))
    # sanity: reconstruct p(M) from inner contour
    pM2 = cauchy_matrix(pv_in, z_in, zp_in, M)
    rec_err = opnorm(pM - pM2)
    L, lmins, mass_err = dlp_min_eig_integral(M, z_out, zp_out)
    bound = (1 - L / 2) + np.sqrt((1 - L / 2) ** 2 + c_val)
    print(f"[{label}] K={K:.6f}  c={c_val:.6f}  2L={2*L:.6f}  "
          f"(*) {'HOLDS' if c_val <= 2*L else 'FAILS'}  bound={bound:.6f}  "
          f"|<pMx0,x0>|={diag_orth:.2e}  mass_err={mass_err:.2e}  rec_err={rec_err:.2e}")
    return dict(K=K, c=c_val, L=L, bound=bound, diag=diag_orth,
                mass_err=mass_err, rec_err=rec_err)

def find_extremal_poly(M, eps, deg, rng, restarts=30, N=256):
    """Maximize ||p(M)|| / max_{Omega_eps}|p| over poly coeffs."""
    from scipy.optimize import minimize
    _, z_out, _ = offset_curve(M, eps, N)

    def neg_ratio(x):
        c = x[:deg + 1] + 1j * x[deg + 1:]
        m = np.max(np.abs(poly_z(c, z_out)))
        if m < 1e-14:
            return 0.0
        return -opnorm(poly_A(list(c), M)) / m

    best, bestc = 0.0, None
    for _ in range(restarts):
        x0 = rng.standard_normal(2 * (deg + 1))
        res = minimize(neg_ratio, x0, method="L-BFGS-B",
                       options={"maxiter": 300, "eps": 1e-8})
        if -res.fun > best:
            best = -res.fun
            bestc = res.x[:deg + 1] + 1j * res.x[deg + 1:]
    return list(bestc), best

if __name__ == "__main__":
    rng = np.random.default_rng(42)
    from crouzeix import crabb_matrix

    # Case 1: Crabb k=2 (W = unit disk), p = z^2. Expect c ~ 0, K ~ 2/(1+eps)^2.
    M = crabb_matrix(2)
    analyze(M, [0, 0, 1], eps=0.05, label="crabb2 z^2 eps=.05")
    analyze(M, [0, 0, 1], eps=0.01, label="crabb2 z^2 eps=.01")

    # Case 2: 2x2 elliptic W: M = [[0,1],[0,1]].
    M = np.array([[0, 1], [0, 1]], dtype=complex)
    c0, r = find_extremal_poly(M, 0.05, 4, rng)
    analyze(M, c0, eps=0.05, label=f"ell2x2 deg4 (r={r:.4f})")

    # Case 3: random 3x3, optimized polynomial.
    for t in range(3):
        M = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        c0, r = find_extremal_poly(M, 0.05, 5, rng)
        analyze(M, c0, eps=0.05, label=f"rand3x3-{t} deg5 (r={r:.4f})")
