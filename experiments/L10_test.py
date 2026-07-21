"""Test scalar candidate L10 on synthetic configurations.

Omega: ellipse x = a cos t, y = b sin t (smooth convex).
w: Cauchy transform of point masses:  w(z) = sum_j a_j / (2 pi i (z - z_j)),
   sum a_j = 1 (a_j complex allowed), z_j inside Omega.
   Positivity required: dnu = 2 Re(w dsigma) >= 0 on the boundary (checked on grid).
p: polynomial, normalized ||p||_dOmega = 1, constrained sum_j a_j p(z_j) = 0
   (i.e. \oint p w = 0), enforced by subtracting a multiple of a correction poly.

Quantities:
  LHS = 2 \oint |p|^2 dnu - |\oint conj(p) w dsigma|^2
  c   = |sum_j a_j g(z_j) p(z_j)|,  g = Phi(p) (interior values via quadrature)
  RHS = (2 - c/2)^2
L10: LHS <= RHS ?
"""
import numpy as np

rng = np.random.default_rng(0)

def ellipse(a, b, N):
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    z = a * np.cos(t) + 1j * b * np.sin(t)
    zp = -a * np.sin(t) + 1j * b * np.cos(t)  # dz/dt
    return z, zp

def poly_z(c, z):
    out = np.full(np.shape(z), c[-1], dtype=complex)
    for k in range(len(c) - 2, -1, -1):
        out = out * z + c[k]
    return out

def phi_at(pvals, z, zp, pts):
    """g(z0) = (1/2pi i) oint conj(p(sigma))/(sigma - z0) dsigma for z0 strictly inside."""
    N = len(z)
    return np.array([np.sum(np.conj(pvals) * zp / (z - z0)) for z0 in pts]) \
        * (2 * np.pi / N) / (2j * np.pi)

def check_config(a_ell, b_ell, zs, ams, cp, N=4096, verbose=False):
    z, zp = ellipse(a_ell, b_ell, N)
    dt = 2 * np.pi / N
    w = np.zeros(N, dtype=complex)
    for zj, aj in zip(zs, ams):
        w += aj / (2j * np.pi * (z - zj))
    nu = 2 * np.real(w * zp) * dt          # dnu at nodes (>= 0 required)
    if nu.min() < -1e-10 * max(1.0, np.abs(nu).max()):
        return None  # positivity fails -> outside admissible class
    # normalize p on boundary
    pv = poly_z(cp, z)
    m = np.max(np.abs(pv))
    cp = [ci / m for ci in cp]
    pv = pv / m
    pz = poly_z(cp, np.array(zs))
    # constraint sum a_j p(z_j) = 0 must hold (caller enforces); verify
    con = abs(np.sum(np.array(ams) * pz))
    mass = np.sum(nu)
    lhs = 2 * np.sum(np.abs(pv) ** 2 * nu) - abs(np.sum(np.conj(pv) * w * zp) * dt) ** 2
    g_at = phi_at(pv, z, zp, np.array(zs))
    c_val = abs(np.sum(np.array(ams) * g_at * pz))
    rhs = (2 - c_val / 2) ** 2
    return dict(lhs=lhs, rhs=rhs, c=c_val, con=con, mass=mass,
                viol=lhs - rhs, numin=nu.min())

def project_constraint(cp, zs, ams, deg):
    """Adjust p |-> p - lambda*q so that sum a_j p(z_j) = 0, using q(z)=z^k with
    largest |sum a_j z_j^k| (k>=0)."""
    zsa = np.array(zs); ama = np.array(ams)
    s = np.sum(ama * poly_z(cp, zsa))
    best_k, best_m = 0, 0
    for k in range(deg + 1):
        mk = abs(np.sum(ama * zsa ** k))
        if mk > best_m:
            best_m, best_k = mk, k
    if best_m < 1e-12:
        return cp
    lam = s / np.sum(ama * zsa ** best_k)
    cp = list(cp)
    cp[best_k] -= lam
    return cp

def random_trial(deg=6, npts=3, shape=(1.0, 0.6)):
    a_ell, b_ell = shape
    # random interior points (inside a margin)
    zs = []
    while len(zs) < npts:
        x = rng.uniform(-a_ell, a_ell); y = rng.uniform(-b_ell, b_ell)
        if (x / a_ell) ** 2 + (y / b_ell) ** 2 < 0.7 ** 2:
            zs.append(x + 1j * y)
    kind = rng.integers(0, 2)
    if kind == 0:
        am = rng.uniform(0.1, 1, npts)
        am = am / am.sum()
    else:
        am = rng.uniform(0.1, 1, npts) * np.exp(1j * rng.uniform(-0.9, 0.9, npts))
        am = am / am.sum()
    cp = list(rng.standard_normal(deg + 1) + 1j * rng.standard_normal(deg + 1))
    cp = project_constraint(cp, zs, am, deg)
    return check_config(a_ell, b_ell, zs, am, cp)

if __name__ == "__main__":
    # sanity: Crabb-in-scalar-form (disk, point mass at 0, p=z^2) -> equality
    r = check_config(1.0, 1.0, [0.0], [1.0], [0, 0, 1.0])
    print(f"crabb-scalar: lhs={r['lhs']:.8f} rhs={r['rhs']:.8f} viol={r['viol']:+.2e}")

    worst = None
    tried = accepted = 0
    for shape in [(1.0, 0.9), (1.0, 0.6), (1.0, 0.3), (1.0, 0.15)]:
        for _ in range(400):
            tried += 1
            r = random_trial(deg=rng.integers(2, 9), npts=rng.integers(1, 5), shape=shape)
            if r is None or r["con"] > 1e-8:
                continue
            accepted += 1
            if worst is None or r["viol"] > worst["viol"]:
                worst = dict(r, shape=shape)
    print(f"accepted {accepted}/{tried}")
    print("worst:", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in worst.items()})
