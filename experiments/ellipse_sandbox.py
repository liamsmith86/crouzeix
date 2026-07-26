"""Exact 2x2 sandbox with Omega = W(A): P2 falsification attempt.

Canonical family: A_h = [[-1, 2h],[0, 1]], W(A_h) = ellipse, foci +-1,
semi-axes a = sqrt(1+h^2), b = h. Every 2x2 matrix is affinely equivalent to
some A_h (affine invariance of the Crouzeix problem), so this sweep covers ALL
2x2 with the critical domain Omega = int W(A).

Riemann map ellipse -> D (classical): phi(z) = sqrt(k) sn((2K/pi) arcsin z, k),
modulus from K'(k)/K(k) = 4 xi0/pi, xi0 = artanh(b/a). Self-calibration check:
|phi| = 1 on boundary. Eigenvalue images phi(+-1) = +-sqrt(k).

Extremal f0 = b_alpha o phi (Blaschke deg <= 1, Crouzeix 04). For 2x2,
f(A_h) = [[f(-1), h(f(1)-f(-1))],[0, f(1)]] exactly -- no quadrature anywhere.

g0 = Phi(f0) = 1/f0 - r1/(z - z1) (principal part at the zero z1 = phi^{-1}(alpha)),
r1 = 1/f0'(z1) = psi-side: 1/(B'(alpha) phi'(z1)).

Quantities: K, C = <g0(A)f0(A)x0, x0>, W, beta, q, s_phase = 2 + ReC/2 - q.
"""

import numpy as np
from mpmath import (
    mp,
    mpf,
    mpc,
    ellipk,
    ellipf,
    ellipfun,
    asin,
    sin,
    sqrt,
    atanh,
)

mp.dps = 30

sn = ellipfun("sn")
cn = ellipfun("cn")
dn = ellipfun("dn")


def modulus_for(xi0, ratio_const=4.0):
    """Solve K(1-m)/K(m) = ratio_const*xi0/pi for m in (0,1)."""
    target = ratio_const * xi0 / mp.pi

    def f(m):
        return ellipk(1 - m) / ellipk(m) - target

    lo, hi = mpf("1e-12"), 1 - mpf("1e-12")
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


class EllipseMap:
    def __init__(self, h, ratio_const=4.0):
        self.h = mpf(h)
        self.a = sqrt(1 + self.h**2)
        self.b = self.h
        self.xi0 = atanh(self.b / self.a)
        self.m = modulus_for(self.xi0, ratio_const)
        self.K = ellipk(self.m)
        self.sk = sqrt(sqrt(self.m))  # sqrt(k), k = sqrt(m)

    def phi(self, z):
        u = (2 * self.K / mp.pi) * asin(z)
        return self.sk * sn(u, self.m)

    def dphi(self, z):
        u = (2 * self.K / mp.pi) * asin(z)
        return (
            self.sk
            * cn(u, self.m)
            * dn(u, self.m)
            * (2 * self.K / mp.pi)
            / sqrt(1 - z**2)
        )

    def phi_inv(self, w):
        """z with phi(z) = w (w in D)."""
        u = ellipf(asin(w / self.sk), self.m)
        return sin(mp.pi * u / (2 * self.K))

    def boundary_check(self, n=24):
        errs = []
        for t in np.linspace(0, 2 * np.pi, n, endpoint=False):
            z = self.a * mp.cos(t) + 1j * self.b * mp.sin(t)
            errs.append(abs(abs(self.phi(z)) - 1))
        return max(errs)


def f_of_A(fm1, fp1, h):
    """[[f(-1), h(f(1)-f(-1))],[0, f(1)]] as numpy complex 2x2."""
    return np.array(
        [
            [complex(fm1), complex(h) * (complex(fp1) - complex(fm1))],
            [0.0, complex(fp1)],
        ],
        dtype=complex,
    )


def mobius(w, al):
    return (w - al) / (1 - np.conj(al) * w)


def sandbox_eval(em, alpha):
    """All quantities for f0 = b_alpha o phi on A_h."""
    h = float(em.h)
    sk = complex(em.sk)
    al = complex(alpha)
    if abs(al) >= 0.9999:
        return None
    # f0 at eigenvalues
    f_m1 = mobius(-sk, al)
    f_p1 = mobius(sk, al)
    F = f_of_A(f_m1, f_p1, h)
    K_norm = np.linalg.norm(F, 2)
    # zero of f0: z1 = phi^{-1}(alpha)
    z1 = em.phi_inv(mpc(al))
    dphi_z1 = em.dphi(z1)
    Bp = 1.0 / (1 - abs(al) ** 2)  # b_alpha'(alpha)
    r1 = 1.0 / (Bp * complex(dphi_z1))
    z1c = complex(z1)
    # g0 at eigenvalues: 1/f0 - r1/(z - z1)
    g_m1 = 1.0 / f_m1 - r1 / (-1.0 - z1c)
    g_p1 = 1.0 / f_p1 - r1 / (1.0 - z1c)
    Gm = f_of_A(g_m1, g_p1, h)
    # svd
    U, s, Vh = np.linalg.svd(F)
    x0 = Vh[0].conj()
    diag = abs(np.vdot(x0, F @ x0))
    C = x0.conj() @ (Gm @ F @ x0)
    v = F @ x0 + Gm.conj().T @ x0
    W = np.linalg.norm(v)
    beta = x0.conj() @ v
    q = np.sqrt(max(W**2 - abs(beta) ** 2, 0.0))
    s_phase = 2 + C.real / 2 - q
    return dict(
        K=K_norm, C=C, W=W, beta=beta, q=q, s_phase=s_phase, diag=diag, z1=z1c, alpha=al
    )


def find_extremal_alpha(em, restarts=12, seed=0):
    from scipy.optimize import minimize

    rng = np.random.default_rng(seed)

    def negK(t):
        al = t[0] + 1j * t[1]
        al = al / np.sqrt(1 + abs(al) ** 2)  # squash to D
        r = sandbox_eval(em, al)
        return -(r["K"]) if r else 0.0

    best, bx = 0.0, None
    starts = [np.array([0.0, 0.0]), np.array([0.5, 0.0]), np.array([0.0, 0.5])]
    while len(starts) < restarts:
        starts.append(rng.standard_normal(2))
    for x0 in starts:
        res = minimize(
            negK,
            x0,
            method="Nelder-Mead",
            options={"maxiter": 3000, "fatol": 1e-14, "xatol": 1e-12},
        )
        if -res.fun > best:
            best, bx = -res.fun, res.x
    al = bx[0] + 1j * bx[1]
    return al / np.sqrt(1 + abs(al) ** 2)


if __name__ == "__main__":
    # calibrate the Riemann-map constant on a moderate ellipse
    for rc in (4.0, 2.0, 8.0):
        em = EllipseMap(0.8, ratio_const=rc)
        print(f"ratio_const={rc}: boundary defect = {float(em.boundary_check()):.2e}")
    print()
    RC = None
    best_def = 1e9
    for rc in (4.0, 2.0, 8.0):
        d = float(EllipseMap(0.8, ratio_const=rc).boundary_check())
        if d < best_def:
            best_def, RC = d, rc
    print(f"using ratio_const={RC} (defect {best_def:.2e})\n")

    print(
        f"{'h':>6} {'K':>9} {'ReC':>10} {'ImC':>10} {'q':>9} {'s_phase':>10} {'diag':>8}"
    )
    worst = None
    for h in [0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.2, 1.6, 2.2, 3.0]:
        em = EllipseMap(h, ratio_const=RC)
        al = find_extremal_alpha(em)
        r = sandbox_eval(em, al)
        print(
            f"{h:6.2f} {r['K']:9.6f} {r['C'].real:+10.6f} {r['C'].imag:+10.6f} "
            f"{r['q']:9.6f} {r['s_phase']:+10.6f} {r['diag']:8.1e}  alpha={r['alpha']:.4f}"
        )
        if worst is None or r["s_phase"] < worst["s_phase"]:
            worst = dict(r, h=h)
    print(f"\nworst s_phase over sweep: {worst['s_phase']:+.6f} at h={worst['h']}")
