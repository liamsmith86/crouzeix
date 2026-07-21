"""Targeted L10 stress test: doubly-symmetric ellipse, eta = (delta_{x} + delta_{-x})/2
Cauchy kernels, p odd with ||p||=1 on boundary and near-unimodular.

By symmetry: constraint oint p w = 0 automatic; m_p = 0; L10 demands
   2I <= 4 - 2|C| + |C|^2/4,  C = p(x) * Phi(p)(x)   (using g odd).
If I -> 2 (p near-unimodular) and C != 0 -> L10 violated.

p: odd polynomial approximating the Riemann map ellipse->disk (unimodular),
obtained by maximizing min |p| on boundary subject to max |p| = 1.
"""
import numpy as np
from scipy.optimize import minimize

from L10_test import ellipse, poly_z, phi_at
from scalar_adv import eval_config

def best_odd_unimodular(b_ell, deg_odd, N=2048, restarts=8, seed=3):
    """Maximize min|p| / max|p| over odd polynomials p(z) = sum_k c_k z^{2k+1}."""
    rng = np.random.default_rng(seed)
    z, _ = ellipse(1.0, b_ell, N)
    nodd = deg_odd // 2 + 1  # number of odd coefficients

    def to_poly(x):
        co = x[:nodd] + 1j * x[nodd:]
        cp = [0.0] * (2 * nodd)
        for k in range(nodd):
            cp[2 * k + 1] = co[k]
        return cp

    def obj(x):
        pv = np.abs(poly_z(to_poly(x), z))
        mx = pv.max()
        if mx < 1e-12:
            return 1.0
        return -(pv.min() / mx)

    best, bx = 1.0, None
    for _ in range(restarts):
        x0 = rng.standard_normal(2 * nodd)
        x0[0] = 1.0
        res = minimize(obj, x0, method="Nelder-Mead",
                       options={"maxiter": 8000, "fatol": 1e-12})
        if res.fun < best:
            best, bx = res.fun, res.x
    cp = to_poly(bx)
    pv = np.abs(poly_z(cp, z))
    cp = [c / pv.max() for c in cp]
    return cp, pv.min() / pv.max()

if __name__ == "__main__":
    for b_ell in [0.6, 0.3]:
        for deg_odd in [5, 9, 13]:
            cp, ratio = best_odd_unimodular(b_ell, deg_odd)
            z, zp = ellipse(1.0, b_ell, 4096)
            for x in [0.3, 0.5, 0.7, 0.85]:
                r = eval_config(b_ell, [x, -x], [0.5, 0.5], [0.0, 0.0], cp, N=8192)
                if r is None:
                    print(f"b={b_ell} deg={deg_odd} x={x}: inadmissible (nu < 0)")
                    continue
                # also report C decomposition
                pts = np.array([x, -x], dtype=complex)
                pvb = poly_z(cp, z)
                g_at = phi_at(pvb / np.max(np.abs(pvb)), z, zp, pts)
                pz = poly_z(cp, pts)
                print(f"b={b_ell} deg={deg_odd} minmod={ratio:.4f} x={x}: "
                      f"viol={r['viol']:+.6f} lhs={r['lhs']:.5f} rhs={r['rhs']:.5f} "
                      f"c={r['c']:.5f} |p(x)|={abs(pz[0]):.4f} |g(x)|={abs(g_at[0]):.5f} "
                      f"numin={r['numin']:.2e} con={r['con']:.2e}")
