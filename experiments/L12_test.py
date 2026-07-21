"""Test L12: at true extremal pairs (K>1), is <(f0 h)(A) x0, x0> = 0 for all h in A(Omega)?
(h=1 is the known orthogonality E2. h = g0 gives c = 0.)

Also: refine extremal precision and watch whether c ~ 0 tightens.
Domain: psi = w + eps w^2, A = tau I + s Crabb(m) (+ optional perturbation).
"""
import numpy as np
from scipy.optimize import minimize

from crouzeix import crabb_matrix
from extremal_pullback import Pullback, find_extremal_blaschke, blaschke
from epoch2_neardisk import max_fill

def refine_blaschke(pb, A, al0, iters=3):
    T = pb.resolvent_stack(A)
    deg = len(al0)

    def to_x(al):
        pts = [a / np.sqrt(1 - abs(a) ** 2 + 1e-15) for a in al]
        return np.array([p.real for p in pts] + [p.imag for p in pts])

    def unpack(x):
        pts = x[:deg] + 1j * x[deg:]
        return [p / np.sqrt(1 + abs(p) ** 2) for p in pts]

    def obj(x):
        return -np.linalg.norm(pb.calc(blaschke(pb.w, unpack(x)), A, T), 2)

    x = to_x(al0)
    for _ in range(iters):
        res = minimize(obj, x, method="Nelder-Mead",
                       options={"maxiter": 20000, "fatol": 1e-16, "xatol": 1e-14})
        x = res.x
    return unpack(x)

def test_L12(pb, A, alphas):
    T = pb.resolvent_stack(A)
    d = pb.extremal_data(A, alphas, T)
    x0 = d["x0"]
    Bv = blaschke(pb.w, alphas)
    rng = np.random.default_rng(1)
    out = []
    tests = [("w^%d" % k, pb.w ** k) for k in range(0, 5)]
    hr = np.polyval(rng.standard_normal(4) + 1j * rng.standard_normal(4), pb.w)
    tests.append(("rand-poly", hr / np.max(np.abs(hr))))
    for name, hv in tests:
        val = x0.conj() @ (pb.calc(Bv * hv, A, T) @ x0)
        out.append((name, abs(val)))
    return d, out

if __name__ == "__main__":
    for m, eps in [(2, 0.15), (2, 0.25), (3, 0.25)]:
        psi = lambda w, e=eps: w + e * w * w
        dpsi = lambda w, e=eps: 1 + 2 * e * w
        pb = Pullback(psi, dpsi, N=2048)
        tau, s, A = max_fill(pb, crabb_matrix(m), np.linspace(0, 0.3, 7))
        al, K1 = find_extremal_blaschke(pb, A, m, restarts=24, seed=3)
        al = refine_blaschke(pb, A, al)
        d, out = test_L12(pb, A, al)
        print(f"m={m} eps={eps}: K={d['K']:.8f} diag={d['diag']:.2e} c={d['c']:.3e} "
              f"G={d['G']:.6f} |b|={abs(d['beta']):.6f} slack={d['slack']:+.6f}")
        print("   L12 <(f0 h)(A)x0,x0>:", "  ".join(f"{n}:{v:.2e}" for n, v in out))
        # zeros of the extremal Blaschke
        print("   zeros:", np.round(np.array(al), 6))
