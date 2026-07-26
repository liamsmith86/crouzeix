"""Verify the candidate L7' violation from advL7_n3 and test L7 at the TRUE
extremal pair of the same (M, Omega).

Steps:
1. Load best record; re-evaluate at N=4096 with strict certificates.
2. Enforce <p(M)x0,x0> = 0 EXACTLY: fix p, adjust x0 by moving within the span
   of x0 and a correction direction (small rotation), re-evaluate.
3. Compute the true extremal polynomial f0 (deg 12, strong multistart) for the
   same (M, Omega); test L7 at (f0, x0=top-sv).
"""
import json

import numpy as np
from scipy.optimize import minimize

from crouzeix import poly_z, poly_A
from adversarial_L7 import quantities
from refined_test import find_extremal_poly, refined_quantities
from star_inequality import offset_curve

EPS = 0.05

recs = [json.loads(line) for line in open("advL7_n3_d4_s12.jsonl")]
r = max(recs, key=lambda q: q["viol"])
M = np.array(r["A_re"]) + 1j * np.array(r["A_im"])
cp = list(np.array(r["c_re"]) + 1j * np.array(r["c_im"]))
x0 = np.array(r["x0_re"]) + 1j * np.array(r["x0_im"])
x0 = x0 / np.linalg.norm(x0)

q = quantities(M, cp, x0, EPS, N=4096, tol=1e-9)
print("step1 strict re-eval:", "REJECTED" if q is None else
      f"lhs={q[0]:.6f} rhs={q[1]:.6f} viol={q[0]-q[1]:+.6f} c={q[2]:.6f} diag={q[3]:.3e}")

# step 2: exact orthogonality. Solve for x in span{x0, v} with <p(M)x,x>=0,
# staying close to x0. Use v = eigvec direction reducing the form.
pM = poly_A([c for c in cp], M)  # note: quantities renormalizes internally; do same
_, z_out, _ = offset_curve(M, EPS, 4096)
m = np.max(np.abs(poly_z(cp, z_out)))
cpn = [c / m for c in cp]
pM = poly_A(cpn, M)

def orth_adjust(pM, x0):
    """Find unit x near x0 with <pM x, x> = 0 via optimizing in full space."""
    n = len(x0)
    def f(t):
        x = t[:n] + 1j * t[n:]
        nx = np.linalg.norm(x)
        if nx < 1e-9:
            return 1e6
        x = x / nx
        return abs(np.vdot(x, pM @ x)) ** 2 * 1e4 + np.linalg.norm(x - x0) ** 2
    t0 = np.concatenate([x0.real, x0.imag])
    res = minimize(f, t0, method="Nelder-Mead",
                   options={"maxiter": 20000, "fatol": 1e-16, "xatol": 1e-12})
    x = res.x[:len(x0)] + 1j * res.x[len(x0):]
    return x / np.linalg.norm(x)

x1 = orth_adjust(pM, x0)
q1 = quantities(M, cp, x1, EPS, N=4096, tol=1e-9)
print("step2 exact-orth  :", "REJECTED" if q1 is None else
      f"lhs={q1[0]:.6f} rhs={q1[1]:.6f} viol={q1[0]-q1[1]:+.6f} c={q1[2]:.6f} diag={q1[3]:.3e}")

# step 3: true extremal for (M, Omega)
rng = np.random.default_rng(99)
c0, K = find_extremal_poly(M, EPS, 12, rng, restarts=80, N=384)
qq = refined_quantities(M, c0, eps=EPS, N=4096)
print(f"step3 true-extremal deg12: K={qq['K']:.6f} W={qq['W']:.6f} |b|={abs(qq['beta']):.6f} "
      f"c={qq['c']:.6f} diag={qq['diag']:.3e}")
print(f"      L7 at extremal: {qq['lhs']:.6f} <= {qq['rhs']:.6f} "
      f"{'HOLDS' if qq['lhs'] <= qq['rhs'] + 1e-9 else 'FAILS'} (slack {qq['rhs']-qq['lhs']:+.6f})")
print(f"      L6 check: K^2={qq['K']**2:.6f} <= K*lhs+c={qq['K']*qq['lhs']+qq['c']:.6f} "
      f"(+|b|diag corr {abs(qq['beta'])*qq['diag']:.2e})")
