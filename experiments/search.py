"""Track C: nonsmooth maximization of the Crouzeix ratio (Greenbaum-Overton style).

Maximize R_inner(A, p) by L-BFGS with finite differences, many restarts.
R_inner overestimates R, so no candidate is missed; candidates are re-verified
with ratio_outer at high resolution before being logged as near-misses.

Usage: search.py [n] [deg] [restarts] [seed]
Writes JSON lines to results_n{n}_d{deg}_s{seed}.jsonl
"""
import json
import sys

import numpy as np
from scipy.optimize import minimize

from crouzeix import ratio_inner, ratio_outer

def pack(A, c):
    return np.concatenate([A.real.ravel(), A.imag.ravel(), np.real(c), np.imag(c)])

def unpack(x, n, d):
    n2 = n * n
    A = (x[:n2] + 1j * x[n2:2 * n2]).reshape(n, n)
    c = x[2 * n2:2 * n2 + d + 1] + 1j * x[2 * n2 + d + 1:]
    return A, c

def negR(x, n, d, ntheta):
    A, c = unpack(x, n, d)
    if np.max(np.abs(c[1:])) < 1e-12:  # constant polynomial => R=1 trivially
        return -1.0
    r = ratio_inner(A, c, ntheta=ntheta)
    if not np.isfinite(r):
        return 0.0
    return -r

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    restarts = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    rng = np.random.default_rng(seed)
    fname = f"results_n{n}_d{d}_s{seed}.jsonl"
    best_overall = 0.0
    with open(fname, "w") as fh:
        for trial in range(restarts):
            # varied initializations: random, near-Crabb perturbations, nilpotent-ish
            kind = trial % 3
            if kind == 0:
                A0 = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
            elif kind == 1:
                A0 = np.diag(rng.standard_normal(n - 1) + 1.0, 1).astype(complex)
                A0 += 0.3 * (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
            else:
                A0 = np.triu(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)), 1)
                A0 += 0.1 * rng.standard_normal((n, n))
            c0 = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
            x0 = pack(A0, c0)
            res = minimize(negR, x0, args=(n, d, 64), method="L-BFGS-B",
                           options={"maxiter": 400, "eps": 1e-7})
            r64 = -res.fun
            A, c = unpack(res.x, n, d)
            # verify at high resolution, both directions
            ri = ratio_inner(A, c, ntheta=1024)
            ro = ratio_outer(A, c, ntheta=1024)
            rec = {"trial": trial, "kind": kind, "r64": r64, "r_inner1024": ri,
                   "r_outer1024": ro,
                   "A_re": A.real.tolist(), "A_im": A.imag.tolist(),
                   "c_re": np.real(c).tolist(), "c_im": np.imag(c).tolist()}
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if ri > best_overall:
                best_overall = ri
                print(f"[n={n} d={d} s={seed}] trial {trial}: R_inner={ri:.8f} R_outer={ro:.8f}", flush=True)
    print(f"DONE n={n} d={d} seed={seed}: best R_inner(1024) = {best_overall:.10f}", flush=True)

if __name__ == "__main__":
    main()
