"""Adversarial H-r falsification: minimize ReC over A (n x n) at near-critical
domain Omega = W(A) + t-offset (t small). Random-direction descent with restarts.

Usage: Hr_adversarial.py n t iters seed
Certificates per eval: theodorsen convergence, DLP positivity, unitality, diag < 1e-5, K > 1.02.
"""
import json
import sys

import numpy as np

from extremal_pullback import find_extremal_blaschke
from theodorsen import theodorsen_map, GeneralPullback

def eval_A(A, t, seed, N=768, restarts=14):
    A = A / np.linalg.norm(A, 2)
    try:
        z, zp, terr = theodorsen_map(A, N=N, inflate=t)
    except Exception:
        return None
    if terr > 1e-8:
        return None
    pb = GeneralPullback(z, zp)
    lmin, merr = pb.dlp_certificate(A)
    if lmin < -1e-9 or merr > 1e-6 or pb.unitality_certificate() > 1e-8:
        return None
    al, _ = find_extremal_blaschke(pb, A, A.shape[0] - 1, restarts=restarts, seed=seed)
    d = pb.extremal_data(A, al)
    if d["K"] < 1.02 or d["diag"] > 1e-5:
        return None
    return d

def main():
    n = int(sys.argv[1]); t = float(sys.argv[2])
    iters = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    from crouzeix import crabb_matrix
    if seed % 2 == 0:
        A = crabb_matrix(n - 1) + 0.4 * (rng.standard_normal((n, n))
                                         + 1j * rng.standard_normal((n, n)))
    else:
        A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    cur = eval_A(A, t, seed)
    while cur is None:
        A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        cur = eval_A(A, t, seed)
    best_r = cur["C"].real
    step = 0.2
    fname = f"Hr_adv_n{n}_t{t}_s{seed}.jsonl"
    with open(fname, "w") as fh:
        rec = dict(it=-1, ReC=best_r, K=cur["K"], s_phase=cur["s_phase"],
                   A_re=A.real.tolist(), A_im=A.imag.tolist())
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        print(f"[n={n} t={t} s={seed}] init: ReC={best_r:+.6f} K={cur['K']:.5f}", flush=True)
        for it in range(iters):
            E = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
            E = E / np.linalg.norm(E, 2)
            A2 = A + step * E
            r2 = eval_A(A2, t, seed + 1000 + it)
            if r2 is not None and r2["C"].real < best_r:
                A, best_r, cur = A2, r2["C"].real, r2
                step = min(step * 1.25, 0.4)
                rec = dict(it=it, ReC=best_r, K=cur["K"], s_phase=cur["s_phase"],
                           A_re=A.real.tolist(), A_im=A.imag.tolist())
                fh.write(json.dumps(rec) + "\n"); fh.flush()
                print(f"[n={n} t={t} s={seed}] it={it}: ReC={best_r:+.6f} "
                      f"K={cur['K']:.5f} s_phase={cur['s_phase']:+.5f}", flush=True)
            else:
                step = max(step * 0.85, 0.02)
    print(f"DONE n={n} t={t} s={seed}: min ReC = {best_r:+.8f} "
          f"(K={cur['K']:.5f}, s_phase={cur['s_phase']:+.5f})", flush=True)

if __name__ == "__main__":
    main()
