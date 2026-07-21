"""Epoch 3: maximize rho = c / (2(2-K)) over matrices at TRUE extremal pairs.
Also track X = 2ReC + G^2 - |beta|^2 and the exact L7@ext slack.

Outer loop: adversarial ascent over matrix entries A (n=3 default) inside a fixed
domain Omega_eps = psi(D), psi = w + eps w^2. Inner loop: exact extremal pair via
Blaschke-zero optimization (deg n-1). A is scale-adjusted to keep W(A) inside Omega
with margin. Certificates: DLP positivity (in scale_into), diag ~ 0 at extremal.

Usage: ratio_adversarial.py eps n outer_iters seed
"""
import json
import sys

import numpy as np

from extremal_pullback import Pullback, find_extremal_blaschke, scale_into

def eval_A(pb, A0, deg, seed, restarts=16):
    A, lmin = scale_into(pb, A0, margin=1e-3)
    if A is None:
        return None
    al, _ = find_extremal_blaschke(pb, A, deg, restarts=restarts, seed=seed)
    d = pb.extremal_data(A, al)
    if d["K"] <= 1.005 or d["diag"] > 1e-5:
        return None
    X = 2 * d["C"].real + d["G"] ** 2 - abs(d["beta"]) ** 2
    rho = d["c"] / (2 * (2 - d["K"])) if d["K"] < 2 else np.inf
    return dict(A=A, alphas=al, K=d["K"], c=d["c"], X=X, rho=rho,
                slack=d["slack"], diag=d["diag"])

def main():
    eps = float(sys.argv[1]); n = int(sys.argv[2])
    iters = int(sys.argv[3]); seed = int(sys.argv[4])
    rng = np.random.default_rng(seed)
    psi = lambda w, e=eps: w + e * w * w
    dpsi = lambda w, e=eps: 1 + 2 * e * w
    pb = Pullback(psi, dpsi, N=768)

    from crouzeix import crabb_matrix
    A_cur = crabb_matrix(n - 1) + 0.2 * (rng.standard_normal((n, n))
                                         + 1j * rng.standard_normal((n, n)))
    r_cur = eval_A(pb, A_cur, n - 1, seed)
    while r_cur is None:
        A_cur = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        r_cur = eval_A(pb, A_cur, n - 1, seed)
    step = 0.25
    best = r_cur
    fname = f"ratio_adv_eps{eps}_n{n}_s{seed}.jsonl"
    with open(fname, "w") as fh:
        for it in range(iters):
            E = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
            E = E / np.linalg.norm(E, 2)
            r_new = eval_A(pb, r_cur and (best["A"] / np.linalg.norm(best["A"], 2) + step * E), n - 1, seed + it)
            if r_new is not None and r_new["rho"] > best["rho"]:
                best = r_new
                step = min(step * 1.3, 0.5)
                rec = {k: (v.tolist() if isinstance(v, np.ndarray) else
                           (list(map(complex, v)) and [[a.real, a.imag] for a in v]
                            if isinstance(v, list) else v))
                       for k, v in best.items() if k != "A"}
                rec["A_re"] = best["A"].real.tolist()
                rec["A_im"] = best["A"].imag.tolist()
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                print(f"[eps={eps} n={n} s={seed}] it={it}: rho={best['rho']:.5f} "
                      f"K={best['K']:.5f} c={best['c']:.5f} X={best['X']:+.5f} "
                      f"slack={best['slack']:+.5f}", flush=True)
            else:
                step = max(step * 0.85, 0.02)
    print(f"DONE eps={eps} n={n} s={seed}: max rho={best['rho']:.6f} "
          f"(K={best['K']:.5f} c={best['c']:.5f} slack={best['slack']:+.5f})", flush=True)

if __name__ == "__main__":
    main()
