"""Near-disk regime probe: Omega_eps = psi_eps(D), psi = w + eps w^2.
A = tau*I + s*(Crabb(m) + t*E): scaled to nearly fill Omega (max s with DLP margin).
At the true Blaschke extremal: K, c, ReC, G, L7@ext slack, and K + c/2 (<= 2?).

Also fits: how do (2-K) and c scale with eps?
"""

import sys

import numpy as np

from crouzeix import crabb_matrix
from extremal_pullback import Pullback, find_extremal_blaschke


def max_fill(pb, base, tau_grid, margin=2e-3):
    """Maximize s over tau in tau_grid such that W(tau I + s*base) inside Omega."""
    best = (None, 0.0, None)
    n = base.shape[0]
    for tau in tau_grid:
        lo, hi = 0.0, 3.0
        for _ in range(40):
            s = (lo + hi) / 2
            A = tau * np.eye(n) + s * base
            lmin, _ = pb.dlp_certificate(A)
            if lmin > margin:
                lo = s
            else:
                hi = s
        if lo > best[1]:
            A = tau * np.eye(n) + lo * base
            best = (tau, lo, A)
    return best


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    tpert = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    rng = np.random.default_rng(seed)
    print(f"# seed={seed} perturbation t={tpert}")
    for m in [2, 3]:
        n = m + 1
        E = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        E = E / np.linalg.norm(E, 2)
        base = crabb_matrix(m) + tpert * E
        for eps in [0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25]:
            def psi(w, e=eps):
                return w + e * w * w

            def dpsi(w, e=eps):
                return 1 + 2 * e * w

            pb = Pullback(psi, dpsi, N=1024)
            tau, s, A = max_fill(pb, base, np.linspace(-0.3, 0.3, 13))
            if A is None:
                continue
            al, K1 = find_extremal_blaschke(pb, A, m, restarts=24, seed=7 * m + seed)
            d = pb.extremal_data(A, al)
            print(
                f"m={m} eps={eps:.2f} tau={tau:+.3f} s={s:.4f}: K={d['K']:.6f} "
                f"c={d['c']:.6f} ReC={d['C'].real:+.6f} ImC={d['C'].imag:+.6f} "
                f"G={d['G']:.6f} |b|={abs(d['beta']):.6f} diag={d['diag']:.2e} "
                f"L7slack={d['slack']:+.6f} K+c/2={d['K'] + d['c'] / 2:.6f}",
                flush=True,
            )


if __name__ == "__main__":
    main()
