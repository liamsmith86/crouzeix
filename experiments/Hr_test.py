"""Test H-r: ReC >= 0 at extremal pairs when Omega -> int W(A), for n = 3.

For each matrix A (n=3, generic nonnormal, smooth strictly-convex W(A)):
  - map Omega_t = W(A) offset by t for t in {0.05, 0.02, 0.01, 0.005, 0.002}
  - Theodorsen conformal map, GeneralPullback, exact Blaschke extremal (deg 2)
  - record ReC(t), s_phase(t), K(t): does ReC increase as t -> 0? sign at smallest t?
Certificates: DLP positivity >= -1e-9 (skip if fails), mass err, unitality, diag.
"""
import numpy as np

from extremal_pullback import find_extremal_blaschke
from theodorsen import theodorsen_map, GeneralPullback

def run_matrix(A, label, ts=(0.05, 0.02, 0.01, 0.005, 0.002), N=1024, seed=0):
    print(f"--- {label}")
    for t in ts:
        try:
            z, zp, terr = theodorsen_map(A, N=N, inflate=t)
        except Exception as e:
            print(f"  t={t}: theodorsen failed: {e}")
            continue
        pb = GeneralPullback(z, zp)
        uerr = pb.unitality_certificate()
        lmin, merr = pb.dlp_certificate(A)
        if lmin < -1e-9 or merr > 1e-7 or uerr > 1e-9:
            print(f"  t={t}: certificates failed (lmin={lmin:.1e} mass={merr:.1e} unit={uerr:.1e})")
            continue
        al, _ = find_extremal_blaschke(pb, A, A.shape[0] - 1, restarts=20, seed=seed)
        d = pb.extremal_data(A, al)
        print(f"  t={t:5.3f}: K={d['K']:.6f} ReC={d['C'].real:+.6f} ImC={d['C'].imag:+.6f} "
              f"q={d['q']:.6f} s_phase={d['s_phase']:+.6f} diag={d['diag']:.1e} "
              f"[lmin={lmin:.1e} theod={terr:.1e}]", flush=True)

if __name__ == "__main__":
    rng = np.random.default_rng(17)
    from crouzeix import crabb_matrix

    # matrices that produced ReC < 0 in psi-domain probes: perturbed Crabb
    E = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    E = E / np.linalg.norm(E, 2)
    run_matrix(crabb_matrix(2) + 0.3 * E, "crabb2+0.3E (seed17)")

    E2 = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    E2 = E2 / np.linalg.norm(E2, 2)
    run_matrix(crabb_matrix(2) + 0.5 * E2, "crabb2+0.5E2")

    A = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    run_matrix(A, "random 3x3")
