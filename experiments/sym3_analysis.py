"""Post-sweep analysis: verify across sym3_sweep records
  (1) rho three-factor formula residual,
  (2) pi0 := <P0 x0, x0> (from closed-form frame) — is Re(pi0) ~ 1/2? >= 0 always?
  (3) g_e >= g_0 always?
  (4) reconstruct rho from factors; check H-r margins.
"""
import glob
import json

import numpy as np

recs = []
for fn in glob.glob("sym3_sweep_s*.jsonl"):
    for line in open(fn):
        recs.append(json.loads(line))
print(f"{len(recs)} records")

bad = 0
for r in recs:
    alpha, tau, kappa = r["alpha"], r["tau"], r["kappa"]
    if alpha < 1e-3:
        continue
    v = (tau ** 2 - alpha ** 2) / (1 - alpha ** 2 * tau ** 2)
    m = alpha ** 2 + v
    sk = np.sqrt(kappa ** 2 - 1)
    M = np.array([[-alpha ** 2, -m * sk], [0, v]])
    U, s, Vh = np.linalg.svd(M)
    x0 = Vh[0]
    E = np.array([[1.0, sk], [0, 0]])
    pi0 = x0.conj() @ (E @ x0)
    rho_formula = alpha ** 2 * (r["ge"] - r["g0"]) * pi0.real
    ok_g = r["ge"] >= r["g0"] - 1e-12
    resid = abs(rho_formula - r["rho"])
    if not ok_g or pi0.real < 0 or resid > 5e-3 * max(abs(r["rho"]), 1e-6):
        bad += 1
        print(f"ANOMALY: K={r['K']:.4f} kappa={kappa:.2f} rho={r['rho']:+.2e} "
              f"formula={rho_formula:+.2e} pi0={pi0.real:+.4f} ge-g0={r['ge']-r['g0']:+.2e}")
    else:
        print(f"K={r['K']:.4f} kap={kappa:6.2f} rho={r['rho']:+.2e} "
              f"pi0={pi0.real:+.5f} ge-g0={r['ge'] - r['g0']:+.5f} resid={resid:.1e}")
print(f"\nanomalies: {bad}/{len(recs)}")
