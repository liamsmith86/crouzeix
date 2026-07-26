"""Sweep the doubly-symmetric 3x3 family; record closed-form data per matrix.

For (a,b,c,d) grid/random: e, kappa, tau, alpha*, z1, psi'(alpha), K, rho,
g0(0), g0(e). Output jsonl for landscape analysis (rho=0 locus, positivity margins,
correlation of extremal zeros with eigenvalue-Blaschke critical points).

Usage: sym3_sweep.py count seed
"""

import json
import sys

import numpy as np

from crouzeix import nr_support
from theodorsen import theodorsen_map, GeneralPullback
from minkowski_test import best_extremal
from zero_geometry import phi_of_points


def one(a, b, c, dd, seed):
    A = np.array([[0, a, 0], [b, 0, c], [0, dd, 0]], dtype=complex)
    e2 = a * b + c * dd
    if e2 <= 1e-6:
        return None
    e = float(np.sqrt(e2))
    th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    hs, _ = nr_support(A, th)
    if (hs + hs[np.arange(64) - 32]).min() < 0.25:
        return None
    chi = np.array([c, 0, -b], dtype=complex)
    rl = np.array([dd, 0, -a], dtype=complex)
    den = abs(rl.conj() @ chi)
    if den < 1e-9:
        return None
    kappa = float(np.linalg.norm(chi) * np.linalg.norm(rl) / den)
    try:
        z, zp, terr = theodorsen_map(A, N=1024, inflate=0.005)
    except Exception:
        return None
    pb = GeneralPullback(z, zp)
    lmin, merr = pb.dlp_certificate(A)
    if lmin < -1e-9 or merr > 1e-6:
        return None
    al, d = best_extremal(pb, A, 2, seed)
    if d["diag"] > 1e-6 or d["K"] < 1.02:
        return None
    tau = float(phi_of_points(pb, [e])[0].real)
    alpha = float(abs(al[0]))
    dw = 1j * pb.w * (2 * np.pi / pb.N)
    z1 = float((np.sum(pb.z / (pb.w - alpha) * dw) / (2j * np.pi)).real)
    dpsi = float((np.sum(pb.z / (pb.w - alpha) ** 2 * dw) / (2j * np.pi)).real)
    v = (tau**2 - alpha**2) / (1 - alpha**2 * tau**2)
    Bp = 2 * alpha / (1 - alpha**4)
    r1 = dpsi / Bp
    g0v = -1.0 / alpha**2 + 2 * r1 / z1 if alpha > 1e-9 else 0.0
    gev = 1.0 / v - 2 * r1 * z1 / (e**2 - z1**2)
    return dict(
        a=a,
        b=b,
        c=c,
        d=dd,
        e=e,
        kappa=kappa,
        tau=tau,
        alpha=alpha,
        z1=z1,
        dpsi=dpsi,
        K=d["K"],
        rho=d["C"].real,
        g0=g0v,
        ge=gev,
        diag=d["diag"],
        zeros=[[x.real, x.imag] for x in al],
    )


def main():
    count = int(sys.argv[1])
    seed = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    with open(f"sym3_sweep_s{seed}.jsonl", "w") as fh:
        got = 0
        tries = 0
        while got < count and tries < 20 * count:
            tries += 1
            a = rng.uniform(0.8, 2.2)
            c = rng.uniform(0.5, 2.0)
            b = rng.uniform(-0.6, 0.6)
            dd = rng.uniform(-0.6, 0.6)
            r = one(a, b, c, dd, seed + tries)
            if r is None:
                continue
            got += 1
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(
                f"[{got}/{count}] K={r['K']:.4f} rho={r['rho']:+.6f} kappa={r['kappa']:.2f} "
                f"tau={r['tau']:.3f} alpha={r['alpha']:.3f} g0={r['g0']:+.4f} ge={r['ge']:+.4f}",
                flush=True,
            )
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
