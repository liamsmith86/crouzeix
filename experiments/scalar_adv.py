"""Adversarial maximization of L10 violation (scalar problem).

Configuration: ellipse (1, b); w = sum_j [ a_j/(2pi i (z-z_j)) + b_j/(2pi i (z-z_j)^2) ]
with sum a_j = 1 (complex allowed); positivity dnu = 2Re(w dsigma) >= 0 enforced
by penalty + final rejection. p = poly deg d, constraint oint p w = 0 by projection.

viol = LHS - RHS,  LHS = 2 oint |p|^2 dnu - |oint conj(p) w dsigma|^2,
RHS = (2 - c/2)^2, c = |oint Phi(p) p w dsigma| (= residue sum).

Usage: scalar_adv.py b npts deg restarts seed [second_order]
"""
import json
import sys

import numpy as np
from scipy.optimize import minimize

from L10_test import ellipse, poly_z, phi_at

def eval_config(b_ell, zs, ams, bms, cp, N=1024, want_parts=False):
    z, zp = ellipse(1.0, b_ell, N)
    dt = 2 * np.pi / N
    w = np.zeros(N, dtype=complex)
    for zj, aj, bj in zip(zs, ams, bms):
        w += aj / (2j * np.pi * (z - zj)) + bj / (2j * np.pi * (z - zj) ** 2)
    nu = 2 * np.real(w * zp) * dt
    numin = nu.min()
    pv = poly_z(cp, z)
    m = np.max(np.abs(pv))
    if m < 1e-13 or not np.isfinite(m):
        return None
    cpn = [ci / m for ci in cp]
    pv = pv / m
    # oint p w and oint conj(p) w via residues / quadrature
    zsa = np.array(zs)
    pz = poly_z(cpn, zsa)
    dpz = np.array([np.polyval(np.polyder(np.array(cpn[::-1])), zj) for zj in zsa])
    con = abs(np.sum(np.array(ams) * pz) + np.sum(np.array(bms) * dpz))
    beta_bar = np.sum(np.conj(pv) * w * zp) * dt
    lhs = 2 * np.sum(np.abs(pv) ** 2 * nu) - abs(beta_bar) ** 2
    g_at = phi_at(pv, z, zp, zsa)
    # residue of g*p*w: first-order part a_j g(z_j)p(z_j); second-order part
    # b_j * d/dz[g p](z_j). g' at z_j via quadrature.
    gp_at = np.array([np.sum(np.conj(pv) * zp / (z - zj) ** 2) for zj in zsa]) \
        * dt / (2j * np.pi)
    d_gp = gp_at * pz + g_at * dpz
    c_val = abs(np.sum(np.array(ams) * g_at * pz) + np.sum(np.array(bms) * d_gp))
    rhs = (2 - c_val / 2) ** 2
    out = dict(lhs=lhs, rhs=rhs, c=c_val, con=con, viol=lhs - rhs, numin=numin,
               mass=float(np.sum(nu)))
    return out

def main():
    b_ell = float(sys.argv[1]); npts = int(sys.argv[2]); deg = int(sys.argv[3])
    restarts = int(sys.argv[4]); seed = int(sys.argv[5])
    second = len(sys.argv) > 6 and sys.argv[6] == "1"
    rng = np.random.default_rng(seed)
    NP = npts

    def unpack(x):
        o = 0
        zs = x[o:o + NP] * 0.9 + 1j * x[o + NP:o + 2 * NP] * 0.9 * b_ell; o += 2 * NP
        am = x[o:o + NP] + 1j * x[o + NP:o + 2 * NP]; o += 2 * NP
        s = np.sum(am)
        if abs(s) < 1e-10:
            return None
        am = am / s
        if second:
            bm = x[o:o + NP] + 1j * x[o + NP:o + 2 * NP]; o += 2 * NP
        else:
            bm = np.zeros(NP, dtype=complex); o += 0
        cp = list(x[o:o + deg + 1] + 1j * x[o + deg + 1:o + 2 * (deg + 1)])
        # points must be strictly inside ellipse
        u = np.real(zs) ** 2 + (np.imag(zs) / b_ell) ** 2
        if np.any(u > 0.92 ** 2):
            return None
        return zs, am, bm, cp

    nvar = (6 if second else 4) * NP + 2 * (deg + 1)

    def obj(x):
        u = unpack(x)
        if u is None:
            return 5.0
        zs, am, bm, cp = u
        r = eval_config(b_ell, zs, am, bm, cp, N=512)
        if r is None:
            return 5.0
        pen = 3e3 * max(0.0, -r["numin"]) + 30.0 * r["con"] ** 2
        return -(r["viol"]) + pen

    tag = f"b{b_ell}_k{npts}_d{deg}_s{seed}{'_2nd' if second else ''}"
    best = -99
    with open(f"scalar_adv_{tag}.jsonl", "w") as fh:
        for t in range(restarts):
            x0 = rng.standard_normal(nvar) * 0.5
            res = minimize(obj, x0, method="Nelder-Mead",
                           options={"maxiter": 4000, "fatol": 1e-10, "xatol": 1e-9})
            u = unpack(res.x)
            if u is None:
                continue
            zs, am, bm, cp = u
            r = eval_config(b_ell, zs, am, bm, cp, N=8192)
            if r is None or r["numin"] < -1e-9 or r["con"] > 1e-6:
                continue
            rec = dict(trial=t, **{k: (v if not isinstance(v, complex) else [v.real, v.imag])
                                   for k, v in r.items()},
                       zs=[[z.real, z.imag] for z in zs],
                       am=[[a.real, a.imag] for a in am],
                       bm=[[a.real, a.imag] for a in bm],
                       cp=[[a.real, a.imag] for a in cp])
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if r["viol"] > best:
                best = r["viol"]
                print(f"[{tag}] t={t}: viol={r['viol']:+.6f} lhs={r['lhs']:.5f} "
                      f"rhs={r['rhs']:.5f} c={r['c']:.5f} numin={r['numin']:.2e}", flush=True)
    print(f"DONE {tag}: max certified viol = {best:+.8f}", flush=True)

if __name__ == "__main__":
    main()
