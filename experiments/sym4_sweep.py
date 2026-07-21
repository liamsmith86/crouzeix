"""sym4 sweep: gather collapse diagnostics for the 4x4 zero-diag tridiagonal family.

Per matrix: weights, eigen-pairs (e1,e2), extremal zeros (verify {0,±alpha}),
K, rho, diag; block data: N1 = A-A+ (even), N2 = A+A- (odd), obliquities;
values f0(e_j) = B(tau_j), g0(e_j); sublattice weights q_j = <Q_j x0, x0>
(spectral projections of the block containing x0); test rho = Re sum B(tau_j) g0(e_j) q_j.

Usage: sym4_sweep.py count seed
"""
import json
import sys

import numpy as np

from crouzeix import nr_support
from theodorsen import theodorsen_map, GeneralPullback
from minkowski_test import best_extremal
from extremal_pullback import blaschke

def one(ws, seed):
    a1, b1, a2, b2, a3, b3 = ws
    A = np.zeros((4, 4), dtype=complex)
    A[0, 1], A[1, 0] = a1, b1
    A[1, 2], A[2, 1] = a2, b2
    A[2, 3], A[3, 2] = a3, b3
    eigs = np.linalg.eigvals(A)
    ers = np.sort(np.abs(eigs))
    # want two distinct REAL pairs for the main sector
    if np.max(np.abs(eigs.imag)) > 1e-9:
        return None
    e1, e2 = ers[3], ers[1]
    if e1 - e2 < 0.05 or e2 < 0.02:
        return None
    th64 = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    hs, _ = nr_support(A, th64)
    if (hs + hs[np.arange(64) - 32]).min() < 0.3:
        return None
    try:
        z, zp, terr = theodorsen_map(A, N=1024, inflate=0.005)
    except Exception:
        return None
    pb = GeneralPullback(z, zp)
    lmin, merr = pb.dlp_certificate(A)
    if lmin < -1e-9 or merr > 1e-6:
        return None
    al, d = best_extremal(pb, A, 3, seed)
    if d["diag"] > 1e-6 or d["K"] < 1.05:
        return None
    # verify odd-zero structure {0, ±alpha}: find the near-zero zero and the pair
    als = np.array(al)
    i0 = np.argmin(np.abs(als))
    rest = np.delete(als, i0)
    if abs(als[i0]) > 0.02 or abs(rest[0] + rest[1]) > 0.02:
        struct = "OTHER"
        alpha = np.nan
    else:
        struct = "ODD"
        alpha = complex(rest[0])
    from zero_geometry import phi_of_points
    taus = phi_of_points(pb, [e1, e2])
    # g0 at e1, e2 via extremal_data pieces: use pb.extremal_data internals — recompute g0 values
    T = pb.resolvent_stack(A)
    Bv = blaschke(pb.w, al)
    # g0 boundary values: 1/Bv - PP; reuse extremal_data by evaluating g0 at points via Cauchy:
    from extremal_pullback import blaschke_prime_at_zero
    alz = list(al)
    dw = 1j * pb.w * (2 * np.pi / pb.N)
    zi, ri = [], []
    for i, aa in enumerate(alz):
        zi.append(np.sum(pb.z / (pb.w - aa) * dw) / (2j * np.pi))
        dpsi = np.sum(pb.z / (pb.w - aa) ** 2 * dw) / (2j * np.pi)
        ri.append(dpsi / blaschke_prime_at_zero(alz, i))
    g0v = 1.0 / Bv - sum(r / (pb.z - z0) for r, z0 in zip(ri, zi))
    def g0_at(p):
        # g0 analytic in Omega: Cauchy from boundary values
        return np.sum(g0v * pb.dz / (pb.z - p)) * (2 * np.pi / pb.N) / (2j * np.pi)
    g0e = [g0_at(e1), g0_at(e2)]
    f0e = [complex(blaschke(np.array([t]), al)[0]) for t in taus]
    # sublattice weights: x0 lives in odd or even sublattice
    x0 = d["x0"]
    odd_mass = abs(x0[0]) ** 2 + abs(x0[2]) ** 2
    sub = "odd" if odd_mass > 0.99 else ("even" if odd_mass < 0.01 else "MIXED")
    # block containing x0 and its spectral projections
    if sub == "odd":
        N = np.array([[a1 * b1, a1 * a2], [b2 * b1, b2 * a2 + a3 * b3]], dtype=complex)
        xb = np.array([x0[0], x0[2]])
    else:
        N = np.array([[b1 * a1 + a2 * b2, a2 * a3], [b3 * b2, b3 * a3]], dtype=complex)
        xb = np.array([x0[1], x0[3]])
    lam, V = np.linalg.eig(N)
    Vi = np.linalg.inv(V)
    order = np.argsort(-np.abs(lam))
    qs = []
    for k in order:
        Qk = np.outer(V[:, k], Vi[k, :])
        qs.append(complex(xb.conj() @ (Qk @ xb)))
    rho_formula = sum((f0e[j] * g0e[j] * qs[j]).real for j in range(2))
    return dict(ws=list(ws), e1=e1, e2=e2, K=d["K"], rho=d["C"].real, diag=d["diag"],
                struct=struct, alpha=[alpha.real, alpha.imag] if struct == "ODD" else None,
                taus=[[t.real, t.imag] for t in taus], sub=sub,
                f0e=[[v.real, v.imag] for v in f0e], g0e=[[v.real, v.imag] for v in g0e],
                qs=[[q.real, q.imag] for q in qs], rho_formula=rho_formula)

def main():
    count = int(sys.argv[1]); seed = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    got = tries = 0
    with open(f"sym4_sweep_s{seed}.jsonl", "w") as fh:
        while got < count and tries < 40 * count:
            tries += 1
            ws = [rng.uniform(0.8, 2.2), rng.uniform(-0.4, 0.4),
                  rng.uniform(0.6, 2.0), rng.uniform(-0.4, 0.4),
                  rng.uniform(0.8, 2.2), rng.uniform(-0.4, 0.4)]
            r = one(ws, seed + tries)
            if r is None:
                continue
            got += 1
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            err = abs(r["rho_formula"] - r["rho"])
            print(f"[{got}] K={r['K']:.4f} rho={r['rho']:+.6f} formula={r['rho_formula']:+.6f} "
                  f"(err {err:.1e}) {r['struct']} sub={r['sub']} "
                  f"q1={r['qs'][0][0]:+.3f} q2={r['qs'][1][0]:+.3f}", flush=True)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
