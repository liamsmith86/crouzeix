"""Test L7: sqrt(W^2 - |beta|^2) <= 2 - c/2 at near-extremal configurations.

Quantities (all x0-localized, computable without the extremal measure):
  p0 : near-extremal polynomial, ||p0||_{dOmega} = 1, K = ||p0(M)||
  x0 : top right singular vector of p0(M)
  g0 = Phi(p0)  (conjugate Cauchy transform)
  W    = || p0(M)x0 + g0(M)* x0 ||
  beta = < (p0(M) + g0(M)*) x0, x0 >
  c    = | < (g0*p0)(M) x0, x0 > |
  diag = | < p0(M)x0, x0 > |   (extremality certificate; -> 0 at true extremal)
L6 check: K^2 <= K*sqrt(W^2-|beta|^2) + c + |beta|*diag  (algebraic, must hold).
L7 test: sqrt(W^2-|beta|^2) <= 2 - c/2 ?
"""
import json
import sys

import numpy as np
from scipy.optimize import minimize

from crouzeix import poly_z, poly_A, opnorm, crabb_matrix
from star_inequality import offset_curve, cauchy_matrix, phi_on_inner, dlp_min_eig_integral

def find_extremal_poly(M, eps, deg, rng, restarts=60, N=256, structured=True):
    _, z_out, _ = offset_curve(M, eps, N)

    def neg_ratio(x):
        c = x[:deg + 1] + 1j * x[deg + 1:]
        m = np.max(np.abs(poly_z(c, z_out)))
        if m < 1e-14:
            return 0.0
        return -opnorm(poly_A(list(c), M)) / m

    inits = []
    if structured:
        for k in range(1, deg + 1):
            e = np.zeros(2 * (deg + 1))
            e[k] = 1.0
            inits.append(e)
    while len(inits) < restarts:
        inits.append(rng.standard_normal(2 * (deg + 1)))
    best, bestc = 0.0, None
    for x0 in inits:
        res = minimize(neg_ratio, x0, method="L-BFGS-B",
                       options={"maxiter": 500, "eps": 1e-8})
        if -res.fun > best:
            best = -res.fun
            bestc = res.x[:deg + 1] + 1j * res.x[deg + 1:]
    return list(bestc), best

def refined_quantities(M, c_poly, eps=0.05, N=2048):
    _, z_out, zp_out = offset_curve(M, eps, N)
    _, z_in, zp_in = offset_curve(M, eps / 2, N)
    pv_out = poly_z(c_poly, z_out)
    m = np.max(np.abs(pv_out))
    c_poly = [ci / m for ci in c_poly]
    pM = poly_A(c_poly, M)
    K = opnorm(pM)
    U, s, Vh = np.linalg.svd(pM)
    x0 = Vh[0].conj()
    diag = abs(np.vdot(x0, pM @ x0))
    phi_in = phi_on_inner(poly_z(c_poly, z_out), z_out, zp_out, z_in)
    pv_in = poly_z(c_poly, z_in)
    g0M = cauchy_matrix(phi_in, z_in, zp_in, M)
    prodM = cauchy_matrix(phi_in * pv_in, z_in, zp_in, M)
    v = pM @ x0 + g0M.conj().T @ x0
    W = np.linalg.norm(v)
    beta = np.vdot(x0, v)
    c_val = abs(np.vdot(x0, prodM @ x0))
    s2 = max(W**2 - abs(beta)**2, 0.0)
    lhs = np.sqrt(s2)
    rhs = 2 - c_val / 2
    L6_lhs = K**2
    L6_rhs = K * lhs + c_val + abs(beta) * diag
    L, _, _ = dlp_min_eig_integral(M, z_out, zp_out)
    return dict(K=K, W=W, beta=beta, c=c_val, diag=diag, lhs=lhs, rhs=rhs,
                L=L, L6_ok=L6_lhs <= L6_rhs + 1e-9, L6_gap=L6_rhs - L6_lhs)

def run_case(name, M, deg, eps, rng):
    c0, r = find_extremal_poly(M, eps, deg, rng)
    q = refined_quantities(M, c0, eps=eps)
    print(f"[{name}] K={q['K']:.6f} W={q['W']:.6f} |b|={abs(q['beta']):.6f} "
          f"c={q['c']:.6f} diag={q['diag']:.2e} | L7: {q['lhs']:.6f} <= {q['rhs']:.6f} "
          f"{'HOLDS' if q['lhs'] <= q['rhs'] + 1e-9 else 'FAILS'} "
          f"(slack {q['rhs']-q['lhs']:+.6f}) | L6gap={q['L6_gap']:+.3e} | 2L={2*q['L']:.4f}")
    return q

if __name__ == "__main__":
    rng = np.random.default_rng(7)

    run_case("crabb2 (disk)", crabb_matrix(2), 8, 0.03, rng)

    E = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    run_case("crabb2+0.2E", crabb_matrix(2) + 0.2 * E, 8, 0.03, rng)
    run_case("crabb2+0.4E", crabb_matrix(2) + 0.4 * E, 8, 0.03, rng)

    E4 = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
    run_case("crabb3+0.2E", crabb_matrix(3) + 0.2 * E4, 8, 0.03, rng)

    run_case("2x2 [[0,1],[0,1]]", np.array([[0, 1], [0, 1]], dtype=complex), 8, 0.03, rng)

    # best matrix from Track C search n=3 d=3
    try:
        best = None
        for line in open("results_n3_d3_s1.jsonl"):
            rec = json.loads(line)
            if best is None or rec["r_inner1024"] > best["r_inner1024"]:
                best = rec
        A = np.array(best["A_re"]) + 1j * np.array(best["A_im"])
        A = A / opnorm(A)
        run_case("trackC-n3-best", A, 8, 0.03, rng)
    except FileNotFoundError:
        pass

    for t in range(2):
        M = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        run_case(f"rand3x3-{t}", M, 8, 0.05, rng)
