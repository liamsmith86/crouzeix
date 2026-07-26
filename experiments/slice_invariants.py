"""Level-4 closure, structured attack: K and the stationarity law in frame invariants.

Active block M(alpha) = A+ F(N), N = A- A+ (2x2, eigs e1^2 > e2^2), F = F1 Q1 + F2 Q2,
F_j = B(tau_j)/e_j, B = odd Blaschke with zeros {0, +-alpha}.

Closed form:  det M = delta * F1 F2   (delta = det A+; det(F1Q1+F2Q2) = F1F2)
              T := ||M||_F^2 = p11 F1^2 + 2 p12 F1 F2 + p22 F2^2,  p_ij = tr(Q_i^T G Q_j),
              G = A+^T A+.
              K^2 = (T + S)/2,  S = sqrt(T^2 - 4 delta^2 F1^2 F2^2).

Stationarity dK^2/dalpha = 0 with dB/dalpha(w) = -2 a w (1-w^4)/(1-a^2 w^2)^2:
    sum_j (T_j + S_j) * Bp_j / e_j = 0,
    T_j = dT/dF_j, S_j = (T T_j - 4 delta^2 F_j F_other^2)/S, Bp_j = dB/dalpha(tau_j).

rho = B1 g1 q1 + B2 g2 q2, q_j = x0^T Q_j x0 (x0 = top right singular vec; q1+q2=1).

Verifies: (V1) K^2 closed form vs direct;  (V2) stationarity law residual at optimizer alpha*;
          (V3) sign table for rho terms;   (V4) Landen-type formula comparisons.
"""

import numpy as np
from mpmath import mp, sqrt as msqrt, pi as mpi, ellipk
from slice_exact import Slice

mp.dps = 25


def invariants(S):
    Ap, Am = S.blocks()
    N = Am @ Ap
    lam, V = np.linalg.eig(N)
    Vi = np.linalg.inv(V)
    E = [float(S.e1**2), float(S.e2**2)]
    # order projections: Q1 <-> e1^2 (outer), Q2 <-> e2^2
    idx1 = 0 if abs(lam[0] - E[0]) < abs(lam[1] - E[0]) else 1
    Q = [np.outer(V[:, i], Vi[i, :]) for i in (idx1, 1 - idx1)]
    G = Ap.T @ Ap
    p = np.array(
        [[np.trace(Q[i].T @ G @ Q[j]) for j in range(2)] for i in range(2)]
    ).real
    delta = np.linalg.det(Ap)
    return Ap, Q, p, delta


def B_and_prime(w, a):
    B = w * (w * w - a * a) / (1 - a * a * w * w)
    Bp = -2 * a * w * (1 - w**4) / (1 - a * a * w * w) ** 2  # dB/dalpha
    return B, Bp


def closed_K2(F1, F2, p, delta):
    T = p[0, 0] * F1**2 + 2 * p[0, 1] * F1 * F2 + p[1, 1] * F2**2
    Ssq = T * T - 4 * delta**2 * F1**2 * F2**2
    return (T + np.sqrt(Ssq)) / 2, T, np.sqrt(Ssq)


def stationarity_residual(S, alpha, p, delta):
    t1, t2 = complex(S.tau[0]).real, complex(S.tau[1]).real
    e1, e2 = complex(S.e1).real, complex(S.e2).real
    B1, Bp1 = B_and_prime(t1, alpha)
    B2, Bp2 = B_and_prime(t2, alpha)
    F1, F2 = B1 / e1, B2 / e2
    K2, T, Ssqrt = closed_K2(F1, F2, p, delta)
    T1 = 2 * (p[0, 0] * F1 + p[0, 1] * F2)
    T2 = 2 * (p[1, 1] * F2 + p[0, 1] * F1)
    S1 = (T * T1 - 4 * delta**2 * F1 * F2**2) / Ssqrt
    S2 = (T * T2 - 4 * delta**2 * F2 * F1**2) / Ssqrt
    lhs = (T1 + S1) * Bp1 / e1 + (T2 + S2) * Bp2 / e2
    scale = abs((T1 + S1) * Bp1 / e1) + abs((T2 + S2) * Bp2 / e2)
    return lhs / scale, K2, (F1, F2, B1, B2)


def rho_terms(S, alpha):
    """Return (B_j g_j q_j) pieces using slice_exact's exact g-values."""
    Kv, rho, g1, g2, diag = S.rho_at(alpha)
    Ap, Q, p, delta = invariants(S)
    t1, t2 = complex(S.tau[0]).real, complex(S.tau[1]).real
    B1, _ = B_and_prime(t1, alpha)
    B2, _ = B_and_prime(t2, alpha)
    F1, F2 = B1 / complex(S.e1).real, B2 / complex(S.e2).real
    F = F1 * Q[0] + F2 * Q[1]
    M = Ap @ F
    U, sv, Vh = np.linalg.svd(M)
    x0 = Vh[0].conj()
    q = [float(np.real(x0.conj() @ Q[j] @ x0)) for j in range(2)]
    terms = [B1 * g1.real * q[0], B2 * g2.real * q[1]]
    return Kv, rho, q, (B1, B2), (g1.real, g2.real), terms, p, delta


if __name__ == "__main__":
    cases = [
        (2.0, 1.2, 1.6, 0.08),
        (1.8, 1.0, 1.4, 0.15),
        (1.5, 2.0, 0.9, 0.12),
        (2.2, 0.9, 1.3, 0.10),
        (1.2, 1.7, 2.1, 0.06),
    ]
    for a1, a2, a3, c in cases:
        S = Slice(a1, a2, a3, c)
        al, Kv = S.extremal()
        Ap, Q, p, delta = invariants(S)
        res, K2c, FB = stationarity_residual(S, al, p, delta)
        Kdirect = S.eval_alpha(al)[0]
        Kv2, rho, q, Bs, gs, terms, _, _ = rho_terms(S, al)
        print(f"a=({a1},{a2},{a3}) c={c}: alpha*={al:.10f}")
        print(
            f"  V1 closed-K: {np.sqrt(K2c):.12f} vs direct {Kdirect:.12f}  "
            f"diff={abs(np.sqrt(K2c) - Kdirect):.1e}"
        )
        print(f"  V2 stationarity-law residual (normalized): {res:+.2e}")
        print(
            f"  V3 rho={rho:+.6e}  q=({q[0]:+.6f},{q[1]:+.6f})  "
            f"B=({Bs[0]:+.6f},{Bs[1]:+.6f})  g=({gs[0]:+.6f},{gs[1]:+.6f})"
        )
        print(f"     terms: B1*g1*q1={terms[0]:+.6e}  B2*g2*q2={terms[1]:+.6e}")
        # V4: Landen-type comparisons
        k = float(S.k)
        kp = float(msqrt(1 - S.m))
        k1 = (1 - kp) / (1 + kp)
        land1 = 1 - float(mpi / (2 * ellipk(k1 * k1)))
        land0 = 1 - float(mpi / (2 * ellipk(S.m)))
        print(
            f"  V4 rho={rho:+.6e} vs 1-pi/2K(k1^2)={land1:+.6e} vs 1-pi/2K(m)={land0:+.6e}"
        )
        print(
            f"     p11={p[0, 0]:.6f} p12={p[0, 1]:.6f} p22={p[1, 1]:.6f} delta={delta:.6f}"
        )
