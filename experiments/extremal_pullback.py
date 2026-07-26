"""Exact extremal-pair machinery on domains Omega = psi(D), psi explicit univalent.

Key facts used:
- (Crouzeix 04) extremal f0 = B o psi^{-1}, B Blaschke deg <= n-1. So optimize over
  Blaschke zeros only -> TRUE extremal pair.
- Pullback quadrature: for h in A(Omega),
    h(A) = (1/2pi i) oint_{|w|=1} h(psi(w)) (psi(w) - A)^{-1} psi'(w) dw.
- g0 = Phi(f0) explicitly: on the boundary conj(f0) = 1/f0 (unimodular), so
    g0(z) = 1/f0(z) - PP(z),  PP(z) = sum_i r_i/(z - z_i),
  z_i = psi(alpha_i) zeros of f0 (alpha_i = Blaschke zeros, assumed simple),
  r_i = psi'(alpha_i)/B'(alpha_i).
- Certificate: P(sigma) >= 0 along boundary (W(A) inside Omega), Phi unitality via
  the same PP mechanism is exact, mass = 2I check.

Diagnostics at extremal pair: K, c=|C|, C, G=||g0(A)*x0||, W, beta, L7@ext.
"""

import numpy as np
from scipy.optimize import minimize


def blaschke(w, alphas):
    out = np.ones_like(w, dtype=complex)
    for a in alphas:
        out *= (w - a) / (1 - np.conj(a) * w)
    return out


def blaschke_prime_at_zero(a_list, i):
    """B'(alpha_i) for Blaschke with simple zeros a_list:
    b_i'(alpha_i) * prod_{j != i} b_j(alpha_i), b_i'(alpha_i) = 1/(1-|alpha_i|^2)."""
    ai = a_list[i]
    out = 1.0 / (1 - abs(ai) ** 2)
    for j, a in enumerate(a_list):
        if j != i:
            out *= (ai - a) / (1 - np.conj(a) * ai)
    return out


class Pullback:
    def __init__(self, psi, dpsi, N=2048):
        self.psi, self.dpsi, self.N = psi, dpsi, N
        th = np.linspace(0, 2 * np.pi, N, endpoint=False)
        self.w = np.exp(1j * th)
        self.z = psi(self.w)  # boundary of Omega
        self.dz = dpsi(self.w) * 1j * self.w  # dz/dtheta

    def resolvent_stack(self, A):
        """Precompute T_j = dz_j (z_j - A)^{-1} dtheta/(2 pi i); h(A) = sum h_j T_j."""
        n = A.shape[0]
        identity = np.eye(n, dtype=complex)
        T = np.empty((self.N, n, n), dtype=complex)
        for j in range(self.N):
            T[j] = (self.dz[j] * (2 * np.pi / self.N) / (2j * np.pi)) * np.linalg.inv(
                self.z[j] * identity - A
            )
        return T

    def calc(self, hvals, A, T=None):
        """h(A) from boundary values hvals of h in A(Omega) (on self.z nodes)."""
        if T is None:
            T = self.resolvent_stack(A)
        return np.einsum("j,jkl->kl", hvals, T)

    def dlp_certificate(self, A):
        """min over nodes of lambda_min(P), and ||mass - 2I||."""
        n = A.shape[0]
        identity = np.eye(n, dtype=complex)
        S = np.zeros((n, n), dtype=complex)
        lmin = np.inf
        for j in range(self.N):
            B = (self.dz[j] / (2j * np.pi)) * np.linalg.inv(
                self.z[j] * identity - A
            )
            Q = B + B.conj().T
            lmin = min(lmin, np.linalg.eigvalsh(Q)[0])
            S += Q
        S *= 2 * np.pi / self.N
        return lmin, np.linalg.norm(S - 2 * identity, 2)

    def extremal_data(self, A, alphas, T=None):
        """Full diagnostic at the Blaschke extremal with given zeros."""
        if T is None:
            T = self.resolvent_stack(A)
        alphas = list(alphas)
        # residue formula assumes simple zeros; nudge near-coincident ones
        for i in range(len(alphas)):
            for j in range(i):
                if abs(alphas[i] - alphas[j]) < 1e-5:
                    alphas[i] += 3e-5 * np.exp(2j * np.pi * (i + 1) / 7)
        Bv = blaschke(self.w, alphas)  # f0 on boundary (pullback)
        f0A = self.calc(Bv, A, T)
        K = np.linalg.norm(f0A, 2)
        U, s, Vh = np.linalg.svd(f0A)
        x0 = Vh[0].conj()
        u0 = U[:, 0]
        diag = abs(np.vdot(x0, f0A @ x0))
        # g0 explicitly. Note r_i = residue of 1/f0 at z_i in the z-plane:
        # 1/f0(z) ~ 1/(f0'(z_i)(z-z_i)), f0'(z_i) = B'(alpha_i)/psi'(alpha_i).
        zi = self.psi(np.array(alphas, dtype=complex))
        ri = np.array(
            [
                self.dpsi(a) / blaschke_prime_at_zero(alphas, i)
                for i, a in enumerate(alphas)
            ],
            dtype=complex,
        )
        g0v = 1.0 / Bv - sum(r / (self.z - z0) for r, z0 in zip(ri, zi))
        g0A = self.calc(g0v, A, T)
        C = np.vdot(x0, g0A @ (f0A @ x0))  # <g0(A)f0(A)x0, x0>... careful order
        # np.vdot(a,b) = conj(a).b -> C = <(g0 f0)(A) x0, x0>? We want <Xx0,x0> = x0* X x0
        C = x0.conj() @ (g0A @ f0A @ x0)
        c_val = abs(C)
        Gv = g0A.conj().T @ x0
        G = np.linalg.norm(Gv)
        v = f0A @ x0 + Gv
        W = np.linalg.norm(v)
        beta = x0.conj() @ v
        lhs = np.sqrt(max(W**2 - abs(beta) ** 2, 0))
        rhs = 2 - c_val / 2
        return dict(
            K=K,
            c=c_val,
            C=C,
            G=G,
            W=W,
            beta=beta,
            diag=diag,
            lhs=lhs,
            rhs=rhs,
            slack=rhs - lhs,
            x0=x0,
            u0=u0,
        )


def find_extremal_blaschke(pb, A, deg, restarts=40, seed=0):
    """Maximize ||(B o phi)(A)|| over Blaschke zeros in D^deg."""
    rng = np.random.default_rng(seed)
    T = pb.resolvent_stack(A)

    def unpack(x):
        pts = x[:deg] + 1j * x[deg:]
        # map C -> D smoothly
        return [p / np.sqrt(1 + abs(p) ** 2) for p in pts]

    def obj(x):
        al = unpack(x)
        Bv = blaschke(pb.w, al)
        return -np.linalg.norm(pb.calc(Bv, A, T), 2)

    best, bx = 0.0, None
    for t in range(restarts):
        x0 = rng.standard_normal(2 * deg) * (0.7 if t % 2 else 1.5)
        res = minimize(
            obj,
            x0,
            method="Nelder-Mead",
            options={"maxiter": 6000, "fatol": 1e-13, "xatol": 1e-11},
        )
        if -res.fun > best:
            best, bx = -res.fun, res.x
    return unpack(bx), best


def scale_into(pb, A, margin=1e-3):
    """Scale A so that W(A) fits inside Omega with positive DLP margin."""
    A = A / np.linalg.norm(A, 2)
    s = 1.0
    for _ in range(60):
        lmin, _ = pb.dlp_certificate(s * A)
        if lmin > margin:
            return s * A, lmin
        s *= 0.85
    return None, None


if __name__ == "__main__":
    rng = np.random.default_rng(5)
    for eps in [0.0, 0.05, 0.1, 0.2]:
        def psi(w, e=eps):
            return w + e * w * w

        def dpsi(w, e=eps):
            return 1 + 2 * e * w

        pb = Pullback(psi, dpsi, N=1024)
        for trial in range(4):
            n = 3 if trial % 2 else 2
            A0 = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
            A, lmin = scale_into(pb, A0)
            if A is None:
                continue
            deg = n - 1
            al, K = find_extremal_blaschke(pb, A, deg, restarts=30, seed=trial)
            d = pb.extremal_data(A, al)
            print(
                f"eps={eps} n={n} t={trial}: K={d['K']:.6f} c={d['c']:.6f} "
                f"ReC={d['C'].real:+.6f} G={d['G']:.6f} W={d['W']:.6f} "
                f"|b|={abs(d['beta']):.6f} diag={d['diag']:.2e} "
                f"L7slack={d['slack']:+.6f} K+c/2={d['K'] + d['c'] / 2:.6f}",
                flush=True,
            )
