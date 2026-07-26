"""Numerical Riemann map D -> int W(A) via Theodorsen's iteration, packaged as a
Pullback-compatible boundary parameterization.

int W(A) is convex (hence starlike w.r.t. any interior point). Boundary in polar
form r = rho(theta) about the centroid. Theodorsen: find boundary correspondence
theta(t) (t = angle on unit circle) with
    log rho(theta(t)) = Re log(psi(e^{it})/e^{it})-conjugate relation:
    theta(t) = t + H[log rho(theta(.))](t)   (H = periodic Hilbert transform),
iterated to fixed point; then psi(e^{it}) = rho(theta(t)) e^{i theta(t)} and psi
extends holomorphically with psi(0) = center, psi'(0) > 0.

We only need boundary values z(t) = psi(e^{it}) and dz/dt (FFT differentiation)
to drive the Pullback machinery (resolvent calculus, extremal Blaschke search).

Certificates: convex tangent winding, DLP positivity for the matrix, mass = 2I.
"""

import numpy as np

from crouzeix import nr_support


def hilbert_periodic(f):
    """Conjugate function (periodic Hilbert transform) via FFT: H[e^{ikt}] = -i sgn(k) e^{ikt}."""
    N = len(f)
    F = np.fft.fft(f)
    k = np.fft.fftfreq(N, d=1.0 / N)
    mult = -1j * np.sign(k)
    mult[0] = 0.0
    return np.real(np.fft.ifft(F * mult))


def boundary_polar(A, N=4096, inflate=0.0):
    """rho(theta) of the boundary of W(A) (optionally offset outward by `inflate`)
    about the centroid; returns interpolant via dense sampling."""
    thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
    _, zs = nr_support(A, thetas)
    zs = zs + inflate * np.exp(1j * thetas)
    zc = np.mean(zs)  # interior point (convexity)
    w = zs - zc
    ang = np.unwrap(np.angle(w))
    rad = np.abs(w)
    # ensure increasing angle parameterization
    if ang[-1] < ang[0]:
        ang, rad = ang[::-1], rad[::-1]
    return zc, ang, rad


def rho_of(theta, ang, rad):
    """Interpolate radius at polar angle theta (periodic)."""
    t0 = ang[0]
    per = 2 * np.pi
    x = np.mod(theta - t0, per) + t0
    return np.interp(x, ang, rad, period=per)


def theodorsen_map(A, N=1024, inflate=0.0, iters=200, tol=1e-13):
    """Return boundary z(t), z'(t) (t-grid uniform on [0,2pi)) of a conformal map
    D -> int(W(A) + inflate-offset), plus convergence diagnostics."""
    zc, ang, rad = boundary_polar(A, N=8 * N, inflate=inflate)
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    theta = t.copy()
    prev = None
    for it in range(iters):
        lr = np.log(rho_of(theta, ang, rad))
        theta_new = t + hilbert_periodic(lr)
        err = np.max(np.abs(theta_new - theta))
        theta = 0.5 * theta + 0.5 * theta_new  # damped for robustness (convex => mild)
        if prev is not None and err < tol:
            break
        prev = err
    z = zc + rho_of(theta, ang, rad) * np.exp(1j * theta)
    k = np.fft.fftfreq(N, d=1.0 / N)
    zp = np.fft.ifft(1j * k * np.fft.fft(z))
    return z, zp, err


class GeneralPullback:
    """Pullback-machinery clone driven by numerical boundary data z(t), z'(t)."""

    def __init__(self, z, zp):
        self.N = len(z)
        self.w = np.exp(1j * np.linspace(0, 2 * np.pi, self.N, endpoint=False))
        self.z = z
        self.dz = zp

    # duck-typed methods matching extremal_pullback.Pullback
    def resolvent_stack(self, A):
        n = A.shape[0]
        identity = np.eye(n, dtype=complex)
        T = np.empty((self.N, n, n), dtype=complex)
        for j in range(self.N):
            T[j] = (self.dz[j] * (2 * np.pi / self.N) / (2j * np.pi)) * np.linalg.inv(
                self.z[j] * identity - A
            )
        return T

    def calc(self, hvals, A, T=None):
        if T is None:
            T = self.resolvent_stack(A)
        return np.einsum("j,jkl->kl", hvals, T)

    def dlp_certificate(self, A):
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

    def unitality_certificate(self):
        """calc of h=1 must equal I for any inside point: use Cauchy of 1 at centroid."""
        zc = np.mean(self.z)
        val = np.sum(self.dz / (self.z - zc)) * (2 * np.pi / self.N) / (2j * np.pi)
        return abs(val - 1.0)

    def extremal_data(self, A, alphas, T=None):
        # zeros of f0 in z-plane and residues: need psi(alpha) and psi'(alpha) for
        # interior alpha: evaluate by Cauchy integrals of psi from boundary.
        if T is None:
            T = self.resolvent_stack(A)
        from extremal_pullback import blaschke, blaschke_prime_at_zero

        alphas = list(alphas)
        for i in range(len(alphas)):
            for j in range(i):
                if abs(alphas[i] - alphas[j]) < 1e-5:
                    alphas[i] += 3e-5 * np.exp(2j * np.pi * (i + 1) / 7)
        Bv = blaschke(self.w, alphas)
        f0A = self.calc(Bv, A, T)
        K = np.linalg.norm(f0A, 2)
        U, s, Vh = np.linalg.svd(f0A)
        x0 = Vh[0].conj()
        u0 = U[:, 0]
        diag = abs(np.vdot(x0, f0A @ x0))
        # psi(alpha), psi'(alpha) via Cauchy integrals over the unit circle in w:
        # psi(a) = (1/2pi i) oint z(t) / (w - a) dw,  dw = i w dt
        dw = 1j * self.w * (2 * np.pi / self.N)
        zi, dpsi_i = [], []
        for a in alphas:
            zi.append(np.sum(self.z / (self.w - a) * dw) / (2j * np.pi))
            dpsi_i.append(np.sum(self.z / (self.w - a) ** 2 * dw) / (2j * np.pi))
        zi = np.array(zi)
        dpsi_i = np.array(dpsi_i)
        ri = np.array(
            [dpsi_i[i] / blaschke_prime_at_zero(alphas, i) for i in range(len(alphas))]
        )
        g0v = 1.0 / Bv - sum(r / (self.z - z0) for r, z0 in zip(ri, zi))
        g0A = self.calc(g0v, A, T)
        C = x0.conj() @ (g0A @ f0A @ x0)
        Gv = g0A.conj().T @ x0
        v = f0A @ x0 + Gv
        W = np.linalg.norm(v)
        beta = x0.conj() @ v
        q = np.sqrt(max(W**2 - abs(beta) ** 2, 0.0))
        s_phase = 2 + C.real / 2 - q
        return dict(
            K=K,
            c=abs(C),
            C=C,
            G=np.linalg.norm(Gv),
            W=W,
            beta=beta,
            diag=diag,
            q=q,
            s_phase=s_phase,
            lhs=q,
            rhs=2 - abs(C) / 2,
            slack=2 - abs(C) / 2 - q,
            x0=x0,
            u0=u0,
        )
