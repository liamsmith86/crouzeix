"""Test Lemma M: rho(s) along the Minkowski family Omega_s = (1-s) W(A) + s D_R.

Support functions add under Minkowski sum: h_s(theta) = (1-s) h_W(theta) + s h_D(theta),
h_D(theta) = Re(c e^{-i theta}) + R for disk D(c, R). Boundary point of a smooth
convex body from its support function: z(theta) = e^{i theta}(h(theta) + i h'(theta)).

We drive Theodorsen with this polar boundary (about the centroid), then the
GeneralPullback extremal machinery. rho(1) = 0 expected (disk); Lemma M predicts
rho(0+) >= 0 with rho monotone-ish decreasing in s.

Usage: minkowski_test.py [seed] [n]
"""
import sys

import numpy as np

from crouzeix import nr_support, crabb_matrix
from extremal_pullback import find_extremal_blaschke
from theodorsen import hilbert_periodic, GeneralPullback

def support_W(A, N):
    th = np.linspace(0, 2 * np.pi, N, endpoint=False)
    hs, _ = nr_support(A, th)
    return th, hs

def boundary_from_support(th, h):
    """z(theta) = e^{i theta} (h + i h') with h' by FFT."""
    N = len(th)
    k = np.fft.fftfreq(N, d=1.0 / N)
    hp = np.real(np.fft.ifft(1j * k * np.fft.fft(h)))
    return np.exp(1j * th) * (h + 1j * hp)

def theodorsen_from_boundary(zs, N):
    """Theodorsen iteration given dense boundary samples zs (any convex curve)."""
    zc = np.mean(zs)
    w = zs - zc
    ang = np.unwrap(np.angle(w))
    rad = np.abs(w)
    if ang[-1] < ang[0]:
        ang, rad = ang[::-1], rad[::-1]
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    theta = t.copy()
    per = 2 * np.pi
    t0 = ang[0]

    def rho_of(x):
        return np.interp(np.mod(x - t0, per) + t0, ang, rad, period=per)

    err = 1.0
    for _ in range(400):
        lr = np.log(rho_of(theta))
        theta_new = t + hilbert_periodic(lr)
        err = np.max(np.abs(theta_new - theta))
        theta = 0.5 * theta + 0.5 * theta_new
        if err < 1e-13:
            break
    z = zc + rho_of(theta) * np.exp(1j * theta)
    k = np.fft.fftfreq(N, d=1.0 / N)
    zp = np.fft.ifft(1j * k * np.fft.fft(z))
    return z, zp, err

def best_extremal(pb, A, maxdeg, seed):
    """Strong extremal search: per-degree multistart + Nelder-Mead polish."""
    from scipy.optimize import minimize
    from extremal_pullback import blaschke
    T = pb.resolvent_stack(A)
    best = None
    for deg in range(1, maxdeg + 1):
        def unpack(x, dg=deg):
            pts = x[:dg] + 1j * x[dg:]
            return [p / np.sqrt(1 + abs(p) ** 2) for p in pts]

        def obj(x):
            return -np.linalg.norm(pb.calc(blaschke(pb.w, unpack(x)), A, T), 2)

        rng = np.random.default_rng(seed + 100 * deg)
        cand = None
        for t in range(24):
            x0 = rng.standard_normal(2 * deg) * (0.6 if t % 2 else 1.6)
            res = minimize(obj, x0, method="Nelder-Mead",
                           options={"maxiter": 4000, "fatol": 1e-14, "xatol": 1e-12})
            if cand is None or res.fun < cand.fun:
                cand = res
        # polish: gradient method then tight Nelder-Mead
        res = minimize(obj, cand.x, method="L-BFGS-B",
                       options={"maxiter": 2000, "eps": 1e-9, "ftol": 1e-16, "gtol": 1e-12})
        if res.fun > cand.fun:
            res = cand
        res2 = minimize(obj, res.x, method="Nelder-Mead",
                        options={"maxiter": 20000, "fatol": 1e-16, "xatol": 1e-14})
        if res2.fun < res.fun:
            res = res2
        al = unpack(res.x)
        d = pb.extremal_data(A, al, T)
        if best is None or d["K"] > best[1]["K"]:
            best = (al, d)
    return best

def run(A, label, svals, N=1024, seed=0, margin=1e-3):
    n = A.shape[0]
    thd, hW = support_W(A, 8 * N)
    # circumscribed-ish disk: center = centroid of W-boundary, R = max dist + margin
    _, zsW = nr_support(A, np.linspace(0, 2 * np.pi, 512, endpoint=False))
    c = np.mean(zsW)
    R = float(np.max(np.abs(zsW - c))) * 1.02
    hD = np.real(c * np.exp(-1j * thd)) + R
    print(f"--- {label} (disk c={c:.3f}, R={R:.3f})")
    for s in svals:
        h = (1 - s) * hW + s * hD + margin * (1 + 3 * (s == 0))
        zs = boundary_from_support(thd, h)
        z, zp, terr = theodorsen_from_boundary(zs, N)
        pb = GeneralPullback(z, zp)
        lmin, merr = pb.dlp_certificate(A)
        if lmin < -1e-9 or merr > 1e-6:
            print(f"  s={s}: cert fail lmin={lmin:.1e} mass={merr:.1e}")
            continue
        al, d = best_extremal(pb, A, n - 1, seed)
        print(f"  s={s:5.3f}: K={d['K']:.6f} rho={d['C'].real:+.6f} ImC={d['C'].imag:+.6f} "
              f"q={d['q']:.6f} s_phase={d['s_phase']:+.6f} diag={d['diag']:.1e}", flush=True)

if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    rng = np.random.default_rng(seed)
    svals = [0.0, 0.1, 0.25, 0.45, 0.7, 0.9, 1.0]
    E = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    E /= np.linalg.norm(E, 2)
    run(crabb_matrix(n - 1) + 0.4 * E, f"crabb{n-1}+0.4E s{seed}", svals, seed=seed)
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    run(A, f"random {n}x{n} s{seed}", svals, seed=seed)
