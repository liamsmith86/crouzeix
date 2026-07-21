"""Epoch 5 opener: doubly-symmetric 4x4 zero-diag tridiagonal — BEYOND known classes.

A = [[0,a1,0,0],[b1,0,a2,0],[0,b2,0,a3],[0,0,b3,0]]; sigma(A) = {±e1, ±e2};
A^4 = pA^2 - qI. W(A) generally NOT elliptic (first non-covered family).

Probe at near-critical domain: extremal Blaschke (deg <= 3) structure (parity of
zeros?), rho >= 0?, and collapse pattern (is f0 odd: zeros {0, ±alpha}?).
Also test: is the squared boundary {z^2} an ellipse (residual)? Expect NO.
"""
import numpy as np

from crouzeix import nr_support
from theodorsen import theodorsen_map, GeneralPullback
from minkowski_test import best_extremal

def probe(ws, label, t=0.005, seed=0):
    a1, b1, a2, b2, a3, b3 = ws
    A = np.zeros((4, 4), dtype=complex)
    A[0, 1], A[1, 0] = a1, b1
    A[1, 2], A[2, 1] = a2, b2
    A[2, 3], A[3, 2] = a3, b3
    eigs = np.linalg.eigvals(A)
    th = np.linspace(0, 2 * np.pi, 2000, endpoint=False)
    _, zs = nr_support(A, th)
    # ellipse-fit residual of squared boundary
    zeta = zs ** 2
    x, y = zeta.real, zeta.imag
    M = np.column_stack([x * x, x * y, y * y, x, y])
    coef, *_ = np.linalg.lstsq(M, np.ones_like(x), rcond=None)
    resid = np.abs(M @ coef - 1).max()
    wmin = (lambda hs: (hs + hs[np.arange(len(hs)) - len(hs) // 2]).min())(
        nr_support(A, np.linspace(0, 2 * np.pi, 64, endpoint=False))[0])
    if wmin < 0.25:
        print(f"{label}: too thin, skip")
        return
    try:
        z, zp, terr = theodorsen_map(A, N=1024, inflate=t)
    except Exception as ex:
        print(f"{label}: map fail {ex}")
        return
    pb = GeneralPullback(z, zp)
    lmin, merr = pb.dlp_certificate(A)
    if lmin < -1e-9 or merr > 1e-6:
        print(f"{label}: cert fail")
        return
    al, d = best_extremal(pb, A, 3, seed)
    print(f"{label}: eigs={np.round(np.sort_complex(eigs), 3)}")
    print(f"  sq-ellipse-resid={resid:.2e}  K={d['K']:.6f} rho={d['C'].real:+.6f} "
          f"ImC={d['C'].imag:+.2e} diag={d['diag']:.1e}")
    print(f"  zeros={np.round(np.array(al), 4)}")

if __name__ == "__main__":
    cases = [
        ((2.0, 0.1, 1.2, 0.2, 1.6, 0.1), "sym4-A"),
        ((1.8, 0.3, 1.0, -0.2, 1.4, 0.2), "sym4-B"),
        ((1.5, 0.2, 2.0, 0.1, 0.9, -0.3), "sym4-C"),
    ]
    for i, (ws, label) in enumerate(cases):
        probe(ws, label, seed=i)
