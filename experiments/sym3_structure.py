"""Structure probe: doubly-symmetric 3x3 family A = [[0,a,0],[b,0,c],[0,d,0]] (real).
W(A) symmetric about both axes; A ~ -A and A real => extremal Blaschke zeros closed
under alpha -> conj(alpha) and alpha -> -alpha.

At near-critical domain (offset t): compute extremal (deg <= 2), report:
 - zero pattern of B, K, rho, per-zero capacity terms  T_i = r_i <h_i(A)x0,x0>
   (rho = 1 - Re sum T_i), and the 2x2-style ratio test.
"""
import numpy as np

from theodorsen import theodorsen_map, GeneralPullback
from minkowski_test import best_extremal

def probe(A, label, t=0.005, N=1024, seed=0):
    pb = None
    try:
        z, zp, terr = theodorsen_map(A, N=N, inflate=t)
        pb = GeneralPullback(z, zp)
    except Exception as e:
        print(f"{label}: map failed {e}")
        return
    lmin, merr = pb.dlp_certificate(A)
    if lmin < -1e-9 or merr > 1e-6:
        print(f"{label}: cert fail lmin={lmin:.1e}")
        return
    al, d = best_extremal(pb, A, 2, seed)
    # per-zero capacity terms via the g0 decomposition pieces
    from extremal_pullback import blaschke, blaschke_prime_at_zero
    T = pb.resolvent_stack(A)
    Bv = blaschke(pb.w, al)
    f0A = pb.calc(Bv, A, T)
    U, s, Vh = np.linalg.svd(f0A)
    x0 = Vh[0].conj()
    dw = 1j * pb.w * (2 * np.pi / pb.N)
    terms = []
    for i, a in enumerate(al):
        zi = np.sum(pb.z / (pb.w - a) * dw) / (2j * np.pi)
        dpsi = np.sum(pb.z / (pb.w - a) ** 2 * dw) / (2j * np.pi)
        ri = dpsi / blaschke_prime_at_zero(al, i)
        hi = Bv / (pb.z - zi)
        hiA = pb.calc(hi, A, T)
        Ti = ri * (x0.conj() @ (hiA @ x0))
        terms.append(Ti)
    rho = d["C"].real
    print(f"{label}: K={d['K']:.6f} rho={rho:+.6f} diag={d['diag']:.1e} "
          f"zeros={np.round(np.array(al), 4)}")
    print(f"    terms T_i = {[f'{T_.real:+.5f}{T_.imag:+.5f}j' for T_ in terms]}  "
          f"1-Re(sum)={1 - sum(T_.real for T_ in terms):+.6f}")

if __name__ == "__main__":
    from crouzeix import nr_support

    def thickness(A):
        th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        hs, _ = nr_support(A, th)
        widths = hs + hs[np.arange(64) - 32]
        return widths.min(), widths.max()

    cases = [
        (2.0, 0.0, 1.4, 0.0),    # pure shift weights (Crabb-like)
        (2.0, 0.2, 1.4, 0.1),    # weakly two-sided
        (1.8, -0.4, 0.9, 0.2),
        (1.5, 0.3, 2.0, -0.2),
    ]
    for i, (a, b, c, dd) in enumerate(cases):
        A = np.array([[0, a, 0], [b, 0, c], [0, dd, 0]], dtype=complex)
        wmin, wmax = thickness(A)
        if wmin < 0.25:
            print(f"tri({a},{b},{c},{dd}): too thin (wmin={wmin:.3f}), skip")
            continue
        probe(A, f"tri({a},{b},{c},{dd})", seed=i)
