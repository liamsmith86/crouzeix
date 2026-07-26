"""Probe: geometry of extremal Blaschke zeros vs eigenvalue images at critical domains.

For matrices A (n=3) at near-critical domain: compute w_j = phi(eig_j) (images of
eigenvalues under the Riemann map, via Cauchy integral of the boundary correspondence)
and extremal zeros alpha_i. Test candidate laws:
  (a) alpha_i ~ critical points of prod_j b_{w_j}(w) (hyperbolic Fekete/median)?
  (b) alpha_i ~ zeros of B such that Pick matrix of (w_j -> B... ) singular (consistency).
Report raw data for pattern hunting.
"""

import numpy as np

from theodorsen import theodorsen_map, GeneralPullback
from minkowski_test import best_extremal


def phi_of_points(pb, pts):
    """phi(p) for interior p via boundary correspondence: phi = psi^{-1};
    psi maps w-circle to z-boundary. Evaluate w(p) by Cauchy integral of the
    inverse boundary function: w(p) = (1/2pi i) oint w(t) psi'(t)/(psi(t)-p) dt...
    i.e. since w -> psi(w) conformal, phi(p) = (1/2pi i) oint_{|w|=1} w psi'(w)/(psi(w)-p) dw."""
    dw = 1j * pb.w * (2 * np.pi / pb.N)
    dpsi_dw = pb.dz / (1j * pb.w)  # dz/dtheta = psi'(w) * i w
    out = []
    for p in pts:
        val = np.sum(pb.w * dpsi_dw / (pb.z - p) * dw) / (2j * np.pi)
        out.append(val)
    return np.array(out)


def hyperbolic_critical_points(wjs):
    """Critical points of the Blaschke product with zeros at w_j (candidate law:
    extremal zeros = critical points of the 'eigenvalue Blaschke')."""

    # B(w) = prod (w - w_j)/(1 - conj(w_j) w); B'/B = sum [1/(w-w_j) + conj(w_j)/(1-conj(w_j)w)]
    # critical points: sum_j (1-|w_j|^2) / ((w-w_j)(1-conj(w_j)w)) = 0
    # clear denominators -> polynomial; roots inside D are the hyperbolic critical points.
    n = len(wjs)
    num = np.zeros(1, dtype=complex)
    for j in range(n):
        term = np.array([1.0], dtype=complex) * (1 - abs(wjs[j]) ** 2)
        for k in range(n):
            if k != j:
                fac1 = np.array([-wjs[k], 1.0])  # (w - w_k)
                term = np.polymul(
                    np.polymul(term, fac1[::-1])[::-1]
                    if False
                    else np.polymul(term[::-1], fac1[::-1])[::-1],
                    np.array([1.0]),
                )
        num = num  # placeholder — do numerically instead

    # numeric root find on the rational equation instead:
    def g(w):
        return sum(
            (1 - abs(wj) ** 2) / ((w - wj) * (1 - np.conj(wj) * w)) for wj in wjs
        )

    # sample grid, Newton polish
    roots = []
    for r0 in np.linspace(0.05, 0.9, 8):
        for t0 in np.linspace(0, 2 * np.pi, 16, endpoint=False):
            w = r0 * np.exp(1j * t0)
            for _ in range(60):
                h = 1e-7
                dg = (g(w + h) - g(w - h)) / (2 * h)
                if abs(dg) < 1e-14:
                    break
                w2 = w - g(w) / dg
                if abs(w2) > 1.5:
                    break
                if abs(w2 - w) < 1e-13:
                    w = w2
                    break
                w = w2
            if abs(w) < 0.999 and abs(g(w)) < 1e-9:
                if not any(abs(w - r) < 1e-6 for r in roots):
                    roots.append(w)
    return np.array(roots)


if __name__ == "__main__":
    rng = np.random.default_rng(17)
    from crouzeix import crabb_matrix

    mats = []
    E = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    E /= np.linalg.norm(E, 2)
    mats.append(("crabb2+0.3E", crabb_matrix(2) + 0.3 * E))
    mats.append(
        (
            "tri(2,.2,1.4,.1)",
            np.array([[0, 2.0, 0], [0.2, 0, 1.4], [0, 0.1, 0]], dtype=complex),
        )
    )
    A = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    mats.append(("random", A))
    for label, A in mats:
        try:
            z, zp, terr = theodorsen_map(A, N=1024, inflate=0.005)
        except Exception as e:
            print(label, "map fail", e)
            continue
        pb = GeneralPullback(z, zp)
        al, d = best_extremal(pb, A, 2, 3)
        eigs = np.linalg.eigvals(A)
        wj = phi_of_points(pb, eigs)
        crit = hyperbolic_critical_points(wj)
        print(f"{label}: K={d['K']:.5f} rho={d['C'].real:+.6f} diag={d['diag']:.1e}")
        print(f"   eig->D: {np.round(wj, 4)}")
        print(f"   extremal zeros: {np.round(np.array(al), 4)}")
        print(f"   B'-critical pts of eigenvalue-Blaschke: {np.round(crit, 4)}")
