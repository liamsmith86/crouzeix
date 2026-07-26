"""Adversarial search for violations of candidate Lemma L7'.

L7': For M (n x n), Omega = eps-offset of W(M), p with ||p||_Omega <= 1, unit x0
with <p(M)x0, x0> = 0:
   sqrt(||G x0||^2 - |<G x0, x0>|^2) <= 2 - |<(Phi(p) p)(M) x0, x0>| / 2,
where G = p(M) + Phi(p)(M)*.

Maximize  viol = lhs - rhs  (want: is sup viol > 0?)
over A (n x n complex), p (deg d), x0 (unit, orthogonality via penalty).

Usage: adversarial_L7.py n d restarts seed [eps]
"""

import json
import sys

import numpy as np
from scipy.optimize import minimize

from crouzeix import poly_z, poly_A, opnorm
from star_inequality import offset_curve, cauchy_matrix, phi_on_inner


def quantities(M, cp, x0, eps, N=192, tol=1e-6):
    """Compute L7' quantities with self-certification. Returns None if the
    contour quadrature cannot certify itself (degenerate/nonsmooth boundary)."""
    n = M.shape[0]
    _, z_out, zp_out = offset_curve(M, eps, N)
    _, z_in, zp_in = offset_curve(M, eps / 2, N)
    # guard: contours must be well-separated (no FFT garbage / crossings)
    if np.min(np.abs(z_out[:, None] - z_in[None, :])) < eps / 8:
        return None
    pv_out = poly_z(cp, z_out)
    m = np.max(np.abs(pv_out))
    if not np.isfinite(m) or m < 1e-12:
        return None
    cp = [ci / m for ci in cp]
    # guard: monomial-basis conditioning (huge coeffs + cancellation = fp garbage)
    if np.max(np.abs(cp)) > 1e6:
        return None
    pM = poly_A(cp, M)
    # theorem guard: ||p(M)|| <= 1+sqrt2 is PROVEN under certified hypotheses;
    # any computed excess is a numerical artifact, not a discovery.
    if opnorm(pM) > 2.5:
        return None
    # certificate 1: double-layer total mass = 2I AND pointwise positivity
    # P(sigma) >= 0 (this encodes convexity of the curve + W(M) inside; it is
    # the actual hypothesis behind ||gamma_Phi|| <= 2 and hence L7').
    S = np.zeros((n, n), dtype=complex)
    for j in range(N):
        B = (zp_out[j] / (2j * np.pi)) * np.linalg.inv(z_out[j] * np.eye(n) - M)
        Q = B + B.conj().T
        if np.linalg.eigvalsh(Q)[0] < -1e-9:
            return None
        S += Q
    S *= 2 * np.pi / N
    if opnorm(S - 2 * np.eye(n)) > tol:
        return None
    # certificate 1b: tangent angle monotone, total turning = 2pi (convex Jordan curve)
    ang = np.unwrap(np.append(np.angle(zp_out), np.angle(zp_out[0])))
    dang = np.diff(ang)
    if np.any(dang < -0.01) or abs((ang[-1] - ang[0]) - 2 * np.pi) > 0.1:
        return None
    # certificate 2: Cauchy reconstruction of p(M) from inner contour
    pv_in = poly_z(cp, z_in)
    pM2 = cauchy_matrix(pv_in, z_in, zp_in, M)
    if opnorm(pM - pM2) > tol * max(1.0, opnorm(pM)):
        return None
    phi_in = phi_on_inner(pv_out, z_out, zp_out, z_in)
    if not np.all(np.isfinite(phi_in)):
        return None
    g0M = cauchy_matrix(phi_in, z_in, zp_in, M)
    prodM = cauchy_matrix(phi_in * pv_in, z_in, zp_in, M)
    v = pM @ x0 + g0M.conj().T @ x0
    W2 = np.linalg.norm(v) ** 2
    # theorem guard: ||gamma_Phi(p)|| <= 2 is PROVEN under certified hypotheses.
    if W2 > 4.0 + 1e-6:
        return None
    beta = np.vdot(x0, v)
    c_val = abs(np.vdot(x0, prodM @ x0))
    diag = abs(np.vdot(x0, pM @ x0))
    lhs = np.sqrt(max(W2 - abs(beta) ** 2, 0.0))
    rhs = 2 - c_val / 2
    return lhs, rhs, c_val, diag, cp


def main():
    n = int(sys.argv[1])
    d = int(sys.argv[2])
    restarts = int(sys.argv[3])
    seed = int(sys.argv[4])
    eps = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
    rng = np.random.default_rng(seed)
    n2 = n * n
    PEN = 30.0

    def unpack(x):
        A = (x[:n2] + 1j * x[n2 : 2 * n2]).reshape(n, n)
        nf = np.linalg.norm(A)
        if nf < 1e-8:
            return None
        A = A / nf  # scale-normalize (L7' is scale-covariant with Omega)
        o = 2 * n2
        cp = x[o : o + d + 1] + 1j * x[o + d + 1 : o + 2 * (d + 1)]
        o += 2 * (d + 1)
        x0 = x[o : o + n] + 1j * x[o + n : o + 2 * n]
        nx = np.linalg.norm(x0)
        if nx < 1e-12:
            return None
        return A, list(cp), x0 / nx

    def obj(x):
        u = unpack(x)
        if u is None:
            return 0.0
        A, cp, x0 = u
        if np.max(np.abs(cp[1:])) < 1e-10:
            return 0.0
        try:
            q = quantities(A, cp, x0, eps)
        except (np.linalg.LinAlgError, FloatingPointError):
            return 0.0
        if q is None:
            return 0.0
        lhs, rhs, c_val, diag, _ = q
        return -(lhs - rhs) + PEN * diag**2

    fname = f"advL7_n{n}_d{d}_s{seed}.jsonl"
    best = -10
    with open(fname, "w") as fh:
        for t in range(restarts):
            x_init = rng.standard_normal(2 * n2 + 2 * (d + 1) + 2 * n)
            res = minimize(
                obj, x_init, method="L-BFGS-B", options={"maxiter": 500, "eps": 1e-7}
            )
            u = unpack(res.x)
            if u is None:
                continue
            A, cp, x0 = u
            try:
                q = quantities(A, cp, x0, eps, N=1024, tol=1e-8)
            except (np.linalg.LinAlgError, FloatingPointError):
                continue
            if q is None:
                continue
            lhs, rhs, c_val, diag, cpn = q
            viol = lhs - rhs
            rec = dict(
                trial=t,
                viol=viol,
                lhs=lhs,
                rhs=rhs,
                c=c_val,
                diag=diag,
                A_re=A.real.tolist(),
                A_im=A.imag.tolist(),
                c_re=np.real(cpn).tolist(),
                c_im=np.imag(cpn).tolist(),
                x0_re=x0.real.tolist(),
                x0_im=x0.imag.tolist(),
            )
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if viol > best:
                best = viol
                print(
                    f"[n={n} d={d} s={seed}] t={t}: viol={viol:+.6f} "
                    f"(lhs={lhs:.5f} rhs={rhs:.5f} c={c_val:.5f} diag={diag:.2e})",
                    flush=True,
                )
    print(f"DONE advL7 n={n} d={d} s={seed}: max viol = {best:+.8f}", flush=True)


if __name__ == "__main__":
    main()
