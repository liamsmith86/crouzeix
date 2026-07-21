"""Sanity checks for crouzeix.py: known equality cases and normal matrices."""
import numpy as np
from crouzeix import (ratio_inner, ratio_outer, crabb_matrix, opnorm, poly_A,
                      inner_boundary_points)

rng = np.random.default_rng(0)

def report(name, A, c):
    ri = ratio_inner(A, c, ntheta=1024)
    ro = ratio_outer(A, c, ntheta=1024)
    print(f"{name:35s} R_inner={ri:.8f}  R_outer={ro:.8f}")

# 1. Jordan 2x2, p = z. Expect exactly 2.
A = np.array([[0, 2], [0, 0]], dtype=complex)
report("Jordan2 p=z", A, [0, 1])

# 2. Crabb matrices, p = z^k. Expect 2.
for k in [2, 3, 4, 5]:
    A = crabb_matrix(k)
    c = [0] * k + [1]
    # check numerical radius ~ 1
    zs = inner_boundary_points(A, 2048)
    print(f"  crabb k={k}: ||A^k||={opnorm(poly_A(c, A)):.6f}  w(A)~{np.max(np.abs(zs)):.8f}")
    report(f"Crabb k={k} p=z^{k}", A, c)

# 3. Random normal matrices: R <= 1.
for t in range(3):
    D = np.diag(rng.standard_normal(5) + 1j * rng.standard_normal(5))
    Q, _ = np.linalg.qr(rng.standard_normal((5, 5)) + 1j * rng.standard_normal((5, 5)))
    A = Q @ D @ Q.conj().T
    c = rng.standard_normal(4) + 1j * rng.standard_normal(4)
    report(f"normal-{t}", A, list(c))

# 4. Random nonnormal: expect R < 2 (conjecture), often well below.
best = 0
for t in range(200):
    n = rng.integers(3, 6)
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    c = list(rng.standard_normal(5) + 1j * rng.standard_normal(5))
    r = ratio_inner(A, c, ntheta=128)
    best = max(best, r)
print(f"random nonnormal max R_inner over 200 trials: {best:.6f}")
