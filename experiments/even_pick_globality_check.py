"""Regression checks for the global even-midpoint theorem.

The analytic proof is ``proof/even_pick_globality.md``.  This script checks the determinant
identity and searches a deterministic complex automorphism grid for an orbit value larger than
the predicted maximum.
"""

from __future__ import annotations

import numpy as np


def automorphism_matrix(d: float, h: float, a: complex) -> np.ndarray:
    involution = np.array([[1.0, 2.0 * h], [0.0, -1.0]], dtype=complex)
    operator = d * involution
    identity = np.eye(2, dtype=complex)
    return (operator + a * identity) @ np.linalg.inv(identity + np.conj(a) * operator)


def determinant_formula(d: float, h: float, a: complex) -> tuple[float, float]:
    midpoint = automorphism_matrix(d, h, 0j)
    lam = np.linalg.norm(midpoint, 2) ** 2
    transformed = automorphism_matrix(d, h, a)
    direct = np.linalg.det(lam * np.eye(2) - transformed.conj().T @ transformed).real

    x = a.real
    r = abs(a) ** 2
    j_value = (
        r * ((lam + 1) * (d**4 - 1) * r + 2 * (1 + d * d) * (lam - d * d))
        - 4 * d * d * (lam - 1) * x * x
    )
    denominator = abs(1 + np.conj(a) * d) ** 2 * abs(1 - np.conj(a) * d) ** 2
    predicted = (lam - 1) * j_value / denominator
    return direct, predicted


def main() -> None:
    max_determinant_error = 0.0
    max_orbit_excess = -np.inf
    tested_superunit_cases = 0

    for d in (0.2, 0.5, 0.8):
        for h in (0.1, 0.5, 2.0, 10.0):
            midpoint_norm = np.linalg.norm(automorphism_matrix(d, h, 0j), 2)
            if midpoint_norm <= 1:
                continue
            tested_superunit_cases += 1
            for radius in np.linspace(0.025, 0.975, 39):
                for angle in np.linspace(0, 2 * np.pi, 80, endpoint=False):
                    a = radius * np.exp(1j * angle)
                    direct, predicted = determinant_formula(d, h, a)
                    relative_error = abs(direct - predicted) / (1 + abs(direct))
                    max_determinant_error = max(max_determinant_error, relative_error)
                    orbit_norm = np.linalg.norm(automorphism_matrix(d, h, a), 2)
                    max_orbit_excess = max(max_orbit_excess, orbit_norm - midpoint_norm)

    assert tested_superunit_cases > 0
    assert max_determinant_error < 2e-12
    assert max_orbit_excess < 0
    print(f"superunit midpoint cases: {tested_superunit_cases}")
    print(f"max determinant relerr:   {max_determinant_error:.3e}")
    print(f"max orbit excess:         {max_orbit_excess:.3e}")


if __name__ == "__main__":
    main()
