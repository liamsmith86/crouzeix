"""Audit parity and the shifted-Mobius phase on the exact elliptic 4x4 slice.

The historical notes inferred that ``A ~ -A`` forces an extremal of definite parity.  Symmetry
does not imply that in general, and this script gives a robust slice configuration where a
shifted degree-one Blaschke product beats the computed odd and the rigorously characterized even
sectors.  The full phase-classification issue is documented in ``proof/slice_closed_form.md``.

This is high-precision numerical evidence, not a certified global proof over the odd sector.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

from slice_exact import Slice


WEIGHTS = (0.8, 2.4, 1.1, 0.3)


@dataclass(frozen=True)
class NodalModel:
    slice_data: Slice
    matrix: np.ndarray
    eigenvalues: np.ndarray
    disk_nodes: np.ndarray
    eigenvectors: np.ndarray
    eigenvectors_inv: np.ndarray

    @classmethod
    def build(cls, weights: tuple[float, float, float, float]) -> "NodalModel":
        data = Slice(*weights)
        a1, a2, a3 = (float(value) for value in data.a)
        c = float(data.c)
        matrix = np.array(
            [
                [0, a1, 0, 0],
                [c * a1, 0, a2, 0],
                [0, c * a2, 0, a3],
                [0, 0, c * a3, 0],
            ],
            dtype=float,
        )
        eigenvalues, eigenvectors = np.linalg.eig(matrix)
        eigenvectors_inv = np.linalg.inv(eigenvectors)
        positive_eigenvalues = np.array([float(data.e1), float(data.e2)])
        positive_nodes = np.array(
            [float(complex(data.tau[0]).real), float(complex(data.tau[1]).real)]
        )
        disk_nodes = []
        for eigenvalue in eigenvalues:
            index = np.argmin(abs(abs(eigenvalue) - positive_eigenvalues))
            disk_nodes.append(np.sign(eigenvalue.real) * positive_nodes[index])
        return cls(
            data,
            matrix,
            eigenvalues,
            np.asarray(disk_nodes),
            eigenvectors,
            eigenvectors_inv,
        )

    def calculus(self, nodal_values: np.ndarray) -> np.ndarray:
        return self.eigenvectors @ np.diag(nodal_values) @ self.eigenvectors_inv

    def norm(self, nodal_values: np.ndarray) -> float:
        return float(np.linalg.norm(self.calculus(nodal_values), 2))


def disk_parameter(raw: np.ndarray) -> complex:
    value = raw[0] + 1j * raw[1]
    return value / np.sqrt(1 + abs(value) ** 2)


def shifted_mobius(model: NodalModel) -> tuple[float, float, float, float]:
    def norm_at(beta: float) -> float:
        values = (model.disk_nodes - beta) / (1 - beta * model.disk_nodes)
        return model.norm(values)

    result = minimize_scalar(
        lambda beta: -norm_at(beta),
        bounds=(-0.999, 0.999),
        method="bounded",
        options={"xatol": 1e-15},
    )
    beta = float(result.x)
    values = (model.disk_nodes - beta) / (1 - beta * model.disk_nodes)
    operator = model.calculus(values)
    _, singular_values, vh = np.linalg.svd(operator)
    vector = vh[0].conj()
    diagonal = abs(vector.conj() @ operator @ vector)

    zero, inverse_derivative = model.slice_data.psi_data(beta)
    residue = float(inverse_derivative) * (1 - beta * beta)
    g_values = 1 / values - residue / (model.eigenvalues - float(zero))
    product_operator = model.calculus(values * g_values)
    rho = float((vector.conj() @ product_operator @ vector).real)
    return beta, float(singular_values[0]), float(diagonal), rho


def even_sector_global_value(model: NodalModel) -> float:
    positive_nodes = np.sort(np.unique(abs(model.disk_nodes)))
    u2, u1 = positive_nodes**2
    midpoint = np.tanh((np.atanh(u1) + np.atanh(u2)) / 2)
    values = (model.disk_nodes**2 - midpoint) / (1 - midpoint * model.disk_nodes**2)
    return max(1.0, model.norm(values))


def odd_sector_numerical_value(model: NodalModel) -> float:
    squared_nodes = model.disk_nodes**2

    def objective(raw: np.ndarray) -> float:
        parameter = disk_parameter(raw)
        factor = (squared_nodes - parameter) / (1 - np.conj(parameter) * squared_nodes)
        return -model.norm(model.disk_nodes * factor)

    result = differential_evolution(
        objective,
        [(-100, 100), (-100, 100)],
        seed=20260721,
        popsize=24,
        maxiter=500,
        tol=1e-11,
        polish=True,
    )
    # The boundary constant factor gives f(w)=w and can dominate all interior odd products.
    return max(-float(result.fun), model.norm(model.disk_nodes))


def main() -> None:
    model = NodalModel.build(WEIGHTS)
    beta, shifted_norm, diagonal, rho = shifted_mobius(model)
    even_norm = even_sector_global_value(model)
    odd_norm = odd_sector_numerical_value(model)

    assert shifted_norm > odd_norm + 0.05
    assert shifted_norm > even_norm + 0.05
    assert diagonal < 1e-7
    assert rho > 0

    print(f"weights:                 {WEIGHTS}")
    print(f"shifted Mobius beta:     {beta:.12f}")
    print(f"shifted Mobius norm:     {shifted_norm:.12f}")
    print(f"odd-sector numerical:    {odd_norm:.12f}")
    print(f"even-sector global:      {even_norm:.12f}")
    print(f"extremal diagonal:       {diagonal:.3e}")
    print(f"rho at shifted phase:    {rho:+.12f}")


if __name__ == "__main__":
    main()
