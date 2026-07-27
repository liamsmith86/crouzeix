#!/usr/bin/env python3
"""Independently finite-difference the exact Gau--Wu normal Hessian."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import sympy as sp

from gau_wu_transverse_hessian import expected_optimized_hessian
from theodorsen import GeneralPullback, theodorsen_map


@dataclass(frozen=True)
class GauWuNumericHessianRecord:
    """Independent true-map comparison with the exact normal Hessian."""

    model_parameter: float
    finite_difference_step: float
    map_resolution: int
    tangent_rank: int
    maximum_map_residual: float
    maximum_hessian_error: float
    numerical_eigenvalues: tuple[float, ...]
    exact_eigenvalues: tuple[float, ...]
    maximum_eigenvalue_error: float
    all_checks_passed: bool


def gau_wu_matrix(parameter: float) -> np.ndarray:
    """Return the floating three-dimensional Gau--Wu matrix."""

    edge = np.sqrt(2 * (1 - parameter**2))
    return np.array(
        [[0, edge, -2 * parameter], [0, parameter, edge], [0, 0, 0]],
        dtype=complex,
    )


def real_vectorize(matrix: np.ndarray) -> np.ndarray:
    """Vectorize into the interlaced real coordinates of the exact checker."""

    return np.array(
        [
            coordinate
            for entry in matrix.flat
            for coordinate in (entry.real, entry.imag)
        ]
    )


def real_matrix(vector: np.ndarray) -> np.ndarray:
    """Invert ``real_vectorize``."""

    return np.array(
        [
            vector[2 * index] + 1j * vector[2 * index + 1]
            for index in range(9)
        ]
    ).reshape(3, 3)


def equality_tangent(parameter: float) -> np.ndarray:
    """Return thirteen independent equality-tangent columns."""

    matrix = gau_wu_matrix(parameter)
    tangents: list[np.ndarray] = []
    for index in (0, 1):
        generator = np.zeros((3, 3), dtype=complex)
        generator[index, index] = 1j
        tangents.append(generator @ matrix - matrix @ generator)
    for row in range(3):
        for column in range(row + 1, 3):
            real_generator = np.zeros((3, 3), dtype=complex)
            real_generator[row, column] = 1
            real_generator[column, row] = -1
            tangents.append(
                real_generator @ matrix - matrix @ real_generator
            )
            imaginary_generator = np.zeros((3, 3), dtype=complex)
            imaginary_generator[row, column] = 1j
            imaginary_generator[column, row] = 1j
            tangents.append(
                imaginary_generator @ matrix - matrix @ imaginary_generator
            )
    identity = np.eye(3)
    tangents.extend((identity, 1j * identity, matrix, 1j * matrix))
    edge = np.sqrt(2 * (1 - parameter**2))
    radial = np.array(
        [
            [0, -2 * parameter / edge, -2],
            [0, 1, -2 * parameter / edge],
            [0, 0, 0],
        ],
        dtype=complex,
    )
    tangents.append(radial)
    return np.stack([real_vectorize(tangent) for tangent in tangents], axis=1)


def slice_basis() -> np.ndarray:
    """Return L332's five bottom-row complement columns."""

    matrices: list[np.ndarray] = []
    for row, column, imaginary in (
        (2, 0, False),
        (2, 0, True),
        (2, 1, False),
        (2, 1, True),
        (2, 2, False),
    ):
        matrix = np.zeros((3, 3), dtype=complex)
        matrix[row, column] = 1j if imaginary else 1
        matrices.append(matrix)
    return np.stack([real_vectorize(matrix) for matrix in matrices], axis=1)


def unpack_zeros(parameters: np.ndarray) -> np.ndarray:
    """Map two unconstrained complex coordinates into the unit disk."""

    points = parameters[:2] + 1j * parameters[2:]
    return np.array(
        [point / np.sqrt(1 + abs(point) ** 2) for point in points]
    )


def pack_zeros(zeros: np.ndarray) -> np.ndarray:
    """Invert ``unpack_zeros``."""

    points = zeros / np.sqrt(1 - abs(zeros) ** 2)
    return np.concatenate((points.real, points.imag))


def blaschke_at_matrix(
    matrix: np.ndarray,
    zeros: np.ndarray,
) -> np.ndarray:
    """Evaluate a degree-two Blaschke product directly."""

    identity = np.eye(matrix.shape[0])
    value = np.eye(matrix.shape[0], dtype=complex)
    for zero in zeros:
        value = (
            value
            @ (matrix - zero * identity)
            @ np.linalg.inv(identity - np.conj(zero) * matrix)
        )
    return value


def normalized_operator(
    matrix: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, float]:
    """Compute the Riemann-pulled operator by Theodorsen and Cauchy."""

    boundary, derivative, residual = theodorsen_map(
        matrix,
        N=resolution,
        inflate=0.0,
        iters=400,
        tol=3e-14,
    )
    pullback = GeneralPullback(boundary, derivative)
    resolvents = pullback.resolvent_stack(matrix)
    operator = pullback.calc(pullback.w, matrix, resolvents)
    return operator, float(residual)


def optimal_value(
    matrix: np.ndarray,
    initial_parameters: np.ndarray,
    resolution: int,
) -> tuple[float, np.ndarray, float]:
    """Optimize both degree-two Blaschke zeros locally."""

    operator, map_residual = normalized_operator(matrix, resolution)

    def objective(parameters: np.ndarray) -> float:
        zeros = unpack_zeros(parameters)
        return -float(np.linalg.norm(blaschke_at_matrix(operator, zeros), 2))

    result = minimize(
        objective,
        initial_parameters,
        method="BFGS",
        options={"maxiter": 300, "gtol": 1e-9},
    )
    if not result.success or np.linalg.norm(result.jac) > 2e-5:
        result = minimize(
            objective,
            result.x,
            method="Nelder-Mead",
            options={
                "maxiter": 1500,
                "xatol": 3e-11,
                "fatol": 2e-13,
            },
        )
    return -float(result.fun), np.asarray(result.x), map_residual


def directional_coefficient(
    matrix: np.ndarray,
    direction: np.ndarray,
    step: float,
    initial_parameters: np.ndarray,
    resolution: int,
) -> tuple[float, float]:
    """Return the symmetric second coefficient of the optimized norm."""

    plus, _, plus_residual = optimal_value(
        matrix + step * direction,
        initial_parameters,
        resolution,
    )
    minus, _, minus_residual = optimal_value(
        matrix - step * direction,
        initial_parameters,
        resolution,
    )
    coefficient = (plus + minus - 4) / (2 * step**2)
    return coefficient, max(plus_residual, minus_residual)


def audit_record(
    parameter: float,
    step: float,
    resolution: int,
) -> GauWuNumericHessianRecord:
    """Construct and compare the complete five-normal matrix."""

    matrix = gau_wu_matrix(parameter)
    tangent = equality_tangent(parameter)
    tangent_vectors, singular_values, _ = np.linalg.svd(
        tangent,
        full_matrices=True,
    )
    tangent_rank = int(np.sum(singular_values > 1e-10))
    if tangent_rank != 13:
        raise RuntimeError("unexpected Gau--Wu tangent rank")
    normal = tangent_vectors[:, tangent_rank:]
    complement = slice_basis()
    coordinate_change = np.linalg.solve(
        np.column_stack((tangent, complement)),
        normal,
    )[-5:, :]

    q_value = (1 - np.sqrt(1 - parameter**2)) / parameter
    symbol = sp.symbols("q", real=True)
    slice_hessian = np.array(
        expected_optimized_hessian(symbol).subs(symbol, q_value).evalf(),
        dtype=float,
    )
    exact_hessian = coordinate_change.T @ slice_hessian @ coordinate_change

    initial_parameters = pack_zeros(
        np.array([0, parameter], dtype=complex)
    )
    numerical_hessian = np.zeros((5, 5))
    maximum_map_residual = 0.0
    normal_matrices = [
        real_matrix(normal[:, index]) for index in range(normal.shape[1])
    ]
    for row, direction in enumerate(normal_matrices):
        coefficient, residual = directional_coefficient(
            matrix,
            direction,
            step,
            initial_parameters,
            resolution,
        )
        numerical_hessian[row, row] = coefficient
        maximum_map_residual = max(maximum_map_residual, residual)
    for row in range(5):
        for column in range(row + 1, 5):
            coefficient, residual = directional_coefficient(
                matrix,
                normal_matrices[row] + normal_matrices[column],
                step,
                initial_parameters,
                resolution,
            )
            entry = (
                coefficient
                - numerical_hessian[row, row]
                - numerical_hessian[column, column]
            ) / 2
            numerical_hessian[row, column] = entry
            numerical_hessian[column, row] = entry
            maximum_map_residual = max(maximum_map_residual, residual)

    numerical_eigenvalues = np.linalg.eigvalsh(numerical_hessian)
    exact_eigenvalues = np.linalg.eigvalsh(exact_hessian)
    maximum_hessian_error = float(
        np.max(abs(numerical_hessian - exact_hessian))
    )
    maximum_eigenvalue_error = float(
        np.max(abs(numerical_eigenvalues - exact_eigenvalues))
    )
    checks = (
        maximum_map_residual < 2e-10
        and maximum_hessian_error < 3e-4
        and maximum_eigenvalue_error < 3e-4
        and numerical_eigenvalues[-1] < -1e-3
    )
    if not checks:
        raise RuntimeError(
            "independent numerical Hessian audit failed: "
            f"map={maximum_map_residual}, "
            f"matrix={maximum_hessian_error}, "
            f"eigen={maximum_eigenvalue_error}"
        )
    return GauWuNumericHessianRecord(
        model_parameter=parameter,
        finite_difference_step=step,
        map_resolution=resolution,
        tangent_rank=tangent_rank,
        maximum_map_residual=maximum_map_residual,
        maximum_hessian_error=maximum_hessian_error,
        numerical_eigenvalues=tuple(map(float, numerical_eigenvalues)),
        exact_eigenvalues=tuple(map(float, exact_eigenvalues)),
        maximum_eigenvalue_error=maximum_eigenvalue_error,
        all_checks_passed=True,
    )


def write_record(
    record: GauWuNumericHessianRecord,
    output: Path,
) -> str:
    """Write one JSON record atomically and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(asdict(record), sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter", type=float, default=0.5)
    parser.add_argument("--step", type=float, default=0.003)
    parser.add_argument("--resolution", type=int, default=384)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_transverse_hessian_numeric_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and report the independent finite-difference audit."""

    arguments = parse_args()
    record = audit_record(
        arguments.parameter,
        arguments.step,
        arguments.resolution,
    )
    digest = write_record(record, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": record.all_checks_passed,
                "record_count": 1,
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
