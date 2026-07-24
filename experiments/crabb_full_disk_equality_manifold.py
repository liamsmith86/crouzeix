#!/usr/bin/env python3
"""Regenerate the full-disk Hardy equality manifold.

The finite Hardy residual has ``(L-1)^2`` real transverse coordinates.
Its differential at the Crabb point is the diagonal-difference map

    (D Psi(E))_(r,c)
      = E_(c+1,L-1-r) - E_(c,L-2-r).

On Hermitian coefficients its kernel is exactly the Hermitian
Toeplitz space.  This checker verifies that formula and rank, builds
the unique zero-Toeplitz formal equality graph through order five,
and independently solves the nonlinear residual equations.  At the
nonlinear solutions it checks the canonical condition square, the
matching characteristic-Blaschke norm, and full ambient first-order
stationarity after the Riemann correction.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, solve_discrete_lyapunov
from scipy.optimize import least_squares
import sympy as sp

from crabb_circular_normal_quadratic_exact import (
    characteristic_series,
    crabb_normal_directions,
    directional_gradient_series,
    disk_model_from_hermitian_series,
)
from crabb_disk_toeplitz_quartic import extend
from crabb_equality_ambient_stationarity import (
    endpoint_derivatives,
    first_riemann_correction,
)
from crabb_full_disk_base_jet_exact import (
    endpoint_delta_from_disk_series,
)
from crabb_full_disk_hardy_geometry import (
    formal_equality_jet,
    linearized_hardy_residual,
    reversal_matrix,
)
from crabb_full_disk_sixth_hardy_factor import (
    hardy_residual_from_disk_series,
    matrix_is_zero,
)


@dataclass(frozen=True)
class DifferentialRecord:
    """One unrestricted-symbol residual differential audit."""

    record_kind: str
    dimension: int
    length: int
    hermitian_variable_count: int
    residual_real_rank: int
    expected_rank: int
    kernel_dimension: int
    expected_toeplitz_dimension: int
    closed_formula_verified: bool
    reversal_hermitian_verified: bool
    all_identities_verified: bool


@dataclass(frozen=True)
class FormalEqualityRecord:
    """One exact formal equality-graph audit."""

    record_kind: str
    dimension: int
    length: int
    jet_order: int
    residual_zero_through_order: int
    endpoint_zero_through_order: int
    circular_normal_response_zero_through_order: int
    nontrivial_higher_correction: bool
    all_identities_verified: bool


@dataclass(frozen=True)
class NonlinearEqualityRecord:
    """One independent floating nonlinear equality audit."""

    record_kind: str
    dimension: int
    length: int
    solver_success: bool
    residual_norm: str
    transverse_correction_norm: str
    minimum_hermitian_eigenvalue: str
    canonical_condition_square: str
    characteristic_blaschke_norm_square: str
    blaschke_top_vector_overlap: str
    maximum_pulled_ambient_derivative: str
    all_identities_verified: bool


def symbolic_hermitian(
    length: int,
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...]]:
    """Return an unrestricted symbolic Hermitian matrix."""

    matrix = sp.zeros(length)
    variables: list[sp.Symbol] = []
    for row in range(length):
        diagonal = sp.symbols(f"hardy_diagonal_{length}_{row}", real=True)
        matrix[row, row] = diagonal
        variables.append(diagonal)
        for column in range(row + 1, length):
            real, imaginary = sp.symbols(
                f"hardy_real_{length}_{row}_{column} "
                f"hardy_imaginary_{length}_{row}_{column}",
                real=True,
            )
            value = real + sp.I * imaginary
            matrix[row, column] = value
            matrix[column, row] = sp.conjugate(value)
            variables.extend((real, imaginary))
    return matrix, tuple(variables)


def differential_record(length: int) -> DifferentialRecord:
    """Audit the residual differential and its exact real rank."""

    dimension = length + 1
    coefficient, variables = symbolic_hermitian(length)
    hermitian = (
        extend(sp.eye(length) / 2),
        extend(coefficient),
    )
    operator, _ = disk_model_from_hermitian_series(list(hermitian))
    residuals, _ = hardy_residual_from_disk_series(
        hermitian,
        tuple(operator),
    )
    residual = residuals[0]
    predicted = linearized_hardy_residual(coefficient)
    formula_verified = matrix_is_zero(residual - predicted)
    reversed_residual = reversal_matrix(length - 1) * residual
    reversal_hermitian = matrix_is_zero(
        reversed_residual - reversed_residual.conjugate().T
    )
    real_equations = []
    for entry in residual:
        expanded = sp.expand_complex(entry)
        real_equations.extend((sp.re(expanded), sp.im(expanded)))
    rank = sp.Matrix(real_equations).jacobian(variables).rank()
    expected_rank = (length - 1) ** 2
    kernel_dimension = len(variables) - rank
    expected_kernel = 2 * length - 1
    verified = bool(
        formula_verified
        and reversal_hermitian
        and rank == expected_rank
        and kernel_dimension == expected_kernel
    )
    if not verified:
        raise RuntimeError(
            f"Hardy differential audit failed in length {length}"
        )
    return DifferentialRecord(
        record_kind="symbolic_differential",
        dimension=dimension,
        length=length,
        hermitian_variable_count=len(variables),
        residual_real_rank=rank,
        expected_rank=expected_rank,
        kernel_dimension=kernel_dimension,
        expected_toeplitz_dimension=expected_kernel,
        closed_formula_verified=formula_verified,
        reversal_hermitian_verified=reversal_hermitian,
        all_identities_verified=verified,
    )


def formal_equality_record(length: int, order: int) -> FormalEqualityRecord:
    """Build and audit one exact non-palindromic equality jet."""

    direction = (
        sp.Integer(0),
        *[
            sp.Rational(index + 1, 43 * length)
            + sp.I
            * (-1) ** index
            * sp.Rational(index + 2, 47 * length)
            for index in range(1, length)
        ],
    )
    hermitian, metric, operator = formal_equality_jet(direction, order)
    residuals, support_verified = hardy_residual_from_disk_series(
        hermitian,
        operator,
    )
    residual_zero = all(matrix_is_zero(value) for value in residuals)

    dimension = length + 1
    doubled_hermitian = (
        *hermitian,
        *[sp.zeros(dimension) for _ in range(order)],
    )
    doubled_operator, doubled_metric = (
        disk_model_from_hermitian_series(list(doubled_hermitian))
    )
    endpoint = endpoint_delta_from_disk_series(
        doubled_hermitian,
        doubled_metric,
        doubled_operator,
    )
    endpoint_zero = all(value == 0 for value in endpoint)

    characteristic = characteristic_series(list(operator))
    response_zero = True
    for mode in range(3, length + 2):
        for normal in crabb_normal_directions([metric[0]], mode):
            response = directional_gradient_series(
                list(operator),
                list(metric),
                characteristic,
                normal,
            )
            response_zero &= all(
                sp.simplify(value) == 0 for value in response
            )
    nontrivial = any(
        not matrix_is_zero(coefficient)
        for coefficient in hermitian[2:]
    )
    verified = bool(
        support_verified
        and residual_zero
        and endpoint_zero
        and response_zero
        and nontrivial
    )
    if not verified:
        raise RuntimeError("the formal equality-graph audit failed")
    return FormalEqualityRecord(
        record_kind="exact_formal_equality_jet",
        dimension=dimension,
        length=length,
        jet_order=order,
        residual_zero_through_order=order,
        endpoint_zero_through_order=2 * order,
        circular_normal_response_zero_through_order=order,
        nontrivial_higher_correction=nontrivial,
        all_identities_verified=verified,
    )


def unpack_hermitian(values: np.ndarray, length: int) -> np.ndarray:
    """Unpack real coordinates into one complex Hermitian matrix."""

    matrix = np.zeros((length, length), dtype=complex)
    cursor = 0
    for index in range(length):
        matrix[index, index] = values[cursor]
        cursor += 1
    for row in range(length):
        for column in range(row + 1, length):
            value = values[cursor] + 1j * values[cursor + 1]
            matrix[row, column] = value
            matrix[column, row] = np.conjugate(value)
            cursor += 2
    return matrix


def pack_hermitian(matrix: np.ndarray) -> np.ndarray:
    """Pack one complex Hermitian matrix into real coordinates."""

    length = len(matrix)
    values = [float(np.real(matrix[index, index])) for index in range(length)]
    for row in range(length):
        for column in range(row + 1, length):
            values.extend(
                (
                    float(np.real(matrix[row, column])),
                    float(np.imag(matrix[row, column])),
                )
            )
    return np.asarray(values)


def toeplitz_coordinates(matrix: np.ndarray) -> np.ndarray:
    """Return the real Hermitian-Toeplitz projection coordinates."""

    length = len(matrix)
    values = [float(np.real(np.mean(np.diag(matrix))))]
    for offset in range(1, length):
        mean = np.mean(np.diag(matrix, offset))
        values.extend((float(np.real(mean)), float(np.imag(mean))))
    return np.asarray(values)


def numerical_disk_model(
    hermitian: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct the floating coefficient-gauge disk model."""

    length = len(hermitian)
    dimension = length + 1
    extended = np.zeros((dimension, dimension), dtype=complex)
    extended[:length, :length] = hermitian
    shift = np.diag(np.ones(length), 1).astype(complex)
    metric = extended + shift.conj().T @ extended @ shift
    operator = np.linalg.solve(metric, 2 * extended @ shift)
    return extended, metric, operator, shift


def numerical_hardy_residual(hermitian: np.ndarray) -> np.ndarray:
    """Return the finite nonlinear Hardy residual by Markov parameters."""

    length = len(hermitian)
    _, _, operator, shift = numerical_disk_model(hermitian)
    endpoint = np.eye(length + 1, dtype=complex)[:, -1]
    power = np.eye(length + 1, dtype=complex)
    residual = np.zeros((length - 1, length - 1), dtype=complex)
    for degree in range(1, length):
        previous = power
        power = power @ operator
        coefficient = (
            power @ endpoint - shift @ previous @ endpoint
        )
        residual[degree - 1, :] = coefficient[1:length]
    return residual


def pack_hermitian_residual(residual: np.ndarray) -> np.ndarray:
    """Pack ``J residual`` as Hermitian real coordinates."""

    hermitian = np.fliplr(np.eye(len(residual))) @ residual
    values = [
        float(np.real(hermitian[index, index]))
        for index in range(len(hermitian))
    ]
    for row in range(len(hermitian)):
        for column in range(row + 1, len(hermitian)):
            values.extend(
                (
                    float(np.real(hermitian[row, column])),
                    float(np.imag(hermitian[row, column])),
                )
            )
    return np.asarray(values)


def nonlinear_equations(
    values: np.ndarray,
    length: int,
    target_toeplitz: np.ndarray,
) -> np.ndarray:
    """Return residual-zero and fixed-Toeplitz equations."""

    hermitian = unpack_hermitian(values, length)
    return np.concatenate(
        (
            pack_hermitian_residual(
                numerical_hardy_residual(hermitian)
            ),
            toeplitz_coordinates(hermitian) - target_toeplitz,
        )
    )


def matrix_polynomial(
    coefficients: np.ndarray,
    operator: np.ndarray,
) -> np.ndarray:
    """Evaluate a descending-coefficient polynomial."""

    value = np.zeros_like(operator)
    identity = np.eye(len(operator), dtype=complex)
    for coefficient in coefficients:
        value = value @ operator + coefficient * identity
    return value


def format_float(value: float) -> str:
    """Format deterministic compact floating audit output."""

    return f"{value:.12e}"


def nonlinear_equality_record(length: int) -> NonlinearEqualityRecord:
    """Solve and independently audit one nonlinear equality point."""

    initial = np.eye(length, dtype=complex) / 2
    for offset in range(1, length):
        value = (
            0.018 * (offset + 1) / length
            + 1j * 0.011 * (length - offset) / length
        )
        for row in range(length - offset):
            initial[row, row + offset] = value
            initial[row + offset, row] = np.conjugate(value)
    target = toeplitz_coordinates(initial)
    solution = least_squares(
        nonlinear_equations,
        pack_hermitian(initial),
        args=(length, target),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=300,
    )
    hermitian = unpack_hermitian(solution.x, length)
    residual_norm = float(
        np.linalg.norm(numerical_hardy_residual(hermitian))
    )
    correction_norm = float(np.linalg.norm(hermitian - initial))
    minimum_eigenvalue = float(np.linalg.eigvalsh(hermitian)[0])

    extended, metric, operator, _ = numerical_disk_model(hermitian)
    defect = extended[:, 0]
    stein_metric = solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(defect, defect.conj()),
    )
    endpoint_values, endpoint_vectors = eigh(stein_metric, metric)
    condition_square = float(endpoint_values[-1] / endpoint_values[0])

    characteristic = np.poly(operator)
    numerator = characteristic[:-1]
    denominator = np.conjugate(numerator[::-1])
    blaschke = (
        matrix_polynomial(numerator, operator)
        @ np.linalg.inv(matrix_polynomial(denominator, operator))
    )
    blaschke_values, blaschke_vectors = eigh(
        blaschke.conj().T @ metric @ blaschke,
        metric,
    )
    blaschke_norm_square = float(blaschke_values[-1])
    endpoint = np.eye(length + 1, dtype=complex)[:, -1]
    overlap = float(
        abs(
            np.vdot(
                blaschke_vectors[:, -1],
                metric @ endpoint,
            )
        )
        / np.sqrt(np.real(endpoint.conj() @ metric @ endpoint))
    )

    pulled_derivatives = []
    for scalar in (1, 1j):
        for row in range(length + 1):
            for column in range(length + 1):
                perturbation = np.zeros_like(operator)
                perturbation[row, column] = scalar
                pulled = perturbation - first_riemann_correction(
                    operator,
                    metric,
                    perturbation,
                )
                lower, upper = endpoint_derivatives(
                    operator,
                    stein_metric,
                    pulled,
                    endpoint_vectors,
                )
                pulled_derivatives.append(
                    upper / endpoint_values[0]
                    - endpoint_values[-1]
                    * lower
                    / endpoint_values[0] ** 2
                )
    maximum_derivative = float(
        np.max(np.abs(np.asarray(pulled_derivatives)))
    )
    verified = bool(
        solution.success
        and residual_norm < 1e-10
        and correction_norm > 1e-6
        and minimum_eigenvalue > 0
        and abs(condition_square - 4) < 1e-10
        and abs(blaschke_norm_square - 4) < 1e-10
        and abs(overlap - 1) < 1e-10
        and maximum_derivative < 1e-9
    )
    if not verified:
        raise RuntimeError(
            f"nonlinear equality audit failed in length {length}"
        )
    return NonlinearEqualityRecord(
        record_kind="independent_nonlinear_equality",
        dimension=length + 1,
        length=length,
        solver_success=solution.success,
        residual_norm=format_float(residual_norm),
        transverse_correction_norm=format_float(correction_norm),
        minimum_hermitian_eigenvalue=format_float(minimum_eigenvalue),
        canonical_condition_square=format_float(condition_square),
        characteristic_blaschke_norm_square=format_float(
            blaschke_norm_square
        ),
        blaschke_top_vector_overlap=format_float(overlap),
        maximum_pulled_ambient_derivative=format_float(
            maximum_derivative
        ),
        all_identities_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-differential-length", type=int, default=7)
    parser.add_argument("--formal-length", type=int, default=6)
    parser.add_argument("--formal-order", type=int, default=5)
    parser.add_argument("--maximum-numerical-length", type=int, default=6)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and serialize every equality-manifold audit."""

    args = parse_args()
    if args.minimum_length < 3:
        raise ValueError("the minimum length must be at least three")
    records: list[
        DifferentialRecord | FormalEqualityRecord | NonlinearEqualityRecord
    ] = []
    for length in range(
        args.minimum_length,
        args.maximum_differential_length + 1,
    ):
        records.append(differential_record(length))
        print(f"verified Hardy differential length {length}", flush=True)
    records.append(
        formal_equality_record(args.formal_length, args.formal_order)
    )
    print(
        f"verified formal equality jet through order {args.formal_order}",
        flush=True,
    )
    for length in range(
        args.minimum_length,
        args.maximum_numerical_length + 1,
    ):
        records.append(nonlinear_equality_record(length))
        print(f"verified nonlinear equality length {length}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(asdict(record), sort_keys=True) for record in records]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
