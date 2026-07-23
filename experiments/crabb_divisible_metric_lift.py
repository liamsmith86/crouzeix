#!/usr/bin/env python3
"""Audit the exact metric lift for divisible Crabb grades.

This is an independent floating regression for the analytic theorem in
``proof/crabb_divisible_dickson_descent.md``.  It constructs the full
physical operator, the degree-``k`` Chebyshev--Blaschke product, the
outer-critical multiplier from L130, and the explicit lifted rank-one
defect.  No full-dimensional defect optimization is used.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, solve_discrete_lyapunov

from blaschke_stein_composition import model_functions
from crabb_central_dual_lift import (
    critical_multiplier,
    normalized_residual,
)
from crabb_palindromic_elliptic_face import (
    disk_matrix,
    ellipse_pullback,
    rank_one_envelope,
)
from crabb_touching_gradient import chebyshev_blaschke_zeros


DEFAULT_DEGREES = (2, 3)
DEFAULT_QUOTIENTS = (3, 4)
DEFAULT_PARAMETERS = (0.2, 0.35)
DEFAULT_AMPLITUDE = 0.04


@dataclass(frozen=True)
class DivisibleMetricLiftRecord:
    degree: int
    quotient: int
    dimension: int
    outer_dimension: int
    ellipse_parameter: float
    descended_parameter: float
    amplitude: float
    outer_image_residual: float
    critical_factorization_residual: float
    model_multiplier_residual: float
    scalar_compression_residual: float
    metric_cross_residual: float
    outer_metric_residual: float
    constructed_condition: float
    outer_condition: float
    condition_gap: float
    inner_minimum_margin: float
    inner_maximum_margin: float


def palindromic_coefficients(
    length: int,
    first_grade: int,
) -> np.ndarray:
    """Return the real phase-one coefficient pair."""

    coefficients = np.zeros(length - 1)
    coefficients[first_grade - 1] = 1
    coefficients[length - first_grade - 1] = 1
    return coefficients


def coordinate_gramian(
    length: int,
    first_grade: int,
    amplitude: float,
) -> np.ndarray:
    """Return the real L123 coordinate Gramian."""

    dimension = length + 1
    toeplitz = np.zeros((dimension, dimension))
    toeplitz[:length, :length] = np.eye(length) / 2
    coefficients = palindromic_coefficients(length, first_grade)
    for grade, coefficient in enumerate(coefficients, start=1):
        for row in range(length - grade):
            toeplitz[row, row + grade] = amplitude * coefficient
            toeplitz[row + grade, row] = amplitude * coefficient

    shift = np.zeros((dimension, dimension))
    for column in range(1, dimension):
        shift[column - 1, column] = 1
    return toeplitz + shift.T @ toeplitz @ shift


def matrix_power(matrix: np.ndarray, exponent: float) -> np.ndarray:
    """Return a real symmetric matrix power."""

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (
        eigenvectors * eigenvalues**exponent
    ) @ eigenvectors.T


def physical_operator(
    length: int,
    first_grade: int,
    amplitude: float,
    ellipse_parameter: float,
) -> np.ndarray:
    """Return the normalized physical equality/ellipse operator."""

    coefficients = palindromic_coefficients(length, first_grade)
    disk_operator = disk_matrix(
        length + 1,
        amplitude,
        coefficients,
    )
    return np.asarray(
        ellipse_pullback(disk_operator, ellipse_parameter),
        dtype=complex,
    )


def optimized_outer_metric(
    operator: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return condition, defect, and rank-one outer metric."""

    condition, tail = rank_one_envelope(
        np.real_if_close(operator).real
    )
    defect = np.concatenate(([1.0], tail)).astype(complex)
    metric = solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(defect, defect.conj()),
    )
    return condition, defect, metric


def make_record(
    degree: int,
    quotient: int,
    ellipse_parameter: float,
    amplitude: float,
) -> DivisibleMetricLiftRecord:
    """Construct and audit one explicit divisible metric lift."""

    length = degree * quotient
    dimension = length + 1
    outer_dimension = quotient + 1
    descended_parameter = ellipse_parameter**degree
    operator = physical_operator(
        length,
        degree,
        amplitude,
        ellipse_parameter,
    )
    outer_operator = physical_operator(
        quotient,
        1,
        amplitude,
        descended_parameter,
    )
    outer_condition, outer_defect, outer_metric = (
        optimized_outer_metric(outer_operator)
    )

    zeros = tuple(
        complex(zero)
        for zero in chebyshev_blaschke_zeros(
            degree + 1,
            ellipse_parameter,
        )
    )
    functions, descended_operator = model_functions(operator, zeros)

    full_gramian = coordinate_gramian(
        length,
        degree,
        amplitude,
    )
    outer_gramian = coordinate_gramian(
        quotient,
        1,
        amplitude,
    )
    outer_indices = [
        residue * degree
        for residue in range(quotient + 1)
    ]
    embedding = np.eye(dimension)[:, outer_indices]
    outer_map = (
        matrix_power(full_gramian, 0.5)
        @ embedding
        @ matrix_power(outer_gramian, -0.5)
    ).astype(complex)
    if (
        np.linalg.norm(
            outer_map.conj().T @ outer_map
            - np.eye(outer_dimension),
            2,
        )
        > 3e-13
    ):
        raise AssertionError("the outer map is not isometric")

    outer_image = (
        outer_map.conj().T @ descended_operator @ outer_map
    )
    phase = np.vdot(outer_operator, outer_image) / np.vdot(
        outer_operator,
        outer_operator,
    )
    outer_image_residual = normalized_residual(
        outer_image - phase * outer_operator,
        np.linalg.norm(outer_operator, 2),
    )

    (
        multiplier,
        multiplier_coefficients,
        critical_factorization_residual,
        model_multiplier_residual,
    ) = critical_multiplier(operator, zeros, functions)
    multiplier_norm_square = float(
        np.vdot(
            multiplier_coefficients,
            multiplier_coefficients,
        ).real
    )
    explicit_lift = np.linalg.solve(
        multiplier.conj().T,
        outer_map,
    )
    scalar_compression_residual = max(
        normalized_residual(
            outer_map.conj().T
            @ function.conj().T
            @ explicit_lift
            - (
                np.conjugate(coefficient)
                / multiplier_norm_square
            )
            * np.eye(outer_dimension),
        )
        for function, coefficient in zip(
            functions,
            multiplier_coefficients,
        )
    )

    full_defect = explicit_lift @ outer_defect
    full_metric = solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(full_defect, full_defect.conj()),
    )
    inner_map = null_space(outer_map.conj().T)
    metric_cross_residual = normalized_residual(
        outer_map.conj().T @ full_metric @ inner_map,
        np.linalg.norm(full_metric, 2),
    )
    outer_compression = (
        outer_map.conj().T @ full_metric @ outer_map
    )
    outer_scale = np.vdot(
        outer_metric,
        outer_compression,
    ) / np.vdot(outer_metric, outer_metric)
    outer_metric_residual = normalized_residual(
        outer_compression - outer_scale * outer_metric,
        np.linalg.norm(outer_compression, 2),
    )

    full_values = np.linalg.eigvalsh(full_metric)
    constructed_condition = float(
        full_values[-1] / full_values[0]
    )
    scaled_outer_values = np.linalg.eigvalsh(
        outer_scale * outer_metric
    )
    inner_values = np.linalg.eigvalsh(
        inner_map.conj().T @ full_metric @ inner_map
    )
    inner_minimum_margin = float(
        inner_values[0] - scaled_outer_values[0]
    )
    inner_maximum_margin = float(
        scaled_outer_values[-1] - inner_values[-1]
    )

    record = DivisibleMetricLiftRecord(
        degree=degree,
        quotient=quotient,
        dimension=dimension,
        outer_dimension=outer_dimension,
        ellipse_parameter=ellipse_parameter,
        descended_parameter=descended_parameter,
        amplitude=amplitude,
        outer_image_residual=outer_image_residual,
        critical_factorization_residual=(
            critical_factorization_residual
        ),
        model_multiplier_residual=model_multiplier_residual,
        scalar_compression_residual=scalar_compression_residual,
        metric_cross_residual=metric_cross_residual,
        outer_metric_residual=outer_metric_residual,
        constructed_condition=constructed_condition,
        outer_condition=outer_condition,
        condition_gap=constructed_condition - outer_condition,
        inner_minimum_margin=inner_minimum_margin,
        inner_maximum_margin=inner_maximum_margin,
    )
    if (
        max(
            outer_image_residual,
            critical_factorization_residual,
            model_multiplier_residual,
            scalar_compression_residual,
            metric_cross_residual,
            outer_metric_residual,
            abs(record.condition_gap),
        )
        > 2e-9
        or inner_minimum_margin <= 0
        or inner_maximum_margin <= 0
    ):
        raise AssertionError(
            f"the divisible metric lift failed: {record}"
        )
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--degrees",
        type=int,
        nargs="+",
        default=DEFAULT_DEGREES,
    )
    parser.add_argument(
        "--quotients",
        type=int,
        nargs="+",
        default=DEFAULT_QUOTIENTS,
    )
    parser.add_argument(
        "--ellipse-parameters",
        type=float,
        nargs="+",
        default=DEFAULT_PARAMETERS,
    )
    parser.add_argument(
        "--amplitude",
        type=float,
        default=DEFAULT_AMPLITUDE,
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic grid and optionally persist JSONL."""

    args = parse_args()
    if (
        any(degree < 1 for degree in args.degrees)
        or any(quotient < 2 for quotient in args.quotients)
    ):
        raise ValueError("invalid degree or quotient")

    records = [
        make_record(
            degree,
            quotient,
            ellipse_parameter,
            args.amplitude,
        )
        for degree in args.degrees
        for quotient in args.quotients
        for ellipse_parameter in args.ellipse_parameters
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
