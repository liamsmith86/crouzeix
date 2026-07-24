#!/usr/bin/env python3
"""Falsify exact disk-model complementarity away from equality.

For a small Toeplitz disk-chart point, construct

* the coefficient-gauge operator ``A`` and metric ``K``;
* its characteristic finite Blaschke product ``B=g/g#``; and
* L145's orbit-complement rank-one Stein metric.

The model identity makes the upper metric complementary to ``B(A)``
on one endpoint line.  It does *not* force equality of the full
condition number and ``||B(A)||_K^2``.  This checker records:

1. an exact rational size-four counterexample; and
2. a binary64 grid showing the small but positive model gap.

The exact counterexample, rather than a floating-point tolerance,
is the load-bearing falsification.  A second exact record audits the
failure of a proposed two-sided comparison between raw and prepared
reflected grades.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh, null_space, solve_discrete_lyapunov
import sympy as sp

from crabb_disk_exact import (
    characteristic_factor,
    coefficient_model as exact_coefficient_model,
    evaluate_polynomial as evaluate_exact_polynomial,
    solve_symmetric_stein,
)

DEFAULT_SEED = 70_223


@dataclass(frozen=True)
class ExactCounterexampleRecord:
    """Exact obstruction to full model complementarity."""

    dimension: int
    coefficients: tuple[str, ...]
    blaschke_rank: int
    model_kernel_annihilates_top_line: bool
    lower_endpoint_is_generalized_eigenvector: bool
    lower_endpoint_residual: tuple[str, ...]
    blaschke_norm_square: str


@dataclass(frozen=True)
class DiskModelGapRecord:
    """Floating-point size scan of the strict model gap."""

    length: int
    dimension: int
    sample_count: int
    minimum_positive_gap: float
    maximum_gap: float
    maximum_condition_square: float
    maximum_second_singular_value: float
    minimum_coordinate_eigenvalue: float


@dataclass(frozen=True)
class ReesAliasRecord:
    """Exact nonlinear high-to-low characteristic-grade alias."""

    length: int
    dimension: int
    raw_offset: int
    raw_reflected_grade: int
    aliased_characteristic_grade: int
    characteristic_coefficient: str
    physical_scaling: str
    raw_norm_order: int
    prepared_norm_order: int
    two_sided_comparison_holds: bool


def matrix_polynomial(
    coefficients: np.ndarray,
    matrix: np.ndarray,
) -> np.ndarray:
    """Evaluate descending scalar coefficients by Horner's rule."""

    result = np.zeros_like(matrix, dtype=complex)
    identity = np.eye(matrix.shape[0], dtype=complex)
    for coefficient in coefficients:
        result = result @ matrix + coefficient * identity
    return result


def coefficient_model(
    length: int,
    coefficients: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``A,K`` for one Hermitian Toeplitz disk-chart point."""

    dimension = length + 1
    toeplitz = np.zeros((dimension, dimension), dtype=complex)
    toeplitz[:length, :length] = np.eye(length) / 2
    for offset, coefficient in enumerate(coefficients, start=1):
        for row in range(length - offset):
            toeplitz[row, row + offset] = coefficient
            toeplitz[row + offset, row] = coefficient.conjugate()
    shift = np.zeros((dimension, dimension), dtype=complex)
    for row in range(length):
        shift[row, row + 1] = 1
    coordinate = toeplitz + shift.conjugate().T @ toeplitz @ shift
    operator = 2 * np.linalg.solve(coordinate, toeplitz @ shift)
    return operator, coordinate


def generalized_blaschke_data(
    operator: np.ndarray,
    coordinate: np.ndarray,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Return the dual square, second singular value, vector, denominator."""

    characteristic = np.poly(operator)
    factor = characteristic[:-1]
    denominator_coefficients = factor[::-1].conjugate()
    numerator = matrix_polynomial(factor, operator)
    denominator = matrix_polynomial(
        denominator_coefficients,
        operator,
    )
    blaschke = numerator @ np.linalg.inv(denominator)

    eigenvalues, eigenvectors = np.linalg.eigh(coordinate)
    square_root = (
        eigenvectors * np.sqrt(eigenvalues)
    ) @ eigenvectors.conjugate().T
    inverse_square_root = (
        eigenvectors * (1 / np.sqrt(eigenvalues))
    ) @ eigenvectors.conjugate().T
    physical_blaschke = (
        square_root @ blaschke @ inverse_square_root
    )
    _, singular_values, right_adjoint = np.linalg.svd(
        physical_blaschke
    )
    physical_right = right_adjoint.conjugate().T[:, 0]
    coefficient_right = inverse_square_root @ physical_right
    return (
        float(singular_values[0] ** 2),
        float(singular_values[1]),
        coefficient_right,
        denominator,
    )


def model_condition_square(
    operator: np.ndarray,
    coordinate: np.ndarray,
    right_vector: np.ndarray,
    denominator: np.ndarray,
) -> float:
    """Construct and evaluate the orbit-complement Stein metric."""

    length = operator.shape[0] - 1
    denominator_inverse = np.linalg.inv(denominator)
    orbit = [
        denominator_inverse @ np.linalg.matrix_power(operator, power)
        @ right_vector
        for power in range(length)
    ]
    orbit_matrix = np.column_stack(orbit)
    defect_space = null_space(orbit_matrix.conjugate().T)
    if defect_space.shape[1] != 1:
        raise AssertionError("the model orbit lost codimension one")
    defect = defect_space[:, 0]
    metric = solve_discrete_lyapunov(
        operator.conjugate().T,
        np.outer(defect, defect.conjugate()),
    )
    levels = eigvalsh(metric, coordinate)
    return float(levels[-1] / levels[0])


def exact_counterexample() -> ExactCounterexampleRecord:
    """Return the rational ``L=3`` counterexample used in the proof note."""

    length = 3
    dimension = length + 1
    coefficients = (sp.Rational(1, 20), sp.Rational(1, 30))
    operator, coordinate, _, _ = exact_coefficient_model(coefficients)

    variable = sp.symbols("xi")
    factor = sp.Poly(
        characteristic_factor(operator, variable),
        variable,
    )
    numerator_coefficients = factor.all_coeffs()
    denominator_coefficients = list(reversed(numerator_coefficients))
    denominator = evaluate_exact_polynomial(
        list(reversed(denominator_coefficients)),
        operator,
    )
    numerator = evaluate_exact_polynomial(
        list(reversed(numerator_coefficients)),
        operator,
    )
    blaschke = sp.simplify(numerator * denominator.inv())
    top_right = sp.eye(dimension)[:, -1]
    orbit = sp.Matrix.hstack(
        *[
            denominator.inv() * operator**power * top_right
            for power in range(length)
        ]
    )
    defect = orbit.T.nullspace()[0]
    metric = solve_symmetric_stein(operator, defect)

    lower = sp.eye(dimension)[:, 0]
    lower_level = sp.simplify(
        (lower.T * metric * lower)[0]
        / (lower.T * coordinate * lower)[0]
    )
    lower_residual = sp.simplify(
        metric * lower - lower_level * coordinate * lower
    )
    model_kernel = sp.simplify(
        metric - blaschke.T * metric * blaschke
    )
    norm_square = next(
        eigenvalue
        for eigenvalue in (
            coordinate.inv()
            * blaschke.T
            * coordinate
            * blaschke
        ).eigenvals()
        if eigenvalue != 0
    )

    if (
        blaschke.rank() != 1
        or model_kernel * top_right != sp.zeros(dimension, 1)
        or lower_residual == sp.zeros(dimension, 1)
    ):
        raise AssertionError("the exact model-gap counterexample failed")
    return ExactCounterexampleRecord(
        dimension=dimension,
        coefficients=tuple(map(str, coefficients)),
        blaschke_rank=blaschke.rank(),
        model_kernel_annihilates_top_line=True,
        lower_endpoint_is_generalized_eigenvector=False,
        lower_endpoint_residual=tuple(map(str, lower_residual)),
        blaschke_norm_square=str(norm_square),
    )


def exact_rees_alias() -> ReesAliasRecord:
    """Return the exact ``L=6`` nonlinear grade-alias obstruction."""

    length = 6
    parameter = sp.symbols("t", real=True)
    coefficients = (
        parameter,
        *[sp.Integer(0) for _ in range(length - 2)],
    )
    operator, _, _, _ = exact_coefficient_model(coefficients)
    variable = sp.symbols("xi")
    factor = sp.Poly(
        characteristic_factor(operator, variable),
        variable,
    )
    aliased_coefficient = sp.factor(factor.nth(3))
    expected = (
        32
        * parameter**7
        / (
            (6 * parameter**2 - 1)
            * (16 * parameter**4 - 14 * parameter**2 + 1)
        )
    )
    if sp.cancel(aliased_coefficient - expected) != 0:
        raise AssertionError("the exact Rees grade alias changed")
    return ReesAliasRecord(
        length=length,
        dimension=length + 1,
        raw_offset=1,
        raw_reflected_grade=5,
        aliased_characteristic_grade=3,
        characteristic_coefficient=str(aliased_coefficient),
        physical_scaling="c=t^5",
        raw_norm_order=52,
        prepared_norm_order=44,
        two_sided_comparison_holds=False,
    )


def make_record(
    length: int,
    sample_count: int,
    rng: np.random.Generator,
) -> DiskModelGapRecord:
    """Run one dimension's deterministic finite-amplitude sample grid."""

    positive_gaps: list[float] = []
    maximum_condition = 0.0
    maximum_second_singular = 0.0
    minimum_coordinate = float("inf")
    for _ in range(sample_count):
        direction = (
            rng.normal(size=length - 1)
            + 1j * rng.normal(size=length - 1)
        )
        direction /= np.linalg.norm(direction)
        amplitude = rng.uniform(0.12, 0.22)
        operator, coordinate = coefficient_model(
            length,
            amplitude * direction,
        )
        minimum_coordinate = min(
            minimum_coordinate,
            float(np.linalg.eigvalsh(coordinate)[0]),
        )
        (
            dual_square,
            second_singular,
            right_vector,
            denominator,
        ) = generalized_blaschke_data(operator, coordinate)
        condition_square = model_condition_square(
            operator,
            coordinate,
            right_vector,
            denominator,
        )
        gap = condition_square - dual_square
        if gap > 1e-11:
            positive_gaps.append(gap)
        maximum_condition = max(maximum_condition, condition_square)
        maximum_second_singular = max(
            maximum_second_singular,
            second_singular,
        )

    if (
        not positive_gaps
        or maximum_condition > 4 + 5e-6
        or maximum_second_singular > 5e-10
        or minimum_coordinate <= 0
    ):
        raise AssertionError("the disk model-gap probe failed")
    return DiskModelGapRecord(
        length=length,
        dimension=length + 1,
        sample_count=sample_count,
        minimum_positive_gap=min(positive_gaps),
        maximum_gap=max(positive_gaps),
        maximum_condition_square=maximum_condition,
        maximum_second_singular_value=maximum_second_singular,
        minimum_coordinate_eigenvalue=minimum_coordinate,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-length", type=int, default=8)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run and optionally persist the exact and floating audits."""

    args = parse_args()
    rng = np.random.default_rng(args.seed)
    records: list[
        ExactCounterexampleRecord | ReesAliasRecord | DiskModelGapRecord
    ] = [
        exact_counterexample(),
        exact_rees_alias(),
        *[
            make_record(length, args.samples, rng)
            for length in range(3, args.maximum_length + 1)
        ],
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    print("\n".join(lines))
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
