#!/usr/bin/env python3
"""Regenerate one-grade primal/dual Blaschke complementarity.

Let ``B=N/D`` be the corrected finite Blaschke product and ``x`` its
simple top right singular vector.  The model space ``K_B`` is

    {p/D : degree(p) < L}.

Choose the rank-one Stein defect ``q`` orthogonal to

    D(T)^(-1) span{x, Tx, ..., T^(L-1)x}.

The model-kernel identity then makes the von Neumann lower bound
complementary to the Stein metric on ``x``.  This checker constructs
the defect tangent over exact rational power series and verifies that
the primal condition Hessian equals L142--L143's dual coefficient
``-64 c^(2k)`` through the complete one-grade face.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

from crabb_faber_blaschke_formal import (
    FormalDualCalculation,
    formal_dual_calculation,
    invert_constant_amplitude_matrix,
)
from crabb_palindromic_elliptic_hessian import (
    Matrix,
    amplitude_power,
    axis_metric_and_defect,
    coordinate_metric_tangent,
    endpoint_condition_hessian,
    matrix_multiply,
    matrix_scale,
    solve_series_system,
)
from exact_truncated_series import Series


DEFAULT_CASES = (
    (3, 1),
    (4, 1),
    (4, 2),
    (5, 2),
    (6, 3),
    (7, 3),
    (8, 4),
    (9, 4),
)


@dataclass(frozen=True)
class BlaschkeComplementarityRecord:
    """One exact complementary defect audit."""

    length: int
    dimension: int
    grade: int
    central: bool
    target_degree: int
    model_orthogonality_vanishes: bool
    primal_first_nonzero_degree: int
    primal_coefficient: str
    dual_coefficient: str
    primal_dual_gap_vanishes: bool
    defect_tangent_prefix: tuple[tuple[int, str], ...]


def matrix_vector_multiply(
    matrix: Matrix,
    vector: list[Series],
) -> list[Series]:
    """Multiply a series matrix and vector."""

    order = vector[0].order
    return [
        sum(
            (
                matrix[row][column] * vector[column]
                for column in range(len(vector))
            ),
            Series.constant(0, order),
        )
        for row in range(len(matrix))
    ]


def vector_inner(
    left: list[Series],
    right: list[Series],
) -> Series:
    """Return the real formal inner product used by the checker."""

    return sum(
        (
            left_entry * right_entry
            for left_entry, right_entry in zip(
                left,
                right,
                strict=True,
            )
        ),
        Series.constant(0, left[0].order),
    )


def add_vectors(
    left: list[Series],
    right: list[Series],
) -> list[Series]:
    """Add two series vectors."""

    return [
        left_entry + right_entry
        for left_entry, right_entry in zip(left, right, strict=True)
    ]


def complementary_defect_tangent(
    length: int,
    grade: int,
    target_degree: int,
    *,
    calculation_override: FormalDualCalculation | None = None,
    coordinate_tangent_override: Matrix | None = None,
) -> tuple[
    list[Series],
    Matrix,
    list[Series],
    list[Fraction],
    bool,
    Series,
    Series,
]:
    """Construct the model-space complementary defect through the face."""

    if (
        calculation_override is None
    ) != (
        coordinate_tangent_override is None
    ):
        raise ValueError(
            "supply both complementary-calculation overrides together"
        )
    calculation = (
        calculation_override
        if calculation_override is not None
        else formal_dual_calculation(
            length,
            (grade,),
            target_degree,
        )
    )
    order = target_degree + 1
    dimension = length + 1
    denominator_inverse = invert_constant_amplitude_matrix(
        calculation.denominator[0],
        target_degree,
    )
    denominator_inverse_tangent = matrix_scale(
        -1,
        matrix_multiply(
            denominator_inverse,
            matrix_multiply(
                calculation.denominator[1],
                denominator_inverse,
            ),
        ),
    )
    right_base = [
        Series.constant(1 if index == length else 0, order)
        for index in range(dimension)
    ]
    right_tangent = list(calculation.first_coupling)

    orbit_base: list[list[Series]] = []
    orbit_tangent: list[list[Series]] = []
    for degree in range(length):
        power = amplitude_power(calculation.operator, degree)
        powered_base = matrix_vector_multiply(power[0], right_base)
        powered_tangent = add_vectors(
            matrix_vector_multiply(power[1], right_base),
            matrix_vector_multiply(power[0], right_tangent),
        )
        orbit_base.append(
            matrix_vector_multiply(
                denominator_inverse,
                powered_base,
            )
        )
        orbit_tangent.append(
            add_vectors(
                matrix_vector_multiply(
                    denominator_inverse_tangent,
                    powered_base,
                ),
                matrix_vector_multiply(
                    denominator_inverse,
                    powered_tangent,
                ),
            )
        )

    axis_metric, base_defect, coordinate_diagonal = (
        axis_metric_and_defect(
            dimension,
            calculation.operator[0],
            order,
        )
    )
    equations = [*orbit_base]
    equations.append(
        [
            Series.constant(1 if index == 0 else 0, order)
            for index in range(dimension)
        ]
    )
    right_hand_side = [
        -vector_inner(base_defect, tangent)
        for tangent in orbit_tangent
    ]
    right_hand_side.append(Series.constant(0, order))
    defect_tangent = solve_series_system(
        equations,
        right_hand_side,
    )

    orthogonality_vanishes = all(
        vector_inner(base_defect, base).valuation() == order
        and (
            vector_inner(defect_tangent, base)
            + vector_inner(base_defect, tangent)
        ).valuation()
        == order
        for base, tangent in zip(
            orbit_base,
            orbit_tangent,
            strict=True,
        )
    )
    if not orthogonality_vanishes:
        raise AssertionError("model-space complementarity failed")

    if coordinate_tangent_override is None:
        direction = [0] * (length - 1)
        direction[grade - 1] = 1
        direction[length - grade - 1] = 1
        coordinate_tangent = coordinate_metric_tangent(
            dimension,
            direction,
            order,
        )
    else:
        coordinate_tangent = coordinate_tangent_override
    primal_hessian, _ = endpoint_condition_hessian(
        calculation.operator,
        axis_metric,
        base_defect,
        coordinate_diagonal,
        coordinate_tangent,
        defect_tangent,
        order,
    )
    return (
        defect_tangent,
        axis_metric,
        base_defect,
        coordinate_diagonal,
        orthogonality_vanishes,
        primal_hessian,
        calculation.quadratic_eigenvalue,
    )


def make_record(
    length: int,
    grade: int,
) -> BlaschkeComplementarityRecord:
    """Construct and validate one complementary primal/dual face."""

    if not 1 <= grade <= length / 2:
        raise ValueError("require 1 <= k <= L/2")
    target_degree = 2 * grade
    (
        defect_tangent,
        _,
        _,
        _,
        orthogonality_vanishes,
        primal,
        dual,
    ) = complementary_defect_tangent(
        length,
        grade,
        target_degree,
    )
    gap = primal - dual
    if (
        primal.valuation() != target_degree
        or primal.coefficient(target_degree) != -64
        or gap.valuation() != gap.order
    ):
        raise AssertionError("primal/dual face complementarity failed")

    return BlaschkeComplementarityRecord(
        length=length,
        dimension=length + 1,
        grade=grade,
        central=2 * grade == length,
        target_degree=target_degree,
        model_orthogonality_vanishes=orthogonality_vanishes,
        primal_first_nonzero_degree=primal.valuation(),
        primal_coefficient=str(primal.coefficient(target_degree)),
        dual_coefficient=str(dual.coefficient(target_degree)),
        primal_dual_gap_vanishes=gap.valuation() == gap.order,
        defect_tangent_prefix=tuple(
            (
                index,
                str(entry.coefficient(entry.valuation())),
            )
            for index, entry in enumerate(defect_tangent)
            if entry.valuation() < entry.order
        ),
    )


def parse_cases(values: list[str]) -> tuple[tuple[int, int], ...]:
    """Parse repeated ``L,k`` values."""

    if not values:
        return DEFAULT_CASES
    cases = []
    for value in values:
        length, grade = value.split(",", maxsplit=1)
        cases.append((int(length), int(grade)))
    return tuple(cases)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="L,K",
        help="length and one low grade; may be repeated",
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the selected exact complementarity audits."""

    args = parse_args()
    records = [
        make_record(length, grade)
        for length, grade in parse_cases(args.case)
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    print("\n".join(lines), flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
