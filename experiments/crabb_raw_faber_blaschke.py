#!/usr/bin/env python3
"""Regenerate the raw disk-flat Faber--Blaschke principal face.

For a raw Toeplitz offset ``j`` in a chain of length ``L``, put
``k=L-j``.  The corrected scalar factor is

    P_L + 2 a (1+c^j) P_k.

This script prepares its finite Blaschke product, evaluates it on the
full coefficient-gauge operator, and constructs L145's canonical
orbit-complement rank-one Stein defect.  Exact rational series verify
that both the dual norm and the feasible primal condition square have
first amplitude-Hessian coefficient

    -64 c^(2k).

It also checks the complete sparse linear Blaschke recurrence through
grade ``k`` and the quadratic top entry ``-16 c^(2k)`` used in the
all-size proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from crabb_disk_flat_elliptic_face import (
    full_disk_operator_expansion,
)
from crabb_faber_blaschke_complementarity import (
    complementary_defect_tangent,
)
from crabb_faber_blaschke_formal import (
    FormalDualCalculation,
    amplitude_polynomial_coefficients,
    add_series_polynomials,
    composed_dickson_polynomials,
    formal_dual_calculation,
    prepare_weierstrass_numerator,
    scale_series_polynomial,
)
from crabb_palindromic_elliptic_hessian import (
    coordinate_metric_tangent,
    zero_matrix,
)
from exact_truncated_series import Series


DEFAULT_CASES = (
    (2, 1),
    (3, 1),
    (3, 2),
    (4, 1),
    (4, 2),
    (4, 3),
    (5, 1),
    (5, 2),
    (5, 3),
    (5, 4),
    (6, 4),
)
DEFAULT_CROSS_CASES = (
    (5, 4, 3),
    (5, 4, 2),
    (5, 3, 2),
    (6, 5, 2),
)


@dataclass(frozen=True)
class RawFaberBlaschkeRecord:
    """One exact raw dual/primal face audit."""

    length: int
    dimension: int
    offset: int
    reflected_grade: int
    target_degree: int
    sparse_linear_recurrence: bool
    linear_top_column_vanishes: bool
    top_right_coupling_vanishes: bool
    quadratic_top_entry: str
    dual_first_nonzero_degree: int
    dual_coefficient: str
    primal_first_nonzero_degree: int
    primal_coefficient: str
    primal_dual_gap_vanishes: bool
    model_orthogonality_vanishes: bool


@dataclass(frozen=True)
class RawFaberBlaschkeCrossRecord:
    """One exact mixed raw-grade audit."""

    length: int
    dimension: int
    offsets: tuple[int, int]
    reflected_grades: tuple[int, int]
    target_degree: int
    dual_cross_vanishes: bool
    primal_cross_vanishes: bool
    combined_model_orthogonality_vanishes: bool
    mixed_primal_dual_gap_vanishes: bool


def raw_prepared_numerator(
    length: int,
    offsets: tuple[int, ...],
    maximum_c_degree: int,
) -> dict:
    """Prepare a sum of corrected raw factors."""

    series_order = maximum_c_degree + 1
    dickson = composed_dickson_polynomials(
        length,
        series_order,
    )
    linear_part = []
    for offset in offsets:
        grade = length - offset
        correction = Series.monomial(offset, series_order)
        linear_part = add_series_polynomials(
            linear_part,
            scale_series_polynomial(
                2,
                add_series_polynomials(
                    dickson[grade],
                    scale_series_polynomial(
                        correction,
                        dickson[grade],
                    ),
                ),
            ),
        )
    prepared_input = amplitude_polynomial_coefficients(
        (dickson[length], linear_part),
        series_order,
    )
    return prepare_weierstrass_numerator(
        length,
        maximum_c_degree,
        prepared_input,
    )


def raw_dual_calculation(
    length: int,
    offsets: tuple[int, ...],
    maximum_c_degree: int,
) -> tuple[FormalDualCalculation, list[list[Series]]]:
    """Return the full-gauge raw dual calculation and coordinate metric."""

    order = maximum_c_degree + 1
    operator, coefficients = full_disk_operator_expansion(
        length,
        offsets,
        order,
    )
    coordinate_tangent = coordinate_metric_tangent(
        length + 1,
        coefficients,
        order,
    )
    calculation = formal_dual_calculation(
        length,
        (),
        maximum_c_degree,
        prepared_override=raw_prepared_numerator(
            length,
            offsets,
            maximum_c_degree,
        ),
        operator_override=operator,
        coordinate_linear_override=coordinate_tangent,
    )
    return calculation, coordinate_tangent


def expected_linear_coefficient(
    length: int,
    offset: int,
    degree: int,
    order: int,
) -> list[list[Series]]:
    """Return the raw sparse Blaschke coefficient at one grade."""

    grade = length - offset
    expected = zero_matrix(length + 1, length + 1, order)
    if degree < grade:
        expected[degree][grade - degree] += Series.constant(
            4,
            order,
        )
    else:
        for row in range(grade, length + 1):
            expected[row][row - grade] += Series.constant(
                4,
                order,
            )
    if degree >= offset:
        expected[degree - offset][
            grade - degree + offset
        ] -= Series.constant(4, order)
    return expected


def make_record(length: int, offset: int) -> RawFaberBlaschkeRecord:
    """Validate one raw offset through its complete principal face."""

    if not 1 <= offset < length:
        raise ValueError("require 1 <= j < L")
    grade = length - offset
    target_degree = 2 * grade
    order = target_degree + 1
    calculation, coordinate_tangent = raw_dual_calculation(
        length,
        (offset,),
        target_degree,
    )

    sparse_recurrence = True
    for degree in range(grade + 1):
        actual = [
            [
                Series.constant(
                    calculation.blaschke[1][row][column].coefficient(
                        degree
                    ),
                    order,
                )
                for column in range(length + 1)
            ]
            for row in range(length + 1)
        ]
        sparse_recurrence &= actual == expected_linear_coefficient(
            length,
            offset,
            degree,
            order,
        )
    top_coupling_vanishes = all(
        entry.valuation() == order
        for entry in calculation.first_coupling
    )
    linear_top_column_vanishes = all(
        calculation.blaschke[1][row][length].valuation() == order
        for row in range(length + 1)
    )
    quadratic_top = calculation.blaschke[2][0][length]
    quadratic_top_ok = (
        all(
            quadratic_top.coefficient(degree) == 0
            for degree in range(target_degree)
        )
        and quadratic_top.coefficient(target_degree) == -16
    )
    dual = calculation.quadratic_eigenvalue
    dual_ok = (
        dual.valuation() == target_degree
        and dual.coefficient(target_degree) == -64
    )
    if not (
        sparse_recurrence
        and linear_top_column_vanishes
        and top_coupling_vanishes
        and quadratic_top_ok
        and dual_ok
    ):
        raise AssertionError("the raw corrected dual face failed")

    (
        primal,
        complementary_dual,
        orthogonality_vanishes,
    ) = complementary_primal(
        length,
        grade,
        target_degree,
        calculation,
        coordinate_tangent,
    )
    gap = primal - complementary_dual
    gap_vanishes = gap.valuation() == order
    if (
        not orthogonality_vanishes
        or primal.valuation() != target_degree
        or primal.coefficient(target_degree) != -64
        or not gap_vanishes
    ):
        raise AssertionError("the raw model complement did not match")

    return RawFaberBlaschkeRecord(
        length=length,
        dimension=length + 1,
        offset=offset,
        reflected_grade=grade,
        target_degree=target_degree,
        sparse_linear_recurrence=sparse_recurrence,
        linear_top_column_vanishes=linear_top_column_vanishes,
        top_right_coupling_vanishes=top_coupling_vanishes,
        quadratic_top_entry=str(
            quadratic_top.coefficient(target_degree)
        ),
        dual_first_nonzero_degree=dual.valuation(),
        dual_coefficient=str(dual.coefficient(target_degree)),
        primal_first_nonzero_degree=primal.valuation(),
        primal_coefficient=str(primal.coefficient(target_degree)),
        primal_dual_gap_vanishes=gap_vanishes,
        model_orthogonality_vanishes=orthogonality_vanishes,
    )


def complementary_primal(
    length: int,
    grade: int,
    target_degree: int,
    calculation: FormalDualCalculation,
    coordinate_tangent: list[list[Series]],
) -> tuple[Series, Series, bool]:
    """Return the canonical primal, dual, and orthogonality flag."""

    (
        _,
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
        calculation_override=calculation,
        coordinate_tangent_override=coordinate_tangent,
    )
    return primal, dual, orthogonality_vanishes


def make_cross_record(
    length: int,
    first_offset: int,
    second_offset: int,
) -> RawFaberBlaschkeCrossRecord:
    """Polarize two distinct raw reflected grades."""

    if (
        not 1 <= first_offset < length
        or not 1 <= second_offset < length
        or first_offset == second_offset
    ):
        raise ValueError("require two distinct raw offsets")
    offsets = (first_offset, second_offset)
    grades = tuple(length - offset for offset in offsets)
    target_degree = sum(grades)

    single_data = [
        raw_dual_calculation(
            length,
            (offset,),
            target_degree,
        )
        for offset in offsets
    ]
    combined_calculation, combined_coordinate = raw_dual_calculation(
        length,
        offsets,
        target_degree,
    )
    dual_cross = combined_calculation.quadratic_eigenvalue
    for calculation, _ in single_data:
        dual_cross -= calculation.quadratic_eigenvalue
    dual_cross_vanishes = all(
        dual_cross.coefficient(degree) == 0
        for degree in range(target_degree + 1)
    )

    single_complements = [
        complementary_primal(
            length,
            grade,
            target_degree,
            calculation,
            coordinate,
        )
        for grade, (calculation, coordinate) in zip(
            grades,
            single_data,
            strict=True,
        )
    ]
    (
        combined_primal,
        combined_dual,
        combined_orthogonality,
    ) = complementary_primal(
        length,
        min(grades),
        target_degree,
        combined_calculation,
        combined_coordinate,
    )
    primal_cross = combined_primal
    for primal, _, _ in single_complements:
        primal_cross -= primal
    primal_cross_vanishes = all(
        primal_cross.coefficient(degree) == 0
        for degree in range(target_degree + 1)
    )
    mixed_gap = combined_primal - combined_dual
    for primal, dual, _ in single_complements:
        mixed_gap -= primal - dual
    mixed_gap_vanishes = (
        mixed_gap.valuation() == target_degree + 1
    )
    if not (
        dual_cross_vanishes
        and primal_cross_vanishes
        and combined_orthogonality
        and mixed_gap_vanishes
    ):
        raise AssertionError("the raw mixed principal face failed")

    return RawFaberBlaschkeCrossRecord(
        length=length,
        dimension=length + 1,
        offsets=offsets,
        reflected_grades=grades,
        target_degree=target_degree,
        dual_cross_vanishes=dual_cross_vanishes,
        primal_cross_vanishes=primal_cross_vanishes,
        combined_model_orthogonality_vanishes=combined_orthogonality,
        mixed_primal_dual_gap_vanishes=mixed_gap_vanishes,
    )


def parse_cases(values: list[str]) -> tuple[tuple[int, int], ...]:
    """Parse repeated ``L,j`` values."""

    if not values:
        return DEFAULT_CASES
    cases = []
    for value in values:
        length, offset = value.split(",", maxsplit=1)
        cases.append((int(length), int(offset)))
    return tuple(cases)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="L,J",
        help="length and raw offset; may be repeated",
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact raw face grid."""

    args = parse_args()
    records = [
        *(
            make_record(length, offset)
            for length, offset in parse_cases(args.case)
        ),
        *(
            make_cross_record(length, first_offset, second_offset)
            for length, first_offset, second_offset in DEFAULT_CROSS_CASES
        ),
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
