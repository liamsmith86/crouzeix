#!/usr/bin/env python3
"""Audit the explicit all-size offset-one elliptic Stein certificate.

For the first phase-palindromic disk-equality offset, the defect
tangent needed through first elliptic order is

    2 e_1 + 2 e_(L-1) - 8 c e_3,

where repeated/missing coordinates are handled by the short-chain
convention below.  This script inserts that defect into the exact
truncated Stein recurrence from
``crabb_palindromic_elliptic_hessian.py``.  It verifies the endpoint
metric pair ``(48, 128)``, the condition coefficient ``-64``, and the
stable Stein-forcing pattern used in the all-size proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import (
    Series,
    axis_metric_and_defect,
    coordinate_metric_tangent,
    endpoint_condition_hessian,
    matrix_add,
    matrix_multiply,
    matrix_scale,
    matrix_transpose,
    monomial,
    operator_expansion,
    stein_gramian_expansion,
    zero,
)


SERIES_ORDER = 3
STABLE_MINIMUM_LENGTH = 7


@dataclass(frozen=True)
class OffsetOneRecord:
    dimension: int
    length: int
    defect_tangent: tuple[tuple[int, str, str], ...]
    condition_coefficients: tuple[str, str, str]
    lower_endpoint_metric_coefficient: str
    upper_endpoint_metric_coefficient: str
    endpoint_linear_coupling_valuation: int
    stationarity_residual_valuation: int | None
    quadratic_stein_forcing_diagonal: tuple[str, ...]
    weighted_forcing_sum: str
    stable_forcing_pattern: bool


def coefficient_string(value: Fraction) -> str:
    """Serialize an exact rational coefficient."""

    return str(value)


def explicit_defect_tangent(dimension: int) -> list[Series]:
    """Return the universal offset-one defect tangent through ``c^1``."""

    length = dimension - 1
    tangent = [zero(SERIES_ORDER) for _ in range(dimension)]
    tangent[1] += 2
    if length > 2:
        tangent[length - 1] += 2
    if length >= 3:
        tangent[3] -= 8 * monomial(1, SERIES_ORDER)
    return tangent


def constant_crabb_part(
    operator: tuple[list[list[Series]], ...],
) -> list[list[Series]]:
    """Extract the ``a^0 c^0`` Crabb shift."""

    dimension = len(operator[0])
    return [
        [
            Series.constant(
                operator[0][row][column].coefficient(0),
                SERIES_ORDER,
            )
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]


def stable_forcing_diagonal(length: int) -> tuple[Fraction, ...]:
    """Return the separated-end forcing pattern for ``L>=7``."""

    result = [Fraction(0) for _ in range(length + 1)]
    for index, coefficient in (
        (0, 48),
        (1, -96),
        (2, -32),
        (3, 16),
        (length - 2, 16),
        (length - 1, -32),
        (length, 64),
    ):
        result[index] += Fraction(coefficient)
    return tuple(result)


def make_record(
    dimension: int,
    audit_stationarity: bool,
) -> OffsetOneRecord:
    """Construct and validate one exact offset-one certificate."""

    if dimension < 3:
        raise ValueError("the dimension must be at least three")

    length = dimension - 1
    operator, coefficients = operator_expansion(
        dimension,
        1,
        SERIES_ORDER,
    )
    axis_metric, base_defect, coordinate_diagonal = (
        axis_metric_and_defect(
            dimension,
            operator[0],
            SERIES_ORDER,
        )
    )
    coordinate_tangent = coordinate_metric_tangent(
        dimension,
        coefficients,
        SERIES_ORDER,
    )
    defect_tangent = explicit_defect_tangent(dimension)
    condition_hessian, quadratic_metric = endpoint_condition_hessian(
        operator,
        axis_metric,
        base_defect,
        coordinate_diagonal,
        coordinate_tangent,
        defect_tangent,
        SERIES_ORDER,
    )
    stationarity_residual_valuation = None
    if audit_stationarity:
        gradient = []
        for index in range(1, dimension):
            positive = defect_tangent.copy()
            negative = defect_tangent.copy()
            positive[index] += 1
            negative[index] -= 1
            positive_value = endpoint_condition_hessian(
                operator,
                axis_metric,
                base_defect,
                coordinate_diagonal,
                coordinate_tangent,
                positive,
                SERIES_ORDER,
            )[0]
            negative_value = endpoint_condition_hessian(
                operator,
                axis_metric,
                base_defect,
                coordinate_diagonal,
                coordinate_tangent,
                negative,
                SERIES_ORDER,
            )[0]
            gradient.append((positive_value - negative_value) / 2)
        stationarity_residual_valuation = min(
            entry.valuation()
            for entry in gradient
        )

    full_gramian = stein_gramian_expansion(
        operator,
        base_defect,
        defect_tangent,
    )
    crabb = constant_crabb_part(operator)
    quadratic_forcing = matrix_add(
        full_gramian[2],
        matrix_scale(
            -1,
            matrix_multiply(
                matrix_transpose(crabb),
                matrix_multiply(full_gramian[2], crabb),
            ),
        ),
    )
    forcing_diagonal = tuple(
        quadratic_forcing[index][index].coefficient(2)
        for index in range(dimension)
    )
    weighted_forcing_sum = (
        4 * forcing_diagonal[0] + sum(forcing_diagonal[1:])
    )

    lower_residuals = [
        full_gramian[1][0][index]
        - 2 * coordinate_tangent[0][index]
        for index in range(dimension)
    ]
    upper_residuals = [
        full_gramian[1][length][index]
        - 8 * coordinate_tangent[length][index]
        for index in range(dimension)
    ]
    endpoint_linear_coupling_valuation = min(
        residual.valuation()
        for residual in (*lower_residuals, *upper_residuals)
    )

    condition_coefficients = tuple(
        condition_hessian.coefficient(degree)
        for degree in range(SERIES_ORDER)
    )
    lower_coefficient = quadratic_metric[0][0].coefficient(2)
    upper_coefficient = quadratic_metric[-1][-1].coefficient(2)
    stable_pattern = (
        length < STABLE_MINIMUM_LENGTH
        or forcing_diagonal == stable_forcing_diagonal(length)
    )
    if (
        condition_coefficients != (0, 0, -64)
        or lower_coefficient != 48
        or upper_coefficient != 128
        or endpoint_linear_coupling_valuation < 2
        or (
            stationarity_residual_valuation is not None
            and stationarity_residual_valuation < 2
        )
        or weighted_forcing_sum != 128
        or not stable_pattern
    ):
        raise AssertionError("the offset-one certificate failed")

    serialized_defect = tuple(
        (
            index,
            coefficient_string(entry.coefficient(0)),
            coefficient_string(entry.coefficient(1)),
        )
        for index, entry in enumerate(defect_tangent)
        if entry.valuation() < 2
    )
    return OffsetOneRecord(
        dimension=dimension,
        length=length,
        defect_tangent=serialized_defect,
        condition_coefficients=tuple(
            coefficient_string(value)
            for value in condition_coefficients
        ),
        lower_endpoint_metric_coefficient=coefficient_string(
            lower_coefficient
        ),
        upper_endpoint_metric_coefficient=coefficient_string(
            upper_coefficient
        ),
        endpoint_linear_coupling_valuation=(
            endpoint_linear_coupling_valuation
        ),
        stationarity_residual_valuation=(
            stationarity_residual_valuation
        ),
        quadratic_stein_forcing_diagonal=tuple(
            coefficient_string(value)
            for value in forcing_diagonal
        ),
        weighted_forcing_sum=coefficient_string(
            weighted_forcing_sum
        ),
        stable_forcing_pattern=stable_pattern,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=16)
    parser.add_argument(
        "--stationarity-maximum-size",
        type=int,
        default=8,
        help=(
            "audit the exact defect gradient through this size; "
            "use a value below the minimum size to disable"
        ),
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact dimension grid and optionally persist JSONL."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")

    records = [
        make_record(
            dimension,
            dimension <= args.stationarity_maximum_size,
        )
        for dimension in range(
            args.minimum_size,
            args.maximum_size + 1,
        )
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
