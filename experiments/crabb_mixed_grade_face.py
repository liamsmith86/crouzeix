#!/usr/bin/env python3
"""Audit distinct-grade polarization on the first elliptic Newton face.

For two real phase-one coefficient pairs at grades ``k<l``, the
optimized amplitude Hessian is a quadratic form.  This checker
reconstructs it exactly for the two basis directions and their sum,
then polarizes:

    cross = (H(k+l)-H(k)-H(l))/2.

The Hardy/Faber target predicts that every coefficient through the
candidate mixed face ``c^(k+l)`` vanishes.  By default the audit stops
there.  Optional pair and order arguments permit deeper falsification
tests without changing the theorem-level assertion.  Every finite audit
is evidence for the all-size no-alias argument; it is not itself that
proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import cache
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import (
    Series,
    optimized_hessian_from_operator,
    operator_expansion_from_coefficients,
    palindromic_direction,
)


@dataclass(frozen=True)
class MixedGradeFaceRecord:
    dimension: int
    length: int
    first_grade: int
    second_grade: int
    target_degree: int
    audit_order: int
    cross_valuation_lower_bound: int
    cross_coefficients: tuple[str, ...]
    first_face_cross_vanishes: bool


def combined_coefficients(
    length: int,
    grades: tuple[int, ...],
) -> tuple[int, ...]:
    """Return the sum of real phase-one coefficient pairs."""

    result = [0] * (length - 1)
    for grade in grades:
        direction = palindromic_direction(length, grade)
        result = [
            left + right
            for left, right in zip(result, direction, strict=True)
        ]
    return tuple(result)


@cache
def optimized_amplitude_hessian(
    dimension: int,
    coefficients: tuple[int, ...],
    order: int,
) -> Series:
    """Return one exact optimized amplitude Hessian."""

    operator, coefficient_list = operator_expansion_from_coefficients(
        dimension,
        list(coefficients),
        order,
    )
    return optimized_hessian_from_operator(
        operator,
        coefficient_list,
        order,
    )[0]


def make_record(
    dimension: int,
    first_grade: int,
    second_grade: int,
    audit_order: int | None = None,
) -> MixedGradeFaceRecord:
    """Construct and validate one exact mixed-grade polarization."""

    length = dimension - 1
    if not 1 <= first_grade < second_grade <= length // 2:
        raise ValueError("require 1 <= first < second <= floor(L/2)")

    target_degree = first_grade + second_grade
    minimum_order = target_degree + 1
    order = minimum_order if audit_order is None else audit_order
    if order < minimum_order:
        raise ValueError(
            "the audit order must include the first mixed face"
        )
    first_coefficients = combined_coefficients(
        length,
        (first_grade,),
    )
    second_coefficients = combined_coefficients(
        length,
        (second_grade,),
    )
    mixed_coefficients = combined_coefficients(
        length,
        (first_grade, second_grade),
    )
    first = optimized_amplitude_hessian(
        dimension,
        first_coefficients,
        order,
    )
    second = optimized_amplitude_hessian(
        dimension,
        second_coefficients,
        order,
    )
    mixed = optimized_amplitude_hessian(
        dimension,
        mixed_coefficients,
        order,
    )
    cross = (mixed - first - second) / 2
    coefficients = tuple(
        cross.coefficient(degree)
        for degree in range(order)
    )
    first_face_cross_vanishes = not any(
        coefficients[:minimum_order]
    )
    if not first_face_cross_vanishes:
        raise AssertionError(
            "a distinct-grade term appeared on or below the first face"
        )

    return MixedGradeFaceRecord(
        dimension=dimension,
        length=length,
        first_grade=first_grade,
        second_grade=second_grade,
        target_degree=target_degree,
        audit_order=order,
        # A zero truncated series reports ``order`` as its valuation.
        # Record that fact as a lower bound rather than as an equality.
        cross_valuation_lower_bound=cross.valuation(),
        cross_coefficients=tuple(
            str(value)
            for value in coefficients
        ),
        first_face_cross_vanishes=first_face_cross_vanishes,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=5)
    parser.add_argument("--maximum-size", type=int, default=7)
    parser.add_argument(
        "--audit-order",
        type=int,
        help=(
            "optional common truncation order; by default each pair "
            "stops immediately after its first mixed face"
        ),
    )
    parser.add_argument("--first-grade", type=int)
    parser.add_argument("--second-grade", type=int)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def selected_grade_pairs(
    length: int,
    first_grade: int | None,
    second_grade: int | None,
) -> tuple[tuple[int, int], ...]:
    """Return either every independent pair or one validated selection."""

    if (first_grade is None) != (second_grade is None):
        raise ValueError(
            "--first-grade and --second-grade must be supplied together"
        )
    if first_grade is None or second_grade is None:
        return tuple(
            (first, second)
            for first in range(1, length // 2 + 1)
            for second in range(first + 1, length // 2 + 1)
        )
    if not 1 <= first_grade < second_grade <= length // 2:
        raise ValueError(
            "the selected pair is unavailable in a requested size"
        )
    return ((first_grade, second_grade),)


def main() -> None:
    """Run every independent real grade pair in the size range."""

    args = parse_args()
    if args.minimum_size < 5 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 5 <= minimum <= maximum")

    records = []
    for dimension in range(args.minimum_size, args.maximum_size + 1):
        length = dimension - 1
        pairs = selected_grade_pairs(
            length,
            args.first_grade,
            args.second_grade,
        )
        for first_grade, second_grade in pairs:
            record = make_record(
                dimension,
                first_grade,
                second_grade,
                args.audit_order,
            )
            records.append(record)
            print(
                json.dumps(asdict(record), sort_keys=True),
                flush=True,
            )

    if args.output is not None:
        args.output.write_text(
            "".join(
                f"{json.dumps(asdict(record), sort_keys=True)}\n"
                for record in records
            ),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
