#!/usr/bin/env python3
"""Audit the first new coprime grade-four elliptic Hessian.

The general exact Hessian checker originally stopped at grade three.
This focused regression tests ``(L,k)=(9,4)``, the smallest new
grade-four pair that is both noncentral and coprime.  It records the
complete optimizer prefix needed through the candidate face
``c^(2k)`` and compares the gauge-invariant endpoint combination with
the central Dickson prediction.

This is finite exact evidence for principal-face locality.  It is not
an all-size proof of that locality statement.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from math import gcd
from pathlib import Path

from crabb_palindromic_elliptic_hessian import (
    axis_metric_and_defect,
    make_record,
    operator_expansion,
)


DEFAULT_LENGTH = 9
DEFAULT_GRADE = 4
DEFAULT_ORDER = 10


@dataclass(frozen=True)
class PrincipalFaceLocalityRecord:
    """Serializable focused grade-four audit."""

    length: int
    dimension: int
    grade: int
    greatest_common_divisor: int
    central_length: int
    central_prediction: int
    first_nonzero_degree: int
    leading_coefficient: int
    lower_endpoint_metric_coefficient: str
    upper_endpoint_metric_coefficient: str
    endpoint_invariant: int
    transport_difference: tuple[tuple[int, int, str], ...]
    correction_disk_energy: str
    transported_face_coefficient: str
    optimized_defect_prefix: tuple[tuple[str, ...], ...]


def transported_defect_prefix(
    length: int,
    grade: int,
    order: int,
) -> tuple[tuple[Fraction, ...], ...]:
    """Return ``2 U(u)d_c`` through the recorded face order."""

    operator, coefficients = operator_expansion(
        length + 1,
        grade,
        order,
    )
    _, axis_defect, _ = axis_metric_and_defect(
        length + 1,
        operator[0],
        order,
    )
    return tuple(
        tuple(
            sum(
                2
                * coefficient
                * axis_defect[coordinate - coefficient_grade].coefficient(
                    degree
                )
                for coefficient_grade, coefficient in enumerate(
                    coefficients,
                    start=1,
                )
                if coefficient and coordinate >= coefficient_grade
            )
            for degree in range(grade + 1)
        )
        for coordinate in range(1, length + 1)
    )


def make_locality_record(
    length: int = DEFAULT_LENGTH,
    grade: int = DEFAULT_GRADE,
    order: int = DEFAULT_ORDER,
) -> PrincipalFaceLocalityRecord:
    """Run and validate one exact noncentral locality audit."""

    if length <= 2 * grade:
        raise ValueError("the focused audit requires a noncentral grade")
    base = make_record(length + 1, grade, order)
    lower = int(base.lower_endpoint_metric_coefficient)
    upper = int(base.upper_endpoint_metric_coefficient)
    endpoint_invariant = upper - 4 * lower
    central_prediction = -64
    transported = transported_defect_prefix(length, grade, order)
    transport_difference = tuple(
        (
            coordinate,
            degree,
            str(
                Fraction(base.optimized_defect_prefix[coordinate - 1][degree])
                - transported[coordinate - 1][degree]
            ),
        )
        for coordinate in range(1, length + 1)
        for degree in range(grade + 1)
        if (
            Fraction(base.optimized_defect_prefix[coordinate - 1][degree])
            != transported[coordinate - 1][degree]
        )
    )
    if any(degree < grade for _, degree, _ in transport_difference):
        raise AssertionError(
            "the transported defect failed below the reflected grade"
        )
    correction_disk_energy = sum(
        (
            Fraction(4)
            if coordinate < length
            else Fraction(8, 3)
        )
        * Fraction(coefficient) ** 2
        for coordinate, degree, coefficient in transport_difference
        if degree == grade
    )
    transported_face = (
        Fraction(base.leading_coefficient) + correction_disk_energy
    )
    if (
        base.leading_coefficient != central_prediction
        or endpoint_invariant != central_prediction
    ):
        raise AssertionError(
            "the noncentral face disagreed with the central prediction"
        )
    return PrincipalFaceLocalityRecord(
        length=length,
        dimension=length + 1,
        grade=grade,
        greatest_common_divisor=gcd(length, grade),
        central_length=2 * grade,
        central_prediction=central_prediction,
        first_nonzero_degree=base.first_nonzero_degree,
        leading_coefficient=base.leading_coefficient,
        lower_endpoint_metric_coefficient=(
            base.lower_endpoint_metric_coefficient
        ),
        upper_endpoint_metric_coefficient=(
            base.upper_endpoint_metric_coefficient
        ),
        endpoint_invariant=endpoint_invariant,
        transport_difference=transport_difference,
        correction_disk_energy=str(correction_disk_energy),
        transported_face_coefficient=str(transported_face),
        optimized_defect_prefix=base.optimized_defect_prefix,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=DEFAULT_LENGTH)
    parser.add_argument("--grade", type=int, default=DEFAULT_GRADE)
    parser.add_argument("--order", type=int, default=DEFAULT_ORDER)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the focused exact audit and optionally persist it."""

    args = parse_args()
    record = make_locality_record(args.length, args.grade, args.order)
    line = json.dumps(asdict(record), sort_keys=True)
    print(line, flush=True)
    if args.output is not None:
        args.output.write_text(f"{line}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
