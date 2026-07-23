#!/usr/bin/env python3
"""Audit divisible Dickson descent for a nonreal coefficient phase.

L132's published checker uses the phase-one pair ``(1,1)``.  The proof
complexifies by replacing it with ``(u, conjugate(u))``.  This
regression uses the independent representative ``u=i`` and verifies
the exact outer pencil, polynomial cross-block cancellation,
subcritical compression, and Hermitian coordinate-Gramian reduction.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_central_dickson_descent import coefficient_crabb
from crabb_divisible_dickson_descent import (
    dickson_values,
    reversal,
)


@dataclass(frozen=True)
class ComplexDivisibleDescentRecord:
    degree: int
    quotient: int
    dimension: int
    exact_outer_pencil: bool
    polynomial_cross_blocks_vanish: bool
    subcritical_compression_vanishes: bool
    coordinate_gramian_reduces: bool
    exact_outer_coordinate_gramian: bool


def complex_coefficients(
    length: int,
    first_grade: int,
    phase: sp.Expr,
) -> list[sp.Expr]:
    """Return the one-pair coefficients ``(u, conjugate(u))``."""

    coefficients: list[sp.Expr] = [sp.Integer(0)] * (length - 1)
    coefficients[first_grade - 1] = phase
    coefficients[length - first_grade - 1] = sp.conjugate(phase)
    return coefficients


def complex_equality_pencil(
    length: int,
    first_grade: int,
    parameter: sp.Expr,
    amplitude: sp.Expr,
    phase: sp.Expr,
) -> sp.Matrix:
    """Return the complex one-pair equality/ellipse pencil."""

    dimension = length + 1
    flip = reversal(dimension)
    crabb = coefficient_crabb(dimension)
    axis = crabb + parameter * flip * crabb * flip
    disk_direction = sp.zeros(dimension)
    coefficients = complex_coefficients(
        length,
        first_grade,
        phase,
    )
    for grade, coefficient in enumerate(coefficients, start=1):
        disk_direction[0, grade + 1] = 2 * coefficient
        disk_direction[length, grade + 1] = -2 * coefficient
    direction = (
        disk_direction
        + parameter
        * flip
        * sp.conjugate(disk_direction)
        * flip
    )
    return axis + amplitude * direction


def complex_coordinate_gramian(
    length: int,
    first_grade: int,
    amplitude: sp.Expr,
    phase: sp.Expr,
) -> sp.Matrix:
    """Return the Hermitian Toeplitz coordinate Gramian."""

    dimension = length + 1
    coefficients = complex_coefficients(
        length,
        first_grade,
        phase,
    )
    toeplitz = sp.zeros(dimension)
    for index in range(length):
        toeplitz[index, index] = sp.Rational(1, 2)
    for grade, coefficient in enumerate(coefficients, start=1):
        for row in range(length - grade):
            toeplitz[row, row + grade] = amplitude * coefficient
            toeplitz[row + grade, row] = (
                amplitude * sp.conjugate(coefficient)
            )

    shift = sp.zeros(dimension)
    for column in range(1, dimension):
        shift[column - 1, column] = 1
    return sp.expand(toeplitz + shift.T * toeplitz * shift)


def make_record(
    degree: int,
    quotient: int,
) -> ComplexDivisibleDescentRecord:
    """Construct and validate one exact complex descent."""

    if degree < 1 or quotient < 3:
        raise ValueError("require degree >= 1 and quotient >= 3")

    parameter, amplitude = sp.symbols("c a", real=True)
    phase = sp.I
    length = degree * quotient
    dimension = length + 1
    outer_dimension = quotient + 1
    pencil = complex_equality_pencil(
        length,
        degree,
        parameter,
        amplitude,
        phase,
    )
    value = dickson_values(pencil, parameter, degree)[-1]
    expected_outer = complex_equality_pencil(
        quotient,
        1,
        parameter**degree,
        amplitude,
        phase,
    )
    outer_indices = [
        residue * degree
        for residue in range(quotient + 1)
    ]
    inner_indices = [
        index
        for index in range(dimension)
        if index not in outer_indices
    ]

    full_gramian = complex_coordinate_gramian(
        length,
        degree,
        amplitude,
        phase,
    )
    expected_outer_gramian = complex_coordinate_gramian(
        quotient,
        1,
        amplitude,
        phase,
    )
    subcritical_values = dickson_values(
        pencil,
        parameter,
        max(0, degree - 1),
    )

    record = ComplexDivisibleDescentRecord(
        degree=degree,
        quotient=quotient,
        dimension=dimension,
        exact_outer_pencil=(
            sp.expand(
                value.extract(outer_indices, outer_indices)
                - expected_outer
            )
            == sp.zeros(outer_dimension)
        ),
        polynomial_cross_blocks_vanish=(
            value.extract(outer_indices, inner_indices)
            == sp.zeros(outer_dimension, len(inner_indices))
            and value.extract(inner_indices, outer_indices)
            == sp.zeros(len(inner_indices), outer_dimension)
        ),
        subcritical_compression_vanishes=all(
            subcritical.extract(outer_indices, outer_indices)
            == sp.zeros(outer_dimension)
            for subcritical in subcritical_values[1:]
        ),
        coordinate_gramian_reduces=(
            full_gramian.extract(outer_indices, inner_indices)
            == sp.zeros(outer_dimension, len(inner_indices))
            and full_gramian.extract(inner_indices, outer_indices)
            == sp.zeros(len(inner_indices), outer_dimension)
        ),
        exact_outer_coordinate_gramian=(
            sp.expand(
                full_gramian.extract(outer_indices, outer_indices)
                - expected_outer_gramian
            )
            == sp.zeros(outer_dimension)
        ),
    )
    if not all(
        (
            record.exact_outer_pencil,
            record.polynomial_cross_blocks_vanish,
            record.subcritical_compression_vanishes,
            record.coordinate_gramian_reduces,
            record.exact_outer_coordinate_gramian,
        )
    ):
        raise AssertionError(f"complex descent failed: {record}")
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-degree", type=int, default=5)
    parser.add_argument("--maximum-quotient", type=int, default=6)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact grid and optionally persist JSONL."""

    args = parse_args()
    if args.maximum_degree < 1 or args.maximum_quotient < 3:
        raise ValueError("invalid degree/quotient bounds")

    records = [
        make_record(degree, quotient)
        for degree in range(1, args.maximum_degree + 1)
        for quotient in range(3, args.maximum_quotient + 1)
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
