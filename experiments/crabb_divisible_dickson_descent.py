#!/usr/bin/env python3
"""Audit exact Dickson descent whenever the circle grade divides the chain.

Write ``L=qk``.  For the phase-one equality direction supported at
grades ``k`` and ``L-k``, the degree-``k`` Dickson polynomial of the
size-``L+1`` coefficient pencil has an exact reducing residue-zero
block.  That block is the complete size-``q+1`` offset-one
equality/ellipse pencil at parameter ``c^k``.  The coordinate Gramian
reduces the same block, and every subcritical Dickson grade compresses
to zero there.

This is the algebraic input for the divisible-grade metric descent.  It
strictly generalizes the central case ``q=2`` checked by
``crabb_central_dickson_descent.py``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_central_dickson_descent import (
    coefficient_crabb,
    dickson_value_and_derivative,
)


@dataclass(frozen=True)
class DivisibleDescentRecord:
    degree: int
    quotient: int
    length: int
    dimension: int
    outer_dimension: int
    exact_outer_pencil: bool
    polynomial_cross_blocks_vanish: bool
    exact_amplitude_linearity: bool
    coordinate_gramian_reduces: bool
    exact_outer_coordinate_gramian: bool
    subcritical_compression_vanishes: bool
    inactive_apex_shifts: bool


def reversal(dimension: int) -> sp.Matrix:
    """Return coordinate reversal."""

    return sp.eye(dimension)[:, ::-1]


def palindromic_coefficients(
    length: int,
    first_grade: int,
) -> list[int]:
    """Return the phase-one pair at ``k`` and ``L-k``."""

    coefficients = [0] * (length - 1)
    coefficients[first_grade - 1] = 1
    coefficients[length - first_grade - 1] = 1
    return coefficients


def equality_direction(
    length: int,
    first_grade: int,
) -> sp.Matrix:
    """Return the coefficient-coordinate disk equality tangent."""

    dimension = length + 1
    coefficients = palindromic_coefficients(length, first_grade)
    direction = sp.zeros(dimension)
    for grade, coefficient in enumerate(coefficients, start=1):
        direction[0, grade + 1] = 2 * coefficient
        direction[length, grade + 1] = -2 * coefficient
    return direction


def equality_pencil(
    length: int,
    first_grade: int,
    parameter: sp.Expr,
    amplitude: sp.Symbol,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return the axis, physical-adjoint tangent, and full pencil."""

    dimension = length + 1
    flip = reversal(dimension)
    crabb = coefficient_crabb(dimension)
    axis = crabb + parameter * flip * crabb * flip
    disk_direction = equality_direction(length, first_grade)
    direction = (
        disk_direction
        + parameter * flip * disk_direction * flip
    )
    return axis, direction, axis + amplitude * direction


def dickson_values(
    matrix: sp.Matrix,
    parameter: sp.Expr,
    maximum_degree: int,
) -> tuple[sp.Matrix, ...]:
    """Return Dickson values through ``maximum_degree``."""

    values = [2 * sp.eye(matrix.rows)]
    if maximum_degree == 0:
        return tuple(values)
    values.append(matrix)
    for _ in range(2, maximum_degree + 1):
        values.append(
            sp.expand(
                matrix * values[-1] - parameter * values[-2]
            )
        )
    return tuple(values)


def coordinate_gramian(
    length: int,
    first_grade: int,
    amplitude: sp.Symbol,
) -> sp.Matrix:
    """Return L123's Toeplitz coordinate Gramian."""

    dimension = length + 1
    coefficients = palindromic_coefficients(length, first_grade)
    toeplitz = sp.zeros(dimension)
    for index in range(length):
        toeplitz[index, index] = sp.Rational(1, 2)
    for grade, coefficient in enumerate(coefficients, start=1):
        for row in range(length - grade):
            toeplitz[row, row + grade] = amplitude * coefficient
            toeplitz[row + grade, row] = amplitude * coefficient

    shift = sp.zeros(dimension)
    for column in range(1, dimension):
        shift[column - 1, column] = 1
    return sp.expand(toeplitz + shift.T * toeplitz * shift)


def make_record(degree: int, quotient: int) -> DivisibleDescentRecord:
    """Construct and audit one exact divisible-grade descent."""

    if degree < 1:
        raise ValueError("degree must be positive")
    if quotient < 2:
        raise ValueError("quotient must be at least two")

    parameter, amplitude = sp.symbols("c a")
    length = degree * quotient
    dimension = length + 1
    outer_dimension = quotient + 1
    axis, direction, pencil = equality_pencil(
        length,
        degree,
        parameter,
        amplitude,
    )
    value, derivative = dickson_value_and_derivative(
        axis,
        direction,
        parameter,
        degree,
    )
    pencil_value = dickson_values(
        pencil,
        parameter,
        degree,
    )[-1]

    descended_parameter = parameter**degree
    _, _, expected_outer = equality_pencil(
        quotient,
        1,
        descended_parameter,
        amplitude,
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
    outer_pencil = pencil_value.extract(
        outer_indices,
        outer_indices,
    )
    polynomial_cross_blocks_vanish = (
        pencil_value.extract(outer_indices, inner_indices)
        == sp.zeros(outer_dimension, len(inner_indices))
        and pencil_value.extract(inner_indices, outer_indices)
        == sp.zeros(len(inner_indices), outer_dimension)
    )

    full_gramian = coordinate_gramian(
        length,
        degree,
        amplitude,
    )
    expected_outer_gramian = coordinate_gramian(
        quotient,
        1,
        amplitude,
    )
    coordinate_cross_blocks_vanish = (
        full_gramian.extract(outer_indices, inner_indices)
        == sp.zeros(outer_dimension, len(inner_indices))
        and full_gramian.extract(inner_indices, outer_indices)
        == sp.zeros(len(inner_indices), outer_dimension)
    )

    subcritical_values = dickson_values(
        pencil,
        parameter,
        max(0, degree - 1),
    )
    subcritical_compression_vanishes = all(
        subcritical_value.extract(outer_indices, outer_indices)
        == sp.zeros(outer_dimension)
        for subcritical_value in subcritical_values[1:]
    )

    apex_value = value.subs(parameter, 0)
    inactive_apex_shifts = True
    expected_inactive = coefficient_crabb(quotient)
    expected_inactive[0, 1] = 1
    for residue in range(1, degree):
        residue_indices = [
            residue + step * degree
            for step in range(quotient)
        ]
        inactive_apex_shifts = (
            inactive_apex_shifts
            and apex_value.extract(
                residue_indices,
                residue_indices,
            )
            == expected_inactive
        )

    record = DivisibleDescentRecord(
        degree=degree,
        quotient=quotient,
        length=length,
        dimension=dimension,
        outer_dimension=outer_dimension,
        exact_outer_pencil=(
            sp.expand(outer_pencil - expected_outer)
            == sp.zeros(outer_dimension)
        ),
        polynomial_cross_blocks_vanish=(
            polynomial_cross_blocks_vanish
        ),
        exact_amplitude_linearity=(
            sp.expand(
                pencil_value - value - amplitude * derivative
            )
            == sp.zeros(dimension)
        ),
        coordinate_gramian_reduces=coordinate_cross_blocks_vanish,
        exact_outer_coordinate_gramian=(
            sp.expand(
                full_gramian.extract(outer_indices, outer_indices)
                - expected_outer_gramian
            )
            == sp.zeros(outer_dimension)
        ),
        subcritical_compression_vanishes=(
            subcritical_compression_vanishes
        ),
        inactive_apex_shifts=inactive_apex_shifts,
    )
    if not all(
        (
            record.exact_outer_pencil,
            record.polynomial_cross_blocks_vanish,
            record.exact_amplitude_linearity,
            record.coordinate_gramian_reduces,
            record.exact_outer_coordinate_gramian,
            record.subcritical_compression_vanishes,
            record.inactive_apex_shifts,
        )
    ):
        raise AssertionError(
            "the exact divisible-grade Dickson descent failed: "
            f"{record}"
        )
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-degree", type=int, default=6)
    parser.add_argument("--maximum-quotient", type=int, default=6)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact grid and optionally persist JSONL."""

    args = parse_args()
    if args.maximum_degree < 1 or args.maximum_quotient < 2:
        raise ValueError("invalid degree/quotient bounds")

    records = [
        make_record(degree, quotient)
        for degree in range(1, args.maximum_degree + 1)
        for quotient in range(2, args.maximum_quotient + 1)
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
