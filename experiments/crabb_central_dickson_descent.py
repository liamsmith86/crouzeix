#!/usr/bin/env python3
"""Regenerate the exact central Dickson descent from L126.

For ``p=2k+1``, the degree-``k`` Dickson polynomial of the coefficient
Crabb ellipse has a reducing outer three-dimensional block.  The same
is true of its Fréchet derivative in the central palindromic equality
direction.  This script verifies both polynomial identities exactly.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class DicksonDescentRecord:
    degree: int
    dimension: int
    outer_axis_identity: bool
    outer_derivative_identity: bool
    full_axis_formula: bool
    full_derivative_formula: bool
    axis_cross_blocks_vanish: bool
    derivative_cross_blocks_vanish: bool


def coefficient_crabb(dimension: int) -> sp.Matrix:
    """Return L123's coefficient-coordinate Crabb shift."""

    matrix = sp.zeros(dimension)
    matrix[0, 1] = 2
    for column in range(2, dimension):
        matrix[column - 1, column] = 1
    return matrix


def dickson_value_and_derivative(
    matrix: sp.Matrix,
    direction: sp.Matrix,
    parameter: sp.Symbol,
    degree: int,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Evaluate ``P_degree`` and its Fréchet derivative by recurrence."""

    dimension = matrix.rows
    previous_value = 2 * sp.eye(dimension)
    value = matrix
    previous_derivative = sp.zeros(dimension)
    derivative = direction
    if degree == 1:
        return value, derivative

    for _ in range(2, degree + 1):
        next_value = sp.expand(
            matrix * value - parameter * previous_value
        )
        next_derivative = sp.expand(
            direction * value
            + matrix * derivative
            - parameter * previous_derivative
        )
        previous_value, value = value, next_value
        previous_derivative, derivative = derivative, next_derivative
    return value, derivative


def make_record(degree: int) -> DicksonDescentRecord:
    """Construct and exactly verify one central descent."""

    if degree < 1:
        raise ValueError("the Dickson degree must be positive")
    parameter = sp.symbols("c")
    dimension = 2 * degree + 1
    crabb = coefficient_crabb(dimension)
    reversal = sp.eye(dimension)[:, ::-1]
    axis = crabb + parameter * reversal * crabb * reversal

    endpoint_difference = sp.zeros(dimension, 1)
    endpoint_difference[0] = 1
    endpoint_difference[-1] = -1
    central_row = (
        sp.eye(dimension)[:, degree + 1].T
        - parameter * sp.eye(dimension)[:, degree - 1].T
    )
    direction = 2 * endpoint_difference * central_row

    value, derivative = dickson_value_and_derivative(
        axis,
        direction,
        parameter,
        degree,
    )
    outer_indices = [0, degree, 2 * degree]
    inner_indices = [
        index
        for index in range(dimension)
        if index not in outer_indices
    ]
    descended_parameter = parameter**degree
    expected_axis = sp.Matrix(
        [
            [0, 2, 0],
            [descended_parameter, 0, 1],
            [0, 2 * descended_parameter, 0],
        ]
    )
    expected_derivative = sp.Matrix(
        [
            [-2 * descended_parameter, 0, 2],
            [0, 0, 0],
            [2 * descended_parameter, 0, -2],
        ]
    )
    expected_full_axis = sp.zeros(dimension)
    expected_full_derivative = sp.zeros(dimension)
    expected_full_axis[degree, 0] = descended_parameter
    expected_full_derivative[0, 0] = -2 * descended_parameter
    expected_full_derivative[2 * degree, 0] = (
        2 * descended_parameter
    )
    for column in range(1, degree):
        expected_full_axis[degree - column, column] = (
            parameter ** (degree - column)
        )
        expected_full_axis[degree + column, column] = (
            descended_parameter
        )
        expected_full_derivative[column, column] = (
            -2 * descended_parameter
        )
        expected_full_derivative[
            2 * degree - column,
            column,
        ] = 2 * parameter ** (degree - column)
    expected_full_axis[0, degree] = 2
    expected_full_axis[2 * degree, degree] = (
        2 * descended_parameter
    )
    for column in range(degree + 1, 2 * degree):
        expected_full_axis[column - degree, column] = 1
        expected_full_axis[3 * degree - column, column] = (
            parameter ** (2 * degree - column)
        )
        expected_full_derivative[column, column] = -2
        expected_full_derivative[
            2 * degree - column,
            column,
        ] = 2 * parameter ** (2 * degree - column)
    expected_full_axis[degree, 2 * degree] = 1
    expected_full_derivative[0, 2 * degree] = 2
    expected_full_derivative[2 * degree, 2 * degree] = -2

    axis_outer = value.extract(outer_indices, outer_indices)
    derivative_outer = derivative.extract(outer_indices, outer_indices)
    axis_cross_vanishes = (
        value.extract(outer_indices, inner_indices)
        == sp.zeros(3, dimension - 3)
        and value.extract(inner_indices, outer_indices)
        == sp.zeros(dimension - 3, 3)
    )
    derivative_cross_vanishes = (
        derivative.extract(outer_indices, inner_indices)
        == sp.zeros(3, dimension - 3)
        and derivative.extract(inner_indices, outer_indices)
        == sp.zeros(dimension - 3, 3)
    )
    record = DicksonDescentRecord(
        degree=degree,
        dimension=dimension,
        outer_axis_identity=axis_outer == expected_axis,
        outer_derivative_identity=derivative_outer == expected_derivative,
        full_axis_formula=value == expected_full_axis,
        full_derivative_formula=derivative == expected_full_derivative,
        axis_cross_blocks_vanish=axis_cross_vanishes,
        derivative_cross_blocks_vanish=derivative_cross_vanishes,
    )
    if not all(
        (
            record.outer_axis_identity,
            record.outer_derivative_identity,
            record.full_axis_formula,
            record.full_derivative_formula,
            record.axis_cross_blocks_vanish,
            record.derivative_cross_blocks_vanish,
        )
    ):
        raise AssertionError("the exact central Dickson descent failed")
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-degree", type=int, default=12)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic symbolic regeneration."""

    args = parse_args()
    records = [
        make_record(degree)
        for degree in range(1, args.maximum_degree + 1)
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
