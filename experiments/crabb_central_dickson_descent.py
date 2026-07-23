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
    exact_amplitude_linearity: bool
    coordinate_gramian_reduces: bool
    outer_coordinate_gramian_identity: bool
    inner_two_layer_identity: bool
    inner_coordinate_gramian_identity: bool
    inner_involution_identity: bool
    subcritical_dickson_compression: bool
    subcritical_moment_compression: bool
    subcritical_boundary_annihilation: bool
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
    parameter, amplitude = sp.symbols("c a")
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

    axis_dickson_values = [2 * sp.eye(dimension), axis]
    pencil_dickson_values = [
        2 * sp.eye(dimension),
        axis + amplitude * direction,
    ]
    for current_degree in range(2, degree):
        axis_dickson_values.append(
            sp.expand(
                axis * axis_dickson_values[-1]
                - parameter * axis_dickson_values[-2]
            )
        )
        pencil_dickson_values.append(
            sp.expand(
                (axis + amplitude * direction)
                * pencil_dickson_values[-1]
                - parameter * pencil_dickson_values[-2]
            )
        )

    value, derivative = dickson_value_and_derivative(
        axis,
        direction,
        parameter,
        degree,
    )
    pencil_value, _ = dickson_value_and_derivative(
        axis + amplitude * direction,
        sp.zeros(dimension),
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

    coordinate_gramian = sp.diag(
        sp.Rational(1, 2),
        *(sp.Integer(1) for _ in range(dimension - 2)),
        sp.Rational(1, 2),
    )
    for row in range(degree):
        column = row + degree
        coordinate_gramian[row, column] += amplitude
        coordinate_gramian[column, row] += amplitude
        coordinate_gramian[row + 1, column + 1] += amplitude
        coordinate_gramian[column + 1, row + 1] += amplitude
    expected_outer_coordinate_gramian = sp.Matrix(
        [
            [sp.Rational(1, 2), amplitude, 0],
            [amplitude, 1, amplitude],
            [0, amplitude, sp.Rational(1, 2)],
        ]
    )
    inner_two_layer_identity = True
    inner_coordinate_gramian_identity = True
    inner_involution_identity = True
    if degree > 1:
        layer_dimension = degree - 1
        layer_reversal = sp.eye(layer_dimension)[:, ::-1]
        layer_diagonal = sp.diag(
            *(
                parameter ** (degree - index)
                for index in range(1, degree)
            )
        )
        layer_operator = layer_reversal * layer_diagonal
        inner_indices_ordered = [
            *range(1, degree),
            *range(degree + 1, 2 * degree),
        ]
        identity = sp.eye(layer_dimension)
        expected_inner_pencil = sp.Matrix.vstack(
            sp.Matrix.hstack(
                layer_operator - 2 * amplitude * descended_parameter * identity,
                identity + 2 * amplitude * layer_operator,
            ),
            sp.Matrix.hstack(
                descended_parameter * identity
                + 2 * amplitude * layer_operator,
                layer_operator - 2 * amplitude * identity,
            ),
        )
        expected_inner_gramian = sp.Matrix.vstack(
            sp.Matrix.hstack(identity, 2 * amplitude * identity),
            sp.Matrix.hstack(2 * amplitude * identity, identity),
        )
        inner_two_layer_identity = (
            pencil_value.extract(
                inner_indices_ordered,
                inner_indices_ordered,
            )
            == expected_inner_pencil
        )
        inner_coordinate_gramian_identity = (
            coordinate_gramian.extract(
                inner_indices_ordered,
                inner_indices_ordered,
            )
            == expected_inner_gramian
        )
        inner_involution_identity = (
            layer_operator**2
            == descended_parameter * identity
        )

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
    outer_embedding = sp.eye(dimension)[:, outer_indices]
    subcritical_dickson_compression = all(
        value.extract(outer_indices, outer_indices) == sp.zeros(3)
        for value in pencil_dickson_values[1:degree]
    )
    subcritical_moment_compression = all(
        ((axis + amplitude * direction) ** moment).extract(
            outer_indices,
            outer_indices,
        )
        == (
            sp.binomial(moment, moment // 2)
            * parameter ** (moment // 2)
            * sp.eye(3)
            if moment % 2 == 0
            else sp.zeros(3)
        )
        for moment in range(degree)
    )
    endpoint_column = 2 * endpoint_difference
    subcritical_boundary_annihilation = all(
        (
            central_row
            * axis_dickson_values[current_degree]
            * outer_embedding
            == sp.zeros(1, 3)
        )
        for current_degree in range(max(0, degree - 1))
    ) and all(
        (
            outer_embedding.T
            * axis_dickson_values[current_degree]
            * endpoint_column
            == sp.zeros(3, 1)
        )
        for current_degree in range(1, degree)
    )
    record = DicksonDescentRecord(
        degree=degree,
        dimension=dimension,
        outer_axis_identity=axis_outer == expected_axis,
        outer_derivative_identity=derivative_outer == expected_derivative,
        full_axis_formula=value == expected_full_axis,
        full_derivative_formula=derivative == expected_full_derivative,
        exact_amplitude_linearity=(
            pencil_value == value + amplitude * derivative
        ),
        coordinate_gramian_reduces=(
            coordinate_gramian.extract(outer_indices, inner_indices)
            == sp.zeros(3, dimension - 3)
            and coordinate_gramian.extract(inner_indices, outer_indices)
            == sp.zeros(dimension - 3, 3)
        ),
        outer_coordinate_gramian_identity=(
            coordinate_gramian.extract(outer_indices, outer_indices)
            == expected_outer_coordinate_gramian
        ),
        inner_two_layer_identity=inner_two_layer_identity,
        inner_coordinate_gramian_identity=(
            inner_coordinate_gramian_identity
        ),
        inner_involution_identity=inner_involution_identity,
        subcritical_dickson_compression=(
            subcritical_dickson_compression
        ),
        subcritical_moment_compression=(
            subcritical_moment_compression
        ),
        subcritical_boundary_annihilation=(
            subcritical_boundary_annihilation
        ),
        axis_cross_blocks_vanish=axis_cross_vanishes,
        derivative_cross_blocks_vanish=derivative_cross_vanishes,
    )
    if not all(
        (
            record.outer_axis_identity,
            record.outer_derivative_identity,
            record.full_axis_formula,
            record.full_derivative_formula,
            record.exact_amplitude_linearity,
            record.coordinate_gramian_reduces,
            record.outer_coordinate_gramian_identity,
            record.inner_two_layer_identity,
            record.inner_coordinate_gramian_identity,
            record.inner_involution_identity,
            record.subcritical_dickson_compression,
            record.subcritical_moment_compression,
            record.subcritical_boundary_annihilation,
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
