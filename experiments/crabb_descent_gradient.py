#!/usr/bin/env python3
"""Regenerate the exact algebra behind the L120 descent gradient.

The checker uses symbolic Dickson recurrences.  For each size it verifies:

* the Chebyshev descent to the reversal matrix;
* the polynomial conditional expectation on endpoint compressions;
* the differential fibre-quadrature identity for an arbitrary matrix
  direction; and
* the finite-path power of the elliptic parameter in the endpoint
  differential.

It also compares the resulting leading chain-rule coefficients with three
independently derived support/Riemann/Stein series.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from general_crabb_weighted_series import (
    crabb_matrix,
    weighted_condition_derivative,
)


@dataclass(frozen=True)
class DescentRecord:
    dimension: int
    checked_monomial_degree: int
    first_endpoint_power: int
    required_endpoint_power: int
    polynomial_reversal: bool
    conditional_expectation: bool
    differential_quadrature: bool


def dickson_values(
    matrix: sp.Matrix,
    length: int,
    parameter: sp.Symbol,
) -> list[sp.Matrix]:
    """Return ``D_n(matrix,parameter)`` through degree ``length``."""

    values = [2 * sp.eye(matrix.rows), matrix]
    for _ in range(2, length + 1):
        values.append(
            (
                matrix * values[-1] - parameter * values[-2]
            ).applyfunc(sp.expand)
        )
    return values[: length + 1]


def dickson_tangent(
    matrix: sp.Matrix,
    perturbation: sp.Matrix,
    values: list[sp.Matrix],
    parameter: sp.Symbol,
) -> sp.Matrix:
    """Differentiate the last matrix Dickson value."""

    length = len(values) - 1
    if length == 1:
        return perturbation
    previous = sp.zeros(matrix.rows)
    current = perturbation
    for degree in range(2, length + 1):
        previous, current = (
            current,
            (
                perturbation * values[degree - 1]
                + matrix * current
                - parameter * previous
            ).applyfunc(sp.expand),
        )
    return current


def dickson_at_matrix(
    matrix: sp.Matrix,
    degree: int,
    parameter: sp.Expr,
) -> sp.Matrix:
    """Evaluate one scalar-parameter Dickson polynomial at a matrix."""

    if degree == 0:
        return 2 * sp.eye(matrix.rows)
    previous = 2 * sp.eye(matrix.rows)
    current = matrix
    for _ in range(2, degree + 1):
        previous, current = (
            current,
            (matrix * current - parameter * previous).applyfunc(
                sp.expand
            ),
        )
    return current


def fibre_average_of_monomial(
    degree: int,
    length: int,
    parameter: sp.Symbol,
    quotient: sp.Matrix,
) -> sp.Matrix:
    """Evaluate the fibre average of ``z**degree`` at ``quotient``."""

    retained: dict[int, sp.Expr] = {}
    for lower_steps in range(degree + 1):
        exponent = degree - 2 * lower_steps
        if exponent % length == 0:
            quotient_power = exponent // length
            retained[quotient_power] = (
                retained.get(quotient_power, 0)
                + sp.binomial(degree, lower_steps)
                * parameter**lower_steps
            )

    result = retained.get(0, 0) * sp.eye(quotient.rows)
    descended_parameter = parameter**length
    for quotient_power in sorted(
        power for power in retained if power > 0
    ):
        positive = retained[quotient_power]
        negative = retained.get(-quotient_power, 0)
        expected_negative = positive * descended_parameter**quotient_power
        if sp.expand(negative - expected_negative) != 0:
            raise RuntimeError("the fibre average lost reflection symmetry")
        result += (
            positive
            * dickson_at_matrix(
                quotient, quotient_power, descended_parameter
            )
        )
    return result.applyfunc(sp.expand)


def audit_dimension(dimension: int) -> DescentRecord:
    """Run every exact identity in one dimension."""

    length = dimension - 1
    parameter = sp.symbols("c", positive=True)
    crabb = crabb_matrix(dimension)
    axis = crabb + parameter * crabb.T
    values = dickson_values(axis, length, parameter)
    descended = values[-1]

    diagonal = sp.diag(
        *[
            parameter ** sp.Rational(index, 2)
            for index in range(dimension)
        ]
    )
    reversal = sp.zeros(dimension)
    for index in range(dimension):
        reversal[index, length - index] = 1
    expected = (
        2
        * parameter ** sp.Rational(length, 2)
        * diagonal
        * reversal
        * diagonal.inv()
    )
    polynomial_reversal = descended.equals(expected)

    endpoint_indices = [0, length]
    quotient = descended.extract(endpoint_indices, endpoint_indices)
    checked_monomial_degree = 3 * length
    conditional_expectation = True
    for degree in range(checked_monomial_degree + 1):
        compression = (axis**degree).extract(
            endpoint_indices, endpoint_indices
        )
        expected_compression = fibre_average_of_monomial(
            degree, length, parameter, quotient
        )
        if not compression.equals(expected_compression):
            conditional_expectation = False
            break

    coordinates = sp.symbols(
        f"y0:{dimension * dimension}", complex=True
    )
    perturbation = sp.Matrix(dimension, dimension, coordinates)
    tangent = dickson_tangent(
        axis, perturbation, values, parameter
    )

    left_coefficients = {
        -1: sp.expand(tangent[length, 0] / 2),
        0: sp.expand(
            (tangent[0, 0] + tangent[length, length]) / 2
        ),
        1: sp.expand(tangent[0, length] / 2),
    }
    right_coefficients = {-1: 0, 0: 0, 1: 0}
    endpoint_factor = 1 / sp.sqrt(2)
    support_weights = [
        endpoint_factor,
        *[sp.Integer(1) for _ in range(length - 1)],
        endpoint_factor,
    ]
    for row in range(dimension):
        for column in range(dimension):
            coordinate = coordinates[row * dimension + column]
            for lower_steps in range(length):
                exponent = (
                    column
                    - row
                    + length
                    - 1
                    - 2 * lower_steps
                )
                if exponent % length != 0:
                    continue
                quotient_power = exponent // length
                if quotient_power not in right_coefficients:
                    continue
                right_coefficients[quotient_power] += (
                    support_weights[row]
                    * support_weights[column]
                    * coordinate
                    * parameter**lower_steps
                )
    differential_quadrature = all(
        sp.expand(
            left_coefficients[power]
            - right_coefficients[power]
        )
        == 0
        for power in left_coefficients
    )

    endpoint_functional = sp.expand(
        tangent[length, 0]
        - parameter**length * tangent[0, length]
    )
    endpoint_polynomial = sp.Poly(
        endpoint_functional, parameter
    )
    first_endpoint_power = min(
        monomial[0]
        for monomial, coefficient in endpoint_polynomial.terms()
        if coefficient != 0
    )

    if not (
        polynomial_reversal
        and conditional_expectation
        and differential_quadrature
        and first_endpoint_power >= length // 2
    ):
        raise RuntimeError(f"descent audit failed in size {dimension}")

    return DescentRecord(
        dimension=dimension,
        checked_monomial_degree=checked_monomial_degree,
        first_endpoint_power=first_endpoint_power,
        required_endpoint_power=length // 2,
        polynomial_reversal=polynomial_reversal,
        conditional_expectation=conditional_expectation,
        differential_quadrature=differential_quadrature,
    )


def audit_sparse_coefficients() -> dict[str, str]:
    """Compare L120's chain coefficients with independent exact series."""

    checks: dict[str, tuple[sp.Expr, sp.Expr]] = {}

    dimension = 3
    direction = crabb_matrix(dimension).T
    direct = weighted_condition_derivative(
        dimension, direction, extra_order=2
    )
    parameter = next(
        symbol
        for symbol in direct.free_symbols
        if symbol.name == "ellipse_parameter"
    )
    checks["p3_axis"] = (
        direct,
        -64 * parameter**5,
    )

    dimension = 4
    direction = sp.zeros(dimension)
    direction[3, 0] = 1
    direct = weighted_condition_derivative(
        dimension, direction, extra_order=1
    )
    parameter = next(iter(direct.free_symbols))
    checks["p4_bottom"] = (
        direct,
        -16 * parameter**7,
    )

    dimension = 5
    direction = sp.zeros(dimension)
    direction[3, 0] = 1 / sp.sqrt(2)
    direction[4, 1] = 1 / sp.sqrt(2)
    direct = weighted_condition_derivative(
        dimension, direction, extra_order=2
    )
    parameter = next(iter(direct.free_symbols))
    checks["p5_penultimate"] = (
        direct,
        -32 * parameter**10,
    )

    for name, (direct, predicted) in checks.items():
        if sp.expand(direct - predicted) != 0:
            raise RuntimeError(f"chain coefficient mismatch: {name}")
    return {
        name: str(predicted)
        for name, (_, predicted) in checks.items()
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSONL path for audit records",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = [
        audit_dimension(dimension)
        for dimension in range(3, 11)
    ]
    coefficient_record = {
        "independent_sparse_coefficients": audit_sparse_coefficients()
    }
    lines = [
        *[
            json.dumps(asdict(record), sort_keys=True)
            for record in records
        ],
        json.dumps(coefficient_record, sort_keys=True),
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
