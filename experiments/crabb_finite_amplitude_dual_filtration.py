#!/usr/bin/env python3
"""Probe the reflected filtration at finite equality amplitude.

L143--L144 expand only to second order in equality amplitudes.  This
checker instead fixes exact nonzero rational equality coefficients and
expands the complete corrected finite-Blaschke lower bound in the
ellipse parameter ``c``.  If the smallest active equality grade is
``k``, the candidate reflected-ideal theorem predicts

    ||B_c(T_c)||^2 - 4 = O(c^(2k)).

The calculation includes every amplitude degree because the
phase-palindromic companion is substituted before the one-variable
``c`` expansion.  It is an exact finite falsification test, not the
all-size proof of the filtration.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from crabb_palindromic_elliptic_hessian import (
    direct_map_coefficients,
    inverse_map_coefficients,
)


DEFAULT_CASES = (
    (2, (Fraction(1, 20),)),
    (3, (Fraction(1, 20),)),
    (4, (Fraction(0), Fraction(1, 20))),
    (5, (Fraction(0), Fraction(1, 20))),
    (
        6,
        (Fraction(0), Fraction(1, 20), Fraction(1, 30)),
    ),
    (
        7,
        (Fraction(0), Fraction(1, 20), Fraction(-1, 30)),
    ),
    (
        8,
        (
            Fraction(0),
            Fraction(0),
            Fraction(1, 20),
            Fraction(1, 30),
        ),
    ),
)


@dataclass(frozen=True)
class FiniteAmplitudeDualRecord:
    """One exact finite-amplitude filtration record."""

    length: int
    dimension: int
    low_grade_coefficients: tuple[str, ...]
    minimum_grade: int
    target_degree: int
    first_nonzero_degree: int
    leading_coefficient: str
    leading_coefficient_negative: bool
    lower_coefficients_vanish: bool
    next_coefficient: str


def truncate(
    expression: sp.Expr,
    parameter: sp.Symbol,
    order: int,
) -> sp.Expr:
    """Truncate a polynomial expression in ``parameter``."""

    expanded = sp.expand(expression)
    return sp.expand(
        sum(
            expanded.coeff(parameter, degree) * parameter**degree
            for degree in range(order)
        )
    )


def truncate_matrix(
    matrix: sp.Matrix,
    parameter: sp.Symbol,
    order: int,
) -> sp.Matrix:
    """Truncate every entry of a matrix polynomial."""

    return matrix.applyfunc(
        lambda entry: truncate(entry, parameter, order)
    )


def coefficient_matrix(
    matrix: sp.Matrix,
    parameter: sp.Symbol,
    degree: int,
) -> sp.Matrix:
    """Extract one coefficient from every matrix entry."""

    return matrix.applyfunc(
        lambda entry: sp.expand(entry).coeff(parameter, degree)
    )


def series_to_expression(
    series: object,
    parameter: sp.Symbol,
    order: int,
) -> sp.Expr:
    """Convert the shared rational ``Series`` type to a SymPy polynomial."""

    return sum(
        (
            sp.Rational(coefficient.numerator, coefficient.denominator)
            * parameter**degree
            for degree, coefficient in enumerate(
                series.coefficients[:order]
            )
        ),
        sp.Integer(0),
    )


def full_coefficients(
    length: int,
    low_coefficients: tuple[Fraction, ...],
) -> dict[int, sp.Rational]:
    """Lift low real equality coefficients by phase-one reflection."""

    result: dict[int, sp.Rational] = {}
    for grade, value in enumerate(low_coefficients, start=1):
        coefficient = sp.Rational(value.numerator, value.denominator)
        result[grade] = coefficient
        result[length - grade] = coefficient
    return result


def reversal_matrix(dimension: int) -> sp.Matrix:
    """Return coordinate reversal."""

    result = sp.zeros(dimension)
    for index in range(dimension):
        result[index, dimension - 1 - index] = 1
    return result


def companion_pencil(
    length: int,
    coefficients: dict[int, sp.Rational],
    parameter: sp.Symbol,
) -> sp.Matrix:
    """Return the exact equality/physical-adjoint coefficient pencil."""

    dimension = length + 1
    crabb = sp.zeros(dimension)
    crabb[0, 1] = 2
    for column in range(2, dimension):
        crabb[column - 1, column] = 1

    equality_tangent = sp.zeros(dimension)
    for column in range(2, dimension):
        value = 2 * coefficients.get(column - 1, 0)
        equality_tangent[0, column] = value
        equality_tangent[length, column] = -value

    reversal = reversal_matrix(dimension)
    equality_companion = crabb + equality_tangent
    return (
        equality_companion
        + parameter * reversal * equality_companion * reversal
    )


def ellipse_pullback(
    pencil: sp.Matrix,
    parameter: sp.Symbol,
    order: int,
) -> sp.Matrix:
    """Evaluate the complete direct ellipse map on the pencil."""

    result = sp.zeros(pencil.rows)
    for scalar_degree, coefficient in enumerate(
        direct_map_coefficients(order - 1, order)
    ):
        scalar = series_to_expression(
            coefficient,
            parameter,
            order,
        )
        result = truncate_matrix(
            result
            + scalar * pencil ** (2 * scalar_degree + 1),
            parameter,
            order,
        )
    return result


def coordinate_gramian(
    length: int,
    coefficients: dict[int, sp.Rational],
) -> sp.Matrix:
    """Return L123's exact equality coordinate Gramian."""

    gramian = sp.diag(
        sp.Rational(1, 2),
        *([sp.Integer(1)] * (length - 1)),
        sp.Rational(1, 2),
    )
    for row in range(length):
        for column in range(row + 1, length):
            value = coefficients.get(column - row, 0)
            gramian[row, column] += value
            gramian[column, row] += value
            gramian[row + 1, column + 1] += value
            gramian[column + 1, row + 1] += value
    return gramian


def prepared_numerator(
    length: int,
    coefficients: dict[int, sp.Rational],
    parameter: sp.Symbol,
    variable: sp.Symbol,
    order: int,
) -> sp.Expr:
    """Prepare the complete finite-amplitude numerator in ``c``."""

    work_order = order + 2 * (order - 1) + 4
    inverse_coefficients = inverse_map_coefficients(
        order - 1,
        work_order,
    )
    inverse_map = sum(
        (
            series_to_expression(
                coefficient,
                parameter,
                order,
            )
            * variable ** (2 * scalar_degree + 1)
            for scalar_degree, coefficient in enumerate(
                inverse_coefficients
            )
        ),
        sp.Integer(0),
    )
    inverse_map = truncate(inverse_map, parameter, order)

    dickson = [sp.Integer(2), inverse_map]
    for _ in range(2, length + 1):
        dickson.append(
            truncate(
                inverse_map * dickson[-1]
                - parameter * dickson[-2],
                parameter,
                order,
            )
        )
    composed = truncate(
        dickson[length]
        + sum(
            (
                2 * coefficients.get(degree, 0) * dickson[degree]
                for degree in range(1, length)
            ),
            sp.Integer(0),
        ),
        parameter,
        order,
    )

    equality_numerator = sp.expand(composed.coeff(parameter, 0))
    numerator_coefficients = [equality_numerator]
    outer_coefficients = [sp.Integer(1)]
    for degree in range(1, order):
        residual = sp.expand(
            composed.coeff(parameter, degree)
            - sum(
                (
                    numerator_coefficients[left_degree]
                    * outer_coefficients[degree - left_degree]
                    for left_degree in range(1, degree)
                ),
                sp.Integer(0),
            )
        )
        quotient, remainder = sp.div(
            residual,
            equality_numerator,
            variable,
        )
        numerator_coefficients.append(sp.expand(remainder))
        outer_coefficients.append(sp.expand(quotient))
    return sum(
        (
            coefficient * parameter**degree
            for degree, coefficient in enumerate(
                numerator_coefficients
            )
        ),
        sp.Integer(0),
    )


def reversed_polynomial(
    numerator: sp.Expr,
    length: int,
    parameter: sp.Symbol,
    variable: sp.Symbol,
    order: int,
) -> sp.Expr:
    """Reverse every numerator ``c`` coefficient to degree ``length``."""

    result = sp.Integer(0)
    for parameter_degree in range(order):
        coefficient = numerator.coeff(parameter, parameter_degree)
        for (scalar_degree,), value in sp.Poly(
            coefficient,
            variable,
        ).terms():
            result += (
                value
                * variable ** (length - scalar_degree)
                * parameter**parameter_degree
            )
    return sp.expand(result)


def matrix_polynomial_evaluation(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    matrix: sp.Matrix,
    parameter: sp.Symbol,
    order: int,
) -> sp.Matrix:
    """Evaluate a scalar bivariate polynomial on a matrix series."""

    poly = sp.Poly(polynomial, variable)
    powers = [sp.eye(matrix.rows)]
    for _ in range(poly.degree()):
        powers.append(
            truncate_matrix(
                powers[-1] * matrix,
                parameter,
                order,
            )
        )
    result = sp.zeros(matrix.rows)
    for (degree,), coefficient in poly.terms():
        result = truncate_matrix(
            result + coefficient * powers[degree],
            parameter,
            order,
        )
    return result


def matrix_series_inverse(
    matrix: sp.Matrix,
    parameter: sp.Symbol,
    order: int,
) -> sp.Matrix:
    """Invert a matrix power series whose constant term is nonsingular."""

    coefficients = [
        coefficient_matrix(matrix, parameter, degree)
        for degree in range(order)
    ]
    inverse_coefficients = [coefficients[0].inv()]
    for degree in range(1, order):
        convolution = sum(
            (
                coefficients[left_degree]
                * inverse_coefficients[degree - left_degree]
                for left_degree in range(1, degree + 1)
            ),
            sp.zeros(matrix.rows),
        )
        inverse_coefficients.append(
            -inverse_coefficients[0] * convolution
        )
    return sum(
        (
            coefficient * parameter**degree
            for degree, coefficient in enumerate(
                inverse_coefficients
            )
        ),
        sp.zeros(matrix.rows),
    )


def simple_top_eigenvalue_series(
    matrix: sp.Matrix,
    coordinate_metric: sp.Matrix,
    reverse_defect: sp.Matrix,
    parameter: sp.Symbol,
    order: int,
) -> tuple[sp.Expr, ...]:
    """Lift the simple eigenvalue four of a rank-one base matrix."""

    generalized = truncate_matrix(
        coordinate_metric.inv() * matrix,
        parameter,
        order,
    )
    coefficients = [
        coefficient_matrix(generalized, parameter, degree)
        for degree in range(order)
    ]
    dimension = matrix.rows
    right_vectors = [sp.eye(dimension)[:, dimension - 1]]
    left_vector = (2 * reverse_defect).T
    eigenvalues = [sp.Integer(4)]
    for degree in range(1, order):
        forcing = sum(
            (
                coefficients[source_degree]
                * right_vectors[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.zeros(dimension, 1),
        )
        eigenvalue = sp.factor((left_vector * forcing)[0])
        eigenvalues.append(eigenvalue)
        residual = -forcing + sum(
            (
                eigenvalues[source_degree]
                * right_vectors[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.zeros(dimension, 1),
        )
        correction = -residual / 4
        if sp.simplify((left_vector * correction)[0]) != 0:
            raise AssertionError("the eigenvector normalization drifted")
        right_vectors.append(correction)
    return tuple(eigenvalues)


def make_record(
    length: int,
    low_coefficients: tuple[Fraction, ...],
) -> FiniteAmplitudeDualRecord:
    """Run one exact complete-amplitude dual calculation."""

    minimum_grade = next(
        grade
        for grade, coefficient in enumerate(low_coefficients, start=1)
        if coefficient
    )
    target_degree = 2 * minimum_grade
    order = target_degree + 2
    parameter, variable = sp.symbols("c w", real=True)
    coefficients = full_coefficients(length, low_coefficients)
    pencil = companion_pencil(length, coefficients, parameter)
    operator = ellipse_pullback(pencil, parameter, order)
    metric = coordinate_gramian(length, coefficients)
    numerator = prepared_numerator(
        length,
        coefficients,
        parameter,
        variable,
        order,
    )
    denominator = reversed_polynomial(
        numerator,
        length,
        parameter,
        variable,
        order,
    )
    numerator_value = matrix_polynomial_evaluation(
        numerator,
        variable,
        operator,
        parameter,
        order,
    )
    denominator_value = matrix_polynomial_evaluation(
        denominator,
        variable,
        operator,
        parameter,
        order,
    )
    blaschke = truncate_matrix(
        numerator_value
        * matrix_series_inverse(
            denominator_value,
            parameter,
            order,
        ),
        parameter,
        order,
    )
    gram = truncate_matrix(
        blaschke.T * metric * blaschke,
        parameter,
        order,
    )

    defect = sp.zeros(length + 1, 1)
    defect[0] = sp.Rational(1, 2)
    for degree in range(1, length):
        defect[degree] = coefficients.get(degree, 0)
    reverse_defect = reversal_matrix(length + 1) * defect
    eigenvalues = simple_top_eigenvalue_series(
        gram,
        metric,
        reverse_defect,
        parameter,
        order,
    )
    first_nonzero = next(
        degree
        for degree, coefficient in enumerate(eigenvalues[1:], start=1)
        if coefficient
    )
    leading = eigenvalues[target_degree]
    lower_vanish = all(
        eigenvalues[degree] == 0
        for degree in range(1, target_degree)
    )
    if (
        first_nonzero != target_degree
        or not lower_vanish
        or leading >= 0
    ):
        raise AssertionError(
            "the finite-amplitude dual filtration failed"
        )
    return FiniteAmplitudeDualRecord(
        length=length,
        dimension=length + 1,
        low_grade_coefficients=tuple(
            str(coefficient) for coefficient in low_coefficients
        ),
        minimum_grade=minimum_grade,
        target_degree=target_degree,
        first_nonzero_degree=first_nonzero,
        leading_coefficient=str(leading),
        leading_coefficient_negative=bool(leading < 0),
        lower_coefficients_vanish=lower_vanish,
        next_coefficient=str(eigenvalues[target_degree + 1]),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the fixed exact case grid."""

    args = parse_args()
    records = [
        make_record(length, coefficients)
        for length, coefficients in DEFAULT_CASES
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
