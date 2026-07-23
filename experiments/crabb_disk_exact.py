"""Shared exact algebra for the Toeplitz disk chart.

The research checkers in this directory use the same coefficient
realization repeatedly.  This module keeps the rational matrix,
polynomial, Fréchet, and Stein constructions in one place.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def coefficient_model(
    coefficients: Sequence[sp.Expr],
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return exact ``A,K,H,R`` for one Toeplitz disk-chart point."""

    length = len(coefficients) + 1
    dimension = length + 1
    toeplitz = sp.zeros(dimension)
    for index in range(length):
        toeplitz[index, index] = sp.Rational(1, 2)
    for offset, coefficient in enumerate(coefficients, start=1):
        for row in range(length - offset):
            toeplitz[row, row + offset] = coefficient
            toeplitz[row + offset, row] = sp.conjugate(coefficient)
    shift = sp.zeros(dimension)
    for row in range(length):
        shift[row, row + 1] = 1
    coordinate = toeplitz + shift.T.conjugate() * toeplitz * shift
    operator = sp.simplify(
        2 * coordinate.inv() * toeplitz * shift
    )
    return operator, coordinate, toeplitz, shift


def characteristic_factor(
    operator: sp.Matrix,
    variable: sp.Symbol,
) -> sp.Expr:
    """Return ``g`` from ``det(variable I-A)=variable*g``."""

    characteristic = operator.charpoly(variable).as_expr()
    return sp.cancel(characteristic / variable)


def ascending_coefficients(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    degree: int,
) -> list[sp.Expr]:
    """Return coefficients of powers zero through ``degree``."""

    poly = sp.Poly(polynomial, variable)
    return [poly.nth(power) for power in range(degree + 1)]


def evaluate_polynomial(
    coefficients: Sequence[sp.Expr],
    matrix: sp.Matrix,
) -> sp.Matrix:
    """Evaluate an ascending-order polynomial at a matrix."""

    result = sp.zeros(matrix.rows)
    power = sp.eye(matrix.rows)
    for coefficient in coefficients:
        result += coefficient * power
        power *= matrix
    return sp.simplify(result)


def frechet_polynomial(
    coefficients: Sequence[sp.Expr],
    matrix: sp.Matrix,
    tangent: sp.Matrix,
) -> sp.Matrix:
    """Evaluate the polynomial Fréchet derivative."""

    result = sp.zeros(matrix.rows)
    for exponent, coefficient in enumerate(coefficients):
        for left_power in range(exponent):
            result += (
                coefficient
                * matrix**left_power
                * tangent
                * matrix ** (exponent - 1 - left_power)
            )
    return sp.simplify(result)


def solve_symmetric_stein(
    operator: sp.Matrix,
    defect: sp.Matrix,
) -> sp.Matrix:
    """Solve ``P-A.T*P*A=defect*defect.T`` exactly over the reals."""

    dimension = operator.rows
    variables: list[sp.Symbol] = []
    metric = sp.zeros(dimension)
    for row in range(dimension):
        for column in range(row, dimension):
            entry = sp.symbols(f"p_{row}_{column}")
            variables.append(entry)
            metric[row, column] = entry
            metric[column, row] = entry
    stein = (
        metric
        - operator.T * metric * operator
        - defect * defect.T
    )
    solution = sp.solve(
        [
            stein[row, column]
            for row in range(dimension)
            for column in range(row, dimension)
        ],
        variables,
        dict=True,
    )
    if len(solution) != 1:
        raise AssertionError("the exact Stein equation was not unique")
    return sp.simplify(metric.subs(solution[0]))
