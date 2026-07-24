#!/usr/bin/env python3
"""Exact size-four counterexample to ambient ``O(Q)`` division.

At the Crabb apex, take a generic Toeplitz disk ray ``z=x h`` and the
fixed ambient matrix unit ``E_(1,2)``.  This script computes the first
inverse-Riemann correction and the characteristic-Blaschke Frechet
derivative over exact Gaussian rationals.  The ambient gradient starts
at order ``x^2``, while ``Q(xh)`` starts at order ``x^4``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json

import sympy as sp


DIMENSION = 4
LENGTH = DIMENSION - 1
SERIES_ORDER = 3


@dataclass(frozen=True)
class ExactAmbientGradientRecord:
    """The exact apex counterexample coefficients."""

    dimension: int
    length: int
    ambient_direction: str
    disk_direction: tuple[str, ...]
    quartic_defect_coefficient: str
    gradient_quadratic_coefficient: str
    gradient_series: str
    strict_order_separation: bool


def truncate(value: sp.Expr, parameter: sp.Symbol) -> sp.Expr:
    """Truncate a scalar modulo the third parameter power."""

    return sp.expand(
        sp.series(
            value,
            parameter,
            0,
            SERIES_ORDER,
        ).removeO()
    )


def truncate_matrix(
    matrix: sp.Matrix,
    parameter: sp.Symbol,
) -> sp.Matrix:
    """Truncate every matrix entry."""

    return matrix.applyfunc(lambda value: truncate(value, parameter))


def disk_model(
    parameter: sp.Symbol,
    direction: tuple[sp.Expr, ...],
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return the exact size-four disk operator and coefficient metric."""

    toeplitz = sp.zeros(DIMENSION)
    for index in range(LENGTH):
        toeplitz[index, index] = sp.Rational(1, 2)
    for offset in range(1, LENGTH):
        for row in range(LENGTH - offset):
            value = parameter * direction[offset]
            toeplitz[row, row + offset] = value
            toeplitz[row + offset, row] = sp.conjugate(value)

    shift = sp.zeros(DIMENSION)
    for row in range(LENGTH):
        shift[row, row + 1] = 1
    metric = toeplitz + shift.T * toeplitz * shift
    operator = 2 * metric.inv() * toeplitz * shift
    return (
        truncate_matrix(operator, parameter),
        truncate_matrix(metric, parameter),
    )


def circle_star(
    expression: sp.Expr,
    parameter: sp.Symbol,
    boundary: sp.Symbol,
) -> sp.Expr:
    """Conjugate a Laurent polynomial on the unit circle."""

    return sp.expand(
        sp.conjugate(expression).subs(
            {
                sp.conjugate(boundary): 1 / boundary,
                sp.conjugate(parameter): parameter,
            }
        )
    )


def laurent_coefficients(
    expression: sp.Expr,
    boundary: sp.Symbol,
) -> dict[int, sp.Expr]:
    """Return the finite Laurent coefficient dictionary."""

    coefficients: dict[int, sp.Expr] = {}
    for term in sp.Add.make_args(sp.expand(expression)):
        power = int(term.as_powers_dict().get(boundary, 0))
        coefficients[power] = sp.simplify(
            coefficients.get(power, 0) + term / boundary**power
        )
    return coefficients


def riemann_pulled_direction(
    operator: sp.Matrix,
    metric: sp.Matrix,
    direction: sp.Matrix,
    parameter: sp.Symbol,
) -> sp.Matrix:
    """Compute the exact first pulled ambient direction through order two."""

    boundary = sp.symbols("boundary")
    monomials = sp.Matrix(
        [boundary**index for index in range(DIMENSION)]
    )
    adjoint_monomials = sp.Matrix(
        [[boundary ** (-index) for index in range(DIMENSION)]]
    )
    denominator = sp.expand(
        (adjoint_monomials * metric * monomials)[0]
    )
    numerator = sp.expand(
        boundary**-1
        * (adjoint_monomials * metric * direction * monomials)[0]
    )
    support = truncate(
        (
            numerator
            + circle_star(numerator, parameter, boundary)
        )
        / (2 * denominator),
        parameter,
    )
    fourier = laurent_coefficients(support, boundary)

    correction = sp.zeros(DIMENSION)
    power = sp.eye(DIMENSION)
    maximum_mode = max(power for power in fourier if power >= 0)
    for degree in range(1, maximum_mode + 2):
        power = truncate_matrix(power * operator, parameter)
        coefficient = (
            fourier.get(0, 0)
            if degree == 1
            else 2 * fourier.get(degree - 1, 0)
        )
        correction = truncate_matrix(
            correction + coefficient * power,
            parameter,
        )
    return truncate_matrix(direction - correction, parameter)


def polynomial_pair(
    coefficients: list[sp.Expr],
    operator: sp.Matrix,
    direction: sp.Matrix,
    parameter: sp.Symbol,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Evaluate a polynomial and its Frechet derivative by Horner steps."""

    value = sp.zeros(DIMENSION)
    derivative = sp.zeros(DIMENSION)
    identity = sp.eye(DIMENSION)
    for coefficient in coefficients:
        derivative = truncate_matrix(
            derivative * operator + value * direction,
            parameter,
        )
        value = truncate_matrix(
            value * operator + coefficient * identity,
            parameter,
        )
    return value, derivative


def make_record() -> ExactAmbientGradientRecord:
    """Derive and validate the exact order separation."""

    parameter = sp.symbols("x", real=True)
    disk_direction = (
        sp.Integer(0),
        sp.Rational(37, 1000) + sp.I * sp.Rational(2, 125),
        sp.Rational(11, 250) + sp.I * sp.Rational(31, 1000),
    )
    operator, metric = disk_model(parameter, disk_direction)
    ambient = sp.zeros(DIMENSION)
    ambient[1, 2] = 1
    pulled = riemann_pulled_direction(
        operator,
        metric,
        ambient,
        parameter,
    )

    characteristic = [
        truncate(coefficient, parameter)
        for coefficient in operator.charpoly().all_coeffs()
    ]
    numerator_coefficients = characteristic[:-1]
    denominator_coefficients = [
        sp.conjugate(coefficient).subs(
            sp.conjugate(parameter),
            parameter,
        )
        for coefficient in numerator_coefficients[::-1]
    ]
    numerator, numerator_derivative = polynomial_pair(
        numerator_coefficients,
        operator,
        pulled,
        parameter,
    )
    denominator, denominator_derivative = polynomial_pair(
        denominator_coefficients,
        operator,
        pulled,
        parameter,
    )
    denominator_inverse = truncate_matrix(
        denominator.inv(),
        parameter,
    )
    blaschke = truncate_matrix(
        numerator * denominator_inverse,
        parameter,
    )
    blaschke_derivative = truncate_matrix(
        numerator_derivative * denominator_inverse
        - blaschke
        * denominator_derivative
        * denominator_inverse,
        parameter,
    )

    defect = metric * sp.eye(DIMENSION)[:, 0]
    endpoint = sp.eye(DIMENSION)[:, -1]
    gradient = truncate(
        8
        * sp.re(
            (
                sp.conjugate(defect).T
                * blaschke_derivative
                * endpoint
            )[0]
        ),
        parameter,
    )
    gradient_coefficient = sp.factor(
        gradient.coeff(parameter, 2)
    )

    coefficients = disk_direction[1:]
    norm_square = sum(
        sp.conjugate(value) * value
        for value in coefficients
    )
    reversal_pairing = sum(
        disk_direction[index] * disk_direction[LENGTH - index]
        for index in range(1, LENGTH)
    )
    quartic = sp.factor(
        norm_square**2
        - sp.conjugate(reversal_pairing) * reversal_pairing
    )
    if gradient_coefficient != -sp.Rational(424, 15625):
        raise AssertionError("the exact ambient gradient changed")
    if quartic != sp.Rational(25281, 15625000000):
        raise AssertionError("the exact quartic defect changed")

    return ExactAmbientGradientRecord(
        dimension=DIMENSION,
        length=LENGTH,
        ambient_direction="E_(1,2)",
        disk_direction=tuple(map(str, disk_direction[1:])),
        quartic_defect_coefficient=str(quartic),
        gradient_quadratic_coefficient=str(gradient_coefficient),
        gradient_series=str(sp.factor(gradient)),
        strict_order_separation=bool(
            gradient_coefficient != 0 and quartic > 0
        ),
    )


def main() -> None:
    """Print the exact deterministic record."""

    print(json.dumps(asdict(make_record()), sort_keys=True))


if __name__ == "__main__":
    main()
