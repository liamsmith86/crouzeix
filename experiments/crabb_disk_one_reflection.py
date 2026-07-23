#!/usr/bin/env python3
"""Audit the first elliptic reflection over a general disk point.

At ``c=0`` let ``g`` be the characteristic factor of L122's
Toeplitz disk operator.  The exact first ellipse pullback is

    T(c) = A + c (A^dagger_K - A^3) + O(c^2).

For the corrected Faber numerator, monic preparation gives

    N_c = g + c n_1 + O(c^2),
    n_1 = rem_g(
        w^3 g'(w) - (g'(w)-g'(0))/w
        + a_(L-1) w^(L-1)
    ).

The last term is L151's sharp raw grade-one Hardy correction.
This checker evaluates the resulting derivative of the simple top
Blaschke norm exactly over the rationals.

It falsifies exact one-reflection stationarity away from the equality
cone.  It also probes the corrected target suggested by the data:
the derivative vanishes quadratically in the disk-normal distance and
is therefore bounded by the Gram determinant ``Q``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_disk_exact import (
    ascending_coefficients,
    characteristic_factor,
    coefficient_model,
    evaluate_polynomial,
    frechet_polynomial,
)


@dataclass(frozen=True)
class OneReflectionCounterexample:
    """One exact nonzero general-disk first derivative."""

    length: int
    dimension: int
    coefficients: tuple[str, ...]
    derivative: str
    derivative_float: float
    derivative_is_nonzero: bool


@dataclass(frozen=True)
class NormalQuadraticRecord:
    """Normal blow-up ratios for one equality anchor."""

    length: int
    dimension: int
    equality_coefficients: tuple[str, ...]
    normal_direction: tuple[str, ...]
    normal_steps: tuple[str, ...]
    derivative_over_quartic: tuple[float, ...]
    maximum_absolute_ratio: float
    equality_derivative_is_zero: bool


def reflected_quartic(coefficients: tuple[sp.Expr, ...]) -> sp.Expr:
    """Return ``||z||^4-|z^T Jz|^2`` for real coefficients."""

    norm_square = sum(value**2 for value in coefficients)
    pairing = sum(
        coefficients[index] * coefficients[-1 - index]
        for index in range(len(coefficients))
    )
    return sp.factor(norm_square**2 - pairing**2)


def first_reflection_derivative(
    coefficients: tuple[sp.Expr, ...],
) -> sp.Expr:
    """Return the exact first derivative of the corrected dual square."""

    operator, coordinate, _, _ = coefficient_model(coefficients)
    length = len(coefficients) + 1
    dimension = length + 1
    variable = sp.symbols("w")
    factor = characteristic_factor(operator, variable)
    numerator_coefficients = ascending_coefficients(
        factor,
        variable,
        length,
    )

    natural_tangent = (
        variable**3 * sp.diff(factor, variable)
        - (
            sp.diff(factor, variable)
            - sp.diff(factor, variable).subs(variable, 0)
        )
        / variable
    )
    corrected_tangent = (
        natural_tangent
        + numerator_coefficients[length - 1]
        * variable ** (length - 1)
    )
    prepared_tangent = sp.Poly(
        corrected_tangent,
        variable,
    ).rem(sp.Poly(factor, variable)).as_expr()
    tangent_coefficients = ascending_coefficients(
        prepared_tangent,
        variable,
        length,
    )

    denominator_coefficients = list(
        reversed(numerator_coefficients)
    )
    tangent_denominator_coefficients = list(
        reversed(tangent_coefficients)
    )
    numerator = evaluate_polynomial(
        numerator_coefficients,
        operator,
    )
    denominator = evaluate_polynomial(
        denominator_coefficients,
        operator,
    )
    blaschke = sp.simplify(numerator * denominator.inv())

    operator_tangent = sp.simplify(
        coordinate.inv()
        * operator.T
        * coordinate
        - operator**3
    )
    numerator_tangent = (
        frechet_polynomial(
            numerator_coefficients,
            operator,
            operator_tangent,
        )
        + evaluate_polynomial(tangent_coefficients, operator)
    )
    denominator_tangent = (
        frechet_polynomial(
            denominator_coefficients,
            operator,
            operator_tangent,
        )
        + evaluate_polynomial(
            tangent_denominator_coefficients,
            operator,
        )
    )
    blaschke_tangent = sp.simplify(
        numerator_tangent * denominator.inv()
        - blaschke * denominator_tangent * denominator.inv()
    )

    top_right = sp.eye(dimension)[:, -1]
    numerator_scalar = (
        (blaschke * top_right).T
        * coordinate
        * (blaschke_tangent * top_right)
    )[0]
    denominator_scalar = (
        top_right.T * coordinate * top_right
    )[0]
    return sp.factor(
        2 * sp.re(numerator_scalar) / denominator_scalar
    )


def counterexample_record() -> OneReflectionCounterexample:
    """Return a compact exact failure of general stationarity."""

    coefficients = (sp.Rational(1, 20), sp.Rational(1, 30))
    derivative = first_reflection_derivative(coefficients)
    if derivative == 0:
        raise AssertionError("the one-reflection counterexample vanished")
    return OneReflectionCounterexample(
        length=3,
        dimension=4,
        coefficients=tuple(map(str, coefficients)),
        derivative=str(derivative),
        derivative_float=float(derivative),
        derivative_is_nonzero=True,
    )


def normal_record(
    equality_coefficients: tuple[sp.Expr, ...],
    normal_direction: tuple[sp.Expr, ...],
) -> NormalQuadraticRecord:
    """Probe exact quadratic normal vanishing at one equality point."""

    length = len(equality_coefficients) + 1
    equality_derivative = first_reflection_derivative(
        equality_coefficients
    )
    steps = tuple(
        sp.Rational(1, denominator)
        for denominator in (200, 400, 800, 1600)
    )
    ratios: list[float] = []
    for step in steps:
        coefficients = tuple(
            equality + step * normal
            for equality, normal in zip(
                equality_coefficients,
                normal_direction,
                strict=True,
            )
        )
        quartic = reflected_quartic(coefficients)
        derivative = first_reflection_derivative(coefficients)
        ratios.append(float(derivative / quartic))
    if equality_derivative != 0 or max(map(abs, ratios)) > 20:
        raise AssertionError("the normal quadratic probe failed")
    return NormalQuadraticRecord(
        length=length,
        dimension=length + 1,
        equality_coefficients=tuple(
            map(str, equality_coefficients)
        ),
        normal_direction=tuple(map(str, normal_direction)),
        normal_steps=tuple(map(str, steps)),
        derivative_over_quartic=tuple(ratios),
        maximum_absolute_ratio=max(map(abs, ratios)),
        equality_derivative_is_zero=True,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run and optionally persist the exact audit."""

    args = parse_args()
    records = [
        counterexample_record(),
        normal_record(
            (sp.Rational(1, 30), sp.Rational(1, 30)),
            (sp.Integer(1), sp.Integer(-1)),
        ),
        normal_record(
            (
                sp.Rational(1, 40),
                sp.Rational(-1, 50),
                sp.Rational(1, 40),
            ),
            (sp.Integer(1), sp.Integer(0), sp.Integer(-1)),
        ),
        normal_record(
            (
                sp.Rational(1, 50),
                sp.Rational(-1, 60),
                sp.Rational(-1, 60),
                sp.Rational(1, 50),
            ),
            (
                sp.Integer(1),
                sp.Integer(0),
                sp.Integer(0),
                sp.Integer(-1),
            ),
        ),
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    print("\n".join(lines))
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
