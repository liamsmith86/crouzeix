#!/usr/bin/env python3
"""Regenerate the lowest Newton edge of the inverse ellipse map.

For ``0 < c < 1``, let ``Psi_c`` be the normalized conformal map from
the unit disk onto

    {w + c conjugate(w) : |w| < 1}.

Writing

    Psi_c(z) = sum_{n >= 0} b_n(c) z^(2n+1),

L125 proves ``b_n(c) = Catalan(n) c^n + O(c^(n+2))``.  This script
derives the coefficients from the exact differential equation for
``Psi_c``.  Jacobi theta series are used for the elliptic modulus, so
all arithmetic in the regeneration is rational and exact.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class InverseMapCoefficient:
    index: int
    catalan_coefficient: int
    computed_leading_coefficient: int
    next_coefficient: int


def truncate(expression: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    """Return an expanded Taylor polynomial through degree ``order - 1``."""

    return sp.series(expression, variable, 0, order).removeO().expand()


def theta_data(
    c: sp.Symbol,
    order: int,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return exact series for the modulus, period factor, and linear term.

    The ellipse nome is ``q=c^2``.  If ``k`` is Jacobi's modulus and
    ``alpha=pi/(2K)``, the theta identities give

        k = 4 c a(c)^2 / theta_3(c^2)^2,
        alpha = 1 / theta_3(c^2)^2,
        b_0 = 1 / (a(c) theta_3(c^2)).
    """

    theta_two_factor = sum(
        c ** (2 * index * (index + 1))
        for index in range(order)
        if 2 * index * (index + 1) < order
    )
    theta_three = 1 + 2 * sum(
        c ** (2 * index * index)
        for index in range(1, order)
        if 2 * index * index < order
    )
    modulus = truncate(
        4 * c * theta_two_factor**2 / theta_three**2,
        c,
        order,
    )
    period_factor = truncate(1 / theta_three**2, c, order)
    linear_coefficient = truncate(
        1 / (theta_two_factor * theta_three),
        c,
        order,
    )
    return modulus, period_factor, linear_coefficient


def convolution(
    coefficients: list[sp.Expr],
    degree: int,
    c: sp.Symbol,
    order: int,
    derivative_weights: bool = False,
) -> sp.Expr:
    """Return one coefficient of ``h^2`` or ``(h+2x h')^2``."""

    total = 0
    for left_index in range(degree + 1):
        right_index = degree - left_index
        if left_index >= len(coefficients) or right_index >= len(coefficients):
            continue
        left = coefficients[left_index]
        right = coefficients[right_index]
        if derivative_weights:
            left *= 2 * left_index + 1
            right *= 2 * right_index + 1
        total += left * right
    return truncate(total, c, order)


def inverse_map_coefficients(
    maximum_degree: int,
    audit_order: int,
) -> tuple[sp.Symbol, sp.Expr, sp.Expr, list[sp.Expr]]:
    """Solve the exact inverse-map ODE coefficient by coefficient."""

    c = sp.symbols("c")
    # Each recurrence divides by a series with valuation one.  Extra
    # internal precision prevents the requested tail from being polluted.
    work_order = audit_order + 2 * maximum_degree + 4
    modulus, period_factor, linear_coefficient = theta_data(c, work_order)
    coefficients = [linear_coefficient]

    for degree in range(1, maximum_degree + 1):
        padded = [*coefficients, sp.Integer(0)]
        current_derivative_square = convolution(
            padded,
            degree,
            c,
            work_order,
            derivative_weights=True,
        )
        previous_derivative_square = convolution(
            coefficients,
            degree - 1,
            c,
            work_order,
            derivative_weights=True,
        )
        second_previous_derivative_square = (
            convolution(
                coefficients,
                degree - 2,
                c,
                work_order,
                derivative_weights=True,
            )
            if degree >= 2
            else sp.Integer(0)
        )
        previous_square = convolution(
            coefficients,
            degree - 1,
            c,
            work_order,
        )

        numerator = truncate(
            (1 + modulus**2) * previous_derivative_square
            - modulus * second_previous_derivative_square
            - period_factor**2 * previous_square
            - modulus * current_derivative_square,
            c,
            work_order,
        )
        denominator = truncate(
            2
            * modulus
            * linear_coefficient
            * (2 * degree + 1),
            c,
            work_order,
        )
        coefficients.append(
            truncate(numerator / denominator, c, work_order)
        )

    return c, modulus, period_factor, [
        truncate(coefficient, c, audit_order)
        for coefficient in coefficients
    ]


def validate_differential_equation(
    c: sp.Symbol,
    modulus: sp.Expr,
    period_factor: sp.Expr,
    coefficients: list[sp.Expr],
    audit_order: int,
) -> None:
    """Check the defining ODE modulo the regenerated truncation."""

    x = sp.symbols("x")
    h = sum(
        coefficient * x**index
        for index, coefficient in enumerate(coefficients)
    )
    derivative_factor = h + 2 * x * sp.diff(h, x)
    residual = (
        derivative_factor**2
        * (modulus - (1 + modulus**2) * x + modulus * x**2)
        - period_factor**2 * (4 * c - x * h**2)
    )
    for degree in range(len(coefficients)):
        coefficient = sp.expand(residual).coeff(x, degree)
        if truncate(coefficient, c, audit_order - 1) != 0:
            raise AssertionError(
                f"inverse-map ODE failed at x-degree {degree}: {coefficient}"
            )


def make_records(
    maximum_degree: int,
    audit_order: int,
) -> list[InverseMapCoefficient]:
    """Regenerate, validate, and summarize the Newton-edge coefficients."""

    c, modulus, period_factor, coefficients = inverse_map_coefficients(
        maximum_degree,
        audit_order,
    )
    validate_differential_equation(
        c,
        modulus,
        period_factor,
        coefficients,
        audit_order,
    )

    records = []
    for index, coefficient in enumerate(coefficients):
        quotient = sp.cancel(coefficient / c**index)
        leading = sp.expand(quotient).coeff(c, 0)
        next_coefficient = sp.expand(quotient).coeff(c, 2)
        catalan = sp.catalan(index)
        if leading != catalan:
            raise AssertionError(
                f"degree {index}: expected Catalan {catalan}, got {leading}"
            )
        for power in range(index):
            if sp.expand(coefficient).coeff(c, power) != 0:
                raise AssertionError(
                    f"degree {index}: unexpected lower power c^{power}"
                )
        for power in range(audit_order):
            if (power - index) % 2 and sp.expand(coefficient).coeff(c, power):
                raise AssertionError(
                    f"degree {index}: parity failed at c^{power}"
                )
        records.append(
            InverseMapCoefficient(
                index=index,
                catalan_coefficient=int(catalan),
                computed_leading_coefficient=int(leading),
                next_coefficient=int(next_coefficient),
            )
        )
    return records


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-degree", type=int, default=7)
    parser.add_argument("--audit-order", type=int, default=18)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact regeneration and optionally persist JSONL output."""

    args = parse_args()
    if args.maximum_degree < 0:
        raise ValueError("maximum degree must be nonnegative")
    if args.audit_order <= args.maximum_degree + 2:
        raise ValueError("audit order is too small for the requested degrees")

    records = make_records(args.maximum_degree, args.audit_order)
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
