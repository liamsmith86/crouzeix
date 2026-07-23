#!/usr/bin/env python3
"""Regenerate the L121 endpoint selection on the disk-flat space.

For the Dickson descent polynomial ``P_{L,c}``, form

    F_c(Y) = (DP_{L,c}(A_c)[Y])_{L,0}
             - c**L (DP_{L,c}(A_c)[Y])_{0,L}.

L120 proves that the first possible power of ``c`` is ``floor(L/2)``.
This checker verifies symbolically that its coefficient is exactly one of
L65's two bottom-mode obstructions (or, in size three, the circular
tangent constraint).  It then imposes the corresponding flat condition
and confirms the improved divisibility by
``c**(floor(L/2)+1)``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_descent_gradient import (
    crabb_matrix,
    dickson_tangent,
    dickson_values,
)


@dataclass(frozen=True)
class FlatSelectionRecord:
    dimension: int
    length: int
    unrestricted_first_power: int
    flat_first_power: int
    required_flat_power: int
    leading_bottom_mode: str
    exact_endpoint_formula: bool
    leading_coefficient_identity: bool


def endpoint_formula(
    direction: sp.Matrix,
    parameter: sp.Symbol,
) -> sp.Expr:
    """Return the closed root-of-unity formula for ``F_c(direction)``."""

    dimension = direction.rows
    length = dimension - 1
    endpoint_weight = 1 / sp.sqrt(2)
    weights = [
        endpoint_weight,
        *[sp.Integer(1) for _ in range(dimension - 2)],
        endpoint_weight,
    ]
    result = 0
    for row in range(dimension):
        for column in range(dimension):
            grade = column - row - 1
            if grade % 2:
                continue
            coefficient = 2 * weights[row] * weights[column]
            if grade <= -2:
                power = length + grade // 2
                result += (
                    coefficient
                    * parameter**power
                    * direction[row, column]
                )
            elif grade >= 0:
                power = length + grade // 2
                result -= (
                    coefficient
                    * parameter**power
                    * direction[row, column]
                )
    return sp.expand(result)


def first_power(expression: sp.Expr, parameter: sp.Symbol) -> int:
    """Return the least nonzero parameter power in a polynomial."""

    polynomial = sp.Poly(sp.expand(expression), parameter)
    return min(
        monomial[0]
        for monomial, coefficient in polynomial.terms()
        if coefficient != 0
    )


def audit_dimension(dimension: int) -> FlatSelectionRecord:
    """Verify the exact endpoint and flat-space selection identities."""

    length = dimension - 1
    parameter = sp.symbols("c", positive=True)
    coordinates = sp.symbols(
        f"y0:{dimension * dimension}",
        complex=True,
    )
    direction = sp.Matrix(dimension, dimension, coordinates)
    crabb = crabb_matrix(dimension)
    axis = crabb + parameter * crabb.T
    values = dickson_values(axis, length, parameter)
    tangent = dickson_tangent(axis, direction, values, parameter)
    endpoint = sp.expand(
        tangent[length, 0]
        - parameter**length * tangent[0, length]
    )

    closed_endpoint = endpoint_formula(direction, parameter)
    exact_formula = sp.expand(endpoint - closed_endpoint) == 0
    unrestricted_power = first_power(endpoint, parameter)
    expected_unrestricted_power = length // 2

    if length % 2:
        bottom_entry = direction[length, 0]
        expected_leading = bottom_entry
        flat_endpoint = endpoint.subs(bottom_entry, 0)
        bottom_mode = "E[L,0]"
    else:
        first_entry = direction[length - 1, 0]
        second_entry = direction[length, 1]
        expected_leading = sp.sqrt(2) * (first_entry + second_entry)
        flat_endpoint = endpoint.subs(second_entry, -first_entry)
        bottom_mode = "sqrt(2)*(E[L-1,0]+E[L,1])"

    leading = sp.Poly(endpoint, parameter).coeff_monomial(
        parameter**expected_unrestricted_power
    )
    leading_identity = sp.expand(leading - expected_leading) == 0
    flat_power = first_power(flat_endpoint, parameter)
    required_flat_power = expected_unrestricted_power + 1

    if not (
        exact_formula
        and leading_identity
        and unrestricted_power == expected_unrestricted_power
        and flat_power >= required_flat_power
    ):
        raise RuntimeError(
            f"flat endpoint selection failed in size {dimension}"
        )

    return FlatSelectionRecord(
        dimension=dimension,
        length=length,
        unrestricted_first_power=unrestricted_power,
        flat_first_power=flat_power,
        required_flat_power=required_flat_power,
        leading_bottom_mode=bottom_mode,
        exact_endpoint_formula=exact_formula,
        leading_coefficient_identity=leading_identity,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=16)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the symbolic audit and optionally persist JSONL records."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")

    records = [
        audit_dimension(dimension)
        for dimension in range(args.minimum_size, args.maximum_size + 1)
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
