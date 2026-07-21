#!/usr/bin/env python3
"""Exact certificates for the sharp nome bound and the rank-one node face.

The proof is recorded in ``proof/slice_coupled_defects.md``.  This checker
regenerates the finite Bernstein certificates used there:

* ``k <= 4c/(1+4c^2)`` on ``0 <= c <= 1/2``;
* the ``p=0`` block-energy inequality on ``0 <= c <= 1/2`` in exact charts
  centered at its curved equality ridge and its degenerate origin;
* the same inequality on ``1/2 <= c <= 1`` after the elementary replacement
  ``k <= 1``.

All coefficients are exact SymPy rationals.  No floating-point table is
loaded or trusted.
"""

from __future__ import annotations

from itertools import product
from typing import TypeAlias

import sympy as sp


Index: TypeAlias = tuple[int, ...]
BernsteinTensor: TypeAlias = dict[Index, sp.Rational]


def bernstein_coefficients(
    polynomial: sp.Poly,
) -> BernsteinTensor:
    """Return tensor-product Bernstein coefficients on the unit cube."""

    variables = polynomial.gens
    degrees = polynomial.degree_list()
    indices = tuple(product(*(range(degree + 1) for degree in degrees)))
    power = {
        index: polynomial.coeff_monomial(
            sp.prod(variable**exponent for variable, exponent in zip(variables, index))
        )
        for index in indices
    }
    coefficients = power
    for axis, degree in enumerate(degrees):
        transformed: dict[tuple[int, ...], sp.Rational] = {}
        for index in indices:
            transformed[index] = sum(
                sp.Rational(
                    sp.binomial(index[axis], exponent), sp.binomial(degree, exponent)
                )
                * coefficients[
                    tuple(
                        exponent if coordinate == axis else index[coordinate]
                        for coordinate in range(len(variables))
                    )
                ]
                for exponent in range(index[axis] + 1)
            )
        coefficients = transformed
    return coefficients


def split_bernstein_tensor(
    coefficients: BernsteinTensor,
    degrees: tuple[int, ...],
    axis: int,
) -> tuple[BernsteinTensor, BernsteinTensor]:
    """Bisect one tensor-product Bernstein axis exactly at one half."""

    left: BernsteinTensor = {}
    right: BernsteinTensor = {}
    other_axes = tuple(
        coordinate for coordinate in range(len(degrees)) if coordinate != axis
    )
    for other_index in product(
        *(range(degrees[coordinate] + 1) for coordinate in other_axes)
    ):

        def index_at(position: int) -> Index:
            index = [0] * len(degrees)
            index[axis] = position
            for coordinate, value in zip(other_axes, other_index):
                index[coordinate] = value
            return tuple(index)

        levels = [
            [coefficients[index_at(position)] for position in range(degrees[axis] + 1)]
        ]
        for _ in range(degrees[axis]):
            previous = levels[-1]
            levels.append(
                [
                    (previous[position] + previous[position + 1]) / 2
                    for position in range(len(previous) - 1)
                ]
            )
        for level in range(degrees[axis] + 1):
            left[index_at(level)] = levels[level][0]
            right[index_at(degrees[axis] - level)] = levels[level][-1]
    return left, right


def logarithmic_size(value: sp.Rational) -> int:
    """Return a cheap exact proxy for log2(abs(value))."""

    if value == 0:
        return -(10**9)
    rational = sp.Rational(value)
    return int(rational.p).bit_length() - int(rational.q).bit_length()


def certify_by_subdivision(polynomial: sp.Poly, label: str) -> None:
    """Certify a polynomial by exact adaptive Bernstein subdivision."""

    degrees = polynomial.degree_list()
    root = bernstein_coefficients(polynomial)
    stack: list[tuple[BernsteinTensor, tuple[int, ...]]] = [(root, (0,) * len(degrees))]
    splits = 0
    proved = 0
    maximum_depth = 0
    while stack:
        coefficients, depths = stack.pop()
        if min(coefficients.values()) >= 0:
            proved += 1
            continue
        scores: list[int] = []
        for axis, degree in enumerate(degrees):
            largest = -(10**9)
            for index, value in coefficients.items():
                if index[axis] == degree:
                    continue
                neighbour = list(index)
                neighbour[axis] += 1
                largest = max(
                    largest,
                    logarithmic_size(coefficients[tuple(neighbour)] - value),
                )
            scores.append(largest - depths[axis])
        axis = max(range(len(degrees)), key=scores.__getitem__)
        left, right = split_bernstein_tensor(coefficients, degrees, axis)
        child_depths = list(depths)
        child_depths[axis] += 1
        child_depths_tuple = tuple(child_depths)
        maximum_depth = max(maximum_depth, child_depths[axis])
        stack.extend(((left, child_depths_tuple), (right, child_depths_tuple)))
        splits += 1
        if splits > 1_000:
            raise AssertionError(f"{label}: subdivision did not terminate")
    print(
        f"{label}: exact Bernstein subdivision "
        f"({proved} leaves, {splits} splits, depth {maximum_depth})"
    )


def assert_nonnegative_bernstein(polynomial: sp.Poly, label: str) -> None:
    coefficients = bernstein_coefficients(polynomial)
    negative = [coefficient for coefficient in coefficients.values() if coefficient < 0]
    if negative:
        raise AssertionError(f"{label}: found {len(negative)} negative coefficients")
    print(
        f"{label}: exact Bernstein certificate "
        f"({len(coefficients)} coefficients, {sum(value == 0 for value in coefficients.values())} zero)"
    )


def sharp_nome_certificate() -> None:
    """Prove the two-factor Jacobi-product estimate exactly."""

    q, coordinate = sp.symbols("q coordinate", nonnegative=True)
    denominator = (1 + q) ** 4 * (1 + q**3) ** 4
    numerator = (1 + 4 * q) * (1 + q**2) ** 4 * (1 + q**4) ** 4
    residual = sp.cancel((denominator - numerator) / q**2)
    scaled = sp.Poly(sp.expand(residual.subs(q, coordinate / 4)), coordinate)
    assert_nonnegative_bernstein(scaled, "sharp nome residual")


def parameter_functions(parameter: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    schur = 3 * parameter / (4 - parameter**2)
    weight = 2 * (1 - parameter**2) / (4 - parameter**2)
    return schur, weight


def rank_one_node_residual(
    c: sp.Expr,
    odd_parameter: sp.Symbol,
    even_parameter: sp.Symbol,
    modulus: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the numerator and positive denominator of the p=0 residual."""

    odd_schur, odd_weight = parameter_functions(odd_parameter)
    even_schur, even_weight = parameter_functions(even_parameter)
    norm_sum = modulus * (1 + c**2) / c
    energy_numerator = (
        odd_weight * even_weight * norm_sum
        + odd_schur**2 * even_weight**2 * modulus**2
        + even_schur**2 * odd_weight**2 * norm_sum**2 / 4
    )
    residual = (1 - odd_schur * even_schur * modulus) ** 2 - energy_numerator
    return sp.cancel(residual).as_numer_denom()


def high_c_rank_one_node_certificate() -> None:
    """Prove the p=0 block-energy target for c >= 1/2."""

    coordinate, odd_parameter, even_parameter = sp.symbols(
        "coordinate odd_parameter even_parameter", nonnegative=True
    )
    c = (1 + coordinate) / 2
    # The residual decreases as k increases, so k <= 1 permits k=1.
    residual, denominator = rank_one_node_residual(
        c, odd_parameter, even_parameter, sp.Integer(1)
    )
    expected_denominator = (
        4
        * (odd_parameter - 2) ** 2
        * (odd_parameter + 2) ** 2
        * (even_parameter - 2) ** 2
        * (even_parameter + 2) ** 2
        * (coordinate + 1) ** 2
    )
    if sp.factor(denominator - expected_denominator) != 0:
        raise AssertionError("unexpected denominator in the p=0 residual")
    polynomial = sp.Poly(sp.expand(residual), coordinate, odd_parameter, even_parameter)
    assert_nonnegative_bernstein(polynomial, "p=0, c>=1/2 residual")


def transformed_polynomial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol, sp.Symbol],
    substitutions: dict[sp.Symbol, sp.Expr],
    divisor: sp.Expr = sp.Integer(1),
) -> sp.Poly:
    """Apply a simultaneous chart substitution and divide a known factor."""

    transformed = sp.cancel(
        sp.expand(expression.subs(substitutions, simultaneous=True)) / divisor
    )
    if sp.denom(transformed) != 1:
        raise AssertionError("chart division did not produce a polynomial")
    return sp.Poly(sp.expand(transformed), *variables)


def low_c_rank_one_node_certificate() -> None:
    """Prove the p=0 block-energy target for c <= 1/2."""

    coordinate, odd_parameter, even_parameter = sp.symbols(
        "coordinate odd_parameter even_parameter", nonnegative=True
    )
    variables = (coordinate, odd_parameter, even_parameter)
    c = coordinate / 2
    sharp_modulus_bound = 2 * coordinate / (1 + coordinate**2)
    residual, denominator = rank_one_node_residual(
        c, odd_parameter, even_parameter, sharp_modulus_bound
    )
    expected_denominator = (
        (odd_parameter - 2) ** 2
        * (odd_parameter + 2) ** 2
        * (even_parameter - 2) ** 2
        * (even_parameter + 2) ** 2
        * (coordinate**2 + 1) ** 2
    )
    if sp.factor(denominator - expected_denominator) != 0:
        raise AssertionError("unexpected denominator in the low-c p=0 residual")
    if sp.Poly(residual, *variables).degree_list() != (4, 4, 4):
        raise AssertionError("unexpected degree in the low-c p=0 residual")

    sixteenth = sp.Rational(1, 16)
    fifteen_sixteenths = sp.Rational(15, 16)
    charts: tuple[tuple[str, dict[sp.Symbol, sp.Expr], sp.Expr], ...] = (
        # b >= 15/16, split at the exact ridge a=coordinate.
        (
            "low-c ridge below",
            {
                even_parameter: fifteen_sixteenths + even_parameter / 16,
                odd_parameter: coordinate * odd_parameter,
            },
            sp.Integer(1),
        ),
        (
            "low-c ridge above",
            {
                even_parameter: fifteen_sixteenths + even_parameter / 16,
                odd_parameter: coordinate + (1 - coordinate) * odd_parameter,
            },
            sp.Integer(1),
        ),
        # b <= 15/16 and coordinate >= 1/16.
        (
            "low-c compact",
            {
                coordinate: sixteenth + fifteen_sixteenths * coordinate,
                even_parameter: fifteen_sixteenths * even_parameter,
            },
            sp.Integer(1),
        ),
        # coordinate <= 1/16, with b or a bounded away from zero.
        (
            "low-c b-middle",
            {
                coordinate: coordinate / 16,
                even_parameter: sixteenth + sp.Rational(14, 16) * even_parameter,
            },
            sp.Integer(1),
        ),
        (
            "low-c a-middle",
            {
                coordinate: coordinate / 16,
                even_parameter: even_parameter / 16,
                odd_parameter: sixteenth + fifteen_sixteenths * odd_parameter,
            },
            sp.Integer(1),
        ),
        # In the remaining tiny cube choose its largest coordinate.  The
        # residual has a common square of that leading coordinate.
        (
            "low-c origin x-chart",
            {
                coordinate: coordinate / 16,
                odd_parameter: coordinate * odd_parameter / 16,
                even_parameter: coordinate * even_parameter / 16,
            },
            coordinate**2,
        ),
        (
            "low-c origin a-chart",
            {
                coordinate: odd_parameter * coordinate / 16,
                odd_parameter: odd_parameter / 16,
                even_parameter: odd_parameter * even_parameter / 16,
            },
            odd_parameter**2,
        ),
        (
            "low-c origin b-chart",
            {
                coordinate: even_parameter * coordinate / 16,
                odd_parameter: even_parameter * odd_parameter / 16,
                even_parameter: even_parameter / 16,
            },
            even_parameter**2,
        ),
    )
    for label, substitutions, divisor in charts:
        polynomial = transformed_polynomial(residual, variables, substitutions, divisor)
        certify_by_subdivision(polynomial, label)


def main() -> None:
    sharp_nome_certificate()
    low_c_rank_one_node_certificate()
    high_c_rank_one_node_certificate()
    print("sharp nome and rank-one node certificates: PASS")


if __name__ == "__main__":
    main()
