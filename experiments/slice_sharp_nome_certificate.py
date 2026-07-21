#!/usr/bin/env python3
"""Exact certificates for the sharp nome bound and the rank-one node face.

The proof is recorded in ``proof/slice_coupled_defects.md``.  This checker
regenerates the two finite Bernstein certificates used there:

* ``k <= 4c/(1+4c^2)`` on ``0 <= c <= 1/2``;
* the ``p=0`` block-energy inequality on ``1/2 <= c <= 1`` after the
  elementary replacement ``k <= 1``.

All coefficients are exact SymPy rationals.  No floating-point table is
loaded or trusted.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


def bernstein_coefficients(
    polynomial: sp.Poly,
) -> dict[tuple[int, ...], sp.Rational]:
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


def main() -> None:
    sharp_nome_certificate()
    high_c_rank_one_node_certificate()
    print("sharp nome and rank-one node certificates: PASS")


if __name__ == "__main__":
    main()
