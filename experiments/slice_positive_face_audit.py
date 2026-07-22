#!/usr/bin/env python3
"""Exact audit for the sharp positive-determinant ``a=b=1`` face.

The two projective determinant charts become manifest even powers on this
face.  Their first coefficients in ``A=1-a`` and ``B=1-b`` also factor with a
fixed nonnegative sign.  This script regenerates both identities directly
from the integer core records and audits the l1 norm of all terms of total
``(A,B)`` degree at least two.
"""

from __future__ import annotations

from math import comb, prod
from typing import TypeAlias

from slice_projective_core import load_records
from slice_projective_interval_certificate import ChartTable, final_chart_tables


ExactPolynomial: TypeAlias = dict[tuple[int, ...], int]
VARIABLE_COUNT = 8  # X, R, Y, c, k, s, gamma_lower, gamma_width
EXPECTED_REMAINDER_L1 = 41_235_531_913


def variable(index: int) -> ExactPolynomial:
    exponent = [0] * VARIABLE_COUNT
    exponent[index] = 1
    return {tuple(exponent): 1}


def add(*polynomials: ExactPolynomial) -> ExactPolynomial:
    output: ExactPolynomial = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            output[monomial] = output.get(monomial, 0) + coefficient
    return {
        monomial: coefficient for monomial, coefficient in output.items() if coefficient
    }


def scale(polynomial: ExactPolynomial, scalar: int) -> ExactPolynomial:
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in polynomial.items()
        if scalar * coefficient
    }


def multiply(*polynomials: ExactPolynomial) -> ExactPolynomial:
    output: ExactPolynomial = {(0,) * VARIABLE_COUNT: 1}
    for polynomial in polynomials:
        product: ExactPolynomial = {}
        for left_monomial, left_coefficient in output.items():
            for right_monomial, right_coefficient in polynomial.items():
                monomial = tuple(
                    left + right for left, right in zip(left_monomial, right_monomial)
                )
                product[monomial] = (
                    product.get(monomial, 0) + left_coefficient * right_coefficient
                )
        output = {
            monomial: coefficient
            for monomial, coefficient in product.items()
            if coefficient
        }
    return output


def power(polynomial: ExactPolynomial, exponent: int) -> ExactPolynomial:
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    output: ExactPolynomial = {(0,) * VARIABLE_COUNT: 1}
    base = polynomial
    while exponent:
        if exponent & 1:
            output = multiply(output, base)
        exponent //= 2
        if exponent:
            base = multiply(base, base)
    return output


def expected_face(chart: int) -> tuple[ExactPolynomial, ExactPolynomial]:
    """Return the face and common first-deficit coefficient for one chart."""

    one: ExactPolynomial = {(0,) * VARIABLE_COUNT: 1}
    x, ratio, y, c, k, s, lower, width = (
        variable(index) for index in range(VARIABLE_COUNT)
    )
    q = add(lower, multiply(y, width))
    one_minus_k = add(one, scale(k, -1))
    if chart == 0:
        z = x
        conformal = add(s, multiply(x, q))
        bracket = add(
            ratio,
            multiply(
                add(one, scale(multiply(ratio, x), -1)),
                power(conformal, 2),
            ),
        )
    else:
        z = multiply(ratio, x)
        conformal = add(s, multiply(z, q))
        bracket = add(
            one,
            multiply(ratio, add(one, scale(x, -1)), power(conformal, 2)),
        )
    one_minus_kz = add(one, scale(multiply(k, z), -1))
    face = scale(
        multiply(
            power(c, 9),
            power(one_minus_k, 2),
            power(one_minus_kz, 2),
            power(bracket, 4),
        ),
        81,
    )
    positive_linear_factor = add(
        scale(one, 2),
        scale(k, 3),
        scale(multiply(z, k), 3),
        scale(multiply(z, power(k, 2)), -8),
    )
    first_deficit = scale(
        multiply(
            power(c, 9),
            one_minus_k,
            one_minus_kz,
            positive_linear_factor,
            power(bracket, 4),
        ),
        54,
    )
    return face, first_deficit


def collected_deficit(table: ChartTable, a_power: int, b_power: int) -> ExactPolynomial:
    """Collect one coefficient after ``a=1-A`` and ``b=1-B``."""

    output: ExactPolynomial = {}
    sign = (-1) ** (a_power + b_power)
    for monomial, coefficient in table.records.items():
        if monomial[3] < a_power or monomial[4] < b_power:
            continue
        index = (monomial[0], monomial[1], monomial[2], *monomial[5:])
        value = (
            coefficient * comb(monomial[3], a_power) * comb(monomial[4], b_power) * sign
        )
        output[index] = output.get(index, 0) + value
    return {
        monomial: coefficient for monomial, coefficient in output.items() if coefficient
    }


def evaluate_modular(
    polynomial: ExactPolynomial, values: tuple[int, ...], modulus: int
) -> int:
    """Evaluate an integer sparse polynomial in one finite field."""

    return (
        sum(
            coefficient
            * prod(
                pow(value, exponent, modulus)
                for value, exponent in zip(values, monomial)
            )
            for monomial, coefficient in polynomial.items()
        )
        % modulus
    )


def audit_deficit_reconstruction(
    table: ChartTable,
    deficit_maps: dict[tuple[int, int], ExactPolynomial],
) -> None:
    """Cross-check the full binomial transform at two finite-field points."""

    original_values = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
    factor_values = (*original_values[:3], *original_values[5:])
    for modulus in (1_000_000_007, 2_147_483_647):
        expected = evaluate_modular(table.records, original_values, modulus)
        a_deficit = (1 - original_values[3]) % modulus
        b_deficit = (1 - original_values[4]) % modulus
        actual = (
            sum(
                evaluate_modular(records, factor_values, modulus)
                * pow(a_deficit, a_power, modulus)
                * pow(b_deficit, b_power, modulus)
                for (a_power, b_power), records in deficit_maps.items()
            )
            % modulus
        )
        if actual != expected:
            raise AssertionError(
                f"{table.label}: deficit reconstruction failed modulo {modulus}"
            )


def audit_table(table: ChartTable, chart: int) -> None:
    face, first_deficit = expected_face(chart)
    deficit_maps = {
        (a_power, b_power): collected_deficit(table, a_power, b_power)
        for a_power in range(5)
        for b_power in range(5)
    }
    actual_face = deficit_maps[0, 0]
    actual_a = deficit_maps[1, 0]
    actual_b = deficit_maps[0, 1]
    if actual_face != face:
        raise AssertionError(f"{table.label}: face factorization failed")
    if actual_a != first_deficit or actual_b != first_deficit:
        raise AssertionError(f"{table.label}: first-deficit factorization failed")

    remainder_l1 = 0
    for a_power in range(5):
        for b_power in range(5):
            if a_power + b_power >= 2:
                remainder_l1 += sum(
                    abs(coefficient)
                    for coefficient in deficit_maps[a_power, b_power].values()
                )
    if remainder_l1 != EXPECTED_REMAINDER_L1:
        raise AssertionError(
            f"{table.label}: remainder l1 {remainder_l1} != {EXPECTED_REMAINDER_L1}"
        )
    audit_deficit_reconstruction(table, deficit_maps)
    print(
        f"PASS {table.label}: face_terms={len(face)}, "
        f"first_terms={len(first_deficit)}, remainder_l1={remainder_l1}"
    )


def main() -> None:
    tables = final_chart_tables(load_records(), 1)[:2]
    for chart, table in enumerate(tables):
        audit_table(table, chart)


if __name__ == "__main__":
    main()
