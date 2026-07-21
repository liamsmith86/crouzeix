#!/usr/bin/env python3
"""Exact completion certificate for the singular low-nome center.

This checker strengthens L35 in two ways.  It proves both cubic-envelope
endpoint residuals positive on

    0 < c <= 1/20,  |(p-p_star)/c**2| <= 50,

and it proves the narrower bridge tube

    0 < c <= 1/12,  |p-p_star| <= 4*c**4.

The first assertion is split into a deep Taylor tube
``p=p_star+c**4*x, |x|<=16`` and the annulus
``16*c**2<=|(p-p_star)/c**2|<=50``.  All arithmetic is rational and every
large polynomial is regenerated from the compact L32 residual formula.
The full audit is intentionally slow; no opaque expansion file is trusted.
"""

from __future__ import annotations

import gc
import time

import sympy as sp
from sympy import QQ
from sympy.polys.rings import PolyElement, ring

from slice_odd_tiny_edge_certificate import (
    assert_active_variables,
    bernstein_rectangle_lower,
    common_series_substitutions,
    compose,
    divide_by_c_power,
    homogenized_residual,
    minimum_c_power,
    theta_bounds,
)


SMALL_EDGE = QQ(1, 20)
BRIDGE_EDGE = QQ(1, 12)
UPPER_SCALE = QQ(9_895_604_649_984)
LOWER_SCALE = QQ(316_659_348_799_488)
LOWER_ANNULUS_SCALE = QQ(712_483_534_798_848)


def polynomial_ring() -> tuple[object, ...]:
    return ring("c,g,s,p,x,G,S,A,B,C,D,E,F,M,N,P,Q,U,V,W,Z", QQ)


def negative_term_lower_bound(
    polynomial: sp.Poly,
    edge: sp.Rational,
) -> sp.Rational:
    result = sp.Integer(0)
    for (degree,), coefficient in polynomial.terms():
        if degree == 0:
            result += coefficient
        elif coefficient < 0:
            result += coefficient * edge**degree
    return sp.factor(result)


def certify_positive_rational(
    expression: sp.Expr,
    c: sp.Symbol,
    vanishing_order: int,
    edge: sp.Rational,
) -> tuple[sp.Rational, sp.Rational]:
    numerator, denominator = sp.cancel(expression / c**vanishing_order).as_numer_denom()
    numerator_poly = sp.Poly(numerator, c)
    denominator_poly = sp.Poly(denominator, c)
    sample = edge / 2
    if numerator_poly.eval(sample) < 0:
        numerator_poly = -numerator_poly
    if denominator_poly.eval(sample) < 0:
        denominator_poly = -denominator_poly
    numerator_lower = negative_term_lower_bound(numerator_poly, edge)
    denominator_lower = negative_term_lower_bound(denominator_poly, edge)
    if numerator_lower <= 0 or denominator_lower <= 0:
        raise AssertionError((numerator_lower, denominator_lower))
    return numerator_lower, denominator_lower


def certify_theta_remainders() -> list[tuple[sp.Rational, sp.Rational]]:
    """Prove all deep theta bounds on the larger bridge edge."""

    c = sp.symbols("c", positive=True)
    bounds = theta_bounds(c)
    g_through_c10 = 2 - 4 * c**2 + 10 * c**4 - 20 * c**6 + 36 * c**8 - 64 * c**10
    s_through_c8 = 1 - 4 * c**2 + 12 * c**4 - 32 * c**6 + 76 * c**8
    g_through_c14 = g_through_c10 + 110 * c**12 - 180 * c**14
    s_through_c14 = s_through_c8 - 168 * c**10 + 352 * c**12 - 704 * c**14
    residuals = (
        (bounds.g_lower - (g_through_c10 - 120 * c**12), 12),
        (g_through_c10 + 120 * c**12 - bounds.g_upper, 12),
        (bounds.s_lower - (s_through_c8 - 200 * c**10), 10),
        (s_through_c8 + 200 * c**10 - bounds.s_upper, 10),
        (bounds.g_lower - (g_through_c14 - 300 * c**16), 16),
        (g_through_c14 + 300 * c**16 - bounds.g_upper, 16),
        (bounds.s_lower - (s_through_c14 - 1500 * c**16), 16),
        (s_through_c14 + 1500 * c**16 - bounds.s_upper, 16),
        (2 - bounds.g_upper, 2),
    )
    return [
        certify_positive_rational(value, c, order, BRIDGE_EDGE)
        for value, order in residuals
    ]


def truncated_expression(polynomial: PolyElement, maximum_c_power: int) -> sp.Expr:
    """Extract an even c truncation and rewrite it in t=c^2."""

    truncated = polynomial.ring.zero
    for monomial, coefficient in polynomial.terms():
        if monomial[0] <= maximum_c_power:
            if monomial[0] % 2:
                raise AssertionError(monomial)
            truncated[monomial] = coefficient
    t, x = sp.symbols("t x")
    return sp.expand(
        sp.sympify(truncated.as_expr()).subs(
            {sp.Symbol("c"): sp.sqrt(t), sp.Symbol("x"): x}
        )
    )


def bernstein_range_minimum(
    expression: sp.Expr,
    edge: sp.Rational,
    lower: int,
    upper: int,
) -> sp.Rational:
    """Lower-bound a two-variable polynomial on one x strip."""

    t, x = sp.symbols("t x")
    values = []
    left = QQ(lower)
    while left < upper:
        right = min(left + QQ(1, 8), QQ(upper))
        values.append(
            bernstein_rectangle_lower(
                expression,
                t,
                (QQ(0), edge**2),
                x,
                (left, right),
            )
        )
        left = right
    return min(values)


def symmetric_tail_bound(
    polynomial: PolyElement,
    first_power: int,
    edge: sp.Rational,
    x_bound: int,
    remainder_bounds: dict[int, int],
) -> object:
    """Absolute coefficient bound for a normalized even-c tail."""

    result = QQ(0)
    for monomial, coefficient in polynomial.terms():
        if monomial[0] < first_power:
            continue
        term = abs(coefficient) * edge ** monomial[0] * QQ(x_bound) ** monomial[4]
        for index, bound in remainder_bounds.items():
            term *= QQ(bound) ** monomial[index]
        result += term
    return result


def certify_upper_tubes(polynomial: PolyElement) -> tuple[tuple[float, ...], ...]:
    expression = truncated_expression(polynomial, 6)

    small_minima = tuple(
        bernstein_range_minimum(expression, SMALL_EDGE, lower, upper)
        for lower, upper in ((-2, 2), (-16, -2), (2, 16))
    )
    small_tail = symmetric_tail_bound(
        polynomial, 8, SMALL_EDGE, 16, {19: 300, 20: 1500}
    )
    if min(small_minima) <= QQ(4, 5) * UPPER_SCALE or small_tail >= QQ(1, 5) * UPPER_SCALE:
        raise AssertionError((small_minima, small_tail / UPPER_SCALE))

    bridge_minima = tuple(
        bernstein_range_minimum(expression, BRIDGE_EDGE, lower, upper)
        for lower, upper in ((-2, 2), (-4, -2), (2, 4))
    )
    bridge_central_tail = symmetric_tail_bound(
        polynomial, 8, BRIDGE_EDGE, 2, {19: 300, 20: 1500}
    )
    bridge_outer_tail = symmetric_tail_bound(
        polynomial, 8, BRIDGE_EDGE, 4, {19: 300, 20: 1500}
    )
    if (
        bridge_minima[0] <= QQ(1, 2) * UPPER_SCALE
        or min(bridge_minima[1:]) <= QQ(7, 5) * UPPER_SCALE
        or bridge_central_tail >= QQ(1, 2) * UPPER_SCALE
        or bridge_outer_tail >= QQ(6, 5) * UPPER_SCALE
    ):
        raise AssertionError(
            (
                bridge_minima,
                bridge_central_tail / UPPER_SCALE,
                bridge_outer_tail / UPPER_SCALE,
            )
        )
    return (
        tuple(float(value / UPPER_SCALE) for value in small_minima)
        + (float(small_tail / UPPER_SCALE),),
        tuple(float(value / UPPER_SCALE) for value in bridge_minima)
        + (
            float(bridge_central_tail / UPPER_SCALE),
            float(bridge_outer_tail / UPPER_SCALE),
        ),
    )


def certify_lower_tubes(polynomial: PolyElement) -> tuple[tuple[float, ...], ...]:
    expression = truncated_expression(polynomial, 4)

    small_minimum = bernstein_range_minimum(expression, SMALL_EDGE, -16, 16)
    small_tail = symmetric_tail_bound(
        polynomial, 6, SMALL_EDGE, 16, {14: 200, 15: 120}
    )
    if small_minimum <= QQ(3, 4) * LOWER_SCALE or small_tail >= QQ(1, 10) * LOWER_SCALE:
        raise AssertionError((small_minimum / LOWER_SCALE, small_tail / LOWER_SCALE))

    bridge_minimum = bernstein_range_minimum(expression, BRIDGE_EDGE, -4, 4)
    bridge_tail = symmetric_tail_bound(
        polynomial, 6, BRIDGE_EDGE, 4, {14: 200, 15: 120}
    )
    if bridge_minimum <= QQ(1, 2) * LOWER_SCALE or bridge_tail >= QQ(1, 4) * LOWER_SCALE:
        raise AssertionError((bridge_minimum / LOWER_SCALE, bridge_tail / LOWER_SCALE))
    return (
        (float(small_minimum / LOWER_SCALE), float(small_tail / LOWER_SCALE)),
        (float(bridge_minimum / LOWER_SCALE), float(bridge_tail / LOWER_SCALE)),
    )


def certify_upper_annulus(polynomial: PolyElement) -> sp.Rational:
    """Prove the upper endpoint when 16c^2<=|y|<=50."""

    base_keys = {(0, 2), (2, 1), (4, 0)}
    base_coefficients: dict[tuple[int, int], object] = {}
    correction = QQ(0)
    for monomial, coefficient in polynomial.terms():
        y_power = monomial[4]
        c_power = monomial[0] - 2 * y_power + 4
        if c_power < 0 or c_power % 2:
            raise AssertionError(monomial)
        if (
            (c_power, y_power) in base_keys
            and not monomial[19]
            and not monomial[20]
        ):
            base_coefficients[c_power, y_power] = coefficient
            continue
        half_power = c_power // 2
        term = (
            abs(coefficient)
            * QQ(300) ** monomial[19]
            * QQ(1500) ** monomial[20]
        )
        if y_power >= 2:
            term *= SMALL_EDGE ** (2 * half_power) * QQ(50) ** (y_power - 2)
        elif y_power == 1:
            if half_power < 1:
                raise AssertionError(monomial)
            term *= SMALL_EDGE ** (2 * (half_power - 1)) / 16
        else:
            if half_power < 2:
                raise AssertionError(monomial)
            term *= SMALL_EDGE ** (2 * (half_power - 2)) / 256
        correction += term

    expected_base = {
        (0, 2): 2 * UPPER_SCALE,
        (2, 1): 4 * UPPER_SCALE,
        (4, 0): 3 * UPPER_SCALE,
    }
    if base_coefficients != expected_base:
        raise AssertionError(base_coefficients)

    # K(2y^2+4c^2y+3c^4)=K(2(y+c^2)^2+c^4)
    # is at least (225/128)K y^2 on this annulus.
    if correction >= QQ(3, 2) * UPPER_SCALE or QQ(225, 128) <= QQ(3, 2):
        raise AssertionError(correction / UPPER_SCALE)
    return correction / UPPER_SCALE


def certify_lower_annulus(polynomial: PolyElement) -> sp.Rational:
    """Prove the lower endpoint when 16c^2<=|y|<=50."""

    correction = QQ(0)
    leading_coefficient = QQ(0)
    positive_constant = QQ(0)
    for monomial, coefficient in polynomial.terms():
        y_power = monomial[4]
        c_power = monomial[0] - 2 * y_power + 2
        if c_power < 0 or c_power % 2:
            raise AssertionError(monomial)
        if (
            (c_power, y_power) == (0, 2)
            and not monomial[14]
            and not monomial[15]
        ):
            leading_coefficient += coefficient
            continue
        # The c^2 constant is positive and may be discarded.
        if (
            (c_power, y_power) == (2, 0)
            and coefficient > 0
            and not monomial[14]
            and not monomial[15]
        ):
            positive_constant += coefficient
            continue
        half_power = c_power // 2
        term = (
            abs(coefficient)
            * QQ(200) ** monomial[14]
            * QQ(120) ** monomial[15]
        )
        if y_power >= 2:
            term *= SMALL_EDGE ** (2 * half_power) * QQ(50) ** (y_power - 2)
        elif y_power == 1:
            if half_power < 1:
                raise AssertionError(monomial)
            term *= SMALL_EDGE ** (2 * (half_power - 1)) / 16
        else:
            if half_power < 2:
                raise AssertionError(monomial)
            term *= SMALL_EDGE ** (2 * (half_power - 2)) / 256
        correction += term
    if (
        leading_coefficient != LOWER_ANNULUS_SCALE
        or positive_constant <= 0
        or correction >= LOWER_ANNULUS_SCALE
    ):
        raise AssertionError(correction / LOWER_ANNULUS_SCALE)
    return correction / LOWER_ANNULUS_SCALE


def generate_upper(data: tuple[object, ...], started: float) -> PolyElement:
    (
        _, c, _, _, _, _, _, _, _, _, c_rem, d_rem, e_rem, f_rem,
        m_rem, n_rem, p_rem, q_rem, u_rem, v_rem, w_rem, z_rem,
    ) = data
    polynomial = homogenized_residual(data, "upper")
    print(f"upper homogenized: {len(polynomial.terms())} terms", flush=True)
    polynomial = common_series_substitutions(data, polynomial, started)
    for variable, value in (
        (c_rem, -20 + c**2 * e_rem),
        (d_rem, -32 + c**2 * f_rem),
        (e_rem, 36 + c**2 * m_rem),
        (f_rem, 76 + c**2 * n_rem),
    ):
        polynomial = compose(polynomial, variable, value, started)
    if minimum_c_power(polynomial) != 10:
        raise AssertionError("unexpected upper vanishing order")
    polynomial = divide_by_c_power(polynomial, 10)
    for variable, value in (
        (m_rem, -64 + c**2 * p_rem),
        (n_rem, -168 + c**2 * q_rem),
        (p_rem, 110 + c**2 * u_rem),
        (q_rem, 352 + c**2 * v_rem),
        (u_rem, -180 + c**2 * w_rem),
        (v_rem, -704 + c**2 * z_rem),
    ):
        polynomial = compose(polynomial, variable, value, started)
    assert_active_variables(polynomial, {0, 4, 19, 20})
    return polynomial


def generate_lower(data: tuple[object, ...], started: float) -> PolyElement:
    (
        _, c, _, _, _, _, _, _, _, _, c_rem, d_rem, e_rem, f_rem,
        m_rem, n_rem, p_rem, *_
    ) = data
    polynomial = homogenized_residual(data, "lower")
    print(f"lower homogenized: {len(polynomial.terms())} terms", flush=True)
    polynomial = common_series_substitutions(data, polynomial, started)
    for variable, value in (
        (c_rem, -20 + c**2 * e_rem),
        (d_rem, -32 + c**2 * f_rem),
        (e_rem, 36 + c**2 * m_rem),
        (f_rem, 76 + c**2 * n_rem),
        (m_rem, -64 + c**2 * p_rem),
    ):
        polynomial = compose(polynomial, variable, value, started)
    if minimum_c_power(polynomial) != 8:
        raise AssertionError("unexpected lower vanishing order")
    polynomial = divide_by_c_power(polynomial, 8)
    assert_active_variables(polynomial, {0, 4, 14, 15})
    return polynomial


def main() -> None:
    started = time.monotonic()
    theta_checks = certify_theta_remainders()
    print("deep theta remainder bounds: exact", flush=True)
    data = polynomial_ring()

    upper = generate_upper(data, started)
    print(f"upper deep series: {len(upper.terms())} terms", flush=True)
    upper_tubes = certify_upper_tubes(upper)
    upper_annulus = certify_upper_annulus(upper)
    del upper
    gc.collect()

    lower = generate_lower(data, started)
    print(f"lower deep series: {len(lower.terms())} terms", flush=True)
    lower_tubes = certify_lower_tubes(lower)
    lower_annulus = certify_lower_annulus(lower)

    print("centered low-nome completion certificate: exact")
    print(f"  upper small/bridge tube data: {upper_tubes}")
    print(f"  lower small/bridge tube data: {lower_tubes}")
    print(f"  annulus correction ratios: ({float(upper_annulus):.12g}, {float(lower_annulus):.12g})")
    print(f"  theta rational checks: {len(theta_checks)}")
    print(f"  elapsed seconds: {time.monotonic() - started:.1f}")


if __name__ == "__main__":
    main()
