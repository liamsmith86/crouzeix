#!/usr/bin/env python3
"""Exact certificate for the tiny-nome ridge in the upper odd block.

This checker proves the two cubic-envelope endpoint inequalities from
``proof/slice_odd_block_reduction.md`` when

    0 < c <= 1/20,  p = p_star + c**4*x,  |x| <= 4.

Only rational polynomial arithmetic is used.  The large intermediate
polynomials are generated from the compact residual formula and are never
stored as opaque certificate data.  A full run takes several minutes.
"""

from __future__ import annotations

import gc
import time
from dataclasses import dataclass

import sympy as sp
from sympy import QQ
from sympy.polys.rings import PolyElement, ring


EDGE = sp.Rational(1, 20)


@dataclass(frozen=True)
class ThetaBounds:
    g_lower: sp.Expr
    g_upper: sp.Expr
    s_lower: sp.Expr
    s_upper: sp.Expr


def theta_bounds(c: sp.Symbol) -> ThetaBounds:
    """Return rational enclosures for g=2R/T and s=T^-2 on the edge."""

    theta_lower = 1 + 2 * c**2 + 2 * c**8
    theta_upper = theta_lower + 2 * c**18 / (1 - c**14)
    r_lower = 1 + c**4 + c**12
    r_upper = r_lower + c**24 / (1 - c**16)
    return ThetaBounds(
        g_lower=sp.factor(2 * r_lower / theta_upper),
        g_upper=sp.factor(2 * r_upper / theta_lower),
        s_lower=sp.factor(1 / theta_upper**2),
        s_upper=sp.factor(1 / theta_lower**2),
    )


def negative_term_lower_bound(poly: sp.Poly, endpoint: sp.Rational) -> sp.Rational:
    """Lower-bound a polynomial by retaining its constant and negative terms."""

    result = sp.Integer(0)
    for (degree,), coefficient in poly.terms():
        if degree == 0:
            result += coefficient
        elif coefficient < 0:
            result += coefficient * endpoint**degree
    return sp.factor(result)


def certify_positive_rational(
    expression: sp.Expr,
    c: sp.Symbol,
    vanishing_order: int,
) -> tuple[sp.Rational, sp.Rational]:
    """Prove an edge residual positive by exact coefficient domination."""

    numerator, denominator = sp.cancel(expression / c**vanishing_order).as_numer_denom()
    numerator_poly = sp.Poly(numerator, c)
    denominator_poly = sp.Poly(denominator, c)
    sample = EDGE / 2
    if numerator_poly.eval(sample) < 0:
        numerator_poly = -numerator_poly
    if denominator_poly.eval(sample) < 0:
        denominator_poly = -denominator_poly
    numerator_lower = negative_term_lower_bound(numerator_poly, EDGE)
    denominator_lower = negative_term_lower_bound(denominator_poly, EDGE)
    if numerator_lower <= 0 or denominator_lower <= 0:
        raise AssertionError((numerator_lower, denominator_lower))
    return numerator_lower, denominator_lower


def certify_theta_remainders() -> list[tuple[sp.Rational, sp.Rational]]:
    """Prove the four remainder bounds used by each endpoint expansion."""

    c = sp.symbols("c", positive=True)
    bounds = theta_bounds(c)

    g_through_c6 = 2 - 4 * c**2 + 10 * c**4 - 20 * c**6
    s_through_c4 = 1 - 4 * c**2 + 12 * c**4
    g_through_c12 = g_through_c6 + 36 * c**8 - 64 * c**10 + 110 * c**12
    s_through_c12 = (
        s_through_c4 - 32 * c**6 + 76 * c**8 - 168 * c**10 + 352 * c**12
    )

    residuals = (
        (bounds.g_lower - (g_through_c6 - 37 * c**8), 8),
        (g_through_c6 + 37 * c**8 - bounds.g_upper, 8),
        (bounds.s_lower - (s_through_c4 - 33 * c**6), 6),
        (s_through_c4 + 33 * c**6 - bounds.s_upper, 6),
        (bounds.g_lower - (g_through_c12 - 200 * c**14), 14),
        (g_through_c12 + 200 * c**14 - bounds.g_upper, 14),
        (bounds.s_lower - (s_through_c12 - 750 * c**14), 14),
        (s_through_c12 + 750 * c**14 - bounds.s_upper, 14),
        (2 - bounds.g_upper, 2),
    )
    return [certify_positive_rational(value, c, order) for value, order in residuals]


def residual_numerator(endpoint: str) -> sp.Expr:
    """Build the numerator of r*F1*F2-g^2*D^2 at a cubic endpoint."""

    c, g, s, p = sp.symbols("c g s p")
    k = c * g**2
    d = k * (1 - p**2) / (1 - k**2 * p**2)
    lam = d * (k * p + 4 * c) / (2 * g)
    f1 = -4 * c**2 * d + 2 * c * g * (p + 1) - d * g**2 * p
    f2 = c**2 * d * g**2 * p - 2 * c * g * (p + 1) + 4 * d
    if endpoint == "upper":
        r = s * p + (1 - s) * p**3
        gap = (1 - s) * p * (1 - p**2)
    elif endpoint == "lower":
        a3 = s * (1 + k**2 - s**2) / 6
        r = s * p + a3 * p**3
        gap = p * (1 - s - a3 * p**2)
    else:
        raise ValueError(endpoint)
    small_d = gap + c**2 * (1 - p * r) - lam * (1 - r)
    residual = r * f1 * f2 - g**2 * small_d**2
    numerator, denominator = sp.together(residual).as_numer_denom()
    denominator_scale = 4 if endpoint == "upper" else 144
    expected_denominator = (
        denominator_scale
        * (c * g**2 * p - 1) ** 2
        * (c * g**2 * p + 1) ** 2
    )
    if sp.factor(denominator - expected_denominator) != 0:
        raise AssertionError(sp.factor(denominator))
    return numerator


def polynomial_ring() -> tuple[object, ...]:
    """Create the common exact ring used by both endpoint certificates."""

    return ring("c,g,s,p,x,G,S,A,B,C,D,E,F,M,N,P,Q,U,V", QQ)


def homogenized_residual(
    polynomial_ring_data: tuple[object, ...],
    endpoint: str,
) -> PolyElement:
    """Substitute p=p_star+c^4*x and clear its twelfth-power denominator."""

    polynomial_ring_object, c, g, s, p, x, *_ = polynomial_ring_data
    source = polynomial_ring_object.from_expr(residual_numerator(endpoint))
    d0 = 2 * g**2 * (1 - 2 * c**2 * g)
    p_numerator = g**2 - 4 * c**2 + c**4 * x * d0
    result = polynomial_ring_object.zero
    for monomial, coefficient in source.terms():
        c_power, g_power, s_power, p_power = monomial[:4]
        if any(monomial[4:]):
            raise AssertionError(monomial)
        result += (
            coefficient
            * c**c_power
            * g**g_power
            * s**s_power
            * p_numerator**p_power
            * d0 ** (12 - p_power)
        )
    return result


def compose(
    polynomial: PolyElement,
    variable: PolyElement,
    value: PolyElement,
    started: float,
) -> PolyElement:
    """Compose once and emit a compact progress record for this long audit."""

    result = polynomial.compose(variable, value)
    print(
        f"  {variable} substituted: {len(result.terms())} terms "
        f"({time.monotonic() - started:.1f}s)",
        flush=True,
    )
    return result


def minimum_c_power(polynomial: PolyElement) -> int:
    return min(monomial[0] for monomial, _ in polynomial.terms())


def divide_by_c_power(polynomial: PolyElement, power: int) -> PolyElement:
    """Divide a polynomial known to be termwise divisible by c**power."""

    result = polynomial.ring.zero
    for monomial, coefficient in polynomial.terms():
        if monomial[0] < power:
            raise AssertionError(monomial)
        result[(monomial[0] - power,) + monomial[1:]] = coefficient
    return result


def assert_active_variables(
    polynomial: PolyElement,
    active_indices: set[int],
) -> None:
    for monomial, _ in polynomial.terms():
        if any(monomial[index] for index in range(len(monomial)) if index not in active_indices):
            raise AssertionError(monomial)


def absolute_remainder_bound(
    polynomial: PolyElement,
    leading_power: int,
    first_remainder_power: int,
    limits: dict[int, object],
) -> object:
    """Bound a polynomial tail on the supplied symmetric box."""

    result = QQ(0)
    for monomial, coefficient in polynomial.terms():
        c_power = monomial[0]
        if c_power < first_remainder_power:
            continue
        term = abs(coefficient) * QQ(1, 20) ** (c_power - leading_power)
        for index, limit in limits.items():
            term *= QQ.convert(limit) ** monomial[index]
        result += term
    return result


def bernstein_rectangle_lower(
    expression: sp.Expr,
    first: sp.Symbol,
    first_interval: tuple[sp.Rational, sp.Rational],
    second: sp.Symbol,
    second_interval: tuple[sp.Rational, sp.Rational],
) -> sp.Rational:
    """Return the least exact tensor-product Bernstein coefficient."""

    u, v = sp.symbols("u v")
    first_lower, first_upper = first_interval
    second_lower, second_upper = second_interval
    transformed = sp.Poly(
        sp.expand(
            expression.subs(
                {
                    first: first_lower + (first_upper - first_lower) * u,
                    second: second_lower + (second_upper - second_lower) * v,
                }
            )
        ),
        u,
        v,
    )
    first_degree = transformed.degree(u)
    second_degree = transformed.degree(v)
    power_coefficients = {
        (i, j): transformed.coeff_monomial(u**i * v**j)
        for i in range(first_degree + 1)
        for j in range(second_degree + 1)
    }
    coefficients = []
    for k in range(first_degree + 1):
        for ell in range(second_degree + 1):
            coefficient = sp.Integer(0)
            for i in range(k + 1):
                for j in range(ell + 1):
                    coefficient += (
                        power_coefficients[i, j]
                        * sp.Rational(sp.binomial(k, i), sp.binomial(first_degree, i))
                        * sp.Rational(sp.binomial(ell, j), sp.binomial(second_degree, j))
                    )
            coefficients.append(coefficient)
    return min(coefficients)


def truncated_upper_lower_bounds(expression: sp.Expr) -> tuple[sp.Rational, ...]:
    """Bernstein-certify the upper leading polynomial on three x ranges."""

    t, x = sp.symbols("t x")
    ranges = ((-2, 2), (-4, -2), (2, 4))
    minima = []
    for lower, upper in ranges:
        box_minima = []
        left = sp.Rational(lower)
        while left < upper:
            right = min(left + sp.Rational(1, 8), sp.Rational(upper))
            box_minima.append(
                bernstein_rectangle_lower(
                    expression,
                    t,
                    (sp.Integer(0), EDGE**2),
                    x,
                    (left, right),
                )
            )
            left = right
        minima.append(min(box_minima))
    return tuple(minima)


def common_series_substitutions(
    data: tuple[object, ...],
    polynomial: PolyElement,
    started: float,
) -> PolyElement:
    _, c, g, s, _, _, big_g, big_s, a, b, c_rem, d_rem, *_ = data
    substitutions = (
        (g, 2 + c**2 * big_g),
        (s, 1 + c**2 * big_s),
        (big_g, -4 + c**2 * a),
        (big_s, -4 + c**2 * b),
        (a, 10 + c**2 * c_rem),
        (b, 12 + c**2 * d_rem),
    )
    for variable, value in substitutions:
        polynomial = compose(polynomial, variable, value, started)
    return polynomial


def certify_upper_endpoint(
    data: tuple[object, ...],
    started: float,
) -> tuple[tuple[float, ...], float, float]:
    """Generate and certify the upper cubic-envelope endpoint polynomial."""

    (
        _, c, _, _, _, x, _, _, _, _, c_rem, d_rem, e_rem, f_rem,
        m_rem, n_rem, p_rem, q_rem, u_rem, v_rem,
    ) = data
    polynomial = homogenized_residual(data, "upper")
    print(f"upper homogenized: {len(polynomial.terms())} terms", flush=True)
    polynomial = common_series_substitutions(data, polynomial, started)
    substitutions = (
        (c_rem, -20 + c**2 * e_rem),
        (d_rem, -32 + c**2 * f_rem),
        (e_rem, 36 + c**2 * m_rem),
        (f_rem, 76 + c**2 * n_rem),
    )
    for variable, value in substitutions:
        polynomial = compose(polynomial, variable, value, started)
    if minimum_c_power(polynomial) != 10:
        raise AssertionError("unexpected upper vanishing order")
    polynomial = divide_by_c_power(polynomial, 10)
    polynomial = compose(polynomial, m_rem, -64 + c**2 * p_rem, started)
    polynomial = compose(polynomial, n_rem, -168 + c**2 * q_rem, started)
    polynomial = compose(polynomial, p_rem, 110 + c**2 * u_rem, started)
    polynomial = compose(polynomial, q_rem, 352 + c**2 * v_rem, started)
    assert_active_variables(polynomial, {0, 4, 17, 18})

    scale = QQ(9_895_604_649_984)
    expected_truncated = (
        scale * (2 * x**2 + 4 * x + 3)
        + c**2
        * (
            -QQ(2_275_989_069_496_320) * x**2
            - QQ(4_690_516_604_092_416) * x
            - QQ(3_225_967_115_894_784)
        )
        + c**4
        * (
            -QQ(13_194_139_533_312) * x**3
            + QQ(132_911_164_588_818_432) * x**2
            + QQ(280_995_589_640_945_664) * x
            + QQ(177_393_556_757_938_176)
        )
    )
    actual_truncated = polynomial.ring.zero
    for monomial, coefficient in polynomial.terms():
        if monomial[0] <= 4:
            actual_truncated[monomial] = coefficient
    if actual_truncated != expected_truncated:
        raise AssertionError(actual_truncated)

    t_symbol, x_symbol = sp.symbols("t x")
    truncated_expression = sp.sympify(expected_truncated.as_expr()).subs(
        {sp.Symbol("c"): sp.sqrt(t_symbol), sp.Symbol("x"): x_symbol}
    )
    truncated_expression = sp.expand(truncated_expression)
    truncated_minima = truncated_upper_lower_bounds(truncated_expression)
    if not (
        truncated_minima[0] > sp.Rational(4, 5) * int(scale)
        and truncated_minima[1] > 2 * int(scale)
        and truncated_minima[2] > 14 * int(scale)
    ):
        raise AssertionError(truncated_minima)

    central = absolute_remainder_bound(
        polynomial, 0, 6, {4: 2, 17: 200, 18: 750}
    )
    outer = absolute_remainder_bound(
        polynomial, 0, 6, {4: 4, 17: 200, 18: 750}
    )
    if central >= QQ(1, 10) * scale or outer >= QQ(1, 4) * scale:
        raise AssertionError((central / scale, outer / scale))
    return (
        tuple(float(value / int(scale)) for value in truncated_minima),
        float(central / scale),
        float(outer / scale),
    )


def certify_lower_endpoint(data: tuple[object, ...], started: float) -> float:
    """Generate and certify the lower cubic-envelope endpoint polynomial."""

    _, c, _, _, _, _, _, _, _, _, c_rem, _, e_rem, *_ = data
    polynomial = homogenized_residual(data, "lower")
    print(f"lower homogenized: {len(polynomial.terms())} terms", flush=True)
    polynomial = common_series_substitutions(data, polynomial, started)
    polynomial = compose(polynomial, c_rem, -20 + c**2 * e_rem, started)
    if minimum_c_power(polynomial) != 8:
        raise AssertionError("unexpected lower vanishing order")
    assert_active_variables(polynomial, {0, 4, 10, 11})

    scale = QQ(316_659_348_799_488)
    lead_terms = [
        (monomial, coefficient)
        for monomial, coefficient in polynomial.terms()
        if monomial[0] == 8
    ]
    expected_monomial = (8,) + (0,) * (len(data) - 2)
    if lead_terms != [(expected_monomial, scale)]:
        raise AssertionError(lead_terms)
    remainder = absolute_remainder_bound(polynomial, 8, 10, {4: 4, 10: 33, 11: 37})
    if remainder >= QQ(1, 2) * scale:
        raise AssertionError(remainder / scale)
    return float(remainder / scale)


def main() -> None:
    started = time.monotonic()
    remainder_checks = certify_theta_remainders()
    print("theta remainder bounds: exact", flush=True)
    data = polynomial_ring()
    upper_minima, upper_central, upper_outer = certify_upper_endpoint(data, started)
    gc.collect()
    lower = certify_lower_endpoint(data, started)
    print("tiny-edge endpoint certificate: exact")
    print(f"  upper truncated minima/scale: {upper_minima}")
    print(f"  upper |x|<=2 tail/scale <= {upper_central:.12g}")
    print(f"  upper |x|<=4 tail/scale <= {upper_outer:.12g}")
    print(f"  lower |x|<=4 remainder/lead <= {lower:.12g}")
    print(f"  theta rational checks: {len(remainder_checks)}")
    print(f"  elapsed seconds: {time.monotonic() - started:.1f}")


if __name__ == "__main__":
    main()
