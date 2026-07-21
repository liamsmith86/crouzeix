#!/usr/bin/env python3
"""Exact certificate for the full tiny-nome odd-block range.

Together with ``slice_odd_tiny_edge_certificate.py``, this checker proves
the two L32 endpoint residuals positive on their discriminant branch for

    0 < c <= 1/20,  0 <= p <= 1.

The proof partitions the p interval into two boundary regions, a regular
middle region, and the singular central band.  This file certifies the
boundary and regular-middle estimates and proves that the branch cannot
enter the two excluded boundary strips.  The companion checker certifies
the central band in coordinates centered at the exact vertex p_star.

All arithmetic is rational.  Every polynomial is regenerated from the
compact residual and branch formulas; no saved coefficient table is used.
"""

from __future__ import annotations

from collections import defaultdict
from math import comb
import time

import sympy as sp
from sympy import QQ
from sympy.polys.rings import PolyElement, ring

from slice_odd_tiny_edge_certificate import (
    EDGE,
    certify_positive_rational,
    residual_numerator,
    theta_bounds,
)


W_BOUND = QQ(300)
Z_BOUND = QQ(1500)
LOWER_P = QQ(1, 4)
UPPER_P = QQ(2, 3)
BRANCH_TOP = QQ(3, 4)
MIDDLE_SEPARATION = QQ(96)


def polynomial_ring() -> tuple[object, ...]:
    """Return the exact ring used throughout this checker."""

    return ring("c,g,s,p,q,z,W,Z", QQ)


def deep_theta_series(c: PolyElement, w: PolyElement, z: PolyElement) -> tuple[PolyElement, PolyElement]:
    """Theta quotients through c^14, with bounded c^16 remainders."""

    g_series = (
        2
        - 4 * c**2
        + 10 * c**4
        - 20 * c**6
        + 36 * c**8
        - 64 * c**10
        + 110 * c**12
        - 180 * c**14
        + c**16 * w
    )
    s_series = (
        1
        - 4 * c**2
        + 12 * c**4
        - 32 * c**6
        + 76 * c**8
        - 168 * c**10
        + 352 * c**12
        - 704 * c**14
        + c**16 * z
    )
    return g_series, s_series


def certify_deep_theta_bounds() -> list[tuple[sp.Rational, sp.Rational]]:
    """Prove the remainder and coarse bounds used by coefficient domination."""

    c = sp.symbols("c", positive=True)
    bounds = theta_bounds(c)
    g_polynomial = (
        2
        - 4 * c**2
        + 10 * c**4
        - 20 * c**6
        + 36 * c**8
        - 64 * c**10
        + 110 * c**12
        - 180 * c**14
    )
    s_polynomial = (
        1
        - 4 * c**2
        + 12 * c**4
        - 32 * c**6
        + 76 * c**8
        - 168 * c**10
        + 352 * c**12
        - 704 * c**14
    )
    residuals = (
        (bounds.g_lower - (g_polynomial - 300 * c**16), 16),
        (g_polynomial + 300 * c**16 - bounds.g_upper, 16),
        (bounds.s_lower - (s_polynomial - 1500 * c**16), 16),
        (s_polynomial + 1500 * c**16 - bounds.s_upper, 16),
        (2 - bounds.g_upper, 2),
        (bounds.g_lower - sp.Rational(19, 10), 0),
        (1 - bounds.s_upper, 2),
    )
    results = [certify_positive_rational(value, c, order) for value, order in residuals]

    # These coarse estimates give 1/2 <= p_star <= 1/2 + 2c^2:
    # (p_star-1/2)/c^2=(g^3-2)/(g^2(1-2c^2g)).
    coefficient_upper = sp.Rational(6) / (
        sp.Rational(19, 10) ** 2 * sp.Rational(99, 100)
    )
    if coefficient_upper >= 2 or sp.Rational(19, 10) ** 3 <= 2:
        raise AssertionError(coefficient_upper)
    return results


def series_residual(data: tuple[object, ...], endpoint: str) -> PolyElement:
    """Substitute the certified theta series into one endpoint numerator."""

    polynomial_ring_object, c, g, s, _, _, _, w, z = data
    source = polynomial_ring_object.from_expr(residual_numerator(endpoint))
    g_series, s_series = deep_theta_series(c, w, z)
    result = source.compose(g, g_series).compose(s, s_series)
    if any(monomial[index] for monomial, _ in result.terms() for index in (1, 2, 4, 5)):
        raise AssertionError("inactive variable survived residual substitution")
    return result


def branch_numerator(endpoint: str) -> tuple[sp.Expr, sp.Expr]:
    """Return the numerator and factored denominator of 2A-B."""

    c, g, s, p = sp.symbols("c g s p")
    k = c * g**2
    d = k * (1 - p**2) / (1 - k**2 * p**2)
    a3 = s * (1 + k**2 - s**2) / 6
    if endpoint == "upper":
        r = s * p + (1 - s) * p**3
    elif endpoint == "lower":
        r = s * p + a3 * p**3
    else:
        raise ValueError(endpoint)
    a_form = g * p * (2 * d * (1 + c**2 * r) - c * g * (1 + r)) / (c * (1 + r))
    b_form = (
        2 * c**2 * g * p * r
        - 2 * c**2 * g
        - c * d * g**2 * p * (1 + r)
        + 4 * c * d * (1 + r)
        + 2 * g * p
        - 2 * g * r
    ) / (c * (1 + r))
    numerator, denominator = sp.together(2 * a_form - b_form).as_numer_denom()
    return numerator, sp.factor(denominator)


def series_branch(data: tuple[object, ...], endpoint: str) -> tuple[PolyElement, sp.Expr]:
    """Substitute the theta series into the branch numerator."""

    polynomial_ring_object, c, g, s, _, _, _, w, z = data
    numerator, denominator = branch_numerator(endpoint)
    source = polynomial_ring_object.from_expr(numerator)
    g_series, s_series = deep_theta_series(c, w, z)
    return source.compose(g, g_series).compose(s, s_series), denominator


def divide_by_c_power(polynomial: PolyElement, power: int) -> PolyElement:
    """Divide a polynomial known termwise divisible by c**power."""

    result = polynomial.ring.zero
    for monomial, coefficient in polynomial.terms():
        if monomial[0] < power:
            raise AssertionError(monomial)
        result[(monomial[0] - power,) + monomial[1:]] = coefficient
    return result


def bounded_coefficient(monomial: tuple[int, ...], coefficient: object) -> object:
    """Absorb the two certified theta-remainder variables."""

    return abs(coefficient) * W_BOUND ** monomial[6] * Z_BOUND ** monomial[7]


def certify_lower_residual_region(
    data: tuple[object, ...],
    endpoint: str,
    polynomial: PolyElement,
) -> sp.Rational:
    """Prove positivity for c/4 <= p <= 1/4."""

    polynomial_ring_object, c, _, _, p, *_ = data
    scale = QQ(1) if endpoint == "upper" else QQ(36)
    leading = 64 * scale * p * (p + 1) ** 2 * (2 * p - 1) ** 2 * (3 - 4 * p)
    actual_leading = polynomial_ring_object.zero
    for monomial, coefficient in polynomial.terms():
        if monomial[0] == 2:
            actual_leading[monomial] = coefficient
    if actual_leading != c**2 * leading:
        raise AssertionError((endpoint, actual_leading.as_expr()))

    # On [0,1/4], the leading coefficient of p is at least 50.
    correction = QQ(0)
    for monomial, coefficient in polynomial.terms():
        c_power, p_power = monomial[0], monomial[3]
        if c_power == 2:
            continue
        if c_power < 4 or c_power % 2:
            raise AssertionError(monomial)
        term = bounded_coefficient(monomial, coefficient)
        if p_power:
            term *= EDGE ** (c_power - 2) * LOWER_P ** (p_power - 1)
        else:
            # p >= c/4 turns a p-independent term into a multiple of p.
            term *= 4 * EDGE ** (c_power - 3)
        correction += term
    lower_coefficient = 50 * scale
    if correction >= lower_coefficient:
        raise AssertionError((endpoint, correction / lower_coefficient))
    return correction / lower_coefficient


def certify_lower_branch_exclusion(
    data: tuple[object, ...],
    endpoint: str,
) -> sp.Rational:
    """Prove 2A-B<0 when 0 <= p < c/4."""

    _, c, _, _, p, _, z_coordinate, *_ = data
    polynomial, denominator = series_branch(data, endpoint)
    polynomial = polynomial.compose(p, c * z_coordinate)
    power = min(monomial[0] for monomial, _ in polynomial.terms())
    if power != 2:
        raise AssertionError((endpoint, power))
    polynomial = divide_by_c_power(polynomial, power)

    lead_scale = QQ(1) if endpoint == "upper" else QQ(6)
    expected_lead = lead_scale * (24 * z_coordinate - 12)
    actual_lead = polynomial.ring.zero
    for monomial, coefficient in polynomial.terms():
        if monomial[0] == 0:
            actual_lead[monomial] = coefficient
    if actual_lead != expected_lead:
        raise AssertionError((endpoint, actual_lead.as_expr()))

    remainder = QQ(0)
    for monomial, coefficient in polynomial.terms():
        if monomial[0] == 0:
            continue
        remainder += (
            bounded_coefficient(monomial, coefficient)
            * EDGE ** monomial[0]
            * LOWER_P ** monomial[5]
        )
    margin = 6 if endpoint == "upper" else 36
    if remainder >= QQ(margin, 12):
        raise AssertionError((endpoint, remainder, denominator))
    return remainder / margin


def certify_upper_residual_region(
    data: tuple[object, ...],
    endpoint: str,
    polynomial: PolyElement,
) -> tuple[sp.Rational, sp.Rational]:
    """Prove positivity for 2/3 <= p <= 3/4."""

    polynomial_ring_object, c, _, _, p, q_coordinate, *_ = data
    polynomial = divide_by_c_power(polynomial.compose(p, BRANCH_TOP - q_coordinate), 2)
    scale = QQ(1) if endpoint == "upper" else QQ(36)
    leading_q = QQ(12800, 243) * scale

    actual_leading = polynomial_ring_object.zero
    for monomial, coefficient in polynomial.terms():
        if monomial[0] == 0:
            actual_leading[monomial] = coefficient
    source_p = BRANCH_TOP - q_coordinate
    expected_leading = (
        64
        * scale
        * source_p
        * (source_p + 1) ** 2
        * (2 * source_p - 1) ** 2
        * (3 - 4 * source_p)
    )
    if actual_leading != expected_leading:
        raise AssertionError((endpoint, actual_leading.as_expr()))

    base_c = QQ(0)
    q_correction = QQ(0)
    c_correction = QQ(0)
    for monomial, coefficient in polynomial.terms():
        c_power, q_power = monomial[0], monomial[4]
        if c_power == 0:
            continue
        if c_power == 2 and q_power == 0 and not monomial[6] and not monomial[7]:
            base_c += coefficient
            continue
        term = bounded_coefficient(monomial, coefficient)
        if q_power:
            q_correction += term * EDGE**c_power * QQ(1, 12) ** (q_power - 1)
        else:
            if c_power < 2:
                raise AssertionError(monomial)
            c_correction += term * EDGE ** (c_power - 2)

    expected_base = QQ(68943, 256) * scale
    if base_c != expected_base:
        raise AssertionError((endpoint, base_c))
    if q_correction >= QQ(3, 4) * leading_q or c_correction >= QQ(1, 10) * base_c:
        raise AssertionError((endpoint, q_correction / leading_q, c_correction / base_c))
    return q_correction / leading_q, c_correction / base_c


def certify_upper_branch_exclusion(
    data: tuple[object, ...],
    endpoint: str,
) -> tuple[sp.Rational, sp.Rational]:
    """Prove 2A-B<0 when 3/4 <= p <= 1."""

    polynomial_ring_object, c, _, _, p, q_coordinate, *_ = data
    polynomial, _ = series_branch(data, endpoint)
    polynomial = polynomial.compose(p, BRANCH_TOP + q_coordinate)
    power = min(monomial[0] for monomial, _ in polynomial.terms())
    if power != 1:
        raise AssertionError((endpoint, power))
    polynomial = divide_by_c_power(polynomial, power)

    q_margin = QQ(42) if endpoint == "upper" else QQ(252)
    c_margin = QQ(105, 16) if endpoint == "upper" else QQ(315, 8)
    branch_scale = QQ(1) if endpoint == "upper" else QQ(6)
    actual_c0 = polynomial_ring_object.zero
    c1_constant = QQ(0)
    for monomial, coefficient in polynomial.terms():
        if monomial[0] == 0:
            actual_c0[monomial] = coefficient
        elif monomial[0] == 1 and monomial[4] == 0 and not monomial[6] and not monomial[7]:
            c1_constant += coefficient
    expected_c0 = -branch_scale * (32 * q_coordinate**3 + 80 * q_coordinate**2 + 42 * q_coordinate)
    if actual_c0 != expected_c0 or c1_constant != -c_margin:
        raise AssertionError((endpoint, actual_c0.as_expr(), c1_constant))
    q_positive = QQ(0)
    c_positive = QQ(0)
    for monomial, coefficient in polynomial.terms():
        has_uncertain_remainder = bool(monomial[6] or monomial[7])
        if coefficient <= 0 and not has_uncertain_remainder:
            continue
        c_power, q_power = monomial[0], monomial[4]
        term = abs(coefficient) * W_BOUND ** monomial[6] * Z_BOUND ** monomial[7]
        if q_power:
            q_positive += term * EDGE**c_power * LOWER_P ** (q_power - 1)
        else:
            if c_power < 1:
                raise AssertionError(monomial)
            c_positive += term * EDGE ** (c_power - 1)
    if q_positive >= QQ(9, 100) * q_margin or c_positive >= QQ(7, 10) * c_margin:
        raise AssertionError((endpoint, q_positive / q_margin, c_positive / c_margin))
    return q_positive / q_margin, c_positive / c_margin


def bernstein_absolute_bound(
    coefficients: dict[int, object],
    lower: sp.Rational,
    upper: sp.Rational,
) -> sp.Rational:
    """Bound a univariate polynomial by its exact Bernstein coefficients."""

    degree = max(coefficients, default=0)
    power_coefficients = [QQ(0)] * (degree + 1)
    for source_degree, coefficient in coefficients.items():
        for target_degree in range(source_degree + 1):
            power_coefficients[target_degree] += (
                coefficient
                * comb(source_degree, target_degree)
                * lower ** (source_degree - target_degree)
                * (upper - lower) ** target_degree
            )
    bernstein_coefficients = []
    for degree_index in range(degree + 1):
        bernstein_coefficients.append(
            sum(
                power_coefficients[index]
                * QQ(comb(degree_index, index), comb(degree, index))
                for index in range(degree_index + 1)
            )
        )
    return max(abs(value) for value in bernstein_coefficients)


def middle_polynomial_bound(coefficients: dict[int, object]) -> sp.Rational:
    """Bound one p-polynomial on [1/4,2/3] by exact subdivisions."""

    result = QQ(0)
    left = LOWER_P
    while left < UPPER_P:
        right = min(left + QQ(1, 96), UPPER_P)
        result = max(result, bernstein_absolute_bound(coefficients, left, right))
        left = right
    return result


def certify_middle_region(
    data: tuple[object, ...],
    endpoint: str,
    polynomial: PolyElement,
) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
    """Prove positivity when |2p-1| >= 96c^2 in the middle interval."""

    groups: defaultdict[tuple[int, int, int], dict[int, object]] = defaultdict(dict)
    for monomial, coefficient in polynomial.terms():
        key = monomial[0], monomial[6], monomial[7]
        p_power = monomial[3]
        groups[key][p_power] = groups[key].get(p_power, QQ(0)) + coefficient

    bounds = {key: middle_polynomial_bound(coefficients) for key, coefficients in groups.items()}
    p_symbol = sp.symbols("p")
    c4_expression = sum(
        sp.Rational(coefficient.numerator, coefficient.denominator) * p_symbol**p_power
        for (c_power, w_power, z_power), coefficients in groups.items()
        if c_power == 4 and not w_power and not z_power
        for p_power, coefficient in coefficients.items()
    )
    quotient, remainder = sp.div(sp.Poly(c4_expression, p_symbol), sp.Poly(2 * p_symbol - 1, p_symbol))
    if remainder.as_expr() != 0:
        raise AssertionError((endpoint, remainder))
    first_bound = middle_polynomial_bound(
        {degree[0]: QQ.convert(coefficient) for degree, coefficient in quotient.terms()}
    )

    higher_bound = QQ(0)
    for (c_power, w_power, z_power), bound in bounds.items():
        if c_power < 6:
            continue
        higher_bound += (
            bound
            * W_BOUND**w_power
            * Z_BOUND**z_power
            * EDGE ** (c_power - 6)
        )

    scale = QQ(1) if endpoint == "upper" else QQ(36)
    leading_coefficient = QQ(25, 3) * scale
    ratio = (
        first_bound / (leading_coefficient * MIDDLE_SEPARATION)
        + higher_bound / (leading_coefficient * MIDDLE_SEPARATION**2)
    )
    if ratio >= QQ(7, 10):
        raise AssertionError((endpoint, first_bound, higher_bound, ratio))
    return first_bound, higher_bound, ratio


def certify_branch_denominators() -> None:
    """Audit the exact denominator factorizations used in both exclusions."""

    c, g, s, p = sp.symbols("c g s p")
    expected_upper = (
        c
        * (p + 1)
        * (c * g**2 * p - 1)
        * (c * g**2 * p + 1)
        * (p**2 * s - p**2 - p * s + p - 1)
    )
    expected_lower = (
        c
        * (c * g**2 * p - 1)
        * (c * g**2 * p + 1)
        * (-c**2 * g**4 * p**3 * s + p**3 * s**3 - p**3 * s - 6 * p * s - 6)
    )
    for endpoint, expected in (("upper", expected_upper), ("lower", expected_lower)):
        _, actual = branch_numerator(endpoint)
        if sp.factor(actual - expected) != 0:
            raise AssertionError((endpoint, actual))


def certify_elementary_leading_bounds() -> None:
    """Audit the three elementary lower bounds used in the partition."""

    p = sp.symbols("p")
    lower_coefficient = 64 * (p + 1) ** 2 * (1 - 2 * p) ** 2 * (3 - 4 * p)
    expected_derivative = -128 * (p + 1) * (2 * p - 1) * (20 * p**2 - 6 * p - 5)
    if sp.factor(sp.diff(lower_coefficient, p) - expected_derivative) != 0:
        raise AssertionError(sp.factor(sp.diff(lower_coefficient, p)))
    if lower_coefficient.subs(p, LOWER_P) != 50:
        raise AssertionError(lower_coefficient.subs(p, LOWER_P))

    # On [1/4,2/3], p>=1/4, (p+1)^2>=25/16, and 3-4p>=1/3.
    if 64 * LOWER_P * QQ(25, 16) * QQ(1, 3) != QQ(25, 3):
        raise AssertionError("middle leading bound")
    # On [2/3,3/4], the same factors and (3-4p)=4q give this q coefficient.
    if 64 * UPPER_P * QQ(25, 9) * QQ(1, 9) * 4 != QQ(12800, 243):
        raise AssertionError("upper leading bound")


def main() -> None:
    started = time.monotonic()
    theta_checks = certify_deep_theta_bounds()
    certify_branch_denominators()
    certify_elementary_leading_bounds()
    print("deep theta and branch-denominator bounds: exact", flush=True)

    data = polynomial_ring()
    for endpoint in ("upper", "lower"):
        polynomial = series_residual(data, endpoint)
        print(f"{endpoint} residual series: {len(polynomial.terms())} terms", flush=True)
        lower_ratio = certify_lower_residual_region(data, endpoint, polynomial)
        upper_ratios = certify_upper_residual_region(data, endpoint, polynomial)
        middle = certify_middle_region(data, endpoint, polynomial)
        lower_branch = certify_lower_branch_exclusion(data, endpoint)
        upper_branch = certify_upper_branch_exclusion(data, endpoint)
        print(
            f"  {endpoint}: lower correction={float(lower_ratio):.12g}, "
            f"upper corrections=({float(upper_ratios[0]):.12g},"
            f" {float(upper_ratios[1]):.12g})",
            flush=True,
        )
        print(
            f"  {endpoint}: middle ratio={float(middle[2]):.12g}, "
            f"branch ratios=({float(lower_branch):.12g},"
            f" {float(upper_branch[0]):.12g}, {float(upper_branch[1]):.12g})",
            flush=True,
        )

    print("small-edge regular-region certificate: exact")
    print(f"  theta rational checks: {len(theta_checks)}")
    print("  central remainder: |p-p_star|/c^2 < 50 (companion checker)")
    print(f"  elapsed seconds: {time.monotonic() - started:.1f}")


if __name__ == "__main__":
    main()
