#!/usr/bin/env python3
"""Rigorous certificate for the upper modal-block theorem.

The proof in ``proof/slice_upper_block_theorem.md`` reduces the only
transcendental step to three scalar inequalities in the ellipse nome.  This
script checks them with outward-rounded interval arithmetic.  It does not
sample the elliptic functions or solve for the modulus numerically.

Two elementary enclosures make the certificate finite.  With ``q = c**2``,

    theta_3(q) = 1 + 2 sum(q**(n*n)),
    ell = product(((1 + q**(2*n-1))/(1 + q**(2*n)))**2).

The omitted theta tail is bounded by a geometric series.  The logarithm of
the omitted ell product is at most another geometric series, and
``exp(x) <= 1/(1-x)`` supplies a rational upper endpoint.  All interval
operations after that are algebraic.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from typing import Iterable

from mpmath import iv
import sympy as sp


iv.dps = 50

SMALL_END = sp.Rational(1, 50)
INTERVAL_END = Decimal("0.65")
TAIL_TERMS = 12


@dataclass(frozen=True)
class IntervalResult:
    i1: object
    i2: object
    i3: object


def _interval_union(lower: object, upper: object) -> object:
    """Return the interval whose endpoints are the supplied point intervals."""

    return iv.mpf([lower, upper])


def _theta_and_ell_enclosures(c: object, terms: int) -> tuple[object, object]:
    q = c * c

    theta = iv.mpf(1)
    for n in range(1, terms + 1):
        theta += 2 * q ** (n * n)
    theta_tail = 2 * q ** ((terms + 1) ** 2) / (1 - q ** (2 * terms + 3))
    theta = _interval_union(theta.a, (theta + theta_tail).b)

    ell_product = iv.mpf(1)
    for n in range(1, terms + 1):
        numerator = 1 + q ** (2 * n - 1)
        denominator = 1 + q ** (2 * n)
        ell_product *= (numerator / denominator) ** 2

    log_tail = 2 * q ** (2 * terms + 1) / (1 - q * q)
    if not log_tail.b < 1:
        raise AssertionError("ell tail bound requires log_tail < 1")
    ell_upper = ell_product / (1 - log_tail)
    ell = _interval_union(ell_product.a, ell_upper.b)
    return theta, ell


def scalar_intervals(lower: str, upper: str) -> IntervalResult:
    """Enclose the three scalar residuals on a decimal interval for c."""

    c = iv.mpf([lower, upper])
    theta, ell = _theta_and_ell_enclosures(c, TAIL_TERMS)
    s = 1 / (theta * theta)
    k = 4 * c / (ell * ell)
    m = k * k

    k0 = (c + ell) * (1 - c * c) / (ell * (1 + c * ell))
    a3 = s * (1 + m - s * s) / 6
    r1 = (c + ell) * (1 - c * ell) / ((1 + c * ell) * (ell - c))

    return IntervalResult(
        i1=k0 - s,
        i2=k0 * (1 + 2 * c / ell) - s - 3 * a3,
        i3=s + a3 - r1,
    )


def _negative_term_lower_bound(poly: sp.Poly, endpoint: sp.Rational) -> sp.Rational:
    """Lower-bound a polynomial on [0, endpoint] by discarding positive terms."""

    coefficients = {monomial[0]: coefficient for monomial, coefficient in poly.terms()}
    lower = coefficients.get(0, sp.Integer(0))
    for degree, coefficient in coefficients.items():
        if degree and coefficient < 0:
            lower += coefficient * endpoint**degree
    return sp.factor(lower)


def _certify_positive_rational(
    expression: sp.Expr,
    variable: sp.Symbol,
    vanishing_order: int,
    endpoint: sp.Rational,
) -> tuple[sp.Rational, sp.Rational]:
    """Prove a rational expression positive by exact coefficient domination."""

    numerator, denominator = sp.cancel(expression).as_numer_denom()
    normalized_numerator = sp.Poly(sp.cancel(numerator / variable**vanishing_order), variable)
    denominator_poly = sp.Poly(denominator, variable)
    sample = endpoint / 2

    numerator_value = normalized_numerator.eval(sample)
    denominator_value = denominator_poly.eval(sample)
    if numerator_value * denominator_value <= 0:
        raise AssertionError("the rational residual is not positive at the sample point")

    if numerator_value < 0:
        normalized_numerator = -normalized_numerator
    if denominator_value < 0:
        denominator_poly = -denominator_poly

    numerator_lower = _negative_term_lower_bound(normalized_numerator, endpoint)
    denominator_lower = _negative_term_lower_bound(denominator_poly, endpoint)
    if not numerator_lower > 0 or not denominator_lower > 0:
        raise AssertionError((numerator_lower, denominator_lower))
    return numerator_lower, denominator_lower


def certify_small_range() -> list[tuple[sp.Rational, sp.Rational]]:
    """Certify 0 < c <= 1/50 using explicit rational majorants/minorants."""

    c = sp.symbols("c", positive=True)
    q = c * c

    ell_lower = ((1 + q) / (1 + q * q)) ** 2
    log_tail = 2 * q**3 / (1 - q * q)
    if not log_tail.subs(c, SMALL_END) < 1:
        raise AssertionError("small-range ell tail is not below one")
    ell_upper = sp.factor(ell_lower / (1 - log_tail))

    theta_upper = 1 + 2 * q + 2 * q**4 / (1 - q**5)
    s_lower = sp.factor(1 / theta_upper**2)
    s_upper = sp.factor(1 / (1 + 2 * q) ** 2)
    m_lower = sp.factor(16 * c * c / ell_upper**4)
    m_upper = sp.factor(16 * c * c / ell_lower**4)

    def k0(ell: sp.Expr) -> sp.Expr:
        return (c + ell) * (1 - c * c) / (ell * (1 + c * ell))

    def r1(ell: sp.Expr) -> sp.Expr:
        return (c + ell) * (1 - c * ell) / ((1 + c * ell) * (ell - c))

    # K0 decreases with ell.  The other substitutions use the monotonicities
    # recorded explicitly in the proof note.
    residuals = (
        k0(ell_upper) - s_upper,
        k0(ell_upper) * (1 + 2 * c / ell_upper)
        - s_upper * (3 + m_upper - s_upper**2) / 2,
        s_lower * (7 + m_lower - s_lower**2) / 6 - r1(ell_lower),
    )
    orders = (2, 1, 3)
    return [
        _certify_positive_rational(residual, c, order, SMALL_END)
        for residual, order in zip(residuals, orders, strict=True)
    ]


def _decimal_intervals() -> Iterable[tuple[str, str]]:
    ranges = (
        (Decimal("0.02"), Decimal("0.1"), Decimal("0.00001")),
        (Decimal("0.1"), INTERVAL_END, Decimal("0.0001")),
    )
    for start, stop, step in ranges:
        lower = start
        while lower < stop:
            upper = min(lower + step, stop)
            yield str(lower), str(upper)
            lower = upper


def certify_interval_range() -> tuple[int, list[tuple[float, tuple[str, str]]]]:
    """Certify 1/50 <= c <= 13/20 by outward-rounded intervals."""

    minima = [(float("inf"), ("", "")) for _ in range(3)]
    count = 0
    for lower, upper in _decimal_intervals():
        result = scalar_intervals(lower, upper)
        for index, residual in enumerate((result.i1, result.i2, result.i3)):
            if not residual.a > 0:
                raise AssertionError((index + 1, lower, upper, residual))
            reported_lower = float(residual.a)
            if reported_lower < minima[index][0]:
                minima[index] = (reported_lower, (lower, upper))
        count += 1
    return count, minima


def certify_trivial_range() -> sp.Rational:
    """Prove c*ell >= 1 for 13/20 <= c < 1 from two product factors."""

    c = sp.symbols("c", positive=True)
    q = c * c
    product = sp.prod(
        ((1 + q ** (2 * n - 1)) / (1 + q ** (2 * n))) ** 2
        for n in range(1, 3)
    )
    numerator = sp.cancel(c * product - 1).as_numer_denom()[0]
    quotient = sp.Poly(sp.cancel(numerator / (1 - c)), c)
    derivative = sp.Poly(sp.diff(quotient.as_expr(), c), c)
    if any(coefficient < 0 for coefficient in derivative.all_coeffs()):
        raise AssertionError("the two-factor quotient is not increasing")
    value = sp.factor(quotient.eval(sp.Rational(13, 20)))
    if not value > 0:
        raise AssertionError(value)
    return value


def audit_determinant_factorization() -> None:
    """Symbolically check the two roots used in the modal determinant proof."""

    c, ell, r, p = sp.symbols("c ell r p", positive=True)
    h = 4 / ell**2
    matrix = sp.sqrt(h) / (1 + r) * sp.Matrix(
        [
            [sp.sqrt(r) * (1 + p), -(p - r) / c],
            [c * (1 - p * r), sp.sqrt(r) * (1 + p)],
        ]
    )
    determinant = sp.factor((4 * sp.eye(2) - matrix.T * matrix).det())
    p_minus = ell * (c**2 - c * ell * (1 + r) + r) / (
        c**2 * ell * r - c * (1 + r) + ell
    )
    p_plus = ell * (c**2 + c * ell * (1 + r) + r) / (
        c**2 * ell * r + c * (1 + r) + ell
    )
    positive_factor = (
        16
        * (c**2 * ell * r - c * (1 + r) + ell)
        * (c**2 * ell * r + c * (1 + r) + ell)
        / (c**2 * ell**4 * (1 + r) ** 2)
    )
    if sp.factor(determinant - positive_factor * (p - p_minus) * (p_plus - p)) != 0:
        raise AssertionError("modal determinant factorization failed")


def main() -> None:
    audit_determinant_factorization()
    small_bounds = certify_small_range()
    count, interval_minima = certify_interval_range()
    trivial_value = certify_trivial_range()

    print("determinant factorization: exact")
    print("small-range exact coefficient bounds:")
    for index, (numerator, denominator) in enumerate(small_bounds, start=1):
        print(f"  I{index}: numerator >= {float(numerator):.12g}, "
              f"denominator magnitude >= {float(denominator):.12g}")
    print(f"interval boxes certified: {count}")
    for index, (lower, box) in enumerate(interval_minima, start=1):
        print(f"  I{index}: lower endpoint >= {lower:.12g} on {box}")
    print(f"trivial-range product quotient at c=13/20: {float(trivial_value):.12g}")


if __name__ == "__main__":
    main()
