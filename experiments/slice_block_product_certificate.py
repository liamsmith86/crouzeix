#!/usr/bin/env python3
"""Exact Bernstein certificate for the elliptic-slice block-product theorem.

For the modal off-diagonal blocks ``B`` and ``C`` of ``phi(A)``, the
certificate proves ``||B|| ||C|| <= 2``.  The proof uses only integer
polynomial arithmetic after two elementary nome relaxations.  Every
polynomial is regenerated; there are no stored coefficient tables and no
floating-point decisions.
"""

from __future__ import annotations

import gc
import math
import time
from fractions import Fraction
from math import comb

import numpy as np
import sympy as sp


LOW_C_NUMERATOR = 37
LOW_C_DENOMINATOR = 125
HIGH_C_NUMERATOR = 59
HIGH_C_DENOMINATOR = 200


def audit_scalar_reduction() -> None:
    """Check the trace compression and the singular-value majorant exactly."""

    c, delta, modal_trace, w = sp.symbols(
        "c delta modal_trace w", positive=True
    )
    diagonal_squares = modal_trace - (1 + delta**2) * w
    trace_b = sp.expand(
        diagonal_squares + delta**2 * w / c**2 + c**2 * w
    )
    trace_c = sp.expand(
        diagonal_squares + w / c**2 + c**2 * delta**2 * w
    )
    expected_b = modal_trace - (1 - c**2) * (1 - delta**2 / c**2) * w
    expected_c = modal_trace + (1 - c**2) * (1 / c**2 - delta**2) * w
    assert sp.factor(trace_b - expected_b) == 0
    assert sp.factor(trace_c - expected_c) == 0

    large, small = sp.symbols("large small", nonnegative=True)
    upper_trace_bound = large + small - large * small / (large + small)
    assert sp.factor(upper_trace_bound - large - small**2 / (large + small)) == 0

    q = sp.symbols("q", positive=True)
    theta_lower = (1 - q**2) ** 2 / (1 + 2 * q - q**2) ** 2
    displacement_scale = 4 * (1 + q - q**2) / (1 + 2 * q - q**2) ** 2
    assert sp.factor((1 - theta_lower) / q - displacement_scale) == 0


def bernstein_value(
    coefficients: np.ndarray,
    coordinates: tuple[Fraction, Fraction, Fraction, Fraction],
) -> Fraction:
    """Evaluate a small integer Bernstein tensor exactly for machinery audits."""

    total = Fraction(0)
    for index, coefficient in np.ndenumerate(coefficients):
        basis = Fraction(1)
        for axis, bernstein_index in enumerate(index):
            degree = coefficients.shape[axis] - 1
            coordinate = coordinates[axis]
            basis *= (
                comb(degree, bernstein_index)
                * coordinate**bernstein_index
                * (1 - coordinate) ** (degree - bernstein_index)
            )
        total += int(coefficient) * basis
    return total


def audit_bernstein_machinery() -> None:
    """Check conversion, scaling, and rational subdivision on a test polynomial."""

    variables = sp.symbols("x0:4")
    x0, x1, x2, x3 = variables
    polynomial = sp.Poly(
        2 + 3 * x0 - 5 * x1 + 7 * x0**2 * x1**2 - 11 * x2 * x3,
        *variables,
        domain=sp.ZZ,
    )
    coefficients = power_to_bernstein(polynomial, (2, 3))
    degrees = polynomial.degree_list()
    denominator = 3 ** degrees[0]
    for degree in degrees:
        common = 1
        for index in range(degree + 1):
            common = math.lcm(common, comb(degree, index))
        denominator *= common

    point = (Fraction(1, 4), Fraction(2, 7), Fraction(1, 3), Fraction(3, 5))
    expected = Fraction(
        polynomial.as_expr().subs(
            dict(zip(variables, (Fraction(1, 6), *point[1:]), strict=True))
        )
    )
    assert bernstein_value(coefficients, point) / denominator == expected

    left, right = split_bernstein(coefficients, 1, 2, 5)
    child_denominator = denominator * 5 ** degrees[1]
    left_point = list(point)
    left_point[1] *= Fraction(2, 5)
    right_point = list(point)
    right_point[1] = Fraction(2, 5) + Fraction(3, 5) * right_point[1]
    assert bernstein_value(left, point) / child_denominator == Fraction(
        polynomial.as_expr().subs(
            dict(
                zip(
                    variables,
                    (Fraction(1, 6), *left_point[1:]),
                    strict=True,
                )
            )
        )
    )
    assert bernstein_value(right, point) / child_denominator == Fraction(
        polynomial.as_expr().subs(
            dict(
                zip(
                    variables,
                    (Fraction(1, 6), *right_point[1:]),
                    strict=True,
                )
            )
        )
    )


def build_base_polynomials() -> tuple[sp.Poly, sp.Poly]:
    """Build the regular low- and high-nome sufficient numerators."""

    c, p, eta, omega = sp.symbols("c p eta omega")
    q = c**2
    displacement_numerator = 4 * (1 + q - q**2)
    displacement_denominator = (1 + 2 * q - q**2) ** 2
    ratio_numerator = (
        displacement_denominator - q * displacement_numerator * eta
    )

    def poly(expression: sp.Expr | int) -> sp.Poly:
        return sp.Poly(expression, c, p, eta, omega, domain=sp.ZZ)

    one = poly(1)
    q_poly = poly(q)
    p_poly = poly(p)
    eta_poly = poly(eta)
    omega_poly = poly(omega)
    displacement_poly = poly(displacement_numerator)
    one_plus_r = poly(displacement_denominator + p * ratio_numerator)
    one_minus_pr = poly(displacement_denominator - p**2 * ratio_numerator)

    one_plus_r_squared = one_plus_r * one_plus_r
    one_plus_r_fourth = one_plus_r_squared * one_plus_r_squared
    one_minus_pr_squared = one_minus_pr * one_minus_pr
    displacement_squared = displacement_poly * displacement_poly

    x_numerator = (one + p_poly * p_poly) * one_plus_r_squared - (
        (one - q_poly)
        * (
            one_minus_pr_squared
            - p_poly
            * p_poly
            * q_poly
            * displacement_squared
            * eta_poly
            * eta_poly
        )
        * omega_poly
    )
    y_numerator = q_poly * (one + p_poly * p_poly) * one_plus_r_squared + (
        (one - q_poly)
        * (
            one_minus_pr_squared
            - p_poly
            * p_poly
            * q_poly**3
            * displacement_squared
            * eta_poly
            * eta_poly
        )
        * omega_poly
    )

    x_defect = x_numerator * x_numerator - p_poly * p_poly * one_plus_r_fourth
    y_defect = (
        y_numerator * y_numerator
        - p_poly * p_poly * q_poly * q_poly * one_plus_r_fourth
    )

    low_numerator = (
        16 * x_defect * y_defect
        - 4
        * (one + q_poly) ** 4
        * x_numerator
        * y_numerator
        * one_plus_r_fourth
    )
    high_numerator = (
        x_defect * y_defect
        - 4
        * q_poly
        * x_numerator
        * y_numerator
        * one_plus_r_fourth
    )
    return low_numerator, high_numerator


def centered_polynomial(
    polynomial: sp.Poly,
    chart: str,
    sign: int = 1,
) -> sp.Poly:
    """Pull the low-nome numerator into one of the three blow-up charts."""

    coefficients: dict[tuple[int, int, int, int], int] = {}
    for (c_degree, p_degree, eta_degree, omega_degree), coefficient in (
        polynomial.terms()
    ):
        for omega_scale_degree in range(omega_degree + 1):
            output_degrees = (
                range(omega_scale_degree + 1) if chart != "omega" else (0,)
            )
            for output_degree in output_degrees:
                if chart == "c":
                    monomial = (
                        c_degree + p_degree + omega_scale_degree,
                        p_degree,
                        eta_degree,
                        output_degree,
                    )
                    extra_factor = (
                        comb(omega_scale_degree, output_degree)
                        * 2**output_degree
                        * (-1) ** (omega_scale_degree - output_degree)
                    )
                elif chart == "p":
                    monomial = (
                        c_degree + p_degree + omega_scale_degree,
                        c_degree,
                        eta_degree,
                        output_degree,
                    )
                    extra_factor = (
                        comb(omega_scale_degree, output_degree)
                        * 2**output_degree
                        * (-1) ** (omega_scale_degree - output_degree)
                    )
                elif chart == "omega":
                    monomial = (
                        c_degree + p_degree + omega_scale_degree,
                        c_degree,
                        p_degree,
                        eta_degree,
                    )
                    extra_factor = sign**omega_scale_degree
                else:
                    raise ValueError(f"unknown centered chart: {chart}")

                # Multiplication by 16 clears the powers of 1/2 coming from
                # omega = 1/2 + scale.
                value = (
                    int(coefficient)
                    * comb(omega_degree, omega_scale_degree)
                    * 2 ** (4 - omega_degree + omega_scale_degree)
                    * extra_factor
                )
                coefficients[monomial] = coefficients.get(monomial, 0) + value

    variables = sp.symbols("x0:4")
    return sp.Poly.from_dict(coefficients, variables, domain=sp.ZZ)


def power_to_bernstein(
    polynomial: sp.Poly,
    first_coordinate_scale: tuple[int, int] = (1, 1),
) -> np.ndarray:
    """Return a common-denominator integer Bernstein tensor on [0,1]^4."""

    degrees = tuple(polynomial.degree_list())
    coefficients = np.zeros(tuple(degree + 1 for degree in degrees), dtype=object)
    scale_numerator, scale_denominator = first_coordinate_scale
    first_degree = degrees[0]
    scale_factors = [
        scale_numerator**degree * scale_denominator ** (first_degree - degree)
        for degree in range(first_degree + 1)
    ]
    for monomial, coefficient in polynomial.terms():
        coefficients[monomial] = int(coefficient) * scale_factors[monomial[0]]

    for axis in range(4):
        degree = coefficients.shape[axis] - 1
        power_coefficients = np.moveaxis(coefficients, axis, 0)
        bernstein_coefficients = np.empty_like(power_coefficients)
        common_denominator = 1
        for index in range(degree + 1):
            common_denominator = math.lcm(
                common_denominator,
                comb(degree, index),
            )
        weights = [
            common_denominator // comb(degree, index)
            for index in range(degree + 1)
        ]
        for bernstein_index in range(degree + 1):
            value = np.zeros_like(power_coefficients[0])
            for power_index in range(bernstein_index + 1):
                value += power_coefficients[power_index] * (
                    comb(bernstein_index, power_index) * weights[power_index]
                )
            bernstein_coefficients[bernstein_index] = value
        coefficients = np.moveaxis(bernstein_coefficients, 0, axis)
    return coefficients


def split_bernstein(
    coefficients: np.ndarray,
    axis: int,
    numerator: int,
    denominator: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Split an integer Bernstein tensor at an exact rational coordinate."""

    working = np.moveaxis(coefficients, axis, 0)
    degree = working.shape[0] - 1
    left = np.empty_like(working)
    right = np.empty_like(working)
    denominator_powers = [1]
    for _ in range(degree):
        denominator_powers.append(denominator_powers[-1] * denominator)
    left[0] = working[0] * denominator_powers[degree]
    right[degree] = working[degree] * denominator_powers[degree]
    for level in range(1, degree + 1):
        working = (
            (denominator - numerator) * working[:-1]
            + numerator * working[1:]
        )
        scale = denominator_powers[degree - level]
        left[level] = working[0] * scale
        right[degree - level] = working[-1] * scale
    return np.moveaxis(left, 0, axis), np.moveaxis(right, 0, axis)


def certify_nonpositive(label: str, coefficients: np.ndarray) -> None:
    maximum = max(coefficients.flat)
    if maximum > 0:
        raise AssertionError(f"{label} has a positive Bernstein coefficient")
    zero_count = sum(coefficient == 0 for coefficient in coefficients.flat)
    print(
        f"{label}: exact, {coefficients.size} coefficients, "
        f"{zero_count} zero, none positive"
    )


def certify_low_outer_boxes(low_numerator: sp.Poly) -> None:
    coefficients = power_to_bernstein(
        low_numerator,
        (LOW_C_NUMERATOR, LOW_C_DENOMINATOR),
    )
    low_p, high_p = split_bernstein(
        coefficients,
        1,
        LOW_C_NUMERATOR,
        LOW_C_DENOMINATOR,
    )
    high_p_left, high_p_right = split_bernstein(high_p, 1, 1, 2)
    low_omega, _ = split_bernstein(low_p, 3, 51, 250)
    _, high_omega = split_bernstein(low_p, 3, 199, 250)
    certify_nonpositive("low outer p-middle", high_p_left)
    certify_nonpositive("low outer p-high", high_p_right)
    certify_nonpositive("low outer omega-low", low_omega)
    certify_nonpositive("low outer omega-high", high_omega)


def certify_low_center(low_numerator: sp.Poly) -> None:
    charts = (
        ("center c-dominant", centered_polynomial(low_numerator, "c")),
        ("center p-dominant", centered_polynomial(low_numerator, "p")),
        (
            "center omega-positive",
            centered_polynomial(low_numerator, "omega", 1),
        ),
        (
            "center omega-negative",
            centered_polynomial(low_numerator, "omega", -1),
        ),
    )
    for label, polynomial in charts:
        print(f"{label}: regenerating {len(polynomial.terms())} power coefficients")
        coefficients = power_to_bernstein(
            polynomial,
            (LOW_C_NUMERATOR, LOW_C_DENOMINATOR),
        )
        certify_nonpositive(label, coefficients)
        del coefficients
        gc.collect()


def certify_high_range(high_numerator: sp.Poly) -> None:
    coefficients = power_to_bernstein(high_numerator)
    _, coefficients = split_bernstein(
        coefficients,
        0,
        HIGH_C_NUMERATOR,
        HIGH_C_DENOMINATOR,
    )
    boxes = [coefficients]
    for axis in range(4):
        children = []
        for box in boxes:
            children.extend(split_bernstein(box, axis, 1, 2))
        boxes = children
    for index, box in enumerate(boxes):
        certify_nonpositive(f"high box {index + 1:02d}", box)


def main() -> None:
    started = time.monotonic()
    audit_scalar_reduction()
    audit_bernstein_machinery()
    low_numerator, high_numerator = build_base_polynomials()
    print(
        "base numerators regenerated: "
        f"low={len(low_numerator.terms())}, high={len(high_numerator.terms())}"
    )
    certify_low_outer_boxes(low_numerator)
    certify_low_center(low_numerator)
    certify_high_range(high_numerator)
    print(f"block-product certificate: exact ({time.monotonic() - started:.1f}s)")


if __name__ == "__main__":
    main()
