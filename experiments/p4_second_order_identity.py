#!/usr/bin/env python3
"""Regenerate the exact p=4 second-order L21 identity.

All 32 real coordinates of a general complex perturbation ``E`` are kept
symbolic.  The script reconstructs the physical conformal coefficients and
independently eliminates the free first metric variation before checking the
claimed negative sum of squares as a zero polynomial.
"""

from __future__ import annotations

import sympy as sp

from crabb_second_order_symbolic import (
    conformal_coefficients,
    symbolic_perturbation,
)


def main() -> None:
    z = sp.symbols("z", nonzero=True)
    perturbation, real_parts, imaginary_parts = symbolic_perturbation(4)
    root_two = sp.sqrt(2)
    root_three = sp.sqrt(3)
    root_six = sp.sqrt(6)
    base = sp.Matrix(
        [
            [0, root_two, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, root_two],
            [0, 0, 0, 0],
        ]
    )
    metric = sp.diag(1, 2, 2, 4)

    top_vector = sp.Matrix(
        [1 / root_six, z / root_three, z**2 / root_three, z**3 / root_six]
    )
    half_vector = sp.Matrix(
        [-1 / root_three, -z / root_six, z**2 / root_six, z**3 / root_three]
    )
    minus_half_vector = sp.Matrix(
        [1 / root_three, -z / root_six, -(z**2) / root_six, z**3 / root_three]
    )
    bottom_vector = sp.Matrix(
        [-1 / root_six, z / root_three, -(z**2) / root_three, z**3 / root_six]
    )
    tangent, second_order = conformal_coefficients(
        base,
        perturbation,
        top_vector,
        [
            (sp.Rational(1, 2), half_vector),
            (-sp.Rational(1, 2), minus_half_vector),
            (-sp.Integer(1), bottom_vector),
        ],
        z,
    )

    first_stein_forcing = sp.simplify(
        tangent.conjugate().T * metric * base + base.conjugate().T * metric * tangent
    )
    x_real, x_imag, y_real, y_imag, z_real, z_imag = sp.symbols(
        "xr xi yr yi zr zi",
        real=True,
    )
    x = x_real + sp.I * x_imag
    y = y_real + sp.I * y_imag
    free_z = z_real + sp.I * z_imag
    metric_tangent = sp.zeros(4)
    metric_tangent[0, 1] = x
    metric_tangent[1, 0] = sp.conjugate(x)
    metric_tangent[0, 2] = y
    metric_tangent[2, 0] = sp.conjugate(y)
    metric_tangent[0, 3] = free_z
    metric_tangent[3, 0] = sp.conjugate(free_z)
    metric_tangent[1, 1] = first_stein_forcing[1, 1]
    metric_tangent[2, 2] = first_stein_forcing[1, 1] + first_stein_forcing[2, 2]
    metric_tangent[1, 2] = root_two * x + first_stein_forcing[1, 2]
    metric_tangent[2, 1] = sp.conjugate(metric_tangent[1, 2])
    metric_tangent[1, 3] = 2 * y + first_stein_forcing[1, 3]
    metric_tangent[3, 1] = sp.conjugate(metric_tangent[1, 3])
    metric_tangent[2, 3] = (
        2 * x + root_two * first_stein_forcing[1, 2] + first_stein_forcing[2, 3]
    )
    metric_tangent[3, 2] = sp.conjugate(metric_tangent[2, 3])

    second_stein_forcing = sp.simplify(
        second_order.conjugate().T * metric * base
        + base.conjugate().T * metric * second_order
        + tangent.conjugate().T * metric * tangent
        + tangent.conjugate().T * metric_tangent * base
        + base.conjugate().T * metric_tangent * tangent
    )
    first_stein_variation = sp.simplify(
        metric_tangent
        - base.conjugate().T * metric_tangent * base
        - first_stein_forcing
    )
    kernel_to_range = first_stein_variation[1:4, 0]
    contraction_penalty = sp.expand(
        2
        * (
            second_stein_forcing[1, 1]
            + kernel_to_range[0] * sp.conjugate(kernel_to_range[0])
        )
        + 2
        * (
            second_stein_forcing[2, 2]
            + kernel_to_range[1] * sp.conjugate(kernel_to_range[1])
        )
        + second_stein_forcing[3, 3]
        + kernel_to_range[2] * sp.conjugate(kernel_to_range[2])
    )
    lower_penalty = sp.expand(
        x * sp.conjugate(x) + y * sp.conjugate(y) + free_z * sp.conjugate(free_z) / 3
    )
    upper_penalty = sp.expand(
        metric_tangent[3, 0] * metric_tangent[0, 3] / 3
        + metric_tangent[3, 1] * metric_tangent[1, 3] / 2
        + metric_tangent[3, 2] * metric_tangent[2, 3] / 2
    )
    objective = sp.expand(4 * lower_penalty + contraction_penalty + upper_penalty)

    free_variables = (x_real, x_imag, y_real, y_imag, z_real, z_imag)
    zero_free = {variable: 0 for variable in free_variables}
    constant = objective.subs(zero_free)
    linear_coefficients = [
        sp.diff(objective, variable).subs(zero_free) for variable in free_variables
    ]
    quadratic_remainder = sp.expand(
        objective
        - constant
        - sum(
            coefficient * variable
            for coefficient, variable in zip(linear_coefficients, free_variables)
        )
    )
    expected_quadratic = sp.expand(
        8 * (x_real**2 + x_imag**2 + y_real**2 + y_imag**2)
        + sp.Rational(8, 3) * (z_real**2 + z_imag**2)
    )
    if sp.simplify(quadratic_remainder - expected_quadratic) != 0:
        raise AssertionError("unexpected p=4 metric-variable Hessian")
    optimum = sp.expand(
        constant
        - sum(coefficient**2 for coefficient in linear_coefficients[:4]) / 32
        - 3 * sum(coefficient**2 for coefficient in linear_coefficients[4:]) / 32
    )

    def entry(row: int, column: int) -> tuple[sp.Symbol, sp.Symbol]:
        index = 4 * row + column
        return real_parts[index], imaginary_parts[index]

    e00r, e00i = entry(0, 0)
    e01r, _ = entry(0, 1)
    e02r, e02i = entry(0, 2)
    e11r, e11i = entry(1, 1)
    e12r, _ = entry(1, 2)
    e13r, e13i = entry(1, 3)
    e20r, e20i = entry(2, 0)
    e22r, e22i = entry(2, 2)
    e23r, _ = entry(2, 3)
    e30r, e30i = entry(3, 0)
    e31r, e31i = entry(3, 1)
    e33r, e33i = entry(3, 3)
    complex_trace_form_real = (
        e00r + root_two * e02r + 2 * e11r - root_two * e13r - 2 * e22r - e33r
    )
    complex_trace_form_imag = (
        e00i - root_two * e02i + 2 * e11i + root_two * e13i - 2 * e22i - e33i
    )
    target = sp.expand(
        -sp.Rational(1, 2) * (complex_trace_form_real**2 + complex_trace_form_imag**2)
        - 4 * (e01r - e23r) ** 2
        - sp.Rational(4, 9) * (e01r - 2 * root_two * e12r + e23r) ** 2
        - sp.Rational(2, 9) * ((e20r + e31r) ** 2 + (e20i + e31i) ** 2)
        - sp.Rational(52, 9) * (e30r**2 + e30i**2)
    )
    residual = sp.simplify(sp.expand_complex(optimum - target))
    if residual != 0:
        raise AssertionError(f"nonzero p=4 residual: {residual}")
    print("PASS p=4 second-order identity")
    print("e(E) is the five-term negative sum of squares in L64")


if __name__ == "__main__":
    main()
