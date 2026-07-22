#!/usr/bin/env python3
"""Exact one-parameter Riemann-map series for a real 3x3 matrix path.

The routines here work in the physical gauge used by the Crabb local analysis:
the inverse Riemann map ``Psi_e`` is normalized by ``Psi_e(0)=0`` and a positive
derivative at zero.  They are intentionally symbolic and small-dimensional;
this is a certificate helper, not a general numerical conformal-map package.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp


def laurent_modes(expression: sp.Expr, variable: sp.Symbol) -> dict[int, sp.Expr]:
    """Return the nonzero Laurent coefficients of ``expression``."""

    modes: dict[int, sp.Expr] = {}
    for term in sp.Add.make_args(sp.expand(expression)):
        degree = int(term.as_powers_dict().get(variable, 0))
        coefficient = sp.simplify(term / variable**degree)
        modes[degree] = sp.simplify(modes.get(degree, 0) + coefficient)
    return {degree: value for degree, value in modes.items() if value != 0}


def compose_laurent(
    expression: sp.Expr,
    variable: sp.Symbol,
    argument: sp.Expr,
) -> sp.Expr:
    """Compose a finite Laurent polynomial with ``argument``."""

    return sp.expand(
        sum(
            coefficient * argument**degree
            for degree, coefficient in laurent_modes(expression, variable).items()
        )
    )


def matrix_series_product(
    left: Sequence[sp.Matrix],
    right: Sequence[sp.Matrix],
    order: int,
) -> list[sp.Matrix]:
    """Multiply two square-matrix power series through ``order``."""

    dimension = left[0].rows
    product = [sp.zeros(dimension) for _ in range(order + 1)]
    for degree in range(order + 1):
        for left_degree in range(degree + 1):
            right_degree = degree - left_degree
            if left_degree < len(left) and right_degree < len(right):
                product[degree] += left[left_degree] * right[right_degree]
    return product


def matrix_series_power(
    series: Sequence[sp.Matrix],
    exponent: int,
    order: int,
) -> list[sp.Matrix]:
    """Raise a square-matrix power series to a nonnegative integer power."""

    dimension = series[0].rows
    value = [sp.eye(dimension), *[sp.zeros(dimension) for _ in range(order)]]
    for _ in range(exponent):
        value = matrix_series_product(value, series, order)
    return value


def polynomial_at_matrix_series(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    series: Sequence[sp.Matrix],
    order: int,
) -> list[sp.Matrix]:
    """Evaluate an ordinary polynomial at a matrix power series."""

    dimension = series[0].rows
    value = [sp.zeros(dimension) for _ in range(order + 1)]
    for exponent, coefficient in laurent_modes(polynomial, variable).items():
        if exponent < 0:
            raise ValueError("the Riemann-map series must be analytic at zero")
        power = matrix_series_power(series, exponent, order)
        for degree in range(order + 1):
            value[degree] += coefficient * power[degree]
    return value


def inverse_riemann_series(
    base: sp.Matrix,
    perturbation: sp.Matrix,
    order: int,
    epsilon: sp.Symbol,
    variable: sp.Symbol,
) -> tuple[list[sp.Expr], list[sp.Matrix]]:
    """Derive ``Psi_e`` and ``phi_e(base+eE)`` through ``order``.

    ``base`` must be a real 3x3 Crabb-normalized matrix whose top support
    eigenvalue is one.  The boundary is first represented in its normal-angle
    parameter, then formally reparameterized until every nonanalytic Fourier
    mode vanishes.
    """

    if base.shape != (3, 3) or perturbation.shape != (3, 3):
        raise ValueError("this exact helper is restricted to 3x3 paths")
    if any(entry.has(sp.I) for entry in (*base, *perturbation)):
        raise ValueError("this helper expects real symbolic matrix paths")

    matrix = base + epsilon * perturbation
    support = (matrix / variable + variable * matrix.T) / 2
    spectral_parameter = sp.symbols("lambda")
    support_value = sp.Integer(1)
    for degree in range(1, order + 1):
        unknown = sp.symbols(f"support_{degree}")
        determinant = (spectral_parameter * sp.eye(3) - support).det()
        equation = sp.expand(
            sp.series(
                determinant.subs(
                    spectral_parameter,
                    support_value + unknown * epsilon**degree,
                ),
                epsilon,
                0,
                degree + 1,
            ).removeO()
        ).coeff(epsilon, degree)
        solution = sp.solve(equation, unknown)
        if len(solution) != 1:
            raise AssertionError("the simple top support root did not lift uniquely")
        support_value = sp.expand(support_value + solution[0] * epsilon**degree)

    boundary_normal = sp.expand(
        variable * (support_value - variable * sp.diff(support_value, variable))
    )

    def boundary_star(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.conjugate(expression).subs(
                {
                    sp.conjugate(variable): 1 / variable,
                    sp.conjugate(epsilon): epsilon,
                }
            )
        )

    angle_shift = sp.Integer(0)
    inverse_map = [variable]
    for degree in range(1, order + 1):
        normal_argument = variable * sp.series(
            sp.exp(sp.I * angle_shift),
            epsilon,
            0,
            degree + 1,
        ).removeO()
        unmatched = compose_laurent(
            boundary_normal,
            variable,
            normal_argument,
        )
        unmatched = sp.expand(
            sp.series(unmatched, epsilon, 0, degree + 1).removeO()
        ).coeff(epsilon, degree)

        shift_coefficient = sp.Integer(0)
        for mode, coefficient in laurent_modes(unmatched, variable).items():
            if mode <= 0:
                analytic_half = sp.I * coefficient * variable ** (mode - 1)
                shift_coefficient += analytic_half + boundary_star(analytic_half)
        angle_shift = sp.expand(
            angle_shift + epsilon**degree * shift_coefficient
        )

        normal_argument = variable * sp.series(
            sp.exp(sp.I * angle_shift),
            epsilon,
            0,
            degree + 1,
        ).removeO()
        analytic_boundary = compose_laurent(
            boundary_normal,
            variable,
            normal_argument,
        )
        coefficient = sp.expand(
            sp.series(analytic_boundary, epsilon, 0, degree + 1).removeO()
        ).coeff(epsilon, degree)
        bad_modes = {
            mode: value
            for mode, value in laurent_modes(coefficient, variable).items()
            if mode <= 0 and sp.simplify(value) != 0
        }
        if bad_modes:
            raise AssertionError(f"nonanalytic boundary modes remain: {bad_modes}")
        inverse_map.append(sp.collect(coefficient, variable))

    operator = [base]
    for degree in range(1, order + 1):
        trial = [*operator, sp.zeros(3)]
        known = sp.zeros(3)
        for map_degree in range(1, degree + 1):
            evaluated = polynomial_at_matrix_series(
                inverse_map[map_degree],
                variable,
                trial,
                degree - map_degree,
            )
            known += evaluated[degree - map_degree]
        target = perturbation if degree == 1 else sp.zeros(3)
        operator.append(sp.simplify(target - known))

    return inverse_map, operator
