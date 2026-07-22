#!/usr/bin/env python3
"""Regenerate the exact sixth-order certificate on the flat p=3 Crabb mode.

After quotienting the affine-unitary equality orbit, one complex mode-one
direction remains.  Circle symmetry makes its coefficient real, giving the
path ``C_3 + epsilon*E`` below.  This script derives the Riemann pullback
through order six and constructs an explicit rank-one-defect Stein metric.
Its condition number is ``4 - 171*epsilon**6/4096 + O(epsilon**7)``.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from crabb_second_order_symbolic import conformal_coefficients
from formal_riemann_series import inverse_riemann_series


ORDER = 6


def matrix_power_series(
    coefficients: Sequence[sp.Matrix],
    epsilon: sp.Symbol,
) -> sp.Matrix:
    """Assemble a matrix power series as a symbolic matrix."""

    value = sp.zeros(coefficients[0].rows)
    for degree, coefficient in enumerate(coefficients):
        value += epsilon**degree * coefficient
    return value


def stein_gramian_series(
    operator: Sequence[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
) -> list[sp.Matrix]:
    """Solve ``P-T.T*P*T=defect*defect.T`` through ``ORDER``."""

    base = operator[0]
    defect_square = sp.expand(defect * defect.T)
    forcing = [
        defect_square.applyfunc(lambda entry: entry.coeff(epsilon, degree))
        for degree in range(ORDER + 1)
    ]
    gramian: list[sp.Matrix] = []
    for degree in range(ORDER + 1):
        right_hand_side = forcing[degree]
        for left_degree in range(degree + 1):
            for right_degree in range(degree + 1 - left_degree):
                metric_degree = degree - left_degree - right_degree
                if metric_degree < len(gramian):
                    right_hand_side += (
                        operator[left_degree].T
                        * gramian[metric_degree]
                        * operator[right_degree]
                    )
        coefficient = sp.zeros(3)
        for power in range(3):
            coefficient += (
                (base.T**power) * right_hand_side * (base**power)
            )
        gramian.append(sp.simplify(coefficient))
    return gramian


def simple_eigenvalue_series(
    matrix: sp.Matrix,
    base_eigenvalue: int,
    epsilon: sp.Symbol,
) -> sp.Expr:
    """Lift a simple eigenvalue of the constant matrix through ``ORDER``."""

    spectral_parameter = sp.symbols("spectral_parameter")
    series = sp.Integer(base_eigenvalue)
    characteristic = (spectral_parameter * sp.eye(3) - matrix).det()
    for degree in range(1, ORDER + 1):
        unknown = sp.symbols(f"eigenvalue_{base_eigenvalue}_{degree}")
        equation = sp.expand(
            sp.series(
                characteristic.subs(
                    spectral_parameter,
                    series + unknown * epsilon**degree,
                ),
                epsilon,
                0,
                degree + 1,
            ).removeO()
        ).coeff(epsilon, degree)
        solution = sp.solve(equation, unknown)
        if len(solution) != 1:
            raise AssertionError("a simple Gramian eigenvalue did not lift uniquely")
        series = sp.expand(series + solution[0] * epsilon**degree)
    return series


def check_second_order_regression(
    base: sp.Matrix,
    perturbation: sp.Matrix,
    operator: Sequence[sp.Matrix],
    variable: sp.Symbol,
) -> None:
    """Compare the new formal map with the independent L63 derivation."""

    root_two = sp.sqrt(2)
    top_vector = sp.Matrix([sp.Rational(1, 2), variable / root_two, variable**2 / 2])
    zero_vector = sp.Matrix([1 / root_two, 0, -(variable**2) / root_two])
    bottom_vector = sp.Matrix(
        [sp.Rational(1, 2), -variable / root_two, variable**2 / 2]
    )
    tangent, second_order = conformal_coefficients(
        base,
        perturbation,
        top_vector,
        [(sp.Integer(0), zero_vector), (-sp.Integer(1), bottom_vector)],
        variable,
    )
    if operator[1] != tangent or operator[2] != second_order:
        raise AssertionError("formal map disagrees with the independent L63 expansion")


def main() -> None:
    epsilon = sp.symbols("epsilon", real=True)
    variable = sp.symbols("z", nonzero=True)
    root_two = sp.sqrt(2)
    base = sp.Matrix([[0, root_two, 0], [0, 0, root_two], [0, 0, 0]])
    perturbation = sp.diag(
        -sp.Rational(1, 3),
        sp.Rational(2, 3),
        -sp.Rational(1, 3),
    )
    perturbation[0, 2] = 1

    inverse_map, operator = inverse_riemann_series(
        base,
        perturbation,
        ORDER,
        epsilon,
        variable,
    )
    check_second_order_regression(base, perturbation, operator, variable)
    if inverse_map[1] != 5 * variable**2 / 12:
        raise AssertionError("unexpected first inverse-map coefficient")
    if inverse_map[2] != 73 * variable**3 / 576 - 55 * variable / 576:
        raise AssertionError("unexpected second inverse-map coefficient")

    a1, a3, a5 = sp.symbols("a1 a3 a5", real=True)
    b2, b4, b6 = sp.symbols("b2 b4 b6", real=True)
    defect = sp.Matrix(
        [
            1,
            a1 * epsilon + a3 * epsilon**3 + a5 * epsilon**5,
            b2 * epsilon**2 + b4 * epsilon**4 + b6 * epsilon**6,
        ]
    )
    gramian = stein_gramian_series(operator, defect, epsilon)
    full_gramian = matrix_power_series(gramian, epsilon)
    lower = simple_eigenvalue_series(full_gramian, 1, epsilon)
    upper = simple_eigenvalue_series(full_gramian, 4, epsilon)
    condition = sp.expand(
        sp.series(upper / lower, epsilon, 0, ORDER + 1).removeO()
    )

    second = sp.factor(condition.coeff(epsilon, 2))
    expected_second = (288 * a1**2 + 24 * root_two * a1 + 1) / 36
    if sp.simplify(second - expected_second) != 0:
        raise AssertionError("unexpected quadratic condition-number coefficient")
    a1_star = -root_two / 24
    if second.subs(a1, a1_star) != 0:
        raise AssertionError("the second-order coefficient did not vanish")

    fourth = sp.factor(condition.coeff(epsilon, 4).subs(a1, a1_star))
    expected_fourth = (576 * b2 + 193) ** 2 / 124416
    if sp.simplify(fourth - expected_fourth) != 0:
        raise AssertionError("unexpected fourth-order condition-number coefficient")
    b2_star = -sp.Rational(193, 576)
    if fourth.subs(b2, b2_star) != 0:
        raise AssertionError("the fourth-order coefficient did not vanish")

    sixth = sp.factor(
        condition.coeff(epsilon, 6).subs({a1: a1_star, b2: b2_star})
    )
    expected_sixth = (
        10616832 * a3**2 - 354816 * root_two * a3 - 49475
    ) / 1327104
    if sp.simplify(sixth - expected_sixth) != 0:
        raise AssertionError("unexpected sixth-order condition-number coefficient")
    a3_star = 77 * root_two / 4608
    final_coefficient = sp.factor(sixth.subs(a3, a3_star))
    if final_coefficient != -sp.Rational(171, 4096):
        raise AssertionError("the sixth-order certificate did not have the claimed sign")

    print("PASS p=3 Crabb residual mode-one sixth-order certificate")
    print(f"a1 = {a1_star}, b2 = {b2_star}, a3 = {a3_star}")
    print(f"condition(P) = 4 + ({final_coefficient})*epsilon^6 + O(epsilon^7)")


if __name__ == "__main__":
    main()
