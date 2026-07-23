#!/usr/bin/env python3
"""Prove that the L70 equality-center jet is an exact disk-matrix curve."""

from __future__ import annotations

import sympy as sp


def polynomial_remainder(
    expression: sp.Expr,
    relation: sp.Expr,
    variable: sp.Symbol,
) -> sp.Expr:
    """Reduce a rational identity modulo a polynomial relation."""

    numerator = sp.factor(sp.together(expression).as_numer_denom()[0])
    return sp.factor(
        sp.rem(sp.Poly(numerator, variable), sp.Poly(relation, variable)).as_expr()
    )


def main() -> None:
    epsilon, scale = sp.symbols("epsilon scale", real=True)
    x, y, z, spectral = sp.symbols("x y z spectral", real=True)
    root_two = sp.sqrt(2)
    matrix = sp.Matrix(
        [
            [-epsilon / 3, root_two, epsilon],
            [
                3 * root_two * epsilon**2 * scale / 64,
                2 * epsilon / 3,
                root_two,
            ],
            [
                -9 * epsilon**3 * scale / 64,
                3 * root_two * epsilon**2 * scale / 64,
                -epsilon / 3,
            ],
        ]
    )
    disk_relation = (
        81 * scale**2 * epsilon**4
        + (1152 * epsilon**2 - 4096) * scale
        + 4096
    )

    squared_parameter = sp.symbols("squared_parameter", real=True)
    analytic_scale = (
        64
        * (
            32
            - 9 * squared_parameter
            - 8 * sp.sqrt(16 - 9 * squared_parameter)
        )
        / (81 * squared_parameter**2)
    )
    analytic_scale_in_epsilon = analytic_scale.subs(
        squared_parameter,
        epsilon**2,
    )
    expected_scale_series = (
        1
        + sp.Rational(9, 32) * squared_parameter
        + sp.Rational(405, 4096) * squared_parameter**2
        + sp.Rational(5103, 131072) * squared_parameter**3
    )
    actual_scale_series = sp.series(
        analytic_scale,
        squared_parameter,
        0,
        4,
    ).removeO()
    if sp.expand(actual_scale_series - expected_scale_series) != 0:
        raise AssertionError("unexpected analytic disk-center branch")
    if (
        sp.factor(
            disk_relation.subs(
                {
                    epsilon**2: squared_parameter,
                    scale: analytic_scale,
                }
            )
        )
        != 0
    ):
        raise AssertionError("the analytic scale did not solve the disk relation")

    characteristic = sp.Poly(
        sp.expand((spectral * sp.eye(3) - matrix).det()),
        spectral,
    )
    discriminant = sp.factor(sp.discriminant(characteristic.as_expr(), spectral))
    expected_discriminant = (
        -9
        * scale
        * epsilon**6
        * (3 * scale + 64 * epsilon**2 - 128)
        * disk_relation
        / 4194304
    )
    if sp.expand(discriminant - expected_discriminant) != 0:
        raise AssertionError("unexpected characteristic discriminant")

    depressed_quadratic = -characteristic.coeff_monomial(spectral)
    depressed_constant = characteristic.coeff_monomial(1)
    center = sp.factor(3 * depressed_constant / (2 * depressed_quadratic))
    centered_characteristic = (spectral - center) ** 2 * (
        spectral + 2 * center
    )
    if (
        polynomial_remainder(
            characteristic.as_expr() - centered_characteristic,
            disk_relation,
            scale,
        )
        != 0
    ):
        raise AssertionError("the disk center was not a double eigenvalue")

    hermitian_pencil = (
        z * sp.eye(3)
        + x * (matrix + matrix.T) / 2
        + y * (matrix - matrix.T) / (2 * sp.I)
    )
    kippenhahn = sp.Poly(
        sp.expand(hermitian_pencil.det()),
        x,
        y,
        z,
    )

    def coefficient(x_degree: int, y_degree: int, z_degree: int) -> sp.Expr:
        return sp.factor(
            kippenhahn.coeff_monomial(
                x**x_degree * y**y_degree * z**z_degree
            )
        )

    qx = coefficient(2, 0, 1)
    qy = coefficient(0, 2, 1)
    cx = coefficient(3, 0, 0)
    cy = coefficient(1, 2, 0)
    line_coordinate = 2 * center
    identities = (
        cy + line_coordinate * qy,
        cx + line_coordinate * qx + line_coordinate**3,
        qx - qy + sp.Rational(3, 4) * line_coordinate**2,
    )
    if any(
        polynomial_remainder(identity, disk_relation, scale) != 0
        for identity in identities
    ):
        raise AssertionError("the Kippenhahn pencil did not collapse to a circle")

    radius_squared = -qy
    disk_factorization = (z - 2 * center * x) * (
        (z + center * x) ** 2 - radius_squared * (x**2 + y**2)
    )
    if (
        polynomial_remainder(
            kippenhahn.as_expr() - disk_factorization,
            disk_relation,
            scale,
        )
        != 0
    ):
        raise AssertionError("the complete Kippenhahn factorization failed")
    if sp.limit(radius_squared.subs(scale, analytic_scale).subs(epsilon**2, squared_parameter), squared_parameter, 0) != 1:
        raise AssertionError("the disk radius did not converge to one")
    if (
        sp.limit(
            center.subs(scale, analytic_scale_in_epsilon) / epsilon,
            epsilon,
            0,
        )
        != sp.Rational(5, 12)
    ):
        raise AssertionError("the disk center had the wrong tangent")

    print("PASS p=3 Crabb disk-center factorization")
    print(f"disk relation = {disk_relation}")
    print(f"analytic scale = {analytic_scale}")
    print(f"scale series = {actual_scale_series}")
    print("Kippenhahn curve = circle(center) union one interior point")
    print("the circle center is a double eigenvalue")


if __name__ == "__main__":
    main()
