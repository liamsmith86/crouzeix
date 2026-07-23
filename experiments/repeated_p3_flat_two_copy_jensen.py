#!/usr/bin/env python3
"""Prove the quantitative two-copy trace-splitting Jensen gap."""

from __future__ import annotations

import sympy as sp


def constant_laurent_coefficient(
    expression: sp.Expr,
    boundary: sp.Symbol,
    shift: int = 8,
) -> sp.Expr:
    """Return the constant coefficient of a bounded Laurent polynomial."""

    shifted = sp.expand(sp.cancel(expression) * boundary**shift)
    return shifted.coeff(boundary, shift)


def main() -> None:
    boundary = sp.symbols("boundary", nonzero=True)
    edge = sp.symbols("edge", real=True, nonnegative=True)
    diagonal_real, diagonal_imaginary = sp.symbols(
        "diagonal_real diagonal_imaginary",
        real=True,
    )
    scalar_real, scalar_imaginary = sp.symbols(
        "scalar_real scalar_imaginary",
        real=True,
    )
    diagonal = diagonal_real + sp.I * diagonal_imaginary
    scalar = scalar_real + sp.I * scalar_imaginary
    traceless = sp.Matrix([[diagonal, edge], [0, -diagonal]])
    traceless_adjoint = traceless.conjugate().T

    mean_part = sp.Rational(5, 64) * (
        scalar * traceless_adjoint
        + sp.conjugate(scalar) * traceless
    )
    angular_part = -sp.Rational(3, 64) * (
        boundary**2 * scalar * traceless
        + boundary**-2
        * sp.conjugate(scalar)
        * traceless_adjoint
    )

    def vector_inner(first: sp.Matrix, second: sp.Matrix) -> sp.Expr:
        """Return the Euclidean inner product of Pauli vectors."""

        return sp.trace(first * second) / 2

    mean_norm_square = sp.factor(vector_inner(mean_part, mean_part))
    angular_norm_square_mean = sp.factor(
        constant_laurent_coefficient(
            vector_inner(angular_part, angular_part),
            boundary,
        )
    )
    mixed_inner = sp.expand(vector_inner(mean_part, angular_part))
    mixed_square_mean = sp.factor(
        constant_laurent_coefficient(
            mixed_inner**2,
            boundary,
        )
    )
    perpendicular_energy = sp.factor(
        angular_norm_square_mean
        - mixed_square_mean / mean_norm_square
    )

    frobenius_square = (
        2 * diagonal * sp.conjugate(diagonal) + edge**2
    )
    square_trace = 2 * diagonal**2
    base_norm_square = (
        scalar * sp.conjugate(scalar) * frobenius_square
        + (
            sp.conjugate(scalar) ** 2 * square_trace
            + scalar**2 * sp.conjugate(square_trace)
        )
        / 2
    )
    nonnormal_discriminant = sp.factor(
        frobenius_square**2
        - square_trace * sp.conjugate(square_trace)
    )
    expected_discriminant = edge**2 * (
        4 * diagonal * sp.conjugate(diagonal) + edge**2
    )
    expected_perpendicular_energy = (
        9
        * (scalar * sp.conjugate(scalar)) ** 2
        * nonnormal_discriminant
        / (8192 * base_norm_square)
    )

    if sp.simplify(
        mean_norm_square
        - sp.Rational(25, 4096) * base_norm_square
    ) != 0:
        raise AssertionError("the mean Pauli-vector norm failed")
    if sp.simplify(
        nonnormal_discriminant - expected_discriminant
    ) != 0:
        raise AssertionError("the nonnormal discriminant failed")
    if sp.simplify(
        perpendicular_energy - expected_perpendicular_energy
    ) != 0:
        raise AssertionError("the perpendicular Jensen energy failed")

    print("PASS repeated p=3 two-copy trace-splitting Jensen gap")
    print("Z = scalar*I + [[diagonal,edge],[0,-diagonal]]")
    print(
        "mean ||C_perp||^2 = "
        "9*|scalar|^4*edge^2*(4|diagonal|^2+edge^2)"
        "/(8192*base_norm_square)"
    )
    print("the expression is strict exactly off scalar=0 or edge=0")


if __name__ == "__main__":
    main()
