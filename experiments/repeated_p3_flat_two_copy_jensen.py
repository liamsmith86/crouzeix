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

    common_real, common_imaginary = sp.symbols(
        "common_real common_imaginary",
        real=True,
    )
    common = common_real + sp.I * common_imaginary
    common_linear_coefficient = sp.sqrt(2) * (
        boundary**-1 * common
        - boundary**3 * sp.conjugate(common)
    ) / 16
    common_linear_adjoint = sp.sqrt(2) * (
        boundary * sp.conjugate(common)
        - boundary**-3 * common
    ) / 16
    common_support = (
        common_linear_coefficient * traceless
        + common_linear_adjoint * traceless_adjoint
    )
    common_cross = sp.expand(common_support[0, 1])
    expected_common_cross = (
        sp.sqrt(2)
        * edge
        * (
            boundary**-1 * common
            - boundary**3 * sp.conjugate(common)
        )
        / 16
    )
    if sp.simplify(common_cross - expected_common_cross) != 0:
        raise AssertionError("the common-mode cross coefficient failed")
    common_support_mean = common_support.applyfunc(
        lambda entry: constant_laurent_coefficient(entry, boundary)
    )
    if common_support_mean != sp.zeros(2):
        raise AssertionError("the common-mode traceless mean did not vanish")
    even_support = mean_part + angular_part
    if sp.simplify(
        even_support.subs(boundary, -boundary) - even_support
    ) != sp.zeros(2):
        raise AssertionError("the trace-splitting support lost even parity")
    if sp.simplify(
        common_support.subs(boundary, -boundary) + common_support
    ) != sp.zeros(2):
        raise AssertionError("the common-mode support lost odd parity")

    print("PASS repeated p=3 two-copy trace-splitting Jensen gap")
    print("Z = scalar*I + [[diagonal,edge],[0,-diagonal]]")
    print(
        "mean ||C_perp||^2 = "
        "9*|scalar|^4*edge^2*(4|diagonal|^2+edge^2)"
        "/(8192*base_norm_square)"
    )
    print("the expression is strict exactly off scalar=0 or edge=0")
    print("on scalar=0, common w gives endpoint <= -sqrt(2)*edge*|w|")
    print("even/odd pairing prevents common-w cancellation for scalar != 0")


if __name__ == "__main__":
    main()
