#!/usr/bin/env python3
"""Exact algebra audit for the two extreme-orientation transfer faces.

The elementary sign proof is recorded in
``proof/slice_core_projective_reduction.md`` section 4.  This checker
regenerates its scalar determinant factor and the identities used after the
``q = |ab|`` reduction.  It uses exact SymPy arithmetic only.
"""

from __future__ import annotations

import sympy as sp


def assert_zero(expression: sp.Expr, label: str) -> None:
    """Raise unless an exact symbolic identity vanishes."""

    if sp.factor(expression) != 0:
        raise AssertionError(label)


def main() -> None:
    c, x, a, b, q, radius = sp.symbols("c x a b q radius", nonnegative=True)
    modal, p = sp.symbols("modal p", nonnegative=True)
    odd_defect = 1 - a**2
    even_defect = 1 - b**2

    odd_block = even_defect * (4 - a**2) + odd_defect * (4 * b**2 - 1) * c * x
    even_block = odd_defect * (4 - b**2) + even_defect * (4 * a**2 - 1) * x / c
    coupling_square = x / c * (a * even_defect + b * c * odd_defect) ** 2
    scalar_determinant = odd_block * even_block - 9 * coupling_square

    g_factor = (
        4 * a**2 * b**2 * c**2 * x
        - 16 * a**2 * b**2 * c * x**2
        - a**2 * b**2 * c
        + 4 * a**2 * b**2 * x
        - 4 * (a**2 + b**2) * c**2 * x
        + 4 * (a**2 + b**2) * c * x**2
        + 4 * (a**2 + b**2) * c
        - 4 * (a**2 + b**2) * x
        + 18 * a * b * c * x
        + 4 * c**2 * x
        - c * x**2
        - 16 * c
        + 4 * x
    )
    h_factor = (
        4 * a**2 * b**2 * c * x
        - a**2 * b**2
        - a**2 * c * x
        + a**2
        - 4 * b**2 * c * x
        + 4 * b**2
        + c * x
        - 4
    )
    assert_zero(odd_block + h_factor, "odd diagonal factor failed")
    assert_zero(
        scalar_determinant + odd_defect * even_defect * g_factor / c,
        "scalar determinant factor failed",
    )

    constant = -(4 * c - x) * (4 - c * x)
    linear = 4 * (c - x) * (1 - c * x)
    quadratic = 4 * c**2 * x - 16 * c * x**2 - c + 4 * x
    reduced = constant + linear * (a**2 + b**2) + quadratic * (a * b) ** 2
    reduced += 18 * a * b * c * x
    assert_zero(g_factor - reduced, "symmetric scalar reduction failed")

    boundary = sp.factor(g_factor.subs({a: 1, b: q}))
    expected_boundary = -3 * c * (2 * q * x - q + x - 2) * (2 * q * x + q - x - 2)
    assert_zero(boundary - expected_boundary, "boundary factor failed")

    diagonal_upper = constant + 2 * linear * q + quadratic * q**2
    diagonal_upper += 18 * c * x * q
    expected_diagonal = (
        -9 * c * (1 - x) ** 2
        + 6 * c * (4 * x + 1) * (x - 1) * radius
        + quadratic * radius**2
    )
    assert_zero(
        diagonal_upper.subs(q, 1 - radius) - expected_diagonal,
        "diagonal upper bound identity failed",
    )

    corners = tuple(
        sp.factor(h_factor.subs({a: odd, b: even}))
        for odd, even in ((0, 0), (1, 0), (0, 1), (1, 1))
    )
    if corners != (c * x - 4, -3, -3 * c * x, 0):
        raise AssertionError(f"unexpected H corners: {corners}")

    def substitute_x(expression: sp.Expr, value: sp.Expr) -> sp.Expr:
        return expression.subs(x, value)

    def endpoint_core(first: sp.Expr, second: sp.Expr) -> sp.Matrix:
        """Return the endpoint core after a positive diagonal congruence."""

        first_square = first**2
        second_square = second**2
        odd_first = substitute_x(odd_block, first_square)
        odd_second = substitute_x(odd_block, second_square)
        even_first = c * substitute_x(even_block, first_square)
        even_second = c * substitute_x(even_block, second_square)
        common = a * even_defect + b * c * odd_defect
        coupling_first = -3 * first * common
        coupling_second = -3 * second * common
        return sp.Matrix(
            [
                [odd_first, 0, coupling_first, 0],
                [0, odd_second, 0, coupling_second],
                [coupling_first, 0, even_first, 0],
                [0, coupling_second, 0, even_second],
            ]
        )

    def audit_endpoint(first: sp.Expr, second: sp.Expr, label: str) -> None:
        core = endpoint_core(first, second)
        first_square = first**2
        second_square = second**2
        leading_expected = (
            odd_defect
            * even_defect
            * substitute_x(h_factor, second_square)
            * substitute_x(g_factor, first_square)
        )
        determinant_expected = (
            odd_defect**2
            * even_defect**2
            * substitute_x(g_factor, first_square)
            * substitute_x(g_factor, second_square)
        )
        assert_zero(
            core[:3, :3].det() - leading_expected,
            f"{label} leading minor factor failed",
        )
        assert_zero(
            core.det() - determinant_expected,
            f"{label} determinant factor failed",
        )

    # At the two extreme orientations the modal amplitudes are respectively
    # (sqrt(k), p*sqrt(k)) and their interchange.  ``modal`` represents
    # sqrt(k), avoiding radicals while checking the complete 4-by-4 core.
    audit_endpoint(modal, modal * p, "orientation zero")
    audit_endpoint(modal * p, modal, "orientation one")

    print("extreme-orientation factor identities: exact")


if __name__ == "__main__":
    main()
