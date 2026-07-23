#!/usr/bin/env python3
"""Regenerate the leading weighted certificate on the p=3 Crabb slice.

The residual mode two has weight two and the negative slice directions have
weight three relative to the mode-one parameter.  Low ordinary orders on a
small set of integer paths determine every invariant in the weighted
sixth-order form.
"""

from __future__ import annotations

import sympy as sp

from p3_sparse_series import (
    optimized_rank_one_condition,
)


def base_and_directions() -> tuple[sp.Matrix, ...]:
    """Return ``C3`` and the four real local-slice generators."""

    root_two = sp.sqrt(2)
    base = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    mode_one = sp.diag(
        -sp.Rational(1, 3),
        sp.Rational(2, 3),
        -sp.Rational(1, 3),
    )
    mode_one[0, 2] = 1
    mode_two = sp.zeros(3)
    mode_two[1, 0] = mode_two[2, 1] = 1
    superdiagonal_difference = sp.zeros(3)
    superdiagonal_difference[0, 1] = 1
    superdiagonal_difference[1, 2] = -1
    bottom = sp.zeros(3)
    bottom[2, 0] = 1
    return base, mode_one, mode_two, superdiagonal_difference, bottom


def optimized_branch(perturbation: sp.Matrix, order: int) -> tuple[sp.Expr, ...]:
    """Return optimized condition coefficients on one straight path."""

    base, *_ = base_and_directions()
    coefficients, _ = optimized_rank_one_condition(
        [base, perturbation],
        order,
    )
    return coefficients


def main() -> None:
    _, mode_one, mode_two, superdiagonal_difference, bottom = (
        base_and_directions()
    )
    root_two = sp.sqrt(2)

    pure_mode_one = optimized_branch(mode_one, 6)
    pure_mode_two = optimized_branch(mode_two, 4)
    mixed_modes = optimized_branch(mode_one + mode_two, 4)
    constant = pure_mode_one[6]
    mode_cross = sp.factor(
        mixed_modes[4]
        - pure_mode_one[4]
        - pure_mode_two[4]
    )

    fifth_one = optimized_branch(mode_one + mode_two, 5)[5]
    fifth_two = optimized_branch(mode_one + 2 * mode_two, 5)[5]
    mode_linear = sp.factor((8 * fifth_one - fifth_two) / 6)
    mode_cubic = sp.factor(fifth_one - mode_linear)

    bottom_plus = optimized_branch(mode_one + bottom, 4)[4]
    bottom_minus = optimized_branch(mode_one - bottom, 4)[4]
    bottom_linear = sp.factor((bottom_plus - bottom_minus) / 2)
    bottom_quadratic = optimized_branch(bottom, 2)[2]

    coupling_plus = optimized_branch(
        mode_one + mode_two + bottom,
        3,
    )[3]
    coupling_minus = optimized_branch(
        mode_one + mode_two - bottom,
        3,
    )[3]
    mode_bottom_coupling = sp.factor((coupling_plus - coupling_minus) / 2)
    superdiagonal_quadratic = optimized_branch(
        superdiagonal_difference,
        2,
    )[2]

    expected = {
        "constant": -sp.Rational(171, 4096),
        "mode_cross": -sp.Rational(31, 8),
        "mode_linear": -sp.Rational(123, 256) * root_two,
        "mode_cubic": sp.Rational(23, 8) * root_two,
        "bottom_linear": -sp.Rational(117, 128),
        "bottom_quadratic": -sp.Rational(21, 4),
        "mode_bottom_coupling": -6 * root_two,
        "superdiagonal_quadratic": -8,
    }
    actual = {
        "constant": constant,
        "mode_cross": mode_cross,
        "mode_linear": mode_linear,
        "mode_cubic": mode_cubic,
        "bottom_linear": bottom_linear,
        "bottom_quadratic": bottom_quadratic,
        "mode_bottom_coupling": mode_bottom_coupling,
        "superdiagonal_quadratic": superdiagonal_quadratic,
    }
    if actual != expected:
        raise AssertionError(f"unexpected weighted coefficients: {actual}")

    u, v, s = sp.symbols("u v s", real=True)
    leading = (
        constant
        + mode_cross * u**2
        + mode_linear * u
        + bottom_linear * v
        + mode_bottom_coupling * u * v
        + bottom_quadratic * v**2
        + superdiagonal_quadratic * s**2
    )
    completed = (
        -8 * s**2
        - sp.Rational(21, 4)
        * (v + sp.Rational(39, 448) + 4 * root_two * u / 7) ** 2
        - sp.Rational(25, 56)
        * (u - 3 * root_two / 64) ** 2
    )
    if sp.expand(leading - completed) != 0:
        raise AssertionError("the weighted form did not complete to squares")

    print("PASS p=3 Crabb weighted-slice certificate")
    print(f"ordinary fifth linear coefficient = {mode_linear}")
    print(f"ordinary fifth cubic coefficient = {mode_cubic}")
    print("weighted sixth-order form:")
    print(sp.factor(completed))
    print("unique real center: u=3*sqrt(2)/64, v=-9/64, s=0")


if __name__ == "__main__":
    main()
