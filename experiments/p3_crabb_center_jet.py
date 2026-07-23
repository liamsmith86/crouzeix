#!/usr/bin/env python3
"""Regenerate the flat equality-center jet following the L70 leading form."""

from __future__ import annotations

import sympy as sp

from p3_sparse_series import optimized_rank_one_condition


def base_and_mode_one() -> tuple[sp.Matrix, sp.Matrix]:
    """Return the canonical block and normalized residual mode one."""

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
    return base, mode_one


def mode_two(value: sp.Expr) -> sp.Matrix:
    """Return the residual mode-two matrix with coefficient ``value``."""

    result = sp.zeros(3)
    result[1, 0] = result[2, 1] = value
    return result


def bottom(value: sp.Expr) -> sp.Matrix:
    """Return the bottom-left slice direction with coefficient ``value``."""

    result = sp.zeros(3)
    result[2, 0] = value
    return result


def assert_zeros(coefficients: tuple[sp.Expr, ...], stop: int) -> None:
    """Assert vanishing of all positive orders strictly below ``stop``."""

    if any(coefficients[degree] != 0 for degree in range(1, stop)):
        raise AssertionError(f"the center jet was not flat below order {stop}")


def main() -> None:
    base, mode_one = base_and_mode_one()
    root_two = sp.sqrt(2)
    u0 = 3 * root_two / 64
    v0 = -sp.Rational(9, 64)

    center_path = [base, mode_one, mode_two(u0), bottom(v0)]
    center_ten, _ = optimized_rank_one_condition(center_path, 10)
    assert_zeros(center_ten, 10)
    fixed_tenth = -sp.Rational(13851, 4194304)
    if center_ten[10] != fixed_tenth:
        raise AssertionError("unexpected fixed-center tenth-order coefficient")

    x = sp.symbols("x", real=True)
    mode_path = [base, mode_one, mode_two(x), bottom(v0)]
    bottom_path = [base, mode_one, mode_two(u0), bottom(x)]
    mode_eight, _ = optimized_rank_one_condition(mode_path, 8)
    bottom_eight, _ = optimized_rank_one_condition(bottom_path, 8)
    expected_mode_eight = -(
        4194304 * x**4
        - 3014656 * root_two * x**3
        - 909312 * x**2
        + 263232 * root_two * x
        - 19521
    ) / 1048576
    expected_bottom_eight = (
        13 * (64 * x + 9) * (1984 * x + 117) / 524288
    )
    if sp.expand(mode_eight[8] - expected_mode_eight) != 0:
        raise AssertionError("unexpected mode-coordinate eighth-order polynomial")
    if sp.expand(bottom_eight[8] - expected_bottom_eight) != 0:
        raise AssertionError("unexpected bottom-coordinate eighth-order polynomial")

    gradient = sp.Matrix(
        [
            sp.diff(mode_eight[8], x).subs(x, u0),
            sp.diff(bottom_eight[8], x).subs(x, v0),
        ]
    )
    expected_gradient = sp.Matrix(
        [-sp.Rational(1107, 8192) * root_two, -sp.Rational(1053, 4096)]
    )
    if gradient != expected_gradient:
        raise AssertionError("unexpected eighth-order center gradient")

    leading_quadratic = sp.Matrix(
        [
            [sp.Rational(31, 8), 3 * root_two],
            [3 * root_two, sp.Rational(21, 4)],
        ]
    )
    center_shift = sp.simplify(leading_quadratic.inv() * gradient / 2)
    expected_shift = sp.Matrix(
        [sp.Rational(27, 2048) * root_two, -sp.Rational(81, 2048)]
    )
    if center_shift != expected_shift:
        raise AssertionError("unexpected order-two center shift")
    completion = sp.factor(
        (gradient.T * leading_quadratic.inv() * gradient)[0] / 4
    )
    if sp.factor(fixed_tenth + completion) != 0:
        raise AssertionError("the tenth-order center completion did not cancel")

    corrected_path = [
        base,
        mode_one,
        mode_two(u0),
        bottom(v0),
        mode_two(center_shift[0]),
        bottom(center_shift[1]),
    ]
    corrected_twelve, _ = optimized_rank_one_condition(corrected_path, 12)
    assert_zeros(corrected_twelve, 13)

    reversal = sp.zeros(3)
    reversal[0, 2] = reversal[1, 1] = reversal[2, 0] = 1
    superdiagonal_difference = sp.zeros(3)
    superdiagonal_difference[0, 1] = 1
    superdiagonal_difference[1, 2] = -1
    if any(
        reversal * coefficient.T * reversal != coefficient
        for coefficient in corrected_path
    ):
        raise AssertionError("the corrected center lost reversal-transpose symmetry")
    if (
        reversal * superdiagonal_difference.T * reversal
        != -superdiagonal_difference
    ):
        raise AssertionError("the last transverse direction did not reverse sign")

    print("PASS p=3 Crabb equality-center jet")
    print(f"fixed-center tenth coefficient = {fixed_tenth}")
    print(f"eighth-order gradient = {list(gradient)}")
    print(f"order-two center shift = {list(center_shift)}")
    print(f"completion = {completion}")
    print("corrected certificate is flat through order twelve")
    print("reversal-transpose symmetry makes the remaining transverse term -8*s^2")


if __name__ == "__main__":
    main()
