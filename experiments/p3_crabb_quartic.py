#!/usr/bin/env python3
"""Regenerate the quartic rank-one Stein certificate on the p=3 quotient.

Circle symmetry and conjugation leave only three real quartic invariants in
the residual coordinates ``(z,w)``.  Exact calculations at ``(1,0)``,
``(0,1)``, and ``(1,1)`` therefore determine the complete invariant form.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from formal_riemann_series import inverse_riemann_series
from rank_one_stein_series import gramian_condition_series


ORDER = 4


@dataclass(frozen=True)
class PathCertificate:
    first_defect: tuple[sp.Expr, sp.Expr]
    second_defect: tuple[sp.Expr, sp.Expr]
    cubic: sp.Expr
    quartic: sp.Expr


def residual_perturbation(mode_one: int, mode_two: int) -> sp.Matrix:
    """Return the real residual representative ``R1(z)+R2(w)``."""

    z = sp.Integer(mode_one)
    w = sp.Integer(mode_two)
    perturbation = sp.diag(-z / 3, 2 * z / 3, -z / 3)
    perturbation[0, 2] = z
    perturbation[1, 0] = w
    perturbation[2, 1] = w
    return perturbation


def certify_path(mode_one: int, mode_two: int) -> PathCertificate:
    """Optimize the formal Stein defect through fourth order on one path."""

    epsilon = sp.symbols("epsilon", real=True)
    variable = sp.symbols("boundary_variable", nonzero=True)
    root_two = sp.sqrt(2)
    base = sp.Matrix([[0, root_two, 0], [0, 0, root_two], [0, 0, 0]])
    _, operator = inverse_riemann_series(
        base,
        residual_perturbation(mode_one, mode_two),
        ORDER,
        epsilon,
        variable,
    )

    a1, a2, a3 = sp.symbols("a1 a2 a3", real=True)
    b1, b2, b3 = sp.symbols("b1 b2 b3", real=True)
    defect = sp.Matrix(
        [
            1,
            a1 * epsilon + a2 * epsilon**2 + a3 * epsilon**3,
            b1 * epsilon + b2 * epsilon**2 + b3 * epsilon**3,
        ]
    )
    condition = gramian_condition_series(operator, defect, epsilon, ORDER)

    quadratic = sp.factor(condition.coeff(epsilon, 2))
    first_solution = sp.solve(
        [sp.diff(quadratic, a1), sp.diff(quadratic, b1)],
        [a1, b1],
        dict=True,
    )
    if len(first_solution) != 1:
        raise AssertionError("the first defect correction was not unique")
    first = first_solution[0]
    if sp.simplify(quadratic.subs(first)) != 0:
        raise AssertionError("the quotient path was not second-order flat")

    cubic = sp.factor(condition.coeff(epsilon, 3).subs(first))
    quartic_before_second = sp.factor(
        condition.coeff(epsilon, 4).subs(first)
    )
    second_solution = sp.solve(
        [
            sp.diff(quartic_before_second, a2),
            sp.diff(quartic_before_second, b2),
        ],
        [a2, b2],
        dict=True,
    )
    if len(second_solution) != 1:
        raise AssertionError("the second defect correction was not unique")
    second = second_solution[0]
    quartic = sp.factor(quartic_before_second.subs(second))
    return PathCertificate(
        first_defect=(first[a1], first[b1]),
        second_defect=(second[a2], second[b2]),
        cubic=cubic,
        quartic=quartic,
    )


def main() -> None:
    root_two = sp.sqrt(2)
    mode_one = certify_path(1, 0)
    mode_two = certify_path(0, 1)
    mixed = certify_path(1, 1)

    expected = {
        "mode_one": PathCertificate(
            (-root_two / 24, 0),
            (0, -sp.Rational(193, 576)),
            sp.Integer(0),
            sp.Integer(0),
        ),
        "mode_two": PathCertificate(
            (0, -2 * root_two),
            (0, 0),
            sp.Integer(0),
            -sp.Integer(4),
        ),
        "mixed": PathCertificate(
            (-root_two / 24, -2 * root_two),
            (sp.Rational(1, 2), -sp.Rational(193, 576)),
            sp.Integer(0),
            -sp.Rational(63, 8),
        ),
    }
    actual = {"mode_one": mode_one, "mode_two": mode_two, "mixed": mixed}
    if actual != expected:
        raise AssertionError(f"unexpected path certificates: {actual}")

    mode_one_fourth = mode_one.quartic
    mode_two_fourth = mode_two.quartic
    mixed_fourth = mixed.quartic
    cross_fourth = sp.factor(
        mixed_fourth - mode_one_fourth - mode_two_fourth
    )
    if cross_fourth != -sp.Rational(31, 8):
        raise AssertionError("unexpected mixed quartic invariant")

    print("PASS p=3 Crabb residual quartic certificate")
    print("cubic(z,w) = 0")
    print("quartic(z,w) = -4*|w|^4 - 31/8*|z|^2*|w|^2")


if __name__ == "__main__":
    main()
