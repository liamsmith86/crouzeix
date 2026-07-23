#!/usr/bin/env python3
"""Derive the second support function on the L74 repeated-C3 cross quotient."""

from __future__ import annotations

import sympy as sp


def generators() -> tuple[tuple[sp.Matrix, sp.Matrix], ...]:
    """Return the three complex L74 cross-generator pairs."""

    e02 = sp.zeros(3)
    e02[0, 2] = 1
    difference = sp.zeros(3)
    difference[0, 1] = -1
    difference[1, 2] = 1
    return (
        (
            sp.diag(
                -sp.Rational(3, 4),
                sp.Rational(1, 4),
                -sp.Rational(3, 4),
            ),
            e02,
        ),
        (difference, difference),
        (-sp.Rational(4, 3) * e02, sp.diag(1, -sp.Rational(1, 3), 1)),
    )


def main() -> None:
    boundary = sp.symbols("boundary", nonzero=True)
    root_two = sp.sqrt(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    support_operator = (crabb / boundary + boundary * crabb.T) / 2
    reduced_resolvent = (
        sp.eye(3)
        - sp.Rational(1, 4) * support_operator
        - sp.Rational(3, 4) * support_operator**2
    )
    support_vector = sp.Matrix([1 / boundary, root_two, boundary]) / 2
    conjugate_support_vector = sp.Matrix(
        [boundary, root_two, 1 / boundary]
    ) / 2

    alpha = sp.symbols("alpha_0:3")
    alpha_conjugate = sp.symbols("alpha_conjugate_0:3")
    cross = sp.zeros(3)
    cross_adjoint = sp.zeros(3)
    for index, (first, second) in enumerate(generators()):
        cross += alpha_conjugate[index] * (
            first / boundary + boundary * second.T
        ) / 2
        cross_adjoint += alpha[index] * (
            boundary * first.T + second / boundary
        ) / 2

    selected = sp.factor(
        (
            conjugate_support_vector.T
            * cross_adjoint
            * reduced_resolvent
            * cross
            * support_vector
        )[0]
    )
    orthogonal = sp.factor(
        (
            conjugate_support_vector.T
            * cross
            * reduced_resolvent
            * cross_adjoint
            * support_vector
        )[0]
    )
    average = sp.factor((selected + orthogonal) / 2)
    difference = sp.factor(selected - orthogonal)

    expected_average = (
        sp.Rational(5, 128) * alpha[0] * alpha_conjugate[0]
        + sp.Rational(1, 4) * alpha[1] * alpha_conjugate[1]
        + sp.Rational(5, 72) * alpha[2] * alpha_conjugate[2]
        + sp.Rational(1, 32)
        * (
            alpha[0] * alpha_conjugate[2] * boundary**2
            + alpha[2] * alpha_conjugate[0] / boundary**2
        )
    )
    coupling = (
        3 * alpha[0] * alpha_conjugate[1]
        + 4 * alpha[1] * alpha_conjugate[2]
    )
    expected_difference = root_two / 24 * (
        coupling * boundary
        + (
            3 * alpha[1] * alpha_conjugate[0]
            + 4 * alpha[2] * alpha_conjugate[1]
        )
        / boundary
    )
    if sp.factor(average - expected_average) != 0:
        raise AssertionError("unexpected mean second support coefficient")
    if sp.simplify(difference - expected_difference) != 0:
        raise AssertionError("unexpected split second support coefficient")

    if sp.simplify(
        reduced_resolvent * (sp.eye(3) - support_operator)
        * reduced_resolvent
        - reduced_resolvent
    ) != sp.zeros(3):
        raise AssertionError("the reduced support resolvent identity failed")

    print("PASS repeated p=3 second support reduction")
    print("effective support matrix is diagonal")
    print("top coefficient = mean + abs(difference)/2")
    print(f"mean = {expected_average}")
    print(f"difference = {expected_difference}")
    print("only Fourier modes 0 and 2 survive before the absolute-value split")
    print("the first Fourier coefficient is zero, so K(C3)=kappa_hat(0)*C3")


if __name__ == "__main__":
    main()
