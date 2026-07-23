#!/usr/bin/env python3
"""Prove strict cross descent under a repeated-C3 radial copy gap."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators


def main() -> None:
    root_two = sp.sqrt(2)
    alpha = sp.symbols("alpha_0:3")
    alpha_conjugate = sp.symbols("alpha_conjugate_0:3")
    gap, free_zero, free_two = sp.symbols(
        "gap free_zero free_two",
        positive=True,
    )
    adjoint_substitution = dict(
        zip(alpha, alpha_conjugate, strict=True)
    )
    adjoint_substitution.update(
        zip(alpha_conjugate, alpha, strict=True)
    )

    def formal_adjoint(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.T.xreplace(adjoint_substitution)

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb)
    metric = sp.diag(1, 2, 4, 1, 2, 4)
    identity = sp.eye(6)

    cross_to_selected = sp.zeros(3)
    cross_from_selected = sp.zeros(3)
    for index, (first, second) in enumerate(generators()):
        cross_to_selected += alpha_conjugate[index] * first
        cross_from_selected += alpha[index] * second
    cross_direction = sp.zeros(6)
    cross_direction[3:6, 0:3] = cross_to_selected
    cross_direction[0:3, 3:6] = cross_from_selected

    # The physical first perturbation grows the selected Crabb copy and
    # shrinks the other one by the same positive radial amount.
    radial_direction = sp.diag(gap * crabb, -gap * crabb)
    physical_direction = radial_direction + cross_direction

    # The numerical-range radius is 1 + gap*epsilon at first order.  After
    # pullback the selected diagonal tangent vanishes and the other diagonal
    # tangent is -2*gap*C3.
    pulled_first = physical_direction - gap * base

    support_mean = (
        sp.Rational(5, 128) * alpha[0] * alpha_conjugate[0]
        + sp.Rational(1, 4) * alpha[1] * alpha_conjugate[1]
        + sp.Rational(5, 72) * alpha[2] * alpha_conjugate[2]
    )
    support_coupling = (
        3 * alpha[0] * alpha_conjugate[1]
        + 4 * alpha[1] * alpha_conjugate[2]
    )
    first_fourier_operator = root_two * support_coupling * crabb**2 / 24
    support_map_at_crabb = support_mean * crabb + first_fourier_operator
    pulled_second = sp.zeros(6)
    pulled_second[0:3, 0:3] = -support_map_at_crabb
    pulled_second[3:6, 3:6] = (
        2 * gap**2 * crabb - support_map_at_crabb
    )
    pulled_second[3:6, 0:3] = -gap * cross_to_selected
    pulled_second[0:3, 3:6] = -gap * cross_from_selected

    # The two free cross entries are the complete level-0/1 tangent freedom.
    # Active cross Stein equations force the remaining entries.
    cross_metric = sp.Matrix(
        [
            [0, free_zero * alpha[0], 0],
            [
                free_two * alpha[2],
                -2 * root_two * alpha[1],
                (2 * free_zero + 3 * root_two / 2) * alpha[0],
            ],
            [0, (2 * free_two - 2 * root_two) * alpha[2], 0],
        ]
    )
    metric_tangent = sp.zeros(6)
    metric_tangent[0:3, 3:6] = cross_metric
    metric_tangent[3:6, 0:3] = formal_adjoint(cross_metric)
    metric_tangent[3:6, 3:6] = sp.diag(gap, 0, -gap)

    first_contraction = sp.simplify(
        metric_tangent
        - formal_adjoint(base) * metric_tangent * base
        - formal_adjoint(pulled_first) * metric * base
        - formal_adjoint(base) * metric * pulled_first
    )
    selected_contraction = (1, 2)
    orthogonal_contraction = (4, 5)
    if first_contraction.extract(
        selected_contraction,
        selected_contraction,
    ) != sp.zeros(2):
        raise AssertionError("the selected first contraction face moved")
    if first_contraction.extract(
        selected_contraction,
        orthogonal_contraction,
    ) != sp.zeros(2):
        raise AssertionError("the first cross contraction did not vanish")
    expected_orthogonal_slack = sp.diag(6 * gap, 15 * gap)
    if sp.simplify(
        first_contraction.extract(
            orthogonal_contraction,
            orthogonal_contraction,
        )
        - expected_orthogonal_slack
    ) != sp.zeros(2):
        raise AssertionError("the orthogonal first contraction was not strict")

    # Only the selected copy remains active at second order.  Eliminate its
    # lower, Stein, and upper endpoints exactly.
    lower_kernel = (0,)
    lower_range = (1, 2, 4, 5)
    upper_kernel = (2,)
    upper_range = (0, 1, 3, 4)
    contraction_kernel = (1, 2)
    contraction_range = (0, 3)
    lower_base = metric - identity
    upper_base = 4 * identity - metric
    lower_penalty = sp.simplify(
        metric_tangent.extract(lower_kernel, lower_range)
        * lower_base.extract(lower_range, lower_range).inv()
        * metric_tangent.extract(lower_range, lower_kernel)
    )[0]
    upper_penalty = sp.simplify(
        metric_tangent.extract(upper_kernel, upper_range)
        * upper_base.extract(upper_range, upper_range).inv()
        * metric_tangent.extract(upper_range, upper_kernel)
    )[0]
    contraction_penalty = sp.simplify(
        first_contraction.extract(
            contraction_kernel,
            contraction_range,
        )
        * first_contraction.extract(
            contraction_range,
            contraction_kernel,
        )
    )
    forcing = sp.simplify(
        formal_adjoint(pulled_second) * metric * base
        + formal_adjoint(base) * metric * pulled_second
        + formal_adjoint(pulled_first) * metric * pulled_first
        + formal_adjoint(pulled_first) * metric_tangent * base
        + formal_adjoint(base) * metric_tangent * pulled_first
    )
    second_level_one = sp.simplify(
        2 * lower_penalty
        + forcing[1, 1]
        + contraction_penalty[0, 0]
    )
    second_level_two = sp.simplify(
        2 * second_level_one
        + forcing[2, 2]
        + contraction_penalty[1, 1]
    )
    endpoint = sp.factor(second_level_two + upper_penalty)
    expected_endpoint = (
        (
            96 * free_zero**2
            - 48 * root_two * free_zero
            - 63
        )
        * alpha[0]
        * alpha_conjugate[0]
        / 24
        - 8 * alpha[1] * alpha_conjugate[1]
        + (
            96 * free_two**2
            - 256 * root_two * free_two
            + 208
        )
        * alpha[2]
        * alpha_conjugate[2]
        / 24
    )
    if sp.simplify(endpoint - expected_endpoint) != 0:
        raise AssertionError("the radial-gap endpoint formula failed")

    optimum = sp.simplify(
        endpoint.subs(
            {
                free_zero: root_two / 4,
                free_two: 4 * root_two / 3,
            }
        )
    )
    expected_optimum = (
        -sp.Rational(25, 8) * alpha[0] * alpha_conjugate[0]
        - 8 * alpha[1] * alpha_conjugate[1]
        - sp.Rational(50, 9) * alpha[2] * alpha_conjugate[2]
    )
    if sp.simplify(optimum - expected_optimum) != 0:
        raise AssertionError("the optimized radial-gap sign failed")

    print("PASS repeated p=3 radial-gap cross theorem")
    print("orthogonal first contraction slack = diag(6*gap, 15*gap)")
    print("optimized selected second endpoint =")
    print("  -25*abs(alpha_0)^2/8 - 8*abs(alpha_1)^2")
    print("  -50*abs(alpha_2)^2/9")


if __name__ == "__main__":
    main()
