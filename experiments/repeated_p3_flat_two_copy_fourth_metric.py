#!/usr/bin/env python3
"""Prove cancellation of the transverse fourth metric endpoint at L100."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_third_metric import tight_fourth_endpoint


def assert_zero_matrix(matrix: sp.Matrix, message: str) -> None:
    """Prove a small rational-function matrix is zero entry by entry."""

    for entry in matrix:
        if sp.factor(sp.cancel(entry)) != 0:
            raise AssertionError(message)


def polynomial_derivative(
    base: sp.Matrix,
    direction: sp.Matrix,
    power: int,
) -> sp.Matrix:
    """Return the Fréchet derivative of ``base**power``."""

    return sum(
        (
            base**left_power
            * direction
            * base ** (power - 1 - left_power)
            for left_power in range(power)
        ),
        sp.zeros(*base.shape),
    )


def cubic_second_derivative(
    base: sp.Matrix,
    first: sp.Matrix,
    second: sp.Matrix,
) -> sp.Matrix:
    """Return the bilinear second derivative of the matrix cube."""

    derivative = sp.zeros(*base.shape)
    factors = (base, base, base)
    for first_position in range(3):
        for second_position in range(3):
            if first_position == second_position:
                continue
            word = list(factors)
            word[first_position] = first
            word[second_position] = second
            derivative += word[0] * word[1] * word[2]
    return derivative


def flat_direction(copy: sp.Matrix) -> sp.Matrix:
    """Return L88's physical flat perturbation for a copy matrix."""

    generator_first, generator_second = generators()[0]
    return (
        sp.kronecker_product(
            copy.conjugate().T,
            generator_first,
        )
        + sp.kronecker_product(copy, generator_second)
    )


def loop_metric(value: sp.Expr) -> sp.Matrix:
    """Return the canonical first metric for one flat diagonal loop."""

    root_two = sp.sqrt(2)
    upper = sp.Matrix(
        [
            [0, -3 * root_two * value / 8, 0],
            [0, 0, 3 * root_two * value / 4],
            [0, 0, 0],
        ]
    )
    return upper + upper.conjugate().T


def main() -> None:
    root_two = sp.sqrt(2)
    edge, diagonal = sp.symbols(
        "edge diagonal",
        real=True,
    )
    cubic_zero, cubic_two, cubic_four = sp.symbols(
        "cubic_zero cubic_two cubic_four",
        real=True,
    )
    fourth_zero, fourth_one = sp.symbols(
        "fourth_zero fourth_one",
        real=True,
    )
    free_zero_zero, free_zero_one = sp.symbols(
        "free_zero_zero free_zero_one",
        real=True,
    )
    free_one_zero, free_one_one = sp.symbols(
        "free_one_zero free_one_one",
        real=True,
    )

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb)
    metric = sp.diag(1, 2, 4, 1, 2, 4)
    normal_copy = sp.diag(diagonal, -diagonal)
    edge_copy = sp.Matrix([[0, 1], [0, 0]])
    normal = flat_direction(normal_copy)
    transverse = flat_direction(edge_copy)
    first_operator = normal + edge * transverse

    center = 3 * diagonal**2 / (8 * root_two)
    common_block = center * sp.Matrix(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    physical_second = sp.diag(common_block, common_block)

    second_mean = 5 * diagonal**2 / 64
    second_cubic = 9 * diagonal**2 / 64
    second_operator = physical_second - second_mean * base
    normal_third_operator = (
        -second_mean * normal
        - second_cubic * polynomial_derivative(base, normal, 3)
    )
    transverse_third_operator = (
        -second_mean * transverse
        - second_cubic
        * polynomial_derivative(base, transverse, 3)
        - cubic_zero * base
    )
    third_operator = (
        normal_third_operator
        + edge * transverse_third_operator
    )

    normal_fourth_operator = (
        -second_mean * physical_second
        - second_cubic
        * polynomial_derivative(base, physical_second, 3)
        - second_cubic
        * cubic_second_derivative(base, normal, normal)
        / 2
        + second_mean**2 * base
        - fourth_zero * base
        - fourth_one * base**2
    )
    cubic_response_derivative = (
        cubic_zero * normal
        + 2
        * cubic_two
        * polynomial_derivative(base, normal, 3)
        + 2
        * cubic_four
        * polynomial_derivative(base, normal, 5)
    )
    transverse_fourth_operator = (
        -second_cubic
        * cubic_second_derivative(base, normal, transverse)
        - cubic_response_derivative
    )
    fourth_operator = (
        normal_fourth_operator
        + edge * transverse_fourth_operator
    )
    operators = (
        base,
        first_operator,
        second_operator,
        third_operator,
        fourth_operator,
    )

    edge_metric = sp.Matrix(
        [
            [0, -3 * root_two / 8, 0],
            [0, 0, 3 * root_two / 4],
            [0, 0, 0],
        ]
    )
    metric_tangent = sp.zeros(6)
    metric_tangent[0:3, 0:3] = loop_metric(diagonal)
    metric_tangent[3:6, 3:6] = -loop_metric(diagonal)
    metric_tangent[0:3, 3:6] = edge * edge_metric
    metric_tangent[3:6, 0:3] = edge * edge_metric.T

    zero_free_endpoint = tight_fourth_endpoint(
        operators,
        metric,
        metric_tangent,
    )
    zero_free_derivative = zero_free_endpoint.diff(edge).subs(
        edge,
        0,
    )
    expected_zero_free_derivative = (
        25
        * diagonal**3
        / 8
        * sp.Matrix([[0, 1], [1, 0]])
    )
    assert_zero_matrix(
        zero_free_derivative - expected_zero_free_derivative,
        "the uncancelled transverse fourth endpoint changed",
    )

    general_free_block = edge * sp.Matrix(
        [
            [free_zero_zero, free_zero_one],
            [free_one_zero, free_one_one],
        ]
    )
    general_free_endpoint = tight_fourth_endpoint(
        operators,
        metric,
        metric_tangent,
        third_free_block=general_free_block,
    )
    general_free_derivative = general_free_endpoint.diff(edge).subs(
        edge,
        0,
    )
    expected_general_derivative = (
        5
        * diagonal
        / 8
        * (
            5 * diagonal**2
            + 8
            * root_two
            * (free_zero_one - free_one_zero)
        )
        * sp.Matrix([[0, 1], [1, 0]])
    )
    assert_zero_matrix(
        general_free_derivative - expected_general_derivative,
        "the general third-metric freedom formula failed",
    )

    cancellation_scale = 5 * diagonal**2 / (16 * root_two)
    cancelling_free_block = (
        edge
        * cancellation_scale
        * sp.Matrix([[0, -1], [1, 0]])
    )
    cancelled_endpoint = tight_fourth_endpoint(
        operators,
        metric,
        metric_tangent,
        third_free_block=cancelling_free_block,
    )
    assert_zero_matrix(
        cancelled_endpoint.diff(edge).subs(edge, 0),
        "the transverse fourth metric endpoint did not cancel",
    )
    normal_endpoint = cancelled_endpoint.subs(edge, 0)
    assert_zero_matrix(
        normal_endpoint
        - sp.trace(normal_endpoint) * sp.eye(2) / 2,
        "the normal fourth metric endpoint was not scalar",
    )

    print("PASS repeated p=3 weighted normal fourth metric cancellation")
    print("without the free block: E4'=(25*d^3/8)*[[0,1],[1,0]]")
    print("the third-metric skew block cancels the full transverse derivative")
    print("all abstract third/fourth conformal coefficients cancel")


if __name__ == "__main__":
    main()
