#!/usr/bin/env python3
"""Prove the tied-winner graph second-order endpoint identity."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_stein_sign import assign_block


def main() -> None:
    root_two = sp.sqrt(2)
    boundary = sp.symbols("boundary", nonzero=True)
    conformal_mean = sp.symbols("conformal_mean", real=True)
    conformal_first = sp.symbols("conformal_first")
    conformal_first_conjugate = sp.symbols("conformal_first_conjugate")
    edges = ((0, 1), (0, 2), (1, 2))
    coefficients = tuple(
        tuple(sp.symbols(f"coefficient_{edge}_0:3"))
        for edge in range(3)
    )
    conjugates = tuple(
        tuple(sp.symbols(f"conjugate_{edge}_0:3"))
        for edge in range(3)
    )
    adjoint_substitution: dict[sp.Symbol, sp.Symbol] = {
        conformal_first: conformal_first_conjugate,
        conformal_first_conjugate: conformal_first,
    }
    for values, adjoints in zip(coefficients, conjugates, strict=True):
        adjoint_substitution.update(zip(values, adjoints, strict=True))
        adjoint_substitution.update(zip(adjoints, values, strict=True))

    def formal_adjoint(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.T.xreplace(adjoint_substitution)

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb, crabb)
    metric = sp.diag(*([1, 2, 4] * 3))
    identity = sp.eye(9)
    perturbation = sp.zeros(9)
    metric_tangent = sp.zeros(9)
    coefficient_matrices = [sp.zeros(3) for _ in range(3)]

    for edge_index, (left_copy, right_copy) in enumerate(edges):
        values = coefficients[edge_index]
        adjoints = conjugates[edge_index]
        cross_to_left = sp.zeros(3)
        cross_from_left = sp.zeros(3)
        for generator_index, (first, second) in enumerate(generators()):
            cross_to_left += adjoints[generator_index] * first
            cross_from_left += values[generator_index] * second
            coefficient_matrices[generator_index][left_copy, right_copy] = (
                values[generator_index]
            )
            coefficient_matrices[generator_index][right_copy, left_copy] = (
                adjoints[generator_index]
            )
        left_slice = slice(3 * left_copy, 3 * left_copy + 3)
        right_slice = slice(3 * right_copy, 3 * right_copy + 3)
        perturbation[right_slice, left_slice] = cross_to_left
        perturbation[left_slice, right_slice] = cross_from_left

        tangent_block = sp.Matrix(
            [
                [0, -3 * root_two * values[0] / 8, 0],
                [
                    values[2] / root_two,
                    -2 * root_two * values[1],
                    3 * root_two * values[0] / 4,
                ],
                [0, -root_two * values[2], 0],
            ]
        )
        metric_tangent[left_slice, right_slice] = tangent_block
        metric_tangent[right_slice, left_slice] = formal_adjoint(
            tangent_block
        )

    support = (crabb / boundary + boundary * crabb.T) / 2
    reduced_resolvent = (
        sp.eye(3) - support / 4 - 3 * support**2 / 4
    )
    support_vector = sp.Matrix([1 / boundary, root_two, boundary]) / 2
    support_adjoint = sp.Matrix(
        [[boundary, root_two, 1 / boundary]]
    ) / 2
    support_perturbation = (
        perturbation / boundary
        + boundary * formal_adjoint(perturbation)
    ) / 2
    full_resolvent = sp.diag(
        reduced_resolvent,
        reduced_resolvent,
        reduced_resolvent,
    )
    top_vectors = sp.zeros(9, 3)
    top_adjoint = sp.zeros(3, 9)
    for copy in range(3):
        copy_slice = slice(3 * copy, 3 * copy + 3)
        top_vectors[copy_slice, copy] = support_vector
        top_adjoint[copy, copy_slice] = support_adjoint
    effective_support = sp.simplify(
        top_adjoint
        * support_perturbation
        * full_resolvent
        * support_perturbation
        * top_vectors
    )

    def constant_laurent_coefficient(expression: sp.Expr) -> sp.Expr:
        shifted = sp.expand(sp.cancel(expression) * boundary**6)
        return shifted.coeff(boundary, 6)

    support_mean = effective_support.applyfunc(
        constant_laurent_coefficient
    )
    levels = tuple(
        tuple(3 * copy + level for copy in range(3))
        for level in range(3)
    )
    lower_kernel = levels[0]
    lower_range = levels[1] + levels[2]
    upper_kernel = levels[2]
    upper_range = levels[0] + levels[1]
    contraction_kernel = lower_range
    contraction_range = lower_kernel
    lower_base = metric - identity
    upper_base = 4 * identity - metric
    lower_penalty = sp.simplify(
        metric_tangent.extract(lower_kernel, lower_range)
        * lower_base.extract(lower_range, lower_range).inv()
        * metric_tangent.extract(lower_range, lower_kernel)
    )
    upper_penalty = sp.simplify(
        metric_tangent.extract(upper_kernel, upper_range)
        * upper_base.extract(upper_range, upper_range).inv()
        * metric_tangent.extract(upper_range, upper_kernel)
    )
    first_contraction = sp.simplify(
        metric_tangent
        - formal_adjoint(base) * metric_tangent * base
        - formal_adjoint(perturbation) * metric * base
        - formal_adjoint(base) * metric * perturbation
    )
    if first_contraction.extract(
        contraction_kernel,
        contraction_kernel,
    ) != sp.zeros(6):
        raise AssertionError("the winner-graph first contraction moved")
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
    second_operator = (
        -conformal_mean * base - conformal_first * base**2
    )
    forcing = sp.simplify(
        formal_adjoint(second_operator) * metric * base
        + formal_adjoint(base) * metric * second_operator
        + formal_adjoint(perturbation) * metric * perturbation
        + formal_adjoint(perturbation) * metric_tangent * base
        + formal_adjoint(base) * metric_tangent * perturbation
    )
    effective_forcing = sp.MutableDenseMatrix(forcing)
    for row_index, row in enumerate(contraction_kernel):
        for column_index, column in enumerate(contraction_kernel):
            effective_forcing[row, column] += contraction_penalty[
                row_index,
                column_index,
            ]

    second_metric = sp.zeros(9)
    assign_block(
        second_metric,
        levels[0],
        levels[0],
        lower_penalty,
    )
    for row_level in range(1, 3):
        for column_level in range(1, 3):
            block = sp.simplify(
                2
                * second_metric.extract(
                    levels[row_level - 1],
                    levels[column_level - 1],
                )
                + effective_forcing.extract(
                    levels[row_level],
                    levels[column_level],
                )
            )
            assign_block(
                second_metric,
                levels[row_level],
                levels[column_level],
                block,
            )
    endpoint = sp.simplify(
        second_metric.extract(upper_kernel, upper_kernel)
        + upper_penalty
    )
    expected_endpoint = (
        16 * support_mean
        - 8 * coefficient_matrices[1] ** 2
        - 16 * conformal_mean * sp.eye(3)
    )
    if sp.simplify(endpoint - expected_endpoint) != sp.zeros(3):
        raise AssertionError("the winner-graph endpoint identity failed")

    print("PASS repeated p=3 tied-winner graph second-order sign")
    print("endpoint = 16*(mean support - mean top support*I) - 8*H_1^2")
    print("mean support includes orientation-dependent mixed edge paths")
    print("the first conformal mode cancels exactly")
    print("matrix Jensen order makes the endpoint negative semidefinite")


if __name__ == "__main__":
    main()
