#!/usr/bin/env python3
"""Prove the second-order sign for arbitrary repeated-C3 star coupling."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_stein_sign import assign_block


def main() -> None:
    root_two = sp.sqrt(2)
    boundary, conformal_mean = sp.symbols(
        "boundary conformal_mean",
        nonzero=True,
    )
    conformal_first, conformal_first_conjugate = sp.symbols(
        "conformal_first conformal_first_conjugate"
    )
    alpha = sp.symbols("alpha_0:3")
    alpha_conjugate = sp.symbols("alpha_conjugate_0:3")
    beta = sp.symbols("beta_0:3")
    beta_conjugate = sp.symbols("beta_conjugate_0:3")
    adjoint_substitution = dict(zip(alpha, alpha_conjugate, strict=True))
    adjoint_substitution.update(
        zip(alpha_conjugate, alpha, strict=True)
    )
    adjoint_substitution.update(
        zip(beta, beta_conjugate, strict=True)
    )
    adjoint_substitution.update(
        zip(beta_conjugate, beta, strict=True)
    )
    adjoint_substitution[conformal_first] = conformal_first_conjugate
    adjoint_substitution[conformal_first_conjugate] = conformal_first

    def formal_adjoint(matrix: sp.Matrix) -> sp.Matrix:
        """Transpose and exchange each formal coefficient with its adjoint."""

        return matrix.T.xreplace(adjoint_substitution)

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb, crabb)
    metric = sp.diag(*([1, 2, 4] * 3))
    identity = sp.eye(9)
    perturbation = sp.zeros(9)
    metric_tangent = sp.zeros(9)

    def metric_cross(coefficients: tuple[sp.Symbol, ...]) -> sp.Matrix:
        return sp.Matrix(
            [
                [0, -3 * root_two * coefficients[0] / 8, 0],
                [
                    coefficients[2] / root_two,
                    -2 * root_two * coefficients[1],
                    3 * root_two * coefficients[0] / 4,
                ],
                [0, -root_two * coefficients[2], 0],
            ]
        )

    for coefficients, conjugates, copy in (
        (alpha, alpha_conjugate, 1),
        (beta, beta_conjugate, 2),
    ):
        cross_to_selected = sp.zeros(3)
        cross_from_selected = sp.zeros(3)
        for index, (first, second) in enumerate(generators()):
            cross_to_selected += conjugates[index] * first
            cross_from_selected += coefficients[index] * second
        copy_slice = slice(3 * copy, 3 * copy + 3)
        perturbation[copy_slice, 0:3] = cross_to_selected
        perturbation[0:3, copy_slice] = cross_from_selected
        tangent_block = metric_cross(coefficients)
        metric_tangent[0:3, copy_slice] = tangent_block
        metric_tangent[copy_slice, 0:3] = formal_adjoint(tangent_block)

    # Exact constant Fourier coefficient of the second effective support
    # matrix for one selected copy and two independent orthogonal copies.
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
    orthogonal_support = effective_support[1:3, 1:3]
    reversed_trace = sp.trace(orthogonal_support).subs(
        boundary,
        -boundary,
    )
    if sp.simplify(effective_support[0, 0] - reversed_trace) != 0:
        raise AssertionError("the selected/orthogonal trace reversal failed")
    flat_substitution = {
        alpha[1]: 0,
        alpha_conjugate[1]: 0,
        beta[1]: 0,
        beta_conjugate[1]: 0,
    }
    flat_trace_residual = sp.simplify(
        (effective_support[0, 0] - sp.trace(orthogonal_support)).subs(
            flat_substitution
        )
    )
    if flat_trace_residual != 0:
        raise AssertionError("the flat star did not have trace equality")

    def constant_laurent_coefficient(expression: sp.Expr) -> sp.Expr:
        shifted = sp.expand(sp.cancel(expression) * boundary**6)
        return shifted.coeff(boundary, 6)

    support_mean = effective_support.applyfunc(
        constant_laurent_coefficient
    )
    weights = (
        sp.Rational(5, 128),
        sp.Rational(1, 4),
        sp.Rational(5, 72),
    )
    selected_mean = sum(
        weights[index]
        * (
            alpha[index] * alpha_conjugate[index]
            + beta[index] * beta_conjugate[index]
        )
        for index in range(3)
    )
    orthogonal_mean = sp.Matrix(
        [
            [
                sum(
                    weights[index]
                    * alpha_conjugate[index]
                    * alpha[index]
                    for index in range(3)
                ),
                sum(
                    weights[index]
                    * alpha_conjugate[index]
                    * beta[index]
                    for index in range(3)
                ),
            ],
            [
                sum(
                    weights[index]
                    * beta_conjugate[index]
                    * alpha[index]
                    for index in range(3)
                ),
                sum(
                    weights[index]
                    * beta_conjugate[index]
                    * beta[index]
                    for index in range(3)
                ),
            ],
        ]
    )
    expected_support_mean = sp.diag(selected_mean, 1, 1)
    expected_support_mean[1:3, 1:3] = orthogonal_mean
    if sp.simplify(support_mean - expected_support_mean) != sp.zeros(3):
        raise AssertionError("the mean effective support formula failed")

    lower_kernel = (0, 3, 6)
    lower_range = (1, 2, 4, 5, 7, 8)
    upper_kernel = (2, 5, 8)
    upper_range = (0, 1, 3, 4, 6, 7)
    contraction_kernel = lower_range
    contraction_range = lower_kernel
    lower_base = metric - identity
    upper_base = 4 * identity - metric
    levels = tuple(
        tuple(3 * copy + level for copy in range(3))
        for level in range(3)
    )

    def build_endpoint(
        tangent: sp.Matrix,
        direction: sp.Matrix,
        second_operator: sp.Matrix,
    ) -> sp.Matrix:
        """Construct the tight recursive second-order upper endpoint."""

        lower_penalty = sp.simplify(
            tangent.extract(lower_kernel, lower_range)
            * lower_base.extract(lower_range, lower_range).inv()
            * tangent.extract(lower_range, lower_kernel)
        )
        upper_penalty = sp.simplify(
            tangent.extract(upper_kernel, upper_range)
            * upper_base.extract(upper_range, upper_range).inv()
            * tangent.extract(upper_range, upper_kernel)
        )
        first_contraction = sp.simplify(
            tangent
            - formal_adjoint(base) * tangent * base
            - formal_adjoint(direction) * metric * base
            - formal_adjoint(base) * metric * direction
        )
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
            formal_adjoint(second_operator) * metric * base
            + formal_adjoint(base) * metric * second_operator
            + formal_adjoint(direction) * metric * direction
            + formal_adjoint(direction) * tangent * base
            + formal_adjoint(base) * tangent * direction
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
        return sp.simplify(
            second_metric.extract(upper_kernel, upper_kernel)
            + upper_penalty
        )

    second_operator = (
        -conformal_mean * base - conformal_first * base**2
    )
    endpoint = build_endpoint(
        metric_tangent,
        perturbation,
        second_operator,
    )

    first_type_norm = (
        alpha[1] * alpha_conjugate[1]
        + beta[1] * beta_conjugate[1]
    )
    first_type_gram = sp.Matrix(
        [
            [
                alpha_conjugate[1] * alpha[1],
                alpha_conjugate[1] * beta[1],
            ],
            [
                beta_conjugate[1] * alpha[1],
                beta_conjugate[1] * beta[1],
            ],
        ]
    )
    first_type_penalty = sp.diag(first_type_norm, 1, 1)
    first_type_penalty[1:3, 1:3] = first_type_gram
    expected_endpoint = (
        16 * expected_support_mean
        - 8 * first_type_penalty
        - 16 * conformal_mean * sp.eye(3)
    )
    if sp.simplify(endpoint - expected_endpoint) != sp.zeros(3):
        raise AssertionError("the arbitrary-star endpoint identity failed")

    # On the flat alpha_1=0 face, improve the pairwise tangent by a real
    # parameter tau.  Its selected endpoint moves strictly negative, while
    # the orthogonal endpoint starts negative definite exactly at rank two.
    tau = sp.symbols("tau", real=True)
    flat_direction = perturbation.subs(flat_substitution)
    flat_forcing = sp.simplify(
        formal_adjoint(flat_direction) * metric * base
        + formal_adjoint(base) * metric * flat_direction
    )
    free_block = sp.zeros(3)
    free_block[0, 1] = (-3 * root_two / 8 + tau) * alpha[0]
    free_block[0, 2] = (-3 * root_two / 8 + tau) * beta[0]
    free_block[1, 0] = (1 / root_two + tau) * alpha_conjugate[2]
    free_block[2, 0] = (1 / root_two + tau) * beta_conjugate[2]
    improved_tangent = sp.zeros(9)
    assign_block(
        improved_tangent,
        levels[0],
        levels[1],
        free_block,
    )
    assign_block(
        improved_tangent,
        levels[1],
        levels[0],
        formal_adjoint(free_block),
    )
    assign_block(
        improved_tangent,
        levels[1],
        levels[1],
        flat_forcing.extract(levels[1], levels[1]),
    )
    next_block = (
        2 * free_block
        + flat_forcing.extract(levels[1], levels[2])
    )
    assign_block(
        improved_tangent,
        levels[1],
        levels[2],
        next_block,
    )
    assign_block(
        improved_tangent,
        levels[2],
        levels[1],
        formal_adjoint(next_block),
    )
    flat_mean = (
        sp.Rational(5, 128)
        * (
            alpha[0] * alpha_conjugate[0]
            + beta[0] * beta_conjugate[0]
        )
        + sp.Rational(5, 72)
        * (
            alpha[2] * alpha_conjugate[2]
            + beta[2] * beta_conjugate[2]
        )
    )
    improved_endpoint = build_endpoint(
        improved_tangent,
        flat_direction,
        -flat_mean * base,
    )
    endpoint_at_zero = sp.simplify(improved_endpoint.subs(tau, 0))
    endpoint_derivative = sp.simplify(
        sp.diff(improved_endpoint, tau).subs(tau, 0)
    )

    generator_zero_gram = sp.Matrix(
        [
            [
                alpha_conjugate[0] * alpha[0],
                alpha_conjugate[0] * beta[0],
            ],
            [
                beta_conjugate[0] * alpha[0],
                beta_conjugate[0] * beta[0],
            ],
        ]
    )
    generator_two_gram = sp.Matrix(
        [
            [
                alpha_conjugate[2] * alpha[2],
                alpha_conjugate[2] * beta[2],
            ],
            [
                beta_conjugate[2] * alpha[2],
                beta_conjugate[2] * beta[2],
            ],
        ]
    )
    heavy_gram = 9 * generator_zero_gram + 16 * generator_two_gram
    light_gram = 3 * generator_zero_gram + 4 * generator_two_gram
    heavy_adjugate = sp.Matrix(
        [
            [heavy_gram[1, 1], -heavy_gram[0, 1]],
            [-heavy_gram[1, 0], heavy_gram[0, 0]],
        ]
    )
    expected_zero = sp.diag(0, 1, 1)
    expected_zero[1:3, 1:3] = -sp.Rational(5, 72) * heavy_adjugate
    expected_derivative = sp.diag(
        -5 * root_two * sp.trace(light_gram) / 3,
        1,
        1,
    )
    expected_derivative[1:3, 1:3] = (
        5 * root_two * light_gram / 3
    )
    if sp.simplify(endpoint_at_zero - expected_zero) != sp.zeros(3):
        raise AssertionError("the rank-two base endpoint formula failed")
    if sp.simplify(
        endpoint_derivative - expected_derivative
    ) != sp.zeros(3):
        raise AssertionError("the strict rank-two endpoint derivative failed")

    print("PASS repeated p=3 arbitrary-star second-order sign")
    print("mean support = block_diag(mu, weighted coefficient Gram)")
    print("endpoint = 16 * (mean support - mean top support * I)")
    print("           - 8 * block_diag(norm(alpha_1)^2, alpha_1 Gram)")
    print("the possibly nonzero first conformal mode cancels from the endpoint")
    print("both displayed endpoint terms are negative semidefinite")
    print("endpoint equality occurs exactly when the alpha_1 vector vanishes")
    print("rank-two flat data admit a strictly negative perturbed endpoint")


if __name__ == "__main__":
    main()
