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
    contraction_penalty = sp.simplify(
        first_contraction.extract(contraction_kernel, contraction_range)
        * first_contraction.extract(contraction_range, contraction_kernel)
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
    levels = tuple(
        tuple(3 * copy + level for copy in range(3))
        for level in range(3)
    )
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

    print("PASS repeated p=3 arbitrary-star second-order sign")
    print("mean support = block_diag(mu, weighted coefficient Gram)")
    print("endpoint = 16 * (mean support - mean top support * I)")
    print("           - 8 * block_diag(norm(alpha_1)^2, alpha_1 Gram)")
    print("the possibly nonzero first conformal mode cancels from the endpoint")
    print("both displayed endpoint terms are negative semidefinite")


if __name__ == "__main__":
    main()
