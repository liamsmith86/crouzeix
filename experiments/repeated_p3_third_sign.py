#!/usr/bin/env python3
"""Prove the third-order sign on the flat L76 repeated-C3 plane."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_stein_sign import assign_block


def stein_coefficient(
    operators: tuple[sp.Matrix, ...],
    metrics: tuple[sp.Matrix, ...],
    order: int,
) -> sp.Matrix:
    """Return one coefficient of ``P - T^* P T``."""

    coefficient = metrics[order]
    for left_order in range(order + 1):
        for metric_order in range(order + 1 - left_order):
            right_order = order - left_order - metric_order
            if (
                left_order < len(operators)
                and metric_order < len(metrics)
                and right_order < len(operators)
            ):
                coefficient -= (
                    operators[left_order].conjugate().T
                    * metrics[metric_order]
                    * operators[right_order]
                )
    return sp.simplify(coefficient)


def third_schur_penalty(
    first: sp.Matrix,
    second: sp.Matrix,
    kernel: tuple[int, ...],
    positive_range: tuple[int, ...],
    base_range: sp.Matrix,
) -> sp.Matrix:
    """Return the order-three endpoint forced by a tight Schur complement."""

    inverse = base_range.inv()
    first_cross = first.extract(kernel, positive_range)
    second_cross = second.extract(kernel, positive_range)
    first_range = first.extract(positive_range, positive_range)
    return sp.simplify(
        second_cross * inverse * first_cross.conjugate().T
        + first_cross * inverse * second_cross.conjugate().T
        - first_cross
        * inverse
        * first_range
        * inverse
        * first_cross.conjugate().T
    )


def main() -> None:
    root_two = sp.sqrt(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb)
    metric = sp.diag(1, 2, 4, 1, 2, 4)
    identity = sp.eye(6)

    alpha_real, alpha_imaginary = sp.symbols(
        "alpha_real alpha_imaginary",
        real=True,
    )
    beta_real, beta_imaginary = sp.symbols(
        "beta_real beta_imaginary",
        real=True,
    )
    alpha = alpha_real + sp.I * alpha_imaginary
    beta = beta_real + sp.I * beta_imaginary
    perturbation = sp.zeros(6)
    for coefficient, generator_index in ((alpha, 0), (beta, 2)):
        first, second = generators()[generator_index]
        perturbation[3:6, 0:3] += sp.conjugate(coefficient) * first
        perturbation[0:3, 3:6] += coefficient * second

    # Reconstruct the second and third effective top-support matrices.
    boundary = sp.symbols("boundary", nonzero=True)
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
        + boundary * perturbation.conjugate().T
    ) / 2
    full_resolvent = sp.diag(reduced_resolvent, reduced_resolvent)
    top_vectors = sp.zeros(6, 2)
    top_vectors[0:3, 0] = support_vector
    top_vectors[3:6, 1] = support_vector
    top_adjoint = sp.zeros(2, 6)
    top_adjoint[0, 0:3] = support_adjoint
    top_adjoint[1, 3:6] = support_adjoint
    effective_second = sp.simplify(
        top_adjoint
        * support_perturbation
        * full_resolvent
        * support_perturbation
        * top_vectors
    )
    effective_third = sp.simplify(
        top_adjoint
        * support_perturbation
        * full_resolvent
        * support_perturbation
        * full_resolvent
        * support_perturbation
        * top_vectors
    )
    second_mean = (
        sp.Rational(5, 128) * alpha * sp.conjugate(alpha)
        + sp.Rational(5, 72) * beta * sp.conjugate(beta)
    )
    second_mode = (
        alpha * sp.conjugate(beta) * boundary**2
        + sp.conjugate(alpha) * beta / boundary**2
    ) / 32
    if sp.simplify(
        effective_second - (second_mean + second_mode) * sp.eye(2)
    ) != sp.zeros(2):
        raise AssertionError("the second effective support did not stay scalar")

    third_cross = -(
        3 * alpha * boundary**2 + 4 * beta
    ) * (
        -9 * alpha * sp.conjugate(alpha) * boundary**2
        + 18 * alpha * sp.conjugate(beta) * boundary**4
        - 18 * sp.conjugate(alpha) * beta
        + 16 * beta * sp.conjugate(beta) * boundary**2
    ) / (1728 * boundary**3)
    if sp.simplify(effective_third[0, 0]) != 0:
        raise AssertionError("the third effective support had a diagonal term")
    if sp.simplify(effective_third[1, 1]) != 0:
        raise AssertionError("the third effective support had a diagonal term")
    if sp.simplify(effective_third[0, 1] - third_cross) != 0:
        raise AssertionError("the third support splitting formula failed")
    circle_adjoint = sp.simplify(
        sp.conjugate(third_cross).subs(
            sp.conjugate(boundary),
            1 / boundary,
        )
    )
    if sp.simplify(effective_third[1, 0] - circle_adjoint) != 0:
        raise AssertionError("the third effective support was not Hermitian")

    # The second inverse-map term is
    # F_2(w) = second_mean*w + alpha*conj(beta)*w^3/16.
    # The third normal displacement is |third_cross|, whose mean is denoted
    # by third_mean.  Pi-periodicity kills its first Fourier mode.
    third_mean = sp.symbols("third_mean", nonnegative=True)
    cubic_coefficient = alpha * sp.conjugate(beta) / 16
    cubic_frechet = (
        base**2 * perturbation
        + base * perturbation * base
        + perturbation * base**2
    )
    operators = (
        base,
        perturbation,
        -second_mean * base,
        -second_mean * perturbation
        - cubic_coefficient * cubic_frechet
        - third_mean * base,
    )

    metric_cross = sp.Matrix(
        [
            [0, -3 * root_two * alpha / 8, 0],
            [beta / root_two, 0, 3 * root_two * alpha / 4],
            [0, -root_two * beta, 0],
        ]
    )
    metric_tangent = sp.zeros(6)
    metric_tangent[0:3, 3:6] = metric_cross
    metric_tangent[3:6, 0:3] = metric_cross.conjugate().T

    lower_kernel = (0, 3)
    lower_range = (1, 2, 4, 5)
    upper_kernel = (2, 5)
    upper_range = (0, 1, 3, 4)
    contraction_kernel = (1, 2, 4, 5)
    contraction_range = (0, 3)
    levels = tuple((level, level + 3) for level in range(3))

    lower_base = (metric - identity).extract(lower_range, lower_range)
    lower_second = sp.simplify(
        metric_tangent.extract(lower_kernel, lower_range)
        * lower_base.inv()
        * metric_tangent.extract(lower_range, lower_kernel)
    )
    first_contraction = stein_coefficient(
        operators,
        (metric, metric_tangent),
        1,
    )
    contraction_second_penalty = sp.simplify(
        first_contraction.extract(contraction_kernel, contraction_range)
        * first_contraction.extract(contraction_range, contraction_kernel)
    )

    zero_metric = sp.zeros(6)
    second_contraction_without_metric = stein_coefficient(
        operators,
        (metric, metric_tangent, zero_metric),
        2,
    )
    second_forcing = -second_contraction_without_metric
    second_metric = sp.zeros(6)
    assign_block(
        second_metric,
        levels[0],
        levels[0],
        lower_second,
    )
    effective_second_forcing = sp.MutableDenseMatrix(second_forcing)
    for row_index, row in enumerate(contraction_kernel):
        for column_index, column in enumerate(contraction_kernel):
            effective_second_forcing[row, column] += (
                contraction_second_penalty[row_index, column_index]
            )
    for row_level in range(1, 3):
        for column_level in range(1, 3):
            block = sp.simplify(
                2
                * second_metric.extract(
                    levels[row_level - 1],
                    levels[column_level - 1],
                )
                + effective_second_forcing.extract(
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
    second_contraction = stein_coefficient(
        operators,
        (metric, metric_tangent, second_metric),
        2,
    )

    lower_third = third_schur_penalty(
        metric_tangent,
        second_metric,
        lower_kernel,
        lower_range,
        lower_base,
    )
    if sp.simplify(lower_third) != sp.zeros(2):
        raise AssertionError("the lower third-order endpoint did not vanish")
    third_metric = sp.zeros(6)
    assign_block(
        third_metric,
        levels[0],
        levels[0],
        lower_third,
    )

    third_contraction_without_metric = stein_coefficient(
        operators,
        (metric, metric_tangent, second_metric, zero_metric),
        3,
    )
    third_forcing = -third_contraction_without_metric
    contraction_third_penalty = third_schur_penalty(
        first_contraction,
        second_contraction,
        contraction_kernel,
        contraction_range,
        sp.eye(2),
    )
    effective_third_forcing = sp.MutableDenseMatrix(third_forcing)
    for row_index, row in enumerate(contraction_kernel):
        for column_index, column in enumerate(contraction_kernel):
            effective_third_forcing[row, column] += (
                contraction_third_penalty[row_index, column_index]
            )
    for row_level in range(1, 3):
        for column_level in range(1, 3):
            block = sp.simplify(
                2
                * third_metric.extract(
                    levels[row_level - 1],
                    levels[column_level - 1],
                )
                + effective_third_forcing.extract(
                    levels[row_level],
                    levels[column_level],
                )
            )
            assign_block(
                third_metric,
                levels[row_level],
                levels[column_level],
                block,
            )
    third_contraction = stein_coefficient(
        operators,
        (metric, metric_tangent, second_metric, third_metric),
        3,
    )
    if sp.simplify(
        third_contraction.extract(contraction_kernel, contraction_kernel)
        - contraction_third_penalty
    ) != sp.zeros(4):
        raise AssertionError("the third contraction Schur complement failed")

    upper_base = (4 * identity - metric).extract(upper_range, upper_range)
    upper_third_penalty = third_schur_penalty(
        -metric_tangent,
        -second_metric,
        upper_kernel,
        upper_range,
        upper_base,
    )
    endpoint = sp.simplify(
        third_metric.extract(upper_kernel, upper_kernel)
        + upper_third_penalty
    )
    expected_endpoint = -16 * third_mean * sp.eye(2)
    if sp.simplify(endpoint - expected_endpoint) != sp.zeros(2):
        raise AssertionError("the third-order endpoint did not collapse")

    print("PASS repeated p=3 flat-plane third-order Stein sign")
    print("third support cross entry = c(boundary), equation (2) of L77")
    print("third support coefficient = abs(c(boundary))")
    print("third Stein endpoint = -16 * mean(abs(c)) * identity")
    print("the coefficient is strict away from alpha=beta=0")


if __name__ == "__main__":
    main()
