#!/usr/bin/env python3
"""Eliminate the finite L75 repeated-C3 second-order Stein certificate."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators


def assign_block(
    target: sp.Matrix,
    rows: tuple[int, ...],
    columns: tuple[int, ...],
    block: sp.Matrix,
) -> None:
    """Assign a dense block using explicit indices."""

    for row_index, row in enumerate(rows):
        for column_index, column in enumerate(columns):
            target[row, column] = block[row_index, column_index]


def main() -> None:
    root_two = sp.sqrt(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb)
    metric = sp.diag(1, 2, 4, 1, 2, 4)
    identity = sp.eye(6)
    conformal_mean = sp.symbols("conformal_mean", real=True)

    alpha_real = sp.symbols("alpha_real_0:3", real=True)
    alpha_imaginary = sp.symbols("alpha_imaginary_0:3", real=True)
    alpha = tuple(
        alpha_real[index] + sp.I * alpha_imaginary[index]
        for index in range(3)
    )
    perturbation = sp.zeros(6)
    for coefficient, (first, second) in zip(
        alpha,
        generators(),
        strict=True,
    ):
        perturbation[3:6, 0:3] += sp.conjugate(coefficient) * first
        perturbation[0:3, 3:6] += coefficient * second
    second_operator = -conformal_mean * base

    metric_cross = sp.Matrix(
        [
            [0, -3 * root_two * alpha[0] / 8, 0],
            [alpha[2] / root_two, -2 * root_two * alpha[1], 3 * root_two * alpha[0] / 4],
            [0, -root_two * alpha[2], 0],
        ]
    )
    metric_tangent = sp.zeros(6)
    metric_tangent[0:3, 3:6] = metric_cross
    metric_tangent[3:6, 0:3] = metric_cross.conjugate().T

    contraction_tangent = sp.simplify(
        metric_tangent
        - base.conjugate().T * metric_tangent * base
        - perturbation.conjugate().T * metric * base
        - base.conjugate().T * metric * perturbation
    )
    lower_kernel = (0, 3)
    lower_range = (1, 2, 4, 5)
    upper_kernel = (2, 5)
    upper_range = (0, 1, 3, 4)
    contraction_kernel = (1, 2, 4, 5)
    contraction_range = (0, 3)
    if metric_tangent.extract(lower_kernel, lower_kernel) != sp.zeros(2):
        raise AssertionError("the first lower metric face did not vanish")
    if metric_tangent.extract(upper_kernel, upper_kernel) != sp.zeros(2):
        raise AssertionError("the first upper metric face did not vanish")
    if contraction_tangent.extract(
        contraction_kernel,
        contraction_kernel,
    ) != sp.zeros(4):
        raise AssertionError("the first contraction face did not vanish")

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
    expected_lower_penalty = sp.diag(
        sp.Rational(9, 32) * alpha[0] * sp.conjugate(alpha[0]),
        sp.Rational(1, 2) * alpha[2] * sp.conjugate(alpha[2]),
    )
    expected_upper_penalty = sp.diag(
        alpha[2] * sp.conjugate(alpha[2]),
        sp.Rational(9, 16) * alpha[0] * sp.conjugate(alpha[0]),
    )
    if sp.simplify(lower_penalty - expected_lower_penalty) != sp.zeros(2):
        raise AssertionError("the lower endpoint penalty was not canonical")
    if sp.simplify(upper_penalty - expected_upper_penalty) != sp.zeros(2):
        raise AssertionError("the upper endpoint penalty was not canonical")

    forcing = sp.simplify(
        second_operator.conjugate().T * metric * base
        + base.conjugate().T * metric * second_operator
        + perturbation.conjugate().T * metric * perturbation
        + perturbation.conjugate().T * metric_tangent * base
        + base.conjugate().T * metric_tangent * perturbation
    )
    contraction_penalty = sp.simplify(
        contraction_tangent.extract(contraction_kernel, contraction_range)
        * contraction_tangent.extract(contraction_range, contraction_kernel)
    )

    # Construct the complete second metric.  On positive Crabb levels the
    # Stein equality is a two-dimensional weighted-shift recurrence.  Blocks
    # touching level zero away from the lower endpoint are free and are set to
    # zero before the recurrence.
    second_metric = sp.zeros(6)
    level_indices = tuple((level, level + 3) for level in range(3))
    assign_block(
        second_metric,
        level_indices[0],
        level_indices[0],
        lower_penalty,
    )
    effective_forcing = sp.MutableDenseMatrix(forcing)
    for row_index, row in enumerate(contraction_kernel):
        for column_index, column in enumerate(contraction_kernel):
            effective_forcing[row, column] += contraction_penalty[
                row_index,
                column_index,
            ]
    weights = (root_two, root_two)
    for row_level in range(1, 3):
        for column_level in range(1, 3):
            previous = second_metric.extract(
                level_indices[row_level - 1],
                level_indices[column_level - 1],
            )
            current_forcing = effective_forcing.extract(
                level_indices[row_level],
                level_indices[column_level],
            )
            block = sp.simplify(
                weights[row_level - 1]
                * weights[column_level - 1]
                * previous
                + current_forcing
            )
            assign_block(
                second_metric,
                level_indices[row_level],
                level_indices[column_level],
                block,
            )
    if sp.simplify(second_metric - second_metric.conjugate().T) != sp.zeros(6):
        raise AssertionError("the reconstructed second metric is not Hermitian")

    second_contraction = sp.simplify(
        second_metric
        - base.conjugate().T * second_metric * base
        - forcing
    )
    contraction_schur = sp.simplify(
        second_contraction.extract(contraction_kernel, contraction_kernel)
        - contraction_penalty
    )
    if contraction_schur != sp.zeros(4):
        raise AssertionError("the second contraction Schur complement failed")

    endpoint_matrix = sp.simplify(
        second_metric.extract(upper_kernel, upper_kernel) + upper_penalty
    )

    magnitude_squared = tuple(
        sp.expand(coefficient * sp.conjugate(coefficient))
        for coefficient in alpha
    )
    raw_scalar = (
        sp.Rational(5, 8) * magnitude_squared[0]
        - 4 * magnitude_squared[1]
        + sp.Rational(10, 9) * magnitude_squared[2]
        - 16 * conformal_mean
    )
    if sp.simplify(endpoint_matrix - raw_scalar * sp.eye(2)) != sp.zeros(2):
        raise AssertionError("the endpoint certificate did not collapse to a scalar")
    lower_schur = sp.simplify(
        second_metric.extract(lower_kernel, lower_kernel) - lower_penalty
    )
    upper_schur = sp.simplify(raw_scalar * sp.eye(2) - endpoint_matrix)
    if lower_schur != sp.zeros(2) or upper_schur != sp.zeros(2):
        raise AssertionError("an endpoint metric Schur complement failed")

    coupling_magnitude = sp.symbols("coupling_magnitude", nonnegative=True)
    exact_conformal_mean = (
        sp.Rational(5, 128) * magnitude_squared[0]
        + sp.Rational(1, 4) * magnitude_squared[1]
        + sp.Rational(5, 72) * magnitude_squared[2]
        + root_two * coupling_magnitude / (12 * sp.pi)
    )
    final_coefficient = sp.simplify(
        raw_scalar.subs(conformal_mean, exact_conformal_mean)
    )
    expected_final = (
        -8 * magnitude_squared[1]
        - 4 * root_two * coupling_magnitude / (3 * sp.pi)
    )
    if sp.simplify(final_coefficient - expected_final) != 0:
        raise AssertionError("the conformal mean did not force nonpositivity")

    print("PASS repeated p=3 second-order Stein sign")
    print(f"raw endpoint scalar = {raw_scalar}")
    print(f"conformal certificate = {expected_final}")
    print("equality requires alpha_1=0 (and then the coupling also vanishes)")


if __name__ == "__main__":
    main()
