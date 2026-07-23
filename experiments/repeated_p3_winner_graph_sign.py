#!/usr/bin/env python3
"""Prove the tied-winner graph second-order endpoint identity."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_stein_sign import assign_block


def main() -> None:
    root_two = sp.sqrt(2)
    boundary = sp.symbols("boundary", nonzero=True)
    top_support_mean = sp.symbols("top_support_mean", real=True)
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
    loop_zero = sp.symbols("loop_zero_0:3")
    loop_zero_conjugate = sp.symbols("loop_zero_conjugate_0:3")
    loop_one = sp.symbols("loop_one_0:3", real=True)
    common_two, common_two_conjugate = sp.symbols(
        "common_two common_two_conjugate"
    )
    common_three, common_three_conjugate = sp.symbols(
        "common_three common_three_conjugate"
    )
    adjoint_substitution: dict[sp.Symbol, sp.Symbol] = {
        conformal_first: conformal_first_conjugate,
        conformal_first_conjugate: conformal_first,
    }
    for values, adjoints in zip(coefficients, conjugates, strict=True):
        adjoint_substitution.update(zip(values, adjoints, strict=True))
        adjoint_substitution.update(zip(adjoints, values, strict=True))
    adjoint_substitution.update(
        zip(loop_zero, loop_zero_conjugate, strict=True)
    )
    adjoint_substitution.update(
        zip(loop_zero_conjugate, loop_zero, strict=True)
    )
    adjoint_substitution.update(
        {
            common_two: common_two_conjugate,
            common_two_conjugate: common_two,
            common_three: common_three_conjugate,
            common_three_conjugate: common_three,
        }
    )

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
    generator_pairs = generators()
    common_motion = common_two * sp.Matrix(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    common_motion[2, 0] += common_three

    def metric_cross(values: tuple[sp.Expr, sp.Expr, sp.Expr]) -> sp.Matrix:
        """Return the canonical L76 first-metric block."""

        return sp.Matrix(
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

    for edge_index, (left_copy, right_copy) in enumerate(edges):
        values = coefficients[edge_index]
        adjoints = conjugates[edge_index]
        cross_to_left = sp.zeros(3)
        cross_from_left = sp.zeros(3)
        for generator_index, (first, second) in enumerate(generator_pairs):
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

        tangent_block = metric_cross(values)
        metric_tangent[left_slice, right_slice] = tangent_block
        metric_tangent[right_slice, left_slice] = formal_adjoint(
            tangent_block
        )

    # L84's canonical diagonal support-kernel quotient consists of a complex
    # generator-zero loop and a real generator-one loop.  The latter is split
    # evenly between X_1 and Y_1 so its physical coefficient is loop_one.
    for copy in range(3):
        copy_slice = slice(3 * copy, 3 * copy + 3)
        loop_values = (
            loop_zero[copy],
            loop_one[copy] / 2,
            sp.Integer(0),
        )
        loop_adjoints = (
            loop_zero_conjugate[copy],
            loop_one[copy] / 2,
            sp.Integer(0),
        )
        loop_direction = sp.zeros(3)
        for index, (first, second) in enumerate(generator_pairs):
            loop_direction += loop_adjoints[index] * first
            loop_direction += loop_values[index] * second
        perturbation[copy_slice, copy_slice] = (
            loop_direction + common_motion
        )

        loop_block = metric_cross(loop_values)
        metric_tangent[copy_slice, copy_slice] = (
            loop_block + formal_adjoint(loop_block)
        )
        coefficient_matrices[1][copy, copy] = loop_one[copy]

    support = (crabb / boundary + boundary * crabb.T) / 2
    reduced_resolvent = (
        sp.eye(3) - support / 4 - 3 * support**2 / 4
    )
    support_vector = sp.Matrix([1 / boundary, root_two, boundary]) / 2
    support_adjoint = sp.Matrix(
        [[boundary, root_two, 1 / boundary]]
    ) / 2
    common_first_support = sp.factor(
        (
            support_adjoint
            * (
                common_motion / boundary
                + boundary * formal_adjoint(common_motion)
            )
            / 2
            * support_vector
        )[0]
    )
    expected_first_support = (
        root_two
        * (
            common_two / boundary**2
            + common_two_conjugate * boundary**2
        )
        / 4
        + (
            common_three / boundary**3
            + common_three_conjugate * boundary**3
        )
        / 8
    )
    if sp.simplify(common_first_support - expected_first_support) != 0:
        raise AssertionError("the common first-support formula failed")

    first_boundary_schwarz = (
        common_two_conjugate * boundary**2 / root_two
        + common_three_conjugate * boundary**3 / 4
    )
    adjoint_boundary_schwarz = first_boundary_schwarz.xreplace(
        adjoint_substitution
    ).subs(boundary, 1 / boundary)
    imaginary_boundary = (
        first_boundary_schwarz - adjoint_boundary_schwarz
    ) / (2 * sp.I)
    support_derivative = (
        sp.I * boundary * sp.diff(common_first_support, boundary)
    )
    normal_shift = sp.expand(imaginary_boundary - support_derivative)
    normal_shift_square = sp.expand(normal_shift**2 * boundary**6)
    normal_shift_mean = sp.simplify(
        normal_shift_square.coeff(boundary, 6)
    )
    expected_normal_shift_mean = (
        sp.Rational(9, 4) * common_two * common_two_conjugate
        + sp.Rational(1, 2) * common_three * common_three_conjugate
    )
    if sp.simplify(normal_shift_mean - expected_normal_shift_mean) != 0:
        raise AssertionError("the common normal-angle mean failed")
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
    first_map_frechet = (common_two_conjugate / root_two) * (
        perturbation * base**2
        + base * perturbation * base
        + base**2 * perturbation
    ) + (common_three_conjugate / 4) * (
        base * perturbation * base**2
        + base**2 * perturbation * base
    )
    second_normal_mean = top_support_mean - normal_shift_mean / 2
    second_operator = (
        -first_map_frechet
        - second_normal_mean * base
        - conformal_first * base**2
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
        - 16 * top_support_mean * sp.eye(3)
        - sp.Rational(21, 4)
        * common_three
        * common_three_conjugate
        * sp.eye(3)
    )
    if sp.simplify(endpoint - expected_endpoint) != sp.zeros(3):
        raise AssertionError("the winner-graph endpoint identity failed")

    print("PASS repeated p=3 complete tied-winner second-order sign")
    print("endpoint = 16*(mean Q - mean lambda_max(Q)*I) - 8*H_1^2")
    print("           - 21*abs(common mode 3)^2*I/4")
    print("mean support includes orientation-dependent mixed edge paths")
    print("the first conformal mode and common flat mode cancel exactly")
    print("matrix Jensen order makes the endpoint negative semidefinite")


if __name__ == "__main__":
    main()
