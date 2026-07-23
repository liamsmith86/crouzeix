#!/usr/bin/env python3
"""Prove the repeated-C3 multiwinner/loser cross Gram penalty."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_stein_sign import assign_block


def main() -> None:
    root_two = sp.sqrt(2)
    boundary = sp.symbols("boundary", nonzero=True)
    gap = sp.symbols("gap", positive=True)
    top_support_mean = sp.symbols("top_support_mean", real=True)
    conformal_first, conformal_first_conjugate = sp.symbols(
        "conformal_first conformal_first_conjugate"
    )
    alpha = sp.symbols("alpha_0:3")
    alpha_conjugate = sp.symbols("alpha_conjugate_0:3")
    beta = sp.symbols("beta_0:3")
    beta_conjugate = sp.symbols("beta_conjugate_0:3")
    winner_edge = sp.symbols("winner_edge_0:3")
    winner_edge_conjugate = sp.symbols("winner_edge_conjugate_0:3")
    adjoint_substitution = {
        conformal_first: conformal_first_conjugate,
        conformal_first_conjugate: conformal_first,
    }
    for values, adjoints in (
        (alpha, alpha_conjugate),
        (beta, beta_conjugate),
        (winner_edge, winner_edge_conjugate),
    ):
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
    physical_direction = sp.diag(gap * crabb, gap * crabb, -gap * crabb)
    metric_tangent = sp.zeros(9)

    def optimized_cross_metric(values: tuple[sp.Symbol, ...]) -> sp.Matrix:
        """Return L81's optimized winner/loser metric block."""

        return sp.Matrix(
            [
                [0, root_two * values[0] / 4, 0],
                [
                    4 * root_two * values[2] / 3,
                    -2 * root_two * values[1],
                    2 * root_two * values[0],
                ],
                [0, 2 * root_two * values[2] / 3, 0],
            ]
        )

    def tied_cross_metric(values: tuple[sp.Symbol, ...]) -> sp.Matrix:
        """Return L76's canonical tied-winner metric block."""

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

    for winner, (values, adjoints) in enumerate(
        (
            (alpha, alpha_conjugate),
            (beta, beta_conjugate),
        )
    ):
        cross_to_winner = sp.zeros(3)
        cross_from_winner = sp.zeros(3)
        for index, (first, second) in enumerate(generators()):
            cross_to_winner += adjoints[index] * first
            cross_from_winner += values[index] * second
        winner_slice = slice(3 * winner, 3 * winner + 3)
        loser_slice = slice(6, 9)
        physical_direction[loser_slice, winner_slice] = cross_to_winner
        physical_direction[winner_slice, loser_slice] = cross_from_winner
        cross_metric = optimized_cross_metric(values)
        metric_tangent[winner_slice, loser_slice] = cross_metric
        metric_tangent[loser_slice, winner_slice] = formal_adjoint(
            cross_metric
        )

    internal_to_left = sp.zeros(3)
    internal_from_left = sp.zeros(3)
    for index, (first, second) in enumerate(generators()):
        internal_to_left += winner_edge_conjugate[index] * first
        internal_from_left += winner_edge[index] * second
    physical_direction[3:6, 0:3] = internal_to_left
    physical_direction[0:3, 3:6] = internal_from_left
    internal_metric = tied_cross_metric(winner_edge)
    metric_tangent[0:3, 3:6] = internal_metric
    metric_tangent[3:6, 0:3] = formal_adjoint(internal_metric)

    metric_tangent[6:9, 6:9] = sp.diag(gap, 0, -gap)
    pulled_first = physical_direction - gap * base
    pulled_second = (
        -gap * physical_direction
        + gap**2 * base
        - top_support_mean * base
        - conformal_first * base**2
    )

    support = (crabb / boundary + boundary * crabb.T) / 2
    reduced_resolvent = sp.eye(3) - support / 4 - 3 * support**2 / 4
    support_vector = sp.Matrix([1 / boundary, root_two, boundary]) / 2
    support_adjoint = sp.Matrix(
        [[boundary, root_two, 1 / boundary]]
    ) / 2
    support_perturbation = (
        physical_direction / boundary
        + boundary * formal_adjoint(physical_direction)
    ) / 2
    top_vectors = sp.zeros(9, 2)
    top_adjoint = sp.zeros(2, 9)
    for winner in range(2):
        winner_slice = slice(3 * winner, 3 * winner + 3)
        top_vectors[winner_slice, winner] = support_vector
        top_adjoint[winner, winner_slice] = support_adjoint
    effective_support = sp.simplify(
        top_adjoint
        * support_perturbation
        * sp.diag(
            reduced_resolvent,
            reduced_resolvent,
            reduced_resolvent,
        )
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
    winner_levels = tuple(level[:2] for level in levels)
    lower_kernel = winner_levels[0]
    lower_range = levels[1] + levels[2]
    upper_kernel = winner_levels[2]
    upper_range = levels[0] + levels[1]
    contraction_kernel = winner_levels[1] + winner_levels[2]
    contraction_range = levels[0]
    losing_contraction = (levels[1][2], levels[2][2])
    lower_base = metric - identity
    upper_base = 4 * identity - metric

    first_contraction = sp.simplify(
        metric_tangent
        - formal_adjoint(base) * metric_tangent * base
        - formal_adjoint(pulled_first) * metric * base
        - formal_adjoint(base) * metric * pulled_first
    )
    if first_contraction.extract(
        contraction_kernel,
        contraction_kernel,
    ) != sp.zeros(4):
        raise AssertionError("the winner contraction face moved")
    if first_contraction.extract(
        contraction_kernel,
        losing_contraction,
    ) != sp.zeros(4, 2):
        raise AssertionError("the winner/loser first contraction did not split")
    if sp.simplify(
        first_contraction.extract(
            losing_contraction,
            losing_contraction,
        )
        - sp.diag(6 * gap, 15 * gap)
    ) != sp.zeros(2):
        raise AssertionError("the losing first contraction was not strict")

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
        winner_levels[0],
        winner_levels[0],
        lower_penalty,
    )
    for row_level in range(1, 3):
        for column_level in range(1, 3):
            block = sp.simplify(
                2
                * second_metric.extract(
                    winner_levels[row_level - 1],
                    winner_levels[column_level - 1],
                )
                + effective_forcing.extract(
                    winner_levels[row_level],
                    winner_levels[column_level],
                )
            )
            assign_block(
                second_metric,
                winner_levels[row_level],
                winner_levels[column_level],
                block,
            )
    endpoint = sp.simplify(
        second_metric.extract(upper_kernel, upper_kernel) + upper_penalty
    )

    coefficient_columns = (
        sp.Matrix([[alpha[index]], [beta[index]]])
        for index in range(3)
    )
    coefficient_grams = tuple(
        column * formal_adjoint(column) for column in coefficient_columns
    )
    expected_endpoint = (
        16 * (support_mean - top_support_mean * sp.eye(2))
        - 8
        * sp.Matrix(
            [
                [0, winner_edge[1]],
                [winner_edge_conjugate[1], 0],
            ]
        )
        ** 2
        - sp.Rational(25, 8) * coefficient_grams[0]
        - 8 * coefficient_grams[1]
        - sp.Rational(50, 9) * coefficient_grams[2]
    )
    if sp.simplify(endpoint - expected_endpoint) != sp.zeros(2):
        raise AssertionError("the multiwinner Gram endpoint failed")

    print("PASS repeated p=3 multiwinner/loser cross penalty")
    print("losing first contraction slack = diag(6*gap, 15*gap)")
    print("winner endpoint = 16*(mean Q - mean lambda_max(Q)*I) - 8*H_1^2")
    print("  -25*A_0*A_0^*/8 - 8*A_1*A_1^* - 50*A_2*A_2^*/9")
    print("the radial gap and first conformal mode cancel exactly")


if __name__ == "__main__":
    main()
