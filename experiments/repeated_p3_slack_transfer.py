#!/usr/bin/env python3
"""Prove how normal Stein slack transfers to the weighted terminal endpoint."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_third_metric import (
    third_schur_penalty,
    tight_metrics_through_third,
)


def assert_zero_matrix(matrix: sp.Matrix, message: str) -> None:
    """Prove a symbolic matrix is zero entry by entry."""

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


def weighted_problem(
    diagonal: sp.Expr,
    operator_edge: sp.Expr,
    domain_edge: sp.Expr,
    third_mean: sp.Expr,
) -> tuple[
    tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix],
    sp.Matrix,
    sp.Matrix,
    sp.Matrix,
    sp.Matrix,
]:
    """Build the weighted-center operator, metric, and inherited slack.

    ``operator_edge`` is the edge present in the evaluated matrix.  The
    possibly different ``domain_edge`` records the edge that enlarged the
    frozen Riemann-map domain.  Setting the former to zero models ``f(N)``;
    setting both equal models ``f(A)``.
    """

    root_two = sp.sqrt(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb)
    metric = sp.diag(1, 2, 4, 1, 2, 4)
    copy = sp.Matrix([[diagonal, operator_edge], [0, -diagonal]])
    first_operator = flat_direction(copy)

    center = 3 * diagonal**2 / (8 * root_two)
    common_block = center * sp.Matrix(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    physical_second = sp.diag(common_block, common_block)
    second_mean = 5 * (2 * diagonal**2 + domain_edge**2) / 128
    second_cubic = 9 * diagonal**2 / 64
    operators = (
        base,
        first_operator,
        physical_second - second_mean * base,
        -second_mean * first_operator
        - second_cubic
        * polynomial_derivative(base, first_operator, 3)
        - third_mean * base,
    )

    edge_metric = sp.Matrix(
        [
            [0, -3 * root_two * operator_edge / 8, 0],
            [0, 0, 3 * root_two * operator_edge / 4],
            [0, 0, 0],
        ]
    )
    metric_tangent = sp.zeros(6)
    metric_tangent[0:3, 0:3] = loop_metric(diagonal)
    metric_tangent[3:6, 3:6] = -loop_metric(diagonal)
    metric_tangent[0:3, 3:6] = edge_metric
    metric_tangent[3:6, 0:3] = edge_metric.T

    transverse_second_mean = 5 * domain_edge**2 / 128
    second_slack = sp.diag(
        4 * transverse_second_mean,
        8 * transverse_second_mean,
        4 * transverse_second_mean,
        8 * transverse_second_mean,
    )
    mixed_slack = (
        15
        * sp.sqrt(2)
        * domain_edge**2
        * diagonal
        / 256
    )
    third_slack = sp.Matrix(
        [
            [4 * third_mean, mixed_slack, 0, 0],
            [mixed_slack, 8 * third_mean, 0, 0],
            [0, 0, 4 * third_mean, -mixed_slack],
            [0, 0, -mixed_slack, 8 * third_mean],
        ]
    )
    return (
        operators,
        metric,
        metric_tangent,
        second_slack,
        third_slack,
    )


def upper_endpoint_coefficients(
    metric: sp.Matrix,
    metric_tangent: sp.Matrix,
    second_metric: sp.Matrix,
    third_metric: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return the order-two and order-three upper Schur endpoints."""

    upper_kernel = (2, 5)
    upper_range = (0, 1, 3, 4)
    upper_base = (4 * sp.eye(6) - metric).extract(
        upper_range,
        upper_range,
    )
    second_penalty = (
        metric_tangent.extract(upper_kernel, upper_range)
        * upper_base.inv()
        * metric_tangent.extract(upper_range, upper_kernel)
    )
    third_penalty = third_schur_penalty(
        -metric_tangent,
        -second_metric,
        upper_kernel,
        upper_range,
        upper_base,
    )
    return (
        sp.simplify(
            second_metric.extract(upper_kernel, upper_kernel)
            + second_penalty
        ),
        sp.simplify(
            third_metric.extract(upper_kernel, upper_kernel)
            + third_penalty
        ),
    )


def weighted_endpoints(
    diagonal: sp.Expr,
    operator_edge: sp.Expr,
    domain_edge: sp.Expr,
    third_mean: sp.Expr,
    slack_fraction: sp.Expr,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return the quadratic and cubic weighted-center upper endpoints."""

    (
        operators,
        metric,
        metric_tangent,
        second_slack,
        third_slack,
    ) = weighted_problem(
        diagonal,
        operator_edge,
        domain_edge,
        third_mean,
    )
    second_metric, third_metric, *_ = tight_metrics_through_third(
        operators,
        metric,
        metric_tangent,
        second_stein_slack=slack_fraction * second_slack,
        third_stein_slack=slack_fraction * third_slack,
    )
    return upper_endpoint_coefficients(
        metric,
        metric_tangent,
        second_metric,
        third_metric,
    )


def main() -> None:
    diagonal, edge, third_mean, slack_fraction = sp.symbols(
        "diagonal edge third_mean slack_fraction",
        real=True,
    )
    actual_second, actual_third = weighted_endpoints(
        diagonal,
        edge,
        edge,
        third_mean,
        slack_fraction,
    )
    frozen_second, frozen_third = weighted_endpoints(
        diagonal,
        0,
        edge,
        third_mean,
        slack_fraction,
    )
    expected_third = (
        -16
        * (1 - slack_fraction)
        * third_mean
        * sp.eye(2)
    )
    assert_zero_matrix(
        actual_second
        - 5 * edge**2 * slack_fraction * sp.eye(2) / 8,
        "the actual quadratic Stein-slack transfer changed",
    )
    assert_zero_matrix(
        frozen_second
        + 5
        * edge**2
        * (1 - slack_fraction)
        * sp.eye(2)
        / 8,
        "the frozen-normal quadratic Stein-slack transfer changed",
    )
    assert_zero_matrix(
        actual_third - expected_third,
        "the actual weighted Stein-slack transfer changed",
    )
    assert_zero_matrix(
        frozen_third - expected_third,
        "the frozen-normal Stein-slack transfer changed",
    )

    baseline_problem = weighted_problem(diagonal, 0, 0, 0)
    frozen_problem = weighted_problem(
        diagonal,
        0,
        edge,
        third_mean,
    )
    baseline_metrics = tight_metrics_through_third(
        baseline_problem[0],
        baseline_problem[1],
        baseline_problem[2],
    )
    frozen_metrics = tight_metrics_through_third(
        frozen_problem[0],
        frozen_problem[1],
        frozen_problem[2],
        second_stein_slack=frozen_problem[3],
        third_stein_slack=frozen_problem[4],
    )
    assert_zero_matrix(
        frozen_metrics[0] - baseline_metrics[0],
        "the inherited second metric did not stay fixed",
    )
    assert_zero_matrix(
        frozen_metrics[1] - baseline_metrics[1],
        "the inherited third metric did not stay fixed",
    )

    print("PASS repeated p=3 weighted terminal Stein-slack transfer")
    print("E2(A)=5*a^2*theta/8; E2(N)=-5*a^2*(1-theta)/8")
    print("E3=-16*(1-theta)*m3*I for slack fraction theta")
    print("zero slack cancels the normal/transverse quadratic endpoints")
    print("f(A) and frozen f(N) then have the same negative cubic endpoint")
    print("full slack exactly retains the unperturbed normal metric jet")
    print("inherited full slack is flat; zero slack recovers L100")


if __name__ == "__main__":
    main()
