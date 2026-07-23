"""Shared tight third-order metric propagation for two repeated C3 blocks."""

from __future__ import annotations

import sympy as sp

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


def tight_metrics_through_third(
    operators: tuple[sp.Matrix, ...],
    metric: sp.Matrix,
    metric_tangent: sp.Matrix,
    second_free_block: sp.Matrix | None = None,
    third_free_block: sp.Matrix | None = None,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return tight metrics and Stein coefficients through order three.

    ``second_free_block`` is the optional level-zero/level-one copy block of
    the second metric.  It is needed when a first-order flat direction itself
    changes at the next weighted order.  ``third_free_block`` is the analogous
    block at order three; it first affects the order-four endpoint.
    """

    identity = sp.eye(6)
    lower_kernel = (0, 3)
    lower_range = (1, 2, 4, 5)
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
    second_metric = sp.zeros(6)
    assign_block(
        second_metric,
        levels[0],
        levels[0],
        lower_second,
    )
    if second_free_block is not None:
        assign_block(
            second_metric,
            levels[0],
            levels[1],
            second_free_block,
        )
        assign_block(
            second_metric,
            levels[1],
            levels[0],
            second_free_block.conjugate().T,
        )
    effective_second_forcing = sp.MutableDenseMatrix(
        -second_contraction_without_metric
    )
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
    third_metric = sp.zeros(6)
    assign_block(
        third_metric,
        levels[0],
        levels[0],
        lower_third,
    )
    if third_free_block is not None:
        assign_block(
            third_metric,
            levels[0],
            levels[1],
            third_free_block,
        )
        assign_block(
            third_metric,
            levels[1],
            levels[0],
            third_free_block.conjugate().T,
        )

    third_contraction_without_metric = stein_coefficient(
        operators,
        (metric, metric_tangent, second_metric, zero_metric),
        3,
    )
    contraction_third_penalty = third_schur_penalty(
        first_contraction,
        second_contraction,
        contraction_kernel,
        contraction_range,
        sp.eye(2),
    )
    effective_third_forcing = sp.MutableDenseMatrix(
        -third_contraction_without_metric
    )
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

    return (
        second_metric,
        third_metric,
        first_contraction,
        second_contraction,
        third_contraction,
    )


def tight_third_endpoint(
    operators: tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix],
    metric: sp.Matrix,
    metric_tangent: sp.Matrix,
    second_free_block: sp.Matrix | None = None,
) -> sp.Matrix:
    """Propagate the tight metrics and return the order-three endpoint."""

    identity = sp.eye(6)
    upper_kernel = (2, 5)
    upper_range = (0, 1, 3, 4)
    second_metric, third_metric, _, _, _ = (
        tight_metrics_through_third(
            operators,
            metric,
            metric_tangent,
            second_free_block,
        )
    )
    upper_base = (4 * identity - metric).extract(
        upper_range,
        upper_range,
    )
    upper_third_penalty = third_schur_penalty(
        -metric_tangent,
        -second_metric,
        upper_kernel,
        upper_range,
        upper_base,
    )
    return sp.simplify(
        third_metric.extract(upper_kernel, upper_kernel)
        + upper_third_penalty
    )


def fourth_schur_penalty(
    first: sp.Matrix,
    second: sp.Matrix,
    third: sp.Matrix,
    kernel: tuple[int, ...],
    positive_range: tuple[int, ...],
    base_range: sp.Matrix,
) -> sp.Matrix:
    """Return the order-four endpoint forced by a tight Schur complement."""

    inverse = base_range.inv()
    first_cross = first.extract(kernel, positive_range)
    second_cross = second.extract(kernel, positive_range)
    third_cross = third.extract(kernel, positive_range)
    first_range = first.extract(positive_range, positive_range)
    second_range = second.extract(positive_range, positive_range)
    return sp.simplify(
        third_cross * inverse * first_cross.conjugate().T
        + first_cross * inverse * third_cross.conjugate().T
        + second_cross * inverse * second_cross.conjugate().T
        - second_cross
        * inverse
        * first_range
        * inverse
        * first_cross.conjugate().T
        - first_cross
        * inverse
        * first_range
        * inverse
        * second_cross.conjugate().T
        - first_cross
        * inverse
        * second_range
        * inverse
        * first_cross.conjugate().T
        + first_cross
        * inverse
        * first_range
        * inverse
        * first_range
        * inverse
        * first_cross.conjugate().T
    )


def tight_fourth_endpoint(
    operators: tuple[
        sp.Matrix,
        sp.Matrix,
        sp.Matrix,
        sp.Matrix,
        sp.Matrix,
    ],
    metric: sp.Matrix,
    metric_tangent: sp.Matrix,
    second_free_block: sp.Matrix | None = None,
    third_free_block: sp.Matrix | None = None,
) -> sp.Matrix:
    """Propagate the tight metrics and return the order-four endpoint."""

    identity = sp.eye(6)
    lower_kernel = (0, 3)
    lower_range = (1, 2, 4, 5)
    upper_kernel = (2, 5)
    upper_range = (0, 1, 3, 4)
    contraction_kernel = (1, 2, 4, 5)
    contraction_range = (0, 3)
    levels = tuple((level, level + 3) for level in range(3))
    zero_metric = sp.zeros(6)

    (
        second_metric,
        third_metric,
        first_contraction,
        second_contraction,
        third_contraction,
    ) = tight_metrics_through_third(
        operators,
        metric,
        metric_tangent,
        second_free_block,
        third_free_block,
    )

    lower_base = (metric - identity).extract(
        lower_range,
        lower_range,
    )
    lower_fourth = fourth_schur_penalty(
        metric_tangent,
        second_metric,
        third_metric,
        lower_kernel,
        lower_range,
        lower_base,
    )
    fourth_metric = sp.zeros(6)
    assign_block(
        fourth_metric,
        levels[0],
        levels[0],
        lower_fourth,
    )

    fourth_contraction_without_metric = stein_coefficient(
        operators,
        (
            metric,
            metric_tangent,
            second_metric,
            third_metric,
            zero_metric,
        ),
        4,
    )
    contraction_fourth_penalty = fourth_schur_penalty(
        first_contraction,
        second_contraction,
        third_contraction,
        contraction_kernel,
        contraction_range,
        sp.eye(2),
    )
    effective_fourth_forcing = sp.MutableDenseMatrix(
        -fourth_contraction_without_metric
    )
    for row_index, row in enumerate(contraction_kernel):
        for column_index, column in enumerate(contraction_kernel):
            effective_fourth_forcing[row, column] += (
                contraction_fourth_penalty[
                    row_index,
                    column_index,
                ]
            )
    for row_level in range(1, 3):
        for column_level in range(1, 3):
            block = sp.simplify(
                2
                * fourth_metric.extract(
                    levels[row_level - 1],
                    levels[column_level - 1],
                )
                + effective_fourth_forcing.extract(
                    levels[row_level],
                    levels[column_level],
                )
            )
            assign_block(
                fourth_metric,
                levels[row_level],
                levels[column_level],
                block,
            )
    fourth_contraction = stein_coefficient(
        operators,
        (
            metric,
            metric_tangent,
            second_metric,
            third_metric,
            fourth_metric,
        ),
        4,
    )
    if sp.simplify(
        fourth_contraction.extract(
            contraction_kernel,
            contraction_kernel,
        )
        - contraction_fourth_penalty
    ) != sp.zeros(4):
        raise AssertionError(
            "the fourth contraction Schur complement failed"
        )

    upper_base = (4 * identity - metric).extract(
        upper_range,
        upper_range,
    )
    upper_fourth_penalty = fourth_schur_penalty(
        -metric_tangent,
        -second_metric,
        -third_metric,
        upper_kernel,
        upper_range,
        upper_base,
    )
    return sp.simplify(
        fourth_metric.extract(upper_kernel, upper_kernel)
        + upper_fourth_penalty
    )
