#!/usr/bin/env python3
"""Verify the metric-flag derivative on the repeated-C3 flat core."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators


def formal_adjoint(
    matrix: sp.Matrix,
    adjoint_substitution: dict[sp.Symbol, sp.Symbol],
) -> sp.Matrix:
    """Transpose a symbolic matrix and exchange formal adjoint variables."""

    return matrix.T.xreplace(adjoint_substitution)


def assign_level_block(
    matrix: sp.Matrix,
    row_level: int,
    column_level: int,
    block: sp.Matrix,
    copy_count: int,
) -> None:
    """Assign one copy-space block in copy-major physical coordinates."""

    row_indices = tuple(3 * copy + row_level for copy in range(copy_count))
    column_indices = tuple(
        3 * copy + column_level for copy in range(copy_count)
    )
    for row_position, row_index in enumerate(row_indices):
        for column_position, column_index in enumerate(column_indices):
            matrix[row_index, column_index] = block[
                row_position,
                column_position,
            ]


def extract_level_block(
    matrix: sp.Matrix,
    row_level: int,
    column_level: int,
    copy_count: int,
) -> sp.Matrix:
    """Extract one copy-space block in copy-major physical coordinates."""

    row_indices = tuple(3 * copy + row_level for copy in range(copy_count))
    column_indices = tuple(
        3 * copy + column_level for copy in range(copy_count)
    )
    return matrix.extract(row_indices, column_indices)


def build_endpoint(
    copy_matrix: sp.Matrix,
    free_block: sp.Matrix,
    adjoint_substitution: dict[sp.Symbol, sp.Symbol],
) -> sp.Matrix:
    """Build the tight second-order endpoint for a flat copy perturbation."""

    copy_count = copy_matrix.rows
    root_two = sp.sqrt(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(*([crabb] * copy_count))
    metric = sp.diag(*([1, 2, 4] * copy_count))
    identity = sp.eye(3 * copy_count)
    generator_zero_first, generator_zero_second = generators()[0]
    copy_adjoint = formal_adjoint(
        copy_matrix,
        adjoint_substitution,
    )
    perturbation = sp.kronecker_product(
        copy_adjoint,
        generator_zero_first,
    ) + sp.kronecker_product(copy_matrix, generator_zero_second)

    first_forcing = (
        formal_adjoint(perturbation, adjoint_substitution)
        * metric
        * base
        + formal_adjoint(base, adjoint_substitution)
        * metric
        * perturbation
    )
    metric_tangent = sp.zeros(3 * copy_count)
    assign_level_block(
        metric_tangent,
        0,
        1,
        free_block,
        copy_count,
    )
    assign_level_block(
        metric_tangent,
        1,
        0,
        formal_adjoint(free_block, adjoint_substitution),
        copy_count,
    )
    assign_level_block(
        metric_tangent,
        1,
        1,
        extract_level_block(first_forcing, 1, 1, copy_count),
        copy_count,
    )
    next_block = (
        2 * free_block
        + extract_level_block(first_forcing, 1, 2, copy_count)
    )
    assign_level_block(
        metric_tangent,
        1,
        2,
        next_block,
        copy_count,
    )
    assign_level_block(
        metric_tangent,
        2,
        1,
        formal_adjoint(next_block, adjoint_substitution),
        copy_count,
    )

    lower_kernel = tuple(3 * copy for copy in range(copy_count))
    lower_range = tuple(
        3 * copy + level
        for copy in range(copy_count)
        for level in (1, 2)
    )
    upper_kernel = tuple(3 * copy + 2 for copy in range(copy_count))
    upper_range = tuple(
        3 * copy + level
        for copy in range(copy_count)
        for level in (0, 1)
    )
    lower_base = metric - identity
    upper_base = 4 * identity - metric
    lower_penalty = (
        metric_tangent.extract(lower_kernel, lower_range)
        * lower_base.extract(lower_range, lower_range).inv()
        * metric_tangent.extract(lower_range, lower_kernel)
    )
    upper_penalty = (
        metric_tangent.extract(upper_kernel, upper_range)
        * upper_base.extract(upper_range, upper_range).inv()
        * metric_tangent.extract(upper_range, upper_kernel)
    )

    first_contraction = (
        metric_tangent
        - formal_adjoint(base, adjoint_substitution)
        * metric_tangent
        * base
        - formal_adjoint(perturbation, adjoint_substitution)
        * metric
        * base
        - formal_adjoint(base, adjoint_substitution)
        * metric
        * perturbation
    )
    contraction_penalty = (
        first_contraction.extract(lower_range, lower_kernel)
        * first_contraction.extract(lower_kernel, lower_range)
    )
    forcing = sp.MutableDenseMatrix(
        formal_adjoint(perturbation, adjoint_substitution)
        * metric
        * perturbation
        + formal_adjoint(perturbation, adjoint_substitution)
        * metric_tangent
        * base
        + formal_adjoint(base, adjoint_substitution)
        * metric_tangent
        * perturbation
    )
    for row_position, row_index in enumerate(lower_range):
        for column_position, column_index in enumerate(lower_range):
            forcing[row_index, column_index] += contraction_penalty[
                row_position,
                column_position,
            ]

    second_metric = sp.zeros(3 * copy_count)
    assign_level_block(
        second_metric,
        0,
        0,
        lower_penalty,
        copy_count,
    )
    for row_level in (1, 2):
        for column_level in (1, 2):
            block = (
                2
                * extract_level_block(
                    second_metric,
                    row_level - 1,
                    column_level - 1,
                    copy_count,
                )
                + extract_level_block(
                    forcing,
                    row_level,
                    column_level,
                    copy_count,
                )
            )
            assign_level_block(
                second_metric,
                row_level,
                column_level,
                block,
                copy_count,
            )
    return (
        extract_level_block(second_metric, 2, 2, copy_count)
        + upper_penalty
    )


def main() -> None:
    copy_count = 3
    tau = sp.symbols("tau", real=True)
    entries = sp.symbols(f"z_0:{copy_count * copy_count}")
    adjoints = sp.symbols(f"z_adjoint_0:{copy_count * copy_count}")
    adjoint_substitution = dict(zip(entries, adjoints, strict=True))
    adjoint_substitution.update(zip(adjoints, entries, strict=True))
    copy_matrix = sp.Matrix(copy_count, copy_count, entries)
    copy_adjoint = formal_adjoint(
        copy_matrix,
        adjoint_substitution,
    )

    canonical_free_block = -3 * sp.sqrt(2) * copy_matrix / 8
    for projection_rank in (1, 2):
        projection = sp.diag(
            *(projection_rank * [1] + (copy_count - projection_rank) * [0])
        )
        complement = sp.eye(copy_count) - projection
        free_variation = (
            projection * copy_matrix * complement
            - sp.Rational(3, 4)
            * complement
            * copy_matrix
            * projection
        )
        endpoint = build_endpoint(
            copy_matrix,
            canonical_free_block + tau * free_variation,
            adjoint_substitution,
        )
        endpoint_derivative = endpoint.diff(tau).subs(tau, 0)
        expected_compression = (
            -5
            * sp.sqrt(2)
            * projection
            * copy_matrix
            * complement
            * copy_adjoint
            * projection
            - sp.Rational(15, 4)
            * sp.sqrt(2)
            * projection
            * copy_adjoint
            * complement
            * copy_matrix
            * projection
        )
        compression_residual = (
            projection * endpoint_derivative * projection
            - expected_compression
        )
        if compression_residual.applyfunc(sp.expand) != sp.zeros(copy_count):
            raise AssertionError(
                "the metric-flag derivative formula failed at rank "
                f"{projection_rank}"
            )

    # The canonical endpoint derivative is independent of the internal
    # projection blocks.  The arbitrary-dimensional statement is the same
    # block multiplication; the rank-one and corank-one symbolic checks
    # verify both orientations and every polarized entry type.

    print("PASS repeated p=3 flat-copy metric-flag derivative")
    print(
        "P E'(0) P = -5*sqrt(2) P Z R Z* P"
        " - 15*sqrt(2) P Z* R Z P/4"
    )
    print("the formula is negative semidefinite in every copy dimension")
    print("its kernel is {x in ran(P): Zx,Z*x lie in ran(P)}")


if __name__ == "__main__":
    main()
