"""Reusable exact geometry for the full-disk Hardy residual."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from crabb_circular_normal_quadratic_exact import (
    disk_model_from_hermitian_series,
)
from crabb_disk_toeplitz_quartic import extend
from crabb_full_disk_base_jet_exact import toeplitz_direction
from crabb_full_disk_sixth_hardy_factor import (
    hardy_residual_from_disk_series,
    matrix_is_zero,
)


def reversal_matrix(size: int) -> sp.Matrix:
    """Return the order-reversing permutation matrix."""

    matrix = sp.zeros(size)
    for index in range(size):
        matrix[index, size - 1 - index] = 1
    return matrix


def skew_part(matrix: sp.Matrix) -> sp.Matrix:
    """Return the transpose-skew part of a complex matrix."""

    return ((matrix - matrix.T) / 2).applyfunc(sp.expand)


def weighted_projection(
    residual: sp.Matrix,
    length: int,
    mode: int,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the two weighted grade anti-diagonal projections."""

    lower = sum(
        (
            (right - left) * residual[left - 1, right - 1]
            for left in range(1, length)
            for right in range(left + 1, length)
            if left + right == length - mode
        ),
        sp.Integer(0),
    )
    upper = sum(
        (
            (right - left) * residual[left - 1, right - 1]
            for left in range(1, length)
            for right in range(left + 1, length)
            if left + right == length + mode
        ),
        sp.Integer(0),
    )
    return sp.expand(lower), sp.expand(upper)


def linearized_hardy_residual(
    hermitian_coefficient: sp.Matrix,
) -> sp.Matrix:
    """Return the Crabb differential of the finite Hardy residual.

    If the coefficient is an ``L``-square Hermitian matrix, the
    residual has size ``L-1``.  Rows are ordered by negative Hardy
    power, as in ``hardy_residual_from_disk_series``.
    """

    if hermitian_coefficient.rows != hermitian_coefficient.cols:
        raise ValueError("the Hermitian coefficient must be square")
    length = hermitian_coefficient.rows
    residual = sp.zeros(length - 1)
    for row in range(length - 1):
        reflected = length - 2 - row
        for column in range(length - 1):
            residual[row, column] = sp.expand(
                hermitian_coefficient[column + 1, reflected + 1]
                - hermitian_coefficient[column, reflected]
            )
    return residual


def zero_toeplitz_residual_lift(residual: sp.Matrix) -> sp.Matrix:
    """Lift a structured residual to a zero-Toeplitz Hermitian matrix.

    The residual is assumed to satisfy the disk-chart symmetry that
    ``J residual`` is Hermitian.  The returned coefficient has zero
    mean on every upper diagonal and is the unique such lift under
    ``linearized_hardy_residual``.
    """

    if residual.rows != residual.cols:
        raise ValueError("the Hardy residual must be square")
    size = residual.rows
    length = size + 1
    diagonal_difference = reversal_matrix(size) * residual
    coefficient = sp.zeros(length)
    for offset in range(length):
        diagonal_length = length - offset
        initial = -sum(
            (
                (diagonal_length - 1 - index)
                * diagonal_difference[index + offset, index]
                for index in range(diagonal_length - 1)
            ),
            sp.Integer(0),
        ) / diagonal_length
        for row in range(diagonal_length):
            value = initial + sum(
                (
                    diagonal_difference[index + offset, index]
                    for index in range(row)
                ),
                sp.Integer(0),
            )
            coefficient[row, row + offset] = sp.expand(value)
            coefficient[row + offset, row] = sp.expand(
                sp.conjugate(value)
            )
    if not matrix_is_zero(
        linearized_hardy_residual(coefficient) - residual
    ):
        raise RuntimeError("the zero-Toeplitz residual lift failed")
    return coefficient


def extended_disk_series(
    coefficients: Sequence[sp.Matrix],
) -> tuple[
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
]:
    """Extend Hermitian coefficients and construct their disk model."""

    if not coefficients:
        raise ValueError("the Hermitian series cannot be empty")
    hermitian = tuple(extend(coefficient) for coefficient in coefficients)
    operator, metric = disk_model_from_hermitian_series(
        list(hermitian)
    )
    return hermitian, tuple(metric), tuple(operator)


def formal_equality_jet(
    direction: Sequence[sp.Expr],
    order: int,
) -> tuple[
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
    tuple[sp.Matrix, ...],
]:
    """Recenter a Toeplitz direction onto the formal residual-zero graph."""

    if order < 1:
        raise ValueError("the formal equality jet needs positive order")
    length = len(direction)
    coefficients = [
        sp.eye(length) / 2,
        toeplitz_direction(direction),
        *[sp.zeros(length) for _ in range(order - 1)],
    ]
    for degree in range(2, order + 1):
        hermitian, _, operator = extended_disk_series(
            coefficients[: degree + 1]
        )
        residuals, _ = hardy_residual_from_disk_series(
            hermitian,
            operator,
        )
        coefficients[degree] = zero_toeplitz_residual_lift(
            -residuals[degree - 1]
        )
        hermitian, _, operator = extended_disk_series(
            coefficients[: degree + 1]
        )
        residuals, _ = hardy_residual_from_disk_series(
            hermitian,
            operator,
        )
        if not all(matrix_is_zero(residual) for residual in residuals):
            raise RuntimeError(
                f"formal equality recentering failed in degree {degree}"
            )
    return extended_disk_series(coefficients)
