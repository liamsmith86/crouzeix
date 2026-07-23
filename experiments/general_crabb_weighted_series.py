#!/usr/bin/env python3
"""Sparse exact weighted series near an arbitrary Crabb block.

This module extends the scalar Laurent-series machinery used by the ``p=3``
proof to arbitrary fixed dimension.  The simple top support eigenpair is
lifted by a reduced-resolvent recurrence, avoiding a symbolic determinant.
It is intended for weighted transverse-jet discovery and exact checks.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from p3_sparse_series import (
    Laurent,
    Series,
    clean,
    lp_add,
    lp_mul,
    lp_scale,
    lp_shift,
    lp_star,
    polynomial_at_matrix_series,
    series_exp,
)


def crabb_matrix(dimension: int) -> sp.Matrix:
    """Return the exact symbolic ``dimension x dimension`` Crabb block."""

    if dimension < 3:
        raise ValueError("the Crabb dimension must be at least three")
    length = dimension - 1
    matrix = sp.zeros(dimension)
    for index in range(length):
        matrix[index, index + 1] = (
            sp.sqrt(2) if index in (0, length - 1) else 1
        )
    return matrix


def laurent_matrix_vector(
    matrix: Sequence[Sequence[Laurent]], vector: Sequence[Laurent]
) -> list[Laurent]:
    """Multiply a Laurent-polynomial matrix by a Laurent vector."""

    result: list[Laurent] = []
    for row in matrix:
        value: Laurent = {}
        for entry, component in zip(row, vector, strict=True):
            value = lp_add(value, lp_mul(entry, component))
        result.append(value)
    return result


def constant_matrix_vector(
    matrix: sp.Matrix, vector: Sequence[Laurent]
) -> list[Laurent]:
    """Multiply a constant symbolic matrix by a Laurent vector."""

    result: list[Laurent] = []
    for row in range(matrix.rows):
        value: Laurent = {}
        for column in range(matrix.cols):
            value = lp_add(
                value,
                lp_scale(vector[column], matrix[row, column]),
            )
        result.append(value)
    return result


def constant_dot(vector: sp.Matrix, value: Sequence[Laurent]) -> Laurent:
    """Return the bilinear dot product of a constant and Laurent vector."""

    result: Laurent = {}
    for index, component in enumerate(value):
        result = lp_add(result, lp_scale(component, vector[index]))
    return result


def support_series(
    path: Sequence[sp.Matrix], order: int
) -> tuple[Series, list[list[Laurent]]]:
    """Lift the top support eigenvalue and intermediate-normalized eigenvector.

    ``path[j]`` is the coefficient of the real matrix path at the series
    parameter's degree ``j``.
    """

    if not path:
        raise ValueError("a nonempty matrix path is required")
    dimension = path[0].rows
    if any(coefficient.shape != (dimension, dimension) for coefficient in path):
        raise ValueError("all path coefficients must have the same square shape")
    if any(entry.is_real is not True for matrix in path for entry in matrix):
        raise ValueError("the sparse support recurrence expects a real path")

    length = dimension - 1
    endpoint = 1 / sp.sqrt(2)
    top = sp.Matrix(
        [
            (endpoint if index in (0, length) else 1) / sp.sqrt(length)
            for index in range(dimension)
        ]
    )
    base = path[0]
    support_base = (base + base.T) / 2
    reduced_inverse = (
        sp.eye(dimension) - support_base + top * top.T
    ).inv() - top * top.T

    support: list[list[list[Laurent]]] = []
    for degree in range(order + 1):
        coefficient = path[degree] if degree < len(path) else sp.zeros(dimension)
        matrix: list[list[Laurent]] = []
        for row in range(dimension):
            matrix_row: list[Laurent] = []
            for column in range(dimension):
                grade = column - row - 1
                value: Laurent = {}
                if coefficient[row, column] != 0:
                    value[grade] = coefficient[row, column] / 2
                if coefficient[column, row] != 0:
                    reverse_grade = grade + 2
                    value[reverse_grade] = (
                        value.get(reverse_grade, 0)
                        + coefficient[column, row] / 2
                    )
                matrix_row.append(value)
            matrix.append(matrix_row)
        support.append(matrix)

    eigenvalue: Series = [{0: sp.Integer(1)}, *[{} for _ in range(order)]]
    eigenvector: list[list[Laurent]] = [
        [{0: top[index]} for index in range(dimension)]
    ]
    for degree in range(1, order + 1):
        forcing: list[Laurent] = [{} for _ in range(dimension)]
        for source_degree in range(1, degree + 1):
            if source_degree >= len(support):
                continue
            image = laurent_matrix_vector(
                support[source_degree],
                eigenvector[degree - source_degree],
            )
            forcing = [
                lp_add(left, right)
                for left, right in zip(forcing, image, strict=True)
            ]
        eigenvalue[degree] = constant_dot(top, forcing)
        for source_degree in range(1, degree + 1):
            product = [
                lp_mul(
                    eigenvalue[source_degree],
                    eigenvector[degree - source_degree][index],
                )
                for index in range(dimension)
            ]
            forcing = [
                lp_add(left, lp_scale(right, -1))
                for left, right in zip(forcing, product, strict=True)
            ]
        eigenvector.append(constant_matrix_vector(reduced_inverse, forcing))
    return eigenvalue, eigenvector


def inverse_riemann_series(
    path: Sequence[sp.Matrix], order: int
) -> tuple[list[Laurent], list[sp.Matrix]]:
    """Return the inverse Riemann map and pullback for a Crabb-based path."""

    support_value, _ = support_series(path, order)
    boundary: Series = []
    for coefficient in support_value:
        boundary.append(
            {
                mode + 1: clean((1 - mode) * value)
                for mode, value in coefficient.items()
                if clean((1 - mode) * value) != 0
            }
        )

    angle_shift: Series = [{}, *[{} for _ in range(order)]]
    inverse_map: list[Laurent] = [{1: sp.Integer(1)}]

    def reparameterized_coefficient(degree: int) -> Laurent:
        result: Laurent = {}
        modes = {mode for coefficient in boundary for mode in coefficient}
        for mode in modes:
            exponent = [lp_scale(value, sp.I * mode) for value in angle_shift]
            phase = series_exp(exponent, degree)
            for boundary_degree in range(degree + 1):
                value = boundary[boundary_degree].get(mode, 0)
                if value != 0:
                    result = lp_add(
                        result,
                        lp_shift(
                            lp_scale(
                                phase[degree - boundary_degree],
                                value,
                            ),
                            mode,
                        ),
                    )
        return result

    for degree in range(1, order + 1):
        unmatched = reparameterized_coefficient(degree)
        correction: Laurent = {}
        for mode, coefficient in unmatched.items():
            if mode <= 0:
                analytic_half = {mode - 1: sp.I * coefficient}
                correction = lp_add(correction, analytic_half)
                correction = lp_add(correction, lp_star(analytic_half))
        angle_shift[degree] = correction
        analytic = reparameterized_coefficient(degree)
        bad = {
            mode: value
            for mode, value in analytic.items()
            if mode <= 0 and value != 0
        }
        if bad:
            raise AssertionError((degree, bad))
        inverse_map.append(analytic)

    dimension = path[0].rows
    operator = [path[0]]
    for degree in range(1, order + 1):
        trial = [*operator, sp.zeros(dimension)]
        known = sp.zeros(dimension)
        for map_degree in range(1, degree + 1):
            evaluated = polynomial_at_matrix_series(
                inverse_map[map_degree],
                trial,
                degree - map_degree,
            )
            known += evaluated[degree - map_degree]
        target = path[degree] if degree < len(path) else sp.zeros(dimension)
        operator.append(target - known)
    return inverse_map, [matrix.applyfunc(clean) for matrix in operator]


def axis_metric_series(
    dimension: int, order: int, parameter: sp.Symbol
) -> tuple[list[sp.Expr], list[sp.Matrix]]:
    """Expand L117's exact diagonal metric at the disk endpoint."""

    length = dimension - 1

    def hyperbolic_secant(index: int) -> sp.Expr:
        absolute_index = abs(index)
        if absolute_index == 0:
            return sp.Integer(1)
        return (
            2
            * parameter**absolute_index
            / (1 + parameter ** (2 * absolute_index))
        )

    periodized = []
    # Terms beyond this range start after the requested truncation.
    radius = order // (2 * length) + 3
    for index in range(dimension):
        periodized.append(
            sum(
                hyperbolic_secant(index + 2 * length * shift)
                for shift in range(-radius, radius + 1)
            )
        )
    expressions = [
        sp.series(
            periodized[index]
            / periodized[0]
            / parameter**index,
            parameter,
            0,
            order + 1,
        )
        .removeO()
        .expand()
        for index in range(dimension)
    ]
    coefficients = [
        sp.diag(
            *[
                expression.coeff(parameter, degree)
                for expression in expressions
            ]
        )
        for degree in range(order + 1)
    ]
    return expressions, coefficients


def triple_series_coefficient(
    left: Sequence[sp.Matrix],
    middle: Sequence[sp.Matrix],
    right: Sequence[sp.Matrix],
    degree: int,
) -> sp.Matrix:
    """Return the coefficient of ``left.T * middle * right``."""

    dimension = left[0].rows
    result = sp.zeros(dimension)
    for left_degree in range(degree + 1):
        for middle_degree in range(degree - left_degree + 1):
            right_degree = degree - left_degree - middle_degree
            result += (
                left[left_degree].T
                * middle[middle_degree]
                * right[right_degree]
            )
    return result


def weighted_condition_derivative(
    dimension: int,
    perturbation: sp.Matrix,
    extra_order: int = 0,
) -> sp.Expr:
    """Differentiate the rank-one envelope along ``c^L perturbation``.

    The result is a power series in ``c``.  The defect vector is held on the
    exact L117 optimizer, which is legitimate by the envelope theorem.
    """

    if perturbation.shape != (dimension, dimension):
        raise ValueError("the perturbation has the wrong shape")
    if any(entry.is_real is not True for entry in perturbation):
        raise ValueError("the weighted checker expects a real perturbation")

    length = dimension - 1
    order = 2 * length + extra_order
    parameter, transverse = sp.symbols(
        "ellipse_parameter transverse", real=True
    )
    base = crabb_matrix(dimension)
    path = [sp.zeros(dimension) for _ in range(length + 1)]
    path[0] = base
    path[1] = base.T
    path[length] += transverse * perturbation

    _, operator = inverse_riemann_series(path, order)
    axis_operator = [
        coefficient.subs(transverse, 0) for coefficient in operator
    ]
    operator_tangent = [
        coefficient.diff(transverse).subs(transverse, 0)
        for coefficient in operator
    ]
    metric_expressions, metric = axis_metric_series(
        dimension, order, parameter
    )

    forcing = [
        (
            triple_series_coefficient(
                operator_tangent,
                metric,
                axis_operator,
                degree,
            )
            + triple_series_coefficient(
                axis_operator,
                metric,
                operator_tangent,
                degree,
            )
        ).applyfunc(sp.expand)
        for degree in range(order + 1)
    ]

    metric_tangent: list[sp.Matrix] = []
    for degree in range(order + 1):
        right_hand_side = forcing[degree]
        for left_degree in range(degree + 1):
            for metric_degree in range(degree - left_degree + 1):
                right_degree = degree - left_degree - metric_degree
                if metric_degree < len(metric_tangent):
                    right_hand_side += (
                        axis_operator[left_degree].T
                        * metric_tangent[metric_degree]
                        * axis_operator[right_degree]
                    )
        coefficient = sp.zeros(dimension)
        for power in range(dimension):
            coefficient += (
                (base.T**power)
                * right_hand_side
                * (base**power)
            )
        metric_tangent.append(coefficient.applyfunc(sp.simplify))

    lower_derivative = sum(
        coefficient[0, 0] * parameter**degree
        for degree, coefficient in enumerate(metric_tangent)
    )
    upper_derivative = sum(
        coefficient[length, length] * parameter**degree
        for degree, coefficient in enumerate(metric_tangent)
    )
    derivative = upper_derivative - metric_expressions[length] * lower_derivative
    return sp.series(
        derivative,
        parameter,
        0,
        order + 1,
    ).removeO().expand()
