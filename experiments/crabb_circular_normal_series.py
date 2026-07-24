#!/usr/bin/env python3
"""Reusable exact weighted-series tools for circular-normal Crabb faces."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from rank_one_stein_series import (
    diagonal_gramian_condition_series,
    simple_diagonal_eigenvalue_coefficients,
    stein_gramian_series,
)


def ordinary_triple_series_coefficient(
    left: Sequence[sp.Matrix],
    middle: Sequence[sp.Matrix],
    right: Sequence[sp.Matrix],
    degree: int,
) -> sp.Matrix:
    """Return one coefficient of an ordinary three-factor matrix product."""

    dimension = left[0].rows
    result = sp.zeros(dimension)
    for left_degree in range(degree + 1):
        for middle_degree in range(degree - left_degree + 1):
            right_degree = degree - left_degree - middle_degree
            if (
                left_degree < len(left)
                and middle_degree < len(middle)
                and right_degree < len(right)
            ):
                result += (
                    left[left_degree]
                    * middle[middle_degree]
                    * right[right_degree]
                )
    return sp.simplify(result)


def ordinary_matrix_inverse_series(
    coefficients: Sequence[sp.Matrix],
) -> list[sp.Matrix]:
    """Invert an ordinary matrix power series."""

    dimension = coefficients[0].rows
    inverse = [coefficients[0].inv()]
    for degree in range(1, len(coefficients)):
        convolution = sum(
            (
                coefficients[source_degree]
                * inverse[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.zeros(dimension),
        )
        inverse.append(sp.simplify(-inverse[0] * convolution))
    return inverse


def scalar_series_quotient(
    numerator: Sequence[sp.Expr],
    denominator: Sequence[sp.Expr],
) -> list[sp.Expr]:
    """Divide two scalar power series with nonzero denominator constant."""

    quotient: list[sp.Expr] = []
    for degree in range(len(numerator)):
        known = sum(
            (
                quotient[source_degree]
                * denominator[degree - source_degree]
                for source_degree in range(degree)
            ),
            sp.Integer(0),
        )
        quotient.append(
            sp.simplify(
                (numerator[degree] - known) / denominator[0]
            )
        )
    return quotient


def inverse_square_root_from_series(
    gramian: Sequence[sp.Matrix],
) -> list[sp.Matrix]:
    """Solve ``G(epsilon) K(epsilon) G(epsilon) = I`` coefficientwise.

    The constant coefficient must be positive diagonal.  Later
    coefficients may be arbitrary symmetric matrices.  This is the
    general recurrence underlying :func:`inverse_square_root_series`.
    """

    if not gramian:
        raise ValueError("the Gramian series cannot be empty")
    constant = gramian[0]
    if constant != sp.diag(*constant.diagonal()):
        raise ValueError("the constant Gramian coefficient must be diagonal")
    dimension = constant.rows
    square_roots = [
        sp.sqrt(constant[index, index])
        for index in range(dimension)
    ]
    coefficients = [
        sp.diag(*[1 / value for value in square_roots])
    ]

    for degree in range(1, len(gramian)):
        known = sp.zeros(dimension)
        for left_degree in range(degree + 1):
            for middle_degree in range(degree - left_degree + 1):
                right_degree = degree - left_degree - middle_degree
                is_left_unknown = (
                    left_degree == degree
                    and middle_degree == 0
                    and right_degree == 0
                )
                is_right_unknown = (
                    left_degree == 0
                    and middle_degree == 0
                    and right_degree == degree
                )
                if is_left_unknown or is_right_unknown:
                    continue
                if (
                    left_degree < len(coefficients)
                    and right_degree < len(coefficients)
                ):
                    known += (
                        coefficients[left_degree]
                        * gramian[middle_degree]
                        * coefficients[right_degree]
                    )

        coefficient = sp.zeros(dimension)
        for row in range(dimension):
            for column in range(dimension):
                coefficient[row, column] = sp.simplify(
                    -known[row, column]
                    / (square_roots[row] + square_roots[column])
                )
        coefficients.append(coefficient)

    for degree in range(len(gramian)):
        residual = ordinary_triple_series_coefficient(
            coefficients,
            gramian,
            coefficients,
            degree,
        )
        target = sp.eye(dimension) if degree == 0 else sp.zeros(dimension)
        if residual != target:
            raise AssertionError(
                "the inverse-square-root recurrence did not regenerate"
            )
    return coefficients


def inverse_square_root_series(
    constant: sp.Matrix,
    tangent: sp.Matrix,
    order: int,
) -> list[sp.Matrix]:
    """Solve an affine ``G(epsilon) K(epsilon) G(epsilon) = I`` series."""

    dimension = constant.rows
    gramian = [
        constant,
        tangent,
        *[sp.zeros(dimension) for _ in range(order - 1)],
    ]
    return inverse_square_root_from_series(gramian)


def physical_reflected_path(
    dimension: int,
    equality_grade: int,
    strong_parameter: sp.Expr,
    strong_direction: sp.Matrix,
    strong_degree: int,
    order: int,
    *,
    amplitude: sp.Expr = sp.Integer(1),
    ellipse: sp.Expr = sp.Integer(1),
) -> list[sp.Matrix]:
    """Construct a disk-equality/ellipse/strong-normal weighted path."""

    length = dimension - 1
    if not 1 <= equality_grade <= length - 1:
        raise ValueError("the equality grade must lie between 1 and L-1")
    if not 1 <= strong_degree <= order:
        raise ValueError("the strong degree must lie between 1 and the order")
    if strong_direction.shape != (dimension, dimension):
        raise ValueError("the strong direction has the wrong shape")

    shift = sp.zeros(dimension)
    for index in range(length):
        shift[index, index + 1] = 1

    base_toeplitz = sp.zeros(dimension)
    tangent_toeplitz = sp.zeros(dimension)
    for index in range(length):
        base_toeplitz[index, index] = sp.Rational(1, 2)
    for offset in {equality_grade, length - equality_grade}:
        for row in range(length - offset):
            tangent_toeplitz[row, row + offset] = amplitude
            tangent_toeplitz[row + offset, row] = amplitude

    base_coordinate_gramian = (
        base_toeplitz
        + shift.T * base_toeplitz * shift
    )
    tangent_coordinate_gramian = (
        tangent_toeplitz
        + shift.T * tangent_toeplitz * shift
    )
    inverse_square_root = inverse_square_root_series(
        base_coordinate_gramian,
        tangent_coordinate_gramian,
        order,
    )
    toeplitz_series = [
        base_toeplitz,
        tangent_toeplitz,
        *[sp.zeros(dimension) for _ in range(order - 1)],
    ]
    shifted_inverse = [
        shift * coefficient
        for coefficient in inverse_square_root
    ]
    disk_operator = [
        2
        * ordinary_triple_series_coefficient(
            inverse_square_root,
            toeplitz_series,
            shifted_inverse,
            degree,
        )
        for degree in range(order + 1)
    ]

    path = []
    for degree in range(order + 1):
        coefficient = disk_operator[degree]
        if degree >= 1:
            coefficient += ellipse * disk_operator[degree - 1].T
        if degree == strong_degree:
            coefficient += strong_parameter * strong_direction
        path.append(sp.simplify(coefficient))
    return path


def optimized_defect_jets(
    operator: Sequence[sp.Matrix],
    epsilon: sp.Symbol,
    jet_count: int,
) -> sp.Matrix:
    """Solve successive rank-one defect stationarity jets exactly."""

    required_order = 2 * jet_count
    if len(operator) <= required_order:
        raise ValueError("the operator series is too short for the requested jets")

    dimension = operator[0].rows
    defect = sp.zeros(dimension, 1)
    defect[0] = 1
    for jet in range(1, jet_count + 1):
        variables = sp.symbols(
            f"defect_{jet}_1:{dimension}",
            real=True,
        )
        trial = (
            defect
            + epsilon**jet * sp.Matrix([0, *variables])
        )
        condition = diagonal_gramian_condition_series(
            operator,
            trial,
            epsilon,
            2 * jet,
        )
        coefficient = sp.expand(condition).coeff(epsilon, 2 * jet)
        solutions = sp.solve(
            [
                sp.diff(coefficient, variable)
                for variable in variables
            ],
            variables,
            dict=True,
            simplify=False,
        )
        if len(solutions) != 1:
            raise RuntimeError(
                f"defect jet {jet} was not uniquely stationary"
            )
        defect = sp.simplify(trial.subs(solutions[0]))
    return defect


def reciprocal_reversal_defect_jets(
    operator: Sequence[sp.Matrix],
    epsilon: sp.Symbol,
    jet_count: int,
) -> sp.Matrix:
    """Solve real persymmetric stationary defects from L164 self-duality.

    At jet ``j``, impose ``J P^(-1) J = alpha P`` only through degree
    ``j``.  Unlike :func:`optimized_defect_jets`, this requires operator
    and Stein data only through ``j`` rather than through ``2j``.
    """

    if len(operator) <= jet_count:
        raise ValueError("the operator series is too short for the requested jets")

    dimension = operator[0].rows
    reversal = sp.zeros(dimension)
    for index in range(dimension):
        reversal[index, dimension - 1 - index] = 1

    if any(
        sp.simplify(reversal * coefficient.T * reversal - coefficient)
        != sp.zeros(dimension)
        for coefficient in operator[: jet_count + 1]
    ):
        raise ValueError("reciprocal-reversal jets require a persymmetric path")

    defect = sp.eye(dimension)[:, 0]
    for jet in range(1, jet_count + 1):
        variables = sp.symbols(
            f"selfdual_{jet}_1:{dimension}",
            real=True,
        )
        trial = defect + epsilon**jet * sp.Matrix([0, *variables])
        gramian = stein_gramian_series(
            operator,
            trial,
            epsilon,
            jet,
        )
        inverse = ordinary_matrix_inverse_series(gramian)
        reversed_inverse = [
            reversal * coefficient * reversal
            for coefficient in inverse
        ]
        scalar = scalar_series_quotient(
            [coefficient[0, 0] for coefficient in reversed_inverse],
            [coefficient[0, 0] for coefficient in gramian],
        )
        residual = sp.simplify(
            reversed_inverse[jet]
            - sum(
                (
                    scalar[source_degree]
                    * gramian[jet - source_degree]
                    for source_degree in range(jet + 1)
                ),
                sp.zeros(dimension),
            )
        )
        equations = [
            residual[row, column]
            for row in range(dimension)
            for column in range(row, dimension)
            if residual[row, column] != 0
        ]
        solutions = sp.solve(
            equations,
            variables,
            dict=True,
            simplify=False,
        )
        if len(solutions) != 1:
            raise RuntimeError(
                f"self-dual defect jet {jet} was not uniquely determined"
            )
        defect = sp.simplify(trial.subs(solutions[0]))
    return defect


def real_circular_normal_direction(
    dimension: int,
    mode: int,
) -> sp.Matrix:
    """Return the real Riesz representative of one disk-normal support mode."""

    if not 2 <= mode <= dimension:
        raise ValueError("a circular normal mode must lie between 2 and p")
    weights = [
        1 / sp.sqrt(2) if index in (0, dimension - 1) else sp.Integer(1)
        for index in range(dimension)
    ]
    direction = sp.zeros(dimension)
    for row in range(dimension):
        for column in range(dimension):
            grade = column - row - 1
            if grade in (mode, -mode):
                direction[row, column] = weights[row] * weights[column]
    return direction


def endpoint_condition_coefficients(
    operator: Sequence[sp.Matrix],
    defect: sp.Matrix,
    epsilon: sp.Symbol,
    order: int,
) -> tuple[list[sp.Expr], list[sp.Expr], list[sp.Expr]]:
    """Return the two Stein endpoint series and their quotient."""

    gramian = stein_gramian_series(
        operator,
        defect,
        epsilon,
        order,
    )
    lower = simple_diagonal_eigenvalue_coefficients(
        gramian,
        endpoint=0,
        base_eigenvalue=1,
    )
    upper = simple_diagonal_eigenvalue_coefficients(
        gramian,
        endpoint=gramian[0].rows - 1,
        base_eigenvalue=4,
    )
    ratio = [sp.Integer(4), *[sp.Integer(0) for _ in range(order)]]
    for degree in range(1, order + 1):
        ratio[degree] = sp.expand(
            upper[degree]
            - sum(
                lower[source_degree] * ratio[degree - source_degree]
                for source_degree in range(1, degree + 1)
            )
        )
    return lower, upper, ratio
