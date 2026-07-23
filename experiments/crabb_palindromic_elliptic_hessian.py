#!/usr/bin/env python3
"""Exact finite-size audit of the palindromic elliptic Hessian.

This checker expands the locally optimized rank-one Stein condition
square to second order in an equality amplitude ``a`` and exactly in a
truncated rational series in the ellipse parameter ``c``.

For every tested phase-one coefficient supported at offsets ``k`` and
``L-k``, it verifies that the optimized amplitude Hessian has first
term

    -64 c^(2k).

The computation is exact in every coefficient below the declared
audit cutoff.  It is finite-size evidence for A84, not an all-size
proof of its recurrence or a uniform analytic remainder theorem.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import cache
import json
from pathlib import Path
from typing import TypeAlias

from exact_truncated_series import Series


Matrix: TypeAlias = list[list[Series]]
AmplitudeMatrix: TypeAlias = tuple[Matrix, Matrix, Matrix]
ConditionEvaluator: TypeAlias = Callable[
    [list[Series]],
    tuple[Series, Matrix],
]

DEFAULT_SERIES_ORDER = 10


@dataclass(frozen=True)
class HessianRecord:
    dimension: int
    first_offset: int
    target_degree: int
    first_nonzero_degree: int
    leading_coefficient: int
    predicted_coefficient: int
    lower_coefficients_vanish: bool
    lower_endpoint_metric_coefficient: str
    upper_endpoint_metric_coefficient: str
    optimized_defect_prefix: tuple[tuple[str, ...], ...]


def zero(order: int) -> Series:
    """Return the zero series."""

    return Series.constant(0, order)


def one(order: int) -> Series:
    """Return the unit series."""

    return Series.constant(1, order)


def monomial(degree: int, order: int, coefficient: int = 1) -> Series:
    """Return a rational monomial."""

    return Series.monomial(degree, order, coefficient)


def zero_matrix(rows: int, columns: int, order: int) -> Matrix:
    """Return a rectangular zero matrix."""

    return [[zero(order) for _ in range(columns)] for _ in range(rows)]


def identity_matrix(dimension: int, order: int) -> Matrix:
    """Return an identity matrix."""

    result = zero_matrix(dimension, dimension, order)
    for index in range(dimension):
        result[index][index] = one(order)
    return result


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    """Add two matrices."""

    return [
        [
            left[row][column] + right[row][column]
            for column in range(len(left[0]))
        ]
        for row in range(len(left))
    ]


def matrix_scale(scalar: int | Fraction | Series, matrix: Matrix) -> Matrix:
    """Scale a matrix."""

    return [
        [scalar * entry for entry in row]
        for row in matrix
    ]


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    """Multiply two matrices over the truncated-series ring."""

    rows = len(left)
    inner = len(right)
    columns = len(right[0])
    order = left[0][0].order
    result = zero_matrix(rows, columns, order)
    for row in range(rows):
        for index in range(inner):
            if left[row][index].valuation() == order:
                continue
            for column in range(columns):
                result[row][column] += (
                    left[row][index] * right[index][column]
                )
    return result


def matrix_transpose(matrix: Matrix) -> Matrix:
    """Transpose a matrix."""

    return [list(row) for row in zip(*matrix, strict=True)]


def matrix_outer(left: list[Series], right: list[Series]) -> Matrix:
    """Form a rank-one outer product."""

    return [
        [left_entry * right_entry for right_entry in right]
        for left_entry in left
    ]


def diagonal_matrix(diagonal: list[Series]) -> Matrix:
    """Construct a diagonal matrix."""

    result = zero_matrix(len(diagonal), len(diagonal), diagonal[0].order)
    for index, value in enumerate(diagonal):
        result[index][index] = value
    return result


def matrices_equal(left: Matrix, right: Matrix) -> bool:
    """Test exact equality of two truncated matrices."""

    return all(
        left[row][column] == right[row][column]
        for row in range(len(left))
        for column in range(len(left[0]))
    )


def amplitude_add(
    left: AmplitudeMatrix,
    right: AmplitudeMatrix,
) -> AmplitudeMatrix:
    """Add matrices truncated after amplitude degree two."""

    return tuple(
        matrix_add(left[degree], right[degree])
        for degree in range(3)
    )  # type: ignore[return-value]


def amplitude_scale(
    scalar: int | Fraction | Series,
    matrix: AmplitudeMatrix,
) -> AmplitudeMatrix:
    """Scale an amplitude-truncated matrix."""

    return tuple(
        matrix_scale(scalar, matrix[degree])
        for degree in range(3)
    )  # type: ignore[return-value]


def amplitude_multiply(
    left: AmplitudeMatrix,
    right: AmplitudeMatrix,
) -> AmplitudeMatrix:
    """Multiply matrices modulo amplitude degree three."""

    constant = matrix_multiply(left[0], right[0])
    linear = matrix_add(
        matrix_multiply(left[0], right[1]),
        matrix_multiply(left[1], right[0]),
    )
    quadratic = matrix_add(
        matrix_add(
            matrix_multiply(left[0], right[2]),
            matrix_multiply(left[1], right[1]),
        ),
        matrix_multiply(left[2], right[0]),
    )
    return constant, linear, quadratic


def amplitude_power(matrix: AmplitudeMatrix, exponent: int) -> AmplitudeMatrix:
    """Raise a matrix to a nonnegative power modulo ``a^3``."""

    dimension = len(matrix[0])
    order = matrix[0][0][0].order
    result = (
        identity_matrix(dimension, order),
        zero_matrix(dimension, dimension, order),
        zero_matrix(dimension, dimension, order),
    )
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = amplitude_multiply(result, base)
        base = amplitude_multiply(base, base)
        remaining //= 2
    return result


def scalar_convolution(
    coefficients: list[Series],
    degree: int,
    derivative_weights: bool = False,
) -> Series:
    """Return one coefficient of ``h^2`` or ``(h+2xh')^2``."""

    order = coefficients[0].order
    total = zero(order)
    for left_index in range(degree + 1):
        right_index = degree - left_index
        if (
            left_index >= len(coefficients)
            or right_index >= len(coefficients)
        ):
            continue
        factor = 1
        if derivative_weights:
            factor = (2 * left_index + 1) * (2 * right_index + 1)
        total += (
            factor
            * coefficients[left_index]
            * coefficients[right_index]
        )
    return total


def theta_data(order: int) -> tuple[Series, Series, Series]:
    """Return ``k``, ``alpha``, and the inverse-map linear coefficient."""

    theta_two_factor = zero(order)
    index = 0
    while 2 * index * (index + 1) < order:
        theta_two_factor += monomial(
            2 * index * (index + 1),
            order,
        )
        index += 1

    theta_three = one(order)
    index = 1
    while 2 * index * index < order:
        theta_three += 2 * monomial(2 * index * index, order)
        index += 1

    modulus = (
        4
        * monomial(1, order)
        * theta_two_factor**2
        / theta_three**2
    )
    period_factor = one(order) / theta_three**2
    linear_coefficient = one(order) / (
        theta_two_factor * theta_three
    )
    return modulus, period_factor, linear_coefficient


@cache
def inverse_map_coefficients(
    maximum_degree: int,
    work_order: int,
) -> tuple[Series, ...]:
    """Solve L125's inverse-map ODE with an explicit internal guard."""

    modulus, period_factor, linear_coefficient = theta_data(work_order)
    coefficients = [linear_coefficient]
    for degree in range(1, maximum_degree + 1):
        padded = [*coefficients, zero(work_order)]
        numerator = (
            (one(work_order) + modulus**2)
            * scalar_convolution(
                coefficients,
                degree - 1,
                derivative_weights=True,
            )
            - modulus
            * (
                scalar_convolution(
                    coefficients,
                    degree - 2,
                    derivative_weights=True,
                )
                if degree >= 2
                else zero(work_order)
            )
            - period_factor**2
            * scalar_convolution(coefficients, degree - 1)
            - modulus
            * scalar_convolution(
                padded,
                degree,
                derivative_weights=True,
            )
        )
        denominator = (
            2
            * modulus
            * linear_coefficient
            * (2 * degree + 1)
        )
        coefficients.append(numerator / denominator)
    return tuple(coefficients)


def polynomial_convolution(
    left: list[Series],
    right: list[Series],
    maximum_degree: int,
) -> list[Series]:
    """Multiply scalar polynomials in the functional-calculus variable."""

    order = left[0].order
    result = [zero(order) for _ in range(maximum_degree + 1)]
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            degree = left_degree + right_degree
            if degree <= maximum_degree:
                result[degree] += left_coefficient * right_coefficient
    return result


def polynomial_power(
    polynomial: list[Series],
    exponent: int,
    maximum_degree: int,
) -> list[Series]:
    """Raise a scalar polynomial to a nonnegative power."""

    order = polynomial[0].order
    result = [one(order), *(
        zero(order) for _ in range(maximum_degree)
    )]
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = polynomial_convolution(
                result,
                base,
                maximum_degree,
            )
        base = polynomial_convolution(base, base, maximum_degree)
        remaining //= 2
    return result


@cache
def direct_map_coefficients(
    maximum_degree: int,
    output_order: int,
) -> tuple[Series, ...]:
    """Revert the inverse series without contaminating terminal c-jets."""

    # The inverse-map recurrence divides by a valuation-one modulus at
    # every scalar degree.  L125's guard calculation shows that this
    # work order preserves every requested output coefficient.
    work_order = output_order + 2 * maximum_degree + 4
    inverse_coefficients = inverse_map_coefficients(
        maximum_degree,
        work_order,
    )
    coefficients: list[Series] = []
    maximum_polynomial_degree = 2 * maximum_degree + 1
    for index in range(maximum_degree + 1):
        direct_polynomial = [
            zero(work_order)
            for _ in range(maximum_polynomial_degree + 1)
        ]
        for known_index, coefficient in enumerate(coefficients):
            direct_polynomial[2 * known_index + 1] = coefficient

        known_term = zero(work_order)
        for inverse_index in range(index + 1):
            power = polynomial_power(
                direct_polynomial,
                2 * inverse_index + 1,
                maximum_polynomial_degree,
            )
            known_term += (
                inverse_coefficients[inverse_index]
                * power[2 * index + 1]
            )

        if index == 0:
            coefficients.append(
                one(work_order) / inverse_coefficients[0]
            )
        else:
            coefficients.append(
                -known_term / inverse_coefficients[0]
            )

    truncated = tuple(
        Series.from_coefficients(
            coefficient.coefficients,
            output_order,
        )
        for coefficient in coefficients
    )
    for index, coefficient in enumerate(truncated):
        if any(
            coefficient.coefficient(degree)
            for degree in range(index)
        ):
            raise AssertionError(
                "the direct ellipse map violated L125's filtration"
            )
        if (
            index < output_order
            and coefficient.coefficient(index) != (-1) ** index
        ):
            raise AssertionError(
                "the direct ellipse map violated its Newton-edge sign"
            )
    return truncated


def reverse_matrix(matrix: Matrix) -> Matrix:
    """Conjugate a real matrix by coordinate reversal."""

    return [row[::-1] for row in matrix[::-1]]


def palindromic_direction(length: int, first_offset: int) -> list[int]:
    """Return a phase-one coefficient vector on one reversed pair."""

    coefficients = [0] * (length - 1)
    coefficients[first_offset - 1] = 1
    coefficients[length - first_offset - 1] = 1
    return coefficients


def operator_expansion(
    dimension: int,
    first_offset: int,
    order: int,
) -> tuple[AmplitudeMatrix, list[int]]:
    """Return ``T_0+aT_1+a^2T_2`` in L123 coefficient coordinates."""

    length = dimension - 1
    coefficients = palindromic_direction(length, first_offset)
    return operator_expansion_from_coefficients(
        dimension,
        coefficients,
        order,
    )


def operator_expansion_from_coefficients(
    dimension: int,
    coefficients: list[int],
    order: int,
) -> tuple[AmplitudeMatrix, list[int]]:
    """Return the operator jet for an arbitrary real coefficient vector."""

    length = dimension - 1
    if len(coefficients) != length - 1:
        raise ValueError("the coefficient vector has the wrong length")

    crabb = zero_matrix(dimension, dimension, order)
    crabb[0][1] = Series.constant(2, order)
    for column in range(2, dimension):
        crabb[column - 1][column] = one(order)

    tangent = zero_matrix(dimension, dimension, order)
    for column in range(2, dimension):
        value = 2 * coefficients[column - 2]
        tangent[0][column] = Series.constant(value, order)
        tangent[length][column] = Series.constant(-value, order)

    ellipse_parameter = monomial(1, order)
    base = matrix_add(
        crabb,
        matrix_scale(ellipse_parameter, reverse_matrix(crabb)),
    )
    linear = matrix_add(
        tangent,
        matrix_scale(ellipse_parameter, reverse_matrix(tangent)),
    )
    zero_block = zero_matrix(dimension, dimension, order)
    pencil = base, linear, zero_block

    result: AmplitudeMatrix = (
        zero_matrix(dimension, dimension, order),
        zero_matrix(dimension, dimension, order),
        zero_matrix(dimension, dimension, order),
    )
    for index, coefficient in enumerate(
        direct_map_coefficients(order - 1, order)
    ):
        if coefficient.valuation() == order:
            continue
        result = amplitude_add(
            result,
            amplitude_scale(
                coefficient,
                amplitude_power(pencil, 2 * index + 1),
            ),
        )
    return result, coefficients


def sech_series(exponent: int, order: int) -> Series:
    """Expand ``sech(exponent * ell)`` for ``c=exp(-ell)``."""

    if exponent == 0:
        return one(order)
    exponent = abs(exponent)
    return (
        2
        * monomial(exponent, order)
        / (one(order) + monomial(2 * exponent, order))
    )


def axis_metric_and_defect(
    dimension: int,
    operator: Matrix,
    audit_order: int,
) -> tuple[Matrix, list[Series], list[Fraction]]:
    """Construct the exact L117 axis metric and its rank-one defect."""

    order = operator[0][0].order
    length = dimension - 1
    periodized = []
    for index in range(dimension):
        value = zero(order)
        for alias in range(-order, order + 1):
            exponent = index + 2 * length * alias
            if abs(exponent) < order:
                value += sech_series(exponent, order)
        periodized.append(value)

    # Form S_m/c^m directly.  Dividing an already truncated S_m would
    # lose precisely the high-order endpoint coefficients under audit.
    divided_periodized = []
    for index in range(dimension):
        value = zero(order)
        for alias in range(-order, order + 1):
            exponent = abs(index + 2 * length * alias)
            shifted_degree = exponent - index
            if exponent == 0:
                if index == 0:
                    value += one(order)
            elif shifted_degree < order:
                value += (
                    2
                    * monomial(shifted_degree, order)
                    / (one(order) + monomial(2 * exponent, order))
                )
        divided_periodized.append(value)

    physical_diagonal = [
        value / periodized[0]
        for value in divided_periodized
    ]
    coordinate_diagonal = [
        Fraction(1, 2),
        *(Fraction(1) for _ in range(dimension - 2)),
        Fraction(1, 2),
    ]
    # Scale the metric by two so the defect starts at e_0.
    metric = diagonal_matrix(
        [
            2 * coordinate_diagonal[index] * physical_diagonal[index]
            for index in range(dimension)
        ]
    )
    defect_matrix = matrix_add(
        metric,
        matrix_scale(
            -1,
            matrix_multiply(
                matrix_transpose(operator),
                matrix_multiply(metric, operator),
            ),
        ),
    )
    defect = [zero(order) for _ in range(dimension)]
    defect[0] = defect_matrix[0][0].square_root_unit()
    for index in range(1, dimension):
        defect[index] = defect_matrix[index][0] / defect[0]

    rank_one_residual = matrix_add(
        defect_matrix,
        matrix_scale(-1, matrix_outer(defect, defect)),
    )
    if any(
        entry.valuation() < audit_order
        for row in rank_one_residual
        for entry in row
    ):
        raise AssertionError("the axis Stein defect did not remain rank one")
    return metric, defect, coordinate_diagonal


def coordinate_metric_tangent(
    dimension: int,
    coefficients: list[int],
    order: int,
) -> Matrix:
    """Return the derivative of ``K=H+R*HR`` along the Toeplitz curve."""

    length = dimension - 1
    toeplitz = zero_matrix(dimension, dimension, order)
    for offset, coefficient in enumerate(coefficients, start=1):
        for row in range(length - offset):
            toeplitz[row][row + offset] = Series.constant(
                coefficient,
                order,
            )
            toeplitz[row + offset][row] = Series.constant(
                coefficient,
                order,
            )

    shift = zero_matrix(dimension, dimension, order)
    for column in range(1, dimension):
        shift[column - 1][column] = one(order)
    return matrix_add(
        toeplitz,
        matrix_multiply(
            matrix_transpose(shift),
            matrix_multiply(toeplitz, shift),
        ),
    )


def stein_gramian_expansion(
    operator: AmplitudeMatrix,
    base_defect: list[Series],
    defect_tangent: list[Series],
) -> AmplitudeMatrix:
    """Solve the Stein equation modulo ``a^3`` by exact iteration."""

    dimension = len(operator[0])
    order = operator[0][0][0].order
    forcing: AmplitudeMatrix = (
        matrix_outer(base_defect, base_defect),
        matrix_add(
            matrix_outer(base_defect, defect_tangent),
            matrix_outer(defect_tangent, base_defect),
        ),
        matrix_outer(defect_tangent, defect_tangent),
    )
    gramian: AmplitudeMatrix = (
        zero_matrix(dimension, dimension, order),
        zero_matrix(dimension, dimension, order),
        zero_matrix(dimension, dimension, order),
    )
    transpose = tuple(
        matrix_transpose(operator[degree])
        for degree in range(3)
    )
    for _ in range(12 * dimension * order):
        updated = amplitude_add(
            forcing,
            amplitude_multiply(
                transpose,  # type: ignore[arg-type]
                amplitude_multiply(gramian, operator),
            ),
        )
        if all(
            matrices_equal(gramian[degree], updated[degree])
            for degree in range(3)
        ):
            return updated
        gramian = updated
    raise RuntimeError("the formal Stein iteration did not stabilize")


def endpoint_condition_hessian(
    operator: AmplitudeMatrix,
    axis_metric: Matrix,
    base_defect: list[Series],
    coordinate_diagonal: list[Fraction],
    coordinate_tangent: Matrix,
    defect_tangent: list[Series],
    audit_order: int,
) -> tuple[Series, Matrix]:
    """Return the amplitude Hessian of the endpoint condition ratio."""

    dimension = len(operator[0])
    gramian = stein_gramian_expansion(
        operator,
        base_defect,
        defect_tangent,
    )
    if any(
        (
            gramian[0][row][column]
            - axis_metric[row][column]
        ).valuation()
        < audit_order
        for row in range(dimension)
        for column in range(dimension)
    ):
        raise AssertionError("the regenerated axis metric drifted")

    axis_eigenvalues = [
        gramian[0][index][index] / coordinate_diagonal[index]
        for index in range(dimension)
    ]

    def eigenvalue_variations(endpoint: int) -> tuple[Series, Series]:
        eigenvalue = axis_eigenvalues[endpoint]
        first_variation = (
            gramian[1][endpoint][endpoint]
            - eigenvalue * coordinate_tangent[endpoint][endpoint]
        ) / coordinate_diagonal[endpoint]
        second_variation = (
            gramian[2][endpoint][endpoint]
            / coordinate_diagonal[endpoint]
            - first_variation
            * coordinate_tangent[endpoint][endpoint]
            / coordinate_diagonal[endpoint]
        )
        for index in range(dimension):
            if index == endpoint:
                continue
            coupling = (
                gramian[1][endpoint][index]
                - eigenvalue * coordinate_tangent[endpoint][index]
            )
            denominator = (
                gramian[0][index][index]
                - eigenvalue * coordinate_diagonal[index]
            )
            second_variation -= (
                coupling**2
                / (coordinate_diagonal[endpoint] * denominator)
            )
        return first_variation, second_variation

    lower_linear, lower_hessian = eigenvalue_variations(0)
    upper_linear, upper_hessian = eigenvalue_variations(dimension - 1)
    condition_hessian = (
        upper_hessian / axis_eigenvalues[0]
        - upper_linear
        * lower_linear
        / axis_eigenvalues[0] ** 2
        + axis_eigenvalues[-1]
        * lower_linear**2
        / axis_eigenvalues[0] ** 3
        - axis_eigenvalues[-1]
        * lower_hessian
        / axis_eigenvalues[0] ** 2
    )
    return condition_hessian, gramian[2]


def solve_series_system(
    matrix: list[list[Series]],
    right_hand_side: list[Series],
) -> list[Series]:
    """Solve a square system whose leading matrix is nonsingular."""

    dimension = len(right_hand_side)
    augmented = [
        [*matrix[row], right_hand_side[row]]
        for row in range(dimension)
    ]
    for column in range(dimension):
        pivot = next(
            (
                row
                for row in range(column, dimension)
                if augmented[row][column].valuation() == 0
            ),
            None,
        )
        if pivot is None:
            raise RuntimeError("the formal Hessian lost its leading pivot")
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        inverse = augmented[column][column].inverse_unit()
        augmented[column] = [
            entry * inverse
            for entry in augmented[column]
        ]
        for row in range(dimension):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][index]
                - factor * augmented[column][index]
                for index in range(dimension + 1)
            ]
    return [augmented[row][-1] for row in range(dimension)]


def optimized_hessian(
    dimension: int,
    first_offset: int,
    order: int,
) -> tuple[Series, list[Series], Matrix]:
    """Reconstruct and minimize the exact quadratic in the defect tangent."""

    operator, coefficients = operator_expansion(
        dimension,
        first_offset,
        order,
    )
    return optimized_hessian_from_operator(
        operator,
        coefficients,
        order,
    )


def defect_quadratic_data(
    operator: AmplitudeMatrix,
    coefficients: list[int],
    order: int,
) -> tuple[
    Series,
    list[Series],
    list[list[Series]],
    ConditionEvaluator,
]:
    """Reconstruct the exact quadratic in the free defect tangent."""

    dimension = len(operator[0])
    if len(coefficients) != dimension - 2:
        raise ValueError("the coefficient vector has the wrong length")

    audit_order = order
    axis_metric, base_defect, coordinate_diagonal = (
        axis_metric_and_defect(
            dimension,
            operator[0],
            audit_order,
        )
    )
    coordinate_tangent = coordinate_metric_tangent(
        dimension,
        coefficients,
        order,
    )

    zero_tangent = [zero(order) for _ in range(dimension)]

    def evaluate(tangent: list[Series]) -> tuple[Series, Matrix]:
        return endpoint_condition_hessian(
            operator,
            axis_metric,
            base_defect,
            coordinate_diagonal,
            coordinate_tangent,
            tangent,
            audit_order,
        )

    constant, _ = evaluate(zero_tangent)
    positive_values = []
    negative_values = []
    for variable in range(1, dimension):
        tangent = [zero(order) for _ in range(dimension)]
        tangent[variable] = one(order)
        positive_values.append(evaluate(tangent)[0])
        tangent[variable] = -one(order)
        negative_values.append(evaluate(tangent)[0])

    variable_count = dimension - 1
    linear = [zero(order) for _ in range(variable_count)]
    quadratic = [
        [zero(order) for _ in range(variable_count)]
        for _ in range(variable_count)
    ]
    for index in range(variable_count):
        linear[index] = (
            positive_values[index] - negative_values[index]
        ) / 2
        quadratic[index][index] = (
            positive_values[index] + negative_values[index]
        ) / 2 - constant

    for left in range(variable_count):
        for right in range(left + 1, variable_count):
            tangent = [zero(order) for _ in range(dimension)]
            tangent[left + 1] = one(order)
            tangent[right + 1] = one(order)
            value, _ = evaluate(tangent)
            mixed = (
                value
                - constant
                - linear[left]
                - linear[right]
                - quadratic[left][left]
                - quadratic[right][right]
            ) / 2
            quadratic[left][right] = mixed
            quadratic[right][left] = mixed

    return constant, linear, quadratic, evaluate


def optimized_hessian_from_operator(
    operator: AmplitudeMatrix,
    coefficients: list[int],
    order: int,
) -> tuple[Series, list[Series], Matrix]:
    """Optimize the defect tangent for a supplied exact operator jet."""

    constant, linear, quadratic, evaluate = defect_quadratic_data(
        operator,
        coefficients,
        order,
    )
    variable_count = len(linear)
    optimizer = solve_series_system(
        quadratic,
        [-entry / 2 for entry in linear],
    )
    optimized = constant
    optimized += sum(
        linear[index] * optimizer[index]
        for index in range(variable_count)
    )
    optimized += sum(
        optimizer[left]
        * quadratic[left][right]
        * optimizer[right]
        for left in range(variable_count)
        for right in range(variable_count)
    )
    _, optimized_metric_quadratic = evaluate(
        [zero(order), *optimizer]
    )
    return optimized, optimizer, optimized_metric_quadratic


def coefficient_string(value: Fraction) -> str:
    """Serialize a rational coefficient without losing exactness."""

    return str(value)


def make_record(
    dimension: int,
    first_offset: int,
    order: int,
) -> HessianRecord:
    """Run and validate one exact amplitude-Hessian audit."""

    target_degree = 2 * first_offset
    if target_degree >= order:
        raise ValueError(
            "increase the series order so the target lies below the audit cutoff"
        )
    optimized, defect_tangent, metric_quadratic = optimized_hessian(
        dimension,
        first_offset,
        order,
    )
    lower_vanishes = all(
        optimized.coefficient(degree) == 0
        for degree in range(target_degree)
    )
    leading = optimized.coefficient(target_degree)
    if not lower_vanishes or leading != -64:
        raise AssertionError(
            "the exact elliptic Hessian did not have the predicted first term"
        )

    prefix_order = first_offset + 1
    defect_prefix = tuple(
        tuple(
            coefficient_string(entry.coefficient(degree))
            for degree in range(prefix_order)
        )
        for entry in defect_tangent
    )
    return HessianRecord(
        dimension=dimension,
        first_offset=first_offset,
        target_degree=target_degree,
        first_nonzero_degree=optimized.valuation(),
        leading_coefficient=int(leading),
        predicted_coefficient=-64,
        lower_coefficients_vanish=lower_vanishes,
        lower_endpoint_metric_coefficient=coefficient_string(
            metric_quadratic[0][0].coefficient(target_degree)
        ),
        upper_endpoint_metric_coefficient=coefficient_string(
            metric_quadratic[-1][-1].coefficient(target_degree)
        ),
        optimized_defect_prefix=defect_prefix,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=7)
    parser.add_argument(
        "--series-order",
        type=int,
        default=DEFAULT_SERIES_ORDER,
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic exact audit and optionally persist JSONL."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    if args.series_order < 5:
        raise ValueError("the series order must be at least five")

    records = []
    for dimension in range(args.minimum_size, args.maximum_size + 1):
        length = dimension - 1
        for first_offset in range(1, length // 2 + 1):
            if 2 * first_offset >= args.series_order:
                continue
            record = make_record(
                dimension,
                first_offset,
                args.series_order,
            )
            records.append(record)
            print(json.dumps(asdict(record), sort_keys=True), flush=True)

    if args.output is not None:
        args.output.write_text(
            "".join(
                f"{json.dumps(asdict(record), sort_keys=True)}\n"
                for record in records
            ),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
