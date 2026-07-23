"""Exact complex matrix algebra over rational truncated power series.

``exact_truncated_series.Series`` deliberately stores rational
coefficients.  For the complex Stein checkers, represent a matrix by
the pair ``(real_part, imaginary_part)`` instead of weakening that
invariant or converting Gaussian rationals to floating point.

The helpers here mirror only the degree-two amplitude algebra used by
the research checkers and keep the real matrix primitives in
``crabb_palindromic_elliptic_hessian`` as the single implementation.
"""

from __future__ import annotations

from fractions import Fraction
from typing import TypeAlias

from crabb_palindromic_elliptic_hessian import (
    Matrix,
    Series,
    matrix_add,
    matrix_multiply,
    matrix_outer,
    matrix_scale,
    matrix_transpose,
    one,
    zero_matrix,
)


ComplexMatrix: TypeAlias = tuple[Matrix, Matrix]
ComplexAmplitudeMatrix: TypeAlias = tuple[
    ComplexMatrix,
    ComplexMatrix,
    ComplexMatrix,
]
ComplexVector: TypeAlias = tuple[list[Series], list[Series]]


def zero_complex_matrix(
    rows: int,
    columns: int,
    order: int,
) -> ComplexMatrix:
    """Return a complex zero matrix."""

    return (
        zero_matrix(rows, columns, order),
        zero_matrix(rows, columns, order),
    )


def add(
    left: ComplexMatrix,
    right: ComplexMatrix,
) -> ComplexMatrix:
    """Add exact complex matrices."""

    return (
        matrix_add(left[0], right[0]),
        matrix_add(left[1], right[1]),
    )


def scale(
    scalar: int | Fraction | Series,
    matrix: ComplexMatrix,
) -> ComplexMatrix:
    """Scale an exact complex matrix by a real series."""

    return (
        matrix_scale(scalar, matrix[0]),
        matrix_scale(scalar, matrix[1]),
    )


def multiply(
    left: ComplexMatrix,
    right: ComplexMatrix,
) -> ComplexMatrix:
    """Multiply exact complex matrices."""

    return (
        matrix_add(
            matrix_multiply(left[0], right[0]),
            matrix_scale(
                -1,
                matrix_multiply(left[1], right[1]),
            ),
        ),
        matrix_add(
            matrix_multiply(left[0], right[1]),
            matrix_multiply(left[1], right[0]),
        ),
    )


def adjoint(matrix: ComplexMatrix) -> ComplexMatrix:
    """Return the conjugate transpose."""

    return (
        matrix_transpose(matrix[0]),
        matrix_scale(-1, matrix_transpose(matrix[1])),
    )


def identity(dimension: int, order: int) -> ComplexMatrix:
    """Return the exact complex identity."""

    real = zero_matrix(dimension, dimension, order)
    for index in range(dimension):
        real[index][index] = one(order)
    return real, zero_matrix(dimension, dimension, order)


def amplitude_add(
    left: ComplexAmplitudeMatrix,
    right: ComplexAmplitudeMatrix,
) -> ComplexAmplitudeMatrix:
    """Add complex matrices modulo amplitude degree three."""

    return tuple(
        add(left[degree], right[degree])
        for degree in range(3)
    )  # type: ignore[return-value]


def amplitude_scale(
    scalar: int | Fraction | Series,
    matrix: ComplexAmplitudeMatrix,
) -> ComplexAmplitudeMatrix:
    """Scale a complex amplitude expansion."""

    return tuple(
        scale(scalar, matrix[degree])
        for degree in range(3)
    )  # type: ignore[return-value]


def amplitude_multiply(
    left: ComplexAmplitudeMatrix,
    right: ComplexAmplitudeMatrix,
) -> ComplexAmplitudeMatrix:
    """Multiply complex matrices modulo amplitude degree three."""

    constant = multiply(left[0], right[0])
    linear = add(
        multiply(left[0], right[1]),
        multiply(left[1], right[0]),
    )
    quadratic = add(
        add(
            multiply(left[0], right[2]),
            multiply(left[1], right[1]),
        ),
        multiply(left[2], right[0]),
    )
    return constant, linear, quadratic


def amplitude_power(
    matrix: ComplexAmplitudeMatrix,
    exponent: int,
) -> ComplexAmplitudeMatrix:
    """Raise a complex amplitude expansion to a nonnegative power."""

    dimension = len(matrix[0][0])
    order = matrix[0][0][0][0].order
    result = (
        identity(dimension, order),
        zero_complex_matrix(dimension, dimension, order),
        zero_complex_matrix(dimension, dimension, order),
    )
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = amplitude_multiply(result, base)
        base = amplitude_multiply(base, base)
        remaining //= 2
    return result


def vector_cross(
    left: ComplexVector,
    right: ComplexVector,
) -> ComplexMatrix:
    """Return ``left right*``."""

    left_real, left_imaginary = left
    right_real, right_imaginary = right
    return (
        matrix_add(
            matrix_outer(left_real, right_real),
            matrix_outer(left_imaginary, right_imaginary),
        ),
        matrix_add(
            matrix_outer(left_imaginary, right_real),
            matrix_scale(
                -1,
                matrix_outer(left_real, right_imaginary),
            ),
        ),
    )


def vector_outer(vector: ComplexVector) -> ComplexMatrix:
    """Return ``vector vector*``."""

    return vector_cross(vector, vector)


def defect_forcing(
    base: ComplexVector,
    tangent: ComplexVector,
) -> ComplexAmplitudeMatrix:
    """Return ``(base+a tangent)(base+a tangent)*`` through degree two."""

    return (
        vector_outer(base),
        add(vector_cross(base, tangent), vector_cross(tangent, base)),
        vector_outer(tangent),
    )


def stein_gramian_expansion(
    operator: ComplexAmplitudeMatrix,
    base_defect: ComplexVector,
    defect_tangent: ComplexVector,
) -> ComplexAmplitudeMatrix:
    """Solve the complex Stein equation modulo amplitude degree three."""

    dimension = len(operator[0][0])
    order = operator[0][0][0][0].order
    forcing = defect_forcing(base_defect, defect_tangent)
    gramian: ComplexAmplitudeMatrix = (
        zero_complex_matrix(dimension, dimension, order),
        zero_complex_matrix(dimension, dimension, order),
        zero_complex_matrix(dimension, dimension, order),
    )
    operator_adjoint = tuple(
        adjoint(operator[degree])
        for degree in range(3)
    )
    for _ in range(12 * dimension * order):
        updated = amplitude_add(
            forcing,
            amplitude_multiply(
                operator_adjoint,  # type: ignore[arg-type]
                amplitude_multiply(gramian, operator),
            ),
        )
        if gramian == updated:
            return updated
        gramian = updated
    raise RuntimeError("the complex formal Stein iteration did not stabilize")


def squared_modulus(value: tuple[Series, Series]) -> Series:
    """Return the exact squared modulus of a complex series."""

    return value[0] ** 2 + value[1] ** 2
