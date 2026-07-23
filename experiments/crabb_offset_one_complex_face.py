#!/usr/bin/env python3
"""Audit the imaginary phase of the universal offset-one Stein face.

The real phase is handled by ``crabb_offset_one_face.py``.  Conjugation
symmetry kills the mixed real/imaginary term in the amplitude Hessian,
so it remains to prove the same coefficient for the pure-imaginary
phase.  This checker performs that calculation with exact rational
real/imaginary matrix pairs; it never converts Gaussian rationals to
floating point.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import (
    Matrix,
    Series,
    axis_metric_and_defect,
    direct_map_coefficients,
    matrix_add,
    matrix_multiply,
    matrix_scale,
    matrix_transpose,
    monomial,
    one,
    reverse_matrix,
    zero,
    zero_matrix,
)
from crabb_offset_one_face import (
    STABLE_MINIMUM_LENGTH,
    constant_crabb_part,
    stable_forcing_diagonal,
)
from exact_complex_matrix_series import (
    ComplexAmplitudeMatrix,
    ComplexMatrix,
    ComplexVector,
    amplitude_add,
    amplitude_power,
    amplitude_scale,
    squared_modulus,
    stein_gramian_expansion as complex_stein_gramian_expansion,
    zero_complex_matrix as complex_zero_matrix,
)

SERIES_ORDER = 3


@dataclass(frozen=True)
class ComplexOffsetOneRecord:
    dimension: int
    length: int
    defect_tangent: tuple[tuple[int, str, str], ...]
    condition_coefficients: tuple[str, str, str]
    lower_endpoint_metric_coefficient: str
    upper_endpoint_metric_coefficient: str
    endpoint_linear_coupling_valuation: int
    stationarity_residual_valuation: int | None
    quadratic_stein_forcing_diagonal: tuple[str, ...]
    weighted_forcing_sum: str
    stable_forcing_pattern: bool


def imaginary_coefficients(length: int) -> list[int]:
    """Return coefficients for ``u_1=i, u_(L-1)=-i``."""

    if length < 3:
        raise ValueError("the imaginary phase needs two distinct offsets")
    coefficients = [0] * (length - 1)
    coefficients[0] = 1
    coefficients[-1] = -1
    return coefficients


def imaginary_operator_expansion(
    dimension: int,
    order: int,
) -> tuple[ComplexAmplitudeMatrix, list[int]]:
    """Return the exact operator jet for the pure-imaginary phase."""

    length = dimension - 1
    crabb = zero_matrix(dimension, dimension, order)
    crabb[0][1] = 2 * one(order)
    for column in range(2, dimension):
        crabb[column - 1][column] = one(order)

    coefficients = imaginary_coefficients(length)
    imaginary_tangent = zero_matrix(dimension, dimension, order)
    for column in range(2, dimension):
        value = 2 * coefficients[column - 2]
        imaginary_tangent[0][column] = (
            Series.constant(value, order)
        )
        imaginary_tangent[length][column] = (
            Series.constant(-value, order)
        )

    parameter = monomial(1, order)
    base = matrix_add(
        crabb,
        matrix_scale(parameter, reverse_matrix(crabb)),
    )
    # If E=iF, then E+cJ conjugate(E) J=i(F-cJFJ).
    linear_imaginary = matrix_add(
        imaginary_tangent,
        matrix_scale(
            -parameter,
            reverse_matrix(imaginary_tangent),
        ),
    )
    real_zero = zero_matrix(dimension, dimension, order)
    pencil: ComplexAmplitudeMatrix = (
        (base, real_zero),
        (real_zero, linear_imaginary),
        complex_zero_matrix(dimension, dimension, order),
    )

    result: ComplexAmplitudeMatrix = (
        complex_zero_matrix(dimension, dimension, order),
        complex_zero_matrix(dimension, dimension, order),
        complex_zero_matrix(dimension, dimension, order),
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


def imaginary_coordinate_tangent(
    dimension: int,
    coefficients: list[int],
    order: int,
) -> ComplexMatrix:
    """Return the Hermitian coordinate-Gramian tangent."""

    length = dimension - 1
    imaginary_toeplitz = zero_matrix(dimension, dimension, order)
    for offset, coefficient in enumerate(coefficients, start=1):
        for row in range(length - offset):
            value = Series.constant(coefficient, order)
            imaginary_toeplitz[row][row + offset] = value
            imaginary_toeplitz[row + offset][row] = -value

    shift = zero_matrix(dimension, dimension, order)
    for column in range(1, dimension):
        shift[column - 1][column] = one(order)
    imaginary = matrix_add(
        imaginary_toeplitz,
        matrix_multiply(
            matrix_transpose(shift),
            matrix_multiply(imaginary_toeplitz, shift),
        ),
    )
    return zero_matrix(dimension, dimension, order), imaginary


def imaginary_defect_tangent(dimension: int) -> ComplexVector:
    """Return the universal imaginary-phase defect jet."""

    length = dimension - 1
    real = [zero(SERIES_ORDER) for _ in range(dimension)]
    imaginary = [zero(SERIES_ORDER) for _ in range(dimension)]
    imaginary[1] -= 2
    imaginary[length - 1] += 2
    imaginary[3] += 8 * monomial(1, SERIES_ORDER)
    return real, imaginary


def stein_gramian_expansion(
    operator: ComplexAmplitudeMatrix,
    base_defect: list[Series],
    defect_tangent: ComplexVector,
) -> ComplexAmplitudeMatrix:
    """Solve the complex Stein equation with a real axis defect."""

    dimension = len(operator[0][0])
    base: ComplexVector = (
        base_defect,
        [zero(SERIES_ORDER) for _ in range(dimension)],
    )
    return complex_stein_gramian_expansion(
        operator,
        base,
        defect_tangent,
    )


def endpoint_condition_hessian(
    operator: ComplexAmplitudeMatrix,
    axis_metric: Matrix,
    base_defect: list[Series],
    coordinate_diagonal: list[Fraction],
    coordinate_tangent: ComplexMatrix,
    defect_tangent: ComplexVector,
) -> tuple[Series, ComplexMatrix]:
    """Return the imaginary-phase amplitude Hessian."""

    dimension = len(operator[0][0])
    gramian = stein_gramian_expansion(
        operator,
        base_defect,
        defect_tangent,
    )
    if gramian[0][1] != zero_matrix(
        dimension,
        dimension,
        SERIES_ORDER,
    ):
        raise AssertionError("the axis metric acquired an imaginary part")
    if any(
        (
            gramian[0][0][row][column]
            - axis_metric[row][column]
        ).valuation()
        < SERIES_ORDER
        for row in range(dimension)
        for column in range(dimension)
    ):
        raise AssertionError("the complex axis metric drifted")

    axis_eigenvalues = [
        gramian[0][0][index][index] / coordinate_diagonal[index]
        for index in range(dimension)
    ]

    def variations(endpoint: int) -> tuple[Series, Series]:
        eigenvalue = axis_eigenvalues[endpoint]
        first_pair = (
            gramian[1][0][endpoint][endpoint]
            - eigenvalue * coordinate_tangent[0][endpoint][endpoint],
            gramian[1][1][endpoint][endpoint]
            - eigenvalue * coordinate_tangent[1][endpoint][endpoint],
        )
        if first_pair[1].valuation() < SERIES_ORDER:
            raise AssertionError("a Hermitian diagonal became imaginary")
        first = first_pair[0] / coordinate_diagonal[endpoint]
        second = (
            gramian[2][0][endpoint][endpoint]
            / coordinate_diagonal[endpoint]
            - first
            * coordinate_tangent[0][endpoint][endpoint]
            / coordinate_diagonal[endpoint]
        )
        for index in range(dimension):
            if index == endpoint:
                continue
            coupling = (
                gramian[1][0][endpoint][index]
                - eigenvalue
                * coordinate_tangent[0][endpoint][index],
                gramian[1][1][endpoint][index]
                - eigenvalue
                * coordinate_tangent[1][endpoint][index],
            )
            denominator = (
                gramian[0][0][index][index]
                - eigenvalue * coordinate_diagonal[index]
            )
            second -= (
                squared_modulus(coupling)
                / (coordinate_diagonal[endpoint] * denominator)
            )
        return first, second

    lower_linear, lower_quadratic = variations(0)
    upper_linear, upper_quadratic = variations(dimension - 1)
    condition = (
        upper_quadratic / axis_eigenvalues[0]
        - upper_linear
        * lower_linear
        / axis_eigenvalues[0] ** 2
        + axis_eigenvalues[-1]
        * lower_linear**2
        / axis_eigenvalues[0] ** 3
        - axis_eigenvalues[-1]
        * lower_quadratic
        / axis_eigenvalues[0] ** 2
    )
    return condition, gramian[2]


def make_record(
    dimension: int,
    audit_stationarity: bool,
) -> ComplexOffsetOneRecord:
    """Construct and validate one imaginary-phase certificate."""

    if dimension < 4:
        raise ValueError("the imaginary phase requires dimension at least four")

    operator, coefficients = imaginary_operator_expansion(
        dimension,
        SERIES_ORDER,
    )
    axis_metric, base_defect, coordinate_diagonal = (
        axis_metric_and_defect(
            dimension,
            operator[0][0],
            SERIES_ORDER,
        )
    )
    coordinate_tangent = imaginary_coordinate_tangent(
        dimension,
        coefficients,
        SERIES_ORDER,
    )
    defect_tangent = imaginary_defect_tangent(dimension)
    condition, quadratic_metric = endpoint_condition_hessian(
        operator,
        axis_metric,
        base_defect,
        coordinate_diagonal,
        coordinate_tangent,
        defect_tangent,
    )
    gramian = stein_gramian_expansion(
        operator,
        base_defect,
        defect_tangent,
    )
    crabb = constant_crabb_part(
        (
            operator[0][0],
            zero_matrix(dimension, dimension, SERIES_ORDER),
            zero_matrix(dimension, dimension, SERIES_ORDER),
        )
    )
    quadratic_forcing = matrix_add(
        gramian[2][0],
        matrix_scale(
            -1,
            matrix_multiply(
                matrix_transpose(crabb),
                matrix_multiply(gramian[2][0], crabb),
            ),
        ),
    )
    forcing_diagonal = tuple(
        quadratic_forcing[index][index].coefficient(2)
        for index in range(dimension)
    )
    weighted_forcing_sum = (
        4 * forcing_diagonal[0] + sum(forcing_diagonal[1:])
    )
    stable_pattern = (
        dimension - 1 < STABLE_MINIMUM_LENGTH
        or forcing_diagonal
        == stable_forcing_diagonal(dimension - 1)
    )

    lower_residuals = [
        (
            gramian[1][0][0][index]
            - 2 * coordinate_tangent[0][0][index],
            gramian[1][1][0][index]
            - 2 * coordinate_tangent[1][0][index],
        )
        for index in range(dimension)
    ]
    upper_residuals = [
        (
            gramian[1][0][-1][index]
            - 8 * coordinate_tangent[0][-1][index],
            gramian[1][1][-1][index]
            - 8 * coordinate_tangent[1][-1][index],
        )
        for index in range(dimension)
    ]
    endpoint_coupling_valuation = min(
        part.valuation()
        for residual in (*lower_residuals, *upper_residuals)
        for part in residual
    )

    stationarity_residual_valuation = None
    if audit_stationarity:
        gradient = []
        for index in range(1, dimension):
            unit = [zero(SERIES_ORDER) for _ in range(dimension)]
            unit[index] = one(SERIES_ORDER)
            positive = (
                defect_tangent[0].copy(),
                [
                    defect_tangent[1][entry] + unit[entry]
                    for entry in range(dimension)
                ],
            )
            negative = (
                defect_tangent[0].copy(),
                [
                    defect_tangent[1][entry] - unit[entry]
                    for entry in range(dimension)
                ],
            )
            positive_value = endpoint_condition_hessian(
                operator,
                axis_metric,
                base_defect,
                coordinate_diagonal,
                coordinate_tangent,
                positive,
            )[0]
            negative_value = endpoint_condition_hessian(
                operator,
                axis_metric,
                base_defect,
                coordinate_diagonal,
                coordinate_tangent,
                negative,
            )[0]
            gradient.append((positive_value - negative_value) / 2)
        stationarity_residual_valuation = min(
            entry.valuation()
            for entry in gradient
        )

    coefficients_tuple = tuple(
        condition.coefficient(degree)
        for degree in range(SERIES_ORDER)
    )
    lower_coefficient = quadratic_metric[0][0][0].coefficient(2)
    upper_coefficient = quadratic_metric[0][-1][-1].coefficient(2)
    if (
        coefficients_tuple != (0, 0, -64)
        or lower_coefficient != 48
        or upper_coefficient != 128
        or endpoint_coupling_valuation < 2
        or weighted_forcing_sum != 128
        or not stable_pattern
        or (
            stationarity_residual_valuation is not None
            and stationarity_residual_valuation < 2
        )
    ):
        raise AssertionError("the imaginary offset-one certificate failed")

    return ComplexOffsetOneRecord(
        dimension=dimension,
        length=dimension - 1,
        defect_tangent=tuple(
            (
                index,
                str(entry.coefficient(0)),
                str(entry.coefficient(1)),
            )
            for index, entry in enumerate(defect_tangent[1])
            if entry.valuation() < 2
        ),
        condition_coefficients=tuple(
            str(value)
            for value in coefficients_tuple
        ),
        lower_endpoint_metric_coefficient=str(lower_coefficient),
        upper_endpoint_metric_coefficient=str(upper_coefficient),
        endpoint_linear_coupling_valuation=(
            endpoint_coupling_valuation
        ),
        stationarity_residual_valuation=(
            stationarity_residual_valuation
        ),
        quadratic_stein_forcing_diagonal=tuple(
            str(value)
            for value in forcing_diagonal
        ),
        weighted_forcing_sum=str(weighted_forcing_sum),
        stable_forcing_pattern=stable_pattern,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=4)
    parser.add_argument("--maximum-size", type=int, default=12)
    parser.add_argument(
        "--stationarity-maximum-size",
        type=int,
        default=7,
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact dimension grid and optionally persist JSONL."""

    args = parse_args()
    if args.minimum_size < 4 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 4 <= minimum <= maximum")

    records = [
        make_record(
            dimension,
            dimension <= args.stationarity_maximum_size,
        )
        for dimension in range(
            args.minimum_size,
            args.maximum_size + 1,
        )
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
