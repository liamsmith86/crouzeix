#!/usr/bin/env python3
"""Regenerate the noncentral corrected Faber--Blaschke Newton face.

This checker works entirely over truncated rational power series.  It
constructs the Weierstrass numerator of

    G_(a,c)(Psi_c(w))

through amplitude degree two, evaluates the resulting finite Blaschke
product on the exact phase-palindromic matrix pencil, and lifts its
simple top generalized singular pair.  For every noncentral
``1 <= k < L/2`` it tests, the first nonzero amplitude Hessian term is

    -64 c^(2k).

The implementation also checks the sparse first-face tangent and the
only possible endpoint alias.  The finite runs are regeneration tests
for L143's coefficient induction; they are not, by themselves, a
substitute for that all-size induction.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import TypeAlias

from crabb_palindromic_elliptic_hessian import (
    AmplitudeMatrix,
    Matrix,
    amplitude_multiply,
    coordinate_metric_tangent,
    diagonal_matrix,
    identity_matrix,
    inverse_map_coefficients,
    matrix_add,
    matrix_multiply,
    matrix_scale,
    matrix_transpose,
    operator_expansion_from_coefficients,
    zero_matrix,
)
from exact_truncated_series import Series


FractionPolynomial: TypeAlias = list[Fraction]
SeriesPolynomial: TypeAlias = list[Series]

DEFAULT_CASES = (
    (3, 1),
    (4, 1),
    (5, 2),
    (6, 2),
    (7, 3),
    (8, 3),
    (9, 4),
)


@dataclass(frozen=True)
class FormalFaberBlaschkeRecord:
    """One exact noncentral dual-Hessian audit."""

    length: int
    dimension: int
    grade: int
    target_degree: int
    first_nonzero_degree: int
    leading_coefficient: str
    predicted_coefficient: str
    lower_coefficients_vanish: bool
    linear_eigenvalue_vanishes: bool
    first_face_tangent_matches: bool
    endpoint_alias_degree: int
    endpoint_alias_active: bool
    first_singular_coupling_valuation: int
    direct_rayleigh_coefficient: str
    schur_coefficient: str


@dataclass(frozen=True)
class FormalDualCalculation:
    """Reusable exact dual Hessian for one or several grades."""

    operator: AmplitudeMatrix
    denominator: AmplitudeMatrix
    blaschke: AmplitudeMatrix
    coordinate_constant: Matrix
    coordinate_linear: Matrix
    linear_eigenvalue: Series
    first_coupling: tuple[Series, ...]
    direct_rayleigh: Series
    schur: Series
    quadratic_eigenvalue: Series


def zero(order: int) -> Series:
    """Return the zero series."""

    return Series.constant(0, order)


def one(order: int) -> Series:
    """Return the unit series."""

    return Series.constant(1, order)


def trim_series_polynomial(
    polynomial: SeriesPolynomial,
) -> SeriesPolynomial:
    """Remove trailing zero coefficients."""

    while (
        polynomial
        and polynomial[-1].valuation() == polynomial[-1].order
    ):
        polynomial.pop()
    return polynomial


def add_series_polynomials(
    left: SeriesPolynomial,
    right: SeriesPolynomial,
) -> SeriesPolynomial:
    """Add polynomials in the scalar variable with series coefficients."""

    order = (left or right)[0].order
    result = [zero(order) for _ in range(max(len(left), len(right)))]
    for index, coefficient in enumerate(left):
        result[index] += coefficient
    for index, coefficient in enumerate(right):
        result[index] += coefficient
    return trim_series_polynomial(result)


def scale_series_polynomial(
    scalar: int | Fraction | Series,
    polynomial: SeriesPolynomial,
) -> SeriesPolynomial:
    """Scale a series-coefficient polynomial."""

    return trim_series_polynomial(
        [scalar * coefficient for coefficient in polynomial]
    )


def multiply_series_polynomials(
    left: SeriesPolynomial,
    right: SeriesPolynomial,
) -> SeriesPolynomial:
    """Multiply scalar polynomials over the truncated-series ring."""

    if not left or not right:
        return []
    order = left[0].order
    result = [zero(order) for _ in range(len(left) + len(right) - 1)]
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] += (
                left_coefficient * right_coefficient
            )
    return trim_series_polynomial(result)


def trim_fraction_polynomial(
    polynomial: FractionPolynomial,
) -> FractionPolynomial:
    """Remove trailing zero rational coefficients."""

    while polynomial and not polynomial[-1]:
        polynomial.pop()
    return polynomial


def add_fraction_polynomials(
    left: FractionPolynomial,
    right: FractionPolynomial,
) -> FractionPolynomial:
    """Add ordinary rational polynomials."""

    result = [Fraction(0) for _ in range(max(len(left), len(right)))]
    for index, coefficient in enumerate(left):
        result[index] += coefficient
    for index, coefficient in enumerate(right):
        result[index] += coefficient
    return trim_fraction_polynomial(result)


def scale_fraction_polynomial(
    scalar: int | Fraction,
    polynomial: FractionPolynomial,
) -> FractionPolynomial:
    """Scale an ordinary rational polynomial."""

    return trim_fraction_polynomial(
        [Fraction(scalar) * coefficient for coefficient in polynomial]
    )


def multiply_fraction_polynomials(
    left: FractionPolynomial,
    right: FractionPolynomial,
) -> FractionPolynomial:
    """Multiply ordinary rational polynomials."""

    if not left or not right:
        return []
    result = [
        Fraction(0)
        for _ in range(len(left) + len(right) - 1)
    ]
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] += (
                left_coefficient * right_coefficient
            )
    return trim_fraction_polynomial(result)


def composed_dickson_polynomials(
    maximum_degree: int,
    series_order: int,
) -> tuple[SeriesPolynomial, ...]:
    """Return ``P_j(Psi_c(w))`` through the requested ``c`` order."""

    maximum_c_degree = series_order - 1
    work_order = series_order + 2 * maximum_c_degree + 4
    inverse_coefficients = inverse_map_coefficients(
        maximum_c_degree,
        work_order,
    )
    inverse_map = [
        zero(series_order)
        for _ in range(2 * maximum_c_degree + 2)
    ]
    for index, coefficient in enumerate(inverse_coefficients):
        inverse_map[2 * index + 1] = Series.from_coefficients(
            coefficient.coefficients[:series_order],
            series_order,
        )

    polynomials: list[SeriesPolynomial] = [
        [Series.constant(2, series_order)],
        inverse_map,
    ]
    parameter = Series.monomial(1, series_order)
    for _ in range(2, maximum_degree + 1):
        polynomials.append(
            add_series_polynomials(
                multiply_series_polynomials(
                    inverse_map,
                    polynomials[-1],
                ),
                scale_series_polynomial(
                    -parameter,
                    polynomials[-2],
                ),
            )
        )
    return tuple(polynomials)


def weierstrass_numerator_coefficients(
    length: int,
    grades: tuple[int, ...],
    maximum_c_degree: int,
) -> dict[tuple[int, int], FractionPolynomial]:
    """Prepare a corrected multi-grade factor through amplitude degree two."""

    series_order = maximum_c_degree + 1
    dickson = composed_dickson_polynomials(length, series_order)
    constant_part = dickson[length]
    linear_part: SeriesPolynomial = []
    for grade in grades:
        if not 1 <= grade <= length / 2:
            raise ValueError("grades must be distinct low representatives")
        if 2 * grade == length:
            grade_part = scale_series_polynomial(
                2,
                dickson[grade],
            )
            reflected_part = scale_series_polynomial(
                2 * Series.monomial(grade, series_order),
                dickson[grade],
            )
        else:
            grade_part = scale_series_polynomial(
                2,
                add_series_polynomials(
                    dickson[grade],
                    dickson[length - grade],
                ),
            )
            reflected_part = scale_series_polynomial(
                2 * Series.monomial(grade, series_order),
                dickson[length - grade],
            )
        linear_part = add_series_polynomials(
            linear_part,
            add_series_polynomials(
                grade_part,
                reflected_part,
            ),
        )

    prepared_input: dict[
        tuple[int, int],
        FractionPolynomial,
    ] = {}
    for amplitude_degree, polynomial in enumerate(
        (constant_part, linear_part)
    ):
        for c_degree in range(series_order):
            coefficient = trim_fraction_polynomial(
                [
                    entry.coefficient(c_degree)
                    for entry in polynomial
                ]
            )
            if coefficient:
                prepared_input[amplitude_degree, c_degree] = coefficient

    numerator: dict[tuple[int, int], FractionPolynomial] = {
        (0, 0): [Fraction(0)] * length + [Fraction(1)]
    }
    outer: dict[tuple[int, int], FractionPolynomial] = {
        (0, 0): [Fraction(1)]
    }
    for total_degree in range(1, maximum_c_degree + 3):
        for amplitude_degree in range(3):
            c_degree = total_degree - amplitude_degree
            if not 0 <= c_degree <= maximum_c_degree:
                continue
            remainder = prepared_input.get(
                (amplitude_degree, c_degree),
                [],
            )
            for left_amplitude in range(amplitude_degree + 1):
                for left_c in range(c_degree + 1):
                    left_index = (left_amplitude, left_c)
                    right_index = (
                        amplitude_degree - left_amplitude,
                        c_degree - left_c,
                    )
                    if (
                        left_index in ((0, 0), (amplitude_degree, c_degree))
                        or right_index == (0, 0)
                    ):
                        continue
                    remainder = add_fraction_polynomials(
                        remainder,
                        scale_fraction_polynomial(
                            -1,
                            multiply_fraction_polynomials(
                                numerator.get(left_index, []),
                                outer.get(right_index, []),
                            ),
                        ),
                    )

            low = trim_fraction_polynomial(remainder[:length])
            high = trim_fraction_polynomial(remainder[length:])
            if low:
                numerator[amplitude_degree, c_degree] = low
            if high:
                outer[amplitude_degree, c_degree] = high

    for amplitude_degree in range(3):
        for c_degree in range(series_order):
            reconstructed: FractionPolynomial = []
            for left_amplitude in range(amplitude_degree + 1):
                for left_c in range(c_degree + 1):
                    reconstructed = add_fraction_polynomials(
                        reconstructed,
                        multiply_fraction_polynomials(
                            numerator.get(
                                (left_amplitude, left_c),
                                [],
                            ),
                            outer.get(
                                (
                                    amplitude_degree - left_amplitude,
                                    c_degree - left_c,
                                ),
                                [],
                            ),
                        ),
                    )
            residual = add_fraction_polynomials(
                reconstructed,
                scale_fraction_polynomial(
                    -1,
                    prepared_input.get(
                        (amplitude_degree, c_degree),
                        [],
                    ),
                ),
            )
            if residual:
                raise AssertionError("Weierstrass preparation did not close")
    return numerator


def matrix_polynomial_coefficients(
    prepared: dict[tuple[int, int], FractionPolynomial],
    length: int,
    maximum_c_degree: int,
    reverse: bool,
) -> list[list[Series]]:
    """Convert prepared coefficients to numerator or reversed denominator."""

    values = [
        [
            [Fraction(0) for _ in range(maximum_c_degree + 1)]
            for _ in range(length + 1)
        ]
        for _ in range(3)
    ]
    for (amplitude_degree, c_degree), polynomial in prepared.items():
        for scalar_degree, coefficient in enumerate(polynomial):
            target_degree = (
                length - scalar_degree
                if reverse
                else scalar_degree
            )
            values[amplitude_degree][target_degree][c_degree] += (
                coefficient
            )
    return [
        [
            Series.from_coefficients(
                values[amplitude_degree][scalar_degree],
                maximum_c_degree + 1,
            )
            for scalar_degree in range(length + 1)
        ]
        for amplitude_degree in range(3)
    ]


def amplitude_matrix_powers(
    matrix: AmplitudeMatrix,
    maximum_degree: int,
) -> tuple[AmplitudeMatrix, ...]:
    """Build consecutive powers once for both rational factors."""

    dimension = len(matrix[0])
    order = matrix[0][0][0].order
    powers: list[AmplitudeMatrix] = [
        (
            identity_matrix(dimension, order),
            zero_matrix(dimension, dimension, order),
            zero_matrix(dimension, dimension, order),
        )
    ]
    for _ in range(maximum_degree):
        powers.append(amplitude_multiply(powers[-1], matrix))
    return tuple(powers)


def evaluate_amplitude_polynomial(
    coefficients: list[list[Series]],
    powers: tuple[AmplitudeMatrix, ...],
) -> AmplitudeMatrix:
    """Evaluate an amplitude-dependent scalar polynomial."""

    dimension = len(powers[0][0])
    order = powers[0][0][0][0].order
    result = [
        zero_matrix(dimension, dimension, order)
        for _ in range(3)
    ]
    for scalar_degree, power in enumerate(powers):
        for amplitude_degree in range(3):
            for coefficient_degree in range(amplitude_degree + 1):
                result[amplitude_degree] = matrix_add(
                    result[amplitude_degree],
                    matrix_scale(
                        coefficients[coefficient_degree][scalar_degree],
                        power[amplitude_degree - coefficient_degree],
                    ),
                )
    return result[0], result[1], result[2]


def invert_constant_amplitude_matrix(
    matrix: Matrix,
    maximum_c_degree: int,
) -> Matrix:
    """Invert a matrix with identity constant term by a Neumann series."""

    dimension = len(matrix)
    order = maximum_c_degree + 1
    identity = identity_matrix(dimension, order)
    perturbation = matrix_add(matrix, matrix_scale(-1, identity))
    result = identity
    power = identity
    for degree in range(1, maximum_c_degree + 1):
        power = matrix_multiply(power, perturbation)
        result = matrix_add(
            result,
            matrix_scale(-1 if degree % 2 else 1, power),
        )
    product = matrix_multiply(matrix, result)
    if any(
        (
            product[row][column]
            - (1 if row == column else 0)
        ).valuation()
        != order
        for row in range(dimension)
        for column in range(dimension)
    ):
        raise AssertionError("the matrix denominator inverse failed")
    return result


def add_matrices(*matrices: Matrix) -> Matrix:
    """Add a nonempty sequence of matrices."""

    result = matrices[0]
    for matrix in matrices[1:]:
        result = matrix_add(result, matrix)
    return result


def blaschke_amplitude_coefficients(
    numerator: AmplitudeMatrix,
    denominator: AmplitudeMatrix,
    maximum_c_degree: int,
) -> AmplitudeMatrix:
    """Divide the two amplitude-truncated matrix series."""

    inverse_constant = invert_constant_amplitude_matrix(
        denominator[0],
        maximum_c_degree,
    )
    inverse_linear = matrix_scale(
        -1,
        matrix_multiply(
            inverse_constant,
            matrix_multiply(denominator[1], inverse_constant),
        ),
    )
    inverse_quadratic = matrix_scale(
        -1,
        matrix_multiply(
            inverse_constant,
            add_matrices(
                matrix_multiply(
                    denominator[1],
                    inverse_linear,
                ),
                matrix_multiply(
                    denominator[2],
                    inverse_constant,
                ),
            ),
        ),
    )
    constant = matrix_multiply(numerator[0], inverse_constant)
    linear = add_matrices(
        matrix_multiply(numerator[0], inverse_linear),
        matrix_multiply(numerator[1], inverse_constant),
    )
    quadratic = add_matrices(
        matrix_multiply(numerator[0], inverse_quadratic),
        matrix_multiply(numerator[1], inverse_linear),
        matrix_multiply(numerator[2], inverse_constant),
    )
    return constant, linear, quadratic


def weighted_gram(
    left: Matrix,
    metric: Matrix,
    right: Matrix,
) -> Matrix:
    """Return ``left^T metric right``."""

    return matrix_multiply(
        matrix_transpose(left),
        matrix_multiply(metric, right),
    )


def expected_first_face_tangent(
    length: int,
    grade: int,
    order: int,
) -> Matrix:
    """Return L143's sparse noncentral first-face formula."""

    expected = zero_matrix(length + 1, length + 1, order)
    expected[0][length - grade] = Series.constant(-4, order)
    for row in range(grade, length + 1):
        expected[row][row - grade] += Series.constant(4, order)
    expected[grade][length - 2 * grade] += Series.constant(4, order)
    return expected


def formal_dual_calculation(
    length: int,
    grades: tuple[int, ...],
    maximum_c_degree: int,
) -> FormalDualCalculation:
    """Return the exact quadratic dual series for a multi-grade direction."""

    if (
        not grades
        or len(set(grades)) != len(grades)
        or any(not 1 <= grade <= length / 2 for grade in grades)
    ):
        raise ValueError("grades must be distinct low representatives")
    order = maximum_c_degree + 1
    prepared = weierstrass_numerator_coefficients(
        length,
        grades,
        maximum_c_degree,
    )
    numerator_coefficients = matrix_polynomial_coefficients(
        prepared,
        length,
        maximum_c_degree,
        reverse=False,
    )
    denominator_coefficients = matrix_polynomial_coefficients(
        prepared,
        length,
        maximum_c_degree,
        reverse=True,
    )
    direction = [0] * (length - 1)
    for grade in grades:
        direction[grade - 1] = 1
        direction[length - grade - 1] = 1
    operator, _ = operator_expansion_from_coefficients(
        length + 1,
        direction,
        order,
    )
    powers = amplitude_matrix_powers(operator, length)
    numerator = evaluate_amplitude_polynomial(
        numerator_coefficients,
        powers,
    )
    denominator = evaluate_amplitude_polynomial(
        denominator_coefficients,
        powers,
    )
    blaschke = blaschke_amplitude_coefficients(
        numerator,
        denominator,
        maximum_c_degree,
    )

    coordinate_constant = diagonal_matrix(
        [
            Series.constant(
                Fraction(1, 2)
                if index in (0, length)
                else 1,
                order,
            )
            for index in range(length + 1)
        ]
    )
    coordinate_linear = coordinate_metric_tangent(
        length + 1,
        direction,
        order,
    )
    constant_gram = weighted_gram(
        blaschke[0],
        coordinate_constant,
        blaschke[0],
    )
    linear_gram = add_matrices(
        weighted_gram(
            blaschke[1],
            coordinate_constant,
            blaschke[0],
        ),
        weighted_gram(
            blaschke[0],
            coordinate_constant,
            blaschke[1],
        ),
        weighted_gram(
            blaschke[0],
            coordinate_linear,
            blaschke[0],
        ),
    )
    quadratic_gram = add_matrices(
        weighted_gram(
            blaschke[2],
            coordinate_constant,
            blaschke[0],
        ),
        weighted_gram(
            blaschke[0],
            coordinate_constant,
            blaschke[2],
        ),
        weighted_gram(
            blaschke[1],
            coordinate_constant,
            blaschke[1],
        ),
        weighted_gram(
            blaschke[1],
            coordinate_linear,
            blaschke[0],
        ),
        weighted_gram(
            blaschke[0],
            coordinate_linear,
            blaschke[1],
        ),
    )
    if any(
        constant_gram[row][column].valuation() != order
        for row in range(length + 1)
        for column in range(length + 1)
        if row != column
    ):
        raise AssertionError("the elliptic-axis Gramian was not diagonal")

    top = length
    linear_eigenvalue = (
        linear_gram[top][top] - 4 * coordinate_linear[top][top]
    ) / coordinate_constant[top][top]
    first_coupling = [
        zero(order)
        for _ in range(length + 1)
    ]
    for row in range(length):
        right_side = -(
            linear_gram[row][top]
            - 4 * coordinate_linear[row][top]
            - linear_eigenvalue * coordinate_constant[row][top]
        )
        gap = (
            constant_gram[row][row]
            - 4 * coordinate_constant[row][row]
        )
        first_coupling[row] = right_side / gap

    direct_rayleigh = (
        quadratic_gram[top][top]
        - linear_eigenvalue * coordinate_linear[top][top]
    ) / coordinate_constant[top][top]
    schur = zero(order)
    for index in range(length + 1):
        schur += (
            linear_gram[top][index]
            - 4 * coordinate_linear[top][index]
            - linear_eigenvalue * coordinate_constant[top][index]
        ) * first_coupling[index]
    schur /= coordinate_constant[top][top]
    return FormalDualCalculation(
        operator=operator,
        denominator=denominator,
        blaschke=blaschke,
        coordinate_constant=coordinate_constant,
        coordinate_linear=coordinate_linear,
        linear_eigenvalue=linear_eigenvalue,
        first_coupling=tuple(first_coupling),
        direct_rayleigh=direct_rayleigh,
        schur=schur,
        quadratic_eigenvalue=direct_rayleigh + schur,
    )


def make_record(length: int, grade: int) -> FormalFaberBlaschkeRecord:
    """Construct and validate one exact noncentral face."""

    if not 1 <= grade < length / 2:
        raise ValueError("the formal checker requires a noncentral grade")
    target_degree = 2 * grade
    order = target_degree + 1
    calculation = formal_dual_calculation(
        length,
        (grade,),
        target_degree,
    )
    blaschke = calculation.blaschke

    first_face = [
        [
            Series.constant(
                blaschke[1][row][column].coefficient(grade),
                order,
            )
            for column in range(length + 1)
        ]
        for row in range(length + 1)
    ]
    first_face_tangent_matches = (
        first_face
        == expected_first_face_tangent(length, grade, order)
    )
    if not first_face_tangent_matches:
        raise AssertionError("the sparse first-face tangent failed")

    first_coupling = calculation.first_coupling
    direct_rayleigh = calculation.direct_rayleigh
    schur = calculation.schur
    quadratic_eigenvalue = calculation.quadratic_eigenvalue
    linear_eigenvalue_vanishes = (
        calculation.linear_eigenvalue.valuation() == order
    )

    first_nonzero = quadratic_eigenvalue.valuation()
    leading = quadratic_eigenvalue.coefficient(target_degree)
    if (
        not linear_eigenvalue_vanishes
        or first_nonzero != target_degree
        or leading != -64
        or schur.coefficient(target_degree) != 0
    ):
        raise AssertionError("the corrected dual face was not -64")

    endpoint_alias_degree = length - grade
    endpoint_alias_active = endpoint_alias_degree <= target_degree
    if endpoint_alias_active:
        alias = blaschke[1][endpoint_alias_degree][length].coefficient(
            endpoint_alias_degree
        )
        if alias != -2:
            raise AssertionError("the endpoint alias coefficient failed")

    return FormalFaberBlaschkeRecord(
        length=length,
        dimension=length + 1,
        grade=grade,
        target_degree=target_degree,
        first_nonzero_degree=first_nonzero,
        leading_coefficient=str(leading),
        predicted_coefficient="-64",
        lower_coefficients_vanish=all(
            quadratic_eigenvalue.coefficient(degree) == 0
            for degree in range(target_degree)
        ),
        linear_eigenvalue_vanishes=linear_eigenvalue_vanishes,
        first_face_tangent_matches=first_face_tangent_matches,
        endpoint_alias_degree=endpoint_alias_degree,
        endpoint_alias_active=endpoint_alias_active,
        first_singular_coupling_valuation=min(
            coupling.valuation()
            for coupling in first_coupling
        ),
        direct_rayleigh_coefficient=str(
            direct_rayleigh.coefficient(target_degree)
        ),
        schur_coefficient=str(schur.coefficient(target_degree)),
    )


def parse_cases(values: list[str]) -> tuple[tuple[int, int], ...]:
    """Parse repeated ``L,k`` command-line values."""

    if not values:
        return DEFAULT_CASES
    cases = []
    for value in values:
        length_text, grade_text = value.split(",", maxsplit=1)
        cases.append((int(length_text), int(grade_text)))
    return tuple(cases)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="L,K",
        help="noncentral length and grade; may be repeated",
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the selected exact audits."""

    args = parse_args()
    records = [
        make_record(length, grade)
        for length, grade in parse_cases(args.case)
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    print("\n".join(lines), flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
