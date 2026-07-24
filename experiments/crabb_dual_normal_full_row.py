#!/usr/bin/env python3
"""Regenerate complete frozen-dual circular-normal rows exactly.

The checker retains two nilpotent first-order parameters: the equality
amplitude and the strong circular normal.  It constructs the physical
inverse-Riemann pullback, the frozen prepared Blaschke rational
function, and the moving simple top singular value through the full
target weight.

This is an end-to-end adversarial audit of L172.  It is deliberately
independent of the residue proof and retains the nonzero grade-one
case as a discriminator.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import asdict
import argparse
import json
from pathlib import Path
from typing import Generic, TypeVar

import sympy as sp

from crabb_circular_normal_series import (
    physical_reflected_path,
    real_circular_normal_direction,
)
from crabb_faber_blaschke_formal import (
    weierstrass_numerator_coefficients,
)
from general_crabb_weighted_series import inverse_riemann_series


T = TypeVar("T")


DEFAULT_CASES = (
    (3, 1),
    (4, 1),
    (5, 1),
    (4, 2),
    (5, 2),
    (6, 3),
)


@dataclass(frozen=True)
class Jet(Generic[T]):
    """Coefficient of ``1, a, s, a*s`` with ``a^2=s^2=0``."""

    constant: T
    amplitude: T
    strong: T
    mixed: T


@dataclass(frozen=True)
class FullDualNormalRowRecord:
    """One complete exact mixed frozen-dual row."""

    dimension: int
    length: int
    reflected_grade: int
    target_weight: int
    mixed_top_singular_square: str
    predicted_mixed_top_singular_square: str
    prediction_verified: bool
    grade_one_discriminator: bool


def zero_matrix(rows: int, columns: int | None = None) -> sp.Matrix:
    """Return a symbolic zero matrix."""

    return sp.zeros(rows, rows if columns is None else columns)


def matrix_jet_constant(value: sp.Matrix) -> Jet[sp.Matrix]:
    """Embed a constant matrix in the bilinear jet algebra."""

    zero = zero_matrix(value.rows, value.cols)
    return Jet(value, zero, zero, zero)


def scalar_jet_constant(value: sp.Expr) -> Jet[sp.Expr]:
    """Embed a scalar in the bilinear jet algebra."""

    return Jet(value, sp.Integer(0), sp.Integer(0), sp.Integer(0))


def jet_add(left: Jet[T], right: Jet[T]) -> Jet[T]:
    """Add two jets."""

    return Jet(
        left.constant + right.constant,
        left.amplitude + right.amplitude,
        left.strong + right.strong,
        left.mixed + right.mixed,
    )


def jet_scale(scalar: sp.Expr, value: Jet[T]) -> Jet[T]:
    """Scale a jet by a constant scalar."""

    return Jet(
        scalar * value.constant,
        scalar * value.amplitude,
        scalar * value.strong,
        scalar * value.mixed,
    )


def jet_multiply(left: Jet, right: Jet) -> Jet:
    """Multiply jets, retaining only the bilinear quotient."""

    return Jet(
        left.constant * right.constant,
        (
            left.amplitude * right.constant
            + left.constant * right.amplitude
        ),
        left.strong * right.constant + left.constant * right.strong,
        (
            left.mixed * right.constant
            + left.amplitude * right.strong
            + left.strong * right.amplitude
            + left.constant * right.mixed
        ),
    )


def jet_transpose(value: Jet[sp.Matrix]) -> Jet[sp.Matrix]:
    """Transpose a matrix jet."""

    return Jet(
        value.constant.T,
        value.amplitude.T,
        value.strong.T,
        value.mixed.T,
    )


def zero_matrix_jet(rows: int, columns: int | None = None) -> Jet[sp.Matrix]:
    """Return a zero matrix jet."""

    return matrix_jet_constant(zero_matrix(rows, columns))


def series_add(left: list[Jet], right: list[Jet]) -> list[Jet]:
    """Add equally truncated jet series."""

    return [jet_add(a, b) for a, b in zip(left, right, strict=True)]


def series_multiply(left: list[Jet], right: list[Jet]) -> list[Jet]:
    """Multiply equally truncated jet series."""

    order = min(len(left), len(right)) - 1
    result = []
    for degree in range(order + 1):
        value = None
        for source in range(degree + 1):
            term = jet_multiply(left[source], right[degree - source])
            value = term if value is None else jet_add(value, term)
        if value is None:
            raise AssertionError("empty series convolution")
        result.append(value)
    return result


def matrix_series_identity(dimension: int, order: int) -> list[Jet[sp.Matrix]]:
    """Return the constant identity matrix series."""

    return [
        matrix_jet_constant(sp.eye(dimension)),
        *[
            zero_matrix_jet(dimension)
            for _ in range(order)
        ],
    ]


def matrix_series_inverse(
    value: list[Jet[sp.Matrix]],
) -> list[Jet[sp.Matrix]]:
    """Invert a matrix jet series whose constant term is ordinary."""

    order = len(value) - 1
    if any(
        component != sp.zeros(value[0].constant.rows)
        for component in (
            value[0].amplitude,
            value[0].strong,
            value[0].mixed,
        )
    ):
        raise ValueError(
            "the inverse recurrence requires an ordinary degree-zero term"
        )
    inverse_constant = value[0].constant.inv()
    result = [
        matrix_jet_constant(inverse_constant),
        *[zero_matrix_jet(inverse_constant.rows) for _ in range(order)],
    ]
    for degree in range(1, order + 1):
        convolution = zero_matrix_jet(inverse_constant.rows)
        for source in range(1, degree + 1):
            convolution = jet_add(
                convolution,
                jet_multiply(value[source], result[degree - source]),
            )
        result[degree] = jet_scale(
            -1,
            jet_multiply(result[0], convolution),
        )
    return result


def scalar_matrix_series_multiply(
    scalar: list[Jet[sp.Expr]],
    matrix: list[Jet[sp.Matrix]],
) -> list[Jet[sp.Matrix]]:
    """Multiply a scalar jet series into a matrix jet series."""

    return series_multiply(scalar, matrix)


def evaluate_polynomial(
    coefficients: list[list[Jet[sp.Expr]]],
    matrix: list[Jet[sp.Matrix]],
) -> list[Jet[sp.Matrix]]:
    """Evaluate a scalar-series polynomial on a matrix-series jet."""

    order = len(matrix) - 1
    dimension = matrix[0].constant.rows
    powers = [matrix_series_identity(dimension, order)]
    for _ in range(1, len(coefficients)):
        powers.append(series_multiply(powers[-1], matrix))
    result = [zero_matrix_jet(dimension) for _ in range(order + 1)]
    for coefficient, power in zip(coefficients, powers, strict=True):
        result = series_add(
            result,
            scalar_matrix_series_multiply(coefficient, power),
        )
    return result


def prepared_scalar_coefficients(
    length: int,
    grade: int,
    order: int,
) -> tuple[list[list[Jet[sp.Expr]]], list[list[Jet[sp.Expr]]]]:
    """Return numerator and reversed-denominator scalar jet series."""

    prepared = weierstrass_numerator_coefficients(
        length,
        (grade,),
        order,
    )
    numerator = []
    for scalar_degree in range(length + 1):
        coefficient_series = []
        for physical_degree in range(order + 1):
            constant_polynomial = prepared.get(
                (0, physical_degree),
                [],
            )
            constant = (
                sp.Rational(constant_polynomial[scalar_degree])
                if scalar_degree < len(constant_polynomial)
                else sp.Integer(0)
            )
            # The formal preparation treats the equality amplitude as
            # independent of ``c``.  The physical path inserts it as
            # ``epsilon*amplitude``, hence this one-degree shift.
            amplitude_polynomial = prepared.get(
                (1, physical_degree - 1),
                [],
            )
            amplitude = (
                sp.Rational(amplitude_polynomial[scalar_degree])
                if scalar_degree < len(amplitude_polynomial)
                else sp.Integer(0)
            )
            coefficient_series.append(
                Jet(
                    constant,
                    amplitude,
                    sp.Integer(0),
                    sp.Integer(0),
                )
            )
        numerator.append(coefficient_series)
    denominator = list(reversed(numerator))
    return numerator, denominator


def operator_jet(
    length: int,
    grade: int,
    order: int,
) -> list[Jet[sp.Matrix]]:
    """Return the physical normalized operator as a bilinear jet."""

    dimension = length + 1
    face_degree = grade + 1
    amplitude, strong = sp.symbols("amplitude strong", real=True)
    normal = real_circular_normal_direction(
        dimension,
        length + 2 - grade,
    )
    path = physical_reflected_path(
        dimension,
        grade,
        strong,
        normal,
        face_degree,
        order,
        amplitude=amplitude,
        ellipse=sp.Integer(1),
    )
    _, operator = inverse_riemann_series(path, order)
    result = []
    for coefficient in operator:
        constant = coefficient.subs({amplitude: 0, strong: 0})
        amplitude_part = sp.diff(coefficient, amplitude).subs(
            {amplitude: 0, strong: 0}
        )
        strong_part = sp.diff(coefficient, strong).subs(
            {amplitude: 0, strong: 0}
        )
        mixed = sp.diff(coefficient, amplitude, strong).subs(
            {amplitude: 0, strong: 0}
        )
        result.append(
            Jet(
                constant,
                amplitude_part,
                strong_part,
                mixed,
            )
        )
    return result


def simple_top_eigenpair(
    hermitian: list[Jet[sp.Matrix]],
) -> tuple[list[Jet[sp.Expr]], list[Jet[sp.Matrix]]]:
    """Lift the simple top eigenvalue with endpoint normalization."""

    order = len(hermitian) - 1
    dimension = hermitian[0].constant.rows
    endpoint = dimension - 1
    if any(
        component != sp.zeros(dimension)
        for component in (
            hermitian[0].amplitude,
            hermitian[0].strong,
            hermitian[0].mixed,
        )
    ):
        raise ValueError(
            "the eigenpair recurrence requires an ordinary degree-zero term"
        )
    if any(
        hermitian[0].constant[row, column] != 0
        for row in range(dimension)
        for column in range(dimension)
        if row != column
    ):
        raise ValueError("the degree-zero Gram matrix must be diagonal")
    base_value = hermitian[0].constant[endpoint, endpoint]
    base_vector = sp.eye(dimension)[:, endpoint]
    eigenvalues = [scalar_jet_constant(base_value)]
    eigenvectors = [matrix_jet_constant(base_vector)]
    reduced_inverse = sp.zeros(dimension)
    for index in range(dimension - 1):
        if hermitian[0].constant[index, index] == base_value:
            raise ValueError("the top eigenvalue must be simple")
        reduced_inverse[index, index] = 1 / (
            hermitian[0].constant[index, index] - base_value
        )

    for degree in range(1, order + 1):
        raw = zero_matrix_jet(dimension, 1)
        for source in range(1, degree + 1):
            raw = jet_add(
                raw,
                jet_multiply(
                    hermitian[source],
                    eigenvectors[degree - source],
                ),
            )
        for source in range(1, degree):
            raw = jet_add(
                raw,
                jet_scale(
                    -1,
                    jet_multiply(
                        eigenvalues[source],
                        eigenvectors[degree - source],
                    ),
                ),
            )
        eigenvalue = Jet(
            raw.constant[endpoint],
            raw.amplitude[endpoint],
            raw.strong[endpoint],
            raw.mixed[endpoint],
        )
        eigenvalues.append(eigenvalue)
        residual = jet_add(
            raw,
            jet_scale(
                -1,
                jet_multiply(eigenvalue, eigenvectors[0]),
            ),
        )
        eigenvectors.append(
            jet_scale(
                -1,
                jet_multiply(
                    matrix_jet_constant(reduced_inverse),
                    residual,
                ),
            )
        )
    return eigenvalues, eigenvectors


def simple_top_eigenvalue(
    hermitian: list[Jet[sp.Matrix]],
) -> list[Jet[sp.Expr]]:
    """Return only the eigenvalue part of :func:`simple_top_eigenpair`."""

    return simple_top_eigenpair(hermitian)[0]


def mixed_dual_row(length: int, grade: int) -> sp.Expr:
    """Return the exact target mixed coefficient."""

    order = 2 * (grade + 1)
    operator = operator_jet(length, grade, order)
    numerator_coefficients, denominator_coefficients = (
        prepared_scalar_coefficients(length, grade, order)
    )
    numerator = evaluate_polynomial(numerator_coefficients, operator)
    denominator = evaluate_polynomial(denominator_coefficients, operator)
    blaschke = series_multiply(
        numerator,
        matrix_series_inverse(denominator),
    )
    gram = series_multiply(
        [jet_transpose(value) for value in blaschke],
        blaschke,
    )
    eigenvalue = simple_top_eigenvalue(gram)
    return sp.factor(eigenvalue[order].mixed)


def make_record(length: int, grade: int) -> FullDualNormalRowRecord:
    """Return one full exact row and check its predicted value."""

    if length < 3 or not 1 <= grade <= length // 2:
        raise ValueError("require L>=3 and 1<=k<=floor(L/2)")
    actual = mixed_dual_row(length, grade)
    predicted = (
        -4 * sp.Rational(5 * length - 1, length)
        if grade == 1
        else sp.Integer(0)
    )
    verified = sp.simplify(actual - predicted) == 0
    if not verified:
        raise AssertionError("the complete frozen-dual row changed")
    return FullDualNormalRowRecord(
        dimension=length + 1,
        length=length,
        reflected_grade=grade,
        target_weight=2 * (grade + 1),
        mixed_top_singular_square=str(actual),
        predicted_mixed_top_singular_square=str(predicted),
        prediction_verified=verified,
        grade_one_discriminator=(grade == 1 and actual != 0),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic exact cases."""

    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for length, grade in DEFAULT_CASES:
            record = make_record(length, grade)
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
