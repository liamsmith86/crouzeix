#!/usr/bin/env python3
"""Exact quadratic response of disk tangents against Crabb normals.

For a homogeneous Toeplitz disk ray ``z=x h``, the frozen
characteristic-Blaschke ambient gradient can start at order ``x^2``.
This checker pairs it with every coercive circular-normal support
covector at the Crabb point and compares the result with the closed
exterior-square formula used in L173.  In particular it preserves the
nonzero length-six response that disproved the earlier all-size
cancellation shortcut.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from formal_riemann_series import laurent_modes


SERIES_ORDER = 2


MatrixSeries = list[sp.Matrix]
ScalarSeries = list[sp.Expr]


@dataclass(frozen=True)
class CircularNormalQuadraticRecord:
    """Exact normal-gradient data for one dimension and support mode."""

    dimension: int
    length: int
    support_mode: int
    polarization: str
    disk_direction: tuple[str, ...]
    gradient_coefficients: tuple[str, ...]
    quadratic_coefficient: str
    predicted_quadratic_coefficient: str
    formula_residual: str
    formula_verified: bool
    quadratic_cancellation: bool


def clean(matrix: sp.Matrix) -> sp.Matrix:
    """Expand every matrix entry without expensive global simplification."""

    return matrix.applyfunc(sp.expand)


def series_add(left: MatrixSeries, right: MatrixSeries) -> MatrixSeries:
    """Add equal-length matrix series."""

    if len(left) != len(right):
        raise ValueError("matrix series must have equal lengths")
    return [
        clean(left_coefficient + right_coefficient)
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def series_subtract(left: MatrixSeries, right: MatrixSeries) -> MatrixSeries:
    """Subtract equal-length matrix series."""

    if len(left) != len(right):
        raise ValueError("matrix series must have equal lengths")
    return [
        clean(left_coefficient - right_coefficient)
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def series_multiply(left: MatrixSeries, right: MatrixSeries) -> MatrixSeries:
    """Multiply equal-length matrix series through their common order."""

    if len(left) != len(right):
        raise ValueError("matrix series must have equal lengths")
    order = len(left) - 1
    dimension = left[0].rows
    result = [sp.zeros(dimension) for _ in range(order + 1)]
    for degree in range(order + 1):
        for left_degree in range(degree + 1):
            result[degree] += (
                left[left_degree] * right[degree - left_degree]
            )
        result[degree] = clean(result[degree])
    return result


def scalar_matrix_series(
    coefficients: ScalarSeries,
    dimension: int,
) -> MatrixSeries:
    """Embed a scalar series as multiples of the identity."""

    identity = sp.eye(dimension)
    return [
        sp.expand(coefficient) * identity
        for coefficient in coefficients
    ]


def series_inverse(coefficients: MatrixSeries) -> MatrixSeries:
    """Invert a matrix series whose constant term is nonsingular."""

    inverse = [coefficients[0].inv()]
    for degree in range(1, len(coefficients)):
        convolution = sum(
            (
                coefficients[source_degree]
                * inverse[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.zeros(coefficients[0].rows),
        )
        inverse.append(clean(-inverse[0] * convolution))
    return inverse


def series_power(coefficients: MatrixSeries, exponent: int) -> MatrixSeries:
    """Raise a matrix series to a nonnegative integer power."""

    dimension = coefficients[0].rows
    result = [
        sp.eye(dimension),
        *[sp.zeros(dimension) for _ in coefficients[1:]],
    ]
    base = coefficients
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = series_multiply(result, base)
        base = series_multiply(base, base)
        remaining //= 2
    return result


def polynomial_frechet_series(
    coefficients: Sequence[ScalarSeries],
    operator: MatrixSeries,
    direction: MatrixSeries,
) -> tuple[MatrixSeries, MatrixSeries]:
    """Evaluate a varying polynomial and its matrix Frechet derivative."""

    dimension = operator[0].rows
    zero = [sp.zeros(dimension) for _ in operator]
    value = zero
    derivative = [matrix.copy() for matrix in zero]
    for coefficient in coefficients:
        previous_value = value
        derivative = series_add(
            series_multiply(derivative, operator),
            series_multiply(previous_value, direction),
        )
        value = series_add(
            series_multiply(previous_value, operator),
            scalar_matrix_series(coefficient, dimension),
        )
    return value, derivative


def deterministic_disk_direction(length: int) -> tuple[sp.Expr, ...]:
    """Return a generic Gaussian-rational Toeplitz direction."""

    values = [sp.Integer(0)]
    for offset in range(1, length):
        values.append(
            sp.Rational(2 * offset + 1, 20 * length)
            + sp.I * sp.Rational(offset + 1, 25 * length)
        )
    return tuple(values)


def closed_quadratic_response(
    direction: Sequence[sp.Expr],
    mode: int,
) -> sp.Expr:
    """Return L173's complex quadratic support-normal response.

    Coordinates ``1,...,L-1`` are the nonconstant Toeplitz
    coefficients and ``J`` reverses those coordinates.  The summand is
    one coefficient of ``h wedge J conjugate(h)``.
    """

    length = len(direction)
    response = sp.Integer(0)
    for left in range(1, length):
        right = length + mode - left
        if not left < right < length:
            continue
        wedge = (
            direction[left] * sp.conjugate(direction[length - right])
            - direction[right]
            * sp.conjugate(direction[length - left])
        )
        response += (right - left) * wedge
    return sp.expand(8 * (4 * mode - 1) * response / length**2)


def disk_model_series(
    direction: Sequence[sp.Expr],
    order: int = SERIES_ORDER,
    correction: sp.Matrix | None = None,
) -> tuple[MatrixSeries, MatrixSeries]:
    """Construct a coefficient-gauge disk operator and metric series.

    ``correction`` is the quadratic Hermitian coefficient in
    ``H=I/2+sZ+s^2 correction``.  The default reproduces L173's linear
    Toeplitz ray.
    """

    if order < 0:
        raise ValueError("series order must be nonnegative")
    length = len(direction)
    dimension = length + 1

    base_toeplitz = sp.zeros(dimension)
    tangent_toeplitz = sp.zeros(dimension)
    for index in range(length):
        base_toeplitz[index, index] = sp.Rational(1, 2)
    for offset in range(1, length):
        for row in range(length - offset):
            tangent_toeplitz[row, row + offset] = direction[offset]
            tangent_toeplitz[row + offset, row] = sp.conjugate(
                direction[offset]
            )

    hermitian = [
        base_toeplitz,
        *[sp.zeros(dimension) for _ in range(order)],
    ]
    if order >= 1:
        hermitian[1] = tangent_toeplitz
    if correction is not None:
        if correction.shape != (length, length):
            raise ValueError("correction has the wrong shape")
        if order < 2:
            raise ValueError("a quadratic correction requires order >= 2")
        hermitian[2][:length, :length] = correction

    return disk_model_from_hermitian_series(hermitian)


def disk_model_from_hermitian_series(
    hermitian: MatrixSeries,
) -> tuple[MatrixSeries, MatrixSeries]:
    """Construct the disk operator from an arbitrary Hermitian series."""

    if not hermitian:
        raise ValueError("the Hermitian series cannot be empty")
    dimension = hermitian[0].rows
    order = len(hermitian) - 1
    if any(
        coefficient.shape != (dimension, dimension)
        for coefficient in hermitian
    ):
        raise ValueError("Hermitian series coefficients have wrong shapes")
    nilpotent_shift = sp.zeros(dimension)
    for row in range(dimension - 1):
        nilpotent_shift[row, row + 1] = 1

    metric = [
        coefficient
        + nilpotent_shift.T * coefficient * nilpotent_shift
        for coefficient in hermitian
    ]
    metric_inverse = series_inverse(metric)
    operator = series_multiply(
        series_multiply(metric_inverse, hermitian),
        [
            2 * nilpotent_shift,
            *[sp.zeros(dimension) for _ in range(order)],
        ],
    )
    return operator, metric


def circle_star(expression: sp.Expr, boundary: sp.Symbol) -> sp.Expr:
    """Conjugate a Laurent polynomial on the unit circle."""

    return sp.expand(
        sp.conjugate(expression).subs(
            sp.conjugate(boundary),
            1 / boundary,
        )
    )


def support_series(
    metric: MatrixSeries,
    direction: sp.Matrix,
    boundary: sp.Symbol,
) -> ScalarSeries:
    """Return the first-support variation through the supplied order."""

    dimension = direction.rows
    monomials = sp.Matrix(
        [boundary**index for index in range(dimension)]
    )
    adjoint_monomials = sp.Matrix(
        [[boundary ** (-index) for index in range(dimension)]]
    )
    denominator = [
        sp.expand(
            (adjoint_monomials * metric[degree] * monomials)[0]
        )
        for degree in range(len(metric))
    ]
    analytic_numerator = [
        sp.expand(
            boundary**-1
            * (
                adjoint_monomials
                * metric[degree]
                * direction
                * monomials
            )[0]
        )
        for degree in range(len(metric))
    ]
    numerator = [
        sp.expand(
            (
                analytic_numerator[degree]
                + circle_star(analytic_numerator[degree], boundary)
            )
            / 2
        )
        for degree in range(len(metric))
    ]

    constant = sp.expand(denominator[0])
    if boundary in constant.free_symbols:
        raise AssertionError("the Crabb support denominator is not constant")
    denominator_inverse = [1 / constant]
    for degree in range(1, len(denominator)):
        convolution = sum(
            (
                denominator[source_degree]
                * denominator_inverse[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.Integer(0),
        )
        denominator_inverse.append(
            sp.expand(-convolution / constant)
        )
    return [
        sp.expand(
            sum(
                numerator[source_degree]
                * denominator_inverse[degree - source_degree]
                for source_degree in range(degree + 1)
            )
        )
        for degree in range(len(denominator))
    ]


def polynomial_at_operator_series(
    scalar_coefficients: ScalarSeries,
    operator: MatrixSeries,
    boundary: sp.Symbol,
) -> MatrixSeries:
    """Evaluate the inverse-Riemann correction on an operator series."""

    if len(scalar_coefficients) != len(operator):
        raise ValueError("scalar and operator series must have equal lengths")
    order = len(operator) - 1
    dimension = operator[0].rows
    result = [sp.zeros(dimension) for _ in range(order + 1)]
    for scalar_degree, expression in enumerate(scalar_coefficients):
        for mode, coefficient in laurent_modes(
            expression,
            boundary,
        ).items():
            if mode < 0:
                raise AssertionError(
                    "the inverse-Riemann correction is not analytic"
                )
            power = series_power(operator, mode)
            for matrix_degree in range(
                order - scalar_degree + 1
            ):
                result[scalar_degree + matrix_degree] += (
                    coefficient * power[matrix_degree]
                )
    return [clean(matrix) for matrix in result]


def pulled_direction_series(
    operator: MatrixSeries,
    metric: MatrixSeries,
    direction: sp.Matrix,
) -> MatrixSeries:
    """Subtract the exact first inverse-Riemann correction."""

    boundary = sp.symbols("boundary")
    support = support_series(metric, direction, boundary)
    inverse_map = [sp.Integer(0) for _ in operator]
    for degree in range(len(operator)):
        modes = laurent_modes(support[degree], boundary)
        inverse_map[degree] = sp.expand(
            modes.get(0, 0) * boundary
            + 2
            * sum(
                (
                    coefficient * boundary ** (mode + 1)
                    for mode, coefficient in modes.items()
                    if mode >= 1
                ),
                sp.Integer(0),
            )
        )
    correction = polynomial_at_operator_series(
        inverse_map,
        operator,
        boundary,
    )
    raw = [
        direction,
        *[sp.zeros(direction.rows) for _ in operator[1:]],
    ]
    return series_subtract(raw, correction)


def characteristic_series(operator: MatrixSeries) -> list[ScalarSeries]:
    """Return descending characteristic coefficients through series order."""

    order = len(operator) - 1
    parameter = sp.symbols("parameter", real=True)
    matrix = sum(
        (
            parameter**degree * operator[degree]
            for degree in range(order + 1)
        ),
        sp.zeros(operator[0].rows),
    )
    coefficients = matrix.charpoly().all_coeffs()
    return [
        [
            sp.expand(coefficient).coeff(parameter, degree)
            for degree in range(order + 1)
        ]
        for coefficient in coefficients
    ]


def directional_gradient_series(
    operator: MatrixSeries,
    metric: MatrixSeries,
    characteristic: Sequence[ScalarSeries],
    direction: sp.Matrix,
) -> ScalarSeries:
    """Return the frozen characteristic-Blaschke norm derivative."""

    numerator_coefficients = characteristic[:-1]
    denominator_coefficients = [
        [sp.conjugate(value) for value in coefficient]
        for coefficient in numerator_coefficients[::-1]
    ]
    pulled = pulled_direction_series(operator, metric, direction)
    numerator, numerator_derivative = polynomial_frechet_series(
        numerator_coefficients,
        operator,
        pulled,
    )
    denominator, denominator_derivative = polynomial_frechet_series(
        denominator_coefficients,
        operator,
        pulled,
    )
    denominator_inverse = series_inverse(denominator)
    blaschke = series_multiply(numerator, denominator_inverse)
    blaschke_derivative = series_subtract(
        series_multiply(numerator_derivative, denominator_inverse),
        series_multiply(
            series_multiply(blaschke, denominator_derivative),
            denominator_inverse,
        ),
    )

    endpoint = sp.eye(operator[0].rows)[:, -1]
    defect = [
        metric[degree] * sp.eye(operator[0].rows)[:, 0]
        for degree in range(len(metric))
    ]
    result = []
    for degree in range(len(operator)):
        value = sum(
            (
                (
                    sp.conjugate(defect[left_degree]).T
                    * blaschke_derivative[degree - left_degree]
                    * endpoint
                )[0]
                for left_degree in range(degree + 1)
            ),
            sp.Integer(0),
        )
        result.append(sp.simplify(8 * sp.re(value)))
    return result


def crabb_normal_directions(
    metric: MatrixSeries,
    mode: int,
) -> tuple[sp.Matrix, sp.Matrix]:
    """Return real and imaginary Riesz representatives of one support mode."""

    dimension = metric[0].rows
    boundary = sp.symbols("boundary")
    real_values = [sp.Integer(0)] * (2 * dimension**2)
    imaginary_values = [sp.Integer(0)] * (2 * dimension**2)
    for imaginary_basis in (False, True):
        for entry in range(dimension**2):
            direction = sp.zeros(dimension)
            direction[entry] = sp.I if imaginary_basis else 1
            constant_support = support_series(
                metric,
                direction,
                boundary,
            )[0]
            coefficient = laurent_modes(
                constant_support,
                boundary,
            ).get(mode, 0)
            index = entry + (
                dimension**2 if imaginary_basis else 0
            )
            real_values[index] = sp.re(coefficient)
            imaginary_values[index] = sp.im(coefficient)

    def to_matrix(values: Sequence[sp.Expr]) -> sp.Matrix:
        matrix = sp.zeros(dimension)
        for entry in range(dimension**2):
            matrix[entry] = (
                values[entry]
                + sp.I * values[dimension**2 + entry]
            )
        return matrix

    return to_matrix(real_values), to_matrix(imaginary_values)


def records_for_length(length: int) -> list[CircularNormalQuadraticRecord]:
    """Regenerate every coercive normal polarization in one dimension."""

    direction = deterministic_disk_direction(length)
    operator, metric = disk_model_series(direction)
    characteristic = characteristic_series(operator)
    records = []
    for mode in range(3, length + 2):
        predicted_complex = closed_quadratic_response(direction, mode)
        normal_directions = crabb_normal_directions(metric, mode)
        for polarization, normal in zip(
            ("real", "imaginary"),
            normal_directions,
            strict=True,
        ):
            gradient = directional_gradient_series(
                operator,
                metric,
                characteristic,
                normal,
            )
            quadratic = sp.simplify(gradient[2])
            predicted = sp.simplify(
                sp.re(predicted_complex)
                if polarization == "real"
                else sp.im(predicted_complex)
            )
            residual = sp.simplify(quadratic - predicted)
            if residual != 0:
                raise AssertionError(
                    "closed quadratic response formula failed for "
                    f"length={length}, mode={mode}, "
                    f"polarization={polarization}: {residual}"
                )
            records.append(
                CircularNormalQuadraticRecord(
                    dimension=length + 1,
                    length=length,
                    support_mode=mode,
                    polarization=polarization,
                    disk_direction=tuple(
                        str(value) for value in direction[1:]
                    ),
                    gradient_coefficients=tuple(
                        str(sp.factor(value)) for value in gradient
                    ),
                    quadratic_coefficient=str(quadratic),
                    predicted_quadratic_coefficient=str(predicted),
                    formula_residual=str(residual),
                    formula_verified=True,
                    quadratic_cancellation=bool(quadratic == 0),
                )
            )
    return records


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=7)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the exact covariant quadratic audit."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("length range must satisfy 3 <= minimum <= maximum")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        ):
            for record in records_for_length(length):
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")


if __name__ == "__main__":
    main()
