#!/usr/bin/env python3
"""Exact jet audit for the candidate one-reflection ``O(Q)`` bound.

This checker uses a truncated-series ring over exact SymPy rationals.
It verifies the two finite algebraic inputs required by L155:

1. at a phase-palindromic equality anchor, the one-reflection
   coefficient and its anti-palindromic normal derivative vanish;
2. at the Crabb apex, the same coefficient has no disk-amplitude
   term below degree four.

It also audits the proposed endpoint route to the first condition:
at one generic complex equality anchor in each length, the cleared
resolvent numerator vanishes and its normalized Gram-block gradient
is supported only on the two forbidden diagonal-normalization
directions.  Polynomial interpolation in the resolvent variable and
a full matrix-gradient calculation replace random tangent sampling.

Finally, it checks the four inverse-free endpoint recurrences used in
the all-size proof of that gradient identity.  Those recurrences
separate an arbitrary complex disk tangent into its direct and
reverse-conjugate coefficient polynomials.

As a downstream audit of L157, it constructs the normalized
orbit-complement defect and verifies that its value and every real
disk-coordinate derivative agree with the canonical disk defect at a
phase-one equality anchor.

The calculation implements the root-free tangent (2)--(7) from
``proof/crabb_disk_one_reflection.md``.  Finite verification is not an
all-size proof of the companion recurrence.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import TypeAlias

import sympy as sp


ScalarSeries: TypeAlias = list[sp.Expr]
MatrixSeries: TypeAlias = list[sp.Matrix]


@dataclass(frozen=True)
class NormalJetRecord:
    """One exact equality-normal first jet."""

    length: int
    dimension: int
    equality_coefficients: tuple[str, ...]
    normal_direction: tuple[str, ...]
    value: str
    normal_derivative: str


@dataclass(frozen=True)
class ApexJetRecord:
    """One exact Crabb-apex amplitude jet."""

    length: int
    dimension: int
    direction: tuple[int, ...]
    coefficients_through_four: tuple[str, ...]
    vanishes_below_degree_four: bool


@dataclass(frozen=True)
class ResolventSquareRecord:
    """Full real first jet of L149's cleared endpoint residual."""

    length: int
    dimension: int
    equality_coefficients: tuple[str, ...]
    interpolation_points: tuple[int, ...]
    off_diagonal_tangent_dimension: int
    numerator_polynomial_zero: bool
    normalized_gradient_endpoint_support: bool


@dataclass(frozen=True)
class EndpointRecurrenceRecord:
    """Inverse-free audit of the all-size endpoint recurrence."""

    length: int
    dimension: int
    equality_coefficients: tuple[str, ...]
    tangent_coefficients: tuple[str, ...]
    tangent_decomposition: bool
    right_endpoint_recurrence: bool
    left_endpoint_recurrence: bool
    bilinear_cancellation: bool


@dataclass(frozen=True)
class ModelDefectJetRecord:
    """Full coordinate first jet of the normalized model defect."""

    length: int
    dimension: int
    equality_coefficients: tuple[str, ...]
    real_tangent_dimension: int
    normalized_defect_matches_canonical_value: bool
    normalized_defect_matches_canonical_gradient: bool


class SeriesRing:
    """Small exact truncated-series algebra."""

    def __init__(self, order: int) -> None:
        self.order = order

    @staticmethod
    def zero_like(value: sp.Expr | sp.Matrix) -> sp.Expr | sp.Matrix:
        """Return an additive zero with the same scalar/matrix type."""

        if isinstance(value, sp.MatrixBase):
            return sp.zeros(*value.shape)
        return sp.Integer(0)

    def constant(
        self,
        value: sp.Expr | sp.Matrix,
    ) -> ScalarSeries | MatrixSeries:
        """Embed a scalar or matrix as a constant series."""

        return [
            value,
            *[
                self.zero_like(value)
                for _ in range(self.order)
            ],
        ]

    @staticmethod
    def add(
        left: ScalarSeries | MatrixSeries,
        right: ScalarSeries | MatrixSeries,
    ) -> ScalarSeries | MatrixSeries:
        """Add two equally truncated series."""

        return [
            left[index] + right[index]
            for index in range(len(left))
        ]

    @staticmethod
    def subtract(
        left: ScalarSeries | MatrixSeries,
        right: ScalarSeries | MatrixSeries,
    ) -> ScalarSeries | MatrixSeries:
        """Subtract two equally truncated series."""

        return [
            left[index] - right[index]
            for index in range(len(left))
        ]

    def multiply(
        self,
        left: ScalarSeries | MatrixSeries,
        right: ScalarSeries | MatrixSeries,
    ) -> ScalarSeries | MatrixSeries:
        """Multiply series by truncated convolution."""

        zero = self.zero_like(left[0] * right[0])
        return [
            sum(
                (
                    left[left_degree]
                    * right[degree - left_degree]
                    for left_degree in range(degree + 1)
                ),
                zero,
            )
            for degree in range(self.order + 1)
        ]

    @staticmethod
    def scale(
        series: ScalarSeries | MatrixSeries,
        scalar: sp.Expr,
    ) -> ScalarSeries | MatrixSeries:
        """Multiply every coefficient by a scalar."""

        return [scalar * coefficient for coefficient in series]

    @staticmethod
    def adjoint(series: MatrixSeries) -> MatrixSeries:
        """Take the Hermitian adjoint coefficientwise."""

        return [coefficient.conjugate().T for coefficient in series]

    @staticmethod
    def conjugate(series: ScalarSeries) -> ScalarSeries:
        """Conjugate a scalar series in a real perturbation variable."""

        return [sp.conjugate(coefficient) for coefficient in series]

    @staticmethod
    def trace(series: MatrixSeries) -> ScalarSeries:
        """Take the trace coefficientwise."""

        return [sp.trace(coefficient) for coefficient in series]

    def matrix_inverse(self, series: MatrixSeries) -> MatrixSeries:
        """Invert a matrix series recursively."""

        inverse = [
            series[0].inv(),
            *[
                sp.zeros(series[0].rows)
                for _ in range(self.order)
            ],
        ]
        for degree in range(1, self.order + 1):
            convolution = sum(
                (
                    series[source_degree]
                    * inverse[degree - source_degree]
                    for source_degree in range(1, degree + 1)
                ),
                sp.zeros(series[0].rows),
            )
            inverse[degree] = -inverse[0] * convolution
        return inverse

    def scalar_inverse(self, series: ScalarSeries) -> ScalarSeries:
        """Invert a scalar series recursively."""

        inverse = [
            1 / series[0],
            *[sp.Integer(0) for _ in range(self.order)],
        ]
        for degree in range(1, self.order + 1):
            inverse[degree] = -inverse[0] * sum(
                series[source_degree]
                * inverse[degree - source_degree]
                for source_degree in range(1, degree + 1)
            )
        return inverse

    def matrix_power(
        self,
        series: MatrixSeries,
        exponent: int,
    ) -> MatrixSeries:
        """Raise a matrix series to a nonnegative power."""

        result = self.constant(sp.eye(series[0].rows))
        for _ in range(exponent):
            result = self.multiply(result, series)
        return result

    def evaluate_polynomial(
        self,
        coefficients: list[ScalarSeries],
        matrix: MatrixSeries,
    ) -> MatrixSeries:
        """Evaluate a series-coefficient polynomial."""

        result = self.constant(sp.zeros(matrix[0].rows))
        power = self.constant(sp.eye(matrix[0].rows))
        for coefficient in coefficients:
            result = self.add(
                result,
                self.multiply(coefficient, power),
            )
            power = self.multiply(power, matrix)
        return result

    def frechet_polynomial(
        self,
        coefficients: list[ScalarSeries],
        matrix: MatrixSeries,
        tangent: MatrixSeries,
    ) -> MatrixSeries:
        """Evaluate a polynomial Fréchet derivative over the ring."""

        result = self.constant(sp.zeros(matrix[0].rows))
        powers = [
            self.matrix_power(matrix, exponent)
            for exponent in range(len(coefficients))
        ]
        for exponent, coefficient in enumerate(coefficients):
            for left_power in range(exponent):
                term = self.multiply(
                    powers[left_power],
                    self.multiply(
                        tangent,
                        powers[exponent - 1 - left_power],
                    ),
                )
                result = self.add(
                    result,
                    self.multiply(coefficient, term),
                )
        return result

    def characteristic_coefficients(
        self,
        matrix: MatrixSeries,
    ) -> list[ScalarSeries]:
        """Return descending characteristic coefficients by Faddeev."""

        dimension = matrix[0].rows
        identity = self.constant(sp.eye(dimension))
        work = identity
        coefficients = [self.constant(sp.Integer(1))]
        for degree in range(1, dimension + 1):
            product = self.multiply(matrix, work)
            coefficient = self.scale(
                self.trace(product),
                -sp.Rational(1, degree),
            )
            coefficients.append(coefficient)
            work = self.add(
                product,
                [
                    scalar * sp.eye(dimension)
                    for scalar in coefficient
                ],
            )
        return coefficients

    def polynomial_remainder(
        self,
        dividend: list[ScalarSeries],
        monic_divisor: list[ScalarSeries],
    ) -> list[ScalarSeries]:
        """Divide by a monic polynomial over the series ring."""

        remainder = [coefficient[:] for coefficient in dividend]
        zero = self.constant(sp.Integer(0))
        while len(remainder) >= len(monic_divisor):
            leading = remainder[-1]
            shift = len(remainder) - len(monic_divisor)
            for index, coefficient in enumerate(monic_divisor):
                remainder[shift + index] = self.subtract(
                    remainder[shift + index],
                    self.multiply(leading, coefficient),
                )
            remainder.pop()
        while len(remainder) < len(monic_divisor):
            remainder.append(zero)
        return remainder


def forward_shift(dimension: int) -> sp.Matrix:
    """Return the superdiagonal forward shift."""

    shift = sp.zeros(dimension)
    for row in range(dimension - 1):
        shift[row, row + 1] = 1
    return shift


def hermitian_toeplitz(
    coefficients: tuple[sp.Expr, ...],
    *,
    diagonal: sp.Expr = sp.Rational(1, 2),
) -> sp.Matrix:
    """Extend one Hermitian Toeplitz block by a final zero coordinate."""

    length = len(coefficients) + 1
    dimension = length + 1
    toeplitz = sp.zeros(dimension)
    for index in range(length):
        toeplitz[index, index] = diagonal
    for offset, value in enumerate(coefficients, start=1):
        for row in range(length - offset):
            toeplitz[row, row + offset] = value
            toeplitz[row + offset, row] = sp.conjugate(value)
    return toeplitz


def first_reflection_series(
    base: tuple[sp.Expr, ...],
    direction: tuple[sp.Expr, ...],
    order: int,
) -> ScalarSeries:
    """Return the disk-amplitude series of the first reflection."""

    ring = SeriesRing(order)
    length = len(base) + 1
    dimension = length + 1
    toeplitz_base = hermitian_toeplitz(base)
    toeplitz_tangent = hermitian_toeplitz(
        direction,
        diagonal=sp.Integer(0),
    )
    shift = forward_shift(dimension)

    toeplitz = [
        toeplitz_base,
        toeplitz_tangent,
        *[sp.zeros(dimension) for _ in range(order - 1)],
    ]
    shift_series = ring.constant(shift)
    coordinate = ring.add(
        toeplitz,
        ring.multiply(
            ring.adjoint(shift_series),
            ring.multiply(toeplitz, shift_series),
        ),
    )
    operator = ring.scale(
        ring.multiply(
            ring.matrix_inverse(coordinate),
            ring.multiply(toeplitz, shift_series),
        ),
        sp.Integer(2),
    )

    characteristic = ring.characteristic_coefficients(operator)
    factor = [
        characteristic[dimension - (power + 1)]
        for power in range(length + 1)
    ]
    tangent_polynomial = [
        ring.constant(sp.Integer(0))
        for _ in range(length + 3)
    ]
    for power in range(1, length + 1):
        tangent_polynomial[power + 2] = ring.add(
            tangent_polynomial[power + 2],
            ring.scale(factor[power], sp.Integer(power)),
        )
    for power in range(2, length + 1):
        tangent_polynomial[power - 2] = ring.subtract(
            tangent_polynomial[power - 2],
            ring.scale(factor[power], sp.Integer(power)),
        )
    tangent_polynomial[length - 1] = ring.add(
        tangent_polynomial[length - 1],
        factor[length - 1],
    )
    prepared_tangent = ring.polynomial_remainder(
        tangent_polynomial,
        factor,
    )

    denominator_factor = [
        ring.conjugate(coefficient)
        for coefficient in reversed(factor)
    ]
    denominator_tangent = [
        ring.conjugate(coefficient)
        for coefficient in reversed(prepared_tangent)
    ]
    denominator = ring.evaluate_polynomial(
        denominator_factor,
        operator,
    )
    blaschke = ring.multiply(
        ring.evaluate_polynomial(factor, operator),
        ring.matrix_inverse(denominator),
    )
    operator_tangent = ring.subtract(
        ring.multiply(
            ring.matrix_inverse(coordinate),
            ring.multiply(
                ring.adjoint(operator),
                coordinate,
            ),
        ),
        ring.matrix_power(operator, 3),
    )
    numerator_derivative = ring.add(
        ring.frechet_polynomial(
            factor,
            operator,
            operator_tangent,
        ),
        ring.evaluate_polynomial(prepared_tangent, operator),
    )
    denominator_derivative = ring.add(
        ring.frechet_polynomial(
            denominator_factor,
            operator,
            operator_tangent,
        ),
        ring.evaluate_polynomial(
            denominator_tangent,
            operator,
        ),
    )
    inverse_denominator = ring.matrix_inverse(denominator)
    blaschke_derivative = ring.subtract(
        ring.multiply(numerator_derivative, inverse_denominator),
        ring.multiply(
            blaschke,
            ring.multiply(
                denominator_derivative,
                inverse_denominator,
            ),
        ),
    )

    endpoint = sp.eye(dimension)[:, -1]
    endpoint_series = ring.constant(endpoint)
    image = ring.multiply(blaschke, endpoint_series)
    image_derivative = ring.multiply(
        blaschke_derivative,
        endpoint_series,
    )
    numerator = ring.multiply(
        ring.adjoint(image),
        ring.multiply(coordinate, image_derivative),
    )
    denominator = ring.multiply(
        ring.adjoint(endpoint_series),
        ring.multiply(coordinate, endpoint_series),
    )
    numerator_scalar = [entry[0] for entry in numerator]
    denominator_scalar = [entry[0] for entry in denominator]
    return ring.scale(
        ring.multiply(
            numerator_scalar,
            ring.scalar_inverse(denominator_scalar),
        ),
        sp.Integer(2),
    )


def complex_equality_data(length: int) -> tuple[sp.Expr, ...]:
    """Return a small generic complex equality point in the phase-one gauge."""

    size = length - 1
    equality: list[sp.Expr | None] = [None] * size
    for index in range((size + 1) // 2):
        partner = size - 1 - index
        real_part = sp.Rational(index + 1, 211 + 7 * index)
        if index == partner:
            value = real_part
        else:
            value = real_part + sp.I * sp.Rational(
                index + 2,
                307 + 11 * index,
            )
        equality[index] = value
        equality[partner] = sp.conjugate(value)
    return tuple(value for value in equality if value is not None)


def canonical_gaussian_rational(value: sp.Expr) -> sp.Expr:
    """Canonicalize an exact scalar known to lie in the Gaussian rationals."""

    return sp.QQ_I.to_sympy(sp.QQ_I.from_sympy(sp.expand(value)))


def gaussian_polynomial_is_zero(
    value: sp.Expr,
    variable: sp.Symbol,
) -> bool:
    """Test a univariate Gaussian-rational polynomial exactly."""

    return all(
        canonical_gaussian_rational(coefficient) == 0
        for coefficient in sp.Poly(sp.expand(value), variable).all_coeffs()
    )


def gaussian_matrix_is_zero(matrix: sp.Matrix) -> bool:
    """Test a Gaussian-rational matrix exactly after expansion."""

    return all(
        canonical_gaussian_rational(entry) == 0
        for entry in matrix
    )


def cleared_resolvent_gradient_at(
    toeplitz: sp.Matrix,
    variable_value: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Matrix]:
    """Evaluate the cleared numerator and normalized Gram-block gradient."""

    dimension = toeplitz.rows
    shift = forward_shift(dimension)
    shift_adjoint = shift.T
    coordinate = toeplitz + shift_adjoint * toeplitz * shift
    pencil = variable_value * coordinate - 2 * toeplitz * shift
    sharp_pencil = coordinate - 2 * variable_value * shift_adjoint * toeplitz
    pencil_inverse = pencil.inv()
    sharp_pencil_inverse = sharp_pencil.inv()
    pencil_determinant = pencil.det()
    sharp_pencil_determinant = sharp_pencil.det()

    left_endpoint = sp.eye(dimension)[0, :]
    right_endpoint = sp.eye(dimension)[:, -1]
    left = left_endpoint * coordinate
    right = coordinate * right_endpoint
    right_solution = pencil_inverse * right
    left_solution = left * pencil_inverse
    transfer = (left_solution * right)[0]
    numerator = canonical_gaussian_rational(
        pencil_determinant * transfer
        - sharp_pencil_determinant,
    )
    coordinate_gradient = (
        right_solution * left_endpoint
        + right_endpoint * left_solution
        - variable_value * right_solution * left_solution
        + transfer * variable_value * pencil_inverse
        - transfer * sharp_pencil_inverse
    )
    gradient = (
        coordinate_gradient
        + shift * coordinate_gradient * shift_adjoint
        + shift
        * (
            2 * right_solution * left_solution
            - 2 * transfer * pencil_inverse
        )
        + 2
        * transfer
        * variable_value
        * sharp_pencil_inverse
        * shift_adjoint
    )
    return numerator, transfer, gradient


def equality_data(
    length: int,
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    """Return deterministic phase-one equality and normal vectors."""

    size = length - 1
    equality: list[sp.Expr | None] = [None] * size
    normal: list[sp.Expr | None] = [None] * size
    for index in range((size + 1) // 2):
        partner = size - 1 - index
        value = sp.Rational(index + 1, 31 + 5 * index)
        equality[index] = value
        equality[partner] = value
        if index == partner:
            normal[index] = sp.Integer(0)
        else:
            tangent = sp.Rational(index + 2, 23 + 3 * index)
            normal[index] = tangent
            normal[partner] = -tangent
    return (
        tuple(value for value in equality if value is not None),
        tuple(value for value in normal if value is not None),
    )


def normal_record(length: int) -> NormalJetRecord:
    """Return one exact equality-normal record."""

    equality, normal = equality_data(length)
    coefficients = first_reflection_series(equality, normal, 1)
    if coefficients != [0, 0]:
        raise AssertionError("the equality-normal jet did not vanish")
    return NormalJetRecord(
        length=length,
        dimension=length + 1,
        equality_coefficients=tuple(map(str, equality)),
        normal_direction=tuple(map(str, normal)),
        value="0",
        normal_derivative="0",
    )


def apex_record(length: int) -> ApexJetRecord:
    """Return one exact Crabb-apex amplitude record."""

    direction = tuple(
        ((3 * index + length) % 7) - 3
        for index in range(length - 1)
    )
    if not any(direction):
        raise AssertionError("the deterministic direction vanished")
    coefficients = first_reflection_series(
        tuple(sp.Integer(0) for _ in direction),
        tuple(map(sp.Integer, direction)),
        4,
    )
    if any(coefficients[:4]):
        raise AssertionError("a subquartic apex jet survived")
    return ApexJetRecord(
        length=length,
        dimension=length + 1,
        direction=direction,
        coefficients_through_four=tuple(map(str, coefficients)),
        vanishes_below_degree_four=True,
    )


def resolvent_record(length: int) -> ResolventSquareRecord:
    """Certify the full first jet at one complex equality anchor."""

    equality = complex_equality_data(length)
    toeplitz = hermitian_toeplitz(equality)
    interpolation_points = tuple(range(2, length + 3))
    for variable_value in interpolation_points:
        value, transfer, gradient = cleared_resolvent_gradient_at(
            toeplitz,
            sp.Integer(variable_value),
        )
        if value != 0:
            raise AssertionError("the endpoint numerator did not vanish")
        expected = sp.zeros(length)
        expected[0, 0] = transfer
        expected[-1, -1] = transfer
        gradient_block = gradient[:length, :length]
        difference = gradient_block - expected
        if any(
            canonical_gaussian_rational(entry) != 0
            for entry in difference
        ):
            raise AssertionError(
                "the endpoint gradient lost boundary-diagonal support",
            )
    return ResolventSquareRecord(
        length=length,
        dimension=length + 1,
        equality_coefficients=tuple(map(str, equality)),
        interpolation_points=interpolation_points,
        off_diagonal_tangent_dimension=length * (length - 1),
        numerator_polynomial_zero=True,
        normalized_gradient_endpoint_support=True,
    )


def endpoint_recurrence_record(length: int) -> EndpointRecurrenceRecord:
    """Audit the direct/reverse endpoint identities without inverses."""

    dimension = length + 1
    equality = complex_equality_data(length)
    tangent = tuple(
        sp.Rational(offset + 2, 101 + 3 * offset)
        + sp.I * sp.Rational(2 * offset + 1, 137 + 5 * offset)
        for offset in range(1, length)
    )
    tangent_sharp = tuple(
        sp.conjugate(tangent[length - 1 - offset])
        for offset in range(1, length)
    )
    toeplitz = hermitian_toeplitz(equality)
    toeplitz_tangent = hermitian_toeplitz(
        tangent,
        diagonal=sp.Integer(0),
    )
    shift = forward_shift(dimension)
    coordinate = toeplitz + shift.T * toeplitz * shift
    coordinate_tangent = (
        toeplitz_tangent
        + shift.T * toeplitz_tangent * shift
    )

    companion = shift.copy()
    companion[0, 1] += 1
    for offset, coefficient in enumerate(equality, start=1):
        companion[0, offset + 1] += 2 * coefficient
        companion[-1, offset + 1] -= 2 * coefficient
    if not gaussian_matrix_is_zero(
        coordinate * companion - 2 * toeplitz * shift
    ):
        raise AssertionError("the equality companion identity failed")

    direct_row = sp.zeros(1, dimension)
    equality_row = sp.zeros(1, dimension)
    sharp_row = sp.zeros(1, dimension)
    normal_vector = sp.zeros(dimension, 1)
    for offset in range(1, length):
        direct_row[0, offset + 1] = tangent[offset - 1]
        equality_row[0, offset + 1] = equality[offset - 1]
        sharp_row[0, offset + 1] = tangent_sharp[offset - 1]
        normal_vector[offset, 0] = (
            tangent[length - offset - 1]
            - sp.conjugate(tangent[offset - 1])
        )
    left_endpoint = sp.eye(dimension)[:, 0]
    right_endpoint = sp.eye(dimension)[:, -1]
    actual_tangent = (
        2 * toeplitz_tangent * shift
        - coordinate_tangent * companion
    )
    expected_tangent = (
        left_endpoint * direct_row
        + 2 * normal_vector * equality_row
        - right_endpoint * sharp_row
    )
    tangent_decomposition = gaussian_matrix_is_zero(
        actual_tangent - expected_tangent
    )
    if not tangent_decomposition:
        raise AssertionError("the tangent decomposition failed")

    variable = sp.symbols("xi")
    equality_polynomial = sum(
        coefficient * variable**offset
        for offset, coefficient in enumerate(equality, start=1)
    )
    direct_polynomial = sum(
        coefficient * variable**offset
        for offset, coefficient in enumerate(tangent, start=1)
    )
    sharp_polynomial = sum(
        coefficient * variable**offset
        for offset, coefficient in enumerate(tangent_sharp, start=1)
    )
    factor = variable**length + 2 * equality_polynomial
    reversed_factor = 1 + 2 * equality_polynomial
    pencil = variable * sp.eye(dimension) - companion
    right_vector = sp.Matrix(
        [
            reversed_factor + 1,
            *[variable**power for power in range(1, length + 1)],
        ]
    )
    left_vector = sp.Matrix(
        [[
            *[
                variable ** (length - index)
                for index in range(length)
            ],
            reversed_factor + 1,
        ]]
    )
    right_residual = (
        pencil * right_vector
        - variable
        * factor
        * sp.eye(dimension)[:, -1]
    )
    left_residual = (
        left_vector * coordinate * pencil
        - variable
        * factor
        * sp.eye(dimension)[0, :]
        * coordinate
    )
    right_endpoint_recurrence = all(
        gaussian_polynomial_is_zero(entry, variable)
        for entry in right_residual
    )
    left_endpoint_recurrence = all(
        gaussian_polynomial_is_zero(entry, variable)
        for entry in left_residual
    )
    if not right_endpoint_recurrence:
        raise AssertionError("the right endpoint recurrence failed")
    if not left_endpoint_recurrence:
        raise AssertionError("the left endpoint recurrence failed")

    bilinear_residual = (
        (left_vector * actual_tangent * right_vector)[0]
        - variable
        * (
            direct_polynomial * factor
            - 2 * reversed_factor * sharp_polynomial
        )
    )
    bilinear_cancellation = gaussian_polynomial_is_zero(
        bilinear_residual,
        variable,
    )
    if not bilinear_cancellation:
        raise AssertionError("the endpoint bilinear cancellation failed")

    return EndpointRecurrenceRecord(
        length=length,
        dimension=dimension,
        equality_coefficients=tuple(map(str, equality)),
        tangent_coefficients=tuple(map(str, tangent)),
        tangent_decomposition=tangent_decomposition,
        right_endpoint_recurrence=right_endpoint_recurrence,
        left_endpoint_recurrence=left_endpoint_recurrence,
        bilinear_cancellation=bilinear_cancellation,
    )


def model_defect_jet_holds(
    equality: tuple[sp.Expr, ...],
    direction: tuple[sp.Expr, ...],
) -> tuple[bool, bool]:
    """Check one normalized orbit-defect value and directional jet."""

    length = len(equality) + 1
    dimension = length + 1
    toeplitz = hermitian_toeplitz(equality)
    toeplitz_tangent = hermitian_toeplitz(
        direction,
        diagonal=sp.Integer(0),
    )
    shift = forward_shift(dimension)
    coordinate = toeplitz + shift.T * toeplitz * shift
    coordinate_tangent = (
        toeplitz_tangent
        + shift.T * toeplitz_tangent * shift
    )
    operator = 2 * coordinate.inv() * toeplitz * shift
    operator_tangent = coordinate.inv() * (
        2 * toeplitz_tangent * shift
        - coordinate_tangent * operator
    )

    ring = SeriesRing(1)
    operator_series = [operator, operator_tangent]
    characteristic = ring.characteristic_coefficients(operator_series)
    factor = [
        characteristic[dimension - (power + 1)]
        for power in range(length + 1)
    ]
    denominator_factor = [
        ring.conjugate(coefficient)
        for coefficient in reversed(factor)
    ]
    denominator = ring.evaluate_polynomial(
        denominator_factor,
        operator_series,
    )
    denominator_inverse = ring.matrix_inverse(denominator)
    endpoint = ring.constant(sp.eye(dimension)[:, -1])

    orbit_constant: list[sp.Matrix] = []
    orbit_tangent: list[sp.Matrix] = []
    for power in range(length):
        operator_power = ring.matrix_power(operator_series, power)
        orbit_column = ring.multiply(
            denominator_inverse,
            ring.multiply(operator_power, endpoint),
        )
        orbit_constant.append(orbit_column[0])
        orbit_tangent.append(orbit_column[1])
    orbit = sp.Matrix.hstack(*orbit_constant)
    orbit_derivative = sp.Matrix.hstack(*orbit_tangent)

    normalization_row = sp.eye(dimension)[0, :]
    system = sp.Matrix.vstack(
        orbit.conjugate().T,
        normalization_row,
    )
    system_derivative = sp.Matrix.vstack(
        orbit_derivative.conjugate().T,
        sp.zeros(1, dimension),
    )
    target = sp.Matrix(
        [*[sp.Integer(0) for _ in range(length)], sp.Rational(1, 2)]
    )
    normalized_defect = system.inv() * target
    normalized_defect_tangent = -system.inv() * (
        system_derivative * normalized_defect
    )
    canonical_defect = toeplitz[:, 0]
    canonical_defect_tangent = toeplitz_tangent[:, 0]
    return (
        gaussian_matrix_is_zero(
            normalized_defect - canonical_defect
        ),
        gaussian_matrix_is_zero(
            normalized_defect_tangent
            - canonical_defect_tangent
        ),
    )


def model_defect_record(length: int) -> ModelDefectJetRecord:
    """Audit every real disk-coordinate direction at one equality point."""

    equality, _ = equality_data(length)
    value_matches = True
    gradient_matches = True
    for offset in range(length - 1):
        for phase in (sp.Integer(1), sp.I):
            direction = tuple(
                phase if index == offset else sp.Integer(0)
                for index in range(length - 1)
            )
            value_holds, gradient_holds = model_defect_jet_holds(
                equality,
                direction,
            )
            value_matches &= value_holds
            gradient_matches &= gradient_holds
    if not value_matches:
        raise AssertionError(
            "the normalized model defect missed its equality value"
        )
    if not gradient_matches:
        raise AssertionError(
            "the normalized model defect lost a coordinate first jet"
        )
    return ModelDefectJetRecord(
        length=length,
        dimension=length + 1,
        equality_coefficients=tuple(map(str, equality)),
        real_tangent_dimension=2 * (length - 1),
        normalized_defect_matches_canonical_value=True,
        normalized_defect_matches_canonical_gradient=True,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-length", type=int, default=10)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run and optionally persist the exact jet grid."""

    args = parse_args()
    records: list[
        NormalJetRecord
        | ApexJetRecord
        | ResolventSquareRecord
        | EndpointRecurrenceRecord
        | ModelDefectJetRecord
    ] = []
    for length in range(3, args.maximum_length + 1):
        records.append(normal_record(length))
        records.append(apex_record(length))
        records.append(resolvent_record(length))
        records.append(endpoint_recurrence_record(length))
        records.append(model_defect_record(length))
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    print("\n".join(lines))
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
