#!/usr/bin/env python3
"""Regenerate first normal stationarity along Crabb equality branches.

For a phase-one real equality polynomial

    g(z) = z^L + 2 sum_(j=1)^(L-1) u_j z^j,
    u_j = u_(L-j),

let ``A`` be L123's exact companion realization and ``K`` its
coordinate Gramian.  The first ellipse-pullback tangent is

    E = A^dagger_K - A^3.

This checker prepares the first numerator tangent of
``G_c(Psi_c)`` and verifies exactly that the simple top singular
value of the corrected inner evaluation has zero first derivative.
It also separates the two cancellations:

1. varying the finite inner function at fixed ``A``; and
2. varying ``A`` in the tangent ``E`` at fixed inner function.

The symbolic cases use independent equality coefficients through
``L=4``.  Exact rational cases cover two deterministic directions
for every ``2 <= L <= 10``.  These finite checks regenerate L148's
endpoint recurrence; they are not a replacement for its all-size
induction.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class NormalStationarityRecord:
    """One exact first-normal derivative calculation."""

    length: int
    dimension: int
    case_name: str
    symbolic: bool
    prepared_remainder_degree: int
    equality_blaschke_identity: bool
    inner_tangent_derivative: str
    operator_tangent_derivative: str
    total_derivative: str


def polynomial_evaluation(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    matrix: sp.Matrix,
) -> sp.Matrix:
    """Evaluate a nonzero scalar polynomial on a square matrix."""

    if polynomial == 0:
        return sp.zeros(matrix.rows)
    poly = sp.Poly(polynomial, variable)
    powers = [sp.eye(matrix.rows)]
    for _ in range(poly.degree()):
        powers.append(powers[-1] * matrix)
    return sum(
        (
            coefficient * powers[degree]
            for (degree,), coefficient in poly.terms()
        ),
        sp.zeros(matrix.rows),
    )


def polynomial_frechet_derivative(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    matrix: sp.Matrix,
    tangent: sp.Matrix,
) -> sp.Matrix:
    """Evaluate the Fréchet derivative of a scalar polynomial."""

    result = sp.zeros(matrix.rows)
    for (degree,), coefficient in sp.Poly(
        polynomial,
        variable,
    ).terms():
        for left_degree in range(degree):
            result += (
                coefficient
                * matrix**left_degree
                * tangent
                * matrix ** (degree - 1 - left_degree)
            )
    return result


def reversed_polynomial(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    length: int,
) -> sp.Expr:
    """Reverse a real polynomial into a degree-``length`` denominator."""

    return sp.expand(
        sum(
            coefficient * variable ** (length - degree)
            for (degree,), coefficient in sp.Poly(
                polynomial,
                variable,
            ).terms()
        )
    )


def phase_palindromic_coefficients(
    length: int,
    low_coefficients: tuple[sp.Expr, ...],
) -> dict[int, sp.Expr]:
    """Lift low real coefficients to the full phase-one vector."""

    coefficients: dict[int, sp.Expr] = {}
    for grade, coefficient in enumerate(low_coefficients, start=1):
        coefficients[grade] = coefficient
        coefficients[length - grade] = coefficient
    return coefficients


def companion_data(
    length: int,
    coefficients: dict[int, sp.Expr],
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return ``H``, ``K``, the companion ``A``, and its ellipse tangent."""

    dimension = length + 1
    toeplitz = sp.zeros(dimension)
    shift = sp.zeros(dimension)
    for row in range(length):
        toeplitz[row, row] = sp.Rational(1, 2)
        for column in range(row + 1, length):
            value = coefficients[column - row]
            toeplitz[row, column] = value
            toeplitz[column, row] = value
        shift[row, row + 1] = 1

    coordinate_gramian = (
        toeplitz + shift.T * toeplitz * shift
    )
    companion = sp.simplify(
        2 * coordinate_gramian.inv() * toeplitz * shift
    )
    tangent = sp.simplify(
        coordinate_gramian.inv()
        * companion.T
        * coordinate_gramian
        - companion**3
    )
    return toeplitz, coordinate_gramian, companion, tangent


def first_prepared_numerator(
    length: int,
    coefficients: dict[int, sp.Expr],
    variable: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the equality numerator and first prepared remainder."""

    ellipse_parameter = sp.symbols("c", real=True)
    scalar = sp.symbols("z", real=True)
    dickson = [sp.Integer(2), scalar]
    for _ in range(2, length + 1):
        dickson.append(
            sp.expand(
                scalar * dickson[-1]
                - ellipse_parameter * dickson[-2]
            )
        )
    faber = dickson[length] + sum(
        (
            2 * coefficients[degree] * dickson[degree]
            for degree in range(1, length)
        ),
        sp.Integer(0),
    )
    inverse_map = variable + ellipse_parameter * variable**3
    composed = sp.expand(faber.subs(scalar, inverse_map))
    equality_numerator = sp.expand(
        composed.subs(ellipse_parameter, 0)
    )
    first_variation = sp.expand(
        sp.diff(composed, ellipse_parameter).subs(
            ellipse_parameter,
            0,
        )
    )
    _, remainder = sp.div(
        first_variation,
        equality_numerator,
        variable,
    )
    return equality_numerator, sp.expand(remainder)


def blaschke_derivatives(
    numerator: sp.Expr,
    numerator_tangent: sp.Expr,
    variable: sp.Symbol,
    length: int,
    matrix: sp.Matrix,
    matrix_tangent: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Evaluate an inner quotient and its function/operator derivatives."""

    denominator = reversed_polynomial(
        numerator,
        variable,
        length,
    )
    denominator_tangent = reversed_polynomial(
        numerator_tangent,
        variable,
        length,
    )
    numerator_value = polynomial_evaluation(
        numerator,
        variable,
        matrix,
    )
    denominator_value = polynomial_evaluation(
        denominator,
        variable,
        matrix,
    )
    numerator_inner_derivative = polynomial_evaluation(
        numerator_tangent,
        variable,
        matrix,
    )
    denominator_inner_derivative = polynomial_evaluation(
        denominator_tangent,
        variable,
        matrix,
    )
    numerator_operator_derivative = polynomial_frechet_derivative(
        numerator,
        variable,
        matrix,
        matrix_tangent,
    )
    denominator_operator_derivative = polynomial_frechet_derivative(
        denominator,
        variable,
        matrix,
        matrix_tangent,
    )
    denominator_inverse = denominator_value.inv()
    blaschke = sp.simplify(
        numerator_value * denominator_inverse
    )
    inner_derivative = sp.simplify(
        (
            numerator_inner_derivative
            - blaschke * denominator_inner_derivative
        )
        * denominator_inverse
    )
    operator_derivative = sp.simplify(
        (
            numerator_operator_derivative
            - blaschke * denominator_operator_derivative
        )
        * denominator_inverse
    )
    return blaschke, inner_derivative, operator_derivative


def exact_zero(value: sp.Expr) -> sp.Expr:
    """Normalize a rational or low-size symbolic expression."""

    return sp.factor(sp.cancel(value))


def make_record(
    length: int,
    case_name: str,
    low_coefficients: tuple[sp.Expr, ...],
    *,
    symbolic: bool,
) -> NormalStationarityRecord:
    """Run one exact normal-stationarity calculation."""

    variable = sp.symbols("w", real=True)
    coefficients = phase_palindromic_coefficients(
        length,
        low_coefficients,
    )
    numerator, numerator_tangent = first_prepared_numerator(
        length,
        coefficients,
        variable,
    )
    (
        toeplitz,
        coordinate_gramian,
        companion,
        operator_tangent,
    ) = companion_data(length, coefficients)

    (
        blaschke,
        inner_derivative,
        operator_derivative,
    ) = blaschke_derivatives(
        numerator,
        numerator_tangent,
        variable,
        length,
        companion,
        operator_tangent,
    )
    total_derivative = inner_derivative + operator_derivative

    defect = toeplitz[:, 0]
    reversal = sp.zeros(length + 1)
    for index in range(length + 1):
        reversal[index, length - index] = 1
    reverse_defect = reversal * defect
    right_vector = coordinate_gramian.inv() * reverse_defect
    expected_blaschke = (
        4 * sp.eye(length + 1)[:, 0] * reverse_defect.T
    )
    equality_identity = (
        sp.simplify(blaschke - expected_blaschke)
        == sp.zeros(length + 1)
    )
    if not equality_identity:
        raise AssertionError("the equality Blaschke identity failed")

    def endpoint_scalar(derivative: sp.Matrix) -> sp.Expr:
        return exact_zero(
            (defect.T * derivative * right_vector)[0]
        )

    inner_scalar = endpoint_scalar(inner_derivative)
    total_scalar = endpoint_scalar(total_derivative)
    operator_scalar = exact_zero(total_scalar - inner_scalar)
    if any(
        value != 0
        for value in (
            inner_scalar,
            operator_scalar,
            total_scalar,
        )
    ):
        raise AssertionError("the first normal derivative did not vanish")

    return NormalStationarityRecord(
        length=length,
        dimension=length + 1,
        case_name=case_name,
        symbolic=symbolic,
        prepared_remainder_degree=sp.Poly(
            numerator_tangent,
            variable,
        ).degree(),
        equality_blaschke_identity=equality_identity,
        inner_tangent_derivative=str(inner_scalar),
        operator_tangent_derivative=str(operator_scalar),
        total_derivative=str(total_scalar),
    )


def symbolic_records() -> list[NormalStationarityRecord]:
    """Return generic symbolic checks through length four."""

    records = []
    for length in range(2, 5):
        coefficients = sp.symbols(
            f"u1:{length // 2 + 1}",
            real=True,
        )
        records.append(
            make_record(
                length,
                "generic",
                tuple(coefficients),
                symbolic=True,
            )
        )
    return records


def rational_records(
    maximum_length: int,
) -> list[NormalStationarityRecord]:
    """Return two deterministic exact rational checks per length."""

    records = []
    for length in range(2, maximum_length + 1):
        for variant in range(2):
            low_coefficients = tuple(
                sp.Rational(
                    (-1) ** (grade + variant)
                    * (grade + variant + 1),
                    100 + 3 * length,
                )
                for grade in range(1, length // 2 + 1)
            )
            records.append(
                make_record(
                    length,
                    f"rational_{variant}",
                    low_coefficients,
                    symbolic=False,
                )
            )
    return records


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-length", type=int, default=10)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run and optionally persist all exact records."""

    args = parse_args()
    if args.maximum_length < 4:
        raise ValueError("maximum length must be at least four")
    records = [
        *symbolic_records(),
        *rational_records(args.maximum_length),
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
