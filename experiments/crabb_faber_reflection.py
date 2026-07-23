#!/usr/bin/env python3
"""Audit the endpoint reflection law for palindromic Crabb tangents.

Let ``P_m`` be the Dickson polynomials and

    S_0 = C + c J C J.

For the elementary equality coefficient at grade ``j``, the compatible
physical-adjoint tangent is

    H_j = E_j + c J E_(L-j) J,
    E_j = 2 (e_0-e_L) e_(j+1)^*.

The all-size algebra proved in the companion note says that the two
endpoint rows of ``D P_L(S_0)[H_j]`` vanish.  Consequently the endpoint
rows of the Faber transform of the equality polynomial contain its
analytic coefficients at the top and their c-weighted reflection at
the bottom.

The checker uses exact symbolic polynomial arithmetic.  It tests both
real and imaginary components of every phase-one coefficient pair.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class ReflectionRecord:
    length: int
    tested_real_directions: int
    tested_imaginary_directions: int
    maximum_polynomial_degree: int


def reversal(dimension: int) -> sp.Matrix:
    """Return coordinate reversal."""

    result = sp.zeros(dimension)
    for index in range(dimension):
        result[index, dimension - 1 - index] = 1
    return result


def crabb_shift(length: int) -> sp.Matrix:
    """Return the size-``length+1`` Crabb shift."""

    dimension = length + 1
    result = sp.zeros(dimension)
    result[0, 1] = 2
    for column in range(2, dimension):
        result[column - 1, column] = 1
    return result


def dickson_matrices(
    length: int,
    parameter: sp.Symbol,
) -> tuple[sp.Matrix, ...]:
    """Return ``P_0(S_0),...,P_L(S_0)`` exactly."""

    dimension = length + 1
    shift = crabb_shift(length)
    flip = reversal(dimension)
    operator = shift + parameter * flip * shift * flip
    values = [2 * sp.eye(dimension), operator]
    for degree in range(2, length + 1):
        values.append(
            sp.expand(
                operator * values[-1]
                - parameter * values[-2]
            )
        )
    return tuple(values)


def elementary_direction(length: int, grade: int) -> sp.Matrix:
    """Return ``E_j=2(e_0-e_L)e_(j+1)^*``."""

    result = sp.zeros(length + 1)
    result[0, grade + 1] = 2
    result[length, grade + 1] = -2
    return result


def frechet_dickson_values(
    base_values: tuple[sp.Matrix, ...],
    tangent: sp.Matrix,
    parameter: sp.Symbol,
) -> tuple[sp.Matrix, ...]:
    """Return every ``D P_m(S_0)[tangent]`` by recurrence."""

    operator = base_values[1]
    values = [sp.zeros(operator.rows), tangent]
    for degree in range(2, len(base_values)):
        values.append(
            sp.expand(
                tangent * base_values[degree - 1]
                + operator * values[-1]
                - parameter * values[-2]
            )
        )
    return tuple(values)


def assert_zero(expression: sp.Expr, message: str) -> None:
    """Raise if an exact polynomial did not vanish."""

    if sp.expand(expression) != 0:
        raise AssertionError(f"{message}: {sp.factor(expression)}")


def audit_length(length: int) -> ReflectionRecord:
    """Audit every independent real and imaginary coefficient direction."""

    dimension = length + 1
    parameter = sp.symbols("c", real=True)
    flip = reversal(dimension)
    values = dickson_matrices(length, parameter)

    for degree, value in enumerate(values):
        for column in range(dimension):
            assert_zero(
                value[0, column] - (2 if column == degree else 0),
                "top endpoint Dickson row failed",
            )
            assert_zero(
                value[length, column]
                - (
                    2 * parameter**degree
                    if column == length - degree
                    else 0
                ),
                "bottom endpoint Dickson row failed",
            )

    real_count = 0
    imaginary_count = 0
    for first_grade in range(1, length // 2 + 1):
        phases = (1,) if 2 * first_grade == length else (1, -1)
        for phase in phases:
            coefficients = [0] * (length - 1)
            coefficients[first_grade - 1] = 1
            coefficients[length - first_grade - 1] = phase

            raw_tangent = sp.zeros(dimension)
            for grade, coefficient in enumerate(coefficients, start=1):
                raw_tangent += (
                    coefficient
                    * elementary_direction(length, grade)
                )

            # After factoring out i, complex conjugation changes the
            # sign of an imaginary direction.
            tangent = (
                raw_tangent
                + phase * parameter * flip * raw_tangent * flip
            )
            derivatives = frechet_dickson_values(
                values,
                tangent,
                parameter,
            )
            derivative = derivatives[-1]
            for row in (0, length):
                for column in range(dimension):
                    assert_zero(
                        derivative[row, column],
                        "endpoint Frechet row failed",
                    )

            # Equation (4b) in the proof note exposes the complete
            # endpoint cancellation, not only its final zero.
            atomic = elementary_direction(length, first_grade)
            atomic_tangent = (
                atomic
                + parameter
                * flip
                * elementary_direction(
                    length,
                    length - first_grade,
                )
                * flip
            )
            atomic_derivatives = frechet_dickson_values(
                values,
                atomic_tangent,
                parameter,
            )
            for degree in range(1, length + 1):
                expected = sp.zeros(1, dimension)
                if first_grade + degree <= length:
                    expected[0, first_grade + degree] += 2
                else:
                    expected[
                        0,
                        2 * length - first_grade - degree,
                    ] += (
                        2
                        * parameter
                        ** (first_grade + degree - length)
                    )
                if degree <= first_grade:
                    expected[0, first_grade - degree] -= (
                        2 * parameter**degree
                    )
                else:
                    expected[0, degree - first_grade] -= (
                        2 * parameter**first_grade
                    )
                for column in range(dimension):
                    assert_zero(
                        atomic_derivatives[degree][0, column]
                        - expected[0, column],
                        "explicit top endpoint formula failed",
                    )
                    assert_zero(
                        atomic_derivatives[degree][length, column]
                        + expected[0, column],
                        "explicit bottom endpoint formula failed",
                    )

            faber_derivative = derivative.copy()
            for grade, coefficient in enumerate(coefficients, start=1):
                faber_derivative += (
                    2 * coefficient * values[grade]
                )

            for column in range(dimension):
                top_expected = (
                    4 * coefficients[column - 1]
                    if 1 <= column < length
                    else 0
                )
                bottom_expected = sum(
                    4
                    * coefficient
                    * parameter**grade
                    for grade, coefficient in enumerate(
                        coefficients,
                        start=1,
                    )
                    if length - grade == column
                )
                assert_zero(
                    faber_derivative[0, column] - top_expected,
                    "top Faber reflection row failed",
                )
                assert_zero(
                    faber_derivative[length, column]
                    - bottom_expected,
                    "bottom Faber reflection row failed",
                )

            if phase == 1:
                real_count += 1
            else:
                imaginary_count += 1

    maximum_degree = max(
        (
            sp.Poly(entry, parameter).degree()
            for value in values
            for entry in value
            if entry != 0
        ),
        default=0,
    )
    return ReflectionRecord(
        length=length,
        tested_real_directions=real_count,
        tested_imaginary_directions=imaginary_count,
        maximum_polynomial_degree=maximum_degree,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=2)
    parser.add_argument("--maximum-length", type=int, default=14)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact audit and optionally persist JSONL."""

    args = parse_args()
    if (
        args.minimum_length < 2
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError(
            "length range must satisfy 2 <= minimum <= maximum"
        )

    records = [
        audit_length(length)
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
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
