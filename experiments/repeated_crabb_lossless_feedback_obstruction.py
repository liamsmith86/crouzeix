#!/usr/bin/env python3
"""Audit the exact rank obstruction to a one-port lossless delay lift."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp

from crabb_palindromic_elliptic_hessian import direct_map_coefficients


Matrix = sp.Matrix
Series = list[Matrix]


@dataclass(frozen=True)
class LosslessFeedbackObstructionRecord:
    """The exact grade-two monomial obstruction."""

    full_first_active_degree: int
    tail_first_active_degree: int
    full_earlier_coefficients_vanish: bool
    tail_earlier_coefficients_vanish: bool
    full_face: list[list[str]]
    tail_face: list[list[str]]
    full_face_rank: int
    tail_face_rank: int
    full_face_trace: str
    tail_face_trace: str
    all_checks_passed: bool


def zero_series(
    row_dimension: int,
    maximum_degree: int,
    column_dimension: int | None = None,
) -> Series:
    """Return a zero matrix series."""

    columns = row_dimension if column_dimension is None else column_dimension
    return [
        sp.zeros(row_dimension, columns)
        for _ in range(maximum_degree + 1)
    ]


def identity_series(dimension: int, maximum_degree: int) -> Series:
    """Return the constant identity series."""

    result = zero_series(dimension, maximum_degree)
    result[0] = sp.eye(dimension)
    return result


def add_series(left: Series, right: Series) -> Series:
    """Add two equal-length series."""

    return [
        left_coefficient + right_coefficient
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def subtract_series(left: Series, right: Series) -> Series:
    """Subtract two equal-length series."""

    return [
        left_coefficient - right_coefficient
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def multiply_series(
    left: Series,
    right: Series,
    maximum_degree: int,
) -> Series:
    """Multiply two matrix series through one degree."""

    result = zero_series(
        left[0].rows,
        maximum_degree,
        right[0].cols,
    )
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            degree = left_degree + right_degree
            if degree <= maximum_degree:
                result[degree] += left_coefficient * right_coefficient
    return result


def adjoint_series(series: Series) -> Series:
    """Transpose the real exact matrix coefficients."""

    return [coefficient.T for coefficient in series]


def inverse_series(series: Series) -> Series:
    """Invert a series with an invertible constant coefficient."""

    maximum_degree = len(series) - 1
    dimension = series[0].rows
    result = zero_series(dimension, maximum_degree)
    result[0] = series[0].inv()
    for degree in range(1, maximum_degree + 1):
        convolution = sp.zeros(dimension)
        for positive_degree in range(1, degree + 1):
            convolution += (
                series[positive_degree]
                * result[degree - positive_degree]
            )
        result[degree] = -result[0] * convolution
    return result


def power_series(
    series: Series,
    exponent: int,
    maximum_degree: int,
) -> Series:
    """Raise a matrix series to a nonnegative integer power."""

    result = identity_series(series[0].rows, maximum_degree)
    for _ in range(exponent):
        result = multiply_series(result, series, maximum_degree)
    return result


def ellipse_operator_series(
    partial: Matrix,
    initial: Matrix,
    final: Matrix,
    maximum_degree: int,
) -> Series:
    """Return the exact physical ellipse-map series."""

    dimension = partial.rows
    identity = sp.eye(dimension)
    pencil = zero_series(dimension, maximum_degree)
    pencil[0] = partial
    pencil[1] = (
        (identity + final)
        * partial.T
        * (identity + initial)
    )
    result = zero_series(dimension, maximum_degree)
    for index, scalar_series in enumerate(
        direct_map_coefficients(maximum_degree, maximum_degree + 1),
    ):
        power = power_series(
            pencil,
            2 * index + 1,
            maximum_degree,
        )
        for scalar_degree in range(maximum_degree + 1):
            coefficient = scalar_series.coefficient(scalar_degree)
            scalar = sp.Rational(
                coefficient.numerator,
                coefficient.denominator,
            )
            for word_degree in range(
                maximum_degree - scalar_degree + 1,
            ):
                result[scalar_degree + word_degree] += (
                    scalar * power[word_degree]
                )
    return result


def boundary_metric_series(
    partial: Matrix,
    initial: Matrix,
    final: Matrix,
    maximum_degree: int,
) -> Series:
    """Return L219's exact boundary metric."""

    dimension = partial.rows
    result = identity_series(dimension, maximum_degree)
    for degree in range(2, maximum_degree + 1, 2):
        order = degree // 2
        coefficient = (
            partial**order
            * final
            * (partial.T**order)
        )
        for divisor in range(1, order + 1):
            if order % divisor:
                continue
            coefficient += (
                (-1) ** (order // divisor)
                * (partial.T**divisor)
                * initial
                * (partial**divisor)
            )
        result[degree] = coefficient
    return result


def compress_series(
    series: Series,
    left: Matrix,
    right: Matrix,
) -> Series:
    """Apply constant left and right compressions."""

    return [
        left * coefficient * right
        for coefficient in series
    ]


def closed_return_defect_series(
    partial: Matrix,
    initial: Matrix,
    final: Matrix,
    active_degree: int,
) -> Series:
    """Return L258's edge-deleted closed-return defect."""

    dimension = partial.rows
    identity = sp.eye(dimension)
    retained = identity - final
    operator = ellipse_operator_series(
        partial,
        initial,
        final,
        active_degree,
    )
    metric = boundary_metric_series(
        partial,
        initial,
        final,
        active_degree,
    )
    metric[active_degree] = sp.zeros(dimension)

    output_gram = multiply_series(
        multiply_series(
            operator,
            inverse_series(metric),
            active_degree,
        ),
        adjoint_series(operator),
        active_degree,
    )
    retained_gram = compress_series(
        output_gram,
        retained,
        retained,
    )
    entrance = compress_series(
        output_gram,
        retained,
        final,
    )
    exit_ = compress_series(
        output_gram,
        final,
        retained,
    )
    loop = compress_series(
        output_gram,
        final,
        final,
    )
    loop_denominator = subtract_series(
        identity_series(dimension, active_degree),
        loop,
    )
    renewal = add_series(
        retained_gram,
        multiply_series(
            multiply_series(
                entrance,
                inverse_series(loop_denominator),
                active_degree,
            ),
            exit_,
            active_degree,
        ),
    )
    retained_metric = compress_series(
        metric,
        retained,
        retained,
    )
    return subtract_series(
        [
            retained,
            *(
                sp.zeros(dimension)
                for _ in range(active_degree)
            ),
        ],
        multiply_series(
            retained_metric,
            renewal,
            active_degree,
        ),
    )


def string_matrix(matrix: Matrix) -> list[list[str]]:
    """Serialize an exact matrix without losing rational values."""

    return [
        [str(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def exact_record() -> LosslessFeedbackObstructionRecord:
    """Compute and verify the exact monomial rank obstruction."""

    full_partial = sp.Matrix([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
    ])
    full_identity = sp.eye(3)
    full_initial = (
        full_identity
        - full_partial.T * full_partial
    )
    full_final = (
        full_identity
        - full_partial * full_partial.T
    )
    full_series = closed_return_defect_series(
        full_partial,
        full_initial,
        full_final,
        4,
    )
    full_face = full_series[4][1:, 1:]

    tail_partial = sp.Matrix([
        [0, 0],
        [1, 0],
    ])
    tail_identity = sp.eye(2)
    tail_initial = (
        tail_identity
        - tail_partial.T * tail_partial
    )
    tail_final = (
        tail_identity
        - tail_partial * tail_partial.T
    )
    tail_series = closed_return_defect_series(
        tail_partial,
        tail_initial,
        tail_final,
        2,
    )
    tail_face = tail_series[2]

    expected_full = 2 * sp.eye(2)
    expected_tail = sp.diag(0, 4)
    full_earlier_vanish = all(
        coefficient == sp.zeros(3)
        for coefficient in full_series[:4]
    )
    tail_earlier_vanish = all(
        coefficient == sp.zeros(2)
        for coefficient in tail_series[:2]
    )
    verified = bool(
        full_earlier_vanish
        and tail_earlier_vanish
        and full_face == expected_full
        and tail_face == expected_tail
        and full_face.trace() == tail_face.trace() == 4
        and full_face.rank() == 2
        and tail_face.rank() == 1
    )
    if not verified:
        raise RuntimeError("the exact lossless-feedback obstruction failed")
    return LosslessFeedbackObstructionRecord(
        full_first_active_degree=4,
        tail_first_active_degree=2,
        full_earlier_coefficients_vanish=full_earlier_vanish,
        tail_earlier_coefficients_vanish=tail_earlier_vanish,
        full_face=string_matrix(full_face),
        tail_face=string_matrix(tail_face),
        full_face_rank=full_face.rank(),
        tail_face_rank=tail_face.rank(),
        full_face_trace=str(full_face.trace()),
        tail_face_trace=str(tail_face.trace()),
        all_checks_passed=verified,
    )


def write_record(
    record: LosslessFeedbackObstructionRecord,
    output: Path,
) -> str:
    """Write one deterministic record atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(asdict(record), sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_lossless_feedback_obstruction_s70225.jsonl",
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact obstruction."""

    args = parse_args()
    record = exact_record()
    digest = write_record(record, args.output)
    print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
