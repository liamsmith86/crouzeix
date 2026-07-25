#!/usr/bin/env python3
"""Extract the exact length-five elliptic-axis defect-frame jet.

The computation combines the proved periodized-sech metric from the
all-size Crabb axis with the exact theta/ODE direct Riemann recurrence.
It demonstrates that the exact lower-tight frame differs from the
zero-reflection frame already at order four, although the first
transfer coefficient is ``B_5``.
"""

from __future__ import annotations

import argparse
import json
from math import ceil
from pathlib import Path

import sympy as sp

from crabb_palindromic_elliptic_hessian import direct_map_coefficients


LENGTH = 5
DEGREE = 2 * LENGTH


def truncate(
    expression: sp.Expr,
    variable: sp.Symbol,
    degree: int,
) -> sp.Expr:
    """Return a polynomial Taylor truncation."""

    return sp.series(
        expression,
        variable,
        0,
        degree + 1,
    ).removeO().expand()


def truncate_matrix(
    matrix: sp.Matrix,
    variable: sp.Symbol,
    degree: int,
) -> sp.Matrix:
    """Taylor-truncate every matrix entry."""

    return matrix.applyfunc(
        lambda value: truncate(value, variable, degree)
    )


def periodized_sech(
    index: int,
    length: int,
    variable: sp.Symbol,
    degree: int,
) -> sp.Expr:
    """Return the exact periodized-sech series needed at one level."""

    maximum_frequency = degree + index
    period_bound = ceil(maximum_frequency / (2 * length)) + 1
    total = sp.Integer(0)
    for period in range(-period_bound, period_bound + 1):
        frequency = abs(index + 2 * length * period)
        if frequency > maximum_frequency:
            continue
        if frequency == 0:
            term = sp.Integer(1)
        else:
            term = (
                2
                * variable**frequency
                / (1 + variable ** (2 * frequency))
            )
        total += term
    return truncate(total, variable, maximum_frequency)


def exact_axis_data() -> dict[str, object]:
    """Construct and verify the length-five exact axis jet."""

    c = sp.symbols("c", real=True)
    dimension = LENGTH + 1
    base_metric = sp.diag(1, *([2] * (LENGTH - 1)), 4)
    metric_root = sp.diag(
        1,
        *([sp.sqrt(2)] * (LENGTH - 1)),
        2,
    )
    partial = sp.zeros(dimension)
    for index in range(LENGTH):
        partial[index, index + 1] = 1
    operator = metric_root.inv() * partial * metric_root
    pencil = operator + c * operator.T

    pullback = sp.zeros(dimension)
    for scalar_index, coefficient_series in enumerate(
        direct_map_coefficients(DEGREE, DEGREE + 1)
    ):
        coefficient = sum(
            sp.Rational(value.numerator, value.denominator) * c**index
            for index, value in enumerate(
                coefficient_series.coefficients
            )
        )
        pullback += coefficient * pencil ** (2 * scalar_index + 1)
    pullback = truncate_matrix(pullback, c, DEGREE)

    sech_sums = [
        periodized_sech(index, LENGTH, c, DEGREE)
        for index in range(dimension)
    ]
    metric_diagonal = [
        truncate(
            sech_sums[index] / (sech_sums[0] * c**index),
            c,
            DEGREE,
        )
        for index in range(dimension)
    ]
    metric = sp.diag(*metric_diagonal)
    defect_gram = truncate_matrix(
        metric - pullback.T * metric * pullback,
        c,
        DEGREE,
    )
    first_entry = truncate(
        sp.sqrt(defect_gram[0, 0]),
        c,
        DEGREE,
    )
    defect = sp.Matrix(
        [
            truncate(
                defect_gram[index, 0] / first_entry,
                c,
                DEGREE,
            )
            for index in range(dimension)
        ]
    )
    balanced_defect = truncate_matrix(
        metric_root.inv() * defect,
        c,
        DEGREE,
    )
    rank_one_residual = truncate_matrix(
        defect_gram - defect * defect.T,
        c,
        DEGREE,
    )

    right = sp.eye(dimension)[:, 0]
    zero_reflection = (
        right
        + 2
        * sum(
            (
                (-c) ** orbit
                * (partial.T ** (2 * orbit))
                * right
                for orbit in range(1, LENGTH + 1)
            ),
            sp.zeros(dimension, 1),
        )
    ) / (1 + 2 * c**2)
    zero_reflection = truncate_matrix(
        zero_reflection,
        c,
        DEGREE,
    )
    difference = truncate_matrix(
        balanced_defect - zero_reflection,
        c,
        DEGREE,
    )

    expected_leading_difference = sp.zeros(dimension, 1)
    expected_leading_difference[4] = -2 * c**4
    expected_leading_difference[2] = 2 * c**5
    leading_error = truncate_matrix(
        difference - expected_leading_difference,
        c,
        5,
    )
    verified = bool(
        rank_one_residual == sp.zeros(dimension)
        and leading_error == sp.zeros(dimension, 1)
    )
    if not verified:
        raise RuntimeError(
            "the exact length-five axis gauge audit failed"
        )

    return {
        "length": LENGTH,
        "degree": DEGREE,
        "metric_diagonal": [
            str(value) for value in metric_diagonal
        ],
        "balanced_defect": [
            str(value) for value in balanced_defect
        ],
        "zero_reflection_difference": [
            str(value) for value in difference
        ],
        "leading_difference": (
            "-2*c**4*(S @ W @ B5) "
            "+ 2*c**5*(S**3 @ W @ B5)"
        ),
        "rank_one_stein_residual_exactly_zero": True,
        "leading_difference_verified": True,
        "all_checks_passed": verified,
    }


def write_record(record: dict[str, object], output: Path) -> None:
    """Write one deterministic JSON record atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(record, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_axis_gauge_obstruction_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and print the exact axis-gauge audit."""

    args = parse_args()
    record = exact_axis_data()
    write_record(record, args.output)
    print(json.dumps(record, sort_keys=True))


if __name__ == "__main__":
    main()
