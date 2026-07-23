#!/usr/bin/env python3
"""Audit the first weighted transverse gradients of the L118 envelope.

For a strong direction ``Y`` at a ``p x p`` Crabb block, this script forms

    C_p + c C_p^* + s c^(p-1) Y

and differentiates the exact rank-one certificate with respect to ``s``.
The dangerous weighted range ends at degree ``2(p-1)``.  Every tested
derivative vanishes through that degree; two low-size cases also expose the
first nonzero coefficient beyond it.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from general_crabb_weighted_series import weighted_condition_derivative


@dataclass(frozen=True)
class GradientRecord:
    dimension: int
    direction: str
    checked_order: int
    first_nonzero_degree: int | None
    first_nonzero_coefficient: str | None


def coefficient_records(expression: sp.Expr) -> list[tuple[int, sp.Expr]]:
    """Return the nonzero coefficients of a univariate expression."""

    symbols = list(expression.free_symbols)
    if not symbols:
        if expression != 0:
            raise ValueError("a nonzero constant derivative was not expected")
        return []
    if len(symbols) != 1:
        raise ValueError("a univariate derivative series was expected")
    parameter = symbols[0]
    polynomial = sp.Poly(expression, parameter)
    return sorted(
        (
            int(degree[0]),
            sp.simplify(coefficient),
        )
        for degree, coefficient in polynomial.terms()
        if coefficient != 0
    )


def cases() -> list[tuple[int, str, sp.Matrix, int, tuple[int, sp.Expr] | None]]:
    """Return deterministic low-size weighted-gradient cases."""

    dimension_four = sp.zeros(4)
    dimension_four[3, 0] = 1

    dimension_five = sp.zeros(5)
    dimension_five[3, 0] = 1 / sp.sqrt(2)
    dimension_five[4, 1] = 1 / sp.sqrt(2)

    dimension_six_mode_four = sp.zeros(6)
    dimension_six_mode_four[3, 0] = 1

    dimension_six_bottom = sp.zeros(6)
    dimension_six_bottom[5, 0] = 1

    return [
        (4, "bottom_grade_4", dimension_four, 1, (7, sp.Integer(-16))),
        (
            5,
            "penultimate_grade_4",
            dimension_five,
            2,
            (10, sp.Integer(-32)),
        ),
        (6, "paired_grade_4", dimension_six_mode_four, 0, None),
        (6, "bottom_grade_6", dimension_six_bottom, 0, None),
    ]


def run_checks() -> list[GradientRecord]:
    """Run every exact check and return its audit records."""

    records: list[GradientRecord] = []
    for dimension in (3, 4, 5):
        coordinates = sp.symbols(
            f"y0:{dimension * dimension}",
            real=True,
        )
        perturbation = sp.Matrix(dimension, dimension, coordinates)
        derivative = weighted_condition_derivative(dimension, perturbation)
        if derivative != 0:
            raise RuntimeError(
                f"universal weighted cancellation failed in size {dimension}"
            )
        record = GradientRecord(
            dimension=dimension,
            direction="all_real_coordinates",
            checked_order=2 * (dimension - 1),
            first_nonzero_degree=None,
            first_nonzero_coefficient=None,
        )
        records.append(record)
        print(json.dumps(asdict(record), sort_keys=True), flush=True)

    for dimension, name, perturbation, extra_order, expected in cases():
        derivative = weighted_condition_derivative(
            dimension,
            perturbation,
            extra_order=extra_order,
        )
        coefficients = coefficient_records(derivative)
        dangerous_order = 2 * (dimension - 1)
        if any(degree <= dangerous_order for degree, _ in coefficients):
            raise RuntimeError(
                f"dangerous weighted gradient survived for {dimension}:{name}"
            )
        if expected is not None:
            if not coefficients or coefficients[0] != expected:
                raise RuntimeError(
                    f"unexpected first coefficient for {dimension}:{name}"
                )

        first = coefficients[0] if coefficients else None
        record = GradientRecord(
            dimension=dimension,
            direction=name,
            checked_order=dangerous_order + extra_order,
            first_nonzero_degree=first[0] if first else None,
            first_nonzero_coefficient=str(first[1]) if first else None,
        )
        records.append(record)
        print(json.dumps(asdict(record), sort_keys=True), flush=True)
    return records


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSONL path for the regenerated audit records",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = run_checks()
    if args.output is not None:
        args.output.write_text(
            "".join(
                f"{json.dumps(asdict(record), sort_keys=True)}\n"
                for record in records
            ),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
