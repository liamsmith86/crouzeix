#!/usr/bin/env python3
"""Regenerate the grade-three covariant-recurrence obstruction.

The complete grade-three physical cross vanishes in the finite checks,
but its prepared matrix jet has a nonzero ordinary L65 Hessian pairing
with the eligible circular normal.  Thus higher mixed normalization
jets are essential; a proof cannot retain only the first matrix column.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_series import (
    physical_reflected_path,
    real_circular_normal_direction,
)
from general_crabb_weighted_series import inverse_riemann_series


GRADE = 3
WEIGHT = GRADE + 1


@dataclass(frozen=True)
class CovariantGateRecord:
    """One exact grade-three base-Hessian audit."""

    dimension: int
    length: int
    normal_mode: int
    prepared_reduced_coordinate: str
    normal_reduced_coordinate: str
    base_hessian_cross_coefficient: str
    predicted_cross_coefficient: str


def reduced_coordinate(matrix: sp.Matrix, mode: int) -> sp.Expr:
    """Return L65's length-one paired-mode coordinate."""

    entries = sp.Matrix(
        [
            sp.conjugate(matrix[0, mode + 1]),
            matrix[mode - 1, 0],
            matrix[mode, 1],
            matrix[mode + 1, 2],
        ]
    )
    aggregator = sp.Matrix([[1, sp.sqrt(2), 2, sp.sqrt(2)]])
    return sp.factor((aggregator * entries)[0])


def make_record(dimension: int) -> CovariantGateRecord:
    """Derive the effective compact jet and its exact L65 pairing."""

    if dimension < 7:
        raise ValueError("grade three requires dimension at least seven")

    length = dimension - 1
    mode = length - 1
    amplitude = sp.symbols("amplitude", real=True)
    zero = sp.zeros(dimension)
    path = physical_reflected_path(
        dimension=dimension,
        equality_grade=GRADE,
        strong_parameter=sp.Integer(0),
        strong_direction=zero,
        strong_degree=1,
        order=WEIGHT,
        amplitude=amplitude,
        ellipse=sp.Integer(1),
    )
    _, operator = inverse_riemann_series(path, WEIGHT)
    prepared = sp.simplify(
        sp.diff(operator[WEIGHT], amplitude).subs(amplitude, 0)
    )
    normal = real_circular_normal_direction(dimension, mode)
    prepared_reduced = reduced_coordinate(prepared, mode)
    normal_reduced = reduced_coordinate(normal, mode)

    kernel = sp.Rational(
        mode * (mode - 1) * (mode - 2),
        6 * (mode + 1) ** 2,
    )
    cross = sp.factor(
        -2 * prepared_reduced * kernel * normal_reduced
    )
    predicted_prepared = (
        sp.Integer(-4)
        if length == 10
        else 4 * (-2 + sp.sqrt(2))
    )
    predicted = sp.factor(
        -2 * predicted_prepared * kernel * sp.Rational(9, 2)
    )

    if sp.simplify(
        prepared_reduced - predicted_prepared
    ) != 0:
        raise AssertionError("the prepared grade-three coordinate changed")
    if sp.simplify(normal_reduced - sp.Rational(9, 2)) != 0:
        raise AssertionError("the grade-three normal coordinate changed")
    if (
        sp.simplify(cross - predicted) != 0
        or sp.N(cross) <= 0
    ):
        raise AssertionError("the nonzero Hessian obstruction changed")

    return CovariantGateRecord(
        dimension=dimension,
        length=length,
        normal_mode=mode,
        prepared_reduced_coordinate=str(prepared_reduced),
        normal_reduced_coordinate=str(normal_reduced),
        base_hessian_cross_coefficient=str(cross),
        predicted_cross_coefficient=str(predicted),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=7)
    parser.add_argument("--maximum-size", type=int, default=12)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic exact audit."""

    args = parse_args()
    if args.minimum_size < 7 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 7 <= minimum <= maximum")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for dimension in range(
            args.minimum_size,
            args.maximum_size + 1,
        ):
            record = make_record(dimension)
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
