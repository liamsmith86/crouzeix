#!/usr/bin/env python3
"""Regenerate the circular-normal/defect-column coefficient.

For reflected grade ``k >= 2``, the only eligible circular normal has
mode ``m = L + 2 - k``.  At the Crabb point, its first normalized
matrix jet has a unique optimal rank-one-defect response in coordinate
``m``.  This script regenerates that response exactly and compares it
with the closed all-size formula used by L163.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_series import (
    optimized_defect_jets,
    real_circular_normal_direction,
)
from general_crabb_weighted_series import (
    crabb_matrix,
    inverse_riemann_series,
)
from rank_one_stein_series import stein_gramian_series


@dataclass(frozen=True)
class NormalDefectColumnRecord:
    """One exact circular-normal defect-column audit."""

    dimension: int
    length: int
    reflected_grade: int
    normal_mode: int
    exact_coefficient: str
    predicted_coefficient: str
    residual: str
    stein_column_formula_verified: bool
    other_defect_coordinates_vanish: bool


def predicted_coefficient(length: int, grade: int) -> sp.Expr:
    """Return the closed defect-column coefficient ``gamma_(L,k)``."""

    if grade == 2:
        return sp.Integer(2)
    if grade == 3:
        return sp.sqrt(2) * (13 * length - 9) / (8 * length)
    return sp.sqrt(2) * (
        sp.Rational(2 * grade - 3, 2)
        - sp.Rational((grade - 2) ** 2, length)
    )


def make_record(length: int, grade: int) -> NormalDefectColumnRecord:
    """Regenerate one coefficient from the normalized first-order path."""

    if not 2 <= grade <= length // 2:
        raise ValueError("the grade must satisfy 2 <= k <= floor(L/2)")

    dimension = length + 1
    mode = length + 2 - grade
    epsilon = sp.symbols("epsilon", real=True)
    base = crabb_matrix(dimension)
    normal = real_circular_normal_direction(dimension, mode)
    _, operator = inverse_riemann_series(
        [base, normal, sp.zeros(dimension)],
        2,
    )
    defect = optimized_defect_jets(operator, epsilon, jet_count=1)
    exact = sp.factor(-sp.expand(defect[mode]).coeff(epsilon))
    predicted = sp.factor(predicted_coefficient(length, grade))
    residual = sp.factor(exact - predicted)
    other_coordinates_vanish = all(
        sp.expand(defect[index]).coeff(epsilon) == 0
        for index in range(1, dimension)
        if index != mode
    )

    if residual != 0:
        raise AssertionError("the defect-column formula changed")
    if not other_coordinates_vanish:
        raise AssertionError("the normal activated another defect coordinate")

    base_defect = sp.eye(dimension)[:, 0]
    normal_column = stein_gramian_series(
        operator,
        base_defect,
        epsilon,
        1,
    )[1]
    defect_direction = sp.MutableDenseMatrix(base_defect)
    defect_direction[mode] += epsilon
    defect_column = stein_gramian_series(
        [base, sp.zeros(dimension)],
        defect_direction,
        epsilon,
        1,
    )[1]
    column_residual = sp.simplify(
        normal_column - predicted * defect_column
    )
    predicted_residual = sp.zeros(dimension)
    if grade >= 3:
        if grade == 3:
            edge = sp.Rational(5 * length - 9, 8 * length)
            interior_step = None
        else:
            edge = (
                sp.Rational(2 * grade - 5, 2)
                - sp.Rational((grade - 2) ** 2, length)
            )
            interior_step = (
                2 - sp.Rational(2 * (grade - 2), length)
            )
        predicted_residual[0, mode] = -sp.sqrt(2) * edge
        predicted_residual[grade - 2, length] = (
            2 * sp.sqrt(2) * edge
        )
        if interior_step is not None:
            for index in range(1, grade - 2):
                predicted_residual[index, mode + index] = (
                    (2 * index - (grade - 2)) * interior_step
                )
        predicted_residual += predicted_residual.T
    stein_column_residual = sp.simplify(
        column_residual - predicted_residual
    )
    if stein_column_residual != sp.zeros(dimension):
        raise AssertionError("the residual Stein chain formula changed")

    return NormalDefectColumnRecord(
        dimension=dimension,
        length=length,
        reflected_grade=grade,
        normal_mode=mode,
        exact_coefficient=str(exact),
        predicted_coefficient=str(predicted),
        residual=str(residual),
        stein_column_formula_verified=True,
        other_defect_coordinates_vanish=other_coordinates_vanish,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=4)
    parser.add_argument("--maximum-length", type=int, default=12)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic exact audit."""

    args = parse_args()
    if (
        args.minimum_length < 4
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("length range must satisfy 4 <= minimum <= maximum")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        ):
            for grade in range(2, length // 2 + 1):
                record = make_record(length, grade)
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")


if __name__ == "__main__":
    main()
