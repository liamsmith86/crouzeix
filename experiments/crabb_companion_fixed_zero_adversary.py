#!/usr/bin/env python3
"""Falsify a fixed-zero shortcut for the circular-normal telescope.

The grade-two normalized normal variation preserves the determinant
through its target face, suggesting that every grade at least two might
retain a zero eigenvalue.  Grade three disproves that suggestion exactly:
its determinant derivative is already nonzero below the condition face,
although L163's checked condition derivative still vanishes.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
from pathlib import Path
from typing import Sequence

import sympy as sp

from crabb_circular_normal_series import (
    physical_reflected_path,
    real_circular_normal_direction,
)
from general_crabb_weighted_series import inverse_riemann_series


CASES = ((4, 1), (5, 2), (7, 3))


@dataclass(frozen=True)
class FixedZeroRecord:
    """One exact determinant-variation audit."""

    dimension: int
    reflected_grade: int
    target_weight: int
    leading_characteristic_derivative: str
    nonzero_determinant_derivatives: tuple[tuple[int, str], ...]
    preserves_determinant_through_target: bool


def scalar_series_product(
    left: Sequence[sp.Expr],
    right: Sequence[sp.Expr],
    order: int,
) -> list[sp.Expr]:
    """Multiply two scalar series through ``order``."""

    product = [sp.Integer(0) for _ in range(order + 1)]
    for left_degree, left_value in enumerate(left):
        if left_value == 0:
            continue
        maximum_right = order - left_degree
        for right_degree, right_value in enumerate(
            right[: maximum_right + 1]
        ):
            if right_value != 0:
                product[left_degree + right_degree] += (
                    left_value * right_value
                )
    return [sp.expand(value) for value in product]


def permutation_sign(permutation: Sequence[int]) -> int:
    """Return the sign of a permutation."""

    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant_derivative_series(
    matrix: Sequence[sp.Matrix],
    variation: Sequence[sp.Matrix],
) -> list[sp.Expr]:
    """Return ``D det(matrix)[variation]`` as a truncated series.

    Leibniz differentiation avoids a symbolic determinant in the
    bookkeeping parameter and keeps the quadratic-field calculation
    exact.
    """

    order = min(len(matrix), len(variation)) - 1
    dimension = matrix[0].rows
    result = [sp.Integer(0) for _ in range(order + 1)]

    for permutation in itertools.permutations(range(dimension)):
        sign = permutation_sign(permutation)
        for varied_row in range(dimension):
            product = [
                sp.Integer(1),
                *[sp.Integer(0) for _ in range(order)],
            ]
            for row in range(dimension):
                source = variation if row == varied_row else matrix
                entry = [
                    source[degree][row, permutation[row]]
                    for degree in range(order + 1)
                ]
                if not any(entry):
                    break
                product = scalar_series_product(product, entry, order)
            else:
                for degree, value in enumerate(product):
                    result[degree] += sign * value

    return [sp.factor(value) for value in result]


def make_record(dimension: int, grade: int) -> FixedZeroRecord:
    """Construct one exact normalized normal variation."""

    length = dimension - 1
    face_degree = grade + 1
    target_weight = 2 * face_degree
    normal = real_circular_normal_direction(
        dimension,
        length + 2 - grade,
    )

    operators = []
    for strong_value in (0, 1, -1):
        path = physical_reflected_path(
            dimension=dimension,
            equality_grade=grade,
            strong_parameter=sp.Integer(strong_value),
            strong_direction=normal,
            strong_degree=face_degree,
            order=target_weight,
        )
        _, operator = inverse_riemann_series(path, target_weight)
        operators.append(operator)

    base = operators[0]
    variation = [
        sp.simplify((plus - minus) / 2)
        for plus, minus in zip(
            operators[1],
            operators[2],
            strict=True,
        )
    ]
    derivative = determinant_derivative_series(base, variation)
    spectral_parameter = sp.symbols("spectral_parameter")
    characteristic_derivative = sp.factor(
        -sp.trace(
            (
                spectral_parameter * sp.eye(dimension)
                - base[0]
            ).adjugate()
            * variation[face_degree]
        )
    )
    predicted_characteristic_derivative = (
        -grade * spectral_parameter ** (grade - 1)
    )
    if (
        sp.simplify(
            characteristic_derivative
            - predicted_characteristic_derivative
        )
        != 0
    ):
        raise AssertionError(
            "the leading characteristic monomial changed"
        )
    nonzero = tuple(
        (degree, str(value))
        for degree, value in enumerate(derivative)
        if value != 0
    )

    expected = {
        (4, 1): (
            (2, "-1"),
            (4, "-6*(-3 + 2*sqrt(2))"),
        ),
        (5, 2): (),
        (7, 3): (
            (5, "-4"),
            (7, "-2*(-17 + 8*sqrt(2))"),
            (8, "-4*(-17 + 12*sqrt(2))"),
        ),
    }[(dimension, grade)]
    if nonzero != expected:
        raise AssertionError(
            "the determinant-variation adversary changed"
        )

    return FixedZeroRecord(
        dimension=dimension,
        reflected_grade=grade,
        target_weight=target_weight,
        leading_characteristic_derivative=str(
            characteristic_derivative
        ),
        nonzero_determinant_derivatives=nonzero,
        preserves_determinant_through_target=not nonzero,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact adversarial audit."""

    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for dimension, grade in CASES:
            record = make_record(dimension, grade)
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
