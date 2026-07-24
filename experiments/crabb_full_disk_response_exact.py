#!/usr/bin/env python3
"""Regenerate L177's complete quadratic full-disk response formula.

The exact characteristic-Blaschke/Riemann series engine from L173 is
evaluated on every real and imaginary matrix unit.  The resulting
quadratic response is compared with twice L176's curvature covector:

* nested interval-flux covectors in every active paired circle mode;
* the exact weighted-shift covector in grade zero; and
* zero in every inactive mode.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_quadratic_exact import (
    characteristic_series,
    deterministic_disk_direction,
    directional_gradient_series,
    disk_model_series,
)


@dataclass(frozen=True)
class FullDiskResponseRecord:
    """One exact complete-coordinate response audit."""

    dimension: int
    length: int
    disk_direction: tuple[str, ...]
    active_paired_coordinate_count: int
    inactive_coordinate_count: int
    terminal_coordinate_count: int
    paired_response_verified: bool
    inactive_response_verified: bool
    terminal_response_verified: bool


def coefficient_to_physical(direction: sp.Matrix) -> sp.Matrix:
    """Map a coefficient-gauge direction to physical Crabb gauge."""

    dimension = direction.rows
    square_root = [
        1 / sp.sqrt(2),
        *([sp.Integer(1)] * (dimension - 2)),
        1 / sp.sqrt(2),
    ]
    return sp.Matrix(
        dimension,
        dimension,
        lambda row, column: (
            square_root[row]
            * direction[row, column]
            / square_root[column]
        ),
    )


def reduced_mode_vector(
    physical: sp.Matrix,
    length: int,
    mode: int,
) -> sp.Matrix:
    """Return L65's paired reduced vector exactly."""

    size = length - mode
    positive = [
        sp.conjugate(physical[index, index + mode + 1])
        for index in range(size)
    ]
    negative = [
        physical[index + mode - 1, index]
        for index in range(size + 2)
    ]
    if size == 1:
        return sp.Matrix(
            [
                positive[0]
                + sp.sqrt(2) * negative[0]
                + 2 * negative[1]
                + sp.sqrt(2) * negative[2]
            ]
        )

    endpoint = 1 / sp.sqrt(2) if mode == 1 else sp.Integer(1)
    reduced = [
        (
            positive[0]
            + endpoint * negative[0]
            + sp.sqrt(2) * negative[1]
        )
    ]
    reduced.extend(
        positive[index] + negative[index + 1]
        for index in range(1, size - 1)
    )
    reduced.append(
        positive[-1]
        + sp.sqrt(2) * negative[-2]
        + endpoint * negative[-1]
    )
    return sp.Matrix(reduced)


def plucker_matrix(direction: tuple[sp.Expr, ...]) -> sp.Matrix:
    """Return ``z wedge J conjugate(z)`` on nonconstant coordinates."""

    coefficients = sp.Matrix(direction[1:])
    count = len(coefficients)
    reversal = sp.Matrix(
        [sp.conjugate(coefficients[count - 1 - index]) for index in range(count)]
    )
    return coefficients * reversal.T - reversal * coefficients.T


def paired_inverse_image(
    direction: tuple[sp.Expr, ...],
    mode: int,
) -> sp.Matrix:
    """Return ``K_mode s(W)`` from L176's interval formula."""

    length = len(direction)
    anti_diagonal = length - 2 - mode
    size = length - mode
    plucker = plucker_matrix(direction)
    null = sp.Matrix(
        [
            1 / sp.sqrt(2),
            *([sp.Integer(1)] * (size - 2)),
            1 / sp.sqrt(2),
        ]
    )
    normalizer = sp.Matrix(
        [
            (mode + 2) / sp.sqrt(2),
            *([sp.Integer(1)] * (size - 2)),
            (mode + 2) / sp.sqrt(2),
        ]
    )
    image = sp.zeros(size, 1)
    for left in range((anti_diagonal + 1) // 2):
        right = anti_diagonal - left
        indicator = sp.Matrix(
            [
                1
                if left + 1 <= row <= anti_diagonal - left
                else 0
                for row in range(size)
            ]
        )
        coefficient = (
            (normalizer.T * indicator)[0]
            / (normalizer.T * null)[0]
        )
        inverse_image = 16 * (indicator - coefficient * null)
        image += plucker[left, right] * inverse_image
    return image


def predicted_paired_response(
    raw_direction: sp.Matrix,
    disk_direction: tuple[sp.Expr, ...],
    mode: int,
) -> sp.Expr:
    """Return L177's paired response on one raw direction."""

    reduced = reduced_mode_vector(
        coefficient_to_physical(raw_direction),
        len(disk_direction),
        mode,
    )
    inverse_image = paired_inverse_image(disk_direction, mode)
    return sp.simplify(
        2 * sp.re((sp.conjugate(inverse_image).T * reduced)[0])
    )


def terminal_correction(
    direction: tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Return L176's terminal diagonal Hermitian correction."""

    length = len(direction)
    anti_diagonal = length - 2
    plucker = plucker_matrix(direction)
    correction = sp.zeros(length, 1)
    for left in range((anti_diagonal + 1) // 2):
        right = anti_diagonal - left
        for row in range(length):
            pulse = 2 * (
                1
                if left + 1 <= row <= anti_diagonal - left
                else 0
            )
            mean = sp.Rational(
                2 * (anti_diagonal - 2 * left),
                anti_diagonal + 2,
            )
            correction[row] += (
                pulse - mean
            ) * plucker[left, right]
    return sp.simplify(correction)


def weight_tangent(length: int) -> tuple[sp.Matrix, sp.Matrix]:
    """Return the diagonal-H to physical-weight tangent and Crabb weights."""

    tangent = sp.zeros(length)
    weights = []
    for row in range(length):
        weight = (
            sp.sqrt(2)
            if row in (0, length - 1)
            else sp.Integer(1)
        )
        weights.append(weight)
        for column in range(length):
            value = 2 * int(column == row)
            for metric_index in (row, row + 1):
                metric_value = (
                    sp.Rational(1, 2)
                    if metric_index in (0, length)
                    else sp.Integer(1)
                )
                metric_derivative = (
                    int(column == metric_index)
                    + int(column == metric_index - 1)
                )
                value -= (
                    sp.Rational(1, 2)
                    * metric_derivative
                    / metric_value
                )
            tangent[row, column] = sp.simplify(weight * value)
    return tangent, sp.Matrix(weights)


def terminal_curvature_covector(
    direction: tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Return the physical grade-zero curvature covector."""

    length = len(direction)
    correction = terminal_correction(direction)
    tangent, weights = weight_tangent(length)
    difference = sp.zeros(length - 1, length)
    for row in range(length - 1):
        difference[row, row] = -1
        difference[row, row + 1] = 1
    right_side = 4 * difference.T * difference * correction
    unknowns = sp.Matrix(sp.symbols(f"q0:{length}"))
    equations = list(tangent.T * unknowns - right_side)
    equations.append((weights.T * unknowns)[0])
    solution = sp.solve(
        equations,
        list(unknowns),
        dict=True,
    )
    if len(solution) != 1:
        raise RuntimeError("the terminal curvature covector is not unique")
    return sp.Matrix([solution[0][unknown] for unknown in unknowns])


def terminal_physical_weight(
    raw_direction: sp.Matrix,
    row: int,
) -> sp.Expr:
    """Return one physical first-superdiagonal variation."""

    physical = coefficient_to_physical(raw_direction)
    return physical[row, row + 1]


def audit_length(length: int) -> FullDiskResponseRecord:
    """Audit every real ambient coordinate in one length."""

    disk_direction = deterministic_disk_direction(length)
    operator, metric = disk_model_series(disk_direction)
    characteristic = characteristic_series(operator)
    dimension = length + 1
    terminal_covector = terminal_curvature_covector(disk_direction)
    active_count = 0
    inactive_count = 0
    terminal_count = 0
    paired_verified = True
    inactive_verified = True
    terminal_verified = True

    for row in range(dimension):
        for column in range(dimension):
            grade = column - row - 1
            for scalar in (sp.Integer(1), sp.I):
                raw = sp.zeros(dimension)
                raw[row, column] = scalar
                actual = sp.simplify(
                    directional_gradient_series(
                        operator,
                        metric,
                        characteristic,
                        raw,
                    )[2]
                )
                if grade == 0:
                    predicted = sp.simplify(
                        2
                        * sp.re(
                            terminal_covector[row]
                            * terminal_physical_weight(raw, row)
                        )
                    )
                    terminal_count += 1
                    terminal_verified &= sp.simplify(actual - predicted) == 0
                    continue

                mode = abs(grade)
                if 1 <= mode <= length - 3:
                    predicted = predicted_paired_response(
                        raw,
                        disk_direction,
                        mode,
                    )
                    active_count += 1
                    paired_verified &= sp.simplify(actual - predicted) == 0
                else:
                    inactive_count += 1
                    inactive_verified &= actual == 0

    if not (paired_verified and inactive_verified and terminal_verified):
        raise RuntimeError(
            f"the complete response audit failed in length {length}"
        )
    return FullDiskResponseRecord(
        dimension=dimension,
        length=length,
        disk_direction=tuple(str(value) for value in disk_direction[1:]),
        active_paired_coordinate_count=active_count,
        inactive_coordinate_count=inactive_count,
        terminal_coordinate_count=terminal_count,
        paired_response_verified=paired_verified,
        inactive_response_verified=inactive_verified,
        terminal_response_verified=terminal_verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=6)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the complete exact response audit."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("require 3 <= minimum length <= maximum length")
    records = [
        audit_length(length)
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        )
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
