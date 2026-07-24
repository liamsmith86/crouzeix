#!/usr/bin/env python3
"""Regenerate L176's exact full-disk Pluecker correction isometry.

The checker uses exact SymPy arithmetic.  For each Pluecker
anti-diagonal it verifies:

* the interval-pulse Hermitian correction maps under the differentiated
  general-H disk chart to the claimed L65 reduced flux vectors;
* those flux vectors have Gram matrix ``64 I`` in every paired L65
  curvature block; and
* the terminal diagonal correction has Gram matrix ``32 I`` in the
  exact weighted-shift curvature.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class CorrectionIsometryRecord:
    """One exact dimension audit."""

    dimension: int
    length: int
    paired_mode_count: int
    tangent_reduction_identity: bool
    paired_curvature_isometry: bool
    diagonal_curvature_isometry: bool


def correction_matrix(anti_diagonal: int) -> sp.Matrix:
    """Return the interval-pulse correction on one anti-diagonal."""

    size = anti_diagonal + 2
    column_count = (anti_diagonal + 1) // 2
    return sp.Matrix(
        size,
        column_count,
        lambda row, column: (
            2
            * (
                1
                if column + 1
                <= row
                <= anti_diagonal - column
                else 0
            )
            - sp.Rational(
                2 * (anti_diagonal - 2 * column),
                anti_diagonal + 2,
            )
        ),
    )


def plucker_correction(
    direction: tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Assemble L176's Hermitian correction for one Toeplitz direction."""

    length = len(direction)
    coefficients = sp.Matrix(direction[1:])
    coefficient_count = len(coefficients)
    reversal = sp.Matrix(
        [
            sp.conjugate(coefficients[coefficient_count - 1 - index])
            for index in range(coefficient_count)
        ]
    )
    plucker = coefficients * reversal.T - reversal * coefficients.T
    correction = sp.zeros(length)
    for anti_diagonal in range(1, coefficient_count):
        offset = coefficient_count - 1 - anti_diagonal
        pulses = correction_matrix(anti_diagonal)
        intrinsic = sp.Matrix(
            [
                sp.conjugate(
                    plucker[left, anti_diagonal - left]
                )
                for left in range(pulses.cols)
            ]
        )
        values = pulses * intrinsic
        for row, value in enumerate(values):
            correction[row, row + offset] = sp.expand(value)
            if offset:
                correction[row + offset, row] = sp.conjugate(
                    correction[row, row + offset]
                )
    return correction


def tangent_reduction(size: int) -> sp.Matrix:
    """Return the physical disk-tangent to L65 reduced-mode map."""

    reduction = sp.zeros(size)
    reduction[0, 0] = sp.sqrt(2)
    reduction[0, 1] = -sp.sqrt(2)
    reduction[-1, -2] = -sp.sqrt(2)
    reduction[-1, -1] = sp.sqrt(2)
    for row in range(1, size - 1):
        reduction[row, row - 1] = -1
        reduction[row, row] = 2
        reduction[row, row + 1] = -1
    return reduction


def interval_indicators(size: int) -> sp.Matrix:
    """Return the nested interval columns behind the correction."""

    column_count = (size - 1) // 2
    return sp.Matrix(
        size,
        column_count,
        lambda row, column: (
            1
            if column + 1 <= row <= size - 2 - column
            else 0
        ),
    )


def crabb_path_curvature(size: int) -> sp.Matrix:
    """Return L65's singular path matrix ``R_size`` exactly."""

    support = sp.zeros(size)
    for index in range(size - 1):
        weight = (
            sp.sqrt(2)
            if index in (0, size - 2)
            else sp.Integer(1)
        )
        support[index, index + 1] = weight / 2
        support[index + 1, index] = weight / 2
    return (sp.eye(size) - support) / 4


def flux_vectors(size: int) -> sp.Matrix:
    """Return the universal L176 flux columns."""

    return 16 * crabb_path_curvature(size) * interval_indicators(size)


def paired_isometry(size: int, mode: int) -> bool:
    """Verify the exact L65 inverse-curvature Gram identity."""

    path = crabb_path_curvature(size)
    indicators = interval_indicators(size)
    flux = 16 * path * indicators
    null = sp.Matrix(
        [
            sp.sqrt(2) / 2,
            *([sp.Integer(1)] * (size - 2)),
            sp.sqrt(2) / 2,
        ]
    )
    normalizer = sp.Matrix(
        [
            (mode + 2) * sp.sqrt(2) / 2,
            *([sp.Integer(1)] * (size - 2)),
            (mode + 2) * sp.sqrt(2) / 2,
        ]
    )
    corrected_columns = []
    for column in range(indicators.cols):
        indicator = indicators[:, column]
        coefficient = (
            (normalizer.T * indicator)[0]
            / (normalizer.T * null)[0]
        )
        corrected = 16 * (indicator - coefficient * null)
        if path * corrected != flux[:, column]:
            return False
        if (normalizer.T * corrected)[0] != 0:
            return False
        corrected_columns.append(corrected)
    inverse_images = sp.Matrix.hstack(*corrected_columns)
    gram = sp.simplify(flux.T * inverse_images)
    return gram == 64 * sp.eye(indicators.cols)


def diagonal_isometry(length: int) -> bool:
    """Verify the terminal weighted-shift curvature exactly."""

    anti_diagonal = length - 2
    correction = correction_matrix(anti_diagonal)
    difference = sp.zeros(length - 1, length)
    for row in range(length - 1):
        difference[row, row] = -1
        difference[row, row + 1] = 1
    gram = sp.simplify(
        4 * correction.T * difference.T * difference * correction
    )
    return gram == 32 * sp.eye(correction.cols)


def audit_length(length: int) -> CorrectionIsometryRecord:
    """Audit all anti-diagonals in one length."""

    tangent_checks = []
    paired_checks = []
    for anti_diagonal in range(1, length - 2):
        size = anti_diagonal + 2
        mode = length - 2 - anti_diagonal
        correction = correction_matrix(anti_diagonal)
        tangent_checks.append(
            tangent_reduction(size) * correction
            == flux_vectors(size)
        )
        paired_checks.append(paired_isometry(size, mode))
    return CorrectionIsometryRecord(
        dimension=length + 1,
        length=length,
        paired_mode_count=max(0, length - 3),
        tangent_reduction_identity=all(tangent_checks),
        paired_curvature_isometry=all(paired_checks),
        diagonal_curvature_isometry=diagonal_isometry(length),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=30)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the exact isometry audit."""

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
    if any(
        not (
            record.tangent_reduction_identity
            and record.paired_curvature_isometry
            and record.diagonal_curvature_isometry
        )
        for record in records
    ):
        raise RuntimeError("the exact correction isometry audit failed")

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
