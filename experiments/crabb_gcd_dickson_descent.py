#!/usr/bin/env python3
"""Audit exact Dickson descent by the common divisor of length and grade.

Write ``L=dq`` and ``k=ds``.  Degree-``d`` Dickson descent reduces the
residue-zero coordinates of the size-``L+1`` grade-``k`` equality
pencil to the complete size-``q+1`` grade-``s`` pencil at parameter
``c^d``.  Taking ``d=gcd(L,k)`` leaves the reduced pair ``(q,s)``
coprime.

The checker verifies the polynomial pencil, coordinate Gramian,
subcritical compression, and inactive apex fibers exactly.  It is a
regression for the all-size residue proof recorded in L136.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from math import gcd
import json
from pathlib import Path

from crabb_divisible_dickson_descent import (
    audit_residue_descent,
)


@dataclass(frozen=True)
class GcdDescentRecord:
    common_divisor: int
    reduced_length: int
    reduced_grade: int
    length: int
    grade: int
    dimension: int
    outer_dimension: int
    reduced_pair_coprime: bool
    exact_outer_pencil: bool
    polynomial_cross_blocks_vanish: bool
    exact_amplitude_linearity: bool
    coordinate_gramian_reduces: bool
    exact_outer_coordinate_gramian: bool
    subcritical_compression_vanishes: bool
    inactive_apex_shifts: bool


def make_record(
    common_divisor: int,
    reduced_length: int,
    reduced_grade: int,
) -> GcdDescentRecord:
    """Construct and verify one exact common-divisor reduction."""

    if common_divisor < 1:
        raise ValueError("the common divisor must be positive")
    if not 1 <= reduced_grade <= reduced_length // 2:
        raise ValueError(
            "require 1 <= reduced grade <= floor(reduced length / 2)"
        )

    length = common_divisor * reduced_length
    grade = common_divisor * reduced_grade
    audit = audit_residue_descent(
        common_divisor,
        reduced_length,
        reduced_grade,
    )

    record = GcdDescentRecord(
        common_divisor=common_divisor,
        reduced_length=reduced_length,
        reduced_grade=reduced_grade,
        length=length,
        grade=grade,
        dimension=audit.dimension,
        outer_dimension=audit.outer_dimension,
        reduced_pair_coprime=(
            gcd(reduced_length, reduced_grade) == 1
        ),
        exact_outer_pencil=audit.exact_outer_pencil,
        polynomial_cross_blocks_vanish=(
            audit.polynomial_cross_blocks_vanish
        ),
        exact_amplitude_linearity=audit.exact_amplitude_linearity,
        coordinate_gramian_reduces=audit.coordinate_gramian_reduces,
        exact_outer_coordinate_gramian=(
            audit.exact_outer_coordinate_gramian
        ),
        subcritical_compression_vanishes=(
            audit.subcritical_compression_vanishes
        ),
        inactive_apex_shifts=audit.inactive_apex_shifts,
    )
    if not all(
        (
            record.exact_outer_pencil,
            record.polynomial_cross_blocks_vanish,
            record.exact_amplitude_linearity,
            record.coordinate_gramian_reduces,
            record.exact_outer_coordinate_gramian,
            record.subcritical_compression_vanishes,
            record.inactive_apex_shifts,
        )
    ):
        raise AssertionError(
            "the exact common-divisor descent failed: "
            f"{record}"
        )
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-divisor", type=int, default=5)
    parser.add_argument("--maximum-reduced-length", type=int, default=8)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run every coprime reduced pair in the requested grid."""

    args = parse_args()
    if args.maximum_divisor < 2 or args.maximum_reduced_length < 2:
        raise ValueError("invalid divisor or reduced-length bound")

    records = [
        make_record(common_divisor, reduced_length, reduced_grade)
        for common_divisor in range(2, args.maximum_divisor + 1)
        for reduced_length in range(2, args.maximum_reduced_length + 1)
        for reduced_grade in range(1, reduced_length // 2 + 1)
        if gcd(reduced_length, reduced_grade) == 1
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
