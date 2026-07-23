#!/usr/bin/env python3
"""Audit mixed-grade polarization of the corrected Blaschke face.

For two distinct low representatives ``k<ell<=L/2``, compare the exact
quadratic norm series along the combined direction with the sum of the
two one-grade series.  L144 predicts that their difference vanishes
through the complete principal range, in particular at the only
possible mixed face ``c^(k+ell)``.

The computation reuses L143's rational Weierstrass preparation and
simple generalized singular-pair lift.  It includes the two possible
terminal resonances ``L=2k+ell`` and ``L=k+2ell``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import cache
import json
from pathlib import Path

from crabb_faber_blaschke_formal import formal_dual_calculation
from exact_truncated_series import Series


DEFAULT_CASES = (
    (5, 1, 2),
    (6, 1, 2),
    (7, 1, 2),
    (7, 1, 3),
    (7, 2, 3),
    (8, 2, 3),
)


@dataclass(frozen=True)
class MixedFaberBlaschkeRecord:
    """One exact mixed-grade cancellation audit."""

    length: int
    dimension: int
    first_grade: int
    second_grade: int
    target_degree: int
    audit_degree: int
    cross_coefficient: str
    cross_series_vanishes: bool
    first_diagonal_coefficient: str
    second_diagonal_coefficient: str
    left_alias_resonance: bool
    right_alias_resonance: bool


@cache
def quadratic_series(
    length: int,
    grades: tuple[int, ...],
    audit_degree: int,
) -> Series:
    """Return one cached exact quadratic generalized singular series."""

    return formal_dual_calculation(
        length,
        grades,
        audit_degree,
    ).quadratic_eigenvalue


def make_record(
    length: int,
    first_grade: int,
    second_grade: int,
) -> MixedFaberBlaschkeRecord:
    """Construct and validate one mixed-grade polarization."""

    if not (
        1
        <= first_grade
        < second_grade
        <= length / 2
    ):
        raise ValueError("require 1 <= k < ell <= L/2")
    target_degree = first_grade + second_grade
    audit_degree = 2 * second_grade
    combined = quadratic_series(
        length,
        (first_grade, second_grade),
        audit_degree,
    )
    first = quadratic_series(
        length,
        (first_grade,),
        audit_degree,
    )
    second = quadratic_series(
        length,
        (second_grade,),
        audit_degree,
    )
    cross = combined - first - second
    cross_vanishes = cross.valuation() == cross.order
    if (
        not cross_vanishes
        or combined.coefficient(2 * first_grade) != -64
        or combined.coefficient(2 * second_grade) != -64
    ):
        raise AssertionError("mixed-grade diagonalization failed")

    return MixedFaberBlaschkeRecord(
        length=length,
        dimension=length + 1,
        first_grade=first_grade,
        second_grade=second_grade,
        target_degree=target_degree,
        audit_degree=audit_degree,
        cross_coefficient=str(cross.coefficient(target_degree)),
        cross_series_vanishes=cross_vanishes,
        first_diagonal_coefficient=str(
            combined.coefficient(2 * first_grade)
        ),
        second_diagonal_coefficient=str(
            combined.coefficient(2 * second_grade)
        ),
        left_alias_resonance=(
            length == 2 * first_grade + second_grade
        ),
        right_alias_resonance=(
            length == first_grade + 2 * second_grade
        ),
    )


def parse_cases(values: list[str]) -> tuple[tuple[int, int, int], ...]:
    """Parse repeated ``L,k,ell`` values."""

    if not values:
        return DEFAULT_CASES
    cases = []
    for value in values:
        length, first, second = value.split(",", maxsplit=2)
        cases.append((int(length), int(first), int(second)))
    return tuple(cases)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="L,K,ELL",
        help="length and two distinct low grades; may be repeated",
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the selected exact mixed-grade audits."""

    args = parse_args()
    records = [
        make_record(length, first, second)
        for length, first, second in parse_cases(args.case)
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    print("\n".join(lines), flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
