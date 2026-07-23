#!/usr/bin/env python3
"""Audit L137's Fourier Newton edge for the elliptic-axis defect."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import (
    axis_metric_and_defect,
    operator_expansion_from_coefficients,
)


DEFAULT_MAXIMUM_LENGTH = 12
DEFAULT_ORDER = 8


@dataclass(frozen=True)
class AxisDefectFourierRecord:
    """One exact finite axis-defect audit."""

    length: int
    dimension: int
    audit_order: int
    checked_even_modes: int
    odd_modes_vanish: bool
    lower_coefficients_vanish: bool
    leading_coefficients_match: bool
    parity_gap_holds: bool


def make_record(length: int, order: int) -> AxisDefectFourierRecord:
    """Construct and validate one exact defect series."""

    operator, _ = operator_expansion_from_coefficients(
        length + 1,
        [0] * (length - 1),
        order,
    )
    _, defect, _ = axis_metric_and_defect(
        length + 1,
        operator[0],
        order,
    )
    odd_modes_vanish = all(
        defect[index].valuation() == order
        for index in range(1, length + 1, 2)
    )
    checked = tuple(
        first_half_grade
        for first_half_grade in range(1, length // 2 + 1)
        if (
            2 * first_half_grade < length
            and first_half_grade + 1 < order
        )
    )
    lower_vanish = all(
        defect[2 * grade].coefficient(degree) == 0
        for grade in checked
        for degree in range(grade)
    )
    leading_match = all(
        defect[2 * grade].coefficient(grade)
        == 4 * (-1) ** grade
        for grade in checked
    )
    parity_gap = all(
        defect[2 * grade].coefficient(grade + 1) == 0
        for grade in checked
    )
    if not all(
        (
            odd_modes_vanish,
            lower_vanish,
            leading_match,
            parity_gap,
        )
    ):
        raise AssertionError(
            "the axis defect violated its reciprocal-dn Fourier edge"
        )
    return AxisDefectFourierRecord(
        length=length,
        dimension=length + 1,
        audit_order=order,
        checked_even_modes=len(checked),
        odd_modes_vanish=odd_modes_vanish,
        lower_coefficients_vanish=lower_vanish,
        leading_coefficients_match=leading_match,
        parity_gap_holds=parity_gap,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--maximum-length",
        type=int,
        default=DEFAULT_MAXIMUM_LENGTH,
    )
    parser.add_argument("--order", type=int, default=DEFAULT_ORDER)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact finite Fourier regression."""

    args = parse_args()
    if args.maximum_length < 3 or args.order < 4:
        raise ValueError("require maximum length >= 3 and order >= 4")
    records = [
        make_record(length, args.order)
        for length in range(3, args.maximum_length + 1)
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
