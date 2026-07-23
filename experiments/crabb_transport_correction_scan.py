#!/usr/bin/env python3
"""Falsify coordinate-invariance of the first defect correction.

The completed-square target is invariant only after quotienting defect
scale and lower-weight transport.  This focused exact scan compares the
raw optimized defect with ``2 U(u)d_c`` through weight ``k``.  It
records that the visible one-coordinate correction can move into the
transported part as the terminal length changes, so that raw prefixes
must not be promoted to an all-size formula.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import optimized_hessian
from crabb_principal_face_locality import transported_defect_prefix


DEFAULT_PAIRS = (
    (8, 3),
    (9, 3),
    (10, 3),
    (9, 4),
    (10, 4),
    (11, 4),
    (12, 4),
)


@dataclass(frozen=True)
class TransportCorrectionRecord:
    """One exact raw-coordinate transport comparison."""

    length: int
    grade: int
    audit_order: int
    transport_difference: tuple[tuple[int, int, str], ...]


def make_record(length: int, grade: int) -> TransportCorrectionRecord:
    """Compare optimized and transported defect prefixes exactly."""

    order = grade + 1
    _, optimized, _ = optimized_hessian(length + 1, grade, order)
    transported = transported_defect_prefix(length, grade, order)
    difference = tuple(
        (
            coordinate,
            degree,
            str(
                optimized[coordinate - 1].coefficient(degree)
                - transported[coordinate - 1][degree]
            ),
        )
        for coordinate in range(1, length + 1)
        for degree in range(order)
        if (
            optimized[coordinate - 1].coefficient(degree)
            != transported[coordinate - 1][degree]
        )
    )
    if any(degree < grade for _, degree, _ in difference):
        raise AssertionError(
            "the raw transport comparison changed below its first face"
        )
    return TransportCorrectionRecord(
        length=length,
        grade=grade,
        audit_order=order,
        transport_difference=difference,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic grade-three/four scan."""

    args = parse_args()
    records = [
        make_record(length, grade)
        for length, grade in DEFAULT_PAIRS
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
