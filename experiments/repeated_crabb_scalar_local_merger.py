#!/usr/bin/env python3
"""Audit the scalar master-square used in the repeated local merger."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path


@dataclass(frozen=True)
class ScalarLocalMergerRecord:
    """One exact master-square and retained-margin audit."""

    disk_reserve: str
    normal_curvature: str
    response_constant: str
    joint_defect: str
    normal_coordinate: str
    direct_expression: str
    completed_upper_bound: str
    retained_margin_residual: str
    all_checks_passed: bool


def audit_case(
    disk_reserve: Fraction,
    normal_curvature: Fraction,
    response_constant: Fraction,
    defect_fraction: Fraction,
    normal_coordinate: Fraction,
) -> ScalarLocalMergerRecord:
    """Verify one exact completion below the safe defect radius."""

    safe_radius = (
        disk_reserve
        * normal_curvature
        / (4 * response_constant**2)
    )
    joint_defect = defect_fraction * safe_radius
    response = response_constant * joint_defect

    direct = (
        -disk_reserve * joint_defect
        - normal_curvature * normal_coordinate**2 / 2
        + 2 * normal_coordinate * response
    )
    completed_upper = (
        -disk_reserve * joint_defect
        + 2 * response**2 / normal_curvature
    )
    retained_residual = (
        -disk_reserve * joint_defect / 2
        - completed_upper
    )

    verified = bool(
        direct <= completed_upper
        and completed_upper <= -disk_reserve * joint_defect / 2
        and retained_residual >= 0
    )
    if not verified:
        raise RuntimeError("scalar local merger audit failed")

    return ScalarLocalMergerRecord(
        disk_reserve=str(disk_reserve),
        normal_curvature=str(normal_curvature),
        response_constant=str(response_constant),
        joint_defect=str(joint_defect),
        normal_coordinate=str(normal_coordinate),
        direct_expression=str(direct),
        completed_upper_bound=str(completed_upper),
        retained_margin_residual=str(retained_residual),
        all_checks_passed=verified,
    )


def standard_records() -> list[ScalarLocalMergerRecord]:
    """Return the deterministic exact audit grid."""

    return [
        audit_case(
            disk_reserve,
            normal_curvature,
            response_constant,
            defect_fraction,
            normal_coordinate,
        )
        for disk_reserve in (Fraction(1, 3), Fraction(2, 3))
        for normal_curvature in (Fraction(1, 2), Fraction(3, 2))
        for response_constant in (Fraction(1), Fraction(5, 2))
        for defect_fraction in (Fraction(1, 4), Fraction(1))
        for normal_coordinate in (
            Fraction(-1, 3),
            Fraction(0),
            Fraction(2, 5),
        )
    ]


def write_records(
    records: list[ScalarLocalMergerRecord],
    output: Path,
) -> str:
    """Write JSON Lines atomically and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_scalar_local_merger_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the exact master-square audit."""

    arguments = parse_args()
    records = standard_records()
    digest = write_records(records, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "record_count": len(records),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
