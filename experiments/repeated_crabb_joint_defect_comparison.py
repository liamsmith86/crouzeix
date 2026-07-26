#!/usr/bin/env python3
"""Audit the metric comparison used for the joint scalar defect."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class JointDefectComparisonRecord:
    """One exact subspace-distance comparison."""

    ambient_dimension: int
    subspace_dimension: int
    disk_square: str
    actual_distance_square: str
    anchor_distance_square: str
    actual_to_anchor_residual: str
    anchor_to_actual_residual: str
    completion_residual: str
    all_checks_passed: bool


def rational_vector(dimension: int, seed: int) -> sp.Matrix:
    """Return a deterministic rational column."""

    return sp.Matrix(
        [
            sp.Rational(
                ((seed + 2) * (index + 3)) % 11 - 5,
                seed + index + 7,
            )
            for index in range(dimension)
        ]
    )


def squared_distance(
    vector: sp.Matrix,
    subspace_dimension: int,
) -> sp.Expr:
    """Return distance squared to the first coordinate subspace."""

    return sp.simplify(
        sum(
            vector[index] ** 2
            for index in range(subspace_dimension, vector.rows)
        )
    )


def audit_case(
    ambient_dimension: int,
    subspace_dimension: int,
    seed: int,
) -> JointDefectComparisonRecord:
    """Check the two Young bounds and the small-neighbourhood cost."""

    anchor = rational_vector(ambient_dimension, seed)
    perturbation = rational_vector(ambient_dimension, seed + 13)
    scale = sp.Rational(1, seed + 17)
    error = sp.simplify(scale * perturbation)
    actual = sp.simplify(anchor + error)

    disk_square = sp.simplify(error.dot(error))
    actual_distance = squared_distance(actual, subspace_dimension)
    anchor_distance = squared_distance(anchor, subspace_dimension)

    actual_to_anchor = sp.simplify(
        2 * anchor_distance + 2 * disk_square - actual_distance
    )
    anchor_to_actual = sp.simplify(
        2 * actual_distance + 2 * disk_square - anchor_distance
    )

    joint_square = sp.simplify(anchor_distance + disk_square)
    radius = sp.Rational(1, 2)
    scaled_joint_square = sp.simplify(
        joint_square / (1 + 4 * joint_square)
    )
    completion_residual = sp.simplify(
        radius * scaled_joint_square - scaled_joint_square**2
    )

    verified = bool(
        actual_to_anchor >= 0
        and anchor_to_actual >= 0
        and completion_residual >= 0
    )
    if not verified:
        raise RuntimeError("joint defect comparison audit failed")

    return JointDefectComparisonRecord(
        ambient_dimension=ambient_dimension,
        subspace_dimension=subspace_dimension,
        disk_square=str(disk_square),
        actual_distance_square=str(actual_distance),
        anchor_distance_square=str(anchor_distance),
        actual_to_anchor_residual=str(actual_to_anchor),
        anchor_to_actual_residual=str(anchor_to_actual),
        completion_residual=str(completion_residual),
        all_checks_passed=verified,
    )


def standard_records() -> list[JointDefectComparisonRecord]:
    """Return the deterministic exact audit grid."""

    return [
        audit_case(ambient_dimension, subspace_dimension, seed)
        for ambient_dimension in range(3, 7)
        for subspace_dimension in range(1, ambient_dimension)
        for seed in (2, 5)
    ]


def write_records(
    records: list[JointDefectComparisonRecord],
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
            "experiments/repeated_crabb_joint_defect_comparison_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the exact comparison audit."""

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
