#!/usr/bin/env python3
"""Audit the joint disk-endpoint and scalar-channel gap."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class JointDiskGapRecord:
    """One exact two-variable endpoint-angle audit."""

    beta_lower: str
    beta_upper: str
    top_cosine: str
    input_cosine: str
    matrix_trace: str
    matrix_determinant: str
    lower_bound_residual: str
    disk_allocation_residual: str
    all_checks_passed: bool


def pythagorean_pairs() -> list[tuple[Fraction, Fraction]]:
    """Return rational unit-circle pairs in the first quadrant."""

    return [
        (Fraction(1), Fraction(0)),
        (Fraction(4, 5), Fraction(3, 5)),
        (Fraction(12, 13), Fraction(5, 13)),
        (Fraction(20, 29), Fraction(21, 29)),
        (Fraction(0), Fraction(1)),
    ]


def audit_case(
    beta_lower: Fraction,
    beta_upper: Fraction,
    angle_pair: tuple[Fraction, Fraction],
    input_pair: tuple[Fraction, Fraction],
) -> JointDiskGapRecord:
    """Verify one exact spectral-gap and angle specialization."""

    t, t_perp = angle_pair
    a, b = input_pair
    q = beta_lower / (1 + beta_lower)
    mu = beta_upper / (4 - beta_upper)

    matrix = sp.Matrix(
        [
            [q * t_perp**2, -q * t * t_perp],
            [
                -q * t * t_perp,
                q * t**2 + mu / 2,
            ],
        ]
    )
    expected_trace = q + mu / 2
    expected_determinant = q * mu * t_perp**2 / 2
    trace_residual = sp.simplify(matrix.trace() - expected_trace)
    determinant_residual = sp.simplify(
        matrix.det() - expected_determinant
    )

    angle_value = sp.simplify(
        (sp.Matrix([a, b]).T * matrix * sp.Matrix([a, b]))[0]
    )
    angle_constant = q * mu / (2 * q + mu)
    lower_residual = sp.simplify(
        angle_value - angle_constant * t_perp**2
    )

    disk = Fraction(7, 20)
    disk_maximum = Fraction(1, 2)
    disk_constant = min(Fraction(1), mu / (2 * disk_maximum))
    disk_value = disk * a**2 + mu * b**2 / 2
    disk_residual = sp.simplify(
        disk_value - disk_constant * disk
    )

    verified = bool(
        trace_residual == 0
        and determinant_residual == 0
        and lower_residual >= 0
        and disk_residual >= 0
    )
    if not verified:
        raise RuntimeError("joint disk/channel gap audit failed")

    return JointDiskGapRecord(
        beta_lower=str(beta_lower),
        beta_upper=str(beta_upper),
        top_cosine=str(t),
        input_cosine=str(a),
        matrix_trace=str(matrix.trace()),
        matrix_determinant=str(matrix.det()),
        lower_bound_residual=str(lower_residual),
        disk_allocation_residual=str(disk_residual),
        all_checks_passed=verified,
    )


def standard_records() -> list[JointDiskGapRecord]:
    """Return the deterministic exact audit grid."""

    return [
        audit_case(beta_lower, beta_upper, angle, input_pair)
        for beta_lower in (Fraction(1, 3), Fraction(1, 2), Fraction(1))
        for beta_upper in (Fraction(1, 3), Fraction(1, 2), Fraction(1))
        for angle in pythagorean_pairs()
        for input_pair in pythagorean_pairs()
    ]


def write_records(
    records: list[JointDiskGapRecord],
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
            "experiments/"
            "repeated_crabb_joint_disk_scalar_gap_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic audit."""

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
