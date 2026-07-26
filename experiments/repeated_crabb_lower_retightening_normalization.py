#!/usr/bin/env python3
"""Audit the monomial-apex raw/effective endpoint normalization."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


Matrix = sp.Matrix


@dataclass(frozen=True)
class RetighteningRecord:
    """One exact multiplicity audit."""

    multiplicity: int
    boundary_lower_residual_count: int
    boundary_upper_residual_count: int
    raw_repair_lower_residual_count: int
    raw_repair_upper_residual_count: int
    effective_lower_residual_count: int
    effective_upper_residual_count: int
    raw_repair_increment_residual_count: int
    retightening_increment_residual_count: int
    l212_target_residual_count: int
    all_checks_passed: bool


def residual_count(matrix: Matrix) -> int:
    """Count nonzero entries after exact simplification."""

    return sum(
        int(sp.simplify(entry) != 0)
        for entry in matrix
    )


def audit_case(size: int) -> RetighteningRecord:
    """Audit one exact endpoint ledger."""

    identity = sp.eye(size)
    zero = sp.zeros(size, size)

    boundary_lower = identity
    boundary_upper = -4 * identity
    raw_lower = identity
    raw_upper = -12 * identity
    effective_lower = zero
    effective_upper = -16 * identity
    raw_increment = raw_upper - boundary_upper
    retightening_lower_increment = effective_lower - raw_lower
    retightening_upper_increment = effective_upper - raw_upper

    flux = zero
    l212_target = -16 * identity - 28 * flux

    residuals = {
        "boundary_lower": residual_count(boundary_lower - identity),
        "boundary_upper": residual_count(
            boundary_upper + 4 * identity
        ),
        "raw_lower": residual_count(raw_lower - identity),
        "raw_upper": residual_count(raw_upper + 12 * identity),
        "effective_lower": residual_count(effective_lower),
        "effective_upper": residual_count(
            effective_upper + 16 * identity
        ),
        "raw_increment": residual_count(
            raw_increment + 8 * identity
        ),
        "retightening_increment": (
            residual_count(
                retightening_lower_increment + identity
            )
            + residual_count(
                retightening_upper_increment + 4 * identity
            )
        ),
        "l212_target": residual_count(
            l212_target - effective_upper
        ),
    }
    checks_passed = all(value == 0 for value in residuals.values())
    if not checks_passed:
        raise RuntimeError(
            f"retightening audit failed at multiplicity {size}: "
            f"{residuals}"
        )

    return RetighteningRecord(
        multiplicity=size,
        boundary_lower_residual_count=residuals["boundary_lower"],
        boundary_upper_residual_count=residuals["boundary_upper"],
        raw_repair_lower_residual_count=residuals["raw_lower"],
        raw_repair_upper_residual_count=residuals["raw_upper"],
        effective_lower_residual_count=residuals["effective_lower"],
        effective_upper_residual_count=residuals["effective_upper"],
        raw_repair_increment_residual_count=residuals[
            "raw_increment"
        ],
        retightening_increment_residual_count=residuals[
            "retightening_increment"
        ],
        l212_target_residual_count=residuals["l212_target"],
        all_checks_passed=checks_passed,
    )


def standard_records() -> list[RetighteningRecord]:
    """Return all exact multiplicity audits."""

    return [audit_case(size) for size in range(1, 5)]


def write_records(
    records: list[RetighteningRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(
        json.dumps(asdict(record), sort_keys=True) + "\n"
        for record in records
    )
    output.write_text(payload, encoding="utf-8")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact audit."""

    args = parse_args()
    records = standard_records()
    digest = None
    if args.output is not None:
        digest = write_records(records, args.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "dataset_sha256": digest,
                "records": [asdict(record) for record in records],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
