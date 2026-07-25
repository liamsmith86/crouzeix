#!/usr/bin/env python3
"""Audit the candidate delayed trace law beyond the exact word range.

The proposed argument has an unresolved post-Schur endpoint-filtration
step.  This checker is an independent adversarial audit on generic
noncommuting tails.  It uses
the full theta/ODE direct-map series, boundary metric, Stein slack, and
right-defect Schur complement, then verifies

    tr [c**(2*k)] K = 2 ||B_k||_F**2

and the equivalent physical effective trace ``-16 ||B_k||_F**2``.
Grades one through five duplicate the exact word audit; grades six
through eight extend the floating falsification range.  Floating
agreement does not prove the missing filtration step.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from repeated_crabb_boundary_slack_deflation import (
    BoundarySlackDeflationRecord,
    audit_case,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
)


@dataclass(frozen=True)
class DelayedTraceLawRecord:
    """One deterministic full-pipeline trace audit."""

    grade: int
    state_dimension: int
    defect_dimension: int
    maximum_earlier_transfer_norm: str
    maximum_earlier_slack_schur_norm: str
    slack_face_trace_error: str
    reconstructed_upper_trace_error: str
    full_matrix_covariance_error: str
    all_checks_passed: bool


def trace_record(
    record: BoundarySlackDeflationRecord,
) -> DelayedTraceLawRecord:
    """Retain the trace-law fields from one full covariance audit."""

    return DelayedTraceLawRecord(
        grade=record.grade,
        state_dimension=record.state_dimension,
        defect_dimension=record.defect_dimension,
        maximum_earlier_transfer_norm=(
            record.maximum_earlier_transfer_norm
        ),
        maximum_earlier_slack_schur_norm=(
            record.maximum_earlier_slack_schur_norm
        ),
        slack_face_trace_error=record.slack_face_trace_error,
        reconstructed_upper_trace_error=(
            record.reconstructed_upper_trace_error
        ),
        full_matrix_covariance_error=record.deflated_face_error,
        all_checks_passed=record.all_checks_passed,
    )


def standard_records(
    maximum_grade: int,
) -> list[DelayedTraceLawRecord]:
    """Return deterministic generic-tail audits."""

    records: list[DelayedTraceLawRecord] = []
    for grade in range(1, maximum_grade + 1):
        multiplicity = 2 if grade <= 5 else 1
        if grade == 1:
            partial, right, left = random_partial_isometry(
                8,
                multiplicity,
                np.random.default_rng(246_001),
            )
        else:
            partial, right, left, _ = inflated_case(
                grade + 3,
                multiplicity,
                grade,
                multiplicity,
                246_000 + grade,
            )
        records.append(
            trace_record(
                audit_case(grade, partial, right, left)
            )
        )
    return records


def write_records(
    records: list[DelayedTraceLawRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

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
    parser.add_argument("--maximum-grade", type=int, default=8)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_delayed_trace_law_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the delayed trace-law audits."""

    args = parse_args()
    records = standard_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
