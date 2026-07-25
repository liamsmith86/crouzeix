#!/usr/bin/env python3
"""Audit the exact metric-edge part of the delayed dual Schur trace."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_one_image_generator import (
    boundary_metric_coefficient,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class DualMetricTraceSplitRecord:
    """One deterministic delayed metric-edge audit."""

    grade: int
    state_dimension: int
    defect_dimension: int
    maximum_earlier_transfer_norm: str
    lower_endpoint_error: str
    upper_endpoint_error: str
    metric_edge_trace_error: str
    all_checks_passed: bool


def audit_case(
    grade: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> DualMetricTraceSplitRecord:
    """Audit one completely delayed metric edge."""

    coefficient = boundary_metric_coefficient(
        partial,
        right,
        left,
        2 * grade,
    )
    active = transfer_coefficient(
        partial,
        right,
        left,
        grade,
    )
    frobenius_square = float(np.linalg.norm(active) ** 2)
    right_gram = active.conj().T @ active
    left_gram = active @ active.conj().T

    lower = right.conj().T @ coefficient @ right
    upper = left.conj().T @ coefficient @ left
    lower_error = float(np.linalg.norm(lower - right_gram))
    upper_error = float(np.linalg.norm(upper + left_gram))

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    left_projection = left @ left.conj().T
    complement = identity - left_projection
    metric_edge = (
        complement
        @ (
            -coefficient
            + partial @ coefficient @ partial.conj().T
        )
        @ complement
    )
    trace_error = abs(
        float(np.trace(metric_edge).real)
        + 2 * frobenius_square
    )
    earlier_transfer = max(
        (
            float(
                np.linalg.norm(
                    transfer_coefficient(
                        partial,
                        right,
                        left,
                        degree,
                    )
                )
            )
            for degree in range(1, grade)
        ),
        default=0.0,
    )

    tolerance = 8e-10
    verified = max(
        earlier_transfer,
        lower_error,
        upper_error,
        trace_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            "the dual metric trace-split audit failed: "
            f"grade={grade}, lower={lower_error:.3e}, "
            f"upper={upper_error:.3e}, trace={trace_error:.3e}"
        )
    return DualMetricTraceSplitRecord(
        grade=grade,
        state_dimension=dimension,
        defect_dimension=right.shape[1],
        maximum_earlier_transfer_norm=format_float(earlier_transfer),
        lower_endpoint_error=format_float(lower_error),
        upper_endpoint_error=format_float(upper_error),
        metric_edge_trace_error=format_float(float(trace_error)),
        all_checks_passed=verified,
    )


def standard_records(
    maximum_grade: int,
) -> list[DualMetricTraceSplitRecord]:
    """Return generic complete-delay audits."""

    records: list[DualMetricTraceSplitRecord] = []
    for grade in range(1, maximum_grade + 1):
        multiplicity = 2 if grade <= 5 else 1
        if grade == 1:
            partial, right, left = random_partial_isometry(
                8,
                multiplicity,
                np.random.default_rng(247_001),
            )
        else:
            partial, right, left, _ = inflated_case(
                grade + 3,
                multiplicity,
                grade,
                multiplicity,
                247_000 + grade,
            )
        records.append(audit_case(grade, partial, right, left))
    return records


def write_records(
    records: list[DualMetricTraceSplitRecord],
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
            "repeated_crabb_dual_metric_trace_split_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist all audits."""

    args = parse_args()
    records = standard_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
