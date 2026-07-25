#!/usr/bin/env python3
"""Audit the first delayed face of L219's boundary-layer metric.

If ``B_1 = ... = B_(k-1) = 0``, the upper Schur complement of the
physical boundary metric has first coefficient

    -4 B_k B_k*

at order ``c^(2k)``.  This checker constructs the metric coefficients
directly from the two defect-orbit sums; it does not use an elliptic
Riemann-map or Stein-factorization approximation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_delayed_jet import schur_cross_coefficient
from repeated_crabb_one_image_generator import (
    boundary_metric_coefficient,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


Matrix = np.ndarray


@dataclass(frozen=True)
class BoundaryMetricFaceRecord:
    """One all-grade upper-face audit."""

    construction_kind: str
    grade: int
    state_dimension: int
    defect_dimension: int
    maximum_earlier_transfer_norm: str
    active_transfer_norm: str
    maximum_earlier_upper_face_norm: str
    leading_face_error: str
    leading_face_maximum_eigenvalue: str
    all_checks_passed: bool


def metric_coefficients(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    degree: int,
) -> list[Matrix]:
    """Return physical L219 metric coefficients through ``degree``."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    base = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    root = np.asarray(sqrtm(base), dtype=complex)
    coefficients = [base]
    for order in range(1, degree + 1):
        balanced = boundary_metric_coefficient(
            partial,
            right,
            left,
            order,
        )
        coefficients.append(root @ balanced @ root)
    return coefficients


def upper_face(
    coefficients: list[Matrix],
    left: Matrix,
    order: int,
) -> Matrix:
    """Return one upper generalized-endpoint coefficient."""

    result = (
        left.conj().T @ coefficients[order] @ left
        - schur_cross_coefficient(
            coefficients,
            left,
            4,
            order,
        )
    )
    return (result + result.conj().T) / 2


def audit_case(
    construction_kind: str,
    grade: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> BoundaryMetricFaceRecord:
    """Audit one fully delayed boundary-metric face."""

    degree = 2 * grade
    coefficients = metric_coefficients(
        partial,
        right,
        left,
        degree,
    )
    transfer = [
        transfer_coefficient(partial, right, left, index)
        for index in range(grade + 1)
    ]
    earlier_transfer = max(
        (
            float(np.linalg.norm(transfer[index]))
            for index in range(1, grade)
        ),
        default=0.0,
    )
    active = transfer[grade]

    earlier_faces = [
        float(np.linalg.norm(upper_face(coefficients, left, order)))
        for order in range(1, degree)
    ]
    face = upper_face(coefficients, left, degree)
    expected = -4 * active @ active.conj().T
    face_error = float(np.linalg.norm(face - expected))
    face_maximum = float(np.linalg.eigvalsh(face)[-1])

    tolerance = 3e-8
    verified = bool(
        earlier_transfer < tolerance
        and max(earlier_faces, default=0.0) < tolerance
        and face_error < tolerance
        and face_maximum < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the boundary-metric face audit failed: "
            f"kind={construction_kind}, grade={grade}, "
            f"transfer={earlier_transfer:.3e}, "
            f"preface={max(earlier_faces, default=0.0):.3e}, "
            f"face={face_error:.3e}"
        )

    return BoundaryMetricFaceRecord(
        construction_kind=construction_kind,
        grade=grade,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        maximum_earlier_transfer_norm=format_float(earlier_transfer),
        active_transfer_norm=format_float(float(np.linalg.norm(active))),
        maximum_earlier_upper_face_norm=format_float(
            max(earlier_faces, default=0.0)
        ),
        leading_face_error=format_float(face_error),
        leading_face_maximum_eigenvalue=format_float(face_maximum),
        all_checks_passed=verified,
    )


def standard_records() -> list[BoundaryMetricFaceRecord]:
    """Return general, rank-changing, and apex audits."""

    records = []
    for grade in range(1, 8):
        if grade == 1:
            partial, right, left = random_partial_isometry(
                9,
                2,
                np.random.default_rng(107_101),
            )
        else:
            partial, right, left, _ = inflated_case(
                8,
                2,
                grade,
                2,
                107_100 + grade,
            )
        records.append(
            audit_case(
                "unstructured_fully_delayed",
                grade,
                partial,
                right,
                left,
            )
        )

        channel_lengths = (grade, grade + 1, grade + 3)
        records.append(
            audit_case(
                "rank_changing_shift",
                grade,
                *heterogeneous_shift(channel_lengths),
            )
        )
        records.append(
            audit_case(
                "repeated_monomial_apex",
                grade,
                *heterogeneous_shift((grade, grade)),
            )
        )
    return records


def write_records(
    records: list[BoundaryMetricFaceRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_boundary_metric_face_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
