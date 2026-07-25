#!/usr/bin/env python3
"""Audit the two Hardy analysis frames behind the delayed slack.

The right- and left-defect orbit resolutions define two isometric
Hardy analysis maps.  Their cross Gram is the block Hankel matrix of
the transfer coefficients.  The boundary metric is the pullback of
fixed scalar diagonal weights in these two frames, and L228's proposed
leading anticommutator is the Hermitian lift of the first nonzero
Hankel cell.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_elliptic_cokernel import (
    dual_stein_inverse,
    stein_inverse,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class HardyTwoFrameRecord:
    """One audit of the two-frame/Hankel realization."""

    construction_kind: str
    grade: int
    state_dimension: int
    defect_dimension: int
    right_isometry_error: str
    left_isometry_error: str
    maximum_shift_intertwining_error: str
    maximum_hankel_cell_error: str
    maximum_earlier_hankel_cell_norm: str
    first_active_hankel_cell_error: str
    boundary_metric_pullback_error: str
    target_cell_lift_error: str
    all_checks_passed: bool


def analysis_rows(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    row_count: int,
) -> tuple[list[Matrix], list[Matrix]]:
    """Return finite prefixes of the two Hardy analysis maps."""

    right_rows = [
        right.conj().T @ np.linalg.matrix_power(partial, degree)
        for degree in range(row_count)
    ]
    left_rows = [
        left.conj().T
        @ np.linalg.matrix_power(partial.conj().T, degree)
        for degree in range(row_count)
    ]
    return right_rows, left_rows


def boundary_metric_from_frames(
    right_rows: list[Matrix],
    left_rows: list[Matrix],
    parameter: float,
) -> Matrix:
    """Return a finite weighted pullback of the boundary metric."""

    dimension = right_rows[0].shape[1]
    metric = np.eye(dimension, dtype=complex)
    q = parameter**2
    for degree in range(1, len(right_rows)):
        right_weight = -(q**degree) / (1 + q**degree)
        left_weight = q**degree
        metric += (
            right_weight
            * right_rows[degree].conj().T
            @ right_rows[degree]
            + left_weight
            * left_rows[degree].conj().T
            @ left_rows[degree]
        )
    return metric


def boundary_metric_direct(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter: float,
    row_count: int,
) -> Matrix:
    """Return the same finite orbit sum directly in state space."""

    identity = np.eye(len(partial), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    q = parameter**2
    metric = identity.copy()
    for degree in range(1, row_count):
        metric -= (
            (q**degree)
            / (1 + q**degree)
            * np.linalg.matrix_power(partial.conj().T, degree)
            @ right_projection
            @ np.linalg.matrix_power(partial, degree)
        )
        metric += (
            q**degree
            * np.linalg.matrix_power(partial, degree)
            @ left_projection
            @ np.linalg.matrix_power(partial.conj().T, degree)
        )
    return metric


def audit_case(
    construction_kind: str,
    grade: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> HardyTwoFrameRecord:
    """Audit one two-frame realization and active Hankel cell."""

    dimension, multiplicity = right.shape
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    right_gramian = stein_inverse(partial, right_projection)
    left_gramian = dual_stein_inverse(partial, left_projection)
    identity = np.eye(dimension, dtype=complex)
    right_isometry_error = float(np.linalg.norm(right_gramian - identity))
    left_isometry_error = float(np.linalg.norm(left_gramian - identity))

    row_count = max(8, grade + 4)
    right_rows, left_rows = analysis_rows(
        partial,
        right,
        left,
        row_count,
    )
    shift_error = 0.0
    for degree in range(row_count - 1):
        shift_error = max(
            shift_error,
            float(
                np.linalg.norm(
                    right_rows[degree] @ partial
                    - right_rows[degree + 1]
                )
            ),
            float(
                np.linalg.norm(
                    left_rows[degree] @ partial.conj().T
                    - left_rows[degree + 1]
                )
            ),
        )

    hankel_error = 0.0
    earlier_hankel = 0.0
    for right_degree in range(row_count // 2):
        for left_degree in range(row_count // 2):
            transfer_degree = right_degree + left_degree
            cell = (
                right_rows[right_degree]
                @ left_rows[left_degree].conj().T
            )
            transfer = transfer_coefficient(
                partial,
                right,
                left,
                transfer_degree,
            )
            hankel_error = max(
                hankel_error,
                float(np.linalg.norm(cell - transfer.conj().T)),
            )
            if 1 <= transfer_degree < grade:
                earlier_hankel = max(
                    earlier_hankel,
                    float(np.linalg.norm(cell)),
                )

    active_cell = right_rows[1] @ left_rows[grade - 1].conj().T
    active_transfer = transfer_coefficient(
        partial,
        right,
        left,
        grade,
    )
    active_cell_error = float(
        np.linalg.norm(active_cell - active_transfer.conj().T)
    )

    parameter = 0.31
    frame_metric = boundary_metric_from_frames(
        right_rows,
        left_rows,
        parameter,
    )
    direct_metric = boundary_metric_direct(
        partial,
        right,
        left,
        parameter,
        row_count,
    )
    metric_error = float(np.linalg.norm(frame_metric - direct_metric))

    first_right_orbit = (
        right_rows[1].conj().T @ right_rows[1]
    )
    active_left_orbit = (
        left_rows[grade - 1].conj().T
        @ left_rows[grade - 1]
    )
    target = (
        first_right_orbit @ active_left_orbit
        + active_left_orbit @ first_right_orbit
    )
    cell_lift = (
        right_rows[1].conj().T
        @ active_transfer.conj().T
        @ left_rows[grade - 1]
    )
    cell_lift += cell_lift.conj().T
    target_cell_lift_error = float(np.linalg.norm(target - cell_lift))

    tolerance = 2e-11
    verified = max(
        right_isometry_error,
        left_isometry_error,
        shift_error,
        hankel_error,
        earlier_hankel,
        active_cell_error,
        metric_error,
        target_cell_lift_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            "the two-frame audit failed: "
            f"kind={construction_kind}, grade={grade}"
        )
    return HardyTwoFrameRecord(
        construction_kind=construction_kind,
        grade=grade,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        right_isometry_error=format_float(right_isometry_error),
        left_isometry_error=format_float(left_isometry_error),
        maximum_shift_intertwining_error=format_float(shift_error),
        maximum_hankel_cell_error=format_float(hankel_error),
        maximum_earlier_hankel_cell_norm=format_float(earlier_hankel),
        first_active_hankel_cell_error=format_float(active_cell_error),
        boundary_metric_pullback_error=format_float(metric_error),
        target_cell_lift_error=format_float(target_cell_lift_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[HardyTwoFrameRecord]:
    """Return deterministic general and delayed audits."""

    records: list[HardyTwoFrameRecord] = []
    generator = np.random.default_rng(72524)
    for dimension, multiplicity in ((7, 2), (9, 2), (12, 3)):
        records.append(
            audit_case(
                "unstructured_grade_one",
                1,
                *random_partial_isometry(
                    dimension,
                    multiplicity,
                    generator,
                ),
            )
        )
    for grade in (2, 3, 4, 5):
        for multiplicity in (2, 3):
            partial, right, left, _ = inflated_case(
                7 + multiplicity,
                multiplicity,
                grade,
                multiplicity,
                72540 + 10 * grade + multiplicity,
            )
            records.append(
                audit_case(
                    "inflated_complete_delay",
                    grade,
                    partial,
                    right,
                    left,
                )
            )
    return records


def write_records(
    records: list[HardyTwoFrameRecord],
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
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_hardy_two_frame_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
