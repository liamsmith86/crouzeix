#!/usr/bin/env python3
"""Audit the prepared sixth-order face on ``B_1=B_2=0``.

The grade-three gauge adds ``+2 Q S W B_3`` at frame order three and
its negative at order five.  The resulting raw upper face is the
desired channel base plus two explicit future-row endpoint terms.
One polynomial preparation cancels those terms, and L212 leaves
``-16 B_3 B_3*``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_delayed_jet import (
    construct_delayed_metric_jet,
    endpoint_coefficient,
)
from repeated_crabb_elliptic_cokernel import endpoint_motion
from repeated_crabb_grade_two_face import gauge_case
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    transfer_channel,
    transfer_coefficient,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift


@dataclass(frozen=True)
class GradeThreeFaceRecord:
    """One complete sixth-order delayed-face audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    spectral_radius: str
    maximum_earlier_transfer_norm: str
    third_transfer_norm: str
    maximum_lower_schur_error: str
    maximum_pre_sixth_upper_schur_error: str
    raw_endpoint_identity_error: str
    prepared_endpoint_identity_error: str
    final_endpoint_identity_error: str
    apex_total_correction_norm: str
    all_checks_passed: bool


def grade_three_adjustment(
    degree: int,
    partial: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
) -> np.ndarray:
    """Return the endpoint-null gauge needed before the sixth face."""

    if degree not in {3, 5}:
        return np.zeros_like(right)
    projection = np.eye(len(partial)) - right @ right.conj().T
    third = transfer_coefficient(partial, right, left, 3)
    sign = 1 if degree == 3 else -1
    return sign * 2 * projection @ partial @ left @ third


def audit_case(
    construction_kind: str,
    partial: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
) -> GradeThreeFaceRecord:
    """Audit one partial-isometry colligation with two delays."""

    degree = 6
    jet = construct_delayed_metric_jet(
        partial,
        right,
        left,
        degree,
        grade_three_adjustment,
    )
    coefficients = [
        transfer_coefficient(partial, right, left, index)
        for index in range(6)
    ]
    third = coefficients[3]
    left_gram = third @ third.conj().T
    channel = transfer_channel(
        partial,
        right,
        left,
        third.conj().T @ third,
    )
    projection = np.eye(len(partial)) - right @ right.conj().T
    future_row = np.zeros_like(right)
    for offset in (1, 2):
        future_row += (
            np.linalg.matrix_power(partial.conj().T, offset)
            @ right
            @ coefficients[3 + offset].conj().T
            @ third
        )
    future_row = projection @ future_row
    future_endpoint = endpoint_motion(
        jet.operator,
        right,
        left,
        jet.metric_root @ future_row,
    )

    raw_upper = endpoint_coefficient(jet, degree, upper=True)
    raw_expected = (
        12 * left_gram - 28 * channel + 2 * future_endpoint
    )
    raw_error = float(np.linalg.norm(raw_upper - raw_expected))

    preparation_column = -2 * future_row
    prepared_upper = raw_upper + endpoint_motion(
        jet.operator,
        right,
        left,
        jet.metric_root @ preparation_column,
    )
    prepared_expected = 12 * left_gram - 28 * channel
    prepared_error = float(
        np.linalg.norm(prepared_upper - prepared_expected)
    )

    l212_column = (
        -3.5
        * projection
        @ np.linalg.matrix_power(partial, 3)
        @ left
        @ third
    )
    final_upper = prepared_upper + endpoint_motion(
        jet.operator,
        right,
        left,
        jet.metric_root @ l212_column,
    )
    final_error = float(
        np.linalg.norm(final_upper + 16 * left_gram)
    )

    lower_errors = [
        float(
            np.linalg.norm(
                endpoint_coefficient(jet, index, upper=False)
            )
        )
        for index in range(1, degree + 1)
    ]
    earlier_upper_errors = [
        float(
            np.linalg.norm(
                endpoint_coefficient(jet, index, upper=True)
            )
        )
        for index in range(1, degree)
    ]
    earlier_transfer_norm = max(
        float(np.linalg.norm(coefficient))
        for coefficient in coefficients[1:3]
    )
    apex_correction_norm = 0.0
    if construction_kind == "repeated_grade_three_apex":
        apex_correction_norm = float(
            np.linalg.norm(preparation_column + l212_column)
        )

    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    tolerance = 4e-9
    verified = bool(
        spectral_radius < 1
        and earlier_transfer_norm < tolerance
        and max(lower_errors) < tolerance
        and max(earlier_upper_errors) < tolerance
        and raw_error < tolerance
        and prepared_error < tolerance
        and final_error < tolerance
        and apex_correction_norm < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the grade-three face audit failed: "
            f"earlier={earlier_transfer_norm:.3e}, "
            f"lower={max(lower_errors):.3e}, "
            f"upper<6={max(earlier_upper_errors):.3e}, "
            f"raw={raw_error:.3e}, prepared={prepared_error:.3e}, "
            f"final={final_error:.3e}, apex={apex_correction_norm:.3e}"
        )

    return GradeThreeFaceRecord(
        construction_kind=construction_kind,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        spectral_radius=format_float(spectral_radius),
        maximum_earlier_transfer_norm=format_float(
            earlier_transfer_norm
        ),
        third_transfer_norm=format_float(
            float(np.linalg.norm(third))
        ),
        maximum_lower_schur_error=format_float(max(lower_errors)),
        maximum_pre_sixth_upper_schur_error=format_float(
            max(earlier_upper_errors)
        ),
        raw_endpoint_identity_error=format_float(raw_error),
        prepared_endpoint_identity_error=format_float(prepared_error),
        final_endpoint_identity_error=format_float(final_error),
        apex_total_correction_norm=format_float(apex_correction_norm),
        all_checks_passed=verified,
    )


def standard_records() -> list[GradeThreeFaceRecord]:
    """Return deterministic unstructured, delayed, and apex audits."""

    records: list[GradeThreeFaceRecord] = []
    generator = np.random.default_rng(71524)
    for index, (dimension, multiplicity) in enumerate(
        ((8, 2), (10, 2), (12, 3), (15, 3), (16, 4))
    ):
        for repetition in range(2):
            records.append(
                audit_case(
                    "inflated_unstructured_partial_isometry",
                    *inflated_case(
                        dimension,
                        multiplicity,
                        3,
                        multiplicity,
                        72_000 + 10 * index + repetition,
                    )[:3],
                )
            )

    for lengths in (
        (3, 3),
        (3, 4, 5),
        (3, 5, 7),
        (3, 3, 6, 8),
        (4, 6, 8),
    ):
        records.append(
            audit_case(
                "gauged_heterogeneous_shift",
                *gauge_case(
                    *heterogeneous_shift(lengths),
                    generator,
                ),
            )
        )

    for multiplicity in (1, 2, 3):
        records.append(
            audit_case(
                "repeated_grade_three_apex",
                *gauge_case(
                    *heterogeneous_shift((3,) * multiplicity),
                    generator,
                ),
            )
        )
    return records


def write_records(
    records: list[GradeThreeFaceRecord],
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
            "experiments/repeated_crabb_grade_three_face_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete sixth-order audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
