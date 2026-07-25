#!/usr/bin/env python3
"""Audit the first delayed Schur defect of the boundary-metric slack.

Let ``R(c)`` be L219's balanced boundary-layer metric and let

    H(c) = R(c) - A(c)* R(c) A(c)

be its Stein slack for the balanced elliptic pullback.  Schur-compress
``H(c)`` away from the right defect.  On a fully delayed grade-``k``
colligation, the first nonzero coefficient is conjectured to be the
literal lift of the grade-one coefficient on the deflated colligation.

The trace of the grade-one coefficient is proved separately in the
companion note.  This checker tests the all-grade covariance and the
two exact boundary-metric endpoint faces.  It is a falsification audit,
not a proof of the covariance statement.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_delayed_jet import (
    ellipse_operator_coefficients,
    inverse_series,
    series_multiply,
)
from repeated_crabb_one_image_generator import boundary_metric_coefficient
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class BoundarySlackDeflationRecord:
    """One delayed boundary-slack covariance audit."""

    grade: int
    state_dimension: int
    defect_dimension: int
    maximum_earlier_transfer_norm: str
    maximum_earlier_slack_schur_norm: str
    deflated_face_error: str
    slack_face_trace_error: str
    lower_boundary_face_error: str
    upper_boundary_face_error: str
    reconstructed_upper_trace_error: str
    all_checks_passed: bool


def adjoint_series(series: list[Matrix]) -> list[Matrix]:
    """Return the coefficientwise adjoint series."""

    return [coefficient.conj().T for coefficient in series]


def physical_metric_root(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, Matrix]:
    """Return the repeated equality metric root and its inverse."""

    identity = np.eye(len(partial), dtype=complex)
    metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    root = np.asarray(sqrtm(metric), dtype=complex)
    return root, np.linalg.inv(root)


def slack_schur_coefficients(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    maximum_degree: int,
) -> list[Matrix]:
    """Return the right-defect Schur complement of the Stein slack."""

    root, inverse_root = physical_metric_root(partial, right, left)
    physical_operator = inverse_root @ partial @ root
    physical_series = ellipse_operator_coefficients(
        physical_operator,
        maximum_degree,
    )
    balanced_series = [
        root @ coefficient @ inverse_root
        for coefficient in physical_series
    ]
    metric_series = [np.eye(len(partial), dtype=complex)] + [
        boundary_metric_coefficient(
            partial,
            right,
            left,
            degree,
        )
        for degree in range(1, maximum_degree + 1)
    ]
    pulled_metric = series_multiply(
        series_multiply(
            adjoint_series(balanced_series),
            metric_series,
            maximum_degree,
        ),
        balanced_series,
        maximum_degree,
    )
    slack = [
        metric - pulled
        for metric, pulled in zip(
            metric_series,
            pulled_metric,
            strict=True,
        )
    ]

    identity = np.eye(len(partial), dtype=complex)
    right_projection = right @ right.conj().T
    complement = identity - right_projection
    pivot = [
        right.conj().T @ coefficient @ right
        for coefficient in slack
    ]
    cross = [
        complement @ coefficient @ right
        for coefficient in slack
    ]
    corner = [
        complement @ coefficient @ complement
        for coefficient in slack
    ]
    cross_square = series_multiply(
        series_multiply(
            cross,
            inverse_series(pivot),
            maximum_degree,
        ),
        adjoint_series(cross),
        maximum_degree,
    )
    return [
        corner_coefficient - square_coefficient
        for corner_coefficient, square_coefficient in zip(
            corner,
            cross_square,
            strict=True,
        )
    ]


def deflate_complete_delay(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    grade: int,
) -> tuple[Matrix, Matrix, Matrix, Matrix]:
    """Remove the first ``grade - 1`` complete left wandering layers."""

    if grade == 1:
        identity = np.eye(len(partial), dtype=complex)
        return identity, partial, right, left

    removed = np.hstack(
        [
            np.linalg.matrix_power(partial, degree) @ left
            for degree in range(grade - 1)
        ]
    )
    retained = null_space(removed.conj().T)
    deflated = retained.conj().T @ partial @ retained
    deflated_right = retained.conj().T @ right
    promoted_left = (
        retained.conj().T
        @ np.linalg.matrix_power(partial, grade - 1)
        @ left
    )
    return retained, deflated, deflated_right, promoted_left


def audit_case(
    grade: int,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> BoundarySlackDeflationRecord:
    """Audit one fully delayed first face."""

    face_degree = 2 * grade
    slack_schur = slack_schur_coefficients(
        partial,
        right,
        left,
        face_degree,
    )
    retained, deflated, deflated_right, deflated_left = (
        deflate_complete_delay(
            partial,
            right,
            left,
            grade,
        )
    )
    deflated_face = slack_schur_coefficients(
        deflated,
        deflated_right,
        deflated_left,
        2,
    )[2]
    lifted_face = retained @ deflated_face @ retained.conj().T
    face = slack_schur[face_degree]

    coefficients = [
        transfer_coefficient(
            partial,
            right,
            left,
            degree,
        )
        for degree in range(grade + 1)
    ]
    active = coefficients[grade]
    right_gram = active.conj().T @ active
    left_gram = active @ active.conj().T
    frobenius_square = float(np.linalg.norm(active) ** 2)

    boundary_face = boundary_metric_coefficient(
        partial,
        right,
        left,
        face_degree,
    )
    lower_face = right.conj().T @ boundary_face @ right
    upper_face = 4 * left.conj().T @ boundary_face @ left

    earlier_transfer = max(
        (
            float(np.linalg.norm(coefficient))
            for coefficient in coefficients[1:grade]
        ),
        default=0.0,
    )
    earlier_slack = max(
        (
            float(np.linalg.norm(coefficient))
            for coefficient in slack_schur[:face_degree]
        ),
        default=0.0,
    )
    covariance_error = float(np.linalg.norm(face - lifted_face))
    trace_error = abs(
        float(np.trace(face).real) - 2 * frobenius_square
    )
    lower_error = float(np.linalg.norm(lower_face - right_gram))
    upper_error = float(np.linalg.norm(upper_face + 4 * left_gram))

    # If a correction kills the slack Schur face and re-tightens the
    # lower endpoint, the dual trace identity gives
    #
    #   delta_upper - 4 delta_lower = -4 tr(face).
    #
    # Insert delta_lower = -||B_k||_F^2 and add the boundary upper
    # trace -4||B_k||_F^2.
    reconstructed_upper_trace = (
        float(np.trace(upper_face).real)
        + 4 * (-frobenius_square)
        - 4 * float(np.trace(face).real)
    )
    reconstructed_error = abs(
        reconstructed_upper_trace + 16 * frobenius_square
    )

    tolerance = 2e-7
    verified = bool(
        earlier_transfer < tolerance
        and earlier_slack < tolerance
        and covariance_error < tolerance
        and trace_error < tolerance
        and lower_error < tolerance
        and upper_error < tolerance
        and reconstructed_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "boundary-slack deflation audit failed: "
            f"grade={grade}, transfer={earlier_transfer:.3e}, "
            f"earlier={earlier_slack:.3e}, "
            f"covariance={covariance_error:.3e}, "
            f"trace={trace_error:.3e}, lower={lower_error:.3e}, "
            f"upper={upper_error:.3e}, "
            f"reconstructed={reconstructed_error:.3e}"
        )
    return BoundarySlackDeflationRecord(
        grade=grade,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        maximum_earlier_transfer_norm=format_float(earlier_transfer),
        maximum_earlier_slack_schur_norm=format_float(earlier_slack),
        deflated_face_error=format_float(covariance_error),
        slack_face_trace_error=format_float(trace_error),
        lower_boundary_face_error=format_float(lower_error),
        upper_boundary_face_error=format_float(upper_error),
        reconstructed_upper_trace_error=format_float(
            reconstructed_error
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[BoundarySlackDeflationRecord]:
    """Return deterministic unstructured delayed records."""

    records: list[BoundarySlackDeflationRecord] = []
    for grade in range(1, 6):
        if grade == 1:
            partial, right, left = random_partial_isometry(
                8,
                2,
                np.random.default_rng(107_001),
            )
        else:
            partial, right, left, _ = inflated_case(
                grade + 4,
                2,
                grade,
                2,
                107_000 + grade,
            )
        records.append(audit_case(grade, partial, right, left))
    return records


def write_records(
    records: list[BoundarySlackDeflationRecord],
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
            "repeated_crabb_boundary_slack_deflation_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
