#!/usr/bin/env python3
"""Audit the exact one-delay block and terminal-crossing formulas.

Removing the first left wandering layer from a repeated-Crabb partial
isometry puts the balanced ellipse pencil into a two-by-two arrowhead
form.  Its retained resolvent is a tail resolvent with one explicit
rank-``m`` terminal-crossing insertion.  This checker verifies that
normal form, the shifted transfer coefficients, and the corresponding
boundary-metric blocks on deterministic delayed colligations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_slack_deflation import (
    slack_schur_coefficients,
)
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_one_image_generator import (
    boundary_metric_coefficient,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class OneDelayBlockRecord:
    """One numerical audit of the one-delay normal form."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    tail_dimension: int
    block_operator_error: str
    tail_partial_isometry_error: str
    tail_defect_orthogonality_error: str
    maximum_transfer_shift_error: str
    maximum_pencil_block_error: str
    maximum_resolvent_error: str
    maximum_boundary_block_error: str
    leading_residual_covariance_error: str
    degree_six_whole_series_gap: str
    all_checks_passed: bool


def one_delay_tail(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, Matrix, Matrix, Matrix]:
    """Remove the first left wandering layer."""

    retained = null_space(left.conj().T)
    tail = retained.conj().T @ partial @ retained
    tail_right = retained.conj().T @ right
    tail_left = retained.conj().T @ partial @ left
    return retained, tail, tail_right, tail_left


def ordered_basis(
    left: Matrix,
    retained: Matrix,
) -> Matrix:
    """Return the removed-plus-retained unitary basis."""

    return np.hstack((left, retained))


def direct_pencil(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter: float,
) -> Matrix:
    """Return the balanced inverse-ellipse pencil."""

    identity = np.eye(len(partial), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    return partial + parameter * (
        (identity + left_projection)
        @ partial.conj().T
        @ (identity + right_projection)
    )


def predicted_boundary_coefficient(
    tail: Matrix,
    tail_right: Matrix,
    tail_left: Matrix,
    left: Matrix,
    half_order: int,
) -> Matrix:
    """Return the one-delay block formula at ``q**half_order``."""

    multiplicity = left.shape[1]
    tail_dimension = len(tail)
    top_left = np.zeros((multiplicity, multiplicity), dtype=complex)
    top_right = np.zeros(
        (multiplicity, tail_dimension),
        dtype=complex,
    )
    bottom_right = (
        np.linalg.matrix_power(tail, half_order - 1)
        @ tail_left
        @ tail_left.conj().T
        @ np.linalg.matrix_power(
            tail.conj().T,
            half_order - 1,
        )
    )
    for divisor in range(1, half_order + 1):
        if half_order % divisor:
            continue
        sign = (-1) ** (half_order // divisor)
        transfer = transfer_coefficient(
            tail,
            tail_right,
            tail_left,
            divisor - 1,
        )
        right_orbit = (
            np.linalg.matrix_power(tail.conj().T, divisor)
            @ tail_right
        )
        top_left += sign * transfer @ transfer.conj().T
        top_right += sign * transfer @ right_orbit.conj().T
        bottom_right += sign * right_orbit @ right_orbit.conj().T
    return np.block(
        [
            [top_left, top_right],
            [top_right.conj().T, bottom_right],
        ]
    )


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> OneDelayBlockRecord:
    """Audit one complete first delay."""

    dimension, multiplicity = right.shape
    retained, tail, tail_right, tail_left = one_delay_tail(
        partial,
        right,
        left,
    )
    basis = ordered_basis(left, retained)
    zero = np.zeros((multiplicity, multiplicity), dtype=complex)
    zero_top_right = np.zeros(
        (multiplicity, len(tail)),
        dtype=complex,
    )
    block_operator = np.block(
        [
            [zero, zero_top_right],
            [tail_left, tail],
        ]
    )
    block_operator_error = float(
        np.linalg.norm(basis.conj().T @ partial @ basis - block_operator)
    )

    tail_right_projection = tail_right @ tail_right.conj().T
    tail_left_projection = tail_left @ tail_left.conj().T
    tail_identity = np.eye(len(tail), dtype=complex)
    tail_partial_isometry_error = max(
        float(
            np.linalg.norm(
                tail.conj().T @ tail
                - (tail_identity - tail_right_projection)
            )
        ),
        float(
            np.linalg.norm(
                tail @ tail.conj().T
                - (tail_identity - tail_left_projection)
            )
        ),
    )
    tail_defect_orthogonality_error = float(
        np.linalg.norm(tail_right.conj().T @ tail_left)
    )

    transfer_shift_error = 0.0
    for degree in range(1, 6):
        full_transfer = transfer_coefficient(
            partial,
            right,
            left,
            degree + 1,
        )
        tail_transfer = transfer_coefficient(
            tail,
            tail_right,
            tail_left,
            degree,
        )
        transfer_shift_error = max(
            transfer_shift_error,
            float(np.linalg.norm(full_transfer - tail_transfer)),
        )

    pencil_block_error = 0.0
    resolvent_error = 0.0
    for parameter in (0.03, 0.17, 0.41):
        full_pencil = direct_pencil(
            partial,
            right,
            left,
            parameter,
        )
        core_pencil = (
            tail
            + parameter
            * tail.conj().T
            @ (tail_identity + tail_right_projection)
        )
        pencil_block = np.block(
            [
                [
                    zero,
                    2 * parameter * tail_left.conj().T,
                ],
                [tail_left, core_pencil],
            ]
        )
        pencil_block_error = max(
            pencil_block_error,
            float(
                np.linalg.norm(
                    basis.conj().T @ full_pencil @ basis
                    - pencil_block
                )
            ),
        )

        for spectral_parameter in (2.1 + 0.2j, -2.4 + 0.3j):
            full_resolvent = np.linalg.inv(
                spectral_parameter
                * np.eye(dimension, dtype=complex)
                - full_pencil
            )
            retained_resolvent = (
                retained.conj().T @ full_resolvent @ retained
            )
            predicted_resolvent = np.linalg.inv(
                spectral_parameter * tail_identity
                - core_pencil
                - (2 * parameter / spectral_parameter)
                * tail_left_projection
            )
            resolvent_error = max(
                resolvent_error,
                float(
                    np.linalg.norm(
                        retained_resolvent - predicted_resolvent
                    )
                ),
            )

    boundary_block_error = 0.0
    for half_order in range(1, 6):
        coefficient = boundary_metric_coefficient(
            partial,
            right,
            left,
            2 * half_order,
        )
        predicted = predicted_boundary_coefficient(
            tail,
            tail_right,
            tail_left,
            left,
            half_order,
        )
        boundary_block_error = max(
            boundary_block_error,
            float(
                np.linalg.norm(
                    basis.conj().T @ coefficient @ basis
                    - predicted
                )
            ),
        )

    full_residual = slack_schur_coefficients(
        partial,
        right,
        left,
        6,
    )
    tail_residual = slack_schur_coefficients(
        tail,
        tail_right,
        tail_left,
        4,
    )
    leading_residual_covariance_error = float(
        np.linalg.norm(
            retained.conj().T @ full_residual[4] @ retained
            - tail_residual[2]
        )
    )
    degree_six_whole_series_gap = float(
        np.linalg.norm(
            retained.conj().T @ full_residual[6] @ retained
            - tail_residual[4]
        )
    )

    tolerance = 2e-11
    verified = max(
        block_operator_error,
        tail_partial_isometry_error,
        tail_defect_orthogonality_error,
        transfer_shift_error,
        pencil_block_error,
        resolvent_error,
        boundary_block_error,
        leading_residual_covariance_error,
    ) < tolerance
    verified = verified and degree_six_whole_series_gap > 1e-6
    if not verified:
        raise RuntimeError(
            "the one-delay block audit failed: "
            f"kind={construction_kind}, n={dimension}"
        )
    return OneDelayBlockRecord(
        construction_kind=construction_kind,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        tail_dimension=len(tail),
        block_operator_error=format_float(block_operator_error),
        tail_partial_isometry_error=format_float(
            tail_partial_isometry_error
        ),
        tail_defect_orthogonality_error=format_float(
            tail_defect_orthogonality_error
        ),
        maximum_transfer_shift_error=format_float(
            transfer_shift_error
        ),
        maximum_pencil_block_error=format_float(
            pencil_block_error
        ),
        maximum_resolvent_error=format_float(resolvent_error),
        maximum_boundary_block_error=format_float(
            boundary_block_error
        ),
        leading_residual_covariance_error=format_float(
            leading_residual_covariance_error
        ),
        degree_six_whole_series_gap=format_float(
            degree_six_whole_series_gap
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[OneDelayBlockRecord]:
    """Return deterministic delayed audits."""

    records: list[OneDelayBlockRecord] = []
    generator = np.random.default_rng(72324)
    for dimension, multiplicity in ((7, 2), (9, 2), (12, 3), (14, 4)):
        records.append(
            audit_case(
                "unstructured_complete_first_delay",
                *delayed_random_partial_isometry(
                    dimension,
                    multiplicity,
                    generator,
                ),
            )
        )
    for dimension, multiplicity in ((7, 2), (9, 3)):
        for seed in (72401, 72402):
            partial, right, left, _ = inflated_case(
                dimension,
                multiplicity,
                2,
                multiplicity,
                seed + dimension + multiplicity,
            )
            records.append(
                audit_case(
                    "inflated_complete_first_delay",
                    partial,
                    right,
                    left,
                )
            )
    return records


def write_records(
    records: list[OneDelayBlockRecord],
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
            "repeated_crabb_one_delay_block_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
