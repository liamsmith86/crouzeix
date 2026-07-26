#!/usr/bin/env python3
"""Audit the affine copy recurrence and its partial-scale Gram reserve."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_copy_energy_quotient import (
    PARAMETER,
    canonical_moving_pair,
    copy_closure,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray
THETA_VALUES = (0.2, 0.5, 0.8)


@dataclass(frozen=True)
class AffineCopyRecurrenceRecord:
    """One audit of recursive cancellation and partial scaling."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    correction_grades: str
    parameter_scale: str
    preceding_residual_norm: str
    recursive_column_norm: str
    affine_state_reconstruction_error: str
    affine_copy_reconstruction_error: str
    maximum_state_scaling_error: str
    maximum_copy_scaling_error: str
    minimum_state_reserve_eigenvalue: str
    minimum_copy_reserve_eigenvalue: str
    all_checks_passed: bool


def weighted_direction(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    correction_grades: tuple[int, ...],
) -> tuple[Matrix, Matrix]:
    """Return one simultaneous L298 metric/frame correction."""

    metric = np.zeros_like(operator)
    frame = np.zeros_like(right)
    for grade in correction_grades:
        direction = oriented_retightening_direction(
            operator,
            right,
            left,
            grade,
        )
        weight = PARAMETER ** (2 * grade)
        metric += weight * direction.metric
        frame += weight * direction.frame
    return metric, frame


def successor(
    operator_motion: Matrix,
    frame_motion: Matrix,
    operator: Matrix,
    metric: Matrix,
    frame: Matrix,
) -> Matrix:
    """Return the six-term moving successor."""

    return -(
        operator_motion.conj().T @ metric @ operator
        + operator.conj().T @ metric @ operator_motion
        + operator_motion.conj().T
        @ metric
        @ operator_motion
        + frame_motion @ frame.conj().T
        + frame @ frame_motion.conj().T
        + frame @ frame.conj().T
    )


def copy_successor(
    operator_motion: Matrix,
    frame_motion: Matrix,
    operator: Matrix,
    right: Matrix,
    metric: Matrix,
    frame: Matrix,
) -> Matrix:
    """Return the initial-copy successor."""

    operator_part = (
        operator_motion.conj().T @ metric @ operator
        + operator.conj().T @ metric @ operator_motion
        + operator_motion.conj().T
        @ metric
        @ operator_motion
    )
    return -(
        copy_closure(operator_part, operator, right)
        + frame.conj().T @ frame_motion
        + frame_motion.conj().T @ frame
        + frame.conj().T @ frame
    )


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
    correction_grades: tuple[int, ...],
) -> AffineCopyRecurrenceRecord:
    """Audit one genuinely affine recursive correction."""

    moving_operator, moving_frame, factor_error = (
        canonical_moving_pair(operator, right, left)
    )
    operator_motion = moving_operator - operator
    frame_motion = moving_frame - right

    first_metric, first_frame = weighted_direction(
        operator,
        right,
        left,
        correction_grades,
    )
    preceding_residual = successor(
        operator_motion,
        frame_motion,
        operator,
        first_metric,
        first_frame,
    )

    projection = right @ right.conj().T
    recursive_frame = (
        0.3
        * (np.eye(len(operator)) - projection)
        @ (first_frame + operator.conj().T @ first_frame)
    )
    affine_forcing = -preceding_residual
    recursive_metric = stein_inverse(
        operator,
        (
            affine_forcing
            + right @ recursive_frame.conj().T
            + recursive_frame @ right.conj().T
        ),
    )

    moving_correction_residual = (
        recursive_metric
        - moving_operator.conj().T
        @ recursive_metric
        @ moving_operator
        - moving_frame @ recursive_frame.conj().T
        - recursive_frame @ moving_frame.conj().T
        - recursive_frame @ recursive_frame.conj().T
    )
    new_successor = successor(
        operator_motion,
        frame_motion,
        operator,
        recursive_metric,
        recursive_frame,
    )
    state_error = float(
        np.linalg.norm(
            moving_correction_residual
            - affine_forcing
            - new_successor
        )
    )

    full_copy_value = (
        copy_closure(
            recursive_metric
            - moving_operator.conj().T
            @ recursive_metric
            @ moving_operator,
            operator,
            right,
        )
        - recursive_frame.conj().T @ moving_frame
        - moving_frame.conj().T @ recursive_frame
        - recursive_frame.conj().T @ recursive_frame
    )
    forcing_copy_value = copy_closure(
        affine_forcing,
        operator,
        right,
    )
    new_copy_successor = copy_successor(
        operator_motion,
        frame_motion,
        operator,
        right,
        recursive_metric,
        recursive_frame,
    )
    copy_error = float(
        np.linalg.norm(
            full_copy_value
            - forcing_copy_value
            - new_copy_successor
        )
    )

    maximum_state_scaling_error = 0.0
    maximum_copy_scaling_error = 0.0
    minimum_state_reserve = np.inf
    minimum_copy_reserve = np.inf
    state_gram = recursive_frame @ recursive_frame.conj().T
    copy_gram = recursive_frame.conj().T @ recursive_frame
    for theta in THETA_VALUES:
        scaled_state = successor(
            operator_motion,
            frame_motion,
            operator,
            theta * recursive_metric,
            theta * recursive_frame,
        )
        expected_state = (
            theta * new_successor
            + theta * (1 - theta) * state_gram
        )
        maximum_state_scaling_error = max(
            maximum_state_scaling_error,
            float(np.linalg.norm(scaled_state - expected_state)),
        )

        scaled_copy = copy_successor(
            operator_motion,
            frame_motion,
            operator,
            right,
            theta * recursive_metric,
            theta * recursive_frame,
        )
        expected_copy = (
            theta * new_copy_successor
            + theta * (1 - theta) * copy_gram
        )
        maximum_copy_scaling_error = max(
            maximum_copy_scaling_error,
            float(np.linalg.norm(scaled_copy - expected_copy)),
        )
        minimum_state_reserve = min(
            minimum_state_reserve,
            float(
                np.linalg.eigvalsh(
                    theta * (1 - theta) * state_gram
                )[0]
            ),
        )
        minimum_copy_reserve = min(
            minimum_copy_reserve,
            float(
                np.linalg.eigvalsh(
                    theta * (1 - theta) * copy_gram
                )[0]
            ),
        )

    tolerance = 5e-8
    verified = bool(
        factor_error < tolerance
        and state_error < tolerance
        and copy_error < tolerance
        and maximum_state_scaling_error < tolerance
        and maximum_copy_scaling_error < tolerance
        and minimum_state_reserve > -tolerance
        and minimum_copy_reserve > -tolerance
    )
    if not verified:
        raise RuntimeError(
            "affine copy recurrence audit failed: "
            f"{construction_kind=}, {factor_error=}, "
            f"{state_error=}, {copy_error=}, "
            f"{maximum_state_scaling_error=}, "
            f"{maximum_copy_scaling_error=}, "
            f"{minimum_state_reserve=}, "
            f"{minimum_copy_reserve=}"
        )

    return AffineCopyRecurrenceRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        correction_grades=",".join(
            str(grade) for grade in correction_grades
        ),
        parameter_scale=format_float(parameter_scale),
        preceding_residual_norm=format_float(
            float(np.linalg.norm(preceding_residual))
        ),
        recursive_column_norm=format_float(
            float(np.linalg.norm(recursive_frame))
        ),
        affine_state_reconstruction_error=format_float(state_error),
        affine_copy_reconstruction_error=format_float(copy_error),
        maximum_state_scaling_error=format_float(
            maximum_state_scaling_error
        ),
        maximum_copy_scaling_error=format_float(
            maximum_copy_scaling_error
        ),
        minimum_state_reserve_eigenvalue=format_float(
            minimum_state_reserve
        ),
        minimum_copy_reserve_eigenvalue=format_float(
            minimum_copy_reserve
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[AffineCopyRecurrenceRecord]:
    """Return deterministic general, rank-chain, and delayed cases."""

    records = []
    for multiplicity in range(1, 5):
        operator, right, left = random_partial_isometry(
            3 * multiplicity + 5,
            multiplicity,
            np.random.default_rng(143_000 + multiplicity),
        )
        records.append(
            audit_case(
                "unstructured",
                operator,
                right,
                left,
                1,
                (1, 2, 3),
            )
        )

    for multiplicity in (2, 3, 4):
        for parameter_scale in (1.0, 0.1, 0.01):
            operator, right, left, _, _ = rank_chain_case(
                multiplicity,
                143_100 + multiplicity,
                parameter_scale,
            )
            records.append(
                audit_case(
                    "rank_chain",
                    operator,
                    right,
                    left,
                    parameter_scale,
                    (1, 2),
                )
            )

    generator = np.random.default_rng(143_200)
    for multiplicity in (2, 3, 4):
        operator, right, left = delayed_random_partial_isometry(
            4 * multiplicity,
            multiplicity,
            generator,
        )
        records.append(
            audit_case(
                "complete_first_delay",
                operator,
                right,
                left,
                1,
                (2, 3),
            )
        )
    return records


def write_records(
    records: list[AffineCopyRecurrenceRecord],
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
            "repeated_crabb_affine_copy_recurrence_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the affine recurrence audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
