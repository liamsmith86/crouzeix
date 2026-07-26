#!/usr/bin/env python3
"""Audit the copy-space energy quotient of the moving retightening defect.

L305 makes every bridge-word endpoint a rooted transfer term plus a
bounded response.  The rooted class of the complete L299 defect has a
coordinate-free copy-space representative:

    -Psi(Delta* X S + S* X Delta + Delta* X Delta)
    -(C* H + H* C + C* C).

Here ``Psi(Z)=V* G_S^vee(Z) V`` is the dual Stein closure.  The
checker constructs the actual canonical moving operator and defect
frame, verifies this quotient identity, and confirms that the
remaining endpoint lies in the perpendicular response range.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_canonical_quintic_preimage import (
    canonical_factor_series,
    canonical_metric_slack_and_operator,
)
from repeated_crabb_cubic_markov_flux import minimum_response_column
from repeated_crabb_elliptic_cokernel import (
    dual_stein_inverse,
    stein_inverse,
)
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)
from repeated_crabb_rooted_bridge_flux import transfer_channel
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray
DEGREE = 4
PARAMETER = 0.06
THETA = 0.5


@dataclass(frozen=True)
class CopyEnergyQuotientRecord:
    """One audit of the canonical copy-energy quotient."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    correction_grades: str
    parameter_scale: str
    active_transfer_norm: str
    canonical_factor_error: str
    fixed_base_residual_norm: str
    compact_defect_reconstruction_error: str
    copy_quotient_reconstruction_error: str
    endpoint_norm: str
    retained_quotient_norm: str
    response_synthesis_error: str
    response_column_norm: str
    response_column_over_active_transfer: str
    all_checks_passed: bool


def evaluate_series(series: list[Matrix], parameter: float) -> Matrix:
    """Evaluate one finite matrix series."""

    result = np.zeros_like(series[0], dtype=complex)
    power = 1.0
    for coefficient in series:
        result += power * coefficient
        power *= parameter
    return result


def copy_closure(
    forcing: Matrix,
    operator: Matrix,
    right: Matrix,
) -> Matrix:
    """Return ``Psi(forcing)=V* G_S^vee(forcing) V``."""

    return (
        right.conj().T
        @ dual_stein_inverse(operator, forcing)
        @ right
    )


def canonical_moving_pair(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, Matrix, float]:
    """Return the degree-four canonical operator and frame at ``c``."""

    _, slack, operator_series, _ = (
        canonical_metric_slack_and_operator(
            operator,
            right,
            left,
            DEGREE,
        )
    )
    frame_series, factor_error = canonical_factor_series(
        slack,
        right,
    )
    return (
        evaluate_series(operator_series, PARAMETER),
        evaluate_series(frame_series, PARAMETER),
        factor_error,
    )


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
    correction_grades: tuple[int, ...],
) -> CopyEnergyQuotientRecord:
    """Audit one physical moving-pair quotient."""

    moving_operator, moving_frame, factor_error = (
        canonical_moving_pair(operator, right, left)
    )
    directions = [
        oriented_retightening_direction(
            operator,
            right,
            left,
            grade,
        )
        for grade in correction_grades
    ]
    metric_correction = THETA * sum(
        (
            PARAMETER ** (2 * grade) * direction.metric
            for grade, direction in zip(
                correction_grades,
                directions,
                strict=True,
            )
        ),
        np.zeros_like(operator),
    )
    frame_correction = THETA * sum(
        (
            PARAMETER ** (2 * grade) * direction.frame
            for grade, direction in zip(
                correction_grades,
                directions,
                strict=True,
            )
        ),
        np.zeros_like(right),
    )

    fixed_base_residual = (
        metric_correction
        - operator.conj().T @ metric_correction @ operator
        - right @ frame_correction.conj().T
        - frame_correction @ right.conj().T
    )
    full_residual = (
        metric_correction
        - moving_operator.conj().T
        @ metric_correction
        @ moving_operator
        - moving_frame @ frame_correction.conj().T
        - frame_correction @ moving_frame.conj().T
        - frame_correction @ frame_correction.conj().T
    )
    moving_defect = full_residual - fixed_base_residual

    operator_motion = moving_operator - operator
    frame_motion = moving_frame - right
    compact_defect = -(
        operator_motion.conj().T
        @ metric_correction
        @ operator
        + operator.conj().T
        @ metric_correction
        @ operator_motion
        + operator_motion.conj().T
        @ metric_correction
        @ operator_motion
        + frame_motion @ frame_correction.conj().T
        + frame_correction @ frame_motion.conj().T
        + frame_correction @ frame_correction.conj().T
    )

    full_copy_quotient = (
        copy_closure(
            metric_correction
            - moving_operator.conj().T
            @ metric_correction
            @ moving_operator,
            operator,
            right,
        )
        - frame_correction.conj().T @ moving_frame
        - moving_frame.conj().T @ frame_correction
        - frame_correction.conj().T @ frame_correction
    )
    base_copy_quotient = (
        copy_closure(
            metric_correction
            - operator.conj().T
            @ metric_correction
            @ operator,
            operator,
            right,
        )
        - frame_correction.conj().T @ right
        - right.conj().T @ frame_correction
    )
    moving_copy_quotient = full_copy_quotient - base_copy_quotient
    compact_copy_quotient = -(
        copy_closure(
            operator_motion.conj().T
            @ metric_correction
            @ operator
            + operator.conj().T
            @ metric_correction
            @ operator_motion
            + operator_motion.conj().T
            @ metric_correction
            @ operator_motion,
            operator,
            right,
        )
        + frame_correction.conj().T @ frame_motion
        + frame_motion.conj().T @ frame_correction
        + frame_correction.conj().T @ frame_correction
    )

    endpoint = (
        left.conj().T
        @ stein_inverse(operator, moving_defect)
        @ left
    )
    retained_endpoint = transfer_channel(
        moving_copy_quotient,
        operator,
        right,
        left,
    )
    response_target = endpoint - retained_endpoint
    if np.linalg.norm(response_target) < 1e-12:
        response_column = np.zeros_like(frame_correction)
        response_error = float(np.linalg.norm(response_target))
    else:
        response_column, response_error = minimum_response_column(
            operator,
            right,
            left,
            response_target,
        )

    active_transfer_norm = sum(
        float(
            np.linalg.norm(
                left.conj().T
                @ np.linalg.matrix_power(
                    operator.conj().T,
                    grade,
                )
                @ right
            )
        )
        for grade in correction_grades
    )
    column_norm = float(np.linalg.norm(response_column))
    if active_transfer_norm > 1e-12:
        normalized_column = column_norm / active_transfer_norm
    else:
        normalized_column = 0.0

    fixed_norm = float(np.linalg.norm(fixed_base_residual))
    defect_error = float(
        np.linalg.norm(moving_defect - compact_defect)
    )
    quotient_error = float(
        np.linalg.norm(
            moving_copy_quotient - compact_copy_quotient
        )
    )
    endpoint_norm = float(np.linalg.norm(endpoint))
    retained_norm = float(np.linalg.norm(retained_endpoint))

    tolerance = 3e-8
    verified = bool(
        factor_error < tolerance
        and fixed_norm < tolerance
        and defect_error < tolerance
        and quotient_error < tolerance
        and response_error < tolerance
        and (
            active_transfer_norm > 1e-12
            or max(endpoint_norm, retained_norm, column_norm)
            < tolerance
        )
    )
    if not verified:
        raise RuntimeError(
            "the copy-energy quotient audit failed: "
            f"{construction_kind=}, {parameter_scale=}, "
            f"{factor_error=}, {fixed_norm=}, {defect_error=}, "
            f"{quotient_error=}, {response_error=}, "
            f"{active_transfer_norm=}, {endpoint_norm=}, "
            f"{retained_norm=}, {column_norm=}"
        )

    return CopyEnergyQuotientRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        correction_grades=",".join(
            str(grade) for grade in correction_grades
        ),
        parameter_scale=format_float(parameter_scale),
        active_transfer_norm=format_float(active_transfer_norm),
        canonical_factor_error=format_float(factor_error),
        fixed_base_residual_norm=format_float(fixed_norm),
        compact_defect_reconstruction_error=format_float(
            defect_error
        ),
        copy_quotient_reconstruction_error=format_float(
            quotient_error
        ),
        endpoint_norm=format_float(endpoint_norm),
        retained_quotient_norm=format_float(retained_norm),
        response_synthesis_error=format_float(response_error),
        response_column_norm=format_float(column_norm),
        response_column_over_active_transfer=format_float(
            normalized_column
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[CopyEnergyQuotientRecord]:
    """Return deterministic general, rank-chain, and delayed cases."""

    records = []
    for multiplicity in range(1, 5):
        operator, right, left = random_partial_isometry(
            3 * multiplicity + 5,
            multiplicity,
            np.random.default_rng(142_000 + multiplicity),
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
                142_100 + multiplicity,
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

    generator = np.random.default_rng(142_200)
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
    records: list[CopyEnergyQuotientRecord],
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
            "repeated_crabb_copy_energy_quotient_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the copy-energy quotient audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
