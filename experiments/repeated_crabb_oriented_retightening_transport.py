#!/usr/bin/env python3
"""Audit the all-grade oriented lower-retightening transport.

For one transfer coefficient ``B_k``, combine the negative
right-orbit Gram with one seventh of L212's endpoint column.  The
resulting metric/frame direction has physical endpoint pair

    (-B_k* B_k, -4 B_k B_k*).

Adding it to L283's completely delayed raw repair therefore changes
``(+B_k*B_k, -12B_kB_k*)`` to ``(0, -16B_kB_k*)``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import (
    all_grade_column,
    gauged_shift,
    stable_random_partial_isometry,
)
from repeated_crabb_boundary_slack_deflation import (
    physical_metric_root,
)
from repeated_crabb_edge_deleted_balanced_volume import delayed_case
from repeated_crabb_elliptic_cokernel import (
    endpoint_motion,
    stein_inverse,
)
from repeated_crabb_transfer_flag import (
    transfer_channel,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class RetighteningTransportRecord:
    """One deterministic transport audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    grade: int
    spectral_radius: str
    maximum_earlier_transfer_norm: str
    contamination_norm: str
    stein_frame_residual: str
    lower_endpoint_error: str
    upper_endpoint_error: str
    orbit_channel_error: str
    l212_component_error: str
    raw_to_effective_error: str
    apex_scalar_rescaling_error: str
    all_checks_passed: bool


@dataclass(frozen=True)
class RetighteningDirection:
    """The component matrices of one oriented transport."""

    metric: Matrix
    frame: Matrix
    right_gram: Matrix
    left_gram: Matrix
    orbit_gram: Matrix
    perpendicular_column: Matrix
    perpendicular_metric: Matrix
    contamination_norm: float


def validate_exact_rational_case() -> None:
    """Check the full identity on a noncommuting rational delay sum."""

    lengths = (2, 3)
    multiplicity = len(lengths)
    dimension = sum(length + 1 for length in lengths)
    operator = sp.zeros(dimension)
    right_base = sp.zeros(dimension, multiplicity)
    left_base = sp.zeros(dimension, multiplicity)
    offset = 0
    for channel, length in enumerate(lengths):
        right_base[offset, channel] = 1
        left_base[offset + length, channel] = 1
        for level in range(1, length + 1):
            operator[offset + level - 1, offset + level] = 1
        offset += length + 1

    right_gauge = sp.Matrix(
        [
            [sp.Rational(3, 5), sp.Rational(-4, 5)],
            [sp.Rational(4, 5), sp.Rational(3, 5)],
        ]
    )
    left_gauge = sp.Matrix(
        [
            [sp.Rational(5, 13), sp.Rational(-12, 13)],
            [sp.Rational(12, 13), sp.Rational(5, 13)],
        ]
    )
    right = right_base * right_gauge
    left = left_base * left_gauge
    projection = sp.eye(dimension) - right * right.T

    def coefficient(grade: int) -> sp.Matrix:
        return left.T * (operator.T**grade) * right

    def stein_sum(forcing: sp.Matrix) -> sp.Matrix:
        result = sp.zeros(dimension)
        power = sp.eye(dimension)
        for _ in range(dimension + 1):
            result += power.T * forcing * power
            power *= operator
        return sp.simplify(result)

    def channel(copy_matrix: sp.Matrix) -> sp.Matrix:
        return sp.simplify(
            left.T
            * stein_sum(right * copy_matrix * right.T)
            * left
        )

    for grade in range(1, 5):
        active = coefficient(grade)
        right_gram = active.T * active
        left_gram = active * active.T
        core = operator**grade * left * active
        for earlier_grade in range(1, grade):
            core += (
                operator.T ** (grade - earlier_grade)
                * right
                * active.T
                * coefficient(earlier_grade)
            )
        perpendicular = -sp.Rational(1, 2) * projection * core
        orbit = stein_sum(
            right * right_gram * right.T
        )
        perpendicular_metric = stein_sum(
            right * perpendicular.T
            + perpendicular * right.T
        )
        metric_direction = -orbit + perpendicular_metric
        frame_direction = (
            -sp.Rational(1, 2) * right * right_gram
            + perpendicular
        )
        residuals = (
            metric_direction
            - operator.T * metric_direction * operator
            - right * frame_direction.T
            - frame_direction * right.T,
            right.T * metric_direction * right + right_gram,
            left.T * metric_direction * left + left_gram,
            left.T * orbit * left - channel(right_gram),
            (
                left.T * perpendicular_metric * left
                - channel(right_gram)
                + left_gram
            ),
        )
        if not all(residual.is_zero_matrix for residual in residuals):
            raise RuntimeError(
                "the exact rational retightening guard failed: "
                f"{grade=}"
            )


def oriented_retightening_direction(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    grade: int,
) -> RetighteningDirection:
    """Return the balanced metric/frame direction and its two Grams."""

    l212_column, coefficient, contamination = all_grade_column(
        operator,
        right,
        left,
        grade,
    )
    right_gram = coefficient.conj().T @ coefficient
    left_gram = coefficient @ coefficient.conj().T
    orbit_gram = stein_inverse(
        operator,
        right @ right_gram @ right.conj().T,
    )
    perpendicular_column = l212_column / 7
    perpendicular_metric = stein_inverse(
        operator,
        (
            right @ perpendicular_column.conj().T
            + perpendicular_column @ right.conj().T
        ),
    )
    metric_direction = -orbit_gram + perpendicular_metric
    frame_direction = (
        -right @ right_gram / 2
        + perpendicular_column
    )
    return RetighteningDirection(
        metric=metric_direction,
        frame=frame_direction,
        right_gram=right_gram,
        left_gram=left_gram,
        orbit_gram=orbit_gram,
        perpendicular_column=perpendicular_column,
        perpendicular_metric=perpendicular_metric,
        contamination_norm=contamination,
    )


def maximum_earlier_transfer_norm(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    grade: int,
) -> float:
    """Return the largest transfer norm before ``grade``."""

    return max(
        (
            float(
                np.linalg.norm(
                    transfer_coefficient(
                        operator,
                        right,
                        left,
                        earlier,
                    )
                )
            )
            for earlier in range(1, grade)
        ),
        default=0.0,
    )


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    grade: int,
    *,
    require_complete_delay: bool = False,
    require_apex: bool = False,
) -> RetighteningTransportRecord:
    """Audit one all-grade transport identity."""

    direction = oriented_retightening_direction(
        operator,
        right,
        left,
        grade,
    )
    root, inverse_root = physical_metric_root(
        operator,
        right,
        left,
    )
    physical_operator = inverse_root @ operator @ root
    physical_metric = root @ direction.metric @ root

    stein_forcing = (
        right @ direction.frame.conj().T
        + direction.frame @ right.conj().T
    )
    stein_residual = float(
        np.linalg.norm(
            direction.metric
            - operator.conj().T @ direction.metric @ operator
            - stein_forcing
        )
    )
    lower_error = float(
        np.linalg.norm(
            right.conj().T @ physical_metric @ right
            + direction.right_gram
        )
    )
    upper_error = float(
        np.linalg.norm(
            left.conj().T @ physical_metric @ left
            + 4 * direction.left_gram
        )
    )

    channel_value = transfer_channel(
        operator,
        right,
        left,
        direction.right_gram,
    )
    physical_orbit = root @ (-direction.orbit_gram) @ root
    orbit_error = float(
        np.linalg.norm(
            left.conj().T @ physical_orbit @ left
            + 4 * channel_value
        )
    )

    l212_component = endpoint_motion(
        physical_operator,
        right,
        left,
        root @ direction.perpendicular_column,
    )
    l212_error = float(
        np.linalg.norm(
            l212_component
            - 4 * (channel_value - direction.left_gram)
        )
    )

    correction_lower = (
        right.conj().T @ physical_metric @ right
    )
    correction_upper = left.conj().T @ physical_metric @ left
    raw_to_effective_error = float(
        np.linalg.norm(direction.right_gram + correction_lower)
        + np.linalg.norm(
            -12 * direction.left_gram
            + correction_upper
            + 16 * direction.left_gram
        )
    )
    apex_error = (
        float(
            np.linalg.norm(
                direction.metric + np.eye(len(operator))
            )
        )
        if require_apex
        else 0.0
    )
    earlier_norm = maximum_earlier_transfer_norm(
        operator,
        right,
        left,
        grade,
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(operator)))
    )
    tolerance = 3e-9
    verified = bool(
        spectral_radius < 1
        and stein_residual < tolerance
        and lower_error < tolerance
        and upper_error < tolerance
        and orbit_error < tolerance
        and l212_error < tolerance
        and raw_to_effective_error < tolerance
        and (not require_complete_delay or earlier_norm < tolerance)
        and (not require_apex or apex_error < tolerance)
    )
    if not verified:
        raise RuntimeError(
            "oriented retightening transport failed: "
            f"{construction_kind=}, {grade=}"
        )

    return RetighteningTransportRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        grade=grade,
        spectral_radius=format_float(spectral_radius),
        maximum_earlier_transfer_norm=format_float(earlier_norm),
        contamination_norm=format_float(direction.contamination_norm),
        stein_frame_residual=format_float(stein_residual),
        lower_endpoint_error=format_float(lower_error),
        upper_endpoint_error=format_float(upper_error),
        orbit_channel_error=format_float(orbit_error),
        l212_component_error=format_float(l212_error),
        raw_to_effective_error=format_float(raw_to_effective_error),
        apex_scalar_rescaling_error=format_float(apex_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[RetighteningTransportRecord]:
    """Return unstructured, delayed, and apex audits."""

    validate_exact_rational_case()
    records: list[RetighteningTransportRecord] = []
    generator = np.random.default_rng(729_801)
    for dimension, multiplicity in ((7, 2), (9, 3), (12, 3)):
        data = stable_random_partial_isometry(
            dimension,
            multiplicity,
            generator,
        )
        for grade in range(1, 7):
            records.append(
                audit_case(
                    "unstructured_partial_isometry",
                    *data,
                    grade,
                )
            )

    for grade in range(1, 7):
        records.append(
            audit_case(
                "completely_delayed_colligation",
                *delayed_case(grade),
                grade,
                require_complete_delay=True,
            )
        )

    for index, (length, multiplicity) in enumerate(
        ((2, 2), (3, 3), (5, 2))
    ):
        data = gauged_shift(
            (length,) * multiplicity,
            729_900 + index,
        )
        records.append(
            audit_case(
                "repeated_monomial_apex",
                *data,
                length,
                require_complete_delay=True,
                require_apex=True,
            )
        )
    return records


def write_records(
    records: list[RetighteningTransportRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(
        json.dumps(asdict(record), sort_keys=True) + "\n"
        for record in records
    )
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run every deterministic audit."""

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
                "record_count": len(records),
                "records": [asdict(record) for record in records],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
