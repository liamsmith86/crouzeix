#!/usr/bin/env python3
"""Audit the universal polynomial preimage for every transfer grade.

For a balanced pure partial isometry ``S`` with defect frames ``V,W``,
L212 defines

    C_hat[k] = -(7/2) Q (
        S^k W B_k
        + sum_{j<k} (S*)^(k-j) V B_k* B_j
    ).

The lower-grade terms retain exactly the contamination removed by the
flag projection in L208.  The resulting physical endpoint is

    28 (channel(B_k* B_k) - B_k B_k*).

This checker tests each grade and weighted sums on structured and
unstructured partial-isometry colligations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_elliptic_cokernel import endpoint_motion
from repeated_crabb_elliptic_selection import (
    haar_unitary,
    random_partial_isometry,
)
from repeated_crabb_transfer_deflation import heterogeneous_shift
from repeated_crabb_transfer_flag import (
    transfer_channel,
    transfer_coefficient,
)


@dataclass(frozen=True)
class AllGradePreimageRecord:
    """One all-grade polynomial-preimage audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    maximum_grade: int
    spectral_radius: str
    maximum_perpendicular_error: str
    maximum_single_grade_error: str
    maximum_contamination_norm: str
    weighted_sum_error: str
    apex_correction_norm: str
    all_checks_passed: bool


def all_grade_column(
    operator: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    grade: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return L212's balanced column, coefficient, and contamination norm."""

    coefficient = transfer_coefficient(
        operator,
        right,
        left,
        grade,
    )
    column, contamination_norm = generalized_grade_column(
        operator,
        right,
        left,
        grade,
        coefficient,
        scale=-3.5,
    )
    return column, coefficient, contamination_norm


def generalized_grade_column(
    operator: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    grade: int,
    multiplier: np.ndarray,
    scale: float = 1.0,
) -> tuple[np.ndarray, float]:
    """Return the polarized L212 column for one copy multiplier."""

    dimension = len(operator)
    projection = np.eye(dimension) - right @ right.conj().T
    column_core = (
        np.linalg.matrix_power(operator, grade)
        @ left
        @ multiplier
    )
    contamination = np.zeros_like(column_core)
    for earlier_grade in range(1, grade):
        earlier = transfer_coefficient(
            operator,
            right,
            left,
            earlier_grade,
        )
        contamination += (
            np.linalg.matrix_power(
                operator.conj().T,
                grade - earlier_grade,
            )
            @ right
            @ multiplier.conj().T
            @ earlier
        )
    column = scale * projection @ (column_core + contamination)
    return column, float(np.linalg.norm(contamination))


def gauged_shift(
    channel_lengths: tuple[int, ...],
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an independently gauged heterogeneous delay sum."""

    operator, right, left = heterogeneous_shift(channel_lengths)
    generator = np.random.default_rng(seed)
    state_gauge = haar_unitary(len(operator), generator)
    right_gauge = haar_unitary(right.shape[1], generator)
    left_gauge = haar_unitary(left.shape[1], generator)
    return (
        state_gauge @ operator @ state_gauge.conj().T,
        state_gauge @ right @ right_gauge,
        state_gauge @ left @ left_gauge,
    )


def stable_random_partial_isometry(
    dimension: int,
    multiplicity: int,
    generator: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an unstructured partial isometry with radius below 0.97."""

    while True:
        data = random_partial_isometry(
            dimension,
            multiplicity,
            generator,
        )
        if np.max(np.abs(np.linalg.eigvals(data[0]))) < 0.97:
            return data


def audit_case(
    construction_kind: str,
    operator: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    maximum_grade: int,
) -> AllGradePreimageRecord:
    """Audit all individual grades and one weighted aggregate."""

    dimension = len(operator)
    multiplicity = right.shape[1]
    identity = np.eye(dimension)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    metric = 2 * identity - right_projection + 2 * left_projection
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    physical = np.linalg.inv(metric_root) @ operator @ metric_root

    columns: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    perpendicular_error = 0.0
    single_grade_error = 0.0
    contamination_norm = 0.0
    for grade in range(1, maximum_grade + 1):
        balanced_column, coefficient, contamination = all_grade_column(
            operator,
            right,
            left,
            grade,
        )
        physical_column = metric_root @ balanced_column
        target = 28 * (
            transfer_channel(
                operator,
                right,
                left,
                coefficient.conj().T @ coefficient,
            )
            - coefficient @ coefficient.conj().T
        )
        actual = endpoint_motion(
            physical,
            right,
            left,
            physical_column,
        )
        perpendicular_error = max(
            perpendicular_error,
            float(np.linalg.norm(right.conj().T @ physical_column)),
        )
        single_grade_error = max(
            single_grade_error,
            float(np.linalg.norm(actual - target)),
        )
        contamination_norm = max(contamination_norm, contamination)
        columns.append(physical_column)
        targets.append(target)

    reflection_parameter = 0.23
    weights = [
        reflection_parameter ** (2 * grade)
        for grade in range(1, maximum_grade + 1)
    ]
    aggregate_column = sum(
        weight * column
        for weight, column in zip(weights, columns, strict=True)
    )
    aggregate_target = sum(
        weight * target
        for weight, target in zip(weights, targets, strict=True)
    )
    weighted_error = float(
        np.linalg.norm(
            endpoint_motion(
                physical,
                right,
                left,
                aggregate_column,
            )
            - aggregate_target
        )
    )

    apex_norm = 0.0
    if construction_kind == "repeated_monomial_shift":
        apex_norm = max(float(np.linalg.norm(column)) for column in columns)

    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(operator)))
    )
    tolerance = 3e-10
    verified = bool(
        spectral_radius < 1
        and perpendicular_error < tolerance
        and single_grade_error < tolerance
        and weighted_error < tolerance
        and (
            construction_kind != "unstructured_partial_isometry"
            or contamination_norm > 1e-5
        )
        and (
            construction_kind != "repeated_monomial_shift"
            or apex_norm < tolerance
        )
    )
    if not verified:
        raise RuntimeError("the all-grade preimage audit failed")

    return AllGradePreimageRecord(
        construction_kind=construction_kind,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        maximum_grade=maximum_grade,
        spectral_radius=format_float(spectral_radius),
        maximum_perpendicular_error=format_float(perpendicular_error),
        maximum_single_grade_error=format_float(single_grade_error),
        maximum_contamination_norm=format_float(contamination_norm),
        weighted_sum_error=format_float(weighted_error),
        apex_correction_norm=format_float(apex_norm),
        all_checks_passed=verified,
    )


def standard_records() -> list[AllGradePreimageRecord]:
    """Return deterministic structured and unstructured audits."""

    records: list[AllGradePreimageRecord] = []
    generator = np.random.default_rng(71224)
    for dimension, multiplicity in ((7, 2), (9, 3), (12, 3), (14, 4)):
        for _ in range(2):
            data = stable_random_partial_isometry(
                dimension,
                multiplicity,
                generator,
            )
            records.append(
                audit_case(
                    "unstructured_partial_isometry",
                    *data,
                    maximum_grade=6,
                )
            )

    for index, lengths in enumerate(
        ((2, 3, 5), (2, 4, 6, 7), (3, 5, 6, 8))
    ):
        data = gauged_shift(lengths, 71324 + index)
        records.append(
            audit_case(
                "heterogeneous_delay_sum",
                *data,
                maximum_grade=max(lengths),
            )
        )

    for index, (length, multiplicity) in enumerate(
        ((2, 2), (3, 3), (5, 2))
    ):
        data = gauged_shift(
            (length,) * multiplicity,
            71424 + index,
        )
        records.append(
            audit_case(
                "repeated_monomial_shift",
                *data,
                maximum_grade=length,
            )
        )
    return records


def write_records(
    records: list[AllGradePreimageRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines records atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_all_grade_preimage_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the all-grade polynomial-preimage audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
