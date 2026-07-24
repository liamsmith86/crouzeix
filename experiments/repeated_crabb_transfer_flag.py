#!/usr/bin/env python3
"""Audit the higher transfer-intertwining and flagged coboundary identity.

For the genuine transfer coefficients ``B_k = W* (S*)^k V``, L208
identifies the complete failure of an arbitrary dual Gramian to
intertwine ``B_k``.  On the common left kernel of the earlier
coefficients, the lower-grade terms vanish.  Consequently every
compressed target of the form

    eta P (C(B_k* P B_k) - B_k B_k*) P

has the explicit endpoint-map preimage

    C_hat = -(eta/8) Q S^k W P B_k.

The checker tests the general intertwining identity on unstructured
partial isometries and the nontrivial flagged identity on delayed,
independently gauged matrix-inner channels.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, solve_discrete_lyapunov, sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_elliptic_cokernel import endpoint_motion
from repeated_crabb_elliptic_selection import (
    haar_unitary,
    random_partial_isometry,
)


@dataclass(frozen=True)
class TransferFlagRecord:
    """One higher-coefficient flag audit."""

    anchor_kind: str
    state_dimension: int
    defect_dimension: int
    grade: int
    flag_dimension: int
    spectral_radius: str
    intertwining_defect_error: str
    earlier_flag_annihilation_error: str
    perpendicular_column_error: str
    flagged_channel_preimage_error: str
    all_checks_passed: bool


def transfer_coefficient(
    balanced: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    grade: int,
) -> np.ndarray:
    """Return ``B_grade = W* (S*)^grade V``."""

    return (
        left.conj().T
        @ np.linalg.matrix_power(balanced.conj().T, grade)
        @ right
    )


def transfer_channel(
    balanced: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    copy_matrix: np.ndarray,
) -> np.ndarray:
    """Apply the completely positive transfer channel."""

    gramian = solve_discrete_lyapunov(
        balanced.conj().T,
        right @ copy_matrix @ right.conj().T,
    )
    return left.conj().T @ gramian @ left


def common_left_kernel_projection(
    coefficients: list[np.ndarray],
    multiplicity: int,
) -> np.ndarray:
    """Project onto the common kernel of the coefficient adjoints."""

    if not coefficients:
        return np.eye(multiplicity, dtype=complex)
    stacked = np.vstack(
        [coefficient.conj().T for coefficient in coefficients]
    )
    kernel = null_space(stacked, rcond=1e-10)
    return kernel @ kernel.conj().T


def delayed_partial_isometry(
    delays: tuple[int, ...],
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return independently gauged direct-sum delay channels."""

    multiplicity = len(delays)
    dimension = sum(delay + 1 for delay in delays)
    balanced = np.zeros((dimension, dimension), dtype=complex)
    right_base = np.zeros((dimension, multiplicity), dtype=complex)
    left_base = np.zeros((dimension, multiplicity), dtype=complex)
    offset = 0
    for channel, delay in enumerate(delays):
        right_base[offset, channel] = 1
        left_base[offset + delay, channel] = 1
        for level in range(delay):
            balanced[offset + level, offset + level + 1] = 1
        offset += delay + 1
    right = right_base @ haar_unitary(multiplicity, rng)
    left = left_base @ haar_unitary(multiplicity, rng)
    return balanced, right, left


def deterministic_hermitian(
    size: int,
    grade: int,
) -> np.ndarray:
    """Return a deterministic dense Hermitian endpoint test."""

    matrix = np.zeros((size, size), dtype=complex)
    for row in range(size):
        for column in range(size):
            matrix[row, column] = (
                (1 + 2 * row + 3 * column + grade) / (4 * size)
                + 1j
                * (row - 2 * column + grade)
                / (5 * size)
            )
    return matrix + matrix.conj().T


def audit_grade(
    balanced: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    grade: int,
    *,
    anchor_kind: str,
) -> TransferFlagRecord:
    """Audit the general and flagged identities at one grade."""

    dimension = len(balanced)
    multiplicity = right.shape[1]
    identity = np.eye(dimension, dtype=complex)
    projection = identity - right @ right.conj().T
    coefficients = [
        transfer_coefficient(
            balanced,
            right,
            left,
            current_grade,
        )
        for current_grade in range(grade + 1)
    ]

    endpoint_test = deterministic_hermitian(multiplicity, grade)
    dual = solve_discrete_lyapunov(
        balanced,
        left @ endpoint_test @ left.conj().T,
    )
    right_action = right.conj().T @ dual @ right
    residual = projection @ dual @ right
    expected_defect = (
        left.conj().T
        @ np.linalg.matrix_power(balanced.conj().T, grade)
        @ residual
    )
    for earlier_grade in range(1, grade):
        expected_defect += (
            coefficients[earlier_grade]
            @ residual.conj().T
            @ np.linalg.matrix_power(
                balanced.conj().T,
                grade - earlier_grade,
            )
            @ right
        )
    intertwining_error = float(
        np.linalg.norm(
            endpoint_test @ coefficients[grade]
            - coefficients[grade] @ right_action
            - expected_defect
        )
    )

    flag = common_left_kernel_projection(
        coefficients[1:grade],
        multiplicity,
    )
    flag_dimension = int(round(float(np.trace(flag).real)))
    earlier_error = max(
        [
            float(np.linalg.norm(flag @ coefficient))
            for coefficient in coefficients[1:grade]
        ]
        or [0.0]
    )

    metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    operator = np.linalg.inv(metric_root) @ balanced @ metric_root
    eta = 1.375
    coefficient = coefficients[grade]
    balanced_column = (
        -eta
        * projection
        @ np.linalg.matrix_power(balanced, grade)
        @ left
        @ flag
        @ coefficient
        / 8
    )
    physical_column = metric_root @ balanced_column
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ physical_column)
    )
    achieved = endpoint_motion(
        operator,
        right,
        left,
        physical_column,
    )
    right_gram = coefficient.conj().T @ flag @ coefficient
    channel = transfer_channel(
        balanced,
        right,
        left,
        right_gram,
    )
    target = (
        eta
        * flag
        @ (channel - coefficient @ coefficient.conj().T)
        @ flag
    )
    flagged_error = float(
        np.linalg.norm(flag @ achieved @ flag - target)
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(balanced)))
    )
    verified = bool(
        spectral_radius < 1
        and intertwining_error < 2e-9
        and earlier_error < 2e-9
        and perpendicular_error < 2e-9
        and flagged_error < 2e-9
    )
    if not verified:
        raise RuntimeError("the higher transfer-flag audit failed")

    return TransferFlagRecord(
        anchor_kind=anchor_kind,
        state_dimension=dimension,
        defect_dimension=multiplicity,
        grade=grade,
        flag_dimension=flag_dimension,
        spectral_radius=format_float(spectral_radius),
        intertwining_defect_error=format_float(intertwining_error),
        earlier_flag_annihilation_error=format_float(earlier_error),
        perpendicular_column_error=format_float(perpendicular_error),
        flagged_channel_preimage_error=format_float(flagged_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[TransferFlagRecord]:
    """Return unstructured and nontrivially flagged standard records."""

    rng = np.random.default_rng(70224)
    records: list[TransferFlagRecord] = []
    for dimension, multiplicity in ((6, 2), (8, 2), (9, 3), (12, 3)):
        accepted = 0
        while accepted < 2:
            balanced, right, left = random_partial_isometry(
                dimension,
                multiplicity,
                rng,
            )
            if np.max(np.abs(np.linalg.eigvals(balanced))) >= 0.97:
                continue
            for grade in range(1, 6):
                records.append(
                    audit_grade(
                        balanced,
                        right,
                        left,
                        grade,
                        anchor_kind="unstructured_partial_isometry",
                    )
                )
            accepted += 1

    for delays in (
        (1, 2, 3),
        (2, 3, 4),
        (1, 3, 5, 5),
        (2, 2, 4, 6),
    ):
        balanced, right, left = delayed_partial_isometry(delays, rng)
        for grade in range(1, max(delays) + 1):
            records.append(
                audit_grade(
                    balanced,
                    right,
                    left,
                    grade,
                    anchor_kind="gauged_delay_flag",
                )
            )
    return records


def write_records(
    records: list[TransferFlagRecord],
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
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_transfer_flag_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the higher transfer-flag audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
