#!/usr/bin/env python3
"""Audit wandering-chain deflation of a matrix-inner transfer flag.

Let ``S`` be a finite pure partial isometry with right and left defect
frames ``V,W`` and coefficients ``B_n=W*(S*)^nV``.  If a left copy
subspace ``U`` satisfies ``B_j*U=0`` for ``j<k``, then

    WU, S WU, ..., S^(k-2) WU

is an orthonormal co-invariant delay line.  Compressing it away gives
another partial isometry.  Its new left defect replaces ``WU`` by
``S^(k-1)WU``, and its first transfer coefficient on that block is the
old ``B_k``.

The checker uses heterogeneous direct sums of finite shifts, then
applies independent random state and defect-frame unitaries so that
the tested flag subspaces are not coordinate subspaces.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from crabb_block_hardy_equality import format_float


@dataclass(frozen=True)
class TransferDeflationRecord:
    """One exact delay-line deflation audit."""

    channel_lengths: tuple[int, ...]
    deflation_grade: int
    flag_dimension: int
    state_dimension_before: int
    state_dimension_after: int
    lower_coefficient_kernel_error: str
    delay_isometry_error: str
    delay_orthogonality_error: str
    coinvariance_error: str
    compressed_right_defect_error: str
    compressed_left_defect_error: str
    first_transfer_promotion_error: str
    all_checks_passed: bool


def haar_unitary(
    size: int,
    generator: np.random.Generator,
) -> np.ndarray:
    """Return a deterministic-seed Haar unitary."""

    matrix = (
        generator.standard_normal((size, size))
        + 1j * generator.standard_normal((size, size))
    )
    unitary, triangular = np.linalg.qr(matrix)
    diagonal = np.diag(triangular)
    phases = np.ones(size, dtype=complex)
    nonzero = np.abs(diagonal) > 0
    phases[nonzero] = (
        np.conjugate(diagonal[nonzero])
        / np.abs(diagonal[nonzero])
    )
    return unitary @ np.diag(phases)


def heterogeneous_shift(
    channel_lengths: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return direct-sum shifts and their endpoint defect frames."""

    multiplicity = len(channel_lengths)
    dimension = sum(length + 1 for length in channel_lengths)
    operator = np.zeros((dimension, dimension), dtype=complex)
    right = np.zeros((dimension, multiplicity), dtype=complex)
    left = np.zeros((dimension, multiplicity), dtype=complex)
    offset = 0
    for channel, length in enumerate(channel_lengths):
        right[offset, channel] = 1
        left[offset + length, channel] = 1
        for level in range(1, length + 1):
            operator[offset + level - 1, offset + level] = 1
        offset += length + 1
    return operator, right, left


def transfer_coefficient(
    operator: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    degree: int,
) -> np.ndarray:
    """Return ``W*(S*)^degree V``."""

    return (
        left.conj().T
        @ np.linalg.matrix_power(operator.conj().T, degree)
        @ right
    )


def make_record(
    channel_lengths: tuple[int, ...],
    grade: int,
    seed: int,
) -> TransferDeflationRecord:
    """Audit one non-coordinate delayed copy flag."""

    base, base_right, base_left = heterogeneous_shift(channel_lengths)
    dimension = len(base)
    multiplicity = len(channel_lengths)
    generator = np.random.default_rng(seed)

    state_unitary = haar_unitary(dimension, generator)
    right_unitary = haar_unitary(multiplicity, generator)
    left_unitary = haar_unitary(multiplicity, generator)
    operator = state_unitary @ base @ state_unitary.conj().T
    right = state_unitary @ base_right @ right_unitary
    left = state_unitary @ base_left @ left_unitary

    selected = [
        index
        for index, length in enumerate(channel_lengths)
        if length >= grade
    ]
    if not selected:
        raise ValueError("the requested flag is empty")
    selector = np.eye(multiplicity, dtype=complex)[:, selected]
    flag = left_unitary.conj().T @ selector
    flag_complement = null_space(flag.conj().T)
    flag_dimension = flag.shape[1]

    coefficients = [
        transfer_coefficient(operator, right, left, degree)
        for degree in range(grade + 1)
    ]
    kernel_error = max(
        float(
            np.linalg.norm(
                coefficients[degree].conj().T @ flag
            )
        )
        for degree in range(grade)
    )

    delay_columns = [
        np.linalg.matrix_power(operator, degree) @ left @ flag
        for degree in range(grade)
    ]
    delay_isometry_error = max(
        float(
            np.linalg.norm(
                column.conj().T @ column
                - np.eye(flag_dimension)
            )
        )
        for column in delay_columns
    )
    delay_orthogonality_error = 0.0
    for row, first in enumerate(delay_columns):
        for column, second in enumerate(delay_columns):
            if row == column:
                continue
            delay_orthogonality_error = max(
                delay_orthogonality_error,
                float(np.linalg.norm(first.conj().T @ second)),
            )

    removed = np.hstack(delay_columns[:-1])
    removed_projection = removed @ removed.conj().T
    retained_projection = np.eye(dimension) - removed_projection
    coinvariance_error = float(
        np.linalg.norm(
            retained_projection
            @ operator.conj().T
            @ removed
        )
    )

    retained_basis = null_space(removed.conj().T)
    compressed = (
        retained_basis.conj().T @ operator @ retained_basis
    )
    compressed_right = retained_basis.conj().T @ right
    promoted_left = np.hstack(
        (
            left @ flag_complement,
            delay_columns[-1],
        )
    )
    compressed_left = retained_basis.conj().T @ promoted_left
    retained_identity = np.eye(len(compressed), dtype=complex)
    right_error = float(
        np.linalg.norm(
            retained_identity
            - compressed.conj().T @ compressed
            - compressed_right @ compressed_right.conj().T
        )
    )
    left_error = float(
        np.linalg.norm(
            retained_identity
            - compressed @ compressed.conj().T
            - compressed_left @ compressed_left.conj().T
        )
    )

    promoted_first = (
        compressed_left.conj().T
        @ compressed.conj().T
        @ compressed_right
    )
    expected_first = np.vstack(
        (
            flag_complement.conj().T @ coefficients[1],
            flag.conj().T @ coefficients[grade],
        )
    )
    promotion_error = float(
        np.linalg.norm(promoted_first - expected_first)
    )

    tolerance = 2e-10
    verified = bool(
        kernel_error < tolerance
        and delay_isometry_error < tolerance
        and delay_orthogonality_error < tolerance
        and coinvariance_error < tolerance
        and right_error < tolerance
        and left_error < tolerance
        and promotion_error < tolerance
    )
    if not verified:
        raise RuntimeError("the transfer-deflation audit failed")

    return TransferDeflationRecord(
        channel_lengths=channel_lengths,
        deflation_grade=grade,
        flag_dimension=flag_dimension,
        state_dimension_before=dimension,
        state_dimension_after=len(compressed),
        lower_coefficient_kernel_error=format_float(kernel_error),
        delay_isometry_error=format_float(delay_isometry_error),
        delay_orthogonality_error=format_float(
            delay_orthogonality_error
        ),
        coinvariance_error=format_float(coinvariance_error),
        compressed_right_defect_error=format_float(right_error),
        compressed_left_defect_error=format_float(left_error),
        first_transfer_promotion_error=format_float(promotion_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[TransferDeflationRecord]:
    """Return deterministic heterogeneous-delay audits."""

    specifications = (
        ((2, 3), 2),
        ((2, 3), 3),
        ((2, 4, 5), 2),
        ((2, 4, 5), 4),
        ((3, 5, 6), 3),
        ((3, 5, 6), 5),
        ((2, 3, 5, 7), 3),
        ((2, 3, 5, 7), 5),
    )
    return [
        make_record(lengths, grade, 70224 + index)
        for index, (lengths, grade) in enumerate(specifications)
    ]


def write_records(
    records: list[TransferDeflationRecord],
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
            "experiments/"
            "repeated_crabb_transfer_deflation_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the standard transfer-deflation audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
