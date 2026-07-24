#!/usr/bin/env python3
"""Audit transfer-channel covariance under wandering-chain deflation.

L210 proves that wandering-chain deflation shifts the entire delayed
transfer row:

    B_tilde[n] = [U_perp* B[n]; U* B[n+k-1]].

This checker strengthens that result at the completely positive
transfer-channel level and audits the grade-one left Gram.  It uses
gauged delay sums and unstructured partial isometries enlarged by
non-coordinate wandering chains.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, solve_discrete_lyapunov

from crabb_block_hardy_equality import format_float
from repeated_crabb_transfer_deflation import (
    haar_unitary,
    heterogeneous_shift,
)


@dataclass(frozen=True)
class TransferChannelCovarianceRecord:
    """One all-coefficient and channel-covariance audit."""

    construction_kind: str
    state_dimension_before: int
    state_dimension_after: int
    defect_dimension: int
    flag_dimension: int
    deflation_grade: int
    coefficient_count: int
    transfer_coefficient_error: str
    positive_transfer_error: str
    reflected_transfer_error: str
    channel_covariance_error: str
    flagged_left_gram_error: str
    all_checks_passed: bool


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


def transfer_value(
    operator: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    value: complex,
) -> np.ndarray:
    """Evaluate ``W*(I-zS*)^(-1)V``."""

    identity = np.eye(len(operator), dtype=complex)
    return (
        left.conj().T
        @ np.linalg.solve(
            identity - value * operator.conj().T,
            right,
        )
    )


def transfer_channel(
    operator: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    copy_matrix: np.ndarray,
) -> np.ndarray:
    """Return ``sum B_n K B_n*`` by a Stein solve."""

    gramian = solve_discrete_lyapunov(
        operator.conj().T,
        right @ copy_matrix @ right.conj().T,
    )
    return left.conj().T @ gramian @ left


def random_partial_isometry(
    state_dimension: int,
    defect_dimension: int,
    generator: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an unstructured partial isometry with orthogonal defects."""

    identity = np.eye(state_dimension, dtype=complex)
    right = identity[:, :defect_dimension]
    left = identity[:, -defect_dimension:]
    domain = identity[:, defect_dimension:]
    range_frame = identity[:, : state_dimension - defect_dimension]
    bridge = haar_unitary(
        state_dimension - defect_dimension,
        generator,
    )
    operator = range_frame @ bridge @ domain.conj().T
    return operator, right, left


def inflate_wandering_chain(
    base_operator: np.ndarray,
    base_right: np.ndarray,
    base_left: np.ndarray,
    grade: int,
    flag_dimension: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Attach ``grade-1`` input stages to part of the left defect."""

    base_dimension, multiplicity = base_right.shape
    added_dimension = (grade - 1) * flag_dimension
    dimension = base_dimension + added_dimension
    operator = np.zeros((dimension, dimension), dtype=complex)
    operator[:base_dimension, :base_dimension] = base_operator

    promoted = base_left[:, multiplicity - flag_dimension :]
    delay_frames = []
    for level in range(grade - 1):
        frame = np.zeros((dimension, flag_dimension), dtype=complex)
        start = base_dimension + level * flag_dimension
        frame[start : start + flag_dimension, :] = np.eye(flag_dimension)
        delay_frames.append(frame)

    promoted_full = np.zeros((dimension, flag_dimension), dtype=complex)
    promoted_full[:base_dimension, :] = promoted
    chain = [*delay_frames, promoted_full]
    for source, target in zip(chain[:-1], chain[1:], strict=True):
        operator += target @ source.conj().T

    right = np.zeros((dimension, multiplicity), dtype=complex)
    right[:base_dimension, :] = base_right
    unpromoted = np.zeros(
        (dimension, multiplicity - flag_dimension),
        dtype=complex,
    )
    unpromoted[:base_dimension, :] = base_left[
        :, : multiplicity - flag_dimension
    ]
    left = np.hstack((unpromoted, delay_frames[0]))
    selector = np.eye(multiplicity, dtype=complex)[
        :, multiplicity - flag_dimension :
    ]
    return operator, right, left, selector


def gauged_delay_case(
    channel_lengths: tuple[int, ...],
    grade: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return one independently gauged heterogeneous-delay flag."""

    operator, right, left = heterogeneous_shift(channel_lengths)
    dimension = len(operator)
    multiplicity = len(channel_lengths)
    generator = np.random.default_rng(seed)
    state_gauge = haar_unitary(dimension, generator)
    right_gauge = haar_unitary(multiplicity, generator)
    left_gauge = haar_unitary(multiplicity, generator)
    operator = state_gauge @ operator @ state_gauge.conj().T
    right = state_gauge @ right @ right_gauge
    left = state_gauge @ left @ left_gauge
    selected = [
        index
        for index, length in enumerate(channel_lengths)
        if length >= grade
    ]
    selector = np.eye(multiplicity, dtype=complex)[:, selected]
    flag = left_gauge.conj().T @ selector
    return operator, right, left, flag


def inflated_case(
    state_dimension: int,
    multiplicity: int,
    grade: int,
    flag_dimension: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return an unstructured colligation with an attached delay flag."""

    generator = np.random.default_rng(seed)
    base = random_partial_isometry(
        state_dimension,
        multiplicity,
        generator,
    )
    operator, right, left, flag = inflate_wandering_chain(
        *base,
        grade,
        flag_dimension,
    )

    dimension = len(operator)
    state_gauge = haar_unitary(dimension, generator)
    right_gauge = haar_unitary(multiplicity, generator)
    left_gauge = haar_unitary(multiplicity, generator)
    operator = state_gauge @ operator @ state_gauge.conj().T
    right = state_gauge @ right @ right_gauge
    left = state_gauge @ left @ left_gauge
    flag = left_gauge.conj().T @ flag
    return operator, right, left, flag


def make_record(
    construction_kind: str,
    operator: np.ndarray,
    right: np.ndarray,
    left: np.ndarray,
    flag: np.ndarray,
    grade: int,
    seed: int,
) -> TransferChannelCovarianceRecord:
    """Deflate one flag and audit every covariance identity."""

    dimension, multiplicity = right.shape
    flag_dimension = flag.shape[1]
    flag_complement = null_space(flag.conj().T)
    delay_columns = [
        np.linalg.matrix_power(operator, degree) @ left @ flag
        for degree in range(grade)
    ]
    removed = np.hstack(delay_columns[:-1])
    retained_basis = null_space(removed.conj().T)
    deflated = retained_basis.conj().T @ operator @ retained_basis
    deflated_right = retained_basis.conj().T @ right
    deflated_left = retained_basis.conj().T @ np.hstack(
        (
            left @ flag_complement,
            delay_columns[-1],
        )
    )

    coefficient_count = 9
    original_coefficients = [
        transfer_coefficient(operator, right, left, degree)
        for degree in range(coefficient_count + grade)
    ]
    deflated_coefficients = [
        transfer_coefficient(
            deflated,
            deflated_right,
            deflated_left,
            degree,
        )
        for degree in range(coefficient_count)
    ]
    coefficient_error = 0.0
    for degree, coefficient in enumerate(deflated_coefficients):
        expected = np.vstack(
            (
                flag_complement.conj().T
                @ original_coefficients[degree],
                flag.conj().T
                @ original_coefficients[degree + grade - 1],
            )
        )
        coefficient_error = max(
            coefficient_error,
            float(np.linalg.norm(coefficient - expected)),
        )

    bottom = np.eye(multiplicity, dtype=complex)[
        multiplicity - flag_dimension :,
        :,
    ]
    positive_error = 0.0
    reflected_error = 0.0
    ellipse = 0.19
    for index in range(13):
        value = 0.63 * np.exp(2j * np.pi * (index + 0.17) / 13)
        original_value = transfer_value(
            operator,
            right,
            left,
            value,
        )
        deflated_value = transfer_value(
            deflated,
            deflated_right,
            deflated_left,
            value,
        )
        positive_error = max(
            positive_error,
            float(
                np.linalg.norm(
                    flag.conj().T @ original_value
                    - value ** (grade - 1)
                    * bottom
                    @ deflated_value
                )
            ),
        )

        zeta = np.exp(2j * np.pi * (index + 0.29) / 13)
        reflected_value = ellipse / zeta
        original_reflection = transfer_value(
            operator,
            right,
            left,
            reflected_value,
        )
        deflated_reflection = transfer_value(
            deflated,
            deflated_right,
            deflated_left,
            reflected_value,
        )
        reflected_error = max(
            reflected_error,
            float(
                np.linalg.norm(
                    flag.conj().T @ original_reflection
                    - ellipse ** (grade - 1)
                    * zeta ** (1 - grade)
                    * bottom
                    @ deflated_reflection
                )
            ),
        )

    generator = np.random.default_rng(seed + 10_000)
    factor = (
        generator.standard_normal((multiplicity, multiplicity))
        + 1j
        * generator.standard_normal((multiplicity, multiplicity))
    )
    copy_matrix = factor @ factor.conj().T
    original_channel = transfer_channel(
        operator,
        right,
        left,
        copy_matrix,
    )
    deflated_channel = transfer_channel(
        deflated,
        deflated_right,
        deflated_left,
        copy_matrix,
    )
    channel_error = float(
        np.linalg.norm(
            bottom @ deflated_channel @ bottom.conj().T
            - flag.conj().T @ original_channel @ flag
        )
    )

    original_grade = original_coefficients[grade]
    deflated_first = deflated_coefficients[1]
    gram_error = float(
        np.linalg.norm(
            bottom
            @ deflated_first
            @ deflated_first.conj().T
            @ bottom.conj().T
            - flag.conj().T
            @ original_grade
            @ original_grade.conj().T
            @ flag
        )
    )

    tolerance = 3e-10
    verified = bool(
        coefficient_error < tolerance
        and positive_error < tolerance
        and reflected_error < tolerance
        and channel_error < tolerance
        and gram_error < tolerance
    )
    if not verified:
        raise RuntimeError("the transfer-covariance audit failed")
    return TransferChannelCovarianceRecord(
        construction_kind=construction_kind,
        state_dimension_before=dimension,
        state_dimension_after=len(deflated),
        defect_dimension=multiplicity,
        flag_dimension=flag_dimension,
        deflation_grade=grade,
        coefficient_count=coefficient_count,
        transfer_coefficient_error=format_float(coefficient_error),
        positive_transfer_error=format_float(positive_error),
        reflected_transfer_error=format_float(reflected_error),
        channel_covariance_error=format_float(channel_error),
        flagged_left_gram_error=format_float(gram_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[TransferChannelCovarianceRecord]:
    """Return deterministic structured and unstructured audits."""

    records: list[TransferChannelCovarianceRecord] = []
    delay_specs = (
        ((2, 3, 5), 2),
        ((2, 3, 5), 3),
        ((2, 4, 6, 7), 4),
        ((3, 5, 6, 8), 5),
    )
    for index, (lengths, grade) in enumerate(delay_specs):
        data = gauged_delay_case(lengths, grade, 71024 + index)
        records.append(
            make_record(
                "gauged_delay_sum",
                *data,
                grade,
                71024 + index,
            )
        )

    inflated_specs = (
        (8, 2, 2, 1),
        (10, 3, 2, 2),
        (10, 3, 3, 1),
        (12, 4, 4, 2),
        (12, 4, 5, 1),
    )
    for index, (dimension, multiplicity, grade, flag_dimension) in enumerate(
        inflated_specs
    ):
        data = inflated_case(
            dimension,
            multiplicity,
            grade,
            flag_dimension,
            71124 + index,
        )
        records.append(
            make_record(
                "inflated_unstructured",
                *data,
                grade,
                71124 + index,
            )
        )
    return records


def write_records(
    records: list[TransferChannelCovarianceRecord],
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
            "repeated_crabb_transfer_channel_covariance_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the standard transfer-covariance audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
