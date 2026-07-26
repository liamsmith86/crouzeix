#!/usr/bin/env python3
"""Audit Hankel-sandwich normal ordering into the leakage ideal."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_hardy_two_frame import (
    analysis_rows,
    block_backward_shift,
)
from repeated_crabb_toeplitz_window_energy import (
    toeplitz_window,
    transfer_coefficients,
)
from repeated_crabb_transfer_deflation import (
    haar_unitary,
    heterogeneous_shift,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class HankelLeakageIdealRecord:
    """One finite exact-window audit of the sandwich identities."""

    channel_lengths: tuple[int, ...]
    defect_dimension: int
    hardy_rows: int
    right_isometry_error: str
    left_isometry_error: str
    hankel_intertwining_error: str
    shifted_model_projection_error: str
    shifted_complement_error: str
    maximum_sandwich_error: str
    maximum_ideal_remainder_error: str
    all_checks_passed: bool


def compressed_norm(matrix: Matrix, retained: int) -> float:
    """Return the norm of the safe leading square compression."""

    return float(np.linalg.norm(matrix[:retained, :retained]))


def audit_case(
    channel_lengths: tuple[int, ...],
    seed: int,
) -> HankelLeakageIdealRecord:
    """Audit one independently gauged finite inner transfer."""

    partial, right, left = heterogeneous_shift(channel_lengths)
    multiplicity = len(channel_lengths)
    generator = np.random.default_rng(seed)
    right = right @ haar_unitary(multiplicity, generator)
    left = left @ haar_unitary(multiplicity, generator)

    maximum_length = max(channel_lengths)
    row_count = maximum_length + 9
    right_rows, left_rows = analysis_rows(
        partial,
        right,
        left,
        row_count,
    )
    right_analysis = np.vstack(right_rows)
    left_analysis = np.vstack(left_rows)
    state_identity = np.eye(len(partial), dtype=complex)
    right_isometry_error = float(
        np.linalg.norm(
            right_analysis.conj().T @ right_analysis - state_identity
        )
    )
    left_isometry_error = float(
        np.linalg.norm(
            left_analysis.conj().T @ left_analysis - state_identity
        )
    )

    backward = block_backward_shift(row_count, multiplicity)
    forward = backward.conj().T
    hardy_identity = np.eye(row_count * multiplicity, dtype=complex)
    hankel = right_analysis @ left_analysis.conj().T
    shifted_left = backward @ left_analysis
    shifted_hankel = hankel @ forward
    shifted_model = shifted_left @ shifted_left.conj().T

    coefficients = transfer_coefficients(
        partial,
        right,
        left,
        row_count - 1,
    )
    toeplitz = toeplitz_window(coefficients, row_count - 1)
    leakage = toeplitz @ toeplitz.conj().T

    safe_rows = row_count - maximum_length - 3
    retained = safe_rows * multiplicity
    intertwining_error = compressed_norm(
        backward @ shifted_hankel - shifted_hankel @ forward,
        retained,
    )
    shifted_projection_error = compressed_norm(
        shifted_hankel.conj().T @ shifted_hankel - shifted_model,
        retained,
    )
    complement_error = compressed_norm(
        shifted_model + leakage - hardy_identity,
        retained,
    )

    maximum_sandwich_error = 0.0
    maximum_ideal_error = 0.0
    for forward_power in range(5):
        for backward_power in range(5):
            middle = (
                np.linalg.matrix_power(forward, forward_power)
                @ np.linalg.matrix_power(backward, backward_power)
            )
            sandwich = (
                shifted_hankel.conj().T
                @ middle
                @ shifted_hankel
            )
            scalar_left = np.linalg.matrix_power(
                backward,
                forward_power,
            )
            scalar_right = np.linalg.matrix_power(
                forward,
                backward_power,
            )
            expected = scalar_left @ shifted_model @ scalar_right
            scalar_background = scalar_left @ scalar_right
            ideal_remainder = (
                -scalar_left @ leakage @ scalar_right
            )
            maximum_sandwich_error = max(
                maximum_sandwich_error,
                compressed_norm(sandwich - expected, retained),
            )
            maximum_ideal_error = max(
                maximum_ideal_error,
                compressed_norm(
                    sandwich - scalar_background - ideal_remainder,
                    retained,
                ),
            )

    tolerance = 2e-11
    verified = bool(
        max(
            right_isometry_error,
            left_isometry_error,
            intertwining_error,
            shifted_projection_error,
            complement_error,
            maximum_sandwich_error,
            maximum_ideal_error,
        )
        < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the Hankel leakage-ideal audit failed: "
            f"lengths={channel_lengths}"
        )
    return HankelLeakageIdealRecord(
        channel_lengths=channel_lengths,
        defect_dimension=multiplicity,
        hardy_rows=row_count,
        right_isometry_error=format_float(right_isometry_error),
        left_isometry_error=format_float(left_isometry_error),
        hankel_intertwining_error=format_float(intertwining_error),
        shifted_model_projection_error=format_float(
            shifted_projection_error
        ),
        shifted_complement_error=format_float(complement_error),
        maximum_sandwich_error=format_float(maximum_sandwich_error),
        maximum_ideal_remainder_error=format_float(maximum_ideal_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[HankelLeakageIdealRecord]:
    """Return deterministic scalar and block-channel audits."""

    cases = (
        ((2,), 272_001),
        ((2, 4), 272_002),
        ((3, 4, 6), 272_003),
        ((2, 3, 5, 7), 272_004),
    )
    return [
        audit_case(channel_lengths, seed)
        for channel_lengths, seed in cases
    ]


def write_records(
    records: list[HankelLeakageIdealRecord],
    output: Path,
) -> str:
    """Write deterministic records atomically and return their hash."""

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
            "repeated_crabb_hankel_leakage_ideal_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
