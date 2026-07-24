#!/usr/bin/env python3
"""Probe scalar equality channels in the repeated Crabb inner transfer.

For a matrix-inner transfer

    B(z) = sum B_n z^n,

scalar Crouzeix equality requires unit vectors ``v,u`` for which every
``B_n v`` lies on the same line ``C u``.  Equivalently,

    sum_n |u* B_n v|^2 = 1.

This script maximizes that product-state energy by deterministic
alternating eigenvector iteration.  It checks exact score one at the
repeated Crabb apex and a strict numerical gap at the standard
genuinely noncommuting anchors.  The optimization is a probe, not a
certificate that no other product state attains one.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_inner_faber_transfer import (
    canonical_transfer_data,
    strengthened_inverse_toeplitz,
    transfer_coefficients,
)


@dataclass(frozen=True)
class ScalarChannelRecord:
    """One scalar-channel product-state audit."""

    length: int
    multiplicity: int
    toeplitz_strength: str
    coefficient_count: int
    parseval_error: str
    apex_channel_error: str
    noncommuting_channel_score: str
    noncommuting_channel_gap: str
    maximizing_leakage: str
    all_checks_passed: bool


def product_state_step(
    coefficients: np.ndarray,
    right: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Perform one right/left alternating eigenvector step."""

    images = np.einsum("kij,j->ki", coefficients, right)
    left_gram = images.T @ images.conj()
    _, left_vectors = np.linalg.eigh(left_gram)
    left = left_vectors[:, -1]

    adjoint_images = np.einsum(
        "kji,j->ki",
        coefficients.conj(),
        left,
    )
    right_gram = adjoint_images.T @ adjoint_images.conj()
    eigenvalues, right_vectors = np.linalg.eigh(right_gram)
    return left, right_vectors[:, -1], float(eigenvalues[-1])


def maximize_product_state_energy(
    coefficients: tuple[np.ndarray, ...],
    restarts: int = 48,
    iterations: int = 80,
) -> tuple[float, float]:
    """Return a deterministic alternating maximum and its leakage."""

    coefficient_array = np.asarray(coefficients)
    multiplicity = coefficient_array.shape[1]
    generator = np.random.default_rng(70224)
    best_score = -np.inf
    best_pair: tuple[np.ndarray, np.ndarray] | None = None

    starts = [np.eye(multiplicity, dtype=complex)[:, index]
              for index in range(multiplicity)]
    for _ in range(restarts):
        vector = (
            generator.standard_normal(multiplicity)
            + 1j * generator.standard_normal(multiplicity)
        )
        starts.append(vector / np.linalg.norm(vector))

    for start in starts:
        right = start
        score = 0.0
        left = np.zeros(multiplicity, dtype=complex)
        for _ in range(iterations):
            left, following, score = product_state_step(
                coefficient_array,
                right,
            )
            overlap = abs(np.vdot(following, right))
            right = following
            if 1 - overlap < 1e-14:
                break
        if score > best_score:
            best_score = score
            best_pair = (left, right)

    if best_pair is None:
        raise RuntimeError("the product-state iteration produced no pair")
    left, right = best_pair
    images = np.einsum("kij,j->ki", coefficient_array, right)
    scalar_parts = images @ left.conj()
    leakage = images - scalar_parts[:, None] * left[None, :]
    return best_score, float(np.linalg.norm(leakage))


def coefficient_tuple(
    length: int,
    multiplicity: int,
    strength: float,
    count: int,
) -> tuple[tuple[np.ndarray, ...], float]:
    """Return transfer coefficients for one noncommuting anchor."""

    inverse, actual_strength = strengthened_inverse_toeplitz(
        length,
        multiplicity,
        strength,
    )
    data = canonical_transfer_data(
        np.linalg.inv(inverse),
        length,
        multiplicity,
    )
    return transfer_coefficients(data, count)[1:], actual_strength


def make_record(
    length: int,
    multiplicity: int,
    strength: float,
    count: int,
) -> ScalarChannelRecord:
    """Build one scalar-channel audit record."""

    coefficients, actual_strength = coefficient_tuple(
        length,
        multiplicity,
        strength,
        count,
    )
    identity = np.eye(multiplicity, dtype=complex)
    parseval = sum(
        coefficient.conj().T @ coefficient
        for coefficient in coefficients
    )
    parseval_error = float(np.linalg.norm(parseval - identity))
    score, leakage = maximize_product_state_energy(coefficients)

    apex = canonical_transfer_data(
        np.eye(length * multiplicity, dtype=complex) / 2,
        length,
        multiplicity,
    )
    apex_coefficients = transfer_coefficients(apex, count)[1:]
    apex_score, apex_leakage = maximize_product_state_energy(
        apex_coefficients
    )
    apex_error = max(abs(1 - apex_score), apex_leakage)
    gap = 1 - score

    verified = bool(
        parseval_error < 3e-10
        and apex_error < 3e-10
        and gap > 1e-7
        and leakage > 1e-5
    )
    if not verified:
        raise RuntimeError("the scalar-channel audit failed")
    return ScalarChannelRecord(
        length=length,
        multiplicity=multiplicity,
        toeplitz_strength=format_float(actual_strength),
        coefficient_count=count - 1,
        parseval_error=format_float(parseval_error),
        apex_channel_error=format_float(apex_error),
        noncommuting_channel_score=format_float(score),
        noncommuting_channel_gap=format_float(gap),
        maximizing_leakage=format_float(leakage),
        all_checks_passed=verified,
    )


def write_records(
    records: list[ScalarChannelRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines records."""

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
            "experiments/repeated_crabb_scalar_channel_s70224.jsonl"
        ),
    )
    parser.add_argument("--strength", type=float, default=18.0)
    parser.add_argument("--coefficient-count", type=int, default=120)
    return parser.parse_args()


def main() -> None:
    """Run the standard scalar-channel audit."""

    args = parse_args()
    records = [
        make_record(
            length,
            multiplicity,
            args.strength,
            args.coefficient_count,
        )
        for multiplicity in (2, 3)
        for length in range(2, 6)
    ]
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
