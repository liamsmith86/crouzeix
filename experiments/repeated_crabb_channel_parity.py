#!/usr/bin/env python3
"""Audit transfer-channel parity in the doubled-Hardy pencil.

Changing the Hankel channel ``H`` to ``-H`` conjugates the complete
physical pencil by ``diag(I,-I)``.  The checker propagates this exact
symmetry through functional calculus, resolvents, metric defects,
endpoint embeddings, and orientation Schur complements.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float


Matrix = np.ndarray
PARAMETER = 0.23


@dataclass(frozen=True)
class ChannelParityRecord:
    """One doubled-Hardy sign-covariance audit."""

    hardy_rows: int
    multiplicity: int
    pencil_covariance_error: str
    polynomial_covariance_error: str
    resolvent_covariance_error: str
    diagonal_even_error: str
    off_diagonal_odd_error: str
    right_endpoint_covariance_error: str
    left_endpoint_covariance_error: str
    right_sandwich_even_error: str
    left_sandwich_even_error: str
    right_schur_even_error: str
    left_schur_even_error: str
    all_checks_passed: bool


def block_backward_shift(
    row_count: int,
    multiplicity: int,
) -> Matrix:
    """Return the backward shift on a finite block Hardy window."""

    dimension = row_count * multiplicity
    shift = np.zeros((dimension, dimension), dtype=complex)
    identity = np.eye(multiplicity, dtype=complex)
    for row in range(row_count - 1):
        start = row * multiplicity
        shift[
            start : start + multiplicity,
            start + multiplicity : start + 2 * multiplicity,
        ] = identity
    return shift


def row_zero_inclusion(
    row_count: int,
    multiplicity: int,
) -> Matrix:
    """Return the row-zero Hardy inclusion."""

    inclusion = np.zeros(
        (row_count * multiplicity, multiplicity),
        dtype=complex,
    )
    inclusion[:multiplicity] = np.eye(multiplicity)
    return inclusion


def physical_pencil(
    hankel: Matrix,
    backward: Matrix,
    first: Matrix,
) -> Matrix:
    """Return L255's doubled physical ellipse pencil."""

    dimension = len(backward)
    identity = np.eye(dimension, dtype=complex)
    forward = backward.conj().T
    right_right = backward + PARAMETER * (
        forward @ (identity + first)
        + hankel
        @ first
        @ backward
        @ hankel.conj().T
        @ first
    )
    right_left = PARAMETER * (
        -forward @ hankel @ first
        + hankel @ first @ backward
    )
    left_right = (
        -forward @ hankel.conj().T @ first
        + PARAMETER
        * (identity + first)
        @ backward
        @ hankel.conj().T
        @ first
    )
    left_left = (
        forward
        + PARAMETER * (identity + first) @ backward
    )
    return np.block(
        [
            [right_right, right_left],
            [left_right, left_left],
        ]
    )


def odd_polynomial(matrix: Matrix) -> Matrix:
    """Return a nontrivial odd polynomial functional calculus."""

    square = matrix @ matrix
    cube = square @ matrix
    fifth = cube @ square
    return matrix - 0.17 * cube + 0.031 * fifth


def orientation_schur(
    matrix: Matrix,
    dimension: int,
    keep_right: bool,
) -> Matrix:
    """Schur-compress one whole Hardy orientation."""

    right_right = matrix[:dimension, :dimension]
    right_left = matrix[:dimension, dimension:]
    left_right = matrix[dimension:, :dimension]
    left_left = matrix[dimension:, dimension:]
    if keep_right:
        return (
            right_right
            - right_left @ np.linalg.solve(left_left, left_right)
        )
    return (
        left_left
        - left_right @ np.linalg.solve(right_right, right_left)
    )


def audit_case(
    row_count: int,
    multiplicity: int,
    seed: int,
) -> ChannelParityRecord:
    """Audit one random noncommuting channel."""

    generator = np.random.default_rng(seed)
    dimension = row_count * multiplicity
    backward = block_backward_shift(row_count, multiplicity)
    inclusion = row_zero_inclusion(row_count, multiplicity)
    first = inclusion @ inclusion.conj().T
    hankel = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    hankel *= 0.11 / max(1.0, np.linalg.norm(hankel))
    hankel[:multiplicity, :multiplicity] = 0

    pencil = physical_pencil(hankel, backward, first)
    negative_pencil = physical_pencil(-hankel, backward, first)
    identity = np.eye(dimension, dtype=complex)
    sign = np.block(
        [
            [identity, np.zeros_like(identity)],
            [np.zeros_like(identity), -identity],
        ]
    )
    pencil_error = float(
        np.linalg.norm(negative_pencil - sign @ pencil @ sign)
    )

    polynomial = odd_polynomial(pencil)
    negative_polynomial = odd_polynomial(negative_pencil)
    polynomial_error = float(
        np.linalg.norm(
            negative_polynomial - sign @ polynomial @ sign
        )
    )

    spectral_parameter = 1.7 + 0.2j
    doubled_identity = np.eye(2 * dimension, dtype=complex)
    resolvent = np.linalg.inv(
        spectral_parameter * doubled_identity - pencil
    )
    negative_resolvent = np.linalg.inv(
        spectral_parameter * doubled_identity - negative_pencil
    )
    resolvent_error = float(
        np.linalg.norm(
            negative_resolvent - sign @ resolvent @ sign
        )
    )

    right_weights = np.repeat(
        np.linspace(1.0, 0.73, row_count),
        multiplicity,
    )
    left_weights = np.repeat(
        np.linspace(0.0, 0.41, row_count),
        multiplicity,
    )
    metric = np.diag(np.concatenate((right_weights, left_weights)))
    closed = (
        2.3 * doubled_identity
        + metric
        + polynomial.conj().T @ metric @ polynomial
        + 0.07 * resolvent.conj().T @ resolvent
    )
    negative_closed = (
        2.3 * doubled_identity
        + metric
        + negative_polynomial.conj().T
        @ metric
        @ negative_polynomial
        + 0.07
        * negative_resolvent.conj().T
        @ negative_resolvent
    )
    covariance = negative_closed - sign @ closed @ sign
    diagonal_error = max(
        float(np.linalg.norm(covariance[:dimension, :dimension])),
        float(np.linalg.norm(covariance[dimension:, dimension:])),
    )
    off_diagonal_error = max(
        float(
            np.linalg.norm(
                negative_closed[:dimension, dimension:]
                + closed[:dimension, dimension:]
            )
        ),
        float(
            np.linalg.norm(
                negative_closed[dimension:, :dimension]
                + closed[dimension:, :dimension]
            )
        ),
    )

    right_endpoint = np.vstack(
        (inclusion, hankel.conj().T @ inclusion)
    )
    negative_right_endpoint = np.vstack(
        (inclusion, -hankel.conj().T @ inclusion)
    )
    left_endpoint = np.vstack(
        (hankel @ inclusion, inclusion)
    )
    negative_left_endpoint = np.vstack(
        (-hankel @ inclusion, inclusion)
    )
    right_endpoint_error = float(
        np.linalg.norm(
            negative_right_endpoint - sign @ right_endpoint
        )
    )
    left_endpoint_error = float(
        np.linalg.norm(
            negative_left_endpoint + sign @ left_endpoint
        )
    )
    right_sandwich = right_endpoint.conj().T @ closed @ right_endpoint
    negative_right_sandwich = (
        negative_right_endpoint.conj().T
        @ negative_closed
        @ negative_right_endpoint
    )
    left_sandwich = left_endpoint.conj().T @ closed @ left_endpoint
    negative_left_sandwich = (
        negative_left_endpoint.conj().T
        @ negative_closed
        @ negative_left_endpoint
    )
    right_sandwich_error = float(
        np.linalg.norm(negative_right_sandwich - right_sandwich)
    )
    left_sandwich_error = float(
        np.linalg.norm(negative_left_sandwich - left_sandwich)
    )

    right_schur = orientation_schur(
        closed,
        dimension,
        keep_right=True,
    )
    negative_right_schur = orientation_schur(
        negative_closed,
        dimension,
        keep_right=True,
    )
    left_schur = orientation_schur(
        closed,
        dimension,
        keep_right=False,
    )
    negative_left_schur = orientation_schur(
        negative_closed,
        dimension,
        keep_right=False,
    )
    right_schur_error = float(
        np.linalg.norm(negative_right_schur - right_schur)
    )
    left_schur_error = float(
        np.linalg.norm(negative_left_schur - left_schur)
    )

    tolerance = 3e-10
    verified = bool(
        pencil_error < tolerance
        and polynomial_error < tolerance
        and resolvent_error < tolerance
        and diagonal_error < tolerance
        and off_diagonal_error < tolerance
        and right_endpoint_error < tolerance
        and left_endpoint_error < tolerance
        and right_sandwich_error < tolerance
        and left_sandwich_error < tolerance
        and right_schur_error < tolerance
        and left_schur_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "channel-parity audit failed: "
            f"{row_count=}, {multiplicity=}, {pencil_error=}, "
            f"{polynomial_error=}, {resolvent_error=}, "
            f"{diagonal_error=}, {off_diagonal_error=}, "
            f"{right_endpoint_error=}, {left_endpoint_error=}, "
            f"{right_sandwich_error=}, {left_sandwich_error=}, "
            f"{right_schur_error=}, {left_schur_error=}"
        )
    return ChannelParityRecord(
        hardy_rows=row_count,
        multiplicity=multiplicity,
        pencil_covariance_error=format_float(pencil_error),
        polynomial_covariance_error=format_float(polynomial_error),
        resolvent_covariance_error=format_float(resolvent_error),
        diagonal_even_error=format_float(diagonal_error),
        off_diagonal_odd_error=format_float(off_diagonal_error),
        right_endpoint_covariance_error=format_float(
            right_endpoint_error
        ),
        left_endpoint_covariance_error=format_float(
            left_endpoint_error
        ),
        right_sandwich_even_error=format_float(
            right_sandwich_error
        ),
        left_sandwich_even_error=format_float(left_sandwich_error),
        right_schur_even_error=format_float(right_schur_error),
        left_schur_even_error=format_float(left_schur_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[ChannelParityRecord]:
    """Return deterministic parity audits."""

    return [
        audit_case(
            row_count,
            multiplicity,
            145_000 + 100 * multiplicity + row_count,
        )
        for multiplicity in (1, 2, 3)
        for row_count in (5, 6, 7, 8)
    ]


def write_records(
    records: list[ChannelParityRecord],
    output: Path,
) -> str:
    """Write sorted JSON lines and return the SHA-256 digest."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = "\n".join(
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ) + "\n"
    output.write_text(payload)
    return hashlib.sha256(payload.encode()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_channel_parity_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    print(
        json.dumps(
            {
                "records": len(records),
                "sha256": digest,
                "maximum_pencil_covariance_error": max(
                    float(record.pencil_covariance_error)
                    for record in records
                ),
                "maximum_endpoint_even_error": max(
                    max(
                        float(record.right_sandwich_even_error),
                        float(record.left_sandwich_even_error),
                    )
                    for record in records
                ),
                "maximum_schur_even_error": max(
                    max(
                        float(record.right_schur_even_error),
                        float(record.left_schur_even_error),
                    )
                    for record in records
                ),
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
