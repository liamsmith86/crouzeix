#!/usr/bin/env python3
"""Audit the two lossless Redheffer defect identities from L267."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class RedhefferDefectRecord:
    """One deterministic noncommutative feedback audit."""

    external_dimension: int
    state_dimension: int
    maximum_initial_error: str
    maximum_final_error: str
    all_checks_passed: bool


def random_unitary(
    dimension: int,
    generator: np.random.Generator,
) -> np.ndarray:
    """Return a deterministic Haar-style unitary from a complex QR factor."""

    matrix = generator.normal(size=(dimension, dimension))
    matrix = matrix + 1j * generator.normal(size=matrix.shape)
    unitary, triangular = np.linalg.qr(matrix)
    phases = np.diag(triangular)
    return unitary * (phases / np.abs(phases)).conj()


def audit_dimensions(
    external_dimension: int,
    state_dimension: int,
) -> RedhefferDefectRecord:
    """Audit both identities on twenty noncommuting feedback operators."""

    generator = np.random.default_rng(
        267_000 + 100 * external_dimension + state_dimension,
    )
    maximum_initial_error = 0.0
    maximum_final_error = 0.0
    for _ in range(20):
        total_dimension = external_dimension + state_dimension
        unitary = random_unitary(total_dimension, generator)
        split = external_dimension
        a = unitary[:split, :split]
        b = unitary[:split, split:]
        c = unitary[split:, :split]
        d = unitary[split:, split:]

        feedback = generator.normal(size=(state_dimension, state_dimension))
        feedback = feedback + 1j * generator.normal(size=feedback.shape)
        feedback *= 0.2 / np.linalg.norm(feedback)

        identity = np.eye(state_dimension, dtype=complex)
        right_inverse = np.linalg.inv(identity - d @ feedback)
        phi = a + b @ feedback @ right_inverse @ c

        initial = (
            c.conj().T
            @ np.linalg.inv(identity - feedback.conj().T @ d.conj().T)
            @ (identity - feedback.conj().T @ feedback)
            @ right_inverse
            @ c
        )
        final = (
            b
            @ np.linalg.inv(identity - feedback @ d)
            @ (identity - feedback @ feedback.conj().T)
            @ np.linalg.inv(identity - d.conj().T @ feedback.conj().T)
            @ b.conj().T
        )
        external_identity = np.eye(external_dimension, dtype=complex)
        maximum_initial_error = max(
            maximum_initial_error,
            float(np.linalg.norm(external_identity - phi.conj().T @ phi - initial)),
        )
        maximum_final_error = max(
            maximum_final_error,
            float(np.linalg.norm(external_identity - phi @ phi.conj().T - final)),
        )

    verified = max(maximum_initial_error, maximum_final_error) < 2e-14
    if not verified:
        raise RuntimeError("the lossless Redheffer defect audit failed")
    return RedhefferDefectRecord(
        external_dimension=external_dimension,
        state_dimension=state_dimension,
        maximum_initial_error=f"{maximum_initial_error:.3e}",
        maximum_final_error=f"{maximum_final_error:.3e}",
        all_checks_passed=verified,
    )


def write_records(
    records: list[RedhefferDefectRecord],
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
            "experiments/lossless_redheffer_defect_transport_s70225.jsonl",
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic audit."""

    args = parse_args()
    records = [
        audit_dimensions(external_dimension, state_dimension)
        for external_dimension, state_dimension in ((1, 2), (2, 3), (3, 2))
    ]
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
