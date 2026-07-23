#!/usr/bin/env python3
"""Audit the best-phase normal form of the palindromic equality cone."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class NormalFormRecord:
    dimension: int
    samples: int
    maximum_involution_error: float
    maximum_reconstruction_error: float
    maximum_real_orthogonality_error: float
    maximum_norm_identity_error: float
    maximum_quartic_identity_error: float
    maximum_distance_identity_error: float
    maximum_equality_quartic: float


def reversal_conjugate(vector: np.ndarray) -> np.ndarray:
    """Apply ``J conjugate(.)``."""

    return vector.conjugate()[::-1]


def best_phase_split(
    vector: np.ndarray,
) -> tuple[complex, np.ndarray, np.ndarray]:
    """Return the phase, equality component, and normal component."""

    pairing = vector @ vector[::-1]
    phase = pairing / abs(pairing) if pairing != 0 else 1.0 + 0.0j
    reflected = phase * reversal_conjugate(vector)
    equality_component = (vector + reflected) / 2
    normal_component = (vector - reflected) / 2
    return phase, equality_component, normal_component


def audit_dimension(
    dimension: int,
    samples: int,
    generator: np.random.Generator,
) -> NormalFormRecord:
    """Audit all identities on random and exact-equality vectors."""

    error_names = (
        "involution",
        "reconstruction",
        "orthogonality",
        "norm",
        "quartic",
        "distance",
    )
    errors = {name: 0.0 for name in error_names}
    maximum_equality_quartic = 0.0

    vectors: list[np.ndarray] = []
    for _ in range(samples):
        vector = (
            generator.standard_normal(dimension)
            + 1j * generator.standard_normal(dimension)
        )
        vectors.append(vector)

        phase = np.exp(1j * generator.uniform(-np.pi, np.pi))
        seed = (
            generator.standard_normal(dimension)
            + 1j * generator.standard_normal(dimension)
        )
        equality_vector = (
            seed + phase * reversal_conjugate(seed)
        ) / 2
        vectors.append(equality_vector)

    for index, vector in enumerate(vectors):
        phase, equality_component, normal_component = best_phase_split(
            vector
        )
        involution_equality = (
            phase * reversal_conjugate(equality_component)
        )
        involution_normal = (
            phase * reversal_conjugate(normal_component)
        )
        errors["involution"] = max(
            errors["involution"],
            float(
                np.linalg.norm(
                    involution_equality - equality_component
                )
            ),
            float(
                np.linalg.norm(
                    involution_normal + normal_component
                )
            ),
        )
        errors["reconstruction"] = max(
            errors["reconstruction"],
            float(
                np.linalg.norm(
                    vector - equality_component - normal_component
                )
            ),
        )
        errors["orthogonality"] = max(
            errors["orthogonality"],
            abs(
                float(
                    np.vdot(
                        equality_component,
                        normal_component,
                    ).real
                )
            ),
        )

        norm_square = float(np.vdot(vector, vector).real)
        pairing = vector @ vector[::-1]
        equality_norm_square = float(
            np.vdot(equality_component, equality_component).real
        )
        normal_norm_square = float(
            np.vdot(normal_component, normal_component).real
        )
        errors["norm"] = max(
            errors["norm"],
            abs(
                equality_norm_square
                - (norm_square + abs(pairing)) / 2
            ),
            abs(
                normal_norm_square
                - (norm_square - abs(pairing)) / 2
            ),
        )

        quartic = norm_square**2 - abs(pairing) ** 2
        errors["quartic"] = max(
            errors["quartic"],
            abs(
                quartic
                - 4 * equality_norm_square * normal_norm_square
            ),
        )
        errors["distance"] = max(
            errors["distance"],
            abs(
                normal_norm_square
                - (norm_square - abs(pairing)) / 2
            ),
        )
        if index % 2:
            maximum_equality_quartic = max(
                maximum_equality_quartic,
                abs(quartic),
            )

    tolerance = 2e-10 * max(1, dimension)
    if any(error > tolerance for error in errors.values()):
        raise RuntimeError(
            f"normal-form audit failed in dimension {dimension}: {errors}"
        )
    if maximum_equality_quartic > tolerance:
        raise RuntimeError(
            "the constructed equality vector left the quartic null"
        )

    return NormalFormRecord(
        dimension=dimension,
        samples=samples,
        maximum_involution_error=errors["involution"],
        maximum_reconstruction_error=errors["reconstruction"],
        maximum_real_orthogonality_error=errors["orthogonality"],
        maximum_norm_identity_error=errors["norm"],
        maximum_quartic_identity_error=errors["quartic"],
        maximum_distance_identity_error=errors["distance"],
        maximum_equality_quartic=maximum_equality_quartic,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-dimension", type=int, default=20)
    parser.add_argument("--samples", type=int, default=20)
    parser.add_argument("--seed", type=int, default=70223)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the numerical audits and optionally persist JSONL output."""

    args = parse_args()
    if args.maximum_dimension < 1 or args.samples < 1:
        raise ValueError("dimensions and samples must be positive")
    generator = np.random.default_rng(args.seed)
    records = [
        audit_dimension(dimension, args.samples, generator)
        for dimension in range(1, args.maximum_dimension + 1)
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
