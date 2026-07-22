#!/usr/bin/env python3
"""Stress L21 near repeated Crabb equality blocks.

The completely bounded target ``t_*(phi(A)) <= 4`` is sharp on Crabb disk
blocks.  Direct sums create a multiplicity locus where the numerical-range
support eigenvalue is repeated.  This probe applies small deterministic
perturbations in three transverse directions:

* full dense perturbations;
* couplings between rotated copies of the same Crabb block;
* noncommuting operator weights and diagonal blocks in level-major form.

Every case is normalized by its numerical range and passed through the full
map-resolution, Cauchy, double-layer, and primal/dual SDP gates in
``general_similarity_sdp``.  A dense support-eigenvalue-gap diagnostic records
whether the perturbed boundary has separated from the multiplicity locus.
The results are numerical falsification evidence, never a proof or a certified
counterexample.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Iterator

import numpy as np

from crouzeix import crabb_matrix
from general_similarity_sdp import evaluate_case, normalize_numerical_range


PERTURBATION_FAMILIES = ("full", "cross", "operator_weight")
COPY_PHASE_STEP = 0.37


def normalized_direction(matrix: np.ndarray) -> np.ndarray:
    """Remove the scalar part and normalize one nonzero perturbation."""

    dimension = matrix.shape[0]
    centered = matrix - np.trace(matrix) * np.eye(dimension) / dimension
    norm = float(np.linalg.norm(centered, 2))
    if norm <= 1e-14:
        raise ValueError("perturbation direction is numerically zero")
    return centered / norm


def random_complex(
    generator: np.random.Generator,
    shape: tuple[int, int],
) -> np.ndarray:
    """Return a standard complex Gaussian matrix."""

    return generator.standard_normal(shape) + 1j * generator.standard_normal(shape)


def block_diagonal(blocks: list[np.ndarray]) -> np.ndarray:
    """Assemble a dense block diagonal matrix without SciPy."""

    dimension = sum(block.shape[0] for block in blocks)
    output = np.zeros((dimension, dimension), dtype=complex)
    start = 0
    for block in blocks:
        stop = start + block.shape[0]
        output[start:stop, start:stop] = block
        start = stop
    return output


def copy_major_base(block_size: int, multiplicity: int) -> np.ndarray:
    """Return rotated copies of a Crabb equality block."""

    crabb = crabb_matrix(block_size - 1)
    phases = np.exp(COPY_PHASE_STEP * 1j * np.arange(multiplicity))
    return block_diagonal([phase * crabb for phase in phases])


def cross_direction(
    block_size: int,
    multiplicity: int,
    generator: np.random.Generator,
) -> np.ndarray:
    """Couple distinct copy-major Crabb blocks in both directions."""

    dimension = block_size * multiplicity
    output = np.zeros((dimension, dimension), dtype=complex)
    for left in range(multiplicity):
        for right in range(left + 1, multiplicity):
            left_slice = slice(left * block_size, (left + 1) * block_size)
            right_slice = slice(right * block_size, (right + 1) * block_size)
            output[left_slice, right_slice] = random_complex(
                generator,
                (block_size, block_size),
            )
            output[right_slice, left_slice] = random_complex(
                generator,
                (block_size, block_size),
            )
    return normalized_direction(output)


def operator_weight_case(
    block_size: int,
    multiplicity: int,
    generator: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a level-major equality block and a noncommuting perturbation."""

    crabb = crabb_matrix(block_size - 1)
    identity = np.eye(multiplicity)
    base = np.kron(crabb, identity)
    dimension = block_size * multiplicity
    perturbation = np.zeros((dimension, dimension), dtype=complex)

    for level in range(block_size):
        level_slice = slice(level * multiplicity, (level + 1) * multiplicity)
        perturbation[level_slice, level_slice] = random_complex(
            generator,
            (multiplicity, multiplicity),
        )
    for level in range(block_size - 1):
        left = slice(level * multiplicity, (level + 1) * multiplicity)
        right = slice((level + 1) * multiplicity, (level + 2) * multiplicity)
        perturbation[left, right] = random_complex(
            generator,
            (multiplicity, multiplicity),
        )
        perturbation[right, left] = random_complex(
            generator,
            (multiplicity, multiplicity),
        )
    return base, normalized_direction(perturbation)


def generated_directions(
    block_sizes: list[int],
    multiplicities: list[int],
    families: list[str],
    directions: int,
    seed: int,
) -> Iterator[tuple[int, int, str, int, np.ndarray, np.ndarray]]:
    """Yield deterministic equality bases and normalized perturbations."""

    for block_size in block_sizes:
        for multiplicity in multiplicities:
            for family in families:
                for direction in range(directions):
                    family_index = PERTURBATION_FAMILIES.index(family)
                    generator = np.random.default_rng(
                        np.random.SeedSequence(
                            [seed, block_size, multiplicity, family_index, direction]
                        )
                    )
                    if family == "operator_weight":
                        base, perturbation = operator_weight_case(
                            block_size,
                            multiplicity,
                            generator,
                        )
                    else:
                        base = copy_major_base(block_size, multiplicity)
                        dimension = base.shape[0]
                        if family == "full":
                            perturbation = normalized_direction(
                                random_complex(generator, (dimension, dimension))
                            )
                        elif family == "cross":
                            perturbation = cross_direction(
                                block_size,
                                multiplicity,
                                generator,
                            )
                        else:
                            raise ValueError(f"unknown family: {family}")
                    yield (
                        block_size,
                        multiplicity,
                        family,
                        direction,
                        base,
                        perturbation,
                    )


def support_eigenvalue_gap(matrix: np.ndarray, resolution: int) -> float:
    """Return the sampled minimum gap between the two top support eigenvalues."""

    if matrix.shape[0] < 2:
        return float("inf")
    minimum = float("inf")
    angles = np.linspace(0.0, 2.0 * np.pi, resolution, endpoint=False)
    adjoint = matrix.conj().T
    for angle in angles:
        phase = np.exp(-1j * angle)
        support = (phase * matrix + np.conj(phase) * adjoint) / 2
        eigenvalues = np.linalg.eigvalsh(support)
        minimum = min(minimum, float(eigenvalues[-1] - eigenvalues[-2]))
    return minimum


def matrix_sha256(matrix: np.ndarray) -> str:
    """Return a stable fingerprint of a contiguous complex128 matrix."""

    canonical = np.ascontiguousarray(matrix, dtype=np.complex128)
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def evaluate_normalized_matrix(
    normalized: np.ndarray,
    block_size: int,
    multiplicity: int,
    family: str,
    direction: int,
    delta: float,
    inflate: float,
    seed: int,
    resolution: int,
    maximum_resolution: int,
    support_resolution: int,
    support_gap: float,
    fingerprint: str,
) -> dict[str, object]:
    """Run all general-map gates and attach equality-locus diagnostics."""

    label = f"crabb{block_size}_m{multiplicity}_{family}"
    record = evaluate_case(
        normalized.shape[0],
        label,
        direction,
        normalized,
        seed,
        inflate,
        resolution,
        maximum_resolution,
    )
    output: dict[str, object] = asdict(record)
    output.update(
        {
            "block_size": block_size,
            "multiplicity": multiplicity,
            "perturbation_family": family,
            "direction": direction,
            "delta": delta,
            "matrix_sha256": fingerprint,
            "support_gap": support_gap,
            "support_resolution": support_resolution,
        }
    )
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-sizes", nargs="+", type=int, default=[3, 4])
    parser.add_argument("--multiplicities", nargs="+", type=int, default=[2, 3])
    parser.add_argument(
        "--families",
        nargs="+",
        choices=PERTURBATION_FAMILIES,
        default=list(PERTURBATION_FAMILIES),
    )
    parser.add_argument("--directions", type=int, default=1)
    parser.add_argument(
        "--deltas", nargs="+", type=float, default=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
    )
    parser.add_argument("--inflates", nargs="+", type=float, default=[0.0, 0.000625])
    parser.add_argument("--seed", type=int, default=9173401)
    parser.add_argument(
        "--include-base",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="include one exact repeated-Crabb calibration per size",
    )
    parser.add_argument("--resolution", type=int, default=256)
    parser.add_argument("--max-resolution", type=int, default=2048)
    parser.add_argument("--support-resolution", type=int, default=4096)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.directions < 1:
        raise ValueError("directions must be positive")
    if any(block_size < 2 for block_size in args.block_sizes):
        raise ValueError("block sizes must be at least two")
    if any(multiplicity < 2 for multiplicity in args.multiplicities):
        raise ValueError("multiplicities must be at least two")
    if any(delta <= 0 for delta in args.deltas):
        raise ValueError("deltas must be positive")
    if any(inflate < 0 for inflate in args.inflates):
        raise ValueError("inflates must be nonnegative")
    if args.support_resolution < 64:
        raise ValueError("support resolution must be at least 64")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_stream = args.output.open("w")

    def record_case(
        matrix: np.ndarray,
        block_size: int,
        multiplicity: int,
        family: str,
        direction: int,
        delta: float,
    ) -> None:
        normalized = normalize_numerical_range(matrix)
        support_gap = support_eigenvalue_gap(normalized, args.support_resolution)
        fingerprint = matrix_sha256(normalized)
        for inflate in args.inflates:
            record = evaluate_normalized_matrix(
                normalized,
                block_size,
                multiplicity,
                family,
                direction,
                delta,
                inflate,
                args.seed,
                args.resolution,
                args.max_resolution,
                args.support_resolution,
                support_gap,
                fingerprint,
            )
            line = json.dumps(record, allow_nan=True, sort_keys=True)
            print(line, flush=True)
            output_stream.write(line + "\n")
            output_stream.flush()

    try:
        if args.include_base:
            for block_size in args.block_sizes:
                for multiplicity in args.multiplicities:
                    base = np.kron(crabb_matrix(block_size - 1), np.eye(multiplicity))
                    record_case(base, block_size, multiplicity, "base", 0, 0.0)

        for (
            block_size,
            multiplicity,
            family,
            direction,
            base,
            perturbation,
        ) in generated_directions(
            args.block_sizes,
            args.multiplicities,
            args.families,
            args.directions,
            args.seed,
        ):
            for delta in args.deltas:
                record_case(
                    base + delta * perturbation,
                    block_size,
                    multiplicity,
                    family,
                    direction,
                    delta,
                )
    finally:
        output_stream.close()


if __name__ == "__main__":
    main()
