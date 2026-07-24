#!/usr/bin/env python3
"""Audit the repeated-Crabb circular-normal Jensen reduction.

On a maximal common-winner copy space, Fourier independence forces every
true circular-normal coefficient to be scalar.  The scalar L188 response
can then be tested on each copy-space vector.  This checker verifies the
block reflection identity, the resulting statewise absorption, sharp
flux-limit examples, and compatibility with a prescribed residual kernel.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_inverse_gram_kernel import (
    deterministic_copy_unitary,
)


@dataclass(frozen=True)
class CircularJensenRecord:
    """One Fourier-scalarization or residual-absorption audit."""

    record_kind: str
    length: int
    multiplicity: int
    winner_dimension: int
    active_mode_count: int
    reflection_error: str
    maximum_absorption_ratio: str
    minimum_absorption_margin: str
    winner_scalarization_error: str
    winner_cross_error: str
    residual_kernel_error: str
    response_kernel_error: str
    normal_kernel_curvature: str
    all_checks_passed: bool


def normal_curvature(length: int, mode: int) -> float:
    """Return L173's positive curvature in one complex normal mode."""

    remainder = length - mode
    flux = (
        (4 * mode - 1) ** 2
        * remainder
        * (remainder - 1)
        * (remainder - 2)
        / 24
    )
    null_lift = (
        2
        * mode
        * (mode - 1)
        * (mode - 2)
        * (remainder + 1 / 4) ** 2
        / 3
    )
    return (flux + null_lift) / length**4


def response_blocks(
    residual: np.ndarray,
    length: int,
) -> list[np.ndarray]:
    """Return L188's complex circular-normal response matrices."""

    responses = []
    for mode in range(3, length - 2):
        projection = np.zeros_like(residual[0, 0])
        for left in range(length - 1):
            for right in range(left + 1, length - 1):
                if (left + 1) + (right + 1) != length - mode:
                    continue
                skew = (
                    residual[left, right]
                    - residual[right, left]
                ) / 2
                projection += (right - left) * skew
        scale = 4 * (4 * mode - 1) / length**2
        responses.append(scale * projection)
    return responses


def residual_gram(residual: np.ndarray) -> np.ndarray:
    """Return the sum of all right block Gramians."""

    multiplicity = residual[0, 0].shape[0]
    gram = np.zeros((multiplicity, multiplicity), dtype=complex)
    for block in residual.flat:
        gram += block.conj().T @ block
    return (gram + gram.conj().T) / 2


def reflection_error(residual: np.ndarray) -> float:
    """Return the error in ``F_ab = F_(Jb,Ja)^*``."""

    size = residual.shape[0]
    worst = 0.0
    for left in range(size):
        for right in range(size):
            reflected = residual[size - 1 - right, size - 1 - left]
            worst = max(
                worst,
                float(
                    np.linalg.norm(
                        residual[left, right] - reflected.conj().T
                    )
                ),
            )
    return worst


def reflected_residual(
    length: int,
    multiplicity: int,
    seed: int,
    active: np.ndarray | None = None,
) -> np.ndarray:
    """Return a deterministic block residual with Hardy reflection."""

    rng = np.random.default_rng(seed)
    size = length - 1
    residual = np.empty((size, size), dtype=object)
    completed: set[tuple[int, int]] = set()
    for left in range(size):
        for right in range(size):
            if (left, right) in completed:
                continue
            reflected = (size - 1 - right, size - 1 - left)
            block = (
                rng.standard_normal((multiplicity, multiplicity))
                + 1j
                * rng.standard_normal((multiplicity, multiplicity))
            )
            if active is not None:
                block = active @ block @ active
            if reflected == (left, right):
                block = (block + block.conj().T) / 2
            residual[left, right] = block
            residual[reflected] = block.conj().T
            completed.add((left, right))
            completed.add(reflected)
    return residual


def sharp_residual(
    length: int,
    multiplicity: int,
    mode: int,
) -> np.ndarray:
    """Return a flux-Cauchy equality residual in one active mode."""

    size = length - 1
    identity = np.eye(multiplicity, dtype=complex)
    residual = np.empty((size, size), dtype=object)
    for index in np.ndindex(size, size):
        residual[index] = np.zeros_like(identity)
    for left in range(size):
        for right in range(left + 1, size):
            if (left + 1) + (right + 1) != length - mode:
                continue
            weight = right - left
            residual[left, right] = weight * identity
            residual[right, left] = -weight * identity
            reflected = (size - 1 - right, size - 1 - left)
            residual[reflected] = weight * identity
            residual[reflected[::-1]] = -weight * identity
    return residual


def deterministic_states(
    multiplicity: int,
    count: int,
    seed: int,
) -> list[np.ndarray]:
    """Return coordinate and dense deterministic unit copy vectors."""

    states = [
        np.eye(multiplicity, dtype=complex)[:, index]
        for index in range(multiplicity)
    ]
    rng = np.random.default_rng(seed)
    for _ in range(count):
        vector = (
            rng.standard_normal(multiplicity)
            + 1j * rng.standard_normal(multiplicity)
        )
        states.append(vector / np.linalg.norm(vector))
    return states


def absorption_statistics(
    residual: np.ndarray,
    length: int,
    states: Sequence[np.ndarray],
) -> tuple[float, float]:
    """Return the largest completed-gain ratio and smallest margin."""

    gram = residual_gram(residual)
    responses = response_blocks(residual, length)
    curvatures = [
        normal_curvature(length, mode)
        for mode in range(3, length - 2)
    ]
    maximum_ratio = 0.0
    minimum_margin = float("inf")
    for vector in states:
        base = 4 * float(np.real(vector.conj() @ gram @ vector))
        gain = sum(
            abs(vector.conj() @ response @ vector) ** 2
            / (4 * curvature)
            for response, curvature in zip(responses, curvatures)
        )
        ratio = gain / base if base > 0 else 0.0
        maximum_ratio = max(maximum_ratio, float(ratio))
        minimum_margin = min(minimum_margin, float(base - gain))
    return maximum_ratio, minimum_margin


def fourier_scalarization_errors(
    length: int,
    multiplicity: int,
    winner_dimension: int,
) -> tuple[float, float]:
    """Audit coefficientwise scalarization on a common-winner block."""

    unitary = deterministic_copy_unitary(multiplicity)
    winner = unitary[:, :winner_dimension]
    loser = unitary[:, winner_dimension:]
    scalar_error = 0.0
    cross_error = 0.0
    for mode in range(3, length + 1):
        for polarization in range(2):
            scalar = (2 * mode + polarization + 1) / (7 * length)
            coefficient = scalar * np.eye(multiplicity, dtype=complex)
            loser_dimension = loser.shape[1]
            if loser_dimension:
                seed = np.zeros(
                    (loser_dimension, loser_dimension),
                    dtype=complex,
                )
                for row in range(loser_dimension):
                    seed[row, row] = (1 + row) / (11 * length)
                coefficient += loser @ seed @ loser.conj().T
            winner_block = winner.conj().T @ coefficient @ winner
            cross = loser.conj().T @ coefficient @ winner
            scalar_error = max(
                scalar_error,
                float(
                    np.linalg.norm(
                        winner_block
                        - np.trace(winner_block)
                        / winner_dimension
                        * np.eye(winner_dimension)
                    )
                ),
            )
            cross_error = max(cross_error, float(np.linalg.norm(cross)))
    return scalar_error, cross_error


def make_record(
    kind: str,
    length: int,
    multiplicity: int,
    winner_dimension: int,
    seed: int,
) -> CircularJensenRecord:
    """Build one scalarization, absorption, and kernel audit record."""

    unitary = deterministic_copy_unitary(multiplicity)
    kernel = unitary[:, :winner_dimension]
    kernel_projection = kernel @ kernel.conj().T
    active_projection = np.eye(multiplicity) - kernel_projection
    if kind == "sharp_flux":
        residual = sharp_residual(length, multiplicity, mode=3)
    elif kind == "kernel":
        residual = reflected_residual(
            length,
            multiplicity,
            seed,
            active_projection,
        )
    else:
        residual = reflected_residual(
            length,
            multiplicity,
            seed,
        )

    states = deterministic_states(multiplicity, 64, seed + 1)
    ratio, margin = absorption_statistics(residual, length, states)
    gram = residual_gram(residual)
    responses = response_blocks(residual, length)
    residual_kernel = float(np.linalg.norm(gram @ kernel))
    response_kernel = max(
        (
            float(
                np.linalg.norm(response @ kernel)
                + np.linalg.norm(response.conj().T @ kernel)
            )
            for response in responses
        ),
        default=0.0,
    )
    scalar_error, cross_error = fourier_scalarization_errors(
        length,
        multiplicity,
        winner_dimension,
    )
    curvatures = [
        normal_curvature(length, mode)
        for mode in range(3, length - 2)
    ]
    normal_kernel_curvature = min(curvatures, default=0.0)

    expected_kernel = kind == "kernel"
    verified = bool(
        reflection_error(residual) < 1e-12
        and ratio < 1 - 1e-9
        and margin > 1e-10
        and scalar_error < 1e-12
        and cross_error < 1e-12
        and (
            not expected_kernel
            or (
                residual_kernel < 1e-10
                and response_kernel < 1e-10
            )
        )
        and normal_kernel_curvature > 0
    )
    if not verified:
        raise RuntimeError(
            "the repeated circular-normal Jensen audit failed "
            f"for {kind}, L={length}, m={multiplicity}"
        )
    return CircularJensenRecord(
        record_kind=kind,
        length=length,
        multiplicity=multiplicity,
        winner_dimension=winner_dimension,
        active_mode_count=max(0, length - 5),
        reflection_error=format_float(reflection_error(residual)),
        maximum_absorption_ratio=format_float(ratio),
        minimum_absorption_margin=format_float(margin),
        winner_scalarization_error=format_float(scalar_error),
        winner_cross_error=format_float(cross_error),
        residual_kernel_error=format_float(residual_kernel),
        response_kernel_error=format_float(response_kernel),
        normal_kernel_curvature=format_float(normal_kernel_curvature),
        all_checks_passed=verified,
    )


def write_records(
    path: Path,
    records: Sequence[CircularJensenRecord],
) -> None:
    """Write deterministic JSON Lines output atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Run the repeated circular-normal Jensen audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=6)
    parser.add_argument("--maximum-length", type=int, default=14)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(2, 3, 4),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_circular_jensen_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records: list[CircularJensenRecord] = []
    for multiplicity in args.multiplicities:
        winner_dimension = max(1, multiplicity // 2)
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        ):
            for kind in ("generic", "sharp_flux", "kernel"):
                record = make_record(
                    kind,
                    length,
                    multiplicity,
                    winner_dimension,
                    seed=70_224 + 101 * length + 17 * multiplicity,
                )
                records.append(record)
                print(
                    "verified repeated circular-normal Jensen face "
                    f"{kind}, length {length}, "
                    f"multiplicity {multiplicity}",
                    flush=True,
                )
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
