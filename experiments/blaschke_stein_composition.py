#!/usr/bin/env python3
"""Audit the finite-Blaschke Stein composition identities from L127."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_discrete_lyapunov


DEFAULT_ZEROS = (-0.31 + 0.08j, 0.12 - 0.17j, 0.4 + 0.03j, -0.22j)


@dataclass(frozen=True)
class BlaschkeSteinRecord:
    dimension: int
    degree: int
    transfer_stein_residual: float
    gramian_composition_residual: float
    minimum_lifted_eigenvalue: float


def blaschke_factor(matrix: np.ndarray, zero: complex) -> np.ndarray:
    """Evaluate one normalized Blaschke factor, omitting a unit phase."""

    identity = np.eye(matrix.shape[0], dtype=complex)
    return (
        (matrix - zero * identity)
        @ np.linalg.inv(identity - np.conjugate(zero) * matrix)
    )


def model_functions(
    matrix: np.ndarray,
    zeros: tuple[complex, ...],
) -> tuple[list[np.ndarray], np.ndarray]:
    """Evaluate a Takenaka--Malmquist basis and its Blaschke product."""

    identity = np.eye(matrix.shape[0], dtype=complex)
    prefix = identity
    functions: list[np.ndarray] = []
    for zero in zeros:
        denominator = np.linalg.inv(
            identity - np.conjugate(zero) * matrix
        )
        functions.append(
            np.sqrt(1 - abs(zero) ** 2) * denominator @ prefix
        )
        prefix = prefix @ blaschke_factor(matrix, zero)
    return functions, prefix


def transfer(
    functions: list[np.ndarray],
    matrix: np.ndarray,
) -> np.ndarray:
    """Apply the positive finite-model-space transfer map."""

    return sum(
        (function.conj().T @ matrix @ function for function in functions),
        np.zeros_like(matrix, dtype=complex),
    )


def make_record(dimension: int, degree: int) -> BlaschkeSteinRecord:
    """Audit one deterministic dimension/degree pair."""

    rng = np.random.default_rng(10_000 + 100 * dimension + degree)
    raw = (
        rng.standard_normal((dimension, dimension))
        + 1j * rng.standard_normal((dimension, dimension))
    )
    matrix = 0.35 * raw / np.linalg.norm(raw, 2)
    factor = (
        rng.standard_normal((dimension, dimension))
        + 1j * rng.standard_normal((dimension, dimension))
    )
    positive = factor.conj().T @ factor + np.eye(dimension)
    zeros = tuple(DEFAULT_ZEROS[:degree])
    functions, blaschke_value = model_functions(matrix, zeros)

    lifted = transfer(functions, positive)
    transfer_defect = (
        lifted
        - matrix.conj().T @ lifted @ matrix
        - positive
        + blaschke_value.conj().T @ positive @ blaschke_value
    )

    defect = (
        rng.standard_normal(dimension)
        + 1j * rng.standard_normal(dimension)
    )
    forcing = np.outer(defect, np.conjugate(defect))
    direct_gramian = solve_discrete_lyapunov(
        matrix.conj().T,
        forcing,
    )
    descended_gramian = solve_discrete_lyapunov(
        blaschke_value.conj().T,
        forcing,
    )
    composed_gramian = transfer(functions, descended_gramian)

    record = BlaschkeSteinRecord(
        dimension=dimension,
        degree=degree,
        transfer_stein_residual=float(np.linalg.norm(transfer_defect, 2)),
        gramian_composition_residual=float(
            np.linalg.norm(direct_gramian - composed_gramian, 2)
        ),
        minimum_lifted_eigenvalue=float(
            np.linalg.eigvalsh(lifted)[0]
        ),
    )
    if record.transfer_stein_residual > 2e-12:
        raise AssertionError("the Blaschke transfer identity failed")
    if record.gramian_composition_residual > 2e-11:
        raise AssertionError("the Stein Gramian composition failed")
    if record.minimum_lifted_eigenvalue <= 0:
        raise AssertionError("the positive transfer lost definiteness")
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-dimension", type=int, default=3)
    parser.add_argument("--maximum-dimension", type=int, default=7)
    parser.add_argument("--maximum-degree", type=int, default=4)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic audit grid."""

    args = parse_args()
    records = [
        make_record(dimension, degree)
        for dimension in range(
            args.minimum_dimension,
            args.maximum_dimension + 1,
        )
        for degree in range(1, args.maximum_degree + 1)
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
