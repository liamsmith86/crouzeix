#!/usr/bin/env python3
"""Audit the support-weight-to-port-energy identity at Gau--Wu models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import support_jet
from gau_wu_finite_model import (
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    random_interior_zeros,
    real_matrix,
)
from gau_wu_second_support_gram import boundary_angular_derivative


@dataclass(frozen=True)
class SupportPortEnergyRecord:
    """One complete support-factor and port-energy audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    normal_dimension: int
    maximum_support_factor_residual: float
    maximum_null_state_residual: float
    maximum_null_norm_residual: float
    maximum_port_gram_residual: float
    all_checks_passed: bool


@dataclass(frozen=True)
class SupportPortData:
    """Canonical model shift and Ando support-port matrices."""

    shift: np.ndarray
    initial: np.ndarray
    terminal: np.ndarray
    unitary: np.ndarray
    d_matrix: np.ndarray
    e_matrix: np.ndarray
    null_scaling: np.ndarray


def build_support_port_data(matrix: np.ndarray) -> SupportPortData:
    """Build the canonical Gau--Wu support-port realization."""

    dimension = matrix.shape[0]
    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    similarity = np.diag(scaling)
    shift = np.diag(1 / scaling) @ matrix @ similarity
    initial = np.eye(dimension)[:, 0]
    terminal = np.eye(dimension)[:, -1]
    unitary = shift + np.outer(terminal, initial)

    d_values = np.full(dimension, 1 / np.sqrt(2))
    d_values[0] = 0
    d_values[-1] = 1
    e_values = np.full(dimension, 1 / np.sqrt(2))
    e_values[0] = 1
    e_values[-1] = 0
    middle_scaling = np.full(dimension, np.sqrt(2))
    middle_scaling[0] = 1
    middle_scaling[-1] = 1
    return SupportPortData(
        shift=shift,
        initial=initial,
        terminal=terminal,
        unitary=unitary,
        d_matrix=np.diag(d_values),
        e_matrix=np.diag(e_values),
        null_scaling=np.diag(middle_scaling),
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> SupportPortEnergyRecord:
    """Audit the Ando support port at one random Gau--Wu model."""

    rng = np.random.default_rng(seed)
    interior = random_interior_zeros(dimension, rng)
    matrix = gau_wu_model(interior)
    zeros = extremal_zeros(interior)
    normal, _, _ = normal_basis(interior)
    directions = [
        real_matrix(normal[:, index], dimension)
        for index in range(normal.shape[1])
    ]
    _, second_support, _ = support_jet(
        matrix,
        directions,
        angle_count,
    )

    port = build_support_port_data(matrix)

    angles = np.linspace(0, 2 * np.pi, angle_count, endpoint=False)
    weight = boundary_angular_derivative(zeros, angles)
    expected_gram = 2 * np.mean(
        weight[None, None, :] * second_support,
        axis=2,
    )
    actual_gram = np.zeros_like(expected_gram)
    maximum_factor_residual = 0.0
    maximum_null_residual = 0.0
    maximum_norm_residual = 0.0

    identity = np.eye(dimension)
    for angle_index, angle in enumerate(angles):
        phase = np.exp(1j * angle)
        support_factor = (
            port.unitary @ port.d_matrix - phase * port.e_matrix
        )
        support_slack = identity - (
            np.conj(phase) * matrix
            + phase * matrix.conj().T
        ) / 2
        maximum_factor_residual = max(
            maximum_factor_residual,
            np.linalg.norm(
                support_factor.conj().T @ support_factor
                - support_slack
            ),
        )

        model_kernel = np.linalg.solve(
            identity - np.conj(phase) * port.shift,
            port.terminal,
        )
        null_state = port.null_scaling @ model_kernel
        maximum_null_residual = max(
            maximum_null_residual,
            np.linalg.norm(support_factor @ null_state),
        )
        maximum_norm_residual = max(
            maximum_norm_residual,
            abs(np.linalg.norm(null_state) ** 2 - 2 * weight[angle_index]),
        )

        inverse_adjoint = np.linalg.pinv(
            support_factor.conj().T,
            rcond=1e-12,
        )
        responses = np.empty(
            (dimension, len(directions)),
            dtype=complex,
        )
        for index, direction in enumerate(directions):
            support_direction = (
                np.conj(phase) * direction
                + phase * direction.conj().T
            ) / 2
            responses[:, index] = (
                inverse_adjoint @ support_direction @ null_state
            )
        actual_gram += (responses.conj().T @ responses).real / angle_count

    maximum_gram_residual = float(
        np.max(abs(actual_gram - expected_gram))
    )
    checks = (
        maximum_factor_residual < 2e-11
        and maximum_null_residual < 2e-11
        and maximum_norm_residual < 2e-10
        and maximum_gram_residual < 3e-10
    )
    if not checks:
        raise RuntimeError(
            "support port audit failed: "
            f"n={dimension}, factor={maximum_factor_residual}, "
            f"null={maximum_null_residual}, "
            f"norm={maximum_norm_residual}, "
            f"gram={maximum_gram_residual}"
        )
    return SupportPortEnergyRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        normal_dimension=len(directions),
        maximum_support_factor_residual=float(maximum_factor_residual),
        maximum_null_state_residual=float(maximum_null_residual),
        maximum_null_norm_residual=float(maximum_norm_residual),
        maximum_port_gram_residual=maximum_gram_residual,
        all_checks_passed=True,
    )


def write_records(
    records: list[SupportPortEnergyRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_dimensions(value: str) -> tuple[int, ...]:
    """Parse a comma-separated dimension list."""

    dimensions = tuple(int(item) for item in value.split(","))
    if not dimensions or any(dimension < 3 for dimension in dimensions):
        raise argparse.ArgumentTypeError("dimensions must be at least three")
    return dimensions


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dimensions",
        type=parse_dimensions,
        default=(3, 4, 5, 6, 7, 8),
    )
    parser.add_argument("--samples", type=int, default=2)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_support_port_energy_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested support-port audits."""

    arguments = parse_args()
    records: list[SupportPortEnergyRecord] = []
    for dimension in arguments.dimensions:
        for sample in range(arguments.samples):
            seed = arguments.seed + 1009 * dimension + sample
            record = audit_model(
                dimension,
                sample,
                seed,
                arguments.angle_count,
            )
            records.append(record)
            print(
                json.dumps(
                    {
                        "dimension": dimension,
                        "sample": sample,
                        "maximum_null_norm_residual": (
                            record.maximum_null_norm_residual
                        ),
                        "maximum_port_gram_residual": (
                            record.maximum_port_gram_residual
                        ),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
    digest = write_records(records, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "record_count": len(records),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
