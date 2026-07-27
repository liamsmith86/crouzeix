#!/usr/bin/env python3
"""Audit the dimension-free negative second-support Gram at Gau--Wu models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import (
    blaschke_image_jet,
    matrix_polynomial,
    schwarz_map_coefficients,
    support_jet,
)
from gau_wu_finite_model import (
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    random_interior_zeros,
    real_matrix,
)


@dataclass(frozen=True)
class SecondSupportGramRecord:
    """One complete normal-basis audit of the support-energy identity."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    normal_dimension: int
    minimum_boundary_weight: float
    maximum_boundary_weight: float
    minimum_gram_eigenvalue: float
    maximum_gram_eigenvalue: float
    maximum_entry_residual: float
    all_checks_passed: bool


def blaschke_derivative_at_matrix(
    matrix: np.ndarray,
    zeros: np.ndarray,
) -> np.ndarray:
    """Return f'(A) by differentiating the matrix Blaschke product."""

    dimension = matrix.shape[0]
    identity = np.eye(dimension)
    zero = np.zeros_like(matrix)
    velocities = np.zeros(len(zeros), dtype=complex)
    return blaschke_image_jet(
        matrix,
        identity,
        zero,
        zeros,
        velocities,
    )[1]


def boundary_angular_derivative(
    zeros: np.ndarray,
    angles: np.ndarray,
) -> np.ndarray:
    """Return zeta f'(zeta)/f(zeta) as its positive Poisson sum."""

    circle = np.exp(1j * angles)
    return sum(
        (1 - abs(zero) ** 2) / abs(circle - zero) ** 2
        for zero in zeros
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> SecondSupportGramRecord:
    """Audit the complete polarized identity on one random model."""

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
    angles = np.linspace(0, 2 * np.pi, angle_count, endpoint=False)
    weight = boundary_angular_derivative(zeros, angles)
    expected_gram = 2 * np.mean(
        weight[None, None, :] * second_support,
        axis=2,
    )

    derivative = blaschke_derivative_at_matrix(matrix, zeros)
    actual_gram = np.empty_like(expected_gram)
    for row in range(normal.shape[1]):
        for column in range(row, normal.shape[1]):
            second_conformal = matrix_polynomial(
                schwarz_map_coefficients(second_support[row, column]),
                matrix,
            )
            contribution = np.real(
                np.vdot(
                    np.eye(dimension)[:, 0],
                    derivative
                    @ second_conformal
                    @ np.eye(dimension)[:, -1],
                )
            )
            actual_gram[row, column] = contribution
            actual_gram[column, row] = contribution

    eigenvalues = np.linalg.eigvalsh(expected_gram)
    maximum_residual = float(np.max(abs(actual_gram - expected_gram)))
    checks = (
        np.min(weight) > 0
        and eigenvalues[0] > 1e-8
        and maximum_residual < 2e-10
    )
    if not checks:
        raise RuntimeError(
            "second-support Gram audit failed: "
            f"n={dimension}, weight={np.min(weight)}, "
            f"gram={eigenvalues[0]}, residual={maximum_residual}"
        )
    return SecondSupportGramRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        normal_dimension=normal.shape[1],
        minimum_boundary_weight=float(np.min(weight)),
        maximum_boundary_weight=float(np.max(weight)),
        minimum_gram_eigenvalue=float(eigenvalues[0]),
        maximum_gram_eigenvalue=float(eigenvalues[-1]),
        maximum_entry_residual=maximum_residual,
        all_checks_passed=True,
    )


def write_records(
    records: list[SecondSupportGramRecord],
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
        raise argparse.ArgumentTypeError("dimensions must be integers at least 3")
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
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_second_support_gram_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested second-support Gram audit grid."""

    arguments = parse_args()
    records: list[SecondSupportGramRecord] = []
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
                        "minimum_gram_eigenvalue": (
                            record.minimum_gram_eigenvalue
                        ),
                        "maximum_entry_residual": (
                            record.maximum_entry_residual
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
