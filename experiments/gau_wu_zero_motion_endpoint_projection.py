#!/usr/bin/env python3
"""Audit the exact zero-motion endpoint projection identity."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_model import (
    blaschke_at_matrix,
    extremal_zeros,
    gau_wu_model,
    random_interior_zeros,
)


@dataclass(frozen=True)
class ZeroMotionProjectionRecord:
    """One fixed-model audit across several inner perturbations."""

    dimension: int
    sample: int
    seed: int
    boundary_count: int
    movement_count: int
    minimum_endpoint_defect: float
    maximum_endpoint_gap: float
    maximum_projection_residual: float
    all_checks_passed: bool


def blaschke_boundary(
    zeros: np.ndarray,
    circle: np.ndarray,
) -> np.ndarray:
    """Evaluate a normalized finite Blaschke product on the circle."""

    value = np.ones_like(circle)
    for zero in zeros:
        value *= (circle - zero) / (1 - np.conj(zero) * circle)
    return value


def hardy_projection_norm_squared(
    base_phi_zeros: np.ndarray,
    moving_zeros: np.ndarray,
    boundary_count: int,
) -> float:
    """Return ||P_(phi H2) g||² by an independent Fourier projection."""

    circle = np.exp(
        2j * np.pi * np.arange(boundary_count) / boundary_count
    )
    phi = blaschke_boundary(base_phi_zeros, circle)
    moving = blaschke_boundary(moving_zeros, circle)
    quotient = np.conj(phi) * moving
    coefficients = np.fft.fft(quotient) / boundary_count
    nonnegative = coefficients[: boundary_count // 2]
    return float(np.sum(abs(nonnegative) ** 2))


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    boundary_count: int,
    movement_count: int,
) -> ZeroMotionProjectionRecord:
    """Audit both endpoint defects at one fixed model shift."""

    rng = np.random.default_rng(seed)
    interior = random_interior_zeros(dimension, rng)
    matrix = gau_wu_model(interior)
    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    similarity = np.diag(scaling)
    shift = np.diag(1 / scaling) @ matrix @ similarity

    base_f_zeros = extremal_zeros(interior)
    base_phi_zeros = np.asarray((0, *base_f_zeros), dtype=complex)
    minimum_defect = np.inf
    maximum_gap = 0.0
    maximum_projection_residual = 0.0

    for movement in range(movement_count):
        scale = 0.025 + 0.025 * movement
        for _ in range(10_000):
            displacement = scale * (
                rng.normal(size=dimension - 1)
                + 1j * rng.normal(size=dimension - 1)
            )
            moving_zeros = base_f_zeros + displacement
            if np.max(abs(moving_zeros)) < 0.9:
                break
        else:
            raise RuntimeError("could not generate interior moving zeros")

        image = blaschke_at_matrix(shift, moving_zeros)
        right = 1 - np.linalg.norm(image[:, -1]) ** 2
        left = 1 - np.linalg.norm(image[0, :]) ** 2
        projection = hardy_projection_norm_squared(
            base_phi_zeros,
            moving_zeros,
            boundary_count,
        )
        minimum_defect = min(minimum_defect, right, left)
        maximum_gap = max(maximum_gap, abs(right - left))
        maximum_projection_residual = max(
            maximum_projection_residual,
            abs(right - projection),
            abs(left - projection),
        )

    checks = (
        minimum_defect > -2e-12
        and maximum_gap < 2e-11
        and maximum_projection_residual < 2e-10
    )
    if not checks:
        raise RuntimeError(
            "zero-motion endpoint projection failed: "
            f"n={dimension}, defect={minimum_defect}, "
            f"gap={maximum_gap}, projection={maximum_projection_residual}"
        )
    return ZeroMotionProjectionRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        boundary_count=boundary_count,
        movement_count=movement_count,
        minimum_endpoint_defect=float(minimum_defect),
        maximum_endpoint_gap=float(maximum_gap),
        maximum_projection_residual=float(maximum_projection_residual),
        all_checks_passed=True,
    )


def write_records(
    records: list[ZeroMotionProjectionRecord],
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
    parser.add_argument("--movements", type=int, default=3)
    parser.add_argument("--boundary-count", type=int, default=16384)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "gau_wu_zero_motion_endpoint_projection_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested zero-motion projection audit."""

    arguments = parse_args()
    records: list[ZeroMotionProjectionRecord] = []
    for dimension in arguments.dimensions:
        for sample in range(arguments.samples):
            seed = arguments.seed + 1009 * dimension + sample
            record = audit_model(
                dimension,
                sample,
                seed,
                arguments.boundary_count,
                arguments.movements,
            )
            records.append(record)
            print(
                json.dumps(
                    {
                        "dimension": dimension,
                        "sample": sample,
                        "maximum_endpoint_gap": record.maximum_endpoint_gap,
                        "maximum_projection_residual": (
                            record.maximum_projection_residual
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
