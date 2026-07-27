#!/usr/bin/env python3
"""Audit the Gau--Wu cyclic corner as the top weighted support mode."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_disk_chart_recenter import polynomial_support_frame
from gau_wu_finite_hessian_jet import (
    blaschke_image_jet,
    normalized_operator_jet,
    support_jet,
)
from gau_wu_second_support_gram import boundary_angular_derivative
from gau_wu_two_sided_endpoint_defect import build_endpoint_form_audit


@dataclass(frozen=True)
class CyclicCornerShapeRecord:
    """One highest-mode/cyclic-corner identity audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    frame_endpoint_error: float
    minimum_abs_inner_derivative_at_zero: float
    maximum_weighted_support_tail: float
    highest_mode_coefficient_error: float
    cyclic_corner_formula_error: float
    zero_velocity_corner_error: float
    all_checks_passed: bool


def polynomial_weight(
    zeros: np.ndarray,
    boundary_points: np.ndarray,
) -> np.ndarray:
    """Return L343's positive polynomial support weight."""

    phi_zeros = np.asarray((0, *zeros), dtype=complex)
    determinant = np.prod(
        1
        - phi_zeros[:, None]
        * np.conj(boundary_points)[None, :],
        axis=0,
    )
    angular_derivative = boundary_angular_derivative(
        zeros,
        np.angle(boundary_points),
    )
    return abs(determinant) ** 2 * angular_derivative


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> CyclicCornerShapeRecord:
    """Audit the exact corner formula on one nondegenerate model."""

    endpoint = build_endpoint_form_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    directions = list(endpoint.directions)
    frame = polynomial_support_frame(endpoint.matrix, endpoint.zeros)
    endpoint_error = max(
        np.linalg.norm(frame[:, 0] - np.eye(dimension)[:, 0]),
        np.linalg.norm(frame[:, -1] - np.eye(dimension)[:, -1]),
    )

    first_support, _, _ = support_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    boundary_points = np.exp(
        2j * np.pi * np.arange(angle_count) / angle_count
    )
    weight = polynomial_weight(endpoint.zeros, boundary_points)
    weighted_support = weight[None, :] * first_support
    weighted_fourier = np.fft.fft(weighted_support, axis=1) / angle_count
    tail = weighted_fourier[
        :,
        dimension + 1 : angle_count - dimension,
    ]
    maximum_tail = float(np.max(abs(tail))) if tail.size else 0.0

    coefficient_corner = np.asarray(
        [
            (
                frame.conj().T @ direction @ frame
            )[-1, 0]
            / 4
            for direction in directions
        ]
    )
    highest_mode = weighted_fourier[:, -dimension]
    highest_mode_error = float(
        np.linalg.norm(highest_mode - coefficient_corner)
        / max(np.linalg.norm(coefficient_corner), 1e-30)
    )

    first_operators, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    zero_count = dimension - 1
    zero_velocities = np.zeros(zero_count, dtype=complex)
    response_corner = np.asarray(
        [
            blaschke_image_jet(
                endpoint.matrix,
                first_operator,
                np.zeros_like(endpoint.matrix),
                endpoint.zeros,
                zero_velocities,
            )[1][-1, 0]
            for first_operator in first_operators
        ]
    )
    inner_derivative_at_zero = np.prod(-endpoint.zeros[:-1])
    expected_corner = 4 * inner_derivative_at_zero * highest_mode
    corner_error = float(
        np.linalg.norm(response_corner - expected_corner)
        / max(np.linalg.norm(expected_corner), 1e-30)
    )

    zero_corner_error = 0.0
    zero_operator = np.zeros_like(endpoint.matrix)
    for index in range(2 * zero_count):
        velocities = np.zeros(zero_count, dtype=complex)
        velocities[index % zero_count] = 1 if index < zero_count else 1j
        response = blaschke_image_jet(
            endpoint.matrix,
            zero_operator,
            zero_operator,
            endpoint.zeros,
            velocities,
        )[1]
        zero_corner_error = max(
            zero_corner_error,
            abs(response[-1, 0]),
        )

    minimum_derivative = float(abs(inner_derivative_at_zero))
    checks = (
        endpoint_error < 2e-12
        and minimum_derivative > 1e-8
        and maximum_tail < 2e-12
        and highest_mode_error < 2e-11
        and corner_error < 2e-10
        and zero_corner_error < 2e-12
    )
    if not checks:
        raise RuntimeError(
            "cyclic-corner shape audit failed: "
            f"n={dimension}, sample={sample}, "
            f"endpoints={endpoint_error}, "
            f"derivative={minimum_derivative}, "
            f"tail={maximum_tail}, "
            f"highest={highest_mode_error}, "
            f"corner={corner_error}, "
            f"zero={zero_corner_error}"
        )

    return CyclicCornerShapeRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        frame_endpoint_error=float(endpoint_error),
        minimum_abs_inner_derivative_at_zero=minimum_derivative,
        maximum_weighted_support_tail=maximum_tail,
        highest_mode_coefficient_error=highest_mode_error,
        cyclic_corner_formula_error=corner_error,
        zero_velocity_corner_error=float(zero_corner_error),
        all_checks_passed=True,
    )


def write_records(
    records: list[CyclicCornerShapeRecord],
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
            "experiments/gau_wu_cyclic_corner_shape_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic highest-mode/corner audit."""

    arguments = parse_args()
    records: list[CyclicCornerShapeRecord] = []
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
                        "highest_mode": (
                            record.highest_mode_coefficient_error
                        ),
                        "cyclic_corner": (
                            record.cyclic_corner_formula_error
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
