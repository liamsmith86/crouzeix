#!/usr/bin/env python3
"""Audit the Gau--Wu postcomposition/endpoint-trace Euler identity."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_conformal_dual_polarization import first_blaschke_image
from gau_wu_finite_hessian_jet import (
    blaschke_image_jet,
    normalized_operator_jet,
)
from gau_wu_root_tracking_normal_form import (
    inner_derivative_at_simple_zero,
)
from gau_wu_two_sided_endpoint_defect import build_endpoint_form_audit


@dataclass(frozen=True)
class PostcompositionEulerRecord:
    """One postcomposition tangent and Euler-row audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    joint_dimension: int
    postcomposition_first_image_error: float
    endpoint_trace_euler_row_residual: float
    all_checks_passed: bool


def postcomposition_variable(
    value: complex,
    derivatives: np.ndarray,
    physical_dimension: int,
    joint_dimension: int,
) -> np.ndarray:
    """Return the real joint variable for ``b_(t value) o f``."""

    velocities = value / derivatives
    variable = np.zeros(joint_dimension)
    zero_count = len(derivatives)
    variable[
        physical_dimension : physical_dimension + zero_count
    ] = velocities.real
    variable[physical_dimension + zero_count :] = velocities.imag
    return variable


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> PostcompositionEulerRecord:
    """Audit one nondegenerate finite Gau--Wu model."""

    endpoint = build_endpoint_form_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    directions = list(endpoint.directions)
    physical_dimension = len(directions)
    zero_count = dimension - 1
    joint_dimension = dimension**2
    interior_zeros = endpoint.zeros[:-1]
    derivatives = np.asarray(
        [
            *(
                inner_derivative_at_simple_zero(zero, interior_zeros)
                for zero in interior_zeros
            ),
            np.prod(-interior_zeros),
        ]
    )
    real_postcomposition = postcomposition_variable(
        1,
        derivatives,
        physical_dimension,
        joint_dimension,
    )
    imaginary_postcomposition = postcomposition_variable(
        1j,
        derivatives,
        physical_dimension,
        joint_dimension,
    )

    zero_operator = np.zeros_like(endpoint.matrix)
    first_image_error = 0.0
    for value, variable in (
        (1, real_postcomposition),
        (1j, imaginary_postcomposition),
    ):
        velocities = (
            variable[
                physical_dimension : physical_dimension + zero_count
            ]
            + 1j
            * variable[physical_dimension + zero_count :]
        )
        response = blaschke_image_jet(
            endpoint.matrix,
            zero_operator,
            zero_operator,
            endpoint.zeros,
            velocities,
        )[1]
        first_image_error = max(
            first_image_error,
            np.linalg.norm(response + value * np.eye(dimension), 2),
        )

    first_operators, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    identity = np.eye(joint_dimension)
    endpoint_trace_row = np.asarray(
        [
            np.trace(
                first_blaschke_image(
                    identity[index],
                    physical_dimension,
                    first_operators,
                    endpoint.matrix,
                    endpoint.zeros,
                )[[0, -1]][:, [0, -1]]
            )
            for index in range(joint_dimension)
        ]
    )
    expected_trace_row = (4 / 3) * (
        real_postcomposition @ endpoint.sharp_form
        + 1j * imaginary_postcomposition @ endpoint.sharp_form
    )
    row_residual = float(
        np.linalg.norm(endpoint_trace_row - expected_trace_row)
        / np.linalg.norm(endpoint_trace_row)
    )

    checks = (
        joint_dimension == physical_dimension + 2 * zero_count
        and np.min(abs(derivatives)) > 1e-8
        and first_image_error < 2e-11
        and row_residual < 2e-10
    )
    if not checks:
        raise RuntimeError(
            "postcomposition Euler audit failed: "
            f"n={dimension}, sample={sample}, "
            f"first={first_image_error}, row={row_residual}"
        )

    return PostcompositionEulerRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        joint_dimension=joint_dimension,
        postcomposition_first_image_error=float(first_image_error),
        endpoint_trace_euler_row_residual=row_residual,
        all_checks_passed=True,
    )


def write_records(
    records: list[PostcompositionEulerRecord],
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
            "experiments/gau_wu_postcomposition_euler_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic postcomposition Euler audit."""

    arguments = parse_args()
    records: list[PostcompositionEulerRecord] = []
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
                        "first_image": (
                            record.postcomposition_first_image_error
                        ),
                        "euler_row": (
                            record.endpoint_trace_euler_row_residual
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
