#!/usr/bin/env python3
"""Audit the root-tracking normal form of Gau--Wu diagonal responses."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import (
    blaschke_image_jet,
    normalized_operator_jet,
)
from gau_wu_two_sided_endpoint_defect import build_endpoint_form_audit


@dataclass(frozen=True)
class RootTrackingNormalFormRecord:
    """One diagonal/root-mismatch normal-form audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    joint_dimension: int
    maximum_spectral_projection_error: float
    interior_diagonal_normal_form_residual: float
    endpoint_trace_normal_form_residual: float
    cyclic_response_normalization_error: float
    all_checks_passed: bool


def inner_derivative_at_simple_zero(
    zero: complex,
    interior_zeros: np.ndarray,
) -> complex:
    """Return ``f'(zero)`` for the corresponding simple factor."""

    value = zero / (1 - abs(zero) ** 2)
    for other in interior_zeros:
        if other != zero:
            value *= (zero - other) / (1 - np.conj(other) * zero)
    return value


def simple_spectral_projections(
    matrix: np.ndarray,
    interior_zeros: np.ndarray,
) -> tuple[np.ndarray, ...]:
    """Return the Hermite spectral projections at the simple zeros."""

    dimension = matrix.shape[0]
    identity = np.eye(dimension)
    projections: list[np.ndarray] = []
    for zero in interior_zeros:
        numerator = matrix @ matrix
        denominator = zero**2
        for other in interior_zeros:
            if other != zero:
                numerator = numerator @ (matrix - other * identity)
                denominator *= zero - other
        projections.append(numerator / denominator)
    return tuple(projections)


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> RootTrackingNormalFormRecord:
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
    first_operators, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )

    projections = simple_spectral_projections(
        endpoint.matrix,
        interior_zeros,
    )
    zero_projection = np.eye(dimension, dtype=complex) - sum(
        projections,
        np.zeros_like(endpoint.matrix),
    )
    projection_error = 0.0
    for index, (zero, projection) in enumerate(
        zip(interior_zeros, projections, strict=True)
    ):
        projection_error = max(
            projection_error,
            np.linalg.norm(projection @ projection - projection, 2),
            np.linalg.norm(
                endpoint.matrix @ projection - zero * projection,
                2,
            ),
            max(
                (
                    np.linalg.norm(projection @ other_projection, 2)
                    for other_index, other_projection in enumerate(
                        projections
                    )
                    if other_index != index
                ),
                default=0.0,
            ),
        )
    projection_error = max(
        projection_error,
        np.linalg.norm(
            zero_projection @ zero_projection - zero_projection,
            2,
        ),
        np.linalg.norm(
            np.sum(projections, axis=0) + zero_projection - np.eye(dimension),
            2,
        ),
    )

    inner_derivative_zero = np.prod(-interior_zeros)
    simple_derivatives = np.asarray(
        [
            inner_derivative_at_simple_zero(zero, interior_zeros)
            for zero in interior_zeros
        ]
    )

    southwest = np.zeros_like(endpoint.matrix)
    southwest[-1, 0] = 1
    southwest_response = blaschke_image_jet(
        endpoint.matrix,
        southwest,
        np.zeros_like(endpoint.matrix),
        endpoint.zeros,
        np.zeros(zero_count, dtype=complex),
    )[1]
    corner_normalization_error = float(
        abs(southwest_response[-1, 0] - inner_derivative_zero)
        / abs(inner_derivative_zero)
    )

    interior_corner_coefficients = np.asarray(
        [
            (
                southwest_response[index + 1, index + 1]
                - simple_derivatives[index]
                * np.trace(projections[index] @ southwest)
            )
            / inner_derivative_zero
            for index in range(dimension - 2)
        ]
    )
    endpoint_corner_coefficient = (
        southwest_response[0, 0]
        + southwest_response[-1, -1]
        - inner_derivative_zero
        * np.trace(zero_projection @ southwest)
    ) / inner_derivative_zero

    interior_actual: list[np.ndarray] = []
    interior_expected: list[np.ndarray] = []
    endpoint_actual: list[complex] = []
    endpoint_expected: list[complex] = []
    identity = np.eye(joint_dimension)
    for basis_index in range(joint_dimension):
        variable = identity[basis_index]
        physical = variable[:physical_dimension]
        velocities = (
            variable[
                physical_dimension : physical_dimension + zero_count
            ]
            + 1j
            * variable[physical_dimension + zero_count :]
        )
        first_operator = sum(
            (
                physical[index] * first_operators[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(endpoint.matrix),
        )
        response = blaschke_image_jet(
            endpoint.matrix,
            first_operator,
            np.zeros_like(endpoint.matrix),
            endpoint.zeros,
            velocities,
        )[1]
        corner = response[-1, 0]

        actual_diagonal = np.diag(response)[1:-1]
        expected_diagonal = np.asarray(
            [
                simple_derivatives[index]
                * (
                    np.trace(projections[index] @ first_operator)
                    - velocities[index]
                )
                + interior_corner_coefficients[index] * corner
                for index in range(dimension - 2)
            ]
        )
        interior_actual.append(actual_diagonal)
        interior_expected.append(expected_diagonal)

        actual_trace = response[0, 0] + response[-1, -1]
        expected_trace = inner_derivative_zero * (
            np.trace(zero_projection @ first_operator)
            - 2 * velocities[-1]
        ) + endpoint_corner_coefficient * corner
        endpoint_actual.append(actual_trace)
        endpoint_expected.append(expected_trace)

    interior_actual_matrix = np.stack(interior_actual, axis=1)
    interior_expected_matrix = np.stack(interior_expected, axis=1)
    interior_residual = float(
        np.linalg.norm(
            interior_actual_matrix - interior_expected_matrix,
            2,
        )
        / max(np.linalg.norm(interior_actual_matrix, 2), 1e-30)
    )
    endpoint_actual_row = np.asarray(endpoint_actual)
    endpoint_expected_row = np.asarray(endpoint_expected)
    endpoint_residual = float(
        np.linalg.norm(endpoint_actual_row - endpoint_expected_row)
        / max(np.linalg.norm(endpoint_actual_row), 1e-30)
    )

    checks = (
        joint_dimension == physical_dimension + 2 * zero_count
        and projection_error < 5e-9
        and abs(inner_derivative_zero) > 1e-8
        and interior_residual < 2e-11
        and endpoint_residual < 2e-11
        and corner_normalization_error < 2e-12
    )
    if not checks:
        raise RuntimeError(
            "root-tracking normal-form audit failed: "
            f"n={dimension}, sample={sample}, "
            f"projection={projection_error}, "
            f"interior={interior_residual}, "
            f"endpoint={endpoint_residual}, "
            f"corner={corner_normalization_error}"
        )

    return RootTrackingNormalFormRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        joint_dimension=joint_dimension,
        maximum_spectral_projection_error=float(projection_error),
        interior_diagonal_normal_form_residual=interior_residual,
        endpoint_trace_normal_form_residual=endpoint_residual,
        cyclic_response_normalization_error=corner_normalization_error,
        all_checks_passed=True,
    )


def write_records(
    records: list[RootTrackingNormalFormRecord],
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
            "experiments/gau_wu_root_tracking_normal_form_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic root-tracking normal-form audit."""

    arguments = parse_args()
    records: list[RootTrackingNormalFormRecord] = []
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
                        "interior": (
                            record.interior_diagonal_normal_form_residual
                        ),
                        "endpoint": (
                            record.endpoint_trace_normal_form_residual
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
