#!/usr/bin/env python3
"""Audit the signed complex covariance of L342's endpoint residual map."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import normalized_operator_jet
from gau_wu_similarity_hessian import (
    boundary_residual_map,
    build_similarity_hessian_audit,
    gau_wu_metric,
    residual_metric,
)


@dataclass(frozen=True)
class ResidualComplexCovarianceRecord:
    """One signed endpoint-covariance audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    residual_dimension: int
    residual_map_minimum_singular_value: float
    metric_complex_covariance_error: float
    zero_complex_covariance_error: float
    residual_metric_invariance_error: float
    metric_zero_orthogonality_error: float
    all_checks_passed: bool


def standard_complex_structure(complex_dimension: int) -> np.ndarray:
    """Return multiplication by ``i`` in real/imaginary coordinates."""

    zero = np.zeros((complex_dimension, complex_dimension))
    identity = np.eye(complex_dimension)
    return np.block([[zero, -identity], [identity, zero]])


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> ResidualComplexCovarianceRecord:
    """Audit one finite nondegenerate Gau--Wu model."""

    audit = build_similarity_hessian_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    endpoint = audit.endpoint_audit
    directions = list(endpoint.directions)
    physical_dimension = len(directions)
    complex_dimension = dimension - 1
    boundary_dimension = 2 * complex_dimension
    combined_dimension = physical_dimension + 2 * boundary_dimension
    first_operators, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    residual = boundary_residual_map(
        endpoint,
        first_operators,
        gau_wu_metric(dimension),
        combined_dimension,
    )
    metric_columns = residual[
        :,
        physical_dimension : physical_dimension + boundary_dimension,
    ]
    zero_columns = residual[
        :,
        physical_dimension + boundary_dimension :,
    ]

    phase = standard_complex_structure(complex_dimension)
    zero = np.zeros((boundary_dimension, boundary_dimension))
    signed_endpoint_phase = np.block(
        [[-phase, zero], [zero, phase]]
    )
    metric_covariance = float(
        np.linalg.norm(
            signed_endpoint_phase @ metric_columns
            - metric_columns @ phase,
            2,
        )
        / max(np.linalg.norm(metric_columns, 2), 1e-30)
    )
    zero_covariance = float(
        np.linalg.norm(
            signed_endpoint_phase @ zero_columns
            + zero_columns @ phase,
            2,
        )
        / max(np.linalg.norm(zero_columns, 2), 1e-30)
    )

    square_metric = residual_metric(dimension)
    metric_invariance = float(
        np.linalg.norm(
            signed_endpoint_phase.T
            @ square_metric
            @ signed_endpoint_phase
            - square_metric,
            2,
        )
    )
    orthogonality = float(
        np.linalg.norm(
            metric_columns.T @ square_metric @ zero_columns,
            2,
        )
        / max(
            np.linalg.norm(
                metric_columns.T @ square_metric @ metric_columns,
                2,
            ),
            np.linalg.norm(
                zero_columns.T @ square_metric @ zero_columns,
                2,
            ),
            1e-30,
        )
    )
    minimum_singular = float(
        np.linalg.svd(
            np.column_stack((metric_columns, zero_columns)),
            compute_uv=False,
        )[-1]
    )

    checks = (
        residual.shape
        == (2 * boundary_dimension, combined_dimension)
        and minimum_singular > 1e-6
        and metric_covariance < 2e-10
        and zero_covariance < 2e-10
        and metric_invariance < 1e-13
        and orthogonality < 2e-10
    )
    if not checks:
        raise RuntimeError(
            "residual complex-covariance audit failed: "
            f"n={dimension}, sample={sample}, "
            f"metric={metric_covariance}, zero={zero_covariance}, "
            f"invariance={metric_invariance}, "
            f"orthogonality={orthogonality}, "
            f"minimum={minimum_singular}"
        )

    return ResidualComplexCovarianceRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        residual_dimension=2 * boundary_dimension,
        residual_map_minimum_singular_value=minimum_singular,
        metric_complex_covariance_error=metric_covariance,
        zero_complex_covariance_error=zero_covariance,
        residual_metric_invariance_error=metric_invariance,
        metric_zero_orthogonality_error=orthogonality,
        all_checks_passed=True,
    )


def write_records(
    records: list[ResidualComplexCovarianceRecord],
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
            "experiments/gau_wu_residual_complex_covariance_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic signed-covariance audit."""

    arguments = parse_args()
    records: list[ResidualComplexCovarianceRecord] = []
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
                        "metric_covariance": (
                            record.metric_complex_covariance_error
                        ),
                        "zero_covariance": (
                            record.zero_complex_covariance_error
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
