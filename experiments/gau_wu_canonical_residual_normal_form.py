#!/usr/bin/env python3
"""Audit the canonical normal form of the Gau--Wu endpoint residuals."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_canonical_conjugation import canonical_conjugation
from gau_wu_conformal_shape_phase import conformal_shape_map
from gau_wu_finite_hessian_jet import (
    normalized_operator_jet,
    support_jet,
)
from gau_wu_shape_tracking_lagrangian import (
    first_image_rows,
    tracking_rows,
)
from gau_wu_similarity_hessian import (
    boundary_residual_map,
    build_similarity_hessian_audit,
    gau_wu_metric,
)


@dataclass(frozen=True)
class CanonicalResidualNormalFormRecord:
    """One exact homogeneous reduction and physical-factor audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    physical_dimension: int
    corrected_output_dimension: int
    canonical_endpoint_unitarity_error: float
    metric_graph_residual: float
    zero_tracking_normal_form_residual: float
    endpoint_trace_residual: float
    corrected_homogeneous_residual: float
    physical_shape_factor_residual: float
    disk_kernel_residual: float
    observed_physical_response_rank: int
    expected_physical_response_rank: int
    tracking_condition_number: float
    all_checks_passed: bool


def complex_endpoint_blocks(
    residual: np.ndarray,
    complex_dimension: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the lower and upper complex endpoint blocks."""

    lower = (
        residual[:complex_dimension]
        + 1j * residual[complex_dimension : 2 * complex_dimension]
    )
    upper = (
        residual[2 * complex_dimension : 3 * complex_dimension]
        + 1j * residual[3 * complex_dimension :]
    )
    return lower, upper


def real_rank(matrix: np.ndarray, relative_threshold: float = 1e-8) -> int:
    """Return the numerical real rank of a complex-column map."""

    realification = np.vstack((matrix.real, matrix.imag))
    singular_values = np.linalg.svd(realification, compute_uv=False)
    if not len(singular_values) or singular_values[0] == 0:
        return 0
    return int(
        np.sum(singular_values > relative_threshold * singular_values[0])
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> CanonicalResidualNormalFormRecord:
    """Audit one finite nondegenerate Gau--Wu equality model."""

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
    first_support, _, _ = support_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    real_shape = conformal_shape_map(first_support, dimension)
    complex_shape = (
        real_shape[:complex_dimension]
        + 1j * real_shape[complex_dimension:]
    )
    image_rows = first_image_rows(
        endpoint.matrix,
        endpoint.zeros,
        first_operators,
    )
    tracking = tracking_rows(
        endpoint.matrix,
        endpoint.zeros,
        first_operators,
        image_rows,
    )

    residual = boundary_residual_map(
        endpoint,
        first_operators,
        gau_wu_metric(dimension),
        combined_dimension,
    )
    lower, upper = complex_endpoint_blocks(
        residual,
        complex_dimension,
    )
    conjugation = canonical_conjugation(endpoint.matrix)
    endpoint_conjugation = conjugation[
        :complex_dimension,
        1:,
    ]
    conjugated_sum = upper + endpoint_conjugation @ lower.conj()

    physical = slice(0, physical_dimension)
    metric = slice(
        physical_dimension,
        physical_dimension + boundary_dimension,
    )
    zero = slice(physical_dimension + boundary_dimension, None)

    endpoint_unitarity = float(
        np.linalg.norm(
            endpoint_conjugation
            @ endpoint_conjugation.conj().T
            - np.eye(complex_dimension),
            2,
        )
    )
    metric_graph = float(
        np.linalg.norm(conjugated_sum[:, metric], 2)
        / max(np.linalg.norm(residual[:, metric], 2), 1e-30)
    )

    # On pure zero motion, tracking is a complex isomorphism.  The
    # conjugated endpoint sum is antiholomorphic in that coordinate.
    real_zero_tracking = tracking[
        :,
        physical_dimension : physical_dimension + complex_dimension,
    ]
    real_zero_defect = conjugated_sum[
        :,
        physical_dimension
        + boundary_dimension : physical_dimension
        + boundary_dimension
        + complex_dimension,
    ]
    tracking_corrector = np.linalg.solve(
        real_zero_tracking.conj().T,
        real_zero_defect.T,
    ).T

    zero_tracking_normal_form = float(
        np.linalg.norm(
            conjugated_sum[:, zero]
            - tracking_corrector @ tracking[:, physical_dimension:].conj(),
            2,
        )
        / max(np.linalg.norm(conjugated_sum[:, zero], 2), 1e-30)
    )

    embedded_tracking = np.hstack(
        (
            tracking[:, :physical_dimension],
            np.zeros(
                (complex_dimension, boundary_dimension),
                dtype=complex,
            ),
            tracking[:, physical_dimension:],
        )
    )
    endpoint_trace = float(
        np.linalg.norm(
            conjugated_sum[0] - embedded_tracking[0].conj()
        )
        / max(np.linalg.norm(conjugated_sum[0]), 1e-30)
    )
    corrected = (
        conjugated_sum
        - tracking_corrector @ embedded_tracking.conj()
    )
    corrected_homogeneous = float(
        np.linalg.norm(
            corrected[:, physical_dimension:],
            2,
        )
        / max(np.linalg.norm(conjugated_sum[:, physical_dimension:], 2), 1e-30)
    )

    # The following factorization is the remaining physical identity.
    # It is audited numerically but is deliberately not used as a theorem.
    physical_response = corrected[:, physical]
    shape_factor = physical_response @ np.linalg.pinv(complex_shape)
    shape_factor_residual = float(
        np.linalg.norm(
            physical_response - shape_factor @ complex_shape,
            2,
        )
        / max(np.linalg.norm(physical_response, 2), 1e-30)
    )

    _, shape_singular, shape_right = np.linalg.svd(
        real_shape,
        full_matrices=True,
    )
    shape_rank = int(
        np.sum(shape_singular > 1e-8 * shape_singular[0])
    )
    disk_kernel = shape_right[shape_rank:].T
    disk_kernel_residual = float(
        np.linalg.norm(physical_response @ disk_kernel, 2)
        / max(np.linalg.norm(physical_response, 2), 1e-30)
    )
    observed_response_rank = real_rank(physical_response)
    expected_response_rank = 2 * dimension - 4

    checks = (
        endpoint_unitarity < 5e-10
        and metric_graph < 5e-10
        and zero_tracking_normal_form < 5e-10
        and endpoint_trace < 5e-10
        and corrected_homogeneous < 5e-10
        and shape_factor_residual < 5e-9
        and disk_kernel_residual < 5e-9
        and observed_response_rank == expected_response_rank
    )
    if not checks:
        raise RuntimeError(
            "canonical residual normal-form audit failed: "
            f"n={dimension}, sample={sample}, "
            f"unitarity={endpoint_unitarity}, metric={metric_graph}, "
            f"zero={zero_tracking_normal_form}, trace={endpoint_trace}, "
            f"homogeneous={corrected_homogeneous}, "
            f"factor={shape_factor_residual}, "
            f"disk={disk_kernel_residual}, "
            f"rank={observed_response_rank}/{expected_response_rank}"
        )

    return CanonicalResidualNormalFormRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        physical_dimension=physical_dimension,
        corrected_output_dimension=dimension - 2,
        canonical_endpoint_unitarity_error=endpoint_unitarity,
        metric_graph_residual=metric_graph,
        zero_tracking_normal_form_residual=zero_tracking_normal_form,
        endpoint_trace_residual=endpoint_trace,
        corrected_homogeneous_residual=corrected_homogeneous,
        physical_shape_factor_residual=shape_factor_residual,
        disk_kernel_residual=disk_kernel_residual,
        observed_physical_response_rank=observed_response_rank,
        expected_physical_response_rank=expected_response_rank,
        tracking_condition_number=float(
            np.linalg.cond(real_zero_tracking)
        ),
        all_checks_passed=True,
    )


def write_records(
    records: list[CanonicalResidualNormalFormRecord],
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
            "experiments/"
            "gau_wu_canonical_residual_normal_form_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic canonical residual audit."""

    arguments = parse_args()
    records: list[CanonicalResidualNormalFormRecord] = []
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
                        "metric_graph": record.metric_graph_residual,
                        "zero_normal_form": (
                            record.zero_tracking_normal_form_residual
                        ),
                        "physical_factor": (
                            record.physical_shape_factor_residual
                        ),
                        "physical_rank": (
                            record.observed_physical_response_rank
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
