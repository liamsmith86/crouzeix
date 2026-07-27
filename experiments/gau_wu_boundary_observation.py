#!/usr/bin/env python3
"""Identify the observable boundary quotient of the Gau--Wu loss."""

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
    support_jet,
)
from gau_wu_finite_model import real_vectorize
from gau_wu_physical_loss_rank import build_physical_loss_audit


@dataclass(frozen=True)
class BoundaryObservationRecord:
    """One complete loss-kernel/observation-kernel comparison."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    physical_dimension: int
    expected_observation_rank: int
    observed_observation_rank: int
    expected_nullity: int
    observed_nullity: int
    endpoint_observation_rank: int
    support_rank_increment: int
    kernel_projector_residual: float
    maximum_kernel_support_residual: float
    maximum_kernel_endpoint_residual: float
    maximum_kernel_image_residual: float
    all_checks_passed: bool


def real_vector(vector: np.ndarray) -> np.ndarray:
    """Return interlaced real coordinates of one complex vector."""

    return real_vectorize(vector.reshape(-1, 1))


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> BoundaryObservationRecord:
    """Compare the reduced-loss kernel with its boundary observations."""

    audit = build_physical_loss_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    endpoint = audit.endpoint_audit
    directions = list(endpoint.directions)
    physical_dimension = len(directions)
    zero_count = dimension - 1

    first_support, _, _ = support_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    first_operator, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )

    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    similarity = np.diag(scaling)
    inverse_similarity = np.diag(1 / scaling)
    p = np.eye(dimension)[:, 0]
    q = np.eye(dimension)[:, -1]
    zero_matrix = np.zeros_like(endpoint.matrix)

    endpoint_columns: list[np.ndarray] = []
    image_columns: list[np.ndarray] = []
    for index in range(physical_dimension):
        zero_coordinate = audit.zero_optimizer[:, index]
        velocities = (
            zero_coordinate[:zero_count]
            + 1j * zero_coordinate[zero_count:]
        )
        image = blaschke_image_jet(
            endpoint.matrix,
            first_operator[index],
            zero_matrix,
            endpoint.zeros,
            velocities,
        )[1]
        image = inverse_similarity @ image @ similarity
        image_columns.append(real_vectorize(image))
        endpoint_columns.append(
            np.concatenate(
                (
                    real_vector(image @ q),
                    real_vector(image.conj().T @ p),
                )
            )
        )

    endpoint_map = np.stack(endpoint_columns, axis=1)
    image_map = np.stack(image_columns, axis=1)
    observation_map = np.vstack((first_support.T, endpoint_map))

    loss_eigenvalues, loss_eigenvectors = np.linalg.eigh(
        audit.loss_form
    )
    loss_threshold = 1e-7 * loss_eigenvalues[-1]
    loss_kernel = loss_eigenvectors[
        :,
        loss_eigenvalues < loss_threshold,
    ]

    observation_singular_values = np.linalg.svd(
        observation_map,
        compute_uv=False,
    )
    observation_threshold = 1e-8 * observation_singular_values[0]
    observation_rank = int(
        np.sum(observation_singular_values > observation_threshold)
    )
    _, _, observation_right = np.linalg.svd(
        observation_map,
        full_matrices=True,
    )
    observation_kernel = observation_right[observation_rank:].T

    endpoint_singular_values = np.linalg.svd(
        endpoint_map,
        compute_uv=False,
    )
    # The normal basis uses finite-differenced model-zero tangents, so
    # exact zero endpoint modes acquire O(1e-7) singular values in the
    # largest sample.  The active endpoint gap remains above 1e-4.
    endpoint_threshold = 1e-6 * endpoint_singular_values[0]
    endpoint_rank = int(
        np.sum(endpoint_singular_values > endpoint_threshold)
    )

    expected_rank = 6 * dimension - 14
    expected_nullity = (dimension - 4) ** 2
    expected_endpoint_rank = 4 * dimension - 8
    expected_support_increment = 2 * dimension - 6
    projector_residual = float(
        np.linalg.norm(
            loss_kernel @ loss_kernel.T
            - observation_kernel @ observation_kernel.T
        )
    )
    support_residual = float(
        np.linalg.norm(first_support.T @ observation_kernel)
    )
    endpoint_residual = float(
        np.linalg.norm(endpoint_map @ observation_kernel)
    )
    image_residual = float(
        np.linalg.norm(image_map @ observation_kernel)
    )

    checks = (
        observation_rank == expected_rank
        and endpoint_rank == expected_endpoint_rank
        and observation_rank - endpoint_rank
        == expected_support_increment
        and observation_kernel.shape[1] == expected_nullity
        and loss_kernel.shape[1] == expected_nullity
        and projector_residual < 2e-7
        and support_residual < 2e-7
        and endpoint_residual < 2e-7
        and image_residual < 3e-7
    )
    if not checks:
        raise RuntimeError(
            "boundary-observation audit failed: "
            f"n={dimension}, rank={observation_rank}/{expected_rank}, "
            f"null={observation_kernel.shape[1]}/{expected_nullity}, "
            f"projector={projector_residual}, "
            f"support={support_residual}, "
            f"endpoint={endpoint_residual}, image={image_residual}"
        )
    return BoundaryObservationRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        physical_dimension=physical_dimension,
        expected_observation_rank=expected_rank,
        observed_observation_rank=observation_rank,
        expected_nullity=expected_nullity,
        observed_nullity=observation_kernel.shape[1],
        endpoint_observation_rank=endpoint_rank,
        support_rank_increment=observation_rank - endpoint_rank,
        kernel_projector_residual=projector_residual,
        maximum_kernel_support_residual=support_residual,
        maximum_kernel_endpoint_residual=endpoint_residual,
        maximum_kernel_image_residual=image_residual,
        all_checks_passed=True,
    )


def write_records(
    records: list[BoundaryObservationRecord],
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
    if not dimensions or any(dimension < 4 for dimension in dimensions):
        raise argparse.ArgumentTypeError("dimensions must be at least four")
    return dimensions


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dimensions",
        type=parse_dimensions,
        default=(4, 5, 6, 7, 8, 9),
    )
    parser.add_argument("--samples", type=int, default=2)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_boundary_observation_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested boundary-observation probes."""

    arguments = parse_args()
    records: list[BoundaryObservationRecord] = []
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
                        "observation_rank": (
                            record.observed_observation_rank
                        ),
                        "observation_nullity": record.observed_nullity,
                        "projector_residual": (
                            record.kernel_projector_residual
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
