#!/usr/bin/env python3
"""Audit the exact parallel sum of the Gau--Wu endpoint residuals."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_canonical_conjugation import canonical_conjugation
from gau_wu_canonical_residual_normal_form import (
    complex_endpoint_blocks,
)
from gau_wu_finite_hessian_jet import normalized_operator_jet
from gau_wu_shape_tracking_lagrangian import (
    first_image_rows,
    tracking_rows,
)
from gau_wu_similarity_hessian import (
    boundary_residual_map,
    build_similarity_hessian_audit,
    gau_wu_metric,
    residual_metric,
)


@dataclass(frozen=True)
class ResidualParallelSumRecord:
    """One exact endpoint-weight and parallel-sum audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    physical_dimension: int
    endpoint_weight_intertwining_error: float
    metric_lower_minimum_singular_value: float
    parallel_sum_form_residual: float
    scalar_gap_form_residual: float
    optimizer_equation_residual: float
    tracking_condition_number: float
    all_checks_passed: bool


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> ResidualParallelSumRecord:
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
    endpoint_conjugation = canonical_conjugation(endpoint.matrix)[
        :complex_dimension,
        1:,
    ]
    conjugated_sum = upper + endpoint_conjugation @ lower.conj()

    physical = np.arange(physical_dimension)
    metric = np.arange(
        physical_dimension,
        physical_dimension + boundary_dimension,
    )
    zero = np.arange(
        physical_dimension + boundary_dimension,
        combined_dimension,
    )
    joint = np.concatenate((physical, zero))

    lower_weight = np.diag(
        [1.0] * (complex_dimension - 1) + [3.0]
    )
    upper_weight = np.diag(
        [3.0] + [1.0] * (complex_dimension - 1)
    )
    weight_intertwining = float(
        np.linalg.norm(
            endpoint_conjugation.conj().T
            @ upper_weight
            @ endpoint_conjugation
            - lower_weight,
            2,
        )
    )

    lower_metric_realification = np.vstack(
        (lower[:, metric].real, lower[:, metric].imag)
    )
    lower_metric_singular = np.linalg.svd(
        lower_metric_realification,
        compute_uv=False,
    )
    minimum_metric_singular = float(lower_metric_singular[-1])

    square_form = residual.T @ residual_metric(dimension) @ residual
    metric_block = square_form[np.ix_(metric, metric)]
    mixed_block = square_form[np.ix_(joint, metric)]
    joint_block = square_form[np.ix_(joint, joint)]
    optimized_square = (
        joint_block
        - mixed_block @ np.linalg.solve(metric_block, mixed_block.T)
    )

    joint_defect = conjugated_sum[:, joint]
    predicted_square = (
        0.5
        * (
            joint_defect.conj().T
            @ upper_weight
            @ joint_defect
        ).real
    )
    parallel_sum_residual = float(
        np.linalg.norm(optimized_square - predicted_square, 2)
        / max(np.linalg.norm(optimized_square, 2), 1e-30)
    )

    embedded_similarity = np.zeros_like(endpoint.sharp_form)
    embedded_similarity[:physical_dimension, :physical_dimension] = (
        audit.optimized_similarity_form
    )
    scalar_gap = embedded_similarity - 4 * endpoint.sharp_form
    scalar_gap_residual = float(
        np.linalg.norm(scalar_gap - predicted_square, 2)
        / max(np.linalg.norm(scalar_gap, 2), 1e-30)
    )

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
    physical_defect = (
        conjugated_sum[:, :physical_dimension]
        - tracking_corrector
        @ tracking[:, :physical_dimension].conj()
    )

    zero_block = endpoint.sharp_form[
        physical_dimension:,
        physical_dimension:,
    ]
    physical_zero_block = endpoint.sharp_form[
        :physical_dimension,
        physical_dimension:,
    ]
    zero_optimizer = -np.linalg.solve(
        zero_block,
        physical_zero_block.T,
    )
    optimized_tracking = (
        tracking[:, :physical_dimension]
        + tracking[:, physical_dimension:] @ zero_optimizer
    )
    optimizer_equation = float(
        np.linalg.norm(
            tracking_corrector @ optimized_tracking.conj()
            + physical_defect,
            2,
        )
        / max(np.linalg.norm(physical_defect, 2), 1e-30)
    )

    checks = (
        weight_intertwining < 5e-10
        and minimum_metric_singular > 1e-5
        and parallel_sum_residual < 5e-9
        and scalar_gap_residual < 5e-9
        and optimizer_equation < 5e-9
    )
    if not checks:
        raise RuntimeError(
            "residual parallel-sum audit failed: "
            f"n={dimension}, sample={sample}, "
            f"weight={weight_intertwining}, "
            f"metric={minimum_metric_singular}, "
            f"parallel={parallel_sum_residual}, "
            f"gap={scalar_gap_residual}, "
            f"optimizer={optimizer_equation}"
        )

    return ResidualParallelSumRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        physical_dimension=physical_dimension,
        endpoint_weight_intertwining_error=weight_intertwining,
        metric_lower_minimum_singular_value=minimum_metric_singular,
        parallel_sum_form_residual=parallel_sum_residual,
        scalar_gap_form_residual=scalar_gap_residual,
        optimizer_equation_residual=optimizer_equation,
        tracking_condition_number=float(
            np.linalg.cond(real_zero_tracking)
        ),
        all_checks_passed=True,
    )


def write_records(
    records: list[ResidualParallelSumRecord],
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
            "experiments/gau_wu_residual_parallel_sum_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic residual parallel-sum audit."""

    arguments = parse_args()
    records: list[ResidualParallelSumRecord] = []
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
                        "parallel_sum": (
                            record.parallel_sum_form_residual
                        ),
                        "scalar_gap": record.scalar_gap_form_residual,
                        "optimizer": (
                            record.optimizer_equation_residual
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
