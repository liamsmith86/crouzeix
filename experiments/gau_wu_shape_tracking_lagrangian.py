#!/usr/bin/env python3
"""Audit the joint shape/tracking chart and its Hardy-Lagrangian pattern."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_conformal_shape_phase import conformal_shape_map
from gau_wu_finite_hessian_jet import (
    blaschke_image_jet,
    normalized_operator_jet,
    support_jet,
)
from gau_wu_root_tracking_normal_form import simple_spectral_projections
from gau_wu_similarity_hessian import build_similarity_hessian_audit


@dataclass(frozen=True)
class ShapeTrackingLagrangianRecord:
    """One exact-chart and numerical-Lagrangian audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    joint_dimension: int
    expected_active_dimension: int
    observed_coordinate_rank: int
    coordinate_kernel_dimension: int
    tracking_zero_restriction_error: float
    joint_form_condition_number: float
    dual_solve_backward_residual: float
    total_bilinear_isotropy_residual: float
    shape_bilinear_isotropy_residual: float
    shape_tracking_isotropy_residual: float
    tracking_bilinear_isotropy_residual: float
    optimized_tracking_holomorphic_residual: float
    riesz_lower_flag_residual: float
    all_checks_passed: bool


def first_image_rows(
    matrix: np.ndarray,
    zeros: np.ndarray,
    first_operators: list[np.ndarray],
) -> np.ndarray:
    """Return every first-image entry as a row on the joint space."""

    dimension = matrix.shape[0]
    physical_dimension = len(first_operators)
    zero_count = len(zeros)
    joint_dimension = physical_dimension + 2 * zero_count
    rows = np.empty((dimension, dimension, joint_dimension), dtype=complex)
    for index in range(joint_dimension):
        variable = np.eye(joint_dimension)[index]
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
                physical[offset] * first_operators[offset]
                for offset in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        rows[:, :, index] = blaschke_image_jet(
            matrix,
            first_operator,
            np.zeros_like(matrix),
            zeros,
            velocities,
        )[1]
    return rows


def tracking_rows(
    matrix: np.ndarray,
    zeros: np.ndarray,
    first_operators: list[np.ndarray],
    image_rows: np.ndarray,
) -> np.ndarray:
    """Return endpoint-trace and simple-root mismatch rows."""

    physical_dimension = len(first_operators)
    zero_count = len(zeros)
    joint_dimension = physical_dimension + 2 * zero_count
    projections = simple_spectral_projections(matrix, zeros[:-1])
    rows = [image_rows[0, 0] + image_rows[-1, -1]]
    for index, projection in enumerate(projections):
        row = np.zeros(joint_dimension, dtype=complex)
        row[:physical_dimension] = [
            np.trace(projection @ first_operator)
            for first_operator in first_operators
        ]
        row[physical_dimension + index] = -1
        row[physical_dimension + zero_count + index] = -1j
        rows.append(row)
    return np.asarray(rows)


def spectral_velocity_map(
    matrix: np.ndarray,
    zeros: np.ndarray,
    first_operators: list[np.ndarray],
) -> np.ndarray:
    """Return the real map from physical directions to tracked roots."""

    projections = simple_spectral_projections(matrix, zeros[:-1])
    repeated_projection = np.eye(matrix.shape[0]) - sum(
        projections,
        np.zeros_like(matrix),
    )
    complex_map = np.asarray(
        [
            *(
                [
                    np.trace(projection @ first_operator)
                    for first_operator in first_operators
                ]
                for projection in projections
            ),
            [
                np.trace(repeated_projection @ first_operator) / 2
                for first_operator in first_operators
            ],
        ]
    )
    return np.vstack((complex_map.real, complex_map.imag))


def complex_block_residual(
    block: np.ndarray,
    reference: np.ndarray,
) -> float:
    """Normalize one pure complex-bilinear inverse-Hessian block."""

    return float(
        np.linalg.norm(block, 2)
        / max(np.linalg.norm(reference, 2), 1e-30)
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> ShapeTrackingLagrangianRecord:
    """Audit one finite nondegenerate Gau--Wu equality model."""

    audit = build_similarity_hessian_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    endpoint = audit.endpoint_audit
    matrix = endpoint.matrix
    zeros = endpoint.zeros
    directions = list(endpoint.directions)
    physical_dimension = len(directions)
    zero_count = dimension - 1
    joint_dimension = dimension**2
    expected_active = 4 * zero_count

    first_operators, _, _ = normalized_operator_jet(
        matrix,
        directions,
        angle_count,
    )
    first_support, _, _ = support_jet(
        matrix,
        directions,
        angle_count,
    )
    physical_shape = conformal_shape_map(first_support, dimension)
    joint_shape = np.hstack(
        (
            physical_shape,
            np.zeros((2 * zero_count, 2 * zero_count)),
        )
    )
    image_rows = first_image_rows(matrix, zeros, first_operators)
    complex_tracking = tracking_rows(
        matrix,
        zeros,
        first_operators,
        image_rows,
    )
    real_tracking = np.vstack(
        (complex_tracking.real, complex_tracking.imag)
    )
    coordinates = np.vstack((joint_shape, real_tracking))
    normalized_coordinates = coordinates / np.linalg.norm(
        coordinates,
        axis=1,
    )[:, None]
    singular_values = np.linalg.svd(
        normalized_coordinates,
        compute_uv=False,
    )
    rank_threshold = 1e-8 * singular_values[0]
    coordinate_rank = int(np.sum(singular_values > rank_threshold))

    zero_complex = (
        complex_tracking[
            :, physical_dimension : physical_dimension + zero_count
        ]
    )
    zero_imaginary = complex_tracking[:, physical_dimension + zero_count :]
    derivative_at_zero = np.prod(-zeros[:-1])
    expected_zero = np.zeros_like(zero_complex)
    expected_zero[0, -1] = -2 * derivative_at_zero
    for index in range(zero_count - 1):
        expected_zero[index + 1, index] = -1
    zero_restriction_error = float(
        max(
            np.linalg.norm(zero_complex - expected_zero, 2),
            np.linalg.norm(zero_imaginary - 1j * expected_zero, 2),
        )
    )

    shape_x = physical_shape[:zero_count]
    shape_y = physical_shape[zero_count:]
    analytic_shape_source = np.hstack(
        (
            shape_y + 1j * shape_x,
            np.zeros((zero_count, 2 * zero_count)),
        )
    )
    analytic_sources = np.vstack(
        (analytic_shape_source, complex_tracking)
    )
    joint_form = endpoint.sharp_form
    dual_lifts = np.linalg.solve(joint_form, analytic_sources.T)
    backward_residual = float(
        np.linalg.norm(
            joint_form @ dual_lifts - analytic_sources.T,
            2,
        )
        / (
            np.linalg.norm(joint_form, 2) * np.linalg.norm(dual_lifts, 2)
            + np.linalg.norm(analytic_sources, 2)
        )
    )
    pure_pairing = analytic_sources @ dual_lifts
    conjugate_lifts = np.linalg.solve(
        joint_form,
        analytic_sources.conj().T,
    )
    mixed_pairing = analytic_sources @ conjugate_lifts
    total_isotropy = complex_block_residual(pure_pairing, mixed_pairing)
    shape_isotropy = complex_block_residual(
        pure_pairing[:zero_count, :zero_count],
        mixed_pairing[:zero_count, :zero_count],
    )
    cross_isotropy = complex_block_residual(
        pure_pairing[:zero_count, zero_count:],
        mixed_pairing[:zero_count, zero_count:],
    )
    tracking_isotropy = complex_block_residual(
        pure_pairing[zero_count:, zero_count:],
        mixed_pairing[zero_count:, zero_count:],
    )

    physical = joint_form[:physical_dimension, :physical_dimension]
    mixed = joint_form[:physical_dimension, physical_dimension:]
    zero_form = joint_form[physical_dimension:, physical_dimension:]
    optimized_physical = (
        physical - mixed @ np.linalg.solve(zero_form, mixed.T)
    )
    constrained_lift = np.linalg.solve(
        optimized_physical,
        physical_shape.T,
    )
    constrained_lift @= np.linalg.inv(
        physical_shape @ constrained_lift
    )
    optimized_zero = -np.linalg.solve(
        zero_form,
        mixed.T @ constrained_lift,
    )
    mismatch = (
        optimized_zero
        - spectral_velocity_map(matrix, zeros, first_operators)
        @ constrained_lift
    )
    xx = mismatch[:zero_count, :zero_count]
    xy = mismatch[:zero_count, zero_count:]
    yx = mismatch[zero_count:, :zero_count]
    yy = mismatch[zero_count:, zero_count:]
    holomorphic = (xx + yy + 1j * (yx - xy)) / 2
    antiholomorphic = (xx - yy + 1j * (yx + xy)) / 2
    optimizer_holomorphic_residual = float(
        np.linalg.norm(holomorphic, 2)
        / max(np.linalg.norm(antiholomorphic, 2), 1e-30)
    )

    maximum_lower = 0.0
    maximum_image = 0.0
    for lift in dual_lifts.T:
        physical_variable = lift[:physical_dimension]
        velocities = (
            lift[
                physical_dimension : physical_dimension + zero_count
            ]
            + 1j * lift[physical_dimension + zero_count :]
        )
        first_operator = sum(
            (
                physical_variable[index] * first_operators[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        image = blaschke_image_jet(
            matrix,
            first_operator,
            np.zeros_like(matrix),
            zeros,
            velocities,
        )[1]
        maximum_lower = max(
            maximum_lower,
            np.linalg.norm(np.tril(image), 2),
        )
        maximum_image = max(maximum_image, np.linalg.norm(image, 2))
    lower_flag_residual = float(
        maximum_lower / max(maximum_image, 1e-30)
    )

    condition_number = float(np.linalg.cond(joint_form))
    checks = (
        joint_dimension == physical_dimension + 2 * zero_count
        and coordinate_rank == expected_active
        and joint_dimension - coordinate_rank == (dimension - 2) ** 2
        and zero_restriction_error < 2e-10
        and backward_residual < 2e-14
        and total_isotropy < 2e-5
        and shape_isotropy < 2e-5
        and cross_isotropy < 2e-5
        and tracking_isotropy < 2e-5
        and optimizer_holomorphic_residual < 2e-8
        and lower_flag_residual < 2e-5
    )
    if not checks:
        raise RuntimeError(
            "shape/tracking Lagrangian audit failed: "
            f"n={dimension}, sample={sample}, "
            f"rank={coordinate_rank}/{expected_active}, "
            f"zero={zero_restriction_error}, "
            f"isotropy={total_isotropy}, "
            f"optimizer={optimizer_holomorphic_residual}, "
            f"lower={lower_flag_residual}, "
            f"backward={backward_residual}, "
            f"condition={condition_number}"
        )

    return ShapeTrackingLagrangianRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        joint_dimension=joint_dimension,
        expected_active_dimension=expected_active,
        observed_coordinate_rank=coordinate_rank,
        coordinate_kernel_dimension=joint_dimension - coordinate_rank,
        tracking_zero_restriction_error=zero_restriction_error,
        joint_form_condition_number=condition_number,
        dual_solve_backward_residual=backward_residual,
        total_bilinear_isotropy_residual=total_isotropy,
        shape_bilinear_isotropy_residual=shape_isotropy,
        shape_tracking_isotropy_residual=cross_isotropy,
        tracking_bilinear_isotropy_residual=tracking_isotropy,
        optimized_tracking_holomorphic_residual=(
            optimizer_holomorphic_residual
        ),
        riesz_lower_flag_residual=lower_flag_residual,
        all_checks_passed=True,
    )


def write_records(
    records: list[ShapeTrackingLagrangianRecord],
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
            "experiments/gau_wu_shape_tracking_lagrangian_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic joint-chart/Lagrangian audit."""

    arguments = parse_args()
    records: list[ShapeTrackingLagrangianRecord] = []
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
                        "coordinate_rank": record.observed_coordinate_rank,
                        "isotropy": (
                            record.total_bilinear_isotropy_residual
                        ),
                        "tracking_holomorphic": (
                            record.optimized_tracking_holomorphic_residual
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
