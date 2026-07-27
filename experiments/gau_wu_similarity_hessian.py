#!/usr/bin/env python3
"""Reduce and audit the second-order similarity SDP at Gau--Wu models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_discrete_lyapunov

from gau_wu_finite_hessian_jet import (
    blaschke_image_jet,
    normalized_operator_jet,
)
from gau_wu_two_sided_endpoint_defect import (
    EndpointFormAudit,
    build_endpoint_form_audit,
    quadratic_matrix,
)
from general_similarity_second_order_probe import solve_second_order_sdp


@dataclass(frozen=True)
class SimilarityHessianRecord:
    """One complete reduced similarity-Hessian audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    physical_dimension: int
    metric_boundary_dimension: int
    base_metric_error: float
    base_contraction_error: float
    dual_balance_error: float
    maximum_first_stationarity_error: float
    maximum_extension_residual: float
    universal_metric_hessian_error: float
    gap_factorization_error: float
    residual_map_minimum_singular_value: float
    scalar_osculation_error: float
    minimum_similarity_eigenvalue: float
    maximum_similarity_eigenvalue: float
    maximum_full_sdp_error: float
    sign_classification: str
    all_checks_passed: bool


@dataclass(frozen=True)
class SimilarityHessianAudit:
    """Reduced joint and optimized forms for one Gau--Wu model."""

    record: SimilarityHessianRecord
    endpoint_audit: EndpointFormAudit
    joint_metric_form: np.ndarray
    optimized_similarity_form: np.ndarray
    optimized_scalar_form: np.ndarray


def gau_wu_metric(dimension: int) -> np.ndarray:
    """Return the sharp condition-square-four metric."""

    return np.diag([1.0, *([2.0] * (dimension - 2)), 4.0])


def gau_wu_dual_weight(dimension: int) -> np.ndarray:
    """Return the embedded positive Stein dual weight."""

    return np.diag([0.0, *([2.0] * (dimension - 2)), 1.0])


def homogeneous_metric_hessian(dimension: int) -> np.ndarray:
    """Return the universal free-boundary Hessian in real coordinates."""

    one_orientation = [8.0] * (dimension - 2) + [8.0 / 3.0]
    return np.diag([*one_orientation, *one_orientation])


def metric_extension(
    matrix: np.ndarray,
    forcing: np.ndarray,
    boundary_row: np.ndarray,
) -> np.ndarray:
    """Solve the first-order active Stein equation from its boundary row."""

    dimension = matrix.shape[0]
    tail = matrix[1:, 1:]
    leading_row = matrix[0, 1:]
    propagated_row = boundary_row @ tail
    stein_forcing = (
        forcing[1:, 1:]
        + np.outer(leading_row.conj(), propagated_row)
        + np.outer(propagated_row.conj(), leading_row)
    )
    tail_metric = solve_discrete_lyapunov(
        tail.conj().T,
        stein_forcing,
    )
    extension = np.empty((dimension, dimension), dtype=complex)
    extension[0, 0] = 0
    extension[0, 1:] = boundary_row
    extension[1:, 0] = boundary_row.conj()
    extension[1:, 1:] = (tail_metric + tail_metric.conj().T) / 2
    return extension


def reduced_metric_coefficient(
    matrix: np.ndarray,
    metric: np.ndarray,
    dual_weight: np.ndarray,
    first_operator: np.ndarray,
    second_operator: np.ndarray,
    boundary_row: np.ndarray,
) -> tuple[float, float]:
    """Evaluate the exact scalar reduction of the second-order SDP."""

    dimension = matrix.shape[0]
    identity = np.eye(dimension)
    tail_indices = np.arange(1, dimension)
    upper_range = np.arange(dimension - 1)
    contraction_base = metric - matrix.conj().T @ metric @ matrix
    first_forcing = (
        first_operator.conj().T @ metric @ matrix
        + matrix.conj().T @ metric @ first_operator
    )
    metric_tangent = metric_extension(
        matrix,
        first_forcing,
        boundary_row,
    )
    first_defect = (
        metric_tangent
        - matrix.conj().T @ metric_tangent @ matrix
        - first_forcing
    )

    lower_range = (metric - identity)[np.ix_(tail_indices, tail_indices)]
    upper_metric = (4 * identity - metric)[
        np.ix_(upper_range, upper_range)
    ]
    lower_penalty = np.vdot(
        metric_tangent[tail_indices, 0],
        np.linalg.solve(lower_range, metric_tangent[tail_indices, 0]),
    ).real
    upper_penalty = np.vdot(
        metric_tangent[upper_range, -1],
        np.linalg.solve(
            upper_metric,
            metric_tangent[upper_range, -1],
        ),
    ).real

    second_forcing = (
        second_operator.conj().T @ metric @ matrix
        + matrix.conj().T @ metric @ second_operator
        + first_operator.conj().T @ metric @ first_operator
        + first_operator.conj().T @ metric_tangent @ matrix
        + matrix.conj().T @ metric_tangent @ first_operator
    )
    kernel_to_range = first_defect[tail_indices, 0]
    active_weight = dual_weight[np.ix_(tail_indices, tail_indices)]
    contraction_penalty = np.vdot(
        kernel_to_range,
        active_weight @ kernel_to_range,
    ).real / contraction_base[0, 0].real
    weighted_forcing = np.trace(
        active_weight
        @ second_forcing[np.ix_(tail_indices, tail_indices)]
    ).real
    value = (
        4 * lower_penalty
        + upper_penalty
        + contraction_penalty
        + weighted_forcing
    )
    extension_residual = max(
        abs(metric_tangent[-1, -1]),
        np.linalg.norm(
            first_defect[np.ix_(tail_indices, tail_indices)],
            2,
        ),
    )
    return float(value), float(extension_residual)


def optimized_scalar_form(
    endpoint_audit: EndpointFormAudit,
) -> np.ndarray:
    """Maximize the scalar Hessian over all Blaschke-zero velocities."""

    physical_dimension = len(endpoint_audit.directions)
    joint_form = endpoint_audit.sharp_form
    physical = joint_form[:physical_dimension, :physical_dimension]
    mixed = joint_form[:physical_dimension, physical_dimension:]
    zeros = joint_form[physical_dimension:, physical_dimension:]
    return physical - mixed @ np.linalg.solve(zeros, mixed.T)


def boundary_residual_map(
    endpoint_audit: EndpointFormAudit,
    first_operators: list[np.ndarray],
    metric: np.ndarray,
    joint_dimension: int,
) -> np.ndarray:
    """Return the real two-endpoint matching residual map."""

    matrix = endpoint_audit.matrix
    zeros = endpoint_audit.zeros
    dimension = matrix.shape[0]
    physical_dimension = len(endpoint_audit.directions)
    boundary_dimension = 2 * (dimension - 1)
    lower_indices = np.arange(1, dimension)
    upper_indices = np.arange(dimension - 1)

    def residual(variable: np.ndarray) -> np.ndarray:
        physical = variable[:physical_dimension]
        boundary = variable[
            physical_dimension : physical_dimension + boundary_dimension
        ]
        zero_velocity = variable[
            physical_dimension + boundary_dimension :
        ]
        first_operator = sum(
            (
                physical[index] * first_operators[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        first_forcing = (
            first_operator.conj().T @ metric @ matrix
            + matrix.conj().T @ metric @ first_operator
        )
        boundary_row = (
            boundary[: dimension - 1]
            + 1j * boundary[dimension - 1 :]
        )
        metric_tangent = metric_extension(
            matrix,
            first_forcing,
            boundary_row,
        )
        lower_motion = -np.linalg.solve(
            (metric - np.eye(dimension))[
                np.ix_(lower_indices, lower_indices)
            ],
            metric_tangent[lower_indices, 0],
        )
        upper_motion = np.linalg.solve(
            (4 * np.eye(dimension) - metric)[
                np.ix_(upper_indices, upper_indices)
            ],
            metric_tangent[upper_indices, -1],
        )
        velocities = (
            zero_velocity[: dimension - 1]
            + 1j * zero_velocity[dimension - 1 :]
        )
        first_image = blaschke_image_jet(
            matrix,
            first_operator,
            np.zeros_like(matrix),
            zeros,
            velocities,
        )[1]
        lower_residual = (
            first_image[lower_indices, -1] - 2 * lower_motion
        )
        upper_residual = (
            first_image.conj().T[upper_indices, 0] - 2 * upper_motion
        )
        return np.concatenate(
            (
                lower_residual.real,
                lower_residual.imag,
                upper_residual.real,
                upper_residual.imag,
            )
        )

    return np.stack(
        [residual(np.eye(joint_dimension)[index]) for index in range(joint_dimension)],
        axis=1,
    )


def residual_metric(dimension: int) -> np.ndarray:
    """Return the exact metric in the two endpoint-matching squares."""

    lower = [1.0] * (dimension - 2) + [3.0]
    upper = [3.0] + [1.0] * (dimension - 2)
    return np.diag([*lower, *lower, *upper, *upper])


def build_similarity_hessian_audit(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> SimilarityHessianAudit:
    """Build the reduced similarity form and all independent checks."""

    endpoint_audit = build_endpoint_form_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    matrix = endpoint_audit.matrix
    directions = list(endpoint_audit.directions)
    physical_dimension = len(directions)
    boundary_dimension = 2 * (dimension - 1)
    metric = gau_wu_metric(dimension)
    dual_weight = gau_wu_dual_weight(dimension)
    first_operators, second_operators, _ = normalized_operator_jet(
        matrix,
        directions,
        angle_count,
    )

    def coefficient(variable: np.ndarray) -> float:
        physical = variable[:physical_dimension]
        boundary_real = variable[
            physical_dimension : physical_dimension + dimension - 1
        ]
        boundary_imaginary = variable[
            physical_dimension + dimension - 1 :
        ]
        boundary_row = boundary_real + 1j * boundary_imaginary
        first_operator = sum(
            (
                physical[index] * first_operators[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        second_operator = np.einsum(
            "i,j,ijab->ab",
            physical,
            physical,
            second_operators,
        )
        return reduced_metric_coefficient(
            matrix,
            metric,
            dual_weight,
            first_operator,
            second_operator,
            boundary_row,
        )[0]

    joint_dimension = physical_dimension + boundary_dimension
    joint_metric = quadratic_matrix(coefficient, joint_dimension)
    physical_block = joint_metric[
        :physical_dimension,
        :physical_dimension,
    ]
    mixed_block = joint_metric[:physical_dimension, physical_dimension:]
    boundary_block = joint_metric[
        physical_dimension:,
        physical_dimension:,
    ]
    similarity_form = (
        physical_block
        - mixed_block @ np.linalg.solve(boundary_block, mixed_block.T)
    )
    similarity_form = (similarity_form + similarity_form.T) / 2
    scalar_form = optimized_scalar_form(endpoint_audit)
    scalar_joint = endpoint_audit.sharp_form

    combined_dimension = (
        physical_dimension + 2 * boundary_dimension
    )
    gap_form = np.zeros((combined_dimension, combined_dimension))
    metric_indices = np.arange(physical_dimension + boundary_dimension)
    scalar_indices = np.concatenate(
        (
            np.arange(physical_dimension),
            np.arange(
                physical_dimension + boundary_dimension,
                combined_dimension,
            ),
        )
    )
    gap_form[np.ix_(metric_indices, metric_indices)] += joint_metric
    gap_form[np.ix_(scalar_indices, scalar_indices)] -= 4 * scalar_joint
    residual_map = boundary_residual_map(
        endpoint_audit,
        first_operators,
        metric,
        combined_dimension,
    )
    residual_singular_values = np.linalg.svd(
        residual_map,
        compute_uv=False,
    )
    gap_factorization_error = float(
        np.linalg.norm(
            gap_form
            - residual_map.T
            @ residual_metric(dimension)
            @ residual_map,
            2,
        )
    )
    residual_minimum = float(residual_singular_values[-1])

    maximum_stationarity_error = 0.0
    for first_operator in first_operators:
        first_forcing = (
            first_operator.conj().T @ metric @ matrix
            + matrix.conj().T @ metric @ first_operator
        )
        maximum_stationarity_error = max(
            maximum_stationarity_error,
            abs(np.trace(dual_weight @ first_forcing)),
        )

    rng = np.random.default_rng(seed + 7919)
    maximum_extension_residual = 0.0
    maximum_full_sdp_error = 0.0
    check_directions: list[np.ndarray] = []
    for _ in range(2):
        direction = rng.normal(size=physical_dimension)
        check_directions.append(direction / np.linalg.norm(direction))
    soft_values, soft_vectors = np.linalg.eigh(similarity_form)
    check_directions.append(soft_vectors[:, -1])
    for physical in check_directions:
        boundary = -np.linalg.solve(
            boundary_block,
            mixed_block.T @ physical,
        )
        _, extension_residual = reduced_metric_coefficient(
            matrix,
            metric,
            dual_weight,
            sum(
                (
                    physical[index] * first_operators[index]
                    for index in range(physical_dimension)
                ),
                np.zeros_like(matrix),
            ),
            np.einsum(
                "i,j,ijab->ab",
                physical,
                physical,
                second_operators,
            ),
            boundary[: dimension - 1]
            + 1j * boundary[dimension - 1 :],
        )
        maximum_extension_residual = max(
            maximum_extension_residual,
            extension_residual,
        )
        full_solution = solve_second_order_sdp(
            matrix,
            sum(
                (
                    physical[index] * first_operators[index]
                    for index in range(physical_dimension)
                ),
                np.zeros_like(matrix),
            ),
            np.einsum(
                "i,j,ijab->ab",
                physical,
                physical,
                second_operators,
            ),
            metric,
            "CLARABEL",
        )
        reduced_value = float(physical @ similarity_form @ physical)
        maximum_full_sdp_error = max(
            maximum_full_sdp_error,
            abs(full_solution.value - reduced_value),
        )

    base_contraction = metric - matrix.conj().T @ metric @ matrix
    expected_contraction = np.zeros_like(matrix)
    expected_contraction[0, 0] = 1
    expected_dual_balance = np.zeros_like(matrix)
    expected_dual_balance[-1, -1] = 1
    expected_dual_balance[0, 0] = -4
    base_metric_error = max(
        abs(metric[0, 0] - 1),
        abs(metric[-1, -1] - 4),
        max(0.0, 1 - np.linalg.eigvalsh(metric)[0]),
        max(0.0, np.linalg.eigvalsh(metric)[-1] - 4),
    )
    base_contraction_error = float(
        np.linalg.norm(base_contraction - expected_contraction, 2)
    )
    dual_balance_error = float(
        np.linalg.norm(
            dual_weight
            - matrix @ dual_weight @ matrix.conj().T
            - expected_dual_balance,
            2,
        )
    )
    universal_error = float(
        np.linalg.norm(
            boundary_block - homogeneous_metric_hessian(dimension),
            2,
        )
    )
    osculation_error = float(
        np.linalg.norm(similarity_form - 4 * scalar_form, 2)
    )
    minimum_eigenvalue = float(soft_values[0])
    maximum_eigenvalue = float(soft_values[-1])
    sign_classification = (
        "strictly_negative"
        if maximum_eigenvalue < -2e-7
        else "numerically_unresolved_negative"
        if maximum_eigenvalue <= 2e-9
        else "positive_mode"
    )
    checks = (
        physical_dimension == (dimension - 1) ** 2 + 1
        and base_metric_error < 1e-14
        and base_contraction_error < 2e-13
        and dual_balance_error < 2e-13
        and maximum_stationarity_error < 2e-11
        and maximum_extension_residual < 2e-10
        and universal_error < 2e-10
        and gap_factorization_error < 2e-9
        and residual_minimum > 1e-8
        and osculation_error < 1e-7
        and maximum_full_sdp_error < 2e-6
        and maximum_eigenvalue < 2e-9
    )
    if not checks:
        raise RuntimeError(
            "Gau--Wu similarity-Hessian audit failed: "
            f"n={dimension}, metric={base_metric_error}, "
            f"contraction={base_contraction_error}, "
            f"dual={dual_balance_error}, "
            f"stationarity={maximum_stationarity_error}, "
            f"extension={maximum_extension_residual}, "
            f"universal={universal_error}, "
            f"gap={gap_factorization_error}, "
            f"residual_min={residual_minimum}, "
            f"osculation={osculation_error}, "
            f"SDP={maximum_full_sdp_error}, "
            f"lambda_max={maximum_eigenvalue}"
        )

    record = SimilarityHessianRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        physical_dimension=physical_dimension,
        metric_boundary_dimension=boundary_dimension,
        base_metric_error=float(base_metric_error),
        base_contraction_error=base_contraction_error,
        dual_balance_error=dual_balance_error,
        maximum_first_stationarity_error=maximum_stationarity_error,
        maximum_extension_residual=maximum_extension_residual,
        universal_metric_hessian_error=universal_error,
        gap_factorization_error=gap_factorization_error,
        residual_map_minimum_singular_value=residual_minimum,
        scalar_osculation_error=osculation_error,
        minimum_similarity_eigenvalue=minimum_eigenvalue,
        maximum_similarity_eigenvalue=maximum_eigenvalue,
        maximum_full_sdp_error=maximum_full_sdp_error,
        sign_classification=sign_classification,
        all_checks_passed=True,
    )
    return SimilarityHessianAudit(
        record=record,
        endpoint_audit=endpoint_audit,
        joint_metric_form=joint_metric,
        optimized_similarity_form=similarity_form,
        optimized_scalar_form=scalar_form,
    )


def write_records(
    records: list[SimilarityHessianRecord],
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
            "experiments/gau_wu_similarity_hessian_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested similarity-Hessian audits."""

    arguments = parse_args()
    records: list[SimilarityHessianRecord] = []
    for dimension in arguments.dimensions:
        for sample in range(arguments.samples):
            seed = arguments.seed + 1009 * dimension + sample
            audit = build_similarity_hessian_audit(
                dimension,
                sample,
                seed,
                arguments.angle_count,
            )
            records.append(audit.record)
            print(
                json.dumps(
                    {
                        "dimension": dimension,
                        "sample": sample,
                        "maximum_similarity_eigenvalue": (
                            audit.record.maximum_similarity_eigenvalue
                        ),
                        "scalar_osculation_error": (
                            audit.record.scalar_osculation_error
                        ),
                        "maximum_full_sdp_error": (
                            audit.record.maximum_full_sdp_error
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
