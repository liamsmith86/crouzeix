#!/usr/bin/env python3
"""Reduce the Gau--Wu boundary loss to first-order frame data."""

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
    frechet_polynomial,
    matrix_polynomial,
    normalized_operator_jet,
    schwarz_map_coefficients,
    support_jet,
)
from gau_wu_physical_loss_rank import build_physical_loss_audit
from theodorsen import hilbert_periodic


@dataclass(frozen=True)
class FirstJetLossRecord:
    """One exact first-jet loss reduction audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    physical_dimension: int
    second_jet_reduction_error: float
    support_curvature_cancellation_error: float
    endpoint_frame_decomposition_error: float
    first_jet_loss_error: float
    minimum_first_base_reserve_eigenvalue: float
    minimum_frame_only_loss_eigenvalue: float
    minimum_loss_eigenvalue: float
    observed_loss_rank: int
    expected_loss_rank: int
    all_checks_passed: bool


def shifted_data(
    dimension: int,
    matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the shift and the two diagonal similarities."""

    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    similarity = np.diag(scaling)
    inverse_similarity = np.diag(1 / scaling)
    shift = inverse_similarity @ matrix @ similarity
    return shift, similarity, inverse_similarity


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> FirstJetLossRecord:
    """Audit the exact cancellation and first-frame loss formula."""

    audit = build_physical_loss_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    endpoint = audit.endpoint_audit
    matrix = endpoint.matrix
    directions = list(endpoint.directions)
    physical_dimension = len(directions)
    zero_count = dimension - 1
    first_operators, second_operators, _ = normalized_operator_jet(
        matrix,
        directions,
        angle_count,
    )
    first_support, second_support, _ = support_jet(
        matrix,
        directions,
        angle_count,
    )
    shift, similarity, inverse_similarity = shifted_data(
        dimension,
        matrix,
    )
    lower_endpoint = np.eye(dimension, dtype=complex)[:, 0]
    upper_endpoint = np.eye(dimension, dtype=complex)[:, -1]
    right_base = np.eye(dimension) - np.outer(
        lower_endpoint,
        lower_endpoint.conj(),
    )
    left_base = np.eye(dimension) - np.outer(
        upper_endpoint,
        upper_endpoint.conj(),
    )

    def components(variable: np.ndarray) -> tuple[float, ...]:
        raw_direction = sum(
            (
                variable[index] * directions[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        first_operator = sum(
            (
                variable[index] * first_operators[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        second_operator = np.einsum(
            "i,j,ijab->ab",
            variable,
            variable,
            second_operators,
        )
        support_velocity = variable @ first_support
        support_curvature = np.einsum(
            "i,j,ijt->t",
            variable,
            variable,
            second_support,
        )
        velocity_coefficients = schwarz_map_coefficients(
            support_velocity
        )
        frequencies = np.fft.fftfreq(
            angle_count,
            d=1 / angle_count,
        )
        support_derivative = np.fft.ifft(
            1j * frequencies * np.fft.fft(support_velocity)
        ).real
        angle_shift = (
            hilbert_periodic(support_velocity) - support_derivative
        )
        reparametrization_operator = matrix_polynomial(
            schwarz_map_coefficients(angle_shift**2 / 2),
            matrix,
        )
        curvature_operator = matrix_polynomial(
            schwarz_map_coefficients(support_curvature),
            matrix,
        )
        first_curvature = (
            -frechet_polynomial(
                velocity_coefficients,
                matrix,
                first_operator,
            )
            + reparametrization_operator
        )
        second_reduction_error = (
            np.linalg.norm(
                second_operator - first_curvature + curvature_operator,
                2,
            )
            / max(np.linalg.norm(second_operator, 2), 1.0)
        )

        first_shift = (
            inverse_similarity @ first_operator @ similarity
        )
        first_curvature_shift = (
            inverse_similarity @ first_curvature @ similarity
        )
        curvature_shift = (
            inverse_similarity @ curvature_operator @ similarity
        )
        zero_coordinates = audit.zero_optimizer @ variable
        zero_velocities = (
            zero_coordinates[:zero_count]
            + 1j * zero_coordinates[zero_count:]
        )
        first_image = blaschke_image_jet(
            matrix,
            first_operator,
            np.zeros_like(matrix),
            endpoint.zeros,
            zero_velocities,
        )[1]
        first_image = (
            inverse_similarity @ first_image @ similarity
        )
        right_image = first_image @ upper_endpoint
        left_image = first_image.conj().T @ lower_endpoint

        right_forcing = (
            first_shift @ right_base @ shift.conj().T
            + shift @ right_base @ first_shift.conj().T
            - np.outer(right_image, lower_endpoint.conj())
            - np.outer(lower_endpoint, right_image.conj())
        )
        right_frame = solve_discrete_lyapunov(
            shift,
            right_forcing,
        )
        left_forcing = (
            first_shift.conj().T @ left_base @ shift
            + shift.conj().T @ left_base @ first_shift
            - np.outer(left_image, upper_endpoint.conj())
            - np.outer(upper_endpoint, left_image.conj())
        )
        left_frame = solve_discrete_lyapunov(
            shift.conj().T,
            left_forcing,
        )

        right_first_defect = -(
            shift.conj().T @ first_shift
            + first_shift.conj().T @ shift
        )
        left_first_defect = -(
            shift @ first_shift.conj().T
            + first_shift @ shift.conj().T
        )
        frame_term = (
            np.trace(right_first_defect @ right_frame).real
            + np.linalg.norm(right_frame @ lower_endpoint) ** 2
            + np.trace(left_first_defect @ left_frame).real
            + np.linalg.norm(left_frame @ upper_endpoint) ** 2
        )

        def base_term(second_shift: np.ndarray) -> float:
            right_second_defect = -(
                shift.conj().T @ second_shift
                + second_shift.conj().T @ shift
                + first_shift.conj().T @ first_shift
            )
            left_second_defect = -(
                shift @ second_shift.conj().T
                + second_shift @ shift.conj().T
                + first_shift @ first_shift.conj().T
            )
            return float(
                np.trace(right_second_defect @ right_base).real
                + np.trace(left_second_defect @ left_base).real
            )

        first_base = base_term(first_curvature_shift)
        full_base = base_term(
            first_curvature_shift - curvature_shift
        )
        support_contribution = full_base - first_base
        support_energy = float(
            variable @ audit.support_gram @ variable
        )
        first_loss = (
            2 * np.linalg.norm(raw_direction) ** 2
            - first_base
            - frame_term
        )
        physical_remainder = float(
            variable @ audit.physical_remainder @ variable
        )
        tracked_loss = float(variable @ audit.loss_form @ variable)
        return (
            float(second_reduction_error),
            support_contribution - 2 * support_energy,
            full_base + frame_term - physical_remainder,
            first_loss - tracked_loss,
            first_base,
            frame_term,
        )

    basis = np.eye(physical_dimension)
    diagonal = np.asarray(
        [components(basis[index]) for index in range(physical_dimension)]
    )
    matrices = [
        np.diag(diagonal[:, index])
        for index in range(1, 6)
    ]
    maximum_second_error = float(np.max(diagonal[:, 0]))
    for row in range(physical_dimension):
        for column in range(row + 1, physical_dimension):
            values = components(basis[row] + basis[column])
            maximum_second_error = max(
                maximum_second_error,
                values[0],
            )
            for index, matrix_value in enumerate(matrices, start=1):
                entry = (
                    values[index]
                    - diagonal[row, index]
                    - diagonal[column, index]
                ) / 2
                matrix_value[row, column] = entry
                matrix_value[column, row] = entry
    scales = (
        max(np.linalg.norm(2 * audit.support_gram, 2), 1.0),
        max(np.linalg.norm(audit.physical_remainder, 2), 1.0),
        max(np.linalg.norm(audit.loss_form, 2), 1.0),
    )
    scalar_errors = tuple(
        float(np.linalg.norm(matrix_value, 2) / scale)
        for matrix_value, scale in zip(matrices[:3], scales, strict=True)
    )
    errors = (maximum_second_error, *scalar_errors)
    first_base_form = matrices[3]
    frame_form = matrices[4]
    first_base_reserve = 2 * np.eye(physical_dimension) - first_base_form
    minimum_first_base_reserve = float(
        np.linalg.eigvalsh(first_base_reserve)[0]
    )
    minimum_frame_only_loss = float(
        np.linalg.eigvalsh(-frame_form)[0]
    )
    loss_eigenvalues = np.linalg.eigvalsh(audit.loss_form)
    singular_values = np.linalg.svd(
        audit.loss_form,
        compute_uv=False,
    )
    threshold = 1e-7 * singular_values[0]
    observed_rank = int(np.sum(singular_values > threshold))
    expected_rank = 6 * dimension - 14
    checks = (
        max(errors) < 2e-9
        and minimum_first_base_reserve > 1e-8
        and minimum_frame_only_loss < -1e-3
        and loss_eigenvalues[0] > -3e-9
        and observed_rank == expected_rank
    )
    if not checks:
        raise RuntimeError(
            "first-jet loss audit failed: "
            f"n={dimension}, sample={sample}, errors={errors}, "
            f"base={minimum_first_base_reserve}, "
            f"frame={minimum_frame_only_loss}, "
            f"loss={loss_eigenvalues[0]}, "
            f"rank={observed_rank}/{expected_rank}"
        )
    return FirstJetLossRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        physical_dimension=physical_dimension,
        second_jet_reduction_error=errors[0],
        support_curvature_cancellation_error=errors[1],
        endpoint_frame_decomposition_error=errors[2],
        first_jet_loss_error=errors[3],
        minimum_first_base_reserve_eigenvalue=(
            minimum_first_base_reserve
        ),
        minimum_frame_only_loss_eigenvalue=minimum_frame_only_loss,
        minimum_loss_eigenvalue=float(loss_eigenvalues[0]),
        observed_loss_rank=observed_rank,
        expected_loss_rank=expected_rank,
        all_checks_passed=True,
    )


def write_records(
    records: list[FirstJetLossRecord],
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
    """Parse a comma-separated list of dimensions."""

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
        default=(4, 5, 6, 7, 8),
    )
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_first_jet_loss_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic first-jet loss audits."""

    arguments = parse_args()
    records: list[FirstJetLossRecord] = []
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
            print(json.dumps(asdict(record), sort_keys=True), flush=True)
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
