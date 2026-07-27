#!/usr/bin/env python3
"""Audit the two-sided endpoint-defect form at finite Gau--Wu models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import (
    blaschke_image_jet,
    joint_hessian,
    normalized_operator_jet,
)
from gau_wu_finite_model import (
    blaschke_at_matrix,
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    random_interior_zeros,
    real_matrix,
    real_vectorize,
)


@dataclass(frozen=True)
class TwoSidedEndpointRecord:
    """One full joint endpoint-defect audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    normal_dimension: int
    joint_dimension: int
    right_negative_index: int
    left_negative_index: int
    minimum_right_eigenvalue: float
    minimum_left_eigenvalue: float
    minimum_two_sided_eigenvalue: float
    maximum_two_sided_eigenvalue: float
    maximum_endpoint_identity_residual: float
    contraction_defect_rank: int
    characteristic_defect_rank: int
    combined_defect_rank: int
    all_checks_passed: bool


@dataclass(frozen=True)
class EndpointFormAudit:
    """Complete endpoint forms and their public audit record."""

    record: TwoSidedEndpointRecord
    matrix: np.ndarray
    zeros: np.ndarray
    directions: tuple[np.ndarray, ...]
    right_form: np.ndarray
    left_form: np.ndarray
    two_sided_form: np.ndarray
    penalty_form: np.ndarray
    sharp_form: np.ndarray


def quadratic_matrix(
    coefficient,
    dimension: int,
) -> np.ndarray:
    """Polarize a real homogeneous quadratic coefficient."""

    basis = np.eye(dimension)
    diagonal = np.asarray(
        [coefficient(basis[index]) for index in range(dimension)]
    )
    matrix = np.diag(diagonal)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            value = coefficient(basis[row] + basis[column])
            entry = (value - diagonal[row] - diagonal[column]) / 2
            matrix[row, column] = entry
            matrix[column, row] = entry
    return matrix


def build_endpoint_form_audit(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> EndpointFormAudit:
    """Build the complete endpoint forms and their audit record."""

    rng = np.random.default_rng(seed)
    interior = random_interior_zeros(dimension, rng)
    matrix = gau_wu_model(interior)
    zeros = extremal_zeros(interior)
    normal, _, _ = normal_basis(interior)
    directions = [
        real_matrix(normal[:, index], dimension)
        for index in range(normal.shape[1])
    ]
    first_operator, second_operator, _ = normalized_operator_jet(
        matrix,
        directions,
        angle_count,
    )
    physical_dimension = len(directions)
    joint_dimension = dimension**2

    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    similarity = np.diag(scaling)
    inverse_similarity = np.diag(1 / scaling)
    shift = inverse_similarity @ matrix @ similarity
    extremal_at_shift = (
        inverse_similarity
        @ blaschke_at_matrix(matrix, zeros)
        @ similarity
    )

    def joint_jet(variable: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
        physical = variable[:physical_dimension]
        zero_real = variable[
            physical_dimension : physical_dimension + dimension - 1
        ]
        zero_imaginary = variable[
            physical_dimension + dimension - 1 :
        ]
        velocities = zero_real + 1j * zero_imaginary
        first = sum(
            (
                physical[index] * first_operator[index]
                for index in range(physical_dimension)
            ),
            np.zeros_like(matrix),
        )
        second = np.einsum(
            "i,j,ijab->ab",
            physical,
            physical,
            second_operator,
        )
        image = blaschke_image_jet(
            matrix,
            first,
            second,
            zeros,
            velocities,
        )
        first_shift = inverse_similarity @ first @ similarity
        image_shift = [
            inverse_similarity @ value @ similarity for value in image
        ]
        return first_shift, image_shift

    def endpoint_coefficients(
        variable: np.ndarray,
    ) -> tuple[float, float, float, float]:
        _, image = joint_jet(variable)
        first = image[1]
        second = image[2]
        transition = np.real(second[0, -1])
        right = -2 * transition - np.linalg.norm(first[:, -1]) ** 2
        left = -2 * transition - np.linalg.norm(first[0, :]) ** 2
        diagonal_penalty = (
            abs(first[-1, -1]) ** 2 + abs(first[0, 0]) ** 2
        ) / 4
        sharp = -(right + left) / 2 - diagonal_penalty
        return (
            float(right),
            float(left),
            float(diagonal_penalty),
            float(sharp),
        )

    right = quadratic_matrix(
        lambda variable: endpoint_coefficients(variable)[0],
        joint_dimension,
    )
    left = quadratic_matrix(
        lambda variable: endpoint_coefficients(variable)[1],
        joint_dimension,
    )
    penalty = quadratic_matrix(
        lambda variable: endpoint_coefficients(variable)[2],
        joint_dimension,
    )
    sharp = quadratic_matrix(
        lambda variable: endpoint_coefficients(variable)[3],
        joint_dimension,
    )
    direct_sharp, _ = joint_hessian(
        matrix,
        first_operator,
        second_operator,
        zeros,
    )
    identity_residual = float(np.max(abs(sharp - direct_sharp)))

    right_eigenvalues = np.linalg.eigvalsh(right)
    left_eigenvalues = np.linalg.eigvalsh(left)
    two_sided = (right + left + right.T + left.T) / 2
    two_sided_eigenvalues = np.linalg.eigvalsh(two_sided)

    contraction_columns: list[np.ndarray] = []
    characteristic_columns: list[np.ndarray] = []
    for index in range(joint_dimension):
        first_shift, image = joint_jet(np.eye(joint_dimension)[index])
        first_image = image[1]
        contraction = (
            shift.conj().T @ first_shift
            + first_shift.conj().T @ shift
        )
        characteristic = (
            first_shift @ extremal_at_shift + shift @ first_image
        )
        contraction_columns.append(real_vectorize(contraction))
        characteristic_columns.append(real_vectorize(characteristic))
    contraction_map = np.asarray(contraction_columns).T
    characteristic_map = np.asarray(characteristic_columns).T
    maximum_singular = max(
        np.linalg.svd(contraction_map, compute_uv=False)[0],
        np.linalg.svd(characteristic_map, compute_uv=False)[0],
    )
    rank_threshold = 2e-9 * maximum_singular
    contraction_rank = int(
        np.linalg.matrix_rank(contraction_map, rank_threshold)
    )
    characteristic_rank = int(
        np.linalg.matrix_rank(characteristic_map, rank_threshold)
    )
    combined_rank = int(
        np.linalg.matrix_rank(
            np.vstack((contraction_map, characteristic_map)),
            rank_threshold,
        )
    )

    checks = (
        joint_dimension
        == physical_dimension + 2 * (dimension - 1)
        and identity_residual < 2e-11
        and two_sided_eigenvalues[0] > -2e-10
        and combined_rank == joint_dimension
        and np.linalg.eigvalsh(penalty)[0] > -2e-10
    )
    if not checks:
        raise RuntimeError(
            "two-sided endpoint audit failed: "
            f"n={dimension}, identity={identity_residual}, "
            f"two_sided={two_sided_eigenvalues[0]}, "
            f"rank={combined_rank}/{joint_dimension}"
        )

    record = TwoSidedEndpointRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        normal_dimension=physical_dimension,
        joint_dimension=joint_dimension,
        right_negative_index=int(np.sum(right_eigenvalues < -1e-8)),
        left_negative_index=int(np.sum(left_eigenvalues < -1e-8)),
        minimum_right_eigenvalue=float(right_eigenvalues[0]),
        minimum_left_eigenvalue=float(left_eigenvalues[0]),
        minimum_two_sided_eigenvalue=float(two_sided_eigenvalues[0]),
        maximum_two_sided_eigenvalue=float(two_sided_eigenvalues[-1]),
        maximum_endpoint_identity_residual=identity_residual,
        contraction_defect_rank=contraction_rank,
        characteristic_defect_rank=characteristic_rank,
        combined_defect_rank=combined_rank,
        all_checks_passed=True,
    )
    return EndpointFormAudit(
        record=record,
        matrix=matrix,
        zeros=zeros,
        directions=tuple(directions),
        right_form=right,
        left_form=left,
        two_sided_form=two_sided,
        penalty_form=penalty,
        sharp_form=sharp,
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> TwoSidedEndpointRecord:
    """Audit the exact endpoint identity and its remaining sign gate."""

    return build_endpoint_form_audit(
        dimension,
        sample,
        seed,
        angle_count,
    ).record


def write_records(
    records: list[TwoSidedEndpointRecord],
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
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_two_sided_endpoint_defect_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested finite-model endpoint audit."""

    arguments = parse_args()
    records: list[TwoSidedEndpointRecord] = []
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
                        "minimum_two_sided_eigenvalue": (
                            record.minimum_two_sided_eigenvalue
                        ),
                        "right_negative_index": (
                            record.right_negative_index
                        ),
                        "left_negative_index": record.left_negative_index,
                        "combined_defect_rank": (
                            record.combined_defect_rank
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
