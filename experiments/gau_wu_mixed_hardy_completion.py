#!/usr/bin/env python3
"""Audit the exact mixed Hardy completion at finite Gau--Wu models."""

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
from gau_wu_finite_model import (
    blaschke_at_matrix,
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    random_interior_zeros,
    real_matrix,
)


@dataclass(frozen=True)
class MixedHardyCompletionRecord:
    """One complete mixed-coordinate audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    probe_count: int
    maximum_zero_gram_residual: float
    maximum_tangent_output_residual: float
    maximum_mixed_formula_residual: float
    maximum_completion_residual: float
    all_checks_passed: bool


def matrix_blaschke_factor(
    matrix: np.ndarray,
    zero: complex,
) -> np.ndarray:
    """Evaluate one normalized Blaschke factor at a matrix."""

    identity = np.eye(matrix.shape[0])
    return (matrix - zero * identity) @ np.linalg.inv(
        identity - np.conj(zero) * matrix
    )


def hardy_coordinate_at_matrix(
    matrix: np.ndarray,
    zeros: np.ndarray,
    velocities: np.ndarray,
) -> np.ndarray:
    """Evaluate k=sum conj(v_j)/(1-conj(a_j)z) at a matrix."""

    identity = np.eye(matrix.shape[0])
    return sum(
        (
            np.conj(velocity)
            * np.linalg.inv(identity - np.conj(zero) * matrix)
            for zero, velocity in zip(zeros, velocities, strict=True)
        ),
        np.zeros_like(matrix),
    )


def conjugate_coordinate_at_matrix(
    matrix: np.ndarray,
    zeros: np.ndarray,
    velocities: np.ndarray,
) -> np.ndarray:
    """Evaluate J_f k=sum v_j f/(z-a_j) at a matrix."""

    identity = np.eye(matrix.shape[0])
    factors = [
        matrix_blaschke_factor(matrix, zero) for zero in zeros
    ]
    value = np.zeros_like(matrix)
    for index, (zero, velocity) in enumerate(
        zip(zeros, velocities, strict=True)
    ):
        quotient = np.linalg.inv(identity - np.conj(zero) * matrix)
        for factor_index, factor in enumerate(factors):
            if factor_index != index:
                quotient = quotient @ factor
        value += velocity * quotient
    return value


def inner_tangent_at_matrix(
    matrix: np.ndarray,
    zeros: np.ndarray,
    velocities: np.ndarray,
) -> np.ndarray:
    """Evaluate h_k=(zf)k-J_f k at a matrix."""

    inner = blaschke_at_matrix(matrix, zeros)
    coordinate = hardy_coordinate_at_matrix(matrix, zeros, velocities)
    conjugate = conjugate_coordinate_at_matrix(
        matrix,
        zeros,
        velocities,
    )
    return matrix @ inner @ coordinate - conjugate


def block_frechet(
    function,
    matrix: np.ndarray,
    direction: np.ndarray,
) -> np.ndarray:
    """Evaluate a rational Fréchet derivative by a block matrix."""

    dimension = matrix.shape[0]
    block = np.block(
        [
            [matrix, direction],
            [np.zeros_like(matrix), matrix],
        ]
    )
    return function(block)[:dimension, dimension:]


def endpoint_sum(image: list[np.ndarray]) -> float:
    """Return the second coefficient of the two endpoint defects."""

    first = image[1]
    second = image[2]
    return float(
        -4 * np.real(second[0, -1])
        - np.linalg.norm(first[:, -1]) ** 2
        - np.linalg.norm(first[0, :]) ** 2
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
    probe_count: int,
) -> MixedHardyCompletionRecord:
    """Audit the Hardy coordinate and mixed completion on one model."""

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

    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    similarity = np.diag(scaling)
    inverse_similarity = np.diag(1 / scaling)
    shift = inverse_similarity @ matrix @ similarity
    p = np.eye(dimension)[:, 0]
    q = np.eye(dimension)[:, -1]

    maximum_zero_gram_residual = 0.0
    maximum_tangent_output_residual = 0.0
    maximum_mixed_formula_residual = 0.0
    maximum_completion_residual = 0.0

    for _ in range(probe_count):
        physical = rng.normal(size=len(directions))
        velocities = rng.normal(size=len(zeros)) + 1j * rng.normal(
            size=len(zeros)
        )
        first = sum(
            (
                physical[index] * first_operator[index]
                for index in range(len(directions))
            ),
            np.zeros_like(matrix),
        )
        second = np.einsum(
            "i,j,ijab->ab",
            physical,
            physical,
            second_operator,
        )
        zero = np.zeros_like(matrix)
        zero_velocities = np.zeros_like(velocities)

        physical_image = blaschke_image_jet(
            matrix,
            first,
            second,
            zeros,
            zero_velocities,
        )
        zero_image = blaschke_image_jet(
            matrix,
            zero,
            zero,
            zeros,
            velocities,
        )
        joint_image = blaschke_image_jet(
            matrix,
            first,
            second,
            zeros,
            velocities,
        )
        physical_shift_image = [
            inverse_similarity @ value @ similarity
            for value in physical_image
        ]
        zero_shift_image = [
            inverse_similarity @ value @ similarity
            for value in zero_image
        ]
        joint_shift_image = [
            inverse_similarity @ value @ similarity
            for value in joint_image
        ]
        direction = inverse_similarity @ first @ similarity

        coordinate = hardy_coordinate_at_matrix(
            shift,
            zeros,
            velocities,
        )
        conjugate = conjugate_coordinate_at_matrix(
            shift,
            zeros,
            velocities,
        )
        coordinate_vector = coordinate @ q
        conjugate_vector = conjugate @ q
        gram_value = sum(
            (
                velocity_i * np.conj(velocity_j)
                / (1 - zero_i * np.conj(zero_j))
                for zero_i, velocity_i in zip(
                    zeros,
                    velocities,
                    strict=True,
                )
                for zero_j, velocity_j in zip(
                    zeros,
                    velocities,
                    strict=True,
                )
            )
        )
        maximum_zero_gram_residual = max(
            maximum_zero_gram_residual,
            abs(np.linalg.norm(coordinate_vector) ** 2 - gram_value),
            abs(
                endpoint_sum(zero_shift_image)
                - 2 * np.linalg.norm(coordinate_vector) ** 2
            ),
        )
        maximum_tangent_output_residual = max(
            maximum_tangent_output_residual,
            np.linalg.norm(
                zero_shift_image[1] @ q + conjugate_vector
            ),
            np.linalg.norm(
                zero_shift_image[1].conj().T @ p
                + shift @ coordinate_vector
            ),
        )

        first_image = physical_shift_image[1]
        u_value = first_image @ q
        v_value = first_image.conj().T @ p
        phi_frechet = (
            direction
            @ (inverse_similarity @ physical_image[0] @ similarity)
            + shift @ first_image
        )
        tangent_frechet = block_frechet(
            lambda value: inner_tangent_at_matrix(
                value,
                zeros,
                velocities,
            ),
            shift,
            direction,
        )
        transition = np.vdot(p, tangent_frechet @ q)
        transition_from_product = np.vdot(
            p,
            phi_frechet @ coordinate_vector
            - block_frechet(
                lambda value: conjugate_coordinate_at_matrix(
                    value,
                    zeros,
                    velocities,
                ),
                shift,
                direction,
            )
            @ q,
        )
        mixed_expected = (
            -4 * np.real(transition)
            + 2
            * np.real(
                np.vdot(u_value, conjugate_vector)
                + np.vdot(v_value, shift @ coordinate_vector)
            )
        )
        mixed_actual = (
            endpoint_sum(joint_shift_image)
            - endpoint_sum(physical_shift_image)
            - endpoint_sum(zero_shift_image)
        )
        maximum_mixed_formula_residual = max(
            maximum_mixed_formula_residual,
            abs(transition - transition_from_product),
            abs(mixed_actual - mixed_expected),
        )
        completed = (
            endpoint_sum(physical_shift_image)
            + 2 * np.linalg.norm(coordinate_vector) ** 2
            + mixed_expected
        )
        maximum_completion_residual = max(
            maximum_completion_residual,
            abs(endpoint_sum(joint_shift_image) - completed),
        )

    checks = (
        maximum_zero_gram_residual < 2e-10
        and maximum_tangent_output_residual < 2e-10
        and maximum_mixed_formula_residual < 3e-10
        and maximum_completion_residual < 3e-10
    )
    if not checks:
        raise RuntimeError(
            "mixed Hardy completion failed: "
            f"n={dimension}, zero={maximum_zero_gram_residual}, "
            f"tangent={maximum_tangent_output_residual}, "
            f"mixed={maximum_mixed_formula_residual}, "
            f"completion={maximum_completion_residual}"
        )
    return MixedHardyCompletionRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        probe_count=probe_count,
        maximum_zero_gram_residual=float(maximum_zero_gram_residual),
        maximum_tangent_output_residual=float(
            maximum_tangent_output_residual
        ),
        maximum_mixed_formula_residual=float(
            maximum_mixed_formula_residual
        ),
        maximum_completion_residual=float(maximum_completion_residual),
        all_checks_passed=True,
    )


def write_records(
    records: list[MixedHardyCompletionRecord],
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
    parser.add_argument("--probes", type=int, default=4)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_mixed_hardy_completion_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested mixed Hardy audits."""

    arguments = parse_args()
    records: list[MixedHardyCompletionRecord] = []
    for dimension in arguments.dimensions:
        for sample in range(arguments.samples):
            seed = arguments.seed + 1009 * dimension + sample
            record = audit_model(
                dimension,
                sample,
                seed,
                arguments.angle_count,
                arguments.probes,
            )
            records.append(record)
            print(
                json.dumps(
                    {
                        "dimension": dimension,
                        "sample": sample,
                        "maximum_mixed_formula_residual": (
                            record.maximum_mixed_formula_residual
                        ),
                        "maximum_tangent_output_residual": (
                            record.maximum_tangent_output_residual
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
