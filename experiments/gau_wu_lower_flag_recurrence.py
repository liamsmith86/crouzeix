#!/usr/bin/env python3
"""Audit the lower-Schur-flag recurrence for Gau--Wu first responses."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_conformal_dual_polarization import first_blaschke_image
from gau_wu_finite_hessian_jet import normalized_operator_jet
from gau_wu_similarity_hessian import build_similarity_hessian_audit


@dataclass(frozen=True)
class LowerFlagRecurrenceRecord:
    """One first-response commutator and flag-generation audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    joint_dimension: int
    expected_lower_response_rank: int
    observed_lower_response_rank: int
    observed_generator_rank: int
    maximum_commutator_recurrence_residual: float
    lower_row_space_generator_residual: float
    all_checks_passed: bool


def numerical_rank(matrix: np.ndarray) -> int:
    """Return a conservative relative numerical rank."""

    singular_values = np.linalg.svd(matrix, compute_uv=False)
    if not len(singular_values) or singular_values[0] == 0:
        return 0
    return int(np.sum(singular_values > 1e-8 * singular_values[0]))


def row_space_residual(
    matrix: np.ndarray,
    generators: np.ndarray,
) -> float:
    """Measure the part of ``matrix`` outside the generator row space."""

    generator_gram = generators @ generators.conj().T
    projection = (
        matrix
        @ generators.conj().T
        @ np.linalg.pinv(generator_gram)
        @ generators
    )
    return float(
        np.linalg.norm(matrix - projection, 2)
        / np.linalg.norm(matrix, 2)
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> LowerFlagRecurrenceRecord:
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
    zero_count = len(endpoint.zeros)
    joint_dimension = dimension**2
    first_operators, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )

    lower_indices = np.tril_indices(dimension)
    images: list[np.ndarray] = []
    maximum_commutator_residual = 0.0
    identity = np.eye(joint_dimension)
    base_image = np.zeros_like(endpoint.matrix)
    base_image[0, -1] = 2
    for index in range(joint_dimension):
        variable = identity[index]
        image = first_blaschke_image(
            variable,
            physical_dimension,
            first_operators,
            endpoint.matrix,
            endpoint.zeros,
        )
        images.append(image)

        physical = variable[:physical_dimension]
        first_operator = sum(
            (
                physical[direction_index]
                * first_operators[direction_index]
                for direction_index in range(physical_dimension)
            ),
            np.zeros_like(endpoint.matrix),
        )
        commutator = (
            endpoint.matrix @ image
            - image @ endpoint.matrix
            - base_image @ first_operator
            + first_operator @ base_image
        )
        maximum_commutator_residual = max(
            maximum_commutator_residual,
            np.linalg.norm(commutator, 2),
        )

    lower_map = np.stack(
        [image[lower_indices] for image in images],
        axis=1,
    )
    generator_pairs = [
        *((index, index) for index in range(dimension)),
        (dimension - 1, 0),
    ]
    generator_map = np.stack(
        [
            np.asarray([image[row, column] for image in images])
            for row, column in generator_pairs
        ]
    )

    expected_rank = dimension + 1
    lower_rank = numerical_rank(lower_map)
    generator_rank = numerical_rank(generator_map)
    generator_residual = row_space_residual(
        lower_map,
        generator_map,
    )
    checks = (
        zero_count == dimension - 1
        and lower_rank == expected_rank
        and generator_rank == expected_rank
        and maximum_commutator_residual < 2e-10
        and generator_residual < 2e-8
    )
    if not checks:
        raise RuntimeError(
            "lower-flag recurrence audit failed: "
            f"n={dimension}, sample={sample}, "
            f"lower_rank={lower_rank}/{expected_rank}, "
            f"generator_rank={generator_rank}/{expected_rank}, "
            f"commutator={maximum_commutator_residual}, "
            f"row_space={generator_residual}"
        )

    return LowerFlagRecurrenceRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        joint_dimension=joint_dimension,
        expected_lower_response_rank=expected_rank,
        observed_lower_response_rank=lower_rank,
        observed_generator_rank=generator_rank,
        maximum_commutator_recurrence_residual=float(
            maximum_commutator_residual
        ),
        lower_row_space_generator_residual=generator_residual,
        all_checks_passed=True,
    )


def write_records(
    records: list[LowerFlagRecurrenceRecord],
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
            "experiments/gau_wu_lower_flag_recurrence_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic lower-flag recurrence audit."""

    arguments = parse_args()
    records: list[LowerFlagRecurrenceRecord] = []
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
                        "lower_rank": (
                            record.observed_lower_response_rank
                        ),
                        "commutator": (
                            record.maximum_commutator_recurrence_residual
                        ),
                        "row_space": (
                            record.lower_row_space_generator_residual
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
