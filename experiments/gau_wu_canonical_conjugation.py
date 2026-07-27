#!/usr/bin/env python3
"""Audit canonical-conjugation closure of the Gau--Wu endpoint split."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_conformal_dual_polarization import first_blaschke_image
from gau_wu_conformal_shape_phase import conformal_shape_map
from gau_wu_finite_hessian_jet import (
    normalized_operator_jet,
    support_jet,
)
from gau_wu_finite_model import real_vectorize
from gau_wu_similarity_hessian import build_similarity_hessian_audit


@dataclass(frozen=True)
class CanonicalConjugationRecord:
    """One canonical-conjugation covariance audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    canonical_conjugation_error: float
    joint_involution_error: float
    joint_hessian_covariance_error: float
    shape_invariance_error: float
    dual_fixed_point_error: float
    first_image_symmetry_error: float
    endpoint_split_residual: float
    all_checks_passed: bool


def canonical_conjugation(matrix: np.ndarray) -> np.ndarray:
    """Recover the canonical-conjugation matrix with exchanged endpoints."""

    dimension = matrix.shape[0]
    equations: list[np.ndarray] = []
    targets: list[complex] = []

    # The identity A Q = Q A^T expresses C A* C = A for Cx=Q conj(x).
    for row in range(dimension):
        for column in range(dimension):
            equation = np.zeros_like(matrix)
            for index in range(dimension):
                equation[index, column] += matrix[row, index]
                equation[row, index] -= matrix[column, index]
            equations.append(equation.ravel())
            targets.append(0)

    # The model-space conjugation exchanges p=e_0 and q=e_L.
    for column, target_row in ((0, dimension - 1), (dimension - 1, 0)):
        for row in range(dimension):
            equation = np.zeros_like(matrix)
            equation[row, column] = 1
            equations.append(equation.ravel())
            targets.append(complex(row == target_row))

    solution, *_ = np.linalg.lstsq(
        np.asarray(equations),
        np.asarray(targets),
        rcond=None,
    )
    return solution.reshape(dimension, dimension)


def conjugate_operator(
    conjugation: np.ndarray,
    matrix: np.ndarray,
) -> np.ndarray:
    """Evaluate ``C matrix* C`` in ordinary matrix coordinates."""

    return conjugation @ matrix.T @ conjugation.conj().T


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> CanonicalConjugationRecord:
    """Audit one nondegenerate finite Gau--Wu model."""

    audit = build_similarity_hessian_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    endpoint = audit.endpoint_audit
    directions = list(endpoint.directions)
    physical_dimension = len(directions)
    complex_shape_dimension = dimension - 1
    joint_dimension = dimension**2

    conjugation = canonical_conjugation(endpoint.matrix)
    canonical_error = max(
        np.linalg.norm(
            endpoint.matrix @ conjugation
            - conjugation @ endpoint.matrix.T,
            2,
        ),
        np.linalg.norm(conjugation - conjugation.T, 2),
        np.linalg.norm(
            conjugation @ conjugation.conj() - np.eye(dimension),
            2,
        ),
        np.linalg.norm(
            conjugation[:, 0] - np.eye(dimension)[:, -1],
        ),
        np.linalg.norm(
            conjugation[:, -1] - np.eye(dimension)[:, 0],
        ),
    )

    normal = np.stack(
        [real_vectorize(direction) for direction in directions],
        axis=1,
    )
    physical_involution = np.stack(
        [
            normal.T
            @ real_vectorize(conjugate_operator(conjugation, direction))
            for direction in directions
        ],
        axis=1,
    )
    joint_involution = np.block(
        [
            [
                physical_involution,
                np.zeros(
                    (
                        physical_dimension,
                        2 * complex_shape_dimension,
                    )
                ),
            ],
            [
                np.zeros(
                    (
                        2 * complex_shape_dimension,
                        physical_dimension,
                    )
                ),
                np.eye(2 * complex_shape_dimension),
            ],
        ]
    )
    involution_error = float(
        np.linalg.norm(
            joint_involution @ joint_involution
            - np.eye(joint_dimension),
            2,
        )
    )

    joint_form = endpoint.sharp_form
    hessian_covariance = float(
        np.linalg.norm(
            joint_involution.T @ joint_form @ joint_involution
            - joint_form,
            2,
        )
        / np.linalg.norm(joint_form, 2)
    )

    first_support, _, _ = support_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    shape_map = conformal_shape_map(first_support, dimension)
    shape_invariance = float(
        np.linalg.norm(
            shape_map @ physical_involution - shape_map,
            2,
        )
        / np.linalg.norm(shape_map, 2)
    )
    joint_shape_map = np.hstack(
        (
            shape_map,
            np.zeros(
                (
                    2 * complex_shape_dimension,
                    2 * complex_shape_dimension,
                )
            ),
        )
    )
    dual_lifts = np.linalg.solve(joint_form, joint_shape_map.T)
    dual_fixed_point = float(
        np.linalg.norm(
            joint_involution @ dual_lifts - dual_lifts,
            2,
        )
        / np.linalg.norm(dual_lifts, 2)
    )
    analytic_lifts = (
        dual_lifts[:, complex_shape_dimension:]
        + 1j * dual_lifts[:, :complex_shape_dimension]
    )

    first_operators, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    image_symmetry = 0.0
    endpoint_split = 0.0
    for index in range(complex_shape_dimension):
        image = first_blaschke_image(
            analytic_lifts[:, index],
            physical_dimension,
            first_operators,
            endpoint.matrix,
            endpoint.zeros,
        )
        scale = max(np.linalg.norm(image, 2), 1e-30)
        reflected_image = conjugate_operator(conjugation, image)
        image_symmetry = max(
            image_symmetry,
            np.linalg.norm(image - reflected_image, 2) / scale,
        )
        endpoint_split = max(
            endpoint_split,
            abs(image[0, 0] - image[-1, -1]) / scale,
        )

    checks = (
        joint_dimension
        == physical_dimension + 2 * complex_shape_dimension
        and canonical_error < 5e-11
        and involution_error < 5e-11
        and hessian_covariance < 2e-9
        and shape_invariance < 5e-11
        and dual_fixed_point < 5e-10
        and image_symmetry < 5e-10
        and endpoint_split < 5e-11
    )
    if not checks:
        raise RuntimeError(
            "canonical-conjugation audit failed: "
            f"n={dimension}, sample={sample}, "
            f"canonical={canonical_error}, "
            f"involution={involution_error}, "
            f"hessian={hessian_covariance}, "
            f"shape={shape_invariance}, "
            f"dual={dual_fixed_point}, "
            f"image={image_symmetry}, "
            f"split={endpoint_split}"
        )

    return CanonicalConjugationRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        canonical_conjugation_error=float(canonical_error),
        joint_involution_error=involution_error,
        joint_hessian_covariance_error=hessian_covariance,
        shape_invariance_error=shape_invariance,
        dual_fixed_point_error=dual_fixed_point,
        first_image_symmetry_error=float(image_symmetry),
        endpoint_split_residual=float(endpoint_split),
        all_checks_passed=True,
    )


def write_records(
    records: list[CanonicalConjugationRecord],
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
            "experiments/gau_wu_canonical_conjugation_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic canonical-conjugation audit."""

    arguments = parse_args()
    records: list[CanonicalConjugationRecord] = []
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
                        "dual_fixed": record.dual_fixed_point_error,
                        "image_symmetry": (
                            record.first_image_symmetry_error
                        ),
                        "endpoint_split": record.endpoint_split_residual,
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
