#!/usr/bin/env python3
"""Audit the dual Hardy polarization behind Gau--Wu phase covariance."""

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
from gau_wu_similarity_hessian import build_similarity_hessian_audit


@dataclass(frozen=True)
class ConformalDualPolarizationRecord:
    """One joint-Hessian dual-polarization audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    joint_dimension: int
    complex_shape_dimension: int
    joint_form_condition_number: float
    dual_solve_backward_residual: float
    phase_covariance_residual: float
    sharp_isotropy_residual: float
    right_isotropy_residual: float
    left_isotropy_residual: float
    lower_triangular_response_residual: float
    all_checks_passed: bool


def first_blaschke_image(
    variable: np.ndarray,
    physical_dimension: int,
    first_operators: list[np.ndarray],
    matrix: np.ndarray,
    zeros: np.ndarray,
) -> np.ndarray:
    """Return the first image for one joint physical/zero direction."""

    zero_count = len(zeros)
    physical = variable[:physical_dimension]
    zero_real = variable[
        physical_dimension : physical_dimension + zero_count
    ]
    zero_imaginary = variable[physical_dimension + zero_count :]
    first_operator = sum(
        (
            physical[index] * first_operators[index]
            for index in range(physical_dimension)
        ),
        np.zeros_like(matrix),
    )
    velocities = zero_real + 1j * zero_imaginary
    return blaschke_image_jet(
        matrix,
        first_operator,
        np.zeros_like(matrix),
        zeros,
        velocities,
    )[1]


def complex_phase_residual(
    real_form: np.ndarray,
    complex_dimension: int,
) -> float:
    """Measure commutation with simultaneous multiplication by ``i``."""

    zero = np.zeros((complex_dimension, complex_dimension))
    identity = np.eye(complex_dimension)
    phase = np.block([[zero, -identity], [identity, zero]])
    return float(
        np.linalg.norm(
            real_form - phase.T @ real_form @ phase,
            2,
        )
        / np.linalg.norm(real_form, 2)
    )


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> ConformalDualPolarizationRecord:
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
    joint_dimension = dimension**2

    first_support, _, _ = support_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    shape_map = conformal_shape_map(first_support, dimension)
    joint_shape_map = np.hstack(
        (
            shape_map,
            np.zeros(
                (
                    2 * complex_dimension,
                    2 * complex_dimension,
                )
            ),
        )
    )

    joint_form = endpoint.sharp_form
    inverse_shape_form = (
        joint_shape_map
        @ np.linalg.solve(joint_form, joint_shape_map.T)
    )
    phase_residual = complex_phase_residual(
        inverse_shape_form,
        complex_dimension,
    )

    dual_lifts = np.linalg.solve(joint_form, joint_shape_map.T)
    dual_solve_backward = float(
        np.linalg.norm(
            joint_form @ dual_lifts - joint_shape_map.T,
            2,
        )
        / (
            np.linalg.norm(joint_form, 2)
            * np.linalg.norm(dual_lifts, 2)
            + np.linalg.norm(joint_shape_map, 2)
        )
    )
    analytic_lifts = (
        dual_lifts[:, complex_dimension:]
        + 1j * dual_lifts[:, :complex_dimension]
    )
    hermitian_scale = max(
        np.linalg.norm(
            analytic_lifts.conj().T
            @ joint_form
            @ analytic_lifts,
            2,
        ),
        1e-30,
    )

    def isotropy_residual(form: np.ndarray) -> float:
        return float(
            np.linalg.norm(
                analytic_lifts.T @ form @ analytic_lifts,
                2,
            )
            / hermitian_scale
        )

    sharp_isotropy = isotropy_residual(endpoint.sharp_form)
    right_isotropy = isotropy_residual(endpoint.right_form)
    left_isotropy = isotropy_residual(endpoint.left_form)

    first_operators, _, _ = normalized_operator_jet(
        endpoint.matrix,
        directions,
        angle_count,
    )
    real_images = [
        first_blaschke_image(
            dual_lifts[:, index],
            physical_dimension,
            first_operators,
            endpoint.matrix,
            endpoint.zeros,
        )
        for index in range(2 * complex_dimension)
    ]
    lower_residual = 0.0
    for index in range(complex_dimension):
        analytic_image = (
            real_images[complex_dimension + index]
            + 1j * real_images[index]
        )
        lower_residual = max(
            lower_residual,
            np.linalg.norm(np.tril(analytic_image), 2)
            / max(np.linalg.norm(analytic_image, 2), 1e-30),
        )

    condition_number = float(np.linalg.cond(joint_form))
    checks = (
        joint_form.shape == (joint_dimension, joint_dimension)
        and dual_solve_backward < 2e-14
        and phase_residual < 1e-5
        and sharp_isotropy < 1e-5
        and right_isotropy < 1e-5
        and left_isotropy < 1e-5
        and lower_residual < 1e-6
    )
    if not checks:
        raise RuntimeError(
            "conformal dual-polarization audit failed: "
            f"n={dimension}, sample={sample}, "
            f"phase={phase_residual}, "
            f"sharp={sharp_isotropy}, "
            f"right={right_isotropy}, "
            f"left={left_isotropy}, "
            f"lower={lower_residual}, "
            f"backward={dual_solve_backward}, "
            f"condition={condition_number}"
        )

    return ConformalDualPolarizationRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        joint_dimension=joint_dimension,
        complex_shape_dimension=complex_dimension,
        joint_form_condition_number=condition_number,
        dual_solve_backward_residual=dual_solve_backward,
        phase_covariance_residual=phase_residual,
        sharp_isotropy_residual=sharp_isotropy,
        right_isotropy_residual=right_isotropy,
        left_isotropy_residual=left_isotropy,
        lower_triangular_response_residual=float(lower_residual),
        all_checks_passed=True,
    )


def write_records(
    records: list[ConformalDualPolarizationRecord],
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
        default=(4, 5, 6, 7, 8),
    )
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "gau_wu_conformal_dual_polarization_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic dual-polarization audit."""

    arguments = parse_args()
    records: list[ConformalDualPolarizationRecord] = []
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
                        "phase_residual": (
                            record.phase_covariance_residual
                        ),
                        "sharp_isotropy": (
                            record.sharp_isotropy_residual
                        ),
                        "lower_response": (
                            record.lower_triangular_response_residual
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
