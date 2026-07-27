#!/usr/bin/env python3
"""Audit the conformal-Fourier phase structure of the Gau--Wu shape form."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import support_jet
from gau_wu_similarity_hessian import build_similarity_hessian_audit


@dataclass(frozen=True)
class ConformalShapePhaseRecord:
    """One conformal-coordinate rank and phase-covariance audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    physical_dimension: int
    expected_shape_dimension: int
    observed_shape_rank: int
    shape_kernel_dimension: int
    conformal_coordinate_condition_number: float
    phase_covariance_residual: float
    antiholomorphic_block_relative_norm: float
    minimum_hermitian_eigenvalue: float
    maximum_hermitian_eigenvalue: float
    all_checks_passed: bool


def conformal_shape_map(
    first_support: np.ndarray,
    dimension: int,
) -> np.ndarray:
    """Map physical directions to real coordinates of support modes 2 through n."""

    angle_count = first_support.shape[1]
    fourier = np.fft.fft(first_support, axis=1) / angle_count
    # The factor two is the Schwarz-transform convention.  It has no
    # effect on rank or phase covariance, but makes these coordinates
    # the actual coefficients of w H_s(w).
    analytic = 2 * fourier[:, 2 : dimension + 1].T
    return np.vstack((analytic.real, analytic.imag))


def complex_blocks(
    real_form: np.ndarray,
    complex_dimension: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return Hermitian and complex-symmetric parts of a real quadratic form."""

    xx = real_form[:complex_dimension, :complex_dimension]
    xy = real_form[:complex_dimension, complex_dimension:]
    yy = real_form[complex_dimension:, complex_dimension:]
    hermitian = (xx + yy) / 2 + 0.5j * (xy.T - xy)
    symmetric = (xx - yy) / 2 - 0.5j * (xy + xy.T)
    return hermitian, symmetric


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> ConformalShapePhaseRecord:
    """Audit one finite nondegenerate Gau--Wu equality model."""

    audit = build_similarity_hessian_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    endpoint = audit.endpoint_audit
    first_support, _, _ = support_jet(
        endpoint.matrix,
        list(endpoint.directions),
        angle_count,
    )
    shape_map = conformal_shape_map(first_support, dimension)
    singular_values = np.linalg.svd(shape_map, compute_uv=False)
    threshold = 1e-8 * singular_values[0]
    observed_rank = int(np.sum(singular_values > threshold))
    expected_shape_dimension = 2 * dimension - 2

    inverse_form = np.linalg.solve(
        audit.optimized_similarity_form,
        shape_map.T,
    )
    quotient = np.linalg.inv(shape_map @ inverse_form)
    quotient = (quotient + quotient.T) / 2

    complex_dimension = dimension - 1
    zero = np.zeros((complex_dimension, complex_dimension))
    identity = np.eye(complex_dimension)
    phase = np.block([[zero, -identity], [identity, zero]])
    phase_residual = np.linalg.norm(
        quotient - phase.T @ quotient @ phase,
        2,
    ) / np.linalg.norm(quotient, 2)

    hermitian, symmetric = complex_blocks(quotient, complex_dimension)
    antiholomorphic_relative = (
        np.linalg.norm(symmetric, 2) / np.linalg.norm(hermitian, 2)
    )
    eigenvalues = np.linalg.eigvalsh(hermitian)

    physical_dimension = len(endpoint.directions)
    shape_kernel_dimension = physical_dimension - observed_rank
    checks = (
        observed_rank == expected_shape_dimension
        and shape_kernel_dimension == (dimension - 2) ** 2
        and phase_residual < 2e-8
        and antiholomorphic_relative < 2e-8
        and eigenvalues[-1] < 2e-8
    )
    if not checks:
        raise RuntimeError(
            "conformal-shape phase audit failed: "
            f"n={dimension}, sample={sample}, "
            f"rank={observed_rank}/{expected_shape_dimension}, "
            f"kernel={shape_kernel_dimension}/{(dimension - 2) ** 2}, "
            f"phase={phase_residual}, "
            f"antiholomorphic={antiholomorphic_relative}, "
            f"lambda_max={eigenvalues[-1]}"
        )

    return ConformalShapePhaseRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        physical_dimension=physical_dimension,
        expected_shape_dimension=expected_shape_dimension,
        observed_shape_rank=observed_rank,
        shape_kernel_dimension=shape_kernel_dimension,
        conformal_coordinate_condition_number=float(
            singular_values[0] / singular_values[-1]
        ),
        phase_covariance_residual=float(phase_residual),
        antiholomorphic_block_relative_norm=float(
            antiholomorphic_relative
        ),
        minimum_hermitian_eigenvalue=float(eigenvalues[0]),
        maximum_hermitian_eigenvalue=float(eigenvalues[-1]),
        all_checks_passed=True,
    )


def write_records(
    records: list[ConformalShapePhaseRecord],
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
            "experiments/gau_wu_conformal_shape_phase_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic conformal-shape phase audit."""

    arguments = parse_args()
    records: list[ConformalShapePhaseRecord] = []
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
                        "shape_rank": record.observed_shape_rank,
                        "phase_residual": (
                            record.phase_covariance_residual
                        ),
                        "maximum_hermitian_eigenvalue": (
                            record.maximum_hermitian_eigenvalue
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
