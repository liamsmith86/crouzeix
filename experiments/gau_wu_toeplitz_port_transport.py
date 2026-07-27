#!/usr/bin/env python3
"""Transport the Gau--Wu Ando port to the fixed Toeplitz coefficient shift."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_disk_chart_recenter import (
    chart_data,
    polynomial_support_frame,
)
from gau_wu_finite_hessian_jet import support_jet
from gau_wu_finite_model import (
    extremal_zeros,
    gau_wu_model,
    normal_basis,
    random_interior_zeros,
    real_matrix,
)
from gau_wu_second_support_gram import boundary_angular_derivative
from gau_wu_support_port_energy import build_support_port_data


@dataclass(frozen=True)
class ToeplitzPortTransportRecord:
    """One fixed-shift transport audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    normal_dimension: int
    maximum_spectral_factor_residual: float
    maximum_fixed_range_residual: float
    maximum_pointwise_energy_relative_error: float
    maximum_transport_gram_residual: float
    maximum_support_gram_residual: float
    all_checks_passed: bool


def positive_semidefinite_sqrt(matrix: np.ndarray) -> np.ndarray:
    """Return the positive square root, clipping only roundoff at zero."""

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    if eigenvalues[0] < -1e-10 * eigenvalues[-1]:
        raise ValueError("matrix is not positive semidefinite")
    return (
        eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))
    ) @ eigenvectors.conj().T


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> ToeplitzPortTransportRecord:
    """Audit the fixed Toeplitz transport at one Gau--Wu model."""

    rng = np.random.default_rng(seed)
    interior = random_interior_zeros(dimension, rng)
    matrix = gau_wu_model(interior)
    zeros = extremal_zeros(interior)
    normal, _, _ = normal_basis(interior)
    directions = [
        real_matrix(normal[:, index], dimension)
        for index in range(normal.shape[1])
    ]
    _, second_support, _ = support_jet(
        matrix,
        directions,
        angle_count,
    )

    frame = polynomial_support_frame(matrix, zeros)
    _, extended_chart, _ = chart_data(frame, matrix)
    chart_sqrt = positive_semidefinite_sqrt(extended_chart)
    chart_sqrt_pseudoinverse = np.linalg.pinv(
        chart_sqrt,
        rcond=1e-12,
    )
    coefficient_shift = np.diag(np.ones(dimension - 1), 1)
    port = build_support_port_data(matrix)
    identity = np.eye(dimension)

    angles = np.linspace(0, 2 * np.pi, angle_count, endpoint=False)
    angular_derivative = boundary_angular_derivative(zeros, angles)
    expected_gram = 2 * np.mean(
        angular_derivative[None, None, :] * second_support,
        axis=2,
    )
    original_gram = np.zeros_like(expected_gram)
    transported_gram = np.zeros_like(expected_gram)
    maximum_factor_residual = 0.0
    maximum_range_residual = 0.0
    maximum_relative_energy_error = 0.0

    for angle in angles:
        phase = np.exp(1j * angle)
        support_factor = (
            port.unitary @ port.d_matrix - phase * port.e_matrix
        )
        coefficient_factor = chart_sqrt @ (
            identity - np.conj(phase) * coefficient_shift
        )
        maximum_factor_residual = max(
            maximum_factor_residual,
            np.linalg.norm(
                (support_factor @ frame).conj().T
                @ (support_factor @ frame)
                - coefficient_factor.conj().T @ coefficient_factor,
                2,
            ),
        )

        model_kernel = np.linalg.solve(
            identity - np.conj(phase) * port.shift,
            port.terminal,
        )
        null_state = port.null_scaling @ model_kernel
        inverse_support_adjoint = np.linalg.pinv(
            support_factor.conj().T,
            rcond=1e-12,
        )
        original_responses = np.empty(
            (dimension, len(directions)),
            dtype=complex,
        )
        transported_responses = np.empty_like(original_responses)
        for index, direction in enumerate(directions):
            support_direction = (
                np.conj(phase) * direction
                + phase * direction.conj().T
            ) / 2
            image = support_direction @ null_state
            support_motion = (
                np.vdot(null_state, image)
                / np.vdot(null_state, null_state)
            ).real
            orthogonal_image = image - support_motion * null_state
            original = inverse_support_adjoint @ orthogonal_image

            fixed_range_vector = np.linalg.solve(
                identity - phase * coefficient_shift.conj().T,
                frame.conj().T @ orthogonal_image,
            )
            maximum_range_residual = max(
                maximum_range_residual,
                abs(fixed_range_vector[-1]),
            )
            transported = (
                chart_sqrt_pseudoinverse @ fixed_range_vector
            )
            original_responses[:, index] = original
            transported_responses[:, index] = transported
            energy_scale = max(1.0, np.linalg.norm(original) ** 2)
            maximum_relative_energy_error = max(
                maximum_relative_energy_error,
                abs(
                    np.linalg.norm(original) ** 2
                    - np.linalg.norm(transported) ** 2
                )
                / energy_scale,
            )

        original_gram += (
            original_responses.conj().T @ original_responses
        ).real / angle_count
        transported_gram += (
            transported_responses.conj().T @ transported_responses
        ).real / angle_count

    transport_gram_residual = float(
        np.max(abs(original_gram - transported_gram))
    )
    support_gram_residual = float(
        np.max(abs(transported_gram - expected_gram))
    )
    checks = (
        maximum_factor_residual < 2e-10
        and maximum_range_residual < 2e-10
        and maximum_relative_energy_error < 2e-10
        and transport_gram_residual < 5e-9
        and support_gram_residual < 5e-9
    )
    if not checks:
        raise RuntimeError(
            "Toeplitz port transport audit failed: "
            f"n={dimension}, factor={maximum_factor_residual}, "
            f"range={maximum_range_residual}, "
            f"energy={maximum_relative_energy_error}, "
            f"transport={transport_gram_residual}, "
            f"support={support_gram_residual}"
        )
    return ToeplitzPortTransportRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        normal_dimension=len(directions),
        maximum_spectral_factor_residual=float(maximum_factor_residual),
        maximum_fixed_range_residual=float(maximum_range_residual),
        maximum_pointwise_energy_relative_error=float(
            maximum_relative_energy_error
        ),
        maximum_transport_gram_residual=transport_gram_residual,
        maximum_support_gram_residual=support_gram_residual,
        all_checks_passed=True,
    )


def write_records(
    records: list[ToeplitzPortTransportRecord],
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
    parser.add_argument("--angle-count", type=int, default=256)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_toeplitz_port_transport_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic fixed-shift transport audit."""

    arguments = parse_args()
    records: list[ToeplitzPortTransportRecord] = []
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
                        "factor": (
                            record.maximum_spectral_factor_residual
                        ),
                        "transport": (
                            record.maximum_transport_gram_residual
                        ),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
    digest = write_records(records, arguments.output)
    print(f"wrote {len(records)} records to {arguments.output}")
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
