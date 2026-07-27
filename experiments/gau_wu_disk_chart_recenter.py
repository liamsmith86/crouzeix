#!/usr/bin/env python3
"""Recenter the inverse-Toeplitz disk chart at Gau--Wu equality models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import support_jet
from gau_wu_finite_model import real_vectorize
from gau_wu_second_support_gram import boundary_angular_derivative
from gau_wu_similarity_hessian import build_similarity_hessian_audit


@dataclass(frozen=True)
class DiskChartRecenterRecord:
    """One polynomial-frame and disk-tangent quotient audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    support_frame_condition_number: float
    polynomial_null_error: float
    chart_identity_error: float
    inverse_toeplitz_error: float
    minimum_chart_eigenvalue: float
    disk_tangent_dimension: int
    observed_disk_tangent_rank: int
    maximum_disk_hessian_eigenvalue: float
    shape_quotient_dimension: int
    observed_shape_quotient_rank: int
    maximum_shape_polynomial_tail: float
    disk_shape_kernel_projector_error: float
    maximum_shape_schur_eigenvalue: float
    all_checks_passed: bool


def hermitian_sqrt(matrix: np.ndarray) -> np.ndarray:
    """Return the positive square root of a positive Hermitian matrix."""

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return (
        eigenvectors * np.sqrt(eigenvalues)
    ) @ eigenvectors.conj().T


def hermitian_basis(dimension: int) -> tuple[np.ndarray, ...]:
    """Return an orthonormal real basis of Hermitian matrices."""

    basis: list[np.ndarray] = []
    for index in range(dimension):
        diagonal = np.zeros((dimension, dimension), dtype=complex)
        diagonal[index, index] = 1
        basis.append(diagonal)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            real = np.zeros((dimension, dimension), dtype=complex)
            real[row, column] = 1 / np.sqrt(2)
            real[column, row] = 1 / np.sqrt(2)
            basis.append(real)

            imaginary = np.zeros((dimension, dimension), dtype=complex)
            imaginary[row, column] = 1j / np.sqrt(2)
            imaginary[column, row] = -1j / np.sqrt(2)
            basis.append(imaginary)
    return tuple(basis)


def polynomial_support_frame(
    matrix: np.ndarray,
    extremal_zeros: np.ndarray,
) -> np.ndarray:
    """Build the polynomial support frame from the model resolvent."""

    dimension = matrix.shape[0]
    scaling = np.ones(dimension)
    scaling[0] = np.sqrt(2)
    scaling[-1] = 1 / np.sqrt(2)
    similarity = np.diag(scaling)
    shift = np.linalg.solve(similarity, matrix @ similarity)

    state_scaling = np.ones(dimension)
    state_scaling[1:-1] = np.sqrt(2)
    state_scaling = np.diag(state_scaling)
    terminal = np.eye(dimension, dtype=complex)[:, -1]

    phi_zeros = np.asarray((0, *extremal_zeros), dtype=complex)
    determinant = np.asarray((1.0 + 0j,))
    for zero in phi_zeros:
        determinant = np.convolve(determinant, (1, -zero))

    adjugate_columns: list[np.ndarray] = []
    for degree in range(dimension):
        column = sum(
            (
                determinant[index]
                * np.linalg.matrix_power(shift, degree - index)
                @ terminal
                for index in range(degree + 1)
            ),
            np.zeros(dimension, dtype=complex),
        )
        adjugate_columns.append(state_scaling @ column)
    return np.stack(adjugate_columns[::-1], axis=1)


def chart_data(
    frame: np.ndarray,
    matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Recover ``K``, the extended chart matrix ``H``, and coefficient ``A``."""

    dimension = matrix.shape[0]
    metric = frame.conj().T @ frame
    chart = np.zeros_like(metric)
    for row in range(dimension):
        for column in range(dimension):
            predecessor = (
                chart[row - 1, column - 1]
                if row and column
                else 0
            )
            chart[row, column] = metric[row, column] - predecessor
    coefficient_matrix = np.linalg.solve(frame, matrix @ frame)
    return metric, chart, coefficient_matrix


def disk_chart(matrix: np.ndarray) -> np.ndarray:
    """Evaluate the Euclidean disk chart at a positive Hermitian matrix."""

    length = matrix.shape[0]
    dimension = length + 1
    extended = np.zeros((dimension, dimension), dtype=complex)
    extended[:length, :length] = matrix
    shift = np.diag(np.ones(dimension - 1), 1)
    metric = extended + shift.conj().T @ extended @ shift
    square_root = hermitian_sqrt(metric)
    coefficient = 2 * np.linalg.solve(metric, extended @ shift)
    return square_root @ coefficient @ np.linalg.inv(square_root)


def trig_coordinates(values: np.ndarray, degree: int) -> np.ndarray:
    """Return orthonormal real trigonometric coordinates through ``degree``."""

    count = len(values)
    angles = 2 * np.pi * np.arange(count) / count
    rows = [np.ones(count) / np.sqrt(count)]
    for frequency in range(1, degree + 1):
        rows.extend(
            (
                np.sqrt(2 / count) * np.cos(frequency * angles),
                np.sqrt(2 / count) * np.sin(frequency * angles),
            )
        )
    return np.asarray(rows) @ values


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> DiskChartRecenterRecord:
    """Audit one arbitrary finite Gau--Wu equality model."""

    audit = build_similarity_hessian_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    endpoint = audit.endpoint_audit
    matrix = endpoint.matrix
    frame = polynomial_support_frame(matrix, endpoint.zeros)
    metric, extended_chart, coefficient = chart_data(frame, matrix)
    chart = extended_chart[:-1, :-1]
    inverse_chart = np.linalg.inv(chart)
    shift = np.diag(np.ones(dimension - 1), 1)

    angles = np.linspace(0, 2 * np.pi, 17, endpoint=False)
    powers = np.arange(dimension)
    maximum_null_error = 0.0
    for angle in angles:
        point = np.exp(1j * angle)
        vector = frame @ point**powers
        support_slack = np.eye(dimension) - (
            np.conj(point) * matrix
            + point * matrix.conj().T
        ) / 2
        maximum_null_error = max(
            maximum_null_error,
            np.linalg.norm(support_slack @ vector),
        )

    chart_identity_error = max(
        np.linalg.norm(
            metric
            - extended_chart
            - shift.conj().T @ extended_chart @ shift,
            2,
        ),
        np.linalg.norm(
            metric @ coefficient - 2 * extended_chart @ shift,
            2,
        ),
        np.linalg.norm(extended_chart[-1, :]),
        np.linalg.norm(extended_chart[:, -1]),
    )
    inverse_toeplitz_error = max(
        abs(
            inverse_chart[row, column]
            - inverse_chart[row - 1, column - 1]
        )
        for row in range(1, dimension - 1)
        for column in range(1, dimension - 1)
    )

    square_root = hermitian_sqrt(metric)
    unitary = frame @ np.linalg.inv(square_root)
    chart_base = disk_chart(chart)
    base_error = np.linalg.norm(
        matrix - unitary @ chart_base @ unitary.conj().T,
        2,
    )
    chart_identity_error = max(chart_identity_error, base_error)

    normal = np.stack(
        [real_vectorize(direction) for direction in endpoint.directions],
        axis=1,
    )
    step = 2e-6
    tangent_columns: list[np.ndarray] = []
    for direction in hermitian_basis(dimension - 1):
        derivative = (
            disk_chart(chart + step * direction)
            - disk_chart(chart - step * direction)
        ) / (2 * step)
        physical_derivative = unitary @ derivative @ unitary.conj().T
        tangent_columns.append(
            normal.T @ real_vectorize(physical_derivative)
        )
    tangent_map = np.asarray(tangent_columns).T
    tangent_left, tangent_singular, _ = np.linalg.svd(
        tangent_map,
        full_matrices=False,
    )
    expected_tangent_dimension = (dimension - 2) ** 2
    tangent_threshold = 1e-7 * tangent_singular[0]
    observed_tangent_rank = int(
        np.sum(tangent_singular > tangent_threshold)
    )
    tangent_basis = tangent_left[:, :observed_tangent_rank]

    disk_hessian = (
        tangent_basis.T
        @ audit.optimized_similarity_form
        @ tangent_basis
    )
    maximum_disk_hessian = float(
        np.linalg.eigvalsh((disk_hessian + disk_hessian.T) / 2)[-1]
    )

    first_support, _, _ = support_jet(
        matrix,
        list(endpoint.directions),
        angle_count,
    )
    boundary_angles = np.linspace(
        0,
        2 * np.pi,
        angle_count,
        endpoint=False,
    )
    angular_derivative = boundary_angular_derivative(
        endpoint.zeros,
        boundary_angles,
    )
    boundary_points = np.exp(1j * boundary_angles)
    phi_zeros = np.asarray((0, *endpoint.zeros), dtype=complex)
    determinant_values = np.prod(
        1
        - phi_zeros[:, None]
        * np.conj(boundary_points)[None, :],
        axis=0,
    )
    polynomial_weight = (
        abs(determinant_values) ** 2 * angular_derivative
    )
    weighted_support = polynomial_weight[None, :] * first_support
    weighted_fourier = np.fft.fft(weighted_support, axis=1) / angle_count
    polynomial_tail = weighted_fourier[
        :,
        dimension + 1 : angle_count - dimension,
    ]
    maximum_polynomial_tail = float(np.max(abs(polynomial_tail)))
    support_map = np.stack(
        [
            trig_coordinates(
                weighted_support[index],
                dimension,
            )
            for index in range(len(endpoint.directions))
        ],
        axis=1,
    )
    affine_map = np.stack(
        [
            trig_coordinates(polynomial_weight * function, dimension)
            for function in (
                np.ones(angle_count),
                np.cos(boundary_angles),
                np.sin(boundary_angles),
            )
        ],
        axis=1,
    )
    affine_left, affine_singular, _ = np.linalg.svd(
        affine_map,
        full_matrices=True,
    )
    affine_rank = int(
        np.sum(affine_singular > 1e-9 * affine_singular[0])
    )
    shape_map = affine_left[:, affine_rank:].T @ support_map
    _, shape_singular, shape_right = np.linalg.svd(
        shape_map,
        full_matrices=True,
    )
    expected_shape_dimension = 2 * dimension - 2
    shape_threshold = 1e-7 * shape_singular[0]
    observed_shape_rank = int(
        np.sum(shape_singular > shape_threshold)
    )
    shape_kernel = shape_right[observed_shape_rank:].T
    kernel_projector_error = np.linalg.norm(
        tangent_basis @ tangent_basis.T
        - shape_kernel @ shape_kernel.T,
        2,
    )

    form = audit.optimized_similarity_form
    quotient_form = np.linalg.inv(
        shape_map @ np.linalg.solve(form, shape_map.T)
    )
    quotient_form = (quotient_form + quotient_form.T) / 2
    maximum_shape_schur = float(np.linalg.eigvalsh(quotient_form)[-1])

    checks = (
        maximum_null_error < 2e-11
        and chart_identity_error < 5e-10
        and inverse_toeplitz_error < 5e-9
        and np.linalg.eigvalsh(chart)[0] > 1e-6
        and observed_tangent_rank == expected_tangent_dimension
        and maximum_disk_hessian < -1e-4
        and observed_shape_rank == expected_shape_dimension
        and maximum_polynomial_tail < 2e-10
        and kernel_projector_error < 2e-6
        and maximum_shape_schur < 2e-8
    )
    if not checks:
        raise RuntimeError(
            "disk-chart recenter audit failed: "
            f"n={dimension}, null={maximum_null_error}, "
            f"chart={chart_identity_error}, "
            f"Toeplitz={inverse_toeplitz_error}, "
            f"tangent={observed_tangent_rank}/"
            f"{expected_tangent_dimension}, "
            f"disk={maximum_disk_hessian}, "
            f"shape={observed_shape_rank}/{expected_shape_dimension}, "
            f"tail={maximum_polynomial_tail}, "
            f"kernel={kernel_projector_error}, "
            f"Schur={maximum_shape_schur}"
        )
    return DiskChartRecenterRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        support_frame_condition_number=float(np.linalg.cond(frame)),
        polynomial_null_error=float(maximum_null_error),
        chart_identity_error=float(chart_identity_error),
        inverse_toeplitz_error=float(inverse_toeplitz_error),
        minimum_chart_eigenvalue=float(np.linalg.eigvalsh(chart)[0]),
        disk_tangent_dimension=expected_tangent_dimension,
        observed_disk_tangent_rank=observed_tangent_rank,
        maximum_disk_hessian_eigenvalue=maximum_disk_hessian,
        shape_quotient_dimension=expected_shape_dimension,
        observed_shape_quotient_rank=observed_shape_rank,
        maximum_shape_polynomial_tail=maximum_polynomial_tail,
        disk_shape_kernel_projector_error=float(kernel_projector_error),
        maximum_shape_schur_eigenvalue=maximum_shape_schur,
        all_checks_passed=True,
    )


def write_records(
    records: list[DiskChartRecenterRecord],
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
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_disk_chart_recenter_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic polynomial-frame recentering audit."""

    arguments = parse_args()
    records: list[DiskChartRecenterRecord] = []
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
                        "disk_rank": record.observed_disk_tangent_rank,
                        "shape_rank": record.observed_shape_quotient_rank,
                        "disk_max": record.maximum_disk_hessian_eigenvalue,
                        "shape_max": record.maximum_shape_schur_eigenvalue,
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
