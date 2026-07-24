#!/usr/bin/env python3
"""Audit the full general-H disk chart on the Crabb weighted face.

L174 controls the Toeplitz disk slice plus a fixed coercive circular-normal
fibre.  The omitted general-H disk coordinates occur at the same weighted
order as those normal variables.  This script asks the decisive leading-order
question: after allowing *all* L65 strong directions, does the quadratic
completion exceed the Toeplitz quartic deficit, and is its maximizing class
realized by a tangent of L122's exact general-H disk chart?

The calculation is deliberately independent of L173's closed Pluecker
formula.  It reconstructs L65's full real Hessian, extracts the quadratic
ambient gradient of the characteristic-Blaschke lower certificate by
symmetric Richardson extrapolation, and compares the unconstrained maximizer
with the tangent image of arbitrary Hermitian H-variations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crabb_off_equality_dual_gradient import (
    dual_gradient_and_coercive_rows,
    quartic_defect,
)
from general_crabb_second_order_modes import full_quadratic_form


@dataclass(frozen=True)
class FullDiskFaceRecord:
    """One dimension's full weighted-face audit."""

    dimension: int
    length: int
    direction_index: int
    seed: int
    hessian_rank: int
    expected_hessian_rank: int
    disk_tangent_rank: int
    expected_disk_tangent_rank: int
    disk_curvature_rank: int
    expected_disk_curvature_rank: int
    quartic_defect: float
    disk_deficit_coefficient: float
    full_strong_gain: float
    gain_over_disk_deficit: float
    gradient_kernel_residual: float
    disk_tangent_stationarity_residual: float
    explicit_correction_energy_ratio: float
    explicit_curvature_energy_ratio: float
    maximum_hessian_eigenvalue: float
    minimum_positive_curvature: float


@dataclass(frozen=True)
class FullDiskGeometry:
    """Dimension-dependent objects shared by all full-disk face tests."""

    dimension: int
    physical_hessian: np.ndarray
    curvature: np.ndarray
    gauge_map: np.ndarray
    tangent: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    positive: np.ndarray


def deterministic_direction(length: int) -> np.ndarray:
    """Return a normalized generic non-palindromic Toeplitz direction."""

    direction = np.zeros(length, dtype=complex)
    for offset in range(1, length):
        direction[offset] = (
            (2 * offset + 1) / (7 * length)
            + 1j * (offset + 2) / (11 * length)
        )
    norm = np.linalg.norm(direction[1:])
    if norm == 0:
        raise RuntimeError("the deterministic disk direction vanished")
    return direction / norm


def audit_directions(
    length: int,
    count: int,
    seed: int,
) -> list[np.ndarray]:
    """Return one structured direction followed by seeded generic ones."""

    directions = [deterministic_direction(length)]
    rng = np.random.default_rng(seed + length)
    for _ in range(1, count):
        direction = np.zeros(length, dtype=complex)
        direction[1:] = (
            rng.standard_normal(length - 1)
            + 1j * rng.standard_normal(length - 1)
        )
        direction /= np.linalg.norm(direction[1:])
        directions.append(direction)
    return directions


def real_coordinates(matrix: np.ndarray) -> np.ndarray:
    """Use the repository's standard real-then-imaginary matrix ordering."""

    return np.concatenate((matrix.real.ravel(), matrix.imag.ravel()))


def coordinate_gauge_map(dimension: int) -> np.ndarray:
    """Map coefficient-gauge coordinates to physical Crabb coordinates."""

    metric_diagonal = np.ones(dimension)
    metric_diagonal[[0, -1]] = 0.5
    square_root = np.sqrt(metric_diagonal)
    inverse_square_root = 1 / square_root
    coordinate_count = 2 * dimension**2
    transform = np.zeros((coordinate_count, coordinate_count))
    column = 0
    for scalar in (1, 1j):
        for index in range(dimension**2):
            coefficient = np.zeros((dimension, dimension), dtype=complex)
            coefficient.flat[index] = scalar
            physical = (
                square_root[:, None]
                * coefficient
                * inverse_square_root[None, :]
            )
            transform[:, column] = real_coordinates(physical)
            column += 1
    return transform


def hermitian_basis(length: int) -> list[np.ndarray]:
    """Return a real basis of Hermitian ``length x length`` matrices."""

    basis: list[np.ndarray] = []
    for row in range(length):
        diagonal = np.zeros((length, length), dtype=complex)
        diagonal[row, row] = 1
        basis.append(diagonal)
    for row in range(length):
        for column in range(row + 1, length):
            real = np.zeros((length, length), dtype=complex)
            real[row, column] = real[column, row] = 1
            basis.append(real)
            imaginary = np.zeros((length, length), dtype=complex)
            imaginary[row, column] = 1j
            imaginary[column, row] = -1j
            basis.append(imaginary)
    return basis


def hermitian_coefficients(matrix: np.ndarray) -> np.ndarray:
    """Encode a Hermitian matrix in ``hermitian_basis`` ordering."""

    length = len(matrix)
    coefficients = [float(np.real(matrix[index, index])) for index in range(length)]
    for row in range(length):
        for column in range(row + 1, length):
            coefficients.extend(
                (
                    float(np.real(matrix[row, column])),
                    float(np.imag(matrix[row, column])),
                )
            )
    return np.asarray(coefficients)


def plucker_disk_correction(direction: np.ndarray) -> np.ndarray:
    """Return L176's explicit Hermitian disk correction."""

    coefficients = direction[1:]
    coefficient_count = len(coefficients)
    length = len(direction)
    reversal = np.flip(np.conj(coefficients))
    plucker = (
        np.outer(coefficients, reversal)
        - np.outer(reversal, coefficients)
    )
    correction = np.zeros((length, length), dtype=complex)
    for anti_diagonal in range(1, coefficient_count):
        offset = coefficient_count - 1 - anti_diagonal
        for row in range(anti_diagonal + 2):
            value = 0j
            for left in range((anti_diagonal + 1) // 2):
                right = anti_diagonal - left
                if left >= right:
                    continue
                flux = 2 * (
                    left + 1 <= row <= anti_diagonal - left
                )
                mean = 2 * (anti_diagonal - 2 * left) / (
                    anti_diagonal + 2
                )
                value += (flux - mean) * np.conj(
                    plucker[left, right]
                )
            if offset == 0:
                correction[row, row] = np.real(value)
            else:
                correction[row, row + offset] = value
                correction[row + offset, row] = np.conj(value)
    return correction


def inverse_square_root_frechet(
    diagonal: np.ndarray,
    direction: np.ndarray,
) -> np.ndarray:
    """Differentiate ``K -> K^(-1/2)`` at a positive diagonal matrix."""

    values = diagonal ** (-0.5)
    derivative = np.empty_like(direction)
    for row, left in enumerate(diagonal):
        for column, right in enumerate(diagonal):
            if left == right:
                divided_difference = -0.5 * left ** (-1.5)
            else:
                divided_difference = (
                    values[row] - values[column]
                ) / (left - right)
            derivative[row, column] = (
                divided_difference * direction[row, column]
            )
    return derivative


def disk_tangent_map(length: int) -> np.ndarray:
    """Differentiate L122's physical-gauge general-H disk chart."""

    dimension = length + 1
    shift = np.zeros((dimension, dimension), dtype=complex)
    for row in range(length):
        shift[row, row + 1] = 1
    h_zero = np.zeros((dimension, dimension), dtype=complex)
    h_zero[:length, :length] = np.eye(length) / 2
    metric = h_zero + shift.conj().T @ h_zero @ shift
    metric_diagonal = np.real(np.diag(metric))
    inverse_square_root = np.diag(metric_diagonal ** (-0.5))

    columns = []
    for small_direction in hermitian_basis(length):
        direction = np.zeros_like(h_zero)
        direction[:length, :length] = small_direction
        metric_direction = (
            direction + shift.conj().T @ direction @ shift
        )
        operator_direction = (
            2
            * inverse_square_root_frechet(
                metric_diagonal,
                metric_direction,
            )
            @ h_zero
            @ shift
            @ inverse_square_root
            + 2
            * inverse_square_root
            @ direction
            @ shift
            @ inverse_square_root
            + 2
            * inverse_square_root
            @ h_zero
            @ shift
            @ inverse_square_root_frechet(
                metric_diagonal,
                metric_direction,
            )
        )
        columns.append(real_coordinates(operator_direction))
    return np.column_stack(columns)


def quadratic_gradient(
    direction: np.ndarray,
    scale: float,
    resolution: int,
) -> np.ndarray:
    """Extract the even quadratic ambient-gradient coefficient."""

    def even_coefficient(current_scale: float) -> np.ndarray:
        _, positive, _ = dual_gradient_and_coercive_rows(
            current_scale * direction,
            resolution,
        )
        _, negative, _ = dual_gradient_and_coercive_rows(
            -current_scale * direction,
            resolution,
        )
        return (positive + negative) / (2 * current_scale**2)

    coarse = even_coefficient(scale)
    fine = even_coefficient(scale / 2)
    return (4 * fine - coarse) / 3


def full_disk_geometry(length: int) -> FullDiskGeometry:
    """Build the L65 curvature and general-H disk tangent once per size."""

    dimension = length + 1
    physical_hessian = full_quadratic_form(
        dimension,
        max(64, 4 * dimension),
    )
    curvature = -(physical_hessian + physical_hessian.T) / 2
    gauge_map = coordinate_gauge_map(dimension)
    tangent = disk_tangent_map(length)
    eigenvalues, eigenvectors = np.linalg.eigh(curvature)
    tolerance = 1e-8 * max(1.0, float(eigenvalues[-1]))
    positive = eigenvalues > tolerance
    return FullDiskGeometry(
        dimension=dimension,
        physical_hessian=physical_hessian,
        curvature=curvature,
        gauge_map=gauge_map,
        tangent=tangent,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        positive=positive,
    )


def disk_recentering_correction(
    geometry: FullDiskGeometry,
    direction: np.ndarray,
    gradient_scale: float,
    resolution: int,
) -> tuple[np.ndarray, float]:
    """Return the Hermitian second-order disk correction and residual."""

    coefficient_gradient = quadratic_gradient(
        direction,
        gradient_scale,
        resolution,
    )
    physical_gradient = np.linalg.solve(
        geometry.gauge_map.T,
        coefficient_gradient,
    )
    correction = plucker_disk_correction(direction)
    coefficients = hermitian_coefficients(correction)
    stationarity_residual = float(
        np.linalg.norm(
            2
            * geometry.curvature
            @ geometry.tangent
            @ coefficients
            - physical_gradient
        )
        / max(np.finfo(float).eps, np.linalg.norm(physical_gradient))
    )
    return correction, stationarity_residual


def make_record(
    geometry: FullDiskGeometry,
    direction: np.ndarray,
    direction_index: int,
    seed: int,
    scale: float,
    resolution: int,
) -> FullDiskFaceRecord:
    """Construct one full-Hessian/full-disk comparison."""

    dimension = geometry.dimension
    length = dimension - 1
    eigenvalues = geometry.eigenvalues
    eigenvectors = geometry.eigenvectors
    positive = geometry.positive
    rank = int(np.sum(positive))
    expected_rank = dimension * (dimension - 2)

    coefficient_gradient = quadratic_gradient(
        direction,
        scale,
        resolution,
    )
    gradient = np.linalg.solve(
        geometry.gauge_map.T,
        coefficient_gradient,
    )
    range_projection = (
        eigenvectors[:, positive]
        @ (eigenvectors[:, positive].T @ gradient)
    )
    kernel_residual = float(
        np.linalg.norm(gradient - range_projection)
        / max(np.finfo(float).eps, np.linalg.norm(gradient))
    )
    inverse_gradient = (
        eigenvectors[:, positive]
        @ (
            (eigenvectors[:, positive].T @ gradient)
            / eigenvalues[positive]
        )
    )
    gain = float(0.25 * gradient @ inverse_gradient)

    tangent = geometry.tangent
    tangent_rank = int(np.linalg.matrix_rank(tangent, tol=1e-9))
    expected_tangent_rank = length**2 - 1
    disk_curvature_rank = int(
        np.linalg.matrix_rank(
            geometry.curvature @ tangent,
            tol=1e-8,
        )
    )
    expected_disk_curvature_rank = (dimension - 2) ** 2
    explicit_coefficients = hermitian_coefficients(
        plucker_disk_correction(direction)
    )
    stationarity_residual = float(
        np.linalg.norm(
            2
            * geometry.curvature
            @ tangent
            @ explicit_coefficients
            - gradient
        )
        / max(np.finfo(float).eps, np.linalg.norm(gradient))
    )

    defect = quartic_defect(direction)
    disk_deficit = 32 * defect
    explicit_image = tangent @ explicit_coefficients
    explicit_energy = 0.5 * float(gradient @ explicit_image)
    explicit_curvature = float(
        explicit_image @ geometry.curvature @ explicit_image
    )
    return FullDiskFaceRecord(
        dimension=dimension,
        length=length,
        direction_index=direction_index,
        seed=seed,
        hessian_rank=rank,
        expected_hessian_rank=expected_rank,
        disk_tangent_rank=tangent_rank,
        expected_disk_tangent_rank=expected_tangent_rank,
        disk_curvature_rank=disk_curvature_rank,
        expected_disk_curvature_rank=expected_disk_curvature_rank,
        quartic_defect=defect,
        disk_deficit_coefficient=disk_deficit,
        full_strong_gain=gain,
        gain_over_disk_deficit=gain / disk_deficit,
        gradient_kernel_residual=kernel_residual,
        disk_tangent_stationarity_residual=stationarity_residual,
        explicit_correction_energy_ratio=(
            explicit_energy / disk_deficit
        ),
        explicit_curvature_energy_ratio=(
            explicit_curvature / disk_deficit
        ),
        maximum_hessian_eigenvalue=float(
            np.linalg.eigvalsh(geometry.physical_hessian)[-1]
        ),
        minimum_positive_curvature=float(eigenvalues[positive][0]),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=5)
    parser.add_argument("--scale", type=float, default=0.002)
    parser.add_argument("--resolution", type=int, default=4096)
    parser.add_argument("--direction-count", type=int, default=4)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the full-disk weighted-face audit."""

    args = parse_args()
    if (
        args.minimum_length < 2
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("require 2 <= minimum length <= maximum length")
    if args.scale <= 0:
        raise ValueError("scale must be positive")
    if args.direction_count < 1:
        raise ValueError("direction count must be positive")
    if args.resolution < 512 or args.resolution & (args.resolution - 1):
        raise ValueError("resolution must be a power of two at least 512")

    records = []
    for length in range(
        args.minimum_length,
        args.maximum_length + 1,
    ):
        geometry = full_disk_geometry(length)
        directions = audit_directions(
            length,
            args.direction_count,
            args.seed,
        )
        records.extend(
            make_record(
                geometry,
                direction,
                direction_index,
                args.seed,
                args.scale,
                args.resolution,
            )
            for direction_index, direction in enumerate(directions)
        )
    for record in records:
        if record.hessian_rank != record.expected_hessian_rank:
            raise AssertionError("the reconstructed L65 rank changed")
        if record.disk_tangent_rank != record.expected_disk_tangent_rank:
            raise AssertionError("the general-H disk tangent rank changed")
        if (
            record.disk_curvature_rank
            != record.expected_disk_curvature_rank
        ):
            raise AssertionError(
                "the curved general-H disk tangent rank changed"
            )
        if abs(record.gain_over_disk_deficit - 1) > 1e-6:
            raise AssertionError("the full weighted face lost tightness")
        if (
            abs(record.explicit_correction_energy_ratio - 1) > 1e-6
            or abs(record.explicit_curvature_energy_ratio - 1) > 1e-6
        ):
            raise AssertionError("the explicit disk correction changed")
        if (
            record.gradient_kernel_residual > 1e-6
            or record.disk_tangent_stationarity_residual > 1e-6
        ):
            raise AssertionError("the full weighted face lost stationarity")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
