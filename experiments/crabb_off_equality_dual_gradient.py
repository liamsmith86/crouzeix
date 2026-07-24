#!/usr/bin/env python3
"""Probe the off-equality characteristic-Blaschke ambient gradient.

The candidate ``O(Q)`` estimate is tested both in the full ambient
space and after orthogonal projection onto the transported coercive
circular-normal covectors (support modes three through ``p``).
Stable ``sqrt(Q)`` scaling falsifies that estimate as a general route.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm


@dataclass(frozen=True)
class DualGradientRecord:
    """One off-equality gradient-scaling audit."""

    dimension: int
    length: int
    disk_normal_scale: float
    quartic_defect: float
    dual_deficit: float
    full_gradient_norm: float
    coercive_projection_norm: float
    full_over_sqrt_defect: float
    coercive_over_sqrt_defect: float
    coercive_normal_rank: int
    expected_coercive_normal_rank: int
    equality_gradient_control: float


def disk_model(coefficients: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return L122's coefficient-coordinate disk matrix and metric."""

    length = len(coefficients)
    dimension = length + 1
    toeplitz = np.zeros((dimension, dimension), dtype=complex)
    toeplitz[:length, :length] = np.eye(length) / 2
    for offset in range(1, length):
        for row in range(length - offset):
            toeplitz[row, row + offset] = coefficients[offset]
            toeplitz[row + offset, row] = np.conj(
                coefficients[offset]
            )

    shift = np.zeros_like(toeplitz)
    for row in range(length):
        shift[row, row + 1] = 1
    metric = toeplitz + shift.conj().T @ toeplitz @ shift
    operator = 2 * np.linalg.solve(metric, toeplitz @ shift)
    return operator, metric


def matrix_polynomial(
    coefficients: np.ndarray,
    matrix: np.ndarray,
) -> np.ndarray:
    """Evaluate a descending-coefficient polynomial at a matrix."""

    value = np.zeros_like(matrix)
    identity = np.eye(len(matrix), dtype=complex)
    for coefficient in coefficients:
        value = value @ matrix + coefficient * identity
    return value


def polynomial_frechet(
    coefficients: np.ndarray,
    matrix: np.ndarray,
    direction: np.ndarray,
) -> np.ndarray:
    """Evaluate a polynomial Frechet derivative by paired Horner steps."""

    value = np.zeros_like(matrix)
    derivative = np.zeros_like(matrix)
    identity = np.eye(len(matrix), dtype=complex)
    for coefficient in coefficients:
        derivative = derivative @ matrix + value @ direction
        value = value @ matrix + coefficient * identity
    return derivative


def support_data(
    metric: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return boundary angles, monomial vectors, and support denominators."""

    dimension = len(metric)
    angles = 2 * np.pi * np.arange(resolution) / resolution
    boundary = np.exp(1j * angles)
    monomials = boundary[:, None] ** np.arange(dimension)[None, :]
    denominator = np.real(
        np.einsum(
            "ti,ij,tj->t",
            monomials.conj(),
            metric,
            monomials,
        )
    )
    return angles, monomials, denominator


def support_fourier(
    metric: np.ndarray,
    direction: np.ndarray,
    angles: np.ndarray,
    monomials: np.ndarray,
    denominator: np.ndarray,
) -> np.ndarray:
    """Return Fourier coefficients of the first support variation."""

    numerator = np.einsum(
        "ti,ij,tj->t",
        monomials.conj(),
        metric @ direction,
        monomials,
    )
    support = (
        np.real(np.exp(-1j * angles) * numerator)
        / denominator
    )
    return np.fft.fft(support) / len(angles)


def first_riemann_correction(
    operator: np.ndarray,
    support_coefficients: np.ndarray,
    maximum_mode: int,
) -> np.ndarray:
    """Evaluate the analytic first inverse-Riemann correction."""

    ascending = np.zeros(maximum_mode + 2, dtype=complex)
    ascending[1] = support_coefficients[0]
    ascending[2:] = 2 * support_coefficients[1 : maximum_mode + 1]
    return matrix_polynomial(ascending[::-1], operator)


def equality_and_normal(length: int) -> tuple[np.ndarray, np.ndarray]:
    """Return deterministic phase-palindromic and transverse disk data."""

    equality = np.zeros(length, dtype=complex)
    normal = np.zeros(length, dtype=complex)
    for offset in range(1, length):
        partner = length - offset
        if offset < partner:
            equality[offset] = (
                0.04 + 0.015 * offset
                + 1j * (0.01 + 0.007 * offset)
            )
            equality[partner] = np.conj(equality[offset])
            normal[offset] = (
                0.03 + 0.011 * offset
                + 1j * (0.02 + 0.005 * offset)
            )
            normal[partner] = -np.conj(normal[offset])
        elif offset == partner:
            equality[offset] = 0.05
            normal[offset] = 0.04j
    return equality, normal


def quartic_defect(coefficients: np.ndarray) -> float:
    """Return ``||z||^4-|z^T Jz|^2`` on nonconstant coefficients."""

    value = coefficients[1:]
    return float(
        np.vdot(value, value).real**2
        - abs(np.dot(value, np.flip(value))) ** 2
    )


def gradient_record(
    length: int,
    scale: float,
    resolution: int,
    equality_control: float,
) -> DualGradientRecord:
    """Compute one full and coercive-projected dual gradient."""

    equality, normal = equality_and_normal(length)
    coefficients = equality + scale * normal
    operator, metric = disk_model(coefficients)
    dimension = length + 1

    characteristic = np.poly(operator)
    numerator_coefficients = characteristic[:-1]
    denominator_coefficients = np.conj(
        numerator_coefficients[::-1]
    )
    numerator = matrix_polynomial(
        numerator_coefficients,
        operator,
    )
    denominator = matrix_polynomial(
        denominator_coefficients,
        operator,
    )
    denominator_inverse = np.linalg.inv(denominator)
    blaschke = numerator @ denominator_inverse

    metric_sqrt = np.asarray(sqrtm(metric), dtype=complex)
    metric_inverse_sqrt = np.linalg.inv(metric_sqrt)
    physical_blaschke = (
        metric_sqrt @ blaschke @ metric_inverse_sqrt
    )
    left, singular_values, right_adjoint = np.linalg.svd(
        physical_blaschke
    )
    left_vector = left[:, 0]
    right_vector = right_adjoint.conj().T[:, 0]
    singular_value = singular_values[0]

    angles, monomials, support_denominator = support_data(
        metric,
        resolution,
    )
    maximum_mode = min(256, resolution // 4)
    gradient = []
    support_columns = []
    for scalar in (1, 1j):
        for row in range(dimension):
            for column in range(dimension):
                direction = np.zeros_like(operator)
                direction[row, column] = scalar
                support_coefficients = support_fourier(
                    metric,
                    direction,
                    angles,
                    monomials,
                    support_denominator,
                )
                support_columns.append(support_coefficients)
                pulled = direction - first_riemann_correction(
                    operator,
                    support_coefficients,
                    maximum_mode,
                )
                numerator_derivative = polynomial_frechet(
                    numerator_coefficients,
                    operator,
                    pulled,
                )
                denominator_derivative = polynomial_frechet(
                    denominator_coefficients,
                    operator,
                    pulled,
                )
                blaschke_derivative = (
                    numerator_derivative @ denominator_inverse
                    - blaschke
                    @ denominator_derivative
                    @ denominator_inverse
                )
                physical_derivative = (
                    metric_sqrt
                    @ blaschke_derivative
                    @ metric_inverse_sqrt
                )
                gradient.append(
                    2
                    * singular_value
                    * np.real(
                        np.vdot(
                            left_vector,
                            physical_derivative @ right_vector,
                        )
                    )
                )

    gradient_vector = np.asarray(gradient)
    support_columns_array = np.asarray(support_columns)
    normal_rows = []
    for mode in range(3, dimension + 1):
        normal_rows.extend(
            (
                support_columns_array[:, mode].real,
                support_columns_array[:, mode].imag,
            )
        )
    normal_matrix = np.asarray(normal_rows)
    normal_rank = int(np.linalg.matrix_rank(normal_matrix, tol=1e-10))
    projection = (
        normal_matrix.T
        @ np.linalg.pinv(normal_matrix @ normal_matrix.T)
        @ normal_matrix
        @ gradient_vector
    )

    defect = quartic_defect(coefficients)
    full_norm = float(np.linalg.norm(gradient_vector))
    projection_norm = float(np.linalg.norm(projection))
    defect_scale = np.sqrt(defect) if defect > 0 else np.inf
    return DualGradientRecord(
        dimension=dimension,
        length=length,
        disk_normal_scale=scale,
        quartic_defect=defect,
        dual_deficit=float(4 - singular_value**2),
        full_gradient_norm=full_norm,
        coercive_projection_norm=projection_norm,
        full_over_sqrt_defect=full_norm / defect_scale,
        coercive_over_sqrt_defect=projection_norm / defect_scale,
        coercive_normal_rank=normal_rank,
        expected_coercive_normal_rank=2 * dimension - 4,
        equality_gradient_control=equality_control,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=6)
    parser.add_argument("--resolution", type=int, default=4096)
    parser.add_argument(
        "--scales",
        type=float,
        nargs="+",
        default=(0.03, 0.01, 0.003, 0.001),
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic scaling audit."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("length range must satisfy 3 <= minimum <= maximum")
    if args.resolution < 512 or args.resolution & (args.resolution - 1):
        raise ValueError("resolution must be a power of two at least 512")
    if any(scale <= 0 for scale in args.scales):
        raise ValueError("all transverse scales must be positive")

    records = []
    for length in range(
        args.minimum_length,
        args.maximum_length + 1,
    ):
        equality_record = gradient_record(
            length,
            0.0,
            args.resolution,
            equality_control=0.0,
        )
        equality_control = equality_record.full_gradient_norm
        length_records = [
            gradient_record(
                length,
                scale,
                args.resolution,
                equality_control=equality_control,
            )
            for scale in args.scales
        ]
        if equality_control > 1e-9:
            raise AssertionError("equality ambient stationarity failed")
        if any(
            record.coercive_normal_rank
            != record.expected_coercive_normal_rank
            for record in length_records
        ):
            raise AssertionError("coercive normal projection lost rank")
        ratios = [
            record.coercive_over_sqrt_defect
            for record in length_records
        ]
        if max(ratios) / min(ratios) > 1.2:
            raise AssertionError("the projected sqrt(Q) scaling was unstable")
        records.extend(length_records)

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
