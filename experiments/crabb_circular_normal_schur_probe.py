#!/usr/bin/env python3
"""Probe the sharp disk-deficit/circular-normal Schur complement.

The off-equality ambient gradient is only ``O(sqrt(Q))``.  Its square
therefore lives on the same face as the quartic disk deficit and cannot
be discarded by a soft remainder estimate.  This script restricts L65's
negative second variation to the true coercive circular-normal quotient
and compares the resulting completed-square gain with the characteristic
Blaschke deficit.  A phase-palindromic equality anchor is approached with
an independent transverse disk coordinate, and both scales then tend to
the stratified Crabb apex.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_off_equality_dual_gradient import (
    disk_model,
    dual_gradient_and_coercive_rows,
    equality_and_normal,
    quartic_defect,
)
from general_crabb_second_order_modes import reduced_second_order_value


@dataclass(frozen=True)
class CircularNormalSchurRecord:
    """One completed-square ratio on a homogeneous disk-coordinate ray."""

    dimension: int
    length: int
    equality_amplitude: float
    transverse_ratio: float
    quartic_defect: float
    dual_deficit: float
    circular_normal_gradient_norm: float
    gradient_over_sqrt_quartic_defect: float
    gradient_over_amplitude_sqrt_quartic_defect: float
    circular_normal_schur_gain: float
    residual_deficit: float
    gain_over_dual_deficit: float
    residual_over_quartic_defect: float
    normal_dimension: int
    hessian_largest_eigenvalue: float
    hessian_smallest_eigenvalue: float


def real_vector_to_matrix(
    vector: np.ndarray,
    dimension: int,
) -> np.ndarray:
    """Convert the standard real ambient coordinates to a complex matrix."""

    entry_count = dimension**2
    real = vector[:entry_count].reshape((dimension, dimension))
    imaginary = vector[entry_count:].reshape((dimension, dimension))
    return real + 1j * imaginary


def coercive_basis_and_hessian(
    anchor_coefficients: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the transported normal basis and L65's restricted form."""

    length = len(anchor_coefficients)
    dimension = length + 1
    _, _, normal_rows = dual_gradient_and_coercive_rows(
        anchor_coefficients,
        resolution,
    )
    _, singular_values, right_adjoint = np.linalg.svd(
        normal_rows,
        full_matrices=False,
    )
    rank = int(np.sum(singular_values > 1e-10))
    expected_rank = 2 * dimension - 4
    if rank != expected_rank:
        raise RuntimeError(
            f"normal rank {rank} does not match {expected_rank}"
        )
    basis = right_adjoint[:rank].T
    _, metric = disk_model(anchor_coefficients)
    metric_sqrt = np.asarray(sqrtm(metric), dtype=complex)
    metric_inverse_sqrt = np.linalg.inv(metric_sqrt)
    matrices = []
    for index in range(rank):
        coefficient_direction = real_vector_to_matrix(
            basis[:, index],
            dimension,
        )
        matrices.append(
            metric_sqrt
            @ coefficient_direction
            @ metric_inverse_sqrt
        )
    support_resolution = max(64, 4 * dimension)
    hessian = np.zeros((rank, rank))
    diagonal = [
        reduced_second_order_value(matrix, support_resolution)
        for matrix in matrices
    ]
    for row in range(rank):
        hessian[row, row] = diagonal[row]
        for column in range(row):
            combined = reduced_second_order_value(
                matrices[row] + matrices[column],
                support_resolution,
            )
            entry = (
                combined - diagonal[row] - diagonal[column]
            ) / 2
            hessian[row, column] = entry
            hessian[column, row] = entry

    eigenvalues = np.linalg.eigvalsh(hessian)
    if eigenvalues[-1] >= -1e-8:
        raise RuntimeError(
            "L65 is not negative definite on the coercive normal quotient"
        )
    return basis, hessian


def normalized_split(length: int) -> tuple[np.ndarray, np.ndarray]:
    """Return deterministic unit equality and transverse directions."""

    equality, normal = equality_and_normal(length)
    equality_norm = np.linalg.norm(equality)
    normal_norm = np.linalg.norm(normal)
    if equality_norm == 0 or normal_norm == 0:
        raise RuntimeError("the deterministic split vanished")
    return equality / equality_norm, normal / normal_norm


def make_record(
    length: int,
    equality_amplitude: float,
    transverse_ratio: float,
    resolution: int,
    basis: np.ndarray,
    hessian: np.ndarray,
) -> CircularNormalSchurRecord:
    """Evaluate the leading circular-normal completion at one scale."""

    equality, normal = normalized_split(length)
    coefficients = equality_amplitude * (
        equality + transverse_ratio * normal
    )
    dual_norm_square, gradient, _ = dual_gradient_and_coercive_rows(
        coefficients,
        resolution,
    )
    normal_gradient = basis.T @ gradient
    gain = float(
        -0.25
        * normal_gradient
        @ np.linalg.solve(hessian, normal_gradient)
    )
    defect = quartic_defect(coefficients)
    dual_deficit = float(4 - dual_norm_square)
    residual = dual_deficit - gain
    gradient_norm = float(np.linalg.norm(normal_gradient))
    square_root_defect = np.sqrt(defect)
    eigenvalues = np.linalg.eigvalsh(hessian)
    return CircularNormalSchurRecord(
        dimension=length + 1,
        length=length,
        equality_amplitude=equality_amplitude,
        transverse_ratio=transverse_ratio,
        quartic_defect=defect,
        dual_deficit=dual_deficit,
        circular_normal_gradient_norm=gradient_norm,
        gradient_over_sqrt_quartic_defect=(
            gradient_norm / square_root_defect
        ),
        gradient_over_amplitude_sqrt_quartic_defect=(
            gradient_norm
            / (equality_amplitude * square_root_defect)
        ),
        circular_normal_schur_gain=gain,
        residual_deficit=residual,
        gain_over_dual_deficit=gain / dual_deficit,
        residual_over_quartic_defect=residual / defect,
        normal_dimension=basis.shape[1],
        hessian_largest_eigenvalue=float(eigenvalues[-1]),
        hessian_smallest_eigenvalue=float(eigenvalues[0]),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=6)
    parser.add_argument("--resolution", type=int, default=4096)
    parser.add_argument(
        "--transverse-ratios",
        type=float,
        nargs="+",
        default=(0.3, 0.1, 0.03),
    )
    parser.add_argument(
        "--equality-amplitudes",
        type=float,
        nargs="+",
        default=(0.08, 0.04, 0.02),
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic Schur-complement probe."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("length range must satisfy 3 <= minimum <= maximum")
    if args.resolution < 512 or args.resolution & (args.resolution - 1):
        raise ValueError("resolution must be a power of two at least 512")
    if any(scale <= 0 for scale in args.equality_amplitudes):
        raise ValueError("all equality amplitudes must be positive")
    if any(ratio <= 0 for ratio in args.transverse_ratios):
        raise ValueError("all transverse ratios must be positive")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        ):
            equality, _ = normalized_split(length)
            for equality_amplitude in args.equality_amplitudes:
                anchor = equality_amplitude * equality
                basis, hessian = coercive_basis_and_hessian(
                    anchor,
                    args.resolution,
                )
                for transverse_ratio in args.transverse_ratios:
                    record = make_record(
                        length,
                        equality_amplitude,
                        transverse_ratio,
                        args.resolution,
                        basis,
                        hessian,
                    )
                    line = json.dumps(asdict(record), sort_keys=True)
                    print(line, flush=True)
                    output.write(line + "\n")


if __name__ == "__main__":
    main()
