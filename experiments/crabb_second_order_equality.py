#!/usr/bin/env python3
"""Audit the L65 equality kernel against exact affine-unitary motions."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crouzeix import crabb_matrix
from general_crabb_second_order_modes import full_quadratic_form


@dataclass(frozen=True)
class EqualityRecord:
    dimension: int
    support_resolution: int
    quadratic_rank: int
    kernel_dimension: int
    orbit_dimension: int
    residual_dimension: int
    expected_residual_dimension: int
    orbit_quadratic_residual: float


def real_vector(matrix: np.ndarray) -> np.ndarray:
    """Vectorize a complex matrix in the real-form convention of L65 data."""

    return np.concatenate((matrix.real.reshape(-1), matrix.imag.reshape(-1)))


def affine_unitary_generators(base: np.ndarray) -> np.ndarray:
    """Return columns spanning unitary conjugacy, translation, and real scaling."""

    dimension = base.shape[0]
    generators = []
    for row in range(dimension):
        skew = np.zeros_like(base)
        skew[row, row] = 1j
        generators.append(base @ skew - skew @ base)
    for row in range(dimension):
        for column in range(row + 1, dimension):
            real_skew = np.zeros_like(base)
            real_skew[row, column] = 1
            real_skew[column, row] = -1
            generators.append(base @ real_skew - real_skew @ base)

            imaginary_skew = np.zeros_like(base)
            imaginary_skew[row, column] = 1j
            imaginary_skew[column, row] = 1j
            generators.append(base @ imaginary_skew - imaginary_skew @ base)
    generators.extend((np.eye(dimension), 1j * np.eye(dimension), base))
    return np.stack([real_vector(generator) for generator in generators], axis=1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for dimension in range(args.minimum_size, args.maximum_size + 1):
            support_resolution = max(64, 4 * dimension)
            form = full_quadratic_form(dimension, support_resolution)
            eigenvalues = np.linalg.eigvalsh(form)
            quadratic_rank = int(np.count_nonzero(eigenvalues < -1e-8))
            kernel_dimension = 2 * dimension**2 - quadratic_rank
            orbit = affine_unitary_generators(crabb_matrix(dimension - 1))
            orbit_dimension = int(np.linalg.matrix_rank(orbit, tol=1e-9))
            orbit_residual = float(np.linalg.norm(form @ orbit, 2))
            record = EqualityRecord(
                dimension=dimension,
                support_resolution=support_resolution,
                quadratic_rank=quadratic_rank,
                kernel_dimension=kernel_dimension,
                orbit_dimension=orbit_dimension,
                residual_dimension=kernel_dimension - orbit_dimension,
                expected_residual_dimension=2 * dimension - 2,
                orbit_quadratic_residual=orbit_residual,
            )
            if quadratic_rank != dimension * (dimension - 2):
                raise RuntimeError("unexpected L65 rank")
            if orbit_dimension != dimension**2 + 2:
                raise RuntimeError("unexpected affine-unitary orbit dimension")
            if record.residual_dimension != record.expected_residual_dimension:
                raise RuntimeError("unexpected quotient dimension")
            if orbit_residual > 1e-7:
                raise RuntimeError("an exact-orbit tangent left the L65 kernel")
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
