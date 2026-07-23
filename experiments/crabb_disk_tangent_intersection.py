#!/usr/bin/env python3
"""Audit the L65 kernel against the local circular-range manifold.

Lewis--Overton characterize the tangent space to the centered disk-matrix
manifold at a Crabb block.  After adjoining complex translations, the
condition is that all nonconstant support Fourier modes except mode one
vanish.  This script builds those constraints directly and intersects them
with the independently reconstructed L65 quadratic form.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crabb_second_order_equality import affine_unitary_generators
from crouzeix import crabb_matrix
from general_crabb_second_order_modes import full_quadratic_form, mode_aggregator


@dataclass(frozen=True)
class IntersectionRecord:
    dimension: int
    quadratic_rank: int
    circular_tangent_codimension: int
    kernel_dimension: int
    kernel_tangent_intersection_dimension: int
    affine_unitary_dimension: int
    tangent_quotient_dimension: int
    soft_normal_dimension: int
    orbit_constraint_residual: float
    support_factor_residual: float


def circular_tangent_constraints(dimension: int) -> np.ndarray:
    """Return real rows cutting out the arbitrary-center disk tangent.

    If ``d=(1/sqrt(2),1,...,1,1/sqrt(2))``, the positive support
    Fourier coefficient in mode ``k`` is, up to a harmless scalar,

        sum_{m-j-1=k} d_j d_m E_jm
        + sum_{m-j-1=-k} d_j d_m conjugate(E_jm).

    Modes two through ``dimension`` must vanish.  Modes zero and one
    encode radius and center and are therefore unconstrained.
    """

    weights = np.ones(dimension)
    weights[[0, -1]] = 1 / np.sqrt(2)
    entry_count = dimension**2
    rows: list[np.ndarray] = []
    for mode in range(2, dimension + 1):
        real_coefficients = np.zeros(entry_count)
        imaginary_coefficients = np.zeros(entry_count)
        for row in range(dimension):
            for column in range(dimension):
                grade = column - row - 1
                coefficient = weights[row] * weights[column]
                index = row * dimension + column
                if grade == mode:
                    real_coefficients[index] += coefficient
                    imaginary_coefficients[index] += coefficient
                if grade == -mode:
                    real_coefficients[index] += coefficient
                    imaginary_coefficients[index] -= coefficient
        rows.append(np.concatenate((real_coefficients, np.zeros(entry_count))))
        rows.append(np.concatenate((np.zeros(entry_count), imaginary_coefficients)))
    return np.stack(rows)


def maximum_support_factor_residual(dimension: int) -> float:
    """Check that every paired mode at least two factors through L65's ``t``."""

    last = dimension - 1
    physical_weights = np.ones(dimension)
    physical_weights[[0, -1]] = 1 / np.sqrt(2)
    residuals = []
    for mode in range(2, last):
        length = last - mode
        aggregator = mode_aggregator(length, mode)
        support_row = []
        for index in range(length):
            support_row.append(
                physical_weights[index] * physical_weights[index + mode + 1]
            )
        for index in range(length + 2):
            support_row.append(
                physical_weights[index + mode - 1] * physical_weights[index]
            )
        support_row = np.asarray(support_row)
        factor, *_ = np.linalg.lstsq(aggregator.T, support_row, rcond=None)
        residuals.append(np.linalg.norm(aggregator.T @ factor - support_row))
    return float(max(residuals, default=0.0))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=7)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tolerance = 1e-7
    with args.output.open("w") as output:
        for dimension in range(args.minimum_size, args.maximum_size + 1):
            support_resolution = max(64, 4 * dimension)
            form = full_quadratic_form(dimension, support_resolution)
            constraints = circular_tangent_constraints(dimension)
            orbit = affine_unitary_generators(crabb_matrix(dimension - 1))

            ambient_dimension = 2 * dimension**2
            quadratic_rank = int(np.linalg.matrix_rank(form, tolerance))
            constraint_rank = int(np.linalg.matrix_rank(constraints, tolerance))
            stacked_rank = int(
                np.linalg.matrix_rank(np.vstack((form, constraints)), tolerance)
            )
            kernel_dimension = ambient_dimension - quadratic_rank
            intersection_dimension = ambient_dimension - stacked_rank
            orbit_dimension = int(np.linalg.matrix_rank(orbit, tolerance))
            tangent_quotient_dimension = intersection_dimension - orbit_dimension
            soft_normal_dimension = kernel_dimension - intersection_dimension
            orbit_residual = float(np.linalg.norm(constraints @ orbit, 2))
            factor_residual = maximum_support_factor_residual(dimension)

            record = IntersectionRecord(
                dimension=dimension,
                quadratic_rank=quadratic_rank,
                circular_tangent_codimension=constraint_rank,
                kernel_dimension=kernel_dimension,
                kernel_tangent_intersection_dimension=intersection_dimension,
                affine_unitary_dimension=orbit_dimension,
                tangent_quotient_dimension=tangent_quotient_dimension,
                soft_normal_dimension=soft_normal_dimension,
                orbit_constraint_residual=orbit_residual,
                support_factor_residual=factor_residual,
            )
            expected = {
                "quadratic_rank": dimension * (dimension - 2),
                "circular_tangent_codimension": 2 * dimension - 2,
                "kernel_dimension": dimension * (dimension + 2),
                "kernel_tangent_intersection_dimension": dimension**2
                + 2 * dimension
                - 2,
                "affine_unitary_dimension": dimension**2 + 2,
                "tangent_quotient_dimension": 2 * dimension - 4,
                "soft_normal_dimension": 2,
            }
            for field, value in expected.items():
                if getattr(record, field) != value:
                    raise RuntimeError(
                        f"dimension {dimension}: unexpected {field}="
                        f"{getattr(record, field)}, expected {value}"
                    )
            if orbit_residual > 1e-10 or factor_residual > 1e-12:
                raise RuntimeError(f"dimension {dimension}: tangent identity failed")

            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
