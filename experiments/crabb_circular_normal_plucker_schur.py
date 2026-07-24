#!/usr/bin/env python3
"""Regenerate the Crabb apex Plucker--Schur bound.

For a homogeneous Toeplitz disk direction, the quadratic ambient
gradient in each true circular-normal support mode factors through
``h wedge J conjugate(h)``.  The rows occupy disjoint exterior-square
anti-diagonals.  This checker compares the closed response and L65
curvature formulas with independent floating evaluations, and records
the strict completed-square margin below the disk deficit.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_circular_normal_quadratic_exact import (
    closed_quadratic_response,
)
from crabb_circular_normal_schur_probe import real_vector_to_matrix
from crabb_off_equality_dual_gradient import (
    disk_model,
    dual_gradient_and_coercive_rows,
)
from general_crabb_second_order_modes import reduced_second_order_value


@dataclass(frozen=True)
class PluckerSchurRecord:
    """One all-mode apex reconstruction in a fixed length."""

    dimension: int
    length: int
    active_mode_count: int
    maximum_response_residual: float
    maximum_curvature_residual: float
    maximum_row_norm_residual: float
    maximum_row_inner_product: float
    maximum_schur_eigenvalue: float
    disk_threshold: float
    minimum_curvature_margin: float
    strict_schur_bound: bool


def deterministic_direction(length: int) -> np.ndarray:
    """Return a bounded generic complex Toeplitz direction."""

    direction = np.zeros(length, dtype=complex)
    for offset in range(1, length):
        direction[offset] = (
            (2 * offset + 1) / (20 * length)
            + 1j * (offset + 1) / (25 * length)
        )
    return direction


def quadratic_normal_response(
    direction: np.ndarray,
    base_rows: np.ndarray,
    resolution: int,
    step: float,
) -> np.ndarray:
    """Extract the even quadratic response by Richardson extrapolation."""

    zero_coefficients = np.zeros_like(direction)
    _, zero_gradient, _ = dual_gradient_and_coercive_rows(
        zero_coefficients,
        resolution,
    )
    zero_response = base_rows @ zero_gradient

    def centered(current_step: float) -> np.ndarray:
        responses = []
        for sign in (1, -1):
            _, gradient, _ = dual_gradient_and_coercive_rows(
                sign * current_step * direction,
                resolution,
            )
            responses.append(base_rows @ gradient)
        return (
            responses[0] + responses[1] - 2 * zero_response
        ) / (2 * current_step**2)

    coarse = centered(step)
    fine = centered(step / 2)
    return (4 * fine - coarse) / 3


def predicted_response(
    direction: np.ndarray,
    mode: int,
) -> complex:
    """Evaluate the closed quadratic formula in binary64."""

    return complex(closed_quadratic_response(direction, mode))


def row_norm_square(length: int, mode: int) -> float:
    """Return the squared norm of either real Plucker row."""

    remainder = length - mode
    choose_three = remainder * (remainder - 1) * (remainder - 2) / 6
    return (
        128
        * (4 * mode - 1) ** 2
        * choose_three
        / length**4
    )


def normal_curvature(length: int, mode: int) -> float:
    """Return the positive L65 curvature on a raw support Riesz row."""

    remainder = length - mode
    flux = (
        (4 * mode - 1) ** 2
        * remainder
        * (remainder - 1)
        * (remainder - 2)
        / 24
    )
    lifted_null = (
        2
        * mode
        * (mode - 1)
        * (mode - 2)
        * (remainder + 1 / 4) ** 2
        / 3
    )
    return (flux + lifted_null) / length**4


def exterior_rows(length: int) -> list[np.ndarray]:
    """Return complex coefficient rows on the natural exterior basis."""

    pairs = [
        (left, right)
        for left in range(1, length)
        for right in range(left + 1, length)
    ]
    rows = []
    for mode in range(3, length - 2):
        row = np.zeros(len(pairs), dtype=complex)
        scale = 8 * (4 * mode - 1) / length**2
        for index, (left, right) in enumerate(pairs):
            if left + right == length + mode:
                row[index] = scale * (right - left)
        rows.append(row)
    return rows


def make_record(
    length: int,
    resolution: int,
    step: float,
) -> PluckerSchurRecord:
    """Compare every active support mode with the closed formulas."""

    zero = np.zeros(length, dtype=complex)
    _, _, base_rows = dual_gradient_and_coercive_rows(zero, resolution)
    _, base_metric = disk_model(zero)
    metric_sqrt = np.asarray(sqrtm(base_metric), dtype=complex)
    metric_inverse_sqrt = np.linalg.inv(metric_sqrt)
    direction = deterministic_direction(length)
    actual_response = quadratic_normal_response(
        direction,
        base_rows,
        resolution,
        step,
    )

    response_residuals = []
    curvature_residuals = []
    norm_residuals = []
    ratios = []
    margins = []
    rows = exterior_rows(length)
    for mode in range(3, length - 2):
        row_index = 2 * (mode - 3)
        predicted = predicted_response(direction, mode)
        response_residuals.extend(
            (
                abs(actual_response[row_index] - predicted.real),
                abs(actual_response[row_index + 1] - predicted.imag),
            )
        )

        predicted_curvature = normal_curvature(length, mode)
        for polarization in (0, 1):
            normal = real_vector_to_matrix(
                base_rows[row_index + polarization],
                length + 1,
            )
            physical_normal = (
                metric_sqrt @ normal @ metric_inverse_sqrt
            )
            actual_curvature = -reduced_second_order_value(
                physical_normal,
                max(64, 4 * (length + 1)),
            )
            curvature_residuals.append(
                abs(actual_curvature - predicted_curvature)
            )

        predicted_norm = row_norm_square(length, mode)
        constructed_norm = 2 * float(np.vdot(rows[mode - 3], rows[mode - 3]).real)
        norm_residuals.append(abs(constructed_norm - predicted_norm))
        ratio = predicted_norm / predicted_curvature
        ratios.append(ratio)
        margins.append(predicted_curvature - predicted_norm / 512)

    normalized_inner_products = []
    for left, row in enumerate(rows):
        for other in rows[:left]:
            normalized_inner_products.append(
                abs(np.vdot(row, other))
                / np.sqrt(np.vdot(row, row).real * np.vdot(other, other).real)
            )

    maximum_ratio = max(ratios, default=0)
    minimum_margin = min(margins, default=0)
    return PluckerSchurRecord(
        dimension=length + 1,
        length=length,
        active_mode_count=len(rows),
        maximum_response_residual=max(response_residuals, default=0),
        maximum_curvature_residual=max(curvature_residuals, default=0),
        maximum_row_norm_residual=max(norm_residuals, default=0),
        maximum_row_inner_product=max(
            normalized_inner_products,
            default=0,
        ),
        maximum_schur_eigenvalue=maximum_ratio,
        disk_threshold=512,
        minimum_curvature_margin=minimum_margin,
        strict_schur_bound=bool(
            not ratios
            or (maximum_ratio < 512 and minimum_margin > 0)
        ),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=12)
    parser.add_argument("--resolution", type=int, default=1024)
    parser.add_argument("--step", type=float, default=0.003)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the all-mode Plucker--Schur regeneration."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("length range must satisfy 3 <= minimum <= maximum")
    if args.resolution < 512 or args.resolution & (args.resolution - 1):
        raise ValueError("resolution must be a power of two at least 512")
    if not 0 < args.step < 0.02:
        raise ValueError("step must lie strictly between zero and 0.02")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        ):
            record = make_record(length, args.resolution, args.step)
            if (
                record.maximum_response_residual > 2e-5
                or record.maximum_curvature_residual > 2e-10
                or record.maximum_row_norm_residual > 2e-10
                or not record.strict_schur_bound
            ):
                raise AssertionError(
                    f"Plucker--Schur regeneration failed at length {length}: "
                    f"{record}"
                )
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
