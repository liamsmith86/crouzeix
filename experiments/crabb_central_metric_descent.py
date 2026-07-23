#!/usr/bin/env python3
"""Probe exact-looking metric descent on central Crabb equality families.

L126 proves that the degree-``k`` Chebyshev--Blaschke image of the
``(2k+1)``-dimensional central family has an exact reducing size-three
outer block at parameter ``c**k``.  This script compares the optimal
rank-one Stein envelope and the unrestricted similarity SDP of the full
operator with those of that size-three block.

The equality under test is evidence for the missing upper metric lift;
the floating-point experiment is not a proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, solve_discrete_lyapunov

from crabb_palindromic_elliptic_face import (
    disk_matrix,
    ellipse_pullback,
    palindromic_direction,
    rank_one_envelope,
)
from slice_cb_sdp import solve_similarity_sdp


DEFAULT_DEGREES = (2, 3, 4)
DEFAULT_ELLIPSE_PARAMETERS = (0.15, 0.3, 0.5)
DEFAULT_AMPLITUDES = (0.03, 0.08)


@dataclass(frozen=True)
class CentralMetricDescentRecord:
    degree: int
    dimension: int
    ellipse_parameter: float
    descended_parameter: float
    amplitude: float
    full_rank_one_bound: float
    outer_rank_one_bound: float
    rank_one_gap: float
    full_sdp_bound: float
    outer_sdp_bound: float
    sdp_gap: float
    outer_metric_compression_residual: float
    outer_metric_cross_residual: float
    inner_spectrum_lower_slack: float
    inner_spectrum_upper_slack: float


def central_coordinate_gramian(
    degree: int,
    amplitude: float,
) -> np.ndarray:
    """Return the exact central Toeplitz-chart coordinate Gramian."""

    dimension = 2 * degree + 1
    gramian = np.diag([0.5, *([1.0] * (dimension - 2)), 0.5])
    for row in range(degree):
        column = row + degree
        gramian[row, column] += amplitude
        gramian[column, row] += amplitude
        gramian[row + 1, column + 1] += amplitude
        gramian[column + 1, row + 1] += amplitude
    return gramian


def bound_pair(
    degree: int,
    ellipse_parameter: float,
    amplitude: float,
) -> CentralMetricDescentRecord:
    """Compare one full central family with its size-three outer block."""

    dimension = 2 * degree + 1
    descended_parameter = ellipse_parameter**degree

    def bounds_and_metric(
        current_dimension: int,
        first_offset: int,
        current_parameter: float,
    ) -> tuple[float, float, np.ndarray, np.ndarray]:
        coefficients = palindromic_direction(
            current_dimension - 1,
            first_offset,
        )
        disk_operator = disk_matrix(
            current_dimension,
            amplitude,
            coefficients,
        )
        operator = ellipse_pullback(disk_operator, current_parameter)
        rank_one_bound, tail = rank_one_envelope(operator)
        sdp_bound = solve_similarity_sdp(operator).bound
        defect = np.concatenate(([1.0], tail))
        physical_metric = solve_discrete_lyapunov(
            operator.T,
            np.outer(defect, defect),
        )
        coordinate_gramian = central_coordinate_gramian(
            first_offset,
            amplitude,
        )
        eigenvalues, eigenvectors = np.linalg.eigh(coordinate_gramian)
        square_root = (
            eigenvectors * np.sqrt(eigenvalues)
        ) @ eigenvectors.T
        coordinate_metric = (
            square_root @ physical_metric @ square_root
        )
        coordinate_metric /= coordinate_metric[0, 0]
        return (
            rank_one_bound,
            sdp_bound,
            coordinate_metric,
            coordinate_gramian,
        )

    full_rank_one, full_sdp, full_metric, full_gramian = (
        bounds_and_metric(
            dimension,
            degree,
            ellipse_parameter,
        )
    )
    outer_rank_one, outer_sdp, outer_metric, _ = bounds_and_metric(
        3,
        1,
        descended_parameter,
    )
    outer_indices = [0, degree, 2 * degree]
    inner_indices = [
        index
        for index in range(dimension)
        if index not in outer_indices
    ]
    outer_compression = full_metric[np.ix_(outer_indices, outer_indices)]
    compression_residual = float(
        np.max(np.abs(outer_compression - outer_metric))
    )
    cross_residual = float(
        np.max(np.abs(full_metric[np.ix_(outer_indices, inner_indices)]))
    )
    full_generalized_spectrum = eigh(
        full_metric,
        full_gramian,
        eigvals_only=True,
    )
    outer_generalized_spectrum = eigh(
        outer_compression,
        full_gramian[np.ix_(outer_indices, outer_indices)],
        eigvals_only=True,
    )
    inner_generalized_spectrum = eigh(
        full_metric[np.ix_(inner_indices, inner_indices)],
        full_gramian[np.ix_(inner_indices, inner_indices)],
        eigvals_only=True,
    )
    lower_slack = float(
        inner_generalized_spectrum[0]
        - outer_generalized_spectrum[0]
    )
    upper_slack = float(
        outer_generalized_spectrum[-1]
        - inner_generalized_spectrum[-1]
    )
    record = CentralMetricDescentRecord(
        degree=degree,
        dimension=dimension,
        ellipse_parameter=ellipse_parameter,
        descended_parameter=descended_parameter,
        amplitude=amplitude,
        full_rank_one_bound=full_rank_one,
        outer_rank_one_bound=outer_rank_one,
        rank_one_gap=full_rank_one - outer_rank_one,
        full_sdp_bound=full_sdp,
        outer_sdp_bound=outer_sdp,
        sdp_gap=full_sdp - outer_sdp,
        outer_metric_compression_residual=compression_residual,
        outer_metric_cross_residual=cross_residual,
        inner_spectrum_lower_slack=lower_slack,
        inner_spectrum_upper_slack=upper_slack,
    )
    if abs(record.rank_one_gap) > 1e-9:
        raise AssertionError("the central rank-one descent separated")
    if abs(record.sdp_gap) > 1e-6:
        raise AssertionError("the central SDP descent separated")
    if record.outer_metric_compression_residual > 1e-5:
        raise AssertionError("the outer metric compression separated")
    if record.outer_metric_cross_residual > 1e-5:
        raise AssertionError("the optimized metric did not reduce the outer block")
    if min(lower_slack, upper_slack) < -1e-6:
        raise AssertionError("an inner metric level escaped the outer interval")
    if abs(
        full_generalized_spectrum[-1] / full_generalized_spectrum[0]
        - full_rank_one
    ) > 1e-7:
        raise AssertionError("the regenerated full metric condition drifted")
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--degrees",
        type=int,
        nargs="+",
        default=DEFAULT_DEGREES,
    )
    parser.add_argument(
        "--ellipse-parameters",
        type=float,
        nargs="+",
        default=DEFAULT_ELLIPSE_PARAMETERS,
    )
    parser.add_argument(
        "--amplitudes",
        type=float,
        nargs="+",
        default=DEFAULT_AMPLITUDES,
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic comparison grid."""

    args = parse_args()
    if any(degree < 2 for degree in args.degrees):
        raise ValueError("every descent degree must be at least two")
    records = [
        bound_pair(degree, parameter, amplitude)
        for degree in args.degrees
        for parameter in args.ellipse_parameters
        for amplitude in args.amplitudes
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
