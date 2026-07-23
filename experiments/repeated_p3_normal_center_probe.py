#!/usr/bin/env python3
"""Probe nonlinear similarity descent near L100's weighted normal center.

This is a non-load-bearing regression for the exact L100/L102 asymptotics.
"""

from __future__ import annotations

import json

import numpy as np

from crouzeix import crabb_matrix
from general_similarity_sdp import (
    evaluate_map,
    phase_aligned_error,
)
from slice_cb_sdp import solve_similarity_sdp


def terminal_matrix(normal_scale: float, transverse_scale: float) -> np.ndarray:
    """Return the two-copy terminal chart at the quadratic normal center."""

    root_two = np.sqrt(2.0)
    first = np.diag([-0.75, 0.25, -0.75]).astype(complex)
    second = np.zeros((3, 3), dtype=complex)
    second[0, 2] = 1.0
    common = np.zeros((3, 3), dtype=complex)
    common[1, 0] = common[2, 1] = 1.0

    copy = np.array(
        [
            [normal_scale, transverse_scale],
            [0.0, -normal_scale],
        ],
        dtype=complex,
    )
    perturbation = np.kron(copy.conj().T, first) + np.kron(
        copy, second
    )
    center = 3.0 * normal_scale**2 / (8.0 * root_two)
    return (
        np.kron(np.eye(2), crabb_matrix(2))
        + perturbation
        + center * np.kron(np.eye(2), common)
    )


def evaluate(normal_scale: float, transverse_scale: float) -> dict[str, float]:
    """Evaluate one chart at two map resolutions and solve the metric SDP."""

    matrix = terminal_matrix(normal_scale, transverse_scale)
    coarse = evaluate_map(matrix, 2048, 0.0)
    fine = evaluate_map(matrix, 4096, 0.0)
    stability = phase_aligned_error(coarse.operator, fine.operator)
    certificate = solve_similarity_sdp(fine.operator)
    model = transverse_scale * (
        normal_scale**2 + transverse_scale**2
    )
    return {
        "normal": normal_scale,
        "transverse": transverse_scale,
        "bound": certificate.bound,
        "deficit": 4.0 - certificate.bound,
        "model": model,
        "deficit_over_model": (
            (4.0 - certificate.bound) / model if model else np.nan
        ),
        "dual_gap": certificate.bound - certificate.dual_ratio,
        "contraction_slack": certificate.contraction_slack,
        "map_stability": stability,
        "theodorsen_error": fine.theodorsen_error,
        "matrix_unitality_error": fine.matrix_unitality_error,
        "dlp_minimum": fine.dlp_minimum,
        "dlp_mass_error": fine.dlp_mass_error,
    }


def main() -> None:
    cases = (
        (0.12, 0.03),
        (0.12, 0.01),
        (0.12, 0.003),
        (0.08, 0.02),
        (0.08, 0.006),
        (0.05, 0.0125),
        (0.05, 0.004),
    )
    for normal_scale, transverse_scale in cases:
        print(
            json.dumps(evaluate(normal_scale, transverse_scale)),
            flush=True,
        )


if __name__ == "__main__":
    main()
