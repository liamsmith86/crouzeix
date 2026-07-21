#!/usr/bin/env python3
"""Probe the contraction-similarity route on general complex matrices.

For a normalized matrix ``A`` and a small outer offset ``Omega`` of its
numerical range, Theodorsen's method supplies the boundary values of the
Riemann map ``phi: Omega -> D``.  Cauchy functional calculus then constructs
``T = phi(A)``, and L21's SDP computes

    t_*(T) = min{t: I <= P <= t I, T* P T <= P}.

The experiment is deliberately certificate-heavy: it compares two boundary
resolutions, checks scalar and matrix Cauchy unitality, checks the positive
double-layer kernel and its mass, and compares Cauchy calculus with nodal
interpolation when the eigenvector frame is well conditioned.  Results are
still numerical evidence, not a proof of a uniform bound.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterator

import numpy as np

from crouzeix import crabb_matrix, nr_support
from slice_cb_sdp import solve_similarity_sdp
from theodorsen import GeneralPullback, theodorsen_map
from zero_geometry import phi_of_points


MAP_TOLERANCE = 8e-5
CAUCHY_TOLERANCE = 2e-6
DLP_MASS_TOLERANCE = 2e-3
SDP_GAP_TOLERANCE = 2e-3
SDP_SLACK_TOLERANCE = 2e-7
CASE_FAMILIES = ("dense", "near_crabb", "triangular", "near_normal")


@dataclass(frozen=True)
class MapEvaluation:
    operator: np.ndarray
    theodorsen_error: float
    scalar_unitality_error: float
    matrix_unitality_error: float
    dlp_minimum: float
    dlp_mass_error: float
    nodal_error: float
    eigenvector_condition: float
    spectral_radius: float


@dataclass(frozen=True)
class SweepRecord:
    dimension: int
    family: str
    sample: int
    seed: int
    inflate: float
    coarse_resolution: int
    fine_resolution: int
    map_stability: float
    bound: float
    dual_ratio: float
    contraction_slack: float
    solver: str
    theodorsen_error: float
    scalar_unitality_error: float
    matrix_unitality_error: float
    dlp_minimum: float
    dlp_mass_error: float
    nodal_error: float
    eigenvector_condition: float
    spectral_radius: float
    accepted: bool


def normalize_numerical_range(matrix: np.ndarray) -> np.ndarray:
    """Center and scale a matrix so its sampled numerical radius is one."""

    size = matrix.shape[0]
    centered = matrix - np.trace(matrix) * np.eye(size) / size
    angles = np.linspace(0.0, 2.0 * np.pi, 2048, endpoint=False)
    _, boundary = nr_support(centered, angles)
    radius = float(np.max(abs(boundary)))
    if radius <= 1e-14:
        raise ValueError("the sampled numerical range is degenerate")
    return centered / radius


def evaluate_map(
    matrix: np.ndarray,
    resolution: int,
    inflate: float,
) -> MapEvaluation:
    """Construct ``phi(A)`` and all diagnostics at one resolution."""

    boundary, derivative, theodorsen_error = theodorsen_map(
        matrix,
        N=resolution,
        inflate=inflate,
        iters=500,
        tol=1e-12,
    )
    pullback = GeneralPullback(boundary, derivative)
    resolvents = pullback.resolvent_stack(matrix)
    identity = np.eye(matrix.shape[0])
    matrix_unitality_error = float(
        np.linalg.norm(pullback.calc(np.ones(resolution), matrix, resolvents) - identity, 2)
    )
    operator = pullback.calc(pullback.w, matrix, resolvents)
    dlp_minimum, dlp_mass_error = pullback.dlp_certificate(matrix)

    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    eigenvector_condition = float(np.linalg.cond(eigenvectors))
    disk_nodes = phi_of_points(pullback, eigenvalues)
    spectral_radius = float(np.max(abs(disk_nodes)))
    nodal_error = float("nan")
    if eigenvector_condition < 1e8:
        interpolant = eigenvectors @ np.diag(disk_nodes) @ np.linalg.inv(eigenvectors)
        nodal_error = float(np.linalg.norm(operator - interpolant, 2))

    return MapEvaluation(
        operator=operator,
        theodorsen_error=float(theodorsen_error),
        scalar_unitality_error=float(pullback.unitality_certificate()),
        matrix_unitality_error=matrix_unitality_error,
        dlp_minimum=float(dlp_minimum),
        dlp_mass_error=float(dlp_mass_error),
        nodal_error=nodal_error,
        eigenvector_condition=eigenvector_condition,
        spectral_radius=spectral_radius,
    )


def phase_aligned_error(coarse: np.ndarray, fine: np.ndarray) -> float:
    """Return relative map error after removing an irrelevant disk rotation."""

    correlation = np.vdot(coarse, fine)
    phase = correlation / abs(correlation) if abs(correlation) else 1.0
    scale = max(1.0, float(np.linalg.norm(fine, 2)))
    return float(np.linalg.norm(fine - phase * coarse, 2) / scale)


def map_is_accepted(evaluation: MapEvaluation, stability: float) -> bool:
    """Apply the numerical map-certificate gate."""

    nodal_ok = np.isnan(evaluation.nodal_error) or evaluation.nodal_error <= CAUCHY_TOLERANCE
    return bool(
        evaluation.theodorsen_error <= 2e-10
        and evaluation.scalar_unitality_error <= CAUCHY_TOLERANCE
        and evaluation.matrix_unitality_error <= CAUCHY_TOLERANCE
        and evaluation.dlp_minimum >= -2e-7
        and evaluation.dlp_mass_error <= DLP_MASS_TOLERANCE
        and evaluation.spectral_radius < 1.0 + 2e-6
        and nodal_ok
        and stability <= MAP_TOLERANCE
    )


def generated_cases(
    minimum_dimension: int,
    maximum_dimension: int,
    per_family: int,
    seed: int,
) -> Iterator[tuple[int, str, int, np.ndarray]]:
    """Yield varied dense, near-equality, triangular, and near-normal cases."""

    generator = np.random.default_rng(seed)
    for dimension in range(minimum_dimension, maximum_dimension + 1):
        for sample in range(per_family):
            dense = generator.standard_normal((dimension, dimension)) + 1j * generator.standard_normal(
                (dimension, dimension)
            )
            yield dimension, "dense", sample, normalize_numerical_range(dense)

            perturbation = generator.standard_normal(
                (dimension, dimension)
            ) + 1j * generator.standard_normal((dimension, dimension))
            perturbation /= np.linalg.norm(perturbation, 2)
            amplitude = float(np.exp(generator.uniform(np.log(0.02), np.log(0.3))))
            near_crabb = crabb_matrix(dimension - 1) + amplitude * perturbation
            yield dimension, "near_crabb", sample, normalize_numerical_range(near_crabb)

            diagonal = np.linspace(-1.0, 1.0, dimension) + 0.2j * generator.standard_normal(
                dimension
            )
            triangular = np.diag(diagonal).astype(complex)
            triangular += np.diag(np.geomspace(2.5, 0.25, dimension - 1), 1)
            triangular += 0.03 * np.triu(
                generator.standard_normal((dimension, dimension))
                + 1j * generator.standard_normal((dimension, dimension)),
                2,
            )
            yield dimension, "triangular", sample, normalize_numerical_range(triangular)

            normal_diagonal = generator.uniform(-1.0, 1.0, dimension) + 1j * generator.uniform(
                -1.0, 1.0, dimension
            )
            near_normal = np.diag(normal_diagonal)
            near_normal += 0.08 * (
                generator.standard_normal((dimension, dimension))
                + 1j * generator.standard_normal((dimension, dimension))
            )
            yield dimension, "near_normal", sample, normalize_numerical_range(near_normal)


def parse_case_selector(value: str) -> tuple[int, str, int]:
    """Parse a reproducible ``dimension:family:sample`` case selector."""

    try:
        dimension_text, family, sample_text = value.split(":")
        dimension = int(dimension_text)
        sample = int(sample_text)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "case must have the form dimension:family:sample"
        ) from error
    if family not in CASE_FAMILIES:
        raise argparse.ArgumentTypeError(f"unknown family {family!r}")
    return dimension, family, sample


def evaluate_case(
    dimension: int,
    family: str,
    sample: int,
    matrix: np.ndarray,
    seed: int,
    inflate: float,
    coarse_resolution: int,
    maximum_resolution: int,
) -> SweepRecord:
    coarse = evaluate_map(matrix, coarse_resolution, inflate)
    fine_resolution = coarse_resolution
    fine = coarse
    stability = float("inf")
    accepted = False
    while fine_resolution < maximum_resolution:
        fine_resolution *= 2
        fine = evaluate_map(matrix, fine_resolution, inflate)
        stability = phase_aligned_error(coarse.operator, fine.operator)
        accepted = map_is_accepted(fine, stability)
        if accepted:
            break
        coarse = fine

    bound = float("nan")
    dual_ratio = float("nan")
    contraction_slack = float("nan")
    solver = "MAP_REJECTED"
    if accepted:
        try:
            certificate = solve_similarity_sdp(fine.operator)
        except RuntimeError:
            accepted = False
            solver = "SDP_FAILED"
        else:
            bound = certificate.bound
            dual_ratio = certificate.dual_ratio
            contraction_slack = certificate.contraction_slack
            solver = certificate.solver
            accepted = bool(
                abs(bound - dual_ratio) <= SDP_GAP_TOLERANCE
                and contraction_slack >= -SDP_SLACK_TOLERANCE
            )
    return SweepRecord(
        dimension=dimension,
        family=family,
        sample=sample,
        seed=seed,
        inflate=inflate,
        coarse_resolution=coarse_resolution,
        fine_resolution=fine_resolution,
        map_stability=stability,
        bound=bound,
        dual_ratio=dual_ratio,
        contraction_slack=contraction_slack,
        solver=solver,
        theodorsen_error=fine.theodorsen_error,
        scalar_unitality_error=fine.scalar_unitality_error,
        matrix_unitality_error=fine.matrix_unitality_error,
        dlp_minimum=fine.dlp_minimum,
        dlp_mass_error=fine.dlp_mass_error,
        nodal_error=fine.nodal_error,
        eigenvector_condition=fine.eigenvector_condition,
        spectral_radius=fine.spectral_radius,
        accepted=accepted,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-dimension", type=int, default=3)
    parser.add_argument("--max-dimension", type=int, default=8)
    parser.add_argument("--per-family", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--inflate", type=float, default=0.01)
    parser.add_argument(
        "--inflates",
        type=float,
        nargs="+",
        help="evaluate each selected case at several outer offsets",
    )
    parser.add_argument(
        "--case",
        type=parse_case_selector,
        action="append",
        help="restrict to a reproducible dimension:family:sample case",
    )
    parser.add_argument("--resolution", type=int, default=512)
    parser.add_argument("--max-resolution", type=int, default=4096)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    records: list[SweepRecord] = []
    output = arguments.output.open("w", encoding="utf-8") if arguments.output else None
    try:
        generated = generated_cases(
            arguments.min_dimension,
            arguments.max_dimension,
            arguments.per_family,
            arguments.seed,
        )
        selected_cases = set(arguments.case or ())
        inflates = arguments.inflates or [arguments.inflate]
        for case in generated:
            if selected_cases and case[:3] not in selected_cases:
                continue
            for inflate in inflates:
                record = evaluate_case(
                    *case,
                    seed=arguments.seed,
                    inflate=inflate,
                    coarse_resolution=arguments.resolution,
                    maximum_resolution=arguments.max_resolution,
                )
                records.append(record)
                line = json.dumps(asdict(record), allow_nan=True)
                print(line, flush=True)
                if output:
                    output.write(line + "\n")
                    output.flush()
    finally:
        if output:
            output.close()

    accepted = [record for record in records if record.accepted]
    largest = max(accepted, key=lambda record: record.bound) if accepted else None
    print(
        json.dumps(
            {
                "summary": {
                    "total": len(records),
                    "accepted": len(accepted),
                    "largest_accepted": asdict(largest) if largest else None,
                }
            },
            allow_nan=True,
        )
    )


if __name__ == "__main__":
    main()
