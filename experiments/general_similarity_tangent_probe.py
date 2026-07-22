#!/usr/bin/env python3
"""Evaluate the first-order L21 program at repeated Crabb disk blocks.

For ``A(eps) = A0 + eps E`` with ``W(A0)`` the unit disk, first-order support
perturbation theory gives

    s(theta) = lambda_max(V(theta)* Re(e^-i theta E) V(theta)),

where the columns of ``V(theta)`` span the repeated top support eigenspace.
In the disk gauge fixing zero, the inverse Riemann map has tangent

    F(w) = s_hat(0) w + 2 sum_{k>=1} s_hat(k) w^(k+1).

Since a Crabb block of size p is nilpotent of order p, only modes through
``p-2`` affect ``F(A0)``.  Thus

    phi_eps(A(eps)) = A0 + eps (E - F(A0)) + O(eps^2).

This script constructs that finite conformal tangent, evaluates the closed-form
Jensen formula for the tangent optimum, and independently solves the tangent
cone SDP obtained by linearizing the three L21 primal LMIs at the explicit
condition-two Crabb metric.  The two SDP solvers are regression checks for the
closed form.  Numerical output is not a proof of the assumed conformal-map
expansion or of any second-order neighborhood statement.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterator

import cvxpy as cp
import numpy as np

from crouzeix import crabb_matrix
from general_similarity_equality_probe import (
    COPY_PHASE_STEP,
    PERTURBATION_FAMILIES,
    generated_directions,
    matrix_sha256,
)


@dataclass(frozen=True)
class TangentRecord:
    block_size: int
    multiplicity: int
    perturbation_family: str
    direction: int
    seed: int
    perturbation_sha256: str
    support_resolution: int
    support_minimum: float
    support_maximum: float
    support_compression_gap: float
    base_support_separation: float
    fourier_coefficients: list[list[float]]
    conformal_tangent_norm: float
    base_metric_defect: float
    support_jensen_slack: float
    closed_form_derivative: float
    solver: str
    status: str
    objective_derivative: float
    predicted_drop: float
    lower_tangent_slack: float
    upper_tangent_slack: float
    contraction_tangent_slack: float


@dataclass(frozen=True)
class SupportTangent:
    values: np.ndarray
    compression_gap: float
    base_separation: float


def support_tangent(
    base: np.ndarray,
    perturbation: np.ndarray,
    multiplicity: int,
    resolution: int,
) -> SupportTangent:
    """Sample the exact repeated-eigenvalue directional support formula."""

    values = np.empty(resolution)
    minimum_compression_gap = float("inf")
    minimum_base_separation = float("inf")
    adjoint = base.conj().T
    perturbation_adjoint = perturbation.conj().T
    angles = np.linspace(0.0, 2.0 * np.pi, resolution, endpoint=False)
    for index, angle in enumerate(angles):
        phase = np.exp(-1j * angle)
        support = (phase * base + np.conj(phase) * adjoint) / 2
        perturbation_support = (
            phase * perturbation + np.conj(phase) * perturbation_adjoint
        ) / 2
        eigenvalues, eigenvectors = np.linalg.eigh(support)
        top_space = eigenvectors[:, -multiplicity:]
        compression = top_space.conj().T @ perturbation_support @ top_space
        compression_eigenvalues = np.linalg.eigvalsh(compression)
        values[index] = float(compression_eigenvalues[-1])
        if multiplicity > 1:
            minimum_compression_gap = min(
                minimum_compression_gap,
                float(compression_eigenvalues[-1] - compression_eigenvalues[-2]),
            )
        minimum_base_separation = min(
            minimum_base_separation,
            float(eigenvalues[-multiplicity] - eigenvalues[-multiplicity - 1]),
        )
    return SupportTangent(values, minimum_compression_gap, minimum_base_separation)


def conformal_operator_tangent(
    base: np.ndarray,
    perturbation: np.ndarray,
    block_size: int,
    sampled_support: np.ndarray,
) -> tuple[np.ndarray, list[complex]]:
    """Return ``E-F(A0)`` and the Fourier coefficients used in ``F``."""

    coefficients = np.fft.fft(sampled_support) / len(sampled_support)
    retained = [complex(coefficients[index]) for index in range(block_size - 1)]
    correction = retained[0].real * base
    for frequency in range(1, block_size - 1):
        correction += (
            2 * retained[frequency] * np.linalg.matrix_power(base, frequency + 1)
        )
    return perturbation - correction, retained


def crabb_metric(
    block_size: int,
    multiplicity: int,
    family: str,
) -> np.ndarray:
    """Return the explicit condition-four metric in the family's ordering."""

    crabb = crabb_matrix(block_size - 1)
    squared_weights = abs(np.diag(crabb, 1)) ** 2
    levels = np.concatenate(([1.0], np.cumprod(squared_weights)))
    level_metric = np.diag(levels)
    if family == "operator_weight":
        return np.kron(level_metric, np.eye(multiplicity)).astype(complex)
    return np.kron(np.eye(multiplicity), level_metric).astype(complex)


def level_major_matrix(
    matrix: np.ndarray,
    block_size: int,
    multiplicity: int,
    family: str,
) -> np.ndarray:
    """Conjugate a generated case to the canonical ``C_p tensor I_m`` ordering."""

    if family == "operator_weight":
        return matrix

    dimension = block_size * multiplicity
    unitary = np.zeros((dimension, dimension), dtype=complex)
    phases = np.exp(COPY_PHASE_STEP * 1j * np.arange(multiplicity))
    levels = np.arange(block_size)
    for copy, phase in enumerate(phases):
        copy_slice = slice(copy * block_size, (copy + 1) * block_size)
        unitary[copy_slice, copy_slice] = np.diag(phase ** (-levels))
    copy_major = unitary.conj().T @ matrix @ unitary
    permutation = [
        copy * block_size + level
        for level in range(block_size)
        for copy in range(multiplicity)
    ]
    return copy_major[np.ix_(permutation, permutation)]


def closed_form_tangent_derivative(
    perturbation: np.ndarray,
    tangent: np.ndarray,
    block_size: int,
    multiplicity: int,
    family: str,
    support_mean: float,
) -> tuple[float, float]:
    """Return the dual closed form and its support-function Jensen slack."""

    canonical_perturbation = level_major_matrix(
        perturbation,
        block_size,
        multiplicity,
        family,
    )
    canonical_tangent = level_major_matrix(
        tangent,
        block_size,
        multiplicity,
        family,
    )
    weights = abs(np.diag(crabb_matrix(block_size - 1), 1))
    adjacent_sum = np.zeros((multiplicity, multiplicity), dtype=complex)
    tangent_dual = np.zeros((multiplicity, multiplicity), dtype=complex)
    for level, weight in enumerate(weights):
        left = slice(level * multiplicity, (level + 1) * multiplicity)
        right = slice((level + 1) * multiplicity, (level + 2) * multiplicity)
        perturbation_block = canonical_perturbation[left, right]
        tangent_block = canonical_tangent[left, right]
        adjacent_sum += (perturbation_block + perturbation_block.conj().T) / weight
        tangent_dual += 4 * (tangent_block + tangent_block.conj().T) / weight

    adjacent_maximum = float(np.linalg.eigvalsh(adjacent_sum)[-1])
    jensen_slack = 2 * (block_size - 1) * support_mean - adjacent_maximum
    derivative = -4 * jensen_slack
    tangent_derivative = float(np.linalg.eigvalsh(tangent_dual)[-1])
    if abs(derivative - tangent_derivative) > 2e-10:
        raise RuntimeError("closed-form tangent identities are inconsistent")
    return derivative, jensen_slack


def kernel_basis(matrix: np.ndarray, tolerance: float = 1e-9) -> np.ndarray:
    """Return orthonormal columns spanning the numerical kernel of a PSD matrix."""

    hermitian = (matrix + matrix.conj().T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(hermitian)
    scale = max(1.0, float(np.max(abs(eigenvalues))))
    mask = abs(eigenvalues) <= tolerance * scale
    if not np.any(mask):
        raise ValueError("expected a nontrivial active kernel")
    return eigenvectors[:, mask]


def minimum_eigenvalue(matrix: np.ndarray) -> float:
    """Return the least eigenvalue after explicit Hermitian symmetrization."""

    hermitian = (matrix + matrix.conj().T) / 2
    return float(np.linalg.eigvalsh(hermitian)[0])


def solve_tangent_sdp(
    base: np.ndarray,
    tangent: np.ndarray,
    metric: np.ndarray,
    solver: str,
) -> tuple[str, float, float, float, float]:
    """Minimize the first-order condition-square change over metric tangents."""

    dimension = base.shape[0]
    identity = np.eye(dimension)
    lower_kernel = kernel_basis(metric - identity)
    upper_kernel = kernel_basis(4 * identity - metric)
    contraction_defect = metric - base.conj().T @ metric @ base
    contraction_kernel = kernel_basis(contraction_defect)

    metric_tangent = cp.Variable((dimension, dimension), hermitian=True)
    derivative = cp.Variable()
    contraction_tangent = (
        metric_tangent
        - base.conj().T @ metric_tangent @ base
        - tangent.conj().T @ metric @ base
        - base.conj().T @ metric @ tangent
    )
    lower_compression = lower_kernel.conj().T @ metric_tangent @ lower_kernel
    upper_compression = (
        upper_kernel.conj().T @ (derivative * identity - metric_tangent) @ upper_kernel
    )
    contraction_compression = (
        contraction_kernel.conj().T @ contraction_tangent @ contraction_kernel
    )
    problem = cp.Problem(
        cp.Minimize(derivative),
        [
            lower_compression >> 0,
            upper_compression >> 0,
            contraction_compression >> 0,
        ],
    )
    options: dict[str, float | int | bool] = {"verbose": False}
    if solver == "CLARABEL":
        options.update(
            tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=2000
        )
    elif solver == "SCS":
        options.update(eps=1e-7, max_iters=200000)
    problem.solve(solver=solver, **options)
    if metric_tangent.value is None or problem.value is None:
        raise RuntimeError(
            f"{solver} did not return a tangent solution: {problem.status}"
        )

    tangent_value = np.asarray(metric_tangent.value)
    derivative_value = float(problem.value)
    contraction_value = (
        tangent_value
        - base.conj().T @ tangent_value @ base
        - tangent.conj().T @ metric @ base
        - base.conj().T @ metric @ tangent
    )
    lower_value = lower_kernel.conj().T @ tangent_value @ lower_kernel
    upper_value = (
        upper_kernel.conj().T
        @ (derivative_value * identity - tangent_value)
        @ upper_kernel
    )
    contraction_compressed_value = (
        contraction_kernel.conj().T @ contraction_value @ contraction_kernel
    )
    return (
        str(problem.status),
        derivative_value,
        minimum_eigenvalue(lower_value),
        minimum_eigenvalue(upper_value),
        minimum_eigenvalue(contraction_compressed_value),
    )


def records_for_case(
    block_size: int,
    multiplicity: int,
    family: str,
    direction: int,
    base: np.ndarray,
    perturbation: np.ndarray,
    seed: int,
    support_resolution: int,
    solvers: list[str],
) -> Iterator[TangentRecord]:
    """Construct one conformal tangent and solve it with each requested solver."""

    support = support_tangent(base, perturbation, multiplicity, support_resolution)
    tangent, coefficients = conformal_operator_tangent(
        base,
        perturbation,
        block_size,
        support.values,
    )
    metric = crabb_metric(block_size, multiplicity, family)
    base_defect = metric - base.conj().T @ metric @ base
    base_metric_defect = max(0.0, -minimum_eigenvalue(base_defect))
    serialized_coefficients = [[value.real, value.imag] for value in coefficients]
    closed_derivative, jensen_slack = closed_form_tangent_derivative(
        perturbation,
        tangent,
        block_size,
        multiplicity,
        family,
        coefficients[0].real,
    )
    for solver in solvers:
        status, objective, lower_slack, upper_slack, contraction_slack = (
            solve_tangent_sdp(
                base,
                tangent,
                metric,
                solver,
            )
        )
        yield TangentRecord(
            block_size=block_size,
            multiplicity=multiplicity,
            perturbation_family=family,
            direction=direction,
            seed=seed,
            perturbation_sha256=matrix_sha256(perturbation),
            support_resolution=support_resolution,
            support_minimum=float(np.min(support.values)),
            support_maximum=float(np.max(support.values)),
            support_compression_gap=support.compression_gap,
            base_support_separation=support.base_separation,
            fourier_coefficients=serialized_coefficients,
            conformal_tangent_norm=float(np.linalg.norm(tangent, 2)),
            base_metric_defect=base_metric_defect,
            support_jensen_slack=jensen_slack,
            closed_form_derivative=closed_derivative,
            solver=solver,
            status=status,
            objective_derivative=objective,
            predicted_drop=-objective,
            lower_tangent_slack=lower_slack,
            upper_tangent_slack=upper_slack,
            contraction_tangent_slack=contraction_slack,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-sizes", nargs="+", type=int, default=[3, 4])
    parser.add_argument("--multiplicities", nargs="+", type=int, default=[2])
    parser.add_argument(
        "--families",
        nargs="+",
        choices=PERTURBATION_FAMILIES,
        default=["cross"],
    )
    parser.add_argument("--directions", type=int, default=1)
    parser.add_argument("--seed", type=int, default=9173401)
    parser.add_argument("--support-resolution", type=int, default=65536)
    parser.add_argument(
        "--solvers", nargs="+", choices=("CLARABEL", "SCS"), default=["CLARABEL", "SCS"]
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.support_resolution < 64:
        raise ValueError("support resolution must be at least 64")
    if args.directions < 1:
        raise ValueError("directions must be positive")
    if any(block_size < 2 for block_size in args.block_sizes):
        raise ValueError("block sizes must be at least two")
    if any(multiplicity < 2 for multiplicity in args.multiplicities):
        raise ValueError("multiplicities must be at least two")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for case in generated_directions(
            args.block_sizes,
            args.multiplicities,
            args.families,
            args.directions,
            args.seed,
        ):
            block_size, multiplicity, family, direction, base, perturbation = case
            for record in records_for_case(
                block_size,
                multiplicity,
                family,
                direction,
                base,
                perturbation,
                args.seed,
                args.support_resolution,
                args.solvers,
            ):
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")
                output.flush()


if __name__ == "__main__":
    main()
