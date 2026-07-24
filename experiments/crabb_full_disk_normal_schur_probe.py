#!/usr/bin/env python3
"""Probe the next normal Schur face after full-disk recentering.

The full strong completion of a generic Toeplitz quartic is exactly
borderline: ``crabb_full_disk_weighted_face.py`` finds that its maximizing
class is tangent to L122's exact general-H disk chart.  This script follows
that recentered disk path and measures the *next* competition:

* the optimized rank-one base deficit on the exact disk path; versus
* the gain obtained by completing the coercive circular-normal gradient.

The Riemann pullback derivative and the optimized Stein-envelope derivative
are computed independently of the characteristic-Blaschke gradient used by
the first script.  Results are exploratory finite-scale guards, not a proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import (
    eigh,
    solve_discrete_lyapunov,
)
from scipy.optimize import minimize, root

from crabb_full_disk_weighted_face import (
    FullDiskGeometry,
    audit_directions,
    disk_recentering_correction,
    full_disk_geometry,
)
from crabb_off_equality_dual_gradient import (
    first_riemann_correction,
    support_data,
    support_fourier,
)
@dataclass(frozen=True)
class FullDiskNormalSchurRecord:
    """One finite-scale recentered full-disk normal comparison."""

    dimension: int
    length: int
    direction_index: int
    seed: int
    disk_scale: float
    optimized_condition: float
    base_deficit: float
    full_gradient_norm: float
    full_gradient_kernel_residual: float
    formal_range_schur_gain: float
    formal_range_gain_over_base_deficit: float
    normal_gradient_norm: float
    normal_schur_gain: float
    gain_over_base_deficit: float
    residual_deficit: float
    residual_over_scale_six: float
    normal_dimension: int
    hessian_largest_eigenvalue: float
    optimizer_gradient_norm: float
    objective_gradient_relative_residual: float


def general_disk_model(
    hermitian: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return L122's coefficient-gauge disk operator and metric."""

    length = len(hermitian)
    dimension = length + 1
    extended = np.zeros((dimension, dimension), dtype=complex)
    extended[:length, :length] = hermitian
    shift = np.zeros_like(extended)
    for row in range(length):
        shift[row, row + 1] = 1
    metric = extended + shift.conj().T @ extended @ shift
    operator = 2 * np.linalg.solve(metric, extended @ shift)
    return operator, metric


def defect_from_variables(variables: np.ndarray) -> np.ndarray:
    """Decode a scale-fixed complex rank-one Stein defect."""

    half = len(variables) // 2
    return np.concatenate(
        ([1.0], variables[:half] + 1j * variables[half:])
    )


def rank_one_metric(
    operator: np.ndarray,
    defect: np.ndarray,
) -> np.ndarray:
    """Solve ``P-A*PA=qq*``."""

    return solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(defect, defect.conj()),
    )


def defect_coordinate_directions(dimension: int) -> list[np.ndarray]:
    """Return the real coordinate directions for a normalized defect."""

    directions = []
    for scalar in (1, 1j):
        for index in range(1, dimension):
            direction = np.zeros(dimension, dtype=complex)
            direction[index] = scalar
            directions.append(direction)
    return directions


def rank_one_condition_and_gradient(
    operator: np.ndarray,
    metric: np.ndarray,
    variables: np.ndarray,
    defect_directions: list[np.ndarray],
) -> tuple[float, np.ndarray]:
    """Evaluate the rank-one condition and its exact real gradient."""

    defect = defect_from_variables(variables)
    gramian = rank_one_metric(operator, defect)
    eigenvalues, eigenvectors = eigh(gramian, metric)
    lower_value = eigenvalues[0]
    upper_value = eigenvalues[-1]
    lower_vector = eigenvectors[:, 0]
    upper_vector = eigenvectors[:, -1]
    condition = upper_value / lower_value
    gradient = []
    for direction in defect_directions:
        derivative = solve_discrete_lyapunov(
            operator.conj().T,
            (
                np.outer(direction, defect.conj())
                + np.outer(defect, direction.conj())
            ),
        )
        lower_derivative = np.real(
            np.vdot(lower_vector, derivative @ lower_vector)
        )
        upper_derivative = np.real(
            np.vdot(upper_vector, derivative @ upper_vector)
        )
        gradient.append(
            condition
            * (
                upper_derivative / upper_value
                - lower_derivative / lower_value
            )
        )
    return float(condition), np.asarray(gradient)


def optimized_rank_one_metric(
    operator: np.ndarray,
    metric: np.ndarray,
    initial_defect: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray, float]:
    """Minimize the generalized condition on L118's rank-one branch."""

    dimension = len(operator)
    normalized = initial_defect / initial_defect[0]
    initial = np.concatenate(
        (normalized[1:].real, normalized[1:].imag)
    )

    defect_directions = defect_coordinate_directions(dimension)

    def objective_and_gradient(
        variables: np.ndarray,
    ) -> tuple[float, np.ndarray]:
        return rank_one_condition_and_gradient(
            operator,
            metric,
            variables,
            defect_directions,
        )

    result = minimize(
        objective_and_gradient,
        initial,
        jac=True,
        method="BFGS",
        options={"gtol": 1e-10, "maxiter": 2_000},
    )
    if not np.isfinite(result.fun):
        raise RuntimeError("rank-one envelope optimization failed")
    variables = np.asarray(result.x)
    condition, gradient = objective_and_gradient(variables)
    if np.linalg.norm(gradient) > 1e-9:
        refinement = root(
            lambda value: objective_and_gradient(value)[1],
            variables,
            method="lm",
            options={
                "ftol": 1e-13,
                "xtol": 1e-13,
                "gtol": 1e-13,
                "maxiter": 2_000,
            },
        )
        refined_variables = np.asarray(refinement.x)
        refined_condition, refined_gradient = objective_and_gradient(
            refined_variables
        )
        if (
            np.isfinite(refined_condition)
            and refined_condition <= condition + 1e-10
            and np.linalg.norm(refined_gradient) < np.linalg.norm(gradient)
        ):
            variables = refined_variables
            condition = refined_condition
            gradient = refined_gradient
    defect = defect_from_variables(variables)
    gramian = rank_one_metric(operator, defect)
    gradient_norm = float(np.linalg.norm(gradient))
    return condition, defect, gramian, gradient_norm


def objective_gradient_audit(length: int, seed: int) -> float:
    """Compare the analytic defect gradient with central differences."""

    rng = np.random.default_rng(seed + 10_000 + length)
    perturbation = (
        rng.standard_normal((length, length))
        + 1j * rng.standard_normal((length, length))
    )
    perturbation = (perturbation + perturbation.conj().T) / 2
    perturbation /= np.linalg.norm(perturbation)
    hermitian = np.eye(length, dtype=complex) / 2 + 0.02 * perturbation
    operator, metric = general_disk_model(hermitian)
    variables = 0.05 * rng.standard_normal(2 * length)
    directions = defect_coordinate_directions(length + 1)
    _, analytic = rank_one_condition_and_gradient(
        operator,
        metric,
        variables,
        directions,
    )
    step = 1e-5
    finite_difference = []
    for index in range(len(variables)):
        positive = variables.copy()
        negative = variables.copy()
        positive[index] += step
        negative[index] -= step
        positive_value = rank_one_condition_and_gradient(
            operator,
            metric,
            positive,
            directions,
        )[0]
        negative_value = rank_one_condition_and_gradient(
            operator,
            metric,
            negative,
            directions,
        )[0]
        finite_difference.append(
            (positive_value - negative_value) / (2 * step)
        )
    return float(
        np.linalg.norm(np.asarray(finite_difference) - analytic)
        / np.linalg.norm(analytic)
    )


def upper_gradient_and_normal_rows(
    operator: np.ndarray,
    metric: np.ndarray,
    gramian: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Differentiate the optimized envelope and the disk support modes."""

    dimension = len(operator)
    eigenvalues, eigenvectors = eigh(gramian, metric)
    lower_value = eigenvalues[0]
    upper_value = eigenvalues[-1]
    lower_vector = eigenvectors[:, 0]
    upper_vector = eigenvectors[:, -1]
    condition = upper_value / lower_value

    angles, monomials, denominator = support_data(metric, resolution)
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
                    denominator,
                )
                support_columns.append(support_coefficients)
                pulled = direction - first_riemann_correction(
                    operator,
                    support_coefficients,
                    maximum_mode,
                )
                forcing = (
                    pulled.conj().T @ gramian @ operator
                    + operator.conj().T @ gramian @ pulled
                )
                derivative = solve_discrete_lyapunov(
                    operator.conj().T,
                    forcing,
                )
                lower_derivative = np.real(
                    np.vdot(
                        lower_vector,
                        derivative @ lower_vector,
                    )
                )
                upper_derivative = np.real(
                    np.vdot(
                        upper_vector,
                        derivative @ upper_vector,
                    )
                )
                gradient.append(
                    condition
                    * (
                        upper_derivative / upper_value
                        - lower_derivative / lower_value
                    )
                )

    support_columns_array = np.asarray(support_columns)
    normal_rows = []
    for mode in range(3, dimension + 1):
        normal_rows.extend(
            (
                support_columns_array[:, mode].real,
                support_columns_array[:, mode].imag,
            )
        )
    return np.asarray(gradient), np.asarray(normal_rows)


def recentering_correction(
    geometry: FullDiskGeometry,
    direction: np.ndarray,
    resolution: int,
) -> np.ndarray:
    """Return the generic Toeplitz direction and its tight disk correction."""

    correction, stationarity_residual = disk_recentering_correction(
        geometry,
        direction,
        0.002,
        resolution,
    )
    if stationarity_residual > 5e-7:
        raise RuntimeError(
            "the full-disk weighted recentering lost stationarity"
        )
    return correction


def toeplitz_matrix(direction: np.ndarray) -> np.ndarray:
    """Return the Hermitian Toeplitz matrix with the given diagonals."""

    length = len(direction)
    toeplitz = np.zeros((length, length), dtype=complex)
    for offset in range(1, length):
        for row in range(length - offset):
            toeplitz[row, row + offset] = direction[offset]
            toeplitz[row + offset, row] = np.conj(direction[offset])
    return toeplitz


def make_record(
    length: int,
    direction_index: int,
    seed: int,
    scale: float,
    resolution: int,
    direction: np.ndarray,
    correction: np.ndarray,
    coefficient_hessian: np.ndarray,
    objective_gradient_residual: float,
) -> FullDiskNormalSchurRecord:
    """Evaluate one recentered disk point."""

    hermitian = (
        np.eye(length) / 2
        + scale * toeplitz_matrix(direction)
        + scale**2 * correction
    )
    if np.linalg.eigvalsh(hermitian)[0] <= 0:
        raise ValueError("the requested scale left the positive disk chart")
    operator, metric = general_disk_model(hermitian)
    extended = np.zeros(length + 1, dtype=complex)
    extended[:length] = hermitian[:, 0]
    condition, _, gramian, optimizer_gradient = (
        optimized_rank_one_metric(
            operator,
            metric,
            extended,
        )
    )
    gradient, normal_rows = upper_gradient_and_normal_rows(
        operator,
        metric,
        gramian,
        resolution,
    )
    full_curvature = -(coefficient_hessian + coefficient_hessian.T) / 2
    curvature_eigenvalues, curvature_eigenvectors = np.linalg.eigh(
        full_curvature
    )
    curvature_tolerance = 1e-8 * max(
        1.0,
        float(curvature_eigenvalues[-1]),
    )
    positive = curvature_eigenvalues > curvature_tolerance
    projected_gradient = (
        curvature_eigenvectors[:, positive]
        @ (curvature_eigenvectors[:, positive].T @ gradient)
    )
    full_gradient_kernel_residual = float(
        np.linalg.norm(gradient - projected_gradient)
        / max(np.finfo(float).eps, np.linalg.norm(gradient))
    )
    inverse_gradient = (
        curvature_eigenvectors[:, positive]
        @ (
            (curvature_eigenvectors[:, positive].T @ gradient)
            / curvature_eigenvalues[positive]
        )
    )
    full_gain = float(0.25 * gradient @ inverse_gradient)
    _, singular_values, right_adjoint = np.linalg.svd(
        normal_rows,
        full_matrices=False,
    )
    rank = int(np.sum(singular_values > 1e-10))
    expected_rank = 2 * (length + 1) - 4
    if rank != expected_rank:
        raise RuntimeError(
            f"normal rank {rank} does not match {expected_rank}"
        )
    basis = right_adjoint[:rank].T
    normal_gradient = basis.T @ gradient
    restricted_hessian = basis.T @ coefficient_hessian @ basis
    hessian_eigenvalues = np.linalg.eigvalsh(restricted_hessian)
    if hessian_eigenvalues[-1] >= -1e-7:
        raise RuntimeError("the normal Hessian lost strict negativity")
    gain = float(
        -0.25
        * normal_gradient
        @ np.linalg.solve(restricted_hessian, normal_gradient)
    )
    deficit = 4 - condition
    residual = deficit - gain
    return FullDiskNormalSchurRecord(
        dimension=length + 1,
        length=length,
        direction_index=direction_index,
        seed=seed,
        disk_scale=scale,
        optimized_condition=condition,
        base_deficit=deficit,
        full_gradient_norm=float(np.linalg.norm(gradient)),
        full_gradient_kernel_residual=full_gradient_kernel_residual,
        formal_range_schur_gain=full_gain,
        formal_range_gain_over_base_deficit=full_gain / deficit,
        normal_gradient_norm=float(np.linalg.norm(normal_gradient)),
        normal_schur_gain=gain,
        gain_over_base_deficit=gain / deficit,
        residual_deficit=residual,
        residual_over_scale_six=residual / scale**6,
        normal_dimension=rank,
        hessian_largest_eigenvalue=float(hessian_eigenvalues[-1]),
        optimizer_gradient_norm=optimizer_gradient,
        objective_gradient_relative_residual=(
            objective_gradient_residual
        ),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=3)
    parser.add_argument("--maximum-length", type=int, default=5)
    parser.add_argument("--resolution", type=int, default=4096)
    parser.add_argument(
        "--scales",
        type=float,
        nargs="+",
        default=(0.15, 0.1, 0.075, 0.05),
    )
    parser.add_argument("--direction-count", type=int, default=3)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the recentered full-disk normal probe."""

    args = parse_args()
    if (
        args.minimum_length < 3
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("require 3 <= minimum length <= maximum length")
    if args.resolution < 512 or args.resolution & (args.resolution - 1):
        raise ValueError("resolution must be a power of two at least 512")
    if any(scale <= 0 for scale in args.scales):
        raise ValueError("all scales must be positive")
    if args.direction_count < 1:
        raise ValueError("direction count must be positive")

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
        gradient_residual = objective_gradient_audit(
            length,
            args.seed,
        )
        if gradient_residual > 1e-6:
            raise RuntimeError("the analytic objective gradient audit failed")
        coefficient_hessian = (
            geometry.gauge_map.T
            @ geometry.physical_hessian
            @ geometry.gauge_map
        )
        for direction_index, direction in enumerate(directions):
            correction = recentering_correction(
                geometry,
                direction,
                args.resolution,
            )
            records.extend(
                make_record(
                    length,
                    direction_index,
                    args.seed,
                    scale,
                    args.resolution,
                    direction,
                    correction,
                    coefficient_hessian,
                    gradient_residual,
                )
                for scale in args.scales
            )

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
