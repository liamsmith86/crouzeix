#!/usr/bin/env python3
"""Evaluate the second-order L21 program at a single Crabb block.

For a single Crabb block the first-order value in L61 vanishes in every
direction.  This script constructs the first two conformal coefficients from
the first two support-eigenvalue variations and the Schwarz integral.  It also
samples the Riemann-pulled operator in the smooth Theodorsen gauge at
``epsilon = -2h,-h,0,h,2h`` and independently fits

    T(epsilon) = A0 + epsilon G + epsilon^2 H + O(epsilon^3),

and solves the exact second-order Schur-complement SDP obtained from the three
L21 primal inequalities in both gauges.  Two support/map resolutions and two
SDP solvers provide independent stability checks.  The output tests the
analytic reduction; finite numerical samples do not prove its universal sign.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterator

import cvxpy as cp
import numpy as np

from general_similarity_equality_probe import (
    generated_directions,
    matrix_sha256,
)
from general_similarity_sdp import MapEvaluation, evaluate_map
from general_similarity_tangent_probe import crabb_metric
from slice_cb_sdp import solve_similarity_sdp


SECOND_ORDER_FAMILIES = ("full", "operator_weight")


@dataclass(frozen=True)
class CurveFit:
    evaluations: list[MapEvaluation]
    coefficients: np.ndarray


@dataclass(frozen=True)
class SecondOrderSolution:
    status: str
    value: float
    lower_slack: float
    upper_slack: float
    contraction_slack: float
    first_order_residual: float


@dataclass(frozen=True)
class SecondOrderRecord:
    block_size: int
    perturbation_family: str
    direction: int
    seed: int
    perturbation_sha256: str
    fit_step: float
    coarse_resolution: int
    fine_resolution: int
    support_resolution: int
    base_map_error: float
    tangent_resolution_error: float
    second_order_resolution_error: float
    analytic_second_order_resolution_error: float
    cubic_coefficient_norm: float
    quartic_coefficient_norm: float
    gauge_tangent_residual: float
    gauge_constraint_error: float
    theodorsen_error: float
    scalar_unitality_error: float
    matrix_unitality_error: float
    dlp_minimum: float
    dlp_mass_error: float
    solver: str
    status: str
    predicted_second_order: float
    fitted_predicted_second_order: float
    gauge_value_error: float
    lower_second_order_slack: float
    upper_second_order_slack: float
    contraction_second_order_slack: float
    first_order_residual: float
    nonlinear_bound: float
    nonlinear_dual_ratio: float
    nonlinear_contraction_slack: float
    nonlinear_second_order: float


def minimum_eigenvalue(matrix: np.ndarray) -> float:
    """Return the least eigenvalue after Hermitian symmetrization."""

    hermitian = (matrix + matrix.conj().T) / 2
    return float(np.linalg.eigvalsh(hermitian)[0])


def polynomial_at_matrix(coefficients: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Evaluate a scalar polynomial, stored in increasing degree order."""

    value = np.zeros_like(matrix, dtype=complex)
    power = np.eye(matrix.shape[0], dtype=complex)
    for coefficient in coefficients:
        value += coefficient * power
        power = power @ matrix
    return value


def physical_conformal_coefficients(
    base: np.ndarray,
    perturbation: np.ndarray,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return analytic first/second operator coefficients in the gauge fixing zero."""

    dimension = base.shape[0]
    angles = np.linspace(0.0, 2.0 * np.pi, resolution, endpoint=False)
    first_support = np.empty(resolution)
    second_support = np.empty(resolution)
    base_adjoint = base.conj().T
    perturbation_adjoint = perturbation.conj().T
    for index, angle in enumerate(angles):
        phase = np.exp(-1j * angle)
        base_support = (phase * base + np.conj(phase) * base_adjoint) / 2
        perturbation_support = (
            phase * perturbation + np.conj(phase) * perturbation_adjoint
        ) / 2
        eigenvalues, eigenvectors = np.linalg.eigh(base_support)
        top_vector = eigenvectors[:, -1]
        first_support[index] = float(
            np.real(top_vector.conj() @ perturbation_support @ top_vector)
        )
        lower_couplings = (
            eigenvectors[:, :-1].conj().T @ perturbation_support @ top_vector
        )
        second_support[index] = float(
            np.sum(abs(lower_couplings) ** 2 / (1 - eigenvalues[:-1])).real
        )

    first_fourier = np.fft.fft(first_support) / resolution
    analytic_spectrum = np.zeros(resolution, dtype=complex)
    analytic_spectrum[0] = first_fourier[0].real
    analytic_spectrum[1 : resolution // 2] = 2 * first_fourier[1 : resolution // 2]
    boundary_schwarz = np.fft.ifft(resolution * analytic_spectrum)
    frequencies = np.fft.fftfreq(resolution, d=1 / resolution)
    support_derivative = np.fft.ifft(1j * frequencies * np.fft.fft(first_support)).real
    normal_angle_shift = boundary_schwarz.imag - support_derivative
    second_normal_data = second_support - normal_angle_shift**2 / 2
    second_fourier = np.fft.fft(second_normal_data) / resolution

    # F needs degrees through 2p-1 because its Frechet derivative at A may
    # contain A^r E A^s with r,s < p.  K(A) needs only degrees below p.
    first_coefficients = np.zeros(2 * dimension, dtype=complex)
    first_coefficients[1] = first_fourier[0].real
    for degree in range(2, 2 * dimension):
        first_coefficients[degree] = 2 * first_fourier[degree - 1]
    second_coefficients = np.zeros(dimension, dtype=complex)
    second_coefficients[1] = second_fourier[0].real
    for degree in range(2, dimension):
        second_coefficients[degree] = 2 * second_fourier[degree - 1]

    first_map_at_base = polynomial_at_matrix(
        first_coefficients[:dimension],
        base,
    )
    tangent = perturbation - first_map_at_base
    powers = [np.eye(dimension, dtype=complex)]
    for _ in range(1, dimension):
        powers.append(powers[-1] @ base)
    frechet_term = np.zeros_like(base, dtype=complex)
    for degree in range(1, 2 * dimension):
        for left_degree in range(degree):
            right_degree = degree - 1 - left_degree
            if left_degree < dimension and right_degree < dimension:
                frechet_term += first_coefficients[degree] * (
                    powers[left_degree] @ perturbation @ powers[right_degree]
                )

    derivative_coefficients = np.array(
        [
            degree * first_coefficients[degree]
            for degree in range(1, len(first_coefficients))
        ]
    )
    derivative_product = np.convolve(
        derivative_coefficients,
        first_coefficients,
    )
    second_order = (
        -frechet_term
        + polynomial_at_matrix(derivative_product[:dimension], base)
        - polynomial_at_matrix(second_coefficients, base)
    )
    return tangent, second_order


def fit_operator_curve(
    base: np.ndarray,
    perturbation: np.ndarray,
    step: float,
    resolution: int,
) -> CurveFit:
    """Fit five Taylor coefficients from a symmetric five-point map sample."""

    nodes = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    epsilons = step * nodes
    evaluations = [
        evaluate_map(base + epsilon * perturbation, resolution, 0.0)
        for epsilon in epsilons
    ]
    operators = np.asarray([evaluation.operator for evaluation in evaluations])
    vandermonde = np.vander(nodes, 5, increasing=True)
    scaled_coefficients = np.linalg.solve(
        vandermonde,
        operators.reshape(5, -1),
    ).reshape(5, *base.shape)
    powers = step ** np.arange(5)
    coefficients = scaled_coefficients / powers[:, None, None]
    return CurveFit(evaluations, coefficients)


def infinitesimal_automorphism_residual(
    base: np.ndarray,
    fitted_tangent: np.ndarray,
    physical_tangent: np.ndarray,
) -> tuple[float, float]:
    """Check that two disk gauges differ by ``a0 I + a1 A + a2 A^2``."""

    quadratic = base @ base
    basis = [np.eye(base.shape[0]), base]
    has_quadratic = np.linalg.norm(quadratic, 2) > 1e-14
    if has_quadratic:
        basis.append(quadratic)
    design = np.stack([matrix.reshape(-1) for matrix in basis], axis=1)
    coefficients, *_ = np.linalg.lstsq(
        design,
        (fitted_tangent - physical_tangent).reshape(-1),
        rcond=None,
    )
    fitted_difference = sum(
        coefficient * matrix for coefficient, matrix in zip(coefficients, basis)
    )
    residual = float(
        np.linalg.norm(
            fitted_tangent - physical_tangent - fitted_difference,
            2,
        )
    )
    constraint_error = abs(coefficients[1].real)
    if has_quadratic:
        constraint_error = max(
            constraint_error,
            abs(coefficients[2] + np.conj(coefficients[0])),
        )
    constraint_error = float(constraint_error)
    return residual, constraint_error


def solve_second_order_sdp(
    base: np.ndarray,
    tangent: np.ndarray,
    second_order: np.ndarray,
    metric: np.ndarray,
    solver: str,
) -> SecondOrderSolution:
    """Solve the second-order active-face SDP in the fitted conformal gauge."""

    dimension = base.shape[0]
    identity = np.eye(dimension)
    lower_base = metric - identity
    upper_base = 4 * identity - metric
    contraction_base = metric - base.conj().T @ metric @ base

    lower_kernel = np.array([0])
    lower_range = np.arange(1, dimension)
    upper_kernel = np.array([dimension - 1])
    upper_range = np.arange(dimension - 1)
    contraction_kernel = np.arange(1, dimension)
    contraction_range = np.array([0])

    metric_tangent = cp.Variable((dimension, dimension), hermitian=True)
    metric_second_order = cp.Variable((dimension, dimension), hermitian=True)
    bound_second_order = cp.Variable()

    contraction_tangent = (
        metric_tangent
        - base.conj().T @ metric_tangent @ base
        - tangent.conj().T @ metric @ base
        - base.conj().T @ metric @ tangent
    )
    contraction_second_order = (
        metric_second_order
        - base.conj().T @ metric_second_order @ base
        - second_order.conj().T @ metric @ base
        - base.conj().T @ metric @ second_order
        - tangent.conj().T @ metric @ tangent
        - tangent.conj().T @ metric_tangent @ base
        - base.conj().T @ metric_tangent @ tangent
    )

    lower_lmi = cp.bmat(
        [
            [
                metric_second_order[np.ix_(lower_kernel, lower_kernel)],
                metric_tangent[np.ix_(lower_kernel, lower_range)],
            ],
            [
                metric_tangent[np.ix_(lower_range, lower_kernel)],
                lower_base[np.ix_(lower_range, lower_range)],
            ],
        ]
    )
    upper_lmi = cp.bmat(
        [
            [
                bound_second_order * np.eye(1)
                - metric_second_order[np.ix_(upper_kernel, upper_kernel)],
                -metric_tangent[np.ix_(upper_kernel, upper_range)],
            ],
            [
                -metric_tangent[np.ix_(upper_range, upper_kernel)],
                upper_base[np.ix_(upper_range, upper_range)],
            ],
        ]
    )
    contraction_lmi = cp.bmat(
        [
            [
                contraction_second_order[
                    np.ix_(contraction_kernel, contraction_kernel)
                ],
                contraction_tangent[np.ix_(contraction_kernel, contraction_range)],
            ],
            [
                contraction_tangent[np.ix_(contraction_range, contraction_kernel)],
                contraction_base[np.ix_(contraction_range, contraction_range)],
            ],
        ]
    )
    first_order_kernel = contraction_tangent[
        np.ix_(contraction_kernel, contraction_kernel)
    ]
    problem = cp.Problem(
        cp.Minimize(bound_second_order),
        [
            metric_tangent[0, 0] == 0,
            metric_tangent[-1, -1] == 0,
            first_order_kernel == 0,
            lower_lmi >> 0,
            upper_lmi >> 0,
            contraction_lmi >> 0,
        ],
    )
    options: dict[str, float | int | bool] = {"verbose": False}
    if solver == "CLARABEL":
        options.update(
            tol_gap_abs=1e-10,
            tol_gap_rel=1e-10,
            tol_feas=1e-10,
            max_iter=2000,
        )
    elif solver == "SCS":
        options.update(eps=1e-7, max_iters=200000)
    problem.solve(solver=solver, **options)
    if (
        metric_tangent.value is None
        or metric_second_order.value is None
        or bound_second_order.value is None
    ):
        raise RuntimeError(
            f"{solver} did not return a second-order solution: {problem.status}"
        )

    tangent_value = np.asarray(metric_tangent.value)
    second_order_value = np.asarray(metric_second_order.value)
    bound_value = float(bound_second_order.value)
    contraction_tangent_value = (
        tangent_value
        - base.conj().T @ tangent_value @ base
        - tangent.conj().T @ metric @ base
        - base.conj().T @ metric @ tangent
    )
    contraction_second_order_value = (
        second_order_value
        - base.conj().T @ second_order_value @ base
        - second_order.conj().T @ metric @ base
        - base.conj().T @ metric @ second_order
        - tangent.conj().T @ metric @ tangent
        - tangent.conj().T @ tangent_value @ base
        - base.conj().T @ tangent_value @ tangent
    )
    lower_value = np.block(
        [
            [
                second_order_value[np.ix_(lower_kernel, lower_kernel)],
                tangent_value[np.ix_(lower_kernel, lower_range)],
            ],
            [
                tangent_value[np.ix_(lower_range, lower_kernel)],
                lower_base[np.ix_(lower_range, lower_range)],
            ],
        ]
    )
    upper_value = np.block(
        [
            [
                bound_value * np.eye(1)
                - second_order_value[np.ix_(upper_kernel, upper_kernel)],
                -tangent_value[np.ix_(upper_kernel, upper_range)],
            ],
            [
                -tangent_value[np.ix_(upper_range, upper_kernel)],
                upper_base[np.ix_(upper_range, upper_range)],
            ],
        ]
    )
    contraction_value = np.block(
        [
            [
                contraction_second_order_value[
                    np.ix_(contraction_kernel, contraction_kernel)
                ],
                contraction_tangent_value[
                    np.ix_(contraction_kernel, contraction_range)
                ],
            ],
            [
                contraction_tangent_value[
                    np.ix_(contraction_range, contraction_kernel)
                ],
                contraction_base[np.ix_(contraction_range, contraction_range)],
            ],
        ]
    )
    return SecondOrderSolution(
        status=str(problem.status),
        value=bound_value,
        lower_slack=minimum_eigenvalue(lower_value),
        upper_slack=minimum_eigenvalue(upper_value),
        contraction_slack=minimum_eigenvalue(contraction_value),
        first_order_residual=float(
            np.linalg.norm(
                contraction_tangent_value[
                    np.ix_(contraction_kernel, contraction_kernel)
                ],
                2,
            )
        ),
    )


def map_diagnostics(
    evaluations: list[MapEvaluation],
) -> tuple[float, float, float, float, float]:
    """Aggregate the map diagnostics over a fitted five-point curve."""

    return (
        max(evaluation.theodorsen_error for evaluation in evaluations),
        max(evaluation.scalar_unitality_error for evaluation in evaluations),
        max(evaluation.matrix_unitality_error for evaluation in evaluations),
        min(evaluation.dlp_minimum for evaluation in evaluations),
        max(evaluation.dlp_mass_error for evaluation in evaluations),
    )


def records_for_case(
    block_size: int,
    family: str,
    direction: int,
    base: np.ndarray,
    perturbation: np.ndarray,
    seed: int,
    step: float,
    coarse_resolution: int,
    fine_resolution: int,
    support_resolution: int,
    solvers: list[str],
) -> Iterator[SecondOrderRecord]:
    """Fit one map curve and solve its second-order SDP with each solver."""

    coarse = fit_operator_curve(base, perturbation, step, coarse_resolution)
    fine = fit_operator_curve(base, perturbation, step, fine_resolution)
    fitted_base, tangent, second_order = fine.coefficients[:3]

    physical_tangent, physical_second_order = physical_conformal_coefficients(
        base,
        perturbation,
        support_resolution,
    )
    _, coarse_physical_second_order = physical_conformal_coefficients(
        base,
        perturbation,
        support_resolution // 2,
    )
    gauge_residual, gauge_constraint_error = infinitesimal_automorphism_residual(
        base,
        tangent,
        physical_tangent,
    )
    diagnostics = map_diagnostics(fine.evaluations)
    nonlinear = solve_similarity_sdp(fine.evaluations[3].operator)
    nonlinear_second_order = (nonlinear.bound - 4) / step**2
    metric = crabb_metric(block_size, 1, "operator_weight")

    for solver in solvers:
        solution = solve_second_order_sdp(
            base,
            physical_tangent,
            physical_second_order,
            metric,
            solver,
        )
        fitted_solution = solve_second_order_sdp(
            base,
            tangent,
            second_order,
            metric,
            solver,
        )
        yield SecondOrderRecord(
            block_size=block_size,
            perturbation_family=family,
            direction=direction,
            seed=seed,
            perturbation_sha256=matrix_sha256(perturbation),
            fit_step=step,
            coarse_resolution=coarse_resolution,
            fine_resolution=fine_resolution,
            support_resolution=support_resolution,
            base_map_error=float(np.linalg.norm(fitted_base - base, 2)),
            tangent_resolution_error=float(
                np.linalg.norm(tangent - coarse.coefficients[1], 2)
            ),
            second_order_resolution_error=float(
                np.linalg.norm(second_order - coarse.coefficients[2], 2)
            ),
            analytic_second_order_resolution_error=float(
                np.linalg.norm(
                    physical_second_order - coarse_physical_second_order,
                    2,
                )
            ),
            cubic_coefficient_norm=float(np.linalg.norm(fine.coefficients[3], 2)),
            quartic_coefficient_norm=float(np.linalg.norm(fine.coefficients[4], 2)),
            gauge_tangent_residual=gauge_residual,
            gauge_constraint_error=gauge_constraint_error,
            theodorsen_error=diagnostics[0],
            scalar_unitality_error=diagnostics[1],
            matrix_unitality_error=diagnostics[2],
            dlp_minimum=diagnostics[3],
            dlp_mass_error=diagnostics[4],
            solver=solver,
            status=solution.status,
            predicted_second_order=solution.value,
            fitted_predicted_second_order=fitted_solution.value,
            gauge_value_error=abs(solution.value - fitted_solution.value),
            lower_second_order_slack=solution.lower_slack,
            upper_second_order_slack=solution.upper_slack,
            contraction_second_order_slack=solution.contraction_slack,
            first_order_residual=solution.first_order_residual,
            nonlinear_bound=nonlinear.bound,
            nonlinear_dual_ratio=nonlinear.dual_ratio,
            nonlinear_contraction_slack=nonlinear.contraction_slack,
            nonlinear_second_order=nonlinear_second_order,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-sizes", nargs="+", type=int, default=[3, 4])
    parser.add_argument(
        "--families",
        nargs="+",
        choices=SECOND_ORDER_FAMILIES,
        default=list(SECOND_ORDER_FAMILIES),
    )
    parser.add_argument("--directions", type=int, default=1)
    parser.add_argument("--seed", type=int, default=9173401)
    parser.add_argument("--fit-step", type=float, default=0.005)
    parser.add_argument("--coarse-resolution", type=int, default=2048)
    parser.add_argument("--fine-resolution", type=int, default=4096)
    parser.add_argument("--support-resolution", type=int, default=256)
    parser.add_argument(
        "--solvers",
        nargs="+",
        choices=("CLARABEL", "SCS"),
        default=["CLARABEL", "SCS"],
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if any(block_size < 2 for block_size in args.block_sizes):
        raise ValueError("block sizes must be at least two")
    if args.directions < 1:
        raise ValueError("directions must be positive")
    if args.fit_step <= 0:
        raise ValueError("fit step must be positive")
    if args.coarse_resolution >= args.fine_resolution:
        raise ValueError("coarse resolution must be smaller than fine resolution")
    if args.support_resolution < 8 * max(args.block_sizes):
        raise ValueError(
            "support resolution must be at least eight times the largest block"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for case in generated_directions(
            args.block_sizes,
            [1],
            args.families,
            args.directions,
            args.seed,
        ):
            block_size, _, family, direction, base, perturbation = case
            for record in records_for_case(
                block_size,
                family,
                direction,
                base,
                perturbation,
                args.seed,
                args.fit_step,
                args.coarse_resolution,
                args.fine_resolution,
                args.support_resolution,
                args.solvers,
            ):
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")
                output.flush()


if __name__ == "__main__":
    main()
