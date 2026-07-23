#!/usr/bin/env python3
"""Audit L138's homogeneous spectral normal form.

The checker constructs an elliptic Crabb axis from its DCT-I spectral
data, chooses a deterministic complex homogeneous defect multiplier,
and compares:

* the direct physical Stein solve;
* the exact DCT congruence in L138 (1);
* the endpoint perturbation formula in L138 (4)--(5);
* a centered finite difference of the actual condition number.

The matrix identity and perturbation calculation have direct proofs.
This high-precision-friendly floating audit is an independent regression
against orientation, conjugation, and endpoint-normalization mistakes.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import solve_discrete_lyapunov


DEFAULT_LENGTHS = (3, 5, 8)
DEFAULT_PARAMETERS = (0.1, 0.35, 0.65)


@dataclass(frozen=True)
class NormalFormRecord:
    """One deterministic normal-form audit."""

    length: int
    dimension: int
    ellipse_parameter: float
    axis_offdiagonal_error: float
    axis_endpoint_normalization_error: float
    direct_stein_residual: float
    normal_form_error: float
    endpoint_row_error: float
    first_variation_error: float
    hessian_formula_error: float


def modulus_from_nome(nome: float) -> float:
    """Return Jacobi's modulus from its nome."""

    mp_nome = mp.mpf(nome)
    return float(
        (mp.jtheta(2, 0, mp_nome) / mp.jtheta(3, 0, mp_nome)) ** 2
    )


def dct_one(length: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the orthogonal DCT-I matrix and endpoint factors."""

    dimension = length + 1
    endpoint_factors = np.ones(dimension)
    endpoint_factors[[0, length]] = 1 / np.sqrt(2)
    indices = np.arange(dimension)
    transform = (
        np.sqrt(2 / length)
        * endpoint_factors[:, None]
        * endpoint_factors[None, :]
        * np.cos(np.pi * np.outer(indices, indices) / length)
    )
    return transform, endpoint_factors


def spectral_axis_data(
    length: int,
    ellipse_parameter: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return ``T, R, X, U, D, K0, beta`` for one L117 axis."""

    dimension = length + 1
    transform, endpoint_factors = dct_one(length)
    modulus = modulus_from_nome(ellipse_parameter**2)
    parameter = modulus**2
    complement = np.sqrt(1 - parameter)
    quarter_period = float(mp.ellipk(parameter))

    nodes = np.empty(dimension)
    beta = np.empty(dimension)
    for index in range(dimension):
        argument = 2 * quarter_period * index / length
        cn = float(mp.ellipfun("cn", argument, parameter))
        dn = float(mp.ellipfun("dn", argument, parameter))
        nodes[index] = np.sqrt(modulus) * cn / dn
        beta[index] = endpoint_factors[index] * complement / dn

    coordinate_scale = ellipse_parameter ** (
        np.arange(dimension) / 2
    )
    diagonal_gauge = np.diag(coordinate_scale)
    endpoint_metric = np.diag(endpoint_factors**2)
    right_eigenvectors = (
        np.diag(1 / endpoint_factors)
        @ diagonal_gauge
        @ transform
    )
    operator = (
        right_eigenvectors
        @ np.diag(nodes)
        @ np.linalg.inv(right_eigenvectors)
    )
    return (
        operator,
        right_eigenvectors,
        nodes,
        transform,
        diagonal_gauge,
        endpoint_metric,
        beta,
    )


def normalized_axis_metric(
    transform: np.ndarray,
    diagonal_gauge: np.ndarray,
    beta: np.ndarray,
    nodes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return normalized ``P0`` and the correspondingly scaled defect."""

    cauchy = 1 / (1 - nodes[:, None] * nodes[None, :])
    spectral_gramian = beta[:, None] * cauchy * beta[None, :]
    kernel_metric = transform @ spectral_gramian @ transform.T
    inverse_gauge = np.diag(1 / np.diag(diagonal_gauge))
    axis_metric = inverse_gauge @ kernel_metric @ inverse_gauge
    scale = np.sqrt(axis_metric[0, 0])
    return axis_metric / scale**2, beta / scale


def deterministic_multiplier(length: int) -> np.ndarray:
    """Return a reproducible complex node multiplier."""

    indices = np.arange(length + 1)
    return (
        np.cos(0.37 * indices)
        + 0.3 * np.sin(0.61 * indices)
        + 1j * (
            0.4 * np.sin(0.29 * indices)
            - 0.2 * np.cos(0.53 * indices)
        )
    )


def endpoint_hessian(
    axis_levels: np.ndarray,
    multiplier_matrix: np.ndarray,
) -> tuple[float, float]:
    """Return the relative first shift and L138's quadratic coefficient."""

    first_matrix = (
        multiplier_matrix @ np.diag(axis_levels)
        + np.diag(axis_levels) @ multiplier_matrix.conj().T
    )
    quadratic_diagonal = np.sum(
        axis_levels[None, :] * np.abs(multiplier_matrix) ** 2,
        axis=1,
    )
    second_coefficients = quadratic_diagonal.astype(complex)
    for index, level in enumerate(axis_levels):
        second_coefficients[index] += sum(
            abs(first_matrix[index, other]) ** 2
            / (level - axis_levels[other])
            for other in range(len(axis_levels))
            if other != index
        )
    lower = 0
    upper = len(axis_levels) - 1
    relative_first_shift = (
        first_matrix[upper, upper] / axis_levels[upper]
        - first_matrix[lower, lower] / axis_levels[lower]
    )
    hessian_coefficient = (
        second_coefficients[upper]
        - axis_levels[upper] * second_coefficients[lower]
    )
    return (
        float(np.real(relative_first_shift)),
        float(np.real(hessian_coefficient)),
    )


def condition_number(matrix: np.ndarray) -> float:
    """Return the Hermitian spectral condition number."""

    levels = np.linalg.eigvalsh((matrix + matrix.conj().T) / 2)
    return float(levels[-1] / levels[0])


def make_record(
    length: int,
    ellipse_parameter: float,
    finite_difference_step: float,
) -> NormalFormRecord:
    """Construct and validate one normal-form regression."""

    (
        operator,
        right,
        nodes,
        transform,
        diagonal_gauge,
        endpoint_metric,
        raw_beta,
    ) = spectral_axis_data(length, ellipse_parameter)
    axis_metric, beta = normalized_axis_metric(
        transform,
        diagonal_gauge,
        raw_beta,
        nodes,
    )
    axis_levels = np.real(np.diag(axis_metric))
    axis_offdiagonal_error = float(
        np.max(np.abs(axis_metric - np.diag(axis_levels)))
    )
    axis_endpoint_normalization_error = abs(axis_levels[0] - 1)

    multiplier = deterministic_multiplier(length)
    amplitude = 0.07
    homogeneous = 1 + amplitude * multiplier
    spectral_defect = beta * homogeneous
    forcing = np.linalg.solve(right.conj().T, spectral_defect)
    direct_metric = solve_discrete_lyapunov(
        operator.conj().T,
        np.outer(forcing, forcing.conj()),
    )
    direct_stein_residual = float(
        np.max(
            np.abs(
                direct_metric
                - operator.conj().T @ direct_metric @ operator
                - np.outer(forcing, forcing.conj())
            )
        )
    )

    cauchy = 1 / (1 - nodes[:, None] * nodes[None, :])
    spectral_axis = beta[:, None] * cauchy * beta[None, :]
    kernel_metric = transform @ spectral_axis @ transform.T
    multiplier_operator = (
        transform @ np.diag(homogeneous) @ transform.T
    )
    inverse_gauge = np.diag(1 / np.diag(diagonal_gauge))
    normal_metric = (
        inverse_gauge
        @ multiplier_operator
        @ kernel_metric
        @ multiplier_operator.conj().T
        @ inverse_gauge
    )
    endpoint_inverse_sqrt = np.diag(
        1 / np.sqrt(np.diag(endpoint_metric))
    )
    direct_generalized = (
        endpoint_inverse_sqrt
        @ direct_metric
        @ endpoint_inverse_sqrt
    )
    normal_form_error = float(
        np.max(np.abs(direct_generalized - normal_metric))
    )

    homogeneous_generator = (
        transform @ np.diag(multiplier) @ transform.T
    )
    multiplier_matrix = (
        inverse_gauge @ homogeneous_generator @ diagonal_gauge
    )
    endpoint_row_error = abs(
        homogeneous_generator[0, 0]
        - homogeneous_generator[-1, -1]
    )
    first_variation, predicted_hessian = endpoint_hessian(
        axis_levels,
        multiplier_matrix,
    )

    def perturbed_metric(value: float) -> np.ndarray:
        congruence = np.eye(length + 1) + value * multiplier_matrix
        return congruence @ axis_metric @ congruence.conj().T

    base_condition = condition_number(axis_metric)

    def centered_differences(step: float) -> tuple[float, float]:
        plus_condition = condition_number(perturbed_metric(step))
        minus_condition = condition_number(perturbed_metric(-step))
        first = (plus_condition - minus_condition) / (2 * step)
        second_coefficient = (
            plus_condition + minus_condition - 2 * base_condition
        ) / (2 * step**2)
        return first, second_coefficient

    first_full, hessian_full = centered_differences(
        finite_difference_step
    )
    first_half, hessian_half = centered_differences(
        finite_difference_step / 2
    )
    first_finite_difference = (4 * first_half - first_full) / 3
    finite_hessian = (4 * hessian_half - hessian_full) / 3
    first_variation_error = max(
        abs(first_variation),
        abs(first_finite_difference),
    )
    hessian_formula_error = abs(finite_hessian - predicted_hessian)

    thresholds = {
        "axis offdiagonal": (axis_offdiagonal_error, 1e-8),
        "axis normalization": (
            axis_endpoint_normalization_error,
            2e-12,
        ),
        "Stein residual": (direct_stein_residual, 2e-11),
        "normal form": (normal_form_error, 1e-8),
        "endpoint row": (endpoint_row_error, 2e-12),
        "first variation": (first_variation_error, 2e-7),
        "Hessian formula": (hessian_formula_error, 5e-4),
    }
    failures = [
        f"{name}={value:.3e}>{threshold:.3e}"
        for name, (value, threshold) in thresholds.items()
        if value > threshold
    ]
    if failures:
        raise AssertionError("; ".join(failures))

    return NormalFormRecord(
        length=length,
        dimension=length + 1,
        ellipse_parameter=ellipse_parameter,
        axis_offdiagonal_error=axis_offdiagonal_error,
        axis_endpoint_normalization_error=(
            axis_endpoint_normalization_error
        ),
        direct_stein_residual=direct_stein_residual,
        normal_form_error=normal_form_error,
        endpoint_row_error=float(endpoint_row_error),
        first_variation_error=first_variation_error,
        hessian_formula_error=hessian_formula_error,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lengths",
        nargs="+",
        type=int,
        default=DEFAULT_LENGTHS,
    )
    parser.add_argument(
        "--ellipse-parameters",
        nargs="+",
        type=float,
        default=DEFAULT_PARAMETERS,
    )
    parser.add_argument(
        "--finite-difference-step",
        type=float,
        default=2e-4,
    )
    parser.add_argument("--precision", type=int, default=80)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic normal-form audit grid."""

    args = parse_args()
    if any(length < 2 for length in args.lengths):
        raise ValueError("all lengths must be at least two")
    if any(
        not 0 < parameter < 1
        for parameter in args.ellipse_parameters
    ):
        raise ValueError("ellipse parameters must lie in (0,1)")
    if args.precision < 40:
        raise ValueError("use at least 40 digits for elliptic data")
    if not 0 < args.finite_difference_step < 0.01:
        raise ValueError("finite-difference step must lie in (0, .01)")
    mp.mp.dps = args.precision

    records = [
        make_record(
            length,
            parameter,
            args.finite_difference_step,
        )
        for length in args.lengths
        for parameter in args.ellipse_parameters
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
