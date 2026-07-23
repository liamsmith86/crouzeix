#!/usr/bin/env python3
"""Audit the L119 upper/lower touching gradient on elliptic Crabb axes.

For the exact L117 axis operator, this script differentiates independently:

* the condition number of its rank-one Stein Gramian, with the defect vector
  held fixed; and
* the squared norm of the matching degree-``p-1`` Chebyshev--Blaschke
  product.

The sandwich theorem says the two derivatives agree in every matrix
direction.  No conformal-map differentiation or SDP is used here.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.special import ellipj, ellipk

from crabb_elliptic_axis import (
    conformal_crabb_operator,
    elliptic_modulus_from_nome,
    explicit_all_size_metric,
    predicted_similarity_square,
)


@dataclass(frozen=True)
class TouchingRecord:
    dimension: int
    ellipse_parameter: float
    directions: int
    value_error: float
    reversal_error: float
    polynomial_descent_error: float
    maximum_derivative_error: float


def chebyshev_blaschke_zeros(
    dimension: int, ellipse_parameter: float
) -> np.ndarray:
    """Return the real zeros of the elliptic Chebyshev--Blaschke product."""

    length = dimension - 1
    modulus = elliptic_modulus_from_nome(ellipse_parameter**2)
    quarter_period = float(ellipk(modulus**2))
    arguments = (
        (2 * np.arange(length) + 1) * quarter_period / length
    )
    return np.sqrt(modulus) * ellipj(
        quarter_period - arguments, modulus**2
    )[0]


def blaschke_value_and_derivative(
    operator: np.ndarray,
    perturbation: np.ndarray,
    zeros: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate a real-zero Blaschke product and its Fréchet derivative."""

    dimension = operator.shape[0]
    identity = np.eye(dimension)
    factors: list[np.ndarray] = []
    derivatives: list[np.ndarray] = []
    for zero in zeros:
        resolvent = np.linalg.inv(identity - zero * operator)
        factor = (operator - zero * identity) @ resolvent
        derivative = perturbation @ resolvent
        derivative += (
            (operator - zero * identity)
            @ resolvent
            @ (zero * perturbation)
            @ resolvent
        )
        factors.append(factor)
        derivatives.append(derivative)

    prefixes = [identity]
    for factor in factors:
        prefixes.append(prefixes[-1] @ factor)
    suffixes = [identity for _ in range(len(factors) + 1)]
    for index in range(len(factors) - 1, -1, -1):
        suffixes[index] = factors[index] @ suffixes[index + 1]

    product_derivative = sum(
        (
            prefixes[index]
            @ derivatives[index]
            @ suffixes[index + 1]
            for index in range(len(factors))
        ),
        start=np.zeros_like(operator, dtype=complex),
    )
    return prefixes[-1], product_derivative


def condition_derivative(
    operator: np.ndarray,
    perturbation: np.ndarray,
    metric: np.ndarray,
) -> float:
    """Differentiate the fixed-defect Stein condition number."""

    forcing = (
        perturbation.conj().T @ metric @ operator
        + operator.conj().T @ metric @ perturbation
    )
    metric_derivative = solve_discrete_lyapunov(
        operator.conj().T, forcing
    )
    eigenvalues, eigenvectors = np.linalg.eigh(metric)
    lower = eigenvectors[:, 0]
    upper = eigenvectors[:, -1]
    lower_derivative = np.vdot(
        lower, metric_derivative @ lower
    ).real
    upper_derivative = np.vdot(
        upper, metric_derivative @ upper
    ).real
    return float(
        (
            upper_derivative * eigenvalues[0]
            - eigenvalues[-1] * lower_derivative
        )
        / eigenvalues[0] ** 2
    )


def norm_square_derivative(
    value: np.ndarray, derivative: np.ndarray
) -> float:
    """Differentiate a simple largest singular value squared."""

    left, singular_values, right_adjoint = np.linalg.svd(value)
    right = right_adjoint.conj().T[:, 0]
    return float(
        2
        * singular_values[0]
        * np.vdot(left[:, 0], derivative @ right).real
    )


def audit_case(
    dimension: int,
    ellipse_parameter: float,
    directions: int,
    generator: np.random.Generator,
) -> TouchingRecord:
    """Audit one size/parameter pair."""

    operator, _, _ = conformal_crabb_operator(
        dimension, ellipse_parameter
    )
    diagonal = explicit_all_size_metric(
        dimension, ellipse_parameter
    )
    metric = np.diag(diagonal)
    zeros = chebyshev_blaschke_zeros(
        dimension, ellipse_parameter
    )
    value, _ = blaschke_value_and_derivative(
        operator, np.zeros_like(operator), zeros
    )
    predicted = predicted_similarity_square(
        dimension, ellipse_parameter
    )
    value_error = max(
        abs(np.linalg.norm(value, 2) ** 2 - predicted),
        abs(np.linalg.cond(metric) - predicted),
    )

    length = dimension - 1
    amplitude = np.sqrt(
        elliptic_modulus_from_nome(
            ellipse_parameter ** (2 * length)
        )
    )
    diagonal_similarity = ellipse_parameter ** (
        np.arange(dimension) / 2
    )
    reversal = np.fliplr(np.eye(dimension))
    expected = (
        amplitude
        * diagonal_similarity[:, None]
        * reversal
        / diagonal_similarity[None, :]
    )
    reversal_error = min(
        np.linalg.norm(value - expected, 2),
        np.linalg.norm(value + expected, 2),
    )

    crabb_axis = np.zeros_like(operator)
    root_two = np.sqrt(2.0)
    for index in range(length):
        weight = root_two if index in (0, length - 1) else 1.0
        crabb_axis[index, index + 1] = weight
        crabb_axis[index + 1, index] = ellipse_parameter * weight
    chebyshev_argument = crabb_axis / (
        2 * np.sqrt(ellipse_parameter)
    )
    previous = np.eye(dimension)
    current = chebyshev_argument
    for _ in range(2, length + 1):
        previous, current = (
            current,
            2 * chebyshev_argument @ current - previous,
        )
    chebyshev_value = current if length > 0 else previous
    descended = (
        2 * ellipse_parameter ** (length / 2) * chebyshev_value
    )
    expected_descent = (
        2
        * ellipse_parameter ** (length / 2)
        * diagonal_similarity[:, None]
        * reversal
        / diagonal_similarity[None, :]
    )
    polynomial_descent_error = np.linalg.norm(
        descended - expected_descent, 2
    )

    maximum_derivative_error = 0.0
    for _ in range(directions):
        perturbation = (
            generator.standard_normal((dimension, dimension))
            + 1j
            * generator.standard_normal((dimension, dimension))
        )
        perturbation /= np.linalg.norm(perturbation)
        _, product_derivative = blaschke_value_and_derivative(
            operator, perturbation, zeros
        )
        upper_derivative = condition_derivative(
            operator, perturbation, metric
        )
        lower_derivative = norm_square_derivative(
            value, product_derivative
        )
        scale = max(
            1.0, abs(upper_derivative), abs(lower_derivative)
        )
        maximum_derivative_error = max(
            maximum_derivative_error,
            abs(upper_derivative - lower_derivative) / scale,
        )

    return TouchingRecord(
        dimension=dimension,
        ellipse_parameter=ellipse_parameter,
        directions=directions,
        value_error=float(value_error),
        reversal_error=float(reversal_error),
        polynomial_descent_error=float(polynomial_descent_error),
        maximum_derivative_error=float(maximum_derivative_error),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directions", type=int, default=8)
    parser.add_argument("--seed", type=int, default=70223)
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSONL path for audit records",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generator = np.random.default_rng(args.seed)
    records = [
        audit_case(
            dimension,
            ellipse_parameter,
            args.directions,
            generator,
        )
        for dimension in range(3, 11)
        for ellipse_parameter in (0.15, 0.35, 0.5)
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
