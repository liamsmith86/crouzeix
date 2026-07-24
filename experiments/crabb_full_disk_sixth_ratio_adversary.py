#!/usr/bin/env python3
"""Search the first ``p=7`` sixth-order polynomial ratio.

On the real slice, A126's candidate base polynomial and cubic response
suggest

    P_6(z) >= 9 |C_6(z)|**2.

This binary64 discovery harness evaluates ``P_6`` independently from
the canonical Stein endpoint series, samples complex directions, and
uses multistart BFGS to minimize the homogeneous ratio.  It is an
adversarial guard and sharpness locator, not a sum-of-squares proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from crabb_full_disk_weighted_face import plucker_disk_correction


LENGTH = 6
ORDER = 6


@dataclass(frozen=True)
class SixthRatioAdversaryRecord:
    """One seeded search for the sharp sixth-order ratio."""

    dimension: int
    length: int
    seed: int
    random_count: int
    optimization_start_count: int
    random_minimum: float
    optimized_minimum: float
    conjectured_minimum: float
    optimized_gap: float
    minimizing_direction: tuple[str, ...]


def extend(matrix: np.ndarray) -> np.ndarray:
    """Append one zero row and column."""

    result = np.zeros(
        (matrix.shape[0] + 1, matrix.shape[1] + 1),
        dtype=complex,
    )
    result[:-1, :-1] = matrix
    return result


def toeplitz_matrix(direction: np.ndarray) -> np.ndarray:
    """Return the Hermitian Toeplitz matrix encoded by ``direction``."""

    length = len(direction)
    matrix = np.zeros((length, length), dtype=complex)
    for offset in range(1, length):
        for row in range(length - offset):
            matrix[row, row + offset] = direction[offset]
            matrix[row + offset, row] = np.conj(direction[offset])
    return matrix


def inverse_series(
    coefficients: list[np.ndarray],
) -> list[np.ndarray]:
    """Invert a matrix series through its supplied order."""

    inverse = [np.linalg.inv(coefficients[0])]
    for degree in range(1, len(coefficients)):
        convolution = sum(
            (
                coefficients[source_degree]
                @ inverse[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            np.zeros_like(coefficients[0]),
        )
        inverse.append(-inverse[0] @ convolution)
    return inverse


def lyapunov_series(
    operator: list[np.ndarray],
    defect: list[np.ndarray],
) -> list[np.ndarray]:
    """Solve the rank-one Stein series at the nilpotent Crabb base."""

    dimension = len(operator[0])
    base = operator[0]
    metric = []
    for degree in range(len(operator)):
        forcing = sum(
            (
                np.outer(
                    defect[left_degree],
                    np.conj(defect[degree - left_degree]),
                )
                for left_degree in range(degree + 1)
            ),
            np.zeros_like(operator[0]),
        )
        for left_degree in range(degree + 1):
            for metric_degree in range(
                degree + 1 - left_degree
            ):
                right_degree = (
                    degree - left_degree - metric_degree
                )
                if (
                    left_degree == 0
                    and right_degree == 0
                    and metric_degree == degree
                ):
                    continue
                if metric_degree < len(metric):
                    forcing += (
                        operator[left_degree].conj().T
                        @ metric[metric_degree]
                        @ operator[right_degree]
                    )
        solution = np.zeros_like(operator[0])
        power = np.eye(dimension, dtype=complex)
        for _ in range(dimension):
            solution += power.conj().T @ forcing @ power
            power = power @ base
        metric.append(solution)
    return metric


def endpoint_series(
    stein_metric: list[np.ndarray],
    disk_metric: list[np.ndarray],
    endpoint: int,
    initial_value: float,
) -> list[complex]:
    """Expand one simple generalized endpoint."""

    dimension = len(stein_metric[0])
    eigenvalues: list[complex] = [initial_value]
    vectors = [np.eye(dimension, dtype=complex)[:, endpoint]]
    pencil_base = (
        stein_metric[0] - initial_value * disk_metric[0]
    )
    for degree in range(1, len(stein_metric)):
        known = np.zeros(dimension, dtype=complex)
        for pencil_degree in range(1, degree + 1):
            pencil = stein_metric[pencil_degree].copy()
            for eigenvalue_degree in range(pencil_degree + 1):
                if eigenvalue_degree >= len(eigenvalues):
                    continue
                pencil -= (
                    eigenvalues[eigenvalue_degree]
                    * disk_metric[
                        pencil_degree - eigenvalue_degree
                    ]
                )
            known += pencil @ vectors[degree - pencil_degree]
        eigenvalue = (
            known[endpoint] / disk_metric[0][endpoint, endpoint]
        )
        eigenvalues.append(eigenvalue)
        known -= eigenvalue * disk_metric[0] @ vectors[0]
        vector = np.zeros(dimension, dtype=complex)
        for index in range(dimension):
            if index != endpoint:
                vector[index] = (
                    -known[index] / pencil_base[index, index]
                )
        vectors.append(vector)
    return eigenvalues


def endpoint_delta_six(direction: np.ndarray) -> float:
    """Return ``[s^6](lambda_+ - 4 lambda_-)``."""

    dimension = LENGTH + 1
    nilpotent_shift = np.zeros(
        (dimension, dimension),
        dtype=complex,
    )
    nilpotent_shift[
        np.arange(LENGTH),
        np.arange(1, dimension),
    ] = 1
    hermitian = [
        extend(np.eye(LENGTH) / 2),
        extend(toeplitz_matrix(direction)),
        extend(plucker_disk_correction(direction)),
        *[
            np.zeros((dimension, dimension), dtype=complex)
            for _ in range(ORDER - 2)
        ],
    ]
    disk_metric = [
        coefficient
        + nilpotent_shift.conj().T
        @ coefficient
        @ nilpotent_shift
        for coefficient in hermitian
    ]
    inverse = inverse_series(disk_metric)
    operator = []
    for degree in range(ORDER + 1):
        operator.append(
            2
            * sum(
                (
                    inverse[left_degree]
                    @ hermitian[degree - left_degree]
                    @ nilpotent_shift
                    for left_degree in range(degree + 1)
                ),
                np.zeros((dimension, dimension), dtype=complex),
            )
        )
    stein_metric = lyapunov_series(
        operator,
        [coefficient[:, 0] for coefficient in hermitian],
    )
    lower = endpoint_series(
        stein_metric,
        disk_metric,
        endpoint=0,
        initial_value=0.5,
    )
    upper = endpoint_series(
        stein_metric,
        disk_metric,
        endpoint=dimension - 1,
        initial_value=2,
    )
    return float(np.real(upper[6] - 4 * lower[6]))


def cubic_form(direction: np.ndarray) -> complex:
    """Return equation (7)'s unscaled complex cubic ``C_6``."""

    coefficients = direction[1:]
    reversal = np.flip(np.conj(coefficients))
    plucker = (
        np.outer(coefficients, reversal)
        - np.outer(reversal, coefficients)
    )
    return (
        6 * coefficients[3] * plucker[0, 3]
        - 5 * coefficients[4] * plucker[0, 2]
        + 2 * coefficients[3] * plucker[1, 2]
    )


def decode(variables: np.ndarray) -> np.ndarray:
    """Decode and normalize ten real variables as one complex ray."""

    direction = np.zeros(LENGTH, dtype=complex)
    direction[1:] = variables[:5] + 1j * variables[5:]
    norm = np.linalg.norm(direction[1:])
    if norm == 0:
        direction[1] = 1
    else:
        direction[1:] /= norm
    return direction


def sixth_ratio(variables: np.ndarray) -> float:
    """Return ``P_6/|C_6|^2`` on one normalized ray."""

    direction = decode(variables)
    cubic = cubic_form(direction)
    base = -225 * endpoint_delta_six(direction) / 32
    return float(base / (abs(cubic) ** 2 + 1e-30))


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=70226)
    parser.add_argument("--random-count", type=int, default=2_000)
    parser.add_argument("--start-count", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the seeded random and multistart ratio search."""

    args = parse_args()
    if args.random_count < 1 or args.start_count < 1:
        raise ValueError("search counts must be positive")
    rng = np.random.default_rng(args.seed)
    best_value = np.inf
    best_variables = np.zeros(10)
    starts = []
    for _ in range(args.random_count):
        variables = rng.standard_normal(10)
        value = sixth_ratio(variables)
        if value < best_value:
            best_value = value
            best_variables = variables.copy()
        if len(starts) < args.start_count:
            starts.append(variables)
    random_minimum = best_value
    starts[0] = best_variables
    for variables in starts:
        result = minimize(
            sixth_ratio,
            variables,
            method="BFGS",
            options={"gtol": 1e-10, "maxiter": 1_000},
        )
        if np.isfinite(result.fun) and result.fun < best_value:
            best_value = float(result.fun)
            best_variables = np.asarray(result.x)
    minimizing_direction = decode(best_variables)
    record = SixthRatioAdversaryRecord(
        dimension=LENGTH + 1,
        length=LENGTH,
        seed=args.seed,
        random_count=args.random_count,
        optimization_start_count=args.start_count,
        random_minimum=random_minimum,
        optimized_minimum=best_value,
        conjectured_minimum=9.0,
        optimized_gap=best_value - 9,
        minimizing_direction=tuple(
            str(value) for value in minimizing_direction[1:]
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        f"{json.dumps(asdict(record), sort_keys=True)}\n",
        encoding="utf-8",
    )
    print(json.dumps(asdict(record), sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
