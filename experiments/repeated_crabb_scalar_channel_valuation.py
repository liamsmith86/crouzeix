#!/usr/bin/env python3
"""Audit the leading scalar-channel leakage at a monomial inner apex.

The proof is the covariance perturbation in
``proof/repeated_crabb_scalar_channel_valuation.md``.  This script
reconstructs exact matrix-inner paths from matrix Schur parameters and
uses a deterministic global search on the two-dimensional input sphere
only as an independent diagnostic.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from crabb_block_hardy_equality import format_float
from repeated_crabb_matrix_schur_chart import reconstruct_series


MatrixSeries = list[np.ndarray]


@dataclass(frozen=True)
class ScalarChannelValuationRecord:
    """One matrix-inner channel-valuation audit."""

    construction_kind: str
    length: int
    predicted_leakage_coefficient: str
    smallest_step: str
    scaled_leakage_at_smallest_step: str
    scaled_leakage_error: str
    tangent_coefficient_error: str
    maximum_parseval_error: str
    exact_channel_error: str
    all_checks_passed: bool


def sphere_vector(theta: float, phase: float) -> np.ndarray:
    """Parametrize the complex two-sphere modulo global phase."""

    return np.array(
        [
            np.cos(theta / 2),
            np.exp(1j * phase) * np.sin(theta / 2),
        ],
        dtype=complex,
    )


def optimize_on_sphere(
    objective: Callable[[np.ndarray], float],
    *,
    maximize: bool,
) -> float:
    """Deterministically optimize a smooth real function of one qubit."""

    sign = -1.0 if maximize else 1.0
    best = -np.inf if maximize else np.inf
    for theta in np.linspace(0.0, np.pi, 7):
        for phase in np.linspace(0.0, 2 * np.pi, 13)[:-1]:
            result = minimize(
                lambda coordinates: sign
                * float(
                    objective(
                        sphere_vector(
                            coordinates[0],
                            coordinates[1],
                        )
                    )
                ),
                np.array([theta, phase]),
                method="Nelder-Mead",
                bounds=((0.0, np.pi), (0.0, 2 * np.pi)),
                options={
                    "maxiter": 800,
                    "xatol": 2e-13,
                    "fatol": 2e-15,
                },
            )
            value = sign * float(result.fun)
            best = max(best, value) if maximize else min(best, value)
    return float(best)


def channel_score(coefficients: MatrixSeries) -> float:
    """Compute the diagnostic product-state score for 2 x 2 data."""

    def fixed_input_score(vector: np.ndarray) -> float:
        images = np.asarray(
            [coefficient @ vector for coefficient in coefficients]
        )
        covariance = images.T @ images.conj()
        return float(np.linalg.eigvalsh(covariance)[-1])

    return optimize_on_sphere(fixed_input_score, maximize=True)


def leading_leakage(first_jets: MatrixSeries, length: int) -> float:
    """Compute formula (3) by deterministic sphere minimization."""

    off_monomial = [
        coefficient
        for degree, coefficient in enumerate(first_jets)
        if degree != length
    ]

    def leakage(vector: np.ndarray) -> float:
        projector = np.eye(2, dtype=complex) - np.outer(
            vector,
            vector.conj(),
        )
        return float(
            sum(
                np.linalg.norm(projector @ coefficient @ vector) ** 2
                for coefficient in off_monomial
            )
        )

    return optimize_on_sphere(leakage, maximize=False)


def schur_path(
    length: int,
    directions: dict[int, np.ndarray],
    step: float,
    count: int,
    terminal_generator: np.ndarray | None = None,
) -> MatrixSeries:
    """Reconstruct one exact matrix-inner Schur path."""

    zero = np.zeros((2, 2), dtype=complex)
    parameters = [zero.copy() for _ in range(length)]
    for index, direction in directions.items():
        parameters[index] = step * direction
    terminal = np.eye(2, dtype=complex)
    if terminal_generator is not None:
        terminal = expm(step * terminal_generator)
    return reconstruct_series(parameters, terminal, count)


def numerical_tangent(
    length: int,
    directions: dict[int, np.ndarray],
    count: int,
    terminal_generator: np.ndarray | None = None,
) -> MatrixSeries:
    """Return a centered finite-difference coefficient tangent."""

    step = 2e-7
    plus = schur_path(
        length,
        directions,
        step,
        count,
        terminal_generator,
    )
    minus = schur_path(
        length,
        directions,
        -step,
        count,
        terminal_generator,
    )
    return [
        (positive - negative) / (2 * step)
        for positive, negative in zip(plus, minus, strict=True)
    ]


def expected_tangent(
    length: int,
    directions: dict[int, np.ndarray],
    count: int,
    terminal_generator: np.ndarray | None = None,
) -> MatrixSeries:
    """Return the exact monomial Schur tangent."""

    tangent = [
        np.zeros((2, 2), dtype=complex) for _ in range(count)
    ]
    for index, direction in directions.items():
        tangent[index] += direction
        tangent[2 * length - index] -= direction.conj().T
    if terminal_generator is not None:
        tangent[length] += terminal_generator
    return tangent


def parseval_error(coefficients: MatrixSeries) -> float:
    """Return the truncated matrix-inner Parseval residual."""

    gram = sum(
        coefficient.conj().T @ coefficient
        for coefficient in coefficients
    )
    return float(np.linalg.norm(gram - np.eye(2)))


def make_positive_record(
    kind: str,
    directions: dict[int, np.ndarray],
    predicted: float | None,
) -> ScalarChannelValuationRecord:
    """Audit one path with a positive leading leakage coefficient."""

    length = 3
    count = 80
    steps = (0.16, 0.10, 0.06, 0.04, 0.025)
    exact_tangent = expected_tangent(length, directions, count)
    finite_tangent = numerical_tangent(length, directions, count)
    tangent_error = max(
        float(np.linalg.norm(actual - expected))
        for actual, expected in zip(
            finite_tangent,
            exact_tangent,
            strict=True,
        )
    )
    measured_prediction = leading_leakage(exact_tangent, length)
    if predicted is None:
        predicted = measured_prediction
    prediction_error = abs(measured_prediction - predicted)

    maximum_parseval = 0.0
    scaled_values: list[float] = []
    for step in steps:
        coefficients = schur_path(
            length,
            directions,
            step,
            count,
        )
        maximum_parseval = max(
            maximum_parseval,
            parseval_error(coefficients),
        )
        score = channel_score(coefficients)
        scaled_values.append((1 - score) / step**2)

    scaled_error = abs(scaled_values[-1] - predicted)
    verified = bool(
        predicted > 1e-5
        and prediction_error < 3e-10
        and tangent_error < 3e-8
        and maximum_parseval < 3e-11
        and scaled_error < 4e-4
        and scaled_error < abs(scaled_values[0] - predicted)
    )
    if not verified:
        raise RuntimeError(f"{kind} positive channel audit failed")
    return ScalarChannelValuationRecord(
        construction_kind=kind,
        length=length,
        predicted_leakage_coefficient=format_float(predicted),
        smallest_step=format_float(steps[-1]),
        scaled_leakage_at_smallest_step=format_float(
            scaled_values[-1]
        ),
        scaled_leakage_error=format_float(scaled_error),
        tangent_coefficient_error=format_float(tangent_error),
        maximum_parseval_error=format_float(maximum_parseval),
        exact_channel_error=format_float(0.0),
        all_checks_passed=verified,
    )


def make_exact_channel_record(
    kind: str,
    directions: dict[int, np.ndarray],
    terminal_generator: np.ndarray | None = None,
) -> ScalarChannelValuationRecord:
    """Audit a path which retains an exact scalar channel."""

    length = 3
    count = 80
    steps = (0.16, 0.08, 0.04)
    exact_tangent = expected_tangent(
        length,
        directions,
        count,
        terminal_generator,
    )
    finite_tangent = numerical_tangent(
        length,
        directions,
        count,
        terminal_generator,
    )
    tangent_error = max(
        float(np.linalg.norm(actual - expected))
        for actual, expected in zip(
            finite_tangent,
            exact_tangent,
            strict=True,
        )
    )
    predicted = leading_leakage(exact_tangent, length)

    channel_error = 0.0
    maximum_parseval = 0.0
    for step in steps:
        coefficients = schur_path(
            length,
            directions,
            step,
            count,
            terminal_generator,
        )
        maximum_parseval = max(
            maximum_parseval,
            parseval_error(coefficients),
        )
        channel_error = max(
            channel_error,
            abs(1 - channel_score(coefficients)),
        )

    verified = bool(
        abs(predicted) < 3e-11
        and tangent_error < 3e-8
        and maximum_parseval < 3e-11
        and channel_error < 3e-11
    )
    if not verified:
        raise RuntimeError(f"{kind} exact-channel audit failed")
    return ScalarChannelValuationRecord(
        construction_kind=kind,
        length=length,
        predicted_leakage_coefficient=format_float(predicted),
        smallest_step=format_float(steps[-1]),
        scaled_leakage_at_smallest_step=format_float(0.0),
        scaled_leakage_error=format_float(0.0),
        tangent_coefficient_error=format_float(tangent_error),
        maximum_parseval_error=format_float(maximum_parseval),
        exact_channel_error=format_float(channel_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[ScalarChannelValuationRecord]:
    """Return positive and exactly split deterministic paths."""

    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    pauli_z = np.diag([1, -1]).astype(complex)
    amplitude = 0.25
    generic_one = np.array(
        [[0.20, 0.10j], [-0.05, 0.03]],
        dtype=complex,
    )
    generic_two = np.array(
        [[0.02 + 0.04j, -0.11], [0.07j, -0.15]],
        dtype=complex,
    )
    diagonal_one = np.diag([0.20, -0.12]).astype(complex)
    diagonal_two = np.diag([0.08j, -0.17j]).astype(complex)
    rotation = np.array([[0, 0.21], [-0.21, 0]], dtype=complex)

    return [
        make_positive_record(
            "pauli_transverse",
            {
                1: amplitude * pauli_x,
                2: amplitude * pauli_z,
            },
            2 * amplitude**2,
        ),
        make_positive_record(
            "generic_nonnormal_transverse",
            {1: generic_one, 2: generic_two},
            None,
        ),
        make_exact_channel_record(
            "commuting_diagonal_channel",
            {1: diagonal_one, 2: diagonal_two},
        ),
        make_exact_channel_record(
            "terminal_unitary_rotation",
            {},
            rotation,
        ),
    ]


def write_records(
    records: list[ScalarChannelValuationRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_scalar_channel_valuation_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the scalar-channel valuation audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
