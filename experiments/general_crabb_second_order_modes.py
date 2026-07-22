#!/usr/bin/env python3
"""Regenerate the arbitrary-size single-Crabb second-variation formula.

The script evaluates the finite L62 conformal/metric reduction independently
of the closed mode formula in L65.  It checks one deterministic complex vector
in every nonzero rotational mode and, at the smallest sizes, reconstructs the
entire real quadratic form by polarization.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crouzeix import crabb_matrix
from general_similarity_second_order_probe import physical_conformal_coefficients
from general_similarity_tangent_probe import crabb_metric


@dataclass(frozen=True)
class ValidationRecord:
    dimension: int
    seed: int
    support_resolution: int
    maximum_mode_residual: float
    maximum_relative_mode_residual: float
    minimum_mode_value: float
    grade_zero_value: float
    full_form_rank: int | None
    expected_rank: int
    full_form_maximum_eigenvalue: float | None


def limiting_inverse(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Return ``lim_c->infinity (matrix+c vv*)^-1`` stably."""

    regularized = matrix + np.outer(vector, vector)
    inverse = np.linalg.inv(regularized)
    image = inverse @ vector
    return inverse - np.outer(image, image.conj()) / np.vdot(vector, image)


def mode_kernel(length: int, mode: int) -> np.ndarray:
    """Return the positive kernel ``K_{length,mode}`` from L65."""

    if length == 1:
        if mode < 3:
            return np.zeros((1, 1))
        scalar = mode * (mode - 1) * (mode - 2) / (6 * (mode + 1) ** 2)
        return np.array([[scalar]])

    crabb = crabb_matrix(length - 1)
    support = np.real((crabb + crabb.conj().T) / 2)
    base = (np.eye(length) - support) / 4
    vector = np.concatenate(
        (
            [(mode + 2) / np.sqrt(2)],
            np.ones(length - 2),
            [(mode + 2) / np.sqrt(2)],
        )
    )
    if mode < 3:
        return limiting_inverse(base, vector)
    coefficient = 3 / (2 * mode * (mode - 1) * (mode - 2))
    return np.linalg.inv(base + coefficient * np.outer(vector, vector))


def mode_aggregator(length: int, mode: int) -> np.ndarray:
    """Map the aligned ``(+mode,-mode)`` entries to the reduced vector."""

    if length == 1:
        return np.array([[1, np.sqrt(2), 2, np.sqrt(2)]])

    aggregator = np.zeros((length, 2 * length + 2))
    aggregator[:, :length] = np.eye(length)
    edge = 1 / np.sqrt(2) if mode == 1 else 1
    aggregator[0, length] = edge
    aggregator[0, length + 1] = np.sqrt(2)
    for index in range(1, length - 1):
        aggregator[index, length + index + 1] = 1
    aggregator[-1, 2 * length] = np.sqrt(2)
    aggregator[-1, 2 * length + 1] = edge
    return aggregator


def reduced_second_order_value(
    perturbation: np.ndarray,
    support_resolution: int,
) -> float:
    """Eliminate the general L62 single-block metric program in scalars."""

    dimension = perturbation.shape[0]
    last = dimension - 1
    base = crabb_matrix(last)
    metric = crabb_metric(dimension, 1, "operator_weight")
    weights = np.real(np.diag(base, 1))
    tangent, second_order = physical_conformal_coefficients(
        base,
        perturbation,
        support_resolution,
    )
    first_forcing = tangent.conj().T @ metric @ base + base.conj().T @ metric @ tangent
    contraction_weights = np.concatenate((2 * np.ones(last - 1), [1]))

    def objective(free: np.ndarray) -> float:
        metric_tangent = np.zeros((dimension, dimension), dtype=complex)
        metric_tangent[0, 1:] = free
        metric_tangent[1:, 0] = free.conj()
        for row in range(1, dimension):
            for column in range(1, dimension):
                metric_tangent[row, column] = (
                    weights[row - 1]
                    * weights[column - 1]
                    * metric_tangent[row - 1, column - 1]
                    + first_forcing[row, column]
                )
        if abs(metric_tangent[last, last]) > 2e-9:
            raise RuntimeError("terminal first-order consistency failed")

        second_forcing = (
            second_order.conj().T @ metric @ base
            + base.conj().T @ metric @ second_order
            + tangent.conj().T @ metric @ tangent
            + tangent.conj().T @ metric_tangent @ base
            + base.conj().T @ metric_tangent @ tangent
        )
        first_variation = (
            metric_tangent - base.conj().T @ metric_tangent @ base - first_forcing
        )
        range_coupling = first_variation[1:, 0]
        lower = np.sum(abs(free[:-1]) ** 2) + abs(free[-1]) ** 2 / 3
        upper = (
            abs(metric_tangent[last, 0]) ** 2 / 3
            + np.sum(abs(metric_tangent[last, 1:last]) ** 2) / 2
        )
        contraction = np.sum(
            contraction_weights
            * (np.real(np.diag(second_forcing)[1:]) + abs(range_coupling) ** 2)
        )
        return float(np.real(4 * lower + upper + contraction))

    zero = np.zeros(last, dtype=complex)
    optimum = objective(zero)
    for index in range(last):
        real_basis = np.zeros(last, dtype=complex)
        imaginary_basis = np.zeros(last, dtype=complex)
        real_basis[index] = 1
        imaginary_basis[index] = 1j
        real_linear = (objective(real_basis) - objective(-real_basis)) / 2
        imaginary_linear = (
            objective(imaginary_basis) - objective(-imaginary_basis)
        ) / 2
        hessian = 8 if index < last - 1 else 8 / 3
        optimum -= (real_linear**2 + imaginary_linear**2) / (4 * hessian)
    return optimum


def paired_mode_case(
    dimension: int,
    mode: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float]:
    """Return a random perturbation in one paired mode and its closed value."""

    last = dimension - 1
    length = last - mode
    perturbation = np.zeros((dimension, dimension), dtype=complex)
    for index in range(length):
        perturbation[index, index + mode + 1] = (
            rng.standard_normal() + 1j * rng.standard_normal()
        )
    for index in range(length + 2):
        perturbation[index + mode - 1, index] = (
            rng.standard_normal() + 1j * rng.standard_normal()
        )
    positive = np.conj(
        [perturbation[index, index + mode + 1] for index in range(length)]
    )
    negative = np.array(
        [perturbation[index + mode - 1, index] for index in range(length + 2)]
    )
    reduced = mode_aggregator(length, mode) @ np.concatenate((positive, negative))
    kernel = mode_kernel(length, mode)
    predicted = -float(np.real(np.vdot(reduced, kernel @ reduced)))
    return perturbation, predicted


def bottom_mode_cases(
    dimension: int,
    rng: np.random.Generator,
) -> list[tuple[np.ndarray, float]]:
    """Return the two unpaired bottom modes (L65 equations (19)--(20))."""

    last = dimension - 1
    penultimate = np.zeros((dimension, dimension), dtype=complex)
    values = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    for index, value in enumerate(values):
        penultimate[index + last - 1, index] = value
    coefficient = (last - 1) * (last - 2) / (3 * last)
    penultimate_value = -coefficient * abs(np.sum(values)) ** 2

    bottom = np.zeros((dimension, dimension), dtype=complex)
    bottom[last, 0] = rng.standard_normal() + 1j * rng.standard_normal()
    coefficient = (last**2 + 36 * last - 13) / (6 * last)
    bottom_value = -coefficient * abs(bottom[last, 0]) ** 2
    return [(penultimate, penultimate_value), (bottom, bottom_value)]


def full_quadratic_form(
    dimension: int,
    support_resolution: int,
) -> np.ndarray:
    """Reconstruct the real quadratic form by exact polarization in floats."""

    basis = []
    for imaginary in (False, True):
        for index in range(dimension**2):
            matrix = np.zeros((dimension, dimension), dtype=complex)
            matrix.flat[index] = 1j if imaginary else 1
            basis.append(matrix)
    size = len(basis)
    form = np.zeros((size, size))
    for row, matrix in enumerate(basis):
        form[row, row] = reduced_second_order_value(matrix, support_resolution)
    for row in range(size):
        for column in range(row):
            mixed = reduced_second_order_value(
                basis[row] + basis[column],
                support_resolution,
            )
            entry = (mixed - form[row, row] - form[column, column]) / 2
            form[row, column] = form[column, row] = entry
    return form


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=30)
    parser.add_argument("--full-form-maximum-size", type=int, default=8)
    parser.add_argument("--seed", type=int, default=70221)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    rng = np.random.default_rng(args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for dimension in range(args.minimum_size, args.maximum_size + 1):
            support_resolution = max(64, 4 * dimension)
            residuals = []
            relative_residuals = []
            predicted_values = []
            for mode in range(1, dimension - 1):
                perturbation, predicted = paired_mode_case(dimension, mode, rng)
                actual = reduced_second_order_value(
                    perturbation,
                    support_resolution,
                )
                residuals.append(abs(actual - predicted))
                relative_residuals.append(
                    abs(actual - predicted) / max(1, abs(actual), abs(predicted))
                )
                predicted_values.append(predicted)
            for perturbation, predicted in bottom_mode_cases(dimension, rng):
                actual = reduced_second_order_value(
                    perturbation,
                    support_resolution,
                )
                residuals.append(abs(actual - predicted))
                relative_residuals.append(
                    abs(actual - predicted) / max(1, abs(actual), abs(predicted))
                )
                predicted_values.append(predicted)

            grade_zero = np.zeros((dimension, dimension), dtype=complex)
            grade_zero[np.arange(dimension - 1), np.arange(1, dimension)] = (
                rng.standard_normal(dimension - 1)
                + 1j * rng.standard_normal(dimension - 1)
            )
            grade_zero_value = reduced_second_order_value(
                grade_zero,
                support_resolution,
            )

            rank = None
            maximum_eigenvalue = None
            if dimension <= args.full_form_maximum_size:
                form = full_quadratic_form(dimension, support_resolution)
                eigenvalues = np.linalg.eigvalsh(form)
                rank = int(np.count_nonzero(eigenvalues < -1e-8))
                maximum_eigenvalue = float(eigenvalues[-1])
                if rank != dimension * (dimension - 2):
                    raise RuntimeError("unexpected full-form rank")
                if maximum_eigenvalue > 1e-8:
                    raise RuntimeError("full form is not negative semidefinite")
            if max(relative_residuals) > 1e-10:
                raise RuntimeError("closed mode formula failed its residual gate")
            if max(predicted_values) > 1e-10 or grade_zero_value > 1e-8:
                raise RuntimeError("a predicted second-order value is positive")
            record = ValidationRecord(
                dimension=dimension,
                seed=args.seed,
                support_resolution=support_resolution,
                maximum_mode_residual=max(residuals),
                maximum_relative_mode_residual=max(relative_residuals),
                minimum_mode_value=min(predicted_values),
                grade_zero_value=grade_zero_value,
                full_form_rank=rank,
                expected_rank=dimension * (dimension - 2),
                full_form_maximum_eigenvalue=maximum_eigenvalue,
            )
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
