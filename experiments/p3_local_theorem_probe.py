#!/usr/bin/env python3
"""Numerically stress the L73 rank-one certificate near the p=3 disk curve."""

from __future__ import annotations

import argparse

import numpy as np
from scipy.linalg import solve_discrete_lyapunov
from scipy.optimize import minimize

from general_similarity_sdp import evaluate_map


ROOT_TWO = np.sqrt(2.0)


def disk_scale(radius: float) -> float:
    """Return the analytic root of L71's disk relation near one."""

    squared = radius * radius
    roots = np.roots([81 * squared**2, 1152 * squared - 4096, 4096])
    return float(min(roots, key=lambda value: abs(value - 1)).real)


def slice_matrix(
    z: complex,
    w: complex,
    s: float,
    v: complex,
) -> np.ndarray:
    """Construct the exact L69 slice matrix."""

    return np.array(
        [
            [-z / 3, ROOT_TWO + s, np.conjugate(z)],
            [w, 2 * z / 3, ROOT_TWO - s],
            [v, w, -z / 3],
        ],
        dtype=complex,
    )


def rank_one_condition(operator: np.ndarray, variables: np.ndarray) -> float:
    """Evaluate the condition of the rank-one Stein Gramian."""

    defect = np.array(
        [
            1,
            variables[0] + 1j * variables[1],
            variables[2] + 1j * variables[3],
        ]
    )[:, None]
    metric = solve_discrete_lyapunov(
        operator.conjugate().T,
        defect @ defect.conjugate().T,
    )
    eigenvalues = np.linalg.eigvalsh(metric)
    return float(eigenvalues[-1] / eigenvalues[0])


def complex_normal(generator: np.random.Generator) -> complex:
    """Return one normalized complex Gaussian draw."""

    return complex(*generator.normal(size=2)) / ROOT_TWO


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=18)
    parser.add_argument("--seed", type=int, default=73022)
    parser.add_argument("--resolution", type=int, default=2048)
    args = parser.parse_args()

    generator = np.random.default_rng(args.seed)
    maximum_condition = -np.inf
    maximum_map_error = 0.0
    violations = 0
    for sample in range(args.samples):
        radius = float(np.exp(generator.uniform(np.log(0.03), np.log(0.2))))
        phase = np.exp(1j * generator.uniform(0.0, 2 * np.pi))
        z = radius * phase
        scale = disk_scale(radius)
        disk_w = 3 * ROOT_TWO * z**2 * scale / 64
        disk_v = -9 * z**3 * scale / 64

        family = sample % 3
        if family == 0:
            soft = 0.2 * radius**2 * complex_normal(generator)
            strong = 0.2 * radius**3 * complex_normal(generator)
            transverse = 0.2 * radius**3 * generator.normal()
        elif family == 1:
            soft = 0.01 * radius**2 * complex_normal(generator)
            strong = (
                -4 * ROOT_TWO * z * soft / 7
                + 0.002 * radius**3 * complex_normal(generator)
            )
            transverse = 0.002 * radius**3 * generator.normal()
        else:
            soft = 0.04 * radius * complex_normal(generator)
            strong = 0.04 * radius * complex_normal(generator)
            transverse = 0.04 * radius * generator.normal()

        matrix = slice_matrix(
            z,
            disk_w + soft,
            transverse,
            disk_v + strong,
        )
        evaluation = evaluate_map(matrix, args.resolution, 0.0)
        optimum = minimize(
            lambda variables: rank_one_condition(
                evaluation.operator,
                variables,
            ),
            np.zeros(4),
            method="BFGS",
            options={"gtol": 1e-10, "maxiter": 1000},
        )
        condition = float(optimum.fun)
        maximum_condition = max(maximum_condition, condition)
        maximum_map_error = max(
            maximum_map_error,
            evaluation.theodorsen_error,
            evaluation.matrix_unitality_error,
        )
        if condition > 4 + 5e-8:
            violations += 1

    if violations:
        raise AssertionError(f"{violations} numerical certificate violations")
    print("PASS p=3 local-theorem numerical probe")
    print(f"samples = {args.samples}, seed = {args.seed}")
    print(f"maximum rank-one condition = {maximum_condition:.12f}")
    print(f"maximum map diagnostic = {maximum_map_error:.3e}")


if __name__ == "__main__":
    main()
