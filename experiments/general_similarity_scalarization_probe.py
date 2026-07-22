#!/usr/bin/env python3
"""Test state scalarizations of the Crouzeix--Palencia correction.

After pulling a numerical-range domain back to the disk, write

    theta(f) = f(T),
    R(f) = theta(alpha(f))*.

Crouzeix--Palencia controls ``theta + R``, whereas the Hartz--McCarthy
inequality accepts a scalar correction ``theta + beta I``.  A natural family
of attempted scalarizations is

    beta_Q(f) = tr(Q R(f)),  Q >= 0, tr(Q) = 1.

If ``(theta + beta_Q I) / 2`` were completely contractive, its operator-valued
Toeplitz moment matrices would be positive semidefinite at every order.  This
script maximizes the least eigenvalue of one such matrix over every density
matrix Q.  A negative optimum numerically excludes the entire state family at
that order.  The computation is a falsification probe, not an interval proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from general_similarity_sdp import generated_cases, parse_case_selector
from theodorsen import GeneralPullback, theodorsen_map


@dataclass(frozen=True)
class ScalarizationRecord:
    dimension: int
    family: str
    sample: int
    per_family: int
    seed: int
    matrix_sha256: str
    inflate: float
    resolution: int
    order: int
    solver: str
    status: str
    objective: float
    verified_minimum_eigenvalue: float
    density_minimum_eigenvalue: float
    density_trace_error: float
    theodorsen_error: float
    scalar_unitality_error: float
    matrix_unitality_error: float
    dlp_minimum: float
    dlp_mass_error: float


def select_case(
    selector: tuple[int, str, int],
    seed: int,
    per_family: int,
) -> np.ndarray:
    """Regenerate and return one deterministic general-matrix test case."""

    dimension, family, sample = selector
    if sample >= per_family:
        raise ValueError(f"sample {sample} requires --per-family at least {sample + 1}")
    for case_dimension, case_family, case_sample, matrix in generated_cases(
        3,
        dimension,
        per_family,
        seed,
    ):
        if (case_dimension, case_family, case_sample) == selector:
            return matrix
    raise ValueError(f"case was not generated: {dimension}:{family}:{sample}")


def matrix_sha256(matrix: np.ndarray) -> str:
    """Return a stable fingerprint of a contiguous complex128 matrix."""

    canonical = np.ascontiguousarray(matrix, dtype=np.complex128)
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def cp_moments(
    matrix: np.ndarray,
    resolution: int,
    inflate: float,
    order: int,
) -> tuple[list[np.ndarray], list[np.ndarray], dict[str, float]]:
    """Compute theta(z^k), R(z^k), and conformal-map diagnostics."""

    boundary, derivative, map_error = theodorsen_map(
        matrix,
        N=resolution,
        inflate=inflate,
        iters=500,
        tol=1e-13,
    )
    pullback = GeneralPullback(boundary, derivative)
    resolvents = pullback.resolvent_stack(matrix)
    identity = np.eye(matrix.shape[0], dtype=complex)

    theta = [identity]
    correction = [identity]
    for exponent in range(1, order + 1):
        values = pullback.w**exponent
        theta.append(pullback.calc(values, matrix, resolvents))
        correction.append(pullback.calc(np.conj(values), matrix, resolvents).conj().T)

    calculated_identity = pullback.calc(np.ones(resolution), matrix, resolvents)
    dlp_minimum, dlp_mass_error = pullback.dlp_certificate(matrix)
    diagnostics = {
        "theodorsen_error": float(map_error),
        "scalar_unitality_error": float(pullback.unitality_certificate()),
        "matrix_unitality_error": float(
            np.linalg.norm(calculated_identity - identity, 2)
        ),
        "dlp_minimum": float(dlp_minimum),
        "dlp_mass_error": float(dlp_mass_error),
    }
    return theta, correction, diagnostics


def toeplitz_blocks(
    theta: list[np.ndarray],
    correction: list[np.ndarray],
    density: cp.Variable,
) -> list[list[cp.Expression | np.ndarray]]:
    """Build the unital Toeplitz moment matrix for a variable state."""

    order = len(theta) - 1
    dimension = theta[0].shape[0]
    identity = np.eye(dimension)
    rows: list[list[cp.Expression | np.ndarray]] = []
    for row_index in range(order + 1):
        row: list[cp.Expression | np.ndarray] = []
        for column_index in range(order + 1):
            exponent = column_index - row_index
            if exponent == 0:
                block = identity
            elif exponent > 0:
                beta = cp.trace(density @ correction[exponent])
                block = (theta[exponent] + beta * identity) / 2
            else:
                positive_exponent = -exponent
                beta = cp.trace(density @ correction[positive_exponent])
                block = (
                    theta[positive_exponent].conj().T + cp.conj(beta) * identity
                ) / 2
            row.append(block)
        rows.append(row)
    return rows


def numerical_toeplitz(
    theta: list[np.ndarray],
    correction: list[np.ndarray],
    density: np.ndarray,
) -> np.ndarray:
    """Independently rebuild a numeric Toeplitz matrix from an SDP state."""

    order = len(theta) - 1
    dimension = theta[0].shape[0]
    identity = np.eye(dimension)
    matrix = np.zeros(((order + 1) * dimension,) * 2, dtype=complex)
    beta = [complex(np.trace(density @ value)) for value in correction]
    for row_index in range(order + 1):
        for column_index in range(order + 1):
            exponent = column_index - row_index
            if exponent == 0:
                block = identity
            elif exponent > 0:
                block = (theta[exponent] + beta[exponent] * identity) / 2
            else:
                positive_exponent = -exponent
                block = (
                    theta[positive_exponent].conj().T
                    + np.conj(beta[positive_exponent]) * identity
                ) / 2
            row_slice = slice(row_index * dimension, (row_index + 1) * dimension)
            column_slice = slice(
                column_index * dimension, (column_index + 1) * dimension
            )
            matrix[row_slice, column_slice] = block
    return matrix


def solve_state_probe(
    theta: list[np.ndarray],
    correction: list[np.ndarray],
    solver: str,
) -> tuple[str, float, np.ndarray, float]:
    """Maximize the least Toeplitz eigenvalue over all density matrices."""

    dimension = theta[0].shape[0]
    toeplitz_dimension = len(theta) * dimension
    density = cp.Variable((dimension, dimension), hermitian=True)
    lower_bound = cp.Variable()
    toeplitz = cp.bmat(toeplitz_blocks(theta, correction, density))
    problem = cp.Problem(
        cp.Maximize(lower_bound),
        [
            density >> 0,
            cp.trace(density) == 1,
            toeplitz - lower_bound * np.eye(toeplitz_dimension) >> 0,
        ],
    )
    options: dict[str, float | int | bool] = {"verbose": False}
    if solver == "CLARABEL":
        options.update(tol_gap_abs=1e-9, tol_feas=1e-9, tol_gap_rel=1e-9, max_iter=2000)
    elif solver == "SCS":
        options.update(eps=1e-8, max_iters=200000)
    problem.solve(solver=solver, **options)
    if density.value is None or problem.value is None:
        raise RuntimeError(
            f"{solver} did not return a primal solution: {problem.status}"
        )

    state = np.asarray(density.value)
    rebuilt = numerical_toeplitz(theta, correction, state)
    verified_minimum = float(np.linalg.eigvalsh(rebuilt)[0])
    return str(problem.status), float(problem.value), state, verified_minimum


def evaluate(
    selector: tuple[int, str, int],
    seed: int,
    per_family: int,
    matrix: np.ndarray,
    resolution: int,
    inflate: float,
    order: int,
    solver: str,
) -> ScalarizationRecord:
    """Evaluate one resolution and package all audit data."""

    theta, correction, diagnostics = cp_moments(matrix, resolution, inflate, order)
    status, objective, density, verified_minimum = solve_state_probe(
        theta,
        correction,
        solver,
    )
    dimension, family, sample = selector
    return ScalarizationRecord(
        dimension=dimension,
        family=family,
        sample=sample,
        per_family=per_family,
        seed=seed,
        matrix_sha256=matrix_sha256(matrix),
        inflate=inflate,
        resolution=resolution,
        order=order,
        solver=solver,
        status=status,
        objective=objective,
        verified_minimum_eigenvalue=verified_minimum,
        density_minimum_eigenvalue=float(np.linalg.eigvalsh(density)[0]),
        density_trace_error=float(abs(np.trace(density) - 1)),
        **diagnostics,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=parse_case_selector, default=(3, "dense", 0))
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--per-family", type=int, default=5)
    parser.add_argument("--inflate", type=float, default=0.01)
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument(
        "--resolutions", nargs="+", type=int, default=[512, 1024, 2048, 4096]
    )
    parser.add_argument("--solver", choices=("CLARABEL", "SCS"), default="CLARABEL")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.order < 1:
        raise ValueError("order must be positive")
    if any(resolution < 32 for resolution in args.resolutions):
        raise ValueError("every resolution must be at least 32")

    matrix = select_case(args.case, args.seed, args.per_family)
    records = []
    for resolution in args.resolutions:
        record = evaluate(
            args.case,
            args.seed,
            args.per_family,
            matrix,
            resolution,
            args.inflate,
            args.order,
            args.solver,
        )
        records.append(asdict(record))
        print(json.dumps(records[-1], sort_keys=True), flush=True)

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            "".join(json.dumps(record, sort_keys=True) + "\n" for record in records)
        )


if __name__ == "__main__":
    main()
