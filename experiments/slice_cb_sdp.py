"""Search for a complete-2 similarity certificate on the elliptic 4x4 slice.

For T = phi(A), solve the convex problem

    minimize t  subject to  I <= P <= t I  and  T* P T <= P.

If t <= 4, then P^(1/2) T P^(-1/2) is a contraction and its similarity has condition number
at most 2.  Von Neumann's inequality then proves the complete Crouzeix bound for that matrix.

The SDP output is exploratory numerical evidence.  A rigorous uniform construction of P with
t <= 4 is the proof target; solver output alone is not a certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import cvxpy as cp
import numpy as np

from slice_phase_audit import NodalModel


DEFAULT_CASES = (
    (2.0, 1.2, 1.6, 0.08),
    (1.5, 2.0, 0.9, 0.12),
    (0.8, 2.4, 1.1, 0.3),
    (2.5, 0.4, 2.0, 0.4),
)


@dataclass(frozen=True)
class SimilarityCertificate:
    """Numerical primal/dual output for the similarity SDP."""

    bound: float
    metric: np.ndarray
    dual_witness: np.ndarray
    contraction_slack: float
    dual_ratio: float
    solver: str


def conformal_matrix(weights: tuple[float, float, float, float]) -> np.ndarray:
    model = NodalModel.build(weights)
    return np.real_if_close(model.calculus(model.disk_nodes)).real


def trace_ratio(operator: np.ndarray, witness: np.ndarray) -> float:
    """Return tr(D_-)/tr(D_+) for D = Z - T Z T*.

    The ratio is the exact dual quantity from ``proof/slice_similarity_duality.md``.  Inputs here
    are floating-point solver output, so the result is a regression diagnostic, not a certificate.
    """

    adjoint = operator.conj().T
    difference = witness - operator @ witness @ adjoint
    eigenvalues = np.linalg.eigvalsh((difference + difference.conj().T) / 2)
    positive_trace = float(np.maximum(eigenvalues, 0).sum())
    negative_trace = float(np.maximum(-eigenvalues, 0).sum())
    if positive_trace == 0:
        return np.inf if negative_trace > 0 else 0.0
    return negative_trace / positive_trace


def solve_similarity_sdp(operator: np.ndarray) -> SimilarityCertificate:
    """Solve the SDP, with a robust fallback for nearly degenerate slice corners."""

    size = operator.shape[0]
    identity = np.eye(size)
    is_complex = bool(np.iscomplexobj(operator) and np.max(abs(operator.imag)) > 1e-14)
    metric = (
        cp.Variable((size, size), hermitian=True)
        if is_complex
        else cp.Variable((size, size), symmetric=True)
    )
    bound = cp.Variable()
    adjoint = operator.conj().T
    constraints = [
        metric - identity >> 0,
        bound * identity - metric >> 0,
        metric - adjoint @ metric @ operator >> 0,
    ]
    problem = cp.Problem(
        cp.Minimize(bound),
        constraints,
    )

    solver = "CLARABEL"
    try:
        problem.solve(
            solver=solver,
            tol_gap_abs=1e-9,
            tol_gap_rel=1e-9,
            tol_feas=1e-9,
            max_iter=500,
        )
    except cp.error.SolverError:
        pass
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        solver = "SCS"
        problem.solve(solver=solver, eps=2e-7, max_iters=100_000, verbose=False)
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"similarity SDP failed: {problem.status}")

    metric_value = np.asarray(metric.value)
    dual_witness = np.asarray(constraints[2].dual_value)
    contraction_defect = metric_value - adjoint @ metric_value @ operator
    contraction_slack = np.min(
        np.linalg.eigvalsh((contraction_defect + contraction_defect.conj().T) / 2)
    )
    return SimilarityCertificate(
        bound=float(bound.value),
        metric=metric_value,
        dual_witness=dual_witness,
        contraction_slack=float(contraction_slack),
        dual_ratio=trace_ratio(operator, dual_witness),
        solver=solver,
    )


def similarity_sdp(weights: tuple[float, float, float, float]) -> tuple[float, float, float]:
    """Compatibility wrapper used by the original command-line report."""

    result = solve_similarity_sdp(conformal_matrix(weights))
    return result.bound, float(np.sqrt(result.bound)), result.contraction_slack


def random_cases(count: int, seed: int) -> list[tuple[float, float, float, float]]:
    generator = np.random.default_rng(seed)
    cases = []
    for _ in range(count):
        weights = np.exp(generator.uniform(np.log(0.2), np.log(3.0), 3))
        c = float(np.exp(generator.uniform(np.log(0.015), np.log(0.85))))
        cases.append((*map(float, weights), c))
    return cases


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--random", type=int, default=0, help="number of additional random cases")
    parser.add_argument("--seed", type=int, default=123)
    arguments = parser.parse_args()

    cases = list(DEFAULT_CASES) + random_cases(arguments.random, arguments.seed)
    largest_bound = 0.0
    for weights in cases:
        result = solve_similarity_sdp(conformal_matrix(weights))
        condition = float(np.sqrt(result.bound))
        largest_bound = max(largest_bound, result.bound)
        print(
            f"weights={weights} t={result.bound:.9f} sqrt(t)={condition:.9f} "
            f"dual_ratio={result.dual_ratio:.9f} "
            f"contraction_slack={result.contraction_slack:+.2e} solver={result.solver}"
        )
    print(f"largest t: {largest_bound:.9f}")


if __name__ == "__main__":
    main()
