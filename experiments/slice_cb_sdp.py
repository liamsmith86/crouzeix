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

import cvxpy as cp
import numpy as np

from slice_phase_audit import NodalModel


DEFAULT_CASES = (
    (2.0, 1.2, 1.6, 0.08),
    (1.5, 2.0, 0.9, 0.12),
    (0.8, 2.4, 1.1, 0.3),
    (2.5, 0.4, 2.0, 0.4),
)


def conformal_matrix(weights: tuple[float, float, float, float]) -> np.ndarray:
    model = NodalModel.build(weights)
    return np.real_if_close(model.calculus(model.disk_nodes)).real


def similarity_sdp(weights: tuple[float, float, float, float]) -> tuple[float, float, float]:
    operator = conformal_matrix(weights)
    size = operator.shape[0]
    identity = np.eye(size)
    metric = cp.Variable((size, size), symmetric=True)
    bound = cp.Variable()
    problem = cp.Problem(
        cp.Minimize(bound),
        [
            metric - identity >> 0,
            bound * identity - metric >> 0,
            metric - operator.T @ metric @ operator >> 0,
        ],
    )
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-9,
        tol_gap_rel=1e-9,
        tol_feas=1e-9,
        max_iter=500,
    )
    if problem.status != cp.OPTIMAL:
        raise RuntimeError(f"SDP failed for {weights}: {problem.status}")
    contraction_slack = np.min(
        np.linalg.eigvalsh(metric.value - operator.T @ metric.value @ operator)
    )
    return float(bound.value), float(np.sqrt(bound.value)), float(contraction_slack)


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
        bound, condition, slack = similarity_sdp(weights)
        largest_bound = max(largest_bound, bound)
        print(
            f"weights={weights} t={bound:.9f} sqrt(t)={condition:.9f} "
            f"contraction_slack={slack:+.2e}"
        )
    print(f"largest t: {largest_bound:.9f}")


if __name__ == "__main__":
    main()
