#!/usr/bin/env python3
"""Probe the fixed-weight elliptic Crabb axis in arbitrary dimension.

For the ``p x p`` Crabb block ``C_p`` and ``0 < c < 1``, set

    A_c = C_p + c C_p^*.

Its numerical range is the ellipse with semiaxes ``1+c`` and ``1-c``.
This script constructs the normalized Riemann pullback ``T_c = phi_c(A_c)``
from its exact elliptic functional calculus and compares:

* the unrestricted L21 similarity SDP;
* the same SDP restricted to diagonal metrics;
* the Chebyshev--Blaschke lower bound ``k(c^(2p-2))/c^(p-1)``.

The SDP output is numerical evidence, not a proof of the all-size upper
bound.  The Chebyshev--Blaschke lower bound is an exact theorem once the
standard extremal finite-Blaschke result on an interval is imported.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from math import pi, sqrt
from pathlib import Path

import cvxpy as cp
import numpy as np
from scipy.special import ellipj, ellipk

from crouzeix import crabb_matrix
from slice_cb_sdp import solve_similarity_sdp


DEFAULT_ELLIPSE_PARAMETERS = (0.05, 0.15, 0.4, 0.6)


@dataclass(frozen=True)
class AxisRecord:
    dimension: int
    ellipse_parameter: float
    predicted_bound: float
    unrestricted_bound: float
    diagonal_bound: float
    unrestricted_dual_ratio: float
    unrestricted_prediction_error: float
    diagonal_prediction_error: float
    unrestricted_metric_offdiagonal: float
    diagonal_metric_reflection_error: float
    diagonal_stein_rank_one_error: float
    alternation_reversal_error: float
    explicit_low_size_certificate_error: float | None
    unrestricted_solver: str
    diagonal_solver: str


@dataclass(frozen=True)
class DiagonalCertificate:
    bound: float
    diagonal: np.ndarray
    solver: str


def elliptic_modulus_from_nome(nome: float) -> float:
    """Return Jacobi's modulus ``k(q)`` from its nome using the product.

    The product is preferable to numerical inversion here: the predicted
    all-size value involves nomes as small as ``c**(2*(p-1))``.
    """

    if not 0 < nome < 1:
        raise ValueError("the nome must lie strictly between zero and one")
    product = 1.0
    power = nome
    for _ in range(1, 100_000):
        numerator = 1.0 + power * nome
        denominator = 1.0 + power
        product *= (numerator / denominator) ** 4
        power *= nome**2
        if power < np.finfo(float).eps:
            break
    else:
        raise RuntimeError("Jacobi modulus product did not converge")
    return 4.0 * sqrt(nome) * product


def conformal_crabb_operator(
    dimension: int, ellipse_parameter: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct ``phi(C_p+c C_p*)`` by exact nodal functional calculus."""

    if dimension < 3:
        raise ValueError("the Crabb dimension must be at least three")
    if not 0 < ellipse_parameter < 1:
        raise ValueError("the ellipse parameter must lie in (0,1)")

    crabb = np.asarray(crabb_matrix(dimension - 1).real)
    symmetric_model = sqrt(ellipse_parameter) * (crabb + crabb.T)
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric_model)

    modulus = elliptic_modulus_from_nome(ellipse_parameter**2)
    parameter = modulus**2
    quarter_period = float(ellipk(parameter))
    lobatto_nodes = np.clip(
        eigenvalues / (2.0 * sqrt(ellipse_parameter)), -1.0, 1.0
    )
    arguments = (2.0 * quarter_period / pi) * np.arcsin(lobatto_nodes)
    disk_nodes = sqrt(modulus) * ellipj(arguments, parameter)[0]
    symmetric_pullback = (eigenvectors * disk_nodes) @ eigenvectors.T

    diagonal = ellipse_parameter ** (np.arange(dimension) / 2.0)
    operator = diagonal[:, None] * symmetric_pullback / diagonal[None, :]
    return operator, disk_nodes, eigenvectors


def solve_diagonal_similarity_sdp(operator: np.ndarray) -> DiagonalCertificate:
    """Minimize the similarity square over positive diagonal metrics."""

    dimension = operator.shape[0]
    diagonal = cp.Variable(dimension)
    bound = cp.Variable()
    metric = cp.diag(diagonal)
    constraints = [
        diagonal >= 1.0,
        diagonal <= bound,
        metric - operator.T @ metric @ operator >> 0,
    ]
    problem = cp.Problem(cp.Minimize(bound), constraints)

    solver = "CLARABEL"
    try:
        problem.solve(
            solver=solver,
            tol_gap_abs=1e-10,
            tol_gap_rel=1e-10,
            tol_feas=1e-10,
            max_iter=1_000,
        )
    except cp.error.SolverError:
        pass
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        solver = "SCS"
        problem.solve(solver=solver, eps=1e-8, max_iters=200_000, verbose=False)
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"diagonal similarity SDP failed: {problem.status}")
    return DiagonalCertificate(
        bound=float(bound.value),
        diagonal=np.asarray(diagonal.value),
        solver=solver,
    )


def predicted_similarity_square(dimension: int, ellipse_parameter: float) -> float:
    """Return ``k(c^(2p-2))/c^(p-1)`` without a small-nome inversion."""

    length = dimension - 1
    descended_modulus = elliptic_modulus_from_nome(ellipse_parameter ** (2 * length))
    return descended_modulus / ellipse_parameter**length


def explicit_low_size_metric(
    dimension: int, ellipse_parameter: float
) -> np.ndarray | None:
    """Return the exact-form metric diagonal for dimensions three and four."""

    modulus = elliptic_modulus_from_nome(ellipse_parameter**2)
    if dimension == 3:
        upper_square = modulus / (2.0 * ellipse_parameter)
        lower_square = ellipse_parameter * modulus / 2.0
        discriminant = max(0.0, 1.0 - 4.0 * upper_square * lower_square)
        middle = (1.0 - sqrt(discriminant)) / (2.0 * lower_square)
        return np.array([1.0, middle, middle**2])
    if dimension == 4:
        quarter_period = float(ellipk(modulus**2))
        eta = float(ellipj(quarter_period / 3.0, modulus**2)[0])
        return np.array(
            [
                1.0,
                modulus * eta / ellipse_parameter,
                modulus**2 * eta**3 / ellipse_parameter**2,
                modulus**3 * eta**4 / ellipse_parameter**3,
            ]
        )
    return None


def explicit_certificate_error(
    operator: np.ndarray,
    dimension: int,
    ellipse_parameter: float,
    predicted_bound: float,
) -> float | None:
    """Measure feasibility and endpoint identities for the low-size formulas."""

    diagonal = explicit_low_size_metric(dimension, ellipse_parameter)
    if diagonal is None:
        return None
    metric = np.diag(diagonal)
    defect = metric - operator.T @ metric @ operator
    eigenvalues = np.linalg.eigvalsh((defect + defect.T) / 2.0)
    errors = [
        abs(diagonal[0] - 1.0),
        abs(diagonal[-1] - predicted_bound),
        max(0.0, 1.0 - float(diagonal.min())),
        max(0.0, float(diagonal.max()) - predicted_bound),
        max(0.0, -float(eigenvalues.min())),
        float(np.max(np.abs(eigenvalues[:-1]))),
    ]
    return max(errors)


def alternation_reversal_error(
    dimension: int,
    ellipse_parameter: float,
    eigenvectors: np.ndarray,
) -> float:
    """Check the DCT alternation identity behind the exact lower bound."""

    length = dimension - 1
    amplitude = sqrt(
        elliptic_modulus_from_nome(ellipse_parameter ** (2 * length))
    )
    alternating = (-1.0) ** np.arange(dimension)
    symmetric_alternation = (eigenvectors * alternating) @ eigenvectors.T
    diagonal = ellipse_parameter ** (np.arange(dimension) / 2.0)
    evaluated = (
        amplitude
        * diagonal[:, None]
        * symmetric_alternation
        / diagonal[None, :]
    )
    reversal = np.fliplr(np.eye(dimension))
    expected = amplitude * diagonal[:, None] * reversal / diagonal[None, :]
    scale = max(1.0, np.linalg.norm(expected, 2))
    return float(
        min(
            np.linalg.norm(evaluated - expected, 2),
            np.linalg.norm(evaluated + expected, 2),
        )
        / scale
    )


def make_record(dimension: int, ellipse_parameter: float) -> AxisRecord:
    operator, _, eigenvectors = conformal_crabb_operator(
        dimension, ellipse_parameter
    )
    predicted = predicted_similarity_square(dimension, ellipse_parameter)
    unrestricted = solve_similarity_sdp(operator)
    diagonal = solve_diagonal_similarity_sdp(operator)

    unrestricted_metric = (
        unrestricted.metric + unrestricted.metric.conj().T
    ) / 2.0
    unrestricted_offdiagonal = unrestricted_metric - np.diag(
        np.diag(unrestricted_metric)
    )
    unrestricted_offdiagonal_error = float(
        np.linalg.norm(unrestricted_offdiagonal, 2)
        / max(1.0, np.linalg.norm(unrestricted_metric, 2))
    )

    diagonal_values = diagonal.diagonal
    reflection_products = (
        diagonal_values * diagonal_values[::-1] / diagonal.bound
    )
    reflection_error = float(np.max(np.abs(reflection_products - 1.0)))

    metric = np.diag(diagonal_values)
    stein_defect = metric - operator.T @ metric @ operator
    defect_eigenvalues = np.linalg.eigvalsh((stein_defect + stein_defect.T) / 2)
    rank_one_error = float(
        np.max(np.abs(defect_eigenvalues[:-1]))
        / max(1.0, abs(defect_eigenvalues[-1]))
    )

    return AxisRecord(
        dimension=dimension,
        ellipse_parameter=ellipse_parameter,
        predicted_bound=predicted,
        unrestricted_bound=unrestricted.bound,
        diagonal_bound=diagonal.bound,
        unrestricted_dual_ratio=unrestricted.dual_ratio,
        unrestricted_prediction_error=abs(unrestricted.bound - predicted),
        diagonal_prediction_error=abs(diagonal.bound - predicted),
        unrestricted_metric_offdiagonal=unrestricted_offdiagonal_error,
        diagonal_metric_reflection_error=reflection_error,
        diagonal_stein_rank_one_error=rank_one_error,
        alternation_reversal_error=alternation_reversal_error(
            dimension, ellipse_parameter, eigenvectors
        ),
        explicit_low_size_certificate_error=explicit_certificate_error(
            operator, dimension, ellipse_parameter, predicted
        ),
        unrestricted_solver=unrestricted.solver,
        diagonal_solver=diagonal.solver,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=10)
    parser.add_argument(
        "--ellipse-parameters",
        type=float,
        nargs="+",
        default=DEFAULT_ELLIPSE_PARAMETERS,
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    if any(not 0 < value < 1 for value in args.ellipse_parameters):
        raise ValueError("all ellipse parameters must lie in (0,1)")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    records = []
    with args.output.open("w") as output:
        for dimension in range(args.minimum_size, args.maximum_size + 1):
            for ellipse_parameter in args.ellipse_parameters:
                record = make_record(dimension, ellipse_parameter)
                records.append(record)
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")

    maximum_prediction_error = max(
        max(record.unrestricted_prediction_error, record.diagonal_prediction_error)
        for record in records
    )
    maximum_alternation_error = max(
        record.alternation_reversal_error for record in records
    )
    maximum_explicit_error = max(
        (
            record.explicit_low_size_certificate_error
            for record in records
            if record.explicit_low_size_certificate_error is not None
        ),
        default=0.0,
    )
    if maximum_prediction_error > 2e-5:
        raise RuntimeError(
            f"similarity prediction regression failed: {maximum_prediction_error:.3e}"
        )
    if maximum_alternation_error > 2e-12:
        raise RuntimeError(
            f"alternation identity regression failed: {maximum_alternation_error:.3e}"
        )
    if maximum_explicit_error > 2e-11:
        raise RuntimeError(
            f"explicit low-size certificate failed: {maximum_explicit_error:.3e}"
        )


if __name__ == "__main__":
    main()
