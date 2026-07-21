"""Regression and exploratory grid for the elliptic-slice similarity duality.

The script checks three independent descriptions of ``T = phi(A)``:

* nodal functional calculus for the original 4x4 matrix;
* the three-parameter SVD/modal normal form;
* the primal similarity SDP and its trace-ratio dual witness.

All calculations are floating point.  They test the exact derivations in
``proof/slice_similarity_duality.md`` but do not certify the remaining uniform inequality.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
from math import asin, atan, cos, log, pi, sin, sqrt, tan

import numpy as np
from scipy.optimize import brentq
from scipy.special import ellipj, ellipk

from slice_cb_sdp import DEFAULT_CASES, conformal_matrix, solve_similarity_sdp


PARITY_ORDER = np.array([0, 2, 1, 3])


@dataclass(frozen=True)
class ModalSlice:
    c: float
    eigenvalue_ratio: float
    left_angle: float
    right_angle: float
    modulus: float
    nodes: np.ndarray
    weights: tuple[float, float, float, float]
    operator: np.ndarray
    left_rotation: np.ndarray
    right_rotation: np.ndarray


def rotation(angle: float) -> np.ndarray:
    return np.array([[cos(angle), -sin(angle)], [sin(angle), cos(angle)]])


@lru_cache(maxsize=None)
def elliptic_modulus(c: float) -> float:
    """Invert the ellipse nome relation exp(-pi K'/K) = c^2."""

    if not 0 < c < 1:
        raise ValueError("c must lie strictly between zero and one")
    target = -2 * log(c) / pi

    def residual(modulus: float) -> float:
        return float(ellipk(1 - modulus) / ellipk(modulus) - target)

    epsilon = np.finfo(float).eps
    return float(brentq(residual, epsilon, 1 - epsilon, xtol=1e-14, rtol=1e-14))


def modal_slice(c: float, eigenvalue_ratio: float, left_angle: float) -> ModalSlice:
    """Construct the exact normalized modal form from ``(c, e2/e1, u)``."""

    if not 0 < eigenvalue_ratio < 1:
        raise ValueError("the eigenvalue ratio must lie strictly between zero and one")
    if not 0 < left_angle < pi / 2:
        raise ValueError("the left singular-vector angle must lie in (0, pi/2)")

    right_angle = atan(eigenvalue_ratio * tan(left_angle))
    left_rotation = rotation(left_angle)
    right_rotation = rotation(right_angle)
    singular_values = np.diag([1.0, eigenvalue_ratio])
    bidiagonal = left_rotation @ singular_values @ right_rotation.T
    if abs(bidiagonal[0, 1]) > 2e-13:
        raise AssertionError("the modal angle relation did not produce a bidiagonal matrix")
    weights = (
        float(bidiagonal[0, 0] / sqrt(c)),
        float(bidiagonal[1, 0] / sqrt(c)),
        float(bidiagonal[1, 1] / sqrt(c)),
        c,
    )

    modulus = elliptic_modulus(c)
    k = sqrt(modulus)
    complete_integral = float(ellipk(modulus))
    elliptic_argument = (2 * complete_integral / pi) * asin(eigenvalue_ratio)
    inner_sn = float(ellipj(elliptic_argument, modulus)[0])
    nodes = np.array([sqrt(k), sqrt(k) * inner_sn])
    sigma = np.diag(nodes)
    diagonal = np.diag([1.0, c])
    diagonal_inv = np.diag([1.0, 1.0 / c])
    upper_right = (
        diagonal @ left_rotation @ sigma @ right_rotation.T @ diagonal_inv / sqrt(c)
    )
    lower_left = (
        sqrt(c) * diagonal @ right_rotation @ sigma @ left_rotation.T @ diagonal_inv
    )
    zero = np.zeros((2, 2))
    operator = np.block([[zero, upper_right], [lower_left, zero]])
    return ModalSlice(
        c=c,
        eigenvalue_ratio=eigenvalue_ratio,
        left_angle=left_angle,
        right_angle=right_angle,
        modulus=modulus,
        nodes=nodes,
        weights=weights,
        operator=operator,
        left_rotation=left_rotation,
        right_rotation=right_rotation,
    )


def modal_slice_from_weights(weights: tuple[float, float, float, float]) -> ModalSlice:
    """Extract the scale-free modal parameters from positive slice weights."""

    a1, a2, a3, c = weights
    bidiagonal = sqrt(c) * np.array([[a1, 0.0], [a2, a3]])
    left_vectors, singular_values, right_vectors_t = np.linalg.svd(bidiagonal)
    left_leader = left_vectors[:, 0]
    right_leader = right_vectors_t.T[:, 0]
    if left_leader[0] < 0:
        left_leader = -left_leader
    if right_leader[0] < 0:
        right_leader = -right_leader
    left_angle = atan(float(left_leader[1] / left_leader[0]))
    right_angle = atan(float(right_leader[1] / right_leader[0]))
    eigenvalue_ratio = float(singular_values[1] / singular_values[0])
    result = modal_slice(c, eigenvalue_ratio, left_angle)
    if abs(result.right_angle - right_angle) > 2e-12:
        raise AssertionError("the extracted SVD angles violate tan(v) = r tan(u)")
    return result


def modal_lmi_slacks(data: ModalSlice, metric: np.ndarray, bound: float) -> np.ndarray:
    """Evaluate every 2x2 LMI in the modal metric formulation."""

    odd_metric = metric[:2, :2]
    even_metric = metric[2:, 2:]
    diagonal = np.diag([1.0, data.c])
    odd_change = diagonal @ data.left_rotation
    even_change = diagonal @ data.right_rotation
    odd_gram = odd_change.T @ odd_change
    even_gram = even_change.T @ even_change
    odd_modal = odd_change.T @ odd_metric @ odd_change
    even_modal = even_change.T @ even_metric @ even_change
    sigma = np.diag(data.nodes)
    matrices = (
        odd_modal - odd_gram,
        bound * odd_gram - odd_modal,
        even_modal - even_gram,
        bound * even_gram - even_modal,
        odd_modal - data.c * sigma @ even_modal @ sigma,
        even_modal - sigma @ odd_modal @ sigma / data.c,
    )
    return np.array([np.linalg.eigvalsh(matrix).min() for matrix in matrices])


def run_regression() -> None:
    largest_reconstruction_error = 0.0
    largest_dual_gap = 0.0
    for weights in DEFAULT_CASES:
        data = modal_slice_from_weights(weights)
        direct = conformal_matrix(weights)[np.ix_(PARITY_ORDER, PARITY_ORDER)]
        reconstruction_error = float(np.max(np.abs(data.operator - direct)))
        result = solve_similarity_sdp(data.operator)
        dual_gap = abs(result.bound - result.dual_ratio)
        slacks = modal_lmi_slacks(data, result.metric, result.bound)
        largest_reconstruction_error = max(largest_reconstruction_error, reconstruction_error)
        largest_dual_gap = max(largest_dual_gap, dual_gap)
        print(
            f"weights={weights} modal_error={reconstruction_error:.2e} "
            f"t={result.bound:.9f} dual_ratio={result.dual_ratio:.9f} "
            f"dual_gap={dual_gap:.2e} modal_slack={slacks.min():+.2e}"
        )
        if reconstruction_error > 1e-9:
            raise AssertionError("modal reconstruction disagrees with the nodal calculus")
        if dual_gap > 2e-6:
            raise AssertionError("the primal and dual numerical values disagree")
    print(f"largest modal reconstruction error: {largest_reconstruction_error:.3e}")
    print(f"largest primal-dual trace-ratio gap: {largest_dual_gap:.3e}")


def run_grid(c_count: int, ratio_count: int, angle_count: int) -> None:
    c_values = np.geomspace(1e-3, 0.8, c_count)
    ratio_values = np.linspace(0.03, 0.97, ratio_count)
    angle_values = np.linspace(0.03, pi / 2 - 0.03, angle_count)
    largest_bound = -np.inf
    largest_case: tuple[float, float, float] | None = None
    fallback_count = 0
    for c in c_values:
        for eigenvalue_ratio in ratio_values:
            for left_angle in angle_values:
                data = modal_slice(float(c), float(eigenvalue_ratio), float(left_angle))
                result = solve_similarity_sdp(data.operator)
                fallback_count += result.solver != "CLARABEL"
                if result.bound > largest_bound:
                    largest_bound = result.bound
                    largest_case = (float(c), float(eigenvalue_ratio), float(left_angle))
    total = c_count * ratio_count * angle_count
    print(f"grid points: {total}")
    print(f"largest t: {largest_bound:.9f} at (c,r,u)={largest_case}")
    print(f"SCS fallbacks: {fallback_count}")
    if largest_bound >= 4.002:
        raise AssertionError("the exploratory grid produced a robust-looking t > 4")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", action="store_true", help="run the deterministic modal grid")
    parser.add_argument("--c-count", type=int, default=15)
    parser.add_argument("--ratio-count", type=int, default=9)
    parser.add_argument("--angle-count", type=int, default=9)
    arguments = parser.parse_args()
    run_regression()
    if arguments.grid:
        run_grid(arguments.c_count, arguments.ratio_count, arguments.angle_count)


if __name__ == "__main__":
    main()
