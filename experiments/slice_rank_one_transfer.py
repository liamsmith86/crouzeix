#!/usr/bin/env python3
"""Audit the transfer reduction for the rank-one/rank-one dual face.

The exact reduction is proved in ``proof/slice_coupled_defects.md``.  This
script has two roles:

* verify the direct and cancellation-free transfer formulas with rational
  arithmetic, including the matrix-inner identity of the node function;
* compare the transfer-norm upper bound with the similarity SDP on the
  deterministic default cases.

The optimization is floating point and is evidence only.  None of its output
is used in the exact algebra audit.
"""

from __future__ import annotations

import argparse
from math import exp, log, pi, sqrt

import numpy as np
import sympy as sp
from scipy.optimize import differential_evolution
from scipy.special import ellipk

from slice_cb_sdp import DEFAULT_CASES, solve_similarity_sdp
from slice_similarity_duality import (
    elliptic_modulus,
    modal_slice,
    modal_slice_from_weights,
)


RANK_TOLERANCE = 2e-6
OPTIMIZATION_TOLERANCE = 5e-5


def transfer_explicit(
    upper: np.ndarray,
    lower: np.ndarray,
    odd_parameter: float,
    even_parameter: float,
) -> np.ndarray:
    """Return the boundary-stable transfer matrix from formula (19)."""

    identity = np.eye(upper.shape[0])
    odd_even = upper @ lower
    even_odd = lower @ upper
    product = odd_parameter * even_parameter
    odd_resolvent = np.linalg.inv(identity - product * odd_even)
    even_resolvent = np.linalg.inv(identity - product * even_odd)
    off_diagonal_scale = sqrt(
        max(0.0, 1.0 - odd_parameter**2)
        * max(0.0, 1.0 - even_parameter**2)
    )
    return np.block(
        [
            [
                (odd_parameter * identity - even_parameter * odd_even)
                @ odd_resolvent,
                -off_diagonal_scale * odd_resolvent @ upper,
            ],
            [
                -off_diagonal_scale * even_resolvent @ lower,
                (even_parameter * identity - odd_parameter * even_odd)
                @ even_resolvent,
            ],
        ]
    )


def transfer_direct(
    operator: np.ndarray,
    odd_parameter: float,
    even_parameter: float,
) -> np.ndarray:
    """Return the transfer from its defining linear-fractional formula."""

    block_size = operator.shape[0] // 2
    odd_scale = sqrt(1.0 - odd_parameter**2)
    even_scale = sqrt(1.0 - even_parameter**2)
    parameter = np.diag(
        [odd_parameter] * block_size + [even_parameter] * block_size
    )
    defect = np.diag([odd_scale] * block_size + [even_scale] * block_size)
    identity = np.eye(operator.shape[0])
    return -np.linalg.solve(
        defect,
        (operator - parameter)
        @ np.linalg.solve(identity - parameter @ operator, defect),
    )


def transfer_core_lmi(
    operator: np.ndarray,
    odd_parameter: float,
    even_parameter: float,
) -> np.ndarray:
    """Return the denominator-free core equivalent to norm at most two."""

    block_size = operator.shape[0] // 2
    upper = operator[:block_size, block_size:]
    lower = operator[block_size:, :block_size]
    identity = np.eye(block_size)
    odd_square = odd_parameter**2
    even_square = even_parameter**2
    odd_defect = 1.0 - odd_square
    even_defect = 1.0 - even_square
    odd_block = (
        even_defect * (4.0 - odd_square) * identity
        + odd_defect * (4.0 * even_square - 1.0) * lower.T @ lower
    )
    even_block = (
        odd_defect * (4.0 - even_square) * identity
        + even_defect * (4.0 * odd_square - 1.0) * upper.T @ upper
    )
    coupling = -3.0 * (
        odd_parameter * even_defect * upper
        + even_parameter * odd_defect * lower.T
    )
    return np.block([[odd_block, coupling], [coupling.T, even_block]])


def small_gain_matrix(
    operator: np.ndarray,
    odd_parameter: float,
    even_parameter: float,
) -> np.ndarray:
    """Return the square-completion small-gain matrix from formula (33)."""

    block_size = operator.shape[0] // 2
    upper = operator[:block_size, block_size:]
    lower = operator[block_size:, :block_size]
    identity = np.eye(block_size)
    odd_coefficient = 2.0 * (4.0 - odd_parameter**2) / (
        1.0 - odd_parameter**2
    )
    even_coefficient = 2.0 * (4.0 - even_parameter**2) / (
        1.0 - even_parameter**2
    )
    odd_schur = 3.0 * odd_parameter / (4.0 - odd_parameter**2)
    even_schur = 3.0 * even_parameter / (4.0 - even_parameter**2)
    product = odd_schur * even_schur
    odd_resolvent = np.linalg.inv(identity - product * upper @ lower)
    even_resolvent = np.linalg.inv(identity - product * lower @ upper)
    common_scale = 4.0 / np.sqrt(odd_coefficient * even_coefficient)
    return np.block(
        [
            [
                common_scale * lower @ odd_resolvent,
                4.0
                * odd_schur
                * lower
                @ odd_resolvent
                @ upper
                / even_coefficient,
            ],
            [
                4.0
                * even_schur
                * upper
                @ even_resolvent
                @ lower
                / odd_coefficient,
                common_scale * upper @ even_resolvent,
            ],
        ]
    )


def block_norm_energy(matrix: np.ndarray) -> float:
    """Return the sum of squared operator norms of the four equal blocks."""

    block_size = matrix.shape[0] // 2
    return float(
        sum(
            np.linalg.norm(
                matrix[
                    row : row + block_size,
                    column : column + block_size,
                ],
                2,
            )
            ** 2
            for row in (0, block_size)
            for column in (0, block_size)
        )
    )


def exact_algebra_audit() -> None:
    """Check the transfer and matrix-inner identities symbolically."""

    upper = sp.Matrix(
        [
            [sp.Rational(1, 5), sp.Rational(2, 7)],
            [sp.Rational(-1, 4), sp.Rational(1, 3)],
        ]
    )
    lower = sp.Matrix(
        [
            [sp.Rational(2, 9), sp.Rational(-1, 6)],
            [sp.Rational(3, 8), sp.Rational(1, 7)],
        ]
    )
    zero = sp.zeros(2)
    operator = zero.row_join(upper).col_join(lower.row_join(zero))
    odd_parameter, odd_scale = sp.Rational(3, 5), sp.Rational(4, 5)
    even_parameter, even_scale = sp.Rational(5, 13), sp.Rational(12, 13)
    parameter = sp.diag(
        odd_parameter, odd_parameter, even_parameter, even_parameter
    )
    defect = sp.diag(odd_scale, odd_scale, even_scale, even_scale)
    identity_four = sp.eye(4)
    direct = -defect.inv() * (operator - parameter) * (
        identity_four - parameter * operator
    ).inv() * defect

    identity_two = sp.eye(2)
    odd_even = upper * lower
    even_odd = lower * upper
    product = odd_parameter * even_parameter
    odd_resolvent = (identity_two - product * odd_even).inv()
    even_resolvent = (identity_two - product * even_odd).inv()
    explicit = sp.BlockMatrix(
        [
            [
                (odd_parameter * identity_two - even_parameter * odd_even)
                * odd_resolvent,
                -odd_scale * even_scale * odd_resolvent * upper,
            ],
            [
                -odd_scale * even_scale * even_resolvent * lower,
                (even_parameter * identity_two - odd_parameter * even_odd)
                * even_resolvent,
            ],
        ]
    ).as_explicit()
    if sp.simplify(direct - explicit) != sp.zeros(4):
        raise AssertionError("the explicit transfer formula failed")

    core = (
        identity_four
        - operator.T * operator
        + 3
        * (identity_four - operator.T * parameter)
        * defect.inv() ** 2
        * (identity_four - parameter * operator)
    )
    congruence = (identity_four - parameter * operator).inv() * defect
    if sp.simplify(
        4 * identity_four - direct.T * direct - congruence.T * core * congruence
    ) != sp.zeros(4):
        raise AssertionError("the transfer defect identity failed")

    odd_square = odd_parameter**2
    even_square = even_parameter**2
    polynomial_odd = (
        (1 - even_square) * (4 - odd_square) * identity_two
        + (1 - odd_square)
        * (4 * even_square - 1)
        * lower.T
        * lower
    )
    polynomial_even = (
        (1 - odd_square) * (4 - even_square) * identity_two
        + (1 - even_square)
        * (4 * odd_square - 1)
        * upper.T
        * upper
    )
    polynomial_coupling = -3 * (
        odd_parameter * (1 - even_square) * upper
        + even_parameter * (1 - odd_square) * lower.T
    )
    polynomial = sp.BlockMatrix(
        [
            [polynomial_odd, polynomial_coupling],
            [polynomial_coupling.T, polynomial_even],
        ]
    ).as_explicit()
    if sp.simplify(
        (1 - odd_square) * (1 - even_square) * core - polynomial
    ) != sp.zeros(4):
        raise AssertionError("the polynomial core formula failed")

    first_column = identity_two.col_join(-odd_parameter * upper.T)
    second_column = (-even_parameter * lower.T).col_join(identity_two)
    central_defect = sp.diag(
        identity_two - lower.T * lower,
        identity_two - upper.T * upper,
    )
    defect_square_form = (
        3 * (1 - even_square) * first_column * first_column.T
        + 3 * (1 - odd_square) * second_column * second_column.T
        + (1 - odd_square) * (1 - even_square) * central_defect
    )
    if sp.simplify(polynomial - defect_square_form) != sp.zeros(4):
        raise AssertionError("the two-square transfer decomposition failed")

    odd_coefficient = 2 * (4 - odd_square) / (1 - odd_square)
    even_coefficient = 2 * (4 - even_square) / (1 - even_square)
    odd_linear = 6 * odd_parameter / (1 - odd_square)
    even_linear = 6 * even_parameter / (1 - even_square)
    completion = sp.BlockMatrix(
        [
            [
                sp.sqrt(odd_coefficient) * identity_two,
                -odd_linear * upper / sp.sqrt(odd_coefficient),
            ],
            [
                -even_linear * lower / sp.sqrt(even_coefficient),
                sp.sqrt(even_coefficient) * identity_two,
            ],
        ]
    ).as_explicit()
    loss = sp.diag(
        *([1 / sp.sqrt(even_coefficient)] * 2),
        *([1 / sp.sqrt(odd_coefficient)] * 2),
    ) * sp.diag(lower, upper)
    if sp.simplify(2 * core - completion.T * completion + 16 * loss.T * loss) != sp.zeros(4):
        raise AssertionError("the hyperbolic square completion failed")

    small_gain = 4 * loss * completion.inv()
    odd_schur = odd_linear / odd_coefficient
    even_schur = even_linear / even_coefficient
    schur_product = odd_schur * even_schur
    odd_schur_resolvent = (
        identity_two - schur_product * upper * lower
    ).inv()
    even_schur_resolvent = (
        identity_two - schur_product * lower * upper
    ).inv()
    explicit_small_gain = sp.BlockMatrix(
        [
            [
                4
                * lower
                * odd_schur_resolvent
                / sp.sqrt(odd_coefficient * even_coefficient),
                4
                * odd_schur
                * lower
                * odd_schur_resolvent
                * upper
                / even_coefficient,
            ],
            [
                4
                * even_schur
                * upper
                * even_schur_resolvent
                * lower
                / odd_coefficient,
                4
                * upper
                * even_schur_resolvent
                / sp.sqrt(odd_coefficient * even_coefficient),
            ],
        ]
    ).as_explicit()
    if sp.simplify(small_gain - explicit_small_gain) != sp.zeros(4):
        raise AssertionError("the explicit small-gain formula failed")

    a, b, s, t, z = sp.symbols("a b s t z", nonzero=True)
    denominator = 1 - a * b * z**2
    node_function = sp.Matrix(
        [[a - b * z**2, -s * t * z], [-s * t * z, b - a * z**2]]
    ) / denominator
    paraunitary = node_function.subs(z, 1 / z).T * node_function - sp.eye(2)
    for entry in paraunitary:
        numerator = sp.cancel(entry).as_numer_denom()[0]
        numerator = sp.expand(
            numerator.subs(s**2, 1 - a**2).subs(t**2, 1 - b**2)
        )
        if numerator != 0:
            raise AssertionError("the node function is not paraunitary")

    determinant = sp.factor(
        ((a - b * z**2) * (b - a * z**2) - (1 - a**2) * (1 - b**2) * z**2)
        / denominator**2
    )
    expected = (a * b - z**2) / (1 - a * b * z**2)
    if sp.simplify(determinant - expected) != 0:
        raise AssertionError("the Blaschke determinant identity failed")


def numerical_rank(matrix: np.ndarray) -> int:
    eigenvalues = np.linalg.eigvalsh((matrix + matrix.T) / 2)
    scale = max(1.0, float(np.max(abs(eigenvalues))))
    return int(np.count_nonzero(eigenvalues > RANK_TOLERANCE * scale))


def maximize_transfer(operator: np.ndarray) -> tuple[float, np.ndarray]:
    block_size = operator.shape[0] // 2
    upper = operator[:block_size, block_size:]
    lower = operator[block_size:, :block_size]

    def objective(parameters: np.ndarray) -> float:
        transfer = transfer_explicit(upper, lower, *parameters)
        return -float(np.linalg.norm(transfer, 2))

    result = differential_evolution(
        objective,
        [(-1.0, 1.0), (-1.0, 1.0)],
        seed=20260721,
        popsize=16,
        maxiter=400,
        tol=1e-10,
        polish=True,
    )
    return -float(result.fun), result.x


def run_numerical_regression() -> None:
    for weights in DEFAULT_CASES:
        data = modal_slice_from_weights(weights)
        upper = data.operator[:2, 2:]
        lower = data.operator[2:, :2]
        direct = transfer_direct(data.operator, 0.3, -0.4)
        explicit = transfer_explicit(upper, lower, 0.3, -0.4)
        formula_error = float(np.max(abs(direct - explicit)))
        if formula_error > 2e-12:
            raise AssertionError("the two floating-point transfer formulas disagree")
        core = transfer_core_lmi(data.operator, 0.3, -0.4)
        if (np.linalg.norm(explicit, 2) <= 2.0) != (
            np.linalg.eigvalsh(core).min() >= -2e-12
        ):
            raise AssertionError("the transfer norm and polynomial LMI disagree")
        small_gain = small_gain_matrix(data.operator, 0.3, -0.4)
        if (np.linalg.eigvalsh(core).min() >= -2e-12) != (
            np.linalg.norm(small_gain, 2) <= 1.0 + 2e-12
        ):
            raise AssertionError("the polynomial LMI and small-gain norm disagree")

        maximum, parameters = maximize_transfer(data.operator)
        sdp = solve_similarity_sdp(data.operator)
        witness = sdp.dual_witness
        dual_ranks = (
            numerical_rank(witness[:2, :2]),
            numerical_rank(witness[2:, 2:]),
        )
        transfer_square = maximum**2
        optimized_small_gain = small_gain_matrix(data.operator, *parameters)
        energy = block_norm_energy(optimized_small_gain)
        if transfer_square > 4.0 + OPTIMIZATION_TOLERANCE:
            raise AssertionError("the transfer search found a factor-four violation")
        if dual_ranks == (1, 1) and abs(transfer_square - sdp.bound) > OPTIMIZATION_TOLERANCE:
            raise AssertionError("the transfer and rank-one/rank-one SDP values disagree")
        print(
            f"weights={weights} dual_ranks={dual_ranks} "
            f"transfer_squared={transfer_square:.9f} sdp={sdp.bound:.9f} "
            f"parameters=({parameters[0]:+.6f},{parameters[1]:+.6f}) "
            f"block_energy={energy:.9f} formula_error={formula_error:.2e}"
        )


def maximize_global_block_energy(maximum_iterations: int) -> tuple[float, np.ndarray]:
    """Search the five live real parameters for a violation of (BE)."""

    def objective(parameters: np.ndarray) -> float:
        log_c, eigenvalue_ratio, left_angle, odd_parameter, even_parameter = (
            parameters
        )
        data = modal_slice(exp(log_c), eigenvalue_ratio, left_angle)
        small_gain = small_gain_matrix(
            data.operator, odd_parameter, even_parameter
        )
        return -block_norm_energy(small_gain)

    result = differential_evolution(
        objective,
        [
            (log(0.003), log(0.8)),
            (0.002, 0.998),
            (0.002, pi / 2 - 0.002),
            (-0.999999, 0.999999),
            (-0.999999, 0.999999),
        ],
        seed=20260721,
        popsize=20,
        maxiter=maximum_iterations,
        tol=1e-9,
        polish=True,
    )
    parameters = result.x.copy()
    parameters[0] = exp(parameters[0])
    return -float(result.fun), parameters


def envelope_operator(
    c: float,
    node_ratio: float,
    envelope_coordinate: float,
    left_angle: float,
) -> tuple[np.ndarray, float]:
    """Construct a modal operator at a point of the cubic conformal envelope."""

    modulus = elliptic_modulus(c)
    k = sqrt(modulus)
    s_zero = pi / (2.0 * float(ellipk(modulus)))
    cubic_coefficient = s_zero * (1.0 + modulus - s_zero**2) / 6.0
    lower_ratio = s_zero * node_ratio + cubic_coefficient * node_ratio**3
    upper_ratio = (
        s_zero * node_ratio + (1.0 - s_zero) * node_ratio**3
    )
    eigenvalue_ratio = lower_ratio + envelope_coordinate * (
        upper_ratio - lower_ratio
    )
    tangent = np.tan(left_angle)
    denominator = sqrt(
        (1.0 + tangent**2)
        * (1.0 + eigenvalue_ratio**2 * tangent**2)
    )
    h = k / c
    upper = sqrt(h) / denominator * np.array(
        [
            [
                1.0 + node_ratio * eigenvalue_ratio * tangent**2,
                tangent * (eigenvalue_ratio - node_ratio) / c,
            ],
            [
                c * tangent * (1.0 - node_ratio * eigenvalue_ratio),
                node_ratio + eigenvalue_ratio * tangent**2,
            ],
        ]
    )
    lower = c * sqrt(h) / denominator * np.array(
        [
            [
                1.0 + node_ratio * eigenvalue_ratio * tangent**2,
                tangent * (1.0 - node_ratio * eigenvalue_ratio) / c,
            ],
            [
                c * tangent * (eigenvalue_ratio - node_ratio),
                node_ratio + eigenvalue_ratio * tangent**2,
            ],
        ]
    )
    zero = np.zeros((2, 2))
    return np.block([[zero, upper], [lower, zero]]), eigenvalue_ratio


def maximize_envelope_block_energy(
    maximum_iterations: int,
) -> tuple[float, np.ndarray]:
    """Search the stronger cubic-envelope version of (BE)."""

    def objective(parameters: np.ndarray) -> float:
        (
            log_c,
            node_ratio,
            envelope_coordinate,
            left_angle,
            odd_parameter,
            even_parameter,
        ) = parameters
        operator, _ = envelope_operator(
            exp(log_c), node_ratio, envelope_coordinate, left_angle
        )
        return -block_norm_energy(
            small_gain_matrix(operator, odd_parameter, even_parameter)
        )

    result = differential_evolution(
        objective,
        [
            (log(0.003), log(0.63)),
            (0.002, 0.998),
            (0.0, 1.0),
            (0.002, pi / 2 - 0.002),
            (-0.999999, 0.999999),
            (-0.999999, 0.999999),
        ],
        seed=20260721,
        popsize=22,
        maxiter=maximum_iterations,
        tol=1e-9,
        polish=True,
    )
    parameters = result.x.copy()
    parameters[0] = exp(parameters[0])
    return -float(result.fun), parameters


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--global-energy",
        action="store_true",
        help="run the optional five-parameter numerical search for (BE)",
    )
    parser.add_argument(
        "--envelope-energy",
        action="store_true",
        help="search (BE) on the full rigorous cubic conformal envelope",
    )
    parser.add_argument("--maxiter", type=int, default=500)
    arguments = parser.parse_args()
    exact_algebra_audit()
    print("rank-one transfer algebra: exact")
    run_numerical_regression()
    if arguments.global_energy:
        energy, parameters = maximize_global_block_energy(arguments.maxiter)
        print(
            "global block-energy search: "
            f"energy={energy:.12f} "
            f"(c,r,u,a,b)={tuple(float(value) for value in parameters)}"
        )
        if energy > 1.0 + OPTIMIZATION_TOLERANCE:
            raise AssertionError("the global search found a block-energy violation")
    if arguments.envelope_energy:
        energy, parameters = maximize_envelope_block_energy(arguments.maxiter)
        print(
            "envelope block-energy search: "
            f"energy={energy:.12f} "
            f"(c,p,t,u,a,b)={tuple(float(value) for value in parameters)}"
        )
        if energy > 1.0 + OPTIMIZATION_TOLERANCE:
            raise AssertionError("the envelope search found a block-energy violation")


if __name__ == "__main__":
    main()
