#!/usr/bin/env python3
"""Regression for the odd-Blaschke coupled-face reduction.

All output is floating-point evidence.  The exact identities being audited
are proved in ``proof/slice_odd_block_reduction.md``; the final scalar
nonnegativity statement remains a conjecture.
"""

from __future__ import annotations

from math import asin, cos, pi, sin, sqrt

import numpy as np
from scipy.special import ellipk, ellipkinc

from slice_similarity_duality import elliptic_modulus, modal_slice, rotation


def inverse_ellipse_map(p: float, modulus: float) -> float:
    scale = pi / (2 * ellipk(modulus))
    return float(sin(scale * ellipkinc(asin(p), modulus)))


def cubic_envelope(p: float, modulus: float) -> tuple[float, float]:
    k = sqrt(modulus)
    slope = pi / (2 * ellipk(modulus))
    cubic = slope * (1 + k * k - slope * slope) / 6
    return slope * p + cubic * p**3, slope * p + (1 - slope) * p**3


def odd_values(k: float, p: float, parameter: float) -> tuple[float, float]:
    first = sqrt(k) * (k - parameter) / (1 - parameter * k)
    second_node_square = k * p * p
    second = (
        sqrt(k)
        * p
        * (second_node_square - parameter)
        / (1 - parameter * second_node_square)
    )
    return first, second


def reciprocal_block(
    c: float, r: float, q1: float, q2: float
) -> np.ndarray:
    return np.array(
        [
            [sqrt(r) * (q1 + q2), (r * q1 - q2) / c],
            [c * (q1 - r * q2), sqrt(r) * (q1 + q2)],
        ]
    ) / (sqrt(c) * (1 + r))


def determinant_numerator(
    c: float, r: float, q1: float, q2: float
) -> float:
    return (1 + r) ** 2 * (c * q1 * q1 * q2 * q2 + 16 * c**3) - 4 * (
        (c * c + r) ** 2 * q1 * q1
        - 2 * r * (1 - c * c) ** 2 * q1 * q2
        + (1 + c * c * r) ** 2 * q2 * q2
    )


def direct_block(
    c: float,
    r: float,
    angle: float,
    q1: float,
    q2: float,
) -> np.ndarray:
    right_angle = np.arctan(r * np.tan(angle))
    diagonal = np.diag([1.0, c])
    return (
        diagonal
        @ rotation(angle)
        @ np.diag([q1, q2])
        @ rotation(float(right_angle)).T
        @ np.linalg.inv(diagonal)
        / sqrt(c)
    )


def audit(seed: int = 20260721, count: int = 20_000) -> None:
    generator = np.random.default_rng(seed)
    largest_upper = (0.0, None)
    largest_lower = (0.0, None)
    largest_reciprocal_error = 0.0
    largest_orientation_relative_error = 0.0
    largest_determinant_error = 0.0
    largest_colligation_error = 0.0
    smallest_envelope_slack = np.inf
    smallest_determinant = np.inf
    smallest_endpoint_determinant = np.inf
    smallest_convex_linear_coefficient = np.inf

    for _ in range(count):
        c = float(np.exp(generator.uniform(np.log(1e-3), np.log(0.8))))
        r = float(generator.uniform(1e-3, 1 - 1e-3))
        angle = float(generator.uniform(1e-3, pi / 2 - 1e-3))
        parameter = float(generator.uniform(-1, 1))
        data = modal_slice(c, r, angle)
        k = float(data.nodes[0] ** 2)
        p = float(data.nodes[1] / data.nodes[0])
        q1, q2 = odd_values(k, p, parameter)

        upper = direct_block(c, r, angle, q1, q2)
        lower = (
            sqrt(c)
            * np.diag([1.0, c])
            @ data.right_rotation
            @ np.diag([q1, q2])
            @ data.left_rotation.T
            @ np.diag([1.0, 1.0 / c])
        )
        upper_norm = float(np.linalg.norm(upper, 2))
        lower_norm = float(np.linalg.norm(lower, 2))
        if upper_norm > largest_upper[0]:
            largest_upper = (upper_norm, (c, r, angle, parameter))
        if lower_norm > largest_lower[0]:
            largest_lower = (lower_norm, (c, r, angle, parameter))

        midpoint = reciprocal_block(c, r, q1, q2)
        direct_midpoint = direct_block(c, r, np.arctan(1 / sqrt(r)), q1, q2)
        largest_reciprocal_error = max(
            largest_reciprocal_error,
            float(np.max(abs(midpoint - direct_midpoint))),
        )

        tangent_square = np.tan(angle) ** 2
        denominator_square = (1 + tangent_square) * (1 + r * r * tangent_square)
        direct_polynomial = denominator_square * np.linalg.det(
            4 * np.eye(2) - upper.T @ upper
        )
        leading = (4 - q1 * q1 / c) * (4 - q2 * q2 / c)
        midpoint_polynomial = (1 + 1 / r) * (1 + r) * np.linalg.det(
            4 * np.eye(2) - midpoint.T @ midpoint
        )
        middle = r * (midpoint_polynomial - 2 * leading)
        rebuilt_polynomial = leading * (1 + r * r * tangent_square**2) + (
            middle * tangent_square
        )
        largest_orientation_relative_error = max(
            largest_orientation_relative_error,
            abs(direct_polynomial - rebuilt_polynomial)
            / max(1.0, abs(direct_polynomial), abs(rebuilt_polynomial)),
        )

        numerator = determinant_numerator(c, r, q1, q2)
        determinant_formula = numerator / (c**3 * (1 + r) ** 2)
        largest_determinant_error = max(
            largest_determinant_error,
            abs(
                np.linalg.det(4 * np.eye(2) - midpoint.T @ midpoint)
                - determinant_formula
            ),
        )

        if abs(parameter) < 0.99:
            upper_original = data.operator[:2, 2:]
            lower_original = data.operator[2:, :2]
            vector = generator.normal(size=2)
            sine = sqrt(1 - parameter * parameter)
            product = upper_original @ lower_original
            negative_vector = (product - parameter * np.eye(2)) @ vector / sine
            positive_vector = (
                parameter * lower_original - np.linalg.inv(upper_original)
            ) @ vector / sine
            even_witness = (
                lower_original
                @ np.outer(vector, vector)
                @ lower_original.T
                + np.outer(positive_vector, positive_vector)
            )
            odd_difference = np.outer(vector, vector) - (
                upper_original @ even_witness @ upper_original.T
            )
            even_difference = even_witness - (
                lower_original @ np.outer(vector, vector) @ lower_original.T
            )
            largest_colligation_error = max(
                largest_colligation_error,
                float(
                    np.max(
                        abs(odd_difference + np.outer(negative_vector, negative_vector))
                    )
                ),
                float(
                    np.max(
                        abs(even_difference - np.outer(positive_vector, positive_vector))
                    )
                ),
            )

        exact_r = inverse_ellipse_map(p, data.modulus)
        lower_envelope, upper_envelope = cubic_envelope(p, data.modulus)
        smallest_envelope_slack = min(
            smallest_envelope_slack,
            exact_r - lower_envelope,
            upper_envelope - exact_r,
        )
        smallest_determinant = min(
            smallest_determinant,
            determinant_numerator(c, exact_r, q1, q2),
        )
        smallest_endpoint_determinant = min(
            smallest_endpoint_determinant,
            determinant_numerator(c, lower_envelope, q1, q2),
            determinant_numerator(c, upper_envelope, q1, q2),
        )
        quadratic_coefficient = (
            16 * c**3
            + c * q1 * q1 * q2 * q2
            - 4 * q1 * q1
            - 4 * c**4 * q2 * q2
        )
        linear_coefficient = (
            32 * c**3
            + 2 * c * q1 * q1 * q2 * q2
            - 8 * c * c * (q1 * q1 + q2 * q2)
            + 8 * (1 - c * c) ** 2 * q1 * q2
        )
        if quadratic_coefficient > 1e-14:
            smallest_convex_linear_coefficient = min(
                smallest_convex_linear_coefficient,
                linear_coefficient,
            )

    if largest_reciprocal_error > 2e-11:
        raise AssertionError("the reciprocal block formula failed")
    if largest_orientation_relative_error > 2e-10:
        raise AssertionError("the reciprocal-quadratic orientation identity failed")
    if largest_determinant_error > 2e-7:
        raise AssertionError("the midpoint determinant formula failed")
    if largest_colligation_error > 2e-7:
        raise AssertionError("the orthogonal-colligation construction failed")
    if smallest_envelope_slack < -2e-13:
        raise AssertionError("the two-sided conformal envelope failed numerically")
    if largest_upper[0] > 2 + 2e-9 or largest_lower[0] > 2 + 2e-9:
        raise AssertionError("the odd-block target failed numerically")
    if smallest_determinant < -2e-12:
        raise AssertionError("the midpoint determinant target failed numerically")
    if smallest_endpoint_determinant < -2e-12:
        raise AssertionError("the cubic-envelope endpoint target failed numerically")
    if smallest_convex_linear_coefficient < -2e-12:
        raise AssertionError("the convex-quadratic endpoint reduction failed numerically")

    print(f"random cases: {count}")
    print(f"largest upper odd block: {largest_upper}")
    print(f"largest lower odd block: {largest_lower}")
    print(f"largest reciprocal reconstruction error: {largest_reciprocal_error:.3e}")
    print(
        "largest orientation identity relative error: "
        f"{largest_orientation_relative_error:.3e}"
    )
    print(f"largest determinant identity error: {largest_determinant_error:.3e}")
    print(f"largest colligation identity error: {largest_colligation_error:.3e}")
    print(f"smallest envelope slack: {smallest_envelope_slack:.3e}")
    print(f"smallest determinant numerator: {smallest_determinant:.3e}")
    print(f"smallest envelope-endpoint numerator: {smallest_endpoint_determinant:.3e}")
    print(
        "smallest B_r when A_r is positive: "
        f"{smallest_convex_linear_coefficient:.3e}"
    )


if __name__ == "__main__":
    audit()
