"""Numerical regression for the proved elliptic-slice boundary inequalities.

The proofs are in ``proof/slice_boundary_theorems.md``.  This script only guards the algebraic
normal forms and stresses the inequalities away from exact symbolic simplifications.
"""

from __future__ import annotations

import argparse
from math import pi

import numpy as np

from slice_similarity_duality import modal_slice, rotation


def limit_weight_products(eigenvalue_ratio: float, left_angle: float) -> np.ndarray:
    right_angle = np.arctan(eigenvalue_ratio * np.tan(left_angle))
    bidiagonal = (
        rotation(left_angle)
        @ np.diag([1.0, eigenvalue_ratio])
        @ rotation(float(right_angle)).T
    )
    a, b, d = bidiagonal[0, 0], bidiagonal[1, 0], bidiagonal[1, 1]
    return np.array([2 * a, 2 * b, 2 * d, 4 * a * b, 4 * b * d, 8 * a * b * d])


def run_grid(c_count: int, ratio_count: int, angle_count: int) -> None:
    largest_nome_ratio = 0.0
    largest_lower_norm = 0.0
    largest_limit_product = 0.0
    largest_unproved_upper_norm = 0.0
    for c in np.geomspace(1e-3, 0.8, c_count):
        for eigenvalue_ratio in np.linspace(0.002, 0.998, ratio_count):
            for left_angle in np.linspace(0.002, pi / 2 - 0.002, angle_count):
                data = modal_slice(float(c), float(eigenvalue_ratio), float(left_angle))
                h = data.nodes[0] ** 2 / c
                nome_bound = 4 / (1 + c * c) ** 2
                largest_nome_ratio = max(largest_nome_ratio, h / nome_bound)
                upper = data.operator[:2, 2:]
                lower = data.operator[2:, :2]
                largest_lower_norm = max(largest_lower_norm, np.linalg.norm(lower, 2))
                largest_unproved_upper_norm = max(
                    largest_unproved_upper_norm, np.linalg.norm(upper, 2)
                )
                largest_limit_product = max(
                    largest_limit_product,
                    limit_weight_products(float(eigenvalue_ratio), float(left_angle)).max(),
                )

    print(f"largest nome-bound ratio:       {largest_nome_ratio:.12f}")
    print(f"largest proved lower norm:      {largest_lower_norm:.12f}")
    print(f"largest c=0 weight product:     {largest_limit_product:.12f}")
    print(f"largest upper norm (unproved):  {largest_unproved_upper_norm:.12f}")
    if largest_nome_ratio > 1 + 2e-11:
        raise AssertionError("the nome product bound failed numerically")
    if largest_lower_norm > 2 + 2e-10:
        raise AssertionError("the proved lower-block norm bound failed numerically")
    if largest_limit_product > 2 + 2e-10:
        raise AssertionError("the c=0 weighted-shift product bound failed numerically")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-count", type=int, default=20)
    parser.add_argument("--ratio-count", type=int, default=30)
    parser.add_argument("--angle-count", type=int, default=30)
    arguments = parser.parse_args()
    run_grid(arguments.c_count, arguments.ratio_count, arguments.angle_count)


if __name__ == "__main__":
    main()
