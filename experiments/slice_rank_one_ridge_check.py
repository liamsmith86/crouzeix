#!/usr/bin/env python3
"""High-precision audit of L44 at its sharp small-nome corner.

The floating-point block-energy search approaches equality when

    c, p, a, b -> 0,

where ``p`` is the conformal node ratio, ``u`` the left modal angle, and
``a,b`` the transfer parameters.  This script evaluates the exact conformal
ratio ``r=H(p)`` with mpmath and computes the four block norms in L43 at high
precision.  It also checks the closed boundary formula at ``p=a=b=0``; that
formula is independent of the modal angle ``u``.

This is a falsification/audit tool.  A passing finite scan is not a proof of
L44 on a neighbourhood of the corner.
"""

from __future__ import annotations

import argparse
from itertools import product

import mpmath as mp


def elliptic_data(c: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """Return the elliptic modulus ``k``, parameter ``m``, and ``s0``."""

    nome = c**2
    theta_two = mp.jtheta(2, 0, nome)
    theta_three = mp.jtheta(3, 0, nome)
    modulus = (theta_two / theta_three) ** 2
    parameter = modulus**2
    s_zero = 1 / theta_three**2
    return modulus, parameter, s_zero


def conformal_ratio(c: mp.mpf, node_ratio: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """Return ``H(p)`` and ``sqrt(k/c)`` without nome inversion."""

    modulus, parameter, s_zero = elliptic_data(c)
    ratio = mp.sin(s_zero * mp.ellipf(mp.asin(node_ratio), parameter))
    return ratio, mp.sqrt(modulus / c)


def modal_blocks(
    c: mp.mpf,
    node_ratio: mp.mpf,
    left_angle: mp.mpf,
) -> tuple[mp.matrix, mp.matrix]:
    """Construct the exact envelope-form blocks at the true conformal ratio."""

    eigenvalue_ratio, root_h = conformal_ratio(c, node_ratio)
    tangent = mp.tan(left_angle)
    denominator = mp.sqrt(
        (1 + tangent**2) * (1 + eigenvalue_ratio**2 * tangent**2)
    )
    upper = root_h / denominator * mp.matrix(
        [
            [
                1 + node_ratio * eigenvalue_ratio * tangent**2,
                tangent * (eigenvalue_ratio - node_ratio) / c,
            ],
            [
                c * tangent * (1 - node_ratio * eigenvalue_ratio),
                node_ratio + eigenvalue_ratio * tangent**2,
            ],
        ]
    )
    lower = c * root_h / denominator * mp.matrix(
        [
            [
                1 + node_ratio * eigenvalue_ratio * tangent**2,
                tangent * (1 - node_ratio * eigenvalue_ratio) / c,
            ],
            [
                c * tangent * (eigenvalue_ratio - node_ratio),
                node_ratio + eigenvalue_ratio * tangent**2,
            ],
        ]
    )
    return upper, lower


def operator_norm_squared(matrix: mp.matrix) -> mp.mpf:
    """Return the largest eigenvalue of ``M^T M`` for a real 2-by-2 matrix."""

    trace = sum(matrix[row, column] ** 2 for row in range(2) for column in range(2))
    determinant = mp.det(matrix)
    discriminant = max(mp.mpf("0"), trace**2 - 4 * determinant**2)
    return (trace + mp.sqrt(discriminant)) / 2


def block_energy(
    c: mp.mpf,
    node_ratio: mp.mpf,
    left_angle: mp.mpf,
    odd_parameter: mp.mpf,
    even_parameter: mp.mpf,
) -> mp.mpf:
    """Return L44's sum of four squared block norms."""

    upper, lower = modal_blocks(c, node_ratio, left_angle)
    identity = mp.eye(2)
    odd_coefficient = 2 * (4 - odd_parameter**2) / (1 - odd_parameter**2)
    even_coefficient = 2 * (4 - even_parameter**2) / (1 - even_parameter**2)
    odd_schur = 3 * odd_parameter / (4 - odd_parameter**2)
    even_schur = 3 * even_parameter / (4 - even_parameter**2)
    schur_product = odd_schur * even_schur
    odd_even = upper * lower
    even_odd = lower * upper
    odd_resolvent = (identity - schur_product * odd_even) ** -1
    even_resolvent = (identity - schur_product * even_odd) ** -1
    common_scale = 4 / mp.sqrt(odd_coefficient * even_coefficient)
    blocks = (
        common_scale * lower * odd_resolvent,
        4 * odd_schur * even_odd * even_resolvent / even_coefficient,
        4 * even_schur * odd_even * odd_resolvent / odd_coefficient,
        common_scale * upper * even_resolvent,
    )
    return sum(operator_norm_squared(block) for block in blocks)


def boundary_energy(c: mp.mpf) -> mp.mpf:
    """Return the exact ``p=a=b=0`` energy, valid for every modal angle."""

    modulus, _, _ = elliptic_data(c)
    return modulus * (1 + c**2) / (4 * c)


def centered_node(c: mp.mpf) -> mp.mpf:
    """Return L34's old odd-block ridge center for a staleness cross-check."""

    modulus, _, _ = elliptic_data(c)
    g = mp.sqrt(modulus / c)
    return (g**2 - 4 * c**2) / (2 * g**2 * (1 - 2 * c**2 * g))


def run_scan(digits: int) -> None:
    mp.mp.dps = digits
    coefficient_values = (mp.mpf("0"), mp.mpf("0.5"), mp.mpf("1"), mp.mpf("2"), mp.mpf("4"))
    transfer_values = (mp.mpf("-2"), mp.mpf("0"), mp.mpf("2"))
    for exponent in (2, 4, 8, 12, 20):
        c = mp.power(10, -exponent)
        exact_boundary = boundary_energy(c)
        evaluated_boundaries = [
            block_energy(c, mp.mpf("0"), angle, mp.mpf("0"), mp.mpf("0"))
            for angle in (mp.mpf("0"), mp.mpf("0.1"), mp.mpf("0.4"))
        ]
        if max(abs(exact_boundary - value) for value in evaluated_boundaries) > mp.power(
            10, -(digits - 15)
        ):
            raise AssertionError("the closed ridge formula and block calculation disagree")

        maximum = mp.mpf("-inf")
        maximizer: tuple[mp.mpf, ...] | None = None
        for p_scale, angle_scale, odd_scale, even_scale in product(
            coefficient_values,
            coefficient_values,
            transfer_values,
            transfer_values,
        ):
            energy = block_energy(
                c,
                p_scale * c,
                angle_scale * c,
                odd_scale * c,
                even_scale * c,
            )
            if energy > maximum:
                maximum = energy
                maximizer = (p_scale, angle_scale, odd_scale, even_scale)
        if maximum > 1:
            raise AssertionError("the high-precision ridge scan found L44 > 1")

        old_center = centered_node(c)
        old_center_energies = [
            block_energy(
                c,
                old_center + c**4 * offset,
                c,
                mp.mpf("0"),
                mp.mpf("0"),
            )
            for offset in (mp.mpf("-16"), mp.mpf("0"), mp.mpf("16"))
        ]
        largest_old_center_energy = max(old_center_energies)
        print(
            f"c=1e-{exponent:02d} boundary={mp.nstr(exact_boundary, 30)} "
            f"margin/c^2={mp.nstr((1-exact_boundary)/c**2, 20)} "
            f"grid_max={mp.nstr(maximum, 30)} scales={tuple(map(str, maximizer or ())) } "
            f"old_pstar_margin/c^2="
            f"{mp.nstr((1-largest_old_center_energy)/c**2, 20)}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--digits", type=int, default=80)
    arguments = parser.parse_args()
    if arguments.digits < 50:
        raise ValueError("the sharp-ridge audit requires at least 50 decimal digits")
    run_scan(arguments.digits)


if __name__ == "__main__":
    main()
