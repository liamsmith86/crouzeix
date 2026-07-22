#!/usr/bin/env python3
"""Regenerate the exact p=3 second-order L21 identity from L63.

All 18 real coordinates of a general complex perturbation ``E`` are kept
symbolic.  The script derives the first two support variations, the two Schwarz
coefficients needed by the nilpotent functional calculus, the physical-gauge
operator coefficients ``G,H``, and the scalar elimination of the L62 SDP.  It
then proves that the residual from the claimed negative sum of squares is the
zero polynomial.
"""

from __future__ import annotations

import sympy as sp

from crabb_second_order_symbolic import (
    conformal_coefficients,
    symbolic_perturbation,
)


def main() -> None:
    z = sp.symbols("z", nonzero=True)
    perturbation, real_parts, imaginary_parts = symbolic_perturbation(3)
    root_two = sp.sqrt(2)
    base = sp.Matrix([[0, root_two, 0], [0, 0, root_two], [0, 0, 0]])
    metric = sp.diag(1, 2, 4)

    top_vector = sp.Matrix([sp.Rational(1, 2), z / root_two, z**2 / 2])
    zero_vector = sp.Matrix([1 / root_two, 0, -(z**2) / root_two])
    bottom_vector = sp.Matrix([sp.Rational(1, 2), -z / root_two, z**2 / 2])
    tangent, second_order = conformal_coefficients(
        base,
        perturbation,
        top_vector,
        [(sp.Integer(0), zero_vector), (-sp.Integer(1), bottom_vector)],
        z,
    )

    first_stein_forcing = sp.simplify(
        tangent.conjugate().T * metric * base + base.conjugate().T * metric * tangent
    )
    fixed_metric_tangent = sp.zeros(3)
    fixed_metric_tangent[1, 1] = first_stein_forcing[1, 1]
    fixed_metric_tangent[1, 2] = first_stein_forcing[1, 2]
    fixed_metric_tangent[2, 1] = first_stein_forcing[2, 1]
    second_stein_forcing = sp.simplify(
        second_order.conjugate().T * metric * base
        + base.conjugate().T * metric * second_order
        + tangent.conjugate().T * metric * tangent
        + tangent.conjugate().T * fixed_metric_tangent * base
        + base.conjugate().T * fixed_metric_tangent * tangent
    )

    dual_level_weight = sp.diag(0, 2, 1)
    linear_matrix = base * dual_level_weight * tangent.conjugate().T
    linear_x = sp.expand(
        linear_matrix[1, 0]
        + 2 * linear_matrix[2, 1]
        + sp.conjugate(linear_matrix[0, 1] + 2 * linear_matrix[1, 2])
        + sp.conjugate(first_stein_forcing[1, 2])
        - 2 * first_stein_forcing[1, 0]
    )
    linear_z = sp.expand(
        linear_matrix[2, 0]
        + sp.conjugate(linear_matrix[0, 2])
        - first_stein_forcing[2, 0]
    )
    constant = sp.expand(
        sp.Abs(first_stein_forcing[1, 2]) ** 2 / 2
        + 2 * sp.Abs(first_stein_forcing[1, 0]) ** 2
        + sp.Abs(first_stein_forcing[2, 0]) ** 2
        + 2 * second_stein_forcing[1, 1]
        + second_stein_forcing[2, 2]
    )
    optimum = sp.expand(
        constant - sp.Abs(linear_x) ** 2 / 8 - 3 * sp.Abs(linear_z) ** 2 / 8
    )
    target = sp.expand(
        -2 * (real_parts[1] - real_parts[5]) ** 2
        - sp.Rational(21, 4) * (real_parts[6] ** 2 + imaginary_parts[6] ** 2)
    )
    residual = sp.simplify(sp.expand_complex(optimum - target))
    if residual != 0:
        raise AssertionError(f"nonzero L63 residual: {residual}")
    print("PASS p=3 second-order identity")
    print("e(E) = -2*(Re(E01-E12))^2 - 21/4*|E20|^2")


if __name__ == "__main__":
    main()
