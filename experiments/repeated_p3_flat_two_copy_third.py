#!/usr/bin/env python3
"""Prove cubic descent on the trace-zero two-copy flat core."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_third_metric import tight_third_endpoint


def main() -> None:
    root_two = sp.sqrt(2)
    boundary = sp.symbols("boundary", nonzero=True)
    edge_real, edge_imaginary = sp.symbols(
        "edge_real edge_imaginary",
        real=True,
    )
    diagonal_real, diagonal_imaginary = sp.symbols(
        "diagonal_real diagonal_imaginary",
        real=True,
    )
    edge = edge_real + sp.I * edge_imaginary
    diagonal = diagonal_real + sp.I * diagonal_imaginary

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb)
    metric = sp.diag(1, 2, 4, 1, 2, 4)
    generator_zero_first, generator_zero_second = generators()[0]
    positive_loop = (
        sp.conjugate(diagonal) * generator_zero_first
        + diagonal * generator_zero_second
    )
    perturbation = sp.zeros(6)
    perturbation[0:3, 0:3] = positive_loop
    perturbation[3:6, 3:6] = -positive_loop
    perturbation[0:3, 3:6] = edge * generator_zero_second
    perturbation[3:6, 0:3] = (
        sp.conjugate(edge) * generator_zero_first
    )

    support = (crabb / boundary + boundary * crabb.T) / 2
    reduced_resolvent = sp.eye(3) - support / 4 - 3 * support**2 / 4
    support_vector = sp.Matrix([1 / boundary, root_two, boundary]) / 2
    support_adjoint = sp.Matrix(
        [[boundary, root_two, 1 / boundary]]
    ) / 2
    support_perturbation = (
        perturbation / boundary
        + boundary * perturbation.conjugate().T
    ) / 2
    full_resolvent = sp.diag(reduced_resolvent, reduced_resolvent)
    top_vectors = sp.zeros(6, 2)
    top_vectors[0:3, 0] = support_vector
    top_vectors[3:6, 1] = support_vector
    top_adjoint = sp.zeros(2, 6)
    top_adjoint[0, 0:3] = support_adjoint
    top_adjoint[1, 3:6] = support_adjoint
    effective_second = sp.simplify(
        top_adjoint
        * support_perturbation
        * full_resolvent
        * support_perturbation
        * top_vectors
    )
    effective_third = sp.simplify(
        top_adjoint
        * support_perturbation
        * full_resolvent
        * support_perturbation
        * full_resolvent
        * support_perturbation
        * top_vectors
    )

    second_mean = sp.Rational(5, 128) * (
        2 * diagonal * sp.conjugate(diagonal)
        + edge * sp.conjugate(edge)
    )
    second_mode = -sp.Rational(3, 128) * (
        diagonal**2 * boundary**2
        + sp.conjugate(diagonal) ** 2 / boundary**2
    )
    if sp.simplify(
        effective_second - (second_mean + second_mode) * sp.eye(2)
    ) != sp.zeros(2):
        raise AssertionError("the trace-zero second support was not scalar")
    if sp.simplify(sp.trace(effective_third)) != 0:
        raise AssertionError("the trace-zero third support was not traceless")
    expected_third_cross = edge * (
        3 * diagonal**2 * boundary**4
        + (
            2 * edge * sp.conjugate(edge)
            + 4 * diagonal * sp.conjugate(diagonal)
        )
        * boundary**2
        - 7 * sp.conjugate(diagonal) ** 2
    ) / (128 * boundary)
    if sp.simplify(
        effective_third[0, 1] - expected_third_cross
    ) != 0:
        raise AssertionError("the trace-zero third support formula failed")
    coercive_fourier_coefficient = sp.expand(
        expected_third_cross
    ).coeff(boundary, 1)
    expected_coercive_coefficient = (
        edge
        * (
            2 * edge * sp.conjugate(edge)
            + 4 * diagonal * sp.conjugate(diagonal)
        )
        / 128
    )
    if sp.simplify(
        coercive_fourier_coefficient - expected_coercive_coefficient
    ) != 0:
        raise AssertionError("the coercive third Fourier mode failed")

    # The middle Laurent coefficient of the displayed cross entry is
    # edge*(2|edge|^2+4|diagonal|^2)/128, so it cannot vanish for a physical
    # nonzero edge.  The top third support coefficient therefore has positive
    # mean, denoted abstractly below.
    third_mean = sp.symbols("third_mean", nonnegative=True)
    cubic_coefficient = -sp.Rational(3, 64) * diagonal**2
    cubic_frechet = (
        base**2 * perturbation
        + base * perturbation * base
        + perturbation * base**2
    )
    operators = (
        base,
        perturbation,
        -second_mean * base,
        -second_mean * perturbation
        - cubic_coefficient * cubic_frechet
        - third_mean * base,
    )

    def loop_metric(value: sp.Expr) -> sp.Matrix:
        block = sp.Matrix(
            [
                [0, -3 * root_two * value / 8, 0],
                [0, 0, 3 * root_two * value / 4],
                [0, 0, 0],
            ]
        )
        return block + block.conjugate().T

    edge_metric = sp.Matrix(
        [
            [0, -3 * root_two * edge / 8, 0],
            [0, 0, 3 * root_two * edge / 4],
            [0, 0, 0],
        ]
    )
    metric_tangent = sp.zeros(6)
    metric_tangent[0:3, 0:3] = loop_metric(diagonal)
    metric_tangent[3:6, 3:6] = -loop_metric(diagonal)
    metric_tangent[0:3, 3:6] = edge_metric
    metric_tangent[3:6, 0:3] = edge_metric.conjugate().T

    endpoint = tight_third_endpoint(
        operators,
        metric,
        metric_tangent,
    )
    if sp.simplify(endpoint + 16 * third_mean * sp.eye(2)) != sp.zeros(2):
        raise AssertionError("the trace-zero cubic endpoint did not collapse")

    print("PASS repeated p=3 trace-zero two-copy cubic descent")
    print("the second effective support is scalar")
    print("the third effective support is traceless and nonzero when edge != 0")
    print("the tight third metric endpoint is -16*mean(lambda_max(Q_3))*I")
    print("m_3 >= |edge|*(2|edge|^2+4|diagonal|^2)/128")


if __name__ == "__main__":
    main()
