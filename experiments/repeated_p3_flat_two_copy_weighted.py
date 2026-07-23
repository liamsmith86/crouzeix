#!/usr/bin/env python3
"""Prove weighted cubic descent at a two-copy flat terminal block."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators
from repeated_p3_third_metric import tight_third_endpoint


def constant_laurent_coefficient(
    expression: sp.Expr,
    boundary: sp.Symbol,
    shift: int = 8,
) -> sp.Expr:
    """Return the constant coefficient of a bounded Laurent polynomial."""

    shifted = sp.expand(sp.cancel(expression) * boundary**shift)
    return shifted.coeff(boundary, shift)


def loop_direction(value: sp.Expr) -> sp.Matrix:
    """Return the generator-zero loop associated with one copy scalar."""

    generator_first, generator_second = generators()[0]
    return (
        sp.conjugate(value) * generator_first
        + value * generator_second
    )


def loop_metric(value: sp.Expr) -> sp.Matrix:
    """Return the canonical first metric for one generator-zero loop."""

    root_two = sp.sqrt(2)
    block = sp.Matrix(
        [
            [0, -3 * root_two * value / 8, 0],
            [0, 0, 3 * root_two * value / 4],
            [0, 0, 0],
        ]
    )
    return block + block.conjugate().T


def main() -> None:
    root_two = sp.sqrt(2)
    boundary = sp.symbols("boundary", nonzero=True)
    edge = sp.symbols("edge", real=True, nonnegative=True)
    diagonal_real, diagonal_imaginary = sp.symbols(
        "diagonal_real diagonal_imaginary",
        real=True,
    )
    scalar_real, scalar_imaginary = sp.symbols(
        "scalar_real scalar_imaginary",
        real=True,
    )
    common_real, common_imaginary = sp.symbols(
        "common_real common_imaginary",
        real=True,
    )
    tangent_diagonal_real, tangent_diagonal_imaginary = sp.symbols(
        "tangent_diagonal_real tangent_diagonal_imaginary",
        real=True,
    )
    tangent_edge = sp.symbols("tangent_edge", real=True)
    diagonal = diagonal_real + sp.I * diagonal_imaginary
    scalar = scalar_real + sp.I * scalar_imaginary
    common = common_real + sp.I * common_imaginary
    tangent_diagonal = (
        tangent_diagonal_real + sp.I * tangent_diagonal_imaginary
    )

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.diag(crabb, crabb)
    metric = sp.diag(1, 2, 4, 1, 2, 4)
    generator_first, generator_second = generators()[0]
    first_direction = sp.zeros(6)
    first_direction[0:3, 0:3] = loop_direction(diagonal)
    first_direction[3:6, 3:6] = -loop_direction(diagonal)
    first_direction[0:3, 3:6] = edge * generator_second
    first_direction[3:6, 0:3] = edge * generator_first

    common_direction = common * sp.Matrix(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    second_direction = sp.zeros(6)
    second_direction[0:3, 0:3] = (
        loop_direction(scalar) + common_direction
    )
    second_direction[3:6, 3:6] = (
        loop_direction(scalar) + common_direction
    )
    second_direction[0:3, 0:3] += loop_direction(
        tangent_diagonal
    )
    second_direction[3:6, 3:6] -= loop_direction(
        tangent_diagonal
    )
    second_direction[0:3, 3:6] = tangent_edge * generator_second
    second_direction[3:6, 0:3] = tangent_edge * generator_first

    support = (crabb / boundary + boundary * crabb.T) / 2
    reduced_resolvent = (
        sp.eye(3) - support / 4 - 3 * support**2 / 4
    )
    support_vector = sp.Matrix(
        [1 / boundary, root_two, boundary]
    ) / 2
    support_adjoint = sp.Matrix(
        [[boundary, root_two, 1 / boundary]]
    ) / 2
    first_support = (
        first_direction / boundary
        + boundary * first_direction.conjugate().T
    ) / 2
    second_support = (
        second_direction / boundary
        + boundary * second_direction.conjugate().T
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
        * (
            second_support
            + first_support * full_resolvent * first_support
        )
        * top_vectors
    )
    if sp.simplify(effective_second[0, 1]) != 0:
        raise AssertionError("the weighted second support was not diagonal")
    if sp.simplify(
        effective_second[0, 0] - effective_second[1, 1]
    ) != 0:
        raise AssertionError("the weighted second support was not scalar")

    effective_third = sp.simplify(
        top_adjoint
        * (
            first_support
            * full_resolvent
            * first_support
            * full_resolvent
            * first_support
            + first_support * full_resolvent * second_support
            + second_support * full_resolvent * first_support
        )
        * top_vectors
    )
    traceless_third = (
        effective_third
        - sp.trace(effective_third) * sp.eye(2) / 2
    )
    third_mean_matrix = effective_third.applyfunc(
        lambda entry: constant_laurent_coefficient(entry, boundary)
    )
    traceless_copy = sp.Matrix(
        [[diagonal, edge], [0, -diagonal]]
    )
    tangent_copy = sp.Matrix(
        [[tangent_diagonal, tangent_edge], [0, -tangent_diagonal]]
    )
    scalar_third_mean = sp.Rational(5, 64) * (
        sp.conjugate(diagonal) * tangent_diagonal
        + diagonal * sp.conjugate(tangent_diagonal)
        + edge * tangent_edge
    )
    expected_third_mean = (
        scalar_third_mean * sp.eye(2)
        + sp.Rational(5, 64)
        * (
            sp.conjugate(scalar) * traceless_copy
            + scalar * traceless_copy.conjugate().T
        )
    )
    if sp.simplify(
        third_mean_matrix - expected_third_mean
    ) != sp.zeros(2):
        raise AssertionError("the weighted third mean formula failed")
    coercive_cross_coefficient = sp.expand(
        sp.cancel(traceless_third[0, 1]) * boundary**8
    ).coeff(boundary, 9)
    expected_cross_coefficient = (
        edge
        * (
            2 * edge**2
            + 4 * diagonal * sp.conjugate(diagonal)
        )
        / 128
    )
    if sp.simplify(
        coercive_cross_coefficient - expected_cross_coefficient
    ) != 0:
        raise AssertionError("the weighted coercive Fourier mode failed")

    second_mean = sp.Rational(5, 128) * (
        2 * diagonal * sp.conjugate(diagonal) + edge**2
    )
    inverse_cubic = (
        -sp.Rational(3, 64) * diagonal**2
        + sp.conjugate(common) / root_two
    )
    cubic_frechet = (
        base**2 * first_direction
        + base * first_direction * base
        + first_direction * base**2
    )
    third_mean = sp.symbols("third_mean", real=True)
    conformal_first_real, conformal_first_imaginary = sp.symbols(
        "conformal_first_real conformal_first_imaginary",
        real=True,
    )
    conformal_first = (
        conformal_first_real + sp.I * conformal_first_imaginary
    )
    operators = (
        base,
        first_direction,
        second_direction - second_mean * base,
        -second_mean * first_direction
        - inverse_cubic * cubic_frechet
        - third_mean * base
        - conformal_first * base**2,
    )

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
    metric_tangent[3:6, 0:3] = edge_metric.T

    second_free_block = -3 * root_two * tangent_copy / 8
    endpoint = tight_third_endpoint(
        operators,
        metric,
        metric_tangent,
        second_free_block,
    )
    expected_endpoint = (
        16 * third_mean_matrix - 16 * third_mean * sp.eye(2)
    )
    if sp.simplify(endpoint - expected_endpoint) != sp.zeros(2):
        raise AssertionError("the weighted cubic endpoint failed")

    # If scalar is nonzero, the mean coefficient and the q^2 coefficient
    # cannot share a reducing line unless the Schur block is normal.
    second_fourier_matrix = traceless_third.applyfunc(
        lambda entry: sp.expand(
            sp.cancel(entry) * boundary**8
        ).coeff(boundary, 10)
    )
    expected_second_fourier = (
        -sp.Rational(3, 64) * scalar * traceless_copy
    )
    if sp.simplify(
        second_fourier_matrix - expected_second_fourier
    ) != sp.zeros(2):
        raise AssertionError("the weighted q^2 coefficient failed")

    print("PASS repeated p=3 weighted two-copy terminal descent")
    print("the order-two effective support remains scalar")
    print("the order-three endpoint is 16*(mean Q_3 - mean top(Q_3)*I)")
    print("a second-order trace-zero tangent is absorbed by its free metric")
    print("if scalar != 0, [mean Q_3,hat Q_3(2)] detects nonnormality")
    print("if scalar = 0, the uncancellable q^1 cross mode makes it strict")


if __name__ == "__main__":
    main()
