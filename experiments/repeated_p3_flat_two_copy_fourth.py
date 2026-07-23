#!/usr/bin/env python3
"""Prove fourth-order persistence at L100's weighted normal center."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from repeated_p3_second_support import generators


def support(matrix: sp.Matrix, boundary: sp.Symbol) -> sp.Matrix:
    """Return the Hermitian support matrix on the formal unit circle."""

    return (
        matrix / boundary
        + boundary * matrix.conjugate().T
    ) / 2


def projected_product(
    factors: Sequence[sp.Matrix],
    resolvent: sp.Matrix,
    top_adjoint: sp.Matrix,
    top_vectors: sp.Matrix,
) -> sp.Matrix:
    """Compress a product with one reduced resolvent between each factor."""

    product = factors[0]
    for factor in factors[1:]:
        product = product * resolvent * factor
    return top_adjoint * product * top_vectors


def traceless(matrix: sp.Matrix) -> sp.Matrix:
    """Return the traceless part of a two-by-two matrix."""

    return matrix - sp.trace(matrix) * sp.eye(2) / 2


def assert_zero_matrix(matrix: sp.Matrix, message: str) -> None:
    """Prove a small rational-function matrix is zero entry by entry."""

    for entry in matrix:
        if sp.factor(sp.cancel(entry)) != 0:
            raise AssertionError(message)


def main() -> None:
    root_two = sp.sqrt(2)
    boundary = sp.symbols("boundary", nonzero=True)
    diagonal = sp.symbols("diagonal", real=True)

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    generator_first, generator_second = generators()[0]

    normal_copy = sp.diag(diagonal, -diagonal)
    edge_copy = sp.Matrix([[0, 1], [0, 0]])

    def flat_direction(copy: sp.Matrix) -> sp.Matrix:
        return (
            sp.kronecker_product(
                copy.conjugate().T,
                generator_first,
            )
            + sp.kronecker_product(copy, generator_second)
        )

    common_center = 3 * diagonal**2 / (8 * root_two)
    common_block = common_center * sp.Matrix(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    normal_support = support(
        flat_direction(normal_copy),
        boundary,
    )
    edge_support = support(
        flat_direction(edge_copy),
        boundary,
    )
    second_support = support(
        sp.diag(common_block, common_block),
        boundary,
    )

    base_support = support(crabb, boundary)
    reduced_block = (
        sp.eye(3)
        - base_support / 4
        - 3 * base_support**2 / 4
    )
    reduced_resolvent = sp.diag(reduced_block, reduced_block)

    top_vector = sp.Matrix(
        [1 / boundary, root_two, boundary]
    ) / 2
    top_vector_adjoint = sp.Matrix(
        [[boundary, root_two, 1 / boundary]]
    ) / 2
    top_vectors = sp.zeros(6, 2)
    top_vectors[0:3, 0] = top_vector
    top_vectors[3:6, 1] = top_vector
    top_adjoint = sp.zeros(2, 6)
    top_adjoint[0, 0:3] = top_vector_adjoint
    top_adjoint[1, 3:6] = top_vector_adjoint

    def project(factors: Sequence[sp.Matrix]) -> sp.Matrix:
        return projected_product(
            factors,
            reduced_resolvent,
            top_adjoint,
            top_vectors,
        )

    effective_second = project([second_support]) + project(
        [normal_support, normal_support]
    )
    second_scalar = sp.simplify(effective_second[0, 0])
    assert_zero_matrix(
        effective_second - second_scalar * sp.eye(2),
        "the centered second support was not scalar",
    )

    fourth_sequences = (
        (
            normal_support,
            normal_support,
            normal_support,
            normal_support,
        ),
        (second_support, normal_support, normal_support),
        (normal_support, second_support, normal_support),
        (normal_support, normal_support, second_support),
        (second_support, second_support),
    )
    effective_fourth = sum(
        (project(sequence) for sequence in fourth_sequences),
        sp.zeros(2),
    )
    effective_fourth -= second_scalar * (
        top_adjoint
        * normal_support
        * reduced_resolvent**2
        * normal_support
        * top_vectors
    )
    assert_zero_matrix(
        traceless(effective_fourth),
        "the centered fourth support was not scalar",
    )

    edge_derivative = sp.zeros(2)
    normal_positions = (
        (0, 1, 2, 3),
        (1, 2),
        (0, 2),
        (0, 1),
    )
    for sequence, positions in zip(
        fourth_sequences[:-1],
        normal_positions,
        strict=True,
    ):
        for index in positions:
            differentiated = list(sequence)
            differentiated[index] = edge_support
            edge_derivative += project(differentiated)
    edge_derivative -= second_scalar * (
        top_adjoint
        * (
            edge_support
            * reduced_resolvent**2
            * normal_support
            + normal_support
            * reduced_resolvent**2
            * edge_support
        )
        * top_vectors
    )
    assert_zero_matrix(
        edge_derivative,
        "the full fourth edge derivative did not vanish",
    )

    print("PASS repeated p=3 weighted normal fourth-order persistence")
    print("the fourth effective support is scalar at the L100 center")
    print("the full fourth effective support has zero transverse derivative")


if __name__ == "__main__":
    main()
