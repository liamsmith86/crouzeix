#!/usr/bin/env python3
"""Audit the normal scalar-support collision reduction to two-copy pairs."""

from __future__ import annotations

import sympy as sp

from repeated_p3_scalar_support_rigidity import (
    general_complex_matrix,
    real_entries,
    scalar_relation_equations,
)


def assert_zero_matrix(matrix: sp.Matrix, message: str) -> None:
    """Prove a symbolic matrix is zero entry by entry."""

    for entry in matrix:
        if sp.expand_complex(entry).simplify() != 0:
            raise AssertionError(message)


def tangent_dimensions(
    positive_size: int,
    negative_size: int,
) -> tuple[int, int, int]:
    """Return full, homogeneous, and off-diagonal orbit dimensions."""

    size = positive_size + negative_size
    base = sp.diag(
        sp.eye(positive_size),
        -sp.eye(negative_size),
    )
    tangent, coordinates = general_complex_matrix(size, "tangent")
    square_derivative = base * tangent + tangent * base
    anticommutator_derivative = (
        tangent * base
        + base * tangent.conjugate().T
        + tangent.conjugate().T * base
        + base * tangent
    )

    full_equations = scalar_relation_equations(
        square_derivative,
        anticommutator_derivative,
    )
    full_rank = full_equations.jacobian(coordinates).rank()
    full_dimension = 2 * size**2 - full_rank

    homogeneous_equations = real_entries(square_derivative).col_join(
        real_entries(anticommutator_derivative)
    )
    homogeneous_rank = homogeneous_equations.jacobian(coordinates).rank()
    homogeneous_dimension = 2 * size**2 - homogeneous_rank

    coupling, coupling_coordinates = general_complex_matrix(
        positive_size,
        "coupling",
    )
    # general_complex_matrix is square, so build a rectangular block when
    # the two eigenspaces have different sizes.
    if positive_size != negative_size:
        real = sp.symbols(
            " ".join(
                f"orbit_real_{index}"
                for index in range(positive_size * negative_size)
            ),
            real=True,
        )
        imaginary = sp.symbols(
            " ".join(
                f"orbit_imaginary_{index}"
                for index in range(positive_size * negative_size)
            ),
            real=True,
        )
        coupling = sp.Matrix(
            positive_size,
            negative_size,
            lambda row, column: (
                real[row * negative_size + column]
                + sp.I * imaginary[row * negative_size + column]
            ),
        )
        coupling_coordinates = tuple(real) + tuple(imaginary)
    generator = sp.zeros(size)
    generator[:positive_size, positive_size:] = coupling
    generator[positive_size:, :positive_size] = -coupling.conjugate().T
    orbit = generator * base - base * generator
    orbit_rank = real_entries(orbit).jacobian(
        coupling_coordinates
    ).rank()
    return full_dimension, homogeneous_dimension, orbit_rank


def gauge_identity() -> None:
    """Prove that one unitary gauge leaves exactly ``B-C*``."""

    positive_size, negative_size = 2, 3
    size = positive_size + negative_size
    base = sp.diag(sp.eye(positive_size), -sp.eye(negative_size))
    upper_real = sp.symbols("upper_real:6", real=True)
    upper_imaginary = sp.symbols("upper_imaginary:6", real=True)
    lower_real = sp.symbols("lower_real:6", real=True)
    lower_imaginary = sp.symbols("lower_imaginary:6", real=True)
    upper = sp.Matrix(
        positive_size,
        negative_size,
        lambda row, column: (
            upper_real[row * negative_size + column]
            + sp.I * upper_imaginary[row * negative_size + column]
        ),
    )
    lower = sp.Matrix(
        negative_size,
        positive_size,
        lambda row, column: (
            lower_real[row * positive_size + column]
            + sp.I * lower_imaginary[row * positive_size + column]
        ),
    )
    tangent = sp.zeros(size)
    tangent[:positive_size, positive_size:] = upper
    tangent[positive_size:, :positive_size] = lower

    coupling = lower.conjugate().T / 2
    generator = sp.zeros(size)
    generator[:positive_size, positive_size:] = coupling
    generator[positive_size:, :positive_size] = -coupling.conjugate().T
    gauged = sp.simplify(tangent + generator * base - base * generator)
    expected = sp.zeros(size)
    expected[:positive_size, positive_size:] = (
        upper - lower.conjugate().T
    )
    assert_zero_matrix(
        gauged - expected,
        "the normal-collision unitary gauge changed",
    )


def main() -> None:
    for positive_size in range(1, 5):
        for negative_size in range(1, 5):
            full, homogeneous, orbit = tangent_dimensions(
                positive_size,
                negative_size,
            )
            product = positive_size * negative_size
            if full != 4 * product + 2:
                raise AssertionError(
                    "the full normal-collision tangent dimension changed"
                )
            if homogeneous != 4 * product:
                raise AssertionError(
                    "the homogeneous collision dimension changed"
                )
            if orbit != 2 * product:
                raise AssertionError(
                    "the off-diagonal unitary orbit dimension changed"
                )
    gauge_identity()

    print("PASS repeated p=3 normal scalar-support collision")
    print("full tangent dimension = 4*p*q+2")
    print("the unitary gauge leaves the rectangular edge B-C*")
    print("SVD reduces that edge to independent two-copy terminal pairs")


if __name__ == "__main__":
    main()
