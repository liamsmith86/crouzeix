#!/usr/bin/env python3
"""Audit tangent rigidity of repeated scalar-support copy blocks."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class TangentDimensions:
    """Dimensions of the relation tangent and unitary orbit."""

    full: int
    homogeneous: int
    orbit: int


def general_complex_matrix(
    size: int,
    prefix: str,
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...]]:
    """Return a general complex matrix and its real coordinates."""

    real = sp.symbols(
        " ".join(f"{prefix}_real_{index}" for index in range(size**2)),
        real=True,
    )
    imaginary = sp.symbols(
        " ".join(
            f"{prefix}_imaginary_{index}" for index in range(size**2)
        ),
        real=True,
    )
    matrix = sp.Matrix(
        size,
        size,
        lambda row, column: (
            real[row * size + column]
            + sp.I * imaginary[row * size + column]
        ),
    )
    return matrix, tuple(real) + tuple(imaginary)


def general_hermitian_matrix(
    size: int,
    prefix: str,
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...]]:
    """Return a general Hermitian matrix and its real coordinates."""

    matrix = sp.zeros(size)
    coordinates: list[sp.Symbol] = []
    for row in range(size):
        diagonal = sp.symbols(f"{prefix}_diagonal_{row}", real=True)
        matrix[row, row] = diagonal
        coordinates.append(diagonal)
        for column in range(row + 1, size):
            real, imaginary = sp.symbols(
                f"{prefix}_real_{row}_{column} "
                f"{prefix}_imaginary_{row}_{column}",
                real=True,
            )
            matrix[row, column] = real + sp.I * imaginary
            matrix[column, row] = real - sp.I * imaginary
            coordinates.extend((real, imaginary))
    return matrix, tuple(coordinates)


def real_entries(matrix: sp.Matrix) -> sp.Matrix:
    """Vectorize the real and imaginary parts of a complex matrix."""

    entries: list[sp.Expr] = []
    for entry in matrix:
        entries.extend(
            (
                sp.re(entry).expand(complex=True),
                sp.im(entry).expand(complex=True),
            )
        )
    return sp.Matrix(entries)


def scalar_relation_equations(
    first: sp.Matrix,
    second: sp.Matrix,
) -> sp.Matrix:
    """Return equations saying two matrices are scalar."""

    size = first.rows
    equations: list[sp.Expr] = []
    for row in range(size):
        for column in range(size):
            if row != column:
                equations.extend(
                    (
                        sp.re(first[row, column]).expand(complex=True),
                        sp.im(first[row, column]).expand(complex=True),
                        sp.re(second[row, column]).expand(complex=True),
                        sp.im(second[row, column]).expand(complex=True),
                    )
                )
            elif row:
                first_difference = first[row, row] - first[0, 0]
                equations.extend(
                    (
                        sp.re(first_difference).expand(complex=True),
                        sp.im(first_difference).expand(complex=True),
                        sp.re(
                            second[row, row] - second[0, 0]
                        ).expand(complex=True),
                    )
                )
    return sp.Matrix(equations)


def tangent_dimensions(
    irreducible: sp.Matrix,
    multiplicity: int,
) -> TangentDimensions:
    """Return exact tangent dimensions for one repeated irreducible."""

    base = sp.kronecker_product(irreducible, sp.eye(multiplicity))
    size = base.rows
    tangent, tangent_coordinates = general_complex_matrix(size, "tangent")
    square_derivative = base * tangent + tangent * base
    anticommutator_derivative = (
        tangent * base.conjugate().T
        + base * tangent.conjugate().T
        + tangent.conjugate().T * base
        + base.conjugate().T * tangent
    )

    homogeneous_equations = real_entries(square_derivative).col_join(
        real_entries(anticommutator_derivative)
    )
    homogeneous_rank = homogeneous_equations.jacobian(
        tangent_coordinates
    ).rank()
    homogeneous_dimension = 2 * size**2 - homogeneous_rank

    full_equations = scalar_relation_equations(
        square_derivative,
        anticommutator_derivative,
    )
    full_rank = full_equations.jacobian(tangent_coordinates).rank()
    full_dimension = 2 * size**2 - full_rank

    hermitian, hermitian_coordinates = general_hermitian_matrix(
        size,
        "generator",
    )
    orbit_tangent = (
        sp.I * hermitian * base - base * sp.I * hermitian
    )
    orbit_rank = real_entries(orbit_tangent).jacobian(
        hermitian_coordinates
    ).rank()
    return TangentDimensions(
        full=full_dimension,
        homogeneous=homogeneous_dimension,
        orbit=orbit_rank,
    )


def main() -> None:
    irreducibles = {
        "nilpotent": sp.Matrix([[0, 1], [0, 0]]),
        "invertible": sp.Matrix(
            [[sp.sqrt(2), sp.I], [sp.I, -sp.sqrt(2)]]
        ),
    }
    for name, irreducible in irreducibles.items():
        for multiplicity in range(1, 5):
            dimensions = tangent_dimensions(irreducible, multiplicity)
            expected_orbit = 3 * multiplicity**2
            if dimensions.homogeneous != expected_orbit:
                raise AssertionError(
                    f"{name} homogeneous tangent dimension changed"
                )
            if dimensions.orbit != expected_orbit:
                raise AssertionError(
                    f"{name} unitary orbit dimension changed"
                )
            if dimensions.full != expected_orbit + 3:
                raise AssertionError(
                    f"{name} full relation tangent dimension changed"
                )

    print("PASS repeated p=3 scalar-support tangent rigidity")
    print("homogeneous relation tangents equal unitary-orbit tangents")
    print("the full tangent adds exactly three two-by-two parameters")
    print("checked nilpotent and invertible strata through multiplicity four")


if __name__ == "__main__":
    main()
