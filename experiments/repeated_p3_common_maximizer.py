#!/usr/bin/env python3
"""Classify repeated-C3 common-maximizer directions modulo unitary orbit."""

from __future__ import annotations

import sympy as sp


def real_vector(first: sp.Matrix, second: sp.Matrix) -> sp.Matrix:
    """Flatten a complex matrix pair as grouped real and imaginary entries."""

    entries: list[sp.Expr] = []
    for matrix in (first, second):
        entries.extend(sp.re(entry).expand(complex=True) for entry in matrix)
        entries.extend(sp.im(entry).expand(complex=True) for entry in matrix)
    return sp.Matrix(entries)


def single_real_vector(matrix: sp.Matrix) -> sp.Matrix:
    """Flatten one complex matrix as grouped real and imaginary entries."""

    return sp.Matrix(
        [sp.re(entry).expand(complex=True) for entry in matrix]
        + [sp.im(entry).expand(complex=True) for entry in matrix]
    )


def support_constraint_matrix() -> sp.Matrix:
    """Return the real Fourier constraints for one selected/cross copy pair."""

    boundary = sp.symbols("boundary", nonzero=True)
    root_two = sp.sqrt(2)
    support = sp.Matrix([1 / boundary, root_two, boundary])
    conjugate_support = sp.Matrix([boundary, root_two, 1 / boundary])
    x_real = sp.symbols("x_real_0:9", real=True)
    x_imaginary = sp.symbols("x_imaginary_0:9", real=True)
    y_real = sp.symbols("y_real_0:9", real=True)
    y_imaginary = sp.symbols("y_imaginary_0:9", real=True)
    first = sp.Matrix(
        3,
        3,
        [x_real[index] + sp.I * x_imaginary[index] for index in range(9)],
    )
    second = sp.Matrix(
        3,
        3,
        [y_real[index] + sp.I * y_imaginary[index] for index in range(9)],
    )
    compressed = sp.expand(
        (
            conjugate_support.T
            * (first / boundary + boundary * second.conjugate().T)
            * support
        )[0]
    )
    polynomial = sp.Poly(sp.expand(compressed * boundary**4), boundary)
    equations: list[sp.Expr] = []
    for coefficient in polynomial.all_coeffs():
        equations.extend(
            (
                sp.re(coefficient).expand(complex=True),
                sp.im(coefficient).expand(complex=True),
            )
        )
    variables = (
        *x_real,
        *x_imaginary,
        *y_real,
        *y_imaginary,
    )
    matrix, _ = sp.linear_eq_to_matrix(equations, variables)
    return matrix


def orbit_matrix() -> sp.Matrix:
    """Return cross-copy infinitesimal unitary-orbit directions."""

    root_two = sp.sqrt(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    real_entries = sp.symbols("orbit_real_0:9", real=True)
    imaginary_entries = sp.symbols("orbit_imaginary_0:9", real=True)
    generator = sp.Matrix(
        3,
        3,
        [
            real_entries[index] + sp.I * imaginary_entries[index]
            for index in range(9)
        ],
    )
    first = crabb * generator - generator * crabb
    second = -(
        crabb * generator.conjugate().T
        - generator.conjugate().T * crabb
    )
    variables = (*real_entries, *imaginary_entries)
    matrix, _ = sp.linear_eq_to_matrix(
        list(real_vector(first, second)),
        variables,
    )
    return matrix


def quotient_representatives() -> sp.Matrix:
    """Return the six real canonical residual cross directions."""

    first_generators = (
        sp.diag(-sp.Rational(3, 4), sp.Rational(1, 4), -sp.Rational(3, 4)),
        sp.Matrix([[0, -1, 0], [0, 0, 1], [0, 0, 0]]),
        sp.Matrix([[0, 0, -sp.Rational(4, 3)], [0, 0, 0], [0, 0, 0]]),
    )
    second_generators = (
        sp.Matrix([[0, 0, 1], [0, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, -1, 0], [0, 0, 1], [0, 0, 0]]),
        sp.diag(1, -sp.Rational(1, 3), 1),
    )
    columns: list[sp.Matrix] = []
    for first, second in zip(
        first_generators,
        second_generators,
        strict=True,
    ):
        columns.append(real_vector(first, second))
        columns.append(real_vector(-sp.I * first, sp.I * second))
    return sp.Matrix.hstack(*columns)


def diagonal_support_constraint_matrix() -> sp.Matrix:
    """Return the real Fourier constraints for one diagonal copy block."""

    boundary = sp.symbols("diagonal_boundary", nonzero=True)
    root_two = sp.sqrt(2)
    support = sp.Matrix([1 / boundary, root_two, boundary])
    conjugate_support = sp.Matrix([boundary, root_two, 1 / boundary])
    real_entries = sp.symbols("diagonal_real_0:9", real=True)
    imaginary_entries = sp.symbols("diagonal_imaginary_0:9", real=True)
    direction = sp.Matrix(
        3,
        3,
        [
            real_entries[index] + sp.I * imaginary_entries[index]
            for index in range(9)
        ],
    )
    compressed = sp.expand(
        (
            conjugate_support.T
            * (direction / boundary + boundary * direction.conjugate().T)
            * support
        )[0]
    )
    polynomial = sp.Poly(sp.expand(compressed * boundary**4), boundary)
    equations: list[sp.Expr] = []
    for coefficient in polynomial.all_coeffs():
        equations.extend(
            (
                sp.re(coefficient).expand(complex=True),
                sp.im(coefficient).expand(complex=True),
            )
        )
    matrix, _ = sp.linear_eq_to_matrix(
        equations,
        (*real_entries, *imaginary_entries),
    )
    return matrix


def diagonal_orbit_matrix() -> sp.Matrix:
    """Return within-copy infinitesimal unitary-orbit directions."""

    root_two = sp.sqrt(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    diagonal = sp.symbols("skew_diagonal_0:3", real=True)
    upper_real = sp.symbols("skew_upper_real_0:3", real=True)
    upper_imaginary = sp.symbols("skew_upper_imaginary_0:3", real=True)
    generator = sp.diag(*(sp.I * value for value in diagonal))
    for index, (row, column) in enumerate(((0, 1), (0, 2), (1, 2))):
        value = upper_real[index] + sp.I * upper_imaginary[index]
        generator[row, column] = value
        generator[column, row] = -sp.conjugate(value)
    direction = crabb * generator - generator * crabb
    matrix, _ = sp.linear_eq_to_matrix(
        list(single_real_vector(direction)),
        (*diagonal, *upper_real, *upper_imaginary),
    )
    return matrix


def diagonal_quotient_representatives() -> sp.Matrix:
    """Return the three real canonical diagonal support-kernel directions."""

    first_zero = sp.diag(
        -sp.Rational(3, 4),
        sp.Rational(1, 4),
        -sp.Rational(3, 4),
    )
    second_zero = sp.zeros(3)
    second_zero[0, 2] = 1
    first_one = sp.zeros(3)
    first_one[0, 1] = -1
    first_one[1, 2] = 1
    return sp.Matrix.hstack(
        single_real_vector(first_zero + second_zero),
        single_real_vector(-sp.I * first_zero + sp.I * second_zero),
        single_real_vector(first_one),
    )


def main() -> None:
    support = support_constraint_matrix()
    orbit = orbit_matrix()
    quotient = quotient_representatives()
    if support.rank() != 14:
        raise AssertionError("unexpected common-maximizer support rank")
    if orbit.rank() != 16:
        raise AssertionError("unexpected cross unitary-orbit rank")
    if support * orbit != sp.zeros(support.rows, orbit.cols):
        raise AssertionError("a unitary-orbit tangent changed the support compression")

    quotient_constraints = sp.Matrix.vstack(support, orbit.T)
    if quotient_constraints.rank() != 30:
        raise AssertionError("unexpected cross quotient dimension")
    if quotient.rank() != 6:
        raise AssertionError("the canonical quotient representatives lost rank")
    if quotient_constraints * quotient != sp.zeros(
        quotient_constraints.rows,
        quotient.cols,
    ):
        raise AssertionError("a quotient representative violated its constraints")

    diagonal_support = diagonal_support_constraint_matrix()
    diagonal_orbit = diagonal_orbit_matrix()
    diagonal_quotient = diagonal_quotient_representatives()
    if diagonal_support.rank() != 7:
        raise AssertionError("unexpected diagonal support rank")
    if diagonal_orbit.rank() != 8:
        raise AssertionError("unexpected diagonal unitary-orbit rank")
    if diagonal_support * diagonal_orbit != sp.zeros(
        diagonal_support.rows,
        diagonal_orbit.cols,
    ):
        raise AssertionError("a diagonal unitary tangent moved the support")
    if diagonal_support * diagonal_quotient != sp.zeros(
        diagonal_support.rows,
        diagonal_quotient.cols,
    ):
        raise AssertionError("a diagonal quotient tangent moved the support")
    if diagonal_orbit.T * diagonal_quotient != sp.zeros(
        diagonal_orbit.cols,
        diagonal_quotient.cols,
    ):
        raise AssertionError("a diagonal quotient tangent was not canonical")
    if sp.Matrix.hstack(diagonal_orbit, diagonal_quotient).rank() != 11:
        raise AssertionError("the diagonal support kernel was not exhausted")

    print("PASS repeated p=3 common-maximizer cross quotient")
    print("support kernel: 36 - 14 = 22 real dimensions")
    print("unitary cross orbit: 16 real dimensions")
    print("residual quotient: three complex directions")
    print("for each alpha: X=conj(alpha)*X_j, Y=alpha*Y_j")
    print("diagonal support kernel: 18 - 7 = 11 real dimensions")
    print("diagonal unitary orbit: 8 real dimensions")
    print("diagonal residual quotient: one complex loop and one real loop")


if __name__ == "__main__":
    main()
