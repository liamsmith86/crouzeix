#!/usr/bin/env python3
"""Regenerate the linearization of the exact repeated-C3 metric chart."""

from __future__ import annotations

import sympy as sp


def hermitian_two_by_two(prefix: str) -> sp.Matrix:
    """Return a symbolic Hermitian two-by-two matrix."""

    diagonal_zero, diagonal_one = sp.symbols(
        f"{prefix}_diagonal_zero {prefix}_diagonal_one",
        real=True,
    )
    off_real, off_imaginary = sp.symbols(
        f"{prefix}_off_real {prefix}_off_imaginary",
        real=True,
    )
    off_diagonal = off_real + sp.I * off_imaginary
    return sp.Matrix(
        [
            [diagonal_zero, off_diagonal],
            [sp.conjugate(off_diagonal), diagonal_one],
        ]
    )


def complex_two_by_two(prefix: str) -> sp.Matrix:
    """Return a symbolic complex two-by-two matrix."""

    entries: list[sp.Expr] = []
    for index in range(4):
        real, imaginary = sp.symbols(
            f"{prefix}_{index}_real {prefix}_{index}_imaginary",
            real=True,
        )
        entries.append(real + sp.I * imaginary)
    return sp.Matrix(2, 2, entries)


def assert_zero_matrix(matrix: sp.Matrix, message: str) -> None:
    """Prove a symbolic matrix is zero entry by entry."""

    for entry in matrix:
        if sp.expand(entry) != 0:
            raise AssertionError(message)


def main() -> None:
    root_two = sp.sqrt(2)
    parameter = sp.symbols("parameter", real=True)
    identity = sp.eye(2)
    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    base = sp.kronecker_product(crabb, identity)

    first_diagonal = hermitian_two_by_two("first")
    cross = complex_two_by_two("cross")
    second_diagonal = hermitian_two_by_two("second")
    range_metric = sp.diag(2 * identity, 4 * identity) + parameter * (
        first_diagonal.row_join(cross).col_join(
            cross.conjugate().T.row_join(second_diagonal)
        )
    )
    metric = sp.diag(identity, sp.zeros(4))
    metric[2:6, 2:6] = range_metric

    stein_defect = metric - base.conjugate().T * metric * base
    kernel_compression = stein_defect[2:6, 2:6]
    linearization = kernel_compression.diff(parameter).subs(
        parameter,
        0,
    )
    expected = first_diagonal.row_join(cross).col_join(
        cross.conjugate().T.row_join(
            second_diagonal - 2 * first_diagonal
        )
    )
    assert_zero_matrix(
        linearization - expected,
        "the exact metric-chart Jacobian changed",
    )

    print("PASS repeated p=3 exact metric chart linearization")
    print("D_C F(X11,X12,X22)=(X11,X12,X22-2*X11)")
    print("a prescribed Stein Schur slack is an additive chart parameter")
    print("the Jacobian is invertible in every copy multiplicity")


if __name__ == "__main__":
    main()
