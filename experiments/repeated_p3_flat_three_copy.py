#!/usr/bin/env python3
"""Classify the three-copy nilpotent flat-matrix equality cases."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    boundary = sp.symbols("boundary", nonzero=True)
    first, middle, last = sp.symbols("first middle last")
    first_adjoint, middle_adjoint, last_adjoint = sp.symbols(
        "first_adjoint middle_adjoint last_adjoint"
    )
    adjoint_substitution = {
        first: first_adjoint,
        first_adjoint: first,
        middle: middle_adjoint,
        middle_adjoint: middle,
        last: last_adjoint,
        last_adjoint: last,
    }

    def formal_adjoint(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.T.xreplace(adjoint_substitution)

    copy_matrix = sp.Matrix(
        [[0, first, middle], [0, 0, last], [0, 0, 0]]
    )
    copy_adjoint = formal_adjoint(copy_matrix)
    square = copy_matrix**2
    expected_square = sp.zeros(3)
    expected_square[0, 2] = first * last
    if square != expected_square:
        raise AssertionError("the three-copy nilpotent square failed")

    symmetric_square = (
        copy_matrix * copy_adjoint + copy_adjoint * copy_matrix
    )
    effective_support = (
        5 * symmetric_square
        - 3
        * (
            boundary**2 * copy_matrix**2
            + boundary**-2 * copy_adjoint**2
        )
    )

    # If first*last is nonzero, a common eigenvector must be killed by both
    # E_02 and E_20 and is therefore the middle coordinate.  Its two
    # off-diagonal entries under the constant coefficient vanish iff middle
    # is zero (under the same nonzero hypothesis).
    center = sp.Matrix([0, 1, 0])
    center_residual = sp.simplify(
        symmetric_square * center
        - (first * first_adjoint + last * last_adjoint) * center
    )
    expected_residual = sp.Matrix(
        [middle * last_adjoint, 0, middle_adjoint * first]
    )
    if center_residual != expected_residual:
        raise AssertionError("the common-center obstruction failed")

    # On the path face middle=0, the center is a common top vector.  After
    # removing its eigenvalue, the negative of the remaining two-coordinate
    # block has positive diagonal and determinant 16*|first*last|^2.
    path_support = effective_support.subs(
        {middle: 0, middle_adjoint: 0}
    )
    center_eigenvalue = 5 * (
        first * first_adjoint + last * last_adjoint
    )
    complement = (0, 2)
    top_defect = sp.simplify(
        center_eigenvalue * sp.eye(2)
        - path_support.extract(complement, complement)
    )
    if sp.factor(top_defect.det()) != (
        16 * first * first_adjoint * last * last_adjoint
    ):
        raise AssertionError("the path top-defect determinant failed")
    if (
        sp.factor(top_defect[0, 0]) != 5 * last * last_adjoint
        or sp.factor(top_defect[1, 1]) != 5 * first * first_adjoint
    ):
        raise AssertionError("the path top-defect diagonal failed")

    # If first*last=0 then N^2=0, so the support is constant and its top
    # eigenspace is independent of the boundary angle.
    for zero_substitution in (
        {first: 0, first_adjoint: 0},
        {last: 0, last_adjoint: 0},
    ):
        if sp.diff(
            effective_support.subs(zero_substitution),
            boundary,
        ) != sp.zeros(3):
            raise AssertionError("a square-zero face retained an angle mode")

    print("PASS repeated p=3 three-copy flat equality classification")
    print("N = [[0,a,b],[0,0,c],[0,0,0]]")
    print("a common top branch can occur only when a*b*c = 0")
    print("a*c = 0 gives the square-zero pure-pair stratum")
    print("b = 0 with a*c != 0 gives the rank-two path/star stratum")
    print("a*b*c != 0 has a strict matrix-Jensen gap")


if __name__ == "__main__":
    main()
