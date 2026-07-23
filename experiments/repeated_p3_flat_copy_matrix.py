#!/usr/bin/env python3
"""Derive the copy-matrix formula for the repeated-C3 flat quotient."""

from __future__ import annotations

import sympy as sp

from repeated_p3_second_support import generators


def main() -> None:
    boundary = sp.symbols("boundary", nonzero=True)
    root_two = sp.sqrt(2)
    entries = sp.symbols("entry_0:4")
    adjoint_entries = sp.symbols("adjoint_entry_0:4")
    adjoint_substitution = dict(
        zip(entries, adjoint_entries, strict=True)
    )
    adjoint_substitution.update(
        zip(adjoint_entries, entries, strict=True)
    )

    def formal_adjoint(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.T.xreplace(adjoint_substitution)

    copy_matrix = sp.Matrix(2, 2, entries)
    copy_adjoint = formal_adjoint(copy_matrix)
    generator_zero_first, generator_zero_second = generators()[0]
    perturbation = sp.kronecker_product(
        copy_adjoint,
        generator_zero_first,
    ) + sp.kronecker_product(
        copy_matrix,
        generator_zero_second,
    )

    crabb = sp.Matrix(
        [[0, root_two, 0], [0, 0, root_two], [0, 0, 0]]
    )
    support = (crabb / boundary + boundary * crabb.T) / 2
    reduced_resolvent = sp.eye(3) - support / 4 - 3 * support**2 / 4
    support_vector = sp.Matrix([1 / boundary, root_two, boundary]) / 2
    support_adjoint = sp.Matrix(
        [[boundary, root_two, 1 / boundary]]
    ) / 2
    support_perturbation = (
        perturbation / boundary
        + boundary * formal_adjoint(perturbation)
    ) / 2
    top_vectors = sp.zeros(6, 2)
    top_adjoint = sp.zeros(2, 6)
    for copy in range(2):
        copy_slice = slice(3 * copy, 3 * copy + 3)
        top_vectors[copy_slice, copy] = support_vector
        top_adjoint[copy, copy_slice] = support_adjoint

    first_compression = sp.simplify(
        top_adjoint * support_perturbation * top_vectors
    )
    if first_compression != sp.zeros(2):
        raise AssertionError("the flat copy-matrix compression did not vanish")

    effective_support = sp.simplify(
        top_adjoint
        * support_perturbation
        * sp.diag(reduced_resolvent, reduced_resolvent)
        * support_perturbation
        * top_vectors
    )
    expected_support = (
        sp.Rational(5, 128)
        * (
            copy_matrix * copy_adjoint
            + copy_adjoint * copy_matrix
        )
        - sp.Rational(3, 128)
        * (
            boundary**2 * copy_matrix**2
            + boundary**-2 * copy_adjoint**2
        )
    )
    if sp.simplify(effective_support - expected_support) != sp.zeros(2):
        raise AssertionError("the flat copy-matrix support formula failed")

    print("PASS repeated p=3 flat copy-matrix reduction")
    print("E_flat = Z^* tensor X_0 + Z tensor Y_0")
    print("the first support compression vanishes identically")
    print("Q(q) = 5*(Z*Z^* + Z^*Z)/128")
    print("       - 3*(q^2*Z^2 + q^-2*(Z^*)^2)/128")


if __name__ == "__main__":
    main()
