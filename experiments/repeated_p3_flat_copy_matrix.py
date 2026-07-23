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
    common_two, common_two_adjoint = sp.symbols(
        "common_two common_two_adjoint"
    )
    adjoint_substitution = dict(
        zip(entries, adjoint_entries, strict=True)
    )
    adjoint_substitution.update(
        zip(adjoint_entries, entries, strict=True)
    )
    adjoint_substitution.update(
        {
            common_two: common_two_adjoint,
            common_two_adjoint: common_two,
        }
    )

    def formal_adjoint(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.T.xreplace(adjoint_substitution)

    copy_matrix = sp.Matrix(2, 2, entries)
    copy_adjoint = formal_adjoint(copy_matrix)
    generator_zero_first, generator_zero_second = generators()[0]
    flat_perturbation = sp.kronecker_product(
        copy_adjoint,
        generator_zero_first,
    ) + sp.kronecker_product(
        copy_matrix,
        generator_zero_second,
    )
    common_motion = common_two * sp.Matrix(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    perturbation = flat_perturbation + sp.kronecker_product(
        sp.eye(2),
        common_motion,
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
    expected_first_support = root_two * (
        common_two / boundary**2
        + common_two_adjoint * boundary**2
    ) / 4
    if sp.simplify(
        first_compression - expected_first_support * sp.eye(2)
    ) != sp.zeros(2):
        raise AssertionError("the common first compression was not scalar")

    effective_support = sp.simplify(
        top_adjoint
        * support_perturbation
        * sp.diag(reduced_resolvent, reduced_resolvent)
        * support_perturbation
        * top_vectors
    )
    flat_support = (
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
    common_scalar = (
        2 * common_two * common_two_adjoint
        - boundary**4 * common_two_adjoint**2
        - boundary**-4 * common_two**2
    ) / 16
    common_linear = root_two * (
        boundary**-1 * common_two
        - boundary**3 * common_two_adjoint
    ) / 16
    common_linear_adjoint = root_two * (
        boundary * common_two_adjoint
        - boundary**-3 * common_two
    ) / 16
    expected_support = (
        flat_support
        + common_scalar * sp.eye(2)
        + common_linear * copy_matrix
        + common_linear_adjoint * copy_adjoint
    )
    if sp.simplify(effective_support - expected_support) != sp.zeros(2):
        raise AssertionError("the flat copy-matrix support formula failed")

    print("PASS repeated p=3 flat copy-matrix reduction")
    print("E_flat = Z^* tensor X_0 + Z tensor Y_0")
    print("the Z first compression vanishes; the common-w compression is scalar")
    print("Q(q) = 5*(Z*Z^* + Z^*Z)/128")
    print("       - 3*(q^2*Z^2 + q^-2*(Z^*)^2)/128")
    print("common w adds r_w(q)*I + ell_w(q)*Z + conj(ell_w(q))*Z^*")


if __name__ == "__main__":
    main()
