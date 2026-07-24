#!/usr/bin/env python3
"""Regenerate the exact real ``p=7`` sixth-face Gram certificate.

On the real length-six slice, A126 writes the canonical base polynomial
as ``P_6`` and the first cubic response as ``C_6``.  The Schur face
requires

    P_6 - (1089/290) C_6**2 >= 0.

This checker constructs an exact rational rank-seven Gram matrix on a
degree-three basis for the real phase-palindromic equality ideal.  It
verifies the polynomial identity, a positive-definite pivot core, and
the exact low-rank reconstruction.  No floating solver is used.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_full_disk_real_sos_probe import (
    degree_exponents,
    exact_polynomials,
)


REQUIRED_CONSTANT = sp.Rational(1089, 290)
PIVOT_INDICES = (0, 1, 2, 3, 7, 8, 10)


@dataclass(frozen=True)
class RealSixthCertificateRecord:
    """Exact metadata for the real length-six Gram certificate."""

    dimension: int
    length: int
    required_constant: str
    ideal_basis_dimension: int
    gram_rank: int
    pivot_indices: tuple[int, ...]
    leading_principal_minors: tuple[str, ...]
    equality_branches_annihilated: bool
    equality_ideal_basis_verified: bool
    polynomial_identity_verified: bool
    gram_reconstruction_verified: bool
    pivot_core_positive_definite: bool
    certificate_verified: bool


def equality_ideal_basis(
    variables: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    """Return a cubic basis vanishing on both real equality branches."""

    first, second, middle, fourth, fifth = variables
    return sp.Matrix(
        [
            -(first**2) * fourth + first * second * fifth,
            -first * second * middle + first * middle * fourth,
            -(first**2) * middle + first * middle * fifth,
            -first * second**2 + first * fourth**2,
            -(first**2) * second + first * fourth * fifth,
            -(first**3) + first * fifth**2,
            -first * second * fourth + second**2 * fifth,
            -(second**2) * middle + second * middle * fourth,
            -first * second * middle + second * middle * fifth,
            -(second**3) + second * fourth**2,
            -first * second**2 + second * fourth * fifth,
            -(first**2) * second + second * fifth**2,
            -second * middle**2 + middle**2 * fourth,
            -first * middle**2 + middle**2 * fifth,
            -(second**2) * middle + middle * fourth**2,
            -first * second * middle + middle * fourth * fifth,
            -(first**2) * middle + middle * fifth**2,
            -(second**2) * fourth + fourth**3,
            -first * second * fourth + fourth**2 * fifth,
            -(first**2) * fourth + fourth * fifth**2,
            -(first**2) * fifth + fifth**3,
        ]
    )


def certificate_columns() -> sp.Matrix:
    """Return the sparse ``21 x 7`` column factor of the Gram matrix."""

    columns = sp.zeros(21, 7)
    entries = {
        (0, 0): 100,
        (1, 1): 72,
        (1, 5): sp.Rational(1089, 29),
        (2, 2): sp.Rational(20655, 58),
        (2, 6): sp.Rational(4563, 29),
        (3, 3): 100,
        (3, 6): -100,
        (4, 0): 100,
        (4, 1): 216,
        (4, 5): sp.Rational(3267, 29),
        (6, 3): 100,
        (6, 6): -100,
        (7, 4): sp.Rational(8262, 145),
        (7, 6): 108,
        (8, 1): sp.Rational(1089, 29),
        (8, 5): 450,
        (10, 2): sp.Rational(4563, 29),
        (10, 3): -100,
        (10, 4): 108,
        (10, 6): 424,
        (11, 0): -100,
        (12, 0): 60,
        (13, 3): 150,
        (13, 6): -150,
        (14, 2): sp.Rational(1521, 29),
        (14, 4): sp.Rational(-3042, 145),
        (15, 0): 180,
        (16, 2): sp.Rational(-7605, 58),
        (16, 4): sp.Rational(1521, 29),
        (18, 2): sp.Rational(4563, 29),
        (18, 4): sp.Rational(-9126, 145),
    }
    for (row, column), value in entries.items():
        columns[row, column] = value
    return columns


def branch_annihilation(
    basis: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
) -> bool:
    """Check the two real phase-palindromic equality branches."""

    first, second, middle, fourth, fifth = variables
    plus = {
        fourth: second,
        fifth: first,
    }
    minus = {
        middle: 0,
        fourth: -second,
        fifth: -first,
    }
    return all(
        sp.expand(polynomial.subs(substitution)) == 0
        for polynomial in basis
        for substitution in (plus, minus)
    )


def ideal_basis_verified(
    basis: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
) -> bool:
    """Check independence and the dimension of the cubic ideal slice."""

    source_exponents = degree_exponents(len(variables), 3)
    source_monomials = [
        sp.prod(
            variable**exponent
            for variable, exponent in zip(
                variables,
                exponents,
                strict=True,
            )
        )
        for exponents in source_exponents
    ]
    first, second, middle, fourth, fifth = variables
    branch_variables = sp.symbols(
        "branch_first branch_second branch_middle",
    )
    branch_data = (
        (
            branch_variables,
            {
                first: branch_variables[0],
                second: branch_variables[1],
                middle: branch_variables[2],
                fourth: branch_variables[1],
                fifth: branch_variables[0],
            },
        ),
        (
            branch_variables[:2],
            {
                first: branch_variables[0],
                second: branch_variables[1],
                middle: 0,
                fourth: -branch_variables[1],
                fifth: -branch_variables[0],
            },
        ),
    )
    restriction_rows = []
    for target_variables, substitution in branch_data:
        restricted = [
            sp.Poly(
                monomial.subs(substitution),
                *target_variables,
            )
            for monomial in source_monomials
        ]
        for target_exponents in degree_exponents(
            len(target_variables),
            3,
        ):
            target_monomial = sp.prod(
                variable**exponent
                for variable, exponent in zip(
                    target_variables,
                    target_exponents,
                    strict=True,
                )
            )
            restriction_rows.append(
                [
                    polynomial.coeff_monomial(target_monomial)
                    for polynomial in restricted
                ]
            )
    restriction = sp.Matrix(restriction_rows)
    basis_coefficients = sp.Matrix(
        len(source_exponents),
        basis.rows,
        lambda row, column: sp.Poly(
            basis[column],
            *variables,
        ).coeff_monomial(source_monomials[row]),
    )
    restriction_rank = restriction.rank()
    return bool(
        restriction_rank == 14
        and len(source_exponents) - restriction_rank == basis.rows
        and basis_coefficients.rank() == basis.rows
    )


def build_record() -> RealSixthCertificateRecord:
    """Construct and verify the exact rational certificate."""

    variables, base, cubic = exact_polynomials()
    basis = equality_ideal_basis(variables)
    columns = certificate_columns()
    core = columns.extract(PIVOT_INDICES, range(columns.cols))
    gram = sp.simplify(columns * core.inv() * columns.T)
    leading_minors = tuple(
        sp.factor(core[:size, :size].det()) for size in range(1, core.rows + 1)
    )
    polynomial = sp.expand((basis.T * gram * basis)[0])
    target = sp.expand(base.as_expr() - REQUIRED_CONSTANT * cubic.as_expr() ** 2)
    polynomial_identity = sp.expand(polynomial - target) == 0
    reconstruction = (
        core == core.T
        and gram[:, PIVOT_INDICES] == columns
        and gram.rank() == len(PIVOT_INDICES)
    )
    positive_core = all(minor > 0 for minor in leading_minors)
    equality_branches = branch_annihilation(basis, variables)
    equality_ideal = ideal_basis_verified(basis, variables)
    verified = (
        polynomial_identity
        and reconstruction
        and positive_core
        and equality_branches
        and equality_ideal
    )
    if not verified:
        raise RuntimeError(
            "real sixth-face certificate failed: "
            f"polynomial={polynomial_identity}, "
            f"reconstruction={reconstruction}, "
            f"positive_core={positive_core}, "
            f"branches={equality_branches}, "
            f"ideal_basis={equality_ideal}"
        )
    return RealSixthCertificateRecord(
        dimension=7,
        length=6,
        required_constant=str(REQUIRED_CONSTANT),
        ideal_basis_dimension=basis.rows,
        gram_rank=gram.rank(),
        pivot_indices=PIVOT_INDICES,
        leading_principal_minors=tuple(str(minor) for minor in leading_minors),
        equality_branches_annihilated=equality_branches,
        equality_ideal_basis_verified=equality_ideal,
        polynomial_identity_verified=polynomial_identity,
        gram_reconstruction_verified=reconstruction,
        pivot_core_positive_definite=positive_core,
        certificate_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Regenerate the exact certificate."""

    args = parse_args()
    record = build_record()
    line = json.dumps(asdict(record), sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(f"{line}\n", encoding="utf-8")
    print(line, flush=True)


if __name__ == "__main__":
    main()
