#!/usr/bin/env python3
"""Regenerate the exact complex ``p=7`` sixth-face certificate.

For the length-six recentered full-disk path, write

    [s**6] delta = -(32/225) P_6.

The true-normal Schur face requires

    P_6 - (1089/290) |C_6|**2 >= 0.

This checker reconstructs ``P_6`` from the endpoint recurrence and
verifies an explicit ten-square rational certificate on complex
Pluecker-tensor coordinates.  No floating solver is used.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_full_disk_base_jet_exact import endpoint_delta_series
from crabb_full_disk_real_sixth_certificate import REQUIRED_CONSTANT
from crabb_full_disk_tensor_gram import (
    RankOneFactor,
    gram_from_factors,
    pluecker_tensor_features,
    positive_factors_verified,
)


TENSOR_MAPPING = (
    (0, 0, 1),
    (0, 0, 2),
    (0, 0, 3),
    (0, 1, 2),
    (0, 1, 3),
    (0, 1, 4),
    (0, 2, 4),
    (1, 0, 1),
    (1, 0, 2),
    (1, 0, 3),
    (1, 0, 4),
    (1, 1, 2),
    (1, 1, 4),
    (1, 2, 3),
    (2, 0, 2),
    (2, 0, 3),
    (2, 0, 4),
    (2, 1, 2),
    (2, 1, 3),
    (3, 0, 3),
)


@dataclass(frozen=True)
class ComplexSixthCertificateRecord:
    """Exact metadata for the complex length-six certificate."""

    dimension: int
    length: int
    required_constant: str
    polarized_target_term_count: int
    tensor_feature_count: int
    real_square_count: int
    imaginary_square_count: int
    real_gram_rank: int
    imaginary_gram_rank: int
    target_conjugation_invariant: bool
    polynomial_identity_verified: bool
    gram_symmetry_verified: bool
    positive_rank_one_factors_verified: bool
    certificate_verified: bool


def rank_one_factors() -> tuple[
    tuple[RankOneFactor, ...],
    tuple[RankOneFactor, ...],
]:
    """Return the real- and imaginary-coordinate Gram factors."""

    rational = sp.Rational
    common_scale = rational(7605, 58)
    common_vector = (
        rational(1),
        rational(0),
        rational(-2, 5),
        rational(-1),
        rational(2, 5),
        rational(6, 5),
    )
    real_factors = (
        RankOneFactor(
            (0, 5, 10, 15, 17),
            rational(1),
            tuple(map(rational, (10, 10, -10, 18, 6))),
        ),
        RankOneFactor(
            (1, 7, 4, 9, 11, 14),
            rational(100),
            tuple(map(rational, (0, 1, 1, -1, 0))) + (rational(3, 2),),
        ),
        RankOneFactor(
            (1, 7, 4, 9, 11, 14),
            rational(225),
            (
                rational(1),
                rational(0),
                rational(0),
                rational(6, 5),
                rational(2, 5),
                rational(0),
            ),
        ),
        RankOneFactor(
            (2, 3, 8),
            rational(72),
            tuple(map(rational, (3, 1, 0))),
        ),
        RankOneFactor(
            (2, 3, 8),
            rational(450),
            tuple(map(rational, (0, 0, 1))),
        ),
        RankOneFactor(
            (6, 12, 13, 16, 18, 19),
            common_scale,
            common_vector,
        ),
    )
    imaginary_factors = (
        RankOneFactor(
            (0, 5, 10, 15, 17),
            rational(1),
            tuple(map(rational, (10, -10, 10, -18, -6))),
        ),
        RankOneFactor(
            (1, 7, 4, 9, 11, 14),
            rational(100),
            tuple(map(rational, (0, 1, -1, 1, 0))) + (rational(-3, 2),),
        ),
        RankOneFactor(
            (1, 7, 4, 9, 11, 14),
            rational(225),
            (
                rational(1),
                rational(0),
                rational(0),
                rational(-6, 5),
                rational(-2, 5),
                rational(0),
            ),
        ),
        RankOneFactor(
            (6, 12, 13, 16, 18, 19),
            common_scale,
            common_vector,
        ),
    )
    return real_factors, imaginary_factors


def polarized_target(
    variables: tuple[sp.Symbol, ...],
    conjugate_variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Build the exact polarized Schur-face polynomial."""

    delta = endpoint_delta_series(
        (sp.Integer(0), *variables),
        order=6,
    )[6]
    delta = delta.xreplace(
        {
            sp.conjugate(variable): conjugate_variable
            for variable, conjugate_variable in zip(
                variables,
                conjugate_variables,
                strict=True,
            )
        }
    )
    base = sp.expand(-sp.Rational(225, 32) * delta)

    def pluecker(left: int, right: int) -> sp.Expr:
        return (
            variables[left] * conjugate_variables[4 - right]
            - variables[right] * conjugate_variables[4 - left]
        )

    def conjugate_pluecker(left: int, right: int) -> sp.Expr:
        return (
            conjugate_variables[left] * variables[4 - right]
            - conjugate_variables[right] * variables[4 - left]
        )

    cubic = sp.expand(
        6 * variables[3] * pluecker(0, 3)
        - 5 * variables[4] * pluecker(0, 2)
        + 2 * variables[3] * pluecker(1, 2)
    )
    conjugate_cubic = sp.expand(
        6 * conjugate_variables[3] * conjugate_pluecker(0, 3)
        - 5 * conjugate_variables[4] * conjugate_pluecker(0, 2)
        + 2 * conjugate_variables[3] * conjugate_pluecker(1, 2)
    )
    return sp.expand(base - REQUIRED_CONSTANT * cubic * conjugate_cubic)


def build_record() -> ComplexSixthCertificateRecord:
    """Construct and verify the exact rational certificate."""

    variables = sp.symbols("z0:5")
    conjugate_variables = sp.symbols("zb0:5")
    all_variables = (*variables, *conjugate_variables)
    target = polarized_target(variables, conjugate_variables)
    features, conjugate_features = pluecker_tensor_features(
        TENSOR_MAPPING,
        variables,
        conjugate_variables,
    )
    real_coordinates = (features + conjugate_features) / 2
    imaginary_coordinates = (features - conjugate_features) / (2 * sp.I)
    real_factors, imaginary_factors = rank_one_factors()
    real_gram = gram_from_factors(real_factors, len(TENSOR_MAPPING))
    imaginary_gram = gram_from_factors(
        imaginary_factors,
        len(TENSOR_MAPPING),
    )
    certificate = sp.expand(
        (real_coordinates.T * real_gram * real_coordinates)[0]
        + (imaginary_coordinates.T * imaginary_gram * imaginary_coordinates)[0]
    )
    polynomial_identity = sp.Poly(
        certificate - target,
        *all_variables,
    ).is_zero
    swap = {
        **dict(zip(variables, conjugate_variables, strict=True)),
        **dict(zip(conjugate_variables, variables, strict=True)),
    }
    conjugation_invariant = sp.Poly(
        target.xreplace(swap) - target,
        *all_variables,
    ).is_zero
    gram_symmetry = real_gram == real_gram.T and imaginary_gram == imaginary_gram.T
    positive_factors = positive_factors_verified((*real_factors, *imaginary_factors))
    verified = bool(
        polynomial_identity
        and conjugation_invariant
        and gram_symmetry
        and positive_factors
    )
    if not verified:
        raise RuntimeError("the complex sixth-face certificate failed")
    return ComplexSixthCertificateRecord(
        dimension=7,
        length=6,
        required_constant=str(REQUIRED_CONSTANT),
        polarized_target_term_count=len(sp.Poly(target, *all_variables).terms()),
        tensor_feature_count=len(TENSOR_MAPPING),
        real_square_count=len(real_factors),
        imaginary_square_count=len(imaginary_factors),
        real_gram_rank=real_gram.rank(),
        imaginary_gram_rank=imaginary_gram.rank(),
        target_conjugation_invariant=bool(conjugation_invariant),
        polynomial_identity_verified=bool(polynomial_identity),
        gram_symmetry_verified=gram_symmetry,
        positive_rank_one_factors_verified=positive_factors,
        certificate_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Regenerate and serialize the exact certificate metadata."""

    args = parse_args()
    record = build_record()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(asdict(record), sort_keys=True)
    args.output.write_text(f"{line}\n", encoding="utf-8")
    print(line, flush=True)


if __name__ == "__main__":
    main()
