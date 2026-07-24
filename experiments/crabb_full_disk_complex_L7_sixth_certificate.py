#!/usr/bin/env python3
"""Regenerate the exact complex ``p=8`` sixth-face certificate.

For the length-seven recentered full-disk path, let ``delta_6`` be the
sixth endpoint-excess coefficient and let ``G_3,G_4`` be its two cubic
true-normal responses.  The required Schur face is

    -2 delta_6
    - |G_3|**2 / (4 b_3)
    - |G_4|**2 / (4 b_4) >= 0.

This checker reconstructs ``delta_6`` from the endpoint recurrence and
verifies an explicit fifteen-square rational certificate on complex
Pluecker-tensor coordinates.  No floating solver is used.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_full_disk_base_jet_exact import endpoint_delta_series
from crabb_full_disk_cubic_normal_exact import (
    highest_mode_cubic_response,
    normal_curvature,
)
from crabb_full_disk_tensor_gram import (
    RankOneFactor,
    TensorIndex,
    gram_from_factors,
    pluecker_tensor_features,
    positive_factors_verified,
)


TENSOR_GROUPS: tuple[tuple[TensorIndex, ...], ...] = (
    (
        (0, 0, 4),
        (0, 1, 3),
        (1, 0, 3),
        (1, 1, 2),
        (2, 0, 2),
    ),
    (
        (0, 0, 3),
        (0, 1, 2),
        (0, 1, 4),
        (1, 0, 2),
        (1, 0, 4),
        (1, 1, 3),
        (2, 0, 3),
        (2, 1, 2),
    ),
    (
        (0, 0, 2),
        (0, 1, 5),
        (0, 2, 4),
        (1, 0, 1),
        (1, 0, 5),
        (1, 2, 3),
        (2, 0, 4),
        (2, 1, 3),
        (3, 0, 3),
    ),
    (
        (0, 0, 1),
        (0, 2, 5),
        (0, 3, 4),
        (1, 1, 5),
        (1, 2, 4),
        (2, 0, 5),
        (3, 0, 4),
        (3, 1, 3),
    ),
    (
        (0, 3, 5),
        (1, 2, 5),
        (1, 3, 4),
        (2, 1, 5),
        (3, 0, 5),
        (3, 1, 4),
        (4, 0, 4),
    ),
)
TENSOR_MAPPING = tuple(coordinate for group in TENSOR_GROUPS for coordinate in group)


@dataclass(frozen=True)
class ComplexL7SixthCertificateRecord:
    """Exact metadata for the complex length-seven certificate."""

    dimension: int
    length: int
    active_cubic_modes: tuple[int, ...]
    normal_curvatures: tuple[str, ...]
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


def group_indices(group_index: int) -> tuple[int, ...]:
    """Return global tensor-coordinate indices for one grade block."""

    start = sum(len(group) for group in TENSOR_GROUPS[:group_index])
    return tuple(range(start, start + len(TENSOR_GROUPS[group_index])))


def rank_one_factors() -> tuple[
    tuple[RankOneFactor, ...],
    tuple[RankOneFactor, ...],
]:
    """Return the nine real and six imaginary square factors."""

    rational = sp.Rational

    def factor(
        scale: sp.Rational,
        group: int,
        coefficients: tuple[sp.Rational, ...],
    ) -> RankOneFactor:
        return RankOneFactor(
            group_indices(group),
            scale,
            coefficients,
        )

    real_factors = (
        factor(
            rational(2048, 9),
            0,
            (
                rational(1),
                rational(1, 2),
                rational(0),
                rational(0),
                rational(0),
            ),
        ),
        factor(
            rational(512, 9),
            0,
            (
                rational(0),
                rational(1),
                rational(-1),
                rational(0),
                rational(3, 2),
            ),
        ),
        factor(
            rational(4608, 25),
            0,
            (
                rational(0),
                rational(0),
                rational(1),
                rational(1, 3),
                rational(0),
            ),
        ),
        factor(
            rational(2304, 25),
            1,
            (
                rational(1),
                rational(1, 3),
                rational(0),
                rational(0),
                rational(10, 9),
                rational(5, 9),
                rational(0),
                rational(0),
            ),
        ),
        factor(
            rational(256, 9),
            1,
            (
                rational(0),
                rational(0),
                rational(1),
                rational(3, 2),
                rational(-1),
                rational(0),
                rational(9, 5),
                rational(3, 5),
            ),
        ),
        factor(
            rational(64),
            2,
            (
                rational(1),
                rational(2, 3),
                rational(0),
                rational(0),
                rational(-2, 3),
                rational(0),
                rational(4, 3),
                rational(2, 3),
                rational(0),
            ),
        ),
        factor(
            rational(64),
            2,
            (
                rational(0),
                rational(0),
                rational(1),
                rational(2, 3),
                rational(0),
                rational(-2, 5),
                rational(-1),
                rational(2, 5),
                rational(6, 5),
            ),
        ),
        factor(
            rational(73984, 6957),
            3,
            (
                rational(1),
                rational(3, 2),
                rational(0),
                rational(0),
                rational(0),
                rational(-3, 2),
                rational(2),
                rational(1),
            ),
        ),
        factor(
            rational(1557504, 22525),
            4,
            (
                rational(1),
                rational(1, 3),
                rational(-5, 9),
                rational(-1, 3),
                rational(-1),
                rational(5, 9),
                rational(10, 9),
            ),
        ),
    )
    imaginary_factors = (
        factor(
            rational(2304, 25),
            1,
            (
                rational(1),
                rational(1, 3),
                rational(0),
                rational(0),
                rational(-10, 9),
                rational(-5, 9),
                rational(0),
                rational(0),
            ),
        ),
        factor(
            rational(256, 9),
            1,
            (
                rational(0),
                rational(0),
                rational(1),
                rational(-3, 2),
                rational(-1),
                rational(0),
                rational(9, 5),
                rational(3, 5),
            ),
        ),
        factor(
            rational(64),
            2,
            (
                rational(1),
                rational(-2, 3),
                rational(0),
                rational(0),
                rational(2, 3),
                rational(0),
                rational(-4, 3),
                rational(-2, 3),
                rational(0),
            ),
        ),
        factor(
            rational(64),
            2,
            (
                rational(0),
                rational(0),
                rational(1),
                rational(-2, 3),
                rational(0),
                rational(-2, 5),
                rational(-1),
                rational(2, 5),
                rational(6, 5),
            ),
        ),
        factor(
            rational(73984, 6957),
            3,
            (
                rational(1),
                rational(-3, 2),
                rational(0),
                rational(0),
                rational(0),
                rational(3, 2),
                rational(-2),
                rational(-1),
            ),
        ),
        factor(
            rational(1557504, 22525),
            4,
            (
                rational(1),
                rational(1, 3),
                rational(-5, 9),
                rational(-1, 3),
                rational(-1),
                rational(5, 9),
                rational(10, 9),
            ),
        ),
    )
    return real_factors, imaginary_factors


def polarized_target(
    variables: tuple[sp.Symbol, ...],
    conjugate_variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Build the exact polarized length-seven Schur-face polynomial."""

    coefficient_count = len(variables)
    delta = endpoint_delta_series(
        (sp.Integer(0), *variables),
        order=6,
    )[6]
    conjugate_replacement = {
        sp.conjugate(variable): conjugate_variable
        for variable, conjugate_variable in zip(
            variables,
            conjugate_variables,
            strict=True,
        )
    }
    delta = delta.xreplace(conjugate_replacement)

    def pluecker(left: int, right: int) -> sp.Expr:
        return (
            variables[left] * conjugate_variables[coefficient_count - 1 - right]
            - variables[right] * conjugate_variables[coefficient_count - 1 - left]
        )

    def conjugate_pluecker(left: int, right: int) -> sp.Expr:
        return (
            conjugate_variables[left] * variables[coefficient_count - 1 - right]
            - conjugate_variables[right] * variables[coefficient_count - 1 - left]
        )

    response_three = sp.Rational(176, 147) * (
        3 * variables[5] * pluecker(0, 2)
        - 4 * variables[3] * pluecker(0, 4)
        - 2 * variables[3] * pluecker(1, 3)
        - 2 * conjugate_variables[0] * conjugate_pluecker(0, 1)
    )
    swap = {
        **dict(zip(variables, conjugate_variables, strict=True)),
        **dict(zip(conjugate_variables, variables, strict=True)),
    }
    conjugate_response_three = response_three.xreplace(swap)
    response_four = highest_mode_cubic_response((sp.Integer(0), *variables)).xreplace(
        conjugate_replacement
    )
    conjugate_response_four = response_four.xreplace(swap)
    gain = response_three * conjugate_response_three / (
        4 * normal_curvature(7, 3)
    ) + response_four * conjugate_response_four / (4 * normal_curvature(7, 4))
    return sp.expand(-2 * delta - gain)


def build_record() -> ComplexL7SixthCertificateRecord:
    """Construct and verify the exact rational certificate."""

    variables = sp.symbols("z0:6")
    conjugate_variables = sp.symbols("zb0:6")
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
        raise RuntimeError("the complex length-seven sixth-face certificate failed")
    return ComplexL7SixthCertificateRecord(
        dimension=8,
        length=7,
        active_cubic_modes=(3, 4),
        normal_curvatures=(
            str(normal_curvature(7, 3)),
            str(normal_curvature(7, 4)),
        ),
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
