#!/usr/bin/env python3
"""Regenerate the sparse dual reflected/normal Hessian face.

At the Crabb point, L143's prepared grade-k reflected tangent couples
the top right singular endpoint only to index ``L-k``.  The leading
normalized circular-normal variation of ``C^L`` couples it only to
index ``k-2``.  These indices are disjoint for ``k<=L/2``.

This checker derives both matrix tangents exactly and verifies that
the direct quadratic and singular-vector-coupling parts of the mixed
top-singular-value Hessian vanish.  It also regenerates L172's Cauchy
residue selection for the sole remaining endpoint/cofactor term.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_series import (
    real_circular_normal_direction,
)
from crabb_leading_endpoint_transfer import riemann_coefficient
from general_crabb_weighted_series import crabb_matrix


@dataclass(frozen=True)
class DualNormalAssociatedFaceRecord:
    """One exact associated dual Hessian audit."""

    dimension: int
    length: int
    reflected_grade: int
    normal_mode: int
    reflected_coupling_index: int
    strong_coupling_index: int
    coupling_indices_disjoint: bool
    strong_evaluation_formula_verified: bool
    direct_quadratic_term: str
    singular_vector_coupling_term: str
    quadratic_hessian_cross: str
    reflected_blaschke_modes: tuple[int, int]
    endpoint_residual_mode: int
    endpoint_product_modes: tuple[int, int]
    endpoint_residue_vanishes: bool


def prepared_reflected_tangent(
    length: int,
    grade: int,
) -> sp.Matrix:
    """Return L143's first prepared Blaschke-evaluation face."""

    dimension = length + 1
    tangent = sp.zeros(dimension)
    tangent[0, length - grade] = -4
    for row in range(grade, length + 1):
        tangent[row, row - grade] += 4
    tangent[grade, length - 2 * grade] += 4
    return tangent


def strong_blaschke_evaluation_tangent(
    crabb: sp.Matrix,
    normal: sp.Matrix,
) -> sp.Matrix:
    """Differentiate ``C^L`` in the normalized strong direction."""

    length = crabb.rows - 1
    tangent = sp.zeros(crabb.rows)
    for power in range(length):
        tangent += (
            crabb**power
            * normal
            * crabb ** (length - 1 - power)
        )
    return sp.simplify(tangent)


def predicted_strong_evaluation(
    length: int,
    grade: int,
) -> sp.Matrix:
    """Return the closed sparse formula for the strong tangent."""

    dimension = length + 1
    offset = grade - 2
    predicted = sp.zeros(dimension)
    if grade == 2:
        for index in range(dimension):
            predicted[index, index] = (
                1 if index in (0, length) else 2
            )
        return predicted

    for row in range(dimension - offset):
        predicted[row, row + offset] = (
            (grade - 1) * sp.sqrt(2)
            if row in (0, length - offset)
            else grade
        )
    return predicted


def make_record(
    length: int,
    grade: int,
) -> DualNormalAssociatedFaceRecord:
    """Derive one exact sparse mixed Hessian."""

    if length < 4 or not 2 <= grade <= length // 2:
        raise ValueError(
            "the face requires L>=4 and 2<=k<=floor(L/2)"
        )

    dimension = length + 1
    mode = length + 2 - grade
    crabb = crabb_matrix(dimension)
    normal_riesz = real_circular_normal_direction(
        dimension,
        mode,
    )
    normalized_normal = (
        normal_riesz
        - riemann_coefficient(length, grade)
        * crabb ** (mode + 1)
    )
    strong_tangent = strong_blaschke_evaluation_tangent(
        crabb,
        normalized_normal,
    )
    predicted_strong = predicted_strong_evaluation(
        length,
        grade,
    )
    strong_formula_verified = sp.simplify(
        strong_tangent - predicted_strong
    ) == sp.zeros(dimension)
    if not strong_formula_verified:
        raise AssertionError(
            "the strong Blaschke-evaluation formula changed"
        )

    reflected_tangent = prepared_reflected_tangent(
        length,
        grade,
    )
    base_blaschke = crabb**length
    reflected_gram_tangent = (
        base_blaschke.T * reflected_tangent
        + reflected_tangent.T * base_blaschke
    )
    strong_gram_tangent = (
        base_blaschke.T * strong_tangent
        + strong_tangent.T * base_blaschke
    )
    mixed_gram = (
        reflected_tangent.T * strong_tangent
        + strong_tangent.T * reflected_tangent
    )

    reflected_index = length - grade
    strong_index = grade - 2
    disjoint = reflected_index != strong_index
    if not disjoint:
        raise AssertionError("the two endpoint couplings collided")

    direct = sp.factor(mixed_gram[length, length])
    coupling = sp.factor(
        sum(
            (
                reflected_gram_tangent[length, index]
                * strong_gram_tangent[index, length]
                + strong_gram_tangent[length, index]
                * reflected_gram_tangent[index, length]
            )
            / 4
            for index in range(length)
        )
    )
    cross = sp.factor(direct + coupling)
    if direct != 0 or coupling != 0 or cross != 0:
        raise AssertionError(
            "the sparse dual Hessian selection changed"
        )

    reflected_modes = (length - grade, length + grade)
    residual_mode = 1 - grade
    product_modes = tuple(
        mode + residual_mode
        for mode in reflected_modes
    )
    residue_vanishes = (
        product_modes
        == (length - 2 * grade + 1, length + 1)
        and min(product_modes) >= 1
    )
    if not residue_vanishes:
        raise AssertionError(
            "the direct endpoint residue selection changed"
        )

    return DualNormalAssociatedFaceRecord(
        dimension=dimension,
        length=length,
        reflected_grade=grade,
        normal_mode=mode,
        reflected_coupling_index=reflected_index,
        strong_coupling_index=strong_index,
        coupling_indices_disjoint=disjoint,
        strong_evaluation_formula_verified=strong_formula_verified,
        direct_quadratic_term=str(direct),
        singular_vector_coupling_term=str(coupling),
        quadratic_hessian_cross=str(cross),
        reflected_blaschke_modes=reflected_modes,
        endpoint_residual_mode=residual_mode,
        endpoint_product_modes=product_modes,
        endpoint_residue_vanishes=residue_vanishes,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=4)
    parser.add_argument("--maximum-length", type=int, default=14)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and persist all exact sparse face audits."""

    args = parse_args()
    if (
        args.minimum_length < 4
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("lengths must satisfy 4 <= minimum <= maximum")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        ):
            for grade in range(2, length // 2 + 1):
                record = make_record(length, grade)
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")


if __name__ == "__main__":
    main()
