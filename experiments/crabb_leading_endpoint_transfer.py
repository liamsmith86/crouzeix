#!/usr/bin/env python3
"""Regenerate the leading circular-normal endpoint transfer.

For grade ``k`` and eligible normal mode ``m=L+2-k``, the first
Riemann-normalized variation at a Crabb block changes the endpoint
resolvent by

    D log(e0^T (zI-C)^(-1) eL)
        = k z^(-m) + beta_(L,k) z^m.

For every ``k>=2``, subtracting L163's matching symmetrized defect
response leaves the logarithmic-inner tangent
``k(z^(-m)-z^m)``.  The defect endpoint factor is ``1`` for ``k=2``
and ``sqrt(2)`` for ``k>=3``.
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
from crabb_normal_defect_column import predicted_coefficient
from general_crabb_weighted_series import crabb_matrix


@dataclass(frozen=True)
class EndpointTransferRecord:
    """One exact leading endpoint-transfer identity."""

    dimension: int
    length: int
    reflected_grade: int
    normal_mode: int
    riemann_coefficient: str
    characteristic_derivative: str
    characteristic_inner_derivative: str
    schur_iterate_derivative: str
    schur_grade_discriminator_verified: bool
    endpoint_log_derivative: str
    positive_mode_coefficient: str
    symmetrized_defect_log_derivative: str | None
    defect_coefficient_relation_verified: bool | None
    inner_tangent_verified: bool | None


def riemann_coefficient(length: int, grade: int) -> sp.Expr:
    """Return the first Schwarz coefficient of the normal support mode."""

    if grade <= 2:
        return sp.Integer(0)
    if grade == 3:
        return sp.Rational(9, 4 * length)
    return sp.Rational(2 * (grade - 2), length)


def positive_mode_coefficient(length: int, grade: int) -> sp.Expr:
    """Return the positive endpoint-transfer coefficient."""

    if grade <= 2:
        return sp.Integer(0)
    if grade == 3:
        return sp.Rational(length - 9, 4 * length)
    return (
        grade
        - 3
        - sp.Rational(2 * (grade - 2) ** 2, length)
    )


def make_record(length: int, grade: int) -> EndpointTransferRecord:
    """Derive one transfer identity without a weighted-series expansion."""

    if not 1 <= grade <= length // 2:
        raise ValueError("the grade must satisfy 1 <= k <= floor(L/2)")

    dimension = length + 1
    mode = length + 2 - grade
    spectral_parameter = sp.symbols("spectral_parameter")
    crabb = crabb_matrix(dimension)
    normal = real_circular_normal_direction(dimension, mode)
    coefficient = riemann_coefficient(length, grade)
    normalized = normal - coefficient * crabb ** (mode + 1)

    resolvent = sum(
        (
            crabb**power
            / spectral_parameter ** (power + 1)
            for power in range(dimension)
        ),
        sp.zeros(dimension),
    )
    first_endpoint = sp.eye(dimension)[:, 0]
    last_endpoint = sp.eye(dimension)[:, length]
    base_transfer = (
        first_endpoint.T * resolvent * last_endpoint
    )[0]
    transfer_derivative = (
        first_endpoint.T
        * resolvent
        * normalized
        * resolvent
        * last_endpoint
    )[0]
    logarithmic_derivative = sp.factor(
        transfer_derivative / base_transfer
    )

    characteristic_derivative = sp.factor(
        -sp.trace(
            (
                spectral_parameter * sp.eye(dimension)
                - crabb
            ).adjugate()
            * normalized
        )
    )
    predicted_characteristic = (
        -grade * spectral_parameter ** (grade - 1)
    )
    if sp.simplify(
        characteristic_derivative - predicted_characteristic
    ) != 0:
        raise AssertionError(
            "the leading characteristic derivative changed"
        )

    base_inner = spectral_parameter**dimension
    denominator_derivative = sp.factor(
        spectral_parameter**dimension
        * characteristic_derivative.subs(
            spectral_parameter,
            1 / spectral_parameter,
        )
    )
    characteristic_inner_derivative = sp.factor(
        characteristic_derivative
        - base_inner * denominator_derivative
    )
    feedthrough_derivative = sp.expand(
        characteristic_inner_derivative
    ).coeff(spectral_parameter, 0)
    schur_iterate_derivative = sp.factor(
        (
            characteristic_inner_derivative
            - feedthrough_derivative
        )
        / spectral_parameter
        + feedthrough_derivative
        * spectral_parameter ** (2 * dimension - 1)
    )
    predicted_schur_derivative = (
        sp.Integer(0)
        if grade == 1
        else grade
        * (
            -spectral_parameter ** (grade - 2)
            + spectral_parameter ** (2 * dimension - grade)
        )
    )
    schur_discriminator = sp.simplify(
        schur_iterate_derivative - predicted_schur_derivative
    ) == 0
    if not schur_discriminator:
        raise AssertionError("the Schur grade discriminator changed")

    beta = positive_mode_coefficient(length, grade)
    predicted_transfer = (
        grade * spectral_parameter ** (-mode)
        + beta * spectral_parameter**mode
    )
    if sp.simplify(
        logarithmic_derivative - predicted_transfer
    ) != 0:
        raise AssertionError(
            "the endpoint transfer formula changed"
        )

    symmetrized_defect: sp.Expr | None = None
    relation: bool | None = None
    inner_tangent_verified: bool | None = None
    if grade >= 2:
        defect_coordinate = sp.eye(dimension)[:, mode]
        one_sided_defect = sp.factor(
            (
                defect_coordinate.T
                * resolvent
                * last_endpoint
            )[0]
            / base_transfer
        )
        symmetrized_defect = sp.factor(2 * one_sided_defect)
        gamma = predicted_coefficient(length, grade)
        relation = sp.simplify(
            beta * spectral_parameter**mode
            + grade * spectral_parameter**mode
            - gamma * symmetrized_defect
        ) == 0
        if not relation:
            raise AssertionError(
                "the endpoint/defect coefficient relation changed"
            )
        inner_tangent = sp.factor(
            logarithmic_derivative
            - gamma * symmetrized_defect
        )
        predicted_inner_tangent = grade * (
            spectral_parameter ** (-mode)
            - spectral_parameter**mode
        )
        inner_tangent_verified = sp.simplify(
            inner_tangent - predicted_inner_tangent
        ) == 0
        if not inner_tangent_verified:
            raise AssertionError("the leading inner tangent changed")

    return EndpointTransferRecord(
        dimension=dimension,
        length=length,
        reflected_grade=grade,
        normal_mode=mode,
        riemann_coefficient=str(coefficient),
        characteristic_derivative=str(characteristic_derivative),
        characteristic_inner_derivative=str(
            characteristic_inner_derivative
        ),
        schur_iterate_derivative=str(schur_iterate_derivative),
        schur_grade_discriminator_verified=schur_discriminator,
        endpoint_log_derivative=str(logarithmic_derivative),
        positive_mode_coefficient=str(beta),
        symmetrized_defect_log_derivative=(
            str(symmetrized_defect)
            if symmetrized_defect is not None
            else None
        ),
        defect_coefficient_relation_verified=relation,
        inner_tangent_verified=inner_tangent_verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=4)
    parser.add_argument("--maximum-length", type=int, default=14)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact all-size shadows."""

    args = parse_args()
    if (
        args.minimum_length < 2
        or args.maximum_length < args.minimum_length
    ):
        raise ValueError("length range must satisfy 2 <= minimum <= maximum")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        ):
            for grade in range(1, length // 2 + 1):
                record = make_record(length, grade)
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")


if __name__ == "__main__":
    main()
