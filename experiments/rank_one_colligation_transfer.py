#!/usr/bin/env python3
"""Regenerate the rank-one Stein/colligation transfer identity.

An orthogonal colligation

    U = [[S, b], [a^T, d]]

gives a rank-one Stein pair after an arbitrary positive diagonal
similarity.  This checker verifies exactly that its scalar transfer is
the finite inner determinant quotient

    d + z q^T (I-zT)^(-1) r
      = omega det(zI-T^T) / det(I-zT).

The calculation uses unrelated rational colligations and similarities
in several dimensions; it is an audit of the general algebra, not a
fit to Crabb data.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


RATIONAL_ROTATIONS = (
    (sp.Rational(3, 5), sp.Rational(4, 5)),
    (sp.Rational(5, 13), sp.Rational(12, 13)),
    (sp.Rational(8, 17), sp.Rational(15, 17)),
)


@dataclass(frozen=True)
class ColligationTransferRecord:
    """One exact rational colligation audit."""

    dimension: int
    colligation_determinant: str
    determinant_phase: str
    right_stein_verified: bool
    left_stein_verified: bool
    transfer_determinant_verified: bool
    outside_resolvent_verified: bool


def adjacent_rotation(
    dimension: int,
    index: int,
    cosine: sp.Rational,
    sine: sp.Rational,
) -> sp.Matrix:
    """Return one exact adjacent Givens rotation."""

    rotation = sp.eye(dimension)
    rotation[index, index] = cosine
    rotation[index, index + 1] = sine
    rotation[index + 1, index] = -sine
    rotation[index + 1, index + 1] = cosine
    return rotation


def rational_colligation(dimension: int) -> sp.Matrix:
    """Return a deterministic dense rational orthogonal colligation."""

    size = dimension + 1
    colligation = sp.eye(size)
    for step in range(2 * size + 1):
        cosine, sine = RATIONAL_ROTATIONS[
            step % len(RATIONAL_ROTATIONS)
        ]
        rotation = adjacent_rotation(
            size,
            step % dimension,
            cosine,
            sine,
        )
        colligation = rotation * colligation
    if sp.simplify(
        colligation.T * colligation - sp.eye(size)
    ) != sp.zeros(size):
        raise AssertionError("the rational colligation is not orthogonal")
    return colligation


def make_record(dimension: int) -> ColligationTransferRecord:
    """Verify the transfer theorem in one dimension."""

    if dimension < 2:
        raise ValueError("the state dimension must be at least two")

    colligation = rational_colligation(dimension)
    state = colligation[:dimension, :dimension]
    left_defect = colligation[:dimension, dimension]
    right_defect = colligation[dimension, :dimension].T
    feedthrough = colligation[dimension, dimension]

    similarity = sp.diag(
        *[sp.Integer(index + 1) for index in range(dimension)]
    )
    metric = similarity**2
    operator = similarity.inv() * state * similarity
    defect = similarity * right_defect
    inverse_defect = similarity.inv() * left_defect

    right_stein = sp.simplify(
        metric
        - operator.T * metric * operator
        - defect * defect.T
    ) == sp.zeros(dimension)
    left_stein = sp.simplify(
        metric.inv()
        - operator * metric.inv() * operator.T
        - inverse_defect * inverse_defect.T
    ) == sp.zeros(dimension)
    if not right_stein or not left_stein:
        raise AssertionError("a rank-one Stein identity changed")

    spectral_parameter = sp.symbols("spectral_parameter")
    phase = sp.factor(
        (-1) ** dimension * colligation.det()
    )
    transfer = sp.factor(
        feedthrough
        + spectral_parameter
        * (
            defect.T
            * (
                sp.eye(dimension)
                - spectral_parameter * operator
            ).inv()
            * inverse_defect
        )[0]
    )
    determinant_transfer = sp.factor(
        phase
        * (
            spectral_parameter * sp.eye(dimension)
            - operator.T
        ).det()
        / (
            sp.eye(dimension)
            - spectral_parameter * operator
        ).det()
    )
    transfer_verified = sp.factor(
        transfer - determinant_transfer
    ) == 0
    if not transfer_verified:
        raise AssertionError("the determinant transfer identity changed")

    outside_parameter = sp.symbols("outside_parameter")
    outside_transfer = sp.factor(
        feedthrough
        + (
            defect.T
            * (
                outside_parameter * sp.eye(dimension)
                - operator
            ).inv()
            * inverse_defect
        )[0]
    )
    reciprocal_transfer = sp.factor(
        transfer.subs(
            spectral_parameter,
            1 / outside_parameter,
        )
    )
    outside_verified = sp.factor(
        outside_transfer - reciprocal_transfer
    ) == 0
    if not outside_verified:
        raise AssertionError("the outside resolvent identity changed")

    return ColligationTransferRecord(
        dimension=dimension,
        colligation_determinant=str(
            sp.factor(colligation.det())
        ),
        determinant_phase=str(phase),
        right_stein_verified=right_stein,
        left_stein_verified=left_stein,
        transfer_determinant_verified=transfer_verified,
        outside_resolvent_verified=outside_verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=2)
    parser.add_argument("--maximum-size", type=int, default=7)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact rational audits."""

    args = parse_args()
    if (
        args.minimum_size < 2
        or args.maximum_size < args.minimum_size
    ):
        raise ValueError("sizes must satisfy 2 <= minimum <= maximum")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for dimension in range(
            args.minimum_size,
            args.maximum_size + 1,
        ):
            record = make_record(dimension)
            line = json.dumps(asdict(record), sort_keys=True)
            print(line, flush=True)
            output.write(line + "\n")


if __name__ == "__main__":
    main()
