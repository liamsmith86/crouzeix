#!/usr/bin/env python3
"""Regenerate the raw disk-flat Faber endpoint identity.

For one Hermitian Toeplitz coefficient at offset ``j`` in a Crabb
chain of length ``L``, put ``k=L-j``.  This checker differentiates the
full coefficient-gauge pencil, rather than the phase-palindromic
companion shortcut, and verifies

    d det(xi I-A) = 2 conjugate(z) xi^(k+1),

and

    e_0^* (DP_L[Y] + 2 conjugate(z) P_k) = 4 conjugate(z) e_k^*,
    e_L^* (DP_L[Y] + 2 conjugate(z) P_k)
        = 4 conjugate(z) c^k e_j^*.

The two phases ``z=1`` and ``z=i`` audit the real and imaginary
directions.  The identities are polynomial in the ellipse parameter
``c`` and are checked exactly over the Gaussian rationals.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class RawFaberEndpointRecord:
    """One exact offset/phase endpoint audit."""

    length: int
    dimension: int
    offset: int
    reflected_grade: int
    phase: str
    characteristic_identity: bool
    top_endpoint_identity: bool
    bottom_endpoint_identity: bool


def shift(dimension: int) -> sp.Matrix:
    """Return the nilpotent forward shift."""

    result = sp.zeros(dimension)
    for row in range(dimension - 1):
        result[row, row + 1] = 1
    return result


def toeplitz_tangent(
    length: int,
    offset: int,
    phase: sp.Expr,
) -> sp.Matrix:
    """Return one Hermitian Toeplitz coefficient tangent."""

    result = sp.zeros(length + 1)
    for row in range(length - offset):
        result[row, row + offset] = phase
        result[row + offset, row] = sp.conjugate(phase)
    return result


def dickson_endpoint_audit(
    length: int,
    offset: int,
    phase: sp.Expr,
) -> RawFaberEndpointRecord:
    """Differentiate the coefficient gauge and verify both endpoints."""

    dimension = length + 1
    grade = length - offset
    c, amplitude, xi = sp.symbols("c amplitude xi")
    reversal = shift(dimension)
    base_toeplitz = sp.zeros(dimension)
    for index in range(length):
        base_toeplitz[index, index] = sp.Rational(1, 2)
    tangent_toeplitz = toeplitz_tangent(length, offset, phase)

    coordinate_constant = (
        base_toeplitz
        + reversal.T * base_toeplitz * reversal
    )
    coordinate_tangent = (
        tangent_toeplitz
        + reversal.T * tangent_toeplitz * reversal
    )
    coordinate_inverse = coordinate_constant.inv()
    numerator_constant = (
        base_toeplitz * reversal
        + c * reversal.T * base_toeplitz
    )
    numerator_tangent = (
        tangent_toeplitz * reversal
        + c * reversal.T * tangent_toeplitz
    )
    axis = sp.expand(2 * coordinate_inverse * numerator_constant)
    direction = sp.expand(
        2
        * (
            coordinate_inverse * numerator_tangent
            - coordinate_inverse
            * coordinate_tangent
            * coordinate_inverse
            * numerator_constant
        )
    )

    disk_axis = axis.subs(c, 0)
    disk_direction = direction.subs(c, 0)
    disk_pencil = disk_axis + amplitude * disk_direction
    characteristic_derivative = sp.diff(
        (xi * sp.eye(dimension) - disk_pencil).det(),
        amplitude,
    ).subs(amplitude, 0)
    expected_characteristic = (
        2 * sp.conjugate(phase) * xi ** (grade + 1)
    )
    characteristic_ok = (
        sp.expand(
            characteristic_derivative - expected_characteristic
        )
        == 0
    )

    polynomials = [2 * sp.eye(dimension), axis]
    derivatives = [sp.zeros(dimension), direction]
    for degree in range(2, length + 1):
        polynomials.append(
            sp.expand(
                axis * polynomials[-1]
                - c * polynomials[-2]
            )
        )
        derivatives.append(
            sp.expand(
                direction * polynomials[-2]
                + axis * derivatives[-1]
                - c * derivatives[-2]
            )
        )

    faber_derivative = sp.expand(
        derivatives[length]
        + 2 * sp.conjugate(phase) * polynomials[grade]
    )
    expected_top = sp.zeros(1, dimension)
    expected_top[0, grade] = 4 * sp.conjugate(phase)
    expected_bottom = sp.zeros(1, dimension)
    expected_bottom[0, offset] = (
        4 * sp.conjugate(phase) * c**grade
    )
    top_ok = all(
        sp.expand(
            faber_derivative[0, column]
            - expected_top[0, column]
        )
        == 0
        for column in range(dimension)
    )
    bottom_ok = all(
        sp.expand(
            faber_derivative[length, column]
            - expected_bottom[0, column]
        )
        == 0
        for column in range(dimension)
    )
    if not (characteristic_ok and top_ok and bottom_ok):
        raise AssertionError("the raw Faber endpoint identity failed")

    return RawFaberEndpointRecord(
        length=length,
        dimension=dimension,
        offset=offset,
        reflected_grade=grade,
        phase=str(phase),
        characteristic_identity=characteristic_ok,
        top_endpoint_identity=top_ok,
        bottom_endpoint_identity=bottom_ok,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-length", type=int, default=14)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact offset/phase grid."""

    args = parse_args()
    records = [
        dickson_endpoint_audit(length, offset, phase)
        for length in range(2, args.maximum_length + 1)
        for offset in range(1, length)
        for phase in (sp.Integer(1), sp.I)
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
