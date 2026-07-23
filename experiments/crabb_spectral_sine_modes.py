#!/usr/bin/env python3
"""Audit the exact spectral sine transform of Crabb equality tangents.

The elliptic Crabb axis has the explicit DCT-I eigenvector matrix

    R = K_0^(-1/2) D U.

For the grade-j companion row

    v_j^* = e_(j+1)^* - c e_(j-1)^*,

elementary cosine subtraction gives

    (v_j^* R)_n
      = -2 sqrt(2/L) eps_n c^((j+1)/2)
        sin(j theta_n) sin(theta_n).

This checker verifies the axis diagonalization, the tangent identity,
and the discrete-sine orthogonality at high precision.  The identities
have direct all-size proofs; the finite run is only a regression.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from mpmath import mp


DEFAULT_PARAMETERS = (0.05, 0.2, 0.6)


@dataclass(frozen=True)
class SpectralSineRecord:
    dimension: int
    length: int
    ellipse_parameter: str
    axis_residual: str
    tangent_residual: str
    sine_orthogonality_residual: str


def maximum_absolute(matrix: mp.matrix) -> mp.mpf:
    """Return the largest entry modulus of a matrix."""

    return max(
        (
            abs(matrix[row, column])
            for row in range(matrix.rows)
            for column in range(matrix.cols)
        ),
        default=mp.mpf("0"),
    )


def coordinate_axis(length: int, parameter: mp.mpf) -> mp.matrix:
    """Return ``C+cJCJ`` in the companion coordinates."""

    dimension = length + 1
    result = mp.matrix(dimension)
    result[0, 1] = 2
    for column in range(2, dimension):
        result[column - 1, column] = 1
    result[length, length - 1] += 2 * parameter
    for row in range(1, length):
        result[row, row - 1] += parameter
    return result


def dct_data(
    length: int,
    parameter: mp.mpf,
) -> tuple[mp.matrix, mp.matrix, list[mp.mpf], list[mp.mpf]]:
    """Return the right eigenvectors, nodes, endpoint factors, and angles."""

    dimension = length + 1
    endpoint_factors = [mp.mpf(1) for _ in range(dimension)]
    endpoint_factors[0] = endpoint_factors[-1] = 1 / mp.sqrt(2)
    angles = [
        mp.pi * index / length
        for index in range(dimension)
    ]
    dct = mp.matrix(dimension)
    scale = mp.sqrt(mp.mpf(2) / length)
    for row in range(dimension):
        for column in range(dimension):
            dct[row, column] = (
                scale
                * endpoint_factors[row]
                * endpoint_factors[column]
                * mp.cos(row * angles[column])
            )

    right_eigenvectors = mp.matrix(dimension)
    for row in range(dimension):
        coordinate_scale = (
            parameter ** (mp.mpf(row) / 2)
            / endpoint_factors[row]
        )
        for column in range(dimension):
            right_eigenvectors[row, column] = (
                coordinate_scale * dct[row, column]
            )

    nodes = mp.matrix(dimension)
    for index, angle in enumerate(angles):
        nodes[index, index] = 2 * mp.sqrt(parameter) * mp.cos(angle)
    return right_eigenvectors, nodes, endpoint_factors, angles


def make_record(
    dimension: int,
    parameter: mp.mpf,
    tolerance: mp.mpf,
) -> SpectralSineRecord:
    """Construct and validate one high-precision regression record."""

    length = dimension - 1
    axis = coordinate_axis(length, parameter)
    right, nodes, endpoint_factors, angles = dct_data(
        length,
        parameter,
    )
    axis_residual = maximum_absolute(axis * right - right * nodes)

    tangent_residual = mp.mpf("0")
    scale = mp.sqrt(mp.mpf(2) / length)
    for grade in range(1, length):
        for node in range(dimension):
            actual = (
                right[grade + 1, node]
                - parameter * right[grade - 1, node]
            )
            predicted = (
                -2
                * scale
                * endpoint_factors[node]
                * parameter ** (mp.mpf(grade + 1) / 2)
                * mp.sin(grade * angles[node])
                * mp.sin(angles[node])
            )
            tangent_residual = max(
                tangent_residual,
                abs(actual - predicted),
            )

    sine_orthogonality_residual = mp.mpf("0")
    for first_grade in range(1, length):
        for second_grade in range(1, length):
            inner_product = mp.fsum(
                mp.sin(first_grade * angles[node])
                * mp.sin(second_grade * angles[node])
                for node in range(1, length)
            )
            expected = (
                mp.mpf(length) / 2
                if first_grade == second_grade
                else mp.mpf("0")
            )
            sine_orthogonality_residual = max(
                sine_orthogonality_residual,
                abs(inner_product - expected),
            )

    maximum_residual = max(
        axis_residual,
        tangent_residual,
        sine_orthogonality_residual,
    )
    if maximum_residual > tolerance:
        raise AssertionError(
            "the spectral sine factorization failed: "
            f"{mp.nstr(maximum_residual, 8)}"
        )

    return SpectralSineRecord(
        dimension=dimension,
        length=length,
        ellipse_parameter=mp.nstr(parameter, 20),
        axis_residual=mp.nstr(axis_residual, 8),
        tangent_residual=mp.nstr(tangent_residual, 8),
        sine_orthogonality_residual=mp.nstr(
            sine_orthogonality_residual,
            8,
        ),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=20)
    parser.add_argument(
        "--ellipse-parameters",
        type=str,
        nargs="+",
        default=tuple(str(value) for value in DEFAULT_PARAMETERS),
    )
    parser.add_argument("--precision", type=int, default=80)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic high-precision grid."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    if args.precision < 40:
        raise ValueError("use at least 40 decimal digits")

    mp.dps = args.precision
    parameters = tuple(mp.mpf(value) for value in args.ellipse_parameters)
    if any(not 0 < value < 1 for value in parameters):
        raise ValueError("ellipse parameters must lie in (0,1)")
    tolerance = mp.power(10, -(args.precision // 2))

    records = [
        make_record(dimension, parameter, tolerance)
        for dimension in range(args.minimum_size, args.maximum_size + 1)
        for parameter in parameters
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
