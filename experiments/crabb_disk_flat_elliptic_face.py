#!/usr/bin/env python3
"""Regenerate the elliptic Hessian on general disk-flat coefficients.

The older exact Hessian engine uses the simplified companion pencil
available on a phase-palindromic equality branch.  That simplification
is invalid for a general Toeplitz disk direction.  This checker starts
instead from the full coefficient-gauge formula

    S(a,c) = 2 K(a)^(-1) (H(a)R + c R^*H(a))

and only then applies the exact ellipse map and optimized rank-one
Stein Hessian.

For a raw Toeplitz coefficient at offset ``j=L-k``, the first
elliptic face is

    -64 |z_j|^2 c^(2k).

The checker also polarizes selected distinct raw offsets and verifies
that their first possible mixed face vanishes.  These exact finite
records audit L151's all-size coefficient-gauge induction.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import (
    AmplitudeMatrix,
    Series,
    amplitude_add,
    amplitude_power,
    amplitude_scale,
    direct_map_coefficients,
    diagonal_matrix,
    matrix_add,
    matrix_multiply,
    matrix_scale,
    matrix_transpose,
    matrices_equal,
    monomial,
    one,
    operator_expansion_from_coefficients,
    optimized_hessian_from_operator,
    zero_matrix,
)


DEFAULT_SINGLE_CASES = (
    (2, 1),
    (3, 1),
    (3, 2),
    (4, 1),
    (4, 2),
    (4, 3),
    (5, 1),
    (5, 2),
    (5, 3),
    (6, 1),
    (6, 2),
    (6, 3),
    (7, 1),
    (7, 2),
    (8, 1),
)
DEFAULT_CROSS_CASES = (
    (5, 1, 2),
    (5, 1, 3),
    (5, 2, 3),
)


@dataclass(frozen=True)
class RawEllipticFaceRecord:
    """One exact raw-coordinate Hessian record."""

    length: int
    dimension: int
    audit_order: int
    offsets: tuple[int, ...]
    reflected_grades: tuple[int, ...]
    target_degree: int
    first_nonzero_degree: int | None
    target_coefficient: str
    predicted_coefficient: str
    lower_coefficients_vanish: bool
    cross_face_vanishes: bool | None


def equality_coordinate_blocks(
    length: int,
    offsets: tuple[int, ...],
    order: int,
) -> tuple[
    list[list[Series]],
    list[list[Series]],
    list[list[Series]],
]:
    """Return ``K_0``, ``K_1``, and the exact numerator blocks."""

    dimension = length + 1
    base_toeplitz = zero_matrix(dimension, dimension, order)
    tangent_toeplitz = zero_matrix(dimension, dimension, order)
    shift = zero_matrix(dimension, dimension, order)
    for index in range(length):
        base_toeplitz[index][index] = one(order) / 2
        shift[index][index + 1] = one(order)
    for offset in offsets:
        for row in range(length - offset):
            tangent_toeplitz[row][row + offset] = one(order)
            tangent_toeplitz[row + offset][row] = one(order)

    shift_adjoint = matrix_transpose(shift)
    coordinate_constant = matrix_add(
        base_toeplitz,
        matrix_multiply(
            shift_adjoint,
            matrix_multiply(base_toeplitz, shift),
        ),
    )
    coordinate_tangent = matrix_add(
        tangent_toeplitz,
        matrix_multiply(
            shift_adjoint,
            matrix_multiply(tangent_toeplitz, shift),
        ),
    )
    ellipse_parameter = monomial(1, order)
    numerator_constant = matrix_add(
        matrix_multiply(base_toeplitz, shift),
        matrix_scale(
            ellipse_parameter,
            matrix_multiply(shift_adjoint, base_toeplitz),
        ),
    )
    numerator_tangent = matrix_add(
        matrix_multiply(tangent_toeplitz, shift),
        matrix_scale(
            ellipse_parameter,
            matrix_multiply(shift_adjoint, tangent_toeplitz),
        ),
    )
    return (
        coordinate_constant,
        coordinate_tangent,
        (numerator_constant, numerator_tangent),
    )


def full_disk_operator_expansion(
    length: int,
    offsets: tuple[int, ...],
    order: int,
) -> tuple[AmplitudeMatrix, list[int]]:
    """Return the correct full disk-chart operator amplitude jet."""

    dimension = length + 1
    (
        coordinate_constant,
        coordinate_tangent,
        numerator,
    ) = equality_coordinate_blocks(length, offsets, order)
    coordinate_inverse_constant = diagonal_matrix(
        [
            one(order) / coordinate_constant[index][index]
            for index in range(dimension)
        ]
    )
    coordinate_inverse_linear = matrix_scale(
        -1,
        matrix_multiply(
            coordinate_inverse_constant,
            matrix_multiply(
                coordinate_tangent,
                coordinate_inverse_constant,
            ),
        ),
    )
    coordinate_inverse_quadratic = matrix_multiply(
        coordinate_inverse_constant,
        matrix_multiply(
            coordinate_tangent,
            matrix_multiply(
                coordinate_inverse_constant,
                matrix_multiply(
                    coordinate_tangent,
                    coordinate_inverse_constant,
                ),
            ),
        ),
    )
    numerator_constant, numerator_linear = numerator
    pencil: AmplitudeMatrix = (
        matrix_scale(
            2,
            matrix_multiply(
                coordinate_inverse_constant,
                numerator_constant,
            ),
        ),
        matrix_scale(
            2,
            matrix_add(
                matrix_multiply(
                    coordinate_inverse_constant,
                    numerator_linear,
                ),
                matrix_multiply(
                    coordinate_inverse_linear,
                    numerator_constant,
                ),
            ),
        ),
        matrix_scale(
            2,
            matrix_add(
                matrix_multiply(
                    coordinate_inverse_linear,
                    numerator_linear,
                ),
                matrix_multiply(
                    coordinate_inverse_quadratic,
                    numerator_constant,
                ),
            ),
        ),
    )

    operator: AmplitudeMatrix = (
        zero_matrix(dimension, dimension, order),
        zero_matrix(dimension, dimension, order),
        zero_matrix(dimension, dimension, order),
    )
    for scalar_degree, coefficient in enumerate(
        direct_map_coefficients(order - 1, order)
    ):
        operator = amplitude_add(
            operator,
            amplitude_scale(
                coefficient,
                amplitude_power(pencil, 2 * scalar_degree + 1),
            ),
        )

    coefficient_vector = [0] * (length - 1)
    for offset in offsets:
        coefficient_vector[offset - 1] = 1
    return operator, coefficient_vector


def optimized_series(
    length: int,
    offsets: tuple[int, ...],
    order: int,
) -> Series:
    """Return the exact optimized amplitude Hessian."""

    operator, coefficients = full_disk_operator_expansion(
        length,
        offsets,
        order,
    )
    optimized, _, _ = optimized_hessian_from_operator(
        operator,
        coefficients,
        order,
    )
    return optimized


def audit_palindromic_gauge_bridge() -> None:
    """Compare the full gauge with the known companion gauge."""

    cases = (
        (3, (1, 2)),
        (4, (1, 3)),
        (4, (2,)),
        (5, (1, 4)),
        (5, (2, 3)),
    )
    order = 6
    for length, offsets in cases:
        full_operator, coefficients = full_disk_operator_expansion(
            length,
            offsets,
            order,
        )
        companion_operator, _ = operator_expansion_from_coefficients(
            length + 1,
            coefficients,
            order,
        )
        if any(
            not matrices_equal(full_operator[degree], companion_operator[degree])
            for degree in range(3)
        ):
            raise AssertionError(
                "the full and companion coefficient gauges separated"
            )


def make_single_record(length: int, grade: int) -> RawEllipticFaceRecord:
    """Validate one raw offset ``j=L-grade``."""

    offset = length - grade
    target_degree = 2 * grade
    order = target_degree + 2
    optimized = optimized_series(length, (offset,), order)
    lower_vanish = all(
        optimized.coefficient(degree) == 0
        for degree in range(target_degree)
    )
    coefficient = optimized.coefficient(target_degree)
    if (
        optimized.valuation() != target_degree
        or coefficient != -64
        or not lower_vanish
    ):
        raise AssertionError("the raw elliptic face was not -64")
    return RawEllipticFaceRecord(
        length=length,
        dimension=length + 1,
        audit_order=order,
        offsets=(offset,),
        reflected_grades=(grade,),
        target_degree=target_degree,
        first_nonzero_degree=optimized.valuation(),
        target_coefficient=str(coefficient),
        predicted_coefficient="-64",
        lower_coefficients_vanish=lower_vanish,
        cross_face_vanishes=None,
    )


def make_cross_record(
    length: int,
    first_grade: int,
    second_grade: int,
) -> RawEllipticFaceRecord:
    """Polarize two distinct raw reflected grades."""

    offsets = (
        length - first_grade,
        length - second_grade,
    )
    target_degree = first_grade + second_grade
    order = target_degree + 2
    first = optimized_series(length, (offsets[0],), order)
    second = optimized_series(length, (offsets[1],), order)
    combined = optimized_series(length, offsets, order)
    cross = combined - first - second
    cross_vanishes = all(
        cross.coefficient(degree) == 0
        for degree in range(target_degree + 1)
    )
    if not cross_vanishes:
        raise AssertionError("the raw mixed elliptic face did not vanish")
    return RawEllipticFaceRecord(
        length=length,
        dimension=length + 1,
        audit_order=order,
        offsets=offsets,
        reflected_grades=(first_grade, second_grade),
        target_degree=target_degree,
        first_nonzero_degree=(
            cross.valuation()
            if cross.valuation() < order
            else None
        ),
        target_coefficient=str(cross.coefficient(target_degree)),
        predicted_coefficient="0",
        lower_coefficients_vanish=True,
        cross_face_vanishes=cross_vanishes,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run and optionally persist the exact case grid."""

    args = parse_args()
    audit_palindromic_gauge_bridge()
    records = [
        *(
            make_single_record(length, grade)
            for length, grade in DEFAULT_SINGLE_CASES
        ),
        *(
            make_cross_record(length, first_grade, second_grade)
            for length, first_grade, second_grade in DEFAULT_CROSS_CASES
        ),
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
