#!/usr/bin/env python3
"""Exact audit of the sharp positive determinant corner square.

At ``X=R=Y=1`` and with the exact envelope identity
``gamma_width=1-s-gamma_lower``, both positive determinant charts reduce to
the same polynomial ``c^7 F(A,B,c,k)^2``, where ``A=1-a`` and ``B=1-b``.
This script regenerates that identity directly from every integer chart
record without symbolic factorization or floating-point arithmetic.
"""

from __future__ import annotations

from collections import defaultdict
from math import comb
from typing import TypeAlias

from slice_projective_core import CoreRecords, load_records
from slice_projective_interval_certificate import final_chart_tables


CornerMonomial: TypeAlias = tuple[int, int, int, int, int, int]
CornerPolynomial: TypeAlias = dict[CornerMonomial, int]
RootMonomial: TypeAlias = tuple[int, int, int, int]
RootPolynomial: TypeAlias = dict[RootMonomial, int]


def actual_corner(records: dict[tuple[int, ...], int]) -> CornerPolynomial:
    """Substitute the three unit faces and collect in ``A,B,c,k,s,lower``."""

    output: defaultdict[CornerMonomial, int] = defaultdict(int)
    for monomial, coefficient in records.items():
        (
            _,
            _,
            _,
            a_power,
            b_power,
            c_power,
            k_power,
            s_power,
            lower_power,
            width_power,
        ) = monomial
        for a_deficit_power in range(a_power + 1):
            a_coefficient = comb(a_power, a_deficit_power) * (-1) ** a_deficit_power
            for b_deficit_power in range(b_power + 1):
                b_coefficient = comb(b_power, b_deficit_power) * (-1) ** b_deficit_power
                for extra_s_power in range(width_power + 1):
                    for extra_lower_power in range(width_power - extra_s_power + 1):
                        multinomial = comb(width_power, extra_s_power) * comb(
                            width_power - extra_s_power,
                            extra_lower_power,
                        )
                        key = (
                            a_deficit_power,
                            b_deficit_power,
                            c_power,
                            k_power,
                            s_power + extra_s_power,
                            lower_power + extra_lower_power,
                        )
                        output[key] += (
                            coefficient
                            * a_coefficient
                            * b_coefficient
                            * multinomial
                            * (-1) ** (extra_s_power + extra_lower_power)
                        )
    return {
        monomial: coefficient for monomial, coefficient in output.items() if coefficient
    }


def expected_root() -> RootPolynomial:
    """Return the 30-term square root ``F(A,B,c,k)``."""

    output: RootPolynomial = {}

    def add(
        a_power: int,
        b_power: int,
        terms: tuple[tuple[int, int, int], ...],
    ) -> None:
        for c_power, k_power, coefficient in terms:
            output[a_power, b_power, c_power, k_power] = coefficient

    add(2, 2, ((2, 1, 4), (1, 2, -16), (1, 0, -1), (0, 1, 4)))
    mixed_cubic = ((2, 1, -8), (1, 2, 32), (1, 0, 2), (0, 1, -8))
    add(2, 1, mixed_cubic)
    add(1, 2, mixed_cubic)
    pure_quadratic = ((1, 2, -12), (1, 0, 3))
    add(2, 0, pure_quadratic)
    add(0, 2, pure_quadratic)
    add(
        1,
        1,
        ((2, 1, 16), (1, 2, -64), (1, 1, 18), (1, 0, -4), (0, 1, 16)),
    )
    pure_linear = ((1, 2, 24), (1, 1, -18), (1, 0, -6))
    add(1, 0, pure_linear)
    add(0, 1, pure_linear)
    add(0, 0, ((1, 2, -9), (1, 1, 18), (1, 0, -9)))
    return output


def expected_corner() -> CornerPolynomial:
    """Square the explicit root and multiply by ``c^7`` exactly."""

    root = expected_root()
    output: defaultdict[CornerMonomial, int] = defaultdict(int)
    for left_monomial, left_coefficient in root.items():
        for right_monomial, right_coefficient in root.items():
            a_left, b_left, c_left, k_left = left_monomial
            a_right, b_right, c_right, k_right = right_monomial
            key = (
                a_left + a_right,
                b_left + b_right,
                c_left + c_right + 7,
                k_left + k_right,
                0,
                0,
            )
            output[key] += left_coefficient * right_coefficient
    return {
        monomial: coefficient for monomial, coefficient in output.items() if coefficient
    }


def audit_corner_tables(records: CoreRecords | None = None) -> None:
    """Regenerate and verify the corner square in both determinant charts."""

    expected = expected_corner()
    source = records if records is not None else load_records()
    tables = final_chart_tables(source, 1)[:2]
    actual = [actual_corner(table.records) for table in tables]
    for table, polynomial in zip(tables, actual):
        if polynomial != expected:
            raise AssertionError(f"{table.label}: positive corner square failed")
        print(
            f"PASS {table.label}: corner_terms={len(polynomial)}, "
            f"square_root_terms={len(expected_root())}"
        )
    if actual[0] != actual[1]:
        raise AssertionError("the two positive corner charts do not agree")


def main() -> None:
    audit_corner_tables()


if __name__ == "__main__":
    main()
