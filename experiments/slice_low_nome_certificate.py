#!/usr/bin/env python3
"""Exact low-nome normalization for the elliptic transfer-core charts.

The ordinary projective certificate becomes cancellation-sensitive as the
ellipse nome ``c`` tends to zero.  This companion rewrites

``k=c(4+c^2 K2)``, ``s=1+c^2 S2``,
``gamma_lower=c^2(4+c^2 L2)``, and
``gamma_width=c^4(32+c^2 D2)``

before collecting coefficients.  The determinant then has exact order nine
in ``c`` and every final minor chart has exact order two.  The scalar
deviations below are evaluated without removable divisions at ``c=0``.

The current charts expose the remaining sharp ridge but do not yet certify a
full low-nome interval; a failed command is diagnostic, not a counterexample.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Iterator
from collections import defaultdict
from fractions import Fraction
import gc
from math import comb, inf, nextafter
import time

from flint import arb, ctx, fmpq

from slice_projective_core import load_records
from slice_projective_interval_certificate import (
    ArbTaylorJet,
    Box,
    ChartTable,
    FactorJets,
    IntervalTensor,
    arb_positive_series_tail,
    arb_map_to_tensor,
    as_arb_jet,
    certify,
    collapse_bernstein_axis,
    factor_jets_from_data,
    final_chart_tables,
    interval_tensor_for_chart,
    parse_fraction,
    power_to_bernstein_axis,
    power_to_bernstein_tensor,
)


LOW_CONSTANTS = (4, 1, 4, 32)
LOW_BASE_POWERS = (1, 1, 0, 2, 4)
LOW_DEVIATION_SERIES = (
    (-16, 0, 56, 0, -160, 0, 404, 0, -944, 0, 2072, 0, -4320, 0, 8648, 0),
    (-4, 0, 12, 0, -32, 0, 76, 0, -168, 0, 352, 0, -704, 0, 1356, 0),
    (
        -44,
        0,
        304,
        0,
        -1644,
        0,
        7544,
        0,
        -30_672,
        0,
        113_440,
        0,
        -388_460,
        0,
        1_247_316,
        0,
    ),
    (
        -272,
        0,
        1568,
        0,
        -7376,
        0,
        30_320,
        0,
        -112_736,
        0,
        387_104,
        0,
        -1_244_784,
        0,
        3_786_480,
        0,
    ),
)

SparseCoefficient = fmpq | arb
SparseMap = dict[tuple[int, ...], SparseCoefficient]
ORIENTATION_CORNER_SCALE = Fraction(1, 2)
ZERO_PARAMETER_CORNER_SCALE = Fraction(1, 4)
ZERO_PARAMETER_NOME_ARM_CAP = Fraction(1, 2)
POSITIVE_ZERO_PARAMETER_MAIN_ARM_CAP = Fraction(1)
NEGATIVE_ZERO_PARAMETER_MAIN_ARM_CAP = Fraction(1, 2)
ZERO_PARAMETER_MAIN_ARM_LOWER = Fraction(1, 32)
MINOR_CORNER_SCALE = Fraction(1, 4)
SECOND_DETERMINANT_RADIAL_MODELS = tuple(
    (Fraction(index, 10), Fraction(1, 10), degree)
    for index, degree in enumerate((8, 8, 8, 8, 10, 12, 14, 16, 20, 24), start=1)
)


def low_nome_table(table: ChartTable) -> ChartTable:
    """Apply the exact removable-singularity normalization to one chart."""

    expected_order = 9 if table.label.startswith("det-") else 2
    output: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for monomial, coefficient in table.records.items():
        parameters = monomial[:5]
        factor_exponents = monomial[5:]
        base_power = sum(
            exponent * power
            for exponent, power in zip(factor_exponents, LOW_BASE_POWERS)
        )
        states = {(base_power - expected_order, 0, 0, 0, 0): coefficient}
        for deviation_axis, (exponent, constant) in enumerate(
            zip(factor_exponents[1:], LOW_CONSTANTS), start=1
        ):
            if not exponent:
                continue
            expanded: defaultdict[tuple[int, ...], int] = defaultdict(int)
            for state, value in states.items():
                for deviation_power in range(exponent + 1):
                    key = list(state)
                    key[0] += 2 * deviation_power
                    key[deviation_axis] += deviation_power
                    expanded[tuple(key)] += (
                        value
                        * comb(exponent, deviation_power)
                        * constant ** (exponent - deviation_power)
                    )
            states = expanded
        for factor_key, value in states.items():
            output[parameters + factor_key] += value

    for monomial in [key for key, value in output.items() if not value]:
        del output[monomial]
    actual_order = min(monomial[5] for monomial in output)
    if actual_order != 0:
        raise AssertionError(
            f"{table.label}: normalized c-order is {actual_order}, expected zero"
        )
    return ChartTable(f"{table.label}-low-c{expected_order}", dict(output))


def multiply_series(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Multiply two equally truncated exact integer series."""

    if len(left) != len(right):
        raise ValueError("series truncations do not match")
    return tuple(
        sum(left[index] * right[order - index] for index in range(order + 1))
        for order in range(len(left))
    )


def power_series(series: tuple[int, ...], exponent: int) -> tuple[int, ...]:
    """Raise a truncated exact series to a nonnegative integer power."""

    if exponent < 0:
        raise ValueError(exponent)
    output = (1,) + (0,) * (len(series) - 1)
    base = series
    remaining = exponent
    while remaining:
        if remaining & 1:
            output = multiply_series(output, base)
        remaining //= 2
        if remaining:
            base = multiply_series(base, base)
    return output


def low_parameter_taylor_maps(
    table: ChartTable, degree: int = 14
) -> list[dict[tuple[int, ...], int]]:
    """Collect exact Maclaurin coefficient maps of a normalized low-nome chart."""

    if not 0 <= degree < len(LOW_DEVIATION_SERIES[0]):
        raise ValueError("requested degree exceeds the audited deviation series")
    maps: list[defaultdict[tuple[int, ...], int]] = [
        defaultdict(int) for _ in range(degree + 1)
    ]
    factor_cache: dict[tuple[int, ...], tuple[int, ...]] = {}
    for monomial, coefficient in table.records.items():
        factor_key = monomial[5:]
        factor_series = factor_cache.get(factor_key)
        if factor_series is None:
            c_power, *deviation_powers = factor_key
            if c_power < 0:
                raise AssertionError("low-nome normalization retained a negative power")
            factor_series = (1,) + (0,) * degree
            for exponent, series in zip(deviation_powers, LOW_DEVIATION_SERIES):
                factor_series = multiply_series(
                    factor_series, power_series(series[: degree + 1], exponent)
                )
            factor_series = (0,) * min(c_power, degree + 1) + factor_series
            factor_series = factor_series[: degree + 1]
            factor_series += (0,) * (degree + 1 - len(factor_series))
            factor_cache[factor_key] = factor_series
        for order, value in enumerate(factor_series):
            if value:
                maps[order][monomial[:5]] += coefficient * value
    return [{key: value for key, value in part.items() if value} for part in maps]


def scaled_corner_coefficient(
    maps: list[dict[tuple[int, ...], int]],
    degree: int,
    *,
    second_scale: bool,
) -> dict[tuple[int, ...], int]:
    """Return one exact exceptional coefficient at the sharp determinant corner."""

    output: defaultdict[tuple[int, ...], int] = defaultdict(int)
    ua_weight = 2 if second_scale else 1
    for nome_order, coefficient_map in enumerate(maps[: degree + 1]):
        for (
            u_power,
            v_power,
            y_power,
            a_power,
            b_power,
        ), coefficient in coefficient_map.items():
            for centered_v in range(v_power + 1):
                for centered_a in range(a_power + 1):
                    total = (
                        nome_order
                        + ua_weight * u_power
                        + centered_v
                        + ua_weight * centered_a
                        + b_power
                    )
                    if total != degree:
                        continue
                    key = (u_power, centered_v, y_power, centered_a, b_power)
                    output[key] += (
                        coefficient
                        * comb(v_power, centered_v)
                        * (-1) ** centered_v
                        * comb(a_power, centered_a)
                        * (-1) ** centered_a
                    )
    return {key: value for key, value in output.items() if value}


def audit_determinant_ridge() -> None:
    """Verify all low orders and the determinant exceptional hierarchy."""

    expected_first = {
        (1, 0, 0, 0, 0): 4608,
        (0, 0, 0, 1, 0): 1536,
    }
    core = load_records()
    for sign in (1, -1):
        expected_second = {
            (0, 0, 0, 0, 0): 13_824,
            (1, 0, 0, 0, 0): 4_608,
            (0, 2, 0, 0, 0): 576,
            (0, 0, 0, 1, 0): 1_536,
            (0, 0, 0, 0, 1): -13_824 * sign,
            (0, 0, 0, 0, 2): 3_456,
        }
        for source in final_chart_tables(core, sign)[:2]:
            table = low_nome_table(source)
            maps = low_parameter_taylor_maps(table)
            first = scaled_corner_coefficient(maps, 1, second_scale=False)
            second = scaled_corner_coefficient(maps, 2, second_scale=True)
            if first != expected_first or second != expected_second:
                raise AssertionError(f"unexpected ridge hierarchy in {source.label}")
            if sign == 1:
                ridge = defaultdict(int)
                for nome_order, coefficient_map in enumerate(maps):
                    for (
                        u_power,
                        _v,
                        y_power,
                        _a,
                        b_power,
                    ), coefficient in coefficient_map.items():
                        if u_power == 0 and nome_order + b_power <= 6:
                            ridge[(nome_order + b_power, y_power)] += (
                                coefficient * 2**b_power
                            )
                ridge = {key: value for key, value in ridge.items() if value}
                if ridge != {(4, 0): 46_080, (6, 0): -1_622_016}:
                    raise AssertionError(
                        f"unexpected final ridge series in {source.label}"
                    )
            audit_orientation_hierarchy(table)
            if source.label.startswith("det-0"):
                audit_zero_parameter_hierarchy(table)
            print(
                f"{source.label}: c^9 normalization, first/second exceptional "
                "forms, final ridge and orientation hierarchy PASS",
                flush=True,
            )
    for sign in (1, -1):
        for source in final_chart_tables(core, sign)[2:]:
            table = low_nome_table(source)
            if not low_parameter_taylor_maps(table, degree=0)[0]:
                raise AssertionError(
                    f"the c^2 leading minor vanished in {source.label}"
                )
            audit_minor_corner_hierarchy(table)
            print(
                f"{source.label}: exact c^2 normalization and corner form PASS",
                flush=True,
            )


def audit_minor_corner_hierarchy(table: ChartTable) -> None:
    """Verify the exact quadratic at each final minor's sharp corner."""

    power_map = exact_centered_taylor_map(table, degree=2)
    for axis in (1, 3):
        power_map = affine_sparse_map_axis(power_map, axis, Fraction(1), Fraction(-1))

    if table.label.startswith("minor-secondary"):
        transverse_axes = (0, 1, 2, 3, 5)
        expected: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(256),
            (0, 2, 0, 0, 0, 0): fmpq(64),
            (0, 1, 1, 0, 0, 0): fmpq(64),
            (0, 1, 0, 1, 0, 0): fmpq(64),
            (0, 0, 2, 0, 0, 0): fmpq(48),
            (0, 0, 0, 0, 0, 2): fmpq(48),
        }
    else:
        transverse_axes = (0, 1, 2, 3, 4)
        expected = {
            (2, 0, 0, 0, 0, 0): fmpq(256),
            (0, 2, 0, 0, 0, 0): fmpq(64),
            (0, 1, 1, 0, 0, 0): fmpq(64),
            (0, 0, 2, 0, 0, 0): fmpq(48),
            (0, 0, 1, 1, 0, 0): fmpq(96),
            (0, 0, 0, 2, 0, 0): fmpq(96),
            (0, 0, 0, 1, 1, 0): fmpq(96),
            (0, 0, 0, 0, 2, 0): fmpq(48),
        }
        if table.label.startswith("minor-tertiary-1"):
            expected[(0, 0, 1, 0, 1, 0)] = fmpq(96)

    leading = {
        monomial: coefficient
        for monomial, coefficient in power_map.items()
        if sum(monomial[axis] for axis in transverse_axes) == 2
    }
    if leading != expected:
        raise AssertionError(f"unexpected minor corner form in {table.label}")
    if any(
        coefficient
        for monomial, coefficient in power_map.items()
        if all(monomial[axis] == 0 for axis in transverse_axes)
    ):
        raise AssertionError(f"minor corner did not vanish in {table.label}")

    if table.label.startswith("minor-secondary"):
        ratio_chart = projective_sparse_map(power_map, transverse_axes, 3, order=2)
        final_axes = (0, 1, 2, 5)
        final_weights = (1, 2, 1, 1)
        final_leading = {
            monomial: coefficient
            for monomial, coefficient in ratio_chart.items()
            if sum(
                weight * monomial[axis]
                for axis, weight in zip(final_axes, final_weights)
            )
            == 2
        }
        expected_final: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(256),
            (0, 1, 0, 0, 0, 0): fmpq(64),
            (0, 0, 2, 0, 0, 0): fmpq(48),
            (0, 0, 0, 0, 0, 2): fmpq(48),
        }
        if final_leading != expected_final:
            raise AssertionError(
                f"unexpected secondary-minor final form in {table.label}"
            )
        away = exact_centered_taylor_map(table, degree=2)
        away = affine_sparse_map_axis(away, 1, Fraction(1), Fraction(-1))
        away_axes = (0, 1, 2, 5)
        away_weights = (1, 2, 1, 1)
        away_leading = {
            monomial: coefficient
            for monomial, coefficient in away.items()
            if sum(
                weight * monomial[axis] for axis, weight in zip(away_axes, away_weights)
            )
            == 2
        }
        expected_away: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(256),
            (0, 1, 0, 0, 0, 0): fmpq(64),
            (0, 1, 0, 1, 0, 0): fmpq(-64),
            (0, 0, 2, 0, 0, 0): fmpq(48),
            (0, 0, 0, 0, 0, 2): fmpq(48),
        }
        if away_leading != expected_away:
            raise AssertionError(
                f"unexpected secondary-minor ratio-away form in {table.label}"
            )
        origin = exact_centered_taylor_map(table, degree=2)
        origin_axes = (0, 1, 2, 3, 5)
        origin_weights = (1, 2, 1, 2, 1)
        origin_leading = {
            monomial: coefficient
            for monomial, coefficient in origin.items()
            if sum(
                weight * monomial[axis]
                for axis, weight in zip(origin_axes, origin_weights)
            )
            == 2
        }
        expected_origin: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(256),
            (0, 1, 0, 0, 0, 0): fmpq(64),
            (0, 0, 2, 0, 0, 0): fmpq(48),
            (0, 0, 0, 1, 0, 0): fmpq(64),
            (0, 0, 0, 0, 0, 2): fmpq(48),
        }
        if origin_leading != expected_origin:
            raise AssertionError(
                f"unexpected secondary-minor origin form in {table.label}"
            )

    midpoint = exact_centered_taylor_map(table, degree=2)
    midpoint = affine_sparse_map_axis(midpoint, 1, Fraction(1, 2), Fraction(1))
    midpoint = affine_sparse_map_axis(midpoint, 2, Fraction(1), Fraction(-1))
    if not table.label.startswith("minor-secondary"):
        midpoint = affine_sparse_map_axis(midpoint, 3, Fraction(1), Fraction(-1))
    midpoint_axes, _, _, _, _ = minor_midpoint_configuration(table)
    midpoint_weights = tuple(2 if axis == 2 else 1 for axis in midpoint_axes)
    midpoint_leading = {
        monomial: coefficient
        for monomial, coefficient in midpoint.items()
        if sum(
            weight * monomial[axis]
            for axis, weight in zip(midpoint_axes, midpoint_weights)
        )
        == 2
    }
    sign = 1 if "sign-+1" in table.label else -1
    expected_midpoint: SparseMap = {
        (2, 0, 0, 0, 0, 0): fmpq(720),
        (0, 2, 0, 0, 0, 0): fmpq(384),
        (0, 0, 1, 0, 0, 0): fmpq(160),
    }
    if table.label.startswith("minor-secondary"):
        expected_midpoint = {
            (2, 0, 0, 0, 0, 0): fmpq(576),
            (2, 0, 0, 1, 0, 0): fmpq(144),
            (0, 2, 0, 0, 0, 0): fmpq(384),
            (0, 0, 1, 0, 0, 0): fmpq(64),
            (0, 0, 1, 1, 0, 0): fmpq(96),
            (0, 0, 0, 0, 0, 2): fmpq(144),
            (1, 0, 0, 0, 0, 1): fmpq(-576 * sign),
        }
    elif table.label.startswith("minor-tertiary-0"):
        expected_midpoint.update(
            {
                (0, 0, 0, 2, 0, 0): fmpq(144),
                (1, 0, 0, 1, 0, 0): fmpq(-576 * sign),
                (0, 0, 0, 0, 2, 0): fmpq(144),
                (0, 0, 0, 1, 1, 0): fmpq(288),
                (1, 0, 0, 0, 1, 0): fmpq(-576 * sign),
            }
        )
    else:
        expected_midpoint = {
            (0, 0, 0, 2, 0, 0): fmpq(144),
            (0, 0, 0, 2, 1, 0): fmpq(-144),
            (0, 0, 1, 0, 0, 0): fmpq(160),
            (0, 0, 1, 0, 1, 0): fmpq(-224),
            (0, 0, 1, 0, 2, 0): fmpq(64),
            (0, 2, 0, 0, 0, 0): fmpq(384),
            (0, 2, 0, 0, 1, 0): fmpq(-384),
            (1, 0, 0, 1, 0, 0): fmpq(-576 * sign),
            (1, 0, 0, 1, 1, 0): fmpq(576 * sign),
            (2, 0, 0, 0, 0, 0): fmpq(720),
            (2, 0, 0, 0, 1, 0): fmpq(-576),
        }
    if midpoint_leading != expected_midpoint:
        raise AssertionError(f"unexpected minor midpoint form in {table.label}")

    if table.label.startswith("minor-tertiary-1"):
        ratio_zero = exact_centered_taylor_map(table, degree=2)
        ratio_zero = affine_sparse_map_axis(ratio_zero, 4, Fraction(1), Fraction(-1))
        ratio_zero = projective_sparse_map(
            ratio_zero, (0, 4), 4, order=2, weights=(1, 2)
        )
        for axis in (2, 3):
            ratio_zero = affine_sparse_map_axis(
                ratio_zero, axis, Fraction(1), Fraction(-1)
            )
        ratio_zero = affine_sparse_map_axis(ratio_zero, 1, Fraction(1, 2), Fraction(1))
        ratio_axes = (0, 1, 2, 3)
        ratio_weights = (1, 1, 2, 1)
        ratio_leading = {
            monomial: coefficient
            for monomial, coefficient in ratio_zero.items()
            if sum(
                weight * monomial[axis]
                for axis, weight in zip(ratio_axes, ratio_weights)
            )
            == 2
        }
        expected_ratio: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(144),
            (2, 0, 0, 0, 2, 0): fmpq(576),
            (1, 0, 0, 1, 1, 0): fmpq(-576 * sign),
            (0, 2, 0, 0, 0, 0): fmpq(384),
            (0, 0, 1, 0, 0, 0): fmpq(96),
            (0, 0, 1, 0, 2, 0): fmpq(64),
            (0, 0, 0, 2, 0, 0): fmpq(144),
        }
        if ratio_leading != expected_ratio:
            raise AssertionError(
                f"unexpected tertiary ratio-zero form in {table.label}"
            )

    if table.label.startswith("minor-secondary") and sign == 1:
        ridge = exact_centered_taylor_map(table, degree=4, nome_scale=Fraction(1, 200))
        ridge = affine_sparse_map_axis(ridge, 2, Fraction(1), -(MINOR_CORNER_SCALE**2))
        ridge = affine_sparse_map_axis(ridge, 5, Fraction(0), MINOR_CORNER_SCALE)
        ridge = affine_sparse_map_axis(ridge, 1, Fraction(1, 2), -MINOR_CORNER_SCALE)
        ridge = projective_sparse_map(
            ridge,
            (0, 1, 2, 5),
            0,
            order=2,
            weights=(1, 1, 2, 1),
        )
        ridge = affine_sparse_map_axis(ridge, 5, Fraction(1, 25), Fraction(1))
        ridge_axes = (0, 1, 2, 3, 5)
        ridge_weights = (1, 1, 2, 2, 1)
        ridge_leading = {
            monomial: coefficient
            for monomial, coefficient in ridge.items()
            if sum(
                weight * monomial[axis]
                for axis, weight in zip(ridge_axes, ridge_weights)
            )
            == 2
        }
        expected_ridge: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(3, 12_500_000),
            (0, 2, 0, 0, 0, 0): fmpq(24),
            (0, 0, 1, 0, 0, 0): fmpq(4),
            (0, 0, 0, 1, 0, 0): fmpq(9, 2_500),
            (0, 0, 0, 0, 0, 2): fmpq(9),
        }
        if ridge_leading != expected_ridge:
            raise AssertionError(
                f"unexpected secondary-minor ridge form in {table.label}"
            )


def add_centered_parameter_term(
    output: SparseMap,
    nome_power: int,
    parameters: tuple[int, ...],
    value: SparseCoefficient,
) -> None:
    """Center ``v`` and ``a`` at one in a six-variable sparse power map."""

    u_power, v_power, y_power, a_power, b_power = parameters
    for centered_v in range(v_power + 1):
        for centered_a in range(a_power + 1):
            key = (
                nome_power,
                u_power,
                centered_v,
                y_power,
                centered_a,
                b_power,
            )
            contribution = (
                value
                * comb(v_power, centered_v)
                * (-1) ** centered_v
                * comb(a_power, centered_a)
                * (-1) ** centered_a
            )
            output[key] = output.get(key, fmpq(0)) + contribution


def exact_centered_taylor_map(
    table: ChartTable, degree: int, *, nome_scale: Fraction = Fraction(1)
) -> SparseMap:
    """Return the exact centered Taylor polynomial through ``degree``."""

    output: SparseMap = {}
    for nome_power, coefficient_map in enumerate(
        low_parameter_taylor_maps(table, degree=degree)
    ):
        scale = nome_scale**nome_power
        scale_exact = fmpq(scale.numerator, scale.denominator)
        for parameters, coefficient in coefficient_map.items():
            add_centered_parameter_term(
                output,
                nome_power,
                parameters,
                coefficient * scale_exact,
            )
    return {
        key: value for key, value in output.items() if not coefficient_is_zero(value)
    }


def asymptotic_power_map(
    table: ChartTable,
    c_upper: Fraction,
    *,
    remainder_degree: int = 16,
) -> SparseMap:
    """Build an exact Maclaurin polynomial plus parity-aware interval tails."""

    if remainder_degree != 16:
        raise ValueError(
            "the exact center series are currently audited through order fifteen"
        )
    ctx.dps = 80
    output = exact_centered_taylor_map(table, remainder_degree - 1, nome_scale=c_upper)

    keys = {monomial[5:] for monomial in table.records}
    odd_remainder_degree = remainder_degree + 1
    factors = low_factor_jets((Fraction(0), c_upper), keys, odd_remainder_degree)
    remainders: dict[int, defaultdict[tuple[int, ...], arb]] = {
        remainder_degree: defaultdict(arb),
        odd_remainder_degree: defaultdict(arb),
    }
    for monomial, coefficient in table.records.items():
        factor_key = monomial[5:]
        tail_degree = odd_remainder_degree if factor_key[0] % 2 else remainder_degree
        remainders[tail_degree][monomial[:5]] += (
            coefficient * factors[factor_key][1].coefficients[tail_degree]
        )
    for tail_degree, remainder in remainders.items():
        scale = c_upper**tail_degree
        scale_ball = arb(scale.numerator) / scale.denominator
        for parameters, coefficient in remainder.items():
            add_centered_parameter_term(
                output,
                tail_degree,
                parameters,
                coefficient * scale_ball,
            )
    return {
        key: value for key, value in output.items() if not coefficient_is_zero(value)
    }


def coefficient_is_zero(value: SparseCoefficient) -> bool:
    """Return whether a sparse coefficient is exactly zero."""

    return value.is_zero() if isinstance(value, arb) else value == 0


def coefficient_contains_zero(value: SparseCoefficient) -> bool:
    """Return whether a sparse coefficient may be zero."""

    return value.contains(0) if isinstance(value, arb) else value == 0


def fraction_as_fmpq(value: Fraction) -> fmpq:
    """Convert a standard-library rational without passing through binary64."""

    return fmpq(value.numerator, value.denominator)


def centered_sparse_axis_model(
    power_map: SparseMap,
    axis: int,
    *,
    upper: Fraction,
    width: Fraction,
    degree: int,
) -> SparseMap:
    """Taylor-enclose one power axis on ``[upper-width, upper]``.

    Terms below ``degree`` are shifted exactly.  The remaining binomial tail
    is stored as one Arb coefficient multiplying the final unit-coordinate
    power.  This preserves low-order cancellations without retaining a large
    projective degree.
    """

    if not 0 < width <= upper <= 1 or degree < 1:
        raise ValueError("require 0 < width <= upper <= 1 and positive degree")
    grouped: defaultdict[tuple[int, ...], list[tuple[int, SparseCoefficient]]] = (
        defaultdict(list)
    )
    for monomial, coefficient in power_map.items():
        target = list(monomial)
        target[axis] = 0
        grouped[tuple(target)].append((monomial[axis], coefficient))

    maximum_exponent = max(monomial[axis] for monomial in power_map)
    tail_weights: dict[int, fmpq] = {}
    for exponent in range(maximum_exponent + 1):
        weight = sum(
            (
                Fraction(comb(exponent, order))
                * upper ** (exponent - order)
                * width**order
                for order in range(degree, exponent + 1)
            ),
            Fraction(),
        )
        tail_weights[exponent] = fraction_as_fmpq(weight)

    center = fraction_as_fmpq(upper)
    negative_width = -fraction_as_fmpq(width)
    output: SparseMap = {}
    for monomial, terms in grouped.items():
        for order in range(degree):
            coefficient = (
                sum(
                    (
                        value * comb(exponent, order) * center ** (exponent - order)
                        for exponent, value in terms
                        if exponent >= order
                    ),
                    fmpq(0),
                )
                * negative_width**order
            )
            if coefficient_is_zero(coefficient):
                continue
            target = list(monomial)
            target[axis] = order
            output[tuple(target)] = coefficient

        radius = sum(
            (
                abs(value if isinstance(value, arb) else arb(value))
                * arb(tail_weights[exponent])
                for exponent, value in terms
                if exponent >= degree
            ),
            arb(0),
        )
        if radius.is_zero():
            continue
        target = list(monomial)
        target[axis] = degree
        output[tuple(target)] = (-radius).union(radius)
    return output


def iter_half_axis_models(
    power_map: SparseMap,
    axes: tuple[tuple[int, str], ...],
    prefix: str = "",
) -> Iterator[tuple[str, SparseMap]]:
    """Generate a memory-bounded Cartesian cover by half-interval models."""

    if not axes:
        yield prefix, power_map
        return
    (axis, label), *remaining = axes
    for upper in (Fraction(1, 2), Fraction(1)):
        modeled = centered_sparse_axis_model(
            power_map,
            axis,
            upper=upper,
            width=Fraction(1, 2),
            degree=8,
        )
        suffix = f"{prefix}-{label}-{upper - Fraction(1, 2)}-{upper}"
        yield from iter_half_axis_models(modeled, tuple(remaining), suffix)
        del modeled
        gc.collect()


def projective_sparse_map(
    power_map: SparseMap,
    projective_axes: tuple[int, ...],
    chart_axis: int,
    *,
    order: int,
    weights: tuple[int, ...] | None = None,
) -> SparseMap:
    """Apply one exact largest-coordinate substitution to a sparse power map."""

    if chart_axis not in projective_axes:
        raise ValueError(chart_axis)
    if weights is None:
        weights = (1,) * len(projective_axes)
    if len(weights) != len(projective_axes) or any(weight < 1 for weight in weights):
        raise ValueError("projective weights must be positive and match the axes")
    chart_map: SparseMap = {}
    for monomial, coefficient in power_map.items():
        total = sum(
            weight * monomial[axis] for axis, weight in zip(projective_axes, weights)
        )
        if total < order:
            if not coefficient_contains_zero(coefficient):
                raise AssertionError("claimed projective order did not vanish")
            continue
        target = list(monomial)
        target[chart_axis] = total - order
        target_tuple = tuple(target)
        if target_tuple in chart_map:
            raise AssertionError("projective exponent map unexpectedly collided")
        chart_map[target_tuple] = coefficient
    return chart_map


def multiply_sparse_maps(left: SparseMap, right: SparseMap) -> SparseMap:
    """Multiply two exact sparse maps used by structural audits."""

    output: SparseMap = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                left_power + right_power
                for left_power, right_power in zip(left_monomial, right_monomial)
            )
            contribution = left_coefficient * right_coefficient
            output[monomial] = output.get(monomial, fmpq(0)) + contribution
    return {
        key: value for key, value in output.items() if not coefficient_is_zero(value)
    }


def audit_orientation_hierarchy(table: ChartTable) -> None:
    """Verify the exact factorized orientation and endpoint leading forms."""

    power_map = exact_centered_taylor_map(table, degree=6)
    projective_axes = (0, 1, 2, 4, 5)
    main_chart = projective_sparse_map(power_map, projective_axes, 1, order=1)
    for axis in (1, 2, 4):
        main_chart = affine_sparse_map_axis(main_chart, axis, Fraction(1), Fraction(-1))

    orientation_axes = (0, 1, 4, 5)
    orientation_weights = (1, 2, 1, 1)
    square_form: SparseMap = {
        (2, 0, 0, 0, 0, 0): fmpq(16),
        (0, 0, 0, 0, 2, 0): fmpq(3),
        (0, 0, 0, 0, 0, 2): fmpq(3),
    }
    if table.label.startswith("det-0"):
        orientation_order = 4
        first_factor = dict(square_form)
        first_factor[(0, 1, 1, 0, 0, 0)] = fmpq(4)
        second_factor = dict(square_form)
        second_factor[(0, 1, 0, 0, 0, 0)] = fmpq(4)
        second_factor[(0, 1, 1, 0, 0, 0)] = fmpq(-4)
    else:
        orientation_order = 2
        first_factor = {
            (0, 0, 0, 0, 0, 0): fmpq(1),
            (0, 0, 1, 0, 0, 0): fmpq(-1),
        }
        second_factor = dict(square_form)
        second_factor[(0, 1, 0, 0, 0, 0)] = fmpq(4)
        second_factor[(0, 1, 1, 0, 0, 0)] = fmpq(-4)
    leading = {
        monomial: coefficient
        for monomial, coefficient in main_chart.items()
        if sum(
            weight * monomial[axis]
            for axis, weight in zip(orientation_axes, orientation_weights)
        )
        == orientation_order
    }
    expected_leading = {
        key: (16 if orientation_order == 4 else 64) * value
        for key, value in multiply_sparse_maps(first_factor, second_factor).items()
    }
    if leading != expected_leading:
        raise AssertionError(f"unexpected orientation factorization in {table.label}")

    weighted_main_chart = projective_sparse_map(
        main_chart,
        orientation_axes,
        1,
        order=orientation_order,
        weights=orientation_weights,
    )
    if orientation_order == 2:
        ratio_one = affine_sparse_map_axis(main_chart, 2, Fraction(1), Fraction(-1))
        endpoint_axes = (0, 1, 2, 4, 5)
        endpoint_weights = (1, 2, 2, 1, 1)
        endpoint_leading = {
            monomial: coefficient
            for monomial, coefficient in ratio_one.items()
            if sum(
                weight * monomial[axis]
                for axis, weight in zip(endpoint_axes, endpoint_weights)
            )
            == 4
        }
        endpoint_second_factor = dict(square_form)
        endpoint_second_factor[(0, 1, 0, 0, 0, 0)] = fmpq(4)
        endpoint_second_factor[(0, 0, 1, 0, 0, 0)] = fmpq(4)
        expected_endpoint_leading = {
            key: 16 * value
            for key, value in multiply_sparse_maps(
                square_form, endpoint_second_factor
            ).items()
        }
        if endpoint_leading != expected_endpoint_leading:
            raise AssertionError(
                f"unexpected second-chart endpoint form in {table.label}"
            )

        endpoint_face_leading = {
            monomial: coefficient
            for monomial, coefficient in ratio_one.items()
            if monomial[0] == monomial[4] == monomial[5] == 0
            and 2 * (monomial[1] + monomial[2]) == 6
        }
        first_face_factor: SparseMap = {(0, 1, 0, 0, 0, 0): fmpq(1)}
        second_face_factor: SparseMap = {
            (0, 0, 1, 0, 0, 0): fmpq(4),
            (0, 1, 0, 0, 0, 0): fmpq(7),
        }
        third_face_factor: SparseMap = {
            (0, 0, 1, 0, 0, 0): fmpq(4),
            (0, 1, 0, 0, 0, 0): fmpq(4),
        }
        expected_face_leading = {
            key: 16 * value
            for key, value in multiply_sparse_maps(
                multiply_sparse_maps(first_face_factor, second_face_factor),
                third_face_factor,
            ).items()
        }
        if endpoint_face_leading != expected_face_leading:
            raise AssertionError(
                f"unexpected second-chart endpoint-face form in {table.label}"
            )

        endpoint_main_chart = projective_sparse_map(
            ratio_one,
            endpoint_axes,
            2,
            order=4,
            weights=endpoint_weights,
        )
        transverse_axes = (0, 1, 4, 5)
        transverse_weights = (1, 2, 1, 1)
        transverse_leading = {
            monomial: coefficient
            for monomial, coefficient in endpoint_main_chart.items()
            if sum(
                weight * monomial[axis]
                for axis, weight in zip(transverse_axes, transverse_weights)
            )
            == 2
        }
        expected_transverse: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(1_024),
            (0, 0, 0, 0, 2, 0): fmpq(192),
            (0, 0, 0, 0, 0, 2): fmpq(192),
            (0, 1, 2, 0, 0, 0): fmpq(256),
        }
        if transverse_leading != expected_transverse:
            raise AssertionError(
                f"unexpected second-chart transverse form in {table.label}"
            )

        transverse_main_chart = projective_sparse_map(
            endpoint_main_chart,
            transverse_axes,
            1,
            order=2,
            weights=transverse_weights,
        )
        final_axes = (0, 2, 4, 5)
        final_leading = {
            monomial: coefficient
            for monomial, coefficient in transverse_main_chart.items()
            if sum(monomial[axis] for axis in final_axes) == 2
        }
        expected_final: SparseMap = {
            (2, 0, 0, 0, 0, 0): fmpq(1_024),
            (2, 2, 0, 0, 0, 0): fmpq(1_024),
            (0, 0, 2, 0, 0, 0): fmpq(256),
            (0, 1, 1, 0, 1, 0): fmpq(384),
            (0, 2, 0, 0, 0, 2): fmpq(192),
            (0, 2, 0, 0, 2, 0): fmpq(192),
            (0, 2, 2, 0, 0, 0): fmpq(704),
            (0, 3, 1, 0, 1, 0): fmpq(384),
            (0, 4, 2, 0, 0, 0): fmpq(448),
            (0, 0, 0, 0, 0, 2): fmpq(192),
            (0, 0, 0, 0, 2, 0): fmpq(192),
        }
        if final_leading != expected_final:
            raise AssertionError(
                f"unexpected second-chart final quadratic in {table.label}"
            )
        return
    expected_endpoint: SparseMap = {
        (2, 0, 0, 0, 0, 0): fmpq(1_024),
        (0, 2, 0, 0, 0, 0): fmpq(448),
        (0, 1, 0, 0, 1, 0): fmpq(384),
        (0, 0, 0, 0, 2, 0): fmpq(192),
        (0, 0, 0, 0, 0, 2): fmpq(192),
    }
    for endpoint_map in (
        weighted_main_chart,
        affine_sparse_map_axis(weighted_main_chart, 2, Fraction(1), Fraction(-1)),
    ):
        endpoint_face = {
            monomial: coefficient
            for monomial, coefficient in endpoint_map.items()
            if monomial[2] == 0
            and sum(monomial[axis] for axis in orientation_axes) == 2
        }
        if endpoint_face != expected_endpoint:
            raise AssertionError(
                f"unexpected orientation endpoint form in {table.label}"
            )


def audit_zero_parameter_hierarchy(table: ChartTable) -> None:
    """Verify the exact weighted forms at the remaining ``a=0`` corner."""

    power_map = exact_centered_taylor_map(table, degree=6)
    projective_axes = (0, 1, 2, 4, 5)
    ratio_chart = projective_sparse_map(power_map, projective_axes, 2, order=1)
    for axis in (1, 2, 4):
        ratio_chart = affine_sparse_map_axis(
            ratio_chart, axis, Fraction(1), Fraction(-1)
        )

    weights = (1, 2, 2, 2, 1)
    leading = {
        monomial: coefficient
        for monomial, coefficient in ratio_chart.items()
        if sum(
            weight * monomial[axis] for axis, weight in zip(projective_axes, weights)
        )
        == 4
    }
    square_form: SparseMap = {
        (2, 0, 0, 0, 0, 0): fmpq(16),
        (0, 0, 0, 0, 0, 2): fmpq(3),
    }
    second_factor = dict(square_form)
    second_factor[(0, 1, 0, 0, 0, 0)] = fmpq(4)
    second_factor[(0, 0, 1, 0, 0, 0)] = fmpq(4)
    expected_leading = {
        key: 16 * value
        for key, value in multiply_sparse_maps(square_form, second_factor).items()
    }
    if leading != expected_leading:
        raise AssertionError(f"unexpected zero-parameter leading form in {table.label}")

    main_chart = projective_sparse_map(
        ratio_chart,
        projective_axes,
        1,
        order=4,
        weights=weights,
    )
    main_chart = compress_even_sparse_map_axis(main_chart, 1)
    transverse_axes = (0, 2, 4, 5)
    transverse_weights = (1, 2, 1, 1)
    transverse_leading = {
        monomial: coefficient
        for monomial, coefficient in main_chart.items()
        if sum(
            weight * monomial[axis]
            for axis, weight in zip(transverse_axes, transverse_weights)
        )
        == 2
    }
    expected_transverse: SparseMap = {
        (2, 0, 0, 0, 0, 0): fmpq(1_024),
        (0, 0, 0, 0, 0, 2): fmpq(192),
        (0, 1, 0, 0, 2, 0): fmpq(192),
        (0, 1, 1, 0, 0, 0): fmpq(256),
    }
    if transverse_leading != expected_transverse:
        raise AssertionError(
            f"unexpected zero-parameter transverse form in {table.label}"
        )

    for selected_axis, expected_final in (
        (
            2,
            {
                (2, 0, 0, 0, 0, 0): fmpq(1_024),
                (2, 0, 2, 0, 0, 0): fmpq(1_024),
                (0, 0, 0, 0, 0, 2): fmpq(192),
                (0, 0, 2, 0, 0, 2): fmpq(192),
                (0, 1, 0, 0, 0, 0): fmpq(256),
                (0, 1, 0, 0, 2, 0): fmpq(192),
                (0, 1, 1, 0, 1, 0): fmpq(384),
                (0, 1, 2, 0, 0, 0): fmpq(704),
                (0, 1, 2, 0, 2, 0): fmpq(192),
                (0, 1, 3, 0, 1, 0): fmpq(384),
                (0, 1, 4, 0, 0, 0): fmpq(448),
            },
        ),
        (
            4,
            {
                (2, 0, 0, 0, 0, 0): fmpq(1_024),
                (2, 0, 1, 0, 2, 0): fmpq(1_024),
                (0, 0, 0, 0, 0, 2): fmpq(192),
                (0, 0, 1, 0, 2, 2): fmpq(192),
                (0, 1, 0, 0, 0, 0): fmpq(192),
                (0, 1, 1, 0, 0, 0): fmpq(256),
                (0, 1, 1, 0, 1, 0): fmpq(384),
                (0, 1, 1, 0, 2, 0): fmpq(192),
                (0, 1, 2, 0, 2, 0): fmpq(704),
                (0, 1, 2, 0, 3, 0): fmpq(384),
                (0, 1, 3, 0, 4, 0): fmpq(448),
            },
        ),
    ):
        selected = projective_sparse_map(
            main_chart,
            transverse_axes,
            selected_axis,
            order=2,
            weights=transverse_weights,
        )
        final_axes = (0, 1, 5)
        final_weights = (1, 2, 1)
        final_leading = {
            monomial: coefficient
            for monomial, coefficient in selected.items()
            if sum(
                weight * monomial[axis]
                for axis, weight in zip(final_axes, final_weights)
            )
            == 2
        }
        if final_leading != expected_final:
            raise AssertionError(
                f"unexpected zero-parameter final form in {table.label}"
            )
        if any(
            coefficient
            for monomial, coefficient in selected.items()
            if monomial[0] == monomial[1] == monomial[5] == 0
        ):
            raise AssertionError(
                f"zero-parameter final equality line did not vanish in {table.label}"
            )


def affine_sparse_map_axis(
    power_map: SparseMap,
    axis: int,
    center: Fraction,
    slope: Fraction,
) -> SparseMap:
    """Replace one power variable by ``center + slope*unit`` exactly."""

    center_exact = fmpq(center.numerator, center.denominator)
    slope_exact = fmpq(slope.numerator, slope.denominator)
    if center == 0:
        return {
            monomial: coefficient * slope_exact ** monomial[axis]
            for monomial, coefficient in power_map.items()
        }
    output: SparseMap = {}
    for monomial, coefficient in power_map.items():
        exponent = monomial[axis]
        for centered_power in range(exponent + 1):
            target = list(monomial)
            target[axis] = centered_power
            target_tuple = tuple(target)
            contribution = (
                coefficient
                * comb(exponent, centered_power)
                * center_exact ** (exponent - centered_power)
                * slope_exact**centered_power
            )
            output[target_tuple] = output.get(target_tuple, fmpq(0)) + contribution
    return {
        key: value for key, value in output.items() if not coefficient_is_zero(value)
    }


def compress_even_sparse_map_axis(power_map: SparseMap, axis: int) -> SparseMap:
    """Replace an even power variable by its square on the unit interval."""

    if any(monomial[axis] % 2 for monomial in power_map):
        raise ValueError(f"sparse map is not even on axis {axis}")
    output: SparseMap = {}
    for monomial, coefficient in power_map.items():
        target = list(monomial)
        target[axis] //= 2
        target_tuple = tuple(target)
        if target_tuple in output:
            raise AssertionError("even-axis compression unexpectedly collided")
        output[target_tuple] = coefficient
    return output


def sparse_map_bernstein_tensor(
    power_map: SparseMap, *, axis_order: tuple[int, ...] | None = None
) -> IntervalTensor:
    """Convert one exact/Arb sparse map to a directed Bernstein tensor."""

    dimensions = tuple(max(key[axis] for key in power_map) + 1 for axis in range(6))
    interval_map = {
        key: value if isinstance(value, arb) else arb(value)
        for key, value in power_map.items()
    }
    return power_to_bernstein_tensor(
        arb_map_to_tensor(interval_map, dimensions), axis_order=axis_order
    )


def collapsed_sparse_map_bernstein_tensor(
    power_map: SparseMap,
    primary_axis: int,
    *,
    envelope_axis: int = 3,
    collapse_envelope_early: bool = False,
) -> IntervalTensor:
    """Convert a chart efficiently and enclose its irrelevant envelope axis."""

    axis_order = (primary_axis,) + tuple(
        axis for axis in range(6) if axis != primary_axis
    )
    non_envelope_order = (primary_axis,) + tuple(
        axis for axis in range(6) if axis not in (primary_axis, envelope_axis)
    )
    if collapse_envelope_early:
        dimensions = tuple(max(key[axis] for key in power_map) + 1 for axis in range(6))
        interval_map = {
            key: value if isinstance(value, arb) else arb(value)
            for key, value in power_map.items()
        }
        tensor = arb_map_to_tensor(interval_map, dimensions)
        tensor = power_to_bernstein_axis(tensor, envelope_axis)
        tensor = collapse_bernstein_axis(tensor, envelope_axis)
        for axis in non_envelope_order:
            tensor = power_to_bernstein_axis(tensor, axis)
        return tensor
    return collapse_bernstein_axis(
        sparse_map_bernstein_tensor(power_map, axis_order=axis_order), envelope_axis
    )


def orientation_endpoint_maps(
    power_map: SparseMap,
) -> Iterator[tuple[str, SparseMap, tuple[tuple[float, float], ...]]]:
    """Restrict the two sharp endpoints of the weighted orientation chart."""

    scale = ORIENTATION_CORNER_SCALE
    endpoint_width = scale**2
    parameterizations = (
        ("ratio-zero", Fraction(0), endpoint_width),
        ("ratio-one", Fraction(1), -endpoint_width),
    )
    for label, center, slope in parameterizations:
        restricted = affine_sparse_map_axis(power_map, 2, center, slope)
        for axis in (0, 1, 4, 5):
            restricted = affine_sparse_map_axis(restricted, axis, Fraction(0), scale)
        endpoint_interval = (
            (0.0, float(endpoint_width))
            if center == 0
            else (float(1 - endpoint_width), 1.0)
        )
        box = (
            (0.0, float(scale)),
            (0.0, float(scale)),
            endpoint_interval,
            (0.0, 1.0),
            (0.0, float(scale)),
            (0.0, float(scale)),
        )
        yield label, restricted, box


def restricted_zero_parameter_map(power_map: SparseMap, scale: Fraction) -> SparseMap:
    """Center and scale the five active det0 zero-parameter coordinates."""

    restricted = power_map
    for axis in (0, 5):
        restricted = affine_sparse_map_axis(restricted, axis, Fraction(0), scale)
    for axis in (1, 2, 4):
        restricted = affine_sparse_map_axis(restricted, axis, Fraction(1), -(scale**2))
    return restricted


def widened_zero_parameter_corner_tensors(
    power_map: SparseMap,
) -> Iterator[tuple[str, IntervalTensor]]:
    """Build four overlapping weighted charts around the ``a=0`` corner."""

    scale = ZERO_PARAMETER_CORNER_SCALE
    restricted = restricted_zero_parameter_map(power_map, scale)

    axes = (0, 1, 2, 4, 5)
    weights = (1, 2, 2, 2, 1)
    labels = {0: "nome", 2: "ratio-radius", 4: "a", 5: "abs-b"}
    for selected_axis in labels:
        chart = projective_sparse_map(
            restricted, axes, selected_axis, order=4, weights=weights
        )
        for ratio_axis in axes:
            if ratio_axis != selected_axis:
                chart = affine_sparse_map_axis(
                    chart, ratio_axis, Fraction(0), Fraction(2)
                )
        chart = compress_even_sparse_map_axis(chart, selected_axis)
        yield labels[selected_axis], sparse_map_bernstein_tensor(chart)


def zero_parameter_corner_box() -> Box:
    """Return the top-chart box covered by all five weighted corner charts."""

    scale = float(ZERO_PARAMETER_CORNER_SCALE)
    squared_scale = float(ZERO_PARAMETER_CORNER_SCALE**2)
    upper_defect = (1.0 - squared_scale, 1.0)
    return (
        (0.0, scale),
        upper_defect,
        upper_defect,
        (0.0, 1.0),
        upper_defect,
        (0.0, scale),
    )


def iter_asymptotic_chart_maps(
    power_map: SparseMap,
    chart_axis: int,
    label: str,
    c_upper: Fraction,
) -> Iterator[tuple[str, SparseMap]]:
    """Yield one memory-bounded final chart map at a time."""

    projective_axes = (0, 1, 2, 4, 5)
    if chart_axis == 0:
        first = projective_sparse_map(power_map, projective_axes, chart_axis, order=1)
        ridge_center = 2 * c_upper
        for side_label, slope in (
            ("left", -ridge_center),
            ("right", 1 - ridge_center),
        ):
            centered = affine_sparse_map_axis(first, 5, ridge_center, slope)
            for secondary_axis, secondary_label in zip(
                (0, 1, 4), ("c", "main", "one-minus-a")
            ):
                secondary = projective_sparse_map(
                    centered, (0, 1, 4), secondary_axis, order=1
                )
                prefix = f"{label}-{side_label}-secondary-{secondary_label}"
                if secondary_axis != 0:
                    yield prefix, secondary
                    continue
                for tertiary_axis, tertiary_label in zip(
                    projective_axes,
                    (
                        "c",
                        "main",
                        "one-minus-ratio",
                        "one-minus-a",
                        "ridge-offset",
                    ),
                ):
                    yield (
                        f"{prefix}-weighted-{tertiary_label}",
                        projective_sparse_map(
                            secondary,
                            projective_axes,
                            tertiary_axis,
                            order=2,
                            weights=(1, 2, 1, 2, 1),
                        ),
                    )
        return

    if chart_axis == 1:
        first = projective_sparse_map(power_map, projective_axes, chart_axis, order=1)
        centered = first
        for upper_axis in (1, 2, 4):
            centered = affine_sparse_map_axis(
                centered, upper_axis, Fraction(1), Fraction(-1)
            )
        orientation_axes = (0, 1, 4, 5)
        orientation_weights = (1, 2, 1, 1)
        orientation_order = min(
            sum(
                weight * monomial[axis]
                for axis, weight in zip(orientation_axes, orientation_weights)
            )
            for monomial, coefficient in centered.items()
            if not coefficient_is_zero(coefficient)
        )
        for nested_axis, nested_label in zip(
            orientation_axes, ("c", "one-minus-main", "a", "abs-b")
        ):
            yield (
                f"{label}-orientation-weighted{orientation_order}-{nested_label}",
                projective_sparse_map(
                    centered,
                    orientation_axes,
                    nested_axis,
                    order=orientation_order,
                    weights=orientation_weights,
                ),
            )
        return

    yield (
        label,
        projective_sparse_map(power_map, projective_axes, chart_axis, order=1),
    )


def certify_sparse_chart(
    label: str,
    power_map: SparseMap,
    primary_axis: int,
    *,
    max_depth: int,
    max_leaves: int,
    envelope_axis: int = 3,
    collapse_envelope_early: bool = False,
) -> None:
    """Build, print, and require one envelope-collapsed sparse certificate."""

    tensor = collapsed_sparse_map_bernstein_tensor(
        power_map,
        primary_axis,
        envelope_axis=envelope_axis,
        collapse_envelope_early=collapse_envelope_early,
    )
    result = certify(tensor, max_depth=max_depth, max_leaves=max_leaves)
    print(
        f"{label}: shape={tensor.lower.shape}, pass={result.passed}, "
        f"leaves={result.leaves}, depth={result.depth}, "
        f"lower={result.minimum_lower:.3e}",
        flush=True,
    )
    if not result.passed:
        print(
            f"  failure_box={result.failure_box}, index={result.failure_index}",
            flush=True,
        )
        raise SystemExit(1)


def second_determinant_endpoint_map(power_map: SparseMap) -> SparseMap:
    """Center the det1 orientation endpoint before its order-four blow-up."""

    axes = (0, 1, 2, 4, 5)
    centered = projective_sparse_map(power_map, axes, 1, order=1)
    for axis in (1, 2, 4):
        centered = affine_sparse_map_axis(centered, axis, Fraction(1), Fraction(-1))
    return affine_sparse_map_axis(centered, 2, Fraction(1), Fraction(-1))


def first_determinant_zero_parameter_main_map(
    ratio_chart: SparseMap, *, scale: Fraction = ZERO_PARAMETER_CORNER_SCALE
) -> SparseMap:
    """Build the hard U-dominant chart from the det0 ratio chart."""

    axes = (0, 1, 2, 4, 5)
    ratio_chart = restricted_zero_parameter_map(ratio_chart, scale)
    main_chart = projective_sparse_map(
        ratio_chart,
        axes,
        1,
        order=4,
        weights=(1, 2, 2, 2, 1),
    )
    return compress_even_sparse_map_axis(main_chart, 1)


def certify_first_determinant_zero_parameter_main(
    source: ChartTable,
    c_upper: Fraction,
    *,
    max_depth: int,
    max_leaves: int,
    ratio_chart: SparseMap | None = None,
) -> None:
    """Certify the final det0 U-dominant corner hierarchy."""

    if not source.label.startswith("det-0"):
        raise ValueError("the zero-parameter main chart is specific to det0")
    if c_upper > Fraction(1, 200):
        raise ValueError("the audited det0 main tail requires c_upper <= 1/200")

    started = time.monotonic()
    if ratio_chart is None:
        power_map = asymptotic_power_map(low_nome_table(source), c_upper)
        ratio_chart = projective_sparse_map(power_map, (0, 1, 2, 4, 5), 2, order=1)
        del power_map
        gc.collect()
    main_chart = first_determinant_zero_parameter_main_map(ratio_chart)

    transverse_axes = (0, 2, 4, 5)
    transverse_weights = (1, 2, 1, 1)
    for axis, label in ((0, "nome"), (5, "B")):
        chart = projective_sparse_map(
            main_chart,
            transverse_axes,
            axis,
            order=2,
            weights=transverse_weights,
        )
        for s_upper in (Fraction(1, 2), Fraction(1)):
            s_modeled = centered_sparse_axis_model(
                chart,
                1,
                upper=s_upper,
                width=Fraction(1, 2),
                degree=8,
            )
            for r_upper in (Fraction(1, 2), Fraction(1)):
                modeled = centered_sparse_axis_model(
                    s_modeled,
                    2,
                    upper=r_upper,
                    width=Fraction(1, 2),
                    degree=8,
                )
                certify_sparse_chart(
                    f"{source.label}-det0-main-{label}-"
                    f"S-{s_upper - Fraction(1, 2)}-{s_upper}-"
                    f"R-{r_upper - Fraction(1, 2)}-{r_upper}",
                    modeled,
                    axis,
                    max_depth=max_depth,
                    max_leaves=max_leaves,
                )
                del modeled
                gc.collect()
            del s_modeled
            gc.collect()
        del chart
        gc.collect()

    final_axes = (0, 1, 5)
    final_weights = (1, 2, 1)
    radial_configs = (
        (2, "R", ((2, "R"),)),
        (4, "A", ((4, "A"), (2, "R"))),
    )
    for transverse_axis, transverse_label, model_axes in radial_configs:
        transverse = projective_sparse_map(
            main_chart,
            transverse_axes,
            transverse_axis,
            order=2,
            weights=transverse_weights,
        )
        for final_axis, final_label in zip(final_axes, ("nome", "S", "B")):
            final_chart = projective_sparse_map(
                transverse,
                final_axes,
                final_axis,
                order=2,
                weights=final_weights,
            )
            for suffix, modeled in iter_half_axis_models(final_chart, model_axes):
                certify_sparse_chart(
                    f"{source.label}-det0-{transverse_label}-final-"
                    f"{final_label}{suffix}",
                    modeled,
                    final_axis,
                    max_depth=max_depth,
                    max_leaves=max_leaves,
                )
                del modeled
                gc.collect()
            del final_chart
            gc.collect()
        del transverse
        gc.collect()

    print(
        f"{source.label}: det0 zero-parameter main hierarchy PASS in "
        f"{time.monotonic() - started:.2f}s",
        flush=True,
    )


def first_determinant_zero_parameter_nome_arm_map(
    ratio_chart: SparseMap,
) -> SparseMap:
    """Factor the complete nome-dominant equality arm in the det0 top chart."""

    axes = (0, 1, 2, 4, 5)
    chart = projective_sparse_map(
        restricted_zero_parameter_map(ratio_chart, Fraction(1)),
        axes,
        0,
        order=4,
        weights=(1, 2, 2, 2, 1),
    )
    for axis in (1, 2, 4, 5):
        chart = affine_sparse_map_axis(
            chart, axis, Fraction(0), ZERO_PARAMETER_NOME_ARM_CAP
        )
    return compress_even_sparse_map_axis(chart, 0)


def certify_first_determinant_zero_parameter_nome_arm(
    source: ChartTable,
    ratio_chart: SparseMap,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify a uniform narrow tube around the full nome equality arm."""

    chart = first_determinant_zero_parameter_nome_arm_map(ratio_chart)
    certify_sparse_chart(
        f"{source.label}-det0-nome-arm",
        chart,
        0,
        max_depth=max_depth,
        max_leaves=max_leaves,
        collapse_envelope_early=True,
    )


def iter_first_determinant_zero_parameter_main_arm_maps(
    ratio_chart: SparseMap, cap: Fraction
) -> Iterator[tuple[str, SparseMap, int]]:
    """Yield a finite cover of the full main-dominant equality arm."""

    transverse_axes = (0, 2, 4, 5)
    weights = (1, 2, 1, 1)
    main_chart = first_determinant_zero_parameter_main_map(
        ratio_chart, scale=Fraction(1)
    )
    for axis in transverse_axes:
        main_chart = affine_sparse_map_axis(main_chart, axis, Fraction(0), cap)

    for selected_axis, label in zip(
        transverse_axes, ("nome", "ratio-radius", "a", "abs-b")
    ):
        chart = projective_sparse_map(
            main_chart,
            transverse_axes,
            selected_axis,
            order=2,
            weights=weights,
        )
        model_axes: tuple[tuple[int, str], ...] = (
            (selected_axis, "radial"),
            (2, "R"),
        )
        if selected_axis in (2, 4):
            chart = affine_sparse_map_axis(
                chart,
                1,
                ZERO_PARAMETER_MAIN_ARM_LOWER,
                1 - ZERO_PARAMETER_MAIN_ARM_LOWER,
            )
        if selected_axis == 4:
            model_axes += ((4, "A"),)
        for suffix, modeled in iter_half_axis_models(chart, model_axes):
            yield f"{label}{suffix}", modeled, selected_axis


def certify_first_determinant_zero_parameter_main_arm(
    source: ChartTable,
    ratio_chart: SparseMap,
    cap: Fraction,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify a largest-coordinate region around the full main equality arm."""

    started = time.monotonic()
    for (
        label,
        chart,
        selected_axis,
    ) in iter_first_determinant_zero_parameter_main_arm_maps(ratio_chart, cap):
        certify_sparse_chart(
            f"{source.label}-det0-main-arm-{label}",
            chart,
            selected_axis,
            max_depth=max_depth,
            max_leaves=max_leaves,
            collapse_envelope_early=True,
        )
        del chart
        gc.collect()
    print(
        f"{source.label}: det0 main-arm cover PASS in "
        f"{time.monotonic() - started:.2f}s",
        flush=True,
    )


def rounded_up(value: float) -> float:
    """Round one binary64 diagnostic bound conservatively upward."""

    return nextafter(value, inf)


def rounded_down(value: float) -> float:
    """Round one binary64 diagnostic bound conservatively downward."""

    return nextafter(value, -inf)


def inside_zero_parameter_nome_arm(box: Box) -> bool:
    """Recognize boxes contained in the certified nome-arm tube."""

    c_lower = box[0][0]
    if c_lower <= 0:
        return False
    cap = float(ZERO_PARAMETER_NOME_ARM_CAP)
    linear_limit = rounded_down(cap * c_lower)
    square_limit = rounded_down(cap * c_lower * c_lower)
    if box[5][1] > linear_limit:
        return False
    return all(rounded_up(1 - box[axis][0]) <= square_limit for axis in (1, 2, 4))


def inside_zero_parameter_main_arm(box: Box, cap: Fraction) -> bool:
    """Recognize boxes contained in the certified main-arm region."""

    main_lower = rounded_down(1 - box[1][1])
    if main_lower < float(ZERO_PARAMETER_MAIN_ARM_LOWER):
        return False
    cap_float = float(cap)
    linear_limit = rounded_down(cap_float * main_lower)
    square_limit = rounded_down(cap_float * cap_float * main_lower)
    if rounded_up(box[0][1] * box[0][1]) > square_limit:
        return False
    if rounded_up(box[5][1] * box[5][1]) > square_limit:
        return False
    return all(rounded_up(1 - box[axis][0]) <= linear_limit for axis in (2, 4))


def zero_parameter_main_arm_region(cap: Fraction) -> Callable[[Box], bool]:
    """Return the box predicate associated with one proved arm cap."""

    def contains(box: Box) -> bool:
        return inside_zero_parameter_main_arm(box, cap)

    return contains


def minor_corner_configuration(
    source: ChartTable,
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[int, ...], tuple[int, ...], int]:
    """Return projective axes, labels, upper/lower defects, and envelope axis."""

    if source.label.startswith("minor-secondary"):
        return (
            (0, 1, 2, 5, 3),
            ("nome", "main", "secondary", "a", "ratio"),
            (1, 3),
            (2, 5),
            4,
        )
    if source.label.startswith("minor-tertiary"):
        return (
            (0, 1, 2, 3, 4),
            ("nome", "main", "secondary", "tertiary", "ratio"),
            (1, 3),
            (2, 4),
            5,
        )
    raise ValueError(f"unsupported minor chart {source.label}")


def restricted_minor_corner_map(
    power_map: SparseMap, upper_axes: tuple[int, ...], lower_axes: tuple[int, ...]
) -> SparseMap:
    """Restrict one normalized minor map to its sharp parameter corner."""

    restricted = power_map
    for axis in upper_axes:
        restricted = affine_sparse_map_axis(
            restricted, axis, Fraction(1), -MINOR_CORNER_SCALE
        )
    for axis in lower_axes:
        restricted = affine_sparse_map_axis(
            restricted, axis, Fraction(0), MINOR_CORNER_SCALE
        )
    return restricted


def minor_corner_box(source: ChartTable) -> Box:
    """Return the original minor-coordinate box covered by its local charts."""

    _, _, upper_axes, lower_axes, _ = minor_corner_configuration(source)
    box: list[tuple[float, float]] = [(0.0, 1.0) for _ in range(6)]
    scale = float(MINOR_CORNER_SCALE)
    for axis in upper_axes:
        box[axis] = (1.0 - scale, 1.0)
    for axis in lower_axes:
        box[axis] = (0.0, scale)
    return tuple(box)


def certify_secondary_minor_ratio_chart(
    source: ChartTable,
    ratio_chart: SparseMap,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Resolve the secondary minor's sole equality line in its ratio chart."""

    axes = (0, 1, 2, 5)
    weights = (1, 2, 1, 1)
    for selected_axis, label in zip(axes, ("nome", "main", "secondary", "a")):
        chart = projective_sparse_map(
            ratio_chart, axes, selected_axis, order=2, weights=weights
        )
        certify_sparse_chart(
            f"{source.label}-minor-ratio-final-{label}",
            chart,
            selected_axis,
            max_depth=max_depth,
            max_leaves=max_leaves,
            envelope_axis=4,
            collapse_envelope_early=True,
        )
        del chart
        gc.collect()


def certify_secondary_corner_ratio_away(
    source: ChartTable,
    power_map: SparseMap,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify the secondary corner uniformly away from ratio one."""

    restricted = affine_sparse_map_axis(
        power_map, 1, Fraction(1), -(MINOR_CORNER_SCALE**2)
    )
    for axis in (2, 5):
        restricted = affine_sparse_map_axis(
            restricted, axis, Fraction(0), MINOR_CORNER_SCALE
        )
    restricted = affine_sparse_map_axis(
        restricted, 3, Fraction(0), 1 - MINOR_CORNER_SCALE
    )
    axes = (0, 1, 2, 5)
    weights = (1, 2, 1, 1)
    for selected_axis, label in zip(axes, ("nome", "main", "secondary", "a")):
        chart = projective_sparse_map(
            restricted, axes, selected_axis, order=2, weights=weights
        )
        certify_sparse_chart(
            f"{source.label}-minor-corner-ratio-away-{label}",
            chart,
            selected_axis,
            max_depth=max_depth,
            max_leaves=max_leaves,
            envelope_axis=4,
            collapse_envelope_early=True,
        )
        del chart
        gc.collect()


def secondary_corner_ratio_away_box() -> Box:
    """Return the secondary-corner box with ratio bounded away from one."""

    scale = float(MINOR_CORNER_SCALE)
    squared_scale = float(MINOR_CORNER_SCALE**2)
    return (
        (0.0, 1.0),
        (1.0 - squared_scale, 1.0),
        (0.0, scale),
        (0.0, 1.0 - scale),
        (0.0, 1.0),
        (0.0, scale),
    )


def certify_secondary_origin_corner(
    source: ChartTable,
    power_map: SparseMap,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify the secondary chart's all-lower weighted corner."""

    restricted = power_map
    for axis in (1, 3):
        restricted = affine_sparse_map_axis(
            restricted, axis, Fraction(0), MINOR_CORNER_SCALE**2
        )
    for axis in (2, 5):
        restricted = affine_sparse_map_axis(
            restricted, axis, Fraction(0), MINOR_CORNER_SCALE
        )
    axes = (0, 1, 2, 3, 5)
    weights = (1, 2, 1, 2, 1)
    labels = ("nome", "main", "secondary", "ratio", "a")
    for selected_axis, label in zip(axes, labels):
        chart = projective_sparse_map(
            restricted, axes, selected_axis, order=2, weights=weights
        )
        certify_sparse_chart(
            f"{source.label}-minor-origin-{label}",
            chart,
            selected_axis,
            max_depth=max_depth,
            max_leaves=max_leaves,
            envelope_axis=4,
            collapse_envelope_early=True,
        )
        del chart
        gc.collect()


def secondary_origin_corner_box() -> Box:
    """Return the all-lower secondary-minor corner box."""

    scale = float(MINOR_CORNER_SCALE)
    squared_scale = float(MINOR_CORNER_SCALE**2)
    return (
        (0.0, 1.0),
        (0.0, squared_scale),
        (0.0, scale),
        (0.0, squared_scale),
        (0.0, 1.0),
        (0.0, scale),
    )


def certify_tertiary_one_ratio_zero(
    source: ChartTable,
    power_map: SparseMap,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify tertiary chart one's complete ratio-zero boundary box."""

    restricted = affine_sparse_map_axis(
        power_map, 4, Fraction(1), -(MINOR_CORNER_SCALE**2)
    )
    nome_chart = projective_sparse_map(restricted, (0, 4), 0, order=2, weights=(1, 2))
    certify_sparse_chart(
        f"{source.label}-minor-ratio-zero-nome",
        nome_chart,
        0,
        max_depth=max_depth,
        max_leaves=max_leaves,
        envelope_axis=5,
        collapse_envelope_early=True,
    )
    del nome_chart
    gc.collect()

    ratio_chart = projective_sparse_map(restricted, (0, 4), 4, order=2, weights=(1, 2))
    for axis in (2, 3):
        ratio_chart = affine_sparse_map_axis(
            ratio_chart, axis, Fraction(1), Fraction(-1)
        )
    axes = (0, 1, 2, 3)
    weights = (1, 1, 2, 1)
    labels = ("nome", "main", "secondary", "tertiary")
    for side, slope in (("left", Fraction(-1, 2)), ("right", Fraction(1, 2))):
        centered = affine_sparse_map_axis(ratio_chart, 1, Fraction(1, 2), slope)
        for selected_axis, label in zip(axes, labels):
            chart = projective_sparse_map(
                centered, axes, selected_axis, order=2, weights=weights
            )
            certify_sparse_chart(
                f"{source.label}-minor-ratio-zero-{side}-{label}",
                chart,
                selected_axis,
                max_depth=max_depth,
                max_leaves=max_leaves,
                envelope_axis=5,
                collapse_envelope_early=True,
            )
            del chart
            gc.collect()
        del centered
        gc.collect()


def tertiary_one_ratio_zero_box() -> Box:
    """Return tertiary chart one's complete ratio-zero boundary box."""

    squared_scale = float(MINOR_CORNER_SCALE**2)
    box: list[tuple[float, float]] = [(0.0, 1.0) for _ in range(6)]
    box[4] = (1.0 - squared_scale, 1.0)
    return tuple(box)


def minor_midpoint_configuration(
    source: ChartTable,
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[int, ...], tuple[int, ...], int]:
    """Return the transverse coordinates for a minor's midpoint ridge."""

    if source.label.startswith("minor-secondary"):
        return (0, 1, 2, 5), ("nome", "main", "secondary", "a"), (), (5,), 4
    if source.label.startswith("minor-tertiary-0"):
        return (
            (0, 1, 2, 3, 4),
            ("nome", "main", "secondary", "tertiary", "ratio"),
            (3,),
            (4,),
            5,
        )
    if source.label.startswith("minor-tertiary-1"):
        return (
            (0, 1, 2, 3),
            ("nome", "main", "secondary", "tertiary"),
            (3,),
            (),
            5,
        )
    raise ValueError(f"unsupported minor chart {source.label}")


def iter_minor_midpoint_corner_maps(
    source: ChartTable, power_map: SparseMap
) -> Iterator[tuple[str, SparseMap]]:
    """Build the two signed-main charts around a minor midpoint ridge."""

    _, _, upper_axes, lower_axes, _ = minor_midpoint_configuration(source)
    base = affine_sparse_map_axis(power_map, 2, Fraction(1), -(MINOR_CORNER_SCALE**2))
    for axis in upper_axes:
        base = affine_sparse_map_axis(base, axis, Fraction(1), -MINOR_CORNER_SCALE)
    for axis in lower_axes:
        base = affine_sparse_map_axis(base, axis, Fraction(0), MINOR_CORNER_SCALE)
    if source.label.startswith("minor-tertiary-1"):
        base = affine_sparse_map_axis(base, 4, Fraction(0), 1 - MINOR_CORNER_SCALE**2)
    for side, slope in (
        ("left", -MINOR_CORNER_SCALE),
        ("right", MINOR_CORNER_SCALE),
    ):
        yield side, affine_sparse_map_axis(base, 1, Fraction(1, 2), slope)


def certify_secondary_midpoint_nome_chart(
    source: ChartTable,
    nome_chart: SparseMap,
    *,
    side: str,
    c_upper: Fraction,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Resolve the positive-sign secondary midpoint's ``a=2c`` ridge."""

    ridge_center = 2 * c_upper / MINOR_CORNER_SCALE
    axes = (0, 1, 2, 3, 5)
    weights = (1, 1, 2, 2, 1)
    labels = ("radial", "main", "secondary", "ratio", "a-offset")
    for offset_side, slope in (
        ("left", -ridge_center),
        ("right", 1 - ridge_center),
    ):
        centered = affine_sparse_map_axis(nome_chart, 5, ridge_center, slope)
        for selected_axis, label in zip(axes, labels):
            chart = projective_sparse_map(
                centered, axes, selected_axis, order=2, weights=weights
            )
            certify_sparse_chart(
                f"{source.label}-minor-midpoint-{side}-nome-ridge-"
                f"{offset_side}-{label}",
                chart,
                selected_axis,
                max_depth=max_depth,
                max_leaves=max_leaves,
                envelope_axis=4,
                collapse_envelope_early=True,
            )
            del chart
            gc.collect()
        del centered
        gc.collect()


def certify_minor_midpoint_corner(
    source: ChartTable,
    power_map: SparseMap,
    *,
    c_upper: Fraction,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify the weighted corner at secondary=0 and main=1/2."""

    axes, labels, _, _, envelope_axis = minor_midpoint_configuration(source)
    weights = tuple(2 if axis == 2 else 1 for axis in axes)
    for side, corner_map in iter_minor_midpoint_corner_maps(source, power_map):
        for selected_axis, label in zip(axes, labels):
            chart = projective_sparse_map(
                corner_map, axes, selected_axis, order=2, weights=weights
            )
            if (
                source.label.startswith("minor-secondary")
                and source.label.endswith("sign-+1")
                and label == "nome"
            ):
                certify_secondary_midpoint_nome_chart(
                    source,
                    chart,
                    side=side,
                    c_upper=c_upper,
                    max_depth=max_depth,
                    max_leaves=max_leaves,
                )
                del chart
                gc.collect()
                continue
            certify_sparse_chart(
                f"{source.label}-minor-midpoint-{side}-{label}",
                chart,
                selected_axis,
                max_depth=max_depth,
                max_leaves=max_leaves,
                envelope_axis=envelope_axis,
                collapse_envelope_early=True,
            )
            del chart
            gc.collect()
        del corner_map
        gc.collect()


def minor_midpoint_corner_boxes(source: ChartTable) -> tuple[Box, Box]:
    """Return the two original-coordinate boxes around the midpoint ridge."""

    _, _, upper_axes, lower_axes, _ = minor_midpoint_configuration(source)
    scale = float(MINOR_CORNER_SCALE)
    squared_scale = float(MINOR_CORNER_SCALE**2)
    shared: list[tuple[float, float]] = [(0.0, 1.0) for _ in range(6)]
    shared[2] = (1.0 - squared_scale, 1.0)
    for axis in upper_axes:
        shared[axis] = (1.0 - scale, 1.0)
    for axis in lower_axes:
        shared[axis] = (0.0, scale)
    if source.label.startswith("minor-tertiary-1"):
        shared[4] = (0.0, 1.0 - squared_scale)
    left = list(shared)
    right = list(shared)
    left[1] = (0.5 - scale, 0.5)
    right[1] = (0.5, 0.5 + scale)
    return tuple(left), tuple(right)


def minor_prevalidated_boxes(source: ChartTable) -> tuple[Box, ...]:
    """Return the Cartesian regions covered by a minor's local charts."""

    boxes = [minor_corner_box(source), *minor_midpoint_corner_boxes(source)]
    if source.label.startswith("minor-secondary"):
        boxes.extend((secondary_corner_ratio_away_box(), secondary_origin_corner_box()))
    if source.label.startswith("minor-tertiary-1"):
        boxes.append(tertiary_one_ratio_zero_box())
    return tuple(boxes)


def certify_low_minor_corner(
    source: ChartTable,
    c_upper: Fraction,
    *,
    max_depth: int,
    max_leaves: int,
    power_map: SparseMap | None = None,
) -> None:
    """Certify both projective boundary covers for one low-nome minor."""

    if c_upper > Fraction(1, 200):
        raise ValueError("the audited minor tail requires c_upper <= 1/200")
    axes, labels, upper_axes, lower_axes, envelope_axis = minor_corner_configuration(
        source
    )
    started = time.monotonic()
    if power_map is None:
        power_map = asymptotic_power_map(low_nome_table(source), c_upper)
    corner_map = restricted_minor_corner_map(power_map, upper_axes, lower_axes)
    for selected_axis, label in zip(axes, labels):
        chart = projective_sparse_map(corner_map, axes, selected_axis, order=2)
        if source.label.startswith("minor-secondary") and label == "ratio":
            certify_secondary_minor_ratio_chart(
                source,
                chart,
                max_depth=max_depth,
                max_leaves=max_leaves,
            )
            del chart
            gc.collect()
            continue
        certify_sparse_chart(
            f"{source.label}-minor-corner-{label}",
            chart,
            selected_axis,
            max_depth=max_depth,
            max_leaves=max_leaves,
            envelope_axis=envelope_axis,
            collapse_envelope_early=True,
        )
        del chart
        gc.collect()
    del corner_map
    gc.collect()
    if source.label.startswith("minor-secondary"):
        certify_secondary_corner_ratio_away(
            source,
            power_map,
            max_depth=max_depth,
            max_leaves=max_leaves,
        )
        certify_secondary_origin_corner(
            source,
            power_map,
            max_depth=max_depth,
            max_leaves=max_leaves,
        )
    if source.label.startswith("minor-tertiary-1"):
        certify_tertiary_one_ratio_zero(
            source,
            power_map,
            max_depth=max_depth,
            max_leaves=max_leaves,
        )
    certify_minor_midpoint_corner(
        source,
        power_map,
        c_upper=c_upper,
        max_depth=max_depth,
        max_leaves=max_leaves,
    )
    print(
        f"{source.label}: low-minor local cover PASS in "
        f"{time.monotonic() - started:.2f}s",
        flush=True,
    )


def certify_asymptotic_minor(
    source: ChartTable,
    c_upper: Fraction,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify one complete final minor on a low-nome interval."""

    started = time.monotonic()
    power_map = asymptotic_power_map(low_nome_table(source), c_upper)
    certify_low_minor_corner(
        source,
        c_upper,
        max_depth=max_depth,
        max_leaves=max_leaves,
        power_map=power_map,
    )
    tensor = sparse_map_bernstein_tensor(power_map)
    built = time.monotonic()
    result = certify(
        tensor,
        max_depth=max_depth,
        max_leaves=max_leaves,
        prevalidated_boxes=minor_prevalidated_boxes(source),
    )
    print(
        f"{source.label}-minor-global: shape={tensor.lower.shape}, "
        f"pass={result.passed}, leaves={result.leaves}, depth={result.depth}, "
        f"lower={result.minimum_lower:.3e}, "
        f"chart={time.monotonic() - built:.2f}s",
        flush=True,
    )
    if not result.passed:
        print(
            f"  failure_box={result.failure_box}, index={result.failure_index}",
            flush=True,
        )
        raise SystemExit(1)
    print(
        f"{source.label}: complete low-minor PASS in {time.monotonic() - started:.2f}s",
        flush=True,
    )


def certify_second_determinant_orientation(
    source: ChartTable,
    c_upper: Fraction,
    *,
    max_depth: int,
    max_leaves: int,
) -> None:
    """Certify the exact det1 endpoint hierarchy on the asymptotic nome box."""

    if not source.label.startswith("det-1"):
        raise ValueError("the endpoint hierarchy is specific to det1")
    if c_upper > Fraction(1, 200):
        raise ValueError("the audited det1 tail hierarchy requires c_upper <= 1/200")

    started = time.monotonic()
    power_map = asymptotic_power_map(low_nome_table(source), c_upper)
    endpoint = second_determinant_endpoint_map(power_map)
    endpoint_axes = (0, 1, 2, 4, 5)
    endpoint_weights = (1, 2, 2, 1, 1)
    labels = {0: "nome", 1: "S", 4: "D", 5: "B"}
    for axis, axis_label in labels.items():
        chart = projective_sparse_map(
            endpoint,
            endpoint_axes,
            axis,
            order=4,
            weights=endpoint_weights,
        )
        certify_sparse_chart(
            f"{source.label}-det1-endpoint-{axis_label}",
            chart,
            axis,
            max_depth=max_depth,
            max_leaves=max_leaves,
        )
        del chart
        gc.collect()

    endpoint_main = projective_sparse_map(
        endpoint,
        endpoint_axes,
        2,
        order=4,
        weights=endpoint_weights,
    )
    del endpoint, power_map
    gc.collect()

    transverse_axes = (0, 1, 4, 5)
    transverse_weights = (1, 2, 1, 1)
    transverse_labels = {0: "nome", 1: "S", 4: "D", 5: "B"}
    for axis, axis_label in transverse_labels.items():
        chart = projective_sparse_map(
            endpoint_main,
            transverse_axes,
            axis,
            order=2,
            weights=transverse_weights,
        )
        radial_models = (
            SECOND_DETERMINANT_RADIAL_MODELS[1:]
            if axis == 1
            else SECOND_DETERMINANT_RADIAL_MODELS
        )
        for upper, width, degree in radial_models:
            modeled = centered_sparse_axis_model(
                chart,
                2,
                upper=upper,
                width=width,
                degree=degree,
            )
            certify_sparse_chart(
                f"{source.label}-det1-transverse-{axis_label}-"
                f"P-{upper - width}-{upper}",
                modeled,
                axis,
                max_depth=max_depth,
                max_leaves=max_leaves,
            )
            del modeled
            gc.collect()

        if axis != 1:
            del chart
            gc.collect()
            continue

        low_radial = affine_sparse_map_axis(chart, 2, Fraction(0), Fraction(1, 10))
        final_axes = (0, 2, 4, 5)
        for final_axis, final_label in zip(final_axes, ("nome", "P", "D", "B")):
            final_chart = projective_sparse_map(
                low_radial, final_axes, final_axis, order=2
            )
            for s_upper in (Fraction(1, 2), Fraction(1)):
                s_modeled = centered_sparse_axis_model(
                    final_chart,
                    1,
                    upper=s_upper,
                    width=Fraction(1, 2),
                    degree=8,
                )
                for p_upper in (Fraction(1, 2), Fraction(1)):
                    modeled = centered_sparse_axis_model(
                        s_modeled,
                        2,
                        upper=p_upper,
                        width=Fraction(1, 2),
                        degree=8,
                    )
                    certify_sparse_chart(
                        f"{source.label}-det1-final-{final_label}-"
                        f"S-{s_upper - Fraction(1, 2)}-{s_upper}-"
                        f"P-{p_upper - Fraction(1, 2)}-{p_upper}",
                        modeled,
                        final_axis,
                        max_depth=max_depth,
                        max_leaves=max_leaves,
                    )
                    del modeled
                    gc.collect()
                del s_modeled
                gc.collect()
            del final_chart
            gc.collect()
        del low_radial, chart
        gc.collect()

    print(
        f"{source.label}: det1 endpoint hierarchy PASS in "
        f"{time.monotonic() - started:.2f}s",
        flush=True,
    )


def certify_asymptotic_determinant(
    source: ChartTable,
    c_upper: Fraction,
    *,
    max_depth: int,
    max_leaves: int,
    chart_filter: str | None = None,
) -> None:
    """Run the five exact-order asymptotic determinant charts."""

    started = time.monotonic()
    table = low_nome_table(source)
    normalized = time.monotonic()
    power_map = asymptotic_power_map(table, c_upper)
    expanded = time.monotonic()
    labels = ("c", "main", "one-minus-ratio", "one-minus-a", "abs-b")
    for chart_axis, label in zip((0, 1, 2, 4, 5), labels):
        if chart_filter and label != chart_filter:
            continue
        for chart_label, chart_map in iter_asymptotic_chart_maps(
            power_map, chart_axis, label, c_upper
        ):
            prevalidated_boxes: list[Box] = []
            prevalidated_regions: list[Callable[[Box], bool]] = []
            local_maps: list[
                tuple[str, SparseMap, tuple[tuple[float, float], ...]]
            ] = []
            is_zero_parameter_chart = (
                source.label.startswith("det-0") and chart_label == "one-minus-ratio"
            )
            if chart_label.endswith("orientation-weighted4-one-minus-main"):
                local_maps.extend(orientation_endpoint_maps(chart_map))
            for local_label, local_map, local_box in local_maps:
                local_tensor = sparse_map_bernstein_tensor(local_map)
                local_result = certify(
                    local_tensor,
                    max_depth=max_depth,
                    max_leaves=max_leaves,
                )
                print(
                    f"{source.label}-asymptotic-{chart_label}-{local_label}: "
                    f"shape={local_tensor.lower.shape}, "
                    f"pass={local_result.passed}, "
                    f"leaves={local_result.leaves}, "
                    f"depth={local_result.depth}, "
                    f"lower={local_result.minimum_lower:.3e}",
                    flush=True,
                )
                if not local_result.passed:
                    raise SystemExit(1)
                prevalidated_boxes.append(local_box)
                del local_tensor, local_map
                gc.collect()
            if is_zero_parameter_chart:
                main_arm_cap = (
                    POSITIVE_ZERO_PARAMETER_MAIN_ARM_CAP
                    if source.label.endswith("sign-+1")
                    else NEGATIVE_ZERO_PARAMETER_MAIN_ARM_CAP
                )
                for (
                    corner_label,
                    corner_tensor,
                ) in widened_zero_parameter_corner_tensors(chart_map):
                    corner_result = certify(
                        corner_tensor,
                        max_depth=max_depth,
                        max_leaves=max_leaves,
                    )
                    print(
                        f"{source.label}-asymptotic-{chart_label}-"
                        f"zero-parameter-{corner_label}: "
                        f"shape={corner_tensor.lower.shape}, "
                        f"pass={corner_result.passed}, "
                        f"leaves={corner_result.leaves}, "
                        f"depth={corner_result.depth}, "
                        f"lower={corner_result.minimum_lower:.3e}",
                        flush=True,
                    )
                    if not corner_result.passed:
                        raise SystemExit(1)
                    del corner_tensor
                    gc.collect()
                certify_first_determinant_zero_parameter_main(
                    source,
                    c_upper,
                    max_depth=max_depth,
                    max_leaves=max_leaves,
                    ratio_chart=chart_map,
                )
                certify_first_determinant_zero_parameter_nome_arm(
                    source,
                    chart_map,
                    max_depth=max_depth,
                    max_leaves=max_leaves,
                )
                certify_first_determinant_zero_parameter_main_arm(
                    source,
                    chart_map,
                    main_arm_cap,
                    max_depth=max_depth,
                    max_leaves=max_leaves,
                )
                # One of c, |b|, sqrt(1-u), sqrt(1-v), or sqrt(1-a)
                # is largest, so the four widened charts plus the U chart
                # cover this entire Cartesian corner box.
                prevalidated_boxes.append(zero_parameter_corner_box())
                prevalidated_regions.extend(
                    (
                        inside_zero_parameter_nome_arm,
                        zero_parameter_main_arm_region(main_arm_cap),
                    )
                )
            tensor = sparse_map_bernstein_tensor(chart_map)
            built = time.monotonic()
            result = certify(
                tensor,
                max_depth=max_depth,
                max_leaves=max_leaves,
                prevalidated_boxes=tuple(prevalidated_boxes),
                prevalidated_regions=tuple(prevalidated_regions),
            )
            print(
                f"{source.label}-asymptotic-{chart_label}: "
                f"shape={tensor.lower.shape}, pass={result.passed}, "
                f"leaves={result.leaves}, depth={result.depth}, "
                f"lower={result.minimum_lower:.3e}, "
                f"chart={time.monotonic() - built:.2f}s",
                flush=True,
            )
            if not result.passed:
                print(
                    f"  failure_box={result.failure_box}, index={result.failure_index}",
                    flush=True,
                )
                raise SystemExit(1)
            del tensor, chart_map
            gc.collect()
    print(
        f"{source.label}: normalize={normalized - started:.2f}s, "
        f"asymptotic-map={expanded - normalized:.2f}s, "
        f"total={time.monotonic() - started:.2f}s",
        flush=True,
    )


def arb_low_nome_data(
    midpoint: Fraction, degree: int, *, terms: int = 12
) -> tuple[ArbTaylorJet, ...]:
    """Return center jets for ``(c,K2,S2,L2,D2)``."""

    ctx.dps = 80
    midpoint_ball = arb(midpoint.numerator) / midpoint.denominator
    return _arb_low_nome_data_from_ball(midpoint_ball, degree, terms=terms)


def arb_low_nome_data_box(
    lower: Fraction,
    upper: Fraction,
    degree: int,
    *,
    terms: int = 12,
) -> tuple[ArbTaylorJet, ...]:
    """Return deviation-jet enclosures on a rational nome box."""

    ctx.dps = 80
    lower_ball = arb(lower.numerator) / lower.denominator
    upper_ball = arb(upper.numerator) / upper.denominator
    return _arb_low_nome_data_from_ball(
        lower_ball.union(upper_ball), degree, terms=terms
    )


def _arb_low_nome_data_from_ball(
    c_ball: arb, degree: int, *, terms: int
) -> tuple[ArbTaylorJet, ...]:
    if terms < 2:
        raise ValueError("retain at least two normalized theta terms")
    c = ArbTaylorJet([c_ball, arb(1)] + [arb(0)] * (degree - 1))
    b_series = as_arb_jet(0, degree)
    e_series = as_arb_jet(0, degree)
    for index in range(1, terms + 1):
        b_series += c ** (2 * index * (index + 1) - 4)
        if index >= 2:
            e_series += 2 * c ** (2 * index * index - 8)

    first_omitted = terms + 1
    b_coefficients = list(b_series.coefficients)
    e_coefficients = list(e_series.coefficients)
    for derivative in range(degree + 1):
        b_coefficients[derivative] += arb_positive_series_tail(
            1,
            2 * first_omitted * (first_omitted + 1) - 4,
            derivative,
            c_ball,
            first_omitted,
        )
        e_coefficients[derivative] += arb_positive_series_tail(
            2,
            2 * first_omitted * first_omitted - 8,
            derivative,
            c_ball,
            first_omitted,
        )
    b_series = ArbTaylorJet(b_coefficients)
    e_series = ArbTaylorJet(e_coefficients)

    x = c**2
    a_series = 2 + x**3 * e_series
    theta = 1 + x * a_series
    k2 = (
        4
        * (2 * x * b_series + x**3 * b_series**2 - 2 * a_series - x * a_series**2)
        / theta**2
    )
    s2 = -(2 * a_series + x * a_series**2) / theta**2

    l2_numerator = (
        12 * a_series
        + 8 * x * a_series**2
        + 2 * x**2 * a_series**3
        + 4 * x**2 * e_series
        + 6 * x**3 * e_series * a_series
        + 4 * x**4 * e_series * a_series**2
        + x**5 * e_series * a_series**3
        + 64 * x * b_series
        + 96 * x**3 * b_series**2
        + 64 * x**5 * b_series**3
        + 16 * x**7 * b_series**4
        - 144 * a_series
        - 360 * x * a_series**2
        - 480 * x**2 * a_series**3
        - 360 * x**3 * a_series**4
        - 144 * x**4 * a_series**5
        - 24 * x**5 * a_series**6
    )
    l2 = l2_numerator / (6 * theta**6)

    divided_square_defect = (
        4 * x**2 * e_series
        + x**5 * e_series**2
        - 24 * a_series
        - 60 * x * a_series**2
        - 80 * x**2 * a_series**3
        - 60 * x**3 * a_series**4
        - 24 * x**4 * a_series**5
        - 4 * x**5 * a_series**6
    )
    d2_numerator = (
        8 * x * e_series
        + 48 * divided_square_defect
        + 92 * a_series**3
        - 64 * b_series
        + 83 * x * a_series**4
        + x**2 * (36 * a_series**5 - 96 * b_series**2)
        + 6 * x**3 * a_series**6
        - 64 * x**4 * b_series**3
        - 16 * x**6 * b_series**4
    )
    d2 = d2_numerator / (6 * theta**6)
    return c, k2, s2, l2, d2


def low_factor_jets(
    box: tuple[Fraction, Fraction],
    keys: set[tuple[int, ...]],
    degree: int,
) -> FactorJets:
    midpoint = (box[0] + box[1]) / 2
    center = arb_low_nome_data(midpoint, degree)
    whole = arb_low_nome_data_box(box[0], box[1], degree)
    return factor_jets_from_data(center, whole, keys, degree)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lower", type=parse_fraction, nargs="?")
    parser.add_argument("upper", type=parse_fraction, nargs="?")
    parser.add_argument(
        "--audit-ridge",
        action="store_true",
        help="verify the exact determinant exceptional-divisor hierarchy",
    )
    parser.add_argument(
        "--certify-asymptotic",
        action="store_true",
        help="run the five Taylor/projective determinant remainder charts",
    )
    parser.add_argument(
        "--certify-det1-orientation",
        action="store_true",
        help="run the factorized second-determinant endpoint hierarchy",
    )
    parser.add_argument(
        "--certify-det0-zero-main",
        action="store_true",
        help="run the final first-determinant U-dominant corner hierarchy",
    )
    parser.add_argument(
        "--certify-minor-corners",
        action="store_true",
        help="run the low-nome projective corner covers for the final minors",
    )
    parser.add_argument(
        "--certify-asymptotic-minors",
        action="store_true",
        help="run the integrated low-nome certificates for the final minors",
    )
    parser.add_argument(
        "--asymptotic-chart",
        choices=("c", "main", "one-minus-ratio", "one-minus-a", "abs-b"),
        help="run only one top-level asymptotic chart",
    )
    parser.add_argument("--degree", type=int, default=10)
    parser.add_argument("--sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--chart", help="run only labels containing this text")
    parser.add_argument("--max-depth", type=int, default=40)
    parser.add_argument("--max-leaves", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.audit_ridge:
        audit_determinant_ridge()
        if args.lower is None and args.upper is None:
            return
    if args.lower is None or args.upper is None:
        raise SystemExit("provide both lower and upper, or use --audit-ridge")
    if not 0 <= args.lower < args.upper < 1:
        raise SystemExit("require 0 <= lower < upper < 1")
    if args.degree < 1:
        raise SystemExit("require degree >= 1")
    selected = final_chart_tables(load_records(), args.sign)
    if args.chart:
        selected = [table for table in selected if args.chart in table.label]
    if not selected:
        raise SystemExit("no chart matched")
    if args.certify_det1_orientation:
        if args.lower != 0:
            raise SystemExit("the det1 endpoint charts require lower=0")
        det1_sources = [table for table in selected if table.label.startswith("det-1")]
        if not det1_sources:
            raise SystemExit("no det1 chart matched")
        for source in det1_sources:
            certify_second_determinant_orientation(
                source,
                args.upper,
                max_depth=args.max_depth,
                max_leaves=args.max_leaves,
            )
        return
    if args.certify_det0_zero_main:
        if args.lower != 0:
            raise SystemExit("the det0 zero-parameter charts require lower=0")
        det0_sources = [table for table in selected if table.label.startswith("det-0")]
        if not det0_sources:
            raise SystemExit("no det0 chart matched")
        for source in det0_sources:
            certify_first_determinant_zero_parameter_main(
                source,
                args.upper,
                max_depth=args.max_depth,
                max_leaves=args.max_leaves,
            )
        return
    if args.certify_minor_corners or args.certify_asymptotic_minors:
        if args.lower != 0:
            raise SystemExit("the low-minor charts require lower=0")
        minor_sources = [
            table for table in selected if table.label.startswith("minor-")
        ]
        if not minor_sources:
            raise SystemExit("no minor chart matched")
        certificate = (
            certify_asymptotic_minor
            if args.certify_asymptotic_minors
            else certify_low_minor_corner
        )
        for source in minor_sources:
            certificate(
                source,
                args.upper,
                max_depth=args.max_depth,
                max_leaves=args.max_leaves,
            )
        return
    if args.certify_asymptotic:
        if args.lower != 0:
            raise SystemExit("the asymptotic charts require lower=0")
        for source in selected:
            if not source.label.startswith("det-"):
                continue
            certify_asymptotic_determinant(
                source,
                args.upper,
                max_depth=args.max_depth,
                max_leaves=args.max_leaves,
                chart_filter=args.asymptotic_chart,
            )
        return

    for source in selected:
        started = time.monotonic()
        table = low_nome_table(source)
        normalized = time.monotonic()
        tensor = interval_tensor_for_chart(
            table,
            (args.lower, args.upper),
            args.degree,
            factor_provider=low_factor_jets,
        )
        built = time.monotonic()
        result = certify(tensor, max_depth=args.max_depth, max_leaves=args.max_leaves)
        print(
            f"{table.label}: records={len(table.records)}, shape={tensor.lower.shape}, "
            f"normalize={normalized - started:.2f}s, build={built - normalized:.2f}s, "
            f"pass={result.passed}, leaves={result.leaves}, depth={result.depth}, "
            f"lower={result.minimum_lower:.3e}, total={time.monotonic() - started:.2f}s",
            flush=True,
        )
        if not result.passed:
            print(
                f"  failure_box={result.failure_box}, index={result.failure_index}",
                flush=True,
            )
            raise SystemExit(1)


if __name__ == "__main__":
    main()
