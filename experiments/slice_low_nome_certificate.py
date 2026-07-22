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
from math import comb
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
    factor_jets_from_data,
    final_chart_tables,
    interval_tensor_for_chart,
    parse_fraction,
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
            print(f"{source.label}: exact c^2 normalization PASS", flush=True)


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


def sparse_map_bernstein_tensor(power_map: SparseMap) -> IntervalTensor:
    """Convert one exact/Arb sparse map to a directed Bernstein tensor."""

    dimensions = tuple(max(key[axis] for key in power_map) + 1 for axis in range(6))
    interval_map = {
        key: value if isinstance(value, arb) else arb(value)
        for key, value in power_map.items()
    }
    return power_to_bernstein_tensor(arb_map_to_tensor(interval_map, dimensions))


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


def widened_zero_parameter_corner_tensors(
    power_map: SparseMap,
) -> Iterator[tuple[str, IntervalTensor, int]]:
    """Build four overlapping weighted charts around the ``a=0`` corner."""

    scale = Fraction(1, 4)
    squared_scale = scale**2
    restricted = power_map
    for axis in (0, 5):
        restricted = affine_sparse_map_axis(restricted, axis, Fraction(0), scale)
    for axis in (1, 2, 4):
        restricted = affine_sparse_map_axis(
            restricted, axis, Fraction(1), -squared_scale
        )

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
        yield labels[selected_axis], sparse_map_bernstein_tensor(chart), selected_axis


def widened_zero_parameter_region_contains(box: Box, selected_axis: int) -> bool:
    """Test containment in one factor-two weighted corner chart."""

    scale = Fraction(1, 4)
    ratio_bound = 2
    coordinates = tuple(
        tuple(Fraction.from_float(endpoint) for endpoint in interval)
        for interval in box
    )
    nome_upper = coordinates[0][1]
    abs_b_upper = coordinates[5][1]

    if selected_axis in (0, 5):
        selected_lower = coordinates[selected_axis][0]
        selected_upper = coordinates[selected_axis][1]
        if selected_upper > scale:
            return False
        other_linear_upper = abs_b_upper if selected_axis == 0 else nome_upper
        if other_linear_upper > ratio_bound * selected_lower:
            return False
        squared_selected = selected_lower**2
        return all(
            1 - coordinates[axis][0] <= ratio_bound * squared_selected
            for axis in (1, 2, 4)
        )

    selected_defect_lower = 1 - coordinates[selected_axis][1]
    selected_defect_upper = 1 - coordinates[selected_axis][0]
    if selected_defect_upper > scale**2:
        return False
    if nome_upper**2 > ratio_bound**2 * selected_defect_lower:
        return False
    if abs_b_upper**2 > ratio_bound**2 * selected_defect_lower:
        return False
    return all(
        1 - coordinates[axis][0] <= ratio_bound * selected_defect_lower
        for axis in (1, 2, 4)
        if axis != selected_axis
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
            if source.label == "det-0-sign-+1" and chart_label == "one-minus-ratio":
                for (
                    corner_label,
                    corner_tensor,
                    selected_axis,
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
                    prevalidated_regions.append(
                        lambda box, axis=selected_axis: (
                            widened_zero_parameter_region_contains(box, axis)
                        )
                    )
                    del corner_tensor
                    gc.collect()
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
