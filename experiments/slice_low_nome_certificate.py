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
from collections import defaultdict
from fractions import Fraction
from math import comb
import time

from flint import arb, ctx

from slice_projective_core import load_records
from slice_projective_interval_certificate import (
    ArbTaylorJet,
    ChartTable,
    FactorJets,
    arb_positive_series_tail,
    as_arb_jet,
    certify,
    factor_jets_from_data,
    final_chart_tables,
    interval_tensor_for_chart,
    parse_fraction,
)


LOW_CONSTANTS = (4, 1, 4, 32)
LOW_BASE_POWERS = (1, 1, 0, 2, 4)
LOW_DEVIATION_SERIES = (
    (-16, 0, 56, 0, -160, 0, 404),
    (-4, 0, 12, 0, -32, 0, 76),
    (-44, 0, 304, 0, -1644, 0, 7544),
    (-272, 0, 1568, 0, -7376, 0, 30_320),
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
    table: ChartTable, degree: int = 6
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
            print(
                f"{source.label}: c^9 normalization, first/second exceptional "
                "forms and final ridge PASS",
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
