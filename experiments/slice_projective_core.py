#!/usr/bin/env python3
"""Generate the exact polynomial records for the elliptic transfer core.

The derivation is documented in ``proof/slice_core_projective_reduction.md``.
This module regenerates, rather than checks in, the moderately large sparse
record tables used by the projective Bernstein experiments.  Its three chart
levels are exact integer substitutions:

* the main blow-up of ``(P, 1-o) = (0, 0)``;
* the secondary blow-up of ``(v, 1-|b|) = (0, 0)`` for the leading minor;
* the tertiary blow-up of ``(h, 1-a) = (0, 0)`` in the first secondary chart.

Run with ``--force`` to rebuild the cache from the original 4-by-4 core.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from math import comb
from pathlib import Path
import pickle
import time
from typing import TypeAlias

import sympy as sp


Monomial: TypeAlias = tuple[int, ...]
Records: TypeAlias = dict[Monomial, int]
CACHE_VERSION = 1
CACHE_DIRECTORY = Path(__file__).with_name(".cache")


@dataclass(frozen=True, slots=True)
class CoreRecords:
    """Sparse full-envelope records for the two required core minors."""

    minor3: Records
    determinant: Records


def _clean(records: Mapping[Monomial, int | sp.Integer]) -> Records:
    """Combine a sparse integer mapping and remove exact zeros."""

    return {monomial: int(value) for monomial, value in records.items() if value}


def _divide_exactly(polynomial: sp.Poly, divisor: sp.Poly) -> sp.Poly:
    quotient, remainder = sp.div(polynomial, divisor)
    if not remainder.is_zero:
        raise AssertionError(f"nonzero exact-division remainder: {remainder.as_expr()}")
    return quotient


def _base_minor_polynomials() -> dict[str, sp.Poly]:
    """Return defect-normalized minors before the modal substitution."""

    q, a, b, x, y, z, delta = sp.symbols("q a b x y z delta")
    c = q**2
    odd_defect = 1 - a**2
    even_defect = 1 - b**2
    modal = sp.Matrix([[x, -delta * y], [y, z]])
    scaling = sp.diag(1, c)
    upper = scaling * modal * scaling.inv() / q
    lower = q * scaling * modal.T * scaling.inv()
    identity = sp.eye(2)
    odd_block = (
        even_defect * (4 - a**2) * identity
        + odd_defect * (4 * b**2 - 1) * lower.T * lower
    )
    even_block = (
        odd_defect * (4 - b**2) * identity
        + even_defect * (4 * a**2 - 1) * upper.T * upper
    )
    coupling = a * even_defect * upper + b * odd_defect * lower.T
    core = odd_block.row_join(-3 * coupling).col_join(
        (-3 * coupling.T).row_join(even_block)
    )
    variables = (q, a, b, x, y, z, delta)
    output: dict[str, sp.Poly] = {}
    for name, matrix, defect_power in (
        ("minor3", core[:3, :3], 1),
        ("determinant", core, 2),
    ):
        numerator = sp.together(matrix.det(method="berkowitz")).as_numer_denom()[0]
        polynomial = sp.Poly(numerator, *variables)
        for _ in range(defect_power):
            polynomial = _divide_exactly(polynomial, sp.Poly(odd_defect, *variables))
            polynomial = _divide_exactly(polynomial, sp.Poly(even_defect, *variables))
        output[name] = polynomial
    return output


def _modal_polynomial(base: sp.Poly) -> sp.Poly:
    """Substitute the cancellation-free modal matrix and compactify orientation."""

    c, p, r, tangent, k, orientation, a, b = sp.symbols(
        "c p r tangent k orientation a b"
    )
    x_entry = 1 + p * r * tangent
    y_entry = 1 - p * r
    z_entry = p + r * tangent
    normalization = (1 + tangent) * (1 + r**2 * tangent)
    modal_degree = max(sum(monomial[3:6]) for monomial, _ in base.terms()) // 2
    expression = 0
    for monomial, coefficient in base.terms():
        eq, ea, eb, ex, ey, ez, edelta = monomial
        entry_degree = ex + ey + ez
        if eq % 2 or entry_degree % 2 or ey % 2 or edelta > ey:
            raise AssertionError(f"unexpected modal parity: {monomial}")
        expression += (
            coefficient
            * c ** (eq // 2)
            * a**ea
            * b**eb
            * k ** (entry_degree // 2)
            * tangent ** (ey // 2)
            * (p - r) ** edelta
            * y_entry ** (ey - edelta)
            * x_entry**ex
            * z_entry**ez
            * normalization ** (modal_degree - entry_degree // 2)
        )
    tangent_polynomial = sp.Poly(sp.expand(expression), c, p, r, tangent, k, a, b)
    tangent_degree = tangent_polynomial.degree(tangent)
    compactified = sp.cancel(
        tangent_polynomial.as_expr().subs(tangent, orientation / (1 - orientation))
        * (1 - orientation) ** tangent_degree
    )
    if sp.denom(compactified) != 1:
        raise AssertionError("orientation compactification retained a denominator")
    return sp.Poly(sp.expand(compactified), c, p, r, orientation, k, a, b)


def _full_envelope_records(polynomial: sp.Poly) -> Records:
    """Apply ``r=p(s+gamma P)`` and ``gamma=lower+y*width`` exactly."""

    endpoint: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, coefficient in polynomial.terms():
        ec, ep, er, eo, ek, ea, eb = monomial
        coefficient = int(coefficient)
        for gamma_power in range(er + 1):
            p_power = ep + er + 2 * gamma_power
            if p_power % 2:
                raise AssertionError(
                    f"odd p power after envelope substitution: {monomial}"
                )
            for t_power in range(eo + 1):
                key = (
                    p_power // 2,
                    t_power,
                    ea,
                    eb,
                    ec,
                    ek,
                    er - gamma_power,
                    gamma_power,
                )
                endpoint[key] += (
                    coefficient
                    * comb(er, gamma_power)
                    * comb(eo, t_power)
                    * (-1) ** t_power
                )

    full: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, coefficient in _clean(endpoint).items():
        i_p, i_t, ea, eb, ec, ek, es, egamma = monomial
        for envelope_power in range(egamma + 1):
            key = (
                i_p,
                i_t,
                envelope_power,
                ea,
                eb,
                ec,
                ek,
                es,
                egamma - envelope_power,
                envelope_power,
            )
            full[key] += coefficient * comb(egamma, envelope_power)
    return _clean(full)


def generate_records() -> CoreRecords:
    """Regenerate both exact record tables from the transfer-core definition."""

    base = _base_minor_polynomials()
    minor = _full_envelope_records(_modal_polynomial(base["minor3"]))
    determinant = _full_envelope_records(_modal_polynomial(base["determinant"]))
    return CoreRecords(minor, determinant)


def cache_path() -> Path:
    return CACHE_DIRECTORY / f"slice_projective_core_v{CACHE_VERSION}.pkl"


def load_records(*, force: bool = False) -> CoreRecords:
    """Load the deterministic cache, regenerating it when absent or requested."""

    path = cache_path()
    if path.exists() and not force:
        with path.open("rb") as stream:
            payload = pickle.load(stream)
        if (
            not isinstance(payload, dict)
            or payload.get("version") != CACHE_VERSION
            or not isinstance(payload.get("minor3"), dict)
            or not isinstance(payload.get("determinant"), dict)
        ):
            raise TypeError(f"unexpected cache payload in {path}")
        return CoreRecords(payload["minor3"], payload["determinant"])
    records = generate_records()
    CACHE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        pickle.dump(
            {
                "version": CACHE_VERSION,
                "minor3": records.minor3,
                "determinant": records.determinant,
            },
            stream,
            protocol=5,
        )
    return records


def main_chart(records: Records, chart: int, order: int) -> Records:
    """Apply one of the two largest-coordinate charts in ``(P, 1-o)``."""

    if chart not in (0, 1):
        raise ValueError(chart)
    output: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, coefficient in records.items():
        i_p, i_t = monomial[:2]
        if i_p + i_t < order:
            raise AssertionError("claimed main vanishing order is too large")
        ratio_power = i_t if chart == 0 else i_p
        output[(i_p + i_t - order, ratio_power) + monomial[2:]] += coefficient
    return _clean(output)


def signed_secondary_records(records: Records, sign: int) -> Records:
    """Blow up ``(v, 1-|b|)`` in the second main minor chart."""

    if sign not in (-1, 1):
        raise ValueError(sign)
    output: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, coefficient in records.items():
        i_p, i_t, iy, ea, eb = monomial[:5]
        main_power = i_p + i_t - 3
        coefficient *= sign**eb
        for q_power in range(eb + 1):
            key = (
                main_power,
                i_p,
                q_power,
                iy,
                ea,
            ) + monomial[5:]
            output[key] += coefficient * comb(eb, q_power) * (-1) ** q_power
    cleaned = _clean(output)
    if min(monomial[1] + monomial[2] for monomial in cleaned) != 1:
        raise AssertionError("secondary vanishing order is not one")
    return cleaned


def secondary_chart(records: Records, chart: int) -> Records:
    """Apply one of the two largest-coordinate secondary charts."""

    if chart not in (0, 1):
        raise ValueError(chart)
    output: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, coefficient in records.items():
        main_power, v_power, q_power = monomial[:3]
        ratio_power = q_power if chart == 0 else v_power
        key = (
            main_power,
            v_power + q_power - 1,
            ratio_power,
        ) + monomial[3:]
        output[key] += coefficient
    return _clean(output)


def tertiary_records(records: Records) -> Records:
    """Blow up ``(h, 1-a)`` in the first secondary chart."""

    first_chart = secondary_chart(records, 0)
    centered: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, coefficient in first_chart.items():
        main_power, secondary_power, h_power, iy, ea = monomial[:5]
        for z_power in range(ea + 1):
            key = (
                main_power,
                secondary_power,
                h_power,
                z_power,
                iy,
            ) + monomial[5:]
            centered[key] += coefficient * comb(ea, z_power) * (-1) ** z_power
    cleaned = _clean(centered)
    if min(monomial[2] + monomial[3] for monomial in cleaned) != 1:
        raise AssertionError("tertiary vanishing order is not one")

    output: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, coefficient in cleaned.items():
        main_power, secondary_power, h_power, z_power = monomial[:4]
        for chart in (0, 1):
            ratio_power = z_power if chart == 0 else h_power
            key = (
                chart,
                main_power,
                secondary_power,
                h_power + z_power - 1,
                ratio_power,
            ) + monomial[4:]
            output[key] += coefficient
    return _clean(output)


def audit(records: CoreRecords) -> None:
    """Assert the exact chart orders and print stable record counts."""

    expected = {
        "minor3": (49_448, 3),
        "determinant": (197_563, 4),
    }
    for name, table in (
        ("minor3", records.minor3),
        ("determinant", records.determinant),
    ):
        count, order = expected[name]
        actual_order = min(monomial[0] + monomial[1] for monomial in table)
        if len(table) != count or actual_order != order:
            raise AssertionError(
                f"{name}: got ({len(table)}, {actual_order}), expected ({count}, {order})"
            )
        chart_counts = tuple(len(main_chart(table, chart, order)) for chart in (0, 1))
        print(f"{name}: {len(table)} records, order {order}, charts {chart_counts}")

    for sign in (1, -1):
        secondary = signed_secondary_records(records.minor3, sign)
        secondary_counts = tuple(
            len(secondary_chart(secondary, chart)) for chart in (0, 1)
        )
        tertiary = tertiary_records(secondary)
        tertiary_counts = tuple(
            sum(monomial[0] == chart for monomial in tertiary) for chart in (0, 1)
        )
        print(
            f"minor3 sign {sign:+d}: secondary {len(secondary)} -> "
            f"{secondary_counts}; tertiary {tertiary_counts}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true", help="regenerate the exact cache"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.monotonic()
    records = load_records(force=args.force)
    audit(records)
    print(f"projective-core exact audit: PASS ({time.monotonic() - started:.2f}s)")


if __name__ == "__main__":
    main()
