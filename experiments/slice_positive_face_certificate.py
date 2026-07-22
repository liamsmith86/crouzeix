#!/usr/bin/env python3
"""Certify positive determinant charts in deficit-centered coordinates.

On the face ``a=b=1`` both determinant charts factor positively, and their
linear coefficients in ``A=1-a`` and ``B=1-b`` are nonnegative.  This script
bounds only the terms of total ``(A,B)`` degree at least two.  It preserves
the common nome displacement and converts the surviving ``X,R,Y`` variables
to Bernstein form before taking absolute values.

One independently checkable estimate has the form

    determinant >= face_lower - remainder_bound * epsilon**2

when ``0 <= A,B <= epsilon``.  The final certificate is sharper: it converts
the complete deficit-centered determinant to Bernstein form, so the positive
face and linear terms remain correlated with the higher remainder.  Optional
``X,Y`` localization can strengthen the bound, but the tail certificate uses
the full spatial cube and ``epsilon=1``, hence covers the complete ``a,b``
square rather than only a thin face slab.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
import gc
from itertools import product
from math import comb
import time
from typing import TypeAlias

from flint import arb
import numpy as np

from slice_positive_face_audit import ExactPolynomial, collected_deficit
from slice_projective_core import load_records
from slice_projective_interval_certificate import (
    ArbTaylorJet,
    Box,
    CertificateResult,
    ChartTable,
    IntervalTensor,
    arb_theta_data,
    arb_theta_data_box,
    as_arb_jet,
    certify,
    factor_jets_from_data,
    final_chart_tables,
    parse_fraction,
)
from slice_odd_compact_certificate import round_down, round_up


SpatialIndex: TypeAlias = tuple[int, ...]
ArbPolynomial: TypeAlias = dict[SpatialIndex, arb]


@dataclass(frozen=True, slots=True)
class FaceSlabResult:
    """Rigorous bounds defining one independently certified face slab."""

    table: ChartTable
    face_lower: arb
    remainder_bound: arb | None
    coarse_residual_lower: arb | None
    centered_root_lower: float
    centered_certificate: CertificateResult
    recentered_boxes: tuple[Box, ...]


def rational_arb(value: Fraction) -> arb:
    """Return an Arb point containing one exact rational."""

    return arb(value.numerator) / value.denominator


def _bernstein_evaluate(
    controls: Mapping[SpatialIndex, arb],
    dimensions: Sequence[int],
    point: Sequence[Fraction],
) -> arb:
    """Evaluate a small Arb Bernstein map for machinery audits."""

    output = arb(0)
    for index, value in controls.items():
        basis = Fraction(1)
        for axis, bernstein_index in enumerate(index):
            degree = dimensions[axis] - 1
            coordinate = point[axis]
            basis *= (
                comb(degree, bernstein_index)
                * coordinate**bernstein_index
                * (1 - coordinate) ** (degree - bernstein_index)
            )
        output += value * rational_arb(basis)
    return output


def audit_bernstein_machinery() -> None:
    """Independently check affine restriction and both Bernstein conversions."""

    exact_power = {
        (0, 0, 0): Fraction(2),
        (1, 0, 0): Fraction(-3),
        (2, 1, 0): Fraction(5),
        (0, 1, 2): Fraction(-7),
    }
    dimensions = (3, 2, 3)
    domains = (
        (Fraction(1, 4), Fraction(3, 4)),
        (Fraction(0), Fraction(1)),
        (Fraction(2, 5), Fraction(1)),
    )
    local_point = (Fraction(2, 7), Fraction(3, 5), Fraction(4, 9))
    original_point = tuple(
        lower + (upper - lower) * coordinate
        for (lower, upper), coordinate in zip(domains, local_point)
    )
    expected = sum(
        (
            coefficient
            * original_point[0] ** index[0]
            * original_point[1] ** index[1]
            * original_point[2] ** index[2]
            for index, coefficient in exact_power.items()
        ),
        Fraction(0),
    )
    controls = power_map_to_bernstein(
        {index: rational_arb(value) for index, value in exact_power.items()},
        dimensions,
        domains,
    )
    if not _bernstein_evaluate(controls, dimensions, local_point).contains(
        rational_arb(expected)
    ):
        raise AssertionError("spatial Bernstein audit failed")

    lower, upper = Fraction(1, 5), Fraction(3, 5)
    midpoint = (lower + upper) / 2
    center_coefficients = (
        1 - 2 * midpoint + 3 * midpoint**2 - midpoint**3,
        -2 + 6 * midpoint - 3 * midpoint**2,
        3 - 3 * midpoint,
        Fraction(-1),
    )
    c_controls = c_bernstein_controls(
        [rational_arb(value) for value in center_coefficients],
        arb(0),
        (upper - lower) / 2,
    )
    local_c = Fraction(5, 11)
    c = lower + (upper - lower) * local_c
    expected_c = 1 - 2 * c + 3 * c**2 - c**3
    actual_c = sum(
        (
            control
            * rational_arb(
                Fraction(comb(3, index)) * local_c**index * (1 - local_c) ** (3 - index)
            )
            for index, control in enumerate(c_controls)
        ),
        arb(0),
    )
    if not actual_c.contains(rational_arb(expected_c)):
        raise AssertionError("nome Bernstein audit failed")


def interval_ball(lower: Fraction, upper: Fraction) -> arb:
    """Return the Arb hull of two exact rational endpoints."""

    return rational_arb(lower).union(rational_arb(upper))


def power_map_to_bernstein(
    values: Mapping[SpatialIndex, arb],
    dimensions: Sequence[int],
    domains: Sequence[tuple[Fraction, Fraction]] | None = None,
) -> ArbPolynomial:
    """Restrict a sparse power map and convert it to Bernstein form."""

    output = dict(values)
    if domains is None:
        domains = [(Fraction(0), Fraction(1)) for _ in dimensions]
    if len(domains) != len(dimensions):
        raise ValueError("domain dimension does not match the polynomial")
    for axis, (lower, upper) in enumerate(domains):
        if not 0 <= lower < upper <= 1:
            raise ValueError("Bernstein domains must be nonempty unit intervals")
        if lower == 0 and upper == 1:
            continue
        restricted: ArbPolynomial = {}
        width = upper - lower
        for index, value in output.items():
            power = index[axis]
            for local_power in range(power + 1):
                target = (*index[:axis], local_power, *index[axis + 1 :])
                scalar = (
                    comb(power, local_power)
                    * lower ** (power - local_power)
                    * width**local_power
                )
                restricted[target] = restricted.get(
                    target, arb(0)
                ) + value * rational_arb(scalar)
        output = restricted

    for axis, size in enumerate(dimensions):
        degree = size - 1
        converted: ArbPolynomial = {}
        for index, value in output.items():
            power = index[axis]
            for bernstein_index in range(power, size):
                target = (
                    *index[:axis],
                    bernstein_index,
                    *index[axis + 1 :],
                )
                weight = arb(comb(bernstein_index, power)) / comb(degree, power)
                converted[target] = converted.get(target, arb(0)) + value * weight
        output = converted
    return output


def c_bernstein_controls(
    center_coefficients: Sequence[arb],
    remainder_coefficient: arb,
    radius: Fraction,
) -> tuple[arb, ...]:
    """Enclose a Taylor model using the nome as one shared Bernstein axis."""

    degree = len(center_coefficients)
    if degree < 1:
        raise ValueError("at least one retained Taylor coefficient is required")
    power_coefficients = [arb(0) for _ in range(degree)]
    for order, coefficient in enumerate(center_coefficients):
        for exponent in range(order + 1):
            scalar = (
                radius**order
                * comb(order, exponent)
                * 2**exponent
                * (-1) ** (order - exponent)
            )
            power_coefficients[exponent] += coefficient * rational_arb(scalar)

    bernstein_degree = degree - 1
    controls: list[arb] = []
    for index in range(degree):
        controls.append(
            sum(
                (
                    power_coefficients[exponent]
                    * arb(comb(index, exponent))
                    / comb(bernstein_degree, exponent)
                    for exponent in range(index + 1)
                ),
                arb(0),
            )
        )

    remainder_radius = radius**degree
    if degree % 2:
        displacement = interval_ball(-remainder_radius, remainder_radius)
    else:
        displacement = interval_ball(Fraction(0), remainder_radius)
    error = remainder_coefficient * displacement
    return tuple(control + error for control in controls)


def jet_range(
    center: ArbTaylorJet,
    whole: ArbTaylorJet,
    box: tuple[Fraction, Fraction],
) -> arb:
    """Return the convex hull of a rigorous one-dimensional Taylor model."""

    if center.degree != whole.degree:
        raise ValueError("Taylor degrees do not match")
    controls = c_bernstein_controls(
        center.coefficients[:-1],
        whole.coefficients[-1],
        (box[1] - box[0]) / 2,
    )
    enclosure = controls[0]
    for control in controls[1:]:
        enclosure = enclosure.union(control)
    return enclosure


def face_lower_bound(
    chart: int,
    box: tuple[Fraction, Fraction],
    degree: int,
    x_minimum: Fraction,
    y_minimum: Fraction,
) -> arb:
    """Certify the factored face lower bound for one determinant chart."""

    midpoint = (box[0] + box[1]) / 2
    center_data = arb_theta_data(midpoint, degree)
    whole_data = arb_theta_data_box(*box, degree)
    x = as_arb_jet(x_minimum, degree)
    y = as_arb_jet(y_minimum, degree)

    def candidates(data: tuple[ArbTaylorJet, ...]) -> tuple[ArbTaylorJet, ...]:
        c, k, s, gamma_lower, gamma_width = data
        if chart == 1:
            return (81 * c**9 * (1 - k) ** 4,)
        q_minimum = gamma_lower + y * gamma_width
        common = 81 * c**9 * (1 - k) ** 2
        at_x_minimum = common * (1 - k * x) ** 2 * (s + q_minimum * x) ** 8
        at_one = common * (1 - k) ** 2 * (s + q_minimum) ** 8
        return at_x_minimum, at_one

    center_candidates = candidates(center_data)
    whole_candidates = candidates(whole_data)
    lower_bounds = [
        arb(jet_range(center, whole, box).lower())
        for center, whole in zip(center_candidates, whole_candidates)
    ]
    enclosure = lower_bounds[0]
    for lower_bound in lower_bounds[1:]:
        enclosure = enclosure.union(lower_bound)
    return arb(enclosure.lower())


def higher_deficit_maps(table: ChartTable) -> dict[tuple[int, int], ExactPolynomial]:
    """Collect all coefficients of total deficit degree at least two."""

    return {
        (a_power, b_power): collected_deficit(table, a_power, b_power)
        for a_power in range(5)
        for b_power in range(5)
        if a_power + b_power >= 2
    }


def all_deficit_maps(table: ChartTable) -> dict[tuple[int, int], ExactPolynomial]:
    """Collect the complete determinant in the two deficit variables."""

    return {
        (a_power, b_power): collected_deficit(table, a_power, b_power)
        for a_power in range(5)
        for b_power in range(5)
    }


def scaled_deficit_records(
    table: ChartTable, epsilon: Fraction
) -> dict[tuple[int, ...], Fraction]:
    """Return records after ``A=epsilon*u`` and ``B=epsilon*v``."""

    output: dict[tuple[int, ...], Fraction] = {}
    for (a_power, b_power), records in all_deficit_maps(table).items():
        deficit_scale = epsilon ** (a_power + b_power)
        for monomial, coefficient in records.items():
            index = (
                *monomial[:3],
                a_power,
                b_power,
                *monomial[3:],
            )
            output[index] = output.get(index, Fraction(0)) + coefficient * deficit_scale
    return {index: coefficient for index, coefficient in output.items() if coefficient}


def correlated_polynomial_tensor(
    records: Mapping[tuple[int, ...], Fraction],
    spatial_dimensions: Sequence[int],
    spatial_domains: Sequence[tuple[Fraction, Fraction]],
    box: tuple[Fraction, Fraction],
    degree: int,
) -> IntervalTensor:
    """Build a directed tensor after sharing the factored nome displacement."""

    spatial_count = len(spatial_dimensions)
    factor_keys = {monomial[spatial_count:] for monomial in records}
    midpoint = (box[0] + box[1]) / 2
    center_data = arb_theta_data(midpoint, degree)
    whole_data = arb_theta_data_box(*box, degree)

    factors = factor_jets_from_data(center_data, whole_data, factor_keys, degree)
    center_maps: list[ArbPolynomial] = [dict() for _ in range(degree)]
    remainder_map: ArbPolynomial = {}
    for monomial, rational_coefficient in records.items():
        spatial_index = monomial[:spatial_count]
        center_factor, whole_factor = factors[monomial[spatial_count:]]
        coefficient = rational_arb(rational_coefficient)
        for order in range(degree):
            value = coefficient * center_factor.coefficients[order]
            center_maps[order][spatial_index] = (
                center_maps[order].get(spatial_index, arb(0)) + value
            )
        value = coefficient * whole_factor.coefficients[degree]
        remainder_map[spatial_index] = remainder_map.get(spatial_index, arb(0)) + value

    center_bernstein = [
        power_map_to_bernstein(values, spatial_dimensions, spatial_domains)
        for values in center_maps
    ]
    remainder_bernstein = power_map_to_bernstein(
        remainder_map, spatial_dimensions, spatial_domains
    )
    radius = (box[1] - box[0]) / 2
    shape = (degree, *spatial_dimensions)
    lower = np.empty(shape)
    upper = np.empty(shape)
    for spatial_index in product(*(range(size) for size in spatial_dimensions)):
        controls = c_bernstein_controls(
            [values.get(spatial_index, arb(0)) for values in center_bernstein],
            remainder_bernstein.get(spatial_index, arb(0)),
            radius,
        )
        for c_index, control in enumerate(controls):
            index = (c_index, *spatial_index)
            lower[index] = round_down(float(control.lower()))
            upper[index] = round_up(float(control.upper()))
    return IntervalTensor(lower, upper)


def determinant_slab_certificate(
    table: ChartTable,
    box: tuple[Fraction, Fraction],
    degree: int,
    x_minimum: Fraction,
    y_minimum: Fraction,
    epsilon: Fraction,
    max_depth: int,
    max_leaves: int,
    recenter_depth: int,
    max_recentered_boxes: int,
) -> tuple[float, CertificateResult, tuple[Box, ...]]:
    """Certify the complete centered face slab by directed subdivision."""

    original_dimensions = tuple(
        max(monomial[axis] for monomial in table.records) + 1 for axis in range(3)
    )
    dimensions = (*original_dimensions, 5, 5)
    domains = (
        (x_minimum, Fraction(1)),
        (Fraction(0), Fraction(1)),
        (y_minimum, Fraction(1)),
        (Fraction(0), Fraction(1)),
        (Fraction(0), Fraction(1)),
    )
    records = scaled_deficit_records(table, epsilon)
    tensor = correlated_polynomial_tensor(
        records,
        dimensions,
        domains,
        box,
        degree,
    )
    root_lower = float(np.min(tensor.lower))
    recentered_boxes: list[Box] = []
    while True:
        result = certify(
            tensor,
            max_depth=max_depth,
            max_leaves=max_leaves,
            prevalidated_boxes=tuple(recentered_boxes),
        )
        if result.passed:
            return root_lower, result, tuple(recentered_boxes)
        if len(recentered_boxes) >= max_recentered_boxes:
            return root_lower, result, tuple(recentered_boxes)
        if result.failure_box is None:
            raise AssertionError("failed certificate omitted its failure box")

        recentered = dyadic_ancestor_box(result.failure_box, recenter_depth)
        local_c_box = restrict_rational_interval(box, recentered[0])
        local_spatial_domains = tuple(
            restrict_rational_interval(domain, local)
            for domain, local in zip(domains, recentered[1:])
        )
        local_tensor = correlated_polynomial_tensor(
            records,
            dimensions,
            local_spatial_domains,
            local_c_box,
            degree,
        )
        local_result = certify(
            local_tensor,
            max_depth=max_depth,
            max_leaves=max_leaves,
        )
        print(
            f"  RECENTER {table.label}: index={len(recentered_boxes) + 1}, "
            f"box={recentered}, root={float(np.min(local_tensor.lower)):.3e}, "
            f"pass={local_result.passed}, leaves={local_result.leaves}",
            flush=True,
        )
        if not local_result.passed:
            return root_lower, local_result, tuple(recentered_boxes)
        recentered_boxes.append(recentered)
        del local_tensor
        gc.collect()


def dyadic_ancestor_box(box: Box, maximum_depth: int) -> Box:
    """Enlarge a dyadic subdivision leaf to a bounded-depth ancestor."""

    if maximum_depth < 0:
        raise ValueError("recenter depth must be nonnegative")
    output: list[tuple[float, float]] = []
    for lower_float, upper_float in box:
        lower = Fraction.from_float(lower_float)
        upper = Fraction.from_float(upper_float)
        width = upper - lower
        if (
            width <= 0
            or width.numerator != 1
            or width.denominator & (width.denominator - 1)
        ):
            raise ValueError("certificate subdivision box is not dyadic")
        depth = width.denominator.bit_length() - 1
        ancestor_depth = min(depth, maximum_depth)
        denominator = 2**ancestor_depth
        midpoint = (lower + upper) / 2
        index = min(denominator - 1, int(midpoint * denominator))
        ancestor_lower = Fraction(index, denominator)
        ancestor_upper = Fraction(index + 1, denominator)
        if not ancestor_lower <= lower < upper <= ancestor_upper:
            raise AssertionError("computed dyadic ancestor does not contain its leaf")
        output.append((exact_binary64(ancestor_lower), exact_binary64(ancestor_upper)))
    return tuple(output)


def restrict_rational_interval(
    interval: tuple[Fraction, Fraction], unit_subinterval: tuple[float, float]
) -> tuple[Fraction, Fraction]:
    """Map an exact binary64 unit subinterval into a rational interval."""

    lower, upper = interval
    unit_lower = Fraction.from_float(unit_subinterval[0])
    unit_upper = Fraction.from_float(unit_subinterval[1])
    width = upper - lower
    return lower + width * unit_lower, lower + width * unit_upper


def audit_recenter_machinery() -> None:
    """Check dyadic ancestry and exact affine transport of a failure leaf."""

    leaf: Box = (
        (0.8125, 0.875),
        (0.0, 0.03125),
    )
    ancestor = dyadic_ancestor_box(leaf, 2)
    if ancestor != ((0.75, 1.0), (0.0, 0.25)):
        raise AssertionError("dyadic ancestor audit failed")
    restricted = restrict_rational_interval(
        (Fraction(1, 3), Fraction(2, 3)), ancestor[0]
    )
    if restricted != (Fraction(7, 12), Fraction(2, 3)):
        raise AssertionError("rational recenter transport audit failed")


def remainder_bound(
    table: ChartTable,
    box: tuple[Fraction, Fraction],
    degree: int,
    x_minimum: Fraction,
    y_minimum: Fraction,
) -> arb:
    """Bound the correlated quadratic-and-higher deficit coefficient sum."""

    deficit_maps = higher_deficit_maps(table)
    factor_keys = {
        monomial[3:] for records in deficit_maps.values() for monomial in records
    }
    midpoint = (box[0] + box[1]) / 2
    center_data = arb_theta_data(midpoint, degree)
    whole_data = arb_theta_data_box(*box, degree)

    factors = factor_jets_from_data(center_data, whole_data, factor_keys, degree)
    dimensions = tuple(
        max(monomial[axis] for monomial in table.records) + 1 for axis in range(3)
    )
    spatial_domains = (
        (x_minimum, Fraction(1)),
        (Fraction(0), Fraction(1)),
        (y_minimum, Fraction(1)),
    )
    radius = (box[1] - box[0]) / 2
    total = arb(0)
    for records in deficit_maps.values():
        center_maps: list[ArbPolynomial] = [dict() for _ in range(degree)]
        remainder_map: ArbPolynomial = {}
        for monomial, integer_coefficient in records.items():
            spatial_index = monomial[:3]
            center_factor, whole_factor = factors[monomial[3:]]
            for order in range(degree):
                value = integer_coefficient * center_factor.coefficients[order]
                center_maps[order][spatial_index] = (
                    center_maps[order].get(spatial_index, arb(0)) + value
                )
            value = integer_coefficient * whole_factor.coefficients[degree]
            remainder_map[spatial_index] = (
                remainder_map.get(spatial_index, arb(0)) + value
            )

        center_bernstein = [
            power_map_to_bernstein(values, dimensions, spatial_domains)
            for values in center_maps
        ]
        remainder_bernstein = power_map_to_bernstein(
            remainder_map, dimensions, spatial_domains
        )
        coefficient_bound = arb(0)
        for spatial_index in product(*(range(size) for size in dimensions)):
            controls = c_bernstein_controls(
                [values.get(spatial_index, arb(0)) for values in center_bernstein],
                remainder_bernstein.get(spatial_index, arb(0)),
                radius,
            )
            for control in controls:
                coefficient_bound = arb(
                    coefficient_bound.union(control.abs_upper()).upper()
                )
        total += coefficient_bound
    return total


def face_slab_result(
    table: ChartTable,
    chart: int,
    box: tuple[Fraction, Fraction],
    degree: int,
    x_minimum: Fraction,
    y_minimum: Fraction,
    epsilon: Fraction,
    max_depth: int,
    max_leaves: int,
    recenter_depth: int = 1,
    max_recentered_boxes: int = 64,
    audit_remainder: bool = False,
) -> FaceSlabResult:
    """Build and verify the slab estimate for one chart."""

    face = face_lower_bound(chart, box, degree, x_minimum, y_minimum)
    remainder = None
    coarse_residual = None
    if audit_remainder:
        remainder = remainder_bound(table, box, degree, x_minimum, y_minimum)
        coarse_residual = face - remainder * rational_arb(epsilon**2)
    centered_root_lower, centered_certificate, recentered_boxes = (
        determinant_slab_certificate(
            table,
            box,
            degree,
            x_minimum,
            y_minimum,
            epsilon,
            max_depth,
            max_leaves,
            recenter_depth,
            max_recentered_boxes,
        )
    )
    if not centered_certificate.passed:
        raise ValueError(
            f"{table.label}: centered face slab failed: "
            f"lower={centered_certificate.minimum_lower}, "
            f"box={centered_certificate.failure_box}"
        )
    return FaceSlabResult(
        table,
        face,
        remainder,
        coarse_residual,
        centered_root_lower,
        centered_certificate,
        recentered_boxes,
    )


def exact_binary64(value: Fraction) -> float:
    """Convert a dyadic rational to binary64 without changing its value."""

    converted = float(value)
    if Fraction.from_float(converted) != value:
        raise ValueError(f"region threshold is not exactly binary64: {value}")
    return converted


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lower", type=parse_fraction)
    parser.add_argument("upper", type=parse_fraction)
    parser.add_argument("--degree", type=int, default=10)
    parser.add_argument("--epsilon-power", type=int, default=0)
    parser.add_argument("--x-minimum", type=parse_fraction, default=Fraction(0))
    parser.add_argument("--y-minimum", type=parse_fraction, default=Fraction(0))
    parser.add_argument("--chart", type=int, choices=(0, 1))
    parser.add_argument("--audit-remainder", action="store_true")
    parser.add_argument("--max-depth", type=int, default=40)
    parser.add_argument("--max-leaves", type=int, default=100_000)
    parser.add_argument("--recenter-depth", type=int, default=1)
    parser.add_argument("--max-recentered-boxes", type=int, default=64)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 < args.lower < args.upper < 1:
        raise SystemExit("require 0 < lower < upper < 1")
    if args.degree < 1:
        raise SystemExit("require degree >= 1")
    if args.epsilon_power < 0:
        raise SystemExit("require epsilon power >= 0")
    if args.max_depth < 0 or args.max_leaves < 1:
        raise SystemExit("require nonnegative depth and positive leaf limit")
    if args.recenter_depth < 0 or args.max_recentered_boxes < 0:
        raise SystemExit("require nonnegative recenter limits")
    if not 0 <= args.x_minimum < 1 or not 0 <= args.y_minimum < 1:
        raise SystemExit("require spatial thresholds in [0,1)")
    audit_bernstein_machinery()
    audit_recenter_machinery()

    box = args.lower, args.upper
    epsilon = Fraction(1, 2**args.epsilon_power)
    tables = final_chart_tables(load_records(), 1)[:2]
    if args.chart is not None:
        tables = tables[args.chart : args.chart + 1]

    for table in tables:
        chart = int(table.label.split("-")[1])
        started = time.monotonic()
        result = face_slab_result(
            table,
            chart,
            box,
            args.degree,
            args.x_minimum,
            args.y_minimum,
            epsilon,
            args.max_depth,
            args.max_leaves,
            recenter_depth=args.recenter_depth,
            max_recentered_boxes=args.max_recentered_boxes,
            audit_remainder=args.audit_remainder,
        )
        remainder_text = (
            f"remainder_bound={result.remainder_bound}, "
            f"coarse_residual_lower={result.coarse_residual_lower}, "
            if result.remainder_bound is not None
            else ""
        )
        print(
            f"PASS {table.label}: face_lower={result.face_lower}, "
            f"epsilon=2^-{args.epsilon_power}, "
            f"x_minimum={args.x_minimum}, y_minimum={args.y_minimum}, "
            f"{remainder_text}"
            f"centered_root_lower={result.centered_root_lower:.3e}, "
            f"centered_leaves={result.centered_certificate.leaves}, "
            f"centered_depth={result.centered_certificate.depth}, "
            f"recentered={len(result.recentered_boxes)}, "
            f"build={time.monotonic() - started:.2f}s",
            flush=True,
        )


if __name__ == "__main__":
    main()
