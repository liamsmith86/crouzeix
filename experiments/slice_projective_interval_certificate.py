#!/usr/bin/env python3
"""Directed interval probes for the final elliptic transfer-core charts.

This is the continuous-nome companion to ``slice_projective_core.py``.  It
keeps the Taylor displacement in the nome as a shared Bernstein coordinate;
enclosing Taylor orders independently loses the cancellations on sharp faces.

The current command certifies one rational nome box at a time.  It is not yet
the finite global cover claimed by the open transfer theorem, so its output is
reported as a local certificate only.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from math import comb
import math
import time
from typing import Callable, TypeAlias

import numpy as np
from numpy.typing import NDArray

try:
    from flint import arb, ctx
except ImportError as error:  # pragma: no cover - incomplete environment only
    raise SystemExit(
        "python-flint is required for rigorous center-coefficient cancellation"
    ) from error

from slice_odd_compact_certificate import (
    Interval,
    round_down,
    round_up,
)
from slice_projective_core import (
    CoreRecords,
    Records,
    load_records,
    main_chart,
    secondary_chart,
    signed_secondary_records,
    tertiary_records,
)


FloatArray: TypeAlias = NDArray[np.float64]
Box: TypeAlias = tuple[tuple[float, float], ...]
FactorJets: TypeAlias = dict[tuple[int, ...], tuple["ArbTaylorJet", "ArbTaylorJet"]]
FactorJetProvider: TypeAlias = Callable[
    [tuple[Fraction, Fraction], set[tuple[int, ...]], int], FactorJets
]

# A full finite-difference scan gives better subdivision choices on the compact
# certificate tensors.  Above this size its temporary arrays dominate memory,
# so the large weighted corner charts use a local score at the worst control
# point instead.
FULL_VARIATION_SCORE_MAX_SIZE = 5_000_000


@dataclass(frozen=True, slots=True)
class IntervalTensor:
    lower: FloatArray
    upper: FloatArray


@dataclass(frozen=True, slots=True)
class ChartTable:
    label: str
    records: Records


@dataclass(frozen=True, slots=True)
class CertificateResult:
    passed: bool
    leaves: int
    depth: int
    minimum_lower: float
    failure_box: Box | None = None
    failure_index: tuple[int, ...] | None = None


class ArbTaylorJet:
    """Truncated Taylor series whose coefficients are rigorous Arb balls."""

    __slots__ = ("coefficients",)

    def __init__(self, coefficients: tuple[arb, ...] | list[arb]) -> None:
        self.coefficients = tuple(coefficients)

    @property
    def degree(self) -> int:
        return len(self.coefficients) - 1

    def __add__(self, other: object) -> "ArbTaylorJet":
        rhs = as_arb_jet(other, self.degree)
        return ArbTaylorJet(
            [left + right for left, right in zip(self.coefficients, rhs.coefficients)]
        )

    __radd__ = __add__

    def __neg__(self) -> "ArbTaylorJet":
        return ArbTaylorJet([-coefficient for coefficient in self.coefficients])

    def __sub__(self, other: object) -> "ArbTaylorJet":
        return self + (-as_arb_jet(other, self.degree))

    def __rsub__(self, other: object) -> "ArbTaylorJet":
        return as_arb_jet(other, self.degree) - self

    def __mul__(self, other: object) -> "ArbTaylorJet":
        rhs = as_arb_jet(other, self.degree)
        return ArbTaylorJet(
            [
                sum(
                    (
                        self.coefficients[left] * rhs.coefficients[order - left]
                        for left in range(order + 1)
                    ),
                    arb(0),
                )
                for order in range(self.degree + 1)
            ]
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "ArbTaylorJet":
        if self.coefficients[0].contains(0):
            raise ZeroDivisionError("Arb Taylor constant ball contains zero")
        output = [1 / self.coefficients[0]]
        for order in range(1, self.degree + 1):
            tail = sum(
                (
                    self.coefficients[index] * output[order - index]
                    for index in range(1, order + 1)
                ),
                arb(0),
            )
            output.append(-tail / self.coefficients[0])
        return ArbTaylorJet(output)

    def __truediv__(self, other: object) -> "ArbTaylorJet":
        return self * as_arb_jet(other, self.degree).reciprocal()

    def __rtruediv__(self, other: object) -> "ArbTaylorJet":
        return as_arb_jet(other, self.degree) * self.reciprocal()

    def __pow__(self, exponent: int) -> "ArbTaylorJet":
        if exponent < 0:
            return self.reciprocal() ** (-exponent)
        output = as_arb_jet(1, self.degree)
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                output *= base
            remaining //= 2
            if remaining:
                base *= base
        return output


def as_arb_jet(value: object, degree: int) -> ArbTaylorJet:
    if isinstance(value, ArbTaylorJet):
        if value.degree != degree:
            raise ValueError("Taylor degrees do not match")
        return value
    if isinstance(value, Fraction):
        coefficient = arb(value.numerator) / value.denominator
    else:
        coefficient = arb(value)
    return ArbTaylorJet([coefficient] + [arb(0)] * degree)


def fraction_interval(value: Fraction) -> Interval:
    """Return a binary64 interval containing an exact rational."""

    approximation = float(value)
    represented = Fraction.from_float(approximation)
    lower = approximation if represented <= value else round_down(approximation)
    upper = approximation if represented >= value else round_up(approximation)
    return Interval(lower, upper)


def arb_theta_data(
    midpoint: Fraction, degree: int, *, terms: int = 12
) -> tuple[ArbTaylorJet, ...]:
    """Return rigorous high-precision center jets for the five nome factors."""

    ctx.dps = 80
    midpoint_ball = arb(midpoint.numerator) / midpoint.denominator
    return _arb_theta_data_from_ball(midpoint_ball, degree, terms=terms)


def arb_theta_data_box(
    lower: Fraction,
    upper: Fraction,
    degree: int,
    *,
    terms: int = 12,
) -> tuple[ArbTaylorJet, ...]:
    """Return Arb Taylor-coefficient enclosures over a rational nome box."""

    ctx.dps = 80
    lower_ball = arb(lower.numerator) / lower.denominator
    upper_ball = arb(upper.numerator) / upper.denominator
    return _arb_theta_data_from_ball(lower_ball.union(upper_ball), degree, terms=terms)


def _arb_theta_data_from_ball(
    c_ball: arb, degree: int, *, terms: int
) -> tuple[ArbTaylorJet, ...]:
    c = ArbTaylorJet([c_ball, arb(1)] + [arb(0)] * (degree - 1))
    theta = as_arb_jet(1, degree)
    r_series = as_arb_jet(1, degree)
    for index in range(1, terms + 1):
        theta += 2 * c ** (2 * index * index)
        r_series += c ** (2 * index * (index + 1))

    first_omitted = terms + 1
    theta_coefficients = list(theta.coefficients)
    r_coefficients = list(r_series.coefficients)
    for derivative in range(degree + 1):
        theta_coefficients[derivative] += arb_positive_series_tail(
            2,
            2 * first_omitted * first_omitted,
            derivative,
            c_ball,
            first_omitted,
        )
        r_coefficients[derivative] += arb_positive_series_tail(
            1,
            2 * first_omitted * (first_omitted + 1),
            derivative,
            c_ball,
            first_omitted,
        )
    theta = ArbTaylorJet(theta_coefficients)
    r_series = ArbTaylorJet(r_coefficients)
    s = theta**-2
    g = 2 * r_series / theta
    k = c * g**2
    gamma_lower = s * (1 + k**2 - s**2) / 6
    gamma_width = 1 - s - gamma_lower
    return c, k, s, gamma_lower, gamma_width


def arb_positive_series_tail(
    coefficient: int,
    exponent: int,
    derivative: int,
    c_ball: arb,
    first_index: int,
) -> arb:
    if exponent < derivative:
        raise ValueError("retain enough theta terms to cover the Taylor degree")
    if c_ball.is_zero():
        return arb(0)
    upper_magnitude = c_ball.abs_upper()
    first = (
        coefficient
        * comb(exponent, derivative)
        * upper_magnitude ** (exponent - derivative)
    )
    ratio = 4**derivative * upper_magnitude ** (4 * first_index + 2)
    if not ratio < 1:
        raise ValueError(
            "theta-tail ratio bound is not contractive; retain more series terms"
        )
    upper = (first / (1 - ratio)).upper()
    return arb(upper / 2, upper / 2)


def factor_jets(
    box: tuple[Fraction, Fraction],
    keys: set[tuple[int, ...]],
    degree: int,
) -> FactorJets:
    midpoint = (box[0] + box[1]) / 2
    center = arb_theta_data(midpoint, degree)
    whole = arb_theta_data_box(box[0], box[1], degree)
    return factor_jets_from_data(center, whole, keys, degree)


def factor_jets_from_data(
    center: tuple[ArbTaylorJet, ...],
    whole: tuple[ArbTaylorJet, ...],
    keys: set[tuple[int, ...]],
    degree: int,
) -> FactorJets:
    """Build products of common scalar jets for a sparse factor-key set."""

    output: FactorJets = {}
    for key in keys:
        center_factor = as_arb_jet(1, degree)
        whole_factor = as_arb_jet(1, degree)
        for exponent, center_value, whole_value in zip(key, center, whole):
            center_factor *= center_value**exponent
            whole_factor *= whole_value**exponent
        output[key] = center_factor, whole_factor
    return output


def final_chart_tables(core: CoreRecords, sign: int) -> list[ChartTable]:
    determinant_tables: list[ChartTable] = []
    for chart in (0, 1):
        table = main_chart(core.determinant, chart, 4)
        signed = {
            monomial: coefficient * sign ** monomial[4]
            for monomial, coefficient in table.items()
        }
        determinant_tables.append(ChartTable(f"det-{chart}-sign-{sign:+d}", signed))

    secondary = signed_secondary_records(core.minor3, sign)
    minor_tables = [
        ChartTable(f"minor-secondary-q-sign-{sign:+d}", secondary_chart(secondary, 1))
    ]
    tertiary = tertiary_records(secondary)
    for chart in (0, 1):
        table = {
            monomial[1:]: coefficient
            for monomial, coefficient in tertiary.items()
            if monomial[0] == chart
        }
        minor_tables.append(ChartTable(f"minor-tertiary-{chart}-sign-{sign:+d}", table))
    return determinant_tables + minor_tables


def interval_tensor_for_chart(
    table: ChartTable,
    box: tuple[Fraction, Fraction],
    degree: int,
    *,
    factor_provider: FactorJetProvider = factor_jets,
) -> IntervalTensor:
    keys = {monomial[5:] for monomial in table.records}
    factors = factor_provider(box, keys, degree)
    dimensions = tuple(
        max(monomial[axis] for monomial in table.records) + 1 for axis in range(5)
    )
    center_maps: list[dict[tuple[int, ...], arb]] = [dict() for _ in range(degree)]
    remainder_map: dict[tuple[int, ...], arb] = {}
    for monomial, coefficient in table.records.items():
        index = monomial[:5]
        center_factor, whole_factor = factors[monomial[5:]]
        for order in range(degree):
            value = coefficient * center_factor.coefficients[order]
            center_maps[order][index] = center_maps[order].get(index, arb(0)) + value
        value = coefficient * whole_factor.coefficients[degree]
        remainder_map[index] = remainder_map.get(index, arb(0)) + value

    center = [arb_map_to_tensor(values, dimensions) for values in center_maps]
    remainder = arb_map_to_tensor(remainder_map, dimensions)
    center_bernstein = [power_to_bernstein_tensor(part) for part in center]
    remainder_bernstein = power_to_bernstein_tensor(remainder)
    radius = (box[1] - box[0]) / 2
    c_power = IntervalTensor(
        np.zeros((degree,) + dimensions),
        np.zeros((degree,) + dimensions),
    )
    for order, part in enumerate(center_bernstein):
        for exponent in range(order + 1):
            scalar = (
                radius**order
                * comb(order, exponent)
                * 2**exponent
                * (-1) ** (order - exponent)
            )
            contribution = multiply_scalar(part, fraction_interval(scalar))
            _add_tensor_slice(c_power, exponent, contribution)
    result = power_to_bernstein_axis(c_power, 0)

    remainder_scale = fraction_interval(radius**degree)
    if degree % 2:
        displacement = Interval(-remainder_scale.upper, remainder_scale.upper)
    else:
        displacement = Interval(0, remainder_scale.upper)
    error = multiply_scalar(remainder_bernstein, displacement)
    for index in range(degree):
        _add_tensor_slice(result, index, error)
    return result


def arb_map_to_tensor(
    values: dict[tuple[int, ...], arb], dimensions: tuple[int, ...]
) -> IntervalTensor:
    lower = np.zeros(dimensions)
    upper = np.zeros(dimensions)
    for index, value in values.items():
        if value.is_zero():
            continue
        lower[index] = round_down(float(value.lower()))
        upper[index] = round_up(float(value.upper()))
    return IntervalTensor(lower, upper)


def _add_tensor_slice(
    target: IntervalTensor,
    index: int,
    contribution: IntervalTensor,
) -> None:
    target.lower[index] = directed_add(target.lower[index], contribution.lower, False)
    target.upper[index] = directed_add(target.upper[index], contribution.upper, True)


def directed_add(left: FloatArray, right: FloatArray, upward: bool) -> FloatArray:
    value = left + right
    rounded = np.nextafter(value, math.inf if upward else -math.inf)
    rounded = np.where(right == 0, left, rounded)
    rounded = np.where(left == 0, right, rounded)
    return rounded


def multiply_scalar(tensor: IntervalTensor, scalar: Interval) -> IntervalTensor:
    products = np.stack(
        (
            tensor.lower * scalar.lower,
            tensor.lower * scalar.upper,
            tensor.upper * scalar.lower,
            tensor.upper * scalar.upper,
        )
    )
    lower = np.nextafter(np.min(products, axis=0), -math.inf)
    upper = np.nextafter(np.max(products, axis=0), math.inf)
    exact_zero = (tensor.lower == 0) & (tensor.upper == 0)
    lower[exact_zero] = 0
    upper[exact_zero] = 0
    return IntervalTensor(lower, upper)


def multiply_positive_scalar(
    tensor: IntervalTensor, scalar: Interval
) -> IntervalTensor:
    """Multiply by a nonnegative scalar interval without four-product storage."""

    if scalar.lower < 0:
        raise ValueError("scalar interval is not nonnegative")
    lower_scale = np.where(tensor.lower < 0, scalar.upper, scalar.lower)
    upper_scale = np.where(tensor.upper > 0, scalar.upper, scalar.lower)
    lower = np.nextafter(tensor.lower * lower_scale, -math.inf)
    upper = np.nextafter(tensor.upper * upper_scale, math.inf)
    lower[tensor.lower == 0] = 0
    upper[tensor.upper == 0] = 0
    return IntervalTensor(lower, upper)


def power_to_bernstein_tensor(tensor: IntervalTensor) -> IntervalTensor:
    output = tensor
    for axis in range(tensor.lower.ndim):
        output = power_to_bernstein_axis(output, axis)
    return output


def power_to_bernstein_axis(tensor: IntervalTensor, axis: int) -> IntervalTensor:
    lower = np.moveaxis(tensor.lower, axis, 0)
    upper = np.moveaxis(tensor.upper, axis, 0)
    degree = lower.shape[0] - 1
    output_lower = np.zeros_like(lower)
    output_upper = np.zeros_like(upper)
    for bernstein_index in range(degree + 1):
        for power_index in range(bernstein_index + 1):
            coefficient = Fraction(
                comb(bernstein_index, power_index), comb(degree, power_index)
            )
            term = multiply_positive_scalar(
                IntervalTensor(lower[power_index], upper[power_index]),
                fraction_interval(coefficient),
            )
            output_lower[bernstein_index] = directed_add(
                output_lower[bernstein_index], term.lower, False
            )
            output_upper[bernstein_index] = directed_add(
                output_upper[bernstein_index], term.upper, True
            )
    return IntervalTensor(
        np.moveaxis(output_lower, 0, axis), np.moveaxis(output_upper, 0, axis)
    )


def split_axis(
    tensor: IntervalTensor, axis: int
) -> tuple[IntervalTensor, IntervalTensor]:
    lower = np.moveaxis(tensor.lower, axis, 0)
    upper = np.moveaxis(tensor.upper, axis, 0)
    degree = lower.shape[0] - 1
    left_lower = np.empty_like(lower)
    left_upper = np.empty_like(upper)
    right_lower = np.empty_like(lower)
    right_upper = np.empty_like(upper)
    left_lower[0] = lower[0]
    left_upper[0] = upper[0]
    right_lower[degree] = lower[degree]
    right_upper[degree] = upper[degree]

    work_lower = lower
    work_upper = upper
    for level in range(1, degree + 1):
        work_lower = directed_half_sum(work_lower[:-1], work_lower[1:], False)
        work_upper = directed_half_sum(work_upper[:-1], work_upper[1:], True)
        left_lower[level] = work_lower[0]
        left_upper[level] = work_upper[0]
        right_index = degree - level
        right_lower[right_index] = work_lower[-1]
        right_upper[right_index] = work_upper[-1]

    left = IntervalTensor(
        np.moveaxis(left_lower, 0, axis), np.moveaxis(left_upper, 0, axis)
    )
    right = IntervalTensor(
        np.moveaxis(right_lower, 0, axis), np.moveaxis(right_upper, 0, axis)
    )
    return left, right


def directed_half_sum(left: FloatArray, right: FloatArray, upward: bool) -> FloatArray:
    value = (left + right) * 0.5
    rounded = np.nextafter(value, math.inf if upward else -math.inf)
    exact_zero = (left == 0) & (right == 0)
    rounded[exact_zero] = 0
    return rounded


def split_scores(
    tensor: IntervalTensor,
    depths: tuple[int, ...],
    minimum_index: tuple[int, ...],
) -> list[float]:
    """Rank subdivision axes without large temporaries on oversized tensors."""

    if tensor.lower.size <= FULL_VARIATION_SCORE_MAX_SIZE:
        center = (tensor.lower + tensor.upper) * 0.5
        return [
            (
                float(np.max(np.abs(np.diff(center, axis=axis)))) / 2**depth
                if center.shape[axis] > 1
                else -1.0
            )
            for axis, depth in enumerate(depths)
        ]

    scores: list[float] = []
    for axis, (degree, depth) in enumerate(
        zip((size - 1 for size in tensor.lower.shape), depths)
    ):
        if degree == 0:
            scores.append(-1.0)
            continue
        left_index = list(minimum_index)
        right_index = list(minimum_index)
        left_index[axis] = max(0, minimum_index[axis] - 1)
        right_index[axis] = min(degree, minimum_index[axis] + 1)
        left = (tensor.lower[tuple(left_index)] + tensor.upper[tuple(left_index)]) * 0.5
        right = (
            tensor.lower[tuple(right_index)] + tensor.upper[tuple(right_index)]
        ) * 0.5
        scores.append(float(abs(right - left)) / 2**depth)
    return scores


def certify(
    root: IntervalTensor,
    *,
    max_depth: int,
    max_leaves: int,
    prevalidated_boxes: tuple[Box, ...] = (),
    prevalidated_regions: tuple[Callable[[Box], bool], ...] = (),
) -> CertificateResult:
    """Certify a tensor, optionally reusing independently certified regions."""

    unit_box = tuple((0.0, 1.0) for _ in range(root.lower.ndim))
    if any(len(box) != root.lower.ndim for box in prevalidated_boxes):
        raise ValueError("prevalidated box dimension does not match the tensor")
    stack = [(root, (0,) * root.lower.ndim, unit_box)]
    leaves = 0
    deepest = 0
    while stack:
        tensor, depths, box = stack.pop()
        inside_prevalidated_box = any(
            all(
                lower >= accepted_lower and upper <= accepted_upper
                for (lower, upper), (accepted_lower, accepted_upper) in zip(
                    box, accepted
                )
            )
            for accepted in prevalidated_boxes
        )
        if inside_prevalidated_box or any(
            contains(box) for contains in prevalidated_regions
        ):
            leaves += 1
            continue
        minimum_flat_index = int(np.argmin(tensor.lower))
        minimum = float(tensor.lower.flat[minimum_flat_index])
        if minimum >= 0:
            leaves += 1
            continue
        if max(depths) >= max_depth or leaves + len(stack) >= max_leaves:
            failure_index = tuple(
                int(index)
                for index in np.unravel_index(minimum_flat_index, tensor.lower.shape)
            )
            return CertificateResult(
                False, leaves, deepest, minimum, box, failure_index
            )
        minimum_index = tuple(
            int(index)
            for index in np.unravel_index(minimum_flat_index, tensor.lower.shape)
        )
        scores = split_scores(tensor, depths, minimum_index)
        axis = int(np.argmax(scores))
        left, right = split_axis(tensor, axis)
        new_depths = list(depths)
        new_depths[axis] += 1
        depth_tuple = tuple(new_depths)
        deepest = max(deepest, max(depth_tuple))
        lower, upper = box[axis]
        midpoint = (lower + upper) * 0.5
        left_box = list(box)
        right_box = list(box)
        left_box[axis] = (lower, midpoint)
        right_box[axis] = (midpoint, upper)
        stack.extend(
            (
                (left, depth_tuple, tuple(left_box)),
                (right, depth_tuple, tuple(right_box)),
            )
        )
    return CertificateResult(True, leaves, deepest, 0.0)


def parse_fraction(value: str) -> Fraction:
    try:
        return Fraction(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lower", type=parse_fraction)
    parser.add_argument("upper", type=parse_fraction)
    parser.add_argument("--degree", type=int, default=10)
    parser.add_argument("--sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--chart", help="run only labels containing this text")
    parser.add_argument("--max-depth", type=int, default=40)
    parser.add_argument("--max-leaves", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 < args.lower < args.upper < 1:
        raise SystemExit("require 0 < lower < upper < 1")
    if args.degree < 1:
        raise SystemExit("require degree >= 1")
    selected = final_chart_tables(load_records(), args.sign)
    if args.chart:
        selected = [table for table in selected if args.chart in table.label]
    if not selected:
        raise SystemExit("no chart matched")
    for table in selected:
        started = time.monotonic()
        tensor = interval_tensor_for_chart(table, (args.lower, args.upper), args.degree)
        built = time.monotonic()
        result = certify(tensor, max_depth=args.max_depth, max_leaves=args.max_leaves)
        print(
            f"{table.label}: shape={tensor.lower.shape}, build={built - started:.2f}s, "
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
