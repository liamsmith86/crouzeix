#!/usr/bin/env python3
"""Directed-interval certificate for the compact low-nome odd-block range.

This checker proves the two endpoint residuals in L32 nonnegative on their
discriminant branch for

    1/12 <= c <= 12599/20000,  0 <= p <= 1.

The lower part uses the exact ridge vertex p_star as a coordinate split;
the upper part uses (c,p) directly.  Every accepted box either misses the
branch ``A>0, 2A+B>0, 2A-B>0`` or has a nonnegative interval lower bound
for ``r*F1*F2-g**2*D**2``.  Arithmetic is binary64 with one-ulp outward
rounding after every elementary operation.  Theta tails and their first two
derivatives are enclosed by explicit geometric majorants.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import time
from typing import Literal


NEGATIVE_INFINITY = -math.inf
POSITIVE_INFINITY = math.inf


def round_down(value: float) -> float:
    if math.isnan(value):
        raise ArithmeticError("NaN in directed interval operation")
    if math.isinf(value):
        return value
    return math.nextafter(value, NEGATIVE_INFINITY)


def round_up(value: float) -> float:
    if math.isnan(value):
        raise ArithmeticError("NaN in directed interval operation")
    if math.isinf(value):
        return value
    return math.nextafter(value, POSITIVE_INFINITY)


@dataclass(frozen=True, slots=True)
class Interval:
    lower: float
    upper: float

    def __add__(self, other: object) -> "Interval":
        rhs = as_interval(other)
        return Interval(
            round_down(self.lower + rhs.lower),
            round_up(self.upper + rhs.upper),
        )

    __radd__ = __add__

    def __neg__(self) -> "Interval":
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other: object) -> "Interval":
        return self + (-as_interval(other))

    def __rsub__(self, other: object) -> "Interval":
        return as_interval(other) - self

    def __mul__(self, other: object) -> "Interval":
        rhs = as_interval(other)
        products = (
            self.lower * rhs.lower,
            self.lower * rhs.upper,
            self.upper * rhs.lower,
            self.upper * rhs.upper,
        )
        return Interval(round_down(min(products)), round_up(max(products)))

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "Interval":
        rhs = as_interval(other)
        if rhs.lower <= 0 <= rhs.upper:
            return Interval(NEGATIVE_INFINITY, POSITIVE_INFINITY)
        reciprocal = Interval(round_down(1 / rhs.upper), round_up(1 / rhs.lower))
        return self * reciprocal

    def __rtruediv__(self, other: object) -> "Interval":
        return as_interval(other) / self

    def __pow__(self, exponent: int) -> "Interval":
        if not isinstance(exponent, int) or exponent < 0:
            if isinstance(exponent, int):
                return Interval(1.0, 1.0) / self ** (-exponent)
            raise TypeError(exponent)
        result = Interval(1.0, 1.0)
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = result * base
            remaining //= 2
            if remaining:
                base = base * base
        return result


def as_interval(value: object) -> Interval:
    if isinstance(value, Interval):
        return value
    scalar = float(value)
    return Interval(scalar, scalar)


ZERO = Interval(0.0, 0.0)
ONE = Interval(1.0, 1.0)


@dataclass(frozen=True, slots=True)
class Jet:
    value: Interval
    c: Interval = ZERO
    p: Interval = ZERO
    cc: Interval = ZERO
    cp: Interval = ZERO
    pp: Interval = ZERO

    def __add__(self, other: object) -> "Jet":
        rhs = as_jet(other)
        return Jet(
            self.value + rhs.value,
            self.c + rhs.c,
            self.p + rhs.p,
            self.cc + rhs.cc,
            self.cp + rhs.cp,
            self.pp + rhs.pp,
        )

    __radd__ = __add__

    def __neg__(self) -> "Jet":
        return Jet(-self.value, -self.c, -self.p, -self.cc, -self.cp, -self.pp)

    def __sub__(self, other: object) -> "Jet":
        return self + (-as_jet(other))

    def __rsub__(self, other: object) -> "Jet":
        return as_jet(other) - self

    def __mul__(self, other: object) -> "Jet":
        rhs = as_jet(other)
        return Jet(
            self.value * rhs.value,
            self.c * rhs.value + self.value * rhs.c,
            self.p * rhs.value + self.value * rhs.p,
            self.cc * rhs.value + 2 * self.c * rhs.c + self.value * rhs.cc,
            self.cp * rhs.value
            + self.c * rhs.p
            + self.p * rhs.c
            + self.value * rhs.cp,
            self.pp * rhs.value + 2 * self.p * rhs.p + self.value * rhs.pp,
        )

    __rmul__ = __mul__

    def compose(
        self,
        value: Interval,
        first: Interval,
        second: Interval,
    ) -> "Jet":
        return Jet(
            value,
            first * self.c,
            first * self.p,
            second * self.c * self.c + first * self.cc,
            second * self.c * self.p + first * self.cp,
            second * self.p * self.p + first * self.pp,
        )

    def __pow__(self, exponent: int) -> "Jet":
        value = self.value**exponent
        first = exponent * self.value ** (exponent - 1)
        second = exponent * (exponent - 1) * self.value ** (exponent - 2)
        return self.compose(value, first, second)

    def __truediv__(self, other: object) -> "Jet":
        return self * as_jet(other) ** -1

    def __rtruediv__(self, other: object) -> "Jet":
        return as_jet(other) * self**-1


def as_jet(value: object) -> Jet:
    if isinstance(value, Jet):
        return value
    return Jet(as_interval(value))


def theta_data(c: Jet, terms: int = 6) -> tuple[Jet, Jet, Jet]:
    """Enclose s=T^-2, g=2R/T, and k=c*g^2 through second derivatives."""

    theta = as_jet(1)
    r_series = as_jet(1)
    for n in range(1, terms + 1):
        theta += 2 * c ** (2 * n * n)
        r_series += c ** (2 * n * (n + 1))

    n = terms + 1
    c_value = c.value
    theta_value_tail = 2 * c_value ** (2 * n * n) / (1 - c_value ** (4 * n + 2))
    theta_first_tail = (
        4
        * n
        * n
        * c_value ** (2 * n * n - 1)
        / (1 - 4 * c_value ** (4 * n + 2))
    )
    theta_second_tail = (
        4
        * n
        * n
        * (2 * n * n - 1)
        * c_value ** (2 * n * n - 2)
        / (1 - 16 * c_value ** (4 * n + 2))
    )
    r_value_tail = c_value ** (2 * n * (n + 1)) / (
        1 - c_value ** (4 * n + 4)
    )
    r_first_tail = (
        2
        * n
        * (n + 1)
        * c_value ** (2 * n * (n + 1) - 1)
        / (1 - 4 * c_value ** (4 * n + 4))
    )
    r_second_tail = (
        2
        * n
        * (n + 1)
        * (2 * n * (n + 1) - 1)
        * c_value ** (2 * n * (n + 1) - 2)
        / (1 - 16 * c_value ** (4 * n + 4))
    )
    theta += Jet(
        Interval(0.0, theta_value_tail.upper),
        Interval(0.0, theta_first_tail.upper),
        ZERO,
        Interval(0.0, theta_second_tail.upper),
    )
    r_series += Jet(
        Interval(0.0, r_value_tail.upper),
        Interval(0.0, r_first_tail.upper),
        ZERO,
        Interval(0.0, r_second_tail.upper),
    )
    s = theta**-2
    g = 2 * r_series / theta
    k = c * g * g
    return s, g, k


def endpoint_forms(
    c: Jet,
    p: Jet,
    endpoint: Literal["L", "U"],
    data: tuple[Jet, Jet, Jet] | None = None,
) -> tuple[Jet, ...]:
    """Return the three branch forms and L32's cancellation-free residual."""

    s, g, k = theta_data(c) if data is None else data
    d = k * (1 - p * p) / (1 - k * k * p * p)
    a3 = s * (1 + k * k - s * s) / 6
    r = s * p + (a3 if endpoint == "L" else 1 - s) * p**3
    coefficient_a = g * p * (2 * d * (1 + c * c * r) - c * g * (1 + r)) / (
        c * (1 + r)
    )
    coefficient_b = (
        2 * c * c * g * p * r
        - 2 * c * c * g
        - c * d * g * g * p * (1 + r)
        + 4 * c * d * (1 + r)
        + 2 * g * p
        - 2 * g * r
    ) / (c * (1 + r))
    lam = d * (k * p + 4 * c) / (2 * g)
    f1 = -4 * c * c * d + 2 * c * g * (p + 1) - d * g * g * p
    f2 = c * c * d * g * g * p - 2 * c * g * (p + 1) + 4 * d
    gap = (
        (1 - s) * p * (1 - p * p)
        if endpoint == "U"
        else p * (1 - s - a3 * p * p)
    )
    small_d = gap + c * c * (1 - p * r) - lam * (1 - r)
    residual = r * f1 * f2 - g * g * small_d * small_d
    return (
        coefficient_a,
        2 * coefficient_a + coefficient_b,
        2 * coefficient_a - coefficient_b,
        residual,
    )


def p_star(c: Jet, g: Jet) -> Jet:
    return (g * g - 4 * c * c) / (2 * g * g * (1 - 2 * c * c * g))


def coordinate_forms(
    c_lower: float,
    c_upper: float,
    p_lower: float,
    p_upper: float,
    endpoint: Literal["L", "U"],
    coordinate: Literal["direct", "lower", "upper"],
) -> tuple[Jet, ...]:
    c = Jet(Interval(c_lower, c_upper), ONE)
    parameter = Jet(Interval(p_lower, p_upper), ZERO, ONE)
    data = theta_data(c)
    if coordinate == "direct":
        p = parameter
    else:
        _, g, _ = data
        center = p_star(c, g)
        p = center * parameter if coordinate == "lower" else parameter + (1 - parameter) * center
    return endpoint_forms(c, p, endpoint, data)


def taylor_enclosure(center: Jet, whole: Jet, c_width: float, p_width: float) -> Interval:
    delta_c = Interval(-c_width, c_width)
    delta_p = Interval(-p_width, p_width)
    delta_c_squared = Interval(0.0, round_up(c_width * c_width))
    delta_p_squared = Interval(0.0, round_up(p_width * p_width))
    cross = Interval(-round_up(c_width * p_width), round_up(c_width * p_width))
    return (
        center.value
        + center.c * delta_c
        + center.p * delta_p
        + (
            whole.cc * delta_c_squared
            + 2 * whole.cp * cross
            + whole.pp * delta_p_squared
        )
        / 2
    )


@dataclass(frozen=True, slots=True)
class Box:
    c_lower: float
    c_upper: float
    p_lower: float
    p_upper: float
    endpoint: Literal["L", "U"]
    coordinate: Literal["direct", "lower", "upper"]
    depth: int = 0


def evaluate_box(box: Box) -> Literal["outside", "proved", "split"]:
    c_midpoint = (box.c_lower + box.c_upper) / 2
    p_midpoint = (box.p_lower + box.p_upper) / 2
    whole = coordinate_forms(
        box.c_lower,
        box.c_upper,
        box.p_lower,
        box.p_upper,
        box.endpoint,
        box.coordinate,
    )
    center = coordinate_forms(
        c_midpoint,
        c_midpoint,
        p_midpoint,
        p_midpoint,
        box.endpoint,
        box.coordinate,
    )
    values = tuple(
        taylor_enclosure(
            center_value,
            whole_value,
            (box.c_upper - box.c_lower) / 2,
            (box.p_upper - box.p_lower) / 2,
        )
        for center_value, whole_value in zip(center, whole, strict=True)
    )
    if any(value.upper <= 0 for value in values[:3]):
        return "outside"
    if values[3].lower >= 0:
        return "proved"
    return "split"


def outward_fraction(value: Fraction, upper: bool) -> float:
    converted = float(value)
    return round_up(converted) if upper else round_down(converted)


def initial_boxes(
    c_start: Fraction,
    c_stop: Fraction,
    c_step: Fraction,
    p_step: Fraction,
    coordinates: tuple[Literal["direct", "lower", "upper"], ...],
) -> list[Box]:
    boxes = []
    c_lower_fraction = c_start
    while c_lower_fraction < c_stop:
        c_upper_fraction = min(c_lower_fraction + c_step, c_stop)
        p_lower_fraction = Fraction(0)
        while p_lower_fraction < 1:
            p_upper_fraction = min(p_lower_fraction + p_step, Fraction(1))
            for endpoint in ("L", "U"):
                for coordinate in coordinates:
                    boxes.append(
                        Box(
                            outward_fraction(c_lower_fraction, False),
                            outward_fraction(c_upper_fraction, True),
                            outward_fraction(p_lower_fraction, False),
                            outward_fraction(p_upper_fraction, True),
                            endpoint,
                            coordinate,
                        )
                    )
            p_lower_fraction = p_upper_fraction
        c_lower_fraction = c_upper_fraction
    return boxes


def audit_p_star_range(
    c_start: Fraction,
    c_stop: Fraction,
    c_step: Fraction,
) -> None:
    """Prove that the two ridge-coordinate charts cover the full p interval."""

    lower_fraction = c_start
    while lower_fraction < c_stop:
        upper_fraction = min(lower_fraction + c_step, c_stop)
        c = Jet(
            Interval(
                outward_fraction(lower_fraction, False),
                outward_fraction(upper_fraction, True),
            ),
            ONE,
        )
        _, g, _ = theta_data(c)
        center = p_star(c, g).value
        if center.lower < 0 or center.upper > 1:
            raise AssertionError((lower_fraction, upper_fraction, center))
        lower_fraction = upper_fraction


def split_box(box: Box) -> tuple[Box, Box]:
    relative_c_width = (box.c_upper - box.c_lower) / box.c_lower
    p_width = box.p_upper - box.p_lower
    if relative_c_width > p_width:
        midpoint = math.sqrt(box.c_lower * box.c_upper)
        return (
            Box(
                box.c_lower,
                midpoint,
                box.p_lower,
                box.p_upper,
                box.endpoint,
                box.coordinate,
                box.depth + 1,
            ),
            Box(
                midpoint,
                box.c_upper,
                box.p_lower,
                box.p_upper,
                box.endpoint,
                box.coordinate,
                box.depth + 1,
            ),
        )
    midpoint = (box.p_lower + box.p_upper) / 2
    return (
        Box(
            box.c_lower,
            box.c_upper,
            box.p_lower,
            midpoint,
            box.endpoint,
            box.coordinate,
            box.depth + 1,
        ),
        Box(
            box.c_lower,
            box.c_upper,
            midpoint,
            box.p_upper,
            box.endpoint,
            box.coordinate,
            box.depth + 1,
        ),
    )


def certify(boxes: list[Box], label: str, max_splits: int = 150_000) -> tuple[int, ...]:
    counts = {"outside": 0, "proved": 0, "split": 0}
    splits = 0
    maximum_depth = 0
    started = time.monotonic()
    while boxes:
        box = boxes.pop()
        status = evaluate_box(box)
        counts[status] += 1
        if status != "split":
            continue
        maximum_depth = max(maximum_depth, box.depth)
        if box.depth >= 40 or splits >= max_splits:
            raise AssertionError((label, box, counts, len(boxes)))
        boxes.extend(split_box(box))
        splits += 1
        if splits % 10_000 == 0:
            print(
                f"{label}: {splits} splits, {len(boxes)} pending "
                f"({time.monotonic() - started:.1f}s)",
                flush=True,
            )
    print(
        f"{label}: certified with {splits} splits, depth {maximum_depth}, "
        f"counts {counts} ({time.monotonic() - started:.1f}s)",
        flush=True,
    )
    return counts["outside"], counts["proved"], splits, maximum_depth


def main() -> None:
    audit_p_star_range(Fraction(1, 12), Fraction(1, 2), Fraction(1, 200))
    compact_boxes = initial_boxes(
        Fraction(1, 12),
        Fraction(1, 2),
        Fraction(1, 200),
        Fraction(1, 200),
        ("lower", "upper"),
    )
    high_boxes = initial_boxes(
        Fraction(1, 2),
        Fraction(12599, 20000),
        Fraction(1, 500),
        Fraction(1, 200),
        ("direct",),
    )
    compact = certify(compact_boxes, "ridge-coordinate range")
    high = certify(high_boxes, "direct-coordinate range")
    print("compact low-nome endpoint certificate: proved")
    print(f"  ridge-coordinate summary: {compact}")
    print(f"  direct-coordinate summary: {high}")


if __name__ == "__main__":
    main()
