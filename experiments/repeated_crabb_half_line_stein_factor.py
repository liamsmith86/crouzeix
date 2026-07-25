#!/usr/bin/env python3
"""Audit the exact half-line Stein factor in coisometric word algebra.

For the backward-shift half-line, let ``s a = 1`` and
``E = 1 - a s``.  The right-balanced ellipse pencil is

    Xi = s + c a(1 + E).

The half-line boundary metric has diagonal orbit weights
``1/(1+c^(2n))``.  This script checks, in exact rational truncated
series, that its elliptic Stein slack is the rank-one Gram generated
by the Jacobi/theta Fourier column recorded in the repeated-Crabb
axis calculation.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import direct_map_coefficients


Polynomial = dict[str, Fraction]
Series = list[Polynomial]


@dataclass(frozen=True)
class HalfLineSteinFactorRecord:
    """One exact coefficient audit."""

    degree: int
    slack_word_count: int
    factor_word_count: int
    residual_word_count: int
    all_checks_passed: bool


def reduce_word(word: str) -> str:
    """Reduce a word by the coisometry relation ``s a = 1``."""

    while "sa" in word:
        word = word.replace("sa", "")
    return word


def clean(polynomial: Polynomial) -> Polynomial:
    """Remove zero coefficients."""

    return {
        word: coefficient
        for word, coefficient in polynomial.items()
        if coefficient
    }


def add(*polynomials: Polynomial) -> Polynomial:
    """Add sparse word polynomials."""

    result: defaultdict[str, Fraction] = defaultdict(Fraction)
    for polynomial in polynomials:
        for word, coefficient in polynomial.items():
            result[word] += coefficient
    return clean(dict(result))


def scale(coefficient: Fraction | int, polynomial: Polynomial) -> Polynomial:
    """Scale one sparse polynomial."""

    coefficient = Fraction(coefficient)
    return clean({
        word: coefficient * value
        for word, value in polynomial.items()
    })


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    """Multiply and coisometrically reduce two polynomials."""

    result: defaultdict[str, Fraction] = defaultdict(Fraction)
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            result[reduce_word(left_word + right_word)] += (
                left_coefficient * right_coefficient
            )
    return clean(dict(result))


def adjoint(polynomial: Polynomial) -> Polynomial:
    """Apply the formal adjoint."""

    translation = str.maketrans({"s": "a", "a": "s"})
    return clean({
        reduce_word(word.translate(translation)[::-1]): coefficient
        for word, coefficient in polynomial.items()
    })


IDENTITY = {"": Fraction(1)}
S = {"s": Fraction(1)}
STAR = {"a": Fraction(1)}
E = add(IDENTITY, scale(-1, multiply(STAR, S)))


def zero_series(maximum_degree: int) -> Series:
    """Return a zero series through one degree."""

    return [{} for _ in range(maximum_degree + 1)]


def series_add(*series: Series) -> Series:
    """Add equal-length word series."""

    return [
        add(*(item[degree] for item in series))
        for degree in range(len(series[0]))
    ]


def series_multiply(left: Series, right: Series) -> Series:
    """Multiply equal-length word series."""

    maximum_degree = len(left) - 1
    result = zero_series(maximum_degree)
    for left_degree, left_polynomial in enumerate(left):
        for right_degree, right_polynomial in enumerate(right):
            degree = left_degree + right_degree
            if degree <= maximum_degree:
                result[degree] = add(
                    result[degree],
                    multiply(left_polynomial, right_polynomial),
                )
    return result


def series_adjoint(series: Series) -> Series:
    """Apply the formal adjoint coefficientwise."""

    return [adjoint(polynomial) for polynomial in series]


def series_power(series: Series, exponent: int) -> Series:
    """Raise a word series to a nonnegative integer power."""

    maximum_degree = len(series) - 1
    result = [IDENTITY, *({} for _ in range(maximum_degree))]
    base = series
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = series_multiply(result, base)
        base = series_multiply(base, base)
        remaining //= 2
    return result


def scalar_times_series(
    scalars: list[Fraction],
    series: Series,
) -> Series:
    """Convolve a scalar series with a word series."""

    maximum_degree = len(series) - 1
    result = zero_series(maximum_degree)
    for scalar_degree, scalar in enumerate(scalars):
        if not scalar:
            continue
        for word_degree, polynomial in enumerate(series):
            degree = scalar_degree + word_degree
            if degree <= maximum_degree:
                result[degree] = add(
                    result[degree],
                    scale(scalar, polynomial),
                )
    return result


def operator_series(maximum_degree: int) -> Series:
    """Return the half-line elliptic direct-map series."""

    right_correction = multiply(STAR, add(IDENTITY, E))
    pencil = [
        S,
        right_correction,
        *({} for _ in range(maximum_degree - 1)),
    ]
    result = zero_series(maximum_degree)
    for index, coefficient in enumerate(
        direct_map_coefficients(maximum_degree, maximum_degree + 1)
    ):
        scalars = [
            coefficient.coefficient(degree)
            for degree in range(maximum_degree + 1)
        ]
        result = series_add(
            result,
            scalar_times_series(
                scalars,
                series_power(pencil, 2 * index + 1),
            ),
        )
    return result


def metric_series(maximum_degree: int) -> Series:
    """Return the half-line right-orbit metric."""

    result = [IDENTITY, *({} for _ in range(maximum_degree))]
    for orbit in range(1, maximum_degree // 2 + 1):
        projection = multiply(
            multiply({"a" * orbit: Fraction(1)}, E),
            {"s" * orbit: Fraction(1)},
        )
        multiple = 1
        while 2 * orbit * multiple <= maximum_degree:
            degree = 2 * orbit * multiple
            result[degree] = add(
                result[degree],
                scale((-1) ** multiple, projection),
            )
            multiple += 1
    return result


def inverse_theta_three(maximum_degree: int) -> list[Fraction]:
    """Return ``1/theta_3(c^2)`` through one degree."""

    theta = [Fraction(0) for _ in range(maximum_degree + 1)]
    theta[0] = Fraction(1)
    index = 1
    while 2 * index * index <= maximum_degree:
        theta[2 * index * index] = Fraction(2)
        index += 1
    inverse = [Fraction(0) for _ in range(maximum_degree + 1)]
    inverse[0] = Fraction(1)
    for degree in range(1, maximum_degree + 1):
        inverse[degree] = -sum(
            theta[positive] * inverse[degree - positive]
            for positive in range(1, degree + 1)
        )
    return inverse


def defect_column_series(
    maximum_degree: int,
) -> list[dict[int, Fraction]]:
    """Return the theta/Fourier defect column by orbit exponent."""

    inverse_theta = inverse_theta_three(maximum_degree)
    result: list[defaultdict[int, Fraction]] = [
        defaultdict(Fraction)
        for _ in range(maximum_degree + 1)
    ]
    for degree, coefficient in enumerate(inverse_theta):
        result[degree][0] += coefficient

    for frequency in range(1, maximum_degree + 1):
        denominator_power = 0
        while (
            frequency + 4 * frequency * denominator_power
            <= maximum_degree
        ):
            shift = frequency + 4 * frequency * denominator_power
            denominator_sign = (-1) ** denominator_power
            for theta_degree, theta_coefficient in enumerate(
                inverse_theta
            ):
                degree = shift + theta_degree
                if degree > maximum_degree:
                    break
                result[degree][2 * frequency] += (
                    2
                    * (-1) ** frequency
                    * denominator_sign
                    * theta_coefficient
                )
            denominator_power += 1
    return [dict(coefficient) for coefficient in result]


def factor_series(maximum_degree: int) -> Series:
    """Return the rank-one Gram series of the defect column."""

    column = defect_column_series(maximum_degree)
    result = zero_series(maximum_degree)
    for left_degree, left_orbits in enumerate(column):
        for right_degree, right_orbits in enumerate(column):
            degree = left_degree + right_degree
            if degree > maximum_degree:
                continue
            for left_orbit, left_coefficient in left_orbits.items():
                for right_orbit, right_coefficient in right_orbits.items():
                    projection = multiply(
                        multiply(
                            {"a" * left_orbit: Fraction(1)},
                            E,
                        ),
                        {"s" * right_orbit: Fraction(1)},
                    )
                    result[degree] = add(
                        result[degree],
                        scale(
                            left_coefficient * right_coefficient,
                            projection,
                        ),
                    )
    return result


def exact_records(maximum_degree: int) -> list[HalfLineSteinFactorRecord]:
    """Return exact coefficient records through one degree."""

    operator = operator_series(maximum_degree)
    metric = metric_series(maximum_degree)
    pulled = series_multiply(
        series_multiply(series_adjoint(operator), metric),
        operator,
    )
    slack = series_add(
        metric,
        [scale(-1, coefficient) for coefficient in pulled],
    )
    factor = factor_series(maximum_degree)

    records: list[HalfLineSteinFactorRecord] = []
    for degree, (slack_coefficient, factor_coefficient) in enumerate(
        zip(slack, factor, strict=True)
    ):
        residual = add(
            slack_coefficient,
            scale(-1, factor_coefficient),
        )
        verified = not residual
        if not verified:
            raise RuntimeError(
                "the half-line Stein-factor audit failed: "
                f"degree={degree}, residual={residual}"
            )
        records.append(
            HalfLineSteinFactorRecord(
                degree=degree,
                slack_word_count=len(slack_coefficient),
                factor_word_count=len(factor_coefficient),
                residual_word_count=len(residual),
                all_checks_passed=verified,
            )
        )
    return records


def write_records(
    records: list[HalfLineSteinFactorRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-degree", type=int, default=8)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_half_line_stein_factor_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact half-line audit."""

    args = parse_args()
    records = exact_records(args.maximum_degree)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
