#!/usr/bin/env python3
"""Audit the delayed boundary-slack residual in exact word algebra.

Words use ``s`` for S and ``a`` for S*.  The calculation imposes the
partial-isometry relations, orthogonality of the two defect spaces, and the
complete-delay relations ``E S^j F=0`` and their adjoints.  Through every
requested finite grade, it checks

    [c^(2k)] K_S(c) = E_1 F_(k-1) + F_(k-1) E_1.

This is an exact finite-grade falsification audit, not the arbitrary-grade
induction needed to prove the displayed identity.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from crabb_palindromic_elliptic_hessian import direct_map_coefficients


Word = str
Polynomial = dict[Word, Fraction]
Series = list[Polynomial]


@dataclass(frozen=True)
class DelayedSlackAnticommutatorRecord:
    """One exact complete-delay word-algebra audit."""

    grade: int
    face_degree: int
    face_word_count: int
    unreduced_face_difference_word_count: int
    maximum_earlier_reduced_word_count: int
    reduced_face_difference_word_count: int
    all_checks_passed: bool


@lru_cache(maxsize=500_000)
def reduce_word(word: Word) -> Polynomial:
    """Reduce a word by the partial-isometry and orthogonal-defect rules."""

    for pattern, replacement in (("sas", "s"), ("asa", "a")):
        position = word.find(pattern)
        if position >= 0:
            reduced = (
                word[:position] + replacement + word[position + 3 :]
            )
            return reduce_word(reduced)

    # EF=FE=0 after E=1-as and F=1-sa.
    for pattern in ("assa", "saas"):
        position = word.find(pattern)
        if position < 0:
            continue
        prefix = word[:position]
        suffix = word[position + 4 :]
        return add(
            reduce_word(prefix + "as" + suffix),
            reduce_word(prefix + "sa" + suffix),
            scale(-1, reduce_word(prefix + suffix)),
        )
    return {word: Fraction(1)}


def clean(polynomial: Polynomial) -> Polynomial:
    """Remove zero coefficients."""

    return {
        word: coefficient
        for word, coefficient in polynomial.items()
        if coefficient
    }


def add(*polynomials: Polynomial) -> Polynomial:
    """Add sparse word polynomials."""

    result: Polynomial = {}
    for polynomial in polynomials:
        for word, coefficient in polynomial.items():
            result[word] = result.get(word, Fraction()) + coefficient
    return clean(result)


def scale(coefficient: Fraction | int, polynomial: Polynomial) -> Polynomial:
    """Scale a sparse word polynomial."""

    coefficient = Fraction(coefficient)
    return clean(
        {
            word: coefficient * value
            for word, value in polynomial.items()
        }
    )


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    """Multiply sparse word polynomials."""

    result: Polynomial = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            for word, reduction_coefficient in reduce_word(
                left_word + right_word
            ).items():
                result[word] = (
                    result.get(word, Fraction())
                    + left_coefficient
                    * right_coefficient
                    * reduction_coefficient
                )
    return clean(result)


def adjoint(polynomial: Polynomial) -> Polynomial:
    """Apply the formal adjoint."""

    translation = str.maketrans({"s": "a", "a": "s"})
    return {
        word.translate(translation)[::-1]: coefficient
        for word, coefficient in polynomial.items()
    }


IDENTITY = {"": Fraction(1)}
S = {"s": Fraction(1)}
STAR = {"a": Fraction(1)}
E = add(IDENTITY, scale(-1, multiply(STAR, S)))
F = add(IDENTITY, scale(-1, multiply(S, STAR)))
Q = multiply(STAR, S)
E1 = multiply(multiply(STAR, E), S)


def zero_series(degree: int) -> Series:
    """Return a zero series through one degree."""

    return [{} for _ in range(degree + 1)]


def series_add(*series: Series) -> Series:
    """Add equal-length polynomial series."""

    return [
        add(*(item[degree] for item in series))
        for degree in range(len(series[0]))
    ]


def series_scale(scalars: list[Fraction], series: Series) -> Series:
    """Convolve a scalar series with a polynomial series."""

    degree = len(series) - 1
    result = zero_series(degree)
    for left_degree, scalar in enumerate(scalars):
        if not scalar:
            continue
        for right_degree, polynomial in enumerate(series):
            if left_degree + right_degree <= degree:
                result[left_degree + right_degree] = add(
                    result[left_degree + right_degree],
                    scale(scalar, polynomial),
                )
    return result


def series_multiply(left: Series, right: Series) -> Series:
    """Multiply polynomial series."""

    degree = len(left) - 1
    result = zero_series(degree)
    for left_degree, left_polynomial in enumerate(left):
        for right_degree, right_polynomial in enumerate(right):
            if left_degree + right_degree <= degree:
                result[left_degree + right_degree] = add(
                    result[left_degree + right_degree],
                    multiply(left_polynomial, right_polynomial),
                )
    return result


def series_adjoint(series: Series) -> Series:
    """Apply the formal adjoint coefficientwise."""

    return [adjoint(polynomial) for polynomial in series]


def series_power(series: Series, exponent: int) -> Series:
    """Raise a polynomial series to a nonnegative power."""

    degree = len(series) - 1
    result = [IDENTITY] + [{} for _ in range(degree)]
    for _ in range(exponent):
        result = series_multiply(result, series)
    return result


def operator_series(degree: int) -> Series:
    """Return the balanced elliptic operator series."""

    j_operator = multiply(
        add(IDENTITY, F),
        multiply(STAR, add(IDENTITY, E)),
    )
    pencil = [S, j_operator] + [{} for _ in range(degree - 1)]
    result = zero_series(degree)
    for index, coefficient_series in enumerate(
        direct_map_coefficients(degree, degree + 1)
    ):
        scalars = [
            coefficient_series.coefficient(order)
            for order in range(degree + 1)
        ]
        result = series_add(
            result,
            series_scale(
                scalars,
                series_power(pencil, 2 * index + 1),
            ),
        )
    return result


def boundary_metric_series(degree: int) -> Series:
    """Return L219's boundary metric series."""

    result = [IDENTITY] + [{} for _ in range(degree)]
    for order in range(2, degree + 1, 2):
        half_order = order // 2
        coefficient = multiply(
            multiply({"s" * half_order: Fraction(1)}, F),
            {"a" * half_order: Fraction(1)},
        )
        for divisor in range(1, half_order + 1):
            if half_order % divisor:
                continue
            right_orbit = multiply(
                multiply({"a" * divisor: Fraction(1)}, E),
                {"s" * divisor: Fraction(1)},
            )
            coefficient = add(
                coefficient,
                scale((-1) ** (half_order // divisor), right_orbit),
            )
        result[order] = coefficient
    return result


def schur_residual_series(degree: int) -> Series:
    """Return the lifted right-defect Schur residual."""

    operator = operator_series(degree)
    metric = boundary_metric_series(degree)
    pulled = series_multiply(
        series_multiply(series_adjoint(operator), metric),
        operator,
    )
    slack = series_add(metric, [scale(-1, item) for item in pulled])

    pivot_inverse = zero_series(degree)
    pivot_inverse[0] = E
    for order in range(1, degree + 1):
        convolution: Polynomial = {}
        for positive_degree in range(1, order + 1):
            pivot_coefficient = multiply(
                multiply(E, slack[positive_degree]),
                E,
            )
            convolution = add(
                convolution,
                multiply(
                    pivot_coefficient,
                    pivot_inverse[order - positive_degree],
                ),
            )
        pivot_inverse[order] = scale(-1, convolution)

    square = series_multiply(
        series_multiply(slack, pivot_inverse),
        slack,
    )
    return [
        multiply(
            multiply(Q, add(slack_coefficient, scale(-1, square_coefficient))),
            Q,
        )
        for slack_coefficient, square_coefficient in zip(
            slack,
            square,
            strict=True,
        )
    ]


def target(grade: int) -> Polynomial:
    """Return E_1 F_(grade-1) + F_(grade-1) E_1."""

    delay = grade - 1
    left_orbit = multiply(
        multiply({"s" * delay: Fraction(1)}, F),
        {"a" * delay: Fraction(1)},
    )
    return add(
        multiply(E1, left_orbit),
        multiply(left_orbit, E1),
    )


@lru_cache(maxsize=100_000)
def reduce_delayed_word(word: Word, maximum_delay: int) -> Polynomial:
    """Reduce one base-normal word under E S^j F=0 and its adjoint."""

    base = reduce_word(word)
    if len(base) != 1 or next(iter(base.values())) != 1:
        result: Polynomial = {}
        for base_word, coefficient in base.items():
            result = add(
                result,
                scale(
                    coefficient,
                    reduce_delayed_word(base_word, maximum_delay),
                ),
            )
        return result
    word = next(iter(base))

    for delay in range(1, maximum_delay + 1):
        for letter, other in (("s", "a"), ("a", "s")):
            pattern = other + letter * (delay + 2) + other
            position = word.find(pattern)
            if position < 0:
                continue
            prefix = {word[:position]: Fraction(1)}
            suffix = {
                word[position + len(pattern) :]: Fraction(1)
            }
            replacement = add(
                scale(-1, {letter * delay: Fraction(1)}),
                {letter * (delay + 1) + other: Fraction(1)},
                {other + letter * (delay + 1): Fraction(1)},
            )
            expanded = multiply(
                multiply(prefix, replacement),
                suffix,
            )
            result: Polynomial = {}
            for new_word, coefficient in expanded.items():
                result = add(
                    result,
                    scale(
                        coefficient,
                        reduce_delayed_word(
                            new_word,
                            maximum_delay,
                        ),
                    ),
                )
            return result
    return {word: Fraction(1)}


def reduce_delay(
    polynomial: Polynomial,
    maximum_delay: int,
) -> Polynomial:
    """Reduce a polynomial modulo the complete-delay relations."""

    result: Polynomial = {}
    for word, coefficient in polynomial.items():
        result = add(
            result,
            scale(
                coefficient,
                reduce_delayed_word(word, maximum_delay),
            ),
        )
    return result


def exact_records(
    maximum_grade: int,
) -> list[DelayedSlackAnticommutatorRecord]:
    """Return exact word-algebra audits through one maximum grade."""

    residual = schur_residual_series(2 * maximum_grade)
    records: list[DelayedSlackAnticommutatorRecord] = []
    for grade in range(1, maximum_grade + 1):
        face_degree = 2 * grade
        maximum_earlier = max(
            (
                len(reduce_delay(residual[degree], grade - 1))
                for degree in range(face_degree)
            ),
            default=0,
        )
        difference = add(
            residual[face_degree],
            scale(-1, target(grade)),
        )
        reduced_difference = reduce_delay(difference, grade - 1)
        verified = maximum_earlier == 0 and not reduced_difference
        if not verified:
            raise RuntimeError(
                "the exact delayed-slack audit failed: "
                f"grade={grade}, earlier={maximum_earlier}, "
                f"face={len(reduced_difference)}"
            )
        records.append(
            DelayedSlackAnticommutatorRecord(
                grade=grade,
                face_degree=face_degree,
                face_word_count=len(residual[face_degree]),
                unreduced_face_difference_word_count=len(difference),
                maximum_earlier_reduced_word_count=maximum_earlier,
                reduced_face_difference_word_count=len(
                    reduced_difference
                ),
                all_checks_passed=verified,
            )
        )
    return records


def write_records(
    records: list[DelayedSlackAnticommutatorRecord],
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
    parser.add_argument(
        "--maximum-grade",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_delayed_slack_anticommutator_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact finite-grade audit."""

    args = parse_args()
    records = exact_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
