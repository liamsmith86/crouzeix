#!/usr/bin/env python3
"""Audit the cyclic radial form of the delayed edge-deleted volume.

L252 proves that the radial trace

    4 tr(Q_(k+2) - (k+2) Q_1 + (k+1) I)

equals ``4 ||B_k||_F**2`` under a complete delay.  This checker
constructs L250's edge-deleted normalized mass in exact word algebra
and verifies that its active trace reduces to that radial expression
through finite grades.  It is an exact finite-grade audit, not the
missing arbitrary-grade cyclic proof.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    Series,
    add,
    boundary_metric_series,
    multiply,
    operator_series,
    reduce_delayed_word,
    scale,
    series_add,
    series_adjoint,
    series_multiply,
    zero_series,
)


@dataclass(frozen=True)
class CyclicRadialVolumeRecord:
    """One exact cyclic radial audit."""

    grade: int
    face_degree: int
    raw_face_word_count: int
    cyclic_face_word_count: int
    radial_difference_word_count: int
    transfer_difference_word_count: int
    all_checks_passed: bool


def inverse_series(series: Series) -> Series:
    """Invert a word-polynomial series with identity constant term."""

    maximum_degree = len(series) - 1
    result = zero_series(maximum_degree)
    result[0] = IDENTITY
    for degree in range(1, maximum_degree + 1):
        convolution: Polynomial = {}
        for positive_degree in range(1, degree + 1):
            convolution = add(
                convolution,
                multiply(
                    series[positive_degree],
                    result[degree - positive_degree],
                ),
            )
        result[degree] = scale(-1, convolution)
    return result


@lru_cache(maxsize=500_000)
def cyclic_reduce_word(word: str, maximum_delay: int) -> Polynomial:
    """Reduce one word modulo algebra, delay, and trace cyclicity.

    A cyclic rotation is used as a rewrite only when the existing
    algebra/delay reducer shortens every resulting word.  This makes
    termination explicit and avoids assuming a new trace identity.
    """

    reduced = reduce_delayed_word(word, maximum_delay)
    if reduced != {word: Fraction(1)}:
        result: Polynomial = {}
        for new_word, coefficient in reduced.items():
            result = add(
                result,
                scale(
                    coefficient,
                    cyclic_reduce_word(new_word, maximum_delay),
                ),
            )
        return result

    if not word:
        return IDENTITY

    rotations = [
        word[position:] + word[:position]
        for position in range(len(word))
    ]
    for rotation in rotations:
        rotated_reduction = reduce_delayed_word(
            rotation,
            maximum_delay,
        )
        if rotated_reduction == {rotation: Fraction(1)}:
            continue
        if not all(
            len(new_word) < len(word)
            for new_word in rotated_reduction
        ):
            raise RuntimeError(
                "cyclic reduction did not strictly shorten a word"
            )
        result: Polynomial = {}
        for new_word, coefficient in rotated_reduction.items():
            result = add(
                result,
                scale(
                    coefficient,
                    cyclic_reduce_word(new_word, maximum_delay),
                ),
            )
        return result

    return {min(rotations): Fraction(1)}


def cyclic_reduce(
    polynomial: Polynomial,
    maximum_delay: int,
) -> Polynomial:
    """Reduce a polynomial in the cyclic delayed quotient."""

    result: Polynomial = {}
    for word, coefficient in polynomial.items():
        result = add(
            result,
            scale(
                coefficient,
                cyclic_reduce_word(word, maximum_delay),
            ),
        )
    return result


def radial_target(grade: int) -> Polynomial:
    """Return four times the L252 three-term radial expression."""

    q_one = multiply({"a": Fraction(1)}, {"s": Fraction(1)})
    q_top = multiply(
        {"a" * (grade + 2): Fraction(1)},
        {"s" * (grade + 2): Fraction(1)},
    )
    return add(
        scale(4, q_top),
        scale(-4 * (grade + 2), q_one),
        scale(4 * (grade + 1), IDENTITY),
    )


def transfer_target(grade: int) -> Polynomial:
    """Return four times the state trace representative of ||B_k||²."""

    right_orbit = multiply(
        multiply(
            {"a" * grade: Fraction(1)},
            E,
        ),
        {"s" * grade: Fraction(1)},
    )
    return scale(4, multiply(right_orbit, F))


def mass_face(grade: int) -> Polynomial:
    """Return L250's active edge-deleted mass coefficient."""

    maximum_degree = 2 * grade
    operator = operator_series(maximum_degree)
    metric = boundary_metric_series(maximum_degree)
    metric[maximum_degree] = {}

    final_weight = [F] + [{} for _ in range(maximum_degree)]
    row_gram = series_multiply(
        series_multiply(
            series_adjoint(operator),
            final_weight,
        ),
        operator,
    )
    row_denominator = series_add(
        metric,
        [scale(-1, coefficient) for coefficient in row_gram],
    )

    pulled_metric = series_multiply(
        series_multiply(
            series_adjoint(operator),
            metric,
        ),
        operator,
    )
    stein_slack = series_add(
        metric,
        [scale(-1, coefficient) for coefficient in pulled_metric],
    )
    normalized_mass = series_multiply(
        inverse_series(row_denominator),
        stein_slack,
    )
    return normalized_mass[maximum_degree]


def audit_grade(grade: int) -> CyclicRadialVolumeRecord:
    """Audit one finite grade exactly."""

    face = mass_face(grade)
    maximum_delay = grade - 1
    cyclic_face = cyclic_reduce(face, maximum_delay)
    radial = cyclic_reduce(radial_target(grade), maximum_delay)
    transfer = cyclic_reduce(
        transfer_target(grade),
        maximum_delay,
    )
    radial_difference = add(cyclic_face, scale(-1, radial))
    transfer_difference = add(cyclic_face, scale(-1, transfer))
    verified = not radial_difference and not transfer_difference
    if not verified:
        raise RuntimeError(
            "cyclic radial volume audit failed: "
            f"grade={grade}, radial={len(radial_difference)}, "
            f"transfer={len(transfer_difference)}"
        )
    return CyclicRadialVolumeRecord(
        grade=grade,
        face_degree=2 * grade,
        raw_face_word_count=len(face),
        cyclic_face_word_count=len(cyclic_face),
        radial_difference_word_count=len(radial_difference),
        transfer_difference_word_count=len(transfer_difference),
        all_checks_passed=verified,
    )


def standard_records(
    maximum_grade: int,
) -> list[CyclicRadialVolumeRecord]:
    """Return exact audits through one maximum grade."""

    return [
        audit_grade(grade)
        for grade in range(1, maximum_grade + 1)
    ]


def write_records(
    records: list[CyclicRadialVolumeRecord],
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
    parser.add_argument("--maximum-grade", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_cyclic_radial_volume_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact audits."""

    args = parse_args()
    records = standard_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
