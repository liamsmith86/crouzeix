#!/usr/bin/env python3
"""Audit the all-grade induction for the unweighted return fan."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from repeated_crabb_cyclic_radial_volume import DelayedQuotient
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    S,
    STAR,
)
from repeated_crabb_remote_endpoint_coefficients import (
    RenewalPieces,
    renewal_pieces,
)


@dataclass(frozen=True)
class FullFanRecord:
    """One exact fan-restricted full-minus-tail audit."""

    record_type: str
    grade: int
    first_family_inherited_maximum: int
    second_family_inherited_maximum: int
    new_first_coefficient: int
    new_second_coefficient: int
    deep_coefficient: int
    shifted_deep_coefficient: int
    full_fan_word_count: int
    full_fan_coefficients_are_one: bool
    full_additional_nonradial_word_count: int
    layer_additional_nonradial_word_count: int
    all_checks_passed: bool


@dataclass(frozen=True)
class FanEmbeddingRecord:
    """One exact tail-to-full fan embedding audit."""

    record_type: str
    grade: int
    first_family_embeddings_pass: bool
    second_family_embeddings_pass: bool
    deep_embedding_passes: bool
    all_checks_passed: bool


def coefficient(polynomial: Polynomial, word: str) -> int:
    """Return one integral coefficient."""

    value = polynomial.get(word, Fraction())
    if value.denominator != 1:
        raise RuntimeError(f"nonintegral coefficient for {word}: {value}")
    return value.numerator


def add_pieces(
    quotient: DelayedQuotient,
    pieces: RenewalPieces,
    degree: int,
) -> Polynomial:
    """Add the direct and returned coefficients of one renewal."""

    return quotient.add(
        pieces.direct[degree],
        pieces.immediate_return[degree],
        pieces.later_returns[degree],
    )


def first_fan_word(grade: int, index: int) -> str:
    """Return the indexed first-family fan word."""

    return (
        "a" * index
        + "s" * (grade + 2)
        + "a" * (grade + 2 - index)
    )


def second_fan_word(grade: int, index: int) -> str:
    """Return the indexed second-family fan word."""

    return (
        "s" * (index + 2)
        + "a" * (grade + 2)
        + "s" * (grade - index)
    )


def deep_word(grade: int) -> str:
    """Return the deep fan endpoint."""

    return (
        "s"
        + "a" * (grade + 2)
        + "s" * (grade + 2)
        + "a"
    )


def full_fan_words(grade: int) -> list[str]:
    """Return the complete unweighted fan word list."""

    return [
        *(first_fan_word(grade, index) for index in range(grade + 1)),
        *(second_fan_word(grade, index) for index in range(grade)),
        deep_word(grade),
    ]


def is_radial_word(word: str) -> bool:
    """Return whether a reduced word is an identity, Q, or R power."""

    if not word:
        return True
    half_length, remainder = divmod(len(word), 2)
    if remainder:
        return False
    return word in {
        "a" * half_length + "s" * half_length,
        "s" * half_length + "a" * half_length,
    }


def audit_grade(grade: int) -> FullFanRecord:
    """Audit one fan-restricted associated recursion."""

    quotient = DelayedQuotient(grade)
    full_degree = 2 * grade
    full_pieces = renewal_pieces(
        quotient,
        S,
        E,
        F,
        IDENTITY,
        full_degree,
    )
    full = add_pieces(quotient, full_pieces, full_degree)

    if grade == 1:
        delta = full
    else:
        tail_identity = quotient.add(
            IDENTITY,
            quotient.scale(-1, F),
        )
        tail_partial = quotient.multiply(S, tail_identity)
        tail_final = quotient.multiply(
            quotient.multiply(S, F),
            STAR,
        )
        tail_degree = 2 * (grade - 1)
        tail_pieces = renewal_pieces(
            quotient,
            tail_partial,
            E,
            tail_final,
            tail_identity,
            tail_degree,
        )
        tail = add_pieces(quotient, tail_pieces, tail_degree)
        delta = quotient.add(full, quotient.scale(-1, tail))

    inherited_first = [
        abs(coefficient(delta, first_fan_word(grade, index)))
        for index in range(grade)
    ]
    inherited_second = [
        abs(coefficient(delta, second_fan_word(grade, index)))
        for index in range(1, grade)
    ]
    new_first = coefficient(delta, first_fan_word(grade, grade))
    new_second = coefficient(delta, second_fan_word(grade, 0))
    deep = deep_word(grade)
    deep_value = coefficient(delta, deep)
    shifted_deep_value = coefficient(delta, "s" + deep + "a")
    fan_coefficients = [
        coefficient(full, word)
        for word in full_fan_words(grade)
    ]
    full_allowed = set(full_fan_words(grade))
    full_additional = {
        word
        for word in full
        if word not in full_allowed and not is_radial_word(word)
    }
    if grade == 1:
        layer_allowed = full_allowed
    else:
        layer_allowed = {
            first_fan_word(grade, grade),
            second_fan_word(grade, 0),
            deep,
            "s" + deep + "a",
        }
    layer_additional = {
        word
        for word in delta
        if word not in layer_allowed and not is_radial_word(word)
    }

    if grade == 1:
        expected_inherited_first = 1
        expected_inherited_second = 0
        expected_shifted_deep = 0
    else:
        expected_inherited_first = 0
        expected_inherited_second = 0
        expected_shifted_deep = -1
    first_maximum = max(inherited_first, default=0)
    second_maximum = max(inherited_second, default=0)
    verified = bool(
        first_maximum == expected_inherited_first
        and second_maximum == expected_inherited_second
        and new_first == 1
        and new_second == 1
        and deep_value == 1
        and shifted_deep_value == expected_shifted_deep
        and all(value == 1 for value in fan_coefficients)
        and not full_additional
        and not layer_additional
    )
    if not verified:
        raise RuntimeError(
            "the full-fan audit failed: "
            f"grade={grade}, inherited=({first_maximum},"
            f"{second_maximum}), new=({new_first},{new_second}), "
            f"deep=({deep_value},{shifted_deep_value}), "
            f"fan={fan_coefficients}, "
            f"additional=({sorted(full_additional)},"
            f"{sorted(layer_additional)})"
        )

    return FullFanRecord(
        record_type="physical_fan",
        grade=grade,
        first_family_inherited_maximum=first_maximum,
        second_family_inherited_maximum=second_maximum,
        new_first_coefficient=new_first,
        new_second_coefficient=new_second,
        deep_coefficient=deep_value,
        shifted_deep_coefficient=shifted_deep_value,
        full_fan_word_count=len(fan_coefficients),
        full_fan_coefficients_are_one=all(
            value == 1 for value in fan_coefficients
        ),
        full_additional_nonradial_word_count=len(full_additional),
        layer_additional_nonradial_word_count=len(layer_additional),
        all_checks_passed=verified,
    )


def fan_polynomials(
    quotient: DelayedQuotient,
    partial: Polynomial,
    identity: Polynomial,
    grade: int,
) -> tuple[list[Polynomial], list[Polynomial], Polynomial]:
    """Return both fan families and the deep word for one partial isometry."""

    adjoint = quotient.polynomial_adjoint(partial)
    first = [
        quotient.multiply(
            quotient.multiply(
                quotient.polynomial_power(adjoint, index, identity),
                quotient.polynomial_power(
                    partial,
                    grade + 2,
                    identity,
                ),
            ),
            quotient.polynomial_power(
                adjoint,
                grade + 2 - index,
                identity,
            ),
        )
        for index in range(grade + 1)
    ]
    second = [
        quotient.multiply(
            quotient.multiply(
                quotient.polynomial_power(
                    partial,
                    index + 2,
                    identity,
                ),
                quotient.polynomial_power(
                    adjoint,
                    grade + 2,
                    identity,
                ),
            ),
            quotient.polynomial_power(
                partial,
                grade - index,
                identity,
            ),
        )
        for index in range(grade)
    ]
    deep = quotient.multiply(
        quotient.multiply(
            quotient.multiply(
                partial,
                quotient.polynomial_power(
                    adjoint,
                    grade + 2,
                    identity,
                ),
            ),
            quotient.polynomial_power(
                partial,
                grade + 2,
                identity,
            ),
        ),
        adjoint,
    )
    return first, second, deep


def audit_embedding(grade: int) -> FanEmbeddingRecord:
    """Audit the one-layer embedding of the complete tail fan."""

    if grade < 2:
        raise ValueError("tail fan embedding starts at grade two")

    quotient = DelayedQuotient(grade)
    tail_identity = quotient.add(
        IDENTITY,
        quotient.scale(-1, F),
    )
    tail_partial = quotient.multiply(S, tail_identity)
    tail_first, tail_second, tail_deep = fan_polynomials(
        quotient,
        tail_partial,
        tail_identity,
        grade - 1,
    )
    full_first, full_second, full_deep = fan_polynomials(
        quotient,
        S,
        IDENTITY,
        grade,
    )

    first_passes = tail_first == full_first[:-1]
    second_passes = tail_second == full_second[1:]
    shifted_deep = quotient.multiply(
        quotient.multiply(S, full_deep),
        STAR,
    )
    deep_passes = tail_deep == shifted_deep
    verified = bool(first_passes and second_passes and deep_passes)
    if not verified:
        raise RuntimeError(
            "the fan embedding audit failed: "
            f"grade={grade}, first={first_passes}, "
            f"second={second_passes}, deep={deep_passes}"
        )
    return FanEmbeddingRecord(
        record_type="tail_embedding",
        grade=grade,
        first_family_embeddings_pass=first_passes,
        second_family_embeddings_pass=second_passes,
        deep_embedding_passes=deep_passes,
        all_checks_passed=verified,
    )


def exact_records(
    maximum_physical_grade: int,
    maximum_embedding_grade: int,
) -> list[FullFanRecord | FanEmbeddingRecord]:
    """Return physical coefficient and algebraic embedding audits."""

    return [
        *(audit_grade(grade) for grade in range(1, maximum_physical_grade + 1)),
        *(
            audit_embedding(grade)
            for grade in range(2, maximum_embedding_grade + 1)
        ),
    ]


def write_records(
    records: list[FullFanRecord | FanEmbeddingRecord],
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
    parser.add_argument("--maximum-physical-grade", type=int, default=6)
    parser.add_argument("--maximum-embedding-grade", type=int, default=12)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_full_fan_coefficients_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact full-fan audits."""

    args = parse_args()
    records = exact_records(
        args.maximum_physical_grade,
        args.maximum_embedding_grade,
    )
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
