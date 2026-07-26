#!/usr/bin/env python3
"""Audit the four remote coefficients in the unweighted delay recursion."""

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


@dataclass(frozen=True)
class RemoteEndpointRecord:
    """Exact source coefficients for one delayed grade."""

    grade: int
    first_crossed_direct: int
    first_crossed_immediate_return: int
    first_crossed_later_returns: int
    second_crossed_direct: int
    second_crossed_immediate_return: int
    second_crossed_later_returns: int
    deep_full_direct: int
    deep_full_all_returns: int
    shifted_tail_direct_difference: int
    shifted_tail_all_returns_difference: int
    all_checks_passed: bool


@dataclass(frozen=True)
class RenewalPieces:
    """The three pieces of ``I - Z_ret`` relevant to the audit."""

    direct: list[Polynomial]
    immediate_return: list[Polynomial]
    later_returns: list[Polynomial]


def subtract_series(
    quotient: DelayedQuotient,
    left: list[Polynomial],
    right: list[Polynomial],
) -> list[Polynomial]:
    """Subtract two quotient series coefficientwise."""

    return [
        quotient.add(left_item, quotient.scale(-1, right_item))
        for left_item, right_item in zip(left, right, strict=True)
    ]


def renewal_pieces(
    quotient: DelayedQuotient,
    partial: Polynomial,
    initial: Polynomial,
    final: Polynomial,
    identity: Polynomial,
    active_degree: int,
) -> RenewalPieces:
    """Return direct, zero-loop, and positive-loop unweighted pieces."""

    operator = quotient.operator_series_for(
        partial,
        initial,
        final,
        identity,
    )
    metric = quotient.metric_series_for(
        partial,
        initial,
        final,
        identity,
    )
    metric[active_degree] = {}
    output_gram = quotient.series_multiply(
        quotient.series_multiply(
            operator,
            quotient.inverse_series(metric, identity),
        ),
        quotient.series_adjoint(operator),
    )
    retained = quotient.add(
        identity,
        quotient.scale(-1, final),
    )
    retained_gram = quotient.compress_series(
        output_gram,
        retained,
        retained,
    )
    entrance = quotient.compress_series(
        output_gram,
        retained,
        final,
    )
    exit_ = quotient.compress_series(
        output_gram,
        final,
        retained,
    )
    loop = quotient.compress_series(
        output_gram,
        final,
        final,
    )
    loop_denominator = quotient.series_subtract(
        quotient.constant_series(identity),
        loop,
    )
    all_returns = quotient.series_multiply(
        quotient.series_multiply(
            entrance,
            quotient.inverse_series(loop_denominator, identity),
        ),
        exit_,
    )
    immediate_return = quotient.series_multiply(entrance, exit_)
    later_returns = subtract_series(
        quotient,
        all_returns,
        immediate_return,
    )
    return RenewalPieces(
        direct=[
            quotient.scale(-1, coefficient)
            for coefficient in retained_gram
        ],
        immediate_return=[
            quotient.scale(-1, coefficient)
            for coefficient in immediate_return
        ],
        later_returns=[
            quotient.scale(-1, coefficient)
            for coefficient in later_returns
        ],
    )


def coefficient(
    series: list[Polynomial],
    degree: int,
    word: str,
) -> int:
    """Return one integral word coefficient."""

    value = series[degree].get(word, Fraction())
    if value.denominator != 1:
        raise RuntimeError(f"nonintegral coefficient for {word}: {value}")
    return value.numerator


def audit_grade(grade: int) -> RemoteEndpointRecord:
    """Audit the four universal remote coefficients at one grade."""

    if grade < 2:
        raise ValueError("the remote recursion starts at grade two")

    quotient = DelayedQuotient(grade)
    full_degree = 2 * grade
    full = renewal_pieces(
        quotient,
        S,
        E,
        F,
        IDENTITY,
        full_degree,
    )

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
    tail = renewal_pieces(
        quotient,
        tail_partial,
        E,
        tail_final,
        tail_identity,
        tail_degree,
    )

    first_crossed = (
        "a" * grade
        + "s" * (grade + 2)
        + "a" * 2
    )
    second_crossed = (
        "s" * 2
        + "a" * (grade + 2)
        + "s" * grade
    )
    deep = (
        "s"
        + "a" * (grade + 2)
        + "s" * (grade + 2)
        + "a"
    )
    shifted_deep = "s" + deep + "a"

    sign = (-1) ** grade
    first_values = (
        coefficient(full.direct, full_degree, first_crossed),
        coefficient(
            full.immediate_return,
            full_degree,
            first_crossed,
        ),
        coefficient(full.later_returns, full_degree, first_crossed),
    )
    second_values = (
        coefficient(full.direct, full_degree, second_crossed),
        coefficient(
            full.immediate_return,
            full_degree,
            second_crossed,
        ),
        coefficient(full.later_returns, full_degree, second_crossed),
    )
    deep_values = (
        coefficient(full.direct, full_degree, deep),
        coefficient(full.immediate_return, full_degree, deep)
        + coefficient(full.later_returns, full_degree, deep),
    )
    shifted_values = (
        -coefficient(tail.direct, tail_degree, shifted_deep),
        -coefficient(
            tail.immediate_return,
            tail_degree,
            shifted_deep,
        )
        - coefficient(tail.later_returns, tail_degree, shifted_deep),
    )

    expected_crossed = (1 + 4 * sign, -4 * sign, 0)
    verified = bool(
        first_values == expected_crossed
        and second_values == expected_crossed
        and deep_values == (1, 0)
        and shifted_values == (-1, 0)
    )
    if not verified:
        raise RuntimeError(
            "the remote endpoint audit failed: "
            f"grade={grade}, first={first_values}, "
            f"second={second_values}, deep={deep_values}, "
            f"shifted={shifted_values}"
        )

    return RemoteEndpointRecord(
        grade=grade,
        first_crossed_direct=first_values[0],
        first_crossed_immediate_return=first_values[1],
        first_crossed_later_returns=first_values[2],
        second_crossed_direct=second_values[0],
        second_crossed_immediate_return=second_values[1],
        second_crossed_later_returns=second_values[2],
        deep_full_direct=deep_values[0],
        deep_full_all_returns=deep_values[1],
        shifted_tail_direct_difference=shifted_values[0],
        shifted_tail_all_returns_difference=shifted_values[1],
        all_checks_passed=verified,
    )


def exact_records(maximum_grade: int) -> list[RemoteEndpointRecord]:
    """Return exact audits from grade two through one maximum."""

    return [
        audit_grade(grade)
        for grade in range(2, maximum_grade + 1)
    ]


def write_records(
    records: list[RemoteEndpointRecord],
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
            "repeated_crabb_remote_endpoint_coefficients_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact remote-coefficient audits."""

    args = parse_args()
    records = exact_records(args.maximum_grade)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
