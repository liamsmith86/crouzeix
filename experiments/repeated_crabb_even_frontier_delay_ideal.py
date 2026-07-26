#!/usr/bin/env python3
"""Audit the intact-metric even frontier and its delay-ideal factorization."""

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
    reduce_delay,
)
from repeated_crabb_odd_frontier_delay_ideal import (
    delay_generator,
    product,
)


@dataclass(frozen=True)
class EvenFrontierRecord:
    """One exact physical intact-metric frontier audit."""

    record_type: str
    delay: int
    face_degree: int
    frontier_word_count: int
    formula_difference_word_count: int
    next_delay_remainder_word_count: int
    tail_formula_difference_word_count: int
    associated_layer_word_count: int
    all_checks_passed: bool


@dataclass(frozen=True)
class EvenFactorizationRecord:
    """One algebraic factorization and tail-embedding audit."""

    record_type: str
    delay: int
    frontier_word_count: int
    factorization_difference_word_count: int
    next_delay_remainder_word_count: int
    tail_embedding_difference_word_count: int
    all_checks_passed: bool


def even_frontier_formula(
    quotient: DelayedQuotient,
    delay: int,
    partial: Polynomial,
    identity: Polynomial,
) -> Polynomial:
    """Return the universal five-word intact-metric even frontier."""

    adjoint = quotient.polynomial_adjoint(partial)

    def power(polynomial: Polynomial, exponent: int) -> Polynomial:
        return quotient.polynomial_power(
            polynomial,
            exponent,
            identity,
        )

    return quotient.add(
        quotient.scale(2, identity),
        quotient.scale(
            -2,
            product(quotient, adjoint, partial),
        ),
        quotient.scale(
            -2,
            product(
                quotient,
                power(partial, delay + 1),
                power(adjoint, delay + 1),
            ),
        ),
        product(
            quotient,
            adjoint,
            power(partial, delay + 2),
            power(adjoint, delay + 1),
        ),
        product(
            quotient,
            power(partial, delay + 1),
            power(adjoint, delay + 2),
            partial,
        ),
    )


def even_factorization(
    quotient: DelayedQuotient,
    delay: int,
) -> Polynomial:
    """Return ``D_h (S*)^h + S^h D_h*``."""

    forward = delay_generator(quotient, delay, "s")
    backward = delay_generator(quotient, delay, "a")
    return quotient.add(
        product(
            quotient,
            forward,
            {"a" * delay: Fraction(1)},
        ),
        product(
            quotient,
            {"s" * delay: Fraction(1)},
            backward,
        ),
    )


def tail_data(
    quotient: DelayedQuotient,
) -> tuple[Polynomial, Polynomial, Polynomial]:
    """Return the independently balanced tail partial, final, and identity."""

    identity = quotient.add(
        IDENTITY,
        quotient.scale(-1, F),
    )
    partial = quotient.multiply(S, identity)
    final = quotient.multiply(
        quotient.multiply(S, F),
        STAR,
    )
    return partial, final, identity


def intact_defect_coefficient(
    quotient: DelayedQuotient,
    partial: Polynomial,
    initial: Polynomial,
    final: Polynomial,
    identity: Polynomial,
    degree: int,
    cutoff: int,
) -> Polynomial:
    """Return one coefficient while deleting only the later cutoff face."""

    return quotient.closed_return_defect_for(
        partial,
        initial,
        final,
        identity,
        cutoff,
    )[degree]


def physical_records(
    maximum_delay: int,
) -> list[EvenFrontierRecord]:
    """Return exact physical audits through one finite delay."""

    records: list[EvenFrontierRecord] = []
    for delay in range(1, maximum_delay + 1):
        quotient = DelayedQuotient(delay)
        degree = 2 * delay
        cutoff = degree + 2
        quotient.maximum_degree = cutoff
        actual = intact_defect_coefficient(
            quotient,
            S,
            E,
            F,
            IDENTITY,
            degree,
            cutoff,
        )
        expected = even_frontier_formula(
            quotient,
            delay,
            S,
            IDENTITY,
        )
        formula_difference = quotient.add(
            actual,
            quotient.scale(-1, expected),
        )
        delayed_remainder = reduce_delay(expected, delay)

        tail_difference: Polynomial = {}
        associated_layer: Polynomial = {}
        if delay >= 2:
            tail_partial, tail_final, tail_identity = tail_data(
                quotient
            )
            tail_degree = degree - 2
            actual_tail = intact_defect_coefficient(
                quotient,
                tail_partial,
                E,
                tail_final,
                tail_identity,
                tail_degree,
                cutoff,
            )
            expected_tail = even_frontier_formula(
                quotient,
                delay - 1,
                tail_partial,
                tail_identity,
            )
            tail_difference = quotient.add(
                actual_tail,
                quotient.scale(-1, expected_tail),
            )
            associated_layer = quotient.add(
                actual,
                quotient.scale(-1, actual_tail),
            )

        verified = bool(
            not formula_difference
            and not delayed_remainder
            and not tail_difference
            and not associated_layer
        )
        if not verified:
            raise RuntimeError(
                "even-frontier physical audit failed: "
                f"delay={delay}, formula={formula_difference}, "
                f"delayed={delayed_remainder}, tail={tail_difference}, "
                f"layer={associated_layer}"
            )
        records.append(
            EvenFrontierRecord(
                record_type="physical_frontier",
                delay=delay,
                face_degree=degree,
                frontier_word_count=len(actual),
                formula_difference_word_count=len(formula_difference),
                next_delay_remainder_word_count=len(delayed_remainder),
                tail_formula_difference_word_count=len(tail_difference),
                associated_layer_word_count=len(associated_layer),
                all_checks_passed=verified,
            )
        )
    return records


def factorization_records(
    maximum_delay: int,
) -> list[EvenFactorizationRecord]:
    """Return exact algebraic audits through one larger delay."""

    records: list[EvenFactorizationRecord] = []
    for delay in range(1, maximum_delay + 1):
        quotient = DelayedQuotient(delay)
        quotient.maximum_degree = 2 * delay + 2
        frontier = even_frontier_formula(
            quotient,
            delay,
            S,
            IDENTITY,
        )
        factorization_difference = quotient.add(
            frontier,
            quotient.scale(
                -1,
                even_factorization(quotient, delay),
            ),
        )
        delayed_remainder = reduce_delay(frontier, delay)

        tail_embedding_difference: Polynomial = {}
        if delay >= 2:
            tail_partial, _, tail_identity = tail_data(quotient)
            embedded_tail = even_frontier_formula(
                quotient,
                delay - 1,
                tail_partial,
                tail_identity,
            )
            tail_embedding_difference = quotient.add(
                frontier,
                quotient.scale(-1, embedded_tail),
            )

        verified = bool(
            not factorization_difference
            and not delayed_remainder
            and not tail_embedding_difference
        )
        if not verified:
            raise RuntimeError(
                "even-frontier factorization audit failed: "
                f"delay={delay}, factor={factorization_difference}, "
                f"delayed={delayed_remainder}, "
                f"tail={tail_embedding_difference}"
            )
        records.append(
            EvenFactorizationRecord(
                record_type="algebraic_factorization",
                delay=delay,
                frontier_word_count=len(frontier),
                factorization_difference_word_count=len(
                    factorization_difference
                ),
                next_delay_remainder_word_count=len(delayed_remainder),
                tail_embedding_difference_word_count=len(
                    tail_embedding_difference
                ),
                all_checks_passed=verified,
            )
        )
    return records


def write_records(
    records: list[EvenFrontierRecord | EvenFactorizationRecord],
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
    parser.add_argument("--maximum-physical-delay", type=int, default=4)
    parser.add_argument("--maximum-factor-delay", type=int, default=12)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_even_frontier_delay_ideal_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact audits."""

    args = parse_args()
    records: list[EvenFrontierRecord | EvenFactorizationRecord] = [
        *physical_records(args.maximum_physical_delay),
        *factorization_records(args.maximum_factor_delay),
    ]
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
