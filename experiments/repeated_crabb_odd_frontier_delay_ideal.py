#!/usr/bin/env python3
"""Audit the all-grade odd return frontier and its delay-ideal factorization."""

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
from repeated_crabb_remote_endpoint_coefficients import renewal_pieces


@dataclass(frozen=True)
class OddFrontierRecord:
    """One exact physical odd-frontier audit."""

    record_type: str
    delay: int
    face_degree: int
    frontier_word_count: int
    formula_difference_word_count: int
    next_delay_remainder_word_count: int
    tail_formula_difference_word_count: int
    layer_formula_difference_word_count: int
    all_checks_passed: bool


@dataclass(frozen=True)
class OddFactorizationRecord:
    """One algebraic factorization and tail-recursion audit."""

    record_type: str
    delay: int
    frontier_word_count: int
    factorization_difference_word_count: int
    next_delay_remainder_word_count: int
    tail_recursion_difference_word_count: int
    all_checks_passed: bool


def monomial(letter: str, exponent: int) -> Polynomial:
    """Return one word monomial."""

    return {letter * exponent: Fraction(1)}


def product(
    quotient: DelayedQuotient,
    *polynomials: Polynomial,
) -> Polynomial:
    """Multiply a finite list in the delayed quotient."""

    result = IDENTITY
    for polynomial in polynomials:
        result = quotient.multiply(result, polynomial)
    return result


def odd_frontier_formula(
    quotient: DelayedQuotient,
    delay: int,
    partial: Polynomial,
    identity: Polynomial,
) -> Polynomial:
    """Return the closed twelve-word odd frontier for delay at least two."""

    if delay < 2:
        raise ValueError("the twelve-word formula starts at delay two")

    adjoint = quotient.polynomial_adjoint(partial)
    outer = delay + 2

    def power(polynomial: Polynomial, exponent: int) -> Polynomial:
        return quotient.polynomial_power(
            polynomial,
            exponent,
            identity,
        )

    return quotient.add(
        quotient.scale(
            2,
            product(
                quotient,
                power(adjoint, delay - 1),
                power(partial, outer),
                adjoint,
            ),
        ),
        quotient.scale(
            2,
            product(
                quotient,
                partial,
                power(adjoint, outer),
                power(partial, delay - 1),
            ),
        ),
        product(
            quotient,
            power(partial, delay),
            power(adjoint, outer),
        ),
        product(
            quotient,
            power(partial, outer),
            power(adjoint, delay),
        ),
        product(
            quotient,
            power(partial, delay + 1),
            power(adjoint, outer + 1),
        ),
        product(
            quotient,
            power(partial, outer + 1),
            power(adjoint, delay + 1),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                power(adjoint, 2),
                power(partial, outer),
                power(adjoint, outer),
            ),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                adjoint,
                power(partial, outer),
                power(adjoint, outer + 1),
            ),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                power(partial, outer),
                power(adjoint, outer),
                power(partial, 2),
            ),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                power(partial, outer + 1),
                power(adjoint, outer),
                partial,
            ),
        ),
        quotient.scale(
            -2,
            product(
                quotient,
                partial,
                power(adjoint, outer),
                power(partial, outer),
                power(adjoint, 3),
            ),
        ),
        quotient.scale(
            -2,
            product(
                quotient,
                power(partial, 3),
                power(adjoint, outer),
                power(partial, outer),
                adjoint,
            ),
        ),
    )


def first_frontier_formula(
    quotient: DelayedQuotient,
) -> Polynomial:
    """Return the exceptional degree-three frontier."""

    return quotient.add(
        quotient.scale(3, {"sssa": Fraction(1)}),
        {"ssssaa": Fraction(1)},
        quotient.scale(-3, {"sssaaasssa": Fraction(1)}),
        quotient.scale(3, {"saaa": Fraction(1)}),
        quotient.scale(-1, {"asssaaaa": Fraction(1)}),
        quotient.scale(-3, {"saaasssaaa": Fraction(1)}),
        {"ssaaaa": Fraction(1)},
        quotient.scale(-1, {"ssssaaas": Fraction(1)}),
    )


def delay_generator(
    quotient: DelayedQuotient,
    delay: int,
    letter: str,
) -> Polynomial:
    """Return ``E S^delay F`` or its adjoint in base-normal form."""

    other = "a" if letter == "s" else "s"
    return quotient.add(
        {other + letter * (delay + 2) + other: Fraction(1)},
        {letter * delay: Fraction(1)},
        quotient.scale(
            -1,
            {letter * (delay + 1) + other: Fraction(1)},
        ),
        quotient.scale(
            -1,
            {other + letter * (delay + 1): Fraction(1)},
        ),
    )


def odd_factorization(
    quotient: DelayedQuotient,
    delay: int,
) -> Polynomial:
    """Return the explicit two-sided delay-ideal factorization."""

    forward = delay_generator(quotient, delay, "s")
    backward = delay_generator(quotient, delay, "a")
    if delay == 1:
        return quotient.add(
            quotient.scale(
                -3,
                product(
                    quotient,
                    S,
                    monomial("a", 2),
                    forward,
                    monomial("a", 2),
                ),
            ),
            quotient.scale(
                -3,
                product(
                    quotient,
                    monomial("s", 3),
                    monomial("a", 2),
                    forward,
                ),
            ),
            quotient.scale(
                -1,
                product(
                    quotient,
                    forward,
                    monomial("a", 3),
                ),
            ),
            quotient.scale(
                -3,
                product(
                    quotient,
                    monomial("s", 2),
                    backward,
                    S,
                ),
            ),
            quotient.scale(
                -1,
                product(
                    quotient,
                    monomial("s", 3),
                    backward,
                ),
            ),
        )

    return quotient.add(
        quotient.scale(
            -2,
            product(
                quotient,
                S,
                monomial("a", delay + 1),
                forward,
                monomial("a", 2),
            ),
        ),
        quotient.scale(
            -2,
            product(
                quotient,
                monomial("s", 3),
                monomial("a", delay + 1),
                forward,
            ),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                STAR,
                forward,
                monomial("a", delay + 1),
            ),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                forward,
                monomial("a", delay + 2),
            ),
        ),
        quotient.scale(
            -2,
            product(
                quotient,
                monomial("s", 2),
                backward,
                monomial("s", delay),
            ),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                monomial("s", delay + 1),
                backward,
                S,
            ),
        ),
        quotient.scale(
            -1,
            product(
                quotient,
                monomial("s", delay + 2),
                backward,
            ),
        ),
        quotient.scale(
            2,
            product(
                quotient,
                monomial("a", delay - 2),
                forward,
            ),
        ),
    )


def tail_layer_formula(
    quotient: DelayedQuotient,
    delay: int,
) -> Polynomial:
    """Return the eight new-crossing words in the odd associated layer."""

    outer = delay + 2
    return quotient.add(
        quotient.scale(
            2,
            {"a" * (delay - 1) + "s" * outer + "a": Fraction(1)},
        ),
        quotient.scale(
            -2,
            {"a" * (delay - 2) + "s" * outer + "aa": Fraction(1)},
        ),
        quotient.scale(
            2,
            {"s" + "a" * outer + "s" * (delay - 1): Fraction(1)},
        ),
        quotient.scale(
            -2,
            {"ss" + "a" * outer + "s" * (delay - 2): Fraction(1)},
        ),
        quotient.scale(
            -2,
            {"s" + "a" * outer + "s" * outer + "aaa": Fraction(1)},
        ),
        quotient.scale(
            -2,
            {"sss" + "a" * outer + "s" * outer + "a": Fraction(1)},
        ),
        quotient.scale(
            2,
            {"ss" + "a" * outer + "s" * outer + "aaaa": Fraction(1)},
        ),
        quotient.scale(
            2,
            {"ssss" + "a" * outer + "s" * outer + "aa": Fraction(1)},
        ),
    )


def renewal_coefficient(
    quotient: DelayedQuotient,
    partial: Polynomial,
    initial: Polynomial,
    final: Polynomial,
    identity: Polynomial,
    degree: int,
) -> Polynomial:
    """Return one unweighted closed-renewal coefficient."""

    pieces = renewal_pieces(
        quotient,
        partial,
        initial,
        final,
        identity,
        degree,
    )
    return quotient.add(
        pieces.direct[degree],
        pieces.immediate_return[degree],
        pieces.later_returns[degree],
    )


def physical_records(
    maximum_delay: int,
) -> list[OddFrontierRecord]:
    """Return exact physical audits through one finite delay."""

    records: list[OddFrontierRecord] = []
    for delay in range(1, maximum_delay + 1):
        quotient = DelayedQuotient(delay)
        degree = 2 * delay + 1
        quotient.maximum_degree = degree
        actual = renewal_coefficient(
            quotient,
            S,
            E,
            F,
            IDENTITY,
            degree,
        )
        expected = (
            first_frontier_formula(quotient)
            if delay == 1
            else odd_frontier_formula(
                quotient,
                delay,
                S,
                IDENTITY,
            )
        )
        formula_difference = quotient.add(
            actual,
            quotient.scale(-1, expected),
        )
        delayed_remainder = reduce_delay(expected, delay)

        tail_difference: Polynomial = {}
        layer_difference: Polynomial = {}
        if delay >= 3:
            tail_identity = quotient.add(
                IDENTITY,
                quotient.scale(-1, F),
            )
            tail_partial = quotient.multiply(S, tail_identity)
            tail_final = quotient.multiply(
                quotient.multiply(S, F),
                STAR,
            )
            tail_degree = degree - 2
            actual_tail = renewal_coefficient(
                quotient,
                tail_partial,
                E,
                tail_final,
                tail_identity,
                tail_degree,
            )
            expected_tail = odd_frontier_formula(
                quotient,
                delay - 1,
                tail_partial,
                tail_identity,
            )
            tail_difference = quotient.add(
                actual_tail,
                quotient.scale(-1, expected_tail),
            )
            actual_layer = quotient.add(
                actual,
                quotient.scale(-1, actual_tail),
            )
            expected_layer = tail_layer_formula(quotient, delay)
            layer_difference = quotient.add(
                actual_layer,
                quotient.scale(-1, expected_layer),
            )

        verified = bool(
            not formula_difference
            and not delayed_remainder
            and not tail_difference
            and not layer_difference
        )
        if not verified:
            raise RuntimeError(
                "odd-frontier physical audit failed: "
                f"delay={delay}, formula={formula_difference}, "
                f"delayed={delayed_remainder}, tail={tail_difference}, "
                f"layer={layer_difference}"
            )
        records.append(
            OddFrontierRecord(
                record_type="physical_frontier",
                delay=delay,
                face_degree=degree,
                frontier_word_count=len(actual),
                formula_difference_word_count=len(formula_difference),
                next_delay_remainder_word_count=len(delayed_remainder),
                tail_formula_difference_word_count=len(tail_difference),
                layer_formula_difference_word_count=len(layer_difference),
                all_checks_passed=verified,
            )
        )
    return records


def factorization_records(
    maximum_delay: int,
) -> list[OddFactorizationRecord]:
    """Return exact algebraic audits through one larger delay."""

    records: list[OddFactorizationRecord] = []
    for delay in range(1, maximum_delay + 1):
        quotient = DelayedQuotient(delay)
        quotient.maximum_degree = 2 * delay + 1
        frontier = (
            first_frontier_formula(quotient)
            if delay == 1
            else odd_frontier_formula(
                quotient,
                delay,
                S,
                IDENTITY,
            )
        )
        factorization_difference = quotient.add(
            frontier,
            quotient.scale(
                -1,
                odd_factorization(quotient, delay),
            ),
        )
        delayed_remainder = reduce_delay(frontier, delay)

        recursion_difference: Polynomial = {}
        if delay >= 3:
            tail_identity = quotient.add(
                IDENTITY,
                quotient.scale(-1, F),
            )
            tail_partial = quotient.multiply(S, tail_identity)
            embedded_tail = odd_frontier_formula(
                quotient,
                delay - 1,
                tail_partial,
                tail_identity,
            )
            recursion_difference = quotient.add(
                frontier,
                quotient.scale(-1, embedded_tail),
                quotient.scale(
                    -1,
                    tail_layer_formula(quotient, delay),
                ),
            )

        verified = bool(
            not factorization_difference
            and not delayed_remainder
            and not recursion_difference
        )
        if not verified:
            raise RuntimeError(
                "odd-frontier factorization audit failed: "
                f"delay={delay}, factor={factorization_difference}, "
                f"delayed={delayed_remainder}, "
                f"recursion={recursion_difference}"
            )
        records.append(
            OddFactorizationRecord(
                record_type="algebraic_factorization",
                delay=delay,
                frontier_word_count=len(frontier),
                factorization_difference_word_count=len(
                    factorization_difference
                ),
                next_delay_remainder_word_count=len(delayed_remainder),
                tail_recursion_difference_word_count=len(
                    recursion_difference
                ),
                all_checks_passed=verified,
            )
        )
    return records


def write_records(
    records: list[OddFrontierRecord | OddFactorizationRecord],
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
            "repeated_crabb_odd_frontier_delay_ideal_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact audits."""

    args = parse_args()
    records: list[OddFrontierRecord | OddFactorizationRecord] = [
        *physical_records(args.maximum_physical_delay),
        *factorization_records(args.maximum_factor_delay),
    ]
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
