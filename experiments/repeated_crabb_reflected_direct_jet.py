#!/usr/bin/env python3
"""Audit the universal first-reflection direct-map jet.

After a complete delay of length ``r``, L237 expresses the difference
between the finite terminal self-energy and the half-line self-energy as

    gamma_r - eta = c^r z^(-(2r-1)) + O(c^(r+1)).

This script performs the resulting Laurent-resolvent coefficient
calculation in an exact free-word algebra.  It checks that, for every
tested ``r >= 2``, the first three reflected coefficients are the same
three word polynomials (up to the common sign ``(-1)^r``).

The audit deliberately does not impose partial-isometry relations.  It
therefore checks the scalar Riemann coefficients, Laurent residues, and
noncommutative multiplication order before any state-algebra collapse.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp

from crabb_palindromic_elliptic_hessian import direct_map_coefficients


Polynomial = dict[str, Fraction]


@dataclass(frozen=True)
class ReflectedDirectJetRecord:
    """One exact delay audit."""

    delay: int
    leading_degree: int
    leading_word_count: int
    odd_word_count: int
    next_even_word_count: int
    self_energy_leading_coefficient: str
    self_energy_next_coefficient: str
    symbolic_all_delay_check_passed: bool
    all_checks_passed: bool


@dataclass(frozen=True)
class Letter:
    """One letter in the half-line resolvent expansion."""

    word: str
    c_degree: int
    extra_z_degree: int
    coefficient: Fraction


S_LETTER = Letter("s", 0, 0, Fraction(1))
J_LETTER = Letter("j", 1, 0, Fraction(1))


def clean(polynomial: Polynomial) -> Polynomial:
    """Remove zero coefficients."""

    return {
        word: coefficient
        for word, coefficient in polynomial.items()
        if coefficient
    }


def add(*polynomials: Polynomial) -> Polynomial:
    """Add sparse free-word polynomials."""

    result: defaultdict[str, Fraction] = defaultdict(Fraction)
    for polynomial in polynomials:
        for word, coefficient in polynomial.items():
            result[word] += coefficient
    return clean(dict(result))


def scale(coefficient: Fraction | int, polynomial: Polynomial) -> Polynomial:
    """Scale a free-word polynomial."""

    coefficient = Fraction(coefficient)
    return clean({
        word: coefficient * value
        for word, value in polynomial.items()
    })


def monomial(word: str, coefficient: Fraction | int = 1) -> Polynomial:
    """Return one free-word monomial."""

    coefficient = Fraction(coefficient)
    return {word: coefficient} if coefficient else {}


def series_inverse_one_minus_t(
    coefficients: list[Fraction],
    order: int,
) -> list[Fraction]:
    """Return ``1 / (1 - t*h(t))`` through one order."""

    result = [Fraction(0) for _ in range(order)]
    result[0] = Fraction(1)
    for degree in range(1, order):
        result[degree] = sum(
            coefficients[index] * result[degree - 1 - index]
            for index in range(min(degree, len(coefficients)))
        )
    return result


def catalan_coefficients(order: int) -> list[Fraction]:
    """Return the Catalan fixed-point series through one order."""

    result = [Fraction(0) for _ in range(order)]
    result[0] = Fraction(1)
    for degree in range(1, order):
        result[degree] = sum(
            result[index] * result[degree - 1 - index]
            for index in range(degree)
        )
    return result


def terminal_error_coefficients(
    delay: int,
    extra_order: int,
) -> list[Fraction]:
    """Return coefficients of ``t^(-(r-1))(h_r-C)``."""

    order = delay + extra_order
    h = [Fraction(2), *(
        Fraction(0) for _ in range(order - 1)
    )]
    for _ in range(2, delay + 1):
        h = series_inverse_one_minus_t(h, order)
    catalan = catalan_coefficients(order)
    difference = [
        left - right
        for left, right in zip(h, catalan, strict=True)
    ]
    return difference[delay - 1 : delay - 1 + extra_order]


def weak_compositions(total: int, parts: int) -> list[tuple[int, ...]]:
    """Return weak compositions of ``total`` into ``parts`` slots."""

    if parts == 1:
        return [(total,)]
    result: list[tuple[int, ...]] = []
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            result.append((first, *tail))
    return result


def expected_jet() -> tuple[Polynomial, Polynomial, Polynomial]:
    """Return the three universal unreduced word coefficients."""

    leading = add(*(monomial("s" * left + "f" + "s" * right)
                    for left, right in weak_compositions(1, 2)))

    odd = add(
        monomial("jf"),
        monomial("fj"),
        *(scale(-1, monomial("s" * left + "f" + "s" * right))
          for left, right in weak_compositions(3, 2)),
    )

    next_even = scale(3, leading)
    next_even = add(
        next_even,
        *(monomial("s" * left + "f" + "s" * right)
          for left, right in weak_compositions(5, 2)),
    )
    next_even = add(
        next_even,
        *(scale(
            -1,
            monomial(
                "s" * left
                + "j"
                + "s" * middle
                + "f"
                + "s" * right
            ),
        ) for left, middle, right in weak_compositions(2, 3)),
        *(scale(
            -1,
            monomial(
                "s" * left
                + "f"
                + "s" * middle
                + "j"
                + "s" * right
            ),
        ) for left, middle, right in weak_compositions(2, 3)),
        *(scale(
            -2,
            monomial(
                "s" * left
                + "f"
                + "s" * middle
                + "f"
                + "s" * right
            ),
        ) for left, middle, right in weak_compositions(1, 3)),
    )
    return leading, odd, next_even


def symbolic_stable_jet() -> tuple[
    dict[str, sp.Expr],
    dict[str, sp.Expr],
    dict[str, sp.Expr],
]:
    """Enumerate the stable ``r >= 3`` jet with symbolic delay."""

    delay = sp.symbols("r", integer=True, positive=True)
    letters = (
        Letter("s", 0, 0, Fraction(1)),
        Letter("j", 1, 0, Fraction(1)),
        Letter("f", 1, -1, Fraction(1)),
        Letter("f", 2, -3, Fraction(1)),
        Letter("f", 3, -5, Fraction(2)),
    )
    terminal_coefficients = (
        sp.Integer(1),
        2 * delay - 2,
        2 * delay**2 - delay - 2,
    )
    result: list[defaultdict[str, sp.Expr]] = [
        defaultdict(lambda: sp.Integer(0))
        for _ in range(3)
    ]

    for length in range(6):
        for sequence in itertools.product(letters, repeat=length):
            sequence_c = sum(letter.c_degree for letter in sequence)
            extra_z = sum(letter.extra_z_degree for letter in sequence)
            sequence_scalar = sp.prod(
                int(letter.coefficient)
                for letter in sequence
            )
            for split in range(length + 1):
                word = (
                    "".join(letter.word for letter in sequence[:split])
                    + "f"
                    + "".join(letter.word for letter in sequence[split:])
                )
                for error_index, error_coefficient in enumerate(
                    terminal_coefficients
                ):
                    numerator = (
                        length
                        - extra_z
                        + 2 * error_index
                        - 1
                    )
                    if numerator % 2:
                        continue
                    scalar_index_offset = numerator // 2
                    scalar_index = delay + scalar_index_offset
                    edge_sign = (-1) ** scalar_index_offset
                    for extra_scalar_degree in (0, 2):
                        relative_degree = (
                            sequence_c
                            + error_index
                            + scalar_index_offset
                            + extra_scalar_degree
                        )
                        if not 0 <= relative_degree <= 2:
                            continue
                        scalar_coefficient = edge_sign
                        if extra_scalar_degree == 2:
                            scalar_coefficient *= 2 * scalar_index + 1
                        result[relative_degree][word] += (
                            sequence_scalar
                            * error_coefficient
                            * scalar_coefficient
                        )
    return tuple(
        {
            word: sp.expand(coefficient)
            for word, coefficient in polynomial.items()
            if sp.expand(coefficient) != 0
        }
        for polynomial in result
    )  # type: ignore[return-value]


def resolvent_sequences(maximum_length: int) -> list[tuple[Letter, ...]]:
    """Enumerate the only half-line words relevant through relative order two."""

    eta_letters = [
        Letter("f", index + 1, -(2 * index + 1), coefficient)
        for index, coefficient in enumerate(catalan_coefficients(3))
    ]
    alphabet = [S_LETTER, J_LETTER, *eta_letters]
    return [
        sequence
        for length in range(maximum_length + 1)
        for sequence in itertools.product(alphabet, repeat=length)
    ]


def reflected_coefficients(
    delay: int,
) -> tuple[Polynomial, Polynomial, Polynomial]:
    """Extract the linear first-reflection response through relative order two."""

    maximum_degree = 2 * delay + 2
    scalar_coefficients = direct_map_coefficients(
        delay + 2,
        maximum_degree + 1,
    )
    terminal_error = terminal_error_coefficients(delay, 3)
    sequences = resolvent_sequences(5)
    result = [
        defaultdict(Fraction)
        for _ in range(3)
    ]

    for sequence in sequences:
        sequence_c = sum(letter.c_degree for letter in sequence)
        sequence_z = -(len(sequence) + 2) + sum(
            letter.extra_z_degree for letter in sequence
        )
        sequence_scalar = Fraction(1)
        for letter in sequence:
            sequence_scalar *= letter.coefficient
        for split in range(len(sequence) + 1):
            left = sequence[:split]
            right = sequence[split:]
            word = (
                "".join(letter.word for letter in left)
                + "f"
                + "".join(letter.word for letter in right)
            )
            for error_index, error_coefficient in enumerate(
                terminal_error
            ):
                if not error_coefficient:
                    continue
                error_c = delay + error_index
                error_z = -(2 * delay + 2 * error_index - 1)
                for scalar_index, scalar_series in enumerate(
                    scalar_coefficients
                ):
                    scalar_z = 2 * scalar_index + 1
                    for scalar_c in range(maximum_degree + 1):
                        scalar_coefficient = scalar_series.coefficient(
                            scalar_c
                        )
                        if not scalar_coefficient:
                            continue
                        total_c = (
                            sequence_c
                            + error_c
                            + scalar_c
                        )
                        relative_degree = total_c - 2 * delay
                        if not 0 <= relative_degree <= 2:
                            continue
                        total_z = (
                            sequence_z + error_z + scalar_z
                        )
                        if total_z != -1:
                            continue
                        result[relative_degree][word] += (
                            sequence_scalar
                            * error_coefficient
                            * scalar_coefficient
                        )
    return tuple(clean(dict(item)) for item in result)  # type: ignore[return-value]


def exact_records(maximum_delay: int) -> list[ReflectedDirectJetRecord]:
    """Run the exact audit for delays two through one maximum."""

    expected = expected_jet()
    symbolic = symbolic_stable_jet()
    symbolic_verified = all(
        set(actual) == set(target)
        and all(
            sp.expand(actual[word] - target[word]) == 0
            for word in target
        )
        for actual, target in zip(symbolic, expected, strict=True)
    )
    if not symbolic_verified:
        raise RuntimeError("the symbolic all-delay jet audit failed")

    records: list[ReflectedDirectJetRecord] = []
    for delay in range(2, maximum_delay + 1):
        actual = reflected_coefficients(delay)
        sign = (-1) ** delay
        checks = [
            coefficient == scale(sign, target)
            for coefficient, target in zip(actual, expected, strict=True)
        ]
        terminal_error = terminal_error_coefficients(delay, 3)
        checks.extend([
            terminal_error[0] == 1,
            terminal_error[1] == 2 * delay - 2,
        ])
        verified = all(checks)
        if not verified:
            raise RuntimeError(
                f"reflected direct jet failed at delay {delay}: {checks}"
            )
        records.append(
            ReflectedDirectJetRecord(
                delay=delay,
                leading_degree=2 * delay,
                leading_word_count=len(actual[0]),
                odd_word_count=len(actual[1]),
                next_even_word_count=len(actual[2]),
                self_energy_leading_coefficient=str(terminal_error[0]),
                self_energy_next_coefficient=str(terminal_error[1]),
                symbolic_all_delay_check_passed=symbolic_verified,
                all_checks_passed=verified,
            )
        )
    return records


def write_records(
    records: list[ReflectedDirectJetRecord],
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
    parser.add_argument("--maximum-delay", type=int, default=5)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_reflected_direct_jet_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact reflected-jet audit."""

    args = parse_args()
    records = exact_records(args.maximum_delay)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
