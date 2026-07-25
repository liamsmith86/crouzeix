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

from crabb_palindromic_elliptic_hessian import direct_map_coefficients
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    S,
    Series,
    add,
    multiply,
    reduce_delay,
    reduce_delayed_word,
    scale,
)


@dataclass(frozen=True)
class CyclicRadialVolumeRecord:
    """One exact cyclic radial audit."""

    grade: int
    face_degree: int
    quotient_face_word_count: int
    cyclic_face_word_count: int
    stein_cyclic_word_count: int
    whitening_cyclic_word_count: int
    nonradial_cancellation_word_count: int
    radial_difference_word_count: int
    transfer_difference_word_count: int
    all_checks_passed: bool


class DelayedQuotient:
    """Exact word-series algebra with the active delay imposed early."""

    def __init__(self, grade: int) -> None:
        self.grade = grade
        self.maximum_degree = 2 * grade
        self.maximum_delay = grade - 1

    def reduce(self, polynomial: Polynomial) -> Polynomial:
        """Reduce a polynomial in the active two-sided delay ideal."""

        return reduce_delay(polynomial, self.maximum_delay)

    def add(self, *polynomials: Polynomial) -> Polynomial:
        """Add and immediately reduce quotient polynomials."""

        return self.reduce(add(*polynomials))

    def scale(
        self,
        coefficient: Fraction | int,
        polynomial: Polynomial,
    ) -> Polynomial:
        """Scale and immediately reduce a quotient polynomial."""

        return self.reduce(scale(coefficient, polynomial))

    def multiply(
        self,
        left: Polynomial,
        right: Polynomial,
    ) -> Polynomial:
        """Multiply in the active two-sided delay quotient."""

        return self.reduce(multiply(left, right))

    def series_multiply(self, left: Series, right: Series) -> Series:
        """Multiply two quotient series through the active degree."""

        result = [{} for _ in range(self.maximum_degree + 1)]
        for left_degree, left_polynomial in enumerate(left):
            if not left_polynomial:
                continue
            remaining = self.maximum_degree - left_degree
            for right_degree, right_polynomial in enumerate(right[: remaining + 1]):
                if not right_polynomial:
                    continue
                degree = left_degree + right_degree
                result[degree] = self.add(
                    result[degree],
                    self.multiply(left_polynomial, right_polynomial),
                )
        return result

    def series_adjoint(self, series: Series) -> Series:
        """Apply the formal adjoint coefficientwise."""

        return [self.polynomial_adjoint(polynomial) for polynomial in series]

    def series_power(
        self,
        series: Series,
        exponent: int,
        identity: Polynomial = IDENTITY,
    ) -> Series:
        """Raise a quotient series to a nonnegative integer power."""

        result = [
            identity,
            *({} for _ in range(self.maximum_degree)),
        ]
        base = series
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = self.series_multiply(result, base)
            remaining //= 2
            if remaining:
                base = self.series_multiply(base, base)
        return result

    def polynomial_adjoint(self, polynomial: Polynomial) -> Polynomial:
        """Apply the formal adjoint to one quotient polynomial."""

        translation = str.maketrans({"s": "a", "a": "s"})
        return self.reduce({
            word.translate(translation)[::-1]: coefficient
            for word, coefficient in polynomial.items()
        })

    def operator_series_for(
        self,
        partial: Polynomial,
        initial: Polynomial,
        final: Polynomial,
        identity: Polynomial,
    ) -> Series:
        """Return the balanced ellipse map for one quotient corner."""

        reverse = self.multiply(
            self.add(identity, final),
            self.multiply(
                self.polynomial_adjoint(partial),
                self.add(identity, initial),
            ),
        )
        pencil = [
            partial,
            reverse,
            *({} for _ in range(self.maximum_degree - 1)),
        ]
        result = [{} for _ in range(self.maximum_degree + 1)]
        coefficients = direct_map_coefficients(
            self.maximum_degree,
            self.maximum_degree + 1,
        )
        for index, coefficient_series in enumerate(coefficients):
            power = self.series_power(
                pencil,
                2 * index + 1,
                identity,
            )
            for scalar_degree in range(self.maximum_degree + 1):
                scalar = coefficient_series.coefficient(scalar_degree)
                if not scalar:
                    continue
                remaining = self.maximum_degree - scalar_degree
                for word_degree, polynomial in enumerate(power[: remaining + 1]):
                    if not polynomial:
                        continue
                    degree = scalar_degree + word_degree
                    result[degree] = self.add(
                        result[degree],
                        self.scale(scalar, polynomial),
                    )
        return result

    def operator_series(self) -> Series:
        """Return the balanced ellipse map in the delayed quotient."""

        return self.operator_series_for(S, E, F, IDENTITY)

    def polynomial_power(
        self,
        polynomial: Polynomial,
        exponent: int,
        identity: Polynomial,
    ) -> Polynomial:
        """Raise one quotient polynomial to a nonnegative power."""

        result = identity
        base = polynomial
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = self.multiply(result, base)
            remaining //= 2
            if remaining:
                base = self.multiply(base, base)
        return result

    def metric_series_for(
        self,
        partial: Polynomial,
        initial: Polynomial,
        final: Polynomial,
        identity: Polynomial,
    ) -> Series:
        """Return L219's metric for one quotient corner."""

        result = [
            identity,
            *({} for _ in range(self.maximum_degree)),
        ]
        partial_adjoint = self.polynomial_adjoint(partial)
        for order in range(2, self.maximum_degree + 1, 2):
            half_order = order // 2
            forward = self.polynomial_power(
                partial,
                half_order,
                identity,
            )
            backward = self.polynomial_power(
                partial_adjoint,
                half_order,
                identity,
            )
            coefficient = self.multiply(
                self.multiply(forward, final),
                backward,
            )
            for divisor in range(1, half_order + 1):
                if half_order % divisor:
                    continue
                forward = self.polynomial_power(
                    partial,
                    divisor,
                    identity,
                )
                backward = self.polynomial_power(
                    partial_adjoint,
                    divisor,
                    identity,
                )
                right_orbit = self.multiply(
                    self.multiply(backward, initial),
                    forward,
                )
                coefficient = self.add(
                    coefficient,
                    self.scale(
                        (-1) ** (half_order // divisor),
                        right_orbit,
                    ),
                )
            result[order] = coefficient
        return result

    def metric_series(self) -> Series:
        """Return L219's metric in the active delay quotient."""

        return self.metric_series_for(S, E, F, IDENTITY)

    def inverse_series(
        self,
        series: Series,
        identity: Polynomial = IDENTITY,
    ) -> Series:
        """Invert a quotient series with identity constant term."""

        result = [
            identity,
            *({} for _ in range(self.maximum_degree)),
        ]
        for degree in range(1, self.maximum_degree + 1):
            convolution: Polynomial = {}
            for positive_degree in range(1, degree + 1):
                convolution = self.add(
                    convolution,
                    self.multiply(
                        series[positive_degree],
                        result[degree - positive_degree],
                    ),
                )
            result[degree] = self.scale(-1, convolution)
        return result

    def mass_components_for(
        self,
        partial: Polynomial,
        initial: Polynomial,
        final: Polynomial,
        identity: Polynomial,
        active_degree: int,
    ) -> tuple[Polynomial, Polynomial, Polynomial]:
        """Return Stein, whitening, and total faces for one corner."""

        operator = self.operator_series_for(
            partial,
            initial,
            final,
            identity,
        )
        metric = self.metric_series_for(
            partial,
            initial,
            final,
            identity,
        )
        metric[active_degree] = {}

        final_weight = [
            final,
            *({} for _ in range(self.maximum_degree)),
        ]
        row_gram = self.series_multiply(
            self.series_multiply(
                self.series_adjoint(operator),
                final_weight,
            ),
            operator,
        )
        row_denominator = [
            self.add(left, self.scale(-1, right))
            for left, right in zip(metric, row_gram, strict=True)
        ]

        pulled_metric = self.series_multiply(
            self.series_multiply(
                self.series_adjoint(operator),
                metric,
            ),
            operator,
        )
        stein_slack = [
            self.add(left, self.scale(-1, right))
            for left, right in zip(metric, pulled_metric, strict=True)
        ]
        normalized_mass = self.series_multiply(
            self.inverse_series(row_denominator, identity),
            stein_slack,
        )
        stein_face = stein_slack[active_degree]
        total_face = normalized_mass[active_degree]
        whitening_face = self.add(
            total_face,
            self.scale(-1, stein_face),
        )
        return stein_face, whitening_face, total_face

    def mass_components(
        self,
    ) -> tuple[Polynomial, Polynomial, Polynomial]:
        """Return active Stein, whitening, and total mass faces."""

        return self.mass_components_for(
            S,
            E,
            F,
            IDENTITY,
            self.maximum_degree,
        )


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

    rotations = [word[position:] + word[:position] for position in range(len(word))]
    for rotation in rotations:
        rotated_reduction = reduce_delayed_word(
            rotation,
            maximum_delay,
        )
        if rotated_reduction == {rotation: Fraction(1)}:
            continue
        if not all(len(new_word) < len(word) for new_word in rotated_reduction):
            raise RuntimeError("cyclic reduction did not strictly shorten a word")
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


def audit_grade(grade: int) -> CyclicRadialVolumeRecord:
    """Audit one finite grade exactly."""

    quotient = DelayedQuotient(grade)
    stein_face, whitening_face, face = quotient.mass_components()
    maximum_delay = grade - 1
    cyclic_stein = cyclic_reduce(stein_face, maximum_delay)
    cyclic_whitening = cyclic_reduce(
        whitening_face,
        maximum_delay,
    )
    cyclic_face = cyclic_reduce(face, maximum_delay)
    radial = cyclic_reduce(radial_target(grade), maximum_delay)
    transfer = cyclic_reduce(
        transfer_target(grade),
        maximum_delay,
    )
    radial_words = set(radial)
    nonradial_stein = {
        word: coefficient
        for word, coefficient in cyclic_stein.items()
        if word not in radial_words
    }
    nonradial_whitening = {
        word: coefficient
        for word, coefficient in cyclic_whitening.items()
        if word not in radial_words
    }
    nonradial_cancellation = add(
        nonradial_stein,
        nonradial_whitening,
    )
    radial_difference = add(cyclic_face, scale(-1, radial))
    transfer_difference = add(cyclic_face, scale(-1, transfer))
    verified = bool(
        not nonradial_cancellation and not radial_difference and not transfer_difference
    )
    if not verified:
        raise RuntimeError(
            "cyclic radial volume audit failed: "
            f"grade={grade}, "
            f"nonradial={len(nonradial_cancellation)}, "
            f"radial={len(radial_difference)}, "
            f"transfer={len(transfer_difference)}"
        )
    return CyclicRadialVolumeRecord(
        grade=grade,
        face_degree=2 * grade,
        quotient_face_word_count=len(face),
        cyclic_face_word_count=len(cyclic_face),
        stein_cyclic_word_count=len(cyclic_stein),
        whitening_cyclic_word_count=len(cyclic_whitening),
        nonradial_cancellation_word_count=len(nonradial_cancellation),
        radial_difference_word_count=len(radial_difference),
        transfer_difference_word_count=len(transfer_difference),
        all_checks_passed=verified,
    )


def standard_records(
    maximum_grade: int,
) -> list[CyclicRadialVolumeRecord]:
    """Return exact audits through one maximum grade."""

    return [audit_grade(grade) for grade in range(1, maximum_grade + 1)]


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
    parser.add_argument("--maximum-grade", type=int, default=6)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiments/repeated_crabb_cyclic_radial_volume_s70224.jsonl"),
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
