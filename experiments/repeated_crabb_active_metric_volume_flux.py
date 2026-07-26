#!/usr/bin/env python3
"""Audit the active metric-deletion response and the exact volume flux."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from repeated_crabb_cyclic_radial_volume import (
    DelayedQuotient,
    cyclic_reduce,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    S,
)
from repeated_crabb_even_frontier_delay_ideal import (
    even_frontier_formula,
)
from repeated_crabb_odd_frontier_delay_ideal import product


@dataclass(frozen=True)
class ActivePhysicalRecord:
    """One exact physical active-face audit."""

    record_type: str
    grade: int
    face_degree: int
    maximum_lower_word_count: int
    intact_word_count: int
    response_word_count: int
    edge_deleted_word_count: int
    response_difference_word_count: int
    cyclic_target_difference_word_count: int
    all_checks_passed: bool


@dataclass(frozen=True)
class ActiveAlgebraicRecord:
    """One exact algebraic active-flux audit."""

    record_type: str
    grade: int
    metric_word_count: int
    intact_cyclic_difference_word_count: int
    response_cyclic_difference_word_count: int
    total_cyclic_difference_word_count: int
    all_checks_passed: bool


def orbit(
    quotient: DelayedQuotient,
    left: str,
    middle: Polynomial,
    right: str,
) -> Polynomial:
    """Return one defect orbit."""

    return product(
        quotient,
        {left: Fraction(1)},
        middle,
        {right: Fraction(1)},
    )


def metric_coefficient(
    quotient: DelayedQuotient,
    grade: int,
) -> Polynomial:
    """Return L219's exact coefficient ``X_grade``."""

    terms = [
        orbit(
            quotient,
            "s" * grade,
            F,
            "a" * grade,
        )
    ]
    for divisor in range(1, grade + 1):
        if grade % divisor:
            continue
        terms.append(
            quotient.scale(
                (-1) ** (grade // divisor),
                orbit(
                    quotient,
                    "a" * divisor,
                    E,
                    "s" * divisor,
                ),
            )
        )
    return quotient.add(*terms)


def retained_projection(
    quotient: DelayedQuotient,
) -> Polynomial:
    """Return the retained output projection."""

    return quotient.add(
        IDENTITY,
        quotient.scale(-1, F),
    )


def deletion_response(
    quotient: DelayedQuotient,
    metric: Polynomial,
) -> Polynomial:
    """Return ``P(X-SXS*)P`` for deleting the active metric face."""

    shifted = product(
        quotient,
        S,
        metric,
        quotient.polynomial_adjoint(S),
    )
    retained = retained_projection(quotient)
    return product(
        quotient,
        retained,
        quotient.add(metric, quotient.scale(-1, shifted)),
        retained,
    )


def radial_flux_target(
    quotient: DelayedQuotient,
    grade: int,
    multiplier: int,
) -> Polynomial:
    """Return a multiple of the three-term radial flux."""

    q_one = product(
        quotient,
        {"a": Fraction(1)},
        {"s": Fraction(1)},
    )
    q_top = product(
        quotient,
        {"a" * (grade + 2): Fraction(1)},
        {"s" * (grade + 2): Fraction(1)},
    )
    return quotient.add(
        quotient.scale(multiplier, q_top),
        quotient.scale(-multiplier * (grade + 2), q_one),
        quotient.scale(multiplier * (grade + 1), IDENTITY),
    )


def physical_records(
    maximum_grade: int,
) -> list[ActivePhysicalRecord]:
    """Return exact physical audits through one finite grade."""

    records: list[ActivePhysicalRecord] = []
    for grade in range(1, maximum_grade + 1):
        quotient = DelayedQuotient(grade)
        face_degree = 2 * grade
        cutoff = face_degree + 2
        quotient.maximum_degree = cutoff
        edge_series = quotient.closed_return_defect_for(
            S,
            E,
            F,
            IDENTITY,
            face_degree,
        )
        intact_series = quotient.closed_return_defect_for(
            S,
            E,
            F,
            IDENTITY,
            cutoff,
        )
        edge = edge_series[face_degree]
        intact = intact_series[face_degree]
        metric = metric_coefficient(quotient, grade)
        response = deletion_response(quotient, metric)
        response_difference = quotient.add(
            edge,
            quotient.scale(-1, intact),
            quotient.scale(-1, response),
        )
        cyclic_difference = quotient.add(
            cyclic_reduce(edge, grade - 1),
            quotient.scale(
                -1,
                radial_flux_target(quotient, grade, 4),
            ),
        )
        maximum_lower = max(
            (
                len(edge_series[degree])
                for degree in range(face_degree)
            ),
            default=0,
        )
        verified = bool(
            maximum_lower == 0
            and not response_difference
            and not cyclic_difference
        )
        if not verified:
            raise RuntimeError(
                "active metric response audit failed: "
                f"grade={grade}, lower={maximum_lower}, "
                f"response={response_difference}, "
                f"cyclic={cyclic_difference}"
            )
        records.append(
            ActivePhysicalRecord(
                record_type="physical_active_face",
                grade=grade,
                face_degree=face_degree,
                maximum_lower_word_count=maximum_lower,
                intact_word_count=len(intact),
                response_word_count=len(response),
                edge_deleted_word_count=len(edge),
                response_difference_word_count=len(
                    response_difference
                ),
                cyclic_target_difference_word_count=len(
                    cyclic_difference
                ),
                all_checks_passed=verified,
            )
        )
    return records


def algebraic_records(
    maximum_grade: int,
) -> list[ActiveAlgebraicRecord]:
    """Return exact cyclic audits through one larger grade."""

    records: list[ActiveAlgebraicRecord] = []
    for grade in range(1, maximum_grade + 1):
        quotient = DelayedQuotient(grade)
        quotient.maximum_degree = 2 * grade + 2
        intact = even_frontier_formula(
            quotient,
            grade,
            S,
            IDENTITY,
        )
        metric = metric_coefficient(quotient, grade)
        response = deletion_response(quotient, metric)
        half_target = radial_flux_target(quotient, grade, 2)
        full_target = radial_flux_target(quotient, grade, 4)
        intact_difference = quotient.add(
            cyclic_reduce(intact, grade - 1),
            quotient.scale(-1, half_target),
        )
        response_difference = quotient.add(
            cyclic_reduce(response, grade - 1),
            quotient.scale(-1, half_target),
        )
        total_difference = quotient.add(
            cyclic_reduce(
                quotient.add(intact, response),
                grade - 1,
            ),
            quotient.scale(-1, full_target),
        )
        verified = bool(
            not intact_difference
            and not response_difference
            and not total_difference
        )
        if not verified:
            raise RuntimeError(
                "active algebraic flux audit failed: "
                f"grade={grade}, intact={intact_difference}, "
                f"response={response_difference}, "
                f"total={total_difference}"
            )
        records.append(
            ActiveAlgebraicRecord(
                record_type="algebraic_cyclic_flux",
                grade=grade,
                metric_word_count=len(metric),
                intact_cyclic_difference_word_count=len(
                    intact_difference
                ),
                response_cyclic_difference_word_count=len(
                    response_difference
                ),
                total_cyclic_difference_word_count=len(
                    total_difference
                ),
                all_checks_passed=verified,
            )
        )
    return records


def write_records(
    records: list[ActivePhysicalRecord | ActiveAlgebraicRecord],
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
    parser.add_argument("--maximum-physical-grade", type=int, default=4)
    parser.add_argument("--maximum-algebraic-grade", type=int, default=12)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_active_metric_volume_flux_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact audits."""

    args = parse_args()
    records: list[ActivePhysicalRecord | ActiveAlgebraicRecord] = [
        *physical_records(args.maximum_physical_grade),
        *algebraic_records(args.maximum_algebraic_grade),
    ]
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
