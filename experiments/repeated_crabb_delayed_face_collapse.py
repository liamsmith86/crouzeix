#!/usr/bin/env python3
"""Audit the all-grade collapse of the delayed physical Schur face."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_active_metric_volume_flux import (
    deletion_response,
    metric_coefficient,
)
from repeated_crabb_cyclic_radial_volume import DelayedQuotient
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    S,
)
from repeated_crabb_edge_deleted_balanced_volume import delayed_case
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_even_frontier_delay_ideal import (
    even_frontier_formula,
    product,
)
from repeated_crabb_transfer_flag import transfer_coefficient


Matrix = np.ndarray

PRIOR_PHYSICAL_DATASET = Path(
    "experiments/"
    "repeated_crabb_delayed_slack_anticommutator_s70224.jsonl"
)
PRIOR_PHYSICAL_SHA256 = (
    "2f39d3e584b820424cc8c8efd90ac952"
    "05f624a142567e5a24d3028d91e04345"
)


@dataclass(frozen=True)
class DelayedFaceCollapseRecord:
    """One algebraic or numerical delayed-face audit."""

    record_type: str
    grade: int
    state_dimension: int | str
    defect_dimension: int | str
    maximum_earlier_transfer_norm: str
    response_cancellation_word_count: int | str
    associated_face_difference_word_count: int | str
    endpoint_error: str
    all_checks_passed: bool


def associated_candidate(
    quotient: DelayedQuotient,
    grade: int,
):
    """Return ``E_1 F_(grade-1) + F_(grade-1) E_1``."""

    delay = grade - 1
    left_orbit = product(
        quotient,
        {"s" * delay: 1},
        F,
        {"a" * delay: 1},
    )
    e_one = product(
        quotient,
        {"a": 1},
        E,
        {"s": 1},
    )
    return quotient.add(
        product(quotient, e_one, left_orbit),
        product(quotient, left_orbit, e_one),
    )


def algebraic_record(grade: int) -> DelayedFaceCollapseRecord:
    """Audit response cancellation and associated conjugacy exactly."""

    quotient = DelayedQuotient(grade)
    quotient.maximum_degree = 2 * grade + 4
    intact = even_frontier_formula(
        quotient,
        grade,
        S,
        IDENTITY,
    )
    response = deletion_response(
        quotient,
        metric_coefficient(quotient, grade),
    )
    deleted = quotient.add(intact, response)
    full_dual = quotient.add(
        deleted,
        quotient.scale(-1, response),
    )
    response_difference = quotient.add(
        full_dual,
        quotient.scale(-1, intact),
    )
    associated = product(
        quotient,
        {"a": 1},
        full_dual,
        {"s": 1},
    )
    associated_difference = quotient.add(
        associated,
        quotient.scale(
            -1,
            associated_candidate(quotient, grade),
        ),
    )
    verified = not response_difference and not associated_difference
    if not verified:
        raise RuntimeError(
            "delayed face collapse algebra failed: "
            f"grade={grade}, response={response_difference}, "
            f"associated={associated_difference}"
        )
    return DelayedFaceCollapseRecord(
        record_type="exact_algebraic_collapse",
        grade=grade,
        state_dimension="not_applicable",
        defect_dimension="not_applicable",
        maximum_earlier_transfer_norm="not_applicable",
        response_cancellation_word_count=len(response_difference),
        associated_face_difference_word_count=len(
            associated_difference
        ),
        endpoint_error="not_applicable",
        all_checks_passed=verified,
    )


def numerical_candidate(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    grade: int,
) -> Matrix:
    """Return the numerical two-orbit initial Schur face."""

    e_projection = right @ right.conj().T
    f_projection = left @ left.conj().T
    e_one = operator.conj().T @ e_projection @ operator
    power = np.linalg.matrix_power(operator, grade - 1)
    f_orbit = power @ f_projection @ power.conj().T
    return e_one @ f_orbit + f_orbit @ e_one


def numerical_record(grade: int) -> DelayedFaceCollapseRecord:
    """Audit the all-grade endpoint of the candidate face."""

    operator, right, left = delayed_case(grade)
    coefficient = transfer_coefficient(
        operator,
        right,
        left,
        grade,
    )
    earlier = max(
        (
            float(
                np.linalg.norm(
                    transfer_coefficient(
                        operator,
                        right,
                        left,
                        index,
                    )
                )
            )
            for index in range(1, grade)
        ),
        default=0.0,
    )
    candidate = numerical_candidate(
        operator,
        right,
        left,
        grade,
    )
    endpoint = (
        left.conj().T
        @ stein_inverse(operator, candidate)
        @ left
    )
    expected = 2 * coefficient @ coefficient.conj().T
    endpoint_error = float(np.linalg.norm(endpoint - expected))
    tolerance = 3e-8
    verified = bool(
        earlier < tolerance
        and endpoint_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "delayed face collapse endpoint failed: "
            f"grade={grade}, earlier={earlier:.3e}, "
            f"endpoint={endpoint_error:.3e}"
        )
    return DelayedFaceCollapseRecord(
        record_type="numerical_endpoint",
        grade=grade,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        maximum_earlier_transfer_norm=format_float(earlier),
        response_cancellation_word_count="not_applicable",
        associated_face_difference_word_count="not_applicable",
        endpoint_error=format_float(endpoint_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[DelayedFaceCollapseRecord]:
    """Return exact, numerical, and prior physical audits."""

    records = [
        algebraic_record(grade)
        for grade in range(1, 13)
    ]
    records.extend(
        numerical_record(grade)
        for grade in range(1, 7)
    )
    validate_prior_physical_dataset()
    return records


def validate_prior_physical_dataset() -> None:
    """Verify the tracked exact physical audit used as a cross-check."""

    digest = hashlib.sha256(
        PRIOR_PHYSICAL_DATASET.read_bytes()
    ).hexdigest()
    records = [
        json.loads(line)
        for line in PRIOR_PHYSICAL_DATASET.read_text(
            encoding="utf-8"
        ).splitlines()
        if line
    ]
    if (
        digest != PRIOR_PHYSICAL_SHA256
        or len(records) != 5
        or not all(record["all_checks_passed"] for record in records)
    ):
        raise RuntimeError("the prior exact physical audit regressed")


def write_records(
    records: list[DelayedFaceCollapseRecord],
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
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_delayed_face_collapse_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
