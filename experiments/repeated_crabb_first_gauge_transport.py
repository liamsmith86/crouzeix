#!/usr/bin/env python3
"""Audit the first coupled transport of the normalized preparation gauge."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_canonical_quartic_preimage import (
    evaluate_polynomial,
)
from repeated_crabb_canonical_quintic_preimage import (
    prepared_endpoint_series,
)
from repeated_crabb_cyclic_radial_volume import DelayedQuotient
from repeated_crabb_delay_normalized_preparations import (
    SPECIFICATIONS,
    coefficient_l1,
    format_fraction,
    gauge_data,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    Q,
    S,
    STAR,
    Polynomial,
    add,
    adjoint,
    multiply,
    scale,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_gauge_transport_homology import (
    first_successor_transport,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
WitnessTerms = tuple[tuple[Fraction, str], ...]


@dataclass(frozen=True)
class ExactFirstGaugeTransportRecord:
    """Exact coupled quartic-to-quintic transport certificate."""

    record_type: str
    transport_witness_word_count: int
    transport_witness_coefficient_l1: str
    correction_word_count: int
    correction_coefficient_l1: str
    transported_quintic_word_count: int
    transported_quintic_coefficient_l1: str
    maximum_exact_residual_word_count: int
    correction_double_delay_quotient_word_count: int
    transported_double_delay_quotient_word_count: int
    all_checks_passed: bool


@dataclass(frozen=True)
class NumericalFirstGaugeTransportRecord:
    """One independent matrix audit of the coupled transport."""

    record_type: str
    state_dimension: int
    defect_dimension: int
    first_transfer_norm: str
    second_transfer_norm: str
    witness_stein_error: str
    witness_lower_endpoint_norm: str
    witness_upper_endpoint_norm: str
    full_endpoint_invariance_error: str
    transported_quintic_column_norm: str
    expected_delay_behavior: str
    all_checks_passed: bool


AuditRecord = (
    ExactFirstGaugeTransportRecord | NumericalFirstGaugeTransportRecord
)


TRANSPORT_WITNESS_TERMS: WitnessTerms = (
    (Fraction(-9), "aa"),
    (Fraction(-1), "aaaaaa"),
    (Fraction(-2), "assssa"),
    (Fraction(1), "aaaaaaas"),
    (Fraction(1), "aasssssa"),
    (Fraction(1), "asssssaa"),
    (Fraction(1), "saaaaaaa"),
    (Fraction(-1, 2), "aaassssssa"),
    (Fraction(1, 2), "aassssssaa"),
    (Fraction(-1), "asssaaaaas"),
    (Fraction(-1, 2), "assssssaaa"),
    (Fraction(-1, 2), "assssssssa"),
    (Fraction(-1), "saaaaasssa"),
    (Fraction(-1, 2), "aaaaaaasssaa"),
    (Fraction(1, 2), "aaasssssssaa"),
    (Fraction(-1, 2), "aasssaaaaaaa"),
    (Fraction(1, 2), "aasssssssaaa"),
    (Fraction(-1), "asssaaassssa"),
    (Fraction(-1), "assssaaasssa"),
    (Fraction(-3, 2), "aaasssaaaasssa"),
    (Fraction(-3, 2), "aasssaaaasssaa"),
    (Fraction(1, 2), "aasssaaasssssa"),
    (Fraction(2), "aassssaaassssa"),
    (Fraction(1, 2), "aasssssaaasssa"),
    (Fraction(1, 2), "asssaaaaaasssa"),
    (Fraction(-3, 2), "asssaaaasssaaa"),
    (Fraction(1, 2), "asssaaasssssaa"),
    (Fraction(2), "assssaaassssaa"),
    (Fraction(1, 2), "asssssaaasssaa"),
)


def transport_witness() -> Polynomial:
    """Return the Hermitian metric witness for the first transport."""

    witness: Polynomial = {}
    for coefficient, word in TRANSPORT_WITNESS_TERMS:
        monomial = {word: Fraction(1)}
        witness = add(
            witness,
            scale(coefficient, monomial),
            scale(coefficient, adjoint(monomial)),
        )
    return witness


def first_transport_data() -> tuple[
    Polynomial,
    Polynomial,
    Polynomial,
    Polynomial,
]:
    """Return raw forcing, witness, correction, and new quintic lift."""

    raw_transport = first_successor_transport()
    _, _, quintic_gauge, normalized_quintic = gauge_data(
        SPECIFICATIONS[1]
    )
    raw_forcing = add(
        raw_transport,
        scale(-1, quintic_gauge),
        scale(-1, adjoint(quintic_gauge)),
    )
    witness = transport_witness()
    coboundary = add(
        witness,
        scale(-1, multiply(multiply(STAR, witness), S)),
    )
    correction = multiply(
        multiply(
            Q,
            add(coboundary, scale(-1, raw_forcing)),
        ),
        E,
    )
    transported_quintic = add(normalized_quintic, correction)
    return raw_forcing, witness, correction, transported_quintic


def exact_record() -> ExactFirstGaugeTransportRecord:
    """Verify the complete rational word certificate."""

    raw, witness, correction, transported = first_transport_data()
    coboundary = add(
        witness,
        scale(-1, multiply(multiply(STAR, witness), S)),
    )
    quotient = DelayedQuotient(3)
    correction_quotient = quotient.reduce(correction)
    transported_quotient = quotient.reduce(transported)
    residuals = (
        add(witness, scale(-1, adjoint(witness))),
        multiply(multiply(E, witness), E),
        multiply(multiply(F, witness), F),
        add(
            raw,
            correction,
            adjoint(correction),
            scale(-1, coboundary),
        ),
        multiply(E, correction),
        correction_quotient,
        transported_quotient,
    )
    maximum_residual = max(map(len, residuals))
    verified = maximum_residual == 0
    if not verified:
        raise RuntimeError(
            "the exact first gauge-transport certificate failed: "
            f"residuals={tuple(map(len, residuals))}"
        )
    return ExactFirstGaugeTransportRecord(
        record_type="exact_coupled_transport",
        transport_witness_word_count=len(witness),
        transport_witness_coefficient_l1=format_fraction(
            coefficient_l1(witness)
        ),
        correction_word_count=len(correction),
        correction_coefficient_l1=format_fraction(
            coefficient_l1(correction)
        ),
        transported_quintic_word_count=len(transported),
        transported_quintic_coefficient_l1=format_fraction(
            coefficient_l1(transported)
        ),
        maximum_exact_residual_word_count=maximum_residual,
        correction_double_delay_quotient_word_count=len(
            correction_quotient
        ),
        transported_double_delay_quotient_word_count=len(
            transported_quotient
        ),
        all_checks_passed=verified,
    )


def numerical_record(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> NumericalFirstGaugeTransportRecord:
    """Evaluate the transported series and both full endpoints."""

    raw, witness, correction, transported = first_transport_data()
    _, _, quartic_gauge, _ = gauge_data(SPECIFICATIONS[0])
    _, _, quintic_gauge, _ = gauge_data(SPECIFICATIONS[1])

    raw_matrix = evaluate_polynomial(raw, operator)
    correction_matrix = evaluate_polynomial(correction, operator)
    witness_matrix = evaluate_polynomial(witness, operator)
    reconstructed = stein_inverse(
        operator,
        raw_matrix + correction_matrix + correction_matrix.conj().T,
    )
    witness_error = float(np.linalg.norm(reconstructed - witness_matrix))
    lower_witness = float(
        np.linalg.norm(right.conj().T @ reconstructed @ right)
    )
    upper_witness = float(
        np.linalg.norm(left.conj().T @ reconstructed @ left)
    )

    quartic_change = (
        -evaluate_polynomial(quartic_gauge, operator) @ right
    )
    quintic_change = (
        evaluate_polynomial(
            add(scale(-1, quintic_gauge), correction),
            operator,
        )
        @ right
    )
    original_upper, original_lower, _ = prepared_endpoint_series(
        operator,
        right,
        left,
    )
    changed_upper, changed_lower, _ = prepared_endpoint_series(
        operator,
        right,
        left,
        additional_factor_columns={
            4: quartic_change,
            5: quintic_change,
        },
    )
    endpoint_error = max(
        (
            float(
                np.linalg.norm(changed_upper[degree] - original_upper[degree])
                + np.linalg.norm(
                    changed_lower[degree] - original_lower[degree]
                )
            )
            for degree in range(6)
        ),
        default=0.0,
    )

    transported_column = (
        evaluate_polynomial(transported, operator) @ right
    )
    transported_norm = float(np.linalg.norm(transported_column))
    first_transfer = float(
        np.linalg.norm(
            transfer_coefficient(operator, right, left, 1)
        )
    )
    second_transfer = float(
        np.linalg.norm(
            transfer_coefficient(operator, right, left, 2)
        )
    )

    tolerance = 4e-8
    common_checks = bool(
        witness_error < tolerance
        and lower_witness < tolerance
        and upper_witness < tolerance
        and endpoint_error < tolerance
    )
    if construction_kind == "double_complete_delay":
        behavior = "first_two_transfers_and_new_quintic_zero"
        behavior_check = bool(
            first_transfer < tolerance
            and second_transfer < tolerance
            and transported_norm < tolerance
        )
    else:
        behavior = "global_endpoint_invariance"
        behavior_check = True

    verified = common_checks and behavior_check
    if not verified:
        raise RuntimeError(
            "the numerical first gauge-transport audit failed: "
            f"kind={construction_kind}, "
            f"witness={witness_error:.3e}, "
            f"lower={lower_witness:.3e}, "
            f"upper={upper_witness:.3e}, "
            f"endpoint={endpoint_error:.3e}, "
            f"column={transported_norm:.3e}"
        )
    return NumericalFirstGaugeTransportRecord(
        record_type=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        first_transfer_norm=format_float(first_transfer),
        second_transfer_norm=format_float(second_transfer),
        witness_stein_error=format_float(witness_error),
        witness_lower_endpoint_norm=format_float(lower_witness),
        witness_upper_endpoint_norm=format_float(upper_witness),
        full_endpoint_invariance_error=format_float(endpoint_error),
        transported_quintic_column_norm=format_float(transported_norm),
        expected_delay_behavior=behavior,
        all_checks_passed=verified,
    )


def standard_records() -> list[AuditRecord]:
    """Return exact, unstructured, flag, and delayed audits."""

    records: list[AuditRecord] = [exact_record()]
    operator, right, left = random_partial_isometry(
        12,
        2,
        np.random.default_rng(702_261),
    )
    records.append(
        numerical_record("unstructured", operator, right, left)
    )

    operator, right, left, _, _ = rank_chain_case(
        4,
        702_262,
        0.8,
    )
    records.append(
        numerical_record("rank_changing_flag", operator, right, left)
    )

    for grade, construction_kind, seed in (
        (2, "single_complete_delay", 702_263),
        (3, "double_complete_delay", 702_264),
    ):
        operator, right, left, _ = inflated_case(
            9,
            2,
            grade,
            2,
            seed,
        )
        records.append(
            numerical_record(
                construction_kind,
                operator,
                right,
                left,
            )
        )
    return records


def write_records(records: list[AuditRecord], output: Path) -> str:
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
            "repeated_crabb_first_gauge_transport_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
