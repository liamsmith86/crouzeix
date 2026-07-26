#!/usr/bin/env python3
"""Audit the second coupled transport of the normalized preparation gauge."""

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
from repeated_crabb_canonical_quartic_preimage import evaluate_polynomial
from repeated_crabb_canonical_quintic_preimage import (
    joint_left_kernel,
    prepared_endpoint_series,
)
from repeated_crabb_canonical_sextic_preimage import (
    prepared_exact_metrics,
    sextic_column,
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
    IDENTITY,
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
from repeated_crabb_first_gauge_transport import first_transport_data
from repeated_crabb_sextic_transport_trace_obstruction import (
    transported_sextic_forcing,
)
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
WitnessTerms = tuple[tuple[Fraction, str], ...]


@dataclass(frozen=True)
class ExactSecondGaugeTransportRecord:
    """Exact rational certificate for the coupled sextic transport."""

    record_type: str
    transport_witness_word_count: int
    transport_witness_coefficient_l1: str
    perpendicular_correction_word_count: int
    perpendicular_correction_coefficient_l1: str
    parallel_correction_word_count: int
    parallel_correction_coefficient_l1: str
    transported_sextic_word_count: int
    transported_sextic_coefficient_l1: str
    endpoint_response_word_count: int
    endpoint_response_coefficient_l1: str
    endpoint_factor_word_count: int
    endpoint_factor_coefficient_l1: str
    maximum_exact_residual_word_count: int
    correction_triple_delay_quotient_word_count: int
    transported_triple_delay_quotient_word_count: int
    endpoint_first_delay_quotient_word_count: int
    all_checks_passed: bool


@dataclass(frozen=True)
class NumericalSecondGaugeTransportRecord:
    """One independent matrix audit of the sextic transport."""

    record_type: str
    state_dimension: int
    defect_dimension: int
    first_transfer_norm: str
    second_transfer_norm: str
    surviving_flag_dimension: int
    witness_stein_error: str
    witness_lower_endpoint_norm: str
    witness_upper_endpoint_norm: str
    witness_flag_endpoint_norm: str
    nonlinear_lower_sixth_change_norm: str
    lower_order_upper_invariance_error: str
    sextic_flag_invariance_error: str
    sextic_endpoint_prediction_error: str
    transported_sextic_column_norm: str
    expected_delay_behavior: str
    all_checks_passed: bool


AuditRecord = (
    ExactSecondGaugeTransportRecord | NumericalSecondGaugeTransportRecord
)


# Each listed word represents its Hermitian orbit.  A self-adjoint reduced
# word is inserted once; every other word is inserted together with its
# formal adjoint.
TRANSPORT_WITNESS_TERMS: WitnessTerms = (
    (Fraction(36), ""),
    (Fraction(-36), "as"),
    (Fraction(-36), "sa"),
    (Fraction(-97, 4), "aaaa"),
    (Fraction(9), "aasssa"),
    (Fraction(9), "asssaa"),
    (Fraction(-7), "saaaaa"),
    (Fraction(25, 4), "aasssaaaaaasssaa"),
    (Fraction(-6), "aaaaas"),
    (Fraction(23, 4), "aaaaasssaaaasssa"),
    (Fraction(21, 4), "aaasssaaaasssaaa"),
    (Fraction(-21, 4), "aassssssaaassssa"),
    (Fraction(21, 4), "asssaaaasssaaaaa"),
    (Fraction(-21, 4), "assssssaaassssaa"),
    (Fraction(5), "assssssa"),
    (Fraction(19, 4), "aaasssaaaaaasssa"),
    (Fraction(17, 4), "aaaasssaaaasssaa"),
    (Fraction(17, 4), "aasssaaaasssaaaa"),
    (Fraction(-17, 4), "aassssaaassssssa"),
    (Fraction(17, 4), "asssaaaaaasssaaa"),
    (Fraction(-17, 4), "assssaaassssssaa"),
    (Fraction(-3), "asssaaaasssa"),
    (Fraction(-2), "aaaaaass"),
    (Fraction(-2), "aassssaa"),
    (Fraction(-2), "aaaaasssaaaass"),
    (Fraction(2), "aaassssaaassssaa"),
    (Fraction(2), "aassssaaassssaaa"),
    (Fraction(3, 2), "aaaaaasssa"),
    (Fraction(-3, 2), "aaassssaaasssa"),
    (Fraction(-3, 2), "asssaaassssaaa"),
    (Fraction(3, 2), "aasssaaassssaaas"),
    (Fraction(3, 2), "saaassssaaasssaa"),
    (Fraction(1), "aaaaaaaa"),
    (Fraction(1), "aasssssssa"),
    (Fraction(1), "asssssssaa"),
    (Fraction(1), "aaaaasssaaas"),
    (Fraction(-1), "aaasssaaaaas"),
    (Fraction(1), "aaassssssaaa"),
    (Fraction(-1), "saaaaasssaaa"),
    (Fraction(-1), "aaaaaasssaaass"),
    (Fraction(-1), "aaaasssaaaaass"),
    (Fraction(-1), "aasssaaassssaa"),
    (Fraction(-1), "aassssaaasssaa"),
    (Fraction(-1), "asssaaassssssa"),
    (Fraction(1), "asssssaaassssa"),
    (Fraction(1, 2), "aaaasssaaa"),
    (Fraction(1, 2), "asssaaaaaa"),
    (Fraction(-1, 2), "aaaassssssaa"),
    (Fraction(-1, 2), "aassssssaaaa"),
    (Fraction(-1, 2), "aassssssssaa"),
    (Fraction(-1, 2), "saaaaaasssaa"),
    (Fraction(-1, 2), "saaasssaaaaa"),
    (Fraction(-1, 2), "aaaaaaaasssaaa"),
    (Fraction(1, 2), "aaaasssssssaaa"),
    (Fraction(-1, 2), "aaasssaaaaaaaa"),
    (Fraction(1, 2), "aaasssssssaaaa"),
    (Fraction(-1, 2), "aassssssssaaas"),
    (Fraction(1, 2), "assssaaasssssa"),
    (Fraction(-1, 2), "saaassssaaaaaa"),
    (Fraction(1, 2), "aaaasssaaaaasssa"),
    (Fraction(-1, 2), "aaaasssaaasssaaa"),
    (Fraction(1, 2), "aaasssaaaaasssaa"),
    (Fraction(-1, 2), "aaasssaaasssaaaa"),
    (Fraction(1, 2), "aaasssaaasssssaa"),
    (Fraction(1, 2), "aaasssssaaasssaa"),
    (Fraction(1, 2), "aasssaaaaaaasssa"),
    (Fraction(1, 2), "aasssaaaaasssaaa"),
    (Fraction(1, 2), "aasssaaasssssaaa"),
    (Fraction(1, 2), "aasssssaaasssaaa"),
    (Fraction(-1, 2), "asssaaasssaaaaaa"),
    (Fraction(1, 2), "asssaaasssssssaa"),
    (Fraction(1, 2), "asssssssaaaasssa"),
)


ENDPOINT_FACTOR_TERMS: WitnessTerms = (
    (Fraction(4), "aaaaa"),
    (Fraction(-4), "saaaaaa"),
    (Fraction(2), "aaaaaasssaa"),
    (Fraction(3), "aaaaasssaaa"),
    (Fraction(2), "aaaasssaaaa"),
    (Fraction(1), "aaasssaaaaa"),
    (Fraction(-2), "saaaaaaasssaa"),
    (Fraction(-3), "saaaaaasssaaa"),
    (Fraction(-2), "saaaaasssaaaa"),
    (Fraction(-1), "saaaasssaaaaa"),
)


def transport_witness() -> Polynomial:
    """Return the Hermitian metric witness for the second transport."""

    witness: Polynomial = {}
    for coefficient, word in TRANSPORT_WITNESS_TERMS:
        monomial = {word: Fraction(1)}
        witness = add(witness, scale(coefficient, monomial))
        if adjoint(monomial) != monomial:
            witness = add(
                witness,
                scale(coefficient, adjoint(monomial)),
            )
    return witness


def endpoint_factor() -> Polynomial:
    """Return the sparse left factor of the upper endpoint response."""

    factor: Polynomial = {}
    for coefficient, word in ENDPOINT_FACTOR_TERMS:
        factor = add(
            factor,
            scale(coefficient, {word: Fraction(1)}),
        )
    return factor


def full_upper_endpoint_response(witness: Polynomial) -> Polynomial:
    """Return one minus-quarter of the complete sextic endpoint change."""

    quartic_witness, _, _, _ = gauge_data(SPECIFICATIONS[0])
    metric_two = prepared_exact_metrics()[2]
    upper_complement_inverse = add(
        IDENTITY,
        scale(-1, F),
        scale(Fraction(-2, 3), E),
    )
    metric_four_change = scale(-1, quartic_witness)
    schur_change = add(
        multiply(
            multiply(
                multiply(metric_two, upper_complement_inverse),
                metric_four_change,
            ),
            F,
        ),
        multiply(
            multiply(
                multiply(metric_four_change, upper_complement_inverse),
                metric_two,
            ),
            F,
        ),
    )
    return add(
        multiply(multiply(F, witness), F),
        multiply(F, schur_change),
    )


def second_transport_data() -> tuple[
    Polynomial,
    Polynomial,
    Polynomial,
    Polynomial,
    Polynomial,
]:
    """Return forcing, witness, corrections, and transported sextic lift."""

    forcing, parallel_correction = transported_sextic_forcing()
    witness = transport_witness()
    coboundary = add(
        witness,
        scale(-1, multiply(multiply(STAR, witness), S)),
    )
    perpendicular_correction = multiply(
        multiply(
            Q,
            add(coboundary, scale(-1, forcing)),
        ),
        E,
    )
    _, _, _, normalized_sextic = gauge_data(SPECIFICATIONS[2])
    transported_sextic = add(
        normalized_sextic,
        parallel_correction,
        perpendicular_correction,
    )
    return (
        forcing,
        witness,
        parallel_correction,
        perpendicular_correction,
        transported_sextic,
    )


def exact_record() -> ExactSecondGaugeTransportRecord:
    """Verify the complete rational word certificate."""

    (
        forcing,
        witness,
        parallel,
        correction,
        transported,
    ) = second_transport_data()
    coboundary = add(
        witness,
        scale(-1, multiply(multiply(STAR, witness), S)),
    )
    endpoint_response = full_upper_endpoint_response(witness)
    first_channel = multiply(multiply(E, S), F)
    factored_endpoint = multiply(endpoint_factor(), first_channel)
    triple_delay = DelayedQuotient(4)
    first_delay = DelayedQuotient(2)
    correction_quotient = triple_delay.reduce(correction)
    transported_quotient = triple_delay.reduce(transported)
    endpoint_quotient = first_delay.reduce(endpoint_response)
    residuals = (
        add(witness, scale(-1, adjoint(witness))),
        multiply(multiply(E, witness), E),
        add(
            forcing,
            correction,
            adjoint(correction),
            scale(-1, coboundary),
        ),
        multiply(E, correction),
        correction_quotient,
        transported_quotient,
        endpoint_quotient,
        add(
            endpoint_response,
            scale(
                -1,
                add(factored_endpoint, adjoint(factored_endpoint)),
            ),
        ),
    )
    maximum_residual = max(map(len, residuals))
    verified = maximum_residual == 0
    if not verified:
        raise RuntimeError(
            "the exact second gauge-transport certificate failed: "
            f"residuals={tuple(map(len, residuals))}"
        )
    return ExactSecondGaugeTransportRecord(
        record_type="exact_flag_preserving_transport",
        transport_witness_word_count=len(witness),
        transport_witness_coefficient_l1=format_fraction(
            coefficient_l1(witness)
        ),
        perpendicular_correction_word_count=len(correction),
        perpendicular_correction_coefficient_l1=format_fraction(
            coefficient_l1(correction)
        ),
        parallel_correction_word_count=len(parallel),
        parallel_correction_coefficient_l1=format_fraction(
            coefficient_l1(parallel)
        ),
        transported_sextic_word_count=len(transported),
        transported_sextic_coefficient_l1=format_fraction(
            coefficient_l1(transported)
        ),
        endpoint_response_word_count=len(endpoint_response),
        endpoint_response_coefficient_l1=format_fraction(
            coefficient_l1(endpoint_response)
        ),
        endpoint_factor_word_count=len(endpoint_factor()),
        endpoint_factor_coefficient_l1=format_fraction(
            coefficient_l1(endpoint_factor())
        ),
        maximum_exact_residual_word_count=maximum_residual,
        correction_triple_delay_quotient_word_count=len(
            correction_quotient
        ),
        transported_triple_delay_quotient_word_count=len(
            transported_quotient
        ),
        endpoint_first_delay_quotient_word_count=len(
            endpoint_quotient
        ),
        all_checks_passed=verified,
    )


def numerical_record(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> NumericalSecondGaugeTransportRecord:
    """Evaluate the transported series and its surviving endpoint flag."""

    (
        forcing,
        witness,
        parallel,
        correction,
        transported,
    ) = second_transport_data()
    _, _, quartic_gauge, _ = gauge_data(SPECIFICATIONS[0])
    _, _, quintic_gauge, _ = gauge_data(SPECIFICATIONS[1])
    _, _, sextic_gauge, _ = gauge_data(SPECIFICATIONS[2])
    _, _, quintic_correction, _ = first_transport_data()

    forcing_matrix = evaluate_polynomial(
        add(forcing, correction, adjoint(correction)),
        operator,
    )
    witness_matrix = evaluate_polynomial(witness, operator)
    reconstructed = stein_inverse(operator, forcing_matrix)
    witness_error = float(np.linalg.norm(reconstructed - witness_matrix))
    lower_witness = float(
        np.linalg.norm(right.conj().T @ reconstructed @ right)
    )
    upper_witness = float(
        np.linalg.norm(left.conj().T @ reconstructed @ left)
    )

    first = transfer_coefficient(operator, right, left, 1)
    second = transfer_coefficient(operator, right, left, 2)
    flag = joint_left_kernel(first, second)
    flag_witness = float(
        np.linalg.norm(
            flag.conj().T
            @ left.conj().T
            @ reconstructed
            @ left
            @ flag
        )
    )

    original_sextic = sextic_column(operator, right)
    original_upper, original_lower, _ = prepared_endpoint_series(
        operator,
        right,
        left,
        6,
        additional_factor_columns={6: original_sextic},
    )
    changes = {
        4: -evaluate_polynomial(quartic_gauge, operator) @ right,
        5: evaluate_polynomial(
            add(scale(-1, quintic_gauge), quintic_correction),
            operator,
        )
        @ right,
        6: original_sextic
        + evaluate_polynomial(
            add(
                scale(-1, sextic_gauge),
                parallel,
                correction,
            ),
            operator,
        )
        @ right,
    }
    changed_upper, changed_lower, _ = prepared_endpoint_series(
        operator,
        right,
        left,
        6,
        additional_factor_columns=changes,
    )
    nonlinear_lower_change = float(
        np.linalg.norm(changed_lower[6] - original_lower[6])
    )
    lower_order_upper_error = max(
        float(
            np.linalg.norm(
                changed_upper[degree] - original_upper[degree]
            )
        )
        for degree in range(6)
    )
    sextic_difference = changed_upper[6] - original_upper[6]
    sextic_flag_error = float(
        np.linalg.norm(
            flag.conj().T @ sextic_difference @ flag
        )
    )
    predicted_sextic_difference = (
        -4
        * left.conj().T
        @ evaluate_polynomial(
            full_upper_endpoint_response(witness),
            operator,
        )
        @ left
    )
    endpoint_prediction_error = float(
        np.linalg.norm(
            sextic_difference - predicted_sextic_difference
        )
    )
    transported_column_norm = float(
        np.linalg.norm(
            evaluate_polynomial(transported, operator) @ right
        )
    )
    first_norm = float(np.linalg.norm(first))
    second_norm = float(np.linalg.norm(second))

    tolerance = 2e-7
    common_checks = bool(
        witness_error < tolerance
        and lower_witness < tolerance
        and flag_witness < tolerance
        and lower_order_upper_error < tolerance
        and sextic_flag_error < tolerance
        and endpoint_prediction_error < tolerance
    )
    if construction_kind == "triple_complete_delay":
        behavior = "first_three_normalized_columns_zero"
        behavior_check = bool(
            first_norm < tolerance
            and second_norm < tolerance
            and transported_column_norm < tolerance
        )
    else:
        behavior = "sextic_motion_zero_only_on_surviving_flag"
        behavior_check = True

    verified = common_checks and behavior_check
    if not verified:
        raise RuntimeError(
            "the numerical second gauge-transport audit failed: "
            f"kind={construction_kind}, "
            f"witness={witness_error:.3e}, "
            f"lower_metric={lower_witness:.3e}, "
            f"lower_order={lower_order_upper_error:.3e}, "
            f"flag={sextic_flag_error:.3e}, "
            f"prediction={endpoint_prediction_error:.3e}, "
            f"column={transported_column_norm:.3e}"
        )
    return NumericalSecondGaugeTransportRecord(
        record_type=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        first_transfer_norm=format_float(first_norm),
        second_transfer_norm=format_float(second_norm),
        surviving_flag_dimension=flag.shape[1],
        witness_stein_error=format_float(witness_error),
        witness_lower_endpoint_norm=format_float(lower_witness),
        witness_upper_endpoint_norm=format_float(upper_witness),
        witness_flag_endpoint_norm=format_float(flag_witness),
        nonlinear_lower_sixth_change_norm=format_float(
            nonlinear_lower_change
        ),
        lower_order_upper_invariance_error=format_float(
            lower_order_upper_error
        ),
        sextic_flag_invariance_error=format_float(sextic_flag_error),
        sextic_endpoint_prediction_error=format_float(
            endpoint_prediction_error
        ),
        transported_sextic_column_norm=format_float(
            transported_column_norm
        ),
        expected_delay_behavior=behavior,
        all_checks_passed=verified,
    )


def standard_records() -> list[AuditRecord]:
    """Return exact, unstructured, partial-flag, and delayed audits."""

    records: list[AuditRecord] = [exact_record()]
    operator, right, left = random_partial_isometry(
        12,
        2,
        np.random.default_rng(702_271),
    )
    records.append(
        numerical_record("unstructured", operator, right, left)
    )

    for multiplicity, seed in (
        (3, 702_272),
        (4, 702_273),
        (5, 702_274),
    ):
        operator, right, left, _, _ = rank_chain_case(
            multiplicity,
            seed,
            0.8,
        )
        records.append(
            numerical_record(
                f"rank_changing_flag_{multiplicity}",
                operator,
                right,
                left,
            )
        )

    operator, right, left, _ = inflated_case(
        9,
        2,
        4,
        2,
        702_275,
    )
    records.append(
        numerical_record(
            "triple_complete_delay",
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
            "repeated_crabb_second_gauge_transport_s70225.jsonl"
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
