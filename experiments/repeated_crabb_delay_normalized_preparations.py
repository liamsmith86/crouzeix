#!/usr/bin/env python3
"""Audit endpoint-null delay normalizations of the prepared columns."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_canonical_quartic_preimage import (
    evaluate_polynomial,
    quartic_column_lift,
)
from repeated_crabb_canonical_quintic_preimage import (
    joint_left_kernel,
    prepared_endpoint_series,
    quintic_column_lift,
)
from repeated_crabb_canonical_sextic_preimage import (
    sextic_column_lift,
)
from repeated_crabb_cyclic_radial_volume import DelayedQuotient
from repeated_crabb_boundary_metric_flag import rank_chain_case
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
from repeated_crabb_edge_deleted_balanced_volume import delayed_case
from repeated_crabb_elliptic_cokernel import (
    endpoint_motion,
    stein_inverse,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
WitnessTerms = tuple[tuple[Fraction, str], ...]


@dataclass(frozen=True)
class GaugeSpecification:
    """One prepared column and its endpoint-null witness."""

    order: int
    delay_length: int
    column_lift: Callable[[], Polynomial]
    witness_terms: WitnessTerms


@dataclass(frozen=True)
class DelayNormalizedPreparationRecord:
    """One exact or numerical delay-normalization audit."""

    record_type: str
    order: int
    delay_length: int
    state_dimension: int | str
    defect_dimension: int | str
    witness_coefficient_l1: str
    gauge_column_coefficient_l1: str
    normalized_column_coefficient_l1: str
    maximum_exact_residual_word_count: int | str
    stein_reconstruction_error: str
    lower_endpoint_error: str
    upper_endpoint_error: str
    endpoint_response_difference_error: str
    normalized_delay_column_error: str
    simultaneous_quintic_flag_norm: str
    all_checks_passed: bool


SPECIFICATIONS = (
    GaugeSpecification(
        order=4,
        delay_length=1,
        column_lift=quartic_column_lift,
        witness_terms=(
            (Fraction(1), "aaaa"),
            (Fraction(1, 2), "aaaaaasssa"),
            (Fraction(1, 2), "aaaaasssaa"),
            (Fraction(1, 2), "aaaasssaaa"),
            (Fraction(1, 2), "aaasssaaaa"),
            (Fraction(1, 2), "aasssaaaaa"),
            (Fraction(-1, 2), "aasssssssa"),
            (Fraction(1, 2), "asssaaaaaa"),
            (Fraction(-1, 2), "asssssssaa"),
        ),
    ),
    GaugeSpecification(
        order=5,
        delay_length=2,
        column_lift=quintic_column_lift,
        witness_terms=(
            (Fraction(45, 8), "aa"),
            (Fraction(15, 8), "aaas"),
            (Fraction(15, 8), "saaa"),
            (Fraction(-5), "aaaaaa"),
            (Fraction(-15, 8), "assssa"),
            (Fraction(-1, 2), "aaaaaaas"),
            (Fraction(3, 2), "aaaasssa"),
            (Fraction(3, 2), "aaasssaa"),
            (Fraction(3, 2), "aasssaaa"),
            (Fraction(-3, 2), "aasssssa"),
            (Fraction(3, 2), "asssaaaa"),
            (Fraction(-3, 2), "asssssaa"),
            (Fraction(-1, 2), "saaaaaaa"),
            (Fraction(1, 2), "assssssssa"),
        ),
    ),
    GaugeSpecification(
        order=6,
        delay_length=3,
        column_lift=sextic_column_lift,
        witness_terms=(
            (Fraction(9, 2), "aaaa"),
            (Fraction(11, 4), "aaaaas"),
            (Fraction(11, 4), "saaaaa"),
            (Fraction(5), "aaaaaaaa"),
            (Fraction(-11, 4), "assssssa"),
            (Fraction(-5), "aaaaaasssa"),
            (Fraction(-5), "aaaaasssaa"),
            (Fraction(-5), "aaaasssaaa"),
            (Fraction(-5), "aaasssaaaa"),
            (Fraction(-5), "aasssaaaaa"),
            (Fraction(5), "aasssssssa"),
            (Fraction(-5), "asssaaaaaa"),
            (Fraction(5), "asssssssaa"),
        ),
    ),
)


def coefficient_l1(polynomial: Polynomial) -> Fraction:
    """Return the exact coefficient l1 norm."""

    return sum((abs(value) for value in polynomial.values()), Fraction())


def format_fraction(value: Fraction) -> str:
    """Format one exact rational number deterministically."""

    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def hermitian_witness(specification: GaugeSpecification) -> Polynomial:
    """Return the explicit Hermitian endpoint-null witness."""

    witness: Polynomial = {}
    for coefficient, word in specification.witness_terms:
        monomial = {word: Fraction(1)}
        witness = add(
            witness,
            scale(coefficient, monomial),
            scale(coefficient, adjoint(monomial)),
        )
    return witness


def gauge_data(
    specification: GaugeSpecification,
) -> tuple[Polynomial, Polynomial, Polynomial, Polynomial]:
    """Return witness, forcing, gauge lift, and normalized lift."""

    witness = hermitian_witness(specification)
    forcing = add(
        witness,
        scale(-1, multiply(multiply(STAR, witness), S)),
    )
    gauge_lift = multiply(multiply(Q, forcing), E)
    normalized_lift = add(
        specification.column_lift(),
        scale(-1, gauge_lift),
    )
    return witness, forcing, gauge_lift, normalized_lift


def exact_record(
    specification: GaugeSpecification,
) -> DelayNormalizedPreparationRecord:
    """Verify one gauge certificate in exact word algebra."""

    witness, forcing, gauge_lift, normalized_lift = gauge_data(
        specification
    )
    quotient = DelayedQuotient(specification.delay_length + 1)
    residuals = (
        add(witness, scale(-1, adjoint(witness))),
        multiply(multiply(E, witness), E),
        multiply(multiply(F, witness), F),
        multiply(multiply(Q, forcing), Q),
        multiply(multiply(E, forcing), E),
        add(
            forcing,
            scale(-1, gauge_lift),
            scale(-1, adjoint(gauge_lift)),
        ),
        multiply(E, normalized_lift),
        quotient.reduce(normalized_lift),
    )
    maximum_residual = max(map(len, residuals))
    verified = maximum_residual == 0
    if not verified:
        raise RuntimeError(
            "delay-normalized exact certificate failed: "
            f"order={specification.order}, "
            f"residuals={tuple(map(len, residuals))}"
        )
    return DelayNormalizedPreparationRecord(
        record_type="exact_word_certificate",
        order=specification.order,
        delay_length=specification.delay_length,
        state_dimension="not_applicable",
        defect_dimension="not_applicable",
        witness_coefficient_l1=format_fraction(
            coefficient_l1(witness)
        ),
        gauge_column_coefficient_l1=format_fraction(
            coefficient_l1(gauge_lift)
        ),
        normalized_column_coefficient_l1=format_fraction(
            coefficient_l1(normalized_lift)
        ),
        maximum_exact_residual_word_count=maximum_residual,
        stein_reconstruction_error="not_applicable",
        lower_endpoint_error="not_applicable",
        upper_endpoint_error="not_applicable",
        endpoint_response_difference_error="not_applicable",
        normalized_delay_column_error="not_applicable",
        simultaneous_quintic_flag_norm="not_applicable",
        all_checks_passed=verified,
    )


def numerical_record(
    specification: GaugeSpecification,
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> DelayNormalizedPreparationRecord:
    """Independently evaluate one endpoint-null gauge."""

    witness, forcing, gauge_lift, normalized_lift = gauge_data(
        specification
    )
    witness_matrix = evaluate_polynomial(witness, operator)
    forcing_matrix = evaluate_polynomial(forcing, operator)
    gauge_column = evaluate_polynomial(gauge_lift, operator) @ right
    normalized_column = (
        evaluate_polynomial(normalized_lift, operator) @ right
    )
    original_column = (
        evaluate_polynomial(
            specification.column_lift(),
            operator,
        )
        @ right
    )
    column_forcing = (
        right @ gauge_column.conj().T
        + gauge_column @ right.conj().T
    )
    reconstructed = stein_inverse(operator, column_forcing)

    stein_error = float(
        np.linalg.norm(reconstructed - witness_matrix)
        + np.linalg.norm(column_forcing - forcing_matrix)
    )
    lower_error = float(
        np.linalg.norm(
            right.conj().T @ reconstructed @ right
        )
    )
    upper_error = float(
        np.linalg.norm(
            left.conj().T @ reconstructed @ left
        )
    )
    response_error = float(
        np.linalg.norm(
            endpoint_motion(
                operator,
                right,
                left,
                original_column,
            )
            - endpoint_motion(
                operator,
                right,
                left,
                normalized_column,
            )
        )
    )
    delay_error = (
        float(np.linalg.norm(normalized_column))
        if construction_kind == "complete_delay"
        else 0.0
    )
    tolerance = 4e-8
    verified = bool(
        stein_error < tolerance
        and lower_error < tolerance
        and upper_error < tolerance
        and response_error < tolerance
        and delay_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "delay-normalized numerical audit failed: "
            f"order={specification.order}, "
            f"kind={construction_kind}, "
            f"stein={stein_error:.3e}, "
            f"lower={lower_error:.3e}, "
            f"upper={upper_error:.3e}, "
            f"response={response_error:.3e}, "
            f"delay={delay_error:.3e}"
        )
    return DelayNormalizedPreparationRecord(
        record_type=f"numerical_{construction_kind}",
        order=specification.order,
        delay_length=specification.delay_length,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        witness_coefficient_l1=format_fraction(
            coefficient_l1(witness)
        ),
        gauge_column_coefficient_l1=format_fraction(
            coefficient_l1(gauge_lift)
        ),
        normalized_column_coefficient_l1=format_fraction(
            coefficient_l1(normalized_lift)
        ),
        maximum_exact_residual_word_count="not_applicable",
        stein_reconstruction_error=format_float(stein_error),
        lower_endpoint_error=format_float(lower_error),
        upper_endpoint_error=format_float(upper_error),
        endpoint_response_difference_error=format_float(
            response_error
        ),
        normalized_delay_column_error=format_float(delay_error),
        simultaneous_quintic_flag_norm="not_applicable",
        all_checks_passed=verified,
    )


def simultaneous_scope_guard_record() -> DelayNormalizedPreparationRecord:
    """Reject naive simultaneous reuse of the stored later columns."""

    operator, right, left, _, _ = rank_chain_case(
        4,
        109_204,
        1.0,
    )
    specifications = {
        specification.order: specification
        for specification in SPECIFICATIONS
    }
    additional_columns = {}
    for order in (4, 5):
        _, _, gauge_lift, _ = gauge_data(specifications[order])
        additional_columns[order] = (
            -evaluate_polynomial(gauge_lift, operator) @ right
        )
    upper, _, _ = prepared_endpoint_series(
        operator,
        right,
        left,
        additional_factor_columns=additional_columns,
    )
    first = transfer_coefficient(
        operator,
        right,
        left,
        1,
    )
    second = transfer_coefficient(
        operator,
        right,
        left,
        2,
    )
    flag = joint_left_kernel(first, second)
    flag_norm = float(
        np.linalg.norm(flag.conj().T @ upper[5] @ flag)
    )
    verified = flag_norm > 1e-4
    if not verified:
        raise RuntimeError(
            "the simultaneous-normalization scope guard regressed: "
            f"flag_norm={flag_norm:.3e}"
        )
    return DelayNormalizedPreparationRecord(
        record_type="simultaneous_reuse_scope_guard",
        order=5,
        delay_length=2,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        witness_coefficient_l1="not_applicable",
        gauge_column_coefficient_l1="not_applicable",
        normalized_column_coefficient_l1="not_applicable",
        maximum_exact_residual_word_count="not_applicable",
        stein_reconstruction_error="not_applicable",
        lower_endpoint_error="not_applicable",
        upper_endpoint_error="not_applicable",
        endpoint_response_difference_error="not_applicable",
        normalized_delay_column_error="not_applicable",
        simultaneous_quintic_flag_norm=format_float(flag_norm),
        all_checks_passed=verified,
    )


def standard_records() -> list[DelayNormalizedPreparationRecord]:
    """Return exact, unstructured, and complete-delay audits."""

    records: list[DelayNormalizedPreparationRecord] = []
    for specification in SPECIFICATIONS:
        records.append(exact_record(specification))
        generator = np.random.default_rng(
            702_250 + specification.order
        )
        operator, right, left = random_partial_isometry(
            9 + specification.order,
            2,
            generator,
        )
        records.append(
            numerical_record(
                specification,
                "unstructured",
                operator,
                right,
                left,
            )
        )
        operator, right, left = delayed_case(
            specification.delay_length + 1
        )
        records.append(
            numerical_record(
                specification,
                "complete_delay",
                operator,
                right,
                left,
            )
        )
    records.append(simultaneous_scope_guard_record())
    return records


def write_records(
    records: list[DelayNormalizedPreparationRecord],
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
            "repeated_crabb_delay_normalized_preparations_s70225.jsonl"
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
