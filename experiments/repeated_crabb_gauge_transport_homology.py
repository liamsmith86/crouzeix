#!/usr/bin/env python3
"""Audit the triangular transport equation for endpoint-null gauges."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_canonical_quartic_preimage import (
    evaluate_polynomial,
)
from repeated_crabb_canonical_quintic_preimage import (
    exact_canonical_factor_lifts,
)
from repeated_crabb_cyclic_radial_volume import DelayedQuotient
from repeated_crabb_delay_normalized_preparations import (
    SPECIFICATIONS,
    gauge_data,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    Q,
    Polynomial,
    add,
    adjoint,
    multiply,
    operator_series,
    scale,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class ExactGaugeTransportRecord:
    """Exact first-successor transport and ideal-obstruction audit."""

    record_type: str
    gauge_order: int
    successor_order: int
    recurrence_residual_word_count: int
    hermitian_residual_word_count: int
    lower_corner_word_count: int
    transport_word_count: int
    transport_coefficient_l1: str
    initial_delay_quotient_word_count: int
    initial_delay_quotient_coefficient_l1: str
    initial_delay_q_corner_word_count: int
    initial_delay_q_corner_coefficient_l1: str
    naive_ideal_invariance_holds: bool
    all_checks_passed: bool


@dataclass(frozen=True)
class NumericalGaugeTransportRecord:
    """One independent matrix evaluation of the successor transport."""

    record_type: str
    state_dimension: int
    defect_dimension: int
    first_transfer_norm: str
    second_transfer_norm: str
    stein_residual_norm: str
    lower_endpoint_norm: str
    upper_endpoint_norm: str
    upper_endpoint_trace_absolute: str
    expected_delay_behavior: str
    all_checks_passed: bool


AuditRecord = ExactGaugeTransportRecord | NumericalGaugeTransportRecord


def coefficient_l1(polynomial: Polynomial) -> Fraction:
    """Return the exact coefficient l1 norm."""

    return sum((abs(value) for value in polynomial.values()), Fraction())


def format_fraction(value: Fraction) -> str:
    """Format one exact rational number deterministically."""

    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def first_successor_transport() -> Polynomial:
    """Return the raw order-five transport from the order-four gauge."""

    specification = SPECIFICATIONS[0]
    witness, _, gauge_lift, _ = gauge_data(specification)
    factor_one = exact_canonical_factor_lifts(5)[1]
    operator = operator_series(1)
    factor_transport = add(
        scale(-1, multiply(gauge_lift, adjoint(factor_one))),
        scale(-1, multiply(factor_one, adjoint(gauge_lift))),
    )
    metric_transport = add(
        scale(
            -1,
            multiply(
                multiply(adjoint(operator[1]), witness),
                operator[0],
            ),
        ),
        scale(
            -1,
            multiply(
                multiply(adjoint(operator[0]), witness),
                operator[1],
            ),
        ),
    )
    return add(factor_transport, metric_transport)


def independently_enumerated_successor() -> Polynomial:
    """Recompute the successor by coefficient convolution."""

    specification = SPECIFICATIONS[0]
    witness, _, gauge_lift, _ = gauge_data(specification)
    factors = exact_canonical_factor_lifts(5)
    factor_change = [{} for _ in factors]
    factor_change[specification.order] = scale(-1, gauge_lift)

    slack_change: Polynomial = {}
    successor = specification.order + 1
    for left_degree in range(successor + 1):
        right_degree = successor - left_degree
        slack_change = add(
            slack_change,
            multiply(
                factor_change[left_degree],
                adjoint(factors[right_degree]),
            ),
            multiply(
                factors[left_degree],
                adjoint(factor_change[right_degree]),
            ),
            multiply(
                factor_change[left_degree],
                adjoint(factor_change[right_degree]),
            ),
        )

    operator = operator_series(successor)
    metric_change = [{} for _ in range(successor + 1)]
    metric_change[specification.order] = scale(-1, witness)
    metric_image_change: Polynomial = {}
    for left_degree in range(successor + 1):
        for metric_degree in range(successor):
            right_degree = successor - left_degree - metric_degree
            if 0 <= right_degree <= successor:
                metric_image_change = add(
                    metric_image_change,
                    multiply(
                        multiply(
                            adjoint(operator[left_degree]),
                            metric_change[metric_degree],
                        ),
                        operator[right_degree],
                    ),
                )
    return add(slack_change, metric_image_change)


def exact_record() -> ExactGaugeTransportRecord:
    """Verify the recurrence and reject raw delay-ideal invariance."""

    transport = first_successor_transport()
    independent = independently_enumerated_successor()
    quotient = DelayedQuotient(2)
    reduced_transport = quotient.reduce(transport)
    reduced_q_corner = quotient.reduce(
        multiply(multiply(Q, transport), Q)
    )
    recurrence_residual = add(transport, scale(-1, independent))
    hermitian_residual = add(transport, scale(-1, adjoint(transport)))
    lower_corner = multiply(multiply(E, transport), E)

    verified = bool(
        not recurrence_residual
        and not hermitian_residual
        and not lower_corner
        and reduced_transport
        and reduced_q_corner
    )
    if not verified:
        raise RuntimeError(
            "the exact gauge-transport audit failed: "
            f"recurrence={len(recurrence_residual)}, "
            f"Hermitian={len(hermitian_residual)}, "
            f"lower={len(lower_corner)}, "
            f"delay={len(reduced_transport)}, "
            f"q_corner={len(reduced_q_corner)}"
        )
    return ExactGaugeTransportRecord(
        record_type="exact_first_successor",
        gauge_order=4,
        successor_order=5,
        recurrence_residual_word_count=len(recurrence_residual),
        hermitian_residual_word_count=len(hermitian_residual),
        lower_corner_word_count=len(lower_corner),
        transport_word_count=len(transport),
        transport_coefficient_l1=format_fraction(
            coefficient_l1(transport)
        ),
        initial_delay_quotient_word_count=len(reduced_transport),
        initial_delay_quotient_coefficient_l1=format_fraction(
            coefficient_l1(reduced_transport)
        ),
        initial_delay_q_corner_word_count=len(reduced_q_corner),
        initial_delay_q_corner_coefficient_l1=format_fraction(
            coefficient_l1(reduced_q_corner)
        ),
        naive_ideal_invariance_holds=False,
        all_checks_passed=verified,
    )


def numerical_record(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> NumericalGaugeTransportRecord:
    """Evaluate the first-successor transport on one colligation."""

    transport = evaluate_polynomial(
        first_successor_transport(),
        operator,
    )
    metric = stein_inverse(operator, transport)
    stein_residual = float(
        np.linalg.norm(
            metric
            - operator.conj().T @ metric @ operator
            - transport
        )
    )
    lower = right.conj().T @ metric @ right
    upper = left.conj().T @ metric @ left
    lower_norm = float(np.linalg.norm(lower))
    upper_norm = float(np.linalg.norm(upper))
    upper_trace = float(abs(np.trace(upper)))
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

    tolerance = 3e-8
    common_checks = bool(
        stein_residual < tolerance
        and lower_norm < tolerance
        and upper_trace < tolerance
    )
    if construction_kind == "unstructured":
        behavior = "nonzero_endpoint_allowed"
        behavior_check = upper_norm > 1e-5
    elif construction_kind == "single_complete_delay":
        behavior = "first_transfer_zero_endpoint_nonzero"
        behavior_check = bool(
            first_transfer < tolerance
            and second_transfer > 1e-4
            and upper_norm > 1e-5
        )
    elif construction_kind == "double_complete_delay":
        behavior = "first_two_transfers_zero_endpoint_zero"
        behavior_check = bool(
            first_transfer < tolerance
            and second_transfer < tolerance
            and upper_norm < tolerance
        )
    else:
        raise ValueError(f"unknown construction kind: {construction_kind}")

    verified = common_checks and behavior_check
    if not verified:
        raise RuntimeError(
            "the numerical gauge-transport audit failed: "
            f"kind={construction_kind}, "
            f"stein={stein_residual:.3e}, "
            f"lower={lower_norm:.3e}, "
            f"upper={upper_norm:.3e}, "
            f"trace={upper_trace:.3e}, "
            f"B1={first_transfer:.3e}, "
            f"B2={second_transfer:.3e}"
        )
    return NumericalGaugeTransportRecord(
        record_type=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        first_transfer_norm=format_float(first_transfer),
        second_transfer_norm=format_float(second_transfer),
        stein_residual_norm=format_float(stein_residual),
        lower_endpoint_norm=format_float(lower_norm),
        upper_endpoint_norm=format_float(upper_norm),
        upper_endpoint_trace_absolute=format_float(upper_trace),
        expected_delay_behavior=behavior,
        all_checks_passed=verified,
    )


def standard_records() -> list[AuditRecord]:
    """Return exact and independent matrix audits."""

    records: list[AuditRecord] = [exact_record()]
    operator, right, left = random_partial_isometry(
        12,
        2,
        np.random.default_rng(702_251),
    )
    records.append(
        numerical_record(
            "unstructured",
            operator,
            right,
            left,
        )
    )
    for grade, construction_kind, seed in (
        (2, "single_complete_delay", 702_252),
        (3, "double_complete_delay", 702_253),
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
            "repeated_crabb_gauge_transport_homology_s70225.jsonl"
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
