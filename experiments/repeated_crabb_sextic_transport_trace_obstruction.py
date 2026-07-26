#!/usr/bin/env python3
"""Certify the scalar trace obstruction to endpoint-null sextic transport."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp

from repeated_crabb_canonical_quintic_preimage import (
    exact_canonical_factor_lifts,
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
    Polynomial,
    add,
    adjoint,
    multiply,
    operator_series,
    scale,
)
from repeated_crabb_first_gauge_transport import (
    first_transport_data,
)


@dataclass(frozen=True)
class SexticTransportTraceObstructionRecord:
    """Exact scalar obstruction to a globally endpoint-null successor."""

    record_type: str
    state_dimension: int
    defect_dimension: int
    spectral_radius: str
    first_transfer_coefficient: str
    normalized_forcing_word_count: int
    normalized_forcing_coefficient_l1: str
    lower_parallel_word_count: int
    lower_parallel_coefficient_l1: str
    lower_parallel_triple_delay_quotient_word_count: int
    normalized_lower_corner_word_count: int
    exact_stein_residual_nonzero_entry_count: int
    exact_upper_endpoint: str
    endpoint_null_transport_possible: bool
    all_checks_passed: bool


def transported_sextic_forcing() -> tuple[
    Polynomial,
    Polynomial,
]:
    """Return the lower-normalized forcing and its parallel correction."""

    quartic_witness, _, quartic_gauge, _ = gauge_data(
        SPECIFICATIONS[0]
    )
    _, _, quintic_gauge, _ = gauge_data(SPECIFICATIONS[1])
    _, _, sextic_gauge, _ = gauge_data(SPECIFICATIONS[2])
    _, quintic_transport_witness, quintic_correction, _ = (
        first_transport_data()
    )

    factor_changes = {
        4: scale(-1, quartic_gauge),
        5: add(scale(-1, quintic_gauge), quintic_correction),
    }
    metric_changes = {
        4: scale(-1, quartic_witness),
        5: quintic_transport_witness,
    }
    factors = exact_canonical_factor_lifts(6)
    operator = operator_series(6)

    raw_transport: Polynomial = {}
    order = 6
    for left_degree in range(1, order):
        right_degree = order - left_degree
        left_change = factor_changes.get(left_degree)
        right_change = factor_changes.get(right_degree)
        if left_change is not None:
            raw_transport = add(
                raw_transport,
                multiply(
                    left_change,
                    adjoint(factors[right_degree]),
                ),
            )
        if right_change is not None:
            raw_transport = add(
                raw_transport,
                multiply(
                    factors[left_degree],
                    adjoint(right_change),
                ),
            )
        if left_change is not None and right_change is not None:
            raw_transport = add(
                raw_transport,
                multiply(left_change, adjoint(right_change)),
            )

    for left_degree in range(order + 1):
        for metric_degree, metric_change in metric_changes.items():
            right_degree = order - left_degree - metric_degree
            if 0 <= right_degree <= order:
                raw_transport = add(
                    raw_transport,
                    multiply(
                        multiply(
                            adjoint(operator[left_degree]),
                            metric_change,
                        ),
                        operator[right_degree],
                    ),
                )

    forcing = add(
        raw_transport,
        scale(-1, sextic_gauge),
        scale(-1, adjoint(sextic_gauge)),
    )
    lower_corner = multiply(multiply(E, forcing), E)
    parallel_correction = scale(Fraction(-1, 2), lower_corner)
    lower_normalized = add(forcing, scale(-1, lower_corner))
    return lower_normalized, parallel_correction


def rational_counterexample() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return the exact partial isometry and its two defect frames."""

    operator = sp.Matrix(
        [
            [0, sp.Rational(7, 9), sp.Rational(-4, 9), sp.Rational(-4, 9)],
            [0, sp.Rational(-4, 9), sp.Rational(1, 9), sp.Rational(-8, 9)],
            [0, sp.Rational(-4, 9), sp.Rational(-8, 9), sp.Rational(1, 9)],
            [0, 0, 0, 0],
        ]
    )
    right = sp.eye(4)[:, 0]
    left = sp.eye(4)[:, 3]
    return operator, right, left


def evaluate_polynomial_exact(
    polynomial: Polynomial,
    operator: sp.Matrix,
) -> sp.Matrix:
    """Evaluate one rational word polynomial exactly."""

    result = sp.zeros(operator.rows)
    for word, coefficient in polynomial.items():
        term = sp.eye(operator.rows)
        for letter in word:
            term *= operator if letter == "s" else operator.T
        result += sp.Rational(
            coefficient.numerator,
            coefficient.denominator,
        ) * term
    return sp.simplify(result)


def solve_stein_exact(
    operator: sp.Matrix,
    forcing: sp.Matrix,
) -> sp.Matrix:
    """Solve ``X-S.T X S=forcing`` over the rationals."""

    size = operator.rows
    variables = sp.symbols(f"x0:{size * size}")
    metric = sp.Matrix(size, size, variables)
    equations = list(metric - operator.T * metric * operator - forcing)
    solution = sp.solve(equations, variables, dict=True)
    if len(solution) != 1:
        raise RuntimeError("the exact Stein equation was not unique")
    return sp.simplify(metric.subs(solution[0]))


def exact_record() -> SexticTransportTraceObstructionRecord:
    """Build and verify the rational obstruction certificate."""

    forcing_polynomial, parallel = transported_sextic_forcing()
    quotient = DelayedQuotient(4)
    parallel_quotient = quotient.reduce(parallel)
    lower_corner = multiply(multiply(E, forcing_polynomial), E)

    operator, right, left = rational_counterexample()
    identity = sp.eye(4)
    initial = identity - right * right.T
    final = identity - left * left.T
    partial_errors = (
        operator.T * operator - initial,
        operator * operator.T - final,
        right.T * left,
    )
    forcing = evaluate_polynomial_exact(forcing_polynomial, operator)
    metric = solve_stein_exact(operator, forcing)
    stein_residual = sp.simplify(
        metric - operator.T * metric * operator - forcing
    )
    endpoint = sp.factor((left.T * metric * left)[0])
    first_transfer = sp.factor((right.T * operator * left)[0])
    characteristic = sp.factor(operator.charpoly().as_expr())

    verified = bool(
        all(error == sp.zeros(*error.shape) for error in partial_errors[:2])
        and partial_errors[2] == sp.zeros(1, 1)
        and not lower_corner
        and not parallel_quotient
        and stein_residual == sp.zeros(4)
        and endpoint == sp.Rational(151040, 177147)
        and first_transfer == sp.Rational(-4, 9)
        and characteristic
        == sp.Symbol("lambda") ** 2
        * (3 * sp.Symbol("lambda") + 2) ** 2
        / 9
    )
    if not verified:
        raise RuntimeError(
            "the sextic transport trace obstruction failed: "
            f"endpoint={endpoint}, transfer={first_transfer}, "
            f"characteristic={characteristic}"
        )
    return SexticTransportTraceObstructionRecord(
        record_type="exact_scalar_trace_obstruction",
        state_dimension=4,
        defect_dimension=1,
        spectral_radius="2/3",
        first_transfer_coefficient=str(first_transfer),
        normalized_forcing_word_count=len(forcing_polynomial),
        normalized_forcing_coefficient_l1=format_fraction(
            coefficient_l1(forcing_polynomial)
        ),
        lower_parallel_word_count=len(parallel),
        lower_parallel_coefficient_l1=format_fraction(
            coefficient_l1(parallel)
        ),
        lower_parallel_triple_delay_quotient_word_count=len(
            parallel_quotient
        ),
        normalized_lower_corner_word_count=len(lower_corner),
        exact_stein_residual_nonzero_entry_count=sum(
            int(entry != 0) for entry in stein_residual
        ),
        exact_upper_endpoint=str(endpoint),
        endpoint_null_transport_possible=False,
        all_checks_passed=verified,
    )


def write_records(
    records: list[SexticTransportTraceObstructionRecord],
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
            "repeated_crabb_sextic_transport_trace_obstruction_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the exact obstruction."""

    args = parse_args()
    records = [exact_record()]
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
