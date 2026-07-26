#!/usr/bin/env python3
"""Certify the nonzero raw quintic copy trace exactly."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from repeated_crabb_canonical_cubic_selection import (
    cyclic_trace_classes,
)
from repeated_crabb_canonical_quintic_preimage import (
    canonical_factor_series,
    canonical_metric_slack_and_operator,
    exact_canonical_factor_lifts,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    Q,
    S,
    STAR,
    add,
    adjoint,
    multiply,
    operator_series,
    scale,
)
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)


Matrix = np.ndarray
EXPECTED_TRACE = Fraction(33_264, 15_625)


@dataclass(frozen=True)
class RawQuinticObstructionRecord:
    """Exact and independent numerical audits of the raw quintic."""

    state_dimension: int
    defect_dimension: int
    first_transfer: str
    spectral_radius: str
    exact_cubic_class_count: int
    exact_cyclic_class_count: int
    exact_quintic_trace_numerator: int
    exact_quintic_trace_denominator: int
    direct_quintic_trace: str
    direct_trace_error: str
    cubic_trace: str
    partial_isometry_error: str
    all_checks_passed: bool


def polynomial_power(polynomial: Polynomial, exponent: int) -> Polynomial:
    """Return a nonnegative word-polynomial power."""

    result = IDENTITY
    for _ in range(exponent):
        result = multiply(result, polynomial)
    return result


def retightening_frame_lift(grade: int) -> Polynomial:
    """Return the exact state lift ``F_grade V*`` from L298."""

    forward = polynomial_power(S, grade)
    result = multiply(
        multiply(
            multiply(forward, F),
            polynomial_power(STAR, grade),
        ),
        E,
    )
    for earlier in range(1, grade):
        term = multiply(
            polynomial_power(STAR, grade - earlier),
            E,
        )
        term = multiply(term, forward)
        term = multiply(term, F)
        term = multiply(
            term,
            polynomial_power(STAR, earlier),
        )
        term = multiply(term, E)
        result = add(result, multiply(Q, term))
    return scale(Fraction(-1, 2), result)


def finite_stein_sum(
    forcing: Polynomial,
    maximum_steps: int = 12,
) -> Polynomial:
    """Return a terminating exact Stein sum."""

    result: Polynomial = {}
    term = forcing
    for _ in range(maximum_steps):
        result = add(result, term)
        term = multiply(multiply(STAR, term), S)
        if not term:
            return result
    raise RuntimeError("the exact grade-one Stein sum did not terminate")


def exact_raw_quintic_trace_lift() -> Polynomial:
    """Return a finite trace-equivalent lift of the raw quintic."""

    operators = operator_series(5)
    factors = exact_canonical_factor_lifts(5)
    frame_one = retightening_frame_lift(1)
    frame_two = retightening_frame_lift(2)
    metric_one = finite_stein_sum(
        add(frame_one, adjoint(frame_one))
    )

    dual_solution = scale(
        2,
        add(multiply(S, S), multiply(STAR, STAR)),
    )
    grade_two_forcing = add(frame_two, adjoint(frame_two))
    terms = [
        # Trace of the A_1--R_2 cross after L302's dual telescope.
        multiply(grade_two_forcing, dual_solution),
        multiply(
            multiply(adjoint(operators[3]), metric_one),
            S,
        ),
        multiply(multiply(STAR, metric_one), operators[3]),
        multiply(
            multiply(adjoint(operators[1]), metric_one),
            operators[2],
        ),
        multiply(
            multiply(adjoint(operators[2]), metric_one),
            operators[1],
        ),
        multiply(factors[1], adjoint(frame_two)),
        multiply(frame_two, adjoint(factors[1])),
        multiply(factors[3], adjoint(frame_one)),
        multiply(frame_one, adjoint(factors[3])),
    ]
    return scale(-1, add(*terms))


def exact_raw_cubic_trace_lift() -> Polynomial:
    """Return the exact raw cubic state coefficient."""

    operators = operator_series(3)
    factors = exact_canonical_factor_lifts(3)
    frame_one = retightening_frame_lift(1)
    metric_one = finite_stein_sum(
        add(frame_one, adjoint(frame_one))
    )
    return scale(
        -1,
        add(
            multiply(
                multiply(adjoint(operators[1]), metric_one),
                S,
            ),
            multiply(multiply(STAR, metric_one), operators[1]),
            multiply(factors[1], adjoint(frame_one)),
            multiply(frame_one, adjoint(factors[1])),
        ),
    )


def rational_colligation() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return the stable scalar-copy rational partial isometry."""

    operator = sp.Matrix(
        [
            [0, sp.Rational(3, 5), sp.Rational(-4, 5)],
            [0, sp.Rational(4, 5), sp.Rational(3, 5)],
            [0, 0, 0],
        ]
    )
    right = sp.eye(3)[:, :1]
    left = sp.eye(3)[:, 2:]
    return operator, right, left


def evaluate_polynomial(
    polynomial: Polynomial,
    operator: sp.Matrix,
) -> sp.Matrix:
    """Evaluate an exact word polynomial."""

    dimension = operator.rows
    result = sp.zeros(dimension)
    cache: dict[str, sp.Matrix] = {"": sp.eye(dimension)}
    for word, coefficient in polynomial.items():
        if word not in cache:
            value = sp.eye(dimension)
            for letter in word:
                value *= operator if letter == "s" else operator.T
            cache[word] = value
        result += (
            sp.Rational(coefficient.numerator, coefficient.denominator)
            * cache[word]
        )
    return result


def direct_coefficients(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[float, float, float]:
    """Return independent cubic/quintic traces and factor error."""

    _, slack, operator_series_value, _ = (
        canonical_metric_slack_and_operator(
            operator,
            right,
            left,
            5,
        )
    )
    factors, factor_error = canonical_factor_series(slack, right)
    direction_one = oriented_retightening_direction(
        operator,
        right,
        left,
        1,
    )
    direction_two = oriented_retightening_direction(
        operator,
        right,
        left,
        2,
    )

    cubic = -(
        operator_series_value[1].conj().T
        @ direction_one.metric
        @ operator
        + operator.conj().T
        @ direction_one.metric
        @ operator_series_value[1]
        + factors[1] @ direction_one.frame.conj().T
        + direction_one.frame @ factors[1].conj().T
    )
    quintic = -(
        operator_series_value[1].conj().T
        @ direction_two.metric
        @ operator
        + operator.conj().T
        @ direction_two.metric
        @ operator_series_value[1]
        + operator_series_value[3].conj().T
        @ direction_one.metric
        @ operator
        + operator.conj().T
        @ direction_one.metric
        @ operator_series_value[3]
        + operator_series_value[1].conj().T
        @ direction_one.metric
        @ operator_series_value[2]
        + operator_series_value[2].conj().T
        @ direction_one.metric
        @ operator_series_value[1]
        + factors[1] @ direction_two.frame.conj().T
        + direction_two.frame @ factors[1].conj().T
        + factors[3] @ direction_one.frame.conj().T
        + direction_one.frame @ factors[3].conj().T
    )
    return (
        float(np.trace(cubic).real),
        float(np.trace(quintic).real),
        factor_error,
    )


def standard_record() -> RawQuinticObstructionRecord:
    """Construct the exact obstruction and all independent checks."""

    exact_operator, exact_right, exact_left = rational_colligation()
    identity = sp.eye(3)
    partial_error_exact = (
        exact_operator.T * exact_operator
        - (identity - exact_right * exact_right.T)
    )
    final_error_exact = (
        exact_operator * exact_operator.T
        - (identity - exact_left * exact_left.T)
    )
    first_transfer = (
        exact_left.T * exact_operator.T * exact_right
    )[0, 0]

    trace_lift = exact_raw_quintic_trace_lift()
    cubic_classes = cyclic_trace_classes(
        exact_raw_cubic_trace_lift()
    )
    cyclic_classes = cyclic_trace_classes(trace_lift)
    exact_trace = sp.trace(
        evaluate_polynomial(trace_lift, exact_operator)
    )

    operator = np.array(exact_operator, dtype=complex)
    right = np.array(exact_right, dtype=complex)
    left = np.array(exact_left, dtype=complex)
    cubic_trace, quintic_trace, factor_error = direct_coefficients(
        operator,
        right,
        left,
    )
    expected_float = float(EXPECTED_TRACE)
    trace_error = abs(quintic_trace - expected_float)
    partial_error = float(
        np.linalg.norm(np.array(partial_error_exact, dtype=float))
        + np.linalg.norm(np.array(final_error_exact, dtype=float))
    )
    spectral_radius = max(
        abs(complex(value))
        for value in exact_operator.eigenvals()
    )

    verified = bool(
        partial_error_exact.is_zero_matrix
        and final_error_exact.is_zero_matrix
        and first_transfer == sp.Rational(-4, 5)
        and exact_trace
        == sp.Rational(
            EXPECTED_TRACE.numerator,
            EXPECTED_TRACE.denominator,
        )
        and not cubic_classes
        and len(cyclic_classes) == 8
        and abs(cubic_trace) < 1e-11
        and trace_error < 1e-11
        and factor_error < 1e-11
        and spectral_radius < 1
    )
    if not verified:
        raise RuntimeError(
            "raw quintic obstruction audit failed: "
            f"{first_transfer=}, {exact_trace=}, "
            f"{len(cubic_classes)=}, "
            f"{len(cyclic_classes)=}, {cubic_trace=}, "
            f"{quintic_trace=}, {trace_error=}, "
            f"{factor_error=}, {spectral_radius=}"
        )

    return RawQuinticObstructionRecord(
        state_dimension=3,
        defect_dimension=1,
        first_transfer=str(first_transfer),
        spectral_radius=str(sp.Rational(4, 5)),
        exact_cubic_class_count=len(cubic_classes),
        exact_cyclic_class_count=len(cyclic_classes),
        exact_quintic_trace_numerator=EXPECTED_TRACE.numerator,
        exact_quintic_trace_denominator=EXPECTED_TRACE.denominator,
        direct_quintic_trace=f"{quintic_trace:.12e}",
        direct_trace_error=f"{trace_error:.12e}",
        cubic_trace=f"{cubic_trace:.12e}",
        partial_isometry_error=f"{partial_error:.12e}",
        all_checks_passed=verified,
    )


def write_record(
    record: RawQuinticObstructionRecord,
    output: Path,
) -> str:
    """Write the deterministic record atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(asdict(record), sort_keys=True) + "\n",
        encoding="utf-8",
    )
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
            "repeated_crabb_raw_quintic_obstruction_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the raw-quintic obstruction audit."""

    args = parse_args()
    record = standard_record()
    digest = write_record(record, args.output)
    print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
