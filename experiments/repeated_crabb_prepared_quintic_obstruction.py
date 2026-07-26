#!/usr/bin/env python3
"""Disprove global odd-response parity after the first preparations.

The scalar-copy colligation from L308 has no nonzero endpoint response.
This checker performs the cubic cancellation, the complete quartic
lower neutralization, and the quintic lower neutralization exactly.
The remaining physical upper quintic is nevertheless nonzero.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from repeated_crabb_boundary_metric_flag import full_endpoint_gap_series
from repeated_crabb_canonical_quintic_preimage import (
    canonical_factor_series,
    canonical_metric_slack_and_operator,
    exact_canonical_factor_lifts,
)
from repeated_crabb_delayed_jet import series_multiply
from repeated_crabb_delayed_slack_anticommutator import (
    Polynomial,
    operator_series,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_lower_metric_flag import full_lower_endpoint_series
from repeated_crabb_moving_retightening_defect import exact_stein_inverse
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)
from repeated_crabb_raw_quintic_obstruction import (
    rational_colligation,
    retightening_frame_lift,
)
from repeated_crabb_transfer_channel_covariance import transfer_coefficient


Matrix = np.ndarray
ExactMatrix = sp.Matrix
DEGREE = 5
THETA = sp.Rational(1, 2)
EXPECTED_LOWER_FIVE = sp.Rational(3384, 3125)
EXPECTED_UPPER_FIVE = sp.Rational(66384, 15625)


@dataclass(frozen=True)
class PreparedQuinticObstructionRecord:
    """One exact and independently floating prepared-quintic audit."""

    state_dimension: int
    defect_dimension: int
    first_transfer: str
    second_transfer: str
    spectral_radius: str
    cubic_upper: str
    quartic_lower: str
    quartic_upper: str
    pre_neutral_quintic_lower: str
    post_neutral_quintic_lower: str
    post_neutral_quintic_upper_numerator: int
    post_neutral_quintic_upper_denominator: int
    floating_upper_quintic: str
    floating_exact_error: str
    response_dimension: int
    all_checks_passed: bool


def evaluate_polynomial(
    polynomial: Polynomial,
    operator: ExactMatrix,
) -> ExactMatrix:
    """Evaluate one exact word polynomial."""

    result = sp.zeros(operator.rows)
    cache: dict[str, ExactMatrix] = {"": sp.eye(operator.rows)}
    for word, coefficient in polynomial.items():
        if word not in cache:
            value = sp.eye(operator.rows)
            for letter in word:
                value *= operator if letter == "s" else operator.T
            cache[word] = value
        result += (
            sp.Rational(coefficient.numerator, coefficient.denominator)
            * cache[word]
        )
    return sp.simplify(result)


def exact_product(
    left: list[ExactMatrix],
    right: list[ExactMatrix],
) -> list[ExactMatrix]:
    """Multiply two truncated exact matrix series."""

    degree = min(len(left), len(right)) - 1
    return [
        sp.simplify(
            sum(
                (
                    left[index] * right[order - index]
                    for index in range(order + 1)
                ),
                sp.zeros(left[0].rows, right[0].cols),
            )
        )
        for order in range(degree + 1)
    ]


def exact_metric_from_factor(
    operator: ExactMatrix,
    operator_coefficients: list[ExactMatrix],
    factor_coefficients: list[ExactMatrix],
) -> list[ExactMatrix]:
    """Solve ``M-A.T M A=D D.T`` coefficient by coefficient."""

    degree = len(operator_coefficients) - 1
    slack = exact_product(
        factor_coefficients,
        [coefficient.T for coefficient in factor_coefficients],
    )
    metric = [sp.eye(operator.rows)]
    for order in range(1, degree + 1):
        forcing = slack[order]
        for left_degree in range(order + 1):
            for metric_degree in range(order):
                right_degree = order - left_degree - metric_degree
                if 0 <= right_degree <= degree:
                    forcing += (
                        operator_coefficients[left_degree].T
                        * metric[metric_degree]
                        * operator_coefficients[right_degree]
                    )
        metric.append(exact_stein_inverse(operator, forcing))
    return metric


def exact_correction_metric(
    operator: ExactMatrix,
    operator_coefficients: list[ExactMatrix],
    raw_factor: list[ExactMatrix],
    correction_factor: list[ExactMatrix],
) -> list[ExactMatrix]:
    """Solve the exact moving correction recurrence."""

    prepared_factor = [
        raw + correction
        for raw, correction in zip(
            raw_factor,
            correction_factor,
            strict=True,
        )
    ]
    raw_slack = exact_product(
        raw_factor,
        [coefficient.T for coefficient in raw_factor],
    )
    prepared_slack = exact_product(
        prepared_factor,
        [coefficient.T for coefficient in prepared_factor],
    )
    slack_change = [
        prepared - raw
        for prepared, raw in zip(
            prepared_slack,
            raw_slack,
            strict=True,
        )
    ]
    metric_change: list[ExactMatrix] = []
    degree = len(operator_coefficients) - 1
    for order in range(degree + 1):
        forcing = slack_change[order]
        for left_degree in range(order + 1):
            for metric_degree in range(order):
                right_degree = order - left_degree - metric_degree
                if 0 <= right_degree <= degree:
                    forcing += (
                        operator_coefficients[left_degree].T
                        * metric_change[metric_degree]
                        * operator_coefficients[right_degree]
                    )
        metric_change.append(exact_stein_inverse(operator, forcing))
    return metric_change


def exact_inverse_series(series: list[ExactMatrix]) -> list[ExactMatrix]:
    """Invert a square exact matrix series."""

    inverse = [series[0].inv()]
    for order in range(1, len(series)):
        convolution = sum(
            (
                series[index] * inverse[order - index]
                for index in range(1, order + 1)
            ),
            sp.zeros(series[0].rows),
        )
        inverse.append(sp.simplify(-inverse[0] * convolution))
    return inverse


def exact_schur_endpoint(
    series: list[ExactMatrix],
    endpoint_index: int,
) -> list[sp.Expr]:
    """Short a Hermitian series to one coordinate endpoint."""

    interior_indices = [
        index
        for index in range(series[0].rows)
        if index != endpoint_index
    ]
    endpoint = [
        coefficient[endpoint_index, endpoint_index]
        for coefficient in series
    ]
    cross = [
        coefficient.extract([endpoint_index], interior_indices)
        for coefficient in series
    ]
    interior = [
        coefficient.extract(interior_indices, interior_indices)
        for coefficient in series
    ]
    cross_square = exact_product(
        exact_product(cross, exact_inverse_series(interior)),
        [coefficient.T for coefficient in cross],
    )
    return [
        sp.simplify(value - cross_square[order][0, 0])
        for order, value in enumerate(endpoint)
    ]


def exact_endpoint_differences(
    metric: list[ExactMatrix],
    metric_change: list[ExactMatrix],
    right: ExactMatrix,
    left: ExactMatrix,
) -> tuple[list[sp.Expr], list[sp.Expr]]:
    """Return exact lower and physical-upper prepared-minus-raw gaps."""

    identity = sp.eye(metric[0].rows)
    initial = right * right.T
    final = left * left.T
    equality_metric = 2 * identity - initial + 2 * final
    upper_target = 4 * equality_metric.inv()
    lower_target = equality_metric.inv()
    prepared_metric = [
        base + change
        for base, change in zip(metric, metric_change, strict=True)
    ]
    raw_upper = [upper_target - metric[0]] + [
        -coefficient for coefficient in metric[1:]
    ]
    prepared_upper = [upper_target - prepared_metric[0]] + [
        -coefficient for coefficient in prepared_metric[1:]
    ]
    raw_lower = [metric[0] - lower_target] + metric[1:]
    prepared_lower = [
        prepared_metric[0] - lower_target
    ] + prepared_metric[1:]
    raw_upper_endpoint = exact_schur_endpoint(raw_upper, 2)
    prepared_upper_endpoint = exact_schur_endpoint(prepared_upper, 2)
    raw_lower_endpoint = exact_schur_endpoint(raw_lower, 0)
    prepared_lower_endpoint = exact_schur_endpoint(prepared_lower, 0)
    lower = [
        sp.simplify(prepared - raw)
        for prepared, raw in zip(
            prepared_lower_endpoint,
            raw_lower_endpoint,
            strict=True,
        )
    ]
    upper = [
        sp.simplify(4 * (prepared - raw))
        for prepared, raw in zip(
            prepared_upper_endpoint,
            raw_upper_endpoint,
            strict=True,
        )
    ]
    return lower, upper


def exact_prepared_endpoints() -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    sp.Expr,
    sp.Expr,
    sp.Expr,
]:
    """Construct the lower-neutral preparation exactly."""

    operator, right, left = rational_colligation()
    operator_coefficients = [
        evaluate_polynomial(coefficient, operator)
        for coefficient in operator_series(DEGREE)
    ]
    raw_factor = [
        evaluate_polynomial(coefficient, operator) * right
        for coefficient in exact_canonical_factor_lifts(DEGREE)
    ]
    metric = exact_metric_from_factor(
        operator,
        operator_coefficients,
        raw_factor,
    )
    first_frame = (
        evaluate_polynomial(retightening_frame_lift(1), operator)
        * right
    )
    second_frame = (
        evaluate_polynomial(retightening_frame_lift(2), operator)
        * right
    )
    correction_factor = [sp.zeros(3, 1) for _ in range(DEGREE + 1)]
    correction_factor[2] = THETA * first_frame
    correction_factor[4] = THETA * second_frame

    first_transfer = (left.T * operator.T * right)[0, 0]
    second_transfer = (left.T * (operator.T**2) * right)[0, 0]
    first_gram = first_transfer**2
    quartic_neutralizer = (
        (4 * THETA + THETA**2 / 2) * first_gram
        + (9 * THETA / 2 - 3 * THETA**2 / 4) * first_gram**2
    )
    correction_factor[4] += (
        sp.Rational(1, 2) * right * quartic_neutralizer
    )

    metric_change = exact_correction_metric(
        operator,
        operator_coefficients,
        raw_factor,
        correction_factor,
    )
    lower, upper = exact_endpoint_differences(
        metric,
        metric_change,
        right,
        left,
    )
    pre_neutral_lower_five = lower[5]
    correction_factor[5] = (
        -sp.Rational(1, 2) * right * pre_neutral_lower_five
    )
    metric_change = exact_correction_metric(
        operator,
        operator_coefficients,
        raw_factor,
        correction_factor,
    )
    lower, upper = exact_endpoint_differences(
        metric,
        metric_change,
        right,
        left,
    )
    return (
        lower,
        upper,
        first_transfer,
        second_transfer,
        pre_neutral_lower_five,
    )


def exact_response_values(
    operator: ExactMatrix,
    right: ExactMatrix,
    left: ExactMatrix,
) -> tuple[sp.Expr, ...]:
    """Evaluate the scalar response on a real basis of columns."""

    values: list[sp.Expr] = []
    for state_index in (1, 2):
        for phase in (sp.Integer(1), sp.I):
            column = sp.zeros(3, 1)
            column[state_index, 0] = phase
            forcing = right * column.conjugate().T + column * right.T
            response_metric = exact_stein_inverse(operator, forcing)
            values.append(
                sp.simplify((left.T * response_metric * left)[0, 0])
            )
    return tuple(values)


def numeric_correction_metric(
    operator: Matrix,
    operator_coefficients: list[Matrix],
    raw_factor: list[Matrix],
    correction_factor: list[Matrix],
) -> list[Matrix]:
    """Independently solve the floating moving correction recurrence."""

    prepared_factor = [
        raw + correction
        for raw, correction in zip(
            raw_factor,
            correction_factor,
            strict=True,
        )
    ]
    raw_slack = series_multiply(
        raw_factor,
        [coefficient.conj().T for coefficient in raw_factor],
        DEGREE,
    )
    prepared_slack = series_multiply(
        prepared_factor,
        [coefficient.conj().T for coefficient in prepared_factor],
        DEGREE,
    )
    metric_change: list[Matrix] = []
    for order in range(DEGREE + 1):
        forcing = prepared_slack[order] - raw_slack[order]
        for left_degree in range(order + 1):
            for metric_degree in range(order):
                right_degree = order - left_degree - metric_degree
                if 0 <= right_degree <= DEGREE:
                    forcing += (
                        operator_coefficients[left_degree].conj().T
                        @ metric_change[metric_degree]
                        @ operator_coefficients[right_degree]
                    )
        metric_change.append(stein_inverse(operator, forcing))
    return metric_change


def numeric_endpoint_differences(
    metric: list[Matrix],
    metric_change: list[Matrix],
    equality_metric: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[list[Matrix], list[Matrix]]:
    """Independently form the complete floating Schur endpoints."""

    prepared_metric = [
        base + change
        for base, change in zip(metric, metric_change, strict=True)
    ]
    upper_target = 4 * np.linalg.inv(equality_metric)
    lower_target = np.linalg.inv(equality_metric)
    raw_upper = [upper_target - metric[0]] + [
        -coefficient for coefficient in metric[1:]
    ]
    prepared_upper = [upper_target - prepared_metric[0]] + [
        -coefficient for coefficient in prepared_metric[1:]
    ]
    raw_lower = [metric[0] - lower_target] + metric[1:]
    prepared_lower = [
        prepared_metric[0] - lower_target
    ] + prepared_metric[1:]
    raw_upper_endpoint, _ = full_endpoint_gap_series(raw_upper, left)
    prepared_upper_endpoint, _ = full_endpoint_gap_series(
        prepared_upper,
        left,
    )
    raw_lower_endpoint, _ = full_lower_endpoint_series(raw_lower, right)
    prepared_lower_endpoint, _ = full_lower_endpoint_series(
        prepared_lower,
        right,
    )
    lower = [
        prepared - raw
        for prepared, raw in zip(
            prepared_lower_endpoint,
            raw_lower_endpoint,
            strict=True,
        )
    ]
    upper = [
        4 * (prepared - raw)
        for prepared, raw in zip(
            prepared_upper_endpoint,
            raw_upper_endpoint,
            strict=True,
        )
    ]
    return lower, upper


def floating_upper_quintic(
    exact_lower_five: sp.Expr,
) -> tuple[float, float]:
    """Independently reconstruct the final floating upper quintic."""

    exact_operator, exact_right, exact_left = rational_colligation()
    operator = np.array(exact_operator, dtype=complex)
    right = np.array(exact_right, dtype=complex)
    left = np.array(exact_left, dtype=complex)
    metric, slack, operator_coefficients, equality_metric = (
        canonical_metric_slack_and_operator(
            operator,
            right,
            left,
            DEGREE,
        )
    )
    raw_factor, _ = canonical_factor_series(slack, right)
    first = oriented_retightening_direction(
        operator,
        right,
        left,
        1,
    )
    second = oriented_retightening_direction(
        operator,
        right,
        left,
        2,
    )
    first_transfer = transfer_coefficient(
        operator,
        right,
        left,
        1,
    )
    first_gram = first_transfer.conj().T @ first_transfer
    theta = float(THETA)
    neutralizer = (
        (4 * theta + theta**2 / 2) * first_gram
        + (9 * theta / 2 - 3 * theta**2 / 4)
        * (first_gram @ first_gram)
    )
    correction_factor = [
        np.zeros_like(right) for _ in range(DEGREE + 1)
    ]
    correction_factor[2] = theta * first.frame
    correction_factor[4] = (
        theta * second.frame + 0.5 * right @ neutralizer
    )
    correction_factor[5] = (
        -0.5 * right * float(exact_lower_five)
    )
    metric_change = numeric_correction_metric(
        operator,
        operator_coefficients,
        raw_factor,
        correction_factor,
    )
    lower, upper = numeric_endpoint_differences(
        metric,
        metric_change,
        equality_metric,
        right,
        left,
    )
    return float(lower[5][0, 0].real), float(upper[5][0, 0].real)


def standard_record() -> PreparedQuinticObstructionRecord:
    """Run every exact and independent floating check."""

    (
        lower,
        upper,
        first_transfer,
        second_transfer,
        pre_neutral_lower,
    ) = exact_prepared_endpoints()
    floating_lower, floating_upper = floating_upper_quintic(
        pre_neutral_lower
    )
    floating_error = abs(floating_upper - float(EXPECTED_UPPER_FIVE))
    expected_quartic_lower = -THETA * second_transfer**2
    expected_quartic_upper = (
        4 * THETA * second_transfer**2
        - (8 * THETA + 4 * THETA**2) * first_transfer**2
        - (8 * THETA - 4 * THETA**2) * first_transfer**4
    )
    operator, right, left = rational_colligation()
    identity = sp.eye(3)
    initial_defect_error = (
        operator.T * operator - (identity - right * right.T)
    )
    final_defect_error = (
        operator * operator.T - (identity - left * left.T)
    )
    eigenvalues = operator.eigenvals()
    response_values = exact_response_values(operator, right, left)
    verified = bool(
        initial_defect_error.is_zero_matrix
        and final_defect_error.is_zero_matrix
        and first_transfer == sp.Rational(-4, 5)
        and second_transfer == sp.Rational(9, 25)
        and max(abs(complex(value)) for value in eigenvalues) < 1
        and upper[3] == 0
        and lower[4] == expected_quartic_lower
        and upper[4] == expected_quartic_upper
        and pre_neutral_lower == EXPECTED_LOWER_FIVE
        and lower[5] == 0
        and upper[5] == EXPECTED_UPPER_FIVE
        and all(value == 0 for value in response_values)
        and abs(floating_lower) < 1e-11
        and floating_error < 1e-11
        and EXPECTED_UPPER_FIVE != 0
    )
    if not verified:
        raise RuntimeError(
            "prepared quintic obstruction audit failed: "
            f"{first_transfer=}, {second_transfer=}, "
            f"{lower[4]=}, {upper[3]=}, {upper[4]=}, "
            f"{pre_neutral_lower=}, {lower[5]=}, {upper[5]=}, "
            f"{response_values=}, "
            f"{floating_lower=}, {floating_upper=}"
        )
    return PreparedQuinticObstructionRecord(
        state_dimension=3,
        defect_dimension=1,
        first_transfer=str(first_transfer),
        second_transfer=str(second_transfer),
        spectral_radius=str(sp.Rational(4, 5)),
        cubic_upper=str(upper[3]),
        quartic_lower=str(lower[4]),
        quartic_upper=str(upper[4]),
        pre_neutral_quintic_lower=str(pre_neutral_lower),
        post_neutral_quintic_lower=str(lower[5]),
        post_neutral_quintic_upper_numerator=(
            int(EXPECTED_UPPER_FIVE.p)
        ),
        post_neutral_quintic_upper_denominator=(
            int(EXPECTED_UPPER_FIVE.q)
        ),
        floating_upper_quintic=f"{floating_upper:.12e}",
        floating_exact_error=f"{floating_error:.12e}",
        response_dimension=int(any(value != 0 for value in response_values)),
        all_checks_passed=verified,
    )


def write_record(
    record: PreparedQuinticObstructionRecord,
    output: Path,
) -> str:
    """Write one deterministic record atomically and return its hash."""

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
            "repeated_crabb_prepared_quintic_obstruction_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the audit and persist the record."""

    args = parse_args()
    record = standard_record()
    digest = write_record(record, args.output)
    print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
