#!/usr/bin/env python3
"""Audit the polynomial preparation of the canonical quintic flag.

After the cubic and quartic defect-factor columns have been inserted,
the fifth upper endpoint is generally nonzero.  This checker verifies
an explicit bounded fifth column which moves that endpoint entirely
off

    ker(B_1*) intersect ker(B_2*).

The symbolic part is a finite exact rational-word certificate.  The
numerical part independently reconstructs the prepared positive Stein
factor and both physical Schur endpoints through degree five.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import (
    full_endpoint_gap_series,
    rank_chain_case,
)
from repeated_crabb_boundary_slack_deflation import physical_metric_root
from repeated_crabb_canonical_cubic_preimage import (
    coboundary_witness as cubic_coboundary_witness,
    polynomial_column as cubic_column,
)
from repeated_crabb_canonical_quartic_preimage import (
    evaluate_polynomial,
    quartic_coboundary_witness,
    quartic_column,
    quartic_column_lift,
)
from repeated_crabb_canonical_quartic_trace import raw_slack_series
from repeated_crabb_canonical_repair_flag_obstruction import (
    canonical_repair_metric_series,
)
from repeated_crabb_delayed_jet import (
    ellipse_operator_coefficients,
    series_multiply as matrix_series_multiply,
)
from repeated_crabb_delayed_slack_anticommutator import (
    E,
    F,
    IDENTITY,
    Polynomial,
    S,
    STAR,
    add,
    adjoint,
    boundary_metric_series,
    multiply,
    operator_series,
    scale,
    schur_residual_series,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_lower_metric_flag import full_lower_endpoint_series
from repeated_crabb_transfer_channel_covariance import (
    inflated_case,
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
DEGREE = 5


@dataclass(frozen=True)
class CanonicalQuinticPreimageRecord:
    """One exact-structure and numerical quintic preparation audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    colligation_error: str
    spectral_radius: str
    first_transfer_norm: str
    second_transfer_norm: str
    quintic_endpoint_norm: str
    predicted_endpoint_norm: str
    endpoint_reconstruction_error: str
    left_flag_dimension: int
    left_flag_compression_norm: str
    polynomial_column_norm: str
    perpendicular_column_error: str
    lower_endpoint_motion_error: str
    canonical_factor_error: str
    exact_coboundary_residual_word_count: int
    exact_flag_factor_residual_word_count: int
    exact_witness_hermitian_residual_word_count: int
    exact_column_perpendicular_residual_word_count: int
    all_checks_passed: bool


def cubic_column_lift() -> Polynomial:
    """Return the exact state lift ``C_3 V*``."""

    first_right_gram = multiply(
        multiply(
            multiply(
                multiply(E, S),
                F,
            ),
            STAR,
        ),
        E,
    )
    return scale(
        3,
        multiply({"aa": 1}, first_right_gram),
    )


def quintic_column_bracket() -> Polynomial:
    """Return the twenty-term polynomial inside the left projection."""

    terms: tuple[tuple[Fraction | int, str], ...] = (
        (Fraction(1, 2), "ssssaa"),
        (Fraction(3, 2), "aaasssaa"),
        (Fraction(3, 2), "aasssaaa"),
        (Fraction(-1, 2), "sssssaaa"),
        (3, "aaaassssaa"),
        (2, "aaassssaaa"),
        (Fraction(5, 2), "aaaaaaasssaa"),
        (-3, "aaaaaasssaaa"),
        (-1, "aaaasssaaaaa"),
        (-1, "aaasssaaaaaa"),
        (Fraction(3, 2), "aaasssssssaa"),
        (Fraction(-1, 2), "aasssaaaaaaa"),
        (Fraction(1, 2), "aasssssssaaa"),
        (Fraction(1, 2), "saaassssssaa"),
        (Fraction(-3, 2), "ssaaaaasssaa"),
        (Fraction(3, 2), "ssaaaasssaaa"),
        (Fraction(1, 2), "ssssaaasssaa"),
        (17, "aaasssaaasssaa"),
        (Fraction(-1, 2), "aasssaaaasssaa"),
        (Fraction(15, 2), "aasssaaasssaaa"),
    )
    result: Polynomial = {}
    for coefficient, word in terms:
        result = add(result, scale(coefficient, {word: 1}))
    return result


def quintic_column_lift() -> Polynomial:
    """Return the exact perpendicular state lift ``C_5 V*``."""

    initial_complement = add(IDENTITY, scale(-1, E))
    return multiply(
        multiply(
            initial_complement,
            quintic_column_bracket(),
        ),
        E,
    )


def quintic_column(
    partial: Matrix,
    right: Matrix,
) -> Matrix:
    """Evaluate the explicit fifth defect-factor column."""

    return evaluate_polynomial(quintic_column_lift(), partial) @ right


def quintic_coboundary_witness() -> Polynomial:
    """Return the Hermitian 55-term fifth Stein witness."""

    terms: tuple[tuple[Fraction | int, str, str], ...] = (
        (-4, "aa", "ss"),
        (22, "aaas", "asss"),
        (22, "saaa", "sssa"),
        (2, "aaaaaa", "ssssss"),
        (-1, "aaaass", "aassss"),
        (-18, "assssa", "saaaas"),
        (-4, "ssaaaa", "ssssaa"),
        (-1, "aaaaaaas", "asssssss"),
        (Fraction(13, 2), "aaaasssa", "saaassss"),
        (Fraction(-13, 2), "aaasssaa", "ssaaasss"),
        (Fraction(-29, 2), "aasssaaa", "sssaaass"),
        (Fraction(5, 2), "aasssssa", "saaaaass"),
        (2, "asssaaaa", "ssssaaas"),
        (Fraction(5, 2), "asssssaa", "ssaaaaas"),
        (-1, "saaaaaaa", "sssssssa"),
        (-1, "aaaassssaa", "ssaaaassss"),
        (5, "aaasssaaas", "asssaaasss"),
        (Fraction(1, 2), "aaassssssa", "saaaaaasss"),
        (Fraction(25, 2), "aasssaaaas", "assssaaass"),
        (2, "aassssssaa", "ssaaaaaass"),
        (Fraction(-3, 2), "asssaaaaas", "asssssaaas"),
        (Fraction(1, 2), "assssssaaa", "sssaaaaaas"),
        (Fraction(1, 2), "assssssssa", "saaaaaaaas"),
        (-7, "saaaaasssa", "saaasssssa"),
        (8, "saaaasssaa", "ssaaassssa"),
        (13, "saaasssaaa", "sssaaasssa"),
        (Fraction(1, 2), "aaaaaaasssaa", "ssaaasssssss"),
        (-1, "aaaaasssaaaa", "ssssaaasssss"),
        (Fraction(-7, 2), "aaaasssaaass", "aasssaaassss"),
        (1, "aaasssaaaass", "aassssaaasss"),
        (Fraction(-1, 2), "aaasssssssaa", "ssaaaaaaasss"),
        (Fraction(1, 2), "aasssaaaaaaa", "sssssssaaass"),
        (Fraction(-5, 2), "aasssaaaaass", "aasssssaaass"),
        (Fraction(-5, 2), "aassssssaaas", "asssaaaaaass"),
        (Fraction(-1, 2), "aasssssssaaa", "sssaaaaaaass"),
        (-1, "asssaaassssa", "saaaasssaaas"),
        (Fraction(-13, 2), "assssaaasssa", "saaasssaaaas"),
        (Fraction(-5, 2), "saaassssssaa", "ssaaaaaasssa"),
        (-8, "ssaaaaasssaa", "ssaaasssssaa"),
        (4, "ssaaaasssaaa", "sssaaassssaa"),
        (-1, "aaaasssaaasssa", "saaasssaaassss"),
        (Fraction(1, 2), "aaasssaaaasssa", "saaassssaaasss"),
        (Fraction(-15, 2), "aaasssaaasssaa", "ssaaasssaaasss"),
        (Fraction(5, 2), "aasssaaaaasssa", "saaasssssaaass"),
        (Fraction(-5, 2), "aasssaaaasssaa", "ssaaassssaaass"),
        (Fraction(-15, 2), "aasssaaasssaaa", "sssaaasssaaass"),
        (5, "aasssaaasssssa", "saaaaasssaaass"),
        (-3, "aassssaaassssa", "saaaasssaaaass"),
        (Fraction(-1, 2), "aasssssaaasssa", "saaasssaaaaass"),
        (2, "asssaaaaaasssa", "saaassssssaaas"),
        (Fraction(5, 2), "asssaaaaasssaa", "ssaaasssssaaas"),
        (Fraction(-1, 2), "asssaaaasssaaa", "sssaaassssaaas"),
        (5, "asssaaasssssaa", "ssaaaaasssaaas"),
        (-3, "assssaaassssaa", "ssaaaasssaaaas"),
        (Fraction(-1, 2), "asssssaaasssaa", "ssaaasssaaaaas"),
    )
    result: Polynomial = {}
    for coefficient, word, adjoint_word in terms:
        result = add(
            result,
            scale(
                coefficient,
                {word: 1, adjoint_word: 1},
            ),
        )
    return result


def exact_canonical_factor_lifts() -> list[Polynomial]:
    """Return canonical factor coefficients ``D_j V*`` through degree five."""

    raw = raw_slack_series(DEGREE)
    residual = schur_residual_series(DEGREE)
    slack = [
        add(raw_item, scale(-1, residual_item))
        for raw_item, residual_item in zip(raw, residual, strict=True)
    ]
    initial_complement = multiply(STAR, S)
    factors: list[Polynomial] = [E]
    for order in range(1, DEGREE + 1):
        remainder = slack[order]
        for left_degree in range(1, order):
            remainder = add(
                remainder,
                scale(
                    -1,
                    multiply(
                        factors[left_degree],
                        adjoint(factors[order - left_degree]),
                    ),
                ),
            )
        factor = add(
            multiply(multiply(initial_complement, remainder), E),
            scale(
                Fraction(1, 2),
                multiply(multiply(E, remainder), E),
            ),
        )
        factors.append(factor)

        reconstruction: Polynomial = {}
        for left_degree in range(order + 1):
            reconstruction = add(
                reconstruction,
                multiply(
                    factors[left_degree],
                    adjoint(factors[order - left_degree]),
                ),
            )
        if add(slack[order], scale(-1, reconstruction)):
            raise RuntimeError(
                f"exact canonical factor failed at order {order}"
            )
    return factors


@lru_cache(maxsize=1)
def exact_fifth_components() -> tuple[Polynomial, Polynomial]:
    """Return the fifth Stein forcing and upper Schur cross term."""

    factors = exact_canonical_factor_lifts()
    factors[3] = add(factors[3], cubic_column_lift())
    factors[4] = add(factors[4], quartic_column_lift())
    prepared_slack: list[Polynomial] = []
    for order in range(DEGREE + 1):
        coefficient: Polynomial = {}
        for left_degree in range(order + 1):
            coefficient = add(
                coefficient,
                multiply(
                    factors[left_degree],
                    adjoint(factors[order - left_degree]),
                ),
            )
        prepared_slack.append(coefficient)

    operator = operator_series(DEGREE)
    boundary = boundary_metric_series(DEGREE)
    residual = schur_residual_series(DEGREE)
    metric = [{} for _ in range(DEGREE + 1)]
    metric[0] = IDENTITY
    metric[2] = add(
        boundary[2],
        scale(-1, residual[2]),
    )
    metric[3] = cubic_coboundary_witness()
    metric[4] = quartic_coboundary_witness()

    forcing = prepared_slack[DEGREE]
    for left_degree in range(DEGREE + 1):
        for metric_degree in range(DEGREE):
            right_degree = DEGREE - left_degree - metric_degree
            if 0 <= right_degree <= DEGREE:
                forcing = add(
                    forcing,
                    multiply(
                        multiply(
                            adjoint(operator[left_degree]),
                            metric[metric_degree],
                        ),
                        operator[right_degree],
                    ),
                )

    upper_complement_inverse = add(
        IDENTITY,
        scale(-1, F),
        scale(Fraction(-2, 3), E),
    )
    schur_cross = add(
        multiply(
            multiply(
                multiply(metric[2], upper_complement_inverse),
                metric[3],
            ),
            F,
        ),
        multiply(
            multiply(
                multiply(metric[3], upper_complement_inverse),
                metric[2],
            ),
            F,
        ),
    )
    return forcing, schur_cross


def fifth_endpoint_lift() -> Polynomial:
    """Return one minus-quarter of the physical fifth endpoint."""

    _, schur_cross = exact_fifth_components()
    witness = quintic_coboundary_witness()
    return add(
        multiply(multiply(F, witness), F),
        multiply(F, schur_cross),
    )


def fifth_flag_factor_lift() -> Polynomial:
    """Factor the endpoint through the first two delayed channels."""

    first_left = multiply(
        F,
        add(
            scale(Fraction(7, 2), {"aaaasssaa": 1}),
            scale(Fraction(5, 2), {"aasssaaaa": 1}),
        ),
    )
    second_left = multiply(F, {"aassssaa": 1})
    first_channel = multiply(multiply(E, S), F)
    second_channel = multiply(multiply(E, {"ss": 1}), F)
    first_term = multiply(first_left, first_channel)
    second_term = multiply(second_left, second_channel)
    return add(
        first_term,
        adjoint(first_term),
        second_term,
        adjoint(second_term),
    )


@lru_cache(maxsize=1)
def exact_residuals() -> tuple[int, int, int, int]:
    """Return exact word counts for all certificate residuals."""

    forcing, _ = exact_fifth_components()
    column = quintic_column_lift()
    witness = quintic_coboundary_witness()
    coboundary = add(
        witness,
        scale(
            -1,
            multiply(multiply(STAR, witness), S),
        ),
    )
    forcing_residual = add(
        forcing,
        column,
        adjoint(column),
        scale(-1, coboundary),
    )
    flag_residual = add(
        fifth_endpoint_lift(),
        scale(-1, fifth_flag_factor_lift()),
    )
    hermitian_residual = add(
        witness,
        scale(-1, adjoint(witness)),
    )
    perpendicular_residual = multiply(E, column)
    return (
        len(forcing_residual),
        len(flag_residual),
        len(hermitian_residual),
        len(perpendicular_residual),
    )


def matrix_adjoint_series(series: list[Matrix]) -> list[Matrix]:
    """Return the coefficientwise adjoint of a matrix series."""

    return [coefficient.conj().T for coefficient in series]


def canonical_metric_slack_and_operator(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    degree: int,
) -> tuple[list[Matrix], list[Matrix], list[Matrix], Matrix]:
    """Return the canonical metric, its slack, and balanced operator."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    root, inverse_root = physical_metric_root(partial, right, left)
    physical_operator = inverse_root @ partial @ root
    physical_series = ellipse_operator_coefficients(
        physical_operator,
        degree,
    )
    operator = [
        root @ coefficient @ inverse_root
        for coefficient in physical_series
    ]
    upper_gap, _ = canonical_repair_metric_series(
        partial,
        right,
        left,
        degree,
    )
    equality_metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    upper_target = 4 * np.linalg.inv(equality_metric)
    metric = [upper_target - upper_gap[0]] + [
        -coefficient for coefficient in upper_gap[1:]
    ]
    pulled = matrix_series_multiply(
        matrix_series_multiply(
            matrix_adjoint_series(operator),
            metric,
            degree,
        ),
        operator,
        degree,
    )
    slack = [
        coefficient - image
        for coefficient, image in zip(metric, pulled, strict=True)
    ]
    return metric, slack, operator, equality_metric


def canonical_factor_series(
    slack: list[Matrix],
    right: Matrix,
) -> tuple[list[Matrix], float]:
    """Factor a rank-defect positive series in the Hermitian gauge."""

    projection = right @ right.conj().T
    complement = np.eye(len(projection)) - projection
    factor = [right]
    maximum_error = 0.0
    for order in range(1, len(slack)):
        remainder = slack[order].copy()
        for left_degree in range(1, order):
            remainder -= (
                factor[left_degree]
                @ factor[order - left_degree].conj().T
            )
        coefficient = (
            complement @ remainder @ right
            + 0.5
            * right
            @ (right.conj().T @ remainder @ right)
        )
        factor.append(coefficient)
        reconstruction = np.zeros_like(remainder)
        for left_degree in range(order + 1):
            reconstruction += (
                factor[left_degree]
                @ factor[order - left_degree].conj().T
            )
        maximum_error = max(
            maximum_error,
            float(np.linalg.norm(slack[order] - reconstruction)),
        )
    return factor, maximum_error


def prepared_endpoint_series(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    degree: int = DEGREE,
) -> tuple[list[Matrix], list[Matrix], float]:
    """Return both endpoints after all three polynomial preparations."""

    metric, slack, operator, equality_metric = (
        canonical_metric_slack_and_operator(
            partial,
            right,
            left,
            degree,
        )
    )
    canonical_factor, factor_error = canonical_factor_series(slack, right)
    prepared_factor = [
        coefficient.copy() for coefficient in canonical_factor
    ]
    prepared_factor[3] += cubic_column(partial, right, left)
    prepared_factor[4] += quartic_column(partial, right)
    prepared_factor[5] += quintic_column(partial, right)
    prepared_slack = matrix_series_multiply(
        prepared_factor,
        matrix_adjoint_series(prepared_factor),
        degree,
    )
    slack_change = [
        prepared - canonical
        for prepared, canonical in zip(
            prepared_slack,
            slack,
            strict=True,
        )
    ]

    metric_change: list[Matrix] = []
    for order in range(degree + 1):
        forcing = slack_change[order].copy()
        for left_degree in range(order + 1):
            for metric_degree in range(order):
                right_degree = order - left_degree - metric_degree
                if 0 <= right_degree <= degree:
                    forcing += (
                        operator[left_degree].conj().T
                        @ metric_change[metric_degree]
                        @ operator[right_degree]
                    )
        metric_change.append(stein_inverse(partial, forcing))

    prepared_metric = [
        base + change
        for base, change in zip(metric, metric_change, strict=True)
    ]
    upper_target = 4 * np.linalg.inv(equality_metric)
    upper_gap = [upper_target - prepared_metric[0]] + [
        -coefficient for coefficient in prepared_metric[1:]
    ]
    lower_target = np.linalg.inv(equality_metric)
    lower_gap = [prepared_metric[0] - lower_target] + [
        coefficient for coefficient in prepared_metric[1:]
    ]
    upper_endpoint, upper_error = full_endpoint_gap_series(
        upper_gap,
        left,
    )
    lower_endpoint, lower_error = full_lower_endpoint_series(
        lower_gap,
        right,
    )
    if max(upper_error, lower_error) > 3e-8:
        raise RuntimeError("the prepared endpoint base was inconsistent")
    return (
        [4 * coefficient for coefficient in upper_endpoint],
        lower_endpoint,
        factor_error,
    )


def joint_left_kernel(*matrices: Matrix) -> Matrix:
    """Return an orthonormal frame for the joint kernels of ``B_j*``."""

    stacked = np.vstack([matrix.conj().T for matrix in matrices])
    _, singular_values, adjoint_frame = np.linalg.svd(
        stacked,
        full_matrices=True,
    )
    tolerance = max(
        1e-10,
        1e-8 * singular_values[0] if len(singular_values) else 0,
    )
    rank = int(np.count_nonzero(singular_values > tolerance))
    return adjoint_frame.conj().T[:, rank:]


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
) -> CanonicalQuinticPreimageRecord:
    """Audit one exact quintic polynomial preparation."""

    upper, _, factor_error = prepared_endpoint_series(
        partial,
        right,
        left,
    )
    first = transfer_coefficient(partial, right, left, 1)
    second = transfer_coefficient(partial, right, left, 2)
    flag = joint_left_kernel(first, second)
    flag_compression = flag.conj().T @ upper[5] @ flag

    endpoint_lift = evaluate_polynomial(
        fifth_endpoint_lift(),
        partial,
    )
    predicted_endpoint = (
        -4 * left.conj().T @ endpoint_lift @ left
    )
    endpoint_error = float(
        np.linalg.norm(upper[5] - predicted_endpoint)
    )

    column = quintic_column(partial, right)
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ column)
    )
    column_forcing = (
        right @ column.conj().T
        + column @ right.conj().T
    )
    column_response = stein_inverse(partial, column_forcing)
    lower_motion_error = float(
        np.linalg.norm(
            right.conj().T @ column_response @ right
        )
    )

    identity = np.eye(len(partial), dtype=complex)
    colligation_error = max(
        float(
            np.linalg.norm(
                partial.conj().T @ partial
                - (identity - right @ right.conj().T)
            )
        ),
        float(
            np.linalg.norm(
                partial @ partial.conj().T
                - (identity - left @ left.conj().T)
            )
        ),
        float(np.linalg.norm(right.conj().T @ left)),
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    flag_error = float(np.linalg.norm(flag_compression))
    exact_forcing, exact_flag, exact_hermitian, exact_perpendicular = (
        exact_residuals()
    )

    tolerance = 3e-8
    verified = bool(
        colligation_error < tolerance
        and spectral_radius < 1
        and endpoint_error < tolerance
        and flag_error < tolerance
        and perpendicular_error < tolerance
        and lower_motion_error < tolerance
        and factor_error < tolerance
        and exact_forcing == 0
        and exact_flag == 0
        and exact_hermitian == 0
        and exact_perpendicular == 0
    )
    if not verified:
        raise RuntimeError(
            "the canonical quintic preimage audit failed: "
            f"kind={construction_kind}, "
            f"endpoint={endpoint_error:.3e}, "
            f"flag={flag_error:.3e}, "
            f"perpendicular={perpendicular_error:.3e}, "
            f"lower={lower_motion_error:.3e}, "
            f"factor={factor_error:.3e}, "
            f"exact=({exact_forcing},{exact_flag},"
            f"{exact_hermitian},{exact_perpendicular})"
        )

    return CanonicalQuinticPreimageRecord(
        construction_kind=construction_kind,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        parameter_scale=format_float(parameter_scale),
        colligation_error=format_float(colligation_error),
        spectral_radius=format_float(spectral_radius),
        first_transfer_norm=format_float(float(np.linalg.norm(first))),
        second_transfer_norm=format_float(float(np.linalg.norm(second))),
        quintic_endpoint_norm=format_float(
            float(np.linalg.norm(upper[5]))
        ),
        predicted_endpoint_norm=format_float(
            float(np.linalg.norm(predicted_endpoint))
        ),
        endpoint_reconstruction_error=format_float(endpoint_error),
        left_flag_dimension=flag.shape[1],
        left_flag_compression_norm=format_float(flag_error),
        polynomial_column_norm=format_float(float(np.linalg.norm(column))),
        perpendicular_column_error=format_float(perpendicular_error),
        lower_endpoint_motion_error=format_float(lower_motion_error),
        canonical_factor_error=format_float(factor_error),
        exact_coboundary_residual_word_count=exact_forcing,
        exact_flag_factor_residual_word_count=exact_flag,
        exact_witness_hermitian_residual_word_count=exact_hermitian,
        exact_column_perpendicular_residual_word_count=exact_perpendicular,
        all_checks_passed=verified,
    )


def standard_records() -> list[CanonicalQuinticPreimageRecord]:
    """Return unstructured, rank-changing, and delayed records."""

    records: list[CanonicalQuinticPreimageRecord] = []
    for defect_dimension in range(1, 5):
        for repetition in range(2):
            state_dimension = 4 * defect_dimension + 3 + repetition
            partial, right, left = random_partial_isometry(
                state_dimension,
                defect_dimension,
                np.random.default_rng(
                    116_000 + 100 * defect_dimension + repetition
                ),
            )
            records.append(
                audit_case(
                    "unstructured",
                    partial,
                    right,
                    left,
                    1,
                )
            )

    for defect_dimension in range(3, 7):
        for parameter_scale in (1.0, 0.5, 0.2, 0.05):
            partial, right, left, _, _ = rank_chain_case(
                defect_dimension,
                109_200 + defect_dimension,
                parameter_scale,
            )
            records.append(
                audit_case(
                    "rank_chain",
                    partial,
                    right,
                    left,
                    parameter_scale,
                )
            )

    for defect_dimension in (2, 3, 4):
        partial, right, left, _ = inflated_case(
            4 * defect_dimension + 3,
            defect_dimension,
            3,
            defect_dimension,
            116_800 + defect_dimension,
        )
        records.append(
            audit_case(
                "complete_double_delay",
                partial,
                right,
                left,
                1,
            )
        )
    return records


def write_records(
    records: list[CanonicalQuinticPreimageRecord],
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
            "repeated_crabb_canonical_quintic_preimage_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the complete audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
