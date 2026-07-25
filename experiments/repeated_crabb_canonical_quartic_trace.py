#!/usr/bin/env python3
"""Audit the quartic trace after the canonical cubic preparation.

Perturb the exact canonical slack factor by the L230 column

    c^3 C_3,   C_3 = 3 (S*)^2 V (B_1* B_1),

and complete the metric by the variable Stein inverse.  The prepared
physical upper-gap endpoint has the quartic trace

    12 ||B_2||_F^2
    + 32 ||B_1||_F^2
    + 56 tr((B_1* B_1)^2).

The companion note proves this by exact word reduction.  This checker
regenerates that reduction and audits the full upper/lower Schur
series on unstructured, rank-changing, and completely delayed data.
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
from repeated_crabb_boundary_slack_deflation import (
    physical_metric_root,
)
from repeated_crabb_canonical_cubic_preimage import (
    coboundary_witness,
    polynomial_column,
)
from repeated_crabb_canonical_cubic_selection import (
    cyclic_trace_classes,
)
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
    series_add,
    series_adjoint,
    series_multiply,
)
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_lower_metric_flag import (
    full_lower_endpoint_series,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray
DEGREE = 4


@dataclass(frozen=True)
class CanonicalQuarticTraceRecord:
    """One prepared quartic trace and flag audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    colligation_error: str
    spectral_radius: str
    first_transfer_norm: str
    second_transfer_norm: str
    prepared_cubic_upper_norm: str
    prepared_cubic_lower_norm: str
    quartic_upper_trace: str
    predicted_quartic_upper_trace: str
    quartic_trace_error: str
    quartic_left_flag_dimension: int
    quartic_left_flag_minimum_eigenvalue: str
    quartic_right_flag_dimension: int
    quartic_right_flag_minimum_eigenvalue: str
    delayed_upper_face_error: str
    delayed_lower_face_error: str
    exact_trace_residual_class_count: int
    exact_rank_null_residual_matches: bool
    all_checks_passed: bool


def matrix_adjoint_series(series: list[Matrix]) -> list[Matrix]:
    """Return the coefficientwise adjoint of a matrix series."""

    return [coefficient.conj().T for coefficient in series]


def prepared_endpoint_series(
    partial: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[list[Matrix], list[Matrix]]:
    """Return upper and lower endpoints through the prepared quartic."""

    dimension = len(partial)
    identity = np.eye(dimension, dtype=complex)
    root, inverse_root = physical_metric_root(partial, right, left)
    physical_operator = inverse_root @ partial @ root
    physical_series = ellipse_operator_coefficients(
        physical_operator,
        DEGREE,
    )
    operator = [
        root @ coefficient @ inverse_root
        for coefficient in physical_series
    ]

    upper_gap, _ = canonical_repair_metric_series(
        partial,
        right,
        left,
        DEGREE,
    )
    equality_metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    upper_target = 4 * np.linalg.inv(equality_metric)
    canonical_metric = [upper_target - upper_gap[0]] + [
        -coefficient for coefficient in upper_gap[1:]
    ]

    pulled = matrix_series_multiply(
        matrix_series_multiply(
            matrix_adjoint_series(operator),
            canonical_metric,
            DEGREE,
        ),
        operator,
        DEGREE,
    )
    canonical_slack = [
        metric - image
        for metric, image in zip(
            canonical_metric,
            pulled,
            strict=True,
        )
    ]

    right_projection = right @ right.conj().T
    first_frame = (
        (identity - right_projection)
        @ canonical_slack[1]
        @ right
        + 0.5
        * right
        @ (
            right.conj().T
            @ canonical_slack[1]
            @ right
        )
    )
    cubic_column = polynomial_column(partial, right, left)
    cubic_slack_change = (
        right @ cubic_column.conj().T
        + cubic_column @ right.conj().T
    )
    quartic_slack_change = (
        first_frame @ cubic_column.conj().T
        + cubic_column @ first_frame.conj().T
    )

    cubic_metric_change = stein_inverse(
        partial,
        cubic_slack_change,
    )
    quartic_forcing = (
        quartic_slack_change
        + operator[1].conj().T
        @ cubic_metric_change
        @ partial
        + partial.conj().T
        @ cubic_metric_change
        @ operator[1]
    )
    quartic_metric_change = stein_inverse(
        partial,
        quartic_forcing,
    )

    prepared_metric = [
        coefficient.copy() for coefficient in canonical_metric
    ]
    prepared_metric[3] += cubic_metric_change
    prepared_metric[4] += quartic_metric_change

    upper_gap = [upper_target - prepared_metric[0]] + [
        -coefficient for coefficient in prepared_metric[1:]
    ]
    lower_target = np.linalg.inv(equality_metric)
    lower_gap = [prepared_metric[0] - lower_target] + [
        coefficient for coefficient in prepared_metric[1:]
    ]
    upper_endpoint, upper_constant_error = full_endpoint_gap_series(
        upper_gap,
        left,
    )
    lower_endpoint, lower_constant_error = (
        full_lower_endpoint_series(lower_gap, right)
    )
    if max(upper_constant_error, lower_constant_error) > 3e-8:
        raise RuntimeError("the prepared endpoint base was inconsistent")
    return (
        [4 * coefficient for coefficient in upper_endpoint],
        lower_endpoint,
    )


def raw_slack_series(degree: int) -> list[Polynomial]:
    """Return the exact boundary-metric Stein slack series."""

    operator = operator_series(degree)
    metric = boundary_metric_series(degree)
    pulled = series_multiply(
        series_multiply(series_adjoint(operator), metric),
        operator,
    )
    return series_add(
        metric,
        [scale(-1, coefficient) for coefficient in pulled],
    )


def exact_quartic_components() -> tuple[Polynomial, Polynomial, Polynomial]:
    """Return quartic Stein forcing, second metric, and Schur square."""

    operator = operator_series(DEGREE)
    boundary = boundary_metric_series(DEGREE)
    raw_slack = raw_slack_series(DEGREE)
    schur_residual = schur_residual_series(DEGREE)
    canonical_slack = [
        add(raw, scale(-1, residual))
        for raw, residual in zip(
            raw_slack,
            schur_residual,
            strict=True,
        )
    ]

    prepared_metric = [{} for _ in range(DEGREE + 1)]
    prepared_metric[0] = IDENTITY
    prepared_metric[2] = add(
        boundary[2],
        scale(-1, schur_residual[2]),
    )
    prepared_metric[3] = coboundary_witness()

    lifted_gram = multiply(
        multiply(
            multiply(
                multiply(E, S),
                F,
            ),
            STAR,
        ),
        E,
    )
    first_frame_factor = add(
        IDENTITY,
        scale(Fraction(-1, 2), E),
    )
    one_cross = multiply(
        multiply(
            multiply(
                first_frame_factor,
                raw_slack[1],
            ),
            lifted_gram,
        ),
        {"ss": 1},
    )
    quartic_slack_change = scale(
        3,
        add(one_cross, adjoint(one_cross)),
    )
    prepared_slack_four = add(
        canonical_slack[4],
        quartic_slack_change,
    )

    variable_operator_terms: Polynomial = {}
    for left_degree in range(DEGREE + 1):
        for metric_degree in range(DEGREE):
            right_degree = (
                DEGREE - left_degree - metric_degree
            )
            if 0 <= right_degree <= DEGREE:
                variable_operator_terms = add(
                    variable_operator_terms,
                    multiply(
                        multiply(
                            adjoint(operator[left_degree]),
                            prepared_metric[metric_degree],
                        ),
                        operator[right_degree],
                    ),
                )
    quartic_stein_forcing = add(
        prepared_slack_four,
        variable_operator_terms,
    )

    upper_complement_inverse = add(
        IDENTITY,
        scale(-1, F),
        scale(Fraction(-2, 3), E),
    )
    upper_schur_square = multiply(
        multiply(
            multiply(
                prepared_metric[2],
                upper_complement_inverse,
            ),
            prepared_metric[2],
        ),
        F,
    )
    return (
        quartic_stein_forcing,
        prepared_metric[2],
        upper_schur_square,
    )


def exact_quartic_trace_polynomial() -> Polynomial:
    """Return the state polynomial whose trace is the upper trace."""

    quartic_stein_forcing, _, upper_schur_square = (
        exact_quartic_components()
    )
    return scale(
        -4,
        add(quartic_stein_forcing, upper_schur_square),
    )


def predicted_trace_polynomial() -> Polynomial:
    """Return the state lift of the predicted positive trace."""

    first_gram = multiply(
        multiply(
            multiply(
                multiply(E, S),
                F,
            ),
            STAR,
        ),
        E,
    )
    second_gram = multiply(
        multiply(
            multiply(
                multiply(E, {"ss": 1}),
                F,
            ),
            {"aa": 1},
        ),
        E,
    )
    return add(
        scale(12, second_gram),
        scale(32, first_gram),
        scale(56, multiply(first_gram, first_gram)),
    )


@lru_cache(maxsize=1)
def exact_trace_residual() -> tuple[int, bool]:
    """Audit the sole rank-null trace residual after word reduction."""

    residual_classes = cyclic_trace_classes(
        add(
            exact_quartic_trace_polynomial(),
            scale(-1, predicted_trace_polynomial()),
        )
    )
    expected = {
        "": Fraction(-8),
        "as": Fraction(16),
        "aass": Fraction(-8),
    }
    return len(residual_classes), residual_classes == expected


def stable_kernel(matrix: Matrix) -> Matrix:
    """Return a stable kernel frame for a rank-changing transfer."""

    _, singular_values, adjoint = np.linalg.svd(
        matrix,
        full_matrices=True,
    )
    if not len(singular_values):
        return adjoint.conj().T
    tolerance = max(1e-10, 1e-8 * singular_values[0])
    rank = int(np.count_nonzero(singular_values > tolerance))
    return adjoint.conj().T[:, rank:]


def minimum_flag_eigenvalue(
    coefficient: Matrix,
    kernel: Matrix,
) -> float:
    """Return the least eigenvalue on one flag, or zero if empty."""

    if kernel.shape[1] == 0:
        return 0.0
    compression = kernel.conj().T @ coefficient @ kernel
    compression = (compression + compression.conj().T) / 2
    return float(np.linalg.eigvalsh(compression)[0])


def audit_case(
    construction_kind: str,
    partial: Matrix,
    right: Matrix,
    left: Matrix,
    parameter_scale: float,
) -> CanonicalQuarticTraceRecord:
    """Audit one prepared quartic trace and its two flags."""

    upper, lower = prepared_endpoint_series(partial, right, left)
    identity = np.eye(len(partial), dtype=complex)
    colligation_error = max(
        float(
            np.linalg.norm(
                identity
                - partial.conj().T @ partial
                - right @ right.conj().T
            )
        ),
        float(
            np.linalg.norm(
                identity
                - partial @ partial.conj().T
                - left @ left.conj().T
            )
        ),
        float(np.linalg.norm(right.conj().T @ left)),
    )
    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(partial)))
    )
    first = transfer_coefficient(partial, right, left, 1)
    second = transfer_coefficient(partial, right, left, 2)
    first_right_gram = first.conj().T @ first
    predicted_trace = (
        12 * np.linalg.norm(second) ** 2
        + 32 * np.linalg.norm(first) ** 2
        + 56 * float(np.trace(first_right_gram @ first_right_gram).real)
    )
    trace = float(np.trace(upper[4]).real)
    trace_error = abs(trace - predicted_trace)

    left_kernel = stable_kernel(first.conj().T)
    right_kernel = stable_kernel(first)
    left_minimum = minimum_flag_eigenvalue(
        upper[4],
        left_kernel,
    )
    right_minimum = minimum_flag_eigenvalue(
        lower[4],
        right_kernel,
    )

    delayed = np.linalg.norm(first) < 1e-9
    delayed_upper_error = 0.0
    delayed_lower_error = 0.0
    if delayed:
        delayed_upper_error = float(
            np.linalg.norm(
                upper[4] - 12 * second @ second.conj().T
            )
        )
        delayed_lower_error = float(
            np.linalg.norm(
                lower[4] - second.conj().T @ second
            )
        )

    exact_class_count, exact_match = exact_trace_residual()
    tolerance = 3e-8
    verified = bool(
        colligation_error < tolerance
        and spectral_radius < 1
        and np.linalg.norm(upper[3]) < tolerance
        and np.linalg.norm(lower[3]) < tolerance
        and trace_error < tolerance
        and right_minimum > -tolerance
        and delayed_upper_error < tolerance
        and delayed_lower_error < tolerance
        and exact_class_count == 3
        and exact_match
    )
    if not verified:
        raise RuntimeError(
            "the prepared quartic trace audit failed: "
            f"kind={construction_kind}, "
            f"colligation={colligation_error:.3e}, "
            f"radius={spectral_radius:.6f}, "
            f"trace={trace_error:.3e}, "
            f"right-min={right_minimum:.3e}, "
            f"delayed-upper={delayed_upper_error:.3e}, "
            f"delayed-lower={delayed_lower_error:.3e}"
        )

    return CanonicalQuarticTraceRecord(
        construction_kind=construction_kind,
        state_dimension=len(partial),
        defect_dimension=right.shape[1],
        parameter_scale=format_float(parameter_scale),
        colligation_error=format_float(colligation_error),
        spectral_radius=format_float(spectral_radius),
        first_transfer_norm=format_float(float(np.linalg.norm(first))),
        second_transfer_norm=format_float(float(np.linalg.norm(second))),
        prepared_cubic_upper_norm=format_float(
            float(np.linalg.norm(upper[3]))
        ),
        prepared_cubic_lower_norm=format_float(
            float(np.linalg.norm(lower[3]))
        ),
        quartic_upper_trace=format_float(trace),
        predicted_quartic_upper_trace=format_float(predicted_trace),
        quartic_trace_error=format_float(trace_error),
        quartic_left_flag_dimension=left_kernel.shape[1],
        quartic_left_flag_minimum_eigenvalue=format_float(left_minimum),
        quartic_right_flag_dimension=right_kernel.shape[1],
        quartic_right_flag_minimum_eigenvalue=format_float(right_minimum),
        delayed_upper_face_error=format_float(delayed_upper_error),
        delayed_lower_face_error=format_float(delayed_lower_error),
        exact_trace_residual_class_count=exact_class_count,
        exact_rank_null_residual_matches=exact_match,
        all_checks_passed=verified,
    )


def standard_records() -> list[CanonicalQuarticTraceRecord]:
    """Return deterministic unstructured, flag, and delay records."""

    records: list[CanonicalQuarticTraceRecord] = []
    for defect_dimension in range(1, 5):
        state_dimension = 3 * defect_dimension + 3
        partial, right, left = random_partial_isometry(
            state_dimension,
            defect_dimension,
            np.random.default_rng(112_000 + defect_dimension),
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
        for parameter_scale in (1.0, 0.5, 0.2):
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

    generator = np.random.default_rng(112_800)
    for defect_dimension in (2, 3, 4):
        state_dimension = 4 * defect_dimension
        partial, right, left = delayed_random_partial_isometry(
            state_dimension,
            defect_dimension,
            generator,
        )
        records.append(
            audit_case(
                "complete_delay",
                partial,
                right,
                left,
                1,
            )
        )
    return records


def write_records(
    records: list[CanonicalQuarticTraceRecord],
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
            "repeated_crabb_canonical_quartic_trace_s70224.jsonl"
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
