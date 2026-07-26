#!/usr/bin/env python3
"""Audit the nonlinear defect of L298 under the moving elliptic pair."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_elliptic_cokernel import stein_inverse
from repeated_crabb_oriented_retightening_transport import (
    oriented_retightening_direction,
)
from repeated_crabb_transfer_channel_covariance import (
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class MovingDefectRecord:
    """One exact or numerical cubic moving-defect audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    parameter_scale: str
    first_transfer_rank: int
    lower_corner_norm: str
    endpoint_trace_absolute: str
    endpoint_norm: str
    flag_compression_norm: str
    flag_compression_value: str
    endpoint_minimum_eigenvalue: str
    endpoint_maximum_eigenvalue: str
    all_checks_passed: bool


def ellipse_tangent(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> Matrix:
    """Return the first balanced pulled-ellipse coefficient."""

    identity = np.eye(len(operator), dtype=complex)
    right_projection = right @ right.conj().T
    left_projection = left @ left.conj().T
    reflection = (
        (identity + left_projection)
        @ operator.conj().T
        @ (identity + right_projection)
    )
    return reflection - np.linalg.matrix_power(operator, 3)


def raw_frame_tangent(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    operator_tangent: Matrix,
) -> Matrix:
    """Return the first coefficient of L227's raw slack factor."""

    slack_tangent = -(
        operator_tangent.conj().T @ operator
        + operator.conj().T @ operator_tangent
    )
    pivot_tangent = (
        right.conj().T @ slack_tangent @ right
    )
    return (
        slack_tangent @ right
        - right @ pivot_tangent / 2
    )


def cubic_moving_defect(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, Matrix, Matrix]:
    """Return the cubic state defect, its endpoint, and ``B_1``."""

    operator_tangent = ellipse_tangent(operator, right, left)
    frame_tangent = raw_frame_tangent(
        operator,
        right,
        left,
        operator_tangent,
    )
    direction = oriented_retightening_direction(
        operator,
        right,
        left,
        1,
    )
    state_defect = -(
        operator_tangent.conj().T
        @ direction.metric
        @ operator
        + operator.conj().T
        @ direction.metric
        @ operator_tangent
        + frame_tangent @ direction.frame.conj().T
        + direction.frame @ frame_tangent.conj().T
    )
    state_defect = (state_defect + state_defect.conj().T) / 2
    endpoint = (
        left.conj().T
        @ stein_inverse(operator, state_defect)
        @ left
    )
    endpoint = (endpoint + endpoint.conj().T) / 2
    first_transfer = transfer_coefficient(
        operator,
        right,
        left,
        1,
    )
    return state_defect, endpoint, first_transfer


def left_kernel(matrix: Matrix) -> Matrix:
    """Return an orthonormal basis for ``ker(matrix*)``."""

    gram = matrix @ matrix.conj().T
    eigenvalues, eigenvectors = np.linalg.eigh(
        (gram + gram.conj().T) / 2
    )
    tolerance = max(1e-11, 1e-8 * float(eigenvalues[-1]))
    return eigenvectors[:, eigenvalues < tolerance]


def numerical_record(
    multiplicity: int,
    parameter_scale: float,
    seed: int,
) -> MovingDefectRecord:
    """Audit one rank-deficient noncommuting Schur colligation."""

    operator, right, left, colligation_error, _ = rank_chain_case(
        multiplicity,
        seed,
        parameter_scale,
    )
    state_defect, endpoint, first_transfer = cubic_moving_defect(
        operator,
        right,
        left,
    )
    kernel = left_kernel(first_transfer)
    compression = kernel.conj().T @ endpoint @ kernel
    compression = (compression + compression.conj().T) / 2
    endpoint_eigenvalues = np.linalg.eigvalsh(endpoint)
    lower_norm = float(
        np.linalg.norm(
            right.conj().T @ state_defect @ right
        )
    )
    trace_absolute = abs(float(np.trace(endpoint).real))
    endpoint_norm = float(np.linalg.norm(endpoint))
    compression_norm = float(np.linalg.norm(compression))
    tolerance = 3e-8
    verified = bool(
        colligation_error < tolerance
        and np.linalg.matrix_rank(first_transfer, 1e-8) == 1
        and kernel.shape[1] == multiplicity - 1
        and lower_norm < tolerance
        and trace_absolute < tolerance
        and endpoint_norm > 1e-12
        and compression_norm > 1e-12
        and endpoint_eigenvalues[0] < -1e-12
        and endpoint_eigenvalues[-1] > 1e-12
    )
    if not verified:
        raise RuntimeError(
            "the moving-retightening cubic audit failed: "
            f"{multiplicity=}, {parameter_scale=}"
        )
    return MovingDefectRecord(
        construction_kind="noncommuting_rank_chain",
        state_dimension=len(operator),
        defect_dimension=multiplicity,
        parameter_scale=format_float(parameter_scale),
        first_transfer_rank=1,
        lower_corner_norm=format_float(lower_norm),
        endpoint_trace_absolute=format_float(trace_absolute),
        endpoint_norm=format_float(endpoint_norm),
        flag_compression_norm=format_float(compression_norm),
        flag_compression_value="not_applicable",
        endpoint_minimum_eigenvalue=format_float(
            float(endpoint_eigenvalues[0])
        ),
        endpoint_maximum_eigenvalue=format_float(
            float(endpoint_eigenvalues[-1])
        ),
        all_checks_passed=verified,
    )


def exact_stein_inverse(
    operator: sp.Matrix,
    forcing: sp.Matrix,
) -> sp.Matrix:
    """Solve one exact rational Stein equation."""

    dimension = operator.rows
    variables = sp.symbols(f"x:{dimension * dimension}")
    unknown = sp.Matrix(dimension, dimension, variables)
    solution_set = sp.linsolve(
        list(unknown - operator.T * unknown * operator - forcing),
        variables,
    )
    solution = next(iter(solution_set))
    return sp.simplify(
        unknown.subs(dict(zip(variables, solution, strict=True)))
    )


def exact_rational_record() -> MovingDefectRecord:
    """Return a rigorous rank-one-flag counterexample."""

    dimension = 6
    multiplicity = 2
    cosine_one = sp.Rational(3, 5)
    sine_one = sp.Rational(4, 5)
    cosine_two = sp.Rational(5, 13)
    sine_two = sp.Rational(12, 13)
    rotation_one = sp.eye(4)
    rotation_one[0, 0] = rotation_one[1, 1] = cosine_one
    rotation_one[0, 1] = -sine_one
    rotation_one[1, 0] = sine_one
    rotation_two = sp.eye(4)
    rotation_two[1, 1] = rotation_two[2, 2] = cosine_two
    rotation_two[1, 2] = -sine_two
    rotation_two[2, 1] = sine_two
    bridge = rotation_two * rotation_one

    identity = sp.eye(dimension)
    right = identity[:, :multiplicity]
    left = identity[:, 4:]
    domain = identity[:, 2:]
    range_frame = identity[:, :4]
    operator = range_frame * bridge * domain.T
    right_projection = right * right.T
    left_projection = left * left.T
    reflection = (
        (identity + left_projection)
        * operator.T
        * (identity + right_projection)
    )
    operator_tangent = reflection - operator**3
    slack_tangent = -(
        operator_tangent.T * operator
        + operator.T * operator_tangent
    )
    pivot_tangent = right.T * slack_tangent * right
    frame_tangent = (
        slack_tangent * right
        - right * pivot_tangent / 2
    )

    first_transfer = left.T * operator.T * right
    right_gram = first_transfer.T * first_transfer
    projection = identity - right_projection
    perpendicular = (
        -sp.Rational(1, 2)
        * projection
        * operator
        * left
        * first_transfer
    )
    orbit = exact_stein_inverse(
        operator,
        right * right_gram * right.T,
    )
    perpendicular_metric = exact_stein_inverse(
        operator,
        right * perpendicular.T + perpendicular * right.T,
    )
    metric_direction = -orbit + perpendicular_metric
    frame_direction = (
        -right * right_gram / 2 + perpendicular
    )
    state_defect = -(
        operator_tangent.T * metric_direction * operator
        + operator.T * metric_direction * operator_tangent
        + frame_tangent * frame_direction.T
        + frame_direction * frame_tangent.T
    )
    endpoint = sp.simplify(
        left.T
        * exact_stein_inverse(operator, state_defect)
        * left
    )
    expected_endpoint = sp.Matrix(
        [
            [
                -sp.Rational(2_985_984, 54_865_681),
                -sp.Rational(139_470_336, 3_566_269_265),
            ],
            [
                -sp.Rational(139_470_336, 3_566_269_265),
                sp.Rational(2_985_984, 54_865_681),
            ],
        ]
    )
    flag = sp.Matrix([0, 1])
    flag_value = sp.factor((flag.T * endpoint * flag)[0])
    expected_flag_value = sp.Rational(2_985_984, 54_865_681)
    exact_checks = (
        operator.T * operator == identity - right_projection,
        operator * operator.T == identity - left_projection,
        first_transfer.rank() == 1,
        right.T * state_defect * right == sp.zeros(multiplicity),
        endpoint == expected_endpoint,
        sp.trace(endpoint) == 0,
        flag_value == expected_flag_value,
    )
    if not all(exact_checks):
        raise RuntimeError("the exact moving-defect guard failed")

    endpoint_float = np.asarray(endpoint, dtype=float)
    endpoint_eigenvalues = np.linalg.eigvalsh(endpoint_float)
    return MovingDefectRecord(
        construction_kind="exact_rational_rank_one_flag",
        state_dimension=dimension,
        defect_dimension=multiplicity,
        parameter_scale="not_applicable",
        first_transfer_rank=1,
        lower_corner_norm="0",
        endpoint_trace_absolute="0",
        endpoint_norm=format_float(
            float(np.linalg.norm(endpoint_float))
        ),
        flag_compression_norm=str(expected_flag_value),
        flag_compression_value=str(expected_flag_value),
        endpoint_minimum_eigenvalue=format_float(
            float(endpoint_eigenvalues[0])
        ),
        endpoint_maximum_eigenvalue=format_float(
            float(endpoint_eigenvalues[-1])
        ),
        all_checks_passed=True,
    )


def standard_records() -> list[MovingDefectRecord]:
    """Return the exact guard and nearby noncommuting audits."""

    records = [exact_rational_record()]
    for multiplicity in (2, 3, 4):
        for parameter_scale in (0.1, 0.3, 1.0):
            records.append(
                numerical_record(
                    multiplicity,
                    parameter_scale,
                    550_000 + multiplicity,
                )
            )
    return records


def write_records(
    records: list[MovingDefectRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(
        json.dumps(asdict(record), sort_keys=True) + "\n"
        for record in records
    )
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run every deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = None
    if args.output is not None:
        digest = write_records(records, args.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "dataset_sha256": digest,
                "record_count": len(records),
                "records": [asdict(record) for record in records],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
