#!/usr/bin/env python3
"""Audit exact two-graph Schur transport and finite-jet domination."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from crabb_block_hardy_equality import format_float


Matrix = np.ndarray
MatrixSeries = list[Matrix]


@dataclass(frozen=True)
class ExactTwoGraphTransportRecord:
    """One exact rational audit of the mixed-graph identity."""

    record_type: str
    state_dimension: int
    endpoint_dimension: int
    old_interior_determinant: str
    new_interior_determinant: str
    left_schur_residual_nonzero_entry_count: int
    right_schur_residual_nonzero_entry_count: int
    mixed_graph_residual_nonzero_entry_count: int
    adjoint_mixed_graph_residual_nonzero_entry_count: int
    all_checks_passed: bool


@dataclass(frozen=True)
class SeriesTwoGraphTransportRecord:
    """One coefficientwise audit of the triangular series identity."""

    record_type: str
    state_dimension: int
    endpoint_dimension: int
    maximum_degree: int
    first_perturbation_degree: int
    endpoint_difference_before_perturbation_norm: str
    maximum_mixed_graph_coefficient_error: str
    maximum_endpoint_hermitian_error: str
    all_checks_passed: bool


@dataclass(frozen=True)
class FiniteJetDominationRecord:
    """One deterministic audit of the terminal-Gram tail estimate."""

    record_type: str
    copy_dimension: int
    terminal_grade: int
    terminal_minimum_singular_value: str
    certified_radius: str
    tested_parameter: str
    finite_jet_minimum_eigenvalue: str
    full_endpoint_minimum_eigenvalue: str
    analytic_tail_norm: str
    all_checks_passed: bool


AuditRecord = (
    ExactTwoGraphTransportRecord
    | SeriesTwoGraphTransportRecord
    | FiniteJetDominationRecord
)


def sympy_graph(
    matrix: sp.Matrix,
    endpoint_dimension: int,
) -> sp.Matrix:
    """Return the exact endpoint graph column of one block matrix."""

    cross_adjoint = matrix[endpoint_dimension:, :endpoint_dimension]
    interior = matrix[endpoint_dimension:, endpoint_dimension:]
    lower = -interior.inv() * cross_adjoint
    return sp.eye(endpoint_dimension).col_join(lower)


def sympy_schur(
    matrix: sp.Matrix,
    endpoint_dimension: int,
) -> sp.Matrix:
    """Return the exact Schur complement onto the endpoint block."""

    endpoint = matrix[:endpoint_dimension, :endpoint_dimension]
    cross = matrix[:endpoint_dimension, endpoint_dimension:]
    cross_adjoint = matrix[endpoint_dimension:, :endpoint_dimension]
    interior = matrix[endpoint_dimension:, endpoint_dimension:]
    return sp.simplify(
        endpoint - cross * interior.inv() * cross_adjoint
    )


def nonzero_entry_count(matrix: sp.Matrix) -> int:
    """Count nonzero entries after exact simplification."""

    return sum(int(sp.simplify(entry) != 0) for entry in matrix)


def exact_record() -> ExactTwoGraphTransportRecord:
    """Verify the two orientations of the identity over rationals."""

    old_endpoint = sp.Matrix([[3, 1], [1, 2]])
    old_cross = sp.Matrix([[1, -1, 2], [0, 1, 1]])
    old_interior = sp.Matrix(
        [
            [4, 1, 0],
            [1, 3, 1],
            [0, 1, 5],
        ]
    )
    old = old_endpoint.row_join(old_cross).col_join(
        old_cross.T.row_join(old_interior)
    )

    perturbation = sp.Matrix(
        [
            [1, 0, 1, -1, 0],
            [0, -1, 0, 1, 1],
            [1, 0, 1, 0, 1],
            [-1, 1, 0, 2, -1],
            [0, 1, 1, -1, 1],
        ]
    )
    new = old + sp.Rational(1, 3) * perturbation
    endpoint_dimension = 2

    old_graph = sympy_graph(old, endpoint_dimension)
    new_graph = sympy_graph(new, endpoint_dimension)
    old_schur = sympy_schur(old, endpoint_dimension)
    new_schur = sympy_schur(new, endpoint_dimension)
    difference = sp.simplify(new_schur - old_schur)
    mixed = sp.simplify(new_graph.T * (new - old) * old_graph)
    adjoint_mixed = sp.simplify(
        old_graph.T * (new - old) * new_graph
    )
    old_left_residual = sp.simplify(
        old * old_graph
        - old_schur.col_join(
            sp.zeros(old.rows - endpoint_dimension, endpoint_dimension)
        )
    )
    new_left_residual = sp.simplify(
        new * new_graph
        - new_schur.col_join(
            sp.zeros(new.rows - endpoint_dimension, endpoint_dimension)
        )
    )
    residuals = (
        nonzero_entry_count(old_left_residual),
        nonzero_entry_count(new_left_residual),
        nonzero_entry_count(difference - mixed),
        nonzero_entry_count(difference - adjoint_mixed),
    )
    verified = bool(
        old_interior.det() != 0
        and new[endpoint_dimension:, endpoint_dimension:].det() != 0
        and max(residuals) == 0
    )
    if not verified:
        raise RuntimeError(
            "the exact two-graph certificate failed: "
            f"residuals={residuals}"
        )
    return ExactTwoGraphTransportRecord(
        record_type="exact_rational_two_graph_identity",
        state_dimension=old.rows,
        endpoint_dimension=endpoint_dimension,
        old_interior_determinant=str(old_interior.det()),
        new_interior_determinant=str(
            new[endpoint_dimension:, endpoint_dimension:].det()
        ),
        left_schur_residual_nonzero_entry_count=residuals[0],
        right_schur_residual_nonzero_entry_count=residuals[1],
        mixed_graph_residual_nonzero_entry_count=residuals[2],
        adjoint_mixed_graph_residual_nonzero_entry_count=residuals[3],
        all_checks_passed=verified,
    )


def zero_series(
    rows: int,
    columns: int,
    maximum_degree: int,
) -> MatrixSeries:
    """Return a zero rectangular matrix series."""

    return [
        np.zeros((rows, columns), dtype=complex)
        for _ in range(maximum_degree + 1)
    ]


def series_add(
    left: MatrixSeries,
    right: MatrixSeries,
) -> MatrixSeries:
    """Add two equally shaped matrix series."""

    return [
        left_coefficient + right_coefficient
        for left_coefficient, right_coefficient in zip(
            left,
            right,
            strict=True,
        )
    ]


def series_multiply(
    left: MatrixSeries,
    right: MatrixSeries,
) -> MatrixSeries:
    """Multiply compatible rectangular matrix series."""

    maximum_degree = min(len(left), len(right)) - 1
    result = zero_series(
        left[0].shape[0],
        right[0].shape[1],
        maximum_degree,
    )
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            degree = left_degree + right_degree
            if degree <= maximum_degree:
                result[degree] += (
                    left_coefficient @ right_coefficient
                )
    return result


def series_adjoint(series: MatrixSeries) -> MatrixSeries:
    """Apply the matrix adjoint coefficientwise."""

    return [coefficient.conj().T for coefficient in series]


def inverse_series(series: MatrixSeries) -> MatrixSeries:
    """Invert a square matrix series with invertible constant term."""

    result = [np.linalg.inv(series[0])]
    for degree in range(1, len(series)):
        convolution = np.zeros_like(series[0])
        for positive_degree in range(1, degree + 1):
            convolution += (
                series[positive_degree]
                @ result[degree - positive_degree]
            )
        result.append(-result[0] @ convolution)
    return result


def block_series(
    series: MatrixSeries,
    row_slice: slice,
    column_slice: slice,
) -> MatrixSeries:
    """Extract one fixed rectangular block coefficientwise."""

    return [
        coefficient[row_slice, column_slice]
        for coefficient in series
    ]


def graph_series(
    series: MatrixSeries,
    endpoint_dimension: int,
) -> MatrixSeries:
    """Return the endpoint graph column as a matrix series."""

    state_dimension = series[0].shape[0]
    maximum_degree = len(series) - 1
    endpoint_slice = slice(0, endpoint_dimension)
    interior_slice = slice(endpoint_dimension, state_dimension)
    interior = block_series(
        series,
        interior_slice,
        interior_slice,
    )
    cross_adjoint = block_series(
        series,
        interior_slice,
        endpoint_slice,
    )
    lower = [
        -coefficient
        for coefficient in series_multiply(
            inverse_series(interior),
            cross_adjoint,
        )
    ]
    result = zero_series(
        state_dimension,
        endpoint_dimension,
        maximum_degree,
    )
    result[0][:endpoint_dimension, :] = np.eye(
        endpoint_dimension,
        dtype=complex,
    )
    for degree in range(maximum_degree + 1):
        result[degree][endpoint_dimension:, :] = lower[degree]
    return result


def schur_series(
    series: MatrixSeries,
    endpoint_dimension: int,
) -> MatrixSeries:
    """Return the endpoint Schur-complement series."""

    state_dimension = series[0].shape[0]
    endpoint_slice = slice(0, endpoint_dimension)
    interior_slice = slice(endpoint_dimension, state_dimension)
    endpoint = block_series(
        series,
        endpoint_slice,
        endpoint_slice,
    )
    cross = block_series(
        series,
        endpoint_slice,
        interior_slice,
    )
    cross_adjoint = block_series(
        series,
        interior_slice,
        endpoint_slice,
    )
    interior = block_series(
        series,
        interior_slice,
        interior_slice,
    )
    cross_square = series_multiply(
        series_multiply(cross, inverse_series(interior)),
        cross_adjoint,
    )
    return [
        endpoint_coefficient - cross_coefficient
        for endpoint_coefficient, cross_coefficient in zip(
            endpoint,
            cross_square,
            strict=True,
        )
    ]


def random_hermitian(
    generator: np.random.Generator,
    dimension: int,
) -> Matrix:
    """Return one deterministic complex Hermitian matrix."""

    raw = (
        generator.standard_normal((dimension, dimension))
        + 1j * generator.standard_normal((dimension, dimension))
    )
    return (raw + raw.conj().T) / 2


def series_record() -> SeriesTwoGraphTransportRecord:
    """Audit every coefficient of the formal two-graph identity."""

    generator = np.random.default_rng(702_281)
    state_dimension = 7
    endpoint_dimension = 3
    maximum_degree = 8
    first_perturbation_degree = 4
    old = zero_series(
        state_dimension,
        state_dimension,
        maximum_degree,
    )
    old[0] = np.diag(
        np.linspace(2.0, 4.0, state_dimension)
    ).astype(complex)
    for degree in range(1, maximum_degree + 1):
        old[degree] = (
            0.025 / degree
        ) * random_hermitian(generator, state_dimension)

    perturbation = zero_series(
        state_dimension,
        state_dimension,
        maximum_degree,
    )
    for degree in range(
        first_perturbation_degree,
        maximum_degree + 1,
    ):
        perturbation[degree] = (
            0.018 / degree
        ) * random_hermitian(generator, state_dimension)
    new = series_add(old, perturbation)

    old_schur = schur_series(old, endpoint_dimension)
    new_schur = schur_series(new, endpoint_dimension)
    endpoint_difference = [
        new_coefficient - old_coefficient
        for new_coefficient, old_coefficient in zip(
            new_schur,
            old_schur,
            strict=True,
        )
    ]
    old_graph = graph_series(old, endpoint_dimension)
    new_graph = graph_series(new, endpoint_dimension)
    mixed = series_multiply(
        series_multiply(
            series_adjoint(new_graph),
            perturbation,
        ),
        old_graph,
    )
    before_norm = max(
        np.linalg.norm(endpoint_difference[degree])
        for degree in range(first_perturbation_degree)
    )
    mixed_error = max(
        np.linalg.norm(
            endpoint_difference[degree] - mixed[degree]
        )
        for degree in range(maximum_degree + 1)
    )
    hermitian_error = max(
        np.linalg.norm(
            coefficient - coefficient.conj().T
        )
        for coefficient in endpoint_difference
    )
    tolerance = 2e-12
    verified = bool(
        before_norm < tolerance
        and mixed_error < tolerance
        and hermitian_error < tolerance
    )
    if not verified:
        raise RuntimeError(
            "the series two-graph audit failed: "
            f"before={before_norm:.3e}, "
            f"mixed={mixed_error:.3e}, "
            f"hermitian={hermitian_error:.3e}"
        )
    return SeriesTwoGraphTransportRecord(
        record_type="formal_series_two_graph_identity",
        state_dimension=state_dimension,
        endpoint_dimension=endpoint_dimension,
        maximum_degree=maximum_degree,
        first_perturbation_degree=first_perturbation_degree,
        endpoint_difference_before_perturbation_norm=format_float(
            before_norm
        ),
        maximum_mixed_graph_coefficient_error=format_float(
            mixed_error
        ),
        maximum_endpoint_hermitian_error=format_float(
            hermitian_error
        ),
        all_checks_passed=verified,
    )


def finite_jet_record() -> FiniteJetDominationRecord:
    """Audit the elementary terminal-Gram remainder domination."""

    generator = np.random.default_rng(702_282)
    copy_dimension = 4
    terminal_grade = 5
    terminal = random_hermitian(generator, copy_dimension)
    terminal += 3.5 * np.eye(copy_dimension)
    minimum_singular = float(
        np.linalg.svd(terminal, compute_uv=False)[-1]
    )
    margin = 0.7
    tail_bound = 2.25
    certified_radius = (
        margin * minimum_singular**2 / (2 * tail_bound)
    )
    parameter = min(0.08, certified_radius / 2)

    finite_jet = (
        margin
        * parameter ** (2 * terminal_grade)
        * terminal
        @ terminal.conj().T
    )
    for grade in range(1, terminal_grade):
        factor = (
            0.15
            * parameter**grade
            * random_hermitian(generator, copy_dimension)
        )
        finite_jet += factor @ factor.conj().T

    tail_direction = random_hermitian(generator, copy_dimension)
    tail_direction *= (
        tail_bound / np.linalg.norm(tail_direction, ord=2)
    )
    tail = parameter ** (2 * terminal_grade + 1) * tail_direction
    full = finite_jet + tail
    finite_minimum = float(np.linalg.eigvalsh(finite_jet)[0])
    full_minimum = float(np.linalg.eigvalsh(full)[0])
    tail_norm = float(np.linalg.norm(tail, ord=2))
    predicted_floor = (
        parameter ** (2 * terminal_grade)
        * (
            margin * minimum_singular**2
            - parameter * tail_bound
        )
    )
    tolerance = 2e-12
    verified = bool(
        parameter < certified_radius
        and finite_minimum
        >= (
            margin
            * minimum_singular**2
            * parameter ** (2 * terminal_grade)
            - tolerance
        )
        and full_minimum >= predicted_floor - tolerance
        and full_minimum > 0
    )
    if not verified:
        raise RuntimeError(
            "the finite-jet domination audit failed: "
            f"radius={certified_radius:.3e}, "
            f"parameter={parameter:.3e}, "
            f"minimum={full_minimum:.3e}, "
            f"floor={predicted_floor:.3e}"
        )
    return FiniteJetDominationRecord(
        record_type="terminal_gram_finite_jet_domination",
        copy_dimension=copy_dimension,
        terminal_grade=terminal_grade,
        terminal_minimum_singular_value=format_float(
            minimum_singular
        ),
        certified_radius=format_float(certified_radius),
        tested_parameter=format_float(parameter),
        finite_jet_minimum_eigenvalue=format_float(finite_minimum),
        full_endpoint_minimum_eigenvalue=format_float(full_minimum),
        analytic_tail_norm=format_float(tail_norm),
        all_checks_passed=verified,
    )


def standard_records() -> list[AuditRecord]:
    """Return exact, series, and finite-tail audits."""

    return [
        exact_record(),
        series_record(),
        finite_jet_record(),
    ]


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
            "repeated_crabb_two_graph_transport_s70226.jsonl"
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
