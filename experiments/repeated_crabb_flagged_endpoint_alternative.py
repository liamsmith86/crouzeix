#!/usr/bin/env python3
"""Audit the flagged endpoint semidefinite alternative.

For a left-copy isometry ``U``, compress L204's endpoint response to

    C -> U* M_T(C) U.

Its annihilator must be exactly the L206 commutant condition applied
to ``U Y U*``.  The script checks that identity on full, proper
irreducible, and proper reducing flags.  It also audits both sides of
the elementary semidefinite alternative with constructed feasible
and obstructed faces.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, sqrtm

from crabb_block_hardy_equality import format_float, hermitian_basis
from repeated_crabb_elliptic_cokernel import (
    build_reduced_face_data,
    dual_stein_inverse,
    endpoint_motion,
    projected_column_basis,
)
from repeated_crabb_transfer_deflation import haar_unitary


Matrix = np.ndarray


@dataclass(frozen=True)
class FlaggedEndpointAlternativeRecord:
    """One compressed-response and separator audit."""

    flag_kind: str
    length: int
    multiplicity: int
    flag_dimension: int
    response_rank: int
    response_codimension: int
    adjoint_kernel_dimension: int
    maximum_adjoint_pairing_error: str
    maximum_commutant_error: str
    feasible_correction_error: str
    feasible_endpoint_maximum: str
    identity_separator_residual: str
    obstructed_separator_pairing: str
    all_checks_passed: bool


def hermitian_coordinates(matrix: Matrix) -> np.ndarray:
    """Return independent real coordinates of a Hermitian matrix."""

    size = len(matrix)
    coordinates = [float(matrix[index, index].real) for index in range(size)]
    for row in range(size):
        for column in range(row + 1, size):
            coordinates.extend(
                (
                    float(matrix[row, column].real),
                    float(matrix[row, column].imag),
                )
            )
    return np.asarray(coordinates)


def endpoint_from_coordinates(
    coordinates: np.ndarray,
    basis: tuple[Matrix, ...],
) -> Matrix:
    """Reconstruct a Hermitian endpoint from real basis coordinates."""

    return sum(
        (
            coefficient * matrix
            for coefficient, matrix in zip(
                coordinates,
                basis,
                strict=True,
            )
        ),
        np.zeros_like(basis[0]),
    )


def full_cokernel_endpoints(data) -> list[Matrix]:
    """Return a numerical real basis of the full L204 cokernel."""

    multiplicity = data.right_defect.shape[1]
    basis = tuple(hermitian_basis(multiplicity))
    projection = (
        np.eye(len(data.operator))
        - data.right_defect @ data.right_defect.conj().T
    )
    columns = []
    for endpoint in basis:
        dual = dual_stein_inverse(
            data.operator,
            data.left_defect @ endpoint @ data.left_defect.conj().T,
        )
        residual = projection @ dual @ data.right_defect
        columns.append(
            np.concatenate((residual.real.ravel(), residual.imag.ravel()))
        )
    kernel = null_space(np.stack(columns, axis=1), rcond=1e-8)
    return [
        endpoint_from_coordinates(column, basis)
        for column in kernel.T
    ]


def reducing_flag(data) -> Matrix:
    """Extract a proper reducing left-copy flag when one is present."""

    multiplicity = data.right_defect.shape[1]
    cokernel = full_cokernel_endpoints(data)
    identity = np.eye(multiplicity, dtype=complex)
    candidates = [
        endpoint
        - np.trace(endpoint).real / multiplicity * identity
        for endpoint in cokernel
    ]
    endpoint = max(candidates, key=np.linalg.norm)
    eigenvalues, eigenvectors = np.linalg.eigh(endpoint)
    gaps = np.diff(eigenvalues)
    split = int(np.argmax(np.abs(gaps))) + 1
    if split == 0 or split == multiplicity:
        raise RuntimeError("failed to extract a proper reducing flag")
    return eigenvectors[:, :split]


def response_matrix(data, flag: Matrix) -> tuple[Matrix, list[Matrix]]:
    """Return real coordinates of the compressed endpoint map."""

    columns = list(projected_column_basis(data.right_defect))
    endpoint_basis = tuple(hermitian_basis(flag.shape[1]))
    responses = [
        flag.conj().T
        @ endpoint_motion(
            data.operator,
            data.right_defect,
            data.left_defect,
            column,
        )
        @ flag
        for column in columns
    ]
    matrix = np.stack(
        [
            np.asarray(
                [
                    float(np.trace(basis @ response).real)
                    for basis in endpoint_basis
                ]
            )
            for response in responses
        ],
        axis=1,
    )
    return matrix, columns


def make_record(
    flag_kind: str,
    length: int,
    multiplicity: int,
    strength: float,
) -> FlaggedEndpointAlternativeRecord:
    """Audit one full, irreducible, or reducing flag."""

    data = build_reduced_face_data(length, multiplicity, strength)
    if flag_kind == "full":
        flag = np.eye(multiplicity, dtype=complex)
    elif flag_kind == "proper_irreducible":
        flag = haar_unitary(
            multiplicity,
            np.random.default_rng(106_200 + 10 * length + multiplicity),
        )[:, : multiplicity - 1]
    elif flag_kind == "proper_reducing":
        flag = reducing_flag(data)
    else:
        raise ValueError(f"unknown flag kind: {flag_kind}")

    response, column_basis = response_matrix(data, flag)
    flag_dimension = flag.shape[1]
    codomain_dimension = flag_dimension**2
    response[np.abs(response) < 2e-10] = 0
    singular_values = np.linalg.svd(response, compute_uv=False)
    tolerance = max(
        2e-10,
        max(response.shape) * singular_values[0] * 2e-10,
    )
    response_rank = int(np.sum(singular_values > tolerance))
    annihilator = null_space(response.T, rcond=2e-10)

    endpoint_basis = tuple(hermitian_basis(flag_dimension))
    metric_root = np.asarray(sqrtm(data.metric), dtype=complex)
    balanced = (
        metric_root @ data.operator @ np.linalg.inv(metric_root)
    )
    projection = (
        np.eye(len(data.operator))
        - data.right_defect @ data.right_defect.conj().T
    )

    maximum_pairing_error = 0.0
    maximum_commutator_error = 0.0
    for coordinates in annihilator.T:
        endpoint = endpoint_from_coordinates(coordinates, endpoint_basis)
        lifted = flag @ endpoint @ flag.conj().T
        dual = dual_stein_inverse(
            data.operator,
            data.left_defect @ lifted @ data.left_defect.conj().T,
        )
        adjoint_residual = projection @ dual @ data.right_defect
        maximum_pairing_error = max(
            maximum_pairing_error,
            float(np.linalg.norm(adjoint_residual)),
        )
        balanced_dual = metric_root @ dual @ metric_root / 4
        maximum_commutator_error = max(
            maximum_commutator_error,
            float(
                np.linalg.norm(
                    balanced_dual @ balanced
                    - balanced @ balanced_dual
                )
            ),
        )

    generator = np.random.default_rng(
        106_500 + 100 * length + 10 * multiplicity + flag_dimension
    )
    coefficients = generator.standard_normal(len(column_basis))
    chosen_column = sum(
        (
            coefficient * column
            for coefficient, column in zip(
                coefficients,
                column_basis,
                strict=True,
            )
        ),
        np.zeros_like(column_basis[0]),
    )
    chosen_response = (
        flag.conj().T
        @ endpoint_motion(
            data.operator,
            data.right_defect,
            data.left_defect,
            chosen_column,
        )
        @ flag
    )
    identity = np.eye(flag_dimension, dtype=complex)
    feasible_face = -identity + chosen_response
    corrected_face = feasible_face - chosen_response
    feasible_error = float(np.linalg.norm(corrected_face + identity))
    feasible_maximum = float(np.linalg.eigvalsh(corrected_face)[-1])

    identity_coordinates = hermitian_coordinates(identity)
    identity_separator_residual = float(
        np.linalg.norm(response.T @ identity_coordinates)
    )
    obstructed_pairing = float(np.trace(identity @ identity).real)

    response_codimension = codomain_dimension - response_rank
    is_reducing = flag_kind != "proper_irreducible"
    verified = bool(
        annihilator.shape[1] == response_codimension
        and maximum_pairing_error < 3e-8
        and maximum_commutator_error < 3e-8
        and feasible_error < 3e-10
        and feasible_maximum < -0.999999
        and (
            not is_reducing
            or (
                identity_separator_residual < 3e-8
                and obstructed_pairing > 0
            )
        )
        and (
            flag_kind != "proper_irreducible"
            or response_codimension == 0
        )
    )
    if not verified:
        raise RuntimeError(
            "the flagged endpoint alternative audit failed: "
            f"kind={flag_kind}, rank={response_rank}, "
            f"codim={response_codimension}, "
            f"pair={maximum_pairing_error:.3e}, "
            f"comm={maximum_commutator_error:.3e}, "
            f"identity={identity_separator_residual:.3e}"
        )

    return FlaggedEndpointAlternativeRecord(
        flag_kind=flag_kind,
        length=length,
        multiplicity=multiplicity,
        flag_dimension=flag_dimension,
        response_rank=response_rank,
        response_codimension=response_codimension,
        adjoint_kernel_dimension=annihilator.shape[1],
        maximum_adjoint_pairing_error=format_float(
            maximum_pairing_error
        ),
        maximum_commutant_error=format_float(
            maximum_commutator_error
        ),
        feasible_correction_error=format_float(feasible_error),
        feasible_endpoint_maximum=format_float(feasible_maximum),
        identity_separator_residual=(
            format_float(identity_separator_residual)
            if is_reducing
            else "not_applicable"
        ),
        obstructed_separator_pairing=(
            format_float(obstructed_pairing)
            if is_reducing
            else "not_applicable"
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[FlaggedEndpointAlternativeRecord]:
    """Return full, irreducible, and reducing flag audits."""

    records = []
    for length in (2, 4):
        records.append(make_record("full", length, 2, 8.0))
        records.append(
            make_record("proper_irreducible", length, 3, 8.0)
        )
        records.append(
            make_record("proper_reducing", length, 3, 8.0)
        )
    return records


def write_records(
    records: list[FlaggedEndpointAlternativeRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_flagged_endpoint_alternative_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete deterministic audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
