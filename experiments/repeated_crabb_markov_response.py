#!/usr/bin/env python3
"""Audit the quantum-Markov factorization of the endpoint response."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import (
    gauged_shift,
    generalized_grade_column,
    stable_random_partial_isometry,
)
from repeated_crabb_elliptic_cokernel import (
    dual_stein_inverse,
    endpoint_motion,
)
from repeated_crabb_transfer_flag import (
    transfer_channel,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class MarkovResponseRecord:
    """One polarized-response and Markov-Laplacian audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    maximum_grade: int
    spectral_radius: str
    maximum_polarized_response_error: str
    left_unital_error: str
    right_unital_error: str
    markov_response_error: str
    dirichlet_defect_error: str
    endpoint_column_norm: str
    all_checks_passed: bool


def hermitian_part(matrix: Matrix) -> Matrix:
    """Return the Hermitian part of a square matrix."""

    return (matrix + matrix.conj().T) / 2


def transfer_channel_adjoint(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    copy_matrix: Matrix,
) -> Matrix:
    """Apply the Hilbert--Schmidt adjoint transfer channel."""

    gramian = dual_stein_inverse(
        operator,
        left @ copy_matrix @ left.conj().T,
    )
    return right.conj().T @ gramian @ right


def physical_data(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, Matrix]:
    """Return the physical operator and equality metric root."""

    identity = np.eye(len(operator), dtype=complex)
    metric = (
        2 * identity
        - right @ right.conj().T
        + 2 * left @ left.conj().T
    )
    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    physical = np.linalg.inv(metric_root) @ operator @ metric_root
    return physical, metric_root


def random_hermitian(
    size: int,
    generator: np.random.Generator,
) -> Matrix:
    """Return a normalized deterministic random Hermitian matrix."""

    raw = (
        generator.standard_normal((size, size))
        + 1j * generator.standard_normal((size, size))
    )
    matrix = hermitian_part(raw)
    return matrix / np.linalg.norm(matrix)


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    maximum_grade: int,
    seed: int,
) -> MarkovResponseRecord:
    """Audit one structured or unstructured colligation."""

    multiplicity = right.shape[1]
    generator = np.random.default_rng(seed)
    physical, metric_root = physical_data(operator, right, left)

    maximum_polarized_error = 0.0
    for grade in range(1, min(maximum_grade, 8) + 1):
        coefficient = transfer_coefficient(
            operator,
            right,
            left,
            grade,
        )
        multiplier = (
            generator.standard_normal((multiplicity, multiplicity))
            + 1j
            * generator.standard_normal((multiplicity, multiplicity))
        )
        balanced_column, _ = generalized_grade_column(
            operator,
            right,
            left,
            grade,
            multiplier,
        )
        left_face = hermitian_part(
            coefficient @ multiplier.conj().T
        )
        right_face = hermitian_part(
            multiplier.conj().T @ coefficient
        )
        expected = 8 * (
            left_face
            - transfer_channel(
                operator,
                right,
                left,
                right_face,
            )
        )
        actual = endpoint_motion(
            physical,
            right,
            left,
            metric_root @ balanced_column,
        )
        maximum_polarized_error = max(
            maximum_polarized_error,
            float(np.linalg.norm(actual - expected)),
        )

    identity = np.eye(multiplicity, dtype=complex)
    left_unital_error = float(
        np.linalg.norm(
            transfer_channel(
                operator,
                right,
                left,
                identity,
            )
            - identity
        )
    )
    right_unital_error = float(
        np.linalg.norm(
            transfer_channel_adjoint(
                operator,
                right,
                left,
                identity,
            )
            - identity
        )
    )

    endpoint = random_hermitian(multiplicity, generator)
    channel_adjoint_endpoint = transfer_channel_adjoint(
        operator,
        right,
        left,
        endpoint,
    )
    markov_endpoint = transfer_channel(
        operator,
        right,
        left,
        channel_adjoint_endpoint,
    )

    aggregate_column = np.zeros(
        (len(operator), multiplicity),
        dtype=complex,
    )
    dirichlet_defect = 0.0
    for grade in range(1, maximum_grade + 1):
        coefficient = transfer_coefficient(
            operator,
            right,
            left,
            grade,
        )
        balanced_column, _ = generalized_grade_column(
            operator,
            right,
            left,
            grade,
            endpoint @ coefficient,
        )
        aggregate_column += balanced_column
        intertwining_defect = (
            endpoint @ coefficient
            - coefficient @ channel_adjoint_endpoint
        )
        dirichlet_defect += float(
            np.linalg.norm(intertwining_defect) ** 2
        )

    actual_markov_response = endpoint_motion(
        physical,
        right,
        left,
        metric_root @ aggregate_column,
    )
    expected_markov_response = 8 * (
        endpoint - markov_endpoint
    )
    markov_error = float(
        np.linalg.norm(
            actual_markov_response - expected_markov_response
        )
    )
    dirichlet = float(
        np.trace(
            endpoint @ (endpoint - markov_endpoint)
        ).real
    )
    dirichlet_error = abs(dirichlet_defect - dirichlet)
    column_norm = float(np.linalg.norm(aggregate_column))

    spectral_radius = float(
        np.max(np.abs(np.linalg.eigvals(operator)))
    )
    tolerance = 2e-8
    is_apex = construction_kind == "repeated_monomial_shift"
    verified = bool(
        spectral_radius < 1
        and maximum_polarized_error < tolerance
        and left_unital_error < tolerance
        and right_unital_error < tolerance
        and markov_error < tolerance
        and dirichlet_error < tolerance
        and (not is_apex or column_norm < tolerance)
    )
    if not verified:
        raise RuntimeError(
            "Markov response audit failed: "
            f"kind={construction_kind}, "
            f"polarized={maximum_polarized_error:.3e}, "
            f"left={left_unital_error:.3e}, "
            f"right={right_unital_error:.3e}, "
            f"markov={markov_error:.3e}, "
            f"dirichlet={dirichlet_error:.3e}, "
            f"column={column_norm:.3e}"
        )

    return MarkovResponseRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=multiplicity,
        maximum_grade=maximum_grade,
        spectral_radius=format_float(spectral_radius),
        maximum_polarized_response_error=format_float(
            maximum_polarized_error
        ),
        left_unital_error=format_float(left_unital_error),
        right_unital_error=format_float(right_unital_error),
        markov_response_error=format_float(markov_error),
        dirichlet_defect_error=format_float(dirichlet_error),
        endpoint_column_norm=format_float(column_norm),
        all_checks_passed=verified,
    )


def standard_records() -> list[MarkovResponseRecord]:
    """Return deterministic structured and unstructured audits."""

    generator = np.random.default_rng(280_001)
    records: list[MarkovResponseRecord] = []
    for index, (dimension, multiplicity) in enumerate(
        ((7, 2), (9, 2), (11, 3)),
    ):
        while True:
            data = stable_random_partial_isometry(
                dimension,
                multiplicity,
                generator,
            )
            if np.max(np.abs(np.linalg.eigvals(data[0]))) < 0.78:
                break
        records.append(
            audit_case(
                "unstructured_partial_isometry",
                *data,
                maximum_grade=100,
                seed=280_100 + index,
            )
        )

    apex = gauged_shift((4, 4, 4), 280_200)
    records.append(
        audit_case(
            "repeated_monomial_shift",
            *apex,
            maximum_grade=12,
            seed=280_201,
        )
    )
    return records


def write_records(
    records: list[MarkovResponseRecord],
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
            "repeated_crabb_markov_response_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist all audits."""

    records = standard_records()
    digest = write_records(records, parse_args().output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
