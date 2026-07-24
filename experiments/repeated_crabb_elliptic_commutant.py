#!/usr/bin/env python3
"""Audit the commutant description of the L204 elliptic cokernel.

For a Hermitian left-defect matrix ``Y``, let ``Z_Y`` solve

    Z_Y - T Z_Y T* = W Y W*.

L204 says that ``Y`` lies in the endpoint cokernel exactly when
``(I-VV*) Z_Y V = 0``.  After balancing, the normalized Gramian
``H_Y`` should then commute with the partial isometry ``S`` and
intertwine the two defect actions:

    H_Y W = W Y,  H_Y V = V A,
    Y B_n = B_n A.

The checker finds the complete real cokernel, verifies these identities,
and checks the resulting weighted compatibility with L204's target.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space, sqrtm

from crabb_block_hardy_equality import (
    format_float,
    hermitian_basis,
)
from repeated_crabb_elliptic_cokernel import (
    build_reduced_face_data,
    dual_stein_inverse,
)


@dataclass(frozen=True)
class EllipticCommutantRecord:
    """One complete cokernel-to-commutant audit."""

    length: int
    multiplicity: int
    toeplitz_strength: str
    endpoint_cokernel_dimension: int
    expected_cokernel_dimension: int
    cokernel_residual: str
    state_commutator_error: str
    left_defect_intertwining_error: str
    right_defect_intertwining_error: str
    transfer_intertwining_error: str
    weighted_l203_identity_error: str
    weighted_target_compatibility_error: str
    all_checks_passed: bool


def real_matrix_coordinates(matrix: np.ndarray) -> np.ndarray:
    """Return all real and imaginary entries of a complex matrix."""

    return np.concatenate((matrix.real.ravel(), matrix.imag.ravel()))


def make_record(
    length: int,
    multiplicity: int,
    strength: float,
) -> EllipticCommutantRecord:
    """Audit the complete Hermitian cokernel at one equality anchor."""

    data = build_reduced_face_data(
        length,
        multiplicity,
        strength,
    )
    operator = data.operator
    metric = data.metric
    right = data.right_defect
    left = data.left_defect
    projection = np.eye(len(operator)) - right @ right.conj().T
    endpoint_basis = hermitian_basis(multiplicity)

    dual_metrics = [
        dual_stein_inverse(
            operator,
            left @ endpoint @ left.conj().T,
        )
        for endpoint in endpoint_basis
    ]
    adjoint_columns = np.stack(
        [
            real_matrix_coordinates(
                projection @ dual @ right
            )
            for dual in dual_metrics
        ],
        axis=1,
    )
    cokernel_coordinates = null_space(
        adjoint_columns,
        rcond=1e-8,
    )
    cokernel_dimension = cokernel_coordinates.shape[1]

    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    balanced = (
        metric_root
        @ operator
        @ np.linalg.inv(metric_root)
    )
    coefficients = [
        left.conj().T
        @ np.linalg.matrix_power(balanced.conj().T, degree)
        @ right
        for degree in range(1, 2 * length + 9)
    ]

    cokernel_error = 0.0
    commutator_error = 0.0
    left_error = 0.0
    right_error = 0.0
    transfer_error = 0.0
    weighted_identity_error = 0.0
    compatibility_error = 0.0
    base_forcing = (
        data.fixed_forcing
        + right @ data.parallel_column.conj().T
        + data.parallel_column @ right.conj().T
    )
    for column in cokernel_coordinates.T:
        endpoint = sum(
            coefficient * basis
            for coefficient, basis in zip(
                column,
                endpoint_basis,
                strict=True,
            )
        )
        dual = sum(
            coefficient * basis
            for coefficient, basis in zip(
                column,
                dual_metrics,
                strict=True,
            )
        )
        balanced_dual = metric_root @ dual @ metric_root / 4
        right_action = (
            right.conj().T @ balanced_dual @ right
        )

        cokernel_error = max(
            cokernel_error,
            float(
                np.linalg.norm(
                    projection @ dual @ right
                )
            ),
        )
        commutator_error = max(
            commutator_error,
            float(
                np.linalg.norm(
                    balanced_dual @ balanced
                    - balanced @ balanced_dual
                )
            ),
        )
        left_error = max(
            left_error,
            float(
                np.linalg.norm(
                    balanced_dual @ left - left @ endpoint
                )
            ),
        )
        right_error = max(
            right_error,
            float(
                np.linalg.norm(
                    balanced_dual @ right
                    - right @ right_action
                )
            ),
        )
        for coefficient in coefficients:
            transfer_error = max(
                transfer_error,
                float(
                    np.linalg.norm(
                        endpoint @ coefficient
                        - coefficient @ right_action
                    )
                ),
            )
        weighted_identity_error = max(
            weighted_identity_error,
            float(
                abs(
                    np.trace(dual @ base_forcing)
                    + np.trace(
                        endpoint
                        @ data.normal_corner
                        @ data.normal_corner.conj().T
                    )
                )
            ),
        )
        compatibility_error = max(
            compatibility_error,
            float(
                abs(np.trace(endpoint @ data.target_gap))
            ),
        )

    # The deterministic m=2 line is irreducible; the m=3 generator
    # retains one proper reducing copy block.
    expected_dimension = 1 if multiplicity == 2 else 2
    verified = bool(
        cokernel_dimension == expected_dimension
        and cokernel_error < 2e-9
        and commutator_error < 2e-9
        and left_error < 2e-9
        and right_error < 2e-9
        and transfer_error < 2e-9
        and weighted_identity_error < 2e-9
        and compatibility_error < 2e-9
    )
    if not verified:
        raise RuntimeError("the elliptic commutant audit failed")

    return EllipticCommutantRecord(
        length=length,
        multiplicity=multiplicity,
        toeplitz_strength=format_float(data.actual_strength),
        endpoint_cokernel_dimension=cokernel_dimension,
        expected_cokernel_dimension=expected_dimension,
        cokernel_residual=format_float(cokernel_error),
        state_commutator_error=format_float(commutator_error),
        left_defect_intertwining_error=format_float(left_error),
        right_defect_intertwining_error=format_float(right_error),
        transfer_intertwining_error=format_float(transfer_error),
        weighted_l203_identity_error=format_float(
            weighted_identity_error
        ),
        weighted_target_compatibility_error=format_float(
            compatibility_error
        ),
        all_checks_passed=verified,
    )


def write_records(
    records: list[EllipticCommutantRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines records atomically."""

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
            "repeated_crabb_elliptic_commutant_s70224.jsonl"
        ),
    )
    parser.add_argument(
        "--strengths",
        type=float,
        nargs="+",
        default=(0.1, 8.0),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete cokernel/commutant audit."""

    args = parse_args()
    records = [
        make_record(length, multiplicity, strength)
        for multiplicity in (2, 3)
        for length in range(2, 6)
        for strength in args.strengths
    ]
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
