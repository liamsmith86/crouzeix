#!/usr/bin/env python3
"""Audit the candidate second elliptic face at repeated Crabb equality.

L202 identifies the sole partial-isometry-normal block in the real
ellipse tangent as

    Z = W* dot(S) V.

The correctly oriented upper endpoint acts on the left defect space, so
its candidate second coefficient is ``-Z Z*``.  This is not the same
matrix as the right transfer Gram ``-16 B_1* B_1`` when the copy
coefficients do not commute, although the two matrices have the same
eigenvalues.

This checker builds the complete second-order Stein equation, imposes a
zero lower endpoint and the oriented upper endpoint ``-Z Z*``, and solves
the resulting real linear system.  It also verifies the exact scalar
trace-face identity.  Linear solvability is evidence for the matrix face;
it is not an analytic all-anchor lifting theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import sqrtm

from crabb_block_hardy_equality import (
    format_float,
    hermitian_basis,
    pack_hermitian,
)
from repeated_crabb_elliptic_first_jet import physical_hardy_data
from repeated_crabb_inner_faber_transfer import (
    canonical_transfer_data,
    strengthened_inverse_toeplitz,
    transfer_coefficients,
)


@dataclass(frozen=True)
class EllipticSecondFaceRecord:
    """One oriented second-face linear audit."""

    length: int
    multiplicity: int
    toeplitz_strength: str
    coefficient_commutator_norm: str
    normal_corner_norm: str
    right_transfer_gram_error: str
    trace_face_error: str
    second_stein_equation_error: str
    lower_endpoint_error: str
    oriented_upper_endpoint_error: str
    solution_norm: str
    linear_system_rank: int
    linear_system_row_count: int
    all_checks_passed: bool


def hermitian_coordinates(matrix: np.ndarray) -> np.ndarray:
    """Pack a numerically Hermitian matrix into real coordinates."""

    hermitian = (matrix + matrix.conj().T) / 2
    return pack_hermitian(hermitian)


def second_ellipse_coefficients(
    operator: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the first two coefficients of the exact ellipse pullback.

    The direct ellipse map has expansion

        phi_c(w) = w - c w^3 + c^2(2w + w^5) + O(c^3).

    Substitution of ``w=T+cT*`` gives the two matrices below.
    """

    adjoint = operator.conj().T
    first = adjoint - np.linalg.matrix_power(operator, 3)
    second = (
        2 * operator
        - operator @ operator @ adjoint
        - operator @ adjoint @ operator
        - adjoint @ operator @ operator
        + np.linalg.matrix_power(operator, 5)
    )
    return first, second


def canonical_first_forcing_motion(
    operator: np.ndarray,
    metric: np.ndarray,
    right_defect: np.ndarray,
    first: np.ndarray,
) -> np.ndarray:
    """Factor the first Stein derivative in the canonical defect gauge."""

    derivative = -(
        first.conj().T @ metric @ operator
        + operator.conj().T @ metric @ first
    )
    projection = right_defect @ right_defect.conj().T
    corner = right_defect.conj().T @ derivative @ right_defect
    return (
        (np.eye(len(operator)) - projection)
        @ derivative
        @ right_defect
        + right_defect @ corner / 2
    )


def second_forcing(
    operator: np.ndarray,
    metric: np.ndarray,
    first: np.ndarray,
    second: np.ndarray,
    first_forcing: np.ndarray,
) -> np.ndarray:
    """Return the fixed part of the second Stein equation."""

    return (
        first_forcing @ first_forcing.conj().T
        + second.conj().T @ metric @ operator
        + operator.conj().T @ metric @ second
        + first.conj().T @ metric @ first
    )


def solve_oriented_second_face(
    operator: np.ndarray,
    fixed_forcing: np.ndarray,
    right_defect: np.ndarray,
    left_defect: np.ndarray,
    target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, int, int]:
    """Solve the complete real-linear second Stein and endpoint system."""

    dimension = len(operator)
    multiplicity = right_defect.shape[1]
    metric_basis = hermitian_basis(dimension)
    columns: list[np.ndarray] = []

    for direction in metric_basis:
        stein = (
            direction
            - operator.conj().T @ direction @ operator
        )
        columns.append(
            np.concatenate(
                (
                    hermitian_coordinates(stein),
                    hermitian_coordinates(
                        right_defect.conj().T
                        @ direction
                        @ right_defect
                    ),
                    hermitian_coordinates(
                        left_defect.conj().T
                        @ direction
                        @ left_defect
                    ),
                )
            )
        )

    forcing_directions: list[np.ndarray] = []
    for row in range(dimension):
        for column in range(multiplicity):
            for phase in (1 + 0j, 1j):
                direction = np.zeros(
                    (dimension, multiplicity),
                    dtype=complex,
                )
                direction[row, column] = phase
                forcing_directions.append(direction)
                stein = -(
                    right_defect @ direction.conj().T
                    + direction @ right_defect.conj().T
                )
                columns.append(
                    np.concatenate(
                        (
                            hermitian_coordinates(stein),
                            np.zeros(multiplicity**2),
                            np.zeros(multiplicity**2),
                        )
                    )
                )

    system = np.stack(columns, axis=1)
    right_hand_side = np.concatenate(
        (
            hermitian_coordinates(fixed_forcing),
            np.zeros(multiplicity**2),
            hermitian_coordinates(target),
        )
    )
    solution, _, rank, _ = np.linalg.lstsq(
        system,
        right_hand_side,
        rcond=1e-11,
    )

    metric_count = len(metric_basis)
    metric_direction = sum(
        coefficient * basis
        for coefficient, basis in zip(
            solution[:metric_count],
            metric_basis,
            strict=True,
        )
    )
    second_forcing_motion = np.zeros(
        (dimension, multiplicity),
        dtype=complex,
    )
    for coefficient, direction in zip(
        solution[metric_count:],
        forcing_directions,
        strict=True,
    ):
        second_forcing_motion += coefficient * direction

    residual = float(
        np.linalg.norm(system @ solution - right_hand_side)
    )
    return (
        metric_direction,
        second_forcing_motion,
        residual,
        int(rank),
        int(system.shape[0]),
    )


def make_record(
    length: int,
    multiplicity: int,
    strength: float,
) -> EllipticSecondFaceRecord:
    """Build one noncommutative oriented second-face record."""

    (
        operator,
        metric,
        right,
        left,
        _,
        commutator,
        actual_strength,
    ) = physical_hardy_data(length, multiplicity, strength)
    first, second = second_ellipse_coefficients(operator)
    first_forcing = canonical_first_forcing_motion(
        operator,
        metric,
        right,
        first,
    )
    fixed_forcing = second_forcing(
        operator,
        metric,
        first,
        second,
        first_forcing,
    )

    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    metric_root_inverse = np.linalg.inv(metric_root)
    balanced_first = metric_root @ first @ metric_root_inverse
    normal_corner = left.conj().T @ balanced_first @ right
    oriented_target = -normal_corner @ normal_corner.conj().T

    (
        metric_direction,
        second_forcing_motion,
        equation_error,
        system_rank,
        system_rows,
    ) = solve_oriented_second_face(
        operator,
        fixed_forcing,
        right,
        left,
        oriented_target,
    )

    reconstructed = (
        metric_direction
        - operator.conj().T @ metric_direction @ operator
        - right @ second_forcing_motion.conj().T
        - second_forcing_motion @ right.conj().T
    )
    second_stein_error = float(
        np.linalg.norm(reconstructed - fixed_forcing)
    )
    lower_error = float(
        np.linalg.norm(
            right.conj().T @ metric_direction @ right
        )
    )
    upper_error = float(
        np.linalg.norm(
            left.conj().T @ metric_direction @ left
            - oriented_target
        )
    )

    inverse, _ = strengthened_inverse_toeplitz(
        length,
        multiplicity,
        actual_strength,
    )
    transfer = canonical_transfer_data(
        np.linalg.inv(inverse),
        length,
        multiplicity,
    )
    first_transfer = transfer_coefficients(transfer, 2)[1]
    transfer_error = float(
        np.linalg.norm(
            np.linalg.eigvalsh(
                normal_corner.conj().T @ normal_corner
            )
            - np.linalg.eigvalsh(
                16 * first_transfer.conj().T @ first_transfer
            )
        )
    )

    dual_tangent = 4 * (
        np.linalg.inv(metric) - right @ right.conj().T
    )
    trace_face_error = float(
        abs(
            np.trace(dual_tangent @ fixed_forcing)
            + np.linalg.norm(normal_corner) ** 2
        )
    )
    solution_norm = float(
        np.hypot(
            np.linalg.norm(metric_direction),
            np.linalg.norm(second_forcing_motion),
        )
    )

    verified = bool(
        commutator > 1e-8
        and np.linalg.norm(normal_corner) > 1e-8
        and transfer_error < 3e-10
        and trace_face_error < 3e-10
        and equation_error < 3e-10
        and second_stein_error < 3e-10
        and lower_error < 3e-10
        and upper_error < 3e-10
        and solution_norm < 100
    )
    if not verified:
        raise RuntimeError("the oriented elliptic second-face audit failed")
    return EllipticSecondFaceRecord(
        length=length,
        multiplicity=multiplicity,
        toeplitz_strength=format_float(actual_strength),
        coefficient_commutator_norm=format_float(commutator),
        normal_corner_norm=format_float(np.linalg.norm(normal_corner)),
        right_transfer_gram_error=format_float(transfer_error),
        trace_face_error=format_float(trace_face_error),
        second_stein_equation_error=format_float(second_stein_error),
        lower_endpoint_error=format_float(lower_error),
        oriented_upper_endpoint_error=format_float(upper_error),
        solution_norm=format_float(solution_norm),
        linear_system_rank=system_rank,
        linear_system_row_count=system_rows,
        all_checks_passed=verified,
    )


def write_records(
    records: list[EllipticSecondFaceRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines records."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/repeated_crabb_elliptic_second_face_s70224.jsonl"
        ),
    )
    parser.add_argument("--strength", type=float, default=18.0)
    return parser.parse_args()


def main() -> None:
    """Run the standard noncommutative second-face audit."""

    args = parse_args()
    records = [
        make_record(length, multiplicity, args.strength)
        for multiplicity in (2, 3)
        for length in range(2, 6)
    ]
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
