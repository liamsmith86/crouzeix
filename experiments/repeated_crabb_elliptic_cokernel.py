#!/usr/bin/env python3
"""Reduce the repeated elliptic second face to one copy-space range equation.

L203 writes the second Stein equation as

    X - T* X T = F_2 + V C_2* + C_2 V*.

Because ``T V = 0``, the part of ``C_2`` parallel to ``V`` can be
eliminated explicitly to impose ``V* X V = 0``.  The only remaining
choice is a column ``C`` perpendicular to ``V``.  Its effect on the
left-defect endpoint is the real-linear map

    M_T(C) = W* G_T(V C* + C V*) W,

where ``G_T`` is the stable Stein inverse.  This checker verifies that
reduction, its exact adjoint, and the observed bounded divisibility of
the L203 target as a noncommuting equality anchor approaches Crabb.

The range solves are numerical evidence.  The reduction and adjoint
formulas are finite-dimensional identities.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_discrete_lyapunov, sqrtm

from crabb_block_hardy_equality import (
    format_float,
    hermitian_basis,
)
from repeated_crabb_elliptic_first_jet import physical_hardy_data
from repeated_crabb_elliptic_second_face import (
    canonical_first_forcing_motion,
    hermitian_coordinates,
    second_ellipse_coefficients,
    second_forcing,
)


@dataclass(frozen=True)
class EllipticCokernelRecord:
    """One reduced second-face range and adjoint audit."""

    length: int
    multiplicity: int
    requested_toeplitz_strength: str
    actual_toeplitz_strength: str
    coefficient_commutator_norm: str
    endpoint_map_rank: int
    endpoint_space_dimension: int
    endpoint_cokernel_dimension: int
    target_gap_norm: str
    target_gap_over_strength_squared: str
    correction_norm: str
    correction_over_strength: str
    target_trace_error: str
    range_equation_error: str
    adjoint_identity_error: str
    full_second_stein_error: str
    lower_endpoint_error: str
    upper_endpoint_error: str
    all_checks_passed: bool


@dataclass(frozen=True)
class ReducedFaceData:
    """Matrices defining one reduced L204 endpoint equation."""

    operator: np.ndarray
    metric: np.ndarray
    right_defect: np.ndarray
    left_defect: np.ndarray
    fixed_forcing: np.ndarray
    parallel_column: np.ndarray
    base_metric_direction: np.ndarray
    normal_corner: np.ndarray
    target_endpoint: np.ndarray
    target_gap: np.ndarray
    commutator_norm: float
    actual_strength: float


def checked_discrete_lyapunov(
    coefficient: np.ndarray,
    forcing: np.ndarray,
) -> np.ndarray:
    """Solve one discrete Lyapunov equation with a residual guard.

    SciPy's automatic bilinear branch can be inaccurate for some
    nonnormal real matrices of dimension at least ten.  Retrying the
    same equation in complex arithmetic avoids that implementation
    corner without making every larger solve use the much more
    expensive direct Kronecker method.
    """

    coefficient = np.asarray(coefficient)
    forcing = np.asarray(forcing)

    def solve(
        matrix: np.ndarray,
        right_hand_side: np.ndarray,
    ) -> np.ndarray:
        result = solve_discrete_lyapunov(matrix, right_hand_side)
        return (result + result.conj().T) / 2

    solution = solve(coefficient, forcing)
    residual = (
        coefficient @ solution @ coefficient.conj().T
        - solution
        + forcing
    )
    scale = (
        1
        + np.linalg.norm(solution)
        + np.linalg.norm(forcing)
    )
    if np.linalg.norm(residual) > 1e-10 * scale:
        solution = solve(
            np.asarray(coefficient, dtype=complex),
            np.asarray(forcing, dtype=complex),
        )
        residual = (
            coefficient @ solution @ coefficient.conj().T
            - solution
            + forcing
        )
    if np.linalg.norm(residual) > 1e-9 * scale:
        raise RuntimeError(
            "the discrete Lyapunov solve failed its residual guard"
        )
    return solution


def stein_inverse(operator: np.ndarray, forcing: np.ndarray) -> np.ndarray:
    """Solve ``X - T* X T = forcing`` and remove roundoff skew."""

    return checked_discrete_lyapunov(operator.conj().T, forcing)


def dual_stein_inverse(
    operator: np.ndarray,
    forcing: np.ndarray,
) -> np.ndarray:
    """Solve ``Z - T Z T* = forcing``."""

    return checked_discrete_lyapunov(operator, forcing)


def projected_column_basis(
    right_defect: np.ndarray,
) -> tuple[np.ndarray, ...]:
    """Return a Parseval real frame for columns perpendicular to ``V``."""

    dimension, multiplicity = right_defect.shape
    projection = (
        np.eye(dimension)
        - right_defect @ right_defect.conj().T
    )
    basis: list[np.ndarray] = []
    for row in range(dimension):
        for column in range(multiplicity):
            unit = np.zeros((dimension, multiplicity), dtype=complex)
            unit[row, column] = 1
            projected = projection @ unit
            basis.extend((projected, 1j * projected))
    return tuple(basis)


def endpoint_motion(
    operator: np.ndarray,
    right_defect: np.ndarray,
    left_defect: np.ndarray,
    column: np.ndarray,
) -> np.ndarray:
    """Apply the reduced endpoint map ``M_T``."""

    forcing = (
        right_defect @ column.conj().T
        + column @ right_defect.conj().T
    )
    metric = stein_inverse(operator, forcing)
    endpoint = left_defect.conj().T @ metric @ left_defect
    return (endpoint + endpoint.conj().T) / 2


def build_reduced_face_data(
    length: int,
    multiplicity: int,
    requested_strength: float,
) -> ReducedFaceData:
    """Construct the matrices in the reduced L204 range equation."""

    (
        operator,
        metric,
        right,
        left,
        _,
        commutator,
        actual_strength,
    ) = physical_hardy_data(
        length,
        multiplicity,
        requested_strength,
    )
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

    # The parallel part of C_2 is forced by the lower endpoint.
    lower_corner = (
        right.conj().T @ fixed_forcing @ right
    )
    lower_corner = (lower_corner + lower_corner.conj().T) / 2
    parallel_column = -right @ lower_corner / 2
    base_forcing = (
        fixed_forcing
        + right @ parallel_column.conj().T
        + parallel_column @ right.conj().T
    )
    base_metric = stein_inverse(operator, base_forcing)

    metric_root = np.asarray(sqrtm(metric), dtype=complex)
    balanced_first = (
        metric_root
        @ first
        @ np.linalg.inv(metric_root)
    )
    normal_corner = left.conj().T @ balanced_first @ right
    target = -normal_corner @ normal_corner.conj().T
    base_endpoint = left.conj().T @ base_metric @ left
    target_gap = target - base_endpoint
    target_gap = (target_gap + target_gap.conj().T) / 2

    return ReducedFaceData(
        operator=operator,
        metric=metric,
        right_defect=right,
        left_defect=left,
        fixed_forcing=fixed_forcing,
        parallel_column=parallel_column,
        base_metric_direction=base_metric,
        normal_corner=normal_corner,
        target_endpoint=target,
        target_gap=target_gap,
        commutator_norm=commutator,
        actual_strength=actual_strength,
    )


def reduced_second_face(
    length: int,
    multiplicity: int,
    requested_strength: float,
) -> EllipticCokernelRecord:
    """Build one reduced range equation and its minimum-frame solution."""

    data = build_reduced_face_data(
        length,
        multiplicity,
        requested_strength,
    )
    operator = data.operator
    right = data.right_defect
    left = data.left_defect
    fixed_forcing = data.fixed_forcing
    parallel_column = data.parallel_column
    base_metric = data.base_metric_direction
    target = data.target_endpoint
    target_gap = data.target_gap
    commutator = data.commutator_norm
    actual_strength = data.actual_strength

    column_basis = projected_column_basis(right)
    map_columns = np.stack(
        [
            hermitian_coordinates(
                endpoint_motion(
                    operator,
                    right,
                    left,
                    column,
                )
            )
            for column in column_basis
        ],
        axis=1,
    )
    target_coordinates = hermitian_coordinates(target_gap)
    solution, _, rank, _ = np.linalg.lstsq(
        map_columns,
        target_coordinates,
        rcond=1e-10,
    )
    free_column = sum(
        coefficient * column
        for coefficient, column in zip(
            solution,
            column_basis,
            strict=True,
        )
    )
    range_error = float(
        np.linalg.norm(
            map_columns @ solution - target_coordinates
        )
    )

    free_forcing = (
        right @ free_column.conj().T
        + free_column @ right.conj().T
    )
    free_metric = stein_inverse(operator, free_forcing)
    metric_direction = base_metric + free_metric
    second_column = parallel_column + free_column

    reconstructed = (
        metric_direction
        - operator.conj().T @ metric_direction @ operator
    )
    expected = (
        fixed_forcing
        + right @ second_column.conj().T
        + second_column @ right.conj().T
    )
    full_stein_error = float(np.linalg.norm(reconstructed - expected))
    lower_error = float(
        np.linalg.norm(
            right.conj().T @ metric_direction @ right
        )
    )
    upper_error = float(
        np.linalg.norm(
            left.conj().T @ metric_direction @ left - target
        )
    )

    # Audit <Y,M(C)> = 2 Re tr(C* Z_Y V) on complete real bases.
    adjoint_error = 0.0
    for endpoint_test in hermitian_basis(multiplicity):
        dual_metric = dual_stein_inverse(
            operator,
            left @ endpoint_test @ left.conj().T,
        )
        for column in column_basis:
            primal_value = float(
                np.trace(
                    endpoint_test
                    @ endpoint_motion(
                        operator,
                        right,
                        left,
                        column,
                    )
                ).real
            )
            adjoint_value = float(
                2
                * np.trace(
                    column.conj().T
                    @ dual_metric
                    @ right
                ).real
            )
            adjoint_error = max(
                adjoint_error,
                abs(primal_value - adjoint_value),
            )

    strength = float(actual_strength)
    gap_norm = float(np.linalg.norm(target_gap))
    correction_norm = float(np.linalg.norm(free_column))
    trace_error = float(abs(np.trace(target_gap)))
    endpoint_dimension = multiplicity**2
    cokernel_dimension = endpoint_dimension - int(rank)
    verified = bool(
        commutator > 1e-12
        and strength > 0
        and gap_norm / strength**2 < 0.02
        and correction_norm / strength < 0.08
        and trace_error < 5e-10
        and range_error < 5e-10
        and adjoint_error < 5e-10
        and full_stein_error < 5e-10
        and lower_error < 5e-10
        and upper_error < 5e-10
    )
    if not verified:
        raise RuntimeError("the reduced elliptic cokernel audit failed")

    return EllipticCokernelRecord(
        length=length,
        multiplicity=multiplicity,
        requested_toeplitz_strength=format_float(requested_strength),
        actual_toeplitz_strength=format_float(actual_strength),
        coefficient_commutator_norm=format_float(commutator),
        endpoint_map_rank=int(rank),
        endpoint_space_dimension=endpoint_dimension,
        endpoint_cokernel_dimension=cokernel_dimension,
        target_gap_norm=format_float(gap_norm),
        target_gap_over_strength_squared=format_float(
            gap_norm / strength**2
        ),
        correction_norm=format_float(correction_norm),
        correction_over_strength=format_float(
            correction_norm / strength
        ),
        target_trace_error=format_float(trace_error),
        range_equation_error=format_float(range_error),
        adjoint_identity_error=format_float(adjoint_error),
        full_second_stein_error=format_float(full_stein_error),
        lower_endpoint_error=format_float(lower_error),
        upper_endpoint_error=format_float(upper_error),
        all_checks_passed=verified,
    )


def write_records(
    records: list[EllipticCokernelRecord],
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
            "repeated_crabb_elliptic_cokernel_s70224.jsonl"
        ),
    )
    parser.add_argument(
        "--strengths",
        type=float,
        nargs="+",
        default=(0.03, 0.1, 0.3, 1.0, 3.0, 8.0),
    )
    return parser.parse_args()


def main() -> None:
    """Run the standard rank-jump and divisibility audit."""

    args = parse_args()
    records = [
        reduced_second_face(
            length,
            multiplicity,
            strength,
        )
        for multiplicity in (2, 3)
        for length in range(2, 6)
        for strength in args.strengths
    ]
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
