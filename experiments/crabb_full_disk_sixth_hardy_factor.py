#!/usr/bin/env python3
"""Regenerate the Hardy factor behind the all-size sixth-order face.

For L176's recentered exact disk path, the coefficient of order three
in the canonical Berger/Hardy residual is a skew-symmetric matrix
``R_3``.  This checker verifies:

* the residual vanishes through order two;
* ``R_3`` is skew-symmetric;
* its two circle-grade anti-diagonals reproduce every L182 response;
* the endpoint base deficit is exactly ``8 ||R_3||^2``; and
* the stronger flux-only completed-square inequality is nonnegative.

The transfer identities are checked with unrestricted complex symbols
in small sizes and on exact Gaussian-rational directions thereafter.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_quadratic_exact import disk_model_series
from crabb_disk_toeplitz_quartic import extend, shift
from crabb_full_disk_base_jet_exact import (
    endpoint_delta_series,
    toeplitz_direction,
)
from crabb_full_disk_correction_isometry import plucker_correction
from crabb_full_disk_cubic_normal_exact import varied_disk_direction
from crabb_full_disk_cubic_response_formula import cubic_normal_response
from crabb_full_disk_terminal_sixth_exact import terminal_direction
from formal_riemann_series import laurent_modes


@dataclass(frozen=True)
class SymbolicHardyRecord:
    """One unrestricted-symbol transfer audit."""

    record_kind: str
    dimension: int
    length: int
    variable_count: int
    lower_residual_orders_zero: bool
    cubic_laurent_support_verified: bool
    cubic_residual_skew_symmetric: bool
    active_mode_count: int
    response_projection_verified: bool
    all_identities_verified: bool


@dataclass(frozen=True)
class ExactHardyRecord:
    """One exact base-factor and strong-face audit."""

    record_kind: str
    dimension: int
    length: int
    direction_kind: str
    active_mode_count: int
    base_deficit: str
    hardy_residual_norm_square: str
    strong_flux_gain: str
    strong_flux_margin: str
    lower_residual_orders_zero: bool
    cubic_laurent_support_verified: bool
    cubic_residual_skew_symmetric: bool
    response_projection_verified: bool
    base_hardy_factor_verified: bool
    strong_flux_inequality_verified: bool
    terminal_equality_verified: bool | None
    all_identities_verified: bool


def clean(matrix: sp.Matrix) -> sp.Matrix:
    """Expand every matrix entry."""

    return matrix.applyfunc(sp.expand)


def matrix_is_zero(matrix: sp.Matrix) -> bool:
    """Test exact entrywise zero after expansion."""

    return all(sp.expand(entry) == 0 for entry in matrix)


def hardy_residual_matrices(
    direction: tuple[sp.Expr, ...],
    order: int = 3,
) -> tuple[tuple[sp.Matrix, ...], bool]:
    """Return disk-chart Hardy residual coefficients through ``order``.

    Rows of each returned matrix are negative Fourier powers
    ``w^-1,...,w^-(L-1)``.  Columns are coefficient coordinates
    ``1,...,L-1``.  The orbit coordinate zero is removed by its exact
    moving ``H``-orthogonal projection.
    """

    if order < 2:
        raise ValueError("the recentered Hardy residual needs order >= 2")
    length = len(direction)
    dimension = length + 1
    coefficient_count = length - 1
    boundary = sp.symbols("hardy_boundary")
    correction = plucker_correction(direction)
    operator, _ = disk_model_series(
        direction,
        order=order,
        correction=correction,
    )

    hermitian = [
        extend(sp.eye(length) / 2),
        extend(toeplitz_direction(direction)),
        extend(correction),
        *[sp.zeros(dimension) for _ in range(order - 2)],
    ]
    diagonal = [coefficient[0, 0] for coefficient in hermitian]
    diagonal_inverse = [1 / diagonal[0]]
    for degree in range(1, order + 1):
        convolution = sum(
            (
                diagonal[source_degree]
                * diagonal_inverse[degree - source_degree]
                for source_degree in range(1, degree + 1)
            ),
            sp.Integer(0),
        )
        diagonal_inverse.append(
            sp.expand(-convolution / diagonal[0])
        )

    projected_rows = [
        sum(
            (
                hermitian[source_degree][0, :]
                * diagonal_inverse[degree - source_degree]
                for source_degree in range(degree + 1)
            ),
            sp.zeros(1, dimension),
        )
        for degree in range(order + 1)
    ]
    first = sp.eye(dimension)[:, 0]
    projection = [
        sp.eye(dimension) - first * projected_rows[0],
        *[
            -first * projected_rows[degree]
            for degree in range(1, order + 1)
        ],
    ]

    base_resolvent = sum(
        (
            operator[0] ** power / boundary**power
            for power in range(dimension)
        ),
        sp.zeros(dimension),
    )
    inverse = [base_resolvent]
    for degree in range(1, order + 1):
        forcing = sum(
            (
                operator[source_degree]
                * inverse[degree - source_degree]
                / boundary
                for source_degree in range(1, degree + 1)
            ),
            sp.zeros(dimension),
        )
        inverse.append(clean(base_resolvent * forcing))

    nilpotent_shift = shift(dimension)
    pulled = [
        clean(
            (sp.eye(dimension) - nilpotent_shift / boundary)
            * coefficient
        )
        for coefficient in inverse
    ]
    last = sp.eye(dimension)[:, -1]
    residual_matrices = []
    support_verified = True
    for degree in range(1, order + 1):
        residual = clean(
            sum(
                (
                    projection[source_degree]
                    * pulled[degree - source_degree]
                    * last
                    for source_degree in range(degree + 1)
                ),
                sp.zeros(dimension, 1),
            )
        )
        support_verified &= sp.expand(residual[0]) == 0
        matrix = sp.zeros(coefficient_count)
        for coordinate in range(1, length):
            modes = laurent_modes(
                sp.expand(residual[coordinate]),
                boundary,
            )
            support_verified &= all(
                -coefficient_count <= mode <= -1 or coefficient == 0
                for mode, coefficient in modes.items()
            )
            for power in range(1, length):
                matrix[power - 1, coordinate - 1] = modes.get(
                    -power,
                    0,
                )
        residual_matrices.append(matrix)
    return tuple(residual_matrices), support_verified


def response_projection_verified(
    direction: tuple[sp.Expr, ...],
    residual: sp.Matrix,
) -> bool:
    """Check both anti-diagonal formulas for every active response."""

    length = len(direction)
    for mode in range(3, length - 2):
        flux_coordinate = sp.expand(
            cubic_normal_response(direction, mode)
            * length**2
            / (16 * (4 * mode - 1))
        )
        lower = sum(
            (
                (right - left)
                * residual[left - 1, right - 1]
                for left in range(1, length)
                for right in range(left + 1, length)
                if left + right == length - mode
            ),
            sp.Integer(0),
        )
        upper = sum(
            (
                (right - left)
                * residual[left - 1, right - 1]
                for left in range(1, length)
                for right in range(left + 1, length)
                if left + right == length + mode
            ),
            sp.Integer(0),
        )
        if (
            sp.expand(lower - 4 * flux_coordinate) != 0
            or sp.expand(
                upper - 4 * sp.conjugate(flux_coordinate)
            )
            != 0
        ):
            return False
    return True


def residual_norm_square(residual: sp.Matrix) -> sp.Expr:
    """Return the exterior-square norm of a skew residual."""

    return sp.factor(
        sum(
            (
                residual[left, right]
                * sp.conjugate(residual[left, right])
                for left in range(residual.rows)
                for right in range(left + 1, residual.cols)
            ),
            sp.Integer(0),
        )
    )


def strong_flux_gain(direction: tuple[sp.Expr, ...]) -> sp.Expr:
    """Return the completed gain using only L173's flux curvature."""

    length = len(direction)
    gain = sp.Integer(0)
    for mode in range(3, length - 2):
        response_coordinate = (
            cubic_normal_response(direction, mode)
            * length**2
            / (16 * (4 * mode - 1))
        )
        gain += (
            256
            * response_coordinate
            * sp.conjugate(response_coordinate)
            / sp.binomial(length - mode, 3)
        )
    return sp.factor(sp.simplify(gain))


def symbolic_audit(length: int) -> SymbolicHardyRecord:
    """Audit the transfer identities with unrestricted symbols."""

    variables = sp.symbols(f"z0:{length - 1}")
    direction = (sp.Integer(0), *variables)
    (linear, quadratic, cubic), support = hardy_residual_matrices(
        direction
    )
    lower_zero = matrix_is_zero(linear) and matrix_is_zero(quadratic)
    skew = matrix_is_zero(cubic + cubic.T)
    projection = response_projection_verified(direction, cubic)
    verified = bool(lower_zero and support and skew and projection)
    if not verified:
        raise RuntimeError(
            f"symbolic Hardy transfer audit failed in length {length}"
        )
    return SymbolicHardyRecord(
        record_kind="symbolic_transfer",
        dimension=length + 1,
        length=length,
        variable_count=2 * (length - 1),
        lower_residual_orders_zero=lower_zero,
        cubic_laurent_support_verified=support,
        cubic_residual_skew_symmetric=skew,
        active_mode_count=max(0, length - 5),
        response_projection_verified=projection,
        all_identities_verified=verified,
    )


def exact_audit(
    direction: tuple[sp.Expr, ...],
    direction_kind: str,
    terminal: bool,
) -> ExactHardyRecord:
    """Audit the base norm and flux-only inequality on one direction."""

    length = len(direction)
    (linear, quadratic, cubic), support = hardy_residual_matrices(
        direction
    )
    lower_zero = matrix_is_zero(linear) and matrix_is_zero(quadratic)
    skew = matrix_is_zero(cubic + cubic.T)
    projection = response_projection_verified(direction, cubic)
    norm_square = residual_norm_square(cubic)
    base_deficit = sp.factor(
        -2 * endpoint_delta_series(direction, order=6)[6]
    )
    base_factor = sp.simplify(base_deficit - 8 * norm_square) == 0
    gain = strong_flux_gain(direction)
    margin = sp.factor(sp.simplify(base_deficit - gain))
    strong = bool(margin >= 0)
    terminal_equality = bool(margin == 0) if terminal else None
    verified = bool(
        lower_zero
        and support
        and skew
        and projection
        and base_factor
        and strong
        and terminal_equality is not False
    )
    if not verified:
        raise RuntimeError(
            f"exact Hardy factor audit failed in length {length} "
            f"on {direction_kind}"
        )
    return ExactHardyRecord(
        record_kind="exact_direction",
        dimension=length + 1,
        length=length,
        direction_kind=direction_kind,
        active_mode_count=max(0, length - 5),
        base_deficit=str(base_deficit),
        hardy_residual_norm_square=str(norm_square),
        strong_flux_gain=str(gain),
        strong_flux_margin=str(margin),
        lower_residual_orders_zero=lower_zero,
        cubic_laurent_support_verified=support,
        cubic_residual_skew_symmetric=skew,
        response_projection_verified=projection,
        base_hardy_factor_verified=base_factor,
        strong_flux_inequality_verified=strong,
        terminal_equality_verified=terminal_equality,
        all_identities_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbolic-maximum-length", type=int, default=8)
    parser.add_argument("--exact-maximum-length", type=int, default=10)
    parser.add_argument("--terminal-maximum-length", type=int, default=12)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and serialize the symbolic and exact audits."""

    args = parse_args()
    if (
        args.symbolic_maximum_length < 4
        or args.exact_maximum_length < 4
        or args.terminal_maximum_length < 6
    ):
        raise ValueError("maximum lengths are below their valid minima")

    records: list[SymbolicHardyRecord | ExactHardyRecord] = []
    for length in range(4, args.symbolic_maximum_length + 1):
        records.append(symbolic_audit(length))
        print(f"verified symbolic length {length}", flush=True)
    for length in range(4, args.exact_maximum_length + 1):
        records.append(
            exact_audit(
                varied_disk_direction(length, direction_index=0),
                direction_kind="dense_gaussian_rational",
                terminal=False,
            )
        )
        print(f"verified dense exact length {length}", flush=True)
    for length in range(6, args.terminal_maximum_length + 1):
        records.append(
            exact_audit(
                terminal_direction(length),
                direction_kind="terminal_gaussian_rational",
                terminal=True,
            )
        )
        print(f"verified terminal exact length {length}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(asdict(record), sort_keys=True) for record in records]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
