#!/usr/bin/env python3
"""Certify the full complex cubic response in length seven by polarization.

A direct characteristic expansion in six generic complex variables is
unnecessarily large.  The response coefficient is instead a homogeneous
cubic in the twelve real and imaginary Toeplitz coordinates.  This
checker evaluates the exact characteristic/Riemann recurrence on an
explicit 364-point unisolvent set:

* every coordinate vector;
* ``e_i + e_j`` and ``e_i - e_j`` for every pair;
* ``e_i + e_j + e_k`` for every triple.

The exact degree-three evaluation matrix has full rank 364.  Hence
agreement on these sparse points proves the two proposed complex
response formulas coefficientwise.  The same points have full rank for
homogeneous quadratics and also prove complete quadratic cancellation.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import itertools
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_quadratic_exact import (
    characteristic_series,
    crabb_normal_directions,
    directional_gradient_series,
    disk_model_series,
)
from crabb_full_disk_correction_isometry import plucker_correction
from crabb_full_disk_cubic_normal_exact import highest_mode_cubic_response


REAL_VARIABLE_COUNT = 12
LENGTH = 7
EXPECTED_ACTIVE_MODES = (3, 4)


@dataclass(frozen=True)
class PointAudit:
    """Exact response result on one sparse interpolation point."""

    point_index: int
    quadratic_response_zero: bool
    response_formulas_verified: bool
    observed_active_modes: tuple[int, ...]


@dataclass(frozen=True)
class ComplexL7CubicResponseRecord:
    """Metadata for the exact sparse-polarization response audit."""

    dimension: int
    length: int
    real_variable_count: int
    quadratic_monomial_count: int
    cubic_monomial_count: int
    unisolvent_point_count: int
    quadratic_evaluation_rank: int
    cubic_evaluation_rank: int
    normal_polarization_evaluation_count: int
    active_cubic_modes: tuple[int, ...]
    candidate_homogeneity_verified: bool
    quadratic_response_zero: bool
    response_formulas_verified: bool
    full_complex_response_verified: bool


def homogeneous_exponents(
    variable_count: int,
    degree: int,
) -> tuple[tuple[int, ...], ...]:
    """Return every exponent tuple of one homogeneous degree."""

    return tuple(
        exponents
        for exponents in itertools.product(
            range(degree + 1),
            repeat=variable_count,
        )
        if sum(exponents) == degree
    )


def unisolvent_points(
    variable_count: int,
) -> tuple[tuple[int, ...], ...]:
    """Return the sparse degree-three polarization set."""

    points = []
    for first in range(variable_count):
        point = [0] * variable_count
        point[first] = 1
        points.append(tuple(point))
    for first in range(variable_count):
        for second in range(first + 1, variable_count):
            for sign in (1, -1):
                point = [0] * variable_count
                point[first] = 1
                point[second] = sign
                points.append(tuple(point))
    for first in range(variable_count):
        for second in range(first + 1, variable_count):
            for third in range(second + 1, variable_count):
                point = [0] * variable_count
                point[first] = 1
                point[second] = 1
                point[third] = 1
                points.append(tuple(point))
    return tuple(points)


def exact_evaluation_rank(
    points: tuple[tuple[int, ...], ...],
    degree: int,
) -> tuple[int, int]:
    """Return monomial count and exact evaluation-matrix rank."""

    exponents = homogeneous_exponents(len(points[0]), degree)
    rows = [
        [
            sp.QQ(
                sp.prod(
                    coordinate**exponent
                    for coordinate, exponent in zip(
                        point,
                        monomial_exponents,
                        strict=True,
                    )
                )
            )
            for monomial_exponents in exponents
        ]
        for point in points
    ]
    matrix = sp.polys.matrices.DomainMatrix.from_list(rows, sp.QQ)
    return len(exponents), matrix.rank()


def response_candidates(
    direction: tuple[sp.Expr, ...],
) -> dict[int, sp.Expr]:
    """Return the proposed mode-three and mode-four cubic responses."""

    coefficients = direction[1:]
    coefficient_count = len(coefficients)

    def pluecker(left: int, right: int) -> sp.Expr:
        return coefficients[left] * sp.conjugate(
            coefficients[coefficient_count - 1 - right]
        ) - coefficients[right] * sp.conjugate(
            coefficients[coefficient_count - 1 - left]
        )

    response_three = sp.Rational(176, 147) * (
        3 * coefficients[5] * pluecker(0, 2)
        - 4 * coefficients[3] * pluecker(0, 4)
        - 2 * coefficients[3] * pluecker(1, 3)
        - 2 * sp.conjugate(coefficients[0]) * sp.conjugate(pluecker(0, 1))
    )
    return {
        3: sp.expand(response_three),
        4: highest_mode_cubic_response(direction),
    }


def candidate_homogeneity_verified() -> bool:
    """Check that both proposed responses are homogeneous real cubics."""

    real_variables = sp.symbols(
        f"x0:{REAL_VARIABLE_COUNT}",
        real=True,
    )
    coefficients = tuple(
        real_variables[index] + sp.I * real_variables[6 + index] for index in range(6)
    )
    candidates = response_candidates((sp.Integer(0), *coefficients))
    return all(
        sp.Poly(
            sp.expand_complex(component),
            *real_variables,
        ).is_homogeneous
        and sp.Poly(
            sp.expand_complex(component),
            *real_variables,
        ).total_degree()
        == 3
        for response in candidates.values()
        for component in (sp.re(response), sp.im(response))
    )


def audit_point(
    indexed_point: tuple[int, tuple[int, ...]],
) -> PointAudit:
    """Evaluate every true-normal polarization on one exact point."""

    point_index, point = indexed_point
    coefficients = tuple(
        sp.Integer(point[index]) + sp.I * sp.Integer(point[6 + index])
        for index in range(6)
    )
    direction = (sp.Integer(0), *coefficients)
    correction = plucker_correction(direction)
    operator, metric = disk_model_series(
        direction,
        order=3,
        correction=correction,
    )
    characteristic = characteristic_series(operator)
    candidates = response_candidates(direction)
    quadratic_zero = True
    formulas_verified = True
    active_modes = []
    for mode in range(3, LENGTH + 2):
        responses = [
            directional_gradient_series(
                operator,
                metric,
                characteristic,
                normal,
            )
            for normal in crabb_normal_directions([metric[0]], mode)
        ]
        if len(responses) != 2:
            raise RuntimeError(f"mode {mode} has {len(responses)} polarizations")
        quadratic_zero &= all(sp.simplify(response[2]) == 0 for response in responses)
        cubic_response = sp.expand(responses[0][3] + sp.I * responses[1][3])
        if cubic_response != 0:
            active_modes.append(mode)
        formulas_verified &= sp.simplify(cubic_response - candidates.get(mode, 0)) == 0
    return PointAudit(
        point_index=point_index,
        quadratic_response_zero=quadratic_zero,
        response_formulas_verified=formulas_verified,
        observed_active_modes=tuple(active_modes),
    )


def build_record(workers: int) -> ComplexL7CubicResponseRecord:
    """Run the exact unisolvent response audit."""

    points = unisolvent_points(REAL_VARIABLE_COUNT)
    quadratic_count, quadratic_rank = exact_evaluation_rank(
        points,
        degree=2,
    )
    cubic_count, cubic_rank = exact_evaluation_rank(
        points,
        degree=3,
    )
    indexed_points = tuple(enumerate(points))
    records = []
    if workers == 1:
        for completed, indexed_point in enumerate(indexed_points, 1):
            records.append(audit_point(indexed_point))
            if completed % 25 == 0 or completed == len(points):
                print(
                    f"verified {completed}/{len(points)} points",
                    flush=True,
                )
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(audit_point, indexed_point)
                for indexed_point in indexed_points
            ]
            for completed, future in enumerate(as_completed(futures), 1):
                records.append(future.result())
                if completed % 25 == 0 or completed == len(points):
                    print(
                        f"verified {completed}/{len(points)} points",
                        flush=True,
                    )
    records.sort(key=lambda record: record.point_index)
    active_modes = tuple(
        sorted({mode for record in records for mode in record.observed_active_modes})
    )
    homogeneity = candidate_homogeneity_verified()
    quadratic_zero = all(record.quadratic_response_zero for record in records)
    formulas_verified = all(record.response_formulas_verified for record in records)
    verified = bool(
        quadratic_rank == quadratic_count
        and cubic_rank == cubic_count == len(points)
        and active_modes == EXPECTED_ACTIVE_MODES
        and homogeneity
        and quadratic_zero
        and formulas_verified
    )
    if not verified:
        raise RuntimeError("the complex length-seven response audit failed")
    return ComplexL7CubicResponseRecord(
        dimension=8,
        length=LENGTH,
        real_variable_count=REAL_VARIABLE_COUNT,
        quadratic_monomial_count=quadratic_count,
        cubic_monomial_count=cubic_count,
        unisolvent_point_count=len(points),
        quadratic_evaluation_rank=quadratic_rank,
        cubic_evaluation_rank=cubic_rank,
        normal_polarization_evaluation_count=(len(points) * 2 * (LENGTH - 1)),
        active_cubic_modes=active_modes,
        candidate_homogeneity_verified=homogeneity,
        quadratic_response_zero=quadratic_zero,
        response_formulas_verified=formulas_verified,
        full_complex_response_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and serialize the exact sparse-polarization audit."""

    args = parse_args()
    if args.workers < 1:
        raise ValueError("workers must be positive")
    record = build_record(args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(asdict(record), sort_keys=True)
    args.output.write_text(f"{line}\n", encoding="utf-8")
    print(line, flush=True)


if __name__ == "__main__":
    main()
