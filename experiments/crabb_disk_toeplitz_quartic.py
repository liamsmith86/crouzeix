#!/usr/bin/env python3
"""Regenerate the exact Toeplitz disk chart and its quartic metric jet.

The local disk chart is parametrized by a positive Hermitian ``L x L``
matrix ``H``.  After extending it by one zero row and column and writing
``R`` for the unweighted shift, set

    K = H + R* H R,
    X = 2 K**(-1/2) H R K**(-1/2).

For a Toeplitz perturbation ``H = I/2 + s Z``, this script constructs the
rank-one Stein metric in the similar coefficient coordinates and expands
its two simple endpoint generalized eigenvalues through order four.  It
checks the all-size formula on fully symbolic directions at small sizes
and on exact rational directions thereafter.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class ToeplitzQuarticRecord:
    dimension: int
    direction_kind: str
    disk_factorization: bool
    endpoint_column_recurrence: bool
    lower_endpoint_identity: bool
    upper_endpoint_identity: bool
    quartic_identity: bool
    quartic_value: str


def shift(dimension: int) -> sp.Matrix:
    """Return the unweighted nilpotent superdiagonal shift."""

    matrix = sp.zeros(dimension)
    for index in range(dimension - 1):
        matrix[index, index + 1] = 1
    return matrix


def extend(matrix: sp.Matrix) -> sp.Matrix:
    """Append one zero row and column."""

    result = sp.zeros(matrix.rows + 1)
    result[: matrix.rows, : matrix.cols] = matrix
    return result


def toeplitz_direction(
    length: int,
    symbolic: bool,
) -> tuple[sp.Matrix, list[sp.Expr]]:
    """Return a zero-diagonal Hermitian Toeplitz direction."""

    if symbolic:
        real = sp.symbols(f"ar1:{length}", real=True)
        imaginary = sp.symbols(f"ai1:{length}", real=True)
        coefficients = [
            real[index] + sp.I * imaginary[index]
            for index in range(length - 1)
        ]
    else:
        coefficients = [
            sp.Rational(index + 1, index + 3)
            + sp.I * sp.Rational(index + 2, index + 5)
            for index in range(length - 1)
        ]

    direction = sp.zeros(length)
    for offset, coefficient in enumerate(coefficients, start=1):
        for row in range(length - offset):
            direction[row, row + offset] = coefficient
            direction[row + offset, row] = sp.conjugate(coefficient)
    return direction, coefficients


def inverse_series(
    base: sp.Matrix,
    tangent: sp.Matrix,
    order: int,
) -> list[sp.Matrix]:
    """Expand ``(base+s*tangent)^-1`` through ``order``."""

    coefficients = [base.inv()]
    for _ in range(order):
        coefficients.append(
            (-coefficients[0] * tangent * coefficients[-1]).applyfunc(
                sp.expand
            )
        )
    return coefficients


def lyapunov_series(
    operator: list[sp.Matrix],
    defect: list[sp.Matrix],
    order: int,
) -> list[sp.Matrix]:
    """Solve the rank-one Stein equation coefficient by coefficient."""

    dimension = operator[0].rows
    base = operator[0]
    defect_adjoint = [
        vector.conjugate().T
        for vector in defect
    ]
    metric: list[sp.Matrix] = []
    for degree in range(order + 1):
        forcing = sp.zeros(dimension)
        for left_degree in range(degree + 1):
            forcing += (
                defect[left_degree]
                * defect_adjoint[degree - left_degree]
            )
        for left_degree in range(degree + 1):
            for metric_degree in range(degree + 1 - left_degree):
                right_degree = degree - left_degree - metric_degree
                if (
                    left_degree == 0
                    and right_degree == 0
                    and metric_degree == degree
                ):
                    continue
                if metric_degree < len(metric):
                    forcing += (
                        operator[left_degree].conjugate().T
                        * metric[metric_degree]
                        * operator[right_degree]
                    )

        solution = sp.zeros(dimension)
        power = sp.eye(dimension)
        for _ in range(dimension):
            solution += power.conjugate().T * forcing * power
            power = power * base
        metric.append(solution.applyfunc(sp.expand))
    return metric


def generalized_endpoint_series(
    metric: list[sp.Matrix],
    metric_base: sp.Matrix,
    metric_tangent: sp.Matrix,
    endpoint: int,
    initial_value: sp.Rational,
    order: int,
) -> list[sp.Expr]:
    """Expand one simple generalized eigenvalue of ``(metric,K)``."""

    dimension = metric_base.rows
    eigenvalues = [initial_value]
    vectors = [sp.eye(dimension)[:, endpoint]]
    pencil_base = metric[0] - initial_value * metric_base

    for degree in range(1, order + 1):
        known = (
            metric[degree]
            - eigenvalues[degree - 1] * metric_tangent
        ) * vectors[0]
        for pencil_degree in range(1, degree):
            pencil = (
                metric[pencil_degree]
                - eigenvalues[pencil_degree] * metric_base
                - eigenvalues[pencil_degree - 1] * metric_tangent
            )
            known += pencil * vectors[degree - pencil_degree]

        eigenvalue = sp.cancel(
            known[endpoint] / metric_base[endpoint, endpoint]
        )
        eigenvalues.append(eigenvalue)
        known -= eigenvalue * metric_base * vectors[0]

        vector = sp.zeros(dimension, 1)
        for index in range(dimension):
            if index != endpoint:
                vector[index] = sp.cancel(
                    -known[index] / pencil_base[index, index]
                )
        vectors.append(vector)
    return eigenvalues


def quartic_form(coefficients: list[sp.Expr]) -> sp.Expr:
    """Return the Cauchy--Schwarz quartic ``||z||^4-|z^T Jz|^2``."""

    norm_square = sum(
        coefficient * sp.conjugate(coefficient)
        for coefficient in coefficients
    )
    reversal_pairing = sum(
        coefficients[index] * coefficients[-1 - index]
        for index in range(len(coefficients))
    )
    return sp.expand(
        norm_square**2
        - reversal_pairing * sp.conjugate(reversal_pairing)
    )


def audit_dimension(
    dimension: int,
    symbolic: bool,
) -> ToeplitzQuarticRecord:
    """Audit the disk identity and the endpoint quartic in one size."""

    order = 4
    length = dimension - 1
    direction, coefficients = toeplitz_direction(length, symbolic)
    base_half = sp.eye(length) / 2
    base = extend(base_half)
    tangent = extend(direction)
    nilpotent_shift = shift(dimension)
    metric_base = (
        base + nilpotent_shift.T * base * nilpotent_shift
    )
    metric_tangent = (
        tangent
        + nilpotent_shift.T * tangent * nilpotent_shift
    )

    circle_parameter, chart_parameter = sp.symbols("w s")
    chart_matrix = base + chart_parameter * tangent
    chart_metric = (
        chart_matrix
        + nilpotent_shift.T * chart_matrix * nilpotent_shift
    )
    disk_left = (
        chart_metric
        - (
            1 / circle_parameter
            * chart_matrix
            * nilpotent_shift
            + circle_parameter
            * nilpotent_shift.T
            * chart_matrix
        )
    )
    disk_right = (
        (sp.eye(dimension) - circle_parameter * nilpotent_shift.T)
        * chart_matrix
        * (
            sp.eye(dimension)
            - nilpotent_shift / circle_parameter
        )
    )
    disk_factorization = disk_left.equals(disk_right)

    inverse = inverse_series(metric_base, metric_tangent, order)
    operator: list[sp.Matrix] = []
    for degree in range(order + 1):
        coefficient = inverse[degree] * base * nilpotent_shift
        if degree:
            coefficient += (
                inverse[degree - 1]
                * tangent
                * nilpotent_shift
            )
        operator.append((2 * coefficient).applyfunc(sp.expand))

    zero_vector = sp.zeros(dimension, 1)
    defect = [
        base[:, 0],
        tangent[:, 0],
        *[zero_vector for _ in range(order - 1)],
    ]
    stein_metric = lyapunov_series(operator, defect, order)
    quartic = quartic_form(coefficients)
    last_basis_vector = sp.eye(dimension)[:, -1]
    endpoint_column_recurrence = all(
        entry == 0
        for vector in (
            stein_metric[1] * last_basis_vector
            - 2 * metric_tangent * last_basis_vector,
            stein_metric[2] * last_basis_vector,
            stein_metric[3] * last_basis_vector,
        )
        for entry in vector.applyfunc(sp.expand)
    ) and sp.expand(
        (
            last_basis_vector.conjugate().T
            * stein_metric[4]
            * last_basis_vector
        )[0]
        + 8 * quartic
    ) == 0
    lower = generalized_endpoint_series(
        stein_metric,
        metric_base,
        metric_tangent,
        endpoint=0,
        initial_value=sp.Rational(1, 2),
        order=order,
    )
    upper = generalized_endpoint_series(
        stein_metric,
        metric_base,
        metric_tangent,
        endpoint=dimension - 1,
        initial_value=sp.Integer(2),
        order=order,
    )

    lower_identity = all(
        sp.expand(coefficient) == 0
        for coefficient in lower[1:]
    )
    upper_identity = (
        all(
            sp.expand(coefficient) == 0
            for coefficient in upper[1:4]
        )
        and sp.expand(upper[4] + 16 * quartic) == 0
    )
    quartic_identity = lower_identity and upper_identity

    if not (
        disk_factorization
        and endpoint_column_recurrence
        and lower_identity
        and upper_identity
        and quartic_identity
    ):
        raise RuntimeError(
            f"Toeplitz quartic audit failed in size {dimension}: "
            f"disk={disk_factorization}, "
            f"endpoint_column={endpoint_column_recurrence}, "
            f"lower={lower_identity}, "
            f"upper={upper_identity}, quartic={quartic_identity}, "
            f"lower_coefficients={lower}, upper_coefficients={upper}"
        )

    return ToeplitzQuarticRecord(
        dimension=dimension,
        direction_kind="symbolic" if symbolic else "exact_rational",
        disk_factorization=disk_factorization,
        endpoint_column_recurrence=endpoint_column_recurrence,
        lower_endpoint_identity=lower_identity,
        upper_endpoint_identity=upper_identity,
        quartic_identity=quartic_identity,
        quartic_value=(
            "||z||^4-|z^T Jz|^2"
            if symbolic
            else str(sp.factor(quartic))
        ),
    )


def audit_optimized_size_four() -> dict[str, bool | str]:
    """Verify that defect-vector optimization does not fill the quartic null."""

    dimension = 4
    length = dimension - 1
    order = 4
    direction, coefficients = toeplitz_direction(length, symbolic=True)
    quartic = quartic_form(coefficients)
    base = extend(sp.eye(length) / 2)
    tangent = extend(direction)
    nilpotent_shift = shift(dimension)
    metric_base = base + nilpotent_shift.T * base * nilpotent_shift
    metric_tangent = (
        tangent + nilpotent_shift.T * tangent * nilpotent_shift
    )
    inverse = inverse_series(metric_base, metric_tangent, order)
    operator: list[sp.Matrix] = []
    for degree in range(order + 1):
        coefficient = inverse[degree] * base * nilpotent_shift
        if degree:
            coefficient += (
                inverse[degree - 1]
                * tangent
                * nilpotent_shift
            )
        operator.append((2 * coefficient).applyfunc(sp.expand))

    real = sp.symbols(f"ur1:{dimension}", real=True)
    imaginary = sp.symbols(f"ui1:{dimension}", real=True)
    correction = sp.Matrix(
        [
            0,
            *[
                real[index] + sp.I * imaginary[index]
                for index in range(dimension - 1)
            ],
        ]
    )
    zero_vector = sp.zeros(dimension, 1)
    defect = [
        base[:, 0],
        tangent[:, 0],
        correction,
        zero_vector,
        zero_vector,
    ]
    stein_metric = lyapunov_series(operator, defect, order)
    lower = generalized_endpoint_series(
        stein_metric,
        metric_base,
        metric_tangent,
        endpoint=0,
        initial_value=sp.Rational(1, 2),
        order=order,
    )
    upper = generalized_endpoint_series(
        stein_metric,
        metric_base,
        metric_tangent,
        endpoint=dimension - 1,
        initial_value=sp.Integer(2),
        order=order,
    )
    ratio: list[sp.Expr] = []
    for degree in range(order + 1):
        known = upper[degree] - sum(
            ratio[index] * lower[degree - index]
            for index in range(degree)
        )
        ratio.append(sp.cancel(known / lower[0]))

    correction_square = (
        16
        * sum(
            real[index] ** 2 + imaginary[index] ** 2
            for index in range(2)
        )
        + sp.Rational(32, 3)
        * (real[2] ** 2 + imaginary[2] ** 2)
    )
    optimized_coefficient = -32 * quartic
    expected_raw_coefficient = optimized_coefficient + correction_square
    raw_coefficient = sp.expand(ratio[4])
    coefficient_identity = (
        sp.expand(raw_coefficient - expected_raw_coefficient) == 0
    )
    if not coefficient_identity:
        raise RuntimeError("the derived size-four defect quartic changed")
    gradient_solution = sp.solve(
        [
            sp.diff(raw_coefficient, variable)
            for variable in (*real, *imaginary)
        ],
        (*real, *imaginary),
        dict=True,
    )
    expected_solution = {
        variable: 0
        for variable in (*real, *imaginary)
    }
    if gradient_solution != [expected_solution]:
        raise RuntimeError("the size-four defect correction was not minimized")
    return {
        "dimension": 4,
        "optimized_defect_quartic": str(
            sp.factor(optimized_coefficient)
        ),
        "derived_coefficient_identity": coefficient_identity,
        "quartic_null_is_not_filled": True,
    }


def audit_palindromic_null(maximum_size: int = 8) -> dict[str, object]:
    """Check exact finite points on the quartic null stratum.

    This is supporting evidence for the next equality-stratum problem, not
    an all-size proof that the displayed metric is globally optimal there.
    """

    checked_dimensions: list[int] = []
    parameter = sp.Rational(1, 20)
    eigenvalue = sp.symbols("lambda")
    for dimension in range(3, maximum_size + 1):
        length = dimension - 1
        coefficient_count = length - 1
        half = [
            sp.Rational(index + 1, index + 2)
            for index in range((coefficient_count + 1) // 2)
        ]
        coefficients = [
            half[
                min(index, coefficient_count - 1 - index)
            ]
            for index in range(coefficient_count)
        ]
        direction = sp.zeros(length)
        for offset, coefficient in enumerate(coefficients, start=1):
            for row in range(length - offset):
                direction[row, row + offset] = coefficient
                direction[row + offset, row] = coefficient

        chart = extend(sp.eye(length) / 2 + parameter * direction)
        nilpotent_shift = shift(dimension)
        coordinate_metric = (
            chart + nilpotent_shift.T * chart * nilpotent_shift
        )
        operator = (
            2
            * coordinate_metric.inv()
            * chart
            * nilpotent_shift
        )
        defect = chart[:, 0]

        variables = sp.symbols(f"m0:{dimension * dimension}")
        stein_metric = sp.Matrix(dimension, dimension, variables)
        equation = (
            stein_metric
            - operator.T * stein_metric * operator
            - defect * defect.T
        )
        solution = sp.solve(
            list(equation),
            variables,
            dict=True,
        )
        if len(solution) != 1:
            raise RuntimeError("the exact null-stratum Stein solve failed")
        stein_metric = stein_metric.subs(solution[0])
        characteristic = sp.factor(
            (stein_metric - eigenvalue * coordinate_metric).det()
        )
        expected_roots = (
            (eigenvalue - sp.Rational(1, 2))
            * (eigenvalue - 1) ** (dimension - 2)
            * (eigenvalue - 2)
        )
        quotient, remainder = sp.div(
            sp.Poly(characteristic, eigenvalue),
            sp.Poly(expected_roots, eigenvalue),
        )
        if remainder.as_expr() != 0 or quotient.degree() != 0:
            raise RuntimeError(
                f"null-stratum spectrum failed in size {dimension}"
            )
        checked_dimensions.append(dimension)

    return {
        "palindromic_null_checked_dimensions": checked_dimensions,
        "generalized_spectrum": ["1/2", "1", "2"],
        "finite_exact_evidence_only": True,
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=10)
    parser.add_argument(
        "--symbolic-maximum-size",
        type=int,
        default=6,
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run all exact audits and optionally persist JSONL records."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    records = [
        audit_dimension(
            dimension,
            symbolic=dimension <= args.symbolic_maximum_size,
        )
        for dimension in range(args.minimum_size, args.maximum_size + 1)
    ]
    lines = [
        *[
            json.dumps(asdict(record), sort_keys=True)
            for record in records
        ],
        json.dumps(audit_optimized_size_four(), sort_keys=True),
        json.dumps(audit_palindromic_null(), sort_keys=True),
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
