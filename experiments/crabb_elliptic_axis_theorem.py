#!/usr/bin/env python3
"""Independently audit the all-size elliptic Crabb-axis certificate.

Unlike ``crabb_elliptic_axis.py``, this checker uses explicit DCT-I and
high-precision Jacobi formulas rather than an eigensolver or an SDP.  It
tests every identity used in L117, including the periodized-sech/Poisson
formula and the rank-one Stein defect.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp


DEFAULT_ELLIPSE_PARAMETERS = ("0.01", "0.05", "0.15", "0.4", "0.6", "0.8")


@dataclass(frozen=True)
class AuditRecord:
    dimension: int
    ellipse_parameter: str
    orthogonality_error: str
    jacobi_kernel_error: str
    dct_offdiagonal_error: str
    alias_error: str
    poisson_error: str
    reflection_error: str
    stein_rank_one_error: str
    endpoint_error: str
    minimum_metric_increment: str
    strict_bound_gap: str


def scientific(value: mp.mpf, digits: int = 12) -> str:
    """Serialize a high-precision number without demoting it to binary64."""

    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def modulus_from_nome(nome: mp.mpf) -> mp.mpf:
    """Return Jacobi's modulus from its nome."""

    return (mp.jtheta(2, 0, nome) / mp.jtheta(3, 0, nome)) ** 2


def dct_one(dimension: int) -> mp.matrix:
    """Return the orthogonal DCT-I matrix with endpoint normalization."""

    length = dimension - 1
    endpoints = [
        1 / mp.sqrt(2) if index in (0, length) else mp.mpf(1)
        for index in range(dimension)
    ]
    transform = mp.matrix(dimension)
    for row in range(dimension):
        for column in range(dimension):
            transform[row, column] = (
                mp.sqrt(2 / mp.mpf(length))
                * endpoints[row]
                * endpoints[column]
                * mp.cos(mp.pi * row * column / length)
            )
    return transform


def maximum_entry(matrix: mp.matrix) -> mp.mpf:
    """Return the largest entry modulus of a matrix."""

    return max(
        abs(matrix[row, column])
        for row in range(matrix.rows)
        for column in range(matrix.cols)
    )


def periodized_sech(
    index: int, length: int, logarithmic_parameter: mp.mpf, tolerance: mp.mpf
) -> mp.mpf:
    """Sum ``sech((index+2*length*n)*ell)`` symmetrically."""

    total = 1 / mp.cosh(index * logarithmic_parameter)
    radius = 1
    while True:
        increment = 1 / mp.cosh(
            (index + 2 * length * radius) * logarithmic_parameter
        )
        increment += 1 / mp.cosh(
            (index - 2 * length * radius) * logarithmic_parameter
        )
        total += increment
        if increment <= tolerance * total:
            return total
        radius += 1
        if radius > 1_000_000:
            raise RuntimeError("periodized sech sum did not converge")


def audit_case(
    dimension: int, ellipse_parameter: str, precision: int
) -> AuditRecord:
    """Audit one dimension/ellipse pair at the active precision."""

    length = dimension - 1
    c = mp.mpf(ellipse_parameter)
    ell = -mp.log(c)
    modulus = modulus_from_nome(c**2)
    parameter = modulus**2
    complement = mp.sqrt(1 - parameter)
    quarter_period = mp.ellipk(parameter)
    endpoint_factors = [
        1 / mp.sqrt(2) if index in (0, length) else mp.mpf(1)
        for index in range(dimension)
    ]
    transform = dct_one(dimension)
    identity = mp.eye(dimension)
    orthogonality_error = maximum_entry(transform * transform.T - identity)

    nodes = mp.matrix(dimension, dimension)
    beta = mp.matrix(dimension, 1)
    arguments: list[mp.mpf] = []
    for index in range(dimension):
        argument = 2 * quarter_period * index / length
        arguments.append(argument)
        cn = mp.ellipfun("cn", argument, parameter)
        dn = mp.ellipfun("dn", argument, parameter)
        nodes[index, index] = mp.sqrt(modulus) * cn / dn
        beta[index] = endpoint_factors[index] * complement / dn

    gramian = mp.matrix(dimension)
    addition_gramian = mp.matrix(dimension)

    def jacobi_sum(argument: mp.mpf) -> mp.mpf:
        return mp.ellipfun("dn", argument, parameter) + modulus * mp.ellipfun(
            "cn", argument, parameter
        )

    for row in range(dimension):
        for column in range(dimension):
            gramian[row, column] = beta[row] * beta[column] / (
                1 - nodes[row, row] * nodes[column, column]
            )
            addition_gramian[row, column] = (
                endpoint_factors[row]
                * endpoint_factors[column]
                * (
                    jacobi_sum(arguments[row] - arguments[column])
                    + jacobi_sum(arguments[row] + arguments[column])
                )
                / 2
            )
    jacobi_kernel_error = maximum_entry(gramian - addition_gramian)

    kernel_metric = transform * gramian * transform.T
    kernel_scale = kernel_metric[0, 0]
    offdiagonal_error = max(
        abs(kernel_metric[row, column] / kernel_scale)
        for row in range(dimension)
        for column in range(dimension)
        if row != column
    )

    summation_tolerance = mp.power(10, -precision + 30)
    periodized = [
        periodized_sech(index, length, ell, summation_tolerance)
        for index in range(dimension)
    ]
    weights = [value / periodized[0] for value in periodized]
    alias_error = max(
        abs(kernel_metric[index, index] / kernel_scale - weights[index])
        for index in range(dimension)
    )

    descended_modulus = modulus_from_nome(c ** (2 * length))
    complementary_parameter = 1 - descended_modulus**2
    complementary_period = mp.ellipk(complementary_parameter)
    poisson_weights = [
        mp.ellipfun(
            "dn",
            complementary_period * index / length,
            complementary_parameter,
        )
        for index in range(dimension)
    ]
    poisson_error = max(
        abs(weights[index] - poisson_weights[index])
        for index in range(dimension)
    )
    reflection_error = max(
        abs(weights[index] * weights[length - index] - descended_modulus)
        for index in range(dimension)
    )

    symmetric_pullback = transform * nodes * transform.T
    direct_metric = mp.diag(weights)
    defect = (
        direct_metric
        - symmetric_pullback * direct_metric * symmetric_pullback
    )
    physical_beta = transform * beta / mp.sqrt(kernel_scale)
    stein_error = maximum_entry(defect - physical_beta * physical_beta.T)

    metric_diagonal = [
        weights[index] / c**index for index in range(dimension)
    ]
    predicted_bound = descended_modulus / c**length
    endpoint_error = abs(metric_diagonal[-1] - predicted_bound)
    minimum_increment = min(
        metric_diagonal[index + 1] - metric_diagonal[index]
        for index in range(length)
    )
    strict_bound_gap = 4 - predicted_bound

    tolerance = mp.power(10, -precision // 2)
    errors = (
        orthogonality_error,
        jacobi_kernel_error,
        offdiagonal_error,
        alias_error,
        poisson_error,
        reflection_error,
        stein_error,
        endpoint_error,
    )
    if max(errors) > tolerance:
        raise RuntimeError(
            f"identity audit failed for p={dimension}, c={ellipse_parameter}"
        )
    if minimum_increment <= 0 or strict_bound_gap <= 0:
        raise RuntimeError(
            f"strict inequality audit failed for p={dimension}, c={ellipse_parameter}"
        )

    return AuditRecord(
        dimension=dimension,
        ellipse_parameter=ellipse_parameter,
        orthogonality_error=scientific(orthogonality_error),
        jacobi_kernel_error=scientific(jacobi_kernel_error),
        dct_offdiagonal_error=scientific(offdiagonal_error),
        alias_error=scientific(alias_error),
        poisson_error=scientific(poisson_error),
        reflection_error=scientific(reflection_error),
        stein_rank_one_error=scientific(stein_error),
        endpoint_error=scientific(endpoint_error),
        minimum_metric_increment=scientific(minimum_increment),
        strict_bound_gap=scientific(strict_bound_gap),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=30)
    parser.add_argument(
        "--ellipse-parameters",
        nargs="+",
        default=DEFAULT_ELLIPSE_PARAMETERS,
    )
    parser.add_argument("--precision", type=int, default=260)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    if args.precision < 80:
        raise ValueError("precision must be at least 80 decimal digits")
    if any(not 0 < mp.mpf(value) < 1 for value in args.ellipse_parameters):
        raise ValueError("all ellipse parameters must lie in (0,1)")

    mp.mp.dps = args.precision
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as output:
        for dimension in range(args.minimum_size, args.maximum_size + 1):
            for ellipse_parameter in args.ellipse_parameters:
                record = audit_case(
                    dimension, ellipse_parameter, args.precision
                )
                line = json.dumps(asdict(record), sort_keys=True)
                print(line, flush=True)
                output.write(line + "\n")
                output.flush()


if __name__ == "__main__":
    main()
