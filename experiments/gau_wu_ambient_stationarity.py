#!/usr/bin/env python3
"""Audit ambient first-order stationarity at non-Crabb Gau--Wu models."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp

from gau_wu_disk_model import gau_wu_matrix, parameter_from_schwarz


@dataclass(frozen=True)
class GauWuStationarityRecord:
    """One exact rational direction audit."""

    schwarz_parameter: str
    model_parameter: str
    direction_index: int
    raw_functional_real: str
    support_functional_real: str
    stationarity_residual: str
    all_checks_passed: bool


def frechet_functional(
    matrix: sp.Matrix,
    direction: sp.Matrix,
    parameter: sp.Expr,
) -> sp.Expr:
    """Return e1* Df_a(A)[E] e3 exactly."""

    identity = sp.eye(3)
    resolvent = (identity - parameter * matrix).inv()
    numerator = matrix * (matrix - parameter * identity)
    derivative = (
        (
            direction * (matrix - parameter * identity)
            + matrix * direction
        )
        * resolvent
        + numerator
        * resolvent
        * (parameter * direction)
        * resolvent
    )
    return sp.simplify(derivative[0, 2])


def support_mean_matrix() -> tuple[sp.Symbol, sp.Expr, sp.Matrix]:
    """Derive the free-symbol residue matrix in L331 (13)."""

    schwarz_parameter = sp.symbols("q", real=True, nonzero=True)
    phase = sp.symbols("z", nonzero=True)
    parameter = parameter_from_schwarz(schwarz_parameter)
    edge = (
        sp.sqrt(2)
        * (1 - schwarz_parameter**2)
        / (1 + schwarz_parameter**2)
    )
    matrix = sp.Matrix(
        [
            [0, edge, -2 * parameter],
            [0, parameter, edge],
            [0, 0, 0],
        ]
    )
    support = (matrix / phase + phase * matrix.T) / 2
    middle_eigenvalue = parameter * (phase + phase**-1) / 2
    projection = (
        (support + sp.eye(3))
        * (support - middle_eigenvalue * sp.eye(3))
        / (2 * (1 - middle_eigenvalue))
    )
    poisson = (
        (1 - parameter**2)
        / ((1 - parameter * phase) * (1 - parameter / phase))
    )
    mean_matrix = sp.zeros(3)
    for row in range(3):
        for column in range(3):
            integrand = sp.factor(
                2
                * (1 + poisson)
                * projection[column, row]
                / phase**2
            )
            residue = sp.residue(integrand, phase, 0)
            residue += sp.limit(
                (phase - schwarz_parameter) * integrand,
                phase,
                schwarz_parameter,
            )
            residue += sp.limit(
                (phase - parameter) * integrand,
                phase,
                parameter,
            )
            mean_matrix[row, column] = sp.factor(residue)

    return schwarz_parameter, parameter, mean_matrix


def deterministic_directions() -> list[sp.Matrix]:
    """Return a small spanning-style collection of complex directions."""

    directions: list[sp.Matrix] = []
    for row, column in (
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
        (1, 2),
        (2, 0),
        (2, 1),
        (2, 2),
    ):
        real_direction = sp.zeros(3)
        real_direction[row, column] = 1
        directions.append(real_direction)
        imaginary_direction = sp.zeros(3)
        imaginary_direction[row, column] = sp.I
        directions.append(imaginary_direction)
    directions.append(
        sp.Matrix(
            [
                [1 + sp.I, -2 + sp.I, sp.Rational(1, 3)],
                [2 * sp.I, -1, 3 - sp.I],
                [sp.Rational(2, 5), 1 + 2 * sp.I, -sp.I],
            ]
        )
    )
    return directions


def audit_records() -> list[GauWuStationarityRecord]:
    """Run the symbolic residue proof and rational direction checks."""

    free_q, free_a, mean_matrix = support_mean_matrix()
    free_edge = sp.sqrt(2) * (1 - free_q**2) / (1 + free_q**2)
    expected_mean = sp.Matrix(
        [
            [0, free_edge, -free_a],
            [0, 2 * free_a, free_edge],
            [0, 0, 0],
        ]
    )
    if any(
        sp.simplify(entry) != 0
        for entry in (mean_matrix - expected_mean)
    ):
        raise RuntimeError("free-symbol support residue identity failed")

    records: list[GauWuStationarityRecord] = []
    directions = deterministic_directions()
    for q_value in (
        Fraction(-1, 2),
        Fraction(-1, 5),
        Fraction(1, 5),
        Fraction(1, 2),
        Fraction(2, 3),
    ):
        q_exact = sp.Rational(q_value.numerator, q_value.denominator)
        parameter = sp.simplify(parameter_from_schwarz(q_exact))
        matrix = gau_wu_matrix(parameter)
        specialized_mean = sp.simplify(mean_matrix.subs(free_q, q_exact))
        for index, direction in enumerate(directions):
            raw = frechet_functional(matrix, direction, parameter)
            support_value = sp.trace(
                specialized_mean.T * direction
            )
            residual = sp.simplify(
                sp.re(raw).expand(complex=True)
                - sp.re(support_value).expand(complex=True)
            )
            if residual != 0:
                raise RuntimeError(
                    "stationarity identity failed for "
                    f"q={q_value}, direction={index}"
                )
            records.append(
                GauWuStationarityRecord(
                    schwarz_parameter=str(q_value),
                    model_parameter=str(parameter),
                    direction_index=index,
                    raw_functional_real=str(
                        sp.simplify(sp.re(raw).expand(complex=True))
                    ),
                    support_functional_real=str(
                        sp.simplify(
                            sp.re(support_value).expand(complex=True)
                        )
                    ),
                    stationarity_residual=str(residual),
                    all_checks_passed=True,
                )
            )
    return records


def write_records(
    records: list[GauWuStationarityRecord],
    output: Path,
) -> str:
    """Write JSON Lines atomically and return its SHA-256."""

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
            "experiments/gau_wu_ambient_stationarity_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and report the exact audit."""

    arguments = parse_args()
    records = audit_records()
    digest = write_records(records, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "record_count": len(records),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
