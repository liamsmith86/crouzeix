#!/usr/bin/env python3
"""Audit the explicit non-Crabb Gau--Wu scalar equality family."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class GauWuSharpRecord:
    """One exact nonnilpotent sharp-disk audit."""

    parameter: str
    support_factorization: str
    functional_value: str
    functional_norm_squared: str
    spectrum: tuple[str, ...]
    noncrabb_certificate: bool
    all_checks_passed: bool


def symbolic_identities() -> None:
    """Verify the two load-bearing identities with a free parameter."""

    parameter = sp.symbols("a", real=True)
    phase = sp.symbols("z", nonzero=True)
    eigenvalue = sp.symbols("lambda")
    edge = sp.sqrt(2 * (1 - parameter**2))
    matrix = sp.Matrix(
        [
            [0, edge, -2 * parameter],
            [0, parameter, edge],
            [0, 0, 0],
        ]
    )
    support = (matrix / phase + phase * matrix.T) / 2
    expected_characteristic = (
        (eigenvalue - 1)
        * (eigenvalue + 1)
        * (
            eigenvalue
            - parameter * (phase + phase**-1) / 2
        )
    )
    actual_characteristic = support.charpoly(eigenvalue).as_expr()
    if sp.simplify(actual_characteristic - expected_characteristic) != 0:
        raise RuntimeError("symbolic support factorization failed")

    identity = sp.eye(3)
    functional_value = sp.simplify(
        matrix
        * (matrix - parameter * identity)
        * (identity - parameter * matrix).inv()
    )
    if functional_value != 2 * sp.eye(3)[:, 0] * sp.eye(3)[2, :]:
        raise RuntimeError("symbolic Blaschke functional calculus failed")


def audit_parameter(parameter: Fraction) -> GauWuSharpRecord:
    """Verify one rational member exactly."""

    value = sp.Rational(parameter.numerator, parameter.denominator)
    edge = sp.sqrt(2 * (1 - value**2))
    matrix = sp.Matrix(
        [
            [0, edge, -2 * value],
            [0, value, edge],
            [0, 0, 0],
        ]
    )
    identity = sp.eye(3)
    functional_value = sp.simplify(
        matrix
        * (matrix - value * identity)
        * (identity - value * matrix).inv()
    )
    expected_value = sp.zeros(3)
    expected_value[0, 2] = 2
    gram = sp.simplify(functional_value.T * functional_value)
    norm_squared = max(gram.eigenvals(), key=lambda item: float(item))
    spectrum = tuple(
        sorted(
            (str(item) for item in matrix.eigenvals()),
            key=lambda item: (item != "0", item),
        )
    )
    expected_support = (
        "(lambda - 1)*(lambda + 1)"
        "*(lambda - a*(z + 1/z)/2)"
    )
    checks = (
        functional_value == expected_value
        and norm_squared == 4
        and matrix.eigenvals() == {sp.Integer(0): 2, value: 1}
        and value != 0
    )
    if not checks:
        raise RuntimeError(f"exact Gau--Wu audit failed at a={parameter}")

    return GauWuSharpRecord(
        parameter=str(parameter),
        support_factorization=expected_support,
        functional_value=str(functional_value.tolist()),
        functional_norm_squared=str(norm_squared),
        spectrum=spectrum,
        noncrabb_certificate=True,
        all_checks_passed=True,
    )


def standard_records() -> list[GauWuSharpRecord]:
    """Return the deterministic exact parameter grid."""

    symbolic_identities()
    return [
        audit_parameter(parameter)
        for parameter in (
            Fraction(-4, 5),
            Fraction(-1, 2),
            Fraction(1, 5),
            Fraction(1, 2),
            Fraction(3, 5),
            Fraction(4, 5),
        )
    ]


def write_records(
    records: list[GauWuSharpRecord],
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
            "experiments/gau_wu_noncrabb_sharp_stratum_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and report the exact audit."""

    arguments = parse_args()
    records = standard_records()
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
