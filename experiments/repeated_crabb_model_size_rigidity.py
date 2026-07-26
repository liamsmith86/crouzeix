#!/usr/bin/env python3
"""Audit the power invariant separating one Crabb chain from repeats."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class ModelSizeRigidityRecord:
    """One exact incompatible model-size comparison."""

    crabb_dimension: int
    repeated_multiplicity: int
    model_dimension: int
    separating_power: int
    repeated_power_norm_square: str
    long_model_power_norm_square: str
    squared_norm_gap: str
    all_checks_passed: bool


def crabb_shift(dimension: int) -> sp.Matrix:
    """Return the superdiagonal Crabb weighted shift."""

    matrix = sp.zeros(dimension)
    for index in range(dimension - 1):
        matrix[index, index + 1] = (
            sp.sqrt(2)
            if index in {0, dimension - 2}
            else sp.Integer(1)
        )
    return matrix


def squared_power_norm(matrix: sp.Matrix, exponent: int) -> sp.Expr:
    """Return the exact squared norm of a weighted-shift power."""

    power = matrix**exponent
    gram = sp.simplify(power.conjugate().T * power)
    if not gram.is_diagonal():
        raise RuntimeError("Crabb power Gram is not diagonal")
    return max(sp.simplify(gram[index, index]) for index in range(matrix.rows))


def audit_case(
    crabb_dimension: int,
    repeated_multiplicity: int,
) -> ModelSizeRigidityRecord:
    """Compare the limiting restrictions at the separating power."""

    exponent = crabb_dimension - 1
    repeated = sp.kronecker_product(
        crabb_shift(crabb_dimension),
        sp.eye(repeated_multiplicity),
    )
    model_dimension = crabb_dimension * repeated_multiplicity
    long_model = crabb_shift(model_dimension)

    repeated_norm = squared_power_norm(repeated, exponent)
    long_norm = squared_power_norm(long_model, exponent)
    gap = sp.simplify(repeated_norm - long_norm)
    verified = bool(repeated_norm == 4 and long_norm == 2 and gap == 2)
    if not verified:
        raise RuntimeError("model-size rigidity audit failed")

    return ModelSizeRigidityRecord(
        crabb_dimension=crabb_dimension,
        repeated_multiplicity=repeated_multiplicity,
        model_dimension=model_dimension,
        separating_power=exponent,
        repeated_power_norm_square=str(repeated_norm),
        long_model_power_norm_square=str(long_norm),
        squared_norm_gap=str(gap),
        all_checks_passed=verified,
    )


def standard_records() -> list[ModelSizeRigidityRecord]:
    """Return the deterministic exact audit family."""

    return [
        audit_case(crabb_dimension, repeated_multiplicity)
        for crabb_dimension in range(3, 8)
        for repeated_multiplicity in (2, 3)
    ]


def write_records(
    records: list[ModelSizeRigidityRecord],
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
            "experiments/repeated_crabb_model_size_rigidity_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the exact power-separation audit."""

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
