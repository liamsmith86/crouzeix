#!/usr/bin/env python3
"""Audit the nilpotent specialization of the Gau--Wu model block."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class ScalarModelSplitRecord:
    """One exact finite-model specialization."""

    dimension: int
    length: int
    model_equals_crabb_transpose: bool
    lower_power_norm_squares: tuple[str, ...]
    terminal_power_norm_square: str
    next_power_vanishes: bool
    all_checks_passed: bool


def forward_model_shift(dimension: int) -> sp.Matrix:
    """Return the compression of the forward shift to polynomials."""

    shift = sp.zeros(dimension)
    for column in range(dimension - 1):
        shift[column + 1, column] = 1
    return shift


def crabb_superdiagonal_shift(dimension: int) -> sp.Matrix:
    """Return the repository's superdiagonal Crabb convention."""

    crabb = sp.zeros(dimension)
    for index in range(dimension - 1):
        weight = (
            sp.sqrt(2)
            if index in {0, dimension - 2}
            else sp.Integer(1)
        )
        crabb[index, index + 1] = weight
    return crabb


def squared_operator_norm_of_power(
    matrix: sp.Matrix,
    exponent: int,
) -> sp.Expr:
    """Return the exact squared norm for this weighted-shift power."""

    power = matrix**exponent
    gram = sp.simplify(power.conjugate().T * power)
    if gram == sp.zeros(matrix.rows):
        return sp.Integer(0)
    if not gram.is_diagonal():
        raise RuntimeError("weighted-shift power Gram is not diagonal")
    return max(sp.simplify(gram[index, index]) for index in range(matrix.rows))


def audit_dimension(dimension: int) -> ScalarModelSplitRecord:
    """Verify the Gau--Wu endpoint scaling at one Crabb length."""

    length = dimension - 1
    shift = forward_model_shift(dimension)
    endpoint_scaling = sp.diag(
        1 / sp.sqrt(2),
        *([sp.Integer(1)] * (dimension - 2)),
        sp.sqrt(2),
    )
    model = sp.simplify(
        endpoint_scaling * shift * endpoint_scaling.inv()
    )
    crabb = crabb_superdiagonal_shift(dimension)

    lower_norms = tuple(
        squared_operator_norm_of_power(model, exponent)
        for exponent in range(1, length)
    )
    terminal_norm = squared_operator_norm_of_power(model, length)
    next_power = sp.simplify(model ** (length + 1))

    model_match = bool(sp.simplify(model.T - crabb) == sp.zeros(dimension))
    verified = bool(
        model_match
        and all(value == 2 for value in lower_norms)
        and terminal_norm == 4
        and next_power == sp.zeros(dimension)
    )
    if not verified:
        raise RuntimeError(f"scalar model split audit failed in dimension {dimension}")

    return ScalarModelSplitRecord(
        dimension=dimension,
        length=length,
        model_equals_crabb_transpose=model_match,
        lower_power_norm_squares=tuple(str(value) for value in lower_norms),
        terminal_power_norm_square=str(terminal_norm),
        next_power_vanishes=bool(next_power == sp.zeros(dimension)),
        all_checks_passed=verified,
    )


def standard_records() -> list[ScalarModelSplitRecord]:
    """Return the deterministic exact audit family."""

    return [audit_dimension(dimension) for dimension in range(3, 10)]


def write_records(
    records: list[ScalarModelSplitRecord],
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
            "experiments/repeated_crabb_scalar_model_split_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the exact finite-model audit."""

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
