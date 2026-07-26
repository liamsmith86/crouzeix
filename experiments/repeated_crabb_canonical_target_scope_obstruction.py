#!/usr/bin/env python3
"""Audit the apex mismatch between two repeated-Crabb metric branches."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


Matrix = sp.Matrix


@dataclass(frozen=True)
class ScopeObstructionRecord:
    """One exact terminal-unitary comparison."""

    multiplicity: int
    unitary_kind: str
    left_unital_residual_count: int
    right_unital_residual_count: int
    flux_residual_count: int
    l212_face_residual_count: int
    repair_face_residual_count: int
    comparison_gap_residual_count: int
    comparison_violation: str
    l212_response_residual_count: int
    all_checks_passed: bool


def adjoint(matrix: Matrix) -> Matrix:
    """Return the conjugate transpose."""

    return matrix.conjugate().T


def residual_count(matrix: Matrix) -> int:
    """Count nonzero entries after exact simplification."""

    return sum(
        int(sp.simplify(entry) != 0)
        for entry in matrix
    )


def terminal_unitary(size: int) -> tuple[str, Matrix]:
    """Return a deterministic exact terminal unitary."""

    if size == 1:
        return "scalar_phase", Matrix([[sp.I]])
    unitary = sp.zeros(size, size)
    for column in range(size):
        unitary[(column + 1) % size, column] = (
            -1 if column == 0 else 1
        )
    return "signed_cycle", unitary


def audit_case(size: int) -> ScopeObstructionRecord:
    """Audit one repeated monomial apex exactly."""

    unitary_kind, active = terminal_unitary(size)
    identity = sp.eye(size)

    left_unital = active * adjoint(active) - identity
    right_unital = adjoint(active) * active - identity

    def phi(matrix: Matrix) -> Matrix:
        return active * matrix * adjoint(active)

    active_left_gram = active * adjoint(active)
    active_right_gram = adjoint(active) * active
    flux = phi(active_right_gram) - active_left_gram

    l212_face = (
        12 * active_left_gram
        - 28 * phi(active_right_gram)
    )
    repair_face = -12 * active_left_gram
    comparison_gap = repair_face - l212_face
    expected_l212 = -16 * identity
    expected_repair = -12 * identity
    expected_gap = 4 * identity
    l212_response = 28 * flux

    checks_passed = bool(
        residual_count(left_unital) == 0
        and residual_count(right_unital) == 0
        and residual_count(flux) == 0
        and residual_count(l212_face - expected_l212) == 0
        and residual_count(repair_face - expected_repair) == 0
        and residual_count(comparison_gap - expected_gap) == 0
        and residual_count(l212_response) == 0
    )
    if not checks_passed:
        raise RuntimeError(
            f"scope-obstruction audit failed at multiplicity {size}"
        )

    return ScopeObstructionRecord(
        multiplicity=size,
        unitary_kind=unitary_kind,
        left_unital_residual_count=residual_count(left_unital),
        right_unital_residual_count=residual_count(right_unital),
        flux_residual_count=residual_count(flux),
        l212_face_residual_count=residual_count(
            l212_face - expected_l212
        ),
        repair_face_residual_count=residual_count(
            repair_face - expected_repair
        ),
        comparison_gap_residual_count=residual_count(
            comparison_gap - expected_gap
        ),
        comparison_violation=str(comparison_gap[0, 0]),
        l212_response_residual_count=residual_count(l212_response),
        all_checks_passed=checks_passed,
    )


def standard_records() -> list[ScopeObstructionRecord]:
    """Return all exact terminal-unitary audits."""

    return [audit_case(size) for size in range(1, 5)]


def write_records(
    records: list[ScopeObstructionRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(
        json.dumps(asdict(record), sort_keys=True) + "\n"
        for record in records
    )
    output.write_text(payload, encoding="utf-8")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the exact audit."""

    args = parse_args()
    records = standard_records()
    digest = None
    if args.output is not None:
        digest = write_records(records, args.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "dataset_sha256": digest,
                "records": [asdict(record) for record in records],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
