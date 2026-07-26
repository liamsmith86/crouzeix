#!/usr/bin/env python3
"""Audit the joint channel/disk response factorization."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class JointFactorRecord:
    """One exact invariant-monomial and square-completion audit."""

    flag_length: int
    maximum_degree: int
    invariant_monomial_count: int
    forbidden_monomial_count: int
    young_residual: str
    completion_residual: str
    all_checks_passed: bool


def invariant_exponents(
    variable_count: int,
    maximum_degree: int,
) -> list[tuple[int, ...]]:
    """Enumerate nonconstant monomials invariant under joint sign."""

    exponents: list[tuple[int, ...]] = []

    def visit(prefix: tuple[int, ...], remaining: int) -> None:
        if len(prefix) == variable_count:
            degree = sum(prefix)
            if 0 < degree <= maximum_degree and degree % 2 == 0:
                exponents.append(prefix)
            return
        for exponent in range(remaining + 1):
            visit(prefix + (exponent,), remaining - exponent)

    visit((), maximum_degree)
    return exponents


def audit_case(flag_length: int, maximum_degree: int) -> JointFactorRecord:
    """Verify parity and the two exact completion identities."""

    # Two transfer cross variables and ``flag_length`` accumulated
    # disk-square variables all reverse sign in the defect-frame gauge.
    exponents = invariant_exponents(2 + flag_length, maximum_degree)
    forbidden = [
        exponent
        for exponent in exponents
        if sum(exponent) < 2
    ]

    delta, disk = sp.symbols("delta disk", nonnegative=True)
    young_residual = sp.expand(
        delta + disk - 2 * sp.sqrt(delta * disk)
    )
    # This is (sqrt(delta)-sqrt(disk))^2.
    young_check = sp.simplify(
        young_residual
        - (sp.sqrt(delta) - sp.sqrt(disk)) ** 2
    )

    gamma = sp.Rational(7, 3)
    response = sp.symbols("response", real=True)
    normal = sp.symbols("normal", real=True)
    completed = sp.expand(
        -gamma * normal**2 / 2
        + 2 * normal * response
        - (
            -gamma
            * (normal - 2 * response / gamma) ** 2
            / 2
            + 2 * response**2 / gamma
        )
    )

    verified = bool(
        exponents
        and not forbidden
        and young_check == 0
        and completed == 0
    )
    if not verified:
        raise RuntimeError("joint channel/disk factor audit failed")

    return JointFactorRecord(
        flag_length=flag_length,
        maximum_degree=maximum_degree,
        invariant_monomial_count=len(exponents),
        forbidden_monomial_count=len(forbidden),
        young_residual=str(young_residual),
        completion_residual=str(completed),
        all_checks_passed=verified,
    )


def standard_records() -> list[JointFactorRecord]:
    """Return deterministic finite-flag audits."""

    return [
        audit_case(flag_length, maximum_degree)
        for flag_length in range(1, 6)
        for maximum_degree in (4, 6, 8)
    ]


def write_records(
    records: list[JointFactorRecord],
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
            "experiments/"
            "repeated_crabb_joint_channel_disk_factor_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the deterministic audit."""

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
