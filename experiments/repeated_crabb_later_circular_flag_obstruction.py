#!/usr/bin/env python3
"""Audit the obstruction to a Schur-only later circular-normal flag.

The first joint face, its reducing kernel, the cross-order gain, and
the transported Schur square do not determine a later kernel-block
coefficient.  This checker regenerates the exact polynomial
two-by-two certificate in L319.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import sympy as sp


@dataclass(frozen=True)
class LaterCircularFlagObstructionRecord:
    """One exact Schur-only obstruction audit."""

    record_kind: str
    first_face: str
    cross_order_gain: int
    transported_square: str
    disk_schur_quotient: str
    joint_schur_quotient: str
    disk_determinant: str
    joint_determinant: str
    correction: str
    positive_test_parameters: tuple[str, ...]
    all_checks_passed: bool


def polynomial_order(expression: sp.Expr, parameter: sp.Symbol) -> int:
    """Return the valuation of a nonzero polynomial at the origin."""

    polynomial = sp.Poly(sp.expand(expression), parameter)
    return min(monomial[0] for monomial, _ in polynomial.terms())


def exact_record() -> LaterCircularFlagObstructionRecord:
    """Build and verify L319's exact two-by-two certificate."""

    parameter = sp.symbols("s", real=True)
    active = parameter**2
    cross = parameter**3
    disk_kernel = parameter**4 + parameter**6
    joint_kernel = parameter**4 - parameter**6
    disk = sp.Matrix([[active, cross], [cross, disk_kernel]])
    joint = sp.Matrix([[active, cross], [cross, joint_kernel]])

    transported = sp.cancel(cross**2 / active)
    disk_schur = sp.factor(disk_kernel - transported)
    joint_schur = sp.factor(joint_kernel - transported)
    disk_determinant = sp.factor(disk.det())
    joint_determinant = sp.factor(joint.det())
    correction = sp.simplify(joint - disk)
    first_face = disk.applyfunc(
        lambda entry: sp.limit(entry / parameter**2, parameter, 0)
    )

    test_parameters = (
        sp.Rational(1, 7),
        sp.Rational(1, 3),
        sp.Rational(1, 2),
    )
    signs_hold = all(
        disk_schur.subs(parameter, value) > 0
        and joint_schur.subs(parameter, value) < 0
        for value in test_parameters
    )
    verified = bool(
        first_face == sp.diag(1, 0)
        and polynomial_order(cross, parameter)
        - polynomial_order(active, parameter)
        == 1
        and transported == parameter**4
        and disk_schur == parameter**6
        and joint_schur == -parameter**6
        and disk_determinant == parameter**8
        and joint_determinant == -parameter**8
        and correction == sp.diag(0, -2 * parameter**6)
        and signs_hold
    )
    if not verified:
        raise RuntimeError("the exact later-circular flag obstruction failed")

    return LaterCircularFlagObstructionRecord(
        record_kind="exact_schur_only_later_flag_obstruction",
        first_face=str(first_face),
        cross_order_gain=1,
        transported_square=str(transported),
        disk_schur_quotient=str(disk_schur),
        joint_schur_quotient=str(joint_schur),
        disk_determinant=str(disk_determinant),
        joint_determinant=str(joint_determinant),
        correction=str(correction),
        positive_test_parameters=tuple(
            str(value) for value in test_parameters
        ),
        all_checks_passed=verified,
    )


def write_records(
    path: Path,
    records: Sequence[LaterCircularFlagObstructionRecord],
) -> None:
    """Write deterministic JSON Lines output atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Regenerate the exact obstruction dataset."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_later_circular_flag_obstruction_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    record = exact_record()
    write_records(args.output, [record])
    print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
