#!/usr/bin/env python3
"""Audit the terminal two-coefficient sixth-order Schur face.

For length ``L >= 6``, put ``k=L-3`` and restrict the recentered
full-disk path to Toeplitz directions supported only on its last two
coefficients ``a,b``.  The exact endpoint and characteristic-series
engines are compared with the closed formulas

    base = 256 k^2 |a|^4 |b|^2 / (k+2)^2,

    G = 16 k (4k-1) a^2 conjugate(b) / (L^2 (k+2)),

and the resulting completed-square ratio.  This is an exact
finite-size audit of the all-size formulas recorded in A127.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_full_disk_base_jet_exact import endpoint_delta_series
from crabb_full_disk_cubic_normal_exact import (
    audit_direction,
    highest_mode_cubic_response,
)


@dataclass(frozen=True)
class TerminalSixthRecord:
    """One exact terminal two-coefficient audit."""

    dimension: int
    length: int
    support_mode: int
    penultimate_coefficient: str
    terminal_coefficient: str
    base_deficit: str
    predicted_base_deficit: str
    cubic_response: str
    predicted_cubic_response: str
    schur_ratio: str
    predicted_schur_ratio: str
    symbolic_base_identity_verified: bool | None
    all_identities_verified: bool


def terminal_direction(length: int) -> tuple[sp.Expr, ...]:
    """Return a fixed Gaussian-rational terminal direction."""

    if length < 6:
        raise ValueError("terminal sixth-order audit requires length >= 6")
    direction = [sp.Integer(0)] * length
    direction[-2] = sp.Rational(2, 7) + sp.I * sp.Rational(1, 11)
    direction[-1] = sp.Rational(3, 10) - sp.I * sp.Rational(2, 13)
    return tuple(direction)


def predicted_base_deficit(
    length: int,
    penultimate: sp.Expr,
    terminal: sp.Expr,
) -> sp.Expr:
    """Return the closed terminal sixth-order base deficit."""

    mode = length - 3
    return sp.factor(
        sp.simplify(
            sp.Rational(256 * mode**2, (mode + 2) ** 2)
            * (penultimate * sp.conjugate(penultimate)) ** 2
            * terminal
            * sp.conjugate(terminal)
        )
    )


def predicted_cubic_response(
    length: int,
    penultimate: sp.Expr,
    terminal: sp.Expr,
) -> sp.Expr:
    """Return the closed highest-mode terminal cubic response."""

    mode = length - 3
    return sp.expand(
        sp.Rational(
            16 * mode * (4 * mode - 1),
            length**2 * (mode + 2),
        )
        * penultimate**2
        * sp.conjugate(terminal)
    )


def predicted_ratio(length: int) -> sp.Expr:
    """Return the terminal completed-square gain/base ratio."""

    mode = length - 3
    return sp.factor(
        sp.Rational(
            6 * (4 * mode - 1) ** 2,
            6 * (4 * mode - 1) ** 2 + 169 * mode * (mode - 1) * (mode - 2),
        )
    )


def verify_symbolic_base_identity(length: int) -> bool:
    """Verify the base formula with unrestricted complex symbols."""

    penultimate, terminal = sp.symbols(
        "penultimate terminal",
    )
    direction = [sp.Integer(0)] * length
    direction[-2] = penultimate
    direction[-1] = terminal
    actual_delta = endpoint_delta_series(
        tuple(direction),
        order=6,
    )[6]
    predicted_delta = (
        -predicted_base_deficit(
            length,
            penultimate,
            terminal,
        )
        / 2
    )
    return bool(sp.expand(actual_delta - predicted_delta) == 0)


def make_record(length: int) -> TerminalSixthRecord:
    """Build and verify one exact terminal record."""

    direction = terminal_direction(length)
    penultimate = direction[-2]
    terminal = direction[-1]
    audit = audit_direction(direction, direction_index=0)
    mode = length - 3
    if audit.active_cubic_modes != (mode,):
        raise RuntimeError(
            f"unexpected terminal active modes at length {length}: "
            f"{audit.active_cubic_modes}"
        )
    actual_cubic = sp.factor(
        sp.sympify(audit.cubic_responses[0])
        + sp.I * sp.sympify(audit.cubic_responses[1])
    )
    closed_base = predicted_base_deficit(
        length,
        penultimate,
        terminal,
    )
    closed_cubic = predicted_cubic_response(
        length,
        penultimate,
        terminal,
    )
    factored_cubic = highest_mode_cubic_response(direction)
    closed_ratio = predicted_ratio(length)
    actual_base = sp.sympify(audit.base_deficit_sixth_coefficient)
    actual_ratio = sp.sympify(audit.gain_over_base_deficit)
    symbolic_verified = verify_symbolic_base_identity(length) if length <= 8 else None
    verified = (
        all(
            sp.simplify(residual) == 0
            for residual in (
                actual_base - closed_base,
                actual_cubic - closed_cubic,
                actual_cubic - factored_cubic,
                actual_ratio - closed_ratio,
            )
        )
        and symbolic_verified is not False
    )
    if not verified:
        raise RuntimeError(f"terminal formulas failed at length {length}")
    return TerminalSixthRecord(
        dimension=length + 1,
        length=length,
        support_mode=mode,
        penultimate_coefficient=str(penultimate),
        terminal_coefficient=str(terminal),
        base_deficit=str(sp.factor(actual_base)),
        predicted_base_deficit=str(closed_base),
        cubic_response=str(actual_cubic),
        predicted_cubic_response=str(sp.factor(closed_cubic)),
        schur_ratio=str(sp.factor(actual_ratio)),
        predicted_schur_ratio=str(closed_ratio),
        symbolic_base_identity_verified=symbolic_verified,
        all_identities_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=6)
    parser.add_argument("--maximum-length", type=int, default=15)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run the exact terminal audit."""

    args = parse_args()
    if args.minimum_length < 6 or args.maximum_length < args.minimum_length:
        raise ValueError("require 6 <= minimum length <= maximum length")
    records = [
        make_record(length)
        for length in range(
            args.minimum_length,
            args.maximum_length + 1,
        )
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(asdict(record), sort_keys=True) for record in records]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
