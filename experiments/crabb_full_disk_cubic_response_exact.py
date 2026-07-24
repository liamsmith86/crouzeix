#!/usr/bin/env python3
"""Audit the all-size cubic circular-normal response formula exactly.

For each requested length and deterministic Gaussian-rational
Toeplitz direction, this checker independently runs the complete
characteristic/reversed-Horner/inverse-Riemann series engine through
degree three.  It compares both real normal polarizations with the
closed interval-flux formula and checks that every other circular mode
vanishes.

Finite-size regeneration is an audit of the general path-count proof,
not an interpolation premise for it.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp

from crabb_circular_normal_quadratic_exact import (
    characteristic_series,
    crabb_normal_directions,
    directional_gradient_series,
    disk_model_series,
)
from crabb_full_disk_correction_isometry import plucker_correction
from crabb_full_disk_cubic_response_formula import cubic_normal_response


@dataclass(frozen=True)
class CubicResponseRecord:
    """One exact all-mode response audit."""

    dimension: int
    length: int
    direction_index: int
    active_modes: tuple[int, ...]
    expected_active_modes: tuple[int, ...]
    normal_polarization_count: int
    quadratic_response_zero: bool
    cubic_formula_verified: bool
    inactive_modes_zero: bool
    all_identities_verified: bool


def deterministic_direction(
    length: int,
    direction_index: int,
) -> tuple[sp.Expr, ...]:
    """Return one dense exact complex direction."""

    values = [sp.Integer(0)]
    for offset in range(1, length):
        real_sign = -1 if (direction_index + 1) * offset % 3 == 0 else 1
        imaginary_sign = -1 if (direction_index + offset) % 2 else 1
        values.append(
            real_sign
            * sp.Rational(
                (direction_index + 2) * offset + 1,
                (17 + direction_index) * length,
            )
            + sp.I
            * imaginary_sign
            * sp.Rational(
                (2 * direction_index + 3) * offset + 2,
                (23 + 2 * direction_index) * length,
            )
        )
    return tuple(values)


def audit_direction(
    task: tuple[int, int],
) -> CubicResponseRecord:
    """Compare every circular-normal polarization on one direction."""

    length, direction_index = task
    direction = deterministic_direction(length, direction_index)
    correction = plucker_correction(direction)
    operator, metric = disk_model_series(
        direction,
        order=3,
        correction=correction,
    )
    characteristic = characteristic_series(operator)

    quadratic_zero = True
    formula_verified = True
    inactive_zero = True
    active_modes = []
    polarization_count = 0
    for mode in range(3, length + 2):
        responses = [
            directional_gradient_series(
                operator,
                metric,
                characteristic,
                normal,
            )
            for normal in crabb_normal_directions([metric[0]], mode)
        ]
        if len(responses) != 2:
            raise RuntimeError(
                f"mode {mode} has {len(responses)} normal polarizations"
            )
        polarization_count += len(responses)
        quadratic_zero &= all(
            sp.simplify(response[2]) == 0 for response in responses
        )
        actual = sp.expand(responses[0][3] + sp.I * responses[1][3])
        predicted = cubic_normal_response(direction, mode)
        formula_verified &= sp.simplify(actual - predicted) == 0
        if actual != 0:
            active_modes.append(mode)
        if mode > length - 3:
            inactive_zero &= actual == 0

    expected = tuple(range(3, length - 2))
    verified = bool(
        quadratic_zero
        and formula_verified
        and inactive_zero
        and tuple(active_modes) == expected
    )
    if not verified:
        raise RuntimeError(
            f"cubic response audit failed at L={length}, "
            f"direction={direction_index}"
        )
    return CubicResponseRecord(
        dimension=length + 1,
        length=length,
        direction_index=direction_index,
        active_modes=tuple(active_modes),
        expected_active_modes=expected,
        normal_polarization_count=polarization_count,
        quadratic_response_zero=quadratic_zero,
        cubic_formula_verified=formula_verified,
        inactive_modes_zero=inactive_zero,
        all_identities_verified=verified,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-length", type=int, default=6)
    parser.add_argument("--maximum-length", type=int, default=11)
    parser.add_argument("--direction-count", type=int, default=2)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run and serialize the exact audits."""

    args = parse_args()
    if args.minimum_length < 6 or args.maximum_length < args.minimum_length:
        raise ValueError("require 6 <= minimum length <= maximum length")
    if args.direction_count < 1 or args.workers < 1:
        raise ValueError("direction-count and workers must be positive")

    tasks = [
        (length, direction_index)
        for length in range(args.minimum_length, args.maximum_length + 1)
        for direction_index in range(args.direction_count)
    ]
    records = []
    if args.workers == 1:
        for completed, task in enumerate(tasks, 1):
            records.append(audit_direction(task))
            print(f"verified {completed}/{len(tasks)} directions", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(audit_direction, task) for task in tasks]
            for completed, future in enumerate(as_completed(futures), 1):
                records.append(future.result())
                print(
                    f"verified {completed}/{len(tasks)} directions",
                    flush=True,
                )

    records.sort(key=lambda record: (record.length, record.direction_index))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(asdict(record), sort_keys=True) for record in records]
    args.output.write_text(
        "".join(f"{line}\n" for line in lines),
        encoding="utf-8",
    )
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
