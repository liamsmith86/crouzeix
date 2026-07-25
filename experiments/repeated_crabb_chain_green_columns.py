#!/usr/bin/env python3
"""Audit the exact incoming and outgoing Green columns of a clean delay.

The full block inverse after a complete delay uses only the last row of
the chain resolvent and the chain response to its tail coupling.  Both
have closed continuant formulas.  On the Joukowski contour those formulas
reduce to monomials times ``(1+t**j)/(1+t**r)``.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_multidelay_terminal_resolvent import (
    chain_pencil,
    continuant,
)


@dataclass(frozen=True)
class ChainGreenColumnsRecord:
    """One exact-form chain Green-column audit."""

    delay: int
    maximum_continuant_row_error: str
    maximum_continuant_column_error: str
    maximum_joukowski_row_error: str
    maximum_joukowski_column_error: str
    maximum_self_energy_error: str
    all_checks_passed: bool


def audit_delay(delay: int) -> ChainGreenColumnsRecord:
    """Audit one delay length at deterministic complex contour points."""

    row_error = 0.0
    column_error = 0.0
    joukowski_row_error = 0.0
    joukowski_column_error = 0.0
    self_energy_error = 0.0

    for c in (0.03, 0.17, 0.41):
        chain, balance, last = chain_pencil(delay, 1, c)
        for zeta in (
            1.7 * np.exp(0.23j),
            1.9 * np.exp(1.17j),
            2.2 * np.exp(-0.71j),
        ):
            x = 1 / zeta
            rho = c / zeta
            t = c / zeta**2
            z = zeta + c / zeta
            green = np.linalg.inv(z * np.eye(delay) - chain)
            row = (last.conj().T @ green).reshape(-1)
            column = (
                green @ (c * balance @ last)
            ).reshape(-1)

            delta_delay = continuant(delay, z, c)
            predicted_row = np.array(
                [
                    continuant(index, z, c) / delta_delay
                    for index in range(delay)
                ],
                dtype=complex,
            )
            predicted_column = np.array(
                [
                    (
                        2 * c**delay / delta_delay
                        if index == 0
                        else c ** (delay - index)
                        * continuant(index, z, c)
                        / delta_delay
                    )
                    for index in range(delay)
                ],
                dtype=complex,
            )
            row_error = max(
                row_error,
                float(np.linalg.norm(row - predicted_row)),
            )
            column_error = max(
                column_error,
                float(np.linalg.norm(column - predicted_column)),
            )

            predicted_joukowski_row = np.empty(delay, dtype=complex)
            predicted_joukowski_column = np.empty(
                delay,
                dtype=complex,
            )
            for index in range(delay):
                if index == 0:
                    predicted_joukowski_row[index] = (
                        x**delay / (1 + t**delay)
                    )
                    predicted_joukowski_column[index] = (
                        2 * rho**delay / (1 + t**delay)
                    )
                else:
                    predicted_joukowski_row[index] = (
                        x ** (delay - index)
                        * (1 + t**index)
                        / (1 + t**delay)
                    )
                    predicted_joukowski_column[index] = (
                        rho ** (delay - index)
                        * (1 + t**index)
                        / (1 + t**delay)
                    )
            joukowski_row_error = max(
                joukowski_row_error,
                float(
                    np.linalg.norm(
                        row - predicted_joukowski_row
                    )
                ),
            )
            joukowski_column_error = max(
                joukowski_column_error,
                float(
                    np.linalg.norm(
                        column - predicted_joukowski_column
                    )
                ),
            )
            self_energy_error = max(
                self_energy_error,
                abs(
                    row
                    @ (c * balance @ last).reshape(-1)
                    - (
                        c / zeta
                        + zeta
                        * t**delay
                        * (1 - t)
                        / (1 + t**delay)
                    )
                ),
            )

    tolerance = 8e-11
    verified = max(
        row_error,
        column_error,
        joukowski_row_error,
        joukowski_column_error,
        self_energy_error,
    ) < tolerance
    if not verified:
        raise RuntimeError(
            f"the chain Green-column audit failed at delay {delay}"
        )
    return ChainGreenColumnsRecord(
        delay=delay,
        maximum_continuant_row_error=format_float(row_error),
        maximum_continuant_column_error=format_float(column_error),
        maximum_joukowski_row_error=format_float(
            joukowski_row_error
        ),
        maximum_joukowski_column_error=format_float(
            joukowski_column_error
        ),
        maximum_self_energy_error=format_float(self_energy_error),
        all_checks_passed=verified,
    )


def standard_records(
    maximum_delay: int,
) -> list[ChainGreenColumnsRecord]:
    """Return audits through one maximum delay."""

    return [
        audit_delay(delay)
        for delay in range(1, maximum_delay + 1)
    ]


def write_records(
    records: list[ChainGreenColumnsRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

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
    parser.add_argument("--maximum-delay", type=int, default=8)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_chain_green_columns_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the chain Green-column audits."""

    args = parse_args()
    records = standard_records(args.maximum_delay)
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
