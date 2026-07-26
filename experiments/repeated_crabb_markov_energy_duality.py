#!/usr/bin/env python3
"""Audit the exact Dirichlet-energy dual of Markov endpoint repair."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from crabb_block_hardy_equality import format_float


Matrix = np.ndarray


@dataclass(frozen=True)
class MarkovEnergyDualityRecord:
    """One primal/dual audit for the depolarizing Markov Laplacian."""

    positive_endpoint: str
    negative_endpoint: str
    expected_energy: str
    primal_energy: str
    dual_energy: str
    ratio_energy: str
    primal_error: str
    dual_error: str
    ratio_error: str
    corrected_endpoint_maximum: str
    all_checks_passed: bool


def depolarizing_laplacian(matrix):
    """Return ``matrix - trace(matrix) I / 2``."""

    return matrix - cp.trace(matrix) * np.eye(2) / 2


def audit_face(
    positive_endpoint: float,
    negative_endpoint: float,
) -> MarkovEnergyDualityRecord:
    """Compare the exact energy with independent primal and dual SDPs."""

    if not (
        positive_endpoint > 0
        and positive_endpoint + negative_endpoint < 0
    ):
        raise ValueError("the scalar fixed-point separator must be negative")

    endpoint = np.diag((positive_endpoint, negative_endpoint))
    expected = positive_endpoint**2 / 32

    correction = cp.Variable((2, 2), symmetric=True)
    laplacian_correction = depolarizing_laplacian(correction)
    primal = cp.Problem(
        cp.Minimize(cp.sum_squares(laplacian_correction)),
        [endpoint + 8 * laplacian_correction << 0],
    )
    primal.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-11,
        tol_feas=1e-11,
    )

    separator = cp.Variable((2, 2), symmetric=True)
    laplacian_separator = depolarizing_laplacian(separator)
    dual = cp.Problem(
        cp.Maximize(
            cp.trace(separator @ endpoint)
            - 16 * cp.sum_squares(laplacian_separator)
        ),
        [separator >> 0],
    )
    dual.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-11,
        tol_feas=1e-11,
    )

    primal_value = float(primal.value)
    dual_value = float(dual.value)
    ratio_value = positive_endpoint**2 / (
        64 * 0.5
    )
    corrected = (
        endpoint
        + 8
        * np.asarray(laplacian_correction.value, dtype=float)
    )
    corrected_maximum = float(np.linalg.eigvalsh(corrected)[-1])

    primal_error = abs(primal_value - expected)
    dual_error = abs(dual_value - expected)
    ratio_error = abs(ratio_value - expected)
    tolerance = 2e-8
    verified = bool(
        primal.status == cp.OPTIMAL
        and dual.status == cp.OPTIMAL
        and primal_error < tolerance
        and dual_error < tolerance
        and ratio_error < tolerance
        and corrected_maximum < tolerance
    )
    if not verified:
        raise RuntimeError(
            "Markov energy duality audit failed: "
            f"primal={primal_value:.12e}, "
            f"dual={dual_value:.12e}, "
            f"expected={expected:.12e}, "
            f"endpoint={corrected_maximum:.12e}"
        )

    return MarkovEnergyDualityRecord(
        positive_endpoint=format_float(positive_endpoint),
        negative_endpoint=format_float(negative_endpoint),
        expected_energy=format_float(expected),
        primal_energy=format_float(primal_value),
        dual_energy=format_float(dual_value),
        ratio_energy=format_float(ratio_value),
        primal_error=format_float(primal_error),
        dual_error=format_float(dual_error),
        ratio_error=format_float(ratio_error),
        corrected_endpoint_maximum=format_float(corrected_maximum),
        all_checks_passed=verified,
    )


def standard_records() -> list[MarkovEnergyDualityRecord]:
    """Return deterministic scalar-separator test faces."""

    return [
        audit_face(0.25, -0.75),
        audit_face(0.8, -1.1),
        audit_face(1.7, -2.4),
    ]


def write_records(
    records: list[MarkovEnergyDualityRecord],
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
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_markov_energy_duality_s70225.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete deterministic audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
