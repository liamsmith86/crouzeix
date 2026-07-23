#!/usr/bin/env python3
"""Audit L130's outside-critical-factor fiber trace identity.

For deterministic real finite Blaschke products, factor

    N'D - ND' = kappa Q Q#

with all finite zeros of ``Q`` outside the disk.  The script compares the
unweighted traces of ``p/Q`` on several fibers of ``B=N/D`` and checks their
common value against the independent Hardy-space inner product

    <p/D, Q/D> / ||Q/D||^2.

This is a floating regression for the residue proof in L130.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from blaschke_stein_composition import real_blaschke_critical_factor


DEFAULT_MINIMUM_DEGREE = 2
DEFAULT_MAXIMUM_DEGREE = 8
DEFAULT_QUADRATURE_NODES = 16_384
FIBER_VALUES = (
    0.0 + 0.0j,
    0.13 + 0.07j,
    -0.22 + 0.11j,
    0.31 - 0.19j,
)


@dataclass(frozen=True)
class CriticalTraceRecord:
    degree: int
    minimum_zero_separation: float
    factorization_residual: float
    minimum_outer_critical_modulus: float
    fiber_trace_variation: float
    hardy_pairing_residual: float


def deterministic_zeros(degree: int) -> tuple[complex, ...]:
    """Generate separated real zeros strictly inside the disk."""

    rng = np.random.default_rng(81_000 + degree)
    for _ in range(10_000):
        zeros = np.sort(rng.uniform(-0.78, 0.78, degree))
        if np.min(np.diff(zeros)) > 0.06:
            return tuple(complex(zero) for zero in zeros)
    raise RuntimeError("failed to generate separated deterministic zeros")


def pad_to_same_length(
    first: np.ndarray,
    second: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Pad two ascending coefficient arrays to a common length."""

    length = max(len(first), len(second))
    first_padded = np.zeros(length, dtype=complex)
    second_padded = np.zeros(length, dtype=complex)
    first_padded[: len(first)] = first
    second_padded[: len(second)] = second
    return first_padded, second_padded


def rational_value(
    numerator: np.ndarray,
    denominator: np.ndarray,
    points: np.ndarray,
) -> np.ndarray:
    """Evaluate one ascending-coefficient rational function."""

    return (
        np.polynomial.polynomial.polyval(points, numerator)
        / np.polynomial.polynomial.polyval(points, denominator)
    )


def make_record(
    degree: int,
    quadrature_nodes: int,
) -> CriticalTraceRecord:
    """Audit one deterministic real Blaschke product."""

    zeros = deterministic_zeros(degree)
    factorization = real_blaschke_critical_factor(zeros)
    numerator = factorization.numerator
    denominator = factorization.denominator
    outer_factor = factorization.outer_factor

    rng = np.random.default_rng(82_000 + degree)
    polynomial = (
        rng.standard_normal(degree)
        + 1j * rng.standard_normal(degree)
    )
    numerator_padded, denominator_padded = pad_to_same_length(
        numerator,
        denominator,
    )
    fiber_traces: list[complex] = []
    for fiber_value in FIBER_VALUES:
        fiber_polynomial = (
            numerator_padded - fiber_value * denominator_padded
        )
        fiber_polynomial = np.polynomial.polynomial.polytrim(
            fiber_polynomial,
            tol=1e-13 * np.linalg.norm(fiber_polynomial),
        )
        preimages = np.polynomial.polynomial.polyroots(fiber_polynomial)
        fiber_traces.append(
            complex(
                np.mean(
                    rational_value(
                        polynomial,
                        outer_factor,
                        preimages,
                    )
                )
            )
        )
    fiber_trace_variation = float(
        max(abs(value - fiber_traces[0]) for value in fiber_traces)
    )

    angles = 2 * np.pi * np.arange(quadrature_nodes) / quadrature_nodes
    circle = np.exp(1j * angles)
    function_values = rational_value(
        polynomial,
        denominator,
        circle,
    )
    multiplier_values = rational_value(
        outer_factor,
        denominator,
        circle,
    )
    inner_product = np.mean(
        function_values * np.conjugate(multiplier_values)
    )
    multiplier_norm_square = float(
        np.mean(abs(multiplier_values) ** 2)
    )
    hardy_prediction = inner_product / multiplier_norm_square
    hardy_pairing_residual = float(
        abs(fiber_traces[0] - hardy_prediction)
    )

    outer_critical_points = np.polynomial.polynomial.polyroots(
        outer_factor
    )
    minimum_outer_modulus = (
        float(min(abs(point) for point in outer_critical_points))
        if len(outer_critical_points)
        else float("inf")
    )
    record = CriticalTraceRecord(
        degree=degree,
        minimum_zero_separation=float(
            np.min(np.diff(np.sort(np.real(zeros))))
        ),
        factorization_residual=factorization.factorization_residual,
        minimum_outer_critical_modulus=minimum_outer_modulus,
        fiber_trace_variation=fiber_trace_variation,
        hardy_pairing_residual=hardy_pairing_residual,
    )
    if record.factorization_residual > 2e-10:
        raise AssertionError("the critical Wronskian factorization failed")
    if record.minimum_outer_critical_modulus <= 1.0 + 1e-8:
        raise AssertionError("the selected critical factor was not outer")
    if record.fiber_trace_variation > 2e-10:
        raise AssertionError("the unweighted fiber trace was not constant")
    if record.hardy_pairing_residual > 2e-10:
        raise AssertionError("the Hardy pairing did not recover the trace")
    return record


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--minimum-degree",
        type=int,
        default=DEFAULT_MINIMUM_DEGREE,
    )
    parser.add_argument(
        "--maximum-degree",
        type=int,
        default=DEFAULT_MAXIMUM_DEGREE,
    )
    parser.add_argument(
        "--quadrature-nodes",
        type=int,
        default=DEFAULT_QUADRATURE_NODES,
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the deterministic critical-trace audit."""

    args = parse_args()
    records = [
        make_record(degree, args.quadrature_nodes)
        for degree in range(
            args.minimum_degree,
            args.maximum_degree + 1,
        )
    ]
    lines = [
        json.dumps(asdict(record), sort_keys=True)
        for record in records
    ]
    for line in lines:
        print(line, flush=True)
    if args.output is not None:
        args.output.write_text(
            "".join(f"{line}\n" for line in lines),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
