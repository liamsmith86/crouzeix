#!/usr/bin/env python3
"""Stress-test the candidate uniform Crabb equality-tube estimate.

This is a falsification probe, not a proof.  It samples deterministic
one- and mixed-grade real phase-palindromic equality branches, applies
the exact disk chart and elliptic pullback, and minimizes the local
rank-one Stein envelope.  The candidate L147 half-face estimate is

    4 - Gamma_eq >= 8 c^(2L) + 32 sum_k |u_k|^2 c^(2k).

The reported ratio should therefore stay above one in a sufficiently
small neighbourhood if the candidate reflected-ideal theorem is true.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from mpmath import mp
import numpy as np

from crabb_palindromic_elliptic_face import (
    disk_matrix,
    ellipse_pullback,
    rank_one_envelope,
)


DEFAULT_ELLIPSE_PARAMETERS = (0.08, 0.12, 0.18)
DEFAULT_AMPLITUDES = (0.01, 0.03, 0.06)


@dataclass(frozen=True)
class UniformRemainderRecord:
    """One deterministic candidate-tube sample."""

    dimension: int
    length: int
    direction_name: str
    low_grade_direction: tuple[float, ...]
    ellipse_parameter: float
    amplitude: float
    rank_one_bound: float
    candidate_face_size: float
    domination_ratio: float


def phase_palindromic_coefficients(
    length: int,
    low_grades: np.ndarray,
) -> np.ndarray:
    """Lift low-grade real coordinates to a phase-one coefficient vector."""

    result = np.zeros(length - 1)
    for grade, value in enumerate(low_grades, start=1):
        result[grade - 1] = value
        result[length - grade - 1] = value
    return result


def deterministic_directions(
    length: int,
) -> tuple[tuple[str, np.ndarray], ...]:
    """Return normalized one-grade and mixed-grade directions."""

    grade_count = length // 2
    directions: list[tuple[str, np.ndarray]] = []
    for grade in range(1, grade_count + 1):
        direction = np.zeros(grade_count)
        direction[grade - 1] = 1
        directions.append((f"grade_{grade}", direction))
    if grade_count > 1:
        all_positive = np.ones(grade_count)
        all_positive /= np.linalg.norm(all_positive)
        alternating = np.asarray(
            [(-1.0) ** grade for grade in range(grade_count)]
        )
        alternating /= np.linalg.norm(alternating)
        directions.extend(
            (
                ("mixed_positive", all_positive),
                ("mixed_alternating", alternating),
            )
        )
    return tuple(directions)


def make_records(
    dimension: int,
    ellipse_parameters: tuple[float, ...],
    amplitudes: tuple[float, ...],
) -> list[UniformRemainderRecord]:
    """Evaluate all deterministic branches in one dimension."""

    length = dimension - 1
    records: list[UniformRemainderRecord] = []
    for direction_name, low_direction in deterministic_directions(length):
        coefficients = phase_palindromic_coefficients(
            length,
            low_direction,
        )
        for ellipse_parameter in ellipse_parameters:
            initial: np.ndarray | None = None
            for amplitude in amplitudes:
                operator = ellipse_pullback(
                    disk_matrix(
                        dimension,
                        amplitude,
                        coefficients,
                    ),
                    ellipse_parameter,
                )
                bound, initial = rank_one_envelope(operator, initial)
                weighted_amplitude = sum(
                    abs(amplitude * value) ** 2
                    * ellipse_parameter ** (2 * grade)
                    for grade, value in enumerate(
                        low_direction,
                        start=1,
                    )
                )
                face_size = (
                    8 * ellipse_parameter ** (2 * length)
                    + 32 * weighted_amplitude
                )
                ratio = (4 - bound) / face_size
                if ratio <= 0:
                    raise AssertionError(
                        "the candidate equality tube lost its sign"
                    )
                records.append(
                    UniformRemainderRecord(
                        dimension=dimension,
                        length=length,
                        direction_name=direction_name,
                        low_grade_direction=tuple(
                            float(value) for value in low_direction
                        ),
                        ellipse_parameter=ellipse_parameter,
                        amplitude=amplitude,
                        rank_one_bound=bound,
                        candidate_face_size=face_size,
                        domination_ratio=ratio,
                    )
                )
    return records


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-size", type=int, default=3)
    parser.add_argument("--maximum-size", type=int, default=8)
    parser.add_argument(
        "--ellipse-parameters",
        type=float,
        nargs="+",
        default=DEFAULT_ELLIPSE_PARAMETERS,
    )
    parser.add_argument(
        "--amplitudes",
        type=float,
        nargs="+",
        default=DEFAULT_AMPLITUDES,
    )
    parser.add_argument("--precision", type=int, default=50)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    """Run the probe and optionally persist its JSONL records."""

    args = parse_args()
    if args.minimum_size < 3 or args.maximum_size < args.minimum_size:
        raise ValueError("size range must satisfy 3 <= minimum <= maximum")
    if any(amplitude <= 0 for amplitude in args.amplitudes):
        raise ValueError("amplitudes must be positive")
    if any(
        not 0 < parameter < 1
        for parameter in args.ellipse_parameters
    ):
        raise ValueError("ellipse parameters must lie in (0, 1)")
    mp.dps = args.precision

    records = [
        record
        for dimension in range(
            args.minimum_size,
            args.maximum_size + 1,
        )
        for record in make_records(
            dimension,
            tuple(args.ellipse_parameters),
            tuple(args.amplitudes),
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
