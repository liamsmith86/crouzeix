#!/usr/bin/env python3
"""Audit two-sided channel leakage and scalar response factorization."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_matrix_schur_chart import reconstruct_series


MatrixSeries = list[np.ndarray]


@dataclass(frozen=True)
class ChannelResponseRecord:
    """One exact matrix-inner channel-defect audit."""

    construction_kind: str
    step: str
    state_count: int
    maximum_parseval_error: str
    maximum_column_identity_error: str
    maximum_row_identity_error: str
    maximum_fourier_ratio: str
    maximum_cross_pairing_ratio: str
    maximum_gauge_error: str
    maximum_split_channel_error: str
    all_checks_passed: bool


def orthogonal_vector(vector: np.ndarray) -> np.ndarray:
    """Return the canonical unit vector orthogonal to one qubit."""

    return np.array(
        [-vector[1].conjugate(), vector[0].conjugate()],
        dtype=complex,
    )


def channel_states() -> list[tuple[np.ndarray, np.ndarray]]:
    """Return coordinate and dense deterministic channel frames."""

    root_two = np.sqrt(2)
    right_states = [
        np.array([1, 0], dtype=complex),
        np.array([0, 1], dtype=complex),
        np.array([1, 1j], dtype=complex) / root_two,
        np.array([1, -1], dtype=complex) / root_two,
    ]
    left_states = [
        np.array([1, 0], dtype=complex),
        np.array([0, 1], dtype=complex),
        np.array([1, -1j], dtype=complex) / root_two,
        np.array([1, 1], dtype=complex) / root_two,
    ]
    return list(zip(left_states, right_states, strict=True))


def schur_path(
    kind: str,
    step: float,
    count: int,
) -> MatrixSeries:
    """Return one exact matrix-inner curve through ``z^3 I_2``."""

    zero = np.zeros((2, 2), dtype=complex)
    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    pauli_z = np.diag([1, -1]).astype(complex)
    parameters = [zero.copy() for _ in range(3)]
    if kind == "transverse_pauli":
        parameters[1] = 0.24 * step * pauli_x
        parameters[2] = 0.19 * step * pauli_z
    elif kind == "generic_nonnormal":
        parameters[1] = step * np.array(
            [[0.18, 0.07j], [-0.04, 0.02]],
            dtype=complex,
        )
        parameters[2] = step * np.array(
            [[0.03j, -0.09], [0.05j, -0.13]],
            dtype=complex,
        )
    elif kind == "split_diagonal":
        parameters[1] = (
            step * np.diag([0.18, -0.11]).astype(complex)
        )
        parameters[2] = step * np.diag([0.06j, -0.14j])
    else:
        raise ValueError(f"unknown path kind: {kind}")
    return reconstruct_series(
        parameters,
        np.eye(2, dtype=complex),
        count,
    )


def parseval_error(coefficients: MatrixSeries) -> float:
    """Return the truncated column-Parseval residual."""

    gram = sum(
        coefficient.conj().T @ coefficient
        for coefficient in coefficients
    )
    return float(np.linalg.norm(gram - np.eye(2)))


def fixed_channel_statistics(
    coefficients: MatrixSeries,
    left: np.ndarray,
    right: np.ndarray,
) -> tuple[float, float, float, float, float]:
    """Return identity, Fourier, pairing, and gauge errors."""

    left_perp = orthogonal_vector(left)
    right_perp = orthogonal_vector(right)
    scalar = np.asarray(
        [left.conj() @ coefficient @ right for coefficient in coefficients]
    )
    column_cross = np.asarray(
        [
            left_perp.conj() @ coefficient @ right
            for coefficient in coefficients
        ]
    )
    row_cross = np.asarray(
        [
            left.conj() @ coefficient @ right_perp
            for coefficient in coefficients
        ]
    )
    deficit = float(max(0.0, 1 - np.sum(np.abs(scalar) ** 2)))
    column_mass = float(np.sum(np.abs(column_cross) ** 2))
    row_mass = float(np.sum(np.abs(row_cross) ** 2))
    column_error = abs(column_mass - deficit)
    row_error = abs(row_mass - deficit)

    maximum_fourier_ratio = 0.0
    maximum_pairing_ratio = 0.0
    for shift in range(1, 11):
        autocorrelation = np.vdot(scalar[shift:], scalar[:-shift])
        cross_pairing = np.vdot(
            row_cross[shift:],
            column_cross[:-shift],
        )
        if deficit > 1e-12:
            maximum_fourier_ratio = max(
                maximum_fourier_ratio,
                abs(autocorrelation) / deficit,
            )
            maximum_pairing_ratio = max(
                maximum_pairing_ratio,
                abs(cross_pairing) / deficit,
            )
        elif max(abs(autocorrelation), abs(cross_pairing)) > 3e-11:
            raise RuntimeError("zero-defect Fourier response did not vanish")

    input_gauge = np.column_stack([right, -right_perp])
    output_gauge = np.column_stack([left, -left_perp])
    original_input = np.column_stack([right, right_perp])
    original_output = np.column_stack([left, left_perp])
    gauge_error = 0.0
    for coefficient in coefficients:
        block = (
            original_output.conj().T
            @ coefficient
            @ original_input
        )
        gauged = (
            output_gauge.conj().T
            @ coefficient
            @ input_gauge
        )
        gauge_error = max(
            gauge_error,
            float(
                abs(gauged[0, 0] - block[0, 0])
                + abs(gauged[0, 1] + block[0, 1])
                + abs(gauged[1, 0] + block[1, 0])
                + abs(gauged[1, 1] - block[1, 1])
            ),
        )
    return (
        column_error,
        row_error,
        maximum_fourier_ratio,
        maximum_pairing_ratio,
        gauge_error,
    )


def make_record(
    kind: str,
    step: float,
    count: int,
) -> ChannelResponseRecord:
    """Audit one exact matrix-inner Schur path."""

    coefficients = schur_path(kind, step, count)
    maximum_column_error = 0.0
    maximum_row_error = 0.0
    maximum_fourier_ratio = 0.0
    maximum_pairing_ratio = 0.0
    maximum_gauge_error = 0.0
    for left, right in channel_states():
        (
            column_error,
            row_error,
            fourier_ratio,
            pairing_ratio,
            gauge_error,
        ) = fixed_channel_statistics(
            coefficients,
            left,
            right,
        )
        maximum_column_error = max(
            maximum_column_error,
            column_error,
        )
        maximum_row_error = max(maximum_row_error, row_error)
        maximum_fourier_ratio = max(
            maximum_fourier_ratio,
            fourier_ratio,
        )
        maximum_pairing_ratio = max(
            maximum_pairing_ratio,
            pairing_ratio,
        )
        maximum_gauge_error = max(
            maximum_gauge_error,
            gauge_error,
        )

    split_error = 0.0
    if kind == "split_diagonal":
        for index in range(2):
            vector = np.eye(2, dtype=complex)[:, index]
            (
                column_error,
                row_error,
                _,
                _,
                _,
            ) = fixed_channel_statistics(
                coefficients,
                vector,
                vector,
            )
            scalar_mass = sum(
                abs(vector.conj() @ coefficient @ vector) ** 2
                for coefficient in coefficients
            )
            split_error = max(
                split_error,
                column_error,
                row_error,
                abs(1 - scalar_mass),
            )

    parseval = parseval_error(coefficients)
    verified = bool(
        parseval < 3e-11
        and maximum_column_error < 3e-11
        and maximum_row_error < 3e-11
        and maximum_fourier_ratio <= 1 + 3e-10
        and maximum_pairing_ratio <= 1 + 3e-10
        and maximum_gauge_error < 3e-12
        and split_error < 3e-11
    )
    if not verified:
        raise RuntimeError(f"{kind} channel-response audit failed")
    return ChannelResponseRecord(
        construction_kind=kind,
        step=format_float(step),
        state_count=len(channel_states()),
        maximum_parseval_error=format_float(parseval),
        maximum_column_identity_error=format_float(
            maximum_column_error
        ),
        maximum_row_identity_error=format_float(maximum_row_error),
        maximum_fourier_ratio=format_float(maximum_fourier_ratio),
        maximum_cross_pairing_ratio=format_float(
            maximum_pairing_ratio
        ),
        maximum_gauge_error=format_float(maximum_gauge_error),
        maximum_split_channel_error=format_float(split_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[ChannelResponseRecord]:
    """Return all deterministic channel-response audits."""

    return [
        make_record(kind, step, 100)
        for kind in (
            "transverse_pauli",
            "generic_nonnormal",
            "split_diagonal",
        )
        for step in (0.17, 0.08, 0.035)
    ]


def write_records(
    records: list[ChannelResponseRecord],
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
            "repeated_crabb_scalar_channel_response_factor_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the channel-response factorization audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
