#!/usr/bin/env python3
"""Audit complete circular-normal absorption at a repeated Crabb block.

The scalar L188 estimate cannot be amplified one-sidedly: noncommuting copy
matrices distinguish ``X*X`` from ``XX*``.  The two reflected circular
characters instead produce the Jordan-symmetric gain

    8 / C * (P*P + PP*).

This script checks that this gain is dominated by the operator-valued Hardy
residual Gram, verifies the strict improvement from L173's null lift, and
records a sharp nonnormal example that falsifies the tempting one-sided
amplification.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from crabb_block_hardy_equality import format_float


@dataclass(frozen=True)
class JordanAbsorptionRecord:
    """One random reflection-symmetric block-residual audit."""

    record_kind: str
    length: int
    multiplicity: int
    active_mode_count: int
    reflection_error: str
    flux_residual_minimum_eigenvalue: str
    actual_residual_minimum_eigenvalue: str
    minimum_null_lift_ratio: str
    all_checks_passed: bool


@dataclass(frozen=True)
class OneSidedObstructionRecord:
    """One sharp nonnormal counterexample to one-sided amplification."""

    record_kind: str
    length: int
    multiplicity: int
    mode: int
    symmetric_flux_residual_norm: str
    one_sided_residual_minimum_eigenvalue: str
    actual_residual_minimum_eigenvalue: str
    one_sided_amplification_fails: bool
    all_checks_passed: bool


AuditRecord = JordanAbsorptionRecord | OneSidedObstructionRecord


def block_view(
    matrix: np.ndarray,
    block_count: int,
    multiplicity: int,
) -> np.ndarray:
    """Return a ``block_count x block_count`` array of copy blocks."""

    return matrix.reshape(
        block_count,
        multiplicity,
        block_count,
        multiplicity,
    ).transpose(0, 2, 1, 3)


def reflected_residual(
    length: int,
    multiplicity: int,
    seed: int,
) -> np.ndarray:
    """Return ``F`` with ``(J tensor I) F`` Hermitian."""

    block_count = length - 1
    rng = np.random.default_rng(seed)
    raw = rng.standard_normal(
        (block_count * multiplicity,) * 2
    ) + 1j * rng.standard_normal((block_count * multiplicity,) * 2)
    hermitian = (raw + raw.conj().T) / 2
    reversal = np.kron(
        np.fliplr(np.eye(block_count)),
        np.eye(multiplicity),
    )
    return reversal @ hermitian


def mode_constants(length: int, mode: int) -> tuple[float, float, float]:
    """Return ``alpha``, flux curvature, and actual L173 curvature."""

    remainder = length - mode
    choose_three = remainder * (remainder - 1) * (remainder - 2) / 6
    alpha = 4 * (4 * mode - 1) / length**2
    flux = alpha**2 * choose_three / 64
    null_lift = (
        2
        * mode
        * (mode - 1)
        * (mode - 2)
        * (remainder + 1 / 4) ** 2
        / (3 * length**4)
    )
    return alpha, flux, flux + null_lift


def residual_gram(blocks: np.ndarray) -> np.ndarray:
    """Return the copy-space right Gram of all residual blocks."""

    multiplicity = blocks.shape[-1]
    gram = np.zeros((multiplicity, multiplicity), dtype=complex)
    for row in range(blocks.shape[0]):
        for column in range(blocks.shape[1]):
            block = blocks[row, column]
            gram += block.conj().T @ block
    return gram


def response_row(
    blocks: np.ndarray,
    length: int,
    mode: int,
) -> np.ndarray:
    """Return L188's lower anti-diagonal block response ``P``."""

    multiplicity = blocks.shape[-1]
    response = np.zeros((multiplicity, multiplicity), dtype=complex)
    remainder = length - mode
    for left in range(1, length):
        right = remainder - left
        if left < right < length:
            skew = (
                blocks[left - 1, right - 1]
                - blocks[right - 1, left - 1]
            ) / 2
            response += (right - left) * skew
    return response


def completed_gains(
    blocks: np.ndarray,
    length: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return flux-only and actual Jordan-symmetric completed gains."""

    multiplicity = blocks.shape[-1]
    flux_gain = np.zeros((multiplicity, multiplicity), dtype=complex)
    actual_gain = np.zeros_like(flux_gain)
    ratios = []
    for mode in range(3, length - 2):
        remainder = length - mode
        choose_three = (
            remainder * (remainder - 1) * (remainder - 2) / 6
        )
        response = response_row(blocks, length, mode)
        symmetric_square = (
            response.conj().T @ response
            + response @ response.conj().T
        )
        alpha, flux, actual = mode_constants(length, mode)
        flux_gain += 8 * symmetric_square / choose_three
        actual_gain += alpha**2 * symmetric_square / (8 * actual)
        ratios.append(actual / flux)
    return flux_gain, actual_gain, min(ratios, default=float("inf"))


def random_record(
    length: int,
    multiplicity: int,
    seed: int,
) -> JordanAbsorptionRecord:
    """Check both complete absorption inequalities on random block data."""

    residual = reflected_residual(length, multiplicity, seed)
    block_count = length - 1
    blocks = block_view(residual, block_count, multiplicity)
    reversal = np.kron(
        np.fliplr(np.eye(block_count)),
        np.eye(multiplicity),
    )
    reflection_error = float(
        np.linalg.norm(reversal @ residual - (reversal @ residual).conj().T)
    )
    base_loss = 4 * residual_gram(blocks)
    flux_gain, actual_gain, minimum_ratio = completed_gains(blocks, length)
    flux_residual = (base_loss - flux_gain + (base_loss - flux_gain).conj().T) / 2
    actual_residual = (
        base_loss - actual_gain + (base_loss - actual_gain).conj().T
    ) / 2
    flux_minimum = float(np.linalg.eigvalsh(flux_residual)[0])
    actual_minimum = float(np.linalg.eigvalsh(actual_residual)[0])
    verified = bool(
        reflection_error < 1e-11
        and flux_minimum > -1e-10
        and actual_minimum > -1e-10
        and minimum_ratio > 1
    )
    if not verified:
        raise RuntimeError("the random Jordan-absorption audit failed")
    return JordanAbsorptionRecord(
        record_kind="random_jordan_absorption",
        length=length,
        multiplicity=multiplicity,
        active_mode_count=max(0, length - 5),
        reflection_error=format_float(reflection_error),
        flux_residual_minimum_eigenvalue=format_float(flux_minimum),
        actual_residual_minimum_eigenvalue=format_float(actual_minimum),
        minimum_null_lift_ratio=format_float(minimum_ratio),
        all_checks_passed=verified,
    )


def sharp_obstruction_record(length: int, multiplicity: int) -> OneSidedObstructionRecord:
    """Build a flux-sharp example on which the one-sided lift is false."""

    if length < 6 or multiplicity < 2:
        raise ValueError("the sharp obstruction needs length >= 6 and multiplicity >= 2")
    mode = 3
    remainder = length - mode
    block_count = length - 1
    blocks = np.zeros(
        (block_count, block_count, multiplicity, multiplicity),
        dtype=complex,
    )
    nilpotent = np.zeros((multiplicity, multiplicity), dtype=complex)
    nilpotent[0, 1] = 1
    for left in range(1, length):
        right = remainder - left
        if not (left < right < length):
            continue
        weight = right - left
        first = weight * nilpotent
        second = -first
        blocks[left - 1, right - 1] = first
        blocks[right - 1, left - 1] = second
        reflected_left = length - right
        reflected_right = length - left
        blocks[reflected_left - 1, reflected_right - 1] = first.conj().T
        blocks[reflected_right - 1, reflected_left - 1] = second.conj().T

    response = response_row(blocks, length, mode)
    choose_three = remainder * (remainder - 1) * (remainder - 2) / 6
    base_loss = 4 * residual_gram(blocks)
    symmetric_flux_gain = 8 * (
        response.conj().T @ response
        + response @ response.conj().T
    ) / choose_three
    one_sided_gain = 16 * response.conj().T @ response / choose_three
    alpha, _, actual = mode_constants(length, mode)
    actual_gain = alpha**2 * (
        response.conj().T @ response
        + response @ response.conj().T
    ) / (8 * actual)
    symmetric_residual = base_loss - symmetric_flux_gain
    one_sided_residual = (
        base_loss - one_sided_gain
        + (base_loss - one_sided_gain).conj().T
    ) / 2
    actual_residual = (
        base_loss - actual_gain
        + (base_loss - actual_gain).conj().T
    ) / 2
    symmetric_norm = float(np.linalg.norm(symmetric_residual))
    one_sided_minimum = float(np.linalg.eigvalsh(one_sided_residual)[0])
    actual_minimum = float(np.linalg.eigvalsh(actual_residual)[0])
    verified = bool(
        symmetric_norm < 1e-10
        and one_sided_minimum < -1e-6
        and actual_minimum > -1e-10
    )
    if not verified:
        raise RuntimeError("the one-sided obstruction audit failed")
    return OneSidedObstructionRecord(
        record_kind="sharp_one_sided_obstruction",
        length=length,
        multiplicity=multiplicity,
        mode=mode,
        symmetric_flux_residual_norm=format_float(symmetric_norm),
        one_sided_residual_minimum_eigenvalue=format_float(
            one_sided_minimum
        ),
        actual_residual_minimum_eigenvalue=format_float(actual_minimum),
        one_sided_amplification_fails=True,
        all_checks_passed=verified,
    )


def write_records(path: Path, records: Sequence[AuditRecord]) -> None:
    """Write deterministic JSON Lines output atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True))
            handle.write("\n")
    temporary.replace(path)


def main() -> None:
    """Run random and sharp complete-amplification audits."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--minimum-length", type=int, default=6)
    parser.add_argument("--maximum-length", type=int, default=12)
    parser.add_argument(
        "--multiplicities",
        type=int,
        nargs="+",
        default=(2, 3, 4),
    )
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_circular_jordan_absorption_s70224.jsonl"
        ),
    )
    args = parser.parse_args()

    records: list[AuditRecord] = []
    for multiplicity in args.multiplicities:
        for length in range(args.minimum_length, args.maximum_length + 1):
            record = random_record(
                length,
                multiplicity,
                args.seed + 100 * multiplicity + length,
            )
            records.append(record)
            print(
                "verified Jordan absorption "
                f"length {length}, multiplicity {multiplicity}",
                flush=True,
            )
        obstruction = sharp_obstruction_record(
            args.maximum_length,
            multiplicity,
        )
        records.append(obstruction)
        print(
            "verified one-sided obstruction "
            f"length {args.maximum_length}, multiplicity {multiplicity}",
            flush=True,
        )
    write_records(args.output, records)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
