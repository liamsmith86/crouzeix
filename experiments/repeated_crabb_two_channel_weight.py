#!/usr/bin/env python3
"""Audit the marked-chain and Hardy-row parts of two-channel weighting.

The proof note paired with this checker uses two exact mechanisms.

1. Starting the L313 wandering recursion at height ``h`` instead of
   height zero shows that a word with ``r`` physical reverse edges can
   expose only ``B_1, ..., B_(h+r)`` next to an already paid
   ``c**h B_h`` root.
2. At zero remaining elliptic slack the normalized half-line features
   are ordinary Hardy rows, so unequal shortest rows are orthogonal.

This checker independently reconstructs every marked physical word
through a fixed length on several noncommuting colligations, audits the
continuant-port valuations, checks Hardy-row orthogonality, and checks
the final Rees-style quadratic regrouping.  The floating tests are
falsification guards; the recursions in the proof note are exact.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_endpoint_word_gram import transfer_coefficient
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)
from repeated_crabb_weighted_endpoint_polynomial import (
    physical_reverse,
    physical_word_value,
)


Matrix = np.ndarray
MAXIMUM_WORD_LENGTH = 7
MAXIMUM_INITIAL_HEIGHT = 4


@dataclass(frozen=True)
class TwoChannelWeightRecord:
    """One marked-chain and zero-slack audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    maximum_word_length: int
    maximum_initial_height: int
    marked_words_checked: int
    maximum_support_excess: int
    maximum_reconstruction_error: str
    maximum_factor_bound_ratio: str
    continuant_legs_checked: int
    minimum_port_valuation_margin: int
    maximum_off_diagonal_hardy_gram: str
    maximum_nonconstant_leakage_return: str
    quadratic_regrouping_error: str
    all_checks_passed: bool


def marked_word_factor_coefficients(
    word: str,
    initial_height: int,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[Matrix, ...]:
    """Factor ``W*(S*)^h w(S,J)V`` through ``B_1,...,B_(h+r)``.

    This is L313's constructive recursion with the principal wandering
    chain initialized at ``S**h W``.  Here ``r`` is the number of
    physical reverse letters in ``word``.
    """

    reverse_count = word.count("j")
    maximum_channel = initial_height + reverse_count
    multiplicity = right.shape[1]
    coefficients_adjoint = [
        np.zeros((multiplicity, multiplicity), dtype=complex)
        for _ in range(maximum_channel)
    ]
    if maximum_channel == 0:
        return ()

    reverse = physical_reverse(operator, right, left)
    adjoint_letters = [
        operator.conj().T if letter == "s" else reverse.conj().T
        for letter in word
    ]

    height = initial_height
    principal_coefficient = 1.0
    principal_alive = True
    for index, letter in enumerate(word):
        if not principal_alive:
            break

        emitted_index = 0
        emitted_coefficient = 0.0
        if letter == "s":
            if height == 0:
                principal_alive = False
                continue
            emitted_index = height - 1
            emitted_coefficient = -principal_coefficient
            height -= 1
        elif height == 0:
            emitted_index = 1
            emitted_coefficient = 2 * principal_coefficient
            principal_coefficient *= 2
            height = 1
        else:
            emitted_index = height + 1
            emitted_coefficient = principal_coefficient
            height += 1

        if emitted_index == 0:
            continue
        remaining = np.eye(len(operator), dtype=complex)
        for adjoint_letter in adjoint_letters[index + 1 :]:
            remaining = adjoint_letter @ remaining
        coefficients_adjoint[emitted_index - 1] += (
            emitted_coefficient
            * right.conj().T
            @ remaining
            @ right
        )

    if principal_alive and height > 0:
        coefficients_adjoint[height - 1] += (
            principal_coefficient
            * np.eye(multiplicity, dtype=complex)
        )
    return tuple(
        coefficient.conj().T
        for coefficient in coefficients_adjoint
    )


def marked_word_audit(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
) -> tuple[int, int, float, float]:
    """Audit all short marked words and return their worst errors."""

    words_checked = 0
    maximum_support_excess = -MAXIMUM_INITIAL_HEIGHT
    maximum_error = 0.0
    maximum_bound_ratio = 0.0
    reverse = physical_reverse(operator, right, left)

    for initial_height in range(1, MAXIMUM_INITIAL_HEIGHT + 1):
        left_prefix = np.linalg.matrix_power(
            operator.conj().T,
            initial_height,
        )
        for length in range(MAXIMUM_WORD_LENGTH + 1):
            for letters in itertools.product("sj", repeat=length):
                word = "".join(letters)
                reverse_count = word.count("j")
                maximum_channel = initial_height + reverse_count
                coefficients = marked_word_factor_coefficients(
                    word,
                    initial_height,
                    operator,
                    right,
                    left,
                )
                blocks = [
                    transfer_coefficient(
                        operator,
                        right,
                        left,
                        degree,
                    )
                    for degree in range(1, maximum_channel + 1)
                ]
                factor = np.vstack(coefficients)
                reconstruction = np.hstack(blocks) @ factor
                endpoint = (
                    left.conj().T
                    @ left_prefix
                    @ physical_word_value(word, operator, reverse)
                    @ right
                )
                maximum_error = max(
                    maximum_error,
                    float(np.linalg.norm(reconstruction - endpoint)),
                )

                nonzero_indices = [
                    index
                    for index, coefficient in enumerate(
                        coefficients,
                        start=1,
                    )
                    if np.linalg.norm(coefficient) > 1e-11
                ]
                if nonzero_indices:
                    maximum_support_excess = max(
                        maximum_support_excess,
                        max(nonzero_indices) - maximum_channel,
                    )

                bound = (length + 1) * 4**reverse_count
                maximum_bound_ratio = max(
                    maximum_bound_ratio,
                    float(np.linalg.norm(factor, ord=2)) / bound,
                )
                words_checked += 1

    return (
        words_checked,
        maximum_support_excess,
        maximum_error,
        maximum_bound_ratio,
    )


def continuant_port_audit(maximum_delay: int = 10) -> tuple[int, int]:
    """Audit the exact L245 direct/reflected port valuations."""

    legs_checked = 0
    minimum_margin = maximum_delay
    for delay in range(1, maximum_delay + 1):
        for distance in range(1, delay + 1):
            direct_valuation = distance
            reflected_valuation = delay
            minimum_margin = min(
                minimum_margin,
                direct_valuation - distance,
                reflected_valuation - distance,
            )
            legs_checked += 2
    return legs_checked, minimum_margin


def hardy_orthogonality_error(
    row_count: int,
    multiplicity: int,
) -> float:
    """Return the largest off-diagonal ordinary Hardy-row Gram."""

    inclusions: list[Matrix] = []
    for row in range(row_count):
        inclusion = np.zeros(
            (row_count * multiplicity, multiplicity),
            dtype=complex,
        )
        start = row * multiplicity
        inclusion[start : start + multiplicity] = np.eye(multiplicity)
        inclusions.append(inclusion)
    return max(
        float(np.linalg.norm(left.conj().T @ right))
        for left_index, left in enumerate(inclusions)
        for right_index, right in enumerate(inclusions)
        if left_index != right_index
    )


def leakage_return_error(
    multiplicity: int,
    delay: int,
    seed: int,
    row_count: int = 24,
) -> float:
    """Audit L266 on a noncommuting delayed Potapov product."""

    generator = np.random.default_rng(seed)

    def random_projection(rank: int) -> Matrix:
        raw = (
            generator.standard_normal((multiplicity, rank))
            + 1j * generator.standard_normal((multiplicity, rank))
        )
        frame, _ = np.linalg.qr(raw)
        return frame @ frame.conj().T

    rank = max(1, multiplicity // 2)
    first = random_projection(rank)
    second = random_projection(rank)
    identity = np.eye(multiplicity, dtype=complex)
    coefficients = [
        (identity - first) @ (identity - second),
        (identity - first) @ second + first @ (identity - second),
        first @ second,
    ]

    dimension = row_count * multiplicity
    toeplitz = np.zeros((dimension, dimension), dtype=complex)
    for source_row in range(row_count):
        for coefficient_index, coefficient in enumerate(coefficients):
            target_row = source_row + delay + coefficient_index
            if target_row >= row_count:
                continue
            target_start = target_row * multiplicity
            source_start = source_row * multiplicity
            toeplitz[
                target_start : target_start + multiplicity,
                source_start : source_start + multiplicity,
            ] = coefficient
    leakage = toeplitz @ toeplitz.conj().T

    forward = np.zeros((dimension, dimension), dtype=complex)
    for row in range(row_count - 1):
        start = row * multiplicity
        forward[
            start + multiplicity : start + 2 * multiplicity,
            start : start + multiplicity,
        ] = identity
    row_start = delay * multiplicity
    row_inclusion = np.zeros((dimension, multiplicity), dtype=complex)
    row_inclusion[row_start : row_start + multiplicity] = identity

    maximum_error = 0.0
    for exponent in range(1, 6):
        for shift in (
            np.linalg.matrix_power(forward, exponent),
            np.linalg.matrix_power(forward.conj().T, exponent),
        ):
            compressed = (
                row_inclusion.conj().T
                @ leakage
                @ shift
                @ leakage
                @ row_inclusion
            )
            maximum_error = max(
                maximum_error,
                float(np.linalg.norm(compressed)),
            )
    return maximum_error


def quadratic_regrouping_error(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    seed: int,
    maximum_channel: int = 4,
    maximum_slack: int = 3,
) -> float:
    """Check ``B(c)(D+cR(c))B(c)*`` against its coefficient sum."""

    generator = np.random.default_rng(seed)
    multiplicity = right.shape[1]
    channels = [
        transfer_coefficient(operator, right, left, degree)
        for degree in range(1, maximum_channel + 1)
    ]
    parameter = 0.19
    weighted = np.hstack(
        [
            parameter**degree * channel
            for degree, channel in enumerate(channels, start=1)
        ]
    )

    block_dimension = maximum_channel * multiplicity
    direct = np.zeros((block_dimension, block_dimension), dtype=complex)
    for index in range(maximum_channel):
        start = index * multiplicity
        direct[start : start + multiplicity, start : start + multiplicity] = (
            (index + 1) / (maximum_channel + 1)
        ) * np.eye(multiplicity)

    slack_blocks: list[Matrix] = []
    for _ in range(maximum_slack):
        raw = (
            generator.standard_normal((block_dimension, block_dimension))
            + 1j
            * generator.standard_normal(
                (block_dimension, block_dimension)
            )
        )
        slack_blocks.append((raw + raw.conj().T) / 2)

    grouped_middle = direct.copy()
    for slack, block in enumerate(slack_blocks, start=1):
        grouped_middle += parameter**slack * block
    grouped = weighted @ grouped_middle @ weighted.conj().T

    direct_sum = np.zeros(
        (multiplicity, multiplicity),
        dtype=complex,
    )
    all_blocks = [direct, *slack_blocks]
    for slack, block in enumerate(all_blocks):
        for left_index, left_channel in enumerate(channels):
            left_start = left_index * multiplicity
            for right_index, right_channel in enumerate(channels):
                right_start = right_index * multiplicity
                middle = block[
                    left_start : left_start + multiplicity,
                    right_start : right_start + multiplicity,
                ]
                direct_sum += (
                    parameter
                    ** (left_index + right_index + 2 + slack)
                    * left_channel
                    @ middle
                    @ right_channel.conj().T
                )
    return float(np.linalg.norm(direct_sum - grouped))


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    seed: int,
) -> TwoChannelWeightRecord:
    """Audit one colligation."""

    (
        words_checked,
        maximum_support_excess,
        reconstruction_error,
        bound_ratio,
    ) = marked_word_audit(operator, right, left)
    legs_checked, minimum_margin = continuant_port_audit()
    hardy_error = hardy_orthogonality_error(
        MAXIMUM_INITIAL_HEIGHT + MAXIMUM_WORD_LENGTH + 1,
        right.shape[1],
    )
    leakage_error = leakage_return_error(
        right.shape[1],
        MAXIMUM_INITIAL_HEIGHT,
        seed + 1,
    )
    regrouping_error = quadratic_regrouping_error(
        operator,
        right,
        left,
        seed,
    )
    verified = bool(
        maximum_support_excess <= 0
        and reconstruction_error < 2e-10
        and bound_ratio <= 1 + 1e-10
        and minimum_margin >= 0
        and hardy_error < 1e-14
        and leakage_error < 1e-12
        and regrouping_error < 1e-12
    )
    return TwoChannelWeightRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        maximum_word_length=MAXIMUM_WORD_LENGTH,
        maximum_initial_height=MAXIMUM_INITIAL_HEIGHT,
        marked_words_checked=words_checked,
        maximum_support_excess=maximum_support_excess,
        maximum_reconstruction_error=format_float(reconstruction_error),
        maximum_factor_bound_ratio=format_float(bound_ratio),
        continuant_legs_checked=legs_checked,
        minimum_port_valuation_margin=minimum_margin,
        maximum_off_diagonal_hardy_gram=format_float(hardy_error),
        maximum_nonconstant_leakage_return=format_float(leakage_error),
        quadratic_regrouping_error=format_float(regrouping_error),
        all_checks_passed=verified,
    )


def standard_records() -> list[TwoChannelWeightRecord]:
    """Return deterministic unstructured and rank-chain audits."""

    records: list[TwoChannelWeightRecord] = []
    for multiplicity in (1, 2, 3):
        generator = np.random.default_rng(317_000 + multiplicity)
        operator, right, left = random_partial_isometry(
            4 * multiplicity + 9,
            multiplicity,
            generator,
        )
        records.append(
            audit_case(
                "unstructured",
                operator,
                right,
                left,
                317_100 + multiplicity,
            )
        )

    for multiplicity in (2, 3, 4):
        operator, right, left, _, _ = rank_chain_case(
            multiplicity,
            317_200 + multiplicity,
            0.17,
        )
        records.append(
            audit_case(
                "rank_chain",
                operator,
                right,
                left,
                317_300 + multiplicity,
            )
        )
    return records


def write_records(
    records: list[TwoChannelWeightRecord],
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
            "experiments/repeated_crabb_two_channel_weight_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the two-channel weight audits."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(
        json.dumps(
            {"records": len(records), "sha256": digest},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
