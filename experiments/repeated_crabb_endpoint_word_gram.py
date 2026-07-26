#!/usr/bin/env python3
"""Audit sharp and constructive Gram bounds for endpoint words.

For a balanced pure partial-isometry colligation and a word ``w`` of
length ``d`` in ``S,S*``, put ``K_w=W* w V``.  The lossless-circuit
conjecture predicts

    K_w K_w* <= sum_{j=1}^d B_j B_j*,

where ``B_j=W*(S*)^j V``.  The proved normal-form theorem gives the
same factorization with norm at most ``ceil(d/2)``.  The checker
exhausts all binary words through the requested horizon on
unstructured, rank-changing, and completely delayed colligations.
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
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_rooted_bridge_flux import word_value
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
)


Matrix = np.ndarray
MAXIMUM_HORIZON = 10


@dataclass(frozen=True)
class EndpointWordGramRecord:
    """One exhaustive finite-horizon endpoint-word audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    horizon: int
    words_checked: int
    maximum_forward_energy_identity_error: str
    minimum_gram_slack_eigenvalue: str
    maximum_contracting_factor_norm: str
    maximum_factor_reconstruction_error: str
    maximum_normal_form_factor_norm: str
    maximum_normal_form_bound_ratio: str
    maximum_normal_form_reconstruction_error: str
    maximum_delayed_short_word_norm: str
    all_checks_passed: bool


def transfer_coefficient(
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    delay: int,
) -> Matrix:
    """Return ``B_delay``."""

    return (
        left.conj().T
        @ np.linalg.matrix_power(operator.conj().T, delay)
        @ right
    )


def binary_words(length: int) -> tuple[str, ...]:
    """Return every word of the requested length deterministically."""

    return tuple(
        "".join(letters)
        for letters in itertools.product("sa", repeat=length)
    )


def endpoint_normal_form_coefficients(
    word: str,
    operator: Matrix,
    right: Matrix,
) -> tuple[Matrix, ...]:
    """Return the constructive prefix coefficients of ``W*word*V``.

    If the current word is ``a^p s v``, use

        (S*)^p S = (S*)^(p-1) - (S*)^(p-1) V V*

    at the left endpoint.  Every pass deletes two letters, and every
    emitted right multiplier is a contraction ``V* v V``.
    """

    multiplicity = right.shape[1]
    coefficients = [
        np.zeros((multiplicity, multiplicity), dtype=complex)
        for _ in word
    ]
    current = word
    while current:
        if current[0] == "s":
            break
        initial_adjoint_run = len(current) - len(
            current.lstrip("a")
        )
        if initial_adjoint_run == len(current):
            coefficients[initial_adjoint_run - 1] += np.eye(
                multiplicity
            )
            break

        suffix = current[initial_adjoint_run + 1 :]
        transfer_index = initial_adjoint_run - 1
        if transfer_index >= 1:
            suffix_value = (
                word_value(suffix, operator)
                if suffix
                else np.eye(len(operator), dtype=complex)
            )
            coefficients[transfer_index - 1] -= (
                right.conj().T @ suffix_value @ right
            )
        current = "a" * transfer_index + suffix

    return tuple(coefficients)


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    horizon: int,
    *,
    complete_delay: int = 0,
) -> EndpointWordGramRecord:
    """Audit every word through one horizon."""

    transfer_blocks: list[Matrix] = []
    gram = np.zeros(
        (left.shape[1], left.shape[1]),
        dtype=complex,
    )
    minimum_slack = np.inf
    maximum_factor_norm = 0.0
    maximum_reconstruction_error = 0.0
    maximum_normal_factor_norm = 0.0
    maximum_normal_bound_ratio = 0.0
    maximum_normal_reconstruction_error = 0.0
    maximum_delayed_word_norm = 0.0
    maximum_energy_identity_error = 0.0
    words_checked = 0

    for length in range(1, horizon + 1):
        transfer = transfer_coefficient(
            operator,
            right,
            left,
            length,
        )
        transfer_blocks.append(transfer)
        gram += transfer @ transfer.conj().T
        forward_residual = (
            np.linalg.matrix_power(operator, length + 1) @ left
        )
        forward_energy_gram = (
            np.eye(left.shape[1])
            - forward_residual.conj().T @ forward_residual
        )
        maximum_energy_identity_error = max(
            maximum_energy_identity_error,
            float(np.linalg.norm(gram - forward_energy_gram)),
        )
        transfer_row = np.hstack(transfer_blocks)
        transfer_row_pseudoinverse = np.linalg.pinv(
            transfer_row,
            rcond=1e-11,
        )

        for word in binary_words(length):
            endpoint_word = (
                left.conj().T
                @ word_value(word, operator)
                @ right
            )
            word_gram = endpoint_word @ endpoint_word.conj().T
            slack = gram - word_gram
            slack = (slack + slack.conj().T) / 2
            minimum_slack = min(
                minimum_slack,
                float(np.linalg.eigvalsh(slack)[0]),
            )

            factor = transfer_row_pseudoinverse @ endpoint_word
            reconstruction_error = float(
                np.linalg.norm(transfer_row @ factor - endpoint_word)
            )
            maximum_reconstruction_error = max(
                maximum_reconstruction_error,
                reconstruction_error,
            )
            maximum_factor_norm = max(
                maximum_factor_norm,
                float(np.linalg.norm(factor, ord=2)),
            )

            normal_coefficients = endpoint_normal_form_coefficients(
                word,
                operator,
                right,
            )
            normal_factor = np.vstack(normal_coefficients)
            normal_reconstruction_error = float(
                np.linalg.norm(
                    transfer_row @ normal_factor - endpoint_word
                )
            )
            normal_factor_norm = float(
                np.linalg.norm(normal_factor, ord=2)
            )
            normal_bound = (length + 1) // 2
            maximum_normal_reconstruction_error = max(
                maximum_normal_reconstruction_error,
                normal_reconstruction_error,
            )
            maximum_normal_factor_norm = max(
                maximum_normal_factor_norm,
                normal_factor_norm,
            )
            maximum_normal_bound_ratio = max(
                maximum_normal_bound_ratio,
                normal_factor_norm / normal_bound,
            )

            if length < complete_delay:
                maximum_delayed_word_norm = max(
                    maximum_delayed_word_norm,
                    float(np.linalg.norm(endpoint_word)),
                )
            words_checked += 1

    tolerance = 2e-8
    verified = bool(
        minimum_slack > -tolerance
        and maximum_energy_identity_error < tolerance
        and maximum_factor_norm < 1 + tolerance
        and maximum_reconstruction_error < tolerance
        and maximum_normal_reconstruction_error < tolerance
        and maximum_normal_bound_ratio < 1 + tolerance
        and maximum_delayed_word_norm < tolerance
    )
    if not verified:
        raise RuntimeError(
            "endpoint-word Gram audit failed: "
            f"{construction_kind=}, {horizon=}, "
            f"{minimum_slack=}, {maximum_energy_identity_error=}, "
            f"{maximum_factor_norm=}, "
            f"{maximum_reconstruction_error=}, "
            f"{maximum_normal_factor_norm=}, "
            f"{maximum_normal_bound_ratio=}, "
            f"{maximum_normal_reconstruction_error=}, "
            f"{maximum_delayed_word_norm=}"
        )

    return EndpointWordGramRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        horizon=horizon,
        words_checked=words_checked,
        maximum_forward_energy_identity_error=format_float(
            maximum_energy_identity_error
        ),
        minimum_gram_slack_eigenvalue=format_float(minimum_slack),
        maximum_contracting_factor_norm=format_float(
            maximum_factor_norm
        ),
        maximum_factor_reconstruction_error=format_float(
            maximum_reconstruction_error
        ),
        maximum_normal_form_factor_norm=format_float(
            maximum_normal_factor_norm
        ),
        maximum_normal_form_bound_ratio=format_float(
            maximum_normal_bound_ratio
        ),
        maximum_normal_form_reconstruction_error=format_float(
            maximum_normal_reconstruction_error
        ),
        maximum_delayed_short_word_norm=format_float(
            maximum_delayed_word_norm
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[EndpointWordGramRecord]:
    """Return deterministic exhaustive word audits."""

    records = []
    for multiplicity in (1, 2, 3, 4):
        operator, right, left = random_partial_isometry(
            4 * multiplicity + 8,
            multiplicity,
            np.random.default_rng(146_000 + multiplicity),
        )
        records.append(
            audit_case(
                "unstructured",
                operator,
                right,
                left,
                MAXIMUM_HORIZON,
            )
        )

    for multiplicity in (2, 3, 4):
        operator, right, left, _, _ = rank_chain_case(
            multiplicity,
            146_100 + multiplicity,
            0.001,
        )
        records.append(
            audit_case(
                "rank_chain",
                operator,
                right,
                left,
                MAXIMUM_HORIZON,
            )
        )

    generator = np.random.default_rng(146_200)
    for multiplicity in (2, 3, 4):
        operator, right, left = delayed_random_partial_isometry(
            4 * multiplicity,
            multiplicity,
            generator,
        )
        records.append(
            audit_case(
                "complete_first_delay",
                operator,
                right,
                left,
                MAXIMUM_HORIZON,
                complete_delay=2,
            )
        )

    return records


def write_records(
    records: list[EndpointWordGramRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(
                json.dumps(asdict(record), sort_keys=True) + "\n"
            )
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
            "repeated_crabb_endpoint_word_gram_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist the endpoint-word Gram audits."""

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
