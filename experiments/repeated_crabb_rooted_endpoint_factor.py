#!/usr/bin/env python3
"""Audit the hereditary endpoint factor of every rooted bridge channel.

L305 sends a rooted bridge word to a transfer-channel value

    Phi(B_j* K + K* B_j)

modulo an ideal-valued response.  Polarized L280 converts that channel
value exactly to the literal endpoint factor

    K B_j* + B_j K*

modulo one explicit bounded response column.  This checker validates
the conversion, its norm bound, its flag support, and the important
fact that the converting column can remain nonzero on a complete
delay even though its endpoint response is then zero.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_all_grade_preimage import generalized_grade_column
from repeated_crabb_boundary_metric_flag import rank_chain_case
from repeated_crabb_canonical_quintic_preimage import joint_left_kernel
from repeated_crabb_endpoint_null_gauge import (
    delayed_random_partial_isometry,
)
from repeated_crabb_rooted_bridge_flux import (
    response_endpoint,
    transfer_channel,
)
from repeated_crabb_transfer_channel_covariance import (
    random_partial_isometry,
    transfer_coefficient,
)


Matrix = np.ndarray


@dataclass(frozen=True)
class RootedEndpointFactorRecord:
    """One rooted-channel-to-endpoint-factor audit."""

    construction_kind: str
    state_dimension: int
    defect_dimension: int
    delay: int
    transfer_norm: str
    multiplier_norm: str
    endpoint_factor_norm: str
    conversion_error: str
    perpendicular_error: str
    column_norm: str
    column_bound_slack: str
    transfer_flag_dimension: int
    factor_flag_compression_norm: str
    response_norm: str
    all_checks_passed: bool


def random_multiplier(
    multiplicity: int,
    generator: np.random.Generator,
) -> Matrix:
    """Return one normalized complex copy multiplier."""

    multiplier = (
        generator.standard_normal((multiplicity, multiplicity))
        + 1j
        * generator.standard_normal((multiplicity, multiplicity))
    )
    return multiplier / np.linalg.norm(multiplier)


def audit_case(
    construction_kind: str,
    operator: Matrix,
    right: Matrix,
    left: Matrix,
    delay: int,
    multiplier: Matrix,
    *,
    require_complete_delay: bool = False,
) -> RootedEndpointFactorRecord:
    """Audit one exact polarized conversion numerically."""

    transfer = transfer_coefficient(
        operator,
        right,
        left,
        delay,
    )
    root = (
        transfer.conj().T @ multiplier
        + multiplier.conj().T @ transfer
    )
    channel_value = transfer_channel(
        root,
        operator,
        right,
        left,
    )
    endpoint_factor = (
        transfer @ multiplier.conj().T
        + multiplier @ transfer.conj().T
    )
    column, _ = generalized_grade_column(
        operator,
        right,
        left,
        delay,
        multiplier,
        scale=1.0,
    )
    response = response_endpoint(
        column,
        operator,
        right,
        left,
    )
    conversion_error = float(
        np.linalg.norm(
            channel_value - (endpoint_factor - response)
        )
    )
    perpendicular_error = float(
        np.linalg.norm(right.conj().T @ column)
    )
    column_norm = float(np.linalg.norm(column))
    multiplier_norm = float(np.linalg.norm(multiplier))
    column_bound = delay * multiplier_norm
    flag = joint_left_kernel(transfer)
    flag_compression = (
        flag.conj().T @ endpoint_factor @ flag
        if flag.shape[1]
        else np.zeros((0, 0), dtype=complex)
    )
    transfer_norm = float(np.linalg.norm(transfer))
    response_norm = float(np.linalg.norm(response))
    tolerance = 2e-10
    verified = bool(
        conversion_error < tolerance
        and perpendicular_error < tolerance
        and column_norm <= column_bound + tolerance
        and np.linalg.norm(flag_compression) < tolerance
        and (
            not require_complete_delay
            or (
                transfer_norm < tolerance
                and response_norm < tolerance
                and column_norm > 1e-8
            )
        )
    )
    if not verified:
        raise RuntimeError(
            "rooted endpoint-factor audit failed: "
            f"{construction_kind=}, {delay=}, "
            f"{transfer_norm=}, {conversion_error=}, "
            f"{perpendicular_error=}, {column_norm=}, "
            f"{column_bound=}, {response_norm=}"
        )
    return RootedEndpointFactorRecord(
        construction_kind=construction_kind,
        state_dimension=len(operator),
        defect_dimension=right.shape[1],
        delay=delay,
        transfer_norm=format_float(transfer_norm),
        multiplier_norm=format_float(multiplier_norm),
        endpoint_factor_norm=format_float(
            float(np.linalg.norm(endpoint_factor))
        ),
        conversion_error=format_float(conversion_error),
        perpendicular_error=format_float(perpendicular_error),
        column_norm=format_float(column_norm),
        column_bound_slack=format_float(column_bound - column_norm),
        transfer_flag_dimension=flag.shape[1],
        factor_flag_compression_norm=format_float(
            float(np.linalg.norm(flag_compression))
        ),
        response_norm=format_float(response_norm),
        all_checks_passed=verified,
    )


def standard_records() -> list[RootedEndpointFactorRecord]:
    """Return deterministic generic, rank-defective, and delayed audits."""

    generator = np.random.default_rng(71310)
    records: list[RootedEndpointFactorRecord] = []
    for dimension, multiplicity in ((7, 2), (9, 3), (12, 4)):
        operator, right, left = random_partial_isometry(
            dimension,
            multiplicity,
            generator,
        )
        for delay in range(1, 5):
            records.append(
                audit_case(
                    "unstructured_partial_isometry",
                    operator,
                    right,
                    left,
                    delay,
                    random_multiplier(multiplicity, generator),
                )
            )

    for dimension, multiplicity in ((8, 2), (12, 3), (16, 4)):
        operator, right, left = delayed_random_partial_isometry(
            dimension,
            multiplicity,
            generator,
        )
        records.append(
            audit_case(
                "complete_first_delay",
                operator,
                right,
                left,
                1,
                random_multiplier(multiplicity, generator),
                require_complete_delay=True,
            )
        )
    for multiplicity in (3, 4, 5):
        operator, right, left, _, _ = rank_chain_case(
            multiplicity,
            81400 + multiplicity,
        )
        for delay in (1, 2):
            records.append(
                audit_case(
                    "rank_chain_partial_isometry",
                    operator,
                    right,
                    left,
                    delay,
                    random_multiplier(multiplicity, generator),
                )
            )
    return records


def write_records(
    records: list[RootedEndpointFactorRecord],
    output: Path,
) -> str:
    """Write deterministic records atomically and return their hash."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        "".join(
            json.dumps(asdict(record), sort_keys=True) + "\n"
            for record in records
        ),
        encoding="utf-8",
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
            "repeated_crabb_rooted_endpoint_factor_s70226.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run and persist every audit."""

    args = parse_args()
    records = standard_records()
    digest = write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))
    print(json.dumps({"records": len(records), "sha256": digest}))


if __name__ == "__main__":
    main()
