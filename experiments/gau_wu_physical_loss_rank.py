#!/usr/bin/env python3
"""Probe the low-rank physical loss after the exact Gau--Wu reductions."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from gau_wu_finite_hessian_jet import support_jet
from gau_wu_second_support_gram import boundary_angular_derivative
from gau_wu_two_sided_endpoint_defect import build_endpoint_form_audit


@dataclass(frozen=True)
class PhysicalLossRankRecord:
    """One complete reduced physical-loss audit."""

    dimension: int
    sample: int
    seed: int
    angle_count: int
    physical_dimension: int
    expected_loss_rank: int
    observed_loss_rank: int
    expected_nullity: int
    observed_nullity: int
    minimum_physical_remainder_eigenvalue: float
    minimum_loss_eigenvalue: float
    maximum_loss_eigenvalue: float
    smallest_active_loss_singular_value: float
    largest_null_loss_singular_value: float
    all_checks_passed: bool


def audit_model(
    dimension: int,
    sample: int,
    seed: int,
    angle_count: int,
) -> PhysicalLossRankRecord:
    """Audit the reduced physical remainder and its loss rank."""

    audit = build_endpoint_form_audit(
        dimension,
        sample,
        seed,
        angle_count,
    )
    physical_dimension = len(audit.directions)
    zero_dimension = 2 * (dimension - 1)
    two_sided = audit.two_sided_form
    physical_block = two_sided[:physical_dimension, :physical_dimension]
    mixed_block = two_sided[:physical_dimension, physical_dimension:]
    zero_block = two_sided[physical_dimension:, physical_dimension:]
    physical_remainder = physical_block - mixed_block @ np.linalg.solve(
        zero_block,
        mixed_block.T,
    )

    _, second_support, _ = support_jet(
        audit.matrix,
        list(audit.directions),
        angle_count,
    )
    angles = np.linspace(0, 2 * np.pi, angle_count, endpoint=False)
    weight = boundary_angular_derivative(audit.zeros, angles)
    support_gram = 2 * np.mean(
        weight[None, None, :] * second_support,
        axis=2,
    )

    loss = (
        2 * np.eye(physical_dimension)
        + 2 * support_gram
        - physical_remainder
    )
    loss = (loss + loss.T) / 2
    loss_eigenvalues = np.linalg.eigvalsh(loss)
    loss_singular_values = np.linalg.svd(loss, compute_uv=False)
    expected_rank = 6 * dimension - 14
    expected_nullity = (dimension - 4) ** 2
    rank_threshold = 1e-7 * loss_singular_values[0]
    observed_rank = int(np.sum(loss_singular_values > rank_threshold))
    observed_nullity = physical_dimension - observed_rank
    smallest_active = loss_singular_values[observed_rank - 1]
    largest_null = (
        loss_singular_values[observed_rank]
        if observed_rank < physical_dimension
        else 0.0
    )
    minimum_remainder = float(
        np.linalg.eigvalsh(physical_remainder)[0]
    )

    checks = (
        zero_dimension == two_sided.shape[0] - physical_dimension
        and physical_dimension == (dimension - 1) ** 2 + 1
        and expected_rank <= physical_dimension
        and loss_eigenvalues[0] > -3e-9
        and minimum_remainder > -3e-9
        and observed_rank == expected_rank
        and observed_nullity == expected_nullity
        and largest_null < 2e-6 * smallest_active
    )
    if not checks:
        raise RuntimeError(
            "physical loss-rank audit failed: "
            f"n={dimension}, remainder={minimum_remainder}, "
            f"loss={loss_eigenvalues[0]}, "
            f"rank={observed_rank}/{expected_rank}, "
            f"null={observed_nullity}/{expected_nullity}, "
            f"tail={largest_null}/{smallest_active}"
        )
    return PhysicalLossRankRecord(
        dimension=dimension,
        sample=sample,
        seed=seed,
        angle_count=angle_count,
        physical_dimension=physical_dimension,
        expected_loss_rank=expected_rank,
        observed_loss_rank=observed_rank,
        expected_nullity=expected_nullity,
        observed_nullity=observed_nullity,
        minimum_physical_remainder_eigenvalue=minimum_remainder,
        minimum_loss_eigenvalue=float(loss_eigenvalues[0]),
        maximum_loss_eigenvalue=float(loss_eigenvalues[-1]),
        smallest_active_loss_singular_value=float(smallest_active),
        largest_null_loss_singular_value=float(largest_null),
        all_checks_passed=True,
    )


def write_records(
    records: list[PhysicalLossRankRecord],
    output: Path,
) -> str:
    """Write deterministic JSON Lines atomically and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)
    return hashlib.sha256(output.read_bytes()).hexdigest()


def parse_dimensions(value: str) -> tuple[int, ...]:
    """Parse a comma-separated dimension list."""

    dimensions = tuple(int(item) for item in value.split(","))
    if not dimensions or any(dimension < 4 for dimension in dimensions):
        raise argparse.ArgumentTypeError("dimensions must be at least four")
    return dimensions


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dimensions",
        type=parse_dimensions,
        default=(4, 5, 6, 7, 8, 9),
    )
    parser.add_argument("--samples", type=int, default=2)
    parser.add_argument("--angle-count", type=int, default=512)
    parser.add_argument("--seed", type=int, default=70224)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/gau_wu_physical_loss_rank_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the requested physical loss-rank probes."""

    arguments = parse_args()
    records: list[PhysicalLossRankRecord] = []
    for dimension in arguments.dimensions:
        for sample in range(arguments.samples):
            seed = arguments.seed + 1009 * dimension + sample
            record = audit_model(
                dimension,
                sample,
                seed,
                arguments.angle_count,
            )
            records.append(record)
            print(
                json.dumps(
                    {
                        "dimension": dimension,
                        "sample": sample,
                        "minimum_loss_eigenvalue": (
                            record.minimum_loss_eigenvalue
                        ),
                        "observed_loss_rank": record.observed_loss_rank,
                        "observed_nullity": record.observed_nullity,
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
    digest = write_records(records, arguments.output)
    print(
        json.dumps(
            {
                "all_checks_passed": all(
                    record.all_checks_passed for record in records
                ),
                "record_count": len(records),
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
