#!/usr/bin/env python3
"""Audit the ordered Schur-kernel flag of the repeated transfer.

The one-step matrix Schur identity splits the de Branges--Rovnyak
kernel into one rank-``m`` layer and the shifted kernel of the next
Schur iterate.  Iteration gives ``L`` ordered feature layers.  This
checker verifies that factorization on L201 transfers and at repeated
Crabb apices.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from crabb_block_hardy_equality import format_float
from repeated_crabb_inner_faber_transfer import (
    canonical_transfer_data,
    strengthened_inverse_toeplitz,
    transfer_coefficients,
)
from repeated_crabb_matrix_schur_chart import (
    defect_roots,
    recover_parameters,
    schur_value,
)
@dataclass(frozen=True)
class SchurKernelFlagRecord:
    """One ordered model-kernel factorization audit."""

    construction_kind: str
    length: int
    multiplicity: int
    toeplitz_strength: str
    kernel_factorization_error: str
    feature_spanning_singular_value: str
    apex_monomial_feature_error: str
    all_checks_passed: bool


def stage_values(
    parameters: list[np.ndarray],
    terminal: np.ndarray,
    value: complex,
) -> list[np.ndarray]:
    """Return every nested Schur iterate at one point."""

    stages = [
        np.zeros_like(terminal)
        for _ in range(len(parameters) + 1)
    ]
    stages[-1] = terminal
    identity = np.eye(len(terminal), dtype=complex)
    for index in reversed(range(len(parameters))):
        parameter = parameters[index]
        left, right = defect_roots(parameter)
        shifted = value * stages[index + 1]
        stages[index] = (
            parameter
            + left
            @ shifted
            @ np.linalg.inv(
                identity + parameter.conj().T @ shifted
            )
            @ right
        )
    return stages


def schur_features(
    parameters: list[np.ndarray],
    terminal: np.ndarray,
    value: complex,
) -> list[np.ndarray]:
    """Return the ordered feature factors of the model kernel."""

    stages = stage_values(parameters, terminal, value)
    identity = np.eye(len(terminal), dtype=complex)
    product = identity
    features: list[np.ndarray] = []
    for index, parameter in enumerate(parameters):
        _, right = defect_roots(parameter)
        factor = (
            np.linalg.inv(
                identity
                + parameter.conj().T
                @ (value * stages[index + 1])
            )
            @ right
        )
        product = factor @ product
        features.append(value**index * product)
    return features


def audit_transfer(
    construction_kind: str,
    length: int,
    multiplicity: int,
    inverse_toeplitz: np.ndarray,
    strength: float,
) -> SchurKernelFlagRecord:
    """Audit one L201 transfer and its ordered feature flag."""

    hermitian = np.linalg.inv(inverse_toeplitz)
    data = canonical_transfer_data(
        hermitian,
        length,
        multiplicity,
    )
    count = max(80, 12 * length)
    coefficients = list(transfer_coefficients(data, count))
    parameters, terminal, _ = recover_parameters(
        coefficients,
        length,
    )
    identity = np.eye(multiplicity, dtype=complex)

    maximum_kernel_error = 0.0
    for first_index in range(6):
        first = 0.73 * np.exp(
            2j * np.pi * (first_index + 0.17) / 6
        )
        for second_index in range(5):
            second = 0.69 * np.exp(
                2j * np.pi * (second_index + 0.31) / 5
            )
            first_value = schur_value(
                parameters,
                terminal,
                first,
            )
            second_value = schur_value(
                parameters,
                terminal,
                second,
            )
            direct = (
                identity - first_value.conj().T @ second_value
            ) / (1 - np.conjugate(first) * second)
            first_features = schur_features(
                parameters,
                terminal,
                first,
            )
            second_features = schur_features(
                parameters,
                terminal,
                second,
            )
            factored = sum(
                first_feature.conj().T @ second_feature
                for first_feature, second_feature in zip(
                    first_features,
                    second_features,
                    strict=True,
                )
            )
            maximum_kernel_error = max(
                maximum_kernel_error,
                float(np.linalg.norm(direct - factored)),
            )

    sampled_columns = []
    for sample in range(length):
        value = 0.61 * np.exp(
            2j * np.pi * (sample + 0.23) / length
        )
        sampled_columns.append(
            np.vstack(
                schur_features(parameters, terminal, value)
            )
        )
    feature_matrix = np.hstack(sampled_columns)
    smallest_singular = float(
        np.linalg.svd(feature_matrix, compute_uv=False)[-1]
    )

    apex_error = 0.0
    is_apex = construction_kind == "repeated_crabb_apex"
    if is_apex:
        for sample in range(5):
            value = 0.77 * np.exp(
                2j * np.pi * (sample + 0.29) / 5
            )
            features = schur_features(
                parameters,
                terminal,
                value,
            )
            for index, feature in enumerate(features):
                apex_error = max(
                    apex_error,
                    float(
                        np.linalg.norm(
                            feature - value**index * identity
                        )
                    ),
                )

    tolerance = 3e-8
    verified = bool(
        maximum_kernel_error < tolerance
        and smallest_singular > 1e-4
        and (not is_apex or apex_error < tolerance)
    )
    if not verified:
        raise RuntimeError(
            "Schur kernel flag audit failed: "
            f"kernel={maximum_kernel_error:.3e}, "
            f"span={smallest_singular:.3e}, "
            f"apex={apex_error:.3e}"
        )
    return SchurKernelFlagRecord(
        construction_kind=construction_kind,
        length=length,
        multiplicity=multiplicity,
        toeplitz_strength=format_float(strength),
        kernel_factorization_error=format_float(
            maximum_kernel_error
        ),
        feature_spanning_singular_value=format_float(
            smallest_singular
        ),
        apex_monomial_feature_error=(
            format_float(apex_error)
            if is_apex
            else "not_applicable"
        ),
        all_checks_passed=verified,
    )


def standard_records() -> list[SchurKernelFlagRecord]:
    """Return deterministic apex and noncommuting records."""

    records: list[SchurKernelFlagRecord] = []
    for length in range(2, 6):
        for multiplicity in (2, 3):
            inverse, strength = strengthened_inverse_toeplitz(
                length,
                multiplicity,
                0.22,
            )
            records.append(
                audit_transfer(
                    "noncommuting_inverse_block_toeplitz",
                    length,
                    multiplicity,
                    inverse,
                    strength,
                )
            )
            apex_inverse = 2 * np.eye(
                length * multiplicity,
                dtype=complex,
            )
            records.append(
                audit_transfer(
                    "repeated_crabb_apex",
                    length,
                    multiplicity,
                    apex_inverse,
                    0.0,
                )
            )
    return records


def write_records(
    records: list[SchurKernelFlagRecord],
    output: Path,
) -> None:
    """Write deterministic JSON Lines atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    temporary.replace(output)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "experiments/"
            "repeated_crabb_schur_kernel_flag_s70224.jsonl"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the complete ordered kernel audit."""

    args = parse_args()
    records = standard_records()
    write_records(records, args.output)
    for record in records:
        print(json.dumps(asdict(record), sort_keys=True))


if __name__ == "__main__":
    main()
