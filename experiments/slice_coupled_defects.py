#!/usr/bin/env python3
"""Audit the Stein-defect parameterization of elliptic-slice metrics.

This is a floating-point KKT regression for
``proof/slice_coupled_defects.md``.  The parameterization itself is exact;
the SDP ranks printed here are diagnostics, not proof certificates.
"""

from __future__ import annotations

import numpy as np

from slice_cb_sdp import DEFAULT_CASES, solve_similarity_sdp
from slice_similarity_duality import ModalSlice, modal_slice_from_weights


RANK_TOLERANCE = 2e-6


def modal_metric_and_defects(
    data: ModalSlice, metric: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    diagonal = np.diag([1.0, data.c])
    odd_change = diagonal @ data.left_rotation
    even_change = diagonal @ data.right_rotation
    odd_metric = odd_change.T @ metric[:2, :2] @ odd_change
    even_metric = even_change.T @ metric[2:, 2:] @ even_change
    sigma = np.diag(data.nodes)
    odd_defect = odd_metric - data.c * sigma @ even_metric @ sigma
    even_defect = even_metric - sigma @ odd_metric @ sigma / data.c
    return odd_metric, even_metric, odd_defect, even_defect


def reconstruct_from_defects(
    data: ModalSlice, odd_defect: np.ndarray, even_defect: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    node_products = np.outer(data.nodes, data.nodes)
    kernel = 1.0 / (1.0 - node_products**2)
    sigma = np.diag(data.nodes)
    odd_metric = kernel * (
        odd_defect + data.c * sigma @ even_defect @ sigma
    )
    even_metric = kernel * (
        even_defect + sigma @ odd_defect @ sigma / data.c
    )
    return odd_metric, even_metric


def numerical_rank(matrix: np.ndarray) -> int:
    eigenvalues = np.linalg.eigvalsh((matrix + matrix.T) / 2)
    scale = max(1.0, float(np.max(abs(eigenvalues))))
    return int(np.count_nonzero(eigenvalues > RANK_TOLERANCE * scale))


def run_case(weights: tuple[float, float, float, float]) -> None:
    data = modal_slice_from_weights(weights)
    result = solve_similarity_sdp(data.operator)
    odd_metric, even_metric, odd_defect, even_defect = modal_metric_and_defects(
        data, result.metric
    )
    rebuilt_odd, rebuilt_even = reconstruct_from_defects(
        data, odd_defect, even_defect
    )
    reconstruction_error = max(
        float(np.max(abs(rebuilt_odd - odd_metric))),
        float(np.max(abs(rebuilt_even - even_metric))),
    )

    witness = result.dual_witness
    dual_ranks = (numerical_rank(witness[:2, :2]), numerical_rank(witness[2:, 2:]))
    defect_ranks = (numerical_rank(odd_defect), numerical_rank(even_defect))

    if reconstruction_error > 2e-7:
        raise AssertionError("the defect formula failed to reconstruct the modal metric")
    if min(np.linalg.eigvalsh(odd_defect)) < -2e-6:
        raise AssertionError("the odd contraction defect is not positive semidefinite")
    if min(np.linalg.eigvalsh(even_defect)) < -2e-6:
        raise AssertionError("the even contraction defect is not positive semidefinite")
    for dual_rank, defect_rank in zip(dual_ranks, defect_ranks, strict=True):
        if dual_rank + defect_rank > 2:
            raise AssertionError("complementary rank inequality failed numerically")

    print(
        f"weights={weights} t={result.bound:.9f} "
        f"dual_ranks={dual_ranks} defect_ranks={defect_ranks} "
        f"reconstruction_error={reconstruction_error:.2e}"
    )


def main() -> None:
    for weights in DEFAULT_CASES:
        run_case(weights)


if __name__ == "__main__":
    main()
