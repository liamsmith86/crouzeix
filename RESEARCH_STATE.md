# Research summary

Updated 2026-09-06 for the repository's public-facing research record.

## Original campaign

The July work explored scalar Crouzeix bounds, stronger complete bounds,
structured matrix families and local perturbations of equality cases.
The [README](README.md) introduces the main results and reproduction commands.
The full historical statement and approach records are in [archive/](archive/README.md).

The final local calculation, L355, expresses the metric/zero residual gap as
one weighted square. Its physical shape factorization and common quadratic
term remain open in that approach. The September
[local audit](proof/audit_l352_l355_20260905.md) identifies the optimization
scope correction and the separate local-to-global steps.

## General-proof reconstruction

The September work reconstructed Lorist–Schwenninger's all-powers argument
for every finite complex matrix and scalar polynomial:

- [Finite-matrix proof](proof/CROUZEIX_PROOF.md)
- [Independent reconstruction and audit](proof/audit_global_proof_20260905.md)
- [Source versions and related proofs](proof/literature_refresh_20260905.md)
- [Exact algebra check](experiments/global_power_lemma_exact_audit_20260905.py)
- [Additional contour and hypothesis checks](experiments/all_powers_dilation_audit.py)

The proof audit records a successful verification of that finite-dimensional
argument. The general proof and its decisive lemma are attributed to the
source authors; the earlier local work is preserved separately.
