# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-20 (Epoch 4 open)

## Objective & status
Resolve Crouzeix's conjecture. OPEN in literature (verified 2026-07). No counterexample found by
campaign (all candidate violations were certified numerical artifacts). Campaign has reduced the
conjecture to a single sharply-supported positivity conjecture and proved it for n = 2:

**H-r: ρ := Re⟨(f₀·Φ(f₀))(A)x₀,x₀⟩ ≥ 0 at extremal pairs for Ω = int W(A).**
H-r ⟹ Crouzeix (via SV24 Thm 6.1 K ≤ 1+√(1−ρ), + shrinking; proof/strategy_S.md).
Posed-but-unattacked in literature (SV24 §6 remark); no refutation exists (searched).

## Proved by campaign (see proof/)
- **2×2 theorem**: ρ = 1 − π/(2K(m)) ≥ 0 in closed form at critical (confocal-ellipse) domains;
  equality iff disk. Geometric form φ'(0) ≥ φ(1). [Modulo α=0 symmetry step — numerically certain,
  write-up pending.] Also re-derives K ≤ 2 for 2×2 via K = √k(h+√(1+h²)) ≤ 2. (rho_2x2_theorem.md)
- P1: K² + ρ ≤ Kq, q = √(W²−|β|²) (phase-preserving; refines SV24's K²+ρ≤2K by β-subtraction).
- Ceiling: K² + 2ρ + G² ≤ 4 ⟹ at K=2: ρ ≤ −G²/2; with H-r ⟹ Crabb rigidity at K=2.
- Hadamard variation of Φ (single-layer form); Λ-density formula for frozen dρ; ∮Λds = 0 at disks
  (verified numerically to 6 digits).
- DLP contact degeneracy at Ω = W(A): ker P(σ_θ) = span((σ_θ−A)x_θ) (derivation done, write-up pending).
- Scalar reduction L10 DISPROVED (odd-symmetric mechanism) — operator structure essential.

## Key structural picture (Epoch 3–4)
- Criticality = confocality: Kippenhahn curve of ∂W(A) has foci at σ(A); the 2×2 mechanism
  (confocal ellipse ⟹ ρ ≥ 0 via elliptic integrals) is the n=2 case of a general
  "Pick-problem at eigenvalue images + confocal domain" structure.
- Extremal problem = n-point Pick-boundary problem at wⱼ = φ(σ(A)); extremal zeros ≈ critical
  points of ∏b_{wⱼ} (EXACT at n=2 — hyperbolic midpoint; close at n=3, deviation = non-Hermitian
  weighting). (zero_geometry.py)
- ρ-forms: ρ = K·Re⟨g₀(A)u₀,x₀⟩ = 1 − ReΣrᵢ⟨hᵢ(A)x₀,x₀⟩ ("capacity partition", terms ≈ harmonic
  split; Σ = 1 exactly on disk).
- Empirics at critical domains: ρ ≥ 0 always; ρ ↑ as Ω ↓ W(A) (offsets, Minkowski; textbook-monotone);
  adversarial min-ρ floors: +2e-4 (n=3), +1.7e-2 (n=4), +1.6e-2 (n=5, weak search); n=6 pending.
  ρ = 0 attained at disk-like configs. Hadamard: pair-response dominates frozen term (heavy route).

## Current next actions (Epoch 4)
1. n=6 adversarial result (bg job bj2wx10ia); log_Hr_n6_s42.txt.
2. Semi-closed ρ for the doubly-symmetric 3×3 family (analog of A² = I collapse: use
   Cayley–Hamilton A³ = c₁A for zero-diag tridiagonal + symmetry) — hunt the n=3 elliptic-type
   formula and its positivity mechanism. THE main theory push.
3. Rigor debts: (a) α=0 (2×2) symmetry/uniqueness proof; (b) contact-degeneracy lemma write-up;
   (c) mpmath independent recheck of one near-tight n=3 extremal; (d) improve general-domain
   pipeline accuracy (interp error ~1e-4 found by exact-disk cross-check — spline/higher-N).
4. Read: de Vries PhD thesis (Twente 2025, "C-spectral sets and related estimates") — ρ follow-ups;
   SV25 PAMS (ρ-contractions).
5. Pick-geometry formulation: express ρ entirely in terms of (wⱼ, Pick data, Jordan weights);
   the conjecture becomes a finite-dim inequality per n — target a proof for n = 3.
6. Keep committing after each task (user instruction).

## Files map
proof/: track_A_crouzeix_palencia.md, refined_master_inequality.md, epoch2_extremal_structure.md,
P2_target.md, strategy_S.md, rho_2x2_theorem.md, rho_positivity_program.md (MASTER).
experiments/: crouzeix.py (basics), extremal_pullback.py (ψ-domain exact machinery),
theodorsen.py (+GeneralPullback: arbitrary convex domains), ellipse_sandbox.py (2×2 exact),
minkowski_test.py (best_extremal solver + Minkowski families), Hr_test.py / Hr_adversarial.py,
sym3_structure.py, zero_geometry.py, L12_test.py, sanity.py; searches: search.py, adversarial_L7.py,
scalar stuff (L10_test, scalar_adv, odd_symmetric_test — historical).
Ledgers: LEMMA_LEDGER.md, APPROACH_LEDGER.md (pitfalls P1–P4!), LITERATURE_LEDGER.md,
COUNTEREXAMPLE_SEARCH.md. Audit: chatgpt/FABLE_RESEARCH_AUDIT.md (reconciled 2026-07-20).
