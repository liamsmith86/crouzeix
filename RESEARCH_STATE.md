# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-20 (end of Epoch 2, mid-campaign)

## Objective
Resolve Crouzeix's conjecture (‖p(A)‖ ≤ 2 max_{W(A)}|p|): proof, certified counterexample, or forced-stop checkpoint. OPEN in literature as of 2026-07 (verified); best universal constant 1+√2.

## CampaignSTATUS one-liner
No counterexample found (all "violations" were quadrature artifacts — see APPROACH_LEDGER pitfalls P1–P4).
Main line: refined master inequality (L6, proved) + target inequality at extremal pairs. Current concrete target:
**T: c ≤ 2(2−K) and X-excess control at true extremal pairs ⟹ CONJECTURE.**

## The framework (all in proof/refined_master_inequality.md, proof/epoch2_extremal_structure.md)
Setting: Ω smooth convex ⊇ W(A); γ = calculus; Φ = conjugate Cauchy transform; extremal pair
(f₀,x₀), K = ‖f₀(A)‖ = ‖γ‖ > 1; f₀ = B∘φ Blaschke (Crouzeix 04); u₀ = f₀(A)x₀/K;
g₀ = Φ(f₀) = 1/f₀ − PP (principal parts at zeros of f₀) [D4]; C = ⟨(g₀f₀)(A)x₀,x₀⟩ = ∫g₀f₀dμ,
c = |C|; G = ‖g₀(A)*x₀‖; W = ‖(f₀(A)+g₀(A)*)x₀‖; β = ⟨(f₀(A)+g₀(A)*)x₀,x₀⟩ = conj(∫g₀dμ).
- L6 (PROVED): K² ≤ K√(W²−|β|²) + c. Identity: W² = K² + 2ReC + G².
- L7@ext (target): √(W²−|β|²) ≤ 2−c/2 at extremals ⟹ K ≤ 2. Equivalent form:
  K² + X ≤ (2−c/2)², X := 2ReC + (G²−|β|²).

## Epoch-2 empirical laws at TRUE extremal pairs (exact Blaschke machinery, domains ψ_ε = w+εw²)
(experiments/extremal_pullback.py, epoch2_neardisk.py, L12_test.py; diag certificate ~1e-8)
- E-N1: c = 0 exactly on disk (theorem); c > 0 generically off-disk but TINY: c ~ 0.006 at ε=0.25
  with 0.3-perturbed Crabb. Roughly c = O(ε²)–O(ε³) along families.
- E-N2: G ≈ |β| to ~4 digits (x₀ ≈ eigenvector of g₀(A)* with eigenvalue β̄); orthogonal defect
  same order as c. So X ≈ 2ReC ≤ 0 mostly (ReC predominantly negative off-disk).
- E-N3: L7@ext HOLDS in all probes, slack ≥ 0.01 (Crabb near-fill) to 0.46; c/(2(2−K)) ≤ ~0.01.
- E-N4: L12-strong (f₀dμ ⊥ A(Ω)) FALSE (m=3: ⟨(f₀w)(A)x₀,x₀⟩ ~ 3e-3 ≠ 0). c-smallness is subtler.
- Danger zone K→2 empirically forces ε→0 (disk-like), where c→0 faster. Two-regime structure.

## Disproved / dead (do NOT revisit without new idea)
- L10 scalar reduction — DISPROVED via odd-symmetric construction (ellipse, ±x point masses, odd
  near-unimodular p). Operator compression-loss is essential.
- (⋆) c ≤ 2∮λ_min(P)ds — wrong shape (not ⊕-stable).
- L12-strong annihilation — false (above).
- Relaxing extremality to only ⟨p(A)x₀,x₀⟩=0: L7′ unverified violations only from bad numerics;
  but L10-counterexample warns the relaxed version may be false; not load-bearing either way.

## EPOCH 3 (current) — P2 target after external-audit reconciliation
Audit chatgpt/FABLE_RESEARCH_AUDIT.md adopted after independent verification (see proof/P2_target.md):
- P1 (PROVED): K² + ReC ≤ Kq, q := √(W²−|β|²). [SV24 arXiv:2409.15954 §6 has K²+ρ≤2K, ρ=ReC —
  verified against PDF; our new content = β-subtraction q < 2.]
- **P2 (CENTRAL TARGET): q ≤ 2 + ReC/2 at extremal pairs ⟹ CROUZEIX.** Hard case ReC < 0 (generic).
- Stinespring defect form: 4−q² = 4[‖(I−Q)Uξ‖² + |⟨Uξ,ξ⟩|²]; P2 ⟺ leakage+overlap ≥ −r/2 − r²/16.
Next actions:
1. [numerics] Retarget adversarial search to minimize s_phase := 2 + ReC/2 − q at TRUE extremal
   pairs (ψ-domains, n=2,3; push r < 0, q → 2). experiments/ratio_adversarial.py → p2_adversarial.
2. [numerics] Exact 2×2 sandbox with Ω = W(A) (ellipse; elliptic-function Riemann map; deg-1
   Blaschke): attempt P2 falsification over eccentricity/matrix params. STOP/GO: P2 false at exact
   2×2 extremal ⟹ abandon scalar-tradeoff targets, return to angle form of P1.
3. [theory] If P2 survives: Euler–Lagrange at extremals (Blaschke zero-variations, singular-vector
   stationarity, D1 via conjugate-algebra argument) → Stinespring stability bound.
4. [rigor] mpmath independent reimplementation for near-tight cases; error bars.
5. [coverage] Near-polygon ψ-families (flat-ish sides) once P2 tooling is in place.
6. [ops] Git repo active: commit after each task (user instruction 2026-07-20).

## Key files
- proof/track_A_crouzeix_palencia.md — C–P/S–dV dissection (where 1+√2 enters).
- proof/refined_master_inequality.md — L3–L9.
- proof/epoch2_extremal_structure.md — extremal-pair theory E1–E5, D1–D5, two-regime plan.
- experiments/crouzeix.py — ratio basics; sanity suite (Crabb R=2 ✓).
- experiments/extremal_pullback.py — EXACT extremal machinery on ψ(D) domains (main tool).
- experiments/epoch2_neardisk.py, L12_test.py — Epoch-2 probes (data in this file).
- experiments/adversarial_L7.py + star_inequality.py + refined_test.py — operator L7′ tools (certified).
- experiments/L10_test.py, scalar_adv.py, odd_symmetric_test.py — scalar track (L10 disproof).
- COUNTEREXAMPLE_SEARCH.md, APPROACH_LEDGER.md (pitfalls P1–P4 IMPORTANT), LEMMA_LEDGER.md, LITERATURE_LEDGER.md.
