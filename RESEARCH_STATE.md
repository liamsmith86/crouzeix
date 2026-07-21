# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-20 (Epoch 1, late)

## Objective
Resolve Crouzeix's conjecture (R(A,p) ≤ 2 universally): general proof, certified counterexample, or forced-stop checkpoint.

## Verified literature baseline (2026-07-20)
- OPEN as of July 2026. Universal constant: 1+√2 (C–P 2017). Per-shape: 1+√(1+a(Ω)) with a(Ω)<1 (MMOR IMRN 2025). Per-dimension: C_N < 1+√2 (MMOR 2024, non-constructive).
- a(Ω) → 1 for thin quadrilaterals ⟹ domain-only constants can NEVER reach 2. Operator-coupled structure required.
- R–S 2018: abstract Lemma 1.1 sharp at 1+√2; their example has non-unital α and Ω ⊉ W(T). Question 4.1 (unitality ⟹ 2?) — status unknown, likely open (BGG+20 poses it as open; check).
- Full dissection of C–P/S–dV machinery: proof/track_A_crouzeix_palencia.md.

## THE CURRENT FRONTIER (original work, this campaign)
proof/refined_master_inequality.md. At extremal (f₀,x₀,μ), K>1, with W=‖γ_Φ(f₀)x₀‖, β=⟨γ_Φ(f₀)x₀,x₀⟩, c=|∫Φ(f₀)f₀dμ|:
- L6 (PROVED): K² ≤ K·√(W²−|β|²) + c.  [Strictly sharper than S–dV Thm 5: vector-localized, free scalar subtraction via ⟨γ(f₀)x₀,x₀⟩=0.]
- **L7/L7′ (UNRESOLVED, central target): √(W²−|β|²) ≤ 2 − c/2 ⟹ CROUZEIX.**
  - Sharp at Crabb equality cases (2 ≤ 2). Numerically HOLDS at all near-extremal configs tested (slack +0.11…+0.76, smallest at most-extremal). Adversarial searches running (n=2,3,4).
- L5 (PROVED): W² ≤ 2∮|f₀|²dν, ν := ⟨P(·)x₀,x₀⟩ds (mass 2). W=2 ⟹ |f₀|=1 ν-a.e.
- Saturation rigidity (derived, to write up): W=2 forces P^{1/2}(σ)y = e^{−iψ}f₀(σ)P^{1/2}(σ)x₀ (y = output direction). First contradiction attempts suggest W=2 & c>0 incompatible — quantify into W ≤ 2 − φ(c).
- L8 (heuristic): direct sums (M₁⊕normal junk) kill any global-operator-norm route (e.g. ∮λ_min(P)ds → 0) but NOT the localized quantities. Explains why MMOR is dimension/shape-limited. The localized frame is the right invariance class.

## Track status
- **A (dissect & improve C–P):** primary; produced L3–L8. Next: prove L7 near saturation via rigidity; understand ν vs μ (ν = 2Re η, η = resolvent boundary measure; μ = positive rep. of analytic part).
- **C (counterexample):** 5 search runs done (dims 3–5): max R_inner ≈ 1.9747, no violation; consistent with G–O. Coverage gap: pipeline rejects nonsmooth/flat W(M) (certificates); need smoothed-hull Ω to search those. adversarial_L7 runs: br053wnqg (n2), btutuitrr (n3), bunfhcuya (n4).
- **B, D, E, F:** queued. E partially absorbed into A (double-layer geometry).
- **G (audit):** active — caught: (i) fake violations from degenerate boundaries (fixed via self-certifying quadrature: mass=2I, positivity of P, tangent winding, Cauchy reconstruction); (ii) stale-process contamination of logs (fixed: harness-tracked background jobs). RULE: no claim from a run without process hygiene + certificates.

## Exact next actions
1. Await adversarial L7' results; analyze near-tight configurations (which (M,p,x₀) minimize slack?).
2. Write the rigidity argument (W→2) rigorously; attempt quantitative version W ≤ 2 − φ(c) — candidate route via L5 + phase alignment + mass constraint ∮dν=2 + |∫Φf₀dμ| ≤ 1.
3. Determine status of R–S Question 4.1 (if open, L7′ is strictly stronger and more concrete — good).
4. Read Crouzeix 2004 & 2007 equality analysis + BGG+20 survey (extremal pair structure: is |f₀|=1 on supp μ? is f₀ Blaschke-like of degree ≤ n−1? — needed for saturation analysis).
5. Extend pipeline to nonsmooth W via smoothed support function; rerun searches there (flat portions are conjecturally hard cases).
