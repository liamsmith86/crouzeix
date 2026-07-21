# The ρ-positivity program (H-r) — consolidated theory (end of Epoch 3)

**Central conjecture (H-r).** At every extremal pair for (A, Ω) with Ω = int W(A) (σ(A) ⊂ int W(A),
smooth strictly convex case first): ρ := Re⟨(f₀·Φ(f₀))(A)x₀, x₀⟩ ≥ 0.
**H-r ⟹ Crouzeix** via SV24 Thm 6.1 + shrinking (see strategy_S.md). Status of ρ ≥ 0 in the
literature: posed by SV24 (§6, "interestingly, if..."), no proof or refutation found (searched 2026-07-20).

## Equivalent forms (all verified algebraically)
1. ρ = Re∫f₀g₀ dμ = Re∫f₀g₀ dμ̃ (μ, μ̃ = right/left extremal measures; D2-exact).
2. ρ = K·Re⟨g₀(A)u₀, x₀⟩ — "conjugate transform maps left singular vector toward right with ≥ 0
   correlation" (u₀ ⊥ x₀).
3. ρ = 1 − Re Σᵢ rᵢ⟨hᵢ(A)x₀,x₀⟩, hᵢ = f₀/(·−zᵢ), rᵢ = 1/f₀'(zᵢ) (zᵢ = zeros of f₀ = B∘φ).
   "Capacity partition": terms Tᵢ ≈ harmonic split of 1 (numerics: Tᵢ ≈ 1/2 each for symmetric deg-2).
4. Stinespring: ρ-floor ⟺ leakage/overlap bound (P2-defect) — see P2_target.md.

## Bounds already proved
- Ceiling: K² + 2ρ + G² ≤ 4 (from W ≤ 2 and W² = K² + 2ρ + G², G = ‖g₀(A)*x₀‖).
  At K = 2: ρ ≤ −G²/2 ≤ 0; with H-r ⟹ ρ = 0, G = 0: full Crabb rigidity at K = 2.
- |ρ| ≤ a(Ω) < 1 (SV24/MMOR); ρ trivially ∈ [min, max of Re(f₀g₀) on supp μ].
- P1: K² + ρ ≤ Kq (q = √(W²−|β|²)).

## The 2×2 theorem (proved modulo α=0-symmetry step; rho_2x2_theorem.md)
ρ = 1 − π/(2K(m)) ≥ 0, equality iff disk. Geometric form: φ'(0) ≥ φ(1), i.e.
crad(Ω,0)·tanh(d_hyp,Ω(0, eigenvalue)/2) ≤ 1 for CONFOCAL ellipse Ω (foci = eigenvalues).
**Criticality = confocality for 2×2.** Non-confocal Ω ⊇ σ(A) can violate (observed ρ < 0 for
egg-shaped Ω ⊋ W(A)). The n ≥ 3 proof must encode "Ω = W(A)" the way confocality does for n = 2.

## Where criticality must enter (all identities so far are criticality-blind)
- DLP contact degeneracy at Ω = W(A): for each θ, ker P(σ_θ) ∋ k_θ := (σ_θ − A)x_θ, where x_θ =
  maximizing vector of Re(e^{-iθ}A) [derivation: ⟨P(σ)x,x⟩ = ⟨H_θ Rx, Rx⟩, H_θ = Re(e^{-iθ}(σ−A)),
  R = (σ−A)^{-1}; kernel = span((σ_θ−A)x_θ) at contact — TO WRITE UP RIGOROUSLY].
- Hadamard variation of Φ at fixed f (derived via Cauchy–Pompeiu):
  δΦ(f)(z) = (1/π)∮_{∂Ω} conj(f'(σ))(σ−z)^{-1} δh(σ) ds — single-layer with density conj(f')δh.
  Full dρ/dt needs coupled extremal-pair variation (δf₀, δx₀) — NOT killed by envelope (ρ ≠ K).
  Next: derive δK, δf₀ stationarity system; test dρ/dt sign semi-analytically.

## Empirical laws (certified numerics; see RESEARCH_STATE for data)
- ρ ≥ 0 at critical domains: exact 2×2 (all), n=3 shrink-families (3/3, ρ monotone ↑ as Ω ↓ W(A)),
  adversarial min-ρ floors: +2e-4 (n=3), +1.7e-2 (n=4); n=5,6 running.
- ρ = 0 attained: disks (proved); ρ → 0 also along near-shift/near-disk families.
- Minkowski path W(A)→disk: ρ monotone ↓ to 0 (crabb case textbook; random case needs
  pipeline accuracy fix — flagged).

## Epoch-4 additions (derived 2026-07-20, after sym3 reduction)
- Annihilation counting: exact g₀(σ(A)) = 0 is 2n real conditions vs 2(n−1) zero-freedoms — always
  2 short. K-maximization drives toward the annihilation manifold; ρ is LINEAR in the residuals
  (g₀(λⱼ)-values) with coefficients from the extremal frame ⟹ positivity must come from a
  KKT/multiplier alignment: **conjectured proof shape: ρ = (multiplier ≥ 0)·(slack) via the
  constrained-optimality of the extremal Blaschke.** NEXT: write the full stationarity system
  (∂K/∂αᵢ = 0) and check whether it forces the sign of Σ coeffⱼ·g₀(λⱼ).
- Third measure ρ̂ ≥ 0 (mixed): ⟨q(A)x₀,u₀⟩ = K∫q·conj(f₀)dρ̂ (from additive-perturbation
  stationarity; f̄₀ = 1/f₀ on ∂Ω). Boundary density J := φ·(log B)'∘φ = Σᵢ(1−|αᵢ|²)/|φ−αᵢ|² > 0
  (angular derivative; |f₀'| = J·|φ'| on ∂Ω).
- Löwner structure: under domain inflation, δφ = −φ·H (H Herglotz); frozen-B derivative of K² is
  −2K·Re⟨(B'(φ)φH)(A)x₀,u₀⟩ = −2K²∫J·ReH dρ̂ ≤ 0 automatically (ReH|∂ = delta-type ≥ 0) —
  K-monotonicity constraints carry NO extra info (dead end, recorded to avoid repetition).

## Open determinations / next attacks (Epoch 4 candidates)
1. Rigorous write-up: contact-degeneracy lemma; α=0 for 2×2; the K ≤ 2 (2×2) elliptic re-derivation.
2. Hadamard system: δ(f₀,x₀,K,ρ) under support variation; sign of dρ/dt — first do it numerically
   (finite differences along offset families — cheap, already validated), then symbolically.
3. n=3 confocal analog: which domain family plays "confocal" for a 3×3? Guess: level sets of
   |char-poly|^{2/n}-capacity or the "Blaschke-critical" domains where zeros zᵢ align with... probe:
   compute zᵢ(A) at critical domains and look for invariant characterization (numerics).
4. de Vries dissertation "C-spectral sets and related estimates" (Twente 2025) — read for ρ-related
   follow-ups; also SV25 "sharp bound for functional calculus of ρ-contractions" (PAMS 2025).
5. Formalize the K=2 rigidity endpoint: H-r + ceiling ⟹ any K=2 extremal is Crabb-type; then
   a stability/openness argument could exclude K > 2 without full H-r (alternative route).
