# Level-4 closed-form theory (Epoch 5, 2026-07-21) — nodal reduction, phases, even-phase theorems

Family: 4×4 zero-diag tridiagonal, slice b_j = c·a_j (W(A) = exact ellipse, foci ±e₁ = outer
eigenvalue pair, τ₁ = φ(e₁) = √k). General principle first; slice used for exactness.

## 1. Nodal reduction ⟹ the extremal problem is a 4-point Pick problem with explicit objective
A diagonalizable ⟹ f₀(A) depends only on the nodal values v = (f₀(e))_{e∈σ(A)}. Hence
K(f₀) = Ψ(v) with Ψ explicit (norm of the interpolant matrix), and
  K_max = max{ Ψ(v) : v admissible Pick data for nodes φ(σ(A)) ⊂ D }.
Pick: v admissible iff the 4×4 Pick matrix ⪰ 0; extremals are Blaschke of deg ≤ 3 = rank.
Z₂ symmetry (A ~ −A via D = diag(1,−1,1,−1)) ⟹ extremal f₀ has definite parity:
ODD sector f₀(−z) = −f₀(z) or EVEN sector f₀(−z) = f₀(z); each leaves 2 free values
(v₁,v₂) = (f₀(e₁), f₀(e₂)).

## 2. Master closed form for the objective (VERIFIED 1e-16, slice_invariants.py V1)
Active block (odd sector) M = A₊(F₁Q₁+F₂Q₂), F_j = v_j/e_j, Q_j = spectral projections of
N = A₋A₊. Then, with G = A₊ᵀA₊, p_ij = tr(Q_iᵀGQ_j), δ = det A₊:
  T = p₁₁F₁² + 2p₁₂F₁F₂ + p₂₂F₂²  (= ‖M‖_F²),   det M = δF₁F₂,
  **K² = (T + √(T² − 4δ²F₁²F₂²))/2.**
Even sector: X = v₁Q₁+v₂Q₂ (per diagonal block), same formula with G = I, δ = 1.
The frame enters through THREE scalars (p₁₁,p₁₂,p₂₂) + δ only. This explains the
frame-transcendence finding: no (α,τ,e)-only law exists (PSLQ negative was correct), but the
law is algebraic once p_ij, δ are adjoined.

## 3. Odd phase: explicit stationarity law (VERIFIED 1e-10 at optimizer α*)
B odd Blaschke zeros {0,±α}: ∂_αB(w) = −2αw(1−w⁴)/(1−α²w²)². dK²/dα = 0 ⟺
  **Σ_j (T_j + S_j)·τ_j(1−τ_j⁴) / (e_j(1−α²τ_j²)²) = 0**,
T_j = ∂T/∂F_j, S_j = (T·T_j − 4δ²F_jF_{3−j}²)/√(T²−4δ²F₁²F₂²).
Frame identities (VERIFIED 1e-16): c_j = ∂K²/∂F_j/(2K) = (T_j+S_j)/(4K) and Σ_jF_jc_j = K
(Hellmann–Feynman; c_j = ⟨A₊Q_jx₀,y₀⟩ the frame coefficients from the earlier session).

## 4. EVEN PHASE THEOREMS (new; the sym3 mechanism generalizes exactly)
Even f₀ = c·(φ²−α²)/(1−α²φ²) = c·b_v(φ(z)²), v = α² — a DEGREE-1 Möbius in u = φ(z)²
with nodes u_j = τ_j². On the slice u₁ = τ₁² = k exactly (focus law).
(i) **Midpoint law and globality** [PROVED]: K²(B₁,B₂) is symmetric (κ := ‖Q₁‖_F² = ‖Q₂‖_F² for
    complementary idempotents ⟹ T = κ(B₁²+B₂²) + 2(1−κ)B₁B₂) and even under joint negation;
    the Möbius flow is ∂_v b_v(u) = −(1−b_v(u)²)/(1−v²). Hence at B₂ = −B₁ stationarity holds
    identically: ∂₁K²(B,−B)+∂₂K²(−B,B)-cancellation by symmetry+evenness. So the even extremal
    satisfies **b_v(u₂) = −b_v(u₁) ⟺ v = pseudo-hyperbolic midpoint of (u₁,u₂)** — the
    2×2/sym3 midpoint law, one Landen level up. The missing global step is now proved in
    `proof/even_pick_globality.md`: for the complete two-node Schur/Pick problem,
    `sup ||F(u₁)Q₁+F(u₂)Q₂|| = max(1, d||Q₁−Q₂||)`, and in the norm>1 regime the midpoint
    automorphism is the unique nonconstant global maximizer up to phase.
(ii) **q = 1/2 law** [PROVED, sympy]: at B₂ = −B₁, X = B₁R with R = Q₁−Q₂ a traceless
    involution; for the generic 2×2 traceless involution the top right-singular vector obeys
    x₀ᵀRx₀ = 0 (symbolic identity), hence q_j = x₀ᵀQ_jx₀ = 1/2 exactly.
(iii) **ρ-form**: g₀ = 1/f₀ − 2r₊z₁/(z²−z₁²), r₊ = ψ'(α)/(c·B_e'(α)), B_e'(α) = 2α/(1−α⁴),
    z₁ = ψ(α). With (i)+(ii):
    **ρ = 1 − B₁·r₊·z₁·(e₁²−e₂²)/((e₁²−z₁²)(z₁²−e₂²))** — the sym3-shaped target.
    VALIDATED vs exact machinery to 6-7 digits: (1.5,2.0,0.9,c=.12): ρ=+7.65e-4;
    (1.2,1.7,2.1,c=.06): ρ=+4.19e-5.
(iv) **EL4 (the even-phase inequality, uniformizer form)**. Parametrize by U ∈ (0,K):
    w = √k·sn(U,m), z = f·sin(sU), s = π/(2K), m = k². Node U₂ ∈ (0,K), outer node pinned
    at U₁ = K (focus law τ₁ = √k). v = pseudo-hyp midpoint of (k, k·sn²U₂) [quadratic formula;
    NOT arithmetic in U — tested false], W = sn⁻¹(√(v/k), m), B₁ = (k−v)/(1−kv). Then
    **ρ = 1 − Θ, Θ = B₁·s(1−v²)·sin(sW)cos²(sU₂) / [2k·snW·cnW·dnW·cos(sW)·(sin²sW−sin²sU₂)]**
    Target: Θ ≤ 1 on (k,U₂) ∈ (0,1)×(0,K). **PROVED 2026-07-21** by the
    Schwarzian-D2 comparison theorem plus a positive Fourier expansion for the squared-ellipse
    map; see `proof/el4_schwarzian_theorem.md`. Earlier verification: 50×50 grid (0 negatives,
    min 1−Θ = 2.4e-15 at degenerate corner); edge asymptotics 1−Θ ~ c·ε⁴ at U₂=(1−ε)K
    (checked dps 50/80: +7.9e-19, +7.9e-27 — dps-25 sign flips were cancellation noise, P3);
    **U₂→0 edge = sym3/Landen theorem EXACTLY** (Θ(k,0⁺) = π/(2K(k₁)), 12 digits).
    So EL4 is the 2-parameter generalization of the Landen theorem with a marked point.
(v) Effective-modulus structure: κ_eff (K(κ_eff)=π/(2Θ)) ≈ k₁·cd²((1+k′)U₂/2, m₁) to 4–6
    digits but NOT exact (Û/K₁ − U₂/K ~ 1e-3): the second Landen descent is leading-order
    only — the marked point breaks it. No clean pure-modulus closed form exists.
(vi) **D2-convex FALSIFIED** (important negative result): the coordinate-free form
    ρ = 1 − (δ/2)[1/(ζ₁−ζ₀)+1/(ζ₀−ζ₂)]/h'(ζ₀) with h:Ω̃→D, h(ζ₀)=0, h(ζ±)=±δ suggests the
    "harmonic-mean-of-difference-quotients ≤ derivative" inequality for convex Ω̃. It is FALSE:
    Koebe violates (univalent case, 1+δ²); z+z²/4 violates exactly (V = 1/(1−δ²/16)); 55+/60
    random Herglotz-generated convex maps violate (up to 1.26). Equality holds for ALL Möbius h
    (proved via midpoint relation — disk rigidity ρ=0). The truth requires the criticality
    coupling: outer node PINNED to the focus image (U₁ = K endpoint). Fourth confirmation of
    the "no soft proof" phenomenon (after L10, free-S≤1, L12-strong).

## 5. Phase geography (slice, empirics via full solver best_extremal)
- ODD interior phase (e.g. (2.0,1.2,1.6,.08), (1.8,1.0,1.4,.15)): zeros {0,±α}, α ∈ (τ₂,τ₁);
  ρ = B₁g₁q₁+B₂g₂q₂ with dominant positive outer term.
- EVEN phase (e.g. (1.5,2.0,0.9,.12), (1.2,1.7,2.1,.06) — a₂ > a₁ cases): zeros {±α},
  theory in §4. (In deg-3 solver runs the third zero migrates to ∂D = removable.)
- MÖBIUS phase = α→1 endpoint of odd family, B → −w, f₀ = ±φ (e.g. (2.2,0.9,1.3,.10)):
  K = 1.7750, ρ = +3.02e-2, both nodal terms positive; q = (0.708, 0.292).
- Family-restricted α-searches (slice_exact.extremal) can land on non-global stationary points
  in the even/Möbius phases — always cross-check with best_extremal (this caught cases 3,4,5).
- Pipeline→exact consistency: K(inflate→0) ↗ exact-ellipse K (1.811→1.896 vs 1.8987 at
  inflate .02→.0005), same α. Validates slice_exact against the certified pipeline.

## 5b. OFF-SLICE GENERALITY (2026-07-21, late): D2-crit
The L16 mechanism is domain-independent, and VERIFIED off the elliptic slice:
- L16(i)-(ii) proofs never used the ellipse: per-block norm is symmetric+even in (B₁,B₂) for
  ANY complementary 2×2 idempotent frame (κ₁ = κ₂ always), so the midpoint solution is always
  stationary; q = 1/2 by the involution identity. Any real even deg-2 Blaschke = b_v∘(φ²)
  (b_β·b_{−β} = b_{β²}(w²)), so the even phase is always deg-1 in the collapsed variable.
- Bi-conic check (sweep case ws=(1.235,.091,1.931,.035,1.965,.165), K=1.798): true extremal
  zeros ±0.34525 (deg-3 third zero = boundary junk); midpoint law v_mid = α² to 2.3e-7;
  b_v(u₁) = −b_v(u₂) to 5e-7; the symmetric critical point IS the global max off-slice too.
- D2-form ρ = +1.18e-5 (independent Newton-ψ + finite-diff φ′) vs pipeline ρ = +2.13e-5:
  sign/order agree; exact match precision-limited (1−Θ ~ 1e-5 needs φ′ to 1e-6 relative).
- DATA BUG FLAG: sym4_sweep_s61.jsonl OTHER-branch records have broken taus (τ₂ = 0 exactly)
  and inconsistent f0e — do not use those fields for OTHER records; recompute.
**D2-crit (general even-phase target):** for the even-phase class (any n, active block 2-dim),
Ω = int W(A) critical (Kippenhahn: nodes = real foci of the boundary-generating curve),
ρ = 1 − (δ/2)[1/(ζ₁−ζ₀)+1/(ζ₀−ζ₂)]/h'(ζ₀) ≥ 0 for the midpoint-normalized h = b_v∘χ.
Proved: 2×2 (both nodes at foci), sym3 (Landen; node at focus + center), slice-EL4 verified;
bi-conic verified (sign/order). FALSE without foci-pinning (§4(vi)). The abstract question
"why do Kippenhahn-foci nodes make it true" is the distilled deep quest for the even sector.

## 6. Reduced open targets (ranked)
1. Prove that a global extremal can be chosen with definite parity, or explicitly classify and
   cover symmetry-breaking four-node Pick data. Symmetry of the objective alone is insufficient.
2. Odd-phase positivity: ρ = B₁g₁q₁+B₂g₂q₂ ≥ 0 given §3 stationarity (3-parameter).
3. Möbius-phase positivity (deg-1 f₀ = φ on ellipse domains — may admit a general theorem
   beyond 4×4: only uses domain ellipse + deg-1 extremal + nodal weights).
4. De-symmetrize off the slice (bi-conic Kippenhahn data replaces the exact ellipse).
