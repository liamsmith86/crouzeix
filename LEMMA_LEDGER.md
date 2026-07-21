# LEMMA_LEDGER.md

Candidate lemmas. Status ∈ {PROVED, DISPROVED, NUMERICAL, HEURISTIC, UNRESOLVED}.
Every status change requires: precise hypotheses, tests on small/defective/highly-nonnormal matrices, equality cases, dimension-dependence check.

| ID | Statement | Status | Evidence / notes |
|---|---|---|---|
| L0 | (C–P 2017) For convex Ω ⊇ W(A), f ∈ A(Ω): ‖f(A) + g(A)*‖ ≤ 2‖f‖_{∂Ω}, where g = C(f̄) conjugate Cauchy transform | PROVED (literature) | Positivity of operator double-layer kernel; total mass 2 |
| L1 | (C–P abstract lemma) ‖a+b*‖≤2, plus contractive f↦g structure ⇒ ‖a‖ ≤ 1+√2 | PROVED (literature); SHARP abstractly (R–S) | Cannot improve without concrete structure |
| L2 | Normal A: R(A,p) ≤ 1 | PROVED (classical) | spectral theorem + σ(A) ⊆ W(A) |
| L3 | Exact energy identity K² = ⟨γ_Φ(f₀)x₀, γ(f₀)x₀⟩ − conj(∫Φ(f₀)f₀dμ) at extremal pair | PROVED | Algebra from S–dV Prop 1/3 + multiplicativity. proof/refined_master_inequality.md |
| L4 | Free scalar subtraction: ⟨γ(f₀)x₀,(γ_Φ(f₀)−ω₀)x₀⟩ independent of ω₀ ∈ C | PROVED | Needs K>1 (Prop 1) |
| L5 | W² ≤ 2∮\|f₀\|²dν, ν=⟨Px₀,x₀⟩ds; W=2 ⟹ \|f₀\|=1 ν-a.e. | PROVED | Double Cauchy–Schwarz on PSD operator measure |
| L6 | K² ≤ K√(W²−\|β\|²) + c at extremal pair (K>1) | PROVED | Strictly refines S–dV Thm 5 (verified numerically: gap ≥ 0 always, =0 at exact extremality) |
| L7′ | √(W²−\|β\|²) ≤ 2 − c/2 whenever ⟨p(M)x₀,x₀⟩=0, ‖p‖_Ω ≤ 1, Ω convex smooth ⊇ W(M) | NUMERICAL (no verified violation; adversarial searches produced only quadrature artifacts, all killed by strict re-eval) | Suspect TRUE but possibly not provable without extremality; central variant is L7@extremals |
| L7@ext | √(W²−\|β\|²) ≤ 2 − c/2 at TRUE extremal pairs (f₀,x₀) with K>1 | **UNRESOLVED — CENTRAL TARGET** | ⟹ Crouzeix (with L6). Holds at every computed true extremal so far (slacks +0.11…+0.76) |
| L10 | Scalar reduction: 2∮\|p\|²dν − \|∮p̄wdσ\|² ≤ (2−c/2)² for admissible scalar (w,p) | **DISPROVED (2026-07-20)** | Odd-symmetric mechanism: ellipse b=0.6, η=(δ_x+δ_{−x})/2-Cauchy, p odd near-unimodular deg 5; viol +0.33, clean certificates. Mechanism: Φ(p)(x) ≠ 0 off-disk while L5-slack W² < 2∮\|p\|²dν is REAL operator structure (V*-compression loss). Lesson: scalar data (ν,η) alone insufficient — operator alignment obstruction essential |
| L8 | No ⊕-unstable (global-operator-norm) sufficient condition can prove the conjecture | HEURISTIC | Formalize: for any candidate bound B(M) with B(M₁⊕M₂) ≥ B-degradation... make precise before use |
| L9 | Saturation rigidity: W=2 ⟹ P^{1/2}(σ)y = e^{−iψ}f₀(σ)P^{1/2}(σ)x₀ a.e. (y = normalized γ_Φ(f₀)x₀) | SUPERSEDED by Stinespring defect identities (proof/P2_target.md) | Exact: 4−q² = 4[‖(I−Q)Uξ‖² + \|⟨Uξ,ξ⟩\|²] |
| P1 | K² + ReC ≤ Kq, q = √(W²−\|β\|²), at extremal pairs (K>1) | PROVED (2026-07-20, phase-preserving refinement of L6; verified vs SV24 §6 which has K²+ρ≤2K) | β-subtraction is the new content over arXiv:2409.15954 Thm 6.1 |
| P2 | q ≤ 2 + ReC/2 at extremal pairs | **UNRESOLVED — CENTRAL TARGET (Epoch 3)** | ⟹ Crouzeix via (K−2)(K−r/2) ≤ 0. Difficult case r<0 (generic off-disk). Equivalent to Stinespring leakage bound. Sharp at Crabb |
| L6/L7@ext | (former central target) | SUPERSEDED by P1/P2 (audit-adopted: modulus wasted the ImC phase) | All L7@ext numerical support transfers a fortiori to P2 |
