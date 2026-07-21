# RESEARCH_STATE.md — Crouzeix Conjecture Campaign

**Last updated:** 2026-07-21 (Epoch 5, late — level-4 closed-form theory + D2 landscape)

## NEWEST (2026-07-21 late session) — read proof/slice_closed_form.md + proof/D2_landscape.md
- **L15 PROVED**: level-4 nodal closed form K² = (T+√(T²−4δ²F₁²F₂²))/2 (frame invariants
  p_ij, δ only; verified 1e-16); odd-phase stationarity law explicit (1e-10); Hellmann–Feynman
  frame identities. Extremal problem = 4-point Pick problem with explicit objective.
- **L16 PROVED**: even phase = deg-1 Möbius in collapsed variable; midpoint law (symmetry
  proof) + q = 1/2 (involution identity) — DOMAIN-GENERAL, verified off-slice (2e-7).
- **EL4 (central target)**: even-phase ρ = 1 − Θ(k,U₂) explicit elliptic formula; ≤ 1 verified
  (50×50 grid, ε⁴-edges at dps 80); U₂→0 edge = Landen theorem exactly. Θ-form validated
  vs exact machinery (6-7 digits) and on a bi-conic case (sign/order).
- **Phases on slice**: odd {0,±α} / even {±α} / Möbius (= odd family's α→1 endpoint);
  family-restricted searches can miss the global phase — always cross-check best_extremal.
- **D2 landscape mapped, soft routes ALL DEAD (certified)**: free-convex false (z+tz² exact,
  Herglotz-random ≤1.45); odd+G'-increasing FALSE (linearized step-deficit closed form +
  certified counterexample V−1 = +9.8e-8 at 40 dps, entire-poly step, Noshiro–Warschawski);
  symmetric-node case PROVED (one line); wedge family exactly solvable (γ ≤ 1 ⟺ pass);
  convex trace bound |d/ds log g'| ≤ 2 proved (Poisson). Ellipse/sq-ellipse pass with FREE
  nodes (focus-pinning not needed — mechanism is the map class regularity).
- **Surviving proof route for EL4**: deformation path — V ≡ 1 at disk, show dV/dk ≤ 0 with
  the explicit Möbius-Jacobian kernel functional applied to the elliptic velocity field
  ∂ψ_k/∂k (linear in velocity at each k). See D2_landscape.md end.

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

## MILESTONES (Epoch 5, 2026-07-21) — THE LANDEN THEOREM
- **sym3 ρ-positivity now CLASSICAL-COMPLETE closed form** (proof/landen_theorem.md):
  W(A) ellipse explicit (semi-axes from support function; foci ±e proved); z₁ = e/√2 EXACT
  (quarter-period u₁ = K/2, sn(K/2) = 1/√(1+k′)); r₁ = eπ/(2√2 K k);
  **ρ = 1 − π/((1+k′)K(k)) = 1 − π/(2K(k₁)) ≥ 0** (Landen descent k₁ = (1−k′)/(1+k′)),
  equality iff disk. Verified: independent mpmath 30 dps + Landen identity 1e-31.
- ζ = z²-collapse ≡ Landen transformation ⟹ **Landen-tower conjecture** for the symmetric
  tridiagonal family (induction over levels; needs elliptic range at each level).
- **Elliptic sym4 slice found analytically: b_j = c·a_j** (= c-deformed weighted shifts!);
  W(A) exact ellipse (verified 3e-15). [CORRECTED later: foci = ±e₁ = OUTER eigenvalue pair,
  τ₁ = √k — see proof/landen_theorem.md; an earlier foci² = e₁²+e₂² note was an algebra slip.]
  ρ > 0 on slice; no pure-modulus Landen formula (marked point breaks second descent) — the
  correct closed form is EL4's ρ = 1 − Θ(k,U₂) (see NEWEST above).
- Adversarial H-r floors (all positive): n=3 +2e-4, n=4 +1.7e-2/+7.2e-2, n=5 +1.6e-2/+1.8e-2
  (n=6 running as of 2026-07-21, bash job in experiments/log_Hr_n6_s42.txt). Direct ratio
  searches DONE: n=6 best 1.148, n=7 best 1.465 — nothing near 2.

## MILESTONES (Epoch 4, 2026-07-20)
- **sym3 (= GKL 2018 class, real slice) H-r PROVED-by-reduction**: exact ζ = z² transform onto the
  canonical 2×2 confocal configuration; stationarity = pseudo-hyperbolic-midpoint law; π₀ = 1/2
  exact; ρ = (α²/2)(g₀(e)−g₀(0)) = 1 − S ≥ 0 via the 2×2 elliptic theorem. NEW MECHANISM,
  known class (novelty calibrated vs arXiv:1701.01365). proof/sym3_reduction.md is the master note.
- Rigidity of the inequality: S ≤ 1 FALSE for free odd convex maps (1.019), false on partial-focal
  synthetic loci (up to 1.24) — full criticality essential; no soft proof exists.
- L13 Clark-type transition proved; L14 three-factor formula; collapse theorem v = α².
- **sym4 = FIRST BEYOND-LITERATURE TERRITORY: ρ > 0 confirmed numerically (3 cases); extremal
  Blaschke is ODD ({0,±α} zeros) ⟹ f₀ = z·F(z²) collapse exists; squared boundary NOT elliptic
  (resid 1e-4–1e-3).** Next: close sym4 analytically = first new Crouzeix class in campaign.

## Current next actions (Epoch 5, refreshed 2026-07-21 late — items 1–2 of the old list DONE, see NEWEST)
1. **EL4 via deformation path** (top target): derive dV/dk analytically — kernel functional
   (Möbius-Jacobian two-interval average, formula in proof/D2_landscape.md) applied to the
   elliptic velocity field ∂ψ_k/∂k; prove its sign. Gate check passed (dV/dk < 0 strictly,
   bounded away from 0). Do plain-ellipse family first, then squared-ellipse (= EL4 proper),
   then bi-conic (general level-4 even phase).
2. **Odd phase positivity**: ρ = B₁g₁q₁+B₂g₂q₂ ≥ 0 given the L15 stationarity law
   (3-parameter; interlacing τ₂ < α < τ₁; term-1 dominance observed). Try the same
   deformation/kernel machinery.
3. n=6 adversarial result (running; check experiments/log_Hr_n6_s42.txt; harness will notify).
4. Rigor debts: 2×2 α=0 step; contact-degeneracy write-up; ellipse-squared analytic proof;
   GKL comparison read (arXiv:1701.01365); UW thesis check (elliptic n≥4 novelty — is the
   elliptic sym4 slice class already in literature?).
5. General-n: parity-collapse induction (L16 mechanism is dimension-generic for 2-dim active
   blocks; n=6 even phase = deg-2-in-χ, needs the next rung); de-symmetrization using ρ-margins.
Keep committing+pushing after each task (user instruction).

## Files map (handoff-ready, 2026-07-21)
proof/ — read in this order for the current frontier:
  rho_positivity_program.md (MASTER program), slice_closed_form.md (level-4 theory: L15/L16/EL4/
  D2-crit), D2_landscape.md (soft-class falsifications + deformation route), landen_theorem.md
  (sym3 closed form), sym3_reduction.md (collapse theorem), sym4_program.md (family setup),
  rho_2x2_theorem.md (n=2 base), graded_induction_skeleton.md (general-n plan);
  background: strategy_S.md, P2_target.md, refined_master_inequality.md,
  epoch2_extremal_structure.md, track_A_crouzeix_palencia.md.
experiments/: crouzeix.py (basics: poly_A, nr_support, ratio_inner/outer, crabb_matrix),
  extremal_pullback.py (ψ-domain exact machinery), theodorsen.py (+GeneralPullback: arbitrary
  convex domains, unitality certificate), minkowski_test.py (best_extremal — THE solver; always
  use it to cross-check phases), ellipse_sandbox.py (2×2 exact), slice_exact.py (elliptic sym4
  slice, 25 dps; family-restricted — see pitfall P5), slice_invariants.py (L15 invariant
  formulas + verifications), sym3_structure/sweep/analysis.py, sym4_probe/sweep.py (OTHER-branch
  taus/f0e fields BUGGY — pitfall P5), Hr_test.py / Hr_adversarial.py (certified min-ρ search),
  zero_geometry.py (phi_of_points); historical: L10/L12/scalar/adversarial_L7/search.py.
Data: sym3_sweep_s51.jsonl (40 rec), sym4_sweep_s61.jsonl (20 rec, ρ column trustworthy).
Ledgers: LEMMA_LEDGER.md, APPROACH_LEDGER.md (pitfalls P1–P6 — READ BEFORE ANY SEARCH),
LITERATURE_LEDGER.md, COUNTEREXAMPLE_SEARCH.md. Audit: chatgpt/FABLE_RESEARCH_AUDIT.md
(reconciled 2026-07-20). Restart: checkpoints/RESTART_PACKET.md (paste-ready instruction).
