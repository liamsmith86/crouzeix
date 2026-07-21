# CROUZEIX RESEARCH CHECKPOINT — NOT A FINAL RESULT
(Restart packet per goal.txt forced-stop protocol; written proactively 2026-07-20 while work continues.)

## Exact current frontier
Campaign has reduced Crouzeix's conjecture to **H-r: ρ = Re∫f₀·Φ(f₀)dμ ≥ 0 at extremal pairs on
the critical domain Ω = int W(A)** (⟹ conjecture via SV24 Thm 6.1 + shrinking; proof/strategy_S.md,
proof/rho_positivity_program.md). H-r: PROVED for all 2×2 (closed form ρ = 1 − π/(2K(m)), modulo
the α=0 symmetry step); PROVED-by-reduction for the sym3/GKL 3×3 elliptic class (proof/sym3_reduction.md
— ζ = z² collapse, midpoint law, π₀ = 1/2); NUMERICALLY CONFIRMED beyond all known classes (sym4;
adversarial floors positive for n=3,4,5). No counterexample to the conjecture found anywhere
(all "violations" ever seen were certified numerical artifacts — APPROACH_LEDGER pitfalls P1–P7).

## Strongest proved lemmas
P1 (K²+ρ ≤ Kq), L13 (Clark-type transition ⟨q(A)x₀,u₀⟩ = K∫q·conj(f₀)dμ), L14 + collapse theorem
(sym3: v = α², ρ = (α²/2)(g₀(e)−g₀(0))), **Landen theorem** (sym3 ρ = 1 − π/(2K(k₁))),
**L15** (level-4 nodal/Pick closed form + odd stationarity law + frame identities),
**L16** (even-phase midpoint law + q = 1/2, domain-general), ceiling K²+2ρ+G² ≤ 4,
2×2 closed form; D2 partial results (symmetric-node case, wedge family, convex trace bound).
See LEMMA_LEDGER.md.

## Failed approaches — do not repeat
Scalar reduction L10 (disproved, odd-symmetric mechanism); global-operator conditions like
∮λ_min(P)ds (not ⊕-stable); soft/free-map versions of the confocal inequality S ≤ 1 (false:
1.019 free convex, 1.24 partial-focal); Löwner/K-monotonicity constraints (vacuous);
L12-strong annihilation (false); frozen-Hadamard route (pair-response dominates);
frame-free level-4 laws + PSLQ hunts (frame-coupled — L15 explains); soft D2 classes:
free univalent (Koebe 1+δ²), convex (z+tz² exact), odd+G'-increasing (certified 40-dps
counterexample) — see proof/D2_landscape.md; pure-modulus double-Landen formula for the slice
(marked point breaks second descent); arithmetic-in-U midpoint identity (false, 1e-2 off).

## Unresolved candidate lemmas
H-r general; **EL4** (even-phase level-4 elliptic inequality — verified, central); odd-phase
level-4 positivity; 2×2 α=0 step; ellipse-squared analytic proof (elementary); contact-degeneracy
write-up; D2 sharp-class question (open classical problem, NOT needed for campaign).

## Best numerical assets (experiments/)
extremal_pullback.py + theodorsen.py + minkowski_test.best_extremal = exact extremal-pair machinery
(certified; best_extremal is the phase-safe solver — pitfall P5). ellipse_sandbox.py = exact 2×2.
slice_exact.py = exact elliptic sym4 slice (25 dps; family-restricted — cross-check phases).
slice_invariants.py = L15 invariant formulas (K closed form 1e-16, law 1e-10, frame identities).
sym3_sweep_s51.jsonl (40 rec); sym4_sweep_s61.jsonl (20 rec; OTHER-branch taus/f0e BUGGY, ρ/K/qs
fine). Hr_adversarial floors: n3 +2e-4, n4 +1.7e-2, n5 +1.6e-2 (n6 dense: NOT OBTAINED — job killed
as futile, see pitfall P7 + COUNTEREXAMPLE_SEARCH; use structured families instead).
Direct ratio searches: n=6 → 1.148, n=7 → 1.465.

## Epoch-5 late findings (2026-07-21, after Landen theorem)
- Landen theorem BANKED (proof/landen_theorem.md): sym3 ρ = 1 − π/(2K(k₁)) classical-complete.
- Elliptic sym4 slice b_j = c·a_j: W(A) exact ellipse CONFOCAL WITH OUTER EIGENVALUES ±e₁
  (τ₁ = √k focus-image law persists); inner pair ±e₂ interior nodes.
- Level-4 stationarity: c₁·∂_αB(τ₁)/e₁ + c₂·∂_αB(τ₂)/e₂ = 0 with frame coefficients
  c_j = (h_jx₀)⟨g_j,y₀⟩; NEAR-balance c₁ ≈ −c₂ (⟺ ⟨Ax₀,u₀⟩ ≈ 0 via L13 moment condition,
  0.2–3%) but NOT exact — frame-free law loses 2e-5–6e-4 in K. Exact closed form must carry the
  algebraic frame system (heavier than sym3 but finite).
- Tools: experiments/slice_exact.py (exact elliptic slice machinery, 25 dps).

## Epoch-5 LATE additions (2026-07-21, second work block)
- L15 (nodal/Pick closed form, odd law) + L16 (even-phase midpoint + q=1/2) PROVED;
  EL4 = ρ = 1−Θ(k,U₂) ≤ 1 is the even-phase central target (verified thoroughly; contains
  the Landen theorem at U₂=0). proof/slice_closed_form.md is the master note.
- Slice phase geography: odd/even/Möbius(=odd α→1 endpoint). Sweep OTHER-branch taus/f0e
  fields are BUGGY — recompute, don't trust.
- D2 soft-class hunt CLOSED (all falsified, incl. odd+G'↑ via certified 40-dps
  counterexample); symmetric-node case + wedge family + convex trace bound PROVED.
  proof/D2_landscape.md. Remaining route: deformation-path (dV/dk ≤ 0 with explicit
  kernel functional on elliptic velocity field). GATE-CHECK PASSED: dV/dk < 0 strictly on full grid (−1.5e-3..−3.7e-2, bounded away from 0 in interior) — robust target.

## Next five concrete actions (refreshed 2026-07-21 late; old items 1–2 DONE → L15/L16/EL4)
1. **EL4 deformation-path derivation** (top): compute dV/dk analytically for the ellipse family —
   the kernel functional dV = h'(v) − (Möbius-Jacobian two-interval average of h')
   (proof/D2_landscape.md, last section) applied to the elliptic velocity field h = ∂ψ_k/∂k;
   prove dV/dk ≤ 0 (gate check: strictly negative on full grid). Then lift to the squared-ellipse
   (EL4 proper: proves even-phase slice H-r) and bi-conic domains (general level-4 even phase =
   FIRST NEW Crouzeix class).
2. **Odd-phase level-4 positivity**: ρ = B₁g₁q₁+B₂g₂q₂ ≥ 0 under the L15 stationarity law
   (slice_closed_form.md §3; interlacing τ₂ < α < τ₁, dominant positive outer term observed).
   Try the same kernel/deformation machinery; the Möbius-equality structure should persist.
3. **n=6 H-r floor, redesigned** (old dense run KILLED as futile — pitfall P7): run the
   adversarial min-ρ search on STRUCTURED n=6 families (zero-diag tridiagonal; graded collapse
   makes extremality certifiable via blocks), or upgrade the solver first.
4. Rigor debts: 2×2 α=0 symbolic step; contact-degeneracy write-up; ellipse-squared analytic
   proof; read GKL arXiv:1701.01365 + UW thesis (novelty calibration for elliptic n≥4).
5. General-n: L16 mechanism is dimension-generic (any 2-dim active block) — formulate the
   even-phase tower induction; n=6 even phase = deg-2-in-χ two-zero problem (next rung).
   De-symmetrize using ρ-margins (free-node robustness on ellipse/sq-ellipse = stability budget).

## Paste-ready continuation instruction
"Continue the Crouzeix campaign in /home/liam/Downloads/crouzeix (git repo; commit+push after each
task). Read RESEARCH_STATE.md (NEWEST section first), then proof/slice_closed_form.md and
proof/D2_landscape.md; master program in proof/rho_positivity_program.md. Resume at restart-packet
action 1 (EL4 deformation-path derivation). Respect APPROACH_LEDGER.md pitfalls P1–P7: every
numerical claim needs the certificate battery; treat any apparent violation as artifact until it
survives strict re-evaluation and an independent implementation; cross-check extremal phases with
best_extremal (P5); verify analyticity/univalence of probe map families (P6). Do not re-open dead
ends listed in the ledgers (esp. soft D2 classes — all falsified with certificates). The
stop-condition remains: rigorous general proof or certified counterexample."
