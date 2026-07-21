# CROUZEIX RESEARCH CHECKPOINT — NOT A FINAL RESULT
(Restart packet per goal.txt forced-stop protocol; written proactively 2026-07-20 while work continues.)

## Exact current frontier
Campaign has reduced Crouzeix's conjecture to **H-r: ρ = Re∫f₀·Φ(f₀)dμ ≥ 0 at extremal pairs on
the critical domain Ω = int W(A)** (⟹ conjecture via SV24 Thm 6.1 + shrinking; proof/strategy_S.md,
proof/rho_positivity_program.md). H-r: PROVED for all 2×2 (closed form ρ = 1 − π/(2K(m)), modulo
the α=0 symmetry step); PROVED-by-reduction for the sym3/GKL 3×3 elliptic class (proof/sym3_reduction.md
— ζ = z² collapse, midpoint law, π₀ = 1/2); NUMERICALLY CONFIRMED beyond all known classes (sym4;
adversarial floors positive for n=3,4,5). No counterexample to the conjecture found anywhere
(all "violations" ever seen were certified numerical artifacts — APPROACH_LEDGER P1–P4).

## Strongest proved lemmas
P1 (K²+ρ ≤ Kq), L13 (Clark-type transition ⟨q(A)x₀,u₀⟩ = K∫q·conj(f₀)dμ), L14 + collapse theorem
(sym3: v = α², M = scaled 2×2 Jordan, ρ = (α²/2)(g₀(e)−g₀(0))), ceiling K²+2ρ+G² ≤ 4,
2×2 closed form. See LEMMA_LEDGER.md.

## Failed approaches — do not repeat
Scalar reduction L10 (disproved, odd-symmetric mechanism); global-operator conditions like
∮λ_min(P)ds (not ⊕-stable); soft/free-map versions of the confocal inequality S ≤ 1 (false:
1.019 free convex, 1.24 partial-focal); Löwner/K-monotonicity constraints (vacuous);
L12-strong annihilation (false); frozen-Hadamard route (pair-response dominates).

## Unresolved candidate lemmas
H-r general; sym4 closure (current task); 2×2 α=0 step; ellipse-squared analytic proof (elementary);
contact-degeneracy write-up; Lemma M (Minkowski monotonicity — superseded-ish by direct H-r work).

## Best numerical assets (experiments/)
extremal_pullback.py + theodorsen.py + minkowski_test.best_extremal = exact extremal-pair machinery
(certified). ellipse_sandbox.py = exact 2×2. sym3_sweep_s51.jsonl = 40-record verified dataset.
sym4_probe.py = new-territory probe (ρ > 0: +1.3e-4, +1.9e-3, +4.5e-4). Hr_adversarial floors:
n3 +2e-4, n4 +1.7e-2, n5 +1.6e-2 (n6 job may be incomplete).

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
  kernel functional on elliptic velocity field).

## Next five concrete actions
1. sym4 closure — SHARPENED by final Epoch-5 data (sym4_sweep_s61.jsonl, 20/20 ρ>0):
   NO equioscillation — K = ‖A₊F(N_even)‖ alone (single active 2×2 block, verified 5 digits;
   other block 4–400× smaller). ODD phase (16/20): B odd, zeros {0,±α}, α interlaces (τ₂,τ₁)
   [13/16; check 3 exceptions]. EVEN phase (4/20): sym3-like, q=1/2 exact. So: single-block 2-node
   problem with F_j = B(τ_j)/e_j linked through α; stationarity dK/dα = 0 on the 2×2
   K(α) = ‖A₊(F₁Q₁+F₂Q₂)‖; then ρ = Re[B(τ₁)g₀(e₁)q₁ + B(τ₂)g₀(e₂)q₂] ≥ 0 via
   interlacing sign structure. Beware: the critical domain is NOT elliptic — the ζ-domain conformal
   data (g₀-values) needs the bi-conic Kippenhahn structure (det even in t: t⁴+c₂t²+c₄).
   This would be the FIRST NEW Crouzeix class of the campaign.
2. Sym4 sweep (adapt sym3_sweep.py) to gather (zeros, ρ, block data); verify collapse numerically
   before proving.
3. Rigor debts: 2×2 α=0; ellipse-squared analytic; read GKL arXiv:1701.01365 (compare mechanism —
   their cb-proof vs our ρ-proof) and de Vries thesis.
4. De-symmetrization: perturb off the symmetric slice; use ρ-margin + smoothness to extend H-r
   locally (stability analysis around solved families).
5. General-n parity induction: zero-diag tridiagonal n×n (Greenbaum–Overton's empirical hard
   configurations are exactly weighted-shift-like!): the graded collapse f₀(A) = A·F(A²) is
   dimension-generic — chase the induction.

## Paste-ready continuation instruction
"Continue the Crouzeix campaign in /home/liam/Downloads/crouzeix (git repo; commit+push after each
task). Read RESEARCH_STATE.md, then proof/sym3_reduction.md and proof/rho_positivity_program.md.
Resume at Epoch 5 action 1 (sym4 closure). Respect APPROACH_LEDGER.md pitfalls P1–P4: every
numerical claim needs the certificate battery; treat any apparent violation as artifact until it
survives strict re-evaluation and an independent implementation. Do not re-open dead ends listed
in the ledgers. The stop-condition remains: rigorous general proof or certified counterexample."
