# APPROACH_LEDGER.md

Every attempted approach; why it succeeded/failed. Includes known-failed approaches from literature (do not repeat blindly).

## Known dead ends / pitfalls (from literature — do not repeat)
1. **Improving the abstract C–P lemma alone.** Ransford–Schwenninger: sharp at 1+√2 abstractly. Any improvement MUST inject concrete structure of (A, W(A), Cauchy transform).
2. **Invalid degree reduction via Cayley–Hamilton.** Replacing p by its remainder mod char poly changes max over W(A) (can increase it), so "wlog deg p < n" is FALSE for the sup side. Valid only as: ‖p(A)‖ computable from remainder.
3. **Pointwise-to-operator-norm leaps.** Scalar estimates on W(A) do not transfer to ‖·‖ for nonnormal A.
4. **Proving the completely-bounded version as a shortcut.** Possibly strictly harder; a cb-proof failure says nothing about scalar version.
5. **Von Neumann on a disk containing W(A).** Gives constant depending on disk vs W(A) geometry; radius arguments alone max out well above 2 for eccentric shapes.

## Campaign attempts
| ID | Track | Idea | Status | Outcome / failure condition |
|---|---|---|---|---|
| A1 | A | Reconstruct C–P line-by-line; locate constant entry | DONE | Constant enters as c ≤ 1 in K² ≤ 2dK + c (S–dV Thm 5). Full dissection in proof/track_A_crouzeix_palencia.md |
| C1 | C | Randomized + gradient search max R, dims 3–5, poly deg ≤ 6 | DONE (round 1) | Max R_inner ≈ 1.9747 (n=3,d=3). No violation. Optimizer stalls below 2 at nonsmooth points; landscape matches G–O |
| A2 | A | (⋆): c ≤ 2L with L = ∮λ_min(P)ds (CGL slack) | ABANDONED — WRONG SHAPE | Not ⊕-stable: M⊕(normal junk) sends L→0, keeps c. Global operator quantities cannot work. Lesson → localized quantities |
| A3 | A | L6 refined master inequality (vector-localized) + L7′ target | SUPERSEDED (→ P1/P2 → H-r) | L6 proved. Chain: L7′ → L7@ext → P1/P2 (audit) → H-r (SV24 §6 anchor) — current central target |
| C2 | C | Adversarial optimization of L7′ violation over (M,p,x₀) | DONE (no violation) | All positives were artifacts (P1–P4); machinery reused for H-r searches |
| A4 | A/B | H-r reduction: ρ ≥ 0 at extremal pairs on Ω = int W(A) ⟹ Crouzeix (SV24 Thm 6.1 + shrinking) | ACTIVE — CENTRAL | proof/strategy_S.md, rho_positivity_program.md |
| A5 | B | 2×2 closed form ρ = 1 − π/(2K(m)) on confocal ellipses | PROVED (mod α=0 write-up) | rho_2x2_theorem.md; K>1 hypothesis essential |
| A6 | B | sym3 ζ=z² collapse → Landen theorem ρ = 1 − π/(2K(k₁)) | PROVED | landen_theorem.md; collapse ≡ Landen transformation |
| A7 | B | Level-4 nodal/Pick closed form (L15) + phases + even-phase structure (L16/L18) | PARTIAL | `slice_closed_form.md`; even-sector midpoint globality and EL4 are proved. Parity/symmetry-breaking classification and odd/Möbius positivity remain |
| A12 | B | Globality of the even midpoint over the full two-node Pick body | PROVED | Boundary data reduce to constants or disk automorphisms; exact determinant shows the automorphism orbit is uniquely maximized at its hyperbolic midpoint when K>1. `proof/even_pick_globality.md` |
| A13 | B/G | Assume Z₂ symmetry forces a parity extremal | DEAD — FALSE | Shifted degree-one Möbius phase on slice `(.8,2.4,1.1,.3)` beats both parity sectors. Must include symmetry breaking explicitly; matches the phenomenon in Kenan Li's sym3 phase diagram |
| A14 | F | Complete-2 similarity SDP for the elliptic slice | ACTIVE — SHARP REDUCTION | L21 dualizes the optimum exactly as `sup tr(Z−TZT*)_-/tr(Z−TZT*)_+`; L22 reduces every slice matrix to three parameters and 2×2 parity blocks. A deterministic 1215-point grid has `t<4` (max 3.9997713 near `c→0`). Remaining theorem: prove the trace ratio ≤4 using the coupled conformal nodes; simple diagonal/Gramian/equality ansatzes fail. `proof/slice_similarity_duality.md` |
| A15 | F | Decompose L20 into modal-block and parameter-boundary theorems | PARTIAL — BOTH BLOCKS PROVED | L23 gives the sharp nome bound; L24/L26 prove `||C||,||B||≤2`; L25 proves the full similarity bound on the sharp `c→0` weighted-shift face. L26 uses a Möbius barrier and a finite interval-certified theta lemma. Only genuinely coupled dual witnesses remain. `proof/slice_boundary_theorems.md`, `proof/slice_upper_block_theorem.md` |
| A16 | F | Eliminate the contraction LMIs with Stein defects and classify coupled KKT rank faces | ACTIVE — MATRIX INNER THEOREM LEFT | L27 parameterizes every contraction metric by two PSD 2×2 defects through a Szegő kernel. A17 closes rank-one/full. L40 reduces rank-one/rank-one to the norm of an explicit two-parameter `2×2` Blaschke--Potapov transfer; proving that norm ≤2 is now the only rank face. `proof/slice_coupled_defects.md` |
| A17 | F/B | Convert the rank-one/full coupled face to an odd two-node Blaschke block | **PROVED (2026-07-21)** | L28 identifies the face with an odd Blaschke block and proves the lower block. L30–L32 reduce the upper block to two endpoint residuals; L33 and L36–L39 now certify every nome, including the singular multiscale ridge. Hence L29 and the generic rank-one/full KKT face are closed. `proof/slice_odd_block_reduction.md` |
| A18 | F/B | Convert rank-one/rank-one to a matrix-valued transfer | ACTIVE — EXACT REDUCTION | Two spectral Gram rotations give `R_{a,b}(T)=−S⁻¹(T−A)(I−AT)⁻¹S`. Its modal node function is `2×2` rational inner with determinant `(ab−z²)/(1−abz²)`. L41 converts the norm test to a quartic polynomial block LMI and proves both diagonal blocks PSD. L42 closes `c≥2^(−2/3)`. L43 gives a small-gain form and numerical scalar target L44. L45 gives a second exact form; L47/L48 add `||B||||C||≤2` and close the complete axes `ab=0`. The unresolved transfer is now genuinely interior. `proof/slice_coupled_defects.md`; `experiments/slice_rank_one_transfer.py` |
| A20 | F/B | Couple the two parity-block norms through their common modal orientation | **PROVED MILESTONE (L47/L48)** | Exact traces/determinants reduce `||B||||C||≤2` to one rational scalar inequality. Elementary theta bounds and four blow-up charts give a 13-second exact integer Bernstein certificate. Orthogonal factorization then proves (RT) on `a=0` and `b=0`. This does not yet control `ab≠0`, but supplies the first sharp coupled invariant beyond the separate L24/L26 bounds. `proof/slice_coupled_defects.md` §6.3; `experiments/slice_block_product_certificate.py` |
| A21 | F/B | Resolve the zero-node transfer ridge with a sharper nome estimate | **PROVED BOUNDARY — L49/L50** | Jacobi's first two product factors give `k≤4c/(1+4c²)` for `c≤1/2` by a 24-coefficient exact Bernstein check. On `p=0`, BE reduces to one scalar inequality. The high half follows from `k≤1`; on the low half, ridge-centered and largest-coordinate blow-up charts give a finite exact rational Bernstein cover. Thus BE/RT holds on the complete zero-node face. The same centered ridge remains relevant for a small-`p` interior neighbourhood. `proof/slice_coupled_defects.md` §6.4; `experiments/slice_sharp_nome_certificate.py` |
| A22 | F/B | Prove the full polynomial transfer core by projective principal-minor charts | **ACTIVE — LOW CORE CLOSED; BRIDGE + COMPACT RANGE OPEN** | L55 proves every determinant/minor chart and both signs on `c∈[0,.005]`; scale-free nome/main arm charts remove the artificial det0 cutoff, and its global complements pass in `235/399` leaves. L54 covers `c∈[.01,.020736]`, all charts/signs, and L42 separately covers `c≥2^(−2/3)`. First bridge `.005→.01` by reusing L55. Then extend the directed five-chart cover across the genuinely unproved compact interval `.020736<c<2^(−2/3)`; fixed-`c` tests there are evidence, not a continuous theorem. Exact low-nome terms remain rational through 15 with Arb tails at 16/17. No negative core has appeared. `proof/slice_core_projective_reduction.md` |
| A19 | F/G | Test L21's contraction-similarity route on general matrices before further slice certification | **NUMERICALLY SURVIVES GENERAL GATE** | A certificate-gated sweep accepted 102/120 varied complex matrices (`n=3..8`) and found no `t*>4`; max `3.8399296`. The sharp near-Crabb case rises `3.68956→3.97783` as the outer offset shrinks `.02→.00125`, always from below. Seventeen polygonal near-normal maps and one inaccurate SDP were rejected, not counted. This makes L21 a credible general route but remains finite evidence for a potentially stronger complete-bounded statement. `proof/general_similarity_probe.md` |
| A8 | B | Frame-free level-4 stationarity laws / PSLQ closed-form hunts | DEAD | Law is frame-coupled; algebraic only after adjoining p_ij, δ (L15 explains) |
| A9 | B | Soft D2 classes: free univalent / convex / odd+G'-increasing | ALL FALSIFIED (certified) | D2_landscape.md; Koebe 1+δ²; z+tz² exact; odd+G'↑ via linearized step-deficit + 40-dps counterexample. 5th instance of the no-soft-proof meta-pattern |
| A10 | B | EL4 deformation path: dV/dk ≤ 0 with Möbius-Jacobian kernel on elliptic velocity field | SUPERSEDED — EL4 PROVED | Gate was correct, but Schwarzian comparison is stronger and shorter: `SG≥0 ⇒ D2`; squared-ellipse `SG≥0` is a positive Weierstrass Fourier series. `proof/el4_schwarzian_theorem.md` |
| A11 | B | Schwarzian/Sturm comparison for D2 | PROVED; ACTIVE FOR BI-CONIC | Closes EL4 analytically. Next test whether the off-slice bi-conic collapsed map has `SG≥0` on its critical real interval; if not, use the more general potential ordering directly |

## Methodological pitfalls discovered (Track G)
- P1. Offset-curve construction silently breaks on matrices with nonsmooth W(M) (near-normal, corners/flats): FFT ringing → garbage integrals that still satisfy holomorphic-identity checks. FIX: certify positivity of P(σ) pointwise + convex tangent winding + mass 2I + Cauchy reconstruction; reject otherwise. COVERAGE GAP: such matrices currently unsearchable — needs smoothed support-function domains.
- P2. Background process hygiene: stale processes with old code contaminated logs → phantom "violations". FIX: harness-tracked jobs only; verify claims by re-evaluating saved records with current code.
- P3. Catastrophic cancellation: optimizer drives poly coefficients to ~1e50 with boundary cancellation; consistency certificates fooled (both sides share garbage). FIX: theorem guards — computed ‖p(M)‖ > 1+√2 or ‖γ_Φ(p)x₀‖ > 2 under certified hypotheses PROVES numerics lie (they contradict established theorems) → reject. Plus coefficient-magnitude cap.
- P4. Near-singular layer-potential quadrature (outer→inner contour, gap ε/2): N=1024 insufficient; passed all other certificates while Φ values were wrong (fake viol +0.358 at n=3, killed at N=4096). FIX: unitality certificate Φ(1)=1 through the same kernel + N-vs-2N stability check.
- P5. **Family-restricted optimizers miss the global phase** (found 2026-07-21): scalar α-searches
  inside one Blaschke family (e.g. odd {0,±α}) can converge to non-global stationary points —
  the true extremal may be in another phase (even {±α}, or Möbius = odd family's α→1 endpoint).
  FIX: always cross-check with the unrestricted best_extremal solver (deg 1..n−1, multistart);
  boundary-hugging zeros (|z| ≈ 1) in its output are removable (unimodular factors), not real.
  Related DATA BUG: sym4_sweep_s61.jsonl OTHER-branch records have broken taus (τ₂ ≡ 0) and
  inconsistent f0e — recompute; the ρ/K/qs columns are fine.
- P6. **Inadmissible probe families in map-class falsification** (found 2026-07-21): test maps
  built from tanh/step smoothings can have poles INSIDE the disk (tanh(κ(z²−t₀²)) poles at
  z² = t₀²+iπ/(2κ)) — huge "violations" from non-analytic non-univalent junk. FIX: before
  trusting any V > 1, verify analyticity (pole locations) and univalence (Re G' > 0
  Noshiro–Warschawski, or boundary-curve simplicity). Entire (polynomial) smoothings +
  explicit Re G' bounds give certified counterexamples (D2_landscape.md).
- P7. **Certification capability gap at n ≥ 6 (dense)** (found 2026-07-21): the true-extremal
  criterion diag = |⟨f₀(A)x₀,x₀⟩| < 1e-5 is currently UNREACHABLE for arbitrary dense 6×6
  matrices — both find_extremal_blaschke and best_extremal stall at diag ~1e-3 on the deg-5
  Blaschke landscape (10 params), so any acceptance-gated loop (e.g. Hr_adversarial's uncapped
  `while cur is None`) spins forever at 100% CPU with zero output (n=6 s42 job killed after
  2h04m / ~500 silent rejections). FIXES: (i) cap rejection loops + log every rejection with
  (K, diag); (ii) at n ≥ 6 use STRUCTURED families (zero-diag tridiagonal: graded collapse
  reduces extremality to small blocks — certifiable, and they are the G–O-hard configurations
  anyway); (iii) or upgrade the solver (analytic gradients / phase-aware polish) before any
  dense-n≥6 campaign.
- P8. **A stronger sufficient target needs its own falsification gate.** L20/L21 prove a
  complete-bounded similarity statement, and L44 is stronger still than the exact small-gain
  inequality. Before investing in a sharp-edge certificate, test the similarity constant on
  general matrices and audit L44 at high precision in its equality corner. Failure of either
  stronger target is not evidence against scalar Crouzeix.
- **Meta-pattern (important):** every relaxation of extremality so far admits sharp counterconfigurations: abstract lemma (R–S, non-unital α), domain-only constants (MMOR, thin quadrilaterals), scalar localization (L10, odd-symmetric mechanism — ours), free-map S ≤ 1, soft D2 classes (convex; odd+G'↑). Conjecture-strength inequalities must engage TRUE extremal/critical structure (global maximality of f₀ + positive extremal measure + critical domain), not just its first-order shadows.

## Epoch 2 outcome
2×2/ψ-domain exact extremal machinery built (extremal_pullback.py). Empirical laws at true extremals: c small (0 on disk, ~ε²-ish off), G ≈ |β|, ReC ≤ 0 generic, L7@ext holds broadly. L12-strong false.

## Epochs 4–5 outcome (2026-07-20/21) — see RESEARCH_STATE.md NEWEST + proof/slice_closed_form.md
H-r reduction anchored (A4); 2×2 and sym3 PROVED in closed form (A5, A6 — Landen theorem);
level-4 closed-form theory built (A7: L15 plus conditional L16 structure); EL4 subsequently
proved by A11; soft-class hunts all falsified with certificates (A8, A9). Adversarial H-r floors
positive n=3..5; direct ratio searches n≤7 found no violation.

## Epoch 3 (historical): P2 target, audit-adopted (chatgpt/FABLE_RESEARCH_AUDIT.md — verified before adoption)
P1 (proved): K² + ReC ≤ Kq. P2 (target): q ≤ 2 + ReC/2 at extremals ⟹ Crouzeix. See proof/P2_target.md.
Attacks: (i) falsification-first: exact 2×2, Ω = W(A) ellipse (elliptic Riemann map) + ψ-family adversarial ascent on −s_phase = q − 2 − r/2; (ii) if P2 survives: Euler–Lagrange at extremals (Blaschke-zero variations + singular-vector stationarity) → Stinespring leakage bound; (iii) fallback: pre-Cauchy–Schwarz angle form of P1.
Audit adoptions verified by hand: P1/P2 algebra, Stinespring identities (D), D2-exactness, SV24 §6 baseline (ρ = ReC, K²+ρ≤2K — our β-subtraction is the increment). Audit hygiene items fixed (ledgers synced, 2410.10678 corrected, run tables filled).
