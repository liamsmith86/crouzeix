# LITERATURE_LEDGER.md

Papers, theorems, assumptions, possible gaps. Status verified via web 2026-07-20.

## Core chain
| Ref | Result | Method | Notes / gaps to exploit |
|---|---|---|---|
| Crouzeix 2004 ("Bounds for analytical functions of matrices") | Conjecture stated; proved for 2×2 (constant 2, sharp) | Direct estimates, conformal maps | 2×2 case: W(A) is an ellipse; sharp via Jordan block [[0,2],[0,0]] |
| Crouzeix 2007 | Universal constant 11.08 | Cauchy integral + conformal splitting | Superseded |
| Crouzeix–Palencia 2017 (SIMAX/FoCM) | W(A) is a (1+√2)-spectral set | f(A)+g(A)* via positive double-layer kernel; abstract lemma | THE target to improve. g = Cauchy transform of f̄ conjugate |
| Ransford–Schwenninger 2018 (SIMAX 1708.08633) | The abstract lemma (‖f(A)+g(A)*‖≤2‖f‖, g=C(f̄)) cannot yield better than 1+√2 in abstract setting | Explicit abstract construction | KEY OBSTRUCTION: need concrete structure of the pair (f(A), g(A)) beyond the abstract hypotheses |
| Malman–Mashreghi–O'Loughlin–Ransford 2024 | For each fixed N: C_N < 1+√2 | Compactness + strict inequality analysis of C–P equality conditions | Non-constructive; no uniform bound. Their equality analysis of C–P may reveal extremal structure |
| MMOR 2025 | Configuration constants via Neumann–Poincaré operator; domain-dependent bounds | Double-layer potential spectral analysis | Domain-dependent c(Ω) < 1+√2 for smooth domains? Check exact statement |
| Schwenninger–de Vries, arXiv:2409.15954 "DLP for spectral constants revisited" | §6: ρ(f₀,x₀) := ∫Re(K_Ω(f₀)*f₀)dμ₀ = ReC; **K² + ρ ≤ 2K** (Thm 6.1: K ≤ 1+√(1−ρ)); ρ ≥ −1 (convex); ρ ≥ 0 ⟹ Crouzeix; disk ⟹ (K−1)ρ = 0; \|ρ\| ≤ a(Ω) (Rem 6.2) | Extremal pair + DLP | READ & VERIFIED (pp. 15–17). Closest baseline to campaign's P1; our increment = β-subtraction (q < 2). Also: Appendix A smooth-approximation lemma (useful for Ω ↓ W(A) rigor); Prop A.2 spectral-constant equivalence |
| arXiv:2410.10678 | "The algebraic numerical range as a spectral set in Banach algebras" | — | CORRECTED ID (audit): earlier ledger entry misattributed this ID to a parameterized C–P extension |
| O'Loughlin–Rani, arXiv:2603.15536 (2026) "q-Numerical Ranges and Spectral Sets" | Bound K = t + √(t²+a(Ω)), t = 1+γ(1)/2 | q-numerical range | CAUTION (audit): sentence after Thm 3.3 claiming K<2-type conclusion may be typo/gap — verify before use |

## Special cases proved (constant 2)
- Normal matrices (trivial: R ≤ 1; equality-2 impossible). von Neumann: W(A) ⊇ spectrum.
- All 2×2 matrices (Crouzeix 2004).
- n×n with W(A) a disk (Badea–Crouzeix–Delyon? — Okubo–Ando: numerical radius ≤ 1 ⇒ disk is 2-spectral... actually unit disk is a complete 2-spectral set for w(A)≤1, Okubo–Ando 1975).
- Nearly Jordan blocks (Choi–Greenbaum).
- Weighted shift matrices (Choi; complete 2-spectral via arXiv:2508.12768 Aug 2025).
- Certain 3×3: tridiagonal with elliptic W(A) centered at eigenvalue; 3×3 KLS matrices (2025-26 work).
- Compressions of the shift / model space operators with certain Blaschke conditions (Bickel–Gorkin school).
- Matrices attaining ‖A‖ = numerical radius conditions (Crabb-type equality cases).

## Numerical evidence
- Greenbaum–Overton (2017-2018, "Numerical investigation of Crouzeix's conjecture"): nonsmooth optimization (BFGS/GRANSO) over A and p; global max of R appears to be exactly 2, attained on Crabb-type / Jordan-block configurations. No R > 2 ever found.
- AIM workshop reports (aimath.org/pastworkshops/crouzeixrep.pdf): consensus directions.

## Known equality structure (R = 2 attained/approached)
- A = [[0,2],[0,0]], p = z. R = 2 exactly. W(A) = disk radius 1.
- Crabb matrix family (nilpotent Jordan-like with specific superdiagonal weights √2,1,...,1,√2? — verify): ‖p(A)‖ = 2 with p = z^k, W(A) = unit disk. Attains 2 for higher powers.
- General principle: equality cases known are nilpotent + monomial + disk. Conjecturally all extremals are these (Greenbaum–Overton observation).

## Equivalent / stronger formulations (Track F cautions)
- Completely bounded version: "W(A) is a complete 2-spectral set" — open, possibly strictly stronger; Paulsen theory: cb-version ⇔ similarity to operator with dilation... (K-spectral ⇒ complete K'-spectral with K' possibly larger). Do NOT conflate.
- Equivalent: enough to prove for f holomorphic on int W(A), continuous on closure (Mergelyan/convexity).
- Equivalent: enough for A with W(A) having smooth boundary (perturb A ↦ A + εB dense case? verify continuity argument), and p with ‖p‖_{W(A)} = 1.
- Reduction: can assume p has degree < n (Cayley–Hamilton changes sup — INVALID as direct reduction; but valid: sup over W(A) of the reduced rep may increase, so inequality for reduced rep does NOT imply original. Only the direction ‖·‖ side is fixed.) — see APPROACH_LEDGER pitfalls.

## Read & digested (Epochs 1–3)
- [x] Ransford–Schwenninger arXiv:1708.08633 (full) — abstract lemma, sharpness (non-unital, Ω ⊅ W(T)), Question 4.1.
- [x] Schwenninger–de Vries arXiv:2302.05389 (full) — extremal pairs/measures, Theorem 5, Prop 8/9, Remarks 12–13.
- [x] MMOR arXiv:2407.19049 (pp. 1–8) — configuration constants c_R = c_C, a(Ω) < 1, ellipse formula, thin-domain obstruction, Theorem 6 curvature bound.
- [x] BGG+20 arXiv:2006.04901 (pp. 4–12) — Blaschke extremal structure (Crouzeix 04 Thm 2.1), cancellation Thm 4.1, extremal measure Thm 4.5, compressed shifts.
- [x] Schwenninger–de Vries arXiv:2409.15954 (pp. 15–20) — see table (P1 baseline).
- [x] NIST DLMF §23.8.1 — Weierstrass Fourier expansion used to prove the positive series in
  EL4; combined with `wp(x)=1/sn²(x|m)−(1+m)/3` for half-periods `(K,iK′)`.
- [x] Mashreghi--Moucha--O'Loughlin--Ransford--Roth, arXiv:2506.23831 — Schwarz--Jack
  convexity: the inverse Riemann map of a bi-circular domain (hence an ellipse) is convex on the
  positive radius, so the direct map is concave. Used only for `r≤p` in L24; no quantitative
  two-node estimate is imported.
- [x] Kanas--Sugawa (2006), *On conformal representations of the interior of an ellipse* — the
  centered inverse ellipse map has positive odd Taylor coefficients; an alternative precise
  justification of `r≤p`. L26 uses this coefficient positivity to obtain the cubic minorant
  that closes the `(p−r)/c` upper-block bound; the subsequent Möbius barrier and theta
  certificate are campaign-derived.
- [x] Crouzeix--Greenbaum, arXiv:2508.12768 (2025/26) — for scalar cyclic weighted shifts the
  disk complete bound equals the largest consecutive weight product via an explicit diagonal
  similarity. L25 uses the elementary nilpotent specialization on the `c→0` face; the paper's
  general theorem does not cover the interior block-weighted slice.

## To read next
- [ ] GKL arXiv:1701.01365 (Glader–Kurula–Lindström 2018, 3×3 tridiagonal elliptic W(A)) —
  novelty calibration vs our sym3 ρ-proof (their mechanism is cb/dilation-based; ours is new)
- [x] Kenan Li, 2021 UW thesis, Chapter 2.3 — for the proportional 3×3 elliptic family it
  explicitly identifies three extremal phases: centered degree 1, shifted real Möbius degree 1,
  and even degree 2. This invalidates any symmetry-only parity inference and directly predicts
  the shifted phase now found on the 4×4 slice. The thesis discusses `A=C+bC*` as the general
  construction of elliptic numerical ranges but does not visibly prove the arbitrary 4×4
  weighted-shift slice; a more targeted novelty audit remains.
- [ ] de Vries thesis (SV24 companion) — extremal-pair machinery details
- [ ] MMOR 2024 (C_N < 1+√2) — equality-case analysis of C–P (relevant to P2 stability)
- [ ] SV24 full §§1–5 (their Prop 2.4, 2.8, 3.5 used in §6; Berger–Stampfli generalization)
- [ ] R–S Question 4.1 citation trail — is unitality known to be insufficient? (still undetermined)
- [ ] BCG+23 arXiv:2312.04537 (Crouzeix, compressions of shifts, nilpotent classes)
