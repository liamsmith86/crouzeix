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
| arXiv:2410.10678 (rev. Feb 2025) | Parameterized extension of C–P framework | Scaled version of C–P | TO READ: does the family bottom out above 2? |

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

## To read next
- [ ] arXiv:2410.10678 full text
- [ ] MMOR 2024 (C_N < 1+√2) — extract their equality-case analysis of C–P
- [ ] MMOR 2025 configuration constants
- [ ] Ransford–Schwenninger construction details (what does the sharp abstract example look like — can it be realized by an actual (A, Ω)? If YES → route to counterexample intuition; if NO → the gap is real structure to exploit)
