# Graded-induction skeleton for H-r on weighted-shift-type families (Epoch 5)

Target class: A n×n zero-diagonal tridiagonal real ("weighted shift + weighted down-shift"),
symmetries A ~ −A (D = alternating diag) and A real. Includes/perturbs the Crabb–Choi weighted
shifts = Greenbaum–Overton's empirically-hardest configurations. n = 2: proved (elliptic closed
form). n = 3: proved-by-reduction (= GKL class). n = 4: 20/20 numerical (this file's target).

## Structure theorems (established at n ≤ 4, conjectured generic)
1. Parity: extremal f₀ has definite parity; both parities occur as PHASES of the family
   (n=4: 16 odd / 4 even records). Phase boundary = ties of the two parity-optima.
2. Graded collapse: f₀ even ⟹ f₀(A) = G(A²) block-diagonal on sublattices; f₀ odd ⟹
   f₀(A) = A·F(A²) block-off-diagonal. Sublattice blocks of A² are (⌈n/2⌉, ⌊n/2⌋)-sized with
   common nonzero spectrum {e_j²}.
3. Single-block activity (n=4 verified): K = norm of ONE active block; the extremal problem is
   a ⌈n/2⌉-node interpolation problem over the ζ = z²-image domain — HALF the size. Induction!
4. ζ-domain geometry: Kippenhahn determinant is even in t ⟹ ζ-boundary = projection of a
   half-degree curve. n=3: exact ellipse (proved via det factorization t·(t²+c₂)). n=4: quartic
   oval, NEAR-ellipse (single-conic residual ~7e-5 at near-shift weights).

## Concrete near-term theorem target (sym4 perturbative closure)
Since the ζ-domain is an O(δ)-perturbed confocal ellipse (δ = second-sheet influence,
computable from weights), and the sym3/2×2 machinery gives ρ_ellipse ≥ margin(m) > 0 off-disk:
**prove |ρ − ρ_ellipse-model| ≤ C·δ with explicit C (stability of extremal data under domain
perturbation — Hadamard bounds), yielding rigorous H-r for the sym4 subfamily δ < margin/C.**
First NEW Crouzeix class. Components needed:
(a) quantitative version of the 2×2/sym3 formula giving margin(m) explicitly ✓ (have: 1 − π/(2K(m)));
(b) Lipschitz bound of ρ in Hausdorff-domain-distance at fixed matrix (Hadamard variation of Φ
    is bounded — have the single-layer formula; extremal-pair stability needs a spectral-gap
    assumption on the top singular value — generically valid, quantifiable);
(c) bound δ(weights) — explicit from the Kippenhahn quartic.
None of these require new ideas — "just" careful analysis. HIGH PRIORITY.

## The abstract quest (parallel track, unsolved)
A criticality-based proof of H-r without explicit maps. Known constraints on any such proof:
- must use Ω = W(A) (identities are criticality-blind; ρ < 0 occurs off-critical);
- must use operator structure (scalar L10 false);
- must degenerate correctly on the rich equality manifold (disks AND all Crabb-type + phase
  boundaries where ρ → 0);
- candidate mechanism: KKT alignment of the K-maximization against the annihilation manifold
  g₀(σ(A)) = 0 (codim-2 short), with the DLP contact kernel field (ker P(σ_θ) = span((σ_θ−A)x_θ))
  supplying the sign. Status: shape identified, proof missing.
