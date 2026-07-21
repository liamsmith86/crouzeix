# D2 landscape (Epoch 5, 2026-07-21 late) — the two-node distortion inequality dissected

Object: V(G; a,b) = (δ/2)[1/(G(a)−G(v)) + 1/(G(v)−G(b))]·(1−v²)G'(v) for univalent real-symmetric
G on D, real nodes a > b, v = pseudo-hyperbolic midpoint, δ = b_v(a). D2: V ≤ 1.
ρ (even-phase H-r) = 1 − V at critical configurations. Equality: ALL Möbius G (proved).

## Hyperbolic trace form
g(s) := G(tanh s) (real trace in hyperbolic arclength), V₀ = midpoint (arithmetic in s!),
D = half-distance, Δ± = |g(V₀±D)−g(V₀)|:
  **V = tanh(D)·[1/Δ₊ + 1/Δ₋]·g'(V₀)/2.**
The tanh(D) < D "hyperbolic discount" is what any proof must exploit (the undiscounted
HM-of-difference-quotients statement fails even for half-planes).

## Proved pieces
1. **Wedge family exact** (g' = e^{2γs}, γ = wedge exponent): V = tanh(D)·γcoth(γD) ≤ 1
   ⟺ γ ≤ 1, via γcoth(γD) ≤ coth(D) (monotone in γ) and tanhD·cothD = 1.
   Equality: γ = 1 = half-plane, ALL D (rigidity). Reflex wedges γ ∈ (1,2] violate for large D.
2. **Convex trace bound**: G convex ⟹ |d/ds log g'(s)| ≤ 2, sharp (half-plane ≡ 2). Proof:
   Herglotz 1+zG''/G' = ∫(1+ze^{−iθ})/(1−ze^{−iθ})dν, real slice ⟹ d/ds log g' =
   (1−t²)/t[∫P_t(θ)dν − 1] − 2t with Poisson P_t ∈ [(1−t)/(1+t),(1+t)/(1−t)]-range ⟹ ∈ [−2,2].
   (Note: this bound does NOT decide V — see taxonomy.)

## Violation taxonomy (all numerically certain)
- z + tz², any t > 0 (convex for t ≤ 1/4): V = 1/(1−t²δ²) > 1 at symmetric nodes. So
  CONVEXITY IS INSUFFICIENT (and |μ'| ≤ 2 is insufficient).
- Odd maps with G' decreasing on (0,1) (e.g. G' = (1−z²/4)^{1.05}, univalent by
  Noshiro–Warschawski): V = δ/G(δ) > 1 forced at symmetric nodes. So ODDNESS IS INSUFFICIENT.
- Random convex (Herglotz, off-axis prevertices): up to 1.45. Koebe: 1+δ².

## Survivor class (leading conjecture)
**D2-mono: G odd univalent real-symmetric with t ↦ G'(t) nondecreasing on [0,1) ⟹ V ≤ 1
for ALL real node pairs.** Evidence: ellipse maps (all k), odd Koebe z/(1−z²) (nonconvex!),
z/(1−cz²), artanh, arcsin, positive-exponent odd slit-maps — worst V = 0.99969–0.99992,
tight only at node-coincidence and Möbius degenerations. Also passes: squared-ellipse χ
(not odd itself — the square of an odd map; free nodes verified ⟹ EL4 with U₁ free).
Odd + G'↑ = "prolate orientation": derivative grows along the node axis.

## Application fit (why this suffices for the campaign)
Critical W(A)-domains in the symmetric tower: Ω centrally symmetric (A ~ −A), nodes =
eigenvalues on the focal/major axis = the G'-increasing direction; ζ-collapse gives squares
of such maps. Free-node robustness on ellipse/sq-ellipse ⟹ stability margin for
de-symmetrization. Chain: D2-mono (or D2-sq for squared maps) ⟹ EL4 ⟹ even-phase slice
H-r; D2-crit (bi-conic version) ⟹ general even-phase H-r at level 4.

## Symmetric-node theorem (PROVED)
For G odd with G'↑ on [0,1): symmetric nodes ±δ give V = tanh(D)g'(0)/Δ and
g'(s) = G'(t)(1−t²) ≥ G'(0)sech²(s) ⟹ Δ ≥ g'(0)tanh(D) ⟹ **V ≤ 1**, equality iff
G' ≡ G'(0) on (0,δ). One line. (Equivalently: G(δ) ≥ δG'(0) ⟸ G' increasing.)

## General-node status: quantitative regularity needed
Steep-ramp probe (G' = exp(β tanh(κ(t²−t₀²)))) produced huge V (up to 92) BUT the family is
INADMISSIBLE: tanh poles at z² = t₀² + iπ/(2κ) lie inside D — not analytic, not univalent.
Univalence forces |d/ds log G'| ≤ 4+2|t| ≤ 6 (pre-Schwarzian bound) — exactly the regularity
the ramps violated. Lesson: D2-mono as stated (odd + G'↑ only) is UNPROVEN for general nodes
and likely needs the univalent-derivative bound explicitly in any proof; no admissible
counterexample found (ellipse/odd-Koebe/artanh/arcsin/positive-exponent slit maps all pass).
Sharp-class boundary = open classical question; the CAMPAIGN only needs the critical
configurations (ellipse, squared-ellipse, bi-conic verified) + a stability margin.

## D2-mono FALSIFIED (2026-07-21, certified)
Linearization around the identity: V = 1 + ε·dV with
  **dV = h'(v) − (1/2x)∫₀ˣ[(1+vx)²/(1+vu)²·h'(m_v(u)) + (1−vx)²/(1−vu)²·h'(m_v(−u))]du**
(G = id + εh; Möbius-Jacobian-weighted two-interval average). The functional is linear in h',
so extreme rays = even steps h' = 𝟙(|t| ≥ w). CLOSED FORM: for m_v(−x) < w < v the average
F = (1+vx)/2 + (1−vx)²(v−w)/(2x(1−v²)) < 1 = h'(v) ⟹ dV > 0. Constructed admissible
counterexample: G' = 1 + ε·I_{t²}(8,8) (regularized incomplete beta — ENTIRE polynomial step;
Re G' ≥ 1/2 on D ⟹ univalent by Noshiro–Warschawski; ε = 8.6e-7 forced by Bernstein growth
max|S| = 5.8e5 on |z|=1): **V(0.96, 0.4446) − 1 = +9.784e-8, certified at 40 dps** — matches
the predicted linear-order excess. So odd + G'-increasing is INSUFFICIENT; only the
symmetric-node case survives (proved above). The sharp class is quantitative: h'-variation
must be dominated by the Möbius-Jacobian averaging kernel (BMO/Muckenhoupt-flavored).

## The surviving proof route for EL4 (deformation path)
V ≡ 1 at k = 0 (disk) and worst-V DECREASES in k numerically (0.99998 → 0.99970 over
k = 0.4 → 0.9, free nodes). So EL4-free-nodes ⟺ dV/dk ≤ 0 along the disk→ellipse
deformation, and dV/dk has the SAME kernel structure with h = ∂ψ_k/∂k = explicit elliptic
velocity field. Target: show the ellipse velocity field satisfies the kernel inequality
h'_k(v) ≤ (Jacobian-average of h'_k) for all (v,x,k) — one explicit elliptic-function
inequality, linear in the velocity at each k. Same route for squared-ellipse (EL4 proper)
and bi-conic (level-4 general even phase). This is the concrete replacement for the dead
soft-class hunts.
