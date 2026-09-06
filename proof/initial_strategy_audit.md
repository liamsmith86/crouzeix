# Initial strategy audit — 2026-07-20

**Scope:** I read the current project Markdown files and audited the claims below against their algebra, the project’s reported experiments, and primary sources. I did not modify any existing project file.

## Bottom line

This is a credible, well-audited research programme, but there is not yet evidence that L7 is true or that a proof is close. The move from global/domain-only quantities to true-extremal, vector-localized operator structure is well motivated. The strongest next move is **not** the current absolute-value L7: retain the phase of the interaction term and attack a strictly weaker sufficient inequality.

## Highest-priority correction: replace `c = |C|` by `r = Re C`

Use the project’s notation

- `C := ∫ g₀f₀ dμ`, `r := Re C`, `c := |C|`;
- `q := sqrt(W²-|β|²)`.

L3 and L4, before taking absolute values, give the sharper real inequality

`K² + r ≤ Kq`.                                                     (P1)

Therefore Crouzeix follows from the **phase-sensitive target**

`q ≤ 2 + r/2`.                                                       (P2)

Indeed, (P1)+(P2) gives

`(K-2)(K-r/2) ≤ 0`.

Here `K-r/2>0` because `K>1` and `|r|≤|C|≤1`, hence `K≤2`.

Why this is materially better than L7:

- If `r≥0`, nothing more is needed: `q≤2≤2+r/2` already proves `K≤2`.
- If `r<0`, P2 asks for `q≤2-|r|/2`, whereas L7 asks for the stronger `q≤2-|C|/2`.
- Since `|r|≤|C|`, P2 discards no needed strength and can be much easier when `C` has a large imaginary part.
- At Crabb, `C=0` and `q=2`, so P2 remains sharp.

**Immediate numerical action:** re-score every stored true-extremal record using

`s_phase := 2 + Re(C)/2 - q`,

not only `2-|C|/2-q`. Test/falsify P2 exactly in the 2×2 ellipse sandbox before investing in a general proof.

## This phase observation has a direct literature match

Section 6 of Schwenninger–de Vries, *The double-layer potential for spectral constants revisited* defines

`ρ(f₀,x₀) = Re ∫ g₀f₀ dμ = Re C`

(after translating their `K_Ω(f₀)^*` to this project’s `g₀`) and proves `K²+ρ≤2K`. They explicitly note that `ρ≥0` would imply Crouzeix. This is the `q≤2` predecessor of P1. L6’s subtraction of `β` refines their `2` to `q`; the paper should therefore be treated as the closest baseline, not merely a future optional read.

Primary source: <https://arxiv.org/abs/2409.15954>, especially Section 6.

## Where the needed strictness actually lives

Epoch 2’s E1 says every true extremal is `f₀=B∘φ`, hence `|f₀|=1` on all of `∂Ω`. This agrees with the published extremal-function structure. Consequently L5 becomes

`W²≤2∫|f₀|²dν=4`

at every true extremal. Thus L5’s conclusion “`W=2` forces `|f₀|=1`” supplies no new rigidity there: unimodularity is already automatic.

The useful quantity is the **operator compression/alignment defect**, not scalar unimodularity. Normalize the positive map

`T(h) := (1/2)∮hP ds`, so `T(1)=I` and `γ_Φ(f₀)=2T(f₀)`.

In a Stinespring representation `T(h)=V*π(h)V`, put `U:=π(f₀)` and `ξ:=Vx₀`. Because `|f₀|=1`, `U` is unitary. With `Q:=VV*`, one has the exact identities

`4-W² = 4||(I-Q)Uξ||²`,

`β = 2⟨Uξ,ξ⟩`,

`4-q² = 4[||(I-Q)Uξ||² + |⟨Uξ,ξ⟩|²]`.                (D)

So P2 for the only difficult case `r<0` is equivalent to

`4-q² ≥ -2r-r²/4`.                                      (P2-defect)

This says the negative interaction `Re C` must force either leakage from the Stinespring range or nonzero self-overlap. That is the concrete quantitative version of L9. L10 failed precisely because scalar measures cannot see the leakage term.

**Suggested proof target:** derive P2-defect from the full true-extremal Euler–Lagrange equations—left/right singular vectors plus variations of the zeros of the finite Blaschke product—not from `ν,μ` alone.

## Adversarial audit of the current lemma chain

- **L3/L4:** Algebraically sound under a fixed inner-product convention. The project should write the phase-preserving real form P1 explicitly before any triangle inequality.
- **L5:** Sound, but its scalar rigidity consequence is vacuous at true extremals because E1 already gives boundary unimodularity. Its equality *defect* (D), however, is valuable.
- **L6:** Sound as a localized refinement. Calling it a refinement of the published master route is fair, but novelty still requires a full literature check. Its use of `+|C|` unnecessarily loses phase.
- **L7@ext:** Sufficient and sharp at known equality cases, but stronger than Crouzeix and stronger than P2. Crouzeix could be true while L7 fails. Exact falsification work should precede proof work.
- **L7′:** Do not prioritize the relaxed version. The campaign’s own meta-pattern and L10 show that first-order shadows admit irrelevant counterconfigurations.
- **L8:** The direct-sum diagnosis is sound for the examples discussed. The universal phrase “ANY global-operator-norm condition” is broader than what is proved; keep it heuristic.
- **L9:** The saturation condition is consistent with equality in operator-valued Cauchy–Schwarz. Replace a qualitative write-up by the exact Stinespring defect formula (D), then seek a uniform stability estimate.
- **D1:** The adjoint extremal measure can be justified by applying the extremal-measure proposition to the conjugate uniform algebra and `f↦γ(f)*`; write this carefully rather than relying on informal stationarity.
- **D2:** Correct, and actually an exact equality before absolute values: for `u₀=f₀(A)x₀/K` and `H=(f₀g₀)(A)`, `⟨Hx₀,x₀⟩=⟨Hu₀,u₀⟩` because `f₀(A)x₀=Ku₀`, `f₀(A)*u₀=Kx₀`, and the functional-calculus factors commute.
- **D3:** Correct; it follows by expanding `||f₀(A)x₀+g₀(A)*x₀||²`. Preserve `Re C` rather than immediately replacing it by `|C|`.
- **D4:** Plausible Cauchy-projection/principal-parts decomposition. State multiplicities and normalization of the principal parts explicitly and verify it first for one-zero Blaschke factors and the disk.
- **Two-regime strategy:** Interesting but premature. No far-from-disk estimate currently reaches 2, and the near-disk first variation is unknown. First determine whether P2 survives exact 2×2 extremals.

## Literature corrections and cautions

1. `LITERATURE_LEDGER.md` misidentifies arXiv:2410.10678. That ID is *The algebraic numerical range as a spectral set in Banach algebras*: <https://arxiv.org/abs/2410.10678>.
2. The described parameterized/scaled extension is apparently O’Loughlin–Rani, *q-Numerical Ranges and Spectral Sets*, arXiv:2603.15536: <https://arxiv.org/abs/2603.15536>.
3. Audit that 2026 preprint rather than trusting it blindly. Its displayed bound is `K=t+sqrt(t²+a(Ω))`, `t=1+γ(1)/2`. Algebra gives `K<2` iff `a(Ω)<-2γ(1)`. The sentence after its Theorem 3.3 claims a `K<2` conclusion, but the preceding theorem only bounds `γ(1)` and does not visibly establish this comparison. Treat that sentence as a possible typo/gap until resolved; if it were valid for every strict smooth neighbourhood, shrinking the neighbourhood would essentially settle Crouzeix.
4. The finite-Blaschke extremal statement for general simply connected smooth `Ω` is supported by BGG+20’s account of Crouzeix 2004: <https://arxiv.org/abs/2006.04901>. Cite the precise theorem and distinguish `H∞(Ω)` boundary limits from `A(Ω)` when needed.

## Recommended order of work

1. Translate Section 6 of arXiv:2409.15954 into `(C,r,q,β)` notation and record exactly what L6 adds.
2. Re-evaluate existing exact/near-exact extremals with `s_phase` and the defect components in (D).
3. Compute the true 2×2 ellipse extremal (`B` has degree 1) and attempt to disprove P2 over eccentricity and smooth outer approximations.
4. If P2 survives, derive the finite-Blaschke zero-variation and left/right singular-vector stationarity equations, then seek P2-defect.
5. Only after that, attempt near-disk/far-disk decomposition or more adversarial relaxed searches.

## Research-state hygiene

The ledgers currently disagree: `RESEARCH_STATE.md` still labels itself Epoch 1, `APPROACH_LEDGER.md` still calls relaxed L7′ the central target, `COUNTEREXAMPLE_SEARCH.md` has an empty run table despite reported runs, and literature items already read in Epoch 2 remain under “to read.” This does not affect the mathematics, but it can steer a continuing agent toward stale work. Synchronize them at the next safe checkpoint.

## Independent checks performed for this memo

- Expanded P2+P1 exactly: the residual factors as `(K-2)(K-r/2)`.
- Verified the L7/defect algebra and master factorization in 10,000 randomized scalar checks each.
- Verified D2 and D3 on 600 random commuting matrix-polynomial pairs each.
- Verified the Stinespring compression-defect identity (D) on 600 random finite positive operator-valued measures.
- Read the TeX sources of arXiv:2302.05389, 2409.15954, and 2603.15536 to check the extremal-measure theorem, the `ρ=Re C` connection, and the 2026 bound.

## Stop/go criterion

If P2 fails at an exact true 2×2 extremal, abandon the whole L7-style scalar tradeoff and return to the still-sharper un-Cauchy–Schwarz angle in L3. If P2 holds exactly in 2×2 and its defect tracks `-Re C`, the Stinespring stability route is the most structurally justified continuation.
