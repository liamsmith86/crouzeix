# Candidate reduction of the delayed trace law to endpoint filtration

> **Closure note (2026-07-25).**  The endpoint-filtration route below
> remains conditional as stated, but L279 proves its scalar conclusion
> independently: the delayed dual-Schur trace is
> \(2\|B_k\|_F^2\) and the effective separator trace is
> \(-16\|B_k\|_F^2\) in every grade.  Do not resume the conditional
> support classification merely to recover that scalar result.

## 1. Status (candidate, not a proved lemma, 2026-07-24)

Let \(S\) be a finite pure partial isometry with equal orthogonal
defects

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad EF=0,
\]

and transfer coefficients

\[
B_j=W^*(S^*)^jV.
\]

Let \(\mathcal K_S(c)\) be L225/L228's right-defect Schur residual
of the boundary-metric elliptic Stein slack.  Suppose

\[
B_1=\cdots=B_{k-1}=0,
\]

If the closed endpoint-word filtration stated in Section 2 is proved
for the **post-Schur** residual, then

\[
\boxed{
\operatorname {tr}[c^{2j}]\mathcal K_S(c)=0
\quad(0\le j<k),}                                  \tag{1}
\]

and

\[
\boxed{
\operatorname {tr}[c^{2k}]\mathcal K_S(c)
=2\|B_k\|_F^2.}                                   \tag{2}
\]

Sections 3--4 prove this implication from the filtration, but the
filtration itself still has the explicit gap described below.
Accordingly (1)--(2) remain a candidate arbitrary-grade trace law,
not a theorem.  They are weaker than the matrix identity

\[
[c^{2k}]\mathcal K_S
=E_1F_{k-1}+F_{k-1}E_1,
\]

which remains open beyond the finite exact audit.  Proving the trace
law would nevertheless be sufficient for L222's pointwise
semidefinite alternative.

Indeed L225's exact lower and upper boundary faces, followed by lower
retightening and its dual trace identity, would then give

\[
\boxed{
\operatorname {tr}E_{2k,\mathrm{eff}}
=-16\|B_k\|_F^2}                                   \tag{3}
\]

on every completely delayed reducing summand.  Thus every nonzero
L222 reducing separator would see a strictly negative face at its
first active transfer grade.

## 2. The Newton--Hankel endpoint filtration

We first record the filtration implicit in L125 and made state-level
exact by L236--L245.

> **Closed endpoint-word candidate.**  After subtracting L240's
> coisometric half-line baseline, a scalar state trace of total
> \(c\)-degree \(d\) is a finite sum of copy-space traces of closed
> endpoint words.  Give \(B_j\) and \(B_j^*\) endpoint weight \(j\).
> Every such word has total endpoint weight at most \(d\).  A word
> with no endpoint letter is zero after the right-defect Schur
> quotient.

Here “closed” means invariant under independent changes of the two
defect frames.  The lemma applies coefficientwise formally, so no
convergence or interchange of an infinite trace is involved.

More explicitly, let \(\mathfrak F_d\) be the span of scalar cyclic
words in \(B_j,B_j^*\) whose sum of subscripts is at most \(d\),
together with the empty word.  The assertion is

\[
[c^d]\operatorname {tr}\mathcal K_S\in\mathfrak F_d,             \tag{4}
\]

and its empty-word component is zero.  This formulation prevents a
hidden appeal to commutativity or to a scalar transfer.

To prove it, use the following exact building blocks.

1. **Scalar Newton cost.**  L125 proves

   \[
   [w^{2n+1}]\phi_c(w)
   =c^n h_n(c^2),                                  \tag{5}
   \]

   with analytic \(h_n\).  Thus an endpoint propagation of Hardy
   length \(n\) cannot occur below elliptic weight \(n\), and every
   correction has two additional powers of \(c\).  Expanding
   \(\Xi^{2n+1}\) into nearest-neighbour Hardy paths and closing a
   path with L236's two defect frames turns each endpoint crossing
   into exactly one cell \(B_j\) or \(B_j^*\).  The power of \(c\)
   supplied by (5) and by the reverse letters is at least the sum of
   the crossed Hardy distances.  Hence the scalar functional
   calculus sends coefficient degree \(d\) into
   \(\mathfrak F_d\).

2. **Metric cost.**  L236 writes the metric as two scalar diagonal
   Hardy weights whose cross Gram has cells

   \[
   (\mathcal O_R\mathcal O_L^*)_{a,b}=B_{a+b}^*.
                                                               \tag{6}
   \]

   Its coefficients are powers \(q^a=c^{2a}\).  A raw cross cell may
   lie arbitrarily far along the opposite frame, but it is not a
   scalar word by itself.  Closing it forces the complementary
   Hardy path; the Newton cost in item 1 (or, equivalently, L245's
   exterior column) supplies that remaining distance.  Hence the
   resulting cyclic word, rather than an isolated raw metric block,
   obeys the asserted cost.

3. **Resolvent cost.**  L243 writes the entire retained reflection
   through

   \[
   K_{u^rB}(\rho,x)^{-1},\qquad
   K_{u^rB}
   =\sum_{a<r}t^aI+t^rK_B.                         \tag{7}
   \]

   Expanding
   \(B(\rho)=\sum_{j\ge1}\rho^jB_j\) and its
   coefficientwise adjoint preserves endpoint order.  L243 proves
   that only the constant and first nonconstant inverse-kernel
   sectors can reach the requested window; every later sector has
   strictly larger cost.

4. **Full state lift.**  L245 shows that every omitted chain factor
   is a scalar power of \(x\) or \(\rho\), times
   \((1+t^a)/(1+t^r)\).  It introduces no new copy-space word and no
   negative elliptic valuation.  Thus closing the exterior state
   columns into a state trace preserves the cost from (4)--(7).

5. **Stein and Schur operations.**  Multiplication adds endpoint
   weights.  The right pivot has constant coefficient \(I_m\), so
   its formal inverse also preserves the filtration.  L244 says
   that the word-free half-line slack is the normalized defect Gram
   of a coisometry.  Its Schur quotient is zero; hence a surviving
   word must contain endpoint data.

The first, second, fourth, and fifth facts are exact.  The third does
**not yet prove** the claimed post-Schur valuation.  In L243,
\(\mathcal B^\sharp(x)=\sum_{j\ge1}x^jB_j^*\) carries no visible
power of \(c\).  Consequently later \(B_j^*\) can occur in the
pre-Schur zero/one-kernel sectors before the naive symmetric endpoint
weight predicts.  L244 strongly indicates that these terms form an
updated coisometric defect Gram and are removed by the Schur square,
and every finite audit shows that cancellation, but the state-lifted
Gram identity has not yet been written down.

L258--L259 subsequently identify the algebraic shape of this gap:
the final output row is first resummed into a closed return, and the
unweighted future series is exactly one row of
\({\cal L}_B={\cal T}_B{\cal T}_B^*\).  L261 proves that the
ordinary Hardy norm of that **entire** row is already its diagonal
energy:

\[
\operatorname {tr}(P_k{\cal L}_B^2P_k)
=\operatorname {tr}(P_k{\cal L}_BP_k)
=\|B_k\|_F^2.
\]

L262 proves that every genuine entry/return sandwich through L244's
normalized half-line uses exactly this ordinary unweighted Hardy
metric.  L272 now proves the placement part of the former obligation:

\[
\boxed{\text{place L258's first physical return between L262's
Wold unitary and its adjoint.}}                                    \tag{7a}
\]

namely, every active term has the leakage-sandwich form with a
copy-scalar middle return.  L263 additionally proves that the fixed physical Schur port and the
analytic motion of the background defect graph/Wold basis cannot
change this first new face.  Thus (7a) may be checked in the constant
Hardy compression after the background graph is removed.  What
remains is the full two-orientation physical numerator, not a hidden
metric or oblique-port weight.

L264 permits the complete right-half-line \(D_R\) coefficient to be
restored at the active degree without changing the scalar target,
while the active left orbit, final corner, and cross row remain
deleted.  Hence the numerator in (7a) should be assembled in that
port-isolated right-half-line gauge.  An arbitrary retained
coefficient is not free; its response is
\(-\operatorname {tr}(EX)\).

The remaining content of (7a) is support, not placement: identify the
copy-scalar middle return, its lower vanishings, and its first
bilateral constant modulo trace-null deep divergences.  L273 proves
the retained-metric cancellation of every interior fan path, leaving
only the unweighted full-fan theorem.  L260 still disproves the
shortcut of checking the two formal ellipse orientations separately.

## 3. Classification of the first trace face

Assume \(B_1=\cdots=B_{k-1}=0\) and (7a).  By the closed
endpoint-word filtration,
no nonconstant closed word has degree below \(2k\), proving (1).

At degree \(2k\), a surviving word has exactly two endpoint letters.
Both must have index \(k\); any later index or any additional letter
would exceed the cost.  Therefore

\[
\tau_k(B_k)
:=\operatorname {tr}[c^{2k}]\mathcal K_S(c)        \tag{8}
\]

is a real homogeneous quadratic function of \(B_k\), independent of
all later transfer coefficients.

Changing defect frames independently,

\[
V\longmapsto VU,\qquad W\longmapsto WZ,
\]

does not change the state operator or (8), while

\[
B_k\longmapsto Z^*B_kU.                            \tag{9}
\]

Thus \(\tau_k\) is biunitarily invariant.  The only real homogeneous
quadratic form on \(M_m(\mathbb C)\) invariant under
\(B\mapsto Z^*BU\) is a scalar multiple of the Frobenius square:

\[
\tau_k(B)=\alpha_{k,m}\operatorname {tr}(B^*B).    \tag{10}
\]

For completeness, diagonal phase matrices in (9) kill cross terms
between distinct matrix entries and distinguish real from
non-invariant complex bilinears; row and column permutations make
all surviving \(|B_{ab}|^2\) coefficients equal.  This proves (10)
without invoking an invariant-theory theorem.

## 4. All-grade normalization

Use L241's length-\(k\) monomial channel.  There

\[
B_k=I_m,\qquad F_{k-1}=E_1,
\]

and L241 proves in every grade that

\[
[c^{2k}]\mathcal K_{S_k}=2E_1.
\]

Taking the state trace gives

\[
\tau_k(I_m)=2m.
\]

Equation (10) gives \(\tau_k(I_m)=\alpha_{k,m}m\), so

\[
\alpha_{k,m}=2
\]

for every \(k,m\).  This proves (2), conditional on (7a).

Finally, L225's boundary faces are

\[
\operatorname {tr}E_{\rm up}=-4\|B_k\|_F^2,\qquad
\operatorname {tr}E_{\rm low}=+\|B_k\|_F^2.
\]

Lower retightening has
\(\Delta_{\rm low}=-\|B_k\|_F^2\), while L225's dual trace formula
and (2) give

\[
\Delta_{\rm up}-4\Delta_{\rm low}
=-4\operatorname {tr}[c^{2k}]\mathcal K_S
=-8\|B_k\|_F^2.
\]

Hence \(\Delta_{\rm up}=-12\|B_k\|_F^2\), and adding the boundary
upper face proves (3).

## 5. Adversarial regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_delayed_trace_law.py \
  --output \
  experiments/repeated_crabb_delayed_trace_law_s70224.jsonl
```

The checker uses generic random tails, not monomial channels.  It
evaluates the full theta/ODE direct map, metric, Stein product,
pivot inverse, and Schur square.  Grades one through five overlap
L228's exact noncommutative audit; grades six through eight extend
the independent floating falsification range.  It also records the
stronger matrix covariance error.  These audits do not prove the
remaining support part of (7a).
The tracked SHA-256 is
`1c2307a744c396c10ba8e0937ad0e4f2a3f00cb7a6aa5a8ae1cde7df4c1ea0a0`.

L272 proves that the complete first-active return is copy-scalar
between two leakage projections.  L266 then applies if its shallow
part is a Laurent series \(\Psi_{\rm phys}(L)\), because

\[
P_k{\cal L}_B\Psi_{\rm phys}(L){\cal L}_BP_k
=[\Psi_{\rm phys}]_0P_k{\cal L}_BP_k.
\]

The same conclusion holds for a copy-scalar shift polynomial whose
unilateral boundary depth is at most \(k\); depth \(k+1\) is a genuine
obstruction.  L273 removes the metric half of this depth calculation
by cancelling every interior fan word.  The remaining scalar
calculation is the unweighted full-fan theorem, its lower vanishings,
and hence bilateral constant \(4\).  Neither the finite audits nor
L266 establishes that unweighted support theorem.
