# The grouped prepared endpoint has sharp two-channel weight

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood and through a fixed terminal jet.
> It is not a general or dimension-uniform Crouzeix theorem.

## 1. Result (L317, 2026-07-26)

Retain the balanced pure partial-isometry colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\qquad B_j=W^*(S^*)^jV ,
\]

and the weighted transfer row

\[
\boxed{{\mathbb B}_L(c)
=[cB_1\ \ c^2B_2\ \cdots\ c^LB_L].}             \tag{1}
\]

Apply L298's retightening seeds and then the complete L306--L307
affine preparation recurrence: compulsory lower parallel
neutralization, fixed-base Stein solution, and L311's ideal-valued
response normalization are all included.  Group the two physical
ellipse orientations before taking the nonresponse upper endpoint.

Through every fixed terminal degree \(N\), that grouped endpoint has
the exact form

\[
\boxed{
{\cal U}_{\le N}(c)
={\mathbb B}_L(c)\{D_{\rm dir}+cR_N(c)\}
 {\mathbb B}_L(c)^*
{\cal M}_S(C_{\rm rsp}(c))+O(c^{N+1}).}          \tag{2}
\]

Here:

1. \(D_{\rm dir}\) is block diagonal and consists exactly of the
   principal same-grade direct Grams fixed by L283/L298;
2. \(R_N(c)=R_N(c)^*\) is a finite matrix polynomial whose blocks
   are locally bounded and analytic on the fixed repeated-block
   neighbourhood;
3. \(C_{\rm rsp}\) is L311-normalized, polynomial/rank-stable, and
   vanishes on the matching complete-delay strata; and
4. no transfer inverse, rank projection, or pseudoinverse is used.

Equivalently, before the regrouping in (2), every retained quadratic
block at elliptic degree \(d\) has the form

\[
\boxed{B_jR_{jk}^{(d)}B_k^*,\qquad j+k\le d.}    \tag{3}
\]

If \(j+k=d\), then

\[
\boxed{j=k=d/2,}                                 \tag{4}
\]

and the term is the principal direct Gram.  Every unequal pair and
every nonprincipal diagonal pair has integer slack

\[
\boxed{d-j-k\ge1.}                               \tag{5}
\]

Thus L217's genuine early frame term
\(-2c^4SWB_5\) causes no exception.  When it is paired with the
grade-one root at degree six, the apparent equality pair
\((B_1,B_5)\) belongs to the unequal zero-slack Hardy cross and
vanishes only after the complete physical endpoint is grouped.  It
would be false to prove (3)--(5) for the raw frame, raw operator, or
individual rooted summands separately.

L317 closes the **support and extra-\(c\)** half of the post-L316
gate.  It does not yet prove that the bounded Hermitian polynomial
\(R_N\) is absorbed by the available fixed-\(\theta=1/2\) margins.
That square completion, with L307's \(-C^*C\) reserve left intact,
is the next step.

## 2. Marked wandering-chain lemma

L313's recursion has a useful initial-height form.  Let \(w\) be a
word in

\[
S,\qquad J=(I+F)S^*(I+E),
\]

and let \(r(w)\) be the number of \(J\)'s.  For \(h\ge1\), put

\[
K_{h,w}=W^*(S^*)^h w(S,J)V .
\]

Then there are explicit copy matrices \(C_a^{(h)}(w)\) such that

\[
\boxed{
K_{h,w}
=\sum_{a=1}^{h+r(w)}B_aC_a^{(h)}(w),}            \tag{6}
\]

\[
\boxed{
\left\|\operatorname{col}_a C_a^{(h)}(w)\right\|
\le (|w|+1)4^{r(w)}.}                            \tag{7}
\]

Indeed, apply \(w^*\) to \(S^hW\).  L313 gives, for \(q\ge1\),

\[
\begin{aligned}
S^*S^qW&=S^{q-1}W-VB_{q-1}^*,\\
J^*S^qW&=S^{q+1}W+VB_{q+1}^*,
\end{aligned}                                    \tag{8}
\]

together with

\[
S^*W=0,\qquad J^*W=2SW+2VB_1^*.                 \tag{9}
\]

The surviving chain begins at height \(h\).  Each \(J^*\) can raise
it by only one, while \(S^*\) lowers it.  Every emitted transfer
therefore has index at most \(h+r(w)\).  Applying the remaining
suffix and compressing by \(V^*\) constructs (6).  The same count as
L313 proves (7).

Suppose a displayed root \(c^kB_k\) began at its sharp seed degree
\(2k\), and subsequent finite physical motion has excess elliptic
degree \(e\).  Every physical reverse edge costs one power of \(c\),
so \(r(w)\le e\).  Taking \(h=k\) in (6) gives

\[
a\le k+e,\qquad
k+a\le 2k+e=d.                                   \tag{10}
\]

This proves (3) for every finite physical-word sector.  More
generally, if the root already has slack \(s\), the same argument
starts with the accumulated budget \(k+s\); later multiplication
cannot decrease the slack.

## 3. Why the stable future row does not break the weight

The marked-word lemma must not be applied termwise to a stable Stein
sum.  Such an expansion contains arbitrarily long strings of
\(S,S^*\), and L259 shows that an isolated coefficient can display a
future \(B_a^*\) with no visible scalar weight.  The correct grouped
object is the whole causal leakage row.

L255 puts every finite-colligation dependence into the doubled-Hardy
Hankel channel

\[
{\cal H}_{n,q}=B_{n+q}^* .
\]

L272 normal-orders every closed pair of such crossings into the
Toeplitz leakage projection

\[
{\cal L}_B={\cal T}_B{\cal T}_B^*.
\]

Its exact block formula is

\[
\boxed{
P_a{\cal L}_BP_b
=\sum_{\ell<\min(a,b)}
B_{a-\ell}B_{b-\ell}^*.}                         \tag{11}
\]

On the grade-\(k\) flag quotient
\(B_1=\cdots=B_{k-1}=0\), this becomes L259's full future row

\[
P_k{\cal L}_BP_b=B_kB_b^*,\qquad b\ge k.         \tag{12}
\]

L251 proves that the two analytic chain ports occur in a closed
entry/return pair; a one-sided \(B^\sharp\) coefficient is not an
endpoint term.  In a term rooted at \(B_k\), use its Hardy row as the
row-height origin for the opposite port.  Moving that port outward
from row \(k\) to row \(a>k\) requires at least \(a-k\) forward
Hardy steps.  L243's inverse-model-kernel filtration and L245's
Green columns make the scalar payment exact:

\[
\rho^q=c^q\zeta^{-q},\qquad
\zeta^q\lambda_r(t)(1-t^q),\quad
\operatorname{ord}_c\lambda_r=r\ge q,            \tag{13}
\]

where \(q\) is the row distance.  Moving inward only lowers the
exposed transfer index and needs no estimate.  Therefore an
excess-degree-\(e\) closed root/port pair obeys

\[
a\le k+e,\qquad k+a\le2k+e.                      \tag{14}
\]

This includes the sharp early-frame example L217: the grade-one row
can reach \(B_5\) after exactly four outward steps.  The raw frame
coefficient looks one order early only when the already paid
grade-one root is omitted from the ledger.

All higher inverse-kernel powers in L243 add strictly positive
degree, and L258's renewal inverse has identity constant term.

Finally, L262 identifies the paired half-line metric as the ordinary
Hardy metric.  Hence the complete row in (12), not any one of its
future blocks, returns as

\[
\boxed{
P_k{\cal L}_B^2P_k=P_k{\cal L}_BP_k=B_kB_k^*.}   \tag{15}
\]

This is L261's matrix-inner Parseval identity.  Away from the exact
grade-\(k\) quotient, every additional summand in (11) contains a
strictly earlier channel \(B_{k-\ell}\), \(\ell\ge1\).  It is
therefore assigned to the already retained flag; lowering the
channel index cannot worsen (3).  This proves that stable
Stein/Hardy closure preserves the same weighted cone as (10).

## 4. Zero-slack rigidity

It remains to identify the equality case.  If \(j+k=d\), every
nonnegative valuation used above must be sharp:

1. there is no nonconstant inverse-kernel, renewal, metric, or graph
   correction;
2. every physical reverse edge follows the shortest wandering path;
3. the stable return uses the constant ordinary Hardy metric; and
4. no earlier leakage row from \(\ell\ge1\) in (11) is present.

Let \(k=\min(j,k)\) and pass to its ordered flag quotient.  A
zero-slack path from row \(k\) to the other exposed row has no spare
reflection or direction reversal.  In L262's ordinary Hardy
coordinates its middle return is therefore the Laurent shift
\(L^{j-k}\), or its adjoint if the orientation is reversed.  L266
gives the exact compression

\[
\boxed{
P_k{\cal L}_BL^h{\cal L}_BP_k=0\quad(h\ne0),
\qquad
P_k{\cal L}_B^2P_k=P_k{\cal L}_BP_k.}            \tag{16}
\]

Thus every unequal shortest return vanishes, while the zero shift
retains precisely the diagonal row energy.  This is stronger than
informally saying that distinct Hardy rows are orthogonal: it also
keeps the complete future leakage row on both sides and kills the
intervening shift.

L266 also covers every unilateral boundary word of depth at most
\(k\).  The depth-\(k+1\) divergence discovered in L269 is not
silently discarded: L277--L279 prove its all-grade grouped
delay-ideal cancellation, and L283 gives the resulting complete
matrix face.  That exceptional divergence belongs to the universal
first physical return.  Inserting any already prepared root into an
additional boundary loop adds a nonconstant inverse-kernel or moving
factor and therefore positive slack by L243/L315; it cannot create a
second zero-slack exception.  L255 requires the two orientations to
be combined before (16) is used, and L316
guarantees that the result is a paired same-endpoint expression.
L263 shows that motion of the Wold basis and Schur graph is one order
too late to alter this first new face.

For \(j=k\), equality forces \(d=2k\).  L315 proves that every
successor of a grade-\(k\) root gains at least one degree, while a
product of two already prepared roots has positive slack.
Consequently the only degree-\(2k\) same-row term is L298's seed.
L283/L298 compute its direct Gram and fix its normalization.  This
proves (4)--(5).

This argument also explains why the stronger statement

> degree \(d\) uses only \(B_1,\ldots,B_{\lfloor d/2\rfloor}\)

is false.  For example \(B_1B_3^*\) may occur at degree five:
its weight is four and it has one legal power of slack.  L317
controls total paired weight, not either channel index separately.

## 5. Closure through the affine preparation

The L298 seed formula is

\[
\begin{aligned}
C_kV^*
={}&-\frac12G_kG_k^*\\
&-\frac12Q\left\{
S^kG_k^*
+\sum_{j=1}^{k-1}(S^*)^{k-j}G_kG_j^*
\right\}.                                       \tag{17}
\end{aligned}
\]

At degree \(2k\), its paired indices satisfy \(j\le k\), with
equality only in \(G_kG_k^*\).  Thus the induction starts in the
cone (3)--(5).

L315 already proves that the displayed root and its elliptic degree
survive every operation in the complete recurrence.  The sharper
paired-weight statement follows from Sections 2--4:

1. two-sided finite physical multiplication uses (10);
2. stable Stein and graph closure use (11)--(15);
3. products add nonnegative slack;
4. lower parallel elimination is the row-zero compression
   \(-ER_dE/2\) and cannot increase a channel index; and
5. L311's relative-commutator normalization changes only the
   response representative, keeps the same physical root, and is
   normal-ordered with the grouped endpoint before the next step.

The last item is essential.  Individual words in L311's four-word
defect expansion need not satisfy the sharp bound; L272 and (11)
apply to their grouped closed return.  This is why abstract
two-sided-ideal membership alone would not prove L317.

All bounds are finite through a fixed jet.  L313 supplies the
finite-word constants, L311 supplies linear word-length response
bounds, L243 permits only finitely many reflection powers at a fixed
degree, and the Hardy/Toeplitz projections have norm one.  Stable
Stein completion is controlled by L315's local power moment
\(\Lambda_S\).  The resulting block coefficients are therefore
locally bounded and rank-stable.

## 6. Rees regrouping

Write one retained coefficient as

\[
c^dB_jR_{jk}^{(d)}B_k^*
=(c^jB_j)c^{d-j-k}R_{jk}^{(d)}(c^kB_k)^*.        \tag{18}
\]

By (3), the middle exponent is nonnegative.  Equations (4)--(5) say
that exponent zero occurs only in the principal diagonal block.
Collect those blocks in \(D_{\rm dir}\).  Every other exponent is at
least one, so factor one \(c\) and collect the remaining finite
polynomial in \(R_N(c)\).  Hermitian pairing gives
\(R_N(c)=R_N(c)^*\).  This proves (2).

## 7. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_two_channel_weight.py \
  --output \
  experiments/repeated_crabb_two_channel_weight_s70226.jsonl
```

The checker starts L313's recursion at heights one through four and
reconstructs every physical word through length seven on
unstructured and rank-chain colligations.  It also checks L245's
direct/reflected port valuations, ordinary Hardy-row orthogonality,
L266's nonconstant leakage-shift annihilation on noncommuting
Potapov products, and the final quadratic regrouping.  All 6,120
marked words pass;
the largest reconstruction error is \(4.62\cdot10^{-14}\), and no
coefficient exceeds the sharp support \(h+r(w)\).

The marked recursion, closed-row identities, and valuation induction
above—not the floating audit—prove L317.
The tracked dataset has SHA-256

```text
f3f24239f45beb756b4da21fd5ce7d7ce6fc25e662cefbfd344c395fde52ed51
```
