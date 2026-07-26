# The accumulated response coordinates are paid by the actual joint defect

> **Campaign scope.**  This is a fixed-\(L,m\), curvewise comparison
> on the scalar zero-Jensen branch.  It closes the quantitative
> interface between L324 and L325.  It does not yet assemble every
> circular/elliptic chart into the final neighbourhood theorem.

## 1. Result (L327, 2026-07-26)

Use L325's actual canonical disk metric.  For the top-cluster
direction \(k\) selected by a scalar norming input, put

\[
 \Delta_k=
 \langle D_Xk,k\rangle+1-\sigma_X(k).             \tag{1}
\]

This is the **directionwise** defect retained in L325's proof before
the last minimization over \(k\).  Along a ramified real-analytic
rank stratum, let

\[
 d\asymp\sum_j\|\eta_j\|^2                        \tag{2}
\]

be L324's accumulated L197 quotient loss, and let

\[
 \delta=
 1-\sup_{\|v\|=1}\sum_n|u^*B_nv|^2               \tag{3}
\]

be the scalar-channel leakage of the terminal equality-anchor
transfer in the channel direction corresponding to \(k\).

There are locally uniform constants \(c,C>0\) such that

\[
 \boxed{
 c(\delta+d)\le\Delta_k\le C(\delta+d).}           \tag{4}
\]

In particular, L324's completed response cost obeys

\[
 O((\delta+d)^2)
 \le C_1\Delta_k^2
 \le\varepsilon\,\Delta_k                         \tag{5}
\]

after one neighbourhood shrink.  L325's directionwise proof gives

\[
 4-\|P^{-1/2}f(C)P^{1/2}x\|^2
 \ge c_{L,m}\Delta_k,                              \tag{6}
\]

so a fixed fraction of (6) absorbs every stacked L324 response.

The last minimization
\(\Delta=\min_k\Delta_k\) is useful for a uniform statement but must
not be used to pay a response attached to a different direction.
Equations (4)--(6) retain the actual norming direction throughout and
remove that possible mismatch.

## 2. Both channel defects are distances to the same product-line set

Let

\[
 {\cal O}_C h=(V^*C^nh)_{n\ge0}.
\]

L325's observability Parseval identity says that \({\cal O}_C\) is an
isometry.  In
\(\ell^2(\mathbb N_0)\otimes\mathbb C^m\), let

\[
 {\cal P}=
 \{a\otimes v:a\in\ell^2,\ \|v\|=1\}.
\]

For a unit \(k\), covariance diagonalization gives exactly

\[
 \operatorname{dist}({\cal O}_Ck,{\cal P})^2
 =1-\sup_{\|v\|=1}
   \sum_n|v^*V^*C^nk|^2
 =1-\sigma_X(k).                                  \tag{7}
\]

At the L324 terminal equality anchor, write \(W\) for the left defect
frame and \(u\) for the selected output channel.  L201 gives

\[
 {\cal O}_C(Wu)=(B_n^*u)_{n\ge1}.
\]

Matrix-inner Parseval and the same distance calculation give

\[
 \operatorname{dist}({\cal O}_C(Wu),{\cal P})^2
 =1-\sup_{\|v\|=1}\sum_n|u^*B_nv|^2
 =\delta.                                         \tag{8}
\]

Thus the channel terms in (1) and (3) are not merely two analytic
functions with a common zero set.  They are squared distances to the
same closed product-line set.

## 3. The disk flag measures the change between the two observations

L197's successive triangular congruences are exact least-squares
orthogonalizations.  Transport the selected scalar state through
them and retract every active quotient column to zero while retaining
the terminal tangential/equality data.  This is precisely L324's
split-channel retraction.

The Hardy quotient here is the same one seen by the observability
map, not an unrelated norm with a matching scalar value.  L193's
Berger-dilation formulas (18)--(25) realize \(4I-P\) as the
orthogonal complement of the defect orbit.  L201's characteristic
kernel identity (15) identifies that defect-orbit quotient with the
pure-partial-isometry observability model.  Minimality makes the
identification unitary; changing to L197's ordered quotient columns
adds only its analytic triangular congruences.  Thus transporting the
L197 residual column is exactly transporting the off-terminal part of
\({\cal O}_Ck\).

L325 identifies the scalar disk endpoint on the same transported
state:

\[
 \langle D_Xk,k\rangle
 \asymp\sum_j\|\eta_j\|^2\asymp d.                 \tag{9}
\]

Because the flag is finite and every factored active block has a
positive analytic leading coefficient, the retraction and its inverse
have bounded analytic graph maps on the ramified stratum.  Hence

\[
 \|{\cal O}_Ck-{\cal O}_{C_{\rm eq}}(W_{\rm eq}u)\|
 \le C_0\left(\sum_j\|\eta_j\|^2\right)^{1/2}
\le C_1\sqrt d.                                  \tag{10}
\]

Concretely, one Schur step writes the transported observation as its
projection onto the previously active Hardy range, plus the quotient
column \(\eta_j\), plus the next residual.  Undoing the triangular
congruence multiplies those three pieces by bounded analytic graph
columns.  Iterating the **finite** flag gives

\[
 {\cal O}_Ck-{\cal O}_{C_{\rm eq}}(W_{\rm eq}u)
 =\sum_{j=0}^{r-1}G_j(s)\eta_j,                   \tag{10a}
\]

with \(\sup_s\|G_j(s)\|<\infty\).  Cauchy--Schwarz in the finite
sum proves (10).  This is the observation-vector version of L324's
statement that the accumulated quotient collection contains every
off-channel disk variable seen by the transported scalar state.

Here all tangential block-Toeplitz/model motion is retained in
\(C_{\rm eq},W_{\rm eq},B\); only L197's quotient coordinates are
retracted.  Therefore (10) does not incorrectly charge tangential
motion to the disk loss.

The load-bearing point is that the complete accumulated collection
is used.  A newest-quotient-only version of (10) is false whenever
the state has already passed through an earlier active Hardy range.

## 4. A frame-free distance comparison

Distance to any closed set is one-Lipschitz.  Put

\[
 a={\cal O}_Ck,\qquad
 h={\cal O}_{C_{\rm eq}}(W_{\rm eq}u),\qquad
 e=\|a-h\|.
\]

Equations (7)--(8) and the triangle inequality give

\[
 \left|\sqrt{1-\sigma_X(k)}-\sqrt\delta\right|
 \le e\le C_1\sqrt d.                             \tag{11}
\]

Young's inequality therefore yields both

\[
\begin{aligned}
1-\sigma_X(k)&\le2\delta+2C_1^2d,\\
\delta&\le2\{1-\sigma_X(k)\}+2C_1^2d.              \tag{12}
\end{aligned}
\]

Add \(d\) and use (9).  This proves (4).

This argument also explains why a merely Lipschitz comparison of the
unsquared channel amplitudes is sufficient: both L324 and L325 use
their **squares**.  There is no dangerous unpaid
\(\sqrt{\delta d}\) term; L324 already charges it by Young's
inequality.

## 5. Rank changes and exact equality

The comparison is used inside a contradiction/curve-selection
argument.  Finite ramification fixes:

1. the top spectral multiplicity;
2. the L197 Schur flag and its factored valuations; and
3. a maximizing scalar product-line frame.

After the powers of the curve parameter are factored, every active
coefficient is positive definite and every triangular graph map in
(9)--(10) is bounded.  If a rank drops, the next L197 quotient is the
new stratum; there are at most \(m\) such steps.  Thus (4) is
curvewise uniform exactly in the sense required to exclude an
analytic counterarc.  No discontinuous pseudoinverse or globally
smooth maximizing channel is used.

If \(\Delta_k=0\) and a scalar function actually attains norm two,
L326 supplies the reducing Gau--Wu disk-model summand and lowers copy
multiplicity.  Equality of (4) alone is not used to infer scalar
innerness.

## 6. Consequence

L324's response bound can now be inserted into L325 without the
unproved sentence “the two losses are the same” and without adding
reserves from different anchors.  On one scalar branch the retained
upper-gap expression has the schematic form

\[
 -c_0\Delta_k-y^*\Gamma y
 +2\operatorname{Re}\langle y,R\rangle,
 \qquad
 \|R\|\le C_2(\delta+d).                           \tag{13}
\]

Complete all components of \(y\) at once.  Equations (4)--(6) give

\[
 -c_0\Delta_k+C_3(\delta+d)^2
 \le-\frac{c_0}{2}\Delta_k                         \tag{14}
\]

after shrinking.  L199's circular-normal and L318's elliptic
curvatures can therefore be retained as separate diagonal margins.

The remaining task is the finite analytic merger across the Jensen,
disk, circular-normal, and elliptic charts; the disk/channel
response-size interface itself is closed.

## 7. Audit

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/repeated_crabb_joint_defect_comparison.py \
  --output \
  experiments/repeated_crabb_joint_defect_comparison_s70224.jsonl
```

The exact checker verifies both Young bounds behind (12) for
deterministic rational vectors and coordinate subspaces, plus the
small-neighbourhood absorption in (5).  It audits the frame-free
distance algebra in Section 4.  The analytic observation comparison
(10) is supplied by L197's exact finite Schur graph and L324's
accumulated retraction, not by the finite computation.
