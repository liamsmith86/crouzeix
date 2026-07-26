# The intact-metric even frontier is exactly the next delay defect

## 1. Result (L278, 2026-07-25)

Let

\[
E=I-S^*S,\qquad F=I-SS^*,\qquad EF=0,
\]

and use L258's weighted closed-return defect

\[
\mathfrak D=I-R_PZ_{\rm ret}.
\]

Fix \(h\geq1\).  Impose only the lower delays

\[
ES^jF=0\qquad(1\leq j<h),                        \tag{1}
\]

and keep the physical retained-metric coefficient at degree \(2h\).
Equivalently, any active metric deletion is made strictly after this
face.  Put

\[
\Theta_h=[c^{2h}]\mathfrak D.
\]

Then, modulo (1),

\[
\boxed{
\Theta_h
=(ES^hF)(S^*)^h
+S^h\{F(S^*)^hE\}.}                              \tag{2}
\]

Therefore the next complete delay \(ES^hF=0\) and its adjoint kill
the whole face:

\[
\boxed{\Theta_h=0.}                              \tag{3}
\]

Together with L277's odd delay-ideal factorization, this proves every
lower coefficient of the active grade-\(k\) closed defect vanishes:

\[
\boxed{
[c^d]\mathfrak D_k=0\qquad(0\leq d<2k).}         \tag{4}
\]

Thus the entire lower-face premise of L276 is closed in arbitrary
grade, more strongly than the required two-sided radiality.  The sole
remaining A223 support gate is now:

\[
\boxed{\text{active associated cyclic radiality through }Q_{k+2}.}
\]

L278 does not prove that active statement and does not yet prove the
edge-deleted \(+4\|B_k\|_F^2\) volume coefficient.

## 2. Five-word frontier

Write \(s=S\), \(a=S^*\), \(Q_1=as\), and
\(R_j=s^ja^j\).  Before imposing the next delay, (2) has the universal
five-word normal form

\[
\boxed{
\begin{aligned}
\Theta_h={}&
2I-2Q_1-2R_{h+1}\\
&+as^{h+2}a^{h+1}
+s^{h+1}a^{h+2}s .
\end{aligned}}                                   \tag{5}
\]

This is an arbitrary-grade first-visit induction, not a finite-grade
interpolation.

Let \(P=I-F\), \(T=SP\), and form the independently balanced deflated
tail.  After L258 has resummed every final-row visit, first-visit
grouping gives the intact-metric covariance

\[
\boxed{
\Theta_h(S)=\iota\{\Theta_{h-1}(T)\}
\qquad(h\geq2),}                                 \tag{6}
\]

where \(\iota\) embeds the tail word into the original state space.

There are two cases in this grouping.

### No new outer crossing

Removing the clean entry/exit layer sends a full closed path of
degree \(2h\) bijectively to its tail path of degree \(2h-2\).  L251
preserves the common analytic channel and L245 preserves both exterior
Green columns, so its complete coefficient is unchanged.

### A new outer crossing

At a lower face the physical metric has **not** been edge-deleted.
The local new-crossing sector is therefore the constant-kernel
half-line sector of L243--L245 with its complete retained metric.
L244 proves that this normalized half-line block is exactly a
coisometry.  Consequently its output return and retained-metric
contribution cancel before coefficient extraction:

\[
\bigl[I-R_PZ_{\rm ret}\bigr]_{\rm new\ crossing}=0. \tag{7}
\]

Any nonconstant model-kernel selection belongs to the retained tail
between the marked first and last visits and is carried by the
no-new-local-crossing bijection.  A second independent reflected
selection lies above the first-face valuation by L243.  Thus (7)
exhausts the local outer layer.

For completeness, the assertion that the surviving tail perturbation
is transported with coefficient one is exactly L263's first-face
Schur-graph derivative.  By induction, the tail defect has no earlier
coefficient.  L245's clean entry/exit columns add degree two, and at
the coisometric background L263's quotient map is the constant tail
embedding.  Motion of the defect graph or port is one order too late.
Hence the first tail face at degree \(2h-2\) becomes its embedded full
face at degree \(2h\), with no additional term.  This proves (6).

The cancellation is false after deleting the active metric
coefficient; that deliberate failure is precisely the edge-deleted
volume flux sought at degree \(2k\).

The grade-one physical metric is still intact here.  Substitution of
L256's exact direct-map jet together with the retained \(c^2\) metric
coefficient gives

\[
\Theta_1
=2I-2as-2s^2a^2+as^3a^2+s^2a^3s,               \tag{8}
\]

which is (5) at \(h=1\).  Now use

\[
T^j=S^{j+1}S^*,\qquad
(T^*)^j=S(S^*)^{j+1}
\]

in the five tail words.  Direct partial-isometry reduction sends
(5) at \(h-1\) exactly to (5) at \(h\).  Equations (6) and (8)
therefore prove (5) in every grade.

This is only a first-face associated identity under the delay
filtration.  It does not assert the false whole-series tail equality
ruled out by L235.

## 3. Exact delay-ideal factorization

Define

\[
\begin{aligned}
D_h&=Es^hF
=as^{h+2}a+s^h-s^{h+1}a-as^{h+1},\\
D_h^*&=Fa^hE
=sa^{h+2}s+a^h-a^{h+1}s-sa^{h+1}.
\end{aligned}                                    \tag{9}
\]

The lower delay relations imply

\[
\boxed{
as^{h+1}a^h+s^ha^{h+1}s
=2R_h-2I+2Q_1.}                                  \tag{10}
\]

For \(h=1\), this is the sum of \(EF=0\) and \(FE=0\) after expanding
the two defect projections.  If it holds at \(h-1\), multiply
\(D_{h-1}=0\) on the right by \(a^{h-1}\), multiply its adjoint on
the left by \(s^{h-1}\), and add.  The result is

\[
\begin{aligned}
as^{h+1}a^h+s^ha^{h+1}s
={}&as^ha^{h-1}+s^{h-1}a^hs\\
&+2R_h-2R_{h-1},
\end{aligned}
\]

which proves (10) by induction.

Expand the right side of (2):

\[
\begin{aligned}
D_ha^h+s^hD_h^*
={}&as^{h+2}a^{h+1}+s^{h+1}a^{h+2}s\\
&+2R_h-2R_{h+1}\\
&-\{as^{h+1}a^h+s^ha^{h+1}s\}.
\end{aligned}
\]

Equation (10) turns this into exactly (5), proving the factorization
(2).  Imposing \(D_h=D_h^*=0\) proves (3).

## 4. All lower vanishings

In an active grade-\(k\) edge-deleted system, only the metric
coefficient at degree \(2k\) is removed.  Every lower even coefficient
\(2h<2k\) therefore retains its physical metric and satisfies L278.
The complete delay includes \(ES^hF=0\), so (3) kills it.  Every lower
odd coefficient vanishes by L277, and L258 gives the zero constant
face.  This proves (4).

In particular, no long-axis independence argument is now needed to
turn a merely radial lower face into zero: the operator coefficient
itself vanishes in the delay quotient.

## 5. Independent exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_even_frontier_delay_ideal.py \
  --maximum-physical-delay 4 \
  --maximum-factor-delay 12 \
  --output \
  experiments/repeated_crabb_even_frontier_delay_ideal_s70225.jsonl
```

The checker independently:

1. constructs the full intact-metric L258 defect and verifies the
   five-word formula through delay four;
2. constructs the independently balanced tail and verifies that the
   full-minus-tail associated layer is exactly zero;
3. verifies the factorization (2) and tail embedding (6) through
   delay twelve; and
4. imposes the next delay and obtains the zero polynomial in every
   record.

The finite physical checks audit the implementation; the arbitrary-
grade proof is the first-visit/coisometry recursion and exact
delay-ideal factorization above.  The tracked dataset has SHA-256
`69a362b89743d7e92416cd6aef56796e08b25a34aac883c94f1ef0428f488a3f`.
