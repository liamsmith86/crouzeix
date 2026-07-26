# The full delayed Schur face collapses to the two-orbit anticommutator

## 1. Result (L283, 2026-07-25)

Retain the balanced pure partial isometry

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F
\]

and impose the complete delay

\[
B_1=\cdots=B_{k-1}=0,\qquad
B_j=W^*(S^*)^jV.                                  \tag{1}
\]

Let \(\mathcal L_S(c)\) be L247's full dual, left-defect Schur
residual for the physical boundary metric, and let
\(\mathcal K_S(c)\) be the associated initial/right-defect residual
from L225--L228.  Put

\[
E_1=S^*ES,\qquad
F_r=S^rF(S^*)^r.
\]

Then the full dual face is exactly L278's intact two-orbit frontier:

\[
\boxed{
[c^{2k}]\mathcal L_S
=(ES^kF)(S^*)^k+S^kF(S^*)^kE.}                   \tag{2}
\]

Every lower coefficient vanishes.  L246's associated-defect
conjugacy therefore gives

\[
\boxed{
[c^{2k}]\mathcal K_S
=E_1F_{k-1}+F_{k-1}E_1.}                          \tag{3}
\]

Thus L228's former all-grade candidate is proved.  Its exact Stein
endpoint is

\[
\boxed{
W^*\mathcal G_S([c^{2k}]\mathcal K_S)W
=2B_kB_k^*.}                                      \tag{4}
\]

Consequently L227's canonical exact repair has, at the first active
grade of a **completely delayed stratum**, the physical metric faces

\[
\boxed{
\text{lower}=+B_k^*B_k,\qquad
\text{upper}=-12B_kB_k^*.}                        \tag{5}
\]

The repair is analytic and uses no transfer inverse, flag projection,
or pseudoinverse.  Formula (5) closes the all-grade
**complete-delay** covariance that was open in L228.  It does not by
itself close bounded selection through partial rank changes: A172
proves that when earlier \(B_j\) are nonzero but singular, mixed odd
coefficients of this unmodified repair can be indefinite on the next
flag.  L230--L234 prepare those mixed faces through grade three; the
arbitrary-grade right-ideal recursion in A178 remains open.  L283 is
not a proof of the repeated-elliptic neighbourhood or of the full
Crouzeix conjecture.

## 2. The edge-deleted residual is the closed-return defect

Let \(X_k=[c^{2k}]R\) be the active boundary-metric coefficient.
Use a circle to denote deletion of \(c^{2k}X_k\).

L258 proves that the edge-deleted closed-return defect
\(\mathfrak D^\circ\) is similar to the normalized final-defect Schur
residual:

\[
K_F^\circ
=(R_P^\circ)^{-1/2}\mathfrak D^\circ
 (R_P^\circ)^{1/2}.                               \tag{6}
\]

L277--L278 prove

\[
\mathfrak D^\circ=O(c^{2k}).
\]

The similarity factors in (6) have constant coefficient \(I\), so
they cannot change the first nonzero coefficient:

\[
[c^{2k}]K_F^\circ=[c^{2k}]\mathfrak D^\circ.       \tag{7}
\]

The deleted dual slack is the metric congruence of the same normalized
final defect.  L246's first-face congruence principle therefore gives

\[
\boxed{
[c^{2k}]\mathcal L_S^\circ
=[c^{2k}]\mathfrak D^\circ.}                      \tag{8}
\]

No determinant or trace is used in (6)--(8).  This is the operator
information that was not needed for L279's scalar volume conclusion.

## 3. Reinserting the active metric cancels the deletion response

Put \(P=I-F\).  L279 proves the exact closed-return response

\[
\boxed{
[c^{2k}](\mathfrak D^\circ-\mathfrak D)
=P(X_k-SX_kS^*)P.}                                \tag{9}
\]

L247 independently proves that reinserting the same metric
coefficient into the dual Schur residual contributes the opposite
operator:

\[
\boxed{
[c^{2k}](\mathcal L_S-\mathcal L_S^\circ)
=P(-X_k+SX_kS^*)P.}                               \tag{10}
\]

Add (8)--(10).  The response cancels before taking a trace:

\[
[c^{2k}]\mathcal L_S=[c^{2k}]\mathfrak D.          \tag{11}
\]

Finally L278 identifies the intact face exactly as

\[
[c^{2k}]\mathfrak D
=(ES^kF)(S^*)^k+S^kF(S^*)^kE,
\]

which proves (2).  L277--L278's lower vanishings, (6), and the fact
that the two metrics differ first in degree \(2k\) give the asserted
lower vanishings for \(\mathcal L_S\).

## 4. Associated-defect conjugacy

L246 gives, at the first nonzero face,

\[
[c^{2k}]\mathcal L_S
=S[c^{2k}]\mathcal K_SS^*.                        \tag{12}
\]

Since \(S\) is unitary from its initial space \(I-E\) onto its final
space \(I-F\),

\[
[c^{2k}]\mathcal K_S
=S^*[c^{2k}]\mathcal L_SS.                        \tag{13}
\]

Apply (13) to the two terms in (2).  Complete delay gives

\[
F(S^*)^{k-1}E=0,\qquad ES^{k-1}F=0.
\]

Therefore

\[
\begin{aligned}
S^*(ES^kF(S^*)^k)S
&=(S^*ES)\,S^{k-1}F(S^*)^{k-1}
=E_1F_{k-1},\\
S^*(S^kF(S^*)^kE)S
&=S^{k-1}F(S^*)^{k-1}(S^*ES)
=F_{k-1}E_1.
\end{aligned}
\]

This proves (3).

## 5. Endpoint and canonical repair

For completeness, the endpoint calculation from L228 is short.
Write \(W_r=S^{k-1}W\) and
\(C_r=S^*VB_k^*\).  The delay makes

\[
E_1F_{k-1}=C_rW_r^*,\qquad
F_{k-1}E_1=W_rC_r^*.
\]

Only the Stein orbit \(n=k-1\) reaches the left endpoint, so each
term contributes \(B_kB_k^*\).  This proves (4).

L227's canonical repair subtracts the Stein inverse of the Schur
residual.  Since (3) acts on the initial complement, its repair has
zero lower face.  Equation (4), together with the physical upper
factor four, gives the repair contribution

\[
-8B_kB_k^*.
\]

L223/L225's boundary metric contributes

\[
+B_k^*B_k\quad\text{below},\qquad
-4B_kB_k^*\quad\text{above}.
\]

Adding proves (5).

The full repair

\[
X(c)=-\mathcal G_{\widehat T(c)}
\left(
H-HV(V^*HV)^{-1}V^*H
\right)
\]

is L227's analytic exact contraction repair.  Its pivot has constant
coefficient \(I\); hence it stays analytic through every rank change
of \(B_k\).  This analyticity must not be confused with the upper
metric bound.  On an exact complete delay, (5) gives the required
first face.  Along a partial flag, A172's mixed cubic term survives
even though the preceding Gram vanishes on the compressed kernel, so
one cannot simply iterate (5).  Any all-grade continuation must
combine this exact delayed boundary value with the right-ideal mixed
preparations exemplified by L230--L234, or close L282's alternative
uniform energy estimate.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_delayed_face_collapse.py \
  --output \
  experiments/repeated_crabb_delayed_face_collapse_s70225.jsonl
```

The checker:

1. verifies the response-sign cancellation and (13) in exact delayed
   word algebra through grade twelve;
2. independently checks (4) on numerical complete-delay
   colligations through grade six; and
3. verifies the tracked pre-existing exact physical Schur audit
   through grade five as an independent finite falsification test.

Regenerate that slower audit independently with

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_delayed_slack_anticommutator.py
```

The identities above, not the finite audit, prove L283.
The tracked L283 dataset has SHA-256

```text
e45b478ae10270de59fc443ec076f8715d5dabe7654a26ff6f6df6de1842ed57
```
