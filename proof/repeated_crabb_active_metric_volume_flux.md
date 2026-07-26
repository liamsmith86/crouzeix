# Active metric deletion doubles the intact delay flux

> **Operator closure note (2026-07-25).**  L283 combines this note's
> exact deletion response with L247's opposite dual-metric response
> and L258's first-face similarity.  The cancellation proves the full
> delayed Schur matrix face.  It closes the complete-delay covariance,
> not A172's mixed partial-flag selection obstruction.

## 1. Result (L279, 2026-07-25)

Retain the completely delayed pure partial isometry

\[
E=I-S^*S,\qquad F=I-SS^*,\qquad
ES^jF=0\quad(1\leq j<k),
\]

and write L219's physical metric as

\[
R(c)=I+\sum_{h\geq1}c^{2h}X_h.
\]

Let \(\mathfrak D_k\) be L258's closed-return defect with the physical
metric intact at degree \(2k\), and let
\(\mathfrak D_k^\circ\) use the edge-deleted metric

\[
R^\circ=R-c^{2k}X_k.
\]

L277--L278 prove every coefficient below degree \(2k\) vanishes.  At
the active face, the metric deletion has the exact operator response

\[
\boxed{
[c^{2k}](\mathfrak D_k^\circ-\mathfrak D_k)
=P(X_k-SX_kS^*)P,\qquad P=I-F.}                 \tag{1}
\]

Both the intact face and the deletion response have the same cyclic
representative:

\[
\boxed{
\begin{aligned}
[c^{2k}]\mathfrak D_k
&\overset{\rm tr}{\sim}2H_k,\\
P(X_k-SX_kS^*)P
&\overset{\rm tr}{\sim}2H_k,
\end{aligned}}                                   \tag{2}
\]

where

\[
H_k=Q_{k+2}-(k+2)Q_1+(k+1)I,\qquad
Q_j=(S^*)^jS^j.
\]

Consequently the full edge-deleted face satisfies the exact
all-grade cyclic identity sought in L252:

\[
\boxed{
[c^{2k}]\mathfrak D_k^\circ
\overset{\rm tr}{\sim}
4\{Q_{k+2}-(k+2)Q_1+(k+1)I\}.}                  \tag{3}
\]

Taking the trace and using L252 gives

\[
\boxed{
\operatorname {tr}[c^{2k}]\mathfrak D_k^\circ
=4\|B_k\|_F^2.}                                  \tag{4}
\]

Since every lower face is zero, L258 linearizes the logarithm and
proves L248's formerly open volume law:

\[
\boxed{
[c^{2k}]\log\mathcal V_{\widetilde A^{\lnot X_k}}
=4\|B_k\|_F^2.}                                  \tag{5}
\]

Thus A194's complete-delay volume identity is proved in arbitrary
grade.  Combining (5) with L247's exact metric split also proves the
full dual Schur trace

\[
\boxed{
\operatorname {tr}[c^{2k}]\mathcal L_S
=2\|B_k\|_F^2.}                                  \tag{6}
\]

By L246 duality and L225/L228's exact boundary retightening algebra,
this gives the all-grade effective separator trace

\[
\boxed{
\operatorname {tr}E_{2k,\mathrm{eff}}
=-16\|B_k\|_F^2.}                                \tag{6a}
\]

Thus L222's pointwise semidefinite range obstruction is closed on
every repeated flag.  L279 is not by itself a proof of the full
Crouzeix conjecture: a bounded analytic selection through commutant
rank changes and the global normal/elliptic merger still remain.

## 2. Linear response of the closed renewal

The intact and deleted metrics agree below degree \(2k\).  At the
active coefficient,

\[
\delta R=-X_k,\qquad
\delta R^{-1}=+X_k.                               \tag{7}
\]

Write

\[
Z=A R^{-1}A^*
\]

and use L258's final-row renewal

\[
Z_{\rm ret}
=PZP+PZF(I-FZF)^{-1}FZP.                         \tag{8}
\]

At \(c=0\),

\[
A(0)=S,\qquad Z(0)=SS^*=P,
\]

so the two cross blocks \(PZ(0)F\) and \(FZ(0)P\) vanish.  Equation
(7) therefore gives

\[
\delta Z=SX_kS^*,\qquad
\delta Z_{\rm ret}=PSX_kS^*P.                    \tag{9}
\]

The returned term in (8) has two zero cross factors at the
background, so its first derivative is zero; no final-row loop has
been discarded.

The retained metric changes by

\[
\delta R_P=-PX_kP.
\]

Since \(R_P(0)=Z_{\rm ret}(0)=P\), differentiating

\[
\mathfrak D=I_P-R_PZ_{\rm ret}
\]

gives

\[
\delta\mathfrak D
=PX_kP-PSX_kS^*P
=P(X_k-SX_kS^*)P,
\]

which proves (1).  Positive-degree motion of \(A\), the port, or the
defect graph cannot enter because the metric perturbation already has
degree \(2k\).  This is the closed-renewal counterpart of L247's
dual metric response, with the expected opposite deletion sign.

## 3. The intact face carries two units of flux

Put \(s=S\), \(a=S^*\), and

\[
D_k=Es^kF.
\]

L278 proves the exact intact face

\[
\boxed{
[c^{2k}]\mathfrak D_k
=D_ka^k+s^kD_k^*.}                               \tag{10}
\]

The second word is cyclically the first:

\[
s^kD_k^*
=s^kFa^kE
\overset{\rm tr}{\sim}
Es^kFa^k
=D_ka^k.
\]

Hence

\[
[c^{2k}]\mathfrak D_k
\overset{\rm tr}{\sim}2D_ka^k.                  \tag{11}
\]

The transfer coefficient is

\[
B_k=W^*a^kV,
\]

with \(E=VV^*\), \(F=WW^*\).  Therefore

\[
\operatorname {tr}(D_ka^k)
=\operatorname {tr}(B_k^*B_k)
=\|B_k\|_F^2.                                    \tag{12}
\]

This already proves that the intact face contributes two units of
the desired scalar energy.  The radial cyclic representative is
identified in Section 5.

## 4. The deletion response carries the same two units

L273's exact metric coefficient is

\[
\boxed{
X_k=s^kFa^k
+\sum_{d\mid k}(-1)^{k/d}a^dEs^d.}              \tag{13}
\]

Cyclicity and \(P^2=P\) give

\[
\begin{aligned}
P(X_k-sX_ks^*)P
&\overset{\rm tr}{\sim}
X_k\{P-s^*Ps\}\\
&=X_k(R_1-Q_1)
=X_k(E-F).                                      \tag{14}
\end{aligned}
\]

Consider the left-orbit term \(L_k=s^kFa^k\) in (13).  Its
\(-F\) product is cyclically zero because \(Fs^k=0\), while

\[
L_kE\overset{\rm tr}{\sim}Es^kFa^k=D_ka^k.       \tag{15}
\]

For a right-orbit term \(Y_d=a^dEs^d\), the \(E\) product is
cyclically zero because \(Ea^d=0\), and

\[
-(-1)^{k/d}Y_dF
\overset{\rm tr}{\sim}
-(-1)^{k/d}D_da^d.                               \tag{16}
\]

Every proper divisor \(d<k\) vanishes by the complete delay.
For \(d=k\), the coefficient is \((-1)^1=-1\), so (16) contributes
one further \(D_ka^k\).  Equations (14)--(16) prove

\[
P(X_k-SX_kS^*)P
\overset{\rm tr}{\sim}2D_ka^k.                  \tag{17}
\]

Thus the deletion response has the same cyclic class and the same
trace \(2\|B_k\|_F^2\) as the intact face.  Adding (11) and (17)
already proves (4).

## 5. Radial cyclic representative

It remains to identify \(D_ka^k\) with L252's three-term radial flux.
For \(1\leq j\leq k\),

\[
Q_j-Q_{j+1}=a^jEs^j.
\]

After a cyclic rotation, expand

\[
s^ja^j=s^{j-1}(I-F)a^{j-1}.
\]

The delay \(Es^{j-1}F=0\) (with \(EF=0\) for \(j=1\)) removes the
second term.  Iterating the same reduction through
\(Es^rF=0\), \(0\leq r<j\), leaves \(E\).  Hence

\[
Q_j-Q_{j+1}
\overset{\rm tr}{\sim}E=I-Q_1.                  \tag{18}
\]

Summing (18) gives

\[
\boxed{
Q_{k+1}
\overset{\rm tr}{\sim}
(k+1)Q_1-kI.}                                    \tag{19}
\]

Now use L278's five-word form:

\[
\begin{aligned}
2D_ka^k
&\overset{\rm tr}{\sim}
2I-2Q_1-2Q_{k+1}+2Q_{k+2}\\
&\overset{\rm tr}{\sim}
2\{Q_{k+2}-(k+2)Q_1+(k+1)I\}\\
&=2H_k.                                          \tag{20}
\end{aligned}
\]

Equations (11), (17), and (20) prove (2)--(3).  L252's radial trace
formula turns (3) into (4).

## 6. Volume and dual-Schur consequences

L277--L278 give

\[
\mathfrak D_k^\circ=O(c^{2k}).
\]

Therefore L258's exact identity

\[
\mathcal V=\det(I+\mathfrak D_k^\circ)
\]

has active logarithmic coefficient

\[
[c^{2k}]\log\mathcal V
=\operatorname {tr}[c^{2k}]\mathfrak D_k^\circ.
\]

Equation (4) proves (5), closing L248's open equation (7).

L247 proves

\[
\operatorname {tr}[c^{2k}]
\{\mathcal L_S-\mathcal L_S^{\lnot X_k}\}
=-2\|B_k\|_F^2,
\]

while L246/L248 identify the deleted first trace with (5).  Hence

\[
\operatorname {tr}[c^{2k}]\mathcal L_S
=4\|B_k\|_F^2-2\|B_k\|_F^2,
\]

which proves (6).

L246 transports this left-defect trace to L228's right-defect Schur
face without changing its first trace.  L228 Section 4 then combines
that value with the exact upper/lower boundary faces and the lower
retightening response, giving (6a).  L222's finite-dimensional
semidefinite alternative therefore has strict negative pairing on
every reducing delayed summand with \(B_k\ne0\).  The still-open issue
there is bounded analytic choice of a correcting column as the
commutant rank changes, not pointwise feasibility.

No large corner or Schur-square term has been estimated separately.
The cancellation is performed first by L258's closed renewal, and the
only comparison made here is the exact first response of that
whitened defect volume.

## 7. Independent exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_active_metric_volume_flux.py \
  --maximum-physical-grade 4 \
  --maximum-algebraic-grade 12 \
  --output \
  experiments/repeated_crabb_active_metric_volume_flux_s70225.jsonl
```

The checker independently:

1. constructs intact and edge-deleted L258 defects and verifies the
   operator response (1) through physical grade four;
2. verifies every lower coefficient is zero;
3. cyclically reduces the intact face, deletion response, and their
   sum to \(2H_k,2H_k,4H_k\), respectively; and
4. verifies the algebraic metric formula and cyclic identities
   through grade twelve.

The finite physical checks audit the closed-renewal differentiation.
The arbitrary-grade proof is the exact response calculation and
delay/cyclic algebra above.

The regenerated data file has SHA-256

```text
521e3e3f1129eaa36832738fb0ab86ffd73c766437bcbf75626522c9e28276b0
```
