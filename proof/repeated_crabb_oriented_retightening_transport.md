# One orbit Gram and one-seventh of L212 retighten every transfer grade

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

> **Moving-series update.**  L299 derives the exact six-term
> nonlinear defect created by inserting this fixed-base correction
> into L227's moving operator/frame.  It begins cubically and is not
> automatically a prior-flag factor.  L298 remains the all-grade
> complete-delay boundary value; L299 is the active homology forcing.

## 1. Result (L298, 2026-07-26)

Retain the balanced pure partial-isometry data

\[
 I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad
 V^*W=0,
\]

and put

\[
\begin{aligned}
 Q&=I-E,\\
 P&=2I-E+2F,\\
 T&=P^{-1/2}SP^{1/2},\\
 B_j&=W^*(S^*)^jV,\\
 {\mathfrak C}(A)&=\sum_{n\ge1}B_nAB_n^* .
\end{aligned}                                      \tag{1}
\]

For every grade \(k\), set

\[
 A_k=B_k^*B_k,\qquad L_k=B_kB_k^*                 \tag{2}
\]

and let \(\widehat C_k\) be L212's balanced column.  Define

\[
\begin{aligned}
 D_k
 &=\frac17\widehat C_k\\
 &=-\frac12Q\left\{
 S^kWB_k+
 \sum_{j=1}^{k-1}(S^*)^{k-j}VB_k^*B_j
 \right\},                                        \tag{3}\\
 \Omega(A)
 &=\sum_{n\ge0}(S^*)^nVAV^*S^n,                   \tag{4}\\
 R_k
 &=-\Omega(A_k)
 +{\cal G}_S(VD_k^*+D_kV^*),                      \tag{5}\\
 F_k^{\rm rt}
 &=-\frac12VA_k+D_k.                              \tag{6}
\end{aligned}
\]

Here \({\cal G}_S(H)=\sum_{n\ge0}(S^*)^nHS^n\).
Then the metric direction \(X_k=P^{1/2}R_kP^{1/2}\)
and frame direction \(C_k^{\rm rt}=P^{1/2}F_k^{\rm rt}\)
obey

\[
\boxed{
X_k-T^*X_kT
=VC_k^{{\rm rt}*}+C_k^{\rm rt}V^*.}               \tag{7}
\]

Their physical endpoint pair is exactly

\[
\boxed{
V^*X_kV=-B_k^*B_k,\qquad
W^*X_kW=-4B_kB_k^*.}                              \tag{8}
\]

Thus L283's completely delayed raw canonical-repair face

\[
\left(+B_k^*B_k,\,-12B_kB_k^*\right)
\]

is changed in every grade to

\[
\boxed{\left(0,\,-16B_kB_k^*\right).}             \tag{9}
\]

This is an exact arbitrary-grade first-face construction.  It uses no
transfer inverse, flag projection, pseudoinverse, or isolated-grade
expansion.  It closes the complete-delay lower-tight endpoint debt
left after L214--L215.

By linearity, for every real summable sequence \(r_k\),

\[
\begin{aligned}
R(r)&=\sum_{k\ge1}r_kR_k,&
F^{\rm rt}(r)&=\sum_{k\ge1}r_kF_k^{\rm rt}
\end{aligned}
\]

obeys the same Stein/frame equation and has endpoint pair

\[
\boxed{
\left(
-\sum_{k\ge1}r_kB_k^*B_k,\,
-4\sum_{k\ge1}r_kB_kB_k^*
\right).}                                         \tag{10}
\]

For \(r_k=|c|^{2k}\), these are the two oriented reflected Grams.
Thus the all-series **base-colligation** transport is already closed.

It does **not** justify inserting all \(c^{2k}X_k\) simultaneously
into the full moving elliptic operator and metric.  That substitution
creates nonlinear cross terms absent from the fixed-\(S\) equation.
A179's failed raw superposition and A172's mixed partial-flag
obstruction therefore remain in force.  The live debt is now
sharper: compute the defect made by replacing \(S\) in (7) with the
full pulled operator, and prove that the defect is a bounded
prior-flag factor or can be absorbed by L289's finite graph
recurrence.

## 2. Stein/frame compatibility

The orbit Gram satisfies

\[
\Omega(A)-S^*\Omega(A)S=VAV^*.                    \tag{11}
\]

Therefore (5) gives

\[
\begin{aligned}
R_k-S^*R_kS
&=-VA_kV^*+VD_k^*+D_kV^*\\
&=V\left(-\frac12VA_k+D_k\right)^*
 +\left(-\frac12VA_k+D_k\right)V^*\\
&=VF_k^{{\rm rt}*}+F_k^{\rm rt}V^*.              \tag{12}
\end{aligned}
\]

Since \(P^{1/2}V=V\), conjugating (12) by \(P^{1/2}\)
proves (7).  Hence the retightening direction is compatible with the
rank-\(m\) Stein frame, not merely with the two endpoint
compressions.

At a first active coefficient \(c^{2k}\), adding
\(c^{2k}(X_k,C_k^{\rm rt})\) to L227's analytic raw repair changes
its Stein equation by exactly (7).  Products with positive-order
raw frame coefficients occur only above order \(2k\).  Thus (9) is a
valid first-face statement on every completely delayed stratum.

## 3. Lower endpoint

Because \(SV=0\), the \(n=0\) term is the only term of (4) surviving
between two copies of \(V\):

\[
V^*\Omega(A_k)V=A_k.                              \tag{13}
\]

Also \(V^*D_k=0\), and the same \(SV=0\) argument gives

\[
V^*{\cal G}_S(VD_k^*+D_kV^*)V=0.                 \tag{14}
\]

Equations (5), (13), and (14) prove the first identity in (8).

## 4. Upper endpoint and the cancellation

The orbit Gram carries its right-copy input through the complete
matrix-inner transfer:

\[
\begin{aligned}
W^*\Omega(A_k)W
&=\sum_{n\ge0}
 W^*(S^*)^nVA_kV^*S^nW\\
&={\mathfrak C}(A_k).                              \tag{15}
\end{aligned}
\]

The physical metric contributes a factor four at the left defect
because \(P^{1/2}W=2W\).  Hence the first term of (5) has upper
endpoint

\[
-4{\mathfrak C}(A_k).                             \tag{16}
\]

L212 proves

\[
{\cal M}_T(P^{1/2}\widehat C_k)
=28\{{\mathfrak C}(A_k)-L_k\}.                    \tag{17}
\]

Scaling by \(1/7\), the second term of (5) has physical upper
endpoint

\[
+4\{{\mathfrak C}(A_k)-L_k\}.                    \tag{18}
\]

The channel terms in (16)--(18) cancel exactly, leaving
\(-4L_k\).  This proves the second identity in (8), and adding it to
L283's raw pair proves (9).

This explains why the two large mechanisms must not be estimated
separately: the orbit transport gives
\(-4{\mathfrak C}(A_k)\), while the scaled L212 column returns the
same channel term with the opposite sign.

## 5. Apex and relation to L297

At a repeated monomial apex,

\[
B_j=0\ (j<k),\qquad B_k\ \hbox{unitary}.
\]

L212 gives \(D_k=0\), while the right Wold resolution gives
\(\Omega(I)=I\).  Consequently

\[
R_k=-I,\qquad X_k=-P.                             \tag{19}
\]

Thus L298 reduces exactly to scalar rescaling of the equality metric,
whose endpoint response is \((-I,-4I)\).  This is L297's
normalization ledger.  Away from the apex, (15)--(18) are the
noncommutative replacement for that scalar rescaling.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_oriented_retightening_transport.py \
  --output \
  experiments/repeated_crabb_oriented_retightening_transport_s70226.jsonl
```

The checker audits the Stein/frame equation, both endpoint
compressions, the two channel pieces separately, the raw-to-effective
face, and the scalar-rescaling apex.  An exact rational
noncommuting heterogeneous delay sum checks grades one through four,
including a nonzero lower-grade contamination term.  The 27
deterministic floating records include unstructured partial
isometries through grade six, completely delayed colligations through
grade six, and repeated monomial apices.  The proof is (11)--(18);
the audits are independent falsification checks.  The tracked dataset
has SHA-256

```text
ba5a8c971e377824280ad4afdc63b029d26f108534ec506a3328db4a4c8d1c7c
```
