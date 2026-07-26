# The delayed dual metric face contributes exactly \(-2\|B_k\|_F^2\)

> **Closure note (2026-07-25).**  L279 proves the complementary
> edge-deleted flux is \(+4\|B_k\|_F^2\), so the full delayed
> dual-Schur trace is \(+2\|B_k\|_F^2\) in every grade.

## 1. Result (L247, 2026-07-24)

Retain L225's balanced boundary metric and elliptic transform

\[
R(c)=P_{\rm bl}(c),\qquad A(c)=\widehat T_c,
\]

and form the dual Stein slack

\[
J(c)=R(c)^{-1}-A(c)R(c)^{-1}A(c)^*.              \tag{1}
\]

Let \(\mathcal L_S(c)\) be its Schur residual away from the left
defect \(F=WW^*\).  L246 shows that its first face is conjugate to
L225/L228's right-defect face.

Assume the complete delay

\[
B_1=\cdots=B_{k-1}=0,
\qquad B_j=W^*(S^*)^jV.
\]

Put

\[
X_k=[c^{2k}]R(c).
\]

Hold the complete direct-map series \(A(c)\) and all lower metric
coefficients fixed, but delete \(c^{2k}X_k\) from \(R(c)\).  Let
\(\mathcal L_S^{\lnot X_k}\) denote the resulting formal dual Schur
residual.  Then the metric-edge contribution is

\[
\boxed{
[c^{2k}]\{\mathcal L_S-\mathcal L_S^{\lnot X_k}\}
=(I-F)(-X_k+SX_kS^*)(I-F).}                       \tag{2}
\]

Its trace is exactly

\[
\boxed{
\operatorname {tr}[c^{2k}]
\{\mathcal L_S-\mathcal L_S^{\lnot X_k}\}
=-2\|B_k\|_F^2.}                                  \tag{3}
\]

Consequently the candidate all-grade trace law

\[
\operatorname {tr}[c^{2k}]\mathcal L_S
=2\|B_k\|_F^2
\]

is equivalent to the single operator/reflection flux identity

\[
\boxed{
\operatorname {tr}[c^{2k}]
\mathcal L_S^{\lnot X_k}
=4\|B_k\|_F^2.}                                   \tag{4}
\]

Equation (4) remains open.  L245's remote Green column has leading
amplitude \(2\), so its squared boundary energy has precisely the
normalization \(4\).  L247 separates that still-unproved
coisometric-flux assertion from the now-completed metric algebra.

## 2. Linear response of the inverse metric

Write

\[
R(c)=R_{<}(c)+c^{2k}X_k+O(c^{2k+1}).
\]

Coefficientwise inversion gives

\[
[c^{2k}]\{R^{-1}-(R^{\lnot X_k})^{-1}\}=-X_k,     \tag{5}
\]

because both metrics have constant coefficient \(I\); every lower
coefficient is the same.  Since \(A(0)=S\), substitution in (1)
gives

\[
[c^{2k}]\{J-J^{\lnot X_k}\}
=-X_k+SX_kS^*.                                    \tag{6}
\]

At \(c=0\), the dual slack is \(F\).  The graph injection defining
the left Schur quotient therefore has constant term \(I-F\).
The perturbation in (6) begins only at degree \(2k\), so every
positive-degree graph correction would move it beyond the requested
coefficient.  Taking the induced quotient form proves (2).

This argument does not assume that either full dual slack is positive.
It is a formal first-coefficient statement about two series that agree
below degree \(2k\).

## 3. Trace from the two endpoint faces

Put \(P=I-F\).  Cyclicity and the partial-isometry relations give

\[
\begin{aligned}
\operatorname {tr}P(-X_k+SX_kS^*)P
&=-\operatorname {tr}(PX_k)
  +\operatorname {tr}(S^*PSX_k)\\
&=-\operatorname {tr}\{(I-F)X_k\}
  +\operatorname {tr}\{(I-E)X_k\}\\
&=\operatorname {tr}\{(F-E)X_k\}.                 \tag{7}
\end{aligned}
\]

Here \(S^*PS=S^*(I-F)S=S^*S=I-E\), since \(FS=0\).

L223's exact delayed upper metric face says

\[
W^*X_kW=-B_kB_k^*,
\]

while L225's exact lower face says

\[
V^*X_kV=B_k^*B_k.
\]

Therefore

\[
\operatorname {tr}(FX_k)=-\|B_k\|_F^2,\qquad
\operatorname {tr}(EX_k)=+\|B_k\|_F^2.            \tag{8}
\]

Substituting (8) in (7) proves (3), and (4) follows by subtraction.

## 4. Why this sharpens the live gate

The conditional endpoint-word argument in A192 mixed two logically
different cancellations:

1. a completely explicit boundary-metric contribution; and
2. the state-lifted coisometric reflection flux.

L247 removes the first from the open problem.  The remaining theorem
has a fixed positive scalar target \(4\), no metric coefficient at the
active grade, and may be attacked at the remote left endpoint by
L245.  Later \(B_j^*\) terms in L243 still have to be shown to belong
to the internal coisometric energy balance; L247 does not assume that
cancellation.

## 5. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_dual_metric_trace_split.py \
  --output \
  experiments/repeated_crabb_dual_metric_trace_split_s70224.jsonl
```

The checker uses generic noncommuting complete delays through grade
eight.  It verifies both endpoint compressions, the matrix formula
(2), and the trace \(-2\|B_k\|_F^2\).  These calculations only use
L219's explicit metric coefficient and partial-isometry algebra; no
theta/ODE truncation is involved.  The tracked dataset SHA-256 is
`26cbe2bd4cd09b832e917eb7446374597e80a5796bd6dea7245ac7963e65d625`.
