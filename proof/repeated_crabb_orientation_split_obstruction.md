# The two ellipse orientations do not cancel separately

## 1. Result (L260, 2026-07-25)

There is a natural two-variable refinement of the real ellipse
series.  Write

\[
[w^{2n+1}]\phi_c(w)
=\sum_{\ell\geq0}\alpha_{n,\ell}c^{n+2\ell}.
\]

Replace one scalar monomial by

\[
c^{n+2\ell}\longmapsto a^{n+\ell}b^\ell,
\]

replace the reverse pencil edge \(cJ\) by \(aJ\), let formal adjoint
swap \(a,b\), and replace the metric weight \(q^h=c^{2h}\) by
\(a^hb^h\).  Setting \(a=b=c\) recovers the physical operator,
metric, and volume exactly.

The tempting claim that the bidegree sectors cancel separately is
false already at grade one.  In the exact cyclic partial-isometry
quotient, both total-degree-one sectors are nonzero:

\[
\begin{aligned}
[a^0b^1]\mu^\circ&=-\operatorname {tr}S^2
 \operatorname {tr}(S^*)^2,\\
[a^1b^0]\mu^\circ&=+\operatorname {tr}S^2
 -\operatorname {tr}(S^*)^2.
\end{aligned}
\]

Only their sum vanishes on the physical diagonal.  At total degree
two, all three sectors \((0,2),(1,1),(2,0)\) are again nonzero.
Only their sum is

\[
\boxed{
4\operatorname {tr}\{(S^*)^3S^3-3S^*S+2I\}.}     \tag{1}
\]

Thus neither lower vanishings nor the active leakage trace may be
proved by assigning independent positivity or cancellation to the
two ellipse orientations.  The full theta/ODE sum must be assembled
before L244's coisometric cancellation or L258's closed-return
quotient is taken.

This is a negative structural result, not progress on A194's
remaining coefficient.  It rules out one seemingly clean way to
justify L259's diagonal selection.

## 2. Exact calculation

At grade one, lift every term of the direct map by the rule above,
form

\[
N=R^\circ-A^*FA,\qquad
H=R^\circ-A^*R^\circ A,\qquad
\mu^\circ=-m+\operatorname {tr}(N^{-1}H),
\]

and retain total bidegree at most two.  Exact word reduction gives

\[
\begin{array}{c|l}
(0,1)&-\operatorname {tr}S^2+\operatorname {tr}(S^*)^2\\
(1,0)&+\operatorname {tr}S^2-\operatorname {tr}(S^*)^2
\end{array}
\]

at total degree one.  Their physical sum is zero.

At total degree two the three separate cyclic representatives are

\[
\begin{aligned}
(0,2):\quad&
4\operatorname {tr}S^*S
-\operatorname {tr}(S^*)^3S^3
-\operatorname {tr}(S^*)^4,\\
(1,1):\quad&
-20\operatorname {tr}S^*S+8\operatorname {tr}I
+6\operatorname {tr}(S^*)^3S^3\\
&+\operatorname {tr}S^4+\operatorname {tr}(S^*)^4,\\
(2,0):\quad&
4\operatorname {tr}S^*S
-\operatorname {tr}(S^*)^3S^3
-\operatorname {tr}S^4.
\end{aligned}
\]

The pure \(S^4,(S^*)^4\) terms and the surplus radial terms cancel
only after summing all three rows, leaving (1).

## 3. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_orientation_split_obstruction.py \
  --output \
  experiments/repeated_crabb_orientation_split_obstruction_s70225.jsonl
```

The checker also verifies that setting \(a=b=c\) regenerates the
ordinary direct operator and edge-deleted metric before testing the
cyclic obstruction.

The tracked dataset regenerates with SHA-256

```text
fe19fdc52b50dd1f61a2a1239917673b69c76f5be77948623a4445cc68561cbb
```
