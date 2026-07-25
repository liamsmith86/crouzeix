# The doubled terminal edge cannot be isolated from its boundary balance

## 1. Result (L249, 2026-07-25)

The factor two in L245's remote Green column is suggestive, but it is
not by itself a proof of L248's volume coefficient.  This can already
be seen exactly on the length-three monomial channel.

Let \(S\) be the backward shift on \(\mathbb C^3\), with initial and
final defects

\[
E=I-S^*S,\qquad F=I-SS^*.
\]

It has \(B_1=0\), so the physical delayed residual must vanish through
degree two.  Replace only the terminal factor in the balanced ellipse
pencil by a real multiplier:

\[
\Xi_\lambda(c)
=S+c(I+\lambda F)S^*(I+E),                        \tag{1}
\]

while keeping L219's boundary metric fixed.  The physical pencil is
\(\lambda=1\); its remote reverse edge has weight \(1+\lambda=2\).

Let \(K_{L,\lambda}\) be the final-defect dual Schur residual formed
from

\[
P_{\rm bl}^{-1}
-\phi_c(\Xi_\lambda)P_{\rm bl}^{-1}
 \phi_c(\Xi_\lambda)^*.
\]

In the natural two-dimensional final quotient,

\[
\boxed{
[c^2]K_{L,\lambda}
=
\begin{bmatrix}
1-\lambda^2&0\\
0&2\lambda-2
\end{bmatrix}.}                                   \tag{2}
\]

Consequently

\[
\boxed{[c^2]K_{L,\lambda}=0\quad\Longleftrightarrow\quad
\lambda=1.}                                       \tag{3}
\]

Thus the exact terminal doubling is uniquely locked to the physical
direct-map cancellation already one grade before the active \(B_2\)
face.  A proof of the \(+4\) flux may use the square of the doubled
Green amplitude only **after** placing it inside L248's lossless
final-row whitening.  Varying or estimating that port while freezing
the rest of the scalar network is false.

L249 does not disprove the coisometric-port strategy.  It rules out
the shortcut in which the endpoint amplitude is treated as an
independent perturbation and its square is declared to be the answer.
The correct proof must preserve the physical value \(\lambda=1\) and
use the exact theta/metric balance.

## 2. Exact calculation

L125's direct ellipse map through degree two is

\[
\phi_c(w)=(1+2c^2)w-cw^3+c^2w^5+O(c^3).
\]

Put

\[
J_\lambda=(I+\lambda F)S^*(I+E).
\]

Substitution of \(\Xi_\lambda=S+cJ_\lambda\) gives

\[
\begin{aligned}
A_0&=S,\\
A_1&=J_\lambda-S^3,\\
A_2&=2S-
 (J_\lambda S^2+SJ_\lambda S+S^2J_\lambda)+S^5.
\end{aligned}                                      \tag{4}
\]

The boundary metric has

\[
[c^2]P_{\rm bl}=SFS^*-S^*ES.                      \tag{5}
\]

For this length-three monomial shift, the matrix in (5) happens to be
zero.  Thus the exact example proves rigidity relative to the fixed
physical direct-map setup; it does not by itself isolate a
boundary-metric contribution on a generic tail.

Invert (5) through degree two, insert (4) in the dual Stein slack,
and Schur-compress away from \(F\).  Exact rational symbolic
reduction gives (2).  The second diagonal entry forces
\(\lambda=1\), and the first then vanishes automatically, proving
(3).

## 3. Adversarial implication for A194

The same frozen-multiplier experiment on generic noncommuting delayed
tails shows the same failure more violently: away from \(\lambda=1\),
coefficients below the intended \(c^{2k}\) face reappear.  Those
floating experiments are diagnostic only; the exact monomial
counterexample above is sufficient.

The safe next route is therefore:

1. keep the physical doubled chain and the complete lower metric;
2. insert L243--L245 directly into L248's whitened initial-defect
   mass;
3. use L244's exact coisometry to cancel the zero/one-reflection
   background; and only then
4. identify the surviving closed remote-channel energy.

## 4. Independent regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_terminal_multiplier_rigidity.py \
  --output \
  experiments/repeated_crabb_terminal_multiplier_rigidity_s70224.jsonl
```

The checker regenerates L125's scalar coefficients, constructs the
three exact operator coefficients in (4), inverts the metric, and
performs the final Schur quotient symbolically in \(\lambda\).  No
floating arithmetic is used.

The tracked dataset SHA-256 is

```text
114ba20c89b0d3e1e69db45ed510836a6306caca0cf1b02f0d89aa8f911ec642
```
