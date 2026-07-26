# The four remote return coefficients are fixed in every delayed grade

## 1. Result (L274, 2026-07-25)

Let

\[
 E=I-S^*S,\qquad F=I-SS^*,\qquad EF=0
\]

and assume the complete delay

\[
 ES^jF=0\qquad(1\leq j\leq k-1),\qquad k\geq2.
\]

Use L258's edge-deleted output Gram

\[
 Z=A(R_k^\circ)^{-1}A^*
\]

and its unweighted final-row renewal

\[
 \mathfrak U_k
 =I-\{PZP+PZF(I-FZF)^{-1}FZP\},\qquad P=I-F.
\]

Let \(\mathfrak U_{k-1}^{\rm tail}\) be the independently balanced
renewal of the deflated tail, embedded in the original word algebra.
Put

\[
 \Delta_k^{\mathfrak U}
 =[c^{2k}]\mathfrak U_k
  -[c^{2k-2}]\mathfrak U_{k-1}^{\rm tail}.       \tag{1}
\]

Define

\[
\begin{aligned}
C_k&=(S^*)^kS^{k+2}(S^*)^2,\\
\widetilde C_k&=S^2(S^*)^{k+2}S^k,\\
G_k&=S(S^*)^{k+2}S^{k+2}S^* .
\end{aligned}                                    \tag{2}
\]

Then the four remote coefficients in (1) are exactly

\[
\boxed{
[C_k]\Delta_k^{\mathfrak U}
=[\widetilde C_k]\Delta_k^{\mathfrak U}
=[G_k]\Delta_k^{\mathfrak U}=1,\qquad
[SG_kS^*]\Delta_k^{\mathfrak U}=-1.}             \tag{3}
\]

The two crossed coefficients have a load-bearing cancellation:

\[
\boxed{
\begin{array}{c|ccc}
&-PZP&-PZFZP&
-PZF\{(I-FZF)^{-1}-I\}FZP\\ \hline
C_k,\widetilde C_k&
1+4(-1)^k&-4(-1)^k&0.
\end{array}}                                     \tag{4}
\]

Thus neither the direct Gram nor the immediate final-row return has
the desired coefficient separately.  Their physical sum is one.
Every return with at least one intervening final-row loop misses both
crossed words.

The \(G_k\) coefficient comes only from the full direct Gram.  The
\(-SG_kS^*\) coefficient comes only from subtracting the deflated
tail direct Gram.  No final-row return contributes to either.

L274 proves the four nonradial coefficients exposed by the
full-minus-tail unweighted audits and feeding the L269--L271
telescopes.  It does **not** prove that these are the only nonradial
words, does not determine the radial remainder, and does not prove
the lower closed-return vanishings.  L276 subsequently shows that
literal exclusion of every additional word is unnecessary: bounded
active cyclic radiality and lower radiality suffice.  The four
endpoint normalizations need not be computed again.
L277--L278 subsequently close every lower face; active cyclic
radiality remains open.

## 2. Direct-map and metric inputs

Write \(s=S\), \(a=S^*\).  The first delay gives

\[
FS^*E=0,
\]

so the balanced reverse edge is

\[
(I+F)S^*(I+E)=3a-a^2s-sa^2.                     \tag{5}
\]

L125/L242 gives the Newton-edge filtration

\[
[w^{2n+1}]\phi_c(w)
=(-1)^nc^n+O(c^{n+2}),                           \tag{6}
\]

including the exact exceptional coefficient at the next diagonal.
L273 gives, below the unavailable delay, the top inverse-metric
boundary

\[
[c^{2h}]R^{-1}
=\hbox{lower radial words}-Q_{h+1}+R_{h+1}.      \tag{7}
\]

Only (5)--(7), the partial-isometry reductions

\[
sas=s,\qquad asa=a,
\]

and the available delay relations enter the four coefficient
extractions below.  Later theta diagonals cannot reach the displayed
outer runs: each later scalar diagonal spends four additional word
degrees before the required first unavailable crossing.

## 3. The two crossed words

Expand

\[
Z=\sum_{i+h+j=2k}A_iY_hA_j^*,\qquad
Y_h=[c^h](R_k^\circ)^{-1}.                       \tag{8}
\]

For \(k\geq3\), reading the three alternating runs of \(C_k\) from
the outside inward leaves the following exhaustive source table in
\(PZP\):

\[
\begin{array}{c|c}
(i,h,j)&[C_k]\{PA_iY_hA_j^*P\}\\ \hline
(0,0,2k)&(-1)^k\\
(1,0,2k-1)&(-1)^{k+1}\\
(1,2k-4,3)&-1\\
(2k-3,0,3)&(-1)^{k+1}\\
(2k-2,0,2)&(-1)^{k+1}\\
(2k-2,2,0)&(-1)^{k+1}\\
(2k,0,0)&(-1)^{k+1}.
\end{array}                                      \tag{9}
\]

Indeed, once the central \(s^{k+2}\) run is selected, every further
reverse-edge placement either creates one of the seven rows in (9),
creates an available \(ES^jF\) reduction with \(j<k\), or changes an
outer run and cannot produce \(C_k\).  Equation (6) fixes the signs
in the zero-metric rows, while (7) fixes the two nonzero-metric rows.
The first two rows cancel.  The remaining five sum to

\[
[C_k]PZP=-1+4(-1)^{k+1}.
\]

The minus sign in \(-PZP\) gives the first entry of (4).  For \(k=2\)
the coincident rows in (9) combine to the five contributions

\[
1,\ -3,\ -1,\ -1,\ -1,
\]

whose sum is \(-5\), so the same formula holds.

In the immediate return \(PZFZP\), the final-row relations
\(FS=0\), \(S^*F=0\) force the unique degree split

\[
(2k-1,1).
\]

The two physical endpoint factors are both two, as retained inside
L251's common analytic channel, and the Newton-edge parity is
\((-1)^k\).  Therefore

\[
[C_k]PZFZP=4(-1)^k.                              \tag{10}
\]

An additional \(FZF\) loop inserts another terminal run.  It cannot
be removed without either using an available lower delay or changing
one of the three runs of \(C_k\).  Hence every positive-loop term has
zero \(C_k\) coefficient.  Equations (9)--(10) prove (4) for \(C_k\).
Taking adjoints proves the identical statement for
\(\widetilde C_k\).

Notice that (10) does not detach and square the terminal multiplier
ruled out by L249.  It is the coefficient of the complete physical
entry/exit product after the direct term in (9) has been retained.

## 4. The deep divergence pair

The same run count for \(G_k\) leaves only

\[
\begin{array}{c|c}
(i,h,j)&[G_k]\{PA_iY_hA_j^*P\}\\ \hline
(0,0,2k)&(-1)^{k+1}\\
(1,0,2k-1)&(-1)^k\\
(1,2k-2,1)&-1\\
(2k-1,0,1)&(-1)^k\\
(2k,0,0)&(-1)^{k+1}.
\end{array}                                      \tag{11}
\]

The first pair and last pair cancel, leaving

\[
[G_k]PZP=-1.
\]

Thus \(-PZP\) contributes \(+G_k\).  A final-row entry or exit would
add a forbidden outer run, so all returned terms have zero
\(G_k\) coefficient.

For the independently balanced tail, the corresponding table for
the embedded word \(SG_kS^*\) is

\[
\begin{array}{c|c}
(i,h,j)&[SG_kS^*]\{P_{\rm t}A_i^{\rm t}
Y_h^{\rm t}(A_j^{\rm t})^*P_{\rm t}\}\\ \hline
(0,0,2k-2)&(-1)^k\\
(1,0,2k-3)&(-1)^{k+1}\\
(1,2k-4,1)&-1\\
(2k-3,0,1)&(-1)^{k+1}\\
(2k-2,0,0)&(-1)^k.
\end{array}                                      \tag{12}
\]

Again the parity pairs cancel and the direct output-Gram coefficient
is \(-1\).  Since (1) subtracts the tail unweighted defect, its
contribution to \(\Delta_k^{\mathfrak U}\) is \(-SG_kS^*\).
For \(k=2\), the three merged entries \(1,-3,1\) again sum to
\(-1\).  Final-row support excludes every returned contribution.
Equations (11)--(12) prove the last two identities in (3).

## 5. Independent exact audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_remote_endpoint_coefficients.py \
  --maximum-grade 4 \
  --output \
  experiments/repeated_crabb_remote_endpoint_coefficients_s70225.jsonl
```

The checker independently constructs the full and deflated
edge-deleted output Grams and separates:

1. the direct retained Gram;
2. the immediate entry/exit return; and
3. every return containing at least one final-row loop.

In grades two through four it verifies (4), the two deep direct
coefficients, and the absence of every claimed returned contribution
using exact rational word arithmetic.  The tracked dataset has
SHA-256
`4d360bee00b01b1f4e658c8ac4a634d78f0ea519714cd392b11dc3d22176b8d2`.
