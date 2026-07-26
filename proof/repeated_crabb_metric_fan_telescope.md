# The delayed metric inverse cancels every interior return fan

## 1. Result (L273, 2026-07-25)

Retain the pure partial isometry

\[
E=I-S^*S,\qquad F=I-SS^*,\qquad EF=0
\]

and assume complete delay through grade \(k-1\).  Put

\[
Q_j=(S^*)^jS^j,\qquad R_j=S^j(S^*)^j
\]

and write L219's metric as a \(q=c^2\) series

\[
R(q)=I+\sum_{h\ge1}q^hX_h.
\]

Delete its active coefficient:

\[
R^\circ_k(q)=R(q)-q^kX_k.
\]

Let \(P=I-F\).  Then the active coefficient of the retained inverse
has the exact support form

\[
\boxed{
P[q^k](R^\circ_k)^{-1}P
=V_k(Q_0,\ldots,Q_{k+1},R_1,\ldots,R_{k+1})
-\mathcal I_k,}                                  \tag{1}
\]

where \(V_k\) is radial and the interior fan is

\[
\boxed{
\begin{aligned}
\mathcal I_k={}&
\sum_{j=2}^{k}
(S^*)^jS^{k+2}(S^*)^{k+2-j}\\
&+\sum_{j=0}^{k-2}
S^{j+2}(S^*)^{k+2}S^{k-j}.
\end{aligned}}                                   \tag{2}
\]

Every word in (2) has coefficient exactly one.  Thus the active
metric inverse contributes coefficient \(-1\) to all \(2k-2\)
interior paths and contributes no other nonradial word.

Here is the exact interface with L258.  Write
\(\mathfrak D=I-R_PZ_{\rm ret}\) and
\(\mathfrak U=I-Z_{\rm ret}\).  If the required lower faces of
\(\mathfrak D\) vanish, then \(Z_{\rm ret}=R_P^{-1}\) through every
lower degree, and direct convolution gives

\[
[q^k](\mathfrak D-\mathfrak U)
=P[q^k](R_k^\circ)^{-1}P.                         \tag{1a}
\]

Thus (1) is precisely the metric half of L271's proposed endpoint
telescope, conditional only on the same lower vanishings that the
physical theorem already requires.  If the unweighted L258 return is
proved to contain the full fan

\[
\begin{aligned}
\mathcal F_k={}&
\sum_{j=0}^{k}
(S^*)^jS^{k+2}(S^*)^{k+2-j}\\
&+\sum_{j=0}^{k-1}
S^{j+2}(S^*)^{k+2}S^{k-j}
+G_k,
\end{aligned}                                    \tag{3}
\]

then (1) cancels every term of (3) except

\[
S^*S^{k+2}(S^*)^{k+1},\quad
S^{k+1}(S^*)^{k+2}S,\quad
R_{k+2},\quad G_k,                               \tag{4}
\]

which is exactly L271's four-word packet \(H_k\).

L273 does **not** prove that the physical unweighted return is (3),
does not prove the lower closed-return vanishings, and does not
identify the remaining radial polynomial.  L275 subsequently proves
coefficient one on every word of (3).  The sole remaining nonradial
question is now whether the physical unweighted return contains any
additional word.

## 2. Radial coefficients of the metric

The defect orbits telescope as

\[
(S^*)^dES^d=Q_d-Q_{d+1},\qquad
S^hF(S^*)^h=R_h-R_{h+1}.
\]

Consequently L219 gives

\[
\boxed{
X_h=R_h-R_{h+1}
+\sum_{d\mid h}(-1)^{h/d}(Q_d-Q_{d+1}).}         \tag{5}
\]

In particular \(X_h\) is radial, is supported at indices at most
\(h+1\), and its top boundary is

\[
\boxed{[X_h]_{\rm top}=Q_{h+1}-R_{h+1}.}         \tag{6}
\]

Let

\[
R(q)^{-1}=I+\sum_{j\ge1}q^jY_j.
\]

For \(j<k\), the complete-delay quotient gives

\[
\boxed{
Y_j\ \hbox{radial of index at most }j+1,\qquad
[Y_j]_{\rm top}=-Q_{j+1}+R_{j+1}.}               \tag{7}
\]

This follows inductively from

\[
Y_j=-X_j-\sum_{h=1}^{j-1}X_hY_{j-h}.             \tag{8}
\]

Indeed, same-side products of radial powers remain radial.  A mixed
product \(Q_aR_b\) or \(R_aQ_b\) is reduced by the delay relation
coming from \(ES^{a+b-2}F=0\) whenever \(a+b-2\le k-1\).
For the top parts of two lower coefficients one obtains the explicit
radial identity

\[
\boxed{
(Q_a-R_a)(-Q_b+R_b)
=Q_{\min(a,b)}+R_{\min(a,b)}-2I}                 \tag{9}
\]

whenever that delay relation is available.  Therefore every product
in the sum in (8) is radial and has index below \(j+1\); the top
boundary of \(Y_j\) comes only from \(-X_j\), proving (7).

For clarity, (9) is not a new commutation rule.  Assume \(a\le b\).
Same-side reduction gives
\(Q_aQ_b=Q_b\) and \(R_aR_b=R_b\).  Expanding
\(ES^{a+b-2}F=0\), multiplying on the left and right by the remaining
powers of \(S^*\), and applying the same three-letter reductions gives

\[
Q_aR_b+R_aQ_b
=Q_a+R_a+Q_b+R_b-2I.
\]

Substitution into the four-term product on the left of (9) proves the
display.  The case \(b<a\) is symmetric.  This also shows explicitly
why the required delay index is \(a+b-2\).

## 3. The first unavailable delay creates the fan

Since the active metric coefficient is deleted, inverse convolution
gives

\[
\boxed{
[q^k](R^\circ_k)^{-1}
=-\sum_{h=1}^{k-1}X_hY_{k-h}.}                   \tag{10}
\]

Fix one split \(1\le h\le k-1\).  Equations (6)--(7) show that the
only mixed top products are

\[
\begin{aligned}
Q_{h+1}R_{k-h+1}
&=(S^*)^{h+1}S^{k+2}(S^*)^{k-h+1},\\
R_{h+1}Q_{k-h+1}
&=S^{h+1}(S^*)^{k+2}S^{k-h+1}.                  \tag{11}
\end{aligned}
\]

Their indices sum to \(k+2\), so reducing them would require the
missing relation \(ES^kF=0\).  Every other mixed product has index
sum at most \(k+1\) and is radialized by the assumed delays.  The two
words in (11) both have coefficient \(+1\) in
\(X_hY_{k-h}\).  The minus sign in (10) gives coefficient \(-1\).

As \(h\) runs from \(1\) to \(k-1\), the first line of (11) is the
first sum in (2), with \(j=h+1\), and the second line is the second
sum, with \(j=h-1\).  This proves (1)--(2).

## 4. Exact audit

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_metric_fan_telescope.py \
  --maximum-grade 16 \
  --output \
  experiments/repeated_crabb_metric_fan_telescope_s70225.jsonl
```

For every grade the checker independently constructs L219's metric,
deletes the active coefficient, inverts the exact word series, and
verifies:

1. every lower inverse coefficient is radial with top boundary (7);
2. each split in (10) has exactly the two nonradial words (11); and
3. adding (2) to the active inverse leaves only radial words.

The script audits the arbitrary-grade proof rather than replacing it.
The tracked dataset regenerates with SHA-256
`bbf8afbe36d039f98a3c8c5df6dc257d34f4e354fe70e126265618ea89d7acff`.
