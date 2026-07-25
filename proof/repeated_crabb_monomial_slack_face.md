# The delayed slack multiplier is one on every monomial channel

## 1. Result (L241, 2026-07-24)

Let \(S_k\) be the multiplicity-\(m\) backward shift of length \(k\)
on

\[
\mathbb C^{k+1}\otimes\mathbb C^m.
\]

Its transfer has the single nonzero coefficient

\[
B_k=I_m.
\]

Let \(R_{\rm bl}(c)\) be L219's balanced boundary metric and let
\(R_{\rm ax}(c)\) be the balanced form of L117's exact elliptic-axis
metric.  Put \(q=c^2\).  Then

\[
\boxed{
R_{\rm bl}(c)-R_{\rm ax}(c)
=q^k(3I-2E)+O(q^{k+1}).}                           \tag{1}
\]

All earlier coefficients vanish.  Since \(R_{\rm ax}\) has exact
rank-\(m\) Stein slack, (1) gives

\[
\boxed{
[q^k]\mathcal K_{S_k}
=2E_1,}                                           \tag{2}
\]

where \(E_1=S_k^*ES_k\) and \(\mathcal K\) is the right-defect
Schur residual of the boundary slack.

On this monomial channel,

\[
F_{k-1}=S_k^{k-1}F(S_k^*)^{k-1}=E_1.
\]

Therefore (2) is precisely

\[
\boxed{
[c^{2k}]\mathcal K_{S_k}
=E_1F_{k-1}+F_{k-1}E_1.}                          \tag{3}
\]

Thus L228's scalar multiplier is exactly one in **every** grade.  This
is an all-\(k\) calculation, not a finite-grade fit.

L241 does not yet prove L228 for a general matrix-inner transfer.
It fixes the only possible scalar normalization once the remaining
Hardy-cell locality statement is proved.  The live issue is now
whether the exact L239-reflected state lift can create additional
ordered terms involving a nonunitary/noncommuting active cell.

## 2. Boundary metric on the monomial path

Index the shift coordinates by \(0,\ldots,k\), with

\[
E=\Pi_0,\qquad F=\Pi_k.
\]

The right and left defect orbits in L219 are coordinate projections:

\[
E_n=\Pi_n,\qquad F_j=\Pi_{k-j}.
\]

After restoring the physical equality weights

\[
P=2I-E+2F
=\operatorname {diag}(1,2,\ldots,2,4),
\]

the boundary metric has diagonal entries

\[
\begin{aligned}
p_{{\rm bl},0}&=1+q^k,\\
p_{{\rm bl},n}
&=\frac2{1+q^n}+2q^{k-n}
 \quad(0<n<k),\\
p_{{\rm bl},k}&=\frac4{1+q^k}.                    \tag{4}
\end{aligned}
\]

## 3. First wrap of the exact axis metric

L117 writes the exact physical axis weights as

\[
p_{{\rm ax},n}
=\frac{S_n}{S_0c^n},\qquad
S_n=\sum_{\ell\in\mathbb Z}
\operatorname {sech}((n+2k\ell)(-\log c)).         \tag{5}
\]

Use

\[
\operatorname {sech}(a(-\log c))
=\frac{2c^{|a|}}{1+c^{2|a|}}\qquad(a\ne0).         \tag{6}
\]

For \(0<n<k\), the aliases \(n\), \(n-2k\), and \(n+2k\)
give, after division by \(c^n\),

\[
\frac2{1+q^n}+2q^{k-n}+2q^k+O(q^{k+1}).
\]

Meanwhile,

\[
S_0=1+4q^k+O(q^{2k}).
\]

Dividing by \(S_0\) therefore gives

\[
p_{{\rm ax},n}
=\frac2{1+q^n}+2q^{k-n}-6q^k+O(q^{k+1})
\quad(0<n<k).                                     \tag{7}
\]

At the two endpoints, the exact normalization and the two Nyquist
alias pairs give

\[
\begin{aligned}
p_{{\rm ax},0}&=1,\\
p_{{\rm ax},k}
&=\frac4{1+q^k}-12q^k+O(q^{k+1})
=4-16q^k+O(q^{k+1}).                              \tag{8}
\end{aligned}
\]

Subtracting (7)--(8) from (4) yields the physical first face

\[
[q^k](P_{\rm bl,phys}-P_{\rm ax})
=\operatorname {diag}(1,6,\ldots,6,12).           \tag{9}
\]

Conjugating by \(P^{-1/2}\) divides the first, interior, and last
entries by \(1,2,4\), respectively.  This proves (1).

## 4. Stein and Schur response

Both metrics use the same pulled operator

\[
A_k(c)=\phi_c\!\left(
S_k+c(I+F)S_k^*(I+E)
\right).
\]

Their difference starts at \(q^k\), and \(A_k(0)=S_k\).  Hence the
first difference of their balanced Stein slacks is

\[
\begin{aligned}
Y-S_k^*YS_k,\qquad Y=3I-2E.
\end{aligned}
\]

Since \(S_k^*S_k=I-E\),

\[
\begin{aligned}
Y-S_k^*YS_k
&=3I-2E-3(I-E)+2S_k^*ES_k\\
&=E+2E_1.                                         \tag{10}
\end{aligned}
\]

L117's exact axis slack has rank \(m\), so its right-defect Schur
residual is zero.  At \(c=0\), its pivot is \(E\) and its cross row
vanishes.  The first Schur response to (10) is consequently its
\(Q=I-E\) compression:

\[
Q(E+2E_1)Q=2E_1.                                  \tag{11}
\]

This proves (2).  Finally, the left defect at level \(k\), advanced
\(k-1\) times by the backward shift, is level one:

\[
F_{k-1}=E_1.
\]

Equation (3) follows.

## 5. Independent exact regeneration

Run

```bash
.venv/bin/python -u \
  experiments/repeated_crabb_monomial_slack_face.py \
  --output \
  experiments/repeated_crabb_monomial_slack_face_s70224.jsonl
```

The exact rational checker expands the periodized-\(\operatorname
{sech}\) weights, not a floating elliptic-function approximation.
For grades one through twelve it verifies:

1. every earlier metric difference is zero;
2. the physical face is \(1,6,\ldots,6,12\);
3. the balanced face is \(3I-2E\);
4. the Stein face is \(E+2E_1\); and
5. the Schur face equals the monomial anticommutator.

The tracked SHA-256 is
`7756334084e5c5fdd9786bd7222b8fd6a45e526fb22f55489dcb2b32fe602e2c`.
