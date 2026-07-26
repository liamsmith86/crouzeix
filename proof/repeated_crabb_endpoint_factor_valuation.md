# Bounded endpoint factors are exactly Smith-valuation divisibility

> **Campaign scope.**  This note is used only in the fixed finite
> repeated-Crabb neighbourhood; see the L290--L316 scope guard in
> `LEMMA_LEDGER.md`.  It is not a global or dimension-uniform
> Crouzeix theorem.

## 1. Result (L292, 2026-07-26)

Let

\[
{\mathbb B}(s):\mathbb C^q\longrightarrow\mathbb C^m
\]

be a one-variable analytic matrix germ, and let \(G(s)=G(s)^*\) be
an analytic Hermitian endpoint germ.  Along the real parameter
\(s\), factorization of the type isolated by L291 is

\[
\boxed{
G(s)=X(s){\mathbb B}(s)^*
     +{\mathbb B}(s)X(s)^*.}                      \tag{1}
\]

The one-variable analytic Smith form gives invertible analytic
matrices \(P(s),Q(s)\) such that

\[
\boxed{
{\mathbb B}(s)=P(s)D(s)Q(s),\qquad
D(s)=
\begin{bmatrix}
\operatorname{diag}
(s^{\nu_1},\ldots,s^{\nu_\rho})&0\\
0&0
\end{bmatrix},}                                   \tag{2}
\]

where

\[
0\le\nu_1\le\cdots\le\nu_\rho
\]

and \(\rho\) is the generic rank.  Put

\[
\widehat G=P^{-1}GP^{-*}.                         \tag{3}
\]

Then (1) has a locally bounded solution \(X\) if and only if it has
an analytic solution, and this occurs exactly when the following
finite valuation conditions hold:

\[
\boxed{
\begin{aligned}
\widehat G_{ab}&\equiv0
 &&(a,b>\rho),\\
\operatorname{ord}_s\widehat G_{ab}
&\ge\min(\nu_a,\nu_b)
 &&(a,b\le\rho),\\
\operatorname{ord}_s\widehat G_{ak}
&\ge\nu_a
 &&(a\le\rho<k).
\end{aligned}}                                    \tag{4}
\]

Thus L291's pointwise flag condition is only the first line of (4).
It does **not** imply bounded factorization through a rank collapse.
The scalar example

\[
{\mathbb B}(s)=s^2,\qquad G(s)=s                 \tag{5}
\]

has a pointwise factor for every \(s\ne0\), and both matrices vanish
at the enlarged origin flag, but every factor is

\[
X(s)=\frac1{2s},
\]

which blows up.  A mixed active/kernel example with Smith exponents
\((1,3)\) similarly admits a pointwise cross factor of order
\(s^{-2}\).

L292 converts the remaining bounded-factor part of A178 into exact
order inequalities along the analytic arcs supplied by L197:

> after putting the transfer row and L289 endpoint residual into
> Smith coordinates, prove (4) at every prepared grade.

This is strictly sharper than proving compression zero, and strictly
weaker than forcing the full state witness into L290's hereditary
module.  It does not prove that the physical endpoint residuals
satisfy (4), nor does it prove the retained even face margins.

The Smith normal form and the criterion are standard one-variable
analytic linear algebra.  Their role here is to identify the exact
rank-changing divisibility debt and to provide explicit obstruction
tests before any further preparation is accepted.

## 2. Reduction to diagonal Smith coordinates

The ring of convergent one-variable complex power-series germs is a
discrete valuation ring.  Therefore ordinary Smith reduction gives
(2); unit factors are absorbed into \(P,Q\).  A real-analytic path
may be complexified locally and then restricted back to real \(s\).

Write

\[
Y=P^{-1}XQ^*.                                     \tag{6}
\]

Substitution of (2)--(3) in (1) gives the equivalent diagonal
problem

\[
\boxed{\widehat G=YD^*+DY^*.}                     \tag{7}
\]

Only the first \(\rho\) columns of \(Y\) occur.  If

\[
d_a=s^{\nu_a},
\]

then (7) says

\[
\begin{aligned}
\widehat G_{ab}
&=Y_{ab}d_b+d_a\overline{Y_{ba}}
 &&(a,b\le\rho),\\
\widehat G_{ak}
&=d_a\overline{Y_{ka}}
 &&(a\le\rho<k),\\
\widehat G_{k\ell}&=0
 &&(k,\ell>\rho).
\end{aligned}                                     \tag{8}
\]

If \(Y\) is merely locally bounded, every right side is
\(O(s^{\min(\nu_a,\nu_b)})\) or \(O(s^{\nu_a})\), as appropriate.
Because the corresponding entry of \(\widehat G\) is analytic, that
big-\(O\) bound is exactly the valuation lower bound in (4).  This
proves necessity even for a nonanalytic bounded factor.  Section 3
constructs an analytic factor from the same conditions.

## 3. Explicit analytic factor under the valuation conditions

Assume (4).  Define the first \(\rho\) columns of \(Y\) as follows.
On the diagonal, put

\[
Y_{aa}=\frac{\widehat G_{aa}}{2d_a}.              \tag{9}
\]

For \(a<b\le\rho\):

\[
\begin{cases}
Y_{ab}=0,\quad
Y_{ba}=\overline{\widehat G_{ab}/d_a},
&\nu_a\le\nu_b,\\[1mm]
Y_{ab}=\widehat G_{ab}/d_b,\quad
Y_{ba}=0,
&\nu_b<\nu_a.
\end{cases}                                      \tag{10}
\]

For \(a\le\rho<k\), put

\[
Y_{ka}=\overline{\widehat G_{ak}/d_a}.            \tag{11}
\]

Set every unused entry to zero.  The divisibilities in (4) make
(9)--(11) analytic.  Hermiticity of \(\widehat G\) and direct
substitution prove (8), hence (7).

Finally,

\[
\boxed{X=PYQ^{-*}}                                \tag{12}
\]

is analytic and proves (1).  Since \(P,Q\), and their inverses are
locally bounded, the factor \(X\) is locally bounded as well.

The construction is triangular and uses no singular-value
pseudoinverse.  Its divisions are only by the explicit Smith powers;
(4) certifies every division before it is performed.

## 4. Why pointwise flag zero is insufficient

Away from \(s=0\), the surviving endpoint flag is the last
\(m-\rho\) Smith coordinates.  L291's pointwise criterion is exactly

\[
\widehat G_{k\ell}=0
\qquad(k,\ell>\rho).                              \tag{13}
\]

When all \(\nu_a>0\), the transfer row may vanish entirely at the
origin.  Pointwise flag zero there adds only

\[
\widehat G(0)=0.                                  \tag{14}
\]

Conditions (13)--(14) do not control the relative vanishing rates.
For (5),

\[
\operatorname{ord}_sG=1<2
=\operatorname{ord}_s{\mathbb B},
\]

so (4) fails by one order and the factor grows like \(s^{-1}\).

For a mixed example, take

\[
D(s)=\operatorname{diag}(s,s^3,0)
\]

and let the only nonzero entries of \(G\) be

\[
G_{23}=G_{32}=s.
\]

The generic kernel block and \(G(0)\) are both zero, but (4) requires
\(\operatorname{ord}_sG_{23}\ge3\).  The unique exposed cross factor
has order \(s^{-2}\).

These examples are guardrails against treating L222/L291 pointwise
feasibility as uniform analytic selection.

## 5. Consequence for curve selection and A178

L197 already proves that a local failure of the repeated-block metric
bound yields a real-analytic arc.  Along such an arc, (2) exists and
L292 gives a finite test.

Therefore a sufficient endpoint-factor continuation is:

1. use L285/L194 to construct the finite admissible state jet;
2. use L289 to form the **complete** endpoint residual, not its raw
   metric corner;
3. put the joint earlier transfer row in Smith form;
4. prove all valuation inequalities (4) for that residual;
5. use (9)--(12) to obtain a bounded factor; and
6. prove the retained even Gram minus prior Schur cost has a uniform
   positive margin.

If bounded selection failed, curve selection would produce an arc on
which at least one inequality in (4) fails.  Consequently proving
(4) on every physical analytic arc rules out a hidden
rank-dependent blow-up.

L288's explicit polynomial first-channel factor automatically
satisfies (4) at the banked sextic step.  L292 is not an induction
from that example.  The next genuine theorem must derive the
valuation inequalities from the arbitrary-grade physical recurrence
or from L282's equivalent energy control; another isolated
coefficient calculation is not enough.

## 6. Regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_endpoint_factor_valuation.py \
  --output \
  experiments/repeated_crabb_endpoint_factor_valuation_s70226.jsonl
```

The exact SymPy checker covers:

1. an analytic active/kernel factor with Smith exponents \((1,3)\);
2. an analytic full-rank factor with exponents \((1,2,4)\);
3. the scalar pointwise-only obstruction (5); and
4. the mixed cross obstruction above.

Every exact factorization residual is zero.  At \(s=10^{-3}\), the
two valid factors remain bounded by \(3\) and \(4\), while the false
pointwise-only factors reach \(500\) and \(10^6\).  The tracked
dataset has SHA-256

```text
a9c28eec243b7854bb7f06959f359a1e9bcb4a40894e6a65aeee4c5eab80136d
```
