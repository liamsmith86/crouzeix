# The second normalized gauge transports through the sextic flag

## 1. Result (L288, 2026-07-26)

Retain L286's coupled quartic and quintic columns

\[
\widetilde C_4=C_4-N_4,\qquad
\widetilde C_5=C_5-N_5+K_5,
\]

and L287's lower-normalized sextic forcing \(H_6^0\).  There is an
explicit polynomial correction \(K_6\) such that

\[
\boxed{
\widetilde C_6
=C_6-N_6+C_6^\parallel+K_6
\in\mathcal I_3,}                                 \tag{1}
\]

where \(\mathcal I_3\) is the triple-delay ideal.  Hence
\(\widetilde C_6=0\) whenever \(B_1=B_2=B_3=0\).

The transported sixth metric coefficient is not globally
endpoint-null; L287 proves that this is impossible.  Its complete
sixth upper-endpoint change instead factors through the first
transfer channel:

\[
\boxed{
\Delta{\cal U}_6
=-4W^*\{A_1(ESF)+(A_1ESF)^*\}W.}                 \tag{2}
\]

Consequently

\[
\boxed{
U^*\Delta{\cal U}_6U=0
\quad\text{for every }U\text{ with }B_1^*U=0.}    \tag{3}
\]

This is stronger than the condition needed at L234's surviving flag
\(\ker B_1^*\cap\ker B_2^*\).  Since L286 leaves every coefficient
through degree five unchanged, (3) preserves both L234's direct
sextic compression and its effective sextic face after the quintic
Schur cost.

L288 is a second exact finite recursive step.  It is not an
arbitrary-grade preparation recurrence and gives no summable
all-grade bound.

## 2. Exact transported equation

Write \(a=S^*\), \(s=S\), \(Q=I-E\).  L287 supplies

\[
EH_6^0E=0
\]

and the parallel factor correction

\[
{\cal C}_6^\parallel=-\frac12EH_6E
\in\mathcal I_3.
\]

There is a Hermitian polynomial \(Z_6^{\rm tr}\) such that, with

\[
\begin{aligned}
Y_6^{\rm tr}
&=Z_6^{\rm tr}-S^*Z_6^{\rm tr}S,\\
{\cal K}_6
&=Q(Y_6^{\rm tr}-H_6^0)E,
\end{aligned}                                     \tag{4}
\]

exact word reduction gives

\[
\boxed{
H_6^0+{\cal K}_6+{\cal K}_6^*
=Z_6^{\rm tr}-S^*Z_6^{\rm tr}S,\qquad
EZ_6^{\rm tr}E=0.}                               \tag{5}
\]

Thus the lower metric compression is still zero, while the upper
metric compression is allowed to carry L287's compulsory trace.
The exact sizes are

\[
\begin{array}{c|c|c}
\text{polynomial}&\text{reduced words}&\ell^1\\ \hline
Z_6^{\rm tr}&141&881/2\\
{\cal C}_6^\parallel&40&48\\
{\cal K}_6&116&261\\
\widetilde C_6V^*&142&1021/2.
\end{array}                                      \tag{6}
\]

The 72 Hermitian-orbit coefficients of \(Z_6^{\rm tr}\) are stored
verbatim in the checker.  Self-adjoint reduced words are inserted
once, not doubled; this convention is checked by the zero Hermitian
and forcing residuals.

## 3. Triple-delay divisibility

Exact quotient reduction proves

\[
\boxed{
{\cal C}_6^\parallel\in\mathcal I_3,\qquad
{\cal K}_6\in\mathcal I_3,\qquad
\widetilde C_6V^*\in\mathcal I_3.}                \tag{7}
\]

The first identity is L287.  The second is a new constraint built
into (4), and the third follows by adding it to L284's normalized
sextic lift and the parallel correction.

Direct coefficient counting gives the dimension-free, but currently
large, bound

\[
\boxed{\|\widetilde C_6\|\le1021/2.}              \tag{8}
\]

No claim is made that bounds of this size are summable in L194's
analytic chart.

## 4. Complete upper-endpoint factorization

The raw metric corner \(FZ_6^{\rm tr}F\) is not the whole order-six
endpoint change.  L286's degree-four metric change
\(\Delta M_4=-Z_4\) first enters the nonlinear upper Schur square at
degree six.

Let \(M_2\) be the unchanged degree-two prepared metric and

\[
J=I-F-\frac23E,
\]

the inverse of the constant upper-complement block.  One
minus-quarter of the complete endpoint change is

\[
\begin{aligned}
\Omega_6={}&FZ_6^{\rm tr}F
+FM_2J(-Z_4)F+F(-Z_4)JM_2F.                     \tag{9}
\end{aligned}
\]

Exact word reduction gives 40 words of coefficient
\(\ell^1\)-norm \(96\).  More importantly,

\[
\boxed{\Omega_6=A_1(ESF)+(A_1ESF)^*,}             \tag{10}
\]

where the ten-word factor of coefficient \(\ell^1\)-norm \(24\) is

\[
\begin{aligned}
A_1={}&
4a^5-4sa^6+2a^6s^3a^2+3a^5s^3a^3
+2a^4s^3a^4+a^3s^3a^5\\
&-2sa^7s^3a^2-3sa^6s^3a^3
-2sa^5s^3a^4-sa^4s^3a^5.                        \tag{11}
\end{aligned}
\]

The checker stores (11) as literal words to avoid ambiguity.
Because

\[
ESFWU=V(V^*SW)U=VB_1^*U,
\]

equations (2)--(3) follow immediately.  The quotient of
\(\Omega_6\) by even the first-delay ideal is zero.

This exact full-endpoint calculation is essential.  Factoring only
\(FZ_6^{\rm tr}F\) would omit the first nonlinear Schur contribution
and would not prove preservation of L234's effective face.

## 5. Scope guard: do not claim global two-endpoint invariance

Equation (5) makes the sixth **metric** lower compression zero.
It does not make the separately recomputed nonlinear lower Schur
coefficient globally unchanged: the degree-four metric motion first
enters that shorting at degree six.  The tracked unstructured and
rank-changing samples show nonzero sixth lower changes.

Therefore L288 proves the Stein-factor identity, triple-delay
normalization, and the complete upper-flag statement (2)--(3).  It
does **not** extend L286's phrase “complete lower and upper endpoints
are unchanged” to degree six.  Any final analytic metric must still
be expressed through L194's exact lower-tight chart; lower
feasibility may not be inferred from a raw coefficientwise
endpoint-null slogan.

The numerical checker records this lower movement rather than hiding
or thresholding it.  On a complete triple delay it vanishes, as all
three normalized columns do.

## 6. Consequence and next gate

The first two coupled transports now display a parity pattern:

1. the quintic successor has a globally endpoint-null witness;
2. the sextic successor necessarily carries an even endpoint trace,
   but its complete motion factors through an already removed
   transfer channel.

This is evidence for an odd-homology/even-flag-ideal recurrence, not
an induction.  The next valid advance is to derive that parity rule
in arbitrary grade and control the growth of its selected
right-ideal columns.  Computing an isolated seventh column would not
address the stated A178 debt.

## 7. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_second_gauge_transport.py \
  --output \
  experiments/repeated_crabb_second_gauge_transport_s70225.jsonl
```

The checker proves (5), (7), and (10) with zero rational residual
words.  It independently reconstructs the Stein witness, recomputes
the complete upper endpoint, and tests unstructured, rank-changing
multiplicities three through five, and triple-delay colligations.
An additional untracked stress audit passed 30 rank-changing and five
triple-delay cases.

The tracked dataset has SHA-256

```text
aaf15d56428af913f025be1d0132d256b72a892907f6f69bdba8881781a5dfa1
```
