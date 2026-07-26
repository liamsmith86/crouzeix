# Candidate: every endpoint word is dominated by its transfer prefix

## 1. Conjecture (A259, 2026-07-26)

Let

\[
I-S^*S=VV^*,\qquad I-SS^*=WW^*,\qquad V^*W=0
\]

be a balanced pure partial-isometry colligation, and put

\[
B_j=W^*(S^*)^jV\qquad(j\ge1).
\]

For a word \(w\) of length \(d\) in \(S,S^*\), define

\[
K_w=W^*w(S,S^*)V.
\]

The following inequality is strongly supported numerically but is
**not yet proved**:

\[
\boxed{K_wK_w^*\preceq\sum_{j=1}^dB_jB_j^*.}     \tag{1}
\]

Equivalently, if \({\mathbb B}_d=[B_1\ \cdots\ B_d]\), the Douglas
lemma would give a contraction \(C_w\) such that

\[
\boxed{K_w={\mathbb B}_dC_w,\qquad\|C_w\|\le1.}   \tag{2}
\]

In particular,

\[
B_1=\cdots=B_d=0\Longrightarrow K_w=0,           \tag{3}
\]

which is consistent with L209.  A proof would add the missing
constant-one quantitative estimate for the multipliers
\(K=W^*qpV\) in L310--L311.

## 2. Exact equivalent energy formulation

Defect telescoping gives

\[
\sum_{j=0}^d(S^*)^jVV^*S^j
=I-(S^*)^{d+1}S^{d+1}.
\]

The \(j=0\) term vanishes after compression by \(W\), because
\(V^*W=0\).  Therefore

\[
\boxed{
\sum_{j=1}^dB_jB_j^*
=I-W^*(S^*)^{d+1}S^{d+1}W.}                     \tag{4}
\]

Consequently (1) is exactly

\[
\boxed{
\|V^*w(S,S^*)^*Wy\|^2+\|S^{d+1}Wy\|^2
\le\|y\|^2\qquad(y\in\mathbb C^m).}              \tag{5}
\]

Thus the conjecture says that an arbitrary mixed forward/backward
endpoint observation of length \(d\), together with the untouched
all-forward residual at depth \(d+1\), cannot carry more energy than
the initial port input.  This is the precise lossless-system
statement that must be proved; it is not just a range/kernel claim.

## 3. Why the tempting circuit proof is incomplete

The Julia colligation

\[
{\cal U}=
\begin{bmatrix}S&W\\V^*&0\end{bmatrix}           \tag{6}
\]

is unitary:

\[
{\cal U}^*{\cal U}={\cal U}{\cal U}^*=I.         \tag{7}
\]

A natural attempted proof unrolls copies of (6), treats \(S,S^*\)
as oppositely oriented zero-port traversals, and tries to cut an
arbitrary path through the ordinary impulse outputs

\[
\begin{bmatrix}V^*SW\\V^*S^2W\\\vdots\\V^*S^dW\end{bmatrix}
=
\begin{bmatrix}B_1^*\\B_2^*\\\vdots\\B_d^*\end{bmatrix}.       \tag{8}
\]

That picture explains (5), but it is not presently a proof.  The
missing step is an explicit block recursion showing that after mixed
orientations the selected output row has zero component on the
all-forward residual state.  Merely saying that the gates
"straighten" assumes exactly this finite-causality assertion.

Backtracking paths can revisit a port.  The first genuinely nonlinear
example

\[
w=(S^*)^2S^3(S^*)^2
\]

has

\[
K_w=-(I-B_1B_1^*)B_1=-B_1(I-B_1^*B_1),          \tag{9}
\]

so a naive Dyck-path or independent-edge argument is insufficient.
No lemma number is assigned until the residual block is eliminated
rigorously.

## 4. Relation to the two Hardy frames

L236's right Hardy analysis map gives

\[
({\cal O}_RWy)_n=V^*S^nWy=B_n^*y.
\]

Thus the right side of (1) is exactly the first-\(d\) output energy.
The candidate says that any length-\(d\) backtracking endpoint path
is a contractive readout of that causal history.  It is stronger
than L209's kernel statement.

The proposed constant is sharp.  For \(w=(S^*)^d\), \(K_w=B_d\).
At a complete delay \(B_1=\cdots=B_{d-1}=0\), (1) is equality.

## 5. Conditional consequence for the margin problem

For L311's root \(pG_jq\), set \(d=|p|+|q|\) and
\(K_{p,q}=W^*qpV\).  If (1) is proved, then

\[
K_{p,q}K_{p,q}^*\preceq\sum_{\ell=1}^dB_\ell B_\ell^*.
\]

Hence

\[
\boxed{
\pm t(K_{p,q}B_j^*+B_jK_{p,q}^*)
\preceq\varepsilon B_jB_j^*
+\frac{t^2}{\varepsilon}\sum_{\ell=1}^dB_\ell B_\ell^*.}       \tag{10}
\]

This still would not prove the prepared finite jet positive.  The
remaining issue is weight matching: a term at elliptic order \(c^n\)
must only be charged to prefix Grams with absorbable weights
\(c^{2\ell}\).  Raw word length is too crude.  The required stronger
statement is a weighted causal factorization for the complete
L306/L307 generating series, not isolated words.

## 6. Numerical evidence and regeneration

The deterministic checker exhausts all binary words through length
ten on unstructured, rank-chain, and completely delayed
colligations of copy multiplicity one through four.  All 2,046 words
per colligation satisfy (1), with only roundoff-size negative slack.
Separate untracked stress tests covered:

1. every word through length sixteen on additional general
   colligations;
2. 100 general colligations of multiplicity up to five with random
   words through length twenty; and
3. the nonlinear word (9), including noncommuting copy matrices.

These tests make (1) a serious conjecture, not a proof.

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_endpoint_word_gram.py \
  --output \
  experiments/repeated_crabb_endpoint_word_gram_s70226.jsonl
```
