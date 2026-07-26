# Every endpoint word has a constructive transfer-prefix normal form

## 1. Result (L312, 2026-07-26)

Let

\[
I-S^*S=VV^*,\qquad I-SS^*=WW^*,\qquad V^*W=0
\]

be a balanced partial-isometry colligation, and put

\[
B_j=W^*(S^*)^jV\qquad(j\ge1).
\]

For a word \(w\) of length \(d\) in \(S,S^*\), define

\[
K_w=W^*w(S,S^*)V.
\]

There are explicitly constructed copy matrices
\(C_1(w),\ldots,C_d(w)\) such that

\[
\boxed{K_w=\sum_{j=1}^dB_jC_j(w)}                \tag{1}
\]

and, for the stacked column

\[
{\cal C}_w=[C_1(w)^*\ \cdots\ C_d(w)^*]^*,
\]

\[
\boxed{\|{\cal C}_w\|\le\left\lceil\frac d2\right\rceil.}       \tag{2}
\]

Consequently

\[
\boxed{
K_wK_w^*
\preceq
\left\lceil\frac d2\right\rceil^2
\sum_{j=1}^dB_jB_j^*.}                           \tag{3}
\]

This is a polynomial factorization.  It uses no transfer inverse,
pseudoinverse, constant-rank assumption, purity, or state-dimension
constant.  It strengthens L209's endpoint-kernel statement for a
fixed word and supplies a proved fallback to A259's still-open
constant-one conjecture.

## 2. Recursive construction

Write \(a=S^*\) and \(s=S\).  Define the coefficients by reducing
the current word from its left endpoint.

1. If the current word begins with \(s\), stop and emit nothing.
2. If it is \(a^r\), emit \(C_r\mathrel{+}=I\) and stop.
3. Otherwise write it uniquely as \(a^p s v\), with \(p\ge1\).
   Emit

   \[
   C_{p-1}\mathrel{-}=V^*v(S,S^*)V               \tag{4}
   \]

   when \(p>1\), replace the current word by \(a^{p-1}v\), and
   continue.  When \(p=1\), the omitted coefficient would multiply
   \(B_0=W^*V=0\).

Every pass through step 3 deletes two letters, so the procedure
terminates and emits at most \(\lceil d/2\rceil\) matrices.

## 3. Proof of the identity

The partial-isometry defect identity gives

\[
a^ps=a^{p-1}S^*S=a^{p-1}(I-VV^*).                \tag{5}
\]

Therefore

\[
\begin{aligned}
W^*a^psvV
&=W^*a^{p-1}vV\\
&\quad-W^*a^{p-1}V\,V^*vV\\
&=K_{a^{p-1}v}-B_{p-1}\{V^*vV\}.                \tag{6}
\end{aligned}
\]

Equation (6) is exactly step 3.  A word beginning with \(s\) has
zero endpoint because \(W^*S=0\), while a pure word \(a^r\) has
endpoint \(B_r\).  Induction over the successive two-letter
deletions proves (1), with every multiplication order retained.

## 4. Norm and Gram bounds

Every emitted matrix in (4) is a contraction:

\[
\|V^*v(S,S^*)V\|\le1,                            \tag{7}
\]

because \(S,S^*,V,V^*\) are contractions.  The terminal identity is
also a contraction.  If \(N\le\lceil d/2\rceil\) is the number of
emissions, then for every copy vector \(x\),

\[
\begin{aligned}
\|{\cal C}_wx\|
&=\left(\sum_j\|C_j(w)x\|^2\right)^{1/2}\\
&\le\sum_j\|C_j(w)x\|
\le N\|x\|.
\end{aligned}
\]

This proves (2).  With
\({\mathbb B}_d=[B_1\ \cdots\ B_d]\), equation (1) is
\(K_w={\mathbb B}_d{\cal C}_w\).  Hence

\[
K_wK_w^*
\preceq\|{\cal C}_w\|^2{\mathbb B}_d{\mathbb B}_d^*,
\]

which is (3).

The normal form need not be contractive.  On general colligations,
words first appearing at length seven already give
\(\|{\cal C}_w\|>1\).  Thus this proof does not establish A259.

## 5. Weight information and scope

Every emitted index \(j\) is at most the number of \(S^*\) letters
in the original word.  This is useful for a finite elliptic word
whose coefficient degree pays for each reverse letter.

It does **not** by itself close the all-series weight problem.
Stein/Hardy closure can create arbitrarily long raw words before
the coisometric cancellations of L236--L245 are grouped.  Applying
(1) term by term to that infinite expansion would lose the physical
elliptic valuation.  The live task remains a grouped weighted normal
form for L306/L307, preferably in L220's Schur-orthogonal
coordinates.

## 6. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_endpoint_word_gram.py \
  --output \
  experiments/repeated_crabb_endpoint_word_gram_s70226.jsonl
```

The checker exhausts every word through length ten on ten
unstructured, rank-chain, and completely delayed colligations.  It
independently reconstructs (1), audits (2), and separately tests
A259's sharper conjecture.  The exact recursion (5)--(7), not the
floating audit, proves L312.
