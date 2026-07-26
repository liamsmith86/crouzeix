# Every physical rooted endpoint has an ideal-valued normalization

## 1. Result (L311, 2026-07-26)

Retain L305/L310's balanced pure partial-isometry colligation

\[
I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,
\qquad
B_j=W^*(S^*)^jV,
\]

and write

\[
G_j=ES^jF=VB_j^*W^*.
\]

For one physical rooted word

\[
T=pG_jq,
\]

put

\[
K=W^*qpV.
\]

L310 converts L305's rooted transfer-channel term to the literal
hereditary factor

\[
KB_j^*+B_jK^*
\]

but its standalone polarized column need not vanish when \(B_j=0\).
That column is unnecessary for a physical root.

Define the state cross

\[
\boxed{X=WKV^*=FqpE}                              \tag{1}
\]

and the Hermitian relative commutator

\[
\boxed{
R_{p,q,j}=[G_j,X]+[G_j,X]^*.}                    \tag{2}
\]

Then

\[
\boxed{
\begin{aligned}
W^*{\cal G}_S(R_{p,q,j})W
={}&\Phi(B_j^*K+K^*B_j)\\
&-(KB_j^*+B_jK^*).
\end{aligned}}                                    \tag{3}
\]

Moreover, L305 gives an explicit perpendicular column
\(C_{p,q,j}^{\rm rel}\) such that

\[
\boxed{
W^*{\cal G}_S(R_{p,q,j})W
={\cal M}_S(C_{p,q,j}^{\rm rel}),}                \tag{4}
\]

\[
\boxed{
C_{p,q,j}^{\rm rel}V^*\in{\cal I}_j,\qquad
\|C_{p,q,j}^{\rm rel}\|_F
\le(4|p|+4|q|+16)\|B_j\|_F.}                     \tag{5}
\]

Combining (3)--(4) with L305's original column
\(C_{p,q,j}\) gives the normalized physical split

\[
\boxed{
\begin{aligned}
W^*{\cal G}_S(T+T^*)W
={}&KB_j^*+B_jK^*\\
&+{\cal M}_S(\widehat C_{p,q,j}),
\end{aligned}}                                    \tag{6}
\]

where

\[
\widehat C_{p,q,j}
=C_{p,q,j}+C_{p,q,j}^{\rm rel},
\]

\[
\boxed{
\widehat C_{p,q,j}V^*\in{\cal I}_j,\qquad
\|\widehat C_{p,q,j}\|_F
\le(5|p|+4|q|+17)\|B_j\|_F.}                     \tag{7}
\]

In particular,

\[
\boxed{B_j=0\quad\Longrightarrow\quad
\widehat C_{p,q,j}=0.}                            \tag{8}
\]

Thus every physical rooted channel simultaneously has:

1. a literal bounded hereditary endpoint factor;
2. automatic L292 Smith valuations;
3. a polynomial, ideal-valued response selection; and
4. exact normalization on every complete-delay stratum.

This closes L310's response-normalization debt term by term.  No
pseudoinverse, endpoint-null gauge, or abstract cyclic-excision
theorem is used.

## 2. Endpoint identity

From \(G_j=VB_j^*W^*\) and \(X=WKV^*\),

\[
\begin{aligned}
G_jX&=VB_j^*KV^*,\\
XG_j&=WKB_j^*W^*.
\end{aligned}
\]

Hence

\[
\boxed{
R_{p,q,j}
=V(B_j^*K+K^*B_j)V^*
-W(KB_j^*+B_jK^*)W^*.}                           \tag{9}
\]

The first term in (9) has upper Stein endpoint
\(\Phi(B_j^*K+K^*B_j)\).  The second has endpoint
\(KB_j^*+B_jK^*\): for \(n\ge1\),

\[
W^*S^n=0,
\]

so only the \(n=0\) term survives in
\(W^*{\cal G}_S(WYW^*)W\).  This proves (3).

## 3. Constructive ideal-valued column

Put \(r=qp\).  Since

\[
F=I-SS^*,\qquad E=I-S^*S,
\]

(1) has the four-word expansion

\[
\boxed{
X=r-rS^*S-SS^*r+SS^*rS^*S.}                    \tag{10}
\]

For each monomial \(x\) in (10), expand

\[
[G_j,x]=G_jx-xG_j.
\]

Apply L305 to \(G_jx\) with prefix \(1\), suffix \(x\), and to
\(xG_j\) with prefix \(x\), suffix \(1\).  Their rooted copy
matrices are identical:

\[
B_j^*W^*xV.
\]

They therefore cancel exactly in the difference.  L305 leaves the
explicit response column

\[
C_{1,x,j}-C_{x,1,j},
\]

whose lift lies in \({\cal I}_j\).  Sum these four differences with
the signs in (10), proving (4) and the ideal statement in (5).

For a word \(x\), L305 gives

\[
\|C_{1,x,j}-C_{x,1,j}\|_F
\le(|x|+2)\|B_j\|_F.
\]

The four lengths in (10) are

\[
|r|,\quad |r|+2,\quad |r|+2,\quad |r|+4.
\]

Their bounds sum to

\[
(4|r|+16)\|B_j\|_F,
\]

which is (5).  Add L305's
\((|p|+1)\|B_j\|_F\) bound to obtain (7).

Every summand of both columns contains \(B_j\) or \(B_j^*\)
linearly.  This proves (8) and gives polynomial dependence without
appealing to a limiting argument.

Terms initially rooted at \(G_j^*\) follow by applying (6) to their
adjoints.  Therefore any finite coefficient in the Hermitian bridge
ideal has the same normalized decomposition.  Since L289 requires
only the finite jet through the terminal grade, the linear
word-length constant in (7) is sufficient; no infinite-series
summability claim is needed.

## 4. Scope and next gate

L311 closes the endpoint-factor, Smith-divisibility, response-module,
and complete-delay-normalization parts of the rooted recurrence.  It
does **not** prove positivity of the complete prepared coefficient.

The remaining elliptic gate is now the finite-jet margin assembly:

1. collect the literal factors from (6) in L306/L307's affine
   successor;
2. at odd orders absorb their cumulative-transfer factors by
   previously retained positive margins; and
3. at even orders isolate the new direct Gram and charge only the
   prior-transfer cost to earlier margins.

L308--L309 forbid replacing this task by global odd-response parity.
L307 forbids treating its positive partial-scale Gram without the
affine leftover.  Another isolated fifth/seventh computation would
not establish the required uniform square-completion rule.

## 5. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_rooted_endpoint_normalization.py \
  --output \
  experiments/repeated_crabb_rooted_endpoint_normalization_s70226.jsonl
```

The checker expands (10) independently, compares it with the direct
defect-projection product, verifies (3)--(7), reconstructs every
ideal lift, and tests unstructured, rank-chain, and complete-delay
colligations.  On the three complete-delay records the raw L310
polarized columns have norms approximately \(1.02,1.21,1.44\), while
the normalized physical columns vanish exactly.  All 24 records
pass.

An independent exact-rational SymPy audit used L308's three-state
colligation with \(j=1,p=(S^*)^2,q=1\).  It gives

\[
B_1=-\frac45,\quad K=\frac9{25},\quad
C_{p,q,1}=-\frac{64}{125}W,\quad
C_{p,q,1}^{\rm rel}=+\frac{64}{125}W,
\]

so \(\widehat C=0\) and the literal endpoint is
\(-72/125\), with every matrix identity in (1)--(7) holding exactly.

The tracked dataset SHA-256 is

```text
785aaab8462afccce9c5693a31677c0215c7fd714bcc37ec28fa512bc37ec8c6
```
