# Every rooted channel is a literal hereditary endpoint factor modulo response

## 1. Result (L310, 2026-07-26)

Retain L280/L305's balanced pure partial-isometry colligation,

\[
I-S^*S=VV^*,\qquad I-SS^*=WW^*,
\qquad
B_j=W^*(S^*)^jV,
\]

the transfer channel

\[
\Phi(A)=\sum_{n\ge1}B_nAB_n^*,
\]

and the balanced endpoint response \({\cal M}_S\).

For an arbitrary copy matrix \(K\), define L280's polarized column

\[
\boxed{
{\cal C}_j(K)
=Q\left\{
 S^jWK+
 \sum_{\ell=1}^{j-1}
 (S^*)^{j-\ell}VK^*B_\ell
\right\},\qquad Q=I-VV^*.}                       \tag{1}
\]

Then

\[
V^*{\cal C}_j(K)=0
\]

and the rooted transfer channel has the exact conversion

\[
\boxed{
\begin{aligned}
\Phi(B_j^*K+K^*B_j)
={}&KB_j^*+B_jK^*\\
&-{\cal M}_S({\cal C}_j(K)).
\end{aligned}}                                    \tag{2}
\]

Moreover,

\[
\boxed{\|{\cal C}_j(K)\|_F\le j\|K\|_F.}          \tag{3}
\]

Thus a channel root at delay \(j\) is response-equivalent to the
literal hereditary endpoint factor

\[
\boxed{KB_j^*+B_jK^*.}                           \tag{4}
\]

This factor:

1. vanishes exactly on \(\ker B_j^*\);
2. has an explicit bounded analytic multiplier \(K\);
3. satisfies L292's Smith valuations automatically along every
   analytic rank-changing arc; and
4. requires no pseudoinverse or rank projection.

Combining (2) with L305 gives the rooted-word form needed after L309.
For

\[
T=p(ES^jF)q,\qquad
K_{p,q}=W^*qpV,
\]

L305 and (2) give

\[
\boxed{
\begin{aligned}
W^*{\cal G}_S(T+T^*)W
={}&K_{p,q}B_j^*+B_jK_{p,q}^*\\
&+{\cal M}_S\!\left(
 C_{p,q,j}-{\cal C}_j(K_{p,q})
 \right).
\end{aligned}}                                    \tag{5}
\]

Since words in \(S,S^*\) are contractions,

\[
\|K_{p,q}\|\le1
\]

and L305's bound yields

\[
\boxed{
\|C_{p,q,j}-{\cal C}_j(K_{p,q})\|_F
\le (|p|+1)\|B_j\|_F+j\|K_{p,q}\|_F.}            \tag{6}
\]

For any finite rooted presentation, (5) therefore converts the
complete retained quotient into:

1. an explicit bounded sum of hereditary factors through its
   cumulative transfer row; and
2. a bounded analytic response.

This closes the **endpoint factor/Smith-divisibility** part of the
corrected L309 odd target.  It also applies at even orders, where the
remaining task is the sign and domination of the hereditary factors.

## 2. Proof from polarized L280

L280 proves, in physical normalization,

\[
{\cal M}_T(P_{\rm met}^{1/2}{\cal C}_j(K))
=8\left\{
\operatorname{sym}(B_jK^*)
-\Phi(\operatorname{sym}(K^*B_j))
\right\}.                                        \tag{7}
\]

The physical upper endpoint is four times the balanced endpoint, so

\[
{\cal M}_T(P_{\rm met}^{1/2}C)=4{\cal M}_S(C).
\]

Divide (7) by four and use

\[
2\operatorname{sym}(B_jK^*)=B_jK^*+KB_j^*,
\]

\[
2\operatorname{sym}(K^*B_j)=K^*B_j+B_j^*K.
\]

This gives

\[
{\cal M}_S({\cal C}_j(K))
=B_jK^*+KB_j^*
-\Phi(K^*B_j+B_j^*K),
\]

which is (2).

Perpendicularity follows from \(V^*Q=0\).  Each term in (1) is
bounded by \(\|K\|_F\), because the state words, frames, and every
\(B_\ell\) are contractions.  There are \(j\) terms, proving (3).

L305's copy root for \(p(ES^jF)q\) is exactly

\[
B_j^*K_{p,q}.
\]

Its Hermitian completion is the left side of (2).  Substitute (2)
into L305's endpoint split to obtain (5), and combine (3) with
L305's response bound to obtain (6).

## 3. Exact scope guard: the response normalization is still open

Equation (5) does **not** complete the arbitrary-grade preparation.
The polarized converting column (1) need not belong to L305's
bridge ideal and need not vanish when \(B_j=0\).

The failure is already visible at a complete first delay.  If

\[
B_1=0,
\]

then both the channel value and the endpoint factor in (2) vanish,
and

\[
{\cal M}_S({\cal C}_1(K))=0.
\]

Nevertheless

\[
{\cal C}_1(K)=QS WK
\]

can have norm exactly \(\|K\|_F\).  It is a nonzero endpoint-null
fixed-base direction.  Inserting it blindly would revive the
normalization problem behind L284--L288: later moving successors can
be nonzero even though the current endpoint vanishes.

Therefore L310 separates the remaining recurrence debt sharply:

> the hereditary endpoint quotient and its Smith valuations are
> closed; one must still normalize the **sum** of the polarized
> response columns compatibly with every complete-delay flag, or
> prove that their later affine successors remain bounded hereditary
> factors with sufficient margins.

This is not a request for another isolated quintic/sextic
calculation.  It is the arbitrary-grade normalization law for (1),
coupled to L301 and L307.

## 4. Independent regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_rooted_endpoint_factor.py \
  --output \
  experiments/repeated_crabb_rooted_endpoint_factor_s70226.jsonl
```

The checker verifies (2)--(3) for delays one through four on
unstructured colligations, verifies (4) on nontrivial rank-chain
flags, and audits three complete-first-delay cases where the
converting column has norm one while its transfer, endpoint factor,
and response all vanish.  All 21 records pass.  The tracked dataset
SHA-256 is

```text
376523b537a27d4f30f65705e8e6ab536220bf65256901294d931777a7064ae1
```
