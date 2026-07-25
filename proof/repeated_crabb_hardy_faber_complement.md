# The delayed target is the complementary Hardy--Faber channel

## 1. Result (L254, 2026-07-25)

Retain L236's pure partial isometry and two Hardy analysis maps

\[
({\cal O}_Rx)_n=V^*S^nx,\qquad
({\cal O}_Lx)_j=W^*(S^*)^jx,
\]

with matrix-inner transfer

\[
B(z)=\sum_{h\ge1}B_hz^h.
\]

Let \(L^*\) be the backward shift on \(H^2(\mathbb C^m)\), put

\[
\widehat{\cal O}_L=L^*{\cal O}_L,
\]

and let \({\cal T}_B\) be multiplication by \(B\).  Then the two
output channels are exact complements:

\[
\boxed{
\widehat{\cal O}_L\widehat{\cal O}_L^*
+{\cal T}_B{\cal T}_B^*=I.}                       \tag{1}
\]

Equivalently,

\[
\widehat{\cal O}_L\widehat{\cal O}_L^*
=P_{{\cal K}_B},\qquad
{\cal K}_B=H^2(\mathbb C^m)\ominus
B H^2(\mathbb C^m).                               \tag{2}
\]

Consequently L253's window energy is the missing energy of the
shifted L236 left-Hardy channel:

\[
\boxed{
\|P_k{\cal T}_BP_k\|_{\rm HS}^2
=\operatorname {tr}
\{P_k(I-\widehat{\cal O}_L\widehat{\cal O}_L^*)P_k\}.} \tag{3}
\]

There is also an exact generating bridge to L201's reflected Faber
Gram.  For \(q=|c|^2<1\), define

\[
{\cal R}_c(B)=\sum_{h\ge1}q^hB_h^*B_h
\]

and

\[
{\cal E}_k(B)=\|P_k{\cal T}_BP_k\|_{\rm HS}^2.
\]

Then

\[
\boxed{
\sum_{k\ge1}q^k{\cal E}_k(B)
=\frac{\operatorname {tr}{\cal R}_c(B)}
       {(1-q)^2}.}                                \tag{4}
\]

Under the complete delay \(B_1=\cdots=B_{k-1}=0\),

\[
[q^k]\frac{4\operatorname {tr}{\cal R}_c(B)}
              {(1-q)^2}
=4\|B_k\|_F^2.                                    \tag{5}
\]

Thus A194's target is simultaneously:

1. four times a causal Toeplitz-window energy;
2. four times the complementary energy to the shifted left Hardy
   frame; and
3. the first delayed coefficient of the double Abel transform of the
   reflected Faber Gram.

Equations (1)--(5) are exact in every grade.  They do **not** prove
that L250's physical whitened volume equals this complementary
channel.  They sharpen that remaining statement: L244's coisometry
must remove the
\(\widehat{\cal O}_L\widehat{\cal O}_L^*\) channel and leave its
\({\cal T}_B{\cal T}_B^*\) complement, with L251's physical factor
two retained on both sides of the closed channel.

## 2. Hardy-channel proof

The characteristic function of \(S\), with the present defect-frame
orientation, is

\[
\Theta(z)=zB(z).
\]

The standard characteristic-kernel telescope gives

\[
\frac{I-\Theta(z)\Theta(w)^*}{1-z\overline w}
=\sum_{i,j\ge0}
z^i\overline w^{\,j}
W^*(S^*)^iS^jW.                                   \tag{6}
\]

Hence \({\cal O}_L{\cal O}_L^*\) is the model projection for
\({\cal K}_{zB}\).  Since

\[
{\cal K}_{zB}
=\mathbb C^m\oplus z{\cal K}_B,                   \tag{7}
\]

applying the backward shift on both sides gives

\[
\widehat{\cal O}_L\widehat{\cal O}_L^*
=P_{{\cal K}_B}.
\]

Matrix innerness gives

\[
{\cal T}_B{\cal T}_B^*=I-P_{{\cal K}_B},
\]

which proves (1)--(2).  Causality gives
\(P_k{\cal T}_B=P_k{\cal T}_BP_k\).  Taking the finite trace of
(1) therefore proves (3).

The one-row shift in \(\widehat{\cal O}_L\) is load-bearing.  The
unshifted L236 left frame realizes \({\cal K}_{zB}\), not
\({\cal K}_B\); dropping that shift gives the wrong matrix identity
for non-scalar copy spaces.

## 3. Abel--Faber proof

L253 gives

\[
{\cal E}_k(B)
=\sum_{h=1}^k(k+1-h)\|B_h\|_F^2.                  \tag{8}
\]

Absolute convergence permits reversing the two sums:

\[
\begin{aligned}
\sum_{k\ge1}q^k{\cal E}_k(B)
&=\sum_{h\ge1}\|B_h\|_F^2
  \sum_{k\ge h}(k+1-h)q^k\\
&=\sum_{h\ge1}\frac{q^h}{(1-q)^2}\|B_h\|_F^2\\
&=\frac{\operatorname {tr}{\cal R}_c(B)}
        {(1-q)^2}.
\end{aligned}
\]

This proves (4).  Under complete delay, only the \(h=k\) term can
contribute to coefficient \(q^k\), proving (5).

## 4. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_toeplitz_window_energy.py \
  --output \
  experiments/repeated_crabb_toeplitz_window_energy_s70224.jsonl
```

For general and completely delayed colligations through grade six,
the checker now also verifies:

1. the finite compression of (1), using rows \(1,\ldots,k+1\) of
   \({\cal O}_L\); and
2. coefficientwise identity (4), equivalently that applying
   \((1-q)^2\) to the window-energy series returns the reflected
   Faber coefficient norms.

The tracked dataset regenerates with SHA-256
`0f6c8973dd891d162889ca2b3b780a1fb79144a88819962ecfd244098e3f5249`.
