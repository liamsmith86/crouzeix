# Scalar equality channels in the repeated Crabb inner transfer

## 1. Result (L205, 2026-07-24)

Let \(T\) be one L193 inverse-block-Toeplitz disk anchor, let

\[
P=2I-VV^*+2WW^*,\qquad C=P^{1/2}TP^{-1/2},
\]

and let L201's genuine transfer be

\[
B(z)=W^*(I-zC^*)^{-1}V
=\sum_{n\ge1}B_nz^n.                               \tag{1}
\]

Then

\[
\boxed{
\sup_{\|f\|_{\mathbb D}\le1}\|f(T)\|=2
}
\]

if and only if there are unit vectors
\(v,u\in\mathbb C^m\) and a scalar inner function \(g\), with
\(g(0)=0\), such that

\[
\boxed{
B(z)v=g(z)u.
}                                                   \tag{2}
\]

Thus scalar equality is much more rigid than complete similarity
equality.  Every L193 anchor has complete similarity square four, but a
scalar function attains two only when the matrix-inner transfer contains
a constant one-dimensional scalar inner channel.

Moreover, (2) orthogonally splits the transfer:

\[
B(z)=g(z)\oplus B'(z)                               \tag{3}
\]

in the domain decomposition
\(\mathbb Cv\oplus v^\perp\) and range decomposition
\(\mathbb Cu\oplus u^\perp\).  The pure partial-isometry model \(C\),
the three-level metric \(P\), and the physical disk operator \(T\)
split accordingly.  Iterating extracts every scalar equality block and
leaves a complementary matrix-inner block with scalar norm strictly
below two.

This theorem does not by itself prove a uniform repeated neighbourhood:
when a scalar channel nearly appears, its strict gap can vanish together
with elliptic and transverse coordinates.  It does show that a proof of
the scalar conjecture need not establish L203's stronger matrix endpoint
on a genuinely channel-free anchor.

## 2. Equality in the condition-two chain

For every scalar Schur function \(f\),

\[
\begin{aligned}
\|f(T)x\|
&=\|P^{-1/2}f(C)P^{1/2}x\|\\
&\le \|P^{-1/2}\|\,\|f(C)\|\,\|P^{1/2}x\|
\le2\|x\|.                                         \tag{4}
\end{aligned}
\]

The endpoint eigenspaces of \(P\) are exactly
\(\operatorname{ran}V\) at one and \(\operatorname{ran}W\) at four.
If equality holds in (4), every inequality in the chain is equality.
Therefore, for some unit defect coordinates \(u,v\),

\[
f(C)Wu=Vv.                                         \tag{5}
\]

Conversely, (5) makes (4) an equality, because

\[
P^{1/2}Wu=2Wu,\qquad P^{-1/2}Vv=Vv.
\]

So it remains only to classify (5).

## 3. Parseval and Cauchy--Schwarz rigidity

Write \(f(z)=\sum_{n\ge0}a_nz^n\).  Since \(V^*W=0\), (1) gives

\[
V^*f(C)Wu=\sum_{n\ge1}a_nB_n^*u.                   \tag{6}
\]

Taking the inner product with \(v\) in (5),

\[
1=\left|
\sum_{n\ge1}a_n\overline{u^*B_nv}
\right|.                                           \tag{7}
\]

The scalar Schur bound and L201's matrix Parseval identity give

\[
\sum_{n\ge1}|a_n|^2
\le\sum_{n\ge0}|a_n|^2\le1,                        \tag{8}
\]

and

\[
\sum_{n\ge1}|u^*B_nv|^2
\le\sum_{n\ge1}\|B_nv\|^2
=v^*\left(\sum_{n\ge1}B_n^*B_n\right)v=1.          \tag{9}
\]

Equality in (7) forces equality everywhere in (8)--(9).  Hence

\[
a_0=0,\qquad
(I-uu^*)B_nv=0\quad(n\ge1),                        \tag{10}
\]

and the two scalar coefficient sequences in (7) are parallel.  Put

\[
b_n=u^*B_nv,\qquad g(z)=\sum_{n\ge1}b_nz^n.
\]

Equation (10) gives \(B_n v=b_nu\), hence (2).  Since \(B\) is
matrix inner,

\[
1=\|B(e^{it})v\|^2=|g(e^{it})|^2
\quad\text{a.e.},
\]

so \(g\) is scalar inner.  This proves necessity.

Conversely, suppose (2) holds and take \(f=g\).  Then

\[
v^*V^*g(C)Wu
=\sum_{n\ge1}b_n\overline{b_n}=1.                  \tag{11}
\]

The operator \(g(C)\) is a contraction.  Its compression in (11)
already has norm one, so \(g(C)Wu=Vv\).  Equation (5) and hence scalar
equality follow.

## 4. Reduction of the colligation

On the circle \(B(\zeta)\) is unitary.  If (2) holds and
\(x\perp v\), then

\[
\langle B(\zeta)x,u\rangle
=\langle x,B(\zeta)^*u\rangle
=\langle x,\overline{g(\zeta)}v\rangle=0.
\]

Analytic continuation proves the off-diagonal transfer blocks vanish,
which is (3).

The characteristic function of \(C\) is \(\Theta_C(z)=zB(z)\).
L201's characteristic-kernel identity realizes the state space as the
finite model space with kernel

\[
\frac{I-\Theta_C(w)^*\Theta_C(z)}
     {1-\overline wz}.                              \tag{12}
\]

The block decomposition (3) makes (12) block diagonal.  Its two
reproducing-kernel spans are orthogonal and invariant for the model
operator and its adjoint.  Purity makes the realization minimal, so
they pull back to reducing state subspaces for \(C\).

The defect projections \(VV^*,WW^*\) respect the same split.  Therefore

\[
P=2I-VV^*+2WW^*
\]

and \(T=P^{-1/2}CP^{1/2}\) also reduce.  The scalar channel is precisely
a scalar full-Hardy disk-equality block of the kind controlled by
L187--L192.

## 5. Consequence for strategy

There are now two distinct equality notions.

1. **Complete equality:** every L193 anchor has \(t_*(T)=4\), witnessed
   by its matrix-valued characteristic transfer.
2. **Scalar equality:** only anchors containing a channel (2) attain
   scalar Crouzeix norm two.

Therefore failure to construct L203's full matrix endpoint would block
the stronger completely bounded route, but need not block the scalar
Crouzeix conjecture.  A scalar repeated proof may instead combine:

- the quantitative channel leakage in (9);
- L190 on every extracted scalar channel;
- L195--L197 on transverse disk residuals; and
- L199's support-Jensen circular-normal gate.

The remaining issue is a uniform mixed-scale estimate when the channel
leakage, reflected coefficient, and transverse residual vanish at
comparable orders.

## 6. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_scalar_channel.py \
  --output \
  experiments/repeated_crabb_scalar_channel_s70224.jsonl
```

The checker verifies matrix Parseval, exact score one at the repeated
Crabb apex, and a strict product-state leakage at the standard
noncommuting anchors through lengths five and multiplicities three.
The alternating product-state maximization is a diagnostic, not the
proof of L205; the proof is the equality analysis (4)--(12).  The
standard dataset SHA-256 is
`c282527cc63bf732cea8065e800959568a4d2bdd1f451f107ea8a3fffab3df4c`.
