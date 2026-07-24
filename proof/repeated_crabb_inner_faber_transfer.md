# The genuine matrix-inner Faber transfer on repeated Crabb equality

## 1. Result (L201, 2026-07-24)

L200 disproves the raw matrix-polynomial analogue of the scalar quotient
\(g/g^\sharp\).  L193's canonical Stein metric nevertheless supplies the
correct replacement.

At every nearby inverse-block-Toeplitz equality anchor, let \(M\) be
L193's canonical metric and put

\[
 C=M^{1/2}AM^{-1/2}.                                \tag{1}
\]

Then \(C\) is a finite pure partial isometry with equal defect
multiplicities \(m\).  For orthonormal defect frames

\[
 V:\mathbb C^m\to\ker C,\qquad
 W:\mathbb C^m\to\ker C^*,
\]

the transfer

\[
 \boxed{B_H(z)=W^*(I-zC^*)^{-1}V}                   \tag{2}
\]

is analytic and matrix inner, and

\[
 B_H(0)=0.                                          \tag{3}
\]

Write

\[
 B_H(z)=\sum_{n\ge1}B_nz^n,\qquad
 B_n=W^*(C^*)^nV.                                   \tag{4}
\]

For sufficiently small complex \(c\), its matrix Faber transform has the
exact boundary reflection

\[
 \boxed{
 ({\cal F}_cB_H)(\zeta+c/\zeta)
 =B_H(\zeta)+\sum_{n\ge1}c^nB_n\zeta^{-n}
 =B_H(\zeta)+B_H(c/\zeta).}                         \tag{5}
\]

Thus the genuine reflected copy-space Gram is

\[
 {\cal R}_c(H)
 =\sum_{n\ge1}|c|^{2n}B_n^*B_n\succeq0.             \tag{6}
\]

At the repeated Crabb apex,

\[
 B_H(z)=U z^L,\qquad U^*U=I_m,                      \tag{7}
\]

so

\[
 {\cal R}_c(H_0)=|c|^{2L}I_m.                       \tag{8}
\]

Near the apex \(B_L\) remains invertible, hence

\[
 {\cal R}_c(H)\succeq |c|^{2L}B_L^*B_L\succ0
 \quad(c\ne0).                                      \tag{9}
\]

L201 supplies the correct noncommutative inner/Faber coordinates and a
coercive reflected Gram.  It does **not** yet identify (6) with the actual
prepared similarity endpoint.  That ordered model-complement identity is
the remaining elliptic gate.

## 2. The canonical metric gives a partial isometry

Use L193's notation

\[
 M-A^*MA=QQ^*,\qquad
 Q=HE_0D^{-1/2},\qquad D=E_0^*HE_0.                 \tag{10}
\]

The exact lower endpoint says

\[
 ME_0=HE_0.
\]

Therefore

\[
 M^{-1}Q=E_0D^{-1/2},\qquad
 Q^*M^{-1}Q=I_m.                                    \tag{11}
\]

Put

\[
 V=M^{-1/2}Q.
\]

Equations (10)--(11) give

\[
 I-C^*C
 =M^{-1/2}QQ^*M^{-1/2}
 =VV^*,\qquad V^*V=I_m.                             \tag{12}
\]

Thus the right defect is an orthogonal projection of rank \(m\), so \(C\)
is a partial isometry and \(\operatorname{ran}V=\ker C\).  The nonzero
singular values of \(C^*C\) and \(CC^*\) agree, hence

\[
 I-CC^*=WW^*                                        \tag{13}
\]

for an orthonormal \(m\)-frame \(W\) of \(\ker C^*\).

Every equality anchor considered here is near the nilpotent Crabb point.
Its spectrum is therefore strictly inside the disk, so \(C\) has no
unitary summand and is pure.

## 3. The characteristic function and the forced zero

The characteristic function of the partial isometry \(C\), in the defect
frames \(V,W\), is

\[
\begin{aligned}
 \Theta_C(z)
 &=
 W^*\{-C+z( I-CC^*)^{1/2}
       (I-zC^*)^{-1}(I-C^*C)^{1/2}\}V\\
 &=zW^*(I-zC^*)^{-1}V
 =zB_H(z).                                          \tag{14}
\end{aligned}
\]

Here \(CV=0\), while the two square roots in (14) are the projections in
(12)--(13).  The standard characteristic-kernel identity, obtained by
expanding both resolvents and telescoping \(I-C^*C\), is

\[
\begin{aligned}
 I-\Theta_C(w)^*\Theta_C(z)
  =(1-\overline wz)\,
  V^*(I-\overline wC)^{-1}
       (I-zC^*)^{-1}V.                              \tag{15}
\end{aligned}
\]

For \(|z|=1\), purity permits radial passage to the boundary and the right
side vanishes.  Thus \(\Theta_C\), and therefore \(B_H=\Theta_C/z\), is
matrix inner.

It remains to prove (3), since division by \(z\) alone only makes \(B_H\)
analytic.  Equation (11) also gives

\[
 V=M^{1/2}E_0D^{-1/2}.                              \tag{16}
\]

The disk model has \(AE_1=2E_0\).  Hence

\[
 V=C\left(\frac12M^{1/2}E_1D^{-1/2}\right)
 \in\operatorname{ran}C.                            \tag{17}
\]

Since \(\ker C^*\perp\operatorname{ran}C\),

\[
 B_H(0)=W^*V=0,
\]

which proves (3) and the Taylor formula (4).

## 4. Exact matrix Faber reflection

Let \(P_0=2,P_1=x\), and

\[
 P_n(x;c)=xP_{n-1}(x;c)-cP_{n-2}(x;c).
\]

For scalar \(\zeta,c\),

\[
 P_n(\zeta+c/\zeta;c)
 =\zeta^n+c^n\zeta^{-n}.                            \tag{18}
\]

Near the Crabb point, the spectral radius of \(C\) is uniformly below one.
Consequently (4) converges on a disk of radius greater than one, and for
small \(c\) its coefficientwise Faber series converges normally on a
fixed annulus.  Substituting (18) gives

\[
\begin{aligned}
\sum_{n\ge1}B_nP_n(\zeta+c/\zeta;c)
&=\sum_{n\ge1}B_n\zeta^n
  \sum_{n\ge1}B_nc^n\zeta^{-n}\\
&=B_H(\zeta)+B_H(c/\zeta),
\end{aligned}
\]

which proves (5) without commuting any copy matrices.

Matrix innerness also gives the Parseval identity

\[
 \sum_{n\ge1}B_n^*B_n=I_m.                          \tag{19}
\]

The reflected energy (6) is therefore an oriented, bounded Gram, not a
formal scalar norm.

## 5. The Crabb monomial and coercivity

At \(H_0=I/2\), the partial isometry \(C\) is the length-\(L+1\) model
shift tensored with \(I_m\).  Its right and left defects are the two
endpoint copy spaces.  Hence

\[
 W^*(C^*)^nV=0\quad(n\ne L),\qquad
 W^*(C^*)^LV=U
\]

for a unitary frame change \(U\).  This proves (7)--(8).

All objects in (10)--(14) are analytic after choosing a local analytic
left defect frame.  Therefore \(B_L\) remains invertible on a fixed
equality neighbourhood.  Keeping only its summand in (6) proves (9).

This terminal coercivity is the matrix-inner version of L190's
\(-16|c^L|^2\) leg.  Lower Taylor coefficients encode the noncommuting
equality amplitudes, while the higher rational tail records the
noncommutative spectral-factor correction missing from L200's raw finite
polynomial.

## 6. Remaining endpoint identification

L201 resolves the innerness and multiplication-order problem.  Two
interfaces are still required for a repeated elliptic theorem.

1. Lift L149's total one-reflection identity to the rank-\(m\)
   colligation, including Riemann/Faber preparation, operator motion,
   denominator, defect frames, and endpoint normalization.
2. Lift L150's orbit-complement metric so that its first nonzero endpoint
   is a negative operator multiple of the Gram (6), with no unaccounted
   left-oriented or earlier-flag term.

Only after those identities are proved may (9) be used to dominate the
analytic remainder uniformly through copy-rank changes.

## 7. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_inner_faber_transfer.py \
  --output \
  experiments/repeated_crabb_inner_faber_transfer_s70224.jsonl
```

The checker constructs strengthened noncommuting inverse-block-Toeplitz
anchors and verifies:

1. both defect operators are rank-\(m\) orthogonal projections;
2. the defect frames are orthonormal and \(B_H(0)=0\);
3. \(B_H\) is unitary on the sampled circle;
4. its Taylor coefficients satisfy Parseval;
5. the Dickson/Faber series equals (5);
6. the reflected Gram is positive definite; and
7. the Crabb transfer is exactly one degree-\(L\) unitary monomial.

The standard records cover lengths two through five and multiplicities two
and three.  The tracked dataset regenerates byte for byte with SHA-256

```text
48cbe3f080678bc57a6b326ddbe43fa648c1518cc1bf91ba0dde6b59a05a9330
```

Equations (10)--(19), not the floating audit, prove L201.
