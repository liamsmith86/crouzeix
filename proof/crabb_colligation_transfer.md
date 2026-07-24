# Rank-one colligation transfer for the optimized metric (L169, 2026-07-23)

## 1. General theorem

Let `T` be stable and let `P>0` solve

\[
P-T^*PT=qq^*.                                      \tag{1}
\]

Suppose

\[
P^{-1}-TP^{-1}T^*=rr^*.                            \tag{2}
\]

Then there are a scalar `d` and a unimodular scalar `omega` such that

\[
\boxed{
\theta(\zeta)
=d+\zeta q^*(I-\zeta T)^{-1}r
=\omega{\det(\zeta I-T^*)\over\det(I-\zeta T)}.
}                                                    \tag{3}
\]

The rational function `theta` is inner.  Equivalently, in the
outside resolvent variable,

\[
\boxed{
\theta(1/z)=d+q^*(zI-T)^{-1}r.
}                                                    \tag{4}
\]

Thus the correct covariant endpoint transfer is not the resolvent
term alone.  It includes the scalar feedthrough `d`, which is the
constant term of the determinant quotient up to phase.

## 2. Orthogonal colligation proof

Balance (1) by putting

\[
S=P^{1/2}TP^{-1/2},\qquad
a=P^{-1/2}q.
\]

Then

\[
I-S^*S=aa^*.                                       \tag{5}
\]

Equation (2), with `b=P^(1/2)r`, gives

\[
I-SS^*=bb^*.                                       \tag{6}
\]

The first `p` columns of

\[
\begin{pmatrix}S\\a^*\end{pmatrix}
\]

are orthonormal in dimension `p+1`.  Their one-dimensional
orthogonal complement may be phased so its upper block is `b`.
Therefore there is a scalar `d` for which

\[
U=\begin{pmatrix}S&b\\a^*&d\end{pmatrix}            \tag{7}
\]

is unitary.

The transfer function of (7) is

\[
\theta(\zeta)=d+\zeta a^*(I-\zeta S)^{-1}b.
\]

The similarity factors cancel inside the scalar product, giving the
first expression in (3).

For the determinant formula, define

\[
M(\zeta)=
\begin{pmatrix}
I-\zeta S&-\zeta b\\
a^*&d
\end{pmatrix}.
\]

Schur complementation gives

\[
\det M(\zeta)=\det(I-\zeta S)\theta(\zeta).          \tag{8}
\]

On the other hand,

\[
M(\zeta)U^*
=
\begin{pmatrix}
S^*-\zeta I&a\\
0&1
\end{pmatrix},
\]

so

\[
\det M(\zeta)=\det(U)\det(S^*-\zeta I).              \tag{9}
\]

Combining (8)--(9), using similarity invariance of the determinants,
and setting

\[
\omega=(-1)^p\det U
\]

proves (3).  Since `U` is unitary, `|omega|=1`; the numerator and
denominator in (3) have equal modulus on the unit circle.  Hence
`theta` is inner.

## 3. Specialization to the optimized Crabb branch

For L118's optimized real persymmetric branch, L164 gives

\[
JP^{-1}J=\alpha P.
\]

Comparing the rank-one forcing in the reversed inverse Stein equation
with (1) permits the analytic choice

\[
\boxed{r=\sqrt\alpha\,Jq.}                          \tag{10}
\]

Thus the canonical characteristic transfer is available directly
from the one optimized defect:

\[
\theta(\zeta)
=d+\zeta\sqrt\alpha\,
q^T(I-\zeta T)^{-1}Jq.                              \tag{11}
\]

At the Crabb point,

\[
P=\operatorname{diag}(1,2,\ldots,2,4),\quad
q=e_0,\quad r=\tfrac12e_L,\quad d=0,
\]

and (11) is `theta(zeta)=zeta^p`.  Formula (4) becomes one half of
L168's raw endpoint resolvent.

For an analytic family of optimized metrics, (3) gives on the circle

\[
\boxed{\operatorname{Re}{\dot\theta\over\theta}=0.} \tag{12}
\]

This identity is exact to every order and automatically carries the
moving characteristic polynomial.

## 4. Why the feedthrough matters

At the Crabb point, the phase-normalized determinant formula is

\[
\theta(\zeta)
={\det(\zeta I-T^*)\over\det(I-\zeta T)}.            \tag{13}
\]

L167's leading grade-`k` characteristic variation therefore gives

\[
D\log\theta
=-k(\zeta^{-m}-\zeta^m)
\quad\hbox{after outside-variable identification}. \tag{14}
\]

L168's residual endpoint response is the negative of (14).  This
explains, rather than merely matches, its anti-self-reciprocal form.

The grade-three adversary A115 shows that `d` need not remain zero at
later mixed weights.  Equation (3) includes this motion
automatically.  Any proof based only on
`q^*(zI-T)^(-1)r` would omit the feedthrough and repeat the
fixed-zero error.

## 5. Consequence and remaining identification

L169 does not by itself prove L163.  It proves that the exact
optimized normal/defect system has a canonical all-order inner
transfer.  At this stage the remaining task was an identification,
not a search for a cancellation:

> Show that L166's transferred one-reflection Blaschke cross is the
> real mean of the corresponding marked coefficient of
> `dot(theta)/theta` from (11).

If this identification holds, (12) kills the coefficient in every
grade.  L162 separately transports the zero-reflection/equality
component.  The grade-one exception must appear as the unmatched
feedthrough/endpoint term when the compact reflected coordinate has
no defect coordinate `e_(L+1)`.

This formulation is compatible with L165: its companion bottom row
is the determinant/feedthrough part of (3), while its commutator
boundary term is the moving state-space basis part.

L171 subsequently reduces the only associated one-reflection weight
to a direct endpoint/cofactor kernel identity and proves separately
that every ordinary singular-Hessian term vanishes.  The complete
differentiated L156 norming functional still has to be identified
with the appropriate colligation coefficient; L169 alone does not
close L163.

L172 subsequently supplies the required associated identification in
the form actually needed: it writes the remaining endpoint row as a
Cauchy residue and proves by its two possible powers that the residue
vanishes.

## 6. Exact regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/rank_one_colligation_transfer.py \
  --minimum-size 2 --maximum-size 7 \
  --output experiments/rank_one_colligation_transfer_s70223.jsonl
```

The checker uses dense rational orthogonal colligations and unrelated
positive diagonal similarities.  It verifies both Stein equations,
the determinant quotient (3), and the outside resolvent identity (4)
exactly in every tested dimension.  The all-size proof is the block
determinant calculation above; the records are an independent audit.
