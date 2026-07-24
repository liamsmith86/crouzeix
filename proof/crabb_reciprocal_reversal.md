# Reciprocal-reversal symmetry of the optimized Crabb metric

## 1. Statement (L164, 2026-07-23)

Let `J` reverse `p` real coordinates and let

\[
T=JT^TJ
\]

be a stable real persymmetric matrix near the Crabb block.  On L118's
normalized local defect chart, let `q_*(T)` be the unique stationary
rank-one defect and let

\[
P=P(T,q_*)=\sum_{n\ge0}(T^T)^nq_*q_*^TT^n .
\]

Then there is a positive analytic scalar `alpha(T)` such that

\[
\boxed{JP^{-1}J=\alpha P.}                           \tag{1}
\]

In particular,

\[
\boxed{\alpha=(\det P)^{-2/p}}                       \tag{2}
\]

and the eigenpairs occur in reversed reciprocal pairs.  If
`Pu=lambda u`, then

\[
PJu={1\over\alpha\lambda}Ju.                         \tag{3}
\]

The phase-gauged complex version follows by replacing reversal with
the corresponding antiunitary conjugate reversal.

Equation (1) applies to the optimized metric on every real
phase-palindromic reflected path used in L160--L163.  It is an exact
identity of the analytic branch, not a finite-series conjecture.

It does **not** by itself prove L163.  Reciprocal pairing constrains
the sum of the two relative endpoint derivatives, whereas L163 needs
their difference to vanish.  Its value is that the remaining
covariant recurrence may be written for one endpoint, or equivalently
for the adjoint condition gradient, with the other endpoint supplied
by (1).

## 2. Rank-one inverse-Stein duality

Start with any positive Stein solution

\[
P-T^TPT=qq^T.                                        \tag{4}
\]

Put

\[
C=P^{1/2}TP^{-1/2},\qquad a=P^{-1/2}q .
\]

Then

\[
I-C^TC=aa^T.                                         \tag{5}
\]

Thus `C` is a square contraction with rank-one right defect.  The
matrices `C^TC` and `CC^T` have the same eigenvalues, so

\[
I-CC^T=bb^T                                         \tag{6}
\]

for a vector `b`, with the same sole nonzero eigenvalue as (5).
Consequently

\[
\boxed{
P^{-1}-TP^{-1}T^T=rr^T,\qquad r=P^{-1/2}b .
}                                                     \tag{7}
\]

This also proves positivity and rank one in (7); no matrix-entry
recurrence is required.

Now set

\[
Q=JP^{-1}J.
\]

Persymmetry gives `T^TJ=JT` and `JT^T=TJ`, so (7) implies

\[
\begin{aligned}
Q-T^TQT
 &=J(P^{-1}-TP^{-1}T^T)J\\
 &=(Jr)(Jr)^T.                                      \tag{8}
\end{aligned}
\]

Therefore reversal of the inverse Gramian is another rank-one Stein
Gramian for the **same** operator.  Its condition number is exactly
that of `P`.

## 3. The optimized branch is a fixed point

Near the Crabb point, the last component of `r` in (7) is nonzero.
Choose its sign analytically, reverse it, and normalize its first
component to one.  Equation (8) defines a local analytic map

\[
{\cal D}_T:q\longmapsto q^\vee                       \tag{9}
\]

on the projective defect chart.  Applying the construction twice
returns the original Gramian up to its irrelevant positive scale.
Thus `D_T` is a local analytic involution.

If

\[
F_T(q)=\kappa(P(T,q)),
\]

then (8) gives

\[
F_T({\cal D}_Tq)=F_T(q).                             \tag{10}
\]

Hence `D_T` carries critical defects to critical defects.  L118's
positive defect Hessian gives a unique critical defect near `e_0`,
so

\[
{\cal D}_Tq_*=q_*.
\]

The two Stein Gramians in (8) must therefore differ only by the
positive defect normalization scale.  This proves (1).  Taking
determinants in (1) gives (2), and applying (1) to an eigenvector
gives (3).

## 4. Why the first numerical audit looked imperfect

Suppose the stationary defect series is known only through
`epsilon^j`.  Its Gramian can satisfy (1) only through the same
degree: the omitted defect coefficient at degree `j+1` enters the
Stein forcing linearly.  In the grade-two size-five calculation:

* the three-jet defect satisfies (1) through degree three and first
  leaves a residual at degree four;
* after solving the fourth defect jet, that entire degree-four
  residual vanishes exactly and the first residual moves to degree
  five.

Thus evaluated residuals of size `O(epsilon^(j+1))` are truncation
errors, not evidence against (1), and must not be mistaken for
roundoff-level confirmation through twice the face weight.

## 5. Regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_reciprocal_reversal.py \
  --output experiments/crabb_reciprocal_reversal_s70223.jsonl
```

The checker has two independent exact parts.

1.  An arbitrary rational defect for the unweighted persymmetric
    nilpotent shift verifies that (7) has rank one and that (8)
    solves the reversed Stein equation.
2.  The size-five grade-two reflected path solves four stationary
    defect jets, verifies (1) coefficientwise through degree four,
    and confirms a nonzero residual first at the omitted fifth jet.

The all-size proof is Sections 2--3; the finite calculation audits the
normalization and the important truncation boundary.

## 6. Remaining L163 gate

Let `H` be the derivative of the optimized Gramian in a circular
normal direction.  Differentiating (1) relates the reversed endpoint
response to `H` and to `D alpha[H]`, but it does not force
`D log kappa[H]=0`.  The required statement is still the coefficient
of weight `k+1` in the adjoint envelope gradient:

\[
D\log\kappa(P)[H]
={v_+^THv_+\over\lambda_+}
 -{v_-^THv_-\over\lambda_-}.                         \tag{11}
\]

Equation (3) removes one independent endpoint series from (11).  The
next proof step must show that the complete inverse-Riemann/Stein
normal forcing has zero coefficient in (11) for grades `k>=2`.
The grade-one exception in L160 remains compatible with (1).
