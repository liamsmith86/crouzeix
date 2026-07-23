# Dual-witness construction of the central metric (2026-07-23)

## 1. Status

This note reduces the remaining central upper lift to one exact
fiber-orthogonality identity.  That identity is numerically exact but
is **not yet proved here**, so this note is not a completed lemma.

Let

\[
 T=T_k(a,c),\qquad U=B_{k,c}(T),\qquad
 V=\operatorname{span}\{e_0,e_k,e_{2k}\}
\]

in physical coordinates.  L126 and L129 give

\[
 U=U_3\oplus U_{\rm in},\qquad
 t_*(U)=t_*(U_3)=:t,                                 \tag{1}
\]

where `U_3=T_1(a,c^k)`.

## 2. Lift the outer dual, not the outer metric

Choose an optimal L21 dual witness `Z_3>=0` for `U_3`.  On the local
size-three branch it has rank two and one-dimensional kernel spanned
by the defect vector `q_3`.  Embed it by zero on the inactive summand:

\[
 Z_U=\iota Z_3\iota^*.
\]

For a model-space basis `f_0,...,f_(k-1)` of `K_(B_(k,c))`, use L127's
dual transfer

\[
 Z_T=\sum_{j=0}^{k-1}f_j(T)Z_Uf_j(T)^*.              \tag{2}
\]

Then

\[
 Z_T-TZ_TT^*=Z_U-UZ_UU^*.                            \tag{3}
\]

The right side is supported on `V` and has the same one-positive,
one-negative trace ratio `t` as the size-three witness.  At the apex,
the `2k` transferred range vectors are independent.  Analyticity then
gives

\[
 \operatorname{rank}Z_T=2k
\]

locally.  Let `q` span `ker Z_T`.

## 3. Kernel alignment is automatic

Put

\[
 q_j=f_j(T)^*q,\qquad
 q_j=\alpha_j\iota q_3+w_j,\quad w_j\in V^\perp.      \tag{4}
\]

This decomposition is forced by positivity, not conjectured.  Indeed,

\[
 0=q^*Z_Tq
 =\sum_j\|Z_3^{1/2}\iota^*f_j(T)^*q\|^2,
\]

so every outer component `iota^*q_j` lies in `ker Z_3`.

The Wold/Stein composition formula gives

\[
 {\cal P}_T(q)=\sum_{j=0}^{k-1}{\cal P}_U(q_j).       \tag{5}
\]

Therefore its outer compression is a scalar copy of the size-three
metric once the outer/inner forcing decouples.

## 4. The sole missing algebraic identity

The exact target is

\[
\boxed{
 \sum_{j=0}^{k-1}\alpha_jw_j^*=0.}                   \tag{CF}
\]

Equivalently, if

\[
 Q=\sum_jq_jq_j^*,
\]

then

\[
 \Pi_VQ(I-\Pi_V)=0.                                   \tag{6}
\]

If (CF) holds, the cross Stein equation for the reducing operator
`U_3+U_in` has zero forcing and a unique stable solution.  Hence

\[
\begin{aligned}
 \Pi_V{\cal P}_T(q)(I-\Pi_V)&=0,\\
 \iota^*{\cal P}_T(q)\iota
 &=\left(\sum_j|\alpha_j|^2\right)
   {\cal P}_{U_3}(q_3).                              \tag{7}
\end{aligned}
\]

At the disk apex the remaining generalized metric levels are strictly
between the lower and upper size-three levels (L123).  Rank, metric,
and generalized eigenvalues vary analytically.  Thus (CF), together
with L129, proves locally

\[
\boxed{t_*(T_k(a,c))=t_*(T_1(a,c^k)).}               \tag{8}
\]

In particular the complete central A84 Hessian is exactly the
size-three Hessian at `r=c^k`, so its first term is
`-64a^2c^(2k)`.

## 5. Why (CF) should be tractable

The identity is stronger than a fit:

* it holds to binary64 residual `3.5e-13` on the deterministic
  degree-`2,...,5` grid;
* it holds for unrelated deterministic rank-two positive matrices
  `Z_3`, not only the optimal dual witness; and
* it is independent of the chosen orthonormal model-space basis.

Use a Clark basis for `K_(B_(k,c))`.  In the exterior ellipse
coordinate its nodes are the `k` roots

\[
 \zeta_j^k=\xi.
\]

The central Toeplitz data depend only on `zeta^k`.  L120's
root-of-unity filter should therefore kill every cross grade between
`{0,k,2k}` and the remaining coefficient grades.  The next proof must
write (6) in that Clark basis and apply this two-variable version of
L120's conditional expectation.

## 6. Deterministic regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_central_dual_lift.py \
  --output experiments/crabb_central_dual_lift_s70223.jsonl
```

The checker constructs the dual lift, its kernel-generated Stein
metric, and the forcing (6) independently of the full rank-one
optimizer.  It compares the resulting condition with the size-three
value and audits the metric splitting.
