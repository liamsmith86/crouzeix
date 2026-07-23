# Finite-Blaschke composition of Stein metrics (2026-07-23)

## 1. Exact transfer identity

Let `B` be a finite Blaschke product of degree `d`, and let
`f_0,...,f_(d-1)` be any orthonormal basis of its model space

\[
 {\cal K}_B=H^2\ominus BH^2.
\]

For a matrix `T` with spectrum in the disk, define the positive linear
map

\[
 {\cal L}_{B,T}(X)
 =\sum_{j=0}^{d-1}f_j(T)^*Xf_j(T).                   \tag{1}
\]

Then

\[
\boxed{
 {\cal L}_{B,T}(X)-T^*{\cal L}_{B,T}(X)T
 =X-B(T)^*XB(T).}                                    \tag{2}
\]

Thus an upper contraction certificate for `B(T)` really can be lifted
to one for `T`: if the right side of (2) is positive semidefinite, so
is the left side.  If `X` is positive definite, then so is
`L_(B,T)(X)`, because a Takenaka--Malmquist basis contains an
everywhere-invertible rational factor.

This corrects only the *feasibility direction*.  It does not say that
the lifted metric has the same Euclidean condition number.

## 2. Proof

The model-space reproducing kernel identity is

\[
 \sum_{j=0}^{d-1}f_j(z)\overline {f_j(w)}
 ={1-B(z)\overline {B(w)}\over1-z\overline w}.        \tag{3}
\]

Multiplying by `1-z conjugate(w)` gives

\[
 \sum_j\{f_j(z)\overline {f_j(w)}
        -zf_j(z)\overline {wf_j(w)}\}
 =1-B(z)\overline {B(w)}.                            \tag{4}
\]

Apply the hereditary functional calculus to (4), putting the
antiholomorphic variable on the left of `X` and the holomorphic
variable on the right.  Since every `f_j(T)` commutes with `T`, the
left side becomes

\[
 {\cal L}_{B,T}(X)-T^*{\cal L}_{B,T}(X)T,
\]

and the right side becomes `X-B(T)^*XB(T)`.  This proves (2).

The argument is finite-dimensional and algebraic after the Blaschke
denominators are cleared; no boundary limit or dilation theorem is
needed.

## 3. Rank-one Gramian factorization

For a stable matrix `A` and a vector `q`, write

\[
 {\cal P}_A(q)=\sum_{n\geq0}(A^*)^nqq^*A^n,
\qquad
 {\cal P}_A(q)-A^*{\cal P}_A(q)A=qq^*.               \tag{5}
\]

Take `X=P_(B(T))(q)` in (2).  Uniqueness of the stable Stein equation
gives the exact composition law

\[
\boxed{
 {\cal P}_T(q)
 ={\cal L}_{B,T}\bigl({\cal P}_{B(T)}(q)\bigr).}      \tag{6}
\]

Consequently L118's entire optimized rank-one envelope can be
recomputed on the descended operator and then transferred through the
finite model space.  This is stronger than the scalar implication
that a contraction `T` makes `B(T)` a contraction.

## 4. The condition-number obstruction

The transfer is order preserving.  If

\[
 \alpha X_0\preceq X\preceq\beta X_0,
\]

then

\[
 \alpha{\cal L}_{B,T}(X_0)
 \preceq{\cal L}_{B,T}(X)
 \preceq\beta{\cal L}_{B,T}(X_0).                    \tag{7}
\]

But the comparison metric in (7) is `L_(B,T)(X_0)`, not automatically
the original physical coordinate Gramian.  Therefore (2) does not by
itself prove

\[
 t_*(T)\leq t_*(B(T)).
\]

For L126's central family, the numerical evidence identifies the
missing special fact much more precisely:

1. the lifted full metric reduces the exact outer space
   `span{e_0,e_k,e_(2k)}`;
2. its outer compression is exactly the size-three optimal metric;
3. every inner generalized eigenvalue lies strictly between the outer
   endpoints.

Proving those three statements for the transfer (1), or constructing
the same metric directly, would establish the conjectured exact
central identity

\[
 t_*(T_k(a,c))=t_*(T_1(a,c^k)).                      \tag{8}
\]

## 5. The dual transfer

There is a companion identity in the orientation used by L21's
trace-ratio duality.  Define

\[
 {\cal R}_{B,T}(Z)
 =\sum_{j=0}^{d-1}f_j(T)Zf_j(T)^*.                   \tag{9}
\]

The same kernel calculation gives

\[
\boxed{
 {\cal R}_{B,T}(Z)-T{\cal R}_{B,T}(Z)T^*
 =Z-B(T)ZB(T)^*.}                                    \tag{10}
\]

Thus a positive dual witness for `B(T)` lifts positively to `T` while
preserving its Hermitian Stein difference exactly.  In particular the
positive/negative trace ratio is unchanged.  This is the explicit
dual mechanism behind `t_*(B(T))<=t_*(T)`.

## 6. Deterministic audit

Run

```bash
.venv/bin/python -u experiments/blaschke_stein_composition.py \
  --output experiments/blaschke_stein_composition_s70223.jsonl
```

The checker uses deterministic complex matrices, explicit
Takenaka--Malmquist functions, and independently solved Stein
equations.  It audits (2) and (6) in dimensions three through seven
and Blaschke degrees one through four.
