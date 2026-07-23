# The descended central inner fibers are uniformly inactive (2026-07-23)

## 1. Result

Use L126--L128 with `k>=2`, and let

\[
 U_k(a,c)=B_{k,c}(T_k(a,c)).
\]

For all sufficiently small `|a|+c`, uniformly in `k`, the orthogonal
outer/inner decomposition from L126 satisfies

\[
\boxed{
 t_*(U_k(a,c))
 =t_*(T_1(a,c^k)).}                                  \tag{1}
\]

The right side is the exact size-three outer block.  Thus every
descended inner fiber is separated from the active similarity
endpoints by a fixed gap.  This proves which block is active after
descent; it does not yet lift the bound back from `B(T)` to `T`.

## 2. Uniform inner contraction metric

Put `r=c^k`.  L128 writes the inner polynomial and coordinate Gramian
as

\[
 P_{\rm in}
 =\begin{pmatrix}H&I\\rI&H\end{pmatrix}
 +2a\begin{pmatrix}-rI&H\\H&-I\end{pmatrix},
\qquad
 K_{\rm in}=\begin{pmatrix}I&2aI\\2aI&I\end{pmatrix}, \tag{2}
\]

where

\[
 H^2=rI,\qquad \|H\|=c,\qquad r\leq c^2.             \tag{3}
\]

The descended disk operator is

\[
 U_{\rm in}=\phi_r(P_{\rm in}).                       \tag{4}
\]

At `a=c=0`,

\[
 U_{\rm in}=N:=
 \begin{pmatrix}0&I\\0&0\end{pmatrix},\qquad
 K_{\rm in}=I.
\]

Take

\[
 M_{\rm in}^{(0)}
 =\begin{pmatrix}\frac12I&0\\0&I\end{pmatrix}.        \tag{5}
\]

Then

\[
 M_{\rm in}^{(0)}
 -N^*M_{\rm in}^{(0)}N
 =\frac12I.                                           \tag{6}
\]

This is a strictly feasible contraction metric of generalized
condition square two, with a dimension-independent Stein gap.

## 3. Uniform persistence

Equations (2)--(3) give

\[
 \|P_{\rm in}-N\|=O(|a|+c)
\]

with an absolute constant independent of `k`.  The normalized ellipse
map obeys `phi_r(z)=z+O(r)` uniformly on a fixed neighbourhood of the
relevant spectra, so

\[
 \|U_{\rm in}-N\|=O(|a|+c)                            \tag{7}
\]

uniformly as well.  By (6), the fixed metric (5) therefore remains
Stein-feasible when `|a|+c` is sufficiently small.

The eigenvalues of `K_in` lie in

\[
 [1-2|a|,1+2|a|].
\]

Consequently the generalized condition of (5) relative to `K_in` is
at most

\[
 2\,{1+2|a|\over1-2|a|}<3                            \tag{8}
\]

after shrinking the neighbourhood once.  Hence

\[
 t_*(U_{\rm in})<3.                                   \tag{9}
\]

On the other hand, the size-three outer value is continuous at
`(a,r)=(0,0)`, where it equals four.  Shrinking the same neighbourhood
gives

\[
 t_*(T_1(a,r))>3.                                     \tag{10}
\]

L126 makes the outer/inner splitting orthogonal in the physical
coordinate Gramian.  The similarity optimum of a direct sum is the
maximum of the optima of its summands.  Equations (9)--(10) prove (1).

## 4. Scope

This removes a possible failure mode from the central campaign:
no hidden inner block can overtake the size-three block in the
singular `a,c->0` regime.  Combined with L127, the remaining upper
problem is solely quantitative:

> transfer the active outer metric through the degree-`k` model space
> without increasing its two endpoint levels at order `a^2c^(2k)`.

The inactive block itself no longer needs an exact elliptic SDP.

**Subsequent closure (L130).**  The outside critical factor of the
finite Blaschke product gives an explicit lift whose metric reduces
the outer space.  The strict apex metric gap then traps its inner
levels, proving the exact local central identity.
