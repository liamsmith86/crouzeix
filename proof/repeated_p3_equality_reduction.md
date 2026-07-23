# Equality reduction after the complete repeated-`C3` sign (2026-07-22)

## 1. The endpoint as a sum of negative terms

Use L86's maximal winner space `W`.  Its complete second endpoint is

\[
 \begin{aligned}
 \mathcal E_W={}&-16\,\operatorname{mean}\Delta(q)-8H_1^2
 -\frac{21}{4}|v|^2I\\
 &-\frac{25}{8}A_0A_0^*-8A_1A_1^*
 -\frac{50}{9}A_2A_2^*,                              \tag{1}
 \end{aligned}
\]

where

\[
 \Delta(q)=\lambda_{\max}(Q(q))I-Q(q)\succeq0.       \tag{2}
\]

Every summand in (1) is negative semidefinite.  Hence a unit vector `x` lies
in the endpoint kernel if and only if

\[
 \boxed{
 \begin{gathered}
 v=0,\qquad H_1x=0,\qquad A_j^*x=0\quad(j=0,1,2),\\
 Q(q)x=\lambda_{\max}(Q(q))x
 \quad\text{for almost every }q.
 \end{gathered}}                                      \tag{3}
\]

The last condition follows because the continuous nonnegative function
`x*Delta(q)x` has zero mean.  Conversely, (3) plainly makes `x*E_W x=0`; for a
negative semidefinite matrix this is equivalent to `E_W x=0`.

Thus (3) is an exact equality classification for the complete explicit L86
certificate, not a genericity statement.

## 2. Selected-copy normal form

Apply a constant copy unitary sending `x` to the first coordinate.  This does
not change the repeated base or the similarity value.  Reapply the L74/L84
gauge after the rotation.  Conditions (3) then say:

1. the selected copy has no quotient coupling of any generator type to the
   first-order losing sector;
2. its generator-one loop and every generator-one edge vanish;
3. the common strong single-block coordinate `v` vanishes;
4. the selected vector is a common top branch of the full second effective
   support.

Consequently the only first-order data incident to the selected branch are
the flat generator-zero/generator-two star and the common flat single-block
modes.  All other graph and losing data live in the orthogonal copy sector
and can affect the selected branch only through higher-order paths or through
collapse of the second-support gap.

This is precisely the intersection of the two previously identified boundary
models:

- L77/L80 for the flat repeated star;
- L67--L73 for the flat single-block modes and disk critical manifold.

No generator-one, common-`v`, or winner--loser direction belongs to the
remaining higher-order core.

## 3. Strictness tests and the actual blow-up variables

Equation (1) is strictly negative on a candidate vector as soon as any one of
the following has a nonzero margin:

\[
 \begin{gathered}
 x^*\operatorname{mean}\Delta x,quad \|H_1x\|^2,
 \quad |v|^2,\quad \|A_0^*x\|^2,\quad
 \|A_1^*x\|^2,\quad \|A_2^*x\|^2.                    \tag{4}
 \end{gathered}
\]

The first-order losing mean gap from L86 is an additional strict margin off
the winner sector.  Therefore a uniform repeated-block theorem needs a
weighted analysis only where **all** quantities in (4), the losing mean gap,
and the higher-order flat-star/common-mode margins tend to zero together.

This removes two possible but incorrect next steps:

- there is no need to classify further first-order quotient directions;
- a higher-order expansion on a generic graph is wasteful, because generic
  graphs already have a strict matrix-Jensen margin.

The next calculation should instead use (3) as the center stratum and resolve
the interaction between its flat star and common disk-manifold coordinates,
with the quantities in (4) treated as strong normal variables.

## 4. Audit

The decomposition (1) regenerates from

```bash
.venv/bin/python -u experiments/repeated_p3_winner_graph_sign.py
.venv/bin/python -u experiments/repeated_p3_multiwinner_gap.py
```

The equality inference uses only the elementary fact that a sum of
nonnegative quadratic forms vanishes exactly on the intersection of their
kernels.  No numerical equality threshold or unproved optimizer uniqueness is
used.
