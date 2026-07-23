# The flat Crabb endpoint selection rule (2026-07-23)

## 1. Result

Let `p=L+1>=3`,

\[
 A_c=C_p+cC_p^*,\qquad
 P=P_{L,c},
\]

where `P` is L119--L120's Dickson descent polynomial.  Define the
endpoint functional

\[
 {\cal F}_c(Y)
 =(DP(A_c)[Y])_{L0}-c^L(DP(A_c)[Y])_{0L}.             \tag{1}
\]

L120 proves `F_c(Y)=O(c^floor(L/2))` for arbitrary matrix directions.
On the disk-flat space

\[
 {\cal D}_p=\ker e_p\cap T_{C_p}{\cal C}_p
\]

from L115, there is one additional exact cancellation:

\[
 \boxed{
 Y\in{\cal D}_p
 \quad\Longrightarrow\quad
 {\cal F}_c(Y)
 =O\!\left(c^{\lfloor L/2\rfloor+1}\right)\|Y\|.}     \tag{2}
\]

Combining (2) with L120's exact chain rule gives

\[
 \boxed{
 D\Gamma_p(A_c)[Y]
 =O\!\left(c^{L+\lfloor L/2\rfloor+1}\right)\|Y\|,
 \qquad Y\in\ker e_p.}                                \tag{3}
\]

In particular (3) holds on all of L115's disk-flat quotient
directions.  This is the strict weighted exponent needed for the
linear part of the remaining disk/ellipse merger.

This note does **not** prove that merger.  It removes its dangerous
linear monomial; a quantitative disk-flat anchor and control of the
higher mixed terms are still required.

## 2. Exact endpoint formula

Put

\[
 d_0=d_L=2^{-1/2},\qquad d_j=1\quad(0<j<L),
\]

and assign the circle grade

\[
 q(j,m)=m-j-1
\]

to the entry `Y_jm`.  The root-of-unity calculation in L120 equation
(12) gives the full coefficient formula

\[
\begin{aligned}
 {\cal F}_c(Y)
={}&2\!\!\sum_{\substack{j,m\\q(j,m)\leq-2\\q(j,m)\ {\rm even}}}
 d_jd_m c^{\,L+q(j,m)/2}Y_{jm}\\
 &-2\!\!\sum_{\substack{j,m\\q(j,m)\geq0\\q(j,m)\ {\rm even}}}
 d_jd_m c^{\,L+q(j,m)/2}Y_{jm}.                       \tag{4}
\end{aligned}
\]

For completeness, L120's differential quadrature retains a matrix
entry precisely when

\[
 q(j,m)+L-2s=\pm L,\qquad 0\leq s<L.
\]

The minus sign gives `s=L+q/2` and the first sum in (4).  The plus sign
gives `s=q/2`; its contribution is multiplied by the `-c^L` in (1),
giving the second sum.  This also proves (4) directly from the finite
Dickson recurrence, without an asymptotic argument.

## 3. The first coefficient is an L65 bottom mode

The smallest physical grade is `-L-1`.

### Odd `L`

Here `-L-1` is even, so the first possible power is `(L-1)/2`.  There
is only one entry of this grade.  Formula (4) gives

\[
 [c^{(L-1)/2}]{\cal F}_c(Y)=Y_{L0}.                   \tag{5}
\]

L65 equation (20) is a strictly negative multiple of `|Y_L0|^2`.
Consequently `Y in ker e_p` implies that (5) vanishes.

### Even `L`

The lowest even physical grade is `-L`, so the first possible power is
`L/2`.  Its two entries give

\[
 [c^{L/2}]{\cal F}_c(Y)
 =\sqrt2\,(Y_{L-1,0}+Y_{L,1}).                        \tag{6}
\]

L65 equation (19) is a strictly negative multiple of
`|Y_(L-1,0)+Y_(L,1)|^2`.  For `L>=4`, equality therefore kills (6)
directly.  For `L=2`, that coefficient vanishes; instead, L115's
circular-tangent equation (5) kills exactly the same combination.  The
residual `p=3` merger is also proved independently in L73.

Thus the first unrestricted coefficient vanishes on the relevant
flat quotient in every size.  The next possible integer power is
`floor(L/2)+1`, proving (2).

For `L>=3` (equivalently `p>=4`), the argument actually applies to the
full L65 equality space, before intersecting it with the circular-range
tangent.  Only the exceptional `p=3` case needs the disk-tangent
constraint.

## 4. Consequence for the weighted merger

L120 proves

\[
 D\Gamma_p(A_c)[Y]
 ={\tau'(c^L)\over2}\operatorname{Re}{\cal F}_c(Y),
 \qquad
 {\tau'(c^L)\over2}=-16c^L+O(c^{3L}).
\]

Equation (3) follows immediately.  Write

\[
 m=\lfloor L/2\rfloor+1.
\]

If the disk-flat restriction supplies a coercive quartic margin
`-a||d||^4`, the linear mixed term has size

\[
 O(\|d\|c^{L+m}).
\]

Young's inequality turns this into

\[
 \varepsilon\|d\|^4+C_\varepsilon c^{4(L+m)/3}.       \tag{7}
\]

The strict selection `m>L/2` gives

\[
 {4(L+m)\over3}>2L,                                   \tag{8}
\]

so the second term in (7) is smaller than L117's axis margin
`-16c^(2L)`.  Without the extra cancellation, equality occurs when
`L` is even and this absorption would lose its strict power advantage.

Equations (7)--(8) are a conditional reduction, not a proof that the
needed quartic disk margin exists.

**Subsequent correction (L122).**  The premise is false: on the exact
Toeplitz disk chart the quartic is only
`-32(||z||^4-|z^T Jz|^2)` and has a large phase-palindromic null cone.
Thus (7)--(8) remain a valid conditional exponent calculation but
cannot close the full merger.

## 5. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_flat_endpoint_selection.py \
  --output experiments/crabb_flat_endpoint_selection_s70223.jsonl
```

For `p=3,...,16`, the checker:

1. differentiates the Dickson recurrence with every matrix entry
   independent;
2. verifies the complete coefficient formula (4);
3. identifies the first coefficient exactly with (5) or (6); and
4. imposes the corresponding L65 bottom equality (or the `p=3`
   circular-tangent condition) and verifies exact divisibility by
   `c^(floor(L/2)+1)`.

All comparisons are symbolic SymPy equalities.
