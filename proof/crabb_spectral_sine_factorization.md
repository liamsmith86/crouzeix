# Spectral sine factorization of equality tangents (2026-07-23)

## 1. Exact statement

Put `p=L+1`, and use the coefficient coordinates of L123.  Let

\[
S_0=C+cJCJ,\qquad
K_0=\operatorname{diag}(1/2,1,\ldots,1,1/2).
\]

Set

\[
\epsilon_0=\epsilon_L=2^{-1/2},\quad
\epsilon_m=1\ (0<m<L),\quad
\theta_n={n\pi\over L},
\]

and let

\[
U_{mn}=\sqrt{2/L}\,\epsilon_m\epsilon_n\cos(m\theta_n),
\qquad
D=\operatorname{diag}(1,c^{1/2},\ldots,c^{L/2}).
\]

L117's DCT-I diagonalization, transported to coefficient coordinates,
is

\[
\boxed{
R=K_0^{-1/2}DU,\qquad
S_0R=R\operatorname{diag}(2\sqrt c\cos\theta_n).
}                                                     \tag{1}
\]

For a coefficient grade `1<=j<L`, L131's companion row is

\[
v_j^*=e_{j+1}^*-ce_{j-1}^*,\qquad
H_j=2(e_0-e_L)v_j^*.
\]

Its right spectral factor is exactly

\[
\boxed{
(v_j^*R)_n
=-2\sqrt{2/L}\,\epsilon_n c^{(j+1)/2}
  \sin(j\theta_n)\sin\theta_n.
}                                                     \tag{2}
\]

Thus distinct companion grades are discrete sine modes after removing
the one common endpoint factor `sin(theta_n)`.  In particular,

\[
\sum_{n=1}^{L-1}
\sin(j\theta_n)\sin(k\theta_n)
={L\over2}\delta_{jk}.                               \tag{3}
\]

Equations (1)--(3) are exact for every `L>=2` and `0<c<1`; they are
not asymptotic or finite-size observations.

## 2. Proof

The physical Crabb axis is

\[
K_0^{1/2}S_0K_0^{-1/2}
=D\{\sqrt c(C_p+C_p^*)\}D^{-1}.
\]

The symmetric matrix in braces is diagonalized by the DCT-I matrix
`U`, with eigenvalues `2sqrt(c)cos(theta_n)`.  This proves (1).

The endpoint factors in `K_0^{-1/2}` cancel the row factors
`\epsilon_m` in `U`.  Therefore

\[
e_m^*R
=\sqrt{2/L}\,c^{m/2}
  (\epsilon_n\cos(m\theta_n))_{n=0}^L.               \tag{4}
\]

Subtract the instances `m=j+1` and `m=j-1`, using the coefficient `c`
in `v_j`.  Both acquire the common power `c^((j+1)/2)`, and

\[
\cos((j+1)\theta)-\cos((j-1)\theta)
=-2\sin(j\theta)\sin\theta.
\]

This is (2).  Equation (3) is the standard finite geometric-sum
identity for the DST-I matrix, or follows by applying
`2sin(jt)sin(kt)=cos((j-k)t)-cos((j+k)t)` and summing over the
`L-1` interior roots.

For a general coefficient vector,

\[
\left(\sum_j u_jv_j^*R\right)_n
=-2\sqrt{2/L}\,\epsilon_n\sin\theta_n
  \sum_j u_jc^{(j+1)/2}\sin(j\theta_n).              \tag{5}
\]

Hence there is no terminal alias in this spectral leg: the
coefficient-to-tangent map is an invertible scaled DST-I transform.

## 3. Consequence and remaining gate

This supplies the missing coordinate bridge behind A95.  L117's axis
Stein kernel is Toeplitz-plus-Hankel and DCT-I diagonal; (2) shows that
the equality tangents are the complementary sine modes.  The exact
vanishing of distinct-grade optimized Hessian jets is therefore
consistent with a full spectral mode decomposition, rather than an
accidental first-face cancellation.

This lemma does **not** prove that decomposition.  L118 minimizes a
quadratic in the defect tangent.  In spectral coordinates the raw
quadratic and its defect Schur correction both couple the two spectral
legs.  The next required identity is that, after the common
`sin(theta_n)` factors in (2) are removed, the Schur kernel is itself
Toeplitz-plus-Hankel (equivalently DST-I diagonal).  Its sine
eigenvalues must then be expanded at small `c`; the target first term
in grade `j` is `-64c^(2j)`.

There is a useful exact calibration for that calculation.  L117 writes
the DCT-I eigenvalues of the axis Szegő kernel as

\[
w_m={S_m\over S_0},\qquad
S_m=\sum_{r\in\mathbb Z}\operatorname{sech}
((m+2Lr)(-\log c)).
\]

For the relevant independent grades
`1<=j<=floor(L/2)`, the lower cosine mode `m=j-1` in (2) satisfies

\[
{w_{j-1}\over\epsilon_{j-1}^2}
=2c^{j-1}+O(c^{j+1}).                               \tag{6}
\]

The endpoint case `j=1` is included because
`w_0=1` and `epsilon_0²=1/2`.  Thus the single principal-symbol
identity still to prove can be normalized as

\[
-32c^{j+1}{w_{j-1}\over\epsilon_{j-1}^2}
=-64c^{2j}+O(c^{2j+2}).                             \tag{7}
\]

The multiplier `-32` is already calibrated by L133 at `j=1`.
What is not yet proved is that the full defect Schur operation
transports that same principal multiplier to every DST-I mode.  Once
that mode-translation statement is established, (6)--(7) give both
the nondivisor diagonal coefficient and distinct-grade orthogonality
on the first face at once.

This formulation is stronger and safer than invoking endpoint
Parseval alone: it explicitly includes the metric optimizer whose
omission falsified L131's naive defect transport.

## 4. High-precision regression

Run

```bash
.venv/bin/python -u experiments/crabb_spectral_sine_modes.py \
  --output experiments/crabb_spectral_sine_modes_s70223.jsonl
```

The checker independently constructs `S_0`, `R`, and every tangent
row.  At 80 decimal digits it audits (1)--(3) in sizes 3 through 20
and at `c=.05,.2,.6`.  The finite grid is a regression for the direct
proof above, not its justification.
