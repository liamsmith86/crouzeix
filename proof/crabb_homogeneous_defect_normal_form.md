# Homogeneous spectral normal form for axis defects (L138, 2026-07-23)

## 1. Exact statement

Fix an L117 elliptic Crabb axis.  In coefficient coordinates write

\[
T_c=RXR^{-1},\qquad
R=K_0^{-1/2}DU,
\]

where `X=diag(x_n)`, `D=diag(c^(m/2))`, and `U` is DCT-I.
Let

\[
{\cal C}_{mn}=\frac1{1-x_mx_n},\qquad
G=\operatorname{diag}(\beta){\cal C}
  \operatorname{diag}(\beta)^* ,
\]

with L117's node defect

\[
\beta_n=\epsilon_n\frac{k'}{\operatorname{dn}(v_n\mid k^2)}.
\]

For an arbitrary rank-one forcing vector `q`, put

\[
y=R^*q,\qquad h_n=\frac{y_n}{\beta_n},\qquad
A_h=U\operatorname{diag}(h)U^*.
\]

Then its Stein Gramian

\[
M(q)=\sum_{\nu\ge0}(T_c^*)^\nu qq^*T_c^\nu
\]

has exact generalized physical metric

\[
\boxed{
K_0^{-1/2}M(q)K_0^{-1/2}
=D^{-1}A_hWA_h^*D^{-1},
\qquad W=UGU^*.
}                                                     \tag{1}
\]

Thus homogeneous defect variation is exactly multiplication by a
sampled scalar function in the DCT node basis.  A constant `h`
only rescales the metric and is the null defect-scale direction.
Equation (1) restores the homogeneous coordinate deliberately fixed
away by `x_0=0` in the exact Hessian engine.

## 2. Proof

The spectral Stein sum is entrywise:

\[
\begin{aligned}
M(q)
&=R^{-*}G_yR^{-1},\\
(G_y)_{mn}
&=\frac{y_m\overline y_n}{1-x_mx_n}.
\end{aligned}                                        \tag{2}
\]

Since `y=diag(h) beta`,

\[
G_y=\operatorname{diag}(h)G\operatorname{diag}(h)^*.
\]

Conjugating by `U` gives

\[
UG_yU^*=A_hWA_h^*.
\]

Finally

\[
R^{-1}=U^*D^{-1}K_0^{1/2}.
\]

Substitution in (2), followed by congruence by `K_0^{-1/2`,
is exactly (1).

No optimization, asymptotic expansion, or finite-size inference is
used.

## 3. Exact condition Hessian

Normalize the axis so that

\[
P_0=D^{-1}WD^{-1}
=\operatorname{diag}(p_0,\ldots,p_L),\qquad p_0=1.
\]

For `h=1+a s`, define

\[
H_s=U\operatorname{diag}(s)U^*,\qquad
B=D^{-1}H_sD.
\]

Equation (1) becomes the exact quadratic congruence

\[
\boxed{
P(a)=(I+aB)P_0(I+aB)^*.
}                                                     \tag{3}
\]

The DCT endpoint rows have identical squared entries, so

\[
(H_s)_{00}=(H_s)_{LL}.
\]

Consequently the lower and upper generalized eigenvalues in (3)
have the same relative first derivative; the condition number is
stationary in every homogeneous defect direction.

The complete second derivative is finite.  Put

\[
\begin{aligned}
E_{ij}&=p_jB_{ij}+p_i\overline{B_{ji}},\\
q_i&=\sum_jp_j|B_{ij}|^2
 +\sum_{j\ne i}\frac{|E_{ij}|^2}{p_i-p_j}.
\end{aligned}                                        \tag{4}
\]

Since the axis endpoint levels are simple, ordinary Hermitian
eigenvalue perturbation gives

\[
\boxed{
Q_{\rm def}(s)=q_L-p_Lq_0.
}                                                     \tag{5}
\]

The products of the common first endpoint shift cancel in the ratio.
Equations (3)--(5) are an exact closed formula for the universal pure
defect Hessian; they replace repeated Stein propagation by DCT
multiplication and two endpoint Schur sums.

## 4. Toeplitz-plus-Hankel structure

Product-to-sum gives

\[
(H_s)_{mn}
=\sum_{\nu=0}^L U_{m\nu}U_{n\nu}s_\nu
\]

as a weighted sum of the DCT coefficients at indices `m-n` and
`m+n`.  Hence `H_s` is Toeplitz plus Hankel, with the standard
endpoint weights.  Formula (4) couples only its two endpoint rows
after the diagonal gauge `D`.

This proves the structural part of A98 that was previously only
described as “restore homogeneous defect scale.”  It does **not**
yet prove that the LDL Newton edge of (5) is

\[
\frac{1+c{\cal S}^2}{1-c{\cal S}^2},
\]

nor does it insert L131's operator/Faber forcing.  Those are now
explicit coefficient calculations in (4), rather than an unknown
Stein-kernel transformation.

## 5. Independent regression

The identity can be checked directly by constructing both sides of
(1) from arbitrary complex node multipliers.  The persisted checker

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_homogeneous_defect_normal_form.py \
  --output experiments/crabb_homogeneous_defect_normal_form_s70223.jsonl
```

compares the spectral formula with a direct Stein solve and checks
(4)--(5) by centered finite differences.  The algebra above is the
proof.
