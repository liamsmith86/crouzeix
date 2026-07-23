# The exact arbitrary-size elliptic Crabb axis (2026-07-22)

## 1. Theorem and status

For the `p x p` Crabb block `C_p`, put

\[
 A_c=C_p+cC_p^*,\qquad 0<c<1,\qquad L=p-1.
\]

This loses no complex soft-axis directions: diagonal unitary conjugation
followed by a scalar rotation sends `C_p+zC_p*` to
`C_p+|z|C_p*`.

Let `phi_c` be the centered Riemann map from the interior of the
numerical-range ellipse of `A_c` to the disk, and set `T_c=phi_c(A_c)`.
Then, for every `p>=3`,

\[
 \boxed{t_*(T_c)={k(c^{2L})\over c^L}<4,}               \tag{1}
\]

where `k(q)=(theta_2(q)/theta_3(q))^2` is Jacobi's elliptic
modulus as a function of the nome.  Here `t_*` is L21's optimal
similarity square.  Consequently the numerical range of every `A_c`
is a complete spectral set with constant strictly smaller than two.

L116 proved the matching lower bound in (1).  The all-size diagonal
metric below proves the upper bound and closes the gap explicitly
recorded on page 46 of Kenan Li's 2021 thesis.  This is a theorem on
the one-parameter fixed-weight elliptic Crabb axis, not on arbitrary
matrices or even arbitrary elliptic weighted shifts.

## 2. Exact elliptic functional calculus

Let

\[
 D=\operatorname{diag}(1,c^{1/2},\ldots,c^{L/2}),\qquad
 H_c=\sqrt c\,(C_p+C_p^*).
\]

Then

\[
 A_c=D H_cD^{-1},\qquad
 \sigma(H_c)=\{2\sqrt c\cos(j\pi/L):0\le j\le L\}.       \tag{2}
\]

The numerical range of `A_c` is the ellipse with semiaxes `1+c` and
`1-c` and foci `+-2sqrt(c)`.  Put

\[
 q=c^2,\qquad k=k(q),\qquad K=K(k^2),
\]

and define the endpoint factors and orthogonal DCT-I matrix

\[
 e_0=e_L=2^{-1/2},\quad e_j=1\ (0<j<L),\qquad
 U_{mj}=\sqrt{2/L}\,e_me_j\cos(mj\pi/L).                \tag{3}
\]

At `v_j=2Kj/L`, the disk nodes are

\[
 x_j=\sqrt{k}\operatorname{cd}(v_j\mid k^2)
     =\sqrt{k}\operatorname{sn}(K-v_j\mid k^2).         \tag{4}
\]

Thus, with `X=diag(x_j)` and `F=UXU^T`,

\[
 T_c=D F D^{-1}.                                       \tag{5}
\]

This is an exact finite functional calculus; no boundary quadrature or
fitted conformal map is used.

## 3. Exact Chebyshev--Blaschke lower bound

The normalized degree-`L` Chebyshev--Blaschke product for the interval
`[-sqrt(k(q)),sqrt(k(q))]` satisfies, up to a common sign,

\[
 B_L(x_j)=(-1)^j\sqrt{k(q^L)}.                         \tag{6}
\]

This is Wang's Jacobi-`cd` multiplication formula; Ng--Tsang prove the
least-deviation property, and Müller-Hermes--Szehr record the equivalent
extremal value.  The DCT alternation identity is

\[
 U\operatorname{diag}((-1)^j)U^T=\pm J,                \tag{7}
\]

where `J` reverses coordinates.  Therefore

\[
 B_L(T_c)=\pm\sqrt{k(c^{2L})}\,D J D^{-1},\qquad
 \|B_L(T_c)\|^2={k(c^{2L})\over c^L}.                  \tag{8}
\]

If `S T_c S^{-1}` is a contraction, von Neumann's inequality applied
to `B_L` gives `||B_L(T_c)||<=cond(S)`.  Taking the infimum proves

\[
 t_*(T_c)\ge {k(c^{2L})\over c^L}.                     \tag{9}
\]

## 4. The elliptic Szegő kernel

Write `k'=sqrt(1-k^2)` and set

\[
 \delta_j={k'\over\operatorname{dn}(v_j\mid k^2)},\qquad
 \beta_j=e_j\delta_j,\qquad
 G_{ij}={\beta_i\beta_j\over1-x_ix_j}.                 \tag{10}
\]

By construction,

\[
 G-XGX=\beta\beta^T.                                  \tag{11}
\]

The point is that `G` has an exact DCT diagonalization.  Define

\[
 g(v)=\operatorname{dn}(v\mid k^2)+
      k\operatorname{cn}(v\mid k^2).
\]

Jacobi's addition formulas give

\[
 {\delta(v)\delta(w)\over
  1-k\operatorname{cd}(v)\operatorname{cd}(w)}
 ={g(v-w)+g(v+w)\over2}.                              \tag{12}
\]

Indeed, the right side equals

\[
 {\operatorname{dn}v\operatorname{dn}w+
   k\operatorname{cn}v\operatorname{cn}w
  \over1-k^2\operatorname{sn}^2v\operatorname{sn}^2w},
\]

and rationalizing uses the elementary identity

\[
 (\operatorname{dn}v\operatorname{dn}w)^2
 -k^2(\operatorname{cn}v\operatorname{cn}w)^2
 =k'^2(1-k^2\operatorname{sn}^2v\operatorname{sn}^2w).
\]

Consequently

\[
 G_{ij}={e_ie_j\over2}\{g(v_i-v_j)+g(v_i+v_j)\}.       \tag{13}
\]

This is the standard Toeplitz-plus-Hankel kernel diagonalized by
DCT-I.  If `hat g_m` is the length-`2L` discrete Fourier transform of
the samples `g(2Kr/L)`, then

\[
 UGU^T=\operatorname{diag}(\widehat g_m/2)_{m=0}^L.    \tag{14}
\]

Thus (11) already supplies a positive rank-one Stein defect for a
diagonal metric.  It remains only to identify and bound its entries.

## 5. Periodized-sech and Poisson formulas

Put `ell=-log(c)`.  Combining the standard Fourier series of `dn` and
`cn`, using that their nome is `q=c^2`, gives

\[
 g(v)={\pi\over2K}\left[
  1+2\sum_{r\ge1}\operatorname{sech}(r\ell)
       \cos {r\pi v\over2K}\right].                   \tag{15}
\]

Sampling at `v_j=2Kj/L` and aliasing the length-`2L` Fourier series
yields

\[
 \widehat g_m={\pi L\over K}S_m,\qquad
 S_m=\sum_{n\in\mathbb Z}
       \operatorname{sech}((m+2Ln)\ell).               \tag{16}
\]

Normalize

\[
 w_m={S_m\over S_0},\qquad W=\operatorname{diag}(w_m). \tag{17}
\]

Equations (11), (14), and (16) show that, for a vector `b`,

\[
 W-FWF=bb^T\succeq0.                                  \tag{18}
\]

There is also a closed Jacobi form for the weights.  Poisson summation
with

\[
 \int_{\mathbb R}\operatorname{sech}(\ell x)e^{-i\xi x}\,dx
 ={\pi\over\ell}\operatorname{sech}{\pi\xi\over2\ell}
\]

gives, for real `x`,

\[
 S(x)={\pi\over2L\ell}\left[
  1+2\sum_{r\ge1}
   \operatorname{sech}{\pi^2r\over2L\ell}
   \cos{\pi r x\over L}\right].                       \tag{19}
\]

Let

\[
 k_0=k(c^{2L}),\qquad \widetilde k=\sqrt{1-k_0^2},
 \qquad \widetilde K=K(\widetilde k^2).
\]

The nome of `tilde k` is `exp(-pi^2/(2L ell))`; comparison with the
Fourier series of `dn` therefore gives

\[
 S(x)={\widetilde K\over L\ell}
 \operatorname{dn}(\widetilde Kx/L\mid\widetilde k^2),
 \qquad
 w_m=\operatorname{dn}(\widetilde Km/L\mid\widetilde k^2). \tag{20}
\]

The quarter-period identity

\[
 \operatorname{dn}(\widetilde K-u\mid\widetilde k^2)
 ={k_0\over\operatorname{dn}(u\mid\widetilde k^2)}
\]

now proves the exact reflection law

\[
 w_mw_{L-m}=k_0,\qquad w_0=1,\qquad w_L=k_0.           \tag{21}
\]

## 6. The optimal metric

Define

\[
 P=D^{-1}WD^{-1}
   =\operatorname{diag}(p_m),\qquad p_m={w_m\over c^m}. \tag{22}
\]

By (5) and (18),

\[
 P-T_c^*PT_c=D^{-1}(W-FWF)D^{-1}\succeq0.             \tag{23}
\]

Moreover, for every real `x`,

\[
 \operatorname{sech}(x+\ell)\ge
 e^{-\ell}\operatorname{sech}(x)=c\operatorname{sech}(x),
\]

because `cosh(x+ell)<=e^ell cosh(x)`.  Applying this term by term in
(16) gives `S_{m+1}>cS_m`, hence

\[
 1=p_0<p_1<\cdots<p_L={k_0\over c^L}.                 \tag{24}
\]

Thus

\[
 I\le P\le {k(c^{2L})\over c^L}I,
\]

which proves the upper bound matching (9).  It also explains the
rank-one defect and reflection products observed in every L116 SDP.

Finally, Jacobi's product is

\[
 k(Q)=4\sqrt Q\prod_{r\ge1}
  \left({1+Q^{2r}\over1+Q^{2r-1}}\right)^4.            \tag{25}
\]

Every factor is strictly below one.  With `Q=c^(2L)`,

\[
 {k(c^{2L})\over c^L}<4,\qquad
 {k(c^{2L})\over c^L}=4-16c^{2L}+O(c^{4L}),            \tag{26}
\]

completing (1).

## 7. Independent regeneration

The SDP/eigensolver checker is

```bash
.venv/bin/python -u experiments/crabb_elliptic_axis.py \
  --output experiments/crabb_elliptic_axis_s70222.jsonl
```

For `p=3,...,10` and `c in {0.05,0.15,0.4,0.6}`, the explicit metric
satisfies the rank-one Stein identity to at worst `6e-11` in binary64;
the SDPs independently recover the same optimum.

The separate 260-decimal checker uses explicit DCT-I/Jacobi matrices
and no SDP:

```bash
.venv/bin/python -u experiments/crabb_elliptic_axis_theorem.py \
  --output experiments/crabb_elliptic_axis_theorem_s70222.jsonl
```

It checks (11)--(24) for `p=3,...,30` and
`c in {0.01,0.05,0.15,0.4,0.6,0.8}`.  In particular it compares the
kernel diagonal, the periodized-sech series, and the complementary
Jacobi formula through independently computed expressions.

## 8. Scope and next target

This theorem closes the fixed-weight all-size gap in Li's proposed
elliptic Crabb similarity and completes L115's sole soft normal axis.
It does **not** prove a neighbourhood of `C_p` for `p>=4`: the margin
in (26) tends to zero rapidly as `c->0` or `p` grows, while L65's
coercive transverse complement must still be controlled uniformly
against domain and metric remainders.

The next target is therefore a tubular certificate around this exact
axis, using the circular-range manifold as the anchor and L65 for the
strong transverse directions.  It is not a generic `5 x 5` slice
calculation and not yet the general Crouzeix conjecture.
