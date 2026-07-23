# A complete-`2` neighbourhood of the single `3 x 3` Crabb block (2026-07-22)

## 1. Theorem and scope

Let

\[
 C_3=\begin{bmatrix}0&\sqrt2&0\\0&0&\sqrt2\\0&0&0\end{bmatrix}.
\]

There is a neighbourhood `N` of `C_3` in the space of complex `3 x 3`
matrices such that, for every `A in N`, its numerical range is a complete
`2`-spectral set for `A`.  Equivalently, if `phi_A` maps the interior of
`W(A)` conformally to the unit disk, then

\[
 t_*(\phi_A(A))\leq4.                                      \tag{1}
\]

This is a local theorem around one irreducible Crabb block.  It says nothing
about repeated blocks, larger Crabb blocks, or arbitrary matrices far from this
neighbourhood, and therefore does not resolve the general conjecture.

## 2. Analytic certificate

The top eigenvalue of

\[
 \operatorname{Re}(e^{-i\theta}A)
\]

is simple with a uniform gap at `C_3`, and its support curve is the unit circle.
Both properties persist nearby.  Finite-dimensional analytic perturbation theory
therefore makes the support function, and hence the strictly convex boundary
embedding

\[
 \zeta_A(\theta)=e^{i\theta}(h_A(\theta)+ih_A'(\theta)),    \tag{2}
\]

real analytic jointly in the matrix entries and in `theta` in a fixed analytic
strip.

The normalized Riemann map also depends real analytically on `A` in this
neighbourhood.  One direct local proof writes its boundary equation as

\[
 \Psi_A(e^{it})=\zeta_A(t+\eta_A(t)).                      \tag{3}
\]

Use exponentially weighted Fourier spaces for the analytic map `Psi` and the
real boundary reparameterization `eta`, and fix the three disk-automorphism
normalizations.  At a circle, the linearization of (3) is the standard Fourier
splitting into nonnegative analytic modes and a real tangential displacement.
More explicitly, after division by `e^{it}`, modes `k<=-2` determine the negative
Fourier modes of the real function `eta`; reality determines its positive modes,
and the remaining three real modes are fixed by the disk-automorphism
normalizations.  The analytic part is then forced.  This gives a bounded inverse
in the exponentially weighted Fourier norm.  The real-analytic Banach
implicit-function theorem gives the asserted dependence.  This is also a very
special smooth-near-circle instance of the analytic-dependence theorems of Rodin
(1986) and Wu (1993).

Consequently

\[
 T(A)=\phi_A(A)                                             \tag{4}
\]

is real analytic.  Its spectral radius is strictly below one nearby.

For `q=(1,x,y)^T`, let `P(T,q)` be the unique Stein solution

\[
 P-T^*PT=qq^*.                                             \tag{5}
\]

It is positive definite and real analytic near
`(T,q)=(C_3,e_1)`.  Since the endpoint eigenvalues of
`P_0=diag(1,2,4)` are simple,

\[
 k(T,q)=\frac{\lambda_{\max}P(T,q)}{\lambda_{\min}P(T,q)}  \tag{6}
\]

is real analytic there.  The exact defect calculation gives

\[
 k(C_3,(1,\epsilon x,\epsilon y))-4
 =\epsilon^2\left(8|x|^2+\frac83|y|^2\right)+O(\epsilon^3). \tag{7}
\]

Thus the defect Hessian is positive definite.  The analytic implicit-function
theorem supplies a unique local critical defect `q_*(A)`, which is a strict local
minimum of (6).  Define the feasible-certificate excess

\[
 g(A)=k(T(A),q_*(A))-4.                                    \tag{8}
\]

Every `g(A)<=0` proves (1); global optimality of the defect is not required.
The order-by-order defect optimization in L67--L70 is precisely the Taylor series
of this unique analytic branch.

This local critical value is independent of the chosen Riemann-map gauge.  If
`b` is a disk automorphism, the standard defect identity

\[
 P-b(T)^*Pb(T)=(1-|\alpha|^2)(I-\alpha T^*)^{-1}
 (P-T^*PT)(I-\bar\alpha T)^{-1}                            \tag{9}
\]

(up to the harmless unimodular convention in `b`) preserves defect rank and
the condition of `P`; the inverse automorphism gives the converse.  Unitary
similarity transports `(P,q)` and also preserves (6).

## 3. Slice and exact critical manifold

By L69, affine-unitary invariance reduces every nearby matrix to the exact slice

\[
 A(z,w,s,v)=C_3+
 \begin{bmatrix}
 -z/3&s&\bar z\\w&2z/3&-s\\v&w&-z/3
 \end{bmatrix},                                           \tag{10}
\]

where `z,w,v` are complex and `s` is real.  The circle action is

\[
 (z,w,s,v)\mapsto(e^{i\theta}z,e^{2i\theta}w,s,e^{3i\theta}v). \tag{11}
\]

The L71 disk locus extends equivariantly from the real calculation as

\[
 w_d(z)=\frac{3\sqrt2}{64}z^2R(|z|^2),\qquad
 v_d(z)=-\frac9{64}z^3R(|z|^2),\qquad s_d=0,               \tag{12}
\]

where

\[
 81t^2R^2+(1152t-4096)R+4096=0,\qquad R(0)=1.             \tag{13}
\]

L72 proves more than constancy on this curve: after conformal pullback the
rank-one condition derivative vanishes in every complex ambient direction.
The defect envelope in (8) has the same derivative.  Hence

\[
 g(z,w_d(z),0,v_d(z))=0,\qquad Dg(z,w_d(z),0,v_d(z))=0.    \tag{14}
\]

Put

\[
 U=w-w_d(z),\qquad V=v-v_d(z).                             \tag{15}
\]

L63 gives the Hessian in the three strong real variables `(s,Re V,Im V)` at
the origin:

\[
 -8s^2-\frac{21}{4}|V|^2.                                 \tag{16}
\]

It is negative definite and remains so nearby.  The implicit-function theorem
therefore gives a unique analytic critical pair
`(s_*(z,U),V_*(z,U))`, which is the local maximum in the strong variables.
Set

\[
 H(z,U)=g(z,w_d(z)+U,s_*(z,U),v_d(z)+V_*(z,U)).            \tag{17}
\]

After shrinking the neighbourhood, strict concavity gives

\[
 g(z,w_d+U,s,v_d+V)\le H(z,U).                             \tag{18}
\]

Equation (14) and uniqueness give

\[
 H(z,0)=0,\qquad D_UH(z,0)=0.                              \tag{19}
\]

## 4. The one-soft-coordinate sign

The certificate and both analytic optimizations preserve (11), conjugation, and
reversal-transpose symmetry.  Thus

\[
 H(e^{i\theta}z,e^{2i\theta}U)=H(z,U).                    \tag{20}
\]

Every Taylor monomial `z^a zbar^b U^c Ubar^d` obeys

\[
 a-b+2(c-d)=0.                                             \tag{21}
\]

By (19), every nonzero monomial has `c+d>=2`.  At degrees two through five,
condition (21) allows only

\[
 |U|^2;\quad \varnothing;\quad |z|^2|U|^2,\ |U|^4;\quad
 \bar z^2U^2\bar U\text{ and its conjugate}.              \tag{22}
\]

The degree-two coefficient is zero by L63.  The exact L70 completed square,
after maximizing its strong variables, gives the mixed degree-four coefficient
`-25/56`.  At `z=0`, circle weights two and three forbid every term linear in
`V` that depends only on `U`, while reversal-transpose invariance makes the
certificate even in `s`.  Thus maximizing the strong variables cannot alter the
pure mode-two quartic; L68 gives its coefficient `-4`.  Therefore

\[
 H(z,U)=-\frac{25}{56}|z|^2|U|^2-4|U|^4+\mathcal R(z,U).  \tag{23}
\]

Analyticity, (19), and the monomial restrictions give, on a sufficiently small
polydisk,

\[
 |\mathcal R(z,U)|
 \le C\bigl(|z|^4|U|^2+|z|^2|U|^3
       +|z|^2|U|^4+|U|^5\bigr).                           \tag{24}
\]

Each term on the right is a small multiple of
`|z|^2|U|^2+|U|^4`.  Shrinking once more, (23)--(24) imply

\[
 H(z,U)\le
 -\frac{25}{112}|z|^2|U|^2-2|U|^4\le0.                   \tag{25}
\]

Combining (18) and (25) proves `g<=0` throughout the slice.  L69 and invariance
then prove (1) for every matrix in a full neighbourhood of `C_3`.

Finally, (5) and `k<=4` mean that `T(A)` is similar to a contraction through a
similarity of condition at most two.  Von Neumann's matrix-valued inequality and
the conformal pullback give the stated complete `2`-spectral-set theorem.

## 5. Audit and regeneration

The finite checks used above are regenerated by

```bash
.venv/bin/python -u experiments/p3_disk_morse_bott.py
.venv/bin/python -u experiments/p3_disk_center_tangent.py
.venv/bin/python -u experiments/p3_crabb_weighted_slice.py
.venv/bin/python -u experiments/p3_local_theorem_probe.py
```

The first script proves (7), enumerates every symmetry-allowed low-order soft
monomial, and verifies the `-25/56` strong-variable elimination.  The other two
regenerate L72 and L70.  All decisions in these scripts are exact.
The last command is a non-load-bearing randomized map/Gramian smoke test; its
default 18 mixed, weighted, and superweighted samples reached at most
`3.999999999999`, with maximum map diagnostic below `9.8e-13`.

The only imported analytic fact is local analytic dependence of the normalized
Riemann map under analytic Jordan-boundary perturbations; Section 2 also gives
the near-circle implicit-function reduction needed here.  No finite numerical
sampling or unproved optimizer-globality assertion enters the sign argument.
