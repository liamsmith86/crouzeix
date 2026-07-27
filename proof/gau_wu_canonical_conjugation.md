# Canonical conjugation closes the Gau--Wu endpoint split

> **Status and scope.**  The covariance and endpoint conclusion
> below are exact at every finite nondegenerate Gau--Wu equality
> model, whenever L345's joint Hessian is nondegenerate so that its
> dual lift is defined.  Together with L349 they prove that both
> endpoint diagonal entries of the dual first image vanish.  The
> cyclic coefficient, the simple-root equations, phase covariance,
> and the Hessian sign remain open.

## 1. The equality model has an endpoint-exchanging conjugation

Let
\[
 \phi=zf,\qquad {\cal K}={\cal K}_\phi,\qquad
 S=S_\phi,
 \tag{1}
\]
and write the Gau--Wu decomposition as
\[
 {\cal K}=\mathbb Cp\oplus{\cal M}\oplus\mathbb Cq,
 \qquad p=f,\quad q=1.
 \tag{2}
\]
The canonical model-space conjugation
\[
 ({\cal C}g)(\zeta)
 =\phi(\zeta)\overline{\zeta g(\zeta)}
 \quad (|\zeta|=1)
 \tag{3}
\]
is antiunitary and satisfies
\[
 {\cal C}S^*{\cal C}=S,\qquad
 {\cal C}p=q,\qquad {\cal C}q=p.                 \tag{4}
\]
It preserves \({\cal M}\).

For
\[
 X=\sqrt2P_p+P_{\cal M}+{1\over\sqrt2}P_q,
 \qquad A=XSX^{-1},
 \tag{5}
\]
the endpoint exchange gives
\[
 {\cal C}X{\cal C}=X^{-1}.
 \tag{6}
\]
Consequently
\[
\boxed{{\cal C}A^*{\cal C}=A.}                   \tag{7}
\]
Thus every finite Gau--Wu equality model is complex symmetric under
an antiunitary which exchanges precisely the two singular
endpoints.

## 2. Reflection preserves both the sharp objective and its shape

Define the complex-linear involution on operators
\[
 \rho(T)={\cal C}T^*{\cal C}.                    \tag{8}
\]
It is an isometric anti-automorphism and
\[
 W(\rho(T))=W(T).                                \tag{9}
\]
For every analytic scalar function \(g\),
\[
 g(\rho(T))=\rho(g(T)),\qquad
 \|g(\rho(T))\|=\|g(T)\|.                        \tag{10}
\]

Apply (8) to an arbitrary joint curve \(T_t,g_t\), leaving \(g_t\)
unchanged.  Equation (9) says that the normalized numerical-range
domain and its Riemann map are unchanged, while (10) says that the
sharp scalar objective is unchanged.  At the fixed point (7), this
induces an orthogonal real involution
\[
 R_\phi:X_\phi\longrightarrow X_\phi             \tag{11}
\]
on L345's joint physical/zero tangent space.  It acts as
\(\dot T\mapsto\rho(\dot T)\) on the physical operator tangent and
as the identity on the moving function.

Let \(J_\phi\) be the joint scalar Hessian and \(S_\phi\) the
extended conformal shape map.  Exact invariance of the objective
and of the numerical range gives
\[
\boxed{
 R_\phi^TJ_\phi R_\phi=J_\phi,\qquad
 S_\phi R_\phi=S_\phi.}                          \tag{12}
\]
The second identity can also be read directly: \(T_t\) and
\(\rho(T_t)\) have the same numerical range for every \(t\), hence
the same first support variation and the same Fourier shape
coefficients.

## 3. The conformal dual lift is reflection-fixed

Write \(S_\phi=[S_x;S_y]\) and retain L345's dual columns
\[
 V_x=J_\phi^{-1}S_x^T,\qquad
 V_y=J_\phi^{-1}S_y^T,\qquad
 W_\phi=V_y+iV_x.                                \tag{13}
\]
Because an orthogonal involution is self-adjoint, (12) implies that
\(R_\phi\) commutes with \(J_\phi\) and fixes the range of
\(S_\phi^T\).  Therefore
\[
\boxed{
 R_\phi V_x=V_x,\qquad
 R_\phi V_y=V_y,\qquad
 R_\phi W_\phi=W_\phi.}                          \tag{14}
\]
No sign or phase-covariance assumption enters this conclusion.

Let \(Y_1(v)\) be the complexified first Blaschke-image response.
Functional-calculus equivariance under (8), including an arbitrary
moving analytic function, gives
\[
 Y_1(R_\phi v)=\rho(Y_1(v)).                     \tag{15}
\]
Combining (14) and (15),
\[
\boxed{
 Y_1(W_\phi)={\cal C}Y_1(W_\phi)^*{\cal C}.}      \tag{16}
\]

## 4. Endpoint equality and endpoint vanishing

Choose any orthonormal coordinate frame whose first and last vectors
are \(p,q\).  If \({\cal C}x=Q\overline x\), then \(Q\) is symmetric
unitary, \(Qe_0=e_L\), \(Qe_L=e_0\), and (16) reads
\[
 Y_1(W_\phi)=QY_1(W_\phi)^TQ^*.
 \tag{17}
\]
Taking the two endpoint diagonal entries gives the missing split:
\[
\boxed{
 (Y_1(W_\phi))_{00}=(Y_1(W_\phi))_{LL}.}          \tag{18}
\]

L349 independently proves their sum is zero.  Hence
\[
\boxed{
 (Y_1(W_\phi))_{00}
 =(Y_1(W_\phi))_{LL}=0.}                         \tag{19}
\]
This closes L348's endpoint-difference debt.  L348's repeated-zero
mean equation still additionally uses vanishing of L347's cyclic
corner, because that corner occurs in its spectral trace formula.

The remaining lower-flag obligations are now exactly:

1. L347's cyclic/highest weighted support coefficient; and
2. the \(n-2\) simple-root tracking equations.

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_canonical_conjugation.py
```

The checker reconstructs \(Q\) from \(AQ=QA^T\) and the two endpoint
conditions, verifies that it is symmetric unitary, builds the
induced real joint involution, and independently audits both
identities in (12), the fixed-point relation (14), the image
symmetry (17), and the endpoint split (18).  The tracked data cover
twelve models through dimension eight.  The largest conjugation,
shape-invariance, dual-fixed-point, image-symmetry, and endpoint
residuals are respectively
\(8.06\cdot10^{-15}\), \(4.38\cdot10^{-15}\),
\(3.41\cdot10^{-15}\), \(3.16\cdot10^{-12}\), and
\(4.49\cdot10^{-14}\).  The independently polarized Hessian, which
is the noisiest numerical object, has covariance residual at most
\(7.57\cdot10^{-10}\).

The dataset SHA-256 is
`8e317c98ea8216aae05b979cff0fda8b6766b25980e19b9c3e39e34ee96fed8b`.
