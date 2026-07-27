# Postcomposition closes the Gau--Wu endpoint-trace Euler equation

> **Status and scope.**  The Hessian identity and its consequence for
> L345's dual lift are exact.  They prove that the sum of L346's two
> endpoint diagonal generators vanishes on the dual range.  L350
> subsequently proves that the entries agree, so both endpoints are
> now closed.  The cyclic coefficient, the interior root-tracking
> equations, phase covariance, and the Hessian sign remain open.

## 1. A canonical pure zero-motion tangent

Let

\[
 Y_0=f(A)=2pq^*,\qquad p=e_0,\quad q=e_L.
 \tag{1}
\]

For \(w\in\mathbb C\), postcompose the finite Blaschke product with
the disk automorphism

\[
 \beta_{tw}(z)={z-tw\over1-t\overline w z}.        \tag{2}
\]

This is a degree-preserving inner curve.  If \(a_j\) is any simple
zero of \(f\), including its terminal zero at the origin, its
velocity is

\[
 v_j={w\over f'(a_j)}.                            \tag{3}
\]

Denote the resulting pure zero-motion joint tangent by \(k_w\).
Since \(Y_0^2=0\),

\[
 {d\over dt}\bigg|_{t=0}\beta_{tw}(Y_0)=-wI.      \tag{4}
\]

Thus postcomposition supplies two canonical real zero tangents
\(k_1,k_i\) whose first image responses are \(-I,-iI\).

## 2. Exact Hessian/endpoint-trace identity

Let \(J_\phi\) be L345's real joint Hessian and
\({\cal B}_\phi\) its real symmetric polarization.  For a real joint
direction \(v\), let \(Y_1(v)\) be its first image and put

\[
 \tau(v)
 =\langle p,Y_1(v)p\rangle+
  \langle q,Y_1(v)q\rangle
 =(Y_1(v))_{00}+(Y_1(v))_{LL}.                    \tag{5}
\]

Then

\[
\boxed{
 {\cal B}_\phi(v,k_w)
 ={3\over4}\operatorname {Re}
  \{\overline w\,\tau(v)\}.}                      \tag{6}
\]

Equivalently, after complexifying the real row maps,

\[
\boxed{
 \tau={4\over3}\{k_1^TJ_\phi+i\,k_i^TJ_\phi\}.}    \tag{7}
\]

To prove (6), put \(Y=Y_0+\varepsilon Y_1+O(\varepsilon^2)\).
The mixed \(\varepsilon t\) coefficient in (2) is

\[
 \overline w(Y_0Y_1+Y_1Y_0).                     \tag{8}
\]

Its \((0,L)\) entry is \(2\overline w\,\tau(v)\).
The cross term from the first-row and final-column squares in the
simple singular-value formula is
\(-\tfrac12\operatorname {Re}\{\overline w\tau(v)\}\).
The coefficient of \(\varepsilon t\) is therefore
\(\tfrac32\operatorname {Re}\{\overline w\tau(v)\}\).
The polarized quadratic form is half that coefficient, which proves
(6).

No second zero acceleration enters: the first derivative in every
zero direction vanishes at the extremal.

## 3. Exact consequence for the dual lift

L345's dual range satisfies

\[
 J_\phi W_\phi=S_y^T+iS_x^T.                     \tag{9}
\]

The conformal shape map is zero on pure zero-motion directions.
Hence

\[
 k_1^TJ_\phi W_\phi=k_i^TJ_\phi W_\phi=0.         \tag{10}
\]

Equations (7) and (10) give the promised Euler identity:

\[
\boxed{
 (Y_1(W_\phi))_{00}+(Y_1(W_\phi))_{LL}=0.}        \tag{11}
\]

This conclusion uses neither numerical phase covariance nor the
unproved sign.

Combining (11) with L348 shows that, once L347's cyclic coefficient
vanishes, the terminal Blaschke zero automatically tracks the mean
velocity of the repeated zero cluster:

\[
 v_0={1\over2}\operatorname {tr}(P_0G).           \tag{12}
\]

Together with L350, the live lower-flag debt is reduced from L346's
\(n+1\) generators to:

1. L347's one cyclic/highest-mode equation;
2. the \(n-2\) simple-root tracking equations.

## 4. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_postcomposition_euler.py
```

The checker constructs the zero velocities (3), verifies (4)
against independent functional calculus, and compares the complete
joint endpoint-trace row with the right side of (7).  Twelve models
through dimension eight give maximum first-image and Euler-row
residuals \(2.67\cdot10^{-12}\) and
\(9.97\cdot10^{-12}\).  The largest values occur in the model where
the product \(f'(0)\) makes the zero-velocity coordinates most ill
conditioned.

The dataset SHA-256 is
`e572460e3f284788bcf5a121e93e49aaa0fb7f0439fa7e2d02c232fde483d1d4`.
