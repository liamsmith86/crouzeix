# Root-tracking normal form of the Gau--Wu diagonal flag

> **Status and scope.**  The decomposition below is exact at every
> finite nondegenerate Gau--Wu model.  Conditional on L347's cyclic
> coefficient vanishing, it turns L346's diagonal equations into
> ordinary eigenvalue/Blaschke-zero tracking plus one endpoint
> splitting equation.  It does not prove that the L345 dual lift
> satisfies those equations.

## 1. Spectral root velocities

Retain the upper-triangular Gau--Wu matrix \(A\), with diagonal

\[
 0,b_1,\ldots,b_{n-2},0,
 \tag{1}
\]

and let \(P_j\) be the spectral projection at the simple eigenvalue
\(b_j\).  Put

\[
 P_0=I-\sum_{j=1}^{n-2}P_j.                       \tag{2}
\]

Thus \(P_0\) is the two-dimensional algebraic spectral projection
for the repeated endpoint eigenvalue zero.

For a first normalized operator direction \(G\), define

\[
 \mu_j(G)=\operatorname {tr}(P_jG),\qquad
 \mu_0(G)={1\over2}\operatorname {tr}(P_0G).       \tag{3}
\]

The first number is the ordinary first velocity of the simple
eigenvalue \(b_j\).  The second is the mean first velocity of the
two-member zero cluster; it does not assume that the defective pair
has two analytic eigenvalue branches.

Let \(v_j\) be the velocity of the matching zero of the moving
Blaschke product \(f_\varepsilon\), with \(v_0\) denoting its final
zero at the origin.  Write

\[
 Y_1=Df(A)[G]+\dot f(A).                           \tag{4}
\]

## 2. Exact diagonal formulas modulo the cyclic corner

Let \(Z=qp^*=e_Le_0^*\) be the southwest matrix unit.  Define the
base scalars

\[
\begin{aligned}
 \kappa_j&=
 { (Df(A)[Z])_{jj}
   -f'(b_j)\operatorname {tr}(P_jZ)
  \over f'(0)},\qquad 1\le j\le n-2,\\
 \kappa_0&=
 { (Df(A)[Z])_{00}+(Df(A)[Z])_{LL}
   -f'(0)\operatorname {tr}(P_0Z)
  \over f'(0)} .
\end{aligned}                                     \tag{5}
\]

Nondegeneracy gives \(f'(0)\ne0\).  Then every simultaneous
operator/zero direction satisfies

\[
\boxed{
 (Y_1)_{jj}
 =f'(b_j)\{\mu_j(G)-v_j\}
  +\kappa_j(Y_1)_{L0},
 \quad 1\le j\le n-2,}                            \tag{6}
\]

and

\[
\boxed{
 (Y_1)_{00}+(Y_1)_{LL}
 =2f'(0)\{\mu_0(G)-v_0\}
  +\kappa_0(Y_1)_{L0}.}                           \tag{7}
\]

Consequently, on the cyclic-free face from L347,

\[
\begin{aligned}
 (Y_1)_{jj}=0
 &\quad\Longleftrightarrow\quad v_j=\mu_j(G),
 &&1\le j\le n-2,\\
 (Y_1)_{00}+(Y_1)_{LL}=0
 &\quad\Longleftrightarrow\quad v_0=\mu_0(G).
\end{aligned}                                     \tag{8}
\]

The only diagonal condition not contained in (8) is the endpoint
split

\[
 (Y_1)_{00}-(Y_1)_{LL}=0.                         \tag{9}
\]

Thus L346's \(n+1\) scalar generators have the following exact
organization:

1. one top weighted support coefficient, by L347;
2. \(n-2\) simple-root tracking equations;
3. one mean repeated-root tracking equation; and
4. one endpoint-splitting equation.

## 3. Proof by the triangular commutator gauge

The same southwest recurrence used in L346 gives a first-order
similarity gauge \(X\), with \(X_{L0}=0\), for which

\[
 G=U+cZ+[A,X],                                    \tag{10}
\]

where \(U\) is upper triangular.  The southwest coefficient is
unchanged by the commutator because the last row and first column
of \(A\) vanish.  L347 therefore gives

\[
 (Y_1)_{L0}=f'(0)c.                               \tag{11}
\]

Functional calculus differentiates similarity equivariantly:

\[
 Df(A)[[A,X]]=[f(A),X].                           \tag{12}
\]

Since \(f(A)\) is supported at \((0,L)\), the commutator in (12)
has zero interior diagonal, and its two endpoint diagonal entries
sum to zero.  For the upper-triangular term,

\[
\begin{aligned}
 (Df(A)[U])_{jj}&=f'(b_j)U_{jj},\\
 (Df(A)[U])_{00}+(Df(A)[U])_{LL}
 &=f'(0)(U_{00}+U_{LL}).
\end{aligned}                                     \tag{13}
\]

Spectral traces kill commutators, so

\[
\begin{aligned}
 U_{jj}
 &=\mu_j(G)-c\operatorname {tr}(P_jZ),\\
 U_{00}+U_{LL}
 &=\operatorname {tr}(P_0G)
   -c\operatorname {tr}(P_0Z).
\end{aligned}                                     \tag{14}
\]

Finally, differentiating the simple zero factor gives

\[
 \dot f(b_j)=-f'(b_j)v_j,\qquad
 \dot f(0)=-f'(0)v_0.                             \tag{15}
\]

Substitution of (11)--(15), with the \(cZ\) response collected
through (5), proves (6)--(7).

## 4. Revised Euler target

The numerical L345 dual lift already satisfies, to solver
precision,

\[
 v_j=\mu_j(G),\qquad
 v_0=\mu_0(G),\qquad
 (Y_1)_{00}=(Y_1)_{LL}=0,
 \tag{16}
\]

once its cyclic residual is resolved.  Equation (16) is evidence,
not part of this lemma.

The proof target is now ordered:

1. prove L347's highest negative support coefficient vanishes;
2. derive the root-tracking Euler equations (8) from the
   zero-velocity part of L344/L345;
3. prove the remaining endpoint split (9); and
4. invoke L346 for the rest of the lower flag.

Do not use the generally false ungauged formula
\((Y_1)_{jj}=f'(b_j)(G_{jj}-v_j)\).  A nonzero cyclic corner closes
upper-triangular paths and contributes the explicit correction in
(6).

## 5. Deterministic audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_root_tracking_normal_form.py
```

The checker constructs the Hermite spectral projections
polynomially, evaluates the response to \(Z\), and compares
(6)--(7) on every real joint basis direction.  In twelve models
through dimension eight, the largest relative interior and endpoint
residuals are \(1.87\cdot10^{-15}\) and
\(1.71\cdot10^{-15}\).  The least separated dimension-eight model
has spectral-projector idempotence error
\(1.28\cdot10^{-10}\); this conditioning diagnostic does not enter
the response identities.

The dataset SHA-256 is
`d87d5cc883aefe7ad77ab0a7c4ec81c8e8643d9779dfa2c4160b4c2e3875a04c`.
