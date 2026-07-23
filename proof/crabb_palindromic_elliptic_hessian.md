# Exact finite elliptic Hessians on the disk equality stratum (2026-07-23)

## 1. Status

This note upgrades A84's floating-point evidence to exact rational
formal algebra in every tested size.  It is **not** the all-size
recurrence proof and is not a uniform remainder theorem.

Let `p=L+1`, choose a phase-one palindromic coefficient vector supported
at offsets `k` and `L-k`, and let `Gamma(a,c)` be the locally optimized
rank-one Stein condition square minus four along A84's exact curve.
The exact checker proves, for every tested pair,

\[
 [a^2]\Gamma(a,c)=-64c^{2k}+O(c^{2k+1}).              \tag{1}
\]

All coefficients below `c^(2k)` vanish as rational numbers; the leading
coefficient is exactly `-64`, not a fit.

The deterministic audit covers every first offset in sizes `p=3,...,8`:

\[
\begin{array}{c|c}
p&k\\ \hline
3&1\\
4&1\\
5&1,2\\
6&1,2\\
7&1,2,3\\
8&1,2,3.
\end{array}                                           \tag{2}
\]

The calculation also recovers A84's observed endpoint pattern.  For
`k=1`, the coefficient of `c^2` in the quadratic metric term is `48`
at the lower generalized endpoint and `128` at the upper endpoint.
For the tested `k=2,3`, the coefficients at `c^(2k)` are `-16` and
`-128`.  Since the axis endpoint eigenvalues in the chosen scale are
`2` and `8`, both cases give

\[
 [c^{2k}]\left({\lambda_+\over\lambda_-}\right)
 =M^{(2)}_{LL}-4M^{(2)}_{00}=-64.                    \tag{3}
\]

Equation (3) is a discovered recurrence target, not yet an induction
in arbitrary `k`.

## 2. Coordinate model

No matrix square root is expanded.  In L123's coefficient coordinates,

\[
\begin{aligned}
A(a)&=C+aE,\\
S(a,c)&=A(a)+cJA(a)J,\\
T(a,c)&=\phi_c(S(a,c)).
\end{aligned}                                         \tag{4}
\]

Here `C` has first superdiagonal weight two and all remaining weights
one, while

\[
 E=2(e_0-e_L)
 \sum_{j=2}^L u_{j-1}e_j^*                            \tag{5}
\]

is rank one.  The phase-palindromic identity makes `K(a)` centrosymmetric,
so its `K(a)`-adjoint is exactly `JA(a)J`.  Thus (4) is exact.

L125's inverse-map ODE is solved over rational truncated power series,
then reverted coefficient by coefficient to obtain `phi_c`.  Matrix
powers are retained only through amplitude degree two:

\[
 T=T_0+aT_1+a^2T_2+O(a^3).                            \tag{6}
\]

The checker verifies the direct-map Newton edge

\[
 [w^{2n+1}]\phi_c(w)=(-1)^nc^n+O(c^{n+2})             \tag{7}
\]

as part of the same exact calculation.

## 3. Exact axis anchor

L117 supplies the physical diagonal axis metric

\[
 p_m={S_m\over S_0c^m},\qquad
 S_m=\sum_{r\in\mathbb Z}
 \operatorname{sech}((m+2Lr)(-\log c)).               \tag{8}
\]

Each hyperbolic term is expanded rationally through

\[
 \operatorname{sech}(j(-\log c))
 ={2c^{|j|}\over1+c^{2|j|}}.                          \tag{9}
\]

The implementation forms `S_m/c^m` directly.  This matters: dividing
an already truncated `S_m` would silently erase the endpoint guard
orders.  After multiplying by the diagonal `K(0)` and scaling by two,
the coordinate metric `M_0` has generalized endpoints `2` and
`k(c^(2L))/c^L * 2`.  The checker independently factors

\[
 M_0-T_0^*M_0T_0=d_0d_0^*.                            \tag{10}
\]

The odd-series reversion is performed at the internal order
`output_order+2*maximum_scalar_degree+4`, as required by L125's
valuation-one divisions, and only then truncated.  Thus every
coefficient below the declared output order is audited.  An earlier
implementation used one common truncation order; its terminal scalar
coefficients were harmless through `p=7` but polluted the first
`p=8,k=3` extension.  The guarded implementation reproduces all old
records and gives the correct exact size-eight jet.

## 4. Second-order Stein and endpoint formulas

Fix the defect scale and write

\[
 d=d_0+ax.
\]

Solving

\[
 M-T^*MT=dd^*
\]

through amplitude degree two gives

\[
 M=M_0+aM_1+a^2M_2+O(a^3).                            \tag{11}
\]

The solution is obtained by exact fixed-point iteration in the
truncated ring.  It terminates because the constant operator is the
nilpotent Crabb shift; every additional long path gains either an
amplitude degree or a positive `c` valuation.

Let `K=K_0+aK_1`.  Write the two endpoint expansions as
`lambda_r+a ell_r+a²q_r`.  Direct generalized-eigenvalue perturbation
gives

\[
 q_r
 ={(M_2)_{rr}\over(K_0)_{rr}}
 -{\ell_r(K_1)_{rr}\over(K_0)_{rr}}
 -\sum_{j\ne r}
 {((M_1)_{rj}-\lambda_r(K_1)_{rj})^2
  \over
  (K_0)_{rr}((M_0)_{jj}-\lambda_r(K_0)_{jj})}.         \tag{12}
\]

The coefficient of `a²` in `lambda_L/lambda_0` also includes the two
standard products `-ell_L ell_0/lambda_0²` and
`lambda_L ell_0²/lambda_0³`.  An audit initially omitted these terms;
the checker now retains them and independently confirms that the
reported leading coefficients are unchanged.

Therefore the condition-square Hessian is an exact quadratic polynomial
in the free defect tangent `x`.  The script reconstructs that quadratic
by polarization and solves its stationarity system by Gaussian
elimination over rational power series.  Its leading Hessian is L118's
positive defect Hessian, so every pivot is a unit.

## 5. The exposed optimizer recurrence

The exact optimizer begins with L123's disk-equality defect.  Lower
elliptic terms then propagate by two coordinates.  Representative
prefixes are:

\[
\begin{array}{c|c|c}
(p,k)&\text{nonzero disk terms}&\text{first correction}\\ \hline
(4,1)&x_1=x_2=2&x_3=-8c+\cdots\\
(5,2)&x_2=2-4c^2+\cdots&x_4=-8c+\cdots\\
(7,3)&x_3=2-4c^2+\cdots&x_5=-8c+\cdots .
\end{array}                                           \tag{13}
\]

Those corrections cancel every putative condition term below
`c^(2k)`.  The first residual endpoint contribution is (3).

This identified the remaining all-size task sharply.  L133
subsequently completed it for the real phase-one offset `k=1`:
the explicit defect jet
`2e_1+2e_(L-1)-8ce_3` gives the universal endpoint pair `(48,128)`
and hence `-64c²` in every size.  For general grades the remaining
items are:

1. extend the stationarity/endpoint recurrence beyond the offset-one
   case or replace it by unequal-residue localization;
2. polarize complex phases and different coefficient grades to obtain
   the diagonal Hardy
   sum, not only the one-coordinate cases.

The finite checker does not justify replacing these remaining steps
by pattern extrapolation.

## 6. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_palindromic_elliptic_hessian.py \
  --output experiments/crabb_palindromic_elliptic_hessian_s70223.jsonl

PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_palindromic_elliptic_hessian.py \
  --minimum-size 8 --maximum-size 8 \
  --output experiments/crabb_palindromic_elliptic_hessian_p8_s70223.jsonl
```

The implementation uses only `fractions.Fraction`; no floating-point
optimizer, eigensolver, elliptic-function evaluation, or saved expansion
coefficient enters the audit.
