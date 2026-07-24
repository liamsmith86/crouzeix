# The grade-two circular-normal leading cancellation (L161, 2026-07-23)

## 1. Result and scope

Put `p=L+1` and `L>=4`.  Take the real phase-one grade-two
phase-palindromic disk coordinate with amplitude `a`, and write `c`
for the elliptic parameter.  Its compact reflected coordinate is

\[
 w_2=ac^2.                                           \tag{1}
\]

Among L160's coercive circular normals, circle character permits only
the mode-`L` representative

\[
 N_L={E_{L-1,0}+E_{L,1}\over\sqrt2}                 \tag{2}
\]

to couple to `w_2` on the first possible face.  The apparent leading
coupling vanishes:

\[
\boxed{
 [\,w_2\overline y\,]\Gamma_p=0,
 \qquad yN_L\in N_{\rm coercive}.
}                                                    \tag{3}
\]

Equivalently, give `a` and `c` weight one and `y` weight three.  The
weight-six rank-one-envelope face is

\[
\boxed{
 [\epsilon^6]\Gamma_p
 =-64a^2c^4+e_p(yN_L),
}                                                    \tag{4}
\]

with no linear term in `y`.  Thus the grade-two face is strictly
easier than the exceptional offset-one face L160: its first possible
normal Schur loss is absent.

The `e_p` term in (4) is L65's already optimized matrix Hessian.
Holding the flat defect (6) fixed is valid by the envelope theorem
only for the linear `y` derivative; its raw `y^2` coefficient is not
the optimized Hessian and is deliberately not used.

This is a leading weighted selection theorem.  It does **not** yet
assert the uniform Rees-division estimate needed to absorb every
higher term near the axes `a=0` and `c=0`.  Those axes have separate
margins from L117/L121 and the exact disk theorem, but their merger
still has to be written in one tubular chart.  Higher grades also
remain open.

## 2. Fourier selection

Under L65's circle action, a matrix entry `Y_jm` has character
`m-j-1`.  The real circular-normal quotient has one complex support
normal in each mode `3,...,p`; mode two is the elliptic soft normal.

Insert the grade-two equality pair and two elliptic reflections.
At total reflected weight three, their character pairs only with
normal mode `p-1=L`.  The real Riesz representative of that support
functional is exactly (2): the endpoint weights in L115 equation (4)
are both `1/sqrt(2)`.  Every other coercive normal therefore has zero
weight-six cross by character, before any coefficient calculation.
Conjugation kills the mixed real/imaginary phase, and circle rotation
turns the real calculation below into the complex invariant pairing.

It remains to show that the scalar coefficient on (2) also vanishes.

## 3. Exact weighted path

Let `H_0` be the half-diagonal Toeplitz coordinate Gramian and let
`H_2` have offsets `2` and `L-2` (one copy when they coincide).  Put

\[
\begin{aligned}
H_\epsilon&=H_0+a\epsilon H_2,\\
K_\epsilon&=H_\epsilon+R^*H_\epsilon R,\\
X_\epsilon&=
 2K_\epsilon^{-1/2}H_\epsilon R K_\epsilon^{-1/2},\\
A_{\epsilon,s}
 &=X_\epsilon+c\epsilon X_\epsilon^*
   +s\epsilon^3N_L.                                  \tag{5}
\end{aligned}
\]

The first two terms in (5) lie on the exact physical
disk-equality/ellipse family.  Let `T_{\epsilon,s}` be the inverse
Riemann pullback and solve L118's analytic rank-one defect through
three jets:

\[
q_*(\epsilon)=e_0+\epsilon q_1+\epsilon^2q_2
                  +\epsilon^3q_3+O(\epsilon^4).       \tag{6}
\]

The defect Hessian at the Crabb point is invertible, so the three
successive stationarity systems are uniquely solvable.  The envelope
theorem permits (6) to be held fixed when differentiating in `s`.

Write `lambda_-` and `lambda_+` for the simple endpoint eigenvalues of
the resulting Stein Gramian; their base values are one and four.
Normal stationarity and Fourier selection remove every lower-weight
condition-ratio term.  Therefore

\[
\partial_s[\epsilon^6]{\lambda_+\over\lambda_-}
 =
\partial_s[\epsilon^6](\lambda_+-4\lambda_-).        \tag{7}
\]

## 4. Endpoint cancellation

Use L160 equation (10) to construct `K_epsilon^(-1/2)` without a
symbolic matrix root, then use L62's finite support recurrence and
the intermediate-normalized simple-eigenpair recurrence.

For `L>=7`, every weight-six path is confined to separated endpoint
neighbourhoods after the mode-`L` Fourier filter.  Substitution gives

\[
\partial_s[\epsilon^6]\lambda_-
=\partial_s[\epsilon^6]\lambda_+=0.                  \tag{8}
\]

Increasing `L` only inserts unused middle sites, so (8) is the
all-size stable case.

The short chains `L=4,5` also give two zero endpoint derivatives.
At `L=6` the endpoint paths collide and the individual derivatives
are nonzero:

\[
\begin{aligned}
d_-&=-2c(a^2+2ac-2c^2),\\
d_+&=-8c(a^2+2ac-2c^2)=4d_- .                       \tag{9}
\end{aligned}
\]

Thus (7) is zero even in the sole collision case.  Equations
(8)--(9) prove (3) for every `L>=4`.

On `s=0`, the same recurrence gives

\[
[\epsilon^6]\Gamma_p=-64a^2c^4,                     \tag{10}
\]

including the nondivisor case `L=7`.  This agrees with the prepared
compact face from L149--L150 and proves (4).

## 5. Exact regeneration

Run

```bash
PYTHONPATH=experiments OPENBLAS_NUM_THREADS=1 \
  .venv/bin/python -u \
  experiments/crabb_grade_two_normal_selection.py \
  --output experiments/crabb_grade_two_normal_selection_s70223.jsonl
```

The checker:

1. constructs the exact physical path (5);
2. regenerates the inverse square root and Riemann pullback;
3. solves all three defect jets successively;
4. lifts both simple Stein endpoint eigenpairs through weight six;
5. keeps `a` and `c` independent for `p=5,6,7`;
6. verifies the collision polynomial (9); and
7. checks the first two separated sizes `p=8,9`.

The finite sizes audit the short collisions and the first stable
instances.  The all-size step is the support separation leading to
(8), not extrapolation from the two separated records.

## 6. Next gate

The observed pattern now suggests:

\[
\text{grade }k\ge2:
\quad\hbox{the first character-allowed circular-normal cross vanishes}.
\]

The next efficient test is grade three at weight eight, using the
same reusable series engine.  In parallel, the uniform tubular proof
must convert the extra weighted order in (3) into a bound valid
through the disk and elliptic axes; no further finite-size slice grind
should precede that division audit.
