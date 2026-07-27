# Fixed model-kernel multiples fail for the Gau--Wu shape form

> **Status and scope.**  This is a reproducible numerical
> falsification of one stronger-than-needed candidate inequality for
> the Gau--Wu boundary-shape form.  It does not challenge the observed
> positivity of that form, L342--L350, or the Crouzeix conjecture.

## 1. The rejected candidate

Let \(m=n-1\), let \(M_\phi\) denote the negative Hermitian
conformal-shape matrix suggested by A292, and eliminate its first
coordinate:

\[
 S_\phi=(M_\phi)_{22}
 -(M_\phi)_{21}(M_\phi)_{11}^{-1}(M_\phi)_{12}.
 \tag{1}
\]

If \(g\) is the product of the \(n-2\) nonzero interior Blaschke
factors, let \(G_g\) be the coefficient Gram of

\[
 {1-g(z)\overline{g(w)}\over1-z\overline w}
 \tag{2}
\]

in the monomials \(1,\ldots,z^{n-3}\), and put
\(D=\operatorname {diag}(1,\ldots,n-2)\).
The first attractive proposed estimate was

\[
 \boxed{S_\phi\succeq3DG_g^{-1}D.}                \tag{3}
\]

The constant \(3\) is the terminal weight in L342's exact endpoint
square, and (3) survived the initial generic samples.  It would have
given a short model-kernel proof of the remaining sign.

## 2. A moderate separated witness rejects weight three

Take \(n=4\) and

\[
 b_1={3\over8}+{i\over4},\qquad
 b_2=-{1\over2}-{i\over4}.
 \tag{4}
\]

The roots are distinct, nonzero, well inside the disk, and separated
by more than one.  Reconstructing the complete joint scalar Hessian,
eliminating all zero velocities, passing to the exact conformal
shape coordinate, and taking its Hermitian block gives

\[
 \lambda\!\left(S_\phi-3DG_g^{-1}D\right)
 \approx(-0.11585777,\ 1.17891488).               \tag{5}
\]

Thus (3) is indefinite, not positive semidefinite.  This is not a
soft-sign artefact: the full \(M_\phi\) remains positive definite,
with smallest eigenvalue approximately \(0.06403464\).  Its first
pivot also reproduces the independently observed exact value

\[
 (M_\phi)_{00}=16|b_1b_2|^2=1.015625             \tag{6}
\]

to \(2\cdot10^{-14}\).

The complete matrices and both eigenvalues are unchanged to displayed
precision from 512 through 2048 Fourier nodes.

## 3. A boundary witness rejects even weight one

Weakening three to two is not enough.  More decisively, take the
separated rational roots

\[
 b_1=-{41\over64}-{47i\over64},\qquad
 b_2=-{34\over64}-{52i\over64}.                  \tag{7}
\]

Their moduli are approximately \(0.9745,0.9708\), while their
separation is approximately \(0.1344\).  At this witness,

\[
 \lambda\!\left(S_\phi-DG_g^{-1}D\right)
 \approx(-275.39428827,\ 246481.1030077).         \tag{8}
\]

The full \(M_\phi\) is nevertheless positive definite, with smallest
eigenvalue approximately \(7.65574\).  The first-pivot check again
matches \(16|b_1b_2|^2\), and 512, 1024, and 2048 Fourier nodes give
relative matrix drift below \(4.1\cdot10^{-10}\).

Because \(DG_g^{-1}D\succ0\), (8) numerically rejects every fixed
coefficient at least one, not just the original coefficient three.
No claim is made here about smaller coefficients or a
model-dependent comparison.

## 4. Consequence

Do not try to sign the Gau--Wu shape form by assigning the complete
terminal weight, or even one complete copy, to an isolated
model-kernel Gram after the first pivot.  L342's two endpoint residual
squares must remain coupled to the model-dependent port and
disk-fibre elimination.  A genuinely model-dependent comparison
could still be true, but it must first pass both (4) and (7).

This falsification does **not** refute:

1. the numerical positivity of \(M_\phi\);
2. L345's phase-covariance criterion;
3. L346--L350's lower-flag reductions; or
4. the scalar or completely bounded Crouzeix conjectures.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_kernel_gram_falsification.py
```

The checker reconstructs both candidates independently at 512, 1024,
and 2048 nodes, verifies the first pivots, the positive sign of both
full shape matrices, the negative eigenvalues in (5) and (8), and
resolution stability.
The dataset SHA-256 is
`41352e20a5af43b1d4d1d3d6bfcb6ad193e055abfbb93c902928be08f2cd9d7a`.
