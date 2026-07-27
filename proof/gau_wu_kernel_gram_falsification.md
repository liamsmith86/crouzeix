# A universal model-kernel lower bound is false

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
The attractive proposed estimate was

\[
 \boxed{S_\phi\succeq3DG_g^{-1}D.}                \tag{3}
\]

The constant \(3\) is the terminal weight in L342's exact endpoint
square, and (3) survived the initial generic samples.  It would have
given a short model-kernel proof of the remaining sign.

## 2. A separated rational witness

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
precision at 256, 512, 1024, and in an untracked audit 4096 Fourier
nodes.  The tracked three-resolution drift is below \(2\cdot10^{-11}\).

## 3. Consequence

Do not try to sign the Gau--Wu shape form by assigning the complete
terminal weight \(3\) to the isolated model-kernel Gram after the
first pivot.  L342's two endpoint residual squares must remain
coupled to the model-dependent port and disk-fibre elimination.
A weaker or corrected Gram comparison could still be true, but it
must first pass (4).

This falsification does **not** refute:

1. the numerical positivity of \(M_\phi\);
2. L345's phase-covariance criterion;
3. L346--L350's lower-flag reductions; or
4. the scalar or completely bounded Crouzeix conjectures.

## 4. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_kernel_gram_falsification.py
```

The checker reconstructs the candidate independently at 256, 512,
and 1024 nodes, verifies (6), the positive sign of the full shape
matrix, the negative eigenvalue in (5), and resolution stability.
The dataset SHA-256 is
`2eaf011e0abb51ab6b496e1b9d7b49b44e7b877349450c254cbb4a9342287c26`.
