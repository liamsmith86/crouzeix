# Every polarized endpoint response has a gap-free energy bound

## 1. Result (L296, 2026-07-26)

Retain L280's bistochastic transfer channel

\[
\Phi(K)=\sum_{n\ge1}B_nKB_n^*,\qquad
A_Y=\Phi^*(Y),
\]

and put

\[
K_{\rm M}=I-\Phi\Phi^*.
\]

For any copy matrix \(X\), define the polarized grade-\(k\) flux

\[
\boxed{
{\cal F}_k(X)
=\operatorname {sym}(B_kX^*)
 -\Phi\!\left(\operatorname {sym}(X^*B_k)\right).} \tag{1}
\]

This is one eighth of L280's exact physical endpoint response.  Every
Hermitian test \(Y\) satisfies

\[
\boxed{
\langle Y,{\cal F}_k(X)\rangle
=\operatorname {Re}\operatorname {tr}
X^*(YB_k-B_kA_Y),}                                \tag{2}
\]

and hence

\[
\boxed{
|\langle Y,{\cal F}_k(X)\rangle|
\le
\|X\|_F\sqrt{\langle Y,K_{\rm M}Y\rangle}.}        \tag{3}
\]

More generally, for any square-summable family \(X_k\),

\[
\boxed{
\left|\left\langle
Y,\sum_k{\cal F}_k(X_k)
\right\rangle\right|
\le
\left(\sum_k\|X_k\|_F^2\right)^{1/2}
\sqrt{\langle Y,K_{\rm M}Y\rangle}.}              \tag{4}
\]

No spectral gap, inverse, compactness argument, rank projection, or
pseudoinverse occurs.

The same statements hold on every L282 flag: simply take
\(\widehat Y=UYU^*\) in (2)--(4).  In particular, if the compressed
mixed physical remainder has the form

\[
U^*R_{\rm mix}U
=U^*\left\{\sum_j\beta_j{\cal F}_j(X_j)\right\}U, \tag{5}
\]

and its polarization energy obeys

\[
\boxed{
\sum_j|\beta_j|^2\|X_j\|_F^2
\le C^2\|U^*B_k\|_F^2,}                           \tag{6}
\]

then L282's required estimate follows immediately:

\[
\boxed{
\langle Y,U^*R_{\rm mix}U\rangle_+
\le C\|U^*B_k\|_F
\sqrt{\langle\widehat Y,K_{\rm M}\widehat Y\rangle}.}          \tag{7}
\]

Thus after L295's branch separation, the repair-branch gate is not
to force every remainder into the special quadratic flux \(G_k\).
It is enough, **modulo favorable Grams and prior-flag factors whose
compression is zero**, to produce L280 polarization columns and
prove the finite energy bound (6).

## 2. Exact pairing

For Hermitian \(Y\),

\[
\begin{aligned}
\left\langle Y,\operatorname {sym}(B_kX^*)\right\rangle
&=\operatorname {Re}\operatorname {tr}(X^*YB_k),\\
\left\langle
Y,\Phi(\operatorname {sym}(X^*B_k))
\right\rangle
&=\operatorname {Re}\operatorname {tr}(X^*B_kA_Y).
\end{aligned}
\]

The second identity uses channel adjointness and \(A_Y=A_Y^*\).
Subtracting proves (2).

L280's exact Dirichlet identity is

\[
\langle Y,K_{\rm M}Y\rangle
=\sum_n\|YB_n-B_nA_Y\|_F^2.                       \tag{8}
\]

Cauchy--Schwarz applied to the \(k\)-th summand proves (3).
Summing (2) first and applying Cauchy--Schwarz across all transfer
grades proves (4).

For a flag test \(\widehat Y\), the same proof is literal.  If a
polarization column is left-supported on \(P=UU^*\), only its
supported norm enters.  For example,

\[
P{\cal F}_k(PB_k)P
=P\{B_kB_k^*-\Phi(B_k^*PB_k)\}P
=-G_{k,P},
\]

so L294's flagged estimate is the special choice \(X=PB_k\).

Finally absorb scalar coefficients into
\(\widetilde X_j=\beta_jX_j\), apply (4) with
\(\widehat Y\), and use (5)--(6).  This proves (7).

## 3. Scope and next gate

L296 does not prove (5)--(6) for the physical L285/L289 remainder.  It
identifies the exact quantity that must be controlled:

> the Frobenius energy of the polynomial polarization columns, not a
> Markov inverse and not the norm of separately enormous endpoint
> pieces.

L288's factored sextic transport illustrates a different favorable
class: it is a bounded prior-channel endpoint factor, so its flag
compression is exactly zero and it need not be a Markov response
(indeed a general such factor need not have zero trace).  The
arbitrary-grade problem is to use L289/L290 to discard those
prior-flag factors first, then derive polarization columns only for
the surviving compressed remainder and prove (6) through terminal
grade \(L\).  A failure should be translated to L292's Smith
valuations; it must not be hidden by a pseudoinverse.

## 4. Exact regeneration

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/repeated_crabb_polarized_flux_inequality.py \
  --output \
  experiments/repeated_crabb_polarized_flux_inequality_s70226.jsonl
```

The exact SymPy checker uses noncommuting depolarizing and
random-unitary-type bistochastic channels.  It verifies the
single-grade, summed, and flagged pairing identities, L280's
Dirichlet identity, and every Cauchy--Schwarz slack without floating
tolerances.  The tracked dataset has SHA-256

```text
58381726281a8102eb3a53b334b3d0ac11b7dcfaccb3eae3cc64d9750d344c48
```
