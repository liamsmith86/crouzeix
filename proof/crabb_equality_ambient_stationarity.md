# Ambient stationarity of the Crabb equality ridge (L162, 2026-07-23)

## 1. The theorem

Use the phase-one gauge of L123.  Thus `p=L+1`, the positive
Toeplitz block has

\[
 z_k=\overline {z_{L-k}},
\quad
 K=H+R^*HR,
\quad
 A=2K^{-1}HR,
\quad q=He_0,
\]

and

\[
 M=K-qq^*+2rr^*,\qquad M-A^*MA=qq^* .
\]

For a nearby coefficient-coordinate operator `A+sE`, keep `K,q`
fixed and let `M_s` solve

\[
 M_s-(A+sE)^*M_s(A+sE)=qq^* .                       \tag{1}
\]

Let `U(A+sE)` be the condition number of `(M_s,K)`.  If `h_E` is the
first Schwarz/Riemann correction determined by the support motion of
`W(A+sE)`, then

\[
\boxed{DU(A)[E-h_E(A)]=0\quad\hbox{for every }E\in M_p(\mathbb C).}
                                                               \tag{2}
\]

Consequently every positive phase-palindromic equality anchor is an
ambient critical point of one explicit feasible rank-one certificate
after numerical-range normalization.  This is the all-size version
of L72's `p=3` stationarity theorem.

This is a first-order critical-ridge theorem, not a neighborhood
inequality.  It does not supply the uniform negative normal Hessian
or finish the Crabb local merger.

## 2. Touching upper and lower certificates

L123 gives generalized eigenvalues

\[
 \operatorname{spec}(M,K)=\{1/2,1,\ldots,1,2\}
\]

with normalized endpoint vectors `sqrt(2)e_0,sqrt(2)e_L`.  It also
gives the finite Blaschke product

\[
 B(\xi)=g(\xi)/g^\sharp(\xi),\qquad
 B(A)=4e_0r^*,\qquad \|B(A)\|_K=2.                   \tag{3}
\]

For small real `s`, (1) makes `A+sE` a contraction in the `M_s`
metric.  Hence

\[
 \|B(A+sE)\|_K^2\leq U(A+sE).
\]

Both sides equal four at zero and both active endpoints are simple.
Their derivatives therefore agree.  Singular-value differentiation
in (3) gives

\[
 DU(A)[E]=8\operatorname{Re}J(E),\qquad
 J(E):=q^*DB(A)[E]e_L.                               \tag{4}
\]

This touching argument avoids differentiating the optimized
similarity problem and also shows why fixing `q` in (1) is sufficient.

## 3. The endpoint boundary formula

Put

\[
\begin{aligned}
 G(\xi)&=g^\sharp(\xi)
       =1+2\sum_{k=1}^{L-1}z_k\xi^k,\\
 f(\xi)&=(1,\xi,\ldots,\xi^L)^T .
\end{aligned}
\]

The first key identity is

\[
\boxed{
 J(E)=\int_{\mathbb T}
 {\bar w\,f(w)^*KEf(w)\over |G(w)|^2}\,dm(w).
}                                                     \tag{5}
\]

Here is an all-size residue proof.  L156's endpoint resolvents are

\[
 (wI-A)^{-1}e_L={c(w)\over wg(w)},\qquad
 q^*(wI-A)^{-1}={d(w)K\over wg(w)},
\]

where

\[
 c=f+Ge_0,\qquad
 d=(w^L,w^{L-1},\ldots,1)+Ge_L^* .
\]

Insert these in the Cauchy formula for `DB(A)[E]`.  On the circle,
phase palindromy gives

\[
 \overline {G(w)}=w^{-L}g(w).
\]

The difference between the resulting integrand and (5) has matrix
numerator

\[
 {cd-f(w^L,w^{L-1},\ldots,1)\over G}
 =e_0(w^L,\ldots,1)+fe_L^*+Ge_0e_L^*.                \tag{6}
\]

After division by `w^2g`, every entry of (6) is `O(w^{-2})` at
infinity.  All zeros of `g` are inside the circle by L123, so the
residue theorem makes the difference integral zero.  This proves
(5) without a generic-root assumption.

## 4. The Schwarz transfer identity

Write

\[
\begin{aligned}
 \widetilde G(t)&=1+2\sum_{k=1}^{L-1}\bar z_kt^k
                 =\det(I-tA),\\
 \widetilde N(t)&=L+2\sum_{k=1}^{L-1}(L-k)\bar z_kt^k .
\end{aligned}
\]

Direct substitution of L123's companion columns in (3) gives the
endpoint transfer function

\[
\boxed{
 \sum_{j\geq0}J(A^{j+1})t^j
 =q^*B'(A)A(I-tA)^{-1}e_L
 ={\widetilde N(t)\over\widetilde G(t)} .
}                                                     \tag{7}
\]

Equivalently, the numerator is
`L*widetilde_G-t*widetilde_G'`.  Multiplication by
`widetilde_G` reduces (7) to the companion recurrence, so (7) is a
finite polynomial identity.

For `w` on the circle put

\[
 N(w)=L+2\sum_{k=1}^{L-1}(L-k)z_kw^k,\qquad
 D(w)=f(w)^*Kf(w).
\]

Entrywise summation of the Toeplitz diagonals gives

\[
 D(w)=2\operatorname{Re}N(w)-L.                      \tag{8}
\]

Phase palindromy and a coefficient comparison give the sharper
kernel identity

\[
\boxed{
 2\operatorname{Re}{N(w)\over G(w)}-L
 ={D(w)\over|G(w)|^2}.
}                                                     \tag{9}
\]

One can also obtain (9) by inserting
`N=LG-wG'`; the off-zero Fourier coefficients telescope in reflected
pairs.

## 5. Cancellation after the Riemann pullback

Let

\[
 s_E(w)
 ={\operatorname{Re}\{\bar w f(w)^*KEf(w)\}\over D(w)}
\]

be the first support motion.  With Fourier convention
`\widehat s(k)=int s(w)w^{-k}dm(w)`, the first inverse-Riemann
correction is

\[
 h_E(\xi)=\widehat s(0)\xi
          +2\sum_{k\geq1}\widehat s(k)\xi^{k+1}.       \tag{10}
\]

The roots of `widetilde G` lie outside the closed disk, so (7) may be
inserted termwise in (10).  Taking real parts and using (9) gives

\[
\begin{aligned}
 \operatorname{Re}J(h_E(A))
 &=\int_{\mathbb T}s_E(w)
   \left(2\operatorname{Re}{N(w)\over G(w)}-L\right)dm\\
 &=\int_{\mathbb T}
   {\operatorname{Re}\{\bar w f(w)^*KEf(w)\}
    \over |G(w)|^2}\,dm\\
 &=\operatorname{Re}J(E),                            \tag{11}
\end{aligned}
\]

where the last equality is (5).  Equations (4) and (11) prove (2).

The phase-removal gauge in L123 transports the theorem to every
unimodular phase branch.

## 6. Regeneration and audit scope

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/crabb_equality_ambient_stationarity.py \
  --output experiments/crabb_equality_ambient_stationarity_s70223.jsonl
```

The exact part uses nonreal Gaussian-rational phase-palindromic
anchors and checks through `p=6`:

1. the Stein identity;
2. the finite rational identity (7);
3. the boundary kernel identity (9); and
4. every matrix entry of the residue cancellation (6).

The independent numerical part checks every real and imaginary
matrix unit through `p=8`.  It constructs (10) by a 4096-point FFT
with 512 analytic coefficients and differentiates both generalized
metric endpoints.  The unnormalized ambient gradient is large,
whereas the pulled gradient is at roundoff.

The exact proof is all-size; the finite records audit its algebra and
do not replace the residue/companion arguments.

## 7. Consequence for the campaign

The equality set in L123 now has the first ingredient of a tubular
Morse--Bott argument: an explicit feasible upper certificate with no
linear ambient term in normalized coordinates.  L160 and L161 give
coercive information on initial circular-normal faces.  The next
gate is to combine these facts in a prepared chart and prove a
uniform negative normal estimate through the stratified Crabb apex.
