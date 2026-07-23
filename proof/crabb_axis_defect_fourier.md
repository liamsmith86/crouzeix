# Fourier Newton edge of the elliptic-axis defect (L137, 2026-07-23)

## 1. Statement

Use L117's coefficient coordinates for the length-`L` elliptic Crabb
axis, and scale its rank-one Stein defect

\[
M_c-T_c^*M_cT_c=d_c d_c^*
\]

so that `(d_c)_0=1+O(c^2)`.  Then every odd coordinate of `d_c`
vanishes, and for `1<=r<L/2`,

\[
\boxed{
(d_c)_{2r}=4(-1)^r c^r+O(c^{r+2}).
}                                                     \tag{1}
\]

Equivalently, before terminal aliasing the associated spatial
generating function is

\[
\boxed{
d_{\rm edge}(x)
=1+4\sum_{r\ge1}(-x)^r
=\frac{1-3x}{1+x},\qquad x=c{\cal S}^2.
}                                                     \tag{2}
\]

This is an all-size Fourier calculation, not an extrapolation from
the finite Hessian records.  It proves the axis-defect half of A98's
completed-square constants.  It does not prove the universal defect
Hessian factor or the Faber-row coupling.

## 2. DCT expression for the defect

Retain L117's notation

\[
q=c^2,\quad v_n=\frac{2Kn}{L},\quad
\delta_n=\frac{k'}{\operatorname{dn}(v_n\mid k^2)},\quad
\beta_n=\epsilon_n\delta_n .
\]

In the spectral node basis the elliptic Szegő kernel satisfies

\[
G-XGX=\beta\beta^*.
\]

Conjugating by the DCT-I matrix `U` gives defect vector `b=U beta`.
L135's coefficient eigenvector matrix is

\[
R=K_0^{-1/2}DU,\qquad D_{mm}=c^{m/2}.
\]

Transporting the spectral rank-one defect back to coefficient
coordinates therefore gives, up to one common scalar,

\[
(d_c)_m=\epsilon_m c^{-m/2}b_m.                     \tag{3}
\]

The common scalar cancels in `(d_c)_m/(d_c)_0`.

## 3. Reciprocal-dn Fourier series

The classical reciprocal Jacobi Fourier series (DLMF 22.11.6) is

\[
\frac{k'}{\operatorname{dn}(v\mid k^2)}
=\frac{\pi}{2K}
+\frac{2\pi}{K}\sum_{s\ge1}
 \frac{(-1)^s q^s}{1+q^{2s}}
 \cos\frac{s\pi v}{K}.                              \tag{4}
\]

At the Lobatto samples `v_n=2Kn/L`, (4) contains only the even
DCT-I modes `m=2s`.  Hence every odd component in (3) is zero.

Write `A_0=pi/(2K)` and

\[
A_s=\frac{2\pi}{K}\frac{(-1)^s q^s}{1+q^{2s}}.
\]

Trapezoidal DCT-I orthogonality gives, away from a terminal alias,

\[
b_0=\sqrt L\,A_0,\qquad
b_{2r}=\sqrt{\frac L2}\,A_r+\text{aliases}.           \tag{5}
\]

The endpoint factors in (3) cancel the square roots in (5), yielding

\[
\frac{(d_c)_{2r}}{(d_c)_0}
=c^{-r}\frac{A_r}{A_0}+\text{aliases}
=4(-1)^r\frac{c^r}{1+c^{4r}}+\text{aliases}.         \tag{6}
\]

The next local correction in (6) has at least two extra powers of
`c`.  A sampled alias has Fourier index congruent to `+/-r mod L`;
when `r<L/2`, its exponent after the factor `c^{-r}` is strictly
larger than `r`.  Finally `(d_c)_0=1+O(c^2)`, so (6) proves (1).
Summing its leading coefficients proves (2).

At `2r=L`, the two DCT aliases collide and the endpoint normalization
changes.  A98 does not infer that central case from (1); L130 already
handles it exactly.

## 4. Relevance to the completed square

The exact finite LDL audit in A98 finds the candidate Hessian outer
factor

\[
\ell_{\rm edge}(x)=\frac{1+x}{1-x}.
\]

Combining it with the proved (2) gives

\[
2\ell_{\rm edge}(x)d_{\rm edge}(x)
=2\frac{1-3x}{1-x}
=2-4\sum_{r\ge1}x^r.                                \tag{7}
\]

Thus the transported defect has a whitened first-reflection tail of
magnitude four.  L131's reflected Faber row has magnitude four.  If
the LDL factor and their opposite orientation are proved in the
homogeneous quotient, their gap is eight and the disk Hessian weight
four gives the observed improvement `4*8^2=256`.

## 5. Exact regression

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/crabb_axis_defect_fourier.py \
  --output experiments/crabb_axis_defect_fourier_s70223.jsonl
```

The checker reconstructs the axis metric and defect from the guarded
Stein engine and verifies the odd-mode cancellation and coefficient
(1) on a finite size grid.  The proof is the Fourier calculation
above, not that grid.
