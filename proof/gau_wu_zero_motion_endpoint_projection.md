# Zero motion is an exact two-ended model-space distance

> **Scope.**  This note closes the sector of L336 in which the
> compressed shift is fixed and only the extremal inner function
> moves.  It does not control the mixed operator/zero terms or prove
> L336's full symmetric-flux sign.

## 1. Exact identity (L337/A284, 2026-07-26)

Let \(f\) be inner with \(f(0)=0\), put

\[
 \phi=zf,\qquad {\cal K}={\cal K}_\phi,\qquad
 S=S(\phi),\qquad p=f,\quad q=1.
\]

For every inner function \(g\), define \(Y=g(S)\).  Then

\[
\boxed{
 1-\|Yq\|^2
 =1-\|Y^*p\|^2
 =\|P_{\phi H^2}g\|_{H^2}^2\geq0.}                \tag{1}
\]

In particular, if

\[
 g_\varepsilon=f+\varepsilon h+O(\varepsilon^2)
\]

is any differentiable inner curve, the two endpoint coefficients in
L336 satisfy

\[
\boxed{
 r_2=\ell_2=\|P_{\phi H^2}h\|_{H^2}^2,\qquad
 r_2+\ell_2=2\|P_{\phi H^2}h\|_{H^2}^2.}          \tag{2}
\]

Thus pure zero motion is an exact Hardy projection square.  Any
remaining difficulty in L336 is necessarily physical or mixed.

## 2. Proof

Model functional calculus at the cyclic vector \(q=1\) gives

\[
 g(S)q=P_{\cal K}g.                               \tag{3}
\]

Since \(g\) is inner, \(\|g\|_{H^2}=1\).  The orthogonal
decomposition \(H^2={\cal K}\oplus\phi H^2\) therefore yields

\[
 1-\|Yq\|^2=\|P_{\phi H^2}g\|^2.                 \tag{4}
\]

It remains to identify the other endpoint.  Let \(J_\phi\) be the
canonical conjugation on \({\cal K}_\phi\):

\[
 (J_\phi u)(\zeta)
 =\phi(\zeta)\overline{\zeta u(\zeta)}
 =f(\zeta)\overline{u(\zeta)}
 \quad\text{a.e. on }\mathbb T.                  \tag{5}
\]

Write \(g=u+\phi v\), where \(u=P_{\cal K}g\).  On the circle,

\[
 f\overline g
 =f\overline u+f\overline\phi\,\overline v
 =J_\phi u+\overline\zeta\,\overline v.           \tag{6}
\]

The last summand lies in the strictly negative Hardy space and is
orthogonal to \({\cal K}\).  The compressed functional calculus gives

\[
 g(S)^*p=P_{\cal K}(f\overline g)=J_\phi u.       \tag{7}
\]

Because \(J_\phi\) is antiunitary,

\[
 \|Y^*p\|=\|u\|=\|Yq\|.
\]

Together with (4), this proves (1).

At \(g=f\), one has \(f\in{\cal K}_\phi\), so
\(P_{\phi H^2}f=0\).  Expanding the last member of (1) along the
inner curve gives

\[
 \|P_{\phi H^2}g_\varepsilon\|^2
 =\varepsilon^2\|P_{\phi H^2}h\|^2+O(\varepsilon^3),
\]

which proves (2).  No zero acceleration enters the coefficient.

## 3. Consequence for the live gate

The pure zero block now has a canonical Hilbert metric rather than an
opaque finite Schur matrix.  For a simultaneous physical direction
\(C\) and inner tangent \(h\), the next calculation should express
the mixed term as a Riesz pairing against
\(P_{\phi H^2}h\), complete this projection square, and leave one
purely physical remainder.  That remainder, not (2), is the live
content of L336.

The identity also explains why the endpoint pairing is natural:
\(J_\phi p=q\) and \(J_\phi S J_\phi=S^*\).  A one-ended
calculation discards this exact model conjugation.

## 4. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/gau_wu_zero_motion_endpoint_projection.py \
  --output \
  experiments/gau_wu_zero_motion_endpoint_projection_s70224.jsonl
```

The checker fixes random Gau--Wu model shifts in dimensions
\(3,\ldots,8\), moves every zero of \(g\), and independently compares
both matrix endpoint defects with the Fourier projection norm in
(1).  The 12-record dataset has SHA-256
`233e72245d54a1d71ad1442332d0d26d9dbc5bff207d43c2e42f35390cad8b4f`;
the largest residual is \(1.43\cdot10^{-15}\).
