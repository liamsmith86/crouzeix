# The delayed elliptic axis gauge

## 1. Result (L213, 2026-07-24)

Retain L202--L212's balanced equality data

\[
 I-S^*S=VV^*=:E,\qquad I-SS^*=WW^*=:F,\qquad
 V^*W=0,
\]

and the physical metric and operator

\[
 P=2I-E+2F,\qquad T=P^{-1/2}SP^{1/2}.
\]

Let

\[
 B_1=W^*S^*V.
\]

On the delayed face \(B_1=0\), the complete second elliptic Stein
equation has the explicit balanced metric coefficient

\[
\boxed{\widehat X_{\rm ax}=SFS^*-S^*ES.}            \tag{1}
\]

If \(\widehat F_2\) is L207's fixed balanced second forcing and
\(K=V^*\widehat F_2V\), the corresponding full second defect-frame
column is

\[
\boxed{
 \widehat C_{2,\rm ax}
 =-\frac12VK+2(S^*)^4V.}                            \tag{2}
\]

It obeys

\[
\widehat X_{\rm ax}-S^*\widehat X_{\rm ax}S
=\widehat F_2+
V\widehat C_{2,\rm ax}^*
+\widehat C_{2,\rm ax}V^*.                         \tag{3}
\]

Both physical endpoint coefficients vanish:

\[
\boxed{
 V^*P^{1/2}\widehat X_{\rm ax}P^{1/2}V=0,\qquad
 W^*P^{1/2}\widehat X_{\rm ax}P^{1/2}W=0.}          \tag{4}
\]

Formula (2) is polynomial in the equality colligation.  It uses no
kernel projection or pseudoinverse and is analytic through the
repeated Crabb apex.

On a repeated length-\(L\) monomial shift, (1) is exactly the
balanced \(c^2\) coefficient of the all-size elliptic-axis metric in
`proof/crabb_elliptic_axis.md`.  Thus (1)--(2), rather than the
minimum-norm L207 representative, are the axis-compatible
normalization before exposing \(B_2\).

L213 does **not** prove the fourth-order \(B_2\) endpoint.  It removes
a false gauge obstruction: higher coefficients must be computed
after (1)--(2), because endpoint-null second-frame motions change the
fourth-order base.

## 2. Ordered second-jet identity

L207 writes

\[
\begin{aligned}
 J&=(I+F)S^*(I+E),\\
 \dot S&=J-S^3,\\
 \ddot S&=2S-(S^2J+SJS+JS^2)+S^5,\\
 A&=-S^*\dot SV,\\
 \widehat F_2
 &=AA^*+\ddot S^*S+S^*\ddot S+\dot S^*\dot S .
\end{aligned}                                      \tag{5}
\]

The delayed hypothesis is equivalently

\[
 FS^*E=WB_1V^*=0,\qquad ESF=0.                     \tag{6}
\]

Use only (6) and

\[
\begin{gathered}
S^*S=I-E,\quad SS^*=I-F,\quad EF=FE=0,\\
SE=0,\quad ES^*=0,\quad S^*F=0,\quad FS=0.
\end{gathered}                                      \tag{7}
\]

Direct order-preserving expansion of (5) first gives

\[
\begin{aligned}
J&=S^*+FS^*+S^*E,\\
\dot S&=S^*+FS^*+S^*E-S^3,\\
\ddot S&=-S-S^2FS^*-S^*ES^2+S^5,
\end{aligned}                                      \tag{8}
\]

and hence

\[
\begin{aligned}
\widehat F_2-E\widehat F_2E
={}&-F+SFS^*-S^*ES-2ES^4\\
&+(S^*)^2ES^2-2(S^*)^4E.                          \tag{9}
\end{aligned}
\]

On the other hand,

\[
\widehat X_{\rm ax}-S^*\widehat X_{\rm ax}S
=-F+SFS^*-S^*ES+(S^*)^2ES^2.                     \tag{10}
\]

Subtracting (9) from (10) proves

\[
\boxed{
\begin{aligned}
&\widehat X_{\rm ax}-S^*\widehat X_{\rm ax}S\\
&\quad-\{\widehat F_2-E\widehat F_2E\}
=2\{ES^4+(S^*)^4E\}.
\end{aligned}}                                    \tag{11}
\]

No copy matrices are commuted in (8)--(11).  The parallel column
\(-VK/2\) contributes \(-E\widehat F_2E\).  The perpendicular
column \(2(S^*)^4V\) contributes

\[
2\{ES^4+(S^*)^4E\},
\]

so (11) proves (3).

## 3. Endpoint and axis checks

Since \(SV=0\), \(S^*W=0\), and \(B_1=0\),

\[
\begin{aligned}
V^*\widehat X_{\rm ax}V
 &=(W^*S^*V)^*(W^*S^*V)=0,\\
W^*\widehat X_{\rm ax}W
 &=-(W^*S^*V)(W^*S^*V)^*=0.
\end{aligned}
\]

The endpoint eigenvalues of \(P\) are one on \(V\) and four on
\(W\), so physical conjugation preserves their vanishing.  This
proves (4).

For the repeated monomial shift, \(E\) is level zero and \(F\) is
level \(L\).  Hence (1) is

\[
 -\Pi_1+\Pi_{L-1},
\]

with cancellation when \(L=2\).  The periodized-sech formula for the
exact axis metric gives, for \(0<j<L\),

\[
 p_j(c)=2-2c^{2j}+2c^{2(L-j)}+o(c^2).
\]

After conjugating the physical coefficient by \(P^{-1/2}\), its
\(c^2\) term is exactly \(-\Pi_1+\Pi_{L-1}\).  This proves the stated
axis compatibility in every length and copy multiplicity.

## 4. A universal endpoint-null direction

There is also a useful null column valid without \(B_1=0\):

\[
\boxed{\widehat C_{\rm null}=2W.}                  \tag{9}
\]

Its balanced forcing is \(2(VW^*+WV^*)\).  Because
\(SV=0\) and \(S^*W=0\), the next Stein iterate vanishes, so its
balanced metric response is the forcing itself.  In physical
coordinates,

\[
X_{\rm null}=4(VW^*+WV^*),
\]

whose two diagonal endpoint compressions are zero.  Therefore (9)
can be added to L207 without changing either endpoint.  At the
length-four monomial apex, \(2(S^*)^4V=2W\), explaining the
axis-alignment observed in the first grade-two probe; outside that
special case, (2) is the correct delayed gauge.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_endpoint_null_gauge.py \
  --output \
  experiments/repeated_crabb_endpoint_null_gauge_s70224.jsonl
```

The 37 deterministic records include general stable partial
isometries, unstructured stable partial isometries constrained only
by \(B_1=0\), and repeated monomial shifts of lengths two through
eight with multiplicities one through three.  They audit (2)--(4),
the universal null direction (9), and exact axis normalization.  The
identities above, not the floating audit, prove L213.  The tracked
dataset SHA-256 is
`918dc868b9874a206a7f8bd1af617af9a42deb005d09d6aaa48f5a6e30b6cb87`.
