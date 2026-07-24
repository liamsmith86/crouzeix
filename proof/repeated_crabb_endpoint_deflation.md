# Exact endpoint-response covariance under delay deflation

## 1. Result (L216, 2026-07-24)

Let \(S\) be a spectrally stable finite pure partial isometry with
equal orthogonal defect frames

\[
 I-S^*S=VV^*,\qquad I-SS^*=WW^*,\qquad V^*W=0,
\]

and transfer coefficients

\[
 B_j=W^*(S^*)^jV.
\]

Fix \(r\ge1\) and assume

\[
 B_1=\cdots=B_r=0.                                  \tag{1}
\]

L209 permits removal of the complete \(r\)-stage wandering line

\[
 {\cal L}_r=\bigoplus_{j=0}^{r-1}S^jW\mathbb C^m.
\]

Let \(J:{\cal H}_r={\cal L}_r^\perp\to{\cal H}\) be inclusion and put

\[
 S_r=J^*SJ,\qquad V_r=J^*V,\qquad W_r=J^*S^rW.       \tag{2}
\]

Then \(S_r\) is a partial isometry with defect frames \(V_r,W_r\).
For the stable Stein inverse

\[
 {\cal G}_S(H)=\sum_{n\ge0}(S^*)^nHS^n,
\]

every forcing supported on the retained space satisfies

\[
\boxed{
\begin{aligned}
 J^*{\cal G}_S(JH_rJ^*)J&={\cal G}_{S_r}(H_r),\\
 W^*{\cal G}_S(JH_rJ^*)W
 &=W_r^*{\cal G}_{S_r}(H_r)W_r .
\end{aligned}}                                      \tag{3}
\]

Thus the *linear physical endpoint response* is exactly covariant
under every number of complete delay removals.

More precisely, use the two physical equality metrics and operators

\[
\begin{aligned}
P&=2I-VV^*+2WW^*,&T&=P^{-1/2}SP^{1/2},\\
P_r&=2I-V_rV_r^*+2W_rW_r^*,&
T_r&=P_r^{-1/2}S_rP_r^{1/2}.                       \tag{4}
\end{aligned}
\]

If \(\widehat C=J\widehat C_r\) and
\(V_r^*\widehat C_r=0\), then L204's endpoint maps obey

\[
\boxed{
 {\cal M}_{T}(P^{1/2}\widehat C)
 ={\cal M}_{T_r}(P_r^{1/2}\widehat C_r).
}                                                   \tag{5}
\]

On (1), L212's grade-\(r+1\) column is also the literal lift of
L207's grade-one column for the deflated colligation:

\[
\boxed{
\begin{aligned}
\widehat C_{r+1}
&=-\frac72(I-VV^*)S^{r+1}WB_{r+1},\\
J^*\widehat C_{r+1}
&=-\frac72(I-V_rV_r^*)S_rW_r\widetilde B_1,
\qquad \widetilde B_1=B_{r+1}.
\end{aligned}}                                      \tag{6}
\]

Consequently (5) transports L207's channel coboundary exactly:

\[
\boxed{
{\cal M}_{T}(P^{1/2}\widehat C_{r+1})
=28\left\{
 {\mathfrak C}_S(B_{r+1}^*B_{r+1})
 -B_{r+1}B_{r+1}^*
\right\}.
}                                                   \tag{7}
\]

L216 closes the linear endpoint-response part of the
L209--L211 physical covariance problem and gives a geometric proof of
L212 on every fully delayed stratum.  It does **not** yet transport
the nonlinear raw Riemann/metric base: lower metric coefficients and
both endpoint Schur squares interact before order \(2r+2\).  That
prepared-base covariance is the remaining gate.

## 2. Retained Stein compression

L209 proves from (1) that

\[
W,SW,\ldots,S^rW
\]

are pairwise orthonormal isometries, that \({\cal H}_r\) is invariant
for \(S\), and that (2) gives the complete deflated defects.  Hence

\[
S^nJ=JS_r^n\qquad(n\ge0).                           \tag{8}
\]

For \(H=JH_rJ^*\), insert (8) termwise in the convergent Stein series:

\[
\begin{aligned}
J^*{\cal G}_S(H)J
&=\sum_{n\ge0}(S_r^*)^nH_rS_r^n
={\cal G}_{S_r}(H_r).
\end{aligned}
\]

This proves the first identity in (3).

For the upper endpoint, the first \(r\) orbit columns lie in the
removed line, while

\[
S^{r+n}W=JS_r^nW_r\qquad(n\ge0).                   \tag{9}
\]

Therefore

\[
\begin{aligned}
W^*{\cal G}_S(H)W
&=\sum_{n\ge0}(S^nW)^*H(S^nW)\\
&=\sum_{n\ge0}(S_r^nW_r)^*H_r(S_r^nW_r)\\
&=W_r^*{\cal G}_{S_r}(H_r)W_r,
\end{aligned}
\]

which proves the second identity in (3).  This is an exact orbit
shift, not an asymptotic coefficient argument.

## 3. Return to physical coordinates

Let

\[
H_r=V_r\widehat C_r^*+\widehat C_rV_r^*,\qquad
H=JH_rJ^*.
\]

Since \(V=JV_r\) and \(\widehat C=J\widehat C_r\),

\[
H=V\widehat C^*+\widehat CV^*.
\]

Conjugating the balanced Stein solutions by the corresponding metric
roots in (4) gives the physical Stein solutions for the columns

\[
C=P^{1/2}\widehat C,\qquad C_r=P_r^{1/2}\widehat C_r.
\]

Both left defects are eigenvectors of their physical metrics with
eigenvalue four:

\[
P^{1/2}W=2W,\qquad P_r^{1/2}W_r=2W_r.
\]

Thus both physical upper compressions are four times the balanced
compressions in (3).  This proves (5).  The hypothesis
\(V_r^*\widehat C_r=0\) is the usual lower-preserving normalization;
the covariance identity itself only needs retained support.

## 4. Pullback of the all-grade column

Put \(k=r+1\).  On (1), every contamination term in L212 vanishes, so
its column is the first line of (6).  L210 gives

\[
\widetilde B_1=B_k.
\]

Also

\[
S_rW_r=J^*S^kW,\qquad
(I-V_rV_r^*)J^*=J^*(I-VV^*).                       \tag{10}
\]

Equations (10) prove the second line of (6).  Apply (5), then L207 to
the deflated grade-one system.  L211's full-delay channel identity is

\[
{\mathfrak C}_{S_r}(B_k^*B_k)
={\mathfrak C}_{S}(B_k^*B_k),
\]

in the common copy coordinates.  This proves (7).

The argument explains why the same polynomial channel correction
works at every grade.  It deliberately does not claim that the
complete raw coefficient at order \(2k\) is a delayed copy of the
grade-one raw coefficient; L214--L215 prove that assertion only after
grade-specific lower-order gauges and future-row preparations.

## 5. Regeneration

Run

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -u \
  experiments/repeated_crabb_endpoint_deflation.py \
  --output \
  experiments/repeated_crabb_endpoint_deflation_s70224.jsonl
```

The checker covers one through five removed delay stages, defect
multiplicities one through three, inflated unstructured partial
isometries, independently gauged heterogeneous shifts, and repeated
apices.  It verifies both identities in (3), the physical identity
(5), the column pullback (6), and the endpoint formula (7).  These
tests audit the construction; equations (8)--(10) are the all-size
proof.  The tracked dataset SHA-256 is
`bc1e168e33f4d9b7bb668df60688f07cffe7fe1a6ed5e58e62101bd231c1b5d5`.
