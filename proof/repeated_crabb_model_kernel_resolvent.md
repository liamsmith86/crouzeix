# The delayed state lift is the inverse model kernel

## 1. Result (L243, 2026-07-24)

Use L238--L239's retained pure partial isometry

\[
I-T^*T=VV^*=:E,\qquad I-TT^*=WW^*=:F,\qquad EF=0,
\]

and write

\[
\mathcal B(u)=W^*(I-uT^*)^{-1}V,\qquad
\mathcal B^\sharp(u)=V^*(I-uT)^{-1}W.
\]

On the Joukowski contour put

\[
x=\zeta^{-1},\qquad \rho=c/\zeta,\qquad
t=\rho x=c/\zeta^2.
\]

First absorb the right-defect correction into

\[
\begin{aligned}
G&=(\zeta I-T)(I-\rho T^*),\\
R_R&=(G-cT^*E)^{-1}.                              \tag{1}
\end{aligned}
\]

Then

\[
\boxed{
W^*R_RW
=\frac{x}{1-t}
\{I+\mathcal B(\rho)\mathcal B^\sharp(x)\}.}       \tag{2}
\]

For a removed delay of length \(r\), let

\[
\delta_r=\gamma_r-\rho,\qquad
\mathcal C(u)=u^r\mathcal B(u).
\]

The remaining endpoint denominator factors exactly as

\[
\boxed{
\begin{aligned}
I-\delta_rW^*R_RW
&=\frac{
 I-t^r\mathcal B(\rho)\mathcal B^\sharp(x)
}{1+t^r}\\
&=\frac{
 I-\mathcal C(\rho)\mathcal C^\sharp(x)
}{1+t^r}.                                        \tag{3}
\end{aligned}}
\]

Define the ordered two-variable model kernel

\[
K_{\mathcal C}(\rho,x)
=\frac{
 I-\mathcal C(\rho)\mathcal C^\sharp(x)
}{1-t}.                                          \tag{4}
\]

The complete retained resolvent is therefore

\[
\boxed{
D_r^{-1}
=R_R+\zeta t^rR_RW
 K_{\mathcal C}(\rho,x)^{-1}
 W^*R_R.}                                        \tag{5}
\]

This is the state-lifted form missing from L239: the two exterior
columns are retained without compression, and the middle factor is
exactly the inverse model kernel of the **delayed** transfer.

Moreover L221's delay splitting appears pointwise:

\[
\boxed{
K_{\mathcal C}(\rho,x)
=\sum_{j=0}^{r-1}t^jI+t^rK_{\mathcal B}(\rho,x).} \tag{6}
\]

Thus the physical continued fraction and the abstract Schur-model
delay flag are the same object, not merely analogous scalar weights.

L243 is not L228.  Formula (5) still has to be inserted into the
theta/ODE contour map and combined with the boundary metric and the
Stein right-defect Schur complement.  It does, however, remove the
last unbounded reflection remainder from that task.

## 2. Sequential Woodbury reduction

Apply Woodbury first to the right correction \(cT^*E\).  L238 gives

\[
V^*G^{-1}cT^*V=tI.
\]

Hence

\[
R_R
=G^{-1}
+G^{-1}cT^*V\frac1{1-t}V^*G^{-1}.                \tag{7}
\]

Using L238--L239's other three endpoint compressions,

\[
\begin{aligned}
W^*G^{-1}W
&=\frac{x}{1-t}(I-t\mathcal B\mathcal B^\sharp),\\
W^*G^{-1}cT^*V&=(1+t)\mathcal B,\\
V^*G^{-1}W&=x\mathcal B^\sharp,
\end{aligned}
\]

equation (7) gives

\[
\begin{aligned}
W^*R_RW
&=\frac{x}{1-t}(I-t\mathcal B\mathcal B^\sharp)
 +\frac{x(1+t)}{1-t}\mathcal B\mathcal B^\sharp\\
&=\frac{x}{1-t}(I+\mathcal B\mathcal B^\sharp),
\end{aligned}
\]

which proves (2), including the noncommutative order.

Now

\[
D_r=(G-cT^*E)-\delta_rWW^*.
\]

A second rank-\(m\) Woodbury step gives

\[
D_r^{-1}
=R_R+\delta_rR_RW
(I-\delta_rW^*R_RW)^{-1}W^*R_R.                  \tag{8}
\]

L239 proves

\[
\frac{\delta_rx}{1-t}=\frac{t^r}{1+t^r}.
\]

Substitute this and (2) into (8).  The endpoint denominator becomes

\[
I-\frac{t^r}{1+t^r}
 (I+\mathcal B\mathcal B^\sharp)
=\frac{I-t^r\mathcal B\mathcal B^\sharp}
       {1+t^r},
\]

proving (3).  Also

\[
\delta_r(1+t^r)=\zeta(1-t)t^r.
\]

Using (4) in (8) now proves (5).

Finally,

\[
\begin{aligned}
K_{\mathcal C}
&=\frac{I-t^r\mathcal B\mathcal B^\sharp}{1-t}\\
&=\frac{1-t^r}{1-t}I
 +t^r\frac{I-\mathcal B\mathcal B^\sharp}{1-t},
\end{aligned}
\]

which is (6).

## 3. Exact reflection-power filtration

Because

\[
\mathcal B(\rho)=\rho B_1+O(\rho^2),
\qquad
\mathcal B^\sharp(x)=xB_1^*+O(x^2),
\]

the nonconstant part of

\[
\{I-t^r\mathcal B(\rho)\mathcal B^\sharp(x)\}^{-1}
\]

starts at degree \(c^{r+1}\) on a fixed \(\zeta\)-contour.  Formula
(5) already has the exterior factor \(t^r\), of degree \(c^r\).
Consequently:

\[
\begin{array}{c|c}
\text{power of }t^r\mathcal B\mathcal B^\sharp
&\text{least total }c\text{-degree in (5)}\\ \hline
0&r\\
1&2r+1\\
p\ge2&(p+1)r+p\ge3r+2.
\end{array}                                       \tag{9}
\]

For every \(r\ge1\),

\[
3r+2>2r+2.
\]

Thus through the direct-map window required by L228, only the
zeroth and first kernel-reflection powers can occur.  This is the
filtration that L240 deliberately left open; retaining an infinite
endpoint Neumann tail is no longer necessary.  The first power must
still be combined with the Stein Schur square—it cannot simply be
discarded.

## 4. Independent audit

Run

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=experiments \
  .venv/bin/python -u \
  experiments/repeated_crabb_model_kernel_resolvent.py \
  --output \
  experiments/repeated_crabb_model_kernel_resolvent_s70224.jsonl
```

Six deterministic noncommuting cases, with delays one through six
and defect multiplicities one through three, verify (1)--(6) at
nine complex Joukowski points each.  The audit keeps both exterior
state columns in (5), rather than checking only endpoint
compressions.  The tracked SHA-256 is
`228ae5c557c4d76bbc68e7874243d7b72636f7c92bd62375248bb2c687e44576`.
