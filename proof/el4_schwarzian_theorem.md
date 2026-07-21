# EL4 via Schwarzian comparison (Epoch 6, 2026-07-21)

This note proves the explicit even-phase elliptic inequality `EL4` from
`slice_closed_form.md`.  It also isolates a reusable sufficient condition for the two-node
distortion inequality `D2`.

## 1. A Schwarzian sufficient condition for D2

For a real (C^3) map (G) with (G'>0), write

\[
 S G=\frac{G'''}{G'}-\frac32\left(\frac{G''}{G'}\right)^2.
\]

**Lemma (Schwarzian-D2).**  Let (-1<b<v<a<1), where (v) is the
pseudo-hyperbolic midpoint:

\[
 \delta=\frac{a-v}{1-av}=\frac{v-b}{1-bv}\in(0,1).
\]

If (S G\geq0) on ([b,a]), then

\[
 \frac{\delta}{2}(1-v^2)G'(v)
 \left(\frac1{G(a)-G(v)}+\frac1{G(v)-G(b)}\right)\leq1. \tag{D2}
\]

**Proof.**  Put (a=\tanh(t_0+D)), (v=\tanh t_0), and
(b=\tanh(t_0-D)).  Then (delta=\tanh D).  For
(g(t)=G(\tanh t)), the Schwarzian chain rule gives

\[
 Sg(t)=SG(\tanh t)\operatorname{sech}^4t-2\geq-2.
\]

Translate (t_0) to zero and postcompose (g) by an affine map so that
(g(0)=0), (g'(0)=1).  These operations do not change the desired inequality.  Set
(Q=Sg/2\geq-1).  The standard Schwarzian/ODE pair

\[
 y_1=\frac{g}{\sqrt{g'}},\qquad y_2=\frac1{\sqrt{g'}}
\]

solves (y''+Qy=0).  If (s,c) are the normalized solutions
(s(0)=0,s'(0)=1,c(0)=1,c'(0)=0), then (y_1=s) and (y_2=c+p s) for
some real (p).  Consequently

\[
 \frac1{g(D)}-\frac1{g(-D)}
 =\frac{c(D)}{s(D)}-\frac{c(-D)}{s(-D)}. \tag{1}
\]

The right side of (1) is the minimum of

\[
 \mathcal E_Q[h]=\int_{-D}^{D}(h'^2-Qh^2)\,dt
\]

over (h\in H^1_0(-D,D)) with (h(0)=1).  Indeed, the minimizer solves
(h''+Qh=0) on each half-interval, and integration by parts gives (1).
The positive solution (y_2=1/\sqrt{g'}) gives the ground-state identity, so the
Dirichlet form is positive.  Since (Q\geq-1), pointwise

\[
 \mathcal E_Q[h]\leq\mathcal E_{-1}[h].
\]

Taking the two minima and solving the constant-potential problem yields

\[
 \frac1{g(D)}-\frac1{g(-D)}\leq2\coth D.
\]

Multiplication by (	anh D/2) is exactly (D2).  (square)

The equality model is (SG=0), i.e. a real Möbius map.  This explains the
previously observed Möbius equality without treating it as an isolated identity.

## 2. The squared-ellipse map has nonnegative Schwarzian

Fix (0<k<1), put (m=k^2), (K=K(m)), and (s=\pi/(2K)).  On the real
focus interval define

\[
 Y(U)=k\,\operatorname{sn}^2(U,m),\qquad
 Z(U)=\sin^2(sU),\qquad 0<U<K,
\]

and let (G) be determined by (Z=G\circ Y).  This is the normalized inverse
map from the squared disk coordinate to the squared ellipse coordinate.

Direct differentiation, followed by the Jacobi duplication formula, gives

\[
 SZ=2s^2-\frac{6s^2}{\sin^2(2sU)},\qquad
 SY=2(1+m)-\frac6{\operatorname{sn}^2(2U,m)}. \tag{2}
\]

We need one classical elliptic-function inequality.  For (0<x<2K),

\[
 \frac1{\operatorname{sn}^2(x,m)}-
 \frac{s^2}{\sin^2(sx)}\geq\frac{1+m-s^2}{3}. \tag{3}
\]

Here is a short proof with no numerical input.  Let (wp) have half-periods
((K,iK')), and let (q=e^{-\pi K'/K}).  The standard identities are

\[
 \wp(x)=\frac1{\operatorname{sn}^2(x,m)}-\frac{1+m}{3}
\]

and the Fourier expansion (DLMF 23.8.1)

\[
 \wp(x)+\frac{\eta_1}{K}-s^2\csc^2(sx)
 =-8s^2\sum_{n\geq1}\frac{nq^{2n}}{1-q^{2n}}\cos(2nsx).
\]

Taking (x\to0) determines the constant term.  Subtraction gives

\[
 \frac1{\operatorname{sn}^2(x,m)}-\frac{s^2}{\sin^2(sx)}
 -\frac{1+m-s^2}{3}
 =8s^2\sum_{n\geq1}\frac{nq^{2n}}{1-q^{2n}}(1-\cos(2nsx))\geq0,
\]

which proves (3), strictly for (0<x<2K).

Apply (3) with (x=2U) in (2).  It gives (SZ-SY\geq0).  Since

\[
 SZ=(SG\circ Y)(Y')^2+SY
\]

and (Y'>0) on ((0,K)), we obtain

\[
 \boxed{SG(u)\geq0\quad(0<u<k).} \tag{4}
\]

## 3. EL4

Let (0<U_2<K),

\[
 a=k,\quad b=k\operatorname{sn}^2(U_2,m),
\]

and let (v=k\operatorname{sn}^2(W,m)) be their pseudo-hyperbolic midpoint.
Then (delta=(k-v)/(1-kv)=B_1).  Substitution of

\[
 G(k\operatorname{sn}^2U)=\sin^2(sU),\qquad
 G'(v)=\frac{s\sin(sW)\cos(sW)}{k\operatorname{sn}W\operatorname{cn}W\operatorname{dn}W}
\]

into (D2) gives exactly

\[
 \Theta=\frac{B_1s(1-v^2)\sin(sW)\cos^2(sU_2)}
 {2k\operatorname{sn}W\operatorname{cn}W\operatorname{dn}W\cos(sW)
  (\sin^2(sW)-\sin^2(sU_2))}\leq1.
\]

Therefore the explicit EL4 inequality is proved, and the even-phase formula from
`slice_closed_form.md` gives (ho=1-\Theta\geq0) for every nondegenerate
elliptic-slice configuration satisfying the midpoint hypotheses.  The degenerate edges follow
by continuity (with equality only at the corresponding disk/coalescence limits).

## 4. Scope and remaining logical debt

This proves the **analytic inequality EL4**, not the whole 4-by-4 case and not Crouzeix's
conjecture.  `proof/even_pick_globality.md` subsequently proved that the midpoint automorphism is
the global even-sector maximizer whenever that sector has norm larger than one, so EL4 now closes
the complete even sector of the elliptic slice.  Exclusion or treatment of symmetry-breaking
four-node extremals remains separate.  Odd and degree-one phases also remain.
