# Ambient first-order stationarity of the non-Crabb Gau--Wu family

> **Campaign scope.**  This is the first local calculation at the
> nonnilpotent equality matrices L330 exposed.  It proves exact
> stationarity after numerical-range conformal normalization.  It
> does not determine the second-order sign or prove a neighbourhood
> theorem.

## 1. Result (L331, 2026-07-26)

Fix real \(a\), \(0<|a|<1\), and let \(G_a,f_a\) be L330's disk
equality pair.  For an arbitrary complex matrix direction \(E\), let

\[
 s_E(\theta)
 =x_\theta^*
 \operatorname{Re}(e^{-i\theta}E)x_\theta,        \tag{1}
\]

where \(x_\theta\) is the unit top-support vector of \(G_a\).
Write

\[
 s_E(\theta)=\sum_{k\in\mathbb Z}\widehat s_k e^{ik\theta},
\qquad
 H_s(z)=\widehat s_0+2\sum_{k\ge1}\widehat s_kz^k. \tag{2}
\]

L101's Riemann-map tangent gives the normalized operator direction

\[
 \boxed{\mathcal G_E=E-G_aH_s(G_a).}              \tag{3}
\]

Then the sharp scalar branch is stationary in every ambient
direction:

\[
 \boxed{
 {d\over d\varepsilon}\bigg|_{\varepsilon=0}
 \left\|
 f_a\!\left(
 \Phi_{G_a+\varepsilon E}(G_a+\varepsilon E)
 \right)
 \right\|^2=0.}                                   \tag{4}
\]

Here the Riemann maps are centered and phased as in L101.  Since
\(f_a(G_a)=2E_{13}\), the top singular value is simple and no
nonsmooth norm selection is involved.

L331 generalizes L162's ambient stationarity mechanism from the
monomial Crabb collision to the first genuinely nonnilpotent
Gau--Wu stratum.  The next live calculation is the transverse
second variation modulo the tangent of the exact disk-model
manifold.

## 2. The functional-calculus gradient

Put

\[
 {\cal L}_a(E)
 =e_1^*Df_a(G_a)[E]e_3.                           \tag{5}
\]

Differentiating

\[
 f_a(X)=X(X-aI)(I-aX)^{-1}
\]

gives

\[
\boxed{
 {\cal L}_a(E)
 =r_aE_{12}-aE_{13}+2aE_{22}+r_aE_{23},
 \qquad r_a=\sqrt{2(1-a^2)}.}                    \tag{6}
\]

The squared-norm derivative in (4) is therefore

\[
 4\operatorname{Re}{\cal L}_a(\mathcal G_E).      \tag{7}
\]

The minimal polynomial \(z^2(z-a)\) gives

\[
 {\cal L}_a(G_a)=4,\qquad
 {\cal L}_a(G_a^{k+1})=2a^k\quad(k\ge1).          \tag{8}
\]

Consequently

\[
\operatorname{Re}{\cal L}_a(G_aH_s(G_a))
=4\widehat s_0+
4\operatorname{Re}\sum_{k\ge1}a^k\widehat s_k.   \tag{9}
\]

## 3. The support-projection identity

Let \(\zeta=e^{i\theta}\) and

\[
 K_\theta=\operatorname{Re}(\zeta^{-1}G_a).
\]

L330 gives its eigenvalues \(1,-1,a\cos\theta\).  Hence its
rank-one top projection is exactly

\[
 P_+(\zeta)
 ={(K_\theta+I)(K_\theta-a\cos\theta\,I)
 \over2(1-a\cos\theta)}.                          \tag{10}
\]

Let

\[
 {\cal P}_a(\theta)
 ={1-a^2\over1-2a\cos\theta+a^2}                 \tag{11}
\]

be the Poisson kernel.  Since \(s_E\) is real,

\[
\widehat s_0+
{\cal P}[s_E](a)
=2\widehat s_0+
2\operatorname{Re}\sum_{k\ge1}a^k\widehat s_k.   \tag{12}
\]

The load-bearing contour identity is

\[
\boxed{
 2\,\mathop{\rm mean}_\theta
 \{1+{\cal P}_a(\theta)\}
 \zeta^{-1}P_+(\zeta)^T
 =
 \begin{bmatrix}
 0&r_a&-a\\
 0&2a&r_a\\
 0&0&0
 \end{bmatrix}.}                                 \tag{13}
\]

To evaluate it without nested radicals, put

\[
 a={2q\over1+q^2},\qquad
 r_a={\sqrt2(1-q^2)\over1+q^2},\qquad |q|<1.
\]

Every entry on the left of (13) is rational in \((q,\zeta)\).
Its circle mean is the sum of the residues at
\(\zeta=0,q,a\); entrywise reduction gives the displayed matrix.

Taking the real pairing with \(E\) in (13) proves

\[
\operatorname{Re}{\cal L}_a(E)
=4\widehat s_0+
4\operatorname{Re}\sum_{k\ge1}a^k\widehat s_k.   \tag{14}
\]

Equations (9) and (14) imply
\(\operatorname{Re}{\cal L}_a(\mathcal G_E)=0\).
Substitution into (7) proves (4).

## 4. What remains

First-order stationarity is necessary but not sufficient.  The
neighbourhood theorem now requires:

1. identify the tangent space of the exact \(G_a\) disk/Gau--Wu
   equality manifold, including affine and unitary directions;
2. compute the optimized second variation on a complementary normal
   slice, allowing the Blaschke zero to move; and
3. prove strict negativity off the tangent, or identify the next
   flat mode before attempting higher order.

The uniform support gap \(1-|a|\) makes all these coordinates
ordinary analytic coordinates for fixed \(a\); no repeated-support
ramification is needed.

## 5. Audit

Run

```bash
PYTHONPATH=experiments .venv/bin/python -u \
  experiments/gau_wu_ambient_stationarity.py \
  --output experiments/gau_wu_ambient_stationarity_s70224.jsonl
```

The checker derives (6), (8), and all nine residue entries in (13)
with \(q\) left symbolic, then tests the resulting real stationarity
identity on deterministic rational complex directions.
