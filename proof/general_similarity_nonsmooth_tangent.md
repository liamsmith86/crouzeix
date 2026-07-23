# Nonsmooth conformal tangent at a repeated disk block (2026-07-22)

## 1. Statement

Let `K_epsilon` be convex bodies containing zero.  Suppose that, as
`epsilon` decreases to zero, their support functions satisfy

\[
 h_\epsilon(t)=1+\epsilon s(t)+r_\epsilon(t),\qquad
 \|r_\epsilon\|_\infty=o(\epsilon),                    \tag{1}
\]

where `s` is continuous and real.  Let

\[
 \Psi_\epsilon:\mathbb D\longrightarrow\mathop{\rm int}K_\epsilon
\]

be normalized by `Psi_epsilon(0)=0` and
`Psi_epsilon'(0)>0`.  If

\[
 \widehat s(k)=\frac1{2\pi}\int_0^{2\pi}
 s(t)e^{-ikt}\,dt,
\]

then, locally uniformly on the disk,

\[
 \boxed{
 \Psi_\epsilon(z)=z+\epsilon F_s(z)+o(\epsilon),\qquad
 F_s(z)=\widehat s(0)z+
 2\sum_{k\geq1}\widehat s(k)z^{k+1}.}                 \tag{2}
\]

No differentiability of `s`, strict convexity of `K_epsilon`, or simplicity
of the active support branch is required.  The result is sequentially
uniform when `s` ranges over a compact subset of `C(T)` and the remainder in
(1) is uniform.

## 2. Support functions give the radial variation

Let `rho_epsilon(t)` be the radial function of `K_epsilon`.  Convex duality
gives

\[
 \rho_\epsilon(t)=
 \min_{|u-t|<\pi/2}\frac{h_\epsilon(u)}{\cos(u-t)}.     \tag{3}
\]

The competitor `u=t` and (1) show that the minimum is at most `1+O(epsilon)`.
Every minimizer therefore has `|u-t|=O(sqrt(epsilon))`; outside such a
neighbourhood the denominator alone makes the quotient larger than the
competitor.  Uniform continuity of `s` now gives

\[
 s(u)=s(t)+o(1).
\]

Using `1/cos(u-t)>=1` in (3) for the lower bound and `u=t` for the upper
bound yields

\[
 \boxed{\rho_\epsilon(t)=1+\epsilon s(t)+o(\epsilon)}
                                                               \tag{4}
\]

uniformly in `t`.  The same proof is sequentially uniform for a compact,
hence equicontinuous, family of profiles.

## 3. Logarithmic boundary proof

Convexity and (4) make the boundaries Jordan curves converging uniformly to
the circle.  The normalized Riemann maps extend homeomorphically to the
closed disk and converge uniformly to the identity; this is the standard
Caratheodory--Rado convergence theorem for Jordan domains.

Write the boundary values as

\[
 \Psi_\epsilon(e^{it})
 =\rho_\epsilon(\tau_\epsilon(t))e^{i\tau_\epsilon(t)}.
\]

Uniform convergence gives `tau_epsilon(t)-t=o(1)`.  Since the only zero of
`Psi_epsilon` is at the origin,

\[
 L_\epsilon(z)=\log\frac{\Psi_\epsilon(z)}z
\]

has an analytic branch on the disk.  Choose it so that
`L_epsilon(0)` is real.  Equations (4) and uniform continuity give

\[
 \operatorname{Re}L_\epsilon(e^{it})
 =\log\rho_\epsilon(\tau_\epsilon(t))
 =\epsilon s(t)+o(\epsilon)                            \tag{5}
\]

uniformly.

The analytic function with real boundary data `s` and real value at zero is
its Schwarz integral

\[
 S_s(z)=\widehat s(0)+2\sum_{k\geq1}\widehat s(k)z^k.  \tag{6}
\]

Poisson integration of (5), followed by harmonic conjugation with the value
at zero fixed, gives

\[
 L_\epsilon(z)=\epsilon S_s(z)+o(\epsilon)             \tag{7}
\]

locally uniformly.  Exponentiating (7) proves (2).

The proof after (4) uses only radial Jordan domains.  Consequently the same
argument may be reapplied after a lower-order analytic recentering: if the
remaining radial displacement is `delta s+o(delta)`, its conformal response
is `delta F_s+o(delta)`, even when `s` is a largest-eigenvalue profile with
crossings.

## 4. Repeated Crabb consequence

Let

\[
 A_\epsilon=I_m\otimes C_p+\epsilon E+R_\epsilon,
 \qquad \|R_\epsilon\|=o(\epsilon).
\]

The top support eigenvalue of the base is one with multiplicity `m` and has
a uniform gap to the remaining physical levels.  Uniform finite-dimensional
degenerate perturbation theory therefore gives

\[
 h_\epsilon(t)=1+\epsilon\lambda_{\max}B_E(t)
 +O(\epsilon^2)+O(\|R_\epsilon\|),                    \tag{8}
\]

where

\[
 B_E(t)=V(t)^*\operatorname{Re}(e^{-it}E)V(t).
\]

For an exact linear path, or when `R_epsilon=O(epsilon^2)`, the last two
terms in (8) combine into a uniform `O(epsilon^2)` remainder.  For the
stated general path they give the `o(epsilon)` remainder required by (1).
The matrix `B_E(t)` is a trigonometric matrix polynomial.  Its largest
eigenvalue is continuous and uniformly Lipschitz in `t` on bounded sets of
directions; eigenvalue crossings do not weaken this conclusion.  Apply (2)
with

\[
 s_E(t)=\lambda_{\max}B_E(t).                           \tag{9}
\]

Inverting the locally uniform map expansion on compact subsets and using
holomorphic functional calculus gives

\[
 \boxed{
 \phi_\epsilon(A_\epsilon)
 =A_0+\epsilon\{E-F_{s_E}(A_0)\}+o(\epsilon).}         \tag{10}
\]

Since `A_0^p=0`, only Fourier coefficients `0,...,p-2` in (2) survive in
`F_{s_E}(A_0)`.  Equation (10) is exactly the conformal expansion used by
L61, now with a uniform proof that permits arbitrary support-branch
crossings.

## 5. Scope

This removes the regularity hypothesis from L61 and provides the
shape-derivative tool needed at later weighted blow-up scales.  It does not
by itself prove the repeated-`C3` neighbourhood theorem: one must still show
that the finite hierarchy of negative endpoint forms absorbs all recentered
remainders near the exact normal manifold and then merge the result through
the L93 flag.
