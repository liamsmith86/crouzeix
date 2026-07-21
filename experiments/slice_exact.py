"""Exact machinery for the elliptic sym4 slice A = S_a + c S_a^T (b_j = c a_j).

W(A) = ellipse, semi-axes A2 = sqrt(L(1+c)^2)/2, B2 = sqrt(L(1-c)^2)/2 where
L = lam_max-type = (T + sqrt(T^2 - 4 t1 t3))/2, T = sum a_j^2 — wait: from support
function h^2 = [(t1+t2+t3) + sqrt((sum)^2 - 4 t1 t3)]/2 * (u cos^2 + v sin^2)/1,
t_j = a_j^2, u = (1+c)^2/4, v = (1-c)^2/4. Thus
foci^2 = A2^2 - B2^2 = L*c = e1^2.

The foci are ±e1 (the outer eigenvalue pair), so f^2 = e1^2.  The Riemann map is
phi(z) = sqrt(k) sn(2K asin(z/f)/pi, k), K'(k)/K(k) = 4 xi0/pi.
Extremal odd Blaschke B(w) = w(w^2-alpha^2)/(1-alpha^2 w^2), alpha real (symmetry sector).
f0(A) = A F(A^2): K = max block norm; all 2x2 exact.

g0 = 1/f0 - PP, PP = r0/z + r+/(z-z1) + r-/(z+z1), z1 = psi(alpha),
r0 = 1/f0'(0), r± = 1/f0'(±z1) (f0 odd => f0' even => r+ = -r-? f0'(z) even,
f0'(-z1) = f0'(z1): r- = r+; note PP structure: r+/(z-z1)+r+/(z+z1) = 2r+ z/(z^2-z1^2)).
rho = Re<g0(A) f0(A) x0, x0> via block algebra with g0-values at {±e1, ±e2} (g0 odd).
"""
from mpmath import (mp, mpf, mpc, sqrt, atanh, asin, sin, cos, pi, ellipk, ellipf,
                    ellipfun, findroot, mpmathify)
import numpy as np

mp.dps = 25
sn = ellipfun('sn'); cn = ellipfun('cn'); dn = ellipfun('dn')

class Slice:
    def __init__(self, a1, a2, a3, c):
        self.a = [mpf(a1), mpf(a2), mpf(a3)]
        self.c = mpf(c)
        t1, t2, t3 = [x*x for x in self.a]
        T = t1+t2+t3
        L = (T + sqrt(T*T - 4*t1*t3))/2
        self.A2 = sqrt(L)*(1+self.c)/2
        self.B2 = sqrt(L)*(1-self.c)/2
        p = self.c*T; q = self.c**2*t1*t3
        disc = sqrt(p*p - 4*q)
        self.e1 = sqrt((p+disc)/2); self.e2 = sqrt((p-disc)/2)
        self.f = sqrt(self.A2**2 - self.B2**2)
        xi0 = atanh(self.B2/self.A2)
        target = 4*xi0/pi
        lo, hi = mpf(10)**-22, 1-mpf(10)**-22
        for _ in range(220):
            mid = (lo+hi)/2
            if ellipk(1-mid)/ellipk(mid) > target: lo = mid
            else: hi = mid
        self.m = (lo+hi)/2
        self.k = sqrt(self.m)
        self.K = ellipk(self.m)
        self.tau = [self.phi(self.e1), self.phi(self.e2)]

    def phi(self, z):
        u = (2*self.K/pi)*asin(z/self.f)
        return sqrt(self.k)*sn(u, self.m)

    def psi_data(self, w):
        """z = psi(w) and psi'(w) for real w in (-1,1)."""
        u = ellipf(asin(w/sqrt(self.k)), self.m)
        zn = sin(pi*u/(2*self.K))
        z = self.f*zn
        phip = sqrt(self.k)*cn(u,self.m)*dn(u,self.m)*(2*self.K/pi)/sqrt(1-zn**2)
        return z, self.f/phip*1/self.f*self.f  # psi'(w) = 1/phi'(z); phi'(z)=phip/f

    def blocks(self):
        a1,a2,a3 = self.a; c = self.c
        Ap = np.array([[float(a1),0],[float(c*a2),float(a3)]])
        Am = np.array([[float(c*a1),float(a2)],[0,float(c*a3)]])
        return Ap, Am

    def eval_alpha(self, alpha):
        """K and rho at odd-Blaschke parameter alpha (float pipeline on exact data)."""
        alpha = mpf(alpha)
        B = lambda w: w*(w*w-alpha*alpha)/(1-alpha*alpha*w*w)
        t1, t2 = self.tau
        e1, e2 = self.e1, self.e2
        F1, F2 = B(t1)/e1, B(t2)/e2   # F(e_j^2)
        Ap, Am = self.blocks()
        Ne = Am@Ap  # even-sublattice block of A^2
        lam, V = np.linalg.eig(Ne)
        Vi = np.linalg.inv(V)
        E = [float(e1**2), float(e2**2)]
        vals = [F1 if abs(l-E[0])<abs(l-E[1]) else F2 for l in lam]
        FN = V@np.diag([complex(v) for v in vals])@Vi
        M1 = Ap@FN
        # odd block
        No = Ap@Am
        lam2, V2 = np.linalg.eig(No)
        Vi2 = np.linalg.inv(V2)
        vals2 = [F1 if abs(l-E[0])<abs(l-E[1]) else F2 for l in lam2]
        M2 = Am@(V2@np.diag([complex(v) for v in vals2])@Vi2)
        K1, K2 = np.linalg.norm(M1,2), np.linalg.norm(M2,2)
        return max(K1, K2), (M1 if K1>=K2 else M2), K1>=K2

    def extremal(self):
        from scipy.optimize import minimize_scalar
        res = minimize_scalar(lambda x: -self.eval_alpha(x)[0],
                              bounds=(0.01, 0.99), method='bounded',
                              options={'xatol':1e-14})
        return res.x, -res.fun

    def rho_at(self, alpha):
        alpha = mpf(alpha)
        Kv, M, even_active = self.eval_alpha(alpha)
        # g0 odd: values at e1, e2. g0 = 1/f0 - PP.
        z1, dpsi = self.psi_data(alpha)
        B = lambda w: w*(w*w-alpha*alpha)/(1-alpha*alpha*w*w)
        # f0'(0): f0 = B(phi(z)): f0'(0) = B'(0) phi'(0); B'(0) = -alpha^2... B'(w) at 0:
        # B = w(w^2-a^2)/(1-a^2w^2): B'(0) = -alpha^2. phi'(0) = sqrt(k)(2K/pi)/f.
        phip0 = sqrt(self.k)*(2*self.K/pi)/self.f
        r0 = 1/((-alpha**2)*phip0)
        # f0'(z1): B'(alpha) phi'(z1): B'(alpha) = alpha(alpha^4-1)... compute:
        # B'(w) = [(3w^2-a^2)(1-a^2w^2) + 2a^2w^2(w^2-a^2)]/(1-a^2w^2)^2 at w=a:
        Bpa = (3*alpha**2-alpha**2)*(1-alpha**4)/(1-alpha**4)**2
        Bpa = 2*alpha**2/(1-alpha**4)
        # psi'(alpha) = dpsi
        r1 = dpsi/Bpa
        def g0(z, w):
            return 1/B(w) - r0/z - 2*r1*z/(z*z-z1*z1)
        g1 = g0(self.e1, self.tau[0]); g2 = g0(self.e2, self.tau[1])
        # rho = Re<(g0 f0)(A)x0,x0>: g0f0 even: H(e_j^2) = B(tau_j) g_j
        H1, H2 = B(self.tau[0])*g1, B(self.tau[1])*g2
        Ap, Am = self.blocks()
        N = (Am@Ap) if even_active else (Ap@Am)
        lam, V = np.linalg.eig(N); Vi = np.linalg.inv(V)
        E = [float(self.e1**2), float(self.e2**2)]
        vals = [complex(H1) if abs(l-E[0])<abs(l-E[1]) else complex(H2) for l in lam]
        HN = V@np.diag(vals)@Vi
        U, s, Vh = np.linalg.svd(M)
        x0 = Vh[0].conj()
        rho = float(np.real(x0.conj() @ (HN @ x0)))
        diag = abs(np.vdot(x0, M@x0))
        return Kv, rho, complex(g1), complex(g2), diag

if __name__ == "__main__":
    for (a1,a2,a3,c) in [(2.0,1.2,1.6,0.08),(1.8,1.0,1.4,0.15),(1.5,2.0,0.9,0.12)]:
        S = Slice(a1,a2,a3,c)
        al, Kv = S.extremal()
        K2, rho, g1, g2, diag = S.rho_at(al)
        print(f"a=({a1},{a2},{a3}) c={c}: alpha={al:.8f} K={K2:.6f} rho={rho:+.8f} "
              f"g(e1)={g1.real:+.6f}{g1.imag:+.1e}j g(e2)={g2.real:+.6f}{g2.imag:+.1e}j diag={diag:.1e}")
        print(f"   tau1={complex(S.tau[0]).real:.6f} tau2={complex(S.tau[1]).real:.6f} "
              f"e1={float(S.e1):.4f} e2={float(S.e2):.4f} f={float(S.f):.4f} k={float(S.k):.6f}")
