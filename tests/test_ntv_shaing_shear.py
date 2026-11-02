"""Magnetic shear in the trapped-particle precession of the Shaing NTV route.

The bounce-averaged toroidal precession of a trapped particle is the radial
derivative of its parallel action.  In a large-aspect-ratio circular tokamak,
B = B0(1 - eps cos theta), J = q R0 * loop integral of v_par dtheta, so the
derivative carries q'(r) and the pitch dependence is

    G(kappa^2) = 2E/K - 1 + 4 s (E/K + kappa^2 - 1),   s = (r/q) dq/dr

(Connor, Hastie and Taylor).  The superbanana-plateau resonance of the Shaing
route is x*wB0*G(kappa0^2) = -q*omega_E.  For |s| > 1/2, G is not monotone and
a pitch can resonate twice.

The oracle is independent of the production code: scipy roots of G, bounce
averages of b(theta) by adaptive quadrature and a numerical dG/dkappa^2.  It
shares only the energy nodes of gntv1d (x_j = (j-1/2)*15/101).  The
collisionality is negligible, so the resonant torque is the
superbanana-plateau term alone.
"""
import math
import re
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ellipe, ellipk

from support import ROOT, NativeTest

NATIVE = ROOT / "MarsQ_2FK"
Q, EPS, WB0, QWE, RHO_I, R_AXIS, WTI = 1.7, 0.2, 1.0, -2.0, 1.0, 6.2, 1.0
B_M = {1: 1.0e-3 + 0.0j, 2: 0.0 + 4.0e-4j}       # b_mn for n = 1
SHEARS = (0.0, 0.3, 1.5)

STUB = """module NTV_pre
contains
subroutine MARS_out(nr,n_chi,n_phi,n_prof,keyntv,rho,mui,Z,Zi,knc,m_max,n_max,R_axis,nu_prof,bc,bs,bc2,bs2)
integer nr,n_chi,n_phi,n_prof,keyntv,Z,Zi,m_max,n_max
real rho(nr),mui,knc,R_axis,nu_prof(nr,n_prof),bc(*),bs(*),bc2(*),bs2(*)
stop 'MARS_out is not used'
end subroutine
end module
"""

DRIVER = f"""program shaing_shear
use common_ntv
use NTV
implicit none
integer :: i
real*8 :: s(4)
s = [0d0, 0.25d0, 0.5d0, 0.75d0]
! Fifteen columns: the shear column is ignored where it is not read.
allocate(nu_prof(3,15))
call basic_input(1, 3, 32, 1, 0, 1, s)
nu_prof(:,2) = {Q}d0
nu_prof(:,3) = {EPS}d0
nu_prof(:,4) = {RHO_I}d0
nu_prof(:,5) = 1d-12
nu_prof(:,6) = {WTI}d0
nu_prof(:,7) = {WB0}d0
nu_prof(:,8) = {QWE}d0
nu_prof(:,15) = [{', '.join(f'{x}d0' for x in SHEARS)}]
{chr(10).join(f'bc({m}+m_max+1,1,:) = {v.real!r}d0; bs({m}+m_max+1,1,:) = {v.imag!r}d0'
              for m, v in B_M.items())}
bc2 = bc
bs2 = bs
R_axis = {R_AXIS}d0
call gntv1d
do i = 1, 3
   write(*,'(2es25.16)') dble(T_ntv_res(i)), dble(T_ntv_non(i))
end do
end program shaing_shear
"""


def precession(m, s):
    e = ellipe(m)/ellipk(m)
    return 2.0*e - 1.0 + 4.0*s*(e + m - 1.0)


def bounce_average(m):
    def b(theta):
        return sum(v*np.exp(1j*(mm - Q)*theta) for mm, v in B_M.items())

    def part(phi, k):
        value = b(2.0*math.asin(math.sqrt(m)*math.sin(phi)))/math.sqrt(1.0 - m*math.sin(phi)**2)
        return value.real if k == 0 else value.imag
    re_, im_ = (quad(part, -math.pi/2, math.pi/2, args=(k,), limit=400, epsabs=0,
                     epsrel=1e-11)[0] for k in (0, 1))
    return complex(re_, im_)/(2.0*ellipk(m))


def resonant_torque(s):
    grid = np.concatenate([np.linspace(0.0, 0.99, 2000), 1.0 - np.logspace(-2, -15, 2000)])
    values = precession(grid, s)
    dx = 15.0/101
    lam4 = 0.0
    for j in range(101):
        x = (j + 0.5)*dx
        target = QWE/WB0/x*-1.0
        k_sbp = 0.0
        for a, b, fa, fb in zip(grid, grid[1:], values, values[1:]):
            if (fa > target) == (fb > target):
                continue
            m = brentq(lambda u: precession(u, s) - target, a, b, xtol=1e-15, rtol=1e-14)
            h = 1e-6*min(m, 1.0 - m)
            slope = (precession(m + h, s) - precession(m - h, s))/(2.0*h)
            k_sbp += 4.0*math.pi*ellipk(m)/abs(WB0*x*slope)*abs(bounce_average(m))**2
        lam4 += 0.5*k_sbp*x**2.5*math.exp(-x)*dx
    t0 = -0.5/math.sqrt(2.0)/math.pi**1.5*RHO_I*R_AXIS**2*math.sqrt(EPS)*Q**2*WTI**2
    return t0*lam4*QWE


class ShaingShear(NativeTest):
    def test_resonant_torque_follows_the_sheared_precession(self):
        sources = {name: (NATIVE / name).read_text()
                   for name in ("ntv_lib.f", "splines.f", "common_ntv.f", "ntv.f")}
        exe = self.compile({"ntv_lib.f": sources["ntv_lib.f"], "splines.f": sources["splines.f"],
                            "stub.f90": STUB, "common_ntv.f": sources["common_ntv.f"],
                            "ntv.f": sources["ntv.f"], "main.f90": DRIVER})
        rows = [list(map(float, line.split()))
                for line in self.execute(exe).stdout.splitlines() if line.strip()]
        self.assertEqual(len(rows), 3)
        for s, (native, _) in zip(SHEARS, rows):
            want = resonant_torque(s)
            self.assertLess(abs(native/want - 1.0), 1e-4, (s, native, want))


SOURCE = (NATIVE / "ntv_pre.f").read_text()
LIBRARY = (NATIVE / "ntv_lib.f").read_text()


def routine(source, name):
    match = re.search(r"^ +subroutine +%s\b.*?^ +end subroutine +%s\b" % (name, name),
                      source, re.M | re.S | re.I)
    if match is None:
        raise ValueError("Native routine not found: " + name)
    return match.group(0) + "\n"


PROFILE = """      program shear_profile
      implicit none
      integer, parameter :: nr = 6, np = 15
      real :: r(nr), drdpsi(nr), rp(nr,np), out(nr,np), a(nr)
      real, parameter :: R0 = 6.2, B0 = -5.3, q0 = 1.1, q2 = 4.0
      integer :: i
c     circular large-aspect-ratio surfaces: minor radius a, eps = a/R0,
c     q = q0 + q2 a**2, dpsi/da = a B0/q; the profile variable is r = a**2
      a = [0.4, 0.7, 1.0, 1.3, 1.6, 1.9]
      r = a**2
      rp = 0.0
      rp(:,1) = r
      rp(:,2) = q0 + q2*a**2
      rp(:,3) = a/R0
      rp(:,4) = 5.0
      rp(:,5) = 5.0
      rp(:,10) = a/R0
      rp(:,13) = -1.0
      rp(:,14) = -1.0
      drdpsi = 2.0*a*rp(:,2)/(a*abs(B0))
      call cal_nu(0, nr, r, drdpsi, np, rp, 2.0, R0, B0, 1, 1, 1.17, out)
      do i = 2, nr-1
         write(*,'(3es24.16)') a(i), rp(i,2), out(i,15)
      end do
      end program shear_profile
"""


class ShaingShearProfile(NativeTest):
    def test_cal_nu_shear_is_r_dlnq_dr(self):
        exe = self.compile({"cal_nu.f": routine(SOURCE, "cal_nu"),
                            "cald_1d.f": routine(LIBRARY, "cald_1d"), "main.f": PROFILE})
        rows = [list(map(float, line.split()))
                for line in self.execute(exe).stdout.splitlines() if line.strip()]
        self.assertEqual(len(rows), 4)
        for a, q, shear in rows:
            want = a*(2.0*4.0*a)/q
            self.assertLess(abs(shear/want - 1.0), 1e-12, (a, shear, want))


if __name__ == "__main__":
    unittest.main()
