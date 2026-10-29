"""Native Hamada projection in NTV_pre::cal_xi_liu on a manufactured mode.

The unmodified production routine is compiled with stubbed displacement
input (a flux-constant m = 0 X1 whose stubbed transform is exp(i*k0*theta_j))
a zero perturbed field and a metric with xi_perp.grad B = X1, so that
dB(theta_j) = exp(i*k0*theta_j) on the periodic midpoint grid
theta_j = (j-1/2)*2*pi/n_chi that xgrid() supplies.  An exact periodic
projection returns bmn(k0) = 1 and bmn(k) = 0 for every other |k| <= m_max.
"""
import re

from support import ROOT, NativeTest


def ntv_routine(name, file="ntv_pre.f"):
    source = (ROOT / "MarsQ_2FK" / file).read_text()
    match = re.search(r"^ +subroutine +%s\b.*?^ +end subroutine +%s\b" % (name, name),
                      source, re.M | re.S | re.I)
    if match is None:
        raise ValueError("Native routine not found: " + name)
    return match.group(0)


LIB = """      module NTV_lib
      implicit none
      integer :: k0 = 0
      contains
      subroutine trans_1d(fc, fs, M0, M, n, n_chi, chi, dphi, fr, fi)
      integer :: M0, n, n_chi, M(M0), j
      real :: fc(M0), fs(M0), chi(n_chi), dphi(n_chi), fr(n_chi)
      real :: fi(n_chi), theta
      do j = 1, n_chi
         theta = (j - 0.5d0)*2*acos(-1d0)/n_chi
         fr(j) = fc(1)*cos(k0*theta) - fs(1)*sin(k0*theta)
         fi(j) = fs(1)*cos(k0*theta) + fc(1)*sin(k0*theta)
      end do
      end subroutine
%s
%s
      end module
"""

FIXTURE = """      module NTV_pre
      use NTV_lib
      implicit none
      real*8,parameter::npi=3.1415926535897932
      integer :: np_fixture = 3
      contains
      subroutine get_dimensV(MM1,NP1,n)
      integer MM1,NP1,n
      MM1 = 1
      NP1 = np_fixture
      n = 1
      end subroutine get_dimensV
      subroutine get_disp(M10,NP1,M,X1R,X1I,X2R,X2I,X3R,X3I)
      integer M10,NP1,M(M10)
      real X1R(NP1,M10),X1I(NP1,M10),X2R(NP1,M10),X2I(NP1,M10)
      real X3R(NP1,M10),X3I(NP1,M10)
      M = 0
      X1R = 1; X1I = 0; X2R = 0; X2I = 0; X3R = 0; X3I = 0
      end subroutine get_disp
      subroutine get_field(M10,NP1,B1R,B1I,B2R,B2I,B3R,B3I)
      integer M10,NP1
      real B1R(NP1,M10),B1I(NP1,M10),B2R(NP1,M10),B2I(NP1,M10)
      real B3R(NP1,M10),B3I(NP1,M10)
      B1R = 0; B1I = 0; B2R = 0; B2I = 0; B3R = 0; B3I = 0
      end subroutine get_field
%s
      end module NTV_pre
"""

ORACLE = """program projection
use NTV_lib, only: k0
use NTV_pre, only: cal_xi_liu
implicit none
integer, parameter :: np = 3, n_chi = 257, m_max = 80, n_max = 1
real :: chiH(np, n_chi), dphiH(np, n_chi), dBav(np), gg(np, n_chi, 3)
real :: dBds(np, n_chi), dBdc(np, n_chi), omega(np), dwds(np), err
real :: s(np), q(np), B(np, n_chi), dm(np, n_chi, 7), bth(np), cp(np)
complex :: bmn(2*m_max+1, n_max, np)
integer :: k, failures, modes(4)
chiH = 0; dphiH = 0; gg = 0; gg(:, :, 3) = 1; dBdc = 0; dBds = 1; omega = 0; dwds = 0
s = [0d0, 0.5d0, 1d0]; q = 0; B = 1; dm = 0; dm(:, :, 7) = 1; bth = 1; cp = 0
modes = [0, 1, 25, -80]
failures = 0
do k = 1, size(modes)
   k0 = modes(k)
   call cal_xi_liu(np, n_chi, chiH, 1, dphiH, dBav, gg, dBds, dBdc, &
                   omega, dwds, m_max, n_max, bmn, s, q, B, dm, bth, cp)
   bmn(k0+m_max+1, 1, :) = bmn(k0+m_max+1, 1, :) - 1
   err = maxval(abs(bmn(:, 1, 2:np)))
   print *, k0, err
   if (.not. err <= 1d-12) failures = failures + 1
end do
if (failures /= 0) error stop 1
end program
"""


class NtvProjection(NativeTest):
    def test_manufactured_mode_projects_exactly(self):
        lib = LIB % (ntv_routine("cald", "ntv_lib.f"),
                     ntv_routine("cald_1d", "ntv_lib.f"))
        exe = self.compile({"ntv_lib.f": lib,
                            "ntv_pre.f": FIXTURE % ntv_routine("cal_xi_liu"),
                            "oracle.f90": ORACLE})
        self.execute(exe)
