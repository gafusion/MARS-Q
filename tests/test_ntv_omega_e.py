"""E x B frequency of the Shaing NTV profiles (NTV_pre::cal_nu).

Both branches of cal_nu deliver nu_prof(:,8), the E x B frequency that gntv
adds to the diamagnetic frequency nu_prof(:,9) and that sets the resonance
x_min = nu_prof(:,8)/nu_prof(:,7).  Without PROFWE (KROTWE /= 4) it is the
radial force balance w_phi - knc*qwsTi + qwsi in rad/s; with PROFWE.IN
(KROTWE = 4) it is the toroidal E x B rotation frequency read by marsq.f
into ROTWEI, the same rad/s quantity that MARS-K uses as OMEGAE0.  For a
PROFWE equal to the force-balance value both branches must therefore return
the same frequency, at every safety factor.
"""
import re

from support import ROOT, NativeTest

SRC = ROOT / "MarsQ_2FK" / "ntv_pre.f"


def routine(name):
    match = re.search(r"^ +subroutine +%s\b.*?^ +end subroutine +%s\b" % (name, name),
                      SRC.read_text(), re.M | re.S | re.I)
    if match is None:
        raise ValueError("Native routine not found: " + name)
    return match.group(0) + "\n"


PROGRAM = """      program omega_e
      implicit none
      integer, parameter :: nr = 4, np = 14
      real :: r(nr), drdpsi(nr), rp(nr,np), out0(nr,np), out4(nr,np)
      integer :: i
      r = [0.2, 0.4, 0.6, 0.8]
      drdpsi = [0.9, 0.7, 0.5, 0.4]
      rp = 0.0
      rp(:,1) = r
      rp(:,2) = [1.1, 1.9, 3.2, 4.7]
      rp(:,3) = [0.9, 0.8, 0.6, 0.4]
      rp(:,4) = [9.0, 7.0, 4.0, 1.5]
      rp(:,5) = [8.0, 6.5, 3.5, 1.2]
      rp(:,6) = [2.0e4, -1.5e4, 3.0e3, 8.0e3]
      rp(:,7) = [-0.1, -0.3, -0.6, -1.1]
      rp(:,8) = [-3.0, -6.0, -9.0, -7.0]
      rp(:,9) = [-2.5, -5.0, -8.0, -6.0]
      rp(:,10) = 0.3*r
      rp(:,13) = -1.0
      rp(:,14) = -1.0
      call cal_nu(0, nr, r, drdpsi, np, rp, 2.0, 6.2, 5.3, 1, 1, 1.17, out0)
      rp(:,12) = out0(:,8)
      call cal_nu(4, nr, r, drdpsi, np, rp, 2.0, 6.2, 5.3, 1, 1, 1.17, out4)
      do i = 1, nr
         write(*,'(3es24.16)') rp(i,2), out0(i,8), out4(i,8)
      end do
      end program omega_e
"""


class NtvOmegaE(NativeTest):
    def test_profwe_branch_equals_force_balance_frequency(self):
        exe = self.compile({"cal_nu.f": routine("cal_nu"), "main.f90": PROGRAM})
        rows = [list(map(float, line.split()))
                for line in self.execute(exe).stdout.split("\n") if line.strip()]
        self.assertEqual(len(rows), 4)
        for q, balance, profwe in rows:
            self.assertGreater(abs(balance), 1.0e2, q)
            self.assertLess(abs(profwe - balance), 1.0e-12 * abs(balance), q)


if __name__ == "__main__":
    import unittest
    unittest.main()
