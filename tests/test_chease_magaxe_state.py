#!/usr/bin/env python3
"""Exercise production MAGAXE state with an independent analytic minimum."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import support


INCLUDES = {
    "DECLAR.inc": """      REAL EPNON0,EPSMCH,EZMAG,RC0P,RC1P,RC1M12,RC1M13
      REAL RC1M14,RC2P,RC2PI
      COMMON /COMDEC/ EPNON0,EPSMCH,EZMAG,RC0P,RC1P,RC1M12,RC1M13,RC1M14,RC2P,RC2PI
""",
    "COMDIM.inc": "",
    "COMBLA.inc": "",
    "COMCON.inc": "      INTEGER NSYM\n      COMMON /COMCON/ NSYM\n",
    "COMESH.inc": "      REAL R0,RZ0\n      INTEGER NMAG\n      COMMON /COMESH/ R0,RZ0,NMAG\n",
    "COMNUM.inc": "",
    "COMPHY.inc": "      REAL ASPCT,RMAG,RZMAG,SPSIM\n      COMMON /COMPHY/ ASPCT,RMAG,RZMAG,SPSIM\n",
    "COMSOL.inc": "      REAL RELAX\n      COMMON /COMSOL/ RELAX\n",
}

DRIVER = """program magaxe_state
  implicit none
  include 'DECLAR.inc'
  include 'COMCON.inc'
  include 'COMESH.inc'
  include 'COMPHY.inc'
  include 'COMSOL.inc'
  real :: target_r, target_z
  integer :: levels(6)
  character(len=16) :: mode
  common /analytic_target/ target_r,target_z
  call get_command_argument(1,mode)
  levels = [1,1,1,2,2,3]
  if (trim(mode) == 'reset') levels = [1,2,3,4,4,5]
  rc1m14 = 1.e-10
  rc1m12 = 100.*rc1m14
  ezmag = 1.e-4
  r0 = 1.
  rz0 = 0.
  aspct = .3125
  target_r = 1.03125
  target_z = 1./64.
  rmag = target_r
  rzmag = 0.
  spsim = 0.
  nmag = 0
  nsym = 0
  relax = 0.
  call magaxe
  call check(1,target_z)
  nmag = 0
  call magaxe
  call check(2,target_z)
  target_z = target_z + 1./65536.
  nmag = 0
  call magaxe
  call check(3,target_z)
  target_z = target_z + 1./1024.
  nmag = 0
  call magaxe
  call check(4,target_z)
  nsym = 1
  target_z = 0.
  nmag = 0
  call magaxe
  call check(5,0.)
  target_z = 1./64.
  nmag = 0
  call magaxe
  call check(6,0.)
  print *, 'PASS native MAGAXE state'
contains
  subroutine check(stage,accepted_z)
    integer, intent(in) :: stage
    real, intent(in) :: accepted_z
    if (abs(rmag-target_r) > 4.e-8 .or. &
        abs(rzmag-accepted_z) > 4.e-8 .or. abs(spsim+1.) > 1.e-13) then
      print *, 'axis mismatch',stage,rmag,rzmag,spsim,accepted_z
      stop 1
    end if
    if (abs(relax-.01*levels(stage)) > 1.e-14) then
      print *, 'RELAX mismatch',stage,relax,.01*levels(stage)
      stop 2
    end if
    print '(A,I1,3ES24.15)', 'state ',stage,relax,rmag,rzmag
  end subroutine
end program

subroutine evlate(kcase,pr,pz,pdpdr,pdpdz,psi)
  implicit none
  integer, intent(in) :: kcase
  real, intent(in) :: pr,pz
  real, intent(out) :: pdpdr,pdpdz,psi
  real :: target_r,target_z
  real, parameter :: curvature = 1./256.
  common /analytic_target/ target_r,target_z
  if (kcase /= 1 .and. kcase /= 2) stop 3
  pdpdr = curvature*(pr-target_r)
  pdpdz = curvature*(pz-target_z)
  psi = -1. + .5*curvature*((pr-target_r)**2+(pz-target_z)**2)
end subroutine

subroutine rvar(name,value)
  implicit none
  character(len=*), intent(in) :: name
  real, intent(in) :: value
  print *, 'unexpected MAGAXE failure',name,value
  stop 4
end subroutine
"""


def run_oracle(source_dir=None, expect_reset=False):
    source_dir = Path(source_dir or support.NATIVE).resolve()
    with patch.object(support, "NATIVE", source_dir):
        native = support.routine("MAGAXE")
    with tempfile.TemporaryDirectory(prefix="chease-magaxe-",
                                     dir=os.environ.get("TMPDIR")) as directory:
        work = Path(directory)
        for name, content in INCLUDES.items():
            (work / name).write_text(content)
        (work / "native.f").write_text(native)
        (work / "oracle.f90").write_text(DRIVER)
        # Match the production GNU optimization and initialization semantics.
        flags = ["-O2", "-fdefault-real-8", "-fdefault-double-8",
                 "-ffixed-line-length-none", "-fallow-argument-mismatch",
                 "-finit-local-zero", "-mcmodel=medium", "-std=legacy",
                 "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
                 "-fbacktrace"]
        build = subprocess.run([os.environ.get("FC", "gfortran"), *flags,
                                "native.f", "oracle.f90", "-o", "oracle"],
                               cwd=work, capture_output=True, text=True, timeout=90)
        if build.returncode:
            raise AssertionError(build.stdout + build.stderr)
        command = [str(work / "oracle")]
        if expect_reset:
            command.append("reset")
        outputs = []
        for _ in range(2):
            result = subprocess.run(command, cwd=work, capture_output=True,
                                    text=True, timeout=30)
            if result.returncode:
                raise AssertionError(result.stdout + result.stderr)
            outputs.append(result.stdout)
        if outputs[0] != outputs[1]:
            raise AssertionError("fresh-process MAGAXE state differs\n" + "\n".join(outputs))
        return outputs[0]


class TestCheaseMagaxeState(unittest.TestCase):
    def test_native_axis_state_persists_between_calls(self):
        self.assertIn("PASS native MAGAXE state", run_oracle())


if __name__ == "__main__":
    if "--source-dir" in sys.argv or "--expect-reset" in sys.argv:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--source-dir", type=Path)
        parser.add_argument("--expect-reset", action="store_true")
        args = parser.parse_args()
        print(run_oracle(args.source_dir, args.expect_reset), end="")
    else:
        unittest.main()
