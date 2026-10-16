"""Native GETXYMORE must invert velocity on each component's radial mesh.

Run with python -m unittest discover -s tests -p test_getxymore_halfmesh.py.
Only the native modules and the verbatim GETXYMORE routine are compiled,
using the native real8/double8 precision flags. Nonzero poloidal harmonics
also exercise the Fourier angle step when both tangential components are derived.
The manufactured velocity follows v = (gamma + i*n*Omega)*xi minus the
rotation-shear term; integer and half meshes use the same analytic profile.
"""

from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
DRIVER = r"""
program check_getxymore
  use dimensim
  use globalm
  use rcomdm
  implicit none
  complex*16 :: al0,alnorm,alam,atau,altau,omega
  real*8 :: epspam,epsdet
  integer :: ncase,nitmax,noncon
  common /pamarg/ al0,alnorm,alam,atau,altau,omega,epspam,epsdet, &
                  ncase,nitmax,noncon
  real*8 :: nodes(4),half(3),b2,error(3),r
  complex*16 :: exact1(4,3),exact2(3,3),exact3(3,3),xi_half,doppler
  integer :: i,m,k
  character(32) :: label
  nr=3; nrp1=4; nchi=16; msmax=3
  allocate(rot(nrp1),rotm(nrp1),drotm(nrp1),rm(msmax,2))
  allocate(v1u(nrp1,msmax),x1u(nrp1,msmax))
  allocate(v2u(nr,msmax),x2u(nr,msmax),v3u(nr,msmax),x3u(nr,msmax))
  allocate(rjam(nr,nchi),g22lm(nr,nchi),reqm(nr,nchi),tm(nr),dpsidsm(nr))
  nodes=[0.d0,.17d0,.61d0,1.d0]
  half=(nodes(1:nr)+nodes(2:nrp1))/2.d0
  rm=0.d0; rm(:,2)=[-1.d0,0.d0,2.d0]
  alnorm=cmplx(.4d0,.7d0,kind=8); rntor=-3.d0; rote=1.d0
  do k=1,6
    kxx1=-1; kyx2=-1; kyx3=-1
    do i=1,nrp1
      r=nodes(i)
      rot(i)=.8d0+.3d0*r+.2d0*r*r
      if(k==3)rot(i)=.8d0
      do m=1,msmax
        exact1(i,m)=cmplx(.2d0+.1d0*r+.07d0*m,-.03d0*m,kind=8)
        if(k==1)exact1(i,m)=0.d0
        v1u(i,m)=(alnorm+ci*rntor*rot(i))*exact1(i,m)
      end do
    end do
    do i=1,nr
      r=half(i)
      rotm(i)=.8d0+.3d0*r+.2d0*r*r
      drotm(i)=.3d0+.4d0*r
      if(k==3)then
        rotm(i)=.8d0
        drotm(i)=0.d0
      end if
      tm(i)=1.2d0+.1d0*r; dpsidsm(i)=.7d0+.2d0*r
      rjam(i,:)=1.1d0+.3d0*r
      reqm(i,:)=2.d0+.1d0*r
      g22lm(i,:)=.9d0+.2d0*r
      b2=g22lm(i,1)*dpsidsm(i)**2/rjam(i,1)**2+tm(i)**2/reqm(i,1)**2
      doppler=alnorm+ci*rntor*rotm(i)
      do m=1,msmax
        exact2(i,m)=cmplx(.9d0+.2d0*r+.1d0*m,-.11d0*m,kind=8)
        exact3(i,m)=cmplx(-.3d0+.15d0*r,.2d0+.09d0*m,kind=8)
        xi_half=(exact1(i,m)+exact1(i+1,m))/2.d0
        v2u(i,m)=doppler*exact2(i,m)+dpsidsm(i)*drotm(i)*xi_half/rjam(i,1)
        v3u(i,m)=doppler*exact3(i,m)-tm(i)*drotm(i)*xi_half/b2
      end do
    end do
    x1u=exact1; x2u=exact2; x3u=exact3
    if(k==4)then
      ! Directly solved displacements must be retained verbatim.
      kxx1=1; kyx2=1; kyx3=1
      v1u=0.d0; v2u=0.d0; v3u=0.d0
    else
      x1u=0.d0; x2u=0.d0; x3u=0.d0
    end if
    if(k==5)then
      kyx2=1; x2u=exact2
    else if(k==6)then
      kyx3=1; x3u=exact3
    end if
    call getxymore
    error=[maxval(abs(x1u-exact1)),maxval(abs(x2u-exact2)),maxval(abs(x3u-exact3))]
    select case(k)
    case(1);label='differential_rotation_xi1_zero'
    case(2);label='differential_rotation_with_shear'
    case(3);label='uniform_rotation_control'
    case(4);label='directly_solved_control'
    case(5);label='only_x3_reconstructed'
    case(6);label='only_x2_reconstructed'
    end select
    write(*,'(A,3ES25.16E3)')trim(label),error
  end do
end program
"""


class GetXYMoreHalfMesh(unittest.TestCase):
    def test_native_velocity_displacement_roundtrip(self):
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is required for the native routine test")
        source = ROOT / "MarsQ_2FK"
        text = (source / "marsq.f").read_text()
        match = re.search(r"^      SUBROUTINE GETXYMORE\s*$", text, re.MULTILINE)
        self.assertIsNotNone(match)
        native = text[match.start() : text.index("*DECK DPWALL", match.end())]
        with tempfile.TemporaryDirectory(prefix="getxymore-") as directory:
            work = Path(directory)
            (work / "getxymore.f").write_text(native)
            (work / "driver.f90").write_text(DRIVER)
            cmd = [
                compiler,
                "-O0",
                "-g",
                "-fcheck=all",
                "-ffpe-trap=invalid,zero,overflow",
                "-ffixed-line-length-none",
                "-fdefault-real-8",
                "-fdefault-double-8",
                "-I",
                str(source),
                "-J",
                str(work),
                str(source / "dimensim.f"),
                str(source / "globalm.f"),
                str(source / "rcomdm.f"),
                str(work / "getxymore.f"),
                str(work / "driver.f90"),
                "-o",
                str(work / "check"),
            ]
            built = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            result = subprocess.run(
                [str(work / "check")], cwd=work, capture_output=True, text=True
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        records = [line.split() for line in result.stdout.splitlines()]
        self.assertEqual(len(records), 6, result.stdout)
        for label, *numbers in records:
            with self.subTest(case=label):
                errors = list(map(float, numbers))
                self.assertEqual(len(errors), 3)
                self.assertLess(max(errors), 2.0e-13, result.stdout)


if __name__ == "__main__":
    unittest.main()
