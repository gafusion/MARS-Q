"""Native fast/slow trapped-bounce kernels must honor the same physical drift switch."""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get("MARS_ORACLE_SOURCE", Path(__file__).resolve().parents[1]))
NATIVE = ROOT / "MarsQ_2FK"

DRIVER = """program bounce_parity
use rcomdm
use dimensim
use globalm
use kineticm
use anisotropicm
use ToolBox
implicit none
integer i,j,mode,ndb,grid,signum
real*8 center,halfwidth,umax,error,largest
complex*16 fast(4),slow(4)
nepk=intMeshPot/2+1
nepk2=2*(nepk-1)
if(nepk2/=intMeshPot) error stop 'Odd native energy mesh'
allocate(epk(nepk),zepko(nepk2),zepkn(nepk2),zepk(nepk2))
allocate(zresu(nepk2,3,5,2),rlm(1),slam0(1,2))
allocate(nueff(3,2,2),epsalpha(3,2,2),omegae0(3,2))
allocate(omegasn(3,2,2),omegast(3,2,2),especies_pre(3,2,2))
allocate(depsalphadpsi(3,2,2),cs(3),dpsids(3),dpsidsm(3),t(3),tm(3))
allocate(especies_tem(3,2,2),alphaa1(3,2,2),ispecies_ek(2))
umax=2d0*atan(20d0)/pi
do i=1,nepk
 epk(i)=umax*(i-1)/(nepk-1)
enddo
do i=1,nepk-1
 center=.5d0*(epk(i+1)+epk(i));halfwidth=.5d0*(epk(i+1)-epk(i))
 zepko(2*i-1)=center-halfwidth/sqrt(3d0)
 zepko(2*i)=center+halfwidth/sqrt(3d0)
enddo
zepkn=tan(.5d0*pi*zepko)
intMeshX=zepkn
! Both native kernels integrate the same Maxwellian energy measure.
! The fast denominator contains x**1.5, hence its measure is x**4 exp(-x).
intMeshWeight=umax/(nepk-1)/2d0*.5d0*pi/(cos(.5d0*pi*zepko))**2 &
              *exp(-zepkn)*zepkn**4
zresu=0d0
zresu(:,2,2,1)=-exp(-zepkn)
slam0=-1d0
cs=[0d0,.5d0,1d0]
dpsids=1d0;dpsidsm=1d0;t=1d0;tm=1d0
depsalphadpsi=0d0;especies_tem=1d0;alphaa1=1d0;especies_pre=1d0
especies_z(1:2)=[1d0,-1d0]
especies_m(1:2)=[2d0,.00054463d0]
ispecies_f0(1:2)=0;ispecies_ek=0
pspecies_fowt=0d0
inckin=1;kntv=21;incdphi=0;islsode=.false.
epsalpha=1d0;omegaci0=1d0;nueff=.2d0;omegae0=.37d0
omegasn=0d0;omegast=0d0
rntor=-3d0;b0k=1d0;omegab=.8d0;drift=.3d0;lam=.2d0
call set_test_omega()
largest=0d0
do mode=1,2
 inutype=mode
 do ndb=0,1
  pspecies_ndb=ndb
  do signum=-1,1,2
   rlm=signum
   do grid=1,2
    kfastrun=1
    call ki_bounce(2,grid,1,0,1,fast)
    kfastrun=0
    call ki_bounce(2,grid,1,0,1,slow)
    if(abs(slow(1))<1d-8) error stop 'Vacuous bounce response'
    error=abs(fast(1)-slow(1))/abs(slow(1))
    largest=max(largest,error)
    write(*,'(4I3,1X,ES24.15)') mode,ndb,signum,grid,error
   enddo
  enddo
 enddo
enddo
if(largest>1d-11) error stop 'Fast/slow bounce drift disagreement'
end program
"""


class BounceFastSlow(unittest.TestCase):
    def test_native_trapped_bounce_drift_switch(self):
        kinetic = (NATIVE / "kinetic.f").read_text()
        anisotropic = (NATIVE / "anisotropic.f").read_text()
        fast = re.search(r"^ +SUBROUTINE +KI_BOUNCE\b.*?^ +END *$", kinetic,
                         re.M | re.S | re.I).group()
        slow = re.search(r"^ +SUBROUTINE +KIA_TRAP\b.*?^ +END SUBROUTINE +KIA_TRAP\b",
                         anisotropic, re.M | re.S | re.I).group()
        toolbox = (NATIVE / "ToolBox.f90").read_text().split("    contains", 1)[0]
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
            work = Path(tmp)
            (work / "toolbox_declarations.f90").write_text(toolbox + "end module\n")
            (work / "native.f").write_text(fast + "\n" + slow + "\n" +
                "      SUBROUTINE SET_TEST_OMEGA()\n"
                "      IMPLICIT NONE\n      INCLUDE 'compam.inc'\n"
                "      OMEGA=(0D0,0D0)\n      END\n"
                "      COMPLEX*16 FUNCTION PLASMAY(Z)\n"
                "      COMPLEX*16 Z\n      STOP 'UNEXPECTED ANALYTIC ROUTE'\n      END\n"
                "      SUBROUTINE ADAPTIVE_ENERGY_INTEGRAL_LSODE()\n"
                "      STOP 'UNEXPECTED LSODE ROUTE'\n      END\n")
            (work / "driver.f90").write_text(DRIVER)
            modules = ["globalm.f", "dimensim.f", "rcomdm.f", "kineticm.f", "anisotropicm.f"]
            command = [os.environ.get("FC", "gfortran"), "-O0", "-fcheck=all", "-std=legacy",
                       "-fdefault-real-8", "-fdefault-double-8", "-ffixed-line-length-none",
                       "-ffree-line-length-none", "-fallow-argument-mismatch", "-I"+str(NATIVE),
                       *[str(NATIVE / f) for f in modules], "toolbox_declarations.f90", "native.f",
                       str(NATIVE / "splines.f"), "driver.f90", "-o", "oracle"]
            built = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=90)
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            run = subprocess.run([str(work / "oracle")], cwd=work,
                                 capture_output=True, text=True, timeout=30)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertEqual(len(run.stdout.splitlines()), 16, run.stdout)
            print(run.stdout, end="")
