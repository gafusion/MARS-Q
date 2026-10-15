#!/usr/bin/env python3
"""Exercise the actual KI0 routine with affine responses of known integral.

Requires gfortran. This isolated check uses the native modules and splines;
only the energy response and pitch-independent distribution are manufactured.
Build products go to a temporary directory under the selected build parent.
"""
from pathlib import Path
import argparse
import math
import re
import subprocess
import tempfile
import unittest

DRIVER = '''program test_ki0_mesh
use globalm
use dimensim
use rcomdm
use kineticm
use anisotropicm
implicit none
COMPLEX*16 AL0,ALNORM,ALAM,ATAU,ALTAU,OMEGA
REAL*8 EPSPAM,EPSDET
INTEGER NCASE,NITMAX,NONCON
COMMON /PAMARG/ AL0,ALNORM,ALAM,ATAU,ALTAU,OMEGA,EPSPAM,EPSDET, &
 NCASE,NITMAX,NONCON
integer j,k
real*8 fraction(4),expected_re,expected_im
NR=3; NRP1=4; NSPECIES=1; NEPK=3; NEPK2=5
WK=1d0/sqrt(3d0); OMEGA=(.1d0,.1d0)
RNTOR=3; B0K=1; OMEGACI0=100; DPRM=1
PSPECIES_NTD=1
allocate(NLAMK0(4,2),LAMK0(4,4,2),LAMM(8),LAMTMP(8))
allocate(ZOMEGADT(4,8,2),WFUN(4,2),NUEFF(4,2,1),OMEGAE0(4,2))
NLAMK0=4; WFUN=1; NUEFF=.01d0; OMEGAE0=.02d0
LAMK0(:,1,:)=.8d0; LAMK0(:,2,:)=.9d0
LAMK0(:,3,:)=1.1d0; LAMK0(:,4,:)=1.2d0
fraction=[.5d0,.00625d0,1d-6,1d0-1d-6]
do k=1,4
 SLAMD0=.8d0+.4d0*fraction(k)
 do j=1,8
  LAMM(j)=.8d0+.4d0*(j-1)/7d0
  ZOMEGADT(:,j,:)=LAMM(j)-SLAMD0
 enddo
 call KI0(2,1,0)
 expected_re=.4d0/sqrt(PI)
 expected_im=.8d0/sqrt(PI)
 if (abs(real(VI0(1,1))-expected_re)>1d-12) error stop 'Affine real integral'
 if (abs(aimag(VI0(1,1))-expected_im)>1d-12) error stop 'Affine imaginary integral'
 if (abs(VI0(2,1)-cmplx(.4d0,1.2d0,kind=8)/sqrt(PI))>1d-12) &
  error stop 'Linear integral'
 write(*,'(I2,2ES26.17E3)') k,VI0(1,1)
enddo
end program
subroutine KDISTRIBF_TYPE(js,grid,kp,lambda)
implicit none
integer js,grid,kp
real*8 lambda
end subroutine
subroutine KI_PRECESSION(js,grid,particle,kp,kopt,v)
use kineticm, only:LAM
implicit none
integer js,grid,particle,kp,kopt
complex*16 v(4)
v=cmplx(1d0,2d0,kind=8)
v(2)=LAM*cmplx(1d0,3d0,kind=8)
end subroutine
'''


def run_check(source, build_parent, fpe_trap=True):
    src=source.resolve()/'MarsQ_2FK'
    text=(src/'kinetic.f').read_text()
    start=re.search(r'^      SUBROUTINE KI0\(',text,re.M).start()
    end=re.search(r'^      END\s*$',text[start:],re.M).end()+start
    with tempfile.TemporaryDirectory(prefix='ki0-mesh-',dir=build_parent) as directory:
        work=Path(directory)
        (work/'ki0.f').write_text(text[start:end]+'\n')
        (work/'driver.f90').write_text(DRIVER)
        files=[src/name for name in ['globalm.f','dimensim.f','rcomdm.f',
               'kineticm.f','anisotropicm.f','ToolBox.f90','splines.f']]
        command=['gfortran','-O0','-fdefault-real-8','-fdefault-double-8',
                 '-ffixed-line-length-none','-ffree-line-length-none',
                 '-fallow-argument-mismatch','-ffunction-sections','-fdata-sections',
                 '-I'+str(src)]
        if fpe_trap:
            command.append('-ffpe-trap=zero,invalid')
        command+=list(map(str,files))+['ki0.f','driver.f90','-Wl,--gc-sections','-o','test']
        compiled=subprocess.run(command,cwd=work,capture_output=True,text=True)
        if compiled.returncode:
            raise RuntimeError(compiled.stdout+compiled.stderr)
        result=subprocess.run([str(work/'test')],cwd=work,capture_output=True,text=True)
        print(result.stdout,end='')
        if result.returncode:
            raise RuntimeError(result.stderr)
        values=[[float(x) for x in line.split()[1:]] for line in result.stdout.splitlines()]
        assert len(values)==4 and all(math.isfinite(x) for v in values for x in v)
        print('KI0 constant/linear integral PASS (four valid reversal locations)')

class TestKI0PitchMesh(unittest.TestCase):
    def test_constant_and_linear_responses(self):
        run_check(Path(__file__).resolve().parents[1], Path.cwd())


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--build-parent',type=Path,default=Path.cwd())
    ap.add_argument('--no-fpe-trap',action='store_true')
    args=ap.parse_args()
    run_check(args.source,args.build_parent,not args.no_fpe_trap)


if __name__=='__main__':
    main()
