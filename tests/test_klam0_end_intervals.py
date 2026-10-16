#!/usr/bin/env python3
"""Native KLAM0 against known linear/quadratic roots on valid spline intervals.

Requires gfortran. Numerical KLAM0 and SPLINE1DR/SPLINE/RSPLINT are extracted
unchanged from the selected source. Only unused MPI diagnostic output is stubbed.
"""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

MPI = '''module MPIENV
integer RANK
integer,parameter :: CHMPI=99
contains
subroutine MPI_OPEN_FILE(rank_arg)
integer rank_arg
error stop 'Unexpected diagnostic output'
end subroutine
subroutine MPI_CLOSE_FILE(rank_arg)
integer rank_arg
error stop 'Unexpected diagnostic output'
end subroutine
end module
'''

DRIVER = '''program test_klam0
use globalm
use dimensim
use kineticm
use anisotropicm
implicit none
integer branch,shape,side,failures,count
real*8 root,lo,hi
NRP1=2; NSPECIES=1; MLMAX=1; JS0=1
allocate(WFUN(2,2),NLAMK0(2,2),NLAMK1(2,2),LAMM(8),LAMTMP(8))
allocate(ZOMEGADT(2,8,2),ZOMEGABT(2,8,2),ZOMEGABP(2,8,2))
allocate(SLAM0(1,1),RLM(1),ISPECIES_EK(1),EPSALPHA(2,2,1),OMEGAE0(2,2))
NLAMK0=4; NLAMK1=4; WFUN=1; ISPECIES_EK=1
EPSALPHA=.5d0; OMEGAE0=1; RNTOR=1; RQK=1; B0K=1; OMEGACI0=1
ESPECIES_M=1; ESPECIES_Z=1
PSPECIES_NP=1; PSPECIES_NTB=1; PSPECIES_NTD=1
LAMM=[.8d0,.9d0,1d0,1.1d0,1.2d0,1.3d0,1.4d0,1.5d0]
failures=0; count=0
! Both end intervals and an ordinary interior interval in all four branches.
do branch=0,3
 if (branch==1) then
  LAMM=[0d0,.1d0,.2d0,.3d0,.4d0,.5d0,.6d0,.7d0]
  lo=LAMM(1); hi=LAMM(7)
 else
  LAMM=[.8d0,.9d0,1d0,1.1d0,1.2d0,1.3d0,1.4d0,1.5d0]
  lo=LAMM(2); hi=LAMM(8)
 endif
 do shape=1,2
  do side=1,3
   root=lo+.05d0
   if(side==2) root=(lo+hi)/2d0+.025d0
   if(side==3) root=hi-.05d0
   call probe(branch,shape,root,root)
  enddo
 enddo
 ! No crossing and an unsampled tangent remain absent.
 call probe(branch,3,(lo+hi)/2d0+.025d0,-1d0)
 call probe(branch,4,(lo+hi)/2d0+.025d0,-1d0)
 ! Multiple crossings retain the native first-bracket policy.
 call probe(branch,5,lo+.05d0,lo+.05d0)
 ! A crossing outside the branch's spline domain must not activate a root.
 if(branch==1) then
  call probe(branch,1,.65d0,-1d0)
 else
  call probe(branch,1,.85d0,-1d0)
 endif
 ! True physical endpoints cannot support the interior singularity formulas.
 call probe(branch,1,LAMM(1),-1d0)
 call probe(branch,1,LAMM(8),-1d0)
 WFUN=0
 call probe(branch,1,lo+.05d0,-1d0)
 WFUN=1
enddo
write(*,*) 'KLAM0 checks:',count,'failures:',failures
if(failures/=0) error stop 'Known-root oracle failed'
contains
subroutine probe(branch,shape,root,expected)
integer,intent(in)::branch,shape
real*8,intent(in)::root,expected
integer j,particle
real*8 target,d,actual
COMPLEX*16 AL0,ALNORM,ALAM,ATAU,ALTAU,OMEGA
REAL*8 EPSPAM,EPSDET
INTEGER NCASE,NITMAX,NONCON
COMMON /PAMARG/ AL0,ALNORM,ALAM,ATAU,ALTAU,OMEGA,EPSPAM,EPSDET, &
 NCASE,NITMAX,NONCON
OMEGA=(0d0,0d0)
RLM=0; target=0; particle=2
if(branch==1) then
 RLM=1; target=.5d0; particle=1
elseif(branch==2) then
 RLM=-1; target=1d0; particle=0
elseif(branch==3) then
 RLM=0; target=-2d0; particle=0
endif
count=count+1
! The polynomial oracle specifies its roots independently of the native spline.
do j=1,8
 d=.2d0*(LAMM(j)-root)
 if(shape==2) d=.2d0*(LAMM(j)-root)*(1d0+.25d0*(LAMM(j)-root))
 if(shape==3) d=.2d0+(LAMM(j)-root)**2
 if(shape==4) d=(LAMM(j)-root)**2
 if(shape==5) d=.2d0*(LAMM(j)-root)*(LAMM(j)-root-.275d0)*(LAMM(j)-root-.425d0)
 ZOMEGABP(:,j,:)=target+d
 ZOMEGABT(:,j,:)=target+d
 ZOMEGADT(:,j,:)=target+d
enddo
call KLAM0(2,1,particle)
actual=SLAMD0
if(particle/=2) actual=SLAM0(1,1)
if(abs(actual-expected)>1d-8) then
 failures=failures+1
 write(*,'(A,3I3,3ES24.15)') 'FAIL ',branch,shape,count,root,expected,actual
endif
end subroutine
end program
'''


def run_check(source, build_parent):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'kinetic.f').read_text()
    start = re.search(r'^      SUBROUTINE KLAM0\(', text, re.M).start()
    end = start + re.search(r'^      END\s*$', text[start:], re.M).end()
    with tempfile.TemporaryDirectory(prefix='klam0-', dir=build_parent) as directory:
        work = Path(directory)
        (work / 'klam0.f').write_text(text[start:end] + '\n')
        (work / 'mpi.f90').write_text(MPI)
        (work / 'driver.f90').write_text(DRIVER)
        names = ['globalm.f', 'dimensim.f', 'rcomdm.f', 'kineticm.f',
                 'anisotropicm.f', 'ToolBox.f90', 'splines.f']
        command = ['gfortran', '-O0', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-fallow-argument-mismatch', '-ffunction-sections', '-fdata-sections',
                   '-ffpe-trap=zero,invalid', '-fcheck=bounds', '-I' + str(src)]
        command += [str(src / name) for name in names]
        command += ['mpi.f90', 'klam0.f', 'driver.f90', '-Wl,--gc-sections', '-o', 'test']
        compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if compiled.returncode:
            raise RuntimeError(compiled.stdout + compiled.stderr)
        result = subprocess.run([str(work / 'test')], cwd=work, capture_output=True, text=True)
        print(result.stdout, end='')
        if result.returncode:
            raise RuntimeError(result.stderr)


class TestKLAM0EndIntervals(unittest.TestCase):
    def test_known_roots_and_guards(self):
        run_check(Path(__file__).resolve().parents[1], Path.cwd())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    run_check(args.source, args.build_parent)


if __name__ == '__main__':
    main()
