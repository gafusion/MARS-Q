#!/usr/bin/env python3
"""Native IPDIVB=1 row acting on manufactured divergence-free physical fields.

The m=0 solver unknown is B1/T, as specified by INITXY and GETXY. This test
executes the actual native row block and finite-element integration macros.
"""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

HEADER = '''module fixture
implicit none
integer,parameter:: NR=2,MSMAX=2,NXCOMP=1,NYCOMP=2,KXB1=1,KYB2=1,KYB3=2
integer IPDIVB
real*8 T(3),CSH(2),RM(2,2),PTRAPH
complex*16 FSUBM(4,2,2),GSUBM(4,2,2),DSUBM(4,4,2)
end module
subroutine native_row
use fixture
implicit none
integer I,MSA,LX,LY
complex*16 C1,CMA,CNA
include 'integc.inc'
FSUBM=0;GSUBM=0;DSUBM=0
C1=(1d0,0d0)
DO MSA=1,MSMAX
LX=(MSA-1)*NXCOMP;LY=(MSA-1)*NYCOMP
CMA=cmplx(0d0,RM(MSA,2),kind=8)
CNA=cmplx(0d0,-3d0,kind=8)
'''
FOOTER = '''ENDDO
end subroutine
program regression
use fixture
implicit none
integer mode,profile,field,quad,I,LX,LY,failures,count
real*8 S(3),scale
complex*16 B1(3),B2(2),B3(2),U(3),RES,EXPECTED
IPDIVB=1;CSH=[.2d0,.3d0];S=[.1d0,.3d0,.6d0]
RM(:,1)=0;RM(:,2)=[0d0,2d0];failures=0;count=0
DO quad=1,3
PTRAPH=.9d0
if(quad==1)PTRAPH=0d0
if(quad==3)PTRAPH=1d0
DO profile=1,2
T=1.4d0
if(profile==2)T=[1d0,1.3d0,1.7d0]
call native_row
DO mode=1,2
LX=(mode-1)*NXCOMP;LY=(mode-1)*NYCOMP
DO field=1,2
B1=cmplx(.7d0,.2d0,kind=8)
if(field==2)B1=B1+cmplx(.3d0,-.1d0,kind=8)*S
B2=cmplx(.11d0,.04d0,kind=8)
DO I=1,NR
 B3(I)=-((B1(I+1)-B1(I))/CSH(I)+cmplx(0d0,RM(mode,2),kind=8)*B2(I))/cmplx(0d0,-3d0,kind=8)
ENDDO
U=B1
if(mode==1)U=B1/T
DO I=1,NR
 count=count+1
 RES=FSUBM(LY+KYB3,LX+KXB1,I)*U(I)+GSUBM(LY+KYB3,LX+KXB1,I)*U(I+1) &
     +DSUBM(LY+KYB3,LY+KYB2,I)*B2(I)+DSUBM(LY+KYB3,LY+KYB3,I)*B3(I)
 EXPECTED=(B1(I+1)-B1(I))/CSH(I)+cmplx(0d0,RM(mode,2),kind=8)*B2(I)+cmplx(0d0,-3d0,kind=8)*B3(I)
 scale=max(1d0,abs(B1(I))/CSH(I),abs(B1(I+1))/CSH(I))
 if(abs(EXPECTED)>1d-14*scale)error stop 'Invalid manufactured field'
 if(abs(RES-EXPECTED)>1d-13*scale)then
  failures=failures+1
  write(*,'(A,4I3,2ES24.15)')'FAIL ',quad,profile,mode,field,RES
 endif
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO
write(*,*)'Native divB checks:',count,'failures:',failures
if(failures/=0)error stop 'Physical divergence constraint failed'
end program
'''


def run_check(source, build_parent):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'marsq.f').read_text()
    start = text.index('      IF (IPDIVB.EQ.1) THEN')
    end = text.index('      IF (IPDIVB.EQ.2.AND.', start)
    block = text[start:end]
    with tempfile.TemporaryDirectory(prefix='m0-divb-', dir=build_parent) as directory:
        work = Path(directory)
        # The fixture is free form; only the original native block is fixed form.
        (work / 'prefix.f90').write_text(HEADER[:HEADER.index('subroutine native_row')])
        native = HEADER[HEADER.index('subroutine native_row'):]
        native = native.replace("include 'integc.inc'", '')
        declarations = native[:native.index('FSUBM=0;')]
        setup = native[native.index('FSUBM=0;'):]
        fixed = '\n'.join('      ' + line for line in (declarations + setup).splitlines())
        fixed = fixed.replace('      FSUBM=0;GSUBM=0;DSUBM=0', "      INCLUDE 'integc.inc'\n      FSUBM=0;GSUBM=0;DSUBM=0")
        (work / 'row.f').write_text(fixed + '\n' + block + '      ENDDO\n      END\n')
        (work / 'driver.f90').write_text(FOOTER[FOOTER.index('program regression'):])
        command = ['gfortran', '-O0', '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-fdefault-real-8', '-fdefault-double-8', '-fcheck=bounds',
                   '-ffpe-trap=zero,invalid', '-I' + str(src),
                   'prefix.f90', 'row.f', 'driver.f90', '-o', 'test']
        compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if compiled.returncode:
            raise RuntimeError((compiled.stdout + compiled.stderr)[-8000:])
        result = subprocess.run([str(work / 'test')], cwd=work, capture_output=True, text=True)
        print(result.stdout, end='')
        if result.returncode:
            raise RuntimeError(result.stderr)


class TestM0DivergenceConstraint(unittest.TestCase):
    def test_physical_divergence_free_fields(self):
        run_check(Path(__file__).resolve().parents[1], Path.cwd())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    run_check(args.source, args.build_parent)


if __name__ == '__main__':
    main()
