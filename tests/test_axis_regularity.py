#!/usr/bin/env python3
"""Magnetic-axis rows of BOUNDC and the AXISY reconstruction.

The production decks BOUNDC, AXISY, IAXPWX, IAXPWY, ANNIHX and ANNIHY are compiled
unchanged against small fixture modules. The oracle tables below are written from the
derivation in the note, independently of IAXPWX/IAXPWY.
"""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

import numpy as np

DECKS = ('BOUNDC', 'AXISY', 'IAXPWX', 'IAXPWY', 'ANNIHX', 'ANNIHY')

MODULES = '''module dimensim
implicit none
integer NR,NRP1,NFIT,MSMAX,MXMAX,MYMAX
end module
module globalm
implicit none
integer NXCOMP,NYCOMP,K_BC_N0
integer KXV1,KXB1,KXJ2U,KXJ3,KXJ2L,KXX1,KXPD,KXPED,KXW1,KXJRE,KXB2L,KXB3L, &
        KXJRE2,KXJRE3,KXJRE2L,KXDPHI
integer KYV2,KYV3,KYB2,KYB3,KYJ1,KYPR,KYPPERP,KYPPARA,KYPE,KYPP,KYRHOP, &
        KYX2,KYX3,KYW2,KYW3,KYJRE1
real*8 RNTOR,SLEFT
real*8,allocatable :: RM(:,:),CS(:),CSM(:),CSH(:)
end module
'''

# Leading power p of s at the axis, Y ~ s**p, from the derivation (note, table 2).
# Entries are (p for m=0, p(|m|) for |m|>=1). On a non-circular axis |m| >= 3 takes
# the power of meff = 1 (odd) or 2 (even); B3 and J3 also carry m = 0 through J/s.
def meff(a):
    return a if a < 3 else 2 - a % 2


PY = {'V2': (0, lambda a: meff(a) - 2), 'B2': (1, lambda a: meff(a) - 1),
      'B3': (1, lambda a: 2 if meff(a) == 1 else 1), 'J1': (2, meff), 'PR': (0, meff),
      'V3': (0, meff), 'X2': (0, lambda a: meff(a) - 2), 'X3': (0, meff)}
PX = {'V1': (1, lambda a: meff(a) - 1), 'B1': (2, meff), 'J2U': (1, lambda a: meff(a) - 1),
      'J3': (1, lambda a: 2 if meff(a) == 1 else 1), 'J2L': (2, meff),
      'X1': (1, lambda a: meff(a) - 1), 'PD': (0, meff)}
XNAMES = ['V1', 'B1', 'J2U', 'J3', 'J2L', 'X1', 'PD']
YNAMES = ['V2', 'B2', 'B3', 'J1', 'PR', 'V3', 'X2', 'X3']


def power(table, name, m):
    p0, p = table[name]
    return p0 if m == 0 else p(abs(m))


def fortran_table(name, table, names):
    rows = []
    for name_k in names:
        rows.append(','.join(str(power(table, name_k, m)) for m in range(-3, 4)))
    data = ','.join(rows)
    return f'integer,parameter :: {name}({len(names)},7) = reshape([{data}],[{len(names)},7],order=[2,1])'


DRIVER = '''program axis_rows
use dimensim
use globalm
implicit none
complex*16,allocatable :: AM(:,:,:),BM(:,:,:),CM(:,:,:),DM(:,:,:),EM(:,:,:),FM(:,:,:),GM(:,:,:),HM(:,:,:)
complex*16,allocatable :: AM0(:,:,:),BM0(:,:,:),CM0(:,:,:),DM0(:,:,:),EM0(:,:,:),FM0(:,:,:),GM0(:,:,:),HM0(:,:,:)
complex*16,allocatable :: Y(:,:,:)
complex*16 :: expect,b1n(2),b2,b3,cim
real*8 :: hh0,a,s1,s2,err(2),exact
integer :: ms,k,kk,ly,lx,i,failures,ilev,ip,mabs,mode
@PXTAB@
@PYTAB@
failures=0
NR=8;NRP1=NR+1;NFIT=3;MSMAX=7;NXCOMP=7;NYCOMP=8;K_BC_N0=11
MXMAX=MSMAX*NXCOMP;MYMAX=MSMAX*NYCOMP
KXV1=1;KXB1=2;KXJ2U=3;KXJ3=4;KXJ2L=5;KXX1=6;KXPD=7
KXPED=-1;KXW1=-1;KXJRE=-1;KXB2L=-1;KXB3L=-1;KXJRE2=-1;KXJRE3=-1;KXJRE2L=-1;KXDPHI=-1
KYV2=1;KYB2=2;KYB3=3;KYJ1=4;KYPR=5;KYV3=6;KYX2=7;KYX3=8
KYPE=-1;KYPP=-1;KYPPARA=-1;KYPPERP=-1;KYRHOP=-1;KYW2=-1;KYW3=-1;KYJRE1=-1
RNTOR=-3d0
allocate(RM(MSMAX,2),CS(NRP1+1),CSM(NRP1),CSH(NRP1))
do ms=1,MSMAX
 RM(ms,:)=dble(ms-4)
enddo
allocate(AM(MXMAX,MXMAX,NRP1+1),BM(MXMAX,MXMAX,NRP1+1),CM(MXMAX,MXMAX,NRP1+1),DM(MYMAX,MYMAX,NRP1+1), &
 EM(MXMAX,MYMAX,NRP1+1),FM(MYMAX,MXMAX,NRP1+1),GM(MYMAX,MXMAX,NRP1+1),HM(MXMAX,MYMAX,NRP1+1))
allocate(Y(NYCOMP,MSMAX,NR))

do mode=1,2
 call grid(0.05d0)
 SLEFT=0d0
 if(mode==2)SLEFT=0.1d0
 call fill()
 call BOUNDC(AM,BM,CM,DM,EM,FM,GM,HM)
 do ms=1,MSMAX
  lx=(ms-1)*NXCOMP;ly=(ms-1)*NYCOMP;mabs=abs(nint(RM(ms,2)))
  ! node-1 X rows: identity, plus X(1)=X(2) where the component is finite
  do k=1,NXCOMP
   do kk=1,MXMAX
    expect=0;if(kk==k+lx)expect=1
    call chk(BM(k+lx,kk,1),expect,'X row B')
    expect=0
    if(kk==k+lx.and.mode==1.and.pxt(k,ms)==0)expect=-1
    if(kk==k+lx.and.mode==2.and.k==KXV1.and.mabs==1)expect=-1
    call chk(CM(k+lx,kk,1),expect,'X row C')
    call chk(AM(k+lx,kk,1),(0d0,0d0),'X row A')
   enddo
   do kk=1,MYMAX
    call chk(EM(k+lx,kk,1),(0d0,0d0),'X row E')
    call chk(HM(k+lx,kk,1),(0d0,0d0),'X row H')
   enddo
  enddo
  ! half-cell-1 Y rows
  do k=1,NYCOMP
   if(mode==2)then
    ! inner boundary: unchanged earlier rule, only the diagonal is set
    do kk=1,MYMAX
     expect=DM0(k+ly,kk,1);if(kk==k+ly)expect=1
     call chk(DM(k+ly,kk,1),expect,'inner D')
    enddo
    cycle
   endif
   do kk=1,MYMAX
    expect=0
    if(kk==k+ly)expect=1
    if(k==KYB2.and.mabs>0)then
     expect=0
     if(kk==KYB2+ly)expect=dcmplx(0d0,RM(ms,2))
     if(kk==KYB3+ly)expect=dcmplx(0d0,RNTOR)
    endif
    call chk(DM(k+ly,kk,1),expect,'Y row D')
   enddo
   do kk=1,MXMAX
    expect=0;if(k==KYB2.and.mabs>0.and.kk==KXB1+lx)expect=-1/CSH(1)
    call chk(FM(k+ly,kk,1),expect,'Y row F')
    expect=0;if(k==KYB2.and.mabs>0.and.kk==KXB1+lx)expect=1/CSH(1)
    call chk(GM(k+ly,kk,1),expect,'Y row G')
   enddo
  enddo
 enddo
 ! V1 fit rows 2..NFIT: V1 ~ s**p with the axis power (inner boundary: |m|-1)
 do i=2,NFIT
  do ms=1,MSMAX
   lx=(ms-1)*NXCOMP;mabs=abs(nint(RM(ms,2)))
   ip=pxt(KXV1,ms)
   if(mode==2)ip=mabs-1
   if(mabs==0.and.mode==2)ip=1
   call chk(BM(KXV1+lx,KXV1+lx,i),dcmplx(CS(i)**(-ip),0d0),'V1 fit B')
   call chk(CM(KXV1+lx,KXV1+lx,i),dcmplx(-CS(i+1)**(-ip),0d0),'V1 fit C')
  enddo
 enddo
 ! interior cells are untouched apart from the V1 fit rows 2..NFIT
 do i=2,NR
  do ms=1,MSMAX
   ly=(ms-1)*NYCOMP
   do k=1,NYCOMP
    do kk=1,MYMAX
     call chk(DM(k+ly,kk,i),DM0(k+ly,kk,i),'interior D')
    enddo
   enddo
  enddo
 enddo
enddo
print *,'axis-row checks done, failures:',failures

! Manufactured divergence-free field, b ~ curl of an analytic potential:
! B1 = i m s**|m| (1 + a s**2), B3 = 0.7 s**(|m|+1), B2 from div b = 0.
! The B2 row of the axis cell must reproduce B2 at CSM(1) to O(hh0) (|m|=1: O(hh0**2)).
SLEFT=0d0
a=0.4d0
do ms=1,MSMAX
 mabs=abs(nint(RM(ms,2)))
 if(mabs==0)cycle
 do ilev=1,2
  hh0=0.04d0/2**(ilev-1)
  call grid(hh0)
  call fill()
  call BOUNDC(AM,BM,CM,DM,EM,FM,GM,HM)
  lx=(ms-1)*NXCOMP;ly=(ms-1)*NYCOMP
  b1n(1)=b1(CS(1));b1n(2)=b1(CS(2))
  b3=0
  b2=-(FM(KYB2+ly,KXB1+lx,1)*b1n(1)+GM(KYB2+ly,KXB1+lx,1)*b1n(2)+DM(KYB2+ly,KYB3+ly,1)*b3) &
     /DM(KYB2+ly,KYB2+ly,1)
  err(ilev)=abs(b2-b2x(CSM(1)))/abs(b2x(0.5d0))
 enddo
 print '(A,I3,2ES12.3)',' div-b axis row m, rel err h, h/2',nint(RM(ms,2)),err
 if(err(2)>0.6d0*err(1).or.err(2)>2d-2)then
  failures=failures+1;print *,'FAIL div-b convergence m=',nint(RM(ms,2))
 endif
enddo

! AXISY: every component of half-cell 1 (closed by Y(1)=0 in the solve) is rebuilt
! from half-cell 2 with s**p, since Y(1) is the value at CSM(1). B2 (m/=0) keeps its
! div-b row and takes the rebuilt B3: B2 -> B2 - n/m (B3(1) - 0).
do ms=1,MSMAX
 mabs=abs(nint(RM(ms,2)))
 do k=1,NYCOMP
  ip=pyt(k,ms)
  do ilev=1,2
   hh0=0.04d0/2**(ilev-1)
   call grid(hh0)
   Y=(7d0,7d0)
   Y(:,ms,1)=(0d0,0d0)
   Y(k,ms,2)=fpow(CSM(2),ip)
   if(k==KYB2.and.mabs>0)then
    Y(KYB2,ms,1)=(0.3d0,0.2d0)
    call AXISY(Y)
    exact=0
    if(abs(Y(KYB2,ms,1)-((0.3d0,0.2d0)-RNTOR/RM(ms,2)*Y(KYB3,ms,1)))>1d-12)then
     failures=failures+1;print *,'FAIL AXISY B2 div',nint(RM(ms,2))
    endif
    err(ilev)=0
   else
    call AXISY(Y)
    err(ilev)=abs(Y(k,ms,1)-fpow(CSM(1),ip))/abs(fpow(CSM(1),ip))
   endif
  enddo
  if(err(1)>0)then
   if(err(2)>0.3d0*err(1).or.err(2)>1d-2)then
    failures=failures+1;print *,'FAIL AXISY convergence',k,nint(RM(ms,2)),err
   endif
  endif
 enddo
enddo
print *,'total failures:',failures
if(failures/=0)error stop 1
contains
subroutine grid(hg)
real*8,intent(in)::hg
integer::j
do j=1,NRP1+1
 CS(j)=(j-1)*hg
enddo
CS(1)=1d-3*CS(2)
do j=1,NRP1
 CSM(j)=0.5d0*(CS(j)+CS(j+1));CSH(j)=CS(j+1)-CS(j)
enddo
end subroutine
subroutine fill()
call rnd(AM);call rnd(BM);call rnd(CM);call rnd(DM);call rnd(EM);call rnd(FM);call rnd(GM);call rnd(HM)
AM0=AM;BM0=BM;CM0=CM;DM0=DM;EM0=EM;FM0=FM;GM0=GM;HM0=HM
end subroutine
subroutine rnd(z)
complex*16,intent(inout)::z(:,:,:)
real*8,allocatable::x(:,:,:),yy(:,:,:)
allocate(x(size(z,1),size(z,2),size(z,3)),yy(size(z,1),size(z,2),size(z,3)))
call random_number(x);call random_number(yy)
z=dcmplx(x+0.5d0,yy-0.5d0)
end subroutine
subroutine chk(got,want,what)
complex*16,intent(in)::got,want
character(*),intent(in)::what
if(abs(got-want)>1d-12*max(1d0,abs(want)))then
 failures=failures+1
 if(failures<20)print *,'FAIL ',what,mode,ms,k,kk,got,want
endif
end subroutine
complex*16 function b1(s)
real*8,intent(in)::s
b1=dcmplx(0d0,RM(ms,2))*s**mabs*(1+a*s**2)
end function
complex*16 function b2x(s)
real*8,intent(in)::s
complex*16::db1
db1=dcmplx(0d0,RM(ms,2))*(mabs*s**(mabs-1)*(1+a*s**2)+2*a*s**(mabs+1))
b2x=-(db1+dcmplx(0d0,RNTOR)*0.7d0*s**(mabs+1))/dcmplx(0d0,RM(ms,2))
end function
complex*16 function fpow(s,p)
real*8,intent(in)::s
integer,intent(in)::p
fpow=dcmplx(1d0,-0.5d0)*s**p*(1+0.9d0*s**2)
end function
end program
'''


def deck(text, name):
    start = re.search(r'^\*DECK ' + name + r'\b', text, re.M)
    if start is None:
        raise ValueError('deck not found: ' + name)
    end = text.find('\n*DECK', start.end())
    return text[start.end():end if end >= 0 else len(text)]


def run_check(source, build_parent):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'marsq.f').read_text()
    native = ''.join(deck(text, name) for name in DECKS)
    driver = (DRIVER.replace('@PXTAB@', fortran_table('pxt', PX, XNAMES))
              .replace('@PYTAB@', fortran_table('pyt', PY, YNAMES)))
    with tempfile.TemporaryDirectory(prefix='axis-rows-', dir=build_parent) as directory:
        work = Path(directory)
        (work / 'modules.f90').write_text(MODULES)
        (work / 'native.f').write_text(native)
        (work / 'driver.f90').write_text(driver)
        command = ['gfortran', '-O0', '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-fdefault-real-8', '-fdefault-double-8', '-fcheck=bounds',
                   '-ffpe-trap=zero,invalid', '-I' + str(src),
                   'modules.f90', 'native.f', 'driver.f90', '-o', 'test']
        compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if compiled.returncode:
            raise RuntimeError((compiled.stdout + compiled.stderr)[-8000:])
        result = subprocess.run([str(work / 'test')], cwd=work, capture_output=True, text=True)
        print(result.stdout, end='')
        if result.returncode:
            raise RuntimeError(result.stdout[-4000:] + result.stderr)


def shaped_powers(mmax=3, delta=0.1, kappa=1.5):
    """Axis powers of chi-harmonics of smooth fields on an elliptic axis whose angle chi
    is not the elliptic angle u = chi + delta sin(2 chi), as for equal-arc chi.
    Independent of the tables: smooth Cartesian fields, exact map, FFT in chi."""
    chi = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    u, du = chi + delta * np.sin(2 * chi), 1 + 2 * delta * np.cos(2 * chi)

    def fields(s):
        x, y = s * np.cos(u), kappa * s * np.sin(u)
        xs, ys = np.cos(u), kappa * np.sin(u)
        xc, yc = -s * np.sin(u) * du, kappa * s * np.cos(u) * du
        d = xs * yc - xc * ys
        rr = 3 + x
        xi = (0.7 + 0.2 * x, -0.4 + 0.1 * y)
        return {'X1': (xi[0] * yc - xi[1] * xc) / d,                 # xi . grad s
                'B3': rr * d * (1 + 0.3 * x + 0.2 * y) / rr,         # J b.grad(phi)
                'PR': 1 + 0.3 * x + 0.2 * y + 0.5 * x * y + 0.1 * x * x}

    s1, s2 = 1e-4, 2e-4
    f1, f2 = fields(s1), fields(s2)
    out = {}
    for name in f1:
        c1, c2 = np.fft.fft(f1[name]) / chi.size, np.fft.fft(f2[name]) / chi.size
        a1 = [max(abs(c1[m]), abs(c1[-m])) for m in range(mmax + 1)]
        a2 = [max(abs(c2[m]), abs(c2[-m])) for m in range(mmax + 1)]
        # a harmonic absent from the field has no power
        out[name] = [int(round(np.log2(y / x))) if x > 1e-12 * max(a1) else None
                     for x, y in zip(a1, a2)]
    return out


class TestAxisRegularity(unittest.TestCase):
    def test_axis_rows_and_reconstruction(self):
        run_check(Path(__file__).resolve().parents[1], Path.cwd())

    def test_shaped_axis_powers(self):
        got = shaped_powers()
        for name, table in (('X1', PX), ('B3', PY), ('PR', PY)):
            self.assertEqual(got[name], [power(table, name, m) for m in range(4)], name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    run_check(args.source, args.build_parent)


if __name__ == '__main__':
    main()
