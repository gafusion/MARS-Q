#!/usr/bin/env python3
"""Actual native OUTPUT records preserve plasma-array rows and Fourier modes.

The CLEAR allocations, nineteen OUTPUT calls, GENOUT and GENINP are compiled
verbatim. Independent integer-coded complex coefficients are checked in Python
and after native input roundtrip. NV=0 is a control; NV>0 separates the plasma
and plasma-plus-vacuum leading dimensions. No equilibrium or solver run needed.
"""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


FIELDS = [('V1U', 'plasma_full'), ('V2U', 'plasma_half'), ('V3U', 'plasma_half'),
          ('B1U', 'total_full'), ('B2U', 'total_half'), ('B3U', 'total_half'),
          ('J1U', 'total_half'), ('J2U', 'total_full'), ('J3U', 'total_full'),
          ('PRE', 'plasma_half'), ('PEE', 'plasma_half'), ('PEP', 'plasma_half'),
          ('PDE', 'plasma_full'), ('PED', 'plasma_full'), ('PPERP', 'plasma_half'),
          ('PPARA', 'plasma_half'), ('X1U', 'plasma_full'),
          ('X2U', 'plasma_half'), ('X3U', 'plasma_half')]


def routine(text, name):
    match = re.search(r'^\s*SUBROUTINE\s+' + name + r'\b', text, re.I | re.M)
    if match is None:
        raise AssertionError('Missing native routine: ' + name)
    end = re.search(r'^\s*END\s*$', text[match.start():], re.I | re.M)
    return text[match.start():match.start() + end.end()] + '\n'


def run_check(source, build_parent):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'marsq.f').read_text()
    clear = routine(text, 'CLEAR')
    start = clear.index('      IF (.NOT. ALLOCATED(V1U)) THEN')
    allocate = clear[start:clear.index('      ENDIF', start) + len('      ENDIF')]
    output = routine(text, 'OUTPUT')
    start = output.index('      CALL GENOUT(V1U,')
    calls = output[start:output.index('\nC\n', start)]
    native = ('      SUBROUTINE NATIVE_ALLOCATE\n'
              '      USE DIMENSIM\n      USE GLOBALM\n' + allocate + '\n'
              '      END\n      SUBROUTINE NATIVE_OUTPUT(CHOUTP)\n'
              '      USE DIMENSIM\n      USE GLOBALM\n'
              '      INTEGER CHOUTP\n' + calls + '\n      END\n' +
              routine(text, 'GENOUT') + routine(text, 'GENINP'))
    counts = {'plasma_full': 'NRP1', 'plasma_half': 'NR',
              'total_full': 'NTP1', 'total_half': 'NTP1-1'}
    fill = '\n'.join(f'call fill({name},{tag})'
                     for tag, (name, _) in enumerate(FIELDS, 1))
    read = '\n'.join(f"call geninp({name},size({name},1),'{name:>6}',20,"
                     f"{counts[kind]},MSMAX,2)\n"
                     f"call compare({name},{counts[kind]},{tag})"
                     for tag, (name, kind) in enumerate(FIELDS, 1))
    zero = '\n'.join(name + '=(-999d0,999d0)' for name, _ in FIELDS)
    driver = '''program check
use dimensim
use globalm
implicit none
integer m
character(20) arg
call get_command_argument(1,arg)
read(arg,*)NV
NRP1=5;NR=NRP1-1;NTP1=NRP1+NV;MSMAX=3;MSDIM=8;MEDIM=MSMAX
allocate(RM(MSMAX,2))
RM(:,1)=0d0;RM(:,2)=[-2d0,0d0,3d0];RN=0d0;RN(2)=-3d0
call native_allocate
''' + fill + '''
open(20,file='native.out',status='replace')
call native_output(20)
close(20)
''' + zero + '''
open(20,file='native.out',status='old')
''' + read + '''
close(20)
write(*,*)'PASS native GENINP roundtrip, NV=',NV
contains
subroutine fill(a,tag)
complex*16,intent(out)::a(:,:)
integer,intent(in)::tag
integer i,m,k
do m=1,size(a,2)
 do i=1,size(a,1)
  k=10000*tag+100*m+i
  a(i,m)=cmplx(dble(k),dble(-2*k),kind=8)
 enddo
enddo
end subroutine
subroutine compare(a,np,tag)
complex*16,intent(in)::a(:,:)
integer,intent(in)::np,tag
integer i,m,k
! The native input must not overwrite rows or modes absent from the file.
do m=1,size(a,2)
 do i=1,size(a,1)
  if(m<=MSMAX.and.i<=np)then
   k=10000*tag+100*m+i
   if(a(i,m)/=cmplx(dble(k),dble(-2*k),kind=8))then
    write(*,*)'ROUNDTRIP MISMATCH',tag,i,m,a(i,m),k
    return
   endif
  else
   if(a(i,m)/=cmplx(-999d0,999d0,kind=8))error stop 'GENINP padding changed'
  endif
 enddo
enddo
end subroutine
end program
'''
    checked = failures = 0
    by_nv = {}
    with tempfile.TemporaryDirectory(prefix='genout-stride-', dir=build_parent) as directory:
        work = Path(directory)
        (work / 'native.f').write_text(native)
        (work / 'driver.f90').write_text(driver)
        command = ['gfortran', '-O0', '-g', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffixed-line-length-none', '-fcheck=all',
                   '-ffpe-trap=invalid,zero,overflow', str(src / 'dimensim.f'),
                   str(src / 'globalm.f'), 'native.f', 'driver.f90', '-o', 'check']
        built = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if built.returncode:
            raise RuntimeError((built.stdout + built.stderr)[-5000:])
        for nv in (0, 1, 3):
            result = subprocess.run([str(work / 'check'), str(nv)], cwd=work,
                                    capture_output=True, text=True)
            if result.returncode:
                raise AssertionError(result.stdout + result.stderr)
            before = failures
            rows = [line for line in (work / 'native.out').read_text().splitlines()
                    if line.strip()]
            offset = 0
            for tag, (name, kind) in enumerate(FIELDS, 1):
                np = {'plasma_full': 5, 'plasma_half': 4,
                      'total_full': 5 + nv, 'total_half': 4 + nv}[kind]
                for mode_index, mode in enumerate((-2, 0, 3), 1):
                    header = rows[offset].split()
                    assert header == [name, f'{mode}.', '-3.'], header
                    offset += 1
                    for radial_index in range(1, np + 1):
                        line = rows[offset]
                        assert len(line) == 60, line
                        actual = complex(float(line[:30].replace('D', 'E')),
                                         float(line[30:].replace('D', 'E')))
                        k = 10000 * tag + 100 * mode_index + radial_index
                        checked += 1
                        failures += actual != complex(k, -2*k)
                        offset += 1
            assert offset == len(rows)
            assert ('ROUNDTRIP MISMATCH' not in result.stdout) == (failures == before)
            by_nv[nv] = failures - before
    print(f'Native GENOUT checks: {checked} failures: {failures}; by_NV={by_nv}')
    if failures:
        raise AssertionError('Native OUTPUT row/mode corruption: ' + str(by_nv))


class TestGenoutPlasmaStride(unittest.TestCase):
    def test_native_output_and_input_preserve_plasma_modes(self):
        run_check(Path(__file__).resolve().parents[1], Path.cwd())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    run_check(args.source, args.build_parent)
