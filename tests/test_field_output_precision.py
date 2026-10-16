#!/usr/bin/env python3
"""Native OUTPUT field records preserve binary64 coefficients on readback.

Run with unittest discovery, or directly with --source/--build-parent.
The real OUTPUT writes and format are extracted without translating them.
GNU's -fdefault-real-8 -fdefault-double-8 matches MARS's documented -r8 build.
Python's binary64 values and independently parsed records are the oracle;
fixed/free form list-directed consumers also return their readback as binary.
"""
import argparse
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import unittest


SAMPLES = [0.0, -0.0, 1.0 + 2.0**-40, 1.0 - 2.0**-40,
           math.nextafter(1.0, 2.0), math.nextafter(1.0, 0.0),
           1e-8 * (1.0 + 2.0**-40), -1e-8 * (1.0 - 2.0**-40),
           1e-300, -1e300, 2.0**-1022, math.nextafter(0.0, 1.0),
           float.fromhex('0x1.fffffffffffffp+1023'), math.pi, -math.e]

MODULE = '''module fixture
implicit none
integer,parameter::MSMAX=3,NTP1=5,NRP1=4
integer CHOUTP,KSW,ISW,KJRER
real*8 RM(MSMAX,2),RNTOR,DPSIDS(NRP1),T(NRP1)
complex*16 B1U(NTP1,MSMAX),B2U(NTP1,MSMAX),B3U(NTP1,MSMAX), &
 X1U(NTP1,MSMAX),X2U(NTP1,MSMAX),X3U(NTP1,MSMAX), &
 J1U(NTP1,MSMAX),J2U(NTP1,MSMAX),J3U(NTP1,MSMAX), &
 JRE1(NTP1,MSMAX),JRE2(NTP1,MSMAX),JRE3(NTP1,MSMAX),JRE(NTP1,MSMAX)
end module
'''
DRIVER = '''program check
use fixture
implicit none
integer k
open(8,file='fixture.bin',access='stream',form='unformatted',status='old')
read(8) B1U,B2U,B3U,X1U,X2U,X3U,J1U,J2U,J3U,JRE1,JRE2,JRE3,JRE,DPSIDS,T
close(8)
RM(:,1)=0;RM(:,2)=[-2d0,0d0,3d0];RNTOR=-3d0;CHOUTP=20;KSW=-1;ISW=1
call native_b
call native_x
KJRER=0
call native_j
call consumers('JPLASMA.OUT',6,0,'J0')
KJRER=1
call native_j
call consumers('JPLASMA.OUT',8,0,'J1')
KJRER=6
call native_j
call consumers('JPLASMA.OUT',14,0,'J6')
call consumers('BPLASMA.OUT',6,0,'B')
call consumers('XPLASMA.OUT',6,1,'X')
end program
subroutine consumers(file,nc,profiles,tag)
implicit none
character(*),intent(in)::file,tag
integer,intent(in)::nc,profiles
integer kind
character(64) target
! Keep each original text file before the next native J variant overwrites it.
call copy_file(file,trim(tag)//'.txt')
do kind=1,2
 write(target,'(A,I1)')trim(tag)//'.read',kind
 if(kind==1)call read_fixed(file,nc,profiles,trim(target))
 if(kind==2)call read_free(file,nc,profiles,trim(target))
enddo
end subroutine
subroutine copy_file(source,target)
implicit none
character(*),intent(in)::source,target
character(512) line
integer ios
open(10,file=source,status='old');open(11,file=target,status='replace')
do
 read(10,'(A)',iostat=ios)line
 if(ios/=0)exit
 write(11,'(A)')trim(line)
enddo
close(10);close(11)
end subroutine
'''
READER = '''subroutine read_free(file,nc,profiles,target)
implicit none
character(*),intent(in)::file,target
integer,intent(in)::nc,profiles
integer nm,nr,i,j,ios
real*8 fn,v(14)
open(10,file=file,status='old')
open(11,file=target,access='stream',form='unformatted',status='replace')
read(10,*)nm,nr,fn
if(nm/=3.or.fn/=-3d0)error stop 'header'
do j=1,nm
 read(10,*)v(1:nc)
enddo
if(profiles==1)then
 do i=1,nr
  read(10,*)v(1:nc)
  write(11)v(1:nc)
 enddo
endif
do j=1,nm
 do i=1,nr
  read(10,*)v(1:nc)
  write(11)v(1:nc)
 enddo
enddo
read(10,*,iostat=ios)v(1)
if(ios>=0)error stop 'trailing or malformed record'
close(10);close(11)
end subroutine
'''


def native_blocks(source):
    text = (source / 'MarsQ_2FK/marsq.f').read_text()
    markers = [('j', 'C     --- OUTPUT THE PERTURBED PLASMA CURRENT',
                'C     --- OUTPUT THE PERTURBED MAGNETIC FIELD'),
               ('b', 'C     --- OUTPUT THE PERTURBED MAGNETIC FIELD',
                'C     OUTPUT ALL FIELD DATA (B,J,V)'),
               ('x', "         OPEN(CHOUTP,FILE='XPLASMA.OUT')",
                'C     --- OUTPUT THE PERTURBED DENSITY')]
    blocks = []
    formats = '\n'.join(re.findall(r'^ 117[12] FORMAT\(.*$', text, re.M))
    if len(formats.splitlines()) != 2:
        raise ValueError('OUTPUT field formats not found uniquely')
    for name, start_marker, end_marker in markers:
        start = text.index(start_marker)
        if name == 'x':
            start = text.rindex('      IF (KSW.LE.0) THEN', 0, start)
        end = text.index(end_marker, start)
        block = re.sub(r'^ 117[12] FORMAT\(.*\n', '', text[start:end], flags=re.M)
        blocks.append('      SUBROUTINE NATIVE_' + name + '\n'
                      '      USE FIXTURE\n      IMPLICIT NONE\n'
                      '      INTEGER MS,II\n      CHARACTER*256 FILENAME\n'
                      + block + '\n' + formats + '\n      END\n')
    return '\n'.join(blocks)


def bit_bytes(values):
    return struct.pack('=' + 'd' * len(values), *values)


def run_check(source, build_parent):
    if shutil.which('gfortran') is None:
        raise RuntimeError('gfortran is required for the native field-output check')
    fields = [[SAMPLES[(2 * k + 3 * a + ri) % len(SAMPLES)]
               for ri in range(15) for k in range(2)] for a in range(13)]
    profiles = [SAMPLES[(i + 2) % len(SAMPLES)] for i in range(8)]
    expected = {}
    for tag, arrays, nr, nc in [('B', fields[:3], 5, 6),
                               ('X', fields[3:6], 4, 6),
                               ('J0', fields[6:9], 5, 6),
                               ('J1', fields[6:9] + fields[12:13], 5, 8),
                               ('J6', fields[6:13], 5, 14)]:
        values = []
        if tag == 'X':
            for i in range(4):
                values += [profiles[i]] * 3 + [profiles[4 + i]] * 3
        for ms in range(3):
            for i in range(nr):
                for a, array in enumerate(arrays):
                    values += ([0.0, 0.0] if tag in ('J1', 'J6') and i >= 4 and a >= 3
                               else array[2 * (5 * ms + i):2 * (5 * ms + i) + 2])
        expected[tag] = (values, nr, nc)
    with tempfile.TemporaryDirectory(prefix='field-roundtrip-', dir=build_parent) as tmp:
        work = Path(tmp)
        (work / 'fixture.bin').write_bytes(bit_bytes(sum(fields, []) + profiles))
        (work / 'fixture.f90').write_text(MODULE)
        (work / 'output.f').write_text(native_blocks(source))
        (work / 'driver.f90').write_text(DRIVER)
        (work / 'reader.f90').write_text(READER)
        fixed = READER.replace('read_free', 'read_fixed')
        (work / 'reader.f').write_text('\n'.join('      ' + line for line in fixed.splitlines()) + '\n')
        flags = ['-O0', '-fdefault-real-8', '-fdefault-double-8',
                 '-ffixed-line-length-none', '-ffree-line-length-none', '-fcheck=all']
        subprocess.run(['gfortran', *flags, 'fixture.f90', 'output.f', 'reader.f',
                        'reader.f90', 'driver.f90', '-o', 'check'], cwd=work, check=True,
                       capture_output=True, text=True)
        subprocess.run([str(work / 'check')], cwd=work, check=True,
                       capture_output=True, text=True)
        failures = []
        checks = 0
        for tag, (values, nr, nc) in expected.items():
            lines = (work / (tag + '.txt')).read_text().splitlines()
            prefix = 4 + (4 if tag == 'X' else 0)
            if len(lines) != prefix + 3 * nr:
                raise AssertionError((tag, 'record count changed'))
            header = lines[0].split()
            if len(header) != nc or any(len(line.split()) != nc for line in lines[1:]):
                raise AssertionError((tag, 'column count changed'))
            if [int(header[0]), int(header[1]), float(header[2])] != [3, nr, -3.0]:
                raise AssertionError((tag, 'header changed'))
            if [float(line.split()[0]) for line in lines[1:4]] != [-2.0, 0.0, 3.0]:
                raise AssertionError((tag, 'mode order changed'))
            parsed = [float(token) for line in lines[4:] for token in line.split()]
            for name, data in [('python', bit_bytes(parsed)),
                               ('fixed', (work / (tag + '.read1')).read_bytes()),
                               ('free', (work / (tag + '.read2')).read_bytes())]:
                checks += len(values)
                if data != bit_bytes(values):
                    got = struct.unpack('=' + 'd' * (len(data) // 8), data)
                    count = sum(bit_bytes([a]) != bit_bytes([b]) for a, b in zip(got, values))
                    failures.append((tag, name, count))
        print(f'Native field roundtrip: {checks} values, {len(failures)} failed reader/file pairs')
        if failures:
            raise AssertionError(f'Binary64 field roundtrip lost coefficients: {failures}')


class TestFieldOutputPrecision(unittest.TestCase):
    def test_native_field_roundtrip(self):
        run_check(Path(__file__).resolve().parents[1], Path.cwd())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    run_check(args.source, args.build_parent)
