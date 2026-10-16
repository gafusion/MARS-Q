#!/usr/bin/env python3
"""FMIND must invert a well-conditioned matrix independently of scratch data."""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

DRIVER = '''program regression
implicit none
complex*16 a(2,2),scratch(2,5),pivots(2),det,expected(2,2)
character(16) mode
call get_command_argument(1,mode)
a=0d0;a(1,1)=.125d0;a(2,2)=.25d0
expected=0d0;expected(1,1)=8d0;expected(2,2)=4d0
scratch=0d0
if(trim(mode)=='finite')scratch=cmplx(huge(1d0)/4d0,0d0,kind=8)
call fmind(a,pivots,2,2,scratch,det,1d-199,0,1,1)
if(any(a/=expected))error stop 'Exact diagonal inverse changed'
write(*,*)'PASS native FMIND inverse diag(8,4), workspace=',trim(mode)
end program
'''


def routine(text, name):
    start = re.search(r'^\s*SUBROUTINE\s+' + name + r'\b', text, re.I | re.M).start()
    end = re.search(r'^\s*END\s*$', text[start:], re.I | re.M).end() + start
    return text[start:end] + '\n'


def run_check(source, build_parent):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'pams.f').read_text()
    with tempfile.TemporaryDirectory(prefix='pams-workspace-', dir=build_parent) as directory:
        work = Path(directory)
        (work / 'native.f').write_text(''.join(routine(text, name)
                                               for name in ['FMIND', 'VZEROC', 'SCOPYC']))
        (work / 'driver.f90').write_text(DRIVER)
        command = ['gfortran', '-O0', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-fcheck=all', '-ffpe-trap=zero,invalid,overflow',
                   'native.f', str(src / 'smlpck.f'), 'driver.f90', '-lblas', '-o', 'check']
        built = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if built.returncode:
            raise RuntimeError((built.stdout + built.stderr)[-5000:])
        results = {}
        for mode in ['zero', 'finite']:
            result = subprocess.run([str(work / 'check'), mode], cwd=work,
                                    capture_output=True, text=True)
            print(mode, 'returncode', result.returncode, result.stdout.strip())
            results[mode] = result
        if any(result.returncode != 0 for result in results.values()):
            raise AssertionError({mode: (r.returncode, r.stderr[:500])
                                  for mode, r in results.items() if r.returncode})


class TestPamsWorkspace(unittest.TestCase):
    def test_native_inverse_does_not_consume_unused_workspace(self):
        run_check(Path(__file__).resolve().parents[1], Path.cwd())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    run_check(args.source, args.build_parent)
