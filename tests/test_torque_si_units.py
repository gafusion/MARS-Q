"""Native printed-torque conversions use the major-radius magnetic unit."""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get('MARSQ_ROOT', Path(__file__).resolve().parents[1]))


class TorqueSIUnits(unittest.TestCase):
    def test_major_radius_magnetic_unit_is_independent_of_minor_radius(self):
        source = (ROOT/'MarsQ_2FK/torque.f').read_text()
        assignments = re.findall(r'(?m)^ +TORQFAC\s*=[^\r\n]+', source)
        self.assertTrue(assignments)
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as directory:
            work=Path(directory)
            driver='      PROGRAM UNITS\n      IMPLICIT NONE\n      REAL*8 B0EXP,R0EXP,ASPCT,PI,TORQFAC\n      PI=ACOS(-1.D0)\n'
            cases=[(1.7,2.4,2.),(1.7,2.4,7.),(2.3,3.1,2.)]
            for field,radius,aspect in cases:
                driver+=f'      B0EXP={field}D0\n      R0EXP={radius}D0\n      ASPCT={aspect}D0\n'
                for assignment in assignments:driver+=assignment+'\n      WRITE(*,*) 0.23D0*TORQFAC\n'
            driver+='      END\n';(work/'native.f').write_text(driver)
            build=subprocess.run(['gfortran','-O2','-fdefault-real-8','-fdefault-double-8','-ffixed-line-length-none','native.f','-o','oracle'],cwd=work,capture_output=True,text=True)
            self.assertEqual(build.returncode,0,build.stderr)
            run=subprocess.run([str(work/'oracle')],cwd=work,capture_output=True,text=True);self.assertEqual(run.returncode,0,run.stderr)
            observed=[float(x) for x in run.stdout.split()]
            expected=[.23*radius**3*field**2/(4e-7*math.pi) for field,radius,aspect in cases for assignment in assignments]
            for actual,reference in zip(observed,expected):self.assertTrue(math.isclose(actual,reference,rel_tol=3e-13),(actual,reference))
            self.assertEqual(len(observed),len(expected))


if __name__=='__main__':unittest.main()
