"""The native hot-ion diagnostic must have a defined value without a hot species."""
import os
import re
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get("MARS_ORACLE_SOURCE", Path(__file__).resolve().parents[1]))


class HotIonDiagnostic(unittest.TestCase):
    def test_native_diagnostic_with_two_and_three_species(self):
        source = (ROOT / "MarsQ_2FK/kinetic.f").read_text()
        begin = re.search(r"(?m)^ +LAMH = 0\.", source).start()
        end = begin + re.search(r"(?m)^ +FREQK\(JS,1\)", source[begin:]).start()
        block = source[begin:end]
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
            work = Path(tmp)
            (work / "native.f").write_text(
                "      SUBROUTINE HOTDIAG(NSPECIES,MLMAX,RLM,SLAM0,LAMH)\n"
                "      IMPLICIT NONE\n"
                "      INTEGER NSPECIES,MLMAX,L\n"
                "      REAL*8 RLM(MLMAX),SLAM0(MLMAX,NSPECIES),LAMH\n"
                + block + "      END\n")
            (work / "driver.f90").write_text(
                "program driver\nimplicit none\ninteger ns\n"
                "real*8 h,r(3),s(3,3)\nread(*,*) ns\n"
                "r=[-1d0,0d0,1d0]\ns(:,1)=1d0\ns(:,2)=2d0\n"
                "s(:,3)=[31d0,37d0,41d0]\n"
                "call hotdiag(ns,3,r,s(:,1:ns),h)\n"
                "write(*,*) h\nend program\n")
            build = subprocess.run(
                ["gfortran", "-O0", "-fcheck=all", "-ffixed-line-length-none",
                 "native.f", "driver.f90", "-o", "oracle"],
                cwd=work, capture_output=True, text=True)
            self.assertEqual(build.returncode, 0, build.stderr)
            for species, expected in [(2, 0.), (3, 37.)]:
                with self.subTest(species=species):
                    run = subprocess.run([str(work / "oracle")],
                                         input=str(species)+"\n", cwd=work,
                                         capture_output=True, text=True)
                    self.assertEqual(run.returncode, 0, run.stderr)
                    self.assertEqual(float(run.stdout), expected)


if __name__ == "__main__":
    unittest.main()
