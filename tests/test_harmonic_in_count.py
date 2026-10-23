"""Execute the native HARMONIC.IN reader on generated bounce-harmonic lists.

The file holds a count followed by that many harmonics. Every count >= 1 must
replace the default list, including lists of 100 or more entries; a count
below one must stop with a diagnostic instead of being ignored.
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class HarmonicInCount(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source_root = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT))
        source = (source_root / "MarsQ_2FK/kinetic.f").read_text()
        start = source.index("      INQUIRE(FILE='HARMONIC.IN', EXIST=EXISTS)")
        close = source.index("         CLOSE (HARMOFILE)\n      ENDIF\n", start)
        end = close + len("         CLOSE (HARMOFILE)\n      ENDIF\n")
        block = source[start:end]
        cls.scratch = tempfile.TemporaryDirectory(prefix="harmonic-in-")
        cls.work = Path(cls.scratch.name)
        (cls.work / "native.f").write_text(
            "      INTEGER FUNCTION ASSIGNFREEFILEUNIT()\n"
            "      ASSIGNFREEFILEUNIT = 77\n"
            "      END\n"
            "      SUBROUTINE GETRLM(MLMAX,RLM)\n"
            "      IMPLICIT NONE\n"
            "      INTEGER MLMAX,J,HARMOFILE,HARMONUM,ASSIGNFREEFILEUNIT\n"
            "      REAL*8, ALLOCATABLE :: RLM(:)\n"
            "      LOGICAL EXISTS\n" + block +
            "      END\n")
        (cls.work / "driver.f90").write_text(
            "program driver\nimplicit none\n"
            "interface\nsubroutine getrlm(mlmax,rlm)\ninteger mlmax\n"
            "real*8, allocatable :: rlm(:)\nend subroutine\nend interface\n"
            "integer mlmax,j\nreal*8, allocatable :: rlm(:)\n"
            "mlmax = 5\nallocate(rlm(mlmax))\n"
            "rlm = [(dble(j-3), j=1,mlmax)]\n"
            "call getrlm(mlmax,rlm)\n"
            "write(*,'(I0)') mlmax\n"
            "write(*,'(I0)') (nint(rlm(j)), j=1,mlmax)\nend program\n")
        build = subprocess.run(
            ["gfortran", "-O2", "-ffixed-line-length-none", "-fcheck=all",
             "native.f", "driver.f90", "-o", "oracle"],
            cwd=cls.work, capture_output=True, text=True)
        if build.returncode:
            raise RuntimeError(build.stdout + build.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.scratch.cleanup()

    def read(self, lines):
        with tempfile.TemporaryDirectory(dir=self.work) as case:
            Path(case, "HARMONIC.IN").write_text(
                "".join(f"{x}\n" for x in lines))
            return subprocess.run([str(self.work / "oracle")], cwd=case,
                                  capture_output=True, text=True)

    def harmonics(self, ells):
        run = self.read([len(ells), *ells])
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        values = [int(v) for v in run.stdout.split()
                  if v.lstrip("-").isdigit()]
        self.assertEqual(values[0], len(ells))
        return values[1:]

    def test_short_list_replaces_default(self):
        ells = list(range(-8, 9))
        self.assertEqual(self.harmonics(ells), ells)

    def test_list_of_100_or_more_replaces_default(self):
        for ells in (list(range(-50, 50)), list(range(-60, 61))):
            self.assertEqual(self.harmonics(ells), ells)

    def test_count_below_one_stops(self):
        for count in (0, -17):
            run = self.read([count])
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("HARMONIC.IN COUNT < 1", run.stdout + run.stderr)


if __name__ == "__main__":
    unittest.main()
