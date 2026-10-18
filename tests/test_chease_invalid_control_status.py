"""Run CHEASE's native aggregate control check and its process exit status.

Only the pure control-check block is linked, with no-op output and elementary
integer helper stubs; no equilibrium solve or private fixture is required.
Run: python3 -m unittest discover -s tests -p test_chease_invalid_control_status.py
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class InvalidControlStatus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source_root = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT))
        source = (source_root / "CheaseMerge/chease.f").read_text()
        start = source.index("         CALL RESETI(ILIMIT,140,0)")
        end = source.index("C*DECK C1S07", start)
        block = source[start:end]
        cls.scratch = tempfile.TemporaryDirectory(prefix="chease-status-")
        cls.work = Path(cls.scratch.name)
        (cls.work / "COMDIM.inc").write_bytes(
            (source_root / "CheaseMerge/COMDIM.inc").read_bytes())
        (cls.work / "check.f").write_text(
            "      SUBROUTINE CHECK(NV,NANAL)\n"
            "      IMPLICIT INTEGER(A-Z)\n"
            "      INCLUDE 'COMDIM.inc'\n"
            "      DIMENSION ILIMIT(140)\n"
            "      NS=1\n      NT=2\n      NPSI=1\n      NCHI=2\n"
            "      NISO=1\n      NSMAX=1\n      MSMAX=1\n"
            "      NTURN=1\n      NBLC0=1\n      NPPR=1\n"
            "      NPPF=1\n      NTNOVA=2\n      NMGAUS=1\n"
            "      NBAL=0\n      NCSCAL=1\n      NDIFPS=0\n"
            "      NDIFT=0\n      NFUNC=1\n      NIDEAL=0\n"
            "      NIPR=1\n      NOPT=0\n      NPLOT=0\n"
            "      NBLOPT=0\n      NPP=1\n      NPPFUN=1\n"
            "      NPROFZ=0\n      NPRPSI=0\n      NRSCAL=0\n"
            "      NSTTP=1\n      NSURF=1\n      NSYM=0\n"
            "      NTCASE=0\n      NTEST=0\n      NTMF0=0\n"
            "      NPROPT=1\n      NBSOPT=0\n      NBSTRP=1\n"
            "      NBSFUN=1\n" + block)
        (cls.work / "driver.f90").write_text(
            "program driver\ninteger nv,nanal\nread(*,*) nv,nanal\n"
            "call check(nv,nanal)\nprint *, 'VALID RETURN'\nend program\n"
            "subroutine output(n)\ninteger n\nend subroutine\n"
            "subroutine reseti(a,n,v)\ninteger n,v,a(n)\na=v\nend subroutine\n"
            "integer function issum(n,a,stride)\ninteger n,stride,a(*)\n"
            "issum=sum(a(1:n*stride:stride))\nend function\n")
        subprocess.run(["gfortran", "-O0", "-fcheck=all", "-std=legacy",
                        "-ffixed-line-length-none", "check.f", "driver.f90",
                        "-o", "check"], cwd=cls.work, check=True,
                       capture_output=True, text=True)

    @classmethod
    def tearDownClass(cls):
        cls.scratch.cleanup()

    def run_control(self, nv, nanal):
        return subprocess.run([str(self.work / "check")],
                              input=f"{nv} {nanal}\n", text=True,
                              capture_output=True)

    def test_valid_control_returns_success(self):
        run = self.run_control(1, 0)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertIn("VALID RETURN", run.stdout)

    def test_exceeded_vacuum_capacity_fails(self):
        run = self.run_control(361, 0)
        self.assertIn("RECOMPILE WITH LARGER NPV", run.stdout)
        self.assertNotIn("VALID RETURN", run.stdout)
        self.assertNotEqual(run.returncode, 0, run.stdout + run.stderr)

    def test_invalid_analysis_selector_fails(self):
        run = self.run_control(1, 2)
        self.assertIn("WRONG VALUE FOR NANAL", run.stdout)
        self.assertNotIn("VALID RETURN", run.stdout)
        self.assertNotEqual(run.returncode, 0, run.stdout + run.stderr)


if __name__ == "__main__":
    unittest.main()
