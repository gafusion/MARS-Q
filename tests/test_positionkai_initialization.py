"""Skipped kinetic surfaces must export deterministic position diagnostics."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class PositionKaiInitialization(unittest.TestCase):
    def test_skipped_surfaces_are_zero_and_computed_rows_survive(self):
        source_root = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT))
        source = (source_root / "MarsQ_2FK/kinetic.f").read_text()
        start = source.index(
            "      IF (ISMPIRUN.EQ.0) THEN\n      ALLOCATE( FREQKSURF")
        end = source.index("      ENDIF", start) + len("      ENDIF")
        driver = """      PROGRAM CHECK_ALLOCATION
      IMPLICIT NONE
      INTEGER NRP1,NCHI,ISMPIRUN
      PARAMETER(NRP1=31,NCHI=17)
      REAL*8,ALLOCATABLE::FREQKSURF(:,:,:),POSITIONKAI(:,:)
      ISMPIRUN=0
""" + source[start:end] + """
      IF (ANY(POSITIONKAI.NE.0.0D0)) STOP 7
      POSITIONKAI(4:NRP1,:)=3.25D0
      IF (ANY(POSITIONKAI(1:3,:).NE.0.0D0)) STOP 8
      IF (ANY(POSITIONKAI(4:NRP1,:).NE.3.25D0)) STOP 9
      END
"""
        with tempfile.TemporaryDirectory(prefix="positionkai-") as directory:
            work = Path(directory)
            (work / "native.f").write_text(driver)
            build = subprocess.run(
                ["gfortran", "-O2", "-ffixed-line-length-none", "native.f",
                 "-o", "oracle"], cwd=work, capture_output=True, text=True)
            self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
            # glibc fills newly allocated storage, exposing the original
            # uninitialized read independently of prior allocator contents.
            environment = dict(os.environ, MALLOC_PERTURB_="165",
                               GLIBC_TUNABLES="glibc.malloc.tcache_count=0")
            run = subprocess.run([str(work / "oracle")], cwd=work,
                                 env=environment, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)


if __name__ == "__main__":
    unittest.main()
