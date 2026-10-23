"""FMIND (PAMERA block inversion) must fail the run, with a nonzero exit code, on a singular block."""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

SRC = Path(os.environ.get("MARSQ_ROOT") or Path(__file__).resolve().parents[1]) / "MarsQ_2FK" / "pams.f"

DRIVER = """      PROGRAM DRV
      IMPLICIT NONE
      INTEGER N,NDIM,I,J
      PARAMETER (N=4,NDIM=4)
      COMPLEX*16 A(NDIM,N),AORIG(NDIM,N),D(N),SCR(N,2*N+1),DET,P
      REAL*8 EPS,ERR
      CHARACTER*8 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      A = (0.D0,0.D0)
      DO I=1,N
        A(I,I) = CMPLX(2.D0+I,1.D0)
        IF (I.LT.N) A(I,I+1) = CMPLX(0.5D0,-0.25D0)
      END DO
      IF (ARG.EQ.'singular') A(:,3) = (0.D0,0.D0)
      AORIG = A
      EPS = 1.D-199
      CALL FMIND(A,D,N,NDIM,SCR,DET,EPS,0,1,241)
C     independent oracle: A*INV(A) = 1 for the regular matrix
      ERR = 0.D0
      DO I=1,N
        DO J=1,N
          P = SUM(AORIG(I,:)*A(:,J))
          IF (I.EQ.J) P = P - 1.D0
          ERR = MAX(ERR,ABS(P))
        END DO
      END DO
      PRINT *, 'INVERSE RESIDUAL', ERR
      IF (ERR.GT.1.D-13) ERROR STOP 'bad inverse'
      END
"""


def routine(source, name):
    """One production routine of the fixed-form source, up to its END line."""
    start = re.search(r"^ +SUBROUTINE +" + name + r"\b", source, re.M | re.I).start()
    end = re.search(r"^ +END *$", source[start:], re.M).end()
    return source[start:start + end] + "\n"


class PameraSingularBlock(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR"))
        self.addCleanup(tmp.cleanup)
        self.work = Path(tmp.name)
        source = SRC.read_text()
        code = DRIVER + "".join(routine(source, n) for n in ("FMIND", "VZEROC", "SCOPYC"))
        (self.work / "fmind.f").write_text(code)
        build = subprocess.run([os.environ.get("FC", "gfortran"), "-O0", "-ffixed-line-length-none",
                                "-std=legacy", "-fdefault-real-8", "-fdefault-double-8",
                                "fmind.f", "-llapack", "-o", "fmind"],
                               cwd=self.work, capture_output=True, text=True)
        self.assertEqual(build.returncode, 0, build.stderr)

    def run_driver(self, arg):
        return subprocess.run([str(self.work / "fmind"), arg], cwd=self.work, capture_output=True, text=True)

    def test_regular_block_is_inverted(self):
        result = self.run_driver("regular")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_singular_block_fails_the_run(self):
        result = self.run_driver("singular")
        self.assertIn("ZGETRF ERROR: IPSI,IFAIL", result.stdout)
        self.assertNotEqual(result.returncode, 0, "MARS-F returned 0 after a fatal ZGETRF failure")


if __name__ == "__main__":
    unittest.main()
