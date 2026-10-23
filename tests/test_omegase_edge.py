"""Native electron diamagnetic profiles preserve a constant-gradient edge."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT=Path(os.environ.get('MARSQ_ROOT',Path(__file__).resolve().parents[1]))


class ElectronDiamagneticBoundary(unittest.TestCase):
    def test_analytic_constant_profile_has_zero_boundary_derivative(self):
        source=(ROOT/'MarsQ_2FK/marsq.f').read_text()
        begin=source.index('      DO J=2,NR',source.index('C     COMPUTE ADDITIONAL EQUILIBRIUM PROFILES'))
        end=source.index('      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN',begin)
        driver='''      PROGRAM EDGE
      IMPLICIT NONE
      INTEGER J,NR,NRP1
      PARAMETER(NR=4,NRP1=5)
      REAL*8 H1,H2,ZTEMP,B0K,OMEGACI0,ALPHAP,ETASTI,ETASTE
      REAL*8 CS(5),RHO(5),RHOM(5),PEQ(5),PEQM(5),PPEQ(5),PPEQM(5)
      REAL*8 DPSIDS(5),DPSIDSM(5),DLNRHO(5),DLNRHOM(5)
      REAL*8 OMEGASI(5),OMEGASE(5),OMEGASIM(5),OMEGASEM(5)
      REAL*8 DOMEGASI(5),DOMEGASE(5),DOMEGASIM(5),DOMEGASEM(5)
      CS=(/0D0,.1D0,.35D0,.6D0,1D0/)
      B0K=2D0
      OMEGACI0=4D0
      ALPHAP=.4D0
      ETASTI=1D0
      ETASTE=1D0
      RHO=1D0
      RHOM=1D0
      PEQ=2D0
      PEQM=2D0
      PPEQ=3D0
      PPEQM=3D0
      DPSIDS=1D0
      DPSIDSM=1D0
      OMEGASE=-999D0
'''+source[begin:end]+'''
      IF(ANY(ABS(OMEGASE-.9D0).GT.1D-13)) STOP 1
      IF(ANY(ABS(OMEGASI+.6D0).GT.1D-13)) STOP 2
      IF(ANY(ABS(DOMEGASEM).GT.1D-13)) STOP 3
      IF(ANY(ABS(DOMEGASIM).GT.1D-13)) STOP 4
      END
'''
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as directory:
            work=Path(directory);(work/'native.f').write_text(driver)
            build=subprocess.run(['gfortran','-O0','-fcheck=all','-fdefault-real-8','-fdefault-double-8','-ffixed-line-length-none','native.f','-o','oracle'],cwd=work,capture_output=True,text=True)
            self.assertEqual(build.returncode,0,build.stderr)
            run=subprocess.run([str(work/'oracle')],cwd=work,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)


if __name__=='__main__':unittest.main()
