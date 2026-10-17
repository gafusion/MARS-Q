"""Focused native thermal-rate and legacy DWK species-sum regressions.

Run with ``make -C MarsQ_2FK test-thermal-species``.  MARS_ORACLE_SOURCE
can select an unmodified source tree to reproduce the missing third species.
Only the real native routines/block are extracted; MPI and profile derivatives
are isolated because these checks use constant, manufactured local profiles.
"""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT)).resolve()
if SOURCE.name != "MarsQ_2FK":
    SOURCE /= "MarsQ_2FK"


def routine(text, name):
    start = re.search(r"^      SUBROUTINE " + name + r"(?:\(|\s*$)", text, re.M)
    if start is None:
        raise AssertionError(f"Native {name} not found")
    end = re.search(r"^      END(?: SUBROUTINE " + name + r")?\s*$",
                    text[start.start():], re.M)
    return text[start.start():start.start() + end.end()] + "\n"


def compile_oracle(work, sources):
    command = [os.environ.get("FC", "gfortran"), "-O0", "-fcheck=all",
               "-fdefault-real-8", "-fdefault-double-8",
               "-ffixed-line-length-none", "-ffree-line-length-none",
               "-I" + str(SOURCE), *map(str, sources), "-o", "oracle"]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=60)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)


class NativeOracle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        scratch = Path(os.environ.get("MARS_TEST_SCRATCH", ROOT / "build")).resolve()
        if scratch == Path("/tmp") or Path("/tmp") in scratch.parents:
            raise AssertionError("Test scratch must be outside /tmp")
        scratch.mkdir(parents=True, exist_ok=True)
        temporary = tempfile.TemporaryDirectory(prefix="thermal-species-", dir=scratch)
        cls.addClassCleanup(temporary.cleanup)
        cls.work = Path(temporary.name)

    def execute(self, inputs):
        result = subprocess.run([str(self.work / "oracle")], cwd=self.work,
                                input=" ".join(map(str, inputs)) + "\n",
                                capture_output=True, text=True, timeout=5)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout, result.stderr


RATE_DRIVER = """
program thermal_rate_oracle
 use rcomdm
 use dimensim
 use globalm
 use kineticm
 use mpienv
 use anisotropicm
 implicit none
 integer :: j,k
!INPUT_METADATA
 NR=2;NRP1=3;ISMPIRUN=1;RANK=ROOT
 read(*,*) KNUSPECIES,NSPECIES,NUMODEL,INUTYPE,NPROFIE,NEXPV,(NUEFFSA(k),k=1,3)
 INCKIN=1;KNTV=21
 allocate(CS(3),CSM(3),GNUI(3),GNUIM(3),GNUE(3),GNUEM(3))
 allocate(NUEFF(3,2,3),FFF(3,2),DFFF(3,2))
 allocate(EPSLONCA(3,2,3),EPSALPHA(3,2,3),AAK(3,2,3))
 allocate(DEPSALPHADPSI(3,2,3),DEPSCDPSI2(3,2,3),DEPSALPHADPSI2(3,2,3))
 allocate(OMEGASN(3,2,3),OMEGAST(3,2,3),OMEGASNA(3,2,3),OMEGASAA(3,2,3),OMEGASCA(3,2,3))
 allocate(ESPECIES_DEN(3,2,3),ESPECIES_TEM(3,2,3))
 CS=[.2d0,.6d0,1d0];CSM=[.4d0,.8d0,1.2d0];ASPCT=3d0;ZTAUA0=5d-6
 GNUI=[.01d0,.02d0,.03d0];GNUIM=[.015d0,.025d0,.035d0]
 GNUE=[.1d0,.2d0,.3d0];GNUEM=[.15d0,.25d0,.35d0]
 EPSLONCA=1d0;EPSALPHA=1d0;AAK=1d0
 ESPECIES_DEN=1d0;ESPECIES_DEN(:,:,2)=2d0;ESPECIES_TEM=1d0
 ESPECIES_M=2d0;ESPECIES_M(2)=1d0/1836.15267343d0;ESPECIES_M(3)=3d0
 ESPECIES_Z=1d0;ESPECIES_Z(2)=-1d0
 ISPECIES_F0=0;B0K=1d0;OMEGACI0=50d0;B0EXP=1d0;R0EXP=1d0
 call KDIAMAG
!READBACK
 do j=1,3
  do k=1,2
   print '(2I3,3ES24.16)',j,k,NUEFF(j,k,:)
  enddo
 enddo
end program
subroutine DFFFDPSI(mode)
 use globalm,only:DFFF
 integer :: mode
 DFFF=0d0
end subroutine
"""


class ThermalRates(NativeOracle):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        text = (SOURCE / "kinetic.f").read_text()
        kernel = routine(text, "KDIAMAG")
        metadata, readback = "", ""
        if "SUBROUTINE KVALIDATESPECIESNU" in text:
            kernel += routine(text, "KVALIDATESPECIESNU")
            kernel += routine((SOURCE / "marsq.f").read_text(), "CALC_NUSTAR")
            readback = "if (KNUSPECIES==1) call CALC_NUSTAR"
        else:
            # Old native KDIAMAG stays unchanged; unused fixture inputs permit
            # a behavioral failure rather than an absent-variable compile error.
            metadata = "integer :: KNUSPECIES\nreal*8 :: NUEFFSA(10)"
        (cls.work / "kernel.f").write_text(kernel)
        (cls.work / "driver.f90").write_text(
            RATE_DRIVER.replace("!INPUT_METADATA", metadata).replace("!READBACK", readback))
        (cls.work / "mpi.f90").write_text(
            "module MPIENV\ninteger :: ISMPIRUN,RANK\ninteger,parameter :: ROOT=0\nend module\n")
        modules = ["globalm.f", "dimensim.f", "rcomdm.f", "kineticm.f", "anisotropicm.f", "torquem.f"]
        compile_oracle(cls.work, [*(SOURCE / name for name in modules),
                                  "mpi.f90", "kernel.f", "driver.f90"])

    def test_physical_rates_and_unchanged_legacy_grids(self):
        # NRL thermal self-collision frequency, T in eV, n in m^-3,
        # ln Lambda=15.  These supplied rates do not assert an all-pair model.
        e, eps0, mp = 1.602176634e-19, 8.8541878128e-12, 1.67262192369e-27
        coefficient = 4 * math.sqrt(math.pi) / 3 * e**2.5 / (4 * math.pi * eps0)**2
        def rate(mass):
            return coefficient * 1e19 * 15 / math.sqrt(mp * mass) / 1000**1.5
        rates = [rate(2), rate(1 / 1836.15267343), rate(3)]
        self.assertAlmostEqual(rates[2] / rates[0], math.sqrt(2 / 3), places=14)
        for model in (0, 1):
            for opt_in, species, supplied in ((0, 2, rates), (1, 2, rates),
                                              (1, 3, rates), (1, 3, [0, 0, 0])):
                with self.subTest(model=model, opt_in=opt_in, species=species):
                    out, err = self.execute([opt_in, species, model, 0, 0, 0, *supplied])
                    self.assertNotIn("STOP", err)
                    rows = [list(map(float, line.split())) for line in out.splitlines()]
                    self.assertEqual(len(rows), 6)
                    for row in rows:
                        j, grid = map(int, row[:2])
                        grid = 2 if j == 1 and grid == 1 else grid
                        coords = [.2, .6, 1] if grid == 1 else [.4, .8, 1.2]
                        factor = 1 if model == 0 else 3 / (2 * coords[j - 1])
                        if opt_in:
                            expected = [value * 5e-6 for value in supplied]
                            if species == 2:
                                expected[2] = 0
                        else:
                            ions = [.01, .02, .03] if grid == 1 else [.015, .025, .035]
                            electrons = [.1, .2, .3] if grid == 1 else [.15, .25, .35]
                            expected = [ions[j - 1], electrons[j - 1], 0]
                        for kp in (2, 0, 1):
                            wanted = expected[kp] * factor
                            self.assertAlmostEqual(row[2 + kp], wanted,
                                                   delta=max(1e-15, abs(wanted) * 3e-15))
                    if opt_in:
                        diagnostic = (self.work / "PROFNUSTAR.OUT").read_text()
                        self.assertIn("DIAGNOSTIC OMITTED", diagnostic)
                        rows = [line.split() for line in
                                (self.work / "THERMAL_COLLISION_RATES.OUT").read_text().splitlines()
                                if line and not line.startswith("%")]
                        self.assertEqual(len(rows), species)
                        for kp, row in enumerate(rows):
                            self.assertEqual(int(row[0]), kp + 1)
                            self.assertAlmostEqual(float(row[1]), supplied[kp], places=12)
                            self.assertAlmostEqual(float(row[2]), supplied[kp] * 5e-6, places=15)

    def test_invalid_active_rates_and_unsupported_profiles_reject(self):
        cases = [([-1, 2, 3], 0, 0, 0, "NONNEGATIVE"),
                 (["NaN", 2, 3], 0, 0, 0, "FINITE"),
                 (["Infinity", 2, 3], 0, 0, 0, "FINITE"),
                 ([1, 2, 0], 1, 0, 0, "POSITIVE THERMAL"),
                 ([1e-320, 2, 3], 1, 0, 0, "POSITIVE ACTIVE"),
                 ([1, 2, 3], 0, 4, 0, "NPROFIE=0 AND NEXPV=0"),
                 ([1, 2, 3], 0, 0, 1, "NPROFIE=0 AND NEXPV=0")]
        for rates, inutype, prof, nexp, message in cases:
            with self.subTest(rates=rates, inutype=inutype, prof=prof, nexp=nexp):
                out, err = self.execute([1, 3, 0, inutype, prof, nexp, *rates])
                self.assertEqual(out, "")
                self.assertIn(message, err)


TORQUE_DRIVER = """
      PROGRAM THERMAL_TORQUE_ORACLE
      IMPLICIT NONE
      INTEGER,PARAMETER :: NR=2
      INTEGER NSPECIES,KNTV,KP,IS,I,INDX,FIDS,ASSIGNFREEFILEUNIT
      INTEGER ISPECIES_F0(3),INDXDWKC(3,10)
      REAL*8 PI,RNTOR,CSM(NR),TSPEC,WORK,FACTOR,POLARITY
      REAL*8 TORQUENTVI(NR),TORQUENTVE(NR),TORQUENTV(NR)
      COMPLEX*16 DWPPARY(NR,30),DWPPERY(NR,30)
      READ(*,*) NSPECIES,ISPECIES_F0(3),POLARITY
      ISPECIES_F0(1:2)=0
      KNTV=21
      PI=ACOS(-1D0)
      RNTOR=-4D0*PI*PI*POLARITY
      CSM=(/.25D0,.75D0/)
      INDXDWKC=-1
      DWPPARY=0
      DWPPERY=0
      DO KP=1,3
         INDXDWKC(KP,3)=10*(KP-1)+3
         INDXDWKC(KP,5)=10*(KP-1)+5
         DO IS=1,NR
            FACTOR=10D0**(IS-1)
            DO I=3,5,2
               INDX=INDXDWKC(KP,I)
               WORK=FACTOR*DBLE(KP)*DBLE(I-2)/4D0
               DWPPARY(IS,INDX)=CMPLX(100D0,-.25D0*WORK,8)
               DWPPERY(IS,INDX)=CMPLX(-35D0,-.75D0*WORK,8)
            ENDDO
C           AN ELECTRON INDEX WITHOUT ITS ION INDEX MUST STILL BE
C           SKIPPED BY THE LEGACY COUPLED CYCLE IN THE FIRST-TWO LOOP.
            INDXDWKC(2,1)=11
            DWPPARY(IS,11)=CMPLX(0D0,-99D0,8)
         ENDDO
      ENDDO
@PRODUCTION@
      DO IS=1,NR
         WRITE(*,'(I4,3ES26.17E3)') IS,TORQUENTVI(IS),
     &         TORQUENTVE(IS),TORQUENTV(IS)
      ENDDO
      END PROGRAM
      INTEGER FUNCTION ASSIGNFREEFILEUNIT()
      INTEGER U
      LOGICAL OPENED
      DO U=20,99
         INQUIRE(UNIT=U,OPENED=OPENED)
         IF (.NOT.OPENED) THEN
            ASSIGNFREEFILEUNIT=U
            RETURN
         ENDIF
      ENDDO
      STOP 'NO FREE FIXTURE UNIT'
      END FUNCTION
"""


class ThermalTorque(NativeOracle):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        body = routine((SOURCE / "kinetic.f").read_text(), "CALCDWKCOMP")
        start = body.index("      IF (KNTV.EQ.21) THEN")
        end = body.index("C     OUTPUT THE PROFILES OF ENERGY DENSITY", start)
        (cls.work / "oracle.f").write_text(
            TORQUE_DRIVER.replace("@PRODUCTION@", body[start:end]))
        compile_oracle(cls.work, ["oracle.f"])

    def test_native_dwk_total_and_resolved_thermal_contributions(self):
        for species, third_f0, polarity in ((2, 0, 1), (3, 0, 1), (3, 0, -1), (3, 1, 1)):
            with self.subTest(species=species, third_f0=third_f0, polarity=polarity):
                out, err = self.execute([species, third_f0, polarity])
                self.assertNotIn("STOP", err)
                totals = [list(map(float, line.split())) for line in out.splitlines()]
                self.assertEqual(len(totals), 2)
                for cell, row in enumerate(totals):
                    factor = 2 * polarity * 10**cell
                    accepted = 6 if species == 3 and third_f0 == 0 else 3
                    self.assertEqual(row, [cell + 1, factor, 2 * factor, accepted * factor])
                    if species == 2:
                        self.assertEqual(row[3], row[1] + row[2], "Legacy two-species sum")
                path = self.work / "NTV_TORQUE_SPECIES.OUT"
                self.assertTrue(path.is_file(), "Native per-species torque report missing")
                rows = [list(map(float, line.split())) for line in path.read_text().splitlines()
                        if line.strip() and not line.startswith("#")]
                self.assertEqual(len(rows), 2 * species)
                resolved = {}
                for row in rows:
                    self.assertEqual(len(row), 4)
                    kp, cell = map(int, row[:2])
                    self.assertNotIn((kp, cell), resolved)
                    wanted = 0 if kp == 3 and third_f0 else 2 * polarity * 10**(cell - 1) * kp
                    self.assertEqual(row, [kp, cell, [.25, .75][cell - 1], wanted])
                    resolved[kp, cell] = row[3]
                for cell, total in enumerate(totals, 1):
                    self.assertEqual(sum(resolved[kp, cell] for kp in range(1, species + 1)), total[3])


if __name__ == "__main__":
    unittest.main()
