"""Execute the native tabulated-rotation conversion with independent units.

Physical rad/s input is multiplied by the Alfven time, including a zero axis
with off-axis flow. Normalised input retains its axis-amplitude convention.
The exact production block runs with floating-point traps, before splining.
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


class PhysicalRotationUnits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source_root = Path(os.environ.get("MARS_ORACLE_SOURCE", ROOT))
        source = (source_root / "MarsQ_2FK/marsq.f").read_text()
        start = source.index("C        ZTEMP = MAXVAL(ABS(RPRV))")
        end = source.index("         IF (IRAD.EQ.1) THEN", start)
        block = source[start:end]
        cls.scratch = tempfile.TemporaryDirectory(prefix="rotation-units-")
        cls.work = Path(cls.scratch.name)
        (cls.work / "native.f").write_text(
            "      SUBROUTINE CONVERT(NEXPV,ROTE,RPRV,ZTAUA0)\n"
            "      IMPLICIT NONE\n"
            "      INTEGER NEXPV\n"
            "      REAL*8 ROTE,RPRV(5),ZTAUA0,ZTEMP\n" + block +
            "      END\n")
        (cls.work / "driver.f90").write_text(
            "program driver\nimplicit none\ninteger nexpv\n"
            "real*8 rote,tau,profile(5)\n"
            "read(*,*) nexpv,rote,tau,profile\n"
            "call convert(nexpv,rote,profile,tau)\n"
            "write(*,'(6ES26.17)') rote,profile\nend program\n")
        build = subprocess.run(
            ["gfortran", "-O2", "-fdefault-real-8", "-fdefault-double-8",
             "-ffixed-line-length-none", "-ffpe-trap=invalid,zero,overflow",
             "-fcheck=all", "native.f", "driver.f90", "-o", "oracle"],
            cwd=cls.work, capture_output=True, text=True)
        if build.returncode:
            raise RuntimeError(build.stdout + build.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.scratch.cleanup()

    def convert(self, mode, amplitude, profile, tau=2.7e-6, reject=False):
        data = " ".join(map(str, [mode, amplitude, tau, *profile])) + "\n"
        run = subprocess.run([str(self.work / "oracle")], input=data,
                             capture_output=True, text=True)
        if reject:
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("ZERO AXIS NORMALISATION", run.stdout)
            return
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        result = np.fromstring(run.stdout, sep=" ")
        self.assertEqual(result.shape, (6,))
        self.assertTrue(np.isfinite(result).all())
        return result

    def test_physical_zero_profile_is_finite_zero(self):
        np.testing.assert_array_equal(self.convert(1, 0, [0]*5), np.zeros(6))

    def test_physical_zero_axis_keeps_signed_off_axis_flow(self):
        profile = np.array([0, 200, -350, 70, 0])
        np.testing.assert_allclose(self.convert(1, 99, profile),
                                   np.r_[0, profile*2.7e-6], rtol=1e-15)

    def test_physical_nonzero_axis_retains_units(self):
        profile = np.array([-120, -80, 30, 60, 0])
        np.testing.assert_allclose(self.convert(1, 99, profile),
                                   np.r_[profile[0], profile]*2.7e-6,
                                   rtol=1e-15)

    def test_normalised_nonzero_axis_retains_amplitude(self):
        profile = np.array([-2., -1., .4, .1, 0])
        np.testing.assert_allclose(self.convert(0, -.03, profile),
                                   np.r_[-.03, -.03*profile/profile[0]],
                                   rtol=1e-15)

    def test_normalised_zero_amplitude_with_reference_is_zero(self):
        np.testing.assert_array_equal(self.convert(0, 0, [2, 1, -1, 0, 0]),
                                      np.zeros(6))

    def test_normalised_identically_zero_profile_is_defined(self):
        np.testing.assert_array_equal(self.convert(0, 0, [0]*5), np.zeros(6))

    def test_normalised_zero_reference_cannot_set_finite_amplitude(self):
        self.convert(0, .03, [0]*5, reject=True)

    def test_normalised_zero_reference_cannot_scale_off_axis_flow(self):
        self.convert(0, 0, [0, 1, -2, .1, 0], reject=True)


if __name__ == "__main__":
    unittest.main()
