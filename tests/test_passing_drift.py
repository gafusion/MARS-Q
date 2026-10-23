"""Passing magnetic drift (KPASSDRIFT) in the native KIA_PASS energy integral.

Compile the unchanged production KIA_PASS, KPASSROOTS and KENERGYMESH with the
owning module declarations and native splines. The oracle integrates the
physical passing response independently with scipy,

    sum_sigma int_0^inf N(x) / (sigma*R*sqrt(x) + R_D*Delta*x + n*w_E - w - i*nu(x)) dx,

N(x) = -x**2.5 * (R_D*f_psi + (n*w_E - w)*f_E), nu = H3 (INUTYPE=0) or
1.25*H3*x**-1.5 (INUTYPE=1), in its zero-width limit by principal value plus
pi*delta. Delta is the transit-averaged drift coefficient DRIFT; the drift
term is absent for KPASSDRIFT=0, which must not depend on DRIFT at all.
"""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad

ROOT = Path(os.environ.get('MARSQ_ROOT') or Path(__file__).resolve().parents[1])
NATIVE = ROOT / 'MarsQ_2FK'
# Fixture constants (tests/fixtures/passing_drift.f90).
RNTOR, RQK, ELL, OMEGAB, OMEGACI0 = -3.0, 0.7, 1.0, 0.8, 10.0
TRANSIT = OMEGAB*(RNTOR*RQK+ELL)*math.sqrt(2.0)
RD = RNTOR/OMEGACI0
CELL = 1.0+ELL**2/4.0


def production_routine(name):
    source = (NATIVE / 'anisotropic.f').read_text()
    match = re.search(r'^ +SUBROUTINE +' + name + r'\b.*?^ +END SUBROUTINE +' + name + r'\b',
                      source, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if match is None:
        raise ValueError('Missing production routine ' + name)
    return match.group() + '\n'


def numerator(x, rotation):
    return -x**2.5*(RD*(-0.4)+RNTOR*rotation*(-1.0))*math.exp(-x)


def oracle(flag, inutype, drift, rotation, width):
    a = RD*drift if flag else 0.0
    total = 0.0j
    for sigma in (1.0, -1.0):
        def real_denominator(x):
            return sigma*TRANSIT*math.sqrt(x)+a*x+RNTOR*rotation

        def width_of(x):
            return width if inutype == 0 else CELL*width*x**-1.5

        if width > 0.0:
            def part(x, k):
                value = numerator(x, rotation)/complex(real_denominator(x), -width_of(x))
                return value.real if k == 0 else value.imag
            roots = [u*u for u in np.roots([a, sigma*TRANSIT, RNTOR*rotation]
                                           if a else [sigma*TRANSIT, RNTOR*rotation])
                     if abs(u.imag) < 1e-12 and u.real > 0.0]
            edges = sorted({0.0, *roots, 60.0})
            for k in (0, 1):
                value = sum(quad(part, lo, hi, args=(k,), limit=400, epsabs=1e-13,
                                 epsrel=1e-12)[0] for lo, hi in zip(edges, edges[1:]))
                total += value if k == 0 else 1j*value
            continue
        # Zero width: 1/(D-i0) = PV(1/D) + i*pi*delta(D).
        roots = sorted(u.real**2 for u in np.roots([a, sigma*TRANSIT, RNTOR*rotation]
                                                  if a else [sigma*TRANSIT, RNTOR*rotation])
                       if abs(u.imag) < 1e-12 and u.real > 0.0)
        windows, start = [], 0.0
        for root in roots:
            half = 0.25*min(root, *(abs(root-other) for other in roots if other != root), 1.0)
            windows.append((start, root-half, None))
            windows.append((root-half, root+half, root))
            start = root+half
        windows.append((start, 60.0, None))
        for lo, hi, root in windows:
            if root is None:
                total += quad(lambda x: numerator(x, rotation)/real_denominator(x), lo, hi,
                              limit=400, epsabs=1e-13, epsrel=1e-12)[0]
                continue
            slope = sigma*TRANSIT/(2.0*math.sqrt(root))+a

            def regular(x):
                if abs(x-root) < 1e-9*root:
                    return numerator(x, rotation)/slope
                return numerator(x, rotation)*(x-root)/real_denominator(x)
            total += quad(regular, lo, hi, weight='cauchy', wvar=root, limit=400)[0]
            total += 1j*math.pi*numerator(root, rotation)/abs(slope)
    return total


class PassingDriftTests(unittest.TestCase):
    def test_native_passing_drift_matches_independent_energy_integral(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', ROOT)) as tmp:
            work = Path(tmp)
            routines = work / 'passing.f'
            routines.write_text(production_routine('KIA_PASS') +
                                production_routine('KENERGYMESH') +
                                "      SUBROUTINE SET_TEST_OMEGA()\n"
                                "      IMPLICIT NONE\n"
                                "      INCLUDE 'compam.inc'\n"
                                "      OMEGA=(0.D0,0.D0)\n"
                                "      END\n")
            modules = ['globalm.f', 'dimensim.f', 'rcomdm.f', 'kineticm.f', 'anisotropicm.f']
            command = [os.environ.get('FC', 'gfortran'), '-O0', '-fcheck=all', '-std=legacy',
                       '-fdefault-real-8', '-fdefault-double-8', '-ffixed-line-length-none',
                       '-ffree-line-length-none', '-fallow-argument-mismatch', '-I' + str(NATIVE),
                       *map(str, (NATIVE / name for name in modules)), str(routines),
                       str(NATIVE / 'splines.f'),
                       str(Path(__file__).with_name('fixtures') / 'passing_drift.f90'),
                       '-o', str(work / 'passing')]
            compiled = subprocess.run(command, cwd=work, text=True, capture_output=True, timeout=120)
            self.assertEqual(compiled.returncode, 0, compiled.stdout + compiled.stderr)
            worst = {}
            for nepk, tolerance in ((1601, 5.e-4), (6401, 2.e-5)):
                run = subprocess.run([str(work / 'passing'), str(nepk)], cwd=work, text=True,
                                     capture_output=True, timeout=60)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                worst[nepk], changed = self.check(run.stdout, tolerance)
        # The narrow high-energy drift root (x ~ 16) converges with the native mesh.
        print('passing drift max oracle error NEPK1601 %.3e NEPK6401 %.3e, min drift effect %.3e'
              % (worst[1601], worst[6401], changed))
        self.assertLess(worst[6401], worst[1601])

    def check(self, stdout, tolerance):
        rows = []
        for line in stdout.splitlines():
            flag, inutype, drift, rotation, width, real, imag = line.split()
            rows.append((int(flag), int(inutype), float(drift), float(rotation), float(width),
                         complex(float(real), float(imag))))
        self.assertEqual(len(rows), 2*(4*2*3+4*2))
        off = {}
        errors = []
        for flag, inutype, drift, rotation, width, value in rows:
            self.assertTrue(math.isfinite(abs(value)))
            if flag == 0:
                key = (inutype, rotation, width)
                # KPASSDRIFT=0 ignores DRIFT exactly (legacy path).
                self.assertEqual(off.setdefault(key, value), value)
            reference = oracle(flag, inutype, drift, rotation, width)
            error = abs(value-reference)/abs(reference)
            errors.append(error)
            with self.subTest(flag=flag, inutype=inutype, drift=drift, rotation=rotation,
                              width=width, tolerance=tolerance):
                self.assertLess(error, tolerance, (value, reference))
        changed = min(abs(v-off[(i, r, w)])/abs(off[(i, r, w)])
                      for f, i, d, r, w, v in rows if f == 1 and d != 0.0)
        self.assertGreater(changed, 0.1, 'drift does not change the response')
        return max(errors), changed


if __name__ == '__main__':
    unittest.main()
