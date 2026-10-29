"""Public circular equilibrium: FF' -> prescribed q -> recovered FF'.

The Miller C1 case is specified in gafusion/MARS-Q#65. All observations use
the unmodified executable's native half-mesh NUPLO output and EXPEQ.OUT. Requiring the
recovered NSTTP=1 equilibrium to reproduce q rejects a false convergence
obtained by zeroing an inner coarea.
"""
import math
from pathlib import Path
import subprocess

from support import EXE, NativeTest


def deck(mode):
    return f"""synthetic circular prescribed-q round trip
public analytic Miller boundary
no machine data
native half-mesh NUPLO and EXPEQ observations
 &EQDATA
 NEQDSK=0, NSURF=6, NTCASE=0, NBLOPT=0, NBSOPT=0,
 NCSCAL=4, CSSPEC=0., QSPEC=1.1, NTMF0=0,
 NSTTP={mode}, NFUNC=4, NPPFUN=4, NPROPT=1, NIPR=1,
 NS=32, NT=32, NISO=97, NPSI=128, NCHI=128,
 NV=16, REXT=1., R0W=1., RZ0W=0., NMESHA=0,
 NSMOOTH=0, EPSLON=1e-8, NINMAP=100, NINSCA=100,
 NOPT=0, NPLOT=1, NBAL=0, NIDEAL=0, NEGP=-1,
 MSMAX=12, NRSCAL=0, B0EXP=1., R0EXP=1., RELAX=.2,
 &END
"""


def expeq(s, pprime, function, mode):
    lines = ['.125', '0.', '0.', '257 1 1']
    for j in range(257):
        angle = -math.pi + 2*math.pi*j/256
        lines.append(f'{1+.125*math.cos(angle):.16e} '
                     f'{.125*math.sin(angle):.16e}')
    lines.append(f'{len(s)} {mode}')
    for values in (s, pprime, function):
        lines.extend(f'{value:.16e}' for value in values)
    return '\n'.join(lines)+'\n'


def native_profiles(folder):
    """Read the existing headless NUPLO physical half-mesh profiles.

    WRTPLOT's dimension header precedes title records, flags and
    fixed eight-value records. Unlike the five-digit console dump, its seven
    significant digits retain the q derivatives needed by the current solve.
    NPLOT only writes this data file and opens no viewer.
    """
    lines = iter((folder/'NUPLO').read_text().splitlines())
    metadata = []
    for _ in range(3):
        metadata.extend(map(int, next(lines).split()))
    (insur, _, nchi1, _, npsi1, _, ns1, nt, nt1, ins, inr,
     _, titles, _, _, _, _, _, _, _, _, _, _, _, niso, _, _, _) = metadata
    for _ in range(titles + 3*((npsi1+39)//40)):
        next(lines)

    def take(count):
        values = []
        for _ in range((count+7)//8):
            values.extend(float(v.replace('D', 'E')) for v in next(lines).split())
        if len(values) != count:
            raise AssertionError('Invalid native NUPLO record length')
        return values

    # Native WRTPLOT: scaling/mesh-shape metadata, meshes and plot abscissae.
    for count in ([10]*11 + [nt1, ns1, npsi1, nchi1, niso+1,
                            nt, nt, insur, insur, insur, ins, nt1,
                            nchi1, nt1, niso+1, ns1]):
        take(count)
    s = take(ins)[1:]  # First entry is synthetic axis metadata.
    take(inr)          # Signed real-space abscissae, a different grid.
    q = take(ins)[1:]
    pressure = [-.025*(.2+.8*(1-v*v)**2) for v in s]
    return [s, pressure, q]


def interpolate(x, y, query):
    # Independent piecewise-linear comparison on the two almost identical
    # native meshes. This never replaces the native equilibrium calculation.
    j = 0
    result = []
    for value in query:
        while j < len(x)-2 and x[j+1] < value:
            j += 1
        result.append(y[j]+(value-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    return result


class PrescribedQRoundtrip(NativeTest):
    def run_case(self, name, mode, equilibrium):
        folder = self.work/name
        folder.mkdir()
        (folder/'EXPEQ').write_text(equilibrium)
        result = subprocess.run([str(EXE.resolve())], input=deck(mode),
                                cwd=folder, text=True, capture_output=True,
                                timeout=180)
        (folder/'native.log').write_text(result.stdout+result.stderr)
        self.assertEqual(result.returncode, 0, result.stdout[-2500:]+result.stderr)
        self.assertIn('AFTER MAPPING:', result.stdout)
        self.assertNotIn('NOT CONVERGED', result.stdout)
        profiles = native_profiles(folder)
        self.assertEqual(len(set(map(len, profiles))), 1)
        self.assertLess(profiles[0][0], .01)
        self.assertEqual(profiles[0][-1], 1)
        self.assertTrue(all(b > a for a, b in zip(profiles[0], profiles[0][1:])))
        self.assertTrue(all(math.isfinite(v) for row in profiles for v in row))
        self.assertTrue(all(q > 0 for q in profiles[2]))
        return profiles, (folder/'EXPEQ.OUT').read_text()

    def test_physical_q_and_independent_current_roundtrip(self):
        s = [j/256 for j in range(257)]
        pressure = [-.025*(.2+.8*(1-v*v)**2) for v in s]
        ffprime = [-1.8*(.05+.95*(1-v*v)**2) for v in s]
        target, _ = self.run_case('target', 1, expeq(s, pressure, ffprime, 1))
        fitted, recovered = self.run_case('prescribed', 4,
                                           expeq(*target, 4))
        roundtrip, _ = self.run_case('recovered', 1, recovered)
        for label, values in [('prescribed', fitted), ('recovered', roundtrip)]:
            expected = interpolate(target[0], target[2], values[0])
            errors = [abs(a/b-1) for a, b in zip(values[2], expected)]
            rms = math.sqrt(sum(e*e for e in errors)/len(errors))
            print(f'{label}: max relative q error={max(errors):.9g}, RMS={rms:.9g}')
            # Tighter than the independently measured NS=32/48 q discretisation
            # difference (2.5e-3), with the same fixed 1e-8 native solver tolerance.
            self.assertLess(max(errors), 2e-4)
            self.assertLess(rms, 5e-5)


if __name__ == '__main__':
    import unittest
    unittest.main()
