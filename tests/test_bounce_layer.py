"""Radial integral of a trapped bounce layer (ell != 0) narrower than one cell.

For a bounce harmonic ell the KIA_TRAP denominator at normalised energy x is
ell*omega_b*sqrt(x) + n*omega_D*x + n*omega_E - omega - i*nu(x).  It vanishes
for thermal x only on the radial layer where n*omega_E - omega lies between
zero and -ell*omega_b*sqrt(x) - n*omega_D*x, of width ~ell*omega_b/|n omega_E'|.
Point sampling at the nodes and half-mesh points of the DWK assembly misses or
hits such a layer depending on where it sits in the cell.  KRADLAYER=2
replaces the point value by the mean over the radial cell with omega_E,
omega_b and omega_D linear on each half cell.

The synthetic profiles here are linear in s, so for fixed x the radial
integral of the integrand is a closed-form logarithm and the window sum of
the cell means must equal it up to the energy quadrature.  The test compiles
the executed KIA_TRAP and checks

* a narrow ell=1 layer: KRADLAYER=2 matches the exact radial integral within
  1 % for ten positions of the layer in the cell, and point sampling
  (KRADLAYER=0 and 1) has an O(1) error in the resonant (imaginary) part;
* linear omega_b and omega_D: the cell mean carries their slopes, which move
  the exact result by more than the tolerance;
* bitwise invariance: KRADLAYER=1 leaves ell != 0 equal to KRADLAYER=0, and
  KRADLAYER=2 leaves ell = 0 equal to KRADLAYER=1.
"""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
from functools import lru_cache

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[1]

RTMP1 = -0.06        # RNTOR B0K/OMEGACI0 for RNTOR=-3, OMEGACI0=50
A_N, B_T = 0.2, 0.3  # dF/dpsi = exp(-x) (A_N + (x - 5/2) B_T)
NU = 0.1             # NUEFF; Im D = -NU (1 + ell^2/4) x^-1.5 for INUTYPE=1
TOL = 0.01

# Case A: narrow layer.  ell*omega_b = 1, no precession, n omega_E' = 7500:
# the layer (x < 25, |u| < 5) is 7e-4 wide, 0.03 of a cell.
CASE_A = dict(nr=40, grad=7500.0, ell=1.0, b0=1.0, bs=0.0, p0=0.0, ds=0.0,
              jcell=20, window=3)
# Case B: omega_b and n omega_D linear in s.  The relative change across the
# window stays below 0.8, the slope terms are ~6 % of n omega_E' at x = 4.
CASE_B = dict(nr=40, grad=500.0, ell=1.0, b0=1.0, bs=8.0, p0=0.5, ds=8.0,
              jcell=20, window=3)
# ell = 0 invariance run: the precession layer of test_rotation_zero_layer.
CASE_0 = dict(nr=40, grad=7500.0, ell=0.0, b0=0.0, bs=0.0, p0=1.0, ds=0.0,
              jcell=20, window=3)


def run_kernel(test, case):
    """Compile the executed KIA_TRAP; rows it, mode, kgrid, js, Re, Im."""
    source_root = Path(os.environ.get('MARS_ORACLE_SOURCE', ROOT))
    source = (source_root / 'MarsQ_2FK/anisotropic.f').read_text()
    kernels = []
    for name in ('KIA_TRAP', 'KENERGYMESH', 'KTRAPDRIFTAT'):
        found = re.search(rf'^      SUBROUTINE {name}(?:\(|\s*$)', source, re.M)
        test.assertIsNotNone(found, name)
        start = found.start()
        end = source.index('END SUBROUTINE ' + name, start)
        kernels.append(source[start:end] + 'END SUBROUTINE ' + name + '\n')
    driver = DRIVER.replace('@NU@', repr(NU))
    for key, value in case.items():
        driver = driver.replace('@' + key.upper() + '@', repr(value))
    scratch = ROOT / 'build'
    scratch.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=scratch, prefix='bounce-layer-') as td:
        work = Path(td)
        (work / 'kernels.f').write_text(''.join(kernels))
        (work / 'oracle.f').write_text(
            '\n'.join('      ' + s for s in driver.splitlines()) + '\n')
        mods = ['globalm.f', 'dimensim.f', 'rcomdm.f', 'kineticm.f', 'anisotropicm.f']
        command = ['gfortran', '-O2', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-I' + str(source_root / 'MarsQ_2FK')]
        command += [str(source_root / 'MarsQ_2FK' / f) for f in mods]
        command += ['kernels.f', 'oracle.f', '-o', 'oracle']
        compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
        test.assertEqual(compiled.returncode, 0, compiled.stdout + compiled.stderr)
        output = subprocess.run(['./oracle'], cwd=work, check=True,
                                capture_output=True, text=True).stdout
    return output.splitlines()


def parse(lines):
    return np.array([[float(t) for t in ln.split()] for ln in lines])


def window_sums(rows, it, mode, lo, hi, h):
    """Node sum over nodes lo..hi and half-mesh sum over cells lo..hi-1."""
    sel = rows[(rows[:, 0] == it) & (rows[:, 1] == mode)]
    node = {int(r[3]): r[4] + 1j * r[5] for r in sel if r[2] == 1}
    half = {int(r[3]): r[4] + 1j * r[5] for r in sel if r[2] == 2}
    return (sum(node[j] * h for j in range(lo, hi + 1)),
            sum(half[j] * h for j in range(lo, hi)))


def exact(case, ta, tb, frozen=False):
    """int_{ta}^{tb} dt int_0^inf dx of the KIA_TRAP integrand, t = s - s0.

    D = b(t) sqrt(x) + p(t) x + G t - i NU (1 + ell^2/4) x^-1.5 with
    b = b0 (1 + bs t), p = p0 (1 + ds t); the numerator is linear in t.
    frozen=True keeps b and p at their values at the zero (bs = ds = 0).
    """
    g, b0, p0 = case['grad'], case['b0'], case['p0']
    bs, ds = (0.0, 0.0) if frozen else (case['bs'], case['ds'])
    nuell = NU * (1 + case['ell']**2 / 4)

    def inner(x):
        f1 = np.exp(-x) * (A_N + (x - 2.5) * B_T)
        f2 = -np.exp(-x)
        p, q = -2 * x**2.5 * RTMP1 * f1, -2 * x**2.5 * g * f2
        a = b0 * np.sqrt(x) + p0 * x - 1j * nuell * x**-1.5
        bb = g + b0 * bs * np.sqrt(x) + p0 * ds * x
        return q / bb * (tb - ta) + (p - q * a / bb) / bb * (
            np.log(a + bb * tb) - np.log(a + bb * ta))

    pts = set()
    for t in (ta, tb):       # Re D(x, t) = 0: (p0 (1+ds t)) y^2 + b0 (1+bs t) y + g t
        for y in np.roots([p0 * (1 + ds * t), b0 * (1 + bs * t), g * t]):
            if abs(y.imag) < 1e-12 and 0 < y.real**2 < 60:
                pts.add(float(y.real**2))
    edges = [0.0, *sorted(pts), 60.0]
    out = 0j
    for part, unit in ((np.real, 1), (np.imag, 1j)):
        val = sum(quad(lambda x: part(inner(x)), a, b, limit=400,
                       epsabs=1e-13, epsrel=1e-11)[0]
                  for a, b in zip(edges[:-1], edges[1:]))
        val += quad(lambda x: part(inner(x)), 60, np.inf, epsabs=1e-13)[0]
        out += unit * val
    return out


@lru_cache(maxsize=None)
def kernel_lines(key):
    return tuple(run_kernel(unittest.TestCase(), dict(key)))


def cached(case):
    return list(kernel_lines(tuple(sorted(case.items()))))


class BounceLayer(unittest.TestCase):
    def compare(self, case, modes, label, frozen_check=False):
        rows = parse(cached(case))
        nr, jcell, w = case['nr'], case['jcell'], case['window']
        cs = np.linspace(0.0, 1.0, nr + 1)
        h = cs[1] - cs[0]
        lo, hi = jcell + 1 - w, jcell + 1 + w        # 1-based node indices
        ratios = {m: [] for m in modes}
        for it in range(10):
            s0 = cs[jcell] + it / 10.0 * h
            ex = exact(case, cs[lo - 1] - h / 2 - s0, cs[hi - 1] + h / 2 - s0)
            ey = exact(case, cs[lo - 1] - s0, cs[hi - 1] - s0)
            if frozen_check:
                fx = exact(case, cs[lo - 1] - h / 2 - s0, cs[hi - 1] + h / 2 - s0,
                           frozen=True)
                # the slopes move the exact layer integral beyond the tolerance
                self.assertGreater(abs(fx.imag / ex.imag - 1), 3 * TOL)
            for mode in modes:
                xsum, ysum = window_sums(rows, it, mode, lo, hi, h)
                rx, ry = xsum / ex, ysum / ey
                ratios[mode].append((xsum.imag / ex.imag, ysum.imag / ey.imag,
                                     abs(rx - 1), abs(ry - 1)))
                print(f'{label} theta={it / 10:.1f} KRADLAYER={mode}: node/exact '
                      f'{rx:.6f}, half/exact {ry:.6f}, Im node/exact '
                      f'{xsum.imag / ex.imag:.6f}, Im half/exact {ysum.imag / ey.imag:.6f}')
        return {m: np.array(v) for m, v in ratios.items()}

    def test_narrow_bounce_layer(self):
        r = self.compare(CASE_A, (0, 2), 'narrow ell=1')
        layer = r[2]
        self.assertLess(np.max(np.abs(layer[:, :2] - 1)), TOL)
        self.assertLess(np.max(layer[:, 2:]), TOL)
        # point sampling depends on where the layer sits and misses it by O(1)
        legacy = r[0][:, :2]
        self.assertGreater(np.ptp(legacy), 0.02)
        self.assertGreater(np.max(np.abs(legacy - 1)), 0.5)

    def test_bounce_layer_carries_frequency_slopes(self):
        r = self.compare(CASE_B, (2,), 'sloped ell=1', frozen_check=True)
        layer = r[2]
        self.assertLess(np.max(np.abs(layer[:, :2] - 1)), TOL)
        self.assertLess(np.max(layer[:, 2:]), TOL)

    def test_switch_leaves_other_paths_bitwise_unchanged(self):
        def by_mode(lines, mode):
            return [ln.split(None, 2)[2] for ln in lines if int(ln.split()[1]) == mode]
        bounce = cached(CASE_A)
        self.assertEqual(by_mode(bounce, 0), by_mode(bounce, 1))
        self.assertNotEqual(by_mode(bounce, 1), by_mode(bounce, 2))
        precession = cached(CASE_0)
        self.assertEqual(by_mode(precession, 1), by_mode(precession, 2))
        self.assertNotEqual(by_mode(precession, 0), by_mode(precession, 1))


DRIVER = '''
program bounce_layer_oracle
    use globalm
    use dimensim
    use rcomdm
    use kineticm
    use anisotropicm
    implicit none
    include 'compam.inc'
    integer, parameter :: nrc = @NR@, jcell = @JCELL@
    real*8, parameter :: tgrad = @GRAD@, tell = @ELL@, tb0 = @B0@, tbs = @BS@
    real*8, parameter :: tp0 = @P0@, tds = @DS@
    integer :: it, mode, js, kgrid, i
    real*8 :: h, s0, theta, sp, sq2, wtail
    complex*16 :: zvi
    NR = nrc; NRP1 = nrc + 1
    allocate(RLM(1), ISPECIES_EK(1), CS(NRP1+1), T(NRP1), TM(NRP1))
    allocate(DPSIDS(NRP1), DPSIDSM(NRP1))
    allocate(EPSALPHA(NRP1,2,1), NUEFF(NRP1,2,1), OMEGAE0(NRP1,2))
    allocate(DEPSALPHADPSI(NRP1,2,1), SLAM0(1,1))
    allocate(NLAMK0(NRP1,2), HKMIN(NRP1,2), HKMAX(NRP1,2))
    allocate(LAMK0(NRP1,2,2), ZOMEGADT(NRP1,4,2), ZOMEGABT(NRP1,4,2))
    NLAMK0=2; HKMIN=.4d0; HKMAX=.6d0; WK=.5d0
    LAMK0(:,1,:)=.4d0; LAMK0(:,2,:)=.6d0
    h = 1d0/nrc
    do i = 1, NRP1 + 1
        CS(i) = (i - 1)*h
    end do
    RLM=tell; ISPECIES_EK=0; T=1; TM=1; DPSIDS=1; DPSIDSM=1
    EPSALPHA=1; NUEFF=@NU@; DEPSALPHADPSI=0; SLAM0=0
    RNTOR=-3; B0K=1; OMEGACI0=50; OMEGA=dcmplx(0d0,0d0)
    INUTYPE=1; RQK=0
    ISPECIES_F0=0; PSPECIES_FOWP=0; PSPECIES_FOWT=0; PSPECIES_NDB=1
    ESPECIES_M=2; ESPECIES_Z=1
    ! RTMP2 = ell*sqrt(2*EPSALPHA*M1/MKP) = ell*sqrt(2): OMEGAB = b/(ell*sqrt(2))
    sq2 = sqrt(2d0)
    if (abs(tell) > 0.1d0) sq2 = tell*sqrt(2d0)
    NEPK=801; NEPK2=2*NEPK-1
    allocate(EPK(NEPK),ZEPK(NEPK2),ZEPKO(NEPK2),ZEPKN(NEPK2))
    allocate(ZRESU(NEPK2,3,0:6,1))
    call KENERGYMESH
    ZRESU=0
    do i=1,NEPK2
        ! FLOOR: upstream KIA_TRAP divides dF/dE by dF/dE, which is
        ! exp(-ZEPKN); the tan energy mesh reaches ZEPKN ~ 2400 at the
        ! top point, where exp underflows to zero and 0/0 becomes NaN.
        ! The floor is below every value the oracle resolves.
        wtail = exp(-ZEPKN(i))
        if (wtail .eq. 0d0) wtail = TINY(1d0)
        ZRESU(i,2,1,1)=wtail*(.2d0+(ZEPKN(i)-2.5d0)*.3d0)
        ZRESU(i,2,2,1)=-wtail
    end do
    do it = 0, 9
        theta = it/10d0
        s0 = CS(jcell+1) + theta*h
        do i = 1, NRP1
            sp = 0.5d0*(CS(i) + CS(i+1))
            OMEGAE0(i,1) = tgrad*(CS(i) - s0)/RNTOR
            OMEGAE0(i,2) = tgrad*(sp - s0)/RNTOR
            ZOMEGABT(i,:,1) = tb0*(1d0 + tbs*(CS(i) - s0))/sq2
            ZOMEGABT(i,:,2) = tb0*(1d0 + tbs*(sp - s0))/sq2
            ! n omega_D x = NDB*DRIFT*RTMP1*x = p x
            ZOMEGADT(i,:,1) = tp0*(1d0 + tds*(CS(i) - s0))/(-0.06d0)
            ZOMEGADT(i,:,2) = tp0*(1d0 + tds*(sp - s0))/(-0.06d0)
        end do
        do mode = 0, 2
            KRADLAYER = mode
            do kgrid = 1, 2
                do js = 2 - (kgrid - 1), NR
                    sp = CS(js)
                    if (kgrid == 2) sp = 0.5d0*(CS(js) + CS(js+1))
                    OMEGAB = tb0*(1d0 + tbs*(sp - s0))/sq2
                    DRIFT = tp0*(1d0 + tds*(sp - s0))/(-0.06d0)
                    call KIA_TRAP(js,kgrid,1,0,1,.5d0,0,zvi)
                    print '(4I6,2ES26.17)', it, mode, kgrid, js, real(zvi), aimag(zvi)
                end do
            end do
        end do
    end do
end program
subroutine SPLINE1D
    error stop 'Unexpected pole subtraction'
end subroutine
'''


if __name__ == '__main__':
    unittest.main()
