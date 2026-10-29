"""Radial integral of a precession layer narrower than one MARS cell.

A rotation zero puts the trapped ell=0 resonance x = -n omega_E/(n omega_D)
inside the energy domain only on a radial layer of width ~ omega_D/omega_E'.
Compile the executed KIA_TRAP kernel, place that layer at ten positions inside
one radial cell, assemble the node and half-mesh responses with the radial
weights of the DWK assembly (kinetic.f), and compare with an independent
scipy double integral of the same integrand over the same radial interval.

A second case lets the precession frequency vary linearly in s.  For fixed
energy x the denominator x*omega_D(s) + n*omega_E(s) - i*nu is then linear in
s, the radial integral is a closed-form logarithm, and a narrow layer weighs
1/|n omega_E' + x n omega_D'|.  Freezing omega_D across the cell (weight
1/|n omega_E'|) misses this by the drift slope.
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

NR = 40              # cells on [0, 1]
GRAD = 7500.0        # d(n omega_E)/ds; layer width 15/GRAD = 0.1 cell
NU = 0.1             # Krook rate; Im D = -NU x^-1.5 (n omega_D = 1)
RTMP1 = -0.06        # RNTOR B0K/OMEGACI0 for RNTOR=-3, OMEGACI0=50
A_N, B_T = 0.2, 0.3  # dF/dpsi = exp(-x) (A_N + (x - 5/2) B_T)
JCELL = 20           # the zero moves from node JCELL+1 towards JCELL+2
THETAS = np.arange(10) / 10.0
WINDOW = 4           # nodes/cells on each side of the zero entering the sum


def integrand(x, u):
    """KIA_TRAP thermal ell=0 integrand at u = n omega_E (OMEGA = 0)."""
    f1 = np.exp(-x) * (A_N + (x - 2.5) * B_T)
    f2 = -np.exp(-x)
    return -2 * x**2.5 * (RTMP1 * f1 + u * f2) / (x + u - 1j * NU * x**-1.5)


@lru_cache(maxsize=None)
def surface(u):
    """Exact energy integral at one radius."""
    pts = [-u] if -u > 0 else None
    out = []
    for part in (np.real, np.imag):
        val = quad(lambda x: part(integrand(x, u)), 0, 60, points=pts,
                   limit=400, epsabs=1e-13, epsrel=1e-11)[0]
        val += quad(lambda x: part(integrand(x, u)), 60, np.inf, epsabs=1e-13)[0]
        out.append(val)
    return complex(*out)


def radial(ua, ub):
    """int_{ua}^{ub} surface(u) du / GRAD, the layer split out explicitly."""
    edges = sorted({ua, ub, *[e for e in (-40.0, 0.0, 40.0) if ua < e < ub]})
    return sum(piece(a, b) for a, b in zip(edges[:-1], edges[1:])) / GRAD


@lru_cache(maxsize=None)
def piece(a, b):
    pts = [p for p in np.linspace(-30, 0, 31) if a < p < b] or None
    out = 0j
    for part in (np.real, np.imag):
        out += (1j if part is np.imag else 1) * quad(
            lambda u: part(surface(u)), a, b, points=pts, limit=400,
            epsabs=1e-11, epsrel=1e-9)[0]
    return out


def run_kernel(test, nr, grad, nu, jcell, dslope):
    """Compile the executed KIA_TRAP and return its node/half-mesh rows.

    The precession frequency at s is (1 + dslope*(s - s0)) in units of its
    value at the rotation zero s0, on every pitch of the ZOMEGADT table.
    """
    source_root = Path(os.environ.get('MARS_ORACLE_SOURCE', ROOT))
    source = (source_root / 'MarsQ_2FK/anisotropic.f').read_text()
    kernels = []
    for name in ('KIA_TRAP', 'KENERGYMESH', 'KTRAPDRIFTAT'):
        found = re.search(rf'^      SUBROUTINE {name}(?:\(|\s*$)', source, re.M)
        if found is None and name == 'KTRAPDRIFTAT':
            continue        # sources that freeze omega_D do not read the table
        start = found.start()
        end = source.index('END SUBROUTINE ' + name, start)
        kernels.append(source[start:end] + 'END SUBROUTINE ' + name + '\n')
    scratch = ROOT / 'build'
    scratch.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=scratch, prefix='layer-oracle-') as td:
        work = Path(td)
        (work / 'kernels.f').write_text(''.join(kernels))
        driver = DRIVER.replace('@NR@', str(nr)).replace('@GRAD@', repr(grad))
        driver = driver.replace('@NU@', repr(nu)).replace('@JCELL@', str(jcell))
        driver = driver.replace('@DSLOPE@', repr(float(dslope)))
        (work / 'oracle.f').write_text('\n'.join('      ' + s for s in driver.splitlines()) + '\n')
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
    return np.array([[float(t) for t in ln.split()] for ln in output.splitlines()])


def window_sums(rows, it, mode, lo, hi, h):
    """Node sum over nodes lo..hi and half-mesh sum over cells lo..hi-1."""
    sel = rows[(rows[:, 0] == it) & (rows[:, 1] == mode)]
    node = {int(r[3]): r[4] + 1j * r[5] for r in sel if r[2] == 1}
    half = {int(r[3]): r[4] + 1j * r[5] for r in sel if r[2] == 2}
    # node JS carries (CSH(JS-1)+CSH(JS))/2, half JS carries CSH(JS)
    return (sum(node[j] * h for j in range(lo, hi + 1)),
            sum(half[j] * h for j in range(lo, hi)))


# omega_D(s) linear in s.  GRAD2 = n omega_E', DSLOPE2 = n omega_D'/(n omega_D)
# at the zero; x*DSLOPE2/GRAD2 = 0.2 at x = 4 as in the ITER pedestal, where
# the layer of the whole drive (x < 10, width 10/GRAD2) fits in one cell.  The
# window keeps omega_D within [0.25, 1.75] of its value at the zero.
NR2, GRAD2, DSLOPE2, JCELL2, WINDOW2 = 200, 2000.0, 100.0, 100, 1
THETAS2 = (0, 3, 7)


def drift_layer_exact(ta, tb):
    """int_{ta}^{tb} ds int_0^inf dx of the KIA_TRAP integrand, t = s - s0.

    For fixed x the denominator x*(1 + DSLOPE2 t) + GRAD2 t - i NU x^-1.5 and
    the numerator are linear in t, so the radial integral is a logarithm.
    """
    def inner(x):
        f1 = np.exp(-x) * (A_N + (x - 2.5) * B_T)
        f2 = -np.exp(-x)
        p, q = -2 * x**2.5 * RTMP1 * f1, -2 * x**2.5 * GRAD2 * f2
        a = x - 1j * NU * x**-1.5
        b = x * DSLOPE2 + GRAD2
        return q / b * (tb - ta) + (p - q * a / b) / b * (
            np.log(a + b * tb) - np.log(a + b * ta))
    pts = sorted({-GRAD2 * t / (1 + DSLOPE2 * t) for t in (ta, tb)
                  if 0 < -GRAD2 * t / (1 + DSLOPE2 * t) < 60})
    out = 0j
    for part, unit in ((np.real, 1), (np.imag, 1j)):
        edges = [0.0, *pts, 60.0]
        val = sum(quad(lambda x: part(inner(x)), a, b, limit=400,
                       epsabs=1e-13, epsrel=1e-11)[0]
                  for a, b in zip(edges[:-1], edges[1:]))
        val += quad(lambda x: part(inner(x)), 60, np.inf, epsabs=1e-13)[0]
        out += unit * val
    return out


class RotationZeroLayer(unittest.TestCase):
    def test_layer_integral_is_independent_of_zero_position(self):
        rows = run_kernel(self, NR, GRAD, NU, JCELL, 0.0)
        cs = np.linspace(0.0, 1.0, NR + 1)
        h = cs[1] - cs[0]
        lo, hi = JCELL + 1 - WINDOW, JCELL + 1 + WINDOW     # 1-based node indices
        report = []
        for it, theta in enumerate(THETAS):
            s0 = cs[JCELL] + theta * h
            for mode in (0, 1):
                xsum, ysum = window_sums(rows, it, mode, lo, hi, h)
                xa, xb = cs[lo - 1] - h / 2, cs[hi - 1] + h / 2
                ya, yb = cs[lo - 1], cs[hi - 1]
                if mode == 0:
                    exact_x = radial(GRAD * (xa - s0), GRAD * (xb - s0))
                    exact_y = radial(GRAD * (ya - s0), GRAD * (yb - s0))
                report.append((mode, xsum.imag / exact_x.imag, ysum.imag / exact_y.imag))
                print(f'theta={theta:.1f} KRADLAYER={mode}: node/exact '
                      f'{xsum / exact_x:.6f}, half/exact {ysum / exact_y:.6f}, '
                      f'Im node/exact {xsum.imag / exact_x.imag:.6f}, '
                      f'Im half/exact {ysum.imag / exact_y.imag:.6f}')
                if mode == 1:
                    for got, ref in ((xsum, exact_x), (ysum, exact_y)):
                        self.assertLess(abs(got / ref - 1), 0.01)
                        self.assertLess(abs(got.imag / ref.imag - 1), 0.01)
        # The legacy point sampling is not a layer integral: its resonant
        # (imaginary) part depends on where the zero sits in the cell and
        # misses the layer by far more than the 1 % tolerance above.
        legacy = np.array([r[1:] for r in report if r[0] == 0])
        self.assertGreater(np.ptp(legacy), 0.02)
        self.assertGreater(np.max(np.abs(legacy - 1)), 0.5)

    def test_layer_weight_includes_precession_slope(self):
        rows = run_kernel(self, NR2, GRAD2, NU, JCELL2, DSLOPE2)
        cs = np.linspace(0.0, 1.0, NR2 + 1)
        h = cs[1] - cs[0]
        lo, hi = JCELL2 + 1 - WINDOW2, JCELL2 + 1 + WINDOW2
        for it in THETAS2:
            s0 = cs[JCELL2] + it / 10.0 * h
            xsum, ysum = window_sums(rows, it, 1, lo, hi, h)
            exact_x = drift_layer_exact(cs[lo - 1] - h / 2 - s0, cs[hi - 1] + h / 2 - s0)
            exact_y = drift_layer_exact(cs[lo - 1] - s0, cs[hi - 1] - s0)
            print(f'theta={it / 10:.1f} omega_D slope: node/exact {xsum / exact_x:.6f}, '
                  f'half/exact {ysum / exact_y:.6f}, Im node/exact '
                  f'{xsum.imag / exact_x.imag:.6f}, Im half/exact '
                  f'{ysum.imag / exact_y.imag:.6f}')
            for got, ref in ((xsum, exact_x), (ysum, exact_y)):
                self.assertLess(abs(got / ref - 1), 0.01)
                self.assertLess(abs(got.imag / ref.imag - 1), 0.01)


DRIVER = '''
program layer_oracle
    use globalm
    use dimensim
    use rcomdm
    use kineticm
    use anisotropicm
    implicit none
    include 'compam.inc'
    integer, parameter :: nrc = @NR@, jcell = @JCELL@
    real*8, parameter :: grad = @GRAD@
    integer :: it, mode, js, kgrid, i
    real*8 :: h, s0, theta, sp, wtail
    real*8, parameter :: dslope = @DSLOPE@
    complex*16 :: zvi
    NR = nrc; NRP1 = nrc + 1
    allocate(RLM(1), ISPECIES_EK(1), CS(NRP1+1), T(NRP1), TM(NRP1))
    allocate(DPSIDS(NRP1), DPSIDSM(NRP1))
    allocate(EPSALPHA(NRP1,2,1), NUEFF(NRP1,2,1), OMEGAE0(NRP1,2))
    allocate(DEPSALPHADPSI(NRP1,2,1), SLAM0(1,1))
    allocate(NLAMK0(NRP1,2), HKMIN(NRP1,2), HKMAX(NRP1,2))
    allocate(LAMK0(NRP1,2,2), ZOMEGADT(NRP1,4,2))
    NLAMK0=2; HKMIN=.4d0; HKMAX=.6d0; WK=.5d0
    LAMK0(:,1,:)=.4d0; LAMK0(:,2,:)=.6d0
    h = 1d0/nrc
    do i = 1, NRP1 + 1
        CS(i) = (i - 1)*h
    end do
    RLM=0; ISPECIES_EK=0; T=1; TM=1; DPSIDS=1; DPSIDSM=1
    EPSALPHA=1; NUEFF=@NU@; DEPSALPHADPSI=0; SLAM0=0
    RNTOR=-3; B0K=1; OMEGACI0=50; OMEGA=dcmplx(0d0,0d0)
    INUTYPE=1; OMEGAB=0; RQK=0
    DRIFT = -50d0/3d0
    ISPECIES_F0=0; PSPECIES_FOWP=0; PSPECIES_FOWT=0
    ESPECIES_M=2; ESPECIES_Z=1
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
            OMEGAE0(i,1) = grad*(CS(i) - s0)/RNTOR
            OMEGAE0(i,2) = grad*(0.5d0*(CS(i) + CS(i+1)) - s0)/RNTOR
            ZOMEGADT(i,:,1) = -50d0/3d0*(1d0 + dslope*(CS(i) - s0))
            ZOMEGADT(i,:,2) = -50d0/3d0*(1d0 + dslope*(0.5d0*(CS(i) + CS(i+1)) - s0))
        end do
        do mode = 0, 1
            KRADLAYER = mode
            do kgrid = 1, 2
                do js = 2 - (kgrid - 1), NR
                    sp = CS(js)
                    if (kgrid == 2) sp = 0.5d0*(CS(js) + CS(js+1))
                    DRIFT = -50d0/3d0*(1d0 + dslope*(sp - s0))
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
