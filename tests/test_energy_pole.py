"""Energy integral across a narrow INUTYPE>=1 resonance at moderate NEPK.

With the rate nu*x**-1.5 the Lorentzian of the trapped precession, trapped
bounce and passing transit denominators is narrower than the energy mesh at
low collisionality.  The two-point Gauss sum then aliases the pole and is
3-80 % off at NEPK 401 for the circular benchmark.  Here the executed KIA_TRAP
and KIA_PASS kernels are compiled and compared with exact integrals: the
circular-benchmark pitch sums (fixture), synthetic bounce and transit
resonances, and the unchanged INUTYPE=0 and high-collisionality results.
"""
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).resolve().parent / 'fixtures' / 'circular_energy_pole.json'
NEPK = 401


def numerator(x, fac, rd, base):
    """KIA numerator of the synthetic cases: fac*x**2.5*(rd*dF/dpsi + base*dF/dE)."""
    return fac * x**2.5 * np.exp(-x) * (rd * (.2 + (x - 2.5) * .3) - base)


def exact(num, den, poles):
    """int_0^inf num/den, split at the poles on geometric intervals in their widths."""
    points = {0.0, 60.0}
    for xr, width in poles:
        for k in range(-2, 12):
            for sign in (-1, 1):
                x = xr + sign * width * 10.0**k
                if 0 < x < 60:
                    points.add(x)
    points = sorted(points)
    total = 0j
    for a, b in zip(points[:-1], points[1:]):
        for part in (np.real, np.imag):
            value = quad(lambda x: part(num(x) / den(x)), a, b, limit=400,
                         epsabs=0, epsrel=1e-12)[0]
            total += value if part is np.real else 1j * value
    tail = quad(lambda x: num(x) / den(x).real, 60, np.inf)[0]
    return total + tail


class EnergyPole(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source_root = Path(os.environ.get('MARS_ORACLE_SOURCE', ROOT))
        source = (source_root / 'MarsQ_2FK/anisotropic.f').read_text()
        kernels = []
        for name in re.findall(r'^      SUBROUTINE (K\w+)', source, re.M):
            if name in ('KIA_TRAP', 'KIA_PASS', 'KENERGYMESH', 'KTRAPDRIFTAT'):
                start = re.search(rf'^      SUBROUTINE {name}(?:\(|\s*$)', source, re.M).start()
                end = source.index('END SUBROUTINE ' + name, start)
                kernels.append(source[start:end] + 'END SUBROUTINE ' + name + '\n')
        scratch = ROOT / 'build'
        scratch.mkdir(exist_ok=True)
        cls.tmp = tempfile.TemporaryDirectory(dir=scratch, prefix='energy-pole-')
        cls.work = Path(cls.tmp.name)
        (cls.work / 'kernels.f').write_text(''.join(kernels))
        (cls.work / 'oracle.f').write_text('\n'.join('      ' + s for s in DRIVER.splitlines()) + '\n')
        mods = ['globalm.f', 'dimensim.f', 'rcomdm.f', 'kineticm.f', 'anisotropicm.f']
        command = ['gfortran', '-O2', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-I' + str(source_root / 'MarsQ_2FK')]
        command += [str(source_root / 'MarsQ_2FK' / f) for f in mods]
        command += [str(source_root / 'MarsQ_2FK/splines.f'), 'kernels.f', 'oracle.f', '-o', 'oracle']
        compiled = subprocess.run(command, cwd=cls.work, capture_output=True, text=True)
        if compiled.returncode != 0:
            cls.tmp.cleanup()
            raise RuntimeError(compiled.stdout + compiled.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_sets(self, sets):
        """Each set: header values of the driver and (g, W) pitch rows; returns sum W*ZVIF."""
        lines = [str(len(sets))]
        for head, rows in sets:
            lines.append(' '.join(repr(float(v)) if isinstance(v, float) else str(v)
                                  for v in head) + f' {len(rows)}')
            lines += [f'{g!r} {w!r}' for g, w in rows]
        (self.work / 'sets.txt').write_text('\n'.join(lines) + '\n')
        output = subprocess.run(['./oracle'], cwd=self.work, check=True,
                                capture_output=True, text=True).stdout
        return [complex(*map(float, line.split())) for line in output.splitlines()]

    def test_circular_benchmark_pitch_sums(self):
        data = json.loads(FIXTURE.read_text())
        ids = list(data['cases'])
        sets = []
        for cid in ids:
            c = data['cases'][cid]
            surface = data['surfaces'][c['surface']]
            # kind inutype kpassdrift ell ntor qk omegab drift ndb omegae nu f1 f2
            head = [0, NEPK, 1, 0, 0.0, float(c['n']), 0.0, 0.0, 0.0, 0.0, c['omega_E'],
                    c['nu0'], c['F1'], c['F2']]
            sets.append((head, list(zip(surface['g'], surface['W']))))
        for cid, value in zip(ids, self.run_sets(sets)):
            ratio = value.imag / data['cases'][cid]['exact']
            print(f'{cid} NEPK{NEPK}: MARS/exact = {ratio:.6f}')
            # high collisionality is resolved by the mesh: keep the quadrature accuracy
            tol = 1e-4 if 'c3000' in cid or 'c0p1_aE0' in cid else 1e-2
            with self.subTest(case=cid):
                self.assertLess(abs(ratio - 1), tol)

    def test_synthetic_bounce_and_transit_resonances(self):
        # trapped ell=-1 bounce harmonic with precession (NDB=1) and passing
        # transit with N*QK+ELL = -1.  After the division the denominator is
        # P*sqrt(x)+Q*x+N*OMEGAE-i*1.25*nu*x**-1.5, P = OMEGAB*(N*QK+ELL)*SQRT(2),
        # Q = 1/2 (trapped precession) or 0 (passing, no drift term).
        nu, ntor = 2e-4, -3.0
        cases = [  # kind, OMEGAB, QK, OMEGAE, numerator factor, P
            (1, 2.0, 0.0, -1.0 / 3.0, -2.0, -2.0 * np.sqrt(2.0)),
            (2, 1.2, 0.0, -1.0 / 3.0, -1.0, -1.2 * np.sqrt(2.0)),
        ]
        sets = []
        for kind, omegab, qk, omegae, _, _ in cases:
            head = [kind, NEPK, 1, 1 if kind == 2 else 0, -1.0, ntor, qk, omegab,
                    -1.0 / 6.0, 1.0, omegae, nu, 0.0, 0.0]
            sets.append((head, [(-1.0 / 6.0, 1.0)]))
        for (kind, omegab, qk, omegae, fac, p), value in zip(cases, self.run_sets(sets)):
            base = ntor * omegae
            q = .5 if kind == 1 else 0.0
            den = lambda x, p=p, q=q, base=base: (p * np.sqrt(x) + q * x + base
                                                  - 1j * nu * 1.25 * x**-1.5)
            roots = np.roots([q, p, base]) if q else np.array([-base / p])
            poles = [(s * s, nu * 1.25 * s**-3 / abs(p / (2 * s) + q))
                     for s in roots[np.isreal(roots)].real if s > 0]
            self.assertTrue(poles)
            expected = exact(lambda x: numerator(x, fac, ntor, base), den, poles)
            print(f'kind {kind} P {p:.4f} poles {poles}: MARS={value} exact={expected}')
            # the resonant (imaginary) part is small beside the real part for ell=-1
            with self.subTest(kind=kind, P=p):
                self.assertLess(abs(value.real / expected.real - 1), 1e-3)
                self.assertLess(abs(value.imag / expected.imag - 1), 1e-3)

    def test_energy_independent_rate_unchanged(self):
        # INUTYPE=0 keeps its analytic pole subtraction (constant nu)
        nu, ntor, omegae = 2e-2, -3.0, 1.0 / 3.0
        base = ntor * omegae
        head = [1, NEPK, 0, 0, 0.0, ntor, 0.0, 0.0, -1.0 / 6.0, 1.0, omegae, nu, 0.0, 0.0]
        value, = self.run_sets([(head, [(-1.0 / 6.0, 1.0)])])
        expected = exact(lambda x: numerator(x, -2.0, ntor, base),
                         lambda x: .5 * x + base - 1j * nu, [(2.0, nu / .5)])
        print(f'INUTYPE=0: MARS={value} exact={expected}')
        self.assertLess(abs(value / expected - 1), 1e-4)


DRIVER = '''
program energy_pole
    use globalm
    use dimensim
    use rcomdm
    use kineticm
    use anisotropicm
    implicit none
    include 'compam.inc'
    integer :: nset, iset, kind, npitch, ip, k, kpd
    real*8 :: ell, f1, f2, g, w, rd
    complex*16 :: zvif, total
    NR=2; NRP1=3
    allocate(RLM(1), ISPECIES_EK(1), CS(3), T(3), DPSIDS(3))
    allocate(EPSALPHA(3,2,1), NUEFF(3,2,1), OMEGAE0(3,2))
    allocate(DEPSALPHADPSI(3,2,1), SLAM0(1,1))
    ISPECIES_EK=0; CS=[0d0,.5d0,1d0]; T=1; DPSIDS=1
    EPSALPHA=1; DEPSALPHADPSI=0; SLAM0=0
    B0K=1; OMEGACI0=1; OMEGA=0
    ISPECIES_F0=0; PSPECIES_FOWP=0; PSPECIES_FOWT=0
    ESPECIES_M=2; ESPECIES_Z=1
    open(10, file='sets.txt', status='old')
    read(10,*) nset
    do iset=1,nset
        read(10,*) kind, NEPK, INUTYPE, kpd, ell, RNTOR, RQK, OMEGAB, DRIFT, PSPECIES_NDB(1), OMEGAE0(2,1), NUEFF(2,1,1), f1, f2, npitch
        RLM=ell; NEPK2=2*NEPK-1
        allocate(EPK(NEPK),ZEPK(NEPK2),ZEPKO(NEPK2),ZEPKN(NEPK2))
        allocate(ZRESU(NEPK2,3,0:6,1))
        call KENERGYMESH
        ZRESU=0
        rd=RNTOR
        do k=1,NEPK2
            ! KIA_* form F_psi/F_E (gafusion/MARS-Q#38): keep F_E nonzero.
            ! Circular: F_E = -1e-100 adds a numerator OMEGAN0*F_E far below
            ! the tolerance. Synthetic: F_E floored above underflow.
            if (kind.eq.0) then
                ZRESU(k,2,1,1)=-.5d0*exp(-ZEPKN(k))*(f1+(ZEPKN(k)-2.5d0)*f2)/rd
                ZRESU(k,2,2,1)=-1d-100
            else
                ZRESU(k,:,1,1)=exp(-min(ZEPKN(k),700d0))*(.2d0+(ZEPKN(k)-2.5d0)*.3d0)
                ZRESU(k,:,2,1)=-exp(-min(ZEPKN(k),700d0))
                if (kind.eq.2) then
                    ZRESU(k,3,:,1)=0
                    ZRESU(k,3,2,1)=-1d-100
                end if
            end if
        end do
        total=0
        do ip=1,npitch
            read(10,*) g, w
            DRIFT=g
            if (kind.eq.2) then
                call KIA_PASS(2,1,1,0,1,.5d0,0,zvif)
            else
                call KIA_TRAP(2,1,1,0,1,.5d0,0,zvif)
            end if
            total=total+w*zvif
        end do
        print '(2es26.17)', real(total), aimag(total)
        deallocate(EPK,ZEPK,ZEPKO,ZEPKN,ZRESU)
    end do
end program
subroutine PRECOBS_SCALAR
    error stop 'Unexpected precession observer'
end subroutine
subroutine PRECOBS_SAMPLE
    error stop 'Unexpected precession observer'
end subroutine
subroutine PRECOBS_INTEGRAL
    error stop 'Unexpected precession observer'
end subroutine
subroutine KELLTRACESELECT(js,kgrid,otrace)
    implicit none
    integer :: js,kgrid
    logical :: otrace
    otrace=.false.
end subroutine
subroutine WRITEELLTRACEENERGY
    error stop 'Unexpected trace'
end subroutine
subroutine WRITEELLTRACEPITCH
    error stop 'Unexpected trace'
end subroutine
'''


if __name__ == '__main__':
    unittest.main()
