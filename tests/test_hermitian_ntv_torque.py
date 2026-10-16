#!/usr/bin/env python3
"""KNTV=21 torque from the same-drive (Hermitian) kinetic work.

The NTV torque is proportional to Im(dW_k).  KJPFILL forms the non-adiabatic pressures from
the KH drive Hbar of the solved perturbed field.  Paired with that same drive,
Im(dW_k) is Im(K) rho |Hbar|^2 (rho = 1/omega_b for passing orbits): zero for a
reactive kernel, and independent of the displacement-based test rows of the
weak contraction.  Those rows reproduce the drive only when the discrete field
satisfies ideal induction; otherwise the weak torque, and the torque from its
dissipative-kernel part, move with the test rows.

Requires gfortran.  KJPFILL, the DWK component cache writer/reader and
CALCDWKCOMP with its helpers are extracted unchanged from the selected source
(default: this checkout) and driven by a manufactured case; no plasma inputs
are needed.  The torque coefficient (n or 2n, issues 21/22) is measured from a
field that satisfies ideal induction, so the test does not depend on it.
One orbit harmonic and one poloidal harmonic on a two-cell radial mesh; the
drive coefficients differ between surfaces and grids so that the node and
centre placement of each cell is exercised.
"""
import argparse
import math
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROUTINES = ['KJPFILL', 'KHERMFILL', 'KHERMADD', 'KHERMQUAD', 'SETDWKCOMPMAT',
            'ALLOCATEDWKCOMPMAT', 'DEALLOCATEDWKCOMPMAT',
            'WRITE_SURFACE_QUANTITIES', 'READ_SURFACE_QUANTITIES',
            'READ_SURFACE_HERMITIAN', 'FILLMATDWKCOMP', 'CALCPRECOMP',
            'CALCDWKPROF', 'CALCDWKCOMP']

G = 0.8 + 0.3j                      # parallel pressure projection (KG)
H = [0.5 - 0.2j, 0.3 + 0.4j, -0.6 + 0.1j, 0.2 + 0.5j, 0.4 - 0.3j]  # KH: X1 X2 B1 B2 B3
X1 = [0.0, 1.0 + 0.5j, 0.9 - 0.2j]
X2 = [0.7 - 0.4j, 0.6 + 0.3j]
B1 = [0.3 + 0.2j, 0.5 - 0.1j, 0.4 + 0.6j]
B2 = [0.2 - 0.3j, -0.1 + 0.4j]
B3 = [0.6 + 0.1j, 0.3 - 0.5j]
# Change of the field that ideal induction would give for the same
# displacement; the test rows see B + eps*DB.
DB1 = [0.3 - 0.2j, 0.2 + 0.6j, -0.3 + 0.6j]
DB2 = [0.5 + 0.3j, 0.4 + 0.3j]
DB3 = [-0.5 + 0.6j, 0.6 + 0.3j]
CS = [0.05, 0.4, 1.0]
CSM = [0.2, 0.7]
DPSIDS = [0.0, 0.9, 1.3]
DPSIDSM = [0.6, 1.1]
OMEGAB = 0.37
# Trapped precession harmonic with drift reversal: drive and projection at the
# reversal pitch, and its bounce frequency.
H0 = [0.2 + 0.6j, -0.4 + 0.1j, 0.3 + 0.3j, 0.5 - 0.2j, -0.1 + 0.4j]
G0 = 0.5 - 0.6j
OMEGAB0 = 0.29
# Passing transit harmonic with a logarithmic resonance at pitch SLAM0 (inside
# the passing range 0..HKMIN=1), evaluated at pitch LAM.
HS = [-0.3 + 0.5j, 0.6 + 0.2j, 0.1 - 0.4j, -0.2 - 0.3j, 0.5 + 0.1j]
GS = -0.2 + 0.7j
OMEGABS = 0.41
LAM, SLAM0, SF0V, B0K = 0.3, 0.5, 0.45 - 0.35j, 1.3
PTRAPH = 0.9
RNTOR = -2.0
PI2 = 2.0 * math.pi * math.pi
# T = c Im(dW_k).  The source uses c = -2n/(4 pi^2); issue 21 / PR 22 propose
# -n/(4 pi^2).  The test measures c and accepts either value.
COEFFICIENTS = (-2.0 * RNTOR / (4.0 * math.pi ** 2), -1.0 * RNTOR / (4.0 * math.pi ** 2))


def scale(js, kgrid):
    """Drive coefficient factor on surface js (1-based) of grid kgrid."""
    return 1.0 + 0.1 * js + 0.3j * (kgrid - 1)


def points(cell):
    """(weight, psi', block factor, drive) at lower node, centre, upper node."""
    i = cell                      # 0-based cell, nodes i and i+1
    nr = len(CSM)
    y_mid = [0.5 * (X1[i] + X1[i + 1]), X2[i], 0.5 * (B1[i] + B1[i + 1]), B2[i], B3[i]]
    y_lo = [X1[i], X2[i] * CSM[i] / CS[i], B1[i], B2[i], B3[i] * CS[i] / CSM[i]]
    y_hi = [X1[i + 1], X2[i] * CSM[i] / CS[i + 1], B1[i + 1], B2[i],
            B3[i] * CS[i + 1] / CSM[i]]
    # The axis node and the edge node have no coefficients of their own;
    # the cell centre coefficients stand in for them.
    lo = (DPSIDSM[i], scale(i + 1, 2)) if i == 0 else (DPSIDS[i], scale(i + 1, 1))
    hi = (DPSIDSM[i], scale(i + 1, 2)) if i + 1 == nr else (DPSIDS[i + 1], scale(i + 2, 1))
    return [(0.5 * PTRAPH, lo[0], lo[1], y_lo),
            (1.0 - PTRAPH, DPSIDSM[i], scale(i + 1, 2), y_mid),
            (0.5 * PTRAPH, hi[0], hi[1], y_hi)]


def drive(factor, y):
    return factor * sum(h * v for h, v in zip(H, y))


def same_drive_work(kernel, cell):
    """Analytic same-drive work, in the sign of DWPPARY+DWPPERY."""
    return -PI2 * sum(w * 0.5 * psi * kernel / OMEGAB * abs(drive(f, y)) ** 2
                      for w, psi, f, y in points(cell))


def pressure(kernel, cell):
    return sum(w * kernel * G * drive(f, y) for w, psi, f, y in points(cell))


def expected_torque(kernel):
    return [(-same_drive_work(kernel, c)).imag for c in range(len(CSM))]


def expected_reversal_torque(kernel, kernel0):
    """Trapped orbits, rho = 1/(2 omega_b): the regular pitch term minus its
    value at the reversal pitch, plus the reversal contribution itself."""
    out = []
    for cell in range(len(CSM)):
        w_cell = 0.0
        for w, psi, f, y in points(cell):
            h2 = abs(drive(f, y)) ** 2
            h02 = abs(f * sum(c * v for c, v in zip(H0, y))) ** 2
            w_cell += w * 0.5 * psi * (kernel * h2 / (2 * OMEGAB)
                                       + (kernel0 - kernel) * h02 / (2 * OMEGAB0))
        out.append((PI2 * w_cell).imag)
    return out


def test_row(cell, eps):
    """KYV2/KYPPARA row.  eps=0: the weak contraction equals the same-drive
    work (ideal induction); eps>0: the test drive comes from B + eps*DB."""
    d = same_drive_work(1.0, cell) / (PI2 * X2[cell].conjugate() * pressure(1.0, cell))
    if eps:
        i = cell
        b = [B1[i] + eps * DB1[i], B2[i] + eps * DB2[i], B3[i] + eps * DB3[i],
             B1[i + 1] + eps * DB1[i + 1]]
        y_test = [0.5 * (X1[i] + X1[i + 1]), X2[i], 0.5 * (b[0] + b[3]), b[1], b[2]]
        _, _, f, y = points(cell)[1]
        d *= (drive(f, y_test) / drive(f, y)).conjugate()
    return d


def expected_log_torque(kernel, sf0):
    """Passing orbits, rho = 1/omega_b: regular term, the logarithmic term
    with its value at SLAM0 subtracted, and that value integrated over pitch."""
    c2 = math.log(abs(LAM - SLAM0)) * sf0
    c3 = (SLAM0 * (math.log(SLAM0) - 1) + (1 - SLAM0) * (math.log(1 - SLAM0) - 1)) \
        * sf0 / (B0K * math.sqrt(math.pi))
    out = []
    for cell in range(len(CSM)):
        w_cell = 0.0
        for w, psi, f, y in points(cell):
            h2 = abs(drive(f, y)) ** 2
            hs2 = abs(f * sum(c * v for c, v in zip(HS, y))) ** 2
            w_cell += w * 0.5 * psi * ((kernel + c2) * h2 / OMEGAB
                                       + (c3 - c2) * hs2 / OMEGABS)
        out.append((PI2 * w_cell).imag)
    return out


def fortran_complex(z):
    z = complex(z)
    return '(%.17e,%.17e)' % (z.real, z.imag)


def array(values):
    return '[' + ','.join(fortran_complex(z) for z in values) + ']'


def real_array(values):
    return '[' + ','.join('%.17e' % v for v in values) + ']'


def driver(kernel, eps, kernel0=None, sf0=None):
    """kernel0 is None: passing transit harmonic.  Otherwise: trapped
    precession harmonic with drift reversal and reversal kernel kernel0.
    sf0 given: passing harmonic with the logarithmic resonance at SLAM0."""
    rows = [test_row(cell, eps) for cell in range(2)]
    trapped = kernel0 is not None
    fill = []
    for kgrid in (1, 2):
        for js in (1, 2, 3):
            f = scale(js, kgrid)
            fill.append(f'''kgrid={kgrid}; js={js}
VX1={fortran_complex(f*H[0])}; VX2={fortran_complex(f*H[1])}; VQ1={fortran_complex(f*H[2])}
VQ2={fortran_complex(f*H[3])}; VQ3={fortran_complex(f*H[4])}
VX10={fortran_complex(f*H0[0])}; VX20={fortran_complex(f*H0[1])}; VQ10={fortran_complex(f*H0[2])}
VQ20={fortran_complex(f*H0[3])}; VQ30={fortran_complex(f*H0[4])}
SVX10={fortran_complex(f*HS[0])}; SVX20={fortran_complex(f*HS[1])}; SVQ10={fortran_complex(f*HS[2])}
SVQ20={fortran_complex(f*HS[3])}; SVQ30={fortran_complex(f*HS[4])}
call fill_surface(js,kgrid)''')
    fills = '\n'.join(fill)
    if trapped:
        species = 'PSPECIES_NTD(1)=1; RLMV=0; SLAMD0=0.5d0; KPART=0'
        calls = 'call KJPFILL(js,js,kgrid,1d0,0,1)\ncall KJPFILL(js,js,kgrid,0d0,0,4)'
        vi0 = fortran_complex(kernel0)
    else:
        species = 'PSPECIES_NP(1)=1; RLMV=1; SLAMD0=0; KPART=1'
        calls = 'call KJPFILL(js,js,kgrid,1d0,1,1)'
        vi0 = fortran_complex(0)
    if sf0 is not None:
        species += f'; SLAM0(1,1)={SLAM0}d0; SF0(1,1,0,1)={fortran_complex(sf0)}; HKMIN=1'
        calls += '\ncall KJPFILL(js,js,kgrid,1d0,1,2)\ncall KJPFILL(js,js,kgrid,0d0,1,3)'
    singular = 'SLAM0=0; SF0=0; LAM=%.17e; B0K=%.17e' % (LAM, B0K)
    return f'''program test_ntv
use globalm
use dimensim
use rcomdm
use kineticm
use torquem
implicit none
complex*16, allocatable :: a(:,:,:),b(:,:,:),c(:,:,:),dm(:,:,:),
& e(:,:,:),f(:,:,:),g(:,:,:),h(:,:,:)
integer js,kgrid,is,KPART
real*8 RLMV
NR=2; NRP1=3; MSMAX=1; NCHI=4; NFIT=0; NSPECIES=2; MLMAX=1
NXCOMP=3; NYCOMP=7; MXMAX=NXCOMP*MSMAX; MYMAX=NYCOMP*MSMAX
KXX1=1; KXB1=2; KXV1=3; KYX2=1; KYB2=2; KYB3=3; KYV2=4; KYPR=5
KYPPARA=6; KYPPERP=7
RNTOR={RNTOR}d0; KNTV=21; INCKIN=1; KEFORM=2; IFOWP=0; IFOWT=0
ODWKCOM=.true.; ISWEEP=1; NSWEEP=1; ALPHAD=1; PTRAPH={PTRAPH}d0; CTEDGE=10
INCDPHI=0; IDIAMB=0; IDIAMTE=0
PSPECIES_AP=0; PSPECIES_AT=0; PSPECIES_NP=0; PSPECIES_NTB=0; PSPECIES_NTD=0
allocate(HKMIN(NRP1,2),HKMAX(NRP1,2),WFUN(NRP1,2),RLM(1))
HKMIN=0; HKMAX=1; WFUN=1
allocate(SF0(1,NSPECIES,0:3,4),SLAM0(1,NSPECIES))
allocate(SVPARA0(1,1,NSPECIES),SVPERP0(1,1,NSPECIES),SVDPHI0(1,1,NSPECIES))
allocate(SVX10(1,1,NSPECIES),SVX20(1,1,NSPECIES),SVQ10(1,1,NSPECIES))
allocate(SVQ20(1,1,NSPECIES),SVQ30(1,1,NSPECIES),SVDP0(1,1,NSPECIES))
SVPARA0={fortran_complex(GS)}; SVPERP0=0; SVDPHI0=0; SVDP0=0
{singular}
{species}
RLM=RLMV
OMEGAB={OMEGAB}d0
call set_reversal_frequency({OMEGAB0}d0,{OMEGABS}d0)
allocate(CS(NRP1),CSM(NR),CSH(NR),RM(MSMAX,2),RJAM(NRP1,NCHI))
allocate(DPSIDS(NRP1),DPSIDSM(NR))
CS={real_array(CS)}; CSM={real_array(CSM)}; CSH=[.4d0,.6d0]; RM=2; RJAM=1
DPSIDS={real_array(DPSIDS)}; DPSIDSM={real_array(DPSIDSM)}
allocate(X1U(NRP1,1),X2U(NRP1,1),B1U(NRP1,1),B2U(NRP1,1),B3U(NRP1,1))
X1U(:,1)={array(X1)}; B1U(:,1)={array(B1)}
X2U=0; B2U=0; B3U=0
X2U(1:2,1)={array(X2)}; B2U(1:2,1)={array(B2)}; B3U(1:2,1)={array(B3)}
allocate(TORQUENTV(NRP1),TORQUENTVI(NRP1),TORQUENTVE(NRP1))
allocate(VI(4,1,NSPECIES),VPARA(1,1),VPERP(1,1),VDPHI(1,1))
allocate(VX1(1,1),VX2(1,1),VQ1(1,1),VQ2(1,1),VQ3(1,1),VDP(1,1))
allocate(VPARA0(1,1),VPERP0(1,1),VDPHI0(1,1))
allocate(VX10(1,1),VX20(1,1),VQ10(1,1),VQ20(1,1),VQ30(1,1),VDP0(1,1))
VI=0; VI(1,1,1)={fortran_complex(kernel)}
VI0=0; VI0(1,1)={vi0}
VPARA={fortran_complex(G)}; VPERP=0; VDPHI=0; VDP=0
VPARA0={fortran_complex(G0)}; VPERP0=0; VDPHI0=0; VDP0=0
call alloc3(VX1PARA); call alloc3(VX1PERP); call alloc3(VX2PARA); call alloc3(VX2PERP)
call alloc3(VQ1PARA); call alloc3(VQ1PERP); call alloc3(VQ2PARA); call alloc3(VQ2PERP)
call alloc3(VQ3PARA); call alloc3(VQ3PERP); call alloc3(VDPPARA); call alloc3(VDPPERP)
call alloc3(VX1PARAM); call alloc3(VX1PERPM); call alloc3(VX2PARAM); call alloc3(VX2PERPM)
call alloc3(VQ1PARAM); call alloc3(VQ1PERPM); call alloc3(VQ2PARAM); call alloc3(VQ2PERPM)
call alloc3(VQ3PARAM); call alloc3(VQ3PERPM); call alloc3(VDPPARAM); call alloc3(VDPPERPM)
call ALLOCATEDWKCOMPMAT
call execute_command_line('mkdir -p DATA_STORAGE')
{fills}
call DEALLOCATEDWKCOMPMAT
allocate(a(MXMAX,MXMAX,NRP1),b(MXMAX,MXMAX,NRP1),c(MXMAX,MXMAX,NRP1))
allocate(dm(MYMAX,MYMAX,NRP1),e(MXMAX,MYMAX,NRP1),f(MYMAX,MXMAX,NRP1))
allocate(g(MYMAX,MXMAX,NRP1),h(MXMAX,MYMAX,NRP1))
a=0; b=0; c=0; dm=0; e=0; f=0; g=0; h=0
dm(KYV2,KYPPARA,1:2)={array(rows)}
! Node-1 radial rows are never contracted (the CALCDWKPROF node row uses I=IS+1).
h(KXV1,KYPPARA,1)=1; h(KXV1,KYPPERP,1)=1
call CALCDWKCOMP(a,b,c,dm,e,f,g,h)
do is=1,NR
 write(*,'(A,I2,ES26.17)') 'TORQUE ',is,TORQUENTVI(is)
enddo
contains
subroutine fill_surface(js,kgrid)
integer js,kgrid
VX1PARAC=0; VX1PERPC=0; VX2PARAC=0; VX2PERPC=0; VQ1PARAC=0
VQ1PERPC=0; VQ2PARAC=0; VQ2PERPC=0; VQ3PARAC=0; VQ3PERPC=0
call clear_hermitian()
{calls}
call WRITE_SURFACE_QUANTITIES(js,kgrid)
end subroutine
subroutine alloc3(x)
complex*16, allocatable :: x(:,:,:)
allocate(x(MSMAX,MSMAX,NRP1)); x=0
end subroutine
end program
'''


# Zeroing the same-drive blocks per surface is KJPCOEFF's job; the extracted
# driver does it through this hook, which is empty for sources without them.
CLEAR_WITH = '''subroutine clear_hermitian()
use kineticm
if (allocated(HWKC)) HWKC = 0
end subroutine
subroutine set_reversal_frequency(x,y)
use kineticm
real*8 x,y
OMEGAB0 = x
allocate(SOMEGAB0(1,NSPECIES))
SOMEGAB0 = y
end subroutine
'''
CLEAR_WITHOUT = '''subroutine clear_hermitian()
end subroutine
subroutine set_reversal_frequency(x,y)
real*8 x,y
end subroutine
'''


def extract(text, name):
    match = re.search(r'^      (?:COMPLEX\*16 )?(?:SUBROUTINE|FUNCTION) %s\b' % name,
                      text, re.M)
    if match is None:
        return ''
    end = re.search(r'^      END(\s+(SUBROUTINE|FUNCTION)\s+%s)?\s*$' % name,
                    text[match.start():], re.M)
    return text[match.start():match.start() + end.end()] + '\n'


def run_native(source, build_parent, kernel, eps, kernel0=None, sf0=None):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'kinetic.f').read_text()
    module = (src / 'kineticm.f').read_text()
    with tempfile.TemporaryDirectory(prefix='ntv-herm-', dir=build_parent) as directory:
        work = Path(directory)
        (work / 'routines.f').write_text(''.join(extract(text, n) for n in ROUTINES))
        clear = CLEAR_WITH if 'HWKC' in module else CLEAR_WITHOUT
        (work / 'clear.f90').write_text(clear)
        (work / 'driver.f90').write_text(driver(kernel, eps, kernel0, sf0).replace('\n& ', ' &\n& '))
        names = ['dimensim.f', 'globalm.f', 'rcomdm.f', 'kineticm.f', 'anisotropicm.f',
                 'gijlm.f', 'convolcofm.f', 'torquem.f', 'ToolBox.f90']
        command = ['gfortran', '-O0', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-fallow-argument-mismatch', '-ffunction-sections', '-fdata-sections',
                   '-ffpe-trap=zero,invalid', '-fcheck=bounds', '-w', '-I' + str(src)]
        command += [str(src / n) for n in names]
        command += ['routines.f', 'clear.f90', 'driver.f90', '-Wl,--gc-sections',
                    '-o', 'test']
        compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if compiled.returncode:
            raise RuntimeError(compiled.stdout + compiled.stderr)
        result = subprocess.run([str(work / 'test')], cwd=work, capture_output=True,
                                text=True)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        torque = [float(line.split()[2]) for line in result.stdout.splitlines()
                  if line.split()[:1] == ['TORQUE']]
        if len(torque) != 2:
            raise RuntimeError('missing torque rows:\n' + result.stdout[-2000:])
        forms = None
        path = work / 'NTV_TORQUE_FORMS.OUT'
        if path.exists():
            forms = [[float(x) for x in line.split()[2:4]]
                     for line in path.read_text().splitlines()
                     if line.strip() and not line.startswith('%')]
        return torque, forms


class TestHermitianNTVTorque(unittest.TestCase):
    source = Path(__file__).resolve().parents[1]
    build_parent = Path.cwd()
    _coefficient = None

    def run_case(self, kernel, eps, kernel0=None, sf0=None):
        return run_native(self.source, self.build_parent, kernel, eps, kernel0, sf0)

    def coefficient(self):
        """c in T = c Im(dW_k), from a field that satisfies ideal induction,
        where the weak and the same-drive work coincide."""
        if TestHermitianNTVTorque._coefficient is None:
            kernel = 0.7 + 0.2j
            torque, _ = self.run_case(kernel, 0.0)
            ratios = [t / e for t, e in zip(torque, expected_torque(kernel))]
            self.assertAlmostEqual(ratios[1] / ratios[0], 1.0, places=12)
            self.assertTrue(any(abs(ratios[0] / c - 1.0) < 1e-12 for c in COEFFICIENTS),
                            'unexpected torque coefficient %r' % ratios[0])
            TestHermitianNTVTorque._coefficient = ratios[0]
        return TestHermitianNTVTorque._coefficient

    def test_consistent_drive_reproduces_same_drive_torque(self):
        # Control: passes with the weak contraction as well.
        c = self.coefficient()
        kernel = 0.3 + 0.5j
        torque, _ = self.run_case(kernel, 0.0)
        for actual, expected in zip(torque, expected_torque(kernel)):
            self.assertNotEqual(expected, 0.0)
            self.assertAlmostEqual(actual / (c * expected), 1.0, places=12)

    def test_reactive_kernel_gives_zero_torque(self):
        scale_ = max(abs(t) for t in self.run_case(0.7j, 0.0)[0])
        for eps in (0.0, 0.05):
            torque, _ = self.run_case(0.7 + 0.0j, eps)
            for actual in torque:
                self.assertLessEqual(abs(actual), 1e-13 * scale_)

    def test_torque_invariant_under_inconsistent_test_field(self):
        kernel = 0.7 + 0.2j
        reference, _ = self.run_case(kernel, 0.0)
        shifts = {}
        for eps in (0.02, 0.04):
            torque, forms = self.run_case(kernel, eps)
            for actual, ref in zip(torque, reference):
                self.assertAlmostEqual(actual / ref, 1.0, places=12)
            # The weak form is still written, and it moves linearly with the
            # mismatch of the test field.
            self.assertIsNotNone(forms, 'NTV_TORQUE_FORMS.OUT not written')
            shifts[eps] = [tw - th for th, tw in forms]
            for (th, tw), actual, ref in zip(forms, torque, reference):
                self.assertAlmostEqual(th / actual, 1.0, places=12)
                self.assertGreater(abs(tw - th), 1e-4 * abs(ref))
        for a, b in zip(shifts[0.04], shifts[0.02]):
            self.assertAlmostEqual(a / b, 2.0, places=6)

    def test_trapped_drift_reversal_uses_bounce_frequency_at_each_pitch(self):
        # Reactive regular kernel, dissipative reversal kernel: the torque is
        # the reversal term alone, weighted with 1/(2 omega_b) at that pitch.
        c = self.coefficient()
        cases = [(0.7 + 0.0j, 0.3 + 0.0j), (0.7 + 0.0j, 0.3 + 0.25j),
                 (0.7 + 0.2j, 0.3 + 0.25j)]
        scale_ = max(abs(c * t) for t in expected_reversal_torque(0.7j, 0.25j))
        for kernel, kernel0 in cases:
            torque, _ = self.run_case(kernel, 0.05, kernel0)
            for actual, ref in zip(torque, expected_reversal_torque(kernel, kernel0)):
                self.assertLessEqual(abs(actual - c * ref), 1e-12 * scale_)

    def test_passing_log_resonance_uses_bounce_frequency_at_each_pitch(self):
        c = self.coefficient()
        cases = [(0.7 + 0.0j, 0.45 + 0.0j), (0.7 + 0.0j, SF0V), (0.7 + 0.2j, SF0V)]
        scale_ = max(abs(c * t) for t in expected_log_torque(0.2j, 0.35j))
        for kernel, sf0 in cases:
            torque, _ = self.run_case(kernel, 0.05, sf0=sf0)
            for actual, ref in zip(torque, expected_log_torque(kernel, sf0)):
                self.assertLessEqual(abs(actual - c * ref), 1e-12 * scale_)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build-parent', type=Path, default=Path.cwd())
    args = parser.parse_args()
    TestHermitianNTVTorque.source = args.source
    TestHermitianNTVTorque.build_parent = args.build_parent
    unittest.main(argv=['test'])


if __name__ == '__main__':
    main()
