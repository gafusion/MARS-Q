"""Hamada-gauge |B| spectrum of NTV_pre::cal_xi_liu.

The unmodified production routines of ntv_pre.f, ntv_lib.f and splines.f are
compiled with a fixture displacement (get_dimensV/get_disp).

(a) A purely field-aligned MHD displacement X3 B contributes nothing, and X3
    never changes the spectrum.
(b) On a Hamada base the spectrum equals an independent Python evaluation of
    the Hamada-gauge kernel with dB_L = b.deltaB + xi_perp.grad B and source
    d = -dB_L/B - kappa.xi_perp, including the Park regularisation. Its values
    for CASE_B below are stored in fixtures/ntv_hamada_rmp_torque.json.
(c) Circular large-aspect geometry through the full chain cal_metric,
    cal_divmetric, coord_trans, coord_map, cal_xi_liu with the ideal field
    curl(xi_perp x B) against an independent oracle: exact Hamada chart,
    Cartesian finite differences of xi_perp x B, |B| and B.grad|B|.
(d) The spectrum of a surface uses no displacement of other surfaces.
"""
import cmath
import json
import math
import re

from support import ROOT, NativeTest

SRC = ROOT / "MarsQ_2FK"
FIX = ROOT / "tests" / "fixtures" / "ntv_hamada_rmp_torque.json"


def routine(name, file="ntv_pre.f"):
    source = (SRC / file).read_text()
    match = re.search(r"^ +subroutine +%s\b.*?^ +end subroutine +%s\b" % (name, name),
                      source, re.M | re.S | re.I)
    if match is None:
        raise ValueError("Native routine not found: " + name)
    return match.group(0)


CASE_B = dict(n=2, iota=0.37, dqds=0.4, bth=0.8, b00=1.3, eps=0.2, c0=0.7, c1=-0.25,
              cp=0.3, fr=0.6, fi=-0.35, s=[0.0, 0.25, 0.5, 0.75, 1.0], surface=3, n_chi=128, m_max=8,
              modes=[-3, -2, -1, 0, 1, 2, 3],
              a=[[0.11, -0.05], [-0.31, 0.2], [0.52, 0.13], [0.07, -0.4],
                 [0.45, 0.33], [-0.18, 0.09], [0.06, 0.02]],
              b=[[0.2, 0.1], [0.0, -0.15], [0.3, 0.0], [-0.1, 0.05],
                 [0.12, -0.22], [0.04, 0.3], [-0.03, 0.01]])

CASE_C = dict(n=2, eps=0.1, np=101, n_chi=256, m_max=10, modes=[-3, -2, -1, 0, 1],
              x1=[[0.3, 0.1], [1.0, -0.2], [0.5, 0.4], [-0.2, 0.1], [0.1, 0.05]],
              x2=[[0.05, 0.0], [0.2, 0.1], [-0.1, 0.3], [0.0, 0.1], [0.07, -0.02]],
              surfaces=[31, 51, 63])


def q_c(s):
    return 1.1 + 0.8 * s * s


def x1_c(m, s):
    return s * (1.0 - s) ** 2


def x2_c(m, s):
    return s * (1.0 - s)


def fortran_list(values):
    return "[" + ", ".join("%.17gd0" % v for v in values) + "]"


FIXTURE = """      module NTV_pre
      use NTV_lib
      implicit none
      real*8,parameter::npi=3.1415926535897932
      integer :: nfix = 0, mfix = 1, npfix = 1
      integer, allocatable :: mlist(:)
      real, allocatable :: x1r_(:,:),x1i_(:,:),x2r_(:,:),x2i_(:,:)
      real, allocatable :: x3r_(:,:),x3i_(:,:)
      real, allocatable :: b1r_(:,:),b1i_(:,:),b2r_(:,:),b2i_(:,:)
      real, allocatable :: b3r_(:,:),b3i_(:,:)
      contains
      subroutine get_dimensV(MM1,NP1,n)
      integer MM1,NP1,n
      MM1 = mfix
      NP1 = npfix
      n = nfix
      end subroutine get_dimensV
      subroutine get_disp(M10,NP1,M,X1R,X1I,X2R,X2I,X3R,X3I)
      integer M10,NP1,M(M10)
      real X1R(NP1,M10),X1I(NP1,M10),X2R(NP1,M10),X2I(NP1,M10)
      real X3R(NP1,M10),X3I(NP1,M10)
      M = mlist
      X1R = x1r_; X1I = x1i_; X2R = x2r_; X2I = x2i_
      X3R = x3r_; X3I = x3i_
      end subroutine get_disp
      subroutine get_field(M10,NP1,B1R,B1I,B2R,B2I,B3R,B3I)
      integer M10,NP1
      real B1R(NP1,M10),B1I(NP1,M10),B2R(NP1,M10),B2I(NP1,M10)
      real B3R(NP1,M10),B3I(NP1,M10)
      B1R = b1r_; B1I = b1i_; B2R = b2r_; B2I = b2i_
      B3R = b3r_; B3I = b3i_
      end subroutine get_field
%s
      end module NTV_pre
"""

ORACLE_B = """program hamada_base
use NTV_pre
implicit none
integer, parameter :: np = 5, n_chi = %(n_chi)d, m_max = %(m_max)d, nm = %(nm)d
real :: chiH(np, n_chi), dphiH(np, n_chi), dBav(np), gg(np, n_chi, 3)
real :: dBds(np, n_chi), dBdc(np, n_chi), omega(np), dwds(np)
real :: s(np), q(np), B(np, n_chi), dm(np, n_chi, 7), bth(np), cp(np), th
real :: a(2, nm), bb(2, nm)
complex :: bmn(2*m_max+1, 1, np)
character(len=8) :: mode
integer :: i, j
call get_command_argument(1, mode)
nfix = %(n)d; mfix = nm; npfix = np
allocate(mlist(nm), x1r_(np,nm), x1i_(np,nm), x2r_(np,nm), x2i_(np,nm))
allocate(x3r_(np,nm), x3i_(np,nm))
allocate(b1r_(np,nm), b1i_(np,nm), b2r_(np,nm), b2i_(np,nm), b3r_(np,nm), b3i_(np,nm))
b1r_ = 0; b1i_ = 0; b3r_ = 0; b3i_ = 0
mlist = %(modes)s
a = reshape(%(a)s, [2, nm]); bb = reshape(%(b)s, [2, nm])
s = %(s)s
q = 1d0/%(iota).17gd0 + %(dqds).17gd0*(s - s(%(surface)d))
do i = 1, np
   x1r_(i,:) = a(1,:)*s(i) + bb(1,:)
   x1i_(i,:) = a(2,:)*s(i) + bb(2,:)
   x3r_(i,:) = 0.7d0*bb(2,:) - 0.2d0
   x3i_(i,:) = 0.3d0*a(1,:) + 0.1d0
end do
x2r_ = 0.9d0*x1i_; x2i_ = -0.4d0*x1r_
do i = 1, np
   b2r_(i,:) = %(fr).17gd0*a(2,:)*s(i) - bb(1,:)
   b2i_(i,:) = %(fi).17gd0*a(1,:) + bb(2,:)*s(i)
end do
if (mode == 'a') then
   x1r_ = 0; x1i_ = 0; x2r_ = 0; x2i_ = 0; b2r_ = 0; b2i_ = 0
else if (mode == 'b') then
   x3r_ = 0; x3i_ = 0
end if
dphiH = 0; gg = 0; dm = 0; omega = 0; dwds = 0
do j = 1, n_chi
   th = (j - 0.5d0)*2*npi/n_chi
   chiH(:, j) = th
   B(:, j) = %(b00).17gd0*(1 - %(eps).17gd0*cos(th))
   dBdc(:, j) = %(b00).17gd0*%(eps).17gd0*sin(th)
   dBds(:, j) = %(c0).17gd0 + %(c1).17gd0*cos(th)
end do
gg(:, :, 3) = %(bth).17gd0
dm(:, :, 7) = 1.7d0
dm(:, :, 6) = -1.7d0**2
bth = %(bth).17gd0; cp = %(cp).17gd0
call cal_xi_liu(np, n_chi, chiH, 1, dphiH, dBav, gg, dBds, dBdc, &
                omega, dwds, m_max, 1, bmn, s, q, B, dm, bth, cp)
do j = 1, 2*m_max+1
   print '(i4,2es26.17)', j-m_max-1, bmn(j, 1, %(surface)d)
end do
end program
"""

ORACLE_C = """program circular_chain
use NTV_pre
implicit none
integer, parameter :: np = %(np)d, n_chi = %(n_chi)d, m_max = %(m_max)d, nm = %(nm)d
real, parameter :: eps = %(eps).17gd0
real :: chi(n_chi), R(np, n_chi), Z(np, n_chi), s(np), q(np), F(np), dpsids(np)
real :: Jac(np, n_chi), gg(np, n_chi, 3), B(np, n_chi), dBds(np, n_chi)
real :: dBdc(np, n_chi), chiH(np, n_chi), dphiH(np, n_chi), dm(np, n_chi, 7)
real :: tmp(np, n_chi), bth(np), cp(np), dBav(np), omega(np), dwds(np)
real :: c1(2, nm), c2(2, nm)
complex :: bmn(2*m_max+1, 1, np)
integer :: i, j, k, ios, surf(%(ns)d)
real :: fb(6)
character(len=8) :: mode
call get_command_argument(1, mode)
nfix = %(n)d; mfix = nm; npfix = np
allocate(mlist(nm), x1r_(np,nm), x1i_(np,nm), x2r_(np,nm), x2i_(np,nm))
allocate(x3r_(np,nm), x3i_(np,nm))
allocate(b1r_(np,nm), b1i_(np,nm), b2r_(np,nm), b2i_(np,nm), b3r_(np,nm), b3i_(np,nm))
b1r_ = 0; b1i_ = 0; b2r_ = 0; b2i_ = 0; b3r_ = 0; b3i_ = 0
mlist = %(modes)s
c1 = reshape(%(x1)s, [2, nm]); c2 = reshape(%(x2)s, [2, nm])
surf = %(surfaces)s
open(10, file='field.dat', status='old')
do
   read(10, *, iostat=ios) i, k, fb
   if (ios /= 0) exit
   b1r_(i,k) = fb(1); b1i_(i,k) = fb(2); b2r_(i,k) = fb(3)
   b2i_(i,k) = fb(4); b3r_(i,k) = fb(5); b3i_(i,k) = fb(6)
end do
close(10)
call xgrid(n_chi, 0d0, 2*npi, chi)
do i = 1, np
   s(i) = (i - 1d0)/(np - 1d0)
   q(i) = 1.1d0 + 0.8d0*s(i)**2
   F(i) = 1
   dpsids(i) = eps**2*s(i)/(q(i)*sqrt(1 - (eps*s(i))**2))
   R(i, :) = 1 + eps*s(i)*cos(chi)
   Z(i, :) = eps*s(i)*sin(chi)
   x1r_(i,:) = c1(1,:)*s(i)*(1 - s(i))**2
   x1i_(i,:) = c1(2,:)*s(i)*(1 - s(i))**2
   x2r_(i,:) = c2(1,:)*s(i)*(1 - s(i))
   x2i_(i,:) = c2(2,:)*s(i)*(1 - s(i))
end do
x3r_ = 0.5d0; x3i_ = -0.25d0
if (mode == 'jump') then
   ! neighbours of the tested surfaces carry an unrelated displacement,
   ! as the undefined half-mesh slot X2U(NRP1) does at the boundary
   do k = 1, size(surf)
      do i = surf(k) - 1, surf(k) + 1, 2
         x1r_(i,:) = 3*x1r_(i,:) + 0.1d0; x2r_(i,:) = 0; x2i_(i,:) = 0
      end do
   end do
end if
call cal_metric(np, n_chi, s, chi, R, Z, F, dpsids, q, Jac, gg, B, dBds, dBdc)
call cal_divmetric(np, n_chi, s, chi, R, Z, F, dpsids, Jac, B, dm)
bth = 0
do i = 2, np
   bth(i) = n_chi*dpsids(i)/sum(Jac(i, :))
end do
cp = 0
call coord_trans(np, n_chi, chi, R, q, Jac, chiH, dphiH)
call coord_map(n_chi, chi, np, chiH, B)
call coord_map(n_chi, chi, np, chiH, dBds)
call coord_map(n_chi, chi, np, chiH, dBdc)
do k = 1, 3
   tmp = gg(:, :, k)
   call coord_map(n_chi, chi, np, chiH, tmp)
   gg(:, :, k) = tmp
end do
do k = 1, 7
   tmp = dm(:, :, k)
   call coord_map(n_chi, chi, np, chiH, tmp)
   dm(:, :, k) = tmp
end do
omega = 0; dwds = 0
call cal_xi_liu(np, n_chi, chiH, 1, dphiH, dBav, gg, dBds, dBdc, &
                omega, dwds, m_max, 1, bmn, s, q, B, dm, bth, cp)
do k = 1, size(surf)
   do j = 1, 2*m_max+1
      print '(2i4,2es26.17)', surf(k), j-m_max-1, bmn(j, 1, surf(k))
   end do
end do
end program
"""


def parse(stdout, keys):
    out = {}
    for line in stdout.split("\n"):
        parts = line.split()
        if len(parts) == keys + 2:
            out[tuple(int(p) for p in parts[:keys])] = complex(float(parts[-2]),
                                                               float(parts[-1]))
    return out


# ---------------------------------------------------------------- oracle (c)
def circular_case(case, surface):
    """Explicit circular equilibrium and displacement of case (c) in Cartesian form."""
    eps, n, nchi = case["eps"], case["n"], case["n_chi"]
    s0 = (surface - 1.0) / (case["np"] - 1.0)
    modes = case["modes"]
    c1 = [complex(*c) for c in case["x1"]]
    c2 = [complex(*c) for c in case["x2"]]

    def psip(s):
        return eps * eps * s / (q_c(s) * math.sqrt(1 - (eps * s) ** 2))

    def frame(s, chi, phi):
        R = 1 + eps * s * math.cos(chi)
        es = (eps * math.cos(chi) * math.cos(phi), -eps * math.cos(chi) * math.sin(phi),
              eps * math.sin(chi))
        ec = (-eps * s * math.sin(chi) * math.cos(phi), eps * s * math.sin(chi) * math.sin(phi),
              eps * s * math.cos(chi))
        ep = (-R * math.sin(phi), -R * math.cos(phi), 0.0)
        J = dot(es, cross(ec, ep))
        gs = tuple(c / J for c in cross(ec, ep))
        gp = tuple(c / J for c in cross(es, ec))
        Bv = tuple(psip(s) * a + b for a, b in zip(cross(gp, gs), gp))   # F = 1
        return es, gs, Bv, J

    def coords(x):
        R = math.hypot(x[0], x[1])
        return (math.hypot(R - 1, x[2]) / eps, math.atan2(x[2], R - 1),
                math.atan2(-x[1], x[0]))

    def xi_perp(x):
        s, chi, phi = coords(x)
        es, gs, Bv, J = frame(s, chi, phi)
        ph = cmath.exp(1j * n * phi)
        X1 = sum(c * x1_c(m, s) * cmath.exp(1j * m * chi) for m, c in zip(modes, c1)) * ph
        X2 = sum(c * x2_c(m, s) * cmath.exp(1j * m * chi) for m, c in zip(modes, c2)) * ph
        B2 = dot(Bv, Bv)
        Bs = dot(Bv, es)
        bx = cross(Bv, gs)
        return tuple(X1 * (es[i] - Bs * Bv[i] / B2) + X2 * J * bx[i] / B2 for i in range(3))

    def bmod(x):
        s, chi, phi = coords(x)
        Bv = frame(s, chi, phi)[2]
        return math.sqrt(dot(Bv, Bv))

    def position(s, chi, phi):
        R = 1 + eps * s * math.cos(chi)
        return (R * math.cos(phi), -R * math.sin(phi), eps * s * math.sin(chi))

    return eps, n, nchi, s0, frame, coords, xi_perp, bmod, position


def circular_field(case, surface, modes):
    """MARS harmonics of J deltaB^(s,chi,phi) of the ideal field curl(xi_perp x B)."""
    eps, n, nchi, s0, frame, coords, xi_perp, bmod, position = circular_case(case, surface)
    h = 1e-5 * eps

    def potential(x):
        return cross(xi_perp(x), frame(*coords(x))[2])

    npts = 128
    values = []
    for j in range(npts):
        chi = 2 * math.pi * j / npts
        x = position(s0, chi, 0.0)
        d = []
        for i in range(3):
            xp = list(x); xm = list(x)
            xp[i] += h; xm[i] -= h
            d.append([(u - v) / (2 * h) for u, v in zip(potential(xp), potential(xm))])
        curl = (d[1][2] - d[2][1], d[2][0] - d[0][2], d[0][1] - d[1][0])
        R = 1 + eps * s0 * math.cos(chi)
        es, gs, Bv, J = frame(s0, chi, 0.0)
        ec = (-eps * s0 * math.sin(chi), 0.0, eps * s0 * math.cos(chi))
        ep = (0.0, -R, 0.0)
        gc = tuple(c / J for c in cross(ep, es))
        gp = tuple(c / J for c in cross(es, ec))
        jac = R * eps * eps * s0                    # MARS Jacobian of cal_metric
        values.append([jac * dot(curl, g) for g in (gs, gc, gp)])
    return {m: [sum(values[j][i] * cmath.exp(-1j * m * 2 * math.pi * j / npts)
                    for j in range(npts)) / npts for i in range(3)] for m in modes}


def circular_oracle(case, surface):
    """Independent |B| perturbation in perturbed Hamada coordinates."""
    eps, n, nchi, s0, frame, coords, xi_perp, bmod, position = circular_case(case, surface)

    h = 1e-5 * eps
    e = eps * s0
    q0 = q_c(s0)
    bth = 1 / (q0 * math.sqrt(1 - e * e))
    sigma = abs(n * 1.6 * s0) * (2.0 / (case["np"] - 1)) * 0.5
    d, dL, bgb, xs = [], [], [], []

    def potential(x):
        return cross(xi_perp(x), frame(*coords(x))[2])

    for j in range(nchi):
        th = (j + 0.5) * 2 * math.pi / nchi
        chi = th
        for _ in range(50):
            chi -= (chi + e * math.sin(chi) - th) / (1 + e * math.cos(chi))
        t = math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(chi / 2))
        t += math.pi * math.floor((chi + math.pi) / (2 * math.pi))
        dphi = q0 * (chi + e * math.sin(chi) - 2 * t)
        x = position(s0, chi, -dphi)
        xs.append(x)
        grad, dA = [], []
        for i in range(3):
            xp = list(x); xm = list(x)
            xp[i] += h; xm[i] -= h
            dA.append([(u - v) / (2 * h) for u, v in zip(potential(xp), potential(xm))])
            grad.append((bmod(xp) - bmod(xm)) / (2 * h))
        curl = (dA[1][2] - dA[2][1], dA[2][0] - dA[0][2], dA[0][1] - dA[1][0])
        Bv = frame(*coords(x))[2]
        xv = xi_perp(x)
        # Park 2009 Eq. 7 with the ideal field, and the Hamada source
        # d = -dB_L/B - kappa.xi_perp with kappa.xi_perp = xi_perp.grad B/B (cp = 0)
        dBL = dot(Bv, curl) / bmod(x) + dot(xv, grad)
        d.append(-dBL / bmod(x) - dot(xv, grad) / bmod(x))
        dL.append(dBL)
        bgb.append(dot(Bv, grad))
    kmin = -(nchi // 2)
    lam, lam0 = [0j] * nchi, {}
    for kk in range(nchi):
        k = kmin + kk
        dk = sum(d[j] * cmath.exp(-1j * k * (j + 0.5) * 2 * math.pi / nchi)
                 for j in range(nchi)) / nchi
        xk = k + n * q0
        lk = 1j * dk * xk / (bth * (xk * xk + sigma * sigma))
        lam0[k] = 1j * dk / (bth * xk)
        for j in range(nchi):
            lam[j] += lk * cmath.exp(1j * k * (j + 0.5) * 2 * math.pi / nchi)

    def lam_unregularised(x):
        """sigma = 0 solution as a function of position, via the exact chart."""
        s, chi, phi = coords(x)
        th = chi + e * math.sin(chi)
        t = math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(chi / 2))
        t += math.pi * math.floor((chi + math.pi) / (2 * math.pi))
        zeta = phi + q0 * (chi + e * math.sin(chi) - 2 * t)
        return sum(v * cmath.exp(1j * (k * th + n * zeta)) for k, v in lam0.items())

    # Independent check of the chart: B.grad lam = -div xi_perp along the
    # Cartesian field line (s is a flux label, so the derivative stays on s0).
    for j in range(0, nchi, 37):
        x = xs[j]
        Bv = frame(*coords(x))[2]
        hb = 1e-4 * eps / math.sqrt(dot(Bv, Bv))
        bgl = (lam_unregularised([x[i] + hb * Bv[i] for i in range(3)])
               - lam_unregularised([x[i] - hb * Bv[i] for i in range(3)])) / (2 * hb)
        assert abs(bgl + d[j]) < 1e-5 * max(abs(v) for v in d), (j, bgl, d[j])
    dB = [dL[j] + lam[j] * bgb[j] for j in range(nchi)]
    return {k: sum(dB[j] * cmath.exp(-1j * k * (j + 0.5) * 2 * math.pi / nchi)
                   for j in range(nchi)) / nchi
            for k in range(-case["m_max"], case["m_max"] + 1)}


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


class NtvHamadaGauge(NativeTest):
    def build(self, oracle):
        lib = (SRC / "ntv_lib.f").read_text()
        pre = FIXTURE % "\n".join(routine(r) for r in (
            "cal_metric", "cal_divmetric", "coord_trans", "coord_map", "cal_xi_liu"))
        return self.compile({"splines.f": (SRC / "splines.f").read_text(),
                             "ntv_lib.f": lib, "ntv_pre.f": pre, "oracle.f90": oracle})

    def hamada_base(self):
        c = CASE_B
        flat = lambda rows: fortran_list([v for row in rows for v in row])
        return self.build(ORACLE_B % dict(
            c, nm=len(c["modes"]), modes="[" + ", ".join(map(str, c["modes"])) + "]",
            a=flat(c["a"]), b=flat(c["b"]), s=fortran_list(c["s"])))

    def test_field_aligned_x3_contributes_nothing(self):
        exe = self.hamada_base()
        pure = parse(self.execute(exe, "a").stdout, 1)
        self.assertEqual(len(pure), 2 * CASE_B["m_max"] + 1)
        self.assertEqual(max(abs(v) for v in pure.values()), 0.0)
        without = parse(self.execute(exe, "b").stdout, 1)
        with_x3 = parse(self.execute(exe, "b3").stdout, 1)
        self.assertEqual(without, with_x3)

    def test_hamada_base_matches_rmp_torque_hamada_path(self):
        ref = json.loads(FIX.read_text())
        self.assertEqual(ref["case"], json.loads(json.dumps(CASE_B)))
        got = parse(self.execute(self.hamada_base(), "b").stdout, 1)
        expect = {int(k): complex(*v) for k, v in ref["delta_b_hamada"].items()}
        scale = max(abs(v) for v in expect.values())
        for k, v in expect.items():
            self.assertLess(abs(got[(k,)] - v), 1e-11 * scale, k)

    def circular_chain(self):
        c = CASE_C
        flat = lambda rows: fortran_list([v for row in rows for v in row])
        # the field couples xi's modes to toroidal sidebands up to eps**|dm|
        modes = list(range(-14, 13))
        pad = lambda key: [c[key][c["modes"].index(m)] if m in c["modes"] else [0.0, 0.0]
                           for m in modes]
        with open(self.work / "field.dat", "w") as out:
            for surface in c["surfaces"]:
                field = circular_field(c, surface, modes)
                for k, m in enumerate(modes, 1):
                    out.write("%d %d %s\n" % (surface, k, " ".join(
                        "%.17e %.17e" % (v.real, v.imag) for v in field[m])))
        exe = self.build(ORACLE_C % dict(
            c, nm=len(modes), ns=len(c["surfaces"]),
            modes="[" + ", ".join(map(str, modes)) + "]",
            x1=flat(pad("x1")), x2=flat(pad("x2")),
            surfaces="[" + ", ".join(map(str, c["surfaces"])) + "]"))
        return exe

    def test_circular_chain_matches_cartesian_oracle(self):
        c = CASE_C
        got = parse(self.execute(self.circular_chain()).stdout, 2)
        for surface in c["surfaces"]:
            expect = circular_oracle(c, surface)
            scale = max(abs(v) for v in expect.values())
            err = max(abs(got[(surface, k)] - v) for k, v in expect.items())
            print("surface", surface, "max|bmn|", scale, "max error", err)
            self.assertLess(err, 5e-3 * scale, surface)

    def test_spectrum_is_local_to_the_surface(self):
        exe = self.circular_chain()
        self.assertEqual(parse(self.execute(exe).stdout, 2),
                         parse(self.execute(exe, "jump").stdout, 2))
