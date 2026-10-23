"""ALDLT final pivot: singularity is judged against the pivot's own scale."""
from support import NativeTest, NATIVE, routine

PRECISION = "module prec_const\ninteger,parameter::rkind=kind(1.d0)\nend module\n"

ORACLE = """program final_pivot
use prec_const, only: rkind
implicit none
real(rkind) :: band(2, 3), singular(2, 2), rhs(3), dense(3, 3)
real(rkind), parameter :: exact(3) = [1._rkind, -2._rkind, 3._rkind]
integer :: info, polarity
! An exactly or nearly singular final Schur complement is still rejected.
do polarity = 0, 1
   singular = reshape([1._rkind, 1._rkind, 1._rkind + polarity*1.0e-12_rkind, &
                       0._rkind], [2, 2])
   info = 0
   call aldlt(singular, 1.0e-10_rkind, 2, 2, 2, info)
   if (info /= -1) error stop 'Singular final Schur complement accepted'
end do
! Rescaling the last unknown must not make an independent diagonal singular.
dense = 0
dense(1, 1) = 1
dense(2, 2) = 1
dense(3, 3) = 1.0e-12_rkind
band = 0
band(1, :) = [1._rkind, 1._rkind, 1.0e-12_rkind]
rhs = matmul(dense, exact)
info = 0
call aldlt(band, 1.0e-10_rkind, 3, 2, 2, info)
if (info /= 0) error stop 'Scaled nonsingular matrix rejected'
call lyv(band, rhs, 3, 3, 2, 2)
call dwy(band, rhs, 3, 3, 2, 2)
call ltxw(band, rhs, 3, 3, 2, 2)
if (.not. all(abs(rhs - exact) <= 1.0e-12_rkind*abs(exact))) error stop 'Scaled solve'
end program
"""


class FinalPivot(NativeTest):
    def test_final_pivot_uses_its_own_scale(self):
        files = {"precision.f90": PRECISION, "DECLAR.inc": (NATIVE / "DECLAR.inc").read_text()}
        files["native.f"] = "\n".join(routine(n) for n in ("ALDLT", "LYV", "DWY", "LTXW", "SAXPY", "SDOT"))
        files["oracle.f90"] = ORACLE
        self.execute(self.compile(files))
