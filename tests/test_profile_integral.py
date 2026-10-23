"""Native cubic-profile primitive against independent polynomial integrals."""
from pathlib import Path
from support import NativeTest, NATIVE, routine


class ProfileIntegral(NativeTest):
    def test_cubic_primitive_both_directions(self):
        source = routine("ISOFUN")
        body = source.split("         IF (NTMF0.EQ.0) THEN", 1)[1].split("C           WRITE(*,1001)", 1)[0]
        native = '         SUBROUTINE NATIVE_INTEGRATE(KN,NTMF0,NSTTP,NRFP,\n     &    AT,PSIISO,TTP,D2TTP,TMF,CIPR,CIDQ)\n         IMPLICIT NONE\n         INTEGER KN,NTMF0,NSTTP,NRFP,J7,J7P1,J8,J8P1\n         DOUBLE PRECISION AT(10),PSIISO(KN),TTP(KN),D2TTP(KN),\n     &    TMF(KN),CIPR(KN),CIDQ(KN),CPI\n         CPI=ACOS(-1.D0)\n         IF (NTMF0.EQ.0) THEN' + body + "         RETURN\n         END\n"
        native += routine("SPLINE") + "\n" + routine("TRIDAG")
        native += """
         SUBROUTINE NATIVE_SPLINE(N,X,Y,M)
         INCLUDE 'DECLAR.inc'
         DOUBLE PRECISION X(N),Y(N),M(N),A(N),B(N)
         RC1M14=1.D-14
         CALL SPLINE(X,Y,N,M,A,B)
         END
"""
        files = {"DECLAR.inc": (NATIVE / "DECLAR.inc").read_text(),
                 "CUCCCC.inc": (NATIVE / "CUCCCC.inc").read_text(),
                 "native.f": native,
                 "oracle.f90": Path(__file__).with_name("profile_integral.f90").read_text()}
        exe = self.compile(files)
        for degree in range(4):
            for direction in range(2):
                with self.subTest(degree=degree, direction=direction):
                    self.execute(exe, degree, direction)
