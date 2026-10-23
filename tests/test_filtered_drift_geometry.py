"""Filtered field, inverse field and drift metric must describe one spectrum."""
import os
from pathlib import Path
from support import NativeTest

ROOT = Path(os.environ.get("MARS_ORACLE_SOURCE", Path(__file__).resolve().parents[1]))

HEADER = """      SUBROUTINE FIELD(KSMOOTH,BK,BPK,HK)
      IMPLICIT NONE
      INTEGER KSMOOTH,NR,NRP1,NCHI,KM,JS,J,K,M
      PARAMETER(NR=2,NRP1=3,NCHI=16,KM=1)
      REAL*8 BK(3,16,2),BPK(3,16,2),HK(3,16,2)
      REAL*8 G22L(3,16),G22LM(3,16),RJA(3,16),RJAM(3,16)
      REAL*8 REQ(3,16),REQM(3,16),DPSIDS(3),DPSIDSM(3)
      REAL*8 T(3),TM(3),RCHI(16),BMC(2),BMS(2),RCHIH,PI
      REAL*8 RTMP,B0K,WANT
      PI=ACOS(-1D0)
      RCHIH=2D0*PI/NCHI
      B0K=5D0
      T=1D0
      TM=1D0
      RJA=1.3D0
      RJAM=0.9D0
      G22L=0.7D0
      G22LM=1.2D0
      DPSIDS=1.1D0
      DPSIDSM=1.3D0
      BK=-999D0
      BPK=-999D0
      HK=-999D0
      DO J=1,NCHI
         RCHI(J)=(J-1)*RCHIH
         DO JS=1,NRP1
            WANT=2D0+0.2D0*COS(RCHI(J))+0.05D0*COS(3*RCHI(J))
            REQ(JS,J)=1D0/SQRT(WANT**2-
     &                   0.7D0*1.1D0**2/1.3D0**2)
            WANT=2.4D0+0.2D0*COS(RCHI(J))+0.05D0*COS(3*RCHI(J))
            REQM(JS,J)=1D0/SQRT(WANT**2-
     &                   1.2D0*1.3D0**2/0.9D0**2)
         ENDDO
      ENDDO
"""

ORACLE = """program oracle
implicit none
real*8 b(3,16,2),bp(3,16,2),h(3,16,2),x,want,ratio,scale
integer mode,j,s,k,last,first
do mode=0,1
 call field(mode,b,bp,h)
 do k=1,2
  first=2; last=3; ratio=1.3d0**2/(.7d0*1.1d0**2); scale=2d0
  if(k==2) then
   first=1; last=2; ratio=.9d0**2/(1.2d0*1.3d0**2); scale=2.4d0
  endif
  do j=1,16
   x=(j-1)*2d0*acos(-1d0)/16d0
   want=scale+.2d0*cos(x)
   if(mode==0) want=want+.05d0*cos(3*x)
   do s=first,last
    if(abs(b(s,j,k)-want)>2d-14) error stop 'Fourier field oracle'
    if(abs(h(s,j,k)*want-5d0)>2d-14) error stop 'Inverse field oracle'
    if(abs(bp(s,j,k)*h(s,j,k)-ratio*5d0)>3d-14) &
       error stop 'Drift metric uses a different field spectrum'
   enddo
  enddo
 enddo
 if(any(bp(1,:,1)/=bp(1,:,2))) error stop 'Axis metric continuation'
enddo
end program
"""


class FilteredDriftGeometry(NativeTest):
    def test_native_filtered_spectrum_and_metric(self):
        source = (ROOT / "MarsQ_2FK/kinetic.f").read_text()
        begin = source.index("C     CALCULATE FIELD STRENGTH B AND H=B0/B")
        end = source.index("C     DEFINE THE WEIGHTING FUNCTION", begin)
        files = {"native.f": HEADER + source[begin:end] + "      END\n",
                 "oracle.f90": ORACLE}
        self.execute(self.compile(files))
