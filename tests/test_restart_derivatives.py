"""Restart remapping must reconstruct derivatives with smoothing disabled."""
from support import NativeTest, NATIVE, routine


class RestartDerivatives(NativeTest):
    def test_radial_remap_preserves_quadratic_flux_derivatives(self):
        files = {p.name: p.read_text() for p in NATIVE.glob('*.inc')}
        names = ('GUESS', 'BASIS1', 'SMOOTH', 'MSPLINE', 'MSPLCY',
                 'TRIDAGM', 'TRICYCM', 'RESETI', 'RESETR', 'SSCAL', 'VZERO', 'ISAMIN')
        files['native.f'] = '\n'.join(routine(n) for n in names)
        files['oracle.f'] = '''
      PROGRAM ORACLE
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMBLA.inc'
      INCLUDE 'COMBND.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMESH.inc'
      INCLUDE 'COMINT.inc'
      INCLUDE 'COMNUM.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'COMSOL.inc'
      NS=8
      NS1=NS+1
      NT=8
      NSO=NS
      NTO=NT
      NSMOOTH=0
      NOPT=1
      NSURF=1
      NSTTP=1
      R0=1.
      R0O=R0
      RZ0=0.
      RZ0O=RZ0
      CPI=ACOS(-1.)
      RC2PI=2.*CPI
      DO J=1,NT+1
         CT(J)=RC2PI*(J-1)/NT
         CTO(J)=CT(J)
         RHOS(J)=1.
      ENDDO
      DO K=1,NS1
         CSIGO(K)=REAL(K-1)/NS
         CSIG(K)=CSIGO(K)**1.5
         DO J=1,NT
            I=4*((K-1)*NT+J)-3
            CPSIO(I)=CSIGO(K)**2-1.
            CPSIO(I+1)=2.*CSIGO(K)
            CPSIO(I+2)=0.
            CPSIO(I+3)=0.
            CPSICL(I)=-999.
            CPSICL(I+1)=77.
            CPSICL(I+2)=77.
            CPSICL(I+3)=77.
            NUPDWN((K-1)*NT+J)=(K-1)*NT+J
            IF(K.EQ.1) CPSICL(I+1)=0.
            IF(K.EQ.NS1) CPSICL(I+1)=2.
         ENDDO
      ENDDO
      CALL GUESS(2)
      ERR=0.
      DO K=2,NS
         DO J=1,NT
            I=4*((K-1)*NT+J)-3
            ERR=MAX(ERR,ABS(CPSICL(I)-(CSIG(K)**2-1.)),
     &          ABS(CPSICL(I+1)-2.*CSIG(K)),
     &          ABS(CPSICL(I+2)),ABS(CPSICL(I+3)))
         ENDDO
      ENDDO
      PRINT *, 'quadratic remap error=',ERR
      IF(ERR.GT.1.E-11) STOP 1
      END
      SUBROUTINE BOUND(N,T,B)
      INCLUDE 'DECLAR.inc'
      DIMENSION T(N),B(N)
      B=1.
      END
'''
        for name in ('BNDSPL', 'CINT', 'ISOFIND', 'ISOFUN', 'PRFUNC', 'RVAR'):
            files['oracle.f'] += f'      SUBROUTINE {name}\n      STOP 2\n      END\n'
        self.execute(self.compile(files))
