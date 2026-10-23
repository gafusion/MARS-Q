         SUBROUTINE REPLAY(NS0,NT0,S,T,C,PO,PERM)
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
         INTEGER NS0,NT0,PERM((NS0+1)*NT0)
         DIMENSION S(NS0+1),T(NT0+1),C(4*(NS0+1)*NT0),
     &     PO(4*(NS0+1)*NT0)
         NS=NS0
         NT=NT0
         NS1=NS+1
         NT1=NT+1
         RC2PI=2.D0*ACOS(-1.D0)
         RC1M14=1.D-10
         DO I=1,NS1
            CSIG(I)=S(I)
         ENDDO
         DO I=1,NT1
            CT(I)=T(I)
         ENDDO
         DO I=1,NS1*NT
            NUPDWN(I)=PERM(I)
         ENDDO
         DO I=1,4*NS1*NT
            CPSICL(I)=C(I)
            CPSI(I)=C(I)
         ENDDO
         CALL NATIVE_SMOOTH
         DO I=1,4*NS1*NT
            C(I)=CPSICL(I)
            PO(I)=CPSI(I)
         ENDDO
         END
