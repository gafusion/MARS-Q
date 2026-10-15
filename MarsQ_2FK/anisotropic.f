C KINETIC EFFECTS OF ENERGETIC PARTICLES WITH ANISOTROPIC 
C PITCH ANGLE DISTRIBUTIONS, USING NUMERICAL INTEGRATION
C ALONG PARTICLE ENERGY
C G.Z. HAO, 05-2012
C Y.Q. LIU, 06-2013

C===================================================    
C COMPUTE AND STORE NORMALISED DRIFT FREQUENCIES OF 
C PARTICLES, WHICH ARE INDEPENT OF IF0TYPE:
C   ZOMEGABT : BOUNCE  FREQUENCY OF TRAPPED PARTICLES   
C   ZOMEGABP : TRANSIT FEQUENCY OF PASSING PARTICLES
C   ZOMEGADT : PRECESSION DRIFT FREQUENCY OF TRAPPED
C
C ALSO COMPUTE \HAT\PSI0 TERMS ASSOCIATED WITH 
C FOW CORRECTION OF PASSING PARTICLES (WITH PARTICLE MASS 
C AND CHARGE DEFINED FOR THERMAL IONS, NEED A CORRECTING
C MULTIPLIER (Z_I/Z)(M/M_I) FOR HPSI0 AND ITS DERIVATIVES,
C AND (Z_I/Z)(M/M_I)^2 FOR TPSI0 AND ITS DERIVATIVE)  
C   TPSI0    : \HAT\PSI0 FOR PASSING THERMAL IONS
C   TPSI0DPSI: D TPSI0 /D PSI
C   TPSI0DLAM: D TPSI0/D LAMBDA
C   HPSI0    : \HAT\PSI0 FOR PASSING HOT IONS
C   HPSI0DPSI: D HPSI0 /D PSI
C   HPSI0DLAM: D HPSI0/D LAMBDA
C===================================================
      SUBROUTINE  KDRIFTFREQ
       
      USE KINETICM
      USE GLOBALM
      USE RCOMDM
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER JS,K,J,N,NN,KCHECK
      REAL*8  RTMPFJ,LAMH
           
      KCHECK = 0

      TPSI0     = 0.
      TPSI0DPSI = 0.
      TPSI0DLAM = 0.
      HPSI0     = 0.
      HPSI0DPSI = 0.
      HPSI0DLAM = 0.

      DO K  =1,2                                   
      DO JS =3-K,NR

      CALL KEQUIL(JS,K)     
      
C     PASSING PARTICLES
      CHIL=0.
      CHIU=2.*PI
      CALL KCHI(1)
      CALL KEQUILK(JS,K)

      IF (K.EQ.1) RTMPFJ = T(JS) *REAL(JACOBI(JS,1))/DPSIDS(JS)
      IF (K.EQ.2) RTMPFJ = TM(JS)*REAL(JACOBM(JS,1))/DPSIDSM(JS)

      DO N=1,NLAMK1(JS,K)-1
         DO J=0,1
         LAM = 0.5*((1.+WK)*LAMK1(JS,N+J,K)+(1.-WK)*LAMK1(JS,N+1-J,K))
         CALL KBTIME(JS,K,1)
         OMEGAB = 2.0*PI/RTK(NCHI2+2)
         ZOMEGABP (JS,2*N-1+J+1,K) = OMEGAB
         HPSI0    (JS,2*N-1+J,K) = OMEGAB*RTMPFJ*B0K/OMEGACI0
         TPSI0DLAM(JS,2*N-1+J,K) =-0.5*RTMPFJ*B0K/OMEGACI0/TPSI0F
         HPSI0DLAM(JS,2*N-1+J,K) =-0.25*RTMPFJ*B0K/OMEGACI0*OMEGAB**2
     &                             *HPSI0F/PI
         ENDDO
      ENDDO 
      LAM = 0.
      CALL KBTIME(JS,K,1)
      ZOMEGABP(JS,1,K) = 2.0*PI/RTK(NCHI2+2)
      ZOMEGABP(JS,2*NLAMK1(JS,K),K) = 0.

C     TRAPPED PARTICLES
      DO N=1,NLAMK0(JS,K)-1 
      DO J=0,1
         LAM = 0.5*((1.+WK)*LAMK0(JS,N+J,K)+(1.-WK)*LAMK0(JS,N+1-J,K))
         CALL KTURN(JS,K)
         CALL KCHI(0)
         CALL KEQUILK(JS,K)
         CALL KBTIME(JS,K,0)
         OMEGAB = PI/RTK(NCHI2+2)
         CALL KDRIFT(JS,K,0)
         ZOMEGABT(JS,2*N-1+J+1,K) = OMEGAB
         ZOMEGADT(JS,2*N-1+J+1,K) = DRIFT
         HATJ0   (JS,2*N-1+J,K) = HATJ0F
      ENDDO
      ENDDO      

      ZOMEGABT(JS,1,K) = 0.
      J = 1
      ZOMEGADT(JS,J,K) = ((1+WK)*ZOMEGADT(JS,J+1,K)-
     &                    (1-WK)*ZOMEGADT(JS,J+2,K))/(2.*WK)
      J = 2*NLAMK0(JS,K)
      ZOMEGABT(JS,J,K) = ((1+WK)*ZOMEGABT(JS,J-1,K)-
     &                    (1-WK)*ZOMEGABT(JS,J-2,K))/(2.*WK)
      ZOMEGADT(JS,J,K) = ((1+WK)*ZOMEGADT(JS,J-1,K)-
     &                    (1-WK)*ZOMEGADT(JS,J-2,K))/(2.*WK)

      ENDDO
      ENDDO
      
C     THE FOLLOWING IS ALTERNATIVE WAY OF COMPUTING DRIFT
C     FOR TRAPPED PARTICLES.
C     RESULT STORED IN TPSI0DPSI (TEMPORARILY), AND SHOULD
C     BE EQUAL TO ZOMEGADT
      IF (KCHECK.EQ.2) THEN
         CALL KDARRDPSI(HATJ0,TPSI0DPSI,0)

         WRITE(*,*) 'CHECK KDRIFT NEW: LAM DRIFT'
         DO N=1,NLAMK0(JS0,1)-1
         DO J=0,1
         LAM  =0.5*((1.+WK)*LAMK0(JS0,N+J,1)+(1.-WK)*LAMK0(JS0,N+1-J,1))
         LAMH =TPSI0DPSI(JS0,2*N-1+J,1)*ZOMEGABT(JS0,2*N-1+J+1,1)/PI 
         WRITE(*,110) LAM,LAMH
         ENDDO
         ENDDO
 110     FORMAT(2(E16.8,1X))
      ENDIF

C     COMPUTE REMAINING TPSI0* AND HPSI0* QUANTITIES, 
C     NEEDED ONLY FOR PASSING PARTICLES WITH FOW CORRECTION ON
      IF (IFOWP.EQ.1) THEN
         
C     COMPUTE TPSI0 BY INTEGRATING TPSI0DLAM ALONG [LAM,HMIN]
C     NOTE SPECIAL TREAMENT OF FIRST TWO GAUSSIAN POINTS
      IF (ABS(PSPECIES_FOWP(1)).GT.0.) THEN
      DO K  =1,2                                   
      DO JS =3-K,NR
      DO N  = 1,NLAMK1(JS,K)-1 
         LAMH = (LAMK1(JS,N+1,K)-LAMK1(JS,N,K))/2.
         J=0
         TPSI0(JS,2*N-1+J,K) = TPSI0(JS,2*N-1+J,K) 
     &                        +TPSI0DLAM(JS,2*N-1,K)*LAMH*WK
     &                        +TPSI0DLAM(JS,2*N,K)*LAMH
         J=1
         TPSI0(JS,2*N-1+J,K) = TPSI0(JS,2*N-1+J,K) 
     &                        +TPSI0DLAM(JS,2*N,K)*LAMH*(1.-WK)
         DO J=0,1
            DO NN=N+1,NLAMK1(JS,K)-1
               LAMH = (LAMK1(JS,NN+1,K)-LAMK1(JS,NN,K))/2.
               TPSI0(JS,2*N-1+J,K) = TPSI0(JS,2*N-1+J,K)+LAMH*(
     &                               TPSI0DLAM(JS,2*NN-1,K)+
     &                               TPSI0DLAM(JS,2*NN,K) )
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      
C     COMPUTE (D TPSI0/D PSI) AND (D HPSI0/D PSI)
      IF (ABS(PSPECIES_FOWP(1)).GT.0.) 
     &   CALL KDARRDPSI(TPSI0,TPSI0DPSI,1)

      CALL KDARRDPSI(HPSI0,HPSI0DPSI,1)
 
      ENDIF

      IF (KCHECK.EQ.1) THEN
         J = 2*NLAMK1(JS0,1)-2
         WRITE(*,120) TPSI0(JS0,J,1),TPSI0DPSI(JS0,J,1),
     &                TPSI0DLAM(JS0,J,1)
         WRITE(*,130) HPSI0(JS0,J,1),HPSI0DPSI(JS0,J,1),
     &                HPSI0DLAM(JS0,J,1) 
 120     FORMAT('TPSI0: ',3(E13.5,1X))
 130     FORMAT('HPSI0: ',3(E13.5,1X))
      ENDIF

      RETURN
      END

C===================================================    
C COMPUTE D ARR/D PSI, FOR ARR(PSI,LAMBDA,KGRID)
C AND SAVE TO DARR, WHERE LAMBDA ARE ADAPTIVE ARRAYS
C DEPENDING ON FLUX SURFACES
C PROCEDURES: 
C   FIRST CONVERT ARR TO NORMALISED LAMBDA MESH
C         ULAM=\LAMBDA/HMIN MESH USING SPLINE, 
C         ULAM IS RELATED TO UNIFORM MESH XLAM
C   SECONDLY COMPUTE PSI-DERIVATIVE IN (PSI,XLAM) SPACE, 
C   FINALLY SPLINE BACK TO (PSI,LAMBDA) SPACE 
C YQ LIU, 07-2013
C===================================================
      SUBROUTINE  KDARRDPSI(ARR,DARR,KPARTICLE)
       
      USE KINETICM
      USE GLOBALM
      USE RCOMDM 
      USE DIMENSIM
      IMPLICIT NONE
      
      INTEGER JS,K,J,N,NN,KPARTICLE,KCHECK,NLAMK2
      REAL*8  ARR(NRP1,2*NLAMK,2),DARR(NRP1,2*NLAMK,2)
      REAL*8  H1,H2,H3,F1,F2,F3,RLAM,XLAM(2*NLAMK)
      REAL*8  ULAM(2*NLAMK),UVAL(2*NLAMK),UVALK(2*NLAMK),UTMP(2*NLAMK)
      REAL*8, DIMENSION(:,:,:),ALLOCATABLE::UARR
           
      KCHECK = 0

      RLAM   = 3.5
      NLAMK2 = 2*NLAMK

C     DEFINE UNIFORM NORMALISED LAMBDA ARRAY
      XLAM(1) = 0.
      H3      = 1./DFLOAT(NLAMK2-1)
      DO N=2,NLAMK2
         XLAM(N) = XLAM(N-1)+H3
      ENDDO

      IF (KPARTICLE.EQ.1) ULAM = 1.-10.**(-RLAM*XLAM)
      IF (KPARTICLE.EQ.0) ULAM = 10.**(RLAM*(XLAM-1))

C     CONVERT ARR(PSI,LAMBDA) TO UARR(PSI,ULAM)
      ALLOCATE( UARR(NRP1,NLAMK2,2) )
 
      DO K  =1,2                                   
      DO JS =3-K,NR
         IF (KPARTICLE.EQ.0) THEN 
            NN = NLAMK0(JS,K)-1
            DO N  =1,NN 
            DO J  =0,1
               LAMM(2*N-1+J) = ((1+WK)*LAMK0(JS,N+J,K)+
     &                          (1-WK)*LAMK0(JS,N-J+1,K))/2.
            ENDDO   
            ENDDO
            LAMM  = (LAMM-HKMIN(JS,K))/(HKMAX(JS,K)-HKMIN(JS,K))
         ELSEIF (KPARTICLE.EQ.1) THEN 
            NN = NLAMK1(JS,K)-1
            DO N  =1,NN 
            DO J  =0,1
               LAMM(2*N-1+J) = ((1+WK)*LAMK1(JS,N+J,K)+
     &                          (1-WK)*LAMK1(JS,N-J+1,K))/2.
            ENDDO   
            ENDDO
            LAMM  = LAMM/HKMIN(JS,K)
         ENDIF

         UVALK(1:2*NN) = ARR(JS,1:2*NN,K)
         CALL SPLINE1D(UVAL,ULAM,NLAMK2,UVALK,LAMM,2*NN,UTMP)
         UARR(JS,1:NLAMK2,K) = UVAL
      ENDDO
      ENDDO
      UARR(1,1:NLAMK2,1)   =2.*UARR(1,1:NLAMK2,2)-UARR(2,1:NLAMK2,1)
      UARR(NRP1,1:NLAMK2,1)=2.*UARR(NR,1:NLAMK2,2)-UARR(NR,1:NLAMK2,1) 

C     COMPUTE D UARR/D PSI AT INTEGER RADIAL MESH
      DO JS =2,NR
         H1 = (CS(JS)-CS(JS-1))/2.
         H2 = (CS(JS+1)-CS(JS))/2.
         F1 = ((H1/H2*HKMIN(JS,2)-H2/H1*HKMIN(JS-1,2))/(H1+H2)-
     &        (H1-H2)*HKMIN(JS,1)/H1/H2)/DPSIDS(JS)
         IF (KPARTICLE.EQ.0)
     &   F2 = ((H1/H2*HKMAX(JS,2)-H2/H1*HKMAX(JS-1,2))/(H1+H2)-
     &        (H1-H2)*HKMAX(JS,1)/H1/H2)/DPSIDS(JS)
         DO K=2,NLAMK2-1
            IF (KPARTICLE.EQ.1) F3=-F1/HKMIN(JS,1)*ULAM(K)
     &                              *10.**(RLAM*XLAM(K))/RLAM/LOG(10.)
            IF (KPARTICLE.EQ.0) F3 =-(F1+ULAM(K)*(F2-F1))/
     &                               (HKMAX(JS,1)-HKMIN(JS,1))
     &                               /ULAM(K)/RLAM/LOG(10.)
            DARR(JS,K,1)=((H1/H2*UARR(JS,K,2)-H2/H1*UARR(JS-1,K,2))/
     &                   (H1+H2)-(H1-H2)*UARR(JS,K,1)/H1/H2)/DPSIDS(JS)
     &                   +(UARR(JS,K+1,1)-UARR(JS,K-1,1))/2./H3*F3
         ENDDO
         K=1
         IF (KPARTICLE.EQ.1) F3 =-F1/HKMIN(JS,1)*ULAM(K)
     &                            *10.**(RLAM*XLAM(K))/RLAM/LOG(10.)
         IF (KPARTICLE.EQ.0) F3 =-(F1+ULAM(K)*(F2-F1))/
     &                            (HKMAX(JS,1)-HKMIN(JS,1))
     &                            /ULAM(K)/RLAM/LOG(10.)
         DARR(JS,K,1)=((H1/H2*UARR(JS,K,2)-H2/H1*UARR(JS-1,K,2))/
     &                (H1+H2)-(H1-H2)*UARR(JS,K,1)/H1/H2)/DPSIDS(JS)
     &                +(UARR(JS,K+1,1)-UARR(JS,K,1))/H3*F3
         K=NLAMK2
         IF (KPARTICLE.EQ.1) F3 =-F1/HKMIN(JS,1)*ULAM(K)
     &                            *10.**(RLAM*XLAM(K))/RLAM/LOG(10.)
         IF (KPARTICLE.EQ.0) F3 =-(F1+ULAM(K)*(F2-F1))/
     &                            (HKMAX(JS,1)-HKMIN(JS,1))
     &                            /ULAM(K)/RLAM/LOG(10.)
         DARR(JS,K,1)=((H1/H2*UARR(JS,K,2)-H2/H1*UARR(JS-1,K,2))/
     &                (H1+H2)-(H1-H2)*UARR(JS,K,1)/H1/H2)/DPSIDS(JS)
     &                +(UARR(JS,K,1)-UARR(JS,K-1,1))/H3*F3
      ENDDO

C     COMPUTE D UARR/D PSI AT HALF-INTEGER RADIAL MESH
      DO JS=1,NR
         H1 = CS(JS+1)-CS(JS)
         F1 = (HKMIN(JS+1,1)-HKMIN(JS,1))/H1/DPSIDSM(JS)
         IF (KPARTICLE.EQ.0)
     &   F2 = (HKMAX(JS+1,1)-HKMAX(JS,1))/H1/DPSIDSM(JS)
         DO K=2,NLAMK2-1
            IF (KPARTICLE.EQ.1) F3 =-F1/HKMIN(JS,2)*ULAM(K)
     &                              *10.**(RLAM*XLAM(K))/RLAM/LOG(10.)
            IF (KPARTICLE.EQ.0) F3 =-(F1+ULAM(K)*(F2-F1))/
     &                               (HKMAX(JS,2)-HKMIN(JS,2))
     &                               /ULAM(K)/RLAM/LOG(10.)
            DARR(JS,K,2)=(UARR(JS+1,K,1)-UARR(JS,K,1))/H1/DPSIDSM(JS)
     &                   +(UARR(JS,K+1,2)-UARR(JS,K-1,2))/2./H3*F3
         ENDDO
         K=1
         IF (KPARTICLE.EQ.1) F3 =-F1/HKMIN(JS,2)*ULAM(K)
     &                            *10.**(RLAM*XLAM(K))/RLAM/LOG(10.)
         IF (KPARTICLE.EQ.0) F3 =-(F1+ULAM(K)*(F2-F1))/
     &                            (HKMAX(JS,2)-HKMIN(JS,2))
     &                            /ULAM(K)/RLAM/LOG(10.)
         DARR(JS,K,2)=(UARR(JS+1,K,1)-UARR(JS,K,1))/H1/DPSIDSM(JS)
     &                +(UARR(JS,K+1,2)-UARR(JS,K,2))/H3*F3
         K=NLAMK2
         IF (KPARTICLE.EQ.1) F3 =-F1/HKMIN(JS,2)*ULAM(K)
     &                            *10.**(RLAM*XLAM(K))/RLAM/LOG(10.)
         IF (KPARTICLE.EQ.0) F3 =-(F1+ULAM(K)*(F2-F1))/
     &                            (HKMAX(JS,2)-HKMIN(JS,2))
     &                            /ULAM(K)/RLAM/LOG(10.)
         DARR(JS,K,2)=(UARR(JS+1,K,1)-UARR(JS,K,1))/H1/DPSIDSM(JS)
     &                +(UARR(JS,K,2)-UARR(JS,K-1,2))/H3*F3
      ENDDO

C     SPLINE BACK TO (PSI,LAMBDA)-SPACE
      DO K  =1,2                                   
      DO JS =3-K,NR
         IF (KPARTICLE.EQ.0) THEN 
            NN = NLAMK0(JS,K)-1
            DO N  =1,NN 
            DO J  =0,1
               LAMM(2*N-1+J) = ((1+WK)*LAMK0(JS,N+J,K)+
     &                          (1-WK)*LAMK0(JS,N-J+1,K))/2.
            ENDDO   
            ENDDO
            LAMM  = (LAMM-HKMIN(JS,K))/(HKMAX(JS,K)-HKMIN(JS,K))
         ELSEIF (KPARTICLE.EQ.1) THEN 
            NN = NLAMK1(JS,K)-1
            DO N  =1,NN 
            DO J  =0,1
               LAMM(2*N-1+J) = ((1+WK)*LAMK1(JS,N+J,K)+
     &                          (1-WK)*LAMK1(JS,N-J+1,K))/2.
            ENDDO   
            ENDDO
            LAMM  = LAMM/HKMIN(JS,K)
         ENDIF

         UVAL  = DARR(JS,1:NLAMK2,K)
         CALL SPLINE1D(UVALK,LAMM,2*NN,UVAL,ULAM,NLAMK2,UTMP)
         DARR(JS,1:2*NN,K) = UVALK(1:2*NN)
      ENDDO
      ENDDO

      DEALLOCATE( UARR )
         
      RETURN
      END

C=========================================================
C FOR DISTRIBUTION TYPE IFOTYPE=3:    
C EVALUATE PITCH ANGLE-DEPENDENT FACTORS OF DISTRIBUTION FUNCTION F0 
C AND ITS PARTIAL DERIVATIVES FOR
C   PASSING PARTICLES(ZF1,ZF3)                       
C   TRAPPED PARTICLES(ZF2)                           
C
C   KO    = 0-4: POWER OF (ZETA-ZETA_I)
C   KD    = 0-1: ORDER OF DERIVATIVE OF (ZETA-ZETA_I) IN PSI 
C   KD    = 2  : DERIVATIVE [D^2/D PSI^2(ZETA-ZETA_I)]
C   KD    = 3  : DERIVATIVE [D^2/D PSI^2(ZETA-ZETA_I)]^2
C
C   KNBI  = 0: NORMAL NBI INJECTION MODEL: ZZETA0=0
C           2: SYMMETRIC NBI INJECTIONS (I.E. BOTH CO- AND COUNTER-INJECTIONS)
C          11: CO-TANGENTIAL INJECTION MODEL BY GORELENKOV
C          12: CO-TANGENTIAL INJECTION MODEL ~ KNBI=11, BUT ALWAYS CONTINUOUS 
C          13: CO-TANGENTIAL INJECTION MODEL, ALWAYS CONTINUOUS AND POSITIVE 
C          31: COUNTER-TANGENTIAL INJECTION MODEL BY GORELENKOV
C          32: COUNTER-TANGENTIAL INJECTION MODEL ~ KNBI=11, BUT ALWAYS CONTINUOUS 
C          33: COUNTER-TANGENTIAL INJECTION MODEL, ALWAYS CONTINUOUS AND POSITIVE 
C=========================================================
      SUBROUTINE KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZLAM,KO,KD)
       
      USE  KINETICM
      USE  GLOBALM
      USE  ANISOTROPICM
      IMPLICIT NONE
        
      INTEGER  K,L,M,J,JS,KGRID,KP,KO,KD,KCHECK,KNBJ
      REAL*8   T2,RTMP,RTMS,RTMP0,LAMB,ZLAM,XISP,ZHK,ZHMIN
      REAL*8   SMI2ME,AMA2MI,ZZALPHA1,ZZETAJ,DZETAJ
      REAL*8   ZCIA(6),ZCIB(6),ZCIC(6),ZCSA,ZCSC
      REAL*8   ZCEA(6),ZCEB(6),ZCEC(6)
      REAL*8   LAMBZH,LAMBSH
      REAL*8   RTMP1,RTMP2,RTMP3,RTMP4,RTMP5,RTMP6
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::RTMPP
       
      KCHECK = 0

      ZHK    = HKMAX(JS,KGRID)
      ZHMIN  = HKMIN(JS,KGRID)
      XISP   = SQRT(1.0-ZHMIN/ZHK)
       
      ZEPK   = ZEPKO

      KNBJ   = KNBI(KP)
      ZZETAJ = ZZETA0(KP)
      DZETAJ = DZETA0(KP)

      ALLOCATE( RTMPP(NEPK2,6) )
       
      ZCIA = 0.
      ZCIB = 0.
      ZCIC = 0.
      ZCEA = 0.
      ZCEB = 0.
      ZCEC = 0.
      ZCSA = 0.
      ZCSC = 0.

      IF (KNBJ.EQ.0) ZZETAJ = 0.
      IF (KNBJ.EQ.11.OR.KNBJ.EQ.12.OR.KNBJ.EQ.13) ZZETAJ= ABS(ZZETAJ)
      IF (KNBJ.EQ.31.OR.KNBJ.EQ.32.OR.KNBJ.EQ.33) ZZETAJ=-ABS(ZZETAJ)

      ZCIB(1:4) = 1.
      ZCEB(1)   = ZZETAJ
      ZCEB(2)   = 2.-ZZETAJ
      ZCEB(3)   =-ZZETAJ
      ZCEB(4)   =-2.+ZZETAJ
      ZCEA(1:4) = ZCEB(1:4)
      ZCEC(1:4) = ZCEB(1:4)
      
      IF (KNBJ.EQ.0.OR.KNBJ.EQ.2) THEN
         ZCIA(1:4) = 1. 
         ZCIC(1:4) = 1.
      ENDIF 
       
      IF (KNBJ.EQ.11) THEN
         ZCIA(1:2) = 2. 
         ZCIA(5:6) =-1.
         ZCSA      = 2.
         ZCEA(5)   = ZCSA*XISP-ZZETAJ
         ZCEA(6)   = ZCSA*XISP-2.+ZZETAJ
       
         ZCIC(5:6) = 1.
         ZCSC      =-2. 
         ZCEC(5)   = ZCSC*XISP+ZZETAJ
         ZCEC(6)   = ZCSC*XISP+2.-ZZETAJ
      ENDIF
       
      IF (KNBJ.EQ.12) THEN
         ZCIA(1:2) = 2. 
         ZCIA(3:4) = 1. 
         ZCIA(5:6) =-1. 
         ZCSA      = 2.
         ZCEA(5)   = ZCSA*XISP-ZZETAJ
         ZCEA(6)   = ZCSA*XISP-2.+ZZETAJ
       
         ZCIC(1:2) = 1.
         ZCIC(5:6) = 1.
         ZCSC      =-2.
         ZCEC(5)   = ZCSC*XISP+ZZETAJ
         ZCEC(6)   = ZCSC*XISP+2.-ZZETAJ
      ENDIF
       
      IF (KNBJ.EQ.13) THEN       
         ZCIA(1:4) = 1. 

         ZCIC(1:2) = 1.
         ZCIC(5:6) = 1.
         ZCSC      =-2.
         ZCEC(5)   = ZCSC*XISP+ZZETAJ
         ZCEC(6)   = ZCSC*XISP+2.-ZZETAJ
      ENDIF
       
      IF (KNBJ.EQ.31) THEN
         ZCIA(5:6) = 1.
         ZCSA      = 2.
         ZCEA(5)   = ZCSA*XISP+ZZETAJ
         ZCEA(6)   = ZCSA*XISP+2.-ZZETAJ

         ZCIC(1:2) = 2. 
         ZCIC(5:6) =-1. 
         ZCSC      =-2.
         ZCEC(5)   = ZCSC*XISP-ZZETAJ
         ZCEC(6)   = ZCSC*XISP-2.+ZZETAJ
      ENDIF
       
      IF (KNBJ.EQ.32) THEN
         ZCIA(1:2) = 1.
         ZCIA(5:6) = 1.
         ZCSA      = 2.
         ZCEA(5)   = ZCSA*XISP+ZZETAJ
         ZCEA(6)   = ZCSA*XISP+2.-ZZETAJ

         ZCIC(1:2) = 2. 
         ZCIC(3:4) = 1. 
         ZCIC(5:6) =-1. 
         ZCSC      =-2.
         ZCEC(5)   = ZCSC*XISP-ZZETAJ
         ZCEC(6)   = ZCSC*XISP-2.+ZZETAJ
      ENDIF
       
      IF (KNBJ.EQ.33) THEN
         ZCIA(1:2) = 1.
         ZCIA(5:6) = 1.
         ZCSA      = 2.
         ZCEA(5)   = ZCSA*XISP+ZZETAJ
         ZCEA(6)   = ZCSA*XISP+2.-ZZETAJ

         ZCIC(1:4) = 1. 
      ENDIF

C     T2 IS EPSLONC (NORMALIZED)              

      SMI2ME = SQRT(ESPECIES_M(1)/ESPECIES_M(2))
      AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
      RTMP   = (3.*SQRT(PI)/4.*SMI2ME)**(2./3.)*AMA2MI
      T2     = RTMP*ZZALPHA1
       
      LAMB   = ZLAM
      LAMBZH = ZLAM/ZHK
      LAMBSH = SQRT(1.-LAMBZH)

      ZDTK2  = DZETAJ**2-LOG(ZEPK**1.5*(1.0+T2**1.5)
     &         /(ZEPK**1.5+T2**1.5))/3.0
       
C     FOR PASSING PARTICLES    
      IF(LAMB.GE.0.0.AND.LAMB.LE.ZHMIN)THEN

C     PASSING PARTICLE WITH SIGMA=+1
      RTMPP(:,1) = -(LAMBSH-ZCEA(1))**2/ZDTK2
      RTMPP(:,2) = -(LAMBSH-ZCEA(2))**2/ZDTK2
      RTMPP(:,3) = -(LAMBSH-ZCEA(3))**2/ZDTK2
      RTMPP(:,4) = -(LAMBSH-ZCEA(4))**2/ZDTK2
      RTMPP(:,5) = -(LAMBSH-ZCEA(5))**2/ZDTK2
      RTMPP(:,6) = -(LAMBSH-ZCEA(6))**2/ZDTK2

      IF (KD.EQ.0) THEN
         FESS1(:,KO,KD) = 
     &        (LAMBSH-ZCEA(1))**KO*ZCIA(1)*EXP(RTMPP(:,1))
     &       +(LAMBSH-ZCEA(2))**KO*ZCIA(2)*EXP(RTMPP(:,2))
     &       +(LAMBSH-ZCEA(3))**KO*ZCIA(3)*EXP(RTMPP(:,3))
     &       +(LAMBSH-ZCEA(4))**KO*ZCIA(4)*EXP(RTMPP(:,4))
     &       +(LAMBSH-ZCEA(5))**KO*ZCIA(5)*EXP(RTMPP(:,5))
     &       +(LAMBSH-ZCEA(6))**KO*ZCIA(6)*EXP(RTMPP(:,6))
      ELSEIF (KD.EQ.1) THEN
C        RTMP=D ZETA/D PSI; RTMP0=D ZETA_S+ / D_PSI
         RTMP  = ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP0 =-0.5/XISP*ZDHMINMAXDPSI(JS,KGRID)
         RTMS  = RTMP-ZCSA*RTMP0

         FESS1(:,KO,KD) = 
     &        RTMP*(LAMBSH-ZCEA(1))**KO*ZCIA(1)*EXP(RTMPP(:,1)) 
     &       +RTMP*(LAMBSH-ZCEA(2))**KO*ZCIA(2)*EXP(RTMPP(:,2))
     &       +RTMP*(LAMBSH-ZCEA(3))**KO*ZCIA(3)*EXP(RTMPP(:,3))
     &       +RTMP*(LAMBSH-ZCEA(4))**KO*ZCIA(4)*EXP(RTMPP(:,4))     
     &       +RTMS*(LAMBSH-ZCEA(5))**KO*ZCIA(5)*EXP(RTMPP(:,5))     
     &       +RTMS*(LAMBSH-ZCEA(6))**KO*ZCIA(6)*EXP(RTMPP(:,6))     
      ELSEIF (KD.EQ.2) THEN
         RTMP  = ZDHMAXDPSI2(JS,KGRID)*LAMB/LAMBSH/2. - 
     &           (ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/2.)**2/LAMBSH**3
         RTMP0 =-0.5/XISP*ZDHMINMAXDPSI2(JS,KGRID)
     &           -0.25/XISP**3*ZDHMINMAXDPSI(JS,KGRID)**2
         RTMS  = RTMP-ZCSA*RTMP0
         RTMP1 = (LAMBSH-ZCEA(1))*RTMP
         RTMP2 = (LAMBSH-ZCEA(2))*RTMP
         RTMP3 = (LAMBSH-ZCEA(3))*RTMP
         RTMP4 = (LAMBSH-ZCEA(4))*RTMP
         RTMP5 = (LAMBSH-ZCEA(5))*RTMS
         RTMP6 = (LAMBSH-ZCEA(6))*RTMS
         
         RTMP  = ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP0 =-0.5/XISP*ZDHMINMAXDPSI(JS,KGRID)
         RTMS  = RTMP-ZCSA*RTMP0
         RTMP1 = RTMP1 + RTMP**2
         RTMP2 = RTMP2 + RTMP**2
         RTMP3 = RTMP3 + RTMP**2
         RTMP4 = RTMP4 + RTMP**2
         RTMP5 = RTMP5 + RTMS**2
         RTMP6 = RTMP6 + RTMS**2
         
         FESS1(:,0,KD) = 
     &        RTMP1*ZCIA(1)*EXP(RTMPP(:,1)) 
     &       +RTMP2*ZCIA(2)*EXP(RTMPP(:,2))
     &       +RTMP3*ZCIA(3)*EXP(RTMPP(:,3))
     &       +RTMP4*ZCIA(4)*EXP(RTMPP(:,4))     
     &       +RTMP5*ZCIA(5)*EXP(RTMPP(:,5))     
     &       +RTMP6*ZCIA(6)*EXP(RTMPP(:,6))     
      ELSEIF (KD.EQ.3) THEN
         RTMP  = ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP0 =-0.5/XISP*ZDHMINMAXDPSI(JS,KGRID)
         RTMS  = RTMP-ZCSA*RTMP0
         RTMP1 = ((LAMBSH-ZCEA(1))*RTMP)**2
         RTMP2 = ((LAMBSH-ZCEA(2))*RTMP)**2
         RTMP3 = ((LAMBSH-ZCEA(3))*RTMP)**2
         RTMP4 = ((LAMBSH-ZCEA(4))*RTMP)**2
         RTMP5 = ((LAMBSH-ZCEA(5))*RTMS)**2
         RTMP6 = ((LAMBSH-ZCEA(6))*RTMS)**2

         FESS1(:,0,KD) = 
     &        RTMP1*ZCIA(1)*EXP(RTMPP(:,1)) 
     &       +RTMP2*ZCIA(2)*EXP(RTMPP(:,2))
     &       +RTMP3*ZCIA(3)*EXP(RTMPP(:,3))
     &       +RTMP4*ZCIA(4)*EXP(RTMPP(:,4))     
     &       +RTMP5*ZCIA(5)*EXP(RTMPP(:,5))     
     &       +RTMP6*ZCIA(6)*EXP(RTMPP(:,6))     
      ENDIF
              
C     PASSING PARTICLE WITH SIGMA=-1
      RTMPP(:,1) = -(-LAMBSH-ZCEC(1))**2/ZDTK2
      RTMPP(:,2) = -(-LAMBSH-ZCEC(2))**2/ZDTK2
      RTMPP(:,3) = -(-LAMBSH-ZCEC(3))**2/ZDTK2
      RTMPP(:,4) = -(-LAMBSH-ZCEC(4))**2/ZDTK2
      RTMPP(:,5) = -(-LAMBSH-ZCEC(5))**2/ZDTK2
      RTMPP(:,6) = -(-LAMBSH-ZCEC(6))**2/ZDTK2
                
      IF (KD.EQ.0) THEN
         FESS3(:,KO,KD) = 
     &        (-LAMBSH-ZCEC(1))**KO*ZCIC(1)*EXP(RTMPP(:,1))
     &       +(-LAMBSH-ZCEC(2))**KO*ZCIC(2)*EXP(RTMPP(:,2))
     &       +(-LAMBSH-ZCEC(3))**KO*ZCIC(3)*EXP(RTMPP(:,3))
     &       +(-LAMBSH-ZCEC(4))**KO*ZCIC(4)*EXP(RTMPP(:,4))
     &       +(-LAMBSH-ZCEC(5))**KO*ZCIC(5)*EXP(RTMPP(:,5))
     &       +(-LAMBSH-ZCEC(6))**KO*ZCIC(6)*EXP(RTMPP(:,6))
      ELSEIF (KD.EQ.1) THEN
         RTMP  =-ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP0 =-0.5/XISP*ZDHMINMAXDPSI(JS,KGRID)
         RTMS  = RTMP - ZCSC*RTMP0
      
         FESS3(:,KO,KD) = 
     &        RTMP*(-LAMBSH-ZCEC(1))**KO*ZCIC(1)*EXP(RTMPP(:,1)) 
     &       +RTMP*(-LAMBSH-ZCEC(2))**KO*ZCIC(2)*EXP(RTMPP(:,2))
     &       +RTMP*(-LAMBSH-ZCEC(3))**KO*ZCIC(3)*EXP(RTMPP(:,3))
     &       +RTMP*(-LAMBSH-ZCEC(4))**KO*ZCIC(4)*EXP(RTMPP(:,4))     
     &       +RTMS*(-LAMBSH-ZCEC(5))**KO*ZCIC(5)*EXP(RTMPP(:,5))     
     &       +RTMS*(-LAMBSH-ZCEC(6))**KO*ZCIC(6)*EXP(RTMPP(:,6))     
      ELSEIF (KD.EQ.2) THEN
         RTMP  =-ZDHMAXDPSI2(JS,KGRID)*LAMB/LAMBSH/2. + 
     &          (ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/2.)**2/LAMBSH**3
         RTMP0 =-0.5/XISP*ZDHMINMAXDPSI2(JS,KGRID)
     &           -0.25/XISP**3*ZDHMINMAXDPSI(JS,KGRID)**2
         RTMS  = RTMP - ZCSC*RTMP0
         RTMP1 = (-LAMBSH-ZCEC(1))*RTMP
         RTMP2 = (-LAMBSH-ZCEC(2))*RTMP
         RTMP3 = (-LAMBSH-ZCEC(3))*RTMP
         RTMP4 = (-LAMBSH-ZCEC(4))*RTMP  
         RTMP5 = (-LAMBSH-ZCEC(5))*RTMS  
         RTMP6 = (-LAMBSH-ZCEC(6))*RTMS  

         RTMP   =-ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP0  =-0.5/XISP*ZDHMINMAXDPSI(JS,KGRID)
         RTMS   = RTMP - ZCSC*RTMP0
         RTMP1 = RTMP1 + RTMP**2
         RTMP2 = RTMP2 + RTMP**2
         RTMP3 = RTMP3 + RTMP**2
         RTMP4 = RTMP4 + RTMP**2
         RTMP5 = RTMP5 + RTMS**2
         RTMP6 = RTMP6 + RTMS**2
      
         FESS3(:,0,KD) = 
     &        RTMP1*ZCIC(1)*EXP(RTMPP(:,1)) 
     &       +RTMP2*ZCIC(2)*EXP(RTMPP(:,2))
     &       +RTMP3*ZCIC(3)*EXP(RTMPP(:,3))
     &       +RTMP4*ZCIC(4)*EXP(RTMPP(:,4))     
     &       +RTMP5*ZCIC(5)*EXP(RTMPP(:,5))     
     &       +RTMP6*ZCIC(6)*EXP(RTMPP(:,6))     
      ELSEIF (KD.EQ.3) THEN
         RTMP  =-ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP0 =-0.5/XISP*ZDHMINMAXDPSI(JS,KGRID)
         RTMS  = RTMP - ZCSC*RTMP0
         RTMP1 = ((-LAMBSH-ZCEC(1))*RTMP)**2
         RTMP2 = ((-LAMBSH-ZCEC(2))*RTMP)**2
         RTMP3 = ((-LAMBSH-ZCEC(3))*RTMP)**2
         RTMP4 = ((-LAMBSH-ZCEC(4))*RTMP)**2  
         RTMP5 = ((-LAMBSH-ZCEC(5))*RTMS)**2  
         RTMP6 = ((-LAMBSH-ZCEC(6))*RTMS)**2  
      
         FESS3(:,0,KD) = 
     &        RTMP1*ZCIC(1)*EXP(RTMPP(:,1)) 
     &       +RTMP2*ZCIC(2)*EXP(RTMPP(:,2))
     &       +RTMP3*ZCIC(3)*EXP(RTMPP(:,3))
     &       +RTMP4*ZCIC(4)*EXP(RTMPP(:,4))     
     &       +RTMP5*ZCIC(5)*EXP(RTMPP(:,5))     
     &       +RTMP6*ZCIC(6)*EXP(RTMPP(:,6))     
      ENDIF

C     FOR TRAPPED PARTICLES      
C     DUE TO SYMMETRY, COMPUTE ONLY HALF ORBIT WITH SIGMA=+1
      ELSEIF (LAMB.GT.ZHMIN.AND.LAMB.LE.ZHK) THEN

      RTMPP(:,1)  = -(LAMBSH-ZCEB(1))**2/ZDTK2
      RTMPP(:,2)  = -(LAMBSH-ZCEB(2))**2/ZDTK2
      RTMPP(:,3)  = -(LAMBSH-ZCEB(3))**2/ZDTK2
      RTMPP(:,4)  = -(LAMBSH-ZCEB(4))**2/ZDTK2
             
      IF (KD.EQ.0) THEN
         FESS2(:,KO,KD) = 
     &        (LAMBSH-ZCEB(1))**KO*ZCIB(1)*EXP(RTMPP(:,1))
     &       +(LAMBSH-ZCEB(2))**KO*ZCIB(2)*EXP(RTMPP(:,2))
     &       +(LAMBSH-ZCEB(3))**KO*ZCIB(3)*EXP(RTMPP(:,3))
     &       +(LAMBSH-ZCEB(4))**KO*ZCIB(4)*EXP(RTMPP(:,4))
      ELSEIF (KD.EQ.1) THEN
         RTMP = ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
      
         FESS2(:,KO,KD) = 
     &        RTMP*(LAMBSH-ZCEB(1))**KO*ZCIB(1)*EXP(RTMPP(:,1)) 
     &       +RTMP*(LAMBSH-ZCEB(2))**KO*ZCIB(2)*EXP(RTMPP(:,2))
     &       +RTMP*(LAMBSH-ZCEB(3))**KO*ZCIB(3)*EXP(RTMPP(:,3))
     &       +RTMP*(LAMBSH-ZCEB(4))**KO*ZCIB(4)*EXP(RTMPP(:,4))
      ELSEIF (KD.EQ.2) THEN
         RTMP  = ZDHMAXDPSI2(JS,KGRID)*LAMB/LAMBSH/2. - 
     &           (ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/2.)**2/LAMBSH**3
         RTMP1 = (LAMBSH-ZCEB(1))*RTMP
         RTMP2 = (LAMBSH-ZCEB(2))*RTMP
         RTMP3 = (LAMBSH-ZCEB(3))*RTMP
         RTMP4 = (LAMBSH-ZCEB(4))*RTMP  
      
         RTMP  = ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP1 = RTMP1 + RTMP**2
         RTMP2 = RTMP2 + RTMP**2
         RTMP3 = RTMP3 + RTMP**2
         RTMP4 = RTMP4 + RTMP**2

         FESS2(:,0,KD) = 
     &        RTMP1*ZCIB(1)*EXP(RTMPP(:,1)) 
     &       +RTMP2*ZCIB(2)*EXP(RTMPP(:,2))
     &       +RTMP3*ZCIB(3)*EXP(RTMPP(:,3))
     &       +RTMP4*ZCIB(4)*EXP(RTMPP(:,4))
      ELSEIF (KD.EQ.3) THEN
         RTMP  = ZDHMAXDPSI(JS,KGRID)*LAMB/ZHK**2/LAMBSH/2.
         RTMP1 = ((LAMBSH-ZCEB(1))*RTMP)**2
         RTMP2 = ((LAMBSH-ZCEB(2))*RTMP)**2
         RTMP3 = ((LAMBSH-ZCEB(3))*RTMP)**2
         RTMP4 = ((LAMBSH-ZCEB(4))*RTMP)**2  
      
         FESS2(:,0,KD) = 
     &        RTMP1*ZCIB(1)*EXP(RTMPP(:,1)) 
     &       +RTMP2*ZCIB(2)*EXP(RTMPP(:,2))
     &       +RTMP3*ZCIB(3)*EXP(RTMPP(:,3))
     &       +RTMP4*ZCIB(4)*EXP(RTMPP(:,4))
      ENDIF
      
      ENDIF
     
      DEALLOCATE(RTMPP)

      RETURN
      END SUBROUTINE KF0_TYPE3

C===================================================    
C     A LOOP FOR FINDING THE VALUE OF ALPHAA1 FOR 
C     A GIVEN F0 OF FAST IONS                     
C     USING BI-SECTION METHOD
C     APPLIES FOR IF0TYPE=3,5,7
C===================================================
      SUBROUTINE  KALPHAA_TYPE3(KP)
       
      USE KINETICM
      USE GLOBALM
      USE RCOMDM
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER KP,JS,K,J,N,KCHECK
      REAL*8  RTMP,RTMP1,RTMP2,H2,ZZPAK
      REAL*8  SMI2ME,AMA2MI,DRHO,H1,AX1,AX2,AX3,AF1,AF2,AF3
      REAL*8  ZALPHA3
           
      KCHECK =0

      DO K  =1,2                                   
      DO JS =3-K,NR

      H1   = ESPECIES_PREF(JS,K,KP)/ESPECIES_DENF(JS,K,KP)

      AX1       = 1.E-8
      AF1       = (1.0-ALPHAP)*ZALPHA3(JS,K,KP,AX1)-AX1*H1
      AX2       = 1.E+0
      AF2       = (1.0-ALPHAP)*ZALPHA3(JS,K,KP,AX2)-AX2*H1
      AX3       = 0.5*(AX1+AX2)
      AF3       = (1.0-ALPHAP)*ZALPHA3(JS,K,KP,AX3)-AX3*H1
      
      DO J=1,1000
         IF (AF1*AF3.LE.0.) THEN
            AX2 = AX3
            AF2 = AF3
         ELSEIF (AF2*AF3.LE.0.) THEN
            AX1 = AX3
            AF1 = AF3
         ELSE
            WRITE(*,*) 'KDIAMAG: K=',K,', JS=',JS,', H1=',H1
            WRITE(*,*) 'AX1=',AX1,'  AX2=',AX2,'  AX3=',AX3
            WRITE(*,*) 'AF1=',AF1,'  AF2=',AF2,'  AF3=',AF3
            STOP 'KALPHAA_TYPE3:ALPHAA1'
         ENDIF
         AX3 = (AX1+AX2)*0.5
         
         AF3 = (1.0-ALPHAP)*ZALPHA3(JS,K,KP,AX3)-AX3*H1
         IF (ABS((AX1-AX2)/AX3).LE.1.E-5) EXIT 
      ENDDO 

      ALPHAA1(JS,K,KP) = AX3

      CALL   ZI1I3(JS,K,KP,AX3)

      ALPHAA2(JS,K,KP) = ALPHAA1(JS,K,KP)*ALPHAP/(1.-ALPHAP)       
      ENDDO
      ENDDO
      
      J = NRP1
      ALPHAA1(J,1,KP)  = 2.*ALPHAA1(J-1,2,KP)-ALPHAA1(J-1,1,KP)
      ALPHAA2(J,1,KP)  = 2.*ALPHAA2(J-1,2,KP)-ALPHAA2(J-1,1,KP)
      ALPHAA3(J,1,KP)  = 2.*ALPHAA3(J-1,2,KP)-ALPHAA3(J-1,1,KP)

      J = 1
      ALPHAA1(J,1,KP)  = 2.*ALPHAA1(J,2,KP)-ALPHAA1(J+1,1,KP)
      ALPHAA2(J,1,KP)  = 2.*ALPHAA2(J,2,KP)-ALPHAA2(J+1,1,KP)
      ALPHAA3(J,1,KP)  = 2.*ALPHAA3(J,2,KP)-ALPHAA3(J+1,1,KP)

      SMI2ME = SQRT(ESPECIES_M(1)/ESPECIES_M(2))
      AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
      DRHO   = (3.*SQRT(PI)/4.*SMI2ME)**(2./3.)*AMA2MI
      EPSLONCA(:,:,KP) = DRHO*ALPHAA1(:,:,KP)

      RETURN
      END SUBROUTINE KALPHAA_TYPE3
       
C===================================================    
C CORE PROCEDURE FOR 
C COMPUTING PARTCILE PHASE SPACE INTERGRALS AT 
C A GIVEN VALUE OF ZZALPHA1, FOR
C I1 = SURFACE AVERAGED DENSITY 
C I3 = SURFACE AVERAGED PRESSURE
C===================================================
      SUBROUTINE ZZZI1I3(JS,KGRID,KP,ZZALPHA1)
                    
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER   JS,KGRID,KP,KCHECK
      INTEGER   N,J,K
      REAL*8    HEPK,ZHLAM,H3,
     &          ZF1L,ZF3L,ZF1LP,ZF3LP,
     &          ZF1J,ZF3J,ZF1JP,ZF3JP,
     &          RTMP,RTMP1,RTMP2,
     &          ZZALPHA1,DZETAJ,
     &          ZF2L,ZF2LP,ZF2J,ZF2JP,
     &          ZLAMB,SMI2ME,AMA2MI,T2 
      REAL*8,  DIMENSION(:),ALLOCATABLE::FUNN,FUNP

      KCHECK  = 0

      ZEPK    = ZEPKO

      ZF1L    = 0.
      ZF3L    = 0.
      ZF1LP   = 0.
      ZF3LP   = 0.

      DZETAJ  = DZETA0(KP)

      SMI2ME = SQRT(ESPECIES_M(1)/ESPECIES_M(2))
      AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
      RTMP   = (3.*SQRT(PI)/4.*SMI2ME)**(2./3.)*AMA2MI
      T2     = RTMP*ZZALPHA1

      FEPK  = 1./(ZEPK**1.5+T2**1.5)
      ZDTK2 = DZETAJ**2-LOG(ZEPK**1.5*(1.+T2**1.5)*FEPK)/3.
      ZDTK  = SQRT(ZDTK2)

      ALLOCATE(FUNN(NEPK2), FUNP(NEPK2))
      FUNN = FEPK/ZDTK*ZEPK**0.5
      FUNP = FUNN*ZEPK

C     FIRST GO THROUGH PASSING PARTICLES
      DO N   = 1,NLAMK1(JS,KGRID)-1
      ZHLAM  = LAMK1(JS,N+1,KGRID)-LAMK1(JS,N,KGRID)
      DO  J = 0,1
         ZLAMB  = 0.5*((1.+ZW)*LAMK1(JS,N+J,KGRID)
     &            +(1.-ZW)*LAMK1(JS,N+1-J,KGRID))
         ZF1J   = 0.0
         ZF3J   = 0.0
         ZF1JP  = 0.0
         ZF3JP  = 0.0

         CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZLAMB,0,0)

         DO K=1,NEPK-1
         HEPK  = EPK(K+1)-EPK(K)
         RTMP1 = 0.5*HEPK*( FESS1(2*K-1,0,0)*FUNN(2*K-1) + 
     &                      FESS1(2*K,0,0)  *FUNN(2*K) )
         RTMP2 = 0.5*HEPK*( FESS3(2*K-1,0,0)*FUNN(2*K-1) + 
     &                      FESS3(2*K,0,0)  *FUNN(2*K) )
         ZF1J  =  ZF1J+RTMP1
         ZF3J  =  ZF3J+RTMP2
          
         RTMP1 = 0.5*HEPK*( FESS1(2*K-1,0,0)*FUNP(2*K-1) + 
     &                      FESS1(2*K,0,0)  *FUNP(2*K) )
         RTMP2 = 0.5*HEPK*( FESS3(2*K-1,0,0)*FUNP(2*K-1) + 
     &                      FESS3(2*K,0,0)  *FUNP(2*K) )
         ZF1JP  = ZF1JP+RTMP1
         ZF3JP  = ZF3JP+RTMP2
         ENDDO

         H3     = ZOMEGABP(JS,2*N-1+J+1,KGRID)
          
         RTMP   = 0.5*ZHLAM/H3
         ZF1L   = ZF1L  + ZF1J* RTMP
         ZF3L   = ZF3L  + ZF3J* RTMP
         ZF1LP  = ZF1LP + ZF1JP*RTMP
         ZF3LP  = ZF3LP + ZF3JP*RTMP
      ENDDO
      ENDDO
       
C     NEXT GO THROUGH TRAPPED PARTICLES       
      ZF2L   = 0.
      ZF2LP  = 0.
      DO  N  = 1,NLAMK0(JS,KGRID)-1 
      ZHLAM  = LAMK0(JS,N+1,KGRID)-LAMK0(JS,N,KGRID)
      DO  J  = 0,1
         ZLAMB = 0.5*((1.+ZW)*LAMK0(JS,N+J,KGRID)
     &           +(1.-ZW)*LAMK0(JS,N+1-J,KGRID))
         ZF2J  = 0.0
         ZF2JP = 0.0  

         CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZLAMB,0,0)

         DO K=1,NEPK-1
         HEPK   = EPK(K+1)-EPK(K)
         RTMP1 = 0.5*HEPK*( FESS2(2*K-1,0,0)*FUNN(2*K-1) + 
     &                      FESS2(2*K,0,0)  *FUNN(2*K) )
         ZF2J   = ZF2J+RTMP1
          
         RTMP2 = 0.5*HEPK*( FESS2(2*K-1,0,0)*FUNP(2*K-1) + 
     &                      FESS2(2*K,0,0)  *FUNP(2*K) )
         ZF2JP  = ZF2JP+RTMP2
         ENDDO

         H3     = ZOMEGABT(JS,2*N-1+J+1,KGRID)

         RTMP   = 0.5*ZHLAM/H3
         ZF2L   = ZF2L +(ZF2J)*RTMP
         ZF2LP  = ZF2LP+(ZF2JP)*RTMP
      ENDDO
      ENDDO
      
C     0.5 COMES FROM \hat_g FACTOR FOR TRAPPED PARTICLES
      RTMP2     = 0.5 
      ZZZI1     = ZF1L +2.0*RTMP2*ZF2L +ZF3L
      ZZZI3     = ZF1LP+2.0*RTMP2*ZF2LP+ZF3LP
      ZZZALPHA3 = ZZZI3/ZZZI1*2./3.

      DEALLOCATE(FUNN,FUNP)

      RETURN
      END SUBROUTINE ZZZI1I3

C===================================================    
C CORE PROCEDURE FOR 
C COMPUTING PARTCILE PHASE SPACE INTERGRALS AT 
C A GIVEN VALUE OF ZZALPHA1, FOR
C I1 = SURFACE AVERAGED DENSITY
C I3 = SURFACE AVERAGED PRESSURE
C FOR TYPE5
C===================================================
      SUBROUTINE ZITYPE5(JS,KGRID,KP,ZZALPHA1)
                    
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER   JS,KGRID,KP,KCHECK,K
      REAL*8    HEPK,H1,H2,H3,H4,H5,H6,H7,
     &          ZF1L,ZF3L,RTMP1,RTMP3,ZZALPHA1
      REAL*8,  DIMENSION(:),ALLOCATABLE::FUNN,FUNP,RPR2

      KCHECK  = 0

      ZEPK    = ZEPKO

      ZF1L    = 0.
      ZF3L    = 0.

      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)                           !=P0
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI !=N0
      H6 = 8.1872E-14                                     !=me*c^2
      H5 = H1/H2/H6                                       !=T0/(me*c^2)
      H7 = ESPECIES_TEM(JS,KGRID,2)/ZZALPHA1*H5           !=EPS_A/(me*c^2)

      ALLOCATE(FUNN(NEPK2),FUNP(NEPK2),RPR2(NEPK))
      CALL SPLINE1D(FEPK,ZEPK,NEPK2,ESPECIES_REE(:,2),ESPECIES_REE(:,1),
     &              NEPK,RPR2)

C     FEPK = FEPK*(ZEPK*H7+1.0)**(-1.5)
      FUNN = FEPK*ZEPK**0.5
      FUNP = FUNN*ZEPK

C     GO THROUGH BOTH PASSING AND TRAPPED PARTICLES
      DO K=1,NEPK-1
         HEPK  = EPK(K+1)-EPK(K)
         RTMP1 = 0.5*HEPK*( FUNN(2*K-1) + FUNN(2*K) )
         RTMP3 = 0.5*HEPK*( FUNP(2*K-1) + FUNP(2*K) )
         
         ZF1L  = ZF1L + RTMP1
         ZF3L  = ZF3L + RTMP3
      ENDDO
          
      ZZZI1     = ZF1L 
      ZZZI3     = ZF3L
      ZZZALPHA3 = ZZZI3/ZZZI1*2./3.

      DEALLOCATE(FUNN,FUNP,RPR2)

      RETURN
      END SUBROUTINE ZITYPE5

C===================================================    
C CORE PROCEDURE FOR 
C COMPUTING PARTCILE PHASE SPACE INTERGRALS AT 
C A GIVEN VALUE OF ZZALPHA1, FOR
C I1 = SURFACE AVERAGED PARALLEL CURRENT DENSITY 
C I3 = SURFACE AVERAGED PRESSURE
C FOR TYPE6
C===================================================
      SUBROUTINE ZITYPE6(JS,KGRID,KP,ZZALPHA1)
                    
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      USE REORBITM
      IMPLICIT NONE
      
      INTEGER   JS,KGRID,KP,KCHECK
      INTEGER   N,J,K
      REAL*8    HEPK,ZHLAM,ZF1L,ZF2L,ZF3L,ZF1J,ZF2J,ZF3J,
     &          ZZALPHA1,ZLAMB,AMA2MI,RTMP,RTMP1,RTMP2
      REAL*8,  DIMENSION(:),ALLOCATABLE::FUNN,FUNP,ZREP,ZREG

      KCHECK  = 0

      ALLOCATE(FUNN(NEPK2), FUNP(NEPK2), ZREP(NEPK2), ZREG(NEPK2))

      ZEPK    = ZEPKO
      ZREP    = 1.0/SQRT(1.0/ZEPK-1.0)
      ZREG    = 1.0/SQRT(1.0-ZEPK)
      ZREP(NEPK2) = ZREP(NEPK2-1)
      ZREG(NEPK2) = ZREG(NEPK2-1)

      ZF1L    = 0.
      ZF2L    = 0.
      ZF3L    = 0.

      AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
      RTMP   = AMA2MI*(C_VA)**2/2.0

      FEPK  = EXP(-(ZREP-RE_PMAX)**2/RE_CONST(7))*ZEPK

C     FIRST GO THROUGH PASSING PARTICLES
      DO N   = 1,NLAMK1(JS,KGRID)-1
      ZHLAM  = LAMK1(JS,N+1,KGRID)-LAMK1(JS,N,KGRID)
      DO  J = 0,1
         ZLAMB  = 0.5*((1.+ZW)*LAMK1(JS,N+J,KGRID)
     &            +(1.-ZW)*LAMK1(JS,N+1-J,KGRID))

         ZDTK = ZREP**2/ZREG/RE_CONST(6)
         FUNP = EXP(-ZLAMB*EQ1_H(JS,KGRID)/2.*ZDTK)*FEPK
         FUNN = EXP((-2.+ZLAMB*EQ1_H(JS,KGRID)/2.0)*ZDTK)*FEPK

         ZF1J   = 0.0
         ZF3J   = 0.0

         DO K=1,NEPK-1
            HEPK  = EPK(K+1)-EPK(K)
            RTMP1 = 0.5*HEPK*( FUNP(2*K-1) + FUNP(2*K) )
            RTMP2 = 0.5*HEPK*( FUNN(2*K-1) + FUNN(2*K) )
            ZF1J  =  ZF1J+RTMP1
            ZF3J  =  ZF3J-RTMP2
         ENDDO

         ZF1L   = ZF1L  + ZF1J*ZHLAM*0.5
         ZF3L   = ZF3L  + ZF3J*ZHLAM*0.5
      ENDDO
      ENDDO
       
C     NEXT GO THROUGH TRAPPED PARTICLES       
      DO N   = 1,NLAMK0(JS,KGRID)-1
      ZHLAM  = LAMK0(JS,N+1,KGRID)-LAMK0(JS,N,KGRID)
      DO  J = 0,1
         ZLAMB  = 0.5*((1.+ZW)*LAMK0(JS,N+J,KGRID)
     &            +(1.-ZW)*LAMK0(JS,N+1-J,KGRID))

         ZDTK = ZLAMB*EQ1_H(JS,KGRID)/2.*ZREP**2/ZREG/RE_CONST(6)
         FUNP = EXP(-ZDTK)*FEPK

         ZF2J   = 0.0

         DO K=1,NEPK-1
            HEPK  = EPK(K+1)-EPK(K)
            RTMP1 = 0.5*HEPK*( FUNP(2*K-1) + FUNP(2*K) )
            ZF2J  =  ZF2J+RTMP1
         ENDDO

         ZF2L   = ZF2L  + ZF2J*ZHLAM*0.5
      ENDDO
      ENDDO
       
      
      ZZZI1     = (ZF1L+ZF3L+ZF2L*2.0)*RTMP
      ZZZI3     = 0.0
      ZZZALPHA3 = 1.0

      DEALLOCATE(FUNN,FUNP,ZREP,ZREG)

      RETURN
      END SUBROUTINE ZITYPE6

C===================================================    
C CORE PROCEDURE FOR 
C COMPUTING PARTCILE PHASE SPACE INTERGRALS AT 
C A GIVEN VALUE OF ZZALPHA1, FOR
C I1 = SURFACE AVERAGED DENSITY 
C I3 = SURFACE AVERAGED PRESSURE
C IF0TYPE=7
C===================================================
      SUBROUTINE ZITYPE7(JS,KGRID,KP,ZZALPHA1)
                    
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER   JS,KGRID,KP,KCHECK
      INTEGER   N,J,K
      REAL*8    HEPK,ZHLAM,H3,
     &          ZF1L,ZF3L,ZF1LP,ZF3LP,
     &          ZF1J,ZF3J,ZF1JP,ZF3JP,
     &          RTMP,RTMP1,RTMP2,
     &          ZZALPHA1,
     &          ZF2L,ZF2LP,ZF2J,ZF2JP,
     &          ZLAMB,SMI2ME,AMA2MI,T2 
      REAL*8,  DIMENSION(:),ALLOCATABLE::FUNN,FUNP

      KCHECK  = 0

      ZEPK    = ZEPKO

      ZF1L    = 0.
      ZF3L    = 0.
      ZF1LP   = 0.
      ZF3LP   = 0.

      SMI2ME = SQRT(ESPECIES_M(1)/ESPECIES_M(2))
      AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
      RTMP   = (3.*SQRT(PI)/4.*SMI2ME)**(2./3.)*AMA2MI
      T2     = RTMP*ZZALPHA1
      FEPK  = 1./(ZEPK**1.5+T2**1.5)

      ALLOCATE(FUNN(NEPK2), FUNP(NEPK2))
      FUNN = FEPK*ZEPK**0.5
      FUNP = FUNN*ZEPK

C     FIRST GO THROUGH PASSING PARTICLES
      DO N   = 1,NLAMK1(JS,KGRID)-1
      ZHLAM  = LAMK1(JS,N+1,KGRID)-LAMK1(JS,N,KGRID)
      DO  J = 0,1
         ZLAMB  = 0.5*((1.+ZW)*LAMK1(JS,N+J,KGRID)
     &            +(1.-ZW)*LAMK1(JS,N+1-J,KGRID))
         ZF1J   = 0.0
         ZF3J   = 0.0
         ZF1JP  = 0.0
         ZF3JP  = 0.0

         DO K=1,NEPK-1
         HEPK  = EPK(K+1)-EPK(K)
         RTMP1 = 0.5*HEPK*( FUNN(2*K-1)+FUNN(2*K) )
         ZF1J  =  ZF1J+RTMP1
         ZF3J  =  ZF3J+RTMP1
          
         RTMP1 = 0.5*HEPK*( FUNP(2*K-1)+FUNP(2*K) )
         ZF1JP  = ZF1JP+RTMP1
         ZF3JP  = ZF3JP+RTMP1
         ENDDO

         H3     = ZOMEGABP(JS,2*N-1+J+1,KGRID)
          
         RTMP   = 0.5*ZHLAM*ZLAMB/H3
         ZF1L   = ZF1L  + ZF1J* RTMP
         ZF3L   = ZF3L  + ZF3J* RTMP
         ZF1LP  = ZF1LP + ZF1JP*RTMP
         ZF3LP  = ZF3LP + ZF3JP*RTMP
      ENDDO
      ENDDO
       
C     NEXT GO THROUGH TRAPPED PARTICLES       
      ZF2L   = 0.
      ZF2LP  = 0.
      DO  N  = 1,NLAMK0(JS,KGRID)-1 
      ZHLAM  = LAMK0(JS,N+1,KGRID)-LAMK0(JS,N,KGRID)
      DO  J  = 0,1
         ZLAMB = 0.5*((1.+ZW)*LAMK0(JS,N+J,KGRID)
     &           +(1.-ZW)*LAMK0(JS,N+1-J,KGRID))
         ZF2J  = 0.0
         ZF2JP = 0.0  

         DO K=1,NEPK-1
         HEPK   = EPK(K+1)-EPK(K)
         RTMP1 = 0.5*HEPK*( FUNN(2*K-1)+FUNN(2*K) )
         ZF2J   = ZF2J+RTMP1
          
         RTMP2 = 0.5*HEPK*( FUNP(2*K-1)+FUNP(2*K) )
         ZF2JP  = ZF2JP+RTMP2
         ENDDO

         H3     = ZOMEGABT(JS,2*N-1+J+1,KGRID)

         RTMP   = 0.5*ZHLAM*ZLAMB/H3
         ZF2L   = ZF2L +(ZF2J)*RTMP
         ZF2LP  = ZF2LP+(ZF2JP)*RTMP
      ENDDO
      ENDDO
      
C     0.5 COMES FROM \hat_g FACTOR FOR TRAPPED PARTICLES
      RTMP2     = 0.5 
      ZZZI1     = ZF1L +2.0*RTMP2*ZF2L +ZF3L
      ZZZI3     = ZF1LP+2.0*RTMP2*ZF2LP+ZF3LP
      ZZZALPHA3 = ZZZI3/ZZZI1*2./3.

      DEALLOCATE(FUNN,FUNP)

      RETURN
      END SUBROUTINE ZITYPE7

C===================================================    
C CREATE A FUNCTION FOR COMPUTING ALPHAA3 (AT A GIVEN 
C ALPHAA1) BY PACKAGING THE CORE SUBROUTINE ZZZI1I3 
C APPLIES FOR: IF0TYPE=3,5,7
C===================================================       
      FUNCTION ZALPHA3(JS,KGRID,KP,ZZALPHA1)
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER JS,KGRID,KP
      REAL*8  ZZALPHA1
      REAL*8  ZALPHA3

      IF (ISPECIES_F0(KP).EQ.3) CALL ZZZI1I3(JS,KGRID,KP,ZZALPHA1)
      IF (ISPECIES_F0(KP).EQ.5) CALL ZITYPE5(JS,KGRID,KP,ZZALPHA1)
      IF (ISPECIES_F0(KP).EQ.7) CALL ZITYPE7(JS,KGRID,KP,ZZALPHA1)

      ZALPHA3   = ZZZALPHA3
       
      RETURN
      END FUNCTION ZALPHA3

C===================================================    
C COMPUTE AND STORE INTEGRALS I1 AND I3, AND ALPHAA3
C BY PACKAGING THE CORE SUBROUTINE ZZZI1I3 
C APPLIES FOR: IF0TYPE=3,5,6,7
C===================================================       
      SUBROUTINE ZI1I3(JS,KGRID,KP,ZZALPHA1)
      USE ANISOTROPICM
      USE KINETICM    
      IMPLICIT NONE
      
      INTEGER JS,KGRID,KP
      REAL*8  ZZALPHA1

      IF (ISPECIES_F0(KP).EQ.3) CALL ZZZI1I3(JS,KGRID,KP,ZZALPHA1)
      IF (ISPECIES_F0(KP).EQ.5) CALL ZITYPE5(JS,KGRID,KP,ZZALPHA1)
      IF (ISPECIES_F0(KP).EQ.6) CALL ZITYPE6(JS,KGRID,KP,ZZALPHA1)
      IF (ISPECIES_F0(KP).EQ.7) CALL ZITYPE7(JS,KGRID,KP,ZZALPHA1)

      ZZI1(JS,KGRID,KP)     = ZZZI1
      ZZI3(JS,KGRID,KP)     = ZZZI3
      ALPHAA3(JS,KGRID,KP)  = ZZZALPHA3

      RETURN
      END SUBROUTINE ZI1I3
       

C==================================================
C COFFICIENTS FOR COMPUTING G-FACTORS OF ADIABATIC 
C CONTRIBUTIONS. EVALUATE EQUILIBRIUM QUANTITIES AT
C GAUSSIAN POINTS ALONG CHI BY SPLINE.
C==================================================
      SUBROUTINE ZKEQUILK(JS,KGRID)
      
      USE KINETICM
      USE GLOBALM
      USE DIMENSIM
      USE RCOMDM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER JS,KGRID,J
      
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            ZRQ1(J) = G12L(JS,J)/RJA(JS,J)/BK(JS,J,KGRID)
            ZRQ2(J) = G22L(JS,J)/RJA(JS,J)/BK(JS,J,KGRID)
            ZRQ3(J) = T(JS)/BK(JS,J,KGRID)
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            ZRQ1(J) = G12LM(JS,J)/RJAM(JS,J)/BK(JS,J,KGRID)
            ZRQ2(J) = G22LM(JS,J)/RJAM(JS,J)/BK(JS,J,KGRID)
            ZRQ3(J) = TM(JS)/BK(JS,J,KGRID)
         ENDDO
      ENDIF
      ZRQ1(NCHI+1) = ZRQ1(1)
      ZRQ2(NCHI+1) = ZRQ2(1)
      ZRQ3(NCHI+1) = ZRQ3(1)

      CALL SPLINE1D(ZRQ1K,RCHIK,NCHI2+2,ZRQ1,RCHI,NCHI+1,RCHI2)
      CALL SPLINE1D(ZRQ2K,RCHIK,NCHI2+2,ZRQ2,RCHI,NCHI+1,RCHI2)
      CALL SPLINE1D(ZRQ3K,RCHIK,NCHI2+2,ZRQ3,RCHI,NCHI+1,RCHI2)
      RETURN      
      END SUBROUTINE ZKEQUILK
        
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES. RESULTS STORED IN 
C ZRESU(EPSILONK,1:3,KOPT,KP)   
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C=====================================================
      SUBROUTINE  KDISTRIBF(JS,KGRID,ZZLAMB)

      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP
      REAL*8   ZZLAMB

      DO KP=1,NSPECIES
         CALL KDISTRIBF_TYPE(JS,KGRID,KP,ZZLAMB)
      ENDDO

      RETURN
      END SUBROUTINE KDISTRIBF
                   
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE(JS,KGRID,KP,ZZLAMB)

      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP
      REAL*8   ZZLAMB

      IF (KNUMDISTRIB(KP).EQ.0) ZRESU(:,:,:,KP) = 0.

      IF (ABS(PSPECIES_AP(KP))+ABS(PSPECIES_AT(KP))+ABS(PSPECIES_NP(KP))
     &    +ABS(PSPECIES_NTB(KP))+ABS(PSPECIES_NTD(KP)).GT.0.) THEN

         IF (ISPECIES_F0(KP).EQ.0.AND.KNUMDISTRIB(KP).EQ.0) 
     &      CALL KDISTRIBF_TYPE0(JS,KGRID,KP,ZZLAMB)

         IF ((ISPECIES_F0(KP).EQ.1.OR.ISPECIES_F0(KP).EQ.2).AND.
     &       KNUMDISTRIB(KP).EQ.0) 
     &      CALL KDISTRIBF_TYPE1(JS,KGRID,KP,ZZLAMB)

         IF (ISPECIES_F0(KP).EQ.3) 
     &      CALL KDISTRIBF_TYPE3(JS,KGRID,KP,ZZLAMB)

         IF (ISPECIES_F0(KP).EQ.4) 
     &      CALL KDISTRIBF_TYPE4(JS,KGRID,KP,ZZLAMB)

         IF (ISPECIES_F0(KP).EQ.5) 
     &      CALL KDISTRIBF_TYPE5(JS,KGRID,KP,ZZLAMB)

         IF (ISPECIES_F0(KP).EQ.6) 
     &      CALL KDISTRIBF_TYPE6(JS,KGRID,KP,ZZLAMB)

         IF (ISPECIES_F0(KP).EQ.7) 
     &      CALL KDISTRIBF_TYPE7(JS,KGRID,KP,ZZLAMB)

      ENDIF

      KNUMDISTRIB(KP) = KNUMDISTRIB(KP) + 1

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE
                   
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES FOR IF0TYPE=0 (MAXWELLIAN)
C RESULTS STORED IN ZRESU(EPSILONK,1:3,KOPT,KP)   
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C NOTE THAT FOR THESE OPTIONS, RESULTS INDEPEDENT OF ZZLAMB
C AND BEING THE SAME FOR BOTH PASSING AND TRAPPED PARTICLES
C AND HENCE NEED TO BE CALLED ONLY ONCE (WHEN KNUMDISTRIB=0)
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE0(JS,KGRID,KP,ZZLAMB)
       
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP,KCHECK
       
      REAL*8   ZZLAMB,H1,H2,H3,H6,H7
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::ATMP

      KCHECK  = 0
      ZEPK    = ZEPKN

      ALLOCATE( ATMP(NEPK2,4) )

      H1  = CPSI(JS,KGRID,KP)
      H2  = DCDPSIL(JS,KGRID,KP)
      H3  = DCDPSIL2(JS,KGRID,KP)

      H6  = DEPSALPHADPSI(JS,KGRID,KP)
      H7  = DEPSALPHADPSI2(JS,KGRID,KP)

      FEPK = EXP(-ZEPK)

C     F0: KOPT=0
      ATMP(:,1)       = H1*FEPK

      ZRESU(:,1,0,KP) = ATMP(:,1)
      ZRESU(:,2,0,KP) = ATMP(:,1)
      ZRESU(:,3,0,KP) = ATMP(:,1)
          
C     D F0/D PSI: KOPT=1       
      ATMP(:,2) = H2+(ZEPK-2.5)*H6
      ATMP(:,3) = ATMP(:,1)*ATMP(:,2)
      
      ZRESU(:,1,1,KP) = ATMP(:,3)
      ZRESU(:,2,1,KP) = ATMP(:,3)
      ZRESU(:,3,1,KP) = ATMP(:,3)
       
C     D F0/D (NORMALISED) EPSLONK: KOPT=2
      ATMP(:,4) =-ATMP(:,1)

      ZRESU(:,1,2,KP) = ATMP(:,4)
      ZRESU(:,2,2,KP) = ATMP(:,4)
      ZRESU(:,3,2,KP) = ATMP(:,4)

C     D F0/D LAMBDA: KOPT=3
      ZRESU(:,:,3,KP) = 0.
      
      IF (ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN

C     D^2 F0/D PSI^2: KOPT=4       
      ATMP(:,4) = ATMP(:,2)**2+H3+(ZEPK-2.5)*H7-ZEPK*H6**2
      ATMP(:,4) = ATMP(:,1)*ATMP(:,4)

      ZRESU(:,1,4,KP) = ATMP(:,4)
      ZRESU(:,2,4,KP) = ATMP(:,4)
      ZRESU(:,3,4,KP) = ATMP(:,4)

C     D^2 F0/D EPSILONK/D PSI: KOPT=5
      ATMP(:,4) =-ATMP(:,3)
      ZRESU(:,1,5,KP) = ATMP(:,4)
      ZRESU(:,2,5,KP) = ATMP(:,4)
      ZRESU(:,3,5,KP) = ATMP(:,4)

C     D^2 F0/D LAMBDA/D PSI: KOPT=6
      ZRESU(:,:,6,KP) = 0.

      ENDIF

      DEALLOCATE(ATMP) 

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE0
                   
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES FOR IF0TYPE=1,2
C RESULTS STORED IN ZRESU(EPSILONK,1:3,KOPT,KP)   
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C NOTE THAT FOR THESE OPTIONS, RESULTS INDEPEDENT OF ZZLAMB
C AND BEING THE SAME FOR BOTH PASSING AND TRAPPED PARTICLES
C AND HENCE NEED TO BE CALLED ONLY ONCE (WHEN KNUMDISTRIB=0)
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE1(JS,KGRID,KP,ZZLAMB)
       
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP,KCHECK
       
      REAL*8   ZZLAMB,H1,H2,H3,H4,H5,H6,H7,EPC,EPC2
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::ATMP
      REAL*8,  DIMENSION(:),  ALLOCATABLE::ZEPK2

      KCHECK  = 0
      ZEPK    = ZEPKO

      ALLOCATE( ATMP(NEPK2,4),ZEPK2(NEPK2) )

      EPC  = EPSLONCA(JS,KGRID,KP)
      EPC2 = EPC**1.5
  
      H1  = CPSI(JS,KGRID,KP)
      H2  = DCDPSIL(JS,KGRID,KP)
      H3  = DCDPSIL2(JS,KGRID,KP)

C     H4  = D EPSLONC**1.5/ D PSI = OMEGASCA/H3
      H4  = B0K/OMEGACI0*ESPECIES_TEM(JS,KGRID,KP)
     &      *ESPECIES_Z(1)/ESPECIES_Z(KP)
      H4  = OMEGASCA(JS,KGRID,KP)/H4
      H5  = DEPSCDPSI2(JS,KGRID,KP)

      H6  = DEPSALPHADPSI(JS,KGRID,KP)
      H7  = DEPSALPHADPSI2(JS,KGRID,KP)

      ZEPK2 = ZEPK**1.5
      FEPK  = 1./(ZEPK2+EPC2)
       
C     F0: KOPT=0
      ATMP(:,1)       = H1*FEPK

      ZRESU(:,1,0,KP) = ATMP(:,1)
      ZRESU(:,2,0,KP) = ATMP(:,1)
      ZRESU(:,3,0,KP) = ATMP(:,1)
          
C     D F0/D PSI: KOPT=1       
      ATMP(:,2) = 1.5*ZEPK2*H6-H4
      ATMP(:,3) = H2-2.5*H6+FEPK*ATMP(:,2)
      ATMP(:,4) = ATMP(:,1)*ATMP(:,3)
      
      ZRESU(:,1,1,KP) = ATMP(:,4)
      ZRESU(:,2,1,KP) = ATMP(:,4)
      ZRESU(:,3,1,KP) = ATMP(:,4)
       
C     D F0/D (NORMALISED) EPSLONK: KOPT=2
      ATMP(:,4) =-1.5*ATMP(:,1)*FEPK*SQRT(ZEPK)

      ZRESU(:,1,2,KP) = ATMP(:,4)
      ZRESU(:,2,2,KP) = ATMP(:,4)
      ZRESU(:,3,2,KP) = ATMP(:,4)

C     D F0/D LAMBDA: KOPT=3
      ZRESU(:,:,3,KP) = 0.
      
      IF (ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN

C     D^2 F0/D PSI^2: KOPT=4       
      ATMP(:,4) = ATMP(:,1)*(ATMP(:,3)**2+H3-2.5*H7+
     &            (FEPK*ATMP(:,2))**2+FEPK*(1.5*ZEPK2*
     &            (H7-1.5*H6**2)-H5))
 
      ZRESU(:,1,4,KP) = ATMP(:,4)
      ZRESU(:,2,4,KP) = ATMP(:,4)
      ZRESU(:,3,4,KP) = ATMP(:,4)

C     D^2 F0/D EPSILONK/D PSI: KOPT=5
      ATMP(:,4) =-1.5*ATMP(:,1)*FEPK*SQRT(ZEPK)*(H2-3.*H6+
     &            2.*FEPK*ATMP(:,2))
      ZRESU(:,1,5,KP) = ATMP(:,4)
      ZRESU(:,2,5,KP) = ATMP(:,4)
      ZRESU(:,3,5,KP) = ATMP(:,4)

C     D^2 F0/D LAMBDA/D PSI: KOPT=6
      ZRESU(:,:,6,KP) = 0.

      ENDIF

      DEALLOCATE(ATMP,ZEPK2) 

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE1
                   
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES FOR IF0TYPE=3. 
C RESULTS STORED IN ZRESU(EPSILONK,1:3,KOPT,KP)   
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE3(JS,KGRID,KP,ZZLAMB)
       
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP,KCHECK
       
      REAL*8   ZZLAMB,H1,H2,H3,H4,H5,H6,H7,EPC,EPC2
      REAL*8   ZHMIN,ZHK,ZHK2,ZZALPHA1,DZETAJ
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::ATMP
      REAL*8,  DIMENSION(:),  ALLOCATABLE::ZEPK2

      KCHECK  = 0
      ZEPK    = ZEPKO

      DZETAJ  = DZETA0(KP)

      ALLOCATE( ATMP(NEPK2,7),ZEPK2(NEPK2) )

      ZZALPHA1 = ALPHAA1(JS,KGRID,KP)

      ZHK   = HKMAX(JS,KGRID)
      ZHK2  = ZDHMAXDPSI(JS,KGRID)
      ZHMIN = HKMIN(JS,KGRID)

      EPC  = EPSLONCA(JS,KGRID,KP)
      EPC2 = EPC**1.5
  
      H1  = CPSI(JS,KGRID,KP)
      H2  = DCDPSIL(JS,KGRID,KP)
      H3  = DCDPSIL2(JS,KGRID,KP)

C     H4  = D EPSLONC**1.5/ D PSI = OMEGASCA/H3
      H4  = B0K/OMEGACI0*ESPECIES_TEM(JS,KGRID,KP)
     &      *ESPECIES_Z(1)/ESPECIES_Z(KP)
      H4  = OMEGASCA(JS,KGRID,KP)/H4
      H5  = DEPSCDPSI2(JS,KGRID,KP)

      H6  = DEPSALPHADPSI(JS,KGRID,KP)
      H7  = DEPSALPHADPSI2(JS,KGRID,KP)

      ZEPK2 = ZEPK**1.5
      FEPK  = 1./(ZEPK2+EPC2)
      ZDTK2 = DZETAJ**2-LOG(ZEPK2*(1.+EPC2)*FEPK)/3.
      ZDTK  = SQRT(ZDTK2)
      FEPKC = (ZEPK2-1.)/(EPC2+1.)
       
C     F0: KOPT=0
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,0,0)

      ATMP(:,1)    = H1*FEPK/ZDTK/2./SQRT(PI)
      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
         ZRESU(:,1,0,KP) = ATMP(:,1)*FESS1(:,0,0)
         ZRESU(:,3,0,KP) = ATMP(:,1)*FESS3(:,0,0)       
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
         ZRESU(:,2,0,KP) = ATMP(:,1)*FESS2(:,0,0)
      ENDIF
          
C     D F0/D PSI: KOPT=1       
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,2,0)
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,1,1)

      ATMP(:,2) = H2-2.5*H6+(1.5*ZEPK2-0.25*EPC2/ZDTK2)*
     &            FEPK*H6+(FEPKC/ZDTK2/6.-1.)*FEPK*H4
      ATMP(:,2) = ATMP(:,1)*ATMP(:,2)
      ATMP(:,3) =-ATMP(:,1)*2./ZDTK2
      ATMP(:,4) = ATMP(:,1)*FEPK/ZDTK2**2*(0.5*EPC2*H6
     &           -FEPKC*H4/3.)
      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
         ZRESU(:,1,1,KP) = ATMP(:,2)*FESS1(:,0,0)
     &                    +ATMP(:,3)*FESS1(:,1,1)
     &                    +ATMP(:,4)*FESS1(:,2,0)
         ZRESU(:,3,1,KP) = ATMP(:,2)*FESS3(:,0,0)
     &                    +ATMP(:,3)*FESS3(:,1,1)
     &                    +ATMP(:,4)*FESS3(:,2,0)
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
         ZRESU(:,2,1,KP)  = ATMP(:,2)*FESS2(:,0,0)
     &                     +ATMP(:,3)*FESS2(:,1,1)
     &                     +ATMP(:,4)*FESS2(:,2,0)
      ENDIF
       
C     D F0/D (NORMALISED) EPSLONK: KOPT=2
      ATMP(:,1) = H1*FEPK**2/ZDTK/ZEPK/4./SQRT(PI)
      ATMP(:,2) = ATMP(:,1)*(0.5*EPC2/ZDTK2-3.*ZEPK2)
      ATMP(:,3) =-ATMP(:,1)*EPC2/ZDTK2**2
       
      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
         ZRESU(:,1,2,KP) = ATMP(:,2)*FESS1(:,0,0)+ATMP(:,3)*FESS1(:,2,0)
         ZRESU(:,3,2,KP) = ATMP(:,2)*FESS3(:,0,0)+ATMP(:,3)*FESS3(:,2,0)
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
         ZRESU(:,2,2,KP) = ATMP(:,2)*FESS2(:,0,0)+ATMP(:,3)*FESS2(:,2,0)
      ENDIF
       
C     D F0/D LAMBDA: KOPT=3
      IF (ABS(PSPECIES_AP(KP))+ABS(PSPECIES_AT(KP)).GT.0.) THEN
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,1,0)
         
      ATMP(:,1) = H1*FEPK/ZDTK**3/2./SQRT(PI)/ZHK/SQRT(1.-ZZLAMB/ZHK)  

      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
         ZRESU(:,1,3,KP)  = ATMP(:,1)*FESS1(:,1,0)
         ZRESU(:,3,3,KP)  =-ATMP(:,1)*FESS3(:,1,0)
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
         ZRESU(:,2,3,KP)  = ATMP(:,1)*FESS2(:,1,0)
      ENDIF
      ENDIF
      
      IF (ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN

C     D^2 F0/D PSI^2: KOPT=4       
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,4,0)
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,3,1)
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,0,2)
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,0,3)

C     ATMP1=F1, ATMP2=F2, ATMP3=F3, ATMP4-7 ARE TEMPORARY
      ATMP(:,7) = H1*FEPK/ZDTK/2./SQRT(PI)
      ATMP(:,6) = H2-2.5*H6+(1.5*ZEPK2-0.25*EPC2/ZDTK2)*
     &            FEPK*H6+(FEPKC/ZDTK2/6.-1.)*FEPK*H4
      ATMP(:,1) = ATMP(:,7)*ATMP(:,6)
      ATMP(:,2) =-ATMP(:,7)*2./ZDTK2
      ATMP(:,3) = ATMP(:,7)*FEPK/ZDTK2**2*(0.5*EPC2*H6-FEPKC*H4/3.)

C     ATMP4 IS COEFFICIENT FOR F_00, ATMP5-7 ARE TEMPORARY
      ATMP(:,5) = H3-2.5*H7+(-2.25*ZEPK2*H6-0.25*H4/ZDTK2+0.25*EPC2/
     &            ZDTK2**2*(0.5*EPC2*H6-FEPKC*H4/3.)*FEPK)*FEPK*H6+
     &            (1.5*ZEPK2-0.25*EPC2/ZDTK2)*FEPK*((1.5*ZEPK2*H6-H4)*
     &            FEPK*H6+H7)-((1.5*ZEPK2*H6+FEPKC*H4)/(EPC2+1.)/ZDTK2+
     &            FEPKC*FEPK*(0.5*EPC2*H6-FEPKC*H4/3.)/ZDTK2**2)/6.*
     &            FEPK*H4+(FEPKC/ZDTK2/6.-1.)*FEPK*((1.5*ZEPK2*H6-H4)*
     &            FEPK*H4+H5)
      ATMP(:,4) = ATMP(:,7)*(ATMP(:,6)**2+ATMP(:,5))

C     ATMP5 IS COEFFICIENT FOR F_20, ATMP6-7 ARE TEMPORARY
      ATMP(:,6) = ATMP(:,3)*(H2-2.5*H6+2.*FEPK*(1.5*ZEPK2*H6-H4)-
     &            2.5*FEPK/ZDTK2*(0.5*EPC2*H6-FEPKC*H4/3.))+
     &            ATMP(:,7)*FEPK/ZDTK2**2*(0.5*H6*H4+0.5*EPC2*H7+
     &            (1.5*ZEPK2*H6+FEPKC*H4)*H4/3./(EPC2+1.)-FEPKC*H5/3.)
      ATMP(:,5) = ATMP(:,6)+ATMP(:,1)*FEPK/ZDTK2**2
     &            *(0.5*EPC2*H6-FEPKC*H4/3.)
         
C     ATMP6 IS COEFFICIENT FOR F_11, ATMP7 IS TEMPORARY
      ATMP(:,7) = H2-2.5*H6+FEPK*(1.5*ZEPK2*H6-H4)-1.5*FEPK/ZDTK2
     &            *(0.5*EPC2*H6-FEPKC*H4/3.)
      ATMP(:,6) = ATMP(:,2)*ATMP(:,7)+2.*ATMP(:,3)-2.*ATMP(:,1)/ZDTK2

C     ATMP7 IS COEFFICIENT FOR F_31
      ATMP(:,7) =-2.*ATMP(:,3)/ZDTK2+ATMP(:,2)*FEPK/ZDTK2**2
     &            *(0.5*EPC2*H6-FEPKC*H4/3.)

C     ATMP1 IS COEFFICIENT FOR F_03
      ATMP(:,1) =-2.*ATMP(:,2)/ZDTK2

C     ATMP3 IS COEFFICIENT FOR F_40
      ATMP(:,3) = ATMP(:,3)*FEPK/ZDTK2**2*(0.5*EPC2*H6-FEPKC*H4/3.)

      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
         ZRESU(:,1,4,KP) = ATMP(:,4)*FESS1(:,0,0)
     &                    +ATMP(:,6)*FESS1(:,1,1)
     &                    +ATMP(:,5)*FESS1(:,2,0)
     &                    +ATMP(:,7)*FESS1(:,3,1)
     &                    +ATMP(:,1)*FESS1(:,0,3)
     &                    +ATMP(:,3)*FESS1(:,4,0)
     &                    +ATMP(:,2)*FESS1(:,0,2)
         ZRESU(:,3,4,KP) = ATMP(:,4)*FESS3(:,0,0)
     &                    +ATMP(:,6)*FESS3(:,1,1)
     &                    +ATMP(:,5)*FESS3(:,2,0)
     &                    +ATMP(:,7)*FESS3(:,3,1)
     &                    +ATMP(:,1)*FESS3(:,0,3)
     &                    +ATMP(:,3)*FESS3(:,4,0)
     &                    +ATMP(:,2)*FESS3(:,0,2)
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
         ZRESU(:,2,4,KP) = ATMP(:,4)*FESS2(:,0,0)
     &                    +ATMP(:,6)*FESS2(:,1,1)
     &                    +ATMP(:,5)*FESS2(:,2,0)
     &                    +ATMP(:,7)*FESS2(:,3,1)
     &                    +ATMP(:,1)*FESS2(:,0,3)
     &                    +ATMP(:,3)*FESS2(:,4,0)
     &                    +ATMP(:,2)*FESS2(:,0,2)
      ENDIF
       
C     D^2 F0/D EPSILONK/D PSI: KOPT=5
C     ATMP4=F4, ATMP5=F5, OTHERS ARE TEMPORARY
      ATMP(:,1) = H1*FEPK**2/ZDTK/ZEPK/4./SQRT(PI)
      ATMP(:,4) = ATMP(:,1)*(0.5*EPC2/ZDTK2-3.*ZEPK2)
      ATMP(:,5) =-ATMP(:,1)*EPC2/ZDTK2**2

C     ATMP1 IS COEFFICIENT FOR F_00
      ATMP(:,6) = H2-2.5*H6+(3.*ZEPK2-0.25*EPC2/ZDTK2)*FEPK*H6
     &            +(FEPKC/ZDTK2/6.-2.)*FEPK*H4
      ATMP(:,7) = 0.125*H1*FEPK**2/ZDTK/ZEPK*( EPC2/ZDTK2
     &            *(1.-0.5*EPC2/ZDTK2*FEPK)*H6+3.*ZEPK2*H6
     &            +(1.+EPC2/ZDTK2*FEPK*FEPKC/3.)/ZDTK2*H4 )/SQRT(PI)
      ATMP(:,1) = ATMP(:,4)*ATMP(:,6)+ATMP(:,7)

C     ATMP2 IS COEFFICIENT FOR F_20
      ATMP(:,6) = H2-2.5*H6+(1.+3.*FEPK*ZEPK2-1.25*EPC2/ZDTK2*FEPK)
     &            *H6+(1./EPC2-2.*FEPK+5.*FEPK*FEPKC/ZDTK2/6.)*H4
      ATMP(:,7) = (0.5*EPC2*H6-FEPKC*H4/3.)*FEPK/ZDTK2**2
      ATMP(:,2) = ATMP(:,4)*ATMP(:,7)+ATMP(:,5)*ATMP(:,6)

C     ATMP3 IS COEFFICIENT FOR F_11
      ATMP(:,3) = 2.*(ATMP(:,5)-ATMP(:,4)/ZDTK2)

C     ATMP6 IS COEFFICIENT FOR F_40
      ATMP(:,6) = ATMP(:,5)*ATMP(:,7)

C     ATMP7 IS COEFFICIENT FOR F_31
      ATMP(:,7) =-2.*ATMP(:,5)/ZDTK2

      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
         ZRESU(:,1,5,KP) = ATMP(:,1)*FESS1(:,0,0)
     &                    +ATMP(:,2)*FESS1(:,2,0)
     &                    +ATMP(:,3)*FESS1(:,1,1)
     &                    +ATMP(:,6)*FESS1(:,4,0)
     &                    +ATMP(:,7)*FESS1(:,3,1)
         ZRESU(:,3,5,KP) = ATMP(:,1)*FESS3(:,0,0)
     &                    +ATMP(:,2)*FESS3(:,2,0)
     &                    +ATMP(:,3)*FESS3(:,1,1)
     &                    +ATMP(:,6)*FESS3(:,4,0)
     &                    +ATMP(:,7)*FESS3(:,3,1)
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
         ZRESU(:,2,5,KP) = ATMP(:,1)*FESS2(:,0,0)
     &                    +ATMP(:,2)*FESS2(:,2,0)
     &                    +ATMP(:,3)*FESS2(:,1,1)
     &                    +ATMP(:,6)*FESS2(:,4,0)
     &                    +ATMP(:,7)*FESS2(:,3,1)
      ENDIF

      ENDIF

C     D^2 F0/D LAMBDA/D PSI: KOPT=6
      IF (ABS(PSPECIES_AP(KP))+ABS(PSPECIES_AT(KP)).GT.0..AND.
     &    ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,0,1)
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,3,0)
      CALL KF0_TYPE3(JS,KGRID,KP,ZZALPHA1,ZZLAMB,2,1)

      ATMP(:,1) = H1*FEPK/ZDTK**3/2./SQRT(PI)/ZHK/SQRT(1.-ZZLAMB/ZHK)
      ATMP(:,2) = H2-2.5*H6+1.5*(ZEPK2-0.5*EPC2/ZDTK2)*FEPK*H6
     &            +(0.5*FEPKC/ZDTK2-1.)*FEPK*H4-(1.+0.5*ZZLAMB/ZHK
     &            /(1.-ZZLAMB/ZHK))*ZHK2/ZHK
      ATMP(:,3) = (0.5*EPC2*H6-FEPKC*H4/3.)*FEPK/ZDTK2**2
      ATMP(:,4) =-2./ZDTK2

      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
         ZRESU(:,1,6,KP)  = ATMP(:,1)*(
     &                      ATMP(:,2)*FESS1(:,1,0)+FESS1(:,0,1)+
     &                      ATMP(:,3)*FESS1(:,3,0)+
     &                      ATMP(:,4)*FESS1(:,2,1) )
         ZRESU(:,3,6,KP)  =-ATMP(:,1)*(
     &                      ATMP(:,2)*FESS3(:,1,0)+FESS3(:,0,1)+
     &                      ATMP(:,3)*FESS3(:,3,0)+
     &                      ATMP(:,4)*FESS3(:,2,1) )
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
         ZRESU(:,2,6,KP)  = ATMP(:,1)*(
     &                      ATMP(:,2)*FESS2(:,1,0)+FESS2(:,0,1)+
     &                      ATMP(:,3)*FESS2(:,3,0)+
     &                      ATMP(:,4)*FESS2(:,2,1) )
      ENDIF
      ENDIF

      DEALLOCATE(ATMP,ZEPK2) 

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         H1 = SQRT(1.-ZZLAMB/ZHK)
         H2 = SQRT(1.-ZHMIN/ZHK)
         IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN) THEN
            WRITE(*,110)  H1,FESS1(NEPK,0,0),H2,ZDTK(NEPK)
            WRITE(*,110) -H1,FESS3(NEPK,0,0),H2,ZDTK(NEPK)
         ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK) THEN
            WRITE(*,110)  H1,FESS2(NEPK,0,0),H2,ZDTK(NEPK)
            WRITE(*,110) -H1,FESS2(NEPK,0,0),H2,ZDTK(NEPK)
         ENDIF
 110     FORMAT(4(E13.5,1X))
      ENDIF

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE3
                   
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES FOR IF0TYPE=4 (ICRH MODEL)
C RESULTS STORED IN ZRESU(EPSILONK,1:3,KOPT,KP)   
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE4(JS,KGRID,KP,ZZLAMB)
       
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP,KCHECK,J,K,KPT
       
      REAL*8   ZZLAMB,H1,H2,H3,H4,H6,H7,H9,H11,H13,H14
      REAL*8   HKM,HKMP,HKMPP,HS,HSP,HSL,HSPP,HSLP,
     &         HU,HUP,HUL,HUPP,HULP,H10,H10P,H10L,H10PP,H10LP
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::ATMP
      REAL*8   REPS
      PARAMETER (REPS=1.0E-80)

      KCHECK  = 0
      ZEPK    = ZEPKN

      ALLOCATE( ATMP(NEPK2,8) )

      H1  = CPSI(JS,KGRID,KP)
      H2  = DCDPSIL(JS,KGRID,KP)
      H3  = DCDPSIL2(JS,KGRID,KP)
      H6  = DEPSALPHADPSI(JS,KGRID,KP)
      H7  = DEPSALPHADPSI2(JS,KGRID,KP)
      H9  = ZZLAMB/HTYPE4C+RTYPE4(JS,KGRID)*ABS(1.-ZZLAMB/HTYPE4C)
      H14 = HTYPE4C

      HKM   = HKMIN     (JS,KGRID)
      HKMP  = HMIN3DPSI (JS,KGRID)
      HKMPP = HMIN3DPSI2(JS,KGRID)

      KPT = 0
      IF (ZZLAMB.LT.HKM) KPT=1

      IF (KSMTYPE4.EQ.1.AND.KPT.EQ.1) THEN 
         H4   = (ZZLAMB/HKM)**2
         HS   = (1.-H4)**2
         HSP  = -4.*ZZLAMB**2/HKM*(1.-H4)*HKMP
         HSL  = -4.*ZZLAMB/HKM**2*(1.-H4)
         HSPP = -4.*ZZLAMB**2*((1.-3.*H4)*HKMP**2+(1.-H4)*HKMPP/HKM)
         HSLP = -8.*ZZLAMB/HKM*(1.-2.*H4)*HKMP
      ELSEIF (KSMTYPE4.EQ.2.AND.KPT.EQ.1) THEN
         IF (HKM.GE.1.) STOP 'HMIN>1, ADJUST B0K'
         H4   = (1.-ZZLAMB/HKM)*(1.-2.*ZZLAMB+ZZLAMB/HKM)
         HU   = (1.-ZZLAMB/HKM)/(1.-2.*ZZLAMB+ZZLAMB/HKM)
         HUP  = -2.*ZZLAMB*(1.-ZZLAMB)/H4*HKMP
         HUL  = 2.*(1.-1./HKM)/H4
         HUPP = -HUP*HUL*ZZLAMB**2*HKMP-2.*ZZLAMB*(1.-ZZLAMB)/H4*HKMPP
         HULP = -2.*((1.-ZZLAMB)**2+(ZZLAMB-ZZLAMB/HKM)**2)/H4**2*HKMP
         H4   = 4./3.*LOG(HU/1.4)
         HS   = 0.8*EXP(-(0.5*H4)**2)/HU
         HSP  = -(1.+H4)*HUP*HS
         HSL  = -(1.+H4)*HUL*HS
         HSPP = (((1.+H4)**2-4./3.)*HUP**2-(1.+H4)*HUPP)*HS
         HSLP = (((1.+H4)**2-4./3.)*HUP*HUL-(1.+H4)*HULP)*HS
      ENDIF

      IF (KPT.EQ.1) THEN
      H10   = STYPE4(JS,KGRID)*HS
      H10P  = STYPE4DPSI(JS,KGRID)*HS+STYPE4(JS,KGRID)*HSP
      H10L  = STYPE4(JS,KGRID)*HSL
      H10PP = STYPE4DPSI2(JS,KGRID)*HS+2.*STYPE4DPSI(JS,KGRID)*HSP+
     &        STYPE4(JS,KGRID)*HSPP
      H10LP = STYPE4DPSI(JS,KGRID)*HSL+STYPE4(JS,KGRID)*HSLP
      H11   = 1.+H10
      H13   = 1.-H10
      ENDIF

      FEPK  = EXP(-ZEPK*H9)
       
C     F0: KOPT=0
      ATMP(:,1)       = H1*FEPK

      IF (KPT.EQ.1) ZRESU(:,1,0,KP) = ATMP(:,1)*H11
      IF (KPT.EQ.0) ZRESU(:,2,0,KP) = ATMP(:,1)
      IF (KPT.EQ.1) ZRESU(:,3,0,KP) = ATMP(:,1)*H13
          
C     D F0/D PSI: KOPT=1 
      ATMP(:,3) = ZEPK*(RTYPE4DPSI(JS,KGRID)*ABS(1.-ZZLAMB/H14)-H9*H6)

      IF (KPT.EQ.1) THEN
      ATMP(:,7) = ATMP(:,1)*(H2-2.5*H6+H10P/H11-ATMP(:,3))
      ZRESU(:,1,1,KP) = ATMP(:,7)*H11
      ATMP(:,7) = ATMP(:,1)*(H2-2.5*H6-H10P/H13-ATMP(:,3))
      ZRESU(:,3,1,KP) = ATMP(:,7)*H13
      ELSE
      ATMP(:,7) = ATMP(:,1)*(H2-2.5*H6         -ATMP(:,3))
      ZRESU(:,2,1,KP) = ATMP(:,7)
      ENDIF
       
C     D F0/D (NORMALISED) EPSLONK: KOPT=2
      ATMP(:,7) =-ATMP(:,1)*H9

      IF (KPT.EQ.1) ZRESU(:,1,2,KP) = ATMP(:,7)*H11
      IF (KPT.EQ.0) ZRESU(:,2,2,KP) = ATMP(:,7)
      IF (KPT.EQ.1) ZRESU(:,3,2,KP) = ATMP(:,7)*H13

C     D F0/D EPSILONK CANNOT BE ZERO
      DO K=1,3
      DO J=1,NEPK2
         IF (ABS(ZRESU(J,K,2,KP)).LT.REPS) ZRESU(J,K,2,KP) = REPS
      ENDDO
      ENDDO

C     D F0/D LAMBDA: KOPT=3
      IF (ABS(PSPECIES_AP(KP))+ABS(PSPECIES_AT(KP)).GT.0.) THEN
      IF (ZZLAMB.GE.H14) ATMP(:,5) = ZEPK*(1.+RTYPE4(JS,KGRID))/H14
      IF (ZZLAMB.LE.H14) ATMP(:,5) = ZEPK*(1.-RTYPE4(JS,KGRID))/H14

      IF (KPT.EQ.1) THEN
      ATMP(:,7) = ATMP(:,1)*( H10L/H11-ATMP(:,5))
      ZRESU(:,1,3,KP) = ATMP(:,7)*H11
      ATMP(:,7) = ATMP(:,1)*(-H10L/H13-ATMP(:,5))
      ZRESU(:,3,3,KP) = ATMP(:,7)*H13
      ELSE
      ATMP(:,7) = ATMP(:,1)*(         -ATMP(:,5))
      ZRESU(:,2,3,KP) = ATMP(:,7)
      ENDIF
      ENDIF
      
      IF (ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN

C     D^2 F0/D PSI^2: KOPT=4       
      ATMP(:,6) =-ATMP(:,3)*H6-ZEPK*H9*H7+ZEPK*(RTYPE4DPSI2(JS,KGRID)-
     &            H6*RTYPE4DPSI(JS,KGRID))*ABS(1.-ZZLAMB/H14)

      IF (KPT.EQ.1) THEN
      ATMP(:,7) = ATMP(:,1)*((H2-2.5*H6+H10P/H11-ATMP(:,3))**2
     &            +H3-2.5*H7+( H10PP*H11-H10P**2)/H11**2-
     &            ATMP(:,6))
      ZRESU(:,1,4,KP) = ATMP(:,7)*H11
      ATMP(:,7) = ATMP(:,1)*((H2-2.5*H6-H10P/H13-ATMP(:,3))**2
     &            +H3-2.5*H7+(-H10PP*H13-H10P**2)/H13**2-
     &            ATMP(:,6))
      ZRESU(:,3,4,KP) = ATMP(:,7)*H13
      ELSE
      ATMP(:,7) = ATMP(:,1)*((H2-2.5*H6         -ATMP(:,3))**2
     &            +H3-2.5*H7-
     &            ATMP(:,6))
      ZRESU(:,2,4,KP) = ATMP(:,7)
      ENDIF

C     D^2 F0/D EPSILONK/D PSI: KOPT=5
      IF (KPT.EQ.1) THEN
      ATMP(:,6) = H2-2.5*H6+H10P/H11-ATMP(:,3)
      ATMP(:,7) =-ATMP(:,1)*(RTYPE4DPSI(JS,KGRID)*ABS(1.-ZZLAMB/H14)+
     &            H9*ATMP(:,6))
      ZRESU(:,1,5,KP) = ATMP(:,7)*H11

      ATMP(:,6) = H2-2.5*H6-H10P/H13-ATMP(:,3)
      ATMP(:,7) =-ATMP(:,1)*(RTYPE4DPSI(JS,KGRID)*ABS(1.-ZZLAMB/H14)+
     &            H9*ATMP(:,6))
      ZRESU(:,3,5,KP) = ATMP(:,7)*H13
      ELSE
      ATMP(:,6) = H2-2.5*H6         -ATMP(:,3)
      ATMP(:,7) =-ATMP(:,1)*(RTYPE4DPSI(JS,KGRID)*ABS(1.-ZZLAMB/H14)+
     &            H9*ATMP(:,6))
      ZRESU(:,2,5,KP) = ATMP(:,7)
      ENDIF

C     D^2 F0/D LAMBDA/D PSI: KOPT=6
      IF (ABS(PSPECIES_AP(KP))+ABS(PSPECIES_AT(KP)).GT.0.) THEN
      IF (ZZLAMB.GE.H14) THEN
         ATMP(:,6) = ZEPK*(1.+RTYPE4(JS,KGRID))/H14
         ATMP(:,7) = ZEPK*( RTYPE4DPSI(JS,KGRID)-H6*(1.+
     &               RTYPE4(JS,KGRID)))/H14
      ELSE
         ATMP(:,6) = ZEPK*(1.-RTYPE4(JS,KGRID))/H14
         ATMP(:,7) = ZEPK*(-RTYPE4DPSI(JS,KGRID)-H6*(1.-
     &               RTYPE4(JS,KGRID)))/H14
      ENDIF

      IF (KPT.EQ.1) THEN
      ATMP(:,8) = ATMP(:,1)*((H2-2.5*H6+H10P/H11-ATMP(:,3))*
     &            ( H10L/H11-ATMP(:,6))+( H10LP*H11-
     &            H10P*H10L)/H11**2-ATMP(:,7))
      ZRESU(:,1,6,KP) = ATMP(:,8)*H11

      ATMP(:,8) = ATMP(:,1)*((H2-2.5*H6-H10P/H13-ATMP(:,3))*
     &            (-H10L/H13-ATMP(:,6))+(-H10LP*H13-
     &            H10P*H10L)/H13**2-ATMP(:,7))
      ZRESU(:,3,6,KP) = ATMP(:,8)*H13
      ELSE
      ATMP(:,8) = ATMP(:,1)*((H2-2.5*H6          -ATMP(:,3))*
     &            (         -ATMP(:,6))-ATMP(:,7))
      ZRESU(:,2,6,KP) = ATMP(:,8)
      ENDIF
      ENDIF

      ENDIF

      DEALLOCATE(ATMP) 

      IF (KCHECK.EQ.1.AND.JS.EQ.JS0.AND.KGRID.EQ.1) THEN
         J = INT(NLAMK0(JS,KGRID)/2.)
         IF (ABS(LAMM(J+1)-ZZLAMB).LT.1.0E-13) THEN
            WRITE(*,121) KP,ZZLAMB
            DO J=1,NEPK2
               WRITE(*,122) ZEPK(J),ZRESU(J,2,0,KP),ZRESU(J,2,1,KP),
     &                              ZRESU(J,2,2,KP),ZRESU(J,2,3,KP)
            ENDDO
         ENDIF
 121     FORMAT('TYPE4:KP,ZZLAMB=',I2,1X,E13.5)
 122     FORMAT(5(E13.5,1X))
      ENDIF

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE4
                   
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES FOR IF0TYPE=5
C RESULTS STORED IN ZRESU(EPSILONK,1:3,KOPT,KP)   
C HAT_F01 = CPSI
C     F02 = REL
C     F03 = REE
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE5(JS,KGRID,KP,ZZLAMB)
       
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP,KCHECK
       
      REAL*8   ZZLAMB,H1,H2,H3,H4,H5,H6,H7,EPC,EPC2,REL,DRELDLL
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::ATMP
      REAL*8,  DIMENSION(:),  ALLOCATABLE::ZEPK2,YREE,YREL,
     &                                     REE,DREEDEL,DREEDEL2

      KCHECK  = 0
      ZEPK    = ZEPKO

      ALLOCATE( ATMP(NEPK2,5),ZEPK2(NEPK2), 
     &          REE(NEPK2),DREEDEL(NEPK2),DREEDEL2(NEPK2),
     &          YREE(NEPK),YREL(NLAMK) )

      CALL SPLINE1D(REE,ZEPK,NEPK2,ESPECIES_REE(:,2),
     &              ESPECIES_REE(:,1),NEPK,YREE)
      CALL DFFFDX(DREEDEL,ZEPK,NEPK2,ESPECIES_REE(:,2),
     &            ESPECIES_REE(:,1),NEPK,1)
      CALL DFFFDX(DREEDEL2,ZEPK,NEPK2,DREEDEL,ZEPK,NEPK2,0)

      H1  = ZZLAMB/HKMAX(JS,KGRID)
      CALL SPLINE1D(REL,H1,1,ESPECIES_REL(:,2),
     &              ESPECIES_REL(:,1),NLAMK,YREL)
      CALL DFFFDX(DRELDLL,H1,1,ESPECIES_REL(:,2),
     &            ESPECIES_REL(:,1),NLAMK,1)

      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)                            !=P0
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI  !=N0
      H6 = 8.1872E-14                                      !=me*c^2
      H5 = H1/H2/H6                                        !=T0/(me*c^2)
      H7 = ESPECIES_TEM(JS,KGRID,2)/ALPHAA1(JS,KGRID,KP)*H5!=EPS_A/(me*c^2)

      H1  = CPSI(JS,KGRID,KP)
      H2  = DCDPSIL(JS,KGRID,KP)
      H3  = DCDPSIL2(JS,KGRID,KP)

      H4  = EPSALPHA(JS,KGRID,KP)
      H5  = DEPSALPHADPSI(JS,KGRID,KP)
      H6  = DEPSALPHADPSI2(JS,KGRID,KP)

      FEPK  = ZEPK*H7 + 1.0
       
C     F0: KOPT=0
C     ATMP(:,1)       = H1*REL*REE*FEPK**(-1.5)
      ATMP(:,1)       = H1*REL*REE

      ZRESU(:,1,0,KP) = ATMP(:,1)
      ZRESU(:,2,0,KP) = ATMP(:,1)
      ZRESU(:,3,0,KP) = ATMP(:,1)
          
C     D F0/D PSI: KOPT=1       
      ATMP(:,2) = H2 - DREEDEL*ZEPK*H5/H4       !FACTOR A
      ATMP(:,5) = ATMP(:,1)*ATMP(:,2)
      
      ZRESU(:,1,1,KP) = ATMP(:,5)
      ZRESU(:,2,1,KP) = ATMP(:,5)
      ZRESU(:,3,1,KP) = ATMP(:,5)
       
C     D F0/D (NORMALISED) EPSLONK: KOPT=2
C     ATMP(:,3) =-1.5*H7/FEPK + DREEDEL    !FACTOR B
      ATMP(:,3) = DREEDEL    !FACTOR B
      ATMP(:,5) = ATMP(:,1)*ATMP(:,3)

      ZRESU(:,1,2,KP) = ATMP(:,5)
      ZRESU(:,2,2,KP) = ATMP(:,5)
      ZRESU(:,3,2,KP) = ATMP(:,5)

C     D F0/D LAMBDA: KOPT=3
      ATMP(:,4) = DRELDLL/HKMAX(JS,KGRID)  !FACTOR C
      ATMP(:,5) = ATMP(:,1)*ATMP(:,4)

      ZRESU(:,1,3,KP) = ATMP(:,5)
      ZRESU(:,2,3,KP) = ATMP(:,5)
      ZRESU(:,3,3,KP) = ATMP(:,5)
      
      IF (ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN

C     D^2 F0/D PSI^2: KOPT=4       
      ATMP(:,5) = ATMP(:,1)*(ATMP(:,2)**2+H3+
     &                       DREEDEL2*ZEPK**2*(H5/H4)**2+
     &                       DREEDEL*ZEPK*(2.0*(H5/H4)**2-H6/H4))
 
      ZRESU(:,1,4,KP) = ATMP(:,5)
      ZRESU(:,2,4,KP) = ATMP(:,5)
      ZRESU(:,3,4,KP) = ATMP(:,5)

C     D^2 F0/D EPSILONK/D PSI: KOPT=5
C     ATMP(:,5) = ATMP(:,1)*(ATMP(:,2)*ATMP(:,3)-1.5/FEPK*H7*H5/H4-
C    &                       DREEDEL2*ZEPK*H5/H4)
      ATMP(:,5) = ATMP(:,1)*(ATMP(:,2)*ATMP(:,3)-DREEDEL2*ZEPK*H5/H4)

      ZRESU(:,1,5,KP) = ATMP(:,5)
      ZRESU(:,2,5,KP) = ATMP(:,5)
      ZRESU(:,3,5,KP) = ATMP(:,5)

C     D^2 F0/D LAMBDA/D PSI: KOPT=6
      ATMP(:,5) = ATMP(:,1)*ATMP(:,2)*ATMP(:,4)

      ZRESU(:,1,6,KP) = ATMP(:,5)
      ZRESU(:,2,6,KP) = ATMP(:,5)
      ZRESU(:,3,6,KP) = ATMP(:,5)

      ENDIF

      DEALLOCATE(ATMP,ZEPK2,REE,DREEDEL,DREEDEL2,YREE,YREL) 

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE5
                   
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES FOR IF0TYPE=6
C RESULTS STORED IN ZRESU(EPSILONK,1:3,KOPT,KP)   
C HAT_F01 = CPSI
C     F02 = REL
C     F03 = REE
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE6(JS,KGRID,KP,ZZLAMB)
       
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      USE REORBITM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP,KCHECK
       
      REAL*8   ZZLAMB,H1,H2,H3,H4
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::ATMP
      REAL*8,  DIMENSION(:),  ALLOCATABLE::ZEPK2,ZREP,ZREG

      KCHECK  = 0

      ALLOCATE( ATMP(NEPK2,5),ZEPK2(NEPK2), 
     &          ZREP(NEPK2),ZREG(NEPK2) )

      ZEPK  = ZEPKO
      ZREP  = 1.0/SQRT(1.0/ZEPK-1.0)
      ZREG  = 1.0/SQRT(1.0-ZEPK)
      ZREP(NEPK2) = ZREP(NEPK2-1)
      ZREG(NEPK2) = ZREG(NEPK2-1)

      H1  = CPSI(JS,KGRID,KP)
      H2  = DCDPSIL(JS,KGRID,KP)
      H3  = DCDPSIL2(JS,KGRID,KP)
      H4  = H2**2+H3

      ZEPK2     = ZREP**2/ZREG/RE_CONST(6)
      ATMP(:,1) = EXP(-ZZLAMB*EQ1_H(JS,KGRID)/2.*ZEPK2)
      ATMP(:,2) = EXP((-2.+ZZLAMB*EQ1_H(JS,KGRID)/2.0)*ZEPK2)
      ZEPK2     = EXP(-(ZREP-RE_PMAX)**2/RE_CONST(7))
      ATMP(:,1) = ATMP(:,1)*ZEPK2*H1 
      ATMP(:,2) = ATMP(:,2)*ZEPK2*H1 

C     F0: KOPT=0
      ZRESU(:,1,0,KP) = ATMP(:,1)
      ZRESU(:,2,0,KP) = ATMP(:,1)
      ZRESU(:,3,0,KP) = ATMP(:,2)
          
C     D F0/D PSI: KOPT=1       
      ZRESU(:,1,1,KP) = ATMP(:,1)*H2
      ZRESU(:,2,1,KP) = ATMP(:,1)*H2
      ZRESU(:,3,1,KP) = ATMP(:,2)*H2
       
C     D F0/D (NORMALISED) EPSLONK: KOPT=2
      ZEPK2     = (ZREG+ZREG**3)/2.0/RE_CONST(6) 
      ATMP(:,3) = -ZZLAMB*EQ1_H(JS,KGRID)/2.*ZEPK2
      ATMP(:,4) = (-2.+ZZLAMB*EQ1_H(JS,KGRID)/2.0)*ZEPK2
      ZEPK2     =-(ZREP-RE_PMAX)*ZREG**4/ZREP/RE_CONST(7)
      ATMP(:,3) = ATMP(:,3)+ZEPK2
      ATMP(:,4) = ATMP(:,4)+ZEPK2

      ZRESU(:,1,2,KP) = ATMP(:,1)*ATMP(:,3)
      ZRESU(:,2,2,KP) = ATMP(:,1)*ATMP(:,3)
      ZRESU(:,3,2,KP) = ATMP(:,2)*ATMP(:,4)

C     D F0/D LAMBDA: KOPT=3
      ATMP(:,5) = ZREP**2*EQ1_H(JS,KGRID)/2.0/RE_CONST(6)/ZREG

      ZRESU(:,1,3,KP) = ATMP(:,1)*ATMP(:,5)
      ZRESU(:,2,3,KP) = ATMP(:,1)*ATMP(:,5)
      ZRESU(:,3,3,KP) =-ATMP(:,2)*ATMP(:,5)
      
      IF (ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN

C     D^2 F0/D PSI^2: KOPT=4       
      ZRESU(:,1,4,KP) = ATMP(:,1)*H4
      ZRESU(:,2,4,KP) = ATMP(:,1)*H4
      ZRESU(:,3,4,KP) = ATMP(:,2)*H4

C     D^2 F0/D EPSILONK/D PSI: KOPT=5
      ZEPK2     = ZZLAMB*(ZREG+ZREG**3)*EQD1_H(JS,KGRID)/2.0/RE_CONST(6)
      ATMP(:,3) = H2*ATMP(:,3) - ZEPK2
      ATMP(:,4) = H2*ATMP(:,4) + ZEPK2

      ZRESU(:,1,5,KP) = ATMP(:,1)*ATMP(:,3)
      ZRESU(:,2,5,KP) = ATMP(:,1)*ATMP(:,3)
      ZRESU(:,3,5,KP) = ATMP(:,2)*ATMP(:,4)

C     D^2 F0/D LAMBDA/D PSI: KOPT=6
      ZEPK2     = ZREP**2*EQD1_H(JS,KGRID)/2.0/RE_CONST(6)/ZREG
      ATMP(:,5) = H2*ATMP(:,5) - ZEPK2

      ZRESU(:,1,6,KP) = ATMP(:,1)*ATMP(:,5)
      ZRESU(:,2,6,KP) = ATMP(:,1)*ATMP(:,5)
      ZRESU(:,3,6,KP) =-ATMP(:,2)*ATMP(:,5)

      ENDIF

      DEALLOCATE(ATMP,ZEPK2,ZREP,ZREG)

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE6
                   
C===================================================== 
C EVALUATE F0 AND ITS DERIVATIVES FOR IF0TYPE=7
C RESULTS STORED IN ZRESU(EPSILONK,1:3,KOPT,KP)   
C   KOPT=0   F0                     
C   KOPT=1   D F0/D PSI                    
C   KOPT=2   D F0/D EPSLONK                
C   KOPT=3   D F0/D LAMBDA
C   KOPT=4   D^2 F0/D PSI^2
C   KOPT=5   D^2 F0/D EPSILONK/D PSI
C   KOPT=6   D^2 F0/D LAMBDA/D PSI
C=====================================================
      SUBROUTINE  KDISTRIBF_TYPE7(JS,KGRID,KP,ZZLAMB)
       
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
       
      INTEGER  JS,KGRID,KP,KCHECK
       
      REAL*8   ZZLAMB,H1,H2,H3,H4,H5,H6,H7,EPC,EPC2
      REAL*8,  DIMENSION(:,:),ALLOCATABLE::ATMP
      REAL*8,  DIMENSION(:),  ALLOCATABLE::ZEPK2

      KCHECK  = 0
      ZEPK    = ZEPKO

      ALLOCATE( ATMP(NEPK2,4),ZEPK2(NEPK2) )

      EPC  = EPSLONCA(JS,KGRID,KP)
      EPC2 = EPC**1.5
  
      H1  = CPSI(JS,KGRID,KP)
      H2  = DCDPSIL(JS,KGRID,KP)
      H3  = DCDPSIL2(JS,KGRID,KP)

C     H4  = D EPSLONC**1.5/ D PSI = OMEGASCA/H3
      H4  = B0K/OMEGACI0*ESPECIES_TEM(JS,KGRID,KP)
     &      *ESPECIES_Z(1)/ESPECIES_Z(KP)
      H4  = OMEGASCA(JS,KGRID,KP)/H4
      H5  = DEPSCDPSI2(JS,KGRID,KP)

      H6  = DEPSALPHADPSI(JS,KGRID,KP)
      H7  = DEPSALPHADPSI2(JS,KGRID,KP)

      ZEPK2 = ZEPK**1.5
      FEPK  = 1./(ZEPK2+EPC2)
       
C     F0: KOPT=0
      ATMP(:,1)       = H1*ZZLAMB*FEPK

      ZRESU(:,1,0,KP) = ATMP(:,1)
      ZRESU(:,2,0,KP) = ATMP(:,1)
      ZRESU(:,3,0,KP) = ATMP(:,1)
          
C     D F0/D PSI: KOPT=1       
      ATMP(:,2) = 1.5*ZEPK2*H6-H4
      ATMP(:,3) = H2-2.5*H6+FEPK*ATMP(:,2)
      ATMP(:,4) = ATMP(:,1)*ATMP(:,3)
      
      ZRESU(:,1,1,KP) = ATMP(:,4)
      ZRESU(:,2,1,KP) = ATMP(:,4)
      ZRESU(:,3,1,KP) = ATMP(:,4)
       
C     D F0/D (NORMALISED) EPSLONK: KOPT=2
      ATMP(:,4) =-1.5*ATMP(:,1)*FEPK*SQRT(ZEPK)

      ZRESU(:,1,2,KP) = ATMP(:,4)
      ZRESU(:,2,2,KP) = ATMP(:,4)
      ZRESU(:,3,2,KP) = ATMP(:,4)

C     D F0/D LAMBDA: KOPT=3
      ZRESU(:,:,3,KP) = ZRESU(:,:,0,KP)/ZZLAMB
      
      IF (ABS(PSPECIES_FOWP(KP))+ABS(PSPECIES_FOWT(KP)).GT.0.) THEN

C     D^2 F0/D PSI^2: KOPT=4       
      ATMP(:,4) = ATMP(:,1)*(ATMP(:,3)**2+H3-2.5*H7+
     &            (FEPK*ATMP(:,2))**2+FEPK*(1.5*ZEPK2*
     &            (H7-1.5*H6**2)-H5))
 
      ZRESU(:,1,4,KP) = ATMP(:,4)
      ZRESU(:,2,4,KP) = ATMP(:,4)
      ZRESU(:,3,4,KP) = ATMP(:,4)

C     D^2 F0/D EPSILONK/D PSI: KOPT=5
      ATMP(:,4) =-1.5*ATMP(:,1)*FEPK*SQRT(ZEPK)*(H2-3.*H6+
     &            2.*FEPK*ATMP(:,2))
      ZRESU(:,1,5,KP) = ATMP(:,4)
      ZRESU(:,2,5,KP) = ATMP(:,4)
      ZRESU(:,3,5,KP) = ATMP(:,4)

C     D^2 F0/D LAMBDA/D PSI: KOPT=6
      ZRESU(:,:,6,KP) = ZRESU(:,:,1,KP)/ZZLAMB

      ENDIF

      DEALLOCATE(ATMP,ZEPK2) 

      RETURN
      END SUBROUTINE KDISTRIBF_TYPE7
                   
C==========================================================
C ALLOCATE VARIABLES FOR ANISOTROPIC SUBROUTINES        
C NOTE: INCLUDE TWO END POINTS FOR ZOMEGABP,ZOMEGABT AND
C       ZOMEGADT
C==========================================================
      SUBROUTINE ZALLOCANISO(KALLOC)
      USE GLOBALM 
      USE KINETICM
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      INTEGER  KCHECK,MMLMAX,KALLOC

      NEPK2   = 2*NEPK-1

      IF (KALLOC.EQ.1) THEN

      ALLOCATE(
     & ZZI1(NRP1,2,3:NSPECIES),    ZZI3(NRP1,2,3:NSPECIES),
     & CPSI(NRP1,2,1:NSPECIES),    
     & DCDPSIL(NRP1,2,NSPECIES), DEPSALPHADPSI(NRP1,2,NSPECIES),
     & DCDPSIL2(NRP1,2,NSPECIES),DEPSALPHADPSI2(NRP1,2,NSPECIES),
     & DEPSCDPSI2(NRP1,2,NSPECIES) )

      ALLOCATE(
     & ZOMEGABP (NRP1,2*NLAMK+2,2),ZOMEGABT (NRP1,2*NLAMK+2,2),
     & ZOMEGADT (NRP1,2*NLAMK+2,2),TPSI0    (NRP1,2*NLAMK,2),
     & TPSI0DPSI(NRP1,2*NLAMK,2),  TPSI0DLAM(NRP1,2*NLAMK,2),
     & HPSI0    (NRP1,2*NLAMK,2),  HPSI0DPSI(NRP1,2*NLAMK,2),
     & HPSI0DLAM(NRP1,2*NLAMK,2),  HATJ0    (NRP1,2*NLAMK,2) )

      ALLOCATE(
     & ZDHMAXDPSI(NRP1,2),         ZDHMINMAXDPSI(NRP1,2),
     & ZDHMAXDPSI2(NRP1,2),        ZDHMINMAXDPSI2(NRP1,2),
     & HMIN3DPSI(NRP1,2),          HMIN3DPSI2(NRP1,2),
     & STYPE4DPSI(NRP1,2),         STYPE4DPSI2(NRP1,2), 
     & RTYPE4DPSI(NRP1,2),         RTYPE4DPSI2(NRP1,2) )

      ALLOCATE( EPK(NEPK),ZEPKO(NEPK2),ZEPKN(NEPK2) )

      ZZI1       = 0.
      ZZI3       = 0.
       
      CPSI       = 0.
      DCDPSIL    = 0.
      DCDPSIL2   = 0.
     
      ZOMEGABP   = 0.
      ZOMEGABT   = 0.
      ZOMEGADT   = 0.
      
      ZDHMAXDPSI = 0.
      ZDHMINMAXDPSI =0.

      ENDIF

      ALLOCATE(
     &   ZRQ1(NCHI+1),ZRQ2(NCHI+1),ZRQ3(NCHI+1),
     &   ZRQ1K(2*NCHI0+2),ZRQ2K(2*NCHI0+2),ZRQ3K(2*NCHI0+2),
     &   FESS1(NEPK2,0:4,0:3), FESS2(NEPK2,0:4,0:3), 
     &   FESS3(NEPK2,0:4,0:3) )

      ALLOCATE(ZRESU(NEPK2,3,0:6,NSPECIES), ZEPK(NEPK2))
      
      ALLOCATE(
     &   ZDTK(NEPK2),ZDTK2(NEPK2), 
     &   FEPK(NEPK2),FEPKC(NEPK2))
      
      MMLMAX = M2-M1-(M1-M2) + 1      
      ALLOCATE( 
     &   ZGL0PA(MMLMAX),ZGL1PA(MMLMAX),ZGL2PA(MMLMAX),ZGL3PA(MMLMAX),
     &   ZGL0PE(MMLMAX),ZGL1PE(MMLMAX),ZGL2PE(MMLMAX),ZGL3PE(MMLMAX),
     &   ZGL0DP(MMLMAX),ZGL1DP(MMLMAX),ZGL2DP(MMLMAX),ZGL3DP(MMLMAX))

      IF (IFOWP.EQ.1.OR.IFOWT.EQ.1) ALLOCATE( 
     &  ZGL0PA1(MMLMAX),ZGL1PA1(MMLMAX),ZGL2PA1(MMLMAX),ZGL3PA1(MMLMAX),
     &  ZGL0PE1(MMLMAX),ZGL1PE1(MMLMAX),ZGL2PE1(MMLMAX),ZGL3PE1(MMLMAX),
     &  ZGL0DP1(MMLMAX),ZGL1DP1(MMLMAX),ZGL2DP1(MMLMAX),ZGL3DP1(MMLMAX))

      KCHECK = 0 

      FESS1      = 0.
      FESS2      = 0.
      FESS3      = 0.

      ZRQ1       = 0.
      ZRQ2       = 0.
      ZRQ3       = 0.
    
      ZRQ1K      = 0.
      ZRQ2K      = 0.
      ZRQ3K      = 0.
      
      ZRESU      = 0.

      ZDTK       = 0.
      ZDTK2      = 0. 
      FEPK       = 0.
      FEPKC      = 0.

      ZGL0PA     = 0.
      ZGL1PA     = 0.
      ZGL2PA     = 0.
      ZGL3PA     = 0.
      ZGL0PE     = 0.
      ZGL1PE     = 0.
      ZGL2PE     = 0.
      ZGL3PE     = 0.
      ZGL0DP     = 0.
      ZGL1DP     = 0.
      ZGL2DP     = 0.
      ZGL3DP     = 0.

      RETURN
      END SUBROUTINE ZALLOCANISO 

C==========================================================
C   DEALLOCATE VARIABLES FOR ANISOTROPIC SUBROUTINES      =
C==========================================================
      SUBROUTINE ZDEALLOCANISO(KALLOC)
      USE GLOBALM 
      USE KINETICM
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      INTEGER  KALLOC

      IF (KALLOC.EQ.1) THEN

      DEALLOCATE(
     &   ZZI1,ZZI3,CPSI,DCDPSIL,DCDPSIL2)
      DEALLOCATE(DEPSALPHADPSI,DEPSALPHADPSI2,DEPSCDPSI2)
      DEALLOCATE(ZOMEGABP,ZOMEGABT,ZOMEGADT)
      DEALLOCATE(TPSI0,TPSI0DPSI,TPSI0DLAM,HPSI0,HPSI0DPSI,HPSI0DLAM,
     &           HATJ0)
      DEALLOCATE
     & (ZDHMAXDPSI,ZDHMINMAXDPSI,ZDHMAXDPSI2,ZDHMINMAXDPSI2,
     &  HMIN3DPSI,HMIN3DPSI2,RTYPE4DPSI,RTYPE4DPSI2,
     &  STYPE4DPSI,STYPE4DPSI2)
      DEALLOCATE(EPK,ZEPKO,ZEPKN)

      ENDIF

      DEALLOCATE(
     &   ZRQ1,ZRQ2,ZRQ3,
     &   ZRQ1K,ZRQ2K,ZRQ3K,
     &   FESS1,FESS2,FESS3 )

      DEALLOCATE(ZRESU,ZEPK)
      
      DEALLOCATE(
     &   ZDTK,ZDTK2, 
     &   FEPK,FEPKC)
      
      DEALLOCATE( 
     &   ZGL0PA,ZGL1PA,ZGL2PA,ZGL3PA,
     &   ZGL0PE,ZGL1PE,ZGL2PE,ZGL3PE,
     &   ZGL0DP,ZGL1DP,ZGL2DP,ZGL3DP)

      IF (IFOWP.EQ.1.OR.IFOWT.EQ.1) DEALLOCATE( 
     &   ZGL0PA1,ZGL1PA1,ZGL2PA1,ZGL3PA1,
     &   ZGL0PE1,ZGL1PE1,ZGL2PE1,ZGL3PE1,
     &   ZGL0DP1,ZGL1DP1,ZGL2DP1,ZGL3DP1)

      RETURN
      END SUBROUTINE ZDEALLOCANISO 

C==========================================================
C COMPUTE C(PSI) FOR IF0TYPE=3,7 
C==========================================================
      SUBROUTINE KCPSI_TYPE3(KP)
      
      USE KINETICM 
      USE GLOBALM
      USE RCOMDM  
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER  KP,JS,K,KCHECK
      REAL*8   RTMP,ZDAK,KZJAM0,ZDPSIDS

      KCHECK=0

      DO K=1,2
      DO JS=2-K+1,NR 
         ZDAK = ESPECIES_DENF(JS,K,KP)
            
C        OBTAIN ALPHAA3 I1 I3      
         CALL  ZI1I3(JS,K,KP,ALPHAA1(JS,K,KP))

         RTMP = 4.0*PI*ZDAK*ESPECIES_PRE(JS,K,2)
     &          /ALPHAA1(JS,K,KP)/ZZI1(JS,K,KP)
 
         IF (K.EQ.1) KZJAM0 = REAL(JACOBI(JS,1))     
         IF (K.EQ.2) KZJAM0 = REAL(JACOBM(JS,1))     
 
         IF (K.EQ.1) ZDPSIDS = DPSIDS(JS)
         IF (K.EQ.2) ZDPSIDS = DPSIDSM(JS)
      
         CPSI(JS,K,KP) = RTMP*B0K*KZJAM0/ZDPSIDS
      ENDDO
      ENDDO
      
      CPSI(1,1,KP)    = 2.*CPSI(1,2,KP)-CPSI(2,1,KP)
      CPSI(NRP1,1,KP) = 2.*CPSI(NR,2,KP)-CPSI(NR,1,KP)

C     MAKE SURE THAT CPSII>0
      IF (CPSI(1,1,KP).LE.0.)    CPSI(1,1,KP)    = CPSI(1,2,KP)
      IF (CPSI(NRP1,1,KP).LE.0.) CPSI(NRP1,1,KP) = CPSI(NR,2,KP)

      RETURN
      END SUBROUTINE KCPSI_TYPE3

C==========================================================
C COMPUTE C(PSI) FOR IF0TYPE=4 
C==========================================================
      SUBROUTINE KCPSI_TYPE4(KP)
      
      USE KINETICM 
      USE GLOBALM
      USE RCOMDM  
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER  KP,JS,K,N,J,KCHECK
      REAL*8   RTMP,H2,H3,ZLAMB,ZF1L,ZF1LP,ZF2L,ZF2LP,ZF3L,ZF3LP,ZHLAM,
     &         KZJAM0,ZDPSIDS

      KCHECK=1

C     GO THROUGH EACH FLUX SURFACE
      DO K=1,2
      DO JS=2-K+1,NR 

C     GO THROUGH PITCH ANGLE INTEGRATION FOR PASSING PARTICLES
      ZF1L    = 0.
      ZF1LP   = 0.
      DO N   = 1,NLAMK1(JS,K)-1
      ZHLAM  = LAMK1(JS,N+1,K)-LAMK1(JS,N,K)
      DO  J = 0,1
         ZLAMB  = 0.5*((1.+ZW)*LAMK1(JS,N+J,K)
     &            +(1.-ZW)*LAMK1(JS,N+1-J,K))

         H2     = ZLAMB/HTYPE4C+RTYPE4(JS,K)*ABS(1.-ZLAMB/HTYPE4C)
         H3     = ZOMEGABP(JS,2*N-1+J+1,K)
         RTMP   = 0.5*ZHLAM/H3
         ZF1L   = ZF1L  + RTMP/H2**1.5
         ZF1LP  = ZF1LP + RTMP/H2**2.5
      ENDDO
      ENDDO
      ZF3L  = ZF1L
      ZF3LP = ZF1LP

C     GO THROUGH PITCH ANGLE INTEGRATION FOR TRAPPED PARTICLES
      ZF2L   = 0.
      ZF2LP  = 0.
      DO  N  = 1,NLAMK0(JS,K)-1 
      ZHLAM  = LAMK0(JS,N+1,K)-LAMK0(JS,N,K)
      DO  J  = 0,1
         ZLAMB = 0.5*((1.+ZW)*LAMK0(JS,N+J,K)
     &           +(1.-ZW)*LAMK0(JS,N+1-J,K))

         H2     = ZLAMB/HTYPE4C+RTYPE4(JS,K)*ABS(1.-ZLAMB/HTYPE4C)
         H3     = ZOMEGABT(JS,2*N-1+J+1,K)

         RTMP   = 0.5*ZHLAM/H3
         ZF2L   = ZF2L  + RTMP/H2**1.5
         ZF2LP  = ZF2LP + RTMP/H2**2.5
      ENDDO
      ENDDO
      
C     0.5 COMES FROM \hat_g FACTOR FOR TRAPPED PARTICLES
      RTMP  = 0.5 
      ZZZI1 = (ZF1L +2.0*RTMP*ZF2L +ZF3L)/2.
      ZZZI3 = (ZF1LP+2.0*RTMP*ZF2LP+ZF3LP)/2.

      EPSALPHA(JS,K,KP) = ZZZI1/ZZZI3*ESPECIES_PRE(JS,K,KP)/
     &                    ESPECIES_DEN(JS,K,KP)

      IF (K.EQ.1) KZJAM0 = REAL(JACOBI(JS,1))     
      IF (K.EQ.2) KZJAM0 = REAL(JACOBM(JS,1))     
 
      IF (K.EQ.1) ZDPSIDS = DPSIDS(JS)
      IF (K.EQ.2) ZDPSIDS = DPSIDSM(JS)

      CPSI(JS,K,KP) = 2.*B0K*KZJAM0/ZDPSIDS*ESPECIES_PRE(JS,K,KP)/ZZZI3

      IF (KCHECK.EQ.1.AND.K.EQ.2) THEN
         IF (JS.EQ.1) WRITE(*,*) 'CS DPSIDS/B0/JAC/2 I1 I3 RHO PEQ'
         WRITE(*,110) CSM(JS),ZDPSIDS/B0K/KZJAM0/2.,ZZZI1,ZZZI3,
     &                RHOM(JS),PEQM(JS)
      ENDIF
 110  FORMAT(6(E14.5))

      ENDDO
      ENDDO
      CPSI(1,1,KP)        = 2.*CPSI(1,2,KP)-CPSI(2,1,KP)
      CPSI(NRP1,1,KP)     = 2.*CPSI(NR,2,KP)-CPSI(NR,1,KP)
      EPSALPHA(1,1,KP)    = 2.*EPSALPHA(1,2,KP)-EPSALPHA(2,1,KP)
      EPSALPHA(NRP1,1,KP) = 2.*EPSALPHA(NR,2,KP)-EPSALPHA(NR,1,KP)

C     MAKE SURE THAT CPSII>0, EPSALPHA>0
      IF (CPSI(1,1,KP).LE.0.)    CPSI(1,1,KP)     = CPSI(1,2,KP)
      IF (CPSI(NRP1,1,KP).LE.0.) CPSI(NRP1,1,KP)  = CPSI(NR,2,KP)
      IF (EPSALPHA(1,1,KP).LE.0.) EPSALPHA(1,1,KP)= EPSALPHA(1,2,KP)
      IF (EPSALPHA(NRP1,1,KP).LE.0.) EPSALPHA(NRP1,1,KP) = 
     &                               EPSALPHA(NR,2,KP)

      RETURN
      END SUBROUTINE KCPSI_TYPE4

C==========================================================
C COMPUTE C(PSI) FOR IF0TYPE=5 (RE MODEL-1)
C==========================================================
      SUBROUTINE KCPSI_TYPE5(KP)
      
      USE KINETICM 
      USE GLOBALM
      USE RCOMDM  
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER  KP,JS,K,N,J,KCHECK
      REAL*8   RTMP,H2,H3,ZLAMB,ZF1L,ZF2L,ZF3L,ZHLAM,
     &         KZJAM0,ZDPSIDS,ZZZI2
      REAL*8,DIMENSION(:),ALLOCATABLE::SCREL

      KCHECK=1

C     SPLINE COEFFICIENTS FOR PITCH-ANGLE DISTRIBUTION FUNCTION
      ALLOCATE(SCREL(NLAMK))
      CALL SPLINE(ESPECIES_REL(:,1),ESPECIES_REL(:,2),
     &            NLAMK,RTMP,RTMP,SCREL)

C     GO THROUGH EACH FLUX SURFACE
      DO K=1,2
      DO JS=2-K+1,NR 

C     GO THROUGH PITCH ANGLE INTEGRATION FOR PASSING PARTICLES
      ZF1L    = 0.
      DO N   = 1,NLAMK1(JS,K)-1
      ZHLAM  = LAMK1(JS,N+1,K)-LAMK1(JS,N,K)
      DO  J = 0,1
         ZLAMB  = 0.5*((1.+ZW)*LAMK1(JS,N+J,K)
     &            +(1.-ZW)*LAMK1(JS,N+1-J,K))

         RTMP = ZLAMB/HKMAX(JS,K)
         CALL SPLINT(ESPECIES_REL(:,1),ESPECIES_REL(:,2),
     &               SCREL,NLAMK,RTMP,H2)

         H3     = ZOMEGABP(JS,2*N-1+J+1,K)
         RTMP   = 0.5*ZHLAM*H2/H3
         ZF1L   = ZF1L  + RTMP
      ENDDO
      ENDDO
      ZF3L  = ZF1L

C     GO THROUGH PITCH ANGLE INTEGRATION FOR TRAPPED PARTICLES
      ZF2L   = 0.
      DO  N  = 1,NLAMK0(JS,K)-1 
      ZHLAM  = LAMK0(JS,N+1,K)-LAMK0(JS,N,K)
      DO  J  = 0,1
         ZLAMB = 0.5*((1.+ZW)*LAMK0(JS,N+J,K)
     &           +(1.-ZW)*LAMK0(JS,N+1-J,K))

         RTMP = ZLAMB/HKMAX(JS,K)
         CALL SPLINT(ESPECIES_REL(:,1),ESPECIES_REL(:,2),
     &               SCREL,NLAMK,RTMP,H2)

         H3     = ZOMEGABT(JS,2*N-1+J+1,K)
         RTMP   = 0.5*ZHLAM*H2/H3
         ZF2L   = ZF2L  + RTMP
      ENDDO
      ENDDO
      
C     0.5 COMES FROM \hat_g FACTOR FOR TRAPPED PARTICLES
      RTMP  = 0.5 
      ZZZI2 = (ZF1L +2.0*RTMP*ZF2L +ZF3L)/2.

      CALL ZITYPE5(JS,K,KP,ALPHAA1(JS,K,KP))
      ZZZI1 = ZZZI1*EPSALPHA(JS,K,KP)**1.5

      IF (K.EQ.1) KZJAM0 = REAL(JACOBI(JS,1))     
      IF (K.EQ.2) KZJAM0 = REAL(JACOBM(JS,1))     
 
      IF (K.EQ.1) ZDPSIDS = DPSIDS(JS)
      IF (K.EQ.2) ZDPSIDS = DPSIDSM(JS)

      CPSI(JS,K,KP) = 2.*SQRT(PI)*B0K*KZJAM0/ZDPSIDS*
     &                ESPECIES_DEN(JS,K,KP)/ZZZI1/ZZZI2

      IF (KCHECK.EQ.1.AND.K.EQ.2) THEN
         IF (JS.EQ.1) WRITE(*,*) 'CS DPSIDS/B0/JAC/2 I1 I3 RHO PEQ'
         WRITE(*,110) CSM(JS),ZDPSIDS/B0K/KZJAM0/2.,ZZZI1,ZZZI3,
     &                RHOM(JS),PEQM(JS)
      ENDIF
 110  FORMAT(6(E14.5))

      ENDDO
      ENDDO
      CPSI(1,1,KP)        = 2.*CPSI(1,2,KP)-CPSI(2,1,KP)
      CPSI(NRP1,1,KP)     = 2.*CPSI(NR,2,KP)-CPSI(NR,1,KP)
      EPSALPHA(1,1,KP)    = 2.*EPSALPHA(1,2,KP)-EPSALPHA(2,1,KP)
      EPSALPHA(NRP1,1,KP) = 2.*EPSALPHA(NR,2,KP)-EPSALPHA(NR,1,KP)

C     MAKE SURE THAT CPSII>0, EPSALPHA>0
      IF (CPSI(1,1,KP).LE.0.)    CPSI(1,1,KP)     = CPSI(1,2,KP)
      IF (CPSI(NRP1,1,KP).LE.0.) CPSI(NRP1,1,KP)  = CPSI(NR,2,KP)
      IF (EPSALPHA(1,1,KP).LE.0.) EPSALPHA(1,1,KP)= EPSALPHA(1,2,KP)
      IF (EPSALPHA(NRP1,1,KP).LE.0.) EPSALPHA(NRP1,1,KP) = 
     &                               EPSALPHA(NR,2,KP)

      DEALLOCATE(SCREL)

      RETURN
      END SUBROUTINE KCPSI_TYPE5

C==========================================================
C COMPUTE C(PSI) FOR IF0TYPE=6 (RE MODEL-2)
C==========================================================
      SUBROUTINE KCPSI_TYPE6(KP)
      
      USE KINETICM 
      USE GLOBALM
      USE RCOMDM  
      USE DIMENSIM
      USE ANISOTROPICM
      IMPLICIT NONE
      
      INTEGER  KP,JS,K,KCHECK
      REAL*8   RTMP,ZDAK,KZJAM0

      KCHECK=1

      RTMP = OMEGACI0/B0EXP*SQRT(ESPECIES_M(1)/ESPECIES_M(KP)) !=e/sqrt(M)

      DO K=1,2
      DO JS=2-K+1,NR 
         ZDAK   = ABS(EQJPARA(JS,K))*JRE_EQFRAC
         KZJAM0 = EQ1_H(JS,K)
            
C        OBTAIN I1      
         CALL  ZI1I3(JS,K,KP,ALPHAA1(JS,K,KP))

         CPSI(JS,K,KP) = SQRT(2.0*PI)*ZDAK/KZJAM0/RTMP/ZZI1(JS,K,KP)
      ENDDO
      ENDDO
      
      CPSI(1,1,KP)    = 2.*CPSI(1,2,KP)-CPSI(2,1,KP)
      CPSI(NRP1,1,KP) = 2.*CPSI(NR,2,KP)-CPSI(NR,1,KP)

C     MAKE SURE THAT CPSII>0
      IF (CPSI(1,1,KP).LE.0.)    CPSI(1,1,KP)    = CPSI(1,2,KP)
      IF (CPSI(NRP1,1,KP).LE.0.) CPSI(NRP1,1,KP) = CPSI(NR,2,KP)

      IF (KCHECK.EQ.1) THEN
         WRITE(*,*) 'KCPSI_TYPE6: EQJPARA:'
         DO JS=1,NRP1
            WRITE(*,*) EQJPARA(JS,1)
         ENDDO
      ENDIF     

      RETURN
      END SUBROUTINE KCPSI_TYPE6

C=============================================================
C I-FACTOR FOR ADIABATIC CONTRIBUTIONS, FOR GENERIC ANISOTROPIC
C DISTRIBUTION
C KOPT = 0: VANISHING ORBIT WIDTH CONTRIBUTION
C        1: FOW CORRECTION 1
C        2: FOW CORRECTION 2
C OUTOUT:
C   ZVIA(1,:)=PRESSURE EQ: I_XI(LAMBDA)    
C   ZVIA(2,:)=PRESSURE EQ: I_Q(LAMBDA)
C   ZVIA(3,:)=PRESSURE EQ: I_DPHI(LAMBDA)
C   ZVIA(4,:)=DPHI EQ: I_XI(LAMBDA)    
C   ZVIA(5,:)=DPHI EQ: I_Q(LAMBDA)
C   ZVIA(5,:)=DPHI EQ: I_DPHI(LAMBDA)
C RFAC=LORENTZ FACTOR GAMMA FOR RE
C=============================================================
      SUBROUTINE KIA_ADIABATIC(JS,KGRID,KOPT,ZZLAMB,ZVIA)
      
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      USE DIMENSIM
      USE RCOMDM
      IMPLICIT NONE

      INTEGER   JS,KGRID,K,KCHECK,KOPT,KP,KEXEC
      REAL*8    ZZLAMB,ZDPSIDS,HEPK,ZHK,ZHMIN,ZPP1JX,ZPP2JX,ZPP3JX,
     &          ZPP1JQ,ZPP2JQ,ZPP3JQ,ZPP1JP,ZPP2JP,ZPP3JP,
     &          RTMP1,RTMP2,RTMP3,FAC0,ZZT,
     &          H1,H2,H3,H4,H5,H6,H7,H8,H9,
     &          PSI0L,PSI0DPSIL,PSI0DLAML
      REAL*8,DIMENSION(:,:), ALLOCATABLE:: ZT1,ZT3,ZT4
      REAL*8,DIMENSION(:),   ALLOCATABLE:: EFAC,RFAC
      REAL*8    CTMP1,CTMP2,CTMP3,ZVIA(6,NSPECIES0)
     
      KCHECK = 0

      ZVIA   = 0.

      IF (KOPT.EQ.1.AND.IFOWPSI0.EQ.2) RETURN

      ALLOCATE( ZT1(NEPK2,3),ZT3(NEPK2,3),ZT4(NEPK2,3),
     &          EFAC(NEPK2),RFAC(NEPK2) )

      ZHK   = HKMAX(JS,KGRID)
      ZHMIN = HKMIN(JS,KGRID)

      IF(KGRID.EQ.1) ZDPSIDS = DPSIDS(JS)
      IF(KGRID.EQ.2) ZDPSIDS = DPSIDSM(JS) 
       
      IF(KGRID.EQ.1) ZZT     = T(JS)
      IF(KGRID.EQ.2) ZZT     = TM(JS)

      DO KP=1,NSPECIES

C     EFAC=1 FOR PARTICLES WITH FINITE BIRTH ENERGY, AND
C     EFAC=PI/2/COS^2(ZEPK*PI/2) FOR PARTICLES WITH INFINITE
C     BIRTH ENERGY 
      EFAC = 1.
      RFAC = 1.
      ZEPK = ZEPKO
      IF (ISPECIES_EK(KP).EQ.0) THEN
         EFAC = 0.5*PI/(COS(0.5*PI*ZEPKO))**2
         ZEPK = ZEPKN
      ENDIF

      IF (ISPECIES_F0(KP).EQ.5) THEN
      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)                            !=P0
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI  !=N0
      H6 = 8.1872E-14                                      !=me*c^2
      H5 = H1/H2/H6                                        !=T0/(me*c^2)
      H7 = ESPECIES_TEM(JS,KGRID,2)/ALPHAA1(JS,KGRID,KP)*H5!=EPS_A/(me*c^2)
C     RFAC = ZEPK*H7+1.
      RFAC = 1.
      ENDIF

      H8 = ESPECIES_Z(KP)*OMEGACI0/B0K/EPSALPHA(JS,KGRID,KP)
      H9 = ESPECIES_Z(KP)/EPSALPHA(JS,KGRID,KP)

      
C     PASSING PARTICLE CONTRIBUTION
C     FOW CORRECTIONS ONLY FOR TRUELY ANISOTROPIC DISTRIBUTIONS
      IF (ZZLAMB.GE.0.0.AND.ZZLAMB.LE.ZHMIN.AND.
     &    ABS(PSPECIES_AP(KP)).GT.0.) THEN

      IF (KOPT.EQ.1) THEN
      IF (KP.EQ.1.OR.KP.EQ.2) THEN
         H7 = ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &        (ESPECIES_M(KP)/ESPECIES_M(1))**2
         PSI0L     = H7*TPSI0L
         PSI0DPSIL = H7*TPSI0DPSIL
         PSI0DLAML = H7*TPSI0DLAML
      ELSE
         H7 = ESPECIES_Z(1)/ESPECIES_Z(KP)*ESPECIES_M(KP)/ESPECIES_M(1)
         PSI0L     = H7*HPSI0L
         PSI0DPSIL = H7*HPSI0DPSIL
         PSI0DLAML = H7*HPSI0DLAML
      ENDIF
      ENDIF

      KEXEC = 0

      IF (KOPT.EQ.0.AND.ABS(PSPECIES_AP(KP)).GT.0.) THEN
         IF (.NOT.((ISPECIES_F0(KP).EQ.0.OR.ISPECIES_F0(KP).EQ.1
     &       .OR.ISPECIES_F0(KP).EQ.2).AND.KFASTRUN.EQ.1)) THEN
         ZT1(:,2) = ZEPK**1.5*EFAC/RFAC**(1.5)
         ZT1(:,1) = ZRESU(:,1,1,KP)*ZT1(:,2)
         ZT1(:,3) = ZRESU(:,3,1,KP)*ZT1(:,2)
         ZT3(:,1) = ZRESU(:,1,3,KP)*ZT1(:,2)
         ZT3(:,3) = ZRESU(:,3,3,KP)*ZT1(:,2)
         ZT4(:,1) = ZRESU(:,1,2,KP)*ZT1(:,2)
         ZT4(:,3) = ZRESU(:,3,2,KP)*ZT1(:,2)
         KEXEC = 1
         ENDIF
      ELSEIF (KOPT.EQ.1.AND.ABS(PSPECIES_FOWP(KP)).GT.0..AND.
     &        (ISPECIES_F0(KP).EQ.3.OR.ISPECIES_F0(KP).EQ.4.OR.
     &         ISPECIES_F0(KP).EQ.5).OR.ISPECIES_F0(KP).EQ.6) THEN
         FAC0 = PSPECIES_FOWP(KP)*SQRT(2.*ESPECIES_M(1)/
     &          ESPECIES_M(KP)*EPSALPHA(JS,KGRID,KP))

         ZT1(:,2) = FAC0*PSI0DPSIL*ZEPK**2*EFAC/RFAC**2
         ZT1(:,1) = ZRESU(:,1,1,KP)*ZT1(:,2)
         ZT1(:,3) =-ZRESU(:,3,1,KP)*ZT1(:,2)

         ZT1(:,2) = FAC0*PSI0DLAML*ZEPK**2*EFAC/RFAC**2
         ZT3(:,1) = ZRESU(:,1,1,KP)*ZT1(:,2)
         ZT3(:,3) =-ZRESU(:,3,1,KP)*ZT1(:,2)

         ZT1(:,2) = FAC0*PSI0L*ZEPK*EFAC*0.5/RFAC**2
         ZT4(:,1) = ZRESU(:,1,1,KP)*ZT1(:,2)
         ZT4(:,3) =-ZRESU(:,3,1,KP)*ZT1(:,2)

         ZT1(:,2) = FAC0*PSI0L*ZEPK**2*EFAC/RFAC**2
         ZT1(:,1) = ZT1(:,1) + ZRESU(:,1,4,KP)*ZT1(:,2)
         ZT1(:,3) = ZT1(:,3) - ZRESU(:,3,4,KP)*ZT1(:,2)
         ZT3(:,1) = ZT3(:,1) + ZRESU(:,1,6,KP)*ZT1(:,2)
         ZT3(:,3) = ZT3(:,3) - ZRESU(:,3,6,KP)*ZT1(:,2)
         ZT4(:,1) = ZT4(:,1) + ZRESU(:,1,5,KP)*ZT1(:,2)
         ZT4(:,3) = ZT4(:,3) - ZRESU(:,3,5,KP)*ZT1(:,2)
         KEXEC = 1
      ELSEIF (KOPT.EQ.2.AND.ABS(PSPECIES_FOWP(KP)).GT.0..AND.
     &        (ISPECIES_F0(KP).EQ.3.OR.ISPECIES_F0(KP).EQ.4.OR.
     &         ISPECIES_F0(KP).EQ.5).OR.ISPECIES_F0(KP).EQ.6) THEN
         FAC0 =-PSPECIES_FOWP(KP)*SQRT(2.*ESPECIES_M(1)/
     &          ESPECIES_M(KP)*EPSALPHA(JS,KGRID,KP))*
     &          ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &          ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0

         ZT1(:,2) = FAC0*ZEPK**2*EFAC/RFAC**2
         ZT1(:,1) = ZRESU(:,1,4,KP)*ZT1(:,2)
         ZT1(:,3) =-ZRESU(:,3,4,KP)*ZT1(:,2)
         ZT3(:,1) = ZRESU(:,1,6,KP)*ZT1(:,2)
         ZT3(:,3) =-ZRESU(:,3,6,KP)*ZT1(:,2)
         ZT4(:,1) = ZRESU(:,1,5,KP)*ZT1(:,2)
         ZT4(:,3) =-ZRESU(:,3,5,KP)*ZT1(:,2)
         KEXEC = 1
      ENDIF

C     FOR PRESSURE EQ
      IF (KEXEC.EQ.1) THEN
      ZPP1JX  = 0.
      ZPP3JX  = 0.
      ZPP1JQ  = 0.
      ZPP3JQ  = 0.
      ZPP1JP  = 0.
      ZPP3JP  = 0.

      DO K=1,NEPK-1
      HEPK    = EPK(K+1)-EPK(K)
      RTMP1   = 0.5*HEPK*(ZT1(2*K-1,1)+ZT1(2*K,1))
      RTMP2   = 0.5*HEPK*(ZT1(2*K-1,3)+ZT1(2*K,3))
      ZPP1JX  =  ZPP1JX+RTMP1
      ZPP3JX  =  ZPP3JX+RTMP2 
             
      RTMP1   = 0.5*HEPK*(ZT3(2*K-1,1)+ZT3(2*K,1))
      RTMP2   = 0.5*HEPK*(ZT3(2*K-1,3)+ZT3(2*K,3))
      ZPP1JQ  =  ZPP1JQ+RTMP1
      ZPP3JQ  =  ZPP3JQ+RTMP2                         

      RTMP1   = 0.5*HEPK*(ZT4(2*K-1,1)+ZT4(2*K,1))
      RTMP2   = 0.5*HEPK*(ZT4(2*K-1,3)+ZT4(2*K,3))
      ZPP1JP  =  ZPP1JP+RTMP1
      ZPP3JP  =  ZPP3JP+RTMP2 
      ENDDO     
      
      CTMP1   =-(ZPP1JX+ZPP3JX)*PSPECIES_AP(KP)*ZDPSIDS
      CTMP2   =-(ZPP1JQ+ZPP3JQ)*PSPECIES_AP(KP)*ZZLAMB
      CTMP3   = (ZPP1JP+ZPP3JP)*PSPECIES_AP(KP)*H8
      ZVIA(1,KP) = ZVIA(1,KP) + CTMP1 
      ZVIA(2,KP) = ZVIA(2,KP) + CTMP2
      ZVIA(3,KP) = ZVIA(3,KP) + CTMP3 
      
C     FOR DPHI EQ
      ZPP1JX  = 0.
      ZPP3JX  = 0.
      ZPP1JQ  = 0.
      ZPP3JQ  = 0.
      ZPP1JP  = 0.
      ZPP3JP  = 0.

      ZT1(:,1) = ZT1(:,1)/ZEPK
      ZT1(:,3) = ZT1(:,3)/ZEPK
      ZT3(:,1) = ZT3(:,1)/ZEPK
      ZT3(:,3) = ZT3(:,3)/ZEPK
      ZT4(:,1) = ZT4(:,1)/ZEPK
      ZT4(:,3) = ZT4(:,3)/ZEPK
      DO K=1,NEPK-1
      HEPK    = EPK(K+1)-EPK(K)
      RTMP1   = 0.5*HEPK*(ZT1(2*K-1,1)+ZT1(2*K,1))
      RTMP2   = 0.5*HEPK*(ZT1(2*K-1,3)+ZT1(2*K,3))
      ZPP1JX  =  ZPP1JX+RTMP1
      ZPP3JX  =  ZPP3JX+RTMP2 
             
      RTMP1   = 0.5*HEPK*(ZT3(2*K-1,1)+ZT3(2*K,1))
      RTMP2   = 0.5*HEPK*(ZT3(2*K-1,3)+ZT3(2*K,3))
      ZPP1JQ  =  ZPP1JQ+RTMP1
      ZPP3JQ  =  ZPP3JQ+RTMP2                         

      RTMP1   = 0.5*HEPK*(ZT4(2*K-1,1)+ZT4(2*K,1))
      RTMP2   = 0.5*HEPK*(ZT4(2*K-1,3)+ZT4(2*K,3))
      ZPP1JP  =  ZPP1JP+RTMP1
      ZPP3JP  =  ZPP3JP+RTMP2 
      ENDDO     
      
      CTMP1   =-(ZPP1JX+ZPP3JX)*PSPECIES_AP(KP)*ZDPSIDS*H9
      CTMP2   =-(ZPP1JQ+ZPP3JQ)*PSPECIES_AP(KP)*ZZLAMB*H9
      CTMP3   = (ZPP1JP+ZPP3JP)*PSPECIES_AP(KP)*H8*H9
      ZVIA(4,KP) = ZVIA(4,KP) + CTMP1 
      ZVIA(5,KP) = ZVIA(5,KP) + CTMP2
      ZVIA(6,KP) = ZVIA(6,KP) + CTMP3 
      
      IF (JS.EQ.JS0.AND.KGRID.EQ.1.AND.KCHECK.EQ.1)
     &   WRITE(*,110) KP,KOPT,ZZLAMB,CTMP1,CTMP2

      ENDIF

C     TRAPPED PARTICLE CONTRIBUTION
C     FOW CORRECTION (KOPT=1,2) GETS EXACT CANCELLATION  
      ELSEIF (ZZLAMB.GT.ZHMIN.AND.ZZLAMB.LE.ZHK.AND.
     &        ABS(PSPECIES_AT(KP)).GT.0..AND.KOPT.EQ.0.AND.
     &        (.NOT.((ISPECIES_F0(KP).EQ.0.OR.ISPECIES_F0(KP).EQ.1
     &       .OR.ISPECIES_F0(KP).EQ.2).AND.KFASTRUN.EQ.1))) THEN

      ZT1(:,1) = ZEPK**1.5*EFAC/RFAC**(1.5)
      ZT1(:,2) = ZRESU(:,2,1,KP)*ZT1(:,1)
      ZT3(:,2) = ZRESU(:,2,3,KP)*ZT1(:,1)
      ZT4(:,2) = ZRESU(:,2,2,KP)*ZT1(:,1)

C     FOR PRESSURE EQ
      ZPP2JX  = 0.
      ZPP2JQ  = 0.
      ZPP2JP  = 0.

      DO K=1,NEPK-1
         HEPK   = EPK(K+1)-EPK(K)
         RTMP1  = 0.5*HEPK*(ZT1(2*K-1,2)+ZT1(2*K,2))
         ZPP2JX =  ZPP2JX+RTMP1
          
         RTMP1  = 0.5*HEPK*(ZT3(2*K-1,2)+ZT3(2*K,2))
         ZPP2JQ = ZPP2JQ+RTMP1

         RTMP1  = 0.5*HEPK*(ZT4(2*K-1,2)+ZT4(2*K,2))
         ZPP2JP =  ZPP2JP+RTMP1
      ENDDO

      CTMP1   =-(ZPP2JX+ZPP2JX)*PSPECIES_AT(KP)*ZDPSIDS
      CTMP2   =-(ZPP2JQ+ZPP2JQ)*PSPECIES_AT(KP)*ZZLAMB
      CTMP3   = (ZPP2JP+ZPP2JP)*PSPECIES_AT(KP)*H8
      ZVIA(1,KP) = ZVIA(1,KP) + CTMP1
      ZVIA(2,KP) = ZVIA(2,KP) + CTMP2
      ZVIA(3,KP) = ZVIA(3,KP) + CTMP3

C     FOR DPHI EQ
      ZPP2JX  = 0.
      ZPP2JQ  = 0.
      ZPP2JP  = 0.

      ZT1(:,2) = ZT1(:,2)/ZEPK
      ZT3(:,2) = ZT3(:,2)/ZEPK
      ZT4(:,2) = ZT4(:,2)/ZEPK
      DO K=1,NEPK-1
         HEPK   = EPK(K+1)-EPK(K)
         RTMP1  = 0.5*HEPK*(ZT1(2*K-1,2)+ZT1(2*K,2))
         ZPP2JX =  ZPP2JX+RTMP1
          
         RTMP1  = 0.5*HEPK*(ZT3(2*K-1,2)+ZT3(2*K,2))
         ZPP2JQ = ZPP2JQ+RTMP1

         RTMP1  = 0.5*HEPK*(ZT4(2*K-1,2)+ZT4(2*K,2))
         ZPP2JP =  ZPP2JP+RTMP1
      ENDDO

      CTMP1   =-(ZPP2JX+ZPP2JX)*PSPECIES_AT(KP)*ZDPSIDS*H9
      CTMP2   =-(ZPP2JQ+ZPP2JQ)*PSPECIES_AT(KP)*ZZLAMB*H9
      CTMP3   = (ZPP2JP+ZPP2JP)*PSPECIES_AT(KP)*H8*H9
      ZVIA(4,KP) = ZVIA(4,KP) + CTMP1
      ZVIA(5,KP) = ZVIA(5,KP) + CTMP2
      ZVIA(6,KP) = ZVIA(6,KP) + CTMP3

      IF (JS.EQ.JS0.AND.KGRID.EQ.1.AND.KCHECK.EQ.1)
     &   WRITE(*,110) KP,KOPT,ZZLAMB,CTMP1,CTMP2

      ENDIF

 110  FORMAT('KIA_ADIABATIC:KP,KOPT,LAM,ZVIA1,ZVIA2:',I2,1X,I1,1X,
     &       3(E13.5,1X)) 

      ENDDO

      DEALLOCATE(ZT1,ZT3,ZT4,EFAC,RFAC)

      RETURN
      END SUBROUTINE KIA_ADIABATIC
 
C======================================================================
C   KG FOR ADIABATIC CONTRIBUTIONS
C   G_XI G_Q1 G_Q2 G_Q3 FOR PRESSURE EQ 
C   G_XI G_Q1 G_Q2 G_Q3 FOR DPHI EQ        
C======================================================================
      SUBROUTINE KG_ADIABATIC(JS,KGRID,KPARTICLE)
      
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      IMPLICIT NONE

      INTEGER JS,J,L,I,MMLMAX,KGRID,KPARTICLE,KCHECK
      
      REAL*8      H1,H2,H3,H4,ZDPSIDS,DIFPI
      COMPLEX*16  CTMP0,CTMP1,CTMP2,CTMP3
      COMPLEX*16  Z0FL,Z0FU,Z1FL,Z1FU,Z2FL,Z2FU,Z3FL,Z3FU  
      INTEGER,DIMENSION(:),ALLOCATABLE::ZRL
      REAL*8,DIMENSION(:),ALLOCATABLE::ZRG0PA,ZRG1PA,ZRG2PA,ZRG3PA,
     &                                 ZRG0PE,ZRG1PE,ZRG2PE,ZRG3PE,
     &                                 ZRG0DP,ZRG1DP,ZRG2DP,ZRG3DP
      COMPLEX*16  EPHASE
      KCHECK = 0 

      MMLMAX = M2-M1-(M1-M2) + 1      

      ALLOCATE(
     &         ZRG0PA(NCHI2+2),ZRG1PA(NCHI2+2),
     &         ZRG2PA(NCHI2+2),ZRG3PA(NCHI2+2),
     &         ZRG0PE(NCHI2+2),ZRG1PE(NCHI2+2),
     &         ZRG2PE(NCHI2+2),ZRG3PE(NCHI2+2),
     &         ZRG0DP(NCHI2+2),ZRG1DP(NCHI2+2),
     &         ZRG2DP(NCHI2+2),ZRG3DP(NCHI2+2))
      ALLOCATE(ZRL(MMLMAX))

      ZGL0PA=0.
      ZGL1PA=0.
      ZGL2PA=0.
      ZGL3PA=0.
                  
      ZGL0PE=0.
      ZGL1PE=0.
      ZGL2PE=0.
      ZGL3PE=0.

      ZGL0DP=0.
      ZGL1DP=0.
      ZGL2DP=0.
      ZGL3DP=0.

      DO J=2,NCHI2+1
         H1 = SQRT(1.0-LAM/RHK(J))
         ZRG0PA(J) = H1*RJBK(J)
         ZRG1PA(J) = H1*ZRQ1K(J)
         ZRG2PA(J) = H1*ZRQ2K(J)
         ZRG3PA(J) = H1*ZRQ3K(J)
         
         H2 = 0.5*LAM/RHK(J)/H1
         ZRG0PE(J) = H2*RJBK(J)
         ZRG1PE(J) = H2*ZRQ1K(J)
         ZRG2PE(J) = H2*ZRQ2K(J)
         ZRG3PE(J) = H2*ZRQ3K(J)

         H2 = 0.5/H1
         ZRG0DP(J) = H2*RJBK(J)
         ZRG1DP(J) = H2*ZRQ1K(J)
         ZRG2DP(J) = H2*ZRQ2K(J)
         ZRG3DP(J) = H2*ZRQ3K(J)
      ENDDO
 
      ZRL(1)=M1-M2
      DO L = 2,MMLMAX           
      ZRL(L)=ZRL(L-1)+1
      ENDDO
      
      DO L=1,M2-M1+1
      
      IF(KPARTICLE.EQ.0) THEN

C     AT TURNING POINT LAM/H=1      
      Z0FL=0.5*SQRT(LAM)*RJBK(1)/SQRT(HPL)
      Z0FU=0.5*SQRT(LAM)*RJBK(NCHI2+2)/SQRT(-HPU)
      Z1FL=0.5*SQRT(LAM)*ZRQ1K(1)/SQRT(HPL)
      Z1FU=0.5*SQRT(LAM)*ZRQ1K(NCHI2+2)/SQRT(-HPU)
      Z2FL=0.5*SQRT(LAM)*ZRQ2K(1)/SQRT(HPL)
      Z2FU=0.5*SQRT(LAM)*ZRQ2K(NCHI2+2)/SQRT(-HPU)
      Z3FL=0.5*SQRT(LAM)*ZRQ3K(1)/SQRT(HPL)
      Z3FU=0.5*SQRT(LAM)*ZRQ3K(NCHI2+2)/SQRT(-HPU)

      ENDIF

      DO J=2,NCHI2+1
      
      EPHASE=EXP(CI*ZRL(L)*RCHIK(J))
      
      ZGL0PA(L)=ZGL0PA(L)+ZRG0PA(J)*EPHASE
      ZGL1PA(L)=ZGL1PA(L)+ZRG1PA(J)*EPHASE
      ZGL2PA(L)=ZGL2PA(L)+ZRG2PA(J)*EPHASE
      ZGL3PA(L)=ZGL3PA(L)+ZRG3PA(J)*EPHASE
      
      IF (KPARTICLE.EQ.1) THEN
       
         ZGL0PE(L)=ZGL0PE(L)+ZRG0PE(J)*EPHASE
         ZGL1PE(L)=ZGL1PE(L)+ZRG1PE(J)*EPHASE
         ZGL2PE(L)=ZGL2PE(L)+ZRG2PE(J)*EPHASE
         ZGL3PE(L)=ZGL3PE(L)+ZRG3PE(J)*EPHASE
  
         ZGL0DP(L)=ZGL0DP(L)+ZRG0DP(J)*EPHASE
         ZGL1DP(L)=ZGL1DP(L)+ZRG1DP(J)*EPHASE
         ZGL2DP(L)=ZGL2DP(L)+ZRG2DP(J)*EPHASE
         ZGL3DP(L)=ZGL3DP(L)+ZRG3DP(J)*EPHASE

      ELSEIF (KPARTICLE.EQ.0) THEN
        
C        FOR PRESSURE EQ
         CTMP0=EPHASE*(ZRG0PE(J)-Z0FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z0FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z0FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z0FU*SQRT(DIFPI(CHIU-RCHIK(J))))
        
         CTMP1=EPHASE*(ZRG1PE(J)-Z1FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z1FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z1FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z1FU*SQRT(DIFPI(CHIU-RCHIK(J))))
        
         CTMP2=EPHASE*(ZRG2PE(J)-Z2FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z2FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z2FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z2FU*SQRT(DIFPI(CHIU-RCHIK(J))))
        
         CTMP3=EPHASE*(ZRG3PE(J)-Z3FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z3FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z3FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z3FU*SQRT(DIFPI(CHIU-RCHIK(J))))
       
         ZGL0PE(L)=ZGL0PE(L)+CTMP0
         ZGL1PE(L)=ZGL1PE(L)+CTMP1
         ZGL2PE(L)=ZGL2PE(L)+CTMP2
         ZGL3PE(L)=ZGL3PE(L)+CTMP3

C        FOR DPHI EQ
         CTMP0=EPHASE*(ZRG0DP(J)-Z0FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z0FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z0FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z0FU*SQRT(DIFPI(CHIU-RCHIK(J))))
        
         CTMP1=EPHASE*(ZRG1DP(J)-Z1FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z1FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z1FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z1FU*SQRT(DIFPI(CHIU-RCHIK(J))))
        
         CTMP2=EPHASE*(ZRG2DP(J)-Z2FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z2FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z2FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z2FU*SQRT(DIFPI(CHIU-RCHIK(J))))
        
         CTMP3=EPHASE*(ZRG3DP(J)-Z3FL/SQRT(DIFPI(RCHIK(J)-CHIL))
     &                          -Z3FU/SQRT(DIFPI(CHIU-RCHIK(J)))
     &             -2.*CI*ZRL(L)*Z3FL*SQRT(DIFPI(RCHIK(J)-CHIL))
     &             +2.*CI*ZRL(L)*Z3FU*SQRT(DIFPI(CHIU-RCHIK(J))))
       
         ZGL0DP(L)=ZGL0DP(L)+CTMP0
         ZGL1DP(L)=ZGL1DP(L)+CTMP1
         ZGL2DP(L)=ZGL2DP(L)+CTMP2
         ZGL3DP(L)=ZGL3DP(L)+CTMP3

         IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0.AND.L.EQ.1) THEN
            I = INT(NLAMK0(JS,KGRID)/2.)
            IF (ABS(LAMM(I+1)-LAM).LT.1.0E-13) 
     &      WRITE(*,130) RCHIK(J),ZRG0PA(J),ZRG0PE(J),
     &                   ZRG0PE(J)*EPHASE,CTMP0
 130        FORMAT(7(E13.5,1X)) 
         ENDIF
      ENDIF
      ENDDO

      IF (KPARTICLE.EQ.0) THEN
         CTMP0    = 4.*SQRT(2.*PI+CHIU-CHIL)/RCHIHK
         CTMP1    = EXP(CI*ZRL(L)*CHIU)
         CTMP2    = EXP(CI*ZRL(L)*CHIL)

         ZGL0PE(L)= ZGL0PE(L)+(Z0FL*CTMP1+Z0FU*CTMP2)*CTMP0
         ZGL1PE(L)= ZGL1PE(L)+(Z1FL*CTMP1+Z1FU*CTMP2)*CTMP0
         ZGL2PE(L)= ZGL2PE(L)+(Z2FL*CTMP1+Z2FU*CTMP2)*CTMP0
         ZGL3PE(L)= ZGL3PE(L)+(Z3FL*CTMP1+Z3FU*CTMP2)*CTMP0

         ZGL0DP(L)= ZGL0DP(L)+(Z0FL*CTMP1+Z0FU*CTMP2)*CTMP0
         ZGL1DP(L)= ZGL1DP(L)+(Z1FL*CTMP1+Z1FU*CTMP2)*CTMP0
         ZGL2DP(L)= ZGL2DP(L)+(Z2FL*CTMP1+Z2FU*CTMP2)*CTMP0
         ZGL3DP(L)= ZGL3DP(L)+(Z3FL*CTMP1+Z3FU*CTMP2)*CTMP0
      ENDIF
      
      ENDDO
      
      DO L=1,M2-M1
         ZGL0PA(MMLMAX-L+1) = CONJG(ZGL0PA(L))
         ZGL0PE(MMLMAX-L+1) = CONJG(ZGL0PE(L))
         ZGL0DP(MMLMAX-L+1) = CONJG(ZGL0DP(L))
         ZGL1PA(MMLMAX-L+1) = CONJG(ZGL1PA(L))
         ZGL1PE(MMLMAX-L+1) = CONJG(ZGL1PE(L))
         ZGL1DP(MMLMAX-L+1) = CONJG(ZGL1DP(L))
         ZGL2PA(MMLMAX-L+1) = CONJG(ZGL2PA(L))
         ZGL2PE(MMLMAX-L+1) = CONJG(ZGL2PE(L))
         ZGL2DP(MMLMAX-L+1) = CONJG(ZGL2DP(L))
         ZGL3PA(MMLMAX-L+1) = CONJG(ZGL3PA(L))
         ZGL3PE(MMLMAX-L+1) = CONJG(ZGL3PE(L))
         ZGL3DP(MMLMAX-L+1) = CONJG(ZGL3DP(L))
      ENDDO

      H3=0.25/PI*RCHIHK
      IF(KGRID.EQ.1)ZDPSIDS=DPSIDS(JS)
      IF(KGRID.EQ.2)ZDPSIDS=DPSIDSM(JS)
      H4=H3*ZDPSIDS
            
      ZGL0PA=H3*ZGL0PA
      ZGL1PA=H4*ZGL1PA
      ZGL2PA=H4*ZGL2PA
      ZGL3PA=H3*ZGL3PA
      
      ZGL0PE=H3*ZGL0PE
      ZGL1PE=H4*ZGL1PE
      ZGL2PE=H4*ZGL2PE
      ZGL3PE=H3*ZGL3PE

      ZGL0DP=H3*ZGL0DP
      ZGL1DP=H4*ZGL1DP
      ZGL2DP=H4*ZGL2DP
      ZGL3DP=H3*ZGL3DP

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0.AND.
     &    KPARTICLE.EQ.0) THEN
         I = INT(NLAMK0(JS,KGRID)/2.)
         IF (ABS(LAMM(I+1)-LAM).LT.1.0E-13) THEN
         DO L=1,MMLMAX
            WRITE(*,140) ZRL(L),ZGL0PA(L),ZGL0PE(L)
         ENDDO
 140     FORMAT(I3,1X,4(E13.5,1X)) 
         ENDIF
      ENDIF

      DEALLOCATE(ZRG0PA,ZRG1PA,ZRG2PA,ZRG3PA,
     &           ZRG0PE,ZRG1PE,ZRG2PE,ZRG3PE,
     &           ZRG0DP,ZRG1DP,ZRG2DP,ZRG3DP,ZRL)

      RETURN
      END SUBROUTINE KG_ADIABATIC

C======================================================================
C   KG FOR ADIABATIC CONTRIBUTIONS WITH FOW
C   G_XI G_Q1 G_Q2 G_Q3           
C======================================================================
      SUBROUTINE KG_ADIABATIC1(JS,KGRID,KPARTICLE)
      
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      IMPLICIT NONE

      INTEGER JS,J,L,MMLMAX,KGRID,KPARTICLE,KCHECK
      
      REAL*8      H1,H2,H3,H4,ZDPSIDS
      INTEGER,DIMENSION(:),ALLOCATABLE::ZRL
      REAL*8,DIMENSION(:),ALLOCATABLE::ZRG0PA,ZRG1PA,ZRG2PA,ZRG3PA,
     &                                 ZRG0PE,ZRG1PE,ZRG2PE,ZRG3PE,
     &                                 ZRG0DP,ZRG1DP,ZRG2DP,ZRG3DP
      COMPLEX*16  EPHASE
      KCHECK = 0 

      MMLMAX = M2-M1-(M1-M2) + 1      

      ALLOCATE(
     &         ZRG0PA(NCHI2+2),ZRG1PA(NCHI2+2),
     &         ZRG2PA(NCHI2+2),ZRG3PA(NCHI2+2),
     &         ZRG0PE(NCHI2+2),ZRG1PE(NCHI2+2),
     &         ZRG2PE(NCHI2+2),ZRG3PE(NCHI2+2),
     &         ZRG0DP(NCHI2+2),ZRG1DP(NCHI2+2),
     &         ZRG2DP(NCHI2+2),ZRG3DP(NCHI2+2))
      ALLOCATE(ZRL(MMLMAX))

      ZGL0PA1=0.
      ZGL1PA1=0.
      ZGL2PA1=0.
      ZGL3PA1=0.
                  
      ZGL0PE1=0.
      ZGL1PE1=0.
      ZGL2PE1=0.
      ZGL3PE1=0.

      ZGL0DP1=0.
      ZGL1DP1=0.
      ZGL2DP1=0.
      ZGL3PE1=0.

      DO J=2,NCHI2+1
         H1 = RHK(J)-LAM
         ZRG0PA(J) = H1*RJBK(J)
         ZRG1PA(J) = H1*ZRQ1K(J)
         ZRG2PA(J) = H1*ZRQ2K(J)
         ZRG3PA(J) = H1*ZRQ3K(J)
         
         H2 = 0.5*LAM
         ZRG0PE(J) = H2*RJBK(J)
         ZRG1PE(J) = H2*ZRQ1K(J)
         ZRG2PE(J) = H2*ZRQ2K(J)
         ZRG3PE(J) = H2*ZRQ3K(J)
         
         H2 = 0.5*RHK(J)
         ZRG0DP(J) = H2*RJBK(J)
         ZRG1DP(J) = H2*ZRQ1K(J)
         ZRG2DP(J) = H2*ZRQ2K(J)
         ZRG3DP(J) = H2*ZRQ3K(J)
      ENDDO
 
      ZRL(1)=M1-M2
      DO L = 2,MMLMAX           
      ZRL(L)=ZRL(L-1)+1
      ENDDO
      
      DO L=1,M2-M1+1
      DO J=2,NCHI2+1
         EPHASE=EXP(CI*ZRL(L)*RCHIK(J))
      
         ZGL0PA1(L)=ZGL0PA1(L)+ZRG0PA(J)*EPHASE
         ZGL1PA1(L)=ZGL1PA1(L)+ZRG1PA(J)*EPHASE
         ZGL2PA1(L)=ZGL2PA1(L)+ZRG2PA(J)*EPHASE
         ZGL3PA1(L)=ZGL3PA1(L)+ZRG3PA(J)*EPHASE
      
         ZGL0PE1(L)=ZGL0PE1(L)+ZRG0PE(J)*EPHASE
         ZGL1PE1(L)=ZGL1PE1(L)+ZRG1PE(J)*EPHASE
         ZGL2PE1(L)=ZGL2PE1(L)+ZRG2PE(J)*EPHASE
         ZGL3PE1(L)=ZGL3PE1(L)+ZRG3PE(J)*EPHASE
      
         ZGL0DP1(L)=ZGL0DP1(L)+ZRG0DP(J)*EPHASE
         ZGL1DP1(L)=ZGL1DP1(L)+ZRG1DP(J)*EPHASE
         ZGL2DP1(L)=ZGL2DP1(L)+ZRG2DP(J)*EPHASE
         ZGL3DP1(L)=ZGL3DP1(L)+ZRG3DP(J)*EPHASE
      ENDDO
      ENDDO

      DO L=1,M2-M1
         ZGL0PA1(MMLMAX-L+1) = CONJG(ZGL0PA1(L))
         ZGL0PE1(MMLMAX-L+1) = CONJG(ZGL0PE1(L))
         ZGL0DP1(MMLMAX-L+1) = CONJG(ZGL0DP1(L))
         ZGL1PA1(MMLMAX-L+1) = CONJG(ZGL1PA1(L))
         ZGL1PE1(MMLMAX-L+1) = CONJG(ZGL1PE1(L))
         ZGL1DP1(MMLMAX-L+1) = CONJG(ZGL1DP1(L))
         ZGL2PA1(MMLMAX-L+1) = CONJG(ZGL2PA1(L))
         ZGL2PE1(MMLMAX-L+1) = CONJG(ZGL2PE1(L))
         ZGL2DP1(MMLMAX-L+1) = CONJG(ZGL2DP1(L))
         ZGL3PA1(MMLMAX-L+1) = CONJG(ZGL3PA1(L))
         ZGL3PE1(MMLMAX-L+1) = CONJG(ZGL3PE1(L))
         ZGL3DP1(MMLMAX-L+1) = CONJG(ZGL3DP1(L))
      ENDDO

      H3=0.25/PI*RCHIHK
      IF(KGRID.EQ.1)ZDPSIDS=DPSIDS(JS)
      IF(KGRID.EQ.2)ZDPSIDS=DPSIDSM(JS)
      H4=H3*ZDPSIDS
            
      ZGL0PA1=H3*ZGL0PA1
      ZGL1PA1=H4*ZGL1PA1
      ZGL2PA1=H4*ZGL2PA1
      ZGL3PA1=H3*ZGL3PA1
      
      ZGL0PE1=H3*ZGL0PE1
      ZGL1PE1=H4*ZGL1PE1
      ZGL2PE1=H4*ZGL2PE1
      ZGL3PE1=H3*ZGL3PE1

      ZGL0DP1=H3*ZGL0DP1
      ZGL1DP1=H4*ZGL1DP1
      ZGL2DP1=H4*ZGL2DP1
      ZGL3DP1=H3*ZGL3DP1

      DEALLOCATE(ZRG0PA,ZRG1PA,ZRG2PA,ZRG3PA,
     &           ZRG0PE,ZRG1PE,ZRG2PE,ZRG3PE,
     &           ZRG0DP,ZRG1DP,ZRG2DP,ZRG3DP,ZRL)

      RETURN
      END SUBROUTINE KG_ADIABATIC1

C=========================================================
C I-FACTOR FOR PRECESSION/BOUNCE RESONANCE OF TRAPPED EPS
C WITH GENERIC EQUILIBRIUM ANISOTROPIC DISTRIBUTION
C APPLYING NUMERICAL ENERGY INTEGRATION, WITH ANALYTIC
C EXTRACTION OF SINGULARITY (I.E. THE RESONANCE)
C KOPT = 0: VANISHING ORBIT WIDTH CONTRIBUTION
C        1: =0
C        2: FOW CORRECTION 2 FOR TRAPPED PARTICLES
C        3: FOW CORRECTION 3 FOR TRAPPED PARTICLES
C KDPHI = 1: I-FACTOR FOR PRESSURE EQ WITH DPHI VARIABLE
C         0: I-FACTOR FOR PRESSURE EQ WITH ALL OTHER VARIABLES
C         3: I-FACTOR FOR DPHI EQ WITH DPHI VARIABLE
C         2: I-FACTOR FOR DPHI EQ WITH ALL OTHER VARIABLES
C=========================================================
      SUBROUTINE KIA_TRAP0(JS,KGRID,KP,KOPT,L,ZZLAMB,KDPHI,ZVIF)
      
      USE RCOMDM
      USE DIMENSIM
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KP,KOPT,L,K,KCHECK,KEXEC,KDPHI
      REAL*8     RL,ZZLAMB,REPS,RTMP,RTMP1,RTMP2,HEPK,FTMP,RFAC
      REAL*8     H1,H2,H3,H4,H5,H6,H7,FAC0,OMEGAEDPSI,ZZT,OMEGANRE
      COMPLEX*16 ZVIF,OMEGAN,OMEGAN0,ZTMP0
      COMPLEX*16 FNUM
      PARAMETER  (REPS=1.0E-10)

      ZVIF  = 0.0
      RL    = RLM(L)

      IF (KOPT.EQ.1.OR.ISPECIES_EK(KP).EQ.0) RETURN

      IF (ABS(PSPECIES_FOWT(KP)).LT.REPS.AND.KOPT.GT.0) RETURN

C     FIRST ORDER FOW CONTRIBUTION VANISHES FOR PRECESSIONAL DRIFTS
      IF (ABS(RL).LT.0.1.AND.KOPT.GT.0) RETURN
      
      KCHECK = 0
      KEXEC  = 0
      K      = NEPK2
      HEPK   = ZEPKO(NEPK2)
      RFAC   = 1.

      IF (ISPECIES_F0(KP).EQ.5) THEN
      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)                            !=P0
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI  !=N0
      H6 = 8.1872E-14                                      !=me*c^2
      H5 = H1/H2/H6                                        !=T0/(me*c^2)
      H7 = ESPECIES_TEM(JS,KGRID,2)/ALPHAA1(JS,KGRID,KP)*H5!=EPS_A/(me*c^2)
C     RFAC = HEPK*H7+1.
      ENDIF

      H3 = NUEFF(JS,KGRID,KP)
      H5 = EPSALPHA(JS,KGRID,KP)
C     H4 = H5**(-1.5)
      H4 = 1.0       
      H6 = DEPSALPHADPSI(JS,KGRID,KP)

      OMEGAN0 = RNTOR*OMEGAE0(JS,KGRID)-OMEGA
      OMEGAN  = OMEGAN0-CI*H3

      RTMP1 = RNTOR*B0K/OMEGACI0*H5*ESPECIES_Z(1)/ESPECIES_Z(KP)

      IF (KGRID.EQ.1) THEN
         H1  = (CS(JS)-CS(JS-1))/2.
         H2  = (CS(JS+1)-CS(JS))/2.
         OMEGAEDPSI=((H1/H2*OMEGAE0(JS,2)-H2/H1*OMEGAE0(JS-1,2))
     &              /(H1+H2)-(H1-H2)*OMEGAE0(JS,1)/H1/H2)/DPSIDS(JS)
         ZZT = T(JS)
      ELSE
         H1  = CS(JS+1) - CS(JS)
         OMEGAEDPSI=(OMEGAE0(JS+1,1)-OMEGAE0(JS,1))/H1/DPSIDSM(JS)
         ZZT = TM(JS)
      ENDIF

C     COMPUTE FNUM=WHOLE NUMERATOR OF THE I-FACTOR
C     N0-FACTOR TIMES F_E: F_PSI/F_E IS 0/0 IN THE MAXWELLIAN TAIL
      FNUM = RTMP1*ZRESU(K,2,1,KP)+OMEGAN0*ZRESU(K,2,2,KP)

      IF (KOPT.EQ.0) THEN
         FAC0 =-2.
         FNUM = FAC0*HEPK**2.5*FNUM/RFAC**(1.5)
         KEXEC = 1
      ELSEIF (KOPT.EQ.2.AND.ABS(PSPECIES_FOWT(KP)).GT.0.) THEN
         FAC0 = PSPECIES_FOWT(KP)*
     &          2.*SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &          ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &          ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FTMP = RTMP1*(ZRESU(K,2,4,KP)+H6*ZRESU(K,2,1,KP))+
     &          RNTOR*OMEGAEDPSI*ZRESU(K,2,2,KP)+OMEGAN0*ZRESU(K,2,5,KP)
         FNUM = FAC0*HEPK**3/RFAC**2*(FTMP-H6*FNUM)
         KEXEC = 1
      ELSEIF (KOPT.EQ.3.AND.ABS(PSPECIES_FOWT(KP)).GT.0.) THEN
         FAC0 =-PSPECIES_FOWT(KP)*
     &          2.*SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &          ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &          ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FNUM = FAC0*HEPK**3*FNUM/RFAC**2
         KEXEC = 1
      ENDIF

      IF (KDPHI.EQ.1) THEN
         RTMP = ESPECIES_Z(KP)*OMEGACI0/B0K/H5
         FNUM  = FNUM*RTMP/HEPK
      ENDIF
      IF (KDPHI.EQ.2) THEN
         RTMP = ESPECIES_Z(KP)/H5
         FNUM  = FNUM*RTMP/HEPK
      ENDIF
      IF (KDPHI.EQ.3) THEN
         RTMP = (ESPECIES_Z(KP)/H5)**2*OMEGACI0/B0K
         FNUM  = FNUM*RTMP/HEPK**2
      ENDIF

      IF (ABS(RL).LT.0.1) THEN
         RTMP2 = RTMP1
         RTMP  = DRIFT*RTMP2
         IF (KEXEC.EQ.1.AND.ABS(RTMP).GT.REPS) ZVIF = FNUM/RTMP
      ELSE
         RTMP2 = RL*SQRT(2.*H5*ESPECIES_M(1)/ESPECIES_M(KP))
         RTMP  = OMEGAB*RTMP2
         IF (KEXEC.EQ.1.AND.ABS(RTMP).GT.REPS) ZVIF = 2.*FNUM/RTMP
      ENDIF

      RETURN
      END SUBROUTINE KIA_TRAP0

C=========================================================
C I-FACTOR FOR PRECESSION/BOUNCE RESONANCE OF TRAPPED EPS
C WITH GENERIC EQUILIBRIUM ANISOTROPIC DISTRIBUTION
C APPLYING NUMERICAL ENERGY INTEGRATION, WITH ANALYTIC
C EXTRACTION OF SINGULARITY (I.E. THE RESONANCE)
C KOPT = 0: VANISHING ORBIT WIDTH CONTRIBUTION
C        1: =0
C        2: FOW CORRECTION 2 FOR TRAPPED PARTICLES
C        3: FOW CORRECTION 3 FOR TRAPPED PARTICLES
C KDPHI = 1: I-FACTOR FOR PRESSURE EQ WITH DPHI VARIABLE
C         0: I-FACTOR FOR PRESSURE EQ WITH ALL OTHER VARIABLES
C         3: I-FACTOR FOR DPHI EQ WITH DPHI VARIABLE
C         2: I-FACTOR FOR DPHI EQ WITH ALL OTHER VARIABLES
C=========================================================
      SUBROUTINE KIA_TRAP(JS,KGRID,KP,KOPT,L,ZZLAMB,KDPHI,ZVIF)
      
      USE RCOMDM
      USE DIMENSIM
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KP,KOPT,L,K,KCHECK,KEXEC,KDPHI
      REAL*8     RL,ZZLAMB,REPS,RTMP,RTMP1,RTMP2,RTMP3,HEPK
      REAL*8     H1,H2,H3,H4,H5,H6,H7,FAC0,OMEGAEDPSI,ZZT,OMEGANRE
      COMPLEX*16 ZVIF,OMEGAN,OMEGAN0,FNUM0,ZTMP0
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::FNUM
      REAL*8,    DIMENSION(:),ALLOCATABLE::FTMP,EFAC,RFAC,RFAC2
      PARAMETER  (REPS=1.0E-10)
      REAL*8     ZEPKS,FNUM0R,FNUM0I

      ZVIF  = 0.0
      RL    = RLM(L)

      IF (KOPT.EQ.1) RETURN

      IF (ABS(PSPECIES_FOWT(KP)).LT.REPS.AND.KOPT.GT.0) RETURN

C     FIRST ORDER FOW CONTRIBUTION VANISHES FOR PRECESSIONAL DRIFTS
      IF (ABS(RL).LT.0.1.AND.KOPT.GT.0) RETURN
      
      KCHECK = 0
      KEXEC  = 0

      ALLOCATE( FNUM(NEPK2),FTMP(NEPK2),EFAC(NEPK2),RFAC(NEPK2),
     &          RFAC2(NEPK2) )

      IF (ISPECIES_EK(KP).EQ.0) ZEPK = ZEPKN
      IF (ISPECIES_EK(KP).EQ.1) ZEPK = ZEPKO

      RFAC  = 1.
      RFAC2 = 1.

      IF (ISPECIES_F0(KP).EQ.5) THEN
      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)                            !=P0
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI  !=N0
      H6 = 8.1872E-14                                      !=me*c^2
      H5 = H1/H2/H6                                        !=T0/(me*c^2)
      H7 = ESPECIES_TEM(JS,KGRID,2)/ALPHAA1(JS,KGRID,KP)*H5!=EPS_A/(me*c^2)
C     RFAC  = ZEPK*H7+1.
      IF (H7.LT.0.5) THEN
         RFAC2 = (1.0-2.0*ZEPK*H7)**(-0.5)
      ELSE
         STOP 'RFAC2'
      ENDIF
      ENDIF

      IF (ISPECIES_F0(KP).EQ.6) THEN
         RFAC2 = (1.0-ZEPK)**(-0.5)
         RFAC2(NEPK2)=RFAC2(NEPK2-1)
      ENDIF

      H3 = NUEFF(JS,KGRID,KP)
      H5 = EPSALPHA(JS,KGRID,KP)
C     NOTE NEED TO SET H4=1 TO RECOVER KFASTRUN=1 OPTION 
C     ENERGY DEPENDENT COLLISIONALITY MODEL (INUTYPE.GE.1)
C     AS WELL AS TO RECOVER SHAING'S NTV MODEL
C     H4 = H5**(-1.5)
      H4 = 1.0       
      H6 = DEPSALPHADPSI(JS,KGRID,KP)

      OMEGAN0 = RNTOR*OMEGAE0(JS,KGRID)-OMEGA
      OMEGAN  = OMEGAN0-CI*H3

      RTMP1 = RNTOR*B0K/OMEGACI0*H5*ESPECIES_Z(1)/ESPECIES_Z(KP)

      IF (KGRID.EQ.1) THEN
         H1  = (CS(JS)-CS(JS-1))/2.
         H2  = (CS(JS+1)-CS(JS))/2.
         OMEGAEDPSI=((H1/H2*OMEGAE0(JS,2)-H2/H1*OMEGAE0(JS-1,2))
     &              /(H1+H2)-(H1-H2)*OMEGAE0(JS,1)/H1/H2)/DPSIDS(JS)
         ZZT = T(JS)
      ELSE
         H1  = CS(JS+1) - CS(JS)
         OMEGAEDPSI=(OMEGAE0(JS+1,1)-OMEGAE0(JS,1))/H1/DPSIDSM(JS)
         ZZT = TM(JS)
      ENDIF

C     COMPUTE FNUM=WHOLE NUMERATOR OF THE I-FACTOR
C     N0-FACTOR TIMES F_E: F_PSI/F_E IS 0/0 IN THE MAXWELLIAN TAIL
      FNUM = RTMP1*ZRESU(:,2,1,KP)+OMEGAN0*ZRESU(:,2,2,KP)

      IF (KOPT.EQ.0) THEN
         FAC0 =-2.
         FNUM = FAC0*ZEPK**2.5*FNUM/RFAC**(1.5)
         KEXEC = 1
      ELSEIF (KOPT.EQ.2.AND.ABS(PSPECIES_FOWT(KP)).GT.0.) THEN
         FAC0 = PSPECIES_FOWT(KP)*
     &          2.*SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &          ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &          ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FTMP = RTMP1*(ZRESU(:,2,4,KP)+H6*ZRESU(:,2,1,KP))+
     &          RNTOR*OMEGAEDPSI*ZRESU(:,2,2,KP)+OMEGAN0*ZRESU(:,2,5,KP)
         FNUM = FAC0*ZEPK**3/RFAC**2*(FTMP-H6*FNUM)
         KEXEC = 1
      ELSEIF (KOPT.EQ.3.AND.ABS(PSPECIES_FOWT(KP)).GT.0.) THEN
         FAC0 =-PSPECIES_FOWT(KP)*
     &          2.*SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &          ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &          ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FNUM = FAC0*ZEPK**3*FNUM/RFAC**2
         KEXEC = 1
      ENDIF

      IF (KEXEC.EQ.0) GOTO 1001

C     EFAC=1 FOR PARTICLES WITH FINITE BIRTH ENERGY, AND
C     EFAC=PI/2/COS^2(ZEPK*PI/2) FOR PARTICLES WITH INFINITE BIRTH ENERGY 
      IF (ISPECIES_EK(KP).EQ.0) EFAC = 0.5*PI/(COS(0.5*PI*ZEPKO))**2

      IF (KDPHI.EQ.1) THEN
         RTMP = ESPECIES_Z(KP)*OMEGACI0/B0K/H5
         FNUM  = FNUM*RTMP/ZEPK
      ENDIF
      IF (KDPHI.EQ.2) THEN
         RTMP = ESPECIES_Z(KP)/H5
         FNUM  = FNUM*RTMP/ZEPK
      ENDIF
      IF (KDPHI.EQ.3) THEN
         RTMP = (ESPECIES_Z(KP)/H5)**2*OMEGACI0/B0K
         FNUM  = FNUM*RTMP/ZEPK**2
      ENDIF

C     NOW CONSIDER WHOLE RESONANCE OPERATOR
      IF (ABS(RL).LT.0.1) THEN
         RTMP2 = RTMP1
         RTMP  = DRIFT*RTMP2
      ELSE
         RTMP2 = RL*SQRT(2.*H5*ESPECIES_M(1)/ESPECIES_M(KP))
         RTMP  = OMEGAB*RTMP2 
         RTMP3 = PSPECIES_NDB(KP)*DRIFT*RTMP1/(OMEGAB*RTMP2)
      ENDIF

      IF (ABS(RTMP).LT.REPS) THEN
         IF (ABS(OMEGAN).LT.REPS) OMEGAN=REPS

         FTMP = 0.
         IF (INUTYPE.EQ.1) FTMP = H3*(1.-ZEPK**(-1.5)*H4*(1.+RL**2/4.))
         IF (INUTYPE.EQ.2) FTMP = H3*(1.-ZEPK**(-1.5)*H4)

         IF (ISPECIES_EK(KP).EQ.0) FNUM = FNUM*EFAC
         FNUM = FNUM/(OMEGAN+CI*FTMP)
      ELSE
         FNUM     = FNUM/RTMP
         OMEGAN   = OMEGAN/RTMP
         OMEGANRE = DREAL(OMEGAN)
          
         IF (INUTYPE.EQ.0.AND.(
     &      (ISPECIES_EK(KP).EQ.0.AND.OMEGANRE.LT.0.).OR.
     &      (ISPECIES_EK(KP).EQ.1.AND.OMEGANRE.LT.0.AND.OMEGANRE.GT.-1.)
     &      )) THEN 
C        NOTE THAT THIS OPTION HOLDS ONLY FOR ENERGY-INDEPENDENT COLLISIONALITY
         IF (ISPECIES_EK(KP).EQ.0) THEN
            IF (ABS(RL).LT.0.1) FNUM = FNUM*(ZEPK*RFAC2-OMEGANRE)
            IF (ABS(RL).GT.0.1) FNUM = FNUM*(SQRT(ZEPK)-OMEGAN)*
     &                                 (ZEPK+OMEGANRE**2)
         ENDIF
      
         IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).LT.0.1) 
     &      ZEPKS = 2.*ATAN(-OMEGANRE)/PI
         IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).GT.0.1) 
     &      ZEPKS = 2.*ATAN(OMEGANRE**2)/PI
         IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).LT.0.1) ZEPKS=-OMEGANRE
         IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).GT.0.1) ZEPKS=OMEGANRE**2

         IF (ISPECIES_EK(KP).EQ.0) FNUM(NEPK2) = 0.
 
         IF (ZEPKS.LE.ZEPKO(1)) THEN
            FNUM0 = FNUM(1)*ZEPKS/ZEPKO(1)
         ELSE  
            CALL SPLINE1D(FNUM0R,ZEPKS,1,DREAL(FNUM),ZEPKO,NEPK2,FTMP)
            CALL SPLINE1D(FNUM0I,ZEPKS,1,DIMAG(FNUM),ZEPKO,NEPK2,FTMP)
            FNUM0 = CMPLX(FNUM0R,FNUM0I)
         ENDIF

         IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).LT.0.1) 
     &      FNUM = (FNUM-FNUM0)/(ZEPK+OMEGAN)/(ZEPK-OMEGANRE)*EFAC
         IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).GT.0.1) 
     &      FNUM = (FNUM-FNUM0)/(ZEPK-OMEGAN**2)/(ZEPK+OMEGANRE**2)*EFAC
         IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).LT.0.1) 
     &      FNUM = (FNUM-FNUM0)/(ZEPK+OMEGAN)
         IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).GT.0.1) 
     &      FNUM = (FNUM-FNUM0)/(SQRT(ZEPK)+OMEGAN)

         IF (ABS(DIMAG(OMEGAN)).GT.REPS) THEN
            IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).LT.0.1) 
     &         ZTMP0 = CLOG(-OMEGAN/OMEGANRE)/(OMEGAN+OMEGANRE)
            IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).GT.0.1) 
     &         ZTMP0 = CLOG(-OMEGANRE**2/OMEGAN**2)/
     &                 (OMEGAN**2+OMEGANRE**2)
            IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).LT.0.1) 
     &         ZTMP0 = CLOG((1.+OMEGAN)/OMEGAN)
            IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).GT.0.1) 
     &         ZTMP0 = 2.-2.*OMEGAN*CLOG((1.+OMEGAN)/OMEGAN)
         ELSE
            IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).LT.0.1) 
     &         ZTMP0 =-CI*PI*RTMP/ABS(RTMP)/2./OMEGANRE 
            IF (ISPECIES_EK(KP).EQ.0.AND.ABS(RL).GT.0.1) 
     &         ZTMP0 = CI*PI*RTMP/ABS(RTMP)/2./OMEGANRE**2 
            IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).LT.0.1) 
     &         ZTMP0 =LOG(-(1+OMEGANRE)/OMEGANRE)+CI*PI*RTMP/ABS(RTMP)
            IF (ISPECIES_EK(KP).EQ.1.AND.ABS(RL).GT.0.1) 
     &         ZTMP0 = 2.-2.*OMEGANRE*(
     &                 LOG(-(1+OMEGANRE)/OMEGANRE)+CI*PI*RTMP/ABS(RTMP))
         ENDIF

         ZVIF = ZVIF + ZTMP0*FNUM0

         ELSE

         FTMP = 0.
         IF (INUTYPE.EQ.1) FTMP = H3*(1.-ZEPK**(-1.5)*H4
     &                            *(1.+RL**2/4.))/RTMP
         IF (INUTYPE.EQ.2) FTMP = H3*(1.-ZEPK**(-1.5)*H4)/RTMP

         IF (ISPECIES_EK(KP).EQ.0) FNUM = FNUM*EFAC
         IF (ABS(RL).LT.0.1) FNUM = FNUM/(ZEPK*RFAC2+OMEGAN+CI*FTMP)
         IF (ABS(RL).GT.0.1) FNUM = FNUM/(SQRT(ZEPK)+ZEPK*RTMP3*RFAC2+
     &                                    OMEGAN+CI*FTMP)
         ENDIF

C        SUBTRACT SINGULAR CONTRIBUTION FOR LATER PITCH ANGLE INTEGRATION
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.0) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,1)
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.1) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,2)
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.2) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,3)
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.3) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,4)

      ENDIF
  
      DO K=1,NEPK-1
         HEPK = EPK(K+1)-EPK(K)
         ZVIF = ZVIF + 0.5*HEPK*(FNUM(2*K-1)+FNUM(2*K))
      ENDDO

 1001 CONTINUE
      DEALLOCATE(FNUM,FTMP,EFAC,RFAC,RFAC2)

      IF (KCHECK.EQ.1.AND.JS.EQ.JSOUT) THEN
         WRITE(*,*) 'KIA_TRAP:',JS,KGRID,KP,KOPT,L,ZZLAMB,KDPHI,ZVIF
      ENDIF

      RETURN
      END SUBROUTINE KIA_TRAP
 
C=========================================================
C SF0-FACTOR ASSOCIATED WITH SINGULARITY SUBTRACTION
C FOR TRANSIT RESONANCE OF PASSING EPS
C WITH GENERIC EQUILIBRIUM ANISOTROPIC DISTRIBUTION
C APPLYING NUMERICAL ENERGY INTEGRATION, WITH ANALYTIC
C EXTRACTION OF SINGULARITY (I.E. THE RESONANCE)
C KOPT = 0: VANISHING ORBIT WIDTH CONTRIBUTION
C        1: FOW CORRECTION 1 FOR PASSING PARTICLES
C        2: FOW CORRECTION 2 FOR PASSING PARTICLES
C        3: FOW CORRECTION 3 FOR PASSING PARTICLES
C=========================================================
      SUBROUTINE KIA_PASS0(JS,KGRID,KP,KOPT,L,ZZLAMB,KDPHI,ZVIF)
      
      USE RCOMDM
      USE DIMENSIM
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KP,KOPT,L,K,KCHECK,KEXEC,KDPHI
      REAL*8     RL,ZZLAMB,REPS,RTMP,RTMP1,RTMP2,HEPK,FTMP1,FTMP3
      REAL*8     H1,H2,H3,H4,H5,H6,H7,FAC0,OMEGAEDPSI,ZZT,RFAC,
     &           PSI0L,PSI0DPSIL,PSI0DLAML,OMEGAN1RE,OMEGAN3RE
      COMPLEX*16 ZVIF,OMEGAN0,OMEGAN,ZTMP0,OMEGAN1,OMEGAN3
      COMPLEX*16 FNUM1,FNUM3
      PARAMETER  (REPS=1.0E-10)

      ZVIF  = 0.0
      RL    = RLM(L)

      IF (KOPT.EQ.1.AND.IFOWPSI0.EQ.2) RETURN
      IF (ISPECIES_EK(KP).EQ.0)        RETURN
      IF (ABS(PSPECIES_FOWP(KP)).LT.REPS.AND.KOPT.GT.0) RETURN

      KCHECK = 0
      KEXEC  = 0
      K      = NEPK2
      HEPK   = ZEPKO(NEPK2)
      RFAC   = 1.

      IF (ISPECIES_F0(KP).EQ.5) THEN
      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)                            !=P0
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI  !=N0
      H6 = 8.1872E-14                                      !=me*c^2
      H5 = H1/H2/H6                                        !=T0/(me*c^2)
      H7 = ESPECIES_TEM(JS,KGRID,2)/ALPHAA1(JS,KGRID,KP)*H5!=EPS_A/(me*c^2)
C     RFAC = HEPK*H7+1.
      ENDIF

      IF (KOPT.GT.0) THEN
      IF (KP.EQ.1.OR.KP.EQ.2) THEN
         H7 = ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &        (ESPECIES_M(KP)/ESPECIES_M(1))**2
         PSI0L     = H7*TPSI0L
         PSI0DPSIL = H7*TPSI0DPSIL
         PSI0DLAML = H7*TPSI0DLAML
      ELSE
         H7 = ESPECIES_Z(1)/ESPECIES_Z(KP)*ESPECIES_M(KP)/ESPECIES_M(1)
         PSI0L     = H7*HPSI0L
         PSI0DPSIL = H7*HPSI0DPSIL
         PSI0DLAML = H7*HPSI0DLAML
      ENDIF
      ENDIF

      H3 = NUEFF(JS,KGRID,KP)
      H5 = EPSALPHA(JS,KGRID,KP)
      H4 = H5**(-1.5)
      H6 = DEPSALPHADPSI(JS,KGRID,KP)

      OMEGAN0 = RNTOR*OMEGAE0(JS,KGRID)-OMEGA
      OMEGAN  = OMEGAN0-CI*H3

      RTMP1 = RNTOR*B0K/OMEGACI0*H5*ESPECIES_Z(1)/ESPECIES_Z(KP)

      IF (KGRID.EQ.1) THEN
         H1 = (CS(JS)-CS(JS-1))/2.
         H2 = (CS(JS+1)-CS(JS))/2.
         OMEGAEDPSI=((H1/H2*OMEGAE0(JS,2)-H2/H1*OMEGAE0(JS-1,2))
     &              /(H1+H2)-(H1-H2)*OMEGAE0(JS,1)/H1/H2)/DPSIDS(JS)
         ZZT = T(JS)
      ELSE
         H1 = CS(JS+1) - CS(JS)
         OMEGAEDPSI=(OMEGAE0(JS+1,1)-OMEGAE0(JS,1))/H1/DPSIDSM(JS)
         ZZT = TM(JS)
      ENDIF

C     COMPUTE FNUM=WHOLE NUMERATOR OF THE I-FACTOR
C     N0-FACTOR TIMES F_E: F_PSI/F_E IS 0/0 IN THE MAXWELLIAN TAIL
      FNUM1 = RTMP1*ZRESU(K,1,1,KP)+OMEGAN0*ZRESU(K,1,2,KP)
      FNUM3 = RTMP1*ZRESU(K,3,1,KP)+OMEGAN0*ZRESU(K,3,2,KP)

      IF (KOPT.EQ.0) THEN
         FAC0  =-1.
         FNUM1 = FAC0*HEPK**2.5*FNUM1/RFAC**(1.5)
         FNUM3 = FAC0*HEPK**2.5*FNUM3/RFAC**(1.5)
         KEXEC = 1
      ELSEIF (KOPT.EQ.1.AND.ABS(PSPECIES_FOWP(KP)).GT.0.
     &        .AND.IFOWPSI0.EQ.1) THEN
         FAC0  =-PSPECIES_FOWP(KP)*PSI0L*
     &           SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)
         FTMP1 = RTMP1*(ZRESU(K,1,4,KP)+H6*ZRESU(K,1,1,KP))
     &           +RNTOR*OMEGAEDPSI*ZRESU(K,1,2,KP)
     &           +OMEGAN0*ZRESU(K,1,5,KP)
         FTMP3 = RTMP1*(ZRESU(K,3,4,KP)+H6*ZRESU(K,3,1,KP))
     &           +RNTOR*OMEGAEDPSI*ZRESU(K,3,2,KP)
     &           +OMEGAN0*ZRESU(K,3,5,KP)
         FNUM1 = FAC0*HEPK**3/RFAC**2*
     &           (0.5/RTMP1/HEPK*(FNUM1+OMEGA*ZRESU(K,1,2,KP))*
     &           (2.*RTMP1*HEPK*PSI0DPSIL/PSI0L+OMEGAN0)+
     &           FTMP1-H6*FNUM1)
         FNUM3 =-FAC0*HEPK**3/RFAC**2*
     &           (0.5/RTMP1/HEPK*(FNUM3+OMEGA*ZRESU(K,3,2,KP))*
     &           (2.*RTMP1*HEPK*PSI0DPSIL/PSI0L+OMEGAN0)+
     &           FTMP3-H6*FNUM3)
         KEXEC = 1
      ELSEIF (KOPT.EQ.2.AND.ABS(PSPECIES_FOWP(KP)).GT.0.) THEN
         FAC0  = PSPECIES_FOWP(KP)*
     &           SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &           ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &           ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FTMP1 =RTMP1*(ZRESU(K,1,4,KP)+H6*ZRESU(K,1,1,KP))+
     &          RNTOR*OMEGAEDPSI*ZRESU(K,1,2,KP)+OMEGAN0*ZRESU(K,1,5,KP)
         FTMP3 =RTMP1*(ZRESU(K,3,4,KP)+H6*ZRESU(K,3,1,KP))+
     &          RNTOR*OMEGAEDPSI*ZRESU(K,3,2,KP)+OMEGAN0*ZRESU(K,3,5,KP)
         FNUM1 = FAC0*HEPK**3/RFAC**2*(FTMP1-H6*FNUM1)
         FNUM3 =-FAC0*HEPK**3/RFAC**2*(FTMP3-H6*FNUM3)
         KEXEC = 1
      ELSEIF (KOPT.EQ.3.AND.ABS(PSPECIES_FOWP(KP)).GT.0.) THEN
         FAC0  =-PSPECIES_FOWP(KP)*
     &           SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &           ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &           ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FNUM1 = FAC0*HEPK**3*FNUM1/RFAC**2
         FNUM3 =-FAC0*HEPK**3*FNUM3/RFAC**2
         KEXEC = 1
      ENDIF

      RTMP2 = (RNTOR*RQK+RL)*SQRT(2.*H5*ESPECIES_M(1)/ESPECIES_M(KP))
      RTMP  = OMEGAB*RTMP2

      IF (KEXEC.EQ.1.AND.ABS(RTMP).GT.REPS) THEN
         FNUM1   = FNUM1/RTMP
         FNUM3   = FNUM3/(-RTMP)
         OMEGAN1 = OMEGAN/RTMP
         OMEGAN3 = OMEGAN/(-RTMP)
         OMEGAN1RE = DREAL(OMEGAN1)
         OMEGAN3RE = DREAL(OMEGAN3)

         IF (KDPHI.EQ.1) THEN
            RTMP  = ESPECIES_Z(KP)*OMEGACI0/B0K/H5
            FNUM1 = FNUM1*RTMP/HEPK
            FNUM3 = FNUM3*RTMP/HEPK
         ENDIF
         IF (KDPHI.EQ.2) THEN
            RTMP  = ESPECIES_Z(KP)/H5
            FNUM1 = FNUM1*RTMP/HEPK
            FNUM3 = FNUM3*RTMP/HEPK
         ENDIF
         IF (KDPHI.EQ.3) THEN
            RTMP  = (ESPECIES_Z(KP)/H5)**2*OMEGACI0/B0K
            FNUM1 = FNUM1*RTMP/HEPK**2
            FNUM3 = FNUM3*RTMP/HEPK**2
         ENDIF

         IF (ABS(OMEGAN1RE+1.).LT.1.0E-5) THEN
            ZVIF = 2.*FNUM1
         ELSEIF (ABS(OMEGAN3RE+1.).LT.1.0E-5) THEN
            ZVIF = 2.*FNUM3
         ELSE
            WRITE(*,*) 'JS,KGRID,KP,KOPT,RL,ZZLAMB='
            WRITE(*,*) JS,KGRID,KP,KOPT,RLM(L),ZZLAMB
            WRITE(*,*) OMEGAN1RE,OMEGAN3RE
            STOP 'KIA_PASS0'
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE KIA_PASS0
 
C=========================================================
C I-FACTOR FOR TRANSIT RESONANCE OF PASSING EPS
C WITH GENERIC EQUILIBRIUM ANISOTROPIC DISTRIBUTION
C APPLYING NUMERICAL ENERGY INTEGRATION, WITH ANALYTIC
C EXTRACTION OF SINGULARITY (I.E. THE RESONANCE)
C KOPT = 0: VANISHING ORBIT WIDTH CONTRIBUTION
C        1: FOW CORRECTION 1 FOR PASSING PARTICLES
C        2: FOW CORRECTION 2 FOR PASSING PARTICLES
C        3: FOW CORRECTION 3 FOR PASSING PARTICLES
C=========================================================
      SUBROUTINE KIA_PASS(JS,KGRID,KP,KOPT,L,ZZLAMB,KDPHI,ZVIF)
      
      USE RCOMDM
      USE DIMENSIM
      USE KINETICM
      USE GLOBALM
      USE ANISOTROPICM
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KP,KOPT,L,K,KCHECK,KEXEC,KDPHI
      REAL*8     RL,ZZLAMB,REPS,FNUM0R,RTMP,RTMP1,RTMP2,HEPK,FNUM0I
      REAL*8     H1,H2,H3,H4,H5,H6,H7,FAC0,OMEGAEDPSI,ZZT,
     &           PSI0L,PSI0DPSIL,PSI0DLAML,OMEGAN1RE,ZEPKS,OMEGAN3RE
      COMPLEX*16 ZVIF,OMEGAN0,OMEGAN,FNUM10,FNUM30,ZTMP0,OMEGAN1,OMEGAN3
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::FNUM1,FNUM3
      REAL*8,    DIMENSION(:),ALLOCATABLE::FTMP1,FTMP3,EFAC,RFAC
      PARAMETER  (REPS=1.0E-10)

      ZVIF  = 0.0
      RL    = RLM(L)

      IF (KOPT.EQ.1.AND.IFOWPSI0.EQ.2) RETURN

      IF (ABS(PSPECIES_FOWP(KP)).LT.REPS.AND.KOPT.GT.0) RETURN

      KCHECK = 0
      KEXEC  = 0

      ALLOCATE( FNUM1(NEPK2),FNUM3(NEPK2),FTMP1(NEPK2),FTMP3(NEPK2),
     &          EFAC(NEPK2),RFAC(NEPK2) )

      IF (ISPECIES_EK(KP).EQ.0) ZEPK = ZEPKN
      IF (ISPECIES_EK(KP).EQ.1) ZEPK = ZEPKO

      RFAC = 1.

      IF (ISPECIES_F0(KP).EQ.5) THEN
      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)                            !=P0
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI  !=N0
      H6 = 8.1872E-14                                      !=me*c^2
      H5 = H1/H2/H6                                        !=T0/(me*c^2)
      H7 = ESPECIES_TEM(JS,KGRID,2)/ALPHAA1(JS,KGRID,KP)*H5!=EPS_A/(me*c^2)
C     RFAC = ZEPK*H7+1.
      ENDIF

      IF (KOPT.GT.0) THEN
      IF (KP.EQ.1.OR.KP.EQ.2) THEN
         H7 = ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &        (ESPECIES_M(KP)/ESPECIES_M(1))**2
         PSI0L     = H7*TPSI0L
         PSI0DPSIL = H7*TPSI0DPSIL
         PSI0DLAML = H7*TPSI0DLAML
      ELSE
         H7 = ESPECIES_Z(1)/ESPECIES_Z(KP)*ESPECIES_M(KP)/ESPECIES_M(1)
         PSI0L     = H7*HPSI0L
         PSI0DPSIL = H7*HPSI0DPSIL
         PSI0DLAML = H7*HPSI0DLAML
      ENDIF
      ENDIF

      H3 = NUEFF(JS,KGRID,KP)
      H5 = EPSALPHA(JS,KGRID,KP)
      H4 = H5**(-1.5)
      H6 = DEPSALPHADPSI(JS,KGRID,KP)

      OMEGAN0 = RNTOR*OMEGAE0(JS,KGRID)-OMEGA
      OMEGAN  = OMEGAN0-CI*H3

      RTMP1 = RNTOR*B0K/OMEGACI0*H5*ESPECIES_Z(1)/ESPECIES_Z(KP)

      IF (KGRID.EQ.1) THEN
         H1 = (CS(JS)-CS(JS-1))/2.
         H2 = (CS(JS+1)-CS(JS))/2.
         OMEGAEDPSI=((H1/H2*OMEGAE0(JS,2)-H2/H1*OMEGAE0(JS-1,2))
     &              /(H1+H2)-(H1-H2)*OMEGAE0(JS,1)/H1/H2)/DPSIDS(JS)
         ZZT = T(JS)
      ELSE
         H1 = CS(JS+1) - CS(JS)
         OMEGAEDPSI=(OMEGAE0(JS+1,1)-OMEGAE0(JS,1))/H1/DPSIDSM(JS)
         ZZT = TM(JS)
      ENDIF

C     COMPUTE FNUM=WHOLE NUMERATOR OF THE I-FACTOR
C     N0-FACTOR TIMES F_E: F_PSI/F_E IS 0/0 IN THE MAXWELLIAN TAIL
      FNUM1 = RTMP1*ZRESU(:,1,1,KP)+OMEGAN0*ZRESU(:,1,2,KP)
      FNUM3 = RTMP1*ZRESU(:,3,1,KP)+OMEGAN0*ZRESU(:,3,2,KP)

      IF (KOPT.EQ.0) THEN
         FAC0  =-1.
         FNUM1 = FAC0*ZEPK**2.5*FNUM1/RFAC**(1.5)
         FNUM3 = FAC0*ZEPK**2.5*FNUM3/RFAC**(1.5)
         KEXEC = 1
      ELSEIF (KOPT.EQ.1.AND.ABS(PSPECIES_FOWP(KP)).GT.0.
     &        .AND.IFOWPSI0.EQ.1) THEN
         FAC0  =-PSPECIES_FOWP(KP)*PSI0L*
     &           SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)
         FTMP1 = RTMP1*(ZRESU(:,1,4,KP)+H6*ZRESU(:,1,1,KP))
     &           +RNTOR*OMEGAEDPSI*ZRESU(:,1,2,KP)
     &           +OMEGAN0*ZRESU(:,1,5,KP)
         FTMP3 = RTMP1*(ZRESU(:,3,4,KP)+H6*ZRESU(:,3,1,KP))
     &           +RNTOR*OMEGAEDPSI*ZRESU(:,3,2,KP)
     &           +OMEGAN0*ZRESU(:,3,5,KP)
         FNUM1 = FAC0*ZEPK**3/RFAC**2*
     &           (0.5/RTMP1/ZEPK*(FNUM1+OMEGA*ZRESU(:,1,2,KP))*
     &           (2.*RTMP1*ZEPK*PSI0DPSIL/PSI0L+OMEGAN0)+
     &           FTMP1-H6*FNUM1)
         FNUM3 =-FAC0*ZEPK**3/RFAC**2*
     &           (0.5/RTMP1/ZEPK*(FNUM3+OMEGA*ZRESU(:,3,2,KP))*
     &           (2.*RTMP1*ZEPK*PSI0DPSIL/PSI0L+OMEGAN0)+
     &           FTMP3-H6*FNUM3)
         KEXEC = 1
      ELSEIF (KOPT.EQ.2.AND.ABS(PSPECIES_FOWP(KP)).GT.0.) THEN
         FAC0  = PSPECIES_FOWP(KP)*
     &           SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &           ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &           ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FTMP1 =RTMP1*(ZRESU(:,1,4,KP)+H6*ZRESU(:,1,1,KP))+
     &          RNTOR*OMEGAEDPSI*ZRESU(:,1,2,KP)+OMEGAN0*ZRESU(:,1,5,KP)
         FTMP3 =RTMP1*(ZRESU(:,3,4,KP)+H6*ZRESU(:,3,1,KP))+
     &          RNTOR*OMEGAEDPSI*ZRESU(:,3,2,KP)+OMEGAN0*ZRESU(:,3,5,KP)
         FNUM1 = FAC0*ZEPK**3/RFAC**2*(FTMP1-H6*FNUM1)
         FNUM3 =-FAC0*ZEPK**3/RFAC**2*(FTMP3-H6*FNUM3)
         KEXEC = 1
      ELSEIF (KOPT.EQ.3.AND.ABS(PSPECIES_FOWP(KP)).GT.0.) THEN
         FAC0  =-PSPECIES_FOWP(KP)*
     &           SQRT(2.*ESPECIES_M(1)/ESPECIES_M(KP)*H5)*
     &           ESPECIES_Z(1)/ESPECIES_Z(KP)*
     &           ESPECIES_M(KP)/ESPECIES_M(1)*ZZT/OMEGACI0
         FNUM1 = FAC0*ZEPK**3*FNUM1/RFAC**2
         FNUM3 =-FAC0*ZEPK**3*FNUM3/RFAC**2
         KEXEC = 1
      ENDIF

      IF (KEXEC.EQ.0) GOTO 1001

C     EFAC=1 FOR PARTICLES WITH FINITE BIRTH ENERGY, AND
C     EFAC=PI/2/COS^2(ZEPK*PI/2) FOR PARTICLES WITH INFINITE BIRTH ENERGY 
      IF (ISPECIES_EK(KP).EQ.0) EFAC = 0.5*PI/(COS(0.5*PI*ZEPKO))**2

      IF (KDPHI.EQ.1) THEN
         RTMP  = ESPECIES_Z(KP)*OMEGACI0/B0K/H5
         FNUM1 = FNUM1*RTMP/ZEPK
         FNUM3 = FNUM3*RTMP/ZEPK
      ENDIF
      IF (KDPHI.EQ.2) THEN
         RTMP  = ESPECIES_Z(KP)/H5
         FNUM1 = FNUM1*RTMP/ZEPK
         FNUM3 = FNUM3*RTMP/ZEPK
      ENDIF
      IF (KDPHI.EQ.3) THEN
         RTMP  = (ESPECIES_Z(KP)/H5)**2*OMEGACI0/B0K
         FNUM1 = FNUM1*RTMP/ZEPK**2
         FNUM3 = FNUM3*RTMP/ZEPK**2
      ENDIF

C     NOW CONSIDER WHOLE RESONANCE OPERATOR
      RTMP2 = (RNTOR*RQK+RL)*SQRT(2.*H5*ESPECIES_M(1)/ESPECIES_M(KP))
      RTMP  = OMEGAB*RTMP2

      IF (ABS(RTMP).LT.REPS) THEN
         IF (ABS(OMEGAN).LT.REPS) OMEGAN=REPS
 
         FTMP1 = 0.
         IF (INUTYPE.EQ.1) FTMP1 = H3*(1.-ZEPK**(-1.5)*H4*(1.+RL**2/4.))
         IF (INUTYPE.EQ.2) FTMP1 = H3*(1.-ZEPK**(-1.5)*H4)

         IF (ISPECIES_EK(KP).EQ.0) THEN
            FNUM1 = FNUM1*EFAC
            FNUM3 = FNUM3*EFAC
         ENDIF
         FNUM1 = FNUM1/(OMEGAN+CI*FTMP1)
         FNUM3 = FNUM3/(OMEGAN+CI*FTMP1)
         KEXEC = 1
      ELSE
         FNUM1   = FNUM1/RTMP*SQRT(RFAC)
         FNUM3   = FNUM3/(-RTMP)*SQRT(RFAC)
         OMEGAN1 = OMEGAN/RTMP
         OMEGAN3 = OMEGAN/(-RTMP)
         OMEGAN1RE = DREAL(OMEGAN1)
         OMEGAN3RE = DREAL(OMEGAN3)

C        INTEGRATION FOR CO-PASSING PARTICLES
         IF (INUTYPE.EQ.0.AND.(
     &       (ISPECIES_EK(KP).EQ.0.AND.OMEGAN1RE.LT.0.).OR.
     &       (ISPECIES_EK(KP).EQ.1.AND.OMEGAN1RE.LT.0.AND.
     &        OMEGAN1RE.GT.-1.))) THEN
C        NOTE THAT THIS OPTION HOLDS ONLY FOR ENERGY-INDEPENDENT COLLISIONALITY
         IF (ISPECIES_EK(KP).EQ.0) 
     &      FNUM1 = FNUM1*(SQRT(ZEPK)-OMEGAN1)*
     &              (ZEPK+OMEGAN1RE**2)

         IF (ISPECIES_EK(KP).EQ.0) ZEPKS = 2.*ATAN(OMEGAN1RE**2)/PI
         IF (ISPECIES_EK(KP).EQ.1) ZEPKS = OMEGAN1RE**2
         IF (ISPECIES_EK(KP).EQ.0) FNUM1(NEPK2) = 0.
 
         IF (ZEPKS.LE.ZEPKO(1)) THEN
            FNUM10 = FNUM1(1)*ZEPKS/ZEPKO(1)
         ELSE  
            CALL SPLINE1D(FNUM0R,ZEPKS,1,DREAL(FNUM1),ZEPKO,NEPK2,FTMP1)
            CALL SPLINE1D(FNUM0I,ZEPKS,1,DIMAG(FNUM1),ZEPKO,NEPK2,FTMP1)
            FNUM10 = CMPLX(FNUM0R,FNUM0I)
         ENDIF

         IF (ISPECIES_EK(KP).EQ.0) FNUM1 = (FNUM1-FNUM10)/
     &      (ZEPK-OMEGAN1**2)/(ZEPK+OMEGAN1RE**2)*EFAC
         IF (ISPECIES_EK(KP).EQ.1) FNUM1 = (FNUM1-FNUM10)/
     &      (SQRT(ZEPK)+OMEGAN1)

         IF (ABS(DIMAG(OMEGAN1)).GT.REPS) THEN
            IF (ISPECIES_EK(KP).EQ.0) 
     &         ZTMP0 = CLOG(-OMEGAN1RE**2/OMEGAN1**2)/
     &                 (OMEGAN1**2+OMEGAN1RE**2)
            IF (ISPECIES_EK(KP).EQ.1) 
     &         ZTMP0 = 2.-2.*OMEGAN1*CLOG((1.+OMEGAN1)/OMEGAN1)
         ELSE
            IF (ISPECIES_EK(KP).EQ.0) 
     &         ZTMP0 = CI*PI*RTMP/ABS(RTMP)/2./OMEGAN1RE**2 
            IF (ISPECIES_EK(KP).EQ.1) ZTMP0 = 2.-2.*OMEGAN1RE*(
     &         LOG(-(1+OMEGAN1RE)/OMEGAN1RE)+CI*PI*RTMP/ABS(RTMP))
         ENDIF

         ZVIF = ZVIF + ZTMP0*FNUM10

         ELSE
         
         FTMP1 = 0.
         IF (INUTYPE.EQ.1) FTMP1 = H3*(1.-ZEPK**(-1.5)*H4*
     &                             (1.+RL**2/4.))/RTMP
         IF (INUTYPE.EQ.2) FTMP1 = H3*(1.-ZEPK**(-1.5)*H4)/RTMP
   
         IF (ISPECIES_EK(KP).EQ.0) FNUM1 = FNUM1*EFAC
         FNUM1 = FNUM1/(SQRT(ZEPK/RFAC)+OMEGAN1+CI*FTMP1)

         ENDIF

C        INTEGRATION FOR COUNTER-PASSING PARTICLES
         IF (INUTYPE.EQ.0.AND.(
     &       (ISPECIES_EK(KP).EQ.0.AND.OMEGAN3RE.LT.0.).OR.
     &       (ISPECIES_EK(KP).EQ.1.AND.OMEGAN3RE.LT.0.AND.
     &        OMEGAN3RE.GT.-1.))) THEN
C        NOTE THAT THIS OPTION HOLDS ONLY FOR ENERGY-INDEPENDENT COLLISIONALITY
         IF (ISPECIES_EK(KP).EQ.0) 
     &      FNUM3 = FNUM3*(SQRT(ZEPK)-OMEGAN3)*(ZEPK+OMEGAN3RE**2)

         IF (ISPECIES_EK(KP).EQ.0) ZEPKS = 2.*ATAN(OMEGAN3RE**2)/PI
         IF (ISPECIES_EK(KP).EQ.1) ZEPKS = OMEGAN3RE**2
         IF (ISPECIES_EK(KP).EQ.0) FNUM3(NEPK2) = 0.

         IF (ZEPKS.LE.ZEPKO(1)) THEN
            FNUM30 = FNUM3(1)*ZEPKS/ZEPKO(1)
         ELSE  
            CALL SPLINE1D(FNUM0R,ZEPKS,1,DREAL(FNUM3),ZEPKO,NEPK2,FTMP3)
            CALL SPLINE1D(FNUM0I,ZEPKS,1,DIMAG(FNUM3),ZEPKO,NEPK2,FTMP3)
            FNUM30 = CMPLX(FNUM0R,FNUM0I)
         ENDIF

         IF (ISPECIES_EK(KP).EQ.0) FNUM3 = (FNUM3-FNUM30)/
     &      (ZEPK-OMEGAN3**2)/(ZEPK+OMEGAN3RE**2)*EFAC
         IF (ISPECIES_EK(KP).EQ.1) FNUM3 = (FNUM3-FNUM30)/
     &      (SQRT(ZEPK)+OMEGAN3)

         IF (ABS(DIMAG(OMEGAN3)).GT.REPS) THEN
            IF (ISPECIES_EK(KP).EQ.0) 
     &         ZTMP0 = CLOG(-OMEGAN3RE**2/OMEGAN3**2)/
     &                 (OMEGAN3**2+OMEGAN3RE**2)
            IF (ISPECIES_EK(KP).EQ.1) 
     &         ZTMP0 = 2.-2.*OMEGAN3*CLOG((1.+OMEGAN3)/OMEGAN3)
         ELSE
            IF (ISPECIES_EK(KP).EQ.0) 
     &         ZTMP0 = CI*PI*(-RTMP)/ABS(-RTMP)/2./OMEGAN3RE**2 
            IF (ISPECIES_EK(KP).EQ.1) ZTMP0 = 2.-2.*OMEGAN3RE*(
     &         LOG(-(1+OMEGAN3RE)/OMEGAN3RE)-CI*PI*(-RTMP)/ABS(-RTMP))
         ENDIF

         ZVIF = ZVIF + ZTMP0*FNUM30

         ELSE
         
         FTMP3 = 0.
         IF (INUTYPE.EQ.1) FTMP3 = H3*(1.-ZEPK**(-1.5)*H4*
     &                             (1.+RL**2/4.))/(-RTMP)
         IF (INUTYPE.EQ.2) FTMP3 = H3*(1.-ZEPK**(-1.5)*H4)/(-RTMP)
   
         IF (ISPECIES_EK(KP).EQ.0) FNUM3 = FNUM3*EFAC
         FNUM3 = FNUM3/(SQRT(ZEPK/RFAC)+OMEGAN3+CI*FTMP3)

         ENDIF

C        SUBTRACT SINGULAR CONTRIBUTION FOR LATER PITCH ANGLE INTEGRATION
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.0) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,1)
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.1) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,2)
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.2) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,3)
         IF (SLAM0(L,KP).GT.0..AND.KDPHI.EQ.3) 
     &      ZVIF=ZVIF-LOG(ABS(ZZLAMB-SLAM0(L,KP)))*SF0(L,KP,KOPT,4)

      ENDIF

      DO K=1,NEPK-1
         HEPK = EPK(K+1)-EPK(K)
         ZVIF = ZVIF + 0.5*HEPK*(FNUM1(2*K-1)+FNUM1(2*K))
         ZVIF = ZVIF + 0.5*HEPK*(FNUM3(2*K-1)+FNUM3(2*K))
      ENDDO
  
 1001 CONTINUE
      DEALLOCATE(FNUM1,FNUM3,FTMP1,FTMP3,EFAC,RFAC)

      RETURN
      END SUBROUTINE KIA_PASS
 
C=========================================================
C MESH FOR PARTICLE ENERGY EPK
C DEFINE TWO ENERGY MESHES: 
C    ZEPKO=[0,1]-MESH 
C    ZEPKN=TAN(PI*ZEPKO/2) IS THE [0,INFTY]-MESH                                           
C FOR MODELS WITH FINITE   BIRTH ENERGY: ZEPK=ZEPKO 
C FOR MODELS WITH INFINITE BIRTH ENERGY: ZEPK=ZEPKN
C AND A MULTIPLIER EFAC=PI/2/TAN(ZEPKO*PI/2) IS APPLIED
C FOR NUMERICAL INTEGRATION 
C NOTE THAT THE ENERGY MESH IS EXTENDED TO INCLUDE THE 
C END POINT (BIRTH ENERGY FOR ISPECIES_EK=1 CASE)  
C=========================================================
      SUBROUTINE KENERGYMESH
      
      USE GLOBALM
      USE ANISOTROPICM
      IMPLICIT NONE
      INTEGER  M,K,KCHECK
      REAL*8   HEPK
      
      KCHECK =0

      EPK(1) = 0.0
      HEPK   = 1.0/DFLOAT(NEPK-1)
      DO M=2,NEPK
         EPK(M)=EPK(M-1)+HEPK
      ENDDO

      DO K=1,NEPK-1
         ZEPKO(2*K-1) = 0.5*(EPK(K)*(1.0+ZW)+EPK(K+1)*(1.0-ZW))
         ZEPKO(2*K)   = 0.5*(EPK(K)*(1.0-ZW)+EPK(K+1)*(1.0+ZW))
      ENDDO

      ZEPKO(NEPK2)     = 1.0
      ZEPKN(1:NEPK2-1) = TAN(ZEPKO(1:NEPK2-1)*0.5*ACOS(-1.))
      ZEPKN(NEPK2)     = 0.

      RETURN
      END SUBROUTINE KENERGYMESH

C====================================================================
C COMPUTE ZC1 FOR ISOTROPIC ANGLE LIMIT
C FOR IF0TYPE=1 OR 2
C====================================================================
      SUBROUTINE KJP_ISOEXTRA(KP)
      
      USE DIMENSIM
      USE KINETICM
      USE RCOMDM
      USE GLOBALM
      USE ANISOTROPICM
      
      IMPLICIT NONE

      INTEGER  KP,JS,K,KCHECK
      REAL*8   DPSIDSL

      KCHECK = 0
       
      DO K  =1,2
      DO JS =1-K+2,NR  
      IF (K.EQ.1) DPSIDSL = DPSIDS(JS)
      IF (K.EQ.2) DPSIDSL = DPSIDSM(JS)
      ZC1(JS,K,KP)=-8.*SQRT(2.)/3.*PI*ESPECIES_DENF(JS,K,KP)*DPSIDSL
     &              *ESPECIES_PRE(JS,K,2)/ALPHAA1(JS,K,KP)*AAK(JS,K,KP)
     &              /(1.+EPSLONCA(JS,K,KP)**1.5)
      ENDDO
      ENDDO
      
      ZC1(1,1,KP)    = 2.*ZC1(1,2,KP)-ZC1(2,1,KP)
      ZC1(NRP1,1,KP) = 2.*ZC1(NR,2,KP)-ZC1(NR,1,KP) 

      ZC1(:,:,KP) = ZC1(:,:,KP)*PSPECIES_AT(KP)

      RETURN
      END SUBROUTINE KJP_ISOEXTRA

C===================================================    
C     A LOOP FOR FINDING THE VALUE OF RE_PMAX     
C     RE_PMAX=P_0 AS DEFINED FOR THE RE EQUILIBRIUM MODEL-2
C     CORRESPONDING TO THE MAX. VALUE OF U(P)-FUNCTION
C     FOR GIVEN VALUES OF E/Ec=RE_CONST(5) and RE_TRAD 
C     USING BISECTION METHOD (CAN BE IMPROVED WITH ITP METHOD) 
C     APPLIES FOR ISPECIES_F0=6
C===================================================
      SUBROUTINE  GET_RE_PMAX
       
      USE ANISOTROPICM
      USE REORBITM
      IMPLICIT NONE
      
      INTEGER J,KCHECK
      REAL*8  AX1,AX2,AX3,AF1,AF2,AF3,RE_PMAX2
      REAL*8  RE_UP,RE_DUP  
           
      KCHECK =1

      RE_PMAX   =-1.0
      RE_EMAX   =-1.0

C     FIND ZERO OF dU(p)/dp FIRST
C     CORRESPONDING TO RE_PMAX
      AX1       = 1.E-1
      AF1       = RE_DUP(AX1)
      AX2       = 1.E+3
      AF2       = RE_DUP(AX2)
      AX3       = 0.5*(AX1+AX2)
      AF3       = RE_DUP(AX3)
      
      DO J=1,1000
         IF (AF1*AF3.LE.0.) THEN
            AX2 = AX3
            AF2 = AF3
         ELSEIF (AF2*AF3.LE.0.) THEN
            AX1 = AX3
            AF1 = AF3
         ELSE
            WRITE(*,*) 'AX1=',AX1,'  AX2=',AX2,'  AX3=',AX3
            WRITE(*,*) 'AF1=',AF1,'  AF2=',AF2,'  AF3=',AF3
            STOP 'GET_RE_PMAX:1'
         ENDIF
         AX3 = (AX1+AX2)*0.5
         
         AF3 = RE_DUP(AX3)
         IF (ABS((AX1-AX2)/AX3).LE.1.E-5) EXIT 
      ENDDO 

C     CHECK WHETHER U(p=AX3)>0
      AX1 = AX3
      AF1 = RE_UP(AX1)

      IF (AF1.GT.0.) THEN
      IF (KCHECK.EQ.1) WRITE(*,*) 'GET_RE_PMAX: P0=',AX1
      RE_PMAX   = AX1
      AX2       = 1.E+3
      AF2       = RE_UP(AX2)
      AX3       = 0.5*(AX1+AX2)
      AF3       = RE_UP(AX3)
      
      DO J=1,1000
         IF (AF1*AF3.LE.0.) THEN
            AX2 = AX3
            AF2 = AF3
         ELSEIF (AF2*AF3.LE.0.) THEN
            AX1 = AX3
            AF1 = AF3
         ELSE
            WRITE(*,*) 'AX1=',AX1,'  AX2=',AX2,'  AX3=',AX3
            WRITE(*,*) 'AF1=',AF1,'  AF2=',AF2,'  AF3=',AF3
            STOP 'GET_RE_PMAX:2'
         ENDIF
         AX3 = (AX1+AX2)*0.5
         
         AF3 = RE_UP(AX3)
         IF (ABS((AX1-AX2)/AX3).LE.1.E-5) EXIT 
      ENDDO 

      RE_PMAX2 = AX3
      ENDIF

C     FINAL RESULT
      IF (RE_PMAX.GT.0.) RE_EMAX=RE_PMAX**2/(RE_PMAX**2+1.0)
      
      IF (RE_PMAX.LE.0.) THEN
         STOP 'GET_RE_PMAX:NEED TO INCREASE RE_CONST(5)'
      ENDIF

      IF (KCHECK.EQ.1) THEN
         WRITE(*,*) 'GET_RE_PMAX:',RE_PMAX,RE_EMAX,RE_PMAX2
      ENDIF

      RETURN
      END SUBROUTINE GET_RE_PMAX

C===================================================    
C COMPUTE U(p) FOR RE
C===================================================       
      FUNCTION RE_UP(P)
      USE REORBITM
      IMPLICIT NONE
      
      REAL*8  P,REA,CZT1,CZT2,REGAM
      REAL*8  RE_UP

      CZT1  = (RE_CONST(3)+1.0)/RE_CONST(5)
      CZT2  = CZT1/RE_TRAD
      REGAM = SQRT(P**2+1.0)
      REA   = 2.0*P**2/REGAM/CZT1

      RE_UP = (1.0/REA-1.0/TANH(REA))*(CZT2*REGAM**2/P-RE_CONST(5))-
     &        1.0 - 1.0/P**2 

      RETURN
      END FUNCTION RE_UP  

C===================================================    
C COMPUTE dU(p)/dp FOR RE
C===================================================       
      FUNCTION RE_DUP(P)
      USE REORBITM
      IMPLICIT NONE
      
      REAL*8  P,REA,CZT1,CZT2,REGAM
      REAL*8  RE_DUP

      CZT1  = (RE_CONST(3)+1.0)/RE_CONST(5)
      CZT2  = CZT1/RE_TRAD
      REGAM = SQRT(P**2+1.0)
      REA   = 2.0*P**2/REGAM/CZT1
      
      RE_DUP  = (1.0/REA-1.0/TANH(REA))*CZT2*(1.0-1.0/P**2) +
     &          (-1.0/REA**2+1.0/(SINH(REA)**2))*2.0*P*
     &          (1.0/REGAM+1.0/REGAM**3)/CZT1*
     &          (CZT2*REGAM**2/P-RE_CONST(5)) + 2.0/P**3 

      RETURN
      END FUNCTION RE_DUP  
