C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE TOROIDAL JXB TORQUE                           LIU YQ 04.3.2011
C TORQUE DENSITY COMPUTED AT HALF-INTEGER RADIAL MESH   LIU YQ 04.3.2011
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE TORQJXB
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE RCOMDM  
      USE RESMATM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'
      INTEGER    I,J,MS,KCHECK,II,KCTORQ
      REAL*8     TORQ1,TORQ2,RCHIH,RCHI,TORQFAC
      COMPLEX*16 TMP1,TMP2,TMP
      REAL*8,DIMENSION(:),ALLOCATABLE::TORQTMPM,TORQTMPI,TORQTMP
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::J1LA,J2LA,B1LA,B2LA,
     &                                       J1LB,J2LB,B1LB,B2LB
 
      KCHECK = 0
      KCTORQ = 2 

      IF (NCASE.EQ.1.OR.NCASE.EQ.2) KCHECK = 1

      IF (.NOT.ALLOCATED(TORQUEJXB)) ALLOCATE(TORQUEJXB(NR))
      TORQUEJXB = 0.0
 
      CALL GETXY(X,Y)

C     FIRST COMPUTE FLUX SURFACE AVERAGED TORQUE DENSITY
      IF (KCTORQ.EQ.1) THEN
      DO I = NTORQ,NR
         DO MS = 1,MSMAX
            TMP1 = 0.5*(J2U(I,MS)+J2U(I+1,MS))
            TMP2 = 0.5*(B1U(I,MS)+B1U(I+1,MS))
            TMP  =-CONJG(TMP1)*TMP2 + CONJG(J1U(I,MS))*B2U(I,MS)
            TORQUEJXB(I) = TORQUEJXB(I) + DREAL(TMP) 
         ENDDO
         TORQUEJXB(I) = TORQUEJXB(I)/REAL(JACOBM(I,1))/2.
      ENDDO
      ELSEIF (KCTORQ.EQ.2) THEN
      DO I = NTORQ,NR
         DO MS = 1,MSMAX
            TMP1 = CONJG(J2U(I,MS))*B1U(I,MS)
            TMP2 = CONJG(J2U(I+1,MS))*B1U(I+1,MS)
            TMP  =-0.5*(TMP1+TMP2) + CONJG(J1U(I,MS))*B2U(I,MS)
            TORQUEJXB(I) = TORQUEJXB(I) + DREAL(TMP) 
         ENDDO
         TORQUEJXB(I) = TORQUEJXB(I)/REAL(JACOBM(I,1))/2.
      ENDDO
      ELSEIF (KCTORQ.EQ.3) THEN
      ALLOCATE(TORQTMP(NR),TORQTMPM(NR),TORQTMPI(NRP1))
      TORQTMPM = 0.0
      RCHIH = 2.*PI/NCHI
      DO I=1,NR
      DO J=1,NCHI
         RCHI = RCHIH*(J-1)
         TMP1 = 0.0
         TMP2 = 0.0
         DO MS=1,MSMAX
            TMP1 = TMP1 + 0.5*(B1U(I,MS)+B1U(I+1,MS))*
     &                    EXP(CI*RM(MS,2)*RCHI)
            TMP2 = TMP2 + B3U(I,MS)*EXP(CI*RM(MS,2)*RCHI)
         ENDDO
         TMP = REQM(I,J)**2/RJAM(I,J)*TMP1*CONJG(TMP2)
         TORQTMPM(I) = TORQTMPM(I) + DREAL(TMP)
      ENDDO
      TORQTMPM(I) = TORQTMPM(I)*RCHIH
      ENDDO 
      CALL SPLINE1D(TORQTMPI,CS,NRP1,TORQTMPM,CSM,NR,TORQTMP)
      DO I=1,NR
         TORQUEJXB(I)=(TORQTMPI(I+1)-TORQTMPI(I))/CSH(I)*PI/
     &                REAL(JACOBM(I,1))/PI**2/4.
      ENDDO
      DEALLOCATE(TORQTMP,TORQTMPM,TORQTMPI)
      ENDIF

      DO I=1,NR
         TORQUEJXB(I) = TORQUEJXB(I)*REAL(JACOBM(I,1))
      ENDDO

C     SMOOTHING
      DO I=1,5
         TORQUEJXB(2:NR-1) = 0.1*TORQUEJXB(1:NR-2) + 
     &                       0.8*TORQUEJXB(2:NR-1) +
     &                       0.1*TORQUEJXB(3:NR)
      ENDDO

      DO I=1,NTORQ-1
         TORQUEJXB(I) = 0.0*TORQUEJXB(NTORQ)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) TORQUEJXB(I) = 0.0
      ENDDO

C     TWO WAYS OF COMPUTING TOTAL TOROIDAL TORQUE
      IF (KEYTORQ.LE.1) THEN
      TORQ1 = 0.0
      DO I=1,NR
         TORQ1 = TORQ1 + TORQUEJXB(I)*CSH(I)
      ENDDO
      TORQ1    = TORQ1*PI*PI*4.0
      TTORQJXB = TORQ1

C     RE-NORMALISATION
      TORQFAC = R0EXP**3*B0EXP**2/4.0E-7/PI
      WRITE(*,110) TORQ1,TORQ1*TORQFAC
 110  FORMAT('COMPUTE TOTAL JXB TORQUE:',E13.6,2X,E13.6)

      IF (KCHECK.EQ.1) THEN
     
      TORQ2 = 0.0
      RCHIH = 2*PI/NCHI
      I     = 1
      II    = NR + I
      DO J=1,NCHI
         RCHI = RCHIH*(J-1)
         TMP1 = 0.0
         TMP2 = 0.0
         DO MS=1,MSMAX
            TMP1 = TMP1 + 0.5*(B1U(II,MS)+B1U(II+1,MS))*
     &                    EXP(CI*RM(MS,2)*RCHI)
            TMP2 = TMP2 + B3U(II,MS)*EXP(CI*RM(MS,2)*RCHI)
         ENDDO
         TMP = VRRM(I,J)**2/VRJAM(I,J)*TMP1*CONJG(TMP2)
         TORQ2 = TORQ2 + DREAL(TMP)
      ENDDO
      TORQ2 = TORQ2*RCHIH*PI 
      
      RCHIH = R0EXP**3*B0EXP**2/4.0E-7/PI
      WRITE(*,110) TORQ2,TORQ2*TORQFAC

      OPEN(CHOUTP,FILE='TORQUEJXB.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,120) CSM(I),TORQUEJXB(I)
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(2(E15.8,1X))
 
      ENDIF
      ENDIF

      IF (KEYTORQ.EQ.1) THEN
      OPEN(CHOUTP,FILE='DATATORQJXB.OUT')
      REWIND(CHOUTP)
      DO I=1,NRP1
         WRITE(CHOUTP,130) CS(I),REAL(JACOBI(I,1)),
     &                     (J1U(I,MS),MS=1,MSMAX),
     &                     (J2U(I,MS),MS=1,MSMAX),
     &                     (B1U(I,MS),MS=1,MSMAX),
     &                     (B2U(I,MS),MS=1,MSMAX)
      ENDDO
      CLOSE(CHOUTP)
 130  FORMAT(900(E15.8,1X))
      ENDIF

      IF (KEYTORQ.EQ.2) THEN
      ALLOCATE(TORQTMP(8*MSMAX))
      ALLOCATE(J1LA(NRP1,MSMAX),J2LA(NRP1,MSMAX),
     &         B1LA(NRP1,MSMAX),B2LA(NRP1,MSMAX), 
     &         J1LB(NRP1,MSMAX),J2LB(NRP1,MSMAX), 
     &         B1LB(NRP1,MSMAX),B2LB(NRP1,MSMAX)) 
      ALLOCATE(TORQCROSSAB(NR),TORQCROSSBA(NR))
      TORQCROSSAB=(0.,0.)
      TORQCROSSBA=(0.,0.)
      II = MSMAX*2

C     READ SOLUTION DATA A
      OPEN(CHOUTP,FILE='DATATORQJXB_A.IN')
      REWIND(CHOUTP)
      DO I=1,NRP1
         READ(CHOUTP,*) RCHI,RCHI,
     &                  (TORQTMP(MS),MS=1,8*MSMAX)
         DO MS=1,MSMAX
            J1LA(I,MS) = TORQTMP(2*MS-1+0*II) + CI*TORQTMP(2*MS+0*II)
            J2LA(I,MS) = TORQTMP(2*MS-1+1*II) + CI*TORQTMP(2*MS+1*II)
            B1LA(I,MS) = TORQTMP(2*MS-1+2*II) + CI*TORQTMP(2*MS+2*II)
            B2LA(I,MS) = TORQTMP(2*MS-1+3*II) + CI*TORQTMP(2*MS+3*II)
         ENDDO
      ENDDO
      CLOSE(CHOUTP)

C     READ SOLUTION DATA B
      OPEN(CHOUTP,FILE='DATATORQJXB_B.IN')
      REWIND(CHOUTP)
      DO I=1,NRP1
         READ(CHOUTP,*) RCHI,RCHI,
     &                  (TORQTMP(MS),MS=1,8*MSMAX)
         DO MS=1,MSMAX
            J1LB(I,MS) = TORQTMP(2*MS-1+0*II) + CI*TORQTMP(2*MS+0*II)
            J2LB(I,MS) = TORQTMP(2*MS-1+1*II) + CI*TORQTMP(2*MS+1*II)
            B1LB(I,MS) = TORQTMP(2*MS-1+2*II) + CI*TORQTMP(2*MS+2*II)
            B2LB(I,MS) = TORQTMP(2*MS-1+3*II) + CI*TORQTMP(2*MS+3*II)
         ENDDO
      ENDDO
      CLOSE(CHOUTP)

C     COMPUTE CROSS-COUPLING TORQUE ELEMENTS FROM A TO B        
      DO I = NTORQ,NR
         DO MS = 1,MSMAX
            TMP1 = CONJG(J2LB(I,MS))*B1LA(I,MS)
            TMP2 = CONJG(J2LB(I+1,MS))*B1LA(I+1,MS)
            TMP  =-0.5*(TMP1+TMP2) + CONJG(J1LB(I,MS))*B2LA(I,MS)
            TORQCROSSAB(I) = TORQCROSSAB(I) + TMP 
         ENDDO
         TORQCROSSAB(I) = TORQCROSSAB(I)/2.
      ENDDO
      
C     COMPUTE CROSS-COUPLING TORQUE ELEMENTS FROM B TO A        
      DO I = NTORQ,NR
         DO MS = 1,MSMAX
            TMP1 = CONJG(J2LA(I,MS))*B1LB(I,MS)
            TMP2 = CONJG(J2LA(I+1,MS))*B1LB(I+1,MS)
            TMP  =-0.5*(TMP1+TMP2) + CONJG(J1LA(I,MS))*B2LB(I,MS)
            TORQCROSSBA(I) = TORQCROSSBA(I) + TMP 
         ENDDO
         TORQCROSSBA(I) = TORQCROSSBA(I)/2.
      ENDDO
      
C     SMOOTHING
      DO I=1,5
         TORQCROSSAB(2:NR-1) = 0.1*TORQCROSSAB(1:NR-2) + 
     &                         0.8*TORQCROSSAB(2:NR-1) +
     &                         0.1*TORQCROSSAB(3:NR)
         TORQCROSSBA(2:NR-1) = 0.1*TORQCROSSBA(1:NR-2) + 
     &                         0.8*TORQCROSSBA(2:NR-1) +
     &                         0.1*TORQCROSSBA(3:NR)
      ENDDO

C     PATCHING
      DO I=1,NTORQ-1
         TORQCROSSAB(I) = 0.0*TORQCROSSAB(NTORQ)
         TORQCROSSBA(I) = 0.0*TORQCROSSBA(NTORQ)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) TORQCROSSAB(I) = 0.0
         IF (CSM(I).GT.CTEDGE) TORQCROSSBA(I) = 0.0
      ENDDO

      TORQFAC = R0EXP**3*B0EXP**2/4.0E-7/PI
      TMP = 0.0
      DO I=1,NR
         TMP = TMP + TORQCROSSAB(I)*CSH(I)
      ENDDO
      TMP    = TMP*PI*PI*4.0
      WRITE(*,115) TMP,TMP*TORQFAC
 115  FORMAT('NET TORQUE JXB_AB:',2(E13.6,1X,E13.6,3X))

      TMP = 0.0
      DO I=1,NR
         TMP = TMP + TORQCROSSBA(I)*CSH(I)
      ENDDO
      TMP    = TMP*PI*PI*4.0
      WRITE(*,116) TMP,TMP*TORQFAC
 116  FORMAT('NET TORQUE JXB_BA:',2(E13.6,1X,E13.6,3X))

C     SAVE CROSS-COUPLING TORQUE ELEMENTS INTO FILES
      OPEN(CHOUTP,FILE='TORQUEJXB_AB.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,140) CSM(I),TORQCROSSAB(I)
      ENDDO
      CLOSE(CHOUTP)

      OPEN(CHOUTP,FILE='TORQUEJXB_BA.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,140) CSM(I),TORQCROSSBA(I)
      ENDDO
      CLOSE(CHOUTP)
 140  FORMAT(3(E15.8,1X))
 
      DEALLOCATE(TORQTMP)
      DEALLOCATE(J1LA,J2LA,B1LA,B2LA,J1LB,J2LB,B1LB,B2LB)
      DEALLOCATE(TORQCROSSAB,TORQCROSSBA)
      ENDIF

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE TOROIDAL TORQUE FROM REYNOLDS STRESS         LIU YQ 18.12.2012
C TORQUE DENSITY COMPUTED AT HALF-INTEGER RADIAL MESH  LIU YQ 18.12.2012
C ALSO COMPUTE MHD INDUCED PARTICLE TRANSPORT TERM      LIU YQ 12.2.2016
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE TORQREY
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER NPSTRT
      INTEGER    I,J,MS,KCHECK,II,MSA,MSB,NSA,NSB,MSPL,MSMI,
     &           LXROW,LXCOL,LYROW,LYCOL
      PARAMETER  (NSA=2,NSB=1)
C
      REAL*8     B_2,TORQFAC,H1,H2,H12,H22
      REAL*8,DIMENSION(:),ALLOCATABLE::DTMP
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW1,RW2,RW3
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::RVPHIV1,RVPHIV2,RV2V2,
     &                                       RVPHIV1M,RVPHIV2M,RV2V2M
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::RVVPHI,RVVPHIM,RVV1M,RVV2M
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::RLAPHI,RLAPHIM,RLA1M,RLA2M,
     &                                       RLBPHI,RLBPHIM,RLB1M,RLB2M
      COMPLEX*16 TMP

      INCLUDE 'integc.inc'

      KCHECK = 0
      IF (NCASE.EQ.1.OR.NCASE.EQ.2) KCHECK = 1

C     STEP 1: COMPUTE NECCESARY EQUILIBRIUM QUANTITIES

      ALLOCATE(RW1(NRP1,NCHI), RW2(NRP1,NCHI), RW3(NRP1,NCHI)) 
      ALLOCATE(RVPHIV1(NRP1,MEDIM),RVPHIV1M(NRP1,MEDIM),
     &         RVPHIV2(NRP1,MEDIM),RVPHIV2M(NRP1,MEDIM),
     &         RV2V2(NRP1,MEDIM),  RV2V2M(NRP1,MEDIM))
      ALLOCATE(RVVPHI(NRP1,MSDIM), RVVPHIM(NRP1,MSDIM),
     &         RVV1M(NRP1,MSDIM),  RVV2M(NRP1,MSDIM))
      ALLOCATE(DTMP(NR))
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'TORQREY'
C
      RVPHIV1  = 0.
      RVPHIV2  = 0.
      RV2V2    = 0.
      RVPHIV1M = 0.
      RVPHIV2M = 0.
      RV2V2M   = 0.

      RVVPHI  = 0.
      RVVPHIM = 0.
      RVV1M   = 0.
      RVV2M   = 0.

      DO J=1,NCHI
         DO I=2,NRP1
            B_2     =G22L(I,J)*DPSIDS(I)**2/RJA(I,J)**2 + 
     &               T(I)**2/REQ(I,J)**2
            RW1(I,J)=-DPSIDS(I)*T(I)*G12L(I,J)/RJA(I,J)/B_2
            RW2(I,J)=-DPSIDS(I)*REQ(I,J)**2*G22L(I,J)/RJA(I,J)/B_2
            RW3(I,J)=RJA(I,J)*REQ(I,J)**2/T(I)
         ENDDO
      ENDDO

      DO J=1,NCHI
         RW1(1,J) = RW1(2,J)
         RW2(1,J) = RW2(2,J)
         RW3(1,J) = RW3(2,J)
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER( RW1,  RVPHIV1, FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RVPHIV1  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW2,  RVPHIV2, FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RVPHIV2  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW3,  RV2V2,   FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RV2V2  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  RVPHIV1, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUFFTP, 'RVPHIV1')
        call FFTOUTPT(RW2,  RVPHIV2, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUFFTP, 'RVPHIV2')
        call FFTOUTPT(RW3,  RV2V2,   NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUFFTP, 'RV2V2')
      ENDIF
C
      DO J=1,NCHI
         DO I=1,NR
            B_2     =G22LM(I,J)*DPSIDSM(I)**2/RJAM(I,J)**2 + 
     &               TM(I)**2/REQM(I,J)**2
            RW1(I,J)=-DPSIDSM(I)*TM(I)*G12LM(I,J)/RJAM(I,J)/B_2
            RW2(I,J)=-DPSIDSM(I)*REQM(I,J)**2*G22LM(I,J)/RJAM(I,J)/B_2
            RW3(I,J)=RJAM(I,J)*REQM(I,J)**2/TM(I)
         ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER( RW1,  RVPHIV1M,FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RVPHIV1M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW2,  RVPHIV2M,FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RVPHIV2M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      5,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW3,  RV2V2M,  FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RV2V2M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  RVPHIV1M,NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUFFTP, 'RVPHIV1M')
        call FFTOUTPT(RW2,  RVPHIV2M,NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUFFTP, 'RVPHIV2M')
        call FFTOUTPT(RW3,  RV2V2M,  NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUFFTP, 'RV2V2M')
      ENDIF

C     STEP 2: COMPUTE RVVPHI=V_PHI^*, RVV1M=J*V1, RVV2M=J*R^2/T*V2
C     NOTE THAT V_PHI^* IS COMPUTED AT BOTH INTERGER & HALF-INTEGER 
C     POINTS. V1 & V2 ARE COMPUTED ONLY AT HALF-INTEGER POINTS

      DO 90 MSA=1,MSMAX
      DO 90 MSB=1,MSMAX
      MSPL =  MPLUS(MSA,NSA,MSB,NSB)
      MSMI = MMINUS(MSA,NSA,MSB,NSB)
      IF (MSPL.LT.1) GOTO 60

      LXROW = MSPL
      LXCOL = MSA 
      LYROW = MSPL
      LYCOL = MSA 

      DO 50 I=2,NRP1
      INCLUDE 'tent.inc'

      RVVPHI(I,LXROW) = RVVPHI(I,LXROW)  
     &+FF(RVPHIV1(i,msb),RVPHIV1M(i-1,msb),RVPHIV1M(i,msb))*V1U(I,LXCOL)
     &+FFM(RVPHIV1M(i-1,msb))*V1U(I-1,LXCOL)
     &+FFP(RVPHIV1M(i,msb))*V1U(I+1,LXCOL)
     &+FGM(RVPHIV2(i,msb),RVPHIV2M(i-1,msb))*V2U(I-1,LYCOL)
     &+FGP(RVPHIV2(i,msb),RVPHIV2M(i,msb))*V2U(I,LYCOL)
 50   CONTINUE
 
      DO 55 I=1,NR
      INCLUDE 'tophat.inc'
 
      RVVPHIM(I,LYROW) = RVVPHIM(I,LYROW) 
     &+GG(RVPHIV2M(i,msb),RVPHIV2(i,msb),RVPHIV2(i+1,msb))*V2U(I,LYCOL)
     &+GF(RVPHIV1(i,msb),RVPHIV1M(i,msb))*V1U(I,LXCOL)
     &+GF(RVPHIV1(i+1,msb),RVPHIV1M(i,msb))*V1U(I+1,LXCOL)
      RVV1M(I,LYROW) = RVV1M(I,LYROW) 
     &+GF(JACOBI(i,msb),JACOBM(i,msb))*V1U(I,LXCOL)
     &+GF(JACOBI(i+1,msb),JACOBM(i,msb))*V1U(I+1,LXCOL)
      RVV2M(I,LYROW) = RVV2M(I,LYROW) 
     &+GG(RV2V2M(i,msb),RV2V2(i,msb),RV2V2(i+1,msb))*V2U(I,LYCOL)
 55   CONTINUE
 60   CONTINUE
      IF (MSB.LT.2) GOTO 80
      IF (MSMI.LT.1) GOTO 80
 
      LXROW = MSMI
      LXCOL = MSA 
      LYROW = MSMI
      LYCOL = MSA 

      DO 70 I=2,NRP1
      INCLUDE 'tent.inc'
 
      RVVPHI(I,LXROW) = RVVPHI(I,LXROW)
     &+CONJG(FF(RVPHIV1(i,msb),RVPHIV1M(i-1,msb),RVPHIV1M(i,msb)))
     & *V1U(I,LXCOL)
     &+CONJG(FFM(RVPHIV1M(i-1,msb)))*V1U(I-1,LXCOL)
     &+CONJG(FFP(RVPHIV1M(i,msb)))*V1U(I+1,LXCOL)
     &+CONJG(FGM(RVPHIV2(i,msb),RVPHIV2M(i-1,msb)))*V2U(I-1,LYCOL)
     &+CONJG(FGP(RVPHIV2(i,msb),RVPHIV2M(i,msb)))*V2U(I,LYCOL)
 70   CONTINUE
 
      DO 75 I=1,NR
      INCLUDE 'tophat.inc'
 
      RVVPHIM(I,LYROW) = RVVPHIM(I,LYROW)
     &+CONJG(GG(RVPHIV2M(i,msb),RVPHIV2(i,msb),RVPHIV2(i+1,msb)))
     & *V2U(I,LYCOL)
     &+CONJG(GF(RVPHIV1(i,msb),RVPHIV1M(i,msb)))*V1U(I,LXCOL)
     &+CONJG(GF(RVPHIV1(i+1,msb),RVPHIV1M(i,msb)))*V1U(I+1,LXCOL)
      RVV1M(I,LYROW) = RVV1M(I,LYROW)
     &+CONJG(GF(JACOBI(i,msb),JACOBM(i,msb)))*V1U(I,LXCOL)
     &+CONJG(GF(JACOBI(i+1,msb),JACOBM(i,msb)))*V1U(I+1,LXCOL)
      RVV2M(I,LYROW) = RVV2M(I,LYROW)
     &+CONJG(GG(RV2V2M(i,msb),RV2V2(i,msb),RV2V2(i+1,msb)))
     & *V2U(I,LYCOL)
 75   CONTINUE
 80   CONTINUE
 90   CONTINUE

C     TERMS WITHOUT CONVOLUTION
      DO MS=1,MSMAX
      LXROW = MS
      DO I=2,NRP1
      INCLUDE 'tent.inc'
      RVVPHI(I,LXROW) = RVVPHI(I,LXROW)  
     &+FGM(T(i)*C1,TM(i-1)*C1)*V3U(I-1,LXROW)
     &+FGP(T(i)*C1,TM(i)*C1)*V3U(I,LXROW)
      ENDDO

      DO I=1,NR
      INCLUDE 'tophat.inc'
      RVVPHIM(I,LXROW) =  RVVPHIM(I,LXROW)
     &+GG(TM(i)*C1,T(i)*C1,T(i+1)*C1)*V3U(I,LXROW)
      ENDDO   
      ENDDO

      DO MS=1,MSMAX
      DO I=1,NRP1
         RVVPHI(I,MS)  = CONJG(RVVPHI(I,MS))
         RVVPHIM(I,MS) = CONJG(RVVPHIM(I,MS))
      ENDDO
      ENDDO

C     STEP 3: COMPUTE THE TORQUE DENSITY AT HALF-INTEGER MESH
C             ALSO COMPUTE MHD INDUCED PARTICLE TRANSPORT TERM

      IF (.NOT.ALLOCATED(TORQUEREY)) ALLOCATE(TORQUEREY(NR))
      IF (.NOT.ALLOCATED(DNTRANMHD)) ALLOCATE(DNTRANMHD(NR),
     &                                        DPTRANMHD(NR))
      TORQUEREY = 0.0
      DNTRANMHD = 0.0
      DPTRANMHD = 0.0
      DTMP      = 0.0

      DO MS=1,MSMAX
      DO I=1,NR
         TORQUEREY(I) = TORQUEREY(I) + REAL(
     &   RVV1M(I,MS)*(RVVPHI(I+1,MS)-RVVPHI(I,MS))/CSH(I) - 
     &   RVV2M(I,MS)*RVVPHIM(I,MS)*CI*RM(MS,2) )
         DPTRANMHD(I) = DPTRANMHD(I) + 
     &      0.5*REAL(RVV1M(I,MS)*CONJG(RHOP(I,MS)))/REAL(JACOBM(I,1))
      ENDDO
      ENDDO

      DO I=1,NR
         TORQUEREY(I) = -TORQUEREY(I)*0.5*RHOM(I)/REAL(JACOBM(I,1))
      ENDDO

      DO I=1,NR
         TORQUEREY(I) = TORQUEREY(I)*REAL(JACOBM(I,1))
         DTMP(I)      =-DPTRANMHD(I)*REAL(JACOBM(I,1))
      ENDDO

      DO I=2,NR-1
         H1  = CSM(I) - CSM(I-1)
         H2  = CSM(I+1) - CSM(I)
         H12 = H1*H1
         H22 = H2*H2
         DNTRANMHD(I) = (H12*DTMP(I+1)-H22*DTMP(I-1)+(H22-H12)*DTMP(I))/
     &                  (H12*H2+H22*H1)/REAL(JACOBM(I,1))
      ENDDO
      DNTRANMHD(1)  = DNTRANMHD(2)
      DNTRANMHD(NR) = DNTRANMHD(NR-1)

C     SMOOTHING
      DO I=1,5
         TORQUEREY(2:NR-1) = 0.1*TORQUEREY(1:NR-2) + 
     &                       0.8*TORQUEREY(2:NR-1) +
     &                       0.1*TORQUEREY(3:NR)
      ENDDO

C     SMOOTHING
      DO I=1,5
         DNTRANMHD(2:NR-1) = 0.1*DNTRANMHD(1:NR-2) + 
     &                       0.8*DNTRANMHD(2:NR-1) +
     &                       0.1*DNTRANMHD(3:NR)
      ENDDO

      DO I=1,NTORQ-1
         TORQUEREY(I) = 0.0*TORQUEREY(NTORQ)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) TORQUEREY(I) = 0.0
      ENDDO

      DO I=1,NDNTR-1
         DNTRANMHD(I) = DNTRANMHD(NDNTR)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CDEDGE) DNTRANMHD(I) = 0.0
         IF (CSM(I).GT.CDEDGE) DPTRANMHD(I) = 0.0
      ENDDO

      IF (KEYTORQ.LE.1) THEN
C     NET TORQUE FROM REYNOLDS STRESS
      TTORQREY = 0.0
      DO I=1,NR
         TTORQREY = TTORQREY + TORQUEREY(I)*CSH(I)
      ENDDO
      TTORQREY    = TTORQREY*PI*PI*4.0

C     NET DENSITY RADIAL TRANSPORT DUE TO MHD TERM
      TDNTRMHD = 0.0
      DO I=1,NR
         TDNTRMHD = TDNTRMHD + DNTRANMHD(I)*CSH(I)
      ENDDO
      TDNTRMHD    = TDNTRMHD*PI*PI*4.0

C     RE-NORMALISATION
      TORQFAC = R0EXP**3*B0EXP**2/4.0E-7/PI
      WRITE(*,110) TTORQREY,TTORQREY*TORQFAC
 110  FORMAT('COMPUTE TOTAL REYNOLDS TORQUE:',E13.6,2X,E13.6)
      WRITE(*,112) TDNTRMHD,TDNTRMHD/TOTDENS
 112  FORMAT('COMPUTE TOTAL MHD DENSITY PUMPOUT:',E13.6,2X,E13.6)

      IF (KCHECK.EQ.1) THEN
      OPEN(CHOUTP,FILE='TORQUEREY.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,120) CSM(I),TORQUEREY(I),DPTRANMHD(I)
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(3(E15.8,1X))
      ENDIF
      ENDIF
      
      IF (KEYTORQ.EQ.1) THEN
      OPEN(CHOUTP,FILE='DATATORQREY.OUT')
      REWIND(CHOUTP)
      DO I=1,NRP1
         WRITE(CHOUTP,130) CS(I),REAL(JACOBI(I,1)),
     &                     (RVV1M(I,MS),MS=1,MSMAX),
     &                     (RVV2M(I,MS),MS=1,MSMAX),
     &                     (RVVPHI(I,MS),MS=1,MSMAX),
     &                     (RVVPHIM(I,MS),MS=1,MSMAX)
      ENDDO
      CLOSE(CHOUTP)
 130  FORMAT(900(E15.8,1X))
      ENDIF

C     DEALLOCATE LOCAL VARIABLES
      DEALLOCATE(RW1,RW2,RW3,DTMP)
      DEALLOCATE(RVPHIV1,RVPHIV2,RV2V2,RVPHIV1M,RVPHIV2M,RV2V2M)
      DEALLOCATE(RVVPHI,RVVPHIM,RVV1M,RVV2M)

      IF (KEYTORQ.EQ.2) THEN
      ALLOCATE(DTMP(8*MSMAX))
      ALLOCATE(RLA1M(NRP1,MSMAX), RLA2M(NRP1,MSMAX),
     &         RLAPHI(NRP1,MSMAX),RLAPHIM(NRP1,MSMAX), 
     &         RLB1M(NRP1,MSMAX), RLB2M(NRP1,MSMAX), 
     &         RLBPHI(NRP1,MSMAX),RLBPHIM(NRP1,MSMAX)) 
      ALLOCATE(TORQCROSSAB(NR),TORQCROSSBA(NR))
      TORQCROSSAB=(0.,0.)
      TORQCROSSBA=(0.,0.)
      II = MSMAX*2

C     READ SOLUTION DATA A
      OPEN(CHOUTP,FILE='DATATORQREY_A.IN')
      REWIND(CHOUTP)
      DO I=1,NRP1
         READ(CHOUTP,*) H22,H22,
     &                  (DTMP(MS),MS=1,8*MSMAX)
         DO MS=1,MSMAX
            RLA1M(I,MS)   = DTMP(2*MS-1+0*II) + CI*DTMP(2*MS+0*II)
            RLA2M(I,MS)   = DTMP(2*MS-1+1*II) + CI*DTMP(2*MS+1*II)
            RLAPHI(I,MS)  = DTMP(2*MS-1+2*II) + CI*DTMP(2*MS+2*II)
            RLAPHIM(I,MS) = DTMP(2*MS-1+3*II) + CI*DTMP(2*MS+3*II)
         ENDDO
      ENDDO
      CLOSE(CHOUTP)

C     READ SOLUTION DATA B
      OPEN(CHOUTP,FILE='DATATORQREY_B.IN')
      REWIND(CHOUTP)
      DO I=1,NRP1
         READ(CHOUTP,*) H22,H22,
     &                  (DTMP(MS),MS=1,8*MSMAX)
         DO MS=1,MSMAX
            RLB1M(I,MS)   = DTMP(2*MS-1+0*II) + CI*DTMP(2*MS+0*II)
            RLB2M(I,MS)   = DTMP(2*MS-1+1*II) + CI*DTMP(2*MS+1*II)
            RLBPHI(I,MS)  = DTMP(2*MS-1+2*II) + CI*DTMP(2*MS+2*II)
            RLBPHIM(I,MS) = DTMP(2*MS-1+3*II) + CI*DTMP(2*MS+3*II)
         ENDDO
      ENDDO
      CLOSE(CHOUTP)

C     COMPUTE CROSS-COUPLING TORQUE ELEMENTS BETWEEN A AND B    
      DO I = NTORQ,NR
      DO MS = 1,MSMAX
         TORQCROSSAB(I) = TORQCROSSAB(I) + 
     &      RLA1M(I,MS)*(RLBPHI(I+1,MS)-RLBPHI(I,MS))/CSH(I) - 
     &      RLA2M(I,MS)*RLBPHIM(I,MS)*CI*RM(MS,2) 
         TORQCROSSBA(I) = TORQCROSSBA(I) + 
     &      RLB1M(I,MS)*(RLAPHI(I+1,MS)-RLAPHI(I,MS))/CSH(I) - 
     &      RLB2M(I,MS)*RLAPHIM(I,MS)*CI*RM(MS,2) 
      ENDDO
      TORQCROSSAB(I) =-TORQCROSSAB(I)*0.5*RHOM(I)
      TORQCROSSBA(I) =-TORQCROSSBA(I)*0.5*RHOM(I)
      ENDDO
      
C     SMOOTHING
      DO I=1,5
         TORQCROSSAB(2:NR-1) = 0.1*TORQCROSSAB(1:NR-2) + 
     &                         0.8*TORQCROSSAB(2:NR-1) +
     &                         0.1*TORQCROSSAB(3:NR)
         TORQCROSSBA(2:NR-1) = 0.1*TORQCROSSBA(1:NR-2) + 
     &                         0.8*TORQCROSSBA(2:NR-1) +
     &                         0.1*TORQCROSSBA(3:NR)
      ENDDO

C     PATCHING
      DO I=1,NTORQ-1
         TORQCROSSAB(I) = 0.0*TORQCROSSAB(NTORQ)
         TORQCROSSBA(I) = 0.0*TORQCROSSBA(NTORQ)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) TORQCROSSAB(I) = 0.0
         IF (CSM(I).GT.CTEDGE) TORQCROSSBA(I) = 0.0
      ENDDO

      TORQFAC = R0EXP**3*B0EXP**2/4.0E-7/PI
      TMP = 0.0
      DO I=1,NR
         TMP = TMP + TORQCROSSAB(I)*CSH(I)
      ENDDO
      TMP    = TMP*PI*PI*4.0
      WRITE(*,115) TMP,TMP*TORQFAC
 115  FORMAT('NET TORQUE REY_AB:',2(E13.6,1X,E13.6,3X))

      TMP = 0.0
      DO I=1,NR
         TMP = TMP + TORQCROSSBA(I)*CSH(I)
      ENDDO
      TMP    = TMP*PI*PI*4.0
      WRITE(*,116) TMP,TMP*TORQFAC
 116  FORMAT('NET TORQUE REY_BA:',2(E13.6,1X,E13.6,3X))

C     SAVE CROSS-COUPLING TORQUE ELEMENTS INTO FILES
      OPEN(CHOUTP,FILE='TORQUEREY_AB.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,140) CSM(I),TORQCROSSAB(I)
      ENDDO
      CLOSE(CHOUTP)

      OPEN(CHOUTP,FILE='TORQUEREY_BA.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,140) CSM(I),TORQCROSSBA(I)
      ENDDO
      CLOSE(CHOUTP)
 140  FORMAT(3(E15.8,1X))
 
      DEALLOCATE(DTMP)
      DEALLOCATE(RLA1M,RLA2M,RLAPHI,RLAPHIM,RLB1M,RLB2M,RLBPHI,RLBPHIM)
      DEALLOCATE(TORQCROSSAB,TORQCROSSBA)
      ENDIF

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE TOROIDAL TORQUE FROM FIELD LINE ERGODIZATION  TOMASINA E. 30/10/2024
C TORQUE DENSITY COMPUTED AT HALF-INTEGER RADIAL MESH  LIU YQ 18.12.2012
C ALSO COMPUTE MHD INDUCED PARTICLE TRANSPORT TERM      LIU YQ 12.2.2016
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         SUBROUTINE TORQERGO
C     ==========================================================
         USE DIMENSIM
         USE GLOBALM
         USE RCOMDM  
         USE TORQUEM
         IMPLICIT NONE
         INCLUDE 'comioc.inc'
         INCLUDE 'compam.inc'
         INCLUDE 'comfft.inc'

         INTEGER I, J, M, MS, KCHECK
         REAL*8 TMPR,H1,H2,H12,H22, TORQFAC, ALPHAST, ZQMI,ZQME
         REAL*8,DIMENSION(:),ALLOCATABLE:: DLNTE, DLNTI, DLNR
         REAL*8,DIMENSION(:),ALLOCATABLE:: BABS2, DST, BPOLM 
         REAL*8,DIMENSION(:),ALLOCATABLE:: BPVPHIM
         REAL*8,DIMENSION(:),ALLOCATABLE:: MUI, SIGMAST, SIGMANEO
         REAL*8,DIMENSION(:),ALLOCATABLE:: E1ST, E1NEO, JIST
         COMPLEX*16 TMP
   
         INCLUDE 'integc.inc'
   
         KCHECK = 0
         IF (NCASE.EQ.1.OR.NCASE.EQ.2) KCHECK = 1
   
         IF (.NOT.ALLOCATED(TORQUEERGO)) ALLOCATE(TORQUEERGO(NR))
         IF (.NOT.ALLOCATED(DNTRANERGO)) ALLOCATE(DNTRANERGO(NR),
     &                                        DPTRANERGO(NR))
         ALLOCATE(BPOLM(NR), BPVPHIM(NR))
         ALLOCATE(DLNTE(NRP1), DLNTI(NRP1), DLNR(NRP1))
         ALLOCATE(BABS2(NR), DST(NR))
         ALLOCATE(MUI(NR), SIGMAST(NR), SIGMANEO(NR))
         ALLOCATE(E1ST(NR), E1NEO(NR), JIST(NR))

         INCLUDE 'setfft.inc'
         SUBNAM    = 'TORQERGO'
C
         TORQUEERGO = 0.
         DNTRANERGO = 0.
         DPTRANERGO = 0.

         IF (CTERGO.GT.0.0.OR.CDERGO.GT.0.0) THEN
         BPOLM = 0.
         BABS2 = 0.

         ZQMI = 1.67261e-27 !PROTON MASS    [KG]
         ZQME = 9.109E-31   !ELECTRON MASS    [KG]

C     COMPUTE DERIVATIVES OF KINETIC QUANTITIES
         !ELECTRON TEMPERATURE
         FFF(:,1) = TEMPE(:)
         FFF(:,2) = TEMPEM(:)
         CALL DFFFDPSI(1)
         DLNTE = DFFF(:,2)*DPSIDSM(1:NRP1)
         !ION TEMPERATURE
         FFF(:,1) = TEMPI(:)
         FFF(:,2) = TEMPIM(:)
         CALL DFFFDPSI(1)
         DLNTI = DFFF(:,2)*DPSIDSM(1:NRP1)
         !DENSITY
         FFF(:,1) = RHO(:)
         FFF(:,2) = RHOM(:)
         CALL DFFFDPSI(1)
         DLNR = DFFF(:,2)*DPSIDSM(1:NRP1)

C     COMPUTE VPHI*BPOL 
         DO I = 1 , NR
            DO J = 1, NCHI 
               TMPR  = ABS(DPSIDS(I)) * SQRT(G22L(I,J)) / RJa(I,J) + 
     &                ABS(DPSIDS(I+1)) * SQRT(G22L(I+1,J)) / RJa(I+1,J)    
               TMPR = 0.5*TMPR          
               BPOLM(I) = BPOLM(I) +  TMPR!TMPR/NCHI
            END DO
            BPOLM(I) = BPOLM(I)/NCHI !SURFACE AVERAGE
            BPVPHIM(I) = BPOLM(I) * ROTM(I)
         END DO
            
C     DETERMINE COLLISIONALITY REGIME
      ALPHAST = (1-ALPHACOL) * 0.5 + ALPHACOL * 1.71

C     COMPUTE TORQUE DENSITY AT HALF-INTEGER MESH AND ERGODIZATION DRIVEN PARTICLE TRANSPORT
         DO I = 1, NR
            !COMPUTE DIFFUSION COEFFICENT DUE TO ERGODIZATION (QL APPROXIMATION)
            DO M = 1, MSMAX
               TMP = 0.5*(B1U(I,M)+B1U(I+1,M))
               TMP = ABS(TMP)**2
               BABS2(I) = BABS2(I) + TMP
            END DO
            !BABS2(I) = BABS2(I) / REAL(JACOBM(I,1))/2. 
            DST(I) = 2.*PI*QM(I) * BABS2(I)
         
            !COMPUTE CONDUCTIVITY COEFFICIENTS      
            SIGMAST(I) = OMEGACI0**2*RHOM(I)* DST(I) *
     &                   SQRT(2.*ZQMI/(ZQME*TEMPEM(I)*PI))
 
   ! SELF CONSISTENT ION CONDUCTIVITY -> REQUIRE TO FIX NORMALIZATION
   !          MUI(I) = SQRT(ZQME/ZQMI)*TEMPIM(I)/OMEGACI0**2*(
   !   &               TEMPIM(I)/TEMPEM(I))**(3./2.)/RESISM(I)
  
            MUI(I) = 3./4.*TEMPIM(I)**(5./2.)/OMEGACI0**4/(SQRT(PI)*10.)  
            SIGMANEO(I) = 3./2. * MUI(I)

            !COMPUTE ION RADIAL CURRENT
            E1ST(I) = -TEMPEM(I)/OMEGACI0 * (DLNR(I) + ALPHAST*DLNTE(I))
            E1NEO(I) = TEMPIM(I)/OMEGACI0 * 
     &                 (DLNR(I) + DLNTI(I)) + BPVPHIM(I)  
            JIST(I) =   SIGMAST(I)*SIGMANEO(I)/(SIGMAST(I)+
     &                  SIGMANEO(I)) * (E1ST(I)-E1NEO(I)) 

            JIST(I) = JIST(I)*REAL(JACOBM(I,1))
            !COMPUTE TORQUE DENSITY     
            TORQUEERGO(I) = JIST(I) * BPOLM(I) !IGNORING CONTRIBUTION TO POLOIDAL FLOW
            !COMPUTE STOCHASTIC PARTICLE PUMP-OUT 
            DPTRANERGO(I) = JIST(I)/OMEGACI0 

         END DO

C     HANDLE NaN
      DO I = 1, NR
         IF (TORQUEERGO(I).NE.TORQUEERGO(I)) TORQUEERGO(I) = 0.
         IF (DPTRANERGO(I).NE.DPTRANERGO(I)) DPTRANERGO(I) = 0.
      END DO

C     COMPUTE STOCHASTIC DENSITY PUMP-OUT 
         DO I=2,NR-1
            H1  = CSM(I) - CSM(I-1)
            H2  = CSM(I+1) - CSM(I)
            H12 = H1*H1
            H22 = H2*H2
            DNTRANERGO(I) = (H12*DPTRANERGO(I+1)-H22*
     &                      DPTRANERGO(I-1)+(H22-H12)*
     &                      DPTRANERGO(I))/(H12*H2+H22*H1)
            DNTRANERGO(I) = DNTRANERGO(I)/REAL(JACOBM(I,1))                    
         ENDDO
         DNTRANERGO(1)  = DNTRANERGO(2)
         DNTRANERGO(NR) = DNTRANERGO(NR-1)

C     SMOOTHING
         DO I=1,5
            TORQUEERGO(2:NR-1) = 0.1*TORQUEERGO(1:NR-2) + 
     &                         0.8*TORQUEERGO(2:NR-1) +
     &                         0.1*TORQUEERGO(3:NR)
         ENDDO

C     SMOOTHING
         DO I=1,5
            DNTRANERGO(2:NR-1) = 0.1*DNTRANERGO(1:NR-2) + 
     &                           0.8*DNTRANERGO(2:NR-1) +
     &                           0.1*DNTRANERGO(3:NR)
         ENDDO

C     PATCH VALUES AT AXIS
         DO I=1,NTORQ-1
            TORQUEERGO(I) = TORQUEERGO(NTORQ)
         ENDDO

         IF (KEYTORQ.LE.1) THEN
C     NET TORQUE FROM EDGE ERGODIZATION
         TTORQERGO = 0.0
         DO I=1,NR
            TTORQERGO = TTORQERGO + TORQUEERGO(I)*CSH(I)
         ENDDO
         TTORQERGO = TTORQERGO*PI*PI*4.0
   
C     NET DENSITY RADIAL TRANSPORT DUE TO MHD TERM
         TDNTRERGO = 0.0
         DO I=1,NR
            TDNTRERGO = TDNTRERGO + DNTRANERGO(I)*CSH(I)
         ENDDO
         TDNTRERGO = TDNTRERGO*PI*PI*4.0
   
C     RE-NORMALISATION
         TORQFAC = R0EXP**3*B0EXP**2/4.0E-7/PI
         WRITE(*,110) TTORQERGO,TTORQERGO*TORQFAC
  110  FORMAT('COMPUTE TOTAL STOCHASTIC TORQUE:',E13.6,2X,E13.6)
         WRITE(*,112) TDNTRERGO,TDNTRERGO/TOTDENS
  112  FORMAT('COMPUTE TOTAL STOCHASTIC DENSITY PUMPOUT:'
     &               ,E13.6,2X,E13.6)
   
         IF (KCHECK.EQ.1) THEN
         OPEN(CHOUTP,FILE='TORQUEERGO.OUT')
         REWIND(CHOUTP)
         DO I=1,NR
            WRITE(CHOUTP,120) CSM(I),TORQUEERGO(I),DPTRANERGO(I)
         ENDDO
         CLOSE(CHOUTP)
  120  FORMAT(3(E15.8,1X))
         ENDIF
         ENDIF
         
         ENDIF  !IF (CT.ERGO.GT.0.0.OR.CDERGO.GT.0.0)

C     <TORQUE MATRIX TO BE IMPLEMENTED>
C     DEALLOCATE QUANTITIES
         DEALLOCATE(BPOLM, BPVPHIM)
         DEALLOCATE(DLNTE, DLNTI, DLNR)
         DEALLOCATE(BABS2, DST)
         DEALLOCATE(SIGMAST, SIGMANEO)
         DEALLOCATE(E1ST, E1NEO, JIST)
                  
         RETURN
         END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE TOROIDAL NTV TORQUE                           LIU YQ 03.5.2011
C TORQUE DENSITY COMPUTED AT HALF-INTEGER RADIAL MESH   LIU YQ 03.5.2011
C ALSO COMPUTE NTV INDUCED PARTICLE TRANSPORT TERM      LIU YQ 12.2.2016
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE TORQNTV
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INTEGER KCHECK,I,J,K
      REAL*8  TORQ1,TORQFAC,H1,H2,H12,H22,TORQI,TORQE
      REAL*8,DIMENSION(:),ALLOCATABLE::DVDSM,DTMP
      REAL*8,DIMENSION(:),ALLOCATABLE::DPTRANNTVI,DPTRANNTVE
      COMPLEX*16 TMP

      KCHECK = 0
      IF (NCASE.EQ.1.OR.NCASE.EQ.2) KCHECK = 1

      IF (.NOT.ALLOCATED(TORQUENTV)) ALLOCATE(TORQUENTV(NR),
     &                                        TORQUENTVI(NR),
     &                                        TORQUENTVE(NR),
     &                                        ROTOFFSET(NR))
      IF (.NOT.ALLOCATED(DNTRANNTV)) ALLOCATE(DNTRANNTV(NR),
     &                                        DPTRANNTV(NR))
      IF (.NOT.ALLOCATED(DPTRANNTVI)) ALLOCATE(DPTRANNTVI(NR),
     &                                         DPTRANNTVE(NR))
      TORQUENTV = 0.
      ROTOFFSET = 0.
      DNTRANNTV = 0.
      DPTRANNTV = 0.
      TORQUENTVI = 0.
      TORQUENTVE = 0.
      DPTRANNTVI = 0.
      DPTRANNTVE = 0.

      IF (KEYTORQ.EQ.2) THEN
         ALLOCATE(TORQCROSSAB(NR),TORQCROSSBA(NR))
         TORQCROSSAB=(0.,0.)
         TORQCROSSBA=(0.,0.)
      ENDIF

C     COMPUTE DV/DS AT HALF-INTEGER POINTS
C     REMOVED VOLTOT FACTOR, 07/11/2016, Y.Q.LIU
      ALLOCATE(DVDSM(NR),DTMP(NR))
      DO I=1,NR
         DVDSM(I) = CSVM(I)*(CSV(I+1)-CSV(I))/CSH(I)
      ENDDO
C     DVDSM = DVDSM*VOLTOT*2.0
      DVDSM = DVDSM*2.0
      
      IF (KNTV.EQ.10) 
     &   CALL TORQNTV1(ZCHARGE,NR,NCHI,INT(-RNTOR),M2,KEYTORQ,
     &                 CS,ASPCT,CTEDGE,DVDSM,KCHECK)
      IF (KNTV.EQ.11) THEN 
         CALL TORQNTV1(1,NR,NCHI,INT(-RNTOR),M2,KEYTORQ,
     &                 CS,ASPCT,CTEDGE,DVDSM,KCHECK)
         TORQUENTVI = TORQUENTV
         CALL TORQNTV1(-1,NR,NCHI,INT(-RNTOR),M2,KEYTORQ,
     &                 CS,ASPCT,CTEDGE,DVDSM,KCHECK)
         TORQUENTVE = TORQUENTV
         TORQUENTV  = TORQUENTVI + TORQUENTVE
      ENDIF

      IF ((KNTV.EQ.20.OR.KNTV.EQ.21).AND.INCKIN.GT.0) CALL KDWKDENSITY
      IF (KNTV.EQ.21.AND.INCKIN.GT.0.AND.KXX1.GT.0.AND.KYX2.GT.0) THEN
        WRITE(*,*) 'CALLING CALCDWKCOMP FOR NTV TORQUE'
        CALL CALCDWKCOMP
C        WRITE(*,*) 'CALLING KDWKDENSITY FOR NTV TORQUE'
C        CALL KDWKDENSITY
      ENDIF
C     SET TORQUENTV=0 WITHIN RESISTIVE LAYER
      IF (DELRATS.GT.0..AND.DELRATS.LT.1..AND.NRATSURF.GT.0) THEN
      DO I=1,NR
         H1 = CSM(I)
         DO K=1,NRATSURF
            J = IRATSURF(K)
            IF (ABS(H1-CS(J)).LE.DELRATS/Q(J)) THEN      
               TORQUENTV(I)  = 0.
               TORQUENTVI(I) = 0.
               TORQUENTVE(I) = 0.
            ENDIF
         ENDDO
      ENDDO
      ENDIF

C     SMOOTHING
      DO I=1,5
         TORQUENTV(2:NR-1) = 0.1*TORQUENTV(1:NR-2) + 
     &                       0.8*TORQUENTV(2:NR-1) +
     &                       0.1*TORQUENTV(3:NR)
      ENDDO
      IF (KNTV.EQ.11.OR.KNTV.EQ.21) THEN
      DO I=1,5
         TORQUENTVI(2:NR-1) = 0.1*TORQUENTVI(1:NR-2) + 
     &                       0.8*TORQUENTVI(2:NR-1) +
     &                       0.1*TORQUENTVI(3:NR)
         TORQUENTVE(2:NR-1) = 0.1*TORQUENTVE(1:NR-2) + 
     &                       0.8*TORQUENTVE(2:NR-1) +
     &                       0.1*TORQUENTVE(3:NR)
      ENDDO
      ENDIF

C     COMPUTE PARTICLE FLUX GAMMA FROM NTV TORQUE DENSITY, USING
C     FLUX-FORCE RELATION. DPTRANNTV SHOULD BE INTEPRETED AS GAMMA*M_I
C     FOR BOTH IONS AND ELECTRONS, WHERE M_I IS THE ION MASS.
C     KNTV=10,20: WITH SINGLE PARTICLE SPECIES
      IF (KNTV.EQ.10.OR.KNTV.EQ.20) THEN 
      DO I=1,NR
         DPTRANNTV(I) =-TORQUENTV(I)*B0K/OMEGACI0/DPSIDSM(I)/ZCHARGE
      ENDDO
      IF (ZCHARGE.EQ.1)  DPTRANNTVI = DPTRANNTV
      IF (ZCHARGE.EQ.-1) DPTRANNTVE = DPTRANNTV
      ENDIF
C     KNTV=11,21: WITH BOTH IONS AND ELECTRONS
C     FROM PARTICLE FLUXES WITH TWO SPECIES, CHOOSE THE ONE WITH MINIMUM AMPLITUDE 
C     AS THE AMBIPOLAR CONTRIBUTION.
      IF (KNTV.EQ.11.OR.KNTV.EQ.21) THEN 
      DO I=1,NR
         DPTRANNTVI(I) =-TORQUENTVI(I)*B0K/OMEGACI0/DPSIDSM(I)/(+1)
         DPTRANNTVE(I) =-TORQUENTVE(I)*B0K/OMEGACI0/DPSIDSM(I)/(-1)
         DPTRANNTV(I)  = DPTRANNTVI(I)
         IF (ABS(DPTRANNTVI(I)).GT.ABS(DPTRANNTVE(I))) 
     &      DPTRANNTV(I)  = DPTRANNTVE(I)
      ENDDO
      ENDIF
      DO I=1,NR
         DTMP(I) =-DPTRANNTV(I)*REAL(JACOBM(I,1))
      ENDDO
      DO I=2,NR-1
         H1  = CSM(I) - CSM(I-1)
         H2  = CSM(I+1) - CSM(I)
         H12 = H1*H1
         H22 = H2*H2
         DNTRANNTV(I) = (H12*DTMP(I+1)-H22*DTMP(I-1)+(H22-H12)*DTMP(I))/
     &                  (H12*H2+H22*H1)/REAL(JACOBM(I,1))
      ENDDO
      DNTRANNTV(1)  = DNTRANNTV(2)
      DNTRANNTV(NR) = DNTRANNTV(NR-1)

C     SMOOTHING
      DO I=1,5
         DNTRANNTV(2:NR-1) = 0.1*DNTRANNTV(1:NR-2) + 
     &                       0.8*DNTRANNTV(2:NR-1) +
     &                       0.1*DNTRANNTV(3:NR)
      ENDDO

      DO I=1,NTORQ-1
         TORQUENTV(I) = 0.0*TORQUENTV(NTORQ)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) TORQUENTV(I) = 0.0
      ENDDO

      IF (KNTV.EQ.11.OR.KNTV.EQ.21) THEN
      DO I=1,NTORQ-1
         TORQUENTVI(I) = 0.0*TORQUENTVI(NTORQ)
         TORQUENTVE(I) = 0.0*TORQUENTVE(NTORQ)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) TORQUENTVI(I) = 0.0
         IF (CSM(I).GT.CTEDGE) TORQUENTVE(I) = 0.0
      ENDDO
      ENDIF

      DO I=1,NDNTR-1
         DNTRANNTV(I) = DNTRANNTV(NDNTR)
         DPTRANNTV(I) = DPTRANNTV(NDNTR)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CDEDGE) DNTRANNTV(I) = 0.0
         IF (CSM(I).GT.CDEDGE) DPTRANNTV(I) = 0.0
      ENDDO

      IF (KEYTORQ.LE.1) THEN
C     TOTAL TORQUE AND DENSITY PUMPOUT DUE TO NTV
C     NOTE THAT NTV TORQUE DENSITY SURFACE AVERAGED IN HAMADA COORDINATES
      TORQ1    = 0.0
      TDNTRNTV = 0.0
      DO I=1,NR
         TORQ1    = TORQ1    + TORQUENTV(I)*CSH(I)
         TDNTRNTV = TDNTRNTV + DNTRANNTV(I)*CSH(I)
      ENDDO
      TORQ1    = TORQ1*PI*PI*4.
      TTORQNTV = TORQ1
      TDNTRNTV = TDNTRNTV*PI*PI*4.

      IF (KNTV.EQ.11.OR.KNTV.EQ.21) THEN
      TORQI    = 0.0
      TORQE    = 0.0
      DO I=1,NR
         TORQI    = TORQI    + TORQUENTVI(I)*CSH(I)
         TORQE    = TORQE    + TORQUENTVE(I)*CSH(I)
      ENDDO
      TORQI    = TORQI*PI*PI*4.
      TORQE    = TORQE*PI*PI*4.
      ENDIF

      TORQFAC = R0EXP**3*B0EXP**2/4.0E-7/PI
      WRITE(*,120) TTORQNTV,TTORQNTV*TORQFAC
 120  FORMAT('TOTAL NTV TORQUE:',E13.6,2X,E13.6)  
 
      IF (KNTV.EQ.11.OR.KNTV.EQ.21) THEN
         WRITE(*,121) TORQI,TORQI*TORQFAC
 121     FORMAT('THERMAL ION NTV TORQUE:',E13.6,2X,E13.6)  
         WRITE(*,122) TORQE,TORQE*TORQFAC
 122     FORMAT('THERMAL ELECTRON NTV TORQUE:',E13.6,2X,E13.6)  
      ENDIF
 
      WRITE(*,112) TDNTRNTV,TDNTRNTV/TOTDENS
 112  FORMAT('TOTAL NTV DENSITY PUMPOUT:',E13.6,2X,E13.6)

      IF (KCHECK.EQ.1.AND.KEYTORQ.NE.2) THEN
C     COMPUTE OFFSET TOROIDAL ROTATION FREQUENCY
      ROTOFFSET = (ROTM(1:NR)+TROTM(1:NR)) - ROTOFFSET*ZTAUA0

C     SAVE OFFSET TOROIDAL ROTATION FREQUENCY
      OPEN(CHOUTP,FILE='PROFOFFSET.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,130) CSM(I),ROTOFFSET(I),ROTM(I)+TROTM(I)
      ENDDO
      CLOSE(CHOUTP)
 130  FORMAT(3(E15.8,1X))

C     SAVE NTV TORQUE AND PARTICLE FLUX DENSITY                 
      OPEN(CHOUTP,FILE='TORQUENTV.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,140) CSM(I),TORQUENTV(I),DPTRANNTV(I),
     &                     DPTRANNTVI(I),DPTRANNTVE(I),
     &                     TORQUENTVI(I),TORQUENTVE(I)
      ENDDO
      CLOSE(CHOUTP)
 140  FORMAT(7(E15.8,1X))
      ENDIF
      ENDIF

      IF (KEYTORQ.EQ.2) THEN
C     SMOOTHING
      DO I=1,5
         TORQCROSSAB(2:NR-1) = 0.1*TORQCROSSAB(1:NR-2) + 
     &                         0.8*TORQCROSSAB(2:NR-1) +
     &                         0.1*TORQCROSSAB(3:NR)
         TORQCROSSBA(2:NR-1) = 0.1*TORQCROSSBA(1:NR-2) + 
     &                         0.8*TORQCROSSBA(2:NR-1) +
     &                         0.1*TORQCROSSBA(3:NR)
      ENDDO

C     PATCHING
      DO I=1,NTORQ-1
         TORQCROSSAB(I) = 0.0*TORQCROSSAB(NTORQ)
         TORQCROSSBA(I) = 0.0*TORQCROSSBA(NTORQ)
      ENDDO
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) TORQCROSSAB(I) = (0.0,0.0)
         IF (CSM(I).GT.CTEDGE) TORQCROSSBA(I) = (0.0,0.0)
      ENDDO

      TORQFAC = R0EXP**3*B0EXP**2/4.0E-7/PI
      TMP = 0.0
      DO I=1,NR
         TMP = TMP + TORQCROSSAB(I)*CSH(I)
      ENDDO
      TMP    = TMP*PI*PI*4.0
      WRITE(*,115) TMP,TMP*TORQFAC
 115  FORMAT('NET TORQUE NTV_AB:',2(E13.6,1X,E13.6,3X))

      TMP = 0.0
      DO I=1,NR
         TMP = TMP + TORQCROSSBA(I)*CSH(I)
      ENDDO
      TMP    = TMP*PI*PI*4.0
      WRITE(*,116) TMP,TMP*TORQFAC
 116  FORMAT('NET TORQUE NTV_BA:',2(E13.6,1X,E13.6,3X))

C     SAVE CROSS-COUPLING TORQUE ELEMENTS INTO FILES
      OPEN(CHOUTP,FILE='TORQUENTV_AB.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,150) CSM(I),TORQCROSSAB(I)
      ENDDO
      CLOSE(CHOUTP)

      OPEN(CHOUTP,FILE='TORQUENTV_BA.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,150) CSM(I),TORQCROSSBA(I)
      ENDDO
      CLOSE(CHOUTP)
 150  FORMAT(3(E15.8,1X))
 
      DEALLOCATE(TORQCROSSAB,TORQCROSSBA)
      ENDIF

      DEALLOCATE(DVDSM,DTMP,DPTRANNTVI,DPTRANNTVE)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE TOROIDAL NTV TORQUE                           LIU YQ 03.5.2011
C TORQUE DENSITY COMPUTED AT HALF-INTEGER RADIAL MESH   LIU YQ 03.5.2011
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     
      SUBROUTINE TORQNTV1(Z_MARS,NR_MARS,NC_MARS,N_MARS,
     &                    M_MARS,KEYNTV_MARS,
     &                    S_MARS,ASPCT_MARS,CTEDGE_MARS,
     &                    DVDSM_MARS,KCHECK)
C     ==========================================================
      use NTV           ! for NTV torque calculation
      use common_ntv    ! for common variants needed for NTV calculation
      use TORQUEM
      use NTV_pre
      IMPLICIT NONE
      INCLUDE 'comioc.inc'

      INTEGER KCHECK,I,J
      REAL*8  SH,RC0,BC0,FAC,ZPATCH,ASPCT_MARS,CTEDGE_MARS
      INTEGER Z_MARS,NR_MARS,NC_MARS,N_MARS,M_MARS,KEYNTV_MARS
      REAL*8,DIMENSION(NR_MARS+1)::S_MARS
      REAL*8,DIMENSION(NR_MARS)::DVDSM_MARS
 
      TORQUENTV = 0.0
 
C     COMPUTE FLUX SURFACE AVERAGED TORQUE DENSITY
      call basic_input(Z_MARS,NR_MARS,NC_MARS,N_MARS,M_MARS,KEYNTV_MARS,
     &                 S_MARS)
      call ntv_predata    ! load data from MARS-F output
      call gntv1d         ! calculate NTV torque density
  
      call get_dimens(I,RC0,BC0)

C     NOTE THAT NTV TORQUE COMPUTED IN HAMADA COORDINATES (V,THETA,ZETA)
C     NEED TO MULTIPLY FACTOR DV/DS TO RECOVER TORQUE DENSITY IN MARS-F  
      FAC = 1.0/BC0/BC0*4e-7*npi

      IF (KEYNTV_MARS.NE.2) THEN
      DO I=1,NR_MARS
         TORQUENTV(I)=REAL(T_ntv_non(I)+T_ntv_res(I))*FAC*DVDSM_MARS(I)
      ENDDO
      ROTOFFSET = rot_diff(:,1)
      ELSE
      DO I=1,NR_MARS
         TORQCROSSAB(I)=(T_ntv_non(I)+T_ntv_res(I))*FAC*DVDSM_MARS(I)
         TORQCROSSBA(I)=CONJG(T_ntv_non(I)+T_ntv_res(I))*
     &                  FAC*DVDSM_MARS(I)
      ENDDO
      ENDIF

      IF (KCHECK.EQ.1.AND.KEYNTV_MARS.NE.2) THEN
     
      OPEN(CHOUTP,FILE='TORQUENTV2.OUT')
      REWIND(CHOUTP)
      DO I=1,NR_MARS
         SH = (S_MARS(I+1)+S_MARS(I))*0.5
         ZPATCH = 1.0
         IF (SH.GT.CTEDGE_MARS) ZPATCH = 0.0
         WRITE(CHOUTP,120) SH,TORQUENTV(I)*ZPATCH,
     &                     REAL(T_ntv_non(I))*FAC*DVDSM_MARS(I)*ZPATCH,
     &                     REAL(T_ntv_res(I))*FAC*DVDSM_MARS(I)*ZPATCH
      ENDDO
      CLOSE(CHOUTP)
C L.PIGATTO 2021 - ADD OUTPUT FILE FOR DV/DS 
      OPEN(CHOUTP,FILE='DVDSM_MARS.OUT')
      REWIND(CHOUTP)
      DO I=1,NR_MARS
         SH = (S_MARS(I+1)+S_MARS(I))*0.5
         ZPATCH = 1.0
         IF (SH.GT.CTEDGE_MARS) ZPATCH = 0.0
         WRITE(CHOUTP,120) SH,FAC,
     &                     DVDSM_MARS(I),
     &                     ZPATCH
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(4(E15.8,1X))

      OPEN(CHOUTP,FILE='PROFNTV.OUT')
      REWIND(CHOUTP)
      DO I=1,NR_MARS
         WRITE(CHOUTP,130) (nu_prof(I,J),J=1,n_prof)
      ENDDO
      CLOSE(CHOUTP)
 130  FORMAT(20(E15.8,1X))

      ENDIF

      RETURN
      END

      SUBROUTINE GET_INORMSOL(INORM)
C     ==========================================================
      USE GLOBALM
      IMPLICIT NONE
      INTEGER INORM

      INORM = INORMSOL

      RETURN
      END

      SUBROUTINE GET_ADAPCRIT(TDELTA,NONCON)
C     ==========================================================
      USE GLOBALM
      IMPLICIT NONE
      INTEGER NONCON
      REAL*8  TDELTA

      NONCON = 3
      IF (TDELTA.GT.TDELTAMAX) NONCON = 1
      IF (TDELTA.LT.TDELTAMIN) NONCON = 2
      WRITE(*,*) 'TDELTA=',TDELTA

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--PREPARE LHS COEFFICIENTS FOR MOMENTUM EQUATION       LIU YQ 10.5.2011
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE COEFFMOMENT
C     ==========================================================
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'

      INTEGER NPSTRT
      INTEGER I,J
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW0
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::CTMP
      REAL*8,DIMENSION(:),ALLOCATABLE::T1A,T1B,T1C,T1D,
     &                                 T2CA,T2CB,T2CC,T2DA,T2DB,T2DC

      REAL*8  T0A,T0B,T0D,T0ALF,T0BET
      PARAMETER(T0A=2., T0D=3., T0ALF=1.5, T0BET=2.3, T0B=T0D/T0ALF)
   
      ALLOCATE( RW0(NRP1,NCHI) ) 
      ALLOCATE( CTMP(NRP1,MEDIM) )
      ALLOCATE( T1A(NRP1), T1B(NRP1), T1C(NRP1), T1D(NRP1) )
      ALLOCATE( T2CA(NRP1), T2CB(NRP1), T2CC(NRP1),
     &          T2DA(NRP1), T2DB(NRP1), T2DC(NRP1) )
      IF (.NOT.ALLOCATED(T3A)) 
     &ALLOCATE( T3A(NRP1),T3B(NRP1),T3C(NRP1),
     &          T4A(NRP1),T4B(NRP1),T4C(NRP1), T1E(NRP1) )

      INCLUDE 'setfft.inc'
      SUBNAM    = 'Coeffmoment'

C     COMPUTE SURFACE AVERAGED EQUILIBRIUM QUANTITIES T1A-E
      DO J=1,NCHI
        DO I=2,NRP1
          RW0(I,J)=REQ(I,J)**2*G22L(I,J)/RJA(I,J)
        ENDDO 
        RW0(1,J)=RW0(2,J)
      ENDDO 
C 
      NPSTRT    =  1
      call FFTDRIVER( RW0,  CTMP,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CTMP  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,  CTMP,   NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'CTMPa')
      ENDIF
C
      DO I=1,NRP1
         T1B(I) = REAL(CTMP(I,1))
      ENDDO
 
      DO J=1,NCHI
        DO I=2,NRP1
          RW0(I,J)=REQ(I,J)**2*RJA(I,J)
        ENDDO 
        RW0(1,J)=RW0(2,J)
      ENDDO 
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,  CTMP,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CTMP  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,  CTMP,   NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'CTMPb')
      ENDIF
C
      DO I=2,NRP1
C        T1C(I)=1. + DPSIDS(I)**2/RHO(I)*T1B(I)/REAL(CTMP(I,1))
         T1C(I)=1. 
      ENDDO
      T1C(1)=1.
   
      DO I=2,NRP1
         T1E(I)=RHO(I)*REAL(CTMP(I,1))/REAL(JACOBI(I,1))
      ENDDO
      T1E(1)=T1E(2)

      DO I=2,NRP1
         T1B(I)=TCHIMI(I)*T1B(I)/REAL(JACOBI(I,1))
      ENDDO
      T1B(1)=T1B(2)

      DO J=1,NCHI
        DO I=2,NRP1
          RW0(I,J)=REQ(I,J)*SQRT(G22L(I,J))
        ENDDO 
        RW0(1,J)=RW0(2,J)
      ENDDO 
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,  CTMP,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CTMP  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,  CTMP,   NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'CTMPc')
      ENDIF
C
      DO I=2,NRP1
         T1D(I) = TVPINCHI(I)*REAL(CTMP(I,1))/REAL(JACOBI(I,1))
      ENDDO
      T1D(1)=T1D(2)
 
      DO I=2,NRP1
         T1A(I)=CS(I)*REAL(JACOBI(I,1))/Q(I)/DPSIDS(I)
      ENDDO
      T1A(1)=0.
     
C     SPECIFY A TEST CASE FOR THE MOMENTUM SOLVER
      IF (NCASE.EQ.8) THEN
         DO I=1,NRP1
            T1A(I) = T0A*EXP(T0BET*CS(I))
            T1B(I) = T0B
            T1C(I) = 1.1 + CS(I)**2
            T1D(I) = T0D
            T1E(I) = 1.
         ENDDO
      ENDIF

C     COMPUTE COEFFICIENTS T2*
      DO I=1,NR
         RW0(I,1)=(T1A(I)+T1A(I+1))*0.5      
         RW0(I,2)=(T1B(I)+T1B(I+1))*0.5      
         RW0(I,3)=(T1C(I)+T1C(I+1))*0.5      
         RW0(I,4)=(T1D(I)+T1D(I+1))*0.5      
         RW0(I,5)=(T1A(I+1)-T1A(I))/RW0(I,1)
      ENDDO
      
      DO I=2,NR
         T2CA(I)=CSH(I-1)*RW0(I-1,3)/6.
         T2CB(I)=(CSH(I-1)*RW0(I-1,3)+CSH(I)*RW0(I,3))/3.
         T2CC(I)=CSH(I)*RW0(I,3)/6.
         T2DA(I)=(CSH(I-1)*RW0(I-1,4)*(RW0(I-1,5)/6.-0.5) + 
     &            RW0(I-1,2)*(1.-0.5*RW0(I-1,5)))/CSH(I-1)
         T2DB(I)=(CSH(I-1)*RW0(I-1,4)*(RW0(I-1,5)/3.-0.5) - 
     &            RW0(I-1,2)*(1.-0.5*RW0(I-1,5)))/CSH(I-1) +
     &           (CSH(I)*RW0(I,4)*(RW0(I,5)/3.+0.5) -
     &            RW0(I,2)*(1.+0.5*RW0(I,5)))/CSH(I)
         T2DC(I)=(CSH(I)*RW0(I,4)*(RW0(I,5)/6.+0.5) +
     &            RW0(I,2)*(1.+0.5*RW0(I,5)))/CSH(I)
      ENDDO

      I = 1
      T2CA(I)=0.
      T2CB(I)=CSH(I)*RW0(I,3)/3.
      T2CC(I)=CSH(I)*RW0(I,3)/6.
      T2DA(I)=0.
      T2DB(I)=(CSH(I)*RW0(I,4)*(RW0(I,5)/3.+0.5) -
     &         RW0(I,2)*(1.+0.5*RW0(I,5)))/CSH(I) - T1D(I)
      T2DC(I)=(CSH(I)*RW0(I,4)*(RW0(I,5)/6.+0.5) +
     &         RW0(I,2)*(1.+0.5*RW0(I,5)))/CSH(I)

      I = NRP1
      T2CA(I)=CSH(I-1)*RW0(I-1,3)/6.
      T2CB(I)=CSH(I-1)*RW0(I-1,3)/3.
      T2CC(I)=0.
      T2DA(I)=(CSH(I-1)*RW0(I-1,4)*(RW0(I-1,5)/6.-0.5) + 
     &         RW0(I-1,2)*(1.-0.5*RW0(I-1,5)))/CSH(I-1)
      T2DB(I)=(CSH(I-1)*RW0(I-1,4)*(RW0(I-1,5)/3.-0.5) - 
     &         RW0(I-1,2)*(1.-0.5*RW0(I-1,5)))/CSH(I-1) + T1D(I) 
      T2DC(I)=0.
      
C     COMPUTE COEFFICIENTS T3A-C AND T4A-C
      T3A = CALPHA5*T2CA - CALPHA6*T2DA
      T3B = CALPHA5*T2CB - CALPHA6*T2DB
      T3C = CALPHA5*T2CC - CALPHA6*T2DC
      T4A = CALPHA5*T2CA + (1.-CALPHA6)*T2DA
      T4B = CALPHA5*T2CB + (1.-CALPHA6)*T2DB
      T4C = CALPHA5*T2CC + (1.-CALPHA6)*T2DC

C
      DEALLOCATE(RW0,CTMP,
     &           T1A,T1B,T1C,T1D,T2CA,T2CB,T2CC,T2DA,T2DB,T2DC)

C     PREPARE FOR TIME-STEPPING OF TOROIDAL MOMENTUM EQUATION
      IF (.NOT.ALLOCATED(TMOM)) THEN
         ALLOCATE(TMOM(NRP1), TMOMOLD(NRP1), TITALF(NRP1) )
         TMOMOLD   = 0.
      ENDIF

      TITALF(1) =0.
      DO I=1,NR
         TITALF(I+1)=-T3C(I)/(T3A(I)*TITALF(I)+T3B(I))
      ENDDO

      IF (NCASE.EQ.8) THEN
         IF (.NOT.ALLOCATED(TORQUEJXB)) ALLOCATE(TORQUEJXB(NR))
         IF (.NOT.ALLOCATED(TORQUENTV)) ALLOCATE(TORQUENTV(NR))
         IF (.NOT.ALLOCATED(TORQUEREY)) ALLOCATE(TORQUEREY(NR))
         IF (.NOT.ALLOCATED(DENSPUMP))  ALLOCATE(DENSPUMP(NR))
         IF (.NOT.ALLOCATED(DISPNORM))  ALLOCATE(DISPNORM(NCHI))
         IF (.NOT.ALLOCATED(DPTRANMHD)) ALLOCATE(DPTRANMHD(NR))
         IF (.NOT.ALLOCATED(DPTRANNTV)) ALLOCATE(DPTRANNTV(NR))
         TORQUEJXB = 0.
         TORQUENTV = 0.
         TORQUEREY = 0.
         DENSPUMP = 0.
         DISPNORM = (0.,0.)
         DPTRANMHD= 0.
         DPTRANNTV= 0.
         TMOMOLD = -2.
      ENDIF

      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--TIME-STEPPING MOMENTUM EQUATION                      LIU YQ 09.5.2011
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE SOLVEMOMENT
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    I
      REAL*8     H1,H2
      REAL*8,DIMENSION(:),ALLOCATABLE::TTORQ,TITBET   

      REAL*8  T0T,T0GAM
      PARAMETER(T0T=3.2, T0GAM=1.7)

      ALLOCATE(TTORQ(NRP1),TITBET(NRP1+1))

C     COMPUTE TOTAL TORQUE
      TITBET = 0.
      IF (NCASE.EQ.8) THEN
         DO I=1,NR
            TITBET(I) = T0T*EXP(T0GAM*CSM(I))
         ENDDO
      ELSE
         DO I=1,NR
            TITBET(I) = CTJXB*TORQUEJXB(I) + CTNTV*TORQUENTV(I) +
     &                  CTREY*TORQUEREY(I) + CTERGO*TORQUEERGO(I)
         ENDDO
      ENDIF

C     COMPUTE TOTAL TORQUE WITH FEM DISCRETISATION
      DO I=2,NR     
         TTORQ(I)= (CSH(I-1)*TITBET(I-1)+CSH(I)*TITBET(I))*0.5
      ENDDO

      I = 1
      TTORQ(I) = CSH(I)*TITBET(I)*0.5

      I = NRP1
      TTORQ(I) = CSH(I-1)*TITBET(I-1)*0.5 

      TTORQ = TTORQ*CALPHA7

C     COMPUTE RHS OF DISCRETISED MOMENTUM EQUATION
      DO I=2,NR
         TTORQ(I) = TTORQ(I) + 
     &   T4A(I)*TMOMOLD(I-1) + T4B(I)*TMOMOLD(I) + T4C(I)*TMOMOLD(I+1)
      ENDDO
      I=1
      TTORQ(I) = TTORQ(I) + T4B(I)*TMOMOLD(I) + T4C(I)*TMOMOLD(I+1)
      I=NRP1
      TTORQ(I) = TTORQ(I) + T4A(I)*TMOMOLD(I-1) + T4B(I)*TMOMOLD(I)

C     FORWARD LOOP TO COMPUTE TITBET
      TITBET(1) = 0.
      DO I=1,NRP1
         TITBET(I+1) = (TTORQ(I)-T3A(I)*TITBET(I))/
     &                 (T3A(I)*TITALF(I)+T3B(I))
      ENDDO

      IF (.FALSE.) THEN
      WRITE(*,*) 'TEST SOLVEMOMENTUM: TORQUE'
      DO I=1,NR
      WRITE(*,105) CSM(I),TORQUEJXB(I),TORQUENTV(I),TORQUEREY(I)
      ENDDO

      WRITE(*,*) 'TEST SOLVEMOMENTUM: T3'
      DO I=1,NRP1
      WRITE(*,105) CS(I),T3A(I),T3B(I),T3C(I)
      ENDDO

      WRITE(*,*) 'TEST SOLVEMOMENTUM: T4'
      DO I=1,NRP1
      WRITE(*,105) CS(I),T4A(I),T4B(I),T4C(I)
      ENDDO

      WRITE(*,*) 'TEST SOLVEMOMENTUM: TI'
      DO I=1,NRP1
         WRITE(*,105) CS(I),TTORQ(I),TITALF(I),TITBET(I)
      ENDDO
      WRITE(*,*) TITBET(NRP1+1)
 105  FORMAT(4(E15.7,1X)) 
      ENDIF

C     BACKWARD LOOP TO COMPUTE THE NEW MOMENTUM SOLUTION
      IF (IT_BC.EQ.2) TITBET(NRP1+1) = 0.
      TMOM(NRP1) = TITBET(NRP1+1)
      DO I=NRP1,2,-1
         TMOM(I-1) = TITALF(I)*TMOM(I) + TITBET(I)
      ENDDO
      
C     COMPUTE MODIFIED TOROIDAL ROTATION FREQUENCY
      DO I=1,NR
         TROTI(I) = TMOM(I)/T1E(I)
      ENDDO
      TROTI(NRP1) = TROTI(NR)
      IF (IT_BC.EQ.2) TROTI(NRP1) = 0.

C     SET SATURATION CONDITION WHEN TOTAL ROTATION VANISHES
      IF (ITSATURAT.EQ.1) THEN
      DO I=1,NRP1
         IF (ROT(I)*(ROT(I)+TROTI(I)).LT.0.) THEN
            TROTI(I) = -ROT(I)+1.e-3
            TMOM(I)  = TROTI(I)*T1E(I)
         ENDIF
      ENDDO
      ENDIF

      DO I=1,NR
         TROTM(I) = (TROTI(I)+TROTI(I+1))*0.5
      ENDDO
      TROTM(NRP1) = 0.

      IF (ITSATURAT.EQ.1) THEN
      DO I=1,NR
         IF (ROTM(I)*(ROTM(I)+TROTM(I)).LT.0.) THEN
            TROTM(I) = -ROTM(I)+1.e-3
         ENDIF
      ENDDO
      ENDIF

C     COMPUTE MODIFIED TOROIDAL ROTATION SHEAR
      DO I = 1,NR
        TDROTM(I) = (TROTI(I+1)-TROTI(I))/(CS(I+1) - CS(I))
      END DO
      DO I = 2,NR
         H1 = (CS(I)-CS(I-1))/2
         H2 = (CS(I+1)-CS(I))/2
         TDROTI(I) = (H1/H2*TROTM(I)-H2/H1*TROTM(I-1))/(H1+H2)-
     &             (H1-H2)*TROTI(I)/H1/H2
      END DO
      TDROTI(1)    = 0.0
      TDROTI(NRP1) = (TROTI(NRP1)-TROTM(NR))/(CS(NRP1)-CSM(NR))
      TDROTM(NRP1) = 0.

C     COPY NEW SOLUTION TO OLD ONE
      TMOMOLD = TMOM

      DEALLOCATE(TTORQ, TITBET)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--PREPARE LHS COEFFICIENTS FOR DENSITY TRANSPORT EQ.   LIU YQ 12.2.2016
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE COEFFDNTRAN
C     ==========================================================
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'

      INTEGER NPSTRT
      INTEGER I,J
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW0
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::CTMP
      REAL*8,DIMENSION(:),ALLOCATABLE::T1A,T1B,
     &                                 T2CA,T2CB,T2CC,T2DA,T2DB,T2DC

      ALLOCATE( RW0(NRP1,NCHI) ) 
      ALLOCATE( CTMP(NRP1,MEDIM) )
      ALLOCATE( T1A(NRP1), T1B(NRP1) )
      ALLOCATE( T2CA(NRP1), T2CB(NRP1), T2CC(NRP1),
     &          T2DA(NRP1), T2DB(NRP1), T2DC(NRP1) )
      IF (.NOT.ALLOCATED(D3A)) 
     &ALLOCATE( D3A(NRP1),D3B(NRP1),D3C(NRP1),
     &          D4A(NRP1),D4B(NRP1),D4C(NRP1) )

      INCLUDE 'setfft.inc'
      SUBNAM    = 'Coeffmoment'

C     COMPUTE SURFACE AVERAGED EQUILIBRIUM QUANTITIES T1A-E
      DO J=1,NCHI
        DO I=2,NRP1
          RW0(I,J)=REQ(I,J)**2*G22L(I,J)/RJA(I,J)
        ENDDO 
        RW0(1,J)=RW0(2,J)
      ENDDO 
C 
      NPSTRT    =  1
      call FFTDRIVER( RW0,  CTMP,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: CTMP  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,  CTMP,   NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'CTMPa')
      ENDIF
C
      DO I=1,NRP1
         T1B(I) = REAL(CTMP(I,1))
      ENDDO
 
      DO I=2,NRP1
         T1B(I)=TCHIDI(I)*T1B(I)/REAL(JACOBI(I,1))
      ENDDO
      T1B(1)=T1B(2)

      DO I=2,NRP1
         T1A(I)=REAL(JACOBI(I,1))
      ENDDO
      T1A(1)=0.
     
C     COMPUTE COEFFICIENTS T2*
      DO I=1,NR
         RW0(I,1)=(T1A(I)+T1A(I+1))*0.5      
         RW0(I,2)=(T1B(I)+T1B(I+1))*0.5      
         RW0(I,5)=(T1A(I+1)-T1A(I))/RW0(I,1)
      ENDDO
      
      DO I=2,NR
         T2CA(I)=CSH(I-1)/6.
         T2CB(I)=(CSH(I-1)+CSH(I))/3.
         T2CC(I)=CSH(I)/6.
         T2DA(I)= RW0(I-1,2)*(1.-0.5*RW0(I-1,5))/CSH(I-1)
         T2DB(I)=-RW0(I-1,2)*(1.-0.5*RW0(I-1,5))/CSH(I-1)
     &           -RW0(I,2)*(1.+0.5*RW0(I,5))/CSH(I)
         T2DC(I)=RW0(I,2)*(1.+0.5*RW0(I,5))/CSH(I)
      ENDDO

      I = 1
      T2CA(I)=0.
      T2CB(I)=CSH(I)/3.
      T2CC(I)=CSH(I)/6.
      T2DA(I)=0.
      T2DB(I)=-RW0(I,2)*(1.+0.5*RW0(I,5))/CSH(I)
      T2DC(I)= RW0(I,2)*(1.+0.5*RW0(I,5))/CSH(I)

      I = NRP1
      T2CA(I)=CSH(I-1)/6.
      T2CB(I)=CSH(I-1)/3.
      T2CC(I)=0.
      T2DA(I)= RW0(I-1,2)*(1.-0.5*RW0(I-1,5))/CSH(I-1)
      T2DB(I)=-RW0(I-1,2)*(1.-0.5*RW0(I-1,5))/CSH(I-1) 
      T2DC(I)=0.
      
C     COMPUTE COEFFICIENTS T3A-C AND T4A-C
      D3A = CALPHA5*T2CA - CALPHA6*T2DA
      D3B = CALPHA5*T2CB - CALPHA6*T2DB
      D3C = CALPHA5*T2CC - CALPHA6*T2DC
      D4A = CALPHA5*T2CA + (1.-CALPHA6)*T2DA
      D4B = CALPHA5*T2CB + (1.-CALPHA6)*T2DB
      D4C = CALPHA5*T2CC + (1.-CALPHA6)*T2DC

C
      DEALLOCATE(RW0,CTMP,T1A,T1B,T2CA,T2CB,T2CC,T2DA,T2DB,T2DC)

C     PREPARE FOR TIME-STEPPING OF DENSITY TRANSPORT EQUATION
      IF (.NOT.ALLOCATED(TRHOIOLD)) THEN
         ALLOCATE(TRHOIOLD(NRP1), DITALF(NRP1) )
         TRHOIOLD  = 0.
      ENDIF

      DITALF(1) =0.
      DO I=1,NR
         DITALF(I+1)=-D3C(I)/(D3A(I)*DITALF(I)+D3B(I))
      ENDDO

      RETURN
      END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--TIME-STEPPING RADIAL DENSITY TRANSPORT EQUATION      LIU YQ 12.2.2016
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE SOLVEDNTRAN
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    I
      REAL*8,DIMENSION(:),ALLOCATABLE::TTORQ,TITBET   

      ALLOCATE(TTORQ(NRP1),TITBET(NRP1+1))

C     COMPUTE TOTAL DENSITY 
      TITBET = 0.
      IF (CDMHD.GE.0..AND.CDNTV.GE.0.) THEN
      DO I=1,NR
         TITBET(I) = CDMHD*DNTRANMHD(I) + CDNTV*DNTRANNTV(I) + 
     &               CDERGO*DNTRANERGO(I)
      ENDDO

      ELSE

      DO I=1,NR
         TITBET(I) =-5.0E-5*CSM(I)**20*(1.-CSM(I)**2)
      ENDDO
      ENDIF

C     COMPUTE TOTAL DENSITY WITH FEM DISCRETISATION
      DO I=2,NR     
         TTORQ(I)= (CSH(I-1)*TITBET(I-1)+CSH(I)*TITBET(I))*0.5
      ENDDO

      I = 1
      TTORQ(I) = CSH(I)*TITBET(I)*0.5

      I = NRP1
      TTORQ(I) = CSH(I-1)*TITBET(I-1)*0.5 

      TTORQ = TTORQ*CALPHA7

C     COMPUTE RHS OF DISCRETISED DENSITY TRANSPORT EQUATION
      DO I=2,NR
         TTORQ(I) = TTORQ(I) + 
     &   D4A(I)*TRHOIOLD(I-1)+D4B(I)*TRHOIOLD(I)+D4C(I)*TRHOIOLD(I+1)
      ENDDO
      I=1
      TTORQ(I) = TTORQ(I)+D4B(I)*TRHOIOLD(I)+D4C(I)*TRHOIOLD(I+1)
      I=NRP1
      TTORQ(I) = TTORQ(I)+D4A(I)*TRHOIOLD(I-1)+D4B(I)*TRHOIOLD(I)

C     FORWARD LOOP TO COMPUTE TITBET
      TITBET(1) = 0.
      DO I=1,NRP1
         TITBET(I+1) = (TTORQ(I)-D3A(I)*TITBET(I))/
     &                 (D3A(I)*DITALF(I)+D3B(I))
      ENDDO

C     BACKWARD LOOP TO COMPUTE THE NEW DENSITY SOLUTION
      IF (ID_BC.EQ.2) TITBET(NRP1+1) = 0.
      TRHOI(NRP1) = TITBET(NRP1+1)
      DO I=NRP1,2,-1
         TRHOI(I-1) = DITALF(I)*TRHOI(I) + TITBET(I)
      ENDDO
      
C     ALWAYS SET SATURATION CONDITION WHEN TOTAL DENSITY VANISHES
      DO I=1,NRP1
         IF ((RHO(I)+TRHOI(I)).LT.0.) THEN
            TRHOI(I) = -RHO(I)+5.e-2
         ENDIF
      ENDDO

      DO I=1,NR
         TRHOM(I) = (TRHOI(I)+TRHOI(I+1))*0.5
      ENDDO
      TRHOM(NRP1) = 0.

C     COPY NEW SOLUTION TO OLD ONE
      TRHOIOLD = TRHOI

      DEALLOCATE(TTORQ, TITBET)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE DENSITY PUMP-OUT                             LIU YQ 06.23.2011
C AT HALF-INTEGER RADIAL MESH   
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE CALCDENS
C     ==========================================================
      USE DIMENSIM
      USE RCOMDM  
      USE GLOBALM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'compam.inc'
      INTEGER    I,J,MS,KCHECK
      REAL*8     RCHIH,RTMP2
      COMPLEX*16 CTMP1

      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::RW1,RW2,RW3
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::DENSFI,DENSFM
      REAL*8,DIMENSION(:),ALLOCATABLE::RTMP1

      ALLOCATE( RW1(NR,NCHI), RW2(NR,NCHI), RW3(NR,NCHI),
     &          DENSFI(NRP1), DENSFM(NR),
     &          RTMP1(NRP1) )

      KCHECK = 0

      IF (.NOT.ALLOCATED(DENSPUMP)) THEN
         ALLOCATE(DENSPUMP(NR))
         DENSPUMP  = 0.0
         DENSPUMPA = 0.0
      ENDIF
 
C     COMPUTE LOG(P0/RHO0^GAMMA)
      DO I=1,NR
         RTMP1(I) = LOG(PEQ(I)/RHO(I)**GAMMA)
      ENDDO
      RTMP1(NRP1) = 2*RTMP1(NR)-RTMP1(NR-1)

C     CONVERT PERTURBED PRESSURE, V1U AND X1U IN REAL SPACE
C     AT HALF-INTEGER RADIAL GRID
      RW1 = (0.,0.)
      RW2 = (0.,0.)
      RW3 = (0.,0.)
      RCHIH = 2.*PI/NCHI
      DO I=1,NR
         DO J=1,NCHI
            DO MS=1,MSMAX
               CTMP1 = EXP(CI*RM(MS,2)*(J-1)*RCHIH)
               RW1(I,J) = RW1(I,J) + (PRE(I,MS)+PEE(I,MS))*CTMP1
               RW2(I,J) = RW2(I,J) + (X1U(I,MS)+X1U(I+1,MS))*.5*CTMP1
               RW3(I,J) = RW3(I,J) + (V1U(I,MS)+V1U(I+1,MS))*.5*CTMP1
            ENDDO
         ENDDO       
      ENDDO

C     COMPUTE J*P*V1U IN REAL SPACE
      DO I=1,NR
         DO J=1,NCHI
            RW1(I,J) = RJAM(I,J)*RHOM(I)/GAMMA*(RW1(I,J)/PEQM(I)+
     &         RW2(I,J)*(RTMP1(I+1)-RTMP1(I))/CSH(I))*CONJG(RW3(I,J))
         ENDDO
      ENDDO

C     COMPUTE DENSITY FLUX DENSFM=(J*P*V1U)_M=0
C     AT HALF-INTEGER RADIAL GRID
      DENSFM = (0.,0.)
      DO I=1,NR
         DO J=1,NCHI
            DENSFM(I) = DENSFM(I) + RW1(I,J)
         ENDDO
      ENDDO
      DENSFM = DENSFM*RCHIH*.5/PI

C     COMPUTE DENSITY FLUX AT INTEGER RADIAL GRID
      DENSFI(1) = DENSFM(1)
      DO I=2,NR
         DENSFI(I) = DENSFM(I-1) + (DENSFM(I)-DENSFM(I-1))*
     &               CSH(I-1)/(CSH(I-1)+CSH(I))
      ENDDO
      DENSFI(NRP1) = DENSFM(NR)

C     COMPUTE RADIAL DISTRIBUTION OF PUMPED OUT DENSITY
C     AS FUNCTION OF TIME
      DO I=1,NR
         DENSPUMP(I) = DENSPUMP(I) - REAL((DENSFI(I+1)-DENSFI(I))/
     &                 CSH(I)/JACOBM(I,1)/CALPHA1)
      ENDDO

C     COMPUTE TOTAL DENSITY PUMP-OUT, AS FRACTION OF TOTAL INITIAL DENSITY
      RTMP2 = 0.
      DO I=1,NR
         RTMP2 = RTMP2 + RHOM(I)*JACOBM(I,1)*CSH(I)
      ENDDO
      DENSPUMPA = DENSPUMPA - DENSFI(NRP1)/RTMP2/CALPHA1

      DEALLOCATE(RW1,RW2,RW3,DENSFI,DENSFM,RTMP1)
      
      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE NORMAL DISPLACEMENT AT PLASMA SURFACE        LIU YQ 05.26.2012
C ALONG POLOIDAL ANGLE          
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE CALCDISPNORM
C     ==========================================================
      USE DIMENSIM
      USE RCOMDM  
      USE GLOBALM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'
      INTEGER    ID,J,MS,KCHECK
      REAL*8     RCHIH
      COMPLEX*16 CTMP1

      KCHECK = 0
      IF (NCASE.EQ.1.OR.NCASE.EQ.2) KCHECK = 1

      IF (.NOT.ALLOCATED(DISPNORM)) ALLOCATE(DISPNORM(NCHI))
 
C     FIND THE CTEDGE-SURFACE WHERE THE NORMAL DISPLACEMENT IS COMPUTED
      ID = 1
      DO J=1,NRP1
         IF (CS(J)<CTEDGE) ID = J
      ENDDO

C     COMPUTE NORMAL DISPLACEMENT ALONG POLOIDAL ANGLE
      DISPNORM = (0.,0.)
      RCHIH    = 2.*PI/NCHI
      DO J=1,NCHI
         DO MS=1,MSMAX
            CTMP1       = EXP(CI*RM(MS,2)*(J-1)*RCHIH)
            DISPNORM(J) = DISPNORM(J) + X1U(ID,MS)*CTMP1
         ENDDO
         DISPNORM(J) = DISPNORM(J)*RJA(ID,J)/REQ(ID,J)/SQRT(G22L(ID,J))
      ENDDO

C     SAVE NORMAL DISPLACEMENT
      IF (KCHECK.EQ.1) THEN
      OPEN(CHOUTP,FILE='PROFDISP.OUT')
      REWIND(CHOUTP)
      DO J=1,NCHI
         WRITE(CHOUTP,130) RCHIH*(J-1),DISPNORM(J)
      ENDDO
      CLOSE(CHOUTP)
 130  FORMAT(3(E15.8,1X))
      ENDIF

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C SAVE WHOLE SOLUTION FOR NCASE=6 RUNS                 LIU YQ 06.06.2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SOLSAVE(X,Y,CA1,CA2,CA3,CA4,RA1)
      USE DIMENSIM
      USE GLOBALM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      COMPLEX*16 X(*),Y(*)
      COMPLEX*16 CA1,CA2,CA3,CA4
      REAL*8     RA1
      INTEGER    I,J
 
      OPEN(CHOUTP,FILE='MARSQ.SOL')
      REWIND(CHOUTP)
      WRITE(CHOUTP,130) RA1,CA1,CA2,CA3,CA4
      DO I=1,NRP1
         WRITE(CHOUTP,130) TMOMOLD(I),TROTI(I),TROTM(I),
     &                     TDROTI(I),TDROTM(I),
     &                     TRHOIOLD(I),TRHOI(I),TRHOM(I)
      ENDDO
      DO I=1,NXCOMP*MSMAX*NTP1
         WRITE(CHOUTP,130) X(I)
      ENDDO
      DO I=1,NYCOMP*MSMAX*NTP1
         WRITE(CHOUTP,130) Y(I)
      ENDDO
      CLOSE(CHOUTP)
 130  FORMAT(9(E18.10E3,1X))

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C READ WHOLE SOLUTION FOR NCASE=6 RUNS                 LIU YQ 06.06.2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SOLREAD(X,Y,CA1,CA2,CA3,CA4,RA1)
      USE DIMENSIM
      USE GLOBALM
      USE TORQUEM
      USE FEEDBACKM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      COMPLEX*16 X(*),Y(*)
      COMPLEX*16 CA1,CA2,CA3,CA4
      REAL*8     RA1,TMP1R,TMP1I,TMP2R,TMP2I,TMP3R,TMP3I,TMP4R,TMP4I
      INTEGER    I,J
 
      OPEN(CHOUTP,FILE='MARSQ.SOL')
      REWIND(CHOUTP)
      READ(CHOUTP,*) RA1,TMP1R,TMP1I,TMP2R,TMP2I,TMP3R,TMP3I,TMP4R,TMP4I
      CA1 = CMPLX(TMP1R,TMP1I)
      CA2 = CMPLX(TMP2R,TMP2I)
      CA3 = CMPLX(TMP3R,TMP3I)
      CA4 = CMPLX(TMP4R,TMP4I)
      DO I=1,NRP1
         READ(CHOUTP,*) TMOMOLD(I),TROTI(I),TROTM(I),
     &                  TDROTI(I),TDROTM(I),
     &                  TRHOIOLD(I),TRHOI(I),TRHOM(I)
      ENDDO
      DO I=1,NXCOMP*MSMAX*NTP1
         READ(CHOUTP,*) TMP1R,TMP1I
         X(I) = CMPLX(TMP1R,TMP1I)*THRESHOLD
      ENDDO
      DO I=1,NYCOMP*MSMAX*NTP1
         READ(CHOUTP,*) TMP1R,TMP1I
         Y(I) = CMPLX(TMP1R,TMP1I)*THRESHOLD
      ENDDO
      CLOSE(CHOUTP)

      RETURN
      END


C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C CHECK WHETHER SOLUTION SATISFIES EQUATIONS           LIU YQ 09.27.2022
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SOLTEST(MD,MDY,ND,NXC,NYC,A,B,C,D,E,F,G,H,X,Y)
      USE DIMENSIM
      USE GLOBALM
      IMPLICIT NONE
      INCLUDE 'cdipam.inc'
      INCLUDE 'comioc.inc'

      COMPLEX*16 CV1,CV2,CV3
      REAL*8     CV1M,CV2M,CV3M,CV1X,CV2X,CV3X,CXYM
      INTEGER    I,MROW,MSA,LXROW,LYROW,LXCOL,LYCOL,KX,KY
 
      OPEN(CHOUTP,FILE='SOLTEST.OUT')
      REWIND(CHOUTP)
      WRITE(CHOUTP,100) NR-1,MSMAX,0

      CV1M = 0.0
      CV2M = 0.0
      CV3M = 0.0
      CV1X = 0.0
      CV2X = 0.0
      CV3X = 0.0

      CXYM = (SUM(ABS(X)) + SUM(ABS(Y)))/(MD*(ND+1)+MDY*ND)

      DO I=2,NR
      DO MROW=1,MSMAX
      LXROW = (MROW-1)*NXC
      LYROW = (MROW-1)*NYC

      CV1 = (0.0,0.0)
      CV2 = (0.0,0.0)
      CV3 = (0.0,0.0)

      DO MSA=1,MSMAX
      LXCOL = (MSA -1)*NXC
      LYCOL = (MSA -1)*NYC
         
      DO KX=1,NXC
         CV1 = CV1 + A(KXV1+LXROW,KX+LXCOL,I)*X(KX+LXCOL,I-1)
     &             + B(KXV1+LXROW,KX+LXCOL,I)*X(KX+LXCOL,I)
     &             + C(KXV1+LXROW,KX+LXCOL,I)*X(KX+LXCOL,I+1)
      ENDDO
      DO KY=1,NYC
         CV1 = CV1 + H(KXV1+LXROW,KY+LYCOL,I)*Y(KY+LYCOL,I-1)
     &             + E(KXV1+LXROW,KY+LYCOL,I)*Y(KY+LYCOL,I)
      ENDDO

      DO KY=1,NYC
         CV2 = CV2 + D(KYV2+LYROW,KY+LYCOL,I)*Y(KY+LYCOL,I)
         CV3 = CV3 + D(KYV3+LYROW,KY+LYCOL,I)*Y(KY+LYCOL,I)
      ENDDO
      DO KX=1,NXC
         CV2 = CV2 + F(KYV2+LYROW,KX+LXCOL,I)*X(KX+LXCOL,I)
     &             + G(KYV2+LYROW,KX+LXCOL,I)*X(KX+LXCOL,I+1)
         CV3 = CV3 + F(KYV3+LYROW,KX+LXCOL,I)*X(KX+LXCOL,I)
     &             + G(KYV3+LYROW,KX+LXCOL,I)*X(KX+LXCOL,I+1)
      ENDDO
      ENDDO

      CV1 = CV1/CXYM
      CV2 = CV2/CXYM
      CV3 = CV3/CXYM

      WRITE(CHOUTP,110) ABS(CV1),ABS(CV2),ABS(CV3)

      IF (I.LE.NFIT) CV1 = (0.0,0.0)

      CV1M = CV1M + ABS(CV1)
      CV2M = CV2M + ABS(CV2)
      CV3M = CV3M + ABS(CV3)

      IF (CV1X.LT.ABS(CV1)) CV1X = ABS(CV1)
      IF (CV2X.LT.ABS(CV2)) CV2X = ABS(CV2)
      IF (CV3X.LT.ABS(CV3)) CV3X = ABS(CV3)

      ENDDO
      ENDDO

      CLOSE(CHOUTP)
 100  FORMAT(3(I4,1X))
 110  FORMAT(3(E11.4,1X))

      CV1M = CV1M/(NR-1)/MSMAX
      CV2M = CV2M/(NR-1)/MSMAX
      CV3M = CV3M/(NR-1)/MSMAX

      WRITE(*,'("MEAN ERROR = ",3E12.4)') CV1M,CV2M,CV3M
      WRITE(*,'(" MAX ERROR = ",3E12.4)') CV1X,CV2X,CV3X

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE CHI(THETA) AT UNIFORM GEOMETRIC ANGLE THETA  LIU YQ 19.12.2012
C AT ALL RATIONAL SURFACES, AND SAVE TO A FILE                         $
C ASSUMING RATIONAL SURFACES (IRATSURF) ALREADY COMPUTED               $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE THETA2CHI
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INTEGER    I,J,IR,K,M,MS
      REAL*8     TQHCHI,QTMP1,QTMP2,REQS,ZEQS
      REAL*8,DIMENSION(:),ALLOCATABLE::TQTET,TQCHI,TQTETN,TQCHIN,WORK
 
      ALLOCATE(TQTET(3*NCHI+1),TQCHI(3*NCHI+1),
     &         TQTETN(NCHI+1), TQCHIN(NCHI+1),
     &         WORK(3*NCHI+1))

      TQTET  = 0.
      TQCHI  = 0.
      TQTETN = 0.
      TQCHIN = 0.

C     UNIFORM CHI-MESH FROM CHEASE, EXTENDED TO [-2*PI,4*PI]
      TQHCHI   = 2.*PI/NCHI
      TQCHI(1) = -2.*PI
      DO J=2,3*NCHI+1
         TQCHI(J) = TQCHI(J-1) + TQHCHI
      ENDDO

      TQTETN(1) = 0.
      DO J=2,NCHI+1
         TQTETN(J) = TQTETN(J-1) + TQHCHI
      ENDDO

C     LOOPING THROUGH ALL RATIONAL SURFACES
      OPEN(CHOUTP,FILE='PANGLE.OUT')
      REWIND(CHOUTP)
      WRITE(CHOUTP,*) NRATSURF
      WRITE(CHOUTP,*) NCHI+1
      
      DO IR=1,NRATSURF
      I = IRATSURF(IR)

C     COMPUTE CORRESPONDING THETA-MESH, USING (REQ,ZEQ)
      M = 1
      DO MS=1,MSMAX
         IF (ABS(-QPLS(I)*RNTOR-RM(MS,2)).LT.0.90) M=MS
      ENDDO
      QTMP1 = -RM(M,2)/RNTOR
      QTMP2 = (QTMP1-QPLS(I))/(QPLS(I+1)-QPLS(I))   
      IF (ABS(QTMP2).GT.1.) THEN
         WRITE(*,*) QPLS(I),QPLS(I+1),-RM(M,2),QTMP1,QTMP2
         STOP 'THETA2CHI: INACCURATE MESH'
      ENDIF

      DO J=1,NCHI
         REQS = REQ(I,J) + QTMP2*(REQ(I+1,J)-REQ(I,J))
         ZEQS = ZEQ(I,J) + QTMP2*(ZEQ(I+1,J)-ZEQ(I,J))
         TQTET(NCHI+J) = DATAN2(ZEQS-ZEQ(1,1),REQS-REQ(1,1))
      ENDDO
      TQTET(2*NCHI+1) = TQTET(NCHI+1)

      K=1
      DO J=1,NCHI
         IF (ABS(TQTET(NCHI+J+1)-TQTET(NCHI+J)).GT.PI) K=J
      ENDDO

      IF (K.GT.1) THEN
         DO J=K+1,NCHI+1
            TQTET(NCHI+J) = TQTET(NCHI+J)+2.*PI
         ENDDO
      ENDIF

      DO J=1,NCHI
         TQTET(J) = TQTET(NCHI+J) - 2.*PI
         TQTET(2*NCHI+J+1) = TQTET(NCHI+J+1) + 2.*PI
      ENDDO

C     COMPUTE CHI(THETA) WITH UNIFORM THETA-MESH 
C     USING SPLINE
      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

C     SAVE TQCHIN INTO A FILE
      DO J=1,NCHI+1
         WRITE(CHOUTP,120) TQCHIN(J)
      ENDDO

      ENDDO

 120  FORMAT(E15.8)
      CLOSE(CHOUTP)

      DEALLOCATE(TQTET,TQCHI,TQTETN,TQCHIN,WORK)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C READ IN CHI_PEST(CHI_EQAC) FROM FILES                LIU YQ 19.12.2012
C AT ALL RATIONAL SURFACES, AND COMPUTE CONVERSION MATRICES FROM       $
C B1EQAC TO B1PEST                                                     $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE B1EQAC2B1PEST
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INTEGER    I,J,CHNEW,MSB
      PARAMETER  (CHNEW=121)
      REAL*8     RMSA
      COMPLEX*16 CTMP1,CTMP2
      REAL*8,DIMENSION(:),ALLOCATABLE::TQTET,TQCHI

      ALLOCATE(TQTET(NCHI+1),TQCHI(NCHI+1))
      IF (.NOT.ALLOCATED(TQB1MAT)) THEN 
         ALLOCATE(TQB1MAT(MSMAX,NRATSURF))
         TQB1MAT = 0.
      ENDIF

C     OPEN TWO FILES TO READ, AND GO THROUGH ALL RATIONAL SURFACES
      OPEN(CHOUTP,FILE='PANGLE_PEST.IN',STATUS='OLD')
      REWIND(CHOUTP)
      READ(CHOUTP,*) I       
      READ(CHOUTP,*) I     

      OPEN(CHNEW,FILE='PANGLE_EQAC.IN',STATUS='OLD')
      REWIND(CHNEW)
      READ(CHNEW,*) I       
      READ(CHNEW,*) I   
     
      TQB1MAT = 0.
      DO I=1,NRATSURF

C     READ IN CHI_EQAC=TQTET AND CHI_PEST=TQCHI
      DO J=1,NCHI+1
         READ(CHNEW,*) TQTET(J)
      ENDDO

      DO J=1,NCHI+1
         READ(CHOUTP,*) TQCHI(J)
      ENDDO

C     FOURIER DECOMPOSITION
      RMSA = -QPLS(IRATSURF(I))*RNTOR
      DO MSB=1,MSMAX
         IF (ABS(RMSA-RM(MSB,2)).LT.0.4) RMSA = RM(MSB,2)
      ENDDO
C     RMSA = RMSA + 1
      DO MSB=1,MSMAX
         DO J=1,NCHI
            CTMP1 = EXP(CI*RM(MSB,2)*TQTET(J)-
     &                  CI*RMSA*TQCHI(J))
            CTMP2 = EXP(CI*RM(MSB,2)*TQTET(J+1)-
     &                  CI*RMSA*TQCHI(J+1))
            TQB1MAT(MSB,I) = TQB1MAT(MSB,I) +
     &      (CTMP1+CTMP2)*(TQTET(J+1)-TQTET(J))
         ENDDO
         TQB1MAT(MSB,I) = TQB1MAT(MSB,I)*.25/PI
      ENDDO

      ENDDO
      CLOSE(CHOUTP)
      CLOSE(CHNEW)

      DEALLOCATE(TQTET,TQCHI)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--UPDATE EQUILIBRIUM ROTATION AND DENSITY PROFILES     LIU YQ 12.2.2016
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE QLIN_UPDATE
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      IMPLICIT NONE

      ROT   = ROTEQ   + TROTI
      ROTM  = ROTEQM  + TROTM
      DROT  = DROTEQ  + TDROTI
      DROTM = DROTEQM + TDROTM
      IF (NDENEQ.EQ.1) THEN
         RHO   = RHOEQ   + TRHOI
         RHOM  = RHOEQM  + TRHOM
      ENDIF

      RETURN
      END

