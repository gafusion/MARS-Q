C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C  NEW VERSION OF FEEDBACK FROM JULY 12, 1999. 
C  USING RESISTANCE INSTEAD OF SELF-INDUCTANCE.
C  USING TAUW AS NORMALIZATION FACTOR FOR PROPORTIONAL GAIN.
C-----------------------------------------------------------------------
C  NEW FEATURES FOR MARSCRG
C  ------------------------
C  ADDED NWALL - NUMBER OF WALLS IN THE VACUUM - TO THE NAMELIST
C  ADDED IWALLJ AND TAUWJ TO <GLOBAL.INC>
C  REPLACED OLD IWALL BY IWALLJ AND TAUW BY TAUWJ INMARSCRG.F AND FEEDBACK.F
C  NOW TO EXCLUDE WALL ONE NEED TO SET NWALL = 0!
C-----------------------------------------------------------------------
C  CHANGES, LYQ, MAY 22, 2002
C  --------------------------
C  ADD <OFEEDI, OSENSP, OSENSR> IN COMMOM BLOCK, FOR OUTPUT INFORMATION
C  OFEEDI --- FEEDBACK CURRENT
C  OSENSR --- RADIAL SENSOR SIGNAL AT POSITION ISENS(1)
C  OSENSP --- POLOIDAL SENSOR SIGNAL AT POSITION ISENS(1)-1
C----------------------------------------------------------------------- 
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C----------SET FEEDBACK DEFAULT VALUES--------Y.Q.LIU 13/04/1999--------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FEEDINI
C     ==================
      USE DIMENSIM
      USE FEEDBACKM
      USE GLOBALM
  
C     --- DEFAULT NAMELIST VALUES
      INCFEED = 0
      NCOIL   = 1
      NSENS   = 0
      NCOILT  = 6
      NSENST  = 6
      IFEED   = 0
      KCOILCURR = 1
      RL0     = 1.0
      IPSIF   = 1
      PHISIGN   =-1.0
      PHIPHASE  = 0.0
      DTRAMPUI  =-1.0
      FEEDI   = (0.,0.)
      FEEDIT  = (0.,0.)
      GAINA   = 0.
      GAINP   = 0.
      HGAINA  = 0.
      HGAINP  = 0.
      RMFS    = 1.0
      RLK     = 0.
      FCCHI   = 0.
      FWCHI   = 0.
      SCCHI   = 0.
      SWCHI   = 0.
      BTCHI   = 0.06
      OFEEDI  = 0.
      OSENSR  = 0.
      OSENSP  = 0.
      EXTCURR = 0.
      FNUP    = (0.,0.)
      FNUD    = (0.,0.)
      THRESHOLD = 1.0
      TAUF    = 2.E+4
      IDYNAM  = 0
      KKF     =-1
      KREADECA= 0
      IBOUT   = 2
      IBOUT2  = 3
      RAN_STD =-1.0
      VFMAX   =-1.0
      FINIC   = 1.0

C     --- THE SENSOR COIL WITH NO. J=0  LOCATES AT IFEED
      ISENS    = 0
      ISENS(0) = IFEED

C     --- WHEN NCASE = 10,20
C     --- FOR THE FIRST RUN (ISWEEP=1), ISWITCH = 0
      ISWITCH = 0

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------CHECK VALIDITY OF NAMELIST 'FEEDBACK'-----Y.Q.LIU 14/04/1999----
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FEEDCTRL
C     ===================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE FEEDBACKM
      USE RAN_MODM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'

      INTEGER    J, K
      REAL*8     SUM
      COMPLEX*16 CTMP

      IF (INCFEED.LT.0) THEN
         WRITE(*,90) INCFEED
         STOP 'INCFEED'
      ENDIF
 90   FORMAT('NO SUCH OPTION FOR FEEDBACK, INCFEED=',I2)

      IF (INCFEED.EQ.0) GOTO 1090

      IF (NCOIL.LT.0)   NCOIL   = 0
      IF (NCOILT.LT.0)  NCOILT  = 0
      IF (NSENS.LT.0)   NSENS   = 0
      IF (KKF.LT.-5)    KKF     = 1

      IF (INCFEED.EQ.6.AND.NWALL.EQ.0) THEN
         STOP 'NWALL MUST BE > 0 IF INCFEED = 6'
      ENDIF

      IF (NCOIL.GT.NCOIL0) THEN
         STOP 'NCOIL > NCOIL0'
      ENDIF
      IF (NCOIL.GT.MSMAX) THEN
         NCOIL = MSMAX
         WRITE(*,*) ' *** WARNING *** NCOIL SET TO MSMAX=',NCOIL
      ENDIF

      IF (NCOILT.GT.NCOILT0) THEN
         STOP 'NCOILT > NCOILT0'
      ENDIF

C     --- SET N-HARMONIC FOR ACTIVE COIL CURRENT AS DEFAULT
C     --- FOR BACKWARD COUPLING RUN 
C     --- '-' SIGN DUE TO OPPOSITE PHI-ANGLE IN CARIDDI
      CTMP = PHISIGN*CI*RNTOR*2.*PI/NCOILT
      DO K=1,NCOILT
         FEEDIT(K) = EXP((K-1)*CTMP+CI*RNTOR*PHIPHASE*PI)
      ENDDO

      IF (NSENS.GT.NSENS0) THEN
         NSENS = NSENS0
         WRITE(*,*) ' *** WARNING *** NSENS SET TO NSENS0'
      ENDIF

C     --- CHECK THE FEEDBACK AND SENSOR COILS LOCATION
      IF (INCFEED.NE.20.AND.INCFEED.NE.21) THEN
      IF (IFEED.LT.1.OR.IFEED.GE.NV) THEN
         WRITE(*,100) IFEED
         STOP 'IFEED'
      ENDIF
      IF (IABS(ISENS(0)-IFEED).GT.1) THEN
         WRITE(*,110) 0,ISENS(0)
         STOP 'ISENS(0)'
      ENDIF
      IF (ISENS(0)-IFEED.EQ.1.AND.(INCFEED.EQ.3.OR.INCFEED.EQ.5)) THEN
         WRITE(*,110) 0,ISENS(0)
         STOP 'ISENS(0)'
      ENDIF
      DO J=1,NSENS
         IF (ISENS(J).GE.NV) THEN
            WRITE(*,110) J,ISENS(J)
C           STOP 'ISENS(J)'
         ENDIF
      ENDDO
 100  FORMAT(' INCORRECT LOCATION OF FEEDBACK COIL, IFEED=',I3)
 110  FORMAT(' INCORRECT LOCATION OF SENSOR COIL NO. J=',I2,
     *       ' ISENS(J)=',I3)
      ENDIF

C     --- CHECK FEEDBACK COIL WIDTH FROM NAMELIST <FEEDBACK>
      DO K=1,NCOIL
         IF (FWCHI(K).LT.BTCHI) THEN
            WRITE(*,130) K,FWCHI(K),BTCHI
C           STOP 'FWCHI(K) < BTCHI'
         ENDIF
      ENDDO
 130  FORMAT(' K=',I2,' FWCHI(K)=',F6.3,' ',F6.3)

C     --- CHECK SENSOR COIL WIDTH FROM NAMELIST <FEEDBACK>
      DO K=1,NCOIL
         IF (SWCHI(K).LT.0.0) THEN
            WRITE(*,140) K,SWCHI(K)
            STOP 'SWCHI(K) < 0'
         ENDIF
      ENDDO
 140  FORMAT(' K=',I2,' SWCHI(K)=',F6.3,/,
     *       ' SENSOR COIL WIDTH SHOULD NOT BE NEGATIVE')

C     --- CHECK <IDYNAM>:
      IF (IDYNAM.LT.0.AND.IDYNAM.GT.4) THEN
          STOP '<IDYNAM> SHOULD BE BETWEEN 0 AND 4'
      ENDIF
C     IF (IDYNAM.EQ.3) FEEDI=(0.,0.)

      IF (INCFEED.GT.0.AND.VFMAX.GT.0.0) THEN
         IF (MSMAX.LT.NCOIL) STOP 'MSMAX < NCOIL'
      ENDIF

C     --- CHECK PROPORTIONAL AND DERIVATIVE GAINS
      DO K=1,NCOIL
         GAIN(K) = GAINA(K)*EXP(GAINP(K)*PI*CI)
      ENDDO

C     SPECIAL TREATMENT: SET T=CONST
      IF (INCFEED.EQ.2.OR.INCFEED.EQ.4.OR.INCFEED.EQ.6.OR.
     *    INCFEED.EQ.10.OR.INCFEED.EQ.9.OR.INCFEED.EQ.12.OR.
     *    INCFEED.EQ.18.OR.INCFEED.EQ.22) THEN
      DO J=1,NRP1
         T(J)  = 1.0    
      ENDDO
      DO J=1,NR
         TM(J) = 1.0    
      ENDDO
      ENDIF

 1090 CONTINUE

C     --- INITIALIZATION ARRAYS AND MATRICES
      IF (.NOT. ALLOCATED (AFEEDS)) THEN
         ALLOCATE(  AFEEDS(MSDIM,MSDIM),  AFEEDF(MSDIM,MSDIM))
         ALLOCATE(  CFEED(MSDIM,NCOIL0),  BFEEDF(NCOIL0,MSDIM))
         ALLOCATE(BFEEDSS(NCOIL0,MSDIM),  BFEEDST(NCOIL0,MSDIM))
         ALLOCATE( HCFEED(MSDIM,NCOIL0))
         ALLOCATE( CFEED2(MSDIM,NCOIL0),  HCFEED2(MSDIM,NCOIL0))
      END IF
 
      CFEED   = (0.,0.)
      HCFEED  = (0.,0.)
      CFEED2  = (0.,0.)
      HCFEED2 = (0.,0.)
      BFEEDF  = (0.,0.)
      BFEEDSS = (0.,0.)
      BFEEDST = (0.,0.)
      DFEED   = (0.,0.)
      AFEEDF  = (0.,0.)
      AFEEDS  = (0.,0.)

C     GENERATE ARRAY OF SENSOR NOISE VS. TIME (ITERATION)
C     AS COMPLEX QUANTITY
      IF (.NOT.ALLOCATED(SENSNOISE)) ALLOCATE(SENSNOISE(NITMAX+1))
      SENSNOISE = (0.0,0.0)
      IF (RAN_STD.GT.0.) THEN
         DO J=1,NITMAX+1
            SENSNOISE(J) = RAN_NORMAL(0.0,RAN_STD)
C    &                     *EXP(CI*RAN_NORMAL(0.0,PI))
         ENDDO
      ENDIF
      OPEN(CHOUTP,FILE='SENSNOISE.OUT')
      REWIND(CHOUTP)
      DO J=1,NITMAX+1
         WRITE(CHOUTP,111) REAL(SENSNOISE(J)),IMAG(SENSNOISE(J))
      ENDDO
      CLOSE(CHOUTP)
 111  FORMAT(2(1X,E15.8))

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------CALCULATE ACTIVE COIL'S RESISTANCE ----Y.Q.LIU 13/07/1999-------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETRFK
C     ===================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE FEEDBACKM

      INTEGER K
      REAL*8  CHIA, CHIB
      
      DO K=1,NCOIL
CC         CHIA = FCCHI(K) - 0.5*FWCHI(K) - 0.5*BTCHI
CC         CHIB = FCCHI(K) + 0.5*FWCHI(K) + 0.5*BTCHI
CC         RLK(K) = RL0*(ASPCT*2 + VCS(IFEED)*(CHIB-CHIA+
CC     *          DCOS(CHIA*PI)+DCOS(CHIB*PI)))/(ASPCT+1.0)/2.0
         RLK(K) = RL0 
      ENDDO

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------CALCULATE MATRIX 'C' OF FOURIER COEFF.-----Y.Q.LIU 14/04/1999---
C KKF=-4: ALLOW COMBINATION OF ESC CURRENT AND COIL CURRENT IN THE FORM
C         J_TOT=J_ESC + J_COIL*FEEDI
C         BOTH CURRENTS NEED TO BE SPECIFIED AT THE SURFACE SURFACE
C         SPECIFIED BY IFEED. THE J_ESC IS REPRESNTED BY HCFEED2 AND
C         CFEED2, AND THE J_COIL IS REPRESNTED BY HCFEED AND CFEED.
C         FOR OTHER OPTIONS OF KKF>-4, THE CURRENT J_TOT IS ALWAYS
C         REPRESENTED BY HCFEED AND CFEED.
C KKF=-5: SIMILAR TO KKF=-4, BUT J_COIL IS NOW ALSO REPRESENTED BY ESC
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CALCFEED
C     ===================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE FEEDBACKM

      INTEGER MS, K, NFC, J, KKF2, NFCK
      PARAMETER ( NFC=1025 )
      REAl*8  DCHI,  DCHI2, CHIK, BDET, CHIF, AFK0, 
     &        CHI1K, CHI2K, FEEDK(NFC),GGK(NFC)
      COMPLEX*16 RMI,   TMP,   RMD,       DHH
      INTEGER FILEID,HARMONUM
      REAL*8 AMP_REAL,AMP_IMAG
     
      DHH = CI/PI/(VCS(IFEED+1) - VCS(IFEED-1))

C     --- KKF=0: USING HEVISAID-FUNCTION FOR J3U REPRESENTATION
      IF (KKF.EQ.0) THEN
      DO MS=1,MSMAX
         IF (IM(MS,2).EQ.0) THEN
            DO K=1,NCOIL
               HCFEED(MS,K) = DHH*FWCHI(K)*PI
            ENDDO
         ELSE
            RMI = -CI*RM(MS,2)*PI
            DO K=1,NCOIL
               RMD = CI*RM(MS,2)*BTCHI*PI
               CHI1K = FCCHI(K) - 0.5*FWCHI(K) - 0.5*BTCHI
               CHI2K = FCCHI(K) + 0.5*FWCHI(K) + 0.5*BTCHI
               TMP = EXP(RMI*CHI1K)*(EXP(-RMD)-1.) +
     &              EXP(RMI*CHI2K)*(EXP(RMD)-1.)
               HCFEED(MS,K) = DHH*TMP/RMI/RMD*PI
            ENDDO
         ENDIF
      ENDDO
      ENDIF

C     --- KKF=-1 OR -3: USING GAUSS-FUNCTION FOR J3U REPRESENTATION,
C     J2U(M) CALCULATED ANALITICALLY BY INTEGRATING OVER [-INF INF].
      IF (KKF.EQ.-1.OR.KKF.EQ.-3.OR.KKF.EQ.-4) THEN
      DO K=1,NCOIL
         BDET  = BTCHI
         CHI1K = (FCCHI(K)+0.5*FWCHI(K))*PI
         CHI2K = (FCCHI(K)-0.5*FWCHI(K))*PI
         DO MS=1,MSDIM
            IF (IM(MS,2).EQ.0) THEN
               HCFEED(MS,K) = DHH*(CHI1K-CHI2K)
            ELSE
               RMI = -CI*RM(MS,2)
               TMP = (EXP(RMI*CHI1K)-EXP(RMI*CHI2K))/RM(MS,2)
               HCFEED(MS,K) = DHH*CI*TMP*EXP((RMI*BDET)**2*PI/4)
            ENDIF
         ENDDO
      ENDDO   
      ENDIF

C     --- KKF>0: USING GAUSS-FUNCTION OF ORDER OF KKF FOR J3U 
C     REPRESENTATION, J2U(M) CALCULATED NUMERICALLY BY INTEGRATING 
C     OVER [-PI PI].
      IF (KKF.GT.0) THEN
      KKF2 = KKF + KKF
      DO K=1,NCOIL
         BDET  = BTCHI
         CHIF  = -1.0
         DCHI  = 2.0/(NFC-1)
         DCHI2 = DCHI/2.0*PI
         CHI1K = FCCHI(K) - 0.5*FWCHI(K) - 0.5*BTCHI
         CHI2K = FCCHI(K) + 0.5*FWCHI(K) + 0.5*BTCHI
         DO J=1,NFC
            FEEDK(J) = EXP(-PI*(((CHIF-CHI2K)/BDET)+0.5)**KKF2)
            FEEDK(J) = FEEDK(J)-
     &                 EXP(-PI*(((CHIF-CHI1K)/BDET)-0.5)**KKF2)
            CHIF     = CHIF + DCHI
         ENDDO

         GGK(1) = 0.0
         DO J=1,NFC-1
            GGK(J+1) = GGK(J) + (FEEDK(J)+FEEDK(J+1))*DCHI2
         ENDDO

         WRITE(*,*) 'GGK(NFC) =',GGK(NFC)

         NFCK = (FCCHI(K)+1)/DCHI + 1
         AFK0 = -GGK(NFCK)
         DO MS=1,MSMAX
            HCFEED(MS,K) = (0.0,0.0)
            RMI  = -CI*RM(MS,2)*PI
            CHIF = -1.0 + DCHI
            DO J=2,NFC-1
               HCFEED(MS,K) = HCFEED(MS,K) + GGK(J)*EXP(RMI*CHIF)
               CHIF = CHIF + DCHI
            ENDDO
            HCFEED(MS,K) = -HCFEED(MS,K)*DHH*DCHI*PI/AFK0
         ENDDO
      ENDDO
      ENDIF
      
      IF (KKF.EQ.-5) THEN
         FILEID = 401
         OPEN(FILEID,FILE='CURHARMO_COIL.IN',STATUS='OLD')
         REWIND(FILEID)
         READ (FILEID,*) HARMONUM
         DO J=1,HARMONUM
            READ (FILEID,*) K, AMP_REAL,AMP_IMAG
            DO MS = 1,MSMAX
                IF (RM(MS,2) .EQ. K) THEN
                    HCFEED(MS,1) = CMPLX(AMP_REAL,AMP_IMAG)/RNTOR
                ENDIF
            ENDDO
         ENDDO
         CLOSE (FILEID)
      ENDIF

      IF ((KKF.EQ.-2.OR.KKF.EQ.-4.OR.KKF.EQ.-5).AND.
     &    ABS(RNTOR).GT.1.0E-10) THEN
         FILEID = 401
         OPEN(FILEID,FILE='CURHARMO.IN',STATUS='OLD')
         REWIND(FILEID)
         READ (FILEID,*) HARMONUM
         DO J=1,HARMONUM
            READ (FILEID,*) K, AMP_REAL,AMP_IMAG
            DO MS = 1,MSMAX
                IF (RM(MS,2) .EQ. K) THEN
                    HCFEED2(MS,1) = CMPLX(AMP_REAL,AMP_IMAG)/RNTOR
                ENDIF
            ENDDO
         ENDDO
         CLOSE (FILEID)

         IF (KKF.EQ.-2) HCFEED = HCFEED2
      ENDIF

C     DISTINGUISH BETWEEN N.NE.0 RWM CONTROL AND VERTICAL STABILITY CONTROL
      IF (ABS(RNTOR).LT.1.0E-10) THEN
         DO MS=1,MSMAX
            DO K=1,NCOIL
               CFEED(MS,K) =-HCFEED(MS,K)*RM(MS,2)
            ENDDO
         ENDDO
      ELSE
         CFEED = HCFEED*RNTOR    
         IF (KKF.EQ.-4.OR.KKF.EQ.-5) CFEED2 = HCFEED2*RNTOR
      ENDIF

C     INCLUDE COIL CURRENT PHASE
C     YQLIU, 2009-10-18
C     NEED TO INTRODUCE A NEW NAMELIST VARIABLE FOR THE CURRENT PHASE!
C     DO MS=1,MSMAX
C        DO K=1,NCOIL
C           CFEED(MS,K) = CFEED(MS,K)*FEEDI(K)
C        ENDDO
C     ENDDO

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------CALCULATE MATRIX 'E' OF FOURIER COEFF.-----Y.Q.LIU 14/04/1999---
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     Jmn = CFEED * FEEDI, m=1,MSMAX
C              DFEED
C     CFEED = [     ]
C              FFEED
C     FEEDI = DFEED^-1 * Jmn, m=1,NCOIL
C-----------------------------------------------------------------------    
      SUBROUTINE CALDFEED
C     ===================
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      INCLUDE 'compam.inc'

      INTEGER MS, K, M, KL
      COMPLEX*16 SCR(NCOIL0*(2*NCOIL0+1)), TMP(NCOIL0), DET

      DO M=1,NCOIL0*NCOIL0
         DFEED(M) = (0.,0.)
      ENDDO

      DO K=1,NCOIL
         KL = (K-1)*NCOIL
         DO MS=1,NCOIL
            DO M=1,MSMAX
               DFEED(MS+KL) = DFEED(MS+KL) + CFEED(M,MS)*CFEED(M,K)
            ENDDO
         ENDDO
      ENDDO

      CALL FMIND(DFEED,TMP,NCOIL,NCOIL,SCR,DET,EPSDET,0,1,IFEED)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------CALCULATE MATRIX 'A' FOR VERSION 3-----Y.Q.LIU 21/04/1999-------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     PSI   = BFEED * b^s
C     AFEED = CFEED * BFEED
C-----------------------------------------------------------------------    
      SUBROUTINE CALAFEED
C     ===================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE FEEDBACKM
      INCLUDE 'compam.inc'

      INTEGER MS, K, M, J, NFC, MS1
      PARAMETER ( NFC=1001 )
      REAL*8    DCHI, CHIS, DTK, CHIAB2, SFACT, CHIA, CHIB 
      COMPLEX*16 RMI, TMP, RMD

      IF (IPSIF.EQ.0) THEN
      DO MS=1,MSMAX
         IF (IM(MS,2).EQ.0) THEN
            DO K=1,NCOIL
               BFEEDF(K,MS) = FWCHI(K)*PI
            ENDDO
         ELSE
            RMI = CI*RM(MS,2)*PI
            DO K=1,NCOIL
               CHIA = FCCHI(K) - 0.5*FWCHI(K) - 0.5*BTCHI
               CHIB = FCCHI(K) + 0.5*FWCHI(K) + 0.5*BTCHI
               RMD = CI*RM(MS,2)*BTCHI*PI
               TMP = EXP(RMI*CHIA)*(EXP(RMD)-1.) +
     *              EXP(RMI*CHIB)*(EXP(-RMD)-1.)
               BFEEDF(K,MS) = -TMP/RMI/RMD*PI
            ENDDO
         ENDIF
      ENDDO
      ELSE
      TMP = -(VCS(IFEED+1)-VCS(IFEED-1))*PI*CI
      DO MS=1,MSMAX
         DO K=1,NCOIL
            BFEEDF(K,MS) = 0.0
         ENDDO
         MS1 = 2-MS-M1-M1
         IF (MS1.GE.1.AND.MS1.LE.MSMAX) THEN
            DO K=1,NCOIL
               BFEEDF(K,MS) = TMP*HCFEED(MS1,K)
            ENDDO
         ENDIF
      ENDDO      
      ENDIF

      DO M=1,MSMAX
         DO MS=1,MSMAX
            AFEEDF(MS,M) = (0.,0.)
            DO K=1,NCOIL
               AFEEDF(MS,M) = AFEEDF(MS,M) + 
     *                        CFEED(MS,K)*BFEEDF(K,M)/RLK(K)
            ENDDO
         ENDDO
      ENDDO

      DO MS=1,MSMAX
         IF (IM(MS,2).EQ.0) THEN
            DO K=1,NCOIL
               BFEEDSS(K,MS) = SWCHI(K)*PI
            ENDDO
         ELSE
            RMI = CI*RM(MS,2)*PI
            DO K=1,NCOIL
               TMP=EXP(RMI*(SCCHI(K)+0.5*SWCHI(K)))-
     *             EXP(RMI*(SCCHI(K)-0.5*SWCHI(K)))
               BFEEDSS(K,MS) = TMP/RMI*PI
            ENDDO
         ENDIF
      ENDDO

      IF (INCFEED.EQ.5) THEN
         DO K=1,NCOIL
            DO MS=1,MSMAX
               BFEEDSS(K,MS) = CI*BFEEDSS(K,MS)
            ENDDO
         ENDDO
      ENDIF

      IF (INCFEED.EQ.11.OR.INCFEED.EQ.12) THEN
      DO K=1,NCOIL
         CHIAB2 = SCCHI(K)*PI
         DO MS=1,MSMAX
            BFEEDSS(K,MS) = EXP(CHIAB2*RM(MS,2)*CI)
         ENDDO
      ENDDO
      ENDIF

      DO K=1,NCOIL
         CHIAB2 = SCCHI(K)*PI
         DO MS=1,MSMAX
            BFEEDST(K,MS) = CI*EXP(CHIAB2*RM(MS,2)*CI)
         ENDDO
      ENDDO

      IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.2
     *    .OR.INCFEED.EQ.5.OR.INCFEED.EQ.12) THEN
      DO M=1,MSMAX
         DO MS=1,MSMAX
            AFEEDS(MS,M) = (0.,0.)
            DO K=1,NCOIL
               AFEEDS(MS,M) = AFEEDS(MS,M) + 
     *                        GAIN(K)*CFEED(MS,K)*BFEEDSS(K,M)/RMFS(K)
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      IF (INCFEED.EQ.3.OR.INCFEED.EQ.4.OR.INCFEED.EQ.8.OR.
     &    INCFEED.EQ.10.OR.INCFEED.EQ.9.OR.INCFEED.EQ.18.OR.
     &    INCFEED.EQ.22) THEN
      DO M=1,MSMAX
         DO MS=1,MSMAX
            AFEEDS(MS,M) = (0.,0.)
            DO K=1,NCOIL
               AFEEDS(MS,M) = AFEEDS(MS,M) + 
     *                        GAIN(K)*CFEED(MS,K)*BFEEDST(K,M)/RMFS(K)
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---SET INITIAL VALUE OF J2U OR J3U AT POINT IFEED---Y.Q.LIU 15/04/1999-
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SETJINI
C     ===================
      USE DIMENSIM
      USE GLOBALM
      USE RESMATM
      USE FEEDBACKM

      INTEGER MS, K, I, LX
      REAL*8  RTMP

      I = NR + IFEED

      RTMP = 1.0
      IF (DTRAMPUI.GT.0..AND.KCOILCURR.EQ.2) RTMP = 0.0
      
      DO MS=1,MSMAX
         LX            = (MS-1)*NXCOMP + KXJFB
         X(LX,I) = (0.,0.)   
         DO K=1,NCOIL
            X(LX,I) = X(LX,I) + CFEED2(MS,K)+CFEED(MS,K)*FEEDI(K)*RTMP
         ENDDO
         IF (ABS(RNTOR).LT.1.0E-10) THEN
            J3U(I,MS) = X(LX,I)
         ELSE
            J2U(I,MS) = X(LX,I)
         ENDIF
      ENDDO

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------CONSTRUCT MATRIX DUE TO FEEDBACK-----Y.Q.LIU 14/04/1999---------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FEEDM(
C     =================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM) 
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE GIJLM
      USE FEEDBACKM
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'

      INTEGER MS, I, II, LX, M, LCOL, J, K, L, KXJ, ID, LY, KYBX,
     &        NSA, NSB, MSA, MSB, MSPL, MSMI, LYROW, LXCOL, LYCOL,
     &        NRTEMP, KY1
      COMPLEX*16 TMP, TMP1, CMA, CNA
      PARAMETER  (NSA=2,NSB=1)

      INCLUDE 'integc.inc'

      IF (INCFEED.LT.1) RETURN

      I   = NR + IFEED

C     --- INCFEED=1,11,3,5,9,12
      IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.3
     *    .OR.INCFEED.EQ.5.OR.INCFEED.EQ.9.OR.INCFEED.EQ.12) THEN
      K = I
      TMP = -1.0
      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP
         BSUBM(KXB1+LX,KXJFB+LX,K) = TMP
      ENDDO

C     USE V2U TO REPRESENT CONTROL COIL VOLTAGE/CURRENT
C     THIS REQUIRES MSMAX > NCOIL
      IF (VFMAX.GT.0.0.AND.NCASE.EQ.4) THEN
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP
            DO J=1,NCOIL
               LY = (J-1)*NYCOMP
               ESUBM(KXB1+LX,KYV2+LY,I) = CFEED(MS,J)
            ENDDO
         ENDDO

         DO J=1,NCOIL
            LY = (J-1)*NYCOMP
            DSUBM(KYV2+LY,KYV2+LY,I) = 1.0
         ENDDO
      ENDIF

      IF (IDYNAM.EQ.2) THEN
      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP
         BSUBM(KXB1+LX,KXJFB,K) = BSUBM(KXB1+LX,KXJFB,K) - 
     *                            TMP*CFEED(MS,1)/CFEED(1,1)
      ENDDO
      ENDIF

      IF (IDYNAM.GT.0.AND.IDYNAM.LT.3) THEN
      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP
         DO M=1,MSMAX
            LCOL = (M-1)*NXCOMP
            IF (ABS(RM(M,2)).GT.0.1) THEN
            BSUBM(KXB1+LX,KXB1+LCOL,K) = -TAUF*AL0*AFEEDF(MS,M)
            ELSE
            BSUBM(KXB1+LX,KXB1+LCOL,K) = -TAUF*AL0*AFEEDF(MS,M)*T(NRP1)
            ENDIF
         ENDDO
      ENDDO
      ENDIF

      ENDIF

C     THE FOLLOWING SENSOR BLOCK IN ACTION ONLY IF VFMAX < 0
      IF (VFMAX.LT.0.0) THEN

C     --- SENSOR COILS: INCFEED=1,11,12
      IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.12) THEN
C     --- TERMS FROM SENSOR COILS WITH J=0
      TMP = -(FNUP(0) + AL0*FNUD(0)*TAUF)
      K = I
      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP + KXB1
         DO M=1,MSMAX
            LCOL = (M-1)*NXCOMP + KXB1
            TMP1 = TMP*AFEEDS(MS,M)
            IF (ABS(RM(M,2)).LT.0.1) TMP1 = TMP1*T(NRP1)
            IF (ISENS(0)-IFEED.EQ.-1) ASUBM(LX,LCOL,K) = TMP1
            IF (ISENS(0)-IFEED.EQ.1)  CSUBM(LX,LCOL,K) = TMP1
            IF (ISENS(0)-IFEED.EQ.0)  BSUBM(LX,LCOL,K) = TMP1 +
     *                                BSUBM(LX,LCOL,K)
         ENDDO
      ENDDO

C     --- TERMS FROM SENSOR COILS WITH J=1
      IF (NSENS.GE.1) THEN
         TMP = -(FNUP(1) + AL0*FNUD(1)*TAUF)
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP + KXB1
            DO M=1,MSMAX
               LCOL = (M-1)*NXCOMP + KXV1
               TMP1 = TMP*AFEEDS(MS,M)
               BSUBM(LX,LCOL,K) = TMP1
            ENDDO
         ENDDO
      ENDIF

C     --- USE V1 TO TRANSFER THE VALUE OF B1 FROM ISENS TO IFEED
      IF (NSENS.GE.1.AND.ISENS(1).GT.1) THEN
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP
            J  = ISENS(1) + NR 
            K = J
            BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
            BSUBM(KXV1+LX,KXB1+LX,K) = -1.0
            DO J=2,ISENS(1)-1
               K = J+NR
               BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
               CSUBM(KXV1+LX,KXV1+LX,K) = -1.0
            ENDDO   
            DO J=ISENS(1)+1,NV-1
               K = J+NR
               BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
               ASUBM(KXV1+LX,KXV1+LX,K) = -1.0
            ENDDO
         ENDDO   
      ENDIF
      ENDIF

C     --- SENSOR COILS: INCFEED=3,5,9
      IF (INCFEED.EQ.3.OR.INCFEED.EQ.5.OR.INCFEED.EQ.9) THEN
      K = I
      IF (INCFEED.EQ.3.OR.INCFEED.EQ.9) KYBX = KYB2   
      IF (INCFEED.EQ.5) KYBX = KYB3   
C     --- TERMS FROM SENSOR COILS WITH J=0
      TMP = -(FNUP(0) + AL0*FNUD(0)*TAUF)
      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP + KXB1
         DO M=1,MSMAX
            LCOL = (M-1)*NYCOMP + KYBX
            TMP1 = TMP*AFEEDS(MS,M)
            IF (ISENS(0)-IFEED.EQ.-1) HSUBM(LX,LCOL,K) = TMP1
            IF (ISENS(0)-IFEED.EQ.0)  ESUBM(LX,LCOL,K) = TMP1
         ENDDO
      ENDDO

C     --- TERMS FROM SENSOR COILS WITH J=1
      IF (NSENS.GE.1) THEN
         TMP = -(FNUP(1) + AL0*FNUD(1)*TAUF)
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP + KXB1
            DO M=1,MSMAX
               LCOL = (M-1)*NXCOMP + KXV1
               TMP1 = TMP*AFEEDS(MS,M)
               BSUBM(LX,LCOL,K) = TMP1
            ENDDO
         ENDDO
      ENDIF

C     --- USE V1 TO TRANSFER THE VALUE OF B2 FROM ISENS TO IFEED
      IF (NSENS.GE.1) THEN
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP
            LY = (MS-1)*NYCOMP
            J  = ISENS(1) + NR
            K = J
            BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
            ESUBM(KXV1+LX,KYBX+LY,K) = -1.0
            DO J=2,ISENS(1)-1
               K = J+NR
               BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
               CSUBM(KXV1+LX,KXV1+LX,K) = -1.0
            ENDDO   
            DO J=ISENS(1)+1,NV-1
               K = J+NR
               BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
               ASUBM(KXV1+LX,KXV1+LX,K) = -1.0
            ENDDO
         ENDDO   
      ENDIF
      ENDIF

      ENDIF
C     THE ABOVE SENSOR BLOCK IN ACTION ONLY IF VFMAX < 0

      K = I
C     --- INCFEED=8: 
C     ERROR FIELD AMPLIFICATION
      IF (INCFEED.EQ.8) THEN
      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP
         BSUBM(KXB1+LX,KXJFB+LX,K) = 1.
      ENDDO
      ENDIF

 1011 CONTINUE

      IF (INCFEED.EQ.8.AND.IDYNAM.EQ.3) THEN
      DO K=1,NRATSURF
         J  = IRATSURF(K)
         MS = NINT(-RNTOR*Q(J))
         IF (M1.LE.MS.AND.MS.LE.M2) THEN
            II = MS-M1+1
            CALL ANNIHX(KXV1,II,J,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            LX = (II-1)*NXCOMP + KXV1
            BSUBM(LX,LX,J) = 1.
         ENDIF
      ENDDO
      ENDIF

      IF (INCFEED.EQ.8.AND.IDYNAM.EQ.4) THEN
      DO K=1,NRATSURF
         L  = IRATSURF(K)
         MS = NINT(-RNTOR*Q(L))
         DO MS=1,MSMAX
         IF (M1.LE.MS.AND.MS.LE.M2.AND.L.GT.1.AND.L.LT.NRP1) THEN
            II = MS-M1+1
            LX = (II-1)*NXCOMP + KXV1
            DO J=L-1,L+1
            CALL ANNIHX(KXV1,II,J,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            ASUBM(LX,LX,J) =-1./CSH(J-1)
            BSUBM(LX,LX,J) = 1./CSH(J-1) + 1./CSH(J)
            CSUBM(LX,LX,J) =-1./CSH(J)
            ENDDO
         ENDIF
         ENDDO
      ENDDO
      ENDIF

C     --- INCFEED=2,4,6,10,9,12,18,20,22: 
C     SELF OR MUTUAL INDUCTANCE CALCULATION
      IF (INCFEED.EQ.2.OR.INCFEED.EQ.4.OR.INCFEED.EQ.6.OR.
     *    INCFEED.EQ.10.OR.INCFEED.EQ.9.OR.INCFEED.EQ.12.OR.
     *    INCFEED.EQ.18.OR.INCFEED.EQ.20.OR.INCFEED.EQ.22) THEN
      K = I
      IF (INCFEED.NE.20) THEN
      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP
         BSUBM(KXB1+LX,KXJFB+LX,K) = 1.
      ENDDO
      ENDIF

      IF (INCFEED.EQ.6) THEN
         DO J=1,NWALL
            L = NR + IWALL(J)
            K = L
            DO MS=1,MSMAX
               CALL ANNIHX(KXB1,MS,K,
     *              ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            ENDDO
            DO MS=1,MSMAX
               LX = (MS-1)*NXCOMP
               IF (ABS(RM(MS,2)).LT.0.1.AND.ABS(RNTOR).LT.1.0E-13) 
     *           BSUBM(KXB1+LX,KXB1+LX,K) = 1.
               IF (ABS(RM(MS,2)).GE.0.1.AND.ABS(RNTOR).LT.1.0E-13) 
     *           BSUBM(KXB1+LX,KXJ3+LX,K) = 1.
               IF (ABS(RNTOR).GT.1.0E-13) 
     *           BSUBM(KXB1+LX,KXJFB+LX,K) = 1.
            ENDDO
         ENDDO
      ENDIF
      
      NRTEMP = NR
      IF (INCFEED.EQ.20.AND.NRATSURF.GE.KRATSURF) 
     &   NRTEMP = IRATSURF(KRATSURF)-1
C     WHEN THE GRID REACHES THE AXIS, HALF-CELL 1 KEEPS THE BOUNDC ROWS
C     (J1=0, DIV B, OTHERS 0): REPLACING ITS B2 ROW BY J1=0 WOULD DUPLICATE
C     THE AXIS J1 ROW AND LEAVE B2, B3 WITH ONE EQUATION (SINGULAR BLOCK)
      KY1 = 1
      IF (SLEFT.LT.1.E-10) KY1 = 2
 
      DO MS=1,MSMAX
         DO K=1,NRTEMP+1
            CALL ANNIHX(KXB1,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            CALL ANNIHX(KXV1,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KXX1.GT.0) CALL ANNIHX(KXX1,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KXPD.GT.0) CALL ANNIHX(KXPD,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (ABS(RNTOR).GT.0.1) CALL ANNIHX(KXJ3,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KXW1.GT.0) CALL ANNIHX(KXW1,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KXDPHI.GT.0) CALL ANNIHX(KXDPHI,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
         ENDDO
         DO K=KY1,NRTEMP
            CALL ANNIHY(KYV2,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYV3.GT.0) CALL ANNIHY(KYV3,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            CALL ANNIHY(KYB2,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            CALL ANNIHY(KYB3,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            CALL ANNIHY(KYPR,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYPE.GT.0) CALL ANNIHY(KYPE,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYPPERP.GT.0) CALL ANNIHY(KYPPERP,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYPPARA.GT.0) CALL ANNIHY(KYPPARA,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYX2.GT.0) CALL ANNIHY(KYX2,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYX3.GT.0) CALL ANNIHY(KYX3,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYW2.GT.0) CALL ANNIHY(KYW2,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYW3.GT.0) CALL ANNIHY(KYW3,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            IF (KYRHOP.GT.0) CALL ANNIHY(KYRHOP,MS,K,
     *           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
         ENDDO
      ENDDO

      DO MS=1,MSMAX
         LX = (MS-1)*NXCOMP
         LY = (MS-1)*NYCOMP
         CMA = CI*RM(MS,2)
         CNA = CI*RN(2)
         DO K=2,NRTEMP+1
            IF (INCFEED.EQ.10.OR.INCFEED.EQ.20.OR. 
     &          (ABS(RM(MS,2)).LT.0.1.AND.ABS(RNTOR).LT.1.0E-13)) THEN
              BSUBM(KXB1+LX,KXB1+LX,K) = 1.0
            ELSEIF (ABS(RNTOR).LT.1.0E-13) THEN
              BSUBM(KXB1+LX,KXJ3+LX,K) = 1.0
            ELSE
              BSUBM(KXB1+LX,KXJFB+LX,K) = 1.0
            ENDIF  
            BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
            IF (KXX1.GT.0) BSUBM(KXX1+LX,KXX1+LX,K) = 1.0
            IF (KXPD.GT.0) BSUBM(KXPD+LX,KXPD+LX,K) = 1.0
            IF (KXW1.GT.0) BSUBM(KXW1+LX,KXW1+LX,K) = 1.0
            IF (KXDPHI.GT.0) BSUBM(KXDPHI+LX,KXDPHI+LX,K) = 1.0
            BSUBM(KXJ3+LX,KXJ3+LX,K) = 1.0
         ENDDO
         DO K=KY1,NRTEMP
            IF (INCFEED.EQ.10.OR.INCFEED.EQ.20) THEN
              DSUBM(KYB2+LY,KYB2+LY,K) = 1.0
            ELSE
              IF (ABS(RM(MS,2)).GT.0.1.OR.ABS(RNTOR).GT.1.0E-13)
     &          DSUBM(KYB2+LY,KYJ1+LY,K) = 1.0
            ENDIF  
            DSUBM(KYV2+LY,KYV2+LY,K) = 1.0
            IF (KYV3.GT.0)    DSUBM(KYV3+LY,KYV3+LY,K) = 1.0
            IF (KYX2.GT.0)    DSUBM(KYX2+LY,KYX2+LY,K) = 1.0
            IF (KYX3.GT.0)    DSUBM(KYX3+LY,KYX3+LY,K) = 1.0
            IF (KYW2.GT.0)    DSUBM(KYW2+LY,KYW2+LY,K) = 1.0
            IF (KYW3.GT.0)    DSUBM(KYW3+LY,KYW3+LY,K) = 1.0
            DSUBM(KYPR+LY,KYPR+LY,K) = 1.0
            IF (KYPE.GT.0)    DSUBM(KYPE+LY,KYPE+LY,K) = 1.0
            IF (KYPPERP.GT.0) DSUBM(KYPPERP+LY,KYPPERP+LY,K) = 1.0
            IF (KYPPARA.GT.0) DSUBM(KYPPARA+LY,KYPPARA+LY,K) = 1.0
            IF (KYRHOP.GT.0)  DSUBM(KYRHOP+LY,KYRHOP+LY,K)   = 1.0
            IF (ABS(RM(MS,2)).GT.0.1.OR.ABS(RNTOR).GT.1.0E-13) THEN
            HC = CSH(K)
            ZNORM = 1./HC
            PTRAP = PTRAPI
            IF (ABS(RM(MS,2)).GT.0.1) THEN
              FSUBM(KYB3+LY,KXB1+LX ,K)= -ZNORM
              GSUBM(KYB3+LY,KXB1+LX ,K)=  ZNORM
            ELSE
              FSUBM(KYB3+LY,KXB1+LX ,K)= -ZNORM*T(K)
              GSUBM(KYB3+LY,KXB1+LX ,K)=  ZNORM*T(K+1)
            ENDIF
            DSUBM(KYB3+LY,KYB2+LY ,K)= CMA * GG(C1,C1,C1)
            DSUBM(KYB3+LY,KYB3+LY ,K)= CNA * GG(C1,C1,C1)
            ENDIF
         ENDDO
         K = 1
         BSUBM(KXB1+LX,KXB1+LX,K) = 1.0
         BSUBM(KXV1+LX,KXV1+LX,K) = 1.0
         IF (KXX1.GT.0) BSUBM(KXX1+LX,KXX1+LX,K) = 1.0
         IF (KXPD.GT.0) BSUBM(KXPD+LX,KXPD+LX,K) = 1.0
         IF (KXW1.GT.0) BSUBM(KXW1+LX,KXW1+LX,K) = 1.0
         BSUBM(KXJ3+LX,KXJ3+LX,K) = 1.0
      ENDDO

C     SPECIAL CASE FOR RNTOR=0 AND MROW=0
      IF (ABS(RNTOR).LT.1.0E-13) THEN
      DO MS=1,MSMAX
         IF (ABS(RM(MS,2)).LT.0.1) THEN
            LX = (MS-1)*NXCOMP
            DO K=2,NRTEMP+1
               BSUBM(KXB1+LX,KXB1+LX,K) = 1.0
            ENDDO
         ENDIF
      ENDDO

      DO 140 MSA=1,MSMAX
      DO 130 MSB=1,MSMAX
         MSPL= MPLUS(MSA,NSA,MSB,NSB)
         MSMI=MMINUS(MSA,NSA,MSB,NSB)
         IF (MSPL.LT.1.OR.ABS(RM(MSPL,2)).GT.0.1) GOTO 80

         LYROW=(MSPL-1)*NYCOMP
         LXCOL=(MSA -1)*NXCOMP
         LYCOL=(MSA -1)*NYCOMP
           
         DO  I=KY1,NRTEMP
            INCLUDE 'tophat.inc'
            DSUBM(KYB3+LYROW,KYB3+LYCOL,I)=
     &           GG(DG33LM(I,MSB),DG33L(I,MSB),DG33L(I+1,MSB))

            IF (ABS(RM(MSA,2)).GT.0.1) THEN 
               FSUBM(KYB2+LYROW,KXB1+LXCOL,I)= 
     &              GF(DG12L(I  ,MSB),DG12LM(I,MSB))
               GSUBM(KYB2+LYROW,KXB1+LXCOL,I)= 
     &              GF(DG12L(I+1,MSB),DG12LM(I,MSB))
            ELSE
               FSUBM(KYB2+LYROW,KXB1+LXCOL,I)= 
     &              GF(DG12L(I  ,MSB)*T(I),DG12LM(I,MSB)*TM(I))
               GSUBM(KYB2+LYROW,KXB1+LXCOL,K)= 
     &              GF(DG12L(I+1,MSB)*T(I+1),DG12LM(I,MSB)*TM(I))
            ENDIF

            DSUBM(KYB2+LYROW,KYB2+LYCOL,I)= 
     &           GG(DG22LM(I,MSB),DG22L(I,MSB),DG22L(I+1,MSB))
         ENDDO

 80      CONTINUE
         IF (MSB.LT.2) GOTO 110
         IF (MSMI.LT.1.OR.ABS(RM(MSMI,2)).GT.0.1) GOTO 110

         LYROW=(MSMI-1)*NYCOMP
         LXCOL=(MSA -1)*NXCOMP
         LYCOL=(MSA -1)*NYCOMP

         DO  I=KY1,NRTEMP
            INCLUDE 'tophat.inc'
            DSUBM(KYB3+LYROW,KYB3+LYCOL,I)=CONJG(
     &           GG(DG33LM(I,MSB),DG33L(I,MSB),DG33L(I+1,MSB)))

            IF (ABS(RM(MSA,2)).GT.0.1) THEN 
               FSUBM(KYB2+LYROW,KXB1+LXCOL,I)=CONJG( 
     &              GF(DG12L(I  ,MSB),DG12LM(I,MSB)))
               GSUBM(KYB2+LYROW,KXB1+LXCOL,I)=CONJG( 
     &              GF(DG12L(I+1,MSB),DG12LM(I,MSB)))
            ELSE
               FSUBM(KYB2+LYROW,KXB1+LXCOL,I)=CONJG( 
     &              GF(DG12L(I  ,MSB)*T(I),DG12LM(I,MSB)*TM(I)))
               GSUBM(KYB2+LYROW,KXB1+LXCOL,K)=CONJG( 
     &              GF(DG12L(I+1,MSB)*T(I+1),DG12LM(I,MSB)*TM(I)))
            ENDIF

            DSUBM(KYB2+LYROW,KYB2+LYCOL,I)=CONJG( 
     &           GG(DG22LM(I,MSB),DG22L(I,MSB),DG22L(I+1,MSB)))
         ENDDO
 110     CONTINUE
 130  CONTINUE
 140  CONTINUE

      ENDIF
      ENDIF

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C----------GET 'DX' DUE TO FEEDBACK--------Y.Q.LIU 14/04/1999-----------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FEEDDX(
C     ==================
     &   MD, MDY, ND, R, RY, X, Y, XOLD, YOLD)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE GVACUUMM
      USE FEEDBACKM
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'

      INTEGER         MD, MDY, ND
      COMPLEX*16      R(MD,ND+1), X(MD,ND+1), RY(MDY,ND), Y(MDY,ND)
      COMPLEX*16      XOLD(MD,ND+1), YOLD(MDY,ND)
      COMPLEX*16      TMP,TTMP,CJ20,CJ21,CJ22,CJ30,CJ31,CJ32
      REAL*8          TMPR1,TMPI1,TMPR2,TMPI2,TMPR3,TMPI3,TMPR4,TMPR
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::TMPCC

      INTEGER MS, I, II, LX, J, K, LCOL, MM, KYBX

      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY


      IF (INCFEED.LT.1) RETURN

      I = NR + IFEED
      TTMP = -CALPHA1*CALPHA4/CALPHA2**2

      IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.3
     *    .OR.INCFEED.EQ.5.OR.INCFEED.EQ.9.OR.INCFEED.EQ.12) THEN
         DO MS=1,MSMAX
            LX      = (MS-1)*NXCOMP + KXB1
            R(LX,I) = (0.,0.)
         ENDDO
      
C     --- INCFEED=1,11,12
         IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.12) THEN
            DO MS=1,MSMAX
               LX      = (MS-1)*NXCOMP + KXB1
               IF (NCOUPL.EQ.9) THEN              
                  DO J=1,NCOIL
                     R(LX,I) = R(LX,I) - CFEED(MS,J)*FEEDI(J)
                  ENDDO
               ENDIF
               DO MM=1,MSMAX
                  LCOL    = (MM-1)*NXCOMP + KXB1
                  TMPR1 = 1.0
                  IF (ABS(RM(MM,2)).LT.0.1) TMPR1 = T(NRP1)
                  IF (IDYNAM.GT.0.AND.IDYNAM.LT.3) R(LX,I) = R(LX,I) + 
     *                          TTMP*TMPR1*TAUF*X(LCOL,I)*AFEEDF(MS,MM)
                  IF (VFMAX.LT.0.0) THEN
                  DO J=0,NSENS
                     R(LX,I) = R(LX,I)+TTMP*FNUD(J)*AFEEDS(MS,MM)*TAUF*
     *                    TMPR1*X(LCOL,NR+ISENS(J)) 
                  ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

C     --- INCFEED=3,5,9
         IF (INCFEED.EQ.3.OR.INCFEED.EQ.5.OR.INCFEED.EQ.9) THEN
            IF (INCFEED.EQ.3.OR.INCFEED.EQ.9) KYBX = KYB2   
            IF (INCFEED.EQ.5) KYBX = KYB3   
            DO MS=1,MSMAX
               LX      = (MS-1)*NXCOMP + KXB1
               DO MM=1,MSMAX
                  LCOL    = (MM-1)*NXCOMP + KXB1
                  TMPR1 = 1.0
                  IF (ABS(RM(MM,2)).LT.0.1) TMPR1 = T(NRP1)
                  IF (IDYNAM.GT.0.AND.IDYNAM.LT.3) R(LX,I) = R(LX,I) + 
     *                          TTMP*TMPR1*TAUF*X(LCOL,I)*AFEEDF(MS,MM)
                  LCOL    = (MM-1)*NYCOMP + KYBX
                  IF (VFMAX.LT.0.0) THEN
                  DO J=0,NSENS
                     R(LX,I) = R(LX,I) + TTMP*FNUD(J)*AFEEDS(MS,MM)
     *                         *TAUF*Y(LCOL,NR+ISENS(J))
                  ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDIF

C     CONTROL VOLTAGE/CURRENT SATURATION 
      IF (VFMAX.GT.0.0.AND.NCASE.EQ.4) THEN
      DO K=1,NCOIL
         LY = (K-1)*NYCOMP + KYV2
         TTMP = 0.0
         IF (NIT.GT.1) TTMP = RY(LY,I)
         RY(LY,I) = 0.0

         DO MS=1,MSMAX
         DO J=0,NSENS
            II = NR + ISENS(J)
            IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.12) THEN
               LX = (MS-1)*NXCOMP + KXB1
               RY(LY,I) = RY(LY,I) - GAIN(K)/RMFS(K)*BFEEDSS(K,MS)*
     &                    (FNUP(J)*X(LX,II)+FNUD(J)*TAUF*
     &                     (X(LX,II)-XOLD(LX,II))*REAL(CALPHA1))
            ENDIF

            IF (INCFEED.EQ.3.OR.INCFEED.EQ.5.OR.INCFEED.EQ.9) THEN
               IF (INCFEED.EQ.3.OR.INCFEED.EQ.9) KYBX = KYB2   
               IF (INCFEED.EQ.5) KYBX = KYB3   
               LX = (MS-1)*NYCOMP + KYBX
               RY(LY,I) = RY(LY,I) - GAIN(K)/RMFS(K)*BFEEDST(K,MS)*
     &                    (FNUP(J)*Y(LX,II)+FNUD(J)*TAUF*
     &                     (Y(LX,II)-YOLD(LX,II))*REAL(CALPHA1))
            ENDIF
         ENDDO
         ENDDO
         TMPR = 0.
         DO J=0,NSENS
            TMPR = TMPR + ABS(FNUD(J))
         ENDDO
         IF (TMPR.GT.0.) RY(LY,I) = 0.9*RY(LY,I) + 0.1*TTMP
         IF (ABS(RY(LY,I)).GT.VFMAX) THEN
            RY(LY,I) = RY(LY,I)/ABS(RY(LY,I))*VFMAX
         ENDIF
      ENDDO
      ENDIF

C     ADD SENSOR NOISE TO THE FEEDBACK LOOP
      IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.3
     *    .OR.INCFEED.EQ.5.OR.INCFEED.EQ.9.OR.INCFEED.EQ.12) THEN
         IF (ICALPHA2.EQ.0) TMPR = 1.
         IF (ICALPHA2.EQ.1) TMPR = 1./CALPHA2
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP + KXB1
            DO MM=1,MSMAX
               TMPR1 = 1.0
               IF (ABS(RM(MM,2)).LT.0.1) TMPR1 = T(NRP1)
               DO J=0,NSENS
                  R(LX,I) = R(LX,I)+TMPR*TMPR1*AFEEDS(MS,MM)*
     &            (FNUP(J)*SENSNOISE(NIT)+ FNUD(J)*TAUF*     
     &            (SENSNOISE(NIT+1)-SENSNOISE(NIT))*REAL(CALPHA1))
               ENDDO
            ENDDO
         ENDDO 
      ENDIF

C     --- INCFEED=20,21:          
      IF (INCFEED.EQ.20.OR.INCFEED.EQ.21.OR.INCFEED.EQ.22)  THEN
         DO MS=1,MXMAX
            DO J=1,NTP1
               R(MS,J) = 0.0
            ENDDO
         ENDDO
         DO MS=1,MYMAX
            DO J=1,NTOT
               RY(MS,J) = 0.0
            ENDDO
         ENDDO
      ENDIF

C     --- INCFEED=2,4,8,10,9,12: 
      IF (INCFEED.EQ.2.OR.INCFEED.EQ.4.OR.INCFEED.EQ.8.OR.
     &    INCFEED.EQ.10.OR.INCFEED.EQ.9.OR.INCFEED.EQ.12.OR.
     &    INCFEED.EQ.22) THEN
         IF (ICALPHA2.EQ.0) TMPR = 1.
         IF (ICALPHA2.EQ.1) TMPR = 1./CALPHA2

C        ADD SOURCE CURRENT TIME EVOLUTION      
         TMP = (1.0,0.0)
         IF (KCOILCURR.EQ.1) TMP = (1.0,0.0)
         IF (NCASE.EQ.5.OR.NCASE.EQ.9.OR.NCASE.EQ.10) 
     &      TOTTIME = TOTTIME+1.0/REAL(CALPHA1)
         IF (NCASE.EQ.6.OR.NCASE.EQ.5) THEN
            IF (KCOILCURR.EQ.2) THEN
               TMP = EXP(CI*IMAG(TALPHA1)*TOTTIME)
               IF (TOTTIME.GT.0..AND.TOTTIME.LT.DTRAMPUI) 
     &         TMP = TMP*TOTTIME/DTRAMPUI
            ENDIF

C           STANDING WAVE FOR THE COIL CURRENT
            IF (KCOILCURR.EQ.4) THEN
               TMPR4 = ABS(FEEDI(1))
               IF (TOTTIME.LE.COILCURR(1,1)) THEN
                  TMPR4 = COILCURR(1,2)
               ELSEIF (TOTTIME.GE.COILCURR(NCOILCURR,1)) THEN
                  TMPR4 = COILCURR(NCOILCURR,2)
               ELSE
                  CALL SPLINE1D(TMPR4,TOTTIME,1,COILCURR(:,2),
     &                          COILCURR(:,1),NCOILCURR,COILCURR2)
               ENDIF
C              TMP = TMPR4/ABS(FEEDI(1))
               TMP = TMPR4
            ENDIF

C           TRAVELLING WAVE FOR THE COIL CURRENT
            IF (KCOILCURR.EQ.5) THEN
               ALLOCATE(TMPCC(NCOIL))
               DO J=1,NCOIL
               IF (TOTTIME.LE.COILCURR(1,1)) THEN
                  TMPR1 = COILCURR(1,2*J)
                  TMPI1 = COILCURR(1,2*J+1)
               ELSEIF (TOTTIME.GE.COILCURR(NCOILCURR,1)) THEN
                  TMPR1 = COILCURR(NCOILCURR,2*J)
                  TMPI1 = COILCURR(NCOILCURR,2*J+1)
               ELSE
                  CALL SPLINE1D(TMPR1,TOTTIME,1,COILCURR(:,2*J),
     &                          COILCURR(:,1),NCOILCURR,COILCURR2)
                  CALL SPLINE1D(TMPI1,TOTTIME,1,COILCURR(:,2*J+1),
     &                          COILCURR(:,1),NCOILCURR,COILCURR2)
               ENDIF
               TMPCC(J) = TMPR1+CI*TMPI1
               ENDDO
            ENDIF
         ENDIF

         IF (INCFEED.EQ.4)  THEN
         DO MS=1,MXMAX
            DO J=1,NRP1
               R(MS,J) = 0.0
            ENDDO
         ENDDO
         DO MS=1,MYMAX
            DO J=1,NR
               RY(MS,J) = 0.0
            ENDDO
         ENDDO
         ENDIF

         DO MS=1,MSMAX
            LX      = (MS-1)*NXCOMP + KXB1
            R(LX,I) = (0.,0.)
            DO J=1,NCOIL
               IF (KCOILCURR.EQ.5.AND.ALLOCATED(TMPCC)) TMP = TMPCC(J)
               R(LX,I) = R(LX,I) + CFEED2(MS,J)*TMPR + 
     &                   CFEED(MS,J)*FEEDI(J)*TMPR*TMP
            ENDDO
         ENDDO
        
         IF (ALLOCATED(TMPCC)) DEALLOCATE(TMPCC)
      ENDIF

      IF (NCASE.EQ.1.AND.AOMEGA<1.0e+5) THEN
         DO MS=1,MSMAX
            LX      = (MS-1)*NXCOMP + KXB1
            DO J=1,NCOIL
               R(LX,I) = R(LX,I) + CFEED(MS,J)*FEEDI(J)
            ENDDO
         ENDDO
      ENDIF

      IF (NCASE.EQ.3.AND.KCOILCURR.EQ.3) THEN
         IF (ICALPHA2.EQ.0) TMPR = 1.
         IF (ICALPHA2.EQ.1) TMPR = 1./CALPHA2

         DO MS=1,MSMAX
            LX      = (MS-1)*NXCOMP + KXB1
            R(LX,I) = (0.,0.)

            TMP = 0.0
            IF (DFLOAT(NIT).GE.EXTCURR(1).AND.
     *          DFLOAT(NIT).LE.EXTCURR(2)) THEN
               TMP = EXP(CI*EXTCURR(3)*ATAU*(NIT-EXTCURR(1)))
            ENDIF
            TMP = TMP*ATAU*AOMEGA

            DO J=1,NCOIL
               R(LX,I) = R(LX,I) + CFEED2(MS,J)*TMPR + 
     &                   CFEED(MS,J)*FEEDI(J)*TMP*TMPR    
            ENDDO
         ENDDO
      ENDIF

 1011 CONTINUE

C     --- INCFEED=6: CALCULATIONS WITH ONLY RW CURRENT
      IF (INCFEED.EQ.6) THEN
         DO MS=1, MXMAX
            DO J=1,NTP1
               R(MS,J) = 0.0
            ENDDO
         ENDDO
         DO MS=1,MYMAX
            DO J=1,NTOT
               RY(MS,J) = 0.0
            ENDDO
         ENDDO

         OPEN(CHOUTP,FILE='JRW.IN')
         REWIND(CHOUTP)
         DO J=1,NWALL
            II = NR + IWALL(J)
            DO MS=1,MSMAX
               LX      = (MS-1)*NXCOMP + KXB1
               READ(CHOUTP,*) TMPR,TMPR1,TMPI1,TMPR2,TMPI2
               IF (ABS(RNTOR).GT.1.0E-13) THEN 
                 R(LX,II) = (TMPR1 + TMPI1*CI)
               ELSEIF (ABS(RM(MS,2)).GT.0.1) THEN 
                 R(LX,II) = (TMPR2 + TMPI2*CI)
               ENDIF
            ENDDO
         ENDDO
         CLOSE(CHOUTP)

         IF (ABS(RM(MS,2)).GT.0.1.OR.ABS(RNTOR).GT.1.0E-13) GOTO 1012

         OPEN(CHOUTP,FILE='JPLASMA_BIOT.IN')
         REWIND(CHOUTP)
         READ(CHOUTP,1171) J,II,TMPR,J,J,J

         DO MS=1,MSMAX
            READ(CHOUTP,1172) TMPR,TMPR,TMPR,
     &                        TMPR,TMPR,TMPR
         ENDDO
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP 
            LY = (MS-1)*NYCOMP 

            DO J=1,NRP1
              READ(CHOUTP,1172) TMPR1,TMPI1,TMPR2,TMPI2,TMPR3,TMPI3
            ENDDO

            J = 1
            CJ20 = 0.0
            CJ30 = 0.0

            J = 2
            READ(CHOUTP,1172) TMPR1,TMPI1,TMPR2,TMPI2,TMPR3,TMPI3
            CJ21 = TMPR2 + TMPI2*CI
            CJ31 = TMPR3 + TMPI3*CI

            RY(LY+KYB2,NR+J-1) = VCSH(J-1)*(0.75*CJ30+0.25*CJ31)*0.5
            RY(LY+KYB3,NR+J-1) =-VCSH(J-1)*(0.75*CJ20+0.25*CJ21)*0.5
                
            DO J=3,II-NR
              READ(CHOUTP,1172) TMPR1,TMPI1,TMPR2,TMPI2,TMPR3,TMPI3
              CJ22 = TMPR2 + TMPI2*CI
              CJ32 = TMPR3 + TMPI3*CI
              IF (ABS(RM(MS,2)).LT.0.1.AND.
     &            ABS(RNTOR).LT.1.0E-13) THEN
                RY(LY+KYB2,NR+J-1) = RY(LY+KYB2,NR+J-2) + 
     &                (VCSH(J-2)*CJ30*0.125+(VCSH(J-2)+VCSH(J-1))*
     &                 CJ31*0.375+VCSH(J-1)*CJ32*0.125)
                RY(LY+KYB3,NR+J-1) = RY(LY+KYB3,NR+J-2) - 
     &                (VCSH(J-2)*CJ20*0.125+(VCSH(J-2)+VCSH(J-1))*
     &                 CJ21*0.375+VCSH(J-1)*CJ22*0.125)
                CJ20 = CJ21
                CJ30 = CJ31
                CJ21 = CJ22
                CJ31 = CJ32
              ENDIF
            ENDDO
         ENDDO
         CLOSE(CHOUTP)
 1012    CONTINUE
      ENDIF
      
C     INCFEED=18: BIOT-SAVART FOR PLASMA CURRENT
      IF (INCFEED.EQ.18) THEN
         DO MS=1, MSMAX
            LX      = (MS-1)*NXCOMP + KXB1
            DO J=1,NRP1
               R(LX,J) = 0.0
            ENDDO
         ENDDO
         DO MS=1,MSMAX
            LY      = (MS-1)*NYCOMP
            DO J=1,NR
               RY(LY+KYB2,J) = 0.0
               RY(LY+KYB3,J) = 0.0
            ENDDO
         ENDDO

         OPEN(CHOUTP,FILE='JPLASMA_BIOT.IN')
         REWIND(CHOUTP)
         READ(CHOUTP,1171) J,II,TMPR,J,J,J

         DO MS=1,MSMAX
            READ(CHOUTP,1172) TMPR,TMPR,TMPR,
     &                        TMPR,TMPR,TMPR
         ENDDO
         DO MS=1,MSMAX
            LX = (MS-1)*NXCOMP 
            LY = (MS-1)*NYCOMP 

            J = 1
            READ(CHOUTP,1172) TMPR1,TMPI1,TMPR2,TMPI2,TMPR3,TMPI3
            CJ20 = TMPR2 + TMPI2*CI
            CJ30 = TMPR3 + TMPI3*CI
            IF (ABS(RNTOR).GT.1.0E-13) THEN
              R(LX+KXB1,J) = CJ20
            ELSEIF (ABS(RM(MS,2)).GT.0.1) THEN 
              R(LX+KXB1,J) = CJ30
            ENDIF
            IF (ABS(RM(MS,2)).GT.0.1.OR.ABS(RNTOR).GT.1.0E-13)   
     &        RY(LY+KYB2,J) = (TMPR1 + TMPI1*CI)

            J = 2
            READ(CHOUTP,1172) TMPR1,TMPI1,TMPR2,TMPI2,TMPR3,TMPI3
            CJ21 = TMPR2 + TMPI2*CI
            CJ31 = TMPR3 + TMPI3*CI
            IF (ABS(RNTOR).GT.1.0E-13) THEN
              R(LX+KXB1,J) = CJ21
            ELSEIF (ABS(RM(MS,2)).GT.0.1) THEN 
              R(LX+KXB1,J) = CJ31
            ENDIF
            IF (ABS(RM(MS,2)).GT.0.1.OR.ABS(RNTOR).GT.1.0E-13)
     &        RY(LY+KYB2,J) = (TMPR1 + TMPI1*CI)

            IF (ABS(RM(MS,2)).LT.0.1.AND.ABS(RNTOR).LT.1.0E-13) THEN
              RY(LY+KYB2,J-1) = CSH(J-1)*(0.75*CJ30+0.25*CJ31)*0.5
              RY(LY+KYB3,J-1) =-CSH(J-1)*(0.75*CJ20+0.25*CJ21)*0.5
            ENDIF
                
            DO J=3,II  
              READ(CHOUTP,1172) TMPR1,TMPI1,TMPR2,TMPI2,TMPR3,TMPI3
              CJ22 = TMPR2 + TMPI2*CI
              CJ32 = TMPR3 + TMPI3*CI
              IF (J.LE.NRP1) THEN
              IF (ABS(RNTOR).GT.1.0E-13) THEN
                R(LX+KXB1,J) = CJ22
              ELSEIF (ABS(RM(MS,2)).GT.0.1) THEN 
                R(LX+KXB1,J) = CJ32
              ENDIF
              ENDIF
              IF (J.LE.NR.AND.(ABS(RM(MS,2)).GT.0.1.OR.
     &            ABS(RNTOR).GT.1.0E-13))   
     &          RY(LY+KYB2,J) = (TMPR1 + TMPI1*CI)
              IF (J.LE.NRP1.AND.ABS(RM(MS,2)).LT.0.1.AND.
     &            ABS(RNTOR).LT.1.0E-13) THEN
                RY(LY+KYB2,J-1) = RY(LY+KYB2,J-2) + 
     &                (CSH(J-2)*CJ30*0.125+(CSH(J-2)+CSH(J-1))*
     &                 CJ31*0.375+CSH(J-1)*CJ32*0.125)
                RY(LY+KYB3,J-1) = RY(LY+KYB3,J-2) - 
     &                (CSH(J-2)*CJ20*0.125+(CSH(J-2)+CSH(J-1))*
     &                 CJ21*0.375+CSH(J-1)*CJ22*0.125)
                CJ20 = CJ21
                CJ30 = CJ31
                CJ21 = CJ22
                CJ31 = CJ32
              ENDIF
            ENDDO
         ENDDO
         CLOSE(CHOUTP)
 1171 FORMAT(I5,1X,I5,1X,E9.2,3(1X,I2))
 1172 FORMAT(E16.8,5(1X,E16.8))
      ENDIF
      
      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C----------OUTPUT FEEDBACK RESULTS--------Y.Q.LIU 15/04/1999------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE FEEDOUT
C     ==================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE GIJLM
      USE RCOMDM
      USE FEEDBACKM
      USE TORQUEM
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'

      INTEGER     M, MS, I, K, NCHN, J, II, IFD, ML, MR, NFC
      PARAMETER   (NCHN=128, NFC=1025)
      REAL*8      CHII(NCHN), HCHII, CHIAA, CHIBB, CHICC, CHIF,
     R            CHIA, CHIB, QTMP1, QTMP2,B1MN
      COMPLEX*16  JCHI2(NCHN), JCHI3(NCHN), HEVIS(NCHN), TMP,
     C            JJJ(NCOIL0), PSIF, ALF, ALFS, ALFT, ALFP, TMP1, TMP2,
     C            FDI(NCOIL0)
      COMPLEX*16  BRPICK,BZPICK,BFPICK

      COMPLEX*16  WORK(2*NCHN), SQRTG11, SQRTG22
      INTEGER     NWORK
      REAL*8      TMPR
      COMPLEX*16  G 
      
      COMPLEX*16,DIMENSION(:),ALLOCATABLE:: PSIS,TPSIM 

      I   = NR + IFEED
      IFD = INCFEED

C     --- CALCULATE FEEDI: METHOD IV
      WRITE(*,*) 'TOTAL FEEDBACK CURRENT ...'
      WRITE(*,70)
 70   FORMAT('METHOD',1X,'Re(FEEDI)',3X,'Im(FEEDI)')
      DO K=1,NCOIL
         CHIF     = -1.0
         IF (NFC.GT.1) HCHII    = (FCCHI(K)+1.0)/(NFC-1)
         FDI(K) = (0.0,0.0)
         TMP2     = (0.0,0.0)
         DO J=1,NFC-1
            TMP1 = (0.0,0.0)
            DO MS=1,MSMAX
               TMP  = EXP(CI*RM(MS,2)*CHIF*PI)
               TMP1 = TMP1 + J3U(I,MS)*TMP
            ENDDO
            FDI(K) = FDI(K) + (TMP1+TMP2)
            TMP2     = TMP1
            CHIF     = CHIF + HCHII
         ENDDO
         FDI(K) = FDI(K)*HCHII*PI/2.0
         FDI(K) = -FDI(K)*(VCS(IFEED+1)-VCS(IFEED-1))/2.0
         OFEEDI(K) = FDI(K)
      ENDDO
      WRITE(*,80) (DREAL(FDI(K)),DIMAG(FDI(K)),K=1,NCOIL)
 80   FORMAT('  IV  ',10(E11.4,1X))


C     --- CALCULATE FEEDI: METHOD III
      IF (ABS(RNTOR).GT.1.E-10) THEN
      DO K=1,NCOIL
         FDI(K) = (0.,0.)
         DO MS=1,MSMAX
            TMP      = EXP(CI*RM(MS,2)*FCCHI(K)*PI)
            FDI(K) = FDI(K) + J2U(I,MS)*TMP
         ENDDO
         FDI(K) = FDI(K)*(VCS(IFEED+1)-VCS(IFEED-1))/2.0
         FDI(K) = FDI(K)/CI/RNTOR
      ENDDO
      WRITE(*,90) (DREAL(FDI(K)),DIMAG(FDI(K)),K=1,NCOIL)
 90   FORMAT(' III  ',10(E11.4,1X))
      ENDIF

C     --- CALCULATE FEEDI: METHOD II
      DO K=1,NCOIL
         JJJ(K) = 0.0
         DO MS=1,MSMAX
            IF (ABS(RNTOR).LT.1.0E-10) THEN
               JJJ(K) = JJJ(K) + CFEED(MS,K)*J3U(I,MS)
            ELSE
               JJJ(K) = JJJ(K) + CFEED(MS,K)*J2U(I,MS)
            ENDIF
         ENDDO
      ENDDO
      DO K=1,NCOIL
         FDI(K) = (0.,0.)
         DO MS=1,NCOIL
            FDI(K) = FDI(K) + DFEED(K+(MS-1)*NCOIL)*JJJ(MS)
         ENDDO
      ENDDO
      WRITE(*,100) (DREAL(FDI(K)),DIMAG(FDI(K)),K=1,NCOIL)
 100  FORMAT('  II  ',10(E11.4,1X))

C     --- TEST CURRENT DECOMPOSITION
      WRITE(*,*) 'TEST COIL CURRENT'
      DO MS=1,MSMAX
         TMP = 0.0
         DO K=1,NCOIL
C           TMP = TMP + CFEED(MS,K)*FDI(K)
            TMP = TMP + CFEED(MS,K)
         ENDDO
CYQL     WRITE(*,101) MS,J2U(I,MS),TMP
      ENDDO
 101  FORMAT(I4,4(E11.4,1X))

C     --- OUTPUT FEEDBACK COILS CURRENT
      OPEN(CHOUTP,FILE='FEEDI.OUT')
      REWIND(CHOUTP)
      DO K=1,NCOIL
         WRITE(CHOUTP,110) DREAL(FDI(K)),DIMAG(FDI(K))
      ENDDO
      CLOSE(CHOUTP)


C     --- CALCULATE MAGN. FLUX 'BS' AND 'BCHI' IN THE VACUUM REGION
      OPEN(CHOUTP,FILE='FLUX.OUT')
      REWIND(CHOUTP)
      DO J=1,NV
         DO K=1,NCOIL
            JJJ(K) = 0.0
            DO MS=1,MSMAX
               JJJ(K) = JJJ(K) + BFEEDSS(K,MS)*B1U(NR+J,MS)
            ENDDO
         ENDDO
         WRITE(CHOUTP,102) VCS(J),(DREAL(JJJ(K)),K=1,NCOIL),
     *                             (DIMAG(JJJ(K)),K=1,NCOIL)
      ENDDO
      DO J=1,NV-1
         DO K=1,NCOIL
            JJJ(K) = 0.0
            DO MS=1,MSMAX
               JJJ(K) = JJJ(K) + BFEEDST(K,MS)*B2U(NR+J,MS)
            ENDDO
         ENDDO
         WRITE(CHOUTP,102) VCSM(J),(DREAL(JJJ(K)),K=1,NCOIL),
     *                             (DIMAG(JJJ(K)),K=1,NCOIL)
      ENDDO
      CLOSE(CHOUTP)
 102  FORMAT(E11.4,2X,18(E12.5,1X))

C     --- CALCULATE SELF-INDUCTANCE OF ACTIVE COILS
      WRITE(*,*) ' SELF-INDUCTANCE OF ACTIVE COILS:'
      DO K=1,NCOIL
         PSIF = 0.0
         DO MS=1,MSMAX
            PSIF = PSIF + BFEEDF(K,MS)*B1U(IFEED+NR,MS)
         ENDDO
         WRITE(*,103) IFEED, K, DREAL(PSIF),DIMAG(PSIF)
         IF (IFD.EQ.2.OR.IFD.EQ.4.OR.IFD.EQ.10) RLK(K) = DREAL(ALF)
      ENDDO

C     --- CALCULATE MUTUAL INDUCTANCE BETWEEN ACTIVE COIL
C     --- AND SENSOR LOOP
      WRITE(*,*) ' Br-MUTUAL INDUCTANCE BETWEEN ACT. COIL & SENSOR,'
      WRITE(*,*) ' FLUX SENSOR:'
      DO J=0,NSENS
         DO K=1,NCOIL
            PSIF = 0.0
            DO MS=1,MSMAX
               PSIF = PSIF + BFEEDSS(K,MS)*B1U(ISENS(J)+NR,MS)
            ENDDO
            WRITE(*,504) ISENS(J), K, DREAL(PSIF),DIMAG(PSIF)
            IF (IFD.EQ.2) RMFS(K) = DREAL(ALF)
         ENDDO
      ENDDO

      WRITE(*,*) ' Br-MUTUAL INDUCTANCE BETWEEN ACT. COIL & SENSOR,'
      WRITE(*,*) ' POINT FIELD SENSOR:'
      DO J=0,NSENS
         DO K=1,NCOIL
            PSIF = 0.0
            DO MS=1,MSMAX
               PSIF = PSIF + BFEEDST(K,MS)*B1U(ISENS(J)+NR,MS)/CI
            ENDDO
            WRITE(*,505) ISENS(J), K, DREAL(PSIF),DIMAG(PSIF)
            IF (IFD.EQ.2) RMFS(K) = DREAL(ALF)
         ENDDO
      ENDDO

C     --- WRITE Br-FIELD AROUND THE POLOIDAL ANGLE
      OPEN(CHOUTP,FILE='EFAF.OUT')
      REWIND(CHOUTP)
      TMP = 2.0*PI/200
      DO J=0,200
         TMP1 = TMP*(J-100)
         PSIF = 0.0
         DO MS=1,MSMAX
            TMP2 = EXP(TMP1*RM(MS,2)*CI)
            PSIF = PSIF + TMP2*B1U(ISENS(NSENS)+NR,MS)
         ENDDO
         WRITE(CHOUTP,191) DREAL(TMP1),DREAL(PSIF),DIMAG(PSIF)
      ENDDO
 191  FORMAT(3E16.8)
      CLOSE(CHOUTP)

      J = 0
      IF (IFD.EQ.3) J = 1
      DO K=1,NCOIL
         PSIF = 0.0
         DO MS=1,MSMAX
            PSIF = PSIF + BFEEDST(K,MS)*B1U(ISENS(1)+NR+J,MS)/CI
         ENDDO
         OSENSR(K) = PSIF
      ENDDO

      WRITE(*,*) ' Bt-MUTUAL INDUCTANCE BETWEEN ACT. COIL & SENSOR:'
      DO J=0,NSENS
         DO K=1,NCOIL
            PSIF = 0.0
            DO MS=1,MSMAX
               PSIF = PSIF + BFEEDST(K,MS)*B2U(ISENS(J)+NR,MS)
            ENDDO
            WRITE(*,506) ISENS(J), K, DREAL(PSIF),DIMAG(PSIF)
            
            PSIF = 0.0
            DO MS=1,MSMAX
               PSIF = PSIF + BFEEDST(K,MS)*B2U(ISENS(J)+NR-1,MS)
            ENDDO
            WRITE(*,507) ISENS(J)-1, K, DREAL(PSIF),DIMAG(PSIF)
         ENDDO
      ENDDO
 103  FORMAT('ISENS=',I3,' K=',I2,' ALF=',2(E12.5,1X))
 104  FORMAT('ISENS=',I3,' K=',I2,' ALFF=',2(E12.5,1X))
 504  FORMAT('SENSOR Flux : ',I3,1X,I2,2(1X,E12.5))
 505  FORMAT('SENSOR Bn   : ',I3,1X,I2,2(1X,E12.5))
 506  FORMAT('SENSOR BtExt: ',I3,1X,I2,2(1X,E12.5))
 507  FORMAT('SENSOR BtInt: ',I3,1X,I2,2(1X,E12.5))
      
      J = 0
      IF (IFD.EQ.11) J = 1
      DO K=1,NCOIL
         PSIF = 0.0
         DO MS=1,MSMAX
            PSIF = PSIF + BFEEDST(K,MS)*B2U(ISENS(1)+NR-J,MS)
         ENDDO
         OSENSP(K) = PSIF
      ENDDO

      DO K=1,NCOIL
         WRITE(*,105) K, RLK(K)
      ENDDO
 105  FORMAT(' RLK(',I2,')=',1(E11.4,1X))

      WRITE(*,*) ' Bt/Br(rs) = '
      DO J=0,NSENS
         DO K=1,NCOIL
            PSIF = 0.0
            DO MS=1,MSMAX
               PSIF = PSIF + BFEEDST(K,MS)*B1U(ISENS(J)+NR,MS)/CI
            ENDDO
            ALF = PSIF
            PSIF = 0.0
            DO MS=1,MSMAX
               PSIF = PSIF + BFEEDST(K,MS)*B2U(ISENS(J)+NR-1,MS)
            ENDDO
            IF (ABS(ALF).NE.0.0) ALF = PSIF/ALF
            WRITE(*,103) ISENS(J), K, ABS(ALF),DREAL(ALF),DIMAG(ALF)
         ENDDO
      ENDDO

C     CALCULATE GEOMETRIC FACTORS SQRTG11 & SQRTG22 AT THE SENSOR POSITION
      WRITE(*,*) ' SQRT_G11 AND SQRT_G22 AT THE SENSOR POSITION:'
      DO J=0,NSENS
         SQRTG11 = 0.0
         SQRTG22 = 0.0
         IF (ISENS(J).GE.1) THEN
         DO MS=1,MEDIM
            SQRTG11 = SQRTG11 + VDG11LM(ISENS(J),MS)
            SQRTG22 = SQRTG22 + VDG22LM(ISENS(J),MS)
         ENDDO
         SQRTG11 = SQRT(SQRTG11)
         SQRTG22 = SQRT(SQRTG22)
         WRITE(*,107) J, DREAL(SQRTG11),DIMAG(SQRTG11),
     *                   DREAL(SQRTG22),DIMAG(SQRTG22),
     *                   DREAL(SQRTG11/SQRTG22),DIMAG(SQRTG11/SQRTG22)
         ENDIF
      ENDDO
 107  FORMAT('J=',I1,' SQRT_G11:SQRT_G22=',6(E11.4,1X))

C     --- OUTPUT JCHI(IFEED,CHI) FOR MATLAB
      HCHII   = 2./NCHN
      CHII(1) = -1.
      DO J=2,NCHN
         CHII(J) = CHII(J-1) + HCHII
      ENDDO
      
      DO J=1,NCHN
         HEVIS(J) = 0.
         DO K=1,NCOIL
            CHIA  = FCCHI(K) - 0.5*FWCHI(K) - 0.5*BTCHI
            CHIB  = FCCHI(K) + 0.5*FWCHI(K) + 0.5*BTCHI
            CHIAA = CHIA + BTCHI
            CHIBB = CHIB - BTCHI
            IF (CHIA.LT.CHII(J).AND.CHII(J).LT.CHIAA) THEN
               HEVIS(J) = HEVIS(J) + FDI(K)*(CHII(J)-CHIA)/BTCHI
            ENDIF
            IF (CHIAA.LE.CHII(J).AND.CHII(J).LE.CHIBB) THEN
               HEVIS(J) = HEVIS(J) + FDI(K)
            ENDIF
            IF (CHIBB.LT.CHII(J).AND.CHII(J).LT.CHIB) THEN
               HEVIS(J) = HEVIS(J) + FDI(K)*(CHIB-CHII(J))/BTCHI
            ENDIF
         ENDDO
         IF (NV.GE.2) HEVIS(J) = 
     *                HEVIS(J)*2*CI*RNTOR/(VCS(IFEED+1)-VCS(IFEED-1))
      ENDDO
         
      DO J=1,NCHN
         JCHI2(J) = (0.,0.)
         JCHI3(J) = (0.,0.)
         DO MS=1,MSMAX
            TMP      = EXP(CI*RM(MS,2)*CHII(J)*PI)
            JCHI2(J) = JCHI2(J) + J2U(I,MS)*TMP
            JCHI3(J) = JCHI3(J) + J3U(I,MS)*TMP
         ENDDO
      ENDDO
      OPEN(CHOUTP,FILE='JFD.OUT')
      REWIND(CHOUTP)
      DO J=1,NCHN
         WRITE(CHOUTP,110) CHII(J),ABS(HEVIS(J)),DREAL(JCHI2(J)),
     *        DIMAG(JCHI2(J)),DREAL(JCHI3(J)),DIMAG(JCHI3(J)),
     *        DREAL(HEVIS(J))
      ENDDO
      CLOSE(CHOUTP)
 110  FORMAT(7(E11.4,1X))

C     --- OUTPUT COEFF. <BFEEDF>
      OPEN(CHOUTP,FILE='BFEEDF.OUT')
      REWIND(CHOUTP)
      DO MS=1,MSMAX
         WRITE(CHOUTP,115) IM(MS,2),(DREAL(BFEEDF(K,MS)),K=1,NCOIL),
     *                         (DIMAG(BFEEDF(K,MS)),K=1,NCOIL)
      ENDDO
      CLOSE(CHOUTP)
 115  FORMAT(I3,1X,10(E11.4,1X))

C     --- OUTPUT COEFF. <CFEED>
      OPEN(CHOUTP,FILE='CFEED.OUT')
      REWIND(CHOUTP)
      DO MS=1,MSMAX
         WRITE(CHOUTP,115) IM(MS,2),(DREAL(CFEED(MS,K)),K=1,NCOIL),
     *                         (DIMAG(CFEED(MS,K)),K=1,NCOIL)
      ENDDO
      CLOSE(CHOUTP)

C     --- OUTPUT THE EDDY CURRENTS IN THE RESISTIVE WALL
      OPEN(CHOUTP,FILE='JRW.OUT')
      REWIND(CHOUTP)
      DO J=1,NWALL
         II = NR + IWALL(J)
         DO MS=1,MSMAX
            WRITE(CHOUTP,117) RM(MS,2),J2U(II,MS),J3U(II,MS)
         ENDDO
      ENDDO
      CLOSE(CHOUTP)
 117  FORMAT(F5.1,1X,4E16.8)

C     --- OUTPUT THE CURRENTS IN THE PLASMA SURFACE
      II = NRP1
      OPEN(CHOUTP,FILE='JPS.OUT')
      REWIND(CHOUTP)
      DO MS=1,MSMAX
         WRITE(CHOUTP,117) RM(MS,2),J2U(II,MS),J3U(II,MS)
      ENDDO
      CLOSE(CHOUTP)

C     --- OUTPUT RADIAL MAGNETIC FIELD AT THE PLASMA SURFACE AND
C     --- EDGE PLASMA PERTURBATION VELOCITY
      II = NRP1
      DO J=1,NCHN
         JCHI2(J) = (0.,0.)
         JCHI3(J) = (0.,0.)
         DO MS=1,MSMAX
            TMP      = EXP(CI*RM(MS,2)*CHII(J)*PI)
            JCHI2(J) = JCHI2(J) + B1U(II,MS)*TMP
            JCHI3(J) = JCHI3(J) + V1U(II,MS)*TMP
         ENDDO
      ENDDO
      OPEN(CHOUTP,FILE='B1UP.OUT')
      REWIND(CHOUTP)
      DO J=1,NCHN
         WRITE(CHOUTP,110) CHII(J),DREAL(JCHI2(J)),DIMAG(JCHI2(J)),
     *                             DREAL(JCHI3(J)),DIMAG(JCHI3(J))
      ENDDO
      CLOSE(CHOUTP)

C     --- OUTPUT RADIAL MAGNETIC FIELD AT THE FIRST RESISTIVE WALL
      IF (NWALL.GT.0) THEN
         II = NR + IWALL(1)
         DO J=1,NCHN
            JCHI2(J) = (0.,0.)
            DO MS=1,MSMAX
               TMP      = EXP(CI*RM(MS,2)*CHII(J)*PI)
               JCHI2(J) = JCHI2(J) + B1U(II,MS)*TMP
            ENDDO
         ENDDO
         OPEN(CHOUTP,FILE='B1UW.OUT')
         REWIND(CHOUTP)
         DO J=1,NCHN
            WRITE(CHOUTP,110) CHII(J),DREAL(JCHI2(J)),DIMAG(JCHI2(J))
         ENDDO
         CLOSE(CHOUTP)
      ENDIF

C     --- OUTPUT RADIAL MAGNETIC FIELD AT SENSOR COIL POSITION
      II = NR + ISENS(1)
      DO J=1,NCHN
         JCHI2(J) = (0.,0.)
         DO MS=1,MSMAX
            TMP      = EXP(CI*RM(MS,2)*CHII(J)*PI)
            JCHI2(J) = JCHI2(J) + B1U(II,MS)*TMP
         ENDDO
      ENDDO
      IF (INCFEED.EQ.0.OR.GAINA(1).EQ.0.0) THEN
         OPEN(CHOUTP,FILE='B1US00.OUT')
         REWIND(CHOUTP)
      ENDIF
      IF ((IFD.EQ.1.OR.IFD.EQ.11.OR.IFD.EQ.3.OR.IFD.EQ.5)
     *    .AND.GAINA(1).NE.0.0) THEN
         OPEN(CHOUTP,FILE='B1US.OUT')
         REWIND(CHOUTP)
         WRITE(CHOUTP,110) CHII(1),DREAL(FDI(1)),DIMAG(FDI(1))
      ENDIF
      IF (INCFEED.EQ.2.OR.INCFEED.EQ.4) THEN
         OPEN(CHOUTP,FILE='B1US0.OUT')
         REWIND(CHOUTP)
      ENDIF
      DO J=1,NCHN
         WRITE(CHOUTP,110) CHII(J),DREAL(JCHI2(J)),DIMAG(JCHI2(J))
      ENDDO
      CLOSE(CHOUTP)

C     --- CONTINUE CS UNTIL NTP1
      DO J=NRP1+1,NTP1
         CS(J) = VCS(J-NRP1+1)
      ENDDO

      ML=1
      MR=MSMAX
      
C     --- OUTPUT B2U(I,MS), MS=1,MSMAX JUST INSIDE COIL POSITION
C     --- EXTRAPOLATED FROM OUTSIDE OF THE COIL
C     --- FOR BACKWARD COUPLING ONLY
      OPEN(CHOUTP,FILE='B2U.OUT')
      REWIND(CHOUTP)
      J = NR + IFEED 
      WRITE(*,*) CSM(J-1),CSM(J),CSM(J+1)
      DO MS=ML,MR
         WORK(MS)=B2U(J+1,MS)+(CSM(J-1)-CSM(J+1))/(CSM(J)-CSM(J+1))
     &            *(B2U(J,MS)-B2U(J+1,MS))
      ENDDO
      WRITE(CHOUTP,120) (DREAL(WORK(MS)),MS=ML,MR),
     *                  (DIMAG(WORK(MS)),MS=ML,MR)
      CLOSE(CHOUTP)

C     --- OUTPUT B3U(I,MS), MS=1,MSMAX JUST INSIDE COIL POSITION
      OPEN(CHOUTP,FILE='B3U.OUT')
      REWIND(CHOUTP)
      J = NR + IFEED 
      DO MS=ML,MR
         WORK(MS)=B3U(J+1,MS)+(CSM(J-1)-CSM(J+1))/(CSM(J)-CSM(J+1))
     &            *(B3U(J,MS)-B3U(J+1,MS))
      ENDDO
      WRITE(CHOUTP,120) (DREAL(WORK(MS)),MS=ML,MR),
     *                  (DIMAG(WORK(MS)),MS=ML,MR)
      CLOSE(CHOUTP)

 120  FORMAT(600(E13.6,1X))

C     --- OUTPUT B1U(I,MS),B2U(I-1,MS),B3U(I-1,MS)
C     --- AT A GIVEN RADIAL LOCATION I=IBOUT, WHERE
C     --- 0 < IBOUT < NR+NV+1
      IF (IBOUT.GT.0.AND.IBOUT.LT.NRP1+NV) THEN
      OPEN(CHOUTP,FILE='B123U.OUT')
      REWIND(CHOUTP)
      J = IBOUT
      WRITE(CHOUTP,120) (DREAL(B1U(J,MS)),MS=ML,MR),
     *                  (DIMAG(B1U(J,MS)),MS=ML,MR),
     *                  (DREAL(B2U(J-1,MS)),MS=ML,MR),
     *                  (DIMAG(B2U(J-1,MS)),MS=ML,MR),
     *                  (DREAL(B3U(J-1,MS)),MS=ML,MR),
     *                  (DIMAG(B3U(J-1,MS)),MS=ML,MR)
      CLOSE(CHOUTP)
      ENDIF

C     --- OUTPUT B1U(I,MS),B2U(I-1,MS),B3U(I-1,MS)
C     --- AT A GIVEN RADIAL LOCATION I=IBOUT2, WHERE
C     --- 0 < IBOUT < NR+NV+1
      IF (IBOUT2.GT.0.AND.IBOUT2.LT.NRP1+NV) THEN
      OPEN(CHOUTP,FILE='B123U2.OUT')
      REWIND(CHOUTP)
      J = IBOUT2
      WRITE(CHOUTP,120) (DREAL(B1U(J,MS)),MS=ML,MR),
     *                  (DIMAG(B1U(J,MS)),MS=ML,MR),
     *                  (DREAL(B2U(J-1,MS)),MS=ML,MR),
     *                  (DIMAG(B2U(J-1,MS)),MS=ML,MR),
     *                  (DREAL(B3U(J-1,MS)),MS=ML,MR),
     *                  (DIMAG(B3U(J-1,MS)),MS=ML,MR)
      CLOSE(CHOUTP)
      ENDIF

C     --- OUTPUT J2U(I,MS) ON RESISTIVE WALL & FEEDBACK, MS=1,MSMAX.
      II = NR + IWALL(1)
      OPEN(CHOUTP,FILE='J2UWF.OUT')
      REWIND(CHOUTP)
      DO MS=1,MSMAX
         WRITE(CHOUTP,130) IM(MS,2),DREAL(J2U(II,MS)),DIMAG(J2U(II,MS)),
     *        DREAL(J2U(NR+IFEED,MS)),DIMAG(J2U(NR+IFEED,MS))
      ENDDO
      CLOSE(CHOUTP)
 130  FORMAT(I3,1X,4(E11.4,1X))
 140  FORMAT(I3,1X,140(E11.4,1X))

C     --- OUTPUT J3U(I,MS) ON RESISTIVE WALL & FEEDBACK, MS=1,MSMAX.
      II = NR + IWALL(1)
      OPEN(CHOUTP,FILE='J3UWF.OUT')
      REWIND(CHOUTP)
      DO MS=1,MSMAX
         WRITE(CHOUTP,130) IM(MS,2),DREAL(J3U(II,MS)),DIMAG(J3U(II,MS)),
     *        DREAL(J3U(NR+IFEED,MS)),DIMAG(J3U(NR+IFEED,MS))
      ENDDO
      CLOSE(CHOUTP)

C     --- OUTPUT SENSOR SIGNALS FOR BACKWARD COUPLING SCHEME
      IF (NCOUPL.EQ.-2.AND.INCFEED.GT.0) THEN
         IF (.NOT.ALLOCATED(PSIS)) ALLOCATE( PSIS(NSENST) )
         PSIS = (0.,0.)
         II = NTP1
         DO J=1,NSENST
            DO MS=1,MSMAX
               PSIS(J)=PSIS(J) - (B01S(J,MS)+AL0*B11S(J,MS))*B1U(II,MS)
     &                         + AL0*A11S(J,MS)*B2U(II-1,MS)
            ENDDO
            DO MS=1,NCOILT
               PSIS(J)=PSIS(J) + (C01S(J,MS)+AL0*C11S(J,MS))*FEEDIT(MS)
            ENDDO
         ENDDO
      ENDIF
      IF (NCOUPL.EQ.-3.AND.INCFEED.GT.0) THEN
         IF (.NOT.ALLOCATED(PSIS)) ALLOCATE( PSIS(NSENST) )
         PSIS = (0.,0.)
         II = NTP1
         DO J=1,NSENST
            DO MS=1,MSMAX
               PSIS(J)=PSIS(J) - (B02S(J,MS)+AL0*B12S(J,MS))*B1U(II,MS)
     &                         + AL0*A12S(J,MS)*B3U(II-1,MS)
            ENDDO
            DO MS=1,NCOILT
               PSIS(J)=PSIS(J) + (C02S(J,MS)+AL0*C12S(J,MS))*FEEDIT(MS)
            ENDDO
         ENDDO
      ENDIF
      IF ((NCOUPL.EQ.-2.OR.NCOUPL.EQ.-3).AND.INCFEED.GT.0) THEN
         DO J=1,NSENST
            WRITE(*,555) J,PSIS(J)
         ENDDO
      ENDIF
 555  FORMAT('SENSOR CarMa: ',I3,1X,4(E12.5,1X))

C     CONVERT B1-FIELD HARMONICS TO SFL COORDINATE SYSTEM
      IF (NCASE.EQ.1.OR.NCASE.EQ.2) THEN
      IF (.NOT.ALLOCATED(TPSIM)) ALLOCATE(TPSIM(NRATSURF))
      DO J=1,NRATSURF
         I = IRATSURF(J)
         TPSIM(J) = (0.,0.)
         M = 1
         DO MS=1,MSMAX
            IF (ABS(-Q(I)*RNTOR-RM(MS,2)).LT.0.4) M=MS
         ENDDO
         QTMP1 = -RM(M,2)/RNTOR
         QTMP2 = (QTMP1-Q(I))/(Q(I+1)-Q(I))   
C        IF (ABS(QTMP2).GT.1.) STOP 'FEEDOUT: INACCURATE MESH'
         IF (NCONVB1.EQ.0) THEN
            TPSIM(J) = B1U(I,M) + QTMP2*(B1U(I+1,M)-B1U(I,M))
         ELSEIF (NCONVB1.EQ.2) THEN
            DO MS=1,MSMAX
               TMP1     = B1U(I,MS) + QTMP2*(B1U(I+1,MS)-B1U(I,MS))
               TPSIM(J) = TPSIM(J)+TQB1MAT(MS,J)*TMP1
            ENDDO
         ENDIF
      ENDDO

C     PRINT B^1MN AND ISLAND WIDTH AT RATIONAL SURFACES
C     ISLAND WIDTH IN TERMS OF DELTA_S
C     BELOW: TMPR=MAGNETIC SHEAR
      DO J=1,NRATSURF
         I    = IRATSURF(J)
         TMPR = (Q(I+1)-Q(I-1))/(CS(I+1)-CS(I-1))*CS(I)/Q(I)
         B1MN = 0.0
         IF (KISLAND.EQ.1) B1MN = ABS(TPSIM(J))
         IF (KISLAND.EQ.2) B1MN = 2.0*REAL(TPSIM(J))
         WRITE(*,121) -Q(I)*RNTOR,CS(I),TPSIM(J),
     &                4.*SQRT(ABS(B1MN/(RNTOR*Q(I)*TMPR*
     &                                      DPSIDS(NRP1))))
      ENDDO
 121  FORMAT('OUTPUT TPSIM: ',E10.2,4(1X,E14.6))
      ENDIF         

C     OUTPUT PERTURBED BR,BZ,BPHI COMPONENT AT GIVEN PICKUP COIL LOCATIONS
      WRITE(*,*) 'PICKUP SENSOR DATA'
      DO K=1,NPICK
         TMP1 = CPICK(K)*PI
         ALFS = 0.0
         ALFT = 0.0
         ALFP = 0.0
         DO MS=1,MSMAX
            TMP2 = EXP(TMP1*RM(MS,2)*CI)
            ALFS = ALFS + TMP2*
     &             (B1U(IPICK(K)+NR,MS)+B1U(IPICK(K)+NR+1,MS))/2.
            ALFT = ALFT + TMP2*B2U(IPICK(K)+NR,MS)
            ALFP = ALFP + TMP2*B3U(IPICK(K)+NR,MS)
         ENDDO
         BZPICK = VZSJ(K)*ALFS + VZCJ(K)*ALFT
         BRPICK = VRSJ(K)*ALFS + VRCJ(K)*ALFT
         BFPICK = VR3J(K)*ALFP
         WRITE(*,193) K,IPICK(K),CPICK(K),DREAL(BRPICK),DIMAG(BRPICK),
     &                DREAL(BZPICK),DIMAG(BZPICK),
     &                DREAL(BFPICK),DIMAG(BFPICK)

      ENDDO
 193  FORMAT(I3,1X,I3,1X,10E12.4)

      RETURN
      END
      
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------TIME EVOLUTION PROFILE-----Y.Q. LIU 16/08/2001------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE FEEDEVOL(MD,MDY,ND,NCASE,X,Y,XOLD,YOLD,ATAU,IFLAG)
C     ==================================
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE FEEDBACKM
      USE TORQUEM
      USE REORBITM
      INCLUDE 'comioc.inc'
C
      INTEGER       I,J,MS, M, K, IFLAG
      REAL*8        QTMP1,QTMP2,SOLMAXX,SOLMAXY
      COMPLEX*16    TMP,PSIS,PSISS,CTMP1,CTMP2,ALFS,ALFT,PSISNOISE

      CHARACTER(LEN=1024) FILENAME

      INTEGER       MD,MDY,ND,NCASE
      COMPLEX*16    ATAU,X(MD,ND+1),Y(MDY,ND),XOLD(MD,ND+1),YOLD(MDY,ND)

      COMPLEX*16,DIMENSION(:),ALLOCATABLE::TPSI,BZPICK
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::AVF,AIF,PSIF,PSIFOLD

      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, N, MY

C     --- CALCULATE FEEDBACK CURRENT USING METHOD III
      ALLOCATE(AVF(NCOIL),AIF(NCOIL),PSIF(NCOIL),PSIFOLD(NCOIL))

      AIF = 0.0
      DO K=1,NCOIL
         DO MS=1,MSMAX
            TMP = EXP(CI*RM(MS,2)*FCCHI(K)*PI)
            AIF(K) = AIF(K) + X(KXJ2U+(MS-1)*NXCOMP,NR+IFEED)*TMP
         ENDDO
      ENDDO
      AIF = AIF*(VCS(IFEED+1)-VCS(IFEED-1))/2.0
      AIF = AIF/CI/RNTOR

C     --- CALCULATE FLUX THROUGH THE ACTIVE COIL
      PSIFOLD = 0.0
      PSIF    = 0.0
      DO K=1,NCOIL
         DO MS=1,MSMAX
            PSIF(K)=PSIF(K)+BFEEDF(K,MS)*X(KXB1+(MS-1)*NXCOMP,NR+IFEED)
            PSIFOLD(K) = PSIFOLD(K) + 
     &                BFEEDF(K,MS)*XOLD(KXB1+(MS-1)*NXCOMP,NR+IFEED)
         ENDDO
      ENDDO

C     --- CALCULATE THE FEEDBACK VOLTAGE
      AVF = 0.0
      IF (VFMAX.LT.0.0) THEN
         IF (IDYNAM.EQ.0.OR.IDYNAM.EQ.1)  AVF = AVF + AIF
         IF (IDYNAM.GT.0.AND.IDYNAM.LT.3) AVF = AVF + 
     &      TAUF*(PSIF-PSIFOLD)*DREAL(CALPHA1)/RL0
      ELSE
         DO K=1,NCOIL
            AVF(K) = Y((K-1)*NYCOMP+KYV2,NR+IFEED)
         ENDDO
      ENDIF
  
      
C     --- CALCULATE THE SENSOR SIGNAL
      PSIS      = (0.0,0.0)
      PSISS     = (0.0,0.0)
      PSISNOISE = (0.0,0.0)
      J = NSENS
      IF (INCFEED.GT.0) THEN
         DO K=1,NCOIL
            DO MS=1,MSMAX
               PSISS = PSISS + 
     &                BFEEDSS(K,MS)*X(KXB1+(MS-1)*NXCOMP,ISENS(J)+NR)
            ENDDO
         ENDDO
      ENDIF

      IF (INCFEED.EQ.1) PSIS = PSISS

      IF (INCFEED.EQ.11.OR.INCFEED.EQ.12.OR.
     &    INCFEED.EQ.4.OR.INCFEED.EQ.8.OR.INCFEED.EQ.22) THEN
         DO K=1,NCOIL
            DO MS=1,MSMAX
               PSIS = PSIS + 
     &                BFEEDST(K,MS)*X(KXB1+(MS-1)*NXCOMP,ISENS(J)+NR)/CI
               PSISNOISE = PSISNOISE + BFEEDST(K,MS)*SENSNOISE(NIT)/CI
            ENDDO
         ENDDO
      ENDIF
      IF (INCFEED.EQ.3.OR.INCFEED.EQ.9) THEN
         DO K=1,NCOIL
            DO MS=1,MSMAX
               PSIS = PSIS + 
     &                BFEEDST(K,MS)*Y(KYB2+(MS-1)*NYCOMP,ISENS(J)+NR)
               PSISNOISE = PSISNOISE + BFEEDST(K,MS)*SENSNOISE(NIT)
            ENDDO
         ENDDO
      ENDIF

C     CALCULATE POINT-WISE RADIAL FIELD
C     AND CONVERT B1-FIELD HARMONICS TO SFL COORDINATE SYSTEM
      IF (CALPHA7.GT.0.) THEN
      IF (.NOT.ALLOCATED(TPSI)) ALLOCATE(TPSI(2*NRATSURF))
      DO J=1,NRATSURF
         I = IRATSURF(J)
         K = NRATSURF+J

C        TOTAL FIELD AT RATIONAL SURFACE AT CHI=0 ANGLE
         TPSI(K) = (0.,0.)
         DO MS=1,MSMAX
            TPSI(K) = TPSI(K) + X(KXB1+(MS-1)*NXCOMP,I)
         ENDDO

C        RESONANT B1 FIELD AT RATIONAL SURFACES, IN SFL COORDINATE SYSTEM
         TPSI(J) = (0.,0.)

C        IF (ABS(RNTOR).LE.3) THEN
         M = 1
         DO MS=1,MSMAX
            IF (ABS(-Q(I)*RNTOR-RM(MS,2)).LT.0.4) M=MS
         ENDDO
         QTMP1 = -RM(M,2)/RNTOR
         QTMP2 = (QTMP1-Q(I))/(Q(I+1)-Q(I))   
         IF (ABS(QTMP2).GT.1.) STOP 'FEEDEVOL: INACCURATE MESH'
         IF (NCONVB1.EQ.0) THEN
            CTMP1   = X(KXB1+(M-1)*NXCOMP,I)
            CTMP2   = X(KXB1+(M-1)*NXCOMP,I+1)
            TPSI(J) = CTMP1 + QTMP2*(CTMP2-CTMP1)
         ELSEIF (NCONVB1.EQ.2) THEN
            DO MS=1,MSMAX
               CTMP1   = X(KXB1+(MS-1)*NXCOMP,I)
               CTMP2   = X(KXB1+(MS-1)*NXCOMP,I+1)
               TMP1    = CTMP1 + QTMP2*(CTMP2-CTMP1)
               TPSI(J) = TPSI(J)+TQB1MAT(MS,J)*TMP1
            ENDDO
         ENDIF
C        ENDIF
      ENDDO

C     OUTPUT PERTURBED BR BZ BFI COMPONENT AT GIVEN PICKUP COIL LOCATIONS
      IF (.NOT.ALLOCATED(BZPICK)) ALLOCATE(BZPICK(NPICK))
      BZPICK = (0.,0.)
      DO K=1,NPICK
         QTMP1 = CPICK(K)*PI
         ALFS = (0.,0.)
         ALFT = (0.,0.)
         DO MS=1,MSMAX
            CTMP2 = EXP(QTMP1*RM(MS,2)*CI)
            ALFS = ALFS + CTMP2*
     &             (X(KXB1+(MS-1)*NXCOMP,IPICK(K)+NR)+
     &              X(KXB1+(MS-1)*NXCOMP,IPICK(K)+NR+1))/2.
            ALFT = ALFT + CTMP2*Y(KYB2+(MS-1)*NYCOMP,IPICK(K)+NR)
         ENDDO
         BZPICK(K) = VZSJ(K)*ALFS + VZCJ(K)*ALFT
      ENDDO

      ENDIF

C     --- WRITE THE PROFILE EVOLUTION TO FILE
      TMP=(0.,0.)
      IF (NCASE.GE.3.AND.NCASE.LE.6) TMP=1./CALPHA1
      IF (NCASE.EQ.9.OR.NCASE.EQ.10) TMP=1./CALPHA1
      IF (NCASE.EQ.7.OR.NCASE.EQ.8)  TMP=1./CALPHA5

      IF (NCASE.EQ.3) THEN
         WRITE(CHTIME,90) DREAL(TMP),
     &                 DREAL(PSIS),DIMAG(PSIS),
     &                 (DREAL(AIF(K)),K=1,NCOIL),
     &                 (DIMAG(AIF(K)),K=1,NCOIL)
      ELSEIF (NCASE.EQ.4) THEN
         WRITE(CHTIME,90) DREAL(TMP),
     &                 DREAL(PSIS),DIMAG(PSIS),
     &                 (DREAL(AVF(K)),K=1,NCOIL),
     &                 (DIMAG(AVF(K)),K=1,NCOIL),
     &                 (DREAL(AIF(K)),K=1,NCOIL),
     &                 (DIMAG(AIF(K)),K=1,NCOIL),
     &                 DREAL(PSISNOISE),DIMAG(PSISNOISE)
      ELSEIF (NCASE.EQ.5) THEN
         WRITE(CHTIME,90) DREAL(TMP),
     &                 DREAL(PSISS),DIMAG(PSISS),
     &                 DREAL(PSIS),DIMAG(PSIS),
     &                 (DREAL(AIF(K)),K=1,NCOIL),
     &                 (DIMAG(AIF(K)),K=1,NCOIL)
      ELSEIF (NCASE.EQ.9.OR.NCASE.EQ.10) THEN
         IF (1.EQ.0) THEN
         WRITE(FILENAME,"(A5,I0.4,A4)") "SOLX_",NIT,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
         REWIND(CHOUTP)
         DO MS=1,MSMAX
         M = (MS-1)*NXCOMP
         DO K=1,NTP1
            WRITE(CHOUTP,111) ABS(X(KXB1+M,K)),ABS(X(KXV1+M,K)),
     &                        ABS(X(KXJ2U+M,K)),ABS(X(KXJ3+M,K)),
     &                        ABS(X(KXJ2L+M,K)),ABS(X(KXJRE+M,K))
         ENDDO
         ENDDO
         CLOSE(CHOUTP)
 111     FORMAT(6(1X,E8.1))
         ENDIF

         SOLMAXX = MAXVAL(ABS(X(:,2:NR)))
         SOLMAXY = MAXVAL(ABS(Y(:,2:NR-1)))
         WRITE(CHTIME,90) DREAL(TMP),
     &                 DREAL(TOT_RJA),
     &                 DREAL(TOT_JRE),DIMAG(TOT_JRE),
     &                 DREAL(TOT_JPA),DIMAG(TOT_JPA),
     &                 DREAL(TOT_EPA),
     &                 -RE_E0,SOLMAXX,SOLMAXY,
     &                 DREAL(COR_RJA),
     &                 DREAL(COR_JRE),DIMAG(COR_JRE),
     &                 DREAL(COR_JPA),DIMAG(COR_JPA),
     &                 DREAL(COR_EPA)
      ELSEIF (CALPHA7.GT.0.) THEN
         WRITE(CHTIME,91) NR,NCHI,NRATSURF,NPICK,
     &                    DREAL(TMP),DREAL(DENSPUMPA),
     &                    TTORQREY,TTORQJXB,TTORQNTV,
     &                    (DREAL(TPSI(J)),J=1,2*NRATSURF),
     &                    (DIMAG(TPSI(J)),J=1,2*NRATSURF),
     &                    (TROTM(J),J=1,NR),
     &                    (TORQUEJXB(J),J=1,NR),
     &                    (TORQUENTV(J),J=1,NR),
     &                    (TORQUEREY(J),J=1,NR),
     &                    (DREAL(DISPNORM(J)),J=1,NCHI),
     &                    (DIMAG(DISPNORM(J)),J=1,NCHI),
     &                    TDNTRMHD/TOTDENS,TDNTRNTV/TOTDENS,
     &                    (TRHOM(J),J=1,NR),
     &                    (DPTRANMHD(J),J=1,NR),
     &                    (DPTRANNTV(J),J=1,NR),
     &                    (DREAL(BZPICK(J)),J=1,NPICK),
     &                    (DIMAG(BZPICK(J)),J=1,NPICK),
     &                    DREAL(VDEXIR),DIMAG(VDEXIR),
     &                    DREAL(VDEXIZ),DIMAG(VDEXIZ)
      ENDIF
C    &                    (TORQUEERGO(J),J=1,NR),
C    &                    TDNTRERGO/TOTDENS,
C    &                    (DPTRANERGO(J),J=1,NR),
      CLOSE(CHTIME)
      OPEN(CHTIME,FILE='TIMEEVOL.OUT',POSITION='APPEND')

 90   FORMAT(3000(E14.5E3,1X))
 91   FORMAT(4(I4,1X),6000(E14.5E3,1X))

C     --- CHECK IF TO TURN ON THE FEEDBACK
      IF (ISWITCH.EQ.0.AND.ABS(PSIS).GE.THRESHOLD) THEN
         ISWITCH = 1
         IFLAG   = 1
      ENDIF 

      IF (ALLOCATED(TPSI))   DEALLOCATE(TPSI)
      IF (ALLOCATED(BZPICK)) DEALLOCATE(BZPICK)
      DEALLOCATE(AVF,AIF,PSIF,PSIFOLD)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------RESCALE THE EIGEN-SOLUTION FOR STARTING TIME EVOLUTION--------
C---------Y.Q.LIU, 19/08/2001-------------------------------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE RESCALEXY
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      INTEGER MS,I,J,K
      REAL*8    AMP
      COMPLEX*16 PSIS

C     --- CALCULATE THE SENSOR SIGNAL
      PSIS = 0.0
      J = NSENS
      IF (INCFEED.EQ.1) THEN
         DO K=1,NCOIL
            DO MS=1,MSMAX
               PSIS = PSIS + 
     &                BFEEDSS(K,MS)*B1U(ISENS(J)+NR,MS)
            ENDDO
         ENDDO
      ENDIF
      IF (INCFEED.EQ.11.OR.INCFEED.EQ.12) THEN
         DO K=1,NCOIL
            DO MS=1,MSMAX
               PSIS = PSIS + 
     &                BFEEDST(K,MS)*B1U(ISENS(J)+NR,MS)/CI
            ENDDO
         ENDDO
      ENDIF
      IF (INCFEED.EQ.3.OR.INCFEED.EQ.9) THEN
         DO K=1,NCOIL
            DO MS=1,MSMAX
               PSIS = PSIS + 
     &                BFEEDST(K,MS)*B2U(ISENS(J)+NR,MS)
            ENDDO
         ENDDO
      ENDIF

      IF (ABS(PSIS).LE.1e-16) PSIS=1.0

      AMP = 0.1*ABS(THRESHOLD)/ABS(PSIS)

      DO 100 MS=1,MSMAX
      DO  20 I=1,NTP1
      V1U(I,MS) = V1U(I,MS)*AMP
      X1U(I,MS) = X1U(I,MS)*AMP
      B1U(I,MS) = B1U(I,MS)*AMP
      J2U(I,MS) = J2U(I,MS)*AMP
      J3U(I,MS) = J3U(I,MS)*AMP
      J2L(I,MS) = J2L(I,MS)*AMP
      PDE(I,MS) = PDE(I,MS)*AMP
 20   CONTINUE
      DO  30 I=1,NTOT
      V2U(I,MS) = V2U(I,MS)*AMP
      V3U(I,MS) = V3U(I,MS)*AMP
      X2U(I,MS) = X2U(I,MS)*AMP
      B2U(I,MS) = B2U(I,MS)*AMP
      B3U(I,MS) = B3U(I,MS)*AMP
      J1U(I,MS) = J1U(I,MS)*AMP
      PRE(I,MS) = PRE(I,MS)*AMP
      PEE(I,MS) = PEE(I,MS)*AMP
      PPARA(I,MS) = PPARA(I,MS)*AMP
      PPERP(I,MS) = PPERP(I,MS)*AMP      
      RHOP(I,MS) = RHOP(I,MS)*AMP
 30   CONTINUE
 100  CONTINUE

      RETURN
      END
