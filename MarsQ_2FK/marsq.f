      PROGRAM MARS
C========================================================================
C== MODIFICATIONS FOR CRAY-2:
C== COMMENT @PROCESS
C== COMMENT IBM ERROR ROUTINE
C========================================================================
C MODIFICATIONS  FOR IBM/3090:
C                  @PROCESS DC(...)   ===> DYNAMIC COMMON DEFINITION
C
C MATRICES .SUBM ARE REQUIRE THE MAJOR PART OF MEMORY SPACE
C ON IBM/3090 MEMORY IN EXCESS OF 16M CAN BE ACHIEVED ONLY THROUGH
C THE DYNAMIC COMMON DEFINITION.
C THE MATRICES .SUBM THUS MUST BE DEFINED IN A COMMON (COMMON /MATRIX/)
C
C***********************************************************************
C
C      VERSION WITH RESISTIVE WALL.  SOLVE VACUUM EQUATIONS ON EACH
CMSC   ITERATION. THE EFFECT OF PLASMA ROTATION IS ASSUMED TO GIVE
CMSC   ONLY A LOCAL DOPPLER SHIFT TO THE PLASMA DYNAMIC RESPONSE
CMSC   THE SHEAR ROTATION PROFILE IS SPECIFIED BY ROTE
CMSC   THE ROTATION PROFILES ARE STORED IN ROT, ROTM
CMSC   THE DERIVATIVES OF THE ROTATION PROFILES ARE STORED IN DROT,DROTM
CMSC      THE SIGNS OF THE CORSS MATRIX ELEMENTS ARE CORRECTED ARTIFICIALLY
CMSC      7/8/94 (10 IN TOTAL)
C      DAMPING COEFFICIENT NU INTRODUCED FOR XI-PARALLEL
C
C      B-MATRICES ELIMINATED  01/08/89
C
C      MATRICES (ASUBM - HSUBM) DIMENSIONED IN MAIN
C      PROGRAM AND THEN TRANSFERED AS SUBROUTINE ARGUMENTS WITH VARIABLE
C      DIMENSION.
C
C      NEWPAMS VERSION OF COEFFICIENT MATRIX WITH SOLUTION SPLIT INTO
C      X'S (X1,V1,B1,J2U,J3,J2L,PDE,PED) AND
C      Y'S (X2,V2,V3,B2,B3,J1,PRE,PEE,PEP,PPARA,PPERP).
C          
C      PAMS SOLVES THE SYSTEM
C
C                         A XX = LAMBDA B XX
C
C                  WITH:  A      = BLOCK MATRIX (MHD OPERATOR)
C                         B      = BLOCK MATRIX
C                         XX     = EIGENVECTOR (XX=(X,Y))
C                         LAMBDA = EIGENVALUE
C
C    FORM OF A AND B MATRICES:             FORM OF XX VECTOR
C
C    ------------------                              ------
C    ! B1  ! E1  ! C1 !                              ! X1 !
C    !-----!-----!----!                              !----!
C    !     !     !    !                              !    !
C    ! F1  ! D1  ! G1 !                              ! Y1 !
C    !     !     !    !                              !    !
C    !-----!-----!----!------------                  !----!
C    ! A1  ! H1  ! B2 !  E2 !  C2 !                  ! X2 !
C    ------------!----!-----!-----!         *        !----!
C                !    !     !     !                  !    !
C                ! F2 !  D2 !  G2 !                  ! Y2 !
C                !    !     !     !                  !    !
C                !----!-----!-----!----              !----!
C                ! A2 !  H2 !  B3 ! ...              ! X3 !
C                !---------------------              !----!
C                            ..............            ..
C                            ..............            ..
C
C TO OBTAIN THE MHD OPERATOR WITH THE USUAL SIGN, THE COEFFICIENTS
C OF THE EQUATIONS FOR V1,V2,V3,B1,B2,B3,PRE HAVE BEEN MULTIPLIED
C BY -1 (EXCEPT THE TERMS ARISING FORM THE TIME DERIVATIVES!)
C
C NOW LAMBDA = - IMM(W) + (0.,1.) REAL(W)
C
C-----------------------------------------------------------------------
C
C MARS version modified to seek critical Gain
C 
C AUG 05, 2003, YQLIU: 
C ADDED STUFF FOR POLOIDAL CONVERGENCE IMPROVEMENT 
C-----------------------------------------------------------------------
C IMPORTANT NOTES!!!!!
C 1. DON'T USE AL0 DIRECTLY INSIDE PROCEURES COEFFI(...) OR VACCOE(...)
C BECAUSE CONJG(...) ARE TAKEN FOR ALL COEFFICIENTS FOR M+K COMPONENTS 
C ONE MUST USE SHIFT TO REPLACE AL0 IN A PROPER WAY
C IT TOOK ME THREE DAYS TO FIGURE THIS PROBLEM OUT
C MARCH 09, 2004, YQLIU
C
C 2. COMPUTATIONAL DOMAIN SHOULD NOT EXTEND TO R<0 REGION, BECAUSE JACOBIAN 
C VANISHES WITH R=0.
C-----------------------------------------------------------------------
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE RESMATM
      USE FEEDBACKM
      USE ADAPTIVEM
      USE REORBITM
      USE MPIENV
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'
C
      INTEGER MF,I,J,K,L,INCKRE
      REAL*8  GAINAI(50,NCOIL0), GAINPI(50,NCOIL0), ALCR, ALCR0
      REAL*8  NEWTON_EPSILONG
      COMPLEX*16 ALNORM1,ALNORM2,ALFUN1,ALFUN2,ALTMP
      COMPLEX*16,DIMENSION(:),ALLOCATABLE:: WORKC
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::SHIFTC,SHIFTM
C
      ALCR  = 2.0E-06
      ALCR0 = 1.0E-02

      OPEN(CHNAME,FILE='RUN.IN',FORM='FORMATTED',STATUS='OLD')
      REWIND(CHNAME)
C
      CALL PRESET
      CALL FEEDINI

      NPARAM = 0
10    NPARAM = NPARAM + 1
      CALL RDNAME
      
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("MXMAX,MYMAX= ",2I5)')MXMAX,MYMAX
      ENDIF
      
      IF (IGO.EQ.0) GOTO 900

      IF (NPARAM.EQ.1) CALL READTOR

      CALL COUPLE
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER COUPLE")')
      ENDIF

C     INITIAL TIME FOR NCASE=6 OR 5
      TOTTIME = 0.0
      ALNORM1 = 0.0
      ALNORM2 = 0.0
      ALFUN1  = 0.0
      ALFUN2  = 0.0
      AL0     = CALPHA1/CALPHA2
      ALNORM  = AL0 

      IF (NPARAM.EQ.1)  CALL READVAC
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER READVAC")')
      ENDIF

      CALL COTROL
      
      DO 100 ISWEEP=1,NSWEEP
      DO 100 IADAPS=0,NADAPS

      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*) 'ISWEEP=',ISWEEP,'  IADAPS=',IADAPS
      ENDIF

      IF (ISWEEP.GT.1.AND.NADAPS.EQ.0) GOTO 200

      CALL GETADAPS

      IF (ABS(RNTOR).LT.0.1.AND.NCASE.EQ.10.AND.NWALL.GT.0
     &    .AND.MWALL.GT.0) CALL UPDATE_WALL(ISWEEP)

      CALL CLEAR
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER CLEAR")')
      ENDIF

      CALL Fourier_Equilibrium
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER Fourier_Equilibrium")')
      ENDIF

      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
        IF (NPARAM.EQ.1) CALL PRINTRMZM
      ENDIF
      
      IF (NCASE.EQ.11) THEN
        STOP "GRID COMPUTED, EXITING.."
      ENDIF

      IF (SLEFT.LT.1.e-10) CALL FIXORI
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER FIXORI")')
      ENDIF

      CALL GEOMET
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER GEOMET")')
      ENDIF

      INCKRE = 0
      IF (INCKIN.GT.0) THEN
         DO K=1,NSPECIES
            IF (ISPECIES_F0(K).EQ.5) INCKRE = 1
            IF (ISPECIES_F0(K).EQ.6) INCKRE = 2
         ENDDO
      ENDIF   
      IF (KXJRE.GT.0.OR.KJRER.EQ.6.OR.INCKRE.GT.1) CALL GET_RE_CONST
      IF (INCKRE.EQ.2.AND.RE_PMAX.LT.0.) CALL GET_RE_PMAX

      CALL PROFIL
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER PROFIL")')
      ENDIF

      IF (NPICK.GT.0) CALL CALC_PICK_DATA
 
      CALL PATCH
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER PATCH")')
      ENDIF
C
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      IF (NPARAM.EQ.1.AND.ISWEEP.EQ.1.AND.IADAPS.EQ.0.AND.1.EQ.0) 
     &   CALL OutputT7
      ENDIF

      MF = MXMAX
      IF (MYMAX.LT.MXMAX) MF=MYMAX
      MF = 2 * MF * MF
      IF (MF.LT.MEDIM*4) MF=MEDIM*4
      
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("MF= ",I9)')MF
      ENDIF
C
      IF (NPARAM.EQ.1) CALL GCONTR(MF)
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER GCONTR")')
      ENDIF

C
      IF (NCASE.EQ.4) THEN
         DO K=1,NCOIL
            GAINA(K)  = 0.0
         ENDDO
      ENDIF

      IF (NCASE.EQ.3) THEN
         DO K=1,NCOIL
            GAINA(K)  = 0.0
            HGAINA(K) = 0.0
         ENDDO
      ENDIF

      IF (NCASE.EQ.4) THEN
      DO K=1,NCOIL
      GAINAI(1,K) = GAINA(K)
      GAINPI(1,K) = GAINP(K)
      DO L=2,NSWEEP
         GAINAI(L,K)=GAINAI(L-1,K) + HGAINA(K)
         GAINPI(L,K)=GAINPI(L-1,K) + HGAINP(K)
      ENDDO
      ENDDO
      ENDIF
         
      IF (NCASE.GE.3.AND.NCASE.LE.10.AND.
     &    (ISMPIRUN.EQ.0.OR.RANK.EQ.ROOT)) THEN
         OPEN(CHTIME,FILE='TIMEEVOL.OUT')
         REWIND(CHTIME)
      ENDIF    

      IF (NCASE.EQ.4) THEN
         DO K=1,NCOIL
            GAINAI(3,K) = GAINAI(2,K)
            GAINPI(3,K) = GAINPI(2,K)
            GAINAI(2,K) = 0.0
            GAINPI(2,K) = 0.0
         ENDDO
      ENDIF    

      IF (NPARAM.EQ.1.OR..NOT.KEEPMO) CALL INITPE
      
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*) ' AFTER INITPE'
      ENDIF
  
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'(" MXMAX,MYMAX = " , 2I5)') MXMAX,MYMAX
      ENDIF
      
      IF (.NOT. ALLOCATED (ASUBM)) THEN
        IF (ISMPIRUN.EQ.0.OR.RANK.EQ.ROOT) THEN
             ALLOCATE(                   ASUBM(MXMAX,MXMAX,NTP1+2),
     &                                   BSUBM(MXMAX,MXMAX,NTP1+3),
     &                                   CSUBM(MXMAX,MXMAX,NTP1+2),
     &                                   DSUBM(MYMAX,MYMAX,NTP1+2),
     &                                   ESUBM(MXMAX,MYMAX,NTP1+2),
     &                                   FSUBM(MYMAX,MXMAX,NTP1+2),
     &                                   GSUBM(MYMAX,MXMAX,NTP1+2),
     &                                   HSUBM(MXMAX,MYMAX,NTP1+2))
C
             ASUBM  = (0.0,0.0)
             BSUBM  = (0.0,0.0)
             CSUBM  = (0.0,0.0)
             DSUBM  = (0.0,0.0)
             ESUBM  = (0.0,0.0)
             FSUBM  = (0.0,0.0)
             GSUBM  = (0.0,0.0)
             HSUBM  = (0.0,0.0)
C
        ELSEIF (ISMPIRUN.GE.1.AND.RANK.NE.ROOT) THEN
C            FAKE MATRIX FOR KJP IN SLAVE PROCESS CALCULATING KINETIC EFFECT ON EACH FLUX SURFACE        
             ALLOCATE(                   ASUBM(2,2,2),
     &                                   BSUBM(2,2,2),
     &                                   CSUBM(2,2,2),
     &                                   DSUBM(2,2,2),
     &                                   ESUBM(2,2,2),
     &                                   FSUBM(2,2,2),
     &                                   GSUBM(2,2,2),
     &                                   HSUBM(2,2,2),
     &                                   SHIFTC(2),SHIFTM(2))
             ASUBM  = (0.0,0.0)
             BSUBM  = (0.0,0.0)
             CSUBM  = (0.0,0.0)
             DSUBM  = (0.0,0.0)
             ESUBM  = (0.0,0.0)
             FSUBM  = (0.0,0.0)
             GSUBM  = (0.0,0.0)
             HSUBM  = (0.0,0.0)
             SHIFTC = (0.0,0.0)
             SHIFTM = (0.0,0.0)
        ELSE
             PRINT *,"UNKOWN RUNNING MODE"
             STOP "UNKOWN RUNNING MODE"
        END IF
C
        ALLOCATE( X(NXCOMP*MSMAX,NTP1), Y(NYCOMP*MSMAX,NTP1))
        ALLOCATE(DX(NXCOMP*MSMAX,NTP1),DY(NYCOMP*MSMAX,NTP1))
        X        = (0.0,0.0)
        Y        = (0.0,0.0)
        DX       = (0.0,0.0)
        DY       = (0.0,0.0)
      END IF
C
      IWORK=2*(3*NTP1*(MXMAX+MYMAX)+10*MXMAX**2+2*MXMAX*MYMAX
     &        +4*MYMAX**2+10*MXMAX+7*MYMAX)

      IF (.NOT. ALLOCATED (WORK)) THEN
        IF (ISMPIRUN .EQ. 0 .OR. RANK.EQ.ROOT) THEN
            ALLOCATE(WORK(IWORK))
        END IF
      END IF

      IF (.NOT. ALLOCATED (WORKC)) THEN
        IF (ISMPIRUN .EQ. 0 .OR. RANK.EQ.ROOT) THEN
            ALLOCATE(WORKC(IWORK))
        END IF
      END IF
C
      IF (ISMPIRUN .EQ. 0 .OR. RANK.EQ.ROOT) THEN
         WORK=0.
         WORKC=0.
      END IF

 200  CONTINUE

      CALL COEFFMOMENT
      CALL COEFFDNTRAN

C     READ SAVED SOLUTION FOR NCASE=6 RUNS
      IF (NCASE.EQ.6.AND.KSOLREAD.EQ.1.AND.ISWEEP.EQ.1) THEN
         IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
         WRITE(*,*) 'READ SAVED SOLUTION FROM PREVIOUS RUN...'
         ENDIF
         CALL SOLREAD(X,Y,CALPHA1,CALPHA5,AL0,ALNORM,TOTTIME)
         CALL GETXY(X,Y)
         CALL COEFFMOMENT
         CALL COEFFDNTRAN
      ENDIF

C     READ SAVED SOLUTION FOR NCASE=1 RUNS
      IF ((NCASE.EQ.1.OR.NCASE.EQ.2.OR.NCASE.EQ.9.OR.NCASE.EQ.10).AND.
     &    KSOLREAD.EQ.1.AND.ISWEEP.EQ.1) THEN
         IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
         WRITE(*,*) 'READ SAVED SOLUTION FROM PREVIOUS RUN...'
         ENDIF
         CALL SOLREAD(X,Y,CALPHA1,CALPHA5,AL0,ALNORM,TOTTIME)
         CALL GETXY(X,Y)
      ENDIF

C     UPDATE EQUILIBRIUM ROTATION AND DENSITY PROFILES
      IF ((NCASE.EQ.6.OR.NCASE.EQ.7).AND.ABS(CALPHA2-CALPHA3).LT.1.0E-3)
     &   CALL QLIN_UPDATE

      CALL FTCOEFF
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
         WRITE(*,*) ' AFTER FTCOEF'
      ENDIF

      CALL CALC_NUSTAR

      IF (ISWEEP.GT.1) NPARAM = NPARAM + 1
      
      ICALPHA2 = 1
      IF (NCASE.EQ.6.AND.ISWEEP.EQ.1) ICALPHA2 = 0
      IF (NCASE.EQ.10.AND.ISWEEP.EQ.1) ICALPHA2 = 0

      IF (NCASE.EQ.4.AND.ISWEEP.GT.1) THEN
         CALPHA1 = 10.*DREAL(ALNORM)
         CALPHA2 = 1.
         CALPHA4 = (1.,0.)
         CALPHA7 = 0.
         INORMSOL= 0
         AL0     = CALPHA1/CALPHA2
      ENDIF
      IF (NCASE.EQ.4.AND.ISWEEP.EQ.2) CALL RESCALEXY
      CALL INITXY(X,Y)
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
         WRITE(*,*) ' AFTER INITXY'
      ENDIF
      
      IF (NCASE.EQ.6.OR.NCASE.EQ.10) TOTTIME = TOTTIME+1.0/REAL(CALPHA1)

C.....ERROR FIELD AMPLIFICATION
      IF (INCFEED.EQ.8.OR.INCFEED.EQ.18) THEN
         DO K=1,NCOIL
            GAINA(K) = 0.0
         ENDDO
      ENDIF

      IF (ISWEEP.GT.1.AND.NCASE.EQ.4) THEN
         DO K=1,NCOIL
            GAINA(K) = GAINAI(ISWEEP,K)
            GAINP(K) = GAINPI(ISWEEP,K)
         ENDDO
      ENDIF

      CALL FEEDCTRL
      IF (INCFEED.GT.0) THEN
         CALL GETRFK
         CALL CALCFEED
         CALL CALDFEED
         CALL CALAFEED
         IF (NCASE.NE.3.AND.NCASE.NE.4) CALL SETJINI
      ENDIF   

      KJPKEY = 1
      KPBKEY = 1
      IF (ISWEEP.GT.1) KJPKEY = 2
      IF (IPERTURB.NE.0) KPBKEY=0

      ALNORM2 = ALNORM1
      ALFUN2  = ALFUN1
      ALNORM1 = ALNORM

C CALL KJP IN SLAVE PROCESS
      IF ((ISMPIRUN.EQ.1.OR.ISMPIRUN.EQ.3).AND.RANK.NE.ROOT) THEN
         CALL MPI_OPEN_FILE(RANK)
         WRITE(CHMPI,*) 'NPARAM,ISWEEP,RANK=',NPARAM,ISWEEP,RANK
         CALL MPI_CLOSE_FILE(RANK) 
         CALL KJP(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM,
     &             SHIFTC,SHIFTM,2,2,2)
      ENDIF

      IF (ISMPIRUN.EQ.0.OR.RANK.EQ.ROOT) THEN

      CALL LINEAR(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
      
      WRITE(*,*) ' AFTER LINEAR'

      IF (MSMAX.LE.7) 
     &   CALL PRINTATOH(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM, 
     &                  HSUBM,MXMAX,MYMAX,NTP1,1,1,JSOUT)
 
C
      IF (ABS(NCOUPL).EQ.1.OR.NCOUPL.EQ.4.OR.INCFEED.EQ.20.OR.
     &    INCFEED.EQ.21.OR.INCFEED.EQ.22) THEN
         OPEN(CHOUTP,FILE='BNORM01.IN',FORM='FORMATTED')
         REWIND(CHOUTP)
      ENDIF

C     CHECK MOMENTUM BALANCE
      IF (KSOLTEST.EQ.1) THEN 
         CALL SOLTEST(MXMAX,MYMAX,NTOT,NXCOMP,NYCOMP,
     &        ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM,X,Y)
         GOTO 900
      ENDIF

      ALAM = CALPHA2
      IF (ISMPIRUN.EQ.0.OR.RANK.EQ.ROOT) CALL CALPAM(
     $      MXMAX,MYMAX,NTOT,NXCOMP,NYCOMP,NCASE,NITMAX
     $     ,EPSPAM,EPSDET
     $     ,AL0,ALAM,ALNORM,NONCON
     $     ,ASUBM(1,1,2)
     $     ,BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM
     $     ,HSUBM(1,1,2)
     $     ,DX,DY, X, Y, WORKC, IWORK )
      WRITE(*,*)  ' AFTER PAMS'

C     NEWTON SOR ITERATION
      ALFUN1 = ALNORM - ALNORM1
      IF (ISWEEP.LE.NSWEEP.AND.ISWEEP.GT.1.AND.NCASE.LE.1) THEN
         IF (ABS(ALNORM1).GT.1E-6) THEN
            NEWTON_EPSILON = ABS(ALNORM1)*1E-4
         ELSE
            NEWTON_EPSILON = 1E-6*1E-4
         ENDIF
         IF (ISWEEP.GE.2.AND.INCKIN.EQ.1.AND.
     $      ABS(ALFUN1-ALFUN2).GT.NEWTON_EPSILON) THEN
C           IF (ISWEEP .EQ. 2) THEN
C              ALTMP   = ALNORM1
C              ALNORM1 = ALNORM2
C              ALNORM2 = ALTMP
C              ALTMP   = ALFUN1
C              ALFUN1  = ALFUN2
C              ALFUN2  = ALTMP
C           ENDIF
            ALNORM = ALNORM1 
     $             - ALFUN1*(ALNORM1-ALNORM2)/(ALFUN1-ALFUN2)
            WRITE(*,*)'NEWTON SOR'
            WRITE(*,*)'ISWEEP=',ISWEEP
            WRITE(*,*)'ALNORM=',ALNORM,'ALNORM1=',ALNORM1
            WRITE(*,*)'ALNORM2=',ALNORM2
            WRITE(*,*)'ALFUN1=',ALFUN1,'ALFUN2=',ALFUN2
         ELSE
            WRITE(*,*) 'SOR'
         ENDIF
         IF (ISWEEP .EQ. NSWEEP .AND. INCKIN.NE.1) THEN
            ALNORM = ALNORM
         ELSE
            ALNORM = ALTAU*ALNORM + (1.0-ALTAU)*ALNORM1
         ENDIF
      ENDIF
      
      IF ((NCASE.EQ.6.OR.NCASE.EQ.10).AND.
     &    (ISWEEP.GT.1.OR.KSOLREAD.EQ.1)) THEN
         IF (NONCON.EQ.1) CALPHA1 = CALPHA1/CALPHA8
         IF (NONCON.EQ.2) CALPHA1 = CALPHA1*CALPHA8
         IF (ABS(1./CALPHA1).GT.TDELTALIM) CALPHA1=1./TDELTALIM
         IF (ABS(1./CALPHA1).LT.TDELTALOW) CALPHA1=1./TDELTALOW
         CALPHA5 = CALPHA1
         AL0     = CALPHA1/CALPHA2
         WRITE(*,*) 'CALPHA1=',CALPHA1
         IF ((NONCON.EQ.1.OR.NONCON.EQ.2).AND.NCASE.EQ.6) THEN
            CALL COEFFMOMENT
            CALL COEFFDNTRAN
         ENDIF
      ENDIF

      IF (ABS(NCOUPL).EQ.1.OR.NCOUPL.EQ.4.OR.INCFEED.EQ.20.OR.
     &    INCFEED.EQ.21.OR.INCFEED.EQ.22) CLOSE(CHOUTP)

 50   CONTINUE
C
C     COPY PAMS VARIABLES AND OUTPUT DATA
C
      CALL GETXY(X,Y)
     
      CALL GETXYMORE

      CALL PLOTLP(X,Y)
      IF (NV.GE.2.AND.NWALL.GT.0) THEN
         DO J=1,NWALL 
            IWALLJ = IWALL(J)
            TAUWJ  = TAUW(J)
            CALL DPWALL
         ENDDO
      ENDIF

      CALL LISTMO(X,Y)
C.....TEST FOR VACUUM-PLASMA BOUNDARY
C     CALL TESTV
C
C     PHYSICS DIGNOSTICS
C
      IF (DCONTI.AND.KYV3.GT.0) CALL LOCONT(Y)
      IF (DTERMS) CALL TERMS(
     &     DX, DY, X, Y, WORK, IWORK,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
C
      IF (NCASE.EQ.1.OR.NCASE.EQ.2.OR.NCASE.EQ.6.OR.NCASE.EQ.10) 
     &   CALL OUTPUT_RECTRZ(ISWEEP)

      VDEXIR = (0.,0.)
      VDEXIZ = (0.,0.)
      IF (ABS(RNTOR).LT.0.1) CALL OUTPUT_XIRZ(ISWEEP)

      IF (NOPP.GT.0.AND.ODJPHI.GT.0.0) CALL OUTPUT_BFILAMENT

      IF (NCASE.EQ.2.AND.INCFEED.EQ.4.AND.KKF.EQ.-3)
     &   CALL BS_ESC

      IF (NCASE.EQ.2.AND.INCFEED.EQ.4.AND.KKF.EQ.-1.AND.NSENS.GT.0)
     &   CALL BS_B1

      CALL OUTPUT(ISWEEP)

      IF (ISWEEP.EQ.1) ALCR0=DREAL(ALNORM)

      IF ((NCASE.EQ.0.OR.NCASE.EQ.1.OR.NCASE.EQ.2.OR.NCASE.EQ.6.OR.
     &     NCASE.EQ.9.OR.NCASE.EQ.10).AND.ISWEEP.EQ.NSWEEP) THEN 
C        SAVE SOLUTION BEFORE CALLING ENERGYMAT, WHICH CHANGES IT
         IF ((NCASE.EQ.1.OR.NCASE.EQ.2.OR.NCASE.EQ.6.OR.NCASE.EQ.9
     &        .OR.NCASE.EQ.10).AND.KSOLSAVE.EQ.1) THEN
            WRITE(*,*) 'SAVE SOLUTION TO BE USED FOR NEXT RUN.'
            CALL SOLSAVE(X,Y,CALPHA1,CALPHA5,AL0,ALNORM,TOTTIME)
         ENDIF

         IF (INCKIN.GT.0.AND.NCASE.NE.6.AND.1.EQ.0) CALL KDWKDENSITY
         IF (NCASE.NE.6.AND.NCASE.NE.10.AND.KEFORM.NE.0) THEN
            CALL ENERGYMAT(X,Y,
     &      ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
            WRITE(*,*) 'AFTER ENERGYMAT'
         ENDIF
      ENDIF

      ENDIF
 100  CONTINUE 

      IF (NREORBIT.GT.0) CALL RE_TRACING

      IF (ISMPIRUN.EQ.0.OR.RANK.EQ.ROOT) THEN
         IF (NCASE.GE.3.AND.NCASE.LE.10) CLOSE(CHTIME)
         IF (NCASE.EQ.0.OR.NCASE.EQ.1)  AL0 = ALNORM
      ENDIF

      CALL MPI_EXIT
      GOTO 10

      OPEN (101,FILE='ENERGY_ANALYZE_INPUT.OUT',FORM='FORMATTED')
      WRITE (101,901) RNTOR, GAMMA,DREAL(ALNORM),DIMAG(ALNORM)
 901  FORMAT(E20.10,5X,E20.10,5X,E30.15,5X,E30.15)
      CLOSE(101)

 900  CONTINUE

      STOP 'NORMAL'
      END
*DECK COUPLE
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C----------COUPLE ROUTINE-------------G.VLAD 29/03/1989-----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.......CONSTRUCT COUPLING CONTROLS.....................................
C
      SUBROUTINE COUPLE
C     =================
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'comioc.inc'
      INTEGER   MSA,MSB,NSA,NSB,MSAMAX,MSBMAX,MPL,MMI,NPL,NMI,
     I          MSPL,MSMI,NSPL,NSMI,NSCH,MSPMAX,MSMMAX,MS,NS,NSEQ,NSPE
      PARAMETER (NSEQ = 1, NSPE = 2)
C
C     ALLOCATE VARIABLES
C
      IF (.NOT. ALLOCATED (RM)) THEN
         ALLOCATE(RM(MEDIM,NSMAX),IM(MEDIM,NSMAX))
         ALLOCATE( MPLUS(MEDIM,NSMAX,MEDIM,NSMAX))
         ALLOCATE(MMINUS(MEDIM,NSMAX,MEDIM,NSMAX))
      END IF
C
C..INITIALIZATION
C
      DO 20 NSA=1,NSMAX
      DO 20 NSB=1,NSMAX
C
      NPLUS(NSA,NSB)  = 0
      NMINUS(NSA,NSB) = 0
C
      DO 10 MSA=1,MEDIM
      DO 10 MSB=1,MEDIM
C
      MPLUS(MSA,NSA,MSB,NSB)  = 0
      MMINUS(MSA,NSA,MSB,NSB) = 0
C
 10   CONTINUE
 20   CONTINUE
C
C
      DO 30 MS=1,MEDIM
      IM(MS,NSEQ) = MS-1
 30   CONTINUE
C
      RN(NSPE) = RNTOR
      IF (M2-M1.NE.MSMAX-1) STOP 'M2-M1'
      DO 40 MS = 1,MSMAX
 40   IM(MS,NSPE) = M1 + MS-1
C
      DO 50 MS=1,MEDIM
      RM(MS,1) = DFLOAT(MS-1)
      RN(1) = 0.
 50   CONTINUE
C
      NS=2
      DO 60 MS=1,MSMAX
      RM(MS,NS)=DFLOAT(IM(MS,NS))
 60   CONTINUE
C
C     MODIFY THIS TO MAKE RNTOR=0 WORK
C     YQL, 2007-10-01
C
      NSA = 2
      NSB = 1
C
      MSAMAX=MSMAX
      DO 180 MSA=1,MSAMAX
C
      MSBMAX=MSMAX
      DO 170 MSB=1,MSBMAX
C
      MPL=IM(MSA,NSA)+IM(MSB,NSB)
      MMI=IM(MSA,NSA)-IM(MSB,NSB)
C
      MSPMAX=MSMAX
      DO 80 MSPL=1,MSPMAX
      IF(IM(MSPL,NSA).EQ.MPL) GOTO 100
C
 80   CONTINUE
C
C
      GOTO 110
C
 100  MPLUS(MSA,NSA,MSB,NSB)=MSPL
C
 110  CONTINUE
C
      MSMMAX=MSMAX
      DO 120 MSMI=1,MSMMAX
      IF (IM(MSMI,NSA).EQ.MMI) GOTO 140
C
 120  CONTINUE
C
      GOTO 150
C
 140  MMINUS(MSA,NSA,MSB,NSB)=MSMI
C
 150  CONTINUE
C
C..TO AVOID DOUBLE COUNTING OF THE MODES
C
 170   CONTINUE
 180   CONTINUE
C
      RETURN
      END
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c------- Calculate Fourier Componets -------- D. Liu 31.07.97 ----------
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      Subroutine Fourier_Equilibrium
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE GIJLM
      USE RESMATM
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      real*8,dimension(:,:),allocatable::RW1,RW2,RW3,RW4,RW5,RW10

      ALLOCATE( RW1(nrp1,nchi), RW2(nrp1,nchi), RW3(nrp1,nchi),
     $          RW4(nrp1,nchi), RW5(nrp1,nchi), RW10(nrp1,nchi))
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'Equilibrium'

c     ------------------------------------------------------E Vac
      do j=1,nchi
         do i=1,nrp1
            rw1(i,j)=Req(i,j) 
            rw2(i,j)=Zeq(i,j)
            rw5(i,j)= DPsiDs(i)*PPeq(i)
         end do
         do i=1,nr
            rw3(i,j)=Reqm(i,j) 
            rw4(i,j)=Zeqm(i,j)
         end do
         rw3(nrp1,j) = rw3(nr,j)
         rw4(nrp1,j) = rw4(nr,j)
         do i = 1,nr
            rw10(i,j)= DPsiDsM(i)*PPeqM(i)
         end do
      end do
C
      NPSTRT    =  1
      call FFTDRIVER( RW1,  RPF,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RPF   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW2,  ZPF,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: ZPF   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW3,  RPFM,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)

      call FFTDRIVER( RW4,  ZPFM,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)

      call FFTDRIVER( RW5,  DPeDs,  FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)

      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for  DPeDs  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RJA,  JACOBI, FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for JACOBI  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW10, DPeDsM, FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DPeDsM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      5,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RJAM, JACOBM,  FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for JACOBM in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  RPF,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'RPF')
        call FFTOUTPT(RW2,  ZPF,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'ZPF')
        call FFTOUTPT(RW3,  RPFM,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'RPFM')
        call FFTOUTPT(RW4,  ZPFM,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'ZPFM')
        call FFTOUTPT(RW5,  DPeDs,  NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'DPeDs')
        call FFTOUTPT(RJA,  JACOBI, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'JACOBI')
        call FFTOUTPT(RW10, DPeDsM, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'DPeDsM')
        call FFTOUTPT(RJAM, JACOBM, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,  NCHI,  KUFFTP, 'JACOBM')
      ENDIF

      do j=1,nchi
         rw1(1,j) = 0.
         rw2(1,j) = 0.
         rw3(1,j) = 0.
         rw4(1,j) = 0.
c 
         do i=2,nrp1
            rw1(i,j)=G11L(i,j) / RJa(i,j)
            rw2(i,j)=G22L(i,j) / RJa(i,j)
            rw3(i,j)=G33L(i,j) / RJa(i,j)
            rw4(i,j)=G12L(i,j) / RJa(i,j)
         end do
      end do
      
C
      NPSTRT    =  1
      call FFTDRIVER( RW1,  DG11L,  FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG11L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      7,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW2,  DG22L,  FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG22L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW3,  DG33L,  FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG33L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      9,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW4,  DG12L,  FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG12L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     10,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  DG11L,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG11L')
        call FFTOUTPT(RW2,  DG22L,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG22L')
        call FFTOUTPT(RW3,  DG33L,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG33L')
        call FFTOUTPT(RW4,  DG12L,    NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG12L')
      ENDIF
c
      do j=1,nchi
         do i=1,nr
            rw1(i,j)=G11LM(i,j) / RJaM(i,j)
            rw2(i,j)=G22LM(i,j) / RJaM(i,j)
            rw3(i,j)=G33LM(i,j) / RJaM(i,j)
            rw4(i,j)=G12LM(i,j) / RJaM(i,j)
         end do
      end do
C
      NPSTRT    =  1
      call FFTDRIVER( RW1,  DG11LM, FORWD, NRP1,  NR,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG11LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     15,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW2,  DG22LM, FORWD, NRP1,  NR,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG22LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     16,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW3,  DG33LM, FORWD, NRP1,  NR,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG33LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     17,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW4,  DG12LM, FORWD, NRP1,  NR,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error  for DG12LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     18,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  DG11LM,   NRP1,  NR,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG11LM')
        call FFTOUTPT(RW2,  DG22LM,   NRP1,  NR,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG22LM')
        call FFTOUTPT(RW3,  DG33LM,   NRP1,  NR,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG33LM')
        call FFTOUTPT(RW4,  DG12LM,   NRP1,  NR,   NPSTRT
     &                     ,MEDIM,    NCHI,  KUFFTP, 'DG12LM')
      ENDIF
C
      DEALLOCATE( RW1,RW2,RW3,RW4,RW5,RW10)
C
      End   !!!---{Fourier_Equilibrium}
      
c     _________________________________________________________________
c     _________________________________________________________________

      Subroutine Out_Fourier_New
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      INCLUDE 'comioc.inc'
      if(1.eq.11) then
      OPEN(FOURIER_OUT,FILE='fourier_new.OUT',FORM='FORMATTED')
      REWIND(FOURIER_OUT)

      WRITE(FOURIER_OUT,1000) NRP1,MEDIM,NSMAX


      CALL GENOUT(DPeDs,NRP1,' DPeDs', FOURIER_OUT,NRP1 ,MEDIM,1)
      CALL GENOUT(JACOBI,NRP1,'JACOBI',FOURIER_OUT,NRP1 ,MEDIM,1)
      CALL GENOUT(DG11L,NRP1,' DG11L', FOURIER_OUT,NRP1 ,MEDIM,1)
      CALL GENOUT(DG22L,NRP1,' DG22L', FOURIER_OUT,NRP1 ,MEDIM,1)
      CALL GENOUT(DG33L,NRP1,' DG33L', FOURIER_OUT,NRP1, MEDIM,1)
      CALL GENOUT(DG12L,NRP1,' DG12L', FOURIER_OUT,NRP1, MEDIM,1)

      CALL GENOUT(DPeDsM,NRP1,'DPeDsM',FOURIER_OUT,NR ,MEDIM,1)
      CALL GENOUT(JACOBM,NRP1,'JACOBM',FOURIER_OUT,NR ,MEDIM,1)
      CALL GENOUT(DG11LM,NRP1,'DG11LM',FOURIER_OUT,NR ,MEDIM,1)
      CALL GENOUT(DG22LM,NRP1,'DG22LM',FOURIER_OUT,NR ,MEDIM,1)
      CALL GENOUT(DG33LM,NRP1,'DG33LM',FOURIER_OUT,NR, MEDIM,1)
      CALL GENOUT(DG12LM,NRP1,'DG12LM',FOURIER_OUT,NR, MEDIM,1)

      CLOSE(FOURIER_OUT)
      endif

 1000 FORMAT(3I20)
      
      End   !!!---{Out_Fourier_New}
      
      
c     _________________________________________________________________

      Subroutine Fourier_Vacuum
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE GIJLM
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      real*8,dimension(:,:),allocatable::RW1,RW2,RW3,RW4,RW6,RW7,
     &     RW8,RW9
      integer i,j
C
      ALLOCATE(RW1(nveq1,nchi),RW2(nveq1,nchi),RW3(nveq1,nchi),
     $                         RW4(nveq1,nchi),RW6(nveq1,nchi),
     $         RW7(nveq1,nchi),RW8(nveq1,nchi),RW9(nveq1,nchi))
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'Vacuum'
C
      do j=1,nchi
         do i=1,nveq1
            rw1(i,j)=VRG11L(i,j)/VRJA(i,j)
            rw2(i,j)=VRG22L(i,j)/VRJA(i,j)
            rw3(i,j)=VRG33L(i,j)/VRJA(i,j)
            rw4(i,j)=VRG12L(i,j)/VRJA(i,j)
         end do
         do i = 1,nveq
            rw6(i,j)=VRG11LM(i,j)/VRJAM(i,j)
            rw7(i,j)=VRG22LM(i,j)/VRJAM(i,j)
            rw8(i,j)=VRG33LM(i,j)/VRJAM(i,j)
            rw9(i,j)=VRG12LM(i,j)/VRJAM(i,j)
         end do
      end do
c
      IF (.NOT. ALLOCATED(VDG11L)) THEN
         ALLOCATE( VDG11L(nveq1,MEDIM), VDG22L(nveq1,MEDIM),
     $             VDG33L(nveq1,MEDIM), VDG12L(nveq1,MEDIM))
         ALLOCATE(VDG11LM(nveq1,MEDIM),VDG22LM(nveq1,MEDIM),
     $            VDG33LM(nveq1,MEDIM),VDG12LM(nveq1,MEDIM))
      END IF
C
C
      NVSTRT    =  1
      call FFTDRIVER( RW1,  VDG11L, FORWD, NVEQ1, NVEQ1,  NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG11L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW2,  VDG22L, FORWD, NVEQ1, NVEQ1,  NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG22L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW3,  VDG33L, FORWD, NVEQ1, NVEQ1,  NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG33L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW4,  VDG12L, FORWD, NVEQ1, NVEQ1,  NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG12L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW6,  VDG11LM,FORWD, NVEQ1, NVEQ,   NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG11LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      5,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW7,  VDG22LM,FORWD, NVEQ1, NVEQ,   NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG22LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW8,  VDG33LM, FORWD,NVEQ1, NVEQ,   NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG33LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      7,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW9,  VDG12LM,FORWD, NVEQ1, NVEQ,   NVSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VDG12LM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  VDG11L, NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG11L')
        call FFTOUTPT(RW2,  VDG22L, NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG22L')
        call FFTOUTPT(RW3,  VDG33L, NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG33L')
        call FFTOUTPT(RW4,  VDG12L, NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG12L')
        call FFTOUTPT(RW6,  VDG11LM,NVEQ1,  NVEQ,    NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG11LM')
        call FFTOUTPT(RW7,  VDG22LM,NVEQ1,  NVEQ,    NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG22LM')
        call FFTOUTPT(RW8,  VDG33LM,NVEQ1,  NVEQ,    NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG33LM')
        call FFTOUTPT(RW9,  VDG12LM,NVEQ1,  NVEQ,    NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'VDG12LM')
      ENDIF

      do j=1,nchi
         do i=1,nveq1
            rw1(i,j)=VRR(i,j)
            rw2(i,j)=VRZ(i,j)
            rw3(i,j)=VRRM(i,j)
            rw4(i,j)=VRZM(i,j)
         end do
      end do
c
      IF (.NOT. ALLOCATED(RVF)) THEN
         ALLOCATE( RVF(nveq1,MEDIM), ZVF(nveq1,MEDIM))
         ALLOCATE(RVFM(nveq1,MEDIM),ZVFM(nveq1,MEDIM))
      END IF
C
      NVSTRT    =  1
      call FFTDRIVER( RW1,   RVF,    FORWD,  NVEQ1, NVEQ1
     &                      ,NVSTRT, MEDIM,  NCHI
     &                      ,KUOUT,  IERSUB, IERPLC,  IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RVF    in ',IERSUB
         call ABORTRUN
     &           (SUBNAM,      9,   MESSAGE
     &           ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
      call FFTDRIVER( RW2,   ZVF,    FORWD,  NVEQ1, NVEQ1
     &                      ,NVSTRT, MEDIM,  NCHI
     &                      ,KUOUT,  IERSUB, IERPLC,  IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: ZVF    in ',IERSUB
         call ABORTRUN
     &           (SUBNAM,     10,   MESSAGE
     &           ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
      call FFTDRIVER( RW3,   RVFM,   FORWD,  NVEQ1, NVEQ1
     &                      ,NVSTRT, MEDIM,  NCHI
     &                      ,KUOUT,  IERSUB, IERPLC,  IERR)
      call FFTDRIVER( RW4,   ZVFM,   FORWD,  NVEQ1, NVEQ1
     &                      ,NVSTRT, MEDIM,  NCHI
     &                      ,KUOUT,  IERSUB, IERPLC,  IERR)
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,  RVF,    NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'RVF')
        call FFTOUTPT(RW2,  ZVF,    NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'ZVF')
        call FFTOUTPT(RW3,  RVFM,   NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'RVFM')
        call FFTOUTPT(RW4,  ZVFM,   NVEQ1,  NVEQ1,   NVSTRT
     &                     ,MEDIM,  NCHI,   KUFFTP, 'ZVFM')
      ENDIF
C
C
      IF (.NOT. ALLOCATED(ZCNDF)) THEN
         ALLOCATE(ZCNDF(NALWALL,MEDIM))
      END IF

      IF (NWALL.GT.0) THEN
        NWSTRT    =  1
        call FFTDRIVER( ZCND,  ZCNDF,  FORWD,  NALWALL, NWALL
     &                        ,NWSTRT, MEDIM,  NCHI
     &                        ,KUOUT,  IERSUB, IERPLC,  IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: ZCNDF  in ',IERSUB
           call ABORTRUN
     &           (SUBNAM,      9,   MESSAGE
     &           ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        IF(KUFFTP .GT. 0) THEN
          call FFTOUTPT(ZCND, ZCNDF,  NALWALL,NWALL,   NWSTRT
     &                       ,MEDIM,  NCHI,   KUFFTP, 'ZCNDF')
        ENDIF
      ENDIF

      if(1.eq.0) then
      CLOSE(FOURIER_OUT)
      OPEN(FOURIER_OUT,FILE='fourier_vac.OUT',FORM='FORMATTED')
      REWIND(FOURIER_OUT)

      WRITE(FOURIER_OUT,1000) NVEQ1,MEDIM,NSMAX


      CALL GENOUT(VDG11L,nveq1,' VDG11L',FOURIER_OUT,NVEQ1,MEDIM,1)
      CALL GENOUT(VDG22L,nveq1,' VDG22L',FOURIER_OUT,NVEQ1,MEDIM,1)
      CALL GENOUT(VDG33L,nveq1,' VDG33L',FOURIER_OUT,NVEQ1,MEDIM,1)
      CALL GENOUT(VDG12L,nveq1,' VDG12L',FOURIER_OUT,NVEQ1,MEDIM,1)

      CALL GENOUT(VDG11LM,nveq1,'VDG11LM',FOURIER_OUT,NVEQ,MEDIM,1)
      CALL GENOUT(VDG22LM,nveq1,'VDG22LM',FOURIER_OUT,NVEQ,MEDIM,1)
      CALL GENOUT(VDG33LM,nveq1,'VDG33LM',FOURIER_OUT,NVEQ,MEDIM,1)
      CALL GENOUT(VDG12LM,nveq1,'VDG12LM',FOURIER_OUT,NVEQ,MEDIM,1)
      CLOSE(FOURIER_OUT)
      endif
C
      DEALLOCATE(RW1,RW2,RW3,RW4,RW6,RW7,RW8,RW9)
C
      RETURN
 1000 FORMAT(3I20)
      
      End   !!!---{Out_Fourier_Vac}
      
*DECK READTOR
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------INPUT ROUTINE-------------A.B.  16/07/1997-------------------
C COMPUTE Q-PROFILE RIGHT AFTER READING EQUILIBRIUM DATA
C STORE Q-PROFILE IN <QPLS>
C YQ LIU, MARCH 20, 2012
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE READTOR
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE GIJLM
      USE MPIENV
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'

      INTEGER I,J,K,L,NRP10,NR0
      REAL*8 X, TTP, TTPM, TMP
      REAL*8,DIMENSION(:),ALLOCATABLE::PSIISO
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RTMP  

      OPEN(NEQBIN,FILE='OUTRMAR',FORM='FORMATTED')

      READ(NEQBIN,*) NRP1, NCHI

C     CALCULATE MEDIM BASED ON TWO CONDITIONS
      MEDIM = MAX(MSMAX+10,INT(NCHI/4))
      IF (MEDIM.GT.INT(NCHI/4)) STOP 'INCREASE NCHI'
      IF (MOUTPUT.GT.MEDIM) MOUTPUT = MEDIM

      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*) 'NRP1, NCHI=, MEDIM=',NRP1,NCHI,MEDIM
      ENDIF

      IF (NRP1 .GT. NCRAY) STOP 'NCRAY'
      IF (NCHI .GT. NRSDIM) STOP 'NRSDIM'
      NR = NRP1 - 1
      READ(NEQBIN,*) X,R0EXP,B0EXP
      ASPCT = X

C     STORAGE FOR 1D and 2D EQUILIBRIUM QUANTITIES
C
      IF (.NOT. ALLOCATED(CSE)) THEN
         ALLOCATE(CSE(NRP1),CSEM(NRP1),T(NRP1),
     &            TM(NRP1),PPEQ(NRP1),PPEQM(NRP1),TP(NRP1),TPM(NRP1),
     &            DPSIDS(NRP1),DPSIDSM(NRP1),
     &            PEQ(NRP1),PEQM(NRP1),
     &            ROTP(NRP1), ROTPM(NRP1),
     &            RHOU(NRP1), RHOUM(NRP1),
     &            DRHOU(NRP1),DRHOUM(NRP1))

C        ELEMENT NRP1 IS READ (TIMES A ZERO EDGE WEIGHT) BUT NEVER SET:
C        START FROM ZERO, HEAP GARBAGE NaN WOULD GIVE 0*NaN = NaN
         DPSIDSM = 0.
         TM      = 0.
         PPEQM   = 0.
         TPM     = 0.
      END IF

      ALLOCATE(PSIISO(2*NRP1),RTMP(NRP1,NCHI))
C

      IF (.NOT. ALLOCATED(REQ)) THEN
         ALLOCATE(REQ(NRP1,NCHI),  REQM(NRP1,NCHI),
     &        ZEQ(NRP1,NCHI),      ZEQM(NRP1,NCHI),
     &        RJA(NRP1,NCHI),      RJAM(NRP1,NCHI),
     &        G11L(NRP1,NCHI),    G11LM(NRP1,NCHI),
     &        G12L(NRP1,NCHI),    G12LM(NRP1,NCHI),
     &        G22L(NRP1,NCHI),    G22LM(NRP1,NCHI),
     &        G33L(NRP1,NCHI),    G33LM(NRP1,NCHI),
     &        RDCDZ(NRP1,NCHI),  RDCDZM(NRP1,NCHI),
     &        RDSDZ(NRP1,NCHI),  RDSDZM(NRP1,NCHI),
     &        RBZ(NRP1,NCHI),    RBZM(NRP1,NCHI)  ,
     &        ROTC(NRP1,NCHI),   ROTCM(NRP1,NCHI) )
      END IF
C
      DO 10 J = 1,NR
      I = J
      READ(NEQBIN,*) CSEM(I),PSIISO(2*I),PEQM(I),TM(I)
      READ(NEQBIN,*) TTPM,PPEQM(I),DPSIDSM(I)
      READ(NEQBIN,*) (REQM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (ZEQM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RJAM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (G11LM(I,L),l=1,NCHI)
      READ(NEQBIN,*) (G22LM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (G33LM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (G12LM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RDCDZM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RDSDZM(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RBZM(I,L),L=1,NCHI)
      TPM(I) = TTPM/TM(I)
C
      I = J + 1 
      READ(NEQBIN,*) CSE(I),PSIISO(2*I-1),PEQ(I),T(I)
      READ(NEQBIN,*) TTP,PPEQ(I),DPSIDS(I)
      READ(NEQBIN,*) (REQ(I,L),L=1,NCHI)
      READ(NEQBIN,*) (ZEQ(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RJA(I,L),L=1,NCHI)
      READ(NEQBIN,*) (G11L(I,L),l=1,NCHI)
      READ(NEQBIN,*) (G22L(I,L),L=1,NCHI)
      READ(NEQBIN,*) (G33L(I,L),L=1,NCHI)
      READ(NEQBIN,*) (G12L(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RDCDZ(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RDSDZ(I,L),L=1,NCHI)
      READ(NEQBIN,*) (RBZ(I,L),L=1,NCHI)
C
      TP(I) = TTP/T(I)
C
 10   CONTINUE
      CLOSE(NEQBIN)
      
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("AFTER READTOR, I = ",I5)')I
      ENDIF
      
C     RE-COMPUTE JACOBIAN
      DO I= 1,NR+1
         DO L=1,NCHI
            RJA(I,L) = SQRT(ABS(G33L(I,L)*(G11L(I,L)*G22L(I,L)-
     &                      G12L(I,L)**2)))
         ENDDO
      ENDDO
      DO I= 1,NR
         DO L=1,NCHI
            RJAM(I,L) = SQRT(ABS(G33LM(I,L)*(G11LM(I,L)*G22LM(I,L)-
     &                       G12LM(I,L)**2)))
         ENDDO
      ENDDO

C     PERFORM OTHER ADJUSTMENTS 
      CSE(1) = 0.
      PEQ(1) = (4.*PEQM(1)-PEQ(2))/3.
      PPEQ(1) = 0.
      T(1)   = (4.*TM(1)  -T(2))  /3.
      TP (1) = 0.
      DPSIDS(1) = 0.
      TM(NRP1) = T(NRP1)
C
      X   = 0.
      TMP = 0.
      DO J=1,NCHI
         X   = X   + REQM(1,J)
         TMP = TMP + ZEQM(1,J)
      ENDDO
      
      DO 40 J = 1,NCHI
        REQ(1,J) = X/NCHI
        ZEQ(1,J) = TMP/NCHI
        RJA(1,J) = 0.
        G11L(1,J) = G11LM(1,J)
        G22L(1,J) = 0.
        G33L(1,J) = G33LM(1,J)
        G12L(1,J) = 0.
 40   CONTINUE

C     REDEFINE ALL INPUT EQUILIBRIUM QUANTITIES ON A NEW RADIAL MESH
C     STARTING FROM SLEFT (DEFAULT SLEFT=0.0)
      IF (SLEFT.GT.0.0.AND.SLEFT.LT.1.0) THEN
         L = 1
         DO I=1,NR
            IF (SLEFT.GE.CSE(I)) L = I
         ENDDO
         L = L+1
         WRITE(*,*) 'SLEFT=',SLEFT,' NLEFT=',L

         NRP10 = NRP1
         NR0   = NR
         NRP1  = NRP10 - L + 1
         NR    = NR0 - L + 1

         PSIISO(1:NRP10) = CSE
         DEALLOCATE(CSE)
         ALLOCATE(CSE(NRP1))   
         CSE = PSIISO(L:NRP10)
         PSIISO(1:NRP10) = CSEM
         DEALLOCATE(CSEM)
         ALLOCATE(CSEM(NRP1))   
         CSEM = PSIISO(L:NRP10)

         PSIISO(1:NRP10) = T
         DEALLOCATE(T)
         ALLOCATE(T(NRP1))   
         T = PSIISO(L:NRP10)
         PSIISO(1:NRP10) = TM
         DEALLOCATE(TM)
         ALLOCATE(TM(NRP1))   
         TM = PSIISO(L:NRP10)

         PSIISO(1:NRP10) = PPEQ
         DEALLOCATE(PPEQ)
         ALLOCATE(PPEQ(NRP1))   
         PPEQ = PSIISO(L:NRP10)
         PSIISO(1:NRP10) = PPEQM
         DEALLOCATE(PPEQM)
         ALLOCATE(PPEQM(NRP1))   
         PPEQM = PSIISO(L:NRP10)

         PSIISO(1:NRP10) = TP
         DEALLOCATE(TP)
         ALLOCATE(TP(NRP1))   
         TP = PSIISO(L:NRP10)
         PSIISO(1:NRP10) = TPM
         DEALLOCATE(TPM)
         ALLOCATE(TPM(NRP1))   
         TPM = PSIISO(L:NRP10)

         PSIISO(1:NRP10) = DPSIDS
         DEALLOCATE(DPSIDS)
         ALLOCATE(DPSIDS(NRP1))   
         DPSIDS = PSIISO(L:NRP10)
         PSIISO(1:NRP10) = DPSIDSM
         DEALLOCATE(DPSIDSM)
         ALLOCATE(DPSIDSM(NRP1))   
         DPSIDSM = PSIISO(L:NRP10)

         PSIISO(1:NRP10) = PEQ
         DEALLOCATE(PEQ)
         ALLOCATE(PEQ(NRP1))   
         PEQ = PSIISO(L:NRP10)
         PSIISO(1:NRP10) = PEQM
         DEALLOCATE(PEQM)
         ALLOCATE(PEQM(NRP1))   
         PEQM = PSIISO(L:NRP10)

         RTMP = REQ
         DEALLOCATE(REQ)
         ALLOCATE(REQ(NRP1,NCHI))   
         REQ = RTMP(L:NRP10,:)
         RTMP = REQM
         DEALLOCATE(REQM)
         ALLOCATE(REQM(NRP1,NCHI))   
         REQM = RTMP(L:NRP10,:)

         RTMP = ZEQ
         DEALLOCATE(ZEQ)
         ALLOCATE(ZEQ(NRP1,NCHI))   
         ZEQ = RTMP(L:NRP10,:)
         RTMP = ZEQM
         DEALLOCATE(ZEQM)
         ALLOCATE(ZEQM(NRP1,NCHI))   
         ZEQM = RTMP(L:NRP10,:)

         RTMP = RJA
         DEALLOCATE(RJA)
         ALLOCATE(RJA(NRP1,NCHI))   
         RJA = RTMP(L:NRP10,:)
         RTMP = RJAM
         DEALLOCATE(RJAM)
         ALLOCATE(RJAM(NRP1,NCHI))   
         RJAM = RTMP(L:NRP10,:)

         RTMP = G11L
         DEALLOCATE(G11L)
         ALLOCATE(G11L(NRP1,NCHI))   
         G11L = RTMP(L:NRP10,:)
         RTMP = G11LM
         DEALLOCATE(G11LM)
         ALLOCATE(G11LM(NRP1,NCHI))   
         G11LM = RTMP(L:NRP10,:)

         RTMP = G12L
         DEALLOCATE(G12L)
         ALLOCATE(G12L(NRP1,NCHI))   
         G12L = RTMP(L:NRP10,:)
         RTMP = G12LM
         DEALLOCATE(G12LM)
         ALLOCATE(G12LM(NRP1,NCHI))   
         G12LM = RTMP(L:NRP10,:)

         RTMP = G22L
         DEALLOCATE(G22L)
         ALLOCATE(G22L(NRP1,NCHI))   
         G22L = RTMP(L:NRP10,:)
         RTMP = G22LM
         DEALLOCATE(G22LM)
         ALLOCATE(G22LM(NRP1,NCHI))   
         G22LM = RTMP(L:NRP10,:)

         RTMP = G33L
         DEALLOCATE(G33L)
         ALLOCATE(G33L(NRP1,NCHI))   
         G33L = RTMP(L:NRP10,:)
         RTMP = G33LM
         DEALLOCATE(G33LM)
         ALLOCATE(G33LM(NRP1,NCHI))   
         G33LM = RTMP(L:NRP10,:)

         RTMP = RDCDZ
         DEALLOCATE(RDCDZ)
         ALLOCATE(RDCDZ(NRP1,NCHI))   
         RDCDZ = RTMP(L:NRP10,:)
         RTMP = RDCDZM
         DEALLOCATE(RDCDZM)
         ALLOCATE(RDCDZM(NRP1,NCHI))   
         RDCDZM = RTMP(L:NRP10,:)

         RTMP = RDSDZ
         DEALLOCATE(RDSDZ)
         ALLOCATE(RDSDZ(NRP1,NCHI))   
         RDSDZ = RTMP(L:NRP10,:)
         RTMP = RDSDZM
         DEALLOCATE(RDSDZM)
         ALLOCATE(RDSDZM(NRP1,NCHI))   
         RDSDZM = RTMP(L:NRP10,:)

         RTMP = RBZ
         DEALLOCATE(RBZ)
         ALLOCATE(RBZ(NRP1,NCHI))   
         RBZ = RTMP(L:NRP10,:)
         RTMP = RBZM
         DEALLOCATE(RBZM)
         ALLOCATE(RBZM(NRP1,NCHI))   
         RBZM = RTMP(L:NRP10,:)
      ENDIF

      IF (.NOT. ALLOCATED(RESIST)) THEN
         ALLOCATE(RESIST(NRP1),RESISM(NRP1),
     &            RHO(NRP1),RHOM(NRP1),ROT(NRP1),ROTM(NRP1),
     &            DROT(NRP1),DROTM(NRP1),GMUNU(NRP1),GMUNUM(NRP1),
     &            GAMARR(NRP1),GAMARM(NRP1),TEMPE(NRP1),TEMPEM(NRP1),
     &            TEMPI(NRP1),TEMPIM(NRP1),GNUI(NRP1),GNUIM(NRP1),
     &            GNUE(NRP1),GNUEM(NRP1),JEQ(NRP1),
     &            GNEOFUNC(NRP1),GNEOFUNCM(NRP1),
     &            DAK(NRP1),DAKM(NRP1),
     &            ZEFFI(NRP1),ZEFFM(NRP1),
     &            DLNRHO(NRP1),DLNRHOM(NRP1),
     &            OMEGASI(NRP1),OMEGASIM(NRP1),
     &            OMEGASE(NRP1),OMEGASEM(NRP1),
     &            DOMEGASI(NRP1),DOMEGASIM(NRP1),
     &            DOMEGASE(NRP1),DOMEGASEM(NRP1),
     &            TCHIMI(NRP1),TCHIMM(NRP1),
     &            TVPINCHI(NRP1),TVPINCHM(NRP1),
     &            TCHIDI(NRP1),TCHIDM(NRP1),
     &            TTCPARAI(NRP1),TTCPARAM(NRP1),
     &            TTCPERPI(NRP1),TTCPERPM(NRP1),
     &            TROTI(NRP1),TROTM(NRP1),
     &            TDROTI(NRP1),TDROTM(NRP1),
     &            TRHOI(NRP1),TRHOM(NRP1),
     &            ROTWEI(NRP1),ROTWEM(NRP1),
     &            ROTEQ(NRP1), ROTEQM(NRP1),
     &            DROTEQ(NRP1),DROTEQM(NRP1),
     &            RHOEQ(NRP1), RHOEQM(NRP1) )
         OMEGASE = 0.
         TROTI  = 0.
         TROTM  = 0.
         TDROTI = 0.
         TDROTM = 0.
         TRHOI  = 0.
         TRHOM  = 0.
      END IF

C     COMPUTE Q-PROFILE
      CALL GETQPLS

C     GENERATE NEW CS (S,THETA,PHI) WITH GEOMETRIC POLOIDAL ANGLE THETA
      IF (NCONVCS.EQ.1) CALL GEOMCS_PLASMA


      IF (.NOT.ALLOCATED(FFF)) THEN
         ALLOCATE( FFF(NRP1,2),DFFF(NRP1,2) )
      ENDIF

      DEALLOCATE(PSIISO,RTMP)

      RETURN
      END
C
      SUBROUTINE INPUT(PS,PSIISO,CPR,TMF,TTP,CPPR,DPSIDS,
     &                 R,Z,JAC,GSS,GTT,GPP,GST,RTC,RTS,RTZ)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'comioc.inc'
C
      INTEGER L
      REAL*8 PS,PSIISO,CPR,TMF,TTP,CPPR,DPSIDS
      REAL*8 R(NRP1,NCHI),Z(NRP1,NCHI),JAC(NRP1,NCHI),
     &     GSS(NRP1,NCHI),GTT(NRP1,NCHI),
     $     GPP(NRP1,NCHI),GST(NRP1,NCHI),
     $     RTC(NRP1,NCHI),RTS(NRP1,NCHI),
     $     RTZ(NRP1,NCHI)
C
      READ(NEQBIN,*) PS,PSIISO,CPR,TMF
      READ(NEQBIN,*) TTP,CPPR,DPSIDS
      READ(NEQBIN,*) (R(1,L),L=1,NCHI)
      READ(NEQBIN,*) (Z(1,L),L=1,NCHI)
      READ(NEQBIN,*) (JAC(1,L),L=1,NCHI)
      READ(NEQBIN,*) (GSS(1,L),l=1,NCHI)
      READ(NEQBIN,*) (GTT(1,L),L=1,NCHI)
      READ(NEQBIN,*) (GPP(1,L),L=1,NCHI)
      READ(NEQBIN,*) (GST(1,L),L=1,NCHI)
      READ(NEQBIN,*) (RTC(1,L),L=1,NCHI)
      READ(NEQBIN,*) (RTS(1,L),L=1,NCHI)
      READ(NEQBIN,*) (RTZ(1,L),L=1,NCHI)
C
      RETURN
      END
c     =================================================================
c     =================================================================

      SUBROUTINE DERCHI(ARR,ARRDCH,NACT,NDIM)
      USE DIMENSIM
      USE GLOBALM
      USE RESMATM
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER NACT,NDIM
      INTEGER MS,I
      REAL*8  ARR(NDIM,NCHI),ARRDCH(NDIM,NCHI)
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::FT
      INTEGER J,J1,J2
      REAL*8 DCHII,DIFF,A,D
C
C     IF (NDIM.NE.NCRAY) STOP 'NDIM_DERCHI'
C
      ALLOCATE(FT(NDIM,MEDIM))
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'DERCHI'
      IERF      = 0
      IERI      = 0
C
      if (.true.) goto 21
        NASTRT    =  1
        call FFTDRIVER(ARR,    FT,     FORWD, NDIM,  NACT,   NASTRT
     &                        ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC
     &                                                     , IERF)
        if(IERF .NE. 0) THEN
           write(MESSAGE,*) 'Error: FT  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERF    ', IERF,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        DO 20 MS = 1,MEDIM
        DO 10 I = 1,NACT
        FT(I,MS) = CMPLX(0.,MS-1.) *FT(I,MS)
 10     CONTINUE
 20     CONTINUE
C
      call FFTDRIVER(ARRDCH, FT,     BCKWD, NDIM,  NACT,  NASTRT
     &                      ,MEDIM,  NCHI, KUOUT, IERSUB, IERPLC
     &                                                 ,  IERI)
      if(IERI .NE. 0) THEN
         write(MESSAGE,*) 'Error: ARRDCH  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERI    ', IERI,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(ARR,   FT,     NDIM,   NACT,    NASTRT
     &                      ,MEDIM,  NCHI,   KUFFTP, 'FT')
        call FFTOUTPT(ARRDCH,FT,     NDIM,   NACT,    NASTRT
     &                      ,MEDIM,  NCHI,   KUFFTP, 'ARRDCH')
      ENDIF
 21   CONTINUE
C
      DCHII= DFLOAT(NCHI)/(4.*PI)
      DIFF = 0.
      DO 40 J = 1,NCHI
      J1 = J-1
      IF (J1.EQ.0) J1 = NCHI
      J2 = J+1
      IF (J.EQ.NCHI) J2 = 1
      DO 30 I = 1,NACT
      A = (ARR(I,J2)-ARR(I,J1))*DCHII
C     D = ABS(A-ARRDCH(I,J))
C     IF (D.LT.DIFF) GOTO 25
C     WRITE(*,*) 'i =',i,' j=',j,' d=',d
C     DIFF = D
 25   ARRDCH(I,J) = A
 30   CONTINUE
 40   CONTINUE
C     WRITE(*,*) ' DERCHI MAXIMUM ERROR IS',DIFF
      DEALLOCATE(FT)
      RETURN
      END 

*DECK READVAC
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------INPUT ROUTINE-------------A.B.  25/10/1997-------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE READVAC
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE GIJLM
      USE FEEDBACKM
      USE MPIENV
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'
      INTEGER I,J,J1,NBPS
C
      
      OPEN(NVACBIN,FILE='OUTVMAR',FORM='FORMATTED')
      REWIND(NVACBIN)
C
      READ(NVACBIN,*) NVEQ1, NCHI, NWBPS
      
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      write(*,*) 'NVEQ1, NCHI, NWBPS=',NVEQ1,NCHI,NWBPS
      ENDIF
      
      IF (NCHI .GT. NRSDIM) STOP 'NRSDIM'
      NVEQ = NVEQ1 - 1
C
      IF (.NOT. ALLOCATED(VRJA)) THEN
         ALLOCATE(  VRJA(nveq1,nchi),   VRJAM(nveq1,nchi),
     &            VRG11L(nveq1,nchi), VRG11LM(nveq1,nchi),
     &            VRG12L(nveq1,nchi), VRG12LM(nveq1,nchi),
     &            VRG22L(nveq1,nchi), VRG22LM(nveq1,nchi),
     &            VRG33L(nveq1,nchi), VRG33LM(nveq1,nchi),
     &               VRR(nveq1,nchi),    VRRM(nveq1,nchi),
     &               VRZ(nveq1,nchi),    VRZM(nveq1,nchi))
      END IF

      NTP1=NV+NRP1
      NTOT=NTP1-1

      IF (.NOT. ALLOCATED(VCS)) THEN
         ALLOCATE(     VCS(NVEQ1),VCSM(NVEQ1),VCSH(0:NVEQ1+1))
         VCS=0.
         VCSM=0.
         VCSH=0.
         ALLOCATE(CS(NTP1+4),CSM(NTP1),CSV(NTP1),CSVM(NTP1),
     &                       CSH(0:NTP1+1),Q(NTP1),QM(NTP1))
      END IF
C
      CS(1:NRP1)=CSE(1:NRP1)
      CSM(1:NR)=CSEM(1:NR)
C
      DO 10 J = 1,NVEQ1
      I = J
      CALL INPUTV(VCS(I),VRJA(I,1),
     &     VRG11L(I,1),VRG22L(I,1),VRG33L(I,1),VRG12L(I,1),
     &     VRR(I,1),VRZ(I,1))
      IF (J .EQ. NVEQ1) GOTO 10
      CALL INPUTV(VCSM(I),VRJAM(I,1),
     &     VRG11LM(I,1),VRG22LM(I,1),VRG33LM(I,1),VRG12LM(I,1),
     &     VRRM(I,1),VRZM(I,1))
C
 10   CONTINUE
C
      CS(NRP1:NTP1) = VCS(1:NTP1-NRP1+1)
      CSM(NR+1:NTP1-1) = VCSM(1:NTP1-NR-1)
      CSH(NRP1+1:NTP1) = VCS(2:NTP1-NRP1+1)-VCS(1:NTP1-NRP1)
      CSH(NTP1+1)      = 0.0

      nalwall=max0(nwbps,nwall)

      IF (.NOT. ALLOCATED(ZCND)) THEN
         ALLOCATE(  ZCND(nalwall,nchi),   ZCNDC(nalwall,nchi))
      END IF

      IF (NALWALL.GT.0) THEN
         ZCND=1.
         ZCNDC=0.
      END IF
      DO J1=2,NWBPS-1
         READ(NVACBIN,*) (ZCND(J1-1,J),J=1,NCHI) 
         READ(NVACBIN,*) (ZCNDC(J1-1,J),J=1,NCHI) 
      ENDDO

      CLOSE(NVACBIN)
      
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      write(*,*) ' READVAC FOURIER_VACUUM'
      ENDIF
      
C     GENERATE NEW CS (S,THETA,PHI) WITH GEOMETRIC POLOIDAL ANGLE THETA
      IF (NCONVCS.EQ.1) CALL GEOMCS_VACUUM
      CALL FOURIER_VACUUM

      RETURN
      END
C
      SUBROUTINE INPUTV(PS,
     &                 JAC,GSS,GTT,GPP,GST,GR,GZ)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'comioc.inc'
C
      INTEGER L
      REAL*8 PS,JAC(NVEQ1,*),
     &     GSS(NVEQ1,*),GTT(NVEQ1,*),
     &     GPP(NVEQ1,*),GST(NVEQ1,*),
     &     GR(NVEQ1,*),GZ(NVEQ1,*)
C
      READ(NVACBIN,*) PS
      READ(NVACBIN,*) (JAC(1,L),L=1,NCHI)
      READ(NVACBIN,*) (GSS(1,L),L=1,NCHI)
      READ(NVACBIN,*) (GTT(1,L),L=1,NCHI)
      READ(NVACBIN,*) (GPP(1,L),L=1,NCHI)
      READ(NVACBIN,*) (GST(1,L),L=1,NCHI)
      READ(NVACBIN,*) (GR(1,L),L=1,NCHI)
      READ(NVACBIN,*) (GZ(1,L),L=1,NCHI)
C
      DO L=1,NCHI
         JAC(1,L) = SQRT(ABS(GSS(1,L)*GTT(1,L)-GST(1,L)**2))*GR(1,L)
      ENDDO

      RETURN
      END

*DECK UPDATE_WALL
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-UPDATE VACUUM MESH DUE TO RELATIVE POSITION SHIFT BETWEEN PLASMA AND
C WALL, AS A RESULT OF VERTICAL MOVEMENT OF PLASMA (RNTOR=0)
C IN ACTION ONLY WITH NCASE=10, RNTOR=0, NWALL>0
C LIU YQ, JUNE 20, 2019
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE UPDATE_WALL(ISW)
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      IMPLICIT NONE

      INTEGER ISW,I,J,K,IWALL0,IWALL1
      REAL*8  VCS0,VCS1,VDRDS,VDRDC,VDZDS,VDZDC,
     &        HCHI,VRR1,VRR2,VRZ1,VRZ2
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RWALN,ZWALN
C
      IF (ISW.EQ.1) THEN
         IF (.NOT.ALLOCATED(RWALL)) THEN
            ALLOCATE(RWALL(NWALL,NCHI), ZWALL(NWALL,NCHI))
         ENDIF
         DO I=1,NWALL
         DO J=1,NCHI
            RWALL(I,J) = VRR(IWALL(I),J)
            ZWALL(I,J) = VRZ(IWALL(I),J)
         ENDDO
         ENDDO
      ELSE
         ALLOCATE(RWALN(0:NWALL,NCHI), ZWALN(0:NWALL,NCHI))
         RWALN(1:NWALL,:) = RWALL - DREAL(VDEXIR)/R0EXP
         ZWALN(1:NWALL,:) = ZWALL - DREAL(VDEXIZ)/R0EXP

         RWALN(0,:) = VRR(1,:)
         ZWALN(0,:) = VRZ(1,:)

C        WORK OUT NEW (VRR,VRZ) FOR VACUUM MESH
C        WITH LINEAR INTERPOLATION BETWEEN WALLS
C        AND LINEAR EXTRAPOLATION OUTSIDE LAST WALL
         DO I=0,NWALL-1
            IF (I.EQ.0) THEN
               IWALL0 = 2
               VCS0   = 1.0
            ELSE
               IWALL0 = IWALL(I)+1
               VCS0 = VCS(IWALL0-1)
            ENDIF
            VCS1 = VCS(IWALL(I+1))
            IF (I.EQ.NWALL-1) THEN
               IWALL1 = NVEQ1
            ELSE
               IWALL1 = IWALL(I+1)
            ENDIF

            DO K=IWALL0,IWALL1
               VRR(K,:) = RWALN(I,:)*(VCS1-VCS(K))/(VCS1-VCS0) + 
     &                    RWALN(I+1,:)*(VCS(K)-VCS0)/(VCS1-VCS0)
               VRZ(K,:) = ZWALN(I,:)*(VCS1-VCS(K))/(VCS1-VCS0) + 
     &                    ZWALN(I+1,:)*(VCS(K)-VCS0)/(VCS1-VCS0)
            ENDDO
         ENDDO

C        COMPUTE NEW METRICS ELEMENTS IN VACUUM
         VRRM(1:NVEQ,:) = (VRR(1:NVEQ,:)+VRR(2:NVEQ1,:))*0.5
         VRZM(1:NVEQ,:) = (VRZ(1:NVEQ,:)+VRZ(2:NVEQ1,:))*0.5

         HCHI = 4.*PI/NCHI

C        HALF-INTEGER MESH
         DO I=1,NVEQ
         DO J=1,NCHI
            IF (J.EQ.1) THEN
               VRR1 = VRRM(I,NCHI)
               VRR2 = VRRM(I,2)
               VRZ1 = VRZM(I,NCHI)
               VRZ2 = VRZM(I,2)
            ELSEIF (J.EQ.NCHI) THEN
               VRR1 = VRRM(I,NCHI-1)
               VRR2 = VRRM(I,1)
               VRZ1 = VRZM(I,NCHI-1)
               VRZ2 = VRZM(I,1)
            ELSE
               VRR1 = VRRM(I,J-1)
               VRR2 = VRRM(I,J+1)
               VRZ1 = VRZM(I,J-1)
               VRZ2 = VRZM(I,J+1)
            ENDIF
            VDRDS = (VRR(I+1,J)-VRR(I,J))/(VCS(I+1)-VCS(I))
            VDZDS = (VRZ(I+1,J)-VRZ(I,J))/(VCS(I+1)-VCS(I))
            VDRDC = (VRR2-VRR1)/HCHI
            VDZDC = (VRZ2-VRZ1)/HCHI
         
            VRG11LM(I,J) = VDRDS**2+VDZDS**2
            VRG12LM(I,J) = VDRDS*VDRDC+VDZDS*VDZDC
            VRG22LM(I,J) = VDRDC**2+VDZDC**2
            VRG33LM(I,J) = VRRM(I,J)**2
            VRJAM(I,J)   = (VDRDS*VDZDC-VDRDC*VDZDS)*VRRM(I,J)
         ENDDO
         ENDDO

C        INTEGER MESH 
         DO I=2,NVEQ
         DO J=1,NCHI
            IF (J.EQ.1) THEN
               VRR1 = VRR(I,NCHI)
               VRR2 = VRR(I,2)
               VRZ1 = VRZ(I,NCHI)
               VRZ2 = VRZ(I,2)
            ELSEIF (J.EQ.NCHI) THEN
               VRR1 = VRR(I,NCHI-1)
               VRR2 = VRR(I,1)
               VRZ1 = VRZ(I,NCHI-1)
               VRZ2 = VRZ(I,1)
            ELSE
               VRR1 = VRR(I,J-1)
               VRR2 = VRR(I,J+1)
               VRZ1 = VRZ(I,J-1)
               VRZ2 = VRZ(I,J+1)
            ENDIF
            VDRDS = (VRR(I+1,J)-VRR(I-1,J))/(VCS(I+1)-VCS(I-1))
            VDZDS = (VRZ(I+1,J)-VRZ(I-1,J))/(VCS(I+1)-VCS(I-1))
            VDRDC = (VRR2-VRR1)/HCHI
            VDZDC = (VRZ2-VRZ1)/HCHI
         
            VRG11L(I,J) = VDRDS**2+VDZDS**2
            VRG12L(I,J) = VDRDS*VDRDC+VDZDS*VDZDC
            VRG22L(I,J) = VDRDC**2+VDZDC**2
            VRG33L(I,J) = VRR(I,J)**2
            VRJA(I,J)   = (VDRDS*VDZDC-VDRDC*VDZDS)*VRR(I,J)
         ENDDO
         ENDDO
         VRG11L(1,:) = VRG11LM(1,:)
         VRG12L(1,:) = VRG12LM(1,:)
         VRG22L(1,:) = VRG22LM(1,:)
         VRG33L(1,:) = VRG33LM(1,:)
         VRJA(1,:)   = VRJAM(1,:)
         VRG11L(NVEQ1,:) = VRG11LM(NVEQ,:)
         VRG12L(NVEQ1,:) = VRG12LM(NVEQ,:)
         VRG22L(NVEQ1,:) = VRG22LM(NVEQ,:)
         VRG33L(NVEQ1,:) = VRG33LM(NVEQ,:)
         VRJA(NVEQ1,:)   = VRJAM(NVEQ,:)
      ENDIF
         

      IF (ALLOCATED(RWALN)) DEALLOCATE(RWALN,ZWALN)

C     GENERATE NEW CS (S,THETA,PHI) WITH GEOMETRIC POLOIDAL ANGLE THETA
      IF (NCONVCS.EQ.1) CALL GEOMCS_VACUUM
      CALL FOURIER_VACUUM

      RETURN
      END
C
*DECK COTROL
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C----- CHECK VALIDITY OF NAMELIST INPUT -------------------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE COTROL
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE RCOMDM
      USE MPIENV
      USE REORBITM
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      LOGICAL  OVERLAP
      INTEGER  I, J, K, NWALLNEW, IWALLNEW(NWALL0), IPASS, ITRAP
      REAL*8   TAUWNEW(NWALL0)
      REAL*8,PARAMETER::CEPS=1.e-10

C     CONTROL NCASE
      IF (NCASE.LT.-3.OR.NCASE.GT.11) STOP 'NCASE' !here

      IF (NCASE.GE.-3.AND.NCASE.LE.-1) THEN
         CALPHA2 = 1.     
         CALPHA3 = 1.      
         CALPHA4 = -1./CALPHA1  
         CALPHA7 = 0.      
         INORMSOL= 0
      ENDIF

      IF (NCASE.GE.0.AND.NCASE.LE.1) THEN
         CALPHA2 = 1.      
         CALPHA3 = 1.     
         CALPHA4 = -1./CALPHA1  
         CALPHA7 = 0.      
         IF (NITMAX.LE.1) STOP 'NITMAX'
         INORMSOL= 1
      ENDIF

      IF (NCASE.EQ.2) THEN
         CALPHA2 = 1.      
         CALPHA3 = 1.     
         CALPHA4 = (0.,0.)     
         CALPHA7 = 0.     
         NSWEEP  = 1
         NITMAX  = 1
         INORMSOL= 0
      ENDIF

      IF (NCASE.EQ.3) THEN
         IF (CALPHA2.LE.0.5) STOP 'TALPHA2'
         CALPHA4 = (1.,0.)      
         CALPHA7 = 0.    
         NSWEEP  = 1
         IF (NITMAX.LE.1) STOP 'NITMAX'
         INORMSOL= 0
      ENDIF

      IF (NCASE.EQ.4) THEN
         CALPHA2 = 1.     
         CALPHA3 = 1.     
         CALPHA4 =-1./CALPHA1   
         CALPHA7 = 0.     
         NSWEEP  = 3
         IF (NITMAX.LE.1) STOP 'NITMAX'
         INORMSOL= 1
      ENDIF

      IF (NCASE.EQ.5.OR.NCASE.EQ.9) THEN
         IF (CALPHA2.LE.0.5) STOP 'TALPHA2'
         CALPHA3 = 0.     
         CALPHA4 = (1.,0.)      
         CALPHA5 = CALPHA1
         IF (CALPHA6.LE.0.5) STOP 'TALPHA6'
         CALPHA7 = 1.     
         NSWEEP  = 1
         IF (NITMAX.LT.1) STOP 'NITMAX'
         INORMSOL= 0
      ENDIF

      IF (NCASE.EQ.6) THEN
         CALPHA1 = REAL(TALPHA1)
         IF (CALPHA2.LE.0.5) STOP 'TALPHA2'
         CALPHA4 = (1.,0.)      
         CALPHA5 = CALPHA1
         IF (CALPHA6.LE.0.4) STOP 'TALPHA6'
         CALPHA7 = 1.       
         IF (CALPHA8.LE.0..OR.CALPHA8.GT.1.) STOP 'TALPHA8'
         IF (NSWEEP.LE.1) STOP 'NSWEEP'
         NITMAX  = 1
         INORMSOL= 0
      ENDIF

      IF (NCASE.EQ.7) THEN
         CALPHA2 = 1.      
         CALPHA3 = 1.      
         CALPHA4 = (0.,0.)      
         IF (CALPHA6.LE.0.5) STOP 'TALPHA6'
         CALPHA7 = 1.       
         IF (CALPHA8.LE.0..OR.CALPHA8.GT.1.) STOP 'TALPHA8'
         IF (NSWEEP.LE.1) STOP 'NSWEEP'
         NITMAX  = 1
         INORMSOL= 0
      ENDIF

      IF (NCASE.EQ.8) THEN
         IF (CALPHA6.LE.0.5) STOP 'TALPHA6'
         CALPHA7 = 1.       
         IF (NITMAX.LE.1) STOP 'NITMAX'
         NSWEEP  = 1
         INORMSOL= 0
      ENDIF

      IF (NCASE.EQ.10) THEN
         CALPHA1 = REAL(TALPHA1)
         IF (CALPHA2.LE.0.5) STOP 'TALPHA2'
         IF (CALPHA3.LT.0.0.OR.CALPHA3.GT.1.0) STOP 'TALPHA3'
         CALPHA4 = (1.,0.)      
         IF (CALPHA8.LE.0..OR.CALPHA8.GT.1.) STOP 'TALPHA8'
         IF (NSWEEP.LE.1) STOP 'NSWEEP'
         NITMAX  = 1
         INORMSOL= 0
      ENDIF

      IF (NSWEEP.LT.1) NSWEEP = 1
      IF (NV.LT.2) NV = 0
      IF (NV.LE.NVEQ1-1.OR.NV.LE.0) GOTO 10
      NV   = NVEQ1 - 1
      WRITE(*,*) ' *** WARNING *** NV SET TO NVEQ=',NV
 10   CONTINUE
      NVP1 = NV + 1
      
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      IF (NV.GT.0) WRITE(*,1000) VCS(NVP1),NV
      IF (NV.EQ.0) WRITE(*,*) ' NO VACUUM, WALL ON PLASMA SURFACE'
      ENDIF
      
      NTOT = NR + NV
      NTP1 = NRP1 + NV
      IF (NTP1.GT.NCRAY) STOP 'NTP1'
C
      NWALLNEW = 0
      DO J=1,NWALL
         IF (IWALL(J).GT.1.AND.IWALL(J).LT.NV) THEN
            OVERLAP = .FALSE.
            DO I=1,NWALLNEW
               IF ( (.NOT.OVERLAP).AND.(IWALL(J).EQ.IWALLNEW(I)) ) THEN
                  OVERLAP = .TRUE.
                  TAUWNEW(I) = TAUWNEW(I) + TAUW(J)
               ENDIF
            ENDDO
            IF (.NOT.OVERLAP) THEN
               NWALLNEW = NWALLNEW + 1
               IWALLNEW(NWALLNEW) = IWALL(J)
               TAUWNEW(NWALLNEW)  = TAUW(J)
            ENDIF
         ENDIF
      ENDDO
      IF (NWALLNEW.LT.NWALL) THEN
         IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
         WRITE(*,2000) NWALL,NWALLNEW
         WRITE(*,*) ' IWALLNEW  TAUWNEW'
         ENDIF
         
         DO I=1,NWALLNEW
            IWALL(I) = IWALLNEW(I)
            TAUW(I)  = TAUWNEW(I)
            
            IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
            WRITE(*,2010) IWALL(I),TAUW(I)
            ENDIF
            
         ENDDO
         NWALL = NWALLNEW
      ENDIF
 2000 FORMAT(//,' *** WARNING *** NWALL = ',I2,', NWALLNEW = ',I2)
 2010 FORMAT(4X,I2,5X,E13.4)

      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*)
      WRITE(*,*) ' MEDIM=',MEDIM,' MSMAX=',MSMAX
      ENDIF
      
 25   IF (MSMAX.LE.MSDIM) GOTO 30
      WRITE(*,*) ' M2 - M1 = ',MSMAX,' > MSDIM = ',MSDIM
      STOP 'MSMAX'
 30   CONTINUE
C
      IF (MSMAX.LE.MEDIM) GOTO 50
      WRITE(*,*) ' # PERTURBED M:',MSMAX,' # EQUILIBRIUM M:',MEDIM
      WRITE(*,*) ' SHOULD HAVE MEDIM .GE. MSMAX (= M2 - M1 + 1)'
      STOP 'MEDIM'
 50   CONTINUE

      IF (NSPECIES.GT.NSPECIES0) STOP 'NSPECIES > NSPECIES0'
      IF (NPICK.GT.NPICK0) STOP 'NPICK > NPICK0'
C
      IF (ABS(PVISC).GT.0..AND.(ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT)) 
     &WRITE(*,*) ' *** WARNING *** PVISC NONZERO, VALUE=',PVISC
 
      IF (NUII .EQ.0.) GOTO 70
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT)
     &   WRITE(*,*) ' *** WARNING *** NUII  NONZERO, VALUE=',NUII
      IF (NUII .LT.0.) STOP 'NEGATIVE'
 70   CONTINUE

      IF (INCKIN.GT.0.AND.GAMMA.GT.0.
     &    .AND.(ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT)) THEN
      WRITE(*,*) 'WARNING: SHOULD SET GAMMA=0  FOR FULL KINETIC DAMPING'
C     GAMMA = 0.0
      ENDIF

      IF (INUTYPE.NE.0) THEN
        IF (ABS(NUEFFIA).LT.1E-15 .OR. ABS(NUEFFEA).LT.1E-15) THEN
            STOP 'NUMERICAL INTEGRATION IS NOT VALID 
     &             FOR COLLISIONLESS CASE'
        ENDIF
      ENDIF

      IF (INCDPHI.GT.0) THEN
         WRITE(*,*) 'WARNING: NEED TO SET V2XKEY=1 OR 3'
         V2XKEY = 1
      ENDIF

C     MUST BE AT LEAST TWO PARTICLE SPECIES
      IF (NSPECIES.LT.2) STOP 'NSPECIES<2'

C     THE FIRST SPECIES MUST BE THERMAL IONS 
      IF (ISPECIES_F0(1).NE.0.OR.ESPECIES_Z(1).LT.0.) 
     &   STOP 'ISPECIES_F0(1)'

C     THE SECOND SPECIES MUST BE THERMAL ELECTRONS
      IF (ISPECIES_F0(2).NE.0.OR.ESPECIES_Z(2).NE.-1.) 
     &   STOP 'ISPECIES_F0(2)'

C     DERIVE IPARTICLE 
      IPASS     = 0
      ITRAP     = 0
      IPARTICLE = 0
      IF (SUM(ABS(PSPECIES_NP)).GT.0.)                         IPASS=1
      IF (SUM(ABS(PSPECIES_NTB))+SUM(ABS(PSPECIES_NTD)).GT.0.) ITRAP=1
      IF (IPASS.EQ.1.AND.ITRAP.EQ.1) IPARTICLE = 1
      IF (IPASS.EQ.1.AND.ITRAP.EQ.0) IPARTICLE = 2
      IF (IPASS.EQ.0.AND.ITRAP.EQ.1) IPARTICLE = 3

C     FOR NUMERICAL COMPUTATION OF ADIABATIC CONTRIBUTIONS 
      KANISOTROPIC  = 0
      DO K=1,NSPECIES
         IF ((ISPECIES_F0(K).EQ.3.OR.ISPECIES_F0(K).EQ.4.OR.
     &        ISPECIES_F0(K).EQ.5.OR.
     &       KFASTRUN.EQ.0).AND.(ABS(PSPECIES_AP(K)).GT.0.
     &       .OR.ABS(PSPECIES_AT(K)).GT.0.)) KANISOTROPIC = 1
      ENDDO

C     DERIVE IFOWP AND IFOWT
      IFOWP = 0
      IFOWT = 0
      IF (SUM(ABS(PSPECIES_FOWP)).GT.0) IFOWP = 1
      IF (SUM(ABS(PSPECIES_FOWT)).GT.0) IFOWT = 1

C     DETERMINE WHETHER PARTICLE HAS FINITE BIRTH ENERGY
C     ISPECIES_EK = 0: FOR IF0TYPE=0,4:   INFINITE BIRTH ENERGY
C     ISPECIES_EK = 1: FOR IF0TYPE=1,3,5,6: FINITE BIRTH ENERGY
C
      IF (.NOT.ALLOCATED(ISPECIES_EK)) 
     &   ALLOCATE( ISPECIES_EK(NSPECIES) )
      DO K=1,NSPECIES
         J = ISPECIES_F0(K)
         IF (J.EQ.0.OR.J.EQ.4)                     ISPECIES_EK(K)=0
         IF (J.EQ.1.OR.J.EQ.3.OR.J.EQ.5.OR.J.EQ.6) ISPECIES_EK(K)=1
      ENDDO

      IF (SLEFT.GE.1.0) STOP 'SLEFT.GE.1'
      IF (SLEFT.LT.0.0) STOP 'SLEFT<0'

      IF (NREORBIT.GT.0.AND.KRE_FLT.EQ.1) THEN
         RE_EFIELD = 0.0
         RE_SAC    = 0.0
         RE_SYNCH  = 0.0
         RE_BREMS  = 0.0
         RE_P0     = 1.0E-2
         RE_LAMBDA0= 0.0
         KRE_STAR  = 0
      ENDIF

      IF (RE_CONST(6).LT.0.) RE_CONST(6)=(RE_CONST(3)+1.)/2./RE_CONST(5)

      RETURN
 1000 FORMAT(/,'  WALL RADIUS =',F8.4,'   # GRID POINTS IN VACUUM =',I5)
      END
*DECK ANNIHX
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------ANNIHILATE X-EQUATION------A. BONDESON 12/05/89-----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE ANNIHX(KROW,MSROW,I,
C     ===============================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'
C
      INTEGER KROW,MSROW,I
      INTEGER MSCOL,LXROW,LXCOL,LYCOL,KX,KY
C
      LXROW = (MSROW-1)*NXCOMP
C
      DO  50 MSCOL = 1,MSMAX
C
      LXCOL = (MSCOL-1)*NXCOMP
C
      DO 30 KX = 1,NXCOMP
       ASUBM(KROW+LXROW,KX+LXCOL,I   ) = 0.
       BSUBM(KROW+LXROW,KX+LXCOL,I   ) = 0.
       CSUBM(KROW+LXROW,KX+LXCOL,I   ) = 0.
C
 30   CONTINUE
C
      LYCOL = (MSCOL-1)*NYCOMP
C
      DO 40 KY = 1,NYCOMP
C
       ESUBM(KROW+LXROW,KY+LYCOL,I   ) = 0.
       HSUBM(KROW+LXROW,KY+LYCOL,I   ) = 0.
C
 40   CONTINUE
 50   CONTINUE
C
      RETURN
      END
*DECK ANNIHY
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------ANNIHILATE Y-EQUATION------A. BONDESON 05/06/89-----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE ANNIHY(KROW,MSROW,I,
C     ===============================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'
C
      INTEGER KROW,MSROW,I
      INTEGER MSCOL,LYROW,LXCOL,LYCOL,KX,KY
C
      LYROW = (MSROW-1)*NYCOMP
C
      DO  50 MSCOL = 1,MSMAX
C
      LXCOL = (MSCOL-1)*NXCOMP
C
      DO 30 KX = 1,NXCOMP
       FSUBM(KROW+LYROW,KX+LXCOL,I   ) = 0.
       GSUBM(KROW+LYROW,KX+LXCOL,I   ) = 0.
C
 30   CONTINUE
C
      LYCOL = (MSCOL-1)*NYCOMP
C
      DO 40 KY = 1,NYCOMP
C
       DSUBM(KROW+LYROW,KY+LYCOL,I   ) = 0.
C
 40   CONTINUE
 50   CONTINUE
C
      RETURN
      END
*DECK GENOUT
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------GENERAL OUTPUT ROUTINE-------------G.VLAD 23/04/1989-----------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....GENERAL OUTPUT ROUTINE:...........................................
C.....                       VECTOR        :  VECTOR....................
C.....                       ALPHA         :  ALPHANUMERIC STRING 6 CH..
C.....                       NUNIT         :  OUTPUT UNIT NUMBER........
C.....                       NCRAY         :  RADIAL VECTOR DIMENSION...
C.....                       NR,MSMAX,NSMAX:  CURRENT VECTOR DIM........
C.....                       RM,RN         :  POLOIDAL AND TOROIDAL # ..
C
C...!!IF NSMAX .GT. 1, MSMAX MUST BE THE TRUE DIMENSION OF THE VECTORS..
C.....TO PRESERVE THE CORRECT ORDERING OF THE VECTOR ELEMENTS...........
C
      SUBROUTINE GENOUT(VECTOR,NDIM,ALPHA,NUNIT,NP,MWRITE,NS)
C     =======================================================
      USE DIMENSIM
      USE GLOBALM
C
      CHARACTER*6      ALPHA
      INTEGER          NP,MWRITE,I,MS,NS,NUNIT,NDIM,M5
      COMPLEX*16          VECTOR(NDIM,MEDIM)
C
C      M5=MIN(5,MWRITE)
        DO 100 MS=1,MWRITE
          WRITE(NUNIT,1010) ALPHA,RM(MS,NS),RN(NS)
          WRITE(NUNIT,1020) (VECTOR(I,MS),I=1,NP)
 100  CONTINUE
C
 1010 FORMAT(//,1X,A,2F20.0)
 1020 FORMAT(2D30.15)
C
      RETURN
      END
*DECK GENINP
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------GENERAL  INPUT ROUTINE-------------G.VLAD 23/04/1989-----------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....GENERAL  INPUT ROUTINE:...........................................
C.....                       VECTOR        :  VECTOR....................
C.....                       ALPHA         :  ALPHANUMERIC STRING 6 CH..
C.....                       NUNIT         :  INPUT  UNIT NUMBER........
C.....                       NCRAY         :  RADIAL VECTOR DIMENSION...
C.....                       NP,MREAD,NREAD:  VECTOR COMPONENTS TO READ.
C.....                       RM,RN         :  POLOIDAL AND TOROIDAL # ..
C
C...!!IF NSMAX .GT. 1, MSMAX MUST BE THE TRUE DIMENSION OF THE VECTORS..
C.....TO PRESERVE THE CORRECT ORDERING OF THE VECTOR ELEMENTS...........
C
      SUBROUTINE GENINP(VECTOR,NDIM,ALPHA,NUNIT,NP,MREAD,NS)
C     ======================================================
      USE DIMENSIM
      USE GLOBALM
C
      CHARACTER*6      ALPHA,ALPHA1
      INTEGER          NP,MREAD,I,MS,NS,NUNIT,NDIM
      REAL*8             XM,XN
      COMPLEX*16          VECTOR(NDIM,MREAD)
C
        DO 100 MS=1,MREAD
          READ(NUNIT,1010,END=300,ERR=150) ALPHA1,XM,XN
          IF(ALPHA1.EQ.ALPHA) GOTO 10
             WRITE(*,*) ' WRONG LABEL:',ALPHA1,' SHOULD BE ',ALPHA
             STOP
  10      IF(ABS(XM-RM(MS,NS)).GT.1.E-5.OR.
     &       ABS(XN-RN(NS)).GT.1.D-5) GOTO 200
C
          READ(NUNIT,1020,END=300,ERR=400)  (VECTOR(I,MS),I=1,NP)
C
 100  CONTINUE
C
C     WRITE(*,*) ' READ M=',XM,' N=',XN,' ARRAY=',ALPHA1
      RETURN
 150  WRITE(*,*) 'ERROR READING LABEL'
 200  WRITE(*,*) 'XM,XN,ALPHA1 ERROR IN LABEL ',XM,XN,ALPHA1
      WRITE(*,*) 'SHOULD BE ',RM(MS,NS),RN(NS),ALPHA
      STOP
C
 300  WRITE(*,*) ' END OF EQUILIBRIUM FILE'
      STOP 'GENINP'
 400  WRITE(*,*) ' ERROR READING ARRAY',XM,XN,ALPHA1
      WRITE(*,*) ' I=',I
      RETURN
 1010 FORMAT(//,1X,A,2F20.0)
 1020 FORMAT(2D30.15)
C
      END
*DECK PRIMAT
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------PRINT ROUTINE FOR DEBUGGING--------A.B.   23/05/1989-----------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE PRIMAT(STRING,C,NROW,NCOL,NDIM)
C     ==========================================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER     NROW,NCOL,NDIM,I,J
      CHARACTER*6 STRING
      COMPLEX*16     C(NDIM,NCOL)
C
      WRITE(*,1010) STRING
      DO 100 I=1,NROW
      WRITE(*,1000) (C(I,J),J=1,NCOL)
100   CONTINUE
      RETURN
1000  FORMAT(6(2(1P,E9.2),2X))
1010  FORMAT(/' PRINTOUT OF MATRIX   ',A,/)
      END
*DECK INITXY
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------INITIALIZE PAMS VECTOR----A.B.   19.07.89---------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE INITXY(X,Y)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      COMPLEX*16    X(NXCOMP,MSMAX,*),Y(NYCOMP,MSMAX,*)
      INTEGER    MS,I,K
C
      DO 100 MS=1,MSMAX
      X(KXV1 ,MS,1:NRP1) =  V1U(1:NRP1,MS)
      IF (KXX1.GT.0) X(KXX1 ,MS,1:NRP1) =  X1U(1:NRP1,MS)
      X(KXB1 ,MS,1:NTP1) =  B1U(1:NTP1,MS)
      X(KXJ2U,MS,1:NTP1) =  J2U(1:NTP1,MS)
      X(KXJ3 ,MS,1:NTP1) =  J3U(1:NTP1,MS)
      IF (KXJ2L.GT.0) X(KXJ2L,MS,1:NTP1) =  J2L(1:NTP1,MS)
      IF (KXPD.GT.0)  X(KXPD ,MS,1:NRP1) =  PDE(1:NRP1,MS)
      IF (KXPED.GT.0) X(KXPED,MS,1:NRP1) =  PED(1:NRP1,MS)
      IF (ABS(RM(MS,2)).LT.0.1) THEN
         X(KXB1,MS,1:NRP1) = B1U(1:NRP1,MS)/T(1:NRP1)
         X(KXB1,MS,NRP1+1:NTP1) = B1U(NRP1+1:NTP1,MS)/T(NRP1)
      ENDIF
      IF (KXJRE.GT.0)   X(KXJRE,MS,1:NRP1)   = JRE(1:NRP1,MS)
      IF (KXB2L.GT.0)   X(KXB2L,MS,1:NTP1)   = 0.0
      IF (KXB3L.GT.0)   X(KXB3L,MS,1:NTP1)   = 0.0
      IF (KXJRE2.GT.0)  X(KXJRE2,MS,1:NRP1)  = JRE2(1:NRP1,MS)
      IF (KXJRE3.GT.0)  X(KXJRE3,MS,1:NRP1)  = JRE3(1:NRP1,MS)
      IF (KXJRE2L.GT.0) X(KXJRE2L,MS,1:NRP1) = JRE2L(1:NRP1,MS)
      IF (KXDPHI.GT.0)  X(KXDPHI,MS,1:NRP1)  = DPHI(1:NRP1,MS)

      Y(KYV2 ,MS,1:NR) = V2U(1:NR,MS)
      IF (KYV3.GT.0)    Y(KYV3,MS,1:NR) = V3U(1:NR,MS)
      IF (KYX2.GT.0)    Y(KYX2,MS,1:NR) = X2U(1:NR,MS)
      IF (KYX3.GT.0)    Y(KYX3,MS,1:NR) = X3U(1:NR,MS)
      Y(KYB2,MS,1:NTOT) = B2U(1:NTOT,MS)
      Y(KYB3,MS,1:NTOT) = B3U(1:NTOT,MS)
      Y(KYJ1,MS,1:NTOT) = J1U(1:NTOT,MS)
      Y(KYPR ,MS,1:NR)  = PRE(1:NR,MS)
      IF (KYPE.GT.0)    Y(KYPE,MS,1:NR)    = PEE(1:NR,MS)
      IF (KYPP.GT.0)    Y(KYPP,MS,1:NR)    = PEP(1:NR,MS)
      IF (KYPPERP.GT.0) Y(KYPPERP,MS,1:NR) = PPERP(1:NR,MS)
      IF (KYPPARA.GT.0) Y(KYPPARA,MS,1:NR) = PPARA(1:NR,MS)
      IF (KYRHOP.GT.0)  Y(KYRHOP,MS,1:NR)  = RHOP(1:NR,MS)
      IF (KYJRE1.GT.0)  Y(KYJRE1,MS,1:NR) = JRE1(1:NR,MS)
 100  CONTINUE

      RETURN
      END
*DECK GETXY
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------COPY       PAMS VECTOR----A.B.   19.07.89---------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE GETXY(X,Y)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      COMPLEX*16    X(NXCOMP,MSMAX,*),Y(NYCOMP,MSMAX,*)
      INTEGER    MS,I
C
      DO 100 MS=1,MSMAX
      V1U(1:NRP1,MS) = X(KXV1,MS,1:NRP1)
      IF (KXX1.GT.0)  X1U(1:NRP1,MS) = X(KXX1,MS,1:NRP1)
      B1U(1:NTP1,MS) = X(KXB1 ,MS,1:NTP1)
      J2U(1:NTP1,MS) = X(KXJ2U,MS,1:NTP1)
      IF (KXJ2L.GT.0) J2L(1:NTP1,MS) = X(KXJ2L,MS,1:NTP1)
      J3U(1:NTP1,MS) = X(KXJ3,MS,1:NTP1)
      IF (KXPD.GT.0)  PDE(1:NRP1,MS) = X(KXPD,MS,1:NRP1)
      IF (KXPED.GT.0) PED(1:NRP1,MS) = X(KXPED,MS,1:NRP1)
      IF (ABS(RM(MS,2)).LT.0.1) THEN
         B1U(1:NRP1,MS) = X(KXB1,MS,1:NRP1)*T(1:NRP1)
         B1U(NRP1+1:NTP1,MS) = X(KXB1,MS,NRP1+1:NTP1)*T(NRP1)
      ENDIF
      IF (KXJRE.GT.0)   JRE(1:NRP1,MS)   = X(KXJRE,MS,1:NRP1)
      IF (KXJRE2.GT.0)  JRE2(1:NRP1,MS)  = X(KXJRE2,MS,1:NRP1)
      IF (KXJRE3.GT.0)  JRE3(1:NRP1,MS)  = X(KXJRE3,MS,1:NRP1)
      IF (KXJRE2L.GT.0) JRE2L(1:NRP1,MS) = X(KXJRE2L,MS,1:NRP1)
      IF (KXDPHI.GT.0)  DPHI(1:NRP1,MS)  = X(KXDPHI,MS,1:NRP1)

      V2U(1:NR,MS) = Y(KYV2,MS,1:NR)
      IF (KYV3.GT.0)    V3U(1:NR,MS)  = Y(KYV3,MS,1:NR)
      IF (KYX2.GT.0)    X2U(1:NR,MS)  = Y(KYX2,MS,1:NR)
      IF (KYX3.GT.0)    X3U(1:NR,MS)  = Y(KYX3,MS,1:NR)
      B2U(1:NTOT,MS) = Y(KYB2,MS,1:NTOT)
      B3U(1:NTOT,MS) = Y(KYB3,MS,1:NTOT)
      J1U(1:NTOT,MS) = Y(KYJ1,MS,1:NTOT)
      PRE(1:NR,MS)   = Y(KYPR,MS,1:NR)
      IF (KYPE.GT.0)    PEE(1:NR,MS)   = Y(KYPE,MS,1:NR)
      IF (KYPP.GT.0)    PEP(1:NR,MS)   = Y(KYPP,MS,1:NR)
      IF (KYPPERP.GT.0) PPERP(1:NR,MS) = Y(KYPPERP,MS,1:NR)
      IF (KYPPARA.GT.0) PPARA(1:NR,MS) = Y(KYPPARA,MS,1:NR)
      IF (KYRHOP.GT.0)  RHOP(1:NR,MS)  = Y(KYRHOP,MS,1:NR)
      IF (KYJRE1.GT.0)  JRE1(1:NR,MS)  = Y(KYJRE1,MS,1:NR)
 100  CONTINUE

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE ADDITIONAL PERTURBED QUANTITIES THAT ARE NOT DIRECTLY SOLVED
C FROM THE EQUATIONS
C VALID ONLY FOR EIGENVALUE PROBLEM
C YQ LIU, 2013-09
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GETXYMORE
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      IMPLICIT NONE
      INCLUDE    'compam.inc'
      INTEGER    I,J,MS,KCHECK
      REAL*8     HCHI
      COMPLEX*16 CTMP1
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::OX1,OX2,OX3

      KCHECK = 0

      ALLOCATE( OX1(NR,NCHI), OX2(NR,NCHI), OX3(NR,NCHI) )

C     COMPUTE X1U
      IF (KXX1.LE.0) THEN
         DO I=1,NRP1
         CTMP1 = ALNORM + RNTOR*ROT(I)*CI
         DO MS=1,MSMAX
            X1U(I,MS) = V1U(I,MS)/CTMP1
         ENDDO
         ENDDO
      ENDIF

C     GET X1U IN REAL SPACE
      IF ((KYX2.LE.0.OR.KYX3.LE.0).AND.ABS(ROTE).GT.0.) THEN
      OX1 = (0.,0.)
      HCHI = 2.*PI/DFLOAT(NCHI)
      DO J=1,NCHI
      DO MS=1,MSMAX
         CTMP1 = EXP(CI*RM(MS,2)*(J-1)*HCHI)
         DO I=1,NR
            OX1(I,J) = OX1(I,J) + (X1U(I,MS)+X1U(I+1,MS))*.5*CTMP1
         ENDDO
      ENDDO
      ENDDO
      ENDIF

C     COMPUTE X2U
      IF (KYX2.LE.0) THEN
      X2U = (0.,0.)
   
      IF (ABS(ROTE).GT.0.) THEN
      DO J=1,NCHI
      DO I=1,NR
         OX2(I,J) = OX1(I,J)/RJAM(I,J)
      ENDDO
      ENDDO
      
      DO J=1,NCHI
      DO MS=1,MSMAX
         CTMP1 = EXP(-CI*RM(MS,2)*(J-1)*HCHI)
         DO I=1,NR
            X2U(I,MS) = X2U(I,MS) + OX2(I,J)*CTMP1
         ENDDO
      ENDDO
      ENDDO

      HCHI = 1./DFLOAT(NCHI)
      DO I=1,NR
      CTMP1 = ALNORM + RNTOR*ROT(I)*CI
      DO MS=1,MSMAX
         X2U(I,MS) = -X2U(I,MS)*DPSIDSM(I)*DROTM(I)*HCHI/CTMP1
      ENDDO
      ENDDO
      ENDIF
      
      DO I=1,NR
      CTMP1 = ALNORM + RNTOR*ROT(I)*CI
      DO MS=1,MSMAX
         X2U(I,MS) = X2U(I,MS) + V2U(I,MS)/CTMP1
      ENDDO
      ENDDO

      ENDIF

C     COMPUTE X3U
      IF (KYX3.LE.0) THEN
      X3U = (0.,0.)
   
      IF (ABS(ROTE).GT.0.) THEN
      DO J=1,NCHI
      DO I=1,NR
         OX3(I,J) = OX1(I,J)/(G22LM(I,J)*DPSIDSM(I)**2/RJAM(I,J)**2 +
     &                        TM(I)**2/REQM(I,J)**2)
      ENDDO
      ENDDO
      
      DO J=1,NCHI
      DO MS=1,MSMAX
         CTMP1 = EXP(-CI*RM(MS,2)*(J-1)*HCHI)
         DO I=1,NR
            X3U(I,MS) = X3U(I,MS) + OX3(I,J)*CTMP1
         ENDDO
      ENDDO
      ENDDO

      HCHI = 1./DFLOAT(NCHI)
      DO I=1,NR
      CTMP1 = ALNORM + RNTOR*ROT(I)*CI
      DO MS=1,MSMAX
         X3U(I,MS) = X3U(I,MS)*TM(I)*DROTM(I)*HCHI/CTMP1
      ENDDO
      ENDDO
      ENDIF
      
      DO I=1,NR
      CTMP1 = ALNORM + RNTOR*ROT(I)*CI
      DO MS=1,MSMAX
         X3U(I,MS) = X3U(I,MS) + V3U(I,MS)/CTMP1
      ENDDO
      ENDDO

      ENDIF

      DEALLOCATE( OX1,OX2,OX3 )

      RETURN
      END

*DECK DPWALL
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------DELTA-PRIMES ON WALL -----------------------------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE DPWALL
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      INTEGER    MS,MS0,MSS,I,IV
      COMPLEX*16    DPRIM,GREXP,DERLOG(3)
C
      I = IWALLJ
      IV = NR + IWALLJ
      WRITE(*,1000)
      DO 10 MS=1,MSMAX
      DPRIM =((B1U(IV+1,MS)-B1U(IV,MS))/(VCS(I+1)-VCS(I))
     &       -(B1U(IV,MS)-B1U(IV-1,MS))/(VCS(I)-VCS(I-1)) )/B1U(IV,MS)
C
C     DELTA-PRIME IN TERMS OF S
C
      GREXP = DPRIM * VCS(I) / TAUWJ
CYQL  WRITE(*,1010) INT(RM(MS,2)),DPRIM,GREXP,ABS(B1U(IV,MS))
 10   CONTINUE
C
      DO 20 MS = 1,MSMAX
      IF (RM(MS,2).EQ.2.) GOTO 30
 20   CONTINUE
      RETURN
 30   MS0 = MS
      DO 50 I = 1,NV-1
      IV = NR + I
      DO 40 MSS = 1,3
      MS = MS0 - 1 + MSS
      DERLOG(MSS) = (B1U(IV+1,MS)-B1U(IV,MS))/(VCS(I+1)-VCS(I))
     &              *2./(B1U(IV,MS)+B1U(IV+1,MS))
 40   CONTINUE
CYQL  WRITE(*,1020) I,(DERLOG(MSS),MSS=1,3)
 50   CONTINUE
C
      RETURN
 1000 FORMAT('  M',4X,'RE(DP)',4X,'IM(DP)',8X,'GAMMA',8X,'OMEGA')
 1010 FORMAT(I3,2F10.3,4X,2E13.4,E14.3)
 1020 FORMAT(I3,3('  ',2F10.4))
      END
*DECK OUTPUT
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------OUTPUT EIGENMODE----------------- 05/04/93 --------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE OUTPUT(ISW)
C     =================
C
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      USE TORQUEM
      USE RCOMDM
      USE GIJLM 
      USE GVACUUMM
      USE ToolBox
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INCLUDE 'newrun.inc'
C
      INTEGER      MS,I,K,ISW,KSW
      CHARACTER*80 LINE
      CHARACTER(LEN=1024) FILENAME
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
      INTEGER RESULTFILE_ID
C
      KSW = 1
      IF (NCASE.NE.6.AND.NCAE.NE.10) KSW = -2
      IF ((NCASE.EQ.6.OR.NCASE.EQ.10).AND.
     &    OSWEEP.GE.1.AND.OSWEEP.LE.NSWEEP) KSW = MOD(ISW,OSWEEP)
      IF ((NCASE.EQ.6.OR.NCASE.EQ.10).AND.ISW.EQ.NSWEEP) KSW = -1
      
      IF (.NOT.OUTDAT) GOTO 20
      OPEN(CHOUTP,FILE='OUTDATA.OUT',FORM='FORMATTED')
      REWIND(CHOUTP)
      WRITE(CHOUTP,BASIC)
      WRITE(CHOUTP,FEEDBACK)
      WRITE(CHOUTP,KINETIC)
      WRITE(CHOUTP,QLIN)
      WRITE(CHOUTP,NUMERIC)
      WRITE(CHOUTP,OUTOPT)
      WRITE(CHOUTP,1000) NRP1,MSMAX,NSMAX
      WRITE(CHOUTP,1010) ALNORM
      DO 10 MS=1,MSMAX
 10   WRITE(CHOUTP,1020) MS,RM(MS,2),RN(2)
      WRITE(CHOUTP,1030) (CSV(I),I=1,NRP1)
C
      CALL GENOUT(V1U,NTP1,'   V1U',CHOUTP,NRP1  ,MSMAX,2)
      CALL GENOUT(V2U,NTP1,'   V2U',CHOUTP,NR    ,MSMAX,2)
      CALL GENOUT(V3U,NTP1,'   V3U',CHOUTP,NR    ,MSMAX,2)
      CALL GENOUT(B1U,NTP1,'   B1U',CHOUTP,NTP1  ,MSMAX,2)
      CALL GENOUT(B2U,NTP1,'   B2U',CHOUTP,NTP1-1,MSMAX,2)
      CALL GENOUT(B3U,NTP1,'   B3U',CHOUTP,NTP1-1,MSMAX,2)
      CALL GENOUT(J1U,NTP1,'   J1U',CHOUTP,NTP1-1,MSMAX,2)
      CALL GENOUT(J2U,NTP1,'   J2U',CHOUTP,NTP1  ,MSMAX,2)
      CALL GENOUT(J3U,NTP1,'   J3U',CHOUTP,NTP1  ,MSMAX,2)
      CALL GENOUT(PRE,NTP1,'   PRE',CHOUTP,NR    ,MSMAX,2)
      CALL GENOUT(PEE,NTP1,'   PEE',CHOUTP,NR    ,MSMAX,2)
      CALL GENOUT(PEP,NTP1,'   PEP',CHOUTP,NR    ,MSMAX,2)
      CALL GENOUT(PDE,NTP1,'   PDE',CHOUTP,NRP1  ,MSMAX,2)
      CALL GENOUT(PED,NTP1,'   PED',CHOUTP,NRP1  ,MSMAX,2)
      CALL GENOUT(PPERP,NTP1,' PPERP',CHOUTP,NR  ,MSMAX,2)
      CALL GENOUT(PPARA,NTP1,' PPARA',CHOUTP,NR  ,MSMAX,2)
      CALL GENOUT(X1U,NTP1,'   X1U',CHOUTP,NRP1  ,MSMAX,2)
      CALL GENOUT(X2U,NTP1,'   X2U',CHOUTP,NR    ,MSMAX,2)
      CALL GENOUT(X3U,NTP1,'   X3U',CHOUTP,NR    ,MSMAX,2)
C
      CLOSE(CHOUTP)
C
C.. OUTPUT DENSITY AND RESISTIVITY PROFILES AGAINST CS AND CSV
C
      IF (.NOT.LPROFI) GOTO 20
      OPEN(CHMAP,FILE='MAP.OUT',FORM='FORMATTED')
      REWIND(CHMAP)
      WRITE(CHMAP,1050) NRP1,(CS(I),CSV(I),RHO(I),RESIST(I),I=1,NRP1)
      CLOSE(CHMAP)
      WRITE(*,1050)     NRP1,(CS(I),CSV(I),RHO(I),RESIST(I),I=1,NRP1)
C
 20   CONTINUE
C----------------------------------------
CYQLIU 15/04/1999
      IF (INCFEED.GE.0) CALL FEEDOUT

CYQLIU 03/05/2011
      IF ((CALPHA7.EQ.0..AND.INCFEED.NE.4).OR.
     &    (INCFEED.EQ.4.AND.KEYTORQ.EQ.2)) THEN
         IF (ABS(CTJXB).GT.0.) CALL TORQJXB
         IF (ABS(CTNTV).GT.0..OR.ABS(CDNTV).GT.0.) CALL TORQNTV
         IF (ABS(CTREY).GT.0..OR.ABS(CDMHD).GT.0.) CALL TORQREY
         IF (ABS(CTERGO).GT.0..OR.ABS(CDERGO).GT.0.) CALL TORQERGO
         CALL CALCDISPNORM
         IF (NADAPS.GT.0) CALL GETSNEWMESH
      ENDIF

CYQLIU 12/05/2003
      IF (ABS(PVISC).GT.0..AND.IVISC.EQ.1) CALL PPVISC
C----------------------------------------
CZRWANG FOR SHELL CODE 11/2010
      RESULTFILE_ID = assignFreeFileUnit ()
      OPEN(RESULTFILE_ID,FILE='SHELLRESULT.OUT',STATUS='REPLACE')
      WRITE (RESULTFILE_ID,*) REAL(ALNORM),DIMAG(ALNORM),NIT
      
      CLOSE (RESULTFILE_ID)
      

      OPEN(CHLIST,FILE='RESULT.OUT',FORM='FORMATTED')
      IF (NPARAM.LE.1) GOTO 40
      DO 30 I = 2,NPARAM
      READ(CHLIST,'(A)',END=40) LINE
 30   CONTINUE
 40   CONTINUE
      IF (INCFEED.EQ.1.OR.INCFEED.EQ.11.OR.INCFEED.EQ.12.OR.
     &    INCFEED.EQ.3.OR.INCFEED.EQ.5.OR.INCFEED.EQ.9) THEN
         WRITE(CHLIST,1045) GAINA(1),GAINP(1),AL0,ALNORM
      ENDIF
      IF (INCFEED.EQ.0) THEN
         WRITE(CHLIST,1045) GAINA(1),DFLOAT(NIT),AL0,ALNORM
      ENDIF
      IF (INCFEED.EQ.2.OR.INCFEED.EQ.4.OR.INCFEED.EQ.10) THEN
         WRITE(CHLIST,1045) RLK(1),RMFS(1),AL0,ALNORM
      ENDIF
      CLOSE(CHLIST)

C     OUTPUT FEEDBACK CURRENTS AND SENSOR SIGNALS
      IF (KSW.LT.0) THEN
      OPEN(CHLIST,FILE='CURRSENS.OUT',FORM='FORMATTED')
      IF (NPARAM.LE.1) GOTO 60
      DO 50 I = 2,NPARAM
         DO K=1,NCOIL
            READ(CHLIST,'(A)',END=60) LINE
         ENDDO
 50   CONTINUE
 60   CONTINUE
      DO K=1,NCOIL
         WRITE(CHLIST,1045) GAINA(K),GAINP(K),OFEEDI(K),
     &                      OSENSR(K),OSENSP(K)
      ENDDO
      CLOSE(CHLIST)
      ENDIF

C     OUTPUT EQUILIBRIUM METRICS 
      IF (KSW.LT.0) THEN
 120  FORMAT(17(E15.8,1X))
      OPEN(CHLIST,FILE='METRICS.OUT')
      DO I=1,NRP1
         WRITE(CHLIST,120) CS(I),
     &        DG11L(I,1),DG12L(I,1),DG22L(I,1),DG33L(I,1),
     &        DG11L(I,2),DG12L(I,2),DG22L(I,2),DG33L(I,2)
      ENDDO
      DO I=1,NVEQ1 
         WRITE(CHLIST,120) VCS(I),
     &        VG11L(I,1),VG12L(I,1),VG22L(I,1),VG33L(I,1),
     &        VG11L(I,2),VG12L(I,2),VG22L(I,2),VG33L(I,2)
      ENDDO
      CLOSE(CHLIST)
      ENDIF

C     --- OUTPUT THE PERTURBED PLASMA CURRENT, FOR CARIDDI INPUT
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='JPLASMA.OUT')
      ELSE
         WRITE(FILENAME,"(A8,I0.4,A4)") "JPLASMA_",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      IF (KJRER.GT.0.AND.KJRER.LT.6) THEN
         WRITE(CHOUTP,1171) MSMAX,NTP1,RNTOR,0,0,0,0,0
         DO MS=1,MSMAX
            WRITE(CHOUTP,1172) RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2),
     &                         RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2)
         ENDDO
         DO MS=1,MSMAX
         DO II=1,NRP1
            WRITE(CHOUTP,1172) REAL(J1U(II,MS)),IMAG(J1U(II,MS)),
     &                         REAL(J2U(II,MS)),IMAG(J2U(II,MS)),
     &                         REAL(J3U(II,MS)),IMAG(J3U(II,MS)),
     &                         REAL(JRE(II,MS)),IMAG(JRE(II,MS))
         ENDDO
         DO II=NRP1+1,NTP1
            WRITE(CHOUTP,1172) REAL(J1U(II,MS)),IMAG(J1U(II,MS)),
     &                         REAL(J2U(II,MS)),IMAG(J2U(II,MS)),
     &                         REAL(J3U(II,MS)),IMAG(J3U(II,MS)),
     &                         0.0,0.0
         ENDDO
         ENDDO
      ELSEIF (KJRER.EQ.6) THEN
         WRITE(CHOUTP,1171) MSMAX,NTP1,RNTOR,0,0,0,0,0,0,0,0,0,0,0
         DO MS=1,MSMAX
            WRITE(CHOUTP,1172) RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2),
     &                         RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2),
     &                         RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2),
     &                         RM(MS,2),RM(MS,2)
         ENDDO
         DO MS=1,MSMAX
         DO II=1,NRP1
            WRITE(CHOUTP,1172) REAL(J1U(II,MS)),  IMAG(J1U(II,MS)),
     &                         REAL(J2U(II,MS)),  IMAG(J2U(II,MS)),
     &                         REAL(J3U(II,MS)),  IMAG(J3U(II,MS)),
     &                         REAL(JRE1(II,MS)), IMAG(JRE1(II,MS)),
     &                         REAL(JRE2(II,MS)), IMAG(JRE2(II,MS)),
     &                         REAL(JRE3(II,MS)), IMAG(JRE3(II,MS)),
     &                         REAL(JRE(II,MS)),  IMAG(JRE(II,MS))
         ENDDO
         DO II=NRP1+1,NTP1
            WRITE(CHOUTP,1172) REAL(J1U(II,MS)), IMAG(J1U(II,MS)),
     &                         REAL(J2U(II,MS)), IMAG(J2U(II,MS)),
     &                         REAL(J3U(II,MS)), IMAG(J3U(II,MS)),
     &                         0.0,0.0,
     &                         0.0,0.0,
     &                         0.0,0.0,
     &                         0.0,0.0
         ENDDO
         ENDDO
      ELSE
         WRITE(CHOUTP,1171) MSMAX,NTP1,RNTOR,0,0,0
         DO MS=1,MSMAX
            WRITE(CHOUTP,1172) RM(MS,2),RM(MS,2),RM(MS,2),
     &                         RM(MS,2),RM(MS,2),RM(MS,2)
         ENDDO
         DO MS=1,MSMAX
         DO II=1,NTP1
            WRITE(CHOUTP,1172) REAL(J1U(II,MS)),IMAG(J1U(II,MS)),
     &                         REAL(J2U(II,MS)),IMAG(J2U(II,MS)),
     &                         REAL(J3U(II,MS)),IMAG(J3U(II,MS))
         ENDDO
         ENDDO
      ENDIF
      CLOSE(CHOUTP)
      ENDIF
 1171 FORMAT(I5,1X,I5,1X,E9.2,11(1X,I2))
 1172 FORMAT(14(E16.8E3,1X))

C     --- OUTPUT THE PERTURBED MAGNETIC FIELD, FOR VISUALIZATION
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='BPLASMA.OUT')
      ELSE
         WRITE(FILENAME,"(A8,I0.4,A4)") "BPLASMA_",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      WRITE(CHOUTP,1171) MSMAX,NTP1,RNTOR,0,0,0
      DO MS=1,MSMAX
         WRITE(CHOUTP,1172) RM(MS,2),RM(MS,2),RM(MS,2),
     &                      RM(MS,2),RM(MS,2),RM(MS,2)
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NTP1
            WRITE(CHOUTP,1172) REAL(B1U(II,MS)),IMAG(B1U(II,MS)),
     &                         REAL(B2U(II,MS)),IMAG(B2U(II,MS)),
     &                         REAL(B3U(II,MS)),IMAG(B3U(II,MS))
         ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF
      
C     OUTPUT ALL FIELD DATA (B,J,V) FOR RE TRACING STUDY     
C     WHEN PERTURBED FIELDS FOR SEVERAL TOROIDAL HARMONICS NEED TO BE
C     SUPERPOSED   
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='FIELD_RE.OUT')
      ELSE
         WRITE(FILENAME,"(A9,I0.4,A4)") "FIELD_RE_",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      WRITE(CHOUTP,1271) NR,NV,M1,M2,RNTOR,REAL(ALNORM),IMAG(ALNORM)
      DO II=1,NR+NV
         WRITE(CHOUTP,1172) CS(II),CSM(II)
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NR+NV
            WRITE(CHOUTP,1172) REAL(B1U(II,MS)),IMAG(B1U(II,MS)),
     &                         REAL(B2U(II,MS)),IMAG(B2U(II,MS)),
     &                         REAL(B3U(II,MS)),IMAG(B3U(II,MS))
         ENDDO
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NR+1
            WRITE(CHOUTP,1172) REAL(J1U(II,MS)),IMAG(J1U(II,MS)),
     &                         REAL(J2U(II,MS)),IMAG(J2U(II,MS)),
     &                         REAL(J3U(II,MS)),IMAG(J3U(II,MS))
         ENDDO
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NR+1
            WRITE(CHOUTP,1172) REAL(V1U(II,MS)),IMAG(V1U(II,MS)),
     &                         REAL(V2U(II,MS)),IMAG(V2U(II,MS)),
     &                         REAL(V3U(II,MS)),IMAG(V3U(II,MS))
         ENDDO
      ENDDO
      CLOSE(CHOUTP)
 1271 FORMAT(4(I5,1X),E9.2,1X,2(E14.7,1X))
      ENDIF
      
C     --- OUTPUT B1U(I,MS), MS=1,MSMAX AT COIL POSITION
      IF (KSW.LT.0) THEN
      OPEN(CHOUTP,FILE='BNORM01.OUT')
      REWIND(CHOUTP)
      II = NR + IFEED
      DO MS=1,MSMAX
         WRITE(CHOUTP,1172) DREAL(B1U(II,MS)),DIMAG(B1U(II,MS))
      ENDDO
      CLOSE(CHOUTP)
      ENDIF

C     --- OUTPUT THE PERTURBED VELOCITY, FOR VISUALIZATION
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='VPLASMA.OUT')
      ELSE
         WRITE(FILENAME,"(A8,I0.4,A4)") "VPLASMA_",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      WRITE(CHOUTP,1171) MSMAX,NRP1,RNTOR,0,0,0
      DO MS=1,MSMAX
C        IF (ABS(RM(MS,2)).GT.0.1.OR.ABS(RNTOR).GT.1.0E-13) THEN
         WRITE(CHOUTP,1172) RM(MS,2),RM(MS,2),RM(MS,2),
     &                      RM(MS,2),RM(MS,2),RM(MS,2)
C        ENDIF
      ENDDO
      DO II=1,NRP1
         WRITE(CHOUTP,1172) DPSIDS(II),DPSIDS(II),DPSIDS(II),
     &                      T(II),T(II),T(II)
      ENDDO
      DO MS=1,MSMAX
C        IF (ABS(RM(MS,2)).GT.0.1.OR.ABS(RNTOR).GT.1.0E-13) THEN
         DO II=1,NRP1
            WRITE(CHOUTP,1172) REAL(V1U(II,MS)),IMAG(V1U(II,MS)),
     &                         REAL(V2U(II,MS)),IMAG(V2U(II,MS)),
     &                         REAL(V3U(II,MS)),IMAG(V3U(II,MS))

         ENDDO
C        ENDIF
      ENDDO
      CLOSE(CHOUTP)
      ENDIF
      
C     --- OUTPUT THE PERTURBED VELOCITY, FOR VISUALIZATION
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='XPLASMA.OUT')
      ELSE
         WRITE(FILENAME,"(A8,I0.4,A4)") "XPLASMA_",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      WRITE(CHOUTP,1171) MSMAX,NRP1,RNTOR,0,0,0
      DO MS=1,MSMAX
         WRITE(CHOUTP,1172) RM(MS,2),RM(MS,2),RM(MS,2),
     &                      RM(MS,2),RM(MS,2),RM(MS,2)
      ENDDO
      DO II=1,NRP1
         WRITE(CHOUTP,1172) DPSIDS(II),DPSIDS(II),DPSIDS(II),
     &                      T(II),T(II),T(II)
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NRP1
            WRITE(CHOUTP,1172) REAL(X1U(II,MS)),IMAG(X1U(II,MS)),
     &                         REAL(X2U(II,MS)),IMAG(X2U(II,MS)),
     &                         REAL(X3U(II,MS)),IMAG(X3U(II,MS))
         ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF

C     --- OUTPUT THE PERTURBED DENSITY, FOR VISUALIZATION
      IF (KSW.LT.0) THEN
      OPEN(CHOUTP,FILE='DPLASMA.OUT')
      REWIND(CHOUTP)
      WRITE(CHOUTP,11710) MSMAX,NRP1,0,0
      DO MS=1,MSMAX
         WRITE(CHOUTP,11720) RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2)
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NRP1
            WRITE(CHOUTP,11720) REAL(RHOP(II,MS)),IMAG(RHOP(II,MS)),
     &                         REAL(PRE(II,MS)),IMAG(PRE(II,MS))
         ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF
11710 FORMAT(I5,1X,I5,2(1X,I2))
11720 FORMAT(E16.8E3,3(1X,E16.8E3))
 
C     --- OUTPUT THE PERTURBED PLASMA PRESSURES, FOR VISUALIZATION
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='PPLASMA.OUT')
      ELSE
         WRITE(FILENAME,"(A8,I0.4,A4)") "PPLASMA_",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      OPEN(CHOUTP,FILE='PPLASMA.OUT')
      REWIND(CHOUTP)
      WRITE(CHOUTP,1175) MSMAX,NRP1,RNTOR,0,0,0,0,0,0,0,0,0
      DO MS=1,MSMAX
         WRITE(CHOUTP,11721) RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2),
     &                       RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2),
     &                       RM(MS,2),RM(MS,2),RM(MS,2),RM(MS,2)
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NRP1
            WRITE(CHOUTP,11721) REAL(PRE(II,MS)),IMAG(PRE(II,MS)),
     &                          REAL(PEE(II,MS)),IMAG(PEE(II,MS)),
     &                          REAL(PEP(II,MS)),IMAG(PEP(II,MS)),
     &                          REAL(PPERP(II,MS)),IMAG(PPERP(II,MS)),
     &                          REAL(PPARA(II,MS)),IMAG(PPARA(II,MS)),
     &                          REAL(DPHI(II,MS)),IMAG(DPHI(II,MS))
         ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF 
 1175 FORMAT(I5,1X,I5,1X,E9.2,9(1X,I2))
11721 FORMAT(F30.15,11(1X,F30.15))

C     --- OUTPUT THE JACOBIAN
      IF (KSW.LT.0) THEN
      OPEN(CHOUTP,FILE='JACOBIAN.OUT')
      REWIND(CHOUTP)
      WRITE(CHOUTP,1171) MSMAX,NRP1
      DO MS=1,MSMAX
         WRITE(CHOUTP,1172) RM(MS,2),RM(MS,2)
      ENDDO
      DO MS=1,MSMAX
         DO II=1,NRP1
            WRITE(CHOUTP,1172) REAL(JACOBI(II,MS)),
     &                         IMAG(JACOBI(II,MS))
         ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF

      RETURN

 1000 FORMAT(3I20)
 1010 FORMAT(///,2D30.15,//)
 1020 FORMAT(I20,2F20.5)
 1030 FORMAT(///,(2D30.15))
 1040 FORMAT(1P,3X,E13.5,3X,2E13.5,3X,2E13.5)
 1045 FORMAT(1P,E10.2,E13.5,2E13.5,2E13.5,2E13.5)
 1050 FORMAT(I5,/,(3F20.10,E20.5))
      END
*DECK BOWALL
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--BOUNDARY COND. FOR CONDUCTING WALL --------- A. BONDESON 24.05.90 ---
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE BOWALL(
C     ==================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'
      INTEGER     MSROW,LXROW
C
C.....SET V1, X1 AND B1 EQUAL TO 0 ON LAST POINT NTP1
C
      DO 10 MSROW = 1,MSMAX

      LXROW = (MSROW-1)*NXCOMP
      CALL ANNIHX(KXV1,MSROW,NTP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXV1+LXROW,KXV1+LXROW,NTP1) = 1.
      IF (KXX1.GT.0) THEN
      CALL ANNIHX(KXX1,MSROW,NTP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXX1+LXROW,KXX1+LXROW,NTP1) = 1.
      ENDIF
      CALL ANNIHX(KXB1,MSROW,NTP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXB1+LXROW,KXB1+LXROW,NTP1) = 1.
      IF (KXW1.GT.0) THEN
      CALL ANNIHX(KXW1,MSROW,NTP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXW1+LXROW,KXW1+LXROW,NTP1) = 1.
      ENDIF
      IF (KXJRE.GT.0) THEN
      CALL ANNIHX(KXJRE,MSROW,NTP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXJRE+LXROW,KXJRE+LXROW,NTP1) = 1.
      ENDIF
      IF (KXB2L.GT.0) THEN
      CALL ANNIHX(KXB2L,MSROW,NTP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXB2L+LXROW,KXB2L+LXROW,NTP1) = 1.
      ENDIF
      IF (KXB3L.GT.0) THEN
      CALL ANNIHX(KXB3L,MSROW,NTP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXB3L+LXROW,KXB3L+LXROW,NTP1) = 1.
      ENDIF

      IF (ABS(RNTOR).LT.1.E-10.AND.
     &    (K_BC_N0.EQ.11.OR.K_BC_N0.EQ.21.OR.K_BC_N0.EQ.31)) THEN
      CALL ANNIHX(KXJ2U,MSROW,NRP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXJ2U+LXROW,KXJ2U+LXROW,NRP1) = 1.
      CALL ANNIHX(KXJ3,MSROW,NRP1,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      BSUBM(KXJ3+LXROW,KXJ3+LXROW,NRP1) = 1.
      ENDIF

 10   CONTINUE

      RETURN
      END
*DECK BOUNDB
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--BOUNDARY COND. FOR BACKWARD COUPLING --------- YQ LIU, AUG.4, 2008 --
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE BOUNDB(
C     ==================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE FEEDBACKM
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
      INTEGER  MSA, MSB, L, MSROW, LXCOL, LYCOL, M16

      M16 = MSMAX*16
      IF (INCFEED.EQ.0) THEN

C     ALLOCATE COUPLING MATRICES
      IF (.NOT. ALLOCATED(A01)) THEN
         ALLOCATE( A01(MSMAX,MSMAX), A02(MSMAX,MSMAX),
     &             A11(MSMAX,MSMAX), A12(MSMAX,MSMAX),
     &             B01(MSMAX,MSMAX), B02(MSMAX,MSMAX),
     &             B11(MSMAX,MSMAX), B12(MSMAX,MSMAX) )
         ALLOCATE( A00(M16) )
      END IF

C     READ IN COUPLING MATRICES 
      OPEN(31,FILE='BWCMAT.IN',FORM='FORMATTED',STATUS='OLD')
      DO MSA=1,MSMAX
         READ(31,*) (A00(L),L=1,M16)
         DO MSB=1,MSMAX
            A01(MSA,MSB)=A00(MSB+0*MSMAX)+A00(MSB+1*MSMAX)*CI
            A11(MSA,MSB)=A00(MSB+2*MSMAX)+A00(MSB+3*MSMAX)*CI
            B01(MSA,MSB)=A00(MSB+4*MSMAX)+A00(MSB+5*MSMAX)*CI
            B11(MSA,MSB)=A00(MSB+6*MSMAX)+A00(MSB+7*MSMAX)*CI
            A02(MSA,MSB)=A00(MSB+8*MSMAX)+A00(MSB+9*MSMAX)*CI
            A12(MSA,MSB)=A00(MSB+10*MSMAX)+A00(MSB+11*MSMAX)*CI
            B02(MSA,MSB)=A00(MSB+12*MSMAX)+A00(MSB+13*MSMAX)*CI
            B12(MSA,MSB)=A00(MSB+14*MSMAX)+A00(MSB+15*MSMAX)*CI
         ENDDO
      ENDDO
      CLOSE(31)

      ELSE

C     ALLOCATE COUPLING MATRICES
      IF (.NOT. ALLOCATED(A01)) THEN
         ALLOCATE( A01(MSMAX,MSMAX),    A02(MSMAX,MSMAX),
     &             A11(MSMAX,MSMAX),    A12(MSMAX,MSMAX),
     &             B01(MSMAX,MSMAX),    B02(MSMAX,MSMAX),
     &             B11(MSMAX,MSMAX),    B12(MSMAX,MSMAX),
     &             C01(MSMAX,NCOILT),   C02(MSMAX,NCOILT),
     &             C11(MSMAX,NCOILT),   C12(MSMAX,NCOILT),
     &             A01S(NSENST,MSMAX),  A02S(NSENST,MSMAX),
     &             A11S(NSENST,MSMAX),  A12S(NSENST,MSMAX),
     &             B01S(NSENST,MSMAX),  B02S(NSENST,MSMAX),
     &             B11S(NSENST,MSMAX),  B12S(NSENST,MSMAX),
     &             C01S(NSENST,NCOILT), C02S(NSENST,NCOILT),
     &             C11S(NSENST,NCOILT), C12S(NSENST,NCOILT) )
         ALLOCATE( A00(M16+NCOILT*8) )
      END IF

C     READ IN COUPLING MATRICES 
      OPEN(31,FILE='BWCMAT.IN',FORM='FORMATTED',STATUS='OLD')
      DO MSA=1,MSMAX
         READ(31,*) (A00(L),L=1,M16+NCOILT*8)
         DO MSB=1,MSMAX
            A01(MSA,MSB)=A00(MSB+0*MSMAX)+A00(MSB+1*MSMAX)*CI
            A11(MSA,MSB)=A00(MSB+2*MSMAX)+A00(MSB+3*MSMAX)*CI
            B01(MSA,MSB)=A00(MSB+4*MSMAX)+A00(MSB+5*MSMAX)*CI
            B11(MSA,MSB)=A00(MSB+6*MSMAX)+A00(MSB+7*MSMAX)*CI
            A02(MSA,MSB)=A00(MSB+8*MSMAX)+A00(MSB+9*MSMAX)*CI
            A12(MSA,MSB)=A00(MSB+10*MSMAX)+A00(MSB+11*MSMAX)*CI
            B02(MSA,MSB)=A00(MSB+12*MSMAX)+A00(MSB+13*MSMAX)*CI
            B12(MSA,MSB)=A00(MSB+14*MSMAX)+A00(MSB+15*MSMAX)*CI
         ENDDO
         DO MSB=1,NCOILT
            C01(MSA,MSB)=A00(MSB+M16+0*NCOILT)+A00(MSB+M16+1*NCOILT)*CI
            C11(MSA,MSB)=A00(MSB+M16+2*NCOILT)+A00(MSB+M16+3*NCOILT)*CI
            C02(MSA,MSB)=A00(MSB+M16+4*NCOILT)+A00(MSB+M16+5*NCOILT)*CI
            C12(MSA,MSB)=A00(MSB+M16+6*NCOILT)+A00(MSB+M16+7*NCOILT)*CI
         ENDDO
      ENDDO
      CLOSE(31)

      OPEN(31,FILE='BWCMATS.IN',FORM='FORMATTED',STATUS='OLD')
      DO MSA=1,NSENST
         READ(31,*) (A00(L),L=1,M16+NCOILT*8)
         DO MSB=1,MSMAX
            A01S(MSA,MSB)=A00(MSB+0*MSMAX)+A00(MSB+1*MSMAX)*CI
            A11S(MSA,MSB)=A00(MSB+2*MSMAX)+A00(MSB+3*MSMAX)*CI
            B01S(MSA,MSB)=A00(MSB+4*MSMAX)+A00(MSB+5*MSMAX)*CI
            B11S(MSA,MSB)=A00(MSB+6*MSMAX)+A00(MSB+7*MSMAX)*CI
            A02S(MSA,MSB)=A00(MSB+8*MSMAX)+A00(MSB+9*MSMAX)*CI
            A12S(MSA,MSB)=A00(MSB+10*MSMAX)+A00(MSB+11*MSMAX)*CI
            B02S(MSA,MSB)=A00(MSB+12*MSMAX)+A00(MSB+13*MSMAX)*CI
            B12S(MSA,MSB)=A00(MSB+14*MSMAX)+A00(MSB+15*MSMAX)*CI
         ENDDO
         DO MSB=1,NCOILT
            C01S(MSA,MSB)=A00(MSB+M16+0*NCOILT)+A00(MSB+M16+1*NCOILT)*CI
            C11S(MSA,MSB)=A00(MSB+M16+2*NCOILT)+A00(MSB+M16+3*NCOILT)*CI
            C02S(MSA,MSB)=A00(MSB+M16+4*NCOILT)+A00(MSB+M16+5*NCOILT)*CI
            C12S(MSA,MSB)=A00(MSB+M16+6*NCOILT)+A00(MSB+M16+7*NCOILT)*CI
         ENDDO
      ENDDO
      CLOSE(31)

      ENDIF

C     FILL IN MATRICES FOR BC AT NTP1
      IF (NCOUPL.EQ.-2) THEN
      DO MSA=1,MSMAX
         MSROW=(MSA-1)*NXCOMP
         DO MSB=1,MSMAX
            LXCOL=(MSB-1)*NXCOMP
            LYCOL=(MSB-1)*NYCOMP
            IF (ABS(RM(MSB,2)).GT.0.1) THEN
            BSUBM(KXB1+MSROW,KXB1+LXCOL,NTP1)=
     &           B01(MSA,MSB)+AL0*B11(MSA,MSB)
            HSUBM(KXB1+MSROW,KYB2+LYCOL,NTP1)=
     &           -A01(MSA,MSB)-AL0*A11(MSA,MSB)
            ELSE
            BSUBM(KXB1+MSROW,KXB1+LXCOL,NTP1)=
     &           (B01(MSA,MSB)+AL0*B11(MSA,MSB))*T(NRP1)
            HSUBM(KXB1+MSROW,KYB3+LYCOL,NTP1)=
     &           -A01(MSA,MSB)-AL0*A11(MSA,MSB)
            ENDIF
         ENDDO
      ENDDO
      ENDIF

      IF (NCOUPL.EQ.-3) THEN
      DO MSA=1,MSMAX
         MSROW=(MSA-1)*NXCOMP
         DO MSB=1,MSMAX
            LXCOL=(MSB-1)*NXCOMP
            LYCOL=(MSB-1)*NYCOMP
            IF (ABS(RM(MSB,2)).GT.0.1) THEN
               BSUBM(KXB1+MSROW,KXB1+LXCOL,NTP1)=
     &              B02(MSA,MSB)+AL0*B12(MSA,MSB)
            ELSE
               BSUBM(KXB1+MSROW,KXB1+LXCOL,NTP1)=
     &             (B02(MSA,MSB)+AL0*B12(MSA,MSB))*T(NRP1)
            ENDIF
            HSUBM(KXB1+MSROW,KYB3+LYCOL,NTP1)=
     &           -A02(MSA,MSB)-AL0*A12(MSA,MSB)
         ENDDO
      ENDDO
      ENDIF

      RETURN
      END
*DECK BOVACU02
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--BOUNDARY COND. FOR COUPLING TO CARIDDI ----- Y.Q. LIU    18.05.2005--
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE BOVACU02(
C     ====================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'
      INCLUDE 'comioc.inc'

      INTEGER     MSROW,LXROW,MSCOL,LXCOL
      REAL*8      TMP1,TMP2
      COMPLEX*16  BNORM(100)
C
C.....SPECIAL BOUNDARY CONDITION FOR B1
C
      WRITE(*,*) 'READ BNORM02'
      OPEN(CHOUTP,FILE='BNORM02.IN',FORM='FORMATTED')
      DO MSROW = 1,MSMAX
         READ(CHOUTP,*) TMP1,TMP2
         BNORM(MSROW) = TMP1 + TMP2*CI
      ENDDO
      CLOSE(CHOUTP)
      
      DO MSROW = 1,MSMAX
      LXROW = (MSROW-1)*NXCOMP
      BSUBM(KXB1+LXROW,KXB1+LXROW,NTP1) = 1.
      DO MSCOL = 1,1
      LXCOL = (MSCOL-1)*NXCOMP
      BSUBM(KXB1+LXROW,KXB1+LXCOL,NTP1) = 
     &   BSUBM(KXB1+LXROW,KXB1+LXCOL,NTP1) - 
     &   BNORM(MSROW)
      ENDDO
      ENDDO
 
      RETURN
      END
*DECK BOVACU03
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--BOUNDARY COND. FOR COUPLING TO CARIDDI ----- Y.Q. LIU    18.05.2005--
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE BOVACU03(
C     ====================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'
      INCLUDE 'comioc.inc'

      INTEGER     MSROW,LXROW,LYCOL
      REAL*8      TMP1,TMP2
      COMPLEX*16  BNORM(100)
C
C.....SPECIAL BOUNDARY CONDITION FOR B1
C
      WRITE(*,*) 'READ BNORM03'
      OPEN(CHOUTP,FILE='BNORM03.IN',FORM='FORMATTED')
      DO MSROW = 1,MSMAX
         READ(CHOUTP,*) TMP1,TMP2
         BNORM(MSROW) = TMP1 + TMP2*CI
      ENDDO
      CLOSE(CHOUTP)
      
      DO MSROW = 1,MSMAX
      LXROW = (MSROW-1)*NXCOMP
      BSUBM(KXB1+LXROW,KXB1+LXROW,NTP1) = 1.
      LYCOL = (MSROW-1)*NYCOMP
      HSUBM(KXB1+LXROW,KYB2+LYCOL,NTP1) = -BNORM(MSROW)
      ENDDO
 
      RETURN
      END
*DECK BOVACU03
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--BOUNDARY COND. FOR COUPLING TO CARIDDI ----- Y.Q. LIU    18.05.2005--
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE BOVACU04(
C     ====================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'
      INCLUDE 'comioc.inc'

      INTEGER     MSROW,LXROW,MSCOL,LYCOL
      REAL*8      TMP1,TMP2
      COMPLEX*16  BNORM(100)
C
C.....SPECIAL BOUNDARY CONDITION FOR B1
C
      WRITE(*,*) 'READ BNORM04'
      OPEN(CHOUTP,FILE='BNORM04.IN',FORM='FORMATTED')
      DO MSROW = 1,MSMAX
         READ(CHOUTP,*) TMP1,TMP2
         BNORM(MSROW) = TMP1 + TMP2*CI
      ENDDO
      CLOSE(CHOUTP)
      
      DO MSROW = 1,MSMAX
      LXROW = (MSROW-1)*NXCOMP
      BSUBM(KXB1+LXROW,KXB1+LXROW,NTP1) = 1.
      DO MSCOL = 1,MSMAX
      LYCOL = (MSCOL-1)*NYCOMP
      HSUBM(KXB1+LXROW,KYB2+LYCOL,NTP1) = -BNORM(MSROW)
      ENDDO
      ENDDO
 
      RETURN
      END
*DECK BOVACU
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--BOUNDARY COND. FOR VACUUM CALCULATION------- A. BONDESON 24.05.90 ---
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE BOVACU(
C     ==================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'
      INTEGER      MS,MM
C
C.....SET B1U=0 ON CONDUCTING WALL
C
C     CALL VZERO(2*MXMAX*MXMAX,BSUBM(1,1,NTP1))
C     CALL VZERO(2*MXMAX*MXMAX,ASUBM(1,1,NTP1))
C     CALL VZERO(2*MXMAX*MYMAX,HSUBM(1,1,NTP1))
C
      DO 10 MS=1,MSMAX
      MM = KXB1 + (MS-1)*NXCOMP
 10   BSUBM(MM,MM,NTP1) = 1.
C
      RETURN
      END
*DECK PRESET
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------PRESET ROUTINE-------------A.B.   31/01/90 --------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....PRESET NAMELIST VARIABLES.......................................
C
      SUBROUTINE PRESET
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE REORBITM
      INCLUDE 'compam.inc'
      INTEGER   NC,J,I
C
      PTRAPI = 1.0
      PTRAPH = 0.9
      TALPHA1= (1.,0.) 
      TALPHA2= 1.
      TALPHA3= 1.
      TALPHA4= (1.,0.) 
      TALPHA5= (1.,0.) 
      TALPHA6= 1.0
      TALPHA7= 0.
      TALPHA8= 0.9
      CTJXB  = 0.
      CTNTV  = 0.
      CTREY  = 0.
      CTERGO = 0.
      CTEDGE = 1.0
      CDMHD  = 0.
      CDNTV  = 0.
      CDERGO = 0.
      CDEDGE = 1.0
      ALPHACOL = 0.
      ALTAU  = CMPLX(0.5,0.)
      ATAU   = CMPLX(1.E+11,0.)
      EPSPAM = 1.E-7
      EPSDET = 1.E-199
      ETA    = 0.
      ETAXIS = 0.
      GAMMA  = 5./3.
      NUII   = 0.
      PVISC  = 0.
      ALPHAR = 2.
      BETAR  = 1.
      ALPHAP = 0.5
      ALPHAD = 1.0
      ALPHAH = 1.0
      FRACPTH= 1.0
      RTAN   = 0.5
      NUEFFIA  = 0.0
      NUEFFEA  = 0.0
      OMEGACI0 = 5.0e+1
      TDELTAMAX= 0.08
      TDELTAMIN= 0.04
      TDELTALIM= 100.
      TDELTALOW= 1.1 
      TCHIM0   = 1.0e-7
      TVPINCH0 = 0.    
      TCHID0   = 1.0e-7
      EQFAC    = 1.
      ETASTI   = 1.
      ETASTE   = 1.
      GAMWID = 1.E9
      GAMQ0  = 0.
      P0     = 0.
      P0OLD  = 0.
      ETACS1 = -1.
      ETACS2 = -1.
      ETACS3 = -1.
      ETACS4 = -1.
      ETAI   = 1.5
      RHO1   = 0.
      RHO2   = 0.
      SIGMA1 = 0.
      SIGMA2 = 0.
      DEXQ   =-2.
      DINERT = 1.
      TPOWVD = 0.
      TPOWDD = 0.
      RNTOR  = -1.
      TAUW   = 2.E+4
      ROTE   = 1.0e-2
      ROTEC  = 0.0e-2
      ROTEP  = 0.0e-2
      OMEGA0 = 0.0
      OMEGA1 = 0.0
      CSROT0 = 1.0
      ADAPSDF= 0.2
      ADAPSD2= 0.0
      ADAPSLF= 0.2
      ADAPSNU= 0.5   
      ZZETA0 = 0.8
      DZETA0 = SQRT(0.015)
      ROTWE0 = 0.0
      ZEFF0  = 1.0
      R0TYPE4= 1.0
      S0TYPE4= 0.0
      HHTYPE4C = 0.5
      DELRATS  = 0.
      TTCCONV0 = 1.
      TTCPARA0 = 0.  
      TTCPERP0 = 0.  
      TTCINERT0= 1.  
      B0K      = 1.0
      SLEFT    = 0.0
      JRE_SEED = 1.0E-5
      JRE_EQFRAC = 1.0
      JRE_DIFF = 0.0
      V1U_DIFF = 0.0
      JRE_EXB  = 1.0
      RE_RP_FAC= 1.0
      EP_DENF0 = 0.05
      EP_PREF0 = 0.5
      ORMIN  = 0.6
      ORMAX  = 1.4
      OZMIN  =-0.4
      OZMAX  = 0.4
      ORBE   = 100.0
      ORBL   = 0.25
      ORBSIGINI      = 1.0

      ESPECIES_Z     = 1.0
      ESPECIES_M     = 2.0
      PSPECIES_AP    = 0.0
      PSPECIES_AT    = 0.0
      PSPECIES_NP    = 0.0
      PSPECIES_NTB   = 0.0
      PSPECIES_NTD   = 0.0
      PSPECIES_NDB   = 0.0
      PSPECIES_FOWP  = 0.0
      PSPECIES_FOWT  = 0.0

      ESPECIES_Z(1)  = 1.0
      ESPECIES_Z(2)  =-1.0
      ESPECIES_M(1)  = 2.0
      ESPECIES_M(2)  = 5.4463e-04
 
      ODJPHI = 0.0
      OCHI0  = 0.0
C
      M1     = 2
      M2     = 2
      NSWEEP = 1
      OSWEEP =-1
      NNSCAN = 0
      NCONVB1= 0
      NCONVCS= 0
      NVACJ  = 0
      NCASE  = 1
      NFIT   = 3
      NTORQ  = 1
      NDNTR  = 1
      NITMAX = 50
      KGAM    = 10
      MOUTPUT = 80
      NPROFN = 0
      NPROFT = 0
      NPROFR = 0
      NPROFRC= 0
      NPROFRP= 0
      NPROFG = 0
      NPROFIE = 0
      NPROFWE = 0
      NPROFUI = 1
      NPROFUE = 1
      NUMODEL = 0
      NPROFK  = 0
      NPROFVD = 0
      NPROFVP = 0
      NPROFDD = 0
      NDENEQ  = 0
      NPROFTTCA = 0
      NPROFTTCE = 0
      NPROFNF = 1
      NPROFZ  = 0
      NRES   = 6
      NRESR  = 6
      NV     = 0
      KJRER  = 0
      KJRE_INIT = 1
      KVSQLIN= .FALSE.
      KVSQL  = .FALSE.
      KVSQL(8)=.TRUE.
      V2U_M0 = .FALSE.
      NWALL  = 0
      MWALL  = 0
      IWO    = 1
      INCKIN = 0
      INCDPHI= 0
      NCOUPL = 0
      IPERTURB  = 1
      KENORM    = 2
      KEFORM    = 2
      IVISC     = 0
      IDIAMV    = 0
      IDIAMB    = 0
      IDIAMTI   = 0
      IDIAMTE   = 0
      FDIAMV    = 1.0
      FDIAMB    = 1.0
      FDIAMTI   = 1.0
      IPDIVB    = 0
      K_BC_N0   = 11
      IGAMMADIVV= 0
      KSMTYPE4  = 1
      JSOUT     = 10
      IOMPNUM   = 1
      INUTYPE   = 0
      NORR      = 0
      NOZZ      = 0
      NOPP      = 0
      IORBIT    = 0
      NREORBIT  = 0
      KRE1      = 1
      KRE2      = 2
      NRE1      = 1
      NRE2      = 1
      KRE_INIT  = 3
      KRE_STEP  =-1
      KRE_STEP_MAX = 10000000
      KRE_TRACE = 10
      KRE_VAC    = 0
      RE_BVAC    = 1.0
      KRE_PERTURB_TIME = 0
      KRE_NMAX  = 0
      KRE_STAR  = 5
      KRE_FO    = 0
      KRE_ODE   = 1
      IFOW      = 0
      NORB      = 1000
      NSPECIES  = 2
      KEPSALPHA = 0
      IFOWPSI0  = 1
      KFASTRUN  = 1
      NPROFR4   = 0
      NPROFS4   = 0
      ISPECIES_F0(1) = 0
      ISPECIES_F0(2) = 0
      IT_BC     = 2
      ID_BC     = 0
      ITSATURAT = 0
      ZCHARGE   = 1
      KRATSURF  = 1
      MRATSURF  = 0
      NADAPS    = 0
      ADAPSOPT  = 1
      KADAPS    = 0
      KNBI      = 13
      NEXPV     = 0
      V2XKEY    = 0
      NLAMK     = 201
      NCHI0     = 120
      NLAMIN    = 100
      NEPK      = 101
      NKL0      = 0
      NKSMOOTHB = 20
      NKSMOOTHR = 10
      NKSINGULAR= 1
      KSOLREAD  = 0
      KSOLSAVE  = 0
      KSOLTEST = 0
      IWALL     = 0
      KNTV      = 10
      KISLAND   = 1
      NFLAGROTA = 0
      INERT1    = .TRUE.
      INERT2    = .TRUE.
      INERT3    = .TRUE.
      ODWKCOM   = .FALSE.
      ISLSODE   = .TRUE.
C
      DCONTI = .FALSE.
      DTERMS = .FALSE.
      OUTDAT = .FALSE.
C
      DPLOLP = .FALSE.
      DPRINT = .FALSE.
C
      DPLOLP( 1) = .TRUE.
      DPLOLP( 4) = .TRUE.

      NPICK = 0
      IPICK = 2
      CPICK = 0.0

      KEYTORQ = 0

      RE_PERTURB = (1.0,0.0)
      RE_EFIELD  = 0.0
      RE_SAC     = 0.0
      RE_SYNCH   = 0.0
      RE_BREMS   = 0.0
      RE_TIME    = 10000.0
      RE_TSTEP   = 1.0E-3
      RE_CONST(1)= 4.0E+20
      RE_CONST(2)= 5.0   
      RE_CONST(3)= 5.0   
      RE_CONST(4)= 1.0    
      RE_CONST(5)= 1.8
      RE_CONST(6)= 2.5
      RE_CONST(7)= 12.5
      RE_PMAX    = 10.0
      RE_SP      = 2.0
      RE_S0      = 0.5
      RE_S0MIN   = 0.0
      RE_S0MAX   = 1.0
      RE_CHI0    = 0.0
      RE_PHI0    = 0.0
      RE_P0      = 1.0E-02 
      RE_LAMBDA0 = 0.0
      KRE_SIGMA0 = 1
      RE_ATOL    = 1.0E-7
      RE_AN      = 5.446278570617179e-04
      RE_CN      =-1.0
      RE_PERTURB_MAX = 1.0
      RE_PEDGE   = 0.97
      RE_EIGVAL  = (1.0,0.0)
      NPHIMAX    = 50
      
C
      RETURN
      END
*DECK RDNAME
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------READ NAMELIST--------------A.B.   31/01/90 --------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE RDNAME
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      USE TORQUEM
      USE REORBITM 
      USE MPIENV
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INCLUDE 'newrun.inc'
      INTEGER    INCKXJ2L,INCKXX1,INCKXPD,INCKXPED,INCKXW1,
     &           INCKXJRE,INCKXB2L,INCKXB3L,INCKXJRE2,INCKXJRE3,
     &           INCKXJRE2L,INCKXDPHI
      INTEGER    INCKYX2,INCKYX3,INCKYPE,INCKYPP,INCKYPPARA,INCKYPPERP,
     &           INCKYRHOP,INCKYW2,INCKYW3,INCKYV3,INCKYJRE1
      
      IGO =0

C     RESET EVERY TIME

      KEEPMO = .FALSE.
      LPROFI = .FALSE.

      READ(CHNAME,BASIC,END=100, ERR=210)
      READ(CHNAME,FEEDBACK,      ERR=220)
      READ(CHNAME,KINETIC,       ERR=230)
      READ(CHNAME,QLIN,          ERR=240)
      READ(CHNAME,NUMERIC,       ERR=250)
      READ(CHNAME,OUTOPT,        ERR=260)
      
      CALL CHECK_MPI_RUN()
      CALL CREATE_WORKDIR()

      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN     
        WRITE(*,1000)
        WRITE(*,BASIC)
        WRITE(*,1000)
        WRITE(*,FEEDBACK)
        WRITE(*,1000)
        WRITE(*,KINETIC)
        WRITE(*,1000)
        WRITE(*,QLIN)
        WRITE(*,1000)
        WRITE(*,NUMERIC)
        WRITE(*,1000)
        WRITE(*,OUTOPT)
      ENDIF

C     DECIDE WHICH PERTURBED VARIABLES TO BE INCLUDED INTO EQUATIONS
      INCKXJ2L   = 0
      INCKXX1    = 0
      INCKXPD    = 0
      INCKXPED   = 0
      INCKXW1    = 0
      INCKXJRE   = 0
      INCKXB2L   = 0
      INCKXB3L   = 0
      INCKXJRE2  = 0
      INCKXJRE3  = 0
      INCKXJRE2L = 0
      INCKXDPHI  = 0
  
      IF (ABS(ETA).GT.0..OR.IWO.EQ.2)                    INCKXJ2L   = 1
      IF ((ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.ABS(ROTE).GT.0.).OR.
     &    V2XKEY.GT.0.OR.(INCKIN.GT.0.AND.IPERTURB.GT.0).OR.
     &    NCASE.EQ.6.OR.INCKIN.EQ.2.OR.
     &    (ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.NPROFRC.GT.0).OR.
     &    (ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.NPROFRP.GT.0))
     &                                                   INCKXX1    = 1
      IF (IDIAMB.EQ.1.OR.IDIAMB.EQ.3)                    INCKXPD    = 1
      IF (IDIAMTI.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3).AND.
     &    INCKIN.GT.0)                                   INCKXPED   = 1
      IF (ABS(TTCPARA0)+ABS(TTCPERP0).GT.0..AND.
     &    (V2XKEY.EQ.0.OR.V2XKEY.EQ.1))                  INCKXW1    = 1
      IF (KJRER.GT.0.AND.KJRER.LE.6)                     INCKXJRE   = 1
      IF (KJRER.EQ.5.AND.JRE_EQFRAC.GT.0.)               INCKXB2L   = 1
      IF (KJRER.EQ.6.AND.JRE_EQFRAC.GT.0.)               INCKXB2L   = 1
      IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10)          INCKXB2L   = 1
      IF (KJRER.EQ.5.AND.JRE_EQFRAC.GT.0.)               INCKXB3L   = 1
      IF (KJRER.EQ.6.AND.JRE_EQFRAC.GT.0.)               INCKXB3L   = 1
      IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10)          INCKXB3L   = 1
      IF (KJRER.EQ.6)                                    INCKXJRE2  = 1
      IF (KJRER.EQ.6)                                    INCKXJRE3  = 1
      IF (KJRER.EQ.6.AND.ABS(ETA).GT.0.)                 INCKXJRE2L = 1
      IF (INCKIN.GT.0.AND.INCDPHI.GT.0)                  INCKXDPHI  = 1

      INCKYV3    = 0
      INCKYX2    = 0
      INCKYX3    = 0
      INCKYPE    = 0
      INCKYPP    = 0
      INCKYPPARA = 0
      INCKYPPERP = 0
      INCKYRHOP  = 0
      INCKYW2    = 0
      INCKYW3    = 0
      INCKYJRE1  = 0
  
      IF (ABS(ROTE).GT.0.OR.ABS(GAMMA).GT.0.)             INCKYV3    = 1
      IF ((ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.ABS(ROTE).GT.0.).OR.
     &    V2XKEY.GT.0.OR.(INCKIN.GT.0.AND.IPERTURB.GT.0).OR.
     &    NCASE.EQ.6.OR.INCKIN.EQ.2.OR.
     &    (ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.NPROFRC.GT.0).OR.
     &    (ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.NPROFRP.GT.0))
     &                                                    INCKYX2    = 1
      IF (V2XKEY.GT.0.AND.INCKYV3.EQ.1.OR.NCASE.EQ.6.OR.
     &    (ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.NPROFRC.GT.0).OR.
     &    (ABS(PVISC).GT.0..AND.IVISC.EQ.1.AND.NPROFRP.GT.0))
     &                                                    INCKYX3    = 1
      IF (IDIAMTI.NE.0.OR.IDIAMTE.NE.0)                   INCKYPE    = 1
      IF (IDIAMB.EQ.1.OR.IDIAMB.EQ.3.OR.IDIAMB.EQ.4)      INCKYPE    = 1
      IF (IDIAMTE.EQ.0.AND.IDIAMB.EQ.1.AND.INCKIN.GT.0)   INCKYPE    = 1
      IF (IDIAMTE.EQ.0.AND.IDIAMB.EQ.3.AND.INCKIN.GT.0)   INCKYPE    = 1
      IF (IDIAMTE.EQ.0.AND.IDIAMB.EQ.1.AND.INCKIN.GT.0)   INCKYPP    = 1
      IF (IDIAMTE.EQ.0.AND.IDIAMB.EQ.3.AND.INCKIN.GT.0)   INCKYPP    = 1
      IF (INCKIN.GT.0)                                    INCKYPPARA = 1
      IF (INCKIN.GT.0)                                    INCKYPPERP = 1
      IF ((INERT3.AND.ABS(ROTE).GT.0.).OR.ABS(CDMHD).GT.0.) INCKYRHOP= 1
      IF (KVSQLIN)                                          INCKYRHOP= 1
      IF (ABS(TTCPARA0)+ABS(TTCPERP0).GT.0..AND.
     &    (V2XKEY.EQ.0.OR.V2XKEY.EQ.1))                   INCKYW2    = 1
      IF (ABS(TTCPARA0)+ABS(TTCPERP0).GT.0..AND.
     &    (V2XKEY.EQ.0.OR.V2XKEY.EQ.1))                   INCKYW3    = 1
      IF (KJRER.EQ.6)                                     INCKYJRE1  = 1

C     DEFINE NXCOMP,NYCOMP AND NUMBERING OF VARIABLES
      KXV1    = 1
      KXB1    = 2
      KXJ2U   = 3
      KXJ3    = 4
      KXJ2L   = -1
      KXX1    = -1
      KXPD    = -1
      KXPED   = -1
      KXW1    = -1
      KXJRE   = -1
      KXB2L   = -1
      KXB3L   = -1
      KXJRE2  = -1
      KXJRE3  = -1
      KXJRE2L = -1
      KXDPHI  = -1

      NXCOMP  = 4
      IF (INCKXJ2L.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXJ2L  = NXCOMP
      ENDIF
      IF (INCKXX1.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXX1   = NXCOMP
      ENDIF
      IF (INCKXPD.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXPD   = NXCOMP
      ENDIF
      IF (INCKXPED.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXPED  = NXCOMP
      ENDIF
      IF (INCKXW1.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXW1   = NXCOMP
      ENDIF
      IF (INCKXJRE.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXJRE  = NXCOMP
      ENDIF
      IF (INCKXB2L.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXB2L  = NXCOMP
      ENDIF
      IF (INCKXB3L.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXB3L  = NXCOMP
      ENDIF
      IF (INCKXJRE2.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXJRE2 = NXCOMP
      ENDIF
      IF (INCKXJRE3.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXJRE3 = NXCOMP
      ENDIF
      IF (INCKXJRE2L.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXJRE2L = NXCOMP
      ENDIF
      IF (INCKXDPHI.EQ.1) THEN
         NXCOMP = NXCOMP + 1
         KXDPHI = NXCOMP
      ENDIF

      KYV2    = 1
      KYB2    = 2
      KYB3    = 3
      KYJ1    = 4
      KYPR    = 5
      KYV3    = -1
      KYX2    = -1
      KYX3    = -1
      KYPE    = -1
      KYPP    = -1
      KYPPARA = -1
      KYPPERP = -1
      KYRHOP  = -1
      KYW2    = -1
      KYW3    = -1
      KYJRE1  = -1
 
      NYCOMP  = 5       
      IF (INCKYV3.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYV3   = NYCOMP
      ENDIF
      IF (INCKYX2.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYX2   = NYCOMP
      ENDIF
      IF (INCKYX3.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYX3   = NYCOMP
      ENDIF
      IF (INCKYPE.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYPE   = NYCOMP
      ENDIF
      IF (INCKYPP.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYPP   = NYCOMP
      ENDIF
      IF (INCKYPPARA.EQ.1) THEN
         NYCOMP  = NYCOMP + 1
         KYPPARA = NYCOMP
      ENDIF
      IF (INCKYPPERP.EQ.1) THEN
         NYCOMP  = NYCOMP + 1
         KYPPERP = NYCOMP
      ENDIF
      IF (INCKYRHOP.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYRHOP = NYCOMP
      ENDIF
      IF (INCKYW2.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYW2   = NYCOMP
      ENDIF
      IF (INCKYW3.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYW3   = NYCOMP
      ENDIF
      IF (INCKYJRE1.EQ.1) THEN
         NYCOMP = NYCOMP + 1
         KYJRE1 = NYCOMP
      ENDIF
  
      IF (ISMPIRUN.EQ.0.OR.RANK.EQ.ROOT) THEN
      WRITE(*,2000) NXCOMP,KXV1,KXB1,KXJ2U,KXJ3,KXJ2L,KXX1,KXPD,KXPED,
     &              KXW1,KXJRE,KXB2L,KXB3L,KXJRE2,KXJRE3,KXJRE2L,KXDPHI
      WRITE(*,2100) NYCOMP,KYV2,KYV3,KYB2,KYB3,KYJ1,KYPR,KYX2,KYX3,   
     &              KYPE,KYPP,KYPPARA,KYPPERP,KYRHOP,KYW2,KYW3,KYJRE1
 2000 FORMAT('NXCOMP=',I2,' KXV1=',I2,' KXB1=',I2,' KXJ2U=',I2,
     &       ' KXJ3=',I2,' KXJ2L=',I2,' KXX1=',I2,' KXPD=',I2,
     &       ' KXPED=',I2,
     &       ' KXW1=',I2,' KXJRE=',I2,' KXB2L=',I2,' KXB3L=',I2,
     &       ' KXJRE2=',I2,' KXJRE3=',I2,' KXJRE2L=',I2,' KXDPHI=',I2)
 2100 FORMAT('NYCOMP=',I2,' KYV2=',I2,' KYV3=',I2,' KYB2=',I2,
     &       ' KYB3=',I2,' KYJ1=',I2,' KYPR=',I2,' KYX2=',I2,
     &       ' KYX3=',I2,' KYPE=',I2,' KYPP=',I2,
     &       ' KYPPARA=',I2,' KYPPERP=',I2,
     &       ' KYRHOP=',I2,' KYW2=',I2,' KYW3=',I2,
     &       ' KYJRE1=',I2)
      ENDIF

C     DERIVED QUANTITIES
      CALPHA1 = TALPHA1
      CALPHA2 = TALPHA2
      CALPHA3 = TALPHA3
      CALPHA4 = TALPHA4
      CALPHA5 = TALPHA5 
      CALPHA6 = TALPHA6
      CALPHA7 = TALPHA7
      CALPHA8 = TALPHA8

      KXJFB   = KXJ2U
      IF (ABS(RNTOR).LT.1.0E-10) THEN
         KXJFB  = KXJ3
      ENDIF

      MSMAX = M2 - M1 + 1
      MSDIM = MSMAX

      MXMAX = MSMAX * NXCOMP
      MYMAX = MSMAX * NYCOMP
      MVMAX = MSMAX * 2
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN     
      WRITE(*,'("MXMAX,MYMAX= ",2I5)')MXMAX,MYMAX
      ENDIF
      
      IGO = 1
      RETURN

 100  IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN     
      WRITE(*,*) ' END AT NAMELIST READ'
      ENDIF
      RETURN

 210  WRITE(*,*) ' ERROR AT NAMELIST READ: BASIC'
      WRITE(*,BASIC)
      RETURN
 220  WRITE(*,*) ' ERROR AT NAMELIST READ: FEEDBACK'
      WRITE(*,FEEDBACK)
      RETURN
 230  WRITE(*,*) ' ERROR AT NAMELIST READ: KINETIC'
      WRITE(*,KINETIC)
      RETURN
 240  WRITE(*,*) ' ERROR AT NAMELIST READ: QLIN'
      WRITE(*,QLIN)
      RETURN
 250  WRITE(*,*) ' ERROR AT NAMELIST READ: NUMERIC'
      WRITE(*,NUMERIC)
      RETURN
 260  WRITE(*,*) ' ERROR AT NAMELIST READ: OUTOPT'
      WRITE(*,OUTOPT)
      RETURN

 1000 FORMAT(//,80('*'),//)
      END

*DECK CLEAR
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------CLEAR VARIABLES------------G.VLAD 05/04/1989-------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....INITIALIZE SOME VARIABLES.......................................
C
      SUBROUTINE CLEAR
C     ================
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE GVACUUMM
      USE MPIENV
      INTEGER          I,MS
C
      NTP1=NRP1+NV
      NTOT =NTP1-1
C
C     ALLOCATE VARIABLES
C
      IF (.NOT. ALLOCATED(V1U)) THEN
         ALLOCATE(V1U(NTP1,MSDIM),V2U(NTP1,MSDIM),V3U(NTP1,MSDIM))
         ALLOCATE(X1U(NRP1,MSDIM),X2U(NRP1,MSDIM),PDE(NRP1,MSDIM))
         ALLOCATE(B1U(NTP1,MSDIM),B2U(NTP1,MSDIM),B3U(NTP1,MSDIM))
         ALLOCATE(J1U(NTP1,MSDIM),J2U(NTP1,MSDIM),J3U(NTP1,MSDIM))
         ALLOCATE(PRE(NRP1,MSDIM),PEE(NRP1,MSDIM),PEP(NRP1,MSDIM))
         ALLOCATE(PED(NRP1,MSDIM),J2L(NTP1,MSDIM))
         ALLOCATE(PPERP(NRP1,MSDIM),PPARA(NRP1,MSDIM))
         ALLOCATE(X3U(NRP1,MSDIM),RHOP(NRP1,MSDIM))
         ALLOCATE(DPHI(NRP1,MSDIM))
      ENDIF
      IF (.NOT. ALLOCATED(JRE).AND.KJRER.GT.0) THEN
         ALLOCATE(JRE(NRP1,MSDIM))
         ALLOCATE(JRE1(NRP1,MSDIM),JRE2(NRP1,MSDIM),JRE3(NRP1,MSDIM))
         ALLOCATE(JRE2L(NRP1,MSDIM))
      ENDIF
 
      V1U=0.
      V2U=0.
      V3U=0.
      X1U=0.
      X2U=0.
      X3U=0.
      B1U=0.
      B2U=0.
      B3U=0.
      J1U=0.
      J2U=0.
      J3U=0.
      J2L=0.
      PRE=0.
      PEE=0.
      PEP=0.
      PDE=0.
      PED=0.
      PPERP=0.
      PPARA=0.
      RHOP=0.
      DPHI=0.
      IF (KJRER.GT.0) THEN
         JRE  =0.
         JRE1 =0.
         JRE2 =0.
         JRE3 =0.
         JRE2L=0.
      ENDIF
C
C
      IF (.NOT. ALLOCATED(RPF)) THEN
      ALLOCATE( RPF(nrp1,MEDIM),    ZPF(nrp1,MEDIM),
     $         RPFM(nrp1,MEDIM),   ZPFM(nrp1,MEDIM),
     $        DPEDS(nrp1,MEDIM), DPEDSM(nrp1,MEDIM))
      ALLOCATE(DG11L(nrp1,MEDIM), DG22L(nrp1,MEDIM),  DG33L(nrp1,MEDIM),
     $       DG12L(nrp1,MEDIM),   GCHDZ(nrp1,MEDIM),
     $       GSDZ(nrp1,MEDIM)    ,GBZ(nrp1,MEDIM)  ,    GBR(nrp1,MEDIM))
      ALLOCATE(DG11LM(nrp1,MEDIM),DG22LM(nrp1,MEDIM),DG33LM(nrp1,MEDIM),
     $       DG12LM(nrp1,MEDIM),  GCHDZM(nrp1,MEDIM),
     $       GSDZM(nrp1,MEDIM)   ,GBZM(nrp1,MEDIM)  ,  GBRM(nrp1,MEDIM))
      ALLOCATE(                   JACOBI(nrp1,MEDIM),JACOBM(nrp1,MEDIM))
      END IF
C
C
      DG11L  = 0.
      DG22L  = 0.
      DG12L  = 0.
      DG33L  = 0.
C
      DG11LM = 0.
      DG22LM = 0.
      DG12LM = 0.
      DG33LM = 0.
C
      JACOBI = 0.
      JACOBM = 0.
C
      RETURN
      END
*DECK INITPE
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------INITILIZE PERTURBATION------- A.B. 31/01/90 -------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....INITIALIZE SOME VARIABLES.......................................
C
      SUBROUTINE INITPE
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      USE RCOMDM
      INCLUDE 'compam.inc'

      INTEGER      I,MS
      COMPLEX*16   ZAMP
      REAL*8       CUM, RAMP, TMP
      PARAMETER    (CUM=3.5)
C
      RAMP = 1.0E-6
      IF (NCASE.EQ.9.OR.NCASE.EQ.10) RAMP = THRESHOLD

      IF (INCFEED.EQ.8.OR.INCFEED.EQ.4.OR.INCFEED.EQ.10.OR.
     &    INCFEED.EQ.18.OR.INCFEED.EQ.20.OR.INCFEED.EQ.21.OR.
     &    INCFEED.EQ.22) RAMP = 0.0

      IF (INCFEED.EQ.8.AND.NCASE.EQ.5) RAMP = THRESHOLD

      DO 10 MS=1,MSMAX
CPPPL         ZAMP = (CUM**2 - RM(MS,2)) * CEXP(CMPLX(0.,100.*MS*MS))
      ZAMP = (CUM**2 - RM(MS,2)**2)
      IF (ABS(RM(MS,2)).GT.CUM) ZAMP = 0.01
      IF (ABS(RM(MS,2)).LT.0.1) ZAMP = 0.0

CLIU  DO  10 I=1,NRP1
CLIU      B1U(I,MS) = ZAMP*DFLOAT((I-1)*(NRP1-I))/DFLOAT(NR*NR)
CLIU  10  CONTINUE
      DO  I=1,NTP1
          B1U(I,MS) = ZAMP*DFLOAT((I-1)*(NTP1-I))/DFLOAT(NTOT*NTOT)
      ENDDO
 10   CONTINUE

C     --- YQLIU, 17/08/2001
C     --- RESCALE INITIAL PERTURBATION TO HAVE THE REAL*8 AMPLITUDE
C     --- EQUAL TO <RAMP>, IMPORTANT FOR TIME EVOLUTION (NCASE=3,4)
      TMP = 0.0
      DO MS=1,MSMAX
      DO  I=1,NTP1
         IF (ABS(B1U(I,MS)).GT.TMP) TMP = ABS(B1U(I,MS))
      ENDDO
      ENDDO

      DO MS=1,MSMAX
      DO I=1,NTP1
         B1U(I,MS) = B1U(I,MS)*RAMP/TMP
      ENDDO
      ENDDO

C     SET UP SEED RE CURRENT FOR M=0 HARMONIC
      IF (KXJRE.GT.0.AND.KJRER.GE.1.AND.KJRER.LE.4) THEN
      DO MS=1,MSMAX
      IF (ABS(RM(MS,2)).LT.0.1) THEN 
         DO I=1,NRP1
            IF (KJRE_INIT.EQ.1) JRE(I,MS) = JRE_SEED
            IF (KJRE_INIT.EQ.2) JRE(I,MS) = JRE_SEED*(1.-CS(I)**2)**4
            IF (KJRE_INIT.EQ.3) JRE(I,MS) = JRE_SEED*
     &         CS(I)**4*(1.-CS(I)**2)/0.1481422009686093
         ENDDO
      ENDIF
      ENDDO
      ENDIF

      RETURN
      END
*DECK GEOMET
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------GEOMETRIC FACTORS --------G. VLAD ------------------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE GEOMET
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE MPIENV
      INCLUDE 'comioc.inc'
      INTEGER          I,J,K
      REAL*8  ZRAD
      REAL*8,DIMENSION(:),ALLOCATABLE::ZVOL
      ALLOCATE(ZVOL(NRP1))
C
C.....EXTENSIONS FOR FREE-BOUNDARY
C
      CSM(NRP1) = CS(NRP1)
      IF (NV.GE.2) CSM(NRP1) = VCSM(1)
      CS(NRP1+1) = 2.*CS(NRP1) - CS(NR)
      IF (NV.GE.2) CS(NRP1+1) = VCS(2)
C
      DO 10 I=1,NR
      CSH(I)=CS(I+1)-CS(I)
 10   CONTINUE
      CSH(0) = 0.
      CSH(NRP1) = 0.
C      IF (NV.GE.2) CSH(NRP1) = VCS(2) - VCS(1)
C
C     THIS IS TO ACHIEVE EQUILIBRIUM J*BEQ = 0 ON AXIS TO ADVANCE
C     B2U AND B3U
C
      IF (NV.GE.2) CALL VACGEO
C
C     CONSTRUCT GRIDS CSV AND CSVM CORRESPONDING TO SQUARE ROOT
C     OF ENCLOSED VOLUME ON CS AND CSM MESHES
C
      ZVOL(1) = 0.
      DO 60  I=1,NR
 60   ZVOL(I+1) = ZVOL(I) + DREAL(JACOBM(I,1))*(CS(I+1)-CS(I))

      VOLTOT = ZVOL(NRP1)

      DO  70 I=1,NRP1
      ZRAD=ZVOL(I)/ZVOL(NRP1)
 70   CSV(I)=SQRT(ZRAD)
      DO 80 I=1,NR
 80   CSVM(I) = (CSV(I)+CSV(I+1))*0.5
C
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*) 'LOCATION OF RATIONAL SURFACES: NRATSURF=',NRATSURF
      WRITE(*,*) 'IRATSURF CS Q  DQDS   DSDR'
      DO J=1,NRATSURF
         K = IRATSURF(J)
         WRITE(*,120) K,CS(K),QPLS(K),(QPLS(K+1)-QPLS(K))/CSH(K),
     &                CSH(K)/(CSV(K+1)-CSV(K))
      ENDDO
      ENDIF
 120  FORMAT(I4,4(1X,E13.5))

      DEALLOCATE(ZVOL)
      RETURN
      END
*DECK LISTMO
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------PRINT PAMS VECTOR--------- A.B.  02/02/90 --------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
      SUBROUTINE LISTMO(X,Y)
C     ======================
C
      USE DIMENSIM
      USE GLOBALM
      COMPLEX*16    X(NXCOMP,MSMAX,*),Y(NYCOMP,MSMAX,*)
      INTEGER       MS,NC,MNX,MNY
C
      DO 10 NC = 1,NXCOMP+NYCOMP
 10   IF (DPRINT(NC)) GOTO 20
      RETURN
C
 20   MNX = MSMAX*NXCOMP
      MNY = MSMAX*NYCOMP
C
      DO 100 MS=1,MSMAX
      WRITE(*,1000) RM(MS,2)
      IF (DPRINT( 1)) CALL PRIMAT('V1    ',X(KXV1,MS,1),1,NRP1,MNX)
      IF (DPRINT( 2)) CALL PRIMAT('V2    ',Y(KYV2,MS,1),1,NR  ,MNY)
      IF (DPRINT( 3).AND.KYV3.GT.0) 
     &                CALL PRIMAT('V3    ',Y(KYV3,MS,1),1,NR  ,MNY)

      IF (DPRINT( 4)) CALL PRIMAT('B1    ',X(KXB1,MS,1),1,NRP1,MNX)
      IF (DPRINT( 5)) CALL PRIMAT('B2    ',Y(KYB2,MS,1),1,NR  ,MNY)
      IF (DPRINT( 6)) CALL PRIMAT('B3    ',Y(KYB3,MS,1),1,NR  ,MNY)

      IF (DPRINT( 7)) CALL PRIMAT('J1    ',Y(KYJ1,MS,1),1,NR  ,MNY)
      IF (DPRINT( 8)) CALL PRIMAT('J2U   ',X(KXJ2U,MS,1),1,NRP1,MNX)
      IF (DPRINT( 9)) CALL PRIMAT('J3    ',X(KXJ3,MS,1),1,NRP1,MNX)
      IF (DPRINT(10).AND.KXJ2L.GT.0) 
     &                CALL PRIMAT('J2L   ',X(KXJ2L,MS,1),1,NRP1,MNX)
      IF (DPRINT(11).AND.KXPD.GT.0) 
     &                CALL PRIMAT('PDE   ',X(KXPD,MS,1),1,NRP1,MNX)
      IF (DPRINT(12)) CALL PRIMAT('PRE   ',Y(KYPR,MS,1),1,NR  ,MNY)
C
      IF (DPRINT(13).AND.KXX1.GT.0) 
     &                CALL PRIMAT('X1    ',X(KXX1,MS,1),1,NRP1,MNX)
      IF (DPRINT(14).AND.KYX2.GT.0) 
     &                CALL PRIMAT('X2    ',Y(KYX2,MS,1),1,NR  ,MNY)
      IF (DPRINT(15).AND.KYPPERP.GT.0) 
     &                CALL PRIMAT('PPERP ',Y(KYPPERP,MS,1),1,NR,MNY)
      IF (DPRINT(16).AND.KYPPARA.GT.0) 
     &                CALL PRIMAT('PPARA ',Y(KYPPARA,MS,1),1,NR,MNY)
      IF (DPRINT(17).AND.KYPE.GT.0) 
     &                CALL PRIMAT('PPE   ',Y(KYPE,MS,1),1,NR,MNY)
      IF (DPRINT(18).AND.KYPP.GT.0) 
     &                CALL PRIMAT('PPE   ',Y(KYPP,MS,1),1,NR,MNY)
      IF (DPRINT(19).AND.KYRHOP.GT.0) 
     &                CALL PRIMAT('RHOP  ',Y(KYRHOP,MS,1),1,NR,MNY)
      IF (DPRINT(20).AND.KYX3.GT.0) 
     &                CALL PRIMAT('X3    ',Y(KYX3,MS,1),1,NR,MNY)
 100  CONTINUE
      RETURN
C
 1000 FORMAT(//,' M = ',F5.1,'  COMPONENTS',/,1X,21('='),/)
      END
*DECK TERMS
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C------ TERMS IN FARADAY ---------- A.B.   04/02/90 --------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.. DIAGNOSTICS (FOR RESISTIVE BALLOONING MODES)
C.. SEPARATE IDEAL AND RESISTIVE TERMS IN DB/DT
C
      SUBROUTINE TERMS(
C     =================
     &     DX, DY, X, Y, WORK, IWORK,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'compam.inc'
      INCLUDE 'specmat.inc'
      INTEGER IWORK
      REAL*8    WORK(IWORK)
      COMPLEX*16 DX(NXCOMP,MSMAX),DY(NYCOMP,MSMAX),
     C         X(NXCOMP,MSMAX), Y(NYCOMP,MSMAX)
C
C.....LOCAL VARIABLES
      COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE::XC,YC
      INTEGER  I,MS,MNX
C
      ALLOCATE( XC(NXCOMP,MSDIM,NTP1),YC(NYCOMP,MSDIM,NTP1))
      MNX = MSMAX*NXCOMP
      CALL GETXY(X,Y)
C
      CALL LINEAR(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
C
C     PUT B TO ZERO TO AVOID SUBTRACTING OFF FULL DB/DT BY SHIFT OF A
C
      DO  1 MS = 1,MSMAX
      DO  1 I = 1,NTP1
      B1U(I,MS) = 0.
      B2U(I,MS) = 0.
 1    B3U(I,MS) = 0.
C
C     COPY BACK INTO TEMPORARY (X,Y) VECTOR (XC,YC)
C
      CALL INITXY(XC,YC)
C
C     TOTAL DB/DT FROM TEMPORARY VECTOR (XC,YC)
C
C     CALL CALPAM(
C    $      MXMAX,MYMAX,NTOT,NXCOMP,NYCOMP, -3  ,NITMAX
C    $     ,EPSPAM,EPSDET
C    $     ,AL0,ALAM,ALNORM,NONCON
C    $     ,ASUBM(1+MXMAX*MXMAX)
C    $     ,BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM
C    $     ,HSUBM(1+MXMAX*MYMAX)
C    $     ,DX,DY, XC, YC
C    $     ,WORK   ,IWORK
C    $                                        )
C
      CALL CALPAM(
     $      MXMAX,MYMAX,NTOT,NXCOMP,NYCOMP, -3  ,NITMAX
     $     ,EPSPAM,EPSDET
     $     ,AL0,ALAM,ALNORM,NONCON
     $     ,ASUBM(1,1,2)
     $     ,BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM
     $     ,HSUBM(1,1,2)
     $     ,DX,DY, XC, YC
     $     ,WORK   ,IWORK
     $                                        )
C
      WRITE(*,1010)
                     CALL PLOTV(DX(KXB1,1),Q ,'B1DTOT',NXCOMP,NRP1)
      IF (DPLOLP(5)) CALL PLOTV(DY(KYB2,1),QM,'B2DTOT',NYCOMP,NR  )
      IF (DPLOLP(6)) CALL PLOTV(DY(KYB3,1),QM,'B3DTOT',NYCOMP,NR  )
C
C     IDEAL PART -- SUBTRACT OFF CURRENTS FROM (XC,YC)
C
      CALL GETXY(X,Y)
      DO 10 MS = 1,MSMAX
      DO 10 I = 1,NTP1
      B1U(I,MS) = 0.
      B2U(I,MS) = 0.
      B3U(I,MS) = 0.
      J1U(I,MS) = 0.
      J2U(I,MS) = 0.
      J3U(I,MS) = 0.
      J2L(I,MS) = 0.
 10   CONTINUE
C
      CALL INITXY(XC,YC)
C
C     CALL CALPAM(
C    $      MXMAX,MYMAX,NTOT,NXCOMP,NYCOMP, -3  ,NITMAX
C    $     ,EPSPAM,EPSDET
C    $     ,AL0,ALAM,ALNORM,NONCON
C    $     ,ASUBM(1+MXMAX*MXMAX)
C    $     ,BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM
C    $     ,HSUBM(1+MXMAX*MYMAX)
C    $     ,DX,DY, XC, YC
C    $     ,WORK   ,IWORK
C    $                                        )
C
      CALL CALPAM(
     $      MXMAX,MYMAX,NTOT,NXCOMP,NYCOMP, -3  ,NITMAX
     $     ,EPSPAM,EPSDET
     $     ,AL0,ALAM,ALNORM,NONCON
     $     ,ASUBM(1,1,2)
     $     ,BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM
     $     ,HSUBM(1,1,2)
     $     ,DX,DY, XC, YC
     $     ,WORK   ,IWORK
     $                                        )
C
      WRITE(*,1010)
                     CALL PLOTV(DX(KXB1,1),Q ,'B1DIDE',NXCOMP,NRP1)
      IF (DPLOLP(5)) CALL PLOTV(DY(KYB2,1),QM,'B2DIDE',NYCOMP,NR  )
      IF (DPLOLP(6)) CALL PLOTV(DY(KYB3,1),QM,'B3DIDE',NYCOMP,NR  )
C
C     NOW DO RESISTIVE PIECE: REFRESH PERTURBATION VECTOR AND ZERO V
C
      CALL GETXY(X,Y)
      DO MS = 1,MSMAX
         B1U(1:NTP1,MS)  = 0.
         B2U(1:NTP1,MS)  = 0.
         B3U(1:NTP1,MS)  = 0.
         X1U(1:NRP1,MS)  = 0.
         X2U(1:NRP1,MS)  = 0.
         X3U(1:NRP1,MS)  = 0.
         V1U(1:NTP1,MS)  = 0.
         V2U(1:NTP1,MS)  = 0.
         V3U(1:NTP1,MS)  = 0.

         IF (KJRER.GT.0) THEN
         JRE(1:NRP1,MS)  = 0.
         JRE1(1:NRP1,MS) = 0.
         JRE2(1:NRP1,MS) = 0.
         JRE3(1:NRP1,MS) = 0.
         JRE2L(1:NRP1,MS)= 0.
         ENDIF

         IF (KXDPHI.GT.0) DPHI(1:NRP1,MS)=0.
      ENDDO
C
      CALL INITXY(XC,YC)
C
      CALL CALPAM(
     $      MXMAX,MYMAX,NTOT,NXCOMP,NYCOMP, -3  ,NITMAX
     $     ,EPSPAM,EPSDET
     $     ,AL0,ALAM,ALNORM,NONCON
     $     ,ASUBM(1,1,2)
     $     ,BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM
     $     ,HSUBM(1,1,2)
     $     ,DX,DY, XC, YC
     $     ,WORK   ,IWORK
     $                                        )
C
      WRITE(*,1010)
                     CALL PLOTV(DX(KXB1,1),Q ,'B1DRES',NXCOMP,NRP1)
      IF (DPLOLP(5)) CALL PLOTV(DY(KYB2,1),QM,'B2DRES',NYCOMP,NR  )
      IF (DPLOLP(6)) CALL PLOTV(DY(KYB3,1),QM,'B3DRES',NYCOMP,NR  )
C
C     RESTORE PERTURBATION VECTOR
C
      CALL GETXY(X,Y)
C
      DEALLOCATE( XC,YC )
      RETURN
 1000 FORMAT(//,' M = ',F5.1,'  COMPONENTS',/,1X,21('='),/)
 1010 FORMAT(1H1)
      END
*DECK VACGEO
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C------ GEOMETRIC FACTORS IN VACUUM ------- A.B. 24.05.90 --------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE VACGEO
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      INCLUDE 'comioc.inc'
      INTEGER          I
C
C..DEFINE VCSH (CELL SIZE)
C
      DO 10 I=1,NV
 10   VCSH(I) = VCS(I+1) - VCS(I)
      VCSH(0) = 0.
      VCSH(NV+1) = 0.
C
      RETURN
      END
*DECK RWALLG
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C------- GET GRAD(S) = SQRT(DG22L * DG33L) ON WALL ---------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE RWALLG
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      REAL*8,DIMENSION(:),ALLOCATABLE::GRADS
      REAL*8,    DIMENSION(:,:),ALLOCATABLE::GRAD2DS
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::FGRAD2DS
      INTEGER    I,J
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'RWALLG'
C
 
      IF (.NOT. ALLOCATED(FGRADS)) THEN
         ALLOCATE(FGRADS(MEDIM))
      END IF
C
      ALLOCATE(GRADS(NCHI))
      ALLOCATE(GRAD2DS(1,NCHI),FGRAD2DS(1,MEDIM))
      I     = IWALLJ

C     GET GRAD(S) AS A FUNCTION OF ANGLE ON THE WALL
      DO J=1,NCHI 
         GRADS(J)=SQRT(ABS(VRG22L(I,J)/(VRG11L(I,J)*VRG22L(I,J)-
     &                                  VRG12L(I,J)**2)))
      GRAD2DS(1,J) = GRADS(J)
      ENDDO    
C
      NDSTRT   =  1
      NDM0     =  1
      NDM1     =  NDSTRT
      call FFTDRIVER(GRAD2DS,FGRAD2DS,FORWD,  NDM0, NDM1, NDSTRT
     &                      ,MEDIM,   NCHI,   KUOUT
     &                      ,IERSUB,  IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: GRADS in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      DO 50 J   = 1,MEDIM
      FGRADS(J) = FGRAD2DS(NDSTRT,J)
 50   CONTINUE
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(GRAD2DS,FGRAD2DS,  NDM0,    NDM1,     NDSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'FGRAD2DS')
      ENDIF
C
      DEALLOCATE(GRADS)
      DEALLOCATE(GRAD2DS,FGRAD2DS)
      RETURN
      END
*DECK VACLIN
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------LINEAR ROUTINE FOR VACUUM--G.VLAD 05/04/1989-------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....CONSTRUCT MATRIX ELEMENTS FOR LINEAR SOLVER.....................
C
      SUBROUTINE VACLIN(
C     ==================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE GVACUUMM
      USE FEEDBACKM
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
C----------------------------------------
CYQLIU  19/04/1999
C      INCLUDE 'feedback.inc'
C----------------------------------------
      INTEGER          MSA,MSB,NSA,NSB,MS
      INTEGER          MSPL,MSMI
      INTEGER          LR,LC,I,J,K,KR,KC,LX,LY,IV
      INTEGER          LXROW,LXCOL,LYROW,LYCOL
      REAL*8           RTMP
      COMPLEX*16       CMROW,CMA,CMB,CNA,SHIFT
      PARAMETER        (NSA=2,NSB=1)

      REAL*8             ZGS
C
      INCLUDE 'integc.inc'
C
      IF (NWALL.LE.0) GO TO 5
C     WRITE(*,*) ' RESISTIVE WALL INCLUDED'
      GO TO 6
 5    CONTINUE
      WRITE(*,*) ' NO RESISTIVE WALL,  NV=',NV
 6    CONTINUE
C
      DO  10 LR=1,MXMAX
        DO  10 LC=1,MXMAX
          DO 10  J=NRP1+1,NTP1
            ASUBM(LR,LC,J)= 0.
            BSUBM(LR,LC,J)= 0.
            CSUBM(LR,LC,J)= 0.
10        CONTINUE

      DO  12 LX=1,MXMAX
        DO  12 LY=1,MYMAX
          DO J=NRP1+1,NTP1
            ESUBM(LX,LY,J)= 0.
            HSUBM(LX,LY,J)= 0.
          ENDDO
          DO J=NRP1,NTP1
            FSUBM(LY,LX,J)= 0.
            GSUBM(LY,LX,J)= 0.
          ENDDO
12      CONTINUE

      DO  13 LR=1,MYMAX
        DO  13 LC=1,MYMAX
          DO 13  J=NRP1,NTP1
            DSUBM(LR,LC,J)= 0.
13        CONTINUE

C
C
C..MSA,NSA     ARE THE INDICES OF PERTURBED   QUANTITIES (NSA>1)
C..MSB,NSB     ARE THE INDICES OF EQUILIBRIUM QUANTITIES (NSB=1)
C
C
        DO 140 MSA=1,MSMAX
C
          DO 130 MSB=1,MSMAX
C      WRITE(*,*) ' VACLIN MSA MSB '
C
C..CONVOLUTIONS WITH GEOMETRICAL FACTORS ARE DONE OVER ALL M (MSB);
C..EQUILIBRIUM QUANTITIES (NSB=1) ARE INITIALIZED = ZERO FOR
C..M > MSDIM (MSMAX IS ASSUMED TO BE > MSDIM)
C
            MSPL= MPLUS(MSA,NSA,MSB,NSB)
            MSMI=MMINUS(MSA,NSA,MSB,NSB)
C
C
C..I(BASE)+LROW IS THE ROW    (EQUATION) INDEX OF THE COEFFICIENTS
C..I(BASE)+LCOL  IS THE COLUMN (VARIABLE) INDEX OF THE COEFFICIENTS
C
C
C..NOTE: NO CHECK ON NPL,NMI IS NEEDED (NSB=1)
C
C      WRITE(*,'(" MSPL ",I5)') MSPL
      IF (MSPL.LT.1) GOTO 80
C
        CMROW =(RM(MSA,NSA) + RM(MSB,NSB))*CI
        CMA   = RM(MSA,NSA)*CI
        CMB   = RM(MSB,NSB)*CI
        CNA   = RN(NSA)    *CI
        SHIFT = AL0
C
        CALL VACCOE(MSPL,MSA ,MSB ,CMROW,CMA,CMB,CNA,SHIFT,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
C
80      CONTINUE
C
              IF (MSB.LT.2) GOTO 110
              IF (MSMI.LT.1) GOTO 110
         CMROW = -(RM(MSA,NSA)-RM(MSB,NSB))*CI
         CMA   = -RM(MSA,NSA)*CI
         CMB   =  RM(MSB,NSB)*CI
         CNA   = -RN(NSA)    *CI
         SHIFT = CONJG(AL0)
C
        CALL VACCOE(MSMI,MSA ,MSB ,CMROW,CMA,CMB,CNA,SHIFT,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
C
        LXROW=(MSMI-1)*NXCOMP
        LYROW=(MSMI-1)*NYCOMP
        LXCOL=(MSA -1)*NXCOMP
        LYCOL=(MSA -1)*NYCOMP
C
        DO 100 KR = 1+LXROW, NXCOMP+LXROW
        DO 100 KC = 1+LXCOL, NXCOMP+LXCOL
         DO 100 J=NRP1+1,NTOT
           ASUBM(KR,KC,J) = CONJG(ASUBM(KR,KC,J))
           BSUBM(KR,KC,J) = CONJG(BSUBM(KR,KC,J))
           CSUBM(KR,KC,J) = CONJG(CSUBM(KR,KC,J))
100      CONTINUE
         DO 101 KR = 1+LXROW, NXCOMP+LXROW
         DO 101 KC = 1+LYCOL, NYCOMP+LYCOL
         DO 101 J=NRP1+1,NTOT
           ESUBM(KR,KC,J) = CONJG(ESUBM(KR,KC,J))
           HSUBM(KR,KC,J) = CONJG(HSUBM(KR,KC,J))
101      CONTINUE
         DO 102 KR = 1+LYROW, NYCOMP+LYROW
         DO 102 KC = 1+LXCOL, NXCOMP+LXCOL
         DO 102 J=NRP1,NTOT
           FSUBM(KR,KC,J) = CONJG(FSUBM(KR,KC,J))
           GSUBM(KR,KC,J) = CONJG(GSUBM(KR,KC,J))
102      CONTINUE
         DO 103 KR=1+LYROW,NYCOMP+LYROW
         DO 103 KC=1+LYCOL,NYCOMP+LYCOL
         DO 103 J=NRP1,NTOT
           DSUBM(KR,KC,J) = CONJG(DSUBM(KR,KC,J))
103      CONTINUE
110      CONTINUE
130     CONTINUE
140   CONTINUE

      DO 200 MSA = 1,MSMAX
C
C..ADD COEFFICIENTS WHICH DO NOT NEED CONVOLUTIONS
C
      LX=(MSA-1)*NXCOMP
      LY=(MSA-1)*NYCOMP
C
      CMA = CI*RM(MSA,NSA)
      CNA = CI*RN(NSA)
C
      IF (NWALL.LE.0) GOTO 150
      DO J=1,NWALL
         IWALLJ = IWALL(J)
         TAUWJ  = TAUW(J)
         I = IWALLJ
         IV = NR + I
C     
         IF (IWO.EQ.0) THEN
           IF (ABS(RM(MSA,2)).GT.0.1) THEN 
              BSUBM(KXB1+LX,KXB1+LX,IV) = -AL0
           ELSE
              BSUBM(KXB1+LX,KXB1+LX,IV) = -AL0*T(NRP1)
           ENDIF
         ENDIF
      ENDDO

 150  CONTINUE
      
      IF (ABS(RNTOR).GT.1.0E-13) THEN
C     DIV B = 0. AS THE EQUATION FOR B3
      DO IV=NRP1,NTOT
      I = IV - NR
      INCLUDE 'vtophat.inc'
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB3+LY,KXB1+LX ,IV)= -ZNORM
      GSUBM(KYB3+LY,KXB1+LX ,IV)=  ZNORM
      ELSE
      FSUBM(KYB3+LY,KXB1+LX ,IV)= -ZNORM*T(NRP1)
      GSUBM(KYB3+LY,KXB1+LX ,IV)=  ZNORM*T(NRP1)
      ENDIF

      DSUBM(KYB3+LY,KYB2+LY ,IV)= CMA * GG(C1,C1,C1)
      DSUBM(KYB3+LY,KYB3+LY ,IV)= CNA * GG(C1,C1,C1)
      ENDDO
      ELSEIF (ABS(RM(MSA,2)).GT.0.1) THEN
C     DIV B = 0. AS THE EQUATION FOR B2
      DO IV=NRP1,NTOT
      I = IV - NR
      INCLUDE 'vtophat.inc'
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB2+LY,KXB1+LX ,IV)= -ZNORM
      GSUBM(KYB2+LY,KXB1+LX ,IV)=  ZNORM
      ELSE
      FSUBM(KYB2+LY,KXB1+LX ,IV)= -ZNORM*T(NRP1)
      GSUBM(KYB2+LY,KXB1+LX ,IV)=  ZNORM*T(NRP1)
      ENDIF

      DSUBM(KYB2+LY,KYB2+LY ,IV)= CMA * GG(C1,C1,C1)
      ENDDO
      ENDIF

 200  CONTINUE
C
C.....SET V1, J2 AND J3 TO 0
C
      DO 210 MS = 1,MSMAX
      LXROW = (MS-1)*NXCOMP
      DO 211 IV = NRP1+1,NTP1

      RTMP = 1.
      DO J=1,NWALL
         IF (IWALL(J)+NR.EQ.IV) RTMP=0.
      ENDDO
      IF (INCFEED.GE.1.AND.INCFEED.NE.20.AND.INCFEED.NE.21) THEN
         IF (IFEED+NR.EQ.IV) RTMP=0.
      ENDIF
      IF (RTMP.GT.0..OR.ABS(RNTOR).LT.1.0E-13) 
     &BSUBM(KXJ2U+LXROW,KXJ2U+LXROW,IV) = 1.

      BSUBM(KXV1 +LXROW,KXV1 +LXROW,IV) = 1.
      IF (KXX1.GT.0)  BSUBM(KXX1 +LXROW,KXX1 +LXROW,IV) = 1.
      IF (KXJ2L.GT.0) BSUBM(KXJ2L+LXROW,KXJ2L+LXROW,IV) = 1.
      BSUBM(KXJ3 +LXROW,KXJ3 +LXROW,IV) = 1.
      IF (KXPD.GT.0)  BSUBM(KXPD +LXROW,KXPD +LXROW,IV) = 1.
      IF (KXPED.GT.0) BSUBM(KXPED +LXROW,KXPED +LXROW,IV) = 1.
      IF (KXW1.GT.0)  BSUBM(KXW1 +LXROW,KXW1 +LXROW,IV) = 1.
      IF (ABS(RNTOR).LT.1.0E-13.AND.ABS(RM(MS,2)).LT.0.1)
     &   BSUBM(KXB1 +LXROW,KXB1 +LXROW,IV) = 1.
      IF (KXJRE.GT.0)   BSUBM(KXJRE  +LXROW,KXJRE  +LXROW,IV) = 1. 
      IF (KXB2L.GT.0)   BSUBM(KXB2L  +LXROW,KXB2L  +LXROW,IV) = 1. 
      IF (KXB3L.GT.0)   BSUBM(KXB3L  +LXROW,KXB3L  +LXROW,IV) = 1. 
      IF (KXJRE2.GT.0)  BSUBM(KXJRE2 +LXROW,KXJRE2 +LXROW,IV) = 1. 
      IF (KXJRE2L.GT.0) BSUBM(KXJRE2L+LXROW,KXJRE2L+LXROW,IV) = 1. 
      IF (KXJRE3.GT.0)  BSUBM(KXJRE3 +LXROW,KXJRE3 +LXROW,IV) = 1. 
      IF (KXDPHI.GT.0)  BSUBM(KXDPHI +LXROW,KXDPHI +LXROW,IV) = 1. 
 211  CONTINUE

C----------------------------------------
CYQLIU 19/04/1999
      LX = LXROW
      DO J=1,NWALL
         I  = IWALL(J)
         IV = NR + I
         INCLUDE 'vtent.inc'
         IF (ABS(RNTOR).LT.1.0E-13) THEN
         BSUBM(KXJ2U+LX,KXJ2U+LX,IV)= FF(C1,C1,C1)
         ASUBM(KXJ2U+LX,KXJ2U+LX,IV)= FFM(C1)
         CSUBM(KXJ2U+LX,KXJ2U+LX,IV)= FFP(C1)
         ENDIF

         IF (ABS(RNTOR).LT.1.0E-13.AND.KXJ2L.GT.0) THEN
         BSUBM(KXJ2L+LX,KXJ2L+LX,IV)= FF(C1,C1,C1)
         ASUBM(KXJ2L+LX,KXJ2L+LX,IV)= FFM(C1)
         CSUBM(KXJ2L+LX,KXJ2L+LX,IV)= FFP(C1)
         ENDIF

         BSUBM(KXJ3 +LX,KXJ3 +LX,IV)= FF(C1,C1,C1)
         ASUBM(KXJ3 +LX,KXJ3 +LX,IV)= FFM(C1)
         CSUBM(KXJ3 +LX,KXJ3 +LX,IV)= FFP(C1)

         IF (ABS(RNTOR).GT.1.0E-13) THEN
         RTMP = RM(MS,2)/RNTOR
         BSUBM(KXJ3 +LX,KXJ2U +LX,IV)= FF(C1,C1,C1)*RTMP
         ASUBM(KXJ3 +LX,KXJ2U +LX,IV)= FFM(C1)*RTMP
         CSUBM(KXJ3 +LX,KXJ2U +LX,IV)= FFP(C1)*RTMP
         ENDIF
      ENDDO

      IF (INCFEED.GE.1.AND.INCFEED.NE.20.AND.INCFEED.NE.21) THEN
         I  = IFEED
         IV = NR + I
         INCLUDE 'vtent.inc'
         IF (ABS(RNTOR).LT.1.0E-13) THEN
         BSUBM(KXJ2U+LX,KXJ2U+LX,IV)= FF(C1,C1,C1)
         ASUBM(KXJ2U+LX,KXJ2U+LX,IV)= FFM(C1)
         CSUBM(KXJ2U+LX,KXJ2U+LX,IV)= FFP(C1)
         ENDIF

         IF (ABS(RNTOR).LT.1.0E-13.AND.KXJ2L.GT.0) THEN
         BSUBM(KXJ2L+LX,KXJ2L+LX,IV)= FF(C1,C1,C1)
         ASUBM(KXJ2L+LX,KXJ2L+LX,IV)= FFM(C1)
         CSUBM(KXJ2L+LX,KXJ2L+LX,IV)= FFP(C1)
         ENDIF

         BSUBM(KXJ3 +LX,KXJ3 +LX,IV)= FF(C1,C1,C1)
         ASUBM(KXJ3 +LX,KXJ3 +LX,IV)= FFM(C1)
         CSUBM(KXJ3 +LX,KXJ3 +LX,IV)= FFP(C1)

         IF (ABS(RNTOR).GT.1.0E-13) THEN
         RTMP = RM(MS,2)/RNTOR
         BSUBM(KXJ3 +LX,KXJ2U +LX,IV)= FF(C1,C1,C1)*RTMP
         ASUBM(KXJ3 +LX,KXJ2U +LX,IV)= FFM(C1)*RTMP
         CSUBM(KXJ3 +LX,KXJ2U +LX,IV)= FFP(C1)*RTMP
         ENDIF
      ENDIF
 210  CONTINUE
C----------------------------------------
C
      DO 220 MS = 1,MSMAX
      LYROW = (MS-1)*NYCOMP
      DO 220 IV = NRP1,NTOT
      DSUBM(KYV2 +LYROW,KYV2 +LYROW,IV) = 1.
      IF (KYV3.GT.0) DSUBM(KYV3 +LYROW,KYV3 +LYROW,IV) = 1.
      IF (KYX2.GT.0) DSUBM(KYX2 +LYROW,KYX2 +LYROW,IV) = 1.
      IF (KYX3.GT.0) DSUBM(KYX3 +LYROW,KYX3 +LYROW,IV) = 1.
      IF (KYW2.GT.0) DSUBM(KYW2 +LYROW,KYW2 +LYROW,IV) = 1.
      IF (KYW3.GT.0) DSUBM(KYW3 +LYROW,KYW3 +LYROW,IV) = 1.
      DSUBM(KYJ1 +LYROW,KYJ1 +LYROW,IV) = 1.
      DSUBM(KYPR +LYROW,KYPR +LYROW,IV) = 1.      
      IF (KYRHOP.GT.0)  DSUBM(KYRHOP +LYROW,KYRHOP +LYROW,IV) = 1.
      IF (KYPE.GT.0)    DSUBM(KYPE   +LYROW,KYPE   +LYROW,IV) = 1.
      IF (KYPP.GT.0)    DSUBM(KYPP   +LYROW,KYPP   +LYROW,IV) = 1.
      IF (KYPPERP.GT.0) DSUBM(KYPPERP+LYROW,KYPPERP+LYROW,IV) = 1.
      IF (KYPPARA.GT.0) DSUBM(KYPPARA+LYROW,KYPPARA+LYROW,IV) = 1.
      IF (KYJRE1.GT.0)  DSUBM(KYJRE1 +LYROW,KYJRE1 +LYROW,IV) = 1.
 220  CONTINUE

C
C.....BOUNDARY CONDITION B1U = 0
C
C     WRITE(*,*) ' BEFORE BOVACU '
         IF (NCOUPL.EQ.0.OR.ABS(NCOUPL).EQ.1.OR.NCOUPL.EQ.4.OR.
     &       NCOUPL.EQ.9) CALL BOVACU(
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
         IF (NCOUPL.EQ.2.OR.NCOUPL.EQ.5) CALL BOVACU02(
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
         IF (NCOUPL.EQ.3.OR.NCOUPL.EQ.6) CALL BOVACU03(
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
         IF (NCOUPL.EQ.8) CALL BOVACU04(
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)

      RETURN
C
      END
*DECK VACCOE
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------MATRIX COEFFICIENTS IN VACUUM ----------------A.B. 23/05/1990--
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE VACCOE(MROW,MSA ,MSB,CMROW,CMA,CMB,CNA,SHIFT,
C     ==================================================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE GVACUUMM
      USE FEEDBACKM
      INCLUDE 'specmat.inc'
C----------------------------------------
CYQLIU  19/04/1999
C      INCLUDE 'feedback.inc'
C----------------------------------------
      INTEGER       MROW,MSA,MSB,NSA,NSB
      PARAMETER     (NSA=2,NSB=1)
      INTEGER       LXCOL,LYCOL,LXROW,LYROW,I,J,K,IV
      INTEGER       IEXB3, IEXJ1
      PARAMETER     (IEXB3=1, IEXJ1=1)
      REAL*8        ZGS,ZFW
      REAL*8        Z1M,Z1P,Z3M,Z3P
      COMPLEX*16    CMROW,CMA,CMB,CNA,SHIFT

C
      INCLUDE 'integc.inc'

C
      LXROW = (MROW-1)*NXCOMP
      LXCOL = (MSA -1)*NXCOMP
      LYROW = (MROW-1)*NYCOMP
      LYCOL = (MSA -1)*NYCOMP
C
C      WRITE (*,'(" VACCOE ")')
      DO 10 IV=NRP1+1,NTOT
      I = IV - NR
      K = IV
      ZFW = 1.
      DO J=1,NWALL
         IF (I.EQ.IWALL(J)) ZFW = 0.
      ENDDO
C----------------------------------------
CYQLIU 19/04/1999
      IF (INCFEED.GE.1.AND.I.EQ.IFEED.AND.
     &    INCFEED.NE.20.AND.INCFEED.NE.21) ZFW = 0.
C----------------------------------------
C-----------------------------------------------------------------------
C.. FIRST EQUATION IN VACUUM : COVARIANT-2-COMP. OF AMPERE
C..                         (KXJ2 = 3, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
      INCLUDE 'vtent.inc'
C
      IF (I.LE.NVEQ1) THEN
      IF (ABS(RNTOR).GT.1.0E-13) THEN
      HSUBM(KXB1+LXROW,KYB2+LYCOL,K)= - CNA * ZFW *
     &     FGM(VG12L(I,MSB),VG12LM(I-1,MSB))
      ESUBM(KXB1+LXROW,KYB2+LYCOL,K)= - CNA * ZFW *
     &     FGP(VG12L(I,MSB),VG12LM(I  ,MSB))
C
      HC = VCSH(I-1)
      HSUBM(KXB1+LXROW,KYB3+LYCOL,K)= - ZFW * (
     &     GG(VG33LM(I-1,MSB),VG33L(I-1,MSB),VG33L(I,MSB))/HC
     &    +FGM(VG33JS(I,MSB),VG33JSM(I-1,MSB)) )
      HC = VCSH(I)
      ESUBM(KXB1+LXROW,KYB3+LYCOL,K)=   ZFW * (
     &     GG(VG33LM(I  ,MSB),VG33L(I,MSB),VG33L(I+1,MSB))/HC
     &    -FGP(VG33JS(I,MSB),VG33JSM(I  ,MSB)) )
C
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)= - CNA * ZFW *
     &     FF(VG11L(I,MSB),VG11LM(I-1,MSB),VG11LM(I,MSB))
      ASUBM(KXB1+LXROW,KXB1+LXCOL,K)=-CNA*ZFW*FFM(VG11LM(I-1,MSB))
      CSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-CNA*ZFW*FFP(VG11LM(I  ,MSB))
      ELSE
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)= - CNA * ZFW * T(NRP1)*
     &     FF(VG11L(I,MSB),VG11LM(I-1,MSB),VG11LM(I,MSB))
      ASUBM(KXB1+LXROW,KXB1+LXCOL,K)=-CNA*ZFW*T(NRP1)*
     &     FFM(VG11LM(I-1,MSB))
      CSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-CNA*ZFW*T(NRP1)*
     &     FFP(VG11LM(I  ,MSB))
      ENDIF
      ELSEIF (ABS(RM(MROW,2)).GT.0.1) THEN
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=ZFW*(
     &      FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     -FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*FF(VDG11L(I,MSB),VDG11LM(I-1,MSB),VDG11LM(I,MSB)))

      ASUBM(KXB1+LXROW,KXB1+LXCOL,K)=ZFW*(
     &      FGM(VDG12L(I-1,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     +CMROW*FFM(VDG11LM(I-1,MSB)))
C
      CSUBM(KXB1+LXROW,KXB1+LXCOL,K)=ZFW*(
     &     -FGP(VDG12L(I+1,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*FFP(VDG11LM(I,MSB)))
      ELSE
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=ZFW*(
     & T(NRP1)*FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &-T(NRP1)*FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &+CMROW*T(NRP1)*FF(VDG11L(I,MSB),VDG11LM(I-1,MSB),VDG11LM(I,MSB)))
C
      ASUBM(KXB1+LXROW,KXB1+LXCOL,K)=ZFW*(
     &      T(NRP1)*FGM(VDG12L(I-1,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     +CMROW*T(NRP1)*FFM(VDG11LM(I-1,MSB)))
C
      CSUBM(KXB1+LXROW,KXB1+LXCOL,K)=ZFW*(
     &     -T(NRP1)*FGP(VDG12L(I+1,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*T(NRP1)*FFP(VDG11LM(I,MSB)))
      ENDIF
C
      HC = VCSH(I-1)
      HSUBM(KXB1+LXROW,KYB2+LYCOL,K)=ZFW*(
     &     GG(VDG22LM(I-1,MSB),VDG22L(I-1,MSB),VDG22L(I,MSB))/VCSH(I-1)
     &  +  CMROW * FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB)))
      HC = VCSH(I)
      ESUBM(KXB1+LXROW,KYB2+LYCOL,K)=ZFW*(
     &  -  GG(VDG22LM(I  ,MSB),VDG22L(I,MSB),VDG22L(I+1,MSB))/VCSH(I)
     &  +  CMROW * FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB)))
      ENDIF
      ENDIF
 10   CONTINUE
C
C     WALL
C
      IF (NWALL.LE.0) GOTO 15
      DO J=1,NWALL
      IWALLJ = IWALL(J)
      TAUWJ  = TAUW(J)
      I = IWALLJ
      IV = NR + I
      IF (IV.GE.NRP1+1.AND.IV.LE.NTOT) THEN
      K = IV
      
      IF (IWO.EQ.0) THEN
C
      CALL RWALLG
      ZGS = (VCS(IWALLJ)/ASPCT)/TAUWJ * FGRADS(MSB)
C
C     --------------------
C     FACTOR OF '2' --- LIU
C     ---------------------
      HSUBM(KXB1+LXROW,KYB2+LYCOL,K) =  (CMA + CMB) * ZGS
      HSUBM(KXB1+LXROW,KYB3+LYCOL,K) =        CNA      * ZGS
      ESUBM(KXB1+LXROW,KYB2+LYCOL,K) = -(CMA + CMB) * ZGS
      ESUBM(KXB1+LXROW,KYB3+LYCOL,K) =       -CNA      * ZGS

      ENDIF

      IF (IWO.EQ.1) THEN
      INCLUDE 'vtent.inc'
      ZGS = VCS(IWALLJ)*(VCS(IWALLJ+1)-VCS(IWALLJ-1))
     &      /TAUWJ/2/ASPCT/ASPCT
      BSUBM(KXB1+LXROW,KXJ3+LXCOL,K)=-CMROW*ZGS*
     &      FF(VG33L(I,MSB),VG33LM(I-1,MSB),VG33LM(I,MSB))
     & +ZGS*FF(VG33JC(I,MSB),VG33JCM(I-1,MSB),VG33JCM(I,MSB))
      ASUBM(KXB1+LXROW,KXJ3+LXCOL,K) = 
     & -CMROW*ZGS*FFM(VG33LM(I-1,MSB))
     &       +ZGS*FFM(VG33JCM(I-1,MSB))
      CSUBM(KXB1+LXROW,KXJ3+LXCOL,K) = 
     & -CMROW*ZGS*FFP(VG33LM(I,MSB))
     &       +ZGS*FFP(VG33JCM(I,MSB))

      Z1M = (VCS(I)/VCSM(I-1))**IEXJ1
      Z1P = (VCS(I)/VCSM(I  ))**IEXJ1
      HSUBM(KXB1+LXROW,KYJ1+LYCOL,K) = CNA*ZGS*
     &          FGM(Z1M*VG12L(I,MSB),VG12LM(I-1,MSB))
      ESUBM(KXB1+LXROW,KYJ1+LYCOL,K) = CNA*ZGS*
     &          FGP(Z1P*VG12L(I,MSB),VG12LM(I  ,MSB))
      BSUBM(KXB1+LXROW,KXJ2U+LXCOL,K) = CNA*ZGS*
     &          FF(VG22L(I,MSB),VG22LM(I-1,MSB),VG22LM(I,MSB))
      ASUBM(KXB1+LXROW,KXJ2U+LXCOL,K) = CNA*ZGS*FFM(VG22LM(I-1,MSB))
      CSUBM(KXB1+LXROW,KXJ2U+LXCOL,K) = CNA*ZGS*FFP(VG22LM(I,MSB))

      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-SHIFT*
     &        FF(VJAC(I,MSB),VJACM(I-1,MSB),VJACM(I,MSB))
      ASUBM(KXB1+LXROW,KXB1+LXCOL,K) =-SHIFT*FFM(VJACM(I-1,MSB))
      CSUBM(KXB1+LXROW,KXB1+LXCOL,K) = -SHIFT*FFP(VJACM(I,MSB))
      ELSE
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-SHIFT*T(NRP1)*
     &        FF(VJAC(I,MSB),VJACM(I-1,MSB),VJACM(I,MSB))
      ASUBM(KXB1+LXROW,KXB1+LXCOL,K) =-SHIFT*T(NRP1)*FFM(VJACM(I-1,MSB))
      CSUBM(KXB1+LXROW,KXB1+LXCOL,K) = -SHIFT*T(NRP1)*FFP(VJACM(I,MSB))
      ENDIF      
      ENDIF      

      IF (IWO.EQ.2) THEN
      INCLUDE 'vtent.inc'
      ZGS = VCS(IWALLJ)*(VCS(IWALLJ+1)-VCS(IWALLJ-1))
     &      /TAUWJ/2/ASPCT/ASPCT

      CALL RWALLG2(J,1)
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-SHIFT*FGRADS(MSB)
      ELSE
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-SHIFT*T(NRP1)*FGRADS(MSB)
      ENDIF

      BSUBM(KXB1+LXROW,KXJ2L+LXCOL,K)=CNA*ZGS*ZCNDF(J,MSB)

      CALL RWALLG2(J,2)
      BSUBM(KXB1+LXROW,KXJ3+LXCOL,K)=-CMROW*ZGS*FGRADS(MSB)

      CALL RWALLG2(J,3)
      BSUBM(KXB1+LXROW,KXJ3+LXCOL,K)=BSUBM(KXB1+LXROW,KXJ3+LXCOL,K)+
     &      ZGS*2.0*FGRADS(MSB)
      ENDIF
         
      IF (IWO.EQ.4) THEN
      INCLUDE 'vtent.inc'
      ZGS = VCS(IWALLJ)*(VCS(IWALLJ+1)-VCS(IWALLJ-1))
     &      /TAUWJ/2/ASPCT/ASPCT

      CALL RWALLG2(J,1)
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-SHIFT*FGRADS(MSB)
      ELSE
      BSUBM(KXB1+LXROW,KXB1+LXCOL,K)=-SHIFT*T(NRP1)*FGRADS(MSB)
      ENDIF

      CALL RWALLG2(J,2)
      BSUBM(KXB1+LXROW,KXJ3+LXCOL,K)=-CMROW*ZGS*FGRADS(MSB)

      CALL RWALLG2(J,3)
      BSUBM(KXB1+LXROW,KXJ3+LXCOL,K)=BSUBM(KXB1+LXROW,KXJ3+LXCOL,K)+
     &      ZGS*2.0*FGRADS(MSB)
      ENDIF
         
      ENDIF         
      ENDDO
 15   CONTINUE
C
      DO  20 IV=NRP1,NTOT
      I = IV - NR
      K = IV
C-----------------------------------------------------------------------
C..SECOND EQUATION IN VACUUM : COVARIANT-1-COMP. OF AMPERE J1U=0
C..                            ( DEFINED ON HALF    MESH)
C-----------------------------------------------------------------------
C
      INCLUDE 'vtophat.inc'
C
      IF (ABS(RNTOR).GT.1.0E-13) THEN
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB2+LYROW,KXB1+LXCOL,K)= CNA *
     &                                GF(VG12L(I  ,MSB),VG12LM(I,MSB))
      GSUBM(KYB2+LYROW,KXB1+LXCOL,K)= CNA *
     &                                GF(VG12L(I+1,MSB),VG12LM(I,MSB))
      ELSE
      FSUBM(KYB2+LYROW,KXB1+LXCOL,K)= CNA * T(NRP1) *
     &                                GF(VG12L(I  ,MSB),VG12LM(I,MSB))
      GSUBM(KYB2+LYROW,KXB1+LXCOL,K)= CNA * T(NRP1) *
     &                                GF(VG12L(I+1,MSB),VG12LM(I,MSB))
      ENDIF
C
      DSUBM(KYB2+LYROW,KYB2+LYCOL,K)= CNA *
     &                  GG(VG22LM(I,MSB),VG22L(I,MSB),VG22L(I+1,MSB))
      DSUBM(KYB2+LYROW,KYB3+LYCOL,K)=
     &           -CMROW*GG(VG33LM(I,MSB),VG33L(I,MSB),VG33L(I+1,MSB))
     &                 +GG(VG33JCM(I,MSB),VG33JC(I,MSB),VG33JC(I+1,MSB))
      ELSEIF (ABS(RM(MROW,2)).LT.0.1) THEN 
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB2+LYROW,KXB1+LXCOL,K)= GF(VDG12L(I  ,MSB),VDG12LM(I,MSB))
      GSUBM(KYB2+LYROW,KXB1+LXCOL,K)= GF(VDG12L(I+1,MSB),VDG12LM(I,MSB))
      ELSE
      FSUBM(KYB2+LYROW,KXB1+LXCOL,K)= T(NRP1) *
     &                                GF(VDG12L(I  ,MSB),VDG12LM(I,MSB))
      GSUBM(KYB2+LYROW,KXB1+LXCOL,K)= T(NRP1) *
     &                                GF(VDG12L(I+1,MSB),VDG12LM(I,MSB))
      ENDIF
C
      DSUBM(KYB2+LYROW,KYB2+LYCOL,K)= 
     &                  GG(VDG22LM(I,MSB),VDG22L(I,MSB),VDG22L(I+1,MSB))
      ENDIF
C-----------------------------------------------------------------------
C
 20   CONTINUE

      IF (ABS(RNTOR).LT.1.0E-13) THEN
      DO  IV=NRP1,NTOT
      I = IV - NR
      K = IV

      INCLUDE 'vtophat.inc'
C
      DSUBM(KYB3+LYROW,KYB3+LYCOL,K)=
     &                  GG(VDG33LM(I,MSB),VDG33L(I,MSB),VDG33L(I+1,MSB))
      ENDDO
      ENDIF

C----------------------------------------
CYQLIU 19/04/1999
      DO J=1,NWALL
      I  = IWALL(J)
      IV = NR + I
      IF (IV.GE.NRP1+1.AND.IV.LE.NTOT) THEN
      K = IV
      INCLUDE 'vtent.inc'
C
      HSUBM(KXJ2U+LXROW,KYB2+LYCOL,K)=-CNA  *
     &                              FGM(VG12L(I,MSB),VG12LM(I-1,MSB))
      ESUBM(KXJ2U+LXROW,KYB2+LYCOL,K)=-CNA *
     &                              FGP(VG12L(I,MSB),VG12LM(I  ,MSB))
C
      HC = VCSH(I-1)
      HSUBM(KXJ2U+LXROW,KYB3+LYCOL,K) =
     & -GG(VG33LM(I-1,MSB),VG33L(I-1,MSB),VG33L(I,MSB))/HC
     & -FGM(VG33JS(I,MSB),VG33JSM(I-1,MSB))
C
      HC = VCSH(I)
      ESUBM(KXJ2U+LXROW,KYB3+LYCOL,K) =
     &  GG(VG33LM(I,MSB),VG33L(I,MSB),VG33L(I+1,MSB))/HC
     & -FGP(VG33JS(I,MSB),VG33JSM(I  ,MSB))
C
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA *
     &           FF(VG11L(I,MSB),VG11LM(I-1,MSB),VG11LM(I,MSB))
      ASUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA * FFM(VG11LM(I-1,MSB))
      CSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA * FFP(VG11LM(I,MSB))
      ELSE
      BSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA * T(NRP1) *
     &           FF(VG11L(I,MSB),VG11LM(I-1,MSB),VG11LM(I,MSB))
      ASUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA*T(NRP1)*FFM(VG11LM(I-1,MSB))
      CSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA*T(NRP1)*FFP(VG11LM(I,MSB))
      ENDIF

      BSUBM(KXJ2U+LXROW,KXJ2U+LXCOL,K)=
     &           FF(VJAC(I,MSB),VJACM(I-1,MSB),VJACM(I,MSB))
      ASUBM(KXJ2U+LXROW,KXJ2U+LXCOL,K)=FFM(VJACM(I-1,MSB))
      CSUBM(KXJ2U+LXROW,KXJ2U+LXCOL,K)=FFP(VJACM(I,MSB))
      ENDIF
      ENDDO

      IF (INCFEED.GE.1.AND.INCFEED.NE.20.AND.INCFEED.NE.21) THEN
      I  = IFEED
      IV = NR + I
      IF (IV.GE.NRP1.AND.IV.LE.NTOT) THEN
      K = IV
      INCLUDE 'vtent.inc'
C
      HSUBM(KXJ2U+LXROW,KYB2+LYCOL,K)=-CNA  *
     &                              FGM(VG12L(I,MSB),VG12LM(I-1,MSB))
      ESUBM(KXJ2U+LXROW,KYB2+LYCOL,K)=-CNA *
     &                              FGP(VG12L(I,MSB),VG12LM(I  ,MSB))
C
      HC = VCSH(I-1)
      HSUBM(KXJ2U+LXROW,KYB3+LYCOL,K) =
     & -GG(VG33LM(I-1,MSB),VG33L(I-1,MSB),VG33L(I,MSB))/HC
     & -FGM(VG33JS(I,MSB),VG33JSM(I-1,MSB))
C
      HC = VCSH(I)
      ESUBM(KXJ2U+LXROW,KYB3+LYCOL,K) =
     &  GG(VG33LM(I,MSB),VG33L(I,MSB),VG33L(I+1,MSB))/HC
     & -FGP(VG33JS(I,MSB),VG33JSM(I  ,MSB))
C
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA *
     &           FF(VG11L(I,MSB),VG11LM(I-1,MSB),VG11LM(I,MSB))
      ASUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA * FFM(VG11LM(I-1,MSB))
      CSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA * FFP(VG11LM(I,MSB))
      ELSE
      BSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA * T(NRP1) *
     &           FF(VG11L(I,MSB),VG11LM(I-1,MSB),VG11LM(I,MSB))
      ASUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA*T(NRP1)*FFM(VG11LM(I-1,MSB))
      CSUBM(KXJ2U+LXROW,KXB1+LXCOL,K)=-CNA*T(NRP1)*FFP(VG11LM(I,MSB))
      ENDIF

      BSUBM(KXJ2U+LXROW,KXJ2U+LXCOL,K)=
     &           FF(VJAC(I,MSB),VJACM(I-1,MSB),VJACM(I,MSB))
      ASUBM(KXJ2U+LXROW,KXJ2U+LXCOL,K)=FFM(VJACM(I-1,MSB))
      CSUBM(KXJ2U+LXROW,KXJ2U+LXCOL,K)=FFP(VJACM(I,MSB))
      ENDIF
      ENDIF
C----------------------------------------
C----------------------------------------
CYQLIU 19/04/1999
      IF (ABS(RNTOR).GT.1.0E-13) GOTO 2013 
      DO J=1,NWALL
         I  = IWALL(J)
         IV = NR + I
         IF (IV.GE.NRP1.AND.IV.LE.NTOT) THEN
         K = IV
         INCLUDE 'vtent.inc'
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &      FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     -FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*FF(VDG11L(I,MSB),VDG11LM(I-1,MSB),VDG11LM(I,MSB))
C
      ASUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &      FGM(VDG12L(I-1,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     +CMROW*FFM(VDG11LM(I-1,MSB))
C
      CSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &     -FGP(VDG12L(I+1,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*FFP(VDG11LM(I,MSB))
      ELSE
      BSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     & T(NRP1)*FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &-T(NRP1)*FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &+CMROW*T(NRP1)*FF(VDG11L(I,MSB),VDG11LM(I-1,MSB),VDG11LM(I,MSB))
C
      ASUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &      T(NRP1)*FGM(VDG12L(I-1,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     +CMROW*T(NRP1)*FFM(VDG11LM(I-1,MSB))
C
      CSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &     -T(NRP1)*FGP(VDG12L(I+1,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*T(NRP1)*FFP(VDG11LM(I,MSB))
      ENDIF
C
      HC = VCSH(I-1)
      HSUBM(KXJ3+LXROW,KYB2+LYCOL,K)=
     &     GG(VDG22LM(I-1,MSB),VDG22L(I-1,MSB),VDG22L(I,MSB))/VCSH(I-1)
     &  +  CMROW * FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))
      HC = VCSH(I)
      ESUBM(KXJ3+LXROW,KYB2+LYCOL,K)=
     &  -  GG(VDG22LM(I  ,MSB),VDG22L(I,MSB),VDG22L(I+1,MSB))/VCSH(I)
     &  +  CMROW * FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))
      ENDIF
      ENDDO

      IF (INCFEED.GE.1.AND.INCFEED.NE.20.AND.INCFEED.NE.21) THEN
         I  = IFEED
         IV = NR + I
         IF (IV.GE.NRP1.AND.IV.LE.NTOT) THEN
         K = IV
         INCLUDE 'vtent.inc'
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      BSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &      FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     -FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*FF(VDG11L(I,MSB),VDG11LM(I-1,MSB),VDG11LM(I,MSB))
C
      ASUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &      FGM(VDG12L(I-1,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     +CMROW*FFM(VDG11LM(I-1,MSB))
C
      CSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &     -FGP(VDG12L(I+1,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*FFP(VDG11LM(I,MSB))
      ELSE
      BSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &      T(NRP1)*FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &     -T(NRP1)*FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*T(NRP1)*FF(VDG11L(I,MSB),VDG11LM(I-1,MSB),
     &                       VDG11LM(I,MSB))
C
      ASUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &      T(NRP1)*FGM(VDG12L(I-1,MSB),VDG12LM(I-1,MSB))/VCSH(I-1)
     &      +CMROW*T(NRP1)*FFM(VDG11LM(I-1,MSB))
C
      CSUBM(KXJ3+LXROW,KXB1+LXCOL,K)=
     &     -T(NRP1)*FGP(VDG12L(I+1,MSB),VDG12LM(I  ,MSB))/VCSH(I)
     &     +CMROW*T(NRP1)*FFP(VDG11LM(I,MSB))
      ENDIF
C
      HC = VCSH(I-1)
      HSUBM(KXJ3+LXROW,KYB2+LYCOL,K)=
     &     GG(VDG22LM(I-1,MSB),VDG22L(I-1,MSB),VDG22L(I,MSB))/VCSH(I-1)
     &  +  CMROW * FGM(VDG12L(I,MSB),VDG12LM(I-1,MSB))
      HC = VCSH(I)
      ESUBM(KXJ3+LXROW,KYB2+LYCOL,K)=
     &  -  GG(VDG22LM(I  ,MSB),VDG22L(I,MSB),VDG22L(I+1,MSB))/VCSH(I)
     &  +  CMROW * FGP(VDG12L(I,MSB),VDG12LM(I  ,MSB))
      ENDIF
      ENDIF

 2013 CONTINUE
C----------------------------------------
C----------------------------------------
CYQLIU 19/04/1999
      IF (KXJ2L.GT.0) THEN
      DO J=1,NWALL
         I  = IWALL(J)
         IV = NR + I
         IF (IV.GE.NRP1.AND.IV.LE.NTOT) THEN
         K = IV
         INCLUDE 'vtent.inc'
         Z1M = (VCS(I)/VCSM(I-1))**IEXJ1
         Z1P = (VCS(I)/VCSM(I  ))**IEXJ1
C
         HSUBM(KXJ2L+LXROW,KYJ1+LYCOL,K) =
     &        - FGM(Z1M*VG12L(I,MSB),VG12LM(I-1,MSB))
         ESUBM(KXJ2L+LXROW,KYJ1+LYCOL,K) =
     &        - FGP(Z1P*VG12L(I,MSB),VG12LM(I  ,MSB))
C     
         BSUBM(KXJ2L+LXROW,KXJ2U+LXCOL,K) =
     &        - FF(VG22L(I,MSB),VG22LM(I-1,MSB),VG22LM(I,MSB))
         ASUBM(KXJ2L+LXROW,KXJ2U+LXCOL,K) = - FFM(VG22LM(I-1,MSB))
         CSUBM(KXJ2L+LXROW,KXJ2U+LXCOL,K) = - FFP(VG22LM(I,MSB))
         ENDIF
      ENDDO
      ENDIF
      
      IF (INCFEED.GE.1.AND.INCFEED.NE.20.AND.INCFEED.NE.21.AND.
     &    KXJ2L.GT.0) THEN
         I  = IFEED
         IV = NR + I
         IF (IV.GE.NRP1.AND.IV.LE.NTOT) THEN
         K = IV
         INCLUDE 'vtent.inc'
         Z1M = (VCS(I)/VCSM(I-1))**IEXJ1
         Z1P = (VCS(I)/VCSM(I  ))**IEXJ1
C
         HSUBM(KXJ2L+LXROW,KYJ1+LYCOL,K) =
     &        - FGM(Z1M*VG12L(I,MSB),VG12LM(I-1,MSB))
         ESUBM(KXJ2L+LXROW,KYJ1+LYCOL,K) =
     &        - FGP(Z1P*VG12L(I,MSB),VG12LM(I  ,MSB))
C     
         BSUBM(KXJ2L+LXROW,KXJ2U+LXCOL,K) =
     &        - FF(VG22L(I,MSB),VG22LM(I-1,MSB),VG22LM(I,MSB))
         ASUBM(KXJ2L+LXROW,KXJ2U+LXCOL,K) = - FFM(VG22LM(I-1,MSB))
         CSUBM(KXJ2L+LXROW,KXJ2U+LXCOL,K) = - FFP(VG22LM(I,MSB))
         ENDIF
      ENDIF
C----------------------------------------
      RETURN
      END
*DECK LOCONT
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C----LOCATE CONTINUUM SINGULARITY---A.B.   10.10.89---------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE LOCONT(Y)
      USE DIMENSIM
      USE GLOBALM
      COMPLEX*16    Y(NYCOMP,MSMAX,*),FACT
      INTEGER    MS,I,IS
      REAL*8,DIMENSION(:),ALLOCATABLE::V3,SING
C
      ALLOCATE( V3(NRP1),SING(NRP1))
C
      WRITE(*,1010)
      FACT = 1.
      IF (DIMAG(Y(KYV3, 1,NR/2))**2.GT.DREAL(Y(KYV3, 1,NR/2))**2)
     &    FACT=(0.,1.)
      DO 200 MS=1,MSMAX
      DO  10 I=1,NR
      Y(KYV3,MS,I) = FACT*Y(KYV3,MS,I)
      V3(I) = DREAL(Y(KYV3 ,MS,I))
 10   CONTINUE
C
      IS = 0
      DO 50 I=2,NR
      IF (V3(I-1)*V3(I).GE.0.) GOTO 50
      IF (I.EQ.2) GOTO 20
      IF (ABS(V3(I-1)).LT.ABS(V3(I-2))) GOTO 50
 20   IF (I.EQ.NR) GOTO 30
      IF (ABS(V3(I)).LT.ABS(V3(I+1))) GOTO 50
 30   IS = IS + 1
      SING(IS) = (CSM(I)*V3(I) - CSM(I-1)*V3(I-1))/(V3(I) - V3(I-1))
 50   CONTINUE
      IF (IS.EQ.0) GOTO 200
      WRITE(*,1000) RM(MS,2),(SING(I),I=1,IS)
 200  CONTINUE
      WRITE(*,1020)
      DEALLOCATE(V3,SING)
      RETURN
 1000 FORMAT(' M = ',F3.6,3X,(4F16.5))
 1010 FORMAT(' CONTINUUM FREQUENCIES',//)
 1020 FORMAT(///)
      END
*DECK FIXORI
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------FIRST GRIDPOINT-----------A.B.   19.07.89---------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE FIXORI
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      INTEGER    MS,M
      REAL*8       ZA,ZF
      PARAMETER  (ZF=0.001)
C
      DO 100 MS = 1,MEDIM
      M = MS - 1
      IF (M.GT.1) M = M-2
      ZA = ZF**M
      CS(1)        = ZF*      CS(2)
C
      DG11L(1,MS)  = ZA/ZF*   DG11L(2,MS)
      DG22L(1,MS)  = ZA*ZF*   DG22L(2,MS)
      DG33L(1,MS)  = ZA/ZF*   DG33L(2,MS)
      DG12L(1,MS)  = ZA*      DG12L(2,MS)
C
      JACOBI(1,MS) = ZA*ZF*   JACOBI(2,MS)
C
 100  CONTINUE
      RETURN
      END
*DECK LINEAR
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------LINEAR ROUTINE-------------G.VLAD 05/04/1989-------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....CONSTRUCT MATRIX ELEMENTS FOR LINEAR SOLVER.....................
C
      SUBROUTINE LINEAR(
C     ==================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE KINETICM
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'

      INTEGER          MSA,MSB,NSA,NSB
      INTEGER          MSPL,MSMI
      INTEGER          LR,LC,I,KR,KC,LX,LY
      INTEGER          LXROW,LXCOL,LYROW,LYCOL
      REAL*8           RTMP
      COMPLEX*16       CMROW,CMA,CMB,CNA,SHIFT
      PARAMETER        (NSA=2,NSB=1)
C
      INCLUDE 'integc.inc'
C
C
C.....CALL VACUUM ROUTINE FIRST BECAUSE IT USES MATRICES
C.....ASUBM TROUGH HSUBM
C
      IF (NV.GE.2.AND.NTP1.GE.NRP1.AND.KJPKEY.GT.0) CALL VACLIN(
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)

      WRITE (*,'("AFTER VACLIN")')

C     BE CAUTIOUS WITH THE OPTION CALPHA2.NE.CALPHA3: ROT(M) ARE USED 
C     IN OTHER PLACES AS WELL
      IF ((NCASE.EQ.6.OR.NCASE.EQ.7).AND.ABS(CALPHA2-CALPHA3).GE.1.0E-3)
     &THEN
         RTMP  = CALPHA3/CALPHA2
         ROT   = ROTEQ   + TROTI*RTMP
         ROTM  = ROTEQM  + TROTM*RTMP
         DROT  = DROTEQ  + TDROTI*RTMP
         DROTM = DROTEQM + TDROTM*RTMP
         IF (NDENEQ.EQ.1) THEN
            RHO   = RHOEQ   + TRHOI*RTMP
            RHOM  = RHOEQM  + TRHOM*RTMP
         ENDIF
      ENDIF
      CALL PLASMALIN(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
      IF ((NCASE.EQ.6.OR.NCASE.EQ.7).AND.ABS(CALPHA2-CALPHA3).GE.1.0E-3)
     &THEN
         ROT   = ROTEQ  
         ROTM  = ROTEQM
         DROT  = DROTEQ
         DROTM = DROTEQM
         RHO   = RHOEQ
         RHOM  = RHOEQM
      ENDIF
C
      CALL BOWALL(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
      CALL BOUNDC(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
      IF (NCOUPL.EQ.-2.OR.NCOUPL.EQ.-3)
     &   CALL BOUNDB(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
C
C----------------------------------------
CYQLIU 14/04/1999
      CALL FEEDM(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
C----------------------------------------

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------PLASMALIN ROUTINE-------------Y.Q. LIU 11/03/2005--------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C.....CONSTRUCT MATRIX ELEMENTS IN PLASMA FOR LINEAR SOLVER.............
C
      SUBROUTINE PLASMALIN(
C     ==================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE RCOMDM
      USE GVACUUMM
      USE CONVOLCOFM
      USE REORBITM
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'

      INTEGER          MSA,MSB,NSA,NSB
      INTEGER          MSPL,MSMI
      INTEGER          LR,LC,I,J,K,KR,KC,LX,LY
      INTEGER          LXROW,LXCOL,LYROW,LYCOL
      COMPLEX*16       CMROW,CMA,CMB,CNA,SHIFT,CTMP1,CTMP2,CTMP3
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::SHIFTC,SHIFTM,
     &                 SHIFTVC,SHIFTVM,SHIFTBC,SHIFTBM
      PARAMETER        (NSA=2,NSB=1)
      integer          iexe,iexv2,iexj1,iexb3
      parameter        (iexv2 = -1, iexj1 = 1, iexb3 = 1)
      real*8           zem,zep,zv2m,zv2p,z1m,z1p,z3m,z3p,zos,zobi,zobe
 
      INCLUDE 'integc.inc'
 
      ALLOCATE(SHIFTC(NRP1),SHIFTM(NRP1),SHIFTVC(NRP1),SHIFTVM(NRP1),
     &         SHIFTBC(NRP1),SHIFTBM(NRP1))

      zos = 0.
      IF (IDIAMV.GE.1) zos = FDIAMV

      DO  10 LR=1,MXMAX
        DO  10 LC=1,MXMAX
          DO 10  J=1,NRP1
            ASUBM(LR,LC,J)= 0.
            BSUBM(LR,LC,J)= 0.
            CSUBM(LR,LC,J)= 0.
10        CONTINUE
      DO  12 LX=1,MXMAX
        DO  12 LY=1,MYMAX
          DO J=1,NRP1
            ESUBM(LX,LY,J)= 0.
            HSUBM(LX,LY,J)= 0.
          ENDDO
          DO J=1,NR
            FSUBM(LY,LX,J)= 0.
            GSUBM(LY,LX,J)= 0.
          ENDDO
12      CONTINUE
      DO  13 LR=1,MYMAX
        DO  13 LC=1,MYMAX
          DO 13  J=1,NR
            DSUBM(LR,LC,J)= 0.
13        CONTINUE
C
C
C..MSA,NSA     ARE THE INDICES OF PERTURBED   QUANTITIES (NSA>1)
C..MSB,NSB     ARE THE INDICES OF EQUILIBRIUM QUANTITIES (NSB=1)
C
C
        DO 140 MSA=1,MSMAX
C
          DO 130 MSB=1,MSMAX
C
C..CONVOLUTIONS WITH GEOMETRICAL FACTORS ARE DONE OVER ALL M (MSB);
C..EQUILIBRIUM QUANTITIES (NSB=1) ARE INITIALIZED = ZERO FOR
C..M > MSDIM (MSMAX IS ASSUMED TO BE > MSDIM)
C
            MSPL= MPLUS(MSA,NSA,MSB,NSB)
            MSMI=MMINUS(MSA,NSA,MSB,NSB)
C
C
C..I(BASE)+LROW IS THE ROW    (EQUATION) INDEX OF THE COEFFICIENTS
C..I(BASE)+LCOL  IS THE COLUMN (VARIABLE) INDEX OF THE COEFFICIENTS
C
C
C..NOTE: NO CHECK ON NPL,NMI IS NEEDED (NSB=1)
C
      IF (MSPL.LT.1) GOTO 80
C
        CMROW =(RM(MSA,NSA) + RM(MSB,NSB))*CI
        CMA   = RM(MSA,NSA)*CI
        CMB   = RM(MSB,NSB)*CI
        CNA   = RN(NSA)    *CI
        SHIFT = - AL0
        DO 228 I=1,NRP1
        SHIFTC(I) =SHIFT-RNTOR*ROT(I)*CI
 228    SHIFTVC(I)=SHIFTC(I)-RNTOR*OMEGASI(I)*CI*zos
        DO 229 I=1,NR
        SHIFTM(I) =SHIFT-RNTOR*ROTM(I)*CI
 229    SHIFTVM(I)=SHIFTM(I)-RNTOR*OMEGASIM(I)*CI*zos
        SHIFTM(NRP1)  = 2.*SHIFTC(NRP1) - SHIFTM(NR)
        SHIFTVM(NRP1) = 2.*SHIFTVC(NRP1) - SHIFTVM(NR)
 
        IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10) THEN
           FHATV(1:NRP1,MSB,1:6)  = FHATVP(1:NRP1,MSB,1:6)
     &                              *CALPHA3/CALPHA2
           FHATVM(1:NRP1,MSB,1:6) = FHATVPM(1:NRP1,MSB,1:6)
     &                              *CALPHA3/CALPHA2
        ENDIF

        CALL COEFFI(MSPL,MSA ,MSB ,CMROW,CMA,CMB,CNA,SHIFT,
     &   SHIFTC,SHIFTM,SHIFTVC,SHIFTVM,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
 
80      CONTINUE
 
        IF (MSB.LT.2) GOTO 110
        IF (MSMI.LT.1) GOTO 110
         CMROW = -(RM(MSA,NSA)-RM(MSB,NSB))*CI
         CMA   = -RM(MSA,NSA)*CI
         CMB   =  RM(MSB,NSB)*CI
         CNA   = -RN(NSA)    *CI
         SHIFT = - CONJG(AL0)
        DO 227 I=1,NRP1
        SHIFTC(I) =SHIFT+RNTOR*ROT(I)*CI
 227    SHIFTVC(I)=SHIFTC(I)+RNTOR*OMEGASI(I)*CI*zos
        DO 226 I=1,NR
        SHIFTM(I) =SHIFT+RNTOR*ROTM(I)*CI
 226    SHIFTVM(I)=SHIFTM(I)+RNTOR*OMEGASIM(I)*CI*zos
        SHIFTM(NRP1)  = 2.*SHIFTC(NRP1) - SHIFTM(NR)
        SHIFTVM(NRP1) = 2.*SHIFTVC(NRP1) - SHIFTVM(NR)
 
        IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10) THEN
           FHATV(1:NRP1,MSB,1:6)  = FHATVN(1:NRP1,MSB,1:6)
     &                              *CALPHA3/CALPHA2
           FHATVM(1:NRP1,MSB,1:6) = FHATVNM(1:NRP1,MSB,1:6)
     &                              *CALPHA3/CALPHA2
        ENDIF

        CALL COEFFI(MSMI,MSA ,MSB ,CMROW,CMA,CMB,CNA,SHIFT,
     &   SHIFTC,SHIFTM,SHIFTVC,SHIFTVM,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
C
        LXROW=(MSMI-1)*NXCOMP
        LYROW=(MSMI-1)*NYCOMP
        LXCOL=(MSA -1)*NXCOMP
        LYCOL=(MSA -1)*NYCOMP
C
         DO 100 KR=1+LXROW,NXCOMP+LXROW
         DO 100 KC=1+LXCOL,NXCOMP+LXCOL
         DO 100 J=1,NRP1
           ASUBM(KR,KC,J) = CONJG(ASUBM(KR,KC,J))
           BSUBM(KR,KC,J) = CONJG(BSUBM(KR,KC,J))
           CSUBM(KR,KC,J) = CONJG(CSUBM(KR,KC,J))
100      CONTINUE
         DO 101 KR=1+LXROW,NXCOMP+LXROW
         DO 101 KC=1+LYCOL,NYCOMP+LYCOL
         DO 101 J=1,NRP1
           ESUBM(KR,KC,J) = CONJG(ESUBM(KR,KC,J))
           HSUBM(KR,KC,J) = CONJG(HSUBM(KR,KC,J))
101      CONTINUE
         DO 102 KR=1+LYROW,NYCOMP+LYROW
         DO 102 KC=1+LXCOL,NXCOMP+LXCOL
         DO 102 J=1,NR
           FSUBM(KR,KC,J) = CONJG(FSUBM(KR,KC,J))
           GSUBM(KR,KC,J) = CONJG(GSUBM(KR,KC,J))
102      CONTINUE
         DO 103 KR=1+LYROW,NYCOMP+LYROW
         DO 103 KC=1+LYCOL,NYCOMP+LYCOL
         DO 103 J=1,NR
           DSUBM(KR,KC,J) = CONJG(DSUBM(KR,KC,J))
103      CONTINUE
110      CONTINUE
130     CONTINUE
140   CONTINUE
C
C
C..ADD COEFFICIENTS WHICH DO NOT NEED CONVOLUTIONS
C

      zobi = 0.
      zobe = 0.
      IF (IDIAMB.EQ.5.OR.IDIAMB.EQ.6) zobe = FDIAMB
      IF (IDIAMB.EQ.6) zobi = FDIAMB
      SHIFT = - AL0
      DO I=1,NRP1
         SHIFTC(I)=SHIFT-RNTOR*ROT(I)*CI
         SHIFTBC(I)=SHIFTC(I)-RNTOR*OMEGASE(I)*CI*zobe
         SHIFTBC(I)=SHIFTBC(I)+RNTOR*OMEGASI(I)*CI*zobi
      ENDDO
      DO I=1,NR
         SHIFTM(I)=SHIFT-RNTOR*ROTM(I)*CI
         SHIFTBM(I)=SHIFTM(I)-RNTOR*OMEGASEM(I)*CI*zobe
         SHIFTBM(I)=SHIFTBM(I)+RNTOR*OMEGASIM(I)*CI*zobi
      ENDDO
      SHIFTM(NRP1) = 2.*SHIFTC(NRP1) - SHIFTM(NR)
      SHIFTBM(NRP1) = 2.*SHIFTBC(NRP1) - SHIFTBM(NR)

      DO 200 MSA = 1,MSMAX
      LX=(MSA-1)*NXCOMP
      LY=(MSA-1)*NYCOMP
 
      CMA = CI*RM(MSA,NSA)
      CNA = CI*RN(NSA)

      DO I=2,NRP1
      INCLUDE 'tent.inc'
      IF (KXJ2L.GT.0) THEN
      BSUBM(KXJ2L+LX,KXJ2L+LX,I)= FF(C1,C1,C1)
      ASUBM(KXJ2L+LX,KXJ2L+LX,I)= FFM(C1)
      CSUBM(KXJ2L+LX,KXJ2L+LX,I)= FFP(C1)
      ENDIF

      IF (KXJRE2L.GT.0) THEN
      BSUBM(KXJRE2L+LX,KXJRE2L+LX,I)= FF(C1,C1,C1)
      ASUBM(KXJRE2L+LX,KXJRE2L+LX,I)= FFM(C1)
      CSUBM(KXJRE2L+LX,KXJRE2L+LX,I)= FFP(C1)
      ENDIF

      IF (KXB2L.GT.0) THEN
      BSUBM(KXB2L+LX,KXB2L+LX,I)= FF(C1,C1,C1)
      ASUBM(KXB2L+LX,KXB2L+LX,I)= FFM(C1)
      CSUBM(KXB2L+LX,KXB2L+LX,I)= FFP(C1)
      ENDIF

      IF (KXB3L.GT.0) THEN
      BSUBM(KXB3L+LX,KXB3L+LX,I)= FF(C1,C1,C1)
      ASUBM(KXB3L+LX,KXB3L+LX,I)= FFM(C1)
      CSUBM(KXB3L+LX,KXB3L+LX,I)= FFP(C1)
      ENDIF

      IF (KXJRE.GT.0) THEN
      BSUBM(KXJRE+LX,KXJRE+LX,I)= FF(SHIFT,SHIFT,SHIFT)
      ASUBM(KXJRE+LX,KXJRE+LX,I)= FFM(SHIFT)
      CSUBM(KXJRE+LX,KXJRE+LX,I)= FFP(SHIFT)

C     ADD RADIAL DIFFUSION FOR KXJRE
      IF ((KJRER.EQ.5.OR.KJRER.EQ.6).AND.I.LT.NRP1) THEN
         CTMP1 = CSH(I-1)*CSH(I)*(CSH(I-1)+CSH(I))/2.
         CTMP2 = JRE_DIFF*ABS(CMA)**2/CTMP1
         BSUBM(KXJRE+LX,KXJRE+LX,I) = BSUBM(KXJRE+LX,KXJRE+LX,I) -
     &                                (CSH(I-1)+CSH(I))*CTMP2
         ASUBM(KXJRE+LX,KXJRE+LX,I) = ASUBM(KXJRE+LX,KXJRE+LX,I) +
     &                                CSH(I)*CTMP2
         CSUBM(KXJRE+LX,KXJRE+LX,I) = CSUBM(KXJRE+LX,KXJRE+LX,I) +
     &                                CSH(I-1)*CTMP2
      ENDIF
      ENDIF

C     ADD RADIAL DIFFUSION FOR KXV1
      IF (V1U_DIFF.GT.1.E-15.AND.I.LT.NRP1.AND.ABS(RNTOR).LT.1.E-9) THEN
         CTMP1 = CSH(I-1)*CSH(I)*(CSH(I-1)+CSH(I))/2.
         CTMP2 = V1U_DIFF/CS(I)/CTMP1
         IF (I.EQ.50.AND.MSA.EQ.25) THEN
            WRITE(*,*) 'V1U_DIFF_0:',BSUBM(KXV1+LX,KXV1+LX,I)
         ENDIF 
         BSUBM(KXV1+LX,KXV1+LX,I) = BSUBM(KXV1+LX,KXV1+LX,I) -
     &                                (CSH(I-1)+CSH(I))*CTMP2
         ASUBM(KXV1+LX,KXV1+LX,I) = ASUBM(KXV1+LX,KXV1+LX,I) +
     &                                CSH(I)*CTMP2
         CSUBM(KXV1+LX,KXV1+LX,I) = CSUBM(KXV1+LX,KXV1+LX,I) +
     &                                CSH(I-1)*CTMP2
         IF (I.EQ.50.AND.MSA.EQ.25) THEN
            WRITE(*,*) 'V1U_DIFF_1:',BSUBM(KXV1+LX,KXV1+LX,I)
         ENDIF 
      ENDIF

      IF (KXPD.GT.0) THEN
      BSUBM(KXPD+LX,KXPD+LX,I)= FF(C1,C1,C1)
      ASUBM(KXPD+LX,KXPD+LX,I)= FFM(C1)
      CSUBM(KXPD+LX,KXPD+LX,I)= FFP(C1)
      ENDIF

      IF (KXPED.GT.0) THEN
      BSUBM(KXPED+LX,KXPED+LX,I)= FF(C1,C1,C1)
      ASUBM(KXPED+LX,KXPED+LX,I)= FFM(C1)
      CSUBM(KXPED+LX,KXPED+LX,I)= FFP(C1)
      ENDIF

      Z3M  = (CSM(I-1)/CS(I))**IEXB3
      Z3P  = (CSM(I)  /CS(I))**IEXB3

      CTMP1 = -DPsiDs(I)*Tp(I)
      CTMP2 = -DPsiDsM(I-1)*TpM(I-1)
      CTMP3 = -DPsiDsM(I)*TpM(I)
      HsubM(kxV1+lx, kyB3+ly,I)=FGM(CTMP1/z3m,CTMP2)
      EsubM(kxV1+lx, kyB3+ly,I)=FGP(CTMP1/z3p,CTMP3)

      CTMP1 = DPsiDs(I)
      CTMP2 = DPsiDsM(I-1)
      CTMP3 = DPsiDsM(I)
      BsubM(kxV1+lx, kxJ3+lx,I)=-FF(CTMP1,CTMP2, CTMP3)
      AsubM(kxV1+lx, kxJ3+lx,I)=-FFM(CTMP2)
      CsubM(kxV1+lx, kxJ3+lx,I)=-FFP(CTMP3)
      ENDDO

      IEXE = -1
      DO I=1,NR
      INCLUDE 'tophat.inc'
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE

      IF (KYV3.GT.0) THEN
      CTMP1 = PPeq(i)*DPsiDs(i)
      CTMP2 = PPeq(i+1)*DPsiDs(i+1)
      CTMP3 = PPeqM(i)*DPsiDsM(i)
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
      FSUBM(KYV3+LY,KXB1+LX,I)=-GF(ZEM*CTMP1,CTMP3)
      GSUBM(KYV3+LY,KXB1+LX,I)=-GF(ZEP*CTMP2,CTMP3)
      ELSE
      FSUBM(KYV3+LY,KXB1+LX,I)=-GF(ZEM*CTMP1*T(I),CTMP3*TM(I))
      GSUBM(KYV3+LY,KXB1+LX,I)=-GF(ZEP*CTMP2*T(I+1),CTMP3*TM(I))
      ENDIF

      DSUBM(KYV3+LY,KYPR+LY,I)=DSUBM(KYV3+LY,KYPR+LY,I)
     &  -CMA*GG(C1*DPsiDsM(I),C1*ZEM*DPsiDs(I),C1*ZEP*DPsiDs(I+1))
      IF (KYPE.GT.0.AND.INCKIN.EQ.0) 
     &DSUBM(KYV3+LY,KYPE+LY,I)=DSUBM(KYV3+LY,KYPE+LY,I)
     &  -CMA*GG(C1*DPsiDsM(I),C1*ZEM*DPsiDs(I),C1*ZEP*DPsiDs(I+1))

      ENDIF

      IF (KYW3.GT.0) 
     &DSUBM(KYW3+LY,KYPR+LY,I)=DSUBM(KYW3+LY,KYPR+LY,I)-CMA*RHOM(I)*
     &       GG(C1*DPsiDsM(I),C1*ZEM*DPsiDs(I),C1*ZEP*DPsiDs(I+1))
      ENDDO

      DO I=2,NRP1
      INCLUDE 'tent.inc'
      CTMP1 = DPsiDs(I)
      CTMP2 = DPsiDsM(I-1)
      CTMP3 = DPsiDsM(I)
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         IF (KVSQL(8)) THEN
         BsubM(kxB1+lx, kxV1+lx,I)=BsubM(kxB1+lx, kxV1+lx,I) +
     &     cma*FF(CTMP1,CTMP2, CTMP3)
         AsubM(kxB1+lx, kxV1+lx,I)=AsubM(kxB1+lx, kxV1+lx,I) +
     &     cma*FFM(CTMP2)
         CsubM(kxB1+lx, kxV1+lx,I)=CsubM(kxB1+lx, kxV1+lx,I) +
     &     cma*FFP(CTMP3)
         ENDIF

         IF (KXJ2L.GT.0) THEN
         BSUBM(KXB1+LX,KXJ2L+LX,I)=
     &     CNA*FF(C1*RESIST(I),C1*RESISM(I-1),C1*RESISM(I))
         ASUBM(KXB1+LX ,KXJ2L+LX,I)= CNA*FFM(C1*RESISM(I-1))
         CSUBM(KXB1+LX ,KXJ2L+LX,I)= CNA*FFP(C1*RESISM(I))
         ENDIF

         IF (KXJRE2L.GT.0) THEN
         BSUBM(KXB1+LX,KXJRE2L+LX,I)=-CNA*JRE_EQFRAC
     &        *FF(C1*RESIST(I),C1*RESISM(I-1),C1*RESISM(I))
         ASUBM(KXB1+LX ,KXJRE2L+LX,I)=-CNA*JRE_EQFRAC
     &        *FFM(C1*RESISM(I-1))
         CSUBM(KXB1+LX ,KXJRE2L+LX,I)=-CNA*JRE_EQFRAC
     &        *FFP(C1*RESISM(I))
         ENDIF
      ELSE
         IF (KXJ2L.GT.0) THEN
         BSUBM(KXB1+LX,KXJ2L+LX,I)=CNA*FF(C1*RESIST(I)/T(I),
     &     C1*RESISM(I-1)/TM(I-1),C1*RESISM(I)/TM(I))
         ASUBM(KXB1+LX,KXJ2L+LX,I)=CNA*FFM(C1*RESISM(I-1)/TM(I-1))
         CSUBM(KXB1+LX,KXJ2L+LX,I)=CNA*FFP(C1*RESISM(I)/TM(I))
         ENDIF

         IF (KXJRE2L.GT.0) THEN
         BSUBM(KXB1+LX,KXJRE2L+LX,I)=-CNA*RE_EQFRAC
     &        *FF(C1*RESIST(I)/T(I),C1*RESISM(I-1)/TM(I-1),
     &            C1*RESISM(I)/TM(I))
         ASUBM(KXB1+LX,KXJRE2L+LX,I)=-CNA*RE_EQFRAC
     &        *FFM(C1*RESISM(I-1)/TM(I-1))
         CSUBM(KXB1+LX,KXJRE2L+LX,I)=-CNA*RE_EQFRAC
     &        *FFP(C1*RESISM(I)/TM(I))
         ENDIF
      ENDIF

      BsubM(kxB1+lx,kxB1+lx,I)=BsubM(kxB1+lx,kxB1+lx,I)+
     &                         FF(SHIFTBC(I),SHIFTBM(I-1),SHIFTBM(I))
      AsubM(kxB1+lx,kxB1+lx,I)=AsubM(kxB1+lx,kxB1+lx,I)+
     &                         FFM(SHIFTBM(I-1))
      CsubM(kxB1+lx,kxB1+lx,I)=CsubM(kxB1+lx,kxB1+lx,I)+
     &                         FFP(SHIFTBM(I))

C     DIAMAGNETIC TERM
      IF (IDIAMB.EQ.1.OR.IDIAMB.EQ.3) THEN
         IF (ABS(RM(MSA,2)).GT.0.1) THEN
            BsubM(kxB1+lx, kxPD+lx,i)=BsubM(kxB1+lx, kxPD+lx,i) 
     $        +cma*FF(T(i)*C1,TM(i-1)*C1,TM(i)*C1)*FDIAMB
            AsubM(kxB1+lx, kxPD+lx,i)=AsubM(kxB1+lx, kxPD+lx,i) 
     &        +cma*FFM(TM(i-1)*C1)*FDIAMB
            CsubM(kxB1+lx, kxPD+lx,i)=CsubM(kxB1+lx, kxPD+lx,i)  
     &        +cma*FFP(TM(i)*C1)*FDIAMB

            IF (KXPED.GT.0) THEN
            BsubM(kxB1+lx, kxPED+lx,i)=BsubM(kxB1+lx, kxPED+lx,i) 
     $        +cma*FF(T(i)*C1,TM(i-1)*C1,TM(i)*C1)*FDIAMB
            AsubM(kxB1+lx, kxPED+lx,i)=AsubM(kxB1+lx, kxPED+lx,i) 
     &        +cma*FFM(TM(i-1)*C1)*FDIAMB
            CsubM(kxB1+lx, kxPED+lx,i)=CsubM(kxB1+lx, kxPED+lx,i)  
     &        +cma*FFP(TM(i)*C1)*FDIAMB
            ENDIF
         ENDIF
      ENDIF
      ENDDO

      DO I=1,NR
      INCLUDE 'tophat.inc'

      IF (IPDIVB.NE.2.OR.(IPDIVB.EQ.2.AND.ABS(RM(MSA,2)).LT.0.1)) THEN
      IF (KVSQL(8)) THEN
      FsubM(kyB2+ly, kxV1+lx,I)=znorm*DPsiDs(i)
      GsubM(kyB2+ly, kxV1+lx,I)=-znorm*DPsiDs(i+1)
      ENDIF

      DsubM(kyB2+ly,kyB2+ly,I)=DsubM(kyB2+ly,kyB2+ly,I)+
     &                         GG(SHIFTBM(I),SHIFTBC(I),SHIFTBC(I+1))

C     DIAMAGNETIC TERM
      IF (IDIAMB.EQ.1.OR.IDIAMB.EQ.3) THEN
         FsubM(kyB2+ly, kxPD+lx,i)=FsubM(kyB2+ly, kxPD+lx,i)  
     $    +znorm*TM(i)*FDIAMB
     $    -GF(DPSIDS(i)*TP(i)*C1,DPSIDSM(i)*TPM(i)*C1)*FDIAMB
         GsubM(kyB2+ly, kxPD+lx,i)=GsubM(kyB2+ly, kxPD+lx,i) 
     $    -znorm*TM(i)*FDIAMB
     $    -GF(DPSIDS(i+1)*TP(i+1)*C1,DPSIDSM(i)*TPM(i)*C1)*FDIAMB

         IF (KXPED.GT.0) THEN
         FsubM(kyB2+ly, kxPED+lx,i)=FsubM(kyB2+ly, kxPED+lx,i)  
     $    +znorm*TM(i)*FDIAMB
     $    -GF(DPSIDS(i)*TP(i)*C1,DPSIDSM(i)*TPM(i)*C1)*FDIAMB
         GsubM(kyB2+ly, kxPED+lx,i)=GsubM(kyB2+ly, kxPED+lx,i) 
     $    -znorm*TM(i)*FDIAMB
     $    -GF(DPSIDS(i+1)*TP(i+1)*C1,DPSIDSM(i)*TPM(i)*C1)*FDIAMB
         ENDIF
      ENDIF
      IF (IDIAMB.EQ.4) THEN
         CTMP1 = B0K/OMEGACI0/RHO(I)*DLNRHO(I)
         CTMP2 = B0K/OMEGACI0/RHO(I+1)*DLNRHO(I+1)
         CTMP3 = B0K/OMEGACI0/RHOM(I)*DLNRHOM(I)
         DsubM(kyB2+ly, kyPE+ly,i)=cna*GG(CTMP3,CTMP1,CTMP2)*FDIAMB
      ENDIF
      ENDIF

      IF (IPDIVB.NE.1) THEN 
      z3m  = (cs(i  )/csm(i))**iexb3
      z3p  = (cs(i+1)/csm(i))**iexb3
      IF (IDIAMB.EQ.4) THEN
         CTMP1 = B0K/OMEGACI0/RHO(I)*DLNRHO(I)
         CTMP2 = B0K/OMEGACI0/RHO(I+1)*DLNRHO(I+1)
         CTMP3 = B0K/OMEGACI0/RHOM(I)*DLNRHOM(I)
         DsubM(kyB3+ly, kyPE+ly,i)=-cma*GG(CTMP3,CTMP1,CTMP2)*FDIAMB
      ENDIF
      DsubM(kyB3+ly, kyB3+ly,I)=DsubM(kyB3+ly, kyB3+ly,I)+
     &   GG(SHIFTBM(I),SHIFTBC(I)*z3m,SHIFTBC(I+1)*z3p)
      ENDIF

      z1m  = (cs(i  )/csm(i))**iexj1
      z1p  = (cs(i+1)/csm(i))**iexj1
      DSUBM(KYJ1+LY,KYJ1+LY,I)=GG(C1,C1*z1m,C1*z1p) 

      IF (KJRER.EQ.6) 
     &DSUBM(KYJRE1+LY,KYJRE1+LY,I)=GG(C1,C1*Z1M,C1*Z1P) 
      ENDDO

C     DIV B = 0. AS THE EQUATION FOR B3
C     YQ LIU, 2010-01-05
      IF (IPDIVB.EQ.1) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB3+LY,KXB1+LX,I)=-ZNORM
      GSUBM(KYB3+LY,KXB1+LX,I)= ZNORM
      ELSE
      FSUBM(KYB3+LY,KXB1+LX,I)=-ZNORM*T(I)
      GSUBM(KYB3+LY,KXB1+LX,I)= ZNORM*T(I)
      ENDIF

      DSUBM(KYB3+LY,KYB2+LY ,I)= CMA * GG(C1,C1,C1)
      DSUBM(KYB3+LY,KYB3+LY ,I)= CNA * GG(C1,C1,C1)
      ENDDO
      ENDIF

      IF (IPDIVB.EQ.2.AND.ABS(RM(MSA,2)).GT.0.1) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'
      FSUBM(KYB2+LY,KXB1+LX ,I)= -ZNORM
      GSUBM(KYB2+LY,KXB1+LX ,I)=  ZNORM
      DSUBM(KYB2+LY,KYB2+LY ,I)= CMA * GG(C1,C1,C1)
      ENDDO
      ENDIF

      DO I=2,NR
      INCLUDE 'tent.inc'
      BSUBM(KXJ2U+LX,KXJ2U+LX,I)=FF(C1,C1,C1)
      ASUBM(KXJ2U+LX,KXJ2U+LX,I)=FFM(C1)
      CSUBM(KXJ2U+LX,KXJ2U+LX,I)=FFP(C1)
      ENDDO

      IF (KJRER.EQ.6) THEN
      DO I=2,NR
      INCLUDE 'tent.inc'
      BSUBM(KXJRE2+LX,KXJRE2+LX,I)=FF(C1,C1,C1)
      ASUBM(KXJRE2+LX,KXJRE2+LX,I)=FFM(C1)
      CSUBM(KXJRE2+LX,KXJRE2+LX,I)=FFP(C1)
      ENDDO
      ENDIF

C     BC FOR J2U and J3U 
      I = NRP1
      HM = CSH(NR)
      HP = CSH(NRP1)
      ZNORM = 2./(HM+VCSH(1))
      BSUBM(KXJ2U+LX,KXJ2U+LX,I)=FF0(C1,C1) +  FF1(C1,C1)
      ASUBM(KXJ2U+LX,KXJ2U+LX,I)=FFM(C1)
      CSUBM(KXJ2U+LX,KXJ2U+LX,I)=FFP(C1)

      IF (KJRER.EQ.6) THEN
      BSUBM(KXJRE2+LX,KXJRE2+LX,I)=FF0(C1,C1) +  FF1(C1,C1)
      ASUBM(KXJRE2+LX,KXJRE2+LX,I)=FFM(C1)
      CSUBM(KXJRE2+LX,KXJRE2+LX,I)=FFP(C1)
      ENDIF

C     FREE BOUNDARY CONDITION FOR J2U WITH N=0
      IF (ABS(RNTOR).LT.1.E-10.AND.
     &    (K_BC_N0.EQ.12.OR.K_BC_N0.EQ.22.OR.K_BC_N0.EQ.32))
     &   ASUBM(KXJ2U+LX,KXJ2U+LX,I)=-FF0(C1,C1)-FF1(C1,C1)

C     ROBIN BC MATCHING SLOPE OF EQUILIBRIUM CURRENT
      IF (ABS(RNTOR).LT.1.E-10.AND.ABS(RM(MSA,2)).LT.0.1.AND.
     &    (K_BC_N0.EQ.13.OR.K_BC_N0.EQ.23.OR.K_BC_N0.EQ.33)) THEN
         CTMP1 = TP(I-1)*DPSIDS(I-1) 
         CTMP2 = TP(I)*DPSIDS(I) 
         CTMP3 = (CTMP2-CTMP1)/CTMP2/CSH(I-1)
         CTMP1 =-1./(1.-CSH(I-1)*CTMP3)
         ASUBM(KXJ2U+LX,KXJ2U+LX,I)=CTMP1*(FF0(C1,C1)+FF1(C1,C1))
      ENDIF

      DO I=2,NR
      INCLUDE 'tent.inc'
      BSUBM(KXJ3+LX,KXJ3+LX,I)=FF(C1,C1,C1)
      ASUBM(KXJ3+LX,KXJ3+LX,I)=FFM(C1)
      CSUBM(KXJ3+LX,KXJ3+LX,I)=FFP(C1)
      ENDDO

      IF (KJRER.EQ.6) THEN
      DO I=2,NR
      INCLUDE 'tent.inc'
      BSUBM(KXJRE3+LX,KXJRE3+LX,I)=FF(C1,C1,C1)
      ASUBM(KXJRE3+LX,KXJRE3+LX,I)=FFM(C1)
      CSUBM(KXJRE3+LX,KXJRE3+LX,I)=FFP(C1)
      ENDDO
      ENDIF

      I = NRP1
      HM = CSH(NR)
      HP = CSH(NRP1)
      ZNORM = 2./(HM+VCSH(1))
      BSUBM(KXJ3+LX,KXJ3+LX,I)=FF0(C1,C1) + FF1(C1,C1)
      ASUBM(KXJ3+LX,KXJ3+LX,I)=FFM(C1)
      CSUBM(KXJ3+LX,KXJ3+LX,I)=FFP(C1)

      IF (KJRER.EQ.6) THEN
      BSUBM(KXJRE3+LX,KXJRE3+LX,I)=FF0(C1,C1) + FF1(C1,C1)
      ASUBM(KXJRE3+LX,KXJRE3+LX,I)=FFM(C1)
      CSUBM(KXJRE3+LX,KXJRE3+LX,I)=FFP(C1)
      ENDIF

C     FREE BOUNDARY CONDITION FOR J3U WITH N=0
      IF (ABS(RNTOR).LT.1.E-10.AND.
     &    (K_BC_N0.EQ.12.OR.K_BC_N0.EQ.22.OR.K_BC_N0.EQ.32))
     &   ASUBM(KXJ3+LX,KXJ3+LX,I)=-FF0(C1,C1)-FF1(C1,C1)

C     ROBIN BC MATCHING SLOPE OF EQUILIBRIUM CURRENT
      IF (ABS(RNTOR).LT.1.E-10.AND.ABS(IM(MSA,2)).LT.MEDIM.AND.
     &    (K_BC_N0.EQ.13.OR.K_BC_N0.EQ.23.OR.K_BC_N0.EQ.33)) THEN
         J = IM(MSA,2)   
         IF (J.GE.0) THEN
            CTMP1 = J3B2(I-1,J+1)
            CTMP2 = J3B2(I,J+1)
         ELSE
            CTMP1 = CONJG(J3B2(I-1,-J+1))
            CTMP2 = CONJG(J3B2(I,-J+1))
         ENDIF
         CTMP3 = (CTMP2-CTMP1)/CTMP2/CSH(I-1)
         CTMP1 =-1./(1.-CSH(I-1)*CTMP3)
         ASUBM(KXJ3+LX,KXJ3+LX,I)=CTMP1*(FF0(C1,C1)+FF1(C1,C1))
         WRITE(*,*) 'K_BC_N0 J3U J CTMP1=',J,CTMP1
      ENDIF

      zos = ALPHAP
      IF (KYPE.LE.0) zos = 1.
      DO I=1,NR
      INCLUDE 'tophat.inc'
      CTMP1 = zos*Peq(i)*DPsiDs(i)
      CTMP2 = zos*Peq(i+1)*DPsiDs(i+1)
      CTMP3 = zos*PeqM(i)*DPsiDsM(i)
      IF (V2XKEY.EQ.0.OR.V2XKEY.EQ.1) THEN
      IF (KYV3.GT.0)
     &DsubM(kyPr+ly, kyV3+ly,I)=DsubM(kyPr+ly, kyV3+ly,I)
     &  -CMA*gamarm(i)*GG(CTMP3,CTMP1,CTMP2) 
      IF (KYW3.GT.0) THEN
      CTMP1 = TTCPARAM(I) + TTCPERPM(I)
      CTMP2 = TTCPARAI(I) + TTCPERPI(I)
      CTMP3 = TTCPARAI(I+1) + TTCPERPI(I+1)
      DSUBM(KYPR+LY,KYW3+LY,I)=DSUBM(KYPR+LY,KYW3+LY,I)+CMA*
     &       GG(C1*DPsiDsM(I)*CTMP1,C1*DPsiDs(I)*CTMP2,
     &          C1*DPsiDs(I+1)*CTMP3)
      ENDIF
      ELSE
      DsubM(kyPr+ly, kyX3+ly,I)=DsubM(kyPr+ly, kyX3+ly,I)
     &  -CMA*gamarm(i)*GG(CTMP3,CTMP1,CTMP2)      
      ENDIF
      ENDDO

      IF (KYRHOP.GT.0.AND.KYV3.GT.0) THEN 
      DO I=1,NR
      INCLUDE 'tophat.inc'
      CTMP1 = RHO(i)*DPsiDs(i)
      CTMP2 = RHO(i+1)*DPsiDs(i+1)
      CTMP3 = RHOM(i)*DPsiDsM(i)
      DSUBM(KYRHOP+LY, KYV3+LY,I)=DSUBM(KYRHOP+LY, KYV3+LY,I)
     &  -CMA*GG(CTMP3,CTMP1,CTMP2) 
      ENDDO
      ENDIF
      
      IF (KXX1.GT.0) THEN
      DO I=2,NRP1
      INCLUDE 'tent.inc'
      BsubM(kxx1+lx,kxx1+lx,i)=BsubM(kxx1+lx,kxx1+lx,i)+
     &                         FF(SHIFTC(I),SHIFTM(I-1),SHIFTM(I))
      AsubM(kxx1+lx,kxx1+lx,i)=AsubM(kxx1+lx,kxx1+lx,i)+
     &                         SHIFTM(I-1)*FFM(C1)
      CsubM(kxx1+lx,kxx1+lx,i)=CsubM(kxx1+lx,kxx1+lx,i)+
     &                         SHIFTM(I)*FFP(C1)
      BsubM(kxx1+lx,kxv1+lx,i)=FF(C1,C1,C1)
      AsubM(kxx1+lx,kxv1+lx,i)=FFM(C1)
      CsubM(kxx1+lx,kxv1+lx,i)=FFP(C1)
      ENDDO
      ENDIF

      IF (IDIAMTE.EQ.2.AND.KYPE.GT.0) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'
      CTMP1 = (1.-ALPHAP)*Peq(i)*DPsiDs(i)
      CTMP2 = (1.-ALPHAP)*Peq(i+1)*DPsiDs(i+1)
      CTMP3 = (1.-ALPHAP)*PeqM(i)*DPsiDsM(i)
      IF (V2XKEY.EQ.0.OR.V2XKEY.EQ.1) THEN
      IF (KYV3.GT.0) 
     &DsubM(kyPe+ly, kyV3+ly,I)=DsubM(kyPe+ly, kyV3+ly,I)
     &  -CMA*gamarm(i)*GG(CTMP3,CTMP1,CTMP2)
      ELSE
      DsubM(kyPe+ly, kyX3+ly,I)=DsubM(kyPe+ly, kyX3+ly,I)
     &  -CMA*gamarm(i)*GG(CTMP3,CTMP1,CTMP2)      
      ENDIF 
      ENDDO
      ENDIF

      IF (IDIAMTE.EQ.1.AND.KYPE.GT.0) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'
      CTMP1 = (1.-ALPHAP)*(PPeq(i)*DPsiDs(i)-Peq(i)*DLNRHO(i))
      CTMP2 = (1.-ALPHAP)*(PPeq(i+1)*DPsiDs(i+1)-Peq(i+1)*DLNRHO(i+1))
      CTMP3 = (1.-ALPHAP)*(PPeqM(i)*DPsiDsM(i)-PeqM(i)*DLNRHOM(i))
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
        FSUBM(KyPe+LY,KxB1+LX,I)=GF(CTMP1,CTMP3)
        GSUBM(KyPe+LY,KxB1+LX,I)=GF(CTMP2,CTMP3)
      ELSE
        FSUBM(KyPe+LY,KxB1+LX,I)=GF(CTMP1*T(i),CTMP3*TM(i))
        GSUBM(KyPe+LY,KxB1+LX,I)=GF(CTMP2*T(i+1),CTMP3*TM(i))
      ENDIF
      DsubM(kyPe+ly, kyPe+ly,I)=DsubM(kyPe+ly, kyPe+ly,I)
     &  +CMA*GG(C1*DpsiDsM(i),C1*DpsiDs(i),C1*DpsiDs(i+1))
      ENDDO
      ENDIF

      IEXE = -1
      IF (KXX1.GT.0.AND.KYX2.GT.0) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      CTMP1 = -DpsiDs(I)*DROT(I)
      CTMP2 = -DpsiDs(I+1)*DROT(I+1)
      CTMP3 = -DpsiDsM(I)*DROTM(I)
      FSUBM(KYX2+LY,KXX1+LX,I)=FSUBM(KYX2+LY,KXX1+LX,I)
     &     + GF(ZEM*CTMP1,CTMP3)
      GSUBM(KYX2+LY,KXX1+LX,I)=GSUBM(KYX2+LY,KXX1+LX,I)
     &     + GF(ZEP*CTMP2,CTMP3)
      ENDDO
      ENDIF
C     GLX------<
      IF (NPROFRP.GT.0) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      CTMP1 = RHOUM(I)*DPsiDsM(I)
      CTMP2 = RHOU(I)*DPsiDs(I)
      CTMP3 = RHOU(I+1)*DPsiDs(I+1)
      IF (kyX2.GT.0) THEN
      DsubM(kyX2+ly, kyX2+ly,I)=DsubM(kyX2+ly, kyX2+ly,I)
     &   - CMA*GG(CTMP1,ZEM*ZV2M*CTMP2,ZEP*ZV2P*CTMP3)
      ENDIF
      IF (kyX3.GT.0) THEN
      DsubM(kyX3+ly, kyX3+ly,I)=DsubM(kyX3+ly, kyX3+ly,I)
     &   - CMA*GG(CTMP1,ZEM*CTMP2,ZEP*CTMP3)
      ENDIF
      ENDDO
      ENDIF
C     --------->

      IF (IPDIVB.NE.1) THEN
      DO 191 I = 1,NR
      INCLUDE 'tophat.inc'
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB3+LY ,KXB1+LX,I)= FSUBM(KYB3+LY ,KXB1+LX,I)+
     &                               GF((DROT(I)+DOMEGASE(I)*zobe-
     &                               DOMEGASI(I)*zobi)*C1,
     &                              (DROTM(I)+DOMEGASEM(I)*zobe-
     &                               DOMEGASIM(I)*zobi)*C1)
      GSUBM(KYB3+LY ,KXB1+LX,I)= GSUBM(KYB3+LY ,KXB1+LX,I)+
     &                               GF((DROT(I+1)+DOMEGASE(I+1)*zobe-
     &                               DOMEGASI(I+1)*zobi)*C1,
     &                              (DROTM(I)+DOMEGASEM(I)*zobe-
     &                               DOMEGASIM(I)*zobi)*C1)
      ELSE
      FSUBM(KYB3+LY,KXB1+LX,I)=FSUBM(KYB3+LY,KXB1+LX,I)+
     &                             GF((DROT(I)+DOMEGASE(I)*zobe-
     &                             DOMEGASI(I)*zobi)*C1*T(I),
     &                            (DROTM(I)+DOMEGASEM(I)*zobe-
     &                             DOMEGASIM(I)*zobi)*C1*TM(I))
      GSUBM(KYB3+LY,KXB1+LX,I)=GSUBM(KYB3+LY,KXB1+LX,I)+
     &                             GF((DROT(I+1)+DOMEGASE(I+1)*zobe-
     &                             DOMEGASI(I+1)*zobi)*C1*T(I+1),
     &                            (DROTM(I)+DOMEGASEM(I)*zobe-
     &                             DOMEGASIM(I)*zobi)*C1*TM(I))
      ENDIF

      IF (KXJ2L.GT.0) THEN
      FSUBM(KYB3+LY ,KXJ2L+LX,I)=znorm*resist(i)
      GSUBM(KYB3+LY ,KXJ2L+LX,I)=-znorm*resist(i+1)
      ENDIF

      IF (KXJRE2L.GT.0) THEN
      FSUBM(KYB3+LY ,KXJRE2L+LX,I)=-znorm*JRE_EQFRAC*resist(i)
      GSUBM(KYB3+LY ,KXJRE2L+LX,I)= znorm*JRE_EQFRAC*resist(i+1)
      ENDIF
 191  CONTINUE
      ENDIF
   
      IF (V2U_M0.AND.ABS(RM(MSA,2)).LT.1.E-3) THEN
         DO I=1,NR
         CALL ANNIHY(KYV2,MSA,I,
     &        ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
         DSUBM(KYV2+LY,KYV2+LY,I) = 1.0
         ENDDO
      ENDIF
         
 200  CONTINUE

      WRITE (*,'("CALLING KJP")')
C
      IF (INCKIN.EQ.1) 
     &   CALL KJP(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM,
     &            SHIFTC,SHIFTM,MXMAX,MYMAX,NRP1)

C     MULTIPLY EQUATIONS INSIDE PLASMA BY EQFAC
C     YQ LIU, 2009-10-01
      DO MSA = 1,MSMAX
         LX=(MSA-1)*NXCOMP
         DO LY=1,MXMAX
            DO I=1,NR
               IF (KXX1.GT.0) THEN
               ASUBM(LX+KXX1,LY,I)=ASUBM(LX+KXX1,LY,I)*EQFAC/EQFAC
               BSUBM(LX+KXX1,LY,I)=BSUBM(LX+KXX1,LY,I)*EQFAC/EQFAC
               CSUBM(LX+KXX1,LY,I)=CSUBM(LX+KXX1,LY,I)*EQFAC/EQFAC
               ENDIF

               ASUBM(LX+KXV1,LY,I)=ASUBM(LX+KXV1,LY,I)*EQFAC/EQFAC
               BSUBM(LX+KXV1,LY,I)=BSUBM(LX+KXV1,LY,I)*EQFAC/EQFAC
               CSUBM(LX+KXV1,LY,I)=CSUBM(LX+KXV1,LY,I)*EQFAC/EQFAC

               ASUBM(LX+KXB1,LY,I)=ASUBM(LX+KXB1,LY,I)*EQFAC
               BSUBM(LX+KXB1,LY,I)=BSUBM(LX+KXB1,LY,I)*EQFAC
               CSUBM(LX+KXB1,LY,I)=CSUBM(LX+KXB1,LY,I)*EQFAC
            ENDDO
         ENDDO
         DO LY=1,MYMAX
            DO I=1,NR
               IF (KXX1.GT.0) THEN
               HSUBM(LX+KXX1,LY,I)=HSUBM(LX+KXX1,LY,I)*EQFAC/EQFAC
               ESUBM(LX+KXX1,LY,I)=ESUBM(LX+KXX1,LY,I)*EQFAC/EQFAC
               ENDIF

               HSUBM(LX+KXV1,LY,I)=HSUBM(LX+KXV1,LY,I)*EQFAC/EQFAC
               ESUBM(LX+KXV1,LY,I)=ESUBM(LX+KXV1,LY,I)*EQFAC/EQFAC

               HSUBM(LX+KXB1,LY,I)=HSUBM(LX+KXB1,LY,I)*EQFAC
               ESUBM(LX+KXB1,LY,I)=ESUBM(LX+KXB1,LY,I)*EQFAC
            ENDDO
         ENDDO

         LX=(MSA-1)*NYCOMP
         DO LY=1,MXMAX
            DO I=1,NR
               IF (KYX2.GT.0) THEN
               FSUBM(LX+KYX2,LY,I)=FSUBM(LX+KYX2,LY,I)*EQFAC/EQFAC
               GSUBM(LX+KYX2,LY,I)=GSUBM(LX+KYX2,LY,I)*EQFAC/EQFAC
               ENDIF

               IF (KYX3.GT.0) THEN
               FSUBM(LX+KYX3,LY,I)=FSUBM(LX+KYX3,LY,I)*EQFAC/EQFAC
               GSUBM(LX+KYX3,LY,I)=GSUBM(LX+KYX3,LY,I)*EQFAC/EQFAC
               ENDIF

               FSUBM(LX+KYV2,LY,I)=FSUBM(LX+KYV2,LY,I)*EQFAC/EQFAC
               GSUBM(LX+KYV2,LY,I)=GSUBM(LX+KYV2,LY,I)*EQFAC/EQFAC

               IF (KYV3.GT.0) THEN
               FSUBM(LX+KYV3,LY,I)=FSUBM(LX+KYV3,LY,I)*EQFAC/EQFAC
               GSUBM(LX+KYV3,LY,I)=GSUBM(LX+KYV3,LY,I)*EQFAC/EQFAC
               ENDIF

               FSUBM(LX+KYB2,LY,I)=FSUBM(LX+KYB2,LY,I)*EQFAC/EQFAC
               GSUBM(LX+KYB2,LY,I)=GSUBM(LX+KYB2,LY,I)*EQFAC/EQFAC

               FSUBM(LX+KYB3,LY,I)=FSUBM(LX+KYB3,LY,I)*EQFAC/EQFAC
               GSUBM(LX+KYB3,LY,I)=GSUBM(LX+KYB3,LY,I)*EQFAC/EQFAC
            ENDDO
         ENDDO
         DO LY=1,MYMAX
            DO I=1,NR
               IF (KYX2.GT.0) 
     &         DSUBM(LX+KYX2,LY,I)=DSUBM(LX+KYX2,LY,I)*EQFAC/EQFAC
               IF (KYX3.GT.0)
     &         DSUBM(LX+KYX3,LY,I)=DSUBM(LX+KYX3,LY,I)*EQFAC/EQFAC
               DSUBM(LX+KYV2,LY,I)=DSUBM(LX+KYV2,LY,I)*EQFAC/EQFAC
               IF (KYV3.GT.0)
     &         DSUBM(LX+KYV3,LY,I)=DSUBM(LX+KYV3,LY,I)*EQFAC/EQFAC
               DSUBM(LX+KYB2,LY,I)=DSUBM(LX+KYB2,LY,I)*EQFAC/EQFAC
               DSUBM(LX+KYB3,LY,I)=DSUBM(LX+KYB3,LY,I)*EQFAC/EQFAC
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(SHIFTC,SHIFTM,SHIFTVC,SHIFTVM,SHIFTBC,SHIFTBM)
C
      RETURN
      END

*DECK PATCH
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--------PATCH FOR GENERALISED FEM ------ A.B. 11.07.90 ----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE PATCH
C     ================
C
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE GIJLM
C
      INTEGER I,MS
C
C     DO 100 MS = 1,MEDIM
CMSC      DPEDS(NRP1+1,MS)=0.
C
C
C100  CONTINUE
C
      RETURN
      END
*DECK GCONTR
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------COMPUTE J*CONTRAV. COMP OF GIJ----------- A.B. 28.11.90 --------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE GCONTR(NFMAX)
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE MPIENV
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER    NFMAX
      REAL*8,DIMENSION(:,:),ALLOCATABLE::R11,R12,R22,R33,WORK1
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::WORK3
      INTEGER    I,J,NFOU2,NFOU
      REAL*8       ZDET,ZSAVE,SUMJ,SUME
C
C     GET MAXIMUM POWER OF 2 ALLOWED BY AVAILABLE SPACE FOR FFT
C
C
      ALLOCATE(R11(NRP1,NFMAX),R12(NRP1,NFMAX),R22(NRP1,NFMAX),
     &         R33(NRP1,NFMAX),WORK1(NRP1,NFMAX),
     &         WORK3(NRP1,NFMAX/2))
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'GCONTR' 
C
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,'("NFMAX,MSDIM= ",2I7)')NFMAX,MSDIM
      ENDIF
C
      NFOU2 = 1
      NFOU = NFMAX
      IF (NFOU.GT.40*MEDIM) NFOU=40*MEDIM
 10   IF (NFOU2*2.GT.NFOU) GOTO 20
      NFOU2 = 2*NFOU2
      GOTO 10
 20   CONTINUE
 
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*) ' GCONTR: # REAL*8 SPACE POINTS FOR FFT = ',NFOU2
      ENDIF
C
C.....TRANSFORM TO REAL*8 SPACE D (=1/J) * GIJ (COVARIANT)
C
      NPSTRT    =  1
      call FFTDRIVER(R11,   DG11L,  BCKWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NFOU2, KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG11L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(R12,   DG12L,  BCKWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NFOU2, KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG12L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(R22,   DG22L,  BCKWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NFOU2, KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG22L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(R33,   DG33L,  BCKWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NFOU2, KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG33L in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(R11,    DG11L,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NFOU2,   KUFFTP, 'R11')
        call FFTOUTPT(R12,    DG12L,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NFOU2,   KUFFTP, 'R12')
        call FFTOUTPT(R22,    DG22L,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NFOU2,   KUFFTP, 'R22')
        call FFTOUTPT(R33,    DG33L,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NFOU2,   KUFFTP, 'R33')
      ENDIF
C
C.....INVERT METRIC TENOR IN CONFIGURATION SPACE
C
      DO 50 I=1,NRP1
      DO 50 J=1,NFOU2
C
      ZDET = 1./(R11(I,J)*R22(I,J) - R12(I,J)*R12(I,J))
      WORK1(I,J) = ZDET/R33(I,J)
C
 50   CONTINUE
C
C     TEST JACOBIAN
C
      NPSTRT    =  1
      call FFTDRIVER(WORK1, WORK3,  FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NFOU2, KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: WORK1 in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      9,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(WORK1,  WORK3,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NFOU2,   KUFFTP, 'WORK3')
      ENDIF
C
      SUMJ=0.
      SUME=0.
      DO 60 J=1,MEDIM
      DO 60 I=1,NRP1
      WORK3(I,J)=WORK3(I,J) - JACOBI(I,J)
      SUMJ=SUMJ + ABS(JACOBI(I,J))
      SUME=SUME + ABS(WORK3(I,J))
 60   CONTINUE
      SUME = SUME/SUMJ
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*) ' CHECK OF INVERSION: RELATIVE ERROR = ',SUME
      ENDIF

      DEALLOCATE(R11,R12,R22,R33)
      DEALLOCATE(WORK1,WORK3)
C
      RETURN
      END
*DECK BOUNDC
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------BOUNDARY CONDITIONS--------A. BONDESON 14.02.91-----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BOUNDC(
C     ==================
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
C
      USE DIMENSIM
      USE GLOBALM
      INCLUDE 'specmat.inc'

C
      INTEGER    NS
      PARAMETER (NS=2)
      INTEGER    MSROW,LXROW,LYROW,ARM,I
      INTEGER    KX,KY
C
      DO 100 MSROW = 1,MSMAX
      ARM = ABS(RM(MSROW,NS))
      IF (ARM.EQ.0) ARM = 2
      LXROW = (MSROW-1)*NXCOMP
      LYROW = (MSROW-1)*NYCOMP

      I = 1
      DO 5 KX = 1,NXCOMP
      CALL ANNIHX(KX,MSROW,I,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)

      BSUBM(KX+LXROW,KX+LXROW,I) = 1.
 5    CONTINUE
      IF (ARM.EQ.1.) CSUBM(KXV1+LXROW,KXV1 +LXROW,I) =-1.

C     FREE BC FOR N=0,M=0 J2U&J3U
      IF (ABS(RNTOR).LT.1.E-10.AND.ARM.LT.1.E-10.AND.
     &    (K_BC_N0.GE.20.AND.K_BC_N0.LE.23)) THEN
         CSUBM(KXJ2U+LXROW,KXJ2U +LXROW,I) =-1. 
         CSUBM(KXJ3+LXROW,KXJ3 +LXROW,I) =-1. 
         IF (KXJ2L.GT.0) CSUBM(KXJ2L+LXROW,KXJ2L +LXROW,I) =-1. 
      ENDIF

      I = 1
      DO 6 KY = 1,NYCOMP
C     CALL ANNIHY(KY,MSROW,I,
C    &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)

      DSUBM(KY+LYROW,KY+LYROW,I) = 1.
 6    CONTINUE

      DO 10 I=2,NFIT
      CALL ANNIHX(KXV1,MSROW,I,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)

      BSUBM(KXV1+LXROW,KXV1 +LXROW,I) = CS (I  )**(-ARM+1)
      CSUBM(KXV1+LXROW,KXV1 +LXROW,I) =-CS (I+1)**(-ARM+1)
 10   CONTINUE

 100  CONTINUE

      RETURN
      END
C
*DECK MUBMAT
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------INERTIAL OPERATOR----------A. BONDESON 14.02.91-----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE MUBMAT
C     =================
     $            (MD,MDY,ND,NXC,NYC,NCASE,NITMAX
     $            ,EPSPAM,EPSDET
     $            ,AL0,ALAM,ALNORM,NONCON
     $            ,A, B, C, D, E, F, G, H
     $            ,DX,DY, X, Y
     $            ,DIX,DIY,R,RY
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2
     $            ,XOLD,YOLD,XPIVOT,YPIVOT
     $            ,TMP1,TMP2,TMPY1,TMPY2
     $                                                    )
C
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE CONVOLCOFM
      USE GVACUUMM
      USE RCOMDM
      USE FEEDBACKM
      USE REORBITM
      INCLUDE 'cdipam.inc'
      INCLUDE 'comioc.inc'
C
      INTEGER       I,J,MS,J1,J2
      INTEGER       MSA,MSB,MSMI,MSPL,NSA,NSB
      INTEGER       LXROW,LYROW,LXCOL,LYCOL
      PARAMETER     (NSA=2,NSB=1)
      INTEGER       IEXV2,IEXJ1,IEXB3,IEXE
      PARAMETER     (IEXV2=-1, IEXJ1=1, IEXB3=1, IEXE=-1)
      REAL*8        ZEM,ZEP,ZV2M,ZV2P,z3m,z3p,AOMEGA1,
     &              RTMP,RTMP1,RTMP2,RTMP3,RTMP4,ZKPAVT
      COMPLEX*16    CTMP,CTMP1,CTMP2,CTMP3,CTMP4,CTMP5,CTMP6,
     &              TTMP,TTMP1,TTMP2,TTMP3
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::BNORM

      REAL*8        RCHI,HCHI,HH1,HH2,HH3,HH4,
     &              REB,REEPARA,RE_GAMMA,RE_Z
      COMPLEX*16    HB1, HB2, HB3, HJ1, HJ2, HJ3,
     &              HB1I,HB2I,HB3I,HJ1I,HJ2I,HJ3I,
     &              HB1M,HB2M,HB3M,HJ1M,HJ2M,HJ3M
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::RW1,RW2,RW3,RW4,
     &                                       RW5,RW6,RW7
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::RB1U,RB2U,RB3U,
     &                                       RJ1U,RJ2U,RJ3U,
     &                                       RV1U,RV2U,RV3U,
     &                                       RRHOE,RPPE
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::FJRERHS,FVSQVGP,FVSQPDV,
     &                                       FVSQRHOV
      COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE::FJRECURL,FVSQJXB,
     &                                         FVSQVGV,FVSQVXB,FVSQJEQ
C
      INCLUDE 'integc.inc'
C
      IF (NCASE.EQ.3.OR.NCASE.EQ.5.OR.NCASE.EQ.9) THEN
         IF (NONCON.EQ.-1) ICALPHA2=0
         IF (NONCON.NE.-1) ICALPHA2=1
      ENDIF

      TTMP = -(CALPHA2-CALPHA3)*CALPHA4/CALPHA2**2
      TTMP1= CONJG(TTMP)

      DO 10 MS=1,MD
      DO 10  I=1,ND+1
 10   R(MS,I) = 0.
      DO 20 MS=1,MDY
      DO 20  I=1,ND
 20   RY(MS,I) = 0.
C
      DO 90 MSA=1,MSMAX
      DO 90 MSB=1,MSMAX
      MSPL =  MPLUS(MSA,NSA,MSB,NSB)
      MSMI = MMINUS(MSA,NSA,MSB,NSB)
      IF (MSPL.LT.1) GOTO 60
C
      LXROW = (MSPL-1)*NXCOMP
      LXCOL = (MSA -1)*NXCOMP
      LYROW = (MSPL-1)*NYCOMP
      LYCOL = (MSA -1)*NYCOMP

      DO 50 I=2,NRP1
C
      TTMP2= -(CALPHA1+(CALPHA2-CALPHA3)*CI*RNTOR*TROTI(I))*
     &        CALPHA4/CALPHA2**2
      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      ZV2P = (CS(I)/CSM(I  ))**IEXV2
      INCLUDE 'tent.inc'
C
      R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP2*(
     &          FF(RGV1G11(i,msb),RGV1G11M(i-1,msb),
     &                  RGV1G11M(i,msb))*X(KXV1+LXCOL,I)
     &        + FFM(RGV1G11M(i-1,msb))*X(KXV1+LXCOL,I-1)
     &        + FFP(RGV1G11M(i,msb))*X(KXV1+LXCOL,I+1)
     &+FGM(RGV1G12(i,msb)*zv2m,RGV1G12M(i-1,msb))
     $     *Y(KYV2+LYCOL,I-1)
     &+FGP(RGV1G12(i,msb)*zv2p,RGV1G12M(i,msb))
     $     *Y(KYV2+LYCOL,I))
 
      IF ((CALPHA2-CALPHA3).NE.0..AND.CALPHA7.GT.0.) THEN
      Z3M  = (CSM(I-1)/CS(I))**IEXB3
      Z3P  = (CSM(I)  /CS(I))**IEXB3
 
      IF (INERT2) THEN
      CTMP1 = IDRXX(I,MSB)*TDROTI(I)
      CTMP2 = IDRXXM(I-1,MSB)*TDROTM(I-1)
      CTMP3 = IDRXXM(I,MSB)*TDROTM(I)
      R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP*(
     &                 FF(CTMP1,CTMP2,CTMP3)*X(KXV1+LXCOL,I)
     &                +FFM(CTMP2)*X(KXV1+LXCOL,I-1)
     &                +FFP(CTMP3)*X(KXV1+LXCOL,I+1))
      ENDIF
 
      IF (INERT1) THEN
      CTMP1 = IRXY(I,MSB)*TROTI(I)
      CTMP2 = IRXYM(I-1,MSB)*TROTM(I-1)
      CTMP3 = IRXYM(I,MSB)*TROTM(I)
      R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP*(
     &                   -FGM(CTMP1*ZV2M,CTMP2)*Y(KYV2+LYCOL,I-1)
     &                   -FGP(CTMP1*ZV2P,CTMP3)*Y(KYV2+LYCOL,I))
 
      CTMP1 = IRXZ(I,MSB)*TROTI(I)
      CTMP2 = IRXZM(I-1,MSB)*TROTM(I-1)
      CTMP3 = IRXZM(I,MSB)*TROTM(I)
      IF (KYV3.GT.0)
     &R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP*(
     &                   FGM(CTMP1,CTMP2)*Y(KYV3+LYCOL,I-1)
     &                 + FGP(CTMP1,CTMP3)*Y(KYV3+LYCOL,I))
      ENDIF
      ENDIF
 
 50   CONTINUE
 
      DO 55 I=1,NR
 
      TTMP2= -(CALPHA1+(CALPHA2-CALPHA3)*CI*RNTOR*TROTM(I))*
     &        CALPHA4/CALPHA2**2
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      INCLUDE 'tophat.inc'
 
      RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) + TTMP2*(
     &   GG(RGV2G22M(i,msb),RGV2G22(i,msb)*zem*zv2m,
     $        RGV2G22(i+1,msb)*zep*zv2p)
     &                             *Y(KYV2+LYCOL,I)
     & + GF(RGV1G12(i,msb)*zem,RGV1G12M(i,msb))
     &                             *X(KXV1+LXCOL,I)
     & + GF(RGV1G12(i+1,msb)*zep,RGV1G12M(i,msb))
     &                             *X(KXV1+LXCOL,I+1))
      IF (KYV3.GT.0)
     &RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP2*(
     &    GG(RGV3G33M(I,MSB),ZEM*RGV3G33(I,MSB),ZEP*RGV3G33(I+1,MSB))
     &                             *Y(KYV3+LYCOL,I))
      IF (KYX2.GT.0)
     &RY(KYX2+LYROW,I) = RY(KYX2+LYROW,I) + TTMP2*(
     &   GG(JACOBM(i,msb),JACOBI(i,msb)*zem*zv2m,
     &        JACOBI(i+1,msb)*zep*zv2p)
     &                             *Y(KYX2+LYCOL,I))
      IF (KYX3.GT.0)
     &RY(KYX3+LYROW,I) = RY(KYX3+LYROW,I) + TTMP2*(
     &   GG(JACOBM(i,msb),JACOBI(i,msb)*zem*zv2m,
     &        JACOBI(i+1,msb)*zep*zv2p)
     &                             *Y(KYX3+LYCOL,I))

      IF ((CALPHA2-CALPHA3).NE.0..AND.CALPHA7.GT.0..AND.
     &    KXX1.GT.0.AND.KYX3.GT.0) THEN
         RY(KYX3+LYROW,I) = RY(KYX3+LYROW,I) + TTMP*(
     &      GF(ZEM*TDROTI(I)*TB2(I,MSB),TDROTM(I)*TB2M(I,MSB))*
     &      X(KXX1+LXCOL,I) 
     &    + GF(ZEP*TDROTI(I+1)*TB2(I+1,MSB),TDROTM(I)*TB2M(I,MSB))*
     &      X(KXX1+LXCOL,I+1)) 
      ENDIF

      IF (V2XKEY.EQ.0.OR.V2XKEY.EQ.1) THEN
         RY(KYPR + LYROW, I) = RY(KYPR + LYROW, I) + TTMP2*TTCINERT0*(
     &       GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) 
     &                                *Y(KYPR+LYCOL,I))
         IF (IDIAMTE.EQ.2.AND.KYPE.GT.0) THEN
         RY(KYPE + LYROW, I) = RY(KYPE + LYROW, I) + TTMP2*(
     &       GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) 
     &                                *Y(KYPE+LYCOL,I))
         ENDIF
      ENDIF

      IF (KYRHOP.GT.0)  
     &RY(KYRHOP + LYROW, I) = RY(KYRHOP + LYROW, I) + TTMP2*(
     &    GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) 
     &                             *Y(KYRHOP+LYCOL,I)) 


      
      IF (IPERTURB.EQ.0.AND.INCKIN.EQ.1) THEN
         IF (V2XKEY.EQ.0.OR.V2XKEY.EQ.2) THEN
            RY(KYPPARA+LYROW,I) = RY(KYPPARA+LYROW,I) + TTMP2*(
     &      GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) 
     &     *Y(KYPPARA+LYCOL,I))
 
            RY(KYPPERP+LYROW,I) = RY(KYPPERP+LYROW,I) + TTMP2*(
     &      GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) 
     &     *Y(KYPPERP+LYCOL,I))

            IF (KYPE.GT.0) 
     &      RY(KYPE+LYROW,I) = RY(KYPE+LYROW,I) + TTMP2*(
     &      GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) 
     &     *Y(KYPE+LYCOL,I))
 
            IF (KYPP.GT.0) 
     &      RY(KYPP+LYROW,I) = RY(KYPP+LYROW,I) + TTMP2*(
     &      GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) 
     &     *Y(KYPP+LYCOL,I))
         ENDIF
      ENDIF

      IF ((CALPHA2-CALPHA3).NE.0..AND.CALPHA7.GT.0.) THEN
      IF (INERT2) THEN
      CTMP4 = IDRYX(I,MSB)*TDROTI(I)
      CTMP5 = IDRYX(I+1,MSB)*TDROTI(I+1)
      CTMP6 = IDRYXM(I,MSB)*TDROTM(I)
      RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) + TTMP*(
     &                   GF(ZEM*CTMP4,CTMP6)*X(KXV1+LXCOL,I)
     &                  +GF(ZEP*CTMP5,CTMP6)*X(KXV1+LXCOL,I+1))
      ENDIF
      IF (INERT1) THEN
      CTMP1 = IRXY(I,MSB)*TROTI(I)
      CTMP2 = IRXY(I+1,MSB)*TROTI(I+1)
      CTMP3 = IRXYM(I,MSB)*TROTM(I)
      CTMP4 = IRYZ(I,MSB)*TROTI(I)
      CTMP5 = IRYZ(I+1,MSB)*TROTI(I+1)
      CTMP6 = IRYZM(I,MSB)*TROTM(I)
      RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) + TTMP*(
     &                   GF(ZEM*CTMP1,CTMP3)*X(KXV1+LXCOL,I)
     &                  +GF(ZEP*CTMP2,CTMP3)*X(KXV1+LXCOL,I+1))
      IF (KYV3.GT.0)
     &RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) - TTMP*
     &         GG(CTMP6,ZEM*CTMP4,ZEP*CTMP5)*Y(KYV3+LYCOL,I)
      ENDIF

      IF (ABS(PVISC).GT.0..AND.IVISC.GT.0.AND.KYV3.GT.0) THEN
      ZKPAVT = ABS(RM(MSA,NSA)/QM(I) + RN(NSA))
     &                 * SQRT(REAL(PEQM(I))/RHOM(I))
      CTMP4 = IRYZ(I,MSB)*TROTI(I)
      CTMP5 = IRYZ(I+1,MSB)*TROTI(I+1)
      CTMP6 = IRYZM(I,MSB)*TROTM(I)
      IF (KYX2.GT.0) 
     &RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP*PVISC*ZKPAVT*0.5*(
     &    GG(CTMP6,ZV2M*ZEM*CTMP4,ZV2P*ZEP*CTMP5)*Y(KYX2+LYCOL,I))
 
      IF (KXX1.GT.0) THEN
      CTMP1 = -RHO(I)*TDROTI(I)*T(I)
      CTMP2 = -RHO(I+1)*TDROTI(I+1)*T(I+1)
      CTMP3 = -RHOM(I)*TDROTM(I)*TM(I)
      CTMP4 = IRXZ(I,MSB)*TROTI(I) 
      CTMP5 = IRXZ(I+1,MSB)*TROTI(I+1) 
      CTMP6 = IRXZM(I,MSB)*TROTM(I) 
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP*PVISC*ZKPAVT*(
     &     (-GF(ZEM*CTMP4,CTMP6)*0.5+GF(ZEM*CTMP1*JACOBI(I,MSB),
     &      CTMP3*JACOBM(I,MSB)))*X(KXX1+LXCOL,I)
     &   + (-GF(ZEP*CTMP5,CTMP6)*0.5+GF(ZEP*CTMP2*JACOBI(I+1,MSB),
     &      CTMP3*JACOBM(I,MSB)))*X(KXX1+LXCOL,I+1))
      ENDIF
      ENDIF
          
      IF (INERT1) THEN
      CTMP4 = IRYZ(I,MSB)*TROTI(I) 
      CTMP5 = IRYZ(I+1,MSB)*TROTI(I+1) 
      CTMP6 = IRYZM(I,MSB)*TROTM(I) 
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP*(
     &    GG(CTMP6,ZV2M*ZEM*CTMP4,ZV2P*ZEP*CTMP5)*Y(KYV2+LYCOL,I))
      CTMP4 = IRXZ(I,MSB)*TROTI(I) 
      CTMP5 = IRXZ(I+1,MSB)*TROTI(I+1) 
      CTMP6 = IRXZM(I,MSB)*TROTM(I) 
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP*(
     &            -GF(ZEM*CTMP4,CTMP6)*X(KXV1+LXCOL,I)
     &            -GF(ZEP*CTMP5,CTMP6)*X(KXV1+LXCOL,I+1))
      ENDIF
      IF (INERT2) THEN
      CTMP1 = -RHO(I)*TDROTI(I)*T(I)
      CTMP2 = -RHO(I+1)*TDROTI(I+1)*T(I+1)
      CTMP3 = -RHOM(I)*TDROTM(I)*TM(I)
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP*(
     &     GF(ZEM*CTMP1*JACOBI(I,MSB),CTMP3*JACOBM(I,MSB))*
     &     X(KXV1+LXCOL,I)
     &    +GF(ZEP*CTMP2*JACOBI(I+1,MSB),CTMP3*JACOBM(I,MSB))*
     &     X(KXV1+LXCOL,I+1))
      ENDIF
      ENDIF

 55   CONTINUE
    
      IF (V2U_M0.AND.ABS(RM(MSPL,2)).LT.1.E-3) RY(KYV2+LYROW,:) = 0.0

 60   CONTINUE

      IF (MSB.LT.2) GOTO 80
      IF (MSMI.LT.1) GOTO 80
 
      LXROW = (MSMI-1)*NXCOMP
      LXCOL = (MSA -1)*NXCOMP
      LYROW = (MSMI-1)*NYCOMP
      LYCOL = (MSA -1)*NYCOMP

      DO 70 I=2,NRP1
 
      TTMP2= -CONJG((CALPHA1+(CALPHA2-CALPHA3)*CI*RNTOR*TROTI(I))*
     &        CALPHA4/CALPHA2**2)
      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      ZV2P = (CS(I)/CSM(I  ))**IEXV2
      INCLUDE 'tent.inc'
 
      R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP2*(
     &    CONJG(FF(RGV1G11(i,msb),RGV1G11M(i-1,msb),RGV1G11M(i,msb)))
     &                                *X(KXV1+LXCOL,I)
     &   +CONJG(FFM(RGV1G11M(i-1,msb)))*X(KXV1+LXCOL,I-1)
     &   +CONJG(FFP(RGV1G11M(i,msb)))*X(KXV1+LXCOL,I+1)
     &   +CONJG(FGM(RGV1G12(i,msb)*zv2m,RGV1G12M(i-1,msb)))
     &                                 *Y(KYV2+LYCOL,I-1)
     &   +CONJG(FGP(RGV1G12(i,msb)*zv2p,RGV1G12M(i,msb)))
     &                                 *Y(KYV2+LYCOL,I))

      IF ((CALPHA2-CALPHA3).NE.0..AND.CALPHA7.GT.0.) THEN
      Z3M  = (CSM(I-1)/CS(I))**IEXB3
      Z3P  = (CSM(I)  /CS(I))**IEXB3
 
      IF (INERT2) THEN
      CTMP1 = IDRXX(I,MSB)*TDROTI(I)
      CTMP2 = IDRXXM(I-1,MSB)*TDROTM(I-1)
      CTMP3 = IDRXXM(I,MSB)*TDROTM(I)
      R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP1*(
     &                 CONJG(FF(CTMP1,CTMP2,CTMP3))*X(KXV1+LXCOL,I)
     &                +CONJG(FFM(CTMP2))*X(KXV1+LXCOL,I-1)
     &                +CONJG(FFP(CTMP3))*X(KXV1+LXCOL,I+1))
      ENDIF
 
      IF (INERT1) THEN
      CTMP1 = IRXY(I,MSB)*TROTI(I)
      CTMP2 = IRXYM(I-1,MSB)*TROTM(I-1)
      CTMP3 = IRXYM(I,MSB)*TROTM(I)
      R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP1*(
     &                   -CONJG(FGM(CTMP1*ZV2M,CTMP2))*Y(KYV2+LYCOL,I-1)
     &                   -CONJG(FGP(CTMP1*ZV2P,CTMP3))*Y(KYV2+LYCOL,I))
 
      CTMP1 = IRXZ(I,MSB)*TROTI(I)
      CTMP2 = IRXZM(I-1,MSB)*TROTM(I-1)
      CTMP3 = IRXZM(I,MSB)*TROTM(I)
      IF (KYV3.GT.0)
     &R (KXV1+LXROW,I) = R (KXV1+LXROW,I) + TTMP1*(
     &                   CONJG(FGM(CTMP1,CTMP2))*Y(KYV3+LYCOL,I-1)
     &                 + CONJG(FGP(CTMP1,CTMP3))*Y(KYV3+LYCOL,I))
      ENDIF
      ENDIF

 70   CONTINUE
 
      DO 75 I=1,NR
 
      TTMP2= -CONJG((CALPHA1+(CALPHA2-CALPHA3)*CI*RNTOR*TROTM(I))*
     &        CALPHA4/CALPHA2**2)
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      INCLUDE 'tophat.inc'
 
      RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) + TTMP2*(
     &   CONJG(GG(RGV2G22M(i,msb),RGV2G22(i,msb)*zem*zv2m,
     $        RGV2G22(i+1,msb)*zep*zv2p))
     &                             *Y(KYV2+LYCOL,I)
     &  + CONJG(GF(RGV1G12(i,msb)*zem,RGV1G12M(i,msb)))
     &                             *X(KXV1+LXCOL,I)
     &  + CONJG(GF(RGV1G12(i+1,msb)*zep,RGV1G12M(i,msb)))
     &                             *X(KXV1+LXCOL,I+1))
      IF (KYV3.GT.0)
     &RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP2*(
     &    CONJG(GG(RGV3G33M(i,msb),RGV3G33(i,msb)*zem,
     &        RGV3G33(i+1,msb)*zep))
     &                             *Y(KYV3+LYCOL,I))
      IF (KYX2.GT.0)
     &RY(KYX2+LYROW,I) = RY(KYX2+LYROW,I) + TTMP2*(
     &   CONJG(GG(JACOBM(i,msb),JACOBI(i,msb)*zem*zv2m,
     &         JACOBI(i+1,msb)*zep*zv2p))
     &                             *Y(KYX2+LYCOL,I))
      IF (KYX3.GT.0)
     &RY(KYX3+LYROW,I) = RY(KYX3+LYROW,I) + TTMP2*(
     &   CONJG(GG(JACOBM(i,msb),JACOBI(i,msb)*zem*zv2m,
     $         JACOBI(i+1,msb)*zep*zv2p))
     &                             *Y(KYX3+LYCOL,I))     

      IF ((CALPHA2-CALPHA3).NE.0..AND.CALPHA7.GT.0..AND.
     &    KXX1.GT.0.AND.KYX3.GT.0) THEN
         RY(KYX3+LYROW,I) = RY(KYX3+LYROW,I) + TTMP*(
     &    CONJG(GF(ZEM*TDROTI(I)*TB2(I,MSB),TDROTM(I)*TB2M(I,MSB)))*
     &    X(KXX1+LXCOL,I) 
     &   +CONJG(GF(ZEP*TDROTI(I+1)*TB2(I+1,MSB),TDROTM(I)*TB2M(I,MSB)))*
     &    X(KXX1+LXCOL,I+1)) 
      ENDIF

      IF (V2XKEY.EQ.0.OR.V2XKEY.EQ.1) THEN
         RY(KYPR + LYROW, I) = RY(KYPR + LYROW, I) + TTMP2*TTCINERT0*(
     &       CONJG(GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) )
     &                                *Y(KYPR+LYCOL,I))
         IF (IDIAMTE.EQ.2.AND.KYPE.GT.0) THEN
         RY(KYPE + LYROW, I) = RY(KYPE + LYROW, I) + TTMP2*(
     &       CONJG(GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) )
     &                                *Y(KYPE+LYCOL,I))
         ENDIF
      ENDIF

      IF (KYRHOP.GT.0)  
     &RY(KYRHOP + LYROW, I) = RY(KYRHOP + LYROW, I) + TTMP2*(
     &    CONJG(GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) )
     &                             *Y(KYRHOP+LYCOL,I))

 
      IF (IPERTURB.EQ.0.AND.INCKIN.EQ.1) THEN
         IF (V2XKEY.EQ.0 .OR. V2XKEY.EQ.2) THEN
            RY(KYPPARA+LYROW,I) = RY(KYPPARA+LYROW,I) + TTMP2*(
     &      CONJG(GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) )
     &     *Y(KYPPARA+LYCOL,I))
 
            RY(KYPPERP+LYROW,I) = RY(KYPPERP+LYROW,I) + TTMP2*(
     &      CONJG(GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) )
     &     *Y(KYPPERP+LYCOL,I))

            IF (KYPE.GT.0)
     &      RY(KYPE+LYROW,I) = RY(KYPE+LYROW,I) + TTMP2*(
     &      CONJG(GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) )
     &     *Y(KYPE+LYCOL,I))

            IF (KYPP.GT.0)
     &      RY(KYPP+LYROW,I) = RY(KYPP+LYROW,I) + TTMP2*(
     &      CONJG(GG(JacobM(i,msb), Jacobi(i,msb), Jacobi(i+1,msb)) )
     &     *Y(KYPP+LYCOL,I))
         ENDIF
      ENDIF
 
      IF ((CALPHA2-CALPHA3).NE.0..AND.CALPHA7.GT.0.) THEN
      IF (INERT2) THEN
      CTMP4 = IDRYX(I,MSB)*TDROTI(I)
      CTMP5 = IDRYX(I+1,MSB)*TDROTI(I+1)
      CTMP6 = IDRYXM(I,MSB)*TDROTM(I)
      RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) + TTMP1*(
     &                   CONJG(GF(ZEM*CTMP4,CTMP6))*X(KXV1+LXCOL,I)
     &                  +CONJG(GF(ZEP*CTMP5,CTMP6))*X(KXV1+LXCOL,I+1))
      ENDIF
      IF (INERT1) THEN
      CTMP1 = IRXY(I,MSB)*TROTI(I)
      CTMP2 = IRXY(I+1,MSB)*TROTI(I+1)
      CTMP3 = IRXYM(I,MSB)*TROTM(I)
      CTMP4 = IRYZ(I,MSB)*TROTI(I)
      CTMP5 = IRYZ(I+1,MSB)*TROTI(I+1)
      CTMP6 = IRYZM(I,MSB)*TROTM(I)
      RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) + TTMP1*(
     &                   CONJG(GF(ZEM*CTMP1,CTMP3))*X(KXV1+LXCOL,I)
     &                  +CONJG(GF(ZEP*CTMP2,CTMP3))*X(KXV1+LXCOL,I+1))
      IF (KYV3.GT.0)
     &RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) - TTMP1*
     &         CONJG(GG(CTMP6,ZEM*CTMP4,ZEP*CTMP5))*Y(KYV3+LYCOL,I)
      ENDIF

      IF (ABS(PVISC).GT.0..AND.IVISC.GT.0.AND.KYV3.GT.0) THEN
      ZKPAVT = ABS(RM(MSA,NSA)/QM(I) + RN(NSA))
     &                 * SQRT(REAL(PEQM(I))/RHOM(I))
      CTMP4 = IRYZ(I,MSB)*TROTI(I)
      CTMP5 = IRYZ(I+1,MSB)*TROTI(I+1)
      CTMP6 = IRYZM(I,MSB)*TROTM(I)
      IF (KYX2.GT.0) 
     &RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP1*PVISC*ZKPAVT*0.5*(
     &   CONJG(GG(CTMP6,ZV2M*ZEM*CTMP4,ZV2P*ZEP*CTMP5))*Y(KYX2+LYCOL,I))
 
      IF (KXX1.GT.0) THEN
      CTMP1 = -RHO(I)*TDROTI(I)*T(I)
      CTMP2 = -RHO(I+1)*TDROTI(I+1)*T(I+1)
      CTMP3 = -RHOM(I)*TDROTM(I)*TM(I)
      CTMP4 = IRXZ(I,MSB)*TROTI(I) 
      CTMP5 = IRXZ(I+1,MSB)*TROTI(I+1) 
      CTMP6 = IRXZM(I,MSB)*TROTM(I) 
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP1*PVISC*ZKPAVT*(
     &     CONJG(-GF(ZEM*CTMP4,CTMP6)*0.5+GF(ZEM*CTMP1*JACOBI(I,MSB),
     &      CTMP3*JACOBM(I,MSB)))*X(KXX1+LXCOL,I)
     &   + CONJG(-GF(ZEP*CTMP5,CTMP6)*0.5+GF(ZEP*CTMP2*JACOBI(I+1,MSB),
     &      CTMP3*JACOBM(I,MSB)))*X(KXX1+LXCOL,I+1))
      ENDIF
      ENDIF
          
      IF (INERT1) THEN
      CTMP4 = IRYZ(I,MSB)*TROTI(I) 
      CTMP5 = IRYZ(I+1,MSB)*TROTI(I+1) 
      CTMP6 = IRYZM(I,MSB)*TROTM(I) 
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP1*(
     &   CONJG(GG(CTMP6,ZV2M*ZEM*CTMP4,ZV2P*ZEP*CTMP5))*Y(KYV2+LYCOL,I))
      CTMP4 = IRXZ(I,MSB)*TROTI(I) 
      CTMP5 = IRXZ(I+1,MSB)*TROTI(I+1) 
      CTMP6 = IRXZM(I,MSB)*TROTM(I) 
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP1*(
     &            -CONJG(GF(ZEM*CTMP4,CTMP6))*X(KXV1+LXCOL,I)
     &            -CONJG(GF(ZEP*CTMP5,CTMP6))*X(KXV1+LXCOL,I+1))
      ENDIF
      IF (INERT2) THEN
      CTMP1 = -RHO(I)*TDROTI(I)*T(I)
      CTMP2 = -RHO(I+1)*TDROTI(I+1)*T(I+1)
      CTMP3 = -RHOM(I)*TDROTM(I)*TM(I)
      RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + TTMP1*(
     &     CONJG(GF(ZEM*CTMP1*JACOBI(I,MSB),CTMP3*JACOBM(I,MSB)))*
     &     X(KXV1+LXCOL,I)
     &    +CONJG(GF(ZEP*CTMP2*JACOBI(I+1,MSB),CTMP3*JACOBM(I,MSB)))*
     &     X(KXV1+LXCOL,I+1))
      ENDIF
      ENDIF

 75   CONTINUE

      IF (V2U_M0.AND.ABS(RM(MSMI,2)).LT.1.E-3) RY(KYV2+LYROW,:) = 0.0

 80   CONTINUE
 90   CONTINUE

C     TERMS WITHOUT CONVOLUTION

C     COMPUTE ROSENBLUTH AVALANCHE TERMS 
C     COMPUTE QUADRATIC TERMS FOR N=0 VERTICAL MOVEMENT
C     COMPUTE ALL TERMS IN REAL SPACE, THEN PERFORM FOURIER DECOMPOSITION
      IF ((KXJRE.GT.0.OR.KVSQLIN).AND.KJRER.LE.4) THEN
         IF (.NOT.ALLOCATED(RB1U)) THEN
            ALLOCATE(RB1U(NRP1,NCHI),RB2U(NRP1,NCHI),RB3U(NRP1,NCHI))
            ALLOCATE(RJ1U(NRP1,NCHI),RJ2U(NRP1,NCHI),RJ3U(NRP1,NCHI))
            ALLOCATE(RV1U(NRP1,NCHI),RV2U(NRP1,NCHI),RV3U(NRP1,NCHI))
            ALLOCATE(RRHOE(NRP1,NCHI),RPPE(NRP1,NCHI),
     &               RW1(NRP1,NCHI),RW2(NRP1,NCHI),RW3(NRP1,NCHI),
     &               RW4(NRP1,NCHI),RW5(NRP1,NCHI),RW6(NRP1,NCHI))
         ENDIF

         RB1U = (0.,0.)
         RB2U = (0.,0.)
         RB3U = (0.,0.)
         RJ1U = (0.,0.)
         RJ2U = (0.,0.)
         RJ3U = (0.,0.)
         RV1U = (0.,0.)
         RV2U = (0.,0.)
         RV3U = (0.,0.)
         RRHOE= (0.,0.)
         RPPE = (0.,0.)
         HCHI = 4.*PI/DFLOAT(NCHI)

         DO J=1,NCHI
         RCHI = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         DO I=2,NRP1
            DO MS=1,MSMAX
               LXROW = (MS-1)*NXCOMP
               LYROW = (MS-1)*NYCOMP
               CTMP1 = EXP(CI*RM(MS,2)*RCHI)

               RB1U(I,J)  = RB1U(I,J) + X(KXB1+LXROW,I)*CTMP1
               RB2U(I,J)  = RB2U(I,J) + Y(KYB2+LYROW,I)*CTMP1
               RB3U(I,J)  = RB3U(I,J) + Y(KYB3+LYROW,I)*CTMP1
               RJ1U(I,J)  = RJ1U(I,J) + Y(KYJ1+LYROW,I)*CTMP1
               RJ2U(I,J)  = RJ2U(I,J) + X(KXJ2U+LXROW,I)*CTMP1
               RJ3U(I,J)  = RJ3U(I,J) + X(KXJ3+LXROW,I)*CTMP1
               RV1U(I,J)  = RV1U(I,J) + X(KXV1+LXROW,I)*CTMP1
               RV2U(I,J)  = RV2U(I,J) + Y(KYV2+LYROW,I)*CTMP1
               RV3U(I,J)  = RV3U(I,J) + Y(KYV3+LYROW,I)*CTMP1
               RRHOE(I,J) = RRHOE(I,J)+ Y(KYRHOP+LYROW,I)*CTMP1
               RPPE(I,J)  = RPPE(I,J) + Y(KYPR+LYROW,I)*CTMP1
            ENDDO
         ENDDO
         ENDDO
      ENDIF         

      IF (KXJRE.GT.0.AND.KJRER.GE.1.AND.KJRER.LE.4) THEN
         IF (.NOT.ALLOCATED(FJRERHS)) 
     &      ALLOCATE(FJRERHS(NRP1,MSMAX),FJRECURL(NRP1,MSMAX,3))

         FJRERHS  = (0.,0.)
         FJRECURL = (0.,0.)

C        CONVERT JRE TO REAL SPACE AND STORE IN RW2
         RW2 = (0.,0.)
         DO J=1,NCHI
         RCHI = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         DO MS=1,MSMAX
         LXROW = (MS-1)*NXCOMP
         CTMP1 = EXP(CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            CTMP2 = X(KXJRE+LXROW,I)*CTMP1 
            RW2(I,J) = RW2(I,J) + CTMP2
         ENDDO
         ENDDO
         ENDDO

C        COMPUTE E-FIELD AND ROSENBLUTH AVALANCHE RHS=RW1 IN REAL SPACE
C        ALSO COMPUTE TOTAL JRE, JPARA AND VOLUME AVERAGED EPARA 
         RW1 = (0.,0.)
         RW4 = (0.,0.)
         TOT_RJA = (0.,0.)
         TOT_JRE = (0.,0.)
         TOT_JPA = (0.,0.)
         TOT_EPA = (0.,0.)
         COR_RJA = (0.,0.)
         COR_JRE = (0.,0.)
         COR_JPA = (0.,0.)
         COR_EPA = (0.,0.)
         DO J=1,NCHI
         DO I=2,NRP1
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))

            HB1 = RB1U(I,J)
            HB2 = (RB2U(I-1,J)*HH1+RB2U(I,J)*HH2)+DPSIDS(I)
            HB3 = (RB3U(I-1,J)*HH1+RB3U(I,J)*HH2)+
     &            T(I)*RJA(I,J)/REQ(I,J)**2
            HJ1 = RJ1U(I-1,J)*HH1+RJ1U(I,J)*HH2
            HJ2 = RJ2U(I,J)-DPSIDS(I)*TP(I) 
            HJ3 = RJ3U(I,J)-RJA(I,J)*(PPEQ(I)+T(I)*TP(I)/REQ(I,J)**2)
            
            REB     = SQRT(G11L(I,J)*ABS(HB1)**2 +
     &                     2.*G12L(I,J)*REAL(HB1*CONJG(HB2)) +
     &                     G22L(I,J)*ABS(HB2)**2 +
     &                     REQ(I,J)**2*ABS(HB3)**2)/RJA(I,J)

            CTMP1 = (G11L(I,J)*HB1*HJ1+G12L(I,J)*(HB1*HJ2+HB2*HJ1)+
     &              G22L(I,J)*HB2*HJ2+REQ(I,J)**2*HB3*HJ3)/RJA(I,J)**2
            REEPARA = RESIST(I)*REAL(CTMP1/REB-RW2(I,J))
            REEPARA =-RE_E0*REEPARA

            IF (I.LT.NRP1) RTMP = CS(I+1) - CS(I-1)
            IF (I.EQ.NRP1) RTMP = 1. - CS(I-1)
            TOT_RJA = TOT_RJA + RJA(I,J)*RTMP
            TOT_JRE = TOT_JRE + RW2(I,J)*RJA(I,J)*RTMP
            TOT_JPA = TOT_JPA + CTMP1/REB*RJA(I,J)*RTMP
            TOT_EPA = TOT_EPA + REEPARA*RJA(I,J)*RTMP

            IF (I.LE.NRES) THEN
            COR_RJA = COR_RJA + RJA(I,J)*RTMP
            COR_JRE = COR_JRE + RW2(I,J)*RJA(I,J)*RTMP
            COR_JPA = COR_JPA + CTMP1/REB*RJA(I,J)*RTMP
            COR_EPA = COR_EPA + REEPARA*RJA(I,J)*RTMP
            ENDIF

            RE_GAMMA = GNEOFUNC(I)
            RE_Z     = RE_CONST(3)
            RTMP     = 1.-1./REEPARA+4.*PI*(RE_Z+1)**2/3./RE_GAMMA/
     &                 (RE_Z+5)/(REEPARA**2+4./RE_GAMMA**2-1.)
            IF (RTMP.LT.1.E-4) RTMP = 1.E-4
            RTMP     = SQRT(PI*RE_GAMMA/3./(RE_Z+5)/RTMP)
            RTMP     = RE_T0/RE_CONST(2)*(REEPARA-1.)*RTMP
            RW1(I,J) = RE_RP_FAC*RTMP*RW2(I,J)

            RW4(I,J) = RESIST(I)*RW2(I,J)/REB
         ENDDO
         ENDDO
         RTMP = PI/DFLOAT(NCHI)
         TOT_RJA = TOT_RJA*RTMP
         TOT_JRE = TOT_JRE*RTMP
         TOT_JPA = TOT_JPA*RTMP
         TOT_EPA = TOT_EPA*RTMP
         COR_RJA = COR_RJA*RTMP
         COR_JRE = COR_JRE*RTMP
         COR_JPA = COR_JPA*RTMP
         COR_EPA = COR_EPA*RTMP

C        FOURIER TRANSFORM OF ROSENBLUTH AVALANCHE RHS TERM FJRERHS
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FJRERHS(I,MS) = FJRERHS(I,MS) + RW1(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FJRERHS = FJRERHS/DFLOAT(NCHI)

C        COMPUTE CURL(JRE...) TERM AND FOURIER TRANSFORM
         IF (KJRER.EQ.1.OR.KJRER.EQ.2) THEN
         RW1 = (0.,0.)
         RW2 = (0.,0.)
         RW3 = (0.,0.)
         DO J=1,NCHI
         IF (J.EQ.1)    J1 = NCHI
         IF (J.GT.1)    J1 = J-1
         IF (J.EQ.NCHI) J2 = 1
         IF (J.LT.NCHI) J2 = J+1
         DO I=2,NR
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))
            HH3  = CSH(I+1)/(CSH(I)+CSH(I+1))
            HH4  = CSH(I)/(CSH(I)+CSH(I+1))

            CTMP  = (RB3U(I-1,J1)*HH1+RB3U(I,J1)*HH2)/
     &               RJA(I,J1)*REQ(I,J1)**2
            CTMP1 = RW4(I,J1)*CTMP
            CTMP  = (RB3U(I-1,J2)*HH1+RB3U(I,J2)*HH2)/
     &               RJA(I,J2)*REQ(I,J2)**2
            CTMP2 = RW4(I,J2)*CTMP
            RW1(I,J) = (CTMP2-CTMP1)/HCHI

            CTMP  = (RB3U(I-1,J)*HH1+RB3U(I,J)*HH2)/
     &               RJA(I,J)*REQ(I,J)**2
            CTMP1 = RW4(I,J)*CTMP
            CTMP  = (RB3U(I,J)*HH3+RB3U(I+1,J)*HH4)/
     &               RJA(I+1,J)*REQ(I+1,J)**2
            CTMP2 = RW4(I+1,J)*CTMP
            RW2(I,J) =-(CTMP2-CTMP1)/CSH(I)

            CTMP  = ((RB2U(I-1,J)*HH1+RB2U(I,J)*HH2))*
     &              G22L(I,J) + G12L(I,J)*RB1U(I,J)
            CTMP1 = RW4(I,J)*CTMP/RJA(I,J)
            CTMP  = ((RB2U(I,J)*HH3+RB2U(I+1,J)*HH4))*
     &              G22L(I+1,J) + G12L(I+1,J)*RB1U(I+1,J)
            CTMP2 = RW4(I+1,J)*CTMP/RJA(I+1,J)
            RW3(I,J) = (CTMP2-CTMP1)/CSH(I)

            CTMP1 = RB2U(I,J1)*G12LM(I,J1) +
     &              (RB1U(I,J1)+RB1U(I+1,J1))/2.*G11LM(I,J1)
            CTMP1 = (RW4(I,J1)+RW4(I+1,J1))/2.*CTMP1/RJAM(I,J1)
            CTMP2 = RB2U(I,J2)*G12LM(I,J2) +
     &              (RB1U(I,J2)+RB1U(I+1,J2))/2.*G11LM(I,J2)
            CTMP2 = (RW4(I,J2)+RW4(I+1,J2))/2.*CTMP2/RJAM(I,J2)
            RW3(I,J) = RW3(I,J) - (CTMP2-CTMP1)/HCHI
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         RW2(NRP1,J) = RW2(NR,J)
         RW3(NRP1,J) = RW3(NR,J)
         ENDDO

         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FJRECURL(I,MS,1) = FJRECURL(I,MS,1) + RW1(I,J)*CTMP1
            FJRECURL(I,MS,2) = FJRECURL(I,MS,2) + RW2(I,J)*CTMP1
            FJRECURL(I,MS,3) = FJRECURL(I,MS,3) + RW3(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FJRECURL = FJRECURL/DFLOAT(NCHI)
         ENDIF
      ENDIF

      IF (KVSQLIN) THEN  
C        FIRST COMPUTE AND STORE RW4=BEQ**2 AT INTEGER GRID
C        AND RW5=HAT_V^2, RW6=HAT_V^3 AT HALF_INTEGER GRID
         RW4 = (0.,0.)
         RW5 = (0.,0.)
         RW6 = (0.,0.)
         DO J=1,NCHI
         DO I=2,NRP1
            RW4(I,J) = G22L(I,J)*DPSIDS(I)**2/RJA(I,J)**2 +
     &                 T(I)**2/REQ(I,J)**2
         ENDDO
         ENDDO

         DO J=1,NCHI
         DO I=2,NR
            CTMP1=(RW4(I,J)+RW4(I+1,J))/2.
            CTMP2=(RV1U(I,J)+RV1U(I+1,J))/2.
            CTMP3=RV3U(I,J)-DPSIDSM(I)*G12LM(I,J)*CTMP2/RJAM(I,J)/CTMP1
            CTMP4=RV2U(I,J)/CTMP1
            RW5(I,J) = CTMP3*DPSIDSM(I)/RJAM(I,J)+CTMP4*TM(I)
            RW6(I,J) = CTMP3*TM(I)/REQM(I,J)**2-
     &                 CTMP4*DPSIDSM(I)*G22LM(I,J)/RJAM(I,J)
         ENDDO
         RW5(NRP1,J) = RW5(NR,J)
         RW6(NRP1,J) = RW6(NR,J)
         ENDDO
      ENDIF

C     COMPUTE QUADRATIC TERM JXB ASSOCIATED WITH VERTICAL INSTABILITY FOR N=0
      IF (KVSQLIN.AND.KVSQL(1)) THEN  
         IF (.NOT.ALLOCATED(FVSQJXB))
     &      ALLOCATE(FVSQJXB(NRP1,MSMAX,3))

         RW1 = (0.,0.)
         RW2 = (0.,0.)
         RW3 = (0.,0.)
         DO J=1,NCHI
         DO I=2,NR
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))

            HB1I = RB1U(I,J)
            HB2I = RB2U(I-1,J)*HH1+RB2U(I,J)*HH2
            HB3I = RB3U(I-1,J)*HH1+RB3U(I,J)*HH2
            HJ1I = RJ1U(I-1,J)*HH1+RJ1U(I,J)*HH2
            HJ2I = RJ2U(I,J)
            HJ3I = RJ3U(I,J)

            HB1M = (RB1U(I,J)+RB1U(I+1,J))/2.
            HB2M = RB2U(I,J)
            HB3M = RB3U(I,J)
            HJ1M = RJ1U(I,J)
            HJ2M = (RJ2U(I,J)+RJ2U(I+1,J))/2.
            HJ3M = (RJ3U(I,J)+RJ3U(I+1,J))/2.

            CTMP2 =-(DPSIDS(I)/RJA(I,J))**2*G12L(I,J)/RW4(I,J)
            CTMP3 =-T(I)*DPSIDS(I)*G12L(I,J)/RJA(I,J)/RW4(I,J)/
     &              REQ(I,J)**2
            RW1(I,J) = (HJ2I*HB3I-HJ3I*HB2I) +
     &                 (HJ3I*HB1I-HJ1I*HB3I)*CTMP2 +
     &                 (HJ1I*HB2I-HJ2I*HB1I)*CTMP3 
                
            CTMP = (RW4(I,J)+RW4(I+1,J))/2.
            CTMP2 = TM(I)/CTMP
            CTMP3 =-DPSIDSM(I)*G22LM(I,J)/RJAM(I,J)/CTMP
            RW2(I,J) = (HJ3M*HB1M-HJ1M*HB3M)*CTMP2 +
     &                 (HJ1M*HB2M-HJ2M*HB1M)*CTMP3

            CTMP2 = DPSIDSM(I)/RJAM(I,J)
            CTMP3 = TM(I)/REQM(I,J)**2
            RW3(I,J) = (HJ3M*HB1M-HJ1M*HB3M)*CTMP2 +
     &                 (HJ1M*HB2M-HJ2M*HB1M)*CTMP3

         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         RW2(NRP1,J) = RW2(NR,J)
         RW3(NRP1,J) = RW3(NR,J)
         ENDDO

         FVSQJXB = (0.,0.)
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQJXB(I,MS,1) = FVSQJXB(I,MS,1) + RW1(I,J)*CTMP1
            FVSQJXB(I,MS,2) = FVSQJXB(I,MS,2) + RW2(I,J)*CTMP1
            FVSQJXB(I,MS,3) = FVSQJXB(I,MS,3) + RW3(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQJXB = FVSQJXB/DFLOAT(NCHI)
      ENDIF

C     COMPUTE QUADRATIC TERM RHO(V,GRAD)V ASSOCIATED WITH N=0
      IF (KVSQLIN.AND.KVSQL(2)) THEN  
         IF (.NOT.ALLOCATED(FVSQVGV))
     &      ALLOCATE(FVSQVGV(NRP1,MSMAX,3),RW7(NRP1,NCHI))

C        COMPUTE RW7=V**2 AT INTEGER GRID
         RW7 = (0.,0.)
         DO J=1,NCHI
         DO I=2,NR
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))

            CTMP2 = RW5(I-1,J)*HH1 + RW5(I,J)*HH2
            CTMP3 = RW6(I-1,J)*HH1 + RW6(I,J)*HH2

            RW7(I,J) = G11L(I,J)*RV1U(I,J)**2 +
     &                 2.*G12L(I,J)*RV1U(I,J)*CTMP2 +
     &                 G22L(I,J)*CTMP2**2 + (REQ(I,J)*CTMP3)**2
         ENDDO
         RW7(NRP1,J) = RW7(NR,J)
         ENDDO

         RW1 = (0.,0.)
         RW2 = (0.,0.)
         RW3 = (0.,0.)
         DO J=1,NCHI
         IF (J.EQ.1)    J1 = NCHI 
         IF (J.GT.1)    J1 = J-1
         IF (J.EQ.NCHI) J2 = 1
         IF (J.LT.NCHI) J2 = J+1
         DO I=2,NR
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))
            
C           FIRST TERM: GRAD(V**2)/2
            CTMP1 = (RW7(I+1,J)-RW7(I-1,J))/(CS(I+1)-CS(I-1))
            IF (J.EQ.1)    CTMP = (RW7(I,2)-RW7(I,NCHI))/HCHI
            IF (J.EQ.NCHI) CTMP = (RW7(I,1)-RW7(I,NCHI-1))/HCHI
            IF (J.GT.1.AND.J.LT.NCHI) 
     &         CTMP = (RW7(I,J+1)-RW7(I,J-1))/HCHI
            CTMP2 =-CTMP*(DPSIDS(I)/RJA(I,J))**2*G12L(I,J)/RW4(I,J)
            RW1(I,J) = (CTMP1+CTMP2)/2.

            CTMP1 = (RW7(I,J1)+RW7(I+1,J1))/2.
            CTMP2 = (RW7(I,J2)+RW7(I+1,J2))/2.
            CTMP = (CTMP2-CTMP1)/HCHI
            CTMP3 = (RW4(I,J)+RW4(I+1,J))/2.
            CTMP4 = CTMP*TM(I)/CTMP3
            RW2(I,J) = CTMP4/2.

            CTMP4 = CTMP*DPSIDSM(I)/RJAM(I,J)
            RW3(I,J) = CTMP4/2.
            
C           SECOND TERM: -VX(CURL V)
            CTMP1 = REQ(I,J1)**2*(RW6(I-1,J1)*HH1+RW6(I,J1)*HH2)
            CTMP2 = REQ(I,J2)**2*(RW6(I-1,J2)*HH1+RW6(I,J2)*HH2)
            CTMP  = (CTMP2-CTMP1)/HCHI
            CTMP2 = RV2U(I-1,J)*HH1+RV2U(I,J)*HH2
            CTMP4 = DPSIDS(I)*G12L(I,J)/RW4(I,J)*CTMP2
            RW1(I,J) = RW1(I,J) - CTMP*CTMP4

            CTMP  =-(REQM(I,J)**2*RW6(I,J)-REQM(I-1,J)**2*RW6(I-1,J))/
     &              (CS(I+1)-CS(I-1))*2.
            CTMP3 = RV3U(I-1,J)*HH1+RV3U(I,J)*HH2
            CTMP4 = DPSIDS(I)*G22L(I,J)/RW4(I,J)*CTMP2 -
     &              RJA(I,J)*T(I)/REQ(I,J)**2*CTMP3
            RW1(I,J) = RW1(I,J) - CTMP*CTMP4

            CTMP4 = RJA(I,J)*T(I)/RW4(I,J)*CTMP2+DPSIDS(I)*CTMP3
            CTMP1 = G12LM(I-1,J)*(RV1U(I-1,J)+RV1U(I,J))/2. +
     &              G22LM(I-1,J)*RW5(I-1,J)
            CTMP2 = G12LM(I,J)*(RV1U(I,J)+RV1U(I+1,J))/2. +
     &              G22LM(I,J)*RW5(I,J)
            CTMP  = (CTMP2-CTMP1)/(CS(I+1)-CS(I-1))*2.
            CTMP1 = G11L(I,J1)*RV1U(I,J1) + 
     &              G12L(I,J1)*(RW5(I-1,J1)*HH1+RW5(I,J1)*HH2)
            CTMP2 = G11L(I,J2)*RV1U(I,J2) + 
     &              G12L(I,J2)*(RW5(I-1,J2)*HH1+RW5(I,J2)*HH2)
            CTMP3 = (CTMP2-CTMP1)/HCHI
            CTMP  = CTMP - CTMP3
            RW1(I,J) = RW1(I,J) - CTMP*CTMP4

            CTMP1 = REQM(I,J1)**2*RW6(I,J1)
            CTMP2 = REQM(I,J2)**2*RW6(I,J2)
            CTMP  = (CTMP2-CTMP1)/HCHI
            CTMP3 = (RV1U(I,J)/RW4(I,J)+RV1U(I+1,J)/RW4(I+1,J))/2.
            CTMP4 = RJAM(I,J)*RV3U(I,J) - DPSIDSM(I)*G12LM(I,J)*CTMP3
            CTMP5 =-RJAM(I,J)*RV2U(I,J)
            RW2(I,J) = RW2(I,J) - CTMP*CTMP4
            RW3(I,J) = RW3(I,J) - CTMP*CTMP5

            CTMP1 = REQM(I-1,J)**2*RW6(I-1,J)
            CTMP2 = REQM(I+1,J)**2*RW6(I+1,J)
            CTMP  =-(CTMP2-CTMP1)/(CSM(I+1)-CSM(I-1))
            CTMP4 =-DPSIDSM(I)*G22LM(I,J)*CTMP3
            CTMP5 = RJAM(I,J)*TM(I)/REQM(I,J)**2*
     &              (RV1U(I,J)+RV1U(I+1,J))/2.
            RW2(I,J) = RW2(I,J) - CTMP*CTMP4
            RW3(I,J) = RW3(I,J) - CTMP*CTMP5

            CTMP1 = G12L(I,J)*RV1U(I,J)
            CTMP2 = G12L(I+1,J)*RV1U(I+1,J)
            CTMP  = (CTMP2-CTMP1)/CSH(I)
            CTMP1 = G22LM(I-1,J)*RW5(I-1,J)
            CTMP2 = G22LM(I+1,J)*RW5(I+1,J)
            CTMP  = CTMP + (CTMP2-CTMP1)/(CSM(I+1)-CSM(I-1))
            CTMP1 = G11LM(I,J1)*(RV1U(I,J1)+RV1U(I+1,J1))/2. +
     &              G12LM(I,J1)*RW5(I,J1)
            CTMP2 = G11LM(I,J2)*(RV1U(I,J2)+RV1U(I+1,J2))/2. +
     &              G12LM(I,J2)*RW5(I,J2)
            CTMP  = CTMP - (CTMP2-CTMP1)/HCHI
            CTMP4 =-RJAM(I,J)*TM(I)*CTMP3
            CTMP5 =-DPSIDSM(I)*(RV1U(I,J)+RV1U(I+1,J))/2.
            RW2(I,J) = RW2(I,J) - CTMP*CTMP4
            RW3(I,J) = RW3(I,J) - CTMP*CTMP5

C           ADD EXTRA FACTOR ASSOCIATED WITH TOTAL DENSITY
            CTMP = RRHOE(I-1,J)*HH1+RRHOE(I,J)*HH2
            RW1(I,J) = RW1(I,J)*(RHO(I)+CTMP)
            RW2(I,J) = RW2(I,J)*(RHOM(I)+RRHOE(I,J))
            RW3(I,J) = RW3(I,J)*(RHOM(I)+RRHOE(I,J))
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         RW2(NRP1,J) = RW2(NR,J)
         RW3(NRP1,J) = RW3(NR,J)
         ENDDO

         FVSQVGV = (0.,0.)
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQVGV(I,MS,1) = FVSQVGV(I,MS,1) + RW1(I,J)*CTMP1
            FVSQVGV(I,MS,2) = FVSQVGV(I,MS,2) + RW2(I,J)*CTMP1
            FVSQVGV(I,MS,3) = FVSQVGV(I,MS,3) + RW3(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQVGV = FVSQVGV/DFLOAT(NCHI)
      ENDIF

C     COMPUTE QUADRATIC TERM CURL(VXB) ASSOCIATED WITH N=0
C     OBSOLETE, SEE NEXT ALTERNATIVE WAY
      IF (KVSQLIN.AND.KVSQL(3).AND.1.EQ.0) THEN  
         IF (.NOT.ALLOCATED(FVSQVXB))
     &      ALLOCATE(FVSQVXB(NRP1,MSMAX,3))

         RW1 = (0.,0.)
         RW2 = (0.,0.)
         RW3 = (0.,0.)
         DO J=1,NCHI
         IF (J.EQ.1)    J1 = NCHI
         IF (J.GT.1)    J1 = J-1
         IF (J.EQ.NCHI) J2 = 1
         IF (J.LT.NCHI) J2 = J+1
         DO I=2,NR
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))
            
C           FIRST TERM: B DOT GRAD(V DOT GRAD F)
C           NB: 1ST AND 2ND TERMS FOR RW1 ARE COMBINED
            CTMP1 = (RB2U(I-1,J1)*HH1+RB2U(I,J1)*HH2)*RV1U(I,J1)
            CTMP2 = (RB2U(I-1,J2)*HH1+RB2U(I,J2)*HH2)*RV1U(I,J2)
            RW1(I,J) = (CTMP2-CTMP1)/HCHI

            CTMP1 = (RW5(I+1,J)-RW5(I-1,J))/(CSM(I+1)-CSM(I-1))
            CTMP2 = CTMP1*(RB1U(I,J)+RB1U(I+1,J))/2.

            IF (J.EQ.1)    CTMP = (RW5(I,2)-RW5(I,NCHI))/HCHI
            IF (J.EQ.NCHI) CTMP = (RW5(I,1)-RW5(I,NCHI-1))/HCHI
            IF (J.GT.1.AND.J.LT.NCHI) CTMP=(RW5(I,J+1)-RW5(I,J-1))/HCHI
            CTMP3 = CTMP*RB2U(I,J)
            RW2(I,J) = CTMP2 + CTMP3

            CTMP1 = (RW6(I+1,J)-RW6(I-1,J))/(CSM(I+1)-CSM(I-1))
            CTMP2 = CTMP1*(RB1U(I,J)+RB1U(I+1,J))/2.

            IF (J.EQ.1)    CTMP = (RW6(I,2)-RW6(I,NCHI))/HCHI
            IF (J.EQ.NCHI) CTMP = (RW6(I,1)-RW6(I,NCHI-1))/HCHI
            IF (J.GT.1.AND.J.LT.NCHI) CTMP=(RW6(I,J+1)-RW6(I,J-1))/HCHI
            CTMP3 = CTMP*RB2U(I,J)
            RW3(I,J) = CTMP2 + CTMP3

C           SECOND TERM: -DIV(*V)
            CTMP1 =-DPSIDS(I)**2*G12L(I,J1)/RJA(I,J1)**2/
     &              RW4(I,J1)*RV1U(I,J1) + T(I)/RW4(I,J1)*
     &              (RV2U(I-1,J1)*HH1+RV2U(I,J1)*HH2) +
     &              DPSIDS(I)/RJA(I,J1)*
     &              (RV3U(I-1,J1)*HH1+RV3U(I,J1)*HH2)
            CTMP1 = CTMP1*RB1U(I,J1)
            CTMP2 =-DPSIDS(I)**2*G12L(I,J2)/RJA(I,J2)**2/
     &              RW4(I,J2)*RV1U(I,J2) + T(I)/RW4(I,J2)*
     &              (RV2U(I-1,J2)*HH1+RV2U(I,J2)*HH2) +
     &              DPSIDS(I)/RJA(I,J2)*
     &              (RV3U(I-1,J2)*HH1+RV3U(I,J2)*HH2)
            CTMP2 = CTMP2*RB1U(I,J2)
            CTMP =  (CTMP2-CTMP1)/HCHI
            RW1(I,J) = RW1(I,J) - CTMP

            CTMP = RB2U(I,J)*(RV1U(I+1,J)-RV1U(I,J))/CSH(I)+
     &             (RV1U(I,J)+RV1U(I+1,J))/2.*
     &             (RB2U(I+1,J)-RB2U(I-1,J))/(CSM(I+1)-CSM(I-1))
            CTMP3 = (RW4(I,J1)+RW4(I+1,J1))/2.
            CTMP1 =-DPSIDSM(I)**2*G12LM(I,J1)/RJAM(I,J1)**2/
     &              CTMP3*(RV1U(I,J1)+RV1U(I+1,J1))/2. + 
     &              TM(I)/CTMP3*RV2U(I,J1) + 
     &              DPSIDSM(I)/RJAM(I,J1)*RV3U(I,J1)
            CTMP1 = CTMP1*RB2U(I,J1)
            CTMP3 = (RW4(I,J2)+RW4(I+1,J2))/2.
            CTMP2 =-DPSIDSM(I)**2*G12LM(I,J2)/RJAM(I,J2)**2/
     &              CTMP3*(RV1U(I,J2)+RV1U(I+1,J2))/2. + 
     &              TM(I)/CTMP3*RV2U(I,J2) + 
     &              DPSIDSM(I)/RJAM(I,J2)*RV3U(I,J2)
            CTMP2 = CTMP2*RB2U(I,J2)
            CTMP = CTMP + (CTMP2-CTMP1)/HCHI
            RW2(I,J) = RW2(I,J) - CTMP

            CTMP = RB3U(I,J)*(RV1U(I+1,J)-RV1U(I,J))/CSH(I)+
     &             (RV1U(I,J)+RV1U(I+1,J))/2.*
     &             (RB3U(I+1,J)-RB3U(I-1,J))/(CSM(I+1)-CSM(I-1))
            CTMP3 = (RW4(I,J1)+RW4(I+1,J1))/2.
            CTMP1 =-DPSIDSM(I)**2*G12LM(I,J1)/RJAM(I,J1)**2/
     &              CTMP3*(RV1U(I,J1)+RV1U(I+1,J1))/2. + 
     &              TM(I)/CTMP3*RV2U(I,J1) + 
     &              DPSIDSM(I)/RJAM(I,J1)*RV3U(I,J1)
            CTMP1 = CTMP1*RB3U(I,J1)
            CTMP3 = (RW4(I,J2)+RW4(I+1,J2))/2.
            CTMP2 =-DPSIDSM(I)**2*G12LM(I,J2)/RJAM(I,J2)**2/
     &              CTMP3*(RV1U(I,J2)+RV1U(I+1,J2))/2. + 
     &              TM(I)/CTMP3*RV2U(I,J2) + 
     &              DPSIDSM(I)/RJAM(I,J2)*RV3U(I,J2)
            CTMP2 = CTMP2*RB3U(I,J2)
            CTMP = CTMP + (CTMP2-CTMP1)/HCHI
            RW3(I,J) = RW3(I,J) - CTMP
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         RW2(NRP1,J) = RW2(NR,J)
         RW3(NRP1,J) = RW3(NR,J)
         ENDDO

         FVSQVXB = (0.,0.)
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQVXB(I,MS,1) = FVSQVXB(I,MS,1) + RW1(I,J)*CTMP1
            FVSQVXB(I,MS,2) = FVSQVXB(I,MS,2) + RW2(I,J)*CTMP1
            FVSQVXB(I,MS,3) = FVSQVXB(I,MS,3) + RW3(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQVXB = FVSQVXB/DFLOAT(NCHI)
      ENDIF

C     COMPUTE QUADRATIC TERM CURL(VXB) ASSOCIATED WITH N=0
C     ALTERNATIVE AND MORE EFFICIENT WAY
      IF (KVSQLIN.AND.KVSQL(3).AND.1.EQ.1) THEN  
         IF (.NOT.ALLOCATED(FVSQVXB))
     &      ALLOCATE(FVSQVXB(NRP1,MSMAX,3))

         RW1 = (0.,0.)
         RW2 = (0.,0.)
         RW3 = (0.,0.)
         DO J=1,NCHI
         IF (J.EQ.1)    J1 = NCHI
         IF (J.GT.1)    J1 = J-1
         IF (J.EQ.NCHI) J2 = 1
         IF (J.LT.NCHI) J2 = J+1
         DO I=2,NR
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))
            HH3  = CSH(I+1)/(CSH(I)+CSH(I+1))
            HH4  = CSH(I)/(CSH(I)+CSH(I+1))
            
            CTMP1 = (RB2U(I-1,J1)*HH1+RB2U(I,J1)*HH2)*RV1U(I,J1)
            CTMP2 = (RB2U(I-1,J2)*HH1+RB2U(I,J2)*HH2)*RV1U(I,J2)
            RW1(I,J) = (CTMP2-CTMP1)/HCHI

            CTMP1 = (RW5(I-1,J1)*HH1+RW5(I,J1)*HH2)*RB1U(I,J1)
            CTMP2 = (RW5(I-1,J2)*HH1+RW5(I,J2)*HH2)*RB1U(I,J2)
            RW1(I,J) = RW1(I,J) - (CTMP2-CTMP1)/HCHI

            CTMP1 = (RB2U(I-1,J)*HH1+RB2U(I,J)*HH2)*RV1U(I,J)
            CTMP2 = (RB2U(I,J)*HH3+RB2U(I+1,J)*HH4)*RV1U(I+1,J)
            RW2(I,J) =-(CTMP2-CTMP1)/CSH(I)

            CTMP1 = (RW5(I-1,J)*HH1+RW5(I,J)*HH2)*RB1U(I,J)
            CTMP2 = (RW5(I,J)*HH3+RW5(I+1,J)*HH4)*RB1U(I+1,J)
            RW2(I,J) = RW2(I,J) + (CTMP2-CTMP1)/CSH(I)

            CTMP1 = (RW6(I-1,J)*HH1+RW6(I,J)*HH2)*RB1U(I,J)
            CTMP2 = (RW6(I,J)*HH3+RW6(I+1,J)*HH4)*RB1U(I+1,J)
            RW3(I,J) = (CTMP2-CTMP1)/CSH(I)

            CTMP1 = (RB3U(I-1,J)*HH1+RB3U(I,J)*HH2)*RV1U(I,J)
            CTMP2 = (RB3U(I,J)*HH3+RB3U(I+1,J)*HH4)*RV1U(I+1,J)
            RW3(I,J) = RW3(I,J) - (CTMP2-CTMP1)/CSH(I)

            CTMP1 = RW5(I,J1)*RB3U(I,J1)-RW6(I,J1)*RB2U(I,J1)
            CTMP2 = RW5(I,J2)*RB3U(I,J2)-RW6(I,J2)*RB2U(I,J2)
            RW3(I,J) = RW3(I,J) - (CTMP2-CTMP1)/HCHI
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         RW2(NRP1,J) = RW2(NR,J)
         RW3(NRP1,J) = RW3(NR,J)
         ENDDO

         FVSQVXB = (0.,0.)
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQVXB(I,MS,1) = FVSQVXB(I,MS,1) + RW1(I,J)*CTMP1
            FVSQVXB(I,MS,2) = FVSQVXB(I,MS,2) + RW2(I,J)*CTMP1
            FVSQVXB(I,MS,3) = FVSQVXB(I,MS,3) + RW3(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQVXB = FVSQVXB/DFLOAT(NCHI)
      ENDIF

C     COMPUTE QUADRATIC TERM (V DOT GRAD P) ASSOCIATED WITH N=0
      IF (KVSQLIN.AND.KVSQL(4)) THEN  
         IF (.NOT.ALLOCATED(FVSQVGP))
     &      ALLOCATE(FVSQVGP(NRP1,MSMAX))

         RW1 = (0.,0.)
         DO J=1,NCHI
         DO I=2,NR
            CTMP1 = (RV1U(I,J)+RV1U(I+1,J))/2.
            CTMP2 = (RPPE(I+1,J)-RPPE(I-1,J))/(CSM(I+1)-CSM(I-1))
            RW1(I,J) = CTMP1*CTMP2*RJAM(I,J)

            IF (J.EQ.1)    CTMP = (RPPE(I,2)-RPPE(I,NCHI))/HCHI
            IF (J.EQ.NCHI) CTMP = (RPPE(I,1)-RPPE(I,NCHI-1))/HCHI
            IF (J.GT.1.AND.J.LT.NCHI) 
     &         CTMP = (RPPE(I,J+1)-RPPE(I,J-1))/HCHI
            CTMP3 = CTMP*RW5(I,J)*RJAM(I,J)
            RW1(I,J) = RW1(I,J) + CTMP3
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         ENDDO

         FVSQVGP = (0.,0.)
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQVGP(I,MS) = FVSQVGP(I,MS) + RW1(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQVGP = FVSQVGP/DFLOAT(NCHI)
      ENDIF

C     COMPUTE QUADRATIC TERM (GAMMA P DIV V) ASSOCIATED WITH N=0
      IF (KVSQLIN.AND.KVSQL(5)) THEN  
         IF (.NOT.ALLOCATED(FVSQPDV))
     &      ALLOCATE(FVSQPDV(NRP1,MSMAX))

         RW1 = (0.,0.)
         DO J=1,NCHI
         IF (J.EQ.1)    J1 = NCHI
         IF (J.GT.1)    J1 = J-1
         IF (J.EQ.NCHI) J2 = 1
         IF (J.LT.NCHI) J2 = J+1
         DO I=2,NR
            CTMP = (RJA(I+1,J)*RV1U(I+1,J)-RJA(I,J)*RV1U(I,J))/CSH(I)
            CTMP3 = (RW4(I,J1)+RW4(I+1,J1))/2.
            CTMP1 =-DPSIDSM(I)**2*G12LM(I,J1)/RJAM(I,J1)/
     &              CTMP3*(RV1U(I,J1)+RV1U(I+1,J1))/2. + 
     &              RJAM(I,J1)*TM(I)/CTMP3*RV2U(I,J1) + 
     &              DPSIDSM(I)*RV3U(I,J1)
            CTMP3 = (RW4(I,J2)+RW4(I+1,J2))/2.
            CTMP2 =-DPSIDSM(I)**2*G12LM(I,J2)/RJAM(I,J2)/
     &              CTMP3*(RV1U(I,J2)+RV1U(I+1,J2))/2. + 
     &              RJAM(I,J2)*TM(I)/CTMP3*RV2U(I,J2) + 
     &              DPSIDSM(I)*RV3U(I,J2)
            CTMP = CTMP + (CTMP2-CTMP1)/HCHI
            RW1(I,J) = CTMP*GAMARM(I)*RPPE(I,J)
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         ENDDO

         FVSQPDV = (0.,0.)
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQPDV(I,MS) = FVSQPDV(I,MS) + RW1(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQPDV = FVSQPDV/DFLOAT(NCHI)
      ENDIF

C     COMPUTE QUADRATIC TERM DIV(RHO_1 V) ASSOCIATED WITH N=0
      IF (KVSQLIN.AND.KVSQL(6)) THEN  
         IF (.NOT.ALLOCATED(FVSQRHOV))
     &      ALLOCATE(FVSQRHOV(NRP1,MSMAX))

         RW1 = (0.,0.)
         DO J=1,NCHI
         IF (J.EQ.1)    J1 = NCHI
         IF (J.GT.1)    J1 = J-1
         IF (J.EQ.NCHI) J2 = 1
         IF (J.LT.NCHI) J2 = J+1
         DO I=2,NR
            CTMP1 = (RJA(I+1,J)*RV1U(I+1,J)-RJA(I,J)*RV1U(I,J))/CSH(I)
            CTMP  = CTMP1*RRHOE(I,J)
            CTMP1 = (RJA(I+1,J)*RV1U(I+1,J)+RJA(I,J)*RV1U(I,J))/2.
            CTMP2 = (RRHOE(I+1,J)-RRHOE(I-1,J))/(CSM(I+1)-CSM(I-1))
            CTMP  = CTMP+CTMP1*CTMP2 
            CTMP3 = (RW4(I,J1)+RW4(I+1,J1))/2.
            CTMP1 =-DPSIDSM(I)**2*G12LM(I,J1)/RJAM(I,J1)*RRHOE(I,J1)/
     &              CTMP3*(RV1U(I,J1)+RV1U(I+1,J1))/2. + 
     &              RJAM(I,J1)*TM(I)/CTMP3*RV2U(I,J1)*RRHOE(I,J1) + 
     &              DPSIDSM(I)*RV3U(I,J1)*RRHOE(I,J1)
            CTMP3 = (RW4(I,J2)+RW4(I+1,J2))/2.
            CTMP2 =-DPSIDSM(I)**2*G12LM(I,J2)/RJAM(I,J2)*RRHOE(I,J2)/
     &              CTMP3*(RV1U(I,J2)+RV1U(I+1,J2))/2. + 
     &              RJAM(I,J2)*TM(I)/CTMP3*RV2U(I,J2)*RRHOE(I,J2) + 
     &              DPSIDSM(I)*RV3U(I,J2)*RRHOE(I,J2)
            RW1(I,J) = CTMP + (CTMP2-CTMP1)/HCHI
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         ENDDO

         FVSQRHOV = (0.,0.)
         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQRHOV(I,MS) = FVSQRHOV(I,MS) + RW1(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQRHOV = FVSQRHOV/DFLOAT(NCHI)
      ENDIF

C     COMPUTE CURL(ETA JEQ) TERM AND FOURIER TRANSFORM
C     AS A PURE EQUILIBRIUM SOURCE/SINK TERM
      IF (KVSQLIN.AND.KVSQL(7)) THEN  
         IF (.NOT.ALLOCATED(FVSQJEQ))
     &      ALLOCATE(FVSQJEQ(NRP1,MSMAX,3))
         FVSQJEQ = (0.,0.)

         RW1 = (0.,0.)
         RW2 = (0.,0.)
         RW3 = (0.,0.)
         DO J=1,NCHI
         IF (J.EQ.1)    J1 = NCHI
         IF (J.GT.1)    J1 = J-1
         IF (J.EQ.NCHI) J2 = 1
         IF (J.LT.NCHI) J2 = J+1
         DO I=2,NR
            CTMP1 =-REQ(I,J1)**2
            CTMP2 =-REQ(I,J2)**2
            RW1(I,J) = PPEQ(I)*RESIST(I)*(CTMP2-CTMP1)/HCHI

            CTMP1 =-RESIST(I)*(REQ(I,J)**2*PPEQ(I)+T(I)*TP(I)) 
            CTMP2 =-RESIST(I+1)*(REQ(I+1,J)**2*PPEQ(I+1)+T(I+1)*TP(I+1))
            RW2(I,J) =-(CTMP2-CTMP1)/CSH(I)

            CTMP1 =-RESIST(I)*TP(I)*DPSIDS(I)*G22L(I,J)/RJA(I,J) 
            CTMP2 =-RESIST(I+1)*TP(I+1)*DPSIDS(I+1)*G22L(I+1,J)/
     &              RJA(I+1,J) 
            RW3(I,J) = (CTMP2-CTMP1)/CSH(I)

            CTMP1 =-G12LM(I,J1)/RJAM(I,J1) 
            CTMP2 =-G12LM(I,J2)/RJAM(I,J2) 
            CTMP  = TPM(I)*DPSIDSM(I)*RESISM(I)
            RW3(I,J) = RW3(I,J) - CTMP*(CTMP2-CTMP1)/HCHI
         ENDDO
         RW1(NRP1,J) = RW1(NR,J)
         RW2(NRP1,J) = RW2(NR,J)
         RW3(NRP1,J) = RW3(NR,J)
         ENDDO

         DO MS=1,MSMAX
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,2)*RCHI)
         DO I=2,NRP1
            FVSQJEQ(I,MS,1) = FVSQJEQ(I,MS,1) + RW1(I,J)*CTMP1
            FVSQJEQ(I,MS,2) = FVSQJEQ(I,MS,2) + RW2(I,J)*CTMP1
            FVSQJEQ(I,MS,3) = FVSQJEQ(I,MS,3) + RW3(I,J)*CTMP1
         ENDDO
         ENDDO
         ENDDO
         FVSQJEQ = FVSQJEQ/DFLOAT(NCHI)
      ENDIF


C     FILL IN RHS WITHOUT FOURIER CONVOLUTION
      TTMP3 = -CALPHA1*CALPHA4/CALPHA2**2
      DO MS=1,MSMAX
      LXROW = (MS-1)*NXCOMP
      LYROW = (MS-1)*NYCOMP
      DO I=2,NR
         INCLUDE 'tent.inc'
         TTMP2= -(CALPHA1+(CALPHA2-CALPHA3)*CI*RNTOR*TROTI(I))*
     &        CALPHA4/CALPHA2**2
         R(KXB1+LXROW,I) = R(KXB1+LXROW,I) + TTMP2*(
     &          FF(C1,C1,C1)*X(KXB1+LXROW,I)
     &        + FFM(C1)*X(KXB1+LXROW,I-1)
     &        + FFP(C1)*X(KXB1+LXROW,I+1))
         IF (KJRER.EQ.1.OR.KJRER.EQ.2) 
     &   R(KXB1+LXROW,I) = R(KXB1+LXROW,I) - FJRECURL(I,MS,1)/CALPHA2
         IF (KVSQLIN.AND.KVSQL(3)) 
     &   R(KXB1+LXROW,I) = R(KXB1+LXROW,I) + TTMP*FVSQVXB(I,MS,1)
         IF (KVSQLIN.AND.KVSQL(7)) 
     &   R(KXB1+LXROW,I) = R(KXB1+LXROW,I) + FVSQJEQ(I,MS,1)/CALPHA2
         IF (KXX1.GT.0) 
     &   R(KXX1+LXROW,I) = R(KXX1+LXROW,I) + TTMP2*(
     &          FF(C1,C1,C1)*X(KXX1+LXROW,I)
     &        + FFM(C1)*X(KXX1+LXROW,I-1)
     &        + FFP(C1)*X(KXX1+LXROW,I+1))
         IF (KXJRE.GT.0) THEN 
            R(KXJRE+LXROW,I) = R(KXJRE+LXROW,I) + TTMP3*(
     &          FF(C1,C1,C1)*X(KXJRE+LXROW,I)
     &        + FFM(C1)*X(KXJRE+LXROW,I-1)
     &        + FFP(C1)*X(KXJRE+LXROW,I+1)) 
            IF (KJRER.GE.1.AND.KJRER.LE.4)
     &      R(KXJRE+LXROW,I) = R(KXJRE+LXROW,I) 
     &        - FJRERHS(I,MS)/CALPHA2
         ENDIF

         IF (KVSQLIN.AND.KVSQL(1)) 
     &   R(KXV1+LXROW,I) = R(KXV1+LXROW,I) - FVSQJXB(I,MS,1)/CALPHA2
         IF (KVSQLIN.AND.KVSQL(2)) 
     &   R(KXV1+LXROW,I) = R(KXV1+LXROW,I) + FVSQVGV(I,MS,1)/CALPHA2
      ENDDO

      I=NRP1
      TTMP2= -(CALPHA1+(CALPHA2-CALPHA3)*CI*RNTOR*TROTI(I))*
     &        CALPHA4/CALPHA2**2
      R(KXB1+LXROW,I) =R(KXB1+LXROW,I) + TTMP2*( 
     &      FF(C1,C1,C1)*X(KXB1+LXROW,I)
     &    + FFM(C1)*X(KXB1+LXROW,I-1))
      IF (KJRER.EQ.1.OR.KJRER.EQ.2) 
     &R(KXB1+LXROW,I) = R(KXB1+LXROW,I) - FJRECURL(I,MS,1)/CALPHA2
      IF (KXX1.GT.0) 
     &R (KXX1+LXROW,I) = R (KXX1+LXROW,I) + TTMP2*(
     &      FF(C1,C1,C1)*X(KXX1+LXROW,I)
     &    + FFM(C1)*X(KXX1+LXROW,I-1))
      IF (KXJRE.GT.0) THEN
         R(KXJRE+LXROW,I) = R(KXJRE+LXROW,I) + TTMP3*(
     &      FF(C1,C1,C1)*X(KXJRE+LXROW,I)
     &    + FFM(C1)*X(KXJRE+LXROW,I-1))
         IF (KJRER.GE.1.AND.KJRER.LE.4) 
     &   R(KXJRE+LXROW,I) = R(KXJRE+LXROW,I) 
     &   - FJRERHS(I,MS)/CALPHA2
      ENDIF

      DO I=1,NR
         TTMP2= -(CALPHA1+(CALPHA2-CALPHA3)*CI*RNTOR*TROTM(I))*
     &        CALPHA4/CALPHA2**2
         z3m  = (cs(i  )/csm(i))**iexb3
         z3p  = (cs(i+1)/csm(i))**iexb3
         INCLUDE 'tophat.inc'
         IF (IPDIVB.NE.2.OR.(IPDIVB.EQ.2.AND.ABS(RM(MS,2)).LT.0.1)) 
     $   THEN
         RY(KYB2+LYROW,I) = RY(KYB2+LYROW,I) 
     $     +GG(C1,C1,C1) * Y(KYB2+LYROW,I) * TTMP2
         IF (KJRER.EQ.1.OR.KJRER.EQ.2) 
     $   RY(KYB2+LYROW,I) = RY(KYB2+LYROW,I)-FJRECURL(I,MS,2)/CALPHA2
         IF (KVSQLIN.AND.KVSQL(3)) 
     $   RY(KYB2+LYROW,I) = RY(KYB2+LYROW,I)+TTMP*FVSQVXB(I,MS,2)
         IF (KVSQLIN.AND.KVSQL(7)) 
     $   RY(KYB2+LYROW,I) = RY(KYB2+LYROW,I)+FVSQJEQ(I,MS,2)/CALPHA2
         ENDIF

         IF (KVSQLIN.AND.KVSQL(4)) 
     $   RY(KYPR+LYROW,I) = RY(KYPR+LYROW,I) + FVSQVGP(I,MS)/CALPHA2
         IF (KVSQLIN.AND.KVSQL(5)) 
     $   RY(KYPR+LYROW,I) = RY(KYPR+LYROW,I) + FVSQPDV(I,MS)/CALPHA2
         IF (KVSQLIN.AND.KVSQL(6)) 
     $   RY(KYRHOP+LYROW,I)=RY(KYRHOP+LYROW,I) + FVSQRHOV(I,MS)/CALPHA2

         IF (IPDIVB.NE.1) THEN
           RY(KYB3+LYROW,I) = RY(KYB3+LYROW,I) +
     $     GG(C1,C1*z3m,C1*z3p)*Y(KYB3+LYROW,I)*TTMP2
           IF (KJRER.EQ.1.OR.KJRER.EQ.2) 
     $     RY(KYB3+LYROW,I)=RY(KYB3+LYROW,I)-FJRECURL(I,MS,3)/CALPHA2
           IF (KVSQLIN.AND.KVSQL(3)) 
     $     RY(KYB3+LYROW,I)=RY(KYB3+LYROW,I)+TTMP*FVSQVXB(I,MS,3)
           IF (KVSQLIN.AND.KVSQL(7)) 
     $     RY(KYB3+LYROW,I)=RY(KYB3+LYROW,I)+FVSQJEQ(I,MS,3)/CALPHA2
         ENDIF
      ENDDO

      IF ((CALPHA2-CALPHA3).NE.0..AND.CALPHA7.GT.0.) THEN
      LXCOL = (MS-1)*NXCOMP
      DO I=1,NR
      INCLUDE 'tophat.inc'

      IF (KXX1.GT.0.AND.KYX2.GT.0) THEN
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      CTMP1 = -DpsiDs(I)*TDROTI(I)
      CTMP2 = -DpsiDs(I+1)*TDROTI(I+1)
      CTMP3 = -DpsiDsM(I)*TDROTM(I)
      RY(KYX2+LYROW,I) = RY(KYX2+LYROW,I) + TTMP*(
     &      GF(ZEM*CTMP1,CTMP3)*X(KXX1+LXCOL,I) 
     &    + GF(ZEP*CTMP2,CTMP3)*X(KXX1+LXCOL,I+1)) 
      ENDIF

      IF (KVSQLIN.AND.KVSQL(1)) THEN
         RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) - FVSQJXB(I,MS,2)/CALPHA2
         RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) - FVSQJXB(I,MS,3)/CALPHA2
      ENDIF
      IF (KVSQLIN.AND.KVSQL(2)) THEN
         RY(KYV2+LYROW,I) = RY(KYV2+LYROW,I) + FVSQVGV(I,MS,2)/CALPHA2
         RY(KYV3+LYROW,I) = RY(KYV3+LYROW,I) + FVSQVGV(I,MS,3)/CALPHA2
      ENDIF

      IF (IPDIVB.NE.1) THEN
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
         RY(KYB3+LYROW,I) = RY(KYB3+LYROW,I) + TTMP*(
     &       GF(TDROTI(I)*C1,TDROTM(I)*C1)*X(KXB1+LXCOL,I)
     &    +GF(TDROTI(I+1)*C1,TDROTM(I)*C1)*X(KXB1+LXCOL,I+1))
      ELSE
         RY(KYB3+LYROW,I) = RY(KYB3+LYROW,I) + TTMP*(
     &       GF(TDROTI(I)*C1*T(I),TDROTM(I)*C1*TM(I))*X(KXB1+LXCOL,I)
     &  +GF(TDROTI(I+1)*C1*T(I+1),TDROTM(I)*C1*TM(I))*X(KXB1+LXCOL,I+1))
      ENDIF
      ENDIF
      ENDDO   
      ENDIF

      ENDDO
C
C     SECTION TO AGREE WITH BOUNDARY CONDITION IMPOSED IN BOUNDC
C
C.....FIRST FIXED-BOUNDARY CASE
C
      IF (NV.GE.2) GOTO 200
      DO 110 MS=1,MSMAX
      R(KXV1+(MS-1)*NXCOMP,NRP1) = 0.
      R(KXB1+(MS-1)*NXCOMP,NRP1) = 0.
 110  CONTINUE
C
C-----------------------------------------------------------------------
C
 200  CONTINUE
      DO 220 MS=1,MSMAX

C.....FIRST POINT ONLY

      R(KXB1+(MS-1)*NXCOMP,1) = 0.

C.....FIRST NFIT POINTS

      DO 210 I=1,NFIT
 210  R(KXV1+(MS-1)*NXCOMP,I) = 0.
 220  CONTINUE
C
C....WALL
C
      TTMP2= -CALPHA1*CALPHA4/CALPHA2**2
      IF (IWO.EQ.0) THEN
      DO MS=1,MSMAX
      IF (ABS(RM(MS,2)).GT.0.1) THEN
         RTMP = 1.0
      ELSE
         RTMP = T(NRP1)
      ENDIF
      DO J=1,NWALL
         I = NR + IWALL(J)
         R (KXB1+(MS-1)*NXCOMP,I) = X(KXB1+(MS-1)*NXCOMP,I)*RTMP*TTMP2
      ENDDO
      ENDDO
      ENDIF

      IF (IWO.EQ.1.OR.IWO.EQ.2.OR.IWO.EQ.4) THEN
      DO J=1,NWALL
      IF (IWO.EQ.1) THEN
         CALL RWALLG2(J,4)
      ELSE
         CALL RWALLG2(J,1)
      ENDIF
      I = NR + IWALL(J)

      DO MSA=1,MSMAX
      DO MSB=1,MSMAX
      MSPL =  MPLUS(MSA,NSA,MSB,NSB)
      MSMI = MMINUS(MSA,NSA,MSB,NSB)
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         RTMP = 1.0
      ELSE
         RTMP = T(NRP1)
      ENDIF
      IF (MSPL.LT.1) GOTO 260
      R (KXB1+(MSPL-1)*NXCOMP,I) = R (KXB1+(MSPL-1)*NXCOMP,I) 
     &        +TTMP2*FGRADS(MSB)*X(KXB1+(MSA-1)*NXCOMP,I)*RTMP
 260  CONTINUE
      IF (MSB.LT.2) GOTO 280
      IF (MSMI.LT.1) GOTO 280
      R (KXB1+(MSMI-1)*NXCOMP,I) = R (KXB1+(MSMI-1)*NXCOMP,I) 
     &     +CONJG(TTMP2*FGRADS(MSB))*X(KXB1+(MSA-1)*NXCOMP,I)*RTMP
 280  CONTINUE
      ENDDO
      ENDDO

      ENDDO
      ENDIF

C     DRIFT KINETIC TERMS 
      IF (INCKIN.EQ.1.AND.IPERTURB.EQ.0.AND.
     &    (V2XKEY.EQ.0.OR.V2XKEY.EQ.2))  
     &   CALL KJPDX(MD, MDY, ND, R, RY, X, Y)

C     MUTIPLY EQUATION INSIDE PLASMA BY EQFAC
C     YQLIU, 2009-10-01
      DO MS=1,MSMAX
         LXROW = (MS-1)*NXCOMP
         DO I=1,NR
            IF (KXX1.GT.0) R(LXROW+KXX1,I)=R(LXROW+KXX1,I)*EQFAC/EQFAC
            R(LXROW+KXV1,I)=R(LXROW+KXV1,I)*EQFAC/EQFAC
            R(LXROW+KXB1,I)=R(LXROW+KXB1,I)*EQFAC
         ENDDO

         LXROW = (MS-1)*NYCOMP
         DO I=1,NR
            IF (KYX2.GT.0) RY(LXROW+KYX2,I)=RY(LXROW+KYX2,I)*EQFAC/EQFAC
            RY(LXROW+KYV2,I)=RY(LXROW+KYV2,I)*EQFAC/EQFAC
            IF (KYV3.GT.0) RY(LXROW+KYV3,I)=RY(LXROW+KYV3,I)*EQFAC/EQFAC
            RY(LXROW+KYB2,I)=RY(LXROW+KYB2,I)*EQFAC/EQFAC
            RY(LXROW+KYB3,I)=RY(LXROW+KYB3,I)*EQFAC/EQFAC
         ENDDO
      ENDDO 

C     FEEDBACK TERMS 
C     YQLIU 14/04/1999
      CALL FEEDDX(MD, MDY, ND, R, RY, X, Y, XOLD, YOLD)

C------------------
C SPECIAL BC AT NV
C------------------
      IF (ABS(NCOUPL).EQ.1.OR.NCOUPL.EQ.4.OR.INCFEED.EQ.20.OR.
     &    INCFEED.EQ.21.OR.INCFEED.EQ.22) THEN

C     WRITE(*,*) 'READ BNORM01'
      IF (.NOT.ALLOCATED(BNORM)) 
     &   ALLOCATE( BNORM(MSMAX) )

      REWIND(CHOUTP)
      DO MS = 1,MSMAX
         READ(CHOUTP,*) ZEM,ZEP
         BNORM(MS) = ZEM + ZEP*CI
      ENDDO

      I = NTP1
      IF (NCOUPL.EQ.-1) I=NR+IFEED
 
      DO MS=1,MSMAX
      IF (ABS(RM(MS,2)).GT.0.1) THEN   
         R(KXB1+(MS-1)*NXCOMP,I) = BNORM(MS)
      ELSE
         R(KXB1+(MS-1)*NXCOMP,I) = BNORM(MS)/T(NRP1)
      ENDIF
      ENDDO
      ENDIF

C----------------------------------
C SPEICIAL BC FOR BACKWARD COUPLING
C----------------------------------
      IF (NCOUPL.EQ.-2.AND.INCFEED.EQ.0) THEN
         I=NTP1
         DO MSA=1,MSMAX
            LXROW=(MSA-1)*NXCOMP
            R(KXB1+LXROW,I)=0.0
            DO MSB=1,MSMAX
               LXCOL=(MSB-1)*NXCOMP
               LYCOL=(MSB-1)*NYCOMP
               IF (ABS(RM(MSB,2)).GT.0.1) THEN
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           -B11(MSA,MSB)*X(KXB1+LXCOL,I)
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           +A11(MSA,MSB)*Y(KYB2+LYCOL,I-1)
               ELSE
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           -B11(MSA,MSB)*X(KXB1+LXCOL,I)*T(NRP1)
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           +A11(MSA,MSB)*Y(KYB3+LYCOL,I-1)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (NCOUPL.EQ.-3.AND.INCFEED.EQ.0) THEN
         I=NTP1
         DO MSA=1,MSMAX
            LXROW=(MSA-1)*NXCOMP
            R(KXB1+LXROW,I)=0.0
            DO MSB=1,MSMAX
               LXCOL=(MSB-1)*NXCOMP
               LYCOL=(MSB-1)*NYCOMP
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           +A12(MSA,MSB)*Y(KYB3+LYCOL,I-1)
               IF (ABS(RM(MSB,2)).GT.0.1) THEN
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           -B12(MSA,MSB)*X(KXB1+LXCOL,I)
               ELSE
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           -B12(MSA,MSB)*X(KXB1+LXCOL,I)*T(NRP1)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (NCOUPL.EQ.-2.AND.INCFEED.GT.0) THEN
         I    = NTP1
         DO MSA=1,MSMAX
            LXROW=(MSA-1)*NXCOMP
            R(KXB1+LXROW,I)=0.0
            DO MSB=1,NCOILT
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           +(C01(MSA,MSB)+AL0*C11(MSA,MSB))*FEEDIT(MSB)
            ENDDO
         ENDDO
      ENDIF

      IF (NCOUPL.EQ.-3.AND.INCFEED.GT.0) THEN
         I    = NTP1
         DO MSA=1,MSMAX
            LXROW=(MSA-1)*NXCOMP
            R(KXB1+LXROW,I)=0.0
            DO MSB=1,NCOILT
               R(KXB1+LXROW,I)=R(KXB1+LXROW,I)
     &           +(C02(MSA,MSB)+AL0*C12(MSA,MSB))*FEEDIT(MSB)
            ENDDO
         ENDDO
      ENDIF

      IF (ALLOCATED(RB1U)) DEALLOCATE(RB1U,RB2U,RB3U,RJ1U,RJ2U,RJ3U,
     &                                RV1U,RV2U,RV3U,RRHOE,RPPE,
     &                                RW1,RW2,RW3,RW4,RW5,RW6)
      IF (ALLOCATED(RW7))     DEALLOCATE(RW7)
      IF (ALLOCATED(FJRERHS)) DEALLOCATE(FJRERHS,FJRECURL)
      IF (ALLOCATED(FVSQJXB)) DEALLOCATE(FVSQJXB)
      IF (ALLOCATED(FVSQVGV)) DEALLOCATE(FVSQVGV)
      IF (ALLOCATED(FVSQVXB)) DEALLOCATE(FVSQVXB)
      IF (ALLOCATED(FVSQVGP)) DEALLOCATE(FVSQVGP)
      IF (ALLOCATED(FVSQPDV)) DEALLOCATE(FVSQPDV)
      IF (ALLOCATED(FVSQRHOV)) DEALLOCATE(FVSQRHOV)
      IF (ALLOCATED(FVSQJEQ)) DEALLOCATE(FVSQJEQ)

      RETURN
      END
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c----- Convolution for Coeffi_New()-------- D. Liu xx.07.97 ------------
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
      SUBROUTINE FTCOEFF
c
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE RESMATM
      USE CONVOLCOFM
      USE REORBITM
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'
      INCLUDE 'comfft.inc'
c     ------------------------------------------------------B
c     Temporary variables for Doing FFTRF
c     ------------------------------------------------------E
c
      real*8,dimension(:,:),allocatable::RW1,RW2,RW3,RW4,RW5,RW6,
     $     RW7,RW8,RW9,RW10,RW11,RW12,B_2,B_2M,B2C,B2CM
C
      real*8,dimension(:,:),allocatable::Vrw1,Vrw2,Vrw3,Vrw4,Vrw5,
     $     Vrw6,Vrw7,Vrw8,Vrw9,Vrw10,Vrw11,Vrw4d,Vrw5d,Vrw9d,
     $     Vrw10d,Brw1,Brw6,Brw3,Brw8,rw4d,rw5d,
     &     rw9d,rw10d
C
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::CV2,CV3,CV4,CV5,CV6,
     &           RV1U,RV2U,RV3U

      integer    I,J
      real*8     TQH1,TQH2,RCHI,HH1,HH2
      COMPLEX*16 CTMP1,CTMP2,CTMP3,CTMP4
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'FTCOEFF'
c
c     =======================================================
c     Fourier transform arrays in equations
c
c     Need to calculate 
c		rho*Gij / Jacob            ------->   RGViGij
c		JiU                        ------->   JiBm  
c		Bi                         ------->   BiJm
c		Jacobi                     ------->   fftfnewto
c		Bi *Jacobi                 ------->   BiGvj
c		Gij                        ------->   GijJm
c		Gij                        ------->   Gijbm
c		Jacobi                     ------->   Gvi
c		dPeds                      ------->   GdPv1

c	and then Transform them into the Fourier space
c	Note: it is Bi=BNi/fbi=dpsids/Jacobi.Be cautious with
c	      equilibrium "f" coeffi.s !!!
c     =======================================================
      ALLOCATE( RW1(nrp1,nchi), RW2(nrp1,nchi), RW3(nrp1,nchi),
     $          RW4(nrp1,nchi), RW5(nrp1,nchi), RW6(nrp1,nchi),
     $          RW7(nrp1,nchi), RW8(nrp1,nchi), RW9(nrp1,nchi),
     $          RW10(nrp1,nchi),RW11(nrp1,nchi),RW12(nrp1,nchi) )
      ALLOCATE( B_2(nrp1,nchi), B_2M(nrp1,nchi) )
      ALLOCATE( B2C(nrp1,nchi), B2CM(nrp1,nchi) )
      ALLOCATE( RV1U(nrp1,nchi),RV2U(nrp1,nchi),RV3U(nrp1,nchi),
     &          CV2(nrp1,nchi), CV3(nrp1,nchi), CV4(nrp1,nchi), 
     &          CV5(nrp1,nchi), CV6(nrp1,nchi) )
C
      ALLOCATE( Vrw1(nveq1,nchi ),Vrw2(nveq1,nchi ),Vrw3(nveq1,nchi ),
     $          Vrw4(nveq1,nchi ),Vrw5(nveq1,nchi ),Vrw6(nveq1,nchi ),
     $          Vrw7(nveq1,nchi ),Vrw8(nveq1,nchi ),Vrw9(nveq1,nchi ),
     $                          Vrw10(nveq1,nchi ),Vrw11(nveq1,nchi ),
     $                          Vrw4d(nveq1,nchi ),Vrw5d(nveq1,nchi ),
     $                         Vrw9d(nveq1,nchi ),Vrw10d(nveq1,nchi ))
      ALLOCATE(             Brw1(2,nchi ),Brw6(1,nchi ),Brw3(2,nchi ),
     $                                                  Brw8(1,nchi ))
C
      ALLOCATE(                     rw4d(nrp1,nchi ),rw5d(nrp1,nchi ),
     &                              rw9d(nrp1,nchi ),rw10d(nrp1,nchi))
C
C     WRITE(*,'("FTCOEF after allocation")')
C
c     ------------------------------------------------------E
!      b2=DPsiDs(i) /RJa(i,j) 
!      b3=T(i)/Req(i,j)**2
!      j2=-DPSIDS(I)*Tp(i) /RJa(i,j) 
!      j3=-(PPeq(i) + T(i)*Tp(i)/REQ(i,j)**2)
c     ------------------------------------------------------B
      do j=1,nchi
        do i=2,nrp1
          B_2(i,j)=g22L(i,j) *DPSIDS(I)**2 /RJa(i,j)**2 +
     $        T(i)**2 /Req(i,j)**2
        end do
        B_2(1,j)=B_2(2,j)
        do i = 1,nr
          B_2M(i,j)=g22LM(i,j) *DPSIDSM(I)**2 /RJaM(i,j)**2 +
     $       TM(i)**2 /ReqM(i,j)**2
        end do
      end do

C     DB^2/DCHI
      CALL DERCHI(B_2,B2C,NRP1,NRP1)
      CALL DERCHI(B_2M,B2CM,NR,NRP1)

c     ------------------------------------------------------E
      do j=1,nchi
         do i=2,nrp1
c----  1st
            rw1(i,j)= rho(i)*RJa(i,j)*( g11L(i,j)-
     $           g12L(i,j)**2*DPSIDS(I)**2 /RJa(i,j)**2 /B_2(i,j) )    
c----  2nd  
            rw2(i,j)= rho(i)*g12L(i,j)*RJa(i,j)
     &         *T(i) /B_2(i,j)
            rw4(i,j)= rho(i)*g22L(i,j) *RJa(i,j)
     &         *G33L(i,j) /B_2(i,j)
c----  3rd  
            rw5(i,j)= rho(i) *B_2(i,j) *RJa(i,j)
            rw11(i,j)=1.0/RJA(i,j)*B2C(i,j)/2.0/B_2(i,j)**2
         end do
         do i = 1,nr
c----  1st
            rw6(i,j)= rhoM(i)*RJaM(i,j)*( g11LM(i,j)-
     $           g12LM(i,j)**2*DPSIDSM(I)**2 /RJaM(i,j)**2 /B_2M(i,j) )
c----  2nd  
            rw7(i,j)= rhoM(i) *g12LM(i,j)*RJaM(i,j)
     $           *TM(i) /B_2M(i,j)
            rw9(i,j)= rhoM(i)*g22LM(i,j) *RJaM(i,j)
     $         *G33LM(i,j) /B_2M(i,j)
c----  3rd  
            rw10(i,j)=rhoM(i) *B_2M(i,j) *RJaM(i,j)
            rw12(i,j)=1.0/RJAM(i,j)*B2CM(i,j)/2.0/B_2M(i,j)**2
         end do
      end do

      do j=1,nchi
        rw1(1,j)=0.
        rw2(1,j)=0.
        rw4(1,j)=0.
        rw5(1,j)=0.
      end do

C
      IF (.NOT. ALLOCATED(RGV1G11)) THEN
         ALLOCATE(RGV1G11(NRP1,MEDIM),RGV1G12(NRP1,MEDIM),
     &            RGV2G22(NRP1,MEDIM),TTCJB4(NRP1,MEDIM),
     &            RGV3G33(NRP1,MEDIM),RGV3G33M(NRP1,MEDIM),
     &            RGV1G11M(NRP1,MEDIM),RGV1G12M(NRP1,MEDIM),
     &            RGV2G22M(NRP1,MEDIM),TTCJB4M(NRP1,MEDIM))
      END IF
      RGV1G11=0.
      RGV1G12=0.
      RGV2G22=0.
      RGV3G33=0.
      TTCJB4 =0.
      RGV3G33M=0.
      RGV1G11M=0.
      RGV1G12M=0.
      RGV2G22M=0.
      TTCJB4M =0.
      
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW1,  RGV1G11, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G11  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW2,  RGV1G12, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G12  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW4,  RGV2G22, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2G22  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW5,  RGV3G33, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3G33  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW11,  TTCJB4, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: TTCJB4  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW6, RGV1G11M, FORWD, NRP1,  NR,     NPSTRT
     &                   ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G11M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      5,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW7, RGV1G12M, FORWD, NRP1,  NR,     NPSTRT
     &                   ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G12M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW9, RGV2G22M, FORWD, NRP1,  NR,     NPSTRT
     &                   ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2G22M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      7,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW10, RGV3G33M, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3G33M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW12, TTCJB4M, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: TTCJB4M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,    RGV1G11,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G11')
        call FFTOUTPT(RW2,    RGV1G12,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G12')
        call FFTOUTPT(RW4,    RGV2G22,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2G22')
        call FFTOUTPT(RW5,    RGV3G33,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3G33')
        call FFTOUTPT(RW11,   TTCJB4,    NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'TTCJB4')
        call FFTOUTPT(RW6,    RGV1G11M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G11M')
        call FFTOUTPT(RW7,    RGV1G12M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G12M')
        call FFTOUTPT(RW9,    RGV2G22M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2G22M')
        call FFTOUTPT(RW10,   RGV3G33M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3G33M')
        call FFTOUTPT(RW12,   TTCJB4M,   NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'TTCJB4M')
      ENDIF

      do j=1,nchi
         do i=2,nrp1
c----  1st
            rw1(i,j)= rho(i)*g22L(i,j) *RJa(i,j)**2/cs(i)
     &         *G33L(i,j) /B_2(i,j)
c----  2nd  
            rw2(i,j)= rho(i)*RJa(i,j)**2/cs(i)*( g11L(i,j)-
     $           g12L(i,j)**2*DPSIDS(I)**2 /RJa(i,j)**2 /B_2(i,j) )

         end do
         do i = 1,nr
c----  1st
            rw3(i,j)= rhoM(i)*g22LM(i,j)*RJaM(i,j)**2/csm(i)
     $           *G33LM(i,j) /B_2M(i,j)
c----  2nd  
            rw4(i,j)= rhoM(i)*RJaM(i,j)**2/csm(i)*( g11LM(i,j)-
     $           g12LM(i,j)**2*DPSIDSM(I)**2 /
     $           RJaM(i,j)**2 /B_2M(i,j) )
         end do
      end do

      do j=1,nchi
        rw1(1,j)=0.
        rw2(1,j)=0.
      end do

C
      IF (.NOT. ALLOCATED( RGV1G22)) THEN
         ALLOCATE(RGV1G22(NRP1,MEDIM),RGV2G11(NRP1,MEDIM),
     &          RGV1G22M(NRP1,MEDIM),RGV2G11M(NRP1,MEDIM))
      END IF
      RGV1G22=0.
      RGV2G11=0.
      RGV1G22M=0.
      RGV2G11M=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW1,  RGV1G22, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G22  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      9,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW2,  RGV2G11, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2G11  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     10,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW3, RGV1G22M, FORWD, NRP1,  NR,     NPSTRT
     &                   ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G22M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     11,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW4, RGV2G11M, FORWD, NRP1,  NR,     NPSTRT
     &                   ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2G11M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     12,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,    RGV1G22,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G22')
        call FFTOUTPT(RW2,    RGV2G11,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2G11')
        call FFTOUTPT(RW3,    RGV1G22M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G22M')
        call FFTOUTPT(RW4,    RGV2G11M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2G11M')
      ENDIF

C----------------------------------------------------------------
C ROTATION SECTION
C---------------------------------------------------------------
      do j=1,nchi
         do i=2,nrp1
            rw1(i,j)= rho(i)*
     $           G12L(I,J)*DPSIDS(I)*T(I)/B_2(I,J)
            rw2(i,j)= rho(i)*
     $           G33L(I,J)*G22L(I,J)*DPSIDS(I)/B_2(I,J)
         end do
         do i = 1,nr
            rw4(i,j)= rhoM(i)*
     $           G12LM(I,J)*DPSIDSM(I)*TM(I)/B_2M(I,J)
            rw5(i,j)= rhoM(i)*
     $           G33LM(I,J)*G22LM(I,J)*DPSIDSM(I)/B_2M(I,J)
         end do
      end do

      do j=1,nchi
        rw1(1,j)=0.
        rw2(1,j)=0.
      end do

C
      IF (.NOT. ALLOCATED (IDRXX)) THEN
         ALLOCATE( IDRXX(NRP1,MEDIM),IDRYX(NRP1,MEDIM),
     &             IDRXXM(NRP1,MEDIM),
     &             IDRYXM(NRP1,MEDIM) )
      END IF
      IDRXX=0.
      IDRYX=0.
      IDRXXM=0.
      IDRYXM=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW1,  IDRXX,  FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IDRXX in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     13,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW2,  IDRYX,  FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IDRYX in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     14,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW4,  IDRXXM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IDRXXM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     15,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW5,  IDRYXM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IDRYXM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     16,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,    IDRXX,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IDRXX')
        call FFTOUTPT(RW2,    IDRYX,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IDRYX')
        call FFTOUTPT(RW4,    IDRXXM,    NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IDRXXM')
        call FFTOUTPT(RW5,    IDRYXM,    NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IDRYXM')
      ENDIF
C
C
C----------------------------------------------------------------
C CENTRIFUGAL FORCE COEFFICIENTS
C----------------------------------------------------------------
      DO J=1, NCHI
         DO I=2, NRP1
            RW1(I,J)=ROT(I)**2*(RJA(I,J)**2*RDCDZ(I,J) 
     $              +DPSIDS(I)**2*G12L(I,J)*RDSDZ(I,J)/B_2(I,J))
            RW2(I,J)=ROT(I)**2*RJA(I,J)**2*T(I)*RDSDZ(I,J)/B_2(I,J)
            RW3(I,J)=ROT(I)**2*RJA(I,J)*DPSIDS(I)*RDSDZ(I,J)
         END DO
         
         DO I=1, NR
            RW4(I,J)=ROTM(I)**2*(RJAM(I,J)**2*RDCDZM(I,J) 
     $              +DPSIDSM(I)**2*G12LM(I,J)*RDSDZM(I,J)/B_2M(I,J))
            RW5(I,J)=ROTM(I)**2*RJAM(I,J)**2*TM(I)*RDSDZM(I,J)/B_2M(I,J)
            RW6(I,J)=ROTM(I)**2*RJAM(I,J)*DPSIDSM(I)*RDSDZM(I,J)
         END DO
      END DO
      IF (.NOT. ALLOCATED (DXDZDSDZ)) THEN
         ALLOCATE( DXDZDSDZ(NRP1,MEDIM),ROTDSDZ(NRP1,MEDIM),
     &             ROT2DSDZ(NRP1,MEDIM),DXDZDSDZM(NRP1,MEDIM),
     &             ROTDSDZM(NRP1,MEDIM),ROT2DSDZM(NRP1,MEDIM))
      END IF
      DXDZDSDZ=0.
      ROTDSDZ=0.
      ROT2DSDZ=0.
      DXDZDSDZM=0.  
      ROTDSDZM=0.
      ROT2DSDZM=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW1,  DXDZDSDZ,  FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DXDZDSDZ in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     17,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW2,  ROTDSDZ,   FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: ROTDSDZ in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     18,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW3,  ROT2DSDZ,  FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: ROT2DSDZ in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     19,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW4,  DXDZDSDZM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DXDZDSDZM in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     20,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW5,  ROTDSDZM,  FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: ROTDSDZM in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     21,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW6,  ROT2DSDZM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: ROT2DSDZM in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     22,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,    DXDZDSDZ,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DXDZDSDZ')
        call FFTOUTPT(RW2,    ROTDSDZ,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'ROTDSDZ')
        call FFTOUTPT(RW3,    ROT2DSDZ,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'ROT2DSDZ')
        call FFTOUTPT(RW4,    DXDZDSDZM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DXDZDDZM')
        call FFTOUTPT(RW5,    ROTDSDZM,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'ROTDSDZM')
        call FFTOUTPT(RW6,    ROT2DSDZM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'ROT2DDZM')
      ENDIF
C

      do j=1,nchi
         do i=2,nrp1
            rw1(i,j)=2.*rho(i)*RJa(i,j)**2*
     $           RBZ(I,J)/B_2(I,J)
            rw3(i,j)=2.*rho(i)*RJa(i,j)**2*
     $           RDCDZ(I,J)*T(I)/G33L(I,J)
            rw5(i,j)=2.*rho(i)*RJa(i,j)**2*RDSDZ(I,J)
         end do
         do i = 1,nr
            rw7(i,j)=2.*rhoM(i)*RJaM(i,j)**2*
     $           RBZM(I,J)/B_2M(I,J)
            rw9(i,j)=2.*rhoM(i)*RJaM(i,j)**2*
     $           RDCDZM(I,J)*TM(I)/G33LM(I,J)
            rw11(i,j)=2.*rhoM(i)*RJaM(i,j)**2*RDSDZM(I,J)
         end do
      end do

      do j=1,nchi
        rw1(1,j)=0.
        rw3(1,j)=0.
        rw5(1,j)=0.
      end do

C
      IF (.NOT. ALLOCATED (IRXY)) THEN
         ALLOCATE( IRXY(NRP1,MEDIM),
     &             IRXZ(NRP1,MEDIM),
     &             IRYZ(NRP1,MEDIM),
     &             IRXYM(NRP1,MEDIM),
     &             IRXZM(NRP1,MEDIM),
     &             IRYZM(NRP1,MEDIM) )
      END IF
      IRXY=0.
      IRXZ=0.
      IRYZ=0.
      IRXYM=0.
      IRXZM=0.
      IRYZM=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW1,  IRXY,   FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IRXY in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     23,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW3,  IRXZ,   FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IRXZ in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     24,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW5,  IRYZ,   FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IRYZ in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     25,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW7,  IRXYM,  FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IRXYM in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     26,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW9,  IRXZM,  FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IRXZM in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     27,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW11, IRYZM,  FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: IRYZM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     28,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW1,    IRXY,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IRXY')
        call FFTOUTPT(RW3,    IRXZ,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IRXZ')
        call FFTOUTPT(RW5,    IRYZ,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IRYZ')
        call FFTOUTPT(RW7,    IRXYM,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IRXYM')
        call FFTOUTPT(RW9,    IRXZM,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IRXZM')
        call FFTOUTPT(RW11,   IRYZM,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'IRYZM')
      ENDIF
C

c===Extra terms in Equ.1,3 (The ones for Equ.2 are three blocks after):
      do j=1,nchi
         do i=2,nrp1
c----  1st
            rw2(i,j)=g12L(i,j) *dPsids(i)**2 /B_2(i,j)/RJa(i,j)
            rw3(i,j)=g12L(i,j) *dPsids(i) /B_2(i,j)
     $           *T(i) /Req(i,j)**2 
         end do
         do i = 1,nr
c----  1st
            rw8(i,j)=g12LM(i,j) *dPsidsM(i)**2 /B_2M(i,j)/RJaM(i,j)
            rw9(i,j)=g12LM(i,j) *dPsidsM(i) /B_2M(i,j)
     $           *TM(i) /ReqM(i,j)**2 
         end do
         rw8(nrp1,j)=0.
         rw9(nrp1,j)=0.
      end do

      do j=1,nchi
        rw2(1,j)=0.
        rw3(1,j)=0.
      end do
C
      IF (.NOT. ALLOCATED (G12B2B2)) THEN
         ALLOCATE( G12B2B2(NRP1,MEDIM),
     &             G12B2B3(NRP1,MEDIM),
     &             G12B2B2M(NRP1,MEDIM),
     &             G12B2B3M(NRP1,MEDIM) )
      END IF
      G12B2B2=0.
      G12B2B3=0.
      G12B2B2M=0.
      G12B2B3M=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW2,  G12B2B2, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G12B2B2  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     29,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW3,  G12B2B3,FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G12B2B3  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     30,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW8, G12B2B2M, FORWD, NRP1,  NR,     NPSTRT
     &                   ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G12B2B2M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     31,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW9, G12B2B3M, FORWD, NRP1,  NR,     NPSTRT
     &                   ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G12B2B3M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     32,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW2,    G12B2B2,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G12B2B2')
        call FFTOUTPT(RW3,    G12B2B3,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G12B2B3')
        call FFTOUTPT(RW8,    G12B2B2M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G12B2B2M')
        call FFTOUTPT(RW9,    G12B2B3M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G12B2B3M')
      ENDIF

      do j=1,nchi
         do i=2,nrp1
c----  1st
           rw4(i,j)= -(PPeq(i) + T(i)*Tp(i)/REQ(i,j)**2)
     $          *RJa(i,j)
c----  2nd  
           rw3(i,j)= RJa(i,j)*(-T(i) *(PPeq(i)+T(i)*Tp(i)/REQ(i,j)**2)
     $        -dPsids(i)**2 *Tp(i) /RJa(i,j)**2
     &        *G22L(i,j) ) /B_2(i,j)
c----  3rd  
         end do
c
         do i = 1,nr
c----  1st
           rw8(i,j)= -(PPeqM(i)+TM(i)*TpM(i)/REQM(i,j)**2) 
     $          *RJaM(i,j)
c----  2nd  
           rw7(i,j)= RJaM(i,j)*( -TM(i) *(PPeqM(i) + 
     &        TM(i)*TpM(i)/REQM(i,j)**2 )
     $        -dPsidsM(i)**2 *TpM(i) /RJaM(i,j)**2 
     &        *G22LM(i,j) ) /B_2M(i,j)
c----  3rd  
         end do
      end do
c
      do j=1,nchi
        rw3(1,j)=0.
        rw4(1,j)=0.
      end do

C
      IF (.NOT. ALLOCATED (J3b1)) THEN
         ALLOCATE( J3b1(NRP1,MEDIM),J3b2(NRP1,MEDIM),
     &             J3b1M(NRP1,MEDIM),J3b2M(NRP1,MEDIM) )
      END IF
      J3b1=0.
      J3b2=0.
      J3b1M=0.
      J3b2M=0.

C     USE J3b1 AND J3b1M AS TEMPORARY VIRIABLES TO COMPUTE EQJRAPA
C     EQJPARA=SURFACE AVERAGED EQUILIBIUM PARALLEL CURRENT DENSITY      
C             WEIGHTED BY JACOBIAN      
C     EQ1_H=SURFACE AVERAGED EQUILIBIUM QUANTITY 1/H=B/B_0
C             WEIGHTED BY JACOBIAN      
      IF (.NOT. ALLOCATED(EQJPARA)) ALLOCATE(EQJPARA(NRP1,2))
      IF (.NOT. ALLOCATED(EQ1_H))   ALLOCATE(EQ1_H(NRP1,2))
      IF (.NOT. ALLOCATED(EQD1_H))  ALLOCATE(EQD1_H(NRP1,2))
      EQJPARA = 0.
      EQ1_H   = 0.
      EQD1_H  = 0.

C     EQJPARA
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J)=(-TP(I)*SQRT(B_2(I,J))-T(I)*PPEQ(I)/SQRT(B_2(I,J)))*
     &            RJA(I,J)     
      ENDDO
      DO I=1,NR
         RW2(I,J)=(-TPM(I)*SQRT(B_2M(I,J))-TM(I)*PPEQM(I)/
     &             SQRT(B_2M(I,J)))*RJAM(I,J)
      ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW1,J3b1,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,J3b1M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

      EQJPARA(2:NRP1,1) = DREAL(J3b1(2:NRP1,1)/JACOBI(2:NRP1,1))
      EQJPARA(1,1)      = EQJPARA(2,1)
      EQJPARA(1:NR,2)   = DREAL(J3b1M(1:NR,1)/JACOBM(1:NR,1))

C     EQ1_H
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J)=SQRT(B_2(I,J))*RJA(I,J)     
      ENDDO
      DO I=1,NR
         RW2(I,J)=SQRT(B_2M(I,J))*RJAM(I,J)
      ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW1,J3b1,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,J3b1M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

      EQ1_H(2:NRP1,1) = DREAL(J3b1(2:NRP1,1)/JACOBI(2:NRP1,1))
      EQ1_H(1,1)      = EQ1_H(2,1)
      EQ1_H(1:NR,2)   = DREAL(J3b1M(1:NR,1)/JACOBM(1:NR,1))

C     D(EQ1_H)/DPSI
      FFF    = EQ1_H
      CALL   DFFFDPSI(0)
      EQD1_H = DFFF

C     RESUME NORMAL WORK FOR J3b1 ETC.
      NPSTRT    =  1
      call FFTDRIVER(RW3,  J3B1,   FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: J3B1  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     33,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW4,   J3B2,   FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: J3B2  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     34,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW7,  J3B1M,  FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: J3B1M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     35,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW8,  J3B2M,  FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: J3B2M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     36,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW3,    J3B1,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'J3B1')
        call FFTOUTPT(RW4,    J3B2,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'J3B2')
        call FFTOUTPT(RW7,    J3B1M,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'J3B1M')
        call FFTOUTPT(RW8,    J3B2M,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'J3B2M')
      ENDIF
C
      do j=1,nchi
         do i=2,nrp1
c----  1st
           rw4(i,j)= 1.0/Req(i,j)**2*RJa(i,j)
c----  2nd  
           rw5(i,j)= T(i) /B_2(i,j) *RJa(i,j)
           rw6(i,j)= G22L(i,j) /B_2(i,j) *dPsids(i) 
c----  3rd  
         end do
         do i = 1,nr
c----  1st
           rw10(i,j)= 1.0/ReqM(i,j)**2*RJaM(i,j)
c----  2nd  
           rw11(i,j)= TM(i) /B_2M(i,j) *RJaM(i,j)
           rw12(i,j)=G22LM(i,j) /B_2M(i,j) *dPsidsM(i) 
c----  3rd  
         end do
      end do
c
      do j=1,nchi
        rw4(1,j)=0.
        rw5(1,j)=0.
        rw6(1,j)=0.
      end do

C
      IF (.NOT. ALLOCATED (B3j2)) THEN
         ALLOCATE( B3j2(NRP1,MEDIM),
     &             TB2(NRP1,MEDIM),G22B2B2(NRP1,MEDIM),
     &             TB2M(NRP1,MEDIM),G22B2B2M(NRP1,MEDIM),
     &             B3j2M(NRP1,MEDIM) )
      END IF
      B3j2=0.
      TB2=0.
      G22B2B2=0.
      TB2M=0.
      G22B2B2M=0.
      B3j2M=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW4,  B3j2,    FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: B3j2  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     37,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW5,  TB2,     FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: TB2   in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     38,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW6,  G22B2B2, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G22B2B2  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     39,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW10, B3j2M,   FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: B3j2M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     40,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW11, TB2M,    FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: TB2M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     41,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW12, G22B2B2M, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G22B2B2M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     42,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW4,    B3j2,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'B3j2')
        call FFTOUTPT(RW5,    TB2,       NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'TB2')
        call FFTOUTPT(RW6,    G22B2B2,   NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G22B2B2')
        call FFTOUTPT(RW10,   B3j2M,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'B3j2M')
        call FFTOUTPT(RW11,   TB2M,      NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'TB2M')
        call FFTOUTPT(RW12,   G22B2B2M,  NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G22B2B2M')
      ENDIF
C
      do j=1,nchi
         rw5(1,j) = 0.
         do i=2,nrp1
           rw5(i,j) = RJa(i,j)
         end do
      end do
      CALL DERCHI(RW5,RW5D,NRP1,NRP1)

      do j=1,nchi
         do i=1,nr
           rw10(i,j)= RJaM(i,j)
         end do
      end do
      CALL DERCHI(RW10,RW10D,NR,NRP1)

      do j=1,nchi
        rw5d(1,j)=0.
      end do

C
      IF (.NOT. ALLOCATED (B3Gv2_3dc)) THEN
         ALLOCATE( B3Gv2_3dc(NRP1,MEDIM),
     &             B3Gv2M_3dc(NRP1,MEDIM) )
      END IF
      B3Gv2_3dc=0.
      B3Gv2M_3dc=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW5d,  B3Gv2_3dc,  FORWD,  NRP1,   NRP1
     &                     ,NPSTRT,     MEDIM,  NCHI
     &                     ,KUOUT,      IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: B3Gv2_3dc in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     43,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW10d, B3Gv2M_3dc, FORWD,  NRP1,   NR
     &                     ,NPSTRT,     MEDIM,  NCHI
     &                     ,KUOUT,      IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: B3Gv2M_3dc in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     44,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW5d,   B3Gv2_3dc, NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'B3Gv23dc')
        call FFTOUTPT(RW10d,  B3Gv2M_3dc,NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'B3Gv2Mdc')
      ENDIF
C
      do j=1,nchi
         rw4(1,j) = 0.
         rw5(1,j) = 0.
         do i=2,nrp1
           rw4(i,j) = G11L(i,j)  /RJa(i,j)
           rw5(i,j) = G12L(i,j)  /RJa(i,j)
         end do
      end do
      CALL DERCHI(RW4,RW4D,NRP1,NRP1)
      CALL DERCHI(RW5,RW5D,NRP1,NRP1)
C
      do j=1,nchi
         do i=1,nr
           rw9(i,j)  = G11LM(i,j)  /RJaM(i,j)
           rw10(i,j) = G12LM(i,j)  /RJaM(i,j)
         end do
      end do
      CALL DERCHI(RW9 , RW9D,NR,NRP1)
      CALL DERCHI(RW10,RW10D,NR,NRP1)

      IF (.NOT. ALLOCATED (G11j1_3dc)) THEN
         ALLOCATE(G11j1_3dc(NRP1,MEDIM),
     &        G12j2_3dc(NRP1,MEDIM),G12j2M_3dc(NRP1,MEDIM),
     &            G11j1M_3dc(NRP1,MEDIM))
      END IF
      G11j1_3dc=0.
      G12j2_3dc=0.
      G12j2M_3dc=0.
      G11j1M_3dc=0.
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW4d,  G11j1_3dc,  FORWD,  NRP1,   NRP1
     &                     ,NPSTRT,     MEDIM,  NCHI
     &                     ,KUOUT,      IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G11j1_3dc in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     45,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW5d,  G12j2_3dc,  FORWD,  NRP1,   NRP1
     &                     ,NPSTRT,     MEDIM,  NCHI
     &                     ,KUOUT,      IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G12j2_3dc in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     46,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW9d,  G11j1M_3dc, FORWD,  NRP1,   NR
     &                     ,NPSTRT,     MEDIM,  NCHI
     &                     ,KUOUT,      IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G11j1M_3dc in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     47,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW10d, G12j2M_3dc, FORWD,  NRP1,   NR
     &                     ,NPSTRT,     MEDIM,  NCHI
     &                     ,KUOUT,      IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: G12j2M_3dc in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     48,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW4d,   G11j1_3dc, NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G11j13dc')
        call FFTOUTPT(RW5d,   G12j2_3dc, NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G12j23dc')
        call FFTOUTPT(RW9d,   G11j1M_3dc,NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G11j1Mdc')
        call FFTOUTPT(RW10d,  G12j2M_3dc,NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'G12j2Mdc')
      ENDIF
C
c     ------------------------------------------------------B Vac

      if(nveq .ge.1) then
  
        IF (NVACJ.EQ.1) THEN
          CALL DERCHI(VRJa,Vrw4,NVEQ1,NVEQ1)
          do j=1,nchi
             do i=1,nveq1
               Vrw4d(i,j) = Vrw4(i,j)*VRR(i,j)/sqrt(abs(VRG11L(i,j)*
     &                      VRG22L(i,j)-VRG12L(i,j)**2))
             end do
          end do

          CALL DERCHI(VRJaM,Vrw9,NVEQ,NVEQ1)
          do j=1,nchi
             do i=1,nveq
               Vrw9d(i,j) = Vrw9(i,j)*VRRM(i,j)/sqrt(abs(VRG11LM(i,j)*
     &                      VRG22LM(i,j)-VRG12LM(i,j)**2))
             end do
          end do

          do j=1,nchi
             do i=1,nveq1
                Vrw1(i,j)= VRG33L(i,j)
                Vrw2(i,j)= VRG12L(i,j)
                Vrw3(i,j)= VRG22L(i,j)
                Vrw4(i,j)= VRG11L(i,j)
                Vrw5(i,j)= VRJA(i,j)
             end do
             do i = 1,nveq
                Vrw6(i,j)= VRG33LM(i,j)
                Vrw7(i,j)= VRG12LM(i,j)
                Vrw8(i,j)= VRG22LM(i,j)
                Vrw9(i,j)= VRG11LM(i,j)
                Vrw10(i,j)=VRJAM(i,j)
             end do
          end do

        ELSE

          Vrw4d = 0.0
          Vrw9d = 0.0

          do j=1,nchi
            do i=1,nveq1
              Vrw1(i,j)= VRG33L(i,j)/VRJA(i,j)
              Vrw2(i,j)= VRG12L(i,j)/VRJA(i,j)
              Vrw3(i,j)= VRG22L(i,j)/VRJA(i,j)
              Vrw4(i,j)= VRG11L(i,j)/VRJA(i,j)
              Vrw5(i,j)= 1.0
            end do
            do i = 1,nveq
              Vrw6(i,j)= VRG33LM(i,j)/VRJAM(i,j)
              Vrw7(i,j)= VRG12LM(i,j)/VRJAM(i,j)
              Vrw8(i,j)= VRG22LM(i,j)/VRJAM(i,j)
              Vrw9(i,j)= VRG11LM(i,j)/VRJAM(i,j)
              Vrw10(i,j)=1.0
           end do
         end do

        ENDIF

        IF (.NOT. ALLOCATED (VG33JS)) THEN
           ALLOCATE( VG33JS(nveq1,MEDIM),
     $               VG11L(nveq1,MEDIM),
     $               VG33L(nveq1,MEDIM),
     $               VG12L(nveq1,MEDIM),
     $               VG22L(nveq1,MEDIM),
     $               VJAC(nveq1,MEDIM),
     $               VG33JC(nveq1,MEDIM) )
        END IF
        VG33JS=0.
        VG11L=0.
        VG33L=0.
        VG12L=0.
        VG22L=0.
        VJAC =0.
        VG33JC=0.                          
C
        IF (.NOT. ALLOCATED (VG33JSM)) THEN
           ALLOCATE( VG33JSM(nveq1,MEDIM),
     $               VG11LM(nveq1,MEDIM),
     $               VG33LM(nveq1,MEDIM),
     $               VG12LM(nveq1,MEDIM),
     $               VG22LM(nveq1,MEDIM),
     $               VJACM(nveq1,MEDIM),
     $               VG33JCM(nveq1,MEDIM) )
        END IF
        VG33JSM=0.
        VG11LM=0.
        VG33LM=0.
        VG12LM=0.
        VG22LM=0.
        VJACM =0.
        VG33JCM=0.
C
        NVSTRT    =  1
        call FFTDRIVER(VRW1,   VG33L,      FORWD,  NVEQ1,  NVEQ1
     &                        ,NVSTRT,     MEDIM,  NCHI
     &                        ,KUOUT,      IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG33L  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     49,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW2,   VG12L,      FORWD,  NVEQ1,  NVEQ1
     &                        ,NVSTRT,     MEDIM,  NCHI
     &                        ,KUOUT,      IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG12L  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     50,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW3,   VG22L,      FORWD,  NVEQ1,  NVEQ1
     &                        ,NVSTRT,     MEDIM,  NCHI
     &                        ,KUOUT,      IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG22L  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     51,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW4,   VG11L,      FORWD,  NVEQ1,  NVEQ1
     &                        ,NVSTRT,     MEDIM,  NCHI
     &                        ,KUOUT,      IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG11L  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     52,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW5,   VJAC,        FORWD,  NVEQ1,  NVEQ1
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VJAC in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     53,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW4d,  VG33JC,      FORWD,  NVEQ1,  NVEQ1
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG33JC  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     54,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW6,   VG33LM,      FORWD,  NVEQ1,  NVEQ
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG33LM  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     55,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW7,   VG12LM,      FORWD,  NVEQ1,  NVEQ
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG12LM  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     56,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW8,   VG22LM,      FORWD,  NVEQ1,  NVEQ
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG22LM  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     57,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW9,   VG11LM,      FORWD,  NVEQ1,  NVEQ
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG11LM  in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     58,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW10,  VJACM,       FORWD,  NVEQ1,  NVEQ
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VJACM in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     59,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        call FFTDRIVER(VRW9d,  VG33JCM,     FORWD,  NVEQ1,  NVEQ
     &                        ,NVSTRT,      MEDIM,  NCHI
     &                        ,KUOUT,       IERSUB, IERPLC, IERR)
        if(IERR .NE. 0) THEN
           write(MESSAGE,*) 'Error: VG33JCM in ',IERSUB
           call ABORTRUN
     &         (SUBNAM,     60,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
        endif
C
        IF(KUFFTP .GT. 0) THEN
          call FFTOUTPT(VRW1,   VG33L,     NVEQ1,   NVEQ1,   NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG33L')
          call FFTOUTPT(VRW2,   VG12L,     NVEQ1,   NVEQ1,   NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG12L')
          call FFTOUTPT(VRW3,   VG22L,     NVEQ1,   NVEQ1,   NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG22L')
          call FFTOUTPT(VRW4,   VG11L,     NVEQ1,   NVEQ1,   NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG11L')
          call FFTOUTPT(VRW5,   VJAC,      NVEQ1,   NVEQ1,   NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VJAC')
          call FFTOUTPT(VRW4d,  VG33JC,    NVEQ1,   NVEQ1,   NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG33JC')
          call FFTOUTPT(VRW6,   VG33LM,    NVEQ1,   NVEQ,    NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG33LM')
          call FFTOUTPT(VRW7,   VG12LM,    NVEQ1,   NVEQ,    NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG12LM')
          call FFTOUTPT(VRW8,   VG22LM,    NVEQ1,   NVEQ,    NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG22LM')
          call FFTOUTPT(VRW9,   VG11LM,    NVEQ1,   NVEQ,    NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG11LM')
          call FFTOUTPT(VRW10,  VJACM,     NVEQ1,   NVEQ,    NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VJACM')
          call FFTOUTPT(VRW9d,  VG33JCM,   NVEQ1,   NVEQ,    NVSTRT
     &                         ,MEDIM,     NCHI,    KUFFTP, 'VG33JCM')
        ENDIF 
      endif
C
c     ------------------------------------------------------B Vac
      if(nveq .gt.1) then
 
        IF (NVACJ.EQ.1) THEN
          DO I=1,NVEQ1
          IF (I.GT.1)     TQH1 = (VCS(I)-VCS(I-1))*0.5     
          IF (I.LT.NVEQ1) TQH2 = (VCS(I+1)-VCS(I))*0.5
          DO J=1,NCHI
             IF (I.EQ.1) THEN
                Vrw4(I,J) = (VRJAM(1,J)-VRJA(1,J))/TQH2
             ELSEIF (I.EQ.NVEQ1) THEN
                Vrw4(I,J) = (VRJA(I,J)-VRJAM(I-1,J))/TQH1
             ELSE
                Vrw4(I,J) = (TQH1**2*(VRJAM(I,J)-VRJA(I,J)) + 
     &                       TQH2**2*(VRJA(I,J)-VRJAM(I-1,J)))/
     &                      (TQH1*TQH2*(TQH1+TQH2))
             ENDIF
          ENDDO
          ENDDO

          do j=1,nchi
             do i=1,nveq1
               Vrw4d(i,j) = Vrw4(i,j)*VRR(i,j)/sqrt(abs(VRG11L(i,j)*
     &                      VRG22L(i,j)-VRG12L(i,j)**2))
             end do
          end do

          do j=1,nchi
             do i=1,nveq
               TQH2       = VCS(i+1)-VCS(i)
               Vrw9d(i,j) = (VRJA(i+1,j)-VRJA(i,j))/TQH2*
     *                      VRRM(i,j)/sqrt(abs(VRG11LM(i,j)*
     &                      VRG22LM(i,j)-VRG12LM(i,j)**2))
             end do
          end do
C
C
          NVSTRT    =  1
          call FFTDRIVER(VRW4d,  VG33JS,     FORWD,  NVEQ1,  NVEQ1
     &                          ,NVSTRT,     MEDIM,  NCHI
     &                          ,KUOUT,      IERSUB, IERPLC, IERR)
          if(IERR .NE. 0) THEN
             write(MESSAGE,*) 'Error: VG33JS  in ',IERSUB
             call ABORTRUN
     &         (SUBNAM,     61,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
          endif
C
          call FFTDRIVER(VRW9d,  VG33JSM,    FORWD,  NVEQ1,  NVEQ
     &                          ,NVSTRT,     MEDIM,  NCHI
     &                          ,KUOUT,      IERSUB, IERPLC, IERR)
          if(IERR .NE. 0) THEN
             write(MESSAGE,*) 'Error: VG33JS  in ',IERSUB
             call ABORTRUN
     &         (SUBNAM,     62,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
          endif
C
          IF(KUFFTP .GT. 0) THEN
            call FFTOUTPT(VRW4d,  VG33JS,    NVEQ1,   NVEQ1,   NVSTRT
     &                           ,MEDIM,     NCHI,    KUFFTP, 'VG33JS')
            call FFTOUTPT(VRW9d,  VG33JSM,   NVEQ1,   NVEQ,    NVSTRT
     &                           ,MEDIM,     NCHI,    KUFFTP, 'VG33JSM')
          ENDIF 
        ENDIF
      endif

      IF (KJRER.GT.0) THEN
      IF (.NOT.ALLOCATED(FRE1B)) THEN
         ALLOCATE( FRE1B(NRP1,MEDIM),    FRE1BM(NRP1,MEDIM),
     &             FREG22JB(NRP1,MEDIM), FREG22JBM(NRP1,MEDIM),  
     &             FREG12JB(NRP1,MEDIM), FREG12JBM(NRP1,MEDIM) ) 
      ENDIF
      FRE1B     = 0.0
      FRE1BM    = 0.0
      FREG22JB  = 0.0
      FREG22JBM = 0.0
      FREG12JB  = 0.0
      FREG12JBM = 0.0

C     1/B
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J)=1./SQRT(B_2(I,J))
      ENDDO
      DO I=1,NR
         RW2(I,J)=1./SQRT(B_2M(I,J))
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRE1B,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRE1BM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     G22L/(JB)
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J)=G22L(I,J)/RJA(I,J)/SQRT(B_2(I,J))
      ENDDO
      DO I=1,NR
         RW2(I,J)=G22LM(I,J)/RJAM(I,J)/SQRT(B_2M(I,J))
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREG22JB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREG22JBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     G12L/(JB)
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J)=G12L(I,J)/RJA(I,J)/SQRT(B_2(I,J))
      ENDDO
      DO I=1,NR
         RW2(I,J)=G12LM(I,J)/RJAM(I,J)/SQRT(B_2M(I,J))
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREG12JB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREG12JBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

      ENDIF

      IF (KJRER.EQ.5.OR.KJRER.EQ.6) THEN
      IF (.NOT.ALLOCATED(FRE1JB)) THEN
         ALLOCATE( FRE1JB(NRP1,MEDIM),  FRE1JBM(NRP1,MEDIM),
     &             FRE1JCB(NRP1,MEDIM), FRE1JCBM(NRP1,MEDIM),  
     &             FRE1R2B(NRP1,MEDIM), FRE1R2BM(NRP1,MEDIM), 
     &             FREJPBS(NRP1,MEDIM), FREJPBSM(NRP1,MEDIM), 
     &             FREJPBC(NRP1,MEDIM), FREJPBCM(NRP1,MEDIM),  
     &             FREJPB(NRP1,MEDIM),  FREJPBM(NRP1,MEDIM),  
     &             FREJPBG12(NRP1,MEDIM),FREJPBG12M(NRP1,MEDIM),  
     &             FREJPBG11(NRP1,MEDIM),FREJPBG11M(NRP1,MEDIM),  
     &             FREJPBG22(NRP1,MEDIM),FREJPBG22M(NRP1,MEDIM) ) 
      ENDIF
      FRE1JB    = 0.0
      FRE1JBM   = 0.0
      FRE1JCB   = 0.0
      FRE1JCBM  = 0.0
      FRE1R2B   = 0.0
      FRE1R2BM  = 0.0
      FREJPBS   = 0.0
      FREJPBSM  = 0.0
      FREJPBC   = 0.0
      FREJPBCM  = 0.0
      FREJPB    = 0.0
      FREJPBM   = 0.0
      FREJPBG12 = 0.0
      FREJPBG12M= 0.0
      FREJPBG11 = 0.0
      FREJPBG11M= 0.0
      FREJPBG22 = 0.0
      FREJPBG22M= 0.0

C     1/(JB)
      DO J=1,NCHI
      DO I=2,NRP1
         RW1(I,J)=1./RJA(I,J)/SQRT(B_2(I,J))
      ENDDO
      RW1(1,J)=RW1(2,J)
      DO I=1,NR
         RW2(I,J)=1./RJAM(I,J)/SQRT(B_2M(I,J))
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRE1JB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRE1JBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     (1/J)(D/DCHI)(1/B)
      DO J=1,NCHI
      DO I=2,NRP1
         RW1(I,J)=-B2C(I,J)/RJA(I,J)/B_2(I,J)**1.5/2.
      ENDDO
      RW1(1,J)=RW1(2,J)
      DO I=1,NR
         RW2(I,J)=-B2CM(I,J)/RJAM(I,J)/B_2M(I,J)**1.5/2.
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRE1JCB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRE1JCBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     1/(R^2B)
      DO J=1,NCHI
      DO I=1,NRP1
         RW1(I,J)=1./REQ(I,J)**2/SQRT(B_2(I,J))
      ENDDO
      DO I=1,NR
         RW2(I,J)=1./REQM(I,J)**2/SQRT(B_2M(I,J))
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRE1R2B,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRE1R2BM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     (1/J)(D/DS)(J_PARA/B)
      DO J=1,NCHI
      DO I=1,NRP1
         RW3(I,J)=-TP(I)-T(I)*PPEQ(I)/B_2(I,J)
      ENDDO
      DO I=1,NR
         RW4(I,J)=-TPM(I)-TM(I)*PPEQM(I)/B_2M(I,J)
      ENDDO
      ENDDO
 
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RW1(I,J)=(H1/H2*RW4(I,J)-H2/H1*RW4(I-1,J))/(H1+H2) -
     &               (H1-H2)*RW3(I,J)/H1/H2
            RW1(I,J)=RW1(I,J)/RJA(I,J)
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RW1(I,J)=(RW3(I-1,J)+3*RW3(I,J)-4*RW4(I-1,J))/H1/RJA(I,J)
         RW1(1,J)=RW1(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RW2(I,J)=(RW3(I+1,J)-RW3(I,J))/H1/RJAM(I,J)
         ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJPBS,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJPBSM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     (1/J)(D/DCHI)(J_PARA/B)
      CALL DERCHI(RW3,RW5,NRP1,NRP1)
      CALL DERCHI(RW4,RW6,NR,NRP1)

      DO J=1,NCHI
      DO I=2,NRP1
         RW1(I,J)=RW5(I,J)/RJA(I,J)
      ENDDO
      RW1(1,J)=RW1(2,J)
      DO I=1,NR
         RW2(I,J)=RW6(I,J)/RJAM(I,J)
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJPBC,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJPBCM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     J_PARA/B
      NPSTRT    =  1
      call FFTDRIVER(RW3,FREJPB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW4,FREJPBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     (J_PARA*G12)/(JB)
      DO J=1,NCHI
      DO I=2,NRP1
         RW1(I,J)=RW3(I,J)*G12L(I,J)/RJA(I,J)
      ENDDO
      RW1(1,J)=RW1(2,J)
      DO I=1,NR
         RW2(I,J)=RW4(I,J)*G12LM(I,J)/RJAM(I,J)
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJPBG12,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJPBG12M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     (J_PARA*G11)/(JB)
      DO J=1,NCHI
      DO I=2,NRP1
         RW1(I,J)=RW3(I,J)*G11L(I,J)/RJA(I,J)
      ENDDO
      RW1(1,J)=RW1(2,J)
      DO I=1,NR
         RW2(I,J)=RW4(I,J)*G11LM(I,J)/RJAM(I,J)
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJPBG11,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJPBG11M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     (J_PARA*G22)/(JB)
      DO J=1,NCHI
      DO I=2,NRP1
         RW1(I,J)=RW3(I,J)*G22L(I,J)/RJA(I,J)
      ENDDO
      RW1(1,J)=RW1(2,J)
      DO I=1,NR
         RW2(I,J)=RW4(I,J)*G22LM(I,J)/RJAM(I,J)
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJPBG22,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJPBG22M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

      ENDIF

      IF (KJRER.EQ.6) THEN
      IF (.NOT.ALLOCATED(FRENRE)) 
     &   ALLOCATE( FRENRE(NRP1,MEDIM),     FRENREM(NRP1,MEDIM),
     &          FRENREJ(NRP1,MEDIM),       FRENREJM(NRP1,MEDIM),
     &          FRENDB(NRP1,MEDIM),        FRENDBM(NRP1,MEDIM),
     &          FRENDB2(NRP1,MEDIM),       FRENDB2M(NRP1,MEDIM),
     &          FRENJDB2(NRP1,MEDIM),      FRENJDB2M(NRP1,MEDIM),
     &          FRENDB2R2(NRP1,MEDIM),     FRENDB2R2M(NRP1,MEDIM),
     &          FRENG12DB2R2(NRP1,MEDIM),  FRENG12DB2R2M(NRP1,MEDIM),
     &          FRENG22DB2(NRP1,MEDIM),    FRENG22DB2M(NRP1,MEDIM),
     &          FRENG11DJB2(NRP1,MEDIM),   FRENG11DJB2M(NRP1,MEDIM),
     &          FRENG12DJB2(NRP1,MEDIM),   FRENG12DJB2M(NRP1,MEDIM),
     &          FRENR2G22DJ2B2(NRP1,MEDIM),FRENR2G22DJ2B2M(NRP1,MEDIM),
     &          FRENR2G12DJ2B2(NRP1,MEDIM),FRENR2G12DJ2B2M(NRP1,MEDIM),
     &          FREJDBR2(NRP1,MEDIM),      FREJDBR2M(NRP1,MEDIM), 
     &          FRE1J(NRP1,MEDIM),         FRE1JM(NRP1,MEDIM),
     &          FRENREJB(NRP1,MEDIM),      FRENREJBM(NRP1,MEDIM),
     &          FREJNREDB(NRP1,MEDIM),     FREJNREDBM(NRP1,MEDIM),
     &          FREJNREDB3(NRP1,MEDIM),    FREJNREDB3M(NRP1,MEDIM),
     &          FREJNREDB3R2(NRP1,MEDIM),  FREJNREDB3R2M(NRP1,MEDIM),
     &          FREJNREG11DJB3(NRP1,MEDIM),FREJNREG11DJB3M(NRP1,MEDIM),
     &          FREJNREG12DJB3(NRP1,MEDIM),FREJNREG12DJB3M(NRP1,MEDIM),
     &          FREJNREG22DJB3(NRP1,MEDIM),FREJNREG22DJB3M(NRP1,MEDIM) )

      FRENRE          = 0.0
      FRENREM         = 0.0
      FRENREJ         = 0.0
      FRENREJM        = 0.0
      FRENDB          = 0.0
      FRENDB2         = 0.0
      FRENDB2M        = 0.0
      FRENJDB2        = 0.0
      FRENJDB2M       = 0.0
      FRENDB2R2       = 0.0
      FRENDB2R2M      = 0.0
      FRENG12DB2R2    = 0.0
      FRENG12DB2R2M   = 0.0
      FRENG22DB2      = 0.0
      FRENG22DB2M     = 0.0
      FRENG11DJB2     = 0.0
      FRENG11DJB2M    = 0.0
      FRENG12DJB2     = 0.0
      FRENG12DJB2M    = 0.0
      FRENR2G22DJ2B2  = 0.0
      FRENR2G22DJ2B2M = 0.0
      FRENR2G12DJ2B2  = 0.0
      FRENR2G12DJ2B2M = 0.0
      FREJDBR2        = 0.0
      FREJDBR2M       = 0.0
      FRE1J           = 0.0
      FRE1JM          = 0.0
      FRENREJB        = 0.0
      FRENREJBM       = 0.0
      FREJNREDB       = 0.0
      FREJNREDBM      = 0.0
      FREJNREDB3      = 0.0
      FREJNREDB3M     = 0.0
      FREJNREDB3R2    = 0.0
      FREJNREDB3R2M   = 0.0
      FREJNREG11DJB3  = 0.0
      FREJNREG11DJB3M = 0.0
      FREJNREG12DJB3  = 0.0
      FREJNREG12DJB3M = 0.0
      FREJNREG22DJB3  = 0.0
      FREJNREG22DJB3M = 0.0

C     J_PARA
      DO J=1,NCHI
      DO I=1,NRP1
         RW3(I,J)=-TP(I)*SQRT(B_2(I,J))-T(I)*PPEQ(I)/SQRT(B_2(I,J))
      ENDDO
      DO I=1,NR
         RW4(I,J)=-TPM(I)*SQRT(B_2M(I,J))-TM(I)*PPEQM(I)/SQRT(B_2M(I,J))
      ENDDO
      ENDDO

C     N_RE   
      RW3 =-1./(OMEGACI0*C_VA)*RW3
      RW4 =-1./(OMEGACI0*C_VA)*RW4
 
      NPSTRT    =  1
      call FFTDRIVER(RW3,FRENRE,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW4,FRENREM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE*J   
      RW1 = RW3*RJA
      RW2 = RW4*RJAM
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENREJ,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENREJM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE/B   
      RW1 = RW3/SQRT(B_2)
      RW2 = RW4/SQRT(B_2M)
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENDB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENDBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE/B^2   
      RW1 = RW3/B_2
      RW2 = RW4/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENDB2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENDB2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE*J/B^2   
      RW1 = RW3*RJA/B_2
      RW2 = RW4*RJAM/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENJDB2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENJDB2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE/(B^2*R^2)   
      RW1 = RW3/B_2/REQ**2
      RW2 = RW4/B_2M/REQM**2
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENDB2R2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENDB2R2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE*G12/(B^2*R^2)   
      RW1 = RW3*G12L/B_2/REQ**2
      RW2 = RW4*G12LM/B_2M/REQM**2
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENG12DB2R2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENG12DB2R2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE*G22/(B^2)   
      RW1 = RW3*G22L/B_2
      RW2 = RW4*G22LM/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENG22DB2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENG22DB2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

      RJA(1,:) = RJA(2,:)

C     N_RE*G11/(J*B^2)   
      RW1 = RW3*G11L/RJA/B_2
      RW2 = RW4*G11LM/RJAM/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENG11DJB2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENG11DJB2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE*G12/(J*B^2)   
      RW1 = RW3*G12L/RJA/B_2
      RW2 = RW4*G12LM/RJAM/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENG12DJB2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENG12DJB2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE*R^2*G22/(J^2*B^2)   
      RW1 = RW3*REQ**2*G22L/RJA**2/B_2
      RW2 = RW4*REQM**2*G22LM/RJAM**2/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENR2G22DJ2B2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENR2G22DJ2B2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     N_RE*R^2*G12/(J^2*B^2)   
      RW1 = RW3*REQ**2*G12L/RJA**2/B_2
      RW2 = RW4*REQM**2*G12LM/RJAM**2/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENR2G12DJ2B2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENR2G12DJ2B2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     J/(B*R^2)   
      RW1 = RJA/SQRT(B_2)/REQ**2
      RW2 = RJAM/SQRT(B_2M)/REQM**2
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJDBR2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJDBR2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     1/J   
      RW1 = 1.0/RJA
      RW2 = 1.0/RJAM
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRE1J,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRE1JM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     JB*NRE 
      RW1 = RW3*RJA*SQRT(B_2)
      RW2 = RW4*RJAM*SQRT(B_2M)
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FRENREJB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FRENREJBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     JRE*NRE/B
      RW1 =-OMEGACI0*C_VA*RW3**2/SQRT(B_2)
      RW2 =-OMEGACI0*C_VA*RW4**2/SQRT(B_2M)
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJNREDB,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJNREDBM,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     JRE*NRE/B^3, NOTE CHANGE OF RW3-4 HERE-AFTER
      RW3 = RW1/B_2
      RW4 = RW2/B_2M
 
      NPSTRT    =  1
      call FFTDRIVER(RW3,FREJNREDB3,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW4,FREJNREDB3M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     JRE*NRE/(B^3*R^2)
      RW1 = RW3/REQ**2
      RW2 = RW4/REQM**2
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJNREDB3R2,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJNREDB3R2M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     JRE*NRE*G11/(JB^3)
      RW1 = RW3*G11L/RJA
      RW2 = RW4*G11LM/RJAM
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJNREG11DJB3,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJNREG11DJB3M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     JRE*NRE*G12/(JB^3)
      RW1 = RW3*G12L/RJA
      RW2 = RW4*G12LM/RJAM
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJNREG12DJB3,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJNREG12DJB3M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

C     JRE*NRE*G22/(JB^3)
      RW1 = RW3*G22L/RJA
      RW2 = RW4*G22LM/RJAM
 
      NPSTRT    =  1
      call FFTDRIVER(RW1,FREJNREG22DJB3,FORWD,NRP1,NRP1,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)
      call FFTDRIVER(RW2,FREJNREG22DJB3M,FORWD,NRP1,NR,NPSTRT,
     &               MEDIM,NCHI,KUOUT,IERSUB,IERPLC,IERR)

      ENDIF

 1101 format(8(E11.4,1X))
C
      IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10) THEN
         IF (.NOT.ALLOCATED(FHATV))
     &      ALLOCATE( FHATV(NRP1,MEDIM,6),  FHATVM(NRP1,MEDIM,6),
     &                FHATVP(NRP1,MEDIM,6), FHATVPM(NRP1,MEDIM,6),
     &                FHATVN(NRP1,MEDIM,6), FHATVNM(NRP1,MEDIM,6))

         FHATV   = (0.,0.)
         FHATVP  = (0.,0.)
         FHATVN  = (0.,0.)
         FHATVM  = (0.,0.)
         FHATVPM = (0.,0.)
         FHATVNM = (0.,0.)

C        PREPARE FOURIER COEFFICIENTS FOR HAT_V^1-6 
C        FOR IMPLICIT TREATMENT OF QUADARTIC TERM CURL(VXB) FOR N=0
         RV1U = (0.,0.)
         RV2U = (0.,0.)
         RV3U = (0.,0.)

         DO J=1,NCHI
         RCHI = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         DO I=1,NR
            DO MS=1,MSMAX
               CTMP1 = EXP(CI*RM(MS,2)*RCHI)
               RV1U(I,J)  = RV1U(I,J) + V1U(I,MS)*CTMP1
               RV2U(I,J)  = RV2U(I,J) + V2U(I,MS)*CTMP1
               RV3U(I,J)  = RV3U(I,J) + V3U(I,MS)*CTMP1
            ENDDO
         ENDDO
         ENDDO

         RV1U(NRP1,1:NCHI) = RV1U(NR,1:NCHI)
         RV2U(NRP1,1:NCHI) = RV2U(NR,1:NCHI)
         RV3U(NRP1,1:NCHI) = RV3U(NR,1:NCHI)
 
         CV2 = (0.,0.)
         CV3 = (0.,0.)
         CV4 = (0.,0.)
         CV5 = (0.,0.)
         CV6 = (0.,0.)

         DO J=1,NCHI
         DO I=2,NR
            CTMP1=(B_2(I,J)+B_2(I+1,J))/2.
            CTMP2=(RV1U(I,J)+RV1U(I+1,J))/2.
            CTMP3=RV3U(I,J)-DPSIDSM(I)*G12LM(I,J)*CTMP2/RJAM(I,J)/CTMP1
            CTMP4=RV2U(I,J)/CTMP1
            CV2(I,J) = CTMP3*DPSIDSM(I)/RJAM(I,J)+CTMP4*TM(I)
            CV3(I,J) = CTMP3*TM(I)/REQM(I,J)**2-
     &                 CTMP4*DPSIDSM(I)*G22LM(I,J)/RJAM(I,J)
            CV4(I,J) = RV1U(I,J)*RJA(I,J)/G22L(I,J)
            CV5(I,J) = CTMP2*G12LM(I,J)+RV2U(I,J)
            CV6(I,J) = RV1U(I,J)*RJA(I,J)/REQ(I,J)**2
         ENDDO
         CV2(NRP1,J) = CV2(NR,J)
         CV3(NRP1,J) = CV3(NR,J)
         CV4(NRP1,J) = CV4(NR,J)
         CV5(NRP1,J) = CV5(NR,J)
         CV6(NRP1,J) = CV6(NR,J)
         ENDDO

C        COMPUTE FOURIER HARMONICS OF RV1U,CV2-6  
C        IN THE RANGE OF EQUILIBRIUM HARMONICS [-MEDIM,MEDIM]
C        NOTE THAT FHATVN* IS COMPLEX CONJUGATE OF FOURIER HARMONICS
C        THE ABOVE ARE NEEDED TO BE COMPATIBLE WITH COEFFI()
C        SIMILAR TREATMENT NEED TO BE FOLLOWED FOR ALL QUADRATIC TERMS
         DO MS=1,MEDIM
         DO J=1,NCHI
         RCHI  = 2.*PI*DFLOAT(J-1)/DFLOAT(NCHI)
         CTMP1 = EXP(-CI*RM(MS,1)*RCHI)
         CTMP2 = CONJG(CTMP1)
         DO I=1,NR
            FHATVP(I,MS,1)  = FHATVP(I,MS,1)  + RV1U(I,J)*CTMP1
            FHATVPM(I,MS,2) = FHATVPM(I,MS,2) + CV2(I,J)*CTMP1
            FHATVPM(I,MS,3) = FHATVPM(I,MS,3) + CV3(I,J)*CTMP1
            FHATVP(I,MS,4)  = FHATVP(I,MS,4)  + CV4(I,J)*CTMP1
            FHATVPM(I,MS,5) = FHATVPM(I,MS,5) + CV5(I,J)*CTMP1
            FHATVP(I,MS,6)  = FHATVP(I,MS,6)  + CV6(I,J)*CTMP1
            FHATVN(I,MS,1)  = FHATVN(I,MS,1)  + CONJG(RV1U(I,J)*CTMP2)
            FHATVNM(I,MS,2) = FHATVNM(I,MS,2) + CONJG(CV2(I,J)*CTMP2)
            FHATVNM(I,MS,3) = FHATVNM(I,MS,3) + CONJG(CV3(I,J)*CTMP2)
            FHATVN(I,MS,4)  = FHATVN(I,MS,4)  + CONJG(CV4(I,J)*CTMP2)
            FHATVNM(I,MS,5) = FHATVNM(I,MS,5) + CONJG(CV5(I,J)*CTMP2)
            FHATVN(I,MS,6)  = FHATVN(I,MS,6)  + CONJG(CV6(I,J)*CTMP2)
         ENDDO
         ENDDO
         FHATVP(NRP1,MS,1)  = FHATVP(NR,MS,1)
         FHATVPM(NRP1,MS,2:3) = FHATVPM(NR,MS,2:3)
         FHATVP(NRP1,MS,4)  = FHATVP(NR,MS,4)
         FHATVPM(NRP1,MS,5) = FHATVPM(NR,MS,5)
         FHATVP(NRP1,MS,6)  = FHATVP(NR,MS,6)
         FHATVN(NRP1,MS,1)  = FHATVN(NR,MS,1)
         FHATVNM(NRP1,MS,2:3) = FHATVNM(NR,MS,2:3)
         FHATVN(NRP1,MS,4)  = FHATVN(NR,MS,4)
         FHATVNM(NRP1,MS,5) = FHATVNM(NR,MS,5)
         FHATVN(NRP1,MS,6)  = FHATVN(NR,MS,6)
         ENDDO
         FHATVP  = FHATVP/DFLOAT(NCHI)
         FHATVPM = FHATVPM/DFLOAT(NCHI)
         FHATVN  = FHATVN/DFLOAT(NCHI)
         FHATVNM = FHATVNM/DFLOAT(NCHI)

         DO MS=1,MEDIM
         DO I=1,NR
            FHATVPM(I,MS,1) = 0.5*(FHATVP(I,MS,1)+FHATVP(I+1,MS,1))
            FHATVPM(I,MS,4) = 0.5*(FHATVP(I,MS,4)+FHATVP(I+1,MS,4))
            FHATVPM(I,MS,6) = 0.5*(FHATVP(I,MS,6)+FHATVP(I+1,MS,6))
            FHATVNM(I,MS,1) = 0.5*(FHATVN(I,MS,1)+FHATVN(I+1,MS,1))
            FHATVNM(I,MS,4) = 0.5*(FHATVN(I,MS,4)+FHATVN(I+1,MS,4))
            FHATVNM(I,MS,6) = 0.5*(FHATVN(I,MS,5)+FHATVN(I+1,MS,6))
         ENDDO
         FHATVPM(NRP1,MS,1) = FHATVP(NRP1,MS,1)
         FHATVPM(NRP1,MS,4) = FHATVP(NRP1,MS,4)
         FHATVPM(NRP1,MS,6) = FHATVP(NRP1,MS,6)
         FHATVNM(NRP1,MS,1) = FHATVN(NRP1,MS,1)
         FHATVNM(NRP1,MS,4) = FHATVN(NRP1,MS,4)
         FHATVNM(NRP1,MS,6) = FHATVN(NRP1,MS,6)
         DO I=2,NR
            HH1  = CSH(I)/(CSH(I-1)+CSH(I))
            HH2  = CSH(I-1)/(CSH(I-1)+CSH(I))
            FHATVP(I,MS,2) = FHATVPM(I-1,MS,2)*HH1+FHATVPM(I,MS,2)*HH2
            FHATVP(I,MS,3) = FHATVPM(I-1,MS,3)*HH1+FHATVPM(I,MS,3)*HH2
            FHATVP(I,MS,5) = FHATVPM(I-1,MS,5)*HH1+FHATVPM(I,MS,5)*HH2
            FHATVN(I,MS,2) = FHATVNM(I-1,MS,2)*HH1+FHATVNM(I,MS,2)*HH2
            FHATVN(I,MS,3) = FHATVNM(I-1,MS,3)*HH1+FHATVNM(I,MS,3)*HH2
            FHATVN(I,MS,5) = FHATVNM(I-1,MS,3)*HH1+FHATVNM(I,MS,5)*HH2
         ENDDO
         FHATVP(NRP1,MS,2:3) = FHATVPM(NR,MS,2:3)
         FHATVP(NRP1,MS,5)   = FHATVPM(NR,MS,5)
         FHATVN(NRP1,MS,2:3) = FHATVNM(NR,MS,2:3)
         FHATVN(NRP1,MS,5)   = FHATVNM(NR,MS,5)
         ENDDO
      ENDIF

      deallocate( RW1,RW2,RW3,RW4,RW5,RW6,
     $            RW7,RW8,RW9,RW10,RW11,RW12,B_2,B_2M,B2C,B2CM,
     $            Vrw1,Vrw2,Vrw3,Vrw4,Vrw5,
     $            Vrw6,Vrw7,Vrw8,Vrw9,Vrw10,Vrw11,Vrw4d,Vrw5d,Vrw9d,
     $            Vrw10d,Brw1,Brw6,Brw3,Brw8,rw4d,rw5d,
     &            rw9d,rw10d,RV1U,RV2U,RV3U,CV2,CV3,CV4,CV5,CV6 )

C
      IF ((IDIAMB.GT.0.AND.IDIAMB.LT.5).OR.IDIAMTI.NE.0.OR.
     &     ABS(TTCPARA0)+ABS(TTCPERP0).GT.0.) CALL DFTCOEFF

      IF (NPROFRC.GT.0) CALL CFTCOEFF

      IF (NPROFRP.GT.0) CALL PFTCOEFF

C     write(*,*) 'BEFORE KFTCOEFF'
      CALL KFTCOEFF

      RETURN
      End   !!!{---FTCOEFF}

C=======================================================================
C PREPARE EQUILIBRIUM FOURIER COEFFICIENTS FOR SUBROUTINE DIAMAGNETIC  = 
C TERM, SHOULD BE CALLED FROM WITHIN FTCOEFF                           =
C YQL, 06-2009                                                         =
C=======================================================================
      SUBROUTINE DFTCOEFF

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      IMPLICIT NONE

      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER NPSTRT
      INTEGER I,J
      REAL*8  H1,H2
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW0,RW0M,RW1,RW1M,
     &                                   B_2,B_2M,B2S,B2SM,B2C,B2CM,      
     &                                   JB2,JB2M,JB2S,JB2SM,JB2C,JB2CM       

      ALLOCATE( RW0(NRP1,NCHI), RW0M(NRP1,NCHI), 
     &          RW1(NRP1,NCHI), RW1M(NRP1,NCHI), 
     &          B_2(NRP1,NCHI), B_2M(NRP1,NCHI),
     &          B2S(NRP1,NCHI), B2SM(NRP1,NCHI),
     &          B2C(NRP1,NCHI), B2CM(NRP1,NCHI), 
     &          JB2(NRP1,NCHI), JB2M(NRP1,NCHI),
     &          JB2S(NRP1,NCHI), JB2SM(NRP1,NCHI),
     &          JB2C(NRP1,NCHI), JB2CM(NRP1,NCHI) )
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'DFTCOEFF'
C
C     B^2
      DO J=1,NCHI
        DO I=2,NRP1
          B_2(I,J)=G22L(I,J)*DPSIDS(I)**2/RJA(I,J)**2 +
     &             T(I)**2/REQ(I,J)**2
        ENDDO 
        B_2(1,J)=T(1)**2/REQ(1,J)**2
        DO I=1,NR
          B_2M(I,J)=G22LM(I,J)*DPSIDSM(I)**2/RJAM(I,J)**2 +
     &             TM(I)**2/REQM(I,J)**2
        ENDDO 
      ENDDO 

C     DB^2/DCHI
      CALL DERCHI(B_2,B2C,NRP1,NRP1)
      CALL DERCHI(B_2M,B2CM,NR,NRP1)

C     DB^2/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            B2S(I,J)=(H1/H2*B_2M(I,J)-H2/H1*B_2M(I-1,J))/(H1+H2) -
     &               (H1-H2)*B_2(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         B2S(I,J)=(B_2(I-1,J)+3*B_2(I,J)-4*B_2M(I-1,J))/H1
         B2S(1,J)=B2S(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            B2SM(I,J)=(B_2(I+1,J)-B_2(I,J))/H1
         ENDDO
      ENDDO

C     1/(JB^2)
      DO J=1,NCHI
        DO I=2,NRP1
          JB2(I,J)=1./(RJA(I,J)*B_2(I,J))
        ENDDO 
        JB2(1,J)=JB2(2,J)
        DO I=1,NR
          JB2M(I,J)=1./(RJAM(I,J)*B_2M(I,J))
        ENDDO 
      ENDDO 

C     D(JB2)/DCHI
      CALL DERCHI(JB2,JB2C,NRP1,NRP1)
      CALL DERCHI(JB2M,JB2CM,NR,NRP1)

C     D(JB2)/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            JB2S(I,J)=(H1/H2*JB2M(I,J)-H2/H1*JB2M(I-1,J))/(H1+H2) -
     &               (H1-H2)*JB2(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         JB2S(I,J)=(JB2(I-1,J)+3*JB2(I,J)-4*JB2M(I-1,J))/H1
         JB2S(1,J)=JB2S(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            JB2SM(I,J)=(JB2(I+1,J)-JB2(I,J))/H1
         ENDDO
      ENDDO

      IF (.NOT.ALLOCATED(DJB2)) THEN
         ALLOCATE( DJB2(NRP1,MEDIM),     DJB2M(NRP1,MEDIM),
     &             DG22J2B2(NRP1,MEDIM), DG22J2B2M(NRP1,MEDIM),  
     &             DG12J2B2(NRP1,MEDIM), DG12J2B2M(NRP1,MEDIM),  
     &             DGJBA(NRP1,MEDIM),    DGJBAM(NRP1,MEDIM),  
     &             DB22(NRP1,MEDIM),     DB22M(NRP1,MEDIM),  
     &             DJFB(NRP1,MEDIM),     DJFBM(NRP1,MEDIM),  
     &             DGJBB(NRP1,MEDIM),    DGJBBM(NRP1,MEDIM),  
     &             DR2B2(NRP1,MEDIM),    DR2B2M(NRP1,MEDIM) )
      ENDIF
      DJB2      = 0.0
      DJB2M     = 0.0
      DG22J2B2  = 0.0
      DG22J2B2M = 0.0
      DG12J2B2  = 0.0
      DG12J2B2M = 0.0
      DGJBA     = 0.0
      DGJBAM    = 0.0
      DB22      = 0.0
      DB22M     = 0.0
      DJFB      = 0.0
      DJFBM     = 0.0
      DGJBB     = 0.0
      DGJBBM    = 0.0
      DR2B2     = 0.0
      DR2B2M    = 0.0

C     1/(JB^2)
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=JB2(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=JB2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,   DJB2,   FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DJB2 in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M,  DJB2M,  FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DJB2M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DJB2,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DJB2')
        call FFTOUTPT(RW0M,   DJB2M,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DJB2M')
      ENDIF
C
C     G22/(J^2B^2)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=G22L(I,J)/RJA(I,J)*JB2(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=G22LM(I,J)/RJAM(I,J)*JB2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DG22J2B2,  FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG22J2B2 in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M, DG22J2B2M, FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG22J2B2M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DG22J2B2,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DG22J2B2')
        call FFTOUTPT(RW0M,   DG22J2B2M, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DG22J2BM')
      ENDIF
C
C     G12/(J^2B^2)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=G12L(I,J)/RJA(I,J)*JB2(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=G12LM(I,J)/RJAM(I,J)*JB2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DG12J2B2,  FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG12J2B2 in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      5,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M, DG12J2B2M, FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DG12J2B2M in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DG12J2B2,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DG12J2B2')
        call FFTOUTPT(RW0M,   DG12J2B2M, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DG12J2BM')
      ENDIF
C
C     Psi'*G22/J*D/DS(1/J/B^2) - Psi'*G12/J*D/DC(1/J/B^2) + J^phi/B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=DPSIDS(I)*G22L(I,J)/RJA(I,J)*JB2S(I,J)
     &              - DPSIDS(I)*G12L(I,J)/RJA(I,J)*JB2C(I,J)
     &              - (PPEQ(I)+T(I)*TP(I)/REQ(I,J)**2)/B_2(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=DPSIDSM(I)*G22LM(I,J)/RJAM(I,J)*JB2SM(I,J)
     &              - DPSIDSM(I)*G12LM(I,J)/RJAM(I,J)*JB2CM(I,J)
     &              - (PPEQM(I)+TM(I)*TPM(I)/REQM(I,J)**2)/B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DGJBA,     FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DGJBA in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      7,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M, DGJBAM,    FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DGJBAM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DGJBA,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DGJBA')
        call FFTOUTPT(RW0M,   DGJBAM   , NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DGJBAM')
      ENDIF
C
C     DB^2/DC/B^4
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=B2C(I,J)/B_2(I,J)**2
        ENDDO
        RW0(1,J) = RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=B2CM(I,J)/B_2M(I,J)**2
        ENDDO
      ENDDO
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DB22,      FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DB22  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      9,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M, DB22M,     FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DBB22M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     10,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DB22,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DB22')
        call FFTOUTPT(RW0M,   DB22M,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DB22M')
      ENDIF
C
C     (JJ^CHI + F/B^2DB^2/DS)/B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=(T(I)/B_2(I,J)*B2S(I,J) - TP(I)*DPSIDS(I))/B_2(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=(TM(I)/B_2M(I,J)*B2SM(I,J) 
     &               - TPM(I)*DPSIDSM(I))/B_2M(I,J)
        ENDDO
      ENDDO
C
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DJFB,      FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)

      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DJFB  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     11,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M, DJFBM ,    FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)

      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DJFBM in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     12,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DJFB,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DJFB')
        call FFTOUTPT(RW0M,   DJFBM,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DJFBM')
      ENDIF
C
C     (JJ^phi - Psi'*G22/(JB^2)*DB^2/DS + Psi'*G12/(JB^2)*DB^2/DC)/B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=(-(PPEQ(I)+T(I)*TP(I)/REQ(I,J)**2)*RJA(I,J)
     &              - DPSIDS(I)*G22L(I,J)*B2S(I,J)*JB2(I,J)
     &              + DPSIDS(I)*G12L(I,J)*B2C(I,J)*JB2(I,J))
     &              /B_2(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=(-(PPEQM(I)+TM(I)*TPM(I)/REQM(I,J)**2)*RJAM(I,J)
     &              - DPSIDSM(I)*G22LM(I,J)*B2SM(I,J)*JB2M(I,J)
     &              + DPSIDSM(I)*G12LM(I,J)*B2CM(I,J)*JB2M(I,J))
     &              /B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DGJBB,     FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
C
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DGJBB in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     13,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M, DGJBBM,    FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
C
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DGJBBM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     14,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DGJBB,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DGJBB')
        call FFTOUTPT(RW0M,   DGJBBM,    NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DGJBBM')
      ENDIF
C
C     1/(R^2*B^2)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=1./REQ(I,J)**2/B_2(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=1./REQM(I,J)**2/B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DR2B2,     FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DR2B2 in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     15,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER(RW0M, DR2B2M,    FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DR2B2M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     16,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DR2B2,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DR2B2')
        call FFTOUTPT(RW0M,   DR2B2M,    NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DR2B2M')
      ENDIF
C
      IF (.NOT.ALLOCATED(DB2CDJB4)) THEN
         ALLOCATE( DB2CDJB4(NRP1,MEDIM), DB2CDJB4M(NRP1,MEDIM) )
      ENDIF
      DB2CDJB4        = 0.0
      DB2CDJB4M       = 0.0
         
C     DB^2/DC/(JB^4)
      DO J=1,NCHI
      DO I=2,NRP1
         RW0(I,J)=B2C(I,J)*JB2(I,J)/B_2(I,J)
      ENDDO
      RW0(1,J) = RW0(2,J)
      DO I=1,NR
         RW0M(I,J)=B2CM(I,J)*JB2M(I,J)/B_2M(I,J)
      ENDDO
      ENDDO
 
      NPSTRT    =  1
      call FFTDRIVER(RW0,  DB2CDJB4,  FORWD,  NRP1,   NRP1
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
 
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DB2CDJB4 in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     13,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
 
      call FFTDRIVER(RW0M, DB2CDJB4M, FORWD,  NRP1,   NR
     &                    ,NPSTRT,    MEDIM,  NCHI
     &                    ,KUOUT,     IERSUB, IERPLC, IERR)
 
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: DB2CDJB4M  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,     14,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
 
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    DB2CDJB4,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DGJBB')
        call FFTOUTPT(RW0M,   DB2CDJB4M, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'DGJBBM')
      ENDIF
 
      DEALLOCATE( RW0,RW0M,B_2,B_2M,B2S,B2SM,B2C,B2CM,
     &                     JB2,JB2M,JB2S,JB2SM,JB2C,JB2CM )
      RETURN
      END

C=======================================================================
C PREPARE EQUILIBRIUM FOURIER COEFFICIENTS FOR SUBROUTINE POLOIDALLY   = 
C VARYING TOROIDAL ROTATION, SHOULD BE CALLED FROM WITHIN FTCOEFF      =
C YQL&LL, 12-2016                                                      =
C=======================================================================
      SUBROUTINE CFTCOEFF

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      IMPLICIT NONE

      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER NPSTRT
      INTEGER I,J
      REAL*8  H1,H2
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW0,RW0M,RW1,RW1M,
     &                                   B_2,B_2M,RCS,RCSM,RCC,RCCM     

      ALLOCATE( RW0(NRP1,NCHI), RW0M(NRP1,NCHI),
     &          RW1(NRP1,NCHI), RW1M(NRP1,NCHI), 
     &          B_2(NRP1,NCHI), B_2M(NRP1,NCHI),
     &          RCS(NRP1,NCHI), RCSM(NRP1,NCHI),
     &          RCC(NRP1,NCHI), RCCM(NRP1,NCHI)) 
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'CFTCOEFF'
C
C     B^2
      DO J=1,NCHI
        DO I=2,NRP1
          B_2(I,J)=G22L(I,J)*DPSIDS(I)**2/RJA(I,J)**2 +
     &             T(I)**2/REQ(I,J)**2
        ENDDO 
        B_2(1,J)=T(1)**2/REQ(1,J)**2
        DO I=1,NR
          B_2M(I,J)=G22LM(I,J)*DPSIDSM(I)**2/RJAM(I,J)**2 +
     &             TM(I)**2/REQM(I,J)**2
        ENDDO 
      ENDDO 

C     DROTC/DCHI = RCC
      CALL DERCHI(ROTC,RCC,NRP1,NRP1)
      CALL DERCHI(ROTCM,RCCM,NR,NRP1)

C     DROTC/DS = RCS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RCS(I,J)=(H1/H2*ROTCM(I,J)-H2/H1*ROTCM(I-1,J))/(H1+H2) -
     &               (H1-H2)*ROTC(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RCS(I,J)=(ROTC(I-1,J)+3*ROTC(I,J)-4*ROTCM(I-1,J))/H1
         RCS(1,J)=RCS(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RCSM(I,J)=(ROTC(I+1,J)-ROTC(I,J))/H1
         ENDDO
      ENDDO

C     LLI
      IF (.NOT.ALLOCATED(RGV3G33C)) THEN
         ALLOCATE( 
     &            RGV3G33C(NRP1,MEDIM),RGV3G33CM(NRP1,MEDIM),
     &            RGROTCJ (NRP1,MEDIM),RGROTCJM (NRP1,MEDIM),
     &            RGROTC  (NRP1,MEDIM),RGROTCM  (NRP1,MEDIM),
     &            RGV2G22C(NRP1,MEDIM),RGV2G22CM(NRP1,MEDIM),
     &            RGV1G12C(NRP1,MEDIM),RGV1G12CM(NRP1,MEDIM),
     &            RGV1G11C(NRP1,MEDIM),RGV1G11CM(NRP1,MEDIM))
      ENDIF

      RGV3G33C  = 0.
      RGV3G33CM = 0.
      RGROTCJ   = 0.
      RGROTCJM  = 0.
      RGV2G22C  = 0.
      RGV2G22CM = 0.
      RGV1G12C  = 0.
      RGV1G12CM = 0.
      RGV1G11C  = 0.
      RGV1G11CM = 0.

C     RHO*J*B^2*ROTC
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RHO(I)*B_2(I,J)*RJA(I,J)*ROTC(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*B_2M(I,J)*RJAM(I,J)*ROTCM(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV3G33C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3G33C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV3G33CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3G33CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV3G33C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3G33C')
        call FFTOUTPT(RW0M,   RGV3G33CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3G33CM')
      ENDIF

C     ROTC*J
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RJA(I,J)*ROTC(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RJAM(I,J)*ROTCM(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGROTCJ, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROTCJ  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGROTCJM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROTCJM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGROTCJ,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGROTCJ')
        call FFTOUTPT(RW0M,   RGROTCJM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGROTCJM')
      ENDIF

C     ROTC
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=ROTC(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=ROTCM(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGROTC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROTC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGROTCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROTCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGROTC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGROTC')
        call FFTOUTPT(RW0M,   RGROTCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGROTCM')
      ENDIF
C
C     LLI
C     RHO*G22L*G33L*J*ROTC/B^2
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J) =RHO(I)*G22L(I,J)*G33L(I,J)*RJA(I,J)*ROTC(I,J) 
     &               /B_2(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*G22LM(I,J)*G33LM(I,J)*RJAM(I,J)*ROTCM(I,J)
     &               /B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2G22C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2G22C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2G22CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2G22CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2G22C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2G22C')
        call FFTOUTPT(RW0M,   RGV2G22CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2G22CM')
      ENDIF

C     LLI
C     RHO*G12L*J*F*ROTC/B^2
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J) =RHO(I)*G12L(I,J)*RJA(I,J)*T(I)*ROTC(I,J) 
     &               /B_2(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*G12LM(I,J)*RJAM(I,J)*TM(I)*ROTCM(I,J)
     &               /B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1G12C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G12C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1G12CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G12CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1G12C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G12C')
        call FFTOUTPT(RW0M,   RGV1G12CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G12CM')
      ENDIF

C     RHO*ROTC*J*(G11L-(G12L*DPSIDS/B/J)**2)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =RHO(I)*ROTC(I,J)*RJA(I,J)*(G11L(I,J)
     &               -G12L(I,J)**2*DPSIDS(I)**2
     &               /B_2(I,J)/RJA(I,J)**2)
        ENDDO
        RW0(1,J) = RW0(2,J) 
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*ROTCM(I,J)*RJAM(I,J)*(G11LM(I,J)
     &               -G12LM(I,J)**2*DPSIDSM(I)**2
     &               /B_2M(I,J)/RJAM(I,J)**2)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1G11C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G11C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1G11CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1G11CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1G11C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G11C')
        call FFTOUTPT(RW0M,   RGV1G11CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1G11CM')
      ENDIF

C    LLI
      IF (.NOT.ALLOCATED(RGIRXYC)) THEN
         ALLOCATE( 
     &            RGIRXYC(NRP1,MEDIM),RGIRXYCM(NRP1,MEDIM),
     &            RGIRXZC(NRP1,MEDIM),RGIRXZCM(NRP1,MEDIM),
     &            RGIRYZC(NRP1,MEDIM),RGIRYZCM(NRP1,MEDIM))
      ENDIF

      RGIRXYC  = 0.
      RGIRXYCM = 0.
      RGIRXZC  = 0.
      RGIRXZCM = 0.
      RGIRYZC  = 0.
      RGIRYZCM = 0.

C    2*rho*ROTC*J^2*RBZ/B^2  
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J) =2.*RHO(I)*ROTC(I,J)*RJA(I,J)**2*
     &               RBZ(I,J)/B_2(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=2.*RHOM(I)*ROTCM(I,J)*RJAM(I,J)**2*
     &               RBZM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGIRXYC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGIRXYC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGIRXYCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGIRXYCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGIRXYC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGIRXYC')
        call FFTOUTPT(RW0M,   RGIRXYCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGIRXYCM')
      ENDIF
C
C    2*rho*ROTC*J^2*RDCDZ*F/R^2  
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J) =2.*RHO(I)*ROTC(I,J)*RJA(I,J)**2*
     &               RDCDZ(I,J)*T(I)/G33L(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=2.*RHOM(I)*ROTCM(I,J)*RJAM(I,J)**2*
     &               RDCDZM(I,J)*TM(I)/G33LM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGIRXZC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGIRXZC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGIRXZCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGIRXZCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGIRXZC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGIRXZC')
        call FFTOUTPT(RW0M,   RGIRXZCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGIRXZCM')
      ENDIF
C
C    2*rho*ROTC*J^2*RDSDZ  
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J) =2.*RHO(I)*ROTC(I,J)*RJA(I,J)**2*RDSDZ(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=2.*RHOM(I)*ROTCM(I,J)*RJAM(I,J)**2*RDSDZM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGIRYZC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGIRYZC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGIRYZCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGIRYZCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGIRYZC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGIRYZC')
        call FFTOUTPT(RW0M,   RGIRYZCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGIRYZCM')
      ENDIF
C
C    LLI
      IF (.NOT.ALLOCATED(RGDXDZDSDZC)) THEN
         ALLOCATE( 
     &            RGDXDZDSDZC(NRP1,MEDIM),RGDXDZDSDZCM(NRP1,MEDIM),
     &            RGROTDSDZC(NRP1,MEDIM), RGROTDSDZCM(NRP1,MEDIM),
     &            RGROT2DSDZC(NRP1,MEDIM),RGROT2DSDZCM(NRP1,MEDIM))
      ENDIF

      RGDXDZDSDZC  = 0.
      RGDXDZDSDZCM = 0.
      RGROTDSDZC   = 0.
      RGROTDSDZCM  = 0.
      RGROT2DSDZC  = 0.
      RGROT2DSDZCM = 0.

C    [(ROT+ROTC)^2-ROT^2]*[J^2*RDCDZ+psi'^2*G12L*RDSDZ/B^2]  
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=((ROT(I)+ROTC(I,J))**2-ROT(I)**2)*
     &               (RJA(I,J)**2*RDCDZ(I,J)
     &                  +DPSIDS(I)**2*G12L(I,J)*RDSDZ(I,J)/B_2(I,J))
        ENDDO
        DO I=1,NR
           RW0M(I,J)=((ROTM(I)+ROTCM(I,J))**2-ROTM(I)**2)*
     &               (RJAM(I,J)**2*RDCDZM(I,J)
     &                  +DPSIDSM(I)**2*G12LM(I,J)*RDSDZM(I,J)/B_2M(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGDXDZDSDZC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGDXDZDSDZC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGDXDZDSDZCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGDXDZDSDZCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGDXDZDSDZC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGDXDZDSDZC')
        call FFTOUTPT(RW0M,   RGDXDZDSDZCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGDXDZDSDZCM')
      ENDIF


C    [(ROT+ROTC)^2-ROT^2]*[J^2*F*RDCDZ/B^2]  
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=((ROT(I)+ROTC(I,J))**2-ROT(I)**2)*
     &              RJA(I,J)**2*T(I)*RDSDZ(I,J)/B_2(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=((ROTM(I)+ROTCM(I,J))**2-ROTM(I)**2)*
     &               RJAM(I,J)**2*TM(I)*RDSDZM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGROTDSDZC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROTDSDZC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGROTDSDZCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROTDSDZCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGROTDSDZC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGROTDSDZC')
        call FFTOUTPT(RW0M,   RGROTDSDZCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGROTDSDZCM')
      ENDIF

C    [(ROT+ROTC)^2-ROT^2]*[J*psi'*RDCDZ]  
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=((ROT(I)+ROTC(I,J))**2-ROT(I)**2)
     &               *RJA(I,J)*DPSIDS(I)*RDSDZ(I,J)
           RW0M(I,J)=((ROTM(I)+ROTCM(I,J))**2-ROTM(I)**2)
     &               *RJAM(I,J)*DPSIDSM(I)*RDSDZM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGROT2DSDZC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROT2DSDZC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGROT2DSDZCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGROT2DSDZCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGROT2DSDZC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGROT2DSDZC')
        call FFTOUTPT(RW0M,   RGROT2DSDZCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGROT2DSDZCM')
      ENDIF

C    LLI
      IF (.NOT.ALLOCATED(RGDRCDSC)) THEN
         ALLOCATE( 
     &            RGDRCDSC(NRP1,MEDIM),RGDRCDSCM(NRP1,MEDIM),
     &            RGDRCDCC(NRP1,MEDIM),RGDRCDCCM(NRP1,MEDIM))
      ENDIF

      RGDRCDSC  = 0.
      RGDRCDSCM = 0.
      RGDRCDCC  = 0.
      RGDRCDCCM = 0.

C   DROTC/DS = RCS 
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RCS(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RCSM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGDRCDSC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGDRCDSC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGDRCDSCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGDRCDSCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGDRCDSC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGDRCDSC')
        call FFTOUTPT(RW0M,   RGDRCDSCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGDRCDSCM')
      ENDIF

C     DROTC/DCHI = RCC
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RCC(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RCCM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGDRCDCC, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGDRCDCC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGDRCDCCM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGDRCDCCM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGDRCDCC,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGDRCDCC')
        call FFTOUTPT(RW0M,   RGDRCDCCM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGDRCDCCM')
      ENDIF

C    LLI
      IF (.NOT.ALLOCATED(RGX2RX1C)) THEN
         ALLOCATE( 
     &            RGX2RX1C(NRP1,MEDIM),RGX2RX1CM(NRP1,MEDIM),
     &            RGX3RX1C(NRP1,MEDIM),RGX3RX1CM(NRP1,MEDIM),
     &            RGX2RX2C(NRP1,MEDIM),RGX2RX2CM(NRP1,MEDIM),
     &            RGX3RX2C(NRP1,MEDIM),RGX3RX2CM(NRP1,MEDIM),
     &            RGX2RX3C(NRP1,MEDIM),RGX2RX3CM(NRP1,MEDIM))
      ENDIF

      RGX2RX1C  = 0.
      RGX2RX1CM = 0.
      RGX3RX1C  = 0.
      RGX3RX1CM = 0.
      RGX2RX2C  = 0.
      RGX2RX2CM = 0.
      RGX3RX2C  = 0.
      RGX3RX2CM = 0.
      RGX2RX3C  = 0.
      RGX2RX3CM = 0.

C     RCS-psi'^2*G12L*RCC/J^2/B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW1(I,J) =RCS(I,J)-DPSIDS(I)**2
     &                 *G12L(I,J)*RCC(I,J)/RJA(I,J)**2/B_2(I,J)
        ENDDO
        RW1(1,J) = RW1(2,J) 
        DO I=1,NR
           RW1M(I,J)=RCSM(I,J)-DPSIDSM(I)**2
     &                 *G12LM(I,J)*RCCM(I,J)/RJAM(I,J)**2/B_2M(I,J)
        ENDDO
      ENDDO

C     psi*RW1
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=DPSIDS(I)*RW1(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=DPSIDSM(I)*RW1M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGX2RX1C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX1C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGX2RX1CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX1CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGX2RX1C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGX2RX1C')
        call FFTOUTPT(RW0M,   RGX2RX1CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGX2RX1CM')
      ENDIF

C     (J*F/B^2)*RW1
      DO J=1,NCHI
        DO I=1,NRP1
          RW0(I,J)=(RJA(I,J)*T(I)/B_2(I,J))*RW1(I,J)
        ENDDO
        DO I=1,NR
          RW0M(I,J)=(RJAM(I,J)*TM(I)/B_2M(I,J))*RW1M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGX3RX1C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX3RX1C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGX3RX1CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX3RX1CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGX3RX1C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGX3RX1C')
        call FFTOUTPT(RW0M,   RGX3RX1CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGX3RX1CM')
      ENDIF

C     psi*RCC*F/B^2
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=DPSIDS(I)*RCC(I,J)*T(I)/B_2(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=DPSIDSM(I)*RCCM(I,J)*TM(I)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGX2RX2C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX2C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGX2RX2CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX2CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGX2RX2C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGX2RX2C')
        call FFTOUTPT(RW0M,   RGX2RX2CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGX2RX2CM')
      ENDIF

C     J*F/B^2*RCC*F/B^2
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RJA(I,J)*T(I)**2*RCC(I,J)/B_2(I,J)**2
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RJAM(I,J)*TM(I)**2*RCCM(I,J)/B_2M(I,J)**2
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGX3RX2C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX3RX2C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGX3RX2CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX3RX2CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGX3RX2C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGX3RX2C')
        call FFTOUTPT(RW0M,   RGX3RX2CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGX3RX2CM')
      ENDIF

C     psi^2*RCC/J
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=DPSIDS(I)**2*RCC(I,J)/RJA(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J) 
        DO I=1,NR
           RW0M(I,J)=DPSIDSM(I)**2*RCCM(I,J)/RJAM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGX2RX3C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX3C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGX2RX3CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX3CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGX2RX3C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGX2RX3C')
        call FFTOUTPT(RW0M,   RGX2RX3CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGX2RX3CM')
      ENDIF

C    LLI
      IF (.NOT.ALLOCATED(RGV1RV1C)) THEN
         ALLOCATE( 
     &            RGV1RV1C(NRP1,MEDIM),RGV1RV1CM(NRP1,MEDIM),
     &            RGV2RV1C(NRP1,MEDIM),RGV2RV1CM(NRP1,MEDIM),
     &            RGV3RV1C(NRP1,MEDIM),RGV3RV1CM(NRP1,MEDIM),
     &            RGV1RV2C(NRP1,MEDIM),RGV1RV2CM(NRP1,MEDIM),
     &            RGV2RV2C(NRP1,MEDIM),RGV2RV2CM(NRP1,MEDIM),
     &            RGV3RV2C(NRP1,MEDIM),RGV3RV2CM(NRP1,MEDIM),
     &            RGV1RV3C(NRP1,MEDIM),RGV1RV3CM(NRP1,MEDIM),
     &            RGV2RV3C(NRP1,MEDIM),RGV2RV3CM(NRP1,MEDIM))
      ENDIF

      RGV1RV1C  = 0.
      RGV1RV1CM = 0.
      RGV2RV1C  = 0.
      RGV2RV1CM = 0.
      RGV3RV1C  = 0.
      RGV3RV1CM = 0.
      RGV1RV2C  = 0.
      RGV1RV2CM = 0.
      RGV2RV2C  = 0.
      RGV2RV2CM = 0.
      RGV3RV2C  = 0.
      RGV3RV2CM = 0.
      RGV1RV3C  = 0.
      RGV1RV3CM = 0.
      RGV2RV3C  = 0.
      RGV2RV3CM = 0.

C     rho*F*psi'*G12L/B^2*RW1
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RHO(I)*T(I)*DPSIDS(I)*G12L(I,J)/B_2(I,J)
     &                 *RW1(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*TM(I)*DPSIDSM(I)*G12LM(I,J)/B_2M(I,J)
     &                 *RW1M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV1C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV1C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV1CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV1CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV1C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV1C')
        call FFTOUTPT(RW0M,   RGV1RV1CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV1CM')
      ENDIF

C     rho*psi'*R^2*G22L/B^2*RW1
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RHO(I)*DPSIDS(I)*G33L(I,J)*G22L(I,J)/B_2(I,J)
     &              *RW1(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*DPSIDSM(I)*G33LM(I,J)*G22LM(I,J)/B_2M(I,J)
     &              *RW1M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV1C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV1C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV1CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV1CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV1C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV1C')
        call FFTOUTPT(RW0M,   RGV2RV1CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV1CM')
      ENDIF

C     rho*F*J*RW1
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RHO(I)*T(I)*RJA(I,J)*RW1(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*TM(I)*RJAM(I,J)*RW1M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV3RV1C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV1C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV3RV1CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV1CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV3RV1C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3RV1C')
        call FFTOUTPT(RW0M,   RGV3RV1CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV3RV1CM')
      ENDIF

C     rho*F*psi'*G12L*RCC*F/B^4
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RHO(I)*T(I)**2*DPSIDS(I)*G12L(I,J)
     &                    *RCC(I,J)/B_2(I,J)**2
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*TM(I)**2*DPSIDSM(I)*G12LM(I,J)
     &                      *RCCM(I,J)/B_2M(I,J)**2
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV2C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV2C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV2CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV2CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV2C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV2C')
        call FFTOUTPT(RW0M,   RGV1RV2CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV2CM')
      ENDIF

C     rho*psi'*R^2*G22L*RCC*F/B^4
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RHO(I)*DPSIDS(I)*G33L(I,J)*G22L(I,J)
     &                    *RCC(I,J)*T(I)/B_2(I,J)**2
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*DPSIDSM(I)*G33LM(I,J)*G22LM(I,J)
     &                      *RCCM(I,J)*TM(I)/B_2M(I,J)**2
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV2C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV2C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV2CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV2CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV2C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV2C')
        call FFTOUTPT(RW0M,   RGV2RV2CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV2CM')
      ENDIF

C     rho*F^2*J*RCC/B^2
      DO J=1,NCHI
        DO I=1,NRP1
           RW0(I,J)=RHO(I)*T(I)**2*RJA(I,J)*RCC(I,J)/B_2(I,J)
        ENDDO
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*TM(I)**2*RJAM(I,J)*RCCM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV3RV2C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV2C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV3RV2CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV2CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV3RV2C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3RV2C')
        call FFTOUTPT(RW0M,   RGV3RV2CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV3RV2CM')
      ENDIF

C     rho*F*psi'^2*G12L/B^2*RCC*/J
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=RHO(I)*T(I)*DPSIDS(I)**2*G12L(I,J)/B_2(I,J)
     &                    *RCC(I,J)/RJA(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J) 
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*TM(I)*DPSIDSM(I)**2*G12LM(I,J)/B_2M(I,J)
     &                    *RCCM(I,J)/RJAM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV3C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV3C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV3CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV3CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV3C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV3C')
        call FFTOUTPT(RW0M,   RGV1RV3CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV3CM')
      ENDIF

C     rho*psi'^2*R^2*G22L/B^2*RCC*/J
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=RHO(I)*DPSIDS(I)**2*G33L(I,J)*G22L(I,J)/B_2(I,J)
     &                    *RCC(I,J)/RJA(I,J)
        ENDDO
        RW0(1,J) = RW0(2,J) 
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*DPSIDSM(I)**2*G33LM(I,J)*G22LM(I,J)
     &                      *RCCM(I,J)/RJAM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV3C, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV3C  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV3CM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV3CM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      8,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV3C,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV3C')
        call FFTOUTPT(RW0M,   RGV2RV3CM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV3CM')
      ENDIF



      DEALLOCATE( RW0,RW0M,RW1,RW1M,B_2,B_2M,RCS,RCSM,RCC,RCCM )

      RETURN
      END

C====================================================================
C PREPARE EQUILIBRIUM FOURIER COEFFICIENTS FOR SUBROUTINE POLOIDAL  =
C ROTATION, SHOULD BE CALLED FROM WITHIN FTCOEFF                    =
C YQL,GLX & LLI------11/2017                                        =
C====================================================================
      SUBROUTINE PFTCOEFF

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      IMPLICIT NONE

      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER NPSTRT
      INTEGER I,J
      REAL*8  H1,H2
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW0,RW0M,RW1,RW1M,RW2,RW2M,
     &                                   RW3,RW3M,RW4,RW4M,
     &                                   RW1D,RW1DM,RW2D,RW2DM,
     &                                   B_2,B_2M,B2S,B2SM,B2C,B2CM   

      ALLOCATE( RW0(NRP1,NCHI),  RW0M(NRP1,NCHI),
     &          RW1(NRP1,NCHI),  RW1M(NRP1,NCHI),
     &          RW2(NRP1,NCHI),  RW2M(NRP1,NCHI),
     &          RW3(NRP1,NCHI),  RW3M(NRP1,NCHI),
     &          RW4(NRP1,NCHI),  RW4M(NRP1,NCHI),
     &          RW1D(NRP1,NCHI), RW1DM(NRP1,NCHI),
     &          RW2D(NRP1,NCHI), RW2DM(NRP1,NCHI),
     &          B_2(NRP1,NCHI),  B_2M(NRP1,NCHI),
     &          B2S(NRP1,NCHI),  B2SM(NRP1,NCHI),
     &          B2C(NRP1,NCHI),  B2CM(NRP1,NCHI) )

      INCLUDE 'setfft.inc'
      SUBNAM    = 'PFTCOEFF'

C     B^2
      DO J=1,NCHI
        DO I=2,NRP1
          B_2(I,J) =G22L(I,J)*DPSIDS(I)**2/RJA(I,J)**2 +
     &              T(I)**2/REQ(I,J)**2
        ENDDO 
        B_2(1,J)=T(1)**2/REQ(1,J)**2
        DO I=1,NR
          B_2M(I,J)=G22LM(I,J)*DPSIDSM(I)**2/RJAM(I,J)**2 +
     &              TM(I)**2/REQM(I,J)**2
        ENDDO 
      ENDDO 

C     DB^2/DCHI
      CALL DERCHI(B_2, B2C, NRP1,NRP1)
      CALL DERCHI(B_2M,B2CM,NR,  NRP1)

C     DB^2/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            B2S(I,J)=(H1/H2*B_2M(I,J)-H2/H1*B_2M(I-1,J))/(H1+H2) -
     &               (H1-H2)*B_2(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         B2S(I,J)=(B_2(I-1,J)+3*B_2(I,J)-4*B_2M(I-1,J))/H1
         B2S(1,J)=B2S(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            B2SM(I,J)=(B_2(I+1,J)-B_2(I,J))/H1
         ENDDO
      ENDDO

C GLX
      IF (.NOT.ALLOCATED(RGX2RX1P)) THEN
      ALLOCATE(   RGX2RX1P(NRP1,MEDIM),RGX2RX1PM(NRP1,MEDIM),
     &            RGX2RX2P(NRP1,MEDIM),RGX2RX2PM(NRP1,MEDIM)  )
      ENDIF

      RGX2RX1P  = 0.
      RGX2RX1PM = 0.
      RGX2RX2P  = 0.
      RGX2RX2PM = 0.

C GLX  PSI'/J 
      DO J=1,NCHI
        DO I=2,NRP1
           RW1(I,J) =DPSIDS(I)/RJA(I,J)
        ENDDO
        RW1(1,J)=0.0
        DO I=1,NR
           RW1M(I,J)=DPSIDSM(I)/RJAM(I,J)
        ENDDO
      ENDDO

C      DRW1/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RW1D(I,J)=(H1/H2*RW1M(I,J)-H2/H1*RW1M(I-1,J))/(H1+H2) -
     &                (H1-H2)*RW1(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RW1D(I,J)=(RW1(I-1,J)+3*RW1(I,J)-4*RW1M(I-1,J))/H1
         RW1D(1,J)=RW1D(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RW1DM(I,J)=(RW1(I+1,J)-RW1(I,J))/H1
         ENDDO
      ENDDO

C GLX  F/R^2 
      DO J=1,NCHI
        DO I=2,NRP1
           RW2(I,J) =T(I)/REQ(I,J)**2
        ENDDO
        RW2(1,J)=RW2(2,J)
        DO I=1,NR
           RW2M(I,J)=TM(I)/REQM(I,J)**2
        ENDDO
      ENDDO

C      DRW2/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RW2D(I,J)=(H1/H2*RW2M(I,J)-H2/H1*RW2M(I-1,J))/(H1+H2) -
     &                (H1-H2)*RW2(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RW2D(I,J)=(RW2(I-1,J)+3*RW2(I,J)-4*RW2M(I-1,J))/H1
         RW2D(1,J)=RW2D(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RW2DM(I,J)=(RW2(I+1,J)-RW2(I,J))/H1
         ENDDO
      ENDDO

C GLX  RHO^-1*U*[J*F/R^2*D(PSI'/J)/DS-PSI'*D(F/R^2)/DS]
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =RHOU(I)*(RJA(I,J)*T(I)/REQ(I,J)**2*RW1D(I,J) -
     &               DPSIDS(I)*RW2D(I,J))
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=RHOUM(I)*(RJAM(I,J)*TM(I)/REQM(I,J)**2*RW1DM(I,J) -
     &               DPSIDSM(I)*RW2DM(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGX2RX1P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX1P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGX2RX1PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX1PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGX2RX1P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGX2RX1P')
        call FFTOUTPT(RW0M,   RGX2RX1PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGX2RX1PM')
      ENDIF

C GLX  DJ/DCHI
      CALL DERCHI(RJA, RW1D, NRP1,NRP1)
      CALL DERCHI(RJAM,RW1DM,NR,  NRP1)

      DO J=1,NCHI
         RW1D(1,J)= 0.0
      ENDDO

C GLX  RHO^-1*U*PSI'/J*DJ/DCHI
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =RHOU(I)*DPSIDS(I)/RJA(I,J)*RW1D(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=RHOUM(I)*DPSIDSM(I)/RJAM(I,J)*RW1DM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGX2RX2P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX2P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGX2RX2PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGX2RX2PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGX2RX2P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGX2RX2P')
        call FFTOUTPT(RW0M,   RGX2RX2PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGX2RX2PM')
      ENDIF

C GLX
      IF (.NOT.ALLOCATED(RGV1RV1P)) THEN
      ALLOCATE(   RGV1RV1P(NRP1,MEDIM),RGV1RV1PM(NRP1,MEDIM),
     &            RGV1RV2P(NRP1,MEDIM),RGV1RV2PM(NRP1,MEDIM),
     &            RGV1RV3P(NRP1,MEDIM),RGV1RV3PM(NRP1,MEDIM),
     &            RGV1RV4P(NRP1,MEDIM),RGV1RV4PM(NRP1,MEDIM),
     &            RGV1RV5P(NRP1,MEDIM),RGV1RV5PM(NRP1,MEDIM),
     &            RGV1RV6P(NRP1,MEDIM),RGV1RV6PM(NRP1,MEDIM),
     &            RGV1RV7P(NRP1,MEDIM),RGV1RV7PM(NRP1,MEDIM)  )
      ENDIF

      RGV1RV1P  = 0.
      RGV1RV1PM = 0.
      RGV1RV2P  = 0.
      RGV1RV2PM = 0.
      RGV1RV3P  = 0.
      RGV1RV3PM = 0.
      RGV1RV4P  = 0.
      RGV1RV4PM = 0.
      RGV1RV5P  = 0.
      RGV1RV5PM = 0.
      RGV1RV6P  = 0.
      RGV1RV6PM = 0.
      RGV1RV7P  = 0.
      RGV1RV7PM = 0.

C GLX  U*PSI'*[G11L-(G12L*PSI'/B/J)^2]
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*DPSIDS(I)*(G11L(I,J)-(G12L(I,J)*
     &               DPSIDS(I)/RJA(I,J))**2/B_2(I,J))
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*DPSIDSM(I)*(G11LM(I,J)-(G12LM(I,J)*
     &               DPSIDSM(I)/RJAM(I,J))**2/B_2M(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV1P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV1P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV1PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV1PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV1P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV1P')
        call FFTOUTPT(RW0M,   RGV1RV1PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV1PM')
      ENDIF

C GLX  U*J*F/R^2*[G11L-(G12L*PSI'/B/J)^2]
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*RJA(I,J)*T(I)/REQ(I,J)**2*(G11L(I,J) -
     &               (G12L(I,J)*DPSIDS(I)/RJA(I,J))**2/B_2(I,J))
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*RJAM(I,J)*TM(I)/REQM(I,J)**2*(G11LM(I,J) -
     &               (G12LM(I,J)*DPSIDSM(I)/RJAM(I,J))**2/B_2M(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV2P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV2P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV2PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV2PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV2P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV2P')
        call FFTOUTPT(RW0M,   RGV1RV2PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV2PM')
      ENDIF

C GLX  U*PSI'*F*G12L/B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*DPSIDS(I)*T(I)*G12L(I,J)/B_2(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*DPSIDSM(I)*TM(I)*G12LM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV3P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV3P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV3PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV3PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV3P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV3P')
        call FFTOUTPT(RW0M,   RGV1RV3PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV3PM')
      ENDIF

C GLX  U*J*F^2*G12L/B^2/R^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*RJA(I,J)*T(I)**2*G12L(I,J)/
     &               B_2(I,J)/REQ(I,J)**2
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*RJAM(I,J)*TM(I)**2*G12LM(I,J)/
     &               B_2M(I,J)/REQM(I,J)**2
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV4P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV4P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV4PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV4PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV4P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV4P')
        call FFTOUTPT(RW0M,   RGV1RV4PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV4PM')
      ENDIF

C GLX  J*F/PSI'/R^2 
      DO J=1,NCHI
        DO I=2,NRP1
           RW1(I,J) =RJA(I,J)*T(I)/DPSIDS(I)/REQ(I,J)**2
        ENDDO
        RW1(1,J) = RW1(2,J) 
        DO I=1,NR
           RW1M(I,J)=RJAM(I,J)*TM(I)/DPSIDSM(I)/REQM(I,J)**2
        ENDDO
      ENDDO

C      DRW1/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RW1D(I,J)=(H1/H2*RW1M(I,J)-H2/H1*RW1M(I-1,J))/(H1+H2) -
     &                (H1-H2)*RW1(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RW1D(I,J)=(RW1(I-1,J)+3*RW1(I,J)-4*RW1M(I-1,J))/H1
         RW1D(1,J)=0.0

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RW1DM(I,J)=(RW1(I+1,J)-RW1(I,J))/H1
         ENDDO
      ENDDO

C GLX  G11L-(G12L*PSI'/B/J)^2 
      DO J=1,NCHI
        DO I=2,NRP1
           RW2(I,J) =G11L(I,J)-(G12L(I,J)*DPSIDS(I)/
     &               RJA(I,J))**2/B_2(I,J)
        ENDDO
        RW2(1,J) = RW2(2,J) 
        DO I=1,NR
           RW2M(I,J)=G11LM(I,J)-(G12LM(I,J)*DPSIDSM(I)/
     &               RJAM(I,J))**2/B_2M(I,J)
        ENDDO
      ENDDO

C      DRW2/DCHI
      CALL DERCHI(RW2, RW2D, NRP1,NRP1)
      CALL DERCHI(RW2M,RW2DM,NR,  NRP1)

C GLX  U*(G12L*PSI'^2*F/B^2/J*DRW1/DS-PSI'*DRW2/DCHI)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*(G12L(I,J)*DPSIDS(I)**2*T(I)/RJA(I,J)/
     &               B_2(I,J)*RW1D(I,J)-DPSIDS(I)*RW2D(I,J))
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*(G12LM(I,J)*DPSIDSM(I)**2*TM(I)/RJAM(I,J)/
     &               B_2M(I,J)*RW1DM(I,J)-DPSIDSM(I)*RW2DM(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV5P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV5P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV5PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV5PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV5P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV5P')
        call FFTOUTPT(RW0M,   RGV1RV5PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV5PM')
      ENDIF

C GLX  ln(J*F/PSI'/R^2) 
      DO J=1,NCHI
        DO I=2,NRP1
           RW1(I,J) =log(RJA(I,J)*T(I)/DPSIDS(I)/REQ(I,J)**2)
        ENDDO
        RW1(1,J) = RW1(2,J) 
        DO I=1,NR
           RW1M(I,J)=log(RJAM(I,J)*TM(I)/DPSIDSM(I)/REQM(I,J)**2)
        ENDDO
      ENDDO

C      DRW1/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RW1D(I,J)=(H1/H2*RW1M(I,J)-H2/H1*RW1M(I-1,J))/(H1+H2) -
     &                (H1-H2)*RW1(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RW1D(I,J)=(RW1(I-1,J)+3*RW1(I,J)-4*RW1M(I-1,J))/H1
         RW1D(1,J)=0.0

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RW1DM(I,J)=(RW1(I+1,J)-RW1(I,J))/H1
         ENDDO
      ENDDO

C GLX  F*G12L/B^2 
      DO J=1,NCHI
        DO I=2,NRP1
           RW2(I,J) =T(I)*G12L(I,J)/B_2(I,J)
        ENDDO
        RW2(1,J) = RW2(2,J) 
        DO I=1,NR
           RW2M(I,J)=TM(I)*G12LM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

C      DRW2/DCHI
      CALL DERCHI(RW2, RW2D, NRP1,NRP1)
      CALL DERCHI(RW2M,RW2DM,NR,  NRP1)

C GLX  U*(-J^2*F'-J^2*F/B^2*P'+G22L*PSI'*F/B^2*DRW1/DS-PSI'*DRW2/DCHI)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*(-RJA(I,J)**2*TP(I)-RJA(I,J)**2*
     &               T(I)/B_2(I,J)*PPEQ(I)+G22L(I,J)*DPSIDS(I)*
     &               T(I)/B_2(I,J)*RW1D(I,J)-DPSIDS(I)*RW2D(I,J))
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*(-RJAM(I,J)**2*TPM(I)-RJAM(I,J)**2*
     &               TM(I)/B_2M(I,J)*PPEQM(I)+G22LM(I,J)*DPSIDSM(I)*
     &               TM(I)/B_2M(I,J)*RW1DM(I,J)-DPSIDSM(I)*RW2DM(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV6P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV6P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV6PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV6PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV6P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV6P')
        call FFTOUTPT(RW0M,   RGV1RV6PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV6PM')
      ENDIF

C GLX  U*(J*DB^2/DS+2*J*DP/DS-G12L*PSI'^2/B^2/J*DB^2/DCHI)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*(RJA(I,J)*B2S(I,J)+2.0*RJA(I,J)*
     &               PPEQ(I)*DPSIDS(I)-G12L(I,J)*DPSIDS(I)**2/
     &               B_2(I,J)/RJA(I,J)*B2C(I,J))
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*(RJAM(I,J)*B2SM(I,J)+2.0*RJAM(I,J)*
     &               PPEQM(I)*DPSIDSM(I)-G12LM(I,J)*DPSIDSM(I)**2/
     &               B_2M(I,J)/RJAM(I,J)*B2CM(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV1RV7P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV7P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV1RV7PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV1RV7PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV1RV7P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV1RV7P')
        call FFTOUTPT(RW0M,   RGV1RV7PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV1RV7PM')
      ENDIF

C GLX
      IF (.NOT.ALLOCATED(RGV2RV1P)) THEN
      ALLOCATE(   RGV2RV1P(NRP1,MEDIM),RGV2RV1PM(NRP1,MEDIM),
     &            RGV2RV2P(NRP1,MEDIM),RGV2RV2PM(NRP1,MEDIM),
     &            RGV2RV3P(NRP1,MEDIM),RGV2RV3PM(NRP1,MEDIM),
     &            RGV2RV4P(NRP1,MEDIM),RGV2RV4PM(NRP1,MEDIM),
     &            RGV2RV5P(NRP1,MEDIM),RGV2RV5PM(NRP1,MEDIM)  )
      ENDIF

      RGV2RV1P  = 0.
      RGV2RV1PM = 0.
      RGV2RV2P  = 0.
      RGV2RV2PM = 0.
      RGV2RV3P  = 0.
      RGV2RV3PM = 0.
      RGV2RV4P  = 0.
      RGV2RV4PM = 0.
      RGV2RV5P  = 0.
      RGV2RV5PM = 0.

C GLX  U*PSI'*G22L*R^2/B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*DPSIDS(I)*G22L(I,J)*
     &               G33L(I,J)/B_2(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*DPSIDSM(I)*G22LM(I,J)*
     &               G33LM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV1P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV1P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV1PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV1PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV1P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV1P')
        call FFTOUTPT(RW0M,   RGV2RV1PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV1PM')
      ENDIF

C GLX  U*J*F*G22L/B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*RJA(I,J)*T(I)*G22L(I,J)/B_2(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*RJAM(I,J)*TM(I)*G22LM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV2P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV2P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV2PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV2PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV2P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV2P')
        call FFTOUTPT(RW0M,   RGV2RV2PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV2PM')
      ENDIF

C GLX  G12L*PSI'*F/B^2/J 
      DO J=1,NCHI
        DO I=2,NRP1
           RW1(I,J) =G12L(I,J)*DPSIDS(I)*T(I)/B_2(I,J)/RJA(I,J)
        ENDDO
        RW1(1,J) = RW1(2,J) 
        DO I=1,NR
           RW1M(I,J)=G12LM(I,J)*DPSIDSM(I)*TM(I)/B_2M(I,J)/RJAM(I,J)
        ENDDO
      ENDDO

C      DRW1/DCHI
      CALL DERCHI(RW1, RW1D, NRP1,NRP1)
      CALL DERCHI(RW1M,RW1DM,NR,  NRP1)

C GLX  U*(-J^2*F'-J^2*F/B^2*P'+J*DRW1/DCHI)
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*(-RJA(I,J)**2*TP(I)-RJA(I,J)**2*
     &               T(I)/B_2(I,J)*PPEQ(I)+RJA(I,J)*RW1D(I,J))
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*(-RJAM(I,J)**2*TPM(I)-RJAM(I,J)**2*
     &               TM(I)/B_2M(I,J)*PPEQM(I)+RJAM(I,J)*RW1DM(I,J))
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV3P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV3P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV3PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV3PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV3P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV3P')
        call FFTOUTPT(RW0M,   RGV2RV3PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV3PM')
      ENDIF

C GLX  G22L*PSI'*R^2/B^2/J 
      DO J=1,NCHI
        DO I=2,NRP1
           RW1(I,J) =G22L(I,J)*DPSIDS(I)*G33L(I,J)/
     &               B_2(I,J)/RJA(I,J)
        ENDDO
        RW1(1,J) = RW1(2,J) 
        DO I=1,NR
           RW1M(I,J)=G22LM(I,J)*DPSIDSM(I)*G33LM(I,J)/
     &               B_2M(I,J)/RJAM(I,J)
        ENDDO
      ENDDO

C      DRW1/DCHI
      CALL DERCHI(RW1, RW1D, NRP1,NRP1)
      CALL DERCHI(RW1M,RW1DM,NR,  NRP1)

C GLX  U*J*DRW1/DCHI
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*RJA(I,J)*RW1D(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*RJAM(I,J)*RW1DM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV4P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV4P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV4PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV4PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV4P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV4P')
        call FFTOUTPT(RW0M,   RGV2RV4PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV4PM')
      ENDIF

C GLX  U*J*F/B^2*DB^2/DCHI
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*RJA(I,J)*T(I)/B_2(I,J)*B2C(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*RJAM(I,J)*TM(I)/B_2M(I,J)*B2CM(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV2RV5P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV5P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV2RV5PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV2RV5PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV2RV5P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV2RV5P')
        call FFTOUTPT(RW0M,   RGV2RV5PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV2RV5PM')
      ENDIF

C GLX
      IF (.NOT.ALLOCATED(RGV3RV1P)) THEN
      ALLOCATE(   RGV3RV1P(NRP1,MEDIM),RGV3RV1PM(NRP1,MEDIM),
     &            RGV3RV2P(NRP1,MEDIM),RGV3RV2PM(NRP1,MEDIM),
     &            RGV3RV3P(NRP1,MEDIM),RGV3RV3PM(NRP1,MEDIM)   )
      ENDIF

      RGV3RV1P  = 0.
      RGV3RV1PM = 0.
      RGV3RV2P  = 0.
      RGV3RV2PM = 0.
      RGV3RV3P  = 0.
      RGV3RV3PM = 0.

C GLX  U*PSI'*B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*DPSIDS(I)*B_2(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*DPSIDSM(I)*B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV3RV1P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV1P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV3RV1PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV1PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV3RV1P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3RV1P')
        call FFTOUTPT(RW0M,   RGV3RV1PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV3RV1PM')
      ENDIF

C GLX  U*J*F/R^2*B^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =ROTP(I)*RJA(I,J)*T(I)/REQ(I,J)**2*B_2(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=ROTPM(I)*RJAM(I,J)*TM(I)/REQM(I,J)**2*B_2M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV3RV2P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV2P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV3RV2PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV2PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV3RV2P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3RV2P')
        call FFTOUTPT(RW0M,   RGV3RV2PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV3RV2PM')
      ENDIF

C GLX
C RHO*J*B^2*D(RHO^-1*U)/DS+U*J*DB^2/DS/2.0-U*PSI'^2*G12L/B^2/J*DB^2/DCHI/2.0 
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =RHO(I)*RJA(I,J)*B_2(I,J)*DRHOU(I) +
     &               ROTP(I)*RJA(I,J)*B2S(I,J)*0.5 -
     &               ROTP(I)*DPSIDS(I)**2*G12L(I,J)/
     &               B_2(I,J)/RJA(I,J)*B2C(I,J)*0.5
        ENDDO
        RW0(1,J)=RW0(2,J)
        DO I=1,NR
           RW0M(I,J)=RHOM(I)*RJAM(I,J)*B_2M(I,J)*DRHOUM(I) +
     &               ROTPM(I)*RJAM(I,J)*B2SM(I,J)*0.5 -
     &               ROTPM(I)*DPSIDSM(I)**2*G12LM(I,J)/
     &               B_2M(I,J)/RJAM(I,J)*B2CM(I,J)*0.5
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGV3RV3P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV3P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGV3RV3PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGV3RV3PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGV3RV3P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGV3RV3P')
        call FFTOUTPT(RW0M,   RGV3RV3PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGV3RV3PM')
      ENDIF

C GLX
      IF (.NOT.ALLOCATED(RGQ1RQ1P)) THEN
      ALLOCATE(   RGQ1RQ1P(NRP1,MEDIM),RGQ1RQ1PM(NRP1,MEDIM),
     &            RGQ1RQ2P(NRP1,MEDIM),RGQ1RQ2PM(NRP1,MEDIM),
     &            RGQ2RQ1P(NRP1,MEDIM),RGQ2RQ1PM(NRP1,MEDIM),
     &            RGQ3RQ1P(NRP1,MEDIM),RGQ3RQ1PM(NRP1,MEDIM)  )
      ENDIF

      RGQ1RQ1P  = 0.
      RGQ1RQ1PM = 0.
      RGQ1RQ2P  = 0.
      RGQ1RQ2PM = 0.
      RGQ2RQ1P  = 0.
      RGQ2RQ1PM = 0.
      RGQ3RQ1P  = 0.
      RGQ3RQ1PM = 0.

C GLX  RHO^-1*U*PSI'/J
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =RHOU(I)*DPSIDS(I)/RJA(I,J)
           RW3(I,J) =RW0(I,J)
        ENDDO
        RW0(1,J)=0.0
        RW3(1,J)=0.0
        DO I=1,NR
           RW0M(I,J)=RHOUM(I)*DPSIDSM(I)/RJAM(I,J)
           RW3M(I,J)=RW0M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGQ1RQ1P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ1RQ1P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGQ1RQ1PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ1RQ1PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGQ1RQ1P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGQ1RQ1P')
        call FFTOUTPT(RW0M,   RGQ1RQ1PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGQ1RQ1PM')
      ENDIF

C GLX  RHO^-1*U*F/R^2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J) =RHOU(I)*T(I)/REQ(I,J)**2
           RW4(I,J) =RW0(I,J)
        ENDDO
        RW0(1,J)=RW0(2,J)
        RW4(1,J)=RW4(2,J)
        DO I=1,NR
           RW0M(I,J)=RHOUM(I)*TM(I)/REQM(I,J)**2
           RW4M(I,J)=RW0M(I,J)
        ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGQ1RQ2P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ1RQ2P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGQ1RQ2PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ1RQ2PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGQ1RQ2P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGQ1RQ2P')
        call FFTOUTPT(RW0M,   RGQ1RQ2PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGQ1RQ2PM')
      ENDIF

C GLX  D(RHO^-1*U*PSI'/J)/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RW0(I,J)=(H1/H2*RW3M(I,J)-H2/H1*RW3M(I-1,J))/(H1+H2) -
     &               (H1-H2)*RW3(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RW0(I,J)=(RW3(I-1,J)+3*RW3(I,J)-4*RW3M(I-1,J))/H1
         RW0(1,J)=RW0(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RW0M(I,J)=(RW3(I+1,J)-RW3(I,J))/H1
         ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGQ2RQ1P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ2RQ1P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGQ2RQ1PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ2RQ1PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGQ2RQ1P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGQ2RQ1P')
        call FFTOUTPT(RW0M,   RGQ2RQ1PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGQ2RQ1PM')
      ENDIF

C GLX  D(RHO^-1*U*F/R^2)/DS
      DO J=1,NCHI
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            RW0(I,J)=(H1/H2*RW4M(I,J)-H2/H1*RW4M(I-1,J))/(H1+H2) -
     &               (H1-H2)*RW4(I,J)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         RW0(I,J)=(RW4(I-1,J)+3*RW4(I,J)-4*RW4M(I-1,J))/H1
         RW0(1,J)=RW0(2,J)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            RW0M(I,J)=(RW4(I+1,J)-RW4(I,J))/H1
         ENDDO
      ENDDO

      NPSTRT    =  1
      call FFTDRIVER(RW0,  RGQ3RQ1P, FORWD, NRP1,  NRP1,   NPSTRT
     &                    ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ2RQ1P  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif

      call FFTDRIVER(RW0M, RGQ3RQ1PM, FORWD, NRP1,  NR,     NPSTRT
     &                    ,MEDIM,    NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: RGQ3RQ1PM  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      6,   MESSAGE
     &          ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RW0,    RGQ3RQ1P,  NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'RGQ3RQ1P')
        call FFTOUTPT(RW0M,   RGQ3RQ1PM, NRP1,    NR,      NPSTRT
     &                       ,MEDIM,    NCHI,    KUFFTP, 'RGQ3RQ1PM')
      ENDIF

      DEALLOCATE(  RW0,RW0M,RW1,RW1M,RW2,RW2M,
     &             RW3,RW3M,RW4,RW4M,
     &             RW1D,RW1DM,RW2D,RW2DM,
     &             B_2,B_2M,B2S,B2SM,B2C,B2CM  )
 
      RETURN
      END
C     GLX------> THE END OF PFTCOEFF

*DECK COEFFI
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C------- MATRIX COEFFICIENTS -------- A.B. 14.02.91 --------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     VERSION FOR GENERALISED FEM WITH INTEGRATION SCHEME THAT IS AN
C     ARBITRARY MIXTURE OF TRAPETZOIDAL AND MID-POINT RULES.  PTRAP
C     IS THE FRACTION OF TRAPETZOIDAL RULE AND THE THREE POINTS IN A
C     CELL ARE WEIGHTED BY PTRAP/2, 1-PTRAP, AND PTRAP/2.
C
C        PTRAP = 0          HYBRID ELEMENTS
C        PTRAP = 1/3        FINITE ELEMENTS
C        PTRAP = 1          FINITE DIFFERENCES (APPROXIMATELY)
C
C***********************************************************************
C
      SUBROUTINE COEFFI(MROW,MSA ,MSB,CMROW,CMA,CMB,CNA,SHIFT,
     &     SHIFTC,SHIFTM,SHIFTVC,SHIFTVM,
     &   ASUBM, BSUBM, CSUBM, DSUBM, ESUBM, FSUBM, GSUBM, HSUBM)
C
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE KINETICM
      USE GIJLM
      USE CONVOLCOFM
      USE REORBITM
      IMPLICIT NONE
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
C
      INTEGER       MROW,MSA,MSB,NSA,NSB
      PARAMETER     (NSA=2,NSB=1)
      INTEGER       LXCOL,LYCOL,LXROW,LYROW,I,K,KP,LXCOL2,LYCOL2
      REAL*8        ZKPAVT,TMP,OCIRC,OTRAP,DCIRC,DTRAP,EPS
      COMPLEX*16    CMROW,CMA,CMB,CNA,SHIFT,CTMP,
     &              CTMP1,CTMP2,CTMP3,CTMP4,CTMP5,CTMP6,
     &              CTMP7,CTMP8,CTMP9,CTMP10,CTMP11,CTMP12
      COMPLEX*16,DIMENSION(NRP1)::SHIFTC,SHIFTM,SHIFTVC,SHIFTVM
      integer       iexe,iexv2,iexj1,iexb3
      parameter     (iexv2 = -1, iexj1 = 1, iexb3 = 1)
      real*8        zem,zep,zv2m,zv2p,z1m,z1p,z3m,z3p
 
C     NOTE: EQUATIONS FOR V2,X2,V3,X3 ARE MULTIPLIED BY LOCAL FACTORS
C     ZEM(ZEP); SOLUTION VARIABLES V2,X2,J1,B3 ARE MULTIPLIED BY
C     LOCAL FACTORS ZV2M(ZV2P),ZV2M(ZV2P),Z1M(Z1P),Z3M(Z3P), RESPECTIVELY.  

      INCLUDE 'integc.inc'
C
      LXROW = (MROW-1)*NXCOMP
      LXCOL = (MSA -1)*NXCOMP
      LYROW = (MROW-1)*NYCOMP
      LYCOL = (MSA -1)*NYCOMP
C
C-----------------------------------------------------------------------
C.. FIRST EQUATION: COVARIANT-S-COMP OF EQ. OF MOTION
C..                              (KXV1 = 1, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
C
      DO 10 I=2,NRP1
      K = I
C
      INCLUDE 'tent.inc'
C
      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      ZV2P = (CS(I)/CSM(I  ))**IEXV2
      Z3M  = (CSM(I-1)/CS(I))**IEXB3
      Z3P  = (CSM(I)  /CS(I))**IEXB3

C
c     ------------------------------------------------------ >
c     Inertial Terms:
      
         BsubM(kxV1+lxrow, kxV1+lxcol,k)=
     $     BsubM(kxV1+lxrow, kxV1+lxcol,k) + DINERT*
     $     FF(SHIFTVC(I)*RGV1G11(i,msb),SHIFTVM(I-1)*RGV1G11M(i-1,msb),
     $     SHIFTVM(I)*RGV1G11M(i,msb))
         AsubM(kxV1+lxrow, kxV1+lxcol,k)=
     $     AsubM(kxV1+lxrow, kxV1+lxcol,k) + DINERT*
     $     SHIFTVM(I-1)*FFM(RGV1G11M(i-1,msb))
         CsubM(kxV1+lxrow, kxV1+lxcol,k)=
     $     CsubM(kxV1+lxrow, kxV1+lxcol,k) + DINERT*
     $     SHIFTVM(I)*FFP(RGV1G11M(i,msb))


         HsubM(kxV1+lxrow, kyV2+lycol,k)=
     $     HsubM(kxV1+lxrow, kyV2+lycol,k) + DINERT*
     $     FGM(SHIFTVC(I)*RGV1G12(i,msb)*zv2m,
     $     SHIFTVM(I-1)*RGV1G12M(i-1,msb))
         EsubM(kxV1+lxrow, kyV2+lycol,k)=
     $     EsubM(kxV1+lxrow, kyV2+lycol,k) + DINERT*
     $     FGP(SHIFTVC(I)*RGV1G12(i,msb)*zv2p,
     $     SHIFTVM(I)*RGV1G12M(i,msb))
c
C     LLI         
      IF (NPROFRC.GT.0) THEN
      BSUBM(KXV1+LXROW,KXV1+LXCOL,K) = BSUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &     -CNA*FF(RGV1G11C(I,MSB),RGV1G11CM(I-1,MSB),RGV1G11CM(I,MSB))
      ASUBM(KXV1+LXROW,KXV1+LXCOL,K) = ASUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &     -CNA*FFM(RGV1G11CM(I-1,MSB))
      CSUBM(KXV1+LXROW,KXV1+LXCOL,K) = CSUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &     -CNA*FFP(RGV1G11CM(I,MSB))


      HSUBM(KXV1+LXROW,KYV2+LYCOL,K) = HSUBM(KXV1+LXROW,KYV2+LYCOL,K)
     &     -CNA*FGM(ZV2M*RGV1G12C(I,MSB),RGV1G12CM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYV2+LYCOL,K) = ESUBM(KXV1+LXROW,KYV2+LYCOL,K)
     &     -CNA*FGP(ZV2P*RGV1G12C(I,MSB),RGV1G12CM(I,MSB))

      ENDIF
c     ------------------------------------------------------ <
c
      HsubM(kxV1+lxrow, kyB2+lycol,k)=
     $     -FGM(J3b2(i,msb),J3b2M(i-1,msb))
      EsubM(kxV1+lxrow, kyB2+lycol,k)=
     $     -FGP(J3b2(i,msb),J3b2M(i,msb))

      IF (ABS(RM(MSA,2)).GT.0.1) THEN
      BsubM(kxV1+lxrow, kxB1+lxcol,k)=
     $     FF(PPeq(i)*G12B2B2(i,msb), PPeqM(i-1)*G12B2B2M(i-1,msb), 
     $     PPeqM(i)*G12B2B2M(i,msb))
      AsubM(kxV1+lxrow, kxB1+lxcol,k)= FFM(PPeqM(i-1)*G12B2B2M(i-1,msb))
      CsubM(kxV1+lxrow, kxB1+lxcol,k)= FFP(PPeqM(i)*G12B2B2M(i,msb))
      ELSE
      BsubM(kxV1+lxrow, kxB1+lxcol,k)=
     $     FF(PPeq(i)*G12B2B2(i,msb)*T(i), 
     $     PPeqM(i-1)*G12B2B2M(i-1,msb)*TM(i-1), 
     $     PPeqM(i)*G12B2B2M(i,msb)*TM(i))
      AsubM(kxV1+lxrow, kxB1+lxcol,k)= 
     $     FFM(PPeqM(i-1)*G12B2B2M(i-1,msb)*TM(i-1))
      CsubM(kxV1+lxrow, kxB1+lxcol,k)= 
     $     FFP(PPeqM(i)*G12B2B2M(i,msb)*TM(i))
      ENDIF

      BsubM(kxV1+lxrow, kxJ2U+lxcol,k)=
     $   FF(B3j2(i,msb)*T(I),B3j2M(i-1,msb)*TM(I-1),B3j2M(i,msb)*TM(I))
      AsubM(kxV1+lxrow, kxJ2U+lxcol,k)= FFM(B3j2M(i-1,msb)*TM(I-1))
      CsubM(kxV1+lxrow, kxJ2U+lxcol,k)= FFP(B3j2M(i,msb)*TM(I))

      HsubM(kxV1+lxrow, kyPr+lycol,k)= znorm * JACOBI(i,msb)
     $     + cma *FGM(G12B2B2(i,msb), G12B2B2M(i-1,msb))
     $     + cna *FGM(G12B2B3(i,msb), G12B2B3M(i-1,msb))
      EsubM(kxV1+lxrow, kyPr+lycol,k)=-znorm * JACOBI(i,msb)
     $     + cma *FGP(G12B2B2(i,msb), G12B2B2M(i,msb)) 
     $     + cna *FGP(G12B2B3(i,msb), G12B2B3M(i,msb))

      IF (KYPE.GT.0.AND.INCKIN.EQ.0) THEN
      HsubM(kxV1+lxrow, kyPe+lycol,k)= znorm * JACOBI(i,msb)
     $     + cma *FGM(G12B2B2(i,msb), G12B2B2M(i-1,msb))
     $     + cna *FGM(G12B2B3(i,msb), G12B2B3M(i-1,msb))
      EsubM(kxV1+lxrow, kyPe+lycol,k)=-znorm * JACOBI(i,msb)
     $     + cma *FGP(G12B2B2(i,msb), G12B2B2M(i,msb)) 
     $     + cna *FGP(G12B2B3(i,msb), G12B2B3M(i,msb))
      ENDIF

 10   CONTINUE
C
C-----------------------------------------------------------------------
C      THE INERITAL TERMS DUE TO PLASMA ROTATION
C-----------------------------------------------------------------------
C
      DO 12 I=2,NRP1
      K = I
      INCLUDE 'tent.inc'
 
      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      ZV2P = (CS(I)/CSM(I  ))**IEXV2
      Z3M  = (CSM(I-1)/CS(I))**IEXB3
      Z3P  = (CSM(I)  /CS(I))**IEXB3
 
      IF (INERT2) THEN
      CTMP1 = IDRXX(I,MSB)*DROT(I)
      CTMP2 = IDRXXM(I-1,MSB)*DROTM(I-1)
      CTMP3 = IDRXXM(I,MSB)*DROTM(I)
      IF (IDIAMV.GE.3) THEN
         CTMP1 = IDRXX(I,MSB)*(DROT(I)+DOMEGASI(I)*FDIAMV)
         CTMP2 = IDRXXM(I-1,MSB)*(DROTM(I-1)+DOMEGASIM(I-1)*FDIAMV)
         CTMP3 = IDRXXM(I,MSB)*(DROTM(I)+DOMEGASIM(I)*FDIAMV)
      ENDIF
      BSUBM(KXV1+LXROW,KXV1+LXCOL,I) = BSUBM(KXV1+LXROW,KXV1+LXCOL,I)
     &                +FF(CTMP1,CTMP2,CTMP3)
      ASUBM(KXV1+LXROW,KXV1+LXCOL,I) = ASUBM(KXV1+LXROW,KXV1+LXCOL,I)
     &                +FFM(CTMP2)
      CSUBM(KXV1+LXROW,KXV1+LXCOL,I) = CSUBM(KXV1+LXROW,KXV1+LXCOL,I)
     &                +FFP(CTMP3)

C     LLI         
      IF (NPROFRC.GT.0) THEN
      BSUBM(KXV1+LXROW,KXV1+LXCOL,K) = BSUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &     +FF(RGV1RV1C(I,MSB),RGV1RV1CM(I-1,MSB),RGV1RV1CM(I,MSB))
      ASUBM(KXV1+LXROW,KXV1+LXCOL,K) = ASUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &     +FFM(RGV1RV1CM(I-1,MSB))
      CSUBM(KXV1+LXROW,KXV1+LXCOL,K) = CSUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &     +FFP(RGV1RV1CM(I,MSB))

      HSUBM(KXV1+LXROW,KYV2+LYCOL,K) = HSUBM(KXV1+LXROW,KYV2+LYCOL,K)
     &     +FGM(ZV2M*RGV1RV2C(I,MSB),RGV1RV2CM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYV2+LYCOL,K) = ESUBM(KXV1+LXROW,KYV2+LYCOL,K)
     &     +FGP(ZV2P*RGV1RV2C(I,MSB),RGV1RV2CM(I,MSB))

      IF (KYV3.GT.0) THEN
      HSUBM(KXV1+LXROW,KYV3+LYCOL,I)=HSUBM(KXV1+LXROW,KYV3+LYCOL,I)
     &     +FGM(RGV1RV3C(I,MSB),RGV1RV3CM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYV3+LYCOL,I)=ESUBM(KXV1+LXROW,KYV3+LYCOL,I)
     &     +FGP(RGV1RV3C(I,MSB),RGV1RV3CM(I,MSB))
      ENDIF
      ENDIF
C     GLX------<         
      IF (NPROFRP.GT.0) THEN
      BSUBM(KXV1+LXROW,KXV1+LXCOL,K) = BSUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &   - CMA*FF(RGV1RV1P(I,MSB),RGV1RV1PM(I-1,MSB),RGV1RV1PM(I,MSB))
     &   - CNA*FF(RGV1RV2P(I,MSB),RGV1RV2PM(I-1,MSB),RGV1RV2PM(I,MSB))
     &   +     FF(RGV1RV5P(I,MSB),RGV1RV5PM(I-1,MSB),RGV1RV5PM(I,MSB))
      ASUBM(KXV1+LXROW,KXV1+LXCOL,K) = ASUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &   - CMA*FFM(RGV1RV1PM(I-1,MSB))
     &   - CNA*FFM(RGV1RV2PM(I-1,MSB))
     &   +     FFM(RGV1RV5PM(I-1,MSB))
      CSUBM(KXV1+LXROW,KXV1+LXCOL,K) = CSUBM(KXV1+LXROW,KXV1+LXCOL,K)
     &   - CMA*FFP(RGV1RV1PM(I,MSB))
     &   - CNA*FFP(RGV1RV2PM(I,MSB))
     &   +     FFP(RGV1RV5PM(I,MSB))

      HSUBM(KXV1+LXROW,KYV2+LYCOL,K) = HSUBM(KXV1+LXROW,KYV2+LYCOL,K)
     &   - CMA*FGM(ZV2M*RGV1RV3P(I,MSB),RGV1RV3PM(I-1,MSB))
     &   - CNA*FGM(ZV2M*RGV1RV4P(I,MSB),RGV1RV4PM(I-1,MSB))
     &   +     FGM(ZV2M*RGV1RV6P(I,MSB),RGV1RV6PM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYV2+LYCOL,K) = ESUBM(KXV1+LXROW,KYV2+LYCOL,K)
     &   - CMA*FGP(ZV2P*RGV1RV3P(I,MSB),RGV1RV3PM(I,MSB))
     &   - CNA*FGP(ZV2P*RGV1RV4P(I,MSB),RGV1RV4PM(I,MSB))
     &   +     FGP(ZV2P*RGV1RV6P(I,MSB),RGV1RV6PM(I,MSB))

      IF (KYV3.GT.0) THEN
      HSUBM(KXV1+LXROW,KYV3+LYCOL,I)=HSUBM(KXV1+LXROW,KYV3+LYCOL,I)
     &   -     FGM(RGV1RV7P(I,MSB),RGV1RV7PM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYV3+LYCOL,I)=ESUBM(KXV1+LXROW,KYV3+LYCOL,I)
     &   -     FGP(RGV1RV7P(I,MSB),RGV1RV7PM(I,MSB))
      ENDIF
      ENDIF
C     --------->
      ENDIF 

      IF (INERT1) THEN
      CTMP1 = IRXY(I,MSB)*ROT(I)
      CTMP2 = IRXYM(I-1,MSB)*ROTM(I-1)
      CTMP3 = IRXYM(I,MSB)*ROTM(I)
      IF (IDIAMV.GE.3) THEN
         CTMP1 = IRXY(I,MSB)*(ROT(I)+OMEGASI(I)*FDIAMV)
         CTMP2 = IRXYM(I-1,MSB)*(ROTM(I-1)+OMEGASIM(I-1)*FDIAMV)
         CTMP3 = IRXYM(I,MSB)*(ROTM(I)+OMEGASIM(I)*FDIAMV)
      ELSEIF (IDIAMV.GE.2) THEN
         CTMP1 = IRXY(I,MSB)*(ROT(I)+OMEGASI(I)*.5*FDIAMV)
         CTMP2 = IRXYM(I-1,MSB)*(ROTM(I-1)+OMEGASIM(I-1)*.5*FDIAMV)
         CTMP3 = IRXYM(I,MSB)*(ROTM(I)+OMEGASIM(I)*.5*FDIAMV)
      ENDIF
      HSUBM(KXV1+LXROW,KYV2+LYCOL,I) = HSUBM(KXV1+LXROW,KYV2+LYCOL,I)
     &                   -FGM(CTMP1*ZV2M,CTMP2)
      ESUBM(KXV1+LXROW,KYV2+LYCOL,I) = ESUBM(KXV1+LXROW,KYV2+LYCOL,I)
     &                   -FGP(CTMP1*ZV2P,CTMP3)

C    LLI
      IF (NPROFRC.GT.0) THEN
      HSUBM(KXV1+LXROW,KYV2+LYCOL,I) = HSUBM(KXV1+LXROW,KYV2+LYCOL,I)
     &                   -FGM(RGIRXYC(I,MSB)*ZV2M,RGIRXYCM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYV2+LYCOL,I) = ESUBM(KXV1+LXROW,KYV2+LYCOL,I)
     &                   -FGP(RGIRXYC(I,MSB)*ZV2P,RGIRXYCM(I,MSB))
      ENDIF

      IF (KYV3.GT.0) THEN
      CTMP1 = IRXZ(I,MSB)*ROT(I)
      CTMP2 = IRXZM(I-1,MSB)*ROTM(I-1)
      CTMP3 = IRXZM(I,MSB)*ROTM(I)
      HSUBM(KXV1+LXROW,KYV3+LYCOL,I)=HSUBM(KXV1+LXROW,KYV3+LYCOL,I) +
     &                               FGM(CTMP1,CTMP2)
      ESUBM(KXV1+LXROW,KYV3+LYCOL,I)=ESUBM(KXV1+LXROW,KYV3+LYCOL,I) + 
     &                               FGP(CTMP1,CTMP3)
C     LLI
      IF (NPROFRC.GT.0) THEN
      HSUBM(KXV1+LXROW,KYV3+LYCOL,I)=HSUBM(KXV1+LXROW,KYV3+LYCOL,I)
     &                +FGM(RGIRXZC(I,MSB),RGIRXZCM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYV3+LYCOL,I)=ESUBM(KXV1+LXROW,KYV3+LYCOL,I)
     &                +FGP(RGIRXZC(I,MSB),RGIRXZCM(I,MSB))
      ENDIF
      ENDIF

      ENDIF

      IF (KYRHOP.GT.0) THEN
      HSUBM(KXV1+LXROW,KYRHOP+LYCOL,I)= HSUBM(KXV1+LXROW,KYRHOP+LYCOL,I)
     &                   +FGM(DXDZDSDZ(I,MSB),DXDZDSDZM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYRHOP+LYCOL,I)= ESUBM(KXV1+LXROW,KYRHOP+LYCOL,I)
     &                   +FGP(DXDZDSDZ(I,MSB),DXDZDSDZM(I,MSB))

C     LLI
      IF (NPROFRC.GT.0) THEN
      HSUBM(KXV1+LXROW,KYRHOP+LYCOL,I)= HSUBM(KXV1+LXROW,KYRHOP+LYCOL,I)
     &                   +FGM(RGDXDZDSDZC(I,MSB),RGDXDZDSDZCM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYRHOP+LYCOL,I)= ESUBM(KXV1+LXROW,KYRHOP+LYCOL,I)
     &                   +FGP(RGDXDZDSDZC(I,MSB),RGDXDZDSDZCM(I,MSB))
      ENDIF
      ENDIF

 12   CONTINUE
 
 15   CONTINUE
C
C-----------------------------------------------------------------------
C.. COEFFICIENTS FOR SECOND EQUATION: COVARIANT-2-COMP. OF EQ. OF MOTION
C..                                  (KYV2 = 1, DEFINED ON HALF MESH)
C-----------------------------------------------------------------------
C
      IEXE = -1

      DO 20 I=1,NR
      K = I
C
      INCLUDE 'tophat.inc'
C
C.....EQUATION MULTIPLIED BY 1/S.
C
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      Z1M  = (CS(I  )/CSM(I))**IEXJ1
      Z1P  = (CS(I+1)/CSM(I))**IEXJ1
C
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
      FSUBM(KYV2+LYROW,KXB1+LXCOL,k) =
     &     GF(ZEM*J3b1(I,  MSB),J3b1M(I,MSB))
      GSUBM(KYV2+LYROW,KXB1+LXCOL,k) =
     &     GF(ZEP*J3b1(I+1,MSB),J3b1M(I,MSB))
      ELSE
      FSUBM(KYV2+LYROW,KXB1+LXCOL,k) =
     &     GF(ZEM*J3b1(I,  MSB)*T(I),J3b1M(I,MSB)*TM(I))
      GSUBM(KYV2+LYROW,KXB1+LXCOL,k) =
     &     GF(ZEP*J3b1(I+1,MSB)*T(I+1),J3b1M(I,MSB)*TM(I))
      ENDIF

      DSUBM(KYV2+LYROW,KYJ1+LYCOL,k)= -GG(JACOBM(I,MSB),
     &     ZEM*Z1M*JACOBI(I,MSB),ZEP*Z1P*JACOBI(I+1,MSB))
      DSUBM(KYV2+LYROW,KYPR+LYCOL,k) =
     &     -  CMA*GG(TB2M(I,MSB),ZEM*TB2(I,MSB),ZEP*TB2(I+1,MSB))
     &     +  CNA*GG(G22B2B2M(I,MSB),ZEM*G22B2B2(I,MSB),
     &        ZEP*G22B2B2(I+1,MSB))
      IF (KYPE.GT.0.AND.INCKIN.EQ.0)
     &DSUBM(KYV2+LYROW,KYPE+LYCOL,k) =
     &     -  CMA*GG(TB2M(I,MSB),ZEM*TB2(I,MSB),ZEP*TB2(I+1,MSB))
     &     +  CNA*GG(G22B2B2M(I,MSB),ZEM*G22B2B2(I,MSB),
     &        ZEP*G22B2B2(I+1,MSB))
C-----------------------------------------------------------------------
C.. INERTIAL TERMS
C-----------------------------------------------------------------------
      FSUBM(KYV2+LYROW,KXV1+LXCOL,k) = FSUBM(KYV2+LYROW,KXV1+LXCOL,k)
     &     + DINERT*GF(SHIFTVC(I)*ZEM*RGV1G12(I,MSB),
     &       SHIFTVM(I)*RGV1G12M(I,MSB))
      GSUBM(KYV2+LYROW,KXV1+LXCOL,k) = GSUBM(KYV2+LYROW,KXV1+LXCOL,k)
     &     + DINERT*GF(SHIFTVC(I+1)*ZEP*RGV1G12(I+1,MSB),
     &       SHIFTVM(I)*RGV1G12M(I,MSB))
      DSUBM(KYV2+LYROW,KYV2+LYCOL,k) = DSUBM(KYV2+LYROW,KYV2+LYCOL,k)  
     &     + DINERT*GG(SHIFTVM(I)*RGV2G22M(I,MSB),
     &     ZEM*ZV2M*SHIFTVC(I)*RGV2G22(I,MSB),
     &          ZEP*ZV2P*SHIFTVC(I+1)*RGV2G22(I+1,MSB))

C     LLI
      IF (NPROFRC.GT.0) THEN 
      FSUBM(KYV2+LYROW,KXV1+LXCOL,k) = FSUBM(KYV2+LYROW,KXV1+LXCOL,k)
     &     - CNA*GF(ZEM*RGV1G12C(I,MSB),RGV1G12CM(I,MSB))
      GSUBM(KYV2+LYROW,KXV1+LXCOL,k) = GSUBM(KYV2+LYROW,KXV1+LXCOL,k)
     &     - CNA*GF(ZEP*RGV1G12C(I+1,MSB),RGV1G12CM(I,MSB))

      DSUBM(KYV2+LYROW,KYV2+LYCOL,k) = DSUBM(KYV2+LYROW,KYV2+LYCOL,k)  
     &     - CNA*GG(RGV2G22CM(I,MSB),ZEM*ZV2M*RGV2G22C(I,MSB),
     &              ZEM*ZV2P*RGV2G22C(I+1,MSB))
      ENDIF
20    CONTINUE
C

C-----------------------------------------------------------------------
C.. INERTIAL TERMS DUE TO PLASMA ROTATION
C-----------------------------------------------------------------------
      DO 22 I=1,NR
      K = I
      INCLUDE 'tophat.inc'
C
C.....EQUATION MULTIPLIED BY 1/S.
C
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      Z1M  = (CS(I  )/CSM(I))**IEXJ1
      Z1P  = (CS(I+1)/CSM(I))**IEXJ1
 
      CTMP1 = IRXY(I,MSB)*ROT(I)
      CTMP2 = IRXY(I+1,MSB)*ROT(I+1)
      CTMP3 = IRXYM(I,MSB)*ROTM(I)
      CTMP4 = IDRYX(I,MSB)*DROT(I)
      CTMP5 = IDRYX(I+1,MSB)*DROT(I+1)
      CTMP6 = IDRYXM(I,MSB)*DROTM(I)
      IF (IDIAMV.GE.3) THEN
         CTMP1 = IRXY(I,MSB)*(ROT(I)+OMEGASI(I)*FDIAMV)
         CTMP2 = IRXY(I+1,MSB)*(ROT(I+1)+OMEGASI(I+1)*FDIAMV)
         CTMP3 = IRXYM(I,MSB)*(ROTM(I)+OMEGASIM(I)*FDIAMV)
         CTMP4 = IDRYX(I,MSB)*(DROT(I)+DOMEGASI(I)*FDIAMV)
         CTMP5 = IDRYX(I+1,MSB)*(DROT(I+1)+DOMEGASI(I+1)*FDIAMV)
         CTMP6 = IDRYXM(I,MSB)*(DROTM(I)+DOMEGASIM(I)*FDIAMV)
      ELSEIF (IDIAMV.GE.2) THEN
         CTMP1 = IRXY(I,MSB)*(ROT(I)+OMEGASI(I)*.5*FDIAMV)
         CTMP2 = IRXY(I+1,MSB)*(ROT(I+1)+OMEGASI(I+1)*.5*FDIAMV)
         CTMP3 = IRXYM(I,MSB)*(ROTM(I)+OMEGASIM(I)*.5*FDIAMV)
      ENDIF

      IF (INERT2) THEN
      FSUBM(KYV2+LYROW,KXV1+LXCOL,I) = FSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &                           +GF(ZEM*CTMP4,CTMP6)
      GSUBM(KYV2+LYROW,KXV1+LXCOL,I) = GSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &                           +GF(ZEP*CTMP5,CTMP6)
C     LLI
      IF (NPROFRC.GT.0) THEN 
      FSUBM(KYV2+LYROW,KXV1+LXCOL,I) = FSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &     +GF(ZEM*RGV2RV1C(I,MSB),RGV2RV1CM(I,MSB))
      GSUBM(KYV2+LYROW,KXV1+LXCOL,I) = GSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &     +GF(ZEP*RGV2RV1C(I+1,MSB),RGV2RV1CM(I,MSB))

      DSUBM(KYV2+LYROW,KYV2+LYCOL,k) = DSUBM(KYV2+LYROW,KYV2+LYCOL,k)  
     &     +GG(RGV2RV2CM(I,MSB),ZEM*ZV2M*RGV2RV2C(I,MSB),
     &                          ZEM*ZV2P*RGV2RV2C(I+1,MSB))

      IF (KYV3.GT.0)
     &DSUBM(KYV2+LYROW,KYV3+LYCOL,I) = DSUBM(KYV2+LYROW,KYV3+LYCOL,I) 
     &     +GG(RGV2RV3CM(I,MSB),ZEM*RGV2RV3C(I,MSB),
     &                          ZEP*RGV2RV3C(I+1,MSB))

      ENDIF

C     GLX------<
      IF (NPROFRP.GT.0) THEN 
      FSUBM(KYV2+LYROW,KXV1+LXCOL,I) = FSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &   - CMA*GF(ZEM*RGV1RV3P(I,MSB),RGV1RV3PM(I,MSB))
     &   - CNA*GF(ZEM*RGV1RV4P(I,MSB),RGV1RV4PM(I,MSB))
     &   -     GF(ZEM*RGV2RV3P(I,MSB),RGV2RV3PM(I,MSB))
      GSUBM(KYV2+LYROW,KXV1+LXCOL,I) = GSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &   - CMA*GF(ZEP*RGV1RV3P(I+1,MSB),RGV1RV3PM(I,MSB))
     &   - CNA*GF(ZEP*RGV1RV4P(I+1,MSB),RGV1RV4PM(I,MSB))
     &   -     GF(ZEP*RGV2RV3P(I+1,MSB),RGV2RV3PM(I,MSB))

      DSUBM(KYV2+LYROW,KYV2+LYCOL,k) = DSUBM(KYV2+LYROW,KYV2+LYCOL,k)  
     &   - CMA*GG(RGV2RV1PM(I,MSB),ZEM*ZV2M*RGV2RV1P(I,MSB),
     &                             ZEP*ZV2P*RGV2RV1P(I+1,MSB))
     &   - CNA*GG(RGV2RV2PM(I,MSB),ZEM*ZV2M*RGV2RV2P(I,MSB),
     &                             ZEP*ZV2P*RGV2RV2P(I+1,MSB))
     &   -     GG(RGV2RV4PM(I,MSB),ZEM*ZV2M*RGV2RV4P(I,MSB),
     &                             ZEP*ZV2P*RGV2RV4P(I+1,MSB))

      IF (KYV3.GT.0)
     &DSUBM(KYV2+LYROW,KYV3+LYCOL,I) = DSUBM(KYV2+LYROW,KYV3+LYCOL,I) 
     &   -     GG(RGV2RV5PM(I,MSB),ZEM*RGV2RV5P(I,MSB),
     &                             ZEP*RGV2RV5P(I+1,MSB))
      ENDIF
C     --------->
      ENDIF

      IF (INERT1) THEN
      FSUBM(KYV2+LYROW,KXV1+LXCOL,I) = FSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &                           +GF(ZEM*CTMP1,CTMP3)
      GSUBM(KYV2+LYROW,KXV1+LXCOL,I) = GSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &                           +GF(ZEP*CTMP2,CTMP3)

C     LLI
      IF (NPROFRC.GT.0) THEN
      FSUBM(KYV2+LYROW,KXV1+LXCOL,I) = FSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &                         +GF(ZEM*RGIRXYC(I,MSB),RGIRXYCM(I,MSB))
      GSUBM(KYV2+LYROW,KXV1+LXCOL,I) = GSUBM(KYV2+LYROW,KXV1+LXCOL,I)
     &                         +GF(ZEP*RGIRXYC(I+1,MSB),RGIRXYCM(I,MSB))
      ENDIF

      CTMP4 = IRYZ(I,MSB)*ROT(I)
      CTMP5 = IRYZ(I+1,MSB)*ROT(I+1)
      CTMP6 = IRYZM(I,MSB)*ROTM(I)
      IF (KYV3.GT.0) THEN
      DSUBM(KYV2+LYROW,KYV3+LYCOL,I) = DSUBM(KYV2+LYROW,KYV3+LYCOL,I) 
     &   - GG(CTMP6,ZEM*CTMP4,ZEP*CTMP5)
C     LLI
      IF (NPROFRC.GT.0) THEN
      DSUBM(KYV2+LYROW,KYV3+LYCOL,I) = DSUBM(KYV2+LYROW,KYV3+LYCOL,I) 
     &   - GG(RGIRYZCM(I,MSB),ZEM*RGIRYZC(I,MSB),ZEP*RGIRYZC(I+1,MSB))
      ENDIF
      ENDIF
      ENDIF
     
      IF (KYRHOP.GT.0) THEN 
      DSUBM(KYV2+LYROW,KYRHOP+LYCOL,I) =DSUBM(KYV2+LYROW,KYRHOP+LYCOL,I)
     &   - GG(ROTDSDZM(I,MSB),ZEM*ROTDSDZ(I,MSB),ZEP*ROTDSDZ(I+1,MSB))

C     LLI
      IF (NPROFRC.GT.0) THEN
      DSUBM(KYV2+LYROW,KYRHOP+LYCOL,I) =DSUBM(KYV2+LYROW,KYRHOP+LYCOL,I)
     &   - GG(RGROTDSDZCM(I,MSB),ZEM*RGROTDSDZC(I,MSB),
     &                           ZEP*RGROTDSDZC(I+1,MSB))
      ENDIF
      ENDIF
 
 22   CONTINUE
 
 
      IF (KJPKEY.EQ.0) GOTO 150

C-----------------------------------------------------------------------
C.. COEFFICIENTS FOR THIRD  EQUATION: COVARIANT-3-COMP. OF EQ. OF MOTION
C..                                  (KYV3 = 2, DEFINED ON HALF MESH)
C-----------------------------------------------------------------------
C
      IEXE = -1
C
      IF (KYV3.GT.0) THEN
      DO 30 I=1,NR
C
C.....EQUATION MULTIPLIED BY 1/S.
C
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      Z1M  = (CS(I  )/CSM(I))**IEXJ1
      Z1P  = (CS(I+1)/CSM(I))**IEXJ1
C
      INCLUDE 'tophat.inc'
C-----------------------------------------------------------------------

      DSUBM(KYV3+LYROW,KYPR+LYCOL,I)=
     &     - CNA*GG(B3j2M(I,MSB)*TM(I),ZEM*B3J2(I,MSB)*T(I),
     &       ZEP*B3J2(I+1,MSB)*T(I+1))
      IF (KYPE.GT.0.AND.INCKIN.EQ.0)
     &DSUBM(KYV3+LYROW,KYPE+LYCOL,I)=
     &     - CNA*GG(B3j2M(I,MSB)*TM(I),ZEM*B3J2(I,MSB)*T(I),
     &       ZEP*B3J2(I+1,MSB)*T(I+1))

      IF (KJRER.EQ.6.AND.ETA.GT.0.) THEN
      CTMP1 = JRE_EXB*OMEGACI0*RESIST(I)
      CTMP2 = JRE_EXB*OMEGACI0*RESIST(I+1)
      CTMP3 = JRE_EXB*OMEGACI0*RESISM(I)
      IF (KXJ2L.GT.0.AND.KXJRE2L.GT.0) THEN
      FSUBM(KYV3+LYROW,KXJ2L+LXCOL,I)=FSUBM(KYV3+LYROW,KXJ2L+LXCOL,I) 
     &     +GF(ZEM*FRENRE(I,MSB)*CTMP1*DPSIDS(I),
     &             FRENREM(I,MSB)*CTMP3*DPSIDSM(I))
      GSUBM(KYV3+LYROW,KXJ2L+LXCOL,I)=GSUBM(KYV3+LYROW,KXJ2L+LXCOL,I)
     &     +GF(ZEP*FRENRE(I+1,MSB)*CTMP2*DPSIDS(I+1),
     &             FRENREM(I,MSB)*CTMP3*DPSIDSM(I))
      FSUBM(KYV3+LYROW,KXJRE2L+LXCOL,I)=
     &FSUBM(KYV3+LYROW,KXJRE2L+LXCOL,I)
     &     -GF(ZEM*FRENRE(I,MSB)*CTMP1*DPSIDS(I),
     &             FRENREM(I,MSB)*CTMP3*DPSIDSM(I))
      GSUBM(KYV3+LYROW,KXJRE2L+LXCOL,I)=
     &GSUBM(KYV3+LYROW,KXJRE2L+LXCOL,I)
     &     -GF(ZEP*FRENRE(I+1,MSB)*CTMP2*DPSIDS(I+1),
     &             FRENREM(I,MSB)*CTMP3*DPSIDSM(I))
      ENDIF

      FSUBM(KYV3+LYROW,KXJ3+LXCOL,I)=FSUBM(KYV3+LYROW,KXJ3+LXCOL,I) 
     &     +GF(ZEM*FRENRE(I,MSB)*CTMP1*T(I),
     &             FRENREM(I,MSB)*CTMP3*TM(I))
      GSUBM(KYV3+LYROW,KXJ3+LXCOL,I)=GSUBM(KYV3+LYROW,KXJ3+LXCOL,I)
     &     +GF(ZEP*FRENRE(I+1,MSB)*CTMP2*T(I+1),
     &             FRENREM(I,MSB)*CTMP3*TM(I))
      FSUBM(KYV3+LYROW,KXJRE3+LXCOL,I)=
     &FSUBM(KYV3+LYROW,KXJRE3+LXCOL,I)
     &     -GF(ZEM*FRENRE(I,MSB)*CTMP1*T(I),
     &             FRENREM(I,MSB)*CTMP3*TM(I))
      GSUBM(KYV3+LYROW,KXJRE3+LXCOL,I)=
     &GSUBM(KYV3+LYROW,KXJRE3+LXCOL,I)
     &     -GF(ZEP*FRENRE(I+1,MSB)*CTMP2*T(I+1),
     &             FRENREM(I,MSB)*CTMP3*TM(I))

      IF (KXJRE.GT.0) THEN
      FSUBM(KYV3+LYROW,KXJRE+LXCOL,I)=
     &FSUBM(KYV3+LYROW,KXJRE+LXCOL,I)
     &     -GF(ZEM*FRENREJB(I,MSB)*CTMP1,FRENREJBM(I,MSB)*CTMP3)
      GSUBM(KYV3+LYROW,KXJRE+LXCOL,I)=
     &GSUBM(KYV3+LYROW,KXJRE+LXCOL,I)
     &     -GF(ZEP*FRENREJB(I+1,MSB)*CTMP2,FRENREJBM(I,MSB)*CTMP3)
      ENDIF

      IF (KXB2L.GT.0) THEN
      FSUBM(KYV3+LYROW,KXB2L+LXCOL,I)=
     &FSUBM(KYV3+LYROW,KXB2L+LXCOL,I)
     &     -GF(ZEM*FREJNREDB(I,MSB)*CTMP1*DPSIDS(I),
     &             FREJNREDBM(I,MSB)*CTMP3*DPSIDSM(I))
      GSUBM(KYV3+LYROW,KXB2L+LXCOL,I)=
     &GSUBM(KYV3+LYROW,KXB2L+LXCOL,I)
     &     -GF(ZEP*FREJNREDB(I+1,MSB)*CTMP2*DPSIDS(I+1),
     &             FREJNREDBM(I,MSB)*CTMP3*DPSIDSM(I))
      ENDIF

      DSUBM(KYV3+LYROW,KYB3+LYCOL,I)=
     &DSUBM(KYV3+LYROW,KYB3+LYCOL,I)
     &     -GG(FREJNREDBM(I,MSB)*CTMP3*TM(I),
     &         ZEM*FREJNREDB(I,MSB)*CTMP1*T(I),
     &         ZEP*FREJNREDB(I+1,MSB)*CTMP2*T(I+1))
      ENDIF
C
C-----------------------------------------------------------------------
C.. INERTIA AND COLLISIONLESS PARALLEL 'VISCOSITY'
C-----------------------------------------------------------------------
C
      DSUBM(KYV3+LYROW,KYV3+LYCOL,I) = (SHIFTM(I)-NUII)
     &       * DINERT*GG(RGV3G33M(I,MSB),ZEM*RGV3G33(I,MSB),
     &                                   ZEP*RGV3G33(I+1,MSB))

      IF (NPROFRC.GT.0) 
     &DSUBM(KYV3+LYROW,KYV3+LYCOL,I) = DSUBM(KYV3+LYROW,KYV3+LYCOL,I) -
     &       CNA*GG(RGV3G33CM(I,MSB),ZEM*RGV3G33C(I,MSB),
     &                               ZEP*RGV3G33C(I+1,MSB))

      IF (ABS(PVISC).GT.0..AND.IVISC.GT.0) THEN
      ZKPAVT = ABS(RM(MSA,NSA)/QM(I) + RN(NSA))
     &                 * SQRT(REAL(PEQM(I))/RHOM(I))
C
C     SQRT(2) IS INCLUDED IN PRESSURE = PI + PE
C
C  THE FOLLOWING MODIFIED TO TAKE INTO ACCOUNT DIFFERENCE BETWEEN 
C     EULERIAN AND LAGRANGIAN FRAME
      DSUBM(KYV3+LYROW,KYV3+LYCOL,I) = DSUBM(KYV3+LYROW,KYV3+LYCOL,I)
     &       - PVISC*ZKPAVT
     &       * GG(RGV3G33M(I,MSB),ZEM*RGV3G33(I,MSB),
     &         ZEP*RGV3G33(I+1,MSB))
C     THE FOLLOWING FOR LAGRANGIAN FRAME
      IF (IVISC.EQ.1) THEN
      CTMP4 = IRYZ(I,MSB)*ROT(I)
      CTMP5 = IRYZ(I+1,MSB)*ROT(I+1)
      CTMP6 = IRYZM(I,MSB)*ROTM(I)
      IF (KYX2.GT.0) THEN
      DSUBM(KYV3+LYROW,KYX2+LYCOL,I) = DSUBM(KYV3+LYROW,KYX2+LYCOL,I)
     &   + PVISC*ZKPAVT*0.5*
     &     GG(CTMP6,ZV2M*ZEM*CTMP4,ZV2P*ZEP*CTMP5)
C     LLI
      IF (NPROFRC.GT.0) THEN
      DSUBM(KYV3+LYROW,KYX2+LYCOL,I)=DSUBM(KYV3+LYROW,KYX2+LYCOL,I) 
     &            -PVISC*ZKPAVT
     &            *GG(RGV3RV2CM(I,MSB),ZV2M*ZEM*RGV3RV2C(I,MSB),
     &                ZV2P*ZEP*RGV3RV2C(I+1,MSB))
     &            +PVISC*ZKPAVT*0.5
     &            *GG(RGIRYZCM(I,MSB),ZV2M*ZEM*RGIRYZC(I,MSB),
     &                ZV2P*ZEP*RGIRYZC(I+1,MSB))
      ENDIF
C     GLX------<
      IF (NPROFRP.GT.0) THEN
      DSUBM(KYV3+LYROW,KYX2+LYCOL,I)=DSUBM(KYV3+LYROW,KYX2+LYCOL,I) 
     &   - PVISC*ZKPAVT
     &          *GG(RGV2RV5PM(I,MSB),ZV2M*ZEM*RGV2RV5P(I,MSB),
     &                               ZV2P*ZEP*RGV2RV5P(I+1,MSB))*0.5
      ENDIF
C     --------->
      ENDIF
 
      IF (KXX1.GT.0) THEN
      CTMP1 = -RHO(I)*DROT(I)*T(I)
      CTMP2 = -RHO(I+1)*DROT(I+1)*T(I+1)
      CTMP3 = -RHOM(I)*DROTM(I)*TM(I)
      CTMP4 = IRXZ(I,MSB)*ROT(I) 
      CTMP5 = IRXZ(I+1,MSB)*ROT(I+1) 
      CTMP6 = IRXZM(I,MSB)*ROTM(I) 
      FSUBM(KYV3+LYROW,KXX1+LXCOL,I)=FSUBM(KYV3+LYROW,KXX1+LXCOL,I)
     &   + PVISC*ZKPAVT*(-GF(ZEM*CTMP4,CTMP6)*0.5
     &   + GF(ZEM*CTMP1*JACOBI(I,MSB),CTMP3*JACOBM(I,MSB)))
      GSUBM(KYV3+LYROW,KXX1+LXCOL,I)=GSUBM(KYV3+LYROW,KXX1+LXCOL,I)
     &   + PVISC*ZKPAVT*(-GF(ZEP*CTMP5,CTMP6)*0.5
     &   + GF(ZEP*CTMP2*JACOBI(I+1,MSB),CTMP3*JACOBM(I,MSB)))

C     LLI
      IF (NPROFRC.GT.0) THEN
      FSUBM(KYV3+LYROW,KXX1+LXCOL,I)=FSUBM(KYV3+LYROW,KXX1+LXCOL,I) 
     &   -PVISC*ZKPAVT*GF(ZEM*RGV3RV1C(I,MSB),RGV3RV1CM(I,MSB))
      GSUBM(KYV3+LYROW,KXX1+LXCOL,I)=GSUBM(KYV3+LYROW,KXX1+LXCOL,I)
     &   -PVISC*ZKPAVT*GF(ZEP*RGV3RV1C(I+1,MSB),RGV3RV1CM(I,MSB))

      FSUBM(KYV3+LYROW,KXX1+LXCOL,I)=FSUBM(KYV3+LYROW,KXX1+LXCOL,I) 
     &   -0.5*PVISC*ZKPAVT*GF(ZEM*RGIRXZC(I,MSB),RGIRXZCM(I,MSB))
      GSUBM(KYV3+LYROW,KXX1+LXCOL,I)=GSUBM(KYV3+LYROW,KXX1+LXCOL,I)
     &   -0.5*PVISC*ZKPAVT*GF(ZEP*RGIRXZC(I+1,MSB),RGIRXZCM(I,MSB))
      ENDIF
C     GLX------<
      IF (NPROFRP.GT.0) THEN
      FSUBM(KYV3+LYROW,KXX1+LXCOL,I)=FSUBM(KYV3+LYROW,KXX1+LXCOL,I) 
     &   - PVISC*ZKPAVT*GF(ZEM*RGV3RV3P(I,MSB),RGV3RV3PM(I,MSB))
      GSUBM(KYV3+LYROW,KXX1+LXCOL,I)=GSUBM(KYV3+LYROW,KXX1+LXCOL,I)
     &   - PVISC*ZKPAVT*GF(ZEP*RGV3RV3P(I+1,MSB),RGV3RV3PM(I,MSB))
      ENDIF
C     --------->
      ENDIF

C     LLI      
      IF (KYX3.GT.0.AND.NPROFRC.GT.0) THEN
      CTMP1 = RHOM(I)*TM(I)*DPSIDSM(I)*RGDRCDCCM(I,MSB)
      CTMP2 = RHO(I)*T(I)*DPSIDS(I)*RGDRCDCC(I,MSB)
      CTMP3 = RHO(I+1)*T(I+1)*DPSIDS(I+1)*RGDRCDCC(I+1,MSB)
      DSUBM(KYV3+LYROW,KYX3+LYCOL,I) =DSUBM(KYV3+LYROW,KYX3+LYCOL,I)
     &   -PVISC*ZKPAVT*GG(CTMP1,ZEM*CTMP2,ZEP*CTMP3)
      ENDIF
C     GLX------<      
      IF (KYX3.GT.0.AND.NPROFRP.GT.0) THEN
      DSUBM(KYV3+LYROW,KYX3+LYCOL,I) =DSUBM(KYV3+LYROW,KYX3+LYCOL,I)
     &   - PVISC*ZKPAVT*CMB*GG(RGV3RV1PM(I,MSB),ZEM*RGV3RV1P(I,MSB),
     &                                    ZEP*RGV3RV1P(I+1,MSB))*0.5
      ENDIF
C     --------->
      ENDIF
      ENDIF
          
 30   CONTINUE
C
C-----------------------------------------------------------------------
C.. INERTIA TERMS DUE TO PLASMA ROTATION
C-----------------------------------------------------------------------
C
      DO 32 I=1,NR
C
C.....EQUATION MULTIPLIED BY 1/S.
C
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      Z1M  = (CS(I  )/CSM(I))**IEXJ1
      Z1P  = (CS(I+1)/CSM(I))**IEXJ1
 
      INCLUDE 'tophat.inc'
C-----------------------------------------------------------------------
      IF (INERT1) THEN
      CTMP4 = IRYZ(I,MSB)*ROT(I) 
      CTMP5 = IRYZ(I+1,MSB)*ROT(I+1) 
      CTMP6 = IRYZM(I,MSB)*ROTM(I) 
      IF (IDIAMV.GE.3) THEN
         CTMP4 = IRYZ(I,MSB)*(ROT(I)+OMEGASI(I)*FDIAMV) 
         CTMP5 = IRYZ(I+1,MSB)*(ROT(I+1)+OMEGASI(I+1)*FDIAMV) 
         CTMP6 = IRYZM(I,MSB)*(ROTM(I)+OMEGASIM(I)*FDIAMV) 
      ELSEIF (IDIAMV.GE.2) THEN
         CTMP4 = IRYZ(I,MSB)*(ROT(I)+OMEGASI(I)*.5*FDIAMV) 
         CTMP5 = IRYZ(I+1,MSB)*(ROT(I+1)+OMEGASI(I+1)*.5*FDIAMV) 
         CTMP6 = IRYZM(I,MSB)*(ROTM(I)+OMEGASIM(I)*.5*FDIAMV) 
      ENDIF
      DSUBM(KYV3+LYROW,KYV2+LYCOL,I) = 
     &   GG(CTMP6,ZV2M*ZEM*CTMP4,ZV2P*ZEP*CTMP5)

C     LLI
      IF (NPROFRC.GT.0) THEN
      DSUBM(KYV3+LYROW,KYV2+LYCOL,I)=DSUBM(KYV3+LYROW,KYV2+LYCOL,I) 
     &                +GG(RGIRYZCM(I,MSB),ZV2M*ZEM*RGIRYZC(I,MSB),
     &                    ZV2P*ZEP*RGIRYZC(I+1,MSB))
      ENDIF
      ENDIF
     
      IF (KYRHOP.GT.0.AND.KYV3.GT.0) THEN
      DSUBM(KYV3+LYROW,KYRHOP+LYCOL,I) =DSUBM(KYV3+LYROW,KYRHOP+LYCOL,I)
     &  -GG(ROT2DSDZM(I,MSB),ZEM*ROT2DSDZ(I,MSB),ZEP*ROT2DSDZ(I+1,MSB))

C     LLI
      IF (NPROFRC.GT.0) THEN
      DSUBM(KYV3+LYROW,KYRHOP+LYCOL,I) =DSUBM(KYV3+LYROW,KYRHOP+LYCOL,I)
     &  -GG(RGROT2DSDZCM(I,MSB),ZEM*RGROT2DSDZC(I,MSB),
     &                          ZEP*RGROT2DSDZC(I+1,MSB))
      ENDIF
      ENDIF
 
      CTMP1 = -RHO(I)*DROT(I)*T(I)
      CTMP2 = -RHO(I+1)*DROT(I+1)*T(I+1)
      CTMP3 = -RHOM(I)*DROTM(I)*TM(I)
      CTMP4 = IRXZ(I,MSB)*ROT(I) 
      CTMP5 = IRXZ(I+1,MSB)*ROT(I+1) 
      CTMP6 = IRXZM(I,MSB)*ROTM(I) 
      IF (IDIAMV.GE.3) THEN
         CTMP1 = -RHO(I)*(DROT(I)+DOMEGASI(I)*FDIAMV)*T(I)
         CTMP2 = -RHO(I+1)*(DROT(I+1)+DOMEGASI(I+1)*FDIAMV)*T(I+1)
         CTMP3 = -RHOM(I)*(DROTM(I)+DOMEGASIM(I)*FDIAMV)*TM(I)
         CTMP4 = IRXZ(I,MSB)*(ROT(I)+OMEGASI(I)*FDIAMV) 
         CTMP5 = IRXZ(I+1,MSB)*(ROT(I+1)+OMEGASI(I+1)*FDIAMV) 
         CTMP6 = IRXZM(I,MSB)*(ROTM(I)+OMEGASIM(I)*FDIAMV) 
      ELSEIF (IDIAMV.GE.2) THEN
         CTMP4 = IRXZ(I,MSB)*(ROT(I)+OMEGASI(I)*.5*FDIAMV) 
         CTMP5 = IRXZ(I+1,MSB)*(ROT(I+1)+OMEGASI(I+1)*.5*FDIAMV) 
         CTMP6 = IRXZM(I,MSB)*(ROTM(I)+OMEGASIM(I)*.5*FDIAMV) 
      ENDIF

      IF (INERT1) THEN
      FSUBM(KYV3+LYROW,KXV1+LXCOL,I)=FSUBM(KYV3+LYROW,KXV1+LXCOL,I) 
     &            -GF(ZEM*CTMP4,CTMP6)
      GSUBM(KYV3+LYROW,KXV1+LXCOL,I)=GSUBM(KYV3+LYROW,KXV1+LXCOL,I)
     &            -GF(ZEP*CTMP5,CTMP6)
C     LLI
      IF (NPROFRC.GT.0) THEN
      FSUBM(KYV3+LYROW,KXV1+LXCOL,I)=FSUBM(KYV3+LYROW,KXV1+LXCOL,I) 
     &            -GF(ZEM*RGIRXZC(I,MSB),RGIRXZCM(I,MSB))
      GSUBM(KYV3+LYROW,KXV1+LXCOL,I)=GSUBM(KYV3+LYROW,KXV1+LXCOL,I)
     &            -GF(ZEP*RGIRXZC(I+1,MSB),RGIRXZCM(I,MSB))
      ENDIF
      ENDIF

      IF (INERT2) THEN
      FSUBM(KYV3+LYROW,KXV1+LXCOL,I)=FSUBM(KYV3+LYROW,KXV1+LXCOL,I) 
     &            +GF(ZEM*CTMP1*JACOBI(I,MSB),CTMP3*JACOBM(I,MSB))
      GSUBM(KYV3+LYROW,KXV1+LXCOL,I)=GSUBM(KYV3+LYROW,KXV1+LXCOL,I)
     &            +GF(ZEP*CTMP2*JACOBI(I+1,MSB),CTMP3*JACOBM(I,MSB))
C     LLI
      IF (NPROFRC.GT.0) THEN
      FSUBM(KYV3+LYROW,KXV1+LXCOL,I)=FSUBM(KYV3+LYROW,KXV1+LXCOL,I) 
     &            -GF(ZEM*RGV3RV1C(I,MSB),RGV3RV1CM(I,MSB))
      GSUBM(KYV3+LYROW,KXV1+LXCOL,I)=GSUBM(KYV3+LYROW,KXV1+LXCOL,I)
     &            -GF(ZEP*RGV3RV1C(I+1,MSB),RGV3RV1CM(I,MSB))

      DSUBM(KYV3+LYROW,KYV2+LYCOL,I)=DSUBM(KYV3+LYROW,KYV2+LYCOL,I) 
     &            -GG(RGV3RV2CM(I,MSB),ZV2M*ZEM*RGV3RV2C(I,MSB),
     &                    ZV2P*ZEP*RGV3RV2C(I+1,MSB))

      CTMP1 = RHOM(I)*TM(I)*DPSIDSM(I)*RGDRCDCCM(I,MSB)
      CTMP2 = RHO(I)*T(I)*DPSIDS(I)*RGDRCDCC(I,MSB)
      CTMP3 = RHO(I+1)*T(I+1)*DPSIDS(I+1)*RGDRCDCC(I+1,MSB)
      DSUBM(KYV3+LYROW,KYV3+LYCOL,I) =DSUBM(KYV3+LYROW,KYV3+LYCOL,I)
     &            -GG(CTMP1,ZEM*CTMP2,ZEP*CTMP3)
      ENDIF
C     GLX------<
      IF (NPROFRP.GT.0) THEN
      CTMP4 = ROTP(I)*JACOBI(I,MSB)*PPEQ(I)*DPSIDS(I) -
     &        RGV3G33(I,MSB)*DRHOU(I)
      CTMP5 = ROTPM(I)*JACOBM(I,MSB)*PPEQM(I)*DPSIDSM(I) -
     &        RGV3G33M(I,MSB)*DRHOUM(I)
      CTMP6 = ROTP(I+1)*JACOBI(I+1,MSB)*PPEQ(I+1)*DPSIDS(I+1) -
     &        RGV3G33(I+1,MSB)*DRHOU(I+1)
      FSUBM(KYV3+LYROW,KXV1+LXCOL,I)=FSUBM(KYV3+LYROW,KXV1+LXCOL,I) 
     &   + GF(ZEM*CTMP4,CTMP5)
      GSUBM(KYV3+LYROW,KXV1+LXCOL,I)=GSUBM(KYV3+LYROW,KXV1+LXCOL,I)
     &   + GF(ZEP*CTMP6,CTMP5)

      DSUBM(KYV3+LYROW,KYV3+LYCOL,I) =DSUBM(KYV3+LYROW,KYV3+LYCOL,I)
     &   - CMROW*GG(RGV3RV1PM(I,MSB),ZEM*RGV3RV1P(I,MSB),
     &                               ZEP*RGV3RV1P(I+1,MSB))
     &   - CNA*  GG(RGV3RV2PM(I,MSB),ZEM*RGV3RV2P(I,MSB),
     &                               ZEP*RGV3RV2P(I+1,MSB))
      ENDIF
C     --------->
      ENDIF
 32   CONTINUE
 
      ENDIF
C
C-----------------------------------------------------------------------
C.. FOURTH EQUATION: CONTRAVARIANT-1-COMP. OF FARADAY
C..                                (KXB1 = 2, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
C
      DO 40 I=2,NRP1
      K = I
C
      INCLUDE 'tent.inc'
C
      IF (ABS(RM(MROW,2)).GT.0.1) THEN
         BsubM(kxB1+lxrow, kxV1+lxcol,k)=
     $     cna*FF(B3j2(i,msb)*T(I),B3j2M(i-1,msb)*TM(I-1), 
     $     B3j2M(i,msb)*TM(I))
         AsubM(kxB1+lxrow, kxV1+lxcol,k)=
     &     cna*FFM(B3j2M(i-1,msb)*TM(I-1))
         CsubM(kxB1+lxrow, kxV1+lxcol,k)=
     &     cna*FFP(B3j2M(i,msb)*TM(I))

         BsubM(kxB1+lxrow, kxJ3+lxcol,k)=
     $        -cmrow*FF(resist(i)*DG33L(i,msb),
     $        resisM(i-1)*DG33LM(i-1,msb), resisM(i)*DG33LM(i,msb))
         AsubM(kxB1+lxrow, kxJ3+lxcol,k)=
     $        -cmrow*FFM(resisM(i-1)*DG33LM(i-1,msb))
         CsubM(kxB1+lxrow, kxJ3+lxcol,k)=
     $        -cmrow*FFP(resisM(i)*DG33LM(i,msb))

         IF (KJRER.EQ.6) THEN
         CTMP1 = JRE_EQFRAC*JRE_EXB
         BsubM(kxB1+lxrow, KXJRE3+lxcol,k)=
     $        +cmrow*CTMP1*FF(resist(i)*DG33L(i,msb),
     $        resisM(i-1)*DG33LM(i-1,msb), resisM(i)*DG33LM(i,msb))
         AsubM(kxB1+lxrow, KXJRE3+lxcol,k)=
     $        +cmrow*CTMP1*FFM(resisM(i-1)*DG33LM(i-1,msb))
         CsubM(kxB1+lxrow, KXJRE3+lxcol,k)=
     $        +cmrow*CTMP1*FFP(resisM(i)*DG33LM(i,msb))
         ENDIF

         IF (KJRER.EQ.1.OR.KJRER.EQ.3.OR.KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &      + CMROW*FF(RESIST(I)*T(I)*FRE1B(I,MSB),
     &                 RESISM(I-1)*TM(I-1)*FRE1BM(I-1,MSB),
     &                 RESISM(I)*TM(I)*FRE1BM(I,MSB))*JRE_EQFRAC
         ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)=ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &      + CMROW*FFM(RESISM(I-1)*TM(I-1)*FRE1BM(I-1,MSB))*JRE_EQFRAC
         CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &      + CMROW*FFP(RESISM(I)*TM(I)*FRE1BM(I,MSB))*JRE_EQFRAC
         ENDIF

         IF (KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &     -CNA*FF(RESIST(I)*DPSIDS(I)*FREG22JB(I,MSB),
     &             RESISM(I-1)*DPSIDSM(I-1)*FREG22JBM(I-1,MSB),
     &             RESISM(I)*DPSIDSM(I)*FREG22JBM(I,MSB))*JRE_EQFRAC
         ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)=ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &     -CNA*FFM(RESISM(I-1)*DPSIDSM(I-1)*FREG22JBM(I-1,MSB))
     &      *JRE_EQFRAC
         CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &     -CNA*FFP(RESISM(I)*DPSIDSM(I)*FREG22JBM(I,MSB))*JRE_EQFRAC

         IF (KXB3L.GT.0) THEN
         BSUBM(KXB1+LXROW,KXB3L+LXCOL,K)=BSUBM(KXB1+LXROW,KXB3L+LXCOL,K)
     &     +CMA*JRE_EQFRAC*FF(RESIST(I)*FREJPB(I,MSB),
     &                        RESISM(I-1)*FREJPBM(I-1,MSB),
     &                        RESISM(I)*FREJPBM(I,MSB))
         ASUBM(KXB1+LXROW,KXB3L+LXCOL,K)=ASUBM(KXB1+LXROW,KXB3L+LXCOL,K)
     &     +CMA*JRE_EQFRAC*FFM(RESISM(I-1)*FREJPBM(I-1,MSB))
         CSUBM(KXB1+LXROW,KXB3L+LXCOL,K)=CSUBM(KXB1+LXROW,KXB3L+LXCOL,K)
     &     +CMA*JRE_EQFRAC*FFP(RESISM(I)*FREJPBM(I,MSB))
         ENDIF

         IF (KXB2L.GT.0) THEN
         BSUBM(KXB1+LXROW,KXB2L+LXCOL,K)=BSUBM(KXB1+LXROW,KXB2L+LXCOL,K)
     &     -CNA*JRE_EQFRAC*FF(RESIST(I)*FREJPB(I,MSB),
     &                        RESISM(I-1)*FREJPBM(I-1,MSB),
     &                        RESISM(I)*FREJPBM(I,MSB))
         ASUBM(KXB1+LXROW,KXB2L+LXCOL,K)=ASUBM(KXB1+LXROW,KXB2L+LXCOL,K)
     &     -CNA*JRE_EQFRAC*FFM(RESISM(I-1)*FREJPBM(I-1,MSB))
         CSUBM(KXB1+LXROW,KXB2L+LXCOL,K)=CSUBM(KXB1+LXROW,KXB2L+LXCOL,K)
     &     -CNA*JRE_EQFRAC*FFP(RESISM(I)*FREJPBM(I,MSB))
         ENDIF
         ENDIF
      ELSE
         BsubM(kxB1+lxrow, kxV1+lxcol,k)= 
     $    cna*FF(B3j2(i,msb),B3j2M(i-1,msb),B3j2M(i,msb))    
         AsubM(kxB1+lxrow, kxV1+lxcol,k)=cna*FFM(B3j2M(i-1,msb))    
         CsubM(kxB1+lxrow, kxV1+lxcol,k)=cna*FFP(B3j2M(i,msb))    

         BsubM(kxB1+lxrow, kxJ3+lxcol,k)=
     $        -cmrow*FF(resist(i)*DG33L(i,msb)/T(I),
     $        resisM(i-1)*DG33LM(i-1,msb)/TM(I-1), 
     $        resisM(i)*DG33LM(i,msb)/TM(I))
         AsubM(kxB1+lxrow, kxJ3+lxcol,k)=
     $        -cmrow*FFM(resisM(i-1)*DG33LM(i-1,msb)/TM(I-1))
         CsubM(kxB1+lxrow, kxJ3+lxcol,k)=
     $        -cmrow*FFP(resisM(i)*DG33LM(i,msb)/TM(I))

         IF (KJRER.EQ.6) THEN
         CTMP1 = JRE_EQFRAC*JRE_EXB
         BsubM(kxB1+lxrow, KXJRE3+lxcol,k)=
     $        +cmrow*CTMP1*FF(resist(i)*DG33L(i,msb)/T(I),
     $        resisM(i-1)*DG33LM(i-1,msb)/TM(I-1), 
     $        resisM(i)*DG33LM(i,msb)/TM(I))
         AsubM(kxB1+lxrow, KXJRE3+lxcol,k)=
     $        +cmrow*CTMP1*FFM(resisM(i-1)*DG33LM(i-1,msb)/TM(I-1))
         CsubM(kxB1+lxrow, KXJRE3+lxcol,k)=
     $        +cmrow*CTMP1*FFP(resisM(i)*DG33LM(i,msb)/TM(I))
         ENDIF

         IF (KJRER.EQ.1.OR.KJRER.EQ.3.OR.KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &      + CMROW*FF(RESIST(I)*FRE1B(I,MSB),
     &                 RESISM(I-1)*FRE1BM(I-1,MSB),
     &                 RESISM(I)*FRE1BM(I,MSB))*JRE_EQFRAC
         ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)=ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &      + CMROW*FFM(RESISM(I-1)*FRE1BM(I-1,MSB))*JRE_EQFRAC
         CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &      + CMROW*FFP(RESISM(I)*FRE1BM(I,MSB))*JRE_EQFRAC
         ENDIF

         IF (KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=BSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &     -CNA*FF(RESIST(I)*DPSIDS(I)/T(I)*FREG22JB(I,MSB),
     &             RESISM(I-1)*DPSIDSM(I-1)/TM(I-1)*FREG22JBM(I-1,MSB),
     &             RESISM(I)*DPSIDSM(I)/TM(I)*FREG22JBM(I,MSB))
     &      *JRE_EQFRAC
         ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)=ASUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &     -CNA*FFM(RESISM(I-1)*DPSIDSM(I-1)/TM(I-1)*FREG22JBM(I-1,MSB))
     &      *JRE_EQFRAC
         CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)=CSUBM(KXB1+LXROW,KXJRE+LXCOL,K)
     &     -CNA*FFP(RESISM(I)*DPSIDSM(I)/TM(I)*FREG22JBM(I,MSB))
     &      *JRE_EQFRAC

         IF (KXB3L.GT.0) THEN
         BSUBM(KXB1+LXROW,KXB3L+LXCOL,K)=BSUBM(KXB1+LXROW,KXB3L+LXCOL,K)
     &     +CMA*JRE_EQFRAC*FF(RESIST(I)/T(I)*FREJPB(I,MSB),
     &                        RESISM(I-1)/TM(I-1)*FREJPBM(I-1,MSB),
     &                        RESISM(I)/TM(I)*FREJPBM(I,MSB))
         ASUBM(KXB1+LXROW,KXB3L+LXCOL,K)=ASUBM(KXB1+LXROW,KXB3L+LXCOL,K)
     &     +CMA*JRE_EQFRAC*FFM(RESISM(I-1)/TM(I-1)*FREJPBM(I-1,MSB))
         CSUBM(KXB1+LXROW,KXB3L+LXCOL,K)=CSUBM(KXB1+LXROW,KXB3L+LXCOL,K)
     &     +CMA*JRE_EQFRAC*FFP(RESISM(I)/TM(I)*FREJPBM(I,MSB))
         ENDIF

         IF (KXB2L.GT.0) THEN
         BSUBM(KXB1+LXROW,KXB2L+LXCOL,K)=BSUBM(KXB1+LXROW,KXB2L+LXCOL,K)
     &     -CNA*JRE_EQFRAC*FF(RESIST(I)/T(I)*FREJPB(I,MSB),
     &                        RESISM(I-1)/TM(I-1)*FREJPBM(I-1,MSB),
     &                        RESISM(I)/TM(I)*FREJPBM(I,MSB))
         ASUBM(KXB1+LXROW,KXB2L+LXCOL,K)=ASUBM(KXB1+LXROW,KXB2L+LXCOL,K)
     &     -CNA*JRE_EQFRAC*FFM(RESISM(I-1)/TM(I-1)*FREJPBM(I-1,MSB))
         CSUBM(KXB1+LXROW,KXB2L+LXCOL,K)=CSUBM(KXB1+LXROW,KXB2L+LXCOL,K)
     &     -CNA*JRE_EQFRAC*FFP(RESISM(I)/TM(I)*FREJPBM(I,MSB))
         ENDIF
         ENDIF
      ENDIF

      IF (NPROFRC.GT.0) THEN
      IF ((ABS(RM(MROW,2)).GT.0.1.AND.ABS(RM(MSA,2)).GT.0.1).OR.
     $    (ABS(RM(MROW,2)).LT.0.1.AND.ABS(RM(MSA,2)).LT.0.1)) THEN
         BsubM(kxB1+lxrow,kxB1+lxcol,k)=BsubM(kxB1+lxrow,kxB1+lxcol,k)
     $    -cna*FF(RGROTC(i,msb),RGROTCM(i-1,msb),RGROTCM(i,msb))
         AsubM(kxB1+lxrow,kxB1+lxcol,k)=AsubM(kxB1+lxrow,kxB1+lxcol,k)
     &    -cna*FFM(RGROTCM(i-1,msb))
         CsubM(kxB1+lxrow,kxB1+lxcol,k)=CsubM(kxB1+lxrow,kxB1+lxcol,k)
     &    -cna*FFP(RGROTCM(i,msb))
       ELSEIF (ABS(RM(MROW,2)).GT.0.1.AND.ABS(RM(MSA,2)).LT.0.1) THEN
         BsubM(kxB1+lxrow,kxB1+lxcol,k)=BsubM(kxB1+lxrow,kxB1+lxcol,k)
     $    -cna*FF(RGROTC(i,msb)*T(i),RGROTCM(i-1,msb)*TM(i-1),
     $            RGROTCM(i,msb)*TM(i))
         AsubM(kxB1+lxrow,kxB1+lxcol,k)=AsubM(kxB1+lxrow,kxB1+lxcol,k)
     &    -cna*FFM(RGROTCM(i-1,msb)*TM(i-1))
         CsubM(kxB1+lxrow,kxB1+lxcol,k)=CsubM(kxB1+lxrow,kxB1+lxcol,k)
     &    -cna*FFP(RGROTCM(i,msb)*TM(i))
       ELSEIF (ABS(RM(MROW,2)).LT.0.1.AND.ABS(RM(MSA,2)).GT.0.1) THEN
         BsubM(kxB1+lxrow,kxB1+lxcol,k)=BsubM(kxB1+lxrow,kxB1+lxcol,k)
     $    -cna*FF(RGROTC(i,msb)/T(i),RGROTCM(i-1,msb)/TM(i-1),
     $            RGROTCM(i,msb)/TM(i))
         AsubM(kxB1+lxrow,kxB1+lxcol,k)=AsubM(kxB1+lxrow,kxB1+lxcol,k)
     &    -cna*FFM(RGROTCM(i-1,msb)/TM(i-1))
         CsubM(kxB1+lxrow,kxB1+lxcol,k)=CsubM(kxB1+lxrow,kxB1+lxcol,k)
     &    -cna*FFP(RGROTCM(i,msb)/TM(i))
       ENDIF
       ENDIF

C     GLX------<      
      IF (NPROFRP.GT.0) THEN
      IF ((ABS(RM(MROW,2)).GT.0.1.AND.ABS(RM(MSA,2)).GT.0.1).OR.
     $    (ABS(RM(MROW,2)).LT.0.1.AND.ABS(RM(MSA,2)).LT.0.1)) THEN
         BsubM(kxB1+lxrow,kxB1+lxcol,k)=BsubM(kxB1+lxrow,kxB1+lxcol,k)
     $     - cmrow*FF(RGQ1RQ1P(i,msb),RGQ1RQ1PM(i-1,msb),
     $                                RGQ1RQ1PM(i,msb))
     $     - cna*FF(RGQ1RQ2P(i,msb),RGQ1RQ2PM(i-1,msb),RGQ1RQ2PM(i,msb))
         AsubM(kxB1+lxrow,kxB1+lxcol,k)=AsubM(kxB1+lxrow,kxB1+lxcol,k)
     &     - cmrow*FFM(RGQ1RQ1PM(i-1,msb))
     &     - cna*FFM(RGQ1RQ2PM(i-1,msb))
         CsubM(kxB1+lxrow,kxB1+lxcol,k)=CsubM(kxB1+lxrow,kxB1+lxcol,k)
     &     - cmrow*FFP(RGQ1RQ1PM(i,msb))
     &     - cna*FFP(RGQ1RQ2PM(i,msb))
       ELSEIF (ABS(RM(MROW,2)).GT.0.1.AND.ABS(RM(MSA,2)).LT.0.1) THEN
         BsubM(kxB1+lxrow,kxB1+lxcol,k)=BsubM(kxB1+lxrow,kxB1+lxcol,k)
     $     - cmrow*FF(RGQ1RQ1P(i,msb)*T(i),RGQ1RQ1PM(i-1,msb)*TM(i-1),
     $                                     RGQ1RQ1PM(i,msb)*TM(i))
     $     - cna*FF(RGQ1RQ2P(i,msb)*T(i),RGQ1RQ2PM(i-1,msb)*TM(i-1),
     $                                   RGQ1RQ2PM(i,msb)*TM(i))
         AsubM(kxB1+lxrow,kxB1+lxcol,k)=AsubM(kxB1+lxrow,kxB1+lxcol,k)
     &     - cmrow*FFM(RGQ1RQ1PM(i-1,msb)*TM(i-1))
     &     - cna*FFM(RGQ1RQ2PM(i-1,msb)*TM(i-1))
         CsubM(kxB1+lxrow,kxB1+lxcol,k)=CsubM(kxB1+lxrow,kxB1+lxcol,k)
     &     - cmrow*FFP(RGQ1RQ1PM(i,msb)*TM(i))
     &     - cna*FFP(RGQ1RQ2PM(i,msb)*TM(i))
       ELSEIF (ABS(RM(MROW,2)).LT.0.1.AND.ABS(RM(MSA,2)).GT.0.1) THEN
         BsubM(kxB1+lxrow,kxB1+lxcol,k)=BsubM(kxB1+lxrow,kxB1+lxcol,k)
     $     - cmrow*FF(RGQ1RQ1P(i,msb)/T(i),RGQ1RQ1PM(i-1,msb)/TM(i-1),
     $                                     RGQ1RQ1PM(i,msb)/TM(i))
     $     - cna*FF(RGQ1RQ2P(i,msb)/T(i),RGQ1RQ2PM(i-1,msb)/TM(i-1),
     $                                   RGQ1RQ2PM(i,msb)/TM(i))
         AsubM(kxB1+lxrow,kxB1+lxcol,k)=AsubM(kxB1+lxrow,kxB1+lxcol,k)
     &     - cmrow*FFM(RGQ1RQ1PM(i-1,msb)/TM(i-1))
     &     - cna*FFM(RGQ1RQ2PM(i-1,msb)/TM(i-1))
         CsubM(kxB1+lxrow,kxB1+lxcol,k)=CsubM(kxB1+lxrow,kxB1+lxcol,k)
     &     - cmrow*FFP(RGQ1RQ1PM(i,msb)/TM(i))
     &     - cna*FFP(RGQ1RQ2PM(i,msb)/TM(i))
       ENDIF
       ENDIF
C    --------->

C     DIAMAGNETIC TERM
      IF (IDIAMB.EQ.1.OR.IDIAMB.EQ.3) THEN
         CTMP4 = DPSIDS(I)*FDIAMB
         CTMP5 = DPSIDSM(I-1)*FDIAMB
         CTMP6 = DPSIDSM(I)*FDIAMB
         IF (ABS(RM(MROW,2)).GT.0.1) THEN
            BsubM(kxB1+lxrow, kxPD+lxcol,k)=
     $       -cna*FF(DG22L(i,msb)*CTMP4,
     $        DG22LM(i-1,msb)*CTMP5,DG22LM(i,msb)*CTMP6)
            AsubM(kxB1+lxrow, kxPD+lxcol,k)= 
     &       -cna*FFM(DG22LM(i-1,msb)*CTMP5)
            CsubM(kxB1+lxrow, kxPD+lxcol,k)= 
     &       -cna*FFP(DG22LM(i,msb)*CTMP6)
         ELSEIF (ABS(RM(MROW,2)).LT.0.1) THEN
            BsubM(kxB1+lxrow, kxPD+lxcol,k)=
     $       -cna*FF(DG22L(i,msb)*CTMP4/T(i),DG22LM(i-1,msb)
     $        *CTMP5/TM(i-1),DG22LM(i,msb)*CTMP6/TM(i))
            AsubM(kxB1+lxrow, kxPD+lxcol,k)= 
     &       -cna*FFM(DG22LM(i-1,msb)*CTMP5/TM(i-1))
            CsubM(kxB1+lxrow, kxPD+lxcol,k)= 
     &       -cna*FFP(DG22LM(i,msb)*CTMP6/TM(i))
         ENDIF

         IF (KXPED.GT.0) THEN
         IF (ABS(RM(MROW,2)).GT.0.1) THEN
            BsubM(kxB1+lxrow, kxPED+lxcol,k)=
     $       -cna*FF(DG22L(i,msb)*CTMP4,
     $        DG22LM(i-1,msb)*CTMP5,DG22LM(i,msb)*CTMP6)
            AsubM(kxB1+lxrow, kxPED+lxcol,k)= 
     &       -cna*FFM(DG22LM(i-1,msb)*CTMP5)
            CsubM(kxB1+lxrow, kxPED+lxcol,k)= 
     &       -cna*FFP(DG22LM(i,msb)*CTMP6)
         ELSEIF (ABS(RM(MROW,2)).LT.0.1) THEN
            BsubM(kxB1+lxrow, kxPED+lxcol,k)=
     $       -cna*FF(DG22L(i,msb)*CTMP4/T(i),DG22LM(i-1,msb)
     $        *CTMP5/TM(i-1),DG22LM(i,msb)*CTMP6/TM(i))
            AsubM(kxB1+lxrow, kxPED+lxcol,k)= 
     &       -cna*FFM(DG22LM(i-1,msb)*CTMP5/TM(i-1))
            CsubM(kxB1+lxrow, kxPED+lxcol,k)= 
     &       -cna*FFP(DG22LM(i,msb)*CTMP6/TM(i))
         ENDIF
         ENDIF
      ENDIF

      IF (IDIAMB.EQ.2.OR.IDIAMB.EQ.3) THEN
         CTMP1 = OMEGASE(I)*DPSIDS(I)*FDIAMB            
         CTMP2 = OMEGASEM(I-1)*DPSIDSM(I-1)*FDIAMB            
         CTMP3 = OMEGASEM(I)*DPSIDSM(I)*FDIAMB            
         CTMP4 = CTMP1*DPSIDS(I)*FDIAMB
         CTMP5 = CTMP2*DPSIDSM(I-1)*FDIAMB
         CTMP6 = CTMP3*DPSIDSM(I)*FDIAMB
         IF (ABS(RM(MROW,2)).GT.0.1.AND.ABS(RM(MSA,2)).GT.0.1) THEN
            BsubM(kxB1+lxrow, kxB1+lxcol,k)=
     $       BsubM(kxB1+lxrow, kxB1+lxcol,k)
     $       +cmrow*FF(DJB2(i,msb)*T(i)*CTMP1,
     $        DJB2M(i-1,msb)*TM(i-1)*CTMP2,DJB2M(i,msb)*TM(i)*CTMP3)
     $       -cna*FF(DG22J2B2(i,msb)*CTMP4,
     $        DG22J2B2M(i-1,msb)*CTMP5,DG22J2B2M(i,msb)*CTMP6)
            AsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     $       AsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       +cmrow*FFM(DJB2M(i-1,msb)*TM(i-1)*CTMP2)
     &       -cna*FFM(DG22J2B2M(i-1,msb)*CTMP5)
            CsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     &       CsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       +cmrow*FFP(DJB2M(i,msb)*TM(i)*CTMP3)
     &       -cna*FFP(DG22J2B2M(i,msb)*CTMP6)
         ELSEIF (ABS(RM(MROW,2)).GT.0.1.AND.ABS(RM(MSA,2)).LT.0.1) THEN
            BsubM(kxB1+lxrow, kxB1+lxcol,k)=
     $       BsubM(kxB1+lxrow, kxB1+lxcol,k)
     $       +cmrow*FF(DJB2(i,msb)*T(i)**2*CTMP1,DJB2M(i-1,msb)
     $        *TM(i-1)**2*CTMP2,DJB2M(i,msb)*TM(i)**2*CTMP3)
     $       -cna*FF(DG22J2B2(i,msb)*CTMP4*T(i),DG22J2B2M(i-1,msb)
     $        *CTMP5*TM(i-1),DG22J2B2M(i,msb)*CTMP6*TM(i))
            AsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     $       AsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       +cmrow*FFM(DJB2M(i-1,msb)*TM(i-1)**2*CTMP2)
     &       -cna*FFM(DG22J2B2M(i-1,msb)*CTMP5*TM(i-1))
            CsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     &       CsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       +cmrow*FFP(DJB2M(i,msb)*TM(i)**2*CTMP3)
     &       -cna*FFP(DG22J2B2M(i,msb)*CTMP6*TM(i))
         ELSEIF (ABS(RM(MROW,2)).LT.0.1.AND.ABS(RM(MSA,2)).GT.0.1) THEN
            BsubM(kxB1+lxrow, kxB1+lxcol,k)=
     $       BsubM(kxB1+lxrow, kxB1+lxcol,k)
     $       -cna*FF(DG22J2B2(i,msb)*CTMP4/T(i),DG22J2B2M(i-1,msb)
     $        *CTMP5/TM(i-1),DG22J2B2M(i,msb)*CTMP6/TM(i))
            AsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     $       AsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       -cna*FFM(DG22J2B2M(i-1,msb)*CTMP5/TM(i-1))
            CsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     &       CsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       -cna*FFP(DG22J2B2M(i,msb)*CTMP6/TM(i))
         ELSEIF (ABS(RM(MROW,2)).LT.0.1.AND.ABS(RM(MSA,2)).LT.0.1) THEN
            BsubM(kxB1+lxrow, kxB1+lxcol,k)=
     $       BsubM(kxB1+lxrow, kxB1+lxcol,k)
     $       -cna*FF(DG22J2B2(i,msb)*CTMP4,
     $        DG22J2B2M(i-1,msb)*CTMP5,DG22J2B2M(i,msb)*CTMP6)
            AsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     $       AsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       -cna*FFM(DG22J2B2M(i-1,msb)*CTMP5)
            CsubM(kxB1+lxrow, kxB1+lxcol,k)= 
     &       CsubM(kxB1+lxrow, kxB1+lxcol,k) 
     &       -cna*FFP(DG22J2B2M(i,msb)*CTMP6)
         ENDIF
      ENDIF

      IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10) THEN
      IF (ABS(RM(MROW,2)).GT.0.1) THEN
         HSUBM(KXB1+LXROW,KYB2+LYCOL,I)=
     &   HSUBM(KXB1+LXROW,KYB2+LYCOL,I) 
     &     +CMROW*FGM(FHATV(I,MSB,1),FHATVM(I-1,MSB,1))
         ESUBM(KXB1+LXROW,KYB2+LYCOL,I)=
     &   ESUBM(KXB1+LXROW,KYB2+LYCOL,I) 
     &     +CMROW*FGP(FHATV(I,MSB,1),FHATVM(I,MSB,1))

      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         BSUBM(KXB1+LXROW,KXB1+LXCOL,I)=
     &   BSUBM(KXB1+LXROW,KXB1+LXCOL,I) 
     &      -CMROW*FF(FHATV(I,MSB,2),FHATVM(I-1,MSB,2),FHATVM(I,MSB,2))
         ASUBM(KXB1+LXROW,KXB1+LXCOL,I)=
     &   ASUBM(KXB1+LXROW,KXB1+LXCOL,I)-CMROW*FFM(FHATVM(I-1,MSB,2))
         CSUBM(KXB1+LXROW,KXB1+LXCOL,I)=
     &   CSUBM(KXB1+LXROW,KXB1+LXCOL,I)-CMROW*FFP(FHATVM(I,MSB,2))
      ELSE
         BSUBM(KXB1+LXROW,KXB1+LXCOL,I)=
     &   BSUBM(KXB1+LXROW,KXB1+LXCOL,I)-CMROW*FF(FHATV(I,MSB,2)*T(I),
     &      FHATVM(I-1,MSB,2)*TM(I-1),FHATVM(I,MSB,2)*TM(I))
         ASUBM(KXB1+LXROW,KXB1+LXCOL,I)=
     &   ASUBM(KXB1+LXROW,KXB1+LXCOL,I) 
     &      -CMROW*FFM(FHATVM(I-1,MSB,2)*TM(I-1))
         CSUBM(KXB1+LXROW,KXB1+LXCOL,I)=
     &   CSUBM(KXB1+LXROW,KXB1+LXCOL,I) 
     &      -CMROW*FFP(FHATVM(I,MSB,2)*TM(I))
      ENDIF
      ENDIF
      ENDIF

 40   CONTINUE
C
C
C-----------------------------------------------------------------------
C.. FIFTH  EQUATION  CONTRAVARIANT-2-COMP. OF FARADAY
C..                              (KYB2 = 3, DEFINED ON HALF    MESH)
C-----------------------------------------------------------------------
C
      IF (IPDIVB.NE.2.OR.(IPDIVB.EQ.2.AND.ABS(RM(MROW,2)).LT.0.1)) THEN
      DO 50 I=1,NR
      K = I
C
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      Z1M  = (CS(I  )/CSM(I))**IEXJ1
      Z1P  = (CS(I+1)/CSM(I))**IEXJ1

C
      INCLUDE 'tophat.inc'
C
C-----------------------------------------------------------------------
C
         DsubM(kyB2+lyrow, kyV2+lycol,k)= 
     $cna*GG(JACOBM(i,msb),JACOBI(i,msb)*zv2m,JACOBI(i+1,msb)*zv2p)

         DsubM(kyB2+lyrow, kyJ1+lycol,k)=
     $   -cna*GG(resisM(i)*DG11LM(i,msb),resist(i)*DG11L(i,msb)*z1m,
     $     resist(i+1)*DG11L(i+1,msb)*z1p)

         FsubM(kyB2+lyrow, kxJ2U+lxcol,k)=
     $    -cna*GF(resist(i)*DG12L(i,msb), resisM(i)*DG12LM(i,msb))
         GsubM(kyB2+lyrow, kxJ2U+lxcol,k)=
     $  -cna*GF(resist(i+1)*DG12L(i+1,msb), resisM(i)*DG12LM(i,msb))

         FsubM(kyB2+lyrow, kxJ3+lxcol,k)= 
     $       -znorm*resist(i)*DG33L(i,msb) 
         GsubM(kyB2+lyrow, kxJ3+lxcol,k)=
     $        znorm*resist(i+1)*DG33L(i+1,msb) 

         IF (KJRER.EQ.6) THEN
         CTMP1 = JRE_EQFRAC*JRE_EXB
         DsubM(kyB2+lyrow, KYJRE1+lycol,k)=CNA*CTMP1
     $     *GG(resisM(i)*DG11LM(i,msb),resist(i)*Z1M*DG11L(i,msb),
     $     resist(i+1)*Z1P*DG11L(i+1,msb))

         FsubM(kyB2+lyrow, KXJRE2+lxcol,k)=CNA*CTMP1
     $     *GF(resist(i)*DG12L(i,msb), resisM(i)*DG12LM(i,msb))
         GsubM(kyB2+lyrow, KXJRE2+lxcol,k)=CNA*CTMP1
     $     *GF(resist(i+1)*DG12L(i+1,msb), resisM(i)*DG12LM(i,msb))

         FsubM(kyB2+lyrow, KXJRE3+lxcol,k)= 
     $        znorm*resist(i)*CTMP1*DG33L(i,msb) 
         GsubM(kyB2+lyrow, KXJRE3+lxcol,k)=
     $       -znorm*resist(i+1)*CTMP1*DG33L(i+1,msb) 
         ENDIF

         IF (KJRER.EQ.1.OR.KJRER.EQ.3.OR.KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         FSUBM(KYB2+LYROW,KXJRE+LXCOL,K)=FSUBM(KYB2+LYROW,KXJRE+LXCOL,K)
     &      + ZNORM*RESIST(I)*T(I)*FRE1B(I,MSB)*JRE_EQFRAC
         GSUBM(KYB2+LYROW,KXJRE+LXCOL,K)=GSUBM(KYB2+LYROW,KXJRE+LXCOL,K)
     &      - ZNORM*RESIST(I+1)*T(I+1)*FRE1B(I+1,MSB)*JRE_EQFRAC
         ENDIF

         IF (KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         FSUBM(KYB2+LYROW,KXJRE+LXCOL,K)=FSUBM(KYB2+LYROW,KXJRE+LXCOL,K)
     &    +CNA*GF(RESIST(I)*DPSIDS(I)*FREG12JB(I,MSB),
     &            RESISM(I)*DPSIDSM(I)*FREG12JBM(I,MSB))*JRE_EQFRAC
         GSUBM(KYB2+LYROW,KXJRE+LXCOL,K)=GSUBM(KYB2+LYROW,KXJRE+LXCOL,K)
     &    +CNA*GF(RESIST(I+1)*DPSIDS(I+1)*FREG12JB(I+1,MSB), 
     &            RESISM(I)*DPSIDSM(I)*FREG12JBM(I,MSB))*JRE_EQFRAC

         IF (KXB3L.GT.0) THEN
         FSUBM(KYB2+LYROW,KXB3L+LXCOL,K)=FSUBM(KYB2+LYROW,KXB3L+LXCOL,K)
     &      + ZNORM*JRE_EQFRAC*RESIST(I)*FREJPB(I,MSB)
         GSUBM(KYB2+LYROW,KXB3L+LXCOL,K)=GSUBM(KYB2+LYROW,KXB3L+LXCOL,K)
     &      - ZNORM*JRE_EQFRAC*RESIST(I+1)*FREJPB(I+1,MSB)
         ENDIF

         IF (ABS(RM(MSA,2)).GT.0.1) THEN
         FSUBM(KYB2+LYROW,KXB1+LXCOL,K)=FSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &    +CNA*JRE_EQFRAC*GF(RESIST(I)*FREJPBG11(I,MSB),
     &                       RESISM(I)*FREJPBG11M(I,MSB))
         GSUBM(KYB2+LYROW,KXB1+LXCOL,K)=GSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &    +CNA*JRE_EQFRAC*GF(RESIST(I+1)*FREJPBG11(I+1,MSB),
     &                       RESISM(I)*FREJPBG11M(I,MSB))
         ELSE
         FSUBM(KYB2+LYROW,KXB1+LXCOL,K)=FSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &    +CNA*JRE_EQFRAC*GF(RESIST(I)*T(I)*FREJPBG11(I,MSB),
     &                       RESISM(I)*TM(I)*FREJPBG11M(I,MSB))
         GSUBM(KYB2+LYROW,KXB1+LXCOL,K)=GSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &    +CNA*JRE_EQFRAC*GF(RESIST(I+1)*T(I+1)*FREJPBG11(I+1,MSB),
     &                       RESISM(I)*TM(I)*FREJPBG11M(I,MSB))
         ENDIF

         DSUBM(KYB2+LYROW,KYB2+LYCOL,k)=DSUBM(KYB2+LYROW,KYB2+LYCOL,K) 
     &     +CNA*JRE_EQFRAC*GG(RESISM(I)*FREJPBG12M(I,MSB),
     &                        RESIST(I)*FREJPBG12(I,MSB),
     &                        RESIST(I+1)*FREJPBG12(I+1,MSB))
         ENDIF

         IF (NPROFRC.GT.0) 
     $   DsubM(kyB2+lyrow,kyB2+lycol,k)=DsubM(kyB2+lyrow,kyB2+lycol,k) 
     $     -cna*GG(RGROTCM(i,msb),RGROTC(i,msb),RGROTC(i+1,msb))

C     DIAMAGNETIC TERM
      IF (IDIAMB.EQ.1.OR.IDIAMB.EQ.3) THEN
         FsubM(kyB2+lyrow, kxPD+lxcol,k)= 
     $    +cna*GF(DG12L(i,msb)*DPSIDS(i),DG12LM(i,msb)*DPSIDSM(i))
     $     *FDIAMB
         GsubM(kyB2+lyrow, kxPD+lxcol,k)=
     $    +cna*GF(DG12L(i+1,msb)*DPSIDS(i+1),DG12LM(i,msb)*DPSIDSM(i))
     $     *FDIAMB

         IF (KXPED.GT.0) THEN
         FsubM(kyB2+lyrow, kxPED+lxcol,k)= 
     $    +cna*GF(DG12L(i,msb)*DPSIDS(i),DG12LM(i,msb)*DPSIDSM(i))
     $     *FDIAMB
         GsubM(kyB2+lyrow, kxPED+lxcol,k)=
     $    +cna*GF(DG12L(i+1,msb)*DPSIDS(i+1),DG12LM(i,msb)*DPSIDSM(i))
     $     *FDIAMB
         ENDIF
      ENDIF

      IF (IDIAMB.EQ.2.OR.IDIAMB.EQ.3) THEN
         CTMP1 = OMEGASE(I)*DPSIDS(I)*FDIAMB            
         CTMP2 = OMEGASE(I+1)*DPSIDS(I+1)*FDIAMB            
         CTMP3 = OMEGASEM(I)*DPSIDSM(I)*FDIAMB            
         CTMP4 = CTMP1*DPSIDS(I)*FDIAMB
         CTMP5 = CTMP2*DPSIDS(I+1)*FDIAMB
         CTMP6 = CTMP3*DPSIDSM(I)*FDIAMB
         IF (ABS(RM(MSA,2)).GT.0.1) THEN
            FsubM(kyB2+lyrow, kxB1+lxcol,k)= 
     $       FsubM(kyB2+lyrow, kxB1+lxcol,k) 
     $       +znorm*DJB2(i,msb)*CTMP1*TM(i) 
     $       +cna*GF(DG12J2B2(i,msb)*CTMP4,DG12J2B2M(i,msb)*CTMP6)
     $       -GF(DJB2(i,msb)*CTMP4*TP(i),DJB2M(i,msb)*CTMP6*TPM(i))
            GsubM(kyB2+lyrow, kxB1+lxcol,k)=
     $       GsubM(kyB2+lyrow, kxB1+lxcol,k)
     $       -znorm*DJB2(i+1,msb)*CTMP2*TM(i) 
     $       +cna*GF(DG12J2B2(i+1,msb)*CTMP5,DG12J2B2M(i,msb)*CTMP6)
     $       -GF(DJB2(i+1,msb)*CTMP5*TP(i+1),DJB2M(i,msb)*CTMP6*TPM(i))
         ELSE
            FsubM(kyB2+lyrow, kxB1+lxcol,k)= 
     $       FsubM(kyB2+lyrow, kxB1+lxcol,k) 
     $       +znorm*DJB2(i,msb)*CTMP1*TM(i)*T(i) 
     $       +cna*GF(DG12J2B2(i,msb)*CTMP4*T(i),
     $        DG12J2B2M(i,msb)*CTMP6*TM(i))
     $       -GF(DJB2(i,msb)*CTMP4*TP(i)*T(i),
     $        DJB2M(i,msb)*CTMP6*TPM(i)*TM(i))
            GsubM(kyB2+lyrow, kxB1+lxcol,k)=
     $       GsubM(kyB2+lyrow, kxB1+lxcol,k)
     $       -znorm*DJB2(i+1,msb)*CTMP2*TM(i)*T(i+1) 
     $       +cna*GF(DG12J2B2(i+1,msb)*CTMP5*T(i+1),
     $        DG12J2B2M(i,msb)*CTMP6*TM(i))
     $       -GF(DJB2(i+1,msb)*CTMP5*TP(i+1)*T(i+1),
     $        DJB2M(i,msb)*CTMP6*TPM(i)*TM(i))
         ENDIF
      ENDIF
C    GLX------< 
      IF (NPROFRP.GT.0) THEN 
      DsubM(kyB2+lyrow,kyB2+lycol,k)=DsubM(kyB2+lyrow,kyB2+lycol,k) 
     &   - cma*GG(RGQ1RQ1PM(i,msb),RGQ1RQ1P(i,msb),RGQ1RQ1P(i+1,msb))
     &   - cna*GG(RGQ1RQ2PM(i,msb),RGQ1RQ2P(i,msb),RGQ1RQ2P(i+1,msb))

      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB2+LYROW,KXB1+LXCOL,K)=FSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &   + GF(RGQ2RQ1P(I,MSB),RGQ2RQ1PM(I,MSB))
      GSUBM(KYB2+LYROW,KXB1+LXCOL,K)=GSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &   + GF(RGQ2RQ1P(I+1,MSB),RGQ2RQ1PM(I,MSB))
      ELSE
      FSUBM(KYB2+LYROW,KXB1+LXCOL,K)=FSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &   + GF(RGQ2RQ1P(I,MSB)*T(I),RGQ2RQ1PM(I,MSB)*TM(I))
      GSUBM(KYB2+LYROW,KXB1+LXCOL,K)=GSUBM(KYB2+LYROW,KXB1+LXCOL,K)
     &   + GF(RGQ2RQ1P(I+1,MSB)*T(I+1),RGQ2RQ1PM(I,MSB)*TM(I))
      ENDIF
      ENDIF

      IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10) THEN
      IF (KXB2L.GT.0) THEN
         FSUBM(KYB2+LYROW, KXB2L+LXCOL,I)= 
     &   FSUBM(KYB2+LYROW, KXB2L+LXCOL,I)+ZNORM*FHATV(I,MSB,4) 
         GSUBM(KYB2+LYROW, KXB2L+LXCOL,I)= 
     &   GSUBM(KYB2+LYROW, KXB2L+LXCOL,I)-ZNORM*FHATV(I+1,MSB,4) 
      ENDIF
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
         FSUBM(KYB2+LYROW, KXB1+LXCOL,I)= 
     &   FSUBM(KYB2+LYROW, KXB1+LXCOL,I)-ZNORM*FHATV(I,MSB,5) 
         GSUBM(KYB2+LYROW, KXB1+LXCOL,I)= 
     &   GSUBM(KYB2+LYROW, KXB1+LXCOL,I)+ZNORM*FHATV(I+1,MSB,5) 
      ELSE
         FSUBM(KYB2+LYROW, KXB1+LXCOL,I)= 
     &   FSUBM(KYB2+LYROW, KXB1+LXCOL,I)-ZNORM*FHATV(I,MSB,5)*T(I) 
         GSUBM(KYB2+LYROW, KXB1+LXCOL,I)= 
     &   GSUBM(KYB2+LYROW, KXB1+LXCOL,I)+ZNORM*FHATV(I+1,MSB,5)*T(I+1) 
      ENDIF
      ENDIF

 50   CONTINUE
      ENDIF
C
      IF (IPDIVB.NE.1) THEN
      DO 60 I=1,NR
      K = I
C-----------------------------------------------------------------------
C.. SIXTH  EQUATION: CONTRAVARIANT-3-COMP. OF FARADAY
C..                                  (KYB3 = 4, DEFINED ON HALF    MESH)
C-----------------------------------------------------------------------
C
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      Z1M  = (CS(I  )/CSM(I))**IEXJ1
      Z1P  = (CS(I+1)/CSM(I))**IEXJ1

      z3m  = (cs(i  )/csm(i))**iexb3
      z3p  = (cs(i+1)/csm(i))**iexb3
      INCLUDE 'tophat.inc'
C-----------------------------------------------------------------------
C
         IF (KVSQL(8)) THEN
         FsubM(kyB3+lyrow, kxV1+lxcol,k)=  znorm*B3j2(i,msb)*T(I)
         GsubM(kyB3+lyrow, kxV1+lxcol,k)= -znorm*B3j2(i+1,msb)*T(I+1)

         DsubM(kyB3+lyrow, kyV2+lycol,k)=
     &     -cma*GG(JACOBM(i,msb),
     $      JACOBI(i,msb)*zv2m, JACOBI(i+1,msb)*zv2p)
     &     -GG(B3Gv2M_3dc(i,msb),
     $      B3Gv2_3dc(i,msb)*zv2m, B3Gv2_3dc(i+1,msb)*zv2p)
         ENDIF

         DsubM(kyB3+lyrow, kyJ1+lycol,k)= 
     $      cma*GG(resisM(i)*DG11LM(i,msb),resist(i)*
     &      DG11L(i,msb)*z1m,resist(i+1)*DG11L(i+1,msb)*z1p)
     $     +GG(resisM(i)*G11j1M_3dc(i,msb),resist(i)*
     &      G11j1_3dc(i,msb)*z1m,resist(i+1)*G11j1_3dc(i+1,msb)*z1p)

         FsubM(kyB3+lyrow, kxJ2U+lxcol,k)=
     & cma*GF(resist(i)*DG12L(i,msb), resisM(i)*DG12LM(i,msb))    
     &    +GF(resist(i)*G12j2_3dc(i,msb), resisM(i)*G12j2M_3dc(i,msb)) 

         GsubM(kyB3+lyrow, kxJ2U+lxcol,k)=
     &   cma*GF(resist(i+1)*DG12L(i+1,msb), resisM(i)*DG12L(i,msb))
     &  +GF(resist(i+1)*G12j2_3dc(i+1,msb), resisM(i)*G12j2M_3dc(i,msb))

         IF (KJRER.EQ.6) THEN
         CTMP1 = JRE_EQFRAC*JRE_EXB
         DsubM(kyB3+lyrow, KYJRE1+lycol,k)= 
     $     -cma*CTMP1*GG(resisM(i)*DG11LM(i,msb),resist(i)*
     &      Z1M*DG11L(i,msb),resist(i+1)*Z1P*DG11L(i+1,msb))
     $     -CTMP1*GG(resisM(i)*G11j1M_3dc(i,msb),resist(i)*
     &      Z1M*G11j1_3dc(i,msb),resist(i+1)*Z1P*G11j1_3dc(i+1,msb))

         FsubM(kyB3+lyrow, KXJRE2+lxcol,k)=
     &     -cma*CTMP1*GF(resist(i)*DG12L(i,msb),
     &                        resisM(i)*DG12LM(i,msb))    
     &     -CTMP1*GF(resist(i)*G12j2_3dc(i,msb), 
     &                   resisM(i)*G12j2M_3dc(i,msb)) 

         GsubM(kyB3+lyrow, KXJRE2+lxcol,k)=
     &     -cma*CTMP1*GF(resist(i+1)*DG12L(i+1,msb), 
     &                        resisM(i)*DG12L(i,msb))
     &     -CTMP1*GF(resist(i+1)*G12j2_3dc(i+1,msb), 
     &                    resisM(i)*G12j2M_3dc(i,msb))
         ENDIF
      
         IF (KJRER.EQ.1.OR.KJRER.EQ.3.OR.KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         FSUBM(KYB3+LYROW,KXJRE+LXCOL,K)=FSUBM(KYB3+LYROW,KXJRE+LXCOL,K)
     &      - ZNORM*RESIST(I)*DPSIDS(I)*FREG22JB(I,MSB)*JRE_EQFRAC
     &      - CMROW*GF(RESIST(I)*DPSIDS(I)*FREG12JB(I,MSB),
     &                 RESISM(I)*DPSIDSM(I)*FREG12JBM(I,MSB))*JRE_EQFRAC
         GSUBM(KYB3+LYROW,KXJRE+LXCOL,K)=GSUBM(KYB3+LYROW,KXJRE+LXCOL,K)
     &      + ZNORM*RESIST(I+1)*DPSIDS(I+1)*FREG22JB(I+1,MSB)*JRE_EQFRAC
     &      - CMROW*GF(RESIST(I+1)*DPSIDS(I+1)*FREG12JB(I+1,MSB),
     &                 RESISM(I)*DPSIDSM(I)*FREG12JBM(I,MSB))*JRE_EQFRAC
         ENDIF

         IF (KJRER.EQ.5.OR.KJRER.EQ.6) THEN
         IF (KXB2L.GT.0) THEN
         FSUBM(KYB3+LYROW,KXB2L+LXCOL,K)=FSUBM(KYB3+LYROW,KXB2L+LXCOL,K)
     &      - ZNORM*JRE_EQFRAC*RESIST(I)*FREJPB(I,MSB)
         GSUBM(KYB3+LYROW,KXB2L+LXCOL,K)=GSUBM(KYB3+LYROW,KXB2L+LXCOL,K)
     &      + ZNORM*JRE_EQFRAC*RESIST(I+1)*FREJPB(I+1,MSB)
         ENDIF

         IF (ABS(RM(MSA,2)).GT.0.1) THEN
         FSUBM(KYB3+LYROW,KXB1+LXCOL,K)=FSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &    -CMROW*JRE_EQFRAC*GF(RESIST(I)*FREJPBG11(I,MSB),
     &                         RESISM(I)*FREJPBG11M(I,MSB))
         GSUBM(KYB3+LYROW,KXB1+LXCOL,K)=GSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &    -CMROW*JRE_EQFRAC*GF(RESIST(I+1)*FREJPBG11(I+1,MSB),
     &                         RESISM(I)*FREJPBG11M(I,MSB))
         ELSE
         FSUBM(KYB3+LYROW,KXB1+LXCOL,K)=FSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &    -CMROW*JRE_EQFRAC*GF(RESIST(I)*T(I)*FREJPBG11(I,MSB),
     &                         RESISM(I)*TM(I)*FREJPBG11M(I,MSB))
         GSUBM(KYB3+LYROW,KXB1+LXCOL,K)=GSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &    -CMROW*JRE_EQFRAC*GF(RESIST(I+1)*T(I+1)*FREJPBG11(I+1,MSB),
     &                         RESISM(I)*TM(I)*FREJPBG11M(I,MSB))
         ENDIF

         DSUBM(KYB3+LYROW,KYB2+LYCOL,k)=DSUBM(KYB3+LYROW,KYB2+LYCOL,K) 
     &     -CMROW*JRE_EQFRAC*GG(RESISM(I)*FREJPBG12M(I,MSB),
     &                          RESIST(I)*FREJPBG12(I,MSB),
     &                          RESIST(I+1)*FREJPBG12(I+1,MSB))
         ENDIF

         IF (NPROFRC.GT.0) 
     $   DsubM(kyB3+lyrow,kyB3+lycol,k)=DsubM(kyB3+lyrow,kyB3+lycol,k) 
     $     -cna*GG(RGROTCM(i,msb),RGROTC(i,msb)*z3m,RGROTC(i+1,msb)*z3p)

C     DIAMAGNETIC TERM
      IF (IDIAMB.EQ.1.OR.IDIAMB.EQ.3) THEN
         CTMP4 = DPSIDS(I)*FDIAMB
         CTMP5 = DPSIDS(I+1)*FDIAMB
         CTMP6 = DPSIDSM(I)*FDIAMB
         FsubM(kyB3+lyrow, kxPD+lxcol,k)= 
     $    -znorm*DG22LM(i,msb)*CTMP6 
     $    -cma*GF(DG12L(i,msb)*CTMP4,DG12LM(i,msb)*CTMP6)
     $    +GF(J3B2(i,msb),J3B2M(i,msb))
         GsubM(kyB3+lyrow, kxPD+lxcol,k)=
     $    +znorm*DG22LM(i,msb)*CTMP6 
     $    -cma*GF(DG12L(i+1,msb)*CTMP5,DG12LM(i,msb)*CTMP6)
     $    +GF(J3B2(i+1,msb),J3B2M(i,msb))

         IF (KXPED.GT.0) THEN
         FsubM(kyB3+lyrow, kxPED+lxcol,k)= 
     $    -znorm*DG22LM(i,msb)*CTMP6 
     $    -cma*GF(DG12L(i,msb)*CTMP4,DG12LM(i,msb)*CTMP6)
     $    +GF(J3B2(i,msb),J3B2M(i,msb))
         GsubM(kyB3+lyrow, kxPED+lxcol,k)=
     $    +znorm*DG22LM(i,msb)*CTMP6 
     $    -cma*GF(DG12L(i+1,msb)*CTMP5,DG12LM(i,msb)*CTMP6)
     $    +GF(J3B2(i+1,msb),J3B2M(i,msb))
         ENDIF
      ENDIF

      IF (IDIAMB.EQ.2.OR.IDIAMB.EQ.3) THEN
         CTMP1 = OMEGASE(I)*DPSIDS(I)*FDIAMB            
         CTMP2 = OMEGASE(I+1)*DPSIDS(I+1)*FDIAMB            
         CTMP3 = OMEGASEM(I)*DPSIDSM(I)*FDIAMB            
         CTMP4 = CTMP1*DPSIDS(I)*FDIAMB
         CTMP5 = CTMP2*DPSIDS(I+1)*FDIAMB
         CTMP6 = CTMP3*DPSIDSM(I)*FDIAMB
         IF (ABS(RM(MSA,2)).GT.0.1) THEN
            FsubM(kyB3+lyrow, kxB1+lxcol,k)= 
     $       FsubM(kyB3+lyrow, kxB1+lxcol,k) 
     $       -znorm*DG22J2B2M(i,msb)*CTMP6 
     $       -cma*GF(DG12J2B2(i,msb)*CTMP4,DG12J2B2M(i,msb)*CTMP6)
     $       +GF(DGJBA(i,msb)*CTMP1,DGJBAM(i,msb)*CTMP3)
            GsubM(kyB3+lyrow, kxB1+lxcol,k)=
     $       GsubM(kyB3+lyrow, kxB1+lxcol,k)
     $       +znorm*DG22J2B2M(i,msb)*CTMP6 
     $       -cma*GF(DG12J2B2(i+1,msb)*CTMP5,DG12J2B2M(i,msb)*CTMP6)
     $       +GF(DGJBA(i+1,msb)*CTMP2,DGJBAM(i,msb)*CTMP3)
         ELSE
            FsubM(kyB3+lyrow, kxB1+lxcol,k)= 
     $       FsubM(kyB3+lyrow, kxB1+lxcol,k) 
     $       -znorm*DG22J2B2M(i,msb)*CTMP6*T(i) 
     $       +GF(DGJBA(i,msb)*CTMP1*T(i),DGJBAM(i,msb)*CTMP3*TM(i))
            GsubM(kyB3+lyrow, kxB1+lxcol,k)=
     $       GsubM(kyB3+lyrow, kxB1+lxcol,k)
     $       +znorm*DG22J2B2M(i,msb)*CTMP6*T(i+1) 
     $       +GF(DGJBA(i+1,msb)*CTMP2*T(i+1),DGJBAM(i,msb)*CTMP3*TM(i))
         ENDIF
      ENDIF

C     LLI
      IF (ABS(RM(MSA,2)).GT.0.1.AND.NPROFRC.GT.0) THEN 
      FSUBM(KYB3+LYROW,KXB1+LXCOL,K)=FSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     $       +GF(RGDRCDSC(I,MSB),RGDRCDSCM(I,MSB))
      GSUBM(KYB3+LYROW,KXB1+LXCOL,K)=GSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &       +GF(RGDRCDSC(I+1,MSB),RGDRCDSCM(I,MSB))
      ELSEIF (NPROFRC.GT.0) THEN
      FSUBM(KYB3+LYROW,KXB1+LXCOL,K)=FSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     $       +GF(RGDRCDSC(I,MSB)*T(I),RGDRCDSCM(I,MSB)*TM(I))
      GSUBM(KYB3+LYROW,KXB1+LXCOL,K)=GSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &       +GF(RGDRCDSC(I+1,MSB)*T(I+1),RGDRCDSCM(I,MSB)*TM(I))
      ENDIF

C     LLI
      IF (NPROFRC.GT.0) 
     $   DsubM(kyB3+lyrow,kyB2+lycol,k)=DsubM(kyB3+lyrow,kyB2+lycol,k) 
     $     +GG(RGDRCDCCM(i,msb),RGDRCDCC(i,msb),RGDRCDCC(i+1,msb))

C     GLX------<
      IF (NPROFRP.GT.0) THEN 
      DsubM(kyB3+lyrow,kyB3+lycol,k)=DsubM(kyB3+lyrow,kyB3+lycol,k) 
     $   - CMROW*GG(RGQ1RQ1PM(i,msb),RGQ1RQ1P(i,msb)*Z3M,
     &                               RGQ1RQ1P(i+1,msb)*Z3P)
     $   - CNA*GG(RGQ1RQ2PM(i,msb),RGQ1RQ2P(i,msb)*Z3M,
     &                             RGQ1RQ2P(i+1,msb)*Z3P)

      DsubM(kyB3+lyrow,kyB2+lycol,k)=DsubM(kyB3+lyrow,kyB2+lycol,k) 
     $   + CMB*GG(RGQ1RQ2PM(i,msb),RGQ1RQ2P(i,msb),RGQ1RQ2P(i+1,msb))

      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
      FSUBM(KYB3+LYROW,KXB1+LXCOL,K)=FSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     $   + GF(RGQ3RQ1P(I,MSB),RGQ3RQ1PM(I,MSB))
      GSUBM(KYB3+LYROW,KXB1+LXCOL,K)=GSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &   + GF(RGQ3RQ1P(I+1,MSB),RGQ3RQ1PM(I,MSB))
      ELSE
      FSUBM(KYB3+LYROW,KXB1+LXCOL,K)=FSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     $   + GF(RGQ3RQ1P(I,MSB)*T(I),RGQ3RQ1PM(I,MSB)*TM(I))
      GSUBM(KYB3+LYROW,KXB1+LXCOL,K)=GSUBM(KYB3+LYROW,KXB1+LXCOL,K)
     &   + GF(RGQ3RQ1P(I+1,MSB)*T(I+1),RGQ3RQ1PM(I,MSB)*TM(I))
      ENDIF
      ENDIF

      IF (KVSQLIN.AND.KVSQL(3).AND.NCASE.EQ.10) THEN
      IF (ABS(RM(MSA,2)).GT.0.1) THEN 
         FSUBM(KYB3+LYROW, KXB1+LXCOL,I)= 
     &   FSUBM(KYB3+LYROW, KXB1+LXCOL,I)-ZNORM*FHATV(I,MSB,3) 
         GSUBM(KYB3+LYROW, KXB1+LXCOL,I)= 
     &   GSUBM(KYB3+LYROW, KXB1+LXCOL,I)+ZNORM*FHATV(I+1,MSB,3) 
      ELSE
         FSUBM(KYB3+LYROW, KXB1+LXCOL,I)= 
     &   FSUBM(KYB3+LYROW, KXB1+LXCOL,I)-ZNORM*FHATV(I,MSB,3)*T(I) 
         GSUBM(KYB3+LYROW, KXB1+LXCOL,I)= 
     &   GSUBM(KYB3+LYROW, KXB1+LXCOL,I)+ZNORM*FHATV(I+1,MSB,3)*T(I+1) 
      ENDIF

      IF (KXB3L.GT.0) THEN 
         FSUBM(KYB3+LYROW, KXB3L+LXCOL,I)= 
     &   FSUBM(KYB3+LYROW, KXB3L+LXCOL,I)+ZNORM*FHATV(I,MSB,6) 
         GSUBM(KYB3+LYROW, KXB3L+LXCOL,I)= 
     &   GSUBM(KYB3+LYROW, KXB3L+LXCOL,I)-ZNORM*FHATV(I+1,MSB,6) 
      ENDIF

      DSUBM(KYB3+LYROW,KYB3+LYCOL,I)=
     &DSUBM(KYB3+LYROW,KYB3+LYCOL,I)-CMROW*GG(FHATVM(I,MSB,2),
     &   FHATV(I,MSB,2),FHATV(I+1,MSB,2))
      DSUBM(KYB3+LYROW,KYB2+LYCOL,I)=
     &DSUBM(KYB3+LYROW,KYB2+LYCOL,I)+CMROW*GG(FHATVM(I,MSB,3),
     &   FHATV(I,MSB,3),FHATV(I+1,MSB,3))
      ENDIF

60    CONTINUE
      ENDIF

C
C-----------------------------------------------------------------------
C.. SEVENTH EQUATION: COVARIANT-1-COMP. OF AMPERE'S LAW J = CURL(B)
C..                              (KYJ1= 4, DEFINED ON HALF    MESH)
C-----------------------------------------------------------------------
C
      DO 70 I=1,NR
      K = I
      z3m  = (cs(i  )/csm(i))**iexb3
      z3p  = (cs(i+1)/csm(i))**iexb3
      z1m  = (cs(i  )/csm(i))**iexj1
      z1p  = (cs(i+1)/csm(i))**iexj1

      INCLUDE 'tophat.inc'

         DsubM(kyJ1+lyrow, kyB3+lycol,k)= -cmrow*GG(DG33LM(i,msb),
     $ DG33L(i,msb)*z3m, DG33L(i+1,msb)*z3p) 

      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         FsubM(kyJ1+lyrow, kxB1+lxcol,k)=
     $        cna*GF(DG12L(i,msb), DG12LM(i,msb))    
         GsubM(kyJ1+lyrow, kxB1+lxcol,k)=
     $        cna*GF(DG12L(i+1,msb), DG12LM(i,msb))
      ELSE
         FsubM(kyJ1+lyrow, kxB1+lxcol,k)=
     $        cna*GF(DG12L(i,msb)*T(I), DG12LM(i,msb)*TM(I))    
         GsubM(kyJ1+lyrow, kxB1+lxcol,k)=
     $        cna*GF(DG12L(i+1,msb)*T(I+1), DG12LM(i,msb)*TM(I))
      ENDIF
         
         DsubM(kyJ1+lyrow, kyB2+lycol,k)= cna*GG(DG22LM(i,msb),
     $        DG22L(i,msb), DG22L(i+1,msb))
 70   CONTINUE
C
C-----------------------------------------------------------------------
C.. EIGHTH EQUATION: COVARIANT-2-COMP. OF AMPERE
C..                              (KXJ2U = 3, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
C
      DO 80 I=2,NR
      K = I
  
      INCLUDE 'tent.inc'
C
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         BsubM(kxJ2U+lxrow, kxB1+lxcol,k)=
     $    - cna*FF(DG11L(i,msb), DG11LM(i-1,msb), DG11LM(i,msb))
         AsubM(kxJ2U+lxrow, kxB1+lxcol,k)= -cna*FFM(DG11LM(i-1,msb))
         CsubM(kxJ2U+lxrow, kxB1+lxcol,k)= -cna*FFP(DG11LM(i,msb))
      ELSE
         BsubM(kxJ2U+lxrow, kxB1+lxcol,k)=
     $    - cna*FF(DG11L(i,msb)*T(I), DG11LM(i-1,msb)*TM(I-1), 
     $      DG11LM(i,msb)*TM(I))
         AsubM(kxJ2U+lxrow, kxB1+lxcol,k)= 
     $    - cna*FFM(DG11LM(i-1,msb)*TM(I-1))
         CsubM(kxJ2U+lxrow, kxB1+lxcol,k)= 
     $    - cna*FFP(DG11LM(i,msb)*TM(I))
      ENDIF

         HsubM(kxJ2U+lxrow, kyB2+lycol,k)=
     $     -cna* FGM(DG12L(i,msb), DG12LM(i-1,msb))
         EsubM(kxJ2U+lxrow, kyB2+lycol,k)=
     $     -cna* FGP(DG12L(i,msb), DG12LM(i,msb))

         hc = csh(i-1)
         z3m = (cs(i-1)/csm(i-1))**iexb3
         z3p = (cs(i  )/csm(i-1))**iexb3
         HsubM(kxJ2U+lxrow, kyB3+lycol,k)=
     $      - GG(DG33LM(i-1,msb),DG33L(i-1,msb)*z3m,
     &        DG33L(i,msb)*z3p) /csh(i-1)
         hc = csh(i)
         z3m = (cs(i)  /csm(i))**iexb3
         z3p = (cs(i+1)/csm(i))**iexb3
         EsubM(kxJ2U+lxrow, kyB3+lycol,k)=
     $        GG(DG33LM(i,msb),DG33L(i,msb)*z3m,
     &        DG33L(i+1,msb)*z3p) /csh(i)
 80   CONTINUE

C
C***********************************************************************
C BOUNDARY POINT NRP1
      I = NRP1
      K = I
C
      HP = VCSH(1)
      HM = CSH(NR)
      ZNORM = 2./(HP + HM)
      PTRAP = PTRAPI

      if (NV.GT.2) then
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         BsubM(kxJ2U+lxrow, kxB1+lxcol,k)=
     $    - cna* (FF0(DG11L(i,msb), DG11LM(i-1,msb)) 
     &         +  FF1(VDG11L(1,msb), VDG11LM(1,msb)) )
         AsubM(kxJ2U+lxrow, kxB1+lxcol,k)= -cna*FFM(DG11LM(i-1,msb))
         CsubM(kxJ2U+lxrow, kxB1+lxcol,k)= -cna*FFP(VDG11LM(1,msb))
      ELSE
         BsubM(kxJ2U+lxrow, kxB1+lxcol,k)=
     $    - cna* (FF0(DG11L(i,msb)*T(I), DG11LM(i-1,msb)*TM(I-1)) 
     &         +  FF1(VDG11L(1,msb)*T(NRP1),VDG11LM(1,msb)*T(NRP1)))
         AsubM(kxJ2U+lxrow, kxB1+lxcol,k)= 
     $    - cna*FFM(DG11LM(i-1,msb)*TM(I-1))
         CsubM(kxJ2U+lxrow, kxB1+lxcol,k)= 
     $    - cna*FFP(VDG11LM(1,msb)*T(NRP1))
      ENDIF

         HsubM(kxJ2U+lxrow, kyB2+lycol,k)=
     $     -cna* FGM(DG12L(i,msb), DG12LM(i-1,msb))
         EsubM(kxJ2U+lxrow, kyB2+lycol,k)=
     $     -cna* FGP(VDG12L(1,msb), VDG12LM(1,msb))

      HC = CSH(NR)
      Z3M = (CS(I-1)/CSM(I-1))**IEXB3
      Z3P = (CS(I  )/CSM(I-1))**IEXB3
         HsubM(kxJ2U+lxrow, kyB3+lycol,k)=
     $      - GG(DG33LM(i-1,msb),DG33L(i-1,msb)*z3m,
     &        DG33L(i,msb)*z3p) /hc

      HC = VCSH(1)
         EsubM(kxJ2U+lxrow, kyB3+lycol,k)=
     $        GG(VDG33LM(1,msb),VDG33L(1,msb),
     $        VDG33L(2,msb))/hc
      endif
 81   CONTINUE

C***********************************************************************
C-----------------------------------------------------------------------
C.. NINTH  EQUATION: COVARIANT-3-COMP. OF AMPERE  J = CURL(B)
C..                            (KXJ3 = 4, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
C
C----------------------------------------
      DO 90 I=2,NR
      K = I
C
      INCLUDE 'tent.inc'
C----------------------------------------
C
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         BsubM(kxJ3+lxrow, kxB1+lxcol,k)=
     $      FGM(DG12L(i,msb),DG12LM(i-1,msb)) / csh(i-1)
     $      - FGP(DG12L(i,msb),  DG12LM(i,msb)) /csh(i)
     $      +cmrow*FF(DG11L(i,msb),DG11LM(i-1,msb), DG11LM(i,msb))
         AsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     &       cmrow*FFM(DG11LM(i-1,msb))
     $       +FGM(DG12L(i-1,msb),  DG12LM(i-1,msb)) / csh(i-1)
         CsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     &        cmrow*FFP(DG11LM(i,msb))
     $        -FGP(DG12L(i+1,msb),  DG12LM(i,msb)) / csh(i)
      ELSE
         BsubM(kxJ3+lxrow, kxB1+lxcol,k)=
     $      FGM(DG12L(i,msb)*T(I),DG12LM(i-1,msb)*TM(I-1)) / csh(i-1)
     $      - FGP(DG12L(i,msb)*T(I),  DG12LM(i,msb)*TM(I)) /csh(i)
     $      +cmrow*FF(DG11L(i,msb)*T(I),DG11LM(i-1,msb)*TM(I-1), 
     $      DG11LM(i,msb)*TM(I))
         AsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     &      cmrow*FFM(DG11LM(i-1,msb)*TM(I-1))
     $      +FGM(DG12L(i-1,msb)*T(I-1),DG12LM(i-1,msb)*TM(I-1))/csh(i-1)
         CsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     &        cmrow*FFP(DG11LM(i,msb)*TM(I))
     $        -FGP(DG12L(i+1,msb)*T(I+1),DG12LM(i,msb)*TM(I))/csh(i)
      ENDIF

         hc = csh(i-1)
         HsubM(kxJ3+lxrow, kyB2+lycol,k)=
     $   GG(DG22LM(i-1,msb),DG22L(i-1,msb),DG22L(i,msb)) /csh(i-1)
     $     + cmrow*FGM(DG12L(i,msb),DG12LM(i-1,msb))
         hc = csh(i)
         EsubM(kxJ3+lxrow, kyB2+lycol,k)=
     $     - GG(DG22LM(i,msb),DG22L(i,msb),DG22L(i+1,msb)) /csh(i)
     $     + cmrow*FGP(DG12L(i,msb),  DG12LM(i,msb))

C     EQUATION FOR JRE WITH N.NE.0
      IF (KXJRE.GT.0.AND.(KJRER.EQ.5.OR.KJRER.EQ.6)) THEN
         TMP = C_VA
         BSUBM(KXJRE+LXROW, KXJRE+LXCOL,K)=
     &     BSUBM(KXJRE+LXROW,KXJRE+LXCOL,K)
     &    -TMP*CMA*FF(DPSIDS(I)*FRE1JB(I,MSB),
     &                DPSIDSM(I-1)*FRE1JBM(I-1,MSB), 
     &                DPSIDSM(I)*FRE1JBM(I,MSB)) 
     &    -TMP*FF(DPSIDS(I)*FRE1JCB(I,MSB),
     &            DPSIDSM(I-1)*FRE1JCBM(I-1,MSB), 
     &            DPSIDSM(I)*FRE1JCBM(I,MSB)) 
     &    -TMP*CNA*FF(T(I)*FRE1R2B(I,MSB),
     &                TM(I-1)*FRE1R2BM(I-1,MSB), 
     &                TM(I)*FRE1R2BM(I,MSB)) 
         ASUBM(KXJRE+LXROW, KXJRE+LXCOL,K)=
     &     ASUBM(KXJRE+LXROW,KXJRE+LXCOL,K)
     &    -TMP*CMA*FFM(DPSIDSM(I-1)*FRE1JBM(I-1,MSB))
     &    -TMP*FFM(DPSIDSM(I-1)*FRE1JCBM(I-1,MSB))
     &    -TMP*CNA*FFM(TM(I-1)*FRE1R2BM(I-1,MSB))
         CSUBM(KXJRE+LXROW, KXJRE+LXCOL,K)=
     &     CSUBM(KXJRE+LXROW,KXJRE+LXCOL,K)
     &    -TMP*CMA*FFP(DPSIDSM(I)*FRE1JBM(I,MSB))
     &    -TMP*FFP(DPSIDSM(I)*FRE1JCBM(I,MSB))
     &    -TMP*CNA*FFP(TM(I)*FRE1R2BM(I,MSB))

         TMP = C_VA
         IF (ABS(RM(MSA,2)).GT.0.1) THEN
            BSUBM(KXJRE+LXROW, KXB1+LXCOL,K)=
     &        BSUBM(KXJRE+LXROW,KXB1+LXCOL,K)
     &       -TMP*FF(FREJPBS(I,MSB),
     &               FREJPBSM(I-1,MSB), 
     &               FREJPBSM(I,MSB)) 
            ASUBM(KXJRE+LXROW, KXB1+LXCOL,K)=
     &        ASUBM(KXJRE+LXROW,KXB1+LXCOL,K)
     &       -TMP*FFM(FREJPBSM(I-1,MSB))
            CSUBM(KXJRE+LXROW, KXB1+LXCOL,K)=
     &        CSUBM(KXJRE+LXROW,KXB1+LXCOL,K)
     &       -TMP*FFP(FREJPBSM(I,MSB))
         ELSE
            BSUBM(KXJRE+LXROW, KXB1+LXCOL,K)=
     &        BSUBM(KXJRE+LXROW,KXB1+LXCOL,K)
     &       -TMP*FF(FREJPBS(I,MSB)*T(I),
     &               FREJPBSM(I-1,MSB)*TM(I-1), 
     &               FREJPBSM(I,MSB)*TM(I)) 
            ASUBM(KXJRE+LXROW, KXB1+LXCOL,K)=
     &        ASUBM(KXJRE+LXROW,KXB1+LXCOL,K)
     &       -TMP*FFM(FREJPBSM(I-1,MSB)*TM(I-1))
            CSUBM(KXJRE+LXROW, KXB1+LXCOL,K)=
     &        CSUBM(KXJRE+LXROW,KXB1+LXCOL,K)
     &       -TMP*FFP(FREJPBSM(I,MSB)*TM(I))
         ENDIF

         HSUBM(KXJRE+LXROW, KYB2+LYCOL,K)=
     &     HSUBM(KXJRE+LXROW, KYB2+LYCOL,K) 
     &    -TMP*FGM(FREJPBC(I,MSB),FREJPBCM(I-1,MSB))
         ESUBM(KXJRE+LXROW, KYB2+LYCOL,K)=
     &     ESUBM(KXJRE+LXROW, KYB2+LYCOL,K) 
     &    -TMP*FGP(FREJPBC(I,MSB),FREJPBCM(I,MSB))
      ENDIF

      IF (KXJRE.GT.0.AND.KJRER.EQ.6) THEN
         Z1M  = (CS(I)/CSM(I-1))**IEXJ1
         Z1P  = (CS(I)/CSM(I))**IEXJ1
         TMP  = C_VA*JRE_EXB
         HC = CSH(I-1)
         HSUBM(KXJRE+LXROW, KYJRE1+LYCOL,I)=+TMP*GG(FRE1JM(I,MSB),
     &     FRE1J(I,MSB)*Z1M,FRE1J(I+1,MSB)*Z1P)/HC
         HC = CSH(I)
         ESUBM(KXJRE+LXROW, KYJRE1+LYCOL,I)=-TMP*GG(FRE1JM(I,MSB),
     &     FRE1J(I,MSB)*Z1M,FRE1J(I+1,MSB)*Z1P)/HC

         BSUBM(KXJRE+LXROW,KXJRE2+LXCOL,I)=-CMA*TMP*FF(FRE1J(I,MSB),
     &     FRE1JM(I-1,MSB),FRE1JM(I,MSB))
         ASUBM(KXJRE+LXROW,KXJRE2+LXCOL,I)=-CMA*TMP*FFM(FRE1JM(I-1,MSB))
         CSUBM(KXJRE+LXROW,KXJRE2+LXCOL,I)=-CMA*TMP*FFP(FRE1JM(I,MSB))

         BSUBM(KXJRE+LXROW,KXJRE3+LXCOL,I)=-CNA*TMP*FF(FRE1J(I,MSB),
     &     FRE1JM(I-1,MSB),FRE1JM(I,MSB))
         ASUBM(KXJRE+LXROW,KXJRE3+LXCOL,I)=-CNA*TMP*FFM(FRE1JM(I-1,MSB))
         CSUBM(KXJRE+LXROW,KXJRE3+LXCOL,I)=-CNA*TMP*FFP(FRE1JM(I,MSB))
      ENDIF

 90   CONTINUE

C
C***********************************************************************
C BOUNDARY POINT NRP1
      I = NRP1
      K = I
C
      HP = VCSH(1)
      HM = CSH(NR)
      ZNORM = 2./(HP + HM)
      PTRAP = PTRAPI

c     bounddddddddddd

      if(NV.GT.2) then
      IF (ABS(RM(MSA,2)).GT.0.1) THEN
         BsubM(kxJ3+lxrow, kxB1+lxcol,k)=
     $      FGM(DG12L(i,msb),DG12LM(i-1,msb)) / csh(i-1)
     $      - FGP(VDG12L(1,msb),  VDG12LM(1,msb)) /Vcsh(1)
     $      +cmrow*(  FF0(DG11L(i,msb),DG11LM(i-1,msb)) 
     &      +          FF1(VDG11L(1,msb),VDG11LM(1,msb)) ) 

         AsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     $       cmrow*FFM(DG11LM(i-1,msb)) 
     $       +FGM(DG12L(i-1,msb), DG12LM(i-1,msb)) / csh(i-1)

         CsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     $        cmrow*FFP(VDG11LM(1,msb)) 
     $        -FGP(VDG12L(2,msb),  VDG12LM(1,msb)) / Vcsh(1)
      ELSE
         BsubM(kxJ3+lxrow, kxB1+lxcol,k)=
     $      FGM(DG12L(i,msb)*T(I),DG12LM(i-1,msb)*TM(I-1))/csh(i-1)
     $      - FGP(VDG12L(1,msb)*T(NRP1),
     $      VDG12LM(1,msb)*T(NRP1))/Vcsh(1)
     $      +cmrow*(  FF0(DG11L(i,msb)*T(I),DG11LM(i-1,msb)*TM(I-1)) 
     &      +FF1(VDG11L(1,msb)*T(NRP1),VDG11LM(1,msb)*T(NRP1))) 

         AsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     $       cmrow*FFM(DG11LM(i-1,msb)*TM(I-1)) 
     $       +FGM(DG12L(i-1,msb)*T(I-1), 
     $       DG12LM(i-1,msb)*TM(I-1)) / csh(i-1)

         CsubM(kxJ3+lxrow, kxB1+lxcol,k)= 
     $        cmrow*FFP(VDG11LM(1,msb)*T(NRP1))  
     $        -FGP(VDG12L(2,msb)*T(NRP1),  
     $        VDG12LM(1,msb)*T(NRP1)) / Vcsh(1)
      ENDIF

      HC = CSH(NR)
         HsubM(kxJ3+lxrow, kyB2+lycol,k)=
     $   GG(DG22LM(i-1,msb),DG22L(i-1,msb),DG22L(i,msb)) /hc
     $     + cmrow*FGM(DG12L(i,msb),DG12LM(i-1,msb))

      HC = VCSH(1)
         EsubM(kxJ3+lxrow, kyB2+lycol,k)=
     $- GG(VDG22LM(1,msb),VDG22L(1,msb),VDG22L(2,msb))/hc
     $+ cmrow*FGP(VDG12L(1,msb),  VDG12LM(1,msb))
      endif   !!! ---{ 1.eq.11 }
 91   CONTINUE

C***********************************************************************
C
C.. TENTH  EQUATION: EXPRESS J2L IN J1U AND J2U
C..                          (KXJ2L= 5, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
C
      IF (KXJ2L.GT.0) THEN
      DO 100 I=2,NRP1
      K = I

      INCLUDE 'tent.inc'
C----------------------------------------
C
         z1m  = (cs(i  )/csm(i-1))**iexj1
         z1p  = (cs(i)/csm(i))**iexj1

         HsubM(kxJ2L+lxrow, kyJ1+lycol,k)=
     $        -FGM(DG12L(i,msb)*z1m,DG12LM(i-1,msb))
         EsubM(kxJ2L+lxrow, kyJ1+lycol,k)=
     $        -FGP(DG12L(i,msb)*z1p, DG12LM(i,msb))

         BsubM(kxJ2L+lxrow, kxJ2U+lxcol,k)=
     $        -FF(DG22L(i,msb),DG22LM(i-1,msb), DG22LM(i,msb))
         AsubM(kxJ2L+lxrow, kxJ2U+lxcol,k)=-FFM(DG22LM(i-1,msb)) 
         CsubM(kxJ2L+lxrow, kxJ2U+lxcol,k)= -FFP(DG22LM(i,msb)) 

 100  CONTINUE
      ENDIF

C***********************************************************************
C   EQUATION FOR KXB2L: EXPRESS B2L IN B1U AND B2U
C   NEEDED FOR OPTION KJRER=5
C-----------------------------------------------------------------------
      IF (KXB2L.GT.0) THEN
      DO I=2,NRP1
         K = I

         INCLUDE 'tent.inc'
 
         HSUBM(KXB2L+LXROW, KYB2+LYCOL,K)=
     $        -FGM(DG22L(I,MSB),DG22LM(I-1,MSB))
         ESUBM(KXB2L+LXROW, KYB2+LYCOL,K)=
     $        -FGP(DG22L(I,MSB),DG22LM(I,MSB))

         IF (ABS(RM(MSA,2)).GT.0.1) THEN
         BSUBM(KXB2L+LXROW, KXB1+LXCOL,K)=
     $        -FF(DG12L(I,MSB),DG12LM(I-1,MSB),DG12LM(I,MSB))
         ASUBM(KXB2L+LXROW, KXB1+LXCOL,K)=-FFM(DG12LM(I-1,MSB)) 
         CSUBM(KXB2L+LXROW, KXB1+LXCOL,K)=-FFP(DG12LM(I,MSB)) 
         ELSE
         BSUBM(KXB2L+LXROW, KXB1+LXCOL,K)=
     $        -FF(T(I)*DG12L(I,MSB),TM(I-1)*DG12LM(I-1,MSB),
     $            TM(I)*DG12LM(I,MSB))
         ASUBM(KXB2L+LXROW, KXB1+LXCOL,K)=-FFM(TM(I-1)*DG12LM(I-1,MSB)) 
         CSUBM(KXB2L+LXROW, KXB1+LXCOL,K)=-FFP(TM(I)*DG12LM(I,MSB)) 
         ENDIF
      ENDDO    
      ENDIF

C***********************************************************************
C   EQUATION FOR KXB3L: EXPRESS B3L IN B3U 
C   NEEDED FOR OPTION KJRER=5
C-----------------------------------------------------------------------
      IF (KXB3L.GT.0) THEN
      DO I=2,NRP1
         K = I

         INCLUDE 'tent.inc'
 
         HSUBM(KXB3L+LXROW, KYB3+LYCOL,K)=
     $        -FGM(DG33L(I,MSB),DG33LM(I-1,MSB))
         ESUBM(KXB3L+LXROW, KYB3+LYCOL,K)=
     $        -FGP(DG33L(I,MSB),DG33LM(I,MSB))
      ENDDO    
      ENDIF

C-----------------------------------------------------------------------
C EQUATION FOR KYJRE1, DEFINED ON HALF-INTEGER MESH
C-----------------------------------------------------------------------
      IF (KJRER.EQ.6) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'
    
      IF (ETA.GT.0.) THEN
      CTMP1 = OMEGACI0*RESIST(I)*T(I)
      CTMP2 = OMEGACI0*RESIST(I+1)*T(I+1)
      CTMP3 = OMEGACI0*RESISM(I)*TM(I)
      FSUBM(KYJRE1+LYROW,KXJRE2L+LXCOL,I)= 
     &    -GF(CTMP1*FRENDB2(I,MSB),CTMP3*FRENDB2M(I,MSB))    
      GSUBM(KYJRE1+LYROW,KXJRE2L+LXCOL,I)= 
     &    -GF(CTMP2*FRENDB2(I+1,MSB),CTMP3*FRENDB2M(I,MSB))    
      FSUBM(KYJRE1+LYROW,KXJ2L+LXCOL,I)= 
     &     GF(CTMP1*FRENDB2(I,MSB),CTMP3*FRENDB2M(I,MSB))    
      GSUBM(KYJRE1+LYROW,KXJ2L+LXCOL,I)= 
     &     GF(CTMP2*FRENDB2(I+1,MSB),CTMP3*FRENDB2M(I,MSB))    
      FSUBM(KYJRE1+LYROW,KXB2L+LXCOL,I)= 
     &    -GF(CTMP1*FREJNREDB3(I,MSB),CTMP3*FREJNREDB3M(I,MSB))    
      GSUBM(KYJRE1+LYROW,KXB2L+LXCOL,I)= 
     &    -GF(CTMP2*FREJNREDB3(I+1,MSB),CTMP3*FREJNREDB3M(I,MSB))    

      CTMP1 = OMEGACI0*RESIST(I)*DPSIDS(I)
      CTMP2 = OMEGACI0*RESIST(I+1)*DPSIDS(I+1)
      CTMP3 = OMEGACI0*RESISM(I)*DPSIDSM(I)
      FSUBM(KYJRE1+LYROW,KXJRE3+LXCOL,I)= 
     &   GF(CTMP1*FRENR2G22DJ2B2(I,MSB),CTMP3*FRENR2G22DJ2B2M(I,MSB))
      GSUBM(KYJRE1+LYROW,KXJRE3+LXCOL,I)= 
     &   GF(CTMP2*FRENR2G22DJ2B2(I+1,MSB),CTMP3*FRENR2G22DJ2B2M(I,MSB))
      FSUBM(KYJRE1+LYROW,KXJ3+LXCOL,I)= 
     &  -GF(CTMP1*FRENR2G22DJ2B2(I,MSB),CTMP3*FRENR2G22DJ2B2M(I,MSB))
      GSUBM(KYJRE1+LYROW,KXJ3+LXCOL,I)= 
     &  -GF(CTMP2*FRENR2G22DJ2B2(I+1,MSB),CTMP3*FRENR2G22DJ2B2M(I,MSB))
      FSUBM(KYJRE1+LYROW,KXB3L+LXCOL,I)= 
     &   GF(CTMP1*FREJNREG22DJB3(I,MSB),CTMP3*FREJNREG22DJB3M(I,MSB))
      GSUBM(KYJRE1+LYROW,KXB3L+LXCOL,I)= 
     &   GF(CTMP2*FREJNREG22DJB3(I+1,MSB),CTMP3*FREJNREG22DJB3M(I,MSB))
      ENDIF

      CTMP1 = OMEGACI0
      IF (I.LE.NRES) CTMP1 = 0.
      FSUBM(KYJRE1+LYROW,KXV1+LXCOL,I)=CTMP1*
     &     GF(FRENREJ(I,MSB),FRENREJM(I,MSB))
      GSUBM(KYJRE1+LYROW,KXV1+LXCOL,I)=CTMP1*
     &     GF(FRENREJ(I+1,MSB),FRENREJM(I,MSB))

      ENDDO
      ENDIF

C-----------------------------------------------------------------------
C EQUATION FOR KXJRE2, DEFINED ON INTEGER MESH
C-----------------------------------------------------------------------
      IF (KJRER.EQ.6) THEN
      DO I=2,NRP1
      INCLUDE 'tent.inc'

      IF (ETA.GT.0.) THEN
      CTMP1 = OMEGACI0*RESIST(I)*T(I)
      CTMP2 = OMEGACI0*RESISM(I)*TM(I)
      CTMP3 = OMEGACI0*RESISM(I-1)*TM(I-1)
      Z1M  = (CS(I)/CSM(I-1))**IEXJ1
      Z1P  = (CS(I)/CSM(I))**IEXJ1
      HSUBM(KXJRE2+LXROW,KYJRE1+LYCOL,I)=
     &     FGM(CTMP1*Z1M*FRENG11DJB2(I,MSB),CTMP3*FRENG11DJB2M(I-1,MSB))
      ESUBM(KXJRE2+LXROW,KYJRE1+LYCOL,I)=
     &     FGP(CTMP1*Z1P*FRENG11DJB2(I,MSB),CTMP2*FRENG11DJB2M(I,MSB))
      HSUBM(KXJRE2+LXROW,KYJ1+LYCOL,I)=
     &    -FGM(CTMP1*Z1M*FRENG11DJB2(I,MSB),CTMP3*FRENG11DJB2M(I-1,MSB))
      ESUBM(KXJRE2+LXROW,KYJ1+LYCOL,I)=
     &    -FGP(CTMP1*Z1P*FRENG11DJB2(I,MSB),CTMP2*FRENG11DJB2M(I,MSB))

      BSUBM(KXJRE2+LXROW, KXJRE2+LXCOL,I)=
     &     FF(CTMP1*FRENG12DJB2(I,MSB),CTMP3*FRENG12DJB2M(I-1,MSB), 
     &        CTMP2*FRENG12DJB2M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXJRE2+LXCOL,I)= 
     &     FFM(CTMP3*FRENG12DJB2M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXJRE2+LXCOL,I)= 
     &     FFP(CTMP2*FRENG12DJB2M(I,MSB))
      BSUBM(KXJRE2+LXROW, KXJ2U+LXCOL,I)=
     &    -FF(CTMP1*FRENG12DJB2(I,MSB),CTMP3*FRENG12DJB2M(I-1,MSB), 
     &        CTMP2*FRENG12DJB2M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXJ2U+LXCOL,I)= 
     &    -FFM(CTMP3*FRENG12DJB2M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXJ2U+LXCOL,I)= 
     &    -FFP(CTMP2*FRENG12DJB2M(I,MSB))

      IF (ABS(RM(MSA,2)).GT.0.1) THEN
      BSUBM(KXJRE2+LXROW, KXB1+LXCOL,I)=
     &    FF(CTMP1*FREJNREG11DJB3(I,MSB),CTMP3*FREJNREG11DJB3M(I-1,MSB),
     &       CTMP2*FREJNREG11DJB3M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXB1+LXCOL,I)= 
     &    FFM(CTMP3*FREJNREG11DJB3M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXB1+LXCOL,I)= 
     &    FFP(CTMP2*FREJNREG11DJB3M(I,MSB))
      ELSE
      BSUBM(KXJRE2+LXROW, KXB1+LXCOL,I)=
     &    FF(CTMP1*T(I)*FREJNREG11DJB3(I,MSB),
     &       CTMP3*TM(I-1)*FREJNREG11DJB3M(I-1,MSB),
     &       CTMP2*TM(I)*FREJNREG11DJB3M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXB1+LXCOL,I)= 
     &    FFM(CTMP3*TM(I-1)*FREJNREG11DJB3M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXB1+LXCOL,I)= 
     &    FFP(CTMP2*TM(I)*FREJNREG11DJB3M(I,MSB))
      ENDIF

      HSUBM(KXJRE2+LXROW,KYB2+LYCOL,I)=
     &   FGM(CTMP1*FREJNREG12DJB3(I,MSB),CTMP3*FREJNREG12DJB3M(I-1,MSB))
      ESUBM(KXJRE2+LXROW,KYB2+LYCOL,I)=
     &   FGP(CTMP1*FREJNREG12DJB3(I,MSB),CTMP2*FREJNREG12DJB3M(I,MSB))

      CTMP1 = OMEGACI0*RESIST(I)*DPSIDS(I)
      CTMP2 = OMEGACI0*RESISM(I)*DPSIDSM(I)
      CTMP3 = OMEGACI0*RESISM(I-1)*DPSIDSM(I-1)
      BSUBM(KXJRE2+LXROW, KXJRE3+LXCOL,I)=
     &   -FF(CTMP1*FRENR2G12DJ2B2(I,MSB),CTMP3*FRENR2G12DJ2B2M(I-1,MSB),
     &       CTMP2*FRENR2G12DJ2B2M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXJRE3+LXCOL,I)= 
     &   -FFM(CTMP3*FRENR2G12DJ2B2M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXJRE3+LXCOL,I)= 
     &   -FFP(CTMP2*FRENR2G12DJ2B2M(I,MSB))
      BSUBM(KXJRE2+LXROW, KXJ3+LXCOL,I)=
     &    FF(CTMP1*FRENR2G12DJ2B2(I,MSB),CTMP3*FRENR2G12DJ2B2M(I-1,MSB),
     &       CTMP2*FRENR2G12DJ2B2M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXJ3+LXCOL,I)= 
     &    FFM(CTMP3*FRENR2G12DJ2B2M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXJ3+LXCOL,I)= 
     &    FFP(CTMP2*FRENR2G12DJ2B2M(I,MSB))

      BSUBM(KXJRE2+LXROW, KXB3L+LXCOL,I)=
     &   -FF(CTMP1*FREJNREG12DJB3(I,MSB),CTMP3*FREJNREG12DJB3M(I-1,MSB),
     &       CTMP2*FREJNREG12DJB3M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXB3L+LXCOL,I)= 
     &   -FFM(CTMP3*FREJNREG12DJB3M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXB3L+LXCOL,I)= 
     &   -FFP(CTMP2*FREJNREG12DJB3M(I,MSB))
      ENDIF

      TMP = 1.0
      IF (I.LE.NRES) TMP = 0.0
      CTMP1 = OMEGACI0*DPSIDS(I)**2*TMP
      CTMP2 = OMEGACI0*DPSIDSM(I)**2*TMP
      CTMP3 = OMEGACI0*DPSIDSM(I-1)**2*TMP
      BSUBM(KXJRE2+LXROW, KXV1+LXCOL,I)=
     &   -FF(CTMP1*FRENG12DJB2(I,MSB),CTMP3*FRENG12DJB2M(I-1,MSB),
     &       CTMP2*FRENG12DJB2M(I,MSB))
      ASUBM(KXJRE2+LXROW, KXV1+LXCOL,I)= 
     &   -FFM(CTMP3*FRENG12DJB2M(I-1,MSB))
      CSUBM(KXJRE2+LXROW, KXV1+LXCOL,I)= 
     &   -FFP(CTMP2*FRENG12DJB2M(I,MSB))

      CTMP1 = OMEGACI0*T(I)*TMP
      CTMP2 = OMEGACI0*TM(I)*TMP
      CTMP3 = OMEGACI0*TM(I-1)*TMP
      HSUBM(KXJRE2+LXROW,KYV2+LYCOL,I)=
     &     FGM(CTMP1*FRENJDB2(I,MSB),CTMP3*FRENJDB2M(I-1,MSB))
      ESUBM(KXJRE2+LXROW,KYV2+LYCOL,I)=
     &     FGP(CTMP1*FRENJDB2(I,MSB),CTMP2*FRENJDB2M(I,MSB))

      ENDDO
      ENDIF

C-----------------------------------------------------------------------
C EQUATION FOR KXJRE3, DEFINED ON INTEGER MESH
C-----------------------------------------------------------------------
      IF (KJRER.EQ.6) THEN
      DO I=2,NRP1
      INCLUDE 'tent.inc'

      IF (ETA.GT.0.) THEN
      CTMP1 = OMEGACI0*RESIST(I)*DPSIDS(I)
      CTMP2 = OMEGACI0*RESISM(I)*DPSIDSM(I)
      CTMP3 = OMEGACI0*RESISM(I-1)*DPSIDSM(I-1)
      Z1M  = (CS(I)/CSM(I-1))**IEXJ1
      Z1P  = (CS(I)/CSM(I))**IEXJ1
      HSUBM(KXJRE3+LXROW,KYJRE1+LYCOL,I)=
     &    -FGM(CTMP1*Z1M*FRENDB2R2(I,MSB),CTMP3*FRENDB2R2M(I-1,MSB))
      ESUBM(KXJRE3+LXROW,KYJRE1+LYCOL,I)=
     &    -FGP(CTMP1*Z1P*FRENDB2R2(I,MSB),CTMP2*FRENDB2R2M(I,MSB))
      HSUBM(KXJRE3+LXROW,KYJ1+LYCOL,I)=
     &     FGM(CTMP1*Z1M*FRENDB2R2(I,MSB),CTMP3*FRENDB2R2M(I-1,MSB))
      ESUBM(KXJRE3+LXROW,KYJ1+LYCOL,I)=
     &     FGP(CTMP1*Z1P*FRENDB2R2(I,MSB),CTMP2*FRENDB2R2M(I,MSB))

      IF (ABS(RM(MSA,2)).GT.0.1) THEN
      BSUBM(KXJRE3+LXROW, KXB1+LXCOL,I)=
     &   -FF(CTMP1*FREJNREDB3R2(I,MSB),CTMP3*FREJNREDB3R2M(I-1,MSB),
     &       CTMP2*FREJNREDB3R2M(I,MSB))
      ASUBM(KXJRE3+LXROW, KXB1+LXCOL,I)= 
     &   -FFM(CTMP3*FREJNREDB3R2M(I-1,MSB))
      CSUBM(KXJRE3+LXROW, KXB1+LXCOL,I)= 
     &   -FFP(CTMP2*FREJNREDB3R2M(I,MSB))
      ELSE
      BSUBM(KXJRE3+LXROW, KXB1+LXCOL,I)=
     &   -FF(CTMP1*T(I)*FREJNREDB3R2(I,MSB),
     &       CTMP3*TM(I-1)*FREJNREDB3R2M(I-1,MSB),
     &       CTMP2*TM(I)*FREJNREDB3R2M(I,MSB))
      ASUBM(KXJRE3+LXROW, KXB1+LXCOL,I)= 
     &   -FFM(CTMP3*TM(I-1)*FREJNREDB3R2M(I-1,MSB))
      CSUBM(KXJRE3+LXROW, KXB1+LXCOL,I)= 
     &   -FFP(CTMP2*TM(I)*FREJNREDB3R2M(I,MSB))
      ENDIF
      ENDIF

      TMP = 1.0
      IF (I.LE.NRES) TMP = 0.0
      CTMP1 = OMEGACI0*DPSIDS(I)*T(I)*TMP
      CTMP2 = OMEGACI0*DPSIDSM(I)*TM(I)*TMP
      CTMP3 = OMEGACI0*DPSIDSM(I-1)*TM(I-1)*TMP
      BSUBM(KXJRE3+LXROW, KXV1+LXCOL,I)=
     &   -FF(CTMP1*FRENG12DB2R2(I,MSB),CTMP3*FRENG12DB2R2M(I-1,MSB),
     &       CTMP2*FRENG12DB2R2M(I,MSB))
      ASUBM(KXJRE3+LXROW, KXV1+LXCOL,I)= 
     &   -FFM(CTMP3*FRENG12DB2R2M(I-1,MSB))
      CSUBM(KXJRE3+LXROW, KXV1+LXCOL,I)= 
     &   -FFP(CTMP2*FRENG12DB2R2M(I,MSB))

      CTMP1 = OMEGACI0*DPSIDS(I)*TMP
      CTMP2 = OMEGACI0*DPSIDSM(I)*TMP
      CTMP3 = OMEGACI0*DPSIDSM(I-1)*TMP
      HSUBM(KXJRE3+LXROW,KYV2+LYCOL,I)=
     &    -FGM(CTMP1*FRENG22DB2(I,MSB),CTMP3*FRENG22DB2M(I-1,MSB))
      ESUBM(KXJRE3+LXROW,KYV2+LYCOL,I)=
     &    -FGP(CTMP1*FRENG22DB2(I,MSB),CTMP2*FRENG22DB2M(I,MSB))

      ENDDO
      ENDIF

C-----------------------------------------------------------------------
C EQUATION FOR KXJRE2L, DEFINED ON INTEGER MESH
C-----------------------------------------------------------------------
      IF (KJRER.EQ.6.AND.KXJRE2L.GT.0) THEN
      DO I=2,NRP1
      INCLUDE 'tent.inc'
 
      Z1M  = (CS(I)/CSM(I-1))**IEXJ1
      Z1P  = (CS(I)/CSM(I))**IEXJ1
      HSUBM(KXJRE2L+LXROW,KYJRE1+LYCOL,I)=
     &     -FGM(DG12L(I,MSB)*Z1M,DG12LM(I-1,MSB))
      ESUBM(KXJRE2L+LXROW,KYJRE1+LYCOL,I)=
     &     -FGP(DG12L(I,MSB)*Z1P,DG12LM(I,MSB))

      BSUBM(KXJRE2L+LXROW,KXJRE2+LXCOL,I)=
     $    -FF(DG22L(I,MSB),DG22LM(I-1,MSB),DG22LM(I,MSB))
      ASUBM(KXJRE2L+LXROW,KXJRE2+LXCOL,I)=-FFM(DG22LM(I-1,MSB)) 
      CSUBM(KXJRE2L+LXROW,KXJRE2+LXCOL,I)=-FFP(DG22LM(I,MSB)) 
     
      ENDDO
      ENDIF

C***********************************************************************
C
C.. INTERMEDIATE VARIABLE PDE                    
C..                          (KXPD, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
C
      IF ((IDIAMB.EQ.1.OR.IDIAMB.EQ.3).AND.KXPD.GT.0) THEN
      DO 102 I=2,NRP1
         K = I

         INCLUDE 'tent.inc'
C----------------------------------------
         TMP   = B0K/OMEGACI0
         CTMP1 = TMP*DPSIDS(I)/RHO(I)     
         CTMP2 = TMP*DPSIDSM(I-1)/RHOM(I-1)     
         CTMP3 = TMP*DPSIDSM(I)/RHOM(I)     
         CTMP4 = TMP*T(I)/RHO(I)     
         CTMP5 = TMP*TM(I-1)/RHOM(I-1)     
         CTMP6 = TMP*TM(I)/RHOM(I)     
         HsubM(KXPD+lxrow, KYPE+lycol,k)=
     $     -cma*FGM(DJB2(i,msb)*CTMP1,DJB2M(i-1,msb)*CTMP2)
     $     -cna*FGM(DR2B2(i,msb)*CTMP4,DR2B2M(i-1,msb)*CTMP5)
         EsubM(KXPD+lxrow, KYPE+lycol,k)=
     $     -cma*FGP(DJB2(i,msb)*CTMP1,DJB2M(i,msb)*CTMP3)
     $     -cna*FGP(DR2B2(i,msb)*CTMP4,DR2B2M(i,msb)*CTMP6)
C        NOTE MINUS SIGN ABOVE SINCE THESE TERMS ARE MOVED TO LHS         
 102  CONTINUE
      ENDIF

C***********************************************************************
C
C.. INTERMEDIATE VARIABLE PED                    
C..                          (KXPED, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
C
      IF ((IDIAMB.EQ.1.OR.IDIAMB.EQ.3).AND.KXPED.GT.0) THEN
      DO I=2,NRP1
         K = I

         INCLUDE 'tent.inc'
C----------------------------------------
         TMP   = B0K/OMEGACI0*0.5
         CTMP1 = TMP*DPSIDS(I)/RHO(I)     
         CTMP2 = TMP*DPSIDSM(I-1)/RHOM(I-1)     
         CTMP3 = TMP*DPSIDSM(I)/RHOM(I)     

         CTMP  = FGM(DB2CDJB4(I,MSB)*CTMP1,DB2CDJB4M(I-1,MSB)*CTMP2)
         HSUBM(KXPED+LXROW, KYPP+LYCOL,K)=-CTMP
         HSUBM(KXPED+LXROW, KYPE+LYCOL,K)= CTMP

         CTMP  = FGP(DB2CDJB4(I,MSB)*CTMP1,DB2CDJB4M(I,MSB)*CTMP3)
         ESUBM(KXPED+LXROW, KYPP+LYCOL,K)=-CTMP
         ESUBM(KXPED+LXROW, KYPE+LYCOL,K)= CTMP
C        NOTE MINUS SIGN ABOVE SINCE THESE TERMS ARE MOVED TO LHS         
      ENDDO   
      ENDIF

 150  CONTINUE
C-----------------------------------------------------------------------
C.. SEVENTEENTH EQUATION: PERTURBED DENSITY FROM DIVERGENCE
C..                            (KYRHOP=11, DEFINED ON HALF  MESH)
C-----------------------------------------------------------------------
      IF (KJPKEY.GT.0.AND.KYRHOP.GT.0) THEN
      DO 115 I=1,NR
         K = I
         INCLUDE 'tophat.inc'
C
         zv2m = (cs(i  )/csm(i))**iexv2
         zv2p = (cs(i+1)/csm(i))**iexv2

         FsubM(KYRHOP+LYROW, KXV1+LXCOL,K)=
     $      RHO(i)*JACOBI(i,msb)*znorm
     $     +cmrow*GF(RHO(i)*G12B2B2(i,msb),RHOM(i)*G12B2B2M(i,msb))
     $     +cna*GF(RHO(i)*G12B2B3(i,msb),RHOM(i)*G12B2B3M(i,msb))

         GsubM(KYRHOP+LYROW, KXV1+LXCOL,K)=  
     $     -RHO(i+1)*JACOBI(i+1,msb)*znorm
     $     +cmrow*GF(RHO(i+1)*G12B2B2(i+1,msb),RHOM(i)*G12B2B2M(i,msb))
     $     +cna*GF(RHO(i+1)*G12B2B3(i+1,msb), RHOM(i)*G12B2B3M(i,msb))

         DsubM(KYRHOP+LYROW, KYV2+LYCOL,K)=
     $      -cmrow*GG(RHOM(i)*TB2M(i,msb),RHO(i)*TB2(i,msb)*zv2m, 
     $       RHO(i+1)*TB2(i+1,msb)*zv2p)
     $      +cna*GG(RHOM(i)*G22B2B2M(i,msb),RHO(i)*G22B2B2(i,msb)*zv2m,
     $       RHO(i+1)*G22B2B2(i+1,msb)*zv2p)

         IF (KYV3.GT.0) 
     $   DsubM(KYRHOP+LYROW, KYV3+LYCOL,K)=
     $      -cna*GG(RHOM(i)*B3j2M(i,msb)*TM(I), 
     $       RHO(i)*B3j2(i,msb)*T(I),RHO(i+1)*B3j2(i+1,msb)*T(I+1)) 

         DSUBM(KYRHOP+LYROW,KYRHOP+LYCOL,K)=SHIFTM(I)*
     $       GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))

         IF (NPROFRC.GT.0) 
     $   DSUBM(KYRHOP+LYROW,KYRHOP+LYCOL,K)=DSUBM(KYRHOP+LYROW,KYRHOP+
     $       LYCOL,K)-CNA*
     $       GG(RGROTCJM(I,MSB),RGROTCJ(I,MSB),RGROTCJ(I+1,MSB))
 115  CONTINUE
      ENDIF
 
C
C-----------------------------------------------------------------------
C.. ELEVENTH EQUATION: THERMAL ION PRESSURE FROM CONVECTION AND DIVERGENCE
C..                            (KYPR=6, DEFINED ON   HALF  MESH)
C-----------------------------------------------------------------------
      IF (KJPKEY.GT.0) THEN
      DO 110 I=1,NR
         K = I
         INCLUDE 'tophat.inc'
C
         zv2m = (cs(i  )/csm(i))**iexv2
         zv2p = (cs(i+1)/csm(i))**iexv2

         TMP = ALPHAP
         IF (KYPE.LE.0) TMP = 1. 

         CTMP1 = 0.
         CTMP2 = 0.
         CTMP3 = 0.
         CTMP4 = DPsiDs(i)*PPeq(i)*TMP
         CTMP5 = DPsiDs(i+1)*PPeq(i+1)*TMP
         CTMP6 = DPsiDsM(i)*PPeqM(i)*TMP

         CTMP7 = Peq(i)
         CTMP8 = Peq(i+1)
         CTMP9 = PeqM(i)
         IF (INCKIN.EQ.1.AND.IGAMMADIVV.EQ.1) THEN
         DO KP=1,NSPECIES
            IF (ABS(PSPECIES_NP(KP))+ABS(PSPECIES_NTB(KP))+
     &          ABS(PSPECIES_NTD(KP)).GT.0.) THEN
               CTMP7 = CTMP7 - ESPECIES_PRE(I,1,KP)
               CTMP8 = CTMP8 - ESPECIES_PRE(I+1,1,KP)
               CTMP9 = CTMP9 - ESPECIES_PRE(I,2,KP)

               CTMP4 = CTMP4 - DPsiDs(i)*  ESPECIES_PREP(I,1,KP)*  TMP
               CTMP5 = CTMP5 - DPsiDs(i+1)*ESPECIES_PREP(I+1,1,KP)*TMP
               CTMP6 = CTMP6 - DPsiDsM(i)* ESPECIES_PREP(I,2,KP)*  TMP
            ENDIF
         ENDDO
         ENDIF

         IF (INCKIN.EQ.1.AND.KFASTRUN.EQ.1.AND.IPERTURB.EQ.0) THEN
         DO KP=1,NSPECIES
            IF (ISPECIES_F0(KP).EQ.0.OR.ISPECIES_F0(KP).EQ.1.OR.
     &          ISPECIES_F0(KP).EQ.2) THEN
            CTMP1 = CTMP1+ESPECIES_PREP(I,1,KP)*PSPECIES_AP(KP)
            CTMP2 = CTMP2+ESPECIES_PREP(I+1,1,KP)*PSPECIES_AP(KP)
            CTMP3 = CTMP3+ESPECIES_PREP(I,2,KP)*PSPECIES_AP(KP)
            ENDIF
         ENDDO
         CTMP1 = DPsiDs(i)*CTMP1*TMP*ALPHAD
         CTMP2 = DPsiDs(i+1)*CTMP2*TMP*ALPHAD
         CTMP3 = DPsiDsM(i)*CTMP3*TMP*ALPHAD

         ELSEIF (INCKIN.NE.1.OR.IPERTURB.NE.0) THEN

         CTMP1 = DPsiDs(i)*PPeq(i)*TMP*TTCCONV0
         CTMP2 = DPsiDs(i+1)*PPeq(i+1)*TMP*TTCCONV0
         CTMP3 = DPsiDsM(i)*PPeqM(i)*TMP*TTCCONV0

         ENDIF
         IF (V2XKEY.EQ.0.OR.V2XKEY.EQ.1) THEN
         FsubM(kyPr+lyrow, kxV1+lxcol,k)=         
     &      -GF(CTMP1*JACOBI(i,msb),CTMP3*JACOBM(i,msb)) +
     &      gamarm(i)*GF(CTMP4*JACOBI(i,msb),CTMP6*JACOBM(i,msb)) +
     &      gamarm(i) *CTMP7*TMP*JACOBI(i,msb) *znorm
     $     +cmrow *gamarm(i)*TMP*GF(CTMP7*G12B2B2(i,msb),
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna *gamarm(i)*TMP*GF(CTMP7*G12B2B3(i,msb),
     &      CTMP9*G12B2B3M(i,msb))

         GsubM(kyPr+lyrow, kxV1+lxcol,k)= 
     &      -GF(CTMP2*JACOBI(i+1,msb),CTMP3*JACOBM(i,msb))+
     &      gamarm(i)*GF(CTMP5*JACOBI(i+1,msb),CTMP6*JACOBM(i,msb))-
     &      gamarm(i)*CTMP8*TMP*JACOBI(i+1,msb)*znorm
     $     +cmrow*gamarm(i)*TMP*GF(CTMP8*G12B2B2(i+1,msb), 
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna*gamarm(i)*TMP*GF(CTMP8*G12B2B3(i+1,msb), 
     &      CTMP9*G12B2B3M(i,msb))

         DsubM(kyPr+lyrow, kyV2+lycol,k)=-cmrow*gamarm(i)*TMP 
     &      *GG(CTMP9*TB2M(i,msb),CTMP7*TB2(i,msb)*zv2m, 
     &       CTMP8*TB2(i+1,msb)*zv2p)+cna *gamarm(i)*TMP
     &      *GG(CTMP9*G22B2B2M(i,msb),CTMP7*G22B2B2(i,msb)*zv2m,
     &       CTMP8*G22B2B2(i+1,msb)*zv2p)

         IF (KYV3.GT.0) 
     &   DsubM(kyPr+lyrow, kyV3+lycol,k)=
     &      -cna*gamarm(i)*TMP*GG(CTMP9*B3j2M(i,msb)*TM(I), 
     &       CTMP7*B3j2(i,msb)*T(I),CTMP8*B3j2(i+1,msb)*T(I+1)) 

         DSUBM(KYPR+LYROW,KYPR+LYCOL,k)=SHIFTM(I)*TTCINERT0*
     &       GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))

         IF (NPROFRC.GT.0)
     &   DSUBM(KYPR+LYROW,KYPR+LYCOL,k)=DSUBM(KYPR+LYROW,KYPR+LYCOL,k)
     &      - CNA*GG(RGROTCJM(I,MSB),RGROTCJ(I,MSB),RGROTCJ(I+1,MSB))

         IF (KXW1.GT.0)
     &   FsubM(kyPr+lyrow, kxW1+lxcol,k)=         
     &     -TTCPERPI(i)*JACOBI(i,msb)*znorm
     &     -cmrow*GF(TTCPERPI(i)*G12B2B2(i,msb),
     &               TTCPERPM(i)*G12B2B2M(i,msb))
     &     -cna  *GF(TTCPERPI(i)*G12B2B3(i,msb),
     &               TTCPERPM(i)*G12B2B3M(i,msb))

         IF (KXW1.GT.0)
     &   GsubM(kyPr+lyrow, kxW1+lxcol,k)= 
     &      TTCPERPI(i+1)*JACOBI(i+1,msb)*znorm
     &     -cmrow*GF(TTCPERPI(i+1)*G12B2B2(i+1,msb), 
     &               TTCPERPM(i)*G12B2B2M(i,msb))
     &     -cna*  GF(TTCPERPI(i+1)*G12B2B3(i+1,msb), 
     &               TTCPERPM(i)*G12B2B3M(i,msb))

         IF (KYW2.GT.0)
     &   DsubM(kyPr+lyrow, kyW2+lycol,k)=
     &      cmrow*GG(TTCPERPM(i)*TB2M(i,msb),
     &               TTCPERPI(i)*TB2(i,msb)*zv2m, 
     &               TTCPERPI(i+1)*TB2(i+1,msb)*zv2p)
     &     -cna  *GG(TTCPERPM(i)*G22B2B2M(i,msb),
     &               TTCPERPI(i)*G22B2B2(i,msb)*zv2m,
     &               TTCPERPI(i+1)*G22B2B2(i+1,msb)*zv2p)

         IF (KYW3.GT.0) THEN
         CTMP10 = TTCPARAM(I) + TTCPERPM(I)
         CTMP11 = TTCPARAI(I) + TTCPERPI(I)
         CTMP12 = TTCPARAI(I+1) + TTCPERPI(I+1)
         DsubM(kyPr+lyrow, kyW3+lycol,k)=
     &      cna*  GG(B3j2M(I,MSB)*TM(I)*CTMP10,
     &               B3J2(I,MSB)*T(I)*CTMP11,
     &               B3J2(I+1,MSB)*T(I+1)*CTMP12)
     &      +     GG(V3PKM(I,MSB)*TTCPARAM(I),
     &               V3PK(I,MSB)*TTCPARAI(I),
     &               V3PK(I+1,MSB)*TTCPARAI(I+1))  
         ENDIF

         IF (ABS(TTCPARA0)+ABS(TTCPERP0).GT.0.) THEN
         IF (ABS(RM(MSA,2)).GT.0.1) THEN
            CTMP10 = TTCPARAI(I)*PPEQ(I)*DPSIDS(I)
            CTMP11 = TTCPARAI(I+1)*PPEQ(I+1)*DPSIDS(I+1)
            CTMP12 = TTCPARAM(I)*PPEQM(I)*DPSIDSM(I)
         ELSE
            CTMP10 = TTCPARAI(I)*PPEQ(I)*DPSIDS(I)*T(I)
            CTMP11 = TTCPARAI(I+1)*PPEQ(I+1)*DPSIDS(I+1)*T(I+1)
            CTMP12 = TTCPARAM(I)*PPEQM(I)*DPSIDSM(I)*TM(I)
         ENDIF
         FsubM(kyPr+lyrow, kxB1+lxcol,k)=
     &      cmrow*GF(DJB2(i,msb)   *CTMP10*DPSIDS(I),
     &               DJB2M(i,msb)  *CTMP12*DPSIDSM(I))    
     &     +cna  *GF(DR2B2(i,msb)  *CTMP10*T(I),
     &               DR2B2M(i,msb) *CTMP12*TM(I))    
     &     +      GF(TTCJB4(i,msb)  *CTMP10*DPSIDS(I),
     &               TTCJB4M(i,msb) *CTMP12*DPSIDSM(I)) 
         GsubM(kyPr+lyrow, kxB1+lxcol,k)=
     &      cmrow*GF(DJB2(i+1,msb) *CTMP11*DPSIDS(I+1), 
     &               DJB2M(i,msb)  *CTMP12*DPSIDSM(I))
     &     +cna  *GF(DR2B2(i+1,msb)*CTMP11*T(I+1),
     &               DR2B2M(i,msb) *CTMP12*TM(I))    
     &     +      GF(TTCJB4(i+1,msb)*CTMP11*DPSIDS(I+1),
     &               TTCJB4M(i,msb) *CTMP12*DPSIDSM(I))    
         ENDIF

         ELSE

         FsubM(kyPr+lyrow, kxX1+lxcol,k)=
     &      -GF(CTMP1*JACOBI(i,msb),CTMP3*JACOBM(i,msb)) +
     &      gamarm(i)*GF(CTMP4*JACOBI(i,msb),CTMP6*JACOBM(i,msb)) +
     &      gamarm(i) *CTMP7*TMP*JACOBI(i,msb) *znorm
     $     +cmrow *gamarm(i)*TMP*GF(CTMP7*G12B2B2(i,msb),
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna *gamarm(i)*TMP*GF(CTMP7*G12B2B3(i,msb),
     &      CTMP9*G12B2B3M(i,msb))

         GsubM(kyPr+lyrow, kxX1+lxcol,k)=  
     &      -GF(CTMP2*JACOBI(i+1,msb),CTMP3*JACOBM(i,msb))+
     &      gamarm(i)*GF(CTMP5*JACOBI(i+1,msb),CTMP6*JACOBM(i,msb))-
     &      gamarm(i)*CTMP8*TMP*JACOBI(i+1,msb)*znorm
     $     +cmrow*gamarm(i)*TMP*GF(CTMP8*G12B2B2(i+1,msb), 
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna*gamarm(i)*TMP*GF(CTMP8*G12B2B3(i+1,msb), 
     &      CTMP9*G12B2B3M(i,msb))

         DsubM(kyPr+lyrow, kyX2+lycol,k)=-cmrow*gamarm(i)*TMP 
     &      *GG(CTMP9*TB2M(i,msb),CTMP7*TB2(i,msb)*zv2m, 
     &       CTMP8*TB2(i+1,msb)*zv2p)+cna *gamarm(i)*TMP
     &      *GG(CTMP9*G22B2B2M(i,msb),CTMP7*G22B2B2(i,msb)*zv2m,
     &       CTMP8*G22B2B2(i+1,msb)*zv2p)

         IF (kyX3.GT.0)
     &   DsubM(kyPr+lyrow, kyX3+lycol,k)=
     &      -cna*gamarm(i)*TMP*GG(CTMP9*B3j2M(i,msb)*TM(I), 
     &       CTMP7*B3j2(i,msb)*T(I),CTMP8*B3j2(i+1,msb)*T(I+1)) 

         DSUBM(KYPR+LYROW,KYPR+LYCOL,k)=
     &      -GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))         
         ENDIF

C        DIAMAGNETIC TERM
         IF (IDIAMTI.EQ.1) THEN
            TMP   = 2.0*gamarm(i)*B0K/OMEGACI0*FDIAMTI
            CTMP1 = TMP*ESPECIES_TEM(I,1,1)
            CTMP2 = TMP*ESPECIES_TEM(I+1,1,1)
            CTMP3 = TMP*ESPECIES_TEM(I,2,1)
            DsubM(kyPr+lyrow, kyPr+lycol,k)=
     &       DsubM(kyPr+lyrow, kyPr+lycol,k)
     &       -cma*GG(DJFBM(i,msb)*CTMP3, 
     &        DJFB(i,msb)*CTMP1,DJFB(i+1,msb)*CTMP2) 
     &       -cna*GG(DGJBBM(i,msb)*CTMP3, 
     &        DGJBB(i,msb)*CTMP1,DGJBB(i+1,msb)*CTMP2) 
         ENDIF
 110  CONTINUE
      ENDIF

C     USE EQUATION FOR KYPR TO REPRESENT DIV(V)
      IF (KJPKEY.EQ.0) THEN
      DO 111 I=1,NR
         K = I
         INCLUDE 'tophat.inc'
C
         zv2m = (cs(i  )/csm(i))**iexv2
         zv2p = (cs(i+1)/csm(i))**iexv2

         FsubM(kyPr+lyrow, kxV1+lxcol,k)=
     &     -JACOBI(i,msb) *znorm
     $     -cmrow *GF(G12B2B2(i,msb),
     &      G12B2B2M(i,msb))
     $     -cna *GF(G12B2B3(i,msb),
     &      G12B2B3M(i,msb))

         GsubM(kyPr+lyrow, kxV1+lxcol,k)=  
     &     +JACOBI(i+1,msb) *znorm
     $     -cmrow *GF(G12B2B2(i+1,msb), 
     &      G12B2B2M(i,msb))
     $     -cna *GF(G12B2B3(i+1,msb), 
     &      G12B2B3M(i,msb))

         DsubM(kyPr+lyrow, kyV2+lycol,k)=cmrow * 
     &       GG(TB2M(i,msb), TB2(i,msb)*zv2m, 
     &       TB2(i+1,msb)*zv2p)
     &      -cna *
     &       GG(G22B2B2M(i,msb), G22B2B2(i,msb)*zv2m,
     &       G22B2B2(i+1,msb)*zv2p)
 111  CONTINUE
      ENDIF

      IF (KJPKEY.EQ.0) GOTO 160
C
C-----------------------------------------------------------------------
C.. EQUATION: THERMAL ELECTRON PRESSURE FROM CONVECTION AND DIVERGENCE
C..                            (KYPE, DEFINED ON   HALF  MESH)
C-----------------------------------------------------------------------
      IF (KYPE.GT.0.AND.INCKIN.EQ.0) THEN
      DO 250 I=1,NR
         K = I
         INCLUDE 'tophat.inc'

         IF (IDIAMTE.EQ.2) THEN

         DSUBM(KYPE+LYROW,KYPE+LYCOL,k)=
     &       GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))

         zv2m = (cs(i  )/csm(i))**iexv2
         zv2p = (cs(i+1)/csm(i))**iexv2

         TMP   = (1.-ALPHAP)

         CTMP1 = 0.
         CTMP2 = 0.
         CTMP3 = 0.
         CTMP4 = DPsiDs(i)*PPeq(i)*TMP
         CTMP5 = DPsiDs(i+1)*PPeq(i+1)*TMP
         CTMP6 = DPsiDsM(i)*PPeqM(i)*TMP

         CTMP7 = Peq(i)
         CTMP8 = Peq(i+1)
         CTMP9 = PeqM(i)
         IF (INCKIN.EQ.1.AND.IGAMMADIVV.EQ.1) THEN
         DO KP=1,NSPECIES
            IF (ABS(PSPECIES_NP(KP))+ABS(PSPECIES_NTB(KP))+
     &          ABS(PSPECIES_NTD(KP)).GT.0.) THEN
               CTMP7 = CTMP7 - ESPECIES_PRE(I,1,KP)
               CTMP8 = CTMP8 - ESPECIES_PRE(I+1,1,KP)
               CTMP9 = CTMP9 - ESPECIES_PRE(I,2,KP)

               CTMP4 = CTMP4 - DPsiDs(i)*  ESPECIES_PREP(I,1,KP)*  TMP
               CTMP5 = CTMP5 - DPsiDs(i+1)*ESPECIES_PREP(I+1,1,KP)*TMP
               CTMP6 = CTMP6 - DPsiDsM(i)* ESPECIES_PREP(I,2,KP)*  TMP
            ENDIF
         ENDDO
         ENDIF

         IF (INCKIN.EQ.1.AND.KFASTRUN.EQ.1.AND.IPERTURB.EQ.0) THEN
         DO KP=1,NSPECIES
            IF (ISPECIES_F0(KP).EQ.0.OR.ISPECIES_F0(KP).EQ.1.OR.
     &          ISPECIES_F0(KP).EQ.2) THEN
            CTMP1 = CTMP1+ESPECIES_PREP(I,1,KP)*PSPECIES_AP(KP)
            CTMP2 = CTMP2+ESPECIES_PREP(I+1,1,KP)*PSPECIES_AP(KP)
            CTMP3 = CTMP3+ESPECIES_PREP(I,2,KP)*PSPECIES_AP(KP)
            ENDIF
         ENDDO
         CTMP1 = DPsiDs(i)*CTMP1*TMP*ALPHAD
         CTMP2 = DPsiDs(i+1)*CTMP2*TMP*ALPHAD
         CTMP3 = DPsiDsM(i)*CTMP3*TMP*ALPHAD

         ELSEIF (INCKIN.NE.1.OR.IPERTURB.NE.0) THEN

         CTMP1 = DPsiDs(i)*PPeq(i)*TMP*TTCCONV0
         CTMP2 = DPsiDs(i+1)*PPeq(i+1)*TMP*TTCCONV0
         CTMP3 = DPsiDsM(i)*PPeqM(i)*TMP*TTCCONV0

         ENDIF
         IF (V2XKEY.EQ.0 .OR. V2XKEY.EQ.1) THEN
         FsubM(kyPe+lyrow, kxV1+lxcol,k)=         
     &      -GF(CTMP1*JACOBI(i,msb),CTMP3*JACOBM(i,msb)) +
     &      gamarm(i)*GF(CTMP4*JACOBI(i,msb),CTMP6*JACOBM(i,msb)) +
     &      gamarm(i) *CTMP7*TMP*JACOBI(i,msb) *znorm
     $     +cmrow *gamarm(i)*TMP*GF(CTMP7*G12B2B2(i,msb),
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna *gamarm(i)*TMP*GF(CTMP7*G12B2B3(i,msb),
     &      CTMP9*G12B2B3M(i,msb))

         GsubM(kyPe+lyrow, kxV1+lxcol,k)= 
     &      -GF(CTMP2*JACOBI(i+1,msb),CTMP3*JACOBM(i,msb))+
     &      gamarm(i)*GF(CTMP5*JACOBI(i+1,msb),CTMP6*JACOBM(i,msb))-
     &      gamarm(i)*CTMP8*TMP*JACOBI(i+1,msb)*znorm
     $     +cmrow*gamarm(i)*TMP*GF(CTMP8*G12B2B2(i+1,msb), 
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna*gamarm(i)*TMP*GF(CTMP8*G12B2B3(i+1,msb), 
     &      CTMP9*G12B2B3M(i,msb))

         DsubM(kyPe+lyrow, kyV2+lycol,k)=-cmrow*gamarm(i)*TMP 
     &      *GG(CTMP9*TB2M(i,msb),CTMP7*TB2(i,msb)*zv2m, 
     &       CTMP8*TB2(i+1,msb)*zv2p)+cna *gamarm(i)*TMP
     &      *GG(CTMP9*G22B2B2M(i,msb),CTMP7*G22B2B2(i,msb)*zv2m,
     &       CTMP8*G22B2B2(i+1,msb)*zv2p)

         IF (KYV3.GT.0) 
     &   DsubM(kyPe+lyrow, kyV3+lycol,k)=
     &      -cna*gamarm(i)*TMP*GG(CTMP9*B3j2M(i,msb)*TM(I), 
     &       CTMP7*B3j2(i,msb)*T(I),CTMP8*B3j2(i+1,msb)*T(I+1)) 

         DSUBM(KYPE+LYROW,KYPE+LYCOL,k)=SHIFTM(I)*
     &       GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))
         IF (NPROFRC.GT.0)
     &   DSUBM(KYPE+LYROW,KYPE+LYCOL,k)=DSUBM(KYPE+LYROW,KYPE+LYCOL,k)
     &       - CNA*GG(RGROTCJM(I,MSB),RGROTCJ(I,MSB),RGROTCJ(I+1,MSB))
         ELSE
         FsubM(kyPe+lyrow, kxX1+lxcol,k)=
     &      -GF(CTMP1*JACOBI(i,msb),CTMP3*JACOBM(i,msb)) +
     &      gamarm(i)*GF(CTMP4*JACOBI(i,msb),CTMP6*JACOBM(i,msb)) +
     &      gamarm(i) *CTMP7*TMP*JACOBI(i,msb) *znorm
     $     +cmrow *gamarm(i)*TMP*GF(CTMP7*G12B2B2(i,msb),
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna *gamarm(i)*TMP*GF(CTMP7*G12B2B3(i,msb),
     &      CTMP9*G12B2B3M(i,msb))

         GsubM(kyPe+lyrow, kxX1+lxcol,k)=  
     &      -GF(CTMP2*JACOBI(i+1,msb),CTMP3*JACOBM(i,msb))+
     &      gamarm(i)*GF(CTMP5*JACOBI(i+1,msb),CTMP6*JACOBM(i,msb))-
     &      gamarm(i)*CTMP8*TMP*JACOBI(i+1,msb)*znorm
     $     +cmrow*gamarm(i)*TMP*GF(CTMP8*G12B2B2(i+1,msb), 
     &      CTMP9*G12B2B2M(i,msb))
     $     +cna*gamarm(i)*TMP*GF(CTMP8*G12B2B3(i+1,msb), 
     &      CTMP9*G12B2B3M(i,msb))

         DsubM(kyPe+lyrow, kyX2+lycol,k)=-cmrow*gamarm(i)*TMP 
     &      *GG(CTMP9*TB2M(i,msb),CTMP7*TB2(i,msb)*zv2m, 
     &       CTMP8*TB2(i+1,msb)*zv2p)+cna *gamarm(i)*TMP
     &      *GG(CTMP9*G22B2B2M(i,msb),CTMP7*G22B2B2(i,msb)*zv2m,
     &       CTMP8*G22B2B2(i+1,msb)*zv2p)

         DsubM(kyPe+lyrow, kyX3+lycol,k)=
     &      -cna*gamarm(i)*TMP*GG(CTMP9*B3j2M(i,msb)*TM(I), 
     &       CTMP7*B3j2(i,msb)*T(I),CTMP8*B3j2(i+1,msb)*T(I+1)) 

         DSUBM(KYPE+LYROW,KYPE+LYCOL,k)=
     &       -GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))         
         ENDIF

         ENDIF

C        DIAMAGNETIC TERM
         IF (IDIAMTE.EQ.1) THEN
            DsubM(kyPe+lyrow, kyPe+lycol,k)=
     &       DsubM(kyPe+lyrow, kyPe+lycol,k)
     &       +cna*GG(B3j2M(i,msb)*TM(i), 
     &        B3j2(i,msb)*T(i),B3j2(i+1,msb)*T(i+1)) 
         ENDIF
 250  CONTINUE
      ENDIF

C-----------------------------------------------------------------------
C.. TWELFTH EQUATION: COVARIANT-S-COMP OF EQ. OF DISPLACEMENT
C..                              (KXX1 = 6, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
      IF (KXX1.GT.0.AND.NPROFRC.GT.0) THEN
      DO I=2,NRP1
      K = I
C
      INCLUDE 'tent.inc'
      
      BsubM(KXX1+lxrow,KXX1+lxcol,k)=BsubM(KXX1+lxrow,KXX1+lxcol,k)
     $  -CNA*FF(RGROTC(i,msb),RGROTCM(i-1,msb),RGROTCM(i,msb))
      AsubM(KXX1+lxrow,KXX1+lxcol,k)=AsubM(KXX1+lxrow,KXX1+lxcol,k)
     $  -CNA*FFM(RGROTCM(i-1,msb))
      CsubM(KXX1+lxrow,KXX1+lxcol,k)=CsubM(KXX1+lxrow,KXX1+lxcol,k)
     $  -CNA*FFP(RGROTCM(i,msb))

      ENDDO
      ENDIF
C     GLX------<
      IF (KXX1.GT.0.AND.NPROFRP.GT.0) THEN
      DO I=2,NRP1
      K = I
      INCLUDE 'tent.inc'  
      BsubM(KXX1+lxrow,KXX1+lxcol,k)=BsubM(KXX1+lxrow,KXX1+lxcol,k)
     $  - CMA*FF(RGQ1RQ1P(I,MSB),RGQ1RQ1PM(I-1,MSB),RGQ1RQ1PM(I,MSB))
     $  - CNA*FF(RGQ1RQ2P(I,MSB),RGQ1RQ2PM(I-1,MSB),RGQ1RQ2PM(I,MSB))
      AsubM(KXX1+lxrow,KXX1+lxcol,k)=AsubM(KXX1+lxrow,KXX1+lxcol,k)
     $  - CMA*FFM(RGQ1RQ1PM(I-1,MSB))
     $  - CNA*FFM(RGQ1RQ2PM(I-1,MSB))
      CsubM(KXX1+lxrow,KXX1+lxcol,k)=CsubM(KXX1+lxrow,KXX1+lxcol,k)
     $  - CMA*FFP(RGQ1RQ1PM(I,MSB))
     $  - CNA*FFP(RGQ1RQ2PM(I,MSB))
      ENDDO
      ENDIF
C     --------->

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C.. COEFFICIENTS FOR THIRTEENTH EQUATION: COVARIANT-2-COMP. OF DISPLACEMENT
C..                                  (KYX2 = 7, DEFINED ON HALF MESH)
C..                                       COVARIANT-3-COMP. OF DISPLACEMENT
C..                                  (KYX3 = 12, DEFINED ON HALF MESH)
C-----------------------------------------------------------------------
C
      IEXE = -1
C
      DO 130 I=1,NR
      K = I
C
      INCLUDE 'tophat.inc'
C
C.....EQUATION MULTIPLIED BY 1/S.
C
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
C
C-----------------------------------------------------------------------
C.. INERTIAL TERMS
C-----------------------------------------------------------------------
      IF (KYX2.GT.0) THEN
      DSUBM(KYX2+LYROW,KYX2+LYCOL,K) = DSUBM(KYX2+LYROW,KYX2+LYCOL,K)  
     &     + GG(SHIFTM(I)*JACOBM(I,MSB),
     &                ZEM*ZV2M*SHIFTC(I)*JACOBI(I,MSB),
     &                ZEP*ZV2P*SHIFTC(I+1)*JACOBI(I+1,MSB))
      DSUBM(KYX2+LYROW,KYV2+LYCOL,K) = DSUBM(KYX2+LYROW,KYV2+LYCOL,K)  
     &     + GG(JACOBM(I,MSB),ZEM*ZV2M*JACOBI(I,MSB),
     &          ZEP*ZV2P*JACOBI(I+1,MSB))
      IF (NPROFRC.GT.0) 
     &DSUBM(KYX2+LYROW,KYX2+LYCOL,K) = DSUBM(KYX2+LYROW,KYX2+LYCOL,K)  
     &     - CNA*GG(RGROTCJM(I,MSB),ZEM*ZV2M*RGROTCJ(I,MSB),
     &              ZEP*ZV2P*RGROTCJ(I+1,MSB))

C     LLI
      IF (NPROFRC.GT.0) THEN
      IF (KXX1.GT.0) THEN  
      FSUBM(KYX2+LYROW,KXX1+LXCOL,K) = FSUBM(KYX2+LYROW,KXX1+LXCOL,K)
     &     - GF(ZEM*RGX2RX1C(I,MSB),RGX2RX1CM(I,MSB))
      GSUBM(KYX2+LYROW,KXX1+LXCOL,K) = GSUBM(KYX2+LYROW,KXX1+LXCOL,K)
     &     - GF(ZEP*RGX2RX1C(I+1,MSB),RGX2RX1CM(I,MSB))
      ENDIF
      
      DSUBM(KYX2+LYROW,KYX2+LYCOL,K) = DSUBM(KYX2+LYROW,KYX2+LYCOL,K)  
     &     - GG(RGX2RX2CM(I,MSB),ZEM*ZV2M*RGX2RX2C(I,MSB),
     &                           ZEP*ZV2P*RGX2RX2C(I+1,MSB))

      IF (KYX3.GT.0) THEN  
      DSUBM(KYX2+LYROW,KYX3+LYCOL,I) = DSUBM(KYX2+LYROW,KYX3+LYCOL,I)
     &     - GG(RGX2RX3CM(I,MSB),ZEM*RGX2RX3C(I,MSB),
     &          ZEP*RGX2RX3C(I+1,MSB))
      ENDIF
      ENDIF
C     GLX------<
      IF (NPROFRP.GT.0) THEN
      IF (KXX1.GT.0) THEN  
      FSUBM(KYX2+LYROW,KXX1+LXCOL,K) = FSUBM(KYX2+LYROW,KXX1+LXCOL,K)
     &   + GF(ZEM*RGX2RX1P(I,MSB),RGX2RX1PM(I,MSB))
      GSUBM(KYX2+LYROW,KXX1+LXCOL,K) = GSUBM(KYX2+LYROW,KXX1+LXCOL,K)
     &   + GF(ZEP*RGX2RX1P(I+1,MSB),RGX2RX1PM(I,MSB))
      ENDIF
      
      CTMP1=RHOUM(i)*TM(i)*B3J2M(i,msb) 
      CTMP2=RHOU(i)*T(i)*B3J2(i,msb) 
      CTMP3=RHOU(i+1)*T(i+1)*B3J2(i+1,msb)
 
      DSUBM(KYX2+LYROW,KYX2+LYCOL,K) = DSUBM(KYX2+LYROW,KYX2+LYCOL,K)  
     &   - CNA*GG(CTMP1,ZEM*ZV2M*CTMP2,ZEP*ZV2P*CTMP3)
     &   -     GG(RGX2RX2PM(I,MSB),ZEM*ZV2M*RGX2RX2P(I,MSB),
     &                             ZEP*ZV2P*RGX2RX2P(I+1,MSB))
      ENDIF
C     --------->
      ENDIF

      IF (KYX3.GT.0.AND.KYV3.GT.0) THEN
      DSUBM(KYX3+LYROW,KYX3+LYCOL,K) = DSUBM(KYX3+LYROW,KYX3+LYCOL,K)  
     &     + GG(SHIFTM(I)*JACOBM(I,MSB),
     &                ZEM*SHIFTC(I)*JACOBI(I,MSB),
     &                ZEP*SHIFTC(I+1)*JACOBI(I+1,MSB))
      DSUBM(KYX3+LYROW,KYV3+LYCOL,K) = DSUBM(KYX3+LYROW,KYV3+LYCOL,K)  
     &     + GG(JACOBM(I,MSB),ZEM*JACOBI(I,MSB),
     &                        ZEP*JACOBI(I+1,MSB))    
      IF (NPROFRC.GT.0)
     &DSUBM(KYX3+LYROW,KYX3+LYCOL,K) = DSUBM(KYX3+LYROW,KYX3+LYCOL,K)  
     &     - CNA*GG(RGROTCJM(I,MSB),ZEM*RGROTCJ(I,MSB),
     &                              ZEP*RGROTCJ(I+1,MSB))
      IF (KXX1.GT.0) THEN
      FSUBM(KYX3+LYROW,KXX1+LXCOL,K) = FSUBM(KYX3+LYROW,KXX1+LXCOL,K)
     &     + GF(ZEM*TB2(I,MSB)*DROT(I),TB2M(I,MSB)*DROTM(I))
      GSUBM(KYX3+LYROW,KXX1+LXCOL,K) = GSUBM(KYX3+LYROW,KXX1+LXCOL,K)
     &     + GF(ZEP*TB2(I+1,MSB)*DROT(I+1),TB2M(I,MSB)*DROTM(I))
      ENDIF

C     LLI
      IF (NPROFRC.GT.0) THEN
      IF (KXX1.GT.0) THEN   
      FSUBM(KYX3+LYROW,KXX1+LXCOL,K) = FSUBM(KYX3+LYROW,KXX1+LXCOL,K)
     &     + GF(ZEM*RGX3RX1C(I,MSB),RGX3RX1CM(I,MSB))
      GSUBM(KYX3+LYROW,KXX1+LXCOL,K) = GSUBM(KYX3+LYROW,KXX1+LXCOL,K)
     &     + GF(ZEP*RGX3RX1C(I+1,MSB),RGX3RX1CM(I,MSB))
      ENDIF

      IF (KYX2.GT.0) THEN   
      DSUBM(KYX3+LYROW,KYX2+LYCOL,K)=DSUBM(KYX3+LYROW,KYX2+LYCOL,K) 
     &     +GG(RGX3RX2CM(I,MSB),ZV2M*ZEM*RGX3RX2C(I,MSB),
     &         ZV2P*ZEP*RGX3RX2C(I+1,MSB))
      ENDIF

      DSUBM(KYX3+LYROW,KYX3+LYCOL,K) = DSUBM(KYX3+LYROW,KYX3+LYCOL,K) 
     &     +GG(RGX2RX2CM(I,MSB),ZEM*RGX2RX2C(I,MSB),
     &         ZEP*RGX2RX2C(I+1,MSB))

      ENDIF
C     GLX------<
      IF (NPROFRP.GT.0) THEN
      IF (KXX1.GT.0) THEN  
      CTMP1=RHOU(I)*(-2.0*V1PK(I,MSB)-TB2(I,MSB)/T(I)*
     &               PPEQ(I)*DPSIDS(I))+JACOBI(I,MSB)*DRHOU(I) 
      CTMP2=RHOUM(I)*(-2.0*V1PKM(I,MSB)-TB2M(I,MSB)/TM(I)*
     &               PPEQM(I)*DPSIDSM(I))+JACOBM(I,MSB)*DRHOUM(I) 
      CTMP3=RHOU(I+1)*(-2.0*V1PK(I+1,MSB)-TB2(I+1,MSB)/T(I+1)*
     &               PPEQ(I+1)*DPSIDS(I+1))+JACOBI(I+1,MSB)*DRHOU(I+1) 
      FSUBM(KYX3+LYROW,KXX1+LXCOL,K) = FSUBM(KYX3+LYROW,KXX1+LXCOL,K)
     &   + GF(ZEM*CTMP1,CTMP2)
      GSUBM(KYX3+LYROW,KXX1+LXCOL,K) = GSUBM(KYX3+LYROW,KXX1+LXCOL,K)
     &   + GF(ZEP*CTMP3,CTMP2)
      ENDIF
      IF (KYX2.GT.0) THEN 
      CTMP4=-2.0*RHOUM(I)*V2PKM(I,MSB)
      CTMP5=-2.0*RHOU(I)*V2PK(I,MSB)
      CTMP6=-2.0*RHOU(I+1)*V2PK(I+1,MSB)
      DSUBM(KYX3+LYROW,KYX2+LYCOL,K)=DSUBM(KYX3+LYROW,KYX2+LYCOL,K) 
     &   + GG(CTMP4,ZV2M*ZEM*CTMP5,ZV2P*ZEP*CTMP6)
      ENDIF

      CTMP1=RHOUM(i)*TM(i)*B3J2M(i,msb) 
      CTMP2=RHOU(i)*T(i)*B3J2(i,msb) 
      CTMP3=RHOU(i+1)*T(i+1)*B3J2(i+1,msb)
      DSUBM(KYX3+LYROW,KYX3+LYCOL,K) = DSUBM(KYX3+LYROW,KYX3+LYCOL,K) 
     &   - CNA*GG(CTMP1,ZEM*CTMP2,ZEP*CTMP3)
      ENDIF
C     --------->
      ENDIF

 130  CONTINUE
C
 
C-----------------------------------------------------------------------
C.. EQUATION FOR KXW1, DEFINED ON INTEGER MESH
C-----------------------------------------------------------------------
      IF (KXW1.GT.0.AND.KYW2.GT.0) THEN
      DO 310 I=2,NRP1
         K = I

         INCLUDE 'tent.inc'

         ZV2M = (CS(I)/CSM(I-1))**IEXV2
         ZV2P = (CS(I)/CSM(I  ))**IEXV2
      
         BsubM(kxW1+lxrow, kxW1+lxcol,k)=
     $     BsubM(kxW1+lxrow, kxW1+lxcol,k) + 
     $     FF(RGV1G11(i,msb),RGV1G11M(i-1,msb),RGV1G11M(i,msb))
         AsubM(kxW1+lxrow, kxW1+lxcol,k)=
     $     AsubM(kxW1+lxrow, kxW1+lxcol,k) + FFM(RGV1G11M(i-1,msb))
         CsubM(kxW1+lxrow, kxW1+lxcol,k)=
     $     CsubM(kxW1+lxrow, kxW1+lxcol,k) + FFP(RGV1G11M(i,msb))
         
         HsubM(kxW1+lxrow, kyW2+lycol,k)=
     $     HsubM(kxW1+lxrow, kyW2+lycol,k) + 
     $     FGM(RGV1G12(i,msb)*zv2m,RGV1G12M(i-1,msb))
         EsubM(kxW1+lxrow, kyW2+lycol,k)=
     $     EsubM(kxW1+lxrow, kyW2+lycol,k) + 
     $     FGP(RGV1G12(i,msb)*zv2p,RGV1G12M(i,msb))

         HsubM(kxW1+lxrow,kyPr+lycol,k) = znorm*RHO(I)*JACOBI(i,msb)
     $     + cma*RHO(I)*FGM(G12B2B2(i,msb), G12B2B2M(i-1,msb))
     $     + cna*RHO(I)*FGM(G12B2B3(i,msb), G12B2B3M(i-1,msb))
         EsubM(kxW1+lxrow,kyPr+lycol,k) =-znorm*RHO(I)*JACOBI(i,msb)
     $     + cma*RHO(I)*FGP(G12B2B2(i,msb), G12B2B2M(i,msb)) 
     $     + cna*RHO(I)*FGP(G12B2B3(i,msb), G12B2B3M(i,msb))
 310  CONTINUE
      ENDIF

C-----------------------------------------------------------------------
C.. EQUATION FOR KYW2, DEFINED ON HALD-INTEGER MESH
C-----------------------------------------------------------------------
      IF (KXW1.GT.0.AND.KYW2.GT.0) THEN
      IEXE = -1
      DO 320 I=1,NR
      K = I
      INCLUDE 'tophat.inc'

C.....EQUATION MULTIPLIED BY 1/S.

      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2

      FSUBM(KYW2+LYROW,KXW1+LXCOL,k) = FSUBM(KYW2+LYROW,KXW1+LXCOL,k)
     &     + GF(ZEM*RGV1G12(I,MSB),RGV1G12M(I,MSB))
      GSUBM(KYW2+LYROW,KXW1+LXCOL,k) = GSUBM(KYW2+LYROW,KXW1+LXCOL,k)
     &     + GF(ZEP*RGV1G12(I+1,MSB),RGV1G12M(I,MSB))
      DSUBM(KYW2+LYROW,KYW2+LYCOL,k) = DSUBM(KYW2+LYROW,KYW2+LYCOL,k)  
     &     + GG(RGV2G22M(I,MSB),ZEM*ZV2M*RGV2G22(I,MSB),
     &          ZEP*ZV2P*RGV2G22(I+1,MSB))
      DSUBM(KYW2+LYROW,KYPR+LYCOL,k) =
     &     -CMA*RHOM(I)*GG(TB2M(I,MSB),ZEM*TB2(I,MSB),ZEP*TB2(I+1,MSB))
     &     +CNA*RHOM(I)*GG(G22B2B2M(I,MSB),ZEM*G22B2B2(I,MSB),
     &                     ZEP*G22B2B2(I+1,MSB))
 320  CONTINUE
      ENDIF

C-----------------------------------------------------------------------
C.. EQUATION FOR KYW3, DEFINED ON HALD-INTEGER MESH
C-----------------------------------------------------------------------
      IF (KYW3.GT.0) THEN
      IEXE = -1
      DO 330 I=1,NR

C.....EQUATION MULTIPLIED BY 1/S.

      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE

      INCLUDE 'tophat.inc'

      DSUBM(KYW3+LYROW,KYW3+LYCOL,I) = 
     &     GG(RGV3G33M(I,MSB),ZEM*RGV3G33(I,MSB),ZEP*RGV3G33(I+1,MSB))
      DSUBM(KYW3+LYROW,KYPR+LYCOL,I)=
     &    -CNA*RHOM(I)*GG(B3j2M(I,MSB)*TM(I),ZEM*B3J2(I,MSB)*T(I),
     &                    ZEP*B3J2(I+1,MSB)*T(I+1))
 330  CONTINUE
      ENDIF

 160  CONTINUE
 
      IF (INCKIN.GT.0) 
     &CALL KCOEFFI(MROW,MSA,MSB,CMROW,CMA,CMB,CNA,SHIFTC,SHIFTM,
     &             ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
 
      RETURN
      END
*DECK PROFIL
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------DENSITY & RESISTIVITY PROFILES ---- A.B.- 05.04.93 -------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
      SUBROUTINE PROFIL
C     =================
C
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE GIJLM
      USE FEEDBACKM
      USE KINETICM
      USE REORBITM
      USE ANISOTROPICM
      USE MPIENV
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INTEGER   I,J,MS,K,NRH,IFAIL,IRAD,NRAD,MMIN,MMAX,MC
      REAL*8    ZRAD,ZEXPD,ZEXPT,ZPRESI,ZRES,ZETA0,ZETA1,Z,ZW0,PSHIFT
      REAL*8    ZROT,ZS,ZF,ZEPSIL,ZRTEPS,X,H1,H2,H3,H4
      REAL*8    PVCON,PVCEN,QMIN,QMAX,ZQE,ZQMI,ZQMU,ZQMU0,ZQEPS0,TMP
      PARAMETER (PVCON=0., PVCEN= 0.)
      INTEGER,DIMENSION(:),ALLOCATABLE::KC
      REAL*8,DIMENSION(:),ALLOCATABLE::RAD,RADM,SR,SRM,
     R          RPRS,RPRV,RPR2,RTMPR,RTMPI
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RPRVV
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::ROTCK,ROTCKM
      REAL*8::ZTEMP
      REAL*8    PARAROTA,DEL,PARANORM,SMI2ME

C     COPY <QPLS> TO <Q>
      DO I=1,NRP1
         Q(I)  = QPLS(I)
         QM(I) = QPLSM(I)
      ENDDO

C     RAD = SQRT(TOROIDAL FLUX)
      ALLOCATE(RAD(NRP1),RADM(NRP1),SR(NRP1),SRM(NRP1))
      CALL RADF2RADP(CS,CSM,Q,RAD,RADM,NRP1)

C     PRESSURE-PROFILE: ADD BASE PRESSURE P0 TO EQUILIBRIUM PRESSURE

      PEQ(1) = (4.*PEQM(1) - PEQ(2))/3.
      PSHIFT = (P0 - P0OLD) * REAL(PEQ(1))
      PEQM(NRP1)=0.
      DO 30 I=1,NRP1
      PEQ (I) = PEQ (I) + PSHIFT
 30   PEQM(I) = PEQM(I) + PSHIFT
      P0OLD = P0

C     EQUILIBRIUM CURRENT PROFILE: PLACE HOLDER

      DO I=2,NRP1
         JEQ(I) = 1./Q(I)
      ENDDO
      JEQ(1) = JEQ(2)

C     GAMMA-PROFILE
C     NPROFG = 0: CONSTANT GAMMA
C              1: GAUSSIAN PROFILE
C              2: SET GAMMA=GAMQ0 NEAR RATIONAL SURFACES
C              3: PEAKING NEAR PLASMA EDGE
      
      IF (NPROFG.EQ.0) THEN
         DO I = 1,NRP1
            GAMARR(I) = GAMMA
            GAMARM(I) = GAMMA
         ENDDO
      ELSEIF (NPROFG.EQ.1) THEN
         DO I = 1,NRP1
            GAMARR(I) = GAMMA*EXP(-((Q(I)-GAMQ0)/GAMWID)**4)
            GAMARM(I) = GAMMA*EXP(-((QM(I)-GAMQ0)/GAMWID)**4)
         ENDDO
      ELSEIF (NPROFG.EQ.2) THEN
         GAMARR = GAMMA
         GAMARM = GAMMA
         DO J=1,NRATSURF
            DO I=IRATSURF(J)-KGAM,IRATSURF(J)+KGAM
               IF (I.GE.1.AND.I.LE.NRP1) THEN
                  GAMARR(I) = GAMQ0
                  GAMARM(I) = GAMQ0
               ENDIF
            ENDDO
         ENDDO
      ELSEIF (NPROFG.EQ.3) THEN
         DO I=1,NRP1
            GAMARR(I) = GAMMA + GAMQ0*CS(I)**GAMWID
            GAMARM(I) = GAMMA + GAMQ0*CSM(I)**GAMWID
         ENDDO
      ELSE 
         STOP 'NPROFG'
      ENDIF

C     DENSITY-PROFILE
C     NPROFN = 0: POLYNOMIAL IN R**2
C              1: POLYNOMIAL IN PSI
C L.P.         11: POWER LAW IN PSI
C              2: POWER OF THE PRESSURE
C              3: POWER OF Q
C              4: EXPERIMENTAL PROFILE

      IF (NPROFN.EQ.0) THEN
         DO I=1,NRP1
            ZRAD=CSV(I)**2
            RHO(I) = 1. + ZRAD*(RHO1 + ZRAD*RHO2)
            ZRAD=CSVM(I)**2
            RHOM(I) = 1. + ZRAD*(RHO1 + ZRAD*RHO2)
            IF (RHO(I).LT.0.) STOP 'RHO'
         ENDDO
      ELSEIF (NPROFN.EQ.1) THEN
         DO I=1,NRP1
            ZRAD=CS(I)*CS(I)
            RHO(I) = 1. + ZRAD*(RHO1 + ZRAD*RHO2)
            ZRAD=CSM(I)*CSM(I)
            RHOM(I) = 1. + ZRAD*(RHO1 + ZRAD*RHO2)
            IF (RHO(I).LT.0.) STOP 'RHO'
         ENDDO
      ELSEIF (NPROFN.EQ.11) THEN
         DO I=1,NRP1
            ZRAD=CS(I)*CS(I)
            RHO(I) = (1. - ZRAD**RHO1)**RHO2
            ZRAD=CSM(I)*CSM(I)
            RHOM(I) = (1. - ZRAD**RHO1)**RHO2
            IF (RHO(I).LT.0) STOP 'RHO'
         ENDDO
      ELSEIF (NPROFN.EQ.2) THEN
         ZEXPD= 1./(1.+ETAI)
         DO I=1,NRP1
            ZPRESI  = PEQ(I) / PEQ(1)
            RHO(I)  = ZPRESI**ZEXPD
            ZPRESI  = PEQM(I) / PEQ(1)
            RHOM(I) = ZPRESI**ZEXPD
         ENDDO
      ELSEIF (NPROFN.EQ.3) THEN
         DO I=1,NRP1
            RHO(I)  = (Q(I)/Q(1))**DEXQ
            RHOM(I) = (QM(I)/Q(1))**DEXQ
         ENDDO
      ELSEIF (NPROFN.EQ.4) THEN
         OPEN(99,FILE='PROFDEN.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV = RPRV/ZTEMP
         IF (NEXPV.EQ.1) ZNE0 = ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFN=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(RHO,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(RHOM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         RHOM(NRP1) = RHO(NRP1)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFN'
      ENDIF

      TOTDENS = 0.0
      DO I=1,NR
         TOTDENS = TOTDENS + RHOM(I)*CSH(I)
      ENDDO
      TOTDENS = TOTDENS*PI*PI*4.0

C     CALCULATE PHYSICAL QUANTITIES AT MAGNETIC AXIS
C     ASSUMING D AS THE BULK ION SPECIES
C     USE NAMELIST VARIABLE FRACPTH = FRACTION OF THERMAL PRESSURE
C     OVER TOTAL EQUILIBRIUM PRESSURE AT MAGNETIC AXIS, TO RECOVER 
C     THE THERMAL TEMPERATURE AT THE MAGNETIC AXIS 
      ZQE   = 1.6021917e-19    !ELECTRON CHARGE [C]
      ZQMI  = 1.67261e-27      !PRONTON MASS    [KG]
      ZQMU  = ESPECIES_M(1)    !THERMAL ION MASS NUMBER
      ZQMU0 = 4.0E-7*PI        ![H/m]
      ZQEPS0= 8.8542E-12       ![F/m]

      IF (NEXPV.EQ.0) 
     &ZNE0     = (OMEGACI0/ZQE/R0EXP)**2*ZQMI*ZQMU/ZQMU0 !ELECTRON DENSITY [1/M^3] 
      IF (NEXPV.EQ.1)
     &OMEGACI0 = ZQE*R0EXP*SQRT(ZNE0/ZQMI/ZQMU*ZQMU0)
      ZTAUA0   = R0EXP*SQRT(ZQMU0*ZQMI*ZQMU*ZNE0)/B0EXP  !ALFVEN TIME [SECOND]

C     NORMALISED DENSITY AND PRESSURE PROFILES FOR HOT IONS
C     ASSUMING GENERALLY MULTIPLE HOT ION SPECIES 
C     NOTE THAT DUMMY DATA ARE READ IN FOR IF0TYPE=2(ISOTROPIC ALPHA-PARTICLE MODEL)
C     THE DENSITY AND PRESSURE FRACTIONS FOR IF0TYPE=2 WILL BE 
C     RECOMPUTED IN <KNPFRACF02>
C     NPROFK = 0: ANALYTIC: CONSTANT
C              1: ANALYTIC PROFILE 1
C              4: EXPERIMENTAL PROFILE WITH AMPLITUDE
C
      IF (.NOT.ALLOCATED(ESPECIES_DENF)) 
     &ALLOCATE( ESPECIES_DENF(NRP1,2,NSPECIES),
     &          ESPECIES_PREF(NRP1,2,NSPECIES),
     &          ESPECIES_DEN (NRP1,2,NSPECIES),
     &          ESPECIES_PRE (NRP1,2,NSPECIES),
     &          ESPECIES_PREP(NRP1,2,NSPECIES),
     &          ESPECIES_TEM (NRP1,2,NSPECIES),
     &          ESPECIES_REL (NLAMK,2), 
     &          ESPECIES_REE (NEPK,2) )

      DO I=1,NLAMK
         ESPECIES_REL(I,1) = DFLOAT(I-1)/DFLOAT(NLAMK-1)
      ENDDO
      DO I=1,NEPK
         ESPECIES_REE(I,1) = DFLOAT(I-1)/DFLOAT(NEPK-1)
      ENDDO
      ESPECIES_REL(:,2) = 0.0
      ESPECIES_REE(:,2) = 0.0

      IF (NSPECIES.GT.2) THEN
      IF (NPROFK.EQ.0) THEN
         ESPECIES_DENF(:,:,3:NSPECIES)  = EP_DENF0
         ESPECIES_PREF(:,:,3:NSPECIES)  = EP_PREF0
         ESPECIES_REL(:,2) = 1.0
         ESPECIES_REE(:,2) = 1.0
      ELSEIF (NPROFK.EQ.1) THEN
         DO K=3,NSPECIES
            ESPECIES_DENF(:,1,K)  = EP_DENF0*(1.-CS(1:NRP1)**2)
            ESPECIES_DENF(:,2,K)  = EP_DENF0*(1.-CSM(1:NRP1)**2)
            ESPECIES_PREF(:,1,K)  = EP_PREF0*(1.-CS(1:NRP1)**4)
            ESPECIES_PREF(:,2,K)  = EP_PREF0*(1.-CSM(1:NRP1)**4)
         ENDDO
         DO I=1,NLAMK
            ZTEMP = DFLOAT(I-1)/DFLOAT(NLAMK-1)
            ESPECIES_REL(I,2) = 1.0 - ZTEMP**2
         ENDDO
         DO I=1,NEPK
            ZTEMP = DFLOAT(I-1)/DFLOAT(NEPK-1)
            ESPECIES_REE(I,2) = ZTEMP*(1.0-ZTEMP)*2.0
         ENDDO
      ELSEIF (NPROFK.EQ.4) THEN
         OPEN(99,FILE='PROFPA.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRVV(NRAD,3:NSPECIES),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),(RPRVV(I,K),K=3,NSPECIES)
         ENDDO
         CLOSE(99)
         IF (NEXPV.EQ.1) FRACPTH = 1.0/(1.0+SUM(RPRVV(1,:)))
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFK=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         DO K=3,NSPECIES
            CALL SPLINE1D(DAK,SR,NRP1,RPRVV(:,K),RPRS,NRAD,RPR2)
            CALL SPLINE1D(DAKM,SRM,NRP1-1,RPRVV(:,K),RPRS,NRAD,RPR2)
            ESPECIES_PREF(:,1,K) = DAK *ALPHAH
            ESPECIES_PREF(:,2,K) = DAKM*ALPHAH
         ENDDO
         DEALLOCATE(RPRS,RPRVV,RPR2)

         OPEN(99,FILE='PROFDA.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRVV(NRAD,3:NSPECIES),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),(RPRVV(I,K),K=3,NSPECIES)
         ENDDO
         CLOSE(99)
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFK=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         DO K=3,NSPECIES
            CALL SPLINE1D(DAK,SR,NRP1,RPRVV(:,K),RPRS,NRAD,RPR2)
            CALL SPLINE1D(DAKM,SRM,NRP1-1,RPRVV(:,K),RPRS,NRAD,RPR2)
            ESPECIES_DENF(:,1,K) = DAK
            ESPECIES_DENF(:,2,K) = DAKM
         ENDDO
         DEALLOCATE(RPRS,RPRVV,RPR2)

         NRH = 0
         DO K=3,NSPECIES
            IF (ISPECIES_F0(K).EQ.5) NRH = 1
         ENDDO
         IF (NRH.EQ.1) THEN
         OPEN(99,FILE='PROFREL.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = MAXVAL(RPRV)
         RPRV = RPRV/ZTEMP
         CALL SPLINE1D(ESPECIES_REL(:,2),ESPECIES_REL(:,1),NLAMK,
     &                 RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)

         OPEN(99,FILE='PROFREE.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = MAXVAL(RPRV)
         RPRV = RPRV/ZTEMP
         CALL SPLINE1D(ESPECIES_REE(:,2),ESPECIES_REE(:,1),NEPK,
     &                 RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
         ENDIF
      ELSE
         STOP 'NPROFK'
      ENDIF
      ENDIF

C     NORMALISED TEMPERATURE PROFILES FOR THERMAL ION & ELECTRONS
C     NOTE THAT TEMPTI/E-ARRAYS ARE ALWAYS NORMALISED TO 1 AT MAGNETIC AXIS
C     NPROFIE = 0: PRESSURE/DENSITY
C               4: EXPERIMENTAL PROFILE WITH AMPLITUDE IN [eV]
      
      IF (NPROFIE.EQ.0) THEN
         ZTEMP = PEQ(1)/RHO(1)
         IF (ZTEMP.EQ.0.0) ZTEMP=1.0
         DO I=1,NRP1
            TEMPI(I)  = PEQ(I)/RHO(I)/ZTEMP
            TEMPE(I)  = TEMPI(I)
         ENDDO
         DO I=1,NR
            TEMPIM(I) = PEQM(I)/RHOM(I)/ZTEMP
            TEMPEM(I) = TEMPIM(I)
         ENDDO
         IF (ABS(RHO(NRP1)).LE.1.0E-8) THEN
            TEMPI(NRP1) = TEMPI(NR)
            TEMPE(NRP1) = TEMPE(NR)
         ENDIF
      ELSEIF (NPROFIE.EQ.4) THEN
         OPEN(99,FILE='PROFTI.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (NEXPV.EQ.1) ZTI0=ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFIE=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(TEMPI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(TEMPIM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)

         OPEN(99,FILE='PROFTE.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (NEXPV.EQ.1) ZTE0=ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFIE=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(TEMPE,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(TEMPEM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFIE'
      ENDIF
      TEMPIM(NRP1) = TEMPI(NRP1)
      TEMPEM(NRP1) = TEMPE(NRP1)

      IF (NEXPV.EQ.1) ALPHAP = ZTI0/(ZTI0+ZTE0)
      IF (NEXPV.EQ.1.AND.NPROFK.NE.4) 
     &   FRACPTH = (ZTI0+ZTE0)/PEQ(1)*ZNE0*ZQE/B0EXP**2*ZQMU0
      IF (NEXPV.EQ.0) THEN
      TMP      = FRACPTH*PEQ(1)/ZNE0/ZQE*B0EXP**2/ZQMU0
      ZTI0     = ALPHAP*TMP                              !BULK ION TEMPERATURE [EV]
      ZTE0     = (1.-ALPHAP)*TMP                         !BULK ELECTRON TEMPERATURE [EV] 
      ENDIF
      
C     RTYPE4=T_PERP/T_PARA FOR IF0TYPE=4 EPS MODEL
C     NPROFR4 = 0: CONSTANT
C               4: EXPERIMENTAL PROFILE
C
      IF (.NOT.ALLOCATED(RTYPE4)) 
     &    ALLOCATE(RTYPE4(NRP1,2))
C
      RTYPE4 = 1.
      IF (NPROFR4.EQ.0) THEN
         RTYPE4 = R0TYPE4
      ELSEIF (NPROFR4.EQ.4) THEN
         OPEN(99,FILE='PROFR4.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (NEXPV.EQ.1) R0TYPE4 = ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFR4=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(ROT,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(ROTM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         ROTM(NRP1) = ROT(NRP1)
         RTYPE4(:,1) = ROT *R0TYPE4
         RTYPE4(:,2) = ROTM*R0TYPE4
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFR4'
      ENDIF

C     CHECK: RTYPE4>0
      IF (MINVAL(RTYPE4(:,1)).LE.0.) STOP 'RTYPE4'

C     STYPE4=ASYMMETRIC FACTOR FOR IF0TYPE=4 PASSING EPS MODEL
C     NPROFS4 = 0: CONSTANT
C               4: EXPERIMENTAL PROFILE
C
      IF (.NOT.ALLOCATED(STYPE4)) 
     &   ALLOCATE(STYPE4(NRP1,2))
C
      STYPE4 = 0.
      IF (NPROFS4.EQ.0) THEN
         STYPE4 = S0TYPE4
      ELSEIF (NPROFS4.EQ.4) THEN
         OPEN(99,FILE='PROFS4.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (NEXPV.EQ.1) S0TYPE4 = ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFS4=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(ROT,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(ROTM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         ROTM(NRP1) = ROT(NRP1)
         STYPE4(:,1) = ROT *S0TYPE4
         STYPE4(:,2) = ROTM*S0TYPE4
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFS4'
      ENDIF

C     CHECK: |STYPE4|<1
      IF (MAXVAL(ABS(STYPE4(:,1))).GT.1.) STOP 'STYPE4'
      
C     READ IN TIME TRACE OF COIL CURRENT FROM <COILCURR.IN>
      IF (KCOILCURR.EQ.4) THEN
         OPEN(99,FILE='COILCURR.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NCOILCURR
         ALLOCATE(COILCURR(NCOILCURR,2),COILCURR2(NCOILCURR))
         DO I=1,NCOILCURR
            READ(99,*) COILCURR(I,1),COILCURR(I,2)
         ENDDO
         CLOSE(99)
      ENDIF

      IF (KCOILCURR.EQ.5) THEN
         OPEN(99,FILE='COILCURR.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NCOILCURR
         ALLOCATE(COILCURR(NCOILCURR,1+2*NCOIL),COILCURR2(NCOILCURR))
         DO I=1,NCOILCURR
            READ(99,*) COILCURR(I,1),(COILCURR(I,2*J),
     &                                COILCURR(I,2*J+1),J=1,NCOIL)
         ENDDO
         CLOSE(99)
      ENDIF

C     COMPUTE DENSITY AND PRESSURE FRACTIONS FOR IF0TYPE=2
C     THE FOLLOWING TWO CALLS ARE NEEDED HERE TO COMPUTE THERMAL
C     PRESSURE, BY SUBTRACTING HOT ION PRESURES FROM THE TOTAL 
C     EQUILIBRIUM PRESSURE
      CALL KNPFRACF02

C     GET (SURFACE AVERAGED) EQUILIBRIUM RADIAL PROFILES OF DENSITY,
C     PRESSURE, AND TEMPERATURE FOR EACH PARTICLE SPECIES
      CALL KEQPROF

      IF (NEXPV.EQ.0.AND.NSPECIES.GE.3) THEN
      FRACPTH  = (ESPECIES_PRE(1,1,1)+ESPECIES_PRE(1,1,2))/PEQ(1)
      TMP      = FRACPTH*PEQ(1)/ZNE0/ZQE*B0EXP**2/ZQMU0
      ZTI0     = ALPHAP*TMP                              !BULK ION TEMPERATURE [EV]
      ZTE0     = (1.-ALPHAP)*TMP                         !BULK ELECTRON TEMPERATURE [EV] 
      ENDIF

C     ZEFF-PROFILE
C     NPROFZ = 0: CONSTANT        
C              4: EXPERIMENTAL PROFILE 
      
      IF (NPROFZ.EQ.0) THEN
         DO I=1,NRP1
            ZEFFI(I) = ZEFF0
            ZEFFM(I) = ZEFF0
         ENDDO
      ELSEIF (NPROFZ.EQ.4) THEN
         OPEN(99,FILE='PROFZEFF.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFZ=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(ZEFFI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(ZEFFM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFZ'
      ENDIF

C     RESISTIVITY-PROFILE
C     NPROFT = 0: POLYNOMIAL IN R**2
C            = 1: POLYNOMIAL IN PSI
C            = 2: POWER OF THE PRESSURE
C            = 3: FROM ELECTRON TEMPERATURE
C            = 5: PIECE-WISE LINEAR MODEL  

      ZETA0 = ETA / ASPCT ** 2
      ZETA1 = ETAXIS / ASPCT ** 2
      IF (NPROFT.EQ.0) THEN
         DO I=1,NRP1
            ZRAD=CSV(I)**2
            RESIST(I) = ZETA0/(1. + ZRAD*(SIGMA1 + ZRAD*SIGMA2))
            IF (RESIST(I).LT.0.) STOP 'RESIST'
         ENDDO
      ELSEIF (NPROFT.EQ.1) THEN
         DO I=1,NRP1
            ZRAD=CS(I)*CS(I)
            RESIST(I) = ZETA0/(1. + ZRAD*(SIGMA1 + ZRAD*SIGMA2))
            IF (RESIST(I).LT.0.) STOP 'RESIST'
         ENDDO
      ELSEIF (NPROFT.EQ.2) THEN
         ZEXPT=-1.5*ETAI/(1.+ETAI)
         DO I=1,NRP1
            ZPRESI    = PEQ(I) / PEQ(1)
            RESIST(I) = ZPRESI**ZEXPT * ZETA0
         ENDDO
      ELSEIF (NPROFT.EQ.3) THEN
         DO I=1,NRP1
            IF (TEMPE(I).NE.0.) THEN
               RESIST(I) = ZETA0*TEMPE(I)**(-1.5)
            ELSE
               WRITE(*,*) 'WARNING: TEMPE=0'
               RESIST(I) = ZETA0*TEMPE(I-1)**(-1.5)
            ENDIF
         ENDDO
      ELSEIF (NPROFT.EQ.5) THEN
         RESIST = 1.0
         IF (NRATSURF.GE.2) THEN
            ZS = CS(IRATSURF(1))
            ZF = CS(IRATSURF(2))
            X  = MIN(SIGMA2,(ZF-ZS)*0.5)
            H1 = (ZS+ZF-X)*0.5
            H2 = (ZS+ZF+X)*0.5
            H3 = (1.0-SIGMA1)/(H1-H2)
            H4 = 1.0-H3*H1
            DO I=1,NRP1
               IF (CS(I).GE.H2) RESIST(I)=SIGMA1
               IF (CS(I).GT.H1.AND.CS(I).LT.H2) RESIST(I)=H3*CS(I)+H4
            ENDDO
            WRITE(*,*) 'RESISTIVITY: ',SIGMA1,X
         ENDIF
         RESIST = ZETA0*RESIST
      ELSE
         STOP 'NPROFT'
      ENDIF

C     MODIFY RESISTIVITY
C     Y.Q.LIU, 2009-08-20

C     ADD ZEFF FACTOR
      RESIST = RESIST*ZEFFI

C     FIX POINTS NEAR AXIS AND PLASMA EDGE
      IF (NFIT.LT.NRES) THEN
      H1 = CS(NFIT)
      H2 = CS(NRES)
      RESIST(1:NFIT) = ZETA1
      DO I=NFIT+1,NRES-1
C        RESIST(I) = ZETA1+(RESIST(NRES)-ZETA1)*
C    &               ((H2-H1)**2-(CS(I)-H2)**2)**2/(H2-H1)**4
         RESIST(I) = ZETA1
      ENDDO
      ENDIF
      IF (NFIT.LT.NRESR) THEN
      H3 = CS(NRP1-NRESR+1)
      H4 = CS(NRP1-NFIT+1)
      RESIST(NRP1-NFIT+1:NRP1) = ZETA1
      DO I=NFIT+1,NRESR-1
         RESIST(NRP1-I+1) = ZETA1+(RESIST(NRP1-NRESR+1)-ZETA1)*
     &               ((H4-H3)**2-(CS(NRP1-I+1)-H3)**2)**2/(H4-H3)**4
      ENDDO
      ENDIF

      IF (NRESR.LT.-1.AND.NFIT.LT.ABS(NRESR)) 
     &   RESIST(NRP1-ABS(NRESR)+1:NRP1) = 0.0

C     FURTHER MODIFICATION OF RESISTIVITY PROFILE NEAR PLASMA EDGE
C     YQLIU, 2009-10-01
      IF (NRESR.EQ.-1) THEN
      DO I=1,NRP1
         IF (RESIST(I).GT.ZETA0*100.0) RESIST(I) = ZETA0*100.0
      ENDDO
      ENDIF

      DO I=1,NR
         RESISM(I) = (RESIST(I)+RESIST(I+1))*0.5
      ENDDO
      RESISM(NRP1) = RESIST(NRP1)
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      WRITE(*,*) 'RESISTIVITY: ',ETA,ZETA0,ASPCT
      ENDIF
      
      IF (.NOT.(0..LE.ETACS1.AND.ETACS1.LT.ETACS2.AND.ETACS2.LT.1.))
     &   GOTO 550
      WRITE(*,*) ' WARNING:  NONCONSTANT RESISTIVITY, ETACS1 =',
     &           ETACS1,' ETACS2 =',ETACS2
      DO 540 I = 1,NRP1
      IF (CS(I).LT.ETACS1) RESIST(I) = 0.
      IF (.NOT.(ETACS1.LE.CS(I).AND.CS(I).LT.ETACS2))
     &   GOTO 540
      Z = (CS(I)-ETACS1)/(ETACS2 - ETACS1)
      RESIST(I) = RESIST(I) * Z*Z*(3.-2.*Z)
 540  CONTINUE
 550  CONTINUE
      IF (.NOT.(0..LT.ETACS3.AND.ETACS3.LT.ETACS4.AND.ETACS4.LE.1.))
     &   GOTO 570
      WRITE(*,*) ' WARNING:  NONCONSTANT RESISTIVITY, ETACS3 =',
     &           ETACS3,' ETACS4 =',ETACS4
      DO 560 I = 1,NRP1
      IF (CS(I).GT.ETACS4) RESIST(I) = 0.
      IF (.NOT.(ETACS3.LE.CS(I).AND.CS(I).LT.ETACS4))
     &   GOTO 560
      Z = (CS(I)-ETACS4)/(ETACS3 - ETACS4)
      RESIST(I) = RESIST(I) * Z*Z*(3.-2.*Z)
 560  CONTINUE
 570  CONTINUE


C     ROTATION-PROFILE
C     NPROFR = 0: UNIFROM ROTATION
C              1: ANALYTIC MODEL FOR DIII-D
C              2: ANALYTIC MODEL FOR JET 
C              3: KEEP ROTATION NEAR RATIONAL SURFACE ONLY
C              4: EXPERIMENTAL PROFILE   
C L.P.         5: ANOTHER ANALYTIC MODEL

      IF (NPROFR.EQ.0) THEN
         DO I=1,NRP1
            ROT(I)  = ROTE
            ROTM(I) = ROTE
         ENDDO
      ELSEIF (NPROFR.EQ.1) THEN
         DO I = 1,NRP1
            X = CS(I)/CSROT0
            IF (X .GT. 1.) ROT(I) = 0.
            IF (X .GE. 0. .AND. X .LE. 1.) 
     &           ROT(I) = (OMEGA0-OMEGA1)*(1.-X*X*(2.-X*X))
            ROT(I) = ROTE*(OMEGA1+ROT(I))
         ENDDO
         DO I = 1,NR
            X = CSM(I)/CSROT0
            IF (X .GT. 1.) ROTM(I) = 0.
            IF (X .GE. 0. .AND. X .LE. 1.) 
     &           ROTM(I) = (OMEGA0-OMEGA1)*(1.-X*X*(2.-X*X))
            ROTM(I) = ROTE*(OMEGA1+ROTM(I))
         ENDDO
      ELSEIF (NPROFR.EQ.2) THEN
         DO I = 1,NRP1
            X = CSV(I)/CSROT0
            ROT(I) = 0.
            IF (X .GE. 0. .AND. X .LT. 1.) 
     &           ROT(I) = (OMEGA0-OMEGA1)*(1.-X**ALPHAR)**BETAR
            ROT(I) = ROTE*(OMEGA1+ROT(I))
         ENDDO
         DO I = 1,NR
            X = CSVM(I)/CSROT0
            ROTM(I) = 0.
            IF (X .GE. 0. .AND. X .LT. 1.) 
     &           ROTM(I) = (OMEGA0-OMEGA1)*(1.-X**ALPHAR)**BETAR
            ROTM(I) = ROTE*(OMEGA1+ROTM(I))
         ENDDO
      ELSEIF (NPROFR.EQ.3) THEN
         ROT  = 0.0    
         ROTM = 0.0    
         DO J=1,NRATSURF
            DO I=IRATSURF(J)-KGAM,IRATSURF(J)+KGAM
               IF (I.GE.1.AND.I.LE.NRP1) THEN
                  ROT(I)  = ROTE  
                  ROTM(I) = ROTE  
               ENDIF
            ENDDO
         ENDDO
      ELSEIF (NPROFR.EQ.4) THEN
         OPEN(99,FILE='PROFROT.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
C        ZTEMP = MAXVAL(ABS(RPRV))
         ZTEMP = RPRV(1)
         IF (NEXPV.EQ.1) ROTE = ZTEMP*ZTAUA0
         RPRV = ROTE*RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFR=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(ROT,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(ROTM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
C  L. PIGATTO =================================================
C  EXTRA ANALYTICAL ROTATION PROFILE USING OMEGA0 AND OMEGA1
      ELSEIF (NPROFR.EQ.5) THEN
         DO I = 1,NRP1
            X = CSV(I)/CSROT0
            ROT(I) = 0.
            IF (X .GE. 0. .AND. X .LT. 1.) 
     &           ROT(I) = (OMEGA0-OMEGA1)*(1.-X*X)
            ROT(I) = ROTE*(OMEGA1+ROT(I))
         ENDDO
         DO I = 1,NR
            X = CSVM(I)/CSROT0
            ROTM(I) = 0.
            IF (X .GE. 0. .AND. X .LT. 1.) 
     &           ROTM(I) = (OMEGA0-OMEGA1)*(1.-X*X)
            ROTM(I) = ROTE*(OMEGA1+ROTM(I))
         ENDDO
      ELSE
         STOP 'NPROFR'
      ENDIF

C  BY LI LI =====================================================
C  LOCAL MODIFICATION OF ROTATION PROFILE NEAR THE FIRST RATIONAL
C  SURFACE. ONLY WORKS FOR OPTION NPROFR=4. USE PARAMETER OMEGA0 TO
C  MIDIFY ROTATION AMPLITUDE AT RATIONAL SURFACE, SUCH THAT THE FINAL
C  ROTATION FREQUENCY AT THE RATIONAL SURFACE IS EQUAL TO ROTE*OMEGA0
C  NFLAGROTA= 0  NO MODIFICATION
C           = 1  METHOD I
C           = 2  METHOD II
C           = 3  METHOD III
              

      IF (NFLAGROTA.GT.0.AND.NPROFR.EQ.4.AND.NRATSURF.GT.0) THEN
         PARANORM = ROT(1)
         DO I = 1,NRP1
            ROT(I)  = ROT(I)/PARANORM
            ROTM(I) = ROTM(I)/PARANORM
         ENDDO
    
         PARAROTA = ROT(IRATSURF(1))
         DO I=1,NRP1
            ROT(I)  = ROT(I) - PARAROTA
            ROTM(I) = ROTM(I) - PARAROTA
         ENDDO

C  METHOD I: ADD OMEGA0 AT RATIONAL POINT
       IF (NFLAGROTA.EQ.1) THEN
         DO I=1,NRP1
            ROT(I)  = ROT(I) + OMEGA0
            ROTM(I) = ROTM(I) + OMEGA0
         ENDDO         

C  METHOD II: EXP NEAR RATIONAL POINT
       ELSEIF (NFLAGROTA.EQ.2) THEN
         DO I = 1,NRP1
            DEL = 1.0 - EXP(-((Q(I)-2.0)/1.0E-3)**2.)
            ROT(I)  = ROT(I)*DEL + OMEGA0
            ROTM(I) = ROTM(I)*DEL + OMEGA0
         ENDDO

C  METHOD III: EXP NEAR RATIONAL POINT, ABS(ROT)         
       ELSEIF (NFLAGROTA.EQ.3) THEN
         DO I = 1,NRP1
            DEL = 1.0 - EXP(-((Q(I)-2.0)/1.0E-3)**2.)
            ROT(I)  = ABS(ROT(I))*DEL + OMEGA0
            ROTM(I) = ABS(ROTM(I))*DEL + OMEGA0
         ENDDO

       ENDIF

       DO I = 1,NRP1
          ROT(I)  = ROT(I)*PARANORM
          ROTM(I) = ROTM(I)*PARANORM
       ENDDO 
       ENDIF    
C  END OF MODIFYING ROTATION PROFILE =============

      IF (INCKIN.GT.0.AND.IPERTURB.EQ.2) THEN
        DO I=1,NRP1
           ROT(I) =0.0
           ROTM(I)=0.0
        ENDDO
      ENDIF


C     ROTATION-SHEAR
      DO I = 1,NR
        DROTM(I) = (ROT(I+1)-ROT(I))/(CS(I+1) - CS(I))
      END DO
      DO I = 2,NR
         H1 = (CS(I)-CS(I-1))/2
         H2 = (CS(I+1)-CS(I))/2
         DROT(I) = (H1/H2*ROTM(I)-H2/H1*ROTM(I-1))/(H1+H2)-
     &             (H1-H2)*ROT(I)/H1/H2
      END DO
C     DROT(1) = 2*DROTM(1) - DROT(2)
C     DROT(NRP1) = 2*DROTM(NR) - DROT(NR)
      DROT(1) = 0.0
      DROT(NRP1) = (ROT(NRP1)-ROTM(NR))/(CS(NRP1)-CSM(NR))
      ROTM(NRP1) = 0.
      DROTM(NRP1) = 0.

C GLX  POLOIDAL ROTATION-PROFILE
C      NPROFRP = 0: NO POLOIDAL ROTATION
C                1: UNIFORM POLOIDAL ROTATION
C                4: EXPERIMENTAL PROFILE   

      ROTP  = 0.0
      ROTPM = 0.0
      IF (NPROFRP.EQ.1 .OR. NPROFRC.EQ.21) THEN
         DO I=1,NRP1
            ROTP(I)  = ROTEP
            ROTPM(I) = ROTEP
         ENDDO
      ELSEIF (NPROFRP.EQ.4 .OR. NPROFRC.EQ.24) THEN
         OPEN(99,FILE='PROFROTP.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)

         ZTEMP = RPRV(1)
         IF (NEXPV.EQ.1) ROTEP = ZTEMP*ZTAUA0
         RPRV = ROTEP*RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFRP=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(ROTP,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(ROTPM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ENDIF

C GLX   RHO^-1*U
        DO I=1,NRP1
          RHOU(I) =1.0/RHO(I)*ROTP(I)
        ENDDO 
        DO I=1,NR
          RHOUM(I)=1.0/RHOM(I)*ROTPM(I)
        ENDDO 

C      D(RHO^-1*U)/DS 
         DO I=2,NR
            H1 = (CS(I)-CS(I-1))/2
            H2 = (CS(I+1)-CS(I))/2
            DRHOU(I)=(H1/H2*RHOUM(I)-H2/H1*RHOUM(I-1))/(H1+H2) -
     &               (H1-H2)*RHOU(I)/H1/H2
         ENDDO
         I  = NRP1
         H1 = CS(I)-CS(I-1)
         DRHOU(I)=(RHOU(I-1)+3*RHOU(I)-4*RHOUM(I-1))/H1
         DRHOU(1)=DRHOU(2)

         DO I=1,NR
            H1 = CS(I+1) - CS(I)
            DRHOUM(I)=(RHOU(I+1)-RHOU(I))/H1
         ENDDO

C     ROTATION-PROFILE FOR ROTC
C     NPROFRC = 0: NO INCLUSION OF ROTC
C               1: ANALYTIC: CONSTANT FOR M=1 HARMONIC
C               21:ANALYTIC: ROTEC*U*F/(RHO*R*R)
C                  WITH UNIFORM PARALLEL ROTATION
C               24:ANALYTIC: ROTEC*U*F/(RHO*R*R)
C                  WITH EXPERIMENTAL PARALLEL ROTATION
C               4: EXPERIMENTAL PROFILE   
C     ROTEC = ROTATION AMPLPLITUDE MULTIPLIER
C             GENERALLY BEING A COMPLEX NUMBER
      IF (NPROFRC.EQ.1) THEN
         MC = 1
         ALLOCATE(KC(MC))
         ALLOCATE(ROTCK(NRP1,MC),ROTCKM(NRP1,MC))
         KC(1)  = 1
         ROTCK  = ROTEC
         ROTCKM = ROTEC
      ELSEIF (NPROFRC.EQ.4) THEN
         OPEN(99,FILE='PROFROTC.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD,MC
         ALLOCATE(KC(MC))
         ALLOCATE(RPRS(NRAD),RPRVV(NRAD,2*MC),RPR2(NRAD))
         ALLOCATE(ROTCK(NRP1,MC),ROTCKM(NRP1,MC))
         ALLOCATE(RTMPR(NRP1),RTMPI(NRP1))

         READ(99,*) (KC(K),K=1,MC)
         DO I=1,NRAD
            READ(99,*) RPRS(I),(RPRVV(I,K),K=1,2*MC)
         ENDDO
         CLOSE(99)
         RPRVV = ROTEC*RPRVV
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFRC=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         DO K=1,MC
            CALL SPLINE1D(RTMPR,SR,NRP1,RPRVV(:,2*K-1),RPRS,NRAD,RPR2)
            CALL SPLINE1D(RTMPI,SR,NRP1,RPRVV(:,2*K),RPRS,NRAD,RPR2)
            ROTCK(:,K) = RTMPR + RTMPI*CI
            CALL SPLINE1D(RTMPR,SRM,NR,RPRVV(:,2*K-1),RPRS,NRAD,RPR2)
            CALL SPLINE1D(RTMPI,SRM,NR,RPRVV(:,2*K),RPRS,NRAD,RPR2)
            ROTCKM(:,K) = RTMPR + RTMPI*CI
         ENDDO
         DEALLOCATE(RPRS,RPRVV,RPR2,RTMPR,RTMPI)
      ENDIF

C     COMPUTE ROTC
      IF(NPROFRC.EQ.1 .OR. NPROFRC.EQ.4) THEN
      ROTC  = 0.
      ROTCM = 0.
      DO J=1,NCHI
      ZTEMP = 2.*PI/DFLOAT(NCHI)*(J-1)
      DO K=1,MC
         ROTC(:,J) = ROTC(:,J) +2.*REAL(ROTCK(:,K) *EXP(CI*ZTEMP*KC(K)))
         ROTCM(:,J)= ROTCM(:,J)+2.*REAL(ROTCKM(:,K)*EXP(CI*ZTEMP*KC(K)))
      ENDDO
      ENDDO
      DEALLOCATE(ROTCK,ROTCKM,KC)
      ENDIF

C GLX   COMPUTE PURE POLOIDAL ROTATION
C       NEED TO SET ROTEC=(-1.0,0.0) 
      IF(NPROFRC.EQ.21 .OR. NPROFRC.EQ.24) THEN
      ROTC  = 0.
      ROTCM = 0.
      DO J=1,NCHI
      DO I=1,NRP1
         ROTC(I,J) = ROTEC*RHOU(I)*T(I)/REQ(I,J)**2
      ENDDO
      DO I=1,NR
         ROTCM(I,J)= ROTEC*RHOUM(I)*TM(I)/REQM(I,J)**2
      ENDDO
      ENDDO
      ENDIF

C     VISCOSITY-PROFILE: NEOCLASSICAL 

      DO 580 I=2,NRP1
      ZEPSIL = CS(I)/ ASPCT
      ZRTEPS = SQRT(ZEPSIL)
      GMUNU(I)=2.*RHO(I)*Q(I)**2/CS(I)**2*ASPCT**2*
     &ZRTEPS * (1.46 -0.46 *ZEPSIL)/(1.-ZRTEPS*(1.46-0.46*ZEPSIL))
     & * PVCON*PVCEN*RHO(I)**2.5/REAL(PEQ(I)+1.e-8)**1.5
      GMUNU(I)=GMUNU(I)/9.
  580 CONTINUE
      GMUNU(1) = GMUNU(2)*2.
      DO I=1,NR
         GMUNUM(I) = (GMUNU(I)+GMUNU(I+1))*0.5
      ENDDO
      GMUNUM(NRP1) = 0.0

C     EFFECTIVE COLLISIONALITY-PROFILE FOR IONS 
C     NPROFUI = 0: ANALYTICAL: CONSTANT      
C               1: ANALYTICAL FORMULA [DNESTROVSKIJ, SUN_NF11, MARS-K]
C               2: ANALYTICAL FORMULA [NRL, MISK]
C               3: ANALYTICAL FORMULA [CLASSICAL]
C               4: EXPERIMENTAL PROFILE 
C     FOR NPROFUI=1,2,3:
C               SET NUEFFIA=1.0 IF DON'T WISH TO SCALE COLLISIONALITY
      
      TMP = 1.0
      IF (NPROFUI.EQ.1) TMP = SQRT(2.0)*PI
      IF (NPROFUI.EQ.2) TMP = 4.0*SQRT(PI)/3.0
      IF (NPROFUI.EQ.3) TMP = PI/SQRT(2.0)

      IF (NPROFUI.EQ.0) THEN
      WRITE(*,*) 'UNIFORM ION COLL., NPROFUI=',NPROFUI
         DO I=1,NRP1
            GNUI(I)  = NUEFFIA
            GNUIM(I) = NUEFFIA
         ENDDO
      ELSEIF (NPROFUI.GE.1.AND.NPROFUI.LE.3) THEN
         WRITE(*,*) 'ANALYTICAL ION COLL., NPROFUI=',NPROFUI
         ZTEMP = TMP*ZQE**2.5/(4.0*PI*ZQEPS0)**2/SQRT(ZQMI)
         ZTEMP = ZTEMP*ESPECIES_Z(1)**4/SQRT(ESPECIES_M(1))
         ZTEMP = ZTEMP*ZNE0/ZTI0**1.5
         ZTEMP = ZTEMP*ZTAUA0*NUEFFIA
         DO I=1,NRP1
            GNUI(I) = ZTEMP*RHO(I)/TEMPI(I)**1.5*
     &                (17.3-0.5*LOG(ZNE0*RHO(I)/1.0E+20)
     &                     +1.5*LOG(ZTI0*TEMPI(I)/1.0E+3))   
         ENDDO
         DO I=1,NR
            GNUIM(I) = ZTEMP*RHOM(I)/TEMPIM(I)**1.5*
     &                 (17.3-0.5*LOG(ZNE0*RHOM(I)/1.0E+20)
     &                      +1.5*LOG(ZTI0*TEMPIM(I)/1.0E+3))   
         ENDDO
         GNUIM(NRP1) = GNUIM(NR)
      ELSEIF (NPROFUI.EQ.4) THEN
         WRITE(*,*) 'READING ION COLL. PROFILE, NPROFUI=',NPROFUI
         OPEN(99,FILE='PROFNUI.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP*NUEFFIA
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFUI=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(GNUI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(GNUIM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFUI'
      ENDIF

C     EFFECTIVE COLLISIONALITY-PROFILE FOR ELECTRONS 
C     NPROFUE = 0: ANALYTICAL      
C               1: ANALYTICAL FORMULA [DNESTROVSKIJ, SUN_NF11]
C               2: ANALYTICAL FORMULA [NRL, MISK]
C               3: ANALYTICAL FORMULA [CLASSICAL]
C               4: EXPERIMENTAL PROFILE 
C     FOR NPROFUE=1,2,3:
C               SET NUEFFEA=1.0 IF DON'T WISH TO SCALE COLLISIONALITY
      
      TMP = 1.0
      IF (NPROFUE.EQ.1) TMP = SQRT(2.0)*PI
      IF (NPROFUE.EQ.2) TMP = 4.0*SQRT(PI)/3.0
      IF (NPROFUE.EQ.3) TMP = PI/SQRT(2.0)

      IF (NPROFUE.EQ.0) THEN
         WRITE(*,*) 'UNIFORM ELECTRON COLL., NPROFUE=',NPROFUE
         DO I=1,NRP1
            GNUE(I)  = NUEFFEA
            GNUEM(I) = NUEFFEA
         ENDDO
      ELSEIF (NPROFUE.GE.1.AND.NPROFUE.LE.3) THEN
         WRITE(*,*) 'ANALYTICAL ELECTRON COLL., NPROFUE=',NPROFUE
         ZTEMP = TMP*ZQE**2.5/(4.0*PI*ZQEPS0)**2/SQRT(ZQMI)
         ZTEMP = ZTEMP*ESPECIES_Z(2)**4/SQRT(ESPECIES_M(2))
         ZTEMP = ZTEMP*ZNE0/ZTE0**1.5
         ZTEMP = ZTEMP*ZTAUA0*NUEFFEA
         DO I=1,NRP1
            GNUE(I) = ZTEMP*RHO(I)/TEMPE(I)**1.5*
     &                (17.3-0.5*LOG(ZNE0*RHO(I)/1.0E+20)
     &                     +1.5*LOG(ZTE0*TEMPE(I)/1.0E+3))   
         ENDDO
         DO I=1,NR
            GNUEM(I) = ZTEMP*RHOM(I)/TEMPEM(I)**1.5*
     &                 (17.3-0.5*LOG(ZNE0*RHOM(I)/1.0E+20)
     &                      +1.5*LOG(ZTE0*TEMPEM(I)/1.0E+3))   
         ENDDO
         GNUEM(NRP1) = GNUEM(NR)
      ELSEIF (NPROFUE.EQ.4) THEN
         WRITE(*,*) 'READING ELECTRON COLL. PROFILE, NPROFUE=', NPROFUE
         OPEN(99,FILE='PROFNUE.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP*NUEFFEA
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFUE=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(GNUE,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(GNUEM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFUE'
      ENDIF

C     NEOCLASSICAL FUNCTION TO BE USED IN ROSENBLUTH-PUTVINSKI 
C     RE AVALANCHE MODEL
C     NPROFNF = 0: CONSTANT 1      
C               1: ANALYTICAL FORMULA (LARGE ASPECT RATIO CIRCULAR)
C               4: NUMRICAL PROFILE (CAN BE COMPUTED BY RUNNING MARS-K)
      IF (NPROFNF.EQ.0) THEN
         GNEOFUNC  = 1.0
         GNEOFUNCM = 1.0
      ELSEIF (NPROFNF.EQ.1) THEN
         DO I=1,NRP1
            ZTEMP        = CS(I)/ASPCT
            GNEOFUNC(I)  = 1./(1.+1.46*SQRT(ZTEMP)+1.72*ZTEMP)
            ZTEMP        = CSM(I)/ASPCT
            GNEOFUNCM(I) = 1./(1.+1.46*SQRT(ZTEMP)+1.72*ZTEMP)
         ENDDO
      ELSEIF (NPROFNF.EQ.4) THEN
         OPEN(99,FILE='PROFNEOFUNC.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFNF=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(GNEOFUNC,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(GNEOFUNCM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFNF'
      ENDIF

C     TOROIDAL MOMEMTUM DIFFUSION COEFFICIENT        
C     NPROFVD = 0-3: ANALYTICAL      
C               4: EXPERIMENTAL PROFILE 
      
      IF (NPROFVD.EQ.0) THEN
         DO I=1,NRP1
            TCHIMI(I)  = 1.0
            TCHIMM(I) = 1.0
         ENDDO
      ELSEIF (NPROFVD.EQ.1) THEN
         DO I=2,NR
            TCHIMI(I) = CS(I)**TPOWVD
            TCHIMM(I) = CSM(I)**TPOWVD
         ENDDO
         TCHIMI(1)    = TCHIMM(2)
         TCHIMM(1)    = TCHIMM(2)
         TCHIMI(NRP1) = TCHIMM(NR)
         TCHIMM(NRP1) = TCHIMM(NR)
      ELSEIF (NPROFVD.EQ.2) THEN
         ZTEMP = TEMPE(1)**(-1.5)
         DO I=2,NR
            TCHIMI(I) = TEMPE(I)**(-1.5)/ZTEMP/CS(I)
            TCHIMM(I) = TEMPEM(I)**(-1.5)/ZTEMP/CSM(I)
         ENDDO
         TCHIMI(1)    = TCHIMM(2)
         TCHIMM(1)    = TCHIMM(2)
         TCHIMI(NRP1) = TCHIMM(NR)
         TCHIMM(NRP1) = TCHIMM(NR)
      ELSEIF (NPROFVD.EQ.3) THEN
         ZTEMP = TEMPE(1)**(-1.5)
         DO I=1,NR
            TCHIMI(I) = TEMPE(I)**(-1.5)/ZTEMP
            TCHIMM(I) = TEMPEM(I)**(-1.5)/ZTEMP
         ENDDO
         TCHIMI(NRP1) = TCHIMM(NR)
         TCHIMM(NRP1) = TCHIMM(NR)
      ELSEIF (NPROFVD.EQ.4) THEN
         OPEN(99,FILE='PROFTVD.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFVD=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(TCHIMI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(TCHIMM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFVD'
      ENDIF
      TCHIMI = TCHIM0*TCHIMI
      TCHIMM = TCHIM0*TCHIMM

C     TOROIDAL MOMEMTUM PINCH COEFFICIENT        
C     NPROFVP = 0: ANALYTICAL      
C               4: EXPERIMENTAL PROFILE 
      
      IF (NPROFVP.EQ.0) THEN
         DO I=1,NRP1
            TVPINCHI(I)  = 1.0
            TVPINCHM(I) = 1.0
         ENDDO
      ELSEIF (NPROFVP.EQ.4) THEN
         OPEN(99,FILE='PROFTVP.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFVP=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(TVPINCHI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(TVPINCHM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFVP'
      ENDIF
      TVPINCHI = TVPINCH0*TVPINCHI
      TVPINCHM = TVPINCH0*TVPINCHM

C     RADIAL PARTICLE DIFFUSION COEFFICIENT        
C     NPROFDD = 0-3: ANALYTICAL      
C               4: EXPERIMENTAL PROFILE 
      
      IF (NPROFDD.EQ.0) THEN
         DO I=1,NRP1
            TCHIDI(I) = 1.0
            TCHIDM(I) = 1.0
         ENDDO
      ELSEIF (NPROFDD.EQ.1) THEN
         DO I=2,NR
            TCHIDI(I) = CS(I)**TPOWDD
            TCHIDM(I) = CSM(I)**TPOWDD
         ENDDO
         TCHIDI(1)    = TCHIDM(2)
         TCHIDM(1)    = TCHIDM(2)
         TCHIDI(NRP1) = TCHIDM(NR)
         TCHIDM(NRP1) = TCHIDM(NR)
      ELSEIF (NPROFDD.EQ.2) THEN
         ZTEMP = TEMPE(1)**(-1.5)
         DO I=2,NR
            TCHIDI(I) = TEMPE(I)**(-1.5)/ZTEMP/CS(I)
            TCHIDM(I) = TEMPEM(I)**(-1.5)/ZTEMP/CSM(I)
         ENDDO
         TCHIDI(1)    = TCHIDM(2)
         TCHIDM(1)    = TCHIDM(2)
         TCHIDI(NRP1) = TCHIDM(NR)
         TCHIDM(NRP1) = TCHIDM(NR)
      ELSEIF (NPROFDD.EQ.3) THEN
         ZTEMP = TEMPE(1)**(-1.5)
         DO I=1,NR
            TCHIDI(I) = TEMPE(I)**(-1.5)/ZTEMP
            TCHIDM(I) = TEMPEM(I)**(-1.5)/ZTEMP
         ENDDO
         TCHIDI(NRP1) = TCHIDM(NR)
         TCHIDM(NRP1) = TCHIDM(NR)
      ELSEIF (NPROFDD.EQ.4) THEN
         OPEN(99,FILE='PROFTDD.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFDD=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(TCHIDI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(TCHIDM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFDD'
      ENDIF
      TCHIDI = TCHID0*TCHIDI
      TCHIDM = TCHID0*TCHIDM

C     TOROIDAL EXB ROTATION FREQUENCY                
C     NPROFWE = 0: CONSTANT        
C     NPROFWE = 1: ANALYTICAL: ROT-OMEGA*      
C               4: EXPERIMENTAL PROFILE 
      
      IF (NPROFWE.EQ.0) THEN
         DO I=1,NRP1
            ROTWEI(I) = ROTWE0
            ROTWEM(I) = ROTWE0
         ENDDO
      ELSEIF (NPROFWE.EQ.1) THEN
C        ROTWE = ROT - OMEGAI*
C        TO BE ADDED LATER
      ELSEIF (NPROFWE.EQ.4) THEN
         OPEN(99,FILE='PROFWE.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         IF (NEXPV.EQ.1) ROTWE0 = ZTEMP*ZTAUA0 
         RPRV  = ROTWE0*RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFWE=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(ROTWEI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(ROTWEM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSEIF (NPROFWE.EQ.5) THEN
         OPEN(99,FILE='PROFWE.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = MAXVAL(RPRV)
         RPRV  = ROTWE0*RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFWE=5'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(ROTWEI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(ROTWEM,SRM,NRP1-1,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFWE'
      ENDIF

      IF (.NOT.ALLOCATED(OMEGAE0)) ALLOCATE(OMEGAE0(NRP1,2))
      DO I=1,NRP1
         OMEGAE0(I,1)=ROTWEI(I)
         OMEGAE0(I,2)=ROTWEM(I)
      ENDDO

C     PARALLEL THERMAL TRANSPORT COEFFICIENT
C     NPROFTTCA = 0-3: ANALYTICAL      
C                   4: EXPERIMENTAL PROFILE 
      IF (NPROFTTCA.EQ.0) THEN
         DO I=1,NRP1
            TTCPARAI(I) = 1.0
            TTCPARAM(I) = 1.0
         ENDDO
      ELSEIF (NPROFTTCA.EQ.4) THEN
         OPEN(99,FILE='PROFTTCPARA.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFTTCA=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(TTCPARAI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(TTCPARAM,SRM,NR,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFTTCA'
      ENDIF
      TTCPARAI = TTCPARA0*TTCPARAI
      TTCPARAM = TTCPARA0*TTCPARAM

C     PERPENDICULAR THERMAL TRANSPORT COEFFICIENT
C     NPROFTTCE = 0-3: ANALYTICAL      
C                   4: EXPERIMENTAL PROFILE 
      IF (NPROFTTCE.EQ.0) THEN
         DO I=1,NRP1
            TTCPERPI(I) = 1.0
            TTCPERPM(I) = 1.0
         ENDDO
      ELSEIF (NPROFTTCE.EQ.4) THEN
         OPEN(99,FILE='PROFTTCPERP.IN',STATUS='OLD',FORM='FORMATTED')
         READ(99,*) NRAD,IRAD
         ALLOCATE(RPRS(NRAD),RPRV(NRAD),RPR2(NRAD))
         DO I=1,NRAD
            READ(99,*) RPRS(I),RPRV(I)
         ENDDO
         CLOSE(99)
         ZTEMP = RPRV(1)
         RPRV  = RPRV/ZTEMP
         IF (IRAD.EQ.1) THEN
            SR  = CS(1:NRP1)
            SRM = CSM(1:NRP1)
         ELSEIF (IRAD.EQ.2) THEN
            ZTEMP = RPRS(NRAD)
            RPRS  = RPRS/ZTEMP
            SR  = RAD
            SRM = RADM
         ELSE
            STOP 'NPROFTTCE=4'
         ENDIF
         SR(NRP1) = 1.-1.E-10
         IF (SRM(NRP1-1).GE.1.) SRM(NRP1-1)=1.-1.E-10
         CALL SPLINE1D(TTCPERPI,SR,NRP1,RPRV,RPRS,NRAD,RPR2)
         CALL SPLINE1D(TTCPERPM,SRM,NR,RPRV,RPRS,NRAD,RPR2)
         DEALLOCATE(RPRS,RPRV,RPR2)
      ELSE
         STOP 'NPROFTTCE'
      ENDIF
      TTCPERPI = TTCPERP0*TTCPERPI
      TTCPERPM = TTCPERP0*TTCPERPM

C     REDEFINE TTCPARA TO FACILITATE CODE IMPLEMENTATION
      TTCPARAI = TTCPARAI - TTCPERPI
      TTCPARAM = TTCPARAM - TTCPERPM

C     COMPUTE ADDITIONAL EQUILIBRIUM PROFILES
C     OMEGASI = OMEGASIN + ETASTI*OMEGASIT
C     OMEGASE = OMEGASEN + ETASTE*OMEGASET
C     NOTE THAT OMEGAS*N INVOLVES DRHO/DPSI=(DRHO/DS)/(DPSI/DS)      
      DO J=2,NR
         H1 = (CS(J)-CS(J-1))/2
         H2 = (CS(J+1)-CS(J))/2
         DLNRHO(J)  = ((H1/H2*RHOM(J)-H2/H1*RHOM(J-1))/(H1+H2)-
     &                (H1-H2)*RHO(J)/H1/H2)/RHO(J)
         ZTEMP      = B0K/OMEGACI0/RHO(J)*(ETASTI*PPEQ(J)+(1.-ETASTI)*
     &                                     DLNRHO(J)*PEQ(J)/DPSIDS(J))
         OMEGASI(J) = -ALPHAP*ZTEMP   
         ZTEMP      = B0K/OMEGACI0/RHO(J)*(ETASTE*PPEQ(J)+(1.-ETASTE)*
     &                                     DLNRHO(J)*PEQ(J)/DPSIDS(J))
         OMEGASE(J) = (1.0-ALPHAP)*ZTEMP   
      ENDDO
      DO J=1,NR
         H1 = CS(J+1)-CS(J)
         DLNRHOM(J)  = (RHO(J+1)-RHO(J))/H1/RHOM(J) 
         ZTEMP       = B0K/OMEGACI0/RHOM(J)*(ETASTI*PPEQM(J)+
     &                 (1.-ETASTI)*DLNRHOM(J)*PEQM(J)/DPSIDSM(J))
         OMEGASIM(J) = -ALPHAP*ZTEMP   
         ZTEMP       = B0K/OMEGACI0/RHOM(J)*(ETASTE*PPEQM(J)+
     &                 (1.-ETASTE)*DLNRHOM(J)*PEQM(J)/DPSIDSM(J))
         OMEGASEM(J) = (1.0-ALPHAP)*ZTEMP   
      ENDDO
      J = 1
      OMEGASI(J) = OMEGASI(J+1)
      OMEGASE(J) = OMEGASE(J+1)
      DLNRHO(J)  = DLNRHO(J+1)
      J = NRP1
      OMEGASI(J) = OMEGASI(J-1)
      DLNRHO(J)  = DLNRHO(J-1)
      OMEGASIM(J)= 0.
      OMEGASEM(J)= 0.
      DLNRHOM(J) = 0.

      DO J=2,NR
         H1 = (CS(J)-CS(J-1))/2
         H2 = (CS(J+1)-CS(J))/2
         DOMEGASI(J) = (H1/H2*OMEGASIM(J)-H2/H1*OMEGASIM(J-1))/(H1+H2)-
     &                 (H1-H2)*OMEGASI(J)/H1/H2
         DOMEGASE(J) = (H1/H2*OMEGASEM(J)-H2/H1*OMEGASEM(J-1))/(H1+H2)-
     &                 (H1-H2)*OMEGASE(J)/H1/H2
      ENDDO
      DO J=1,NR
         H1 = CS(J+1)-CS(J)
         DOMEGASIM(J)  = (OMEGASI(J+1)-OMEGASI(J))/H1
         DOMEGASEM(J)  = (OMEGASE(J+1)-OMEGASE(J))/H1
      ENDDO
      J = 1
      DOMEGASI(J) = DOMEGASIM(J)
      DOMEGASE(J) = DOMEGASEM(J)
      J = NRP1
      DOMEGASI(J) = DOMEGASIM(J-1)
      DOMEGASE(J) = DOMEGASEM(J-1)
      DOMEGASIM(J)= 0.
      DOMEGASEM(J)= 0.

      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
C     PRINT PHYSICAL QUANTITIES
      WRITE(*,*) 'PHYSICAL QUANTITIES :'
      WRITE(*,*) 'R0EXP[M]=      ',R0EXP
      WRITE(*,*) 'B0EXP[T]=      ',B0EXP
      WRITE(*,*) 'NE0[1/M^3]=     ',ZNE0
      WRITE(*,*) 'TAUA0[S]=       ',ZTAUA0
      WRITE(*,*) 'TI0[EV]=       ',ZTI0
      WRITE(*,*) 'TE0[EV]=       ',ZTE0
      WRITE(*,*) 'OMEGA[RAD/S]=  ',ROTE/ZTAUA0
      WRITE(*,*) 'OMEGAE[RAD/S]= ',ROTWE0/ZTAUA0
      WRITE(*,*) 'GNUI(0)= ',GNUI(1)
      WRITE(*,*) 'GNUE(0)= ',GNUE(1)

C     PRINT DIMENSIONLESS QUANTITIES
      WRITE(*,*) 'OMEGACI0=      ',OMEGACI0
      WRITE(*,*) 'ALPHAP=         ',ALPHAP
      WRITE(*,*) 'FRACPTH=        ',FRACPTH
      WRITE(*,*) 'OMEGA=          ',ROTE
      WRITE(*,*) 'OMEGAE=         ',ROTWE0
      ENDIF
      
C     OUTPUT EQUILIBRIUM PROFILES
 120  FORMAT(20(E16.9,1X))
 
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      
      OPEN(99,FILE='PROFEQ.OUT')
      DO I=1,NRP1
         WRITE(99,120) CS(I),Q(I),JEQ(I),PEQ(I),RHO(I),
     &                 ROT(I),RESIST(I),GAMARR(I),GMUNU(I),
     &                 TEMPI(I),TEMPE(I),
     &                 DPSIDS(I),T(I),OMEGASI(I),OMEGASE(I),
     &                 CSV(I),RAD(I),GNUI(I),GNUE(I),ZEFFI(I)
      ENDDO
      CLOSE(99)
 130  FORMAT(' CS',15X,'Q',16X,'PEQ',14X,'RHO',14X,'ROT')
      WRITE(*,130) 
      DO I=1,NRP1
         WRITE(*,120) CS(I),Q(I),PEQ(I),RHO(I),ROT(I)
      ENDDO
      
      ENDIF

C     SAVE EQUILIBRIUM ROTATION AND DENSITY PROFILES
      ROTEQ   = ROT
      ROTEQM  = ROTM
      DROTEQ  = DROT
      DROTEQM = DROTM
      RHOEQ   = RHO
      RHOEQM  = RHOM

      DEALLOCATE(RAD,RADM,SR,SRM)

      RETURN
      END

*DECK RADF2RADP
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C------- CONVERT FROM POLOIDAL FLUX TO TOROIDAL FLUX -------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine RADF2RADP(cs,csm,q,rado,radom,npsi)
      implicit none
      integer::npsi,i
      real*8,dimension(npsi)::rado,q,radom,cs,csm
      real*8,dimension(:),allocatable::psi,phi,psim,phim
      real*8::phiend
c
      allocate(psi(npsi),phi(npsi),psim(npsi),phim(npsi))
      do i=1,npsi
         psi(i)=cs(i)**2
      end do
      do i=1,npsi
         psim(i)=csm(i)**2
      end do
      phi(1)=0.
      do i=2,npsi
         phi(i)=phi(i-1)+(q(i)+q(i-1))/2.*(psi(i)-psi(i-1))
      end do
      phim(1)=(phi(2)+phi(1))/2.
      do i=2,npsi
         phim(i)=phim(i-1)+q(i)*(psim(i)-psim(i-1))
      end do
      phiend=phi(npsi)
      phi=phi/phiend
      phim=phim/phiend
      rado=sqrt(phi)
c      phim(npsi)=phim(npsi-1)
      radom=sqrt(phim)
c      write (*,'(2x,"i",8x,"cs",8x,"psi",11x,"q",8x,
c     &     "rado",8x,"phi",7x,"radom",8x,"phim")') 
c      do i=1,npsi
c         write (*,'(i4,1p7e12.4)')i,cs(i),psi(i),q(i),rado(i),phi(i),
c     & radom(i),phim(i)
c      end do
      deallocate(psi,phi,psim,phim)
      return 
      end

*DECK FTRFC
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C------- SUM FOURIER SERIES FOR PERTURBATION ---------------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE FTRFC(PF, PR, NFOU)
C     =============================
      USE DIMENSIM
      USE GLOBALM
C
      INTEGER   NFOU
      COMPLEX*16   PF(MSMAX), PR(NFOU)
      INTEGER   J, MS
      REAL*8      ZA
C
      DO 10 J = 1,NFOU
 10   PR(J) = 0.
C
      DO 20 MS = 1,MSMAX
      DO 20 J = 1,NFOU
      ZA = 2.*PI*DFLOAT(J-1)/DFLOAT(NFOU)
      ZA = ZA * RM(MS,2)
      PR(J) = PR(J) + PF(MS) * CMPLX(COS(ZA), SIN(ZA))
 20   CONTINUE
C
      RETURN
      END
*DECK PLOTLP
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------PLOT  PAMS VECTOR--------- A.B.  02/02/90 --------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
      SUBROUTINE PLOTLP(X,Y)
C     ======================
C
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      INCLUDE 'comioc.inc'
      INTEGER    I,IV,MS
      REAL*8     TMP,TMPR,TMPI
      COMPLEX*16    X(NXCOMP,MSMAX,*),Y(NYCOMP,MSMAX,*)
      DO 10 I = 1,NV
      IV = NR + I
      IF (I.GT.1) Q(IV) = DFLOAT(I)
      QM(IV)= DFLOAT(I)+0.5
 10   CONTINUE
      IF (NTP1.GT.NRP1) Q(NTP1) = NVP1
C
      DO MS=1,MSMAX
         IF (ABS(RM(MS,2)).LT.0.1) THEN
            DO I=1,NRP1
               X(KXB1,MS,I) = X(KXB1,MS,I)*T(I)
            ENDDO
            DO I=NRP1,NTP1
               X(KXB1,MS,I) = X(KXB1,MS,I)*T(NRP1)
            ENDDO
         ENDIF
      ENDDO

      IF (DPLOLP( 1)) CALL PLOTV(X(KXV1,1,1),Q ,'V1U   ',NXCOMP,NRP1)
      IF (DPLOLP( 2)) CALL PLOTV(Y(KYV2,1,1),QM,'V2U   ',NYCOMP,NR  )
      IF (DPLOLP( 3).AND.KYV3.GT.0) 
     &                CALL PLOTV(Y(KYV3,1,1),QM,'V3U   ',NYCOMP,NR  )

      IF (DPLOLP( 4)) CALL PLOTV(X(KXB1,1,1),Q ,'B1U   ',NXCOMP,NTP1)
      IF (DPLOLP( 5)) CALL PLOTV(Y(KYB2,1,1),QM,'B2U   ',NYCOMP,NTOT)
      IF (DPLOLP( 6)) CALL PLOTV(Y(KYB3,1,1),QM,'B3U   ',NYCOMP,NTOT)

      IF (DPLOLP( 7)) CALL PLOTV(Y(KYJ1,1,1),QM,'J1U   ',NYCOMP,NTOT)
      IF (DPLOLP( 8)) CALL PLOTV(X(KXJ2U,1,1),Q,'J2U   ',NXCOMP,NTP1)
      IF (DPLOLP( 9)) CALL PLOTV(X(KXJ3,1,1),Q ,'J3U   ',NXCOMP,NR  )
      IF (DPLOLP(10).AND.KXJ2L.GT.0) 
     &                CALL PLOTV(X(KXJ2L,1,1),Q,'J2L   ',NXCOMP,NTP1)

      IF (DPLOLP(11)) CALL PLOTV(Y(KYPR,1,1),QM,'PRE   ',NYCOMP,NR  )
      IF (DPLOLP(12).AND.KYPE.GT.0) 
     &                CALL PLOTV(Y(KYPE,1,1),QM,'PEE   ',NYCOMP,NR  )
      IF (DPLOLP(13).AND.KYPP.GT.0) 
     &                CALL PLOTV(Y(KYPP,1,1),QM,'PEP   ',NYCOMP,NR  )
      IF (DPLOLP(14).AND.KXPD.GT.0) 
     &                CALL PLOTV(X(KXPD,1,1),Q ,'PDE   ',NXCOMP,NR  )

      IF (DPLOLP(15).AND.KXX1.GT.0) 
     &                CALL PLOTV(X(KXX1,1,1),Q ,'X1U   ',NXCOMP,NRP1)
      IF (DPLOLP(16).AND.KYX2.GT.0) 
     &                CALL PLOTV(Y(KYX2,1,1),QM,'X2U   ',NYCOMP,NR  )
      IF (DPLOLP(17).AND.KYPPERP.GT.0) 
     &                CALL PLOTV(Y(KYPPERP,1,1),QM,'PPERP ',NYCOMP,NR  )
      IF (DPLOLP(18).AND.KYPPARA.GT.0) 
     &                CALL PLOTV(Y(KYPPARA,1,1),QM,'PPARA ',NYCOMP,NR  )
      IF (DPLOLP(19).AND.KYRHOP.GT.0) 
     &                CALL PLOTV(Y(KYRHOP,1,1),QM,'RHOP ',NYCOMP,NR  )
      IF (DPLOLP(20).AND.KYX3.GT.0) 
     &                CALL PLOTV(Y(KYX3,1,1),QM,'X3U  ',NYCOMP,NR  )

      DO MS=1,MSMAX
         IF (ABS(RM(MS,2)).LT.0.1) THEN
            DO I=1,NRP1
               X(KXB1,MS,I) = X(KXB1,MS,I)/T(I)
            ENDDO
            DO I=NRP1,NTP1
               X(KXB1,MS,I) = X(KXB1,MS,I)/T(NRP1)
            ENDDO
         ENDIF
      ENDDO

C-----YQLIU, 2006-02-23
C-----OUTPUT J3U FOR MANUAL MESH-PACKING
      OPEN(CHOUTP,FILE='J3UMAX.OUT',FORM='FORMATTED')
      REWIND(CHOUTP)
      DO I=1,NR
         TMPR = 0
         TMPI = 0
         DO MS=1,MSMAX
            TMP = ABS(DREAL(X(KXJ3,MS,I)))
            IF (TMPR.LT.TMP) TMPR = TMP
            TMP = ABS(DIMAG(X(KXJ3,MS,I)))
            IF (TMPI.LT.TMP) TMPI = TMP
         ENDDO
         WRITE(CHOUTP,20) Q(I),TMPR,TMPI
      ENDDO
      CLOSE(CHOUTP)
20    FORMAT(3(E14.6))

      RETURN
      END
*DECK PLOTV
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-------PLOT ON LINE-PRINTER-------A. BONDESON 18/10/89-----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
      SUBROUTINE PLOTV(C,QA,STRING,NCOMP,NP)
C     ======================================
C
      USE DIMENSIM
      USE GLOBALM
      CHARACTER*6 STRING
      CHARACTER*132 LINE
      INTEGER     J,IPAGE,NCOMP,I,M,NP
      PARAMETER   (IPAGE=79)
      REAL*8        CMAX,CMIN,CAMAX,Z,QA(NP),ZR,ZI
      COMPLEX*16     C(NCOMP,MSMAX,*)
C
      CMAX  = 0.
      CMIN  = 0.
      CAMAX = 0.
C
      DO 10 M=1,MSMAX
      DO 10 I=1,NP
      ZR = DREAL(C(1,M,I))
      ZI = DIMAG(C(1,M,I))
      CAMAX = DMAX1(CAMAX,ZR**2 + ZI**2)
      CMAX=DMAX1(CMAX,ZR)
      CMAX=DMAX1(CMAX,ZI)
      CMIN=DMIN1(CMIN,ZR)
      CMIN=DMIN1(CMIN,ZI)
10    CONTINUE
      CAMAX = SQRT(CAMAX)
      IF (CMAX.EQ.CMIN) CMAX=CMIN+1.
C
      WRITE(*,1030) STRING,CAMAX,CMAX,CMIN
      WRITE(*,1000) STRING
      DO 100 I=1,NP
      DO 20 J=1,IPAGE
20    LINE(J:J)=' '
      Z = -CMIN/(CMAX-CMIN)
      J= INT( Z*(IPAGE-10)+1.5)
      LINE(J:J) = '0'
C
      DO 30 M=1,MSMAX
C     Z = (DREAL((1.,-1.)*C(1,M,I))-CMIN)/(CMAX-CMIN)
      Z = (DREAL(C(1,M,I))-CMIN)/(CMAX-CMIN)
      J= INT( Z*(IPAGE-10)+1.5)
30    LINE(J:J) = CHAR(ICHAR('0')+M)
      WRITE(*,1020) QA(I),LINE(1:IPAGE-9)
 100  CONTINUE
C
C     IF (.TRUE.) RETURN
      WRITE(*,1010) STRING
C
      DO 200 I=1,NP
      DO 120 J=1,IPAGE
120   LINE(J:J)=' '
      Z = -CMIN/(CMAX-CMIN)
      J= INT( Z*(IPAGE-10)+1.5)
      LINE(J:J) = '0'
C
      DO 130 M=1,MSMAX
      Z = (DIMAG(C(1,M,I))-CMIN)/(CMAX-CMIN)
      J= INT( Z*(IPAGE-10)+1.5)
130   LINE(J:J) = CHAR(ICHAR('0')+M)
      WRITE(*,1020) QA(I),LINE(1:IPAGE-9)
 200  CONTINUE
      RETURN
1000  FORMAT(//,A,' - REAL',//)
1010  FORMAT(//,A,' - IMAG',//)
1020  FORMAT(F10.6,2X,A)
1030  FORMAT(//,A,' PEAK AMPLITUDE',E12.4,' MAX =',E12.4,' MIN =',E12.4)
      END
*DECK GCONTR
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C--COMPUTE POWER DISSIPATION DUE TO PARALLEL VISCOSITY--LIU YQ 12.5.2003
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE PPVISC
C     ==========================================================
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      INCLUDE 'comioc.inc'
      INTEGER    I,MSA,NSA
      REAL*8     TMP,EPS
      COMPLEX*16 TMP1,TMP2,CTMP1,CTMP2,CTMP3
C
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE::CPPVISC
      COMPLEX,DIMENSION(:,:),ALLOCATABLE::DPPVISC
C
      ALLOCATE(CPPVISC(NRP1,MSDIM,2),DPPVISC(NRP1,MSDIM))
C
      NSA = 2
      DO I = 1,NRP1
         DO MSA = 1,MSDIM
            DPPVISC(I,MSA) = 0.0
            CPPVISC(I,MSA,1) = 0.0
            CPPVISC(I,MSA,2) = 0.0
         ENDDO
      ENDDO

      DO I = 1,NR
         DO MSA = 1,MSMAX
            TMP = PVISC
            CPPVISC(I,MSA,1) = TMP*0.5
            CPPVISC(I,MSA,2) = TMP*0.5

            TMP = TMP*ABS(RM(MSA,NSA)/QM(I) + RN(NSA))
     &            * SQRT(REAL(PEQM(I))/RHOM(I))
            TMP1 = CONJG(V3U(I,MSA))*V3U(I,MSA)
            TMP2 = -TMP*TMP1
            DPPVISC(I,MSA) = DPPVISC(I,MSA) + TMP2
          ENDDO
       ENDDO


      OPEN(CHOUTP,FILE='PPVISC.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,120) CSM(I),(CPPVISC(I,MSA,1),MSA=1,MSMAX),
     &        (CPPVISC(I,MSA,2),MSA=1,MSMAX),
     &        (REAL(DPPVISC(I,MSA)),MSA=1,MSMAX)
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(200(E13.6,1X))
      DEALLOCATE(CPPVISC,DPPVISC)
      RETURN
      END


      SUBROUTINE ELINTK(RM,ELINT)
C
C     NOTE: ARGUMENT RM IS K**2, NOT K
C
      IMPLICIT LOGICAL (A-Z)
      REAL*8 ELINT, RM
      REAL*8 Z,SUMA,SUMB
      REAL*8 A0,A1,A2,A3,A4,B0,B1,B2,B3,B4
C
       PARAMETER (A0 = 1.38629 436112D0, B0 = .5D0,
     &            A1 = 0.09666 344259D0, B1 = .12498 593597D0,
     &            A2 = 0.03590 092383D0, B2 = .06880 248576D0,
     &            A3 = 0.03742 563713D0, B3 = .03328 355346D0,
     &            A4 = 0.01451 196212D0, B4 = .00441 787012D0)
C
      IF (RM.GE.0. AND .RM.LT.1.) GOTO 10
      WRITE(*,*) ' ARGUMENT FOR K =',RM
      STOP
 10   Z = 1. - RM
      SUMA = A0 + Z * (A1 + Z * (A2 + Z * (A3 + Z * A4)))
      SUMB = B0 + Z * (B1 + Z * (B2 + Z * (B3 + Z * B4)))
C
      ELINT = SUMA - SUMB * LOG(Z)
C
      RETURN
      END

      subroutine lagrange2(num_order,ne,xx,f,xi,lag,lag1,lag2)
      implicit none
      integer num_order,ne,num_on 
      real*8 xx(ne),f(ne),xi,lag,lag1,lag2
 
      integer i,j,m,i2,j2,j3,ic
      real*8 div,l,l1,l2,l22,mul_a1,mul_a2
      real*8,dimension(:),allocatable::a
c
      allocate(a(100))
c      
c      print'(/"lagrange Order= ",i3)',num_order

      if(xx(1) .gt. xx(ne)) then
         if (xi .ge. xx(1)) then
            num_on=1
            goto 222
         else if( xi .le. xx(ne)) then
            num_on=ne-num_order
            goto 222
         else
            do i=1,ne
               if (xi.le.xx(i) .and. xi.ge.xx(i+1)) then
                  ic=i
                  goto 111
               end if
            end do
         end if
      else
         if (xi .le. xx(1)) then
            num_on=1
            goto 222
         else if( xi .ge. xx(ne)) then
            num_on=ne-num_order
            goto 222
         else
            do i=1,ne
               if (xi.ge.xx(i) .and. xi.le.xx(i+1)) then
                  ic=i
                  goto 111
               end if
            end do
         end if
      end if
      
 111  continue
      if (ic .le. ne-num_order)   num_on=ic
      if (ic .ge. ne-num_order+1) num_on=ne-num_order
      
 222  continue
      
      if(1.eq.11) print'("i=",i4," xx=",f9.5," f=",f9.5)',
     $     (i,xx(i),f(i),i=num_on,num_on+num_order)

c      print*,'num_on=',num_on, '    num_order=',num_order
c      stop
      
      lag=0.0 
      lag1=0.0
      lag2=0.d0
      do i=num_on,num_on+num_order
         div=1.d0 
         m=0  
         l=1.d0
         do j=num_on,num_on+num_order
            if (j .ne.  i) then
               m=m+1
               l=l*(xi-xx(j))/(xx(i)-xx(j))
               a(m)=xi-xx(j)
               div=div*(xx(i)-xx(j))
            end if
         end do
         lag=lag+f(i)*l
         
         l1=0.d0
         l22=0.d0
         do i2=1,num_order
               mul_a1=1.d0
               l2=0.d0
!               print'("i2=",i3/)',i2

               do j2=1,num_order
                  if( j2 .ne. i2 ) then
                     mul_a1=mul_a1*a(j2)
                     mul_a2=1.d0
!                     print'("  j2=",i3/)',j2

                     do j3=1,num_order
                        if( j3 .ne. j2 .and. j3 .ne. i2 ) then
!                           print'("    j3=",i3)',j3

                           mul_a2=mul_a2*a(j3)
                        end if
                     end do   ! j3
                     l2=l2+mul_a2
                  end if 
               end do  ! j2
               l1=l1+mul_a1
               l22=l22+l2
         end do !  i2
         lag1=lag1+f(i)*l1/div
         lag2=lag2+f(i)*l22/div
      end do  ! i
c
      deallocate(a)
c
      return
      End
*DECK RWALLG2
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C------- GET GRAD(S) = D_LOG(ZCND)/D_CHI * DG33L ON WALL ---------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
c     INCLUDE 'process.inc'
      SUBROUTINE RWALLG2(JW,K)
C     =================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      REAL*8,DIMENSION(:),ALLOCATABLE::GRADS,GRADS2,GRADS3,GRADS4
      REAL*8,    DIMENSION(:,:),ALLOCATABLE::RKKU
      REAL*8,    DIMENSION(:,:),ALLOCATABLE::GRAD2DS
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::FGRAD2DS
      INTEGER I,J,JW,K
 
      IF (.NOT. ALLOCATED(FGRADS)) THEN
         ALLOCATE(FGRADS(MEDIM))
      END IF
      ALLOCATE(GRADS(NCHI),GRADS2(NCHI),GRADS3(NCHI),GRADS4(NCHI))
C
      ALLOCATE(RKKU    (NVEQ1,  NCHI))
      ALLOCATE(GRAD2DS (   1   ,NCHI))
      ALLOCATE(FGRAD2DS(   1   ,MEDIM))
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'RWALLG2'
C
      DO J         = 1,NVEQ1
      DO I         = 1,NCHI
      RKKU(J,I)    = 0.0
      ENDDO
      ENDDO
C
      DO I         = 1,NCHI
      GRAD2DS(1,I) = 0.0
      ENDDO
C
      NVSTRT    =  IWALLJ
      NVDIM0    =  NVEQ1
      NVDIM1    =  NVSTRT
C
      call FFTDRIVER(RKKU, VG33JC, BCKWD,  NVDIM0, NVDIM1, NVSTRT
     &                    ,MEDIM,  NCHI,   KUOUT,  IERSUB, IERPLC
     &                                                    ,IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VG33JC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      1,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      DO J      = 1,NCHI
      GRADS2(J) = RKKU(NVSTRT,J)
      ENDDO
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RKKU,   VG33JC,    NVDIM0,  NVDIM1,  NVSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'GRADS2')
      ENDIF
C
      call FFTDRIVER(RKKU, VJAC,   BCKWD,  NVDIM0, NVDIM1, NVSTRT
     &                    ,MEDIM,  NCHI,   KUOUT,  IERSUB, IERPLC
     &                                                    ,IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VJAC  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      2,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      DO J      = 1,NCHI
      GRADS3(J) = RKKU(NVSTRT,J)
      ENDDO
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RKKU,   VJAC,      NVDIM0,  NVDIM1,  NVSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'GRADS3')
      ENDIF
C
      call FFTDRIVER(RKKU, VG33L,  BCKWD,  NVDIM0, NVDIM1, NVSTRT
     &                    ,MEDIM,  NCHI,   KUOUT,  IERSUB, IERPLC
     &                                                    ,IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: VG33L  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      3,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      DO J      = 1,NCHI
      GRADS4(J) = RKKU(NVSTRT,J)
      ENDDO
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(RKKU,   VG33L,     NVDIM0,  NVDIM1,  NVSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'GRADS4')
      ENDIF
 
C     GET GRAD(S) AS A FUNCTION OF ANGLE ON THE WALL
      IF (K.EQ.1) THEN
         DO J=1,NCHI
            GRADS(J) = ZCND(JW,J)**2*GRADS3(J)
         ENDDO
      ELSEIF (K.EQ.2) THEN
         DO J=1,NCHI
            GRADS(J) = ZCND(JW,J)*GRADS4(J)
         ENDDO
      ELSEIF (K.EQ.3) THEN
         DO J=1,NCHI
            GRADS(J) = ZCNDC(JW,J)*GRADS4(J) + ZCND(JW,J)*GRADS2(J)*.5
         ENDDO
      ELSEIF (K.EQ.4) THEN
         DO J=1,NCHI
            GRADS(J) = GRADS3(J)
         ENDDO
      ENDIF
C
C
      NDSTRT    =  1
      NDM0      =  1
      NDM1      =  NDSTRT
      DO J              = 1,NCHI
      GRAD2DS(NDSTRT,J) = GRADS(J)
      ENDDO
C
      call FFTDRIVER(GRAD2DS,FGRAD2DS, FORWD,  NDM0,   NDM1
     &                      ,NDSTRT,   MEDIM,  NCHI
     &                      ,KUOUT,    IERSUB, IERPLC, IERR)

      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error: FGRADS  in ',IERSUB
         call ABORTRUN
     &         (SUBNAM,      4,   MESSAGE
     &         ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      DO J        = 1,MEDIM
      FGRADS(J)   = FGRAD2DS(NDSTRT,J)
      ENDDO
C
      IF(KUFFTP .GT. 0) THEN
        call FFTOUTPT(GRAD2DS,FGRAD2DS,  NDM0,    NDM1,    NDSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'FGRAD2DS')
      ENDIF
C
C    DEALLOCATE SPACE
C
      DEALLOCATE(GRADS,GRADS2,GRADS3,GRADS4)
      DEALLOCATE(RKKU,GRAD2DS,FGRAD2DS)
      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C---------OUTPUT RMZM_F-------------Y.Q.LIU 20/03/2012-----------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PRINTRMZM
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      INCLUDE 'comioc.inc'
      INTEGER I,J,NO4
      REAL*8  ZTEMP(4)
      
C     OUTPUT RZ-COORDINATES IN FOURIER SPACE
         NO4=CHOUTP
         OPEN(UNIT=NO4,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='RMZM_F.OUT')
         WRITE(NO4,1019) MOUTPUT,NRP1,NVEQ1-1,R0EXP
         DO I=1,NRP1
            IF (I.EQ.1) THEN
               WRITE(NO4,1021) CS(I),NRATSURF,QPLS(I),B0EXP
            ELSEIF (I.GT.1.AND.I.LE.NRATSURF+1) THEN
               WRITE(NO4,1021) CS(I),IRATSURF(I-1),QPLS(I),B0EXP
            ELSE
               WRITE(NO4,1020) CS(I),CS(I),QPLS(I),B0EXP
            ENDIF
         ENDDO
         DO I=2,NVEQ1
            WRITE(NO4,1020) VCS(I),VCS(I),VCS(I),VCS(I)
         ENDDO
         DO J=1,MOUTPUT
            DO I=1,NRP1
               ZTEMP(1)=REAL(RPF(I,J))
               ZTEMP(2)=IMAG(RPF(I,J))
               ZTEMP(3)=REAL(ZPF(I,J))
               ZTEMP(4)=IMAG(ZPF(I,J))
               WRITE(NO4,1020) ZTEMP(1),ZTEMP(2),ZTEMP(3),ZTEMP(4)
            ENDDO
            DO I=2,NVEQ1
               ZTEMP(1)=REAL(RVF(I,J))
               ZTEMP(2)=IMAG(RVF(I,J))
               ZTEMP(3)=REAL(ZVF(I,J))
               ZTEMP(4)=IMAG(ZVF(I,J))
               WRITE(NO4,1020) ZTEMP(1),ZTEMP(2),ZTEMP(3),ZTEMP(4)
            ENDDO
         ENDDO
 1019    FORMAT(I5,1X,I5,1X,I5,1X,E15.8)
 1020    FORMAT(E15.8,1X,E15.8,1X,E15.8,1X,E15.8)
 1021    FORMAT(E15.8,1X,I5,1X,E15.8,1X,E15.8)
         CLOSE(UNIT=NO4)

         RETURN
         END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C POSTPROCESSING OUTPUT DATA (PERTURBED FIELD)         LIU YQ 03.01.2013
C ON A RECTANGULAR UNIFORM (R,Z)-GRID                                  $
C MAPPING FROM FLUX SURFACE BASED (R,Z)-GRID TO RECTANGULAR (R,Z)-GRID $
C VIA A TEMPLATE UNIT SQUARE                                           $   
C NEW NAMELIST VARIABLES: ORMIN,ORMAX,OZMIN,OZMAX,NORR,NOZZ            $
C NOTE THAT ORMIN,ORMAX,OZMIN,OZMAX ARE DEFINED IN [M]                 $
C                                                                      $
C POSTPROCESSING OUTPUT DATA (VECTOR POTENTIAL)        LIU YQ 03.02.2014
C USING SPECIAL GUAGE CONDITION OF A_PHI=0 FOR N.NE.0 PERTURBATIONS    $
C ==> A_R=-(I*R/N)*B_Z, A_Z=(I*R/N)*B_R
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE OUTPUT_RECTRZ(ISW)
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      USE RCOMDM  
      USE GVACUUMM
      IMPLICIT NONE
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INTEGER    I,J,MS,KCHECK,IC,JC,KCLOC,KCGLB,INCIRC,ISW,NVMIN
      REAL*8     ODRDS,ODRDC,ODZDS,ODZDC,OJAC,ORADIUS,HCHI,HORR,HOZZ,
     &           AA1,AA2,AA3,AA4,BB1,BB2,BB3,BB4,
     &           OPP,OQQ,OSS,ODD,OX1,OX2,OY1,OY2,OXX,OYY,
     &           WAR,WAZ,WAP,R00,Z00
C
      COMPLEX*16 CTMP1
      REAL*8,DIMENSION(:),ALLOCATABLE::OCHI
      REAL*8,DIMENSION(:,:),ALLOCATABLE::ORM,OZM,ORC,OZC,OCS,OCC
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::OB1S,OB1,OB2,OB3,
     &                                       OBR,OBZ,OBP,
     &                                       OBCR,OBCZ,OBCP,
     &                                       OACR,OACZ,OACP
C
      INTEGER    NDELP,NDELV,KSW
      REAL*8     CDELB
      COMPLEX*16 DELB,BTEMP(6)
      INTEGER    K,L,ICC(10),JCC(10)
      CHARACTER(LEN=64) FILENAME

      KCHECK = 0
      KSW = 1
      IF (NCASE.NE.6) KSW = -2
      IF ((NCASE.EQ.6.OR.NCASE.EQ.10).AND.
     &    OSWEEP.GE.1.AND.OSWEEP.LE.NSWEEP) KSW = MOD(ISW,OSWEEP)
      IF ((NCASE.EQ.6.OR.NCASE.EQ.10).AND.ISW.EQ.NSWEEP) KSW = -1

C     DEFINE (R00,Z00) AS AVERAGE OF (R(1,:),Z(1,:))
      R00 = 0.0
      Z00 = 0.0
      DO J=1,NCHI
         R00 = R00 + REQ(1,J)
         Z00 = Z00 + ZEQ(1,J)
      ENDDO
      R00 = R00/NCHI
      Z00 = Z00/NCHI

      ALLOCATE( OCHI(NCHI+1) )
      DO J=1,NCHI+1
         OCHI(J) = 2.*PI*(J-1)/NCHI
      ENDDO

C     STEP 1: GENERATE FLUX SURFACE BASED (R,Z)-GRID IN WHOLE DOMAIN   
C             USING HALF-INTEGER RADIAL POINTS
      ALLOCATE( ORM(NTOT,NCHI+1), OZM(NTOT,NCHI+1) )

      DO I=1,NR
         DO J=1,NCHI
            ORM(I,J) = REQM(I,J)
            OZM(I,J) = ZEQM(I,J)
         ENDDO
         J = NCHI+1
         ORM(I,J) = REQM(I,1)
         OZM(I,J) = ZEQM(I,1)
      ENDDO
      
      DO I=1,NV
         DO J=1,NCHI
            ORM(NR+I,J) = VRRM(I,J)
            OZM(NR+I,J) = VRZM(I,J)
         ENDDO
         J = NCHI+1
         ORM(NR+I,J) = VRRM(I,1)
         OZM(NR+I,J) = VRZM(I,1)
      ENDDO

C     SMOOTH MAGNETIC FIELD AT PLASMA-VACUUM BOUNDARY FOR VACUUM RUNS
C     THIS PROCEDURE DOES NOT REALLY HELP, BECAUSE THE METRICS TENSOR
C     ALSO HAS JUMPS ACROSS THE INTERFACE
      IF (INCFEED.EQ.4.AND.1.EQ.0) THEN
      NDELP = 1
      DO I=NR,1,-1
         IF (CSM(I).GE.0.99) NDELP = NR-I+1
      ENDDO
      NDELV = IFEED-1
      DO I=1,NWALL
         IF (NDELV.GT.IWALL(I)) NDELV = IWALL(I)
      ENDDO
      CDELB  = 0.5
      WRITE(*,*) 'SMOOTH VACUUM FIELD',NDELP,NDELV,CDELB
      DO MS=1,MSMAX
         DELB = B2U(NRP1,MS)-B2U(NR,MS)
         DO I=NR-NDELP,NR
            B2U(I,MS) = B2U(I,MS) +
     &      (CSM(I)-CSM(NR-NDELP))/(1.-CSM(NR-NDELP))*CDELB*DELB
         ENDDO
         DO I=NR+1,NR+NDELV
            B2U(I,MS) = B2U(I,MS) -
     &      (CSM(NR+NDELV)-CSM(I))/(CSM(NR+NDELV)-1.)*(1.-CDELB)*DELB
         ENDDO

         DELB = B3U(NRP1,MS)-B3U(NR,MS)
         DO I=NR-NDELP,NR
            B3U(I,MS) = B3U(I,MS) +
     &      (CSM(I)-CSM(NR-NDELP))/(1.-CSM(NR-NDELP))*CDELB*DELB
         ENDDO
         DO I=NR+1,NR+NDELV
            B3U(I,MS) = B3U(I,MS) -
     &      (CSM(NR+NDELV)-CSM(I))/(CSM(NR+NDELV)-1.)*(1.-CDELB)*DELB
         ENDDO
      ENDDO
      ENDIF
C
C     RECOMPUTE B3U FROM DIVB=0
C     NO NEED IF IPDIVB=1 
C     IF (IPDIVB.EQ.0.AND.ABS(RNTOR).GT.0.) THEN
      IF (1.EQ.1.AND.ABS(RNTOR).GT.0.) THEN
         DO MS=1,MSMAX
         DO I=1,NTOT
            B3U(I,MS) =-((B1U(I+1,MS)-B1U(I,MS))/(CS(I+1)-CS(I)) +
     &                   CI*RM(MS,2)*B2U(I,MS))/(CI*RNTOR)
         ENDDO
         ENDDO

         DO MS=1,MSMAX
            IF (ABS(RM(MS,2)).GT.0) THEN
               B3U(1,MS) = B3U(2,MS)
               B2U(1,MS) =-((B1U(2,MS)-B1U(1,MS))/(CS(2)-CS(1)) +
     &                   CI*RNTOR*B3U(1,MS))/(CI*RM(MS,2))
            ELSE
               B2U(2,MS) = B2U(3,MS)
               B2U(1,MS) = B2U(3,MS)
            ENDIF
         ENDDO
      ENDIF
C
      
      WRITE(*,*) 'OUTPUT_RECTRZ: AFTER STEP 1'

C     STEP 2: COMPUTE (B1,B2,B3) IN (ORM,OZM)-GRID
      ALLOCATE( OB1(NTOT,NCHI+1), OB2(NTOT,NCHI+1), OB3(NTOT,NCHI+1) )
      ALLOCATE( OB1S(NTOT,NCHI+1) ) 

      OB1 = (0.,0.)
      OB2 = (0.,0.)
      OB3 = (0.,0.)
      OB1S= (0.,0.)

      HCHI = 2.*PI/NCHI
      DO J=1,NCHI+1
         DO MS=1,MSMAX
            CTMP1 = EXP(CI*RM(MS,2)*(J-1)*HCHI)
            DO I=1,NTOT
               OB1(I,J) = OB1(I,J) + (B1U(I,MS)+B1U(I+1,MS))*.5*CTMP1
               OB2(I,J) = OB2(I,J) + B2U(I,MS)*CTMP1
               OB3(I,J) = OB3(I,J) + B3U(I,MS)*CTMP1
               OB1S(I,J)= OB1S(I,J) + B1U(I,MS)*CTMP1
            ENDDO
         ENDDO
      ENDDO

C     RECOMPUTE OB3 FROM DIVB=0
      IF (1.EQ.0.AND.ABS(RNTOR).GT.0.) THEN
         DO I=1,NTOT-1
         DO J=2,NCHI
            OB3(I,J) = -((OB1S(I+1,J)-OB1S(I,J))/(CS(I+1)-CS(I)) +
     &                   (OB2(I,J+1)-OB2(I,J-1))/HCHI/2.)/(CI*RNTOR)
         ENDDO
         OB3(I,1) = -((OB1S(I+1,1)-OB1S(I,1))/(CS(I+1)-CS(I)) +
     &                (OB2(I,2)-OB2(I,NCHI))/HCHI/2.)/(CI*RNTOR)
         OB3(I,NCHI+1) = OB3(I,1)
         ENDDO
      ENDIF

C     FIRST FIND NVMIN, SUCH THAT ALL VRRM>0
      AA1   = B0EXP*R0EXP*R0EXP
C     AA2   = 1.0/(R0EXP*R0EXP)
      AA2   = 1.0
      NVMIN = NV
      DO I=NV,1,-1
      DO J=1,NCHI
         IF (VRRM(I,J).LE.0.0) NVMIN = I
      ENDDO         
      ENDDO         
      NVMIN = NVMIN-1

C     SAVE MARS-F RAW DATA FOR LOCUST
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='BPLASMA_MARSF.OUT')
      ELSE
         WRITE(FILENAME,"(A14,I0.4,A4)") "BPLASMA_MARSF",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      WRITE(CHOUTP,118) INT(RNTOR),MOUTPUT,M1,M2,NR,NVMIN,IFEED,
     &                  ABS(FEEDI(1))*R0EXP*B0EXP/4.0E-7/PI,FINIC
      DO I=1,NR+1
         WRITE(CHOUTP,119) CS(I),CSM(I),Q(I)
      ENDDO
      DO I=NR+2,NR+NVMIN
         WRITE(CHOUTP,119) CS(I),CSM(I),0.0
      ENDDO

      DO J=1,MOUTPUT
      DO I=1,NR
         WRITE(CHOUTP,119) REAL(RPF(I,J))*R0EXP,IMAG(RPF(I,J))*R0EXP,
     &                     REAL(ZPF(I,J))*R0EXP,IMAG(ZPF(I,J))*R0EXP,  
     &                     REAL(RPFM(I,J))*R0EXP,IMAG(RPFM(I,J))*R0EXP,
     &                     REAL(ZPFM(I,J))*R0EXP,IMAG(ZPFM(I,J))*R0EXP
      ENDDO
      DO I=1,NVMIN
         WRITE(CHOUTP,119) REAL(RVF(I,J))*R0EXP,IMAG(RVF(I,J))*R0EXP,
     &                     REAL(ZVF(I,J))*R0EXP,IMAG(ZVF(I,J))*R0EXP,  
     &                     REAL(RVFM(I,J))*R0EXP,IMAG(RVFM(I,J))*R0EXP,
     &                     REAL(ZVFM(I,J))*R0EXP,IMAG(ZVFM(I,J))*R0EXP
      ENDDO
      ENDDO

      DO J=1,MSMAX
      DO I=1,NR+NVMIN
         WRITE(CHOUTP,119) REAL(B1U(I,J))*AA1,IMAG(B1U(I,J))*AA1,
     &                     REAL(B2U(I,J))*AA1,IMAG(B2U(I,J))*AA1,
     &                     REAL(B3U(I,J))*AA1,IMAG(B3U(I,J))*AA1
      ENDDO
      ENDDO

      DO J=1,MSMAX
      DO I=1,NR+1
         WRITE(CHOUTP,119) REAL(X1U(I,J))*AA2,IMAG(X1U(I,J))*AA2
      ENDDO
      ENDDO
      CLOSE(CHOUTP)

 118  FORMAT(7(I4,1X),2(E17.10,1X))
 119  FORMAT(8(E17.10,1X))

      ENDIF

      WRITE(*,*) 'OUTPUT_RECTRZ: AFTER STEP 2'

C     STEP 3: COMPUTE (BR,BZ,BPHI) IN (ORM,OZM)-GRID
C             AND (AR,AZ) IN THE SAME GRID
      ALLOCATE( OBR(NTOT,NCHI+1), OBZ(NTOT,NCHI+1), OBP(NTOT,NCHI+1) )

      DO J=1,NCHI+1
         DO I=1,NTOT
            IF (I.LE.NR.AND.J.LE.NCHI) OJAC = RJAM(I,J)     
            IF (I.LE.NR.AND.J.GT.NCHI) OJAC = RJAM(I,1)     
            IF (I.GT.NR.AND.J.LE.NCHI) OJAC = VRJAM(I-NR,J)     
            IF (I.GT.NR.AND.J.GT.NCHI) OJAC = VRJAM(I-NR,1)     

            IF (I.LE.NR.AND.J.LE.NCHI) ODRDS = 
     &                                 (REQ(I+1,J)-REQ(I,J))/CSH(I)
            IF (I.LE.NR.AND.J.GT.NCHI) ODRDS = 
     &                                 (REQ(I+1,1)-REQ(I,1))/CSH(I)
            IF (I.GT.NR.AND.J.LE.NCHI) ODRDS = 
     &                           (VRR(I+1-NR,J)-VRR(I-NR,J))/VCSH(I-NR)
            IF (I.GT.NR.AND.J.GT.NCHI) ODRDS = 
     &                           (VRR(I+1-NR,1)-VRR(I-NR,1))/VCSH(I-NR)

            IF (I.LE.NR.AND.J.LE.NCHI) ODZDS = 
     &                                 (ZEQ(I+1,J)-ZEQ(I,J))/CSH(I)
            IF (I.LE.NR.AND.J.GT.NCHI) ODZDS = 
     &                                 (ZEQ(I+1,1)-ZEQ(I,1))/CSH(I)
            IF (I.GT.NR.AND.J.LE.NCHI) ODZDS = 
     &                           (VRZ(I+1-NR,J)-VRZ(I-NR,J))/VCSH(I-NR)
            IF (I.GT.NR.AND.J.GT.NCHI) ODZDS = 
     &                           (VRZ(I+1-NR,1)-VRZ(I-NR,1))/VCSH(I-NR)

            IF (2.LE.J.AND.J.LE.NCHI)  ODRDC = 
     &                                 (ORM(I,J+1)-ORM(I,J-1))/HCHI/2.
            IF (J.EQ.1.OR.J.EQ.NCHI+1) ODRDC = 
     &                                 (ORM(I,2)-ORM(I,NCHI))/HCHI/2.

            IF (2.LE.J.AND.J.LE.NCHI)  ODZDC = 
     &                                 (OZM(I,J+1)-OZM(I,J-1))/HCHI/2.
            IF (J.EQ.1.OR.J.EQ.NCHI+1) ODZDC = 
     &                                 (OZM(I,2)-OZM(I,NCHI))/HCHI/2.
C
C           RECOMPUTE OJAC USING DEFINITION
C           TESTED ON AN ITER EQUILIBRIUM ==> DOES NOT AFFECT RESULTS
C           OJAC = (ODRDS*ODZDC-ODRDC*ODZDS)*ORM(I,J)

            OBR(I,J) = (OB1(I,J)*ODRDS+OB2(I,J)*ODRDC)/OJAC
            OBZ(I,J) = (OB1(I,J)*ODZDS+OB2(I,J)*ODZDC)/OJAC
            OBP(I,J) = OB3(I,J)*ORM(I,J)/OJAC
         ENDDO
      ENDDO

      DEALLOCATE( OB1,OB2,OB3, OB1S )

C     SMOOTH B-FIELD AT PLASMA-VACUUM INTERFACE
      IF (1.EQ.1.AND.NV.GT.3) THEN
         DO J=1,NCHI+1

         DO K=1,3
            BTEMP(1:6) = OBR(NR-2:NR+3,J)
            DO I=2,5
               OBR(NR+I-3,J) = (BTEMP(I-1)+BTEMP(I)+BTEMP(I+1))/3.0
            ENDDO
         ENDDO 
 
         DO K=1,3
            BTEMP(1:6) = OBZ(NR-2:NR+3,J)
            DO I=2,5
               OBZ(NR+I-3,J) = (BTEMP(I-1)+BTEMP(I)+BTEMP(I+1))/3.0
            ENDDO
         ENDDO 
 
         DO K=1,3
            BTEMP(1:6) = OBP(NR-2:NR+3,J)
            DO I=2,5
               OBP(NR+I-3,J) = (BTEMP(I-1)+BTEMP(I)+BTEMP(I+1))/3.0
            ENDDO
         ENDDO 

         ENDDO 
      ENDIF
C
C     SAVE B-FIELD AT (S,CHI)-MESH
C     THE DATA CAN BE USED TO DO MAPPING OUTSIDE THE CODE (E.G. USING 
C     GRIDDATA IN MATLAB)
C     THE ITER TEST RESULST SHOW SAME RESULTS AS THAT FROM STEP 5
      IF (1.EQ.0) THEN
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='BPLASMA_PSICHI.OUT')
      ELSE
         WRITE(FILENAME,"(A14,I0.4,A4)") "BPLASMA_PSICHI",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      DO I=1,NTOT
      DO J=1,NCHI+1
C        WRITE(CHOUTP,120) ORM(I,J)*R0EXP,OZM(I,J)*R0EXP,
         WRITE(CHOUTP,120) CSM(I),HCHI*(J-1),
     &                     REAL(OBR(I,J))*B0EXP,IMAG(OBR(I,J))*B0EXP,
     &                     REAL(OBZ(I,J))*B0EXP,IMAG(OBZ(I,J))*B0EXP,
     &                     REAL(OBP(I,J))*B0EXP,IMAG(OBP(I,J))*B0EXP
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF
      ENDIF
C

      WRITE(*,*) 'OUTPUT_RECTRZ: AFTER STEP 3'

      IF (NORR.GT.1.AND.NOZZ.GT.1) THEN

C     STEP 4: CREATING RECTANGULAR (R,Z)-GRID
      ALLOCATE( ORC(NORR,NOZZ), OZC(NORR,NOZZ) )

      HORR = (ORMAX-ORMIN)/DFLOAT(NORR-1)/R0EXP
      HOZZ = (OZMAX-OZMIN)/DFLOAT(NOZZ-1)/R0EXP
      DO J=1,NOZZ
         DO I=1,NORR
            ORC(I,J) = ORMIN/R0EXP + HORR*(I-1)
            OZC(I,J) = OZMIN/R0EXP + HOZZ*(J-1)
         ENDDO
      ENDDO     

C     STEP 5: FOR EACH POINT OF (ORC,OZC)-GRID, 
C             FIND CONTAINING QUADRALATERAL IN (ORM,OZM)-MESH
C             AND INTERPOLATE FIELD
C             USE THE SAME PROCEDURE TO DO INVERSE MAPPING 
C             FROM RECTANGULAR (R,Z)-MESH TO (S,CHI)-MESH
C             I.E. FROM (ORC,OZC) TO (OCS,OCC)
      ALLOCATE( OBCR(NORR,NOZZ), OBCZ(NORR,NOZZ), OBCP(NORR,NOZZ) )
      ALLOCATE( OCS(NORR,NOZZ), OCC(NORR,NOZZ) )

      OBCR = (0.,0.)
      OBCZ = (0.,0.)
      OBCP = (0.,0.)
      OCS  = 0.
      OCC  = 0.

C     COMPUTE AVERAGE RADIUS OF THE INNEREST CIRCLE OF THE (ORM,OZM)-MESH
      INCIRC  = 1
      ORADIUS = 0.
      DO J=1,NCHI
         ORADIUS = ORADIUS + SQRT((ORM(INCIRC,J)-R00)**2+
     &                            (OZM(INCIRC,J)-Z00)**2)
      ENDDO
      ORADIUS = ORADIUS/NCHI

C     LOOPING THROUGH ALL (ORC,OZC)-POINTS
      DO IC=1,NORR
      DO JC=1,NOZZ

C     SPECIAL TREATMENT IF (ORC,OZC)-POINT BELONGS TO INNEREST CIRCLE
      IF (SQRT((ORC(IC,JC)-R00)**2+(OZC(IC,JC)-Z00)**2).LE.
     &   ORADIUS*1.2) THEN
         DO J=1,NCHI
            OBCR(IC,JC) = OBCR(IC,JC) + OBR(INCIRC,J)
            OBCZ(IC,JC) = OBCZ(IC,JC) + OBZ(INCIRC,J)
            OBCP(IC,JC) = OBCP(IC,JC) + OBP(INCIRC,J)
         ENDDO
         OBCR(IC,JC) = OBCR(IC,JC)/NCHI
         OBCZ(IC,JC) = OBCZ(IC,JC)/NCHI
         OBCP(IC,JC) = OBCP(IC,JC)/NCHI
         OCS(IC,JC)  = 0.
         OCC(IC,JC)  = DATAN2(OZC(IC,JC)-Z00,ORC(IC,JC)-R00)
         IF (OCC(IC,JC).LT.0.) OCC(IC,JC) = OCC(IC,JC) + 2.*PI
      ELSE
         KCGLB = 0
         DO I=INCIRC,NTOT-1
         DO J=1,NCHI
            AA1 = ORM(I+1,J)-ORM(I,J) 
            BB1 = OZM(I+1,J)-OZM(I,J) 
            AA2 = ORM(I,J+1)-ORM(I,J)
            BB2 = OZM(I,J+1)-OZM(I,J)
            AA3 = ORM(I,J)+ORM(I+1,J+1)-ORM(I+1,J)-ORM(I,J+1)
            BB3 = OZM(I,J)+OZM(I+1,J+1)-OZM(I+1,J)-OZM(I,J+1)
            AA4 = ORM(I,J)
            BB4 = OZM(I,J)

            OPP = AA3*BB1 - AA1*BB3
            OQQ = AA2*BB1 - AA1*BB2 + (ORC(IC,JC)-AA4)*BB3 - 
     &                                (OZC(IC,JC)-BB4)*AA3
            OSS = (ORC(IC,JC)-AA4)*BB2 - (OZC(IC,JC)-BB4)*AA2
            
            IF (OPP.EQ.0.) THEN
               IF (OQQ.EQ.0.) STOP 'OUTPUT_RECTRZ: 1'
               OX1 = -OSS/OQQ
               OX2 = OX1
            ELSE
               ODD = OQQ*OQQ - 4.*OPP*OSS
               IF (ODD.GT.0.) THEN
                  OX1 = (-OQQ+SQRT(ODD))/OPP/2.
                  OX2 = (-OQQ-SQRT(ODD))/OPP/2.
               ELSE
                  OX1 = -1.
                  OX2 = -1.
               ENDIF
            ENDIF

            OY1 = ((AA2+AA3*OX1)*(ORC(IC,JC)-AA4-AA1*OX1) + 
     &             (BB2+BB3*OX1)*(OZC(IC,JC)-BB4-BB1*OX1))/
     &            ((AA2+AA3*OX1)**2 + (BB2+BB3*OX1)**2)
            OY2 = ((AA2+AA3*OX2)*(ORC(IC,JC)-AA4-AA1*OX2) + 
     &             (BB2+BB3*OX2)*(OZC(IC,JC)-BB4-BB1*OX2))/
     &            ((AA2+AA3*OX2)**2 + (BB2+BB3*OX2)**2)
            
            KCLOC = 0
            IF (0.LE.OX1.AND.OX1.LE.1..AND.0.LE.OY1.AND.OY1.LE.1.) THEN
               OXX   = OX1
               OYY   = OY1
               KCLOC = 1
            ELSEIF (0.LE.OX2.AND.OX2.LE.1..AND.0.LE.OY2.AND.OY2.LE.1.) 
     &             THEN
               OXX    = OX2
               OYY    = OY2
               KCLOC  = 1
            ENDIF

            IF (KCLOC.EQ.1) THEN
               KCGLB = 1
               OBCR(IC,JC) = OBR(I,J)*(1.-OXX)*(1.-OYY) + 
     &                       OBR(I+1,J)*OXX*(1.-OYY) + 
     &                       OBR(I+1,J+1)*OXX*OYY + 
     &                       OBR(I,J+1)*(1.-OXX)*OYY 
               OBCZ(IC,JC) = OBZ(I,J)*(1.-OXX)*(1.-OYY) + 
     &                       OBZ(I+1,J)*OXX*(1.-OYY) + 
     &                       OBZ(I+1,J+1)*OXX*OYY + 
     &                       OBZ(I,J+1)*(1.-OXX)*OYY 
               OBCP(IC,JC) = OBP(I,J)*(1.-OXX)*(1.-OYY) + 
     &                       OBP(I+1,J)*OXX*(1.-OYY) + 
     &                       OBP(I+1,J+1)*OXX*OYY + 
     &                       OBP(I,J+1)*(1.-OXX)*OYY 
               OCS(IC,JC)  = CSM(I)*(1.-OXX) + 
     &                       CSM(I+1)*OXX
               OCC(IC,JC)  = OCHI(J)*(1.-OYY) + 
     &                       OCHI(J+1)*OYY
            ENDIF
         ENDDO
         ENDDO

         IF (KCGLB.EQ.0) THEN
            WRITE(*,*) '(ORC,OZC)=',ORC(IC,JC),OZC(IC,JC)
            STOP 'OUTPUT_RECTRZ: 3'
         ENDIF

      ENDIF

      ENDDO
      ENDDO

C     SMOOTH VACUUM FIELD ALONG PLASMA-VACUUM INTERFACE
      IF (INCFEED.EQ.4.AND.1.EQ.0) THEN
      CDELB = SQRT(HORR**2+HOZZ**2)
      DO J=1,2*NCHI,5
         L = J
         IF (L.GT.NCHI) L = L-NCHI
         I = 0
         DO IC=1,NORR
         DO JC=1,NOZZ
            IF (SQRT((ORC(IC,JC)-REQ(NRP1,L))**2+
     &               (OZC(IC,JC)-ZEQ(NRP1,L))**2).LE.CDELB) THEN  
            I = I+1
            IF (I.GT.10) STOP 'RECTRZ: TOO GLOBAL SMOOTHING'
            ICC(I) = IC
            JCC(I) = JC
            ENDIF
         ENDDO 
         ENDDO 
         IF (I.EQ.0) WRITE(*,*) 'RECTRZ WARNING: J=',J
         IF (I.GT.0) THEN
            DELB = 0.
            DO K=1,I
               DELB = DELB + OBCR(ICC(K),JCC(K))
            ENDDO
            DELB = DELB/DFLOAT(I)
            DO K=1,I
               OBCR(ICC(K),JCC(K)) = DELB
            ENDDO

            DELB = 0.
            DO K=1,I
               DELB = DELB + OBCZ(ICC(K),JCC(K))
            ENDDO
            DELB = DELB/DFLOAT(I)
            DO K=1,I
               OBCZ(ICC(K),JCC(K)) = DELB
            ENDDO

            DELB = 0.
            DO K=1,I
               DELB = DELB + OBCP(ICC(K),JCC(K))
            ENDDO
            DELB = DELB/DFLOAT(I)
            DO K=1,I
               OBCP(ICC(K),JCC(K)) = DELB
            ENDDO
         ENDIF
      ENDDO   
      ENDIF

C     STEP 6: SAVE B-FIELD ON RECTANGULAR (ORC,OZC)-MESH INTO A FILE
C     CONVERT INTO SI-UNITS BEFORE SAVING
      ORC  = ORC*R0EXP
      OZC  = OZC*R0EXP
      OBCR = OBCR*B0EXP
      OBCZ = OBCZ*B0EXP
      OBCP = OBCP*B0EXP

C     CONVERSION RULE FROM (R,Z,PHI) TO (R,PHI,Z):
C     RNTOR -> -RNTOR; B_PHI -> -B_PHI; A_PHI -> -A_PHI
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='BPLASMA_RECTRZ.OUT')
      ELSE
         WRITE(FILENAME,"(A14,I0.4,A4)") "BPLASMA_RECTRZ",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      DO I=1,NORR
      DO J=1,NOZZ
         WRITE(CHOUTP,120) ORC(I,J),OZC(I,J),
     &                     REAL(OBCR(I,J)),IMAG(OBCR(I,J)),
     &                     REAL(OBCZ(I,J)),IMAG(OBCZ(I,J)),
     &                     REAL(-OBCP(I,J)),IMAG(-OBCP(I,J))
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF

      IF (1.EQ.1) THEN
C     SAVE B-FIELD DATA IN KSTAR FLT FORMAT
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='BPLASMA_KSTAR.OUT')
      ELSE
         WRITE(FILENAME,"(A14,I0.4,A4)") "BPLASMA_KSTAR",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      DO I=1,NORR
      DO J=1,NOZZ
         WRITE(CHOUTP,133) (I-1)*NOZZ+J,ORC(I,J),OZC(I,J),
     &                     REAL(OBCR(I,J)),-IMAG(OBCR(I,J)),
     &                     REAL(OBCZ(I,J)),-IMAG(OBCZ(I,J)),
     &                     REAL(-OBCP(I,J)),-IMAG(-OBCP(I,J))
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF
      ENDIF

 133  FORMAT(I7,1X,8(E15.8,1X))
 120  FORMAT(8(E15.8,1X))
 122  FORMAT(2(I4,1X),2(E15.8,1X))

C     SAVE INVERSE MAPPING DATA
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='SCHIMESH_RECTRZ.OUT')
      ELSE
         WRITE(FILENAME,"(A15,I0.4,A4)") "SCHIMESH_RECTRZ",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      WRITE(CHOUTP,122) NORR,NOZZ,R00*R0EXP,Z00*R0EXP
      DO I=1,NORR
      DO J=1,NOZZ
         WRITE(CHOUTP,120) ORC(I,J),OZC(I,J),OCS(I,J),OCC(I,J)
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF

C     STEP 7: COMPUTE A-FIELD FROM B-FIELD ON RECTANGULAR (ORC,OZC)-MESH 
C             USE COMBINED GAUGE CONDITIONS AR=0, AZ=0, AP=0, WITH
C             WEIGHTING FACTORS WAR, WAZ, WAP, RESPECTIVELY
C             SAVE A-FIELD INTO AN OUTPUT FILE  
C     TEST SHOWS THAT THE PURE AR=0 GAUGE COMBINED WITH SMOOTHING OF
C     AP ALONG Z-AXIS GIVES THE BEST ACCURACY
      WAR = 0.0
      WAZ = 0.0
      WAP = 1.0

      IF (ABS(RNTOR).LT.0.1) THEN
         WAR = WAR + WAP
         WAP = 0.
      ENDIF

      IF (ABS(WAR+WAZ+WAP-1.).LT.1.E-10) THEN
         ALLOCATE( OACR(NORR,NOZZ), OACZ(NORR,NOZZ), OACP(NORR,NOZZ) )
         OACR = (0.,0.)
         OACZ = (0.,0.)
         OACP = (0.,0.)
 
         IF (WAR.GT.0.) THEN
            HCHI = WAR*HORR*R0EXP*0.5
            DO I=2,NORR
               OACZ(I,:) = OACZ(I-1,:) + (OBCP(I-1,:)+OBCP(I,:))*HCHI   
               OACP(I,:) = OACP(I-1,:) - (ORC(I-1,:)*OBCZ(I-1,:)+
     &                                    ORC(I,:)*OBCZ(I,:))*HCHI   
            ENDDO
C           OACP(2:NORR,:) = OACP(2:NORR,:)/ORC(2:NORR,:)
            OACP = OACP/ORC

            DO I=1,20
               OACP(:,2:NOZZ-1) = (OACP(:,1:NOZZ-2)+OACP(:,3:NOZZ))/6.
     &                          + OACP(:,2:NOZZ-1)*2./3.
            ENDDO
         ENDIF

         IF (WAZ.GT.0.) THEN
            HCHI = WAZ*HOZZ*R0EXP*0.5
            DO J=2,NOZZ
               OACR(:,J) = OACR(:,J-1) - (OBCP(:,J-1)+OBCP(:,J))*HCHI   
               OACP(:,J) = OACP(:,J-1) + (OBCR(:,J-1)+OBCR(:,J))*HCHI
            ENDDO
         ENDIF

         IF (WAP.GT.0.) THEN
            OACR =OACR + WAP*ORC*OBCZ/RNTOR/CI
            OACZ =OACZ - WAP*ORC*OBCR/RNTOR/CI
         ENDIF

C        CONVERSION RULE FROM (R,Z,PHI) TO (R,PHI,Z):
C        RNTOR -> -RNTOR; B_PHI -> -B_PHI; A_PHI -> -A_PHI
         IF (KSW.LE.0) THEN
         IF (KSW.LT.0) THEN
            OPEN(CHOUTP,FILE='APLASMA_RECTRZ.OUT')
         ELSE
            WRITE(FILENAME,"(A14,I0.4,A4)") "APLASMA_RECTRZ",ISW,".OUT"
            OPEN(CHOUTP,FILE=TRIM(FILENAME))
         ENDIF
         REWIND(CHOUTP)
         DO I=1,NORR
         DO J=1,NOZZ
            WRITE(CHOUTP,120) ORC(I,J),OZC(I,J),
     &                        REAL(OACR(I,J)),IMAG(OACR(I,J)),
     &                        REAL(OACZ(I,J)),IMAG(OACZ(I,J)),
     &                        REAL(-OACP(I,J)),IMAG(-OACP(I,J))
        ENDDO
        ENDDO
        CLOSE(CHOUTP)
        ENDIF
        DEALLOCATE( OACR,OACZ,OACP )
      ENDIF
C
      
      DEALLOCATE( ORC,OZC,OBCR,OBCZ,OBCP,OCS,OCC )

      ENDIF

      DEALLOCATE( ORM,OZM,OBR,OBZ,OBP,OCHI )

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C POSTPROCESSING OUTPUT DATA (PLASMA DISPLACEMENT)     LIU YQ 06.17.2019
C XI_R and XI_Z ON (S,CHI)-GRID
C MAINLY FOR N=0 PLASMA DISPLACEMENT, ALTHOUGH PROCEDURE IS ALSO VALID 
C FOR N.NE.0
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE OUTPUT_XIRZ(ISW)
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      IMPLICIT NONE
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INTEGER    I,J,MS,KCHECK,ISW,KSW
      REAL*8     B_2M,RW1,RW2,RW3
C
      COMPLEX*16 CTMP1
      REAL*8,DIMENSION(:),ALLOCATABLE::OCHI
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::OXI1,OXI2,OXI3,
     &                                       OXIR,OXIZ
C
      CHARACTER(LEN=1024) FILENAME

      KCHECK = 0

      KSW = 1
      IF (NCASE.NE.6) KSW = -2
      IF ((NCASE.EQ.6.OR.NCASE.EQ.10).AND.
     &    OSWEEP.GE.1.AND.OSWEEP.LE.NSWEEP) KSW = MOD(ISW,OSWEEP)
      IF ((NCASE.EQ.6.OR.NCASE.EQ.10).AND.ISW.EQ.NSWEEP) KSW = -1

C     STEP 1: COMPUTE (OXI1,OXI2,OXI3) IN (S,CHI)-MESH
      ALLOCATE( OCHI(NCHI) )
      ALLOCATE( OXI1(NR,NCHI), OXI2(NR,NCHI), OXI3(NR,NCHI),
     &          OXIR(NR,NCHI), OXIZ(NR,NCHI) )

      OCHI = 0.
      OXI1 = (0.,0.)
      OXI2 = (0.,0.)
      OXI3 = (0.,0.)
      OXIR = (0.,0.)
      OXIZ = (0.,0.)

      DO J=1,NCHI
         OCHI(J) = 2.*PI*(J-1)/NCHI
      ENDDO

      DO J=1,NCHI
         DO MS=1,MSMAX
            CTMP1 = EXP(CI*RM(MS,2)*OCHI(J))
            DO I=1,NR
               OXI1(I,J) = OXI1(I,J) + (X1U(I,MS)+X1U(I+1,MS))*.5*CTMP1
               OXI2(I,J) = OXI2(I,J) + X2U(I,MS)*CTMP1
               OXI3(I,J) = OXI3(I,J) + X3U(I,MS)*CTMP1
            ENDDO
         ENDDO
      ENDDO

C     STEP 2: COMPUTE (OXIR,OXIZ)
      DO J=1,NCHI
      DO I=1,NR
         B_2M = G22LM(I,J)*DPSIDSM(I)**2/RJAM(I,J)**2 +
     &          TM(I)**2/REQM(I,J)**2

         RW1 = (RDCDZM(I,J)+(DPSIDSM(I)/RJAM(I,J))**2*G12LM(I,J)/
     &         B_2M*RDSDZM(I,J))*RJAM(I,J)/REQM(I,J)
         RW2 = -RJAM(I,J)*TM(I)/REQM(I,J)/B_2M*RDSDZM(I,J)
         RW3 = -DPSIDSM(I)/REQM(I,J)*RDSDZM(I,J)
         OXIR(I,J) = RW1*OXI1(I,J) + RW2*OXI2(I,J) + RW3*OXI3(I,J)

         RW1 = ((DPSIDSM(I)**2+TM(I)**2*G11LM(I,J))*RDSDZM(I,J)+
     &          TM(I)**2*G12LM(I,J)*RDCDZM(I,J))/REQM(I,J)**2/B_2M
         RW2 = RJAM(I,J)*TM(I)/DPSIDSM(I)/B_2M*RBZM(I,J)
         RW3 = RBZM(I,J)
         OXIZ(I,J) = RW1*OXI1(I,J) + RW2*OXI2(I,J) + RW3*OXI3(I,J)
      ENDDO
      ENDDO

      OXIR = OXIR*R0EXP
      OXIZ = OXIZ*R0EXP

      VDEXIR = (0.,0.)
      VDEXIZ = (0.,0.)
      DO I=1,NR
      DO J=1,NCHI
         VDEXIR = VDEXIR + OXIR(I,J)
         VDEXIZ = VDEXIZ + OXIZ(I,J)
      ENDDO
      ENDDO
      VDEXIR = VDEXIR/NR/NCHI
      VDEXIZ = VDEXIZ/NR/NCHI

C     SAVE (OXIR,OXIZ) ON (S,CHI)-MESH
      IF (KSW.LE.0) THEN
      IF (KSW.LT.0) THEN
         OPEN(CHOUTP,FILE='XPLASMA_XIRZ.OUT')
      ELSE
         WRITE(FILENAME,"(A12,I0.4,A4)") "XPLASMA_XIRZ",ISW,".OUT"
         OPEN(CHOUTP,FILE=TRIM(FILENAME))
      ENDIF
      REWIND(CHOUTP)
      DO I=1,NR
      DO J=1,NCHI
         WRITE(CHOUTP,120) CSM(I),OCHI(J),
     &                     REAL(OXIR(I,J)),IMAG(OXIR(I,J)),
     &                     REAL(OXIZ(I,J)),IMAG(OXIZ(I,J))
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
      ENDIF
 120  FORMAT(6(E15.8,1X))

      DEALLOCATE(OCHI,OXI1,OXI2,OXI3,OXIR,OXIZ)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C FOR A GIVEN NORMAL FIELD AT A SURFACE OUTSIDE OR INSIDE              $
C THE PLASMA, FIND AN EQUIVALENT SURFACE CURRENT, LOCATED AT ANOTHER   $
C SURFACE, THAT PRODUCES THE SAME VACUUM FIELD.                        $
C CALCULATIONS ARE BASED ON BIOT-SAVART LAW, WHICH REQUIRES THAT:      $
C   R1>0, R2>0 FOR ALL POLOIDAL ANGLES                                 $
C LIU YQ 24.01.2013                                                    $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BS_ESC
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE RCOMDM  
      USE FEEDBACKM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'

      INTEGER    MAX_EC, IBNM_EC, MS, J, KS, NP_EC, J1,J2, KCHECK, 
     &           NT_EC,M1_EC,M2_EC
      PARAMETER  (NP_EC=128, NT_EC=12)
      REAL*8     TEMP1,TEMP2,HCHI,HPHI
      COMPLEX*16 CTMP1,CTMP2
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::BNM_EC,J2U_EC,SCR,CTMPARR
      REAL*8,DIMENSION(:),ALLOCATABLE::ECR1,ECZ1,ECR2,ECZ2,
     &                                 ECDRDC,ECDZDC,ECDSDR,ECDSDZ,
     &                                 ECPC,ECPS,ECPCN,ECPSN
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::ECA

      KCHECK = 3

C     STEP 1.1:
C     READ IN NORMAL FIELD IN TERMS OF POLOIDAL FOURIER HARMONICS,
C     AS WELL AS RADIAL POSITION OF THE SURFACE IN TERMS OF GRID POINT
      OPEN(CHOUTP,FILE='INPUT_BNM.IN')
      REWIND(CHOUTP)
      READ(CHOUTP,*) MAX_EC,IBNM_EC
      IF (MAX_EC.NE.MSMAX) STOP 'BS_ESC: 1'
      IF (IBNM_EC.LT.2.OR.IBNM_EC.GE.NTOT) STOP 'BS_ESC: 2'
      ALLOCATE(BNM_EC(MAX_EC))
      DO MS=1,MAX_EC
         READ(CHOUTP,*) TEMP1,TEMP2
         BNM_EC(MS) = CMPLX(TEMP1,TEMP2)
      ENDDO
      CLOSE(CHOUTP)

      ALLOCATE(ECA(MAX_EC,MAX_EC))

C     STEP 1.2: TRUNCATE NUMBER OF POLOIDAL HARMONICS TO [M1_EC,M2_EC]
C     SUCH THAT AT LEAST NT_EC NUMBER OF POLOIDAL POINTS AVAILABLE 
C     FOR EACH HARMONIC
      MS    = INT(NCHI/NT_EC)
      M1_EC = MIN(ABS(M1),MS)
      M2_EC = MIN(ABS(M2),MS)
      IF (M1.LT.0) M1_EC = -M1_EC
      IF (M2.LT.0) M2_EC = -M2_EC
      IF (KCHECK.GE.1) WRITE(*,*) 'M1_EC,M2_EC=',M1_EC,M2_EC
 
      IF (KREADECA.EQ.0) THEN
C     STEP 2:
C     COMPUTE MATRIX ELEMENTS RELATING EQUIVALENT SURFACE CURRENT TO THE
C     GIVEN RADIAL FIELD BNM_EC, IN POLOIDAL FOURIER SPACE,
C     BASED ON BIOT-SAVART LAW IN MARS COORDINATES.
C     IFEED = THE ASSUMED RADIAL LOCATION OF THE EQUIVALENT SURFACE CURRENT
      
C     STEP 2.1: PREPARE NECESSARY METRICS COEFFICIENTS 
      ALLOCATE( ECR1(NCHI),ECZ1(NCHI),ECR2(NCHI),ECZ2(NCHI),
     &          ECDRDC(NCHI),ECDZDC(NCHI),ECDSDR(NCHI),ECDSDZ(NCHI) )

      DO J=1,NCHI
         ECR1(J) = VRR(IFEED,J)         
         ECZ1(J) = VRZ(IFEED,J)
      ENDDO
      IF (IBNM_EC.LE.NRP1) THEN
      DO J=1,NCHI
         ECR2(J) = REQ(IBNM_EC,J)         
         ECZ2(J) = ZEQ(IBNM_EC,J)
      ENDDO
      ELSE
      DO J=1,NCHI
         ECR2(J) = VRR(IBNM_EC-NR,J)         
         ECZ2(J) = VRZ(IBNM_EC-NR,J)
      ENDDO
      ENDIF

C     TESTING...
C     HCHI = 2.*PI/NCHI
C     DO J=1,NCHI
C        TEMP1 = COS((J-1)*HCHI)
C        TEMP2 = SIN((J-1)*HCHI)
C        ECR1(J) = (10.+1.250000000*TEMP1)/10.
C        ECZ1(J) = ( 0.+1.250000000*TEMP2)/10.
C        ECR2(J) = (10.+1.166666667*TEMP1)/10.
C        ECZ2(J) = ( 0.+1.166666667*TEMP2)/10.
C     ENDDO

      OPEN(CHOUTP,FILE='VRZ.OUT')
      REWIND(CHOUTP)
      DO J=1,NCHI
         WRITE(CHOUTP,111) ECR1(J),ECZ1(J),ECR2(J),ECZ2(J)  
      ENDDO
      CLOSE(CHOUTP)
 111  FORMAT(4(1X,E15.8))
      
      HCHI = 2.*PI/NCHI
      DO J=2,NCHI-1
         ECDRDC(J) = (ECR1(J+1)-ECR1(J-1))/HCHI/2.
         ECDZDC(J) = (ECZ1(J+1)-ECZ1(J-1))/HCHI/2.
         ECDSDR(J) = (ECZ2(J+1)-ECZ2(J-1))/HCHI/2.
         ECDSDZ(J) = (ECR2(J+1)-ECR2(J-1))/HCHI/2.
      ENDDO
      J=1
      ECDRDC(J) = (ECR1(2)-ECR1(NCHI))/HCHI/2.
      ECDZDC(J) = (ECZ1(2)-ECZ1(NCHI))/HCHI/2.
      ECDSDR(J) = (ECZ2(2)-ECZ2(NCHI))/HCHI/2.
      ECDSDZ(J) = (ECR2(2)-ECR2(NCHI))/HCHI/2.
      J=NCHI
      ECDRDC(J) = (ECR1(1)-ECR1(NCHI-1))/HCHI/2.
      ECDZDC(J) = (ECZ1(1)-ECZ1(NCHI-1))/HCHI/2.
      ECDSDR(J) = (ECZ2(1)-ECZ2(NCHI-1))/HCHI/2.
      ECDSDZ(J) = (ECR2(1)-ECR2(NCHI-1))/HCHI/2.

      DO J=1,NCHI
         ECDSDR(J) = ECDSDR(J)*ECR2(J)
         ECDSDZ(J) =-ECDSDZ(J)*ECR2(J)
      ENDDO

      HPHI = 2.*PI/NP_EC
      ALLOCATE(ECPC(NP_EC),ECPS(NP_EC),ECPCN(NP_EC),ECPSN(NP_EC))
      DO J=1,NP_EC
         TEMP1    = HPHI*(J-1)-PI
         ECPC(J)  = COS(TEMP1)
         ECPS(J)  = SIN(TEMP1)
         ECPCN(J) = COS(RNTOR*TEMP1)
         ECPSN(J) = SIN(RNTOR*TEMP1)
      ENDDO

C     STEP 2.2: COMPUTE MATRIX ELEMENTS
      IF (KCHECK.GE.1) WRITE(*,*) 'MAX_EC,IBNM_EC=',MAX_EC,IBNM_EC

      ECA = (0.,0.)
        
      DO MS=M1_EC-M1+1,M2_EC-M1+1
      IF (KCHECK.GE.2) WRITE(*,*) 'MS=',MS
      DO KS=M1_EC-M1+1,M2_EC-M1+1
         DO J1=1,NCHI
         DO J2=1,NCHI
            CTMP1 = EXP((RM(KS,2)*(J1-1)-RM(MS,2)*(J2-1))*HCHI*CI)
            TEMP1 = (ECR2(J2)-ECR1(J1))**2 + (ECZ2(J2)-ECZ1(J1))**2
            DO J=1,NP_EC
               TEMP2 = (TEMP1 + 2.*ECR1(J1)*ECR2(J2)*(1.-ECPC(J)))**1.5
               CTMP2 = (ECDSDR(J2)*(ECR1(J1)*(ECZ1(J1)-ECZ2(J2))*
     &                              ECPC(J)*ECPCN(J)*RM(KS,2)/RNTOR + 
     &                              (ECR1(J1)*ECDZDC(J1)+
     &                               (ECZ2(J2)-ECZ1(J1))*ECDRDC(J1))*
     &                              ECPS(J)*ECPSN(J)*CI) + 
     &                  ECDSDZ(J2)*(ECR1(J1)*ECPCN(J)*(-ECR1(J1)+
     &                              ECR2(J2)*ECPC(J))*RM(KS,2)/RNTOR -
     &                              ECR2(J2)*ECPS(J)*ECPSN(J)*
     &                              ECDRDC(J1)*CI))*CTMP1/TEMP2
               ECA(MS,KS) = ECA(MS,KS) + CTMP2
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      TEMP1 = HCHI**2*HPHI/8.0/PI**2
      ECA   =-ECA*TEMP1*(VCS(IFEED+1)-VCS(IFEED-1))*0.5

      DEALLOCATE(ECR1,ECZ1,ECR2,ECZ2,
     &           ECDRDC,ECDZDC,ECDSDR,ECDSDZ,
     &           ECPC,ECPS,ECPCN,ECPSN)

      OPEN(CHOUTP,FILE='ECA.OUT')
      REWIND(CHOUTP)
      DO MS=1,MAX_EC
      DO KS=1,MAX_EC
         WRITE(CHOUTP,*) REAL(ECA(MS,KS)),IMAG(ECA(MS,KS))
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
      
      ELSE

      OPEN(CHOUTP,FILE='ECA.IN')
      REWIND(CHOUTP)
      DO MS=1,MAX_EC
      DO KS=1,MAX_EC
         READ(CHOUTP,*) TEMP1,TEMP2
         ECA(MS,KS) = CMPLX(TEMP1,TEMP2)
      ENDDO
      ENDDO
      CLOSE(CHOUTP)

      ENDIF

C     STEP 3:
C     INVERT MATRIX ECA AND COMPUTE THE EQUILIVALENT SURFACE CURRENT 
C     AS A VECTOR OF POLOIDAL FOURIER HARMONICS
      DO MS=1,M1_EC-M1
         ECA(MS,:)  = (0.,0.)
         ECA(:,MS)  = (0.,0.)
         ECA(MS,MS) = (1.,0.)
      ENDDO
      DO MS=M2_EC-M1+2,MAX_EC
         ECA(MS,:)  = (0.,0.)
         ECA(:,MS)  = (0.,0.)
         ECA(MS,MS) = (1.,0.)
      ENDDO

      ALLOCATE(SCR(MAX_EC*(2*MAX_EC+1)),CTMPARR(MAX_EC))
      CALL FMIND(ECA,CTMPARR,MAX_EC,MAX_EC,SCR,CTMP1,EPSDET,0,1,IFEED)

      DEALLOCATE(SCR,CTMPARR)

      DO MS=1,M1_EC-M1
         ECA(MS,:)  = (0.,0.)
         ECA(:,MS)  = (0.,0.)
      ENDDO
      DO MS=M2_EC-M1+2,MAX_EC
         ECA(MS,:)  = (0.,0.)
         ECA(:,MS)  = (0.,0.)
      ENDDO

      ALLOCATE(J2U_EC(MAX_EC))
      J2U_EC = (0.,0.)
      DO MS=1,MAX_EC
      DO KS=1,MAX_EC
         J2U_EC(MS) = J2U_EC(MS) + ECA(MS,KS)*BNM_EC(KS)
      ENDDO
      ENDDO

      DEALLOCATE(ECA)

C     STEP 4:
C     SAVE EQUIVALENT SURFACE CURRENT INTO A FILE
      OPEN(CHOUTP,FILE='CURHARMO.OUT')
      REWIND(CHOUTP)
      WRITE(CHOUTP,*) MAX_EC,MAX_EC,MAX_EC
      DO MS=1,MAX_EC
         WRITE(CHOUTP,110) IM(MS,2),REAL(J2U_EC(MS)),IMAG(J2U_EC(MS))
      ENDDO
 110  FORMAT(I4,2(1X,E15.8))
      CLOSE(CHOUTP)
      
      DEALLOCATE(BNM_EC,J2U_EC)

      RETURN
      END
      
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE NORMAL FIELD AT A SURFACE JUST OUTSIDE       LIU YQ 27.02.2013
C THE PLASMA, USING TRUE GEOMETRY OF COILS (I.E. DELTA-FUNCTION) FOR   $
C TOROIDAL COMPONENT OF COIL CURRENT.                                  $
C CALCULATIONS ARE BASED ON BIOT-SAVART LAW, WHICH REQUIRES THAT:      $
C   R1>0, R2>0 FOR ALL POLOIDAL ANGLES                                 $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BS_B1
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE FEEDBACKM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'compam.inc'

      INTEGER    IBNM_EC, J, K, NP_EC, J1,J2, KCHECK
      PARAMETER  (NP_EC=128)
      REAL*8     TEMP1,TEMP2,TEMPA,TEMPB,TEMP3,TEMP4,TEMP5,TEMP6,
     &           HCHI,HPHI
      INTEGER,DIMENSION(:),ALLOCATABLE::JK1,JK2
      REAL*8,DIMENSION(:),ALLOCATABLE::ECC,ECR1,ECZ1,ECR2,ECZ2,
     &                                 ECDRDC,ECDZDC,ECDSDR,ECDSDZ,
     &                                 ECPC,ECPS,ECPCN,ECPSN
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::ECB1,ECBR,ECBZ

      KCHECK = 1

C     STEP 1: 
C     RADIAL POSITION WHERE TO COMPUTE B1
      IBNM_EC = ISENS(1)
 
C     STEP 2.1: PREPARE NECESSARY METRICS COEFFICIENTS 
      ALLOCATE( ECC(NCHI),ECR1(NCHI),ECZ1(NCHI),ECR2(NCHI),ECZ2(NCHI),
     &          ECDRDC(NCHI),ECDZDC(NCHI),ECDSDR(NCHI),ECDSDZ(NCHI) )

      HCHI = 2.*PI/NCHI
      DO J=1,NCHI
         ECC(J)  = HCHI*(J-1)
         ECR1(J) = VRR(IFEED,J)         
         ECZ1(J) = VRZ(IFEED,J)
         ECR2(J) = VRR(IBNM_EC,J)         
         ECZ2(J) = VRZ(IBNM_EC,J)
      ENDDO

      OPEN(CHOUTP,FILE='VRZ.OUT')
      REWIND(CHOUTP)
      DO J=1,NCHI
         WRITE(CHOUTP,111) ECC(J),ECR1(J),ECZ1(J),ECR2(J),ECZ2(J)  
      ENDDO
      CLOSE(CHOUTP)
 111  FORMAT(5(1X,E15.8))
      
C     COMPUTE JK1, JK2
C     NOTE THAT CHI IS DEFINED IN [0,2*PI]
      ALLOCATE( JK1(NCOIL), JK2(NCOIL) )
      JK1 = 1
      JK2 = 1
      DO K=1,NCOIL
         TEMP1 = (FCCHI(K) - 0.5*FWCHI(K))*PI     
         TEMP2 = (FCCHI(K) + 0.5*FWCHI(K))*PI     
         IF (TEMP1.LT.0.) TEMP1 = TEMP1 + 2.*PI
         IF (TEMP2.LT.0.) TEMP2 = TEMP2 + 2.*PI
         DO J=2,NCHI
            IF (ABS(ECC(J)-TEMP1).LT.ABS(ECC(JK1(K))-TEMP1)) JK1(K)=J
            IF (ABS(ECC(J)-TEMP2).LT.ABS(ECC(JK2(K))-TEMP2)) JK2(K)=J
         ENDDO
      ENDDO
      
      IF (KCHECK.EQ.1) THEN
         DO K=1,NCOIL
            WRITE(*,*) 'K,JK1,JK2=',K,JK1(K),JK2(K)
         ENDDO
      ENDIF
      
      DO J=2,NCHI-1
         ECDRDC(J) = (ECR1(J+1)-ECR1(J-1))/HCHI/2.
         ECDZDC(J) = (ECZ1(J+1)-ECZ1(J-1))/HCHI/2.
         ECDSDR(J) = (ECZ2(J+1)-ECZ2(J-1))/HCHI/2.
         ECDSDZ(J) = (ECR2(J+1)-ECR2(J-1))/HCHI/2.
      ENDDO
      J=1
      ECDRDC(J) = (ECR1(2)-ECR1(NCHI))/HCHI/2.
      ECDZDC(J) = (ECZ1(2)-ECZ1(NCHI))/HCHI/2.
      ECDSDR(J) = (ECZ2(2)-ECZ2(NCHI))/HCHI/2.
      ECDSDZ(J) = (ECR2(2)-ECR2(NCHI))/HCHI/2.
      J=NCHI
      ECDRDC(J) = (ECR1(1)-ECR1(NCHI-1))/HCHI/2.
      ECDZDC(J) = (ECZ1(1)-ECZ1(NCHI-1))/HCHI/2.
      ECDSDR(J) = (ECZ2(1)-ECZ2(NCHI-1))/HCHI/2.
      ECDSDZ(J) = (ECR2(1)-ECR2(NCHI-1))/HCHI/2.

      DO J=1,NCHI
         ECDSDR(J) = ECDSDR(J)*ECR2(J)
         ECDSDZ(J) =-ECDSDZ(J)*ECR2(J)
      ENDDO

      HPHI = 2.*PI/NP_EC
      ALLOCATE(ECPC(NP_EC),ECPS(NP_EC),ECPCN(NP_EC),ECPSN(NP_EC))
      DO J=1,NP_EC
         TEMP1    = HPHI*(J-1)-PI
         ECPC(J)  = COS(TEMP1)
         ECPS(J)  = SIN(TEMP1)
         ECPCN(J) = COS(RNTOR*TEMP1)
         ECPSN(J) = SIN(RNTOR*TEMP1)
      ENDDO

C     STEP 2.2: COMPUTE B1 VS. POLOIDAL ANGLE 
      ALLOCATE( ECB1(NCHI),ECBR(NCHI),ECBZ(NCHI) )
      ECB1  = (0.,0.)
      ECBR  = (0.,0.)
      ECBZ  = (0.,0.)
  
      DO J2=1,NCHI
      DO K=1,NCOIL
         TEMP5 = 0.
         TEMP6 = 0.
         TEMP1 = (ECR2(J2)-ECR1(JK1(K)))**2 + (ECZ2(J2)-ECZ1(JK1(K)))**2
         TEMP2 = (ECR2(J2)-ECR1(JK2(K)))**2 + (ECZ2(J2)-ECZ1(JK2(K)))**2
         DO J=1,NP_EC
            TEMPA = (TEMP1 + 2.*ECR1(JK1(K))*ECR2(J2)*(1.-ECPC(J)))**1.5
            TEMPB = (TEMP2 + 2.*ECR1(JK2(K))*ECR2(J2)*(1.-ECPC(J)))**1.5
            TEMP3 =-ECR1(JK1(K))*(ECZ1(JK1(K))-ECZ2(J2))/TEMPA
            TEMP4 =-ECR1(JK2(K))*(ECZ1(JK2(K))-ECZ2(J2))/TEMPB
            TEMP5 = TEMP5 + (TEMP4-TEMP3)*ECPC(J)*ECPCN(J)
            TEMP3 = ECR1(JK1(K))*(ECR1(JK1(K))-ECR2(J2)*ECPC(J))/TEMPA
            TEMP4 = ECR1(JK2(K))*(ECR1(JK2(K))-ECR2(J2)*ECPC(J))/TEMPB
            TEMP6 = TEMP6 + (TEMP4-TEMP3)*ECPCN(J)
         ENDDO

         TEMP3 = 0.
         TEMP4 = 0.
         IF (JK2(K).GT.JK1(K)) THEN
         DO J1=JK1(K),JK2(K)-1
         TEMP1 = (ECR2(J2)-ECR1(J1))**2 + (ECZ2(J2)-ECZ1(J1))**2
         DO J=1,NP_EC
            TEMPA = (TEMP1 + 2.*ECR1(J1)*ECR2(J2)*(1.-ECPC(J)))**1.5
            TEMPB = ECPS(J)*ECPSN(J)/TEMPA
            TEMP3 = TEMP3 + TEMPB*(ECR1(J1)*ECDZDC(J1)+
     &                             (ECZ2(J2)-ECZ1(J1))*ECDRDC(J1))
            TEMP4 = TEMP4 + TEMPB*ECR2(J2)*ECDRDC(J1)
         ENDDO
         ENDDO
         ELSE
         DO J1=1,JK2(K)-1
         TEMP1 = (ECR2(J2)-ECR1(J1))**2 + (ECZ2(J2)-ECZ1(J1))**2
         DO J=1,NP_EC
            TEMPA = (TEMP1 + 2.*ECR1(J1)*ECR2(J2)*(1.-ECPC(J)))**1.5
            TEMPB = ECPS(J)*ECPSN(J)/TEMPA
            TEMP3 = TEMP3 + TEMPB*(ECR1(J1)*ECDZDC(J1)+
     &                             (ECZ2(J2)-ECZ1(J1))*ECDRDC(J1))
            TEMP4 = TEMP4 + TEMPB*ECR2(J2)*ECDRDC(J1)
         ENDDO
         ENDDO

         DO J1=JK1(K),NCHI
         TEMP1 = (ECR2(J2)-ECR1(J1))**2 + (ECZ2(J2)-ECZ1(J1))**2
         DO J=1,NP_EC
            TEMPA = (TEMP1 + 2.*ECR1(J1)*ECR2(J2)*(1.-ECPC(J)))**1.5
            TEMPB = ECPS(J)*ECPSN(J)/TEMPA
            TEMP3 = TEMP3 + TEMPB*(ECR1(J1)*ECDZDC(J1)+
     &                             (ECZ2(J2)-ECZ1(J1))*ECDRDC(J1))
            TEMP4 = TEMP4 + TEMPB*ECR2(J2)*ECDRDC(J1)
         ENDDO
         ENDDO
         ENDIF
         TEMP5 = TEMP5 - TEMP3*HCHI*RNTOR
         TEMP6 = TEMP6 + TEMP4*HCHI*RNTOR

         ECB1(J2)=ECB1(J2)+FEEDI(K)*(ECDSDR(J2)*TEMP5+ECDSDZ(J2)*TEMP6)
         ECBR(J2)=ECBR(J2)+FEEDI(K)*TEMP5
         ECBZ(J2)=ECBZ(J2)+FEEDI(K)*TEMP6
      ENDDO
      ENDDO

      ECB1 =-ECB1*HPHI/4./PI
      ECBR =-ECBR*HPHI/4./PI
      ECBZ =-ECBZ*HPHI/4./PI

      OPEN(CHOUTP,FILE='ECB1.OUT')
      REWIND(CHOUTP)
      DO J=1,NCHI
         WRITE(CHOUTP,112) ECC(J),REAL(ECB1(J)),IMAG(ECB1(J)),
     &                            REAL(ECBR(J)),IMAG(ECBR(J)),
     &                            REAL(ECBZ(J)),IMAG(ECBZ(J))
      ENDDO
      CLOSE(CHOUTP)
 112  FORMAT(7(1X,E15.8))

      DEALLOCATE(JK1,JK2)

      DEALLOCATE(ECC,ECR1,ECZ1,ECR2,ECZ2,
     &           ECDRDC,ECDZDC,ECDSDR,ECDSDZ,
     &           ECPC,ECPS,ECPCN,ECPSN)

      DEALLOCATE(ECB1,ECBR,ECBZ)

      RETURN
      END
      

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C GENERATE NEW COORDINATE SYSTEM (S,THETA,CHI), WITH   LIU YQ 13.02.2013
C THETA BEING THE GEOMETRIC POLOIDAL ANGLE, FOR THE PLASMA REGION.     $
C SHOULD BE CALLED RIGHT AFTER READING THE OUTRMAR DATA                $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE GEOMCS_PLASMA
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      IMPLICIT NONE
      INTEGER    I,J,K,NLEFT
      REAL*8     TQHCHI,TQH1,TQH2,TQDRDS,TQDRDT,TQDZDS,TQDZDT,R00,Z00
      REAL*8,DIMENSION(:),ALLOCATABLE::TQTET,TQCHI,TQTETN,TQCHIN,WORK
 
      ALLOCATE(TQTET(3*NCHI+1),TQCHI(3*NCHI+1),
     &         TQTETN(NCHI+1), TQCHIN(NCHI+1),
     &         WORK(3*NCHI+1))

      TQTET  = 0.
      TQCHI  = 0.
      TQTETN = 0.
      TQCHIN = 0.

C     DEFINE UNIFORM THETA-MESH IN [0,2*PI] 
      TQHCHI    = 2.0*PI/NCHI
      TQTETN(1) = 0.
      DO J=2,NCHI+1
         TQTETN(J) = TQTETN(J-1) + TQHCHI
      ENDDO

C     DEFINE (R00,Z00) AS AVERAGE OF (R(1,:),Z(1,:))
      R00 = 0.0
      Z00 = 0.0
      DO J=1,NCHI
         R00 = R00 + REQ(1,J)
         Z00 = Z00 + ZEQ(1,J)
      ENDDO
      R00 = R00/NCHI
      Z00 = Z00/NCHI

C     LOOPING THROUGH ALL INTEGER SURFACES
C     RE-COMPUTING (REQ,ZEQ)-MESH 
      NLEFT = 2
      IF (SLEFT.GT.0.) NLEFT=1

      DO I=NLEFT,NRP1

C     COMPUTE CORRESPONDING THETA-MESH, USING (REQ,ZEQ)
      DO J=1,NCHI
         TQTET(NCHI+J) = DATAN2(ZEQ(I,J)-Z00,REQ(I,J)-R00)
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

C     EXTEND THETA-MESH BY 2*PI IN BOTH DIRECTIONS
      DO J=1,NCHI
         TQTET(J) = TQTET(NCHI+J) - 2.*PI
         TQTET(2*NCHI+J+1) = TQTET(NCHI+J+1) + 2.*PI
      ENDDO

C     COMPUTE REQ IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = REQ(I,J)
         TQCHI(NCHI+J)   = REQ(I,J)
         TQCHI(2*NCHI+J) = REQ(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = REQ(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         REQ(I,J) = TQCHIN(J)
      ENDDO

C     COMPUTE ZEQ IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = ZEQ(I,J)
         TQCHI(NCHI+J)   = ZEQ(I,J)
         TQCHI(2*NCHI+J) = ZEQ(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = ZEQ(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         ZEQ(I,J) = TQCHIN(J)
      ENDDO

      ENDDO

C     LOOPING THROUGH ALL HALF-INTEGER SURFACES
C     RE-COMPUTING (REQM,ZEQM)-MESH 
      DO I=1,NR

C     COMPUTE CORRESPONDING THETA-MESH, USING (REQM,ZEQM)
      DO J=1,NCHI
         TQTET(NCHI+J) = DATAN2(ZEQM(I,J)-Z00,REQM(I,J)-R00)
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

C     EXTEND THETA-MESH BY 2*PI IN BOTH DIRECTIONS
      DO J=1,NCHI
         TQTET(J) = TQTET(NCHI+J) - 2.*PI
         TQTET(2*NCHI+J+1) = TQTET(NCHI+J+1) + 2.*PI
      ENDDO

C     COMPUTE REQM IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = REQM(I,J)
         TQCHI(NCHI+J)   = REQM(I,J)
         TQCHI(2*NCHI+J) = REQM(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = REQM(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         REQM(I,J) = TQCHIN(J)
      ENDDO

C     COMPUTE ZEQM IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = ZEQM(I,J)
         TQCHI(NCHI+J)   = ZEQM(I,J)
         TQCHI(2*NCHI+J) = ZEQM(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = ZEQM(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         ZEQM(I,J) = TQCHIN(J)
      ENDDO

      ENDDO

C     COMPUTE REMAINING 2D EQUILIBRIUM QUANTITIES AT INTEGER SURFACE
      DO I=2,NRP1
      TQH1 = (CSE(I)-CSE(I-1))*0.5     
      IF (I.LT.NRP1) TQH2 = (CSE(I+1)-CSE(I))*0.5
      DO J=1,NCHI
         IF (I.LT.NRP1) THEN
            TQDRDS = (TQH1**2*(REQM(I,J)-REQ(I,J)) + 
     &                TQH2**2*(REQ(I,J)-REQM(I-1,J)))/
     &               (TQH1*TQH2*(TQH1+TQH2))
            TQDZDS = (TQH1**2*(ZEQM(I,J)-ZEQ(I,J)) + 
     &                TQH2**2*(ZEQ(I,J)-ZEQM(I-1,J)))/
     &               (TQH1*TQH2*(TQH1+TQH2))
         ELSE
            TQDRDS = (REQ(I,J)-REQM(I-1,J))/TQH1
            TQDZDS = (ZEQ(I,J)-ZEQM(I-1,J))/TQH1
         ENDIF

         IF (J.EQ.1) THEN
            TQDRDT = (REQ(I,2)-REQ(I,NCHI))*0.5/TQHCHI
            TQDZDT = (ZEQ(I,2)-ZEQ(I,NCHI))*0.5/TQHCHI
         ELSEIF (J.EQ.NCHI) THEN
            TQDRDT = (REQ(I,1)-REQ(I,NCHI-1))*0.5/TQHCHI
            TQDZDT = (ZEQ(I,1)-ZEQ(I,NCHI-1))*0.5/TQHCHI
         ELSE
            TQDRDT = (REQ(I,J+1)-REQ(I,J-1))*0.5/TQHCHI
            TQDZDT = (ZEQ(I,J+1)-ZEQ(I,J-1))*0.5/TQHCHI
         ENDIF

         RJA(I,J)   = ABS(TQDRDS*TQDZDT-TQDRDT*TQDZDS)*REQ(I,J)
         G11L(I,J)  = TQDRDS**2 + TQDZDS**2
         G22L(I,J)  = TQDRDT**2 + TQDZDT**2
         G33L(I,J)  = REQ(I,J)**2
         G12L(I,J)  = TQDRDS*TQDRDT + TQDZDS*TQDZDT
         RDCDZ(I,J) = TQDRDS*REQ(I,J)/RJA(I,J)
         RDSDZ(I,J) =-TQDRDT*REQ(I,J)/RJA(I,J)
         RBZ(I,J)   = TQDZDT*DPSIDS(I)/RJA(I,J)
      ENDDO
      ENDDO

C     COMPUTE REMAINING 2D EQUILIBRIUM QUANTITIES AT HALF-INTEGER SURFACE
      DO I=1,NR
      TQH1 = CSE(I+1)-CSE(I)     
      DO J=1,NCHI
         TQDRDS = (REQ(I+1,J)-REQ(I,J))/TQH1
         TQDZDS = (ZEQ(I+1,J)-ZEQ(I,J))/TQH1

         IF (J.EQ.1) THEN
            TQDRDT = (REQM(I,2)-REQM(I,NCHI))*0.5/TQHCHI
            TQDZDT = (ZEQM(I,2)-ZEQM(I,NCHI))*0.5/TQHCHI
         ELSEIF (J.EQ.NCHI) THEN
            TQDRDT = (REQM(I,1)-REQM(I,NCHI-1))*0.5/TQHCHI
            TQDZDT = (ZEQM(I,1)-ZEQM(I,NCHI-1))*0.5/TQHCHI
         ELSE
            TQDRDT = (REQM(I,J+1)-REQM(I,J-1))*0.5/TQHCHI
            TQDZDT = (ZEQM(I,J+1)-ZEQM(I,J-1))*0.5/TQHCHI
         ENDIF

         RJAM(I,J)   = ABS(TQDRDS*TQDZDT-TQDRDT*TQDZDS)*REQM(I,J)
         G11LM(I,J)  = TQDRDS**2 + TQDZDS**2
         G22LM(I,J)  = TQDRDT**2 + TQDZDT**2
         G33LM(I,J)  = REQM(I,J)**2
         G12LM(I,J)  = TQDRDS*TQDRDT + TQDZDS*TQDZDT
         RDCDZM(I,J) = TQDRDS*REQM(I,J)/RJAM(I,J)
         RDSDZM(I,J) =-TQDRDT*REQM(I,J)/RJAM(I,J)
         RBZM(I,J)   = TQDZDT*DPSIDSM(I)/RJAM(I,J)
      ENDDO
      ENDDO
            
      DO J=1,NCHI
         RJA(1,J)   = 0.
         G11L(1,J)  = G11LM(1,J)
         G22L(1,J)  = 0.
         G33L(1,J)  = REQ(1,J)**2
         G12L(1,J)  = 0.
         RDCDZ(1,J) = RDCDZM(1,J)
         RDSDZ(1,J) = RDSDZM(1,J)
         RBZ(1,J)   = RBZM(1,J)
      ENDDO

      DEALLOCATE(TQTET,TQCHI,TQTETN,TQCHIN,WORK)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C GENERATE NEW COORDINATE SYSTEM (S,THETA,CHI), WITH   LIU YQ 13.02.2013
C THETA BEING THE GEOMETRIC POLOIDAL ANGLE, FOR THE VACUUM REGION.     $
C SHOULD BE CALLED RIGHT AFTER READING THE OUTVMAR DATA                $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE GEOMCS_VACUUM
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      USE GVACUUMM
      IMPLICIT NONE
      INTEGER    I,J,K
      REAL*8     TQHCHI,TQH1,TQH2,TQDRDS,TQDRDT,TQDZDS,TQDZDT,R00,Z00
      REAL*8,DIMENSION(:),ALLOCATABLE::TQTET,TQCHI,TQTETN,TQCHIN,WORK
 
      ALLOCATE(TQTET(3*NCHI+1),TQCHI(3*NCHI+1),
     &         TQTETN(NCHI+1), TQCHIN(NCHI+1),
     &         WORK(3*NCHI+1))

      TQTET  = 0.
      TQCHI  = 0.
      TQTETN = 0.
      TQCHIN = 0.

C     DEFINE UNIFORM THETA-MESH IN [0,2*PI] 
      TQHCHI    = 2.0*PI/NCHI
      TQTETN(1) = 0.
      DO J=2,NCHI+1
         TQTETN(J) = TQTETN(J-1) + TQHCHI
      ENDDO

C     DEFINE (R00,Z00) AS AVERAGE OF (R(1,:),Z(1,:))
      R00 = 0.0
      Z00 = 0.0
      DO J=1,NCHI
         R00 = R00 + REQ(1,J)
         Z00 = Z00 + ZEQ(1,J)
      ENDDO
      R00 = R00/NCHI
      Z00 = Z00/NCHI

C     LOOPING THROUGH ALL INTEGER SURFACES
C     RE-COMPUTING (VRR,VRZ)-MESH 
      DO I=1,NVEQ1

C     COMPUTE CORRESPONDING THETA-MESH, USING (VRR,VRZ)
      DO J=1,NCHI
         TQTET(NCHI+J) = DATAN2(VRZ(I,J)-Z00,VRR(I,J)-R00)
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

C     EXTEND THETA-MESH BY 2*PI IN BOTH DIRECTIONS
      DO J=1,NCHI
         TQTET(J) = TQTET(NCHI+J) - 2.*PI
         TQTET(2*NCHI+J+1) = TQTET(NCHI+J+1) + 2.*PI
      ENDDO

C     COMPUTE VRR IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = VRR(I,J)
         TQCHI(NCHI+J)   = VRR(I,J)
         TQCHI(2*NCHI+J) = VRR(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = VRR(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         VRR(I,J) = TQCHIN(J)
      ENDDO

C     COMPUTE VRZ IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = VRZ(I,J)
         TQCHI(NCHI+J)   = VRZ(I,J)
         TQCHI(2*NCHI+J) = VRZ(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = VRZ(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         VRZ(I,J) = TQCHIN(J)
      ENDDO

      DO K=1,NWALL
      IF (IWALL(K).EQ.I) THEN

C     COMPUTE ZCND IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = ZCND(K,J)
         TQCHI(NCHI+J)   = ZCND(K,J)
         TQCHI(2*NCHI+J) = ZCND(K,J)
      ENDDO
      TQCHI(3*NCHI+1)    = ZCND(K,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         ZCND(K,J) = TQCHIN(J)
      ENDDO

C     COMPUTE ZCNDC IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = ZCNDC(K,J)
         TQCHI(NCHI+J)   = ZCNDC(K,J)
         TQCHI(2*NCHI+J) = ZCNDC(K,J)
      ENDDO
      TQCHI(3*NCHI+1)    = ZCNDC(K,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         ZCNDC(K,J) = TQCHIN(J)
      ENDDO

      ENDIF
      ENDDO

      ENDDO

C     LOOPING THROUGH ALL HALF-INTEGER SURFACES
C     RE-COMPUTING (VRRM,VRZM)-MESH 
      DO I=1,NVEQ1-1

C     COMPUTE CORRESPONDING THETA-MESH, USING (REQM,ZEQM)
      DO J=1,NCHI
         TQTET(NCHI+J) = DATAN2(VRZM(I,J)-Z00,VRRM(I,J)-R00)
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

C     EXTEND THETA-MESH BY 2*PI IN BOTH DIRECTIONS
      DO J=1,NCHI
         TQTET(J) = TQTET(NCHI+J) - 2.*PI
         TQTET(2*NCHI+J+1) = TQTET(NCHI+J+1) + 2.*PI
      ENDDO

C     COMPUTE VRRM IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = VRRM(I,J)
         TQCHI(NCHI+J)   = VRRM(I,J)
         TQCHI(2*NCHI+J) = VRRM(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = VRRM(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         VRRM(I,J) = TQCHIN(J)
      ENDDO

C     COMPUTE VRZM IN UNIFORM THETA-MESH, USING SPLINE       
      DO J=1,NCHI
         TQCHI(J)        = VRZM(I,J)
         TQCHI(NCHI+J)   = VRZM(I,J)
         TQCHI(2*NCHI+J) = VRZM(I,J)
      ENDDO
      TQCHI(3*NCHI+1)    = VRZM(I,1)

      CALL SPLINE1D(TQCHIN,TQTETN,NCHI+1,TQCHI,TQTET,3*NCHI+1,WORK)

      DO J=1,NCHI
         VRZM(I,J) = TQCHIN(J)
      ENDDO

      ENDDO

C     COMPUTE REMAINING 2D EQUILIBRIUM QUANTITIES AT INTEGER SURFACE
      DO I=1,NVEQ1
      IF (I.GT.1)     TQH1 = (VCS(I)-VCS(I-1))*0.5     
      IF (I.LT.NVEQ1) TQH2 = (VCS(I+1)-VCS(I))*0.5
      DO J=1,NCHI
         IF (I.EQ.1) THEN
            TQDRDS = (VRRM(1,J)-VRR(1,J))/TQH2
            TQDZDS = (VRZM(1,J)-VRZ(1,J))/TQH2
         ELSEIF (I.EQ.NVEQ1) THEN
            TQDRDS = (VRR(I,J)-VRRM(I-1,J))/TQH1
            TQDZDS = (VRZ(I,J)-VRZM(I-1,J))/TQH1
         ELSE
            TQDRDS = (TQH1**2*(VRRM(I,J)-VRR(I,J)) + 
     &                TQH2**2*(VRR(I,J)-VRRM(I-1,J)))/
     &               (TQH1*TQH2*(TQH1+TQH2))
            TQDZDS = (TQH1**2*(VRZM(I,J)-VRZ(I,J)) + 
     &                TQH2**2*(VRZ(I,J)-VRZM(I-1,J)))/
     &               (TQH1*TQH2*(TQH1+TQH2))
         ENDIF

         IF (J.EQ.1) THEN
            TQDRDT = (VRR(I,2)-VRR(I,NCHI))*0.5/TQHCHI
            TQDZDT = (VRZ(I,2)-VRZ(I,NCHI))*0.5/TQHCHI
         ELSEIF (J.EQ.NCHI) THEN
            TQDRDT = (VRR(I,1)-VRR(I,NCHI-1))*0.5/TQHCHI
            TQDZDT = (VRZ(I,1)-VRZ(I,NCHI-1))*0.5/TQHCHI
         ELSE
            TQDRDT = (VRR(I,J+1)-VRR(I,J-1))*0.5/TQHCHI
            TQDZDT = (VRZ(I,J+1)-VRZ(I,J-1))*0.5/TQHCHI
         ENDIF

         VRJA(I,J)   = ABS(TQDRDS*TQDZDT-TQDRDT*TQDZDS)*VRR(I,J)
         VRG11L(I,J) = TQDRDS**2 + TQDZDS**2
         VRG22L(I,J) = TQDRDT**2 + TQDZDT**2
         VRG33L(I,J) = VRR(I,J)**2
         VRG12L(I,J) = TQDRDS*TQDRDT + TQDZDS*TQDZDT
      ENDDO
      ENDDO

C     COMPUTE REMAINING 2D EQUILIBRIUM QUANTITIES AT HALF-INTEGER SURFACE
      DO I=1,NVEQ1-1
      TQH1 = VCS(I+1)-VCS(I)     
      DO J=1,NCHI
         TQDRDS = (VRR(I+1,J)-VRR(I,J))/TQH1
         TQDZDS = (VRZ(I+1,J)-VRZ(I,J))/TQH1

         IF (J.EQ.1) THEN
            TQDRDT = (VRRM(I,2)-VRRM(I,NCHI))*0.5/TQHCHI
            TQDZDT = (VRZM(I,2)-VRZM(I,NCHI))*0.5/TQHCHI
         ELSEIF (J.EQ.NCHI) THEN
            TQDRDT = (VRRM(I,1)-VRRM(I,NCHI-1))*0.5/TQHCHI
            TQDZDT = (VRZM(I,1)-VRZM(I,NCHI-1))*0.5/TQHCHI
         ELSE
            TQDRDT = (VRRM(I,J+1)-VRRM(I,J-1))*0.5/TQHCHI
            TQDZDT = (VRZM(I,J+1)-VRZM(I,J-1))*0.5/TQHCHI
         ENDIF

         VRJAM(I,J)   = ABS(TQDRDS*TQDZDT-TQDRDT*TQDZDS)*VRRM(I,J)
         VRG11LM(I,J) = TQDRDS**2 + TQDZDS**2
         VRG22LM(I,J) = TQDRDT**2 + TQDZDT**2
         VRG33LM(I,J) = VRRM(I,J)**2
         VRG12LM(I,J) = TQDRDS*TQDRDT + TQDZDS*TQDZDT
      ENDDO
      ENDDO
            
      DEALLOCATE(TQTET,TQCHI,TQTETN,TQCHIN,WORK)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C READ IN EQUILIBRIUM FROM <FLOW> CODE                 LIU YQ 22.02.2013
C ADD VACUUM REGION TOO                                                $
C THE (R,Z) IS SPECIFIED IN GEOMETRIC ANGLE                            $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE EQUIL_FLOW
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
C     INTEGER    I,J,K
C     REAL*8     TQHCHI,TQH1,TQH2,TQDRDS,TQDRDT,TQDZDS,TQDZDT
C     REAL*8,DIMENSION(:),ALLOCATABLE::TQTET,TQCHI,TQTETN,TQCHIN,WORK
C
C     ALLOCATE(TQTET(3*NCHI+1),TQCHI(3*NCHI+1),
C    &         TQTETN(NCHI+1), TQCHIN(NCHI+1),
C    &         WORK(3*NCHI+1))

C     TQTET  = 0.
C     TQCHI  = 0.
C     TQTETN = 0.
C     TQCHIN = 0.

      INTEGER    NPSI,NTET,I,J,ISHOCK
      REAl*8,DIMENSION(:),  ALLOCATABLE::PSI,FF,GG,FLOW_TOR,FLOW_POL
      REAl*8,DIMENSION(:,:),ALLOCATABLE::RR,ZZ,RHOF  
 
C     NPSI  = NUMBER OF RADIAL POINTS ALONG POLOIDAL FLUX PSI
C     NTET  = NUMBER OF POINTS ALONG A UNIFORM GEOMETRIC POLOIDAL ANGLE, WITH ORIGIN LOCATED AT (RR(1,1),ZZ(1,1))
C     FLOW  = FLOW_POL(PSI)/RHO*B + R^2*FLOW_TOR(PSI)*GRAD(PHI)
C     FF    = FF(PSI) = FREE FUNCTION FOR B_PHI
C     GG    = GG(PSI) = P/RHO^GAMMA
C     RR    = R(PSI,THETA)
C     ZZ    = Z(PSI,THETA)
C     RHOF  = DENSITY AS FUNCTION OF (PSI,THETA)

      OPEN(CHOUTP,FILE='OUTRMAR')
      REWIND(CHOUTP)
      READ(CHOUTP,*) NPSI,NTET,ISHOCK
      IF (ISHOCK.GT.0) NPSI = NPSI + 1

      ALLOCATE( PSI(NPSI), FF(NPSI), GG(NPSI), 
     &          FLOW_TOR(NPSI), FLOW_POL(NPSI) )
      ALLOCATE( RR(NPSI,NTET), ZZ(NPSI,NTET),
     &          RHOF(NPSI,NTET) )

      DO I=1,NPSI
         READ(CHOUTP,*) PSI(I),FF(I),GG(I),FLOW_TOR(I),FLOW_POL(I)
      ENDDO

      DO I=1,NPSI
      DO J=1,NCHI
         READ(CHOUTP,*) RR(I,J),ZZ(I,J),RHOF(I,J)
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
C

      DEALLOCATE(PSI,FF,GG,FLOW_TOR,FLOW_POL,RR,ZZ,RHOF)
      RETURN
      END

      SUBROUTINE CREATE_WORKDIR
      USE MPIENV 
CADT: DELETE FOLLOWING LINE FOR NON MPI RUN AND PGF90
CADT: USE IFPORT
      USE GLOBALM
      IMPLICIT NONE
      LOGICAL(4) RES
         
      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
         IF (ODWKCOM) THEN
            IF (IPERTURB.EQ.0.AND.V2XKEY.NE.1.AND.V2XKEY.NE.3) THEN
               WRITE (*,*) 'ODWKCOM=.T. REQUIRES V2XKEY=1 OR 3'
               STOP 'ODWKCOM=.T. REQUIRES V2XKEY=1 OR 3'
            ENDIF
            WRITE (*,*) 'TURN ON THE DETAILED DWK ANALYSIS.'
            CALL SYSTEM('mkdir DATA_STORAGE')
         ENDIF
      ENDIF
      END SUBROUTINE CREATE_WORKDIR

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE NU* FOR THERMAL IONS AND ELECTRONS           LIU YQ 22.02.2013
C BETTER TO SET NPROFUI(E)=3, NUEFFI(E)A=1.0
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      SUBROUTINE CALC_NUSTAR
C     ======================
C
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE TORQUEM
      USE MPIENV
      IMPLICIT NONE
      INTEGER   I
      REAL*8    NUSI,NUSE,TMP

      IF (ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT) THEN
      
      OPEN(99,FILE='PROFNUSTAR.OUT')
      DO I=2,NR
         TMP  = SQRT(T1E(I)/RHO(I))
         NUSI = GNUI(I)*Q(I)*TMP/SQRT(2.0*ESPECIES_PRE(I,1,1))/
     &          (CS(I)/ASPCT/TMP)**1.5
         NUSE = GNUE(I)*Q(I)*TMP/SQRT(2.0*ESPECIES_PRE(I,1,2))/
     &          (CS(I)/ASPCT/TMP)**1.5/SQRT(ESPECIES_M(1)*1.8361e+3)
         WRITE(99,120) CS(I),NUSI,NUSE
      ENDDO
      CLOSE(99)
 120  FORMAT(3(E16.9,1X))

      ENDIF

      RETURN
      END SUBROUTINE CALC_NUSTAR

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE MAGNETIC FIELD IN (R,PHI,Z) SPACE            LIU YQ 04.01.2016
C THE FIELD IS PRODUCED BY A CURRENT FILAMENT PARALLEL TO THE          $
C EQUILIBRIUM MAGNETIC FIELD LINE INSIDE THE PLASMA.                   $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE OUTPUT_BFILAMENT
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      USE RCOMDM  
      USE GVACUUMM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INTEGER    I,J,K,L,KCHECK,JS
      REAL*8     HCHI,HPHI,OCHI1,OR1,OZ1,OCR,OCZ,OCP,OR3,
     &           ODPHI,ODPHIS,ODPHIC,OTMP
 
      REAL*8,DIMENSION(:),ALLOCATABLE::OPHI1,OPHI2,OQ1,ODR1DC,ODZ1DC
      REAL*8,DIMENSION(:,:),ALLOCATABLE::OR2,OZ2
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE::OBR,OBZ,OBP
 
      KCHECK = 0

C     STEP 1: GENERATE FLUX SURFACE BASED (R,Z)-GRID IN WHOLE DOMAIN   
C             USING HALF-INTEGER RADIAL POINTS
      ALLOCATE( OR2(NTOT,NCHI+1), OZ2(NTOT,NCHI+1), OPHI2(NOPP) )

      DO I=1,NR
         DO J=1,NCHI
            OR2(I,J) = REQM(I,J)
            OZ2(I,J) = ZEQM(I,J)
         ENDDO
         J = NCHI+1
         OR2(I,J) = REQM(I,1)
         OZ2(I,J) = ZEQM(I,1)
      ENDDO
      
      DO I=1,NV
         DO J=1,NCHI
            OR2(NR+I,J) = VRRM(I,J)
            OZ2(NR+I,J) = VRZM(I,J)
         ENDDO
         J = NCHI+1
         OR2(NR+I,J) = VRRM(I,1)
         OZ2(NR+I,J) = VRZM(I,1)
      ENDDO

      HPHI = 2.*PI/NOPP
      OPHI2(1) = 0.0
      DO K=2,NOPP
         OPHI2(K) = OPHI2(K-1) + HPHI
      ENDDO

C     STEP 2: COMPUTE OPHI1 ALONG EQUILIBRIUM MAGNETIC FIELD LINE
C             COMPUTE DR/DCHI AND DZ/DCHI AS WELL
      ALLOCATE( OPHI1(NCHI+1), OQ1(NCHI+1), ODR1DC(NCHI+1),
     &          ODZ1DC(NCHI+1) )
     
      I    = 10
      IF (JSOUT.GT.0) I    = JSOUT
      DO J=1,NCHI
         OQ1(J)   = T(I)*RJA(I,J)/DPSIDS(I)/REQ(I,J)**2
      ENDDO
      OQ1(NCHI+1) = T(I)*RJA(I,1)/DPSIDS(I)/REQ(I,1)**2
      
      HCHI = 2.*PI/NCHI
      DO J=1,NCHI+1
         OCHI1 = HCHI*(J-1)
         OPHI1(J) = OQ1(J)*OCHI1 + OCHI0
      ENDDO

      DO J=2,NCHI
         ODR1DC(J) = (REQ(I,J+1)-REQ(I,J-1))/HCHI/2.
         ODZ1DC(J) = (ZEQ(I,J+1)-ZEQ(I,J-1))/HCHI/2.
      ENDDO
      ODR1DC(1) = (REQ(I,2)-REQ(I,NCHI))/HCHI/2.
      ODZ1DC(1) = (ZEQ(I,2)-ZEQ(I,NCHI))/HCHI/2.
      ODR1DC(NCHI+1) = ODR1DC(1)
      ODZ1DC(NCHI+1) = ODZ1DC(1)

C     STEP 3: COMPUTE MAGNETIC FIELD PRODUCED BY CURRENT FILAMENT
      ALLOCATE( OBR(NTOT,NCHI+1,NOPP), OBZ(NTOT,NCHI+1,NOPP), 
     &          OBP(NTOT,NCHI+1,NOPP) )
     
      OBR = 0.
      OBZ = 0.
      OBP = 0.

      JS = 10
      IF (JSOUT.GT.0) JS = JSOUT

      DO I=1,NTOT
      DO J=1,NCHI+1
      DO K=1,NOPP
         DO L=1,NCHI+1
            IF (L.LE.NCHI) THEN
               OR1  = REQ(JS,L)
               OZ1  = ZEQ(JS,L)
            ELSE
               OR1  = REQ(JS,1)
               OZ1  = ZEQ(JS,1)
            ENDIF
            
            ODPHI  = OPHI1(L) - OPHI2(K)
            ODPHIS = SIN(ODPHI)
            ODPHIC = COS(ODPHI)
            OTMP   = (ODZ1DC(L)+ODR1DC(L)*(OZ2(I,J)-OZ1)/OR1)/OQ1(L)
            OCR    = OTMP*ODPHIS + (OZ2(I,J)-OZ1)*ODPHIC
            OCZ    = OR1 - OR2(I,J)*ODR1DC(L)/OQ1(L)/OR1*ODPHIS -
     &               OR2(I,J)*ODPHIC
            OCP    = OR2(I,J)*ODZ1DC(L)/OQ1(L)/OR1 - OTMP*ODPHIC + 
     &               (OZ2(I,J)-OZ1)*ODPHIS
            OR3    = (OR1**2+OR2(I,J)**2-2*OR1*OR2(I,J)*ODPHIC+
     &                (OZ2(I,J)-OZ1)**2)**1.5
        
            OBR(I,J,K) = OBR(I,J,K) + OCR/OR3
            OBZ(I,J,K) = OBZ(I,J,K) + OCZ/OR3
            OBP(I,J,K) = OBP(I,J,K) + OCP/OR3
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      OTMP = B0EXP*ODJPHI/4./PI*HCHI
      OBR  = OBR*OTMP
      OBZ  = OBZ*OTMP
      OBP  = OBP*OTMP

C     STEP 4: SAVE FIELD INTO A FILE
      OPEN(CHOUTP,FILE="BFILAMENT.OUT")
      REWIND(CHOUTP)
      WRITE(CHOUTP,*) NTOT,NCHI+1,NOPP,NTOT,NCHI+1,NOPP
      DO I=1,NTOT
      DO J=1,NCHI+1
      DO K=1,NOPP
         WRITE(CHOUTP,120) OR2(I,J)*R0EXP,OZ2(I,J)*R0EXP,OPHI2(K),
     &                     OBR(I,J,K),OBZ(I,J,K),OBP(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(6(E11.4,1X))

C     DEALOCATE ARRAYS
      DEALLOCATE( OPHI1,OPHI2,OQ1,ODR1DC,ODZ1DC,OR2,OZ2,OBR,OBZ,OBP )

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE THE FOLLOWING VACUUM METRICS ELEMENTS AT LOCATIONS SPECIFIED $
C BY ARRAYS IPICK AND CPICK IN THE (S,CHI) SPACE:                      $
C    VZSJ=(DZ/DS)/J                                                    $
C    VZCJ=(DZ/DCHI)/J                                                  $
C    VRSJ=(DR/DS)/J                                                    $
C    VRCJ=(DR/DCHI)/J                                                  $
C    VR3J=R/J                                                          $
C Y.Q.LIU, 02/03/2016                                                  $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CALC_PICK_DATA
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE FEEDBACKM
      USE RCOMDM  
      USE GVACUUMM
      IMPLICIT NONE
      INTEGER    I,J,K,KCHECK
      REAL*8,DIMENSION(:),ALLOCATABLE::TMPCHI,TMPCHI2,TMPZ
      REAL*8 OZ1,OZ2,OZ3,OZ4,OJ0,OCHIX,OCHI1,OCHI2,OCHIH,OCHI
      REAL*8 OR1,OR2,OR3,OR4

      KCHECK = 0

      OCHIH = 1.0E-3

      ALLOCATE(TMPCHI(NCHI+1),TMPCHI2(NCHI+1),TMPZ(NCHI+1))
      IF (.NOT.ALLOCATED(VZSJ)) ALLOCATE(VZSJ(NPICK),VZCJ(NPICK))
      IF (.NOT.ALLOCATED(VRSJ)) ALLOCATE(VRSJ(NPICK),VRCJ(NPICK))
      IF (.NOT.ALLOCATED(VR3J)) ALLOCATE(VR3J(NPICK))

      VZSJ = 0.
      VZCJ = 0.
      VRSJ = 0.
      VRCJ = 0.
      VR3J = 0.
      
      TMPCHI  = 0.
      TMPCHI2 = 0.
      TMPZ    = 0.

      OCHI      = 2./NCHI
      DO J=2,NCHI+1
         TMPCHI(J) = TMPCHI(J-1) + OCHI
      ENDDO
      TMPCHI2 = 0.

      DO K=1,NPICK
         I = IPICK(K)  
         IF (I.LE.0.OR.I.GE.NV) STOP 'IPICK INCORRECT'
         OCHIX = CPICK(K)
         OCHI1 = OCHIX - OCHIH
         OCHI2 = OCHIX + OCHIH

         TMPZ(1:NCHI) = VRJAM(I,1:NCHI)
         TMPZ(NCHI+1) = VRJAM(I,1)
         OCHI = OCHIX
         IF (OCHI.LT.0.) OCHI = OCHI + 2.
         CALL SPLINE1D(OJ0,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         TMPZ(1:NCHI) = VRZ(I,1:NCHI)
         TMPZ(NCHI+1) = VRZ(I,1)
         CALL SPLINE1D(OZ1,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         TMPZ(1:NCHI) = VRZ(I+1,1:NCHI)
         TMPZ(NCHI+1) = VRZ(I+1,1)
         CALL SPLINE1D(OZ2,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         TMPZ(1:NCHI) = VRZM(I,1:NCHI)
         TMPZ(NCHI+1) = VRZM(I,1)
         OCHI = OCHI1
         IF (OCHI.LT.0.) OCHI = OCHI + 2.
         CALL SPLINE1D(OZ3,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)
         
         OCHI = OCHI2
         IF (OCHI.LT.0.) OCHI = OCHI + 2.
         IF (OCHI.GT.2.) OCHI = OCHI - 2.
         CALL SPLINE1D(OZ4,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         VZSJ(K) = (OZ2-OZ1)/VCSH(I)/OJ0
         VZCJ(K) = (OZ4-OZ3)/OCHIH/2./PI/OJ0
C ===
         TMPZ(1:NCHI) = VRR(I,1:NCHI)
         TMPZ(NCHI+1) = VRR(I,1)
         CALL SPLINE1D(OR1,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         TMPZ(1:NCHI) = VRR(I+1,1:NCHI)
         TMPZ(NCHI+1) = VRR(I+1,1)
         CALL SPLINE1D(OR2,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         TMPZ(1:NCHI) = VRRM(I,1:NCHI)
         TMPZ(NCHI+1) = VRRM(I,1)
         OCHI = OCHI1
         IF (OCHI.LT.0.) OCHI = OCHI + 2.
         CALL SPLINE1D(OR3,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         OCHI = OCHI2
         IF (OCHI.LT.0.) OCHI = OCHI + 2.
         IF (OCHI.GT.2.) OCHI = OCHI - 2.
         CALL SPLINE1D(OR4,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         VRSJ(K) = (OR2-OR1)/VCSH(I)/OJ0
         VRCJ(K) = (OR4-OR3)/OCHIH/2./PI/OJ0
C ===
         TMPZ(1:NCHI) = VRRM(I,1:NCHI)
         TMPZ(NCHI+1) = VRRM(I,1)
         CALL SPLINE1D(OR1,OCHI,1,TMPZ,TMPCHI,NCHI+1,TMPCHI2)

         VR3J(K) = OR1/OJ0
      ENDDO
         
      DEALLOCATE( TMPCHI,TMPCHI2,TMPZ )

      IF (KCHECK.EQ.1) THEN
         WRITE(*,*) 'CALC_PICK_DATA:'
         DO K=1,NPICK
            WRITE(*,100) K,IPICK(K),CPICK(K),VRSJ(K),VRCJ(K),VZSJ(K),
     &                   VZCJ(K),VR3J(K)
         ENDDO
 100  FORMAT(I3,1X,I3,1X,10E12.4)
      ENDIF

      RETURN
      END

