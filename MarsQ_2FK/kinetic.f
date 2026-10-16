C=======================================================================
C FILL IN MATRICES DUE TO KINETIC PRESSURE TERMS                      
C IMPORTANT NOTE: MPI VERSION IS ONLY VALID FOR ONE MARS-K RUN        
C                 TO SCAN PARAMETER, AN EXTERNAL PROGRAM NEEDS TO BE  
C                 USED                                                
C
C A GENERIC APPROACH TO INCLUDE DRIFT KINETIC EFFECTS INTO MARS-K:    
C ASSUMING AN ARRAY OF PARTICLE SPECIES FOR BOTH EQUILIBRIUM AND      
C PERTURBATIONS, WITH EACH SPECIES HAS ONE OF THE FOLLOWING POSSIBLE  
C TYPES OF EQUILIBRIUM DISTRIBUTION FUNCTIONS:                        
C                                                                     
C THE FOLLOWING QUANTITIES (2D-MATRICES) SPECIFY EQUILIBRIUM RADIAL   
C PROFILES FOR EACH PARTICLE SPECIES. THESE ARE NOT NAMELIST VARIABLES
C                                                                                       
C ESPECIES_DENF = {FRACTION OF DENSITY TO TOTAL ELECTRON DENSITY}     
C ESPECIES_PREF = {FRACTION OF PRESSURE TO THERMAL PRESSURE}          
C ESPECIES_DEN  = {(SURFACE-AVERAGED) RADIAL PROFILE OF DENSITY}      
C EPSECIES_PRE  = {(SURFACE-AVERAGED) RADIAL PROFILE OF PRESSURE}     
C ESPECIES_TEM  = {RADIAL PROFILE OF TEMPERATURE}                     
C EPSECIES_PREP = {(SURFACE-AVERAGED) RADIAL PROFILE OF DP/DPSI}      
C
C ANISOTROPIC NBI MODEL (IF0TYPE=3) IMPLEMENTED BY G.Z.HAO             
C OPEN MP AND MPI IMPLEMENTED BY Z.R.WANG                             
C YQL, 07-2013                                                        
C=======================================================================
      SUBROUTINE KJP(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,
     &               GSUBM,HSUBM,SHIFTC,SHIFTM,
     &               MXMAX_DUMMY,MYMAX_DUMMY,NRP1_DUMMY)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE OMP_LIB
      USE ToolBox
      USE MPIENV
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'compam.inc'
      
      INTEGER MROW,MSA,JS,J,I,LXROW,LYROW,LXCOL,LYCOL
      REAL*8  RTMP,ZV2M,ZV2P,ZB3M,ZB3P
      COMPLEX*16 CTMP,px,plasmax,py,plasmay,pz,plasmaz
      
      INTEGER   IEXV2,IEXB3            
      PARAMETER (IEXV2=-1, IEXB3=1)

      INTEGER KCHECK
      
      INTEGER MXMAX_DUMMY,MYMAX_DUMMY,NRP1_DUMMY
      COMPLEX*16       ASUBM(MXMAX_DUMMY,MXMAX_DUMMY,*),
     &                 BSUBM(MXMAX_DUMMY,MXMAX_DUMMY,*),
     &                 CSUBM(MXMAX_DUMMY,MXMAX_DUMMY,*),
     &                 DSUBM(MYMAX_DUMMY,MYMAX_DUMMY,*),
     &                 ESUBM(MXMAX_DUMMY,MYMAX_DUMMY,*),
     &                 FSUBM(MYMAX_DUMMY,MXMAX_DUMMY,*),
     &                 GSUBM(MYMAX_DUMMY,MXMAX_DUMMY,*),
     &                 HSUBM(MXMAX_DUMMY,MYMAX_DUMMY,*)
      COMPLEX*16,DIMENSION(NRP1_DUMMY)::SHIFTC,SHIFTM
      
C     THE VARIABLES FOR MPI COMPUTATION      
      INTEGER BUFFER_SEND(3),BUFFER_REC(3)
      COMPLEX*16 BUFFER_EIGENVALUE
      INTEGER TAG(4),BUFFERSIZE,COUNTER,INDEX,SLAVENUM
      INTEGER TOTTASK,TASKNUM,MPIJS,MPIKGRID
      COMPLEX*16,ALLOCATABLE:: BUFFER_DATA(:,:,:)
      INTEGER,ALLOCATABLE:: TASKQUEUE(:,:)
      INTEGER STAT(MPI_STATUS_SIZE) 

C     THE FILE CHANEL FOR OUTPUT KINETIC QUANTITY IN 2D PLOT
      INTEGER SURFFILENUM

      INCLUDE 'integc.inc'

      KCHECK=0

      IF (KJPKEY.EQ.0) GOTO 1
      IF (KJPKEY.EQ.2) GOTO 3
      IF ((ISMPIRUN.EQ.1.OR.ISMPIRUN.EQ.3).AND.RANK.NE.ROOT) THEN
          ALLOCATE( 
     &          VX1PARA(MSMAX,MSMAX,1), VX1PARAM(MSMAX,MSMAX,1),
     &          VX1PERP(MSMAX,MSMAX,1), VX1PERPM(MSMAX,MSMAX,1),
     &          VX2PARA(MSMAX,MSMAX,1), VX2PARAM(MSMAX,MSMAX,1),
     &          VX2PERP(MSMAX,MSMAX,1), VX2PERPM(MSMAX,MSMAX,1),
     &          VQ1PARA(MSMAX,MSMAX,1), VQ1PARAM(MSMAX,MSMAX,1),
     &          VQ1PERP(MSMAX,MSMAX,1), VQ1PERPM(MSMAX,MSMAX,1),
     &          VQ2PARA(MSMAX,MSMAX,1), VQ2PARAM(MSMAX,MSMAX,1),
     &          VQ2PERP(MSMAX,MSMAX,1), VQ2PERPM(MSMAX,MSMAX,1),
     &          VQ3PARA(MSMAX,MSMAX,1), VQ3PARAM(MSMAX,MSMAX,1),
     &          VQ3PERP(MSMAX,MSMAX,1), VQ3PERPM(MSMAX,MSMAX,1),
     &          VDPPARA(MSMAX,MSMAX,1), VDPPARAM(MSMAX,MSMAX,1),
     &          VDPPERP(MSMAX,MSMAX,1), VDPPERPM(MSMAX,MSMAX,1) )

          VX1PARA  = 0.     
          VX1PERP  = 0.     
          VX2PARA  = 0.     
          VX2PERP  = 0.     
          VQ1PARA  = 0.     
          VQ1PERP  = 0.     
          VQ2PARA  = 0.     
          VQ2PERP  = 0.     
          VQ3PARA  = 0.     
          VQ3PERP  = 0.     
          VDPPARA  = 0.     
          VDPPERP  = 0.     

          VX1PARAM = 0.     
          VX1PERPM = 0.     
          VX2PARAM = 0.     
          VX2PERPM = 0.     
          VQ1PARAM = 0.     
          VQ1PERPM = 0.     
          VQ2PARAM = 0.     
          VQ2PERPM = 0.     
          VQ3PARAM = 0.     
          VQ3PERPM = 0.     
          VDPPARAM = 0.     
          VDPPERPM = 0.     

          IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &        .AND.INCKIN.GT.0) THEN
          ALLOCATE( 
     &          VX1PARE(MSMAX,MSMAX,1), VX1PAREM(MSMAX,MSMAX,1),
     &          VX1PERE(MSMAX,MSMAX,1), VX1PEREM(MSMAX,MSMAX,1),
     &          VX2PARE(MSMAX,MSMAX,1), VX2PAREM(MSMAX,MSMAX,1),
     &          VX2PERE(MSMAX,MSMAX,1), VX2PEREM(MSMAX,MSMAX,1),
     &          VQ1PARE(MSMAX,MSMAX,1), VQ1PAREM(MSMAX,MSMAX,1),
     &          VQ1PERE(MSMAX,MSMAX,1), VQ1PEREM(MSMAX,MSMAX,1),
     &          VQ2PARE(MSMAX,MSMAX,1), VQ2PAREM(MSMAX,MSMAX,1),
     &          VQ2PERE(MSMAX,MSMAX,1), VQ2PEREM(MSMAX,MSMAX,1),
     &          VQ3PARE(MSMAX,MSMAX,1), VQ3PAREM(MSMAX,MSMAX,1),
     &          VQ3PERE(MSMAX,MSMAX,1), VQ3PEREM(MSMAX,MSMAX,1), 
     &          VDPPARE(MSMAX,MSMAX,1), VDPPAREM(MSMAX,MSMAX,1),
     &          VDPPERE(MSMAX,MSMAX,1), VDPPEREM(MSMAX,MSMAX,1) )

          VX1PARE  = 0.     
          VX1PERE  = 0.     
          VX2PARE  = 0.     
          VX2PERE  = 0.     
          VQ1PARE  = 0.     
          VQ1PERE  = 0.     
          VQ2PARE  = 0.     
          VQ2PERE  = 0.     
          VQ3PARE  = 0.     
          VQ3PERE  = 0.     
          VDPPARE  = 0.     
          VDPPERE  = 0.     

          VX1PAREM = 0.     
          VX1PEREM = 0.     
          VX2PAREM = 0.     
          VX2PEREM = 0.     
          VQ1PAREM = 0.     
          VQ1PEREM = 0.     
          VQ2PAREM = 0.     
          VQ2PEREM = 0.     
          VQ3PAREM = 0.     
          VQ3PEREM = 0.     
          VDPPAREM = 0.     
          VDPPEREM = 0.     
          ENDIF

          IF (INCDPHI.GT.0) THEN
             ALLOCATE( 
     &          VX1DPHI(MSMAX,MSMAX,1), VX1DPHIM(MSMAX,MSMAX,1),
     &          VX2DPHI(MSMAX,MSMAX,1), VX2DPHIM(MSMAX,MSMAX,1),
     &          VQ1DPHI(MSMAX,MSMAX,1), VQ1DPHIM(MSMAX,MSMAX,1),
     &          VQ2DPHI(MSMAX,MSMAX,1), VQ2DPHIM(MSMAX,MSMAX,1),
     &          VQ3DPHI(MSMAX,MSMAX,1), VQ3DPHIM(MSMAX,MSMAX,1),
     &          VDPDPHI(MSMAX,MSMAX,1), VDPDPHIM(MSMAX,MSMAX,1) )

             VX1DPHI  = 0.
             VX2DPHI  = 0.
             VQ1DPHI  = 0.
             VQ2DPHI  = 0.
             VQ3DPHI  = 0.
             VDPDPHI  = 0.

             VX1DPHIM = 0.
             VX2DPHIM = 0.
             VQ1DPHIM = 0.
             VQ2DPHIM = 0.
             VQ3DPHIM = 0.
             VDPDPHIM = 0.
          ENDIF

      ELSE
          ALLOCATE( 
     &          VX1PARA(MSMAX,MSMAX,NRP1), VX1PARAM(MSMAX,MSMAX,NR),
     &          VX1PERP(MSMAX,MSMAX,NRP1), VX1PERPM(MSMAX,MSMAX,NR),
     &          VX2PARA(MSMAX,MSMAX,NRP1), VX2PARAM(MSMAX,MSMAX,NR),
     &          VX2PERP(MSMAX,MSMAX,NRP1), VX2PERPM(MSMAX,MSMAX,NR),
     &          VQ1PARA(MSMAX,MSMAX,NRP1), VQ1PARAM(MSMAX,MSMAX,NR),
     &          VQ1PERP(MSMAX,MSMAX,NRP1), VQ1PERPM(MSMAX,MSMAX,NR),
     &          VQ2PARA(MSMAX,MSMAX,NRP1), VQ2PARAM(MSMAX,MSMAX,NR),
     &          VQ2PERP(MSMAX,MSMAX,NRP1), VQ2PERPM(MSMAX,MSMAX,NR),
     &          VQ3PARA(MSMAX,MSMAX,NRP1), VQ3PARAM(MSMAX,MSMAX,NR),
     &          VQ3PERP(MSMAX,MSMAX,NRP1), VQ3PERPM(MSMAX,MSMAX,NR),
     &          VDPPARA(MSMAX,MSMAX,NRP1), VDPPARAM(MSMAX,MSMAX,NR),
     &          VDPPERP(MSMAX,MSMAX,NRP1), VDPPERPM(MSMAX,MSMAX,NR) )

          VX1PARA  = 0.     
          VX1PERP  = 0.     
          VX2PARA  = 0.     
          VX2PERP  = 0.     
          VQ1PARA  = 0.     
          VQ1PERP  = 0.     
          VQ2PARA  = 0.     
          VQ2PERP  = 0.     
          VQ3PARA  = 0.     
          VQ3PERP  = 0.     
          VDPPARA  = 0.     
          VDPPERP  = 0.     

          VX1PARAM = 0.     
          VX1PERPM = 0.     
          VX2PARAM = 0.     
          VX2PERPM = 0.     
          VQ1PARAM = 0.     
          VQ1PERPM = 0.     
          VQ2PARAM = 0.     
          VQ2PERPM = 0.     
          VQ3PARAM = 0.     
          VQ3PERPM = 0.     
          VDPPARAM = 0.     
          VDPPERPM = 0.     

          IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &        .AND.INCKIN.GT.0) THEN
          ALLOCATE( 
     &          VX1PARE(MSMAX,MSMAX,NRP1), VX1PAREM(MSMAX,MSMAX,NR),
     &          VX1PERE(MSMAX,MSMAX,NRP1), VX1PEREM(MSMAX,MSMAX,NR),
     &          VX2PARE(MSMAX,MSMAX,NRP1), VX2PAREM(MSMAX,MSMAX,NR),
     &          VX2PERE(MSMAX,MSMAX,NRP1), VX2PEREM(MSMAX,MSMAX,NR),
     &          VQ1PARE(MSMAX,MSMAX,NRP1), VQ1PAREM(MSMAX,MSMAX,NR),
     &          VQ1PERE(MSMAX,MSMAX,NRP1), VQ1PEREM(MSMAX,MSMAX,NR),
     &          VQ2PARE(MSMAX,MSMAX,NRP1), VQ2PAREM(MSMAX,MSMAX,NR),
     &          VQ2PERE(MSMAX,MSMAX,NRP1), VQ2PEREM(MSMAX,MSMAX,NR),
     &          VQ3PARE(MSMAX,MSMAX,NRP1), VQ3PAREM(MSMAX,MSMAX,NR),
     &          VQ3PERE(MSMAX,MSMAX,NRP1), VQ3PEREM(MSMAX,MSMAX,NR),
     &          VDPPARE(MSMAX,MSMAX,NRP1), VDPPAREM(MSMAX,MSMAX,NR),
     &          VDPPERE(MSMAX,MSMAX,NRP1), VDPPEREM(MSMAX,MSMAX,NR) )

          VX1PARE  = 0.     
          VX1PERE  = 0.     
          VX2PARE  = 0.     
          VX2PERE  = 0.     
          VQ1PARE  = 0.     
          VQ1PERE  = 0.     
          VQ2PARE  = 0.     
          VQ2PERE  = 0.     
          VQ3PARE  = 0.     
          VQ3PERE  = 0.     
          VDPPARE  = 0.     
          VDPPERE  = 0.     

          VX1PAREM = 0.     
          VX1PEREM = 0.     
          VX2PAREM = 0.     
          VX2PEREM = 0.     
          VQ1PAREM = 0.     
          VQ1PEREM = 0.     
          VQ2PAREM = 0.     
          VQ2PEREM = 0.     
          VQ3PAREM = 0.     
          VQ3PEREM = 0.     
          VDPPAREM = 0.     
          VDPPEREM = 0.     
          ENDIF

          IF (INCDPHI.GT.0) THEN
             ALLOCATE( 
     &          VX1DPHI(MSMAX,MSMAX,NRP1), VX1DPHIM(MSMAX,MSMAX,NRP1),
     &          VX2DPHI(MSMAX,MSMAX,NRP1), VX2DPHIM(MSMAX,MSMAX,NRP1),
     &          VQ1DPHI(MSMAX,MSMAX,NRP1), VQ1DPHIM(MSMAX,MSMAX,NRP1),
     &          VQ2DPHI(MSMAX,MSMAX,NRP1), VQ2DPHIM(MSMAX,MSMAX,NRP1),
     &          VQ3DPHI(MSMAX,MSMAX,NRP1), VQ3DPHIM(MSMAX,MSMAX,NRP1),
     &          VDPDPHI(MSMAX,MSMAX,NRP1), VDPDPHIM(MSMAX,MSMAX,NRP1) )

             VX1DPHI  = 0.
             VX2DPHI  = 0.
             VQ1DPHI  = 0.
             VQ2DPHI  = 0.
             VQ3DPHI  = 0.
             VDPDPHI  = 0.

             VX1DPHIM = 0.
             VX2DPHIM = 0.
             VQ1DPHIM = 0.
             VQ2DPHIM = 0.
             VQ3DPHIM = 0.
             VDPDPHIM = 0.
         ENDIF
      ENDIF

      CALL ALLOCATEDWKCOMPMAT
      
      ALLOCATE( LAMM(2*NLAMK+2), LAMHH(2*NLAMK), LAMTMP(2*NLAMK+2) )

      ALLOCATE( RCHI(NCHI+1), RCHI2(NCHI+1), RW1(NCHI+1),
     &          RW2(NCHI+1),
     &          RJB(NCHI+1),  RX1P(NCHI+1),  RX1B(NCHI+1), 
     &          RX1R(NCHI+1), RX2(NCHI+1),   RQ1(NCHI+1),   
     &          RQ2(NCHI+1),  RBT(NCHI+1),   RPHI(NCHI+1),  
     &          RDMU(NCHI+1), RDB(NCHI+1) )

      ALLOCATE( NLAMK0(NRP1,2), NLAMK1(NRP1,2) )

      ALLOCATE( HKMIN(NRP1,2), HKMAX(NRP1,2), WFUN(NRP1,2) )

      ALLOCATE( BK(NRP1,NCHI,2),      HK(NRP1,NCHI,2),  
     &          BPK(NRP1,NCHI,2),    
     &          LAMK0(NRP1,NLAMK,2), LAMK1(NRP1,NLAMK,2) )

      ALLOCATE( KNUMDISTRIB(NSPECIES) ) 
      ALLOCATE( NUEFF    (NRP1,2,NSPECIES), 
     &          OMEGASN  (NRP1,2,NSPECIES), 
     &          OMEGAST  (NRP1,2,NSPECIES), 
     &          OMEGASNA (NRP1,2,NSPECIES), 
     &          OMEGASAA (NRP1,2,NSPECIES),
     &          OMEGASCA (NRP1,2,NSPECIES), 
     &          ALPHAA1  (NRP1,2,NSPECIES),
     &          ALPHAA2  (NRP1,2,NSPECIES),  
     &          ALPHAA3  (NRP1,2,NSPECIES),
     &          EPSLONCA (NRP1,2,NSPECIES),
     &          EPSALPHA (NRP1,2,NSPECIES),
     &          AAK      (NRP1,2,NSPECIES),
     &          ZC1      (NRP1,2,NSPECIES) )
   
      ALLOCATE( FREQK(NRP1,13) ) 
 
      IF (ISMPIRUN.EQ.0) THEN
      ALLOCATE( FREQKSURF(NRP1,7,NCHI),
     &           POSITIONKAI(NRP1,NCHI) )
      ENDIF

C     NUMBER OF INTEGRATION POINTS ALONG CHI
      NCHIT = 2*NCHI0

      ALLOCATE( RCHIK(NCHIT+2), RPHIK(NCHIT+2), RTK(NCHIT+2),
     &          RHK(NCHIT+2),   RJBK(NCHIT+2),  RX1PK(NCHIT+2),
     &          RX1BK(NCHIT+2), RX1RK(NCHIT+2), RX2K(NCHIT+2),  
     &          RQ1K(NCHIT+2),  RQ2K(NCHIT+2),  RVALK(NCHIT+2) )

      ALLOCATE( RVAK1(NCHIT+2), RVAK2(NCHIT+2), RVAK3(NCHIT+2), 
     &          RVAK4(NCHIT+2), RVAK5(NCHIT+2), RVAK6(NCHIT+2),
     &          RVAK7(NCHIT+2) )

      ALLOCATE( RCHIN(2*NCHIT+2), RVALN(2*NCHIT+2), RHN(2*NCHIT+2) )

      CALL ZALLOCANISO(1)

  
C     SET A SURFACE FOR OUTPUT
C     READ FROM NAMELIST 
      JS0 = 10
      IF (JSOUT.GT.0) JS0 = JSOUT

C     SET INTEGRATION WEIGHTING PARAMETER
      WK = 1.0/SQRT(3.0)

C     SET MESH FOR NUMERICAL PARTICLE ENERGY INTEGRATION
      CALL KENERGYMESH

C     CALCULATE RLM 
      CALL KGETRLM

      ALLOCATE( VPARA(MSMAX,MLMAX), VPERP(MSMAX,MLMAX),
     &          VDPHI(MSMAX,MLMAX), 
     &          VX1(MSMAX,MLMAX),   VX2(MSMAX,MLMAX), 
     &          VQ1(MSMAX,MLMAX),   VQ2(MSMAX,MLMAX), 
     &          VQ3(MSMAX,MLMAX),   VDP(MSMAX,MLMAX),  
     &          VI(4,MLMAX,NSPECIES) ) 

      ALLOCATE( VPARA0(MSMAX,MLMAX), VPERP0(MSMAX,MLMAX),
     &          VDPHI0(MSMAX,MLMAX),
     &          VX10(MSMAX,MLMAX),   VX20(MSMAX,MLMAX), 
     &          VQ10(MSMAX,MLMAX),   VQ20(MSMAX,MLMAX), 
     &          VQ30(MSMAX,MLMAX),   VDP0(MSMAX,MLMAX) )
      ALLOCATE( SVPARA0(MSMAX,MLMAX,NSPECIES),
     &          SVPERP0(MSMAX,MLMAX,NSPECIES), 
     &          SVDPHI0(MSMAX,MLMAX,NSPECIES), 
     &          SVX10  (MSMAX,MLMAX,NSPECIES),
     &          SVX20  (MSMAX,MLMAX,NSPECIES), 
     &          SVQ10  (MSMAX,MLMAX,NSPECIES),
     &          SVQ20  (MSMAX,MLMAX,NSPECIES),
     &          SVQ30  (MSMAX,MLMAX,NSPECIES),
     &          SVDP0  (MSMAX,MLMAX,NSPECIES),
     &          SLAM0  (MLMAX,NSPECIES),       
     &          SF0    (MLMAX,NSPECIES,0:3,4) )

C     IF (IFOWP.EQ.1.OR.IFOWT.EQ.1) THEN
      ALLOCATE( VI1(4,MLMAX,NSPECIES), 
     &          VI2(4,MLMAX,NSPECIES), 
     &          VI3(4,MLMAX,NSPECIES) )
      ALLOCATE( VPARA1(MSMAX,MLMAX), VPERP1(MSMAX,MLMAX),
     &          VDPHI1(MSMAX,MLMAX),
     &          VX11(MSMAX,MLMAX),   VX21(MSMAX,MLMAX), 
     &          VQ11(MSMAX,MLMAX),   VQ21(MSMAX,MLMAX), 
     &          VQ31(MSMAX,MLMAX),   VDP1(MSMAX,MLMAX) )

      ALLOCATE( VPARA01(MSMAX,MLMAX), VPERP01(MSMAX,MLMAX),
     &          VDPHI01(MSMAX,MLMAX),
     &          VX101(MSMAX,MLMAX),   VX201(MSMAX,MLMAX), 
     &          VQ101(MSMAX,MLMAX),   VQ201(MSMAX,MLMAX), 
     &          VQ301(MSMAX,MLMAX),   VDP01(MSMAX,MLMAX) )

      ALLOCATE( SVPARA01(MSMAX,MLMAX,NSPECIES),
     &          SVPERP01(MSMAX,MLMAX,NSPECIES), 
     &          SVDPHI01(MSMAX,MLMAX,NSPECIES), 
     &          SVX101  (MSMAX,MLMAX,NSPECIES),
     &          SVX201  (MSMAX,MLMAX,NSPECIES), 
     &          SVQ101  (MSMAX,MLMAX,NSPECIES),
     &          SVQ201  (MSMAX,MLMAX,NSPECIES),
     &          SVQ301  (MSMAX,MLMAX,NSPECIES),
     &          SVDP01  (MSMAX,MLMAX,NSPECIES) )
      

      ALLOCATE( VX1LNP(MSMAX), VX2LNP(MSMAX), VDPLNP(MSMAX),
     &          VQ1LNP(MSMAX), VQ2LNP(MSMAX), VQ3LNP(MSMAX))
C     ENDIF

C     THE FOLLOWING CALLS NEED TO BE EXECUTED EXACTLY IN ORDER
C     DEFINE LAMBDA-MESH AND EQUILIBRIUM CHI-MESH
      CALL KLAMBDA

C     COMPUTE AND STORE PARTICLE DRIFT FREQUENCIES
C     IT IS DESIRABLE TO PARALLELISE THIS SUBROUTINE FOR EACH SURFACE
      CALL KDRIFTFREQ

C     COMPUTE QUANTITIES ALPHAA1,2,3 ETC. FOR HOT ION SPECIES
C     IT IS DESIRABLE TO PARALLELISE THIS SUBROUTINE FOR EACH SURFACE
      CALL KALPHAA

C     COMPUTE DIAMAGNETIC FREQUECIES
C     SPECIFY COLLISIONALITY COEFFICIENT
      CALL KDIAMAG

C     COMPUTE THERMAL ION DRIFT ORBIT AT SURFACE JSOUT
      IF (IORBIT.GT.0) CALL KORBIT
      
3     CONTINUE

C     MODE FREQUENCY IN PLASMA FRAME
      IF (ABS(ATAU).LT.1.0E+10) THEN
         OMEGA = CI*ATAU
      ELSE
         OMEGA = CI*ALNORM
      ENDIF
      IF (DIMAG(OMEGA).LT.0.0) STOP 'KINETIC:OMEGA'


C     COMPUTE COEFFICIENTS FROM KINETIC PRESSURE TERMS
C     1: INTEGER RADIAL GRID
C     2: HALF-INTEGER RADIAL GRID

      CALL INIT_SURFACE_QUANTITY ()

      FREQK = 0.
      FREQKSURF = 0.
      FREQKSURF = 0.
C     CALL KBEAMLAM

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    MPI  PARALLEL COMPUTATION REGION         CC
CC           09/2012    Z.R.WANG               CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      IF (ISMPIRUN.EQ.1.OR.ISMPIRUN.EQ.3) THEN
        ALLOCATE (BUFFER_DATA(MSMAX,MSMAX,30))
        TAG(1)=1
        TAG(2)=2
        TAG(3)=3
        TAG(4)=4
        BUFFERSIZE = MSMAX*MSMAX*18
        
        IF (RANK.EQ.ROOT) THEN
C       MASTER PROCESS
C       BROADCAST THE NEW EIGENVALUE        
            CALL GET_EIGENVALUE(BUFFER_EIGENVALUE)
        
C       CREATE MPI TASK QUEUE         
            TOTTASK = 2*NR-1   
            ALLOCATE(TASKQUEUE(3,TOTTASK))
            TASKNUM = TOTTASK
            
            SLAVENUM = WORLDSIZE - 1
            IF (SLAVENUM .GT. TOTTASK) THEN
                PRINT *,'MORE SLAVE PROCESS: SLAVENUM=',SLAVENUM,
     &                  'TOTTASK=',TOTTASK
                SLAVENUM = TOTTASK
            ENDIF
            
            TASKQUEUE = 0
            COUNTER = 0
            DO JS=2,NR
                COUNTER = COUNTER + 1
                TASKQUEUE(1,COUNTER)=JS
                TASKQUEUE(2,COUNTER)=1  
            ENDDO
            DO JS=1,NR
                COUNTER = COUNTER + 1
                TASKQUEUE(1,COUNTER)=JS
                TASKQUEUE(2,COUNTER)=2 
            ENDDO
            
C           ASSIGN TASK TO SLAVE PROCESS            
            COUNTER = 1
            DO INDEX=1,SLAVENUM
                BUFFER_SEND(1)=TASKQUEUE(1,COUNTER)
                BUFFER_SEND(2)=TASKQUEUE(2,COUNTER)
                BUFFER_SEND(3)=ISWEEP          
                CALL MPI_SSEND( BUFFER_SEND,3, MPI_INTEGER, 
     $               INDEX,TAG(1), MPI_COMM_WORLD,STAT,IERP)
                CALL CHECK_MPI_ERROR()
                
                CALL MPI_SSEND( BUFFER_EIGENVALUE,1, MPI_DOUBLE_COMPLEX,
     $               INDEX,TAG(2), MPI_COMM_WORLD,STAT,IERP)
                CALL CHECK_MPI_ERROR()

                WRITE(*,*) 'KINETIC SEND TASK: JS,KGRID,PROCESS=',
     &                     BUFFER_SEND(1),BUFFER_SEND(2),INDEX

                COUNTER = COUNTER + 1
            ENDDO
            
            DO WHILE (TASKNUM .GT. 0)
            
                CALL MPI_RECV(BUFFER_REC,2,MPI_INTEGER,
     &           MPI_ANY_SOURCE,TAG(3),MPI_COMM_WORLD,STAT,IERP)
                CALL CHECK_MPI_ERROR()
                
                INDEX = STAT(MPI_SOURCE)
                
                CALL MPI_RECV(BUFFER_DATA,BUFFERSIZE,MPI_DOUBLE_COMPLEX,
     &          INDEX,TAG(4),MPI_COMM_WORLD,STAT,IERP)    
                CALL CHECK_MPI_ERROR()
C               SET MATRIX                
                
                MPIJS = BUFFER_REC(1)
                MPIKGRID = BUFFER_REC(2)
            IF(MPIKGRID.EQ.1) THEN
                VX1PARA(:,:,MPIJS) = BUFFER_DATA(:,:,1) 
                VX2PARA(:,:,MPIJS) = BUFFER_DATA(:,:,2) 
                VQ1PARA(:,:,MPIJS) = BUFFER_DATA(:,:,3)  
                VQ2PARA(:,:,MPIJS) = BUFFER_DATA(:,:,4)  
                VQ3PARA(:,:,MPIJS) = BUFFER_DATA(:,:,5)  
                VDPPARA(:,:,MPIJS) = BUFFER_DATA(:,:,6)  
                VX1PERP(:,:,MPIJS) = BUFFER_DATA(:,:,7) 
                VX2PERP(:,:,MPIJS) = BUFFER_DATA(:,:,8)  
                VQ1PERP(:,:,MPIJS) = BUFFER_DATA(:,:,9)  
                VQ2PERP(:,:,MPIJS) = BUFFER_DATA(:,:,10)   
                VQ3PERP(:,:,MPIJS) = BUFFER_DATA(:,:,11) 
                VDPPERP(:,:,MPIJS) = BUFFER_DATA(:,:,12) 

                IF (INCDPHI.GT.0) THEN
                VX1DPHI(:,:,MPIJS) = BUFFER_DATA(:,:,13) 
                VX2DPHI(:,:,MPIJS) = BUFFER_DATA(:,:,14) 
                VQ1DPHI(:,:,MPIJS) = BUFFER_DATA(:,:,15) 
                VQ2DPHI(:,:,MPIJS) = BUFFER_DATA(:,:,16) 
                VQ3DPHI(:,:,MPIJS) = BUFFER_DATA(:,:,17) 
                VDPDPHI(:,:,MPIJS) = BUFFER_DATA(:,:,18) 
                ENDIF

                IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &              .AND.INCKIN.GT.0) THEN
                VX1PARE(:,:,MPIJS) = BUFFER_DATA(:,:,19) 
                VX2PARE(:,:,MPIJS) = BUFFER_DATA(:,:,20) 
                VQ1PARE(:,:,MPIJS) = BUFFER_DATA(:,:,21)  
                VQ2PARE(:,:,MPIJS) = BUFFER_DATA(:,:,22)  
                VQ3PARE(:,:,MPIJS) = BUFFER_DATA(:,:,23)  
                VDPPARE(:,:,MPIJS) = BUFFER_DATA(:,:,24)  
                VX1PERE(:,:,MPIJS) = BUFFER_DATA(:,:,25) 
                VX2PERE(:,:,MPIJS) = BUFFER_DATA(:,:,26)  
                VQ1PERE(:,:,MPIJS) = BUFFER_DATA(:,:,27)  
                VQ2PERE(:,:,MPIJS) = BUFFER_DATA(:,:,28)   
                VQ3PERE(:,:,MPIJS) = BUFFER_DATA(:,:,29) 
                VDPPERE(:,:,MPIJS) = BUFFER_DATA(:,:,30) 
                ENDIF
            ELSE
                VX1PARAM(:,:,MPIJS) = BUFFER_DATA(:,:,1) 
                VX2PARAM(:,:,MPIJS) = BUFFER_DATA(:,:,2)
                VQ1PARAM(:,:,MPIJS) = BUFFER_DATA(:,:,3)
                VQ2PARAM(:,:,MPIJS) = BUFFER_DATA(:,:,4)
                VQ3PARAM(:,:,MPIJS) = BUFFER_DATA(:,:,5) 
                VDPPARAM(:,:,MPIJS) = BUFFER_DATA(:,:,6) 
                VX1PERPM(:,:,MPIJS) = BUFFER_DATA(:,:,7) 
                VX2PERPM(:,:,MPIJS) = BUFFER_DATA(:,:,8) 
                VQ1PERPM(:,:,MPIJS) = BUFFER_DATA(:,:,9) 
                VQ2PERPM(:,:,MPIJS) = BUFFER_DATA(:,:,10) 
                VQ3PERPM(:,:,MPIJS) = BUFFER_DATA(:,:,11) 
                VDPPERPM(:,:,MPIJS) = BUFFER_DATA(:,:,12) 

                IF (INCDPHI.GT.0) THEN
                VX1DPHIM(:,:,MPIJS) = BUFFER_DATA(:,:,13) 
                VX2DPHIM(:,:,MPIJS) = BUFFER_DATA(:,:,14) 
                VQ1DPHIM(:,:,MPIJS) = BUFFER_DATA(:,:,15) 
                VQ2DPHIM(:,:,MPIJS) = BUFFER_DATA(:,:,16) 
                VQ3DPHIM(:,:,MPIJS) = BUFFER_DATA(:,:,17) 
                VDPDPHIM(:,:,MPIJS) = BUFFER_DATA(:,:,18) 
                ENDIF

                IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &              .AND.INCKIN.GT.0) THEN
                VX1PAREM(:,:,MPIJS) = BUFFER_DATA(:,:,19) 
                VX2PAREM(:,:,MPIJS) = BUFFER_DATA(:,:,20) 
                VQ1PAREM(:,:,MPIJS) = BUFFER_DATA(:,:,21)  
                VQ2PAREM(:,:,MPIJS) = BUFFER_DATA(:,:,22)  
                VQ3PAREM(:,:,MPIJS) = BUFFER_DATA(:,:,23)  
                VDPPAREM(:,:,MPIJS) = BUFFER_DATA(:,:,24)  
                VX1PEREM(:,:,MPIJS) = BUFFER_DATA(:,:,25) 
                VX2PEREM(:,:,MPIJS) = BUFFER_DATA(:,:,26)  
                VQ1PEREM(:,:,MPIJS) = BUFFER_DATA(:,:,27)  
                VQ2PEREM(:,:,MPIJS) = BUFFER_DATA(:,:,28)   
                VQ3PEREM(:,:,MPIJS) = BUFFER_DATA(:,:,29) 
                VDPPEREM(:,:,MPIJS) = BUFFER_DATA(:,:,30) 
                ENDIF
            ENDIF                
                
                TASKNUM = TASKNUM - 1
                
                IF (COUNTER .LE. TOTTASK) THEN
                    BUFFER_SEND(1)=TASKQUEUE(1,COUNTER)
                    BUFFER_SEND(2)=TASKQUEUE(2,COUNTER)
                    BUFFER_SEND(3)=ISWEEP
                    CALL MPI_SSEND( BUFFER_SEND,3, MPI_INTEGER, 
     $              INDEX,TAG(1), MPI_COMM_WORLD,STAT,IERP)
                    CALL CHECK_MPI_ERROR()
                    CALL MPI_SSEND( BUFFER_EIGENVALUE,1,  
     $        MPI_DOUBLE_COMPLEX,INDEX,TAG(2), MPI_COMM_WORLD,STAT,IERP)
                    CALL CHECK_MPI_ERROR()

                    WRITE(*,*) 'KINETIC SEND TASK: JS,KGRID,PROCESS=',
     &                         BUFFER_SEND(1),BUFFER_SEND(2),INDEX

                    COUNTER = COUNTER + 1
                ENDIF
            ENDDO        
            IF (TASKNUM.EQ.0) THEN
               BUFFER_SEND(1)=-1
               DO INDEX=1,SLAVENUM
                  CALL MPI_SSEND(BUFFER_SEND,3,MPI_INTEGER,INDEX,
     &                           TAG(1),MPI_COMM_WORLD,STAT,IERP)
                  CALL CHECK_MPI_ERROR()
               ENDDO
            ENDIF
               
            DEALLOCATE(TASKQUEUE,BUFFER_DATA)        
        ELSE
C       SLAVE PROCESS

        
C       START KINETIC CALCULATION        
        DO WHILE (1.EQ.1)        
            CALL MPI_OPEN_FILE(RANK)
            WRITE(CHMPI,*) 'START RECEIVING TAG 1,RANK=',RANK
            CALL MPI_CLOSE_FILE(RANK)

            CALL MPI_RECV(BUFFER_REC,3,MPI_INTEGER,
     &                ROOT,TAG(1),MPI_COMM_WORLD,STAT,IERP)     
            CALL CHECK_MPI_ERROR()

            CALL MPI_OPEN_FILE(RANK)
            WRITE(CHMPI,*) 'FINISH RECEIVING TAG 1,RANK=',RANK
            CALL MPI_CLOSE_FILE(RANK)

            IF (BUFFER_REC(1).EQ.-1) THEN
               IF (ISMPIRUN.EQ.1.OR.ISMPIRUN.EQ.3)
     &            DEALLOCATE(BUFFER_DATA)
               IF (ISMPIRUN.EQ.1) THEN
                  CALL MPI_OPEN_FILE(RANK)
                  WRITE(CHMPI,*) 'STOP RANK=',RANK
                  CALL MPI_CLOSE_FILE(RANK)
                  CALL MPI_FINALIZE(IERP)
                  STOP
               ENDIF 
            ENDIF 

            CALL MPI_RECV(BUFFER_EIGENVALUE,1,MPI_DOUBLE_COMPLEX,
     &                ROOT,TAG(2),MPI_COMM_WORLD,STAT,IERP)
            CALL CHECK_MPI_ERROR()
            
            CALL SET_EIGENVALUE(BUFFER_EIGENVALUE)
            
            MPIJS = BUFFER_REC(1)
            MPIKGRID = BUFFER_REC(2)
            ISWEEP = BUFFER_REC(3)
            
            CALL MPI_OPEN_FILE(RANK)
            WRITE(CHMPI,*) 'RECEIVE TASK: RANK=',RANK,'JS=',
     &               MPIJS,'KGRID=',MPIKGRID 
            CALL MPI_CLOSE_FILE(RANK)

            CALL INIT_SURFACE_QUANTITY ()
            
            BUFFER_SEND(1) = MPIJS
            BUFFER_SEND(2) = MPIKGRID
            
            CALL KJPCOEFF( MPIJS, 1, MPIKGRID )
            
            IF(BUFFER_SEND(2).EQ.1) THEN
                BUFFER_DATA(:,:,1)  = VX1PARA(:,:,1)
                BUFFER_DATA(:,:,2)  = VX2PARA(:,:,1)
                BUFFER_DATA(:,:,3)  = VQ1PARA(:,:,1)
                BUFFER_DATA(:,:,4)  = VQ2PARA(:,:,1)
                BUFFER_DATA(:,:,5)  = VQ3PARA(:,:,1)
                BUFFER_DATA(:,:,6)  = VDPPARA(:,:,1)
                BUFFER_DATA(:,:,7)  = VX1PERP(:,:,1)
                BUFFER_DATA(:,:,8)  = VX2PERP(:,:,1)
                BUFFER_DATA(:,:,9)  = VQ1PERP(:,:,1)
                BUFFER_DATA(:,:,10) = VQ2PERP(:,:,1)
                BUFFER_DATA(:,:,11) = VQ3PERP(:,:,1)
                BUFFER_DATA(:,:,12) = VDPPERP(:,:,1)

                IF (INCDPHI.GT.0) THEN
                BUFFER_DATA(:,:,13) = VX1DPHI(:,:,1)
                BUFFER_DATA(:,:,14) = VX2DPHI(:,:,1)
                BUFFER_DATA(:,:,15) = VQ1DPHI(:,:,1)
                BUFFER_DATA(:,:,16) = VQ2DPHI(:,:,1)
                BUFFER_DATA(:,:,17) = VQ3DPHI(:,:,1)
                BUFFER_DATA(:,:,18) = VDPDPHI(:,:,1)
                ENDIF

                IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &              .AND.INCKIN.GT.0) THEN
                BUFFER_DATA(:,:,19) = VX1PARE(:,:,1)
                BUFFER_DATA(:,:,20) = VX2PARE(:,:,1)
                BUFFER_DATA(:,:,21) = VQ1PARE(:,:,1)
                BUFFER_DATA(:,:,22) = VQ2PARE(:,:,1)
                BUFFER_DATA(:,:,23) = VQ3PARE(:,:,1)
                BUFFER_DATA(:,:,24) = VDPPARE(:,:,1)
                BUFFER_DATA(:,:,25) = VX1PERE(:,:,1)
                BUFFER_DATA(:,:,26) = VX2PERE(:,:,1)
                BUFFER_DATA(:,:,27) = VQ1PERE(:,:,1)
                BUFFER_DATA(:,:,28) = VQ2PERE(:,:,1)
                BUFFER_DATA(:,:,29) = VQ3PERE(:,:,1)
                BUFFER_DATA(:,:,30) = VDPPERE(:,:,1)
                ENDIF
            ELSE
                BUFFER_DATA(:,:,1)  = VX1PARAM(:,:,1)
                BUFFER_DATA(:,:,2)  = VX2PARAM(:,:,1)
                BUFFER_DATA(:,:,3)  = VQ1PARAM(:,:,1)
                BUFFER_DATA(:,:,4)  = VQ2PARAM(:,:,1)
                BUFFER_DATA(:,:,5)  = VQ3PARAM(:,:,1)
                BUFFER_DATA(:,:,6)  = VDPPARAM(:,:,1)
                BUFFER_DATA(:,:,7)  = VX1PERPM(:,:,1)
                BUFFER_DATA(:,:,8)  = VX2PERPM(:,:,1)
                BUFFER_DATA(:,:,9)  = VQ1PERPM(:,:,1)
                BUFFER_DATA(:,:,10) = VQ2PERPM(:,:,1)
                BUFFER_DATA(:,:,11) = VQ3PERPM(:,:,1)
                BUFFER_DATA(:,:,12) = VDPPERPM(:,:,1)

                IF (INCDPHI.GT.0) THEN
                BUFFER_DATA(:,:,13) = VX1DPHIM(:,:,1)
                BUFFER_DATA(:,:,14) = VX2DPHIM(:,:,1)
                BUFFER_DATA(:,:,15) = VQ1DPHIM(:,:,1)
                BUFFER_DATA(:,:,16) = VQ2DPHIM(:,:,1)
                BUFFER_DATA(:,:,17) = VQ3DPHIM(:,:,1)
                BUFFER_DATA(:,:,18) = VDPDPHIM(:,:,1)
                ENDIF

                IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &              .AND.INCKIN.GT.0) THEN
                BUFFER_DATA(:,:,19) = VX1PAREM(:,:,1)
                BUFFER_DATA(:,:,20) = VX2PAREM(:,:,1)
                BUFFER_DATA(:,:,21) = VQ1PAREM(:,:,1)
                BUFFER_DATA(:,:,22) = VQ2PAREM(:,:,1)
                BUFFER_DATA(:,:,23) = VQ3PAREM(:,:,1)
                BUFFER_DATA(:,:,24) = VDPPAREM(:,:,1)
                BUFFER_DATA(:,:,25) = VX1PEREM(:,:,1)
                BUFFER_DATA(:,:,26) = VX2PEREM(:,:,1)
                BUFFER_DATA(:,:,27) = VQ1PEREM(:,:,1)
                BUFFER_DATA(:,:,28) = VQ2PEREM(:,:,1)
                BUFFER_DATA(:,:,29) = VQ3PEREM(:,:,1)
                BUFFER_DATA(:,:,30) = VDPPEREM(:,:,1)
                ENDIF
            ENDIF
            
            call MPI_SSEND( BUFFER_SEND,2,MPI_INTEGER,ROOT,TAG(3), 
     $                      MPI_COMM_WORLD,STAT,IERP)
            CALL CHECK_MPI_ERROR()
            
            call MPI_SSEND( BUFFER_DATA,BUFFERSIZE, MPI_DOUBLE_COMPLEX, 
     $                 ROOT,TAG(4),MPI_COMM_WORLD,STAT,IERP)
            CALL CHECK_MPI_ERROR()
            
        ENDDO    
            
        ENDIF 
        
         
      ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    PARALLEL COMPUTATION REGION   OPENMP     CC
CC              2010  Z.R.WANG                 CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      CALL OMP_SET_NUM_THREADS(IOMPNUM)
      
      JS=2
           
CC$OMP PARALLEL DEFAULT ( SHARED )
C$OMP PARALLEL 
C$OMP&COPYIN( MLMAX,NCHIT,NCHI1,NCHI2,JS0 )
C$OMP&COPYIN( KJPKEY,KPBKEY )

      CALL ALLOCATEFORPARALLEL()

      DO WHILE (JS .LE. NR)         
C$OMP CRITICAL(PARALLEL_LOOP1)         
         PRIVATEJS = JS
         JS = JS + 1      
C        WRITE (*,*) 'THREAD ID=',OMP_GET_THREAD_NUM(),
C    $   'GRID 1 SURFACE ',PRIVATEJS
c$OMP END CRITICAL(PARALLEL_LOOP1)

         CALL KJPCOEFF( PRIVATEJS,PRIVATEJS, 1 )

      ENDDO
     
      CALL DEALLOCATEFORPARALLEL()
      
C$OMP BARRIER      
C$OMP END PARALLEL

      JS = 1
      
CC$OMP PARALLEL DEFAULT ( SHARED )
C$OMP PARALLEL 
C$OMP&COPYIN( MLMAX,NCHIT,NCHI1,NCHI2,JS0 )
C$OMP&COPYIN( KJPKEY,KPBKEY )

      CALL ALLOCATEFORPARALLEL()

      DO WHILE (JS .LE. NR)
C$OMP CRITICAL(PARALLEL_LOOP2)           
         PRIVATEJS = JS
         JS = JS + 1        
C        WRITE (*,*) 'THREAD ID=',OMP_GET_THREAD_NUM(),
C    $   'GRID 2 SURFACE ',PRIVATEJS
c$OMP END CRITICAL(PARALLEL_LOOP2)
      
         CALL KJPCOEFF(PRIVATEJS,PRIVATEJS,2)

      END DO

      CALL DEALLOCATEFORPARALLEL()

C$OMP BARRIER 
C$OMP END PARALLEL

      END IF
      
      IF (ISMPIRUN.EQ.0) THEN
      
C     OUTPUT ALL DRIFT FREQUENCIES
      OPEN(31,FILE='FREQUENCIES.OUT',FORM='FORMATTED')
      REWIND(31)
      DO JS=2,NRP1-1
      WRITE(31,1001) CS(JS), OMEGAE0(JS,1),  OMEGASN(JS,1,1),
     &                       OMEGASN(JS,1,2), OMEGAST(JS,1,1),
     &                       OMEGAST(JS,1,2), FREQK(JS,1),
     &                       FREQK(JS,2),    FREQK(JS,3),    
     &                       FREQK(JS,4),    FREQK(JS,5),    
     &                       FREQK(JS,6),    FREQK(JS,7),    
     &                       FREQK(JS,8),    FREQK(JS,9),    
     &                       FREQK(JS,10),   FREQK(JS,11),    
     &                       FREQK(JS,12),   FREQK(JS,13)
      ENDDO
 1001 FORMAT(19(E16.8,1X))
      CLOSE(31)
      
 1002 FORMAT ( E16.8,1X,$ )

      SURFFILENUM = assignFreeFileUnit ()
      OPEN(SURFFILENUM, FILE='POSITION_KAI.OUT',FORM='FORMATTED')
      REWIND(SURFFILENUM)
      DO JS = 2,NRP1 - 1
        DO J = 1, NCHI
            WRITE (SURFFILENUM,1002) POSITIONKAI(JS,J)
        ENDDO
            WRITE (SURFFILENUM,*)
      ENDDO
      CLOSE (SURFFILENUM)


      OPEN(SURFFILENUM, FILE='BOUNCESURF_FREQ.OUT',FORM='FORMATTED')
      REWIND(SURFFILENUM)
      DO JS = 2,NRP1 - 1
        DO J = 1, NCHI
            WRITE (SURFFILENUM,1002) FREQKSURF(JS,2,J)
        ENDDO
            WRITE (SURFFILENUM,*)
      ENDDO
      CLOSE (SURFFILENUM)

      OPEN(SURFFILENUM, FILE='PRECESSION_FREQ.OUT',FORM='FORMATTED')
      REWIND(SURFFILENUM)
      DO JS = 2,NRP1 - 1
        DO J = 1, NCHI
            WRITE (SURFFILENUM,1002) FREQKSURF(JS,3,J)
        ENDDO
            WRITE (SURFFILENUM,*)
      ENDDO
      CLOSE (SURFFILENUM)

      OPEN(SURFFILENUM, FILE='PASSINGPART_FRAC.OUT',FORM='FORMATTED')
      REWIND(SURFFILENUM)
      DO JS = 2,NRP1 - 1
        DO J = 1, NCHI
            WRITE (SURFFILENUM,1002) FREQKSURF(JS,6,J)
        ENDDO
            WRITE (SURFFILENUM,*)
      ENDDO
      CLOSE (SURFFILENUM)

      OPEN(SURFFILENUM, FILE='TRAPPEDPART_FRAC.OUT',FORM='FORMATTED')
      REWIND(SURFFILENUM)
      DO JS = 2,NRP1 - 1
        DO J = 1, NCHI
            WRITE (SURFFILENUM,1002) FREQKSURF(JS,7,J)
        ENDDO
            WRITE (SURFFILENUM,*)
      ENDDO
      CLOSE (SURFFILENUM)

      OPEN(SURFFILENUM, FILE='TRANSITPART_FREQ.OUT',FORM='FORMATTED')
      REWIND(SURFFILENUM)
      DO JS = 2,NRP1 - 1
        DO J = 1, NCHI
            WRITE (SURFFILENUM,1002) FREQKSURF(JS,1,J)
        ENDDO
            WRITE (SURFFILENUM,*)
      ENDDO
      CLOSE (SURFFILENUM)

      ENDIF


C     SET DAMPING COEFFICIENTS AT THE MAGNETIC AXIS
      DO MROW=1,MSMAX
         DO MSA=1,MSMAX
            VX1PARA(MROW,MSA,1) = VX1PARAM(MROW,MSA,1)
            VX1PERP(MROW,MSA,1) = VX1PERPM(MROW,MSA,1)
            VX2PARA(MROW,MSA,1) = VX2PARAM(MROW,MSA,1)
            VX2PERP(MROW,MSA,1) = VX2PERPM(MROW,MSA,1)
            VQ1PARA(MROW,MSA,1) = VQ1PARAM(MROW,MSA,1)
            VQ1PERP(MROW,MSA,1) = VQ1PERPM(MROW,MSA,1)
            VQ2PARA(MROW,MSA,1) = VQ2PARAM(MROW,MSA,1)
            VQ2PERP(MROW,MSA,1) = VQ2PERPM(MROW,MSA,1)
            VQ3PARA(MROW,MSA,1) = VQ3PARAM(MROW,MSA,1)
            VQ3PERP(MROW,MSA,1) = VQ3PERPM(MROW,MSA,1)
            VDPPARA(MROW,MSA,1) = VDPPARAM(MROW,MSA,1)
            VDPPERP(MROW,MSA,1) = VDPPERPM(MROW,MSA,1)

            IF (INCDPHI.GT.0) THEN
            VX1DPHI(MROW,MSA,1) = VX1DPHIM(MROW,MSA,1)
            VX2DPHI(MROW,MSA,1) = VX2DPHIM(MROW,MSA,1)
            VQ1DPHI(MROW,MSA,1) = VQ1DPHIM(MROW,MSA,1)
            VQ2DPHI(MROW,MSA,1) = VQ2DPHIM(MROW,MSA,1)
            VQ3DPHI(MROW,MSA,1) = VQ3DPHIM(MROW,MSA,1)
            VDPDPHI(MROW,MSA,1) = VDPDPHIM(MROW,MSA,1)
            ENDIF

            IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &          .AND.INCKIN.GT.0) THEN
            VX1PARE(MROW,MSA,1) = VX1PAREM(MROW,MSA,1)
            VX1PERE(MROW,MSA,1) = VX1PEREM(MROW,MSA,1)
            VX2PARE(MROW,MSA,1) = VX2PAREM(MROW,MSA,1)
            VX2PERE(MROW,MSA,1) = VX2PEREM(MROW,MSA,1)
            VQ1PARE(MROW,MSA,1) = VQ1PAREM(MROW,MSA,1)
            VQ1PERE(MROW,MSA,1) = VQ1PEREM(MROW,MSA,1)
            VQ2PARE(MROW,MSA,1) = VQ2PAREM(MROW,MSA,1)
            VQ2PERE(MROW,MSA,1) = VQ2PEREM(MROW,MSA,1)
            VQ3PARE(MROW,MSA,1) = VQ3PAREM(MROW,MSA,1)
            VQ3PERE(MROW,MSA,1) = VQ3PEREM(MROW,MSA,1)
            VDPPARE(MROW,MSA,1) = VDPPAREM(MROW,MSA,1)
            VDPPERE(MROW,MSA,1) = VDPPEREM(MROW,MSA,1)
            ENDIF
         ENDDO
      ENDDO

C     SET DAMPING COEFFICIENTS AT PLASMA SURFACE    
      I = NRP1
      J = NR
      DO MROW=1,MSMAX
         DO MSA=1,MSMAX
            VX1PARA(MROW,MSA,I) = VX1PARAM(MROW,MSA,J)
            VX1PERP(MROW,MSA,I) = VX1PERPM(MROW,MSA,J)
            VX2PARA(MROW,MSA,I) = VX2PARAM(MROW,MSA,J)
            VX2PERP(MROW,MSA,I) = VX2PERPM(MROW,MSA,J)
            VQ1PARA(MROW,MSA,I) = VQ1PARAM(MROW,MSA,J)
            VQ1PERP(MROW,MSA,I) = VQ1PERPM(MROW,MSA,J)
            VQ2PARA(MROW,MSA,I) = VQ2PARAM(MROW,MSA,J)
            VQ2PERP(MROW,MSA,I) = VQ2PERPM(MROW,MSA,J)
            VQ3PARA(MROW,MSA,I) = VQ3PARAM(MROW,MSA,J)
            VQ3PERP(MROW,MSA,I) = VQ3PERPM(MROW,MSA,J)
            VDPPARA(MROW,MSA,I) = VDPPARAM(MROW,MSA,J)
            VDPPERP(MROW,MSA,I) = VDPPERPM(MROW,MSA,J)

            IF (INCDPHI.GT.0) THEN
            VX1DPHI(MROW,MSA,I) = VX1DPHIM(MROW,MSA,J)
            VX2DPHI(MROW,MSA,I) = VX2DPHIM(MROW,MSA,J)
            VQ1DPHI(MROW,MSA,I) = VQ1DPHIM(MROW,MSA,J)
            VQ2DPHI(MROW,MSA,I) = VQ2DPHIM(MROW,MSA,J)
            VQ3DPHI(MROW,MSA,I) = VQ3DPHIM(MROW,MSA,J)
            VDPDPHI(MROW,MSA,I) = VDPDPHIM(MROW,MSA,J)
            ENDIF

            IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &          .AND.INCKIN.GT.0) THEN
            VX1PARE(MROW,MSA,I) = VX1PAREM(MROW,MSA,J)
            VX1PERE(MROW,MSA,I) = VX1PEREM(MROW,MSA,J)
            VX2PARE(MROW,MSA,I) = VX2PAREM(MROW,MSA,J)
            VX2PERE(MROW,MSA,I) = VX2PEREM(MROW,MSA,J)
            VQ1PARE(MROW,MSA,I) = VQ1PAREM(MROW,MSA,J)
            VQ1PERE(MROW,MSA,I) = VQ1PEREM(MROW,MSA,J)
            VQ2PARE(MROW,MSA,I) = VQ2PAREM(MROW,MSA,J)
            VQ2PERE(MROW,MSA,I) = VQ2PEREM(MROW,MSA,J)
            VQ3PARE(MROW,MSA,I) = VQ3PAREM(MROW,MSA,J)
            VQ3PERE(MROW,MSA,I) = VQ3PEREM(MROW,MSA,J)
            VDPPARE(MROW,MSA,I) = VDPPAREM(MROW,MSA,J)
            VDPPERE(MROW,MSA,I) = VDPPEREM(MROW,MSA,J)
            ENDIF
         ENDDO
      ENDDO

C     SET OVERALL FRACTION ALPHAD FOR KINETIC CONTRIBUTIONS
      VX1PARA  = VX1PARA*ALPHAD
      VX1PERP  = VX1PERP*ALPHAD
      VX2PARA  = VX2PARA*ALPHAD
      VX2PERP  = VX2PERP*ALPHAD
      VQ1PARA  = VQ1PARA*ALPHAD
      VQ1PERP  = VQ1PERP*ALPHAD
      VQ2PARA  = VQ2PARA*ALPHAD
      VQ2PERP  = VQ2PERP*ALPHAD
      VQ3PARA  = VQ3PARA*ALPHAD
      VQ3PERP  = VQ3PERP*ALPHAD
      VDPPARA  = VDPPARA*ALPHAD
      VDPPERP  = VDPPERP*ALPHAD

      VX1PARAM = VX1PARAM*ALPHAD
      VX1PERPM = VX1PERPM*ALPHAD
      VX2PARAM = VX2PARAM*ALPHAD
      VX2PERPM = VX2PERPM*ALPHAD
      VQ1PARAM = VQ1PARAM*ALPHAD
      VQ1PERPM = VQ1PERPM*ALPHAD
      VQ2PARAM = VQ2PARAM*ALPHAD
      VQ2PERPM = VQ2PERPM*ALPHAD
      VQ3PARAM = VQ3PARAM*ALPHAD
      VQ3PERPM = VQ3PERPM*ALPHAD
      VDPPARAM = VDPPARAM*ALPHAD
      VDPPERPM = VDPPERPM*ALPHAD
      
      IF (INCDPHI.GT.0) THEN
      VX1DPHI  = VX1DPHI*ALPHAD
      VX2DPHI  = VX2DPHI*ALPHAD
      VQ1DPHI  = VQ1DPHI*ALPHAD
      VQ2DPHI  = VQ2DPHI*ALPHAD
      VQ3DPHI  = VQ3DPHI*ALPHAD
      VDPDPHI  = VDPDPHI*ALPHAD

      VX1DPHIM = VX1DPHIM*ALPHAD
      VX2DPHIM = VX2DPHIM*ALPHAD
      VQ1DPHIM = VQ1DPHIM*ALPHAD
      VQ2DPHIM = VQ2DPHIM*ALPHAD
      VQ3DPHIM = VQ3DPHIM*ALPHAD
      VDPDPHIM = VDPDPHIM*ALPHAD
      ENDIF

      IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &    .AND.INCKIN.GT.0) THEN
      VX1PARE  = VX1PARE*ALPHAD
      VX1PERE  = VX1PERE*ALPHAD
      VX2PARE  = VX2PARE*ALPHAD
      VX2PERE  = VX2PERE*ALPHAD
      VQ1PARE  = VQ1PARE*ALPHAD
      VQ1PERE  = VQ1PERE*ALPHAD
      VQ2PARE  = VQ2PARE*ALPHAD
      VQ2PERE  = VQ2PERE*ALPHAD
      VQ3PARE  = VQ3PARE*ALPHAD
      VQ3PERE  = VQ3PERE*ALPHAD
      VDPPARE  = VDPPARE*ALPHAD
      VDPPERE  = VDPPERE*ALPHAD

      VX1PAREM = VX1PAREM*ALPHAD
      VX1PEREM = VX1PEREM*ALPHAD
      VX2PAREM = VX2PAREM*ALPHAD
      VX2PEREM = VX2PEREM*ALPHAD
      VQ1PAREM = VQ1PAREM*ALPHAD
      VQ1PEREM = VQ1PEREM*ALPHAD
      VQ2PAREM = VQ2PAREM*ALPHAD
      VQ2PEREM = VQ2PEREM*ALPHAD
      VQ3PAREM = VQ3PAREM*ALPHAD
      VQ3PEREM = VQ3PEREM*ALPHAD
      VDPPAREM = VDPPAREM*ALPHAD
      VDPPEREM = VDPPEREM*ALPHAD
      ENDIF

C     FILL IN MATRICES FOR EQUATIONS FOR PARALLEL AND 
C     PERPENDICULAR KINETIC PRESSURE (RHS ONLY, WHICH
C     DOES NOT REQUIRE CONVOLUTION)
      IF (IPERTURB.EQ.0 .AND. (V2XKEY.EQ.0 .OR. V2XKEY.EQ.2) ) THEN

      DO I=1,NR
         INCLUDE 'tophat.inc'

         ZV2M = (CS(I  )/CSM(I))**IEXV2
         ZV2P = (CS(I+1)/CSM(I))**IEXV2
         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         DO MROW=1,MSMAX
            LYROW = (MROW-1)*NYCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               FSUBM(KYPPARA+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PARA(MROW,MSA,I),VX1PARAM(MROW,MSA,I))

               FSUBM(KYPPERP+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PERP(MROW,MSA,I),VX1PERPM(MROW,MSA,I))

               GSUBM(KYPPARA+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PARA(MROW,MSA,I+1),VX1PARAM(MROW,MSA,I))

               GSUBM(KYPPERP+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PERP(MROW,MSA,I+1),VX1PERPM(MROW,MSA,I))

               DSUBM(KYPPARA+LYROW, KYV2+LYCOL,I)=
     &           GG(VX2PARAM(MROW,MSA,I), VX2PARA(MROW,MSA,I)*ZV2M,
     &           VX2PARA(MROW,MSA,I+1)*ZV2P)

               DSUBM(KYPPERP+LYROW, KYV2+LYCOL,I)=
     &           GG(VX2PERPM(MROW,MSA,I), VX2PERP(MROW,MSA,I)*ZV2M,
     &           VX2PERP(MROW,MSA,I+1)*ZV2P)

               FSUBM(KYPPARA+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I)*VQ1PARA(MROW,MSA,I), 
     &              SHIFTM(I)*VQ1PARAM(MROW,MSA,I))    

               FSUBM(KYPPERP+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I)*VQ1PERP(MROW,MSA,I), 
     &              SHIFTM(I)*VQ1PERPM(MROW,MSA,I))    

               GSUBM(KYPPARA+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VQ1PARA(MROW,MSA,I+1), 
     &              SHIFTM(I)*VQ1PARAM(MROW,MSA,I))    

               GSUBM(KYPPERP+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VQ1PERP(MROW,MSA,I+1), 
     &              SHIFTM(I)*VQ1PERPM(MROW,MSA,I))    

               DSUBM(KYPPARA+LYROW, KYB2+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ2PARAM(MROW,MSA,I), 
     &              SHIFTC(I)*VQ2PARA(MROW,MSA,I),
     &              SHIFTC(I+1)*VQ2PARA(MROW,MSA,I+1))

               DSUBM(KYPPERP+LYROW, KYB2+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ2PERPM(MROW,MSA,I), 
     &              SHIFTC(I)*VQ2PERP(MROW,MSA,I),
     &              SHIFTC(I+1)*VQ2PERP(MROW,MSA,I+1))

               DSUBM(KYPPARA+LYROW, KYB3+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ3PARAM(MROW,MSA,I), 
     &              SHIFTC(I)*VQ3PARA(MROW,MSA,I)*ZB3M,
     &              SHIFTC(I+1)*VQ3PARA(MROW,MSA,I+1)*ZB3P)

               DSUBM(KYPPERP+LYROW, KYB3+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ3PERPM(MROW,MSA,I), 
     &               SHIFTC(I)*VQ3PERP(MROW,MSA,I)*ZB3M,
     &               SHIFTC(I+1)*VQ3PERP(MROW,MSA,I+1)*ZB3P)

               IF (KXDPHI.GT.0) THEN
               FSUBM(KYPPARA+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I)*VDPPARA(MROW,MSA,I), 
     &              SHIFTM(I)*VDPPARAM(MROW,MSA,I))    

               FSUBM(KYPPERP+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I)*VDPPERP(MROW,MSA,I), 
     &              SHIFTM(I)*VDPPERPM(MROW,MSA,I))    

               GSUBM(KYPPARA+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VDPPARA(MROW,MSA,I+1), 
     &              SHIFTM(I)*VDPPARAM(MROW,MSA,I))    

               GSUBM(KYPPERP+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VDPPERP(MROW,MSA,I+1), 
     &              SHIFTM(I)*VDPPERPM(MROW,MSA,I))    
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &    .AND.INCKIN.GT.0) THEN
      DO I=1,NR
         INCLUDE 'tophat.inc'

         ZV2M = (CS(I  )/CSM(I))**IEXV2
         ZV2P = (CS(I+1)/CSM(I))**IEXV2
         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         DO MROW=1,MSMAX
            LYROW = (MROW-1)*NYCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               FSUBM(KYPE+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PARE(MROW,MSA,I),VX1PAREM(MROW,MSA,I))

               FSUBM(KYPP+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PERE(MROW,MSA,I),VX1PEREM(MROW,MSA,I))

               GSUBM(KYPE+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PARE(MROW,MSA,I+1),VX1PAREM(MROW,MSA,I))

               GSUBM(KYPP+LYROW, KXV1+LXCOL,I)=
     &           GF(VX1PERE(MROW,MSA,I+1),VX1PEREM(MROW,MSA,I))

               DSUBM(KYPE+LYROW, KYV2+LYCOL,I)=
     &           GG(VX2PAREM(MROW,MSA,I), VX2PARE(MROW,MSA,I)*ZV2M,
     &           VX2PARE(MROW,MSA,I+1)*ZV2P)

               DSUBM(KYPP+LYROW, KYV2+LYCOL,I)=
     &           GG(VX2PEREM(MROW,MSA,I), VX2PERE(MROW,MSA,I)*ZV2M,
     &           VX2PERE(MROW,MSA,I+1)*ZV2P)

               FSUBM(KYPE+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I)*VQ1PARE(MROW,MSA,I), 
     &              SHIFTM(I)*VQ1PAREM(MROW,MSA,I))    

               FSUBM(KYPP+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I)*VQ1PERE(MROW,MSA,I), 
     &              SHIFTM(I)*VQ1PEREM(MROW,MSA,I))    

               GSUBM(KYPE+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VQ1PARE(MROW,MSA,I+1), 
     &              SHIFTM(I)*VQ1PAREM(MROW,MSA,I))    

               GSUBM(KYPP+LYROW, KXB1+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VQ1PERE(MROW,MSA,I+1), 
     &              SHIFTM(I)*VQ1PEREM(MROW,MSA,I))    

               DSUBM(KYPE+LYROW, KYB2+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ2PAREM(MROW,MSA,I), 
     &              SHIFTC(I)*VQ2PARE(MROW,MSA,I),
     &              SHIFTC(I+1)*VQ2PARE(MROW,MSA,I+1))

               DSUBM(KYPP+LYROW, KYB2+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ2PEREM(MROW,MSA,I), 
     &              SHIFTC(I)*VQ2PERE(MROW,MSA,I),
     &              SHIFTC(I+1)*VQ2PERE(MROW,MSA,I+1))

               DSUBM(KYPE+LYROW, KYB3+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ3PAREM(MROW,MSA,I), 
     &              SHIFTC(I)*VQ3PARE(MROW,MSA,I)*ZB3M,
     &              SHIFTC(I+1)*VQ3PARE(MROW,MSA,I+1)*ZB3P)

               DSUBM(KYPP+LYROW, KYB3+LYCOL,I)=
     &          -GG(SHIFTM(I)*VQ3PEREM(MROW,MSA,I), 
     &               SHIFTC(I)*VQ3PERE(MROW,MSA,I)*ZB3M,
     &               SHIFTC(I+1)*VQ3PERE(MROW,MSA,I+1)*ZB3P)

               IF (KXDPHI.GT.0) THEN
               FSUBM(KYPE+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I)*VDPPARE(MROW,MSA,I), 
     &              SHIFTM(I)*VDPPAREM(MROW,MSA,I))    

               FSUBM(KYPP+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I)*VDPPERE(MROW,MSA,I), 
     &              SHIFTM(I)*VDPPEREM(MROW,MSA,I))    

               GSUBM(KYPE+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VDPPARE(MROW,MSA,I+1), 
     &              SHIFTM(I)*VDPPAREM(MROW,MSA,I))    

               GSUBM(KYPP+LYROW, KXDPHI+LXCOL,I)=
     &          -GF(SHIFTC(I+1)*VDPPERE(MROW,MSA,I+1), 
     &              SHIFTM(I)*VDPPEREM(MROW,MSA,I))    
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      ELSE

      DO I=1,NR
         INCLUDE 'tophat.inc'

         ZV2M = (CS(I  )/CSM(I))**IEXV2
         ZV2P = (CS(I+1)/CSM(I))**IEXV2
         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         DO MROW=1,MSMAX
            LYROW = (MROW-1)*NYCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               FSUBM(KYPPARA+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PARA(MROW,MSA,I),VX1PARAM(MROW,MSA,I))

               FSUBM(KYPPERP+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PERP(MROW,MSA,I),VX1PERPM(MROW,MSA,I))

               GSUBM(KYPPARA+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PARA(MROW,MSA,I+1),VX1PARAM(MROW,MSA,I))

               GSUBM(KYPPERP+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PERP(MROW,MSA,I+1),VX1PERPM(MROW,MSA,I))

               DSUBM(KYPPARA+LYROW, KYX2+LYCOL,I)=
     &           GG(VX2PARAM(MROW,MSA,I), VX2PARA(MROW,MSA,I)*ZV2M,
     &           VX2PARA(MROW,MSA,I+1)*ZV2P)

               DSUBM(KYPPERP+LYROW, KYX2+LYCOL,I)=
     &           GG(VX2PERPM(MROW,MSA,I), VX2PERP(MROW,MSA,I)*ZV2M,
     &           VX2PERP(MROW,MSA,I+1)*ZV2P)

               FSUBM(KYPPARA+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PARA(MROW,MSA,I), VQ1PARAM(MROW,MSA,I))    

               FSUBM(KYPPERP+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PERP(MROW,MSA,I), VQ1PERPM(MROW,MSA,I))    

               GSUBM(KYPPARA+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PARA(MROW,MSA,I+1), VQ1PARAM(MROW,MSA,I))    

               GSUBM(KYPPERP+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PERP(MROW,MSA,I+1), VQ1PERPM(MROW,MSA,I))    

               DSUBM(KYPPARA+LYROW, KYB2+LYCOL,I)=
     &           GG(VQ2PARAM(MROW,MSA,I), VQ2PARA(MROW,MSA,I),
     &           VQ2PARA(MROW,MSA,I+1))

               DSUBM(KYPPERP+LYROW, KYB2+LYCOL,I)=
     &           GG(VQ2PERPM(MROW,MSA,I), VQ2PERP(MROW,MSA,I),
     &           VQ2PERP(MROW,MSA,I+1))

               DSUBM(KYPPARA+LYROW, KYB3+LYCOL,I)=
     &           GG(VQ3PARAM(MROW,MSA,I), VQ3PARA(MROW,MSA,I)*ZB3M,
     &           VQ3PARA(MROW,MSA,I+1)*ZB3P)

               DSUBM(KYPPERP+LYROW, KYB3+LYCOL,I)=
     &           GG(VQ3PERPM(MROW,MSA,I), VQ3PERP(MROW,MSA,I)*ZB3M,
     &           VQ3PERP(MROW,MSA,I+1)*ZB3P)

               IF (KXDPHI.GT.0) THEN
               FSUBM(KYPPARA+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPARA(MROW,MSA,I), VDPPARAM(MROW,MSA,I))    

               FSUBM(KYPPERP+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPERP(MROW,MSA,I), VDPPERPM(MROW,MSA,I))    

               GSUBM(KYPPARA+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPARA(MROW,MSA,I+1), VDPPARAM(MROW,MSA,I))    

               GSUBM(KYPPERP+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPERP(MROW,MSA,I+1), VDPPERPM(MROW,MSA,I))    
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &    .AND.INCKIN.GT.0) THEN
      DO I=1,NR
         INCLUDE 'tophat.inc'

         ZV2M = (CS(I  )/CSM(I))**IEXV2
         ZV2P = (CS(I+1)/CSM(I))**IEXV2
         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         DO MROW=1,MSMAX
            LYROW = (MROW-1)*NYCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               FSUBM(KYPE+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PARE(MROW,MSA,I),VX1PAREM(MROW,MSA,I))

               FSUBM(KYPP+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PERE(MROW,MSA,I),VX1PEREM(MROW,MSA,I))

               GSUBM(KYPE+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PARE(MROW,MSA,I+1),VX1PAREM(MROW,MSA,I))

               GSUBM(KYPP+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PERE(MROW,MSA,I+1),VX1PEREM(MROW,MSA,I))

               DSUBM(KYPE+LYROW, KYX2+LYCOL,I)=
     &           GG(VX2PAREM(MROW,MSA,I), VX2PARE(MROW,MSA,I)*ZV2M,
     &           VX2PARE(MROW,MSA,I+1)*ZV2P)

               DSUBM(KYPP+LYROW, KYX2+LYCOL,I)=
     &           GG(VX2PEREM(MROW,MSA,I), VX2PERE(MROW,MSA,I)*ZV2M,
     &           VX2PERE(MROW,MSA,I+1)*ZV2P)

               FSUBM(KYPE+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PARE(MROW,MSA,I), VQ1PAREM(MROW,MSA,I))    

               FSUBM(KYPP+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PERE(MROW,MSA,I), VQ1PEREM(MROW,MSA,I))    

               GSUBM(KYPE+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PARE(MROW,MSA,I+1), VQ1PAREM(MROW,MSA,I))    

               GSUBM(KYPP+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PERE(MROW,MSA,I+1), VQ1PEREM(MROW,MSA,I))    

               DSUBM(KYPE+LYROW, KYB2+LYCOL,I)=
     &           GG(VQ2PAREM(MROW,MSA,I), VQ2PARE(MROW,MSA,I),
     &           VQ2PARE(MROW,MSA,I+1))

               DSUBM(KYPP+LYROW, KYB2+LYCOL,I)=
     &           GG(VQ2PEREM(MROW,MSA,I), VQ2PERE(MROW,MSA,I),
     &           VQ2PERE(MROW,MSA,I+1))

               DSUBM(KYPE+LYROW, KYB3+LYCOL,I)=
     &           GG(VQ3PAREM(MROW,MSA,I), VQ3PARE(MROW,MSA,I)*ZB3M,
     &           VQ3PARE(MROW,MSA,I+1)*ZB3P)

               DSUBM(KYPP+LYROW, KYB3+LYCOL,I)=
     &           GG(VQ3PEREM(MROW,MSA,I), VQ3PERE(MROW,MSA,I)*ZB3M,
     &           VQ3PERE(MROW,MSA,I+1)*ZB3P)

               IF (KXDPHI.GT.0) THEN
               FSUBM(KYPE+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPARE(MROW,MSA,I), VDPPAREM(MROW,MSA,I))    

               FSUBM(KYPP+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPERE(MROW,MSA,I), VDPPEREM(MROW,MSA,I))    

               GSUBM(KYPE+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPARE(MROW,MSA,I+1), VDPPAREM(MROW,MSA,I))    

               GSUBM(KYPP+LYROW, KXDPHI+LXCOL,I)=
     &           GF(VDPPERE(MROW,MSA,I+1), VDPPEREM(MROW,MSA,I))    
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDIF

C     FILL IN MATRICES FOR KXDPHI FOR NEUTRALITY EQUATION
      IF (KXDPHI.GT.0) THEN
      DO I=2,NRP1
         INCLUDE 'tent.inc'
 
         ZV2M = (CS(I)/CSM(I-1))**IEXV2
         ZV2P = (CS(I)/CSM(I  ))**IEXV2
         ZB3M = (CSM(I-1)/CS(I))**IEXB3
         ZB3P = (CSM(I)  /CS(I))**IEXB3

         DO MROW=1,MSMAX
            LXROW = (MROW-1)*NXCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               BSUBM(KXDPHI+LXROW,KXX1+LXCOL,I) =
     &            FF(VX1DPHI(MROW,MSA,I),VX1DPHIM(MROW,MSA,I-1),
     &               VX1DPHIM(MROW,MSA,I))
               ASUBM(KXDPHI+LXROW,KXX1+LXCOL,I) = 
     &            FFM(VX1DPHIM(MROW,MSA,I-1))
               CSUBM(KXDPHI+LXROW,KXX1+LXCOL,I) = 
     &            FFP(VX1DPHIM(MROW,MSA,I))

               HSUBM(KXDPHI+LXROW,KYX2+LYCOL,I) =
     &            FGM(VX2DPHI(MROW,MSA,I)*ZV2M,VX2DPHIM(MROW,MSA,I-1))
               ESUBM(KXDPHI+LXROW,KYX2+LYCOL,I) =
     &            FGP(VX2DPHI(MROW,MSA,I)*ZV2P,VX2DPHIM(MROW,MSA,I))

               BSUBM(KXDPHI+LXROW,KXB1+LXCOL,I) =
     &            FF(VQ1DPHI(MROW,MSA,I),VQ1DPHIM(MROW,MSA,I-1),
     &               VQ1DPHIM(MROW,MSA,I))
               ASUBM(KXDPHI+LXROW,KXB1+LXCOL,I) = 
     &            FFM(VQ1DPHIM(MROW,MSA,I-1))
               CSUBM(KXDPHI+LXROW,KXB1+LXCOL,I) = 
     &            FFP(VQ1DPHIM(MROW,MSA,I))

               HSUBM(KXDPHI+LXROW,KYB2+LYCOL,I) =
     &            FGM(VQ2DPHI(MROW,MSA,I),VQ2DPHIM(MROW,MSA,I-1))
               ESUBM(KXDPHI+LXROW,KYB2+LYCOL,I) =
     &            FGP(VQ2DPHI(MROW,MSA,I),VQ2DPHIM(MROW,MSA,I))

               HSUBM(KXDPHI+LXROW,KYB3+LYCOL,I) =
     &            FGM(VQ3DPHI(MROW,MSA,I)*ZB3M,VQ3DPHIM(MROW,MSA,I-1))
               ESUBM(KXDPHI+LXROW,KYB3+LYCOL,I) =
     &            FGP(VQ3DPHI(MROW,MSA,I)*ZB3P,VQ3DPHIM(MROW,MSA,I))

               BSUBM(KXDPHI+LXROW,KXDPHI+LXCOL,I) =
     &            FF(VDPDPHI(MROW,MSA,I),VDPDPHIM(MROW,MSA,I-1),
     &               VDPDPHIM(MROW,MSA,I))
               ASUBM(KXDPHI+LXROW,KXDPHI+LXCOL,I)=
     &            FFM(VDPDPHIM(MROW,MSA,I-1))
               CSUBM(KXDPHI+LXROW,KXDPHI+LXCOL,I)=
     &            FFP(VDPDPHIM(MROW,MSA,I))
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      ENDIF
      
      IF (ISMPIRUN.EQ.0) THEN
      IF (KCHECK.EQ.1) THEN
        MROW = 1
        WRITE(*,*) 'KJP:MROW=',MROW
        DO MSA=1,MSMAX
        DO I=JS0,JS0
          WRITE(*,1101) CS(I),
     &                  VX1PARA(MROW,MSA,I),VX1PERP(MROW,MSA,I), 
     &                  VX2PARA(MROW,MSA,I),VX2PERP(MROW,MSA,I), 
     &                  VQ1PARA(MROW,MSA,I),VQ1PERP(MROW,MSA,I), 
     &                  VQ2PARA(MROW,MSA,I),VQ2PERP(MROW,MSA,I), 
     &                  VQ3PARA(MROW,MSA,I),VQ3PERP(MROW,MSA,I),
     &                  VQ3PARAM(MROW,MSA,I),VQ3PERPM(MROW,MSA,I), 
     &                  VDPPARA(MROW,MSA,I),VDPPERP(MROW,MSA,I),
     &                  VDPPARAM(MROW,MSA,I),VDPPERPM(MROW,MSA,I) 
        ENDDO
        ENDDO
 1101   FORMAT(33E16.8)
      ENDIF
      ENDIF
      
 1    CONTINUE

      IF (KJPKEY.EQ.1.OR.KJPKEY.EQ.2) GOTO 2

      DEALLOCATE( RLM )
      DEALLOCATE( NLAMK0, NLAMK1 )
      DEALLOCATE( HKMIN, HKMAX, WFUN )
      DEALLOCATE( KNUMDISTRIB )
      DEALLOCATE( NUEFF,OMEGASN,OMEGAST,
     &            AAK,EPSLONCA,EPSALPHA,OMEGASNA,OMEGASAA,OMEGASCA,
     &            ALPHAA1,ALPHAA2,ALPHAA3,ZC1 )
      DEALLOCATE( FREQK ) 
      IF (ISMPIRUN.EQ.0) THEN
         DEALLOCATE( FREQKSURF,POSITIONKAI )
      ENDIF

      DEALLOCATE( ESPECIES_DENF,ESPECIES_PREF,ESPECIES_PREP,
     &            ESPECIES_DEN,ESPECIES_PRE,ESPECIES_TEM,ISPECIES_EK )
      DEALLOCATE( BK,HK,BPK, LAMK0,LAMK1 )
      DEALLOCATE( VX1PARA, VX1PARAM,
     &            VX1PERP, VX1PERPM,
     &            VX2PARA, VX2PARAM,
     &            VX2PERP, VX2PERPM,
     &            VQ1PARA, VQ1PARAM,
     &            VQ1PERP, VQ1PERPM,
     &            VQ2PARA, VQ2PARAM,
     &            VQ2PERP, VQ2PERPM,
     &            VQ3PARA, VQ3PARAM,
     &            VQ3PERP, VQ3PERPM,
     &            VDPPARA, VDPPARAM,
     &            VDPPERP, VDPPERPM )
      IF (ALLOCATED(VX1PARE))
     &   DEALLOCATE( VX1PARE, VX1PAREM,
     &               VX1PERE, VX1PEREM,
     &               VX2PARE, VX2PAREM,
     &               VX2PERE, VX2PEREM,
     &               VQ1PARE, VQ1PAREM,
     &               VQ1PERE, VQ1PEREM,
     &               VQ2PARE, VQ2PAREM,
     &               VQ2PERE, VQ2PEREM,
     &               VQ3PARE, VQ3PAREM,
     &               VQ3PERE, VQ3PEREM,
     &               VDPPARE, VDPPAREM,
     &               VDPPERE, VDPPEREM )
      IF (INCDPHI.GT.0) DEALLOCATE(
     &            VX1DPHI, VX1DPHIM,
     &            VX2DPHI, VX2DPHIM,
     &            VQ1DPHI, VQ1DPHIM,
     &            VQ2DPHI, VQ2DPHIM,
     &            VQ3DPHI, VQ3DPHIM,
     &            VDPDPHI, VDPDPHIM )
     
      DEALLOCATE( LAMM,LAMHH,LAMTMP )
      DEALLOCATE( RCHI,RCHI2,RW1,RW2,
     &            RJB,RX1P,RX1B,RX1R,RX2,RQ1,RQ2,
     &            RBT,RPHI,RDMU,RDB )
      DEALLOCATE( RCHIK,RPHIK,RTK,
     &            RHK,RJBK,RX1PK,RX1BK,RX1RK,RX2K,
     &            RQ1K,RQ2K,RVALK )
      DEALLOCATE( RVAK1,RVAK2,RVAK3,RVAK4,RVAK5,RVAK6,RVAK7 )
      DEALLOCATE( RCHIN,RVALN,RHN )
      DEALLOCATE( VPARA,VPERP,VDPHI,VX1,VX2,VQ1,VQ2,VQ3,VDP )
      DEALLOCATE( VPARA0,VPERP0,VDPHI0,VX10,VX20,VQ10,VQ20,VQ30,VDP0 )
      DEALLOCATE( SVPARA0,SVPERP0,SVDPHI0,SVX10,SVX20,SVQ10,SVQ20,
     &            SVQ30,SVDP0 )
      DEALLOCATE( SLAM0,SF0 )
      DEALLOCATE( VI )

      CALL DEALLOCATEDWKCOMPMAT
      
C     IF (IFOWP.EQ.1.OR.IFOWT.EQ.1) THEN
      DEALLOCATE( VI1,VI2,VI3 )
      DEALLOCATE( VPARA1,VPERP1,VDPHI1,VX11,VX21,VQ11,VQ21,VQ31,VDP1 )
      DEALLOCATE( VPARA01,VPERP01,VDPHI01,VX101,VX201,VQ101,VQ201,
     &            VQ301,VDP01 )
      DEALLOCATE( SVPARA01,SVPERP01,SVDPHI01,SVX101,SVX201,SVQ101,
     &            SVQ201,SVQ301,SVDP01 )
      DEALLOCATE( VX1LNP,VX2LNP,VQ1LNP,VQ2LNP,VQ3LNP,VDPLNP )
C     ENDIF

      CALL ZDEALLOCANISO(1)

 2    CONTINUE

      RETURN
      END

C=======================================================================
C COMPUTE X'/X FOR FIRST ORDER FOW CORRECTION
C NOTE THAT X'/X ARE COMPUTED USING X FROM PREVIOUS ITERATION IN THE 
C NON-PERTURBATIVE LOOP OVER THE CONVERGENCE OF EIGENVALUE
C ALSO SET UP AN UPPER LIMIT FOR THE AMPLITUDE |X'/X|, ASSOCIATED WITH 
C ESTIMATED |\DELTA\PSI| VALUE 
C YQL, 06-2013 
C=======================================================================
      SUBROUTINE KXLNP(JS,KGRID)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE KINETICM
      IMPLICIT NONE

      INTEGER    JS,KGRID,M,L
      REAL*8     DELTA0,DFRAC0,DTHRESH,RTMP,HJ,HJ1,HHJ,HHJ1,H0
      COMPLEX*16 CX,CX1,CX2,CXP,CXLNP
      PARAMETER  (DFRAC0=100.0)

      INTEGER KCHECK
      KCHECK=0

C     NO NEED TO EXECUTE KXLNP NEAR RATIONAL SURFACES
      IF (DELRATS.GT.0..AND.DELRATS.LT.1.) THEN
         IF (KGRID.EQ.1) H0 = CS (JS)
         IF (KGRID.EQ.2) H0 = CSM(JS)
         DO M=1,NRATSURF
            L = IRATSURF(M)
            IF (ABS(H0-CS(L)).LE.DELRATS/Q(L)) THEN      
               VX1LNP = 0.
               VX2LNP = 0.
               VQ1LNP = 0.
               VQ2LNP = 0.
               VQ3LNP = 0.
               VDPLNP = 0.
               RETURN
            ENDIF
         ENDDO
      ENDIF

C     ESTIMATE DELTA0 BASED ON THERMAL IONS
      IF (KGRID.EQ.1) RTMP = T (JS)
      IF (KGRID.EQ.2) RTMP = TM(JS)
      DELTA0 = 2.*SQRT(2.*ESPECIES_TEM(JS,KGRID,1))*RTMP/OMEGACI0
      IF (DELTA0.LT.1.0E-4) DELTA0=1.0E-4

      DTHRESH = DFRAC0/DELTA0

      IF (KGRID.EQ.1) H0 = 1./DPSIDS (JS)
      IF (KGRID.EQ.2) H0 = 1./DPSIDSM(JS)

C     EVALUATE MESH SIZES
      HJ   = CS(JS+1)-CS(JS)
      HJ1  = CS(JS)-CS(JS-1)
      HHJ  = (CS(JS+2)-CS(JS))/2.
      HHJ1 = (CS(JS+1)-CS(JS-1))/2.

C     FOR V1 AND V2
      IF (V2XKEY.EQ.0 .OR. V2XKEY.EQ.2) THEN
         DO M=1,MSMAX
         IF (KGRID.EQ.1) THEN
            CX    = V1U(JS,M)
            CX1   = V1U(JS-1,M)
            CX2   = V1U(JS+1,M)
            CXP   = ((HJ1**2*CX2-HJ**2*CX1)/(HJ+HJ1)+(HJ-HJ1)*CX)/HJ/HJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: V1U AT (JS,M)=',JS,M
            ENDIF 
            VX1LNP(M) = CXLNP

            CX  = V2U(JS,M)
            CX1 = V2U(JS-1,M)
            CXP = (CX-CX1)/HHJ1
            CX2 = (CX+CX1)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: V2U AT (JS,M)=',JS,M
            ENDIF 
            VX2LNP(M) = CXLNP
         ELSE
            CX    = V1U(JS,M)
            CX1   = V1U(JS+1,M)
            CXP   = (CX1-CX)/HJ
            CX2   = (CX1+CX)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: V1U AT (JS,M)=',JS,M
            ENDIF 
            VX1LNP(M) = CXLNP
            
            CX    = V2U(JS,M)
            CX1   = V2U(JS-1,M)
            CX2   = V2U(JS+1,M)
            CXP   = ((HHJ1**2*CX2-HHJ**2*CX1)/(HHJ+HHJ1)+(HHJ-HHJ1)*CX)
     &              /HHJ/HHJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: V2U AT (JS,M)=',JS,M
            ENDIF 
            VX2LNP(M) = CXLNP
         ENDIF
         ENDDO

C     FOR X1 AND X2
      ELSE
         DO M=1,MSMAX
         IF (KGRID.EQ.1) THEN
            CX    = X1U(JS,M)
            CX1   = X1U(JS-1,M)
            CX2   = X1U(JS+1,M)
            CXP   = ((HJ1**2*CX2-HJ**2*CX1)/(HJ+HJ1)+(HJ-HJ1)*CX)/HJ/HJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: X1U AT (JS,M)=',JS,M
            ENDIF 
            VX1LNP(M) = CXLNP

            CX  = X2U(JS,M)
            CX1 = X2U(JS-1,M)
            CXP = (CX-CX1)/HHJ1
            CX2 = (CX+CX1)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: X2U AT (JS,M)=',JS,M
            ENDIF 
            VX2LNP(M) = CXLNP
         ELSE
            CX    = X1U(JS,M)
            CX1   = X1U(JS+1,M)
            CXP   = (CX1-CX)/HJ
            CX2   = (CX1+CX)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: X1U AT (JS,M)=',JS,M
            ENDIF 
            VX1LNP(M) = CXLNP
            
            CX    = X2U(JS,M)
            CX1   = X2U(JS-1,M)
            CX2   = X2U(JS+1,M)
            CXP   = ((HHJ1**2*CX2-HHJ**2*CX1)/(HHJ+HHJ1)+(HHJ-HHJ1)*CX)
     &              /HHJ/HHJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: X2U AT (JS,M)=',JS,M
            ENDIF 
            VX2LNP(M) = CXLNP
         ENDIF
         ENDDO
      ENDIF

C     FOR Q1,Q2,Q3
      DO M=1,MSMAX
         IF (KGRID.EQ.1) THEN
            CX    = B1U(JS,M)
            CX1   = B1U(JS-1,M)
            CX2   = B1U(JS+1,M)
            CXP   = ((HJ1**2*CX2-HJ**2*CX1)/(HJ+HJ1)+(HJ-HJ1)*CX)/HJ/HJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: B1U AT (JS,M)=',JS,M
            ENDIF 
            VQ1LNP(M) = CXLNP

            CX  = B2U(JS,M)
            CX1 = B2U(JS-1,M)
            CXP = (CX-CX1)/HHJ1
            CX2 = (CX+CX1)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: B2U AT (JS,M)=',JS,M
            ENDIF 
            VQ2LNP(M) = CXLNP

            CX  = B3U(JS,M)
            CX1 = B3U(JS-1,M)
            CXP = (CX-CX1)/HHJ1
            CX2 = (CX+CX1)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: B3U AT (JS,M)=',JS,M
            ENDIF 
            VQ3LNP(M) = CXLNP

            CX    = DPHI(JS,M)
            CX1   = DPHI(JS-1,M)
            CX2   = DPHI(JS+1,M)
            CXP   = ((HJ1**2*CX2-HJ**2*CX1)/(HJ+HJ1)+(HJ-HJ1)*CX)/HJ/HJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: DPHI AT (JS,M)=',JS,M
            ENDIF 
            VDPLNP(M) = CXLNP
         ELSE
            CX    = B1U(JS,M)
            CX1   = B1U(JS+1,M)
            CXP   = (CX1-CX)/HJ
            CX2   = (CX1+CX)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: B1U AT (JS,M)=',JS,M
            ENDIF 
            VQ1LNP(M) = CXLNP
            
            CX    = B2U(JS,M)
            CX1   = B2U(JS-1,M)
            CX2   = B2U(JS+1,M)
            CXP   = ((HHJ1**2*CX2-HHJ**2*CX1)/(HHJ+HHJ1)+(HHJ-HHJ1)*CX)
     &              /HHJ/HHJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: B2U AT (JS,M)=',JS,M
            ENDIF 
            VQ2LNP(M) = CXLNP
            
            CX    = B3U(JS,M)
            CX1   = B3U(JS-1,M)
            CX2   = B3U(JS+1,M)
            CXP   = ((HHJ1**2*CX2-HHJ**2*CX1)/(HHJ+HHJ1)+(HHJ-HHJ1)*CX)
     &              /HHJ/HHJ1
            CXLNP = 0.
            IF (ABS(CX).GT.0.) CXLNP = H0*CXP/CX
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: B3U AT (JS,M)=',JS,M
            ENDIF 
            VQ3LNP(M) = CXLNP

            CX    = DPHI(JS,M)
            CX1   = DPHI(JS+1,M)
            CXP   = (CX1-CX)/HJ
            CX2   = (CX1+CX)/2.
            CXLNP = 0.
            IF (ABS(CX2).GT.0.) CXLNP = H0*CXP/CX2
            RTMP = ABS(CXLNP)
            IF (RTMP.GT.DTHRESH) THEN
               CXLNP = CXLNP/RTMP*DTHRESH
               WRITE(*,*) 'KXLNP WARNING: DPHI AT (JS,M)=',JS,M
            ENDIF 
            VDPLNP(M) = CXLNP
         ENDIF
      ENDDO

      RETURN
      END

C=======================================================================
C COMPUTE RHS DURING INVERSE ITERATION                                 = 
C FOR KINETIC TERMS ASSOCIATED WITH B1,B2,B3,DPHI                           =
C YQL, 03-2008                                                         =
C=======================================================================
      SUBROUTINE KJPDX(MD, MDY, ND, R, RY, X, Y)
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE KINETICM
      IMPLICIT NONE

      INTEGER         MD, MDY, ND
      COMPLEX*16      R(MD,ND+1), X(MD,ND+1), RY(MDY,ND), Y(MDY,ND)
      INTEGER MROW,MSA,I,LYROW,LXCOL,LYCOL
      REAL*8  ZB3M,ZB3P
      COMPLEX*16 TTMP
      INTEGER IEXB3            
      PARAMETER (IEXB3=1)
      INTEGER KCHECK

      INCLUDE 'integc.inc'

      KCHECK=0

      TTMP = -CALPHA1*CALPHA4/CALPHA2**2
      DO I=1,NR
         INCLUDE 'tophat.inc'

         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         DO MROW=1,MSMAX
            LYROW = (MROW-1)*NYCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               RY(KYPPARA+LYROW,I)=RY(KYPPARA+LYROW,I) + TTMP*(
     &          -GF(VQ1PARA(MROW,MSA,I),VQ1PARAM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I)    
     &          -GF(VQ1PARA(MROW,MSA,I+1),VQ1PARAM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I+1)    
     &          -GG(VQ2PARAM(MROW,MSA,I),VQ2PARA(MROW,MSA,I),
     &              VQ2PARA(MROW,MSA,I+1))*Y(KYB2+LYCOL,I)
     &          -GG(VQ3PARAM(MROW,MSA,I),VQ3PARA(MROW,MSA,I)*ZB3M,
     &              VQ3PARA(MROW,MSA,I+1)*ZB3P)*Y(KYB3+LYCOL,I))

               IF (KYPE.GT.0)
     &         RY(KYPE+LYROW,I)=RY(KYPE+LYROW,I) + TTMP*(
     &          -GF(VQ1PARE(MROW,MSA,I),VQ1PAREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I)    
     &          -GF(VQ1PARE(MROW,MSA,I+1),VQ1PAREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I+1)    
     &          -GG(VQ2PAREM(MROW,MSA,I),VQ2PARE(MROW,MSA,I),
     &              VQ2PARE(MROW,MSA,I+1))*Y(KYB2+LYCOL,I)
     &          -GG(VQ3PAREM(MROW,MSA,I),VQ3PARE(MROW,MSA,I)*ZB3M,
     &              VQ3PARE(MROW,MSA,I+1)*ZB3P)*Y(KYB3+LYCOL,I))

               IF (KXDPHI.GT.0) 
     &         RY(KYPPARA+LYROW,I)=RY(KYPPARA+LYROW,I) + TTMP*(
     &          -GF(VDPPARA(MROW,MSA,I),VDPPARAM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I)    
     &          -GF(VDPPARA(MROW,MSA,I+1),VDPPARAM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I+1))    

               IF (KXDPHI.GT.0.AND.KYPE.GT.0) 
     &         RY(KYPE+LYROW,I)=RY(KYPE+LYROW,I) + TTMP*(
     &          -GF(VDPPARE(MROW,MSA,I),VDPPAREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I)    
     &          -GF(VDPPARE(MROW,MSA,I+1),VDPPAREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I+1))    

               RY(KYPPERP+LYROW,I)=RY(KYPPERP+LYROW,I) + TTMP*(
     &          -GF(VQ1PERP(MROW,MSA,I),VQ1PERPM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I)    
     &          -GF(VQ1PERP(MROW,MSA,I+1),VQ1PERPM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I+1)    
     &          -GG(VQ2PERPM(MROW,MSA,I),VQ2PERP(MROW,MSA,I),
     &              VQ2PERP(MROW,MSA,I+1))*Y(KYB2+LYCOL,I)
     &          -GG(VQ3PERPM(MROW,MSA,I),VQ3PERP(MROW,MSA,I)*ZB3M,
     &              VQ3PERP(MROW,MSA,I+1)*ZB3P)*Y(KYB3+LYCOL,I))

               IF (KYPP.GT.0) 
     &         RY(KYPP+LYROW,I)=RY(KYPP+LYROW,I) + TTMP*(
     &          -GF(VQ1PERE(MROW,MSA,I),VQ1PEREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I)    
     &          -GF(VQ1PERE(MROW,MSA,I+1),VQ1PEREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I+1)    
     &          -GG(VQ2PEREM(MROW,MSA,I),VQ2PERE(MROW,MSA,I),
     &              VQ2PERE(MROW,MSA,I+1))*Y(KYB2+LYCOL,I)
     &          -GG(VQ3PEREM(MROW,MSA,I),VQ3PERE(MROW,MSA,I)*ZB3M,
     &              VQ3PERE(MROW,MSA,I+1)*ZB3P)*Y(KYB3+LYCOL,I))

               IF (KXDPHI.GT.0) 
     &         RY(KYPP+LYROW,I)=RY(KYPP+LYROW,I) + TTMP*(
     &          -GF(VDPPERP(MROW,MSA,I),VDPPERPM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I)    
     &          -GF(VDPPERP(MROW,MSA,I+1),VDPPERPM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I+1))   

               IF (KXDPHI.GT.0.AND.KYPP.GT.0) 
     &         RY(KYPP+LYROW,I)=RY(KYPP+LYROW,I) + TTMP*(
     &          -GF(VDPPERE(MROW,MSA,I),VDPPEREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I)    
     &          -GF(VDPPERE(MROW,MSA,I+1),VDPPEREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I+1))   
            ENDDO
         ENDDO
      ENDDO

      IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &    .AND.INCKIN.GT.0) THEN
      DO I=1,NR
         INCLUDE 'tophat.inc'

         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         DO MROW=1,MSMAX
            LYROW = (MROW-1)*NYCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               RY(KYPE+LYROW,I)=RY(KYPE+LYROW,I) + TTMP*(
     &          -GF(VQ1PARE(MROW,MSA,I),VQ1PAREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I)    
     &          -GF(VQ1PARE(MROW,MSA,I+1),VQ1PAREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I+1)    
     &          -GG(VQ2PAREM(MROW,MSA,I),VQ2PARE(MROW,MSA,I),
     &              VQ2PARE(MROW,MSA,I+1))*Y(KYB2+LYCOL,I)
     &          -GG(VQ3PAREM(MROW,MSA,I),VQ3PARE(MROW,MSA,I)*ZB3M,
     &              VQ3PARE(MROW,MSA,I+1)*ZB3P)*Y(KYB3+LYCOL,I))

               IF (KXDPHI.GT.0) 
     &         RY(KYPE+LYROW,I)=RY(KYPE+LYROW,I) + TTMP*(
     &          -GF(VDPPARE(MROW,MSA,I),VDPPAREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I)    
     &          -GF(VDPPARE(MROW,MSA,I+1),VDPPAREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I+1))    

               RY(KYPP+LYROW,I)=RY(KYPP+LYROW,I) + TTMP*(
     &          -GF(VQ1PERE(MROW,MSA,I),VQ1PEREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I)    
     &          -GF(VQ1PERE(MROW,MSA,I+1),VQ1PEREM(MROW,MSA,I))*
     &           X(KXB1+LXCOL,I+1)    
     &          -GG(VQ2PEREM(MROW,MSA,I),VQ2PERE(MROW,MSA,I),
     &              VQ2PERE(MROW,MSA,I+1))*Y(KYB2+LYCOL,I)
     &          -GG(VQ3PEREM(MROW,MSA,I),VQ3PERE(MROW,MSA,I)*ZB3M,
     &              VQ3PERE(MROW,MSA,I+1)*ZB3P)*Y(KYB3+LYCOL,I))

               IF (KXDPHI.GT.0) 
     &         RY(KYPP+LYROW,I)=RY(KYPP+LYROW,I) + TTMP*(
     &          -GF(VDPPERE(MROW,MSA,I),VDPPEREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I)    
     &          -GF(VDPPERE(MROW,MSA,I+1),VDPPEREM(MROW,MSA,I))*
     &           X(KXDPHI+LXCOL,I+1))   
            ENDDO
         ENDDO
      ENDDO
      ENDIF

      RETURN
      END

C=======================================================================
C PREPARE EQUILIBRIUM FOURIER COEFFICIENTS FOR SUBROUTINE KCOEFFI      = 
C SHOULD BE CALLED FROM WITHIN FTCOEFF                                 =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KFTCOEFF

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      USE KINETICM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INCLUDE 'comfft.inc'
C
      INTEGER NPSTRT
      INTEGER I,J
      REAL*8  H1,H2
      REAL*8,DIMENSION(:,:),ALLOCATABLE::RW0,RW0M,
     &                                   B_2,B_2M,B2S,B2SM,B2C,B2CM

      ALLOCATE( RW0(NRP1,NCHI), RW0M(NRP1,NCHI), 
     &          B_2(NRP1,NCHI), B_2M(NRP1,NCHI),
     &          B2S(NRP1,NCHI), B2SM(NRP1,NCHI),
     &          B2C(NRP1,NCHI), B2CM(NRP1,NCHI) )
C
      INCLUDE 'setfft.inc'
      SUBNAM    = 'KFTCOEFF'

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
         FFF(:,1) = B_2(1:NRP1,J)
         FFF(:,2) = B_2M(1:NRP1,J)
         CALL DFFFDPSI(0)
         B2S(1:NRP1,J)  = DFFF(:,1)*DPSIDS(1:NRP1)
         B2SM(1:NRP1,J) = DFFF(:,2)*DPSIDSM(1:NRP1)
      ENDDO

      IF (.NOT.ALLOCATED(V1PK)) THEN
         ALLOCATE( V1PK(NRP1,MEDIM),   V1PKM(NRP1,MEDIM),
     &             V2PK(NRP1,MEDIM),   V2PKM(NRP1,MEDIM),  
     &             V3PK(NRP1,MEDIM),   V3PKM(NRP1,MEDIM),  
     &             VKOX1A(NRP1,MEDIM), VKOX1AM(NRP1,MEDIM),  
     &             VKOX1E(NRP1,MEDIM), VKOX1EM(NRP1,MEDIM),  
     &             VKOQ1(NRP1,MEDIM),  VKOQ1M(NRP1,MEDIM),  
     &             VKOQ2(NRP1,MEDIM),  VKOQ2M(NRP1,MEDIM),  
     &             VKOQ3(NRP1,MEDIM),  VKOQ3M(NRP1,MEDIM) )
      ENDIF

C     V1PK
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=G12L(I,J)*DPSIDS(I)**2*B2C(I,J)/2/RJA(I,J)/
     &              B_2(I,J)**2 - RJA(I,J)*(PPEQ(I)*DPSIDS(I)+
     &              B2S(I,J)/2)/B_2(I,J)
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=G12LM(I,J)*DPSIDSM(I)**2*B2CM(I,J)/2/RJAM(I,J)/
     &              B_2M(I,J)**2 - RJAM(I,J)*(PPEQM(I)*DPSIDSM(I)+
     &              B2SM(I,J)/2)/B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,  V1PK,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  V1PK    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      1,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, V1PKM,   FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  V1PKM   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      2,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    V1PK,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'V1PK')
        CALL FFTOUTPT(RW0M,   V1PKM,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'V1PKM')
      ENDIF
C
C     V2PK
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=-RJA(I,J)*T(I)*B2C(I,J)/2./B_2(I,J)**2
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=-RJAM(I,J)*TM(I)*B2CM(I,J)/2./B_2M(I,J)**2
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,  V2PK,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  V2PK    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      3,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, V2PKM,   FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  V2PKM   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      4,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    V2PK,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'V2PK')
        CALL FFTOUTPT(RW0M,   V2PKM,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'V2PKM')
      ENDIF
C
C     V3PK
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=DPSIDS(I)*B2C(I,J)/2/B_2(I,J)
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=DPSIDSM(I)*B2CM(I,J)/2/B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,  V3PK,    FORWD, NRP1,  NRP1,   NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  V3PK    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      5,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, V3PKM,   FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  V3PKM   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      6,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    V3PK,      NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'V3PK')
        CALL FFTOUTPT(RW0M,   V3PKM,     NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'V3PKM')
      ENDIF
C
C     VKOX1A
      DO J=1,NCHI
        DO I=2,NRP1        
           RW0(I,J)=RJA(I,J)*PEQ(I)/B_2(I,J)*
     &     (
     &              (PCOEFDA-PCOEFKA)*PPEQ(I)*DPSIDS(I)
     &             +(PCOEFDA-PCOEFKA*0.5)*(B2S(I,J)-B2C(I,J)
     &             *DPSIDS(I)**2*G12L(I,J)/RJA(I,J)**2/B_2(I,J))
     &     )
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=RJAM(I,J)*PEQM(I)/B_2M(I,J)*
     &     (
     &              (PCOEFDA-PCOEFKA)*PPEQM(I)*DPSIDSM(I)
     &             +(PCOEFDA-PCOEFKA*0.5)*(B2SM(I,J)-B2CM(I,J)
     &             *DPSIDSM(I)**2*G12LM(I,J)/RJAM(I,J)**2/B_2M(I,J))
     &     )     
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,   VKOX1A, FORWD, NRP1,  NRP1,   NPSTRT
     &                      ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOX1A   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      7,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, VKOX1AM, FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOX1AM  in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      8,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    VKOX1A,    NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOX1A')
        CALL FFTOUTPT(RW0M,   VKOX1AM,   NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOX1AM')
      ENDIF
C
C     VKOX1E
      DO J=1,NCHI
        DO I=2,NRP1        
           RW0(I,J)=RJA(I,J)*PEQ(I)/B_2(I,J)*
     &     (
     &              (PCOEFDE-PCOEFKE)*PPEQ(I)*DPSIDS(I)
     &             +(PCOEFDE-PCOEFKE*0.5)*(B2S(I,J)-B2C(I,J)*
     &              DPSIDS(I)**2*G12L(I,J)/RJA(I,J)**2/B_2(I,J))
     &     )     
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=RJAM(I,J)*PEQM(I)/B_2M(I,J)*
     &     (
     &              (PCOEFDE-PCOEFKE)*PPEQM(I)*DPSIDSM(I)
     &             +(PCOEFDE-PCOEFKE*0.5)*(B2SM(I,J)-B2CM(I,J)*
     &              DPSIDSM(I)**2*G12LM(I,J)/RJAM(I,J)**2/B_2M(I,J))
     &     )     
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,   VKOX1E, FORWD, NRP1,  NRP1,   NPSTRT
     &                      ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOX1E   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,      9,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, VKOX1EM, FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOX1EM  in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,     10,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    VKOX1E,    NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOX1E')
        CALL FFTOUTPT(RW0M,   VKOX1EM,   NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOX1EM')
      ENDIF
C
C     VKOQ1
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=PEQ(I)*DPSIDS(I)*G12L(I,J)/RJA(I,J)/B_2(I,J)
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=PEQM(I)*DPSIDSM(I)*G12LM(I,J)/RJAM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,   VKOQ1,  FORWD, NRP1,  NRP1,   NPSTRT
     &                      ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOQ1    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,     11,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, VKOQ1M,  FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOQ1M   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,     12,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    VKOQ1,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOQ1')
        CALL FFTOUTPT(RW0M,   VKOQ1M,    NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOQ1M')
      ENDIF
C
C     VKOQ2
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=PEQ(I)*DPSIDS(I)*G22L(I,J)/RJA(I,J)/B_2(I,J)
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=PEQM(I)*DPSIDSM(I)*G22LM(I,J)/RJAM(I,J)/B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,   VKOQ2,  FORWD, NRP1,  NRP1,   NPSTRT
     &                      ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOQ2    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,     13,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, VKOQ2M,  FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOQ2M   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,     14,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    VKOQ2,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOQ2')
        CALL FFTOUTPT(RW0M,   VKOQ2M,    NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOQ2M')
      ENDIF
C
C     VKOQ3
      DO J=1,NCHI
        DO I=2,NRP1
           RW0(I,J)=PEQ(I)*T(I)/B_2(I,J)
        ENDDO
        RW0(1,J) = 0.0
        DO I=1,NR
           RW0M(I,J)=PEQM(I)*TM(I)/B_2M(I,J)
        ENDDO
      ENDDO
C
      NPSTRT    =  1
      call FFTDRIVER( RW0,   VKOQ3,  FORWD, NRP1,  NRP1,   NPSTRT
     &                      ,MEDIM,  NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOQ3    in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,     15,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      call FFTDRIVER( RW0M, VKOQ3M,  FORWD, NRP1,  NR,     NPSTRT
     &                     ,MEDIM,   NCHI,  KUOUT, IERSUB, IERPLC, IERR)
      if(IERR .NE. 0) THEN
         write(MESSAGE,*) 'Error for  VKOQ3M   in ',IERSUB
         call ABORTRUN
     &        (SUBNAM,     16,   MESSAGE
     &        ,'IERR    ', IERR,     'IERPLC  ', IERPLC,   -1, KUOUT)
      endif
C
      IF(KUFFTP .GT. 0) THEN
        CALL FFTOUTPT(RW0,    VKOQ3,     NRP1,    NRP1,    NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOQ3')
        CALL FFTOUTPT(RW0M,   VKOQ3M,    NRP1,    NR,      NPSTRT
     &                       ,MEDIM,     NCHI,    KUFFTP, 'VKOQ3M')
      ENDIF
C
      DEALLOCATE( RW0,RW0M,B_2,B_2M,B2S,B2SM,B2C,B2CM )

      RETURN
      END

C=======================================================================
C FILL IN MATRICES DUE TO FLUID PART OF KINETIC PRESSURE               =
C HERE WE INCLUDE:                                                     = 
C  1) CONVOLUTION TERMS FOR LHS (JP) OF KINETIC PRESSURE EQUATIONS     = 
C  2) GRAD(P) DUE TO KINETIC TERMS                                     =
C  3) CGL LIMIT                                                        =
C SHOULD BE CALLED FROM WITHIN COEFFI(...)                             =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KCOEFFI(MROW,MSA,MSB,CMROW,CMA,CMB,CNA,SHIFTC,SHIFTM,
     &           ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE CONVOLCOFM
      USE KINETICM
      IMPLICIT NONE
      INCLUDE 'specmat.inc'

      INTEGER   MROW,MSA,MSB,NSA,NSB
      PARAMETER (NSA=2,NSB=1)
      INTEGER   LXCOL,LYCOL,LXROW,LYROW,I
      REAL*8    RTMP,ZEM,ZEP,ZV2M,ZV2P,ZB3M,ZB3P
      INTEGER   IEXE,IEXV2,IEXB3            
      PARAMETER (IEXV2=-1, IEXB3=1)
      COMPLEX*16 CMROW,CMA,CMB,CNA,CTMP
      COMPLEX*16,DIMENSION(NRP1)::SHIFTC,SHIFTM

      INTEGER KCHECK

      INCLUDE 'integc.inc'

      KCHECK=1

      LXROW = (MROW-1)*NXCOMP
      LXCOL = (MSA -1)*NXCOMP
      LYROW = (MROW-1)*NYCOMP
      LYCOL = (MSA -1)*NYCOMP

      IF (KPBKEY.EQ.0) GOTO 1010

C-----------------------------------------------------------------------
C.. FIRST EQUATION: COVARIANT-S-COMP OF EQ. OF MOTION
C..                              (KXV1 = 1, DEFINED ON INTEGER MESH)
C-----------------------------------------------------------------------
      DO 10 I=2,NRP1
      INCLUDE 'tent.inc'
      
      HSUBM(KXV1+LXROW, KYPPERP+LYCOL,I)= 
     $       ZNORM * JACOBI(I,MSB)
     $     + CMA *FGM(G12B2B2(I,MSB), G12B2B2M(I-1,MSB))
     $     + CNA *FGM(G12B2B3(I,MSB), G12B2B3M(I-1,MSB))
     $     - FGM(V1PK(I,MSB), V1PKM(I-1,MSB))
      ESUBM(KXV1+LXROW, KYPPERP+LYCOL,I)=
     $     - ZNORM * JACOBI(I,MSB)
     $     + CMA *FGP(G12B2B2(I,MSB), G12B2B2M(I,MSB)) 
     $     + CNA *FGP(G12B2B3(I,MSB), G12B2B3M(I,MSB))
     $     - FGP(V1PK(I,MSB), V1PKM(I,MSB))
      HSUBM(KXV1+LXROW,KYPPARA+LYCOL,I)=FGM(V1PK(I,MSB),V1PKM(I-1,MSB))
      ESUBM(KXV1+LXROW,KYPPARA+LYCOL,I)=FGP(V1PK(I,MSB),V1PKM(I,MSB))
 10   CONTINUE

C-----------------------------------------------------------------------
C.. COEFFICIENTS FOR SECOND EQUATION: COVARIANT-2-COMP. OF EQ. OF MOTION
C..                                  (KYV2 = 1, DEFINED ON HALF MESH)
C-----------------------------------------------------------------------
      IEXE = -1
      DO 20 I=1,NR
      INCLUDE 'tophat.inc'

      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE

      DSUBM(KYV2+LYROW,KYPPERP+LYCOL,I) =
     &     -  CMA*GG(TB2M(I,MSB),ZEM*TB2(I,MSB),ZEP*TB2(I+1,MSB))
     &     +  CNA*GG(G22B2B2M(I,MSB),ZEM*G22B2B2(I,MSB),
     &        ZEP*G22B2B2(I+1,MSB))
     &     -  GG(V2PKM(I,MSB),ZEM*V2PK(I,MSB),ZEP*V2PK(I+1,MSB))

      DSUBM(KYV2+LYROW,KYPPARA+LYCOL,I) =
     &        GG(V2PKM(I,MSB),ZEM*V2PK(I,MSB),ZEP*V2PK(I+1,MSB))
 20   CONTINUE

C-----------------------------------------------------------------------
C.. COEFFICIENTS FOR THIRD  EQUATION: COVARIANT-3-COMP. OF EQ. OF MOTION
C..                                  (KYV3 = 2, DEFINED ON HALF MESH)
C-----------------------------------------------------------------------
      IEXE = -1
      IF (KYV3.GT.0) THEN
      DO 30 I=1,NR
      INCLUDE 'tophat.inc'

      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE

      IF (MSB.EQ.1) CTMP  = 1.0
      IF (MSB.NE.1) CTMP  = 0.0

      DSUBM(KYV3+LYROW,KYPPERP+LYCOL,I)=
     &     - GG(V3PKM(I,MSB),ZEM*V3PK(I,MSB),ZEP*V3PK(I+1,MSB))

      DSUBM(KYV3+LYROW,KYPPARA+LYCOL,I)=
     &     - CMA*GG(CTMP*DPSIDSM(I),CTMP*ZEM*DPSIDS(I),
     &       CTMP*ZEP*DPSIDS(I+1))
     &     - CNA*GG(B3J2M(I,MSB)*TM(I),ZEM*B3J2(I,MSB)*T(I),
     &       ZEP*B3J2(I+1,MSB)*T(I+1))
     &     + GG(V3PKM(I,MSB),ZEM*V3PK(I,MSB),ZEP*V3PK(I+1,MSB))
 30   CONTINUE
      ENDIF

 1010 CONTINUE

C-----------------------------------------------------------------------
C.. FIFTEENTH EQUATION: PARALLEL KINETIC PRESSURE 
C..                     (KYPPARA=9, DEFINED ON HALF MESH)
C-----------------------------------------------------------------------
      DO 110 I=1,NR
      INCLUDE 'tophat.inc'

C     LHS OF THE EQUATION
      IF (IPERTURB.EQ.0 .AND. INCKIN.NE.2 
     &    .AND. (V2XKEY.EQ.0 .OR. V2XKEY.EQ.2) ) THEN
      DSUBM(KYPPARA+LYROW,KYPPARA+LYCOL,I)=
     &     + GG(SHIFTM(I)*JACOBM(I,MSB),SHIFTC(I)*JACOBI(I,MSB),
     &          SHIFTC(I+1)*JACOBI(I+1,MSB))
      ELSE
      DSUBM(KYPPARA+LYROW,KYPPARA+LYCOL,I)=
     &     - GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))
      ENDIF

C     CGL LIMIT              
      IF (INCKIN.EQ.2) THEN
         ZV2M = (CS(I  )/CSM(I))**IEXV2
         ZV2P = (CS(I+1)/CSM(I))**IEXV2
         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         FSUBM(KYPPARA+LYROW,KXX1+LXCOL,I)=ALPHAD*
     &       GF(VKOX1A(I,MSB),VKOX1AM(I,MSB))

         GSUBM(KYPPARA+LYROW,KXX1+LXCOL,I)=ALPHAD*  
     &       GF(VKOX1A(I+1,MSB),VKOX1AM(I,MSB))

         DSUBM(KYPPARA+LYROW,KYX2+LYCOL,I)=-ALPHAD*
     &       2*(PCOEFDA-PCOEFKA*0.5)*
     &       GG(PEQM(I)*V2PKM(I,MSB), PEQ(I)*V2PK(I,MSB)*ZV2M,
     &       PEQ(I+1)*V2PK(I+1,MSB)*ZV2P)

         FSUBM(KYPPARA+LYROW,KXB1+LXCOL,I)=ALPHAD*
     &       PCOEFDA*GF(VKOQ1(I,MSB),VKOQ1M(I,MSB))

         GSUBM(KYPPARA+LYROW,KXB1+LXCOL,I)=ALPHAD*  
     &       PCOEFDA*GF(VKOQ1(I+1,MSB),VKOQ1M(I,MSB))
      
         DSUBM(KYPPARA+LYROW,KYB2+LYCOL,I)=ALPHAD*
     &       PCOEFDA*GG(VKOQ2M(I,MSB),VKOQ2(I,MSB),VKOQ2(I+1,MSB))

         DSUBM(KYPPARA+LYROW,KYB3+LYCOL,I)=ALPHAD*
     &       PCOEFDA*GG(VKOQ3M(I,MSB),VKOQ3(I,MSB)*ZB3M,
     &                  VKOQ3(I+1,MSB)*ZB3P)
      ENDIF
 110  CONTINUE

C-----------------------------------------------------------------------
C.. EQUATIONS FOR PARALLEL (KYPE) AND PERPENDICULAR (KYPP) 
C.. KINETIC PRESSURES FOR THERMAL ELECTRONS 
C..                     (DEFINED ON HALF MESH)
C.. CGL LIMIT NOT YET IMPLEMENTED      
C-----------------------------------------------------------------------
      IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &    .AND.INCKIN.GT.0) THEN
      DO I=1,NR
      INCLUDE 'tophat.inc'

C     LHS OF THE EQUATION
      IF (IPERTURB.EQ.0 .AND. INCKIN.NE.2 
     &    .AND. (V2XKEY.EQ.0 .OR. V2XKEY.EQ.2) ) THEN
      DSUBM(KYPE+LYROW,KYPE+LYCOL,I)=
     &     + GG(SHIFTM(I)*JACOBM(I,MSB),SHIFTC(I)*JACOBI(I,MSB),
     &          SHIFTC(I+1)*JACOBI(I+1,MSB))
      DSUBM(KYPP+LYROW,KYPP+LYCOL,I)=
     &     + GG(SHIFTM(I)*JACOBM(I,MSB),SHIFTC(I)*JACOBI(I,MSB),
     &          SHIFTC(I+1)*JACOBI(I+1,MSB))
      ELSE
      DSUBM(KYPE+LYROW,KYPE+LYCOL,I)=
     &     - GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))
      DSUBM(KYPP+LYROW,KYPP+LYCOL,I)=
     &     - GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))
      ENDIF
      ENDDO
      ENDIF

C-----------------------------------------------------------------------
C.. SIXTEENTH EQUATION: PERP. KINETIC PRESSURE 
C..                     (KYPPERP=10, DEFINED ON HALF MESH)
C-----------------------------------------------------------------------
      DO 120 I=1,NR
      INCLUDE 'tophat.inc'

C     LHS OF THE EQUATION
      IF (IPERTURB.EQ.0 .AND. INCKIN.NE.2 
     &    .AND. (V2XKEY.EQ.0.OR.V2XKEY.EQ.2) ) THEN
      DSUBM(KYPPERP+LYROW,KYPPERP+LYCOL,I)=
     &     + GG(SHIFTM(I)*JACOBM(I,MSB),SHIFTC(I)*JACOBI(I,MSB),
     &          SHIFTC(I+1)*JACOBI(I+1,MSB))
      ELSE
      DSUBM(KYPPERP+LYROW,KYPPERP+LYCOL,I)=
     &     - GG(JACOBM(I,MSB),JACOBI(I,MSB),JACOBI(I+1,MSB))
      ENDIF

C     CGL LIMIT             
      IF (INCKIN.EQ.2) THEN
         ZV2M = (CS(I  )/CSM(I))**IEXV2
         ZV2P = (CS(I+1)/CSM(I))**IEXV2
         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         FSUBM(KYPPERP+LYROW,KXX1+LXCOL,I)=ALPHAD*
     &       GF(VKOX1E(I,MSB),VKOX1EM(I,MSB))

         GSUBM(KYPPERP+LYROW,KXX1+LXCOL,I)=ALPHAD*  
     &       GF(VKOX1E(I+1,MSB),VKOX1EM(I,MSB))
      
         DSUBM(KYPPERP+LYROW,KYX2+LYCOL,I)=-ALPHAD*
     &       2.0*(PCOEFDE-PCOEFKE*0.5)
     &      *GG(PEQM(I)*V2PKM(I,MSB), PEQ(I)*V2PK(I,MSB)*ZV2M,
     &       PEQ(I+1)*V2PK(I+1,MSB)*ZV2P)

         FSUBM(KYPPERP+LYROW,KXB1+LXCOL,I)=ALPHAD*
     &       PCOEFDE*GF(VKOQ1(I,MSB),VKOQ1M(I,MSB))

         GSUBM(KYPPERP+LYROW,KXB1+LXCOL,I)=ALPHAD*  
     &       PCOEFDE*GF(VKOQ1(I+1,MSB),VKOQ1M(I,MSB))
      
         DSUBM(KYPPERP+LYROW,KYB2+LYCOL,I)=ALPHAD*
     &       PCOEFDE*GG(VKOQ2M(I,MSB),VKOQ2(I,MSB),VKOQ2(I+1,MSB))

         DSUBM(KYPPERP+LYROW,KYB3+LYCOL,I)=ALPHAD*
     &       PCOEFDE*GG(VKOQ3M(I,MSB),VKOQ3(I,MSB)*ZB3M,
     &                  VKOQ3(I+1,MSB)*ZB3P)
      ENDIF
 120  CONTINUE

      RETURN
      END

C=======================================================================
C CALCULATE RLM:
C   FOR PRECESSIONAL DRIFT RESONANCE ONLY, INCLUDE ONLY L=0 BOUNCE 
C   HARMONIC; OTHERWISE INCLUDE ALL NECCESARY BOUNCE/TRANSIT HARMONICS;
C   READ FROM A FILE IF EXISTS (FEATURE ADDED BY Z.R.WANG)
C YQL, 06-2013 
C=======================================================================
      SUBROUTINE KGETRLM

      USE GLOBALM
      USE KINETICM
      USE ToolBox
      IMPLICIT NONE

      INTEGER KSINGLEL,K,I,J
      REAL*8  RTMP
      INTEGER HARMOFILE
      INTEGER HARMONUM 
      LOGICAL EXISTS

      KSINGLEL = 1
      RTMP     = SUM(ABS(PSPECIES_NP)) + SUM(ABS(PSPECIES_NTB))
      IF (RTMP.GT.0.) KSINGLEL = 0

      IF (KSINGLEL.EQ.0) THEN

C     I = MAX(ABS(M1),ABS(M2))
C     MLMAX = ABS(M1)+ABS(M2)+1+2*I
C     ALLOCATE( RLM(MLMAX) )
C     DO J=1,MLMAX
C        RLM(J) = -I-ABS(M1)+J-1
C     ENDDO
      MLMAX = M2-M1+1+2*NKL0
      ALLOCATE( RLM(MLMAX) )
      DO J=1,MLMAX
         RLM(J) = M1-NKL0+J-1
      ENDDO

      INQUIRE(FILE='HARMONIC.IN', EXIST=EXISTS)
      IF (EXISTS) THEN
         HARMOFILE=ASSIGNFREEFILEUNIT ()
         
         OPEN(HARMOFILE, FILE='HARMONIC.IN',FORM='FORMATTED')
         READ (HARMOFILE,*) HARMONUM
         IF (ABS(HARMONUM).LT.100.0) THEN
            MLMAX=HARMONUM
            DEALLOCATE(RLM)
            ALLOCATE(RLM(MLMAX))
            DO J=1,MLMAX
               READ (HARMOFILE,*) HARMONUM 
               RLM(J)=HARMONUM
            ENDDO
            PRINT *,'RLM(1)=',RLM(1),'RLM(MLMAX)=',RLM(MLMAX)
         ENDIF
         CLOSE (HARMOFILE)
      ENDIF

      ELSE

      MLMAX = 1
      ALLOCATE( RLM(MLMAX) )
      RLM(1) = 0.0

      ENDIF

      RETURN
      END

C=======================================================================
C COMPUTE B,H, AND CREATE LAMBDA MESH                                  =    
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KLAMBDA

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      USE ToolBox
      USE MPIENV

      IMPLICIT NONE

      INTEGER JS,J,M,K
      REAL*8  RTMP,RLAM,MTMP,H1,H2,RTMP1,RTMP2,RTMP3

      INTEGER KCHECK,KSMOOTH,KM

      REAL*8,DIMENSION(:),ALLOCATABLE::BMC,BMS
      REAL*8,ALLOCATABLE:: TMPMESH(:)
      REAL*8 INTCOEF,GRID_STEP,GRID_X,GRID_X1,GRID_X2
      INTEGER TMPMESH_SIZE
      INTEGER IHE(100),NIHE
      
      KCHECK  = 1
      KSMOOTH = 1

      KM      = NKSMOOTHB
      ALLOCATE(BMC(KM+1),BMS(KM+1))

C     DEFINE THE EQUILIBRIUM CHI MESH
C     UNFORTUNATELY CHI IS DEFINED FROM 0 TO 2*PI IN CHEASE
C     THIS CAUSES SOME TROUBLE FOR CHI-INTEGRATION FOR 
C     TRAPPED PARTICLES
      RCHIH = 2*PI/NCHI
      RCHI(1) = 0.0
      DO J=2,NCHI+1
         RCHI(J) = RCHI(J-1) + RCHIH
      ENDDO

C     CALCULATE FIELD STRENGTH B AND H=B0/B 
      DO JS=2,NRP1
         DO J=1,NCHI
            RTMP=G22L(JS,J)*DPSIDS(JS)**2/RJA(JS,J)**2 +
     &           T(JS)**2/REQ(JS,J)**2
            BK(JS,J,1)=SQRT(RTMP)
            BPK(JS,J,1)=RJA(JS,J)**2*BK(JS,J,1)/G22L(JS,J)/
     &                  DPSIDS(JS)**2
         ENDDO
      ENDDO

      DO J=1,NCHI
         BK(1,J,1)=T(1)/REQ(1,J)
      ENDDO

      DO JS=1,NR
         DO J=1,NCHI
            RTMP=G22LM(JS,J)*DPSIDSM(JS)**2/RJAM(JS,J)**2 +
     &           TM(JS)**2/REQM(JS,J)**2
            BK(JS,J,2)=SQRT(RTMP)
            BPK(JS,J,2)=RJAM(JS,J)**2*BK(JS,J,2)/G22LM(JS,J)/
     &                  DPSIDSM(JS)**2
         ENDDO
      ENDDO

      DO J=1,NCHI
         BK(NRP1,J,2) = BK(NRP1,J,1)
      ENDDO

C     BPK SHOULD APPROACH INFINITY AT MAGNETIC AXIS
      DO J=1,NCHI
         BPK(1,J,1)=BPK(1,J,2)
      ENDDO
         
C     SMOOTH EQUILIBRIUM B FIELD FOR BAD EQUILIBRIUM
      IF (KSMOOTH.EQ.1) THEN     
      DO K=1,2
      DO JS=1,NRP1
         BMC(1) = 0.0
         DO J=1,NCHI
            BMC(1) = BMC(1) + BK(JS,J,K)
         ENDDO
         BMC(1) = BMC(1)*RCHIH/2.0/PI
         DO M=1,KM
            BMC(M+1) = 0.0
            BMS(M+1) = 0.0
            DO J=1,NCHI
               BMC(M+1)=BMC(M+1)+BK(JS,J,K)*COS(M*RCHI(J))
               BMS(M+1)=BMS(M+1)+BK(JS,J,K)*SIN(M*RCHI(J))
            ENDDO
            BMC(M+1)=BMC(M+1)*RCHIH/PI
            BMS(M+1)=BMS(M+1)*RCHIH/PI
         ENDDO
         DO J=1,NCHI
            BK(JS,J,K)=BMC(1) 
            DO M=1,KM
               BK(JS,J,K)=BK(JS,J,K)+BMC(M+1)*COS(M*RCHI(J))
               BK(JS,J,K)=BK(JS,J,K)+BMS(M+1)*SIN(M*RCHI(J))
            ENDDO 
         ENDDO
      ENDDO
      ENDDO
      ENDIF

      DO J=1,NCHI
         DO JS=1,NRP1
            HK(JS,J,1)=B0K/BK(JS,J,1)
         ENDDO
         DO JS=1,NR
            HK(JS,J,2)=B0K/BK(JS,J,2)
         ENDDO
      ENDDO

C     DEFINE THE WEIGHTING FUNCTION FOR KINETIC TERMS
C     MAINLY FOR ANNIHILATION OF CONTRIBUTION NEAR RATIONAL SURFACES
C     ACT ONLY ON NON_ADIABATIC CONTRIBUTIONS
      WFUN = 1.
      IF (DELRATS.GT.0..AND.DELRATS.LT.1.) THEN
      DO M=1,2
      DO JS=1,NRP1
         IF (M.EQ.1) H1 = CS (JS)
         IF (M.EQ.2) H1 = CSM(JS)
         DO K=1,NRATSURF
            J = IRATSURF(K)
            IF (ABS(H1-CS(J)).LE.DELRATS/Q(J)) THEN      
               WFUN(JS,M) = 0.
            ENDIF
         ENDDO
      ENDDO
      ENDDO
      ENDIF

C     FIND HMIN AND HMAX
      IF (1.EQ.0) THEN
      DO JS=1,NRP1
         HKMIN(JS,1)=HK(JS,1,1)
         HKMAX(JS,1)=HK(JS,1,1)
         DO J=2,NCHI
            IF (HKMIN(JS,1).GT.HK(JS,J,1)) HKMIN(JS,1)=HK(JS,J,1)
            IF (HKMAX(JS,1).LT.HK(JS,J,1)) HKMAX(JS,1)=HK(JS,J,1)
         ENDDO
      ENDDO

      DO JS=1,NR
         HKMIN(JS,2)=HK(JS,1,2)
         HKMAX(JS,2)=HK(JS,1,2)
         DO J=2,NCHI
            IF (HKMIN(JS,2).GT.HK(JS,J,2)) HKMIN(JS,2)=HK(JS,J,2)
            IF (HKMAX(JS,2).LT.HK(JS,J,2)) HKMAX(JS,2)=HK(JS,J,2)
         ENDDO
      ENDDO
      ENDIF

C     A NEW METHOD TO DEFINE HMIN AND HMAX
C     TO AVOID PROBLEM OF MULTIPLE TRAPPING
      IF (1.EQ.1) THEN
      HKMIN(1,1) = HK(1,1,1)
      HKMAX(1,1) = HK(1,1,1)
      DO M=1,2
      IHE = 0
      DO JS=3-M,NR-M+2
         K=0
         DO J=1,NCHI
            IF (J.GT.1) H1 = HK(JS,J-1,M)
            IF (J.EQ.1) H1 = HK(JS,NCHI,M)
            IF (J.LT.NCHI) H2 = HK(JS,J+1,M)
            IF (J.EQ.NCHI) H2 = HK(JS,1,M)
            RTMP = HK(JS,J,M)
            IF ((RTMP-H1)*(H2-RTMP).LE.0.) THEN
               K = K+1
               IHE(K) = J
            ENDIF
         ENDDO
         NIHE = K
         IF (NIHE.GT.100) THEN
            WRITE(*,*) 'JS,KGRID,NIHE=',JS,M,NIHE
            STOP 'KLAMBDA: NIHE'
         ELSEIF (NIHE.LT.2) THEN
            WRITE(*,*) 'JS,KGRID,NIHE=',JS,M,NIHE
            STOP 'KLAMBDA: NIHE'
         ELSEIF (NIHE.EQ.2) THEN
            H1 = HK(JS,IHE(1),M)  
            H2 = HK(JS,IHE(2),M)  
            IF (H1.LT.H2) THEN
               HKMIN(JS,M) = H1
               HKMAX(JS,M) = H2
            ELSE
               HKMIN(JS,M) = H2
               HKMAX(JS,M) = H1
            ENDIF
         ELSE
            WRITE(*,*) 'WARNING: MULTIPLE TRAPPING NEGLECTED'
            WRITE(*,*) 'AT JS,KGRID,NIHE=',JS,M,NIHE
            HKMIN(JS,M) = HK(JS,IHE(1),M)     
            DO K=2,NIHE
               H1 = HK(JS,IHE(K),M)
               IF (HKMIN(JS,M).GT.H1) HKMIN(JS,M) = H1
            ENDDO
            HKMAX(JS,M) = HK(JS,IHE(1),M)     
            DO K=2,NIHE
               H1 = HK(JS,IHE(K),M)
               IF (HKMAX(JS,M).GT.H1.AND.H1.GT.HKMIN(JS,M))
     &            HKMAX(JS,M) = H1
            ENDDO
         ENDIF
      ENDDO
      ENDDO
      ENDIF 

C     DEFINE ADAPTIVE LAMBDA-MESH
      MTMP = HKMAX(NRP1,1)
      DO JS=1,NRP1
        IF (HKMAX(JS,1).GT.MTMP) MTMP = HKMAX(JS,1)
      ENDDO
      DO JS=1,NR
        IF (HKMAX(JS,2).GT.MTMP) MTMP = HKMAX(JS,2)
      ENDDO
      RTMP = MTMP/DFLOAT(NLAMK-NLAMIN-1)
      
      DO JS=1,NRP1
         NLAMK1(JS,1) = INT(HKMIN(JS,1)/RTMP) + NLAMIN
         NLAMK0(JS,1) = INT((HKMAX(JS,1)-HKMIN(JS,1))/RTMP) + NLAMIN
      ENDDO

      DO JS=1,NR
         NLAMK1(JS,2) = INT(HKMIN(JS,2)/RTMP) + NLAMIN
         NLAMK0(JS,2) = INT((HKMAX(JS,2)-HKMIN(JS,2))/RTMP) + NLAMIN
      ENDDO

C     DEFINE PITCH ANGLE ARRAY FOR PASSING PARTICLES
C     USING LOG-MESH 
      RLAM = 4.0

      DO JS=1,NRP1
         RTMP=(-RLAM-LOG10(HKMIN(JS,1)))/(NLAMK1(JS,1)-2)
         DO J=1,NLAMK1(JS,1)-1
            LAMK1(JS,J,1)=HKMIN(JS,1) - 
     &                    10.**(LOG10(HKMIN(JS,1))+(J-1)*RTMP)
         ENDDO
         LAMK1(JS,NLAMK1(JS,1),1)=HKMIN(JS,1)
      ENDDO

      DO JS=1,NR
         RTMP=(-RLAM-LOG10(HKMIN(JS,2)))/(NLAMK1(JS,2)-2)
         DO J=1,NLAMK1(JS,2)-1
            LAMK1(JS,J,2)=HKMIN(JS,2) - 
     &                    10.**(LOG10(HKMIN(JS,2))+(J-1)*RTMP)
         ENDDO
         LAMK1(JS,NLAMK1(JS,2),2)=HKMIN(JS,2)
      ENDDO

C     DEFINE PITCH ANGLE ARRAY FOR TRAPPED PARTICLES
C     USING LOG-MESH 
      DO JS=1,NRP1
         RTMP = (LOG10(HKMAX(JS,1)-HKMIN(JS,1))+RLAM)/(NLAMK0(JS,1)-2)
         DO J=2,NLAMK0(JS,1)
            LAMK0(JS,J,1)=HKMIN(JS,1) + 10.**(-RLAM+(J-2)*RTMP)
         ENDDO
         LAMK0(JS,1,1)=HKMIN(JS,1)
      ENDDO

      DO JS=1,NR
         RTMP = (LOG10(HKMAX(JS,2)-HKMIN(JS,2))+RLAM)/(NLAMK0(JS,2)-2)
         DO J=2,NLAMK0(JS,2)
            LAMK0(JS,J,2)=HKMIN(JS,2) + 10.**(-RLAM+(J-2)*RTMP)
         ENDDO
         LAMK0(JS,1,2)=HKMIN(JS,2)
      ENDDO

C     COMPUTE ADDITIONAL RADIAL DERIVATIVES FOR ANISOTROPIC HOT IONS
      FFF = HKMAX
      CALL DFFFDPSI(0)
      ZDHMAXDPSI = DFFF

      FFF = HKMIN/HKMAX
      CALL DFFFDPSI(0)
      ZDHMINMAXDPSI = DFFF

      FFF = 1.0/HKMIN
      CALL DFFFDPSI(0)
      HMIN3DPSI = DFFF

      FFF = RTYPE4
      CALL DFFFDPSI(0)
      RTYPE4DPSI = DFFF

      FFF = STYPE4
      CALL DFFFDPSI(0)
      STYPE4DPSI = DFFF

C     COMPUTE SECOND ORDER RADIAL DERIVATIVES FOR ANISOTROPIC HOT IONS
C     WITH FOW CORRECTION
      FFF = ZDHMAXDPSI/HKMAX**2
      CALL DFFFDPSI(0)
      ZDHMAXDPSI2 = DFFF

      FFF = ZDHMINMAXDPSI
      CALL DFFFDPSI(0)
      ZDHMINMAXDPSI2 = DFFF

      FFF = HMIN3DPSI
      CALL DFFFDPSI(0)
      HMIN3DPSI2 = DFFF

      FFF = RTYPE4DPSI
      CALL DFFFDPSI(0)
      RTYPE4DPSI2 = DFFF

      FFF = STYPE4DPSI
      CALL DFFFDPSI(0)
      STYPE4DPSI2 = DFFF

C     COMPUTE HTYPE4C FACTOR FOR IF0TYPE=4
C     HHTYPE4C IS THE NAMELIST VARIABLE
      RTMP1   = MINVAL(HKMIN(:,1))
      RTMP2   = MAXVAL(HKMAX(:,1))
      HTYPE4C = RTMP1 + HHTYPE4C*(RTMP2-RTMP1)

      IF (KCHECK.EQ.1.AND.(ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT)) 
     &   WRITE(*,*) 'RTMP1,RTMP2,HTYPE4C=',RTMP1,RTMP2,HTYPE4C

C     CREATE MESH FOR ENERGY INTEGRAL INCLUDING ENERGY DEPENDENT COLLISIONALITY
      IF (INUTYPE .NE. 0 ) THEN
        TMPMESH_SIZE = intMeshPot/2 + 1
        ALLOCATE(TMPMESH( TMPMESH_SIZE ))
        TMPMESH = 0.0
        CALL createEnergyMesh(TMPMESH(1:1500),
     &                              0.0D0,    1.0D0,    1)
        CALL createEnergyMesh(TMPMESH(1500:4800),
     &                              1.0D0,    25D0,     1)
        CALL createEnergyMesh(TMPMESH(4800:TMPMESH_SIZE),
     &                              25D0,     200D0,     1)        

      INTCOEF = 1.0/sqrt(3.0)
      
      IF (INUTYPE .GE. 1) THEN
C        CALL gaulag(intMeshX,intMeshWeight,intMeshPot,4.0_8)
         DO J=1,TMPMESH_SIZE - 1      
          GRID_STEP = TMPMESH(J+1) - TMPMESH(J)
          GRID_X = 0.5 * ( TMPMESH(J) + TMPMESH(J+1) )
          GRID_X1 = GRID_X -  INTCOEF * GRID_STEP * 0.5
          GRID_X2 = GRID_X +  INTCOEF * GRID_STEP * 0.5
          
          intMeshX(2*J-1) = GRID_X1 
          intMeshWeight(2*J-1)=EXP(-GRID_X1)*GRID_X1**4.0*GRID_STEP*0.5
          
          intMeshX(2*J)   = GRID_X2
          intMeshWeight(2*J) =EXP(-GRID_X2)*GRID_X2**4.0*GRID_STEP*0.5
         ENDDO

      ELSEIF(INUTYPE.EQ.-1) THEN
C        TEST CASE FOR RECOVERING THE RESULT WITH COLLISIONALITY INDEPENDENT OF ENERGY      
C        CALL gaulag(intMeshX,intMeshWeight,intMeshPot,2.5_8)
         DO J=1,TMPMESH_SIZE - 1     
          GRID_STEP = TMPMESH(J+1) - TMPMESH(J)
          GRID_X = 0.5 * ( TMPMESH(J) + TMPMESH(J+1) )
          GRID_X1 = GRID_X -  INTCOEF * GRID_STEP * 0.5
          GRID_X2 = GRID_X +  INTCOEF * GRID_STEP * 0.5
          
          intMeshX(2*J-1) = GRID_X1 
          intMeshWeight(2*J-1)=EXP(-GRID_X1)*GRID_X1**2.5*GRID_STEP*0.5
          
          intMeshX(2*J)   = GRID_X2
          intMeshWeight(2*J) =EXP(-GRID_X2)*GRID_X2**2.5*GRID_STEP*0.5
         ENDDO

      ENDIF
      
      DEALLOCATE(TMPMESH,BMC,BMS)
      
      ENDIF
      
      
      IF (ISMPIRUN.EQ.0) THEN  
      IF (KCHECK.EQ.1) THEN
C        CHECK HMIN AND HMAX
         WRITE(*,*) 'CHECK KLAMBDA: CS HKMIN HKMAX NLAMK1 NLAMK0'
         DO JS=1,NRP1
            WRITE(*,110) CS(JS),HKMIN(JS,1),HKMAX(JS,1),
     &                   NLAMK1(JS,1),NLAMK0(JS,1)
         ENDDO
 110     FORMAT(3(E13.4,2X),2(I3,2X))
      ENDIF

      IF (KCHECK.EQ.2) THEN
C        CHECK PITCH ANGLE ARRAY
         WRITE(*,*) 'CHECK KLAMBDA: LAMK1'
         WRITE(*,120) (LAMK1(JS0,J,1),J=1,NLAMK1(JS0,1))
         WRITE(*,*) 'CHECK KLAMBDA: LAMK0'
         WRITE(*,120) (LAMK0(JS0,J,1),J=1,NLAMK0(JS0,1))
 120     FORMAT(E13.4)
      ENDIF
      
      IF (KCHECK.EQ.3) THEN
C        OUTPUT H(CHI) FOR SURFACE JS=JS0
         WRITE(*,*) 'CHECK KLAMBDA: RCHI HK'
         DO J=1,NCHI
            WRITE(*,130) RCHI(J),HK(JS0,J,2)
         ENDDO
 130     FORMAT(2(E17.8,2X))
      ENDIF
      
      ENDIF
      
      RETURN
      END
C=======================================================================
C COMPUTE DENSITY AND PRESSURE FRACTIONS FOR ALPHA PARTICLES WITH 
C THE IF0TYPE=2 MODEL
C YQL, 06-2013                                                        
C=======================================================================
      SUBROUTINE KNPFRACF02

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      IMPLICIT NONE

      INTEGER JS,K,KP
      REAL*8  H2,DRHO,RTMP,RTMPI,RTMPE,AF1,AF3,SMI2ME,AMA2MI

      INTEGER KCHECK
      KCHECK=0

C     GO THROUGH ALL HOT ION SPECIES
      DO KP=3,NSPECIES

C     COMPUTE FAST ION DENSITY AND PRESSURE FRACTION FOR IF0TYPE=2
C     USING (APPROXIMATE) ANALYTIC ALPHA-PARTICLE MODEL
      IF (ISPECIES_F0(KP).EQ.2) THEN  
         SMI2ME = SQRT(ESPECIES_M(1)/ESPECIES_M(2))
         AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
         DRHO   = (3.*SQRT(PI)/4.*SMI2ME)**(2./3.)*AMA2MI
         DO K=1,2
         DO JS=1,NR
            IF (K.EQ.1) RTMPI = TEMPI(JS)
            IF (K.EQ.2) RTMPI = TEMPIM(JS)
            IF (K.EQ.1) RTMPE = TEMPE(JS)
            IF (K.EQ.2) RTMPE = TEMPEM(JS)
            AF1  = RTMPE*ZTE0/3.5E+6
            RTMP = DRHO*AF1
            H2   = LOG(1.+RTMP**(-1.5))
            AF3  = (1.-RTMP/9.*(SQRT(3.)*PI-6.*SQRT(3.)
     &             *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))
     &             -6.*LOG(1.+SQRT(RTMP))
     &             +3.*LOG(1.-SQRT(RTMP)+RTMP)))/H2
            RTMP = 2.925E-7*(RTMPE*ZTE0/1.E+3)**1.5
     &             *(RTMPI*ZTI0/1.E+3)**2
            H2   = RTMP*AF3/AF1*(1.-ALPHAP)
            ESPECIES_DENF(JS,K,KP) = RTMP
            ESPECIES_PREF(JS,K,KP) = H2
         ENDDO
         ENDDO
         JS = NRP1
         ESPECIES_DENF(JS,1,KP)=2.*ESPECIES_DENF(JS-1,2,KP)-
     &                             ESPECIES_DENF(JS-1,1,KP)
         ESPECIES_PREF(JS,1,KP)=2.*ESPECIES_PREF(JS-1,2,KP)-
     &                             ESPECIES_PREF(JS-1,1,KP)

         IF (KCHECK.EQ.1) THEN
            WRITE(*,*) 'KNPFRACF02: CS,DAK,PAK'
            DO JS=1,NR
               WRITE(*,110) CS(JS),ESPECIES_DENF(JS,1,KP),
     &                      ESPECIES_PREF(JS,1,KP)
            ENDDO
 110        FORMAT(3(E13.5,1X))
        ENDIF
      ENDIF

      ENDDO

      RETURN
      END

C=======================================================================
C GET (SURFACE AVERAGED) EQUILIBRIUM RADIAL PROFILES OF DENSITY,
C PRESSURE, TEMPERATURE, AS WELL AS DP/DPSI FOR EACH PARTICLE SPECIES
C NOTE THAT THE PRESSURE FRACTIONS ARE DEFINED W.R.T. TOTAL THERMAL PRESSURE
C AND THE DENSITY FRACTIONS ARE DEFINED W.R.T. THERMAL ELECTRON DENSITY
C YQL, 06-2013
C=======================================================================
      SUBROUTINE KEQPROF

      USE DIMENSIM
      USE GLOBALM
      IMPLICIT NONE

      INTEGER K,JS

      INTEGER KCHECK
      REAL*8 TI0_TMP,TE0_TMP
      KCHECK=0

C     DENSITY FRACTION FOR THERMAL ELECTRONS IS ALWAYS 1
C     SINCE WE ASSUME EQUILIBRIUM DENSITY RHO=ELECTRON DENSITY
      ESPECIES_DENF(:,:,2) = 1.
      ESPECIES_DEN (:,1,2) = ESPECIES_DENF(:,1,2)*RHO
      ESPECIES_DEN (:,2,2) = ESPECIES_DENF(:,2,2)*RHOM

C     CALCULATE DENSITY FOR EP SPECIES
      DO K=3,NSPECIES
         ESPECIES_DEN(:,1,K) = ESPECIES_DENF(:,1,K)*RHO
         ESPECIES_DEN(:,2,K) = ESPECIES_DENF(:,2,K)*RHOM
      ENDDO

C     CALCULATE THERMAL ION DENSITY (AND FRACTION) FROM THE NUETRALITY 
C     CONDITION: SUM(Z*N)=0
      ESPECIES_DEN(:,:,1) = 0.
      DO K=2,NSPECIES
         ESPECIES_DEN(:,:,1) = ESPECIES_DEN(:,:,1) - 
     &                         ESPECIES_DEN(:,:,K)*ESPECIES_Z(K)
      ENDDO
      ESPECIES_DEN (:,:,1) = ESPECIES_DEN(:,:,1)/ESPECIES_Z(1)
      ESPECIES_DENF(:,1,1) = ESPECIES_DEN(:,1,1)/RHO
      ESPECIES_DENF(:,2,1) = ESPECIES_DEN(:,2,1)/RHOM

      ESPECIES_DEN(NRP1,:,:) = ESPECIES_DEN(NR,:,:)

C     PRESSURE FRACTIONS FOR THERMAL PARTICLES  
      ESPECIES_PREF(:,:,1) = ALPHAP
      ESPECIES_PREF(:,:,2) = 1.-ALPHAP

C     CALCULATE THERMAL PRESSURE
C     AND STORE TEMPORARILY IN ESPECIES_TEM(:,:,1)
      ESPECIES_TEM(:,:,1) = 0.
      DO K=1,NSPECIES
         ESPECIES_TEM(:,:,1) = ESPECIES_TEM(:,:,1) + 
     &                         ESPECIES_PREF(:,:,K)
      ENDDO
      ESPECIES_TEM(:,1,1) = PEQ /ESPECIES_TEM(:,1,1)
      ESPECIES_TEM(:,2,1) = PEQM/ESPECIES_TEM(:,2,1)

C     CALCULATE PRESSURE FOR EACH SPECIES
      DO K=1,NSPECIES
         ESPECIES_PRE(:,:,K) = ESPECIES_PREF(:,:,K)*
     &                         ESPECIES_TEM(:,:,1)
      ENDDO

C     RECOVER OLD VERSION (THIS PART IS ONLY FOR TEST)
      IF (1.EQ.0) THEN
      ESPECIES_DEN(:,:,1) = ESPECIES_DEN(:,:,2)
      ENDIF

C     RESCALE THERMAL PRESSURE WITH THE EXPERIMENTAL TEMPERATURE PROFILES
      IF (NPROFIE.EQ.4.AND.1.EQ.1) THEN
          ESPECIES_PRE(:,:,2)=0.
          DO K=3,NSPECIES
             ESPECIES_PRE(:,:,2) = ESPECIES_PRE(:,:,2)
     &                           + ESPECIES_PRE(:,:,K)
          ENDDO
          ESPECIES_PRE(:,1,1) = PEQ-ESPECIES_PRE(:,1,2)
          ESPECIES_PRE(:,2,1) = PEQM-ESPECIES_PRE(:,2,2)
          TI0_TMP=ALPHAP*ESPECIES_PRE(1,1,1)/ESPECIES_DEN(1,1,1)
          TE0_TMP=(1-ALPHAP)*ESPECIES_PRE(1,1,1)/ESPECIES_DEN(1,1,2)
          ESPECIES_TEM(:,1,1) = TI0_TMP*TEMPI  
          ESPECIES_TEM(:,2,1) = TI0_TMP*TEMPIM  
          ESPECIES_TEM(:,1,2) = TE0_TMP*TEMPE  
          ESPECIES_TEM(:,2,2) = TE0_TMP*TEMPEM
          ESPECIES_TEM(NRP1,:,1:2) = ESPECIES_TEM(NR,:,1:2)

          ESPECIES_PREF(:,:,1:2) = ESPECIES_TEM(:,:,1:2)
     &                           * ESPECIES_DEN(:,:,1:2)
          ESPECIES_PREF(:,:,1) = ESPECIES_PREF(:,:,1)
     &                         / ( ESPECIES_PREF(:,:,1)
     &                         + ESPECIES_PREF(:,:,2) )
          ESPECIES_PREF(:,:,2) = 1.0 
     &                         - ESPECIES_PREF(:,:,1)
          ESPECIES_PRE(:,:,2) = ESPECIES_PRE(:,:,1)
     &                        * ESPECIES_PREF(:,:,2)
          ESPECIES_PRE(:,:,1) = ESPECIES_PRE(:,:,1)
     &                        * ESPECIES_PREF(:,:,1)

      ENDIF

C     CALCULATE TEMPERATURE FOR EACH SPECIES
      ESPECIES_TEM = ESPECIES_PRE/ESPECIES_DEN
      ESPECIES_TEM(NRP1,:,:) = ESPECIES_TEM(NR,:,:)

C     COMPUTE DP/DPSI FOR EACH PARTICLE SPECIES
      DO K=1,NSPECIES
         FFF = ESPECIES_PRE(:,:,K)
         CALL DFFFDPSI(0)
         ESPECIES_PREP(:,:,K) = DFFF
      ENDDO

      RETURN
      END

C=======================================================================
C COMPUTE QUANTITIES ALPHAA1,2,3 ETC. 
C YQL, 06-2013                                                        
C=======================================================================
      SUBROUTINE KALPHAA

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE MPIENV
      USE KINETICM
      USE ANISOTROPICM
      USE REORBITM
      IMPLICIT NONE

      INTEGER JS,K,J,KP,KCHECK,N
      REAL*8  H1,H2,DRHO,RTMP,RTMPI,RTMPE,AX1,AX2,AX3,AF1,AF2,AF3,
     &        SMI2ME,AMA2MI

      KCHECK = 2

      ZC1    = 0.

C     GO THROUGH ALL PARTICLE SPECIES
      DO KP=1,NSPECIES

C     FOR MAXWELLIAN DISTRIBUTION: EPSILON_ALPHA=ESPECIES_TEM(:,:,KP)
C     ==> ALPHAA1,2,3
      IF (ISPECIES_F0(KP).EQ.0) THEN
         EPSALPHA(:,:,KP) = ESPECIES_TEM(:,:,KP) 
         ALPHAA1(:,:,KP)  = ESPECIES_TEM(:,:,2)/EPSALPHA(:,:,KP)
         ALPHAA2(:,:,KP)  = ESPECIES_TEM(:,:,1)/EPSALPHA(:,:,KP)
         ALPHAA3(:,:,KP)  = 1.
         EPSLONCA(:,:,KP) = 1.
         AAK(:,:,KP)      = 1.
         CPSI(:,:,KP)     = ESPECIES_PRE(:,:,KP)
      ENDIF
         
C     FOR ISOTROPIC PITCH ANGLE, SLOWING DOWN ENERGY HOT ION MODEL
C     FIND ALPHAA1 AS ROOT OF A NON-LINEAR EQUATION
      IF (ISPECIES_F0(KP).EQ.1) THEN
         SMI2ME = SQRT(ESPECIES_M(1)/ESPECIES_M(2))
         AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
         DRHO = (3.*SQRT(PI)/4.*SMI2ME)**(2./3.)*AMA2MI
         H1   = 3./8./SQRT(2.)/PI
         DO K=1,2
         DO JS=1,NR
            IF (K.NE.1.OR.JS.NE.1) THEN
            H2 = ESPECIES_PREF(JS,K,KP)/ESPECIES_DENF(JS,K,KP)
            AX1 = 1.E-8 
            AX2 = 1.0E+1
            AX3 = (AX1+AX2)*0.5
            RTMP= DRHO*AX1
            AF1 = (1.-ALPHAP)*(1.-RTMP/9.*(SQRT(3.)*PI-6.*SQRT(3.)
     &            *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))-6.*
     &            LOG(1.+SQRT(RTMP))+3.*LOG(1.-SQRT(RTMP)+RTMP)))
     &            /LOG(1.+RTMP**(-1.5))-AX1*H2
            RTMP= DRHO*AX2
            AF2 = (1.-ALPHAP)*(1.-RTMP/9.*(SQRT(3.)*PI-6.*SQRT(3.)
     &            *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))-6.*
     &            LOG(1.+SQRT(RTMP))+3.*LOG(1.-SQRT(RTMP)+RTMP)))
     &            /LOG(1.+RTMP**(-1.5))-AX2*H2
            RTMP= DRHO*AX3
            AF3 = (1.-ALPHAP)*(1.-RTMP/9.*(SQRT(3.)*PI-6.*SQRT(3.)
     &            *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))-6.*
     &            LOG(1.+SQRT(RTMP))+3.*LOG(1.-SQRT(RTMP)+RTMP)))
     &            /LOG(1.+RTMP**(-1.5))-AX3*H2
     
            LOOP1: DO J=1,1000
               IF (AF1*AF3.LE.0.) THEN
                  AX2 = AX3
                  AF2 = AF3
               ELSEIF (AF2*AF3.LE.0.) THEN
                  AX1 = AX3
                  AF1 = AF3
               ELSE
                  WRITE(*,*) 'KALPHAA:K,JS,J,H2,AX1,AX2,AX3,AF1,AF2,AF3=
     &                       ',K,JS,J,H2,AX1,AX2,AX3,AF1,AF2,AF3
                  IF (ISMPIRUN.EQ.0.AND.KCHECK.EQ.2) THEN
                     DO N=1,101
                        AX2 = 1.0E+1*(N-1.)/100.+1.E-8
                        RTMP= DRHO*AX2
                        AF2 = (1.-ALPHAP)*(1.-RTMP/9.*(SQRT(3.)*PI-
     &                         6.*SQRT(3.)
     &                        *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))-6.*
     &                         LOG(1.+SQRT(RTMP))+3.*LOG(1.-SQRT(RTMP)+
     &                         RTMP)))
     &                        /LOG(1.+RTMP**(-1.5))/AX2-H2
                        WRITE(*,*) AX2,AF2
                     ENDDO
                  ENDIF
                  STOP 'KALPHAA:ALPHAA1'
               ENDIF
               AX3 = (AX1+AX2)*0.5
               RTMP= DRHO*AX3
               AF3 = (1.-ALPHAP)*(1.-RTMP/9.*(SQRT(3.)*PI-6.*SQRT(3.)
     &               *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))-6.*
     &               LOG(1.+SQRT(RTMP))+3.*LOG(1.-SQRT(RTMP)+RTMP)))
     &               /LOG(1.+RTMP**(-1.5))-AX3*H2
               IF (ABS((AX1-AX2)/AX3).LE.1.E-5) EXIT LOOP1
            ENDDO LOOP1

            ALPHAA1(JS,K,KP) = AX3
            RTMP             = DRHO*ALPHAA1(JS,K,KP)
            H2               = LOG(1.+RTMP**(-1.5))
            ALPHAA3(JS,K,KP) = (1.-RTMP/9.*(SQRT(3.)*PI-6.*SQRT(3.)
     &                         *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))
     &                         -6.*LOG(1.+SQRT(RTMP))
     &                         +3.*LOG(1.-SQRT(RTMP)+RTMP)))/H2
            EPSLONCA(JS,K,KP)= RTMP
            ALPHAA2(JS,K,KP) = ALPHAA1(JS,K,KP)*ALPHAP/(1.-ALPHAP)     
            AAK(JS,K,KP)     = H1/H2
            ENDIF
         ENDDO
         ENDDO
         J = 1
         ALPHAA1(J,1,KP)  = 2.*ALPHAA1(J,2,KP)-ALPHAA1(J+1,1,KP)
         ALPHAA2(J,1,KP)  = 2.*ALPHAA2(J,2,KP)-ALPHAA2(J+1,1,KP)
         ALPHAA3(J,1,KP)  = 2.*ALPHAA3(J,2,KP)-ALPHAA3(J+1,1,KP)
         EPSLONCA(J,1,KP) = 2.*EPSLONCA(J,2,KP)-EPSLONCA(J+1,1,KP)
         AAK(J,1,KP)      = 2.*AAK(J,2,KP)-AAK(J+1,1,KP)

         J = NRP1
         ALPHAA1(J,1,KP)  = 2.*ALPHAA1(J-1,2,KP)-ALPHAA1(J-1,1,KP)
         ALPHAA2(J,1,KP)  = 2.*ALPHAA2(J-1,2,KP)-ALPHAA2(J-1,1,KP)
         ALPHAA3(J,1,KP)  = 2.*ALPHAA3(J-1,2,KP)-ALPHAA3(J-1,1,KP)
         EPSLONCA(J,1,KP) = 2.*EPSLONCA(J-1,2,KP)-EPSLONCA(J-1,1,KP)
         AAK(J,1,KP)      = 2.*AAK(J-1,2,KP)-AAK(J-1,1,KP)

         CPSI(:,:,KP) = (2.*PI)**1.5*ESPECIES_DEN(:,:,KP)*AAK(:,:,KP)
     &                  *ESPECIES_TEM(:,:,2)/ALPHAA1(:,:,KP)
         EPSALPHA(:,:,KP) = ESPECIES_TEM(:,:,2)/ALPHAA1(:,:,KP)

         IF (KEPSALPHA.EQ.1) THEN
            EPSALPHA(:,1,KP) = SUM(EPSALPHA(1:NRP1,1,KP))/DFLOAT(NRP1)
            EPSALPHA(:,2,KP) = SUM(EPSALPHA(1:NR,2,KP))/DFLOAT(NR)
         ENDIF
         IF (ABS(PSPECIES_AT(KP)).GT.0.) CALL KJP_ISOEXTRA(KP)    
      ENDIF

C     COMPUTE ALPHAA123 
C     USING (APPROXIMATE) ANALYTIC ALPHA-PARTICLE MODEL
      IF (ISPECIES_F0(KP).EQ.2) THEN  
         SMI2ME = SQRT(ESPECIES_M(1)/ESPECIES_M(2))
         AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
         DRHO = (3.*SQRT(PI)/4.*SMI2ME)**(2./3.)*AMA2MI
         H1   = 3./8./SQRT(2.)/PI
         DO K=1,2
         DO JS=1,NR
            IF (K.EQ.1) RTMPI = TEMPI(JS)
            IF (K.EQ.2) RTMPI = TEMPIM(JS)
            IF (K.EQ.1) RTMPE = TEMPE(JS)
            IF (K.EQ.2) RTMPE = TEMPEM(JS)
            ALPHAA1(JS,K,KP) = RTMPE*ZTE0/3.5E+6
            ALPHAA2(JS,K,KP) = ALPHAA1(JS,K,KP)*ALPHAP/(1.-ALPHAP)  
            RTMP             = DRHO*ALPHAA1(JS,K,KP)
            EPSLONCA(JS,K,KP)= RTMP
            H2               = LOG(1.+RTMP**(-1.5))
            AAK(JS,K,KP)     = H1/H2
            ALPHAA3(JS,K,KP) = (1.-RTMP/9.*(SQRT(3.)*PI-6.*SQRT(3.)
     &                         *ATAN((-2.+SQRT(RTMP))/SQRT(3.*RTMP))
     &                         -6.*LOG(1.+SQRT(RTMP))
     &                         +3.*LOG(1.-SQRT(RTMP)+RTMP)))/H2
         ENDDO
         ENDDO
         J = NRP1
         ALPHAA1(J,1,KP)  = 2.*ALPHAA1(J-1,2,KP)-ALPHAA1(J-1,1,KP)
         ALPHAA2(J,1,KP)  = 2.*ALPHAA2(J-1,2,KP)-ALPHAA2(J-1,1,KP)
         ALPHAA3(J,1,KP)  = 2.*ALPHAA3(J-1,2,KP)-ALPHAA3(J-1,1,KP)
         EPSLONCA(J,1,KP) = 2.*EPSLONCA(J-1,2,KP)-EPSLONCA(J-1,1,KP)
         AAK(J,1,KP)      = 2.*AAK(J-1,2,KP)-AAK(J-1,1,KP)
         CPSI(:,:,KP) = (2.*PI)**1.5*ESPECIES_DEN(:,:,KP)*AAK(:,:,KP)
     &                  *ESPECIES_TEM(:,:,2)/ALPHAA1(:,:,KP)
         EPSALPHA(:,:,KP) = ESPECIES_TEM(:,:,2)/ALPHAA1(:,:,KP)
         IF (ABS(PSPECIES_AT(KP)).GT.0.) CALL KJP_ISOEXTRA(KP)    
      ENDIF

C     FOR ANISOTROPIC MODEL WITH IF0TYPE=3,7
C     FIND ALPHAA1 AS ROOT OF A NONLINEAR EQUATION
      IF (ISPECIES_F0(KP).EQ.3.OR.ISPECIES_F0(KP).EQ.7) THEN
         IF (NNSCAN.EQ.1) THEN
         OPEN(132,FILE='DISTRIBUTION.OUT',FORM='FORMATTED',
     &            STATUS='OLD')
         DO K=1,2
         DO JS=1,NRP1
         READ(132,*)     ZZI1(JS,K,KP),ZZI3(JS,K,KP),CPSI(JS,K,KP),
     &                   ALPHAA1(JS,K,KP),ALPHAA3(JS,K,KP),
     &                   EPSLONCA(JS,K,KP)
         ENDDO
         ENDDO
         DO JS=1,NEPK
         READ(132,*)    EPK(JS)
         ENDDO
         DO JS=1,2*NEPK-2
         READ(132,*)    ZEPKO(JS),ZEPKN(JS)
         ENDDO
         CLOSE(132)

         ELSE

         CALL KALPHAA_TYPE3(KP)
         CALL KCPSI_TYPE3  (KP)

         OPEN(132,FILE='DISTRIBUTION.OUT',FORM='FORMATTED',
     &            STATUS='REPLACE')
         DO K=1,2
         DO JS=1,NRP1
         WRITE(132,*)    ZZI1(JS,K,KP),ZZI3(JS,K,KP),CPSI(JS,K,KP),
     &                   ALPHAA1(JS,K,KP),ALPHAA3(JS,K,KP),
     &                   EPSLONCA(JS,K,KP)
         ENDDO
         ENDDO

         DO JS=1,NEPK
         WRITE(132,*)  EPK(JS)
         ENDDO
         DO JS=1,2*NEPK-2
         WRITE(132,*)  ZEPKO(JS),ZEPKN(JS)
         ENDDO
         CLOSE(132)
         ENDIF

         ALPHAA2(:,:,KP) = ALPHAA1(:,:,KP)*ESPECIES_TEM(:,:,1)
     &                     /ESPECIES_TEM(:,:,2)
         AAK(:,:,KP)     = 1.
         EPSALPHA(:,:,KP)= ESPECIES_TEM(:,:,2)/ALPHAA1(:,:,KP)

         IF (KEPSALPHA.EQ.1) THEN
            EPSALPHA(:,1,KP) = SUM(EPSALPHA(1:NRP1,1,KP))/DFLOAT(NRP1)
            EPSALPHA(:,2,KP) = SUM(EPSALPHA(1:NR,2,KP))/DFLOAT(NR)
         ENDIF
      ENDIF

C     FOR ICRH ANISOTROPIC MODEL
      IF (ISPECIES_F0(KP).EQ.4) THEN
         CALL KCPSI_TYPE4(KP)
         ALPHAA1(:,:,KP)  = ESPECIES_TEM(:,:,2)/EPSALPHA(:,:,KP)
         ALPHAA2(:,:,KP)  = ESPECIES_TEM(:,:,1)/EPSALPHA(:,:,KP)
         ALPHAA3(:,:,KP)  = ESPECIES_TEM(:,:,KP)/EPSALPHA(:,:,KP)
         EPSLONCA(:,:,KP) = 1.
         AAK(:,:,KP)      = 1.
      ENDIF
         
C     FOR RE MODEL-1 WITH IF0TYPE=5
C     FIND ALPHAA1 AS ROOT OF A NONLINEAR EQUATION
C     NOTE THAT USING THE SAME ROOT-FINDING PROCEDURE AS FOR TYPE3 
      IF (ISPECIES_F0(KP).EQ.5) THEN
         CALL KALPHAA_TYPE3(KP)

         ALPHAA2(:,:,KP) = ALPHAA1(:,:,KP)*ESPECIES_TEM(:,:,1)
     &                     /ESPECIES_TEM(:,:,2)
         AAK(:,:,KP)     = 1.
         EPSALPHA(:,:,KP)= ESPECIES_TEM(:,:,2)/ALPHAA1(:,:,KP)

         CALL KCPSI_TYPE5(KP)
      ENDIF

C     FOR RE MODEL-2 WITH IF0TYPE=6
C     ANALYTIC DISTRIBUTION MODEL WITH FIXED EPSALPHA  
      IF (ISPECIES_F0(KP).EQ.6) THEN  
         AMA2MI = ESPECIES_M(KP)/ESPECIES_M(1)
         RTMP   = AMA2MI*C_VA**2/2.0
         EPSALPHA(:,:,KP) = RTMP
         ALPHAA1(:,:,KP)  = ESPECIES_TEM(:,:,2)/RTMP
         ALPHAA2(:,:,KP)  = ESPECIES_TEM(:,:,1)/RTMP
         ALPHAA3(:,:,KP)  = 1.0 !=UPPER BOUND FOR RE
         ESPECIES_TEM(:,:,KP) = ALPHAA3(:,:,KP)*RTMP
         AAK(:,:,KP)      = 1.0
         EPSLONCA(:,:,KP) = 1.0

         CALL KCPSI_TYPE6(KP)
      ENDIF

C     COMPUTE D LN(CPSI)/D PSI
      FFF = CPSI(:,:,KP)
      CALL DFFFDPSI(1)
      DCDPSIL(:,:,KP) = DFFF
      
C     COMPUTE D^2 LN(CPSI)/D PSI^2
      FFF = DCDPSIL(:,:,KP)
      CALL DFFFDPSI(0)
      DCDPSIL2(:,:,KP) = DFFF

      ENDDO

      RETURN
      END

C=======================================================================
C SPECIFY COLLISIONALITY COEFFICIENTS FOR THERMAL PARTICLES
C AND COMPUTE OMEGAS* FREQUENCIES FOR THERMAL IONS, ELECTRONS
C AND FAST IONS                                              
C YQL, 08-2007                                               
C=======================================================================
      SUBROUTINE KDIAMAG

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE MPIENV
      USE ANISOTROPICM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'

      INTEGER JS,K,J,KP
      REAL*8  DRHO,RTMP,RTMPI,RTMPE,AX1,AX2,AX3,AF1,AF2,AF3
      REAL*8  H1,H2,H3,H4,H5,H6

      INTEGER KCHECK
      KCHECK=1

C     SPECIFY EFFECTIVE COLLISIONALITY PROFILES
C     SO FAR SPECIFIED FOR THERMAL PARTICLES ONLY (=0 FOR HOT IONS)
C     ONE CAN USE EITHER ANALYTICAL FORMULA OR GET FROM CHEASE OUTPUT
C     READ AMPLITUDE FROM NAMELIST FILE
      NUEFF = 0.
      DO JS=1,NRP1
         IF (NUMODEL .EQ. 0 ) THEN
             NUEFF(JS,1,1) = GNUI(JS)
             NUEFF(JS,2,1) = GNUIM(JS)
             NUEFF(JS,1,2) = GNUE(JS)
             NUEFF(JS,2,2) = GNUEM(JS)
         ELSE
            NUEFF(JS,1,1) = GNUI(JS)*ASPCT/(CS(JS)*2.0)
            NUEFF(JS,2,1) = GNUIM(JS)*ASPCT/(CSM(JS)*2.0)
            NUEFF(JS,1,2) = GNUE(JS)*ASPCT/(CS(JS)*2.0)
            NUEFF(JS,2,2) = GNUEM(JS)*ASPCT/(CSM(JS)*2.0)
         ENDIF
      ENDDO
      NUEFF(1,1,:) = NUEFF(1,2,:)
      
      DO KP=1,NSPECIES

C     COMPUTE 2ND ORDER DERIVATIVE D^2 (EPSILONCA**1.5)/D PSI^2
C     FIRST COMPUTE THE FIRST ORDER DERIVATIVE, USING 
C     DEPSALPHADPSI AS TEMPORARY VARIABLE
      FFF = EPSLONCA(:,:,KP)**1.5
      CALL DFFFDPSI(0)
      DEPSALPHADPSI(:,:,KP) = DFFF
      
      FFF = DEPSALPHADPSI(:,:,KP)
      CALL DFFFDPSI(0)
      DEPSCDPSI2(:,:,KP) = DFFF

C     COMPUTE LOGARITHMIC DERIVATIVE D LN(EPSILON_ALPHA)/D PSI
      FFF = EPSALPHA(:,:,KP)
      CALL DFFFDPSI(1)
      DEPSALPHADPSI(:,:,KP) = DFFF
      
C     COMPUTE 2ND ORDER LOGARITHMIC DERIVATIVE D^2 LN(EPSILON_ALPHA)/D PSI^2
      FFF = DEPSALPHADPSI(:,:,KP)
      CALL DFFFDPSI(0)
      DEPSALPHADPSI2(:,:,KP) = DFFF
      
      ENDDO

C     COMPUTE DIAMAGNETIC FREQUENCIES ASSOCIATED WITH DENSITY AND 
C     TEMPERATURE GRADIENTS
      OMEGASN = 0.
      OMEGAST = 0.
      DO KP=1,NSPECIES
         FFF = ESPECIES_DEN(:,:,KP)
         CALL DFFFDPSI(1)
         OMEGASN(:,:,KP) =-DFFF*B0K/OMEGACI0*ESPECIES_TEM(:,:,KP)*
     &                     ESPECIES_Z(1)/ESPECIES_Z(KP)

         FFF = ESPECIES_TEM(:,:,KP)
         CALL DFFFDPSI(0)
         OMEGAST(:,:,KP)=-DFFF*B0K/OMEGACI0*ESPECIES_Z(1)/ESPECIES_Z(KP)
      ENDDO

C     COMPUTE ADDITIONAL DIAMAGNETIC FREQUENCIES FOR HOT IONS 
      OMEGASNA = OMEGASN

      DO KP=1,NSPECIES
         FFF = AAK(:,:,KP)
         CALL DFFFDPSI(1)
         OMEGASAA(:,:,KP) =-DFFF*B0K/OMEGACI0*ESPECIES_TEM(:,:,KP)
     &                      *ESPECIES_Z(1)/ESPECIES_Z(KP)
         FFF = EPSLONCA(:,:,KP)**1.5
         CALL DFFFDPSI(0)
         OMEGASCA(:,:,KP) = DFFF*B0K/OMEGACI0*ESPECIES_TEM(:,:,KP)
     &                      *ESPECIES_Z(1)/ESPECIES_Z(KP)
      ENDDO

C     OUTPUT PARAMETRIZED DISTRIBUTION FUNCTION FOR HAGIS
C     DATA IN MARS-K UNITS, WITH THE FOLLOWING NORMALIZATION FACTORS
C     H1=P0=B0^2/MU0                           [Pa]
C     H2=N0=OMEGACI0^2*M0/(MU0*R0^2*Z0^2*e^2)  [1/m^3]
C     H3=Z0=Zi                                 [-]
C     H4=M0=Mi                                 [kg]
C     H5=T0/e=P0/N0/e                          [J/e]=[eV]
C     H6=C0=N0*(M0)^1.5                        [kg^1.5/m^3]
      IF (KCHECK.EQ.1) THEN
      H3 = ESPECIES_Z(1)
      H4 = ESPECIES_M(1)*1.6726E-27
      H1 = B0EXP**2/(4.0E-7*PI)
      H2 = H4*(OMEGACI0/R0EXP/H3/1.6022E-19)**2/4.0E-7/PI
      H5 = H1/H2/1.6022E-19
      H6 = H2*H4**1.5
      OPEN(CHOUTP,FILE='EP_DISTRIB_ANA.OUT',FORM='FORMATTED')
      DO KP=3,NSPECIES
         IF (ISPECIES_F0(KP).GE.1) THEN
            WRITE(CHOUTP,*) '%ISOTROPIC SLOWING DOWN DISTRIBUTION:' 
            WRITE(CHOUTP,*) '%SQRT(PSI_P)[-] EPS_A[eV] C[kg^1.5/m^3]
     & EPS_C[eV]'
            DO JS=1,NRP1
               WRITE(CHOUTP,120) CS(JS),EPSALPHA(JS,1,KP)*H5,
     &         ESPECIES_DEN(JS,1,KP)*AAK(JS,1,KP)*
     &         (ESPECIES_M(KP)/ESPECIES_M(1))**1.5*H6,
     &         EPSLONCA(JS,1,KP)*EPSALPHA(JS,1,KP)*H5
            ENDDO
         ENDIF        
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(4(E13.5,1X))
      ENDIF

C     DIAGNOSTIC OUTPUT
      IF (ISMPIRUN.EQ.0.AND.KCHECK.EQ.1) THEN
         WRITE(*,*) 'KDIAMAG: CS ZDHMAXDPSI ZDHMAXDPSI2 HMIN3DPSI 
     & HMIN3DPSI2 ZDHMINMAXDPSI ZDHMINMAXDPSI2'
         DO JS=1,NRP1
            WRITE(*,110) CS(JS),ZDHMAXDPSI(JS,1),
     &                   ZDHMAXDPSI2(JS,1),HMIN3DPSI(JS,1),
     &                   HMIN3DPSI2(JS,1),ZDHMINMAXDPSI(JS,1),
     &                   ZDHMINMAXDPSI2(JS,1)
         ENDDO

         DO KP=1,NSPECIES
         WRITE(*,*) 'KDIAMAG: KP,IF0TYPE=',KP,ISPECIES_F0(KP)
         WRITE(*,*) 'CS ESPECIES_DEN ESPECIES_PRE ESPECIES_PREP 
     & ALPHAA1 ALPHAA3'
         DO JS=1,NRP1
            WRITE(*,110) CS(JS),ESPECIES_DEN(JS,1,KP),
     &                   ESPECIES_PRE(JS,1,KP),ESPECIES_PREP(JS,1,KP),
     &                   ALPHAA1(JS,1,KP),ALPHAA3(JS,1,KP)
         ENDDO

         WRITE(*,*) 'CS CPSI DCDPSIL DCDPSIL2 EPSALPHA DEPSALPHADPSI
     & DEPSALPHADPSI2'
         DO JS=1,NRP1
            WRITE(*,110) CS(JS),CPSI(JS,1,KP),DCDPSIL(JS,1,KP),
     &                   DCDPSIL2(JS,1,KP),EPSALPHA(JS,1,KP),
     &                   DEPSALPHADPSI(JS,1,KP),DEPSALPHADPSI2(JS,1,KP)
         ENDDO

         WRITE(*,*) 'CS EPSLONCA OMEGASCA DEPSCDPSI2' 
         DO JS=1,NRP1
            WRITE(*,110) CS(JS),EPSLONCA(JS,1,KP),OMEGASCA(JS,1,KP),
     &                   DEPSCDPSI2(JS,1,KP)
         ENDDO

         ENDDO
      ENDIF
 110  FORMAT(9(E13.5,1X))

      RETURN
      END
      
C=======================================================================
C COMPUTE COEFFICIENTS DUE TO KINETIC PRESSURE TERMS, AT GIVEN SURFACE =
C AND FOR ALL FOURIER HARMONICS (MROW, MSA).                           =
C KGRID = 1: INTEGER RADIAL GRID                                       = 
C         2: HALF-INTEGER RADIAL GRID                                  =
C PERFORM NUMERICAL INTEGRATION OVER LAMBDA, USING VARIOUS INTEGRATION =
C SCHEMES ACCORDING TO THE FOLLOWING:                                  =
C WK = 1.0: MIDPOINT                                                   =
C      0.0: TRAPIZOIDAL                                                =
C      1/SQRT(3): GAUSSIAN (THE DEFAULT)                               =
C      1/E: LOG-TYPE                                                   =
C NOTE THE SUBROUTINE KI SHOULD = SUM_{E,I}P*I_L                       =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KJPCOEFF(JS,JS_MAT,KGRID)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      USE KineticExtendAnalysis
      USE MPIENV
      IMPLICIT NONE

      INTEGER    K,M,JS,KP,KGRID,J,L,N,JS_MAT,M00,J00,J2M,J1M,J1P,J2P,
     &           KANISO_ADIA
      REAL*8     LAMH

C     SFD DECLARATION
      INTEGER    KSFD,INCSFD
      REAL*8     AOMEGABPN,AOMEGABTN,AOMEGADIN,AOMEGADEN,AOMEGADAN,
     &           AFRACTP,AFRACTT,ANEOFUNCN
      REAL*8,DIMENSION(:),ALLOCATABLE:: AOMEGABPNSURF,AOMEGABTNSURF,
     &           AOMEGADINSURF,AOMEGADENSURF,
     &           AOMEGADANSURF,AFRACTPSURF,AFRACTTSURF
      REAL*8     ONN,ZRN,energyPartFNum,energyPartFDen
      REAL*8,DIMENSION(:),ALLOCATABLE:: HKN,JACN,HKMN
      REAL*8,DIMENSION(:),ALLOCATABLE:: LAMN,CWN
      type(SurfaceFractionDistribution),POINTER:: SFD
  
      INTEGER KCHECK

      KCHECK = 0

      INCSFD = 0
C     NOTE THAT IF INCSFD=0, SFD WILL NOT BE CALLED
C     AND AVERAGED DRIFT FREQUENCIES WILL BE COMPUTED BASED ON 
C     ANALYTIC FORMULA USING RSS AND RUU
      
C     SFD OPERATION
C     ============================================================
      IF (INCSFD.GT.0) THEN
      SFD   => NULL()
      KSFD  =0
      ENDIF

      AOMEGABPN = 0.
      AOMEGABTN = 0.
      AOMEGADIN = 0.
      AOMEGADEN = 0.
      AOMEGADAN = 0.
      AFRACTP   = 0.
      AFRACTT   = 0.
      ANEOFUNCN = 0.

      IF (ODWKCOM) THEN
      VX1PARAC = 0.
      VX1PERPC = 0.
      VX2PARAC = 0.
      VX2PERPC = 0.
      VQ1PARAC = 0.
      VQ1PERPC = 0.
      VQ2PARAC = 0.
      VQ2PERPC = 0.
      VQ3PARAC = 0.
      VQ3PERPC = 0.
      ENDIF
      
      IF (KGRID.EQ.1) THEN
         ZRN = 0.0
         ONN = 1.0
         IF(.NOT.ALLOCATED(HKN)) ALLOCATE(HKN(NCHI), JACN(NCHI), 
     &                                    HKMN(NCHI))
         DO J=1,NCHI
            HKN(J)  = HK(JS,J,KGRID)
            JACN(J) = RJA(JS,J)
            HKMN(J) = HKMIN(JS,KGRID) 
         ENDDO
         ALLOCATE(AOMEGABPNSURF(NCHI),AOMEGABTNSURF(NCHI),
     &          AOMEGADINSURF(NCHI),AOMEGADENSURF(NCHI),
     &          AOMEGADANSURF(NCHI),AFRACTPSURF(NCHI),AFRACTTSURF(NCHI))
         AOMEGABPNSURF = 0.
         AOMEGABTNSURF = 0.
         AOMEGADINSURF = 0.
         AOMEGADENSURF = 0.
         AOMEGADANSURF = 0.
         AFRACTPSURF = 0.
         AFRACTTSURF = 0.
      ENDIF
C     ============================================================

C     CHECK INPUT PARAMETERS
      IF (KGRID.NE.1.AND.KGRID.NE.2) STOP 'CHECK KJPCOEFF: KGRID'

C     COMPUTE EQUILIBRIUM QUANTITIES ON RCHI GRID
      CALL KEQUIL(JS,KGRID)

C     COMPUTE LOGARITHMIC DERIVATIVE X'/X FOR POLOIDAL
C     FOURIER HARMONICS, TO BE USED FOR FOW CORRECTIONS
C     NOTE HERE X'=DX/DPSI
      IF ((IFOWP.EQ.1.OR.IFOWT.EQ.1).AND.JS.GT.1) CALL KXLNP(JS,KGRID)

C     SET GLOBAL INDEX KNUMDISTRIB FOR COUNTING CALLS OF KDISTRIBF
      KNUMDISTRIB = 0

      IF (IPARTICLE.EQ.1.OR.IPARTICLE.EQ.2.OR.KANISOTROPIC.EQ.1) THEN
C     GO THROUGH PASSING PARTICLES
      IF (JS.EQ.JS0.AND.KGRID.EQ.1
     &    .AND.(ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT)) 
     &   WRITE(*,*) 'PASSING PARTICLES...'

C     DEFINE INTEGRATION POINTS ALONG CHI
      CHIL = 0.0
      CHIU = 2.0*PI

C     DEFINE INTEGRATION POINTS ALONG CHI, BETWEEN CHIL AND CHIU
      CALL KCHI(1)

C     COMPUTE EQUILIBRIUM QUANTITIES ON RCHIK GRID
      CALL KEQUILK(JS,KGRID)

C     COMPUTE PHI(CHI) AT SPECIFIC INTEGRATION POINTS 
      CALL KPHI(JS,KGRID)

      IF (KANISOTROPIC.EQ.1.OR.IFOWP.EQ.1) CALL ZKEQUILK(JS,KGRID)

C     COMPUTE LAMBDA-MESH INCLUDING TWO END POINTS
      DO J=1,NLAMK1(JS,KGRID)-1
      LAMHH(J+J-1) = LAMK1(JS,J+1,KGRID)-LAMK1(JS,J,KGRID)
      LAMHH(J+J)   = LAMHH(J+J-1)
      DO N=0,1
         LAMM(J+J+N) = ((1+WK)*LAMK1(JS,J+N,KGRID)+
     &                  (1-WK)*LAMK1(JS,J-N+1,KGRID))/2.
      ENDDO
      ENDDO
      LAMM(1) = 0.
      LAMM(2*NLAMK1(JS,KGRID)) = HKMIN(JS,KGRID)

      IF (IPARTICLE.EQ.1.OR.IPARTICLE.EQ.2) THEN

C     FIND SINGULAR POINTS IN LAMBDA-SPACE 
      CALL KLAM0(JS,KGRID,1)

C     COMPUTE SF0-FACTOR 
      SF0 = 0.
      IF (NKSINGULAR.EQ.1) THEN
      CALL KSF0(JS,KGRID,1,0)
      IF (IFOWP.EQ.1) THEN
         IF (IFOWPSI0.EQ.1) CALL KSF0(JS,KGRID,1,1)
         CALL KSF0(JS,KGRID,1,2)
         CALL KSF0(JS,KGRID,1,3)
      ENDIF
      ENDIF

C     COMPUTE SG0,SH0-FACTORS FOR SINGULARITY SUBTRACTION IN PITCH ANGLE INTEGRATION
      CALL KSGH0(JS,KGRID,1)

      ENDIF

      DO J=1,2*NLAMK1(JS,KGRID)-2
         LAMH = LAMHH(J)/2./B0K/SQRT(PI)
         LAM = LAMM(J+1)

C        COMPUTE AND STORE F0 AND ITS DERIVATIVES
         CALL KDISTRIBF(JS,KGRID,LAM)

C        COMPUTE (NORMALIZED) TRANSIT TIME AS A FUNCTION OF CHI
C        AT SPECIFIC INTEGRATION POINTS 
         CALL KBTIME(JS,KGRID,1)

C        COMPUTE (NORMALIZED) TRANSIT FREQUENCY
         OMEGAB = 2.0*PI/RTK(NCHI2+2)
         
         IF (IPARTICLE.NE.1.AND.IPARTICLE.NE.2) GOTO 211

C        DEFINE PRECESSION DRIFT WHICH IS NOT IN USE
         DRIFT = 0.0

C        PARTCILE ENERGY INTEGRATION AT ZERO ORBIT WIDTH
         CALL KI(JS,KGRID,1,0)

C        COMPUTE G-FACTORS
         CALL KG(1)

C        COMPUTE H-FACTORS
         CALL KH(JS,KGRID,1)

         IF (IFOWP.EQ.1) THEN
            TPSI0L     = TPSI0    (JS,J,KGRID)
            TPSI0DPSIL = TPSI0DPSI(JS,J,KGRID)
            TPSI0DLAML = TPSI0DLAM(JS,J,KGRID)
            HPSI0L     = HPSI0    (JS,J,KGRID)
            HPSI0DPSIL = HPSI0DPSI(JS,J,KGRID)
            HPSI0DLAML = HPSI0DLAM(JS,J,KGRID)

            CALL KI(JS,KGRID,1,1)
            CALL KI(JS,KGRID,1,2)
            CALL KI(JS,KGRID,1,3)
            CALL KG1(1)
            CALL KH1(JS,KGRID,1)
         ENDIF

         IF (IFOWP.EQ.0) THEN
            CALL KJPFILL (JS,JS_MAT,KGRID,LAMH,1,1)
            CALL KJPFILL (JS,JS_MAT,KGRID,LAMH,1,2)
         ELSE
            CALL KJPFILL1(JS,JS_MAT,KGRID,LAMH,1,1)
            CALL KJPFILL1(JS,JS_MAT,KGRID,LAMH,1,2)
         ENDIF

 211     CONTINUE         
         IF (KGRID.EQ.1) THEN
C           NOTE THE EXTRA FACTOR (B0K*SQRT(PI)) IN LAMH
            AOMEGABPN = AOMEGABPN + OMEGAB*RUU/RSS*LAMH*B0K
            ANEOFUNCN = ANEOFUNCN + LAM/RUU2*LAMH*B0K*SQRT(PI)
         ENDIF

C        ADD ADIABATIC CONTRIBUTION FROM PASSING PARTICLES
         ZKIA = 0.
         IF (KANISOTROPIC.EQ.1) 
     &      CALL KIA_ADIABATIC(JS,KGRID,0,LAM,ZKIA)
         IF (KANISOTROPIC.EQ.1.OR.IFOWP.EQ.1) 
     &      CALL KG_ADIABATIC(JS,KGRID,1)
         IF (IFOWP.EQ.1) THEN
            TPSI0L     = TPSI0    (JS,J,KGRID)
            TPSI0DPSIL = TPSI0DPSI(JS,J,KGRID)
            TPSI0DLAML = TPSI0DLAM(JS,J,KGRID)
            HPSI0L     = HPSI0    (JS,J,KGRID)
            HPSI0DPSIL = HPSI0DPSI(JS,J,KGRID)
            HPSI0DLAML = HPSI0DLAM(JS,J,KGRID)

            CALL KIA_ADIABATIC(JS,KGRID,1,LAM,ZKIA1)
            CALL KIA_ADIABATIC(JS,KGRID,2,LAM,ZKIA2)
            CALL KG_ADIABATIC1(JS,KGRID,1)
         ENDIF
         IF (KANISOTROPIC.EQ.1.AND.IFOWP.EQ.0) 
     &      CALL KJPFILL(JS,JS_MAT,KGRID,LAMH,1,5)
         IF (IFOWP.EQ.1) CALL KJPFILL1(JS,JS_MAT,KGRID,LAMH,1,5)
      ENDDO
      
      IF (IPARTICLE.NE.1.AND.IPARTICLE.NE.2) GOTO 212

C     ADD SPECIAL SINGULAR CONTRIBUTIONS 
      IF (IFOWP.EQ.0) CALL KJPFILL (JS,JS_MAT,KGRID,0.,1,3)
      IF (IFOWP.EQ.1) CALL KJPFILL1(JS,JS_MAT,KGRID,0.,1,3)

C     SFD OPERATION
C     PASSING THERMAL ION BOUNCE FREQUENCY AVERAGE
C     AND PASSING THERMAL/FAST ION FRACTION
C     INDEPENDENT OF PARTICLE ENERGY DISTRIBUTION
      IF (KGRID.EQ.1.AND.INCSFD.GT.0) THEN
         K = NLAMK1(JS,KGRID)*2-2
         ALLOCATE(LAMN(K),CWN(K))
         LAMN = LAMM(2:K+1)
         CWN  = ZOMEGABP(JS,2:K+1,KGRID)
         energyPartFNum = 1.0
         energyPartFDen = SQRT(PI)*0.5
         CALL calcSurfacefractionDistribution(KSFD,SFD,ONN,ONN,
     &        LAMN,CWN,JACN,ZRN,HKMN,ZRN,HKMN,HKN)
         AOMEGABPN = SFD%averfraction*energyPartFNum/energyPartFDen
         AOMEGABPNSURF = SFD%fractionDistribution(:,1)
     &                  *energyPartFNum/energyPartFDen
         POSITIONKAI (JS,:) = SFD%postionKai(:,1)

         energyPartFNum = SQRT(PI)*0.5
         energyPartFDen = SQRT(PI)*0.5
         DO K=1,NLAMK1(JS,KGRID)*2-2
            CWN(K) = 1.0
         ENDDO
         CALL calcSurfacefractionDistribution(KSFD,SFD,ONN,ONN,
     &        LAMN,CWN,JACN,ZRN,HKMN,ZRN,HKN,HKN)
         AFRACTP = SFD%averfraction*energyPartFNum/energyPartFDen
         AFRACTPSURF = SFD%fractionDistribution(:,1)
     &                  *energyPartFNum/energyPartFDen

         DEALLOCATE(LAMN,CWN)
      ENDIF

 212  CONTINUE
      
      ENDIF

      IF (IPARTICLE.EQ.1.OR.IPARTICLE.EQ.3.OR.KANISOTROPIC.EQ.1) THEN
C     GO THROUGH TRAPPED  PARTICLES
C     SPECIAL TREATMENT REQUIRED FOR THE SINGULARITY AT DRIFT=0
C     FOR PARTICLE PRECESSION DRIFT RESONANCE
      IF (JS.EQ.JS0.AND.KGRID.EQ.1
     &    .AND.(ISMPIRUN.EQ.0 .OR. RANK.EQ.ROOT)) 
     &   WRITE(*,*) 'TRAPPED PARTICLES...'
      
C     COMPUTE LAMBDA-MESH INCLUDING TWO END POINTS
C     ALSO SAVED DRIFT(LAM)
      DO J=1,NLAMK0(JS,KGRID)-1
      LAMHH(J+J-1) = LAMK0(JS,J+1,KGRID)-LAMK0(JS,J,KGRID)
      LAMHH(J+J)   = LAMHH(J+J-1)
      DO N=0,1
         LAMM(J+J+N) = ((1+WK)*LAMK0(JS,J+N,KGRID)+
     &                  (1-WK)*LAMK0(JS,J-N+1,KGRID))/2.
      ENDDO
      ENDDO
      LAMM(1) = HKMIN(JS,KGRID)
      LAMM(2*NLAMK0(JS,KGRID)) = HKMAX(JS,KGRID)

      IF (IPARTICLE.EQ.1.OR.IPARTICLE.EQ.3) THEN

C     FIND SINGULAR POINTS IN LAMBDA-SPACE 
      CALL KLAM0(JS,KGRID,0)

C     FIND LAM0 FOR DRIFT=0 USING SPLINE
      IF (SUM(ABS(PSPECIES_NTD)).GT.0) CALL KLAM0(JS,KGRID,2)
      
C     COMPUTE SF0-FACTOR FOR SINGULARITY SUBTRACTION IN PITCH ANGLE INTEGRATION
      SF0 = 0.
      IF (NKSINGULAR.EQ.1) THEN
      CALL KSF0(JS,KGRID,0,0)
      IF (IFOWT.EQ.1) THEN
         CALL KSF0(JS,KGRID,0,2)
         CALL KSF0(JS,KGRID,0,3)
      ENDIF
      ENDIF

C     COMPUTE SG0,SH0-FACTORS FOR SINGULARITY SUBTRACTION IN PITCH ANGLE INTEGRATION
      CALL KSGH0(JS,KGRID,0)

C     INTEGRATE I-FACTOR ON A DENSE LAMBDA MESH FOR RLM(L)=0
      IF (SUM(ABS(PSPECIES_NTD)).GT.0.AND.SLAMD0.GT.0.) THEN   
         CALL KI0(JS,KGRID,0)
         IF (IFOWT.EQ.1) THEN
            CALL KI0(JS,KGRID,2)
            CALL KI0(JS,KGRID,3)
         ENDIF
      ENDIF

      ENDIF

      DO J=1,2*(NLAMK0(JS,KGRID)-1)
         LAMH = LAMHH(J)/2./B0K/SQRT(PI)
         LAM  = LAMM(J+1)

C        COMPUTE AND STORE F0 AND ITS DERIVATIVES
         CALL KDISTRIBF(JS,KGRID,LAM)

C        COMPUTE POLOIDAL ANGLES OF TURNING POINTS
         CALL KTURN(JS,KGRID)

C        DEFINE INTEGRATION POINTS ALONG CHI, BETWEEN CHIL AND CHIU
         CALL KCHI(0)
       
C        COMPUTE EQUILIBRIUM QUANTITIES AT INTEGRATION POINTS ALONG CHI
         CALL KEQUILK(JS,KGRID)
         
         IF (KANISOTROPIC.EQ.1) CALL ZKEQUILK(JS,KGRID)

         IF (IPARTICLE.NE.1.AND.IPARTICLE.NE.3) GOTO 312

C        COMPUTE (NORMALIZED) BOUNCING TIME AS A FUNCTION OF CHI
C        AT SPECIFIC INTEGRATION POINTS 
         CALL KBTIME(JS,KGRID,0)

C        COMPUTE (NORMALIZED) BOUNCING FREQUENCY
C        NOTE RTK IS COMPUTED FOR HALF BOUNCE PERIOD FOR TRAPPED PARTICLES
         OMEGAB = PI/RTK(NCHI2+2)
         
C        GET COEFFICIENT OF DRIFT PRECESSION FREQUENCY
         DRIFT = ZOMEGADT(JS,J+1,KGRID)

C        COMPUTE PHI(CHI) AT SPECIFIC INTEGRATION POINTS 
         CALL KPHI(JS,KGRID)

C        ENERGY INTEGRATION AT ZERO ORBIT WIDTH
         CALL KI(JS,KGRID,0,0)

C        COMPUTE G-FACTORS
         CALL KG(0)

C        COMPUTE H-FACTORS
         CALL KH(JS,KGRID,0)

         IF (IFOWT.EQ.1) THEN
            CALL KI(JS,KGRID,0,2)
            CALL KI(JS,KGRID,0,3)
            CALL KG1(0)
            CALL KH1(JS,KGRID,0)
         ENDIF

         IF ((IPARTICLE.EQ.1.OR.IPARTICLE.EQ.3).AND.IFOWT.EQ.0) THEN
            CALL KJPFILL (JS,JS_MAT,KGRID,LAMH,0,1)
            CALL KJPFILL (JS,JS_MAT,KGRID,LAMH,0,2)
         ENDIF
         IF ((IPARTICLE.EQ.1.OR.IPARTICLE.EQ.3).AND.IFOWT.EQ.1) THEN
            CALL KJPFILL1(JS,JS_MAT,KGRID,LAMH,0,1)
            CALL KJPFILL1(JS,JS_MAT,KGRID,LAMH,0,2)
         ENDIF

 312     CONTINUE

         IF (J.GT.0.AND.KGRID.EQ.1) THEN
C           NOTE THE EXTRA FACTOR (B0K*SQRT(PI)) IN LAMH
            AOMEGABTN = AOMEGABTN + OMEGAB*RUU/RSS*LAMH*B0K
            AOMEGADIN = AOMEGADIN + DRIFT*RUU/RSS*LAMH*0.75*B0K*SQRT(PI)
         ENDIF

C        ADD ADIABATIC CONTRIBUTION FROM ANISOTROPIC DISTRIBUTION OF
C        TRAPPED ENERGETIC PARTICLES
C        NOTE FOW ADIABATIC CONTRIBUTION FROM TRAPPED PARTICLES VANISHES 
         ZKIA = 0.
         IF (KANISOTROPIC.EQ.1) THEN
            CALL KIA_ADIABATIC(JS,KGRID,0,LAM,ZKIA)
            CALL KG_ADIABATIC(JS,KGRID,0)
            CALL KJPFILL(JS,JS_MAT,KGRID,LAMH,0,5)
         ENDIF
      ENDDO

      IF (IPARTICLE.NE.1.AND.IPARTICLE.NE.3) GOTO 314

C     ADD SPECIAL SINGULAR CONTRIBUTION FROM RLM(L)=0
      IF (IFOWT.EQ.0) THEN
         CALL KJPFILL (JS,JS_MAT,KGRID,0.,0,3)
         IF (SLAMD0.GT.0.) CALL KJPFILL (JS,JS_MAT,KGRID,0.,0,4)
      ELSE
         CALL KJPFILL1(JS,JS_MAT,KGRID,0.,0,3)
         IF (SLAMD0.GT.0.) CALL KJPFILL1(JS,JS_MAT,KGRID,0.,0,4)
      ENDIF

C     SFD OPERATION
C     TRAPPED THERMAL ION BOUNCE FREQUENCY AVERAGE
      IF (KGRID.EQ.1.AND.INCSFD.GT.0) THEN
         K = NLAMK0(JS,KGRID)*2-2
         ALLOCATE(LAMN(K),CWN(K))
         LAMN = LAMM(2:K+1)
         CWN  = ZOMEGABT(JS,2:K+1,KGRID)
         energyPartFNum = 1.0
         energyPartFDen = SQRT(PI)*0.5
         CALL calcSurfacefractionDistribution(KSFD,SFD,ONN,ONN,
     &        LAMN,CWN,JACN,HKMN(1),HKN,HKMN(1),HKN,HKN)
         AOMEGABTN = SFD%averfraction*energyPartFNum/energyPartFDen
         AOMEGABTNSURF = SFD%fractionDistribution(:,1)
     &                        *energyPartFNum/energyPartFDen
         POSITIONKAI (JS,:) = SFD%postionKai(:,1)

C        TRAPPED THERMAL ION PRECESSION FREQUENCY AVERAGE
         energyPartFNum = SQRT(PI)*0.75
         energyPartFDen = SQRT(PI)*0.5
         CWN = ZOMEGADT(JS,2:K+1,KGRID)
         CALL calcSurfacefractionDistribution(KSFD,SFD,ONN,ONN,
     &        LAMN,CWN,JACN,HKMN(1),HKN,HKMN(1),HKN,HKN)
         AOMEGADIN = SFD%averfraction*energyPartFNum/energyPartFDen
         AOMEGADINSURF = SFD%fractionDistribution(:,1)
     &                        *energyPartFNum/energyPartFDen

C        TRAPPED THERMAL/FAST ION FRACTION
C        INDEPENDENT OF PARTICLE ENERGY DISTRIBUTION
         energyPartFNum = SQRT(PI)*0.5
         energyPartFDen = SQRT(PI)*0.5
         CWN = 1.
         CALL calcSurfacefractionDistribution(KSFD,SFD,ONN,ONN,
     &        LAMN,CWN,JACN,HKMN(1),HKN,ZRN,HKN,HKN)
         AFRACTT = SFD%averfraction*energyPartFNum/energyPartFDen
         AFRACTTSURF = SFD%fractionDistribution(:,1)
     &                        *energyPartFNum/energyPartFDen
         DEALLOCATE(LAMN,CWN)
      ENDIF

 314  CONTINUE
      ENDIF

C     ADD EXTRA TERM TO ADIABATIC CONTRIBUTION, 
C     ASSOCIATED WITH ISOTROPIC SLOWING DOWN DISTRIBUTION 
      IF (KFASTRUN.EQ.1) CALL KJPFILL(JS,JS_MAT,KGRID,0.,1,6)

C     SFD FINAL CALIBRATION OF BOUNCE & PRECESSION FREQUNCIES
C     FIND RESONANCE CONDITION FOR HOT IONS IN LOCAL SPACE    
C     SAVE ALL FREQUNCIES INTO A FILE
      IF (KGRID.EQ.1.AND.IPARTICLE.GT.0) THEN
         AOMEGABPN = AOMEGABPN*SQRT(2.0*ESPECIES_TEM(JS,1,1))
         AOMEGADAN = 0.0
         IF (NSPECIES.GE.3) THEN
            AOMEGADAN   = AOMEGADIN*ESPECIES_TEM(JS,1,2)*B0K/OMEGACI0
     &                    *ALPHAA3(JS,1,3)/ALPHAA1(JS,1,3)
     &                    *ESPECIES_Z(3)/ESPECIES_Z(1)
            FREQK(JS,6) = AOMEGABTN*SQRT(2.0*ESPECIES_TEM(JS,1,3)
     &                                   /ESPECIES_M(3))
         ENDIF
         IF (NSPECIES.GE.4) THEN
            FREQK(JS,7) = AOMEGADIN*ESPECIES_TEM(JS,1,2)*B0K/OMEGACI0
     &                    *ALPHAA3(JS,1,4)/ALPHAA1(JS,1,4)
            FREQK(JS,8) = AOMEGABTN*SQRT(2.0*ESPECIES_TEM(JS,1,4)
     &                                   /ESPECIES_M(4))
         ENDIF
         AOMEGABTN = AOMEGABTN*SQRT(2.0*ESPECIES_TEM(JS,1,1)
     &                              /ESPECIES_M(1))
         AOMEGADEN =-AOMEGADIN*ESPECIES_TEM(JS,1,2)*B0K/OMEGACI0
         AOMEGADIN = AOMEGADIN*ESPECIES_TEM(JS,1,1)*B0K/OMEGACI0
         ANEOFUNCN = ANEOFUNCN*0.75

         AOMEGABPNSURF = AOMEGABPN*SQRT(2.0*ESPECIES_TEM(JS,1,1))
         AOMEGABTNSURF = AOMEGABTN*SQRT(2.0*ESPECIES_TEM(JS,1,1))
         AOMEGADENSURF =-AOMEGADIN*ESPECIES_TEM(JS,1,2)*B0K/OMEGACI0
         AOMEGADANSURF = 0.0
         IF (NSPECIES.GT.2)
     &   AOMEGADANSURF = AOMEGADIN*ESPECIES_TEM(JS,1,2)*B0K/OMEGACI0
     &                   *ALPHAA3(JS,1,3)/ALPHAA1(JS,1,3)
         AOMEGADINSURF = AOMEGADIN*ESPECIES_TEM(JS,1,1)*B0K/OMEGACI0
         LAMH = 0.
         DO L=1,MLMAX
            IF (ABS(RLM(L)).LT.0.1) LAMH = SLAM0(L,3)
         ENDDO

         FREQK(JS,1)  = AOMEGABPN
         FREQK(JS,2)  = AOMEGABTN
         FREQK(JS,3)  = AOMEGADIN
         FREQK(JS,4)  = AOMEGADEN
         FREQK(JS,5)  = AOMEGADAN
         IF (NSPECIES.EQ.2) THEN
         FREQK(JS,6)  = AFRACTP  
         FREQK(JS,7)  = AFRACTT  
         FREQK(JS,8)  = LAMM(1)
         ENDIF
         FREQK(JS,9)  = LAMH
         FREQK(JS,10) = LAMM(NLAMK0(JS,KGRID)*2)
         FREQK(JS,11) = ZOMEGADT(JS,1,KGRID)
         FREQK(JS,12) = ZOMEGADT(JS,NLAMK0(JS,KGRID)*2,KGRID)
         FREQK(JS,13) = ANEOFUNCN

         IF (ISMPIRUN.EQ.0) THEN
            FREQKSURF(JS,1,:) = AOMEGABPNSURF
            FREQKSURF(JS,2,:) = AOMEGABTNSURF
            FREQKSURF(JS,3,:) = AOMEGADINSURF
            FREQKSURF(JS,4,:) = AOMEGADENSURF
            FREQKSURF(JS,5,:) = AOMEGADANSURF
            FREQKSURF(JS,6,:) = AFRACTPSURF
            FREQKSURF(JS,7,:) = AFRACTTSURF
         ENDIF
      ENDIF



C     SFD RELEASE
      IF (INCSFD.GT.0) CALL releaseSurfaceFractionDistribution(SFD)
      IF (KGRID.EQ.1) THEN
         DEALLOCATE(HKN,JACN,HKMN)
         DEALLOCATE(AOMEGABPNSURF,AOMEGABTNSURF,
     &          AOMEGADINSURF,AOMEGADENSURF,
     &          AOMEGADANSURF,AFRACTPSURF,AFRACTTSURF)


      ENDIF

      IF (INCSFD.GT.0) CALL WRITE_SURFACE_QUANTITIES(JS,KGRID)
      
      RETURN
      END

C=======================================================================
C PUT ALL I,H,G FACTORS TOGETHER AND FILL INTO MATRIX ELEMENTS 
C ICASE = 1: REGULAR NON-ADIABATIC CONTRIBUTION (BOTH PASS&TRAP)
C         2: PITCH ANGLE DEPEDENT SINGULAR PART (BOTH PASS&TRAP)
C         3: PITCH ANGLE INDEPEDNET SINGULAR PART (BOTH PASS&TRAP)
C         4: SPECIAL SINGULAR CONTRIBUTION FROM RLM(L)=0 FOR 
C            PRECESSIONAL DRIFT RESONANCE OF TRAPPED PARTICLES
C         5: ADIABATIC CONTRIBUTION FROM ANISOTROPIC EP (BOTH PASS&TRAP)
C         6: EXTRA TERM TO ADIABATIC CONTRIBUTION, ASSOCIATED WITH 
C            ISOTROPIC SLOWING DOWN DISTRIBUTION 
C ADDED PERTURBED ELECTROSTATIC POTENTIAL DPHI BY Y.Q.LIU IN JUNE 2023 
C NOTE THAT ICASE=5,6 DOES NOT APPLY TO X1PARAE ETC.      
C=======================================================================
      SUBROUTINE KJPFILL(JS,JS_MAT,KGRID,LAMH,KPARTICLE,ICASE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      IMPLICIT NONE

      INTEGER    JS,JS_MAT,KGRID,KPARTICLE,ICASE,R,
     &           K,M,KP,J,L,N,M00,J00,J2M,J1M,J1P,J2P
      REAL*8     LAMH,LAM1,LAM2,H1,H3,H4
      COMPLEX*16 CTMP,CTM2,CTM3,CTM4
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE:: 
     &           X1PARA,X2PARA,Q1PARA,Q2PARA,Q3PARA,DPPARA,
     &           X1PERP,X2PERP,Q1PERP,Q2PERP,Q3PERP,DPPERP,
     &           X1DPHI,X2DPHI,Q1DPHI,Q2DPHI,Q3DPHI,DPDPHI
      INTEGER KCHECK
      KCHECK = 0

      LAM1   = HKMIN(JS,KGRID)
      LAM2   = HKMAX(JS,KGRID)
      H3     = WFUN (JS,KGRID)

      ALLOCATE ( X1PARA(NSPECIES,2),X2PARA(NSPECIES,2),
     &           Q1PARA(NSPECIES,2),Q2PARA(NSPECIES,2),
     &           Q3PARA(NSPECIES,2),X1PERP(NSPECIES,2),
     &           X2PERP(NSPECIES,2),Q1PERP(NSPECIES,2),
     &           Q2PERP(NSPECIES,2),Q3PERP(NSPECIES,2), 
     &           DPPARA(NSPECIES,2),DPPERP(NSPECIES,2), 
     &           X1DPHI(NSPECIES,2),X2DPHI(NSPECIES,2), 
     &           Q1DPHI(NSPECIES,2),Q2DPHI(NSPECIES,2), 
     &           Q3DPHI(NSPECIES,2),DPDPHI(NSPECIES,2) )
      DO K=1,MSMAX
      DO M=1,MSMAX

      X1PARA = 0.
      X1PERP = 0.
      X1DPHI = 0.
      X2PARA = 0.
      X2PERP = 0.
      X2DPHI = 0.
      Q1PARA = 0.
      Q1PERP = 0.
      Q1DPHI = 0.
      Q2PARA = 0.
      Q2PERP = 0.
      Q2DPHI = 0.
      Q3PARA = 0.
      Q3PERP = 0.
      Q3DPHI = 0.
      DPPARA = 0.
      DPPERP = 0.
      DPDPHI = 0.
      
      SELECT CASE (ICASE)
      CASE (1)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
            IF (KPARTICLE.EQ.0 .AND. ABS(RLM(L)).LT.0.1) THEN
               R=2
            ELSE
               R=1
            ENDIF

            CTMP   = VI(1,L,KP)*LAMH*H3
            CTM2   = VI(2,L,KP)*LAMH*H3
            CTM3   = VI(3,L,KP)*LAMH*H3
            CTM4   = VI(4,L,KP)*LAMH*H3
            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA(K,L)*VX1(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP(K,L)*VX1(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI(K,L)*VX1(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA(K,L)*VX2(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP(K,L)*VX2(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI(K,L)*VX2(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA(K,L)*VQ1(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP(K,L)*VQ1(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI(K,L)*VQ1(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA(K,L)*VQ2(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP(K,L)*VQ2(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI(K,L)*VQ2(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA(K,L)*VQ3(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP(K,L)*VQ3(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI(K,L)*VQ3(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA(K,L)*VDP(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP(K,L)*VDP(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI(K,L)*VDP(M,L)

            IF (ABS(RLM(L)).LT.0.1.AND.SLAMD0.GT.0..AND.KPARTICLE.EQ.0) 
     &      THEN
            X1PARA(KP,R) = X1PARA(KP,R) - CTMP*VPARA0(K,L)*VX10(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) - CTMP*VPERP0(K,L)*VX10(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) - CTM3*VDPHI0(K,L)*VX10(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) - CTMP*VPARA0(K,L)*VX20(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) - CTMP*VPERP0(K,L)*VX20(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) - CTM3*VDPHI0(K,L)*VX20(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) - CTMP*VPARA0(K,L)*VQ10(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) - CTMP*VPERP0(K,L)*VQ10(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) - CTM3*VDPHI0(K,L)*VQ10(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) - CTMP*VPARA0(K,L)*VQ20(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) - CTMP*VPERP0(K,L)*VQ20(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) - CTM3*VDPHI0(K,L)*VQ20(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) - CTMP*VPARA0(K,L)*VQ30(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) - CTMP*VPERP0(K,L)*VQ30(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) - CTM3*VDPHI0(K,L)*VQ30(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) - CTM2*VPARA0(K,L)*VDP0(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) - CTM2*VPERP0(K,L)*VDP0(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) - CTM4*VDPHI0(K,L)*VDP0(M,L)
            ENDIF
         ENDDO
         ENDDO
                  
      CASE (2)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
         IF (SLAM0(L,KP).GT.0.) THEN
            IF (KPARTICLE.EQ.0 .AND. ABS(RLM(L)).LT.0.1) THEN
               R=2
            ELSE
               R=1
            ENDIF
            IF (KPARTICLE.EQ.1) THEN
               H4 = PSPECIES_NP(KP)
            ELSEIF (ABS(RLM(L)).GT.0.1) THEN
               H4 = PSPECIES_NTB(KP)
            ELSE
               H4 = PSPECIES_NTD(KP)
            ENDIF
            H1     = SLAM0(L,KP)
            CTMP   = LOG(ABS(LAM-H1))*SF0(L,KP,0,1)*LAMH*H3*H4
            CTM2   = LOG(ABS(LAM-H1))*SF0(L,KP,0,2)*LAMH*H3*H4
            CTM3   = LOG(ABS(LAM-H1))*SF0(L,KP,0,3)*LAMH*H3*H4
            CTM4   = LOG(ABS(LAM-H1))*SF0(L,KP,0,4)*LAMH*H3*H4
            X1PARA(KP,R)=X1PARA(KP,R)+CTMP*VPARA(K,L)*VX1(M,L)
            X1PERP(KP,R)=X1PERP(KP,R)+CTMP*VPERP(K,L)*VX1(M,L)
            X1DPHI(KP,R)=X1DPHI(KP,R)+CTM3*VDPHI(K,L)*VX1(M,L)
            X2PARA(KP,R)=X2PARA(KP,R)+CTMP*VPARA(K,L)*VX2(M,L)
            X2PERP(KP,R)=X2PERP(KP,R)+CTMP*VPERP(K,L)*VX2(M,L)
            X2DPHI(KP,R)=X2DPHI(KP,R)+CTM3*VDPHI(K,L)*VX2(M,L)
            Q1PARA(KP,R)=Q1PARA(KP,R)+CTMP*VPARA(K,L)*VQ1(M,L)
            Q1PERP(KP,R)=Q1PERP(KP,R)+CTMP*VPERP(K,L)*VQ1(M,L)
            Q1DPHI(KP,R)=Q1DPHI(KP,R)+CTM3*VDPHI(K,L)*VQ1(M,L)
            Q2PARA(KP,R)=Q2PARA(KP,R)+CTMP*VPARA(K,L)*VQ2(M,L)
            Q2PERP(KP,R)=Q2PERP(KP,R)+CTMP*VPERP(K,L)*VQ2(M,L)
            Q2DPHI(KP,R)=Q2DPHI(KP,R)+CTM3*VDPHI(K,L)*VQ2(M,L)
            Q3PARA(KP,R)=Q3PARA(KP,R)+CTMP*VPARA(K,L)*VQ3(M,L)
            Q3PERP(KP,R)=Q3PERP(KP,R)+CTMP*VPERP(K,L)*VQ3(M,L)
            Q3DPHI(KP,R)=Q3DPHI(KP,R)+CTM3*VDPHI(K,L)*VQ3(M,L)
            DPPARA(KP,R)=DPPARA(KP,R)+CTM2*VPARA(K,L)*VDP(M,L)
            DPPERP(KP,R)=DPPERP(KP,R)+CTM2*VPERP(K,L)*VDP(M,L)
            DPDPHI(KP,R)=DPDPHI(KP,R)+CTM4*VDPHI(K,L)*VDP(M,L)

            X1PARA(KP,R)=X1PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVX10(M,L,KP)
            X1PERP(KP,R)=X1PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVX10(M,L,KP)
            X1DPHI(KP,R)=X1DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVX10(M,L,KP)
            X2PARA(KP,R)=X2PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVX20(M,L,KP)
            X2PERP(KP,R)=X2PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVX20(M,L,KP)
            X2DPHI(KP,R)=X2DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVX20(M,L,KP)
            Q1PARA(KP,R)=Q1PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVQ10(M,L,KP)
            Q1PERP(KP,R)=Q1PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVQ10(M,L,KP)
            Q1DPHI(KP,R)=Q1DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVQ10(M,L,KP)
            Q2PARA(KP,R)=Q2PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVQ20(M,L,KP)
            Q2PERP(KP,R)=Q2PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVQ20(M,L,KP)
            Q2DPHI(KP,R)=Q2DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVQ20(M,L,KP)
            Q3PARA(KP,R)=Q3PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVQ30(M,L,KP)
            Q3PERP(KP,R)=Q3PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVQ30(M,L,KP)
            Q3DPHI(KP,R)=Q3DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVQ30(M,L,KP)
            DPPARA(KP,R)=DPPARA(KP,R)-CTM2*SVPARA0(K,L,KP)*SVDP0(M,L,KP)
            DPPERP(KP,R)=DPPERP(KP,R)-CTM2*SVPERP0(K,L,KP)*SVDP0(M,L,KP)
            DPDPHI(KP,R)=DPDPHI(KP,R)-CTM4*SVDPHI0(K,L,KP)*SVDP0(M,L,KP)
         ENDIF
         ENDDO
         ENDDO

      CASE (3)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
         IF (SLAM0(L,KP).GT.0.) THEN
            IF (KPARTICLE.EQ.0 .AND. ABS(RLM(L)).LT.0.1) THEN
               R=2
            ELSE
               R=1
            ENDIF
            IF (KPARTICLE.EQ.1) THEN
               H4 = PSPECIES_NP(KP)
            ELSEIF (ABS(RLM(L)).GT.0.1) THEN
               H4 = PSPECIES_NTB(KP)
            ELSE
               H4 = PSPECIES_NTD(KP)
            ENDIF
            H1     = SLAM0(L,KP)
            IF (KPARTICLE.EQ.1) THEN
               CTMP=(H1-0.)*(LOG(H1-0.)-1.)+(LAM1-H1)*(LOG(LAM1-H1)-1.)
            ELSE   
               CTMP=(H1-LAM1)*(LOG(H1-LAM1)-1.)
               CTMP=CTMP+(LAM2-H1)*(LOG(LAM2-H1)-1.)
            ENDIF
            CTM2   = CTMP*SF0(L,KP,0,2)/(B0K*SQRT(PI))*H3*H4
            CTM3   = CTMP*SF0(L,KP,0,3)/(B0K*SQRT(PI))*H3*H4
            CTM4   = CTMP*SF0(L,KP,0,4)/(B0K*SQRT(PI))*H3*H4
            CTMP   = CTMP*SF0(L,KP,0,1)/(B0K*SQRT(PI))*H3*H4
            X1PARA(KP,R)=X1PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVX10(M,L,KP)
            X1PERP(KP,R)=X1PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVX10(M,L,KP)
            X1DPHI(KP,R)=X1DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVX10(M,L,KP)
            X2PARA(KP,R)=X2PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVX20(M,L,KP)
            X2PERP(KP,R)=X2PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVX20(M,L,KP)
            X2DPHI(KP,R)=X2DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVX20(M,L,KP)
            Q1PARA(KP,R)=Q1PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVQ10(M,L,KP)
            Q1PERP(KP,R)=Q1PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVQ10(M,L,KP)
            Q1DPHI(KP,R)=Q1DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVQ10(M,L,KP)
            Q2PARA(KP,R)=Q2PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVQ20(M,L,KP)
            Q2PERP(KP,R)=Q2PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVQ20(M,L,KP)
            Q2DPHI(KP,R)=Q2DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVQ20(M,L,KP)
            Q3PARA(KP,R)=Q3PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVQ30(M,L,KP)
            Q3PERP(KP,R)=Q3PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVQ30(M,L,KP)
            Q3DPHI(KP,R)=Q3DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVQ30(M,L,KP)
            DPPARA(KP,R)=DPPARA(KP,R)+CTM2*SVPARA0(K,L,KP)*SVDP0(M,L,KP)
            DPPERP(KP,R)=DPPERP(KP,R)+CTM2*SVPERP0(K,L,KP)*SVDP0(M,L,KP)
            DPDPHI(KP,R)=DPDPHI(KP,R)+CTM4*SVDPHI0(K,L,KP)*SVDP0(M,L,KP)
         ENDIF
         ENDDO
         ENDDO

      CASE (4)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
         IF (ABS(RLM(L)).LT.0.1.AND.KPARTICLE.EQ.0) THEN
            R=2
            CTMP   = VI0(1,KP)*H3
            CTM2   = VI0(2,KP)*H3
            CTM3   = VI0(3,KP)*H3
            CTM4   = VI0(4,KP)*H3
            X1PARA(KP,R)=X1PARA(KP,R)+CTMP*VPARA0(K,L)*VX10(M,L)
            X1PERP(KP,R)=X1PERP(KP,R)+CTMP*VPERP0(K,L)*VX10(M,L)
            X1DPHI(KP,R)=X1DPHI(KP,R)+CTM3*VDPHI0(K,L)*VX10(M,L)
            X2PARA(KP,R)=X2PARA(KP,R)+CTMP*VPARA0(K,L)*VX20(M,L)
            X2PERP(KP,R)=X2PERP(KP,R)+CTMP*VPERP0(K,L)*VX20(M,L)
            X2DPHI(KP,R)=X2DPHI(KP,R)+CTM3*VDPHI0(K,L)*VX20(M,L)
            Q1PARA(KP,R)=Q1PARA(KP,R)+CTMP*VPARA0(K,L)*VQ10(M,L)
            Q1PERP(KP,R)=Q1PERP(KP,R)+CTMP*VPERP0(K,L)*VQ10(M,L)
            Q1DPHI(KP,R)=Q1DPHI(KP,R)+CTM3*VDPHI0(K,L)*VQ10(M,L)
            Q2PARA(KP,R)=Q2PARA(KP,R)+CTMP*VPARA0(K,L)*VQ20(M,L)
            Q2PERP(KP,R)=Q2PERP(KP,R)+CTMP*VPERP0(K,L)*VQ20(M,L)
            Q2DPHI(KP,R)=Q2DPHI(KP,R)+CTM3*VDPHI0(K,L)*VQ20(M,L)
            Q3PARA(KP,R)=Q3PARA(KP,R)+CTMP*VPARA0(K,L)*VQ30(M,L)
            Q3PERP(KP,R)=Q3PERP(KP,R)+CTMP*VPERP0(K,L)*VQ30(M,L)
            Q3DPHI(KP,R)=Q3DPHI(KP,R)+CTM3*VDPHI0(K,L)*VQ30(M,L)
            DPPARA(KP,R)=DPPARA(KP,R)+CTM2*VPARA0(K,L)*VDP0(M,L)
            DPPERP(KP,R)=DPPERP(KP,R)+CTM2*VPERP0(K,L)*VDP0(M,L)
            DPDPHI(KP,R)=DPDPHI(KP,R)+CTM4*VDPHI0(K,L)*VDP0(M,L)
         ENDIF
         ENDDO
         ENDDO

      CASE (5)   
         R=1
         DO KP=1,NSPECIES
         L = M-K + (M2-M1) + 1
         X1PARA(KP,R) = ZKIA(1,KP)*ZGL0PA(L)*LAMH*H3
         X1PERP(KP,R) = ZKIA(1,KP)*ZGL0PE(L)*LAMH*H3
         X1DPHI(KP,R) = ZKIA(4,KP)*ZGL0DP(L)*LAMH*H3
         Q1PARA(KP,R) = ZKIA(2,KP)*ZGL1PA(L)*LAMH*H3
         Q1PERP(KP,R) = ZKIA(2,KP)*ZGL1PE(L)*LAMH*H3
         Q1DPHI(KP,R) = ZKIA(5,KP)*ZGL1DP(L)*LAMH*H3
         Q2PARA(KP,R) = ZKIA(2,KP)*ZGL2PA(L)*LAMH*H3
         Q2PERP(KP,R) = ZKIA(2,KP)*ZGL2PE(L)*LAMH*H3
         Q2DPHI(KP,R) = ZKIA(5,KP)*ZGL2DP(L)*LAMH*H3
         Q3PARA(KP,R) = ZKIA(2,KP)*ZGL3PA(L)*LAMH*H3
         Q3PERP(KP,R) = ZKIA(2,KP)*ZGL3PE(L)*LAMH*H3
         Q3DPHI(KP,R) = ZKIA(5,KP)*ZGL3DP(L)*LAMH*H3
         DPPARA(KP,R) = ZKIA(3,KP)*ZGL0PA(L)*LAMH*H3
         DPPERP(KP,R) = ZKIA(3,KP)*ZGL0PE(L)*LAMH*H3
         DPDPHI(KP,R) = ZKIA(6,KP)*ZGL0DP(L)*LAMH*H3
         ENDDO
      CASE (6)
         R=1
         DO KP=1,NSPECIES
         IF ((ISPECIES_F0(KP).EQ.1.OR.ISPECIES_F0(KP).EQ.2)
     &       .AND.ABS(PSPECIES_AT(KP)).GT.0.) THEN
         H1 = DEPSALPHADPSI(JS,KGRID,KP)
         IF (KGRID.EQ.1) THEN
            L = K-M 
            IF(L.LT.0.)  X1PARA(KP,R) =-ZC1(JS,1,KP)
     &                                *H1*CONJG(JACOBI(JS,-L+1))*H3
            IF(L.EQ.0)   X1PARA(KP,R) =-ZC1(JS,1,KP)*H1*JACOBI(JS,1)*H3
            IF(L.GT.0.1) X1PARA(KP,R) =-ZC1(JS,1,KP)*H1*JACOBI(JS,L+1)
     &                                  *H3
            X1PERP(KP,R) = X1PARA(KP,R)
         ELSEIF (KGRID.EQ.2) THEN
            L = K - M
            IF(L.LT.0.)  X1PARA(KP,R) =-ZC1(JS,2,KP)*H1
     &                                *CONJG(JACOBM(JS,-L+1))*H3
            IF(L.EQ.0.)  X1PARA(KP,R) =-ZC1(JS,2,KP)*H1*JACOBM(JS,1)*H3
            IF(L.GT.0.1) X1PARA(KP,R) =-ZC1(JS,2,KP)*H1*JACOBM(JS,L+1)
     &                                  *H3
            X1PERP(KP,R) = X1PARA(KP,R)
         ENDIF
         ENDIF
         ENDDO
      END SELECT   

      IF (KGRID.EQ.1) THEN
         VX1PARA(K,M,JS_MAT)=VX1PARA(K,M,JS_MAT)+SUM(SUM(X1PARA,1),1)
         VX1PERP(K,M,JS_MAT)=VX1PERP(K,M,JS_MAT)+SUM(SUM(X1PERP,1),1)
         VX2PARA(K,M,JS_MAT)=VX2PARA(K,M,JS_MAT)+SUM(SUM(X2PARA,1),1)
         VX2PERP(K,M,JS_MAT)=VX2PERP(K,M,JS_MAT)+SUM(SUM(X2PERP,1),1) 
         VQ1PARA(K,M,JS_MAT)=VQ1PARA(K,M,JS_MAT)+SUM(SUM(Q1PARA,1),1)
         VQ1PERP(K,M,JS_MAT)=VQ1PERP(K,M,JS_MAT)+SUM(SUM(Q1PERP,1),1)
         VQ2PARA(K,M,JS_MAT)=VQ2PARA(K,M,JS_MAT)+SUM(SUM(Q2PARA,1),1) 
         VQ2PERP(K,M,JS_MAT)=VQ2PERP(K,M,JS_MAT)+SUM(SUM(Q2PERP,1),1) 
         VQ3PARA(K,M,JS_MAT)=VQ3PARA(K,M,JS_MAT)+SUM(SUM(Q3PARA,1),1) 
         VQ3PERP(K,M,JS_MAT)=VQ3PERP(K,M,JS_MAT)+SUM(SUM(Q3PERP,1),1)
         VDPPARA(K,M,JS_MAT)=VDPPARA(K,M,JS_MAT)+SUM(SUM(DPPARA,1),1) 
         VDPPERP(K,M,JS_MAT)=VDPPERP(K,M,JS_MAT)+SUM(SUM(DPPERP,1),1)

         IF (INCDPHI.GT.0) THEN
         VX1DPHI(K,M,JS_MAT)=VX1DPHI(K,M,JS_MAT)+SUM(SUM(X1DPHI,1),1)
         VX2DPHI(K,M,JS_MAT)=VX2DPHI(K,M,JS_MAT)+SUM(SUM(X2DPHI,1),1)
         VQ1DPHI(K,M,JS_MAT)=VQ1DPHI(K,M,JS_MAT)+SUM(SUM(Q1DPHI,1),1)
         VQ2DPHI(K,M,JS_MAT)=VQ2DPHI(K,M,JS_MAT)+SUM(SUM(Q2DPHI,1),1)
         VQ3DPHI(K,M,JS_MAT)=VQ3DPHI(K,M,JS_MAT)+SUM(SUM(Q3DPHI,1),1)
         VDPDPHI(K,M,JS_MAT)=VDPDPHI(K,M,JS_MAT)+SUM(SUM(DPDPHI,1),1)
         ENDIF

         IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &       .AND.INCKIN.GT.0) THEN
         VX1PARE(K,M,JS_MAT)=VX1PARE(K,M,JS_MAT)+SUM(X1PARA(2,:))
         VX1PERE(K,M,JS_MAT)=VX1PERE(K,M,JS_MAT)+SUM(X1PERP(2,:))
         VX2PARE(K,M,JS_MAT)=VX2PARE(K,M,JS_MAT)+SUM(X2PARA(2,:))
         VX2PERE(K,M,JS_MAT)=VX2PERE(K,M,JS_MAT)+SUM(X2PERP(2,:)) 
         VQ1PARE(K,M,JS_MAT)=VQ1PARE(K,M,JS_MAT)+SUM(Q1PARA(2,:))
         VQ1PERE(K,M,JS_MAT)=VQ1PERE(K,M,JS_MAT)+SUM(Q1PERP(2,:))
         VQ2PARE(K,M,JS_MAT)=VQ2PARE(K,M,JS_MAT)+SUM(Q2PARA(2,:)) 
         VQ2PERE(K,M,JS_MAT)=VQ2PERE(K,M,JS_MAT)+SUM(Q2PERP(2,:)) 
         VQ3PARE(K,M,JS_MAT)=VQ3PARE(K,M,JS_MAT)+SUM(Q3PARA(2,:)) 
         VQ3PERE(K,M,JS_MAT)=VQ3PERE(K,M,JS_MAT)+SUM(Q3PERP(2,:))
         VDPPARE(K,M,JS_MAT)=VDPPARE(K,M,JS_MAT)+SUM(DPPARA(2,:)) 
         VDPPERE(K,M,JS_MAT)=VDPPERE(K,M,JS_MAT)+SUM(DPPERP(2,:))
         ENDIF
      ELSE
         VX1PARAM(K,M,JS_MAT)=VX1PARAM(K,M,JS_MAT)+SUM(SUM(X1PARA,1),1)
         VX1PERPM(K,M,JS_MAT)=VX1PERPM(K,M,JS_MAT)+SUM(SUM(X1PERP,1),1)
         VX2PARAM(K,M,JS_MAT)=VX2PARAM(K,M,JS_MAT)+SUM(SUM(X2PARA,1),1)
         VX2PERPM(K,M,JS_MAT)=VX2PERPM(K,M,JS_MAT)+SUM(SUM(X2PERP,1),1) 
         VQ1PARAM(K,M,JS_MAT)=VQ1PARAM(K,M,JS_MAT)+SUM(SUM(Q1PARA,1),1)
         VQ1PERPM(K,M,JS_MAT)=VQ1PERPM(K,M,JS_MAT)+SUM(SUM(Q1PERP,1),1)
         VQ2PARAM(K,M,JS_MAT)=VQ2PARAM(K,M,JS_MAT)+SUM(SUM(Q2PARA,1),1) 
         VQ2PERPM(K,M,JS_MAT)=VQ2PERPM(K,M,JS_MAT)+SUM(SUM(Q2PERP,1),1) 
         VQ3PARAM(K,M,JS_MAT)=VQ3PARAM(K,M,JS_MAT)+SUM(SUM(Q3PARA,1),1) 
         VQ3PERPM(K,M,JS_MAT)=VQ3PERPM(K,M,JS_MAT)+SUM(SUM(Q3PERP,1),1)
         VDPPARAM(K,M,JS_MAT)=VDPPARAM(K,M,JS_MAT)+SUM(SUM(DPPARA,1),1) 
         VDPPERPM(K,M,JS_MAT)=VDPPERPM(K,M,JS_MAT)+SUM(SUM(DPPERP,1),1)

         IF (INCDPHI.GT.0) THEN
         VX1DPHIM(K,M,JS_MAT)=VX1DPHIM(K,M,JS_MAT)+SUM(SUM(X1DPHI,1),1)
         VX2DPHIM(K,M,JS_MAT)=VX2DPHIM(K,M,JS_MAT)+SUM(SUM(X2DPHI,1),1)
         VQ1DPHIM(K,M,JS_MAT)=VQ1DPHIM(K,M,JS_MAT)+SUM(SUM(Q1DPHI,1),1)
         VQ2DPHIM(K,M,JS_MAT)=VQ2DPHIM(K,M,JS_MAT)+SUM(SUM(Q2DPHI,1),1)
         VQ3DPHIM(K,M,JS_MAT)=VQ3DPHIM(K,M,JS_MAT)+SUM(SUM(Q3DPHI,1),1)
         VDPDPHIM(K,M,JS_MAT)=VDPDPHIM(K,M,JS_MAT)+SUM(SUM(DPDPHI,1),1)
         ENDIF

         IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &       .AND.INCKIN.GT.0) THEN
         VX1PAREM(K,M,JS_MAT)=VX1PAREM(K,M,JS_MAT)+SUM(X1PARA(2,:))
         VX1PEREM(K,M,JS_MAT)=VX1PEREM(K,M,JS_MAT)+SUM(X1PERP(2,:))
         VX2PAREM(K,M,JS_MAT)=VX2PAREM(K,M,JS_MAT)+SUM(X2PARA(2,:))
         VX2PEREM(K,M,JS_MAT)=VX2PEREM(K,M,JS_MAT)+SUM(X2PERP(2,:)) 
         VQ1PAREM(K,M,JS_MAT)=VQ1PAREM(K,M,JS_MAT)+SUM(Q1PARA(2,:))
         VQ1PEREM(K,M,JS_MAT)=VQ1PEREM(K,M,JS_MAT)+SUM(Q1PERP(2,:))
         VQ2PAREM(K,M,JS_MAT)=VQ2PAREM(K,M,JS_MAT)+SUM(Q2PARA(2,:)) 
         VQ2PEREM(K,M,JS_MAT)=VQ2PEREM(K,M,JS_MAT)+SUM(Q2PERP(2,:)) 
         VQ3PAREM(K,M,JS_MAT)=VQ3PAREM(K,M,JS_MAT)+SUM(Q3PARA(2,:)) 
         VQ3PEREM(K,M,JS_MAT)=VQ3PEREM(K,M,JS_MAT)+SUM(Q3PERP(2,:))
         VDPPAREM(K,M,JS_MAT)=VDPPAREM(K,M,JS_MAT)+SUM(DPPARA(2,:)) 
         VDPPEREM(K,M,JS_MAT)=VDPPEREM(K,M,JS_MAT)+SUM(DPPERP(2,:))
         ENDIF
      ENDIF

C     SET ENERGY COMPONENT MATRIX         
      DO KP=1,NSPECIES
      IF (KPARTICLE.EQ.0) THEN
         R=1
C        TRAPPED PARTICLES
         IF (ICASE.EQ.5.OR.ICASE.EQ.6) THEN
C        ADIABATIC COMPONENT
         CALL SETDWKCOMPMAT ( KP,2,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))
         ENDIF         
         IF (ICASE.EQ.1.OR.ICASE.EQ.2.OR.ICASE.EQ.3) THEN
C        BOUNCE COMPONENT             
         CALL SETDWKCOMPMAT ( KP,4,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))
         ENDIF
         IF (ICASE.EQ.1.OR.ICASE.EQ.2.OR.ICASE.EQ.3.OR.ICASE.EQ.4) THEN
C        PRECESSION COMPONENT     
         R=2
         CALL SETDWKCOMPMAT ( KP,5,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))
         ENDIF
      ELSE
C        PASSING PARTICLES
         R=1
         IF (ICASE.EQ.5.OR.ICASE.EQ.6) THEN
C        ADIABATIC COMPONENT
         CALL SETDWKCOMPMAT ( KP,1,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))               
         ENDIF
         IF (ICASE.EQ.1.OR.ICASE.EQ.2.OR.ICASE.EQ.3) THEN
C         TRANSIT COMPONENT
         CALL SETDWKCOMPMAT ( KP,3,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))         
         ENDIF         
      ENDIF
      ENDDO
      
      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         IF (K.EQ.MIN(2-M1+1,M2-M1+1).AND.M.EQ.K) THEN
         IF (ABS(LAMM(2)-LAM).LT.1.0E-13) 
     &      WRITE(*,*) 'KJPFILL: ICASE LAM X1PARA X1PERP'
         WRITE(*,121)ICASE,LAM,X1PARA,X1PERP
         ENDIF
      ENDIF
 121  FORMAT(I2,1X,5(E13.5))

      ENDDO
      ENDDO
      
      DEALLOCATE (X1PARA,X2PARA,Q1PARA,Q2PARA,Q3PARA,DPPARA,
     &            X1PERP,X2PERP,Q1PERP,Q2PERP,Q3PERP,DPPERP,
     &            X1DPHI,X2DPHI,Q1DPHI,Q2DPHI,Q3DPHI,DPDPHI)
      RETURN
      END

C=======================================================================
C PUT ALL I,H,G FACTORS TOGETHER AND FILL INTO MATRIX ELEMENTS 
C INCLUDING THE FIRST ORDER FOW CORRECTIONS
C ICASE = 1: REGULAR NON-ADIABATIC CONTRIBUTION (BOTH PASS&TRAP)
C         2: PITCH ANGLE DEPEDENT SINGULAR PART (BOTH PASS&TRAP)
C         3: PITCH ANGLE INDEPEDNET SINGULAR PART (BOTH PASS&TRAP)
C         4: SPECIAL SINGULAR CONTRIBUTION FROM RLM(L)=0 FOR 
C            PRECESSIONAL DRIFT RESONANCE OF TRAPPED PARTICLES
C         5: ADIABATIC CONTRIBUTION FROM ANISOTROPIC EP (BOTH PASS&TRAP)
C NOTE THAN ICASE=5 DOES NOT APPLY TO X1PARE ETC.      
C=======================================================================
      SUBROUTINE KJPFILL1(JS,JS_MAT,KGRID,LAMH,KPARTICLE,ICASE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      IMPLICIT NONE

      INTEGER    JS,JS_MAT,KGRID,KPARTICLE,ICASE,R,
     &           K,M,KP,J,L,N,M00,J00,J2M,J1M,J1P,J2P
      REAL*8     LAMH,LAM1,LAM2,H1,H2,H3,H4
      COMPLEX*16 CTMP,CTM2,CTM3,CTM4,CTMP1
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE:: 
     &           X1PARA,X2PARA,Q1PARA,Q2PARA,Q3PARA,DPPARA,
     &           X1PERP,X2PERP,Q1PERP,Q2PERP,Q3PERP,DPPERP,
     &           X1DPHI,X2DPHI,Q1DPHI,Q2DPHI,Q3DPHI,DPDPHI
      INTEGER KCHECK
      KCHECK=0

      LAM1   = HKMIN(JS,KGRID)
      LAM2   = HKMAX(JS,KGRID)
      H3     = WFUN (JS,KGRID)
      ALLOCATE ( X1PARA(NSPECIES,2),X2PARA(NSPECIES,2),
     &           Q1PARA(NSPECIES,2),Q2PARA(NSPECIES,2),
     &           Q3PARA(NSPECIES,2),X1PERP(NSPECIES,2),
     &           X2PERP(NSPECIES,2),Q1PERP(NSPECIES,2),
     &           Q2PERP(NSPECIES,2),Q3PERP(NSPECIES,2),
     &           DPPARA(NSPECIES,2),DPPERP(NSPECIES,2), 
     &           X1DPHI(NSPECIES,2),X2DPHI(NSPECIES,2),
     &           Q1DPHI(NSPECIES,2),Q2DPHI(NSPECIES,2),
     &           Q3DPHI(NSPECIES,2),DPDPHI(NSPECIES,2) )
      DO K=1,MSMAX
      DO M=1,MSMAX

      X1PARA = 0.
      X1PERP = 0.
      X1DPHI = 0.
      X2PARA = 0.
      X2PERP = 0.
      X2DPHI = 0.
      Q1PARA = 0.
      Q1PERP = 0.
      Q1DPHI = 0.
      Q2PARA = 0.
      Q2PERP = 0.
      Q2DPHI = 0.
      Q3PARA = 0.
      Q3PERP = 0.
      Q3DPHI = 0.
      DPPARA = 0.
      DPPERP = 0.
      DPDPHI = 0.
      
      SELECT CASE (ICASE)     
      CASE (1)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
            IF (KPARTICLE.EQ.0 .AND. ABS(RLM(L)).LT.0.1) THEN
               R=2
            ELSE
               R=1
            ENDIF
            IF (KPARTICLE.EQ.1) THEN
               CTMP = (VI(1,L,KP)+VI1(1,L,KP))*LAMH*H3
               CTM2 = (VI(2,L,KP)+VI1(2,L,KP))*LAMH*H3
               CTM3 = (VI(3,L,KP)+VI1(3,L,KP))*LAMH*H3
               CTM4 = (VI(4,L,KP)+VI1(4,L,KP))*LAMH*H3
            ENDIF
            IF (KPARTICLE.EQ.0) THEN
               CTMP = VI(1,L,KP)*LAMH*H3
               CTM2 = VI(2,L,KP)*LAMH*H3
               CTM3 = VI(3,L,KP)*LAMH*H3
               CTM4 = VI(4,L,KP)*LAMH*H3
            ENDIF

            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA(K,L)*VX1(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP(K,L)*VX1(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI(K,L)*VX1(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA(K,L)*VX2(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP(K,L)*VX2(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI(K,L)*VX2(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA(K,L)*VQ1(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP(K,L)*VQ1(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI(K,L)*VQ1(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA(K,L)*VQ2(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP(K,L)*VQ2(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI(K,L)*VQ2(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA(K,L)*VQ3(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP(K,L)*VQ3(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI(K,L)*VQ3(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA(K,L)*VDP(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP(K,L)*VDP(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI(K,L)*VDP(M,L)

            CTMP   = VI2(1,L,KP)*LAMH*H3
            CTM2   = VI2(2,L,KP)*LAMH*H3
            CTM3   = VI2(3,L,KP)*LAMH*H3
            CTM4   = VI2(4,L,KP)*LAMH*H3
            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA1(K,L)*VX1(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP1(K,L)*VX1(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI1(K,L)*VX1(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA1(K,L)*VX2(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP1(K,L)*VX2(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI1(K,L)*VX2(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA1(K,L)*VQ1(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP1(K,L)*VQ1(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI1(K,L)*VQ1(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA1(K,L)*VQ2(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP1(K,L)*VQ2(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI1(K,L)*VQ2(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA1(K,L)*VQ3(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP1(K,L)*VQ3(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI1(K,L)*VQ3(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA1(K,L)*VDP(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP1(K,L)*VDP(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI1(K,L)*VDP(M,L)

            CTMP   = VI3(1,L,KP)*LAMH*H3
            CTM2   = VI3(2,L,KP)*LAMH*H3
            CTM3   = VI3(3,L,KP)*LAMH*H3
            CTM4   = VI3(4,L,KP)*LAMH*H3
            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VX11(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VX11(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VX11(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VX21(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VX21(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VX21(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VQ11(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VQ11(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VQ11(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VQ21(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VQ21(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VQ21(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VQ31(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VQ31(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VQ31(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   + CTM2*VPARA(K,L)*VDP1(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   + CTM2*VPERP(K,L)*VDP1(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   + CTM4*VDPHI(K,L)*VDP1(M,L)*VDPLNP(M)

            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VX1(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VX1(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VX1(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VX2(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VX2(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VX2(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VQ1(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VQ1(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VQ1(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VQ2(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VQ2(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VQ2(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VQ3(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VQ3(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VQ3(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   - CTM2*VPARA1(K,L)*VDP(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   - CTM2*VPERP1(K,L)*VDP(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   - CTM4*VDPHI1(K,L)*VDP(M,L)*VDPLNP(M)

            IF (ABS(RLM(L)).LT.0.1.AND.SLAMD0.GT.0.0.AND.KPARTICLE.EQ.0)
     &      THEN
            CTMP   =-VI(1,L,KP)*LAMH*H3
            CTM2   =-VI(2,L,KP)*LAMH*H3
            CTM3   =-VI(3,L,KP)*LAMH*H3
            CTM4   =-VI(4,L,KP)*LAMH*H3
            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA0(K,L)*VX10(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP0(K,L)*VX10(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI0(K,L)*VX10(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA0(K,L)*VX20(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP0(K,L)*VX20(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI0(K,L)*VX20(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA0(K,L)*VQ10(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP0(K,L)*VQ10(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI0(K,L)*VQ10(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA0(K,L)*VQ20(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP0(K,L)*VQ20(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI0(K,L)*VQ20(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA0(K,L)*VQ30(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP0(K,L)*VQ30(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI0(K,L)*VQ30(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA0(K,L)*VDP0(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP0(K,L)*VDP0(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI0(K,L)*VDP0(M,L)

            CTMP   =-VI2(1,L,KP)*LAMH*H3
            CTM2   =-VI2(2,L,KP)*LAMH*H3
            CTM3   =-VI2(3,L,KP)*LAMH*H3
            CTM4   =-VI2(4,L,KP)*LAMH*H3
            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA01(K,L)*VX10(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP01(K,L)*VX10(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI01(K,L)*VX10(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA01(K,L)*VX20(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP01(K,L)*VX20(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI01(K,L)*VX20(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA01(K,L)*VQ10(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP01(K,L)*VQ10(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI01(K,L)*VQ10(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA01(K,L)*VQ20(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP01(K,L)*VQ20(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI01(K,L)*VQ20(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA01(K,L)*VQ30(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP01(K,L)*VQ30(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI01(K,L)*VQ30(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA01(K,L)*VDP0(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP01(K,L)*VDP0(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI01(K,L)*VDP0(M,L)

            CTMP   =-VI3(1,L,KP)*LAMH*H3
            CTM2   =-VI3(2,L,KP)*LAMH*H3
            CTM3   =-VI3(3,L,KP)*LAMH*H3
            CTM4   =-VI3(4,L,KP)*LAMH*H3
            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VX101(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VX101(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VX101(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VX201(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VX201(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VX201(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VQ101(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VQ101(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VQ101(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VQ201(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VQ201(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VQ201(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VQ301(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VQ301(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VQ301(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   + CTM2*VPARA0(K,L)*VDP01(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   + CTM2*VPERP0(K,L)*VDP01(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   + CTM4*VDPHI0(K,L)*VDP01(M,L)*VDPLNP(M)

            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VX10(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VX10(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VX10(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VX20(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VX20(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VX20(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VQ10(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VQ10(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VQ10(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VQ20(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VQ20(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VQ20(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VQ30(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VQ30(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VQ30(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   - CTM2*VPARA01(K,L)*VDP0(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   - CTM2*VPERP01(K,L)*VDP0(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   - CTM4*VDPHI01(K,L)*VDP0(M,L)*VDPLNP(M)
            ENDIF
         ENDDO
         ENDDO
      CASE (2)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
         IF (SLAM0(L,KP).GT.0.) THEN
             IF (KPARTICLE.EQ.0 .AND. ABS(RLM(L)).LT.0.1) THEN
               R=2
            ELSE
               R=1
            ENDIF         
            IF (KPARTICLE.EQ.1) THEN
               H4 = PSPECIES_NP(KP)
            ELSEIF (ABS(RLM(L)).GT.0.1) THEN
               H4 = PSPECIES_NTB(KP)
            ELSE
               H4 = PSPECIES_NTD(KP)
            ENDIF
            H1     = SLAM0(L,KP)
            IF (KPARTICLE.EQ.1) THEN
               H2 = SF0(L,KP,0,1)+SF0(L,KP,1,1)
               CTMP   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
               H2 = SF0(L,KP,0,2)+SF0(L,KP,1,2)
               CTM2   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
               H2 = SF0(L,KP,0,3)+SF0(L,KP,1,3)
               CTM3   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
               H2 = SF0(L,KP,0,4)+SF0(L,KP,1,4)
               CTM4   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
            ENDIF
            IF (KPARTICLE.EQ.0) THEN
               H2 = SF0(L,KP,0,1)
               CTMP   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
               H2 = SF0(L,KP,0,2)
               CTM2   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
               H2 = SF0(L,KP,0,3)
               CTM3   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
               H2 = SF0(L,KP,0,4)
               CTM4   = LOG(ABS(LAM-H1))*H2*LAMH*H3*H4
            ENDIF

            X1PARA(KP,R)=X1PARA(KP,R)+CTMP*VPARA(K,L)*VX1(M,L)
            X1PERP(KP,R)=X1PERP(KP,R)+CTMP*VPERP(K,L)*VX1(M,L)
            X1DPHI(KP,R)=X1DPHI(KP,R)+CTM3*VDPHI(K,L)*VX1(M,L)
            X2PARA(KP,R)=X2PARA(KP,R)+CTMP*VPARA(K,L)*VX2(M,L)
            X2PERP(KP,R)=X2PERP(KP,R)+CTMP*VPERP(K,L)*VX2(M,L)
            X2DPHI(KP,R)=X2DPHI(KP,R)+CTM3*VDPHI(K,L)*VX2(M,L)
            Q1PARA(KP,R)=Q1PARA(KP,R)+CTMP*VPARA(K,L)*VQ1(M,L)
            Q1PERP(KP,R)=Q1PERP(KP,R)+CTMP*VPERP(K,L)*VQ1(M,L)
            Q1DPHI(KP,R)=Q1DPHI(KP,R)+CTM3*VDPHI(K,L)*VQ1(M,L)
            Q2PARA(KP,R)=Q2PARA(KP,R)+CTMP*VPARA(K,L)*VQ2(M,L)
            Q2PERP(KP,R)=Q2PERP(KP,R)+CTMP*VPERP(K,L)*VQ2(M,L)
            Q2DPHI(KP,R)=Q2DPHI(KP,R)+CTM3*VDPHI(K,L)*VQ2(M,L)
            Q3PARA(KP,R)=Q3PARA(KP,R)+CTMP*VPARA(K,L)*VQ3(M,L)
            Q3PERP(KP,R)=Q3PERP(KP,R)+CTMP*VPERP(K,L)*VQ3(M,L)
            Q3DPHI(KP,R)=Q3DPHI(KP,R)+CTM3*VDPHI(K,L)*VQ3(M,L)
            DPPARA(KP,R)=DPPARA(KP,R)+CTM2*VPARA(K,L)*VDP(M,L)
            DPPERP(KP,R)=DPPERP(KP,R)+CTM2*VPERP(K,L)*VDP(M,L)
            DPDPHI(KP,R)=DPDPHI(KP,R)+CTM4*VDPHI(K,L)*VDP(M,L)

            X1PARA(KP,R)=X1PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVX10(M,L,KP)
            X1PERP(KP,R)=X1PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVX10(M,L,KP)
            X1DPHI(KP,R)=X1DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVX10(M,L,KP)
            X2PARA(KP,R)=X2PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVX20(M,L,KP)
            X2PERP(KP,R)=X2PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVX20(M,L,KP)
            X2DPHI(KP,R)=X2DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVX20(M,L,KP)
            Q1PARA(KP,R)=Q1PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVQ10(M,L,KP)
            Q1PERP(KP,R)=Q1PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVQ10(M,L,KP)
            Q1DPHI(KP,R)=Q1DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVQ10(M,L,KP)
            Q2PARA(KP,R)=Q2PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVQ20(M,L,KP)
            Q2PERP(KP,R)=Q2PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVQ20(M,L,KP)
            Q2DPHI(KP,R)=Q2DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVQ20(M,L,KP)
            Q3PARA(KP,R)=Q3PARA(KP,R)-CTMP*SVPARA0(K,L,KP)*SVQ30(M,L,KP)
            Q3PERP(KP,R)=Q3PERP(KP,R)-CTMP*SVPERP0(K,L,KP)*SVQ30(M,L,KP)
            Q3DPHI(KP,R)=Q3DPHI(KP,R)-CTM3*SVDPHI0(K,L,KP)*SVQ30(M,L,KP)
            DPPARA(KP,R)=DPPARA(KP,R)-CTM2*SVPARA0(K,L,KP)*SVDP0(M,L,KP)
            DPPERP(KP,R)=DPPERP(KP,R)-CTM2*SVPERP0(K,L,KP)*SVDP0(M,L,KP)
            DPDPHI(KP,R)=DPDPHI(KP,R)-CTM4*SVDPHI0(K,L,KP)*SVDP0(M,L,KP)

            CTMP   = LOG(ABS(LAM-H1))*SF0(L,KP,2,1)*LAMH*H3*H4
            CTM2   = LOG(ABS(LAM-H1))*SF0(L,KP,2,2)*LAMH*H3*H4
            CTM3   = LOG(ABS(LAM-H1))*SF0(L,KP,2,3)*LAMH*H3*H4
            CTM4   = LOG(ABS(LAM-H1))*SF0(L,KP,2,4)*LAMH*H3*H4
            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA1(K,L)*VX1(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP1(K,L)*VX1(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI1(K,L)*VX1(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA1(K,L)*VX2(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP1(K,L)*VX2(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI1(K,L)*VX2(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA1(K,L)*VQ1(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP1(K,L)*VQ1(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI1(K,L)*VQ1(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA1(K,L)*VQ2(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP1(K,L)*VQ2(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI1(K,L)*VQ2(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA1(K,L)*VQ3(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP1(K,L)*VQ3(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI1(K,L)*VQ3(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA1(K,L)*VDP(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP1(K,L)*VDP(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI1(K,L)*VDP(M,L)

            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   - CTMP*SVPARA01(K,L,KP)*SVX10(M,L,KP)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   - CTMP*SVPERP01(K,L,KP)*SVX10(M,L,KP)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   - CTM3*SVDPHI01(K,L,KP)*SVX10(M,L,KP)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   - CTMP*SVPARA01(K,L,KP)*SVX20(M,L,KP)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   - CTMP*SVPERP01(K,L,KP)*SVX20(M,L,KP)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   - CTM3*SVDPHI01(K,L,KP)*SVX20(M,L,KP)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   - CTMP*SVPARA01(K,L,KP)*SVQ10(M,L,KP)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                    - CTMP*SVPERP01(K,L,KP)*SVQ10(M,L,KP)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                    - CTM3*SVDPHI01(K,L,KP)*SVQ10(M,L,KP)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   - CTMP*SVPARA01(K,L,KP)*SVQ20(M,L,KP)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   - CTMP*SVPERP01(K,L,KP)*SVQ20(M,L,KP)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   - CTM3*SVDPHI01(K,L,KP)*SVQ20(M,L,KP)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   - CTMP*SVPARA01(K,L,KP)*SVQ30(M,L,KP)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   - CTMP*SVPERP01(K,L,KP)*SVQ30(M,L,KP)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   - CTM3*SVDPHI01(K,L,KP)*SVQ30(M,L,KP)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   - CTM2*SVPARA01(K,L,KP)*SVDP0(M,L,KP)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   - CTM2*SVPERP01(K,L,KP)*SVDP0(M,L,KP)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   - CTM4*SVDPHI01(K,L,KP)*SVDP0(M,L,KP)

            CTMP   = LOG(ABS(LAM-H1))*SF0(L,KP,3,1)*LAMH*H3*H4
            CTM2   = LOG(ABS(LAM-H1))*SF0(L,KP,3,2)*LAMH*H3*H4
            CTM3   = LOG(ABS(LAM-H1))*SF0(L,KP,3,3)*LAMH*H3*H4
            CTM4   = LOG(ABS(LAM-H1))*SF0(L,KP,3,4)*LAMH*H3*H4
            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VX11(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VX11(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VX11(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VX21(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VX21(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VX21(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VQ11(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VQ11(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VQ11(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VQ21(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VQ21(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VQ21(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   + CTMP*VPARA(K,L)*VQ31(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   + CTMP*VPERP(K,L)*VQ31(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   + CTM3*VDPHI(K,L)*VQ31(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   + CTM2*VPARA(K,L)*VDP1(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   + CTM2*VPERP(K,L)*VDP1(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   + CTM4*VDPHI(K,L)*VDP1(M,L)*VDPLNP(M)

            X1PARA(KP,R) =X1PARA(KP,R) 
     &                   -CTMP*SVPARA0(K,L,KP)*SVX101(M,L,KP)*VX1LNP(M)
            X1PERP(KP,R) =X1PERP(KP,R) 
     &                   -CTMP*SVPERP0(K,L,KP)*SVX101(M,L,KP)*VX1LNP(M)
            X1DPHI(KP,R) =X1DPHI(KP,R) 
     &                   -CTM3*SVDPHI0(K,L,KP)*SVX101(M,L,KP)*VX1LNP(M)
            X2PARA(KP,R) =X2PARA(KP,R) 
     &                   -CTMP*SVPARA0(K,L,KP)*SVX201(M,L,KP)*VX2LNP(M)
            X2PERP(KP,R) =X2PERP(KP,R) 
     &                   -CTMP*SVPERP0(K,L,KP)*SVX201(M,L,KP)*VX2LNP(M)
            X2DPHI(KP,R) =X2DPHI(KP,R) 
     &                   -CTM3*SVDPHI0(K,L,KP)*SVX201(M,L,KP)*VX2LNP(M)
            Q1PARA(KP,R) =Q1PARA(KP,R) 
     &                   -CTMP*SVPARA0(K,L,KP)*SVQ101(M,L,KP)*VQ1LNP(M)
            Q1PERP(KP,R) =Q1PERP(KP,R) 
     &                   -CTMP*SVPERP0(K,L,KP)*SVQ101(M,L,KP)*VQ1LNP(M)
            Q1DPHI(KP,R) =Q1DPHI(KP,R) 
     &                   -CTM3*SVDPHI0(K,L,KP)*SVQ101(M,L,KP)*VQ1LNP(M)
            Q2PARA(KP,R) =Q2PARA(KP,R) 
     &                   -CTMP*SVPARA0(K,L,KP)*SVQ201(M,L,KP)*VQ2LNP(M)
            Q2PERP(KP,R) =Q2PERP(KP,R) 
     &                   -CTMP*SVPERP0(K,L,KP)*SVQ201(M,L,KP)*VQ2LNP(M)
            Q2DPHI(KP,R) =Q2DPHI(KP,R) 
     &                   -CTM3*SVDPHI0(K,L,KP)*SVQ201(M,L,KP)*VQ2LNP(M)
            Q3PARA(KP,R) =Q3PARA(KP,R) 
     &                   -CTMP*SVPARA0(K,L,KP)*SVQ301(M,L,KP)*VQ3LNP(M)
            Q3PERP(KP,R) =Q3PERP(KP,R) 
     &                   -CTMP*SVPERP0(K,L,KP)*SVQ301(M,L,KP)*VQ3LNP(M)
            Q3DPHI(KP,R) =Q3DPHI(KP,R) 
     &                   -CTM3*SVDPHI0(K,L,KP)*SVQ301(M,L,KP)*VQ3LNP(M)
            DPPARA(KP,R) =DPPARA(KP,R) 
     &                   -CTM2*SVPARA0(K,L,KP)*SVDP01(M,L,KP)*VDPLNP(M)
            DPPERP(KP,R) =DPPERP(KP,R) 
     &                   -CTM2*SVPERP0(K,L,KP)*SVDP01(M,L,KP)*VDPLNP(M)
            DPDPHI(KP,R) =DPDPHI(KP,R) 
     &                   -CTM4*SVDPHI0(K,L,KP)*SVDP01(M,L,KP)*VDPLNP(M)

            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VX1(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VX1(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VX1(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VX2(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VX2(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VX2(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VQ1(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VQ1(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VQ1(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VQ2(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VQ2(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VQ2(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   - CTMP*VPARA1(K,L)*VQ3(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   - CTMP*VPERP1(K,L)*VQ3(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   - CTM3*VDPHI1(K,L)*VQ3(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   - CTM2*VPARA1(K,L)*VDP(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   - CTM2*VPERP1(K,L)*VDP(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   - CTM4*VDPHI1(K,L)*VDP(M,L)*VDPLNP(M)

            X1PARA(KP,R) =X1PARA(KP,R) 
     &                   +CTMP*SVPARA01(K,L,KP)*SVX10(M,L,KP)*VX1LNP(M)
            X1PERP(KP,R) =X1PERP(KP,R) 
     &                   +CTMP*SVPERP01(K,L,KP)*SVX10(M,L,KP)*VX1LNP(M)
            X1DPHI(KP,R) =X1DPHI(KP,R) 
     &                   +CTM3*SVDPHI01(K,L,KP)*SVX10(M,L,KP)*VX1LNP(M)
            X2PARA(KP,R) =X2PARA(KP,R) 
     &                   +CTMP*SVPARA01(K,L,KP)*SVX20(M,L,KP)*VX2LNP(M)
            X2PERP(KP,R) =X2PERP(KP,R) 
     &                   +CTMP*SVPERP01(K,L,KP)*SVX20(M,L,KP)*VX2LNP(M)
            X2DPHI(KP,R) =X2DPHI(KP,R) 
     &                   +CTM3*SVDPHI01(K,L,KP)*SVX20(M,L,KP)*VX2LNP(M)
            Q1PARA(KP,R) =Q1PARA(KP,R) 
     &                   +CTMP*SVPARA01(K,L,KP)*SVQ10(M,L,KP)*VQ1LNP(M)
            Q1PERP(KP,R) =Q1PERP(KP,R) 
     &                   +CTMP*SVPERP01(K,L,KP)*SVQ10(M,L,KP)*VQ1LNP(M)
            Q1DPHI(KP,R) =Q1DPHI(KP,R) 
     &                   +CTM3*SVDPHI01(K,L,KP)*SVQ10(M,L,KP)*VQ1LNP(M)
            Q2PARA(KP,R) =Q2PARA(KP,R) 
     &                   +CTMP*SVPARA01(K,L,KP)*SVQ20(M,L,KP)*VQ2LNP(M)
            Q2PERP(KP,R) =Q2PERP(KP,R) 
     &                   +CTMP*SVPERP01(K,L,KP)*SVQ20(M,L,KP)*VQ2LNP(M)
            Q2DPHI(KP,R) =Q2DPHI(KP,R) 
     &                   +CTM3*SVDPHI01(K,L,KP)*SVQ20(M,L,KP)*VQ2LNP(M)
            Q3PARA(KP,R) =Q3PARA(KP,R) 
     &                   +CTMP*SVPARA01(K,L,KP)*SVQ30(M,L,KP)*VQ3LNP(M)
            Q3PERP(KP,R) =Q3PERP(KP,R) 
     &                   +CTMP*SVPERP01(K,L,KP)*SVQ30(M,L,KP)*VQ3LNP(M)
            Q3DPHI(KP,R) =Q3DPHI(KP,R) 
     &                   +CTM3*SVDPHI01(K,L,KP)*SVQ30(M,L,KP)*VQ3LNP(M)
            DPPARA(KP,R) =DPPARA(KP,R) 
     &                   +CTM2*SVPARA01(K,L,KP)*SVDP0(M,L,KP)*VDPLNP(M)
            DPPERP(KP,R) =DPPERP(KP,R) 
     &                   +CTM2*SVPERP01(K,L,KP)*SVDP0(M,L,KP)*VDPLNP(M)
            DPDPHI(KP,R) =DPDPHI(KP,R) 
     &                   +CTM4*SVDPHI01(K,L,KP)*SVDP0(M,L,KP)*VDPLNP(M)
         ENDIF
         ENDDO
         ENDDO

      CASE (3)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
         IF (SLAM0(L,KP).GT.0.) THEN
            IF (KPARTICLE.EQ.0 .AND. ABS(RLM(L)).LT.0.1) THEN
               R=2
            ELSE
               R=1
            ENDIF         
            IF (KPARTICLE.EQ.1) THEN
               H4 = PSPECIES_NP(KP)
            ELSEIF (ABS(RLM(L)).GT.0.1) THEN
               H4 = PSPECIES_NTB(KP)
            ELSE
               H4 = PSPECIES_NTD(KP)
            ENDIF
            H1 = SLAM0(L,KP)
            IF (KPARTICLE.EQ.1) THEN
               CTMP1=(H1-0.)*(LOG(H1-0.)-1.)+(LAM1-H1)*(LOG(LAM1-H1)-1.)
               CTMP1=CTMP1/(B0K*SQRT(PI))
               CTMP =CTMP1*(SF0(L,KP,0,1)+SF0(L,KP,1,1))*H3*H4
               CTM2 =CTMP1*(SF0(L,KP,0,2)+SF0(L,KP,1,2))*H3*H4
               CTM3 =CTMP1*(SF0(L,KP,0,3)+SF0(L,KP,1,3))*H3*H4
               CTM4 =CTMP1*(SF0(L,KP,0,4)+SF0(L,KP,1,4))*H3*H4
            ELSE
               CTMP1=(H1-LAM1)*(LOG(H1-LAM1)-1.) 
               CTMP1=CTMP1 + (LAM2-H1)*(LOG(LAM2-H1)-1.)
               CTMP1=CTMP1/(B0K*SQRT(PI))
               CTMP =CTMP1*SF0(L,KP,0,1)*H3*H4
               CTM2 =CTMP1*SF0(L,KP,0,2)*H3*H4
               CTM3 =CTMP1*SF0(L,KP,0,3)*H3*H4
               CTM4 =CTMP1*SF0(L,KP,0,4)*H3*H4
            ENDIF
            X1PARA(KP,R)=X1PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVX10(M,L,KP)
            X1PERP(KP,R)=X1PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVX10(M,L,KP)
            X1DPHI(KP,R)=X1DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVX10(M,L,KP)
            X2PARA(KP,R)=X2PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVX20(M,L,KP)
            X2PERP(KP,R)=X2PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVX20(M,L,KP)
            X2DPHI(KP,R)=X2DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVX20(M,L,KP)
            Q1PARA(KP,R)=Q1PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVQ10(M,L,KP)
            Q1PERP(KP,R)=Q1PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVQ10(M,L,KP)
            Q1DPHI(KP,R)=Q1DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVQ10(M,L,KP)
            Q2PARA(KP,R)=Q2PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVQ20(M,L,KP)
            Q2PERP(KP,R)=Q2PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVQ20(M,L,KP)
            Q2DPHI(KP,R)=Q2DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVQ20(M,L,KP)
            Q3PARA(KP,R)=Q3PARA(KP,R)+CTMP*SVPARA0(K,L,KP)*SVQ30(M,L,KP)
            Q3PERP(KP,R)=Q3PERP(KP,R)+CTMP*SVPERP0(K,L,KP)*SVQ30(M,L,KP)
            Q3DPHI(KP,R)=Q3DPHI(KP,R)+CTM3*SVDPHI0(K,L,KP)*SVQ30(M,L,KP)
            DPPARA(KP,R)=DPPARA(KP,R)+CTM2*SVPARA0(K,L,KP)*SVDP0(M,L,KP)
            DPPERP(KP,R)=DPPERP(KP,R)+CTM2*SVPERP0(K,L,KP)*SVDP0(M,L,KP)
            DPDPHI(KP,R)=DPDPHI(KP,R)+CTM4*SVDPHI0(K,L,KP)*SVDP0(M,L,KP)

            CTMP   = CTMP1*SF0(L,KP,2,1)*H3*H4
            CTM2   = CTMP1*SF0(L,KP,2,2)*H3*H4
            CTM3   = CTMP1*SF0(L,KP,2,3)*H3*H4
            CTM4   = CTMP1*SF0(L,KP,2,4)*H3*H4
            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   + CTMP*SVPARA01(K,L,KP)*SVX10(M,L,KP)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   + CTMP*SVPERP01(K,L,KP)*SVX10(M,L,KP)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   + CTM3*SVDPHI01(K,L,KP)*SVX10(M,L,KP)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   + CTMP*SVPARA01(K,L,KP)*SVX20(M,L,KP)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   + CTMP*SVPERP01(K,L,KP)*SVX20(M,L,KP)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   + CTM3*SVDPHI01(K,L,KP)*SVX20(M,L,KP)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   + CTMP*SVPARA01(K,L,KP)*SVQ10(M,L,KP)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   + CTMP*SVPERP01(K,L,KP)*SVQ10(M,L,KP)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   + CTM3*SVDPHI01(K,L,KP)*SVQ10(M,L,KP)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   + CTMP*SVPARA01(K,L,KP)*SVQ20(M,L,KP)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   + CTMP*SVPERP01(K,L,KP)*SVQ20(M,L,KP)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   + CTM3*SVDPHI01(K,L,KP)*SVQ20(M,L,KP)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   + CTMP*SVPARA01(K,L,KP)*SVQ30(M,L,KP)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   + CTMP*SVPERP01(K,L,KP)*SVQ30(M,L,KP)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   + CTM3*SVDPHI01(K,L,KP)*SVQ30(M,L,KP)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   + CTM2*SVPARA01(K,L,KP)*SVDP0(M,L,KP)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   + CTM2*SVPERP01(K,L,KP)*SVDP0(M,L,KP)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   + CTM4*SVDPHI01(K,L,KP)*SVDP0(M,L,KP)

            CTMP   =CTMP1*SF0(L,KP,3,1)*H3*H4
            CTM2   =CTMP1*SF0(L,KP,3,2)*H3*H4
            CTM3   =CTMP1*SF0(L,KP,3,3)*H3*H4
            CTM4   =CTMP1*SF0(L,KP,3,4)*H3*H4
            X1PARA(KP,R) =X1PARA(KP,R) 
     &                   +CTMP*SVPARA0(K,L,KP)*SVX101(M,L,KP)*VX1LNP(M)
            X1PERP(KP,R) =X1PERP(KP,R) 
     &                   +CTMP*SVPERP0(K,L,KP)*SVX101(M,L,KP)*VX1LNP(M)
            X1DPHI(KP,R) =X1DPHI(KP,R) 
     &                   +CTM3*SVDPHI0(K,L,KP)*SVX101(M,L,KP)*VX1LNP(M)
            X2PARA(KP,R) =X2PARA(KP,R) 
     &                   +CTMP*SVPARA0(K,L,KP)*SVX201(M,L,KP)*VX2LNP(M)
            X2PERP(KP,R) =X2PERP(KP,R) 
     &                   +CTMP*SVPERP0(K,L,KP)*SVX201(M,L,KP)*VX2LNP(M)
            X2DPHI(KP,R) =X2DPHI(KP,R) 
     &                   +CTM3*SVDPHI0(K,L,KP)*SVX201(M,L,KP)*VX2LNP(M)
            Q1PARA(KP,R) =Q1PARA(KP,R) 
     &                   +CTMP*SVPARA0(K,L,KP)*SVQ101(M,L,KP)*VQ1LNP(M)
            Q1PERP(KP,R) =Q1PERP(KP,R) 
     &                   +CTMP*SVPERP0(K,L,KP)*SVQ101(M,L,KP)*VQ1LNP(M)
            Q1DPHI(KP,R) =Q1DPHI(KP,R) 
     &                   +CTM3*SVDPHI0(K,L,KP)*SVQ101(M,L,KP)*VQ1LNP(M)
            Q2PARA(KP,R) =Q2PARA(KP,R) 
     &                   +CTMP*SVPARA0(K,L,KP)*SVQ201(M,L,KP)*VQ2LNP(M)
            Q2PERP(KP,R) =Q2PERP(KP,R) 
     &                   +CTMP*SVPERP0(K,L,KP)*SVQ201(M,L,KP)*VQ2LNP(M)
            Q2DPHI(KP,R) =Q2DPHI(KP,R) 
     &                   +CTM3*SVDPHI0(K,L,KP)*SVQ201(M,L,KP)*VQ2LNP(M)
            Q3PARA(KP,R) =Q3PARA(KP,R) 
     &                   +CTMP*SVPARA0(K,L,KP)*SVQ301(M,L,KP)*VQ3LNP(M)
            Q3PERP(KP,R) =Q3PERP(KP,R) 
     &                   +CTMP*SVPERP0(K,L,KP)*SVQ301(M,L,KP)*VQ3LNP(M)
            Q3DPHI(KP,R) =Q3DPHI(KP,R) 
     &                   +CTM3*SVDPHI0(K,L,KP)*SVQ301(M,L,KP)*VQ3LNP(M)
            DPPARA(KP,R) =DPPARA(KP,R) 
     &                   +CTM2*SVPARA0(K,L,KP)*SVDP01(M,L,KP)*VDPLNP(M)
            DPPERP(KP,R) =DPPERP(KP,R) 
     &                   +CTM2*SVPERP0(K,L,KP)*SVDP01(M,L,KP)*VDPLNP(M)
            DPDPHI(KP,R) =DPDPHI(KP,R) 
     &                   +CTM4*SVDPHI0(K,L,KP)*SVDP01(M,L,KP)*VDPLNP(M)

            X1PARA(KP,R) =X1PARA(KP,R) 
     &                   -CTMP*SVPARA01(K,L,KP)*SVX10(M,L,KP)*VX1LNP(M)
            X1PERP(KP,R) =X1PERP(KP,R) 
     &                   -CTMP*SVPERP01(K,L,KP)*SVX10(M,L,KP)*VX1LNP(M)
            X1DPHI(KP,R) =X1DPHI(KP,R) 
     &                   -CTM3*SVDPHI01(K,L,KP)*SVX10(M,L,KP)*VX1LNP(M)
            X2PARA(KP,R) =X2PARA(KP,R) 
     &                   -CTMP*SVPARA01(K,L,KP)*SVX20(M,L,KP)*VX2LNP(M)
            X2PERP(KP,R) =X2PERP(KP,R) 
     &                   -CTMP*SVPERP01(K,L,KP)*SVX20(M,L,KP)*VX2LNP(M)
            X2DPHI(KP,R) =X2DPHI(KP,R) 
     &                   -CTM3*SVDPHI01(K,L,KP)*SVX20(M,L,KP)*VX2LNP(M)
            Q1PARA(KP,R) =Q1PARA(KP,R) 
     &                   -CTMP*SVPARA01(K,L,KP)*SVQ10(M,L,KP)*VQ1LNP(M)
            Q1PERP(KP,R) =Q1PERP(KP,R) 
     &                   -CTMP*SVPERP01(K,L,KP)*SVQ10(M,L,KP)*VQ1LNP(M)
            Q1DPHI(KP,R) =Q1DPHI(KP,R) 
     &                   -CTM3*SVDPHI01(K,L,KP)*SVQ10(M,L,KP)*VQ1LNP(M)
            Q2PARA(KP,R) =Q2PARA(KP,R) 
     &                   -CTMP*SVPARA01(K,L,KP)*SVQ20(M,L,KP)*VQ2LNP(M)
            Q2PERP(KP,R) =Q2PERP(KP,R) 
     &                   -CTMP*SVPERP01(K,L,KP)*SVQ20(M,L,KP)*VQ2LNP(M)
            Q2DPHI(KP,R) =Q2DPHI(KP,R) 
     &                   -CTM3*SVDPHI01(K,L,KP)*SVQ20(M,L,KP)*VQ2LNP(M)
            Q3PARA(KP,R) =Q3PARA(KP,R) 
     &                   -CTMP*SVPARA01(K,L,KP)*SVQ30(M,L,KP)*VQ3LNP(M)
            Q3PERP(KP,R) =Q3PERP(KP,R) 
     &                   -CTMP*SVPERP01(K,L,KP)*SVQ30(M,L,KP)*VQ3LNP(M)
            Q3DPHI(KP,R) =Q3DPHI(KP,R) 
     &                   -CTM3*SVDPHI01(K,L,KP)*SVQ30(M,L,KP)*VQ3LNP(M)
            DPPARA(KP,R) =DPPARA(KP,R) 
     &                   -CTM2*SVPARA01(K,L,KP)*SVDP0(M,L,KP)*VDPLNP(M)
            DPPERP(KP,R) =DPPERP(KP,R) 
     &                   -CTM2*SVPERP01(K,L,KP)*SVDP0(M,L,KP)*VDPLNP(M)
            DPDPHI(KP,R) =DPDPHI(KP,R) 
     &                   -CTM4*SVDPHI01(K,L,KP)*SVDP0(M,L,KP)*VDPLNP(M)
         ENDIF
         ENDDO
         ENDDO

      CASE (4)
         DO KP=1,NSPECIES
         DO L=1,MLMAX
         IF (ABS(RLM(L)).LT.0.1.AND.KPARTICLE.EQ.0) THEN
            R=2
            CTMP   = VI0(1,KP)*H3
            CTM2   = VI0(2,KP)*H3
            CTM3   = VI0(3,KP)*H3
            CTM4   = VI0(4,KP)*H3
            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA0(K,L)*VX10(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP0(K,L)*VX10(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI0(K,L)*VX10(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA0(K,L)*VX20(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP0(K,L)*VX20(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI0(K,L)*VX20(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA0(K,L)*VQ10(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP0(K,L)*VQ10(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI0(K,L)*VQ10(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA0(K,L)*VQ20(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP0(K,L)*VQ20(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI0(K,L)*VQ20(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA0(K,L)*VQ30(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP0(K,L)*VQ30(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI0(K,L)*VQ30(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA0(K,L)*VDP0(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP0(K,L)*VDP0(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI0(K,L)*VDP0(M,L)

            CTMP   = VI02(1,KP)*H3
            CTM2   = VI02(2,KP)*H3
            CTM3   = VI02(3,KP)*H3
            CTM4   = VI02(4,KP)*H3
            X1PARA(KP,R) = X1PARA(KP,R) + CTMP*VPARA01(K,L)*VX10(M,L)
            X1PERP(KP,R) = X1PERP(KP,R) + CTMP*VPERP01(K,L)*VX10(M,L)
            X1DPHI(KP,R) = X1DPHI(KP,R) + CTM3*VDPHI01(K,L)*VX10(M,L)
            X2PARA(KP,R) = X2PARA(KP,R) + CTMP*VPARA01(K,L)*VX20(M,L)
            X2PERP(KP,R) = X2PERP(KP,R) + CTMP*VPERP01(K,L)*VX20(M,L)
            X2DPHI(KP,R) = X2DPHI(KP,R) + CTM3*VDPHI01(K,L)*VX20(M,L)
            Q1PARA(KP,R) = Q1PARA(KP,R) + CTMP*VPARA01(K,L)*VQ10(M,L)
            Q1PERP(KP,R) = Q1PERP(KP,R) + CTMP*VPERP01(K,L)*VQ10(M,L)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) + CTM3*VDPHI01(K,L)*VQ10(M,L)
            Q2PARA(KP,R) = Q2PARA(KP,R) + CTMP*VPARA01(K,L)*VQ20(M,L)
            Q2PERP(KP,R) = Q2PERP(KP,R) + CTMP*VPERP01(K,L)*VQ20(M,L)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) + CTM3*VDPHI01(K,L)*VQ20(M,L)
            Q3PARA(KP,R) = Q3PARA(KP,R) + CTMP*VPARA01(K,L)*VQ30(M,L)
            Q3PERP(KP,R) = Q3PERP(KP,R) + CTMP*VPERP01(K,L)*VQ30(M,L)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) + CTM3*VDPHI01(K,L)*VQ30(M,L)
            DPPARA(KP,R) = DPPARA(KP,R) + CTM2*VPARA01(K,L)*VDP0(M,L)
            DPPERP(KP,R) = DPPERP(KP,R) + CTM2*VPERP01(K,L)*VDP0(M,L)
            DPDPHI(KP,R) = DPDPHI(KP,R) + CTM4*VDPHI01(K,L)*VDP0(M,L)

            CTMP   = VI03(1,KP)*H3
            CTM2   = VI03(2,KP)*H3
            CTM3   = VI03(3,KP)*H3
            CTM4   = VI03(4,KP)*H3
            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VX101(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VX101(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VX101(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VX201(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VX201(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VX201(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VQ101(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VQ101(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VQ101(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VQ201(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VQ201(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VQ201(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   + CTMP*VPARA0(K,L)*VQ301(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   + CTMP*VPERP0(K,L)*VQ301(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   + CTM3*VDPHI0(K,L)*VQ301(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   + CTM2*VPARA0(K,L)*VDP01(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   + CTM2*VPERP0(K,L)*VDP01(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   + CTM4*VDPHI0(K,L)*VDP01(M,L)*VDPLNP(M)

            X1PARA(KP,R) = X1PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VX10(M,L)*VX1LNP(M)
            X1PERP(KP,R) = X1PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VX10(M,L)*VX1LNP(M)
            X1DPHI(KP,R) = X1DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VX10(M,L)*VX1LNP(M)
            X2PARA(KP,R) = X2PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VX20(M,L)*VX2LNP(M)
            X2PERP(KP,R) = X2PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VX20(M,L)*VX2LNP(M)
            X2DPHI(KP,R) = X2DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VX20(M,L)*VX2LNP(M)
            Q1PARA(KP,R) = Q1PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VQ10(M,L)*VQ1LNP(M)
            Q1PERP(KP,R) = Q1PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VQ10(M,L)*VQ1LNP(M)
            Q1DPHI(KP,R) = Q1DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VQ10(M,L)*VQ1LNP(M)
            Q2PARA(KP,R) = Q2PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VQ20(M,L)*VQ2LNP(M)
            Q2PERP(KP,R) = Q2PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VQ20(M,L)*VQ2LNP(M)
            Q2DPHI(KP,R) = Q2DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VQ20(M,L)*VQ2LNP(M)
            Q3PARA(KP,R) = Q3PARA(KP,R) 
     &                   - CTMP*VPARA01(K,L)*VQ30(M,L)*VQ3LNP(M)
            Q3PERP(KP,R) = Q3PERP(KP,R) 
     &                   - CTMP*VPERP01(K,L)*VQ30(M,L)*VQ3LNP(M)
            Q3DPHI(KP,R) = Q3DPHI(KP,R) 
     &                   - CTM3*VDPHI01(K,L)*VQ30(M,L)*VQ3LNP(M)
            DPPARA(KP,R) = DPPARA(KP,R) 
     &                   - CTM2*VPARA01(K,L)*VDP0(M,L)*VDPLNP(M)
            DPPERP(KP,R) = DPPERP(KP,R) 
     &                   - CTM2*VPERP01(K,L)*VDP0(M,L)*VDPLNP(M)
            DPDPHI(KP,R) = DPDPHI(KP,R) 
     &                   - CTM4*VDPHI01(K,L)*VDP0(M,L)*VDPLNP(M)
         ENDIF
         ENDDO
         ENDDO
      CASE (5)   
         R=1
         DO KP=1,NSPECIES
         L = M-K + (M2-M1) + 1
         X1PARA(KP,R) =((ZKIA(1,KP)+ZKIA1(1,KP))*ZGL0PA(L) 
     &                 +ZKIA2(1,KP)*ZGL0PA1(L))*LAMH*H3
         X1PERP(KP,R) =((ZKIA(1,KP)+ZKIA1(1,KP))*ZGL0PE(L) 
     &                 +ZKIA2(1,KP)*ZGL0PE1(L))*LAMH*H3
         X1DPHI(KP,R) =((ZKIA(4,KP)+ZKIA1(4,KP))*ZGL0DP(L) 
     &                 +ZKIA2(4,KP)*ZGL0DP1(L))*LAMH*H3
         Q1PARA(KP,R) =((ZKIA(2,KP)+ZKIA1(2,KP))*ZGL1PA(L) 
     &                 +ZKIA2(2,KP)*ZGL1PA1(L))*LAMH*H3
         Q1PERP(KP,R) =((ZKIA(2,KP)+ZKIA1(2,KP))*ZGL1PE(L) 
     &                 +ZKIA2(2,KP)*ZGL1PE1(L))*LAMH*H3
         Q1DPHI(KP,R) =((ZKIA(5,KP)+ZKIA1(5,KP))*ZGL1DP(L) 
     &                 +ZKIA2(5,KP)*ZGL1DP1(L))*LAMH*H3
         Q2PARA(KP,R) =((ZKIA(2,KP)+ZKIA1(2,KP))*ZGL2PA(L) 
     &                 +ZKIA2(2,KP)*ZGL2PA1(L))*LAMH*H3
         Q2PERP(KP,R) =((ZKIA(2,KP)+ZKIA1(2,KP))*ZGL2PE(L) 
     &                 +ZKIA2(2,KP)*ZGL2PE1(L))*LAMH*H3
         Q2DPHI(KP,R) =((ZKIA(5,KP)+ZKIA1(5,KP))*ZGL2DP(L) 
     &                 +ZKIA2(5,KP)*ZGL2DP1(L))*LAMH*H3
         Q3PARA(KP,R) =((ZKIA(2,KP)+ZKIA1(2,KP))*ZGL3PA(L) 
     &                 +ZKIA2(2,KP)*ZGL3PA1(L))*LAMH*H3
         Q3PERP(KP,R) =((ZKIA(2,KP)+ZKIA1(2,KP))*ZGL3PE(L) 
     &                 +ZKIA2(2,KP)*ZGL3PE1(L))*LAMH*H3
         Q3DPHI(KP,R) =((ZKIA(5,KP)+ZKIA1(5,KP))*ZGL3DP(L) 
     &                 +ZKIA2(5,KP)*ZGL3DP1(L))*LAMH*H3
         DPPARA(KP,R) =((ZKIA(3,KP)+ZKIA1(3,KP))*ZGL0PA(L) 
     &                 +ZKIA2(3,KP)*ZGL0PA1(L))*LAMH*H3
         DPPERP(KP,R) =((ZKIA(3,KP)+ZKIA1(3,KP))*ZGL0PE(L) 
     &                 +ZKIA2(3,KP)*ZGL0PE1(L))*LAMH*H3
         DPDPHI(KP,R) =((ZKIA(6,KP)+ZKIA1(6,KP))*ZGL0DP(L) 
     &                 +ZKIA2(6,KP)*ZGL0DP1(L))*LAMH*H3
         ENDDO
      END SELECT   
      
      IF (KGRID.EQ.1) THEN
         VX1PARA(K,M,JS_MAT)=VX1PARA(K,M,JS_MAT)+SUM(SUM(X1PARA,1),1)
         VX1PERP(K,M,JS_MAT)=VX1PERP(K,M,JS_MAT)+SUM(SUM(X1PERP,1),1)
         VX2PARA(K,M,JS_MAT)=VX2PARA(K,M,JS_MAT)+SUM(SUM(X2PARA,1),1)
         VX2PERP(K,M,JS_MAT)=VX2PERP(K,M,JS_MAT)+SUM(SUM(X2PERP,1),1) 
         VQ1PARA(K,M,JS_MAT)=VQ1PARA(K,M,JS_MAT)+SUM(SUM(Q1PARA,1),1)
         VQ1PERP(K,M,JS_MAT)=VQ1PERP(K,M,JS_MAT)+SUM(SUM(Q1PERP,1),1)
         VQ2PARA(K,M,JS_MAT)=VQ2PARA(K,M,JS_MAT)+SUM(SUM(Q2PARA,1),1) 
         VQ2PERP(K,M,JS_MAT)=VQ2PERP(K,M,JS_MAT)+SUM(SUM(Q2PERP,1),1) 
         VQ3PARA(K,M,JS_MAT)=VQ3PARA(K,M,JS_MAT)+SUM(SUM(Q3PARA,1),1) 
         VQ3PERP(K,M,JS_MAT)=VQ3PERP(K,M,JS_MAT)+SUM(SUM(Q3PERP,1),1)
         VDPPARA(K,M,JS_MAT)=VDPPARA(K,M,JS_MAT)+SUM(SUM(DPPARA,1),1) 
         VDPPERP(K,M,JS_MAT)=VDPPERP(K,M,JS_MAT)+SUM(SUM(DPPERP,1),1)

         IF (INCDPHI.GT.0) THEN
         VX1DPHI(K,M,JS_MAT)=VX1DPHI(K,M,JS_MAT)+SUM(SUM(X1DPHI,1),1)
         VX2DPHI(K,M,JS_MAT)=VX2DPHI(K,M,JS_MAT)+SUM(SUM(X2DPHI,1),1)
         VQ1DPHI(K,M,JS_MAT)=VQ1DPHI(K,M,JS_MAT)+SUM(SUM(Q1DPHI,1),1)
         VQ2DPHI(K,M,JS_MAT)=VQ2DPHI(K,M,JS_MAT)+SUM(SUM(Q2DPHI,1),1)
         VQ3DPHI(K,M,JS_MAT)=VQ3DPHI(K,M,JS_MAT)+SUM(SUM(Q3DPHI,1),1)
         VDPDPHI(K,M,JS_MAT)=VDPDPHI(K,M,JS_MAT)+SUM(SUM(DPDPHI,1),1)
         ENDIF

         IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &       .AND.INCKIN.GT.0) THEN
         VX1PARE(K,M,JS_MAT)=VX1PARE(K,M,JS_MAT)+SUM(X1PARA(2,:))
         VX1PERE(K,M,JS_MAT)=VX1PERE(K,M,JS_MAT)+SUM(X1PERP(2,:))
         VX2PARE(K,M,JS_MAT)=VX2PARE(K,M,JS_MAT)+SUM(X2PARA(2,:))
         VX2PERE(K,M,JS_MAT)=VX2PERE(K,M,JS_MAT)+SUM(X2PERP(2,:)) 
         VQ1PARE(K,M,JS_MAT)=VQ1PARE(K,M,JS_MAT)+SUM(Q1PARA(2,:))
         VQ1PERE(K,M,JS_MAT)=VQ1PERE(K,M,JS_MAT)+SUM(Q1PERP(2,:))
         VQ2PARE(K,M,JS_MAT)=VQ2PARE(K,M,JS_MAT)+SUM(Q2PARA(2,:)) 
         VQ2PERE(K,M,JS_MAT)=VQ2PERE(K,M,JS_MAT)+SUM(Q2PERP(2,:)) 
         VQ3PARE(K,M,JS_MAT)=VQ3PARE(K,M,JS_MAT)+SUM(Q3PARA(2,:)) 
         VQ3PERE(K,M,JS_MAT)=VQ3PERE(K,M,JS_MAT)+SUM(Q3PERP(2,:))
         VDPPARE(K,M,JS_MAT)=VDPPARE(K,M,JS_MAT)+SUM(DPPARA(2,:)) 
         VDPPERE(K,M,JS_MAT)=VDPPERE(K,M,JS_MAT)+SUM(DPPERP(2,:))
         ENDIF
      ELSE
         VX1PARAM(K,M,JS_MAT)=VX1PARAM(K,M,JS_MAT)+SUM(SUM(X1PARA,1),1)
         VX1PERPM(K,M,JS_MAT)=VX1PERPM(K,M,JS_MAT)+SUM(SUM(X1PERP,1),1)
         VX2PARAM(K,M,JS_MAT)=VX2PARAM(K,M,JS_MAT)+SUM(SUM(X2PARA,1),1)
         VX2PERPM(K,M,JS_MAT)=VX2PERPM(K,M,JS_MAT)+SUM(SUM(X2PERP,1),1) 
         VQ1PARAM(K,M,JS_MAT)=VQ1PARAM(K,M,JS_MAT)+SUM(SUM(Q1PARA,1),1)
         VQ1PERPM(K,M,JS_MAT)=VQ1PERPM(K,M,JS_MAT)+SUM(SUM(Q1PERP,1),1)
         VQ2PARAM(K,M,JS_MAT)=VQ2PARAM(K,M,JS_MAT)+SUM(SUM(Q2PARA,1),1) 
         VQ2PERPM(K,M,JS_MAT)=VQ2PERPM(K,M,JS_MAT)+SUM(SUM(Q2PERP,1),1) 
         VQ3PARAM(K,M,JS_MAT)=VQ3PARAM(K,M,JS_MAT)+SUM(SUM(Q3PARA,1),1) 
         VQ3PERPM(K,M,JS_MAT)=VQ3PERPM(K,M,JS_MAT)+SUM(SUM(Q3PERP,1),1)
         VDPPARAM(K,M,JS_MAT)=VDPPARAM(K,M,JS_MAT)+SUM(SUM(DPPARA,1),1) 
         VDPPERPM(K,M,JS_MAT)=VDPPERPM(K,M,JS_MAT)+SUM(SUM(DPPERP,1),1)

         IF (INCDPHI.GT.0) THEN
         VX1DPHIM(K,M,JS_MAT)=VX1DPHIM(K,M,JS_MAT)+SUM(SUM(X1DPHI,1),1)
         VX2DPHIM(K,M,JS_MAT)=VX2DPHIM(K,M,JS_MAT)+SUM(SUM(X2DPHI,1),1)
         VQ1DPHIM(K,M,JS_MAT)=VQ1DPHIM(K,M,JS_MAT)+SUM(SUM(Q1DPHI,1),1)
         VQ2DPHIM(K,M,JS_MAT)=VQ2DPHIM(K,M,JS_MAT)+SUM(SUM(Q2DPHI,1),1)
         VQ3DPHIM(K,M,JS_MAT)=VQ3DPHIM(K,M,JS_MAT)+SUM(SUM(Q3DPHI,1),1)
         VDPDPHIM(K,M,JS_MAT)=VDPDPHIM(K,M,JS_MAT)+SUM(SUM(DPDPHI,1),1)
         ENDIF

         IF (IDIAMTE.EQ.0.AND.(IDIAMB.EQ.1.OR.IDIAMB.EQ.3)
     &       .AND.INCKIN.GT.0) THEN
         VX1PAREM(K,M,JS_MAT)=VX1PAREM(K,M,JS_MAT)+SUM(X1PARA(2,:))
         VX1PEREM(K,M,JS_MAT)=VX1PEREM(K,M,JS_MAT)+SUM(X1PERP(2,:))
         VX2PAREM(K,M,JS_MAT)=VX2PAREM(K,M,JS_MAT)+SUM(X2PARA(2,:))
         VX2PEREM(K,M,JS_MAT)=VX2PEREM(K,M,JS_MAT)+SUM(X2PERP(2,:)) 
         VQ1PAREM(K,M,JS_MAT)=VQ1PAREM(K,M,JS_MAT)+SUM(Q1PARA(2,:))
         VQ1PEREM(K,M,JS_MAT)=VQ1PEREM(K,M,JS_MAT)+SUM(Q1PERP(2,:))
         VQ2PAREM(K,M,JS_MAT)=VQ2PAREM(K,M,JS_MAT)+SUM(Q2PARA(2,:)) 
         VQ2PEREM(K,M,JS_MAT)=VQ2PEREM(K,M,JS_MAT)+SUM(Q2PERP(2,:)) 
         VQ3PAREM(K,M,JS_MAT)=VQ3PAREM(K,M,JS_MAT)+SUM(Q3PARA(2,:)) 
         VQ3PEREM(K,M,JS_MAT)=VQ3PEREM(K,M,JS_MAT)+SUM(Q3PERP(2,:))
         VDPPAREM(K,M,JS_MAT)=VDPPAREM(K,M,JS_MAT)+SUM(DPPARA(2,:)) 
         VDPPEREM(K,M,JS_MAT)=VDPPEREM(K,M,JS_MAT)+SUM(DPPERP(2,:))
         ENDIF
      ENDIF

C     SET ENERGY COMPONENT MATRIX         
      DO KP=1,NSPECIES
      IF (KPARTICLE.EQ.0) THEN
         R=1
C        TRAPPED PARTICLES
         IF (ICASE.EQ.5) THEN
C        ADIABATIC COMPONENT
         CALL SETDWKCOMPMAT ( KP,2,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))
         ENDIF         
         IF (ICASE.EQ.1.OR.ICASE.EQ.2.OR.ICASE.EQ.3) THEN
C        BOUNCE COMPONENT             
         CALL SETDWKCOMPMAT ( KP,4,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))
         ENDIF
         IF (ICASE.EQ.1.OR.ICASE.EQ.2.OR.ICASE.EQ.3.OR.ICASE.EQ.4) THEN
C        PRECESSION COMPONENT     
         R=2
         CALL SETDWKCOMPMAT ( KP,5,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))
         ENDIF
      ELSE
C        PASSING PARTICLES
         R=1
         IF (ICASE.EQ.5) THEN
C        ADIABATIC COMPONENT
         CALL SETDWKCOMPMAT ( KP,1,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))               
         ENDIF
         IF (ICASE.EQ.1.OR.ICASE.EQ.2.OR.ICASE.EQ.3) THEN
C         TRANSIT COMPONENT
         CALL SETDWKCOMPMAT ( KP,3,K,M,
     &   X1PARA(KP,R),X1PERP(KP,R),X2PARA(KP,R),X2PERP(KP,R),
     &   Q1PARA(KP,R),Q1PERP(KP,R),Q2PARA(KP,R),Q2PERP(KP,R),
     &   Q3PARA(KP,R),Q3PERP(KP,R))         
         ENDIF         
      ENDIF
      ENDDO
      
      ENDDO
      ENDDO

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0.AND.ICASE.EQ.2) THEN
         J = INT(NLAMK0(JS,KGRID)/2.)
         IF (ABS(LAMM(J+1)-LAM).LT.1.0E-13) THEN
            K=MIN(2-M1+1,M2-M1+1)
            M=K
            L=1
            WRITE(*,121) LAM,VI(1,L,1),VI2(1,L,1),VI3(1,L,1),
     &                   VPARA(K,L),VPARA1(K,L),VPERP(K,L),VPERP1(K,L),
     &                   VX1(M,L),VX11(M,L),VX1LNP(M),
     &                   VX2(M,L),VX21(M,L),VX2LNP(M),
     &                   VQ1(M,L),VQ11(M,L),VQ1LNP(M),
     &                   VQ2(M,L),VQ21(M,L),VQ2LNP(M),
     &                   VQ3(M,L),VQ31(M,L),VQ3LNP(M),
     &                   VDP(M,L),VDP1(M,L),VDPLNP(M)
         ENDIF
 121     FORMAT('KJPFILL1:LAM,VI,VI2,VI3,VPARA,VPARA1,VPERP,VPERP1,
     &           VX1,VX11,VX1LNP,VX2,VX21,VX2LNP,VQ1,VQ11,VQ1LNP,
     &           VQ2,VQ21,VQ2LNP,VQ3,VQ31,VQ3LNP,VDP,VDP1,VDPLNP',
     &           47(E13.5,1X))
      ENDIF
      DEALLOCATE (X1PARA,X2PARA,Q1PARA,Q2PARA,Q3PARA,DPPARA,
     &            X1PERP,X2PERP,Q1PERP,Q2PERP,Q3PERP,DPPERP,
     &            X1DPHI,X2DPHI,Q1DPHI,Q2DPHI,Q3DPHI,DPDPHI)
      RETURN
      END

C=======================================================================
C COMPUTE POLOIDAL ANGLES OF TURNING POINTS,                           =
C AND H'(CHI) AT TURNING POINTS                                        =
C USING 1D CUBIC SPLINE                                                =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KTURN_OLD(JS,KGRID)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER JS,KGRID,J,K
      REAL*8  CHIL0,CHIU0

      INTEGER KCHECK

      DO J=1,NCHI
         RW1(J) = HK(JS,J,KGRID)
      ENDDO
      RW1(NCHI+1)=RW1(1)

C     FIND INITIAL GUESS
      J=1
 100  IF ((RW1(J)-LAM)*(RW1(J+1)-LAM).GT.0.0) THEN
         J=J+1
         GOTO 100
      ENDIF
      CHIU0=(RCHI(J)+RCHI(J+1))*0.5

      J=J+1
 200  IF ((RW1(J)-LAM)*(RW1(J+1)-LAM).GT.0.0) THEN
         J=J+1
         GOTO 200
      ENDIF
      CHIL0=(RCHI(J)+RCHI(J+1))*0.5

      CALL SPLINE1DR(HPL,CHIL,LAM,CHIL0,RW1,RCHI,NCHI+1,RCHI2)
      CALL SPLINE1DR(HPU,CHIU,LAM,CHIU0,RW1,RCHI,NCHI+1,RCHI2)

      IF (HPU.GT.0.AND.HPL.LT.0.AND.1.EQ.0) THEN
         CHIL0 = CHIL
         CHIL  = CHIU
         CHIU  = CHIL0
         CHIL0 = HPL
         HPL   = HPU
         HPU   = CHIL0
      ENDIF

      KCHECK = 0
      IF (KCHECK.EQ.1) THEN
         IF (JS.EQ.JS0.AND.KGRID.EQ.1) THEN
            WRITE(*,*) 'CHECK KTURN: LAM CHIL CHIU HPL HPU'
            WRITE(*,110) LAM,CHIL,CHIU,HPL,HPU
         ENDIF
 110     FORMAT(5(E12.4,1X))
      ENDIF

      IF (HPL.LE.0.0.OR.HPU.GE.0.0) 
     &   WRITE(*,*) 'HPL,HPU',HPL,HPU,JS,KGRID

      IF (HPL.LE.0.0.OR.HPU.GE.0.0)
     &  STOP 'KINETIC:HPL<0 | HPU>0,TRY KSMOOTH=1'
      RETURN
      END
      
      SUBROUTINE KTURN(JS,KGRID)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER JS,KGRID,J1,J2,K
      REAL*8  CHIL0,CHIU0
      INTEGER STARTPOT
      INTEGER ARRARYSIZE,P
      INTEGER KCHECK
      INTEGER J
      
      STARTPOT = 1
      
      DO J=1,NCHI
         RW1(J) = HK(JS,J,KGRID)
      ENDDO
      RW1(NCHI+1)=RW1(1)

      ARRARYSIZE = NCHI+1
      
      P=STARTPOT - 1
      IF ( P .LT. 0) THEN
        WRITE (*,*) 'STARTPOT IS LESS THAN 1 JS=',JS,'KGRID=',KGRID
        STOP 'STARTPOT IS LESS THAN 1.'
      ENDIF
      J1 = MOD(P,ARRARYSIZE) + 1
      J2 = MOD (P+1,ARRARYSIZE) + 1
      DO WHILE ((RW1(J1)-LAM)*(RW1(J2)-LAM).GT.0.0)         
         P=MOD(P+1,ARRARYSIZE)
         
         IF ( P .EQ. (STARTPOT - 1) ) THEN
            WRITE (*,*) 'NOT FIND CHIU/L 1 JS=',JS,'KGRID=',KGRID
            STOP 'NOT FIND CHIU/L 1'
         ENDIF
         J1 = MOD(P,ARRARYSIZE) + 1
         J2 = MOD (P+1,ARRARYSIZE) + 1         
      ENDDO   
      CHIU0=(RCHI(J1)+RCHI(J2))*0.5

      P= MOD(STARTPOT - 2 + ARRARYSIZE, ARRARYSIZE)
      J1 = MOD(P,ARRARYSIZE) + 1
      J2 = MOD (P+1,ARRARYSIZE) + 1
      DO WHILE ((RW1(J1)-LAM)*(RW1(J2)-LAM).GT.0.0)
         P=MOD(P - 1 + ARRARYSIZE,ARRARYSIZE)
         
         IF ( P .EQ. MOD(STARTPOT - 2 + ARRARYSIZE, ARRARYSIZE) ) THEN
            WRITE (*,*) 'NOT FIND CHIU/L 2 JS=',JS,'KGRID=',KGRID
            STOP 'NOT FIND CHIU/L 2'
         ENDIF
         
         J1 = MOD(P,ARRARYSIZE) + 1
         J2 = MOD (P+1,ARRARYSIZE) + 1
      ENDDO
      
      CHIL0=(RCHI(J1)+RCHI(J2))*0.5

      CALL SPLINE1DR(HPL,CHIL,LAM,CHIL0,RW1,RCHI,NCHI+1,RCHI2)
      CALL SPLINE1DR(HPU,CHIU,LAM,CHIU0,RW1,RCHI,NCHI+1,RCHI2)
      IF (HPU * HPL .GE. 0) THEN
        WRITE (*,*) 'HPU AND HPL ARE UNPHYSICAL JS=',JS,'KGRID=',KGRID
        STOP 'HPU AND HPL ARE UNPHYSICAL'
      ENDIF
C     IF (HPU.GT.0.AND.HPL.LT.0.AND.0.EQ.1) THEN
      IF (HPU.GT.0.AND.HPL.LT.0) THEN
         CHIL0 = CHIL
         CHIL  = CHIU
         CHIU  = CHIL0
         CHIL0 = HPL
         HPL   = HPU
         HPU   = CHIL0
      ENDIF

      IF (CHIU-CHIL.GE.PI) THEN
         WRITE(*,*) 'KTURN: MULTIPLE TRAPPING WITH JS,KGRID,CHIL,CHIU=',
     *                     JS,KGRID,CHIL,CHIU
         STOP 'KTURN: MULTIPLE TRAPPING'
      ENDIF

      KCHECK = 0
      IF (KCHECK.EQ.1) THEN
         IF (JS.EQ.JS0.AND.KGRID.EQ.1) THEN
            WRITE(*,*) 'CHECK KTURN: LAM CHIL CHIU HPL HPU'
            WRITE(*,110) LAM,CHIL,CHIU,HPL,HPU
         ENDIF
 110     FORMAT(5(E12.4,1X))
      ENDIF

      IF (HPL.LE.0.0.OR.HPU.GE.0.0) 
     &   STOP 'KINETIC:HPL<0 | HPU>0,TRY KSMOOTH=1'

      RETURN
      END

C=======================================================================
C COMPUTE INTEGRATION POINTS ALONG CHI, BETWEEN CHIL AND CHIU          =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KCHI(KPARTICLE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER J,KPARTICLE
      REAL*8 DIFFERCHI
      IF (KPARTICLE.EQ.1) THEN
         NCHI1 = NCHI0 
         NCHI2 = 2*NCHI1
         RCHIHK = (CHIU-CHIL)/NCHI1
         RCHIK(1) = CHIL
         DO J=1,NCHI1
            RCHIK(2*J)   = CHIL + (2.0*J-1-WK)*RCHIHK/2.0
            RCHIK(2*J+1) = CHIL + (2.0*J-1+WK)*RCHIHK/2.0
         ENDDO
         RCHIK(NCHI2+2) = CHIU
      ELSEIF (KPARTICLE.EQ.0) THEN
         RCHIHK = 2.0*PI/(NCHI0-9)
         NCHI1  = INT((2*PI-CHIL+CHIU)/RCHIHK) + 8
         NCHI2  = 2*NCHI1
         RCHIHK   = (DIFFERCHI(CHIU,CHIL))/NCHI1
         RCHIK(1) = CHIL
         DO J=1,NCHI1
            RCHIK(2*J)   = CHIL + (2.0*J-1-WK)*RCHIHK/2.0
            RCHIK(2*J+1) = CHIL + (2.0*J-1+WK)*RCHIHK/2.0
         ENDDO
         RCHIK(NCHI2+2) = CHIU
         DO J=2,NCHI2+1
            IF (RCHIK(J).GT.(PI+PI)) RCHIK(J) = RCHIK(J)-PI-PI
         ENDDO
      ENDIF

      RETURN
      END

      
C=======================================================================
C COMPUTE PITCH ANGLE OF BEAM DRIVEN EP ALONG LFS MIDPLANE             =
C USING NORMALISED TANGENCY RADII OF BEAM INJECTION AS INPUT           =
C YQL, 11-2010                                                         =
C=======================================================================
      SUBROUTINE KBEAMLAM

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE MPIENV
      IMPLICIT NONE

      INTEGER JS,J
      REAL*8  VPHIV,VRV,BRB,BPHIB,BLAM
      INTEGER KCHECK

      KCHECK = 0
      J      = 1
      
      IF (ISMPIRUN.EQ.0) THEN
        WRITE(*,*) 'PITCH ANGLE OF BEAM INJECTED EP:'
      ENDIF
      
      DO JS=1,NRP1
         VPHIV = RTAN/REQ(JS,J)
         VRV   = -SQRT(1.-VPHIV**2)
         BRB   = DPSIDS(JS)*RDSDZ(JS,J)/REQ(JS,J)/BK(JS,J,1)
         BPHIB = T(JS)/REQ(JS,J)/BK(JS,J,1)
         BLAM  = B0K*(1-(VRV*BRB+VPHIV*BPHIB)**2)/BK(JS,J,1) 
         IF (ISMPIRUN.EQ.0) THEN
            WRITE(*,100) CS(JS),BLAM
         ENDIF
      ENDDO
 100  FORMAT(E16.8,2X,E16.8)

      RETURN
      END

C=======================================================================
C COMPUTE EQUILIBRIUM QUANTITIES AT INTEGRATION POINTS ALONG CHI       =
C USING 1D CUBIC SPLINE                                                =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KEQUIL(JS,KGRID)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE MPIENV
      IMPLICIT NONE

      INTEGER JS,KGRID,J
      REAL*8  H1,H2,TMP1,TMP2,TMP3
      INTEGER KCHECK
      KCHECK=0

C     Q-VALUE
      IF (KGRID.EQ.1) THEN
         RQK = Q(JS)
         RDPSIK=DPSIDS(JS)
      ELSEIF (KGRID.EQ.2) THEN
         RQK = QM(JS)
         RDPSIK=DPSIDSM(JS)
      ENDIF
      
C     COEFFICIENT FOR G-FACTORS JB
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RJB(J) = RJA(JS,J)*BK(JS,J,KGRID)
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RJB(J) = RJAM(JS,J)*BK(JS,J,KGRID)
         ENDDO
      ENDIF
      RJB(NCHI+1)=RJB(1)

C     NOTE THAT ALL COEFFICIENTS FOR H-FACTORS BELOW ARE DEFINED AS 
C     C*JB/DPSIDS, WHERE C IS DEFINED IN MANUAL

C     COEFFICIENTS FOR XI^1 IN H-FACTOR
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RX1P(J)=2.0*PPEQ(JS)*RJA(JS,J)/BK(JS,J,KGRID)
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RX1P(J)=2.0*PPEQM(JS)*RJAM(JS,J)/BK(JS,J,KGRID)
         ENDDO
      ENDIF
      RX1P(NCHI+1)=RX1P(1)

C     ALSO COMPUTE COEFFICIENTS FOR DRIFT FREQUENCY
      DO J=2,NCHI-1
         RW2(J)=(HK(JS,J+1,KGRID)-HK(JS,J-1,KGRID))/2/RCHIH
      ENDDO
      RW2(1)=(HK(JS,2,KGRID)-HK(JS,NCHI,KGRID))/2/RCHIH
      RW2(NCHI)=(HK(JS,1,KGRID)-HK(JS,NCHI-1,KGRID))/2/RCHIH
      IF (KGRID.EQ.1) THEN
         IF (JS.LT.NRP1) THEN
            H1 = (CS(JS)-CS(JS-1))/2
            H2 = (CS(JS+1)-CS(JS))/2
            DO J=1,NCHI
               RW1(J)=(H1/H2*BK(JS,J,2)-H2/H1*BK(JS-1,J,2))/(H1+H2) -
     &                (H1-H2)*BK(JS,J,1)/H1/H2
            ENDDO
         ELSEIF (JS.EQ.NRP1) THEN
         H1 = CS(JS)-CS(JS-1)
         DO J=1,NCHI
            RW1(J)=(BK(JS-1,J,1)+3*BK(JS,J,1)-4*BK(JS-1,J,2))/H1
         ENDDO
         ENDIF
         DO J=1,NCHI
            RX1B(J)=RW1(J)*RJA(JS,J)/DPSIDS(JS) + 
     &              RW2(J)*DPSIDS(JS)*G12L(JS,J)/RJA(JS,J)/B0K
            RX1R(J)=-RW2(J)*T(JS)*DROT(JS)/B0K
            RDMU(J)=-(RW1(J)+RW2(J)*B0K/HK(JS,J,1)**2*G12L(JS,J)/
     &                G22L(JS,J))/SQRT(HK(JS,J,1))*RJA(JS,J)/
     &                DPSIDS(JS)**2
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         H1 = CS(JS+1) - CS(JS)
         DO J=1,NCHI
            RW1(J)=(BK(JS+1,J,1)-BK(JS,J,1))/H1
            RX1B(J)=RW1(J)*RJAM(JS,J)/DPSIDSM(JS) + 
     &              RW2(J)*DPSIDSM(JS)*G12LM(JS,J)/RJAM(JS,J)/B0K
            RX1R(J)=-RW2(J)*TM(JS)*DROTM(JS)/B0K
            RDMU(J)=-(RW1(J)+RW2(J)*B0K/HK(JS,J,2)**2*G12LM(JS,J)/
     &                G22LM(JS,J))/SQRT(HK(JS,J,2))*RJAM(JS,J)/
     &                DPSIDSM(JS)**2
         ENDDO
      ENDIF
      RX1B(NCHI+1)=RX1B(1)
      RX1R(NCHI+1)=RX1R(1)
      RDMU(NCHI+1)=RDMU(1)

C     CHECK DB/DS
      IF (KCHECK.EQ.2.AND.JS.EQ.JS0.AND.KGRID.EQ.1) THEN
         WRITE(*,*) 'CHECK KEQUIL:DPSIDS=',RDPSIK
         WRITE(*,*) 'CHECK KEQUIL:CHI DB/DS'
         DO J=1,NCHI
            WRITE(*,120) RCHI(J),RW1(J)
         ENDDO
 120     FORMAT(E16.8,2X,E16.8)
      ENDIF
      
C     COEFFICIENTS FOR XI^2 IN H-FACTOR
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RX2(J)=-RW2(J)*RJA(JS,J)*T(JS)/DPSIDS(JS)/B0K
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RX2(J)=-RW2(J)*RJAM(JS,J)*TM(JS)/DPSIDSM(JS)/B0K
         ENDDO
      ENDIF
      RX2(NCHI+1)=RX2(1)

C     COEFFICIENTS FOR Q^1 IN H-FACTOR
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RQ1(J)=G12L(JS,J)/RJA(JS,J)/B0K
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RQ1(J)=G12LM(JS,J)/RJAM(JS,J)/B0K
         ENDDO
      ENDIF
      RQ1(NCHI+1)=RQ1(1)

C     COEFFICIENTS FOR Q^2 IN H-FACTOR
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RQ2(J)=G22L(JS,J)/RJA(JS,J)/B0K
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RQ2(J)=G22LM(JS,J)/RJAM(JS,J)/B0K
         ENDDO
      ENDIF
      RQ2(NCHI+1)=RQ2(1)

C     COEFFICIENTS FOR Q^3 IN H-FACTOR
      IF (KGRID.EQ.1) THEN
         RQ3K=T(JS)/DPSIDS(JS)/B0K
      ELSEIF (KGRID.EQ.2) THEN
         RQ3K=TM(JS)/DPSIDSM(JS)/B0K
      ENDIF

C     COEFFICIENTS FOR COMPUTING BOUNCE TIME RTK
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RBT(J)=RJA(JS,J)*BK(JS,J,KGRID)*SQRT(HK(JS,J,KGRID))/
     &             DPSIDS(JS)
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RBT(J)=RJAM(JS,J)*BK(JS,J,KGRID)*SQRT(HK(JS,J,KGRID))/
     &             DPSIDSM(JS)
         ENDDO
      ENDIF
      RBT(NCHI+1)=RBT(1)

C     COEFFICIENTS FOR COMPUTING PHI(CHI)
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RPHI(J)=T(JS)*RJA(JS,J)/DPSIDS(JS)/REQ(JS,J)**2
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RPHI(J)=TM(JS)*RJAM(JS,J)/DPSIDSM(JS)/REQM(JS,J)**2
         ENDDO
      ENDIF
      RPHI(NCHI+1)=RPHI(1)

C     COEFFICIENT FOR DRIFT FREQUENCY 
C     METHOD I: ORIGNAL
C     CORRECT, INDEPENDENT OF COORDINATE SYSTEM
      DO J=2,NCHI-1
         RW2(J)=(BPK(JS,J+1,KGRID)-BPK(JS,J-1,KGRID))/2/RCHIH
      ENDDO
      RW2(1)=(BPK(JS,2,KGRID)-BPK(JS,NCHI,KGRID))/2/RCHIH
      RW2(NCHI)=(BPK(JS,1,KGRID)-BPK(JS,NCHI-1,KGRID))/2/RCHIH
      IF (KGRID.EQ.1) THEN
         IF (JS.LT.NRP1) THEN
            H1 = (CS(JS)-CS(JS-1))/2
            H2 = (CS(JS+1)-CS(JS))/2
            DO J=1,NCHI
               RW1(J)=(H1/H2*BPK(JS,J,2)-H2/H1*BPK(JS-1,J,2))/(H1+H2) -
     &                (H1-H2)*BPK(JS,J,1)/H1/H2
            ENDDO
         ELSEIF (JS.EQ.NRP1) THEN
         H1 = CS(JS)-CS(JS-1)
         DO J=1,NCHI
            RW1(J)=(BPK(JS-1,J,1)+3*BPK(JS,J,1)-4*BPK(JS-1,J,2))/H1
         ENDDO
         ENDIF
         DO J=1,NCHI
            RDB(J)=2.0*(G22L(JS,J)*RW1(J)/RJA(JS,J) - 
     &             G12L(JS,J)*RW2(J)/RJA(JS,J) -
     &             BPK(JS,J,1)*RJA(JS,J)/DPSIDS(JS)*
     &             (PPEQ(JS)+T(JS)*TP(JS)/REQ(JS,J)**2))
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         H1 = CS(JS+1) - CS(JS)
         DO J=1,NCHI
            RW1(J)=(BPK(JS+1,J,1)-BPK(JS,J,1))/H1
            RDB(J)=2.0*(G22LM(JS,J)*RW1(J)/RJAM(JS,J) - 
     &             G12LM(JS,J)*RW2(J)/RJAM(JS,J) -
     &             BPK(JS,J,2)*RJAM(JS,J)/DPSIDSM(JS)*
     &             (PPEQM(JS)+TM(JS)*TPM(JS)/REQM(JS,J)**2))
         ENDDO
      ENDIF
      RDB(NCHI+1)=RDB(1)
      
C     CHECK INTEGRANT FOR DB
      IF (KCHECK.EQ.1.AND.JS.EQ.JS0.AND.KGRID.EQ.1) THEN
         WRITE(*,*) 'CHECK KEQUIL:CHI DB_INTEGRAND'
         DO J=1,NCHI
            WRITE(*,110) RCHI(J),BPK(JS,J,1),RDB(J)*0.5
         ENDDO
 110     FORMAT(E16.8,2X,E16.8,2X,E16.8)
      ENDIF

 
      RETURN
      END

C=======================================================================
C COMPUTE EQUILIBRIUM QUANTITIES AT INTEGRATION POINTS ALONG CHI       =
C USING 1D CUBIC SPLINE                                                =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KEQUILK(JS,KGRID)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER JS,J,KGRID

C     H=B0/B
      DO J=1,NCHI
         RW1(J) = HK(JS,J,KGRID)
      ENDDO
      RW1(NCHI+1)=RW1(1)
      CALL SPLINE1D(RHK,RCHIK,NCHI2+2,RW1,RCHI,NCHI+1,RCHI2)

C     COEFFICIENT FOR G-FACTORS JB
      CALL SPLINE1D(RJBK,RCHIK,NCHI2+2,RJB,RCHI,NCHI+1,RCHI2)

C     COEFFICIENTS FOR XI^1 IN H-FACTOR
      CALL SPLINE1D(RX1PK,RCHIK,NCHI2+2,RX1P,RCHI,NCHI+1,RCHI2)

      CALL SPLINE1D(RX1BK,RCHIK,NCHI2+2,RX1B,RCHI,NCHI+1,RCHI2)

      CALL SPLINE1D(RX1RK,RCHIK,NCHI2+2,RX1R,RCHI,NCHI+1,RCHI2)
      
C     COEFFICIENTS FOR XI^2 IN H-FACTOR
      CALL SPLINE1D(RX2K,RCHIK,NCHI2+2,RX2,RCHI,NCHI+1,RCHI2)

C     COEFFICIENTS FOR Q^1 IN H-FACTOR
      CALL SPLINE1D(RQ1K,RCHIK,NCHI2+2,RQ1,RCHI,NCHI+1,RCHI2)

C     COEFFICIENTS FOR Q^2 IN H-FACTOR
      CALL SPLINE1D(RQ2K,RCHIK,NCHI2+2,RQ2,RCHI,NCHI+1,RCHI2)

      RETURN
      END

C=======================================================================
C COMPUTE PHI(CHI)=RPHIK AT INTEGRATION POINTS ALONG CHI-ANGLE         =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KPHI(JS,KGRID)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER JS,KGRID,J
      REAL*8  RTMP,PI2,DIFPI

      INTEGER KCHECK
      KCHECK = 0

      PI2 = PI + PI

C     DEFINE NEW INTEGRATION POINTS
      DO J=1,NCHI2+1
         IF (RCHIK(J).LT.RCHIK(J+1)) THEN
            RCHIN(2*J-1) = ((1+WK)*RCHIK(J)+(1-WK)*RCHIK(J+1))/2
            RCHIN(2*J)   = ((1-WK)*RCHIK(J)+(1+WK)*RCHIK(J+1))/2
         ELSE
            RCHIN(2*J-1) = ((1+WK)*RCHIK(J)+(1-WK)*(RCHIK(J+1)+PI2))/2
            RCHIN(2*J)   = ((1-WK)*RCHIK(J)+(1+WK)*(RCHIK(J+1)+PI2))/2
            IF (RCHIN(2*J-1).GT.PI2) RCHIN(2*J-1)=RCHIN(2*J-1)-PI2
            IF (RCHIN(2*J).GT.PI2)   RCHIN(2*J)  =RCHIN(2*J)-PI2
         ENDIF
      ENDDO

C     SPLINE EQUILIBRIUM QUANTITIES
      CALL SPLINE1D(RVALN,RCHIN,2*NCHI2+2,RPHI,RCHI,NCHI+1,RCHI2)
      
C     CUMULATIVE INTEGRATION TO COMPUTE RPHIK
      RTMP = 0.0
      RPHIK(1) = RTMP
      DO J=1,NCHI2+1
         RTMP = RTMP + (RVALN(2*J-1)+RVALN(2*J))*
     &                 DIFPI(RCHIK(J+1)-RCHIK(J))*0.5
         RPHIK(J+1) = RTMP
      ENDDO

      IF (KCHECK.EQ.1) THEN
         IF (JS.EQ.JS0.AND.KGRID.EQ.1) THEN
            WRITE(*,*) 'CHECK KPHI: RCHIK RPHIK'
            DO J=1,NCHI2+2
               WRITE(*,110) RCHIK(J),RPHIK(J)
            ENDDO
 110        FORMAT(E12.4,2X,E12.4)
         ENDIF
      ENDIF

      RETURN
      END

C=======================================================================
C COMPUTE T(CHI)=RTK AT INTEGRATION POINTS ALONG CHI-ANGLE     
C YQL, 08-2007                                                 
C ALSO COMPUTE:
C   TPSI0F=INT[J*SQRT(1-LAM/H)/DPSIDS]DCHI/(2*PI), ASSCOAITED WITH 
C          TPSI0DLAM FOR FOW CORRECTION OF THERMAL PASSING IONS
C   HPSI0F=INT[RBT/(H-LAM)**1.5]DCHI, ASSOCIATED WITH 
C          HPSI0DLAM FOR FOW CORRECTION OF HOT PASSING IONS
C   HATJ0F=INT[RBT*SQRT(H-LAM)/H]DCHI, ASSOCIATED WITH
C          LONGITUDINAL INVARIANT OF PARTICLES
C 
C YQL, 07-2013
C=======================================================================
      SUBROUTINE KBTIME(JS,KGRID,KPARTICLE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER JS,KGRID,KPARTICLE,J,J1,J2
      REAL*8  RTMP,RDPIDS,FL,FU,PI2,DIFPI

      INTEGER KCHECK
      
      REAL*8 DIFFERCHI
      
      PI2 = PI + PI

C     DEFINE NEW INTEGRATION POINTS
      DO J=1,NCHI2+1
         IF (RCHIK(J).LT.RCHIK(J+1)) THEN
            RCHIN(2*J-1) = ((1+WK)*RCHIK(J)+(1-WK)*RCHIK(J+1))/2
            RCHIN(2*J)   = ((1-WK)*RCHIK(J)+(1+WK)*RCHIK(J+1))/2
         ELSE
            RCHIN(2*J-1) = ((1+WK)*RCHIK(J)+(1-WK)*(RCHIK(J+1)+PI2))/2
            RCHIN(2*J)   = ((1-WK)*RCHIK(J)+(1+WK)*(RCHIK(J+1)+PI2))/2
            IF (RCHIN(2*J-1).GT.PI2) RCHIN(2*J-1)=RCHIN(2*J-1)-PI2
            IF (RCHIN(2*J).GT.PI2)   RCHIN(2*J)  =RCHIN(2*J)-PI2
         ENDIF
      ENDDO

C     SPLINE EQUILIBRIUM QUANTITIES
      DO J=1,NCHI
         RW1(J)=HK(JS,J,KGRID)
      ENDDO
      RW1(NCHI+1)=RW1(1)
      CALL SPLINE1D(RHN,RCHIN,2*NCHI2+2,RW1,RCHI,NCHI+1,RCHI2)
      
      CALL SPLINE1D(RVALN,RCHIN,2*NCHI2+2,RBT,RCHI,NCHI+1,RCHI2)
      
C     CUMULATIVE INTEGRATION TO COMPUTE RTK
      RTMP   = 0.
      TPSI0F = 0.
      HPSI0F = 0.
      HATJ0F = 0.
      RTK(1) = RTMP
      IF (KPARTICLE.EQ.1) THEN
         DO J=1,NCHI2+1
            J1 = 2*J - 1
            J2 = J1 + 1
            RTMP = RTMP + (RVALN(J1)/SQRT(RHN(J1)-LAM)+
     &             RVALN(J2)/SQRT(RHN(J2)-LAM))*(RCHIK(J+1)-RCHIK(J))/2.
            RTK(J+1) = RTMP
            TPSI0F   = TPSI0F + (RVALN(J1)*SQRT(RHN(J1)-LAM)+
     &                           RVALN(J2)*SQRT(RHN(J2)-LAM))
     &                          *(RCHIK(J+1)-RCHIK(J))/2.
            HPSI0F   = HPSI0F + (RVALN(J1)/(RHN(J1)-LAM)**1.5+
     &                           RVALN(J2)/(RHN(J2)-LAM)**1.5)
     &                          *(RCHIK(J+1)-RCHIK(J))/2.
            HATJ0F   = HATJ0F + (RVALN(J1)*SQRT(RHN(J1)-LAM)/RHN(J1)+
     &                           RVALN(J2)*SQRT(RHN(J2)-LAM)/RHN(J2))
     &                          *(RCHIK(J+1)-RCHIK(J))/2.
         ENDDO
         TPSI0F = TPSI0F/B0K/PI2
      ELSEIF (KPARTICLE.EQ.0) THEN
         FL = RJBK(1)*SQRT(LAM)/RDPSIK/SQRT(HPL)
         FU = RJBK(NCHI2+2)*SQRT(LAM)/RDPSIK/SQRT(-HPU)
         DO J=1,NCHI2+1
            J1 = 2*J - 1
            J2 = J1 + 1
            RTMP = RTMP + (RVALN(J1)/SQRT(RHN(J1)-LAM)-
     &             FL/SQRT(DIFPI(RCHIN(J1)-CHIL))-
     &             FU/SQRT(DIFPI(CHIU-RCHIN(J1)))+
     &             RVALN(J2)/SQRT(RHN(J2)-LAM)-
     &             FL/SQRT(DIFPI(RCHIN(J2)-CHIL))-
     &             FU/SQRT(DIFPI(CHIU-RCHIN(J2))))*
     &             DIFPI(RCHIK(J+1)-RCHIK(J))/2.
            RTK(J+1) = RTMP + FL*SQRT(DIFPI(RCHIK(J+1)-CHIL))*2. + 
     &                 FU*(SQRT(DIFFERCHI(CHIU,CHIL))-
     &                     SQRT(DIFPI(CHIU-RCHIK(J+1))))*2.
            HATJ0F   = HATJ0F + (RVALN(J1)*SQRT(RHN(J1)-LAM)/RHN(J1)+
     &                           RVALN(J2)*SQRT(RHN(J2)-LAM)/RHN(J2))*
     &                           DIFPI(RCHIK(J+1)-RCHIK(J))/2.
         ENDDO
         HATJ0F = 2.*HATJ0F
      ENDIF
         
C     COMPUTE S-FACTOR
      RSS = 1.0
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RW1(J)=RJA(JS,J)
         ENDDO
         RW1(NCHI+1)=RW1(1)
         CALL SPLINE1D(RVALN,RCHIN,2*NCHI2+2,RW1,RCHI,NCHI+1,RCHI2)

         RTMP = 0.0
         DO J=1,NCHI2+1
            J1 = 2*J - 1
            J2 = J1 + 1
            IF (KPARTICLE.EQ.1)
     &      RTMP = RTMP+(RVALN(J1)+RVALN(J2))*(RCHIK(J+1)-RCHIK(J))/2
            IF (KPARTICLE.EQ.0)
     &      RTMP = RTMP+(RVALN(J1)+RVALN(J2))*
     &                  DIFPI(RCHIK(J+1)-RCHIK(J))/2
         ENDDO
         RSS = RTMP
      ENDIF

C     COMPUTE U-FACTOR
      RUU  = 0.0
      RUU2 = 0.0
      IF (KGRID.EQ.1) THEN
      RTMP = 0.0
      IF (KPARTICLE.EQ.1) THEN
         DO J=1,NCHI2+1
            J1 = 2*J - 1
            J2 = J1 + 1
            RTMP = RTMP + (RVALN(J1)/SQRT(RHN(J1))/SQRT(RHN(J1)-LAM)+
     &                     RVALN(J2)/SQRT(RHN(J2))/SQRT(RHN(J2)-LAM))*
     &                    (RCHIK(J+1)-RCHIK(J))/2
            RUU2 = RUU2 + (SQRT(1.-LAM/RHN(J1))+SQRT(1.-LAM/RHN(J2)))*
     &                    (RCHIK(J+1)-RCHIK(J))/2
         ENDDO
         RUU  = RTMP
         RUU2 = RUU2/2./PI
      ELSEIF (KPARTICLE.EQ.0) THEN
         FL = RJBK(1)*SQRT(LAM)/B0K/SQRT(HPL)
         FU = RJBK(NCHI2+2)*SQRT(LAM)/B0K/SQRT(-HPU)
         DO J=1,NCHI2+1
            J1 = 2*J - 1
            J2 = J1 + 1
            RTMP = RTMP + (RVALN(J1)/SQRT(RHN(J1))/SQRT(RHN(J1)-LAM)-
     &             FL/SQRT(DIFPI(RCHIN(J1)-CHIL))-
     &             FU/SQRT(DIFPI(CHIU-RCHIN(J1)))+
     &             RVALN(J2)/SQRT(RHN(J2))/SQRT(RHN(J2)-LAM)-
     &             FL/SQRT(DIFPI(RCHIN(J2)-CHIL))-
     &             FU/SQRT(DIFPI(CHIU-RCHIN(J2))))*
     &             DIFPI(RCHIK(J+1)-RCHIK(J))/2
         ENDDO
         RUU = RTMP + (FL+FU)*2.0*SQRT(PI2+CHIU-CHIL)
      ENDIF
      ENDIF

      KCHECK = 0
      IF (KCHECK.EQ.1.AND.JS.EQ.JS0.AND.KGRID.EQ.1.AND.
     &   KPARTICLE.EQ.1) THEN
         IF (ABS(LAMM(J+1)-LAM).LT.1.0E-13) THEN
            WRITE(*,*) 'CHECK KBTIME: RCHIK RTK RX1PK RX1BK'
            DO J=1,NCHI2+2
               WRITE(*,110) RCHIK(J),RTK(J),RX1PK(J),RX1BK(J)
            ENDDO
 110        FORMAT(4E17.8)
         ENDIF
      ENDIF

      IF (KCHECK.EQ.2.AND.JS.EQ.JS0.AND.KGRID.EQ.1) THEN
         IF (ABS(LAMM(J+1)-LAM).LT.1.0E-13) 
     &      WRITE(*,*) 'CHECK KBTIME: LAM OMEGAB'            
         IF (KPARTICLE.EQ.1) RTMP=2.0*PI/RTK(NCHI2+2)
         IF (KPARTICLE.EQ.0) RTMP=PI/RTK(NCHI2+2)
         WRITE(*,110) LAM,RTMP
      ENDIF

      RETURN
      END

C=======================================================================
      FUNCTION DIFFERCHI (CHIU,CHIL) RESULT (RES)
      
      USE GLOBALM
      IMPLICIT NONE
      REAL*8,INTENT(IN):: CHIU,CHIL
      REAL*8 RES
C     MODIFIED FOR THE NON UP AND DOWN ASYMMETRY CASE         
      IF (CHIU-CHIL .GT. 0) THEN
        IF (CHIU-CHIL .LE.PI) THEN
            RES = CHIU - CHIL       
        ELSE
            STOP 'VALUE OF CHIU-CHIL MAYBE NOT PHYSICAL.'
        ENDIF
      ELSE
        RES = PI + PI - CHIL + CHIU
      ENDIF  
      
      END FUNCTION DIFFERCHI
      
C=======================================================================
C COMPUTE COEFFICIENT OF TOROIDAL DRIFT PRECESSION FREQUENCY FOR       =
C TRAPPED PARTICLES                                                    =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KDRIFT(JS,KGRID,KPARTICLE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER J,JS,KGRID,KPARTICLE
      REAL*8  DB,DMU,FL,FU,DIFPI,PHASE
      REAL*8  DIFFERCHI
      
      INTEGER KCHECK
      KCHECK = 0

C     FOR TRAPPED PARTICLES
      IF (KPARTICLE.EQ.0) THEN

C     FIRST COMPUTE DB, WHICH DOES NOT INVOLVE SINGULAR INTEGRATION
      CALL SPLINE1D(RVALK,RCHIK,NCHI2+2,RDB,RCHI,NCHI+1,RCHI2)
      
      DB = 0.0
      DO J=2,NCHI2+1
         DB = DB + RVALK(J)*SQRT(1-LAM/RHK(J))
      ENDDO
      DB = DB*RCHIHK/2.*OMEGAB/PI

C     THEN COMPUTE DMU, WHICH INVOLVES SINGULAR INTEGRATION
      CALL SPLINE1D(RVALK,RCHIK,NCHI2+2,RDMU,RCHI,NCHI+1,RCHI2)
      
      DMU = 0.0
      FL = RVALK(1)/SQRT(HPL)
      FU = RVALK(NCHI2+2)/SQRT(-HPU)
      DO J=2,NCHI2+1
         DMU = DMU + RVALK(J)/SQRT(RHK(J)-LAM) - 
     &         FL/SQRT(DIFPI(RCHIK(J)-CHIL)) - 
     &         FU/SQRT(DIFPI(CHIU-RCHIK(J)))
      ENDDO
      DMU = DMU*RCHIHK/2 + (FL+FU)*2.0*SQRT(DIFFERCHI(CHIU,CHIL))
      DMU = DMU*LAM*OMEGAB/PI

      DRIFT = DB + DMU

C     FOR PASSING PARTICLES
      ELSE

      DRIFT = 0.
      DO J=2,NCHI2+1
         PHASE = 1.0/SQRT(1.0-LAM/RHK(J))
         DRIFT = DRIFT + PHASE*((1-LAM/RHK(J))*RX1PK(J)+
     &                         (2-LAM/RHK(J))*RX1BK(J))
      ENDDO
      PHASE = RCHIHK/4*OMEGAB/PI
      IF (KGRID.EQ.1) PHASE = -PHASE/DPSIDS(JS)
      IF (KGRID.EQ.2) PHASE = -PHASE/DPSIDSM(JS)
      DRIFT = DRIFT*PHASE

      ENDIF

      IF (KCHECK.EQ.1.AND.ABS(RQK-Q(JS0)).LT.1.0E-13.AND.
     &    KPARTICLE.EQ.0) THEN
          LAMM(2) = ((1+WK)*LAMK0(JS,1,KGRID)+
     &               (1-WK)*LAMK0(JS,2,KGRID))/2.
         IF (ABS(LAMM(2)-LAM).LT.1.0E-13)
     &      WRITE(*,*) 'CHECK KDRIFT: LAM DMU DB DRIFT'
         WRITE(*,110) LAM,DMU,DB,DRIFT
 110     FORMAT(4(E16.8,1X))
      ENDIF

      RETURN
      END

C=======================================================================
C FIND LAMBDA WHERE ANALYTIC SUBTRACTION OF SINGULARITY IS REQUIRED FOR
C PARTICLE PITCH ANGLE INTEGRATION, USING 1D CUBIC SPLINE                                                
C   KPARTICLE=0: FOR BOUNCE RESONANCE OF TRAPPED PARTICLES
C   KPARTICLE=1: FOR TRANSIT/PRECESSION RESONANCE OF PASSING PARTICLES
C   KPARTICLE=2: FIND DRIFT(LAM0)=0 FOR TRAPPED PARTICLES
C YQL, 08-2013                                                         
C=======================================================================
      SUBROUTINE KLAM0(JS,KGRID,KPARTICLE)

      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      USE MPIENV
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KPARTICLE,NN,J,KP,L
      REAL*8     LAM00,LAM0,RTMP,RTMP1,RL

      INTEGER KCHECK
      KCHECK = 1

C     FIND SLAM0 FOR PASSING PARTICLES
      IF (KPARTICLE.EQ.1) THEN
         SLAM0 =-1.
         IF (ABS(WFUN(JS,KGRID)).GE.1.0E-13) THEN
         NN    = 2*NLAMK1(JS,KGRID)
         DO KP=1,NSPECIES
            IF (ABS(PSPECIES_NP(KP)).GT.0.AND.ISPECIES_EK(KP).EQ.1) THEN
               DO L=1,MLMAX
                  RTMP = ABS(RNTOR*OMEGAE0(JS,KGRID)-DREAL(OMEGA))/
     &                   ABS(RNTOR*RQK+RLM(L))/SQRT(2.*
     &                   EPSALPHA(JS,KGRID,KP)/ESPECIES_M(KP)*
     &                   ESPECIES_M(1))
                  J = 1
 100              IF ((ZOMEGABP(JS,J,KGRID)-RTMP)*
     &                (ZOMEGABP(JS,J+1,KGRID)-RTMP).GT.0.) THEN
                     J = J+1
                     IF (J.LT.NN) GOTO 100
                  ENDIF
                  IF (J.LT.NN-1) THEN
                     LAM00 = (LAMM(J)+LAMM(J+1))*0.5
                     CALL SPLINE1DR(RTMP1,LAM0,RTMP,LAM00,
     &                    ZOMEGABP(JS,:,KGRID),LAMM,NN-1,LAMTMP)
C                    DO NOT ADD AN ENDPOINT LOGARITHMIC SINGULARITY.
                     IF (J.NE.1.OR.ZOMEGABP(JS,1,KGRID).NE.RTMP)
     &                  SLAM0(L,KP) = LAM0
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         ENDIF
      ENDIF

C     FIND SLAM0 FOR TRAPPED PARTICLES 
C     BOUNCE FOR RL.NE.0 AND PRECESSION FOR RL=0 
      IF (KPARTICLE.EQ.0) THEN
         SLAM0 =-1.
         IF (ABS(WFUN(JS,KGRID)).GE.1.0E-13) THEN
         NN    = 2*NLAMK0(JS,KGRID)
         DO KP=1,NSPECIES
            DO L=1,MLMAX
               RL = RLM(L)
               IF (ABS(RL).GT.0.1.AND.ABS(PSPECIES_NTB(KP)).GT.0.AND.
     &             ISPECIES_EK(KP).EQ.1) THEN
                  RTMP =-(RNTOR*OMEGAE0(JS,KGRID)-DREAL(OMEGA))/
     &                   RL/SQRT(2.*EPSALPHA(JS,KGRID,KP)/
     &                   ESPECIES_M(KP)*ESPECIES_M(1))
                  J = 2
 110              IF ((ZOMEGABT(JS,J,KGRID)-RTMP)*
     &                (ZOMEGABT(JS,J+1,KGRID)-RTMP).GT.0.) THEN
                     J = J+1
                     IF (J.LT.NN) GOTO 110
                  ENDIF
                  IF (J.LT.NN) THEN
                     LAM00 = (LAMM(J)+LAMM(J+1))*0.5
                     CALL SPLINE1DR(RTMP1,LAM0,RTMP,LAM00,
     &                    ZOMEGABT(JS,2:NN,KGRID),LAMM(2),NN-1,LAMTMP)
                     IF (J.NE.NN-1.OR.ZOMEGABT(JS,NN,KGRID).NE.RTMP)
     &                  SLAM0(L,KP) = LAM0
                  ENDIF
               ENDIF

               IF (ABS(RL).LT.0.1.AND.ABS(PSPECIES_NTD(KP)).GT.0.AND.
     &             ISPECIES_EK(KP).EQ.1) THEN
                  RTMP =-(RNTOR*OMEGAE0(JS,KGRID)-DREAL(OMEGA))*
     &                   OMEGACI0*ESPECIES_Z(KP)/(RNTOR*B0K*
     &                   ESPECIES_Z(1)*EPSALPHA(JS,KGRID,KP))
                  J = 2
 120              IF ((ZOMEGADT(JS,J,KGRID)-RTMP)*
     &                (ZOMEGADT(JS,J+1,KGRID)-RTMP).GT.0.) THEN
                     J = J+1
                     IF (J.LT.NN) GOTO 120
                  ENDIF
                  IF (J.LT.NN) THEN
                     LAM00 = (LAMM(J)+LAMM(J+1))*0.5
                     CALL SPLINE1DR(RTMP1,LAM0,RTMP/1000.,LAM00,
     &               ZOMEGADT(JS,2:NN,KGRID)/1000.,LAMM(2),NN-1,LAMTMP)
                     IF (J.NE.NN-1.OR.ZOMEGADT(JS,NN,KGRID).NE.RTMP)
     &                  SLAM0(L,KP) = LAM0
                  ENDIF

                  IF (JS.EQ.JS0.AND.KGRID.EQ.1.AND.KCHECK.EQ.2) THEN
                  CALL MPI_OPEN_FILE(RANK)
                  WRITE(CHMPI,*) 'RTMP=',RTMP
                  DO J=1,NN
                     WRITE(CHMPI,*) LAMM(J),ZOMEGADT(JS,J,KGRID)
                  ENDDO
                  CALL MPI_CLOSE_FILE(RANK)
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
         ENDIF
      ENDIF

C     FIND SLAMD0 FOR PRECESSIONAL RESONANCE OF TRAPPED PARTICLES
      IF (KPARTICLE.EQ.2) THEN
         SLAMD0 =-1.
         IF (ABS(WFUN(JS,KGRID)).GE.1.0E-13) THEN
         NN     = 2*NLAMK0(JS,KGRID)

C        THE SPLINE USES LAMM(2:NN), INCLUDING ITS END INTERVALS.
         J=2
 130     IF (ZOMEGADT(JS,J,KGRID)*ZOMEGADT(JS,J+1,KGRID).GT.0.0) THEN
            J=J+1
            IF (J.LT.NN) GOTO 130
         ENDIF
         IF (J.LT.NN) THEN
            LAM00 = (LAMM(J)+LAMM(J+1))*0.5
            CALL SPLINE1DR(RTMP1,LAM0,0.0,LAM00,
     &           ZOMEGADT(JS,2:NN,KGRID)/1000.,LAMM(2),NN-1,LAMTMP)
            IF (J.NE.NN-1.OR.ZOMEGADT(JS,NN,KGRID).NE.0.0) THEN
               SLAMD0 = LAM0
               DPRM   = RTMP1*1000.
            ENDIF
         ENDIF
         ENDIF
      ENDIF

      IF (KCHECK.EQ.1.AND.JS.EQ.JS0.AND.KGRID.EQ.1) THEN
         CALL MPI_OPEN_FILE(RANK)
         IF (KPARTICLE.EQ.2) THEN
            WRITE(CHMPI,*) 'KLAM0:SLAMD0,DPRM=',SLAMD0,DPRM
         ELSEIF (KPARTICLE.EQ.0) THEN
            WRITE(CHMPI,*) 'KLAM0:KPARTICLE=0: KP,RLM(L),SLAM0'
            DO KP=1,NSPECIES
               IF ((ABS(PSPECIES_NTD(KP)).GT.0.OR.ABS(PSPECIES_NTB(KP)) 
     &             .GT.0).AND.ISPECIES_EK(KP).EQ.1) THEN
                  DO L=1,MLMAX
                     WRITE(CHMPI,200) KP,RLM(L),SLAM0(L,KP)
                  ENDDO
               ENDIF
            ENDDO
         ELSEIF (KPARTICLE.EQ.1) THEN
            WRITE(CHMPI,*) 'KLAM0:KPARTICLE=1'
            DO KP=1,NSPECIES
               IF (ABS(PSPECIES_NP(KP)).GT.0.AND.ISPECIES_EK(KP).EQ.1) 
     &            THEN
                  DO L=1,MLMAX
                     WRITE(CHMPI,200) KP,RLM(L),SLAM0(L,KP)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         CALL MPI_CLOSE_FILE(RANK)
      ENDIF
 200  FORMAT(I2,1X,E10.2,1X,E13.5)

      RETURN
      END

C=======================================================================
C ENERGY INTEGRATION FACTOR FOR PRECESSION DRIFT RESONANCE             
C INPUT:
C   JS: SURFACE                                                
C   KGRID=1: INTEGER RADIAL POINT; 
C         2: HALF-INTEGER POINT             
C   KPARTICLE=0: TRAPPED; 1: PASSING  (ONLY TRAPPED CONSIDERED FOR NOW)
C   KP: # IN ISPECIES_F0 LIST 
C   KOPT=0: ZERO ORBIT WIDTH CONTRIBUTION
C        2: FOW CORRECTION 2 FOR TRAPPED PARTICLES
C        3: FOW CORRECTION 3 FOR TRAPPED PARTICLES
C IMPLICITLY ASSUMED INPUT PARAMETERS (VIA KINETICM.F)
C   DRIFT
C   OMEGAB
C   LAM
C OUTPUT:
C   ZVI: COMPUTED ENERGY INTEGRATION VALUE
C NUMERICAL ENERGY INTEGRATION ADDED BY Z.R.WANG
C YQL, 06-2013                                                         
C=======================================================================
      SUBROUTINE KI_PRECESSION(JS,KGRID,KPARTICLE,KP,KOPT,ZVI)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      USE ToolBox
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KPARTICLE,KP,KOPT,L,IF0TYPE,J

      REAL*8     OMEGAE,RTMP,REPS,SIG
      COMPLEX*16 ZVI(4),OMEGAN,OMEGASA,OMEGASB,OMEGASC,
     &           plasmax,plasmaf
      PARAMETER  (REPS=1.0E-13)

      COMPLEX*16 NUMER,DENOM
      INTEGER    COUNTER
      REAL*8     TMP

      INTEGER KCHECK
      KCHECK=0

      IF0TYPE = ISPECIES_F0(KP)
      OMEGAE  = OMEGAE0(JS,KGRID)

      ZVI = 0.

      IF (KPARTICLE.EQ.1) RETURN

C     PRECESSION DRIFT OF TRAPPED PARTICLES WITH IF0TYPE=0
      IF (IF0TYPE.EQ.0.AND.KOPT.EQ.0.AND.KFASTRUN.EQ.1) THEN  
      RTMP = RNTOR*DRIFT*ESPECIES_TEM(JS,KGRID,KP)*B0K/OMEGACI0
     &       *ESPECIES_Z(1)/ESPECIES_Z(KP)

      IF (INUTYPE.EQ.0) THEN

      IF (ABS(RTMP).LT.REPS) THEN
         OMEGAN  = RNTOR*OMEGAE - OMEGA -CI*NUEFF(JS,KGRID,KP)
         OMEGASA = RNTOR*(OMEGASN(JS,KGRID,KP)-1.5*OMEGAST(JS,KGRID,KP))
     &             +CI*NUEFF(JS,KGRID,KP)
         OMEGASB = RNTOR*OMEGAST(JS,KGRID,KP)
         
         IF (ABS(OMEGAN).LT.REPS) OMEGAN=REPS
         ZVI(1) = 15.0/4.0*SQRT(PI)*(1+OMEGASA/OMEGAN) + 
     &            105.0/8.0*SQRT(PI)*OMEGASB/OMEGAN
      ELSE
         OMEGAN  =(RNTOR*OMEGAE - OMEGA -CI*NUEFF(JS,KGRID,KP))/RTMP
         OMEGASA =(RNTOR*(OMEGASN(JS,KGRID,KP)-1.5*OMEGAST(JS,KGRID,KP))
     &            +CI*NUEFF(JS,KGRID,KP))/RTMP
         OMEGASB = RNTOR*OMEGAST(JS,KGRID,KP)/RTMP

         IF (ABS(OMEGAN).LT.1000.0) THEN
         ZVI(1)=15.0*SQRT(PI)/4.0*OMEGASB+4.0*SQRT(PI)*(OMEGAN+OMEGASA-
     &          OMEGAN*OMEGASB)*(3.0/8.0-OMEGAN/4.0+OMEGAN**2/2.0-
     &          OMEGAN**2*plasmax(OMEGAN,-RTMP))
         ELSE
         ZVI(1)=15.0*SQRT(PI)/4.0*OMEGASB+4.0*SQRT(PI)*(OMEGAN+OMEGASA-
     &          OMEGAN*OMEGASB)*(15.0/16.0-105.0/32.0/OMEGAN)/OMEGAN
         ENDIF
      ENDIF

      ELSE

      IF (.NOT.ISLSODE) THEN
      DO COUNTER = 1, intMeshPot
      
         NUMER = RNTOR*(OMEGASN(JS,KGRID,KP) + (intMeshX(COUNTER) 
     &   - 1.5)*OMEGAST(JS,KGRID,KP) + OMEGAE) - OMEGA
         DENOM = RNTOR*OMEGAE + RTMP*intMeshX(COUNTER) - DREAL(OMEGA)
        
         IF (INUTYPE.GE.1) THEN
            TMP   = intMeshX(COUNTER)**1.5
         ELSEIF (INUTYPE.EQ.-1)THEN
            TMP   = 1.0
         ELSE
            STOP 'NO THIS COLLISIONALITY OPTION.'
         ENDIF
        
         DENOM = TMP*DENOM - CI*(NUEFF(JS,KGRID,KP)+DIMAG(OMEGA)*TMP)
      
         ZVI(1) = ZVI(1) + NUMER/DENOM * intMeshWeight(COUNTER) * 2.0
      ENDDO
      ELSE
         RTMP = RNTOR*DRIFT*ESPECIES_TEM(JS,KGRID,KP)*B0K/OMEGACI0
     &       *ESPECIES_Z(1)/ESPECIES_Z(KP)
         CALL ADAPTIVE_ENERGY_INTEGRAL_LSODE (ZVI(1),KPARTICLE,INUTYPE,
     &   RNTOR,OMEGASN(JS,KGRID,KP),OMEGAST(JS,KGRID,KP),OMEGAE,OMEGA,
     &   0.0D0,RTMP,NUEFF(JS,KGRID,KP))
      ENDIF
      ENDIF

      ZVI(1) = ZVI(1)*ESPECIES_PRE(JS,KGRID,KP)

C     SUBTRACT SINGULAR CONTRIBUTION FOR LATER PITCH ANGLE INTEGRATION
      DO L=1,MLMAX
         IF (ABS(RLM(L)).LT.0.1.AND.SLAM0(L,KP).GT.0.) 
     &   ZVI(1) = ZVI(1) - LOG(ABS(LAM-SLAM0(L,KP)))*SF0(L,KP,KOPT,1)
      ENDDO

C     PRECESSION DRIFT OF TRAPPED FAST IONS: IF0TYPE=1,2
      ELSEIF ((IF0TYPE.EQ.1.OR.IF0TYPE.EQ.2).AND.KOPT.EQ.0.AND.
     &        KFASTRUN.EQ.1) THEN   
      RTMP = RNTOR*DRIFT*EPSALPHA(JS,KGRID,KP)*B0K/OMEGACI0
     &       *ALPHAA3(JS,KGRID,KP)*ESPECIES_Z(1)/ESPECIES_Z(KP)

      IF (ABS(RTMP).LT.REPS) THEN
         OMEGAN  = ALPHAA3(JS,KGRID,KP)*(RNTOR*OMEGAE-OMEGA)
         OMEGASA = ALPHAA3(JS,KGRID,KP)*(RNTOR*OMEGAE-OMEGA)
         OMEGASB = RNTOR*(OMEGASNA(JS,KGRID,KP)+OMEGASAA(JS,KGRID,KP)) 
         TMP     = B0K/OMEGACI0*ESPECIES_TEM(JS,KGRID,KP)
     &             *ESPECIES_Z(1)/ESPECIES_Z(KP)
         OMEGASC = RNTOR*(OMEGASCA(JS,KGRID,KP)+
     &             TMP*1.5*EPSLONCA(JS,KGRID,KP)**1.5*
     &             DEPSALPHADPSI(JS,KGRID,KP))
         
         IF (ABS(OMEGAN).LT.REPS) OMEGAN=REPS
         ZVI(1) = plasmaf(OMEGAN,OMEGASA,OMEGASB,OMEGASC,
     &                    EPSLONCA(JS,KGRID,KP),0,1.0)
      ELSE
         OMEGAN  = ALPHAA3(JS,KGRID,KP)*(RNTOR*OMEGAE-OMEGA)/RTMP
         OMEGASA = ALPHAA3(JS,KGRID,KP)*(RNTOR*OMEGAE-OMEGA)/RTMP
         OMEGASB = RNTOR*(OMEGASNA(JS,KGRID,KP)+OMEGASAA(JS,KGRID,KP))/
     &             RTMP 
         TMP     = B0K/OMEGACI0*ESPECIES_TEM(JS,KGRID,KP)
     &             *ESPECIES_Z(1)/ESPECIES_Z(KP)
         OMEGASC = RNTOR*(OMEGASCA(JS,KGRID,KP)+
     &             TMP*1.5*EPSLONCA(JS,KGRID,KP)**1.5*
     &             DEPSALPHADPSI(JS,KGRID,KP))/RTMP
         SIG = RTMP/ALPHAA3(JS,KGRID,KP)
         SIG = SIG/ABS(SIG)
         ZVI(1) = plasmaf(OMEGAN,OMEGASA,OMEGASB,OMEGASC,
     &                   EPSLONCA(JS,KGRID,KP),1,SIG)

      ENDIF
      ZVI(1) = ZVI(1)*ESPECIES_PRE(JS,KGRID,KP)/ALPHAA3(JS,KGRID,KP)
     &         *AAK(JS,KGRID,KP)*(2.*PI)**1.5

C     SUBTRACT SINGULAR CONTRIBUTION FOR LATER PITCH ANGLE INTEGRATION
      DO L=1,MLMAX
         IF (ABS(RLM(L)).LT.0.1.AND.SLAM0(L,KP).GT.0.) 
     &   ZVI(1) = ZVI(1) - LOG(ABS(LAM-SLAM0(L,KP)))*SF0(L,KP,KOPT,1)
      ENDDO

C     ALL OTHER CASES USE NUMERICAL INTEGRATION ALONG PARTICLE ENERGY
      ELSE
         DO L=1,MLMAX
            IF (ABS(RLM(L)).LT.0.1) THEN 
               CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,0,ZVI(1))
               IF (INCDPHI.GT.0) THEN
                  CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,1,ZVI(2))
                  CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,2,ZVI(3))
                  CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,3,ZVI(4))
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         J = INT(NLAMK0(JS,KGRID)/2.)
         IF (ABS(LAMM(J+1)-LAM).LT.1.0E-13) THEN
            WRITE(*,121) KP,KOPT,0.0,LAM,ZVI 
         ENDIF
      ENDIF

      IF (KCHECK.EQ.2.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         WRITE(*,121) KP,KOPT,0.0,LAM,ZVI 
      ENDIF

 121  FORMAT('KI_PRECESSION:KP,KOPT,RL,LAM,ZVI:',2(I2,1X),
     &       7(E13.5,1X))

      RETURN
      END

C=======================================================================
C ENERGY INTEGRATION FACTOR FOR BOUNCE RESONANCE OF TRAPPED PARTICLES             
C INPUT:
C   JS: SURFACE                                                
C   KGRID=1: INTEGER RADIAL POINT; 
C         2: HALF-INTEGER POINT             
C   RL: ORBIT BOUNCE HARMONIC 
C   KP: # IN ISPECIES_F0 LIST 
C IMPLCITLY ASSUMED INPUT PARAMETERS (VIA KINETICM.F)
C   DRIFT
C   OMEGAB
C   LAM
C OUTPUT:
C   ZVI: COMPUTED ENERGY INTEGRATION VALUE
C NUMERICAL ENERGY INTEGRATION ADDED BY Z.R.WANG
C YQL, 06-2013                                                         
C THE CURRENT IMPLEMENTATION ONLY FOR KOPT=0
C=======================================================================
      SUBROUTINE KI_BOUNCE(JS,KGRID,KP,KOPT,L,ZVI)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ToolBox
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KP,KOPT,L,IF0TYPE,J

      REAL*8     RL,OMEGAE,RTMP,REPS
      COMPLEX*16 ZVI(4),OMEGAN,OMEGASA,OMEGASB,OMEGASC,
     &           plasmay
      PARAMETER  (REPS=1.0E-13)

      COMPLEX*16 NUMER,DENOM
      INTEGER    COUNTER
      REAL*8     TMP,RTMP1,TNUEFF

      INTEGER KCHECK
      KCHECK=0

      IF0TYPE = ISPECIES_F0(KP)

      OMEGAE = OMEGAE0(JS,KGRID)
      ZVI    = 0.
      RL     = RLM(L)

C     BOUNCE RESONANCE OF TRAPPED THERMAL PARTICLES WITH IF0TYPE=0
      IF (IF0TYPE.EQ.0.AND.KOPT.EQ.0.AND.KFASTRUN.EQ.1) THEN
      RTMP  = RL*OMEGAB*SQRT(2*ESPECIES_TEM(JS,KGRID,KP)*
     &        ESPECIES_M(1)/ESPECIES_M(KP))

      IF (INUTYPE.EQ.0) THEN

      IF (ABS(RTMP).LT.REPS) THEN
         OMEGAN  = RNTOR*OMEGAE - OMEGA -CI*NUEFF(JS,KGRID,KP)
         OMEGASA = RNTOR*(OMEGASN(JS,KGRID,KP)-1.5*OMEGAST(JS,KGRID,KP))
     &             +CI*NUEFF(JS,KGRID,KP)
         OMEGASB = RNTOR*OMEGAST(JS,KGRID,KP)
         
         IF (ABS(OMEGAN).LT.REPS) OMEGAN=REPS
         ZVI(1) = 15.0/4.0*SQRT(PI)*(1+OMEGASA/OMEGAN) + 
     &            105.0/8.0*SQRT(PI)*OMEGASB/OMEGAN
      ELSE
         OMEGAN  =(RNTOR*OMEGAE - OMEGA -CI*NUEFF(JS,KGRID,KP))/RTMP
         OMEGASA =(RNTOR*(OMEGASN(JS,KGRID,KP)-1.5*OMEGAST(JS,KGRID,KP))
     &            +CI*NUEFF(JS,KGRID,KP))/RTMP
         OMEGASB = RNTOR*OMEGAST(JS,KGRID,KP)/RTMP

         IF (ABS(OMEGAN).LT.100.0) THEN
         ZVI(1) = 4.0*(OMEGAN+OMEGASA+OMEGAN**2*OMEGASB)*
     &         (1.-3.*SQRT(PI)/8.*OMEGAN+0.5*OMEGAN**2-SQRT(PI)/4.*
     &         OMEGAN**3+0.5*OMEGAN**4-SQRT(PI)/2.*OMEGAN**5+
     &         OMEGAN**5*plasmay(OMEGAN)) + 
     &         4.*OMEGASB*(3.-15./16.*SQRT(PI)*OMEGAN)
         ELSE
         ZVI(1) = 4.0*(OMEGAN+OMEGASA+OMEGAN**2*OMEGASB)*
     &            (15./16.*SQRT(PI)-3./OMEGAN+
     &            105./32.*SQRT(PI)/OMEGAN**2)/OMEGAN+
     &            4.*OMEGASB*(3.-15./16.*SQRT(PI)*OMEGAN)
         ENDIF
      ENDIF

      ELSE

C     ENERGY INTERATION INCLDUING THE ENERGY DEPENDENT COLLISIONALITY          
C     THE PRECESSION DRIFT IS ALSO INCLUDED IN THE BOUNCE RESONANT OPERATOR     
C     PRECESSION DRIFT OF TRAPPED PARTICLES
      RTMP1 = RNTOR*DRIFT*ESPECIES_TEM(JS,KGRID,KP)*B0K/OMEGACI0
     &       *ESPECIES_Z(1)/ESPECIES_Z(KP)
     
      IF (.NOT.ISLSODE) THEN
 
      DO COUNTER = 1, intMeshPot
      
        NUMER = RNTOR*(OMEGASN(JS,KGRID,KP) + (intMeshX(COUNTER) 
     &  - 1.5)*OMEGAST(JS,KGRID,KP) + OMEGAE) - OMEGA
     
        DENOM = RNTOR*OMEGAE + RTMP*DSQRT(intMeshX(COUNTER))
     &  - DREAL(OMEGA)
        
        
        IF (INUTYPE.EQ.1) THEN
        
        DENOM = DENOM + RTMP1*intMeshX(COUNTER)
        TMP   = intMeshX(COUNTER)**1.5
        DENOM = TMP*DENOM - CI*(NUEFF(JS,KGRID,KP)*(1+(0.5*RL)*(0.5*RL))
     $        + DIMAG(OMEGA)*TMP)
        
        ELSEIF (INUTYPE.EQ.2) THEN
        
        DENOM = DENOM + RTMP1*intMeshX(COUNTER)
        TMP   = intMeshX(COUNTER)**1.5
        DENOM = TMP*DENOM - CI*(NUEFF(JS,KGRID,KP) + DIMAG(OMEGA)*TMP)
        
        ELSEIF (INUTYPE.EQ.-1) THEN
        
        TMP   = 1.0
        DENOM = TMP*DENOM - CI*(NUEFF(JS,KGRID,KP) + DIMAG(OMEGA)*TMP)
        
        ELSE
            STOP 'NO THIS COLLISIONALITY'
        ENDIF
        
        ZVI(1) = ZVI(1) + NUMER/DENOM * intMeshWeight(COUNTER) * 2.0
      ENDDO
      ELSE
         TNUEFF = NUEFF(JS,KGRID,KP)
         IF (INUTYPE.EQ.1) TNUEFF = TNUEFF * (1+(0.5*RL)*(0.5*RL))
         IF (INUTYPE.EQ.-1) RTMP1 = 0.0
         CALL ADAPTIVE_ENERGY_INTEGRAL_LSODE (ZVI(1),0,INUTYPE,
     &   RNTOR,OMEGASN(JS,KGRID,KP),OMEGAST(JS,KGRID,KP),OMEGAE,OMEGA,
     &   RTMP,RTMP1,TNUEFF)      
      ENDIF
      ENDIF

      ZVI(1) = ZVI(1)*ESPECIES_PRE(JS,KGRID,KP)

C     SUBTRACT SINGULAR CONTRIBUTION FOR LATER PITCH ANGLE INTEGRATION
      IF (SLAM0(L,KP).GT.0.) 
     &ZVI(1) = ZVI(1) - LOG(ABS(LAM-SLAM0(L,KP)))*SF0(L,KP,KOPT,1)

C     ALL OTHER CASES USE NUMERICAL INTEGRATION ALONG PARTICLE ENERGY
      ELSE
         CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,0,ZVI(1))
         IF (INCDPHI.GT.0) THEN
            CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,1,ZVI(2))
            CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,2,ZVI(3))
            CALL KIA_TRAP(JS,KGRID,KP,KOPT,L,LAM,3,ZVI(4))
         ENDIF
      ENDIF

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         J = INT(NLAMK0(JS,KGRID)/2.)
         IF (ABS(LAMM(J+1)-LAM).LT.1.0E-13) THEN
            WRITE(*,121) KP,KOPT,RL,LAM,ZVI
         ENDIF
      ENDIF

      IF (KCHECK.EQ.2.AND.KGRID.EQ.1.AND.JS.EQ.JS0.AND.RL.EQ.1.) THEN
         WRITE(*,121) KP,KOPT,RL,LAM,ZVI
      ENDIF

 121  FORMAT('KI_BOUNCE:KP,KOPT,RL,LAM,ZVI:',2(I2,1X),
     &       7(E13.5,1X))

      RETURN
      END

C=======================================================================
C ENERGY INTEGRATION FACTOR FOR TRANSIT RESONANCE OF PASSING PARTICLES             
C INPUT:
C   JS: SURFACE                                                
C   KGRID=1: INTEGER RADIAL POINT; 
C         2: HALF-INTEGER POINT             
C   RL: ORBIT BOUNCE HARMONIC 
C   KP: # IN ISPECIES_F0 LIST 
C IMPLCITLY ASSUMED INPUT PARAMETERS (VIA KINETICM.F)
C   DRIFT
C   OMEGAB
C   LAM
C OUTPUT:
C   ZVI: COMPUTED ENERGY INTEGRATION VALUE
C NUMERICAL ENERGY INTEGRATION ADDED BY Z.R.WANG
C YQL, 06-2013                    
C THE CURRENT IMPLEMENTATION VALID ONLY FOR KOPT=0                                     
C=======================================================================
      SUBROUTINE KI_TRANSIT(JS,KGRID,KP,KOPT,L,ZVI)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ToolBox
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KP,KOPT,L,IF0TYPE,J

      REAL*8     RL,OMEGAE,RTMP,REPS,RSIGN
      COMPLEX*16 ZVI(4),OMEGAN,OMEGASA,OMEGASB,OMEGASC,
     &           plasmaz
      PARAMETER  (REPS=1.0E-13)

      COMPLEX*16 NUMER,DENOM
      INTEGER    COUNTER
      REAL*8     TMP,TNUEFF

      INTEGER KCHECK
      KCHECK=0

      IF0TYPE = ISPECIES_F0(KP)

      OMEGAE = OMEGAE0(JS,KGRID)
      ZVI    = 0.
      RL     = RLM(L)

C     TRANSIT RESONANCE OF PASSING THERMAL PARTICLES WITH IF0TYPE=0
      IF (IF0TYPE.EQ.0.AND.KOPT.EQ.0.AND.KFASTRUN.EQ.1) THEN
      RTMP  = (RNTOR*RQK+RL)*OMEGAB*SQRT(2*ESPECIES_TEM(JS,KGRID,KP)*
     &        ESPECIES_M(1)/ESPECIES_M(KP))

      IF (INUTYPE.EQ.0) THEN

      IF (ABS(RTMP).LT.REPS) THEN
         OMEGAN  = RNTOR*OMEGAE - OMEGA -CI*NUEFF(JS,KGRID,KP)
         OMEGASA = RNTOR*(OMEGASN(JS,KGRID,KP)-1.5*OMEGAST(JS,KGRID,KP))
     &             +CI*NUEFF(JS,KGRID,KP)
         OMEGASB = RNTOR*OMEGAST(JS,KGRID,KP)
         
         IF (ABS(OMEGAN).LT.REPS) OMEGAN=REPS
         ZVI(1) = 15.0/4.0*SQRT(PI)*(1+OMEGASA/OMEGAN) + 
     &            105.0/8.0*SQRT(PI)*OMEGASB/OMEGAN
      ELSE
         OMEGAN  =(RNTOR*OMEGAE - OMEGA -CI*NUEFF(JS,KGRID,KP))/RTMP
         OMEGASA =(RNTOR*(OMEGASN(JS,KGRID,KP)-1.5*OMEGAST(JS,KGRID,KP))
     &            +CI*NUEFF(JS,KGRID,KP))/RTMP
         OMEGASB = RNTOR*OMEGAST(JS,KGRID,KP)/RTMP

         IF (ABS(OMEGAN).LT.100.0) THEN
         RSIGN = 1.
         IF (RTMP.GT.0.) RSIGN=-1.
         ZVI(1) = -15./4.*SQRT(PI)*OMEGAN*OMEGASB - 2.*SQRT(PI)*OMEGAN*
     &            (OMEGAN+OMEGASA+OMEGAN**2*OMEGASB)*(0.75+
     &            0.5*OMEGAN**2+OMEGAN**4+OMEGAN**5*
     &            RSIGN*plasmaz(OMEGAN*RSIGN))
         ELSE
         ZVI(1) = -15./4.*SQRT(PI)*OMEGAN*OMEGASB + 2.*SQRT(PI)*OMEGAN*
     &            (OMEGAN+OMEGASA+OMEGAN**2*OMEGASB)*(15./8.+
     &            105./16./OMEGAN**2)/OMEGAN**2
         ENDIF
      ENDIF

      ELSE

      IF (.NOT.ISLSODE) THEN
      
      DO COUNTER = 1, intMeshPot
        NUMER = RNTOR*(OMEGASN(JS,KGRID,KP) + (intMeshX(COUNTER) 
     &  - 1.5)*OMEGAST(JS,KGRID,KP) + OMEGAE) - OMEGA
     
        DENOM = RNTOR*OMEGAE - DREAL(OMEGA)
        
        IF (INUTYPE.EQ.1) THEN
        TMP   = intMeshX(COUNTER)**1.5
        DENOM = TMP*DENOM - CI*(NUEFF(JS,KGRID,KP)*(1+(0.5*RL)*(0.5*RL))
     $  +DIMAG(OMEGA)*TMP)
     
        ZVI(1) = ZVI(1) + NUMER * DENOM * 2.0
     $           /(DENOM*DENOM-RTMP*RTMP*intMeshX(COUNTER)**4.0)
     $           *intMeshWeight(COUNTER)
        ELSEIF (INUTYPE.EQ.2) THEN
        TMP   = intMeshX(COUNTER)**1.5
        DENOM = TMP*DENOM - CI*(NUEFF(JS,KGRID,KP) + DIMAG(OMEGA)*TMP)
     
        ZVI(1) = ZVI(1) + NUMER * DENOM * 2.0
     $           /(DENOM*DENOM-RTMP*RTMP*intMeshX(COUNTER)**4.0)
     $           *intMeshWeight(COUNTER)
           
        ELSEIF (INUTYPE.EQ.-1) THEN
        TMP   = 1.0
        DENOM = TMP*DENOM - CI*(NUEFF(JS,KGRID,KP) + DIMAG(OMEGA)*TMP)
        
        ZVI(1) = ZVI(1) + NUMER * DENOM * 2.0
     $           /(DENOM*DENOM-RTMP*RTMP*intMeshX(COUNTER))
     $           *intMeshWeight(COUNTER)   
        ELSE
            STOP 'NO THIS COLLISIONALITY'
        ENDIF
      ENDDO      
      ELSE
         TNUEFF = NUEFF(JS,KGRID,KP)
         IF (INUTYPE.EQ.1) TNUEFF = TNUEFF * (1+(0.5*RL)*(0.5*RL))
         CALL ADAPTIVE_ENERGY_INTEGRAL_LSODE (ZVI(1),1,INUTYPE,
     &   RNTOR,OMEGASN(JS,KGRID,KP),OMEGAST(JS,KGRID,KP),OMEGAE,OMEGA,
     &   RTMP,0.0D0,TNUEFF) 
      ENDIF
      ENDIF

      ZVI(1) = ZVI(1)*ESPECIES_PRE(JS,KGRID,KP)

C     SUBTRACT SINGULAR CONTRIBUTION FOR LATER PITCH ANGLE INTEGRATION
      IF (SLAM0(L,KP).GT.0.) 
     &ZVI(1) = ZVI(1) - LOG(ABS(LAM-SLAM0(L,KP)))*SF0(L,KP,KOPT,1)

C     ALL OTHER CASES USE NUMERICAL INTEGRATION ALONG PARTICLE ENERGY
      ELSE
         CALL KIA_PASS(JS,KGRID,KP,KOPT,L,LAM,0,ZVI(1))
         IF (INCDPHI.GT.0) THEN
            CALL KIA_PASS(JS,KGRID,KP,KOPT,L,LAM,1,ZVI(2))
            CALL KIA_PASS(JS,KGRID,KP,KOPT,L,LAM,2,ZVI(3))
            CALL KIA_PASS(JS,KGRID,KP,KOPT,L,LAM,3,ZVI(4))
         ENDIF
      ENDIF

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         J = INT(NLAMK1(JS,KGRID)/2.)
         IF (ABS(LAMM(J+1)-LAM).LT.1.0E-13) THEN
            WRITE(*,121) KP,KOPT,RL,LAM,ZVI
         ENDIF
      ENDIF

      IF (KCHECK.EQ.2.AND.KGRID.EQ.1.AND.JS.EQ.JS0.AND.RL.EQ.1.) THEN
         WRITE(*,121) KP,KOPT,RL,LAM,ZVI
      ENDIF

 121  FORMAT('KI_TRANSIT:KP,KOPT,RL,LAM,ZVI:',2(I2,1X),E8.1,1X,
     &       5(E13.5,1X))

      RETURN
      END

C=======================================================================
C PARTICLE ENERGY INTEGRATION FROM ALL SPECIES
C KOPT = 0: VANISHING ORBIT WIDTH CONTRIBUTION
C        1: FOW CORRECTION 1 
C        2: FOW CORRECTION 2 
C        3: FOW CORRECTION 3 
C OUTPUT:
C   VI: AN ARRAY FROM ALL BOUNCE HARMONICS 
C YQL, 06-2013 
C=======================================================================
      SUBROUTINE KI(JS,KGRID,KPARTICLE,KOPT)

      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER    JS,KGRID,KPARTICLE,L,KP,KOPT
      REAL*8     RL
      COMPLEX*16 ZVI(4)
      COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE::VVI

      INTEGER KCHECK
      KCHECK=0

C     SET KI-CONTRIBUTIONS VANISH NEAR RATIONAL SURFACES
      IF (ABS(WFUN(JS,KGRID)).LT.1.0E-13) THEN      
         IF (KOPT.EQ.0) VI  = 0.
         IF (KOPT.EQ.1) VI1 = 0.
         IF (KOPT.EQ.2) VI2 = 0.
         IF (KOPT.EQ.3) VI3 = 0.
         RETURN
      ENDIF

      ALLOCATE( VVI(4,MLMAX,NSPECIES) )
      VVI = 0.

      IF (KPARTICLE.EQ.0) THEN
C     CONTRIBUTION FROM TRAPPED PARTICLES

      DO KP=1,NSPECIES   
      DO L=1,MLMAX
         RL = RLM(L)
         IF (ABS(RL).LT.0.1.AND.ABS(PSPECIES_NTD(KP)).GT.0.) THEN
            CALL KI_PRECESSION(JS,KGRID,KPARTICLE,KP,KOPT,ZVI)
            VVI(:,L,KP) = VVI(:,L,KP) + ZVI*PSPECIES_NTD(KP)
         ELSEIF (ABS(RL).GT.0.1.AND.ABS(PSPECIES_NTB(KP)).GT.0.) THEN
            CALL KI_BOUNCE(JS,KGRID,KP,KOPT,L,ZVI)
            VVI(:,L,KP) = VVI(:,L,KP) + ZVI*PSPECIES_NTB(KP)
         ENDIF
      ENDDO
      ENDDO

      ELSEIF (KPARTICLE.EQ.1) THEN
C     CONTRIBUTION FROM PASSING PARTICLES

      DO KP=1,NSPECIES  
      IF (ABS(PSPECIES_NP(KP)).GT.0.) THEN
      DO L=1,MLMAX
         CALL KI_TRANSIT(JS,KGRID,KP,KOPT,L,ZVI)
         VVI(:,L,KP) = VVI(:,L,KP) + ZVI*PSPECIES_NP(KP)
      ENDDO
      ENDIF
      ENDDO

      ENDIF
      
      IF (KOPT.EQ.0) VI  = VVI
      IF (KOPT.EQ.1) VI1 = VVI
      IF (KOPT.EQ.2) VI2 = VVI
      IF (KOPT.EQ.3) VI3 = VVI

      DEALLOCATE(VVI)

      RETURN
      END

C=======================================================================
C ANALYTICAL INTEGRATION OVER PARTICLE ENERGY,                         
C AND NUMERICAL INTEGRATION OVER LAMNBDA                               
C FOR PARTICLE PRECESSION DRIFT RESONANCE                              
C FOR RLM(L)=0                                                         
C USING SPLINE TO OBTAIN DRIFT(LAM) AT A DENSE MESH 
C KOPT = 0: ZERO ORBIT WIDTH CONTRIBUTION
C        2: FOW CORRECTION 2 
C        3: FOW CORRECTION 3                                                             
C YQL, 10-2007                                                         
C=======================================================================
      SUBROUTINE KI0(JS,KGRID,KOPT)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KP,KOPT

      REAL*8     RTMP,REPS
      COMPLEX*16 ZVI(4),OMEGAN,VVI0(4,NSPECIES0),CTMP
      PARAMETER  (REPS=1.0E-14)

      INTEGER    NN,N1,N2,N3,J,K
      PARAMETER  (NN=200)
      REAL*8     LAMN(NN),LAMNG(2*(NN-1)),DRIFTNG(2*(NN-1)),
     &           LAMNH(2*(NN-1)),LAMX(NN)
      REAL*8     LAM0,LAM1,LAM2,LAM3,LAM4,LAM5,LAMD1,LAMD2,LAMD3,LAMA
      PARAMETER  (LAMA=3.0)

      REAL*8 :: HEPK
      REAL*8, DIMENSION(:),ALLOCATABLE:: FNUM,EFAC
      INTEGER KCHECK, NUMSIG
      KCHECK=0
      NUMSIG=3

C     FIRST ORDER FOW CONTRIBUTION VANISHES FOR PRECESSIONAL DRIFTS
      IF (KOPT.EQ.2) THEN
         VI02 = 0.
         RETURN
      ENDIF
      IF (KOPT.EQ.3) THEN
         VI03 = 0.
         RETURN
      ENDIF

C     SET KI-CONTRIBUTIONS VANISH NEAR RATIONAL SURFACES
      IF (ABS(WFUN(JS,KGRID)).LT.1.0E-13) THEN      
         IF (KOPT.EQ.0) VI0  = 0.
         RETURN
      ENDIF

C     DEFINE DENSE LAMBDA MESH, SYMMETRIC AROUND LAM0
C     TRAPPING BOUNDARY AND TO LAM0
      LAM0 = SLAMD0
      LAM1 = LAMK0(JS,1,KGRID)
      LAM2 = LAMK0(JS,NLAMK0(JS,KGRID),KGRID)
      LAMD1= MIN((LAM0-LAM1)*0.5,LAM2-LAM0)*0.9
      LAMD2= LAM0-LAM1-2.*LAMD1
      LAMD3= LAM2-LAM0-LAMD1
      LAM3 = LAM1+LAMD1     
      LAM4 = LAM3+LAMD2
      LAM5 = LAM0+LAMD1

      RTMP = 3.*LAMA*LAMD1+LAMD2+LAMD3
      N1 = INT(LAMA*LAMD1/RTMP*(NN-1))
      N2 = INT(LAMD2/RTMP*(NN-1))
      N3 = NN-1-3*N1-N2

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         WRITE(*,*) 'N1,N2,N3=',N1,N2,N3
         WRITE(*,*) 'LAMM=',LAM1,LAM3,LAM4,LAM0,LAM5,LAM2
         WRITE(*,*) 'LAMD=',LAMD1,LAMD2,LAMD3
      ENDIF

      RTMP = 1.0/DFLOAT(N1)
      LAMX(1) = 0.0
      DO J=2,N1+1
         LAMX(J) = LAMX(J-1) + RTMP                  
      ENDDO
      DO J=2,N1
         LAMX(J) = LAMX(J)**LAMA
      ENDDO

      LAMN(1) = LAM1
      LAMN(N1+1) = LAM3
      LAMN(N1+N2+1) = LAM4
      LAMN(2*N1+N2+1) = LAM0
      LAMN(3*N1+N2+1) = LAM5
      LAMN(NN) = LAM2

      DO J=2,N1
         LAMN(J) = LAM1 + LAMD1*LAMX(J)
         LAMN(N1+N2+J) = LAM0 - LAMD1*LAMX(N1+2-J)
         LAMN(2*N1+N2+J) = LAM0 + LAMD1*LAMX(J)
      ENDDO
      RTMP = LAMD2/DFLOAT(N2)
      DO J=N1+2,N1+N2
         LAMN(J) = LAMN(J-1) + RTMP
      ENDDO
      RTMP = LAMD3/DFLOAT(N3)
      DO J=3*N1+N2+2,NN-1
         LAMN(J) = LAMN(J-1) + RTMP
      ENDDO

C     DEFINE GAUSSIAN POINTS FOR THE DENSE MESH
      DO J=1,NN-1
      RTMP = (LAMN(J+1)-LAMN(J))*0.5/B0K/SQRT(PI)
      DO K=0,1
         LAMNG(J+J+K-1) = ((1+WK)*LAMN(J+K)+(1-WK)*LAMN(J-K+1))*0.5
         LAMNH(J+J+K-1) = RTMP 
      ENDDO
      ENDDO

C     SPLINE DRIFT(LAM)
      CALL SPLINE1D(DRIFTNG,LAMNG,2*(NN-1),ZOMEGADT(JS,:,KGRID),LAMM,
     &              2*NLAMK0(JS,KGRID),LAMTMP)

C     INTEGRATION OVER NEW LAMN
      VVI0 = 0.
      DO J=1,2*(NN-1)
         DRIFT = DRIFTNG(J)
         LAM   = LAMNG(J)
         DO KP=1,NSPECIES
            IF (ABS(PSPECIES_NTD(KP)).GT.0.) THEN
            CALL KDISTRIBF_TYPE(JS,KGRID,KP,LAM)
            CALL KI_PRECESSION(JS,KGRID,0,KP,KOPT,ZVI)
            VVI0(:,KP) = VVI0(:,KP) + ZVI*PSPECIES_NTD(KP)*LAMNH(J)
            ENDIF
         ENDDO
      ENDDO

C     ADD AN IMAGINARY PART DUE TO INTEGRATION ACROSS DRIFT=0
C     NOTE: THIS IMAGINARY PART COMES FROM THE INTEGRATION OVER THE
C           PITCH ANGLE, WHEN OMEGAD CROSSING ZERO. THE ION AND ELECTRON
C           CONTRIBUTIONS DO NOT CANCEL EACH OTHER, BUT ADD UP INSTEAD.
C           THIS IS EQUIVALENT TO THE LANDAU PRESCIPTION. THIS TERM SHOULD 
C           BE PRESENT INDEPENDENT OF ATAU AND OMEGAE.
      IF (NUMSIG.EQ.1.OR.NUMSIG.EQ.3) THEN
C     NUMERICAL TREATMENT OF OMEGAD SINGULARITY FOR ANY SPECIES
      DO KP=1,NSPECIES
         OMEGAN  = RNTOR*OMEGAE0(JS,KGRID)-OMEGA-CI*NUEFF(JS,KGRID,KP)
         IF (ABS(OMEGAN).GE.REPS) CYCLE
         ALLOCATE(FNUM(NEPK2),EFAC(NEPK2))
         EFAC = 1.
         IF (ISPECIES_EK(KP).EQ.0) THEN 
            ZEPK = ZEPKN
            EFAC = 0.5*PI/(COS(0.5*PI*ZEPKO))**2
         ENDIF
         IF (ISPECIES_EK(KP).EQ.1) THEN 
            ZEPK = ZEPKO
            EFAC = 1.
         ENDIF

         FNUM = ZRESU(:,2,1,KP)*ZEPK**1.5
         IF (ISPECIES_EK(KP).EQ.0) FNUM = FNUM*EFAC
         RTMP = 0
         DO J=1,NEPK-1
            HEPK = EPK(J+1)-EPK(J)
            RTMP = RTMP + 0.5*HEPK*(FNUM(2*J-1)+FNUM(2*J))
         ENDDO
         CTMP=CI*RTMP*SQRT(PI)*2/(B0K*DPRM)*ESPECIES_Z(KP)/ESPECIES_Z(1)
         VVI0(:,KP)=VVI0(:,KP)+CTMP
         IF (NUMSIG.EQ.3) THEN
            WRITE(*,*)'KI0 NUMSIG=1 RESIDU=',CTMP,
     $                'JS=',JS,'KGRID=',KGRID,'KP=',KP
         ENDIF
         DEALLOCATE(FNUM,EFAC)
      ENDDO
      ENDIF
      IF (NUMSIG.EQ.0.OR.NUMSIG.EQ.3) THEN
      DO KP=1,NSPECIES
         OMEGAN  = RNTOR*OMEGAE0(JS,KGRID)-OMEGA-CI*NUEFF(JS,KGRID,KP)
         IF (ABS(OMEGAN).LT.REPS) THEN
            RTMP=1.5*PI/DPRM*ESPECIES_DEN(JS,KGRID,KP)*
     &           OMEGACI0/B0K**2*PSPECIES_NTD(KP)
         IF (NUMSIG.EQ.3) THEN
            CTMP=-CI*RTMP*(OMEGASN(JS,KGRID,KP)
     &              +OMEGAST(JS,KGRID,KP)+CI*NUEFF(JS,KGRID,KP)/RNTOR)
            WRITE(*,*)'KI0 NUMSIG=0 RESIDU=',CTMP,
     $                'JS=',JS,'KGRID=',KGRID,'KP=',KP
         ELSE
            VVI0(:,KP)=VVI0(:,KP)-CI*RTMP*(OMEGASN(JS,KGRID,KP)
     &              +OMEGAST(JS,KGRID,KP)+CI*NUEFF(JS,KGRID,KP)/RNTOR)
         ENDIF
         ENDIF
      ENDDO
      ENDIF

      IF (KOPT.EQ.0) VI0  = VVI0

      IF (KCHECK.EQ.1.AND.KGRID.EQ.1.AND.JS.EQ.JS0) THEN
         WRITE(*,120) LAM0,DPRM,VI0
      ENDIF
 120  FORMAT('CHECK KI0: LAM0,DPRM,VI0=',4E14.5)

      RETURN
      END

C=======================================================================
C COMPUTE SF0 FACTORS, FOR ANALYTIC SUBTRACTION OF SINGULARITY
C OF PARTICLE PITCH ANGLE INTERGATION FOR NON-ADIABATIC CONTRIBUTIONS
C SF0 = I-FACTOR AT SPECIAL PITCH ANGLE
C KOPT = 0: SF0(:,:,0,:): VANISHING ORBIT WIDTH CONTRIBUTION
C        1: SF0(:,:,1,:): FOW CORRECTION 1 
C        2: SF0(:,:,2,:): FOW CORRECTION 2 
C        3: SF0(:,:,3,:): FOW CORRECTION 3 
C KDPHI = 1: SF0(:,:,:,2): FOR PRESSURE EQ AND DPHI VARIABLE
C         0: SF0(:,:,:,1): FOR PRESSURE EQ AND ALL OTHER VARIABLES
C         3: SF0(:,:,:,4): FOR DPHI EQ AND DPHI VARIABLE
C         2: SF0(:,:,:,3): FOR DPHI EQ AND ALL OTHER VARIABLES
C YQL, 08-2013 
C=======================================================================
      SUBROUTINE KSF0(JS,KGRID,KPARTICLE,KOPT)

      USE GLOBALM
      USE KINETICM
      USE ANISOTROPICM
      USE MPIENV
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KPARTICLE,L,KP,KOPT,N,K
      REAL*8     RL,H1,H2,H3,H4,H5,H6
      COMPLEX*16 ZVI

      INTEGER KCHECK
      KCHECK = 1
      ZVI    = 0.0

      IF (KPARTICLE.EQ.0) THEN
C     CONTRIBUTION FROM TRAPPED PARTICLES

      DO KP=1,NSPECIES   
      DO L=1,MLMAX
         IF (SLAM0(L,KP).GT.0.) THEN
            LAM = SLAM0(L,KP)
            CALL KDISTRIBF_TYPE(JS,KGRID,KP,LAM)
            OMEGAB =-(RNTOR*OMEGAE0(JS,KGRID)-DREAL(OMEGA))/
     &               RLM(L)/SQRT(2.*EPSALPHA(JS,KGRID,KP)/
     &               ESPECIES_M(KP)*ESPECIES_M(1))
            DRIFT  =-(RNTOR*OMEGAE0(JS,KGRID)-DREAL(OMEGA))*
     &               OMEGACI0*ESPECIES_Z(KP)/(RNTOR*B0K*
     &               ESPECIES_Z(1)*EPSALPHA(JS,KGRID,KP))
            CALL KIA_TRAP0(JS,KGRID,KP,KOPT,L,LAM,0,ZVI)
            SF0(L,KP,KOPT,1) = ZVI
            IF (INCDPHI.GT.0) THEN
               CALL KIA_TRAP0(JS,KGRID,KP,KOPT,L,LAM,1,ZVI)
               SF0(L,KP,KOPT,2) = ZVI
               CALL KIA_TRAP0(JS,KGRID,KP,KOPT,L,LAM,2,ZVI)
               SF0(L,KP,KOPT,3) = ZVI
               CALL KIA_TRAP0(JS,KGRID,KP,KOPT,L,LAM,3,ZVI)
               SF0(L,KP,KOPT,4) = ZVI
            ENDIF
         ENDIF
      ENDDO
      ENDDO

      ELSEIF (KPARTICLE.EQ.1) THEN
C     CONTRIBUTION FROM PASSING PARTICLES

      N = 2*NLAMK1(JS,KGRID)-2
      K = KGRID
      DO KP=1,NSPECIES  
      DO L=1,MLMAX
         IF (SLAM0(L,KP).GT.0.) THEN
            LAM = SLAM0(L,KP)
            CALL KDISTRIBF_TYPE(JS,KGRID,KP,LAM)
            OMEGAB = ABS(RNTOR*OMEGAE0(JS,KGRID)-DREAL(OMEGA))/
     &               ABS(RNTOR*RQK+RLM(L))/SQRT(2.*
     &               EPSALPHA(JS,KGRID,KP)/ESPECIES_M(KP)*
     &               ESPECIES_M(1))
            CALL SPLINE1D(H1,LAM,1,TPSI0    (JS,:,K),LAMM(2),N,LAMTMP)
            CALL SPLINE1D(H2,LAM,1,TPSI0DPSI(JS,:,K),LAMM(2),N,LAMTMP)
            CALL SPLINE1D(H3,LAM,1,TPSI0DLAM(JS,:,K),LAMM(2),N,LAMTMP)
            CALL SPLINE1D(H4,LAM,1,HPSI0    (JS,:,K),LAMM(2),N,LAMTMP)
            CALL SPLINE1D(H5,LAM,1,HPSI0DPSI(JS,:,K),LAMM(2),N,LAMTMP)
            CALL SPLINE1D(H6,LAM,1,HPSI0DLAM(JS,:,K),LAMM(2),N,LAMTMP)
            TPSI0L     = H1 
            TPSI0DPSIL = H2
            TPSI0DLAML = H3
            HPSI0L     = H4
            HPSI0DPSIL = H5
            HPSI0DLAML = H6
            CALL KIA_PASS0(JS,KGRID,KP,KOPT,L,LAM,0,ZVI)
            SF0(L,KP,KOPT,1) = ZVI
            IF (INCDPHI.GT.0) THEN
               CALL KIA_PASS0(JS,KGRID,KP,KOPT,L,LAM,1,ZVI)
               SF0(L,KP,KOPT,2) = ZVI
               CALL KIA_PASS0(JS,KGRID,KP,KOPT,L,LAM,2,ZVI)
               SF0(L,KP,KOPT,3) = ZVI
               CALL KIA_PASS0(JS,KGRID,KP,KOPT,L,LAM,3,ZVI)
               SF0(L,KP,KOPT,4) = ZVI
            ENDIF
         ENDIF
      ENDDO
      ENDDO

      ENDIF
      
      IF (KCHECK.EQ.1.AND.JS.EQ.JS0.AND.KGRID.EQ.1) THEN
         CALL MPI_OPEN_FILE(RANK)
         WRITE(CHMPI,*) 'KSF0: KPARTICLE=',KPARTICLE
         WRITE(CHMPI,*) 'KSF0: KP RL KOPT SF0'
         DO KP=1,NSPECIES
         DO L=1,MLMAX
            RL = RLM(L)
            IF(SLAM0(L,KP).GT.0.) WRITE(CHMPI,110) KP,RL,KOPT,
     &                            SF0(L,KP,KOPT,1),SF0(L,KP,KOPT,2),
     &                            SF0(L,KP,KOPT,3),SF0(L,KP,KOPT,4)
         ENDDO
         ENDDO
         CALL MPI_CLOSE_FILE(RANK)
      ENDIF
 110  FORMAT(I2,1X,E13.5,I2,1X,8(E13.5,1X))

      RETURN
      END

C=======================================================================
C COMPUTE SG0,SH0 FACTORS AT LAM=SLAM0, 
C FOR ANALYTIC SUBTRACTION OF SINGULARITY
C OF PARTICLE PITCH ANGLE INTERGATION FOR NON-ADIABATIC CONTRIBUTIONS
C YQL, 08-2013 
C=======================================================================
      SUBROUTINE KSGH0(JS,KGRID,KPARTICLE)

      USE GLOBALM
      USE KINETICM
      USE MPIENV
      IMPLICIT NONE

      INTEGER    JS,KGRID,KPARTICLE,L,KP

      INTEGER KCHECK
      KCHECK=1

      SVPARA0  = 0.
      SVPERP0  = 0.
      SVDPHI0  = 0.
      SVX10    = 0.
      SVX20    = 0.
      SVQ10    = 0.
      SVQ20    = 0.
      SVQ30    = 0.
      SVDP0    = 0.

      SVPARA01 = 0.
      SVPERP01 = 0.
      SVDPHI01 = 0.
      SVX101   = 0.
      SVX201   = 0.
      SVQ101   = 0.
      SVQ201   = 0.
      SVQ301   = 0.
      SVDP01   = 0.

      VPARA0   = 0.
      VPERP0   = 0.
      VDPHI0   = 0.
      VX10     = 0.
      VX20     = 0.
      VQ10     = 0.
      VQ20     = 0.
      VQ30     = 0.
      VDP0     = 0.
                                         
      VPARA01  = 0.
      VPERP01  = 0.
      VDPHI01  = 0.
      VX101    = 0.
      VX201    = 0.
      VQ101    = 0.
      VQ201    = 0.
      VQ301    = 0.
      VDP01    = 0.

      DO KP=1,NSPECIES   
      DO L=1,MLMAX
         IF (SLAM0(L,KP).GT.0.) THEN
            LAM = SLAM0(L,KP)
            IF (KPARTICLE.EQ.0) THEN
               CALL KTURN(JS,KGRID)
               CALL KCHI(0)
               CALL KEQUILK(JS,KGRID)
            ENDIF
            CALL KBTIME(JS,KGRID,KPARTICLE)
            OMEGAB = 2.*PI/RTK(NCHI2+2)
            IF (KPARTICLE.EQ.0) THEN
               OMEGAB = OMEGAB/2.
               CALL KPHI(JS,KGRID)
            ENDIF

            CALL KG(KPARTICLE)
            SVPARA0(:,L,KP) = VPARA(:,L)
            SVPERP0(:,L,KP) = VPERP(:,L)
            SVDPHI0(:,L,KP) = VDPHI(:,L)

            CALL KH(JS,KGRID,KPARTICLE)
            SVX10(:,L,KP) = VX1(:,L)
            SVX20(:,L,KP) = VX2(:,L)
            SVQ10(:,L,KP) = VQ1(:,L)
            SVQ20(:,L,KP) = VQ2(:,L)
            SVQ30(:,L,KP) = VQ3(:,L)
            SVDP0(:,L,KP) = VDP(:,L)

            IF ((IFOWP.EQ.1.AND.KPARTICLE.EQ.1).OR.
     &          (IFOWT.EQ.1.AND.KPARTICLE.EQ.0)) THEN
               CALL KG1(KPARTICLE)
               SVPARA01(:,L,KP) = VPARA1(:,L)
               SVPERP01(:,L,KP) = VPERP1(:,L)
               SVDPHI01(:,L,KP) = VDPHI1(:,L)

               CALL KH1(JS,KGRID,KPARTICLE)
               SVX101(:,L,KP) = VX11(:,L)
               SVX201(:,L,KP) = VX21(:,L)
               SVQ101(:,L,KP) = VQ11(:,L)
               SVQ201(:,L,KP) = VQ21(:,L)
               SVQ301(:,L,KP) = VQ31(:,L)
               SVDP01(:,L,KP) = VDP1(:,L)
            ENDIF
         ENDIF
      ENDDO
      ENDDO

      IF (KPARTICLE.EQ.0.AND.SLAMD0.GT.0.) THEN
         LAM = SLAMD0
         CALL KTURN(JS,KGRID)
         CALL KCHI(0)
         CALL KEQUILK(JS,KGRID)
         CALL KBTIME(JS,KGRID,0)
         OMEGAB = PI/RTK(NCHI2+2)
         CALL KPHI(JS,KGRID)

         CALL KG(0)
         VPARA0 = VPARA
         VPERP0 = VPERP
         VDPHI0 = VDPHI

         CALL KH(JS,KGRID,0)
         VX10 = VX1
         VX20 = VX2
         VQ10 = VQ1
         VQ20 = VQ2
         VQ30 = VQ3
         VDP0 = VDP

         IF (IFOWT.EQ.1) THEN
            CALL KG1(0)
            VPARA01 = VPARA1
            VPERP01 = VPERP1
            VDPHI01 = VDPHI1

            CALL KH1(JS,KGRID,0)
            VX101 = VX11
            VX201 = VX21
            VQ101 = VQ11
            VQ201 = VQ21
            VQ301 = VQ31
            VDP01 = VDP1
         ENDIF
      ENDIF

      IF (KCHECK.EQ.1.AND.JS.EQ.JS0.AND.KGRID.EQ.1) THEN
         CALL MPI_OPEN_FILE(RANK)
         WRITE(CHMPI,*) 'KSGH0: KPARTICLE=',KPARTICLE
         WRITE(CHMPI,*) 'KSGH0: KP RL SVPARA0'
         DO KP=1,NSPECIES
         DO L=1,MLMAX
            IF(SLAM0(L,KP).GT.0.) WRITE(CHMPI,110) KP,RLM(L),
     &                            SVPARA0(1,L,KP)
         ENDDO
         ENDDO
         CALL MPI_CLOSE_FILE(RANK)
      ENDIF
 110  FORMAT(I2,1X,3(E13.5,1X))

      RETURN
      END

C=======================================================================
C COMPUTE G-FACTORS FOR KINETIC PRESSURES FOR GIVEN:                   =
C   STABILITY HARMONIC RK                                              =
C   ORBIT BOUNCE HARMONIC RL                                           =
C   KPARTICLE = 0: TRAPPED; 1: PASSING                                 =
C PERFORM NUMERICAL INTEGRATION OVER CHI, USING VARIOUS INTEGRATION    =
C SCHEMES ACCORDING TO THE FOLLOWING:                                  =
C WK = 1.0: MIDPOINT                                                   =
C      0.0: TRAPIZOIDAL                                                =
C      1/SQRT(3): GAUSSIAN (THE DEFAULT)                               =
C      1/E: LOG-TYPE                                                   = 
C SINGULAR INTEGRALS FOR TRAPPED PARTICLES ARE TREATED ANALYTICALLY    =
C G-FACTORS STORED AS VPARA, VPERP AND VDPHI                           =
C YQL, 08-2007                                                         =
C=======================================================================
      SUBROUTINE KG(KPARTICLE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER    KPARTICLE,J,K,L
      REAL*8     PHASE,DIFPI
      COMPLEX*16 EPHASE,FL,FU,CTMP
      
      REAL*8 DIFFERCHI

      INTEGER KCHECK
      KCHECK=0

      VPARA = 0.0
      VPERP = 0.0
      VDPHI = 0.0

      DO J=2,NCHI2+1
         RVAK1(J) = RJBK(J)*SQRT(1-LAM/RHK(J))
         RVAK2(J) = 0.5*RJBK(J)*LAM/RHK(J)/SQRT(1-LAM/RHK(J))
         RVAK3(J) = 0.5*RJBK(J)/SQRT(1-LAM/RHK(J))
      ENDDO

      DO L=1,MLMAX
      DO K=1,MSMAX

      IF (KPARTICLE.EQ.0) THEN

C        NOTE THAT RHK(J)=LAM AT BANANA TIPS THUS VDPHI AND VPERP
C        SHARE THE SAME FL AND FU
         FL = RJBK(1)/2*SQRT(LAM)/SQRT(HPL)*EXP(-CI*RM(K,2)*CHIL)
         FU = RJBK(NCHI2+2)/2.*SQRT(LAM)/SQRT(-HPU)*EXP(CI*
     &        (-RNTOR*RPHIK(NCHI2+2)-RM(K,2)*CHIU))*COS(RLM(L)*PI)

         DO J=2,NCHI2+1
            PHASE = - RNTOR*RPHIK(J) - RM(K,2)*RCHIK(J) 
            EPHASE= EXP(CI*PHASE) * COS(RLM(L)*OMEGAB*RTK(J)) 

C        G_PARA DOES NOT REUIRE SINGULARITY TREATMENT
            VPARA(K,L) = VPARA(K,L) + RVAK1(J)*EPHASE

C        G_PERP DOES REQUIRE SINGULARITY TREATMENT
            CTMP = RVAK2(J)*EPHASE - 
     &             FL/SQRT(DIFPI(RCHIK(J)-CHIL)) - 
     &             FU/SQRT(DIFPI(CHIU-RCHIK(J)))
            VPERP(K,L) = VPERP(K,L) + CTMP

            CTMP = RVAK3(J)*EPHASE - 
     &             FL/SQRT(DIFPI(RCHIK(J)-CHIL)) - 
     &             FU/SQRT(DIFPI(CHIU-RCHIK(J)))
            VDPHI(K,L) = VDPHI(K,L) + CTMP
         ENDDO
         VPERP(K,L) = VPERP(K,L) + (FL+FU)*4*
     &                SQRT(DIFFERCHI(CHIU,CHIL))/RCHIHK
         VDPHI(K,L) = VDPHI(K,L) + (FL+FU)*4*
     &                SQRT(DIFFERCHI(CHIU,CHIL))/RCHIHK

      ELSEIF (KPARTICLE.EQ.1) THEN
         DO J=2,NCHI2+1
            PHASE = (RLM(L)+RNTOR*RQK)*OMEGAB*RTK(J)
     &            - RNTOR*RPHIK(J) - RM(K,2)*RCHIK(J) 
            EPHASE= EXP(CI*PHASE) 

C           G_PARA AND G_PERP TERMS DO NOT REQUIRE SINGULARITY TREATMENT FOR 
C           PASSING PARTICLES
            VPARA(K,L) = VPARA(K,L) + RVAK1(J)*EPHASE
            VPERP(K,L) = VPERP(K,L) + RVAK2(J)*EPHASE
            VDPHI(K,L) = VDPHI(K,L) + RVAK3(J)*EPHASE
         ENDDO
      ENDIF
      ENDDO
      ENDDO

      VPARA = VPARA*RCHIHK/4/PI
      VPERP = VPERP*RCHIHK/4/PI
      VDPHI = VDPHI*RCHIHK/4/PI

      IF (KCHECK.EQ.2.AND.
     &   ABS(RQK-Q(JS0)).LT.1.0E-13) THEN
         IF (ABS(LAMM(2)-LAM).LT.1.0E-13)
     &      WRITE(*,*) 'CHECK KG: LAM VPARA_P VPERP_P VPARA_N VPERP_N'
         DO K=1,MSMAX
           WRITE(*,120) LAM,RM(K,2),VPARA(K,1),VPERP(K,1)
         ENDDO
 120     FORMAT(6E17.8)
      ENDIF

      RETURN
      END

C=======================================================================
C G-FACTOR FUE TO FIRST ORDER FOW CORRECTION 
C=======================================================================
      SUBROUTINE KG1(KPARTICLE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE

      INTEGER    KPARTICLE,J,K,L
      REAL*8     PHASE
      COMPLEX*16 EPHASE
      
      INTEGER KCHECK
      KCHECK=0

      VPARA1 = 0.0
      VPERP1 = 0.0
      VDPHI1 = 0.0

      DO J=2,NCHI2+1
         RVAK1(J) = RJBK(J)*(RHK(J)-LAM)
         RVAK2(J) = 0.5*RJBK(J)*LAM
         RVAK3(J) = 0.5*RJBK(J)*RHK(J)
      ENDDO

      DO L=1,MLMAX
      DO K=1,MSMAX

C     TRAPPED PARTICLES DO NOT HAVE SINGULARITY INTEGRATION FOR FOW
      IF (KPARTICLE.EQ.0) THEN
         DO J=2,NCHI2+1
            PHASE = - RNTOR*RPHIK(J) - RM(K,2)*RCHIK(J) 
            EPHASE= EXP(CI*PHASE)*CI*SIN(RLM(L)*OMEGAB*RTK(J)) 
            VPARA1(K,L) = VPARA1(K,L) + RVAK1(J)*EPHASE
            VPERP1(K,L) = VPERP1(K,L) + RVAK2(J)*EPHASE
            VDPHI1(K,L) = VDPHI1(K,L) + RVAK3(J)*EPHASE
         ENDDO
      ELSEIF (KPARTICLE.EQ.1) THEN
         DO J=2,NCHI2+1
            PHASE = (RLM(L)+RNTOR*RQK)*OMEGAB*RTK(J)
     &            - RNTOR*RPHIK(J) - RM(K,2)*RCHIK(J) 
            EPHASE= EXP(CI*PHASE) 
            VPARA1(K,L) = VPARA1(K,L) + RVAK1(J)*EPHASE
            VPERP1(K,L) = VPERP1(K,L) + RVAK2(J)*EPHASE
            VDPHI1(K,L) = VDPHI1(K,L) + RVAK3(J)*EPHASE
         ENDDO
      ENDIF
      ENDDO
      ENDDO

      VPARA1 = VPARA1*RCHIHK/4/PI
      VPERP1 = VPERP1*RCHIHK/4/PI
      VDPHI1 = VDPHI1*RCHIHK/4/PI

      RETURN
      END

C=======================================================================
C COMPUTE H-FACTORS FOR KINETIC PRESSURES FOR GIVEN:                   =
C   STABILITY HARMONIC RMS                                             =
C   ORBIT BOUNCE HARMONIC RL                                           =
C   KPARTICLE = 0: TRAPPED; 1: PASSING                                 =
C PERFORM NUMERICAL INTEGRATION OVER CHI, USING VARIOUS INTEGRATION    =
C SCHEMES ACCORDING TO THE FOLLOWING:                                  =
C WK = 1.0: MIDPOINT                                                   =
C      0.0: TRAPIZOIDAL                                                =
C      1/SQRT(3): GAUSSIAN (THE DEFAULT)                               =
C      1/E: LOG-TYPE                                                   = 
C SINGULAR INTEGRALS FOR TRAPPED PARTICLES ARE TREATED ANALYTICALLY    =
C YQL, 08-2007                                                         =
C
C IMPORTANT NOTES:
C   IN THE CODE IMPLEMENTATION, BOTH SIDES OF EQUATIONS FOR P_PARA AND =
C   P_PERP ARE MULTIPLIED BY (GAMMA+I*N*OMEGAE), WHICH THEN LEADS TO   =
C   THE COUPLING OF P_PARA (_P_PERP) TO V1 & V2, INSTEAD OF X1 AND X2  =
C   THIS EXPLAINS THE ADDITIONAL TERM RX1RK BELOW, ASSOCIATED WITH     = 
C   X1(V1). THIS TRICK IMPROVES NUMERICAL CONVERGENCE OF INVERSE       =
C   ITERATION, WHEN KINETIC EQUATIONS ARE SOLVED TOGETHER WITH FLUID   =
C   EQUATIONS.
C=======================================================================
      SUBROUTINE KH(JS,KGRID,KPARTICLE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KPARTICLE,J,M,L
      REAL*8     PHASE,PHASE0,DIFPI,OMEGAE,DPSIS
      COMPLEX*16 EPHASE,CTMP,CTMPL,CTMPU,CALPHA,
     &           FLX1,FLX2,FLQ1,FLQ2,FLQ3,FLDP,
     &           FUX1,FUX2,FUQ1,FUQ2,FUQ3,FUDP
     
      REAL*8 DIFFERCHI
      
      INTEGER KCHECK

      VX1 = 0.0
      VX2 = 0.0
      VQ1 = 0.0
      VQ2 = 0.0
      VQ3 = 0.0
      VDP = 0.0

      IF (KGRID.EQ.1) THEN
C        OMEGAE = OMEGAE0(JS,1)
         OMEGAE = ROT(JS)
         DPSIS  = DPSIDS(JS)
      ELSEIF (KGRID.EQ.2) THEN
C        OMEGAE = OMEGAE0(JS,2)
         OMEGAE = ROTM(JS)
         DPSIS  = DPSIDSM(JS)
      ENDIF

      CALPHA = CI/(OMEGA-RNTOR*OMEGAE)
      IF (IPERTURB.NE.0) CALPHA = 0.0
      IF (V2XKEY.EQ.1 .OR. V2XKEY.EQ.3) CALPHA = 0.0

      IF (KPARTICLE.EQ.0) PHASE0 = 4*SQRT(DIFFERCHI(CHIU,CHIL))/RCHIHK

      DO J=2,NCHI2+1
         PHASE = 1.0/SQRT(1.0-LAM/RHK(J))
         RVAK1(J) = PHASE*((1-LAM/RHK(J))*RX1PK(J)+
     &                     (2-LAM/RHK(J))*RX1BK(J))
         RVAK6(J) = PHASE*(2-LAM/RHK(J))*RX1RK(J)
         RVAK2(J) = PHASE*(2-LAM/RHK(J))*RX2K(J)
         RVAK3(J) = PHASE*LAM*RQ1K(J)
         RVAK4(J) = PHASE*LAM*RQ2K(J)
         RVAK5(J) = PHASE*LAM*RQ3K
         RVAK7(J) = PHASE*RJBK(J)/DPSIS
      ENDDO

C     PASSING PARTICLE DOES NOT HAVE SINGULAR INTEGRATION
      IF (KPARTICLE.EQ.1) THEN
         DO L=1,MLMAX
         DO M=1,MSMAX
         DO J=2,NCHI2+1
            PHASE = -(RLM(L)+RNTOR*RQK)*OMEGAB*RTK(J) +
     &              RNTOR*RPHIK(J) + RM(M,2)*RCHIK(J)
            EPHASE= EXP(CI*PHASE) 
            VX1(M,L) = VX1(M,L) + EPHASE*
     &                           (RVAK1(J)-RVAK6(J)*CALPHA)
            VX2(M,L) = VX2(M,L) + EPHASE*RVAK2(J)
            VQ1(M,L) = VQ1(M,L) + EPHASE*RVAK3(J)
            VQ2(M,L) = VQ2(M,L) + EPHASE*RVAK4(J)
            VQ3(M,L) = VQ3(M,L) + EPHASE*RVAK5(J)
            VDP(M,L) = VDP(M,L) + EPHASE*RVAK7(J)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF (KPARTICLE.EQ.0) THEN
         DO L=1,MLMAX
         DO M=1,MSMAX
         CTMP = SQRT(LAM)*EXP(CI*RM(M,2)*CHIL)*2.0/SQRT(HPL)
         FLX1 = (RX1BK(1)-RX1RK(1)*CALPHA)*CTMP
         FLX2 = RX2K(1)*CTMP
         FLQ1 = RQ1K(1)*LAM*CTMP
         FLQ2 = RQ2K(1)*LAM*CTMP
         FLQ3 = RQ3K*LAM*CTMP
         FLDP = RJBK(1)/DPSIS*CTMP

         CTMP = SQRT(LAM)*EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &          2.0*COS(RLM(L)*PI)/SQRT(-HPU)
         FUX1 = (RX1BK(NCHI2+2)-RX1RK(NCHI2+2)*CALPHA)*CTMP
         FUX2 = RX2K(NCHI2+2)*CTMP
         FUQ1 = RQ1K(NCHI2+2)*LAM*CTMP
         FUQ2 = RQ2K(NCHI2+2)*LAM*CTMP
         FUQ3 = RQ3K*LAM*CTMP
         FUDP = RJBK(NCHI2+2)/DPSIS*CTMP

         DO J=2,NCHI2+1
            PHASE = RNTOR*RPHIK(J) + RM(M,2)*RCHIK(J)
            EPHASE= EXP(CI*PHASE)*2*COS(RLM(L)*OMEGAB*RTK(J))
            CTMPL=1.0/SQRT(DIFPI(RCHIK(J)-CHIL))
            CTMPU=1.0/SQRT(DIFPI(CHIU-RCHIK(J)))
            VX1(M,L)=VX1(M,L)+EPHASE*
     &             (RVAK1(J)-RVAK6(J)*CALPHA)-CTMPL*FLX1-CTMPU*FUX1 
            VX2(M,L)=VX2(M,L)+EPHASE*RVAK2(J)-CTMPL*FLX2-CTMPU*FUX2
            VQ1(M,L)=VQ1(M,L)+EPHASE*RVAK3(J)-CTMPL*FLQ1-CTMPU*FUQ1
            VQ2(M,L)=VQ2(M,L)+EPHASE*RVAK4(J)-CTMPL*FLQ2-CTMPU*FUQ2
            VQ3(M,L)=VQ3(M,L)+EPHASE*RVAK5(J)-CTMPL*FLQ3-CTMPU*FUQ3
            VDP(M,L)=VDP(M,L)+EPHASE*RVAK7(J)-CTMPL*FLDP-CTMPU*FUDP
         ENDDO
         VX1(M,L) = VX1(M,L) + (FLX1+FUX1)*PHASE0
         VX2(M,L) = VX2(M,L) + (FLX2+FUX2)*PHASE0
         VQ1(M,L) = VQ1(M,L) + (FLQ1+FUQ1)*PHASE0
         VQ2(M,L) = VQ2(M,L) + (FLQ2+FUQ2)*PHASE0
         VQ3(M,L) = VQ3(M,L) + (FLQ3+FUQ3)*PHASE0
         VDP(M,L) = VDP(M,L) + (FLDP+FUDP)*PHASE0
         ENDDO
         ENDDO
      ENDIF

      PHASE = RCHIHK/4*OMEGAB/PI
      VX1 = VX1*PHASE
      VX2 = VX2*PHASE
      VQ1 = VQ1*PHASE
      VQ2 = VQ2*PHASE
      VQ3 = VQ3*PHASE
      VDP = VDP*PHASE

      KCHECK=0
      IF (KCHECK.EQ.1.AND.KPARTICLE.EQ.0.AND.
     &   ABS(RQK-Q(JS0)).LT.1.0E-13) THEN
         PHASE=0.10121399E+01
         IF (ABS(PHASE-LAM).LT.1.0E-05) THEN 
            WRITE(*,*) 'CHECK KH: RCHIK RTK RPHIK VX1 VX2 VQ1 VQ2 VQ3 
     &                  VDP'
            DO J=2,NCHI2+1
            PHASE = RNTOR*RPHIK(J)+RM(1,2)*RCHIK(J)-RLM(1)*OMEGAB*RTK(J)
            EPHASE= EXP(CI*PHASE)
            CTMP=EPHASE/SQRT(1.0-LAM/RHK(J))*OMEGAB/2/PI
            FLX1 = CTMP*((1-LAM/RHK(J))*RX1PK(J)+
     &            (2-LAM/RHK(J))*RX1BK(J))
            FLX2 = CTMP*(2-LAM/RHK(J))*RX2K(J) 
            FLQ1 = CTMP*LAM*RQ1K(J)
            FLQ2 = CTMP*LAM*RQ2K(J)
            FLQ3 = CTMP*LAM*RQ3K 
            FLDP = CTMP*RJBK(J)/DPSIS
            WRITE(*,110) RCHIK(J),RTK(J),RPHIK(J),
     &                   FLX1,FLX2,FLQ1,FLQ2,FLQ3,FLDP
            ENDDO
            DO J=NCHI2+1,2,-1
            PHASE = RNTOR*RPHIK(J)+RM(1,2)*RCHIK(J)+RLM(1)*OMEGAB*RTK(J)
            EPHASE= EXP(CI*PHASE)
            CTMP=EPHASE/SQRT(1.0-LAM/RHK(J))*OMEGAB/2/PI
            FLX1 = CTMP*((1-LAM/RHK(J))*RX1PK(J)+
     &            (2-LAM/RHK(J))*RX1BK(J))
            FLX2 = CTMP*(2-LAM/RHK(J))*RX2K(J) 
            FLQ1 = CTMP*LAM*RQ1K(J)
            FLQ2 = CTMP*LAM*RQ2K(J)
            FLQ3 = CTMP*LAM*RQ3K 
            FLDP = CTMP*RJBK(J)/DPSIS
            WRITE(*,110) RCHIK(J),2*PI/OMEGAB-RTK(J),RPHIK(J),
     &                   FLX1,FLX2,FLQ1,FLQ2,FLQ3,FLDP
            ENDDO
         ENDIF
 110     FORMAT(15E17.8)
      ENDIF

      IF (KCHECK.EQ.2.AND.KPARTICLE.EQ.0.AND.
     &   ABS(RQK-Q(JS0)).LT.1.0E-13) THEN
         PHASE=1.0/1.0099
         IF (ABS(PHASE-LAM).LT.1.0E-13) THEN 
            WRITE(*,*) 'CHECK KH: RCHIK RTK RPHIK VX1 VX2 VQ1 VQ2 VQ3 
     &                  VDP'
            DO J=1,NCHI2+2
            PHASE = RNTOR*RPHIK(J) - RLM(1)*OMEGAB*RTK(J)
            EPHASE= EXP(-CI*PHASE)
            FLX1 = VX1(1,1)*EPHASE
            FLX2 = VX2(1,1)*EPHASE
            FLQ1 = VQ1(1,1)*EPHASE
            FLQ2 = VQ2(1,1)*EPHASE
            FLQ3 = VQ3(1,1)*EPHASE
            FLDP = VDP(1,1)*EPHASE
            WRITE(*,110) RCHIK(J),RTK(J),RPHIK(J),
     &                   FLX1,FLX2,FLQ1,FLQ2,FLQ3,FLDP
            ENDDO
            DO J=NCHI2+1,1,-1
            PHASE = RNTOR*RPHIK(J) + RLM(1)*OMEGAB*RTK(J)
            EPHASE= EXP(-CI*PHASE)
            FLX1 = VX1(1,1)*EPHASE
            FLX2 = VX2(1,1)*EPHASE
            FLQ1 = VQ1(1,1)*EPHASE
            FLQ2 = VQ2(1,1)*EPHASE
            FLQ3 = VQ3(1,1)*EPHASE
            FLDP = VDP(1,1)*EPHASE
            WRITE(*,110) RCHIK(J),2*PI/OMEGAB-RTK(J),RPHIK(J),
     &                   FLX1,FLX2,FLQ1,FLQ2,FLQ3,FLDP
            ENDDO
         ENDIF
      ENDIF

      IF (KCHECK.EQ.3.AND.
     &   ABS(RQK-QM(JS0)).LT.1.0E-13) THEN
         IF (ABS(LAMM(2)-LAM).LT.1.0E-13)
     &      WRITE(*,*) 'CHECK KH: LAM VX1 VX2 VQ1 VQ2 VQ3 VDP'
         WRITE(*,120) LAM,VX1(1,1),VX2(1,1),VQ1(1,1),
     &                VQ2(1,1),VQ3(1,1),VDP(1,1)
 120     FORMAT(13E17.8)
      ENDIF

      IF (KCHECK.EQ.4.AND.ABS(LAM-0.5).LT.1.0E-13.AND.
     &   ABS(RQK-Q(JS0)).LT.1.0E-13) THEN
         WRITE(*,*) 'CHECK KH: M L VX1 VX2 VQ1 VQ2 VQ3 VDP'
         DO M=1,MSMAX
         DO L=1,MLMAX
         WRITE(*,130) RM(M,2),RLM(L),ABS(VX1(M,L)),ABS(VX2(M,L)),
     &                ABS(VQ1(M,L)),ABS(VQ2(M,L)),ABS(VQ3(M,L)),
     &                ABS(VDP(M,L))
         ENDDO
         ENDDO
 130     FORMAT(2F6.1,6E17.8)
      ENDIF

      RETURN
      END

C=======================================================================
C H-FACTOR DUE TO FIRST ORDER FOW CORRECTION
C YQL, 06-2013
C=======================================================================
      SUBROUTINE KH1(JS,KGRID,KPARTICLE)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER    JS,KGRID,KPARTICLE,J,M,L
      REAL*8     PHASE,OMEGAE,DPSIS
      COMPLEX*16 EPHASE,CTMP,CTMPL,CTMPU,CALPHA
     
      INTEGER KCHECK

      VX11 = 0.0
      VX21 = 0.0
      VQ11 = 0.0
      VQ21 = 0.0
      VQ31 = 0.0
      VDP1 = 0.0

      IF (KGRID.EQ.1) THEN
C        OMEGAE = OMEGAE0(JS,1)
         OMEGAE = ROT(JS)
         DPSIS  = DPSIDS(JS)
      ELSEIF (KGRID.EQ.2) THEN
C        OMEGAE = OMEGAE0(JS,2)
         OMEGAE = ROTM(JS)
         DPSIS  = DPSIDSM(JS)
      ENDIF

      CALPHA = CI/(OMEGA-RNTOR*OMEGAE)
      IF (IPERTURB.NE.0) CALPHA = 0.0
      IF (V2XKEY.EQ.1 .OR. V2XKEY.EQ.3) CALPHA = 0.0

      DO J=2,NCHI2+1
         PHASE = RHK(J)
         RVAK1(J) = PHASE*((1-LAM/RHK(J))*RX1PK(J)+
     &                     (2-LAM/RHK(J))*RX1BK(J))
         RVAK6(J) = PHASE*(2-LAM/RHK(J))*RX1RK(J)
         RVAK2(J) = PHASE*(2-LAM/RHK(J))*RX2K(J)
         RVAK3(J) = PHASE*LAM*RQ1K(J)
         RVAK4(J) = PHASE*LAM*RQ2K(J)
         RVAK5(J) = PHASE*LAM*RQ3K
         RVAK7(J) = PHASE*RJBK(J)/DPSIS
      ENDDO

C     PASSING PARTICLE DOES NOT HAVE SINGULAR INTEGRATION
      IF (KPARTICLE.EQ.1) THEN
         DO L=1,MLMAX
         DO M=1,MSMAX
         DO J=2,NCHI2+1
            PHASE = -(RLM(L)+RNTOR*RQK)*OMEGAB*RTK(J) +
     &              RNTOR*RPHIK(J) + RM(M,2)*RCHIK(J)
            EPHASE= EXP(CI*PHASE) 
            VX11(M,L) = VX11(M,L) + EPHASE*
     &                           (RVAK1(J)-RVAK6(J)*CALPHA)
            VX21(M,L) = VX21(M,L) + EPHASE*RVAK2(J)
            VQ11(M,L) = VQ11(M,L) + EPHASE*RVAK3(J)
            VQ21(M,L) = VQ21(M,L) + EPHASE*RVAK4(J)
            VQ31(M,L) = VQ31(M,L) + EPHASE*RVAK5(J)
            VDP1(M,L) = VDP1(M,L) + EPHASE*RVAK7(J)
         ENDDO
         ENDDO
         ENDDO
C     TRAPPED PARTICLES DO NOT HAVE SINGULARITY EITHER
      ELSEIF (KPARTICLE.EQ.0) THEN
         DO L=1,MLMAX
         DO M=1,MSMAX
         DO J=2,NCHI2+1
            PHASE = RNTOR*RPHIK(J) + RM(M,2)*RCHIK(J)
            EPHASE=-EXP(CI*PHASE)*2.*CI*SIN(RLM(L)*OMEGAB*RTK(J))
            VX11(M,L)=VX11(M,L)+EPHASE*(RVAK1(J)-RVAK6(J)*CALPHA)
            VX21(M,L)=VX21(M,L)+EPHASE*RVAK2(J)
            VQ11(M,L)=VQ11(M,L)+EPHASE*RVAK3(J)
            VQ21(M,L)=VQ21(M,L)+EPHASE*RVAK4(J)
            VQ31(M,L)=VQ31(M,L)+EPHASE*RVAK5(J)
            VDP1(M,L)=VDP1(M,L)+EPHASE*RVAK7(J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      PHASE = RCHIHK/4*OMEGAB/PI
      VX11 = VX11*PHASE
      VX21 = VX21*PHASE
      VQ11 = VQ11*PHASE
      VQ21 = VQ21*PHASE
      VQ31 = VQ31*PHASE
      VDP1 = VDP1*PHASE

      RETURN
      END

C=======================================================================
C COMPUTE PERTURBED ENERGY USING THE COMPUTED EIGENFUNCTION AND        =
C THE SYSTEM MATRICES                                                  =
C INCLUDING:                                                           = 
C  1) FLUID POTENTIAL ENERGY                                           = 
C  2) FLIUD KINETIC ENERGY DUE TO INERTIAL AND PLASMA ROTATION         =
C  3) KINETIC ENERGY DUE TO KINETIC PRESSURE PERTURBATION              =
C SHOULD BE CALLED FROM MARS-F AFTER CALPAM(...)                       =
C KENORM = 1: NORMALIZE ENERGY BY TOTAL PLASMA INERTIA ENERGY          =
C          2: NORMALIZE ENERGY BY INERTIA ENERGY ASSOCIATED WITH X1    =
C KEFORM = 1: RAW FORM FOR ALL DW* TERMS
C          2: QUADRATIC FORM FOR DWP,DWJ,DWPPAR,DWPPER
C YQL, 09-2013                                                         =
C=======================================================================
      SUBROUTINE ENERGYMAT(X,Y,
     &                  ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE GIJLM
      USE CONVOLCOFM
      IMPLICIT NONE
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'

      INTEGER    MROW,MSA,MSB,MSMI,MSPL,NSA,NSB
      PARAMETER  (NSA=2,NSB=1)
      INTEGER    LXCOL,LYCOL,LXROW,LYROW,I
      REAL*8     PI2
      COMPLEX*16 X(NXCOMP*MSMAX,*),Y(NYCOMP*MSMAX,*)
      COMPLEX*16 DWP,DWPPAR,DWPPER,DWJ,DWQ,DWV,DWRHO,DWXI,DWX2,DWALL,
     &           DWVAC
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::
     &           DWPX,DWPPARX,DWPPERX,DWJX,DWQX,DWVX,DWXIX,DWX2X,DWRHOX,
     &           DWPY,DWPPARY,DWPPERY,DWJY,DWQY,DWVY,DWXIY,DWX2Y,DWRHOY

      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::R,RY
      INTEGER       IEXV2,IEXE,IEXY
      PARAMETER     (IEXV2 = -1, IEXY=0)
      REAL*8        ZEM,ZEP,ZV2M,ZV2P,z3m,z3p

      COMPLEX*16    AL0SAVE,CTMP
      CHARACTER*80  LINE

      INCLUDE 'integc.inc'

C     FILL IN SYSTEM MATRICES WITH THE CONVERGED EIGENVALUE
      AL0SAVE = AL0
      AL0 = ALNORM
      KJPKEY = 0
      KPBKEY = 1
      CALL LINEAR(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
      AL0 = AL0SAVE
     
      ALLOCATE( DWPX(NTP1), DWPPARX(NTP1), DWPPERX(NTP1), DWJX(NTP1),
     &          DWQX(NTP1), DWVX(NTP1),    DWXIX(NTP1),   DWX2X(NTP1),
     &          DWRHOX(NTP1),
     &          DWPY(NTOT), DWPPARY(NTOT), DWPPERY(NTOT), DWJY(NTOT),
     &          DWQY(NTOT), DWVY(NTOT),    DWXIY(NTOT), DWX2Y(NTOT),
     &          DWRHOY(NTOT)    
     &        )

      DWPX    = 0.
      DWPPARX = 0.
      DWPPERX = 0.
      DWJX    = 0.
      DWQX    = 0.
      DWVX    = 0.
      DWXIX   = 0.
      DWX2X   = 0.
      DWRHOX  = 0.
      DWPY    = 0.
      DWPPARY = 0.
      DWPPERY = 0.
      DWJY    = 0.
      DWQY    = 0.
      DWVY    = 0.
      DWXIY   = 0.
      DWX2Y   = 0.
      DWRHOY  = 0.
      DWVAC   = 0.

      DO MROW=1,MSMAX
         LXROW = (MROW-1)*NXCOMP
         LYROW = (MROW-1)*NYCOMP
         DO MSA=1,MSMAX
            LXCOL = (MSA -1)*NXCOMP
            LYCOL = (MSA -1)*NYCOMP
            DO I=2,NRP1
               IF (KEFORM.EQ.1) THEN
               DWPX(I)=DWPX(I) + CONJG(X1U(I,MROW))*
     &           HSUBM(LXROW+KXV1,LYCOL+KYPR,I)*Y(LYCOL+KYPR,I-1)
               IF (KYPE.GT.0) 
     &         DWPX(I)=DWPX(I) + CONJG(X1U(I,MROW))*
     &           HSUBM(LXROW+KXV1,LYCOL+KYPE,I)*Y(LYCOL+KYPE,I-1)
               IF (KYPPERP.GT.0)
     &         DWPPERX(I)=DWPPERX(I) + CONJG(X1U(I,MROW))*
     &           HSUBM(LXROW+KXV1,LYCOL+KYPPERP,I)*Y(LYCOL+KYPPERP,I-1)
               DWJX(I)=DWJX(I) + CONJG(X1U(I,MROW))*(
     &           HSUBM(LXROW+KXV1,LYCOL+KYJ1,I)*Y(LYCOL+KYJ1,I-1) + 
     &           ASUBM(LXROW+KXV1,LXCOL+KXJ2U,I)*X(LXCOL+KXJ2U,I-1) +
     &           ASUBM(LXROW+KXV1,LXCOL+KXJ3,I)*X(LXCOL+KXJ3,I-1) )
               ENDIF
               IF (KEFORM.EQ.2.AND.INCKIN.GT.0) THEN
               DWPPERX(I)=DWPPERX(I) - CONJG(X1U(I,MROW))*
     &           HSUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*Y(LYCOL+KYPPERP,I-1)
               ENDIF
               IF (KYPPARA.GT.0)
     &         DWPPARX(I)=DWPPARX(I) + CONJG(X1U(I,MROW))*
     &           HSUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*Y(LYCOL+KYPPARA,I-1)
               DWQX(I)=DWQX(I) + CONJG(X1U(I,MROW))*(
     &           ASUBM(LXROW+KXV1,LXCOL+KXB1,I)*X(LXCOL+KXB1,I-1) + 
     &           HSUBM(LXROW+KXV1,LYCOL+KYB2,I)*Y(LYCOL+KYB2,I-1) +
     &           HSUBM(LXROW+KXV1,LYCOL+KYB3,I)*Y(LYCOL+KYB3,I-1) )
               DWVX(I)=DWVX(I) + CONJG(X1U(I,MROW))*(
     &           ASUBM(LXROW+KXV1,LXCOL+KXV1,I)*X(LXCOL+KXV1,I-1) + 
     &           HSUBM(LXROW+KXV1,LYCOL+KYV2,I)*Y(LYCOL+KYV2,I-1) )
               IF (KYV3.GT.0)
     &         DWVX(I)=DWVX(I) + CONJG(X1U(I,MROW))*
     &           HSUBM(LXROW+KXV1,LYCOL+KYV3,I)*Y(LYCOL+KYV3,I-1) 
               IF (KYRHOP.GT.0)
     &         DWRHOX(I)=DWRHOX(I) + CONJG(X1U(I,MROW))*(
     &           HSUBM(LXROW+KXV1,LYCOL+KYRHOP,I)*Y(LYCOL+KYRHOP,I-1) )
               IF (KXX1.GT.0)
     &         DWXIX(I)=DWXIX(I) + CONJG(X1U(I,MROW))*
     &           ASUBM(LXROW+KXV1,LXCOL+KXX1,I)*X(LXCOL+KXX1,I-1)
               IF (KYX2.GT.0)
     &         DWXIX(I)=DWXIX(I) + CONJG(X1U(I,MROW))*
     &           HSUBM(LXROW+KXV1,LYCOL+KYX2,I)*Y(LYCOL+KYX2,I-1) 
            ENDDO

            DO I=1,NR
               DWQX(I)=DWQX(I) + CONJG(X1U(I,MROW))*
     &           CSUBM(LXROW+KXV1,LXCOL+KXB1,I)*X(LXCOL+KXB1,I+1)
               DWVX(I)=DWVX(I) + CONJG(X1U(I,MROW))*
     &           CSUBM(LXROW+KXV1,LXCOL+KXV1,I)*X(LXCOL+KXV1,I+1)
               IF (KXX1.GT.0)
     &         DWXIX(I)=DWXIX(I) + CONJG(X1U(I,MROW))*
     &           CSUBM(LXROW+KXV1,LXCOL+KXX1,I)*X(LXCOL+KXX1,I+1)

               ZEM = CSM(I)**IEXY
               IF (KEFORM.EQ.1) THEN
               DWJX(I)=DWJX(I) + CONJG(X1U(I,MROW))*(
     &           CSUBM(LXROW+KXV1,LXCOL+KXJ2U,I)*X(LXCOL+KXJ2U,I+1) +
     &           CSUBM(LXROW+KXV1,LXCOL+KXJ3,I)*X(LXCOL+KXJ3,I+1) )
               DWPY(I)=DWPY(I) + ZEM*CONJG(X2U(I,MROW))*
     &           DSUBM(LYROW+KYV2,LYCOL+KYPR,I)*Y(LYCOL+KYPR,I)
               IF (KYPE.GT.0)
     &         DWPY(I)=DWPY(I) + ZEM*CONJG(X2U(I,MROW))*
     &           DSUBM(LYROW+KYV2,LYCOL+KYPE,I)*Y(LYCOL+KYPE,I)
               IF (KYPPERP.GT.0)
     &         DWPPERY(I)=DWPPERY(I) + ZEM*CONJG(X2U(I,MROW))*
     &           DSUBM(LYROW+KYV2,LYCOL+KYPPERP,I)*Y(LYCOL+KYPPERP,I)
               DWJY(I)=DWJY(I) + ZEM*CONJG(X2U(I,MROW))*(
     &           DSUBM(LYROW+KYV2,LYCOL+KYJ1,I)*Y(LYCOL+KYJ1,I) +
     &           FSUBM(LYROW+KYV2,LXCOL+KXJ2U,I)*X(LXCOL+KXJ2U,I) +
     &           FSUBM(LYROW+KYV2,LXCOL+KXJ3,I)*X(LXCOL+KXJ3,I) +
     &           GSUBM(LYROW+KYV2,LXCOL+KXJ2U,I)*X(LXCOL+KXJ2U,I+1) +
     &           GSUBM(LYROW+KYV2,LXCOL+KXJ3,I)*X(LXCOL+KXJ3,I+1) )
               ENDIF
               IF (KEFORM.EQ.2) THEN
               DWPY(I)=DWPY(I) + ZEM*Y(LYROW+KYPR,I)*CONJG(
     &          FSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I,MSA)+
     &          GSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I+1,MSA)+
     &          DSUBM(LYROW+KYPR,LYCOL+KYV2,I)*X2U(I,MSA))
               IF (KYPE.GT.0) 
     &         DWPY(I)=DWPY(I) + ZEM*Y(LYROW+KYPE,I)*CONJG(
     &          FSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I,MSA)+
     &          GSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I+1,MSA)+
     &          DSUBM(LYROW+KYPR,LYCOL+KYV2,I)*x2U(I,MSA)) 
               IF (INCKIN.GT.0) THEN
               DWPPERY(I)=DWPPERY(I) - ZEM*CONJG(X2U(I,MROW))*
     &           DSUBM(LYROW+KYV2,LYCOL+KYPPARA,I)*Y(LYCOL+KYPPERP,I)
               DWPPERY(I)=DWPPERY(I) + ZEM*Y(LYROW+KYPPERP,I)*CONJG(
     &          FSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I,MSA)+
     &          GSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I+1,MSA)+
     &          DSUBM(LYROW+KYPR,LYCOL+KYV2,I)*X2U(I,MSA)) 
               ENDIF
               ENDIF
               IF (KYPPARA.GT.0) 
     &         DWPPARY(I)=DWPPARY(I) + ZEM*CONJG(X2U(I,MROW))*
     &           DSUBM(LYROW+KYV2,LYCOL+KYPPARA,I)*Y(LYCOL+KYPPARA,I)
               DWQY(I)=DWQY(I) + ZEM*CONJG(X2U(I,MROW))*(
     &           FSUBM(LYROW+KYV2,LXCOL+KXB1,I)*X(LXCOL+KXB1,I) +
     &           GSUBM(LYROW+KYV2,LXCOL+KXB1,I)*X(LXCOL+KXB1,I+1) +
     &           DSUBM(LYROW+KYV2,LYCOL+KYB2,I)*Y(LYCOL+KYB2,I) +
     &           DSUBM(LYROW+KYV2,LYCOL+KYB3,I)*Y(LYCOL+KYB3,I) )
               DWVY(I)=DWVY(I) + ZEM*CONJG(X2U(I,MROW))*(
     &           FSUBM(LYROW+KYV2,LXCOL+KXV1,I)*X(LXCOL+KXV1,I) +
     &           GSUBM(LYROW+KYV2,LXCOL+KXV1,I)*X(LXCOL+KXV1,I+1) +
     &           DSUBM(LYROW+KYV2,LYCOL+KYV2,I)*Y(LYCOL+KYV2,I) )
               IF (KYV3.GT.0)
     &         DWVY(I)=DWVY(I) + ZEM*CONJG(X2U(I,MROW))*
     &           DSUBM(LYROW+KYV2,LYCOL+KYV3,I)*Y(LYCOL+KYV3,I) 
               IF (KYRHOP.GT.0)
     &         DWRHOY(I)=DWRHOY(I) + ZEM*CONJG(X2U(I,MROW))*(
     &           DSUBM(LYROW+KYV2,LYCOL+KYRHOP,I)*Y(LYCOL+KYRHOP,I) )
               IF (KXX1.GT.0) 
     &         DWXIY(I)=DWXIY(I) + ZEM*CONJG(X2U(I,MROW))*(
     &           FSUBM(LYROW+KYV2,LXCOL+KXX1,I)*X(LXCOL+KXX1,I) +
     &           GSUBM(LYROW+KYV2,LXCOL+KXX1,I)*X(LXCOL+KXX1,I+1) )
               IF (KYX2.GT.0) 
     &         DWXIY(I)=DWXIY(I) + ZEM*CONJG(X2U(I,MROW))*
     &           DSUBM(LYROW+KYV2,LYCOL+KYX2,I)*Y(LYCOL+KYX2,I) 

               IF (KEFORM.EQ.1) THEN
               DWPX(I)=DWPX(I) + CONJG(X1U(I,MROW))*
     &           ESUBM(LXROW+KXV1,LYCOL+KYPR,I)*Y(LYCOL+KYPR,I)
               IF (KYPE.GT.0)
     &         DWPX(I)=DWPX(I) + CONJG(X1U(I,MROW))*
     &           ESUBM(LXROW+KXV1,LYCOL+KYPE,I)*Y(LYCOL+KYPE,I)
               IF (KYPPERP.GT.0)
     &         DWPPERX(I)=DWPPERX(I) + CONJG(X1U(I,MROW))*
     &           ESUBM(LXROW+KXV1,LYCOL+KYPPERP,I)*Y(LYCOL+KYPPERP,I)
               DWJX(I)=DWJX(I) + CONJG(X1U(I,MROW))*(
     &           ESUBM(LXROW+KXV1,LYCOL+KYJ1,I)*Y(LYCOL+KYJ1,I) )
               ENDIF
               IF (KEFORM.EQ.2.AND.INCKIN.GT.0) THEN
               DWPPERX(I)=DWPPERX(I) - CONJG(X1U(I,MROW))*
     &           ESUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*Y(LYCOL+KYPPERP,I)
               ENDIF
               IF (KYPPARA.GT.0) 
     &         DWPPARX(I)=DWPPARX(I) + CONJG(X1U(I,MROW))*
     &           ESUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*Y(LYCOL+KYPPARA,I)
               DWQX(I)=DWQX(I) + CONJG(X1U(I,MROW))*(
     &           ESUBM(LXROW+KXV1,LYCOL+KYB2,I)*Y(LYCOL+KYB2,I) +
     &           ESUBM(LXROW+KXV1,LYCOL+KYB3,I)*Y(LYCOL+KYB3,I) )
               DWVX(I)=DWVX(I) + CONJG(X1U(I,MROW))*
     &           ESUBM(LXROW+KXV1,LYCOL+KYV2,I)*Y(LYCOL+KYV2,I) 
               IF (KYV3.GT.0)
     &         DWVX(I)=DWVX(I) + CONJG(X1U(I,MROW))*
     &           ESUBM(LXROW+KXV1,LYCOL+KYV3,I)*Y(LYCOL+KYV3,I) 
               IF (KYRHOP.GT.0) 
     &         DWRHOX(I)=DWRHOX(I) + CONJG(X1U(I,MROW))*(
     &           ESUBM(LXROW+KXV1,LYCOL+KYRHOP,I)*Y(LYCOL+KYRHOP,I) )
               IF (KYX2.GT.0)
     &         DWXIX(I)=DWXIX(I) + CONJG(X1U(I,MROW))*(
     &           ESUBM(LXROW+KXV1,LYCOL+KYX2,I)*Y(LYCOL+KYX2,I) )
            ENDDO

            DO I=1,NRP1
               IF (KEFORM.EQ.1) THEN
               DWJX(I)=DWJX(I) + CONJG(X1U(I,MROW))*(
     &           BSUBM(LXROW+KXV1,LXCOL+KXJ2U,I)*X(LXCOL+KXJ2U,I) +
     &           BSUBM(LXROW+KXV1,LXCOL+KXJ3,I)*X(LXCOL+KXJ3,I) )
               ENDIF
               DWQX(I)=DWQX(I) + CONJG(X1U(I,MROW))*(
     &           BSUBM(LXROW+KXV1,LXCOL+KXB1,I)*X(LXCOL+KXB1,I) )
               DWVX(I)=DWVX(I) + CONJG(X1U(I,MROW))*(
     &           BSUBM(LXROW+KXV1,LXCOL+KXV1,I)*X(LXCOL+KXV1,I) )
               IF (KXX1.GT.0)
     &         DWXIX(I)=DWXIX(I) + CONJG(X1U(I,MROW))*(
     &           BSUBM(LXROW+KXV1,LXCOL+KXX1,I)*X(LXCOL+KXX1,I) )
            ENDDO
         ENDDO
      ENDDO
      
      IF (KXX1.GT.0.AND.KYX2.GT.0)
     &CALL CALCDWKCOMP(ASUBM,BSUBM,CSUBM,DSUBM,
     &                 ESUBM,FSUBM,GSUBM,HSUBM)

C     COMPUTE X2=(XI,XI^*) USING THE EQUATION FOR V IN MUBMAT
      ALLOCATE( R(MXMAX,NTP1), RY(MYMAX,NTOT) )
      R = 0.
      RY = 0.

C     ANNIHILATE X2 COMPONENT BEFORE COMPUTING INERTIA ENERGY
      IF (KENORM.EQ.2) THEN
        DO MSA=1,MSMAX
        DO I=1,NR
           X2U(I,MSA) = 0.0
        ENDDO
        ENDDO
      ENDIF

      DO 90 MSA=1,MSMAX
      DO 90 MSB=1,MSMAX
      MSPL =  MPLUS(MSA,NSA,MSB,NSB)
      MSMI = MMINUS(MSA,NSA,MSB,NSB)
      IF (MSPL.LT.1) GOTO 60

      DO 50 I=2,NR
      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      ZV2P = (CS(I)/CSM(I  ))**IEXV2
      INCLUDE 'tent.inc'

      R (KXV1+(MSPL-1)*NXCOMP,I) = R (KXV1+(MSPL-1)*NXCOMP,I)
     &        + FF(RGV1G11(i,msb),RGV1G11M(i-1,msb),RGV1G11M(i,msb))
     &                                *X1U(I,MSA)
     &        + FFM(RGV1G11M(i-1,msb))*X1U(I-1,MSA)
     &        + FFP(RGV1G11M(i,msb))  *X1U(I+1,MSA)
     &+FGM(RGV1G12(i,msb)*zv2m,RGV1G12M(i-1,msb))*X2U(I-1,MSA)
     &+FGP(RGV1G12(i,msb)*zv2p,RGV1G12M(i,msb))*X2U(I,MSA)
 50   CONTINUE

      I = NRP1
      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      INCLUDE 'tent.inc'

      R (KXV1+(MSPL-1)*NXCOMP,I) = R (KXV1+(MSPL-1)*NXCOMP,I)
     &        + FF(RGV1G11(i,msb),RGV1G11M(i-1,msb),RGV1G11M(i,msb))
     &                                *X1U(I,MSA)
     &        + FFM(RGV1G11M(i-1,msb))*X1U(I-1,MSA)
     &+FGM(ZV2M*RGV1G12(I,MSB),RGV1G12M(I-1,MSB))*X2U(I-1,MSA)

      IEXE = -1 + IEXY

      DO 55 I=1,NR

      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      INCLUDE 'tophat.inc'

      RY(KYV2+(MSPL-1)*NYCOMP,I) = RY(KYV2+(MSPL-1)*NYCOMP,I)
     & + GG(RGV2G22M(i,msb),RGV2G22(i,msb)*zem*zv2m,
     $        RGV2G22(i+1,msb)*zep*zv2p)*X2U(I,MSA)
     & + GF(RGV1G12(i,msb)*zem,RGV1G12M(i,msb))*X1U(I,MSA)
     & + GF(RGV1G12(i+1,msb)*zep,RGV1G12M(i,msb))*X1U(I+1,MSA)
 55   CONTINUE
 60   CONTINUE
      IF (MSB.LT.2) GOTO 80
      IF (MSMI.LT.1) GOTO 80

      DO 70 I=2,NR

      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      ZV2P = (CS(I)/CSM(I  ))**IEXV2
      INCLUDE 'tent.inc'

      R (KXV1+(MSMI-1)*NXCOMP,I) = R (KXV1+(MSMI-1)*NXCOMP,I)
     &   +CONJG(FF(RGV1G11(i,msb),RGV1G11M(i-1,msb),RGV1G11M(i,msb)))
     &                                 *X1U(I,MSA)
     &   +CONJG(FFM(RGV1G11M(i-1,msb)))*X1U(I-1,MSA)
     &   +CONJG(FFP(RGV1G11M(i,msb)))  *X1U(I+1,MSA)
     &   +CONJG(FGM(RGV1G12(i,msb)*zv2m,RGV1G12M(i-1,msb)))
     &                                 *X2U(I-1,MSA)
     &   +CONJG(FGP(RGV1G12(i,msb)*zv2p,RGV1G12M(i,msb)))
     &                                 *X2U(I,MSA)
 70   CONTINUE

      I = NRP1
      ZV2M = (CS(I)/CSM(I-1))**IEXV2
      INCLUDE 'tent.inc'

      R (KXV1+(MSMI-1)*NXCOMP,I) = R (KXV1+(MSMI-1)*NXCOMP,I)
     &   +CONJG(FF(RGV1G11(i,msb),RGV1G11M(i-1,msb),RGV1G11M(i,msb)))
     &                                 *X1U(I,MSA)
     &   +CONJG(FFM(RGV1G11M(i-1,msb)))*X1U(I-1,MSA)
     &   +CONJG(FGM(RGV1G12(i,msb)*zv2m,RGV1G12M(i-1,msb)))
     &                                 *X2U(I-1,MSA)

      DO 75 I=1,NR

      ZV2M = (CS(I  )/CSM(I))**IEXV2
      ZV2P = (CS(I+1)/CSM(I))**IEXV2
      ZEM  = (CS(I  )/CSM(I))**IEXE
      ZEP  = (CS(I+1)/CSM(I))**IEXE
      INCLUDE 'tophat.inc'

      RY(KYV2+(MSMI-1)*NYCOMP,I) = RY(KYV2+(MSMI-1)*NYCOMP,I)
     &  +CONJG(GG(RGV2G22M(i,msb),RGV2G22(i,msb)*zem*zv2m,
     $        RGV2G22(i+1,msb)*zep*zv2p))*X2U(I,MSA)
     &  + CONJG(GF(RGV1G12(i,msb)*zem,RGV1G12M(i,msb)))
     &                             *X1U(I,MSA)
     &  + CONJG(GF(RGV1G12(i+1,msb)*zep,RGV1G12M(i,msb)))
     &                             *X1U(I+1,MSA)
 75   CONTINUE
 80   CONTINUE
 90   CONTINUE

      IF (NV.GE.2) GOTO 200
      DO 110 MSA=1,MSMAX
      R(KXV1+(MSA-1)*NXCOMP,NRP1) = 0.
 110  CONTINUE
 200  CONTINUE

      DO 220 MSA=1,MSMAX
      DO 210 I=1,NFIT
 210  R(KXV1+(MSA-1)*NXCOMP,I) = 0.
 220  CONTINUE
      
      DO MSA=1,MSMAX
         LXCOL = (MSA -1)*NXCOMP
         LYCOL = (MSA -1)*NYCOMP
         DO I=1,NRP1
            DWX2X(I)=DWX2X(I) + CONJG(X1U(I,MSA))*R(LXCOL+KXV1,I)
         ENDDO
         DO I=1,NR
            DWX2Y(I)=DWX2Y(I) + CONJG(X2U(I,MSA))*RY(LYCOL+KYV2,I)
         ENDDO
      ENDDO

      DO I=1,NR
      WRITE(*,1002) CSM(I),DWPPARY(I),DWPPERY(I),DWPPARX(I),DWPPERX(I)
      ENDDO

      PI2 = PI*PI*2.0
      DO I=1,NR
      DWPY(I)    = PI2*(DWPY(I)+(DWPX(I)+DWPX(I+1))*0.5)*CSH(I)
      DWPPARY(I) = PI2*(DWPPARY(I)+(DWPPARX(I)+DWPPARX(I+1))*0.5)*CSH(I)
      DWPPERY(I) = PI2*(DWPPERY(I)+(DWPPERX(I)+DWPPERX(I+1))*0.5)*CSH(I)
      DWQY(I)    = PI2*(DWQY(I)+(DWQX(I)+DWQX(I+1))*0.5)*CSH(I)
      DWVY(I)    = PI2*(DWVY(I)+(DWVX(I)+DWVX(I+1))*0.5)*CSH(I)
      DWRHOY(I)  = PI2*(DWRHOY(I)+(DWRHOX(I)+DWRHOX(I+1))*0.5)*CSH(I)
      DWXIY(I)   = PI2*(DWXIY(I)+(DWXIX(I)+DWXIX(I+1))*0.5)*CSH(I)
      DWX2Y(I)   = PI2*(DWX2Y(I)+(DWX2X(I)+DWX2X(I+1))*0.5)*CSH(I)
      IF (KEFORM.EQ.1) 
     &DWJY(I)    = PI2*(DWJY(I)+(DWJX(I)+DWJX(I+1))*0.5)*CSH(I)
      ENDDO

      IF (KEFORM.EQ.2) CALL KDWFMAGP(DWJY)
      IF (KEFORM.EQ.2) CALL KDWFMAGV(DWVAC)

      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) THEN
            DWPY(I)    = 0.
            DWPPARY(I) = 0.
            DWPPERY(I) = 0.
            DWJY(I)    = 0.
            DWQY(I)    = 0.
            DWVY(I)    = 0.
            DWRHOY(I)  = 0.
            DWXIY(I)   = 0.
            DWX2Y(I)   = 0.
         ENDIF
      ENDDO

C     COMPUTE INTEGRATED ENERGY OVER RADIUS
      DWP    = 0.
      DWPPAR = 0.
      DWPPER = 0.
      DWJ    = 0.
      DWQ    = 0.
      DWV    = 0.
      DWRHO  = 0.
      DWXI   = 0.
      DWX2   = 0.

      DO I=NFIT+1,NR-0
         DWP    = DWP - DWPY(I)
         DWPPAR = DWPPAR - DWPPARY(I)
         DWPPER = DWPPER - DWPPERY(I)
         DWJ    = DWJ - DWJY(I)
         DWQ    = DWQ - DWQY(I)
         DWV    = DWV - DWVY(I)
         DWRHO  = DWRHO - DWRHOY(I)
         DWXI   = DWXI - DWXIY(I)
         DWX2   = DWX2 + DWX2Y(I)
      ENDDO

C     OUTPUT
      DWALL = DWP+DWPPAR+DWPPER+DWJ+DWQ+DWV+DWXI+DWRHO
      WRITE(*,*) 'PERTURBED ENERGY: DWP DWPPAR DWPPER DWJ' 
      WRITE(*,1001) DWP,DWPPAR,DWPPER,DWJ
      WRITE(*,*) 'PERTURBED ENERGY: DWQ DWV DWXI DWX2 DWALL'
      WRITE(*,1001) DWQ,DWV,DWXI,DWX2,DWALL
      WRITE(*,*) 'SAVE ENERGY DISTRIBUTION TO <EPLASMA>'
      WRITE(*,*) 'CSM DWX2 DWJ+DWQ+DQP DWPPER+DWPPAR DWALL'
      OPEN(CHLIST,FILE='EPLASMA.OUT',FORM='FORMATTED')
      DO I=1,NR
         DWALL = DWPY(I)+DWPPARY(I)+DWPPERY(I)+DWJY(I)+
     &           DWQY(I)+DWVY(I)+DWRHOY(I)+DWXIY(I)
         CTMP = DWX2*CSM(I)*2.0*CSH(I)
         CTMP = 1.0
         WRITE(CHLIST,1002) CSM(I), CSH(I), DWX2Y(I)/CTMP,
     &                 -(DWPY(I))/CTMP,
C    &                 -(DWJY(I)+DWQY(I)+DWPY(I))/CTMP,
     &                 -DWPPARY(I)/CTMP,
     &                 -DWPPERY(I)/CTMP
      ENDDO
      CLOSE(CHLIST)
 1001 FORMAT(5(E11.4,1X,E11.4,2X))
 1002 FORMAT(E12.5,1X,E12.5,2X,4(E11.4,1X,E11.4,2X))

      OPEN(CHLIST,FILE='ENERGY.OUT',FORM='FORMATTED',POSITION='APPEND')
C      IF (NPARAM.LE.1) GOTO 40
C      DO 30 I = 2,NPARAM
C      READ(CHLIST,'(A)',END=40) LINE
C 30   CONTINUE
C 40   CONTINUE
      WRITE(CHLIST,1045) DWX2,DWJ+DWQ+DWP,DWVAC,DWP,DWPPAR+DWPPER
C    &                  ,DWRHO,DWV,DWXI,DWPPAR,DWPPER
C    &                  ,DWJ+DWQ+DWP+DWPPAR+DWPPER+DWV+DWXI+DWRHO
      WRITE(CHLIST,*)
      CLOSE(CHLIST)
 1045 FORMAT(E13.5,$)


      DEALLOCATE(DWPX,DWPPARX,DWPPERX,DWJX,DWQX,DWVX,DWXIX,DWX2X,DWRHOX,
     &           DWPY,DWPPARY,DWPPERY,DWJY,DWQY,DWVY,DWXIY,DWX2Y,DWRHOY)

      DEALLOCATE( R, RY )

      RETURN
      END

      FUNCTION DIFPI(X)
      IMPLICIT NONE
      REAL*8::DIFPI,X

      DIFPI = X
      IF (X.LT.0.0) DIFPI = X + 2.0*ACOS(-1.)
      IF (DIFPI.LT.0.0) STOP 'DIFPI<0!'           
      IF (DIFPI.LT.1.0e-13) DIFPI=1.0e-13             

      RETURN
      END

      SUBROUTINE FILLMATDWKCOMP(I,INDX,
     &                          ASUBM,BSUBM,CSUBM,DSUBM,
     &                          ESUBM,FSUBM,GSUBM,HSUBM)
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE GIJLM
      USE CONVOLCOFM
      IMPLICIT NONE
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      INCLUDE 'integc.inc'
            
      INTEGER I,INDX,J
      INTEGER MROW,MSA

      INTEGER LXCOL,LYCOL,LXROW,LYROW
      REAL*8  ZV2M,ZV2P,ZB3M,ZB3P
      INTEGER   IEXV2,IEXB3            
      PARAMETER (IEXV2=-1, IEXB3=1)      

      COMPLEX*16,DIMENSION(:,:),POINTER::VX1PARAI,VX1PARAMI,VX1PARAI1,
     &                                   VX1PERPI,VX1PERPMI,VX1PERPI1,
     &                                   VX2PARAI,VX2PARAMI,VX2PARAI1,
     &                                   VX2PERPI,VX2PERPMI,VX2PERPI1,
     &                                   VQ1PARAI,VQ1PARAMI,VQ1PARAI1,
     &                                   VQ1PERPI,VQ1PERPMI,VQ1PERPI1,
     &                                   VQ2PARAI,VQ2PARAMI,VQ2PARAI1,
     &                                   VQ2PERPI,VQ2PERPMI,VQ2PERPI1,
     &                                   VQ3PARAI,VQ3PARAMI,VQ3PARAI1,
     &                                   VQ3PERPI,VQ3PERPMI,VQ3PERPI1
      
      VX1PARAI => BUFFER_I(:,:,INDX,1)
      VX1PERPI => BUFFER_I(:,:,INDX,2)
      VX2PARAI => BUFFER_I(:,:,INDX,3)
      VX2PERPI => BUFFER_I(:,:,INDX,4)
      VQ1PARAI => BUFFER_I(:,:,INDX,5)
      VQ1PERPI => BUFFER_I(:,:,INDX,6)
      VQ2PARAI => BUFFER_I(:,:,INDX,7)
      VQ2PERPI => BUFFER_I(:,:,INDX,8)
      VQ3PARAI => BUFFER_I(:,:,INDX,9)
      VQ3PERPI => BUFFER_I(:,:,INDX,10)
      
      VX1PARAMI => BUFFERM_I(:,:,INDX,1)
      VX1PERPMI => BUFFERM_I(:,:,INDX,2)
      VX2PARAMI => BUFFERM_I(:,:,INDX,3)
      VX2PERPMI => BUFFERM_I(:,:,INDX,4)
      VQ1PARAMI => BUFFERM_I(:,:,INDX,5)
      VQ1PERPMI => BUFFERM_I(:,:,INDX,6)
      VQ2PARAMI => BUFFERM_I(:,:,INDX,7)
      VQ2PERPMI => BUFFERM_I(:,:,INDX,8)
      VQ3PARAMI => BUFFERM_I(:,:,INDX,9)
      VQ3PERPMI => BUFFERM_I(:,:,INDX,10)
      
      VX1PARAI1 => BUFFER_I1(:,:,INDX,1)
      VX1PERPI1 => BUFFER_I1(:,:,INDX,2)
      VX2PARAI1 => BUFFER_I1(:,:,INDX,3)
      VX2PERPI1 => BUFFER_I1(:,:,INDX,4)
      VQ1PARAI1 => BUFFER_I1(:,:,INDX,5)
      VQ1PERPI1 => BUFFER_I1(:,:,INDX,6)
      VQ2PARAI1 => BUFFER_I1(:,:,INDX,7)
      VQ2PERPI1 => BUFFER_I1(:,:,INDX,8)
      VQ3PARAI1 => BUFFER_I1(:,:,INDX,9)
      VQ3PERPI1 => BUFFER_I1(:,:,INDX,10)
      
C     FILL IN MATRICES FOR EQUATIONS FOR PARALLEL AND 
C     PERPENDICULAR KINETIC PRESSURE
         INCLUDE 'tophat.inc'
         ZV2M = (CS(I  )/CSM(I))**IEXV2
         ZV2P = (CS(I+1)/CSM(I))**IEXV2
         ZB3M = (CS(I  )/CSM(I))**IEXB3
         ZB3P = (CS(I+1)/CSM(I))**IEXB3

         DO MROW=1,MSMAX
            LYROW = (MROW-1)*NYCOMP
            DO MSA=1,MSMAX
               LXCOL = (MSA-1)*NXCOMP
               LYCOL = (MSA-1)*NYCOMP

               FSUBM(KYPPARA+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PARAI(MROW,MSA),VX1PARAMI(MROW,MSA))

               FSUBM(KYPPERP+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PERPI(MROW,MSA),VX1PERPMI(MROW,MSA))

               GSUBM(KYPPARA+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PARAI1(MROW,MSA),VX1PARAMI(MROW,MSA))

               GSUBM(KYPPERP+LYROW, KXX1+LXCOL,I)=
     &           GF(VX1PERPI1(MROW,MSA),VX1PERPMI(MROW,MSA))

               DSUBM(KYPPARA+LYROW, KYX2+LYCOL,I)=
     &           GG(VX2PARAMI(MROW,MSA), VX2PARAI(MROW,MSA)*ZV2M,
     &           VX2PARAI1(MROW,MSA)*ZV2P)

               DSUBM(KYPPERP+LYROW, KYX2+LYCOL,I)=
     &           GG(VX2PERPMI(MROW,MSA), VX2PERPI(MROW,MSA)*ZV2M,
     &           VX2PERPI1(MROW,MSA)*ZV2P)

               FSUBM(KYPPARA+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PARAI(MROW,MSA), VQ1PARAMI(MROW,MSA))    

               FSUBM(KYPPERP+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PERPI(MROW,MSA), VQ1PERPMI(MROW,MSA))    

               GSUBM(KYPPARA+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PARAI1(MROW,MSA), VQ1PARAMI(MROW,MSA))    

               GSUBM(KYPPERP+LYROW, KXB1+LXCOL,I)=
     &           GF(VQ1PERPI1(MROW,MSA), VQ1PERPMI(MROW,MSA))    

               DSUBM(KYPPARA+LYROW, KYB2+LYCOL,I)=
     &           GG(VQ2PARAMI(MROW,MSA), VQ2PARAI(MROW,MSA),
     &           VQ2PARAI1(MROW,MSA))

               DSUBM(KYPPERP+LYROW, KYB2+LYCOL,I)=
     &           GG(VQ2PERPMI(MROW,MSA), VQ2PERPI(MROW,MSA),
     &           VQ2PERPI1(MROW,MSA))

               DSUBM(KYPPARA+LYROW, KYB3+LYCOL,I)=
     &           GG(VQ3PARAMI(MROW,MSA), VQ3PARAI(MROW,MSA)*ZB3M,
     &           VQ3PARAI1(MROW,MSA)*ZB3P)

               DSUBM(KYPPERP+LYROW, KYB3+LYCOL,I)=
     &           GG(VQ3PERPMI(MROW,MSA), VQ3PERPI(MROW,MSA)*ZB3M,
     &           VQ3PERPI1(MROW,MSA)*ZB3P)
            ENDDO
         ENDDO
      END SUBROUTINE FILLMATDWKCOMP
      
      SUBROUTINE CALCPRECOMP(I,PPARAC,PPERPC,
     &                       ASUBM,BSUBM,CSUBM,DSUBM,
     &                       ESUBM,FSUBM,GSUBM,HSUBM)
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE GIJLM
      USE CONVOLCOFM
      IMPLICIT NONE
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'

      INTEGER    I,J,MROW,MSA,MSB,MS
      INTEGER    LXCOL,LYCOL,LXROW,LYROW
      REAL*8     HCHI
      COMPLEX*16 CTMP1
      COMPLEX*16,DIMENSION(NRP1,MSMAX)::    PPARAC,PPERPC
      COMPLEX*16,DIMENSION(:),ALLOCATABLE:: JPPARA,JPPERP
      COMPLEX*16,DIMENSION(:),ALLOCATABLE:: CTMP2,CTMP3
      
      ALLOCATE (JPPARA(MSMAX),JPPERP(MSMAX),CTMP2(NCHI),CTMP3(NCHI))
      JPPARA=0.0
      JPPERP=0.0

      DO MROW=1,MSMAX
      DO MSA=1,MSMAX
      LXCOL = (MSA -1)*NXCOMP
      LYROW = (MROW-1)*NYCOMP
      LYCOL = (MSA -1)*NYCOMP
C     COMPUTATION OF JPPARA AND JPPERP	  
      JPPARA(MROW)=JPPARA(MROW)+
     &             FSUBM(KYPPARA+LYROW,KXX1+LXCOL,I)*X1U(I,MSA) +
     &             GSUBM(KYPPARA+LYROW,KXX1+LXCOL,I)*X1U(I+1,MSA)+
     &             FSUBM(KYPPARA+LYROW,KXB1+LXCOL,I)*B1U(I,MSA) +
     &             GSUBM(KYPPARA+LYROW,KXB1+LXCOL,I)*B1U(I+1,MSA)+
     &             DSUBM(KYPPARA+LYROW,KYX2+LYCOL,I)*X2U(I,MSA)+
     &             DSUBM(KYPPARA+LYROW,KYB2+LYCOL,I)*B2U(I,MSA)+
     &             DSUBM(KYPPARA+LYROW,KYB3+LYCOL,I)*B3U(I,MSA)
      JPPERP(MROW)=JPPERP(MROW)+
     &             FSUBM(KYPPERP+LYROW,KXX1+LXCOL,I)*X1U(I,MSA) +
     &             GSUBM(KYPPERP+LYROW,KXX1+LXCOL,I)*X1U(I+1,MSA)+
     &             FSUBM(KYPPERP+LYROW,KXB1+LXCOL,I)*B1U(I,MSA) +
     &             GSUBM(KYPPERP+LYROW,KXB1+LXCOL,I)*B1U(I+1,MSA)+
     &             DSUBM(KYPPERP+LYROW,KYX2+LYCOL,I)*X2U(I,MSA)+
     &             DSUBM(KYPPERP+LYROW,KYB2+LYCOL,I)*B2U(I,MSA)+
     &             DSUBM(KYPPERP+LYROW,KYB3+LYCOL,I)*B3U(I,MSA)					
      ENDDO
      ENDDO
	  
C     COMPUTATION OF PPERP AND PPARA  
      CTMP2=0.0
      CTMP3=0.0
      HCHI = 2.*PI/NCHI
      DO J=1,NCHI
         DO MS=1,MSMAX
            CTMP1 = EXP(CI*RM(MS,2)*(J-1)*HCHI)
            CTMP2(J) = CTMP2(J) + JPPARA(MS)*CTMP1
            CTMP3(J) = CTMP3(J) + JPPERP(MS)*CTMP1
         ENDDO
         CTMP2(J) = CTMP2(J) / RJAM(I,J)
         CTMP3(J) = CTMP3(J) / RJAM(I,J)
      ENDDO
      
      DO MS=1,MSMAX
         DO J=1,NCHI
            CTMP1 = EXP(-CI*RM(MS,2)*(J-1)*HCHI)
            PPARAC(I,MS) = PPARAC(I,MS) + CTMP2(J)*CTMP1
            PPERPC(I,MS) = PPERPC(I,MS) + CTMP3(J)*CTMP1
         ENDDO
      ENDDO
      HCHI = 1./DFLOAT(NCHI)
      PPARAC(I,:)=PPARAC(I,:)*HCHI
      PPERPC(I,:)=PPERPC(I,:)*HCHI

      DEALLOCATE (JPPARA,JPPERP,CTMP2,CTMP3)
      END SUBROUTINE CALCPRECOMP
      
      SUBROUTINE CALCDWKPROF (IS,PPARAC,PPERPC,
     &                        DWPPERX,DWPPERY,DWPPARX,DWPPARY,
     &                        ASUBM,BSUBM,CSUBM,DSUBM,
     &                        ESUBM,FSUBM,GSUBM,HSUBM)
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE GIJLM
      USE CONVOLCOFM
      IMPLICIT NONE
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'

      INTEGER IS
      COMPLEX*16,DIMENSION(NRP1,MSMAX)::PPARAC,PPERPC
      COMPLEX*16,DIMENSION(NRP1)::DWPPARX,DWPPERX,
     &                                     DWPPARY,DWPPERY      
      INTEGER    MROW,MSA,MSB,MSMI,MSPL,NSA,NSB,I
      PARAMETER  (NSA=2,NSB=1)
      INTEGER    LXCOL,LYCOL,LXROW,LYROW
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::R,RY
      INTEGER       IEXV2,IEXY
      PARAMETER     (IEXV2 = -1, IEXY=0)
      REAL*8        ZEM,ZEP,ZV2M,ZV2P,z3m,z3p
      
      DO MROW=1,MSMAX
         LXROW = (MROW-1)*NXCOMP
         LYROW = (MROW-1)*NYCOMP
         DO MSA=1,MSMAX
            LXCOL = (MSA -1)*NXCOMP
            LYCOL = (MSA -1)*NYCOMP
            I=IS+1
            IF (KEFORM.EQ.1) THEN
            IF (KYPPERP.GT.0)
     &         DWPPERX(I)=DWPPERX(I) + CONJG(X1U(I,MROW))*
     &         HSUBM(LXROW+KXV1,LYCOL+KYPPERP,I)*PPERPC(I-1,MSA)
            ENDIF
            IF (KEFORM.EQ.2.AND.INCKIN.GT.0) THEN
               DWPPERX(I)=DWPPERX(I) - CONJG(X1U(I,MROW))*
     &         HSUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*PPERPC(I-1,MSA)
            ENDIF
            IF (KYPPARA.GT.0)
     &         DWPPARX(I)=DWPPARX(I) + CONJG(X1U(I,MROW))*
     &         HSUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*PPARAC(I-1,MSA)
            
            I=IS
            ZEM = CSM(I)**IEXY
            IF (KEFORM.EQ.1) THEN
            IF (KYPPERP.GT.0)
     &         DWPPERY(I)=DWPPERY(I) + ZEM*CONJG(X2U(I,MROW))*
     &         DSUBM(LYROW+KYV2,LYCOL+KYPPERP,I)*PPERPC(I,MSA)
            ENDIF
            IF (KEFORM.EQ.2) THEN
            IF (INCKIN.GT.0) THEN
               DWPPERY(I)=DWPPERY(I) - ZEM*CONJG(X2U(I,MROW))*
     &         DSUBM(LYROW+KYV2,LYCOL+KYPPARA,I)*PPERPC(I,MSA)

               DWPPERY(I)=DWPPERY(I) + ZEM*PPERPC(I,MROW)*CONJG(
     &         FSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I,MSA)+
     &         GSUBM(LYROW+KYPR,LXCOL+KXV1,I)*X1U(I+1,MSA)+
     &         DSUBM(LYROW+KYPR,LYCOL+KYV2,I)*X2U(I,MSA)) 
            ENDIF
            ENDIF
            IF (KYPPARA.GT.0) 
     &         DWPPARY(I)=DWPPARY(I) + ZEM*CONJG(X2U(I,MROW))*
     &         DSUBM(LYROW+KYV2,LYCOL+KYPPARA,I)*PPARAC(I,MSA)
            IF (KEFORM.EQ.1) THEN
            IF (KYPPERP.GT.0)
     &         DWPPERX(I)=DWPPERX(I) + CONJG(X1U(I,MROW))*
     &         ESUBM(LXROW+KXV1,LYCOL+KYPPERP,I)*PPERPC(I,MSA)
            ENDIF
            IF (KEFORM.EQ.2.AND.INCKIN.GT.0) THEN
               DWPPERX(I)=DWPPERX(I) - CONJG(X1U(I,MROW))*
     &         ESUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*PPERPC(I,MSA)
            ENDIF
            IF (KYPPARA.GT.0) 
     &         DWPPARX(I)=DWPPARX(I) + CONJG(X1U(I,MROW))*
     &         ESUBM(LXROW+KXV1,LYCOL+KYPPARA,I)*PPARAC(I,MSA)
         ENDDO
      ENDDO
      END SUBROUTINE CALCDWKPROF
      
      SUBROUTINE CALCDWKCOMP(ASUBM,BSUBM,CSUBM,DSUBM,
     $                       ESUBM,FSUBM,GSUBM,HSUBM)
      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE GIJLM
      USE CONVOLCOFM
      USE TORQUEM
      USE ToolBox
      IMPLICIT NONE
      INCLUDE 'specmat.inc'
      INCLUDE 'compam.inc'
      INCLUDE 'comioc.inc'
      
      INTEGER KP,IS,TOTINDX,INDX,I,J,FID
      REAL*8 PI2
      COMPLEX*16,DIMENSION(:),ALLOCATABLE:: DWPPARA,DWPPERP,DWK
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::DWPPARX,DWPPERX,
     &                                       DWPPARY,DWPPERY
      COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE:: PPARAC,PPERPC
      COMPLEX*16,DIMENSION(:,:,:,:),ALLOCATABLE,TARGET:: 
     &                  BUFFER_DATA1,BUFFER_DATA2,BUFFER_DATAM
      COMPLEX*16,DIMENSION(:,:,:,:),POINTER:: TMPPOT

      IF (INCKIN.NE.1) RETURN
      IF (.NOT. ODWKCOM) RETURN
      IF (ISWEEP.NE.NSWEEP) RETURN
      
      CALL ALLOCATEDWKCOMPMAT
      TOTINDX = SIZE(VX1PARAC,3)
      ALLOCATE (DWPPARA(TOTINDX),DWPPERP(TOTINDX),DWK(TOTINDX))
      ALLOCATE (DWPPARX(NRP1,TOTINDX),DWPPERX(NRP1,TOTINDX),
     &          DWPPARY(NRP1,TOTINDX),DWPPERY(NRP1,TOTINDX) )
      ALLOCATE (PPARAC(NRP1,MSMAX,TOTINDX),PPERPC(NRP1,MSMAX,TOTINDX))
      ALLOCATE (BUFFER_DATA1(MSMAX,MSMAX,TOTINDX,30),
     &          BUFFER_DATA2(MSMAX,MSMAX,TOTINDX,30),
     &          BUFFER_DATAM(MSMAX,MSMAX,TOTINDX,30))
      DWPPARX = 0.0
      DWPPERX = 0.0
      DWPPARY = 0.0
      DWPPERY = 0.0
      
      PPARAC  = 0.0
      PPERPC  = 0.0
      
      BUFFER_I  => BUFFER_DATA1
      BUFFER_I1 => BUFFER_DATA2
      BUFFERM_I => BUFFER_DATAM
      
      BUFFERT => BUFFER_I1
      CALL READ_SURFACE_QUANTITIES (1,2)
      DO IS=1,NRP1-1
C     READ MARTIX FROM BINARY FILES ON EACH SURFACE
         BUFFERT => BUFFERM_I
         CALL READ_SURFACE_QUANTITIES (IS,2)
         
         TMPPOT => BUFFER_I
         BUFFER_I => BUFFER_I1
         BUFFER_I1 => TMPPOT
         IF (IS+1.EQ.NRP1) THEN
            BUFFER_I1 = BUFFERM_I
            
         ELSE
            BUFFERT => BUFFER_I1
            CALL READ_SURFACE_QUANTITIES (IS+1,1)
         ENDIF   
         
         DO INDX=1,TOTINDX
C     FILL IN THE GLOBAL MATRIX FOR PRESSURE CALCULATION
            CALL FILLMATDWKCOMP(IS,INDX,
     &                          ASUBM,BSUBM,CSUBM,DSUBM,
     &                          ESUBM,FSUBM,GSUBM,HSUBM)         
C     CALCULATE THE COMPONENTS OF PRESSURE     
            CALL CALCPRECOMP(IS,PPARAC(:,:,INDX),PPERPC(:,:,INDX),
     &                       ASUBM,BSUBM,CSUBM,DSUBM,
     &                       ESUBM,FSUBM,GSUBM,HSUBM)
C     CALCULATE ENERGY PROFILE OF DIFFERENT COMPONENTS
            CALL CALCDWKPROF (IS,PPARAC(:,:,INDX),PPERPC(:,:,INDX),
     &                        DWPPERX(:,INDX),DWPPERY(:,INDX),
     &                        DWPPARX(:,INDX),DWPPARY(:,INDX),
     &                        ASUBM,BSUBM,CSUBM,DSUBM,
     &                        ESUBM,FSUBM,GSUBM,HSUBM)
	  
         ENDDO
      ENDDO

      PI2 = PI*PI*2.0
C     CALCULATE THE TOTAL ENERGY OF EACH COMPONENT
      DO IS=1,NR
         DWPPARY(IS,:) = 
     &   PI2*(DWPPARY(IS,:)+(DWPPARX(IS,:)+DWPPARX(IS+1,:))*0.5)
     
         DWPPERY(IS,:) = 
     &   PI2*(DWPPERY(IS,:)+(DWPPERX(IS,:)+DWPPERX(IS+1,:))*0.5)
      ENDDO
      
      DO IS=1,NR
         IF (CSM(IS).GT.CTEDGE) THEN
            DWPPARY(IS,:)    = 0.
            DWPPERY(IS,:)    = 0.
         ENDIF
      ENDDO

      DWPPARA=0.0
      DWPPERP=0.0
      DWK   =0.0
      
      DO IS=NFIT+1,NR
         DWPPARA = DWPPARA - DWPPARY(IS,:)*CSH(IS)
         DWPPERP = DWPPERP - DWPPERY(IS,:)*CSH(IS)      
      ENDDO
      DWK = DWPPARA+DWPPERP

      IF (KNTV.EQ.21) THEN
C     COMPUTE NTV TORQUE DENSITY FROM DWK
C     FOR IONS AND ELECTRONS SEPARATELY
C     USING: T_NTV = -2*N*IM(DWKA)/(4*PI^2)
      TORQUENTVI = 0.
      TORQUENTVE = 0.
      DO IS=1,NR
      DO I=1,5
         INDX=INDXDWKC(1,I)
         IF (INDX.LT.0) CYCLE
         TORQUENTVI(IS)=TORQUENTVI(IS)+
     &                  IMAG(-DWPPARY(IS,INDX)-DWPPERY(IS,INDX))    

         INDX=INDXDWKC(2,I)
         IF (INDX.LT.0) CYCLE
         TORQUENTVE(IS)=TORQUENTVE(IS)+
     &                  IMAG(-DWPPARY(IS,INDX)-DWPPERY(IS,INDX))    
      ENDDO
      ENDDO
      TORQUENTVI = -2.*RNTOR*TORQUENTVI/(4.*PI*PI)
      TORQUENTVE = -2.*RNTOR*TORQUENTVE/(4.*PI*PI)
      TORQUENTV  = TORQUENTVI + TORQUENTVE
      ENDIF

C     OUTPUT THE PROFILES OF ENERGY DENSITY
      FID=ASSIGNFREEFILEUNIT () 
      OPEN(FID,FILE='DWK_ENERGY_DENSITY.OUT',FORM='FORMATTED',
     &     STATUS='REPLACE')
      DO IS=1,NR
      DO KP=1,NSPECIES
      DO I=1,5
         INDX=INDXDWKC(KP,I)
         IF (INDX.LT.0) CYCLE
         WRITE(FID,1050) KP,I,CSM(IS), CSH(IS),
     &                   -DWPPARY(IS,INDX),-DWPPERY(IS,INDX),
     &                   -(DWPPARY(IS,INDX)+DWPPERY(IS,INDX))
      ENDDO
      ENDDO
      WRITE(FID,1050) 0,0,CSM(IS), CSH(IS),
     &                -SUM(DWPPARY(IS,:)),-SUM(DWPPERY(IS,:)),
     &                -(SUM(DWPPARY(IS,:))+SUM(DWPPERY(IS,:)))
      WRITE(FID,*)
      ENDDO
1050  FORMAT (2I5,8(E13.5),$)

      WRITE(FID,*)'%MATRIX INSTRUCTION'
      WRITE(FID,*)'%EVERY 10 COLUMNS IS A GROUP.'
      WRITE(FID,*)'%THE FORMAT OF EACH GROUP IS'
      WRITE(FID,*)'%KP I CSM CSH DWKPARA_RE DWKPARA_IM '//
     &'DWKPERP_RE DWKPERP_IM DWK_RE DWK_IM'
      WRITE(FID,*)'%KP: THE PARTICLE SPECIE FROM 1 TO NSPECIES'
      WRITE(FID,*)'%    KP=0 INCLUDES ALL PARTICLE SPECIES'
      WRITE(FID,*)'%I=0 TOTAL KINETIC EFFECT'
      WRITE(FID,*)'%  1 ADIABATIC PART OF PASSING PARTICLES
     & (PSPECIES_AP)'
      WRITE(FID,*)'%  2 ADIABATIC PART OF TRAPPED PARTICLES
     & (PSPECIES_TP)'
      WRITE(FID,*)'%  3 NON-ADIABATIC PART OF PASSING
     & PARTICLES (PSPECIES_NP)'
      WRITE(FID,*)'%  4 NON-ADIABATIC PART OF BOUNCE
     & RESONANCE (PSPECIES_NTB)'
      WRITE(FID,*)'%  5 NON-ADIABATIC PART OF PRECESSION
     & RESONACE (PSPECIES_NTB)'
      CLOSE(FID)
      

C     OUTPUT THE ENERGY COMPONENTS      
      FID=ASSIGNFREEFILEUNIT () 
      OPEN(FID,FILE='DWK_COMPONENTS.OUT',FORM='FORMATTED',
     &     STATUS='REPLACE')
      DO KP=1,NSPECIES
      DO I=1,5
         INDX=INDXDWKC(KP,I)
         IF (INDX.LT.0) CYCLE
         WRITE(FID,1100) KP,I,DWPPARA(INDX),DWPPERP(INDX),DWK(INDX)
      ENDDO
      ENDDO
      WRITE(FID,1100) 0,0,SUM(DWPPARA),SUM(DWPPERP),SUM(DWK)
      WRITE(FID,*)
      WRITE(FID,*)
      WRITE(FID,*)'%MATRIX INSTRUCTION'
      WRITE(FID,*)'%COLUMN (1): THE PARTICLE SPECIE FROM KP=1 TO
     & NSPECIES'
      WRITE(FID,*)'%            KP=0 INCLUDES ALL PARTICLE SPECIES'
      WRITE(FID,*)'%COLUMN (2): THE CONTIRBUTION DUE TO DIFFERENT
     & KINETIC EFFECTS'
      WRITE(FID,*)'%            I=0 TOTAL KINETIC EFFECT'
      WRITE(FID,*)'%            I=1 ADIABATIC PART OF PASSING PARTICLES
     & (PSPECIES_AP)'
      WRITE(FID,*)'%            I=2 ADIABATIC PART OF TRAPPED PARTICLES
     & (PSPECIES_TP)'
      WRITE(FID,*)'%            I=3 NON-ADIABATIC PART OF PASSING
     & PARTICLES (PSPECIES_NP)'
      WRITE(FID,*)'%            I=4 NON-ADIABATIC PART OF BOUNCE
     & RESONANCE (PSPECIES_NTB)'
      WRITE(FID,*)'%            I=5 NON-ADIABATIC PART OF PRECESSION
     & RESONACE (PSPECIES_NTD)'
      WRITE(FID,*)'%COLUMN (3,4): DWPPARA(KP,I)'
      WRITE(FID,*)'%COLUMN (5,6): DWPPERP(KP,I)'
      WRITE(FID,*)'%COLUMN (7,8): DWK(KP,I)=DWPPARA(KP,I)+DWPPERP(KP,I)'
      WRITE(FID,*)'%NOTE:INTEGRAL OF ENERGY DENSITY IGNORES FIRST NFIT+1
     & POINTS'
      CLOSE(FID)
1100  FORMAT(2I5,6E13.5)

      WRITE (*,*) SUM(DWPPARA)
      WRITE (*,*) SUM(DWPPERP)
      DEALLOCATE (DWPPARA,DWPPERP,DWK)
      DEALLOCATE (DWPPARX,DWPPERX,DWPPARY,DWPPERY)
      DEALLOCATE (PPARAC,PPERPC)
      DEALLOCATE (BUFFER_DATA1,BUFFER_DATA2,BUFFER_DATAM)

      CALL DEALLOCATEDWKCOMPMAT

      END SUBROUTINE CALCDWKCOMP
      
      SUBROUTINE WRITE_SURFACE_QUANTITIES (IS,KGRID)
      USE KINETICM
      USE GLOBALM
      USE ToolBox
      
      IMPLICIT NONE
      INTEGER IS,KGRID,FID
      CHARACTER*50 FILENAME
      INTEGER SIZEOFARRAY
      
      IF (.NOT. ODWKCOM) RETURN
      IF (ISWEEP.NE.NSWEEP) RETURN
      WRITE(FILENAME,*) IS
      IF (KGRID.EQ.1) FILENAME='PCOEF'//TRIM(ADJUSTL(FILENAME))//'.DAT'
      IF (KGRID.EQ.2) FILENAME='PCOEF'//TRIM(ADJUSTL(FILENAME))//'M.DAT'
C$OMP CRITICAL(WRITE_FILE)        
      FID=ASSIGNFREEFILEUNIT () 
      SIZEOFARRAY = SIZEOF(VX1PARAC)
      FILENAME='./DATA_STORAGE/'//TRIM(ADJUSTL(FILENAME))
      OPEN (UNIT=FID, ACCESS='SEQUENTIAL',FILE=TRIM(FILENAME),
     &FORM='UNFORMATTED',STATUS='REPLACE',ACTION ='WRITE')
      REWIND(FID)         
      WRITE (FID) VX1PARAC*ALPHAD
      WRITE (FID) VX1PERPC*ALPHAD
      WRITE (FID) VX2PARAC*ALPHAD
      WRITE (FID) VX2PERPC*ALPHAD
      WRITE (FID) VQ1PARAC*ALPHAD
      WRITE (FID) VQ1PERPC*ALPHAD
      WRITE (FID) VQ2PARAC*ALPHAD
      WRITE (FID) VQ2PERPC*ALPHAD
      WRITE (FID) VQ3PARAC*ALPHAD
      WRITE (FID) VQ3PERPC*ALPHAD
      FLUSH (FID)
      CLOSE (FID)
C$OMP END CRITICAL(WRITE_FILE)    
C      WRITE (*,*) 'SUCCESS OF WRITE:', FILENAME          
      END SUBROUTINE WRITE_SURFACE_QUANTITIES

      SUBROUTINE READ_SURFACE_QUANTITIES (IS,KGRID)
      USE KINETICM      
      USE GLOBALM
      USE ToolBox
       
      IMPLICIT NONE
      INTEGER IS,KGRID,FID,IVAL
      CHARACTER*50 FILENAME
      INTEGER SIZEOFARRAY,TMPMX,TMPMY,TOTINDX
      COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE::READ_DATA
      TMPMX= SIZE(VX1PARAC,1)
      TMPMY= SIZE(VX1PARAC,2)
      TOTINDX = SIZE(VX1PARAC,3)
      ALLOCATE (READ_DATA(TMPMX,TMPMY,TOTINDX))
      IF (.NOT. ODWKCOM) RETURN
      IF (ISWEEP.NE.NSWEEP) RETURN
      WRITE(FILENAME,*) IS
      IF (KGRID.EQ.1) FILENAME='PCOEF'//TRIM(ADJUSTL(FILENAME))//'.DAT'
      IF (KGRID.EQ.2) FILENAME='PCOEF'//TRIM(ADJUSTL(FILENAME))//'M.DAT'
      FID=ASSIGNFREEFILEUNIT () 
      SIZEOFARRAY=SIZEOF(BUFFERT)
      FILENAME='./DATA_STORAGE/'//TRIM(ADJUSTL(FILENAME))
      OPEN (UNIT=FID, ACCESS='SEQUENTIAL',FILE=TRIM(FILENAME),
     &FORM='UNFORMATTED',STATUS='OLD',ACTION ='READ')
      REWIND(FID) 
      DO IVAL=1,10
         READ (FID) READ_DATA
         BUFFERT(:,:,:,IVAL)=READ_DATA
      ENDDO
      CLOSE (FID)  
      DEALLOCATE(READ_DATA)
C      WRITE (*,*) 'SUCCESS OF READ:', FILENAME          
      END SUBROUTINE READ_SURFACE_QUANTITIES      
      
      SUBROUTINE ALLOCATEDWKCOMPMAT
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      
      IMPLICIT NONE
      
      INTEGER KP,INDX
      
      IF (.NOT. ODWKCOM) RETURN
	
      ALLOCATE (INDXDWKC(NSPECIES,5))
      INDXDWKC = -1
      
      INDX=0
      DO KP=1,NSPECIES
         IF (ABS(PSPECIES_AP(KP)).GT.0) THEN
            INDX=INDX+1
            INDXDWKC(KP,1)=INDX
         ENDIF
         IF (ABS(PSPECIES_AT(KP)).GT.0) THEN
            INDX=INDX+1
            INDXDWKC(KP,2)=INDX
         ENDIF
           IF (ABS(PSPECIES_NP(KP)).GT.0) THEN
            INDX=INDX+1
            INDXDWKC(KP,3)=INDX
         ENDIF
         IF (ABS(PSPECIES_NTB(KP)).GT.0) THEN
            INDX=INDX+1
            INDXDWKC(KP,4)=INDX
         ENDIF
         IF (ABS(PSPECIES_NTD(KP)).GT.0) THEN
            INDX=INDX+1
            INDXDWKC(KP,5)=INDX
         ENDIF
      ENDDO
      
      ALLOCATE ( VX1PARAC(MSMAX,MSMAX,INDX), VX1PERPC(MSMAX,MSMAX,INDX),
     &           VX2PARAC(MSMAX,MSMAX,INDX), VX2PERPC(MSMAX,MSMAX,INDX),
     &           VQ1PARAC(MSMAX,MSMAX,INDX), VQ1PERPC(MSMAX,MSMAX,INDX),
     &           VQ2PARAC(MSMAX,MSMAX,INDX), VQ2PERPC(MSMAX,MSMAX,INDX),
     &           VQ3PARAC(MSMAX,MSMAX,INDX), VQ3PERPC(MSMAX,MSMAX,INDX))
      END SUBROUTINE ALLOCATEDWKCOMPMAT
      
      SUBROUTINE DEALLOCATEDWKCOMPMAT
      USE KINETICM
      USE GLOBALM
      
      IMPLICIT NONE
      
      IF (.NOT. ODWKCOM) RETURN
	IF (.NOT.ALLOCATED(INDXDWKC)) RETURN
	  
      DEALLOCATE ( INDXDWKC )
      DEALLOCATE (VX1PARAC, VX1PERPC, VX2PARAC, VX2PERPC,
     &            VQ1PARAC, VQ1PERPC, VQ2PARAC, VQ2PERPC, 
     &            VQ3PARAC, VQ3PERPC)
      END SUBROUTINE DEALLOCATEDWKCOMPMAT
      
      SUBROUTINE SETDWKCOMPMAT ( KP,COMP,K,M,
     &                           X1PARA,X1PERP,X2PARA,X2PERP,Q1PARA,
     &                           Q1PERP,Q2PARA,Q2PERP,Q3PARA,Q3PERP )
	USE KINETICM
      USE GLOBALM
	  
      IMPLICIT NONE
      
      INTEGER KP,COMP,K,M
      COMPLEX*16 X1PARA,X1PERP,X2PARA,X2PERP,Q1PARA,
     &           Q1PERP,Q2PARA,Q2PERP,Q3PARA,Q3PERP
      INTEGER INDX
      IF (.NOT. ODWKCOM) RETURN      
      IF (ISWEEP.NE.NSWEEP) RETURN
      IF (INDXDWKC(KP,COMP).EQ.-1) RETURN
      
      INDX=INDXDWKC(KP,COMP)
      
      VX1PARAC(K,M,INDX)=VX1PARAC(K,M,INDX)+X1PARA
      VX1PERPC(K,M,INDX)=VX1PERPC(K,M,INDX)+X1PERP
      VX2PARAC(K,M,INDX)=VX2PARAC(K,M,INDX)+X2PARA
      VX2PERPC(K,M,INDX)=VX2PERPC(K,M,INDX)+X2PERP
      VQ1PARAC(K,M,INDX)=VQ1PARAC(K,M,INDX)+Q1PARA
      VQ1PERPC(K,M,INDX)=VQ1PERPC(K,M,INDX)+Q1PERP
      VQ2PARAC(K,M,INDX)=VQ2PARAC(K,M,INDX)+Q2PARA
      VQ2PERPC(K,M,INDX)=VQ2PERPC(K,M,INDX)+Q2PERP
      VQ3PARAC(K,M,INDX)=VQ3PARAC(K,M,INDX)+Q3PARA
      VQ3PERPC(K,M,INDX)=VQ3PERPC(K,M,INDX)+Q3PERP
      
      END SUBROUTINE SETDWKCOMPMAT

	  
      function plasmaz(y)
c     z is defined as the standard plasma dispersion function  
c     according to Fried and Conte
      implicit none
      complex*16::plasmaz,y,y2,z
      integer::k
      real*8::w,w1,w2,ht,t,spi

      integer N0,N
      parameter (N0=3602)
      real*8 tt(N0+1),t0

      y2  = y*y
      spi = sqrt(acos(-1.))

      if (abs(y).lt.3.5) then
      w  = 1.0/sqrt(3.0)
      w1 = (1.-w)*0.5
      w2 = (1.+w)*0.5

      N = 101
      if (abs(y).gt.1.) 
     &   N = 101 + int((abs(y)-1.)**2*400)

      t0 = -4.0
      ht = (0.0-t0)/(N-1)
      do k=1,N
         tt(N+1-k) = 1.0 - 10**(t0+(k-1)*ht)
      enddo
      tt(N+1) = 1.0

      z  = 0.0
      do k=1,N
         ht= tt(k+1)-tt(k)
         t = tt(k)*w1 + tt(k+1)*w2
         z = z + exp(y2*(t*t-1))*ht
         t = tt(k)*w2 + tt(k+1)*w1
         z = z + exp(y2*(t*t-1))*ht
      enddo
      plasmaz = - y*z

      plasmaz = plasmaz + (0.,1.)*spi*exp(-y2)
      else 
         z  = 1./y
         plasmaz = -z
         loop1: do k=1,int(abs(y2))
            z=z/y2*(k-0.5)
            plasmaz=plasmaz-z
            if (abs(z).lt.1.0e-15) exit loop1
         end do loop1
c        the residual part below follows Fried and Conte, not Miyamoto
         if (imag(y).gt.0.) then
           plasmaz = plasmaz + 0.
         elseif (imag(y).lt.0.) then
           plasmaz = plasmaz + (0.,1.)*spi*2.*exp(-y2)
         else
           plasmaz = plasmaz + (0.,1.)*spi*exp(-y2)
         endif
      end if

      return
      end

      function plasmay(y)
c     should be OK for unstable mode
c     need to be extended to stable mode
      implicit none
      complex*16::plasmay,y,y2,z,z1,z2,z3
      integer::k
      real*8::w,ht,t,t2,a,b,p,pi

      integer N0,N
      parameter (N0=3602)
      real*8 tt(N0+1),t0

      pi = acos(-1.)
      y2 = y*y

      if (abs(y).lt.3.5) then

c     compute principle part
      w  = 1.0/sqrt(3.0)
      t0 = -4.0
      N = 101
      if (abs(y).gt.1.) 
     &   N = 101 + int((abs(y)-1.)**2*400)
      ht = (0.0-t0)/(N-1)
      do k=1,N
         tt(N+1-k) = 1.0 - 10**(t0+(k-1)*ht)
      enddo
      tt(N+1) = 1.0
      z1 = 0.0
      z2 = 0.0
      z3 = 0.0
      
      do k=1,N
         ht= tt(k+1)-tt(k)
         t = (tt(k)*(1-w) + tt(k+1)*(1+w))*0.5
         t2 = t*t
         z1 = z1 + exp(y2*(t2-1))*ht
         z2 = z2 + (exp(-abs(y2)*t2-y2)-exp(y2*(t2-1)))/t*ht
         z3 = z3 + exp(-abs(y2)/t2-y2)/t*ht
         t = (tt(k)*(1+w) + tt(k+1)*(1-w))*0.5
         t2 = t*t
         z1 = z1 + exp(y2*(t2-1))*ht
         z2 = z2 + (exp(-abs(y2)*t2-y2)-exp(y2*(t2-1)))/t*ht
         z3 = z3 + exp(-abs(y2)/t2-y2)/t*ht
      enddo

      z = (sqrt(pi)*y*z1 + z2 + z3)*0.5
      plasmay= y*z

c     compute residual part
      a = real(y)
      b = imag(y)
      if (a.gt.0.) then
         p = atan(b/a)
      elseif (a.lt.0.) then
         if (b.gt.0.) then
            p = atan(b/a) + pi
         else
            p = atan(b/a) - pi
         endif
      else
         if (b.gt.0.) then
            p = pi*0.5
         else
            p = -pi*0.5
         endif
      endif
      plasmay = plasmay -(0.,1.)*p*y*exp(-y2)

      else
 
      z1 = sqrt(pi)
      z2 = 1./y
      plasmay = (z1 - z2)*0.5
      loop1: do k=1,16
         z1 = z1*(k-0.5)/y2
         z2 = z2*k/y2
         plasmay = plasmay + (z1 - z2)*0.5
         if (abs(z2).lt.1.0e-14) exit loop1
      end do loop1

      endif
      
      return
      end

      function plasmax(z,a)
c     z has a factor i*gamma/a
c     two branches are defined for gamma>0, with a>0 and a<0
c     with both analytically continued to stable region gamma<0.
c     sign of imaginary part for real and negative z is determined
c     by continuous extension from unstable half-plane
      implicit none
      complex*16::plasmax,plasmaz,z,y
      real*8::a

      if (a.lt.0.) then
         y = sqrt(-z)
      else
         y = -sqrt(-z)
         if (real(z).gt.0.and.abs(imag(z)).lt.1.0e-15) y = sqrt(-z)
      endif

      plasmax = -plasmaz(y)*y/2.

      return
      end


      function plasmaz1(z)
c     this is for abs(imag(z)) less than abs(real(z))
c     I have tested the approximation method, also quite good
c     MSC 6/30/05
      implicit none
      complex*16::plasmaz1,z,zsq,papp,term,s15ddf,pd
      integer::ifail,i
      real*8::absz,pi,sqrt,abspd
c
      zsq=z*z
      absz=abs(z)
      pi=acos(-1.)
c     plasmaz1=(0.,1.)*sqrt(pi)*s15ddf(z,ifail)
      plasmaz1=0.
c
      return
c     the approximation is coded below, note it is
c     the large argument limit is not abolutely convergent
      if (absz.le.100.) then
         papp=(0.,1.)*sqrt(pi)*exp(-zsq)
      end if
      if (absz.le.4.) then
         term=2.*z
         loop1: do i=1,10000
            papp=papp-term
            term=-term*zsq*2./(2.*i+1)
            if (abs(term).le.1.e-10) exit loop1
         end do loop1
      else if (absz.gt.4.) then
         term=1./z
         loop2: do i=1,10
            papp=papp-term
            term=term/zsq*(2*i-1.)/2.
            if (abs(term).le.1.e-10) exit loop2
         end do loop2
      end if
      pd=(plasmaz1-papp)
      abspd=abs(pd)
      write (*,'("   z   ",1p2e12.3)')z
      write (*,'("plasmaz,papp,abspd",1p5e12.3)')plasmaz1,papp,abspd
      return
      end

      function plasmaf_old(cn,ca,cb,cc,rd,kf)
c     numerical integration over particle energy
c     kf = 0: integration for DRIFT=0
c          1: integration for DRIFT<>0
      implicit none
      complex*16::plasmaf_old,cn,ca,cb,cc,c0,c1,z
      real*8::rd,t,t1,w,ht
      integer::kf,kk,k

      integer N0
      parameter (N0=101)
      real*8 tt(N0+1)

c     setup
      w     = 1.0/sqrt(3.0)
      ht    = 1.0/(N0-1)
      tt(1) = 0.0
      do k=2,N0
         tt(k) = tt(k-1) + ht
      enddo

      kk = 0
      if (abs(imag(cn)).lt.1e-14.and.abs(cn).le.1.) kk = 1 

      if (kf.eq.1.and.kk.eq.1) then
         t  = abs(cn)
         t1 = t**1.5 + rd**1.5
         c0 = t**2.5/t1*(1.5*sqrt(t)/t1*ca + cc/t1 + cb)
      endif

c     integration using Gauss quadrature
      z  = 0.0
      do k=1,N0-1
         t = (tt(k)*(1-w) + tt(k+1)*(1+w))*0.5
         t1 = t**1.5 + rd**1.5
         c1 = t**2.5/t1*(1.5*sqrt(t)/t1*ca + cc/t1 + cb)
         if (kf.eq.0) z = z + c1
         if (kf.eq.1.and.kk.eq.0) z = z + c1/(t*t-cn*cn)
         if (kf.eq.1.and.kk.eq.1) z = z + (c1-c0)/(t*t-cn*cn)

         t = (tt(k)*(1+w) + tt(k+1)*(1-w))*0.5
         t1 = t**1.5 + rd**1.5
         c1 = t**2.5/t1*(1.5*sqrt(t)/t1*ca + cc/t1 + cb)
         if (kf.eq.0) z = z + c1
         if (kf.eq.1.and.kk.eq.0) z = z + c1/(t*t-cn*cn)
         if (kf.eq.1.and.kk.eq.1) z = z + (c1-c0)/(t*t-cn*cn)
      enddo

      if (kf.eq.0) z = z/cn*ht
      if (kf.eq.1) z = -z*cn*ht 
      if (kf.eq.1.and.kk.eq.1) 
     &   z = z - c0*(log((1-cn)/(1+cn))+(0.,1.)*acos(-1.))
      if (kf.eq.1.and.kk.eq.0.and.imag(cn).lt.0.0) z=conjg(z)

      plasmaf_old = z

      return
      end

      function plasmaf(cn,ca,cb,cc,rd,kf,sig)
c     numerical integration over particle energy
c     kf = 0: integration for DRIFT=0
c          1: integration for DRIFT<>0
      implicit none
      complex*16::plasmaf,cn,ca,cb,cc,c0,c1,z
      real*8::rd,t,t1,w,ht,sig
      integer::kf,kk,k

      integer N0
      parameter (N0=101)
      real*8 tt(N0+1)

c     setup
      w     = 1.0/sqrt(3.0)
      ht    = 1.0/(N0-1)
      tt(1) = 0.0
      do k=2,N0
         tt(k) = tt(k-1) + ht
      enddo

      kk = 0
      if (abs(imag(cn)).lt.1e-14.and.real(cn).le.0.and.real(cn).ge.-1)
     &    kk = 1 

      if (kf.eq.1.and.kk.eq.1) then
         t  = -cn
         t1 = t**1.5 + rd**1.5
         c0 = t**2.5/t1*(1.5*sqrt(t)/t1*ca + cc/t1 + cb)
      endif

c     integration using Gauss quadrature
      z  = 0.0
      do k=1,N0-1
         t = (tt(k)*(1-w) + tt(k+1)*(1+w))*0.5
         t1 = t**1.5 + rd**1.5
         c1 = t**2.5/t1*(1.5*sqrt(t)/t1*ca + cc/t1 + cb)
         if (kf.eq.0) z = z + c1
         if (kf.eq.1.and.kk.eq.0) z = z + c1/(t+cn)
         if (kf.eq.1.and.kk.eq.1) z = z + (c1-c0)/(t+cn)

         t = (tt(k)*(1+w) + tt(k+1)*(1-w))*0.5
         t1 = t**1.5 + rd**1.5
         c1 = t**2.5/t1*(1.5*sqrt(t)/t1*ca + cc/t1 + cb)
         if (kf.eq.0) z = z + c1
         if (kf.eq.1.and.kk.eq.0) z = z + c1/(t+cn)
         if (kf.eq.1.and.kk.eq.1) z = z + (c1-c0)/(t+cn)
      enddo

      if (kf.eq.0) z = z/cn*ht
      if (kf.eq.1) z = z*ht 
      if (kf.eq.1.and.kk.eq.1) 
     &   z = z + 2.0*c0*(log((1+cn)/(-cn)) + (0.,1.)*acos(-1.)*sig)

      plasmaf = z

      return
      end

C=======================================================================
C COMPUTE D F/D X, FOR GENERIC INPUT AND OUTPUT
C PERFORM NUMERICAL DIFFERENTIATION IN UNIFORM MESH
C KD = 0: DERIVATIVE
C      1: LOGARITHMIC DERIVATIVE
C YQL, 07-2024                                               
C=======================================================================
      SUBROUTINE DFFFDX(YOUT,XOUT,NOUT,YIN,XIN,NIN,KD)

      USE GLOBALM
      IMPLICIT NONE

      INTEGER NOUT,NIN,KD,J,KCHECK
      REAL*8  YOUT(NOUT),XOUT(NOUT),YIN(NIN),XIN(NIN)
      REAL*8  H1
      REAL*8,DIMENSION(:),ALLOCATABLE::PSI,FFI,FFW,DFFI,ZOUT

      KCHECK=0

      ALLOCATE( PSI(NIN),FFI(NIN),FFW(NIN),DFFI(NIN),ZOUT(NOUT) )

C     DEFINE UNIFORM PSI-MESH
      H1     = (XIN(NIN)-XIN(1))/DFLOAT(NIN-1)
      PSI(1) = XIN(1)
      DO J=2,NIN
         PSI(J) = PSI(J-1) + H1
      ENDDO

C     SPLINE YIN ONTO PSI-MESH
      CALL SPLINE1D(FFI,PSI,NIN,YIN,XIN,NIN,FFW)

C     COMPUTE DFFI=D FFI/D PSI AT UNIFORM HALF-INTEGER MESH
      DFFI(2:NIN) = (FFI(2:NIN)-FFI(1:NIN-1))/H1
      
C     EXTEND DFFI TO FULL XIN-RANGE, STORE NEW PSI-MESH INTO FFI
      FFI(1)       = PSI(1)
      FFI(2:NIN-1) = (PSI(1:NIN-2)+PSI(2:NIN-1))*0.5
      FFI(NIN)     = PSI(NIN)
      DFFI(1)      = 1.5*DFFI(2)-0.5*DFFI(3)
      DFFI(NIN)    = 1.5*DFFI(NIN)-0.5*DFFI(NIN-1)

C     SPLINE DFFI BACK ONTO XOUT-MESH
      CALL SPLINE1D(YOUT,XOUT,NOUT,DFFI,FFI,NIN,FFW)

      IF (KD.EQ.1) THEN
         CALL SPLINE1D(ZOUT,XOUT,NOUT,YIN,XIN,NIN,FFW)
         YOUT = YOUT/ZOUT
      ENDIF

      IF (NOUT.GE.3) THEN
         YOUT(1)    = YOUT(2)
         YOUT(NOUT) = YOUT(NOUT-1)

C        SMOOTHING THE ABOVE DERIVATIVE
         H1 = 0.1
         DO J=1,NKSMOOTHR
            YOUT(2:NOUT-1) = H1*YOUT(1:NOUT-2)+(1.-2.*H1)*YOUT(2:NOUT-1)
     &                       +H1*YOUT(3:NOUT)
         ENDDO
      ENDIF

      DEALLOCATE( PSI,FFI,FFW,DFFI,ZOUT )

      RETURN
      END

C=======================================================================
C COMPUTE D F/D PSI, FOR F(JS,KGRID)
C NOTE WE PERFORM NUMERICAL DIFFERENTIATION IN UNIFORM PSI=CS^2 MESH
C YQL, 08-2013                                               
C=======================================================================
      SUBROUTINE DFFFDPSI(KD)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      IMPLICIT NONE

      INTEGER KD,JS,KCHECK
      REAL*8  H1
      REAL*8,DIMENSION(:),ALLOCATABLE::PSI,FFI,FFW,DFFI

      KCHECK=0

      ALLOCATE( PSI(NRP1),FFI(NRP1),FFW(NRP1),DFFI(NRP1) )

C     DEFINE UNIFORM PSI-MESH
      H1     = (CS(NRP1)**2-CS(1)**2)/DFLOAT(NR)
      PSI(1) = CS(1)**2
      DO JS=2,NRP1
         PSI(JS) = PSI(JS-1) + H1
      ENDDO

C     SPLINE FFF(:,1) ONTO PSI-MESH
      DFFI = CS(1:NRP1)**2
      CALL SPLINE1D(FFI,PSI,NRP1,FFF(:,1),DFFI,NRP1,FFW)

C     COMPUTE DFFI=D FFI/D PSI AT UNIFORM HALF-INTEGER MESH
      DFFI(2:NRP1) = (FFI(2:NRP1)-FFI(1:NR))/H1
      
C     EXTEND DFFI TO FULL RADIUS, STORE NEW PSI-MESH INTO FFI
      FFI(1)    = PSI(1)
      FFI(2:NR) = (PSI(1:NR-1)+PSI(2:NR))*0.5
      FFI(NRP1) = PSI(NRP1)
      DFFI(1)   = 1.5*DFFI(2)-0.5*DFFI(3)
      DFFI(NRP1)= 1.5*DFFI(NRP1)-0.5*DFFI(NR)

C     COMPUTE DFFI=D FFI/D S IN NEW PSI-MESH
      DFFI = DFFI*2.*SQRT(FFI)

C     SPLINE DFFI BACK ONTO ORIGINAL CS-MESH
      CALL SPLINE1D(DFFF(:,1),CS(1:NRP1),NRP1,DFFI,SQRT(FFI),NRP1,FFW)
      CALL SPLINE1D(DFFF(1:NR,2),CSM(1:NR),NR,DFFI,SQRT(FFI),NRP1,FFW)

C     FINAL OPERATIONS
      DFFF(2:NR,1) = DFFF(2:NR,1)/DPSIDS (2:NR)
      DFFF(1:NR,2) = DFFF(1:NR,2)/DPSIDSM(1:NR)
    
      IF (KD.EQ.1) THEN
         DFFF(2:NR,1) = DFFF(2:NR,1)/FFF(2:NR,1)
         DFFF(1:NR,2) = DFFF(1:NR,2)/FFF(1:NR,2)
      ENDIF

      DFFF(1,1)    = 2.*DFFF(1,2)  - DFFF(2,1)
      DFFF(NRP1,1) = 2.*DFFF(NR,2) - DFFF(NR,1)

C     SMOOTHING THE ABOVE DERIVATIVE
      H1 = 0.1
      DO JS=1,NKSMOOTHR
         DFFF(2:NR,1) = H1*DFFF(1:NR-1,1)+(1.-2.*H1)*DFFF(2:NR,1)+
     &                  H1*DFFF(3:NRP1,1)
         DFFF(2:NR-1,2) = H1*DFFF(1:NR-2,2)+(1.-2.*H1)*DFFF(2:NR-1,2)+
     &                    H1*DFFF(3:NR,2) 
      ENDDO

      DEALLOCATE( PSI,FFI,FFW,DFFI )

      RETURN
      END

C=======================================================================
C COMPUTE D LN(F)/D PSI, FOR F(JS,KGRID)
C YQL, 07-2013                                               
C=======================================================================
      SUBROUTINE DFFFDPSI_OLD(KD)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      IMPLICIT NONE

      INTEGER KD,K,JS,KCHECK
      REAL*8  H1,H2

      KCHECK=0

      DO JS=2,NR
         H1 = (CS(JS)-CS(JS-1))/2.
         H2 = (CS(JS+1)-CS(JS))/2.
         DFFF(JS,1) = ( (H1/H2*FFF(JS,2)-H2/H1*FFF(JS-1,2))/(H1+H2)-
     &                  (H1-H2)*FFF(JS,1)/H1/H2 )/DPSIDS(JS)
      ENDDO

      DO JS=1,NR
         H1 = CS(JS+1)-CS(JS)
         DFFF(JS,2) = ( FFF(JS+1,1)-FFF(JS,1) )/H1/DPSIDSM(JS)
      ENDDO
      
      IF (KD.EQ.1) THEN
         DFFF(2:NR,1) = DFFF(2:NR,1)/FFF(2:NR,1)
         DFFF(1:NR,2) = DFFF(1:NR,2)/FFF(1:NR,2)
      ENDIF

      DFFF(1,1)    = 2.*DFFF(1,2)  - DFFF(2,1)
      DFFF(NRP1,1) = 2.*DFFF(NR,2) - DFFF(NR,1)

      RETURN
      END

      SUBROUTINE ALLOCATEFORPARALLEL
      USE KINETICM
      USE ANISOTROPICM
      USE DIMENSIM
      USE GLOBALM
      USE OMP_LIB
      IMPLICIT NONE

      IF (OMP_GET_THREAD_NUM() == 0) RETURN

      ALLOCATE( KNUMDISTRIB(NSPECIES) ) 
      ALLOCATE( LAMM(2*NLAMK+2), LAMHH(2*NLAMK),
     &          LAMTMP(2*NLAMK+2) )      

      ALLOCATE( VPARA0(MSMAX,MLMAX), VPERP0(MSMAX,MLMAX),
     &          VDPHI0(MSMAX,MLMAX),
     &          VX10(MSMAX,MLMAX),   VX20(MSMAX,MLMAX), 
     &          VQ10(MSMAX,MLMAX),   VQ20(MSMAX,MLMAX), 
     &          VQ30(MSMAX,MLMAX),   VDP0(MSMAX,MLMAX) )

      ALLOCATE( SVPARA0(MSMAX,MLMAX,NSPECIES),
     &          SVPERP0(MSMAX,MLMAX,NSPECIES), 
     &          SVDPHI0(MSMAX,MLMAX,NSPECIES), 
     &          SVX10  (MSMAX,MLMAX,NSPECIES),
     &          SVX20  (MSMAX,MLMAX,NSPECIES), 
     &          SVQ10  (MSMAX,MLMAX,NSPECIES),
     &          SVQ20  (MSMAX,MLMAX,NSPECIES),
     &          SVQ30  (MSMAX,MLMAX,NSPECIES),
     &          SVDP0  (MSMAX,MLMAX,NSPECIES),
     &          SLAM0  (MLMAX,NSPECIES),       
     &          SF0    (MLMAX,NSPECIES,0:3,4) )

      ALLOCATE( RCHI2(NCHI+1), RW1(NCHI+1),
     &          RW2(NCHI+1),
     &          RJB(NCHI+1),  RX1P(NCHI+1),  RX1B(NCHI+1), 
     &          RX1R(NCHI+1), RX2(NCHI+1),   RQ1(NCHI+1),   
     &          RQ2(NCHI+1),  RBT(NCHI+1),   RPHI(NCHI+1),  
     &          RDMU(NCHI+1), RDB(NCHI+1) )
     
      ALLOCATE( RCHIK(NCHIT+2), RPHIK(NCHIT+2), RTK(NCHIT+2),
     &          RHK(NCHIT+2),   RJBK(NCHIT+2),  RX1PK(NCHIT+2),
     &          RX1BK(NCHIT+2), RX1RK(NCHIT+2), RX2K(NCHIT+2),  
     &          RQ1K(NCHIT+2),  RQ2K(NCHIT+2),  RVALK(NCHIT+2) )
      
      ALLOCATE( RCHIN(2*NCHIT+2), RVALN(2*NCHIT+2), RHN(2*NCHIT+2) )

      ALLOCATE ( VI(4,MLMAX,NSPECIES) )      

      ALLOCATE( VPARA(MSMAX,MLMAX), VPERP(MSMAX,MLMAX),
     &          VDPHI(MSMAX,MLMAX),
     &          VX1(MSMAX,MLMAX),   VX2(MSMAX,MLMAX), 
     &          VQ1(MSMAX,MLMAX),   VQ2(MSMAX,MLMAX), 
     &          VQ3(MSMAX,MLMAX),   VDP(MSMAX,MLMAX) )

      ALLOCATE( RVAK1(NCHIT+2), RVAK2(NCHIT+2), RVAK3(NCHIT+2), 
     &          RVAK4(NCHIT+2), RVAK5(NCHIT+2), RVAK6(NCHIT+2),
     &          RVAK7(NCHIT+2) )
      
      CALL ALLOCATEDWKCOMPMAT
      
C     IF (IFOWP.EQ.1.OR.IFOWT.EQ.1) THEN
         ALLOCATE(VI1(4,MLMAX,NSPECIES), 
     &            VI2(4,MLMAX,NSPECIES), 
     &            VI3(4,MLMAX,NSPECIES)) 
         ALLOCATE(VPARA1(MSMAX,MLMAX), VPERP1(MSMAX,MLMAX),
     &            VDPHI1(MSMAX,MLMAX), 
     &            VX11(MSMAX,MLMAX),   VX21(MSMAX,MLMAX), 
     &            VQ11(MSMAX,MLMAX),   VQ21(MSMAX,MLMAX), 
     &            VQ31(MSMAX,MLMAX),   VDP1(MSMAX,MLMAX) )
         ALLOCATE(VPARA01(MSMAX,MLMAX), VPERP01(MSMAX,MLMAX),
     &            VDPHI01(MSMAX,MLMAX),
     &            VX101(MSMAX,MLMAX),   VX201(MSMAX,MLMAX), 
     &            VQ101(MSMAX,MLMAX),   VQ201(MSMAX,MLMAX), 
     &            VQ301(MSMAX,MLMAX),   VDP01(MSMAX,MLMAX) )
         ALLOCATE(SVPARA01(MSMAX,MLMAX,NSPECIES),
     &            SVPERP01(MSMAX,MLMAX,NSPECIES), 
     &            SVDPHI01(MSMAX,MLMAX,NSPECIES), 
     &            SVX101  (MSMAX,MLMAX,NSPECIES),
     &            SVX201  (MSMAX,MLMAX,NSPECIES), 
     &            SVQ101  (MSMAX,MLMAX,NSPECIES),
     &            SVQ201  (MSMAX,MLMAX,NSPECIES),
     &            SVQ301  (MSMAX,MLMAX,NSPECIES),
     &            SVDP01  (MSMAX,MLMAX,NSPECIES) )
         ALLOCATE(VX1LNP(MSMAX), VX2LNP(MSMAX), VDPLNP(MSMAX), 
     &            VQ1LNP(MSMAX), VQ2LNP(MSMAX), VQ3LNP(MSMAX))
C     ENDIF

      CALL ZALLOCANISO(0)

      END SUBROUTINE ALLOCATEFORPARALLEL
      
      SUBROUTINE DEALLOCATEFORPARALLEL
      USE KINETICM
      USE GLOBALM
      USE OMP_LIB
      
      IMPLICIT NONE
      IF (OMP_GET_THREAD_NUM() == 0) RETURN
      
      DEALLOCATE( KNUMDISTRIB ) 
      DEALLOCATE( LAMM,LAMHH,LAMTMP )
      DEALLOCATE( VPARA0,VPERP0,VDPHI0,VX10,VX20,VQ10,VQ20,VQ30,VDP0 )
      DEALLOCATE( SVPARA0,SVPERP0,SVDPHI0,SVX10,SVX20,SVQ10,SVQ20,
     &            SVQ30,SVDP0 )
      DEALLOCATE( SLAM0,SF0 )
      DEALLOCATE( RCHI2,RW1,RW2,
     &            RJB,RX1P,RX1B,RX1R,RX2,RQ1,RQ2,
     &            RBT,RPHI,RDMU,RDB )
      DEALLOCATE( RCHIK,RPHIK,RTK,
     &            RHK,RJBK,RX1PK,RX1BK,RX1RK,
     &            RX2K,RQ1K,RQ2K,RVALK )
      DEALLOCATE( RCHIN,RVALN,RHN )
      DEALLOCATE( VI )
      DEALLOCATE( VPARA,VPERP,VDPHI,VX1,VX2,VQ1,VQ2,VQ3,VDP )
      DEALLOCATE( RVAK1,RVAK2,RVAK3,RVAK4,RVAK5,RVAK6,RVAK7 )

C     IF (IFOWP.EQ.1.OR.IFOWT.EQ.1) THEN
         DEALLOCATE( VI1,VI2,VI3 ) 
         DEALLOCATE( VPARA1,VPERP1,VDPHI1,VX11,VX21,VQ11,VQ21,VQ31,VDP1)
         DEALLOCATE( VPARA01,VPERP01,VDPHI01,VX101,VX201, 
     &               VQ101,VQ201,VQ301,VDP01 )
         DEALLOCATE( SVPARA01,SVPERP01,SVDPHI01,SVX101,SVX201,
     &               SVQ101,SVQ201,SVQ301,SVDP01 )
         DEALLOCATE( VX1LNP,VX2LNP,VQ1LNP,VQ2LNP,VQ3LNP,VDPLNP ) 
C     ENDIF

      CALL ZDEALLOCANISO(0)

      CALL DEALLOCATEDWKCOMPMAT
      
      END SUBROUTINE DEALLOCATEFORPARALLEL
      
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE RADIAL DENSITY DWK^(A)                       LIU YQ 04.01.2013
C AT HALF-INTEGER RADIAL GRID                                          $
C SUCH THAT TOTAL DWK=INT DWK^(A) DS                                   $
C THIS IS JUST AN ALTERNATIVE WAY FOR COMPUTING THE DRIFT KINETIC ENERGY
C DENSITY IN THE QUADRATIC FORM. IT HAS BEEN TESTED (ON ITER
C EQUILIBRIUM) THAT THIS SUBROUTINE FULLY RECOVERS THAT OF ENERGYMAT WITH
C KEFORM=2.
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE KDWKDENSITY  
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      USE TORQUEM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INTEGER    I,J,MS,KCHECK
      REAL*8     HCHI
      COMPLEX*16 CTMP1,ODWK
      REAL*8,DIMENSION(:,:),ALLOCATABLE::B_2,B_2M,B2C
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::OB1,OB2,OB3,
     &                                       OX1,OX2,OPE,OPA,
     &                                       OFW
      COMPLEX*16,DIMENSION(:),ALLOCATABLE::  DWKA

      KCHECK = 0

C     STEP 1: COMPUTE (B1,B2,B3,X1,X2,PPERP,PPARA) IN (R,Z)-GRID
      ALLOCATE( OB1(NR,NCHI), OB2(NR,NCHI), OB3(NR,NCHI),
     &          OX1(NR,NCHI), OX2(NR,NCHI), OPE(NR,NCHI), OPA(NR,NCHI) )

      OB1 = (0.,0.)
      OB2 = (0.,0.)
      OB3 = (0.,0.)
      OX1 = (0.,0.)
      OX2 = (0.,0.)
      OPE = (0.,0.)
      OPA = (0.,0.)

      HCHI = 2.*PI/NCHI
      DO J=1,NCHI
         DO MS=1,MSMAX
            CTMP1 = EXP(CI*RM(MS,2)*(J-1)*HCHI)
            DO I=1,NR
               OB1(I,J) = OB1(I,J) + (B1U(I,MS)+B1U(I+1,MS))*.5*CTMP1
               OB2(I,J) = OB2(I,J) + B2U(I,MS)*CTMP1
               OB3(I,J) = OB3(I,J) + B3U(I,MS)*CTMP1
               OX1(I,J) = OX1(I,J) + (X1U(I,MS)+X1U(I+1,MS))*.5*CTMP1
               OX2(I,J) = OX2(I,J) + X2U(I,MS)*CTMP1
               OPE(I,J) = OPE(I,J) + PPERP(I,MS)*CTMP1
               OPA(I,J) = OPA(I,J) + PPARA(I,MS)*CTMP1
            ENDDO
         ENDDO
      ENDDO

C     STEP 2: COMPUTE B^2 IN BOTH INTEGER AND HALF-INTEGER RADIAL GRIDS
      ALLOCATE( B_2(NRP1,NCHI), B_2M(NR,NCHI), B2C(NR,NCHI) )

      DO J=1,NCHI
         DO I=2,NRP1
            B_2(I,J)=G22L(I,J)*DPSIDS(I)**2/RJA(I,J)**2 +
     &               T(I)**2/REQ(I,J)**2
         ENDDO 
         B_2(1,J)=T(1)**2/REQ(1,J)**2
         DO I=1,NR
            B_2M(I,J)=G22LM(I,J)*DPSIDSM(I)**2/RJAM(I,J)**2 +
     &             TM(I)**2/REQM(I,J)**2
         ENDDO 
      ENDDO 

C     DB^2/DCHI
      CALL DERCHI(B_2M,B2C,NR,NR)

C     STEP 3: COMPUTE FULL INTEGRAND OFW FOR DWK^(A)
      ALLOCATE( OFW(NR,NCHI) )

      DO J=1,NCHI
      DO I=1,NR
         OFW(I,J) = DPSIDSM(I)*G12LM(I,J)/RJAM(I,J)/B_2M(I,J)*
     &                 CONJG(OB1(I,J))*OPE(I,J) + 
     &              DPSIDSM(I)*G22LM(I,J)/RJAM(I,J)/B_2M(I,J)*
     &                 CONJG(OB2(I,J))*OPE(I,J) + 
     &              TM(I)/B_2M(I,J)*CONJG(OB3(I,J))*OPE(I,J) +
     &              RJAM(I,J)/B_2M(I,J)*PPEQM(I)*DPSIDSM(I)*
     &                 CONJG(OX1(I,J))*OPA(I,J) + 
     &             (RJAM(I,J)/2./B_2M(I,J)*(B_2(I+1,J)-B_2(I,J))/CSH(I)-
     &              DPSIDSM(I)**2*G12LM(I,J)/2./RJAM(I,J)/B_2M(I,J)**2*
     &              B2C(I,J))*CONJG(OX1(I,J))*(OPE(I,J)+OPA(I,J)) +
     &              RJAM(I,J)*TM(I)/2./B_2M(I,J)**2*B2C(I,J)*
     &              CONJG(OX2(I,J))*(OPE(I,J)+OPA(I,J))
      ENDDO
      ENDDO

      DEALLOCATE( OB1,OB2,OB3,OX1,OX2,OPE,OPA,B_2,B_2M,B2C )

C     STEP 4: COMPUTE DWK^(A)
      ALLOCATE( DWKA(NR) )

      DWKA = (0.,0.)

      DO J=1,NCHI
      DO I=1,NR
         DWKA(I) = DWKA(I) + OFW(I,J)
      ENDDO
      ENDDO
      DWKA = DWKA*PI*HCHI
   
      DWKA(1) = 0.0
      DO I=1,NR
         IF (CSM(I).GT.CTEDGE) DWKA(I) = 0.0
      ENDDO

      IF (KNTV.EQ.20) THEN
C     COMPUTE NTV TORQUE DENSITY FROM DWK
C     USING: T_NTV = -2*N*IM(DWKA)/(4*PI^2)
      TORQUENTV = -2.*RNTOR*IMAG(DWKA)/(4.*PI*PI)
 
      ELSE
      
C     STEP 5: COMPUTE TOTAL DWK
      ODWK = (0.,0.)
      DO I=1,NR
         ODWK = ODWK + DWKA(I)*CSH(I)
      ENDDO

C     STEP 6: SAVE ENERGY DENSITY TO A FILE
      OPEN(CHOUTP,FILE='ENERGY_DENSITY.OUT')
      REWIND(CHOUTP)
      DO I=1,NR
         WRITE(CHOUTP,120) CSM(I),
     &                     REAL(DWKA(I)/(4.*PI*PI)),    
     &                     IMAG(DWKA(I)/(4.*PI*PI))    
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(5(E15.8,1X))

      WRITE(*,*) 'TOTAL KINETIC ENERGY = ',ODWK
 
      ENDIF
      
      DEALLOCATE( DWKA,OFW )

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE PERTURBED MAGNETIC ENERGY DENSITY -0.5*INT[|Q|^2]DV          $
C IN THE PLASMA REGION                                                 $
C TO BE CALLED FROM ENERGYMAT WITH KEFORM=2                            $
C YQ LIU, 2013-09
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE KDWFMAGP(DWKA)
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      IMPLICIT NONE
      INTEGER    I,J,MS,KCHECK
      REAL*8     HCHI
      COMPLEX*16 CTMP1,DWKA(*)
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::OB1,OB2,OB3
      REAL*8    ,DIMENSION(:,:),ALLOCATABLE::OFW

      KCHECK = 0

C     STEP 1: COMPUTE (B1,B2,B3) IN (R,Z)-GRID
      ALLOCATE( OB1(NR,NCHI), OB2(NR,NCHI), OB3(NR,NCHI) )

      OB1 = (0.,0.)
      OB2 = (0.,0.)
      OB3 = (0.,0.)

      HCHI = 2.*PI/NCHI
      DO J=1,NCHI
         DO MS=1,MSMAX
            CTMP1 = EXP(CI*RM(MS,2)*(J-1)*HCHI)
            DO I=1,NR
               OB1(I,J) = OB1(I,J) + (B1U(I,MS)+B1U(I+1,MS))*.5*CTMP1
               OB2(I,J) = OB2(I,J) + B2U(I,MS)*CTMP1
               OB3(I,J) = OB3(I,J) + B3U(I,MS)*CTMP1
            ENDDO
         ENDDO
      ENDDO

C     STEP 2: COMPUTE FULL INTEGRAND OFW FOR DWK^(A)
      ALLOCATE( OFW(NR,NCHI) )

      DO J=1,NCHI
      DO I=1,NR
         OFW(I,J) = (ABS(OB1(I,J))**2*G11LM(I,J)+
     &               ABS(OB2(I,J))**2*G22LM(I,J)+
     &               ABS(OB3(I,J))**2*G33LM(I,J)+
     &               2.*DREAL(OB1(I,J)*CONJG(OB2(I,J)))*G12LM(I,J))/
     &               RJAM(I,J)
      ENDDO
      ENDDO

      DEALLOCATE( OB1,OB2,OB3 )

C     STEP 3: COMPUTE DWFMAG DENSITY 
      DWKA(1:NR) = (0.,0.)
      DO I=1,NR
      DO J=1,NCHI
         DWKA(I) = DWKA(I) - OFW(I,J)
      ENDDO
      DWKA(I) = DWKA(I)*CSH(I)*PI*HCHI
      ENDDO
   
      DEALLOCATE(OFW)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C COMPUTE PERTURBED MAGNETIC ENERGY 0.5*INT[|Q|^2]DV                   $
C IN THE VACUUM REGION                                                 $
C TO BE CALLED FROM ENERGYMAT WITH KEFORM=2                            $
C YQ LIU, 2013-09
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE KDWFMAGV(DWV)
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE GVACUUMM
      USE RCOMDM  
      IMPLICIT NONE
      INTEGER    I,J,L,MS,KCHECK
      REAL*8     HCHI
      COMPLEX*16 CTMP1,DWV
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE::OB1,OB2,OB3
      REAL*8    ,DIMENSION(:,:),ALLOCATABLE::OFW
      REAL*8    ,DIMENSION(:)  ,ALLOCATABLE::DWKA

      KCHECK = 0

C     STEP 1: COMPUTE (B1,B2,B3) IN (R,Z)-GRID
      ALLOCATE( OB1(NV,NCHI), OB2(NV,NCHI), OB3(NV,NCHI) )

      OB1 = (0.,0.)
      OB2 = (0.,0.)
      OB3 = (0.,0.)

      HCHI = 2.*PI/NCHI
      DO J=1,NCHI
         DO MS=1,MSMAX
            CTMP1 = EXP(CI*RM(MS,2)*(J-1)*HCHI)
            DO I=1,NV
               L=I+NR
               OB1(I,J) = OB1(I,J) + (B1U(L,MS)+B1U(L+1,MS))*.5*CTMP1
               OB2(I,J) = OB2(I,J) + B2U(L,MS)*CTMP1
               OB3(I,J) = OB3(I,J) + B3U(L,MS)*CTMP1
            ENDDO
         ENDDO
      ENDDO

C     STEP 2: COMPUTE FULL INTEGRAND OFW FOR DWK^(A)
      ALLOCATE( OFW(NV,NCHI) )

      DO J=1,NCHI
      DO I=1,NV
         OFW(I,J) = (ABS(OB1(I,J))**2*VRG11LM(I,J)+
     &               ABS(OB2(I,J))**2*VRG22LM(I,J)+
     &               ABS(OB3(I,J))**2*VRG33LM(I,J)+
     &               2.*DREAL(OB1(I,J)*CONJG(OB2(I,J)))*VRG12LM(I,J))/
     &               VRJAM(I,J)
      ENDDO
      ENDDO

      DEALLOCATE( OB1,OB2,OB3 )

C     STEP 3: COMPUTE DWFMAG DENSITY 
      ALLOCATE( DWKA(NV) )
      DWKA = 0.
      DO I=1,NV
      DO J=1,NCHI
         DWKA(I) = DWKA(I) + OFW(I,J)
      ENDDO
      DWKA(I) = DWKA(I)*VCSH(I)*PI*HCHI
      ENDDO
   
C     STEP 4: COMPUTE TOTAL VACUUM MAGNETIC ENERGY
      DWV = SUM(DWKA)

      DEALLOCATE(OFW,DWKA)

      RETURN
      END

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C THERMAL ION FULL DRIFT TORBIT INTEGRATION NEAR ONE SURFACE           $
C GIVEN BY JSOUT                                                       $
C USING 4TH ORDER RUNGE-KUTTA METHOD AND 2D SPLINE                     $
C LIU YQ 11.06.2013                                                    $
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE KORBIT       
C     ==========================================================
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM  
      USE KINETICM
      IMPLICIT NONE
      INCLUDE 'comioc.inc'
      INTEGER    I,J,K,JS,KSTEP
      REAL*8     ORBS0,ORBC0,ORBP0,ORBS1,ORBC1,ORBP1,
     &           ORBSA,ORBCA,ORBSIG,ORBHK1,ORBTB,ORBWD,ORBWTH,
     &           ORBS2,ORBC2,ORBP2,ORBSK1,ORBCK1,ORBPK1,
     &           ORBSK2,ORBCK2,ORBPK2,ORBSK3,ORBCK3,ORBPK3,
     &           ORBSK4,ORBCK4,ORBPK4,ORBT,ORBDT,ORBR1,ORBZ1,
     &           H1,H2,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMPCP,TMPCB
      REAL*8,DIMENSION(:,:),ALLOCATABLE::ORBES,ORBECP,ORBECM,
     &                                   ORBEPP,ORBEPM,ORBEH,
     &                                   ORBER,ORBEZ,
     &                                   ORBFS,ORBFCP,ORBFCM,
     &                                   ORBFPP,ORBFPM,ORBFH,
     &                                   ORBFR,ORBFZ
      REAL*8,DIMENSION(:,:),ALLOCATABLE::B_2,B_2M,B2S,B2C

C     INPUT PARAMETERS FOR THE THERMAL ION ENERGY AND PITCH ANGLE
C     AS WELL AS TIME STEPPING PARAMETERS
C     ORBE IS NORMALISED (BY T) PARTICLE ENERGY
C     ORBL IS NORMALISED PITCH ANGLE, SUCH THAT 
C     ORBL=[0,0.5] CORRESPONDS TO PASSING PARTICLE, AND 
C     ORBL=[0.5,1] CORRESPONDS TO TRAPPED PARTICLE
C     ORBSIGINI SPECIFIES INITIAL DIRECTION OF PARTICLES
C     ORBE = 1.0
C     ORBL = 0.25
C     NORB = 1000
C     ORBSIGINI = 1.0

      JS   = 10
      IF (JSOUT.GT.0) JS   = JSOUT

C     IFOW = 1: WITH FULL FOW INCLUDED DURING PARTICLE ORBIT INTEGRATION
C            0: WITH VANISHING FOW
      
      IF (ORBL.LT.0.5) THEN
         LAM = 2.*HKMIN(JS,1)*ORBL
      ELSE
         LAM = 2.*(1.-ORBL)*HKMIN(JS,1) + (2.*ORBL-1.)*HKMAX(JS,1)
      ENDIF

C     COMPUTE BOUNCE/TRANSIT TIME AND PRECESSION DRIFT FREQUENCY (TRAPPED)
C     AT VANISHING FOW, FOR BENCHMARKING PURPOSE
C     ALSO CALCULATE ANALYTICALLY BANANA ORBIT WIDTH 
      CALL KEQUIL(JS,1)
      IF (ORBL.LE.0.5) THEN
         CHIL = 0.0
         CHIU = 2.0*PI
         CALL KCHI(1)
         CALL KEQUILK(JS,1)
         CALL KBTIME(JS,1,1)
         OMEGAB = 2.*PI/RTK(NCHI2+2)
         CALL KDRIFT(JS,1,1)
         ORBTB = RTK(NCHI2+2)/SQRT(2.*ESPECIES_PRE(JS,1,1)/RHO(JS)*ORBE)
         ORBWD = DRIFT*ESPECIES_PRE(JS,1,1)/RHO(JS)*B0K/OMEGACI0*ORBE
         ORBWTH= 0.
      ELSE
         CALL KTURN(JS,1)
         CALL KCHI(0)
         CALL KEQUILK(JS,1)
         CALL KBTIME(JS,1,0)
         OMEGAB = PI/RTK(NCHI2+2)
         CALL KDRIFT(JS,1,0)
         ORBTB = 2.*RTK(NCHI2+2)/SQRT(2.*ESPECIES_PRE(JS,1,1)/RHO(JS)*
     &           ORBE)
         ORBWD = DRIFT*ESPECIES_PRE(JS,1,1)/RHO(JS)*B0K/OMEGACI0*ORBE
         ORBWTH=2.*SQRT(2.*ESPECIES_PRE(JS,1,1)/RHO(JS)*ORBE)/OMEGACI0*
     &          T(JS)/DPSIDS(JS)*SQRT(HKMAX(JS,1)*(HKMAX(JS,1)-LAM))
      ENDIF

      ORBDT= ORBTB/NORB*2. 

C     INITIAL POSITION OF PARTICLE AT T=0
      ORBS0 = CS(JS)
      ORBC0 = CHIL
      ORBP0 = 0.
      IF (ORBL.LT.0.5) ORBC0 = PI
      
C     PREPARE EQUILIBRIUM COEFFICIENTS FOR 2D SPLINES
C     IN TOTAL 5 SUCH 2D COEFFICIENTS FOR TRAPPED PARTICLE
      ALLOCATE(ORBES (NRP1,NCHI+1), ORBEH (NRP1,NCHI+1),
     &         ORBECP(NRP1,NCHI+1), ORBECM(NRP1,NCHI+1),
     &         ORBEPP(NRP1,NCHI+1), ORBEPM(NRP1,NCHI+1),
     &         ORBER (NRP1,NCHI+1), ORBEZ (NRP1,NCHI+1),
     &         ORBFS (NRP1,NCHI+1), ORBFH (NRP1,NCHI+1),
     &         ORBFCP(NRP1,NCHI+1), ORBFCM(NRP1,NCHI+1),
     &         ORBFPP(NRP1,NCHI+1), ORBFPM(NRP1,NCHI+1),
     &         ORBFR (NRP1,NCHI+1), ORBFZ (NRP1,NCHI+1))

      ALLOCATE( B_2(NRP1,NCHI),  B_2M(NRP1,NCHI), 
     &          B2S(NRP1,NCHI),  B2C(NRP1,NCHI) )

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

C     DB^2/DS
      DO J=1,NCHI
         FFF(:,1) = B_2(:,J)
         FFF(:,2) = B_2M(:,J)
         CALL DFFFDPSI(0)
         B2S(:,J) = DFFF(:,1)
      ENDDO

      DO J=1,NCHI
      DO I=2,NRP1
         TMP1 = ESPECIES_PRE(I,1,1)/RHO(I)*B0K/OMEGACI0*ORBE
         TMP2 = SQRT(2.*ESPECIES_PRE(I,1,1)/RHO(I)*ORBE)
         TMP3  = 1.-LAM/HK(I,J,1) 
         TMP4  = T(I)/RJA(I,J)/B_2(I,J)
         TMPCP = 2.*TMP3/B_2(I,J)
         TMPCB = (1.-LAM/HK(I,J,1)/2.)/B_2(I,J)

         ORBES(I,J)  = -TMP1*TMP4*TMPCB*B2C(I,J)*IFOW
         TMP5        = TMP1*TMP4*(TMPCP*PPEQ(I)*DPSIDS(I)+
     &                 TMPCB*B2S(I,J))  
C    &                 -TMP4*DPSIDS(I)*OMEGAE0(I,1)
         TMP6        = TMP2*SQRT(ABS(TMP3))*DPSIDS(I)/RJA(I,J)*
     &                 HK(I,J,1)/B0K
         ORBECP(I,J) = TMP5*IFOW + TMP6*ORBSIGINI
         ORBECM(I,J) = TMP5*IFOW - TMP6*ORBSIGINI
         TMP5        = TMP1*DPSIDS(I)/RJA(I,J)**2/B_2(I,J)*
     &                 (G12L(I,J)*TMPCB*B2C(I,J) - 
     &                  G22L(I,J)*(TMPCP*PPEQ(I)*DPSIDS(I)+
     &                             TMPCB*B2S(I,J)))  
C    &                 +G22L(I,J)*(DPSIDS(I)/RJA(I,J))**2/
C    &                 B_2(I,J)*OMEGAE0(I,1)
         TMP6        = TMP2*SQRT(ABS(TMP3))*T(I)*HK(I,J,1)/
     &                 B0K/REQ(I,J)**2
         ORBEPP(I,J) = TMP5*IFOW + TMP6*ORBSIGINI
         ORBEPM(I,J) = TMP5*IFOW - TMP6*ORBSIGINI
      ENDDO

      ORBES(1,J)  = ORBES(2,J)
      ORBECP(1,J) = ORBECP(2,J)
      ORBECM(1,J) = ORBECM(2,J)
      ORBEPP(1,J) = ORBEPP(2,J)
      ORBEPM(1,J) = ORBEPM(2,J)
      ENDDO

      DO J=1,NCHI
      DO I=1,NRP1
         ORBEH(I,J) = HK(I,J,1)
         ORBER(I,J) = REQ(I,J)
         ORBEZ(I,J) = ZEQ(I,J)
      ENDDO
      ENDDO
 
      DO I=1,NRP1
         ORBES(I,NCHI+1)  = ORBES(I,1)
         ORBEH(I,NCHI+1)  = ORBEH(I,1)
         ORBER(I,NCHI+1)  = ORBER(I,1)
         ORBEZ(I,NCHI+1)  = ORBEZ(I,1)
         ORBECP(I,NCHI+1) = ORBECP(I,1)
         ORBECM(I,NCHI+1) = ORBECM(I,1)
         ORBEPP(I,NCHI+1) = ORBEPP(I,1)
         ORBEPM(I,NCHI+1) = ORBEPM(I,1)
      ENDDO

C     STORE 2D SPLINE COEFFICIENTS 
      CALL SPLINE2D(ORBES,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFS)
      CALL SPLINE2D(ORBEH,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFH)
      CALL SPLINE2D(ORBER,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFR)
      CALL SPLINE2D(ORBEZ,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFZ)
      CALL SPLINE2D(ORBECP,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFCP)
      CALL SPLINE2D(ORBECM,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFCM)
      CALL SPLINE2D(ORBEPP,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFPP)
      CALL SPLINE2D(ORBEPM,CS,RCHI,NRP1,NCHI+1,NRP1,ORBFPM)

C     RK4 METHOD FOR TIME STEPPING OF PARTICLE TRAJECTORY
C     FOR SIMPLICITY, ASSUMING UNIFORM TIME STEP
C     CAN BE EXTENDED TO ADAPTIVE TIME STEPPING IN THE FUTURE
      OPEN(CHOUTP,FILE='ORBIT.OUT')
      REWIND(CHOUTP)
      ORBT  = 0.
      ORBS1 = ORBS0
      ORBC1 = ORBC0
      ORBP1 = ORBP0
      ORBSIG= 1.
      CALL SPLINE2DT(ORBR1,ORBS1,ORBC1,1,1,1,
     &     ORBER,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFR)
      CALL SPLINE2DT(ORBZ1,ORBS1,ORBC1,1,1,1,
     &     ORBEZ,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFZ)

      TMP1 = (Q(JS+1)-Q(JS-1))/(CS(JS+1)-CS(JS-1))
      WRITE(CHOUTP,120) CS(JS),Q(JS),TMP1,DFLOAT(IFOW),
     &                  DFLOAT(NORB),0.
      WRITE(CHOUTP,120) ORBE,ORBL,ORBSIGINI,ORBTB,ORBWD,ORBWTH
      WRITE(CHOUTP,120) ORBT,ORBS1,ORBC1,ORBP1,ORBR1,ORBZ1
      DO K=1,NORB
         ORBSA = ORBS1
         ORBCA = ORBC1
         IF (ORBSA.GT.0..AND.ORBSA.LT.1..AND.
     &       ORBCA.GE.0..AND.ORBCA.LE.2.*PI) THEN
            CALL SPLINE2DT(ORBHK1,ORBSA,ORBCA,1,1,1,
     &           ORBEH,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFH)
            IF ((1.-LAM/ORBHK1).LT.0.) ORBSIG = -ORBSIG
            CALL SPLINE2DT(ORBSK1,ORBSA,ORBCA,1,1,1,
     &           ORBES,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFS)
            IF (ORBSIG.GT.0.) THEN
            CALL SPLINE2DT(ORBCK1,ORBSA,ORBCA,1,1,1,
     &           ORBECP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCP)
            CALL SPLINE2DT(ORBPK1,ORBSA,ORBCA,1,1,1,
     &           ORBEPP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPP)
            ELSE
            CALL SPLINE2DT(ORBCK1,ORBSA,ORBCA,1,1,1,
     &           ORBECM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCM)
            CALL SPLINE2DT(ORBPK1,ORBSA,ORBCA,1,1,1,
     &           ORBEPM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPM)
            ENDIF
         ELSE
            ORBSK1 = 0.
            ORBCK1 = 0.
            ORBPK1 = 0.
         ENDIF

         ORBSA = ORBS1 + ORBSK1*ORBDT/2.
         ORBCA = ORBC1 + ORBCK1*ORBDT/2.
         IF (ORBSA.GT.0..AND.ORBSA.LT.1..AND.
     &       ORBCA.GE.0..AND.ORBCA.LE.2.*PI) THEN
            CALL SPLINE2DT(ORBSK2,ORBSA,ORBCA,1,1,1,
     &           ORBES,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFS)
            IF (ORBSIG.GT.0.) THEN
            CALL SPLINE2DT(ORBCK2,ORBSA,ORBCA,1,1,1,
     &           ORBECP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCP)
            CALL SPLINE2DT(ORBPK2,ORBSA,ORBCA,1,1,1,
     &           ORBEPP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPP)
            ELSE
            CALL SPLINE2DT(ORBCK2,ORBSA,ORBCA,1,1,1,
     &           ORBECM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCM)
            CALL SPLINE2DT(ORBPK2,ORBSA,ORBCA,1,1,1,
     &           ORBEPM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPM)
            ENDIF
         ELSE
            ORBSK2 = 0.
            ORBCK2 = 0.
            ORBPK2 = 0.
         ENDIF

         ORBSA = ORBS1 + ORBSK2*ORBDT/2.
         ORBCA = ORBC1 + ORBCK2*ORBDT/2.
         IF (ORBSA.GT.0..AND.ORBSA.LT.1..AND.
     &       ORBCA.GE.0..AND.ORBCA.LE.2.*PI) THEN
            CALL SPLINE2DT(ORBSK3,ORBSA,ORBCA,1,1,1,
     &           ORBES,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFS)
            IF (ORBSIG.GT.0.) THEN
            CALL SPLINE2DT(ORBCK3,ORBSA,ORBCA,1,1,1,
     &           ORBECP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCP)
            CALL SPLINE2DT(ORBPK3,ORBSA,ORBCA,1,1,1,
     &           ORBEPP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPP)
            ELSE
            CALL SPLINE2DT(ORBCK3,ORBSA,ORBCA,1,1,1,
     &           ORBECM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCM)
            CALL SPLINE2DT(ORBPK3,ORBSA,ORBCA,1,1,1,
     &           ORBEPM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPM)
            ENDIF
         ELSE
            ORBSK3 = 0.
            ORBCK3 = 0.
            ORBPK3 = 0.
         ENDIF

         ORBSA = ORBS1 + ORBSK3*ORBDT
         ORBCA = ORBC1 + ORBCK3*ORBDT
         IF (ORBSA.GT.0..AND.ORBSA.LT.1..AND.
     &       ORBCA.GE.0..AND.ORBCA.LE.2.*PI) THEN
            CALL SPLINE2DT(ORBSK4,ORBSA,ORBCA,1,1,1,
     &           ORBES,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFS)
            IF (ORBSIG.GT.0.) THEN
            CALL SPLINE2DT(ORBCK4,ORBSA,ORBCA,1,1,1,
     &           ORBECP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCP)
            CALL SPLINE2DT(ORBPK4,ORBSA,ORBCA,1,1,1,
     &           ORBEPP,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPP)
            ELSE
            CALL SPLINE2DT(ORBCK4,ORBSA,ORBCA,1,1,1,
     &           ORBECM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFCM)
            CALL SPLINE2DT(ORBPK4,ORBSA,ORBCA,1,1,1,
     &           ORBEPM,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFPM)
            ENDIF
         ELSE
            ORBSK4 = 0.
            ORBCK4 = 0.
            ORBPK4 = 0.
         ENDIF

         ORBT =ORBT + ORBDT
         ORBS2=ORBS1+(ORBSK1+2.*ORBSK2+2.*ORBSK3+ORBSK4)*ORBDT/6. 
         ORBC2=ORBC1+(ORBCK1+2.*ORBCK2+2.*ORBCK3+ORBCK4)*ORBDT/6. 
         ORBP2=ORBP1+(ORBPK1+2.*ORBPK2+2.*ORBPK3+ORBPK4)*ORBDT/6. 

C        CHECK AGAIN PARTICLE POSITION
         IF (ORBS2.LT.0.)    ORBS2 = 0.
         IF (ORBS2.GT.1.)    ORBS2 = 1.
         IF (ORBC2.LT.0.)    ORBC2 = ORBC2 + 2.*PI
         IF (ORBC2.GT.2.*PI) ORBC2 = ORBC2 - 2.*PI
         
         CALL SPLINE2DT(ORBR1,ORBS2,ORBC2,1,1,1,
     &        ORBER,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFR)
         CALL SPLINE2DT(ORBZ1,ORBS2,ORBC2,1,1,1,
     &        ORBEZ,CS,RCHI,NRP1,NRP1,NCHI+1,ORBFZ)

         ORBS1 = ORBS2
         ORBC1 = ORBC2
         ORBP1 = ORBP2
         WRITE(CHOUTP,120) ORBT,ORBS1,ORBC1,ORBP1,ORBR1,ORBZ1
      ENDDO
      CLOSE(CHOUTP)
 120  FORMAT(6(E15.8,1X))
      
      DEALLOCATE(ORBES,ORBEH,ORBECP,ORBECM,ORBEPM,ORBER,ORBEZ,
     &           ORBFS,ORBFH,ORBFCP,ORBFCM,ORBFPM,ORBFR,ORBFZ)

      DEALLOCATE(B_2,B_2M,B2S,B2C)
     
      RETURN
      END
