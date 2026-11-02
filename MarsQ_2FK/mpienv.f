! PLEASE NOTE THAT THE STRUCTURE OF MPI COMPUTATION  IS DESIGNED FOR ONE MARS-K RUN NOW
      MODULE MPIENV
      INTEGER ISMPIRUN
      INTEGER RANK
      INTEGER WORLDSIZE
      INTEGER ROOT  
      INTEGER IERP
      INTEGER CHMPI
    
      PARAMETER (ROOT = 0)
      PARAMETER (CHMPI = 99)
    
      CONTAINS
      SUBROUTINE CHECK_MPI_RUN ()
        USE GLOBALM
        USE REORBITM

        IMPLICIT NONE
        INCLUDE 'mpif.h'
        
        CALL MPI_INIT(IERP)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERP)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, WORLDSIZE, IERP)
        
!       CHECK MPI RUN FOR KINETIC CALCULATION
        ISMPIRUN = 0
        IF (IERP.EQ.0.AND.WORLDSIZE.GE.2) THEN
           IF (INCKIN.EQ.1.AND.NREORBIT.LE.1) ISMPIRUN = 1
           IF (INCKIN.EQ.0.AND.NREORBIT.GT.1) ISMPIRUN = 2
           IF (INCKIN.EQ.1.AND.NREORBIT.GT.1) ISMPIRUN = 3
        ENDIF

!        TEST        
!        WORLDSIZE=2
!        ISMPIRUN = 1
!        RANK = 1
!        
        IF (ISMPIRUN.GE.1 .AND. RANK.EQ.ROOT) THEN
           IF (INCKIN.EQ.1) WRITE(*,*) 'MPI RUN WITH KINETIC EFFECTS'
           IF (NREORBIT.GT.1) WRITE(*,*) 'MPI RUN WITH REORBIT'
           WRITE(*,*) 'ISMPIRUN =',ISMPIRUN 
           WRITE(*,*) 'NUMBER OF REQUESTED PROCESSES =',WORLDSIZE
        ENDIF        
        
      END SUBROUTINE CHECK_MPI_RUN
    
      SUBROUTINE INIT_SURFACE_QUANTITY()
        USE KINETICM
        IMPLICIT NONE
                
        VX1PARA = 0.
        VX2PARA = 0.
        VQ1PARA = 0.
        VQ2PARA = 0.
        VQ3PARA = 0.
        VX1PERP = 0.
        VX2PERP = 0.
        VQ1PERP = 0.
        VQ2PERP = 0.
        VQ3PERP = 0.
        
        VX1PARAM= 0.
        VX2PARAM= 0.
        VQ1PARAM= 0.
        VQ2PARAM= 0.
        VQ3PARAM= 0.
        VX1PERPM= 0.
        VX2PERPM= 0.
        VQ1PERPM= 0.
        VQ2PERPM= 0.
        VQ3PERPM= 0.
C       KJPCOEFF ADDS TO EVERY SURFACE COEFFICIENT, SO A NEW SWEEP OR MPI
C       TASK MUST START FROM ZERO IN ALL OF THEM
        VDPPARA = 0.
        VDPPERP = 0.
        VDPPARAM= 0.
        VDPPERPM= 0.
        IF (ALLOCATED(VX1DPHI)) THEN
           VX1DPHI = 0.
           VX2DPHI = 0.
           VQ1DPHI = 0.
           VQ2DPHI = 0.
           VQ3DPHI = 0.
           VDPDPHI = 0.
           VX1DPHIM= 0.
           VX2DPHIM= 0.
           VQ1DPHIM= 0.
           VQ2DPHIM= 0.
           VQ3DPHIM= 0.
           VDPDPHIM= 0.
        ENDIF
        IF (ALLOCATED(VX1PARE)) THEN
           VX1PARE = 0.
           VX1PERE = 0.
           VX2PARE = 0.
           VX2PERE = 0.
           VQ1PARE = 0.
           VQ1PERE = 0.
           VQ2PARE = 0.
           VQ2PERE = 0.
           VQ3PARE = 0.
           VQ3PERE = 0.
           VDPPARE = 0.
           VDPPERE = 0.
           VX1PAREM= 0.
           VX1PEREM= 0.
           VX2PAREM= 0.
           VX2PEREM= 0.
           VQ1PAREM= 0.
           VQ1PEREM= 0.
           VQ2PAREM= 0.
           VQ2PEREM= 0.
           VQ3PAREM= 0.
           VQ3PEREM= 0.
           VDPPAREM= 0.
           VDPPEREM= 0.
        ENDIF
      END SUBROUTINE INIT_SURFACE_QUANTITY
    
      SUBROUTINE CHECK_MPI_ERROR ()
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        
        IF (IERP .EQ. 0) RETURN
        
        CALL MPI_FINALIZE(IERP)
        PRINT *,'ERROR: RANK=',RANK
        STOP 'ERROR: RANK='
      END SUBROUTINE CHECK_MPI_ERROR
    
!    FUNCTION CHECK_QUEUE(QUEUE,QUEUESIZE):: RESULT (ISFINISH)
!        IMPLICIT NONE
!        INTEGER I,ISFINISH
!        ISFINISH = 1
!        DO I=1:QUQUESIZE
!            IF (QUEUE(3,I).EQ.0) THEN
!                ISFINISH = 0
!                RETURN
!            ENDIF
!        ENDDO
!    END FUNCTION CHECK_QUEUE 
    
      SUBROUTINE GET_EIGENVALUE(BUFFER_EIGENVALUE)
        IMPLICIT NONE
        INCLUDE 'compam.inc'
        COMPLEX*16,INTENT(OUT)::BUFFER_EIGENVALUE
        BUFFER_EIGENVALUE = OMEGA
      END SUBROUTINE GET_EIGENVALUE
    
      SUBROUTINE SET_EIGENVALUE(BUFFER_EIGENVALUE)
        IMPLICIT NONE
        INCLUDE 'compam.inc'
        COMPLEX*16,INTENT(IN):: BUFFER_EIGENVALUE
            
        OMEGA = BUFFER_EIGENVALUE   
      END SUBROUTINE SET_EIGENVALUE
    
      SUBROUTINE MPI_EXIT()
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INTEGER INDEX
        INTEGER BUFFER_SEND(3)
        INTEGER STAT(MPI_STATUS_SIZE) 
        IF (ISMPIRUN .EQ. 0) RETURN
        
C       DO INDEX=1,WORLDSIZE-1
C           BUFFER_SEND= -1              
C           CALL MPI_SSEND( BUFFER_SEND,3, MPI_INTEGER,INDEX,
C    &                      1, MPI_COMM_WORLD,STAT,IERP)
C           CALL CHECK_MPI_ERROR()
C       ENDDO

        CALL MPI_FINALIZE(IERP)
        
      END SUBROUTINE MPI_EXIT

      SUBROUTINE MPI_OPEN_FILE(RANK)
C     ==========================================================
      IMPLICIT NONE
      INCLUDE 'comioc.inc'

      INTEGER RANK
      CHARACTER(LEN=1024) FILENAME

      WRITE(FILENAME,"(A7,I0.4,A4)") "LOG_MPI",RANK,".OUT"
      OPEN(CHMPI,FILE=TRIM(FILENAME),POSITION='APPEND')

      RETURN
      END SUBROUTINE MPI_OPEN_FILE

      SUBROUTINE MPI_CLOSE_FILE(RANK)
C     ==========================================================
      IMPLICIT NONE
      INCLUDE 'comioc.inc'

      INTEGER RANK

      CLOSE(CHMPI)

      RETURN
      END SUBROUTINE MPI_CLOSE_FILE

      END MODULE MPIENV
