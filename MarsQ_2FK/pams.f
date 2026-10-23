C-----------------------------------------------------------------------
      SUBROUTINE CALPAM
C     =================
C-----------------------------------------------------------------------
C                                                                       
     $           (                                                      
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,WORK,IWORK                                           
     $                                                    )             
C                                                                       
C
      IMPLICIT LOGICAL (A-Z)
C
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX
C
         REAL*8      EPSPAM,EPSDET
C
         COMPLEX*16   AL0,ALAM,ALNORM
         INTEGER   NONCON
C
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)
C
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)
C
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)
C-----------------------------------------------------------------------
      INTEGER  IDIX,IDIY,IR,IRY
     $        ,IBT1,IBTY1,ICTE1,ICTY1,IDT1,IDTY1,ISC1,ISCY1
     $        ,IBT2,IBTY2,ICTE2,ICTY2,IDT2,IDTY2,ISC2,ISCY2
     $        ,IXOLD,IYOLD,IXPIVO,IYPIVO
     $        ,ITMP1,ITMP2,ITMPY1,ITMPY2
C
      INTEGER IWORKE
      INTEGER IWORK
      COMPLEX*16    WORK(IWORK)
C
      IDIX  =1
      IDIY  =IDIX  +2*MD     *(ND+1)
      IR    =IDIY  +2*MDY    *ND
      IRY   =IR    +2*MD     *(ND+1)
      IBT1  =IRY   +2*MDY    *ND
      IBTY1 =IBT1  +2*MD *MD
      ICTE1 =IBTY1 +2*MD *MD
      ICTY1 =ICTE1 +2*MD *MD
      IDT1  =ICTY1 +2*MD *MDY
      IDTY1 =IDT1  +2*MD
      ISC1  =IDTY1 +2*MDY
      ISCY1 =ISC1  +2*(2*MD+1)*MD
      IBT2  =ISCY1 +2*(2*MDY+1)*MDY
      IBTY2 =IBT2  +2*MD *MD
      ICTE2 =IBTY2 +2*MD *MD
      ICTY2 =ICTE2 +2*MD *MD
      IDT2  =ICTY2 +2*MD *MDY
      IDTY2 =IDT2  +2*MD
      ISC2  =IDTY2 +2*MDY
      ISCY2 =ISC2  +2*(2*MD+1)*MD
      IXOLD =ISCY2 +2*(2*MDY+1)*MDY
      IYOLD =IXOLD +2*MD     *(ND+1)
      IXPIVO=IYOLD +2*MDY    *ND
      IYPIVO=IXPIVO+2*MD
      ITMP1 =IYPIVO+2*MDY
      ITMP2 =ITMP1 +2*MD
      ITMPY1=ITMP2 +2*MD
      ITMPY2=ITMPY1+2*MDY
      IWORKE=ITMPY2+2*MDY-1
C
C  IWORKE = 6*(MD*(ND+1)+MY*N) + 28*MD*MD + 16*MDY*MDY + 26*MD + 18*MDY
C
C   PRESENT VERSION: IWORKE = 2*[ 3*ND*(MD+MDY) + 10*MD*MD+2*MD*MDY
C                                + 4*MDY*MDY + 10*MD + 7*MDY ]
C-----------------------------------------------------------------------
C
      WRITE(*,*)
      WRITE(*,*) ' COMPLEX*16 PAMS PACKAGE CALLED WITH ARGUMENTS'
      WRITE(*,*) '     MD = ',MD
      WRITE(*,*) '    MDY = ',MDY
      WRITE(*,*) '     ND = ',ND
      WRITE(*,*) '  NCASE = ',NCASE
      WRITE(*,*) ' EPSPAM = ',EPSPAM
      WRITE(*,*) ' EPSDET = ',EPSDET
      WRITE(*,*) '    AL0 = ',AL0
C
      WRITE(*,*) '**************************************'
      WRITE(*,*) ' IWORK=',IWORK,' REQUIRED: ',IWORKE
      WRITE(*,*) '**************************************'
      WRITE(*,*)
         IF (IWORK.LT.IWORKE) STOP 'IWORK'
C
      CALL PAMERA( MD,MDY,ND,NXC,NYC,NCASE,NITMAX
     $            ,EPSPAM,EPSDET
     $            ,AL0,ALAM,ALNORM,NONCON
     $            , A, B, C, D, E, F, G, H
     $            ,DX,DY, X, Y
     $            ,WORK(IDIX)  ,WORK(IDIY)  ,WORK(IR)    ,WORK(IRY)
     $            ,WORK(IBT1)  ,WORK(IBTY1) ,WORK(ICTE1) ,WORK(ICTY1)
     $            ,WORK(IDT1)  ,WORK(IDTY1) ,WORK(ISC1)  ,WORK(ISCY1)
     $            ,WORK(IBT2)  ,WORK(IBTY2) ,WORK(ICTE2) ,WORK(ICTY2)
     $            ,WORK(IDT2)  ,WORK(IDTY2) ,WORK(ISC2)  ,WORK(ISCY2)
     $            ,WORK(IXOLD) ,WORK(IYOLD) ,WORK(IXPIVO),WORK(IYPIVO)
     $            ,WORK(ITMP1) ,WORK(ITMP2) ,WORK(ITMPY1),WORK(ITMPY2)
     $                     )
      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE PAMERA(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C..LOCAL VARIABLES
C
         LOGICAL CINVEC
         INTEGER IK
C
C     SET DIMENSIONS OF MATRICES
C
      M = MD
      N = ND
      MY = MDY
C
      NXPNEG = 0
      NYPNEG = 0
      CINVEC = NCASE.GE.0
      IF (NCASE.EQ.-1) GOTO 150
      IF (NCASE.EQ.-2) GOTO 200
      IF (NCASE.EQ.-3) GOTO 300
      WRITE(*,*) ' START OF PAMERA'
       CALL RUNTIM
C
C     0. COPY LEFT-HAND SIDES
C
      WRITE(*,*) ' CALL MAKCOP'
      CALL MAKCOP(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C     2. INITIAL GUESS
C
      WRITE(*,*) ' CALL INITX0'
      IF (NCASE.EQ.0) CALL INITX0(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                   )             

C                                                                       
C
C     3. REDUCTION FOR Y
C
      WRITE(*,*) ' CALL REDUCY'
      CALL REDUCY(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      WRITE(*,*) ' AFTER REDUCY'
      CALL RUNTIM
C
C     4. REDUCTION FOR X
C
      WRITE(*,*) ' CALL REDUCX'
      CALL REDUCX(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      WRITE(*,*) ' AFTER REDUCX'
      CALL RUNTIM
C
C     5. NUMBER OF NEGATIVE EIGENVALUES
C
COMPLEX*16  NEGEIG=NXPNEG+NYPNEG
C
C     6. INVERSE VECTOR ITERATION
C
      WRITE(*,*) ' ITERATE'
      DO  100  IK=1,NITMAX
      NIT=IK
C
C     7. RHS VECTOR U = B * X
C
      NONCON = 0
      IF (NIT.EQ.1.AND.(NCASE.EQ.3.OR.NCASE.EQ.5.OR.NCASE.EQ.9)) 
     &   NONCON=-1

      IF (CINVEC) CALL UBXRHS(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C     8. Y-REDUCTION OF RHS
C
      CALL REDY(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C     9. X-REDUCTION OF RHS
C
         CALL REDX(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C    10. BACKSOLVE FOR X
C
         CALL SOLVEX(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C    11. BACKSOLVE FOR Y
C
         CALL SOLVEY(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
         IF (.NOT.CINVEC) GO TO 120
C
C    12. NORMALIZATION
C
         IF (NCASE.LT.-3.OR.NCASE.GT.9) CALL STPFWD(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
         IF (NCASE.GE.-3.AND.NCASE.LE.10) CALL STPTIME(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C    13. CONVERGENCE CRITERION
C
         IF (NONCON .EQ. 0) GO TO 110
C
  100    CONTINUE
  110    CONTINUE
C
  120    CONTINUE
      WRITE(*,*) ' AFTER ITERATION'
      CALL RUNTIM
C
C    16. PRINT RESULT
C
         CALL OUTEIG(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
        RETURN
 150  CONTINUE
C
C    14. RAYLEIGH QUOTIENT
C
         CALL EIGEN(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
         RETURN
 200  CONTINUE
C
C     CHECK EIGENVALUE
C
      CALL CHECKE(
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             

      RETURN
 300  CONTINUE
C
C     (DX,DY) = A * (X,Y)  (FOR DIAGNOSTIC PURPOSES)
C
      CALL DEQAX(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE DEQAX(
C        ================
C        (DX,DY) = A * (X,Y)
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C     MULTIPLY (R,RY) = A * (X,Y)
C
      CALL MUAMAT(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C      STORE RESULT IN DX AND DY
C
      CALL SCOPYC(MD*(ND+1),R ,1,DX,1)
      CALL SCOPYC(MDY*ND   ,RY,1,DY,1)
C
      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE CHECKE(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
          INTEGER I,J,IFLAG
          REAL*8    ZDIF
C
C---CHECK SECTION: EVALUATE A*X AND COMPARE WITH LAMBDA*B*X
C
C     MULTIPLY DIX = B * X
C
      CALL MUBMAT(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      CALL SCOPYC(MD*(ND+1),R ,1,DIX,1)
      CALL SCOPYC(MDY*ND   ,RY,1,DIY,1)
C
         CALL MUAMAT(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
         WRITE(*,*)
         WRITE(*,*)' -------CHECK SECTION--------'
         WRITE(*,*)
C
         WRITE(*,*)
     & '   I ,   J , A*X -  LAMBDA*B*X(I,J) , X(I,J)   ABS(X-XOLD)'
         WRITE(*,*)
         DO  400  J=1,N+1
         DO  400  I=1,M
         ZDIF = ABS(X(I,J) - XOLD(I,J))
         IFLAG = 0
         IF (ZDIF.GT.EPSPAM) IFLAG = 10
 400     WRITE(*,1400) I,J,(R(I,J)-ALNORM*DIX(I,J)),X(I,J),ZDIF,IFLAG
C
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
     & '   I ,   J , A*Y -  LAMBDA*B*Y(I,J),  Y(I,J),  ABS(Y-YOLD)'
         WRITE(*,*)
         DO  410  J=1,N
         DO  410  I=1,MY
         ZDIF = ABS(Y(I,J) - YOLD(I,J))
         IFLAG = 0
         IF (ZDIF.GT.EPSPAM) IFLAG = 10
 410     WRITE(*,1400) I,J,(RY(I,J)-ALNORM*DIY(I,J)),Y(I,J),ZDIF,IFLAG
C
1400     FORMAT(1X,I4,I5,2(3X,2E11.3),5X,E11.3,5X,I1)
C
      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE MAKCOP(
C        =================
C        COPY LEFT-HAND SIDE AND ORIGINAL SOLUTION VECTOR
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C
         CALL SCOPYC(MD*(ND+1), DX, 1, DIX, 1)
         CALL SCOPYC(MDY*ND   , DY, 1, DIY, 1)
C
         CALL SCOPYC(MD*(ND+1),  X, 1,XOLD, 1)
         CALL SCOPYC(MDY*ND   ,  Y, 1,YOLD, 1)
C
         RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE INITX0(
C        =================
C        INITIALIZE SOLUTION VECTOR BY RANDOM FUNCTION
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C..LOCAL VARIABLES
C
         INTEGER J,K
         COMPLEX*16 TMPLIU
C
         DO  100  K=1,ND+1
         DO  100  J=1,MD
         CALL RANF(TMPLIU)
         X(J,K)=TMPLIU
  100    XOLD(J,K)=X(J,K)
C
         DO  110  K=1,ND
         DO  110  J=1,MDY
         CALL RANF(TMPLIU)
         Y(J,K)=TMPLIU
  110    YOLD(J,K)=Y(J,K)
C
         RETURN
         END
C-----------------------------------------------------------------------
         SUBROUTINE REDUCY(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
      INTEGER         I,J,K,L
C                                                                       
C---------------------------------------------------------------------- 
C
C..LOCAL VARIABLES
C
         COMPLEX*16 DETY1

C
C     ELIMINATE Y
C     ***********
C
C     LOOP FORWARD (TOP->BOTTOM)
C
C     B = B - E * D-1 * F
C
         DO 20  I=1,ND

         CALL FMIND(D(1,1,I),YPIVOT,MY,MDY,SCY1,DETY1,EPSDET,0,1,I)
 10      CALL MXM(E(1,1,I),M,D(1,1,I),MY,CTY1,MY)
         CALL SCOPYC(MD*MDY,CTY1,1,E(1,1,I),1)
         CALL MXM(CTY1,M,F(1,1,I),MY,BTY1,M)
         CALL VMVC(MD*MD,BTY1,B(1,1,I))
C
C     C = C - E * D-1 * G
C
         CALL MXM(CTY1,M,G(1,1,I),MY,BTY1,M)
         CALL VMVC(MD*MD,BTY1,C(1,1,I))
C
C     B = B - H * D-1 * G
C
         CALL MXM(H(1,1,I),M,D(1,1,I),MY,CTY1,MY)
         CALL SCOPYC(MD*MDY,CTY1,1,H(1,1,I),1)
         CALL MXM(CTY1,M,G(1,1,I),MY,BTY1,M)
         CALL VMVC(MD*MD,BTY1,B(1,1,I+1))
C
C     A = A - H * D-1 * F
C
         CALL MXM(CTY1,M,F(1,1,I),MY,BTY1,M)
         CALL VMVC(MD*MD,BTY1,A(1,1,I))
C
 20      CONTINUE
C

         RETURN
         END
C-----------------------------------------------------------------------
         SUBROUTINE REDUCX(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16        ALSNOR,DET1                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
      INTEGER         I,J,K,L
C                                                                       
C---------------------------------------------------------------------- 
C
C     LOOP FORWARD (TOP->BOTTOM)

      DO 20 I=1,ND
         CALL FMIND(B(1,1,I),XPIVOT,M,MD,SC1,DET1,EPSDET,0,1,I)
 10      CALL MXM(A(1,1,I),M,B(1,1,I),M,CTE1,M)
         CALL SCOPYC(MD*MD,CTE1,1,A(1,1,I),1)
         CALL MXM(CTE1,M,C(1,1,I),M,BT1,M)
         CALL VMVC(MD*MD,BT1,B(1,1,I+1))
 20   CONTINUE
C
      I=ND+1
      CALL FMIND(B(1,1,I),XPIVOT,M,MD,SC1,DET1,EPSDET,0,1,I)

      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE UBXRHS(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C     MULTIPLY R=B*X
C
      CALL MUBMAT(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C      FILL IN RHS
C
      CALL SCOPYC(MD*(ND+1),R ,1,DX,1)
      CALL SCOPYC(MDY*ND   ,RY,1,DY,1)

C     COPY NEW SOLUTION INTO OLD ONE
      CALL SCOPYC(MD*(ND+1),X ,1,XOLD,1)
      CALL SCOPYC(MDY*ND   ,Y ,1,YOLD,1)

      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE MUAMAT(
C        =================
C
C        MULTIPLIES SOLUTION VECTOR BY A-MATRIX: (R,RY) = A * (X,Y)
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
      INTEGER         I,J,K,L
C                                                                       
C---------------------------------------------------------------------- 
C
C..LOCAL VARIABLES
C
C
C     MULTIPLY R=A*X
C
      CALL VZEROC(MD*(ND+1),R )
      CALL VZEROC(MDY*ND   ,RY)
C
      CALL MXV(B(1,1,1),M,X(1,1),M,R(1,1))

         DO 120 I=1,ND
         L=I
         CALL MXV(A(1,1,L),M,X(1,I),M,R(1,I+1))
         CALL MXV(B(1,1,L+1),M,X(1,I+1),M,TMP1)
         CALL VPV(M,TMP1,R(1,I+1))
         CALL MXV(C(1,1,L),M,X(1,I+1),M,TMP2)
         CALL VPV(M,TMP2,R(1,I))
C
         CALL MXV(E(1,1,L),M,Y(1,I),MY,TMP1)
         CALL VPV(M,TMP1,R(1,I))
         CALL MXV(H(1,1,L),M,Y(1,I),MY,TMP2)
         CALL VPV(M,TMP2,R(1,I+1))
C
         CALL MXV(F(1,1,L),MY,X(1,I),M,RY(1,I))
         CALL MXV(G(1,1,L),MY,X(1,I+1),M,TMPY1)
         CALL VPV(MY,TMPY1,RY(1,I))
         CALL MXV(D(1,1,L),MY,Y(1,I),MY,TMPY2)
 120     CALL VPV(MY,TMPY2,RY(1,I))
C

      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE REDY(
C        ===============
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
      INTEGER         I,J,K,L
C                                                                       
C---------------------------------------------------------------------- 
C
C     ELIMINATE Y
C     ***********
C
C     LOOP FORWARD (TOP->BOTTOM)
C
C
C     DX = DX - E * D-1 * DY
C
         DO  110  I=1,ND
         L=I
         CALL MXV(E(1,1,L),M,DY(1,I),MY,DTY1)
         CALL VMVC(M,DTY1,DX(1,I))
C
C     DX = DX - H * D-1 * DY
C
         CALL MXV(H(1,1,L),M,DY(1,I),MY,DTY1)
  110    CALL VMVC(M,DTY1,DX(1,I+1))
C

         RETURN
         END
C-----------------------------------------------------------------------
         SUBROUTINE REDX(
C        ===============
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+2),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
      INTEGER         I,J,K,L
C                                                                       
C---------------------------------------------------------------------- 
C
C     REDUCTION FOR X
C     ***************
C
C     LOOP FORWARD (TOP->BOTTOM)
C
         DO 110  I=1,ND
         L=I
         CALL MXV(A(1,1,L),M,DX(1,I),M,DT1)
  110    CALL VMVC(M,DT1,DX(1,I+1))
C
         RETURN
         END
C-----------------------------------------------------------------------
         SUBROUTINE SOLVEX(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+2)              
         COMPLEX*16   G(MDY,MD,ND+1),H(MD,MDY,ND)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
      INTEGER         I,J,K,L
C                                                                       
C---------------------------------------------------------------------- 
C        BACKSOLVE FOR X
C        ***************
C
C        BACKSOLVE THE 2*2 SYSTEM
C
      L = ND+1
      CALL MXV(B(1,1,L),M,DX(1,N+1),M,X(1,N+1))
C
C     RESOLUTION (BOTTOM->TOP)
C
         DO 10 I=ND,1,-1
         L = I
         CALL MXV(C(1,1,L),M,X(1,I+1),M,DT1)
         CALL VMVC(M,DT1,DX(1,I))
 10      CALL MXV(B(1,1,L),M,DX(1,I),M,X(1,I))
C

         RETURN
         END
C-----------------------------------------------------------------------
         SUBROUTINE SOLVEY(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+2),F(MDY,MD,ND+1)              
         COMPLEX*16   G(MDY,MD,ND+1),H(MD,MDY,ND)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY

      INTEGER         I,J,K,L
C                                                                       
C---------------------------------------------------------------------- 
C
C
C     FORWARDSOLVE FOR Y
C     ***************
C     RESOLUTION (TOP->BOTTOM)
C
         DO 110 I=1,ND
         L = I
         CALL MXV(F(1,1,L),MY,X(1,I),M,DTY1)
         CALL VMVC(MY,DTY1,DY(1,I))
         CALL MXV(G(1,1,L),MY,X(1,I+1),M,DTY1)
         CALL VMVC(MY,DTY1,DY(1,I))
  110    CALL MXV(D(1,1,L),MY,DY(1,I),MY,Y(1,I))
C
         RETURN
         END
C-----------------------------------------------------------------------
         COMPLEX*16 FUNCTION CNORM(M0,N0,V1,V2,IXY,NCOMP)
C        =============================================
C-----------------------------------------------------------------------
C
      INTEGER M0,N0,IXY,NCOMP
      COMPLEX*16 V1(M0,N0),V2(M0,N0)
      INTEGER I,J
C
      CNORM=0.
C
      IF (IXY.EQ.3)     GOTO 300
      IF (IXY.EQ.2)     GOTO 200
      IF (IXY.EQ.1)     GOTO 100
      WRITE(*,*) ' IXY=',IXY
      STOP 'IXY'
100   CONTINUE
      DO 110 J=1,M0,NCOMP
      DO 110 I=1,N0
110   CNORM = CNORM+DFLOAT(I-1)/DFLOAT(N0)* CONJG(V1(J,I))*V2(J,I)

      DO 120 J=2,M0,NCOMP
      DO 120 I=1,N0
120   CNORM = CNORM+DFLOAT(I-1)/DFLOAT(N0)* CONJG(V1(J,I))*V2(J,I)

50    CONTINUE
      RETURN
      DO 130 J=1,M0
      DO 130 I=4,N0-3
130   CNORM = CNORM+DFLOAT(I-1)/DFLOAT(N0)* CONJG(V1(J,I))*V2(J,I)
      RETURN
200   CONTINUE
      IF (.TRUE.) RETURN
      DO 210 J=1,M0-NCOMP+1,NCOMP
C        THIS IS A SPECIAL ARRANGEMENT FOR THE CASE WITH OMEGA
      DO 210 I=1,N0
210   CNORM=CNORM+   (DFLOAT(I-1)/DFLOAT(N0))**3
     &      * CONJG(V1(J,I))*V2(J,I)+CONJG(V1(J+1,I))*V2(J+1,I)
300   CONTINUE
 
C     COMPUTE NORM USING ALL EIGENFUNCTIONS,
C     FOR COMPUTING CRITERION OF ADAPTIVE TIME STEPPING

      DO 310 J=1,M0
      DO 310 I=1,N0
310   CNORM = CNORM+DFLOAT(I-1)/DFLOAT(N0)* CONJG(V1(J,I))*V2(J,I)

      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE STPFWD(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+2),E(MD,MDY,ND+1),F(MDY,MD,ND+1)              
         COMPLEX*16   G(MDY,MD,ND+1),H(MD,MDY,ND)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C..LOCAL VARIABLES
C
         COMPLEX*16 XSCXO,ZNORMO
         REAL*8    XSCX,ZSCX,ZNORM,AOMEGA
         INTEGER J,K,INC
         COMPLEX*16 CNORM
C
      AOMEGA = ALAM

      IF (NITMAX.EQ.0) GOTO 500
C     THIS IS FOR USE AS LINEAR RESPONSE CODE
C
C     NORMALIZE VECTOR
C
         XSCX =  DREAL(CNORM(M,(N+1),X,X,   1,NXC)
     &               +CNORM(MY,N   ,Y,Y,   2,NYC))
         XSCXO=       CNORM(M,(N+1),X,XOLD,1,NXC)
     &               +CNORM(MY,N   ,Y,YOLD,2,NYC)
C
      IF (ABS(XSCX).LE.1.0e-16) THEN
         XSCX = 1.0
         XSCXO = 1.0
      ENDIF
      ALSNOR = XSCXO/XSCX
      ALNORM=AL0+ALSNOR
      ZNORM = SQRT(1./XSCX)
      ZNORMO = ZNORM/ALSNOR

      IF (AOMEGA.LT.0.0.OR.AOMEGA.GT.1.0e+10) THEN
         ZNORMO = 1.0
         ZNORM  = 1.0
      ENDIF

      IF (NITMAX.EQ.1) ZNORM=1.0

      WRITE(*,1002) NIT,XSCX,XSCXO,ALNORM
 1002 FORMAT(I4,5E16.8)
C
C     REDEFINE RHS
C
         NXDEV=0
         NYDEV=0
C
      IF (AOMEGA>1.0e+5) THEN
         DO  100  K=1,ND+1
         DO  100  J=1,MD
         XOLD(J,K)=XOLD(J,K)*ZNORMO
 100     X(J,K)=X(J,K)*ZNORM
      ENDIF 
C
      DO 110 K=1,ND+1
      DO 110 J=1,MD,NXC
      INC = 0
      IF (ABS(XOLD(J,K)-X(J,K)) .GT. EPSPAM) INC = 1
      NXDEV=NXDEV+INC
 110  CONTINUE
C
      IF (AOMEGA>1.0e+5) THEN
         DO  200  K=1,ND
         DO  200  J=1,MDY
         YOLD(J,K)=YOLD(J,K)*ZNORMO
 200     Y(J,K)=Y(J,K)*ZNORM
      ENDIF
C
C     CALL CSSCAL(MDY*ND,ZNORM ,Y   ,1)
C     CALL CSCAL (MDY*ND,ZNORMO,YOLD,1)
C
      DO 220 K=1,ND
      DO 210 J=1,MDY,NYC
      INC = 0
      IF (ABS(YOLD(J,K)-Y(J,K)) .GT. EPSPAM) INC = 1
 210  NYDEV=NYDEV+INC
      DO 220 J=2,MDY,NYC
      INC = 0
      IF (ABS(YOLD(J+1,K)-Y(J+1,K)) .GT. EPSPAM) INC = 1
      NYDEV=NYDEV+INC
 220  CONTINUE
C
      NONCON = NXDEV + NYDEV

      IF (NIT.EQ.NITMAX) RETURN
C
      CALL SCOPYC(MD*(ND+1),X ,1,XOLD,1)
      CALL SCOPYC(MDY*ND   ,Y ,1,YOLD,1)
C
      RETURN
C
C     RESPONSE CALCULATION
C
 500  CONTINUE
         XSCX =  DREAL(CNORM(M,(N+1),X,   X,   1,NXC)
     &               +CNORM(MY,N   ,Y,   Y,   2,NYC))
         ZSCX =  DREAL(CNORM(M,(N+1),XOLD,XOLD,1,NXC)
     &               +CNORM(MY,N   ,YOLD,YOLD,2,NYC))
C
      ZNORM = XSCX/ZSCX
      WRITE(*,*) ' RESPONSE =',ZNORM
      ALNORM = ZNORM
C
      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE STPTIME(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+2)                
         COMPLEX*16   D(MDY,MDY,ND+1),E(MD,MDY,ND+1),F(MDY,MD,ND+1)              
         COMPLEX*16   G(MDY,MD,ND+1),H(MD,MDY,ND)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C..LOCAL VARIABLES
C
      COMPLEX*16 XSCXO,COEF1,ZNORM0
      REAL*8     XSCX,ZSCX,ZNORM
      INTEGER    J,K,IFLAG,INORMSOL
      COMPLEX*16 CNORM,TMP
 
C     COMPUTE SOLUTION AT THE NEXT ITERATION
      COEF1  = (1./ALAM-1.)
      
      DO  100  K=1,ND+1
         DO  100  J=1,MD
 100        X(J,K) = X(J,K) - COEF1*XOLD(J,K)
      
      DO  200  K=1,ND
         DO  200  J=1,MDY
 200        Y(J,K) = Y(J,K) - COEF1*YOLD(J,K)

C     COMPUTE NORM OF THE NEW SOLUTION
      XSCX = DREAL(CNORM(M,N+1,X,X,1,NXC)
     &            +CNORM(MY,N ,Y,Y,2,NYC))
      XSCXO= CNORM(M,N+1,X,XOLD,1,NXC)
     &            +CNORM(MY,N,Y,YOLD,2,NYC)

      CALL GET_INORMSOL(INORMSOL)

C     COMPUTE EIGENVALUE AT CURRENT ITERATION
      IF (NCASE.EQ.0.OR.NCASE.EQ.1.OR.(NCASE.EQ.4.AND.INORMSOL.EQ.1)) 
     &THEN
         ALSNOR = XSCXO/XSCX
         ALNORM=AL0+ALSNOR
         WRITE(*,1002) NIT,XSCX,XSCXO,ALNORM
 1002    FORMAT(I4,5E16.8)
      ENDIF

C     NORMALISE THE OLD/NEW SOLUTIONS
      IF (INORMSOL.EQ.1) THEN
         ZNORM  = SQRT(1./XSCX)
         ZNORM0 = ZNORM
         IF (ABS(ALSNOR).GT.0.) ZNORM0 = ZNORM/ALSNOR
         DO  110  K=1,ND+1
         DO  110  J=1,MD
         XOLD(J,K)=XOLD(J,K)*ZNORM0
 110     X(J,K)=X(J,K)*ZNORM

         DO  210  K=1,ND
         DO  210  J=1,MDY
         YOLD(J,K)=YOLD(J,K)*ZNORM0
 210     Y(J,K)=Y(J,K)*ZNORM
      ENDIF 

C     COMPUTE CONVERGENCE INDEX
      NONCON = 3
      IFLAG  = 0

      IF (INORMSOL.EQ.1) THEN

      NXDEV = 0
      DO 120 K=1,ND+1
      DO 130 J=1,MD,NXC
 130  IF (ABS(XOLD(J,K)-X(J,K)).GT.EPSPAM) NXDEV=NXDEV+1
      DO 120 J=2,MD,NXC
      IF (ABS(XOLD(J,K)-X(J,K)).GT.EPSPAM) NXDEV=NXDEV+1
 120  CONTINUE

      NYDEV = 0
      DO 220 K=1,ND
      DO 220 J=1,MDY,NYC
      IF (ABS(YOLD(J,K)-Y(J,K)).GT.EPSPAM) NYDEV=NYDEV+1
 220  CONTINUE
      
C     NONCON = NXDEV + NYDEV
      NONCON = NXDEV
      
      ELSEIF (NCASE.EQ.6.OR.NCASE.EQ.7.OR.NCASE.EQ.10.OR.
     &        (NCASE.EQ.4.AND.INORMSOL.EQ.0)) THEN

C     USE XOLD,YOLD TO TEMPORARILY STORE THE SOLUTION DIFFERENCE BETWEEN
C     TWO TIME STEPS
      DO  140  K=1,ND+1
      DO  140  J=1,MD
 140     XOLD(J,K)=X(J,K)-XOLD(J,K)

      DO  240  K=1,ND
      DO  240  J=1,MDY
 240     YOLD(J,K)=Y(J,K)-YOLD(J,K)

      XSCX = DREAL(CNORM(M,N+1,X,X,3,NXC)
     &            +CNORM(MY,N ,Y,Y,3,NYC))
      ZSCX = DREAL(CNORM(M,N+1,XOLD,XOLD,3,NXC)
     &            +CNORM(MY,N ,YOLD,YOLD,3,NYC))
      ZSCX = SQRT(ZSCX/XSCX)

      CALL GET_ADAPCRIT(ZSCX,NONCON)

C     RECOVER THE OLD SOLUTION XOLD,YOLD
      DO  141  K=1,ND+1
      DO  141  J=1,MD
 141     XOLD(J,K)=X(J,K)-XOLD(J,K)

      DO  241  K=1,ND
      DO  241  J=1,MDY
 241     YOLD(J,K)=Y(J,K)-YOLD(J,K)

      ENDIF
            
C     SOLVE TOROIDAL MOMEMTUM EQUATION
C     AND DENSITY RADIAL TRANSPORT EQUATION
      IF (NCASE.GE.5.AND.NCASE.LE.8) THEN
         IF (NCASE.NE.8) THEN
            CALL TORQJXB
            CALL TORQNTV
            CALL TORQREY
            CALL TORQERGO
            CALL CALCDENS
            CALL CALCDISPNORM
         ENDIF
         CALL SOLVEMOMENT
         CALL SOLVEDNTRAN
      ENDIF

C     SAVE TIME-STEPPING DATA
      IF (NCASE.GE.3.AND.NCASE.LE.10.AND.INORMSOL.EQ.0) THEN
         CALL FEEDEVOL(MD,MDY,ND,NCASE,X,Y,XOLD,YOLD,ATAU,IFLAG)
         IF (IFLAG.EQ.1.AND.NCASE.EQ.4) NONCON=0
      ENDIF
 
      RETURN
      END
C-----------------------------------------------------------------------
         SUBROUTINE OUTEIG(
C        =================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+3),C(MD,MD,ND+1)                
         COMPLEX*16   D(MDY,MDY,ND+1),E(MD,MDY,ND+1),F(MDY,MD,ND+1)              
         COMPLEX*16   G(MDY,MD,ND+1),H(MD,MDY,ND)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
         PRINT 100,AL0,ALNORM
  100    FORMAT (1H1,/1X,'EIGENVALUE SHIFT              =',1P2E15.5,
C    F              //1X,'RAYLEIGH QUOTIENT             =',2E15.5,
     F              //1X,'EIGENVALUE FROM NORMALIZATION =',2E15.5)
         PRINT 110,NITMAX,NIT
C        WRITE(16,110) NITMAX,NIT
  110    FORMAT (//1X,'MAX. NUMBER OF ITERATIONS =',I5,
     F           //1X,'ITERATIONS DONE           =',I5)
COMPLEX*16  PRINT 120,NEGEIG,NXPNEG,NYPNEG
C        WRITE(16,120) NEGEIG,NXPNEG,NYPNEG
C 120    FORMAT (//1X,'NUMBER OF NEGATIVE EIGENVALUES =',I20,
C    F           //1X,'NUMBER OF NEGATIVE X - PIVOTS  =',I20,
C    F           //1X,'NUMBER OF NEGATIVE Y - PIVOTS  =',I20)
         PRINT 130,NXDEV,NYDEV
C        WRITE(16,130) NXDEV,NYDEV
  130    FORMAT (//1X,'NUMBER OF NON CONVERGED X COMPONENTS =',I5,
     F           //1X,'NUMBER OF NON CONVERGED Y COMPONENTS =',I5)
         PRINT 140,EPSPAM
C        WRITE(16,140) EPSPAM
  140    FORMAT (//1X,'CONVERGED IF ERROR SMALLER THAN ',1PE13.5)
         PRINT 150
C        WRITE(16,150)
  150    FORMAT (1H1)
C
         RETURN
         END
C-----------------------------------------------------------------------
         SUBROUTINE EIGEN(
C        ================
C-----------------------------------------------------------------------
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
      IMPLICIT LOGICAL (A-Z)                                            
C     THIS IS TAKEN OUT TO AVOID DOUBLE DECLARATION IN MHDLIB           
C                                                                       
         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX                       
C                                                                       
         REAL*8      EPSPAM,EPSDET                                        
C                                                                       
         COMPLEX*16   AL0,ALAM,ALNORM                                      
         INTEGER   NONCON                                               
C                                                                       
         COMPLEX*16   A(MD,MD,ND+2),B(MD,MD,ND+2),C(MD,MD,ND+1)                
         COMPLEX*16   D(MDY,MDY,ND+1),E(MD,MDY,ND+1),F(MDY,MD,ND+1)              
         COMPLEX*16   G(MDY,MD,ND+1),H(MD,MDY,ND+2)                            
C                                                                       
         COMPLEX*16   DX(MD,ND+1),DY(MDY,ND+1)                               
C                                                                       
         COMPLEX*16   X(MD,ND+1),Y(MDY,ND+1)                                 
C                                                                       
         COMPLEX*16   DIX(MD,ND+1),DIY(MDY,ND)                             
         COMPLEX*16   R(MD,ND+1),RY(MDY,ND)                                
C                                                                       
         COMPLEX*16   BT1(MD,MD),BTY1(MD,MD),CTE1(MD,MD),CTY1(MD,MDY)      
         COMPLEX*16   DT1(MD),DTY1(MDY)                                    
         COMPLEX*16   SC1(MD,2*MD+1),SCY1(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   BT2(MD,MD),BTY2(MD,MD),CTE2(MD,MD),CTY2(MD,MDY)      
         COMPLEX*16   DT2(MD),DTY2(MDY)                                    
         COMPLEX*16   SC2(MD,2*MD+1),SCY2(MDY,2*MDY+1)                     
C                                                                       
         COMPLEX*16   XOLD(MD,ND+1),YOLD(MDY,ND)                           
         COMPLEX*16   XPIVOT(MD),YPIVOT(MDY)                               
C                                                                       
         COMPLEX*16   TMP1(MD),TMP2(MD),TMPY1(MDY),TMPY2(MDY)              
C                                                                       
      REAL*8            RNORM,XNORM                                       
      COMPLEX*16         ALSNOR                                            
      COMMON /COMPAM/ RNORM,XNORM,ALSNOR                                
C                                                                       
      INTEGER         NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
      INTEGER         NEV, M, N, MY
      COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG          
     $               ,NEV, M, N, MY
C                                                                       
C---------------------------------------------------------------------- 
C
C..LOCAL VARIABLES
C
         COMPLEX*16 ZA
         COMPLEX*16 CNORM
C
C     COMPUTE RAYLEIGH QUOTIENT
C
C     MULTIPLY R=A*X
C
      CALL MUAMAT(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C     MULTIPLY ALAM = X(T) * R
C
         ALAM =       CNORM(M,(N+1),X,R,   1,NXC)
     &               +CNORM(MY,N   ,Y,RY,  2,NYC)
C
C     MULTIPLY R = B * X
C
      CALL MUBMAT(
C                                                                       
     $             MD,MDY,ND,NXC,NYC,NCASE,NITMAX                       
     $            ,EPSPAM,EPSDET                                        
     $            ,AL0,ALAM,ALNORM,NONCON                               
     $            , A, B, C, D, E, F, G, H                              
     $            ,DX,DY, X, Y                                          
     $            ,DIX,DIY,R,RY                                         
     $            ,BT1,BTY1,CTE1,CTY1,DT1,DTY1,SC1,SCY1                 
     $            ,BT2,BTY2,CTE2,CTY2,DT2,DTY2,SC2,SCY2                 
     $            ,XOLD,YOLD,XPIVOT,YPIVOT                              
     $            ,TMP1,TMP2,TMPY1,TMPY2                                
     $                                                    )             
C                                                                       
C
C     MULTIPLY ZA = X(T) * R
C
         ZA   =       CNORM(M,(N+1),X,R,   1,NXC)
     &               +CNORM(MY,N   ,Y,RY,  2,NYC)
         ALAM=ALAM/ZA
C
C     CALCULATE ALNORM
C
         ALNORM=AL0+ALSNOR
         ALAM  =AL0+ALAM
C
         PRINT 1000,ALAM,ALNORM
 1000    FORMAT (1H1,/1X,'RAYLEIGH QUOTIENT             =',2E15.5,
     F              //1X,'EIGENVALUE FROM NORMALIZATION =',2E15.5)
C
         RETURN
         END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                 C
C INPUT:                                                          C
C    A    : A GIVEN COMPLEX*16 MATRIX A(LDA,N+M).                    C
C    N    : # OF EQUATIONS.(IE. # OF ROWS OF A*X = B)             C
C           A EQ. A(1:N,1:N). B EQ. A(1:N,N+1:N+M)                C
C    LDA  : LEADING DIMENSION OF MATRIX A.                        C
C    SCR  : A GIVEN WORKING SPACE,SCR(N).                         C
C    EPS  : A GIVEN LOWER LIMIT FOR DETERMINANT PARTIAL PRODUCT.  C
C           INVERSION IS DECLARED SINGULAR IF AT ANY PRODUCTION   C
C           STEP THE PARTIAL PRODUCT OF PIVOT ELEMENTS IS LESS    C
C           THAN OR EQUAL TO EPS IN MAGNITUDE.                    C
C    M    : # OF SOURCE COLUMNS OF B, AND # OF SOLUTIONS OF X.    C
C    MODE : IF MODE = 1, BOTH THE INVERSE(A) AND SOLUTIONS ARE    C
C                        RETURNED.                                C
C           IF MODE = 0, ONLY THE SOLUTIONS ARE RETURNED.         C
C OUTPUT:                                                         C
C    DET  : THE DETERMINANT OF A(1:N,1:N).                        C
C          D(1:N) CONTAINS THE PIVOTS OF A                        C
C    MODE = 1                                                     C
C           (A) A(1:N,1:N) = INVERSE(ORIGINAL(A(1:N,1:N)))        C
C           (B) A(1:N,N+1,N+M) CONTAINS THE SOLUTIONS OF AX=B.    C
C    MODE = 0                                                     C
C           (A) A(1:N,N+1,N+M) CONTAINS THE SOLUTIONS OF AX=B.    C
C           (B) A(1:N,1:N) IS DESTROYED.                          C
C EXIT(1) :                                                       C
C           IF THE DET OF PARTIAL PRODUCT OF PIVOTS IS LESS       C
C           THAN OR EQUAL TO EPS, THEN EXIT(1).                   C
C                                                                 C
C     IPSI SERVES TO IDENTIFY FLUX SURFACE IN CASE OF ZERO-DIVIDE C
C     WHICH GIVES TRACE-BACK ON CRAY-2 CFT77                      C
C     THIS VERSION USES NAG ROUTINES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C-----------------------------------------------------------------------
      SUBROUTINE FMIND(A,D,N,NDIM,SCR,DET,EPS,M,MODE,IPSI)
C     ===================================================
C-----------------------------------------------------------------------
      IMPLICIT LOGICAL (A-Z)
      INTEGER          N,NDIM,M,MODE,IPSI
      COMPLEX*16       A(NDIM,N),D(N),SCR(N,2*N+1),DET
      REAL*8           EPS
      INTEGER          J,K
      REAL*8           DETR,DETI
      INTEGER          ID,IFAIL
      INTEGER,DIMENSION(:),ALLOCATABLE::IPIV
C
      ALLOCATE(IPIV(N))
C
      IF (MODE.NE.1) STOP 'MODE'
      IF (M.NE.0) STOP 'M.NE.0'
C
      CALL VZEROC(N*N,SCR)
C
      DO 10 J=1,N
 10   SCR(J,J)=1.
C
C      CALL F03AHF(N,A,NDIM,DETR,DETI,ID,SCR(1,N+1),IFAIL)
C      CALL F03AHE(N,A,NDIM,DETR,DETI,ID,SCR(1,N+1),IFAIL)

       CALL ZGETRF(N,N,A,NDIM,IPIV,IFAIL)

      IF (IFAIL.EQ.0) GOTO 20
      WRITE(*,*) ' ZGETRF ERROR: IPSI,IFAIL =',IPSI,IFAIL
      ERROR STOP 'FMIND: ZGETRF found a singular block'

C 20   CALL F04AKF(N,N,A,NDIM,SCR(1,N+1),SCR,N)
C 20   CALL F04AKE(N,N,A,NDIM,SCR(1,N+1),SCR,N)
 20   CALL ZGETRS('N',N,N,A,NDIM,IPIV,SCR(1,N+1),NDIM,IFAIL)
      CALL ZGETRS('N',N,N,A,NDIM,IPIV,SCR       ,NDIM,IFAIL)
      CALL SCOPYC(N*N,SCR,1,A,1)
C
      DEALLOCATE(IPIV)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                 C
C INPUT:                                                          C
C    A    : A GIVEN COMPLEX*16 MATRIX A(LDA,N+M).                 C
C    N    : # OF EQUATIONS.(IE. # OF ROWS OF A*X = B)             C
C           A EQ. A(1:N,1:N). B EQ. A(1:N,N+1:N+M)                C
C    LDA  : LEADING DIMENSION OF MATRIX A.                        C
C    SCR  : A GIVEN WORKING SPACE,SCR(N).                         C
C    EPS  : A GIVEN LOWER LIMIT FOR DETERMINANT PARTIAL PRODUCT.  C
C           INVERSION IS DECLARED SINGULAR IF AT ANY PRODUCTION   C
C           STEP THE PARTIAL PRODUCT OF PIVOT ELEMENTS IS LESS    C
C           THAN OR EQUAL TO EPS IN MAGNITUDE.                    C
C    M    : # OF SOURCE COLUMNS OF B, AND # OF SOLUTIONS OF X.    C
C    MODE : IF MODE = 1, BOTH THE INVERSE(A) AND SOLUTIONS ARE    C
C                        RETURNED.                                C
C           IF MODE = 0, ONLY THE SOLUTIONS ARE RETURNED.         C
C OUTPUT:                                                         C
C    DET  : THE DETERMINANT OF A(1:N,1:N).                        C
C          D(1:N) CONTAINS THE PIVOTS OF A                        C
C    MODE = 1                                                     C
C           (A) A(1:N,1:N) = INVERSE(ORIGINAL(A(1:N,1:N)))        C
C           (B) A(1:N,N+1,N+M) CONTAINS THE SOLUTIONS OF AX=B.    C
C    MODE = 0                                                     C
C           (A) A(1:N,N+1,N+M) CONTAINS THE SOLUTIONS OF AX=B.    C
C           (B) A(1:N,1:N) IS DESTROYED.                          C
C EXIT(1) :                                                       C
C           IF THE DET OF PARTIAL PRODUCT OF PIVOTS IS LESS       C
C           THAN OR EQUAL TO EPS, THEN EXIT(1).                   C
C                                                                 C
C     IPSI SERVES TO IDENTIFY FLUX SURFACE IN CASE OF ZERO-DIVIDE C
C     WHICH GIVES TRACE-BACK ON CRAY-2 CFT77                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C-----------------------------------------------------------------------
      SUBROUTINE FMINDR(A,D,N,NDIM,SCR,DET,EPS,M,MODE,IPSI)
C     ===================================================
C-----------------------------------------------------------------------
      IMPLICIT LOGICAL (A-Z)
      INTEGER          N,NDIM,M,MODE,IPSI
      COMPLEX*16          A(NDIM,N),D(N),SCR(N,2*N+1),DET
      REAL*8             EPS
      INTEGER          J
C
      IF (MODE.NE.1) STOP 'MODE'
      IF (M.NE.0) STOP 'M.NE.0'
C
      CALL SCOPYC(N*N,A,1,SCR,1)
      CALL VZEROC(N*N,SCR(1,N+1))
C
      DO 10 J=1,N
 10   SCR(J,J+N)=1.
      CALL GAUSS1(SCR,D,N,N,NDIM,DET,EPS,SCR(1,2*N+1),IPSI)
C     WRITE(*,*) ' IPSI=',IPSI,' LOG(DET)=',DET
      CALL SCOPYC(N*N,SCR(1,N+1),1,A,1)
      RETURN
      END
C----------------------------------------------------------------------
C      SOLVE LINEAR COMPLEX*16 SYSTEM OF EQUATIONS
C      BY GAUSS ELIMINATION WITH ROW PIVOTING
C      M RIGHT-HAND-SIDE VECTORS
C----------------------------------------------------------------------
      SUBROUTINE GAUSS1(A,D,N,M,NDIM,DET,EPS,APIV,IPSI)
      IMPLICIT LOGICAL (A-Z)
      INTEGER   N,M,NDIM,IPSI
      COMPLEX*16   A(NDIM,N+M),D(N),PIV,DET,Z
      REAL*8      APIV(N),EPS
      INTEGER   I,J,II,NXTROW,ICAMAX
C
C     FORWARD SUBSTITUTION
C
      DET=0.
      I=1
10    CONTINUE
C      DO 11 II=I,N
C11    APIV(II)=ABS(A(II,I))
C      MAXI=APIV(I)
C      NXTROW=I
C      DO 12 II=I+1,N
C        IF (APIV(II).LT.MAXI) GOTO 12
C          MAXI=APIV(II)
C          NXTROW=II
C12    CONTINUE
      NXTROW=ICAMAX(N-I+1,A(I,I),1)+I-1
      IF (I.EQ.NXTROW) GOTO 18
        DO 14 II=I,N+M
          PIV=A(I,II)
          A(I,II)=A(NXTROW,II)
          A(NXTROW,II)=PIV
14      CONTINUE
18    CONTINUE
C     WRITE(*,1000) NXTROW-I,ABS(A(I,I))
      D(I)=A(I,I)
      IF (ABS(D(I)).GT.EPS) GOTO 19
      WRITE(*,*) '  GAUSS1, N=',N,' I=',I,' IPSI=',IPSI
     &          ,' DIAG=',ABS(D(I))
      D(I)=EPS
19    CONTINUE
      DET = DET + LOG(D(I))
      Z=1./D(I)
      DO 20 J=I+1,N+M
20      A(I,J)=A(I,J)*Z
C     CALL CSCAL(N+M-I,1./D(I),A(I,I+1),NDIM)
      DO 60 II=I+1,N
C        Z=A(II,I)
C        DO 40 J=I+1,N+M
C40        A(II,J)=A(II,J)-Z*A(I,J)
      CALL CAXPY(N+M-I,-A(II,I),A(I,I+1),NDIM,A(II,I+1),NDIM)
60    CONTINUE
      I=I+1
      IF (I.LT.N) GOTO 10
      D(N)=A(N,N)
      IF (ABS(D(N)).GT.EPS) GOTO 65
      WRITE(*,*) '  GAUSS1, N=',N,' I=',N,' IPSI=',IPSI
     &          ,' DIAG=',ABS(D(N))
      D(N)=EPS
65    CONTINUE
      DET = DET + LOG(D(N))
      Z=1./D(N)
      DO 70 J=N+1,N+M
70    A(N,J)=A(N,J)*Z
C     CALL CSCAL(M,Z,A(N,N+1),NDIM)
C
C
C     BACKWARD SUBSTITUTION
C
      I=N
110   DO 130 II=1,I-1
C      Z=A(II,I)
C      DO 120 J=N+1,N+M
C120   A(II,J)=A(II,J)-Z*A(I,J)
      CALL CAXPY(M,-A(II,I),A(I,N+1),NDIM,A(II,N+1),NDIM)
130   CONTINUE
      I=I-1
      IF (I.GT.1) GOTO 110
      RETURN
1000  FORMAT(' ROW SHIFT',I5,'    ABS(PIVOT)',D12.4)
      END
C-----------------------------------------------------------------------
C.. G.VLAD & A.B.                                              3/4/1990
C..        SUBROUTINE MXM (MATRIX * MATRIX) FOR COMPLEX*16 ELEMENTS
C-----------------------------------------------------------------------
      SUBROUTINE MXM (A, NA, B, NB, C, NC)
C     ====================================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER          I,J,K,NA,NB,NC
      COMPLEX*16          A(NA,NB), B(NB,NC), C(NA,NC)
C
      DO 40 J=1, NC
      DO 10 I=1, NA
10    C(I,J) = 0.
        DO 30 K=1, NB
      CALL CAXPY(NA,B(K,J),A(1,K),1,C(1,J),1)
C          DO 20 I=1, NA
C20        C(I,J)=C(I,J)+A(I,K)*B(K,J)
30      CONTINUE
40    CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C.. G.VLAD  & A.B.                                             3/4/1990
C..        SUBROUTINE MXV (MATRIX * VECTOR) FOR COMPLEX*16 ELEMENTS
C-----------------------------------------------------------------------
      SUBROUTINE MXV (A, NA, B, NB, C)
C     ================================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER          I,K,NA,NB
      COMPLEX*16          A(NA,NB), B(NB), C(NA)
C
      DO 10 I=1, NA
10    C(I) = 0.
        DO 30 K=1, NB
      CALL CAXPY(NA,B(K),A(1,K),1,C,1)
C          DO 20 I=1, NA
C20        C(I)=C(I)+A(I,K)*B(K)
30      CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C.. A. BONDESON          20.10.89
C..        ADD TWO COMPLEX*16 VECTORS  B = B + A
C-----------------------------------------------------------------------
      SUBROUTINE VPV (N,A,B)
C     ======================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER          N,I
      REAL*8             A(*),B(*)
C
      DO 10 I=1,2*N
 10   B(I)=B(I)+A(I)
      RETURN
      END
C-----------------------------------------------------------------------
C.. A. BONDESON          20.10.89
C..    SUBTRACT TWO COMPLEX*16 VECTORS  B = B - A
C-----------------------------------------------------------------------
      SUBROUTINE VMV (N,A,B)
C     ======================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER          N,I
      REAL*8             A(*),B(*)
C
      DO 10 I=1,2*N
 10   B(I)=B(I)-A(I)
      RETURN
      END
      SUBROUTINE VMVC (N,A,B)
C     ======================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER          N,I
      COMPLEX*16       A(*),B(*)
C
      DO 10 I=1,N
 10   B(I)=B(I)-A(I)
      RETURN
      END
C-----------------------------------------------------------------------
C.. A. BONDESON          15.05.90
C..    SET REAL*8 VECTOR TO ZERO
C-----------------------------------------------------------------------
      SUBROUTINE VZERO (N,A)
C     ======================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER          N,I
      REAL*8             A(*)
C
      DO 10 I = 1, N
 10   A(I) = 0.
C
      RETURN
      END
      SUBROUTINE VZEROC (N,A)
C     ======================
C
      IMPLICIT LOGICAL (A-Z)
      INTEGER          N,I
      COMPLEX*16       A(*)
C
      DO 10 I = 1, N
 10   A(I) = 0.
C
      RETURN
      END
         SUBROUTINE RUNTIM
C        #################
C
C U.50   UPDATE CPU TIME (SECS) AND PRINT IT
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         REAL*8
     R     CPTIME,
     R     STIME
         COMMON /COMTIM/ CPTIME,  STIME
         DATA STIME/0/
C
C     SECOND IS A CDC ROUTINE GIVING THE CPU-TIME USED SO FAR (SEC)
C
C        CALL SECOND(CPTIME)
C     CPTIME = SECNDS(0.0)
C
      IF (.TRUE.) RETURN
         CPTIME = CPTIME - STIME
C
         WRITE (*,9900) CPTIME
C
         RETURN
C
 9900    FORMAT(/,1X,'CPU TIME USED SO FAR =',1PE13.4,' SECS')
C
         END
         FUNCTION CVABS(N, X, NX)
C        -----------------------
C
C NORM OF COMPLEX*16 ARRAY X.  X IS INCREMENTED BY NX.
C
	 COMPLEX*16 X
         DIMENSION
     R   X(N*NX)
C
         SD = 0.0
C
         IF (N .LE. 0) GOTO 2
C
         ZR = DREAL(X(1))
	 ZI = DIMAG(X(1))
         SD = ZR * ZR + ZI * ZI
C
         IF (N .EQ. 1) GOTO 2
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
	    ZR = DREAL(X(J))
	    ZI = DIMAG(X(J))
            SD = SD + ZR * ZR + ZI * ZI
   1     CONTINUE
   2     CONTINUE
C
         CVABS = SQRT(SD)
C
         RETURN
         END
C-------------- REPLACEMENTS FOR CRAY LIBRARY BELOW -------
C
         SUBROUTINE CAXPY(N, A, X, NX, Y, NY)
C        ---------- -----
C
C PERFORMS Y = Y + A*X. X AND Y ARE COMPLEX*16 ARRAYS, A IS COMPLEX
C X IS INCREMENTED BY NX AND Y BY NY.
C
	 COMPLEX*16 A, X, Y
         DIMENSION
     R   X(N*NX), Y(N*NY)
C
         IF (A .EQ. 0.0 .OR. N .LE. 0) RETURN
C
         Y(1) = Y(1) + A*X(1)
C
         IF (N .EQ. 1) RETURN
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            Y(J*NY+1) = Y(J*NY+1) + A*X(J*NX+1)
   1     CONTINUE
C
         RETURN
         END
         FUNCTION ICAMAX(N,PV,NX)
C        ------------------------
C
C  FIND ELEMENT WITH MAXIMUM ABSOLUTE VALUE IN COMPLEX*16 ARRAY PV.
C  PV IS INCREMENTED BY NX
C
	 COMPLEX*16   PV
         DIMENSION PV(N*NX)
C
         IF (N .LE. 0) RETURN
C
         ISM = 1
C
         IF (N .EQ. 1) GOTO 2
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            I = J * NX + 1
            IF (ABS(PV(I)) .GT. ABS(PV(ISM))) ISM = I
    1    CONTINUE
    2    CONTINUE
C
         ICAMAX = ISM
C
         RETURN
         END
         SUBROUTINE SCOPY(N, X, NX, Y, NY)
C        ---------------------------------
C
C COPIES REAL*8 ARRAY X INTO REAL*8 ARRAY Y.
C X IS INCREMENTED BY NX AND Y BY NY.
C
C
         REAL*8
     R   X(N*NX), Y(N*NY)
C
         IF (N .LE. 0) RETURN
C
         Y(1) = X(1)
C
         IF (N .EQ. 1) RETURN
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            Y(J*NY+1) = X(J*NX+1)
   1     CONTINUE
C
         RETURN
         END
         SUBROUTINE SCOPYC(N, X, NX, Y, NY)
C        ---------------------------------
C
C COPIES COMPLEX*16 ARRAY X INTO COMPLEX*16 ARRAY Y.
C X IS INCREMENTED BY NX AND Y BY NY.
C
C
         COMPLEX*16
     #   X(N*NX), Y(N*NY)
C
         IF (N .LE. 0) RETURN
C
         Y(1) = X(1)
C
         IF (N .EQ. 1) RETURN
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            Y(J*NY+1) = X(J*NX+1)
   1     CONTINUE
C
         RETURN
         END
C
       SUBROUTINE RANF(TMPLIU)
       INTEGER N
       COMPLEX*16 TMPLIU
       REAL*8     C1,C2,X
       PARAMETER (C1 = 5.1234, C2 = 0.52)
       COMMON /CRANF/ X
       DATA X /0./
C
       X = C1*X + C2
       N = INT(X)
       X = X - DFLOAT(N)
       TMPLIU = X
       RETURN
       END

