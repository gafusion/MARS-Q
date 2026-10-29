C*COMDECK COMSUR
C
***********************************************************************
*                                                                     *
*   C.17  CONSTANT POLOIDAL FLUX SURFACE QUANTITES                    *
*                                                                     *
*                                                                     *
*   ARATIO  : ASPECT RATIO OF FLUX TUBE                               *
*   BETAB   : GENERALIZED BETA ON FLUX SURFACE                        *
*   CDQ     : D(Q) / D(PSI) ON FLUX SURFACE OBTAINED BY INTEGRATION   *
*   CIDI    : TOROIDAL CURRENT INSIDE FLUX SURFACE (NSTTP=4, CINT)    *
*   CIDR    : FLUX SURFACE INTEGRAL USED FOR SPECIFICATION OF         *
*             CURRENT DENSITY (SEE SUBROUTINE CINT)                   *
*   CIDQ    : INTEGRAL FROM CHI=0 TO CHI=2*PI OF J WHERE J IS THE     *
*             JACOBIAN OF THE GENERAL ERATO COORDINATES               *
*   CIDQO   : STORE CIDQ                                              *
*   CID0    : FLUX SURFACE INTEGRAL USED FOR SPECIFICATION OF         *
*             CURRENT DENSITY (SEE SUBROUTINE CINT)                   *
*   CID0O   : STORE CID0                                              *
*   CID2    : FLUX SURFACE INTEGRAL USED FOR SPECIFICATION OF         *
*             CURRENT DENSITY (SEE SUBROUTINE CINT)                   *
*   CID2O   : STORE CID2                                              *
*   CIPR    : I-PRIME                                                 *
*   CDRQ    : RSURF / Q * (D(Q) / D(RSURF))                           *
*   CP      : CONSTANT APPEARING IN J (J IS THE JACOBIAN OF THE       *
*             GENERAL ERATO COORDINATES). CP IS OBTAINED BY NORMA-    *
*             LIZING CHI TO 2*PI                                      *
*   CPDP    : D(CP) / D(PSI) OBTAINED BY INTEGRATION                  *
*   CPPR    : P-PRIME                                                 *
*   CPPRO   : STORE CPPR                                              *
*   CPR     : PRESSURE                                                *
*   DENSTY  : PLASMA DENSITY PROFILE                                  *
*   DPRIME  : DELTA - PRIME (LOWEST A / R0 TERM)                      *
*   D2CIDR  : 2ND DERIVATIVE OF CIDR BY SPLINE INTERPOLATION          *
*   D2CIDQ  : 2ND DERIVATIVE OF CIDQ BY SPLINE INTERPOLATION          *
*   D2CID0  : 2ND DERIVATIVE OF CID0 BY SPLINE INTERPOLATION          *
*   D2CID0O : STORE D2CID0                                            *
*   D2CID2  : 2ND DERIVATIVE OF CID2 BY SPLINE INTERPOLATION          *
*   D2CID2O : STORE D2CID2                                            *
*   D2CIPR  : 2ND DERIVATIVE OF CIPR BY SPLINE INTERPOLATION          *
*   D2CPPR  : 2ND DERIVATIVE OF CPPR BY SPLINE INTERPOLATION          *
*   D2CPPRO : STORE D2CPPR                                            *
*   D2CPR   : 2ND DERIVATIVE OF CPR BY SPLINE INTERPOLATION           *
*   D2TMF   : 2ND DERIVATIVE OF TMF BY SPLINE INTERPOLATION           *
*   D2TMFO  : STORE D2TMF0                                            *
*   D2TTP   : 2ND DERIVATIVE OF TTP BY SPLINE INTERPOLATION           *
*   PSIISO  : CONSTANT FLUX SURFACES USED FOR ERATO MESH              *
*   QPSI    : SAFETY FACTOR                                           *
*   RARE    : POLOIDAL FLUX SURFACE                                   *
*   RB2AV   : AVERAGED TOTAL MAGNETIC FIELD                           *
*   RB2MAX  : MAXIMUM TOTAL MAGNETIC FIELD                            *
*   RDI     : MERCIER CRITERION, LARGE ASPECT RATIO, INCLUDES         *
*             SHAPING TERMS (SEE H.LUETJENS ET AL.,                   *
*                                NUCL.FUSION 32 (1992) P.1625         *     
*   RELL    : ELLIPTICITY OF FLUX TUBE                                *
*   RDEDR   : DERIVATIVE OF ELLIPTICITY OF FLUX TUBE WITH RESPECT TO  *
*             RHO                                                     *
*   RFCIRC  : FRACTION OF CIRCULATING PARTICLES IN HIRSHMAN DEFINITION*
*             OF BOOTSRAP CURRENT                                     *
*   RFPBP   : SURFACE AVERADGED POLOIDAL MAGNETIC FIELD FOR RFP       *
*   RIB2    : INTEGRAL FROM CHI=0 TO CHI=2*PI OF B**2 WHERE J * B IS  *
*             THE MAGNITIC FIELD                                      *
*   RIIE    : INTEGRAL FROM CHI=0 TO CHI=2*PI OF J*(TOROIDAL CURRENT  *
*             DENSITY) / R                                            *
*   RIIR    : INTEGRAL OVER CHI OF J * B_P**2 = C_3                   *
*   RINDUC  : PLASMA INDUCTANCE                                       *
*   RIP     : INTEGRAL FROM CHI=0 TO CHI=2*PI OF J * PRESSURE         *
*   RIPR    : RIIE / CID0                                             *
*   RIP2    : INTEGRAL FROM CHI=0 TO CHI=2*PI OF J * PRESSURE**2      *
*   RIVOL   : CIDR                                                    *
*   RLENG   : ARC-LENGTH OF POLOIDAL FLUX SURFACE                     *
*   RSURF   : GENRALIZED RADIUS RHO                                   *
*   RSY     : MERCIER CRITERION, LARGE ASPECT RATIO, SHAFRANOV-       *
*             YURCHENKO FORMULA                                       *
*   TEMPER  : TEMPERATURE PROFILE                                     *
*   TMF     : TOROIDAL FLUX T                                         *
*   TMFO    : STORE TMF                                               *
*   TTP     : T*T-PRIME                                               *
*   TTPO    : STORE TTP                                               *
*   VSURF   : FLUX AREA INSIDE OF FLUX SURFACE                        *
*                                                                     *
***********************************************************************
C
         COMMON /COMSUR/
     R   ARATIO  ,BETAB   ,CDQ     ,CIDR    ,CIDQ    ,
     R   CID0    ,CID0O   ,CID2    ,CID2O   ,CIPR    ,
     R   CDRQ    ,CP      ,CPDP    ,CPPR    ,CPPRO   ,
     R   CPR     ,DENSTY  ,DPRIME  ,D2CIDR  ,
     R   D2CID0  ,D2CID0O ,D2CID2  ,D2CID2O ,D2CIPR  ,
     R   D2CPPR  ,D2CPPRO ,D2CPR   ,D2TMF   ,D2TMFO  ,
     R   D2TTP   ,PSIISO  ,QPSI    ,RARE    ,RB2AV   ,
     R   RB2MAX  ,RDEDR   ,RDI     ,RELL    ,RFCIRC  ,
     R   RFPBP   ,RIB2    ,RIIE    ,RIIR    ,RINDUC  ,
     R   RIP     ,RIPR    ,RIP2    ,RIVOL   ,RJBSH   ,
     R   RJBSOS,  RJBSR   ,RJDOTB  ,RJPAR   ,RLENG   ,
     R   RLENG1  ,RNUSTAR ,RPEOP   ,RSURF   ,RSY     ,
     R   TEMPER  ,TMF     ,TMFO    ,TTP     ,TTPO    ,
     R   VSURF   ,CIDI

C
         DIMENSION
     R   ARATIO(2*NPISO),  BETAB(NPISO),      CDQ(2*NPISO),      
     R   CIDR(2*NPISO),    CIDQ(2*NPISO),     CID0(2*NPISO),    
     R   CID0O(2*NPISO),   CID2(2*NPISO),     CID2O(2*NPISO),    
     R   CIPR(2*NPISO),    CDRQ(NPISO),       CP(2*NPISO),       
     R   CPDP(2*NPISO),    CPPR(2*NPISO),     CPPRO(2*NPISO),    
     R   CPR(2*NPISO),     DENSTY(2*NPISO),   DPRIME(NPISO),    
     R   D2CIDR(2*NPISO),  D2CID0(2*NPISO),   
     R   D2CID0O(2*NPISO), D2CID2(2*NPISO),   D2CID2O(2*NPISO),  
     R   D2CIPR(2*NPISO),  D2CPPR(2*NPISO),   D2CPPRO(2*NPISO),  
     R   D2CPR(2*NPISO),   D2TMF(2*NPISO),    D2TMFO(2*NPISO),   
     R   D2TTP(2*NPISO),   PSIISO(2*NPISO),   QPSI(2*NPISO), 
     R   RARE(2*NPISO),    RB2AV(2*NPISO),    RB2MAX(2*NPISO),   
     R   RDI(NPISO),       RELL(2*NPISO),     RDEDR(NPISO),      
     R   RFCIRC(2*NPISO),  RFPBP(2*NPISO),    RIB2(2*NPISO),     
     R   RIIE(2*NPISO),    RIIR(2*NPISO),     RINDUC(NPISO),    
     R   RIP(2*NPISO),     RIPR(2*NPISO),     RIP2(2*NPISO),     
     R   RIVOL(2*NPISO),   RJBSH(2*NPISO),    RJBSOS(2*NPISO,2),
     R   RJBSR(2*NPISO),   RJDOTB(2*NPISO),   RJPAR(2*NPISO),    
     R   RLENG(2*NPISO),   RLENG1(2*NPISO),   RNUSTAR(2*NPISO),
     R   RSURF(NPISO),     RSY(NPISO),        
     R   TEMPER(2*NPISO),  TMF(2*NPISO),      TMFO(2*NPISO),     
     R   TTP(2*NPISO),     TTPO(2*NPISO),     VSURF(NPISO),
     R   CIDI(2*NPISO)



