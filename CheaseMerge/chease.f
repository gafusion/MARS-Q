C*DECK CHEASE
         PROGRAM CHEASE
C
C                                        AUTHOR
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
C                      **************************
C                      *                        *
C                      *      C H E A S E       *
C                      *                        *
C                      **************************
C
C
C
C     (C)UBIC (H)ERMITE (E)LEMENT (A)XISYMETIRC (S)TATIC (E)QUILIBRIUM
C
C
C                    CONVERGENCE TESTS ARE SHOWN IN
C
C   [1] H.LUETJENS, A.BONDESON, A.ROY, COMP.PHYS.COMM 69 (1992) P.287-298
C
C   COMPILATION DIRECTIONS: 
C   -----------------------
C   THIS PROGRAM MUST IMPERATIVELY BE COMPILED AT LEAST IN 64BITS
C   PRECISION, I.E. SINGLE PRECISION ON CRAY, DOUBLE-PRECISION
C   ON 32BIT WORKSTATIONS. THIS PROGRAM CONTAINS MANY SUBROUTINES
C   AND THERE ARE COMMANDS UP TO 25 FORTRAN LINES AND DO LOOPS
C   WITH SEVERAL HUNDREDS FORTRAN LINES. THEREFORE, THE USER MUST
C   BE AWARE THAT THE DEFAULT SET-UP OF HIS LOCAL COMPILER CAN BE
C   UNADAPTED FOR THE COMPILATION OF CHEASE. 
C
C   IF MARS EQ'S ARE COMPUTED WITH FFT'S (NFFTOPT=1), THE NAG 
C   LIBRARY MUST BE LINKED. THE CODE USES C06FAE ON 64 BIT MACHINES
C   AND C06FAF ON 32 BITS MACHINES. TO ALLOW THE COMPILATION OF 
C   THE CODE ON COMPUTER'S WITHOUT NAG LIBRARY, A DUMMY C06FAE
C   IS IMPLEMENTED IN THE CODE WHICH STOPS THE RUN AND GIVES A
C   COMPILATION MESSAGE. THE MARS EQ'S CAN BE COMPUTED WITHOUT
C   FFT'S (NFFTOPT=0) WITHOUT THE NAG LIBRARY.
C
C   THE CODE CONTAINS CALLS TO TIME AND DATE ROUTINES WHICH 
C   ARE MACHINE DEPENDENT. PER DEFAULT, THESE CALLS ARE IN
C   COMMENT, AND THE CORRESPONDING FORTRAN LINE STARTS WITH
C   C+DATE. CALLS TO CRAY, SUN AND SILICON GRAPHICS ROUTINES
C   ARE CURRENTLY IMPLEMENTED. IF THE USER REQUIRES INFORMATIONS
C   ABOUT PARTIAL CPU, DATE, ETC..., HE IS ASKED TO REMOVE
C   BY HAND THESE COMMENTS.
C
C   ON CRAY SYSTEMS, CHEASE SHOULD BE COMPILED WITH
C
C   cf77 -Wf"-o agress" -p -l sci 
C
C   ON 32 BITS SUN SPARC-10 WORKSTATIONS WITH
C
C   f77 -Nx300 -Nl30 -r8 -i4 
C
C   ON 32 BITS HP-K200 WORKSTATIONS WITH
C
C   f77 -O +Onolimit +autodblpad 
C   
C   ON 32 BITS IBM-RSK6000 WORKSTATIONS WITH
C
C   xlf -O -qautodbl=dblpad -qaux_size=16384 -qtkq_size=20000 -qst_size=3072
C
C   ON 32 BITS SILICON GRAPHICS WORKSTATIONS WITH
C
C   f77 -O -Nx300 -Nl30 -r8 -i4 chease.f
C
C     TEST CASES : 
C     ------------
C
C     1) SOLOVEV EQUILIBRIUM
C     ----------------------
C        INPUT CHANNEL 5:
C
C          ***
C          ***
C          ***
C          ***
C           $EQDATA
C             NTCASE=1,
C             NCHI=100,NPSI=15,NS= 30,NT= 30,
C           $END
C           $NEWRUN
C           $END 
C
C     2) JET EQUILIBRIUM. J_PHI SPECIFIED WITH TT-PRIME AND P-PRIME
C     -------------------------------------------------------------
C        INPUT CHANNEL 5:
C
C          ***
C          ***
C          ***
C          ***
C           $EQDATA
C             NTCASE=2,
C             NCHI=100,NPSI=15,NS= 30,NT= 30,
C           $END
C           $NEWRUN
C           $END 
C
C     3) NET EQUILIBRIUM. J_PHI SPECIFIED WITH I* AND P-PRIME
C     -------------------------------------------------------
C
C        INPUT CHANNEL 5:
C
C          ***
C          ***
C          ***
C          ***
C           $EQDATA
C             NTCASE=3,
C             NCHI=100,NPSI=15,NS= 30,NT= 30,
C           $END
C           $NEWRUN
C           $END 
C
C     4) ASYMMETRIC EQUILIBRIUM. J_PHI SPECIFIED WITH I* AND P-PRIME
C     --------------------------------------------------------------
C        INPUT CHANNEL 5:
C
C          ***
C          ***
C          ***
C          ***
C           $EQDATA
C             NTCASE=4,
C             NCHI=100,NPSI=15,NS= 30,NT= 30,
C           $END
C           $NEWRUN
C           $END 
C
C
C                            **********
C
C     LIST OF SUBROUTINES
C     ---------------------
C
C     MASTER       CONTROLS THE RUN                        0.01
C
C
C     LABRUN       LABEL THE RUN                           1.01
C     CLEAR        CLEAR ALL COMMONS                       1.02
C     PRESET       SET UP THE DEFAULT CASE                 1.03
C     DATA         READ NAMELIST                           1.04
C     AUXVAL       SET UP AUXILLIARY VALUES                1.05
C     COTROL       CONTROL READ IN PARAMETERS              1.06
C     TCASE        SET UP TEST CASES                       1.07
C
C     STEPON       LEAD THE CALCULATIONS                   2.01
C     BALLIT       LEAD BALLOONING OPTIMZATION AND 
C                  SPECIFICATION OF BOOTSTRAP CURRENT      2.02
C     INITIA       INITIALIZE VERTICAL MATRIX INDEXATION
C                  AND QUADRATURE POINTS FOR EQUILIBRIUM
C                  INTEGRATION                             2.03
C     EQDIM        SET UP SMALL AND FINAL EQUILIBRIUM      2.04
C     MATRIX       LEAD CONSTRUCTION AND LDLT
C                  DECOMPOSITION OF EQ MATRIX              2.05
C     ITIPR        LEAD ITERATION OVER CURRENT PROFILE     2.06
C     NONLIN       LEAD ITERATION OVER NONLINEARITY        2.07
C     CHECK        COMPUTE GLOBAL RESIDU OF A * PSI - B    2.08
C     OLDNEW       STORE  CONVERGED EQUILIBRIUM            2.09
C     OLDEQ        READ CONVERGED EQUILIBRIUM              2.10
C
C     MESH         SET UP DISCRETIZATION MESHES            2.A01
C     PACKME       MESH PACKING WITH LORENTZIANS           2.A02
C     PSVOL        S-MESH PACKING SO THAT D(RHO) / DS = 0  2.A03
C     TETARE       AUTOMATIC THETA-MESH PACKING            2.A04
C     QPLACS       S-MESH PACKING AT PREDEFINED Q-VALUES   2.A05
C     PACKMEP      AS PACKME, BUT FOR 2*PI PERIODIC MESHES 2.A06
C
C     GUESS        INITIALIZE PICARD ITERATION             2.B01
C
C     MAGAXE       FIND MAGNETIC AXIS                      2.C01
C     EVLATE       EVALUATE PSI, D(PSI)/D(R) AND
C                  D(PSI)/D(Z) AT (R,Z)                    2.C02
C
C     SETUPA       CONSTRUCT A                             2.D01
C     SETUPB       CONSTRUCT B                             2.D02
C     LIMITA       IMPOSE BOUNDARY CONDITIONS ON A         2.D03
C     LIMITB       IMPOSE BOUNDARY CONDITIONS ON B         2.D04
C     IDENTA       PERFORM ROW AND COLUMN OPERATIONS IN A  2.D05
C     IDENTB       PERFORM ROW OPERATIONS IN B             2.D06
C     AWAY         REMOVE 1 ROW AND COLUMN IN A            2.D07
C     CENTER       EVALUATE COEFFICENTS REQUIRED TO
C                  IMPOSE BOUNDARY CONDITIONS              2.D08
C
C     SOLVIT       SOLVE GRAD-SHAFRANOV EQUATION           2.E01
C     DIRECT       GAUSS ELIMINATION                       2.E02
C     ERROR        COMPUTE ERROR ON PSI                    2.E03
C     ENERGY       COMPUTE AVERAGED POLOIDAL MAGNETIC 
C                  FIELD ENERGY                            2.E04
C     SMOOTH       BICUBIC SPLINE SMOOTHING OF BICUBIC
C                  HERMITE EQUILIBRIUM SOLUTION            2.E05
C     CONVER       CONVERGENCE TESTS                       2.E06
C
C     NOREPT       EQUILIBRIUM TRANSFORMATIONS             2.F01
C     RSCALE       SCALE EQUILIBRIUM AGAINST R OF
C                  MAGNETIC AXIS                           2.F02
C     TSHIFT       SHIFT TOROIDAL FLUX PROFILE             2.F03
C     PRNORM       SCALE EQUILIBRIUM
C
C     TEST         COMPUTE RELATIVE ERROR FOR SOLOVEV      2.G01
C     SOLOVEV      COMPUTE ANALYTIC SOLOVEV EQUILIBRIUM    2.G02
C
C     MAPPIN       LEAD COMPUTATION OF MAPPINGS            2.M01
C     SURFACE      INTEGRATION OF LOCAL AND GLOBAL FLUX
C                  SURFACE QUANTITIES                      2.M02
C     CHIPSI       INTERPOLATION OF LOCAL FLUX SURFACE
C                  QUANTITES ON ERATO MESH                 2.M03
C     ERDATA       COMPUTE EQ'S FOR ERATO AND LOCAL SHEAR
C                  AND ZERO LINE OF AVERAGED MAGNETIC
C                  FIELD LINE CURVATURE                    2.M04
C     CINT         COMPUTE INTEGRALS NEEDED TO OBTAIN THE
C                  T-TPRIME PROFILE FROM THE I-PRIME AND
C                  THE P-PRIME PROFILE                     2.M05
C     PREMAP       LEAD COMPUTATION OF PROFILES ON S-MESH  2.M06
C     GCHI         INTERPOLATE LOCAL FLUX SURFACE
C                  QUANTITIES ON GAUSS INTEGRATION POINTS  2.M07
C     GIJLIN       COMPUTE LOCAL QUANTITES NEEDED BY MARS  2.M08
C     FOURIER      PERFORM FOURIER TRANSFORM OF QUANTITIES
C                  COMPUTED BY GIJLIN                      2.M09
C     PROFILE      COMPUTE PROFILES ON S-MESH              2.M10
C     BALOON       BALOONING STABILITY AND LOCAL 
C                  INTERCHANGE TESTS                       2.M11
C     GLOQUA       COMPUTE AUXILIARY GLOBAL FLUX SURFACE
C                  QUANTITES                               2.M12
C     VACUMM       COMPUTE VACUUM EQ'S FOR MARS            2.M13
C     VLION        COMPUTE EQ'S FOR LION AT PLASMA SURFACE 2.M14
C     OUTNVW       COMPUTE EQ'S FOR NOVA-W AND PEST        2.M15
C     STCHPS       COMPUTE SIGMA(PSI,CHI) AND
C                          THETA(PSI,CHI) FOR NOVA-W       2.M16
C     JNOVAW       COMPUTE R,Z AND STABILITY MESH 
C                  JACOBIAN AT (SIGMA, THETA)              2.M17
C     TPSI         COMPUTE SIGMA(PSI,THETA-PENN) AND
C                          THETA(PSI,THETA-PENN)           
C                  OR      SIGMA(PSI,THETA-XTOR) AND
C                          THETA(PSI,THETA-XTOR)           2.M18
C     OUTPEN       COMPUTE EQ'S FOR PENN                   2.M19
C     OUTXT        COMPUTE EQ'S FOR XTOR                   2.M20
C
C     FOURFFT      COMPUTE FAST FOURIER TRANSFORMS OF ALL 
C                  EQ'S FOR MARS INSIDE THE PLASMA         2.M22
C     VACUFFT      COMPUTE FAST FOURIER TRANSFORMS OF ALL 
C                  EQV'S FOR MARS IN THE VACUMM            2.M23
C     SPLIFFT      COMPUTES CUBIC SPLINE INTERPOLATION 
C                  AND FAST FOURIER TRANSFORM              2.M24
C
C     CURENT       COMPUTE CURRENT DENSITY FOR GIVEN PSI   2.J01
C                  VALUES
C     ISOFUN       COMPUTE T, T-PRIME, P AND P-PRIME ON
C                  S-MESH                                  2.J02
C
C     PRFUNC       FUNCTIONAL FORM FOR PRESCRIPTION OF
C                  TT'(S),I*(S) OR I_PARA(S)               2.I01
C     ATCOEF       TT'(S),I*(S) OR I_PARA(S) WITH AT'S     2.I02
C     COPYAT       COMPUTE COEFFICIENTS OF POLYNOMIAL
C                  SECTIONS FROM AT'S                      2.I03
C
C     PPRIME       COMPUTE P-PRIME PROFILE                 2.P01
C     BSFUNC       COMPUTE FUNCTION FOR FRACTION OF
C                  BOOTSTRAP CURRENT                       2.P02
C     PPSPLN       SPLINE INTERPOLATION OF P-PRIME         2.P03
C     APCOEF       COMPUTE P-PRIME WITH AP'S               2.P04
C     COPYAP       COMPUTE COEFFICENT OF POLYNOMIAL 
C                  SECTIONS FROM AP'S                      2.P05
C     APCOEF2      ALTERNATIVE VERSION OF APCOEF           2.P06
C     COPYAPP      ALTERNATIVE VERSION OF COPYAPP          2.P07
C     BLTEST       LEAD BALLOONING AND MERCIER STABILITY
C                  CALCULATION FOR BALLOONING OPTIMIZATION 2.P08
C     RESPPR       INITIALIZE PROFILES FOR BALLOONING
C                  OPTIMIZATION                            2.P09
C     PPRM         MODIFICATION P-PRIME PROFILE DURING
C                  BALLOONING OPTIMIZATION                 2.P10
C     PPBSTR       LEAD COMPUTATION OF P-PRIME WHEN
C                  BOOTSTRAP CURRENT DENSITY IS SPECIFIED  2.P11
C
C     POLYNM       COMPUTE POLYNOMIAL COEFFICIENTS OF
C                  DENSITY AND TEMPERATURE WHEN CURRENT
C                  DENSITY IS GIVEN IN TERMS OF DENSITY
C                  AND TEMPERATURE                         2.T01
C     DRHODP       COMPUTE D(RHO)/D(PSI) FOR A GIVEN SET
C                  OF PSI VALUES                           2.T02
C
C     ISOFIND      LEAD TRACING OF CONSTANT FLUX SURFACES  2.U01
C     CUBRT        TRACE CONSTANT FLUX SURFACES            2.U02
C     GAUSS        GAUSS QUADRATURE QUANTITES IN
C                  [0; 1] INTERVAL                         2.U03
C     RMRAD        COMPUTE INTERSECTIONS OF Z OF MAGNETIC
C                  AXIS WITH CONSTANT FLUX SURFACES        2.U04
C
C     BOUND        COMPUTE PLASMA SURFACE                  2.X01
C     BNDINP       READ EXPERIMENTAL BOUNDARY POINTS       2.X02
C     BNDSPL       CUBIC SPLINE INTERPOLATION OF 
C                  EXPERIMENTAL BOUNDARY POINTS            2.X03
C     SUBZ         SHIFT EXPERIMENTAL BOUNDARY VERTICALLY  2.X04
C     RZBOUND      COMPUTE (R,Z) COORDINATES OF BOUNDARY   2.X05
C
C     BASIS1       COMPUTE BASIS FUNCTIONS AT GAUSS 
C                  INTEGRATION POINTS                      2.Y01
C     BASIS2       COMPUTE FIRST DERIVATIVES OF BASIS 
C                  FUNCTIONS AT GAUSS INTEGRATION POINTS   2.Y02
C     BASIS3       COMPUTE FIRST AND 2ND DERIVATIVES OF 
C                  BASIS FUNCTIONS AT GAUSS INTEGRATION 
C                  POINTS                                  2.Y03
C     BASIS4       COMPUTE FIRST DERIVATIVES OF BASIS 
C                  FUNCTIONS IN SIGMA DIRECTION AT GAUSS 
C                  INTEGRATION POINTS                      2.Y04
C     PSICEL       COMPUTE VARIABLES DEFINING THE BICUBIC
C                  EXPANSIONS OF PSI IN KN CELLS           2.Y05
C     PSIBOX       EVALUATE PSI ON A (R,Z) GRID            2.Y06
C
C     OUTPUT       CONTROL INPUT/OUTPUT                    3.A01
C     PRIQQU       PRINT EQUILIBRIUM QUANTITIES Q-VALUES   
C                  SPECIFIED IN QPLACS                     3.A02
C     OUTMKSA      SAVE EQUILIBRIUM QUANTITIES IN MKSA     3.A03
C
C     IODISK       PERFORM DISK FILE OPERATIONS            3.B01
C     WRPLOT       WRITE PLOT QUANTITIES                   3.B02
C     SHAVE        SHAVE AWAY OUTER POLOIDAL FLUX SURFACES 3.B03
C     SURFRZ       SAVE (R,Z)'S OF LAST FLUX SURFACE       3.B04
C     NERAT        COMPUTE ERATO NAMELIST                  3.B05
C     GENOUT       GENERAL OUTPUT ROUTINE USED TO
C                  CONSRUCT INPUT FILE OF LINEAR RESISTIF
C                  CODE                                    3.B06
C
C
C
C                            **********
C
C     AUXILLIARY SUBROUTINES :
C     ------------------------
C
C     ALDLT        DECOMPOSE A = L*D*LT                    MR01
C     LYV          SOLVE  L*Y  = V                         MR02
C     DWY          SOLVE  D*W  = Y                         MR03
C     LTXW         SOLVE  LT*X = W                         MR04
C
C     SPLINE       CUBIC SPLINE INTERPOLATION. BOUNDARY
C                  CONDITIONS BY CUBIC LAGRANGE 
C                  INTERPOLATION                           MSP01
C     MSPLINE      AS SPLINE, BUT DOES M INTERPOLATIONS 
C                  IN PARALLEL                             MSP02
C     SPLCY        CUBIC SPLINE INTERPOLATION. PERIODIC
C                  BOUNDARY CONDITIONS.                    MSP03
C     MSPLCY       AS SPLCY, BUT DOES M INTERPOLATIONS 
C                  IN PARALLEL                             MSP04
C     SPLCYP       CUBIC SPLINE INTERPOLATION OF A PERIODIC
C                  FUNCTION WITH PERIODIC DEFINITION 
C                  INTERVAL                                MSP05
C
C     NTRIDG       LU DECOMPOSE ND2-ND1+1 TRIDIAGONAL
C                  SYSTEMS                                 MRD01
C     TRIDAGM      INVERT M TRIDIAGONAL SYSTEMS IN 
C                  PARALLEL                                MRD02
C     TRICYC       INVERT 1 TRIDIAGONAL SYSTEM WITH
C                  PERIODIC BOUNDARY CONDITIONS            MRD03
C     TRICYCM      INVERT M TRIDIAGONAL SYSTEM WITH
C                  PERIODIC BOUNDARY CONDITIONS IN
C                  PARALLEL                                MRD04
C     TRIDAG       INVERT A TRIDIAGONAL SYSTEM             MRD05
C     SORT         SORTING ROUTINE                         MRD06
C
C     PAGE         JUMP A PAGE                             U1
C     BLINES       JUMP N LINES                            U2
C     MESAGE       WRITE A 48-CHARACTERS STRING            U10
C     RVAR         WRITE NAME AND REAL VALUE               U20
C     RVAR2        WRITE 2*(NAME AND REAL VALUE)           U21
C     RVAR3        WRITE 3*(NAME AND REAL VALUE)           U22
C     IVAR         WRITE NAME AND INTEGER VALUE            U23
C     IVAR2        WRITE 2*(NAME AND INTEGER VALUE)        U24
C     IVAR3        WRITE 3*(NAME AND INTEGER VALUE)        U25
C     HVAR         WRITE NAME AND CHARACTER VALUE          U26
C     LVAR         WRITE NAME AND LOGICAL VALUE            U27
C     RARRAY       WRITE NAME AND REAL ARRAY               U30
C     IARRAY       WRITE NAME AND INTEGER ARRAY            U31
C     LARRAY       WRITE NAME AND LOGICAL ARRAY            U32
C     HARRAY       WRITE NAME AND CHARACTER ARRAY          U33
C     SARRAY       WRITE NAME AND SCALED VALUES OF REAL    U34
C     0ARRAY       WRITE NAME AND REAL ARRAY INTO NCHAN    U35
C     WRTEXT       WRITE REAL VARIABLE NAME AND VALUE INTO
C                  TEXT ARRAY                              U36
C     WITEXT       WRITE INTEGER VARIABLE NAME AND VALUE
C                  INTO TEXT ARRAY                         U37
C     WHTEXT       WRITE TEXT VARIABLE INTI TEXT ARRAY     U38
C     RESETR       RESET REAL ARRAY                        U40
C     RESETI       RESET INTEGER ARRAY                     U41
C     RESETH       RESET CHARACTER ARRAY                   U42
C     RESETL       RESET LOGICAL ARRAY                     U43
C     RESETC       RESET COMPLEX ARRAY                     U44
C     SCOPYR       Y = Y + RF * (X - Y)                    U49
C
C     RUNTIM       UPDATE CPU TIME AND PRINT IT            U50
C     DAYTIM       PRINT OUT DATE AND TIME                 U51
C
C     VZERO        X = 0 (X REAL)                          MAT7
C     CVZERO       X = 0 (X COMPLEX)                       MAT8
C     ACOPY        Y = X (X,Y REAL) WITH INTERMEDIATE
C                        STORAGE                           MAT9
C     ICOPY        Y = X (X,Y INTEGER) WITH INTERMEDIATE
C                        STORAGE                           MAT10
C     CCOPY        Y = X (X,Y COMPLEX) WITH INTERMEDIATE
C                        STORAGE                           MAT11
C
C                            **********
C     COMMENTS (YQL):
C     * IT IS ESSENTIAL TO USE NMESHC=1 TO PACK THE SIGMA-MESH, 
C       IN ORDER TO RESOLVE WELL THE PROFILES NEAR MAGNETIC AXIS.
C
C     CRAY SUBROUTINES AND FUNCTIONS :
C     --------------------------------
C     FOR MORE DETAILS, SEE CRAY LIBRARY MANUAL.
C
C     ISMAX
C     ISMIN
C     ISAMIN
C     ISRCHFGT
C     ISRCHFGE
C     SAXPY  
C     SCOPY  
C     SDOT
C     SSCAL  
C     SSUM
C     ISSUM     
C
C                            **********
C
C     DISK CHANNELS :                     SEE PUBLICATION, TABLE 5
C     ---------------
C     LIST OF COMMONS :                   SEE PUBLICATION, TABLE 8
C     -----------------
C     LIST OF STATEMENT FUNCTIONS DECKS : SEE PUBLICATION, TABLE 7
C     -----------------------------------
C     EQUILIBRIUM NAMELIST VARIABLES :    SEE PUBLICATION, TABLE 9-11
C     --------------------------------
C
C                            **********
C
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C     START THE RUN
C
         CALL MASTER
C
         STOP
         END
C*DECK C0S01
C*CALL PROCESS
         SUBROUTINE MASTER
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C0S01 CONTROLS THE RUN                                              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         COMMON /COMTIM/ STIME
C
C----*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C+DATE   IF (COMPTYP .EQ. 'CRAY') CALL SECOND(STIME)
C
***********************************************************************
*                                                                     *
* 1. SET UP THE CASE                                                  *
*                                                                     *
***********************************************************************
C
         CALL LABRUN
         CALL CLEAR
         CALL PRESET
         CALL DATA
         CALL TCASE
         CALL COTROL
         CALL AUXVAL
C
***********************************************************************
*                                                                     *
* 3. STEP ON THE CALCULATION                                          *
*                                                                     *
***********************************************************************
C
         CALL STEPON
         CALL RUNTIM
C
         RETURN
         END
C*DECK C1S01
C*CALL PROCESS
         SUBROUTINE LABRUN
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C1S01 LABEL THE RUN                                                 *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMLAB.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
***********************************************************************
*                                                                     *
* 1. READ THE LABELS                                                  *
*                                                                     *
***********************************************************************
C
         READ (5,1000) LABEL1
         READ (5,1000) LABEL2
         READ (5,1000) LABEL3
         READ (5,1000) LABEL4
C
***********************************************************************
*                                                                     *
* 2. WRITE THE HEADING                                                *
*                                                                     *
***********************************************************************
C
         CALL PAGE
         CALL BLINES(30)
         WRITE (6,1011)
         WRITE (6,1012)
         WRITE (6,1013)
         WRITE (6,1014)
         WRITE (6,1015)
         WRITE (6,1016)
         WRITE (6,1017)
         WRITE (6,1018)
         WRITE (6,1019)
         WRITE (6,1020)
         WRITE (6,1021)
         WRITE (6,1022)
         WRITE (6,1023)
         CALL BLINES(10)
         WRITE (6,1030)
C
***********************************************************************
*                                                                     *
* 3. WRITE THE LABELS                                                 *
*                                                                     *
***********************************************************************
C
         CALL BLINES(10)
         CALL MESAGE(LABEL1)
         CALL BLINES(1)
         CALL MESAGE(LABEL2)
         CALL BLINES(1)
         CALL MESAGE(LABEL3)
         CALL BLINES(1)
         CALL MESAGE(LABEL4)
         CALL PAGE
C
         RETURN
 1000    FORMAT(A80)
 1011    FORMAT(16X,12HC  CCCCCCCC ,5X,12HH        HH ,5X,
     ,          12HE  EEEEEEEE ,5X,12HA  AAAAAAAA ,5X,
     ,          12HS  SSSSSSSS ,5X,12HE  EEEEEEEE )
 1012    FORMAT(16X,12HCC  CCCCCCCC,5X,12HHH       HHH,5X,
     ,          12HEE  EEEEEEEE,5X,12HAA  AAAAAAAA,5X,
     ,          12HSS  SSSSSSSS,5X,12HEE  EEEEEEEE)
 1013    FORMAT(16X,12HCCC  CCCCCCC,5X,12HHHH      HHH,5X,
     ,          12HEEE  EEEEEEE,5X,12HAAA  AAAAAAA,5X,
     ,          12HSSS  SSSSSSS,5X,12HEEE  EEEEEEE)
 1014    FORMAT(16X,12HCCC      CCC,5X,12HHHH      HHH,5X,
     ,          12HEEE         ,5X,12HAAA      AAA,5X,
     ,          12HSSS         ,5X,12HEEE         )
 1015    FORMAT(16X,12HCCC         ,5X,12HHHH      HHH,5X,
     ,          12HEEE         ,5X,12HAAA      AAA,5X,
     ,          12HSSS         ,5X,12HEEE         )
 1016    FORMAT(16X,12HCCC         ,5X,12HHHHHHHH  HHH,5X,
     ,          12HEEEEEEE     ,5X,12HAAAAAAA  AAA,5X,
     ,          12HSSSSSSSSSSS ,5X,12HEEEEEEE     )
 1017    FORMAT(16X,12HCCC         ,5X,12HHHHHHHHH HHH,5X,
     ,          12HEEEEEEEE    ,5X,12HAAAAAAAA AAA,5X,
     ,          12HSSSSSSSSSSSS,5X,12HEEEEEEEE    )
 1018    FORMAT(16X,12HCCC         ,5X,12HHHHHHHHHHHHH,5X,
     ,          12HEEEEEEEEE   ,5X,12HAAAAAAAAAAAA,5X,
     ,          12HSSSSSSSSSSSS,5X,12HEEEEEEEEE   )
 1019    FORMAT(16X,12HCCC         ,5X,12HHHH      HHH,5X,
     ,          12HEEE         ,5X,12HAAA      AAA,5X,
     ,          12H         SSS,5X,12HEEE         )
 1020    FORMAT(16X,12HCCC      CCC,5X,12HHHH      HHH,5X,
     ,          12HEEE         ,5X,12HAAA      AAA,5X,
     ,          12H         SSS,5X,12HEEE         )
 1021    FORMAT(16X,12HCCCCCCCCCCCC,5X,12HHHH      HHH,5X,
     ,          12HEEEEEEEEEEEE,5X,12HAAA      AAA,5X,
     ,          12HSSSSSSSSSSSS,5X,12HEEEEEEEEEEEE)
 1022    FORMAT(16X,12HCCCCCCCCCCCC,5X,12HHHH      HHH,5X,
     ,          12HEEEEEEEEEEEE,5X,12HAAA      AAA,5X,
     ,          12HSSSSSSSSSSSS,5X,12HEEEEEEEEEEEE)
 1023    FORMAT(16X,12H CCCCCCCCCC ,5X,12H HH       HH,5X,
     ,          12H EEEEEEEEEE ,5X,12H AA       AA,5X,
     ,          12H SSSSSSSSSS ,5X,12H EEEEEEEEEE )
C
 1030    FORMAT(40X,
     ,          'CUBIC HERMIT ELEMENT AXISYMETRIC STATIC EQUILIBRIUM')
C
         END
C*DECK C1S02
C*CALL PROCESS
         SUBROUTINE CLEAR
C        ################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C1S02 CLEAR ALL COMMONS                                             *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMEQD.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMOPT.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMPLO.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
         INCLUDE 'COMVEV.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C         CALL RESETI(NCBAL(1)     ,NBAL1,0 )
         CALL RESETI(NCBAL        ,NPIS0,0 )
         CALL RESETR(ABAL(1,1,1)  ,NBAL2,0.)
C         CALL RESETR(ABAL(1,1,1)  ,NPPSBAL*2*(2*NPTURN*PPCHI+2),0.)
         CALL RESETR(A(1,1)       ,NBLA ,0.)
         CALL RESETR(ALZERO       ,NBND1,0.)
         CALL RESETI(NFOURPB      ,NBND2,0 )
         CALL RESETI(NANAL        ,NICON,0 )
         CALL RESETR(CNR1(1,1)    ,NERA ,0.)
         CALL RESETI(NMAG         ,NESH1,0 )
         CALL RESETR(APLACE(1)    ,NESH2,0.)
         CALL RESETR(ASPCTR       ,NETA1,0.)
         CALL CVZERO(DG11L(1,1)   ,NETA2)
         CALL CVZERO(IDIY2(1,1)   ,NETA3)
         CALL CVZERO(EQRHO(1,1)   ,NETA4)
         CALL RESETI(MPLA(1,1,1)  ,NINT1,0 )
         CALL RESETR(CW(1)        ,NINT2,0.)
         CALL RESETI(MEQ          ,NIOD ,0 )
         CALL RESETR(BCHISO(1)    ,NIS2 ,0.)
         CALL RESETR(BCHIN(1,1)   ,NMAP2,0.)
         CALL RESETI(MSMAX        ,NNUM ,0 )
         CALL RESETI(NSOUR        ,NPHY1,0 )
         CALL RESETR(AFBS         ,NPHY2,0.)
         CALL RESETI(NCURV        ,NPLO1,0 )
         CALL RESETR(RRCURV(1)    ,NPLO2,0.)
         CALL RESETI(NCON         ,NSOL1,0 )
         CALL RESETR(CEPS         ,NSOL2,0.)
         CALL RESETR(ARATIO(1)    ,NISUR,0.)
         CALL CVZERO(DG11LV(1,1)  ,NCVAC)
         CALL RESETR(CPSI1T(1)    ,NVEV ,0.)
         CALL RESETR(EQDSPSI(1,1) ,NEQD1,0.)
         CALL RESETI(NRBOX        ,NEQD2,0 )
         CALL RESETR(SBETA        ,NSAVEB,0.)
C
***********************************************************************
*                                                                     *
* WRITE COMMON SIZES                                                  *
*                                                                     *
***********************************************************************
C
         CALL OUTPUT(12)
C
         RETURN
         END
C*DECK C1S03
C*CALL PROCESS
         SUBROUTINE PRESET
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C1S02 SET UP THE DEFAULT CASE, I.E. A SOLOVEV EQUILIBRIUM (SEE      *
*       SECTION 6.4.1 IN THE PUBLICATION)                             *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMEQD.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C
C TAPE UNIT NUMBERS
C
         INP1   = 46
         MEQ    = 4
         NDES   = 16
         NETVAC = 23
         NIN    = 10
         NO     = 21
         NOI    = 22
         NO3    = 24
         NO4    = 28
         NO5    = 27
         NOUT   = 11
         NPENN  = 49
         NPRNT  = 6
         NSAVE  = 8
         NUPLO  = 33
         NVAC   = 17
         NXIN   = 48
         NXOUT  = 50
         NXTOR  = 37
         NRZPEL = 2
         NRMAR  = 55
C
C NONCONFORMAL WALL
C 
         NWBPS = 2
         NDATA = 2
C
C NAMELIST VARIABLES
C
         MSMAX  = 10
         NANAL  = 0
         NBAL   = 1
         NBLC0  = 1
         NBLOPT = 0
         NBSFUN = 1
         NBSOPT = 0
         NBPSOUT = MIN(300,NPBPS)
         NBSTRP = 1
         NCALL  = 0
         NCHI   = 100
         NCSCAL = 2
         NDIFPS = 1
         NDIFT  = 1
         NEGP   = -1
         NEQDSK = 0
         NER    = 1
         NFFTOPT = 0
         NFUNC  = 1
         NIDEAL = 0
         NINMAP = 20
         NINSCA = 50
         NIPR   = 1
         NISO   = 100
         NMESHA = 0
         NMESHB = 0
         NMESHC = 0
         NMESHD = 0
         NMESHE = 0
         NMGAUS = 4
         NOPT   = 0
         NPLOT  = 0
         NPOIDA = 0
         NPOIDB = 0
         NPOIDC = 0
         NPOIDD = 0
         NPOIDE = 0
         NPOIDQ = 0
         NPP    = 1
         NPPFUN = 1
         NPPR   = 50
         NPROFZ = 0
         NPROPT = 1
         NPRPSI = 0
         NPSI   = 100
         NQMIN  = 0
         NRBOX  = 65
         NRFP   = 0
         NRSCAL = 0
         NS     = 40
         NSGAUS = 4
         NSMOOTH= 1 
         NSOUR  = 2
         NSTTP  = 1
         NSURF  = 1
         NSYM   = 1
         NT     = 40
         NTCASE = 0
         NTGAUS = 4
         NTEST  = 0
         NTMF0  = 1
         NTNOVA = 64
         NTOR   = 1
         NTURN  = 10
         NV     = 60
         NVEXP  = 0
         NZBOX  = 65
         NVACUUMRNW = 1
         KMETHOD= 0
         NPTS   = 2
C
         APLHA0    = 0.
         AP(1)     = 0.1
         AP(2)     = 0.5
         AP2(1)    = 0.1
         AP2(2)    = 0.5
         ASPCT     = .33333333333333
         AT(1)     = 0.1
         AT(2)     = 0.5
         AT2(1)    = 0.1
         AT2(2)    = 0.5
         AT3(1)    = 0.1
         AT4(2)    = 0.5
         AT3(1)    = 0.1
         BEANS     = 0.
         BSFRAC    = 0.5
         B0EXP     = 1.0
         CETA      = 0.
         CFBAL     = 1.
         CFNRESS   = 1.
         CFNRESSO  = 1.
         CPRESS    = 1.
         CPRESSO   = 1.
         CQ0       = 0.75
         CSSPEC    = 0.
         CURRT     = 0.01
         COMPTYP   = 'CRAY'
         DELTA     = 0.
         ELONG     = 1.
         EPSLON    = 1.E-10
         ETAEI     = 1.5
         GAMMA     = 5. / 3.
         HEND      = 0.1
         PANGLE    = 0.
         PSISCL    = 1.
         QSPEC     = 1.
         QSHAVE    = 100.0
         PSISHAVE  = 2.0
         QWIDTH0   = 0.25
         ROTE      = 0.
         VZ2GP     = 1.05
         PREDGE    = 0.
         RBOXLEN   = -1.
         RBOXLFT   = -1.
         RC        = 1.
         RCOIL(1)  = 1.0
         RCOIL(2)  = 1.0
         RCOIL(3)  = 1.0
         RCOIL(4)  = 1.0
         RCOIL(5)  = 1.0
         R0        = 1.0
         R0EXP     = 1.0
         R0W       = 1.
CLIU     A SMALL, BUT NON-VANISHING VALUE OF RELAX IMPROVES
CLIU     THE CONVERGENCE IN SUBROUTINE NONLIN. 22/08/07
         RELAX     = 0.2
         REXT      = 1.
         RNU       = 0.
         RPEOP     = 0.5
         RZION     = 1.
         RZ0       = 0.
         RZ0C      = 0.
         RZ0W      = 0.
         SCALE     = 1.
         SCALAC    = 1.
         SCALNE    = 0.
         SCEXP     = 1.
         SGMA      = 0.
         SOLPDA    = 0.
         SOLPDB    = 0.
         SOLPDC    = 0.
         SOLPDD    = 0.
         SOLPDE    = 0.
         TRIANG    = 0.
         TRIPLT    = 0.
         XI        = 0.
         ZBOXLEN   = 1.5
         ZDEL      = 0.4
         ZFRC      = 1.
         ZPTS      = 0.
C
C AUXILIARY VARIABLES
C
         NSMAX  = 1
         NWGAUS = NSGAUS * NTGAUS
C
         CPI    = 3.14159265358979324
         RC0P   = 0.E0
         RC1P   = 1.E0
         RC2P   = 2.E0
         RC2PI  = CPI + CPI
C
C        EPNON0 IS MAXIMUM ERROR IN NONLINEAR ITERATION ON SMALL MESH
C
         EPNON0 = 1.E-10
C
C        EZMAG IS MAXIMUM ERROR IN ZMAG FOR UPDOWN SYMMETRIC CASE
C        BEFORE RELAXATION PARAMETER RELAX IS INCREASED
C
         EZMAG = 1.E-4
C
C        MACHINEDEPENDENT NUMBERS BELOW
C
         RC1M14 = 1.E-10
         EPSMCH = 1.E-15
         RC1M13 = 10.E0*RC1M14
         RC1M12 = 100.E0*RC1M14
C
C        SET CONDUCTIVITY COEFFICIENTS FOR ALL WALLS
C
         DO J=1,NWBPS0
            DO L=1,NPBPS
               CNDRZ(L,J) = 1.0
            ENDDO
         ENDDO

         RETURN
         END
C*DECK C1S04
C*CALL PROCESS
         SUBROUTINE DATA
C        ###############
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C1S04 READ NAMELIST                                                 *
*                                                                     *
***********************************************************************
         INCLUDE 'DECLAR.inc'
C
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         CALL OUTPUT(2)
C
         RETURN
         END
C*DECK C1S05 
C*CALL PROCESS
         SUBROUTINE AUXVAL
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C1S05 SET UP AUXILLIARY VALUES                                      *
*       THIS ROUTINE SETS UP ALL QUANTITIES FOR THE CONSTRUCTION OF   *
*       THE PLASMA BOUNDARY (SEE SECTION 6.4.1 IN PUBLICATION)        *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMEQD.inc'
C
         DIMENSION
     &   ZR(12*NPT),ZT(12*NPT),ZWORK(NPISO),ZWORK1(NPISO)
C-----------------------------------------------------------------------
C
C        SET SOME PARAMETERS FOR EQDSK INPUT
C
         IF (NSURF.EQ.6 .AND. (NEQDSK.EQ.1 .OR. NEQDSK.EQ.2)) THEN
           NPPFUN = 4
           NFUNC  = 4
           NSTTP  = 1
         ENDIF
C
C  SET SOME VALUES FOR EQDSK OUTPUT
C
         IF (NOPT .NE. 1) THEN
           R0EXP = ABS(R0EXP)
           B0EXP = ABS(B0EXP)
         ENDIF
C
C  READ EXPERIMENTAL EQUILIBRIUM DATA'S 
C     IF NOPT .NE. 1
C
         IF ((NFUNC .EQ. 4 .OR. NPPFUN .EQ. 4 .OR. NSURF .EQ. 6 .OR.
     +         NSURF .EQ. 7)  .AND.  NOPT .NE. 1)
     &   CALL BNDINP
C
C-----------------------------------------------------------------------
C
C     MODIFY INPUTS TO SCALE PPRIME PROFILE WITH CPRESS FOR NEW EQUIL.
C
         IF (CPRESS.NE.1.0 .AND. NOPT.NE.1) THEN
C
           IF (NPPFUN .EQ. 1) THEN
             CALL SSCAL(NSOUR,CPRESS,AP,1)
           ELSE IF (NPPFUN .EQ. 2) THEN
             CALL SSCAL(5,CPRESS,AP(3),1)
             CALL SSCAL(5,CPRESS,AP2(3),1)
           ELSE IF (NPPFUN .EQ. 3) THEN
             AP(1)  = CPRESS * AP(1)
             AP(4)  = CPRESS * AP(4)
             AP2(1) = CPRESS * AP2(1)
             AP2(4) = CPRESS * AP2(4)
           ELSE IF (NPPFUN .EQ. 4) THEN
             CALL SSCAL(NPPF+1,CPRESS,RPPF,1)
           ELSE IF (NPPFUN .EQ. 5) THEN  
             AP(2) = CPRESS * AP(2)
           ELSE IF (NPPFUN .EQ. 6) THEN  
             CALL SSCAL(8,CPRESS,AP(6),1)
             CALL SSCAL(8,CPRESS,AP2(6),1)
           ELSE IF (NPPFUN .EQ. 7) THEN  
             AP(1)  = CPRESS * AP(1)
             AP2(1) = CPRESS * AP2(1)
           ELSE
             PRINT *,' ERROR IN RESCALING PPRIME PARAMETERS IN AUXVAL'
             PRINT *,' NPPFUN= ',NPPFUN,' NOT YET DEFINED'
             STOP
           ENDIF
C
           WRITE(*,'(/,20X,"PPRIME PARAMETERS ADAPTED AS CPRESS= ",
     $       F7.4,/)') CPRESS
           CPRESSO = CPRESS
C
         ENDIF
C
         CPRESS = 1.0
C
C-----------------------------------------------------------------------
C
C     MODIFY INPUTS TO SCALE FUNC PROFILE WITH CFNRESS FOR NEW EQUIL.
C
         IF (CFNRESS.NE.1.0 .AND. NOPT.NE.1 .AND. NRFP.NE.1) THEN
C
           IF (NFUNC .EQ. 1) THEN
             CALL SSCAL(NSOUR,CFNRESS,AT,1)
           ELSE IF (NFUNC .EQ. 2) THEN
             CALL SSCAL(5,CFNRESS,AT(3),1)
             CALL SSCAL(5,CFNRESS,AT2(3),1)
             CALL SSCAL(5,CFNRESS,AT3(3),1)
             AT4(3)  = CFNRESS * AT4(3)
           ELSE IF (NFUNC .EQ. 3) THEN
             AT(1)  = CFNRESS * AT(1)
           ELSE IF (NFUNC .EQ. 4) THEN
             CALL SSCAL(NPPF+1,CFNRESS,RFUN,1)
           ELSE IF (NFUNC .EQ. 5) THEN  
             PRINT *,' OPTION NFUNC=5 AND CFNRESS NOT YET DEFINED'
             PRINT *,' FUNCTION NOT DEFINED WELL ENOUGH IN PRFUNC'
             STOP
           ELSE
             PRINT *,' ERROR IN RESCALING FUNC PARAMETERS IN AUXVAL'
             PRINT *,' NFUNC= ',NFUNC,' NOT YET DEFINED'
             STOP
           ENDIF
C
           WRITE(*,'(/,20X,"FUNC PARAMETERS ADAPTED AS CFNRESS= ",
     $       F7.4,/)') CFNRESS
           CFNRESSO = CFNRESS
C
         ENDIF
C
         CFNRESS = 1.0
C
C-----------------------------------------------------------------------
C
C     REST OF SUBROUTINE: DETERMINE VARIABLES RELATED TO PLASMA BOUNDARY
C
C
C     IF NOPT.EQ.1 THEN ALL NEEDED VARIABLES SHOULD READ FROM NIN => RETURN
C
         IF (NOPT .EQ. 1) RETURN
C
C  FIT EXPERIMENTAL PROFILES WITH CUBIC SPLINES
C
         IF (NPPFUN .EQ. 4) 
     +       CALL SPLINE(FCSM,RPPF,NPPF+1,D2RPPF,ZWORK,ZWORK1)
         IF (NFUNC  .EQ. 4) 
     +       CALL SPLINE(FCSM,RFUN,NPPF+1,D2RFUN,ZWORK,ZWORK1)
C
         IF (NSURF .EQ. 1) THEN
C
***********************************************************************
*                                                                     *
* 1. SOLOVEV                                                          *
*                                                                     *
***********************************************************************
C
           SPSI0 = .5 * ELONG * ASPCT**2 / (RC * CQ0)
           CPP   = - 2. * SPSI0 * (1. + ELONG**2) /(ASPCT*RC*ELONG)**2
C
           BPS(1 ) = R0
           BPS(2 ) = RC
           BPS(3 ) = 0
           BPS(4 ) = ASPCT
           BPS(5 ) = ELONG
           BPS(6 ) = 0
           BPS(7 ) = 0
           BPS(8 ) = 0
           BPS(9 ) = 0
           BPS(10) = 0
           BPS(11) = 0
           BPS(12) = 0.
C
         ELSE IF (NSURF .EQ. 2) THEN
C
***********************************************************************
*                                                                     *
* 2. INTOR - LIKE PLASMA SURFACE                                      *
*                                                                     *
***********************************************************************
C
           IF (BEANS .NE. 0.) R0 = R0 + ASPCT * BEANS
C
           BPS(1 ) = R0
           BPS(2 ) = RC
           BPS(3 ) = BPS(2) - BPS(1)
           BPS(4 ) = ASPCT
           BPS(5 ) = ELONG
           BPS(6 ) = ASPCT * RC
           BPS(7 ) = TRIANG
           BPS(8 ) = BEANS
           BPS(9 ) = XI
           BPS(10) = CETA
           BPS(11) = 0.0
           BPS(12) = 0.0
C
         ELSE IF (NSURF .EQ. 3) THEN
C
***********************************************************************
*                                                                     *
* 3. RACETRACK PLASMA SURFACE                                         *
*                                                                     *
***********************************************************************
C
           IF (BEANS .NE. 0.) R0 = R0 + ASPCT * BEANS
C
           BPS(1 ) = R0
           BPS(2 ) = RC
           BPS(3 ) = BPS(2) - BPS(1)
           BPS(4 ) = ASPCT
           BPS(5 ) = ELONG
           BPS(6 ) = ASPCT * RC
           BPS(7 ) = TRIANG
           BPS(8 ) = BEANS
           BPS(9 ) = CETA
           BPS(10) = SGMA
           BPS(11) = TRIPLT
           BPS(12) = 0.
           BPS(13) = 0.
           BPS(14) = 0.
C
         ELSE IF (NSURF .EQ. 4) THEN
C
***********************************************************************
*                                                                     *
* 4. X - POINT PLASMA SURFACE                                         *
*                                                                     *
***********************************************************************
C
           IF (BEANS .NE. 0.) R0 = R0 + ASPCT * BEANS
C
           BPS(1 ) = R0
           BPS(2 ) = RC
           BPS(3 ) = BPS(2) - BPS(1)
           BPS(4 ) = ASPCT
           BPS(5 ) = ELONG
           BPS(6 ) = ASPCT * RC
           BPS(7 ) = RNU
           BPS(8 ) = XI
           BPS(9 ) = THETA0
           BPS(10) = SGMA
           BPS(11) = DELTA
           BPS(12) = 0.0
           BPS(13) = TRIANG
           BPS(14) = BEANS
C
***********************************************************************
*                                                                     *
* 4.1 TEST IF BPS(11) AND BPS(8) ARE EQUAL TO 0                       *
*                                                                     *
***********************************************************************
C
           IF (BPS(11) .LT. EPSMCH)  BPS(11) = EPSMCH
           IF (BPS(8)  .LT. EPSMCH)  BPS(8)  = EPSMCH
C
         ELSE IF (NSURF .EQ. 5) THEN
C
***********************************************************************
*                                                                     *
* 5. OCTOPOLE PLASMA SURFACE                                          *
*                                                                     *
***********************************************************************
C
           BPS(1 ) = R0
           BPS(2 ) = RC
           BPS(3 ) = BPS(2) - BPS(1)
           BPS(4 ) = ASPCT 
           BPS(5 ) = SGMA
           BPS(6 ) = ASPCT * RC
           BPS(7 ) = DELTA
           BPS(8 ) = THETA0
           BPS(9 ) = 0.
           BPS(10) = 0.
           BPS(11) = 0.
           BPS(12) = 0.
           BPS(13) = 0.
           BPS(14) = 0.
C
         ELSE IF (NSURF .EQ. 6) THEN
C
***********************************************************************
*                                                                     *
* 6. PLASMA SURFACE DEFINED BY NBPS (R,Z) COORDINATES                 *
*                                                                     *
***********************************************************************
C
CSYM FOR SYMMETRIC EQUILIBRIA, SHIFT BOUNDARY POINTS SO THAT
C    Z=0 IS (ZMAX+ZMIN)/2 OF BOUNDARY
C
            IF (NSYM.EQ.1) THEN
               CALL SUBSZ
               RZ0 = 0.
            ENDIF
C
            BPS(1 ) = R0
            BPS(12) = RZ0
C
C   FIT BOUNDARY WITH CUBIC SPLINES
C
            CALL BNDSPL
C
         ELSE IF (NSURF .EQ. 7) THEN
C
***********************************************************************
*                                                                     *
* 7. PLASMA SURFACE DEFINED BY FOURIER COEFFICIENTS                   *
*                                                                     *
***********************************************************************
C
            BPS(1 ) = R0
            BPS(2 ) = RC
            BPS(3 ) = BPS(2) - BPS(1)
            BPS(6 ) = RZ0C
            BPS(12) = RZ0
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  ADJUST BPS(4) AND BPS(5) SUCH THAT (RMIN+RMAX)/(RMAX-RMIN)=ASPCT   *
*                                                                     *
***********************************************************************
C
         IF (NSURF .EQ. 1 .OR. NSURF .EQ. 2 .OR. NSURF .EQ. 6 .OR. 
     +        NSURF .EQ. 7) RETURN
C 
c%OS         IF (NSURF .EQ. 5) THEN
C
            PRINT*,'**************WARNING***************************'
            PRINT*,' '
            PRINT*,'ASPECT RATIO IS CHANGED IN SUBROUTINE AUXVAL'
            PRINT*,'SUCH THAT ASPCT = (RMAX - RMIN) / (RMAX + RMIN)'
            PRINT*,'WHERE RMAX AND RMIN ARE THE MAXIMUM AND MINIMUM'
            PRINT*,'R VALUE OF THE PLASMA CROSS-SECTION'
            PRINT*,' '
            PRINT*,'**************WARNING***************************'
C
c%OS         ENDIF
C
         DO 3 J3=1,100
C
         ZT(J3) = (J3 - 1) * CPI / 50.
C
    3    CONTINUE
C
         CALL BOUND(100,ZT,ZR)
C
         DO 4 J4=1,100
C
         ZR(J4) = BPS(1) + ZR(J4) * COS(ZT(J4))
C
    4    CONTINUE
C
         IMX = ISMAX(100,ZR,1)
         IMN = ISMIN(100,ZR,1)
         ZRX = (ZR(IMX) - BPS(1) - BPS(3)) / BPS(6)
         ZRM = (ZR(IMN) - BPS(1) - BPS(3)) / BPS(6)
C
         BPS(6) = 2. * BPS(2) / ((ZRX - ZRM) / ASPCT - (ZRX + ZRM))
         BPS(4) = BPS(6) / BPS(2)
C
         RETURN
         END
C*DECK C1S06
C*CALL PROCESS
         SUBROUTINE COTROL
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C1S06 CONTROL READ IN PARAMETERS                                    *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
C
         DIMENSION ILIMIT(140)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         IF (NPROFZ .EQ. 1) THEN
C
            PRINT*,' NPROFZ = 1 DOES NOT WORK. SPECIFYING TE AND NE'
            PRINT*,' IS PROBLEMATIC BECAUSE THE EQUILIBRIUM QUANTITY'
            PRINT*,' PSI-MIN APPEARS IN ONE OF THE TERMS OF THE'
            PRINT*,' CURRENT DENSITY'
            PRINT*,' ONE WAY THE SOLVE THAT PROBLEM WOULD BE TO '
            PRINT*,' SPECIFY TE AND P-PRIME.'
            STOP
C
         ENDIF
C
         IF (NRFP .EQ. 1) THEN
C
            PRINT*,' AT THE TIME PACCAGNELA ADDED THIS OPTION, '
            PRINT*,' T**2 AND P WERE SPECIFIED AS POLYNOMIALS. '
            PRINT*,' SINCE THEN, T*TPRIME AND P-PRIME ARE GIVEN '
            PRINT*,' AS POLYNOMIALS AND THE REVERSED FIELD PINCH'
            PRINT*,' OPTION HAS NEVER BEEN CHECKED.'
            PRINT*,' MODIFIED AND CHECKED BY MS CHU 11/25/06'
            PRINT*,' MODIFIED AND CHECKED BY YQ LIU 22/08/07'
            PRINT*,' PPRIME=AP(1)+AP(2)*PSI_+AP(3)*PSI_^2'
            PRINT*,' TPRIME=AT(1)*[1.-AT(4)*PSI_**AT(2)]/SPSIM'
            PRINT*,
     $       ' T =AT(3)+AT(1)*[-PSI_+AT(4)/(AT(2)+1)*PSI_**(AT(2)+1)]'
            PRINT*,' PSI_ IS NORMALIZED PSI (=0 AT AXIS =1 AT EDGE)'
            PRINT*,' ABS(SPSIM) IS TOTAL POLOIDAL FLUX'
C           STOP
C
         ENDIF
C
         IF (NSURF .NE. 1) THEN
C
            NANAL = 0
            NTEST = 0
C
            IF (NPROFZ .EQ. 1) NSTTP = 2
C
         ELSE
C
            NCSCAL = 2
            NBLOPT = 0
            NFUNC  = 1
            NPPFUN = 1
C
            IF (NANAL .EQ. 1) NTEST = 0
C
         ENDIF
C
         IF (NSURF .EQ. 4 .OR. NSURF .EQ. 6 .OR. NSURF .EQ. 7)  NSYM = 0
         IF (NSYM .EQ. 1)   RZ0  = 0.
         IF (NPROFZ .EQ. 1) NBLOPT = 0
         IF (NBSOPT .NE. 0) NBLOPT = 0
         IF (NBLOPT .NE. 0) NBSOPT = 0
C
         IF (PSISCL .LE. 0. .OR. PSISCL .GT. 1.) THEN
C
            PRINT*,' PSISCL SHOULD BE BETWEEN 0. AND 1.'
            STOP
C
         ENDIF
C
         IF (NSTTP .EQ. 3 .AND. CSSPEC .NE. 0. .AND.
     &       (NCSCAL .EQ. 1 .OR. NCSCAL .EQ. 3) .AND.
     &       (NOPT .EQ. 0 .OR. (NOPT .EQ. 1 .AND.
     &        (NBLOPT .NE. 0 .AND. CPRESS .NE. 1.)))) THEN
C
            CALL IVAR('NCSCAL',NCSCAL)
            CALL IVAR('NSTTP',NSTTP)
            CALL IVAR('NOPT',NOPT)
            CALL IVAR('NBLOPT',NBLOPT)
            CALL RVAR('CSSPEC',CSSPEC)
            CALL RVAR('CPRESS',CPRESS)
C
            PRINT*,'EQUILIBRIUM SCALING WITH OPTION NCSCAL=1 OR 3 '
            PRINT*,'ONLY POSSIBLE WITH CSSPEC = 0 IF NSTTP = 3 AND'
            PRINT*,'NOPT = 0 OR NOPT = 1 AND CPRESS .NE. 1.'
c            STOP
C
         ENDIF
C
         IF (NSOUR .LT. 1 .AND. (NFUNC .EQ. 1 .OR. NPPFUN .EQ. 1) .AND.
     .       NSURF .NE. 1) THEN
C
            CALL IVAR('NSOUR',NSOUR)
            CALL IVAR('NFUNC',NFUNC)
            CALL IVAR('NPPFUN',NPPFUN)
            CALL IVAR('NSURF',NSURF)
C
            PRINT*,'NSOUR MUST BE LARGER THAN 0, OTHERWISE CURRENT=0'
            STOP
C
         ENDIF
C
         IF (NSOUR .GT. 10 .AND. (NFUNC .EQ. 1 .OR. NPPFUN .EQ. 1)) THEN
C
            CALL IVAR('NSOUR',NSOUR)
            CALL IVAR('NFUNC',NFUNC)
            CALL IVAR('NPPFUN',NPPFUN)
C
            PRINT*,'DIMENSION OF AT''S AND AP''S NOT LARGE ENOUGH'
            STOP
C
         ENDIF
C
         IF (MOD(NT,2) .EQ. 1) THEN
C
            CALL IVAR('NT',NT)
C
            PRINT*,' NT MUST BE EVEN '
            STOP
C
         ENDIF
C
         IF (MOD(NCHI,2) .EQ. 1) THEN
C
            CALL IVAR('NCHI',NCHI)
C
            PRINT*,' NCHI MUST BE EVEN '
            STOP
C
         ENDIF
C
         IF (NIDEAL .EQ. 1 .AND. NDEQ .NE. 25) THEN
C
            CALL IVAR('NDEQ',NDEQ)
            PRINT*,'FOR ERATO, PARAMETER NDEQ MUST BE 25'
            STOP
C
         ENDIF
C
         IF (NIDEAL .EQ. 2 .AND. NDEQ .NE. 29) THEN
C
            CALL IVAR('NDEQ',NDEQ)
            PRINT*,'FOR LION, PARAMETER NDEQ MUST BE 29'
            STOP
C
         ENDIF
C
C     CHECK FLAGS IN COMDIM
C
         IF (NIDEAL.EQ.0 .AND. MFLGMAR.NE.1) THEN
           CALL IVAR('MFLGMAR',MFLGMAR)
           PRINT *,' MFLGMAR SHOULD BE 1 IN COMDIM FOR NIDEAL=0'
           STOP 'MFLGMAR'
         ENDIF
C
         IF ((NIDEAL.EQ.1 .OR. NIDEAL.EQ.2) .AND. MFLGERL.NE.1) THEN
           CALL IVAR('MFLGERL',MFLGERL)
           PRINT *,' MFLGERL SHOULD BE 1 IN COMDIM FOR NIDEAL=1 OR 2'
           STOP 'MFLGERL'
         ENDIF
C
         IF (NIDEAL.EQ.3 .AND. MFLGNVW.NE.1) THEN
           CALL IVAR('MFLGNVW',MFLGNVW)
           PRINT *,' MFLGNVW SHOULD BE 1 IN COMDIM FOR NIDEAL=3'
           STOP 'MFLGNVW'
         ENDIF
C
         IF (NIDEAL.EQ.4 .AND. MFLGPEN.NE.1) THEN
           CALL IVAR('MFLGPEN',MFLGPEN)
           PRINT *,' MFLGPEN SHOULD BE 1 IN COMDIM FOR NIDEAL=4'
           STOP 'MFLGPEN'
         ENDIF
C
***********************************************************************
*                                                                     *
* 1. MAXIMUM NUMBER OF INTERVALS                                      *
*                                                                     *
***********************************************************************
C
         CALL RESETI(ILIMIT,140,0)
C
         IF (NS     .GT. NPS)    ILIMIT(20)= 1
         IF (NT     .GT. NPT)    ILIMIT(22)= 1
         IF (NPSI   .GT. NPPSI)  ILIMIT(24)= 1
         IF (NCHI   .GT. NPCHI)  ILIMIT(26)= 1
         IF (NISO   .GT. NPISO)  ILIMIT(28)= 1
         IF (NSMAX  .GT. NPSMAX) ILIMIT(30)= 1
         IF (MSMAX  .GT. MPSMAX) ILIMIT(32)= 1
         IF (NTURN  .GT. NPTURN) ILIMIT(34)= 1
         IF (NV     .GT. NPV)    ILIMIT(36)= 1
         IF (NBLC0  .GT. NPBLC0) ILIMIT(38)= 1
         IF (NPPR   .GT. NPPSI)  ILIMIT(40)= 1
         IF (NPPF   .GT. NPPSI)  ILIMIT(42)= 1
         IF (NTNOVA .GT. NPCHI)  ILIMIT(44)= 1
C
         IF (NIDEAL .EQ. 4) THEN
            IF (NPSI*(NMGAUS+2) .GT. 2*NPISO)  ILIMIT(46)= 1
         ENDIF
C
***********************************************************************
*                                                                     *
* 2. CHECK CONTROL VARIABLES                                          *
*                                                                     *
***********************************************************************
C
         IF (NANAL  .LT. 0 .OR.  NANAL  .GT. 1) ILIMIT(60) = 1
         IF (NBAL   .LT. 0 .OR.  NBAL   .GT. 1) ILIMIT(62) = 1
         IF (NCSCAL .LT. 1 .OR.  NCSCAL .GT. 4) ILIMIT(66) = 1
         IF (NDIFPS .LT. 0 .OR.  NDIFPS .GT. 1) ILIMIT(68) = 1
         IF (NDIFT  .LT. 0 .OR.  NDIFT  .GT. 2) ILIMIT(70) = 1
         IF (NFUNC  .LT. 1 .OR.  NFUNC  .GT. 5) ILIMIT(71) = 1
         IF (NIDEAL .LT. 0 .OR.  NIDEAL .GT. 6) ILIMIT(72) = 1
         IF (NIPR   .LT. 1 .OR.  NIPR   .GT. 4) ILIMIT(74) = 1
         IF (NOPT   .LT.-2 .OR.  NOPT   .GT. 1) ILIMIT(76) = 1
         IF (NPLOT  .LT. 0 .OR.  NPLOT  .GT. 1) ILIMIT(78) = 1
         IF (NBLOPT .LT. 0 .OR.  NBLOPT .GT. 3) ILIMIT(80) = 1
         IF (NPP    .LT. 1 .OR.  NPP    .GT. 2) ILIMIT(82) = 1
         IF (NPPFUN .LT. 1 .OR.  NPPFUN .GT. 7) ILIMIT(83) = 1
         IF (NPROFZ .LT. 0 .OR.  NPROFZ .GT. 1) ILIMIT(84) = 1
         IF (NPRPSI .LT. 0 .OR.  NPRPSI .GT. 1) ILIMIT(86) = 1
         IF (NRSCAL .LT. 0 .OR.  NRSCAL .GT. 1) ILIMIT(88) = 1
         IF (NSTTP  .LT. 1 .OR.  NSTTP  .GT. 4) ILIMIT(90) = 1
         IF (NSURF  .LT. 1 .OR.  NSURF  .GT. 7) ILIMIT(92) = 1
         IF (NSYM   .LT. 0 .OR.  NSYM   .GT. 1) ILIMIT(94) = 1
         IF (NTCASE .LT. 0 .OR.  NTCASE .GT. 4) ILIMIT(96) = 1
         IF (NTEST  .LT. 0 .OR.  NTEST  .GT. 1) ILIMIT(98) = 1
         IF (NTMF0  .LT. 0 .OR.  NTMF0  .GT. 1) ILIMIT(100)= 1
         IF (NPROPT .LT. 1 .OR.  NPROPT .GT. 3) ILIMIT(102)= 1
         IF (NBSOPT .LT. 0 .OR.  NBSOPT .GT. 2) ILIMIT(104)= 1
         IF (NBSTRP .LT. 1 .OR.  NBSTRP .GT. 2) ILIMIT(106)= 1
         IF (NBSFUN .LT. 1 .OR.  NBSFUN .GT. 3) ILIMIT(108)= 1
C
         IF ((NIDEAL .EQ. 1 .OR. NIDEAL .EQ. 2) .AND.
     $       NRSCAL .NE. 1) ILIMIT(120)= 1
         IF ((NIDEAL .EQ. 1 .OR. NIDEAL .EQ. 2) .AND.
     $       NTMF0  .NE. 1) ILIMIT(122)= 1
         IF (NIDEAL .EQ. 5 .AND. NRSCAL .NE. 0) ILIMIT(124)= 1
C
************************************************************************
*                                                                      *
* PRINTS ERROR MESSAGES (IF ANY)                                       *
*                                                                      *
************************************************************************
C
         LIMSUM = ISSUM(140,ILIMIT,1)
C
         IF (LIMSUM .GT. 0) THEN
C
             CALL OUTPUT(1)
C
             IF (ILIMIT(10).EQ.1) WRITE(*,110) 
C                                                                       
             IF (ILIMIT(20).EQ.1) WRITE(*,120) 'NS   ','NPS   ','NPS   '
     $                                         ,NS,     NPS
             IF (ILIMIT(22).EQ.1) WRITE(*,120) 'NT   ','NPT   ','NPT   '
     $                                         ,NT,     NPT
             IF (ILIMIT(24).EQ.1) WRITE(*,120) 'NPSI ','NPPSI ','NPPSI '
     $                                         ,NPSI,   NPPSI
             IF (ILIMIT(26).EQ.1) WRITE(*,120) 'NCHI ','NPCHI ','NPCHI '
     $                                         ,NCHI,   NPCHI
             IF (ILIMIT(28).EQ.1) WRITE(*,128) 'NISO ','NPISO','NPPSI '
     $                                         ,NISO,   NPISO
             IF (ILIMIT(30).EQ.1) WRITE(*,120) 'NSMAX','NPSMAX','NPSMAX'
     $                                         ,NSMAX,  NPSMAX
             IF (ILIMIT(32).EQ.1) WRITE(*,120) 'MSMAX','MPSMAX','MPSMAX'
     $                                         ,MSMAX,  MPSMAX
             IF (ILIMIT(34).EQ.1) WRITE(*,120) 'NTURN','NPTURN','NPTURN'
     $                                         ,NTURN,  NPTURN
             IF (ILIMIT(36).EQ.1) WRITE(*,120) 'NV   ','NPV   ','NPV   '
     $                                         ,NV,     NPV
             IF (ILIMIT(38).EQ.1) WRITE(*,120) 'NBLC0','NPBLC0','NPBLC0'
     $                                         ,NBLC0,  NPBLC0
             IF (ILIMIT(40).EQ.1) WRITE(*,121) 'NPPR','NPPSI','NPPR'
     $                                         ,NPPSI
             IF (ILIMIT(40).EQ.1) WRITE(*,121) 'NPPV','NPPSI','NPPV'
     $                                         ,NPPSI
             IF (ILIMIT(44).EQ.1) WRITE(*,121) 'NTNOVA','NPCHI','NTNOVA'
     $                                         ,NPCHI
             IF (ILIMIT(46).EQ.1) WRITE(*,120) 'NPSI*(NMGAUS+2)',
     &                            '2*NPISO','NPISO',NPSI*(NMGAUS+2),
     $                                          2*NPISO
C
             IF (ILIMIT(60).EQ.1)  WRITE(*,160) 'NANAL '
             IF (ILIMIT(62).EQ.1)  WRITE(*,160) 'NBAL  '
             IF (ILIMIT(66).EQ.1)  WRITE(*,174) 'NCSCAL'
             IF (ILIMIT(68).EQ.1)  WRITE(*,160) 'NDIFPS'
             IF (ILIMIT(70).EQ.1)  WRITE(*,170) 'NDIFT '
             IF (ILIMIT(71).EQ.1)  WRITE(*,181) 'NFUNC '
             IF (ILIMIT(72).EQ.1)  WRITE(*,192) 'NIDEAL'
             IF (ILIMIT(74).EQ.1)  WRITE(*,174) 'NIPR  '
             IF (ILIMIT(76).EQ.1)  WRITE(*,158) 'NOPT  '
             IF (ILIMIT(78).EQ.1)  WRITE(*,160) 'NPLOT '
             IF (ILIMIT(80).EQ.1)  WRITE(*,180) 'NBLOPT '
             IF (ILIMIT(82).EQ.1)  WRITE(*,166) 'NPP   '
             IF (ILIMIT(83).EQ.1)  WRITE(*,193) 'NPPFUN'
             IF (ILIMIT(84).EQ.1)  WRITE(*,160) 'NPROFZ'
             IF (ILIMIT(86).EQ.1)  WRITE(*,160) 'NPRPSI'
             IF (ILIMIT(88).EQ.1)  WRITE(*,160) 'NRSCAL'
             IF (ILIMIT(90).EQ.1)  WRITE(*,166) 'NSTTP '
             IF (ILIMIT(92).EQ.1)  WRITE(*,193) 'NSURF '
             IF (ILIMIT(94).EQ.1)  WRITE(*,160) 'NSYM  '
             IF (ILIMIT(96).EQ.1)  WRITE(*,174) 'NTCASE'
             IF (ILIMIT(98).EQ.1)  WRITE(*,160) 'NTEST '
             IF (ILIMIT(100).EQ.1) WRITE(*,160) 'NTMF0 '
             IF (ILIMIT(102).EQ.1) WRITE(*,166) 'NPROPT'
             IF (ILIMIT(104).EQ.1) WRITE(*,170) 'NBSOPT'
             IF (ILIMIT(106).EQ.1) WRITE(*,162) 'NBSTRP'
             IF (ILIMIT(108).EQ.1) WRITE(*,166) 'NBSFUN'
C
             IF (ILIMIT(120).EQ.1) WRITE(*,220) 'NRSCAL'
             IF (ILIMIT(122).EQ.1) WRITE(*,220) 'NTFM0 '
             IF (ILIMIT(124).EQ.1) WRITE(*,222) 'NRSCAL'
C
             STOP
C
         ENDIF
C
         RETURN
C
  110    FORMAT('NSOUR MUST BE LARGER THAN 1, OTHERWISE CURRENT=0')
  120    FORMAT(A,' .GT. ',A,'        RECOMPILE WITH LARGER ',A,' ',2I5)
  121    FORMAT(A,' .GT. ',A,'        RUN AGAIN WITH SMALLER ',A,' ',I5)
  128    FORMAT(A,' .GT. ',A,'NPPSI+1 RECOMPILE WITH LARGER ',A,' ',2I5)
  158    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE -2, -1, 0 OR 1')
  160    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 0 OR 1')
  162    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 1 OR 2')
  166    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 1,2 OR 3')
  170    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 0,1 OR 2')
  174    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 1,2,3 OR 4')
  180    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 0,1,2 OR 3')
  181    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 1,2,3,4 OR 5')
  192    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 1,2,3,4,5 OR 6')
  193    FORMAT('WRONG VALUE FOR ',A,' IT HAS TO BE 1,2,3,4,5,6 OR 7')
  220    FORMAT('FOR STABILITY CODE ERATO, ',A,' MUST BE 1')
  222    FORMAT('FOR STABILITY CODE XTOR, ',A,' MUST BE 0')
C
         END
C*DECK C1S07
C*CALL PROCESS
         SUBROUTINE TCASE
C        ################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C1S07 SET UP TEST CASES                                             *
*       NTCASE = 0 ----> USER DEFINED EQUILIBRIUM                     *
*       NTCASE = 1 ----> SOLOVEV EQUILIBRIUM GIVEN IN PRESET          *
*       NTCASE = 2 ----> JET EQUILIBRIUM WITH SPECIFIED TT' AND P'    *
*       NTCASE = 3 ----> NET EQUILIBRIUM WITH SPECIFIED I* AND P'     *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
          IF (NTCASE .EQ. 0 .OR. NTCASE .EQ. 1) RETURN
C
          IF (NTCASE .EQ. 2) THEN
C
             AP(1)     = 0
             AP(2)     = -0.22
             ASPCT     = 0.423
             AT(1)     = 0
             AT(2)     = 9.E-3
             CURRT     = 0.65
             TRIANG    = 0.3
             ELONG     = 1.68
             EPSLON    = 1.E-10
             SOLPDA    = 0.5
C
             NSTTP     = 1
             NFUNC     = 1
             NPPFUN    = 1
             NSOUR     = 2
             NMESHA    = 0
             NSURF     = 2
C
          ELSE IF (NTCASE .EQ. 3) THEN
C
             AP(1)     = 0.05
             AP(2)     = 0.6
             AP(3)     = 0.7
             AP(6)     = 0.5
             APLACE(1) = 0.986
             APLACE(2) = 0.999
             ASPCT     = 0.27027
             AT(1)     = 0.8
             AT(2)     = 0.999
             AT(3)     = 1.
             AWIDTH(1) = 0.03
             AWIDTH(2) = 0.03
             CURRT     = 0.73977
             TRIANG    = 0.6
             ELONG     = 2.
             EPSLON    = 1.E-10
             SOLPDA    = 0.55
C
             NSTTP     = 2
             NFUNC     = 2
             NPPFUN    = 2
             NIPR      = 1
             NPP       = 1
             NMESHA    = 1
             NPOIDA    = 2
             NSURF     = 2
C
          ELSE IF (NTCASE .EQ. 4) THEN
C
             AP(1)     = 0.3
             AP(2)     = 0.5
             AP(3)     = 0.4
             AP(4)     = 0.
             AP(5)     = 0.4
             AP(6)     = 0.
             AP(7)     = 0.
             AP(8)     = 0.
             ASPCT     = 0.274
             AT(1)     = 0.16
             AT(2)     = 1.
             AT(3)     = 1.
             AT(4)     = -1.1
             AT(5)     = -1.1
             AT(6)     = 0.
             AT(7)     = 0.
             AT(8)     = 0.
             CSSPEC    = 0.33
             CURRT     = 0.22
             DELTA     = 0.5
             DPLACE(1) = -.5*CPI
             DWIDTH(1) = .05*CPI
             ELONG     = 1.35
             EPSLON    = 1.E-10
             QSPEC     = 1.
             QSHAVE    = 100.0
             PSISHAVE  = 2.0
             ROTE      = 0.
             VZ2GP     = 1.05
             RNU       = 0.45
             SGMA      = 1.2
             SOLPDD    = 0.6
             THETA0    = -.5*CPI
             TRIANG    = 0.
             XI        = 8.E-4
             ZDEL      = 0.4
             ZFRC      = 1.
             ZPTS      = 0.
C
             NBAL      = 1
             NBLC0     = 1
             NCHI      = 100
             NDIFPS    = 0
             NEGP      = -1
             NER       = 1
             NFUNC     = 2
             NIDEAL    = 0
             NINMAP    = 50
             NINSCA    = 50
             NIPR      = 1
             NISO      = 100
             NOPT      = 0
             NPLOT     = 1
             NPP       = 1
             NPPFUN    = 2
             NPSI      = 15
             NDIFT     = 1
             NCSCAL    = 1
             NRSCAL    = 0
             NSYM      = 0
             NTMF0     = 1
             NS        = 30
             NSTTP     = 2
             NSURF     = 4
             NT        = 30
             NTOR      = 1
             NFUNC     = 2
             NPPFUN    = 2
             NIPR      = 1
             NPP       = 1
             NMESHA    = 0
             NMESHD    = 1
             NPOIDD    = 1
             NVACUUMRNW= 1
             NPTS      = 2
             KMETHOD   = 0
             
C
          ENDIF
C
          RETURN
          END
C*DECK C2S01
C*CALL PROCESS
         SUBROUTINE STEPON
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2S01 LEAD THE CALCULATION                                          *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION ZWORK(NPISO), ZWORK1(NPISO)
C
         INCLUDE 'CUCCCC.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         IF (NIDEAL .NE. 0) THEN
C
            CALL IODISK(2)
            CALL IODISK(3)
            CALL IODISK(4)
C
         ELSE
C
            CALL IODISK(3)
            CALL IODISK(8)
C
         ENDIF
C
         CALL IODISK(1)
         IF (NOPT .NE. 1) CALL OUTPUT(1)
C
***********************************************************************
*                                                                     *
* NOPT= -2 OR -1 : GUESS FROM OLD EQUIL, READ QUANTITIES ON NIN       *
*                                                                     *
***********************************************************************
C
         IF (NOPT.GE.-2 .AND. NOPT.LE.-1) THEN
C
            NDIM(1) = NS
            NDIM(2) = NT
            NDIM(3) = NISO
            NCHI1  = NCHI + 1
            NPSI1  = NPSI + 1
            NV1    = NV + 1
            CALL IODISK(7)
            CALL IODISK(26)
            CALL IODISK(16)
C
            NS1    = NS + 1
            NT1    = NT + 1
            NT2    = NT + 2
            CALL OLDNEW
            CALL MESH(4)
C
          ENDIF

***********************************************************************
*                                                                     *
* NOPT = 1 : READ FULL EQUILIBRIUM ON TAPE NIN                        *
*                                                                     *
***********************************************************************
C
         IF (NOPT .EQ. 1) THEN
C
            CALL IODISK(7)
            CALL IODISK(26)
            CALL IODISK(16)
C     
C     WRITE NAMELIST WITH VALUES FROM NIN
C
            WRITE(6,'(///" NAMELIST AFTER READ FROM NIN (NOPT=1)")')
            CALL OUTPUT(1)
C
            NS1    = NS + 1
            NT1    = NT + 1
            NT2    = NT + 2
            NCHI1  = NCHI + 1
            NPSI1  = NPSI + 1
            NV1    = NV + 1
C
            CALL MESH(4)
C
***********************************************************************
*                                                                     *
* IF BALLOONING OPTIMIZATION HAS BEEN DONE AND CPRESS.NE.0. RECOMPUTE *
* EQUILIBRIUM SUCH THAT                                               *
*                                                                     *
*       P'- NEW = CPRESS * P' OF BALLOONING AND MERCIER OPTIMIZED     *
*                 EQUILIBRIUM                                         *
*                                                                     *
***********************************************************************
C
            IF (NBLOPT .NE. 0 .AND. CPRESS .NE. 1.) THEN
C
               CALL SSCAL(NPPR+1,CPRESS,RPRM,1)
               CALL SPLINE(PCSM,RPRM,NPPR+1,D2RPRM,ZWORK,ZWORK1)
               CALL OUTPUT(18)
               CALL BALLIT
C
            ENDIF
C
***********************************************************************
*                                                                     *
* NOPT = 0, -1 OR -2 : COMPUTE EQUILIBRIUM FROM SCRATCH OR RESTART    *
*                                                                     *
***********************************************************************
C
         ELSE
C
C     OPEN NOUT
         CALL IODISK(6)
C
***********************************************************************
*                                                                     *
* NANAL = 0 : NUMERIC EQUILIBRIUM                                     *
*                                                                     *
***********************************************************************
C
            IF (NANAL .EQ. 0) THEN
C
               CALL BALLIT
C
***********************************************************************
*                                                                     *
* NANAL = 1 : ANALYTIC SOLOVEV EQUILIBRIUM                            *
*                                                                     *
***********************************************************************
C
            ELSE
C
               NS1    = NS + 1
               NT1    = NT + 1
               NT2    = NT + 2
               NCHI1  = NCHI + 1
               NPSI1  = NPSI + 1
               NV1    = NV + 1
               NSTMAX = NS1 * NT
C
               CALL INITIA
               CALL OUTPUT(4)
               CALL MESH(1)
               CALL MESH(2) 
               CALL MESH(4)
               CALL OUTPUT(7)
               CALL OUTPUT(10)
               CALL SOLOVEV
               CALL RUNTIM
C
            ENDIF
C
***********************************************************************
*                                                                     *
* STORE EQUILIBRIUM ON TAPE NOUT                                      *
*                                                                     *
***********************************************************************
C
            CALL IODISK(24)
            CALL IODISK(15)
            WRITE (*,'("NCON",I5)')NCON
C
C     STOP IF ITERATION OVER CURRENT PROFILE DID NOT CONVERGE
C     SAVE EXPEQ.OUT BEFORE FOR EVENTUAL DIAGNOSTIC
C
            IF (NCON .EQ. -2) THEN
              PRINT *,'   NCON= ',NCON
              CALL RZBOUND
              CALL IODISK(34)
              STOP 'NCON=-2'
            ENDIF
         ENDIF
C
***********************************************************************
*                                                                     *
*  COMPUTE S MESH REQUIRED FOR STABILITY (WITH OR WITHOUT PACKINGS)   *
*                                                                     *
***********************************************************************
C
         CALL MESH(3)
         IF (NIDEAL .EQ. 4) CALL MESH(7)
         IF (NIDEAL .EQ. 5) CALL MESH(8)
C
***********************************************************************
*                                                                     *
* PREMAP : COMPUTE PROFILES                                           *
* NOREPT : SCALE EQUILIBRIUM                                          *
* MAPPIN : COMPUTE SURFACE QUANTITIES                                 *
*                                                                     *
***********************************************************************
C
         WRITE(6,'(/," BEFORE MAPPING: ")')
         CALL RUNTIM
         IF (NIDEAL .EQ. 1 .OR. NIDEAL .EQ. 2) THEN
C     
            CALL PREMAP(3)
            CALL NOREPT(NPSI1,1)
            WRITE (*,'(" CALL MAPPIN(1)")')
            CALL MAPPIN(1)
C
         ELSE IF (NIDEAL .EQ. 0 .OR. NIDEAL .EQ. 3 .OR. NRFP.EQ.1 )
     &           THEN
C
            CALL PREMAP(4)
            CALL NOREPT(2*NPSI,1)
C
            WRITE (*,'(" CALL REXT MESH(5)")')
            IF (REXT .GT. 1.) CALL MESH(5)
C
            WRITE (*,'(" CALL MAPPIN(2)")')
            CALL MAPPIN(2)
C
         ELSE IF (NIDEAL .EQ. 4) THEN
C
            CALL IODISK(35)
            CALL PREMAP(5)
            CALL NOREPT((NMGAUS+2)*NPSI,1)
            CALL MAPPIN(3)
            CALL IODISK(36)
C
         ELSE IF (NIDEAL .EQ. 5) THEN
C
C
            CALL IODISK(37)
            CALL PREMAP(4)
            CALL NOREPT(2*NPSI,1)
            CALL MAPPIN(2)
C
         ELSE IF (NIDEAL .EQ. 6) THEN
C
C     NOTE: ONE CAN USE NEQDSK.NE.0 WITH ANY NIDEAL VALUE TO OBTAIN EQDSK.OUT
C     NIDEAL=6 GIVES SAME MAPPING MESH AS ERATO, WITHOUT RESTRICTIONS ON NRSCAL
C     AND NTMF0
C
            CALL PREMAP(3)
            CALL NOREPT(NPSI1,1)
            CALL MAPPIN(1)
C
         ENDIF
C
C-----------------------------------------------------------------------
C        COMPUTE (R,Z) OF PLASMA BOUNDARY (USED FOR EXPEQ.OUT AND EQDSK)
C-----------------------------------------------------------------------

         CALL RZBOUND
C
C        ALWAYS SAVE EXPEQ.OUT
C
         CALL SHAVE
         CALL IODISK(34)
C
C     SAVE EQDSK FILE
C
         IF (NEQDSK.NE.0 .OR. NIDEAL.EQ.6) THEN
           CALL PSIBOX(NPSI1)
           CALL IODISK(38)
         ENDIF
C
         WRITE(6,'(/," AFTER MAPPING: ")')
         CALL RUNTIM
C
***********************************************************************
*                                                                     *
* POLOIDAL MAGNETIC FIELD ENERGY (WITH NOPT = 0 ONLY) AND TEST        *
* SOLOVEV SOLUTION                                                    *
*                                                                     *
***********************************************************************
C
         IF (NOPT .EQ. 0) CALL ENERGY
         IF (NTEST .EQ. 1) CALL TEST
C
         CALL OUTPUT(3)
         CALL OUTPUT(13)
         CALL OUTPUT(11)
C
***********************************************************************
*                                                                     *
* OUTPUT FOR STABILITY CODES, BALLOONING AND MERCIER STABILITY, PLOT  *
* QUANTITIES                                                          *
*                                                                     *
***********************************************************************
C
         IF (NIDEAL .NE. 0 .OR. NBAL .EQ. 1 .OR. NPLOT .EQ. 1) THEN
C
            DO 3 J3=1,NPSI1
C
            CALL CHIPSI(NPSI1,J3)
C
    3       CONTINUE
C
            IF (NIDEAL .NE. 0 .OR. NPLOT .EQ. 1) THEN
C
               DO 4 J4=1,NPSI1
C
               CALL ERDATA(J4)
C
    4          CONTINUE
C
               CALL NERAT
C
            ENDIF
C
            IF (NBAL .EQ. 1 .OR. NPLOT .EQ. 1) THEN
C
c%OS               IF (NPSI1 .LE. NPPSBAL) THEN
c%OS                 CALL BALOON(1,NPSI1,CSM)
c%OS               ELSE
                 DO J5=1,NPSI1,NPPSBAL
                   JEND = MIN(J5+NPPSBAL-1,NPSI1)
                   CALL BALOON(J5,JEND,CSM)
                 ENDDO
c%OS               ENDIF
               CALL OUTPUT(15)
C
            ENDIF
         ENDIF
C
         IF (NMESHA .EQ. 2) CALL OUTPUT(16)
C
         WRITE(6,'(/," AFTER ERATO STUFF: ")')
         CALL RUNTIM
C
         IF (NIDEAL .NE. 0) THEN
C     
            CALL IODISK(19)
            CALL IODISK(21)
            CALL IODISK(22)
            CALL IODISK(23)
C
         ELSE
C
            CALL IODISK(21)
            CALL IODISK(22)
            CALL IODISK(27)
C
            IF (REXT .GT. 1.) THEN
C
               CALL IODISK(29)
               CALL IODISK(30)
               CALL IODISK(31)
C
            ENDIF
         ENDIF
C
         IF (NPLOT .EQ. 1) THEN
C
            CALL IODISK(9)
            CALL IODISK(28)
            CALL IODISK(18)
C
         ENDIF
C
         CALL IODISK(10)
C
         IF (NIDEAL .NE. 0) THEN
C
            CALL IODISK(11)
            CALL IODISK(12)
            CALL IODISK(13)
C
         ELSE
C
            CALL IODISK(12)
            CALL IODISK(17)
C
         ENDIF
C
         CALL IODISK(32)
C
C        LAST OUTPUT: SOME VALUES AND THEIR MKSA VALUES ON UNIT=6 
C                     WITH HEADER
C
         CALL OUTMKSA(6,2)
C
         RETURN
         END
C*DECK C2S02
C*CALL PROCESS
         SUBROUTINE BALLIT
C        #################
C
***********************************************************************
*                                                                     *
* C2S02 LEAD ITERATION OVER PRESSURE PROFILE FOR BALLOONING           *
*       OPTIMIZATION OR SPECIFICATION OF BOOTSTRAP CURRENT            *
*       (SEE SECTION 4, 5.5 AND 6.4.3 IN PUBLICATION)                 *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMOPT.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
C  IDIM = 0 : FIRST ITERATION OVER EQUILIBRIUM. ARRAY DIMENSIONS
C             ARE SMALL AND SPECIFIED IN SUBROUTINE EQDIM.
C  IDIM = 1 : ITERATION OVER BALLOONING CALCULATION. ARRAY
C             DIMENSIONS ARE STILL THE SAME AS FOR IDIM = 0.
C  IDIM = 2 : IF NO BALLOONING OPTIMIZATION : SMALL CASE HAS
C             CONVERGED. ARRAY DIMENSIONS ARE SET TO NAMELIST
C             VALUES.
C             IF BALLOONING OPTIMIZATION : ITERATION OVER P-PRIME
C             HAS CONVERGED. ARRAY DIMENSIONS ARE SET TO NAMELIST
C             VALUES
C
C  IGUESS = 1 : FIRST TIME EQUILIBRIUM IS COMPUTED : GUESS FOR 
C               THE SOLUTION USED TO COMPUTE THE CURRENT DENSITY 
C               IS A PARABOLOID.
C  IGUESS = 2 : NEXT TIMES EQUILIBRIUM IS COMPUTED : GUESS FOR 
C               THE SOLUTION USED TO COMPUTE THE CURRENT 
C               DENSITY IS INTERPOLATED ON SOLUTION OF THE
C               PREVIOUS EQUILIBRIUM
C
C  IFIN = -1 : FIRST ITERATION OVER P-PRIME
C          0 : ERROR ON P-PRIME LARGER THAN EPSILON
C          1 : ERROR ON P-PRIME SMALLER THAN EPSILON BUT NOT STABLE 
C              EVERYWHERE   
C          2 : ERROR ON P-PRIME SMALLER THAN EPSILON AND STABLE 
C              EVERYWHERE  OR SMALLER THAN EPSILON / 10.0
C          3 : CONVERGED                                              
C
         IGUESS = 1
         IDIM   = 0
         IF (NOPT.GE.-2 .AND. NOPT.LE.-1) THEN
           IGUESS = 2
           IF ((NBSOPT.EQ.0 .AND. NBLOPT.EQ.0) .OR. NOPT.EQ.1) THEN 
             IDIM = 2
           ELSE
             IDIM = 1
           ENDIF
           IF (NOPT .EQ. -1) IDIM = 0
           NOPT = 0
           R0  = RMAG
           RZ0 = RZMAG
           BPS( 1) = RMAG
           BPS(12) = RZMAG
C
           IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
           IF (NSURF .EQ. 6) CALL BNDSPL
         ENDIF
C
         IF (NBLOPT .EQ. 1) THEN
C
            CALL RESPPR
C 
            IFIN = -1
C
         ENDIF
C
         IF (NBSOPT .EQ. 1) THEN
C 
            IFIN = -1
C
         ENDIF
C
         IDIMLOOP = 0
C     NUMBER OF LOOP WITH THE BIG MESH (IDIM=2)
         IDIM2LOOP = 0
C
    1    CONTINUE
C
         IDIMLOOP = IDIMLOOP + 1
         IF (IDIMLOOP .GE. NINSCA .AND. NBLOPT.NE.0) IFIN=3
C
         NCON = 0
C
         CALL EQDIM(IDIM)
         IF (IDIM .EQ. 2) IDIM2LOOP = IDIM2LOOP + 1
C
         IF (NSURF .EQ. 1) IDIM = 2
C
         CALL MESH(4)
         CALL INITIA
         CALL OUTPUT(4)
         CALL MESH(1)
         CALL MESH(2)
C
         IF (NSURF  .NE. 1) CALL GUESS(IGUESS)
         IF (NBLOPT .NE. 0 .OR. NBSOPT .NE. 0) CALL MESH(6)
C
         IGUESS = 2
C
C        CALL OUTPUT(7)
         CALL MATRIX
         CALL OUTPUT(10)
         CALL ITIPR
         CALL OLDNEW
         CALL RUNTIM
C
         IF (IDIM .EQ. 0) THEN
C
            IF ((NBSOPT.EQ.0 .AND. NBLOPT.EQ.0) .OR. NOPT.EQ.1) THEN 
C
               IDIM = 2
C
            ELSE
C
               IDIM = 1
C
            ENDIF
C
            R0  = RMAG
            RZ0 = RZMAG
            BPS( 1) = RMAG
            BPS(12) = RZMAG
C
            WRITE(*,*) 'BNDSPL ...'
            IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
            IF (NSURF .EQ. 6) CALL BNDSPL
C
            GOTO 1
C
         ELSE IF (IDIM .EQ. 1) THEN
C
            CALL PREMAP(2)
C
            IF (NBLOPT .NE. 0) THEN
C
               IF (NBLOPT .EQ. 1) CALL SCOPY(NPPR+1,CPPR,1,RPRM,1)
C
               CALL NOREPT(NPPR+1,1)
               CALL BLTEST
C
C IF IFIN = 3, BALLOONING STABILITY IS TESTED, BUT THE P-PRIME
C PROFILE IS NOT CHANGED.
C
               IF (IFIN .EQ. 3) THEN
C
                  CALL OUTPUT(17)
C
                  IDIM = 2
C
               ELSE
C
                  CALL PPRM(IFIN,IDIMLOOP)
                  CALL OUTPUT(17)
C
               ENDIF
C
               NBLOPT = 2
C
            ENDIF
C
            IF (NBSOPT .NE. 0) THEN
C
               CALL SCOPY(NPPR+1,CPPR,1,RPRM,1)
               CALL NOREPT(NPPR+1,1)
               CALL PPBSTR(IFIN)
C
               IF (IFIN .EQ. 3) THEN
C
                  CALL OUTPUT(17)
C
                  IDIM = 2
C
               ELSE
C
                  CALL OUTPUT(17)
C
               ENDIF
C
               NBSOPT = 2
C
            ENDIF
C
            CALL OLDEQ
C
            R0  = RMAG
            RZ0 = RZMAG
            BPS( 1) = RMAG
            BPS(12) = RZMAG
C
            IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
            IF (NSURF .EQ. 6) CALL BNDSPL
C
            GOTO 1
C
         ELSE IF (IDIM.EQ.2 .AND. (IDIMLOOP.LE.2  .OR.
     +        (NSTTP.EQ.3 .AND. IDIM2LOOP.LE.3)) ) THEN
C
C     THERE MIGHT BE A PROBLEM IF PSIISO(1) WILL BE .LT. CPSICL(1)
C     IN THIS CASE, PROFILE CALLS ISOFIND(K1>1) AND TETPSI(.,1), FOR E.G.,
C     IS NOT COMPUTED IN ISOFIND. THIS HAPPENS IF R0 IS TOO FAR FROM RMAG.
C     THEN SPSIM IS TOO MUCH LOWER THAN CPSICL(1) AND IT
C     LEAVES ROOM FOR PSIISO(1) TO BE IN BETWEEN SPSIM AND CPSICL(1).
C     PSIISO(1) IS TYPICALLY SPSIM*(1.-CSM(1)**2), WHICH FOR AN
C     EQUIDISTANT MESH IS: SPSIM*(1.-(1/2/NPSI)**2).
C
c%OS            ZMACON = SQRT((RMAG-R0)**2+(RZMAG-RZ0)**2)
c%OS            IF (ZMACON .GT. 0.3/ASPCT/DFLOAT(NPSI)) THEN
C
            ZPSITEST = SPSIM*(1. - 1./(2.*NPSI)**2)
            WRITE(6,'(/," IN BALLIT : ZPSITEST= ",1PE15.8,
     +        "  CPSICL(1)= ",E15.8,/)') ZPSITEST, CPSICL(1)
C
C     EXTRA ITERATION ON BIG MESH
            IF (ZPSITEST.LT.CPSICL(1) .OR.
     +        (NIDEAL.EQ.6 .AND. IDIM2LOOP.LE.1)) THEN
               R0  = RMAG
               RZ0 = RZMAG
               BPS( 1) = RMAG
               BPS(12) = RZMAG
C
               IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
               IF (NSURF .EQ. 6) CALL BNDSPL
C
               GOTO 1
            ENDIF
         ENDIF
C
         IF (IDIM.EQ.2 .AND. IDIMLOOP.GT.2) THEN
            ZPSITEST = SPSIM*(1. - 1./(2.*NPSI)**2)
            WRITE(6,'(/," IN BALLIT : ZPSITEST= ",1PE15.8,
     +        "  CPSICL(1)= ",E15.8,/)') ZPSITEST, CPSICL(1)
         ENDIF
C
         CALL CHECK
         CALL SCOPY(NISO,TMF,1,TMFO,1)
C
         RETURN
         END
C*DECK C2SP10
C*CALL PROCESS
         SUBROUTINE PPRM(KFIN,ILOOP)
C        ###########################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP10  MODIFY P' DURING BALLOONING OPTIMIZATION ACCORDING TO THE  *
*          ALGORITHM DESCRIBED IN SECTION 5.5 OF THE PUBLICATION      *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMOPT.inc'
C
         PARAMETER(ZEPS = 1.E-3,EPSDIF = 0.2, CSMIN = 0.05)
C
         DIMENSION ZDIFF(NPISO), ZWORK(NPISO),ZWORK1(NPISO)
C
         ZERR  = 0.
         ISTAB = 10
         CALL GLOQUA(PCSM,PCS,NPPR+1,2)
C
         WRITE(*,1)(PCSM(J),NP4(J),NP3(J),NP2(J),NP1(J),NP0(J),NCBAL(J),
     &     RJBSOS(J,2)/RJPAR(J),SMERCR(J),J=1,NPPR+1)
 1       FORMAT('        CS  NP4  NP3  NP2  NP1  NP0  NCBAL J_BS/J_PAR',
     &    '     -D_RES',/,(F8.3,5I5,I7,2F11.4))
         WRITE(6,100) 'UNSTABLES VALUES : '
C
      WRITE(*,*) 'pplimit'
      WRITE(*,'(/(5f12.4))') 
     &       (PPLIMIT(PCSM(I)**2),I=1,NPPR+1)
         DO 11 J11=1,NPPR+1
C
         XP3(J11) = XP2(J11)
         XP2(J11) = XP1(J11)
         XP1(J11) = XP0(J11)
         XP0(J11) = RPRM(J11)
C
         NP3(J11) = NP2(J11)
         NP2(J11) = NP1(J11)
         NP1(J11) = NP0(J11)
C
         ZBSFRAC = RJBSOS(J11,2)/RJPAR(J11)
         RADVOL = SQRT(VSURF(J11)/VSURF(NPPR+1))
         RPEDESTAL = AT(9)
         IF (RADVOL .GT. RPEDESTAL) GOTO 2
         IF (SMERCI(J11) .GE. 0. 
     &   .AND. ABS(XP0(J11)) .LT. 
     &      CFBAL*PPLIMIT(PCSM(J11)**2)
     &   .AND. NCBAL(J11) .EQ. 0) THEN 
C STABLE
            IF (NP0(J11) .NE. 0) NP0(J11) = 1
C
            IF (XP0(J11).NE.0.0) THEN
C
               IF (ISTAB.EQ.10) ISTAB=1
               IF (ISTAB.EQ.2)  ISTAB=3
C
            ENDIF
C
         ELSE IF (SMERCI(J11) .LT. 0. 
C    &    .OR. SMERCR(J11) .LT. 0.
     &    .OR. ABS(XP0(J11)) .GT. 
     &        CFBAL*PPLIMIT(PCSM(J11)**2)
     &    .OR. NCBAL(J11) .GE. 1) THEN 
C UNSTABLE
C           WRITE(6,200) J11,XP0(J11),ZBSFRAC
C
            IF (NP0(J11) .NE. 0) NP0(J11) = - 1
C
            XP4(J11) = XP0(J11)
            NP4(J11) = NP0(J11)
C
            IF (XP0(J11) .NE. 0.) THEN
C
               IF (ISTAB .EQ. 10) ISTAB=2
               IF (ISTAB .EQ. 1)  ISTAB=3
C
            ENDIF
         ENDIF
         GOTO 11
CJET
CJET New section for pedestal
CJET
 2       CONTINUE
          write(*,*) ' NEW CONDITION USED, J11 =',J11
     &              ,' NPPR=',NPPR
         write(*,*) ' MERCIER=',SMERCI(J11)
         write(*,*) ' ZBSFRAC=',ZBSFRAC
         write(*,*) ' NCBAL  =',NCBAL(J11)
         IF (SMERCI(J11) .GE. 0. 
     &     .AND. ZBSFRAC . LT. 1.
     &     .AND. NCBAL(J11) .EQ. 0) THEN 
C STABLE
            IF (NP0(J11) .NE. 0) NP0(J11) = 1
C
            IF (XP0(J11).NE.0.0) THEN
C
               IF (ISTAB.EQ.10) ISTAB=1
               IF (ISTAB.EQ.2)  ISTAB=3
C
            ENDIF
C
         ELSE IF (SMERCI(J11) .LT. 0. 
     &      .OR. ZBSFRAC . GT. 1.
     &      .OR. NCBAL(J11) .GE. 1) THEN 
C UNSTABLE
C           WRITE(6,200) J11,XP0(J11),ZBSFRAC
C
            IF (NP0(J11) .NE. 0) NP0(J11) = - 1
C
            XP4(J11) = XP0(J11)
            NP4(J11) = NP0(J11)
C
            IF (XP0(J11) .NE. 0.) THEN
C
               IF (ISTAB .EQ. 10) ISTAB=2
               IF (ISTAB .EQ. 1)  ISTAB=3
C
            ENDIF
         ENDIF
C
  11     CONTINUE
C
         WRITE(95,400) (PCSM(J),XP0(J),J=1,NPPR+1)
         WRITE(95,410) (PCSM(J),NP0(J),J=1,NPPR+1)
C
         WRITE(6,*)
C
************************************************************************
*                                                                      *
* ISTAB=1 IF THE PROFILE IS STABLE   EVERYWHERE                        * 
*      =2 IF THE PROFILE IS UNSTABLE EVERYWHERE                        * 
*      =3 OTHERWISE                                                    * 
*                                                                      *
************************************************************************
C
         IF (ISTAB .EQ. 1) THEN
C
            IF (NSRCH .EQ. 2) THEN
C
               DO 12 J12=1,NPPR+1
C
               XPPRMN(J12) = XP0(J12)
               XLAMB(J12)  = 0.
               XP0(J12)    = XP4(J12)
               NP0(J12)    = NP4(J12)
C
  12           CONTINUE
C
                WRITE(*,*) 'MINIMUM FOUND'
                WRITE(*,'(F15.8)') (XPPRMN(J),J=1,NPPR+1)
C
               NSRCH = 2
C
            ENDIF 
C
         ENDIF
C
         IF (ISTAB .EQ. 2) THEN
C
            IF (NSRCH .EQ. 1) THEN
C
               DO 13 J13=1,NPPR+1
C
               XPPRMN(J13) = 0.
               XPPRMX(J13) = XP0(J13)
               XPPRDF(J13) = XPPRMX(J13) - XPPRMN(J13)
               XLAMB(J13)  = 1.
               NP1(J13)    = 0
C
   13          CONTINUE
C
               WRITE(*,*)
               WRITE(*,*) 'MAXIMUM FOUND'
C
               NSRCH = 3
C
            ENDIF
C
         ENDIF
C
************************************************************************
*                                                                      *
* KFIN = -1 : FIRST ITERATION                                          *
*         0 : ERROR LARGER THAN ZEPS                                   *
*         1 : ERROR SMALLER THAN ZEPS BUT NOT STABLE EVERYWHERE        *
*         2 : ERROR SMALLER THAN ZEPS AND STABLE EVERYWHERE  OR        *
*             ERROR SMALLER THAN ZEPS / 10.0                           *
*         3 : CONVERGED                                                *
*                                                                      *
************************************************************************
C
         IF (KFIN .EQ. -1) THEN
C
            DO 14 J14=1,NPPR+1
C
            IF (PCSM(J14) .LT. .1) THEN
C
               XP0(J14) = 0.
               NP0(J14) = 0
C
            ENDIF
C
            XP4(J14) = XP0(J14)
            NP4(J14) = NP0(J14)
C
            XPPRMN(J14) = 0.
            XPPRMX(J14) = XP0(J14)
            XPPRDF(J14) = XP0(J14)
C
   14       CONTINUE
C
            KFIN = 0
C
         ENDIF
C
************************************************************************
*                                                                      *
* IF NSRCH=1, SEARCH A COMPLETELY UNSTABLE PROFILE                     *
*          2,                     STABLE   PROFILE                     *
*          3, SEARCH OPTIMIZED PROFILE (1st time)                      *
*          4,        OPTIMIZED PROFILE                                 *
*                                                                      *
************************************************************************
C
         IF (NSRCH .EQ. 1) THEN
C
            DO 15 J15=1,NPPR+1
C
            IF (NP0(J15) .EQ. 1 .AND. XP0(J15) .NE. 0.) THEN
C
               XLAMB(J15) = 1.1
C
               IF (XP0(J15) .NE. 0.) XP0(J15) = XP0(J15) - 0.02
C
            ELSE
C
               XLAMB(J15)=1.
C
            ENDIF
C
            RPRM(J15) = XLAMB(J15) * XP0(J15)
C
            IF (RPRM(J15) .GT. 0.) RPRM(J15) = 0.
C
            ZDIFF(J15) = RPRM(J15) - XP0(J15)
            AZDIFF = ABS(ZDIFF(J15))
C
            IF (AZDIFF .GT. ZERR) ZERR = AZDIFF
C
   15       CONTINUE
C
            GOTO 9000
C      
         ELSE IF (NSRCH .EQ. 2) THEN
C
            DO 16 J16=1,NPPR+1
C
            IF (NP0(J16) .EQ. -1 .AND. XP0(J16) .NE. 0.) THEN
C
               XLAMB(J16) = .5
C
            ELSE
C 
               XLAMB(J16) = 1.
C
            ENDIF
C
            RPRM(J16) = XLAMB(J16)*XP0(J16) 
C
            IF (RPRM(J16) .GT. 0.) RPRM(J16) = 0.
C
            ZDIFF(J16) = RPRM(J16) - XP0(J16)
            AZDIFF = ABS(ZDIFF(J16))
C
            IF (AZDIFF .GT. ZERR) ZERR = AZDIFF
C
   16       CONTINUE
C      
            GOTO 9000
C
         ELSE IF (NSRCH .EQ. 3) THEN
C
            DO 17 J17=1,NPPR+1
            ZDIFF(J17) = 0.
C
            X2SRCH(J17) = 0.1
            XLAMB(J17)  = 1.
C
            IF (NP0(J17) .EQ. -1) THEN
C
               IF (XP0(J17) .GT. -.005 .AND. XP0(J17-1) .EQ. 0.) THEN
C
                  WRITE(*,*) J17,XP0(J17),' MIS A ZERO'
C
                  XP0(J17)    = 0.
                  XPPRMN(J17) = 0.
                  XPPRDF(J17) = 0.
C
               ENDIF
C
               XLAMB(J17) = XLAMB(J17)  - X2SRCH(J17)
C
            ELSE 
C
               XLAMB(J17) = XLAMB(J17)  + X2SRCH(J17)
C
            ENDIF

C
            RPRM(J17)  = XPPRMN(J17) + XLAMB(J17) * XPPRDF(J17)
C
            IF (RPRM(J17) .GT. 0.) RPRM(J17) = 0.
C
            ZDIFF(J17) = RPRM(J17) - XP0(J17)
            AZDIFF = ABS(ZDIFF(J17))
C
            IF (AZDIFF .GT. ZERR) ZERR = AZDIFF
C
   17       CONTINUE
C
            NSRCH = 4
C
            GOTO 9500
C  
         ELSE IF (NSRCH .EQ. 4) THEN
C
            DO 18 J18=1,NPPR+1
            ZDIFF(J18) = 0.
C
            IF (XPPRDF(J18) .EQ. 0.) GOTO 18
C
            IP01 = NP0(J18) * NP1(J18)
            ISUM = ABS(NP1(J18) + NP2(J18) + NP3(J18))
C
            IF (IP01 .LT. 0) THEN
C 
               X2SRCH(J18) = .5 * ABS(X2SRCH(J18))
C
            ELSE IF (ISUM .EQ. 3 .AND. KFIN .EQ. 0) THEN
C
             X2SRCH(J18) = (1.2 + 3./DFLOAT(10+ILOOP))*ABS(X2SRCH(J18))
C
            ENDIF
            IF (PCSM(J18).LT. CSMIN) THEN
               XP0(J18)    = 0.
               XPPRMN(J18) = 0.
               XPPRDF(J18) = 0.
               RPRM(J18)   = 0.
            ENDIF
C
            IF (NP0(J18) .EQ. -1) THEN
C
               IF (XP0(J18) .GT. -.005 .AND. XP0(J18-1) .EQ. 0.) THEN
C
                  WRITE(*,*) J18,XP0(J18),' MIS A ZERO'
C
                  XP0(J18)    = 0.
                  XPPRMN(J18) = 0.
                  XPPRDF(J18) = 0.
C
               ENDIF
C
               X2SRCH(J18) = - X2SRCH(J18)
C
            ENDIF
C
C
 18    CONTINUE
C
cab
cab   Section added to make the profile prefer to relax to smooth function
cab
      DO 19 J19 = 2,NPPR
      ZDPROF = XP0(J19) - 0.5*(XP0(J19+1)+XP0(J19-1))
cab NOTE THAT P-PRIME IS NEGATIVE
      IF (ZDPROF*X2SRCH(J19).LT.0.) THEN
          IF (X2SRCH(J19-1)*X2SRCH(J19).LT.0..OR. 
     &        X2SRCH(J19)*X2SRCH(J19+1).LT.0.)
     &          X2SRCH(J19) = X2SRCH(J19) * 0.8
      ENDIF
      IF (ZDPROF*X2SRCH(J19).GT.0.) X2SRCH(J19) = 1.1*X2SRCH(J19)
cab
cab   Smoothing change
cab
      IF (X2SRCH(J19-1)*X2SRCH(J19).GT.0..AND.
     &    X2SRCH(J19)*X2SRCH(J19+1).GT.0.)
     &          X2SRCH(J19) = X2SRCH(J19) * 1.1
 19   CONTINUE
      DO 29 J29 = 1,NPPR+1,NPPR
      IF (J29.EQ.1) JJ = 2
      IF (J29.EQ.NPPR+1) JJ = NPPR
cab
cab   Smoothing change
cab
      IF (X2SRCH(JJ)*X2SRCH(J29).LT.0.) 
     &            X2SRCH(J29) = X2SRCH(J29) * 0.9
      IF (X2SRCH(JJ)*X2SRCH(J29).GT.0.) 
     &            X2SRCH(J29) = X2SRCH(J29) * 1.1
C
 29   CONTINUE
      ZERR = 0.
      DO 20 J20 = 1,NPPR+1
      ZDIFF(J20) = X2SRCH(J20) * XPPRDF(J20)
cab
cab   Avoid RPRM < 0
cab
      ZDIFFM = 0.5*ABS(RPRM(J20))
      IF (ZDIFFM.GT.EPSDIF) ZDIFFM = EPSDIF
      ZRED = 1.
      IF (ABS(ZDIFF(J20)).GT.ZDIFFM) ZRED = ZDIFFM/ABS(ZDIFF(J20))
      IF (ABS(ZDIFF(J20)).LT.0.5*ZEPS.AND.ZDIFF(J20).NE.0.) 
     &                               ZRED = .5*ZEPS/ABS(ZDIFF(J20))
      ZDIFF(J20) = ZDIFF(J20)*ZRED
      X2SRCH(J20) = X2SRCH(J20)*ZRED
C
            RPRM(J20) = RPRM(J20) + ZDIFF(J20)
            XLAMB(J20) = XLAMB(J20) + X2SRCH(J20)
            AZDIFF = ABS(ZDIFF(J20))
            IF (AZDIFF .GT. ZERR) ZERR = AZDIFF
 20   CONTINUE
C
            GOTO 9500
C
         ENDIF
C      
************************************************************************
*                                                                      *
*                            OUTPUT                                    *
*                                                                      *
************************************************************************
C
 9000   CONTINUE
C
        WRITE(*,*) 
        WRITE(*,*) 'ISTAB : ',ISTAB,NSRCH
        WRITE(*,*) 
        WRITE(*,*) '   I       CSM       P1P  IP1',
     $             '       P0P  IP0',
     $             '      LAMBDA          PP'
        DO 9100 J910=1,NPPR+1
            WRITE(*,250) J910,PCSM(J910),XP1(J910),NP1(J910),
     $                   XP0(J910),NP0(J910),
     $                   XLAMB(J910),RPRM(J910)
 9100   CONTINUE
         WRITE(*,*) 
         WRITE(*,*) 'ERROR  : ',ZERR,ZEPS
C
         CALL SPLINE(PCSM,RPRM,NPPR+1,D2RPRM,ZWORK,ZWORK1)
C
         RETURN
C
 9500   CONTINUE
C
        WRITE(*,*) 
        WRITE(*,*) 'ISTAB : ',ISTAB,NSRCH-3
        WRITE(*,*) 
        WRITE(*,9550)
 9550   FORMAT(11X,'CSM',11X,'PP0',12X,'PP',' STAB',8X,'LAMBDA',
     &         9X,'PDIFF')
        DO 9600 J960=1,NPPR+1
            WRITE(*,260) PCSM(J960),-XP0(J960),-RPRM(J960),NP0(J960),
     $                   XLAMB(J960),-ZDIFF(J960)
 9600    CONTINUE
C
         WRITE(*,*) 
C
cab
cab    This section stops iteration on pressure profile
cab    if beta changes too little over NSAVEP iterations
cab
      DO 9601 J=NSAVEB,2,-1
 9601 SBETA(J) = SBETA(J-1)
      SBETA(1) = BETAX
      ZBMAX = SBETA(1)
      ZBMIN = SBETA(1)
      DO 9602 J = 2,NSAVEB
      IF (SBETA(J).GT.ZBMAX) ZBMAX = SBETA(J)
      IF (SBETA(J).LT.ZBMIN) ZBMIN = SBETA(J)
 9602 CONTINUE
      WRITE(*,*) 'ITERATION ',ILOOP,'.  TOROIDAL BETA (EXP) OVER ',
     &           NSAVEB,' LAST ITERATIONS:'
      WRITE(*,9650) ZBMAX,ZBMIN,ZBMAX-ZBMIN
 9650 FORMAT(' MAX = ',E14.5,',   MIN = ',E14.5,'  VARIATION =',E14.5)
      IF (ZBMAX-ZBMIN.LT.EPSBET.AND.ZBMIN.GT.0.) THEN
         WRITE(*,9651) EPSBET,NSAVEB
         KFIN = 3
 9651 FORMAT('BETA-VARIATION LESS THAN',E14.5,' OVER',I4,' ITERATIONS')
      ENDIF
C
C
         IF (ZERR .LE. ZEPS) THEN
C
            IF (ZERR .EQ. 0. .AND. ISTAB .NE. 1) STOP 'ZERR=0&ISTAB?1'
C
            WRITE(6,110) 'EPSILON > ERROR > EPSILON / 5. '
            WRITE(6,*) ZEPS,ZERR,ZEPS / 5.
C
            KFIN = 3
C
         ELSE
C
            WRITE(*,*) 'ERROR > EPSILON ',ZERR,ZEPS
C
         ENDIF
C
         CALL SPLINE(PCSM,RPRM,NPPR+1,D2RPRM,ZWORK,ZWORK1)
C
         RETURN
C
  100    FORMAT (A)
  110    FORMAT (A,$)
  120    FORMAT (A,F10.3,A)
C
  200    FORMAT (I8,2X,2F15.8,$)
  240    FORMAT (I5,3(2X,F15.8),I5)
  250    FORMAT (I5,2X,F15.8,2(2X,F15.8,I5),2(2X,F15.8))
  260    FORMAT (3(F14.5),I5,2(F14.5))
C
  300    FORMAT (F6.2) 
  400    FORMAT (E15.8,2X,E15.8)
  410    FORMAT (E15.8,2X,I5)
C
  500    FORMAT ('THE ',F10.4,' % OPT PROFILE STABLE EVERYWHERE')
  510    FORMAT ('THE ',F10.4,' % OPT PROFILE NOT STABLE EVERYWHERE')
  520    FORMAT ('THE ',F10.4,' % OPT PROFILE NOT UNSTABLE EVERYWHERE')
  530    FORMAT ('THE ',F10.4,' % OPT PROFILE UNSTABLE EVERYWHERE')
C
  610    FORMAT (I4) 
  620    FORMAT (1PE15.8) 
C
         END
CJET
CJET  Function added to limit pprime in optimization
CJET  (generalization of CFBAL procedure)
CJET
      FUNCTION PPLIMIT(PSIN)
      INCLUDE 'COMPHY.inc'
      REAL PPLILIT,PSIN,X
      IF (PSIN .GT. 1. .OR. PSIN .LT. 0.) THEN
        WRITE(*,*) 'WARNING IN PPLIMIT: PSIN =',PSIN
      END IF
      IF (PSIN .LT. 0.) PSIN = 0.
      IF (PSIN .GT. 1.) PSIN = 1.
C     X = (1 - PSIN)**2
C     PPLIMIT = X*(0.25 + X*X)
      CALL PPRIM0(1, SPSIM*(1. - PSIN), PPLIMIT)
      PPLIMIT = ABS(PPLIMIT)
      RETURN
      END
C*DECK C2S03
C*CALL PROCESS
         SUBROUTINE INITIA
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2S03 INITIALIZE VERTICAL MATRIX NUMBERING AND GAUSS QUADRATURE     *
*       POINTS FOR THE INTEGRATION OF THE VARIATIONAL FORM OF THE     *
*       EQUILIBRIUM (SEE EQ. (27) IN THE PUBLICATION)                 *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     I   MPLA1(4,4),   MPLA2(4,4),   MPLA3(4,4),   MPLA4(4,4),
     R   ZRACS(NPSGS), ZRACT(NPTGS), ZWGTS(NPSGS), ZWGTT(NPTGS)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         CALL RESETI(MPLA1,16,0)
         CALL RESETI(MPLA2,16,0)
         CALL RESETI(MPLA3,16,0)
         CALL RESETI(MPLA4,16,0)
C
***********************************************************************
*                                                                     *
* 1. DEFINE HORIZONTAL BAND MATRIX ADRESSES                           *
*                                                                     *
***********************************************************************
C
         MPLA1(1,1) = 1
         MPLA1(2,1) = NT1
         MPLA1(3,1) = 2
         MPLA1(4,1) = NT + 2
C
         MPLA1(2,2) = 1
         MPLA1(4,2) = 2
C
         MPLA1(2,3) = NT
         MPLA1(3,3) = 1
         MPLA1(4,3) = NT1
C
         MPLA1(4,4) = 1
C
         MPLA2(1,1) = 1
         MPLA2(2,1) = NT1
         MPLA2(3,1) = 3
         MPLA2(4,1) = NT + 3
C
         MPLA2(2,2) = 1
         MPLA2(4,2) = 3
C
         MPLA2(2,3) = NT - 1
         MPLA2(3,3) = 1
         MPLA2(4,3) = NT1
C
         MPLA2(4,4) = 1
C
         MPLA3(1,1) = 1
         MPLA3(2,1) = NT1
         MPLA3(4,1) = NT - 1
C
         MPLA3(1,3) = 3
         MPLA3(2,3) = NT + 3
         MPLA3(3,3) = 1
         MPLA3(4,3) = NT1
C
         MPLA3(2,2) = 1
C
         MPLA3(2,4) = 3
         MPLA3(4,4) = 1
C
         MPLA4(1,1) = 1
         MPLA4(2,1) = NT1
         MPLA4(4,1) = NT
C
         MPLA4(1,3) = 2
         MPLA4(2,3) = NT + 2
         MPLA4(3,3) = 1
         MPLA4(4,3) = NT1
C
         MPLA4(2,2) = 1
C
         MPLA4(2,4) = 2
         MPLA4(4,4) = 1
C
         DO 1 I=1,4
         DO 1 JT=1,NT
C
         IF (JT .EQ. 1) THEN
C
         MPLA(JT,I,1) = MPLA1(I,1)
         MPLA(JT,I,2) = MPLA1(I,2)
         MPLA(JT,I,3) = MPLA1(I,3)
         MPLA(JT,I,4) = MPLA1(I,4)
C
         ELSE IF (JT .GE. 2 .AND. JT .LE. NT/2) THEN
C
         MPLA(JT,I,1) = MPLA2(I,1)
         MPLA(JT,I,2) = MPLA2(I,2)
         MPLA(JT,I,3) = MPLA2(I,3)
         MPLA(JT,I,4) = MPLA2(I,4)
C
         ELSE IF (JT .GE. NT/2+2 .AND. JT .LE. NT) THEN
C
         MPLA(JT,I,1) = MPLA3(I,1)
         MPLA(JT,I,2) = MPLA3(I,2)
         MPLA(JT,I,3) = MPLA3(I,3)
         MPLA(JT,I,4) = MPLA3(I,4)
C
         ELSE IF (JT .EQ. NT/2+1) THEN
C
         MPLA(JT,I,1) = MPLA4(I,1)
         MPLA(JT,I,2) = MPLA4(I,2)
         MPLA(JT,I,3) = MPLA4(I,3)
         MPLA(JT,I,4) = MPLA4(I,4)
C
         ENDIF
C           
    1    CONTINUE
C
***********************************************************************
*                                                                     *
* 2. DEFINE ISOPARAMETRIC INTEGRATION POINTS FOR GAUSSIAN QUADRATURE  *
*    (16 POINTS FORMULA)                                              *
*                                                                     *
***********************************************************************
C
         CALL GAUSS(NSGAUS,ZRACS,ZWGTS)
         CALL GAUSS(NTGAUS,ZRACT,ZWGTT)
C
         DO 6 J6=1,NSGAUS
C
         ZRACS(J6) = .5 * (1. + ZRACS(J6))
C
    6    CONTINUE
C
         DO 7 J7=1,NTGAUS
C
         ZRACT(J7) = .5 * (1. + ZRACT(J7))
C
    7    CONTINUE
C
***********************************************************************
*                                                                     *
* 3.DZETA CONTAINS VALUES (X ; Y) ON THE 16 GAUSSIAN QUADRATURE       *
*   POINTS. CW CONTAINS INTEGRATION WEIGHTS                           *
*                                                                     *
***********************************************************************
C
         DO 9 J9=1,NTGAUS
C
         DO 8 J8=1,NSGAUS
C
         KPOINT          = (J9 - 1) * NSGAUS + J8
         DZETA(KPOINT,1) = ZRACS(J8)
         DZETA(KPOINT,2) = ZRACT(J9)
         CW(KPOINT)      = ZWGTS(J8) * ZWGTT(J9)
C
    8    CONTINUE
    9    CONTINUE
C
         RETURN
         END
C*DECK C2S04
C*CALL PROCESS
         SUBROUTINE EQDIM(K)
C        ###################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2S04 CONTROL SIZE OF SMALL AND FINAL EQUILIBRIUM                   *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
C
         IF (K .EQ. 0) THEN
C
            IF (NSURF .NE. 1) THEN
C
               NDIM(1) = NS
               NDIM(2) = NT
               NDIM(3) = NISO
C
               NS      = MIN(24,NDIM(1))
               NT      = MIN(24,NDIM(2))
               NISO    = MIN(100,NDIM(3))
C
               CEPS    = MAX(EPNON0,EPSLON)
C
            ENDIF
C
            NS1    = NS + 1
            NT1    = NT + 1
            NT2    = NT + 2
            NPSI1  = NPSI + 1
            NCHI1  = NCHI + 1
            NV1    = NV + 1
            NBAND  = 4 * NT + 12
            NSTMAX = NS1 * NT
            N4NSNT = 4 * NSTMAX
            NWGAUS = NSGAUS * NTGAUS
C
         ELSE IF (K .EQ. 2) THEN
C
            NS     = NDIM(1)
            NT     = NDIM(2)
            NISO   = NDIM(3)
            CEPS   = EPSLON
            NS1    = NS + 1
            NT1    = NT + 1
            NT2    = NT + 2
            NBAND  = 4 * NT + 12
            NSTMAX = NS1 * NT
            N4NSNT = 4 * NSTMAX
            NWGAUS = NSGAUS * NTGAUS
C
         ENDIF
C
         RETURN
         END
C*DECK C2S05
C*CALL PROCESS
         SUBROUTINE MATRIX
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2S05 SET UP MATRIX A AND PERFORMS L * D * LT DECOMPOSITION OF A    *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMNUM.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C
         CALL SETUPA
         CALL ALDLT(A,RC1M14,N4NSNT,NBAND,NPBAND,ISGN)
C
         IF (ISGN .EQ. -1) STOP 'ISGN=-1'
C
         RETURN
         END
C*DECK C2S06
C*CALL PROCESS
         SUBROUTINE ITIPR
C        ################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2S06 LEAD ITERATION OVER CURRENT PROFILE                           *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   ZTTP(2*NPISO), ZCPPR(2*NPISO)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         DO 3 J3=1,NINMAP
C
         CALL NONLIN
C
         IF (NSTTP .EQ. 1) GOTO 4
C
         CALL SCOPY(NISO,TTP,1,ZTTP,1)
         CALL SCOPY(NISO,CPPR,1,ZCPPR,1)
C
         CALL PREMAP(1)
C
         ZNM  = 0.
         ZNDM = 0.
C
         DO 2 J2=1,NISO
C
         ZNM  = ZNM + (CPPR(J2) - ZCPPR(J2))**2 +
     +                (TTP(J2) - ZTTP(J2))**2
         ZNDM = ZNDM + TTP(J2)**2 + CPPR(J2)**2
C
    2    CONTINUE
C
         RESMAP = SQRT(ZNM/ZNDM)
C
C IF THE NSTTP = 3 : RESCALE AT'S AT EVERY ITERATION OVER CURRENT 
C PROFILE IF NO BALLOONING OPTIMIZATION IS USED OR THE OPTIMIZED 
C BALLOONING PRESSURE PROFILE IS RESCALED BY CPRESS. IF BALLOONING
C OPTIMIZATION IS USED, THIS OPERATION IS DONE ONLY AFTER EVERY
C STEP OVER THE PRESSURE PROFILE, SINCE IT HAS BEEN OBSERVED THAT 
C THE CONVERGENCE IS FASTER IN THAT WAY.
C
         IF (NSTTP .EQ. 3 .AND. 
     &       (NBLOPT .EQ. 0 .OR. (NOPT.EQ.1 .AND. 
     &       (NBLOPT.NE.0 .AND. CPRESS.NE.1.)))) THEN
C
            CALL NOREPT(NISO,0)
C
         ENDIF
C
         CALL OUTPUT(8)
         CALL CONVER(2,NCON)
C
         IF (NCON .EQ. 2) GOTO 4
C
    3    CONTINUE
C
C  ITERATION OVER CURRENT PROFILE HAS NOT CONVERGED
C  (WILL STOP AFTER HAVING SAVED THE EQUIL. IN STEPON)
C
         NCON = -2
C
         CALL OUTPUT(14)
C
    4    CONTINUE
C
         RETURN
         END
C*DECK C2S07
C*CALL PROCESS
         SUBROUTINE NONLIN
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2S07 LEAD ITERATION OVER THE NONLINEARITY                          *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         DO 1 J1=1,NINSCA
C
C SET UP RIGHT HAND SIDE B 
C
         CALL SETUPB
C
C SAVE OLD SOLUTION IN CPSIO;
C SOLVE SYSTEM L * D * LT * X = B
C
         CALL SCOPY(N4NSNT,CPSICL,1,CPSIO,1)
         CALL SOLVIT
         call scopyr(relax,n4nsnt,cpsio,1,cpsicl,1)
C
C FIND PSIMIN AND MAGNETIC AXIS
C
         CALL MAGAXE
C
C PRINT OUT :
C          - SPSIM, RMAG AND RZMAG FOR SOLOVEV CASE
C          - SPSIM, RMAG ,RZMAG, RESIDU AND EPSLON
C            FOR OTHER CASES
C
         CALL OUTPUT(5)
C
         IF (NSMOOTH .EQ. 1) CALL SMOOTH
C
         CALL ERROR1(CPSIO,CPSICL)
C
         IF (NSURF .EQ. 1) GOTO 2
C
         CALL OUTPUT(6)
C
C CONVERGENCE TEST
C
         CALL CONVER(1,NCON)
C
         IF (NCON .EQ. 1) GOTO 2
C
    1    CONTINUE
C
C ITERATION OVER NON-LINEARITY HAS NOT CONVERGED
C
         CALL OUTPUT(9)
C
    2    CONTINUE
C
         RETURN
         END
C*DECK C2S08
C*CALL PROCESS
         SUBROUTINE CHECK
C        ################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2S08 : CHECKS THE SOLUTION. IT COMPUTES :                          *
*                                                                     *
*         SQRT(SASUM((A(I,J)*CPSI(J) - B(I))**2) / N4NSNT, I=1,N4NSNT *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R   ZRES(NP4NST)
C
         INCLUDE 'BNDIND.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         CALL VZERO(ZRES,N4NSNT)
         CALL SETUPA
         CALL SETUPB
C
         DO 3 J3=1,N4NSNT
C
         I1 = MAX(1,J3-NBAND+1)
         I2 = MIN(N4NSNT,J3+NBAND-1)
C
         DO 2 J2=I1,I2
C
         ICOL = INDCOL(J2,J3)
         IROW = INDROW(J2,J3)
C
         ZRES(J3) = ZRES(J3) + A(ICOL,IROW) * CPSILI(J2)
C
    2    CONTINUE
C
         ZRES(J3) = (ZRES(J3) - B(J3)) * (ZRES(J3) - B(J3))
C
    3    CONTINUE
C
         SCHECK = SQRT(SSUM(N4NSNT,ZRES,1)) /DFLOAT(N4NSNT)
C
         RETURN
         END
C*DECK C2S09
C*CALL PROCESS
         SUBROUTINE OLDNEW
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2S09 : STORE CONVERGED SOLUTION OF GRAD-SHAFRANOV INTO AUXILLARY   *
*         ARRAYS                                                      *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         CALL SCOPY(NS1,CSIG,1,CSIGO,1)
         CALL SCOPY(NT1,CT,1,CTO,1)
         CALL SCOPY(NISO,CID0,1,CID0O,1)
         CALL SCOPY(NISO,CID2,1,CID2O,1)
         CALL SCOPY(NISO,CPPR,1,CPPRO,1)
         CALL SCOPY(NISO,D2CID0,1,D2CID0O,1)
         CALL SCOPY(NISO,D2CID2,1,D2CID2O,1)
         CALL SCOPY(NISO,D2CPPR,1,D2CPPRO,1)
         CALL SCOPY(NISO,TTP,1,TTPO,1)
         CALL SCOPY(NISO,TMF,1,TMFO,1)
         CALL SCOPY(NISO,D2TMF,1,D2TMFO,1)
         CALL SCOPY(NISO,CSIPR,1,CSIPRO,1)
         CALL SCOPY(N4NSNT,CPSICL,1,CPSIO,1)
         CALL SCOPY(14,BPS,1,BPSO,1)
C
         SPSIMO  = SPSIM
         R0O     = R0
         RZ0O    = RZ0
         R0WO    = R0W
         RZ0WO   = RZ0W
         RMAGO   = RMAG
         RZMAGO  = RZMAG
         NTO     = NT
         NSO     = NS
         NISOO   = NISO
C
         RETURN
         END
C*DECK C2S10
C*CALL PROCESS
         SUBROUTINE OLDEQ
C        ################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2S10 : READ CONVERGED SOLUTION OF GRAD-SHAFRANOV IN AUXILLARY      *
*         ARRAYS                                                      *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         CALL SCOPY(NS1,CSIGO,1,CSIG,1)
         CALL SCOPY(NT1,CTO,1,CT,1)
         CALL SCOPY(NISO,CID0O,1,CID0,1)
         CALL SCOPY(NISO,CID2O,1,CID2,1)
         CALL SCOPY(NISO,CPPRO,1,CPPR,1)
         CALL SCOPY(NISO,D2CID0O,1,D2CID0,1)
         CALL SCOPY(NISO,D2CID2O,1,D2CID2,1)
         CALL SCOPY(NISO,D2CPPRO,1,D2CPPR,1)
         CALL SCOPY(NISO,TTPO,1,TTP,1)
         CALL SCOPY(NISO,TMFO,1,TMF,1)
         CALL SCOPY(NISO,D2TMFO,1,D2TMF,1)
         CALL SCOPY(NISO,CSIPRO,1,CSIPR,1)
         CALL SCOPY(N4NSNT,CPSIO,1,CPSICL,1)
         CALL SCOPY(14,BPSO,1,BPS,1)
C
         SPSIM  = SPSIMO
         R0     = R0O
         RZ0    = RZ0O
         R0W    = R0WO
         RZ0W   = RZ0WO
         RMAG   = RMAGO
         RZMAG  = RZMAGO
         NT     = NTO
         NS     = NSO
         NISO   = NISOO
C
         RETURN
         END
C*DECK C2SA01
C*CALL PROCESS
         SUBROUTINE MESH(K)
C        ##################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SA01 : K = 1 COMPUTE S-MESH FOR I* OR J-PARALLEL.                 *
*                COMPUTE SIGMA & THETA MESHES FOR EQUILIBRIUM         *
*                DISCRETIZATION. COMPUTE PLASMA RADII RHOS(THETA)     *
*                AND D(RHOS(THETA)) / D(THETA).                       *
*                                                                     *
*          K = 2 COMPUTE ON EVERY GAUSSIAN INTEGRATION POINT FOR THE  *
*                EQUILIBRIUM SOLVER:                                  *
*                       - SIGMA & THETA                               *
*                       - RHOS(THETA) & D(RHOS(THETA) / D(THETA)      *
*                       - THE VALUES FB(I), I=1,...,16 OF THE BASIS   *
*                         FUNCTIONS                                   *
*                       - THE RELATION BETWEEN INVERSE CLOCKWISE AND  *
*                         UP-DOWN NUMEROTATION OF THE THETA MESH      *
*                                                                     *
*          K = 3 COMPUTE THE S-MESH FOR ERATO, LION, MARS AND         *
*                BALLOONING                                           *
*                                                                     *
*          K = 4 COMPUTE THE CHI-MESH FOR ERATO, LION, MARS AND       *
*                BALLOONING                                           *
*                                                                     *
*          K = 5 COMPUTE THE VACUUM S-MESH FOR MARS.                  *
*                                                                     *
*          K = 6 COMPUTE THE S-MESH FOR THE BALLOONING OPTIMIZATION   *
*                                                                     *
*          K = 7 COMPUTE THE S-MESH AND THE THETA-MESH FOR PENN       *
*                                                                     *
*          K = 8 THETA MESH FOR XTOR                                  *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMBND.inc'
C
         DIMENSION
     R   ZF(NPT,16),    ZRAC(NPMGS+1),ZRHOS(NTP1),    ZRHOS1(NTP1),  
     R   ZRHOS2(NTP1),  ZRHOS3(NTP1), ZRHOS4(NTP1),   ZS(NPT),       
     R   ZS1(NPT),      ZS2(NPT),     ZT(NTP1+NPCHI), ZT1(NPT),      
     R   ZT2(NPT),      ZT3(NPT),     ZT4(NPT),       ZWGT(NPMGS+1),
     R   CSOLD(NPISO),  ZQOLD(NPISO)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         GOTO (100,200,300,400,500,600,700,800) K
C
  100    CONTINUE
C
         ZDPIPR = 1. /DFLOAT(NISO - 1)
C
         DO 101 J101=1,NISO
 101        CSIPRI(J101) = (J101 - 1.) * ZDPIPR
C
***********************************************************************
*                                                                     *
*        NMESHB = 0 ===> EQUIDISTANT S-MESH                           *
*        NMESHB = 1 ===> WEIGHTED S-MESH                              *
*        NPOIDB = 0 ===> NO WEIGHTING IS POSSIBLE FOR S-MESH          *
*        SOLPDB = 1 ===> NO WEIGHTING IS POSSIBLE FOR S-MESH          *
*                                                                     *
***********************************************************************
C
         IF (NMESHB .NE. 0 .AND. NPOIDB .NE. 0 .AND. 
     +       SOLPDB .NE. 1.) THEN
            CALL PACKME(NISO,NPOIDB,CSIPRI,BPLACE,BWIDTH,SOLPDB)
         ENDIF
C
         DO 102 J102=1,NISO-1
 102        CSIPR(J102) = .5 * (CSIPRI(J102+1) + CSIPRI(J102))
C
         CSIPR(NISO) = 1.
C
***********************************************************************
*                                                                     *
*        FILL IN THETA VALUES : CONSTANT STEP = ZDT                   *
*                                                                     *
***********************************************************************
C
         ZDT = 2. * CPI /DFLOAT(NT)
C
         DO 103 J103=1,NT1
 103        CT(J103) = .5 * (2. * J103 - 3.) * ZDT
C
***********************************************************************
*                                                                     *
*        NMESHD = 0 ===> EQUIDISTANT THETA-MESH                       *
*        NMESHD = 1 ===> WEIGHTED THETA-MESH                          *
*        NPOIDD = 0 ===> NO WEIGHTING IS POSSIBLE FOR THETA-MESH      *
*        SOLPDD = 1 ===> NO WEIGHTING IS POSSIBLE FOR THETA-MESH      *
*     ASSUME DPLACE IN [-PI,+PI] TO EASE SYMMETRIC PACKING            *
*                                                                     *
***********************************************************************
C    
         ISHIFT = 0
C
         IF (NMESHD .NE. 0 .AND. NPOIDD .NE. 0 .AND. 
     +        SOLPDD .NE. 1.) THEN
C
            DO 104 J104=1,NPOIDD
               IF (DPLACE(J104) .GT. CPI) THEN
                 PRINT *,' DPLACE(',J104,')=',DPLACE(J104),
     +             ' SHOULD BE IN [-PI,PI]'
                 STOP
               ENDIF
               IF (DPLACE(J104) .LT. 0.)
     +                DPLACE(J104) = DPLACE(J104) + 2.*CPI
 104        CONTINUE
C
            CALL PACKMEP(NT1,NPOIDD,CT,DPLACE,DWIDTH,SOLPDD)
C
            DO 105 J105=1,NPOIDD
               IF (DPLACE(J105) .GT. CPI)
     +          DPLACE(J105) = DPLACE(J105) - 2.*CPI
 105        CONTINUE
C    
            ISHIFT = 1
         ENDIF
C
***********************************************************************
*                                                                     *
*    THETA MESH WITH CONSTANT AREA BETWEEN CONSTANT THETA SURFACES    *
*                                                                     *
***********************************************************************
C
         IF (NDIFT .NE. 0) THEN
            CALL TETARE(CT,NT)
            ISHIFT = 1
         ENDIF
C    
***********************************************************************
*     *
*     IF THETA MESH HAS BEEN PACKED OR NDIFT = 1 SHIFT THETA MESH FROM *
*     [0 ; 2*PI] TO  [(CT(NT)-2*PI)/2 ; 2*PI-(CT(NT)-2*PI)/2]          *
*     *
***********************************************************************
C     
         IF (ISHIFT .EQ. 1) THEN
            ZT(1) = .5 * (CT(1) + CT(NT) - 2.*CPI)
C     
            DO 106 J106=2,NT1
 106        ZT(J106) = .5 * (CT(J106) + CT(J106-1))
C     
            CALL SCOPY(NT1,ZT,1,CT,1)
         ENDIF
C
***********************************************************************
*                                                                     *
*        FILL IN SIGMA VALUES : CONSTANT STEP = ZDS                   *
*                                                                     *
***********************************************************************
C
         ZDS = 1. /DFLOAT(NS)
C
         DO 107 J107=1,NS1
 107        CSIG(J107) = (J107 - 1) * ZDS
C
***********************************************************************
*                                                                     *
*        NMESHC = 0 ===> EQUIDISTANT SIGMA-MESH                       *
*        NMESHC = 1 ===> WEIGHTED SIGMA-MESH                          *
*        NPOIDC = 0 ===> NO WEIGHTING IS POSSIBLE FOR SIGMA-MESH      *
*        SOLPDC = 1 ===> NO WEIGHTING IS POSSIBLE FOR SIGMA-MESH      *
*                                                                     *
***********************************************************************
C
         IF (NMESHC .NE. 0 .AND. NPOIDC .NE. 0 .AND. 
     +       SOLPDC .NE. 1.) THEN
            CALL PACKME(NS1,NPOIDC,CSIG,CPLACE,CWIDTH,SOLPDC)
         ENDIF
C
***********************************************************************
*                                                                     *
*        COMPUTE PLASMA SURFACE PARAMETERS:                           *
*                                                                     *
*        1) COMPUTE BOUNDARY VECTOR RADIUS : RHOS(CT(J))              *
*                                                                     *
***********************************************************************
C
         CALL BOUND(NT1,CT,RHOS)
C
         ZEPS = 1.E-3
C
         DO 108 J108=1,NT
            ZT1(J108) = CT(J108) - 2. * ZEPS
            ZT2(J108) = CT(J108) -      ZEPS
            ZT3(J108) = CT(J108) +      ZEPS
 108        ZT4(J108) = CT(J108) + 2. * ZEPS
C
         CALL BOUND(NT,ZT1,ZRHOS1)
         CALL BOUND(NT,ZT2,ZRHOS2)
         CALL BOUND(NT,ZT3,ZRHOS3)
         CALL BOUND(NT,ZT4,ZRHOS4)
C
         DO 109 J109=1,NT
 109        DRSDT(J109) = (ZRHOS1(J109)+8*(ZRHOS3(J109)-ZRHOS2(J109))-
     -                     ZRHOS4(J109))/(12.*ZEPS)
C
         DRSDT(NT1) = DRSDT(1)
C
         RETURN
C
 200     CONTINUE
C
***********************************************************************
*                                                                     *
* 4.5.2  COMPUTE RHO-SURFACE AND D(RHO-SURFACE)/DTHETA FOR EACH       *
*        THETA GIVEN BY THE POSITION OF THE 16 GAUSSIAN INTEGRATION   *
*        POINTS IN THE CELL.                                          *
*                                                                     *
***********************************************************************
C
         ZEPS = 1.E-3
C
         DO 204 J204=1,NWGAUS
            DO 201 J201=1,NT
               ZT(J201)  = CT(J201)+(CT(J201+1)-CT(J201))*DZETA(J204,2)
               ZT1(J201) = ZT(J201) - 2. * ZEPS
               ZT2(J201) = ZT(J201) -      ZEPS
               ZT3(J201) = ZT(J201) +      ZEPS
 201           ZT4(J201) = ZT(J201) + 2. * ZEPS
C
            CALL BOUND(NT,ZT,ZRHOS)
            CALL BOUND(NT,ZT1,ZRHOS1)
            CALL BOUND(NT,ZT2,ZRHOS2)
            CALL BOUND(NT,ZT3,ZRHOS3)
            CALL BOUND(NT,ZT4,ZRHOS4)
C
            DO 202 J202=1,NT
               YRST(J202,J204)   = ZRHOS(J202)
 202           YDRSDT(J202,J204) = (ZRHOS1(J202) + 8 * (ZRHOS3(J202) - 
     -                              ZRHOS2(J202)) - ZRHOS4(J202)) / 
     /                             (12. * ZEPS)
 204     CONTINUE
C
***********************************************************************
*                                                                     *
* 4.5.4. COMPUTE SIGMA AND THETA VALUES ON THE 16 INTEGRATION POINTS  *
*                                                                     *
***********************************************************************
C
         DO 207 J207=1,NWGAUS
            DO 205 J205=1,NS
 205           RSINT(J205,J207) = CSIG(J205)+(CSIG(J205+1)-CSIG(J205))*
     *                            DZETA(J207,1)
            DO 206 J206=1,NT
 206           RTINT(J206,J207) = CT(J206) + (CT(J206+1) - CT(J206)) *
     *                            DZETA(J207,2)
 207     CONTINUE
C
         DO 208 J208=1,NT
            ZT1(J208) = CT(J208)
 208        ZT2(J208) = CT(J208+1)
C
         DO 213 J213=1,NS
            CALL RESETR(ZS1,NT,CSIG(J213))
            CALL RESETR(ZS2,NT,CSIG(J213+1))
C
            DO 212 J212=1,NWGAUS
               DO 209 J209=1,NT
                  ZS(J209) = RSINT(J213,J212)
 209              ZT(J209) = RTINT(J209,J212)
C
            CALL BASIS1(NT,NPT,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZF)
C
               DO 211 J211=1,16
                  DO 210 J210=1,NT
 210                 FB(J210,J212,J211,J213) = ZF(J210,J211)
 211              CONTINUE
 212        CONTINUE
 213     CONTINUE
C
***********************************************************************
*                                                                     *
* 4.5.6. COMPUTE UP AND DOWN NUMEROTATION STARTING FROM INVERSE       *
*        CLOCKWISE NUMEROTATION                                       *
*                                                                     *
***********************************************************************
C
         DO 215 J215=1,NS1
            DO 214 J214=1,NT
               I = (J215 - 1) * NT + J214
               IF (J214 .GE. 2 .AND. J214 .LE. NT/2+1) THEN
                  NUPDWN(I) = I + J214 - 2
               ELSE IF (J214 .GE. NT/2+2 .AND. J214 .LE. NT) THEN
                  NUPDWN(I) = (J215 + 1) * NT - 2 * J214 + 3
               ELSE IF (J214 .EQ. 1) THEN
                  NUPDWN(I) = I
               ENDIF
 214        CONTINUE
 215     CONTINUE
C
         RETURN
C
 300     CONTINUE
C
***********************************************************************
*                                                                     *
* 1.     INITIALIZATION  AND  COMPUTING S-MESH                        *
*                                                                     *
* 1.1.   EQUIDISTANT S-MESH : ZDS = CONSTANT STEP IN S-DIRECTION      *
*                                                                     *
***********************************************************************
C
         ZDP = 1. /DFLOAT(NPSI)
C
         DO 301 J301=1,NPSI1
 301        CS(J301) = (J301 - 1) * ZDP
C
***********************************************************************
*                                                                     *
*        NMESHA = 0 ===> EQUIDISTANT S-MESH                           *
*        NMESHA = 1 ===> WEIGHTED S-MESH ON SPECIFIED S VALUES        *
*        NMESHA = 2 ===> WEIGHTED S-MESH ON SPECIFIED Q VALUES        *
*                                                                     *
*        NPOIDA = 0 AND NPOIDQ = 0                                    *
*                   ===> NO WEIGHTING IS POSSIBLE FOR S-MESH          *
*        SOLPDA = 1 ===> NO WEIGHTING IS POSSIBLE FOR S-MESH          *
*                                                                     *
***********************************************************************
C
         IF (NMESHA .NE. 0 .AND. (NPOIDA .NE. 0 .OR. NPOIDQ .NE. 0)
     +       .AND. SOLPDA .NE. 1.) THEN
C
            IF (NMESHA.GE.2) 
     +      CALL QPLACS(CSOLD,ZQOLD)
C
            IF (NPOIDA .NE. 0 .AND. SOLPDA .NE. 1.) THEN
               CALL PACKME(NPSI1,NPOIDA,CS,APLACE,AWIDTH,SOLPDA)
            ENDIF
         ENDIF
C
         DO 302 J302=1,NPSI
 302        CSM(J302) = .5 * (CS(J302+1) + CS(J302))
C
         CSM(NPSI1) = 1.
C
***********************************************************************
*                                                                     *
*    S-MESH WITH CONSTANT FLUX VOLUME BETWEEN SURFACES                *
*                                                                     *
***********************************************************************
C
         IF (NDIFPS .EQ. 1) CALL PSVOL

C YQL, 2005-03-29
C SHIFT S-MESH AROUND RESONANT Q-SURFACES
C        IF (NMESHA.GE.3) CALL QAUTOPACKB(CSOLD,ZQOLD) 
C
         RETURN
C
 400     CONTINUE
C
***********************************************************************
*                                                                     *
* 4. FILL IN EQUIDISTANT CHI-MESH                                     *
*                                                                     *
***********************************************************************
C
         ZDCHI = 2. * CPI /DFLOAT(NCHI)
C
         DO 401 J401=1,NCHI1
 401        CHI(J401) = (J401 - 1.5) * ZDCHI
C
***********************************************************************
*                                                                     *
*        NMESHE = 0 ===> EQUIDISTANT CHI-MESH                         *
*        NMESHE = 1 ===> WEIGHTED CHI-MESH                            *
*        NPOIDE = 0 ===> NO WEIGHTING IS POSSIBLE FOR CHI-MESH        *
*        SOLPDE = 1 ===> NO WEIGHTING IS POSSIBLE FOR CHI-MESH        *
*                                                                     *
***********************************************************************
C
         IF (NMESHE .NE. 0 .AND. NPOIDE .NE. 0 .AND. 
     +           SOLPDE .NE. 1.) THEN
C     
            DO 402 J402=1,NPOIDE
               IF (EPLACE(J402) .GT. CPI) THEN
                 PRINT *,' EPLACE(',J402,')=',EPLACE(J402),
     +             ' SHOULD BE IN [-PI,PI]'
                 STOP
               ENDIF
               IF (EPLACE(J402) .LT. 0.0)
     +                EPLACE(J402) = EPLACE(J402) + 2.*CPI
 402        CONTINUE

            CALL PACKMEP(NCHI1,NPOIDE,CHI,EPLACE,EWIDTH,SOLPDE)
C     
            DO 403 J403=1,NPOIDE
               IF (EPLACE(J403) .GT. CPI)
     +                EPLACE(J403) = EPLACE(J403) - 2.*CPI
 403        CONTINUE
C
            ZCSHFT = .5 * CHI(2)
C
            DO 404 J404=1,NCHI1
 404           CHI(J404) = CHI(J404) - ZCSHFT
         ENDIF
C
         DO 405 J405=1,NCHI
 405        CHIM(J405) = .5 * (CHI(J405) + CHI(J405+1))
C
         RETURN
C
 500     CONTINUE
C
 
         IF (NVEXP.EQ.0.AND.NV.GT.1) THEN
C
***********************************************************************
*                                                                     *
*    EQUIDISTANT VACUUM S-MESH FOR MARS                               *
*                                                                     *
***********************************************************************
C
            DO 501 J501 = 1, NV
 501           CSV(J501) = 1. + (REXT - 1.) *DFLOAT(J501-1) /DFLOAT(NV)
C
            CSV(NV+1) = REXT

            TMP_HV    = (REXT-1.0)/DFLOAT(NV)
            NW(1)     = 1
            NW(NWBPS) = NV + 1
            DO J=2,NWBPS-1
               NW(J) = 1 + INT(0.5+(RW(J)-1.0)/TMP_HV)
            ENDDO

            DO J=1,NWBPS
               WRITE(*,*) 'J=',J,' NW(J)=',NW(J)
            ENDDO

            DO J=1,5
               NCOIL = 1 + INT(0.5+(RCOIL(J)-1.0)/TMP_HV)
               WRITE(*,*) 'J=',J,' NCOIL=',NCOIL
            ENDDO

         ELSE IF ((NVEXP.EQ.1.OR.NVEXP.EQ.2.OR.NVEXP.EQ.5
     &            .OR.NVEXP.EQ.6).AND.NV.GT.1) THEN
C
CYQLIU, NOV 21,2006
CYQLIU, MAR 17,2007: added NVEXP=5,6,50,60
***********************************************************************
*                                                                     *
*    EQUIDISTANT VACUUM S-MESH FOR MARS, EXCEPT MESH POINTS CLOSE TO  *
*    WALLS ARE SHIFTED TO MATCH WALLS POSITIONS                       *
*                                                                     *
***********************************************************************
C
            IF ((NVEXP.EQ.5.OR.NVEXP.EQ.6).AND.NV.GT.6) THEN
              TMP_HV    = (REXT-1.0)/DFLOAT(NV-2)
              DO J501 = 1, 5
                 CSV(J501) = 1. + DFLOAT(J501-1)*TMP_HV*0.5 
              ENDDO
              DO J501 = 6, NV
                 CSV(J501) = CSV(5) + DFLOAT(J501-5)*TMP_HV 
              ENDDO
            ELSE
              DO J501 = 1, NV
                 CSV(J501) = 1.+(REXT-1.)*DFLOAT(J501-1)/DFLOAT(NV)
              ENDDO
            ENDIF
            CSV(NV+1) = REXT

            IF ((NVEXP.EQ.2.OR.NVEXP.EQ.6).AND.NV.GE.3.AND.
     &          (CSV(2)-CSV(1)).GT.(CS(NPSI1)-CS(NPSI))) THEN
               CSV(2) = CSV(1) + (CS(NPSI1)-CS(NPSI))
            ENDIF

            TMP_HV    = (REXT-1.0)/DFLOAT(NV)
            NW(1)     = 1
            NW(NWBPS) = NV + 1
            DO J=2,NWBPS-1
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
            ENDDO
            
            DO J=1,5
               IF (RCOIL(J).GT.1.0) THEN
               JRW = 1
               DO JV=1,NV
                 IF (CSV(JV).LE.RCOIL(J).AND.CSV(JV+1).GT.RCOIL(J)) THEN
                    JRW = JV
                 ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RCOIL(J)-CSV(JRW).LT.CSV(JRW+1)-RCOIL(J)) THEN
                  NCOIL = JRW
                  CSV(JRW) = RCOIL(J)
               ELSE 
                  NCOIL = JRW + 1
                  CSV(JRW+1) = RCOIL(J)
               ENDIF
               WRITE(*,5829) RCOIL(J),TMP1,TMP2,NCOIL
               ENDIF 
            ENDDO

         ELSE IF ((NVEXP.EQ.10.OR.NVEXP.EQ.20.OR.NVEXP.EQ.50
     &            .OR.NVEXP.EQ.60).AND.NV.GT.1) THEN
C
CYQLIU, JAN 03,2007
***********************************************************************
*                                                                     *
*    EQUIDISTANT VACUUM S-MESH FOR MARS, EXCEPT MESH POINTS CLOSE TO  *
*    WALLS ARE SHIFTED TO MATCH WALLS POSITIONS, PLUS SPECIFIC OPTION *
*    FOR ITER DOUBLE WALL WITH UNIFROM WALL THICKNESS                 *
***********************************************************************
C
            IF ((NVEXP.EQ.50.OR.NVEXP.EQ.60).AND.NV.GT.6) THEN
              TMP_HV    = (REXT-1.0)/DFLOAT(NV-2)
              DO J501 = 1, 5
                 CSV(J501) = 1. + DFLOAT(J501-1)*TMP_HV*0.5 
              ENDDO
              DO J501 = 6, NV
                 CSV(J501) = CSV(5) + DFLOAT(J501-5)*TMP_HV 
              ENDDO
            ELSE
              DO J501 = 1, NV
                 CSV(J501) = 1.+(REXT-1.)*DFLOAT(J501-1)/DFLOAT(NV)
              ENDDO
            ENDIF
            CSV(NV+1) = REXT

            IF ((NVEXP.EQ.20.OR.NVEXP.EQ.60).AND.NV.GE.3.AND.
     &          (CSV(2)-CSV(1)).GT.(CS(NPSI1)-CS(NPSI))) THEN
               CSV(2) = CSV(1) + (CS(NPSI1)-CS(NPSI))
            ENDIF

            TMP_HV    = (REXT-1.0)/DFLOAT(NV)
            NW(1)     = 1
            NW(NWBPS) = NV + 1
            DO J=3,3
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               NW(J-1) = NW(J) - 1
               NW(J+1) = NW(J) + 1
               CSV(NW(J-1)) = RW(J-1)
               CSV(NW(J+1)) = RW(J+1)
               WRITE(*,5829) RW(J-1),TMP1,TMP2,NW(J-1)
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
               WRITE(*,5829) RW(J+1),TMP1,TMP2,NW(J+1)
            ENDDO
            DO J=6,6
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               NW(J-1) = NW(J) - 1
               NW(J+1) = NW(J) + 1
               CSV(NW(J-1)) = RW(J-1)
               CSV(NW(J+1)) = RW(J+1)
               WRITE(*,5829) RW(J-1),TMP1,TMP2,NW(J-1)
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
               WRITE(*,5829) RW(J+1),TMP1,TMP2,NW(J+1)
            ENDDO
            DO J=8,NWBPS-1
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
            ENDDO
            
            DO J=1,5
               IF (RCOIL(J).GT.1.0) THEN
               JRW = 1
               DO JV=1,NV
                 IF (CSV(JV).LE.RCOIL(J).AND.CSV(JV+1).GT.RCOIL(J)) THEN
                    JRW = JV
                 ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RCOIL(J)-CSV(JRW).LT.CSV(JRW+1)-RCOIL(J)) THEN
                  NCOIL = JRW
                  CSV(JRW) = RCOIL(J)
               ELSE 
                  NCOIL = JRW + 1
                  CSV(JRW+1) = RCOIL(J)
               ENDIF
               WRITE(*,5829) RCOIL(J),TMP1,TMP2,NCOIL
               ENDIF 
            ENDDO

         ELSE IF (NVEXP.EQ.3.AND.NV.GT.1) THEN
C
***********************************************************************
*                                                                     *
*    EXPONENTIAL VACUUM S-MESH FOR MARS                               *
*                                                                     *
***********************************************************************
C
C.. FIRST DETERMINE ZALPHA TO SATISFY ABOVE RELATIONS
C
            ZQ     = (REXT - 1.) / (CS(NPSI1) - CS(NPSI))
            ZALPHA = 2.
C
            IF (ZQ .LT. NV) ZALPHA = 0.5
C
C...START VALUE FOR ZALPHA
C
            JITER = 0
C
  502       CONTINUE
C
            ZFA    = (ZALPHA**NV - 1.) / (ZALPHA - 1.) - ZQ
            JITER  = JITER + 1
            ZFPRIM = (((NV - 1.) * ZALPHA - NV) * ZALPHA**(NV-1) + 1.) /
     /               (ZALPHA - 1.)**2
            ZFA    = ZFA / ZFPRIM
            ZALPHA = ZALPHA - ZFA
C      
            IF (ABS(ZFA) .LT. RC1M12) GOTO 503
C
            IF (JITER .LT. 500) GOTO 502
C
            WRITE(*,*) ' ITERATION FOR ALPHA DOES NOT CONVERGE'
            STOP 'NO_CONV'
C
  503       CONTINUE
C
            ZD1    = CS(NPSI1) - CS(NPSI)
            CSV(1) = 1.
            ZEXP   = 1.
C      
            DO 504 J504=1,NV
               ZEXP        = ZEXP * ZALPHA
 504           CSV(J504+1) = 1. + ZD1 * (ZEXP - 1.) / (ZALPHA - 1.)

            NW(1) = 1
            NW(NWBPS) = NV + 1
         ELSE IF (NVEXP.EQ.4.AND.NV.GT.1) THEN

***********************************************************************
*                                                                     *
*    EXPONENTIAL VACUUM S-MESH FOR MARS-F, ALSO WALL POSITIONS AND    *
*    THE FIRST VACUUM GRID IS MOVED.                                  *
*                                                                     *
***********************************************************************
            ZQ     = (REXT - 1.) / (CS(NPSI1) - CS(NPSI))
            ZALPHA = 1.

            JITER = 0
 1502       CONTINUE

            ZLOGH  = LOG(HEND/(CS(NPSI1)-CS(NPSI)))
            ZFA    = 1.0 + HEND/(CS(NPSI1)-CS(NPSI)) - ZQ
            ZFPRIM = 0.0
            DO J=2,NV-1
               ZTMP1 = LOG(DFLOAT(J-1)/DFLOAT(NV-1))
               ZTMP2 = EXP(ZALPHA*ZTMP1)
               ZTMP3 = EXP(ZLOGH*ZTMP2)
               ZFA   = ZFA + ZTMP3
               ZFPRIM = ZFPRIM + ZTMP3*ZLOGH*ZTMP2*ZTMP1
            ENDDO
            JITER  = JITER + 1
            ZFA    = ZFA / ZFPRIM
            ZALPHA = ZALPHA - ZFA
C      
            IF (ABS(ZFA) .LT. RC1M12) GOTO 1503
C
            IF (JITER .LT. 500) GOTO 1502
C
            WRITE(*,*) ' ITERATION FOR ALPHA DOES NOT CONVERGE'
            STOP 'NO_CONV'
C
 1503       CONTINUE
C
            ZD1    = CS(NPSI1) - CS(NPSI)
            CSV(1) = 1.
            CSV(2) = CSV(1) + ZD1

            DO J=2,NV
               ZTMP1 = LOG(DFLOAT(J-1)/DFLOAT(NV-1))
               ZTMP2 = EXP(ZALPHA*ZTMP1)
               ZTMP3 = EXP(ZLOGH*ZTMP2)*ZD1
               CSV(J+1) = CSV(J) + ZTMP3
            ENDDO

            IF (ABS(CSV(NV+1)-REXT).GT.RC1M12) THEN
               WRITE(*,*) 'CSV(NV+1)=',CSV(NV+1),' REXT=',REXT
               STOP 'NVEXP=4 REXT'
            ENDIF

            WRITE(*,*) 'NVEXP = ',NVEXP,'   ZALPHA = ',ZALPHA
            WRITE(*,*) 'LAST CELL SIZE = ',CSV(NV+1)-CSV(NV)
            NW(1) = 1
            NW(NWBPS) = NV + 1
            DO J=2,NWBPS-1
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
            ENDDO
            
            DO J=1,5
               IF (RCOIL(J).GT.1.0) THEN
               JRW = 1
               DO JV=1,NV
                 IF (CSV(JV).LE.RCOIL(J).AND.CSV(JV+1).GT.RCOIL(J)) THEN
                    JRW = JV
                 ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RCOIL(J)-CSV(JRW).LT.CSV(JRW+1)-RCOIL(J)) THEN
                  NCOIL = JRW
                  CSV(JRW) = RCOIL(J)
               ELSE 
                  NCOIL = JRW + 1
                  CSV(JRW+1) = RCOIL(J)
               ENDIF
               WRITE(*,5829) RCOIL(J),TMP1,TMP2,NCOIL
               ENDIF 
            ENDDO
            
         ELSE IF ((NVEXP.EQ.8.OR.NVEXP.EQ.80.OR.NVEXP.EQ.800).AND.
     &            NV.GT.4.AND.NWBPS.GT.2) THEN
CYQLIU, 2007-06-07
***********************************************************************
*                                                                     *
*    NEW ROBUST EXPONENTIAL VACUUM S-MESH FOR MARS-F                  *
*    EXPONENTIAL MESH BETWEEN PLASMA SURFACE AND FIRST WALL           *
*    EQUIDISTANT MESH BETWEEN FIRST WALL AND OUTER BOUNDARY           *
*    MESH POSITIONS ARE SHIFTED TO MATCH WALL POSITIONS               *
*    THE FIRST VACUUM GRID SIZE MATCHES THE LAST PLASMA MESH SIZE     *
*    NVEXP = 80:  specific for ITER double wall + EFCC
*    NVEXP = 800: specific for ITER double wall + ELM  
*                                                                     *
***********************************************************************
C           COMPUTE NW FOR FIRST WALL = JITER
            ZTEMP1 = (RW(2) - 1.0)/(CS(NPSI1) - CS(NPSI))
            ZTEMP2 = (REXT - RW(2))/(CS(NPSI1) - CS(NPSI))
            JITER = 3
            ZTEMP3 = (1-ZTEMP1)/(ZTEMP2/(NV+1-JITER)-ZTEMP1)
            ZTEMP3 = ZTEMP3**(JITER-2)
            ZTEMP3 = ABS(ZTEMP3 - ZTEMP2/(NV+1-JITER))
            ZFA = ZTEMP3
            DO J=4,NV
              ZTEMP3 = ((1-ZTEMP1)/(ZTEMP2/(NV+1-J)-ZTEMP1))**(J-2)
              ZTEMP3 = ABS(ZTEMP3 - ZTEMP2/(NV+1-J))
              IF (ZTEMP3.LT.ZFA) THEN
                JITER = J
                ZFA   = ZTEMP3
              ENDIF 
            ENDDO

C           COMPUTE ZALPHA
            ZTEMP3 = (1-ZTEMP1)/(ZTEMP2/(NV+1-JITER)-ZTEMP1)
            DO J=1,10
              ZQ = ZTEMP3**(JITER-1) - ZTEMP1*ZTEMP3 + ZTEMP1 - 1.0
              ZFA = (JITER-1)*ZTEMP3**(JITER-2) - ZTEMP1
              ZALPHA = ZTEMP3 - ZQ/ZFA
              ZTEMP3 = ZALPHA
            ENDDO
            IF (ZALPHA.LE.0.) THEN
               ZD1    = CS(NPSI1) - CS(NPSI)
               ZALPHA = ((REXT-1.)/NV/ZD1)**(1./(JITER-2))
            ENDIF
            
            WRITE(*,*) 'ZALPHA=',ZALPHA

C           COMPUTE VACUUM MESH-S
            ZD1    = CS(NPSI1) - CS(NPSI)
            CSV(1) = 1.
            CSV(2) = CSV(1) + ZD1

            DO J=3,JITER
               ZTMP3 = ZD1*ZALPHA**(J-2)
               CSV(J) = CSV(J-1) + ZTMP3
            ENDDO
      
            ZD1 = (REXT-CSV(JITER))/(NV+1-JITER)
            DO J=JITER+1,NV+1
              CSV(J) = CSV(JITER) + ZD1*(J-JITER)
            ENDDO

C           DEFINE WALL AND COIL POSITIONS
            NW(1) = 1
            NW(NWBPS) = NV + 1
            IF (NVEXP.EQ.8) THEN
            DO J=2,NWBPS-1
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
            ENDDO
            ELSEIF (NVEXP.EQ.80) THEN
            DO J=3,3
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               NW(J-1) = NW(J) - 1
               NW(J+1) = NW(J) + 1
               CSV(NW(J-1)) = RW(J-1)
               CSV(NW(J+1)) = RW(J+1)
               WRITE(*,5829) RW(J-1),TMP1,TMP2,NW(J-1)
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
               WRITE(*,5829) RW(J+1),TMP1,TMP2,NW(J+1)
            ENDDO
            DO J=6,6
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               NW(J-1) = NW(J) - 1
               NW(J+1) = NW(J) + 1
               CSV(NW(J-1)) = RW(J-1)
               CSV(NW(J+1)) = RW(J+1)
               WRITE(*,5829) RW(J-1),TMP1,TMP2,NW(J-1)
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
               WRITE(*,5829) RW(J+1),TMP1,TMP2,NW(J+1)
            ENDDO
            DO J=8,NWBPS-1
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
            ENDDO
            ELSEIF (NVEXP.EQ.800) THEN
            DO J=4,4
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               NW(J-2) = NW(J) - 2
               NW(J-1) = NW(J) - 1
               NW(J+1) = NW(J) + 1
               CSV(NW(J-2)) = RW(J-2)
               CSV(NW(J-1)) = RW(J-1)
               CSV(NW(J+1)) = RW(J+1)
               WRITE(*,5829) RW(J-2),TMP1,TMP2,NW(J-2)
               WRITE(*,5829) RW(J-1),TMP1,TMP2,NW(J-1)
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
               WRITE(*,5829) RW(J+1),TMP1,TMP2,NW(J+1)
            ENDDO
            DO J=7,7
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               NW(J-1) = NW(J) - 1
               NW(J+1) = NW(J) + 1
               CSV(NW(J-1)) = RW(J-1)
               CSV(NW(J+1)) = RW(J+1)
               WRITE(*,5829) RW(J-1),TMP1,TMP2,NW(J-1)
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
               WRITE(*,5829) RW(J+1),TMP1,TMP2,NW(J+1)
            ENDDO
            DO J=9,NWBPS-1
               JRW = 1
               DO JV=1,NV
                  IF (CSV(JV).LE.RW(J).AND.CSV(JV+1).GT.RW(J)) THEN
                     JRW = JV
                  ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RW(J)-CSV(JRW).LT.CSV(JRW+1)-RW(J)) THEN
                  NW(J) = JRW
                  CSV(JRW) = RW(J)
               ELSE 
                  NW(J) = JRW + 1
                  CSV(JRW+1) = RW(J)
               ENDIF
               WRITE(*,5829) RW(J),TMP1,TMP2,NW(J)
            ENDDO
            ENDIF
            IF (NW(2).NE.JITER) WRITE(*,*) 'JITER,NW(2)=',JITER,NW(2)
            
            DO J=1,5
               IF (RCOIL(J).GT.1.0) THEN
               JRW = 1
               DO JV=1,NV
                 IF (CSV(JV).LE.RCOIL(J).AND.CSV(JV+1).GT.RCOIL(J)) THEN
                    JRW = JV
                 ENDIF
               ENDDO
               TMP1 = CSV(JRW)
               TMP2 = CSV(JRW+1)
               IF (RCOIL(J)-CSV(JRW).LT.CSV(JRW+1)-RCOIL(J)) THEN
                  NCOIL = JRW
                  CSV(JRW) = RCOIL(J)
               ELSE 
                  NCOIL = JRW + 1
                  CSV(JRW+1) = RCOIL(J)
               ENDIF
               WRITE(*,5829) RCOIL(J),TMP1,TMP2,NCOIL
               ENDIF 
            ENDDO
         ELSE IF (NV .EQ. 1) THEN
C
            CSV(1) = 1.0
            CSV(NV+1) = REXT
C
         ELSE
            STOP 'INVALID NVEXP OR NV'
         ENDIF
C
         WRITE(*,*) 'VACUUM MESH CSV = '
         DO JV=1,NV+1
            WRITE(*,5830) JV,CSV(JV)
         ENDDO

         DO 505 J505=1,NV
 505        CSMV(J505) = .5 * (CSV(J505+1) + CSV(J505))
C

 5829    FORMAT ('RW=',E12.4,'  CSV=',2E12.4,'  NW=',I3)
 5830    FORMAT (I3,2X,E15.8)

         RETURN
C
 600     CONTINUE
C
***********************************************************************
*                                                                     *
* 5. S-MESH FOR BALLOONING OPTIMIZATION                               *
*                                                                     *
***********************************************************************
C
         ZDPS = 1. / DFLOAT(NPPR)
C
         DO 601 J601=1,NPPR+1
 601        PCS(J601) = (J601 - 1)  * ZDPS
C
         PCS(NPPR+1) = 1.
C
***********************************************************************
*                                                                     *
*        NMESHB = 0 ===> EQUIDISTANT S-MESH                           *
*        NMESHB = 1 ===> WEIGHTED S-MESH                              *
*        NPOIDB = 0 ===> NO WEIGHTING IS POSSIBLE FOR S-MESH          *
*        SOLPDB = 1 ===> NO WEIGHTING IS POSSIBLE FOR S-MESH          *
*                                                                     *
***********************************************************************
C
         IF (NMESHB .NE. 0 .AND. NPOIDB .NE. 0 .AND. 
     +       SOLPDB .NE. 1.) THEN
            CALL PACKME(NPPR+1,NPOIDB,PCS,BPLACE,BWIDTH,SOLPDB)
         ENDIF
C
         DO 603 J603=1,NPPR
 603        PCSM(J603) = .5 * (PCS(J603+1) + PCS(J603))
C
         PCSM(NPPR+1) = 1.
C
         RETURN
C
***********************************************************************
*                                                                     *
* 7.1. CSPEN-MESH                                                     *
*                                                                     *
***********************************************************************
C
 700     CONTINUE
C
         CALL SCOPY(NPSI,CS(2),1,CSPEN(5),5)
         CALL GAUSS(NMGAUS,ZRAC,ZWGT)
C
         DO 702 J702=1,NPSI
            ZADD = CS(J702+1) + CS(J702)
            ZDIF = CS(J702+1) - CS(J702)
C
            I = (J702 - 1) * (NMGAUS + 1)
C
            DO 701 J701=1,NMGAUS
 701           CSPEN(I+J701) = .5 * (ZADD + ZDIF * ZRAC(J701))
 702     CONTINUE
C
***********************************************************************
*                                                                     *
* 7.2. FILL IN EQUIDISTANT CTPEN-MESH                                 *
*                                                                     *
***********************************************************************
C
         ZDT = 2. * CPI / DFLOAT(NCHI)
C
         DO 703 J703=1,NCHI1
 703        ZT(J703) = (J703 - 1.) * ZDT
C
***********************************************************************
*                                                                     *
*        NMESHE = 0 ===> EQUIDISTANT CTPEN-MESH                       *
*        NMESHE = 1 ===> WEIGHTED CTPEN-MESH                          *
*        NPOIDE = 0 ===> NO WEIGHTING IS POSSIBLE FOR CTPEN-MESH      *
*        SOLPDE = 1 ===> NO WEIGHTING IS POSSIBLE FOR CTPEN-MESH      *
*                                                                     *
***********************************************************************
C
         IF (NMESHE .NE. 0 .AND. NPOIDE .NE. 0 .AND. 
     +           SOLPDE .NE. 1.) THEN
C     
            DO 704 J704=1,NPOIDE
               IF (EPLACE(J704) .GT. CPI) THEN
                 PRINT *,' EPLACE(',J704,')=',EPLACE(J704),
     +             ' SHOULD BE IN [-PI,PI]'
                 STOP
               ENDIF
               IF (EPLACE(J704) .LT. 0.0)
     +                EPLACE(J704) = EPLACE(J704) + 2.*CPI

 704           CONTINUE
C     
            CALL PACKMEP(NCHI1,NPOIDE,ZT,EPLACE,EWIDTH,SOLPDE)
C     
            DO 705 J705=1,NPOIDE
               IF (EPLACE(J705) .GT. CPI)
     +                EPLACE(J705) = EPLACE(J705) - 2.*CPI
 705           CONTINUE
C
         ENDIF
C
***********************************************************************
*                                                                     *
* 7.3. CTPEN MESH WITH CONSTANT AREA BETWEEN CONSTANT THETA SURFACES  *
*                                                                     *
***********************************************************************
C
         IF (NDIFT .NE. 0) THEN
            CALL TETARE(ZT,NCHI)
         ENDIF
C
         CALL SCOPY(NCHI1,ZT,1,CTPEN,5)
C
         DO 707 J707=1,NCHI
            ZADD = ZT(J707+1) + ZT(J707)
            ZDIF = ZT(J707+1) - ZT(J707)
C
            I = (J707 - 1) * (NMGAUS + 1)
C
            DO 706 J706=1,NMGAUS
 706           CTPEN(I+J706+1) = .5 * (ZADD + ZDIF * ZRAC(J706))
 707     CONTINUE
C
         RETURN
C
 800     CONTINUE
C
***********************************************************************
*                                                                     *
* 8. FILL IN EQUIDISTANT CTXT-MESH                                    *
*                                                                     *
***********************************************************************
C
         ZDT = 2. * CPI / DFLOAT(NTNOVA)
C
         DO 803 J803=1,NTNOVA
 803        CTXT(J803) = (J803 - 1.) * ZDT
C
         RETURN
         END
C*DECK C2SA02
C*CALL PROCESS
         SUBROUTINE PACKME(KN,KPOID,PMESH,PPLACE,PWIDTH,PSOLPD)
C        =======================================================
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
         INCLUDE 'DECLAR.inc'
         PARAMETER (IM = 401)
C
         DIMENSION
     R   PMESH(KN),   PPLACE(KPOID),   PWIDTH(KPOID),   ZW(IM)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C
***********************************************************************
*     PMESH                                                           *
*      ^                                                              *
*    1-!          ### : LORENTZIANS                          ++       *
*      !                                                +++++         *
*      !          +++ : ZW(I)                       ++++              *
*      !                                         +++                  *
*      !                                        +                     *
*      !                                      ++                      *
*      !                                    ++                        *
*      !                               +++++                          *
*      !                          +++++                               *
*      !                     +++++   !                                *
*      !                +++++        !SLOPE PSOLPD                    *
*      !            ++++                                              *
*      !         +++                                                  *
*      !        +                                                     *
*      !       +###                              ###                  *
*      !      +# ! #                            # ! #                 *
*      !   +++#  !  #                          #  !  #                *
*      !+++ ##   !   ##                      ##   !   ##              *
*    0-!---------!--------------------------------!-----------!-> S'  *
*      0     PPLACE(1)                         PPLACE(2)      1       *
*                                                                     *
***********************************************************************
*                                                                     *
* 1. STEP FOR EQUIDISTANT S'-MESH                                     *
*                                                                     *
***********************************************************************
C
         ZM = 1. / DFLOAT(IM - 1)
C
***********************************************************************
*                                                                     *
* 2. FILL IN DENSITY FUNCTION                                         *
*                                                                     *
***********************************************************************
C
         DO 2 J2=1,IM
C
         ZS     = (J2 - 1) * ZM
         ZW(J2) = 0.
C
         DO 1 J1=1,KPOID
C
         ZW(J2) = ZW(J2) + ATAN((ZS - PPLACE(J1)) / PWIDTH(J1))
     +                   + ATAN(      PPLACE(J1)  / PWIDTH(J1))
C
    1    CONTINUE
    2    CONTINUE
C
***********************************************************************
*                                                                     *
* 3. NORMALIZE IT TO ONE                                              *
*                                                                     *
***********************************************************************
C
         ZC = (1 - PSOLPD) / ZW(IM)
C
         DO 3 J3=1,IM
C
         ZS     = (J3 - 1) * ZM
         ZW(J3) = ZS * PSOLPD + ZC * ZW(J3)
C
    3    CONTINUE
C
***********************************************************************
*                                                                     *
* 4. FIND MESH POSITIONS                                              *
*                                                                     *
***********************************************************************
C
         PMESH( 1) = 0.
         PMESH(KN) = 1.
C
         ZDP = 1. / DFLOAT(KN - 1)
         ZF  = ZDP
         I   = 1
C
         DO 5 J5=2,IM
C
    4    CONTINUE
C
         IF (ZW(J5) .LE. ZF) GOTO 5
C
         I        = I + 1
         ZS       = (J5 - 2) * ZM
         PMESH(I) = ZS + (ZF - ZW(J5-1)) * ZM / (ZW(J5) - ZW(J5-1))
         ZF       = ZF + ZDP
C
         GOTO 4
C
    5    CONTINUE
C
         RETURN
         END
C*DECK C2SA03
C*CALL PROCESS
         SUBROUTINE PSVOL
C        ################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SA03 : DENSIFY STABILITY S-MESH SO THAT THE STABILITY MESH IS     *
*          EQUIDISTANT IN RHO (SEE SECTION 6.4.5 IN PUBLICATION       *
*          AND TABLE 1 FOR THE DEFINITION OF RHO)                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   ZDVDS(NPPSI1), ZPAR(NPPSI1), ZPISO(NPPSI1),   ZS(NPPSI1),    
     R   ZU1(NPPSI1),   ZU2(NPPSI1),  ZVOL(NPPSI1),    ZW(NPPSI1)
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(CID0,2*NPISO)
         CALL VZERO(CIDR,2*NPISO)
         CALL VZERO(CIDQ,2*NPISO)
         CALL VZERO(CID2,2*NPISO)
         CALL VZERO(ZVOL,NPPSI1)
         CALL VZERO(ZU1,NPPSI1)
         CALL VZERO(ZU2,NPPSI1)
C
         DO 1 J1=1,NPSI1
C
         PSIISO(J1) = SPSIM * (1. - CSM(J1) * CSM(J1))
         ZPISO(J1)  = SPSIM * (1. - CS(J1) * CS(J1))
C
    1    CONTINUE
C
         CALL RMRAD(NPSI1,SPSIM,RC0P,RC0P,ZPAR,SIGMAP,TETMAP,NTP2)
C
         IP = ISRCHFGE(NPSI1,PSIISO,1,CPSICL(1))
C
         IF (IP.LT.1)     IP = 1
         IF (IP.GT.NPSI1) IP = NPSI1
C
         CALL ISOFIND(IP,NPSI1,SIGPSI,TETPSI,WGTPSI,SPSIM,RC0P)
C
         DO 2 J2=IP,NPSI1
C
         CALL CINT(J2,SIGPSI(1,J2),TETPSI(1,J2),WGTPSI(1,J2))
C
   2     CONTINUE
C
         IF (IP .GT. 1) THEN
C
            DO 3 J3=1,IP-1
C
            I = J3
C
            IF (J3 .EQ. 1) I = 2
C
            CIDR(J3) = FCCCC0(CIDR(I-1),CIDR(I),CIDR(I+1),CIDR(I+2),
     ,                        PSIISO(I-1),PSIISO(I),PSIISO(I+1),
     ,                        PSIISO(I+2),PSIISO(J3))
C
    3       CONTINUE
C
         ENDIF
C
         ZVOL(1) = 0.
C
         DO 4 J4=2,NPSI1
C
         ZVOL(J4) = ZVOL(J4-1) + CIDR(J4-1) * (ZPISO(J4)-ZPISO(J4-1))
C
    4    CONTINUE
C
         DO 5 J5=1,NPSI1
C
         ZVOL(J5) = SQRT(ZVOL(J5))
C
    5    CONTINUE
C
         DO 6 J6=1,NPSI
C
         ZDVDS(J6) = (ZVOL(J6+1) - ZVOL(J6)) / (CS(J6+1) - CS(J6))
C
    6    CONTINUE
C
         ZU1(1) = 0.
         ZU2(1) = 0.
C
         DO 8 J8=2,NPSI1
C
         ZU1(J8) = ZVOL(J8)
C
         IF (NPOIDA .EQ. 0 .OR. NMESHA .EQ. 0) THEN
C
            ZU2(J8) = ZVOL(J8)
C
         ELSE
C
            ZU2(J8) = ZU2(J8-1)
C
            DO 7 J7=1,NPOIDA
C
            ZU2(J8) = ZU2(J8) + ZDVDS(J8-1) * AWIDTH(J7) /
     /                (AWIDTH(J7)**2 + (APLACE(J7)-CSM(J8-1))**2) *
     *                (CS(J8) - CS(J8-1))
C
    7       CONTINUE
C
         ENDIF
C
    8    CONTINUE
C
***********************************************************************
*                                                                     *
* 3. NORMALIZE IT TO ONE                                              *
*                                                                     *
***********************************************************************
C
         DO 9 J9=1,NPSI1
C
         ZW(J9) = SOLPDA * ZU1(J9) / ZU1(NPSI1) +
     +            (1. - SOLPDA) * ZU2(J9) / ZU2(NPSI1)
C
    9    CONTINUE
C
***********************************************************************
*                                                                     *
* 4. FIND MESH POSITIONS                                              *
*                                                                     *
***********************************************************************
C
         ZS(1)     = 0.
         ZS(NPSI1) = 1.
C
         DO 10 J10=2,NPSI
C
         ZWS = DFLOAT(J10 - 1) / DFLOAT(NPSI)
C
         IS = ISRCHFGE(NPSI1,ZW,1,ZWS)
C
         IF (IS.LT.2)     IS = 2
         IF (IS.GT.NPSI1) IS = NPSI1
C
         ZS(J10) = CS(IS-1) + (CS(IS) - CS(IS-1)) * (ZWS - ZW(IS-1)) /
     /                        (ZW(IS) -  ZW(IS-1))
C
   10    CONTINUE
C
         CALL SCOPY(NPSI1,ZS,1,CS,1)
C
         DO 11 J11=1,NPSI
C
         CSM(J11) = .5 * (CS(J11+1) + CS(J11))
C
  11     CONTINUE
C
         CSM(NPSI1) = 1.
C
         RETURN
         END
C*DECK C2SA04
C*CALL PROCESS
         SUBROUTINE TETARE(PT,KN)
C        ########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SA04 : DENSIFY EQUILIBRIUM THETA-MESH AUTOMATICALLY               *
*          NDIFT= 1: CONSTANT POLOIDAL FLUX BETWEEN 2 SUCESSIVE ANGLES*
*          NDIFT= 2: CONSTANT ARC LENGHT AT PLASMA SURFACE FOR EVERY  *
*                    THETA INTERVAL                                   *
*          (SEE SECTION 6.4.5 IN PUBLICATION)                         *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
C
         PARAMETER (IM = 301)
C
         DIMENSION
     R   PT(*),
     R   ZBND(IM),   ZBND1(IM),   ZBND2(IM),   ZDRSDT(IM),   ZTET(IM),   
     R   ZTET1(IM),  ZTET2(IM),   ZU1(IM),     ZU2(IM),      ZW(IM)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(ZU1,IM)
         CALL VZERO(ZU2,IM)
C
         ZEPS = 1.E-4
         ZDT  = 2.*CPI / DFLOAT(IM - 1)
C
         DO 1 J1=1,IM
C
         ZTET(J1)  = (J1 - 1) * ZDT
         ZTET1(J1) = ZTET(J1) - ZEPS
         ZTET2(J1) = ZTET(J1) + ZEPS
C
    1    CONTINUE
C
         CALL BOUND(IM,ZTET,ZBND)
         CALL BOUND(IM,ZTET1,ZBND1)
         CALL BOUND(IM,ZTET2,ZBND2)
C
         ZU1(1) = 0.
         ZU2(1) = 0.
C
         IF (NDIFT .EQ. 1) THEN
C
         DO 3 J3=2,IM
C
         ZU1(J3) = ZU1(J3-1) + .25 * (ZTET(J3) - ZTET(J3-1)) *
     *             (ZBND(J3)**2 + ZBND(J3-1)**2)
C
         IF (NPOIDD .EQ. 0 .OR. NMESHD .EQ. 0) THEN
C
            ZU2(J3) = ZU1(J3)
C
         ELSE
C
            ZU2(J3) = ZU2(J3-1)
C
            DO 2 J2=1,NPOIDD
C
            ZU2(J3)= ZU2(J3) + .5 * (ZTET(J3) - ZTET(J3-1)) *
     *               (ZBND(J3-1)**2 * DWIDTH(J2) / 
     +               (DWIDTH(J2)**2+SIN((DPLACE(J2)-ZTET(J3-1))/2.)**2)+
     *               ZBND(J3)**2 * DWIDTH(J2) / 
     +               (DWIDTH(J2)**2+SIN((DPLACE(J2)-ZTET(J3))/2.)**2))
C
    2       CONTINUE
C
         ENDIF
C
    3    CONTINUE
C
         ELSE IF (NDIFT .EQ. 2) THEN
C
         DO 4 J4=1,IM
C
         ZDRSDT(J4) = .5 *(ZBND2(J4) - ZBND1(J4)) / ZEPS
C
    4    CONTINUE
C
         DO 6 J6=2,IM
C
         ZU1(J6) = ZU1(J6-1) + .5 * (ZTET(J6) - ZTET(J6-1)) *
     *             (SQRT(ZBND(J6)**2 + ZDRSDT(J6)**2) +
     +              SQRT(ZBND(J6-1)**2 + ZDRSDT(J6-1)**2))
C
         IF (NPOIDD .EQ. 0 .OR. NMESHD .EQ. 0) THEN
C
            ZU2(J6) = ZU1(J6)
C
         ELSE
C
            ZU2(J6) = ZU2(J6-1)
C
            DO 5 J5=1,NPOIDD
C
            ZU2(J6)= ZU2(J6) + .5 * (ZTET(J6) - ZTET(J6-1)) *
     *               (ZBND(J6-1) * DWIDTH(J5) / 
     /               (DWIDTH(J5)**2+SIN((DPLACE(J5)-ZTET(J6-1))/2.)**2)+
     *               ZBND(J6) * DWIDTH(J5) / 
     /               (DWIDTH(J5)**2+SIN((DPLACE(J5)-ZTET(J6))/2.)**2))
C
    5       CONTINUE
C
         ENDIF
C
    6    CONTINUE
C
         ENDIF
C
***********************************************************************
*                                                                     *
* 3. NORMALIZE IT TO ONE                                              *
*                                                                     *
***********************************************************************
C
         DO 7 J7=1,IM
C
         ZW(J7) = 2. * CPI * (SOLPDD * ZU1(J7) / ZU1(IM) +
     +                        (1. - SOLPDD) * ZU2(J7) / ZU2(IM))
C
    7    CONTINUE
C
***********************************************************************
*                                                                     *
* 4. FIND MESH POSITIONS                                              *
*                                                                     *
***********************************************************************
C
         PT(1)   = 0.
         PT(KN+1) = 2.*CPI
C
         ZDP = 2.*CPI / DFLOAT(KN)
         ZF  = ZDP
         I   = 1
C
         DO 9 J9=2,IM
C
    8    CONTINUE
C
         IF (ZW(J9) .LE. ZF) GOTO 9
C
         I     = I + 1
         ZT    = (J9 - 2) * ZDT
         PT(I) = ZT + (ZF - ZW(J9-1)) * ZDT / (ZW(J9) - ZW(J9-1))
         ZF    = ZF + ZDP
C
         GOTO 8
C
    9    CONTINUE
C
         RETURN
         END
C*DECK C2SA05
C*CALL PROCESS
         SUBROUTINE QPLACS(CSOLD,ZQOLD)
C        ##############################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SA05 : DENSIFY STABILITY MESH AT A PREDEFINED SET OF Q-VALUES     *
*          GIVEN IN QPLACS (SEE SECTION 6.4.5 IN PUBLICATION)         *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R   ZQ(NPPSI1),   ZS(NPPSI1), zder(nppsi1),
     R   CSOLD(*),     ZQOLD(*)
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         DO 1 J1=1,NPSI
C
         CSM(J1) = .5 * (CS(J1+1) + CS(J1))
C
    1    CONTINUE
C
         CSM(NPSI1) = 1.
C
         CALL PREMAP(3)
C
         IF (NRFP .EQ. 0) THEN
C
***********************************************************************
*                                                                     *
* TOKAMAK EQUILIBRIUM                                                 *
*                                                                     *
***********************************************************************
C
         IF (NCSCAL .EQ. 1) THEN
C
            DO 2 J2=1,NPSI1
C
            ZS(J2) = SQRT(1. - PSIISO(J2) / SPSIM)
            ZQ(J2) = .5 * TMF(J2) * CIDQ(J2) / CPI
C
    2       CONTINUE
cab
      IF (NQMIN .LE. 0) GOTO 103
cab   csspec is redefined only if nqmin > 0
cab
      do 101 j = npsi,2,-1
      zder(j) = tmf(j)*(tmf(j+1)-tmf(j-1))*cidq(j)**2
     &        + (2.*cpi*qspec)**2*(cidq(j+1)-cidq(j-1))/cidq(j)
      if (zder(j).lt.0.) goto 102
 101  continue
      csspec=0.
      goto 103
 102  csspec=zs(j)-zder(j)*(zs(j+1)-zs(j))/(zder(j+1)-zder(j))
 103  WRITE(*,*) ' NQMIN = ',NQMIN,'   CSSPEC =',CSSPEC
cab
C
            ICS = ISRCHFGE(NPSI1,ZS,1,CSSPEC) - 1
C
            IF (ICS .LT. 2)         ICS = 2
            IF (ICS .GE. NPSI1 - 1) ICS = NPSI1 - 2
C
            QICS = FCCCC0(ZQ(ICS-1),ZQ(ICS),ZQ(ICS+1),ZQ(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS),ZS(ICS+1),ZS(ICS+2),CSSPEC)
            TICS = FCCCC0(TMF(ICS-1),TMF(ICS),TMF(ICS+1),TMF(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS), ZS(ICS+1),ZS(ICS+2),CSSPEC)
C
            ZCSHFT = TICS**2 * ((QSPEC / QICS)**2 - 1.)
C
            DO 3 J3=1,NPSI1
C
            ZQ(J3) = ZQ(J3) * SQRT(1. + ZCSHFT / TMF(J3)**2)
C
    3       CONTINUE
C
         ELSE IF (NCSCAL .EQ. 2) THEN
C
            T0 = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                  PSIISO(1),PSIISO(2),PSIISO(3),PSIISO(4),SPSIM)
C
            DO 4 J4=1,NPSI1
C
            ZQ(J4) = .5 * TMF(J4) * CIDQ(J4) / CPI
C
    4       CONTINUE
C
            IF (NSURF .EQ. 1) THEN
C
               ZCSCAL = 1.
C
            ELSE
C
               ZCSCAL = CURRT / CUROLD
C
            ENDIF
C
            IF (NTMF0 .EQ. 0) THEN
C
               ZCSHFT = 1. - (ZCSCAL * TMF(NPSI1))**2
C
            ELSE IF (NTMF0 .EQ. 1) THEN
C
               ZCSHFT = 1. - (ZCSCAL * T0)**2
C
            ENDIF
C
            DO 5 J5=1,NPSI1
C 
            ZQ(J5) = ZQ(J5) * SQRT(1. + ZCSHFT / (ZCSCAL*TMF(J5))**2)
C
    5       CONTINUE
C
         ELSE IF (NCSCAL .EQ. 3) THEN
C
            ZS(1) = CSM(1)**2 * ABS(SPSIM) * CIDR(1)
            ZQ(1) = .5 * TMF(1) * CIDQ(1) / CPI
C
            DO 7 J7=2,NPSI1
C
            ZS(J7) = ZS(J7-1) + ABS(SPSIM) * (CSM(J7) - CSM(J7-1)) *
     *               (CSM(J7) * CIDR(J7) + CSM(J7-1) * CIDR(J7-1))
            ZQ(J7) = .5 * TMF(J7) * CIDQ(J7) / CPI
C
    7       CONTINUE
C
            DO 8 J8=1,NPSI1
C
            ZS(J8) = SQRT(ZS(J8) / ZS(NPSI1))
C
    8       CONTINUE
C
            ICS = ISRCHFGE(NPSI1,ZS,1,CSSPEC) - 1
C
            IF (ICS .LT. 2)      ICS = 2
            IF (ICS .GE. NPSI1 - 1) ICS = NPSI1 - 2
C
            QICS = FCCCC0(ZQ(ICS-1),ZQ(ICS),ZQ(ICS+1),ZQ(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS),ZS(ICS+1),ZS(ICS+2),CSSPEC)
            TICS = FCCCC0(TMF(ICS-1),TMF(ICS),TMF(ICS+1),TMF(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS), ZS(ICS+1),ZS(ICS+2),CSSPEC)
C
            ZCSHFT = TICS**2 * ((QSPEC / QICS)**2 - 1.)
C
            DO 9 J9=1,NPSI1
C
            ZQ(J9) = ZQ(J9) * SQRT(1. + ZCSHFT / TMF(J9)**2)
C
    9       CONTINUE
C
         ELSE IF (NCSCAL .EQ. 4) THEN
C
            DO 10 J10=1,NPSI1
C
            ZQ(J10) = .5 * TMF(J10) * CIDQ(J10) / CPI
C
   10       CONTINUE
C
         ENDIF
C
***********************************************************************
*                                                                     *
* REVERSED FIELD PINCH EQUILIBRIUM                                    *
*                                                                     *
***********************************************************************
C
         ELSE IF (NRFP .EQ. 1) THEN
C
            DO 11 J11=1,NPSI1
C
            ZQ(J11) = .5 * TMF(J11) * CIDQ(J11) / CPI
C
   11       CONTINUE
C
         ENDIF

         DO I=1,NPSI1
            CSOLD(I)=CSM(I)
            ZQOLD(I)=ZQ(I)
         ENDDO

C YQL, 2005-03-27
C AUTOMATIC PACKING
         IF (NMESHA.GE.3) CALL QAUTOPACKA(ZQ)

         JPLACE = 0
C
         DO 13 J13=1,NPOIDQ
C
         DO 12 J12=1,NPSI
C
         IF ((QPLACE(J13)-ZQ(J12))*(QPLACE(J13)-ZQ(J12+1)) .LE. 0.) THEN
C
            JPLACE = JPLACE + 1
C
            APLACE(JPLACE) = CSM(J12) + (CSM(J12+1) - CSM(J12)) *
     *                       (QPLACE(J13)-ZQ(J12)) / (ZQ(J12+1)-ZQ(J12))
            AWIDTH(JPLACE) = QWIDTH(J13)
CLIU Y.Q.   -----------------------------------------------------------
C           IF (APLACE(JPLACE).LT.0.5) JPLACE = JPLACE - 1
CLIU Y.Q.   -----------------------------------------------------------
C
         ENDIF
C
   12    CONTINUE
   13    CONTINUE
C
         IF (NMESHA.EQ.5.OR.NMESHA.EQ.6) THEN
           JPLACE = JPLACE + 1
           APLACE(JPLACE) = CSM(NPSI1)
           IF (JPLACE.GT.1) THEN
              AWIDTH(JPLACE) = AWIDTH(JPLACE-1)
           ELSE
              AWIDTH(JPLACE) = 1.0/NPSI*QWIDTH0
           ENDIF
         ENDIF  

         NPOIDA = JPLACE
C
         WRITE(6,'(/,A20,/,A20,I2,A,/,(1P10E12.4))') ' AFTER Q PACKING:'
     +     ,'APLACE(1:',NPOIDA,')',(APLACE(J),J=1,NPOIDA)
         WRITE(6,'(A20,I2,A,/,(1P10E12.4))')
     +     'AWIDTH(1:',NPOIDA,')',(AWIDTH(J),J=1,NPOIDA)
C
         RETURN
         END


         SUBROUTINE QAUTOPACKA(ZQ)
C        #########################
C
C                                        AUTHORS:
C                                        Y.Q. LIU
***********************************************************************
*                                                                     *
* AUTOMATIC Q-PACKING WITH PLASMA ROTATION                            *     
*                                                                     *
***********************************************************************
C
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
C
         INTEGER    CHROTD
         PARAMETER  (CHROTD=26)
         REAL       TEMP
         DIMENSION
     R   ZQ(*),          
     R   CSN(NPISO),  RQ(NPISO),
     R   SROT(NPISO), ZROT(NPISO), RROT(NPISO), 
     R   SDEN(NPISO), ZDEN(NPISO), RDEN(NPISO),
     R   ZWORK(NPISO), ZWORK1(NPISO),
     R   D2RALL(NPISO),
     I   IC(NPISO),    IT(NPISO) 
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C-----------------------------------------------------------------------
C FIND MMIN AND MMAX FROM Q
C--------------------------
       QMIN=ZQ(1)
       QMAX=ZQ(NPSI1)
       DO I=1,NPSI1
          IF (ZQ(I).LT.QMIN) QMIN=ZQ(I)
          IF (ZQ(I).GT.QMAX) QMAX=ZQ(I)
       ENDDO
       MMIN=INT(NTOR*QMIN)+1
       MMAX=INT(NTOR*QMAX)
       WRITE(*,*) 
       WRITE(*,*) 'QAUTOPACKA:'

C-----------------------------------------------------------------------
C FIND RATIONAL QPLACE
C---------------------
       IF (NMESHA.EQ.3.OR.NMESHA.EQ.5) THEN
          NPOIDQ=MMAX-MMIN+1
          IF (NPOIDQ.GT.NPPACK) THEN
             WRITE(*,*) 'TOO MANY PACKING SURFACES, NPOIDQ=',NPOIDQ
             NPOIDQ=NPPACK
             MMAX=NPOIDQ+MMIN-1
          ENDIF
          DO I=MMIN,MMAX
             QPLACE(I-MMIN+1)=DFLOAT(I)/NTOR
          ENDDO 
       ENDIF

       IF (NMESHA.EQ.4.OR.NMESHA.EQ.6) THEN
C-----------------------------------------------------------------------
C READ PRORD
C-----------
       OPEN(CHROTD,FILE='PROFROT',FORM='FORMATTED')
       READ (CHROTD,*)NROT,I
       IF (NROT.GE.1.AND.NROT.LE.NPPSI1) THEN
          DO I=1,NROT
             READ (CHROTD,*)SROT(I),ZROT(I)
          ENDDO
          TEMP = ABS(ZROT(1))
          DO I=2,NROT
            IF (TEMP.LT.ABS(ZROT(I))) TEMP=ABS(ZROT(I))
          ENDDO
          DO I=1,NROT
             ZROT(I) = ZROT(I)/TEMP
          ENDDO
       ELSE
          WRITE(*,*) 'NROT=',NROT,' OUT OF [1,',NPPSI1,']'
          STOP 'NROT'
       ENDIF
       CLOSE(CHROTD)
       
       OPEN(CHROTD,FILE='PROFDEN',FORM='FORMATTED')
       READ (CHROTD,*)NDEN,I
       IF (NDEN.GE.1.AND.NDEN.LE.NPPSI1) THEN
          DO I=1,NDEN
             READ (CHROTD,*)SDEN(I),ZDEN(I)
          ENDDO
          TEMP = ZDEN(1)
          DO I=1,NDEN
             ZDEN(I) = ZDEN(I)/TEMP
          ENDDO
       ELSE
          WRITE(*,*) 'NDEN=',NDEN,' OUT OF [1,',NPPSI1,']'
          STOP 'NDEN'
       ENDIF
       CLOSE(CHROTD)

C-----------------------------------------------------------------------
C NEW MESH DENSER AT HIGH Q
C--------------------------
       QSUM=0.0
       DO I=1,NPSI
          QSUM=QSUM+1.0/ZQ(I)**4
       ENDDO
       H0=1.0/QSUM
       CSN(1)=0.0
       DO I=1,NPSI
          CSN(I+1)=CSN(I) + H0/ZQ(I)**4
       ENDDO

C-----------------------------------------------------------------------
C SPLINE PRORD
C-------------
      CALL SPLINE(CSM,ZQ,NPSI1,D2RALL,ZWORK,ZWORK1)
      CALL RESETI(IC,NPSI1,1)
      DO JG=1,NPSI1
         DO JT=1,NPSI1
            IF (IC(JG).EQ.0) GOTO 1
            IT(JG) = JT-1
            IF (CSM(JT).GE.CSN(JG)) IC(JG) = 0
 1          CONTINUE
         ENDDO
      ENDDO
      DO J=1,NPSI1
         ICHI = IT(J)
         IF (ICHI.LT.1) ICHI = 1
         IF (ICHI.GT.NPSI1-1) ICHI = NPSI1-1
         ZH = CSM(ICHI+1) - CSM(ICHI)
         ZA = (CSM(ICHI+1) - CSN(J)) / ZH
         ZB = (CSN(J) - CSM(ICHI)) / ZH
         ZC = (ZA + 1)*(ZA - 1)*ZH*(CSM(ICHI+1)-CSN(J))/6.
         ZD = (ZB + 1)*(ZB - 1)*ZH*(CSN(J)-CSM(ICHI))/6.
         RQ(J) = ZA*ZQ(ICHI) + ZB*ZQ(ICHI+1) + 
     &           ZC*D2RALL(ICHI) + ZD*D2RALL(ICHI+1)
      ENDDO

      CALL SPLINE(SROT,ZROT,NROT,D2RALL,ZWORK,ZWORK1)
      CALL RESETI(IC,NPSI1,1)
      DO JG=1,NPSI1
         DO JT=1,NROT
            IF (IC(JG).EQ.0) GOTO 2
            IT(JG) = JT-1
            IF (SROT(JT).GE.CSN(JG)) IC(JG) = 0
 2          CONTINUE
         ENDDO
      ENDDO
      DO J=1,NPSI1
         ICHI = IT(J)
         IF (ICHI.LT.1) ICHI = 1
         IF (ICHI.GT.NROT-1) ICHI = NROT-1
         ZH = SROT(ICHI+1) - SROT(ICHI)
         ZA = (SROT(ICHI+1) - CSN(J)) / ZH
         ZB = (CSN(J) - SROT(ICHI)) / ZH
         ZC = (ZA + 1)*(ZA - 1)*ZH*(SROT(ICHI+1)-CSN(J))/6.
         ZD = (ZB + 1)*(ZB - 1)*ZH*(CSN(J)-SROT(ICHI))/6.
         RROT(J) = ZA*ZROT(ICHI) + ZB*ZROT(ICHI+1) + 
     &             ZC*D2RALL(ICHI) + ZD*D2RALL(ICHI+1)
      ENDDO

      CALL SPLINE(SDEN,ZDEN,NDEN,D2RALL,ZWORK,ZWORK1)
      CALL RESETI(IC,NPSI1,1)
      DO JG=1,NPSI1
         DO JT=1,NDEN
            IF (IC(JG).EQ.0) GOTO 3
            IT(JG) = JT-1
            IF (SDEN(JT).GE.CSN(JG)) IC(JG) = 0
 3          CONTINUE
         ENDDO
      ENDDO
      DO J=1,NPSI1
         ICHI = IT(J)
         IF (ICHI.LT.1) ICHI = 1
         IF (ICHI.GT.NDEN-1) ICHI = NDEN-1
         ZH = SDEN(ICHI+1) - SDEN(ICHI)
         ZA = (SDEN(ICHI+1) - CSN(J)) / ZH
         ZB = (CSN(J) - SDEN(ICHI)) / ZH
         ZC = (ZA + 1)*(ZA - 1)*ZH*(SDEN(ICHI+1)-CSN(J))/6.
         ZD = (ZB + 1)*(ZB - 1)*ZH*(CSN(J)-SDEN(ICHI))/6.
         RDEN(J) = ZA*ZDEN(ICHI) + ZB*ZDEN(ICHI+1) + 
     &             ZC*D2RALL(ICHI) + ZD*D2RALL(ICHI+1)
      ENDDO
      IF (RDEN(NPSI1).LE.1.0E-15) RDEN(NPSI1)=RDEN(NPSI)*0.1

C-----------------------------------------------------------------------
C RESCALE ROTATION
C-----------------
       DO I=1,NPSI1
          RROT(I)=ROTE*ABS(RROT(I))
       ENDDO
       
C       do i=1,npsi1
C          write(*,*) csn(i),rq(i),rrot(i),rden(i)
C       enddo

C-----------------------------------------------------------------------
C FIND QPLACE FOR GIVEN ROTATION
C-------------------------------
       JPLACE=0
       DO M=MMIN,MMAX
          Q=DFLOAT(M)/NTOR
          DO I=1,NPSI
             IF ((RQ(I)-Q)*(RQ(I+1)-Q).LE.0.0) THEN
                ROT=RROT(I) + (RROT(I+1)-RROT(I))*
     &                        (Q-RQ(I))/(RQ(I+1)-RQ(I))
                ROT=ROT*NTOR
                WRITE(*,*) 'NTOR=',NTOR,' M=',M,' ROT=',ROT
                DO J=1,NPSI
                   FPS1=1.0 + (RQ(J)/(M-1-NTOR*RQ(J)))**2 +
     &                      (RQ(J)/(M+1-NTOR*RQ(J)))**2
                   FPS2=1.0 + (RQ(J+1)/(M-1-NTOR*RQ(J+1)))**2 +
     &                      (RQ(J+1)/(M+1-NTOR*RQ(J+1)))**2
                   F1=(DFLOAT(M)/RQ(J)-NTOR)/SQRT(RDEN(J)*FPS1)
                   F2=(DFLOAT(M)/RQ(J+1)-NTOR)/SQRT(RDEN(J+1)*FPS2)
                   IF ((F1-ROT)*(F2-ROT).LE.0.0) THEN
                      QQ=RQ(J) + (RQ(J+1)-RQ(J))*(ROT-F1)/(F2-F1)
                      IF (QQ.GT.Q-0.5/NTOR.AND.QQ.LT.Q+0.5/NTOR) THEN
                      JPLACE=JPLACE + 1
                      IF (JPLACE.GT.NPPACK) THEN
                         WRITE(*,*) 'TOO MANY PACKING SURFACES, 
     &                               JPLACE=',JPLACE
                         STOP 'JPLACE'
                      ENDIF
                      QPLACE(JPLACE)=QQ
                      ENDIF
                   ENDIF
                   IF ((F1+ROT)*(F2+ROT).LE.0.0) THEN
                      QQ=RQ(J) + (RQ(J+1)-RQ(J))*(-ROT-F1)/(F2-F1)
                      IF (QQ.GT.Q-0.5/NTOR.AND.QQ.LT.Q+0.5/NTOR) THEN
                      JPLACE=JPLACE + 1
                      IF (JPLACE.GT.NPPACK) THEN
                         WRITE(*,*) 'TOO MANY PACKING SURFACES, 
     &                               JPLACE=',JPLACE
                         STOP 'JPLACE'
                      ENDIF
                      QPLACE(JPLACE)=QQ
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       NPOIDQ=JPLACE

       ENDIF

C-----------------------------------------------------------------------
C DETERMINE QWIDTH BY Q-PROFILE
C------------------------------
       DO J=1,NPOIDQ
          QWIDTH(J)=1.0/NPSI/SQRT(QPLACE(J))*QWIDTH0
       ENDDO

C-----------------------------------------------------------------------
C OUTPUT RESULTS
C---------------
       WRITE(*,1000) (QPLACE(J),J=1,NPOIDQ)
       WRITE(*,1100) (QWIDTH(J),J=1,NPOIDQ)
 1000  FORMAT(' QPLACE=',40(F7.4,','))
 1100  FORMAT(' QWIDTH=',40(F7.4,','))
          
       RETURN
       END

         SUBROUTINE QAUTOPACKB(CSOLD,ZQOLD)
C        ##################################
C
C                                        AUTHORS:
C                                        Y.Q. LIU
***********************************************************************
*                                                                     *
* AUTOMATIC SHIFTING OF MESH POINTS CLOSEST TO Q-RESONANT SURFACES    *
*                                                                     *
***********************************************************************
C
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     R   CSOLD(*),    ZQOLD(*),
     R   ZS(NPISO),   RQ(NPISO),
     R   ZWORK(NPISO), ZWORK1(NPISO),
     R   D2RALL(NPISO),
     I   IC(NPISO),    IT(NPISO) 
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
      WRITE(*,*)
      WRITE(*,*) 'QAUTOPACKB:'
C-----------------------------------------------------------------------
C SPLINE ZQOLD ON NEW MESH CS
C----------------------------
      CALL SPLINE(CSOLD,ZQOLD,NPSI1,D2RALL,ZWORK,ZWORK1)
      CALL RESETI(IC,NPSI1,1)
      DO JG=1,NPSI1
         DO JT=1,NPSI1
            IF (IC(JG).EQ.0) GOTO 1
            IT(JG) = JT-1
            IF (CSOLD(JT).GE.CS(JG)) IC(JG) = 0
 1          CONTINUE
         ENDDO
      ENDDO
      DO J=1,NPSI1
         ICHI = IT(J)
         IF (ICHI.LT.1) ICHI = 1
         IF (ICHI.GT.NPSI1-1) ICHI = NPSI1-1
         ZH = CSOLD(ICHI+1) - CSOLD(ICHI)
         ZA = (CSOLD(ICHI+1) - CS(J)) / ZH
         ZB = (CS(J) - CSOLD(ICHI)) / ZH
         ZC = (ZA + 1)*(ZA - 1)*ZH*(CSOLD(ICHI+1)-CS(J))/6.
         ZD = (ZB + 1)*(ZB - 1)*ZH*(CS(J)-CSOLD(ICHI))/6.
         RQ(J) = ZA*ZQOLD(ICHI) + ZB*ZQOLD(ICHI+1) + 
     &           ZC*D2RALL(ICHI) + ZD*D2RALL(ICHI+1)
      ENDDO

C-----------------------------------------------------------------------
C SHIFT CS ACCORDING TO RQ
C-------------------------
      CALL SCOPY(NPSI1,CS,1,ZS,1)
      DO I=1,NPOIDQ
         Q=QPLACE(I)
         DO J=1,NPSI
            IF ((RQ(J)-Q)*(RQ(J+1)-Q).LE.0.0) THEN
               SNEW=(Q*(CS(J+1)-CS(J)) - 
     &               (CS(J+1)*RQ(J)-CS(J)*RQ(J+1)))/(RQ(J+1)-RQ(J))
               DELTAS=SNEW - 0.5*(CS(J)+CS(J+1))
               write(*,*) 'CS=',CS(J),CS(J+1)
               write(*,*) ' Q=',RQ(J),RQ(J+1)
               write(*,*) 'QS=',Q,' DELTAS=',DELTAS
               IF ((CS(J)+DELTAS).GT.ZS(J-1).AND.
     &             (CS(J+1)+DELTAS).LT.ZS(J+2)) THEN
                  ZS(J)=CS(J) + DELTAS
                  ZS(J+1)=CS(J+1) + DELTAS
               ELSE
                  WRITE(*,*) 'WARNING: NOT POSSIBLE TO SHIFT 
     &                        MESH AROUND Q=',Q
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      CALL SCOPY(NPSI1,ZS,1,CS,1)

      DO J=1,NPSI
         CSM(J)=0.5*(CS(J)+CS(J+1))
      ENDDO
      CSM(NPSI1)=1.0
          
       RETURN
       END

C*DECK C2SA06
C*CALL PROCESS
         SUBROUTINE PACKMEP(KN,KPOID,PMESH,PPLACE,PWIDTH,PSOLPD)
C        =======================================================
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMSOL.inc'
         PARAMETER (IM = 401)
C
         DIMENSION
     R   PMESH(KN),   PPLACE(KPOID),   PWIDTH(KPOID),   ZW(IM)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
***********************************************************************
*                                                                     *
* 1. STEP FOR EQUIDISTANT THETA'-MESH                                 *
*                                                                     *
***********************************************************************
C
         ZM = 2. * CPI / DFLOAT(IM - 1)
C
***********************************************************************
*                                                                     *
* 2. FILL IN DENSITY FUNCTION                                         *
*                                                                     *
***********************************************************************
C     
         DO 2 J2=1,IM
C
         ZS     = (J2 - 1) * ZM
         ZW(J2) = 0.
C
         DO 1 J1=1,KPOID
C
         ZZ = SQRT(PWIDTH(J1)**2+1)
C
         ZW(J2) = ZW(J2) + 2./(ZZ*PWIDTH(J1))*(
     &            ATAN2(PWIDTH(J1)*TAN(.25*(ZS-PPLACE(J1))),ZZ+1)+
     &            ATAN2(PWIDTH(J1)*TAN(.25*(ZS-PPLACE(J1))),ZZ-1)+
     &            ATAN2(PWIDTH(J1)*TAN(.25*(   PPLACE(J1))),ZZ+1)+
     &            ATAN2(PWIDTH(J1)*TAN(.25*(   PPLACE(J1))),ZZ-1))
C
    1    CONTINUE
    2    CONTINUE
C     
***********************************************************************
*                                                                     *
* 3. NORMALIZE IT TO ONE                                              *
*                                                                     *
***********************************************************************
C
         ZC = 2.* CPI * (1 - PSOLPD) / ZW(IM)
C
         DO 3 J3=1,IM
C
         ZS     = (J3 - 1) * ZM
         ZW(J3) = ZS * PSOLPD + ZC * ZW(J3)
C
    3    CONTINUE
C
***********************************************************************
*                                                                     *
* 4. FIND MESH POSITIONS                                              *
*                                                                     *
***********************************************************************
C
         PMESH( 1) = 0.
         PMESH(KN) = 2.*CPI
C
         ZDP = 2.*CPI / DFLOAT(KN - 1)
         ZF  = ZDP
         I   = 1
C
         DO 5 J5=2,IM
C
    4    CONTINUE
C
         IF (ZW(J5) .LE. ZF) GOTO 5
C
         I        = I + 1
         ZS       = (J5 - 2) * ZM
         PMESH(I) = ZS + (ZF - ZW(J5-1)) * ZM / (ZW(J5) - ZW(J5-1))
         ZF       = ZF + ZDP
C
         GOTO 4
C
    5    CONTINUE
C
         RETURN
         END
C*DECK C2SB01
C*CALL PROCESS
         SUBROUTINE GUESS(KGUESS)
C        ########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SB01 INITIALIZE PICARD ITERATION:                                 *
*        KGUESS = 1:  PSI(I,J) = -.1 * (1 - CSIG(I,J)**2)             *
*        KGUESS = 2:  INTEPOLATE SOLUTION ON PREVIOUS DISCRETIZATION  *
*                     MESH                                            *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   IC(NPT+NPPSI1), IS0(NPT+NPPSI1), IT0(NPT+NPPSI1),
     R   ZBND(NPT),  ZPCEL(NPT,16),  ZF(NPT,16),    
     R   ZS(NPT),    ZT(NPT),        ZSIGMA(NPT),
     R   ZT1(NPT),   ZT2(NPT),       ZS1(NPT),      ZS2(NPT)
C
         INCLUDE 'CUCCCC.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C     
         IF (KGUESS .EQ. 2) THEN
C
            BPS( 1) = R0O
            BPS(12) = RZ0O
C
            IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
            IF (NSURF .EQ. 6) CALL BNDSPL
C
            DO 6 J6=1,NS1
C
            DO 1 J1=1,NT
C
            ZR = CSIG(J6) * RHOS(J1) * COS(CT(J1)) + R0
            ZZ = CSIG(J6) * RHOS(J1) * SIN(CT(J1)) + RZ0
C
            ZS(J1) = SQRT((ZR - R0O)**2 + (ZZ - RZ0O)**2)
C     
            IF (ZR .EQ. R0O) THEN
               ZT(J1) = CTO(J1)
            ELSE
               ZT(J1) = ATAN2(ZZ - RZ0O,ZR - R0O)
            ENDIF
C 
            IF (ZT(J1) .LT. CTO(1)) ZT(J1) = ZT(J1) + 2. * CPI
C
    1       CONTINUE
C
            CALL BOUND(NT,ZT,ZBND)
C
            CALL RESETI(IC,NT,1)
            DO 2 JS = 1,NSO+1
               DO 2 JG=1,NT
                  IF (IC(JG).EQ.0) GOTO 2
                  ZSIGMA(JG) = ZS(JG) / ZBND(JG)
                  IS0(JG) = JS-1
                  IF (ZSIGMA(JG).LE.CSIGO(JS)) IC(JG)  = 0
 2          CONTINUE
            CALL RESETI(IC,NT,1)
            DO 3 JT = 1,NTO+1
               DO 3 JG=1,NT
                  IF (IC(JG).EQ.0) GOTO 3
                  IT0(JG) = JT-1
                  IF (ZT(JG).LE.CTO(JT)) IC(JG)  = 0
 3          CONTINUE
C     
         DO 4 J4=1,NT
C
         IF (IS0(J4).LT.  1) IS0(J4) = 1
         IF (IS0(J4).GE.NSO) IS0(J4) = NSO
         IF (IT0(J4).LT.  1) IT0(J4) = 1
         IF (IT0(J4).GE.NTO) IT0(J4) = NTO
C     
         ZS1(J4) = CSIGO(IS0(J4))
         ZS2(J4) = CSIGO(IS0(J4)+1)
         ZT1(J4) = CTO(IT0(J4))
         ZT2(J4) = CTO(IT0(J4)+1)
C     
         I = (IS0(J4) - 1) * NTO + IT0(J4)
C
         ZPCEL(J4, 1) = CPSIO(4*I-3)
         ZPCEL(J4, 2) = CPSIO(4*I-2)
         ZPCEL(J4, 3) = CPSIO(4*I-1)
         ZPCEL(J4, 4) = CPSIO(4*I  )
         ZPCEL(J4, 5) = CPSIO(4*(I+NTO)-3)
         ZPCEL(J4, 6) = CPSIO(4*(I+NTO)-2)
         ZPCEL(J4, 7) = CPSIO(4*(I+NTO)-1)
         ZPCEL(J4, 8) = CPSIO(4*(I+NTO)  )
C
         IF (IT0(J4) .NE. NTO) THEN
C
            ZPCEL(J4, 9) = CPSIO(4*I+1)
            ZPCEL(J4,10) = CPSIO(4*I+2)
            ZPCEL(J4,11) = CPSIO(4*I+3)
            ZPCEL(J4,12) = CPSIO(4*I+4)
            ZPCEL(J4,13) = CPSIO(4*(I+NTO)+1)
            ZPCEL(J4,14) = CPSIO(4*(I+NTO)+2)
            ZPCEL(J4,15) = CPSIO(4*(I+NTO)+3)
            ZPCEL(J4,16) = CPSIO(4*(I+NTO)+4)
C
         ELSE
C
            ZPCEL(J4, 9) = CPSIO(4*(I-NTO)+1)
            ZPCEL(J4,10) = CPSIO(4*(I-NTO)+2)
            ZPCEL(J4,11) = CPSIO(4*(I-NTO)+3)
            ZPCEL(J4,12) = CPSIO(4*(I-NTO)+4)
            ZPCEL(J4,13) = CPSIO(4*I+1)
            ZPCEL(J4,14) = CPSIO(4*I+2)
            ZPCEL(J4,15) = CPSIO(4*I+3)
            ZPCEL(J4,16) = CPSIO(4*I+4)
C
         ENDIF
C
 4       CONTINUE
C
         CALL BASIS1(NT,NPT,ZS1,ZS2,ZT1,ZT2,ZSIGMA,ZT,ZF)
C
         DO 5 J5=1,NT
C
         I = 4 * ((J6 - 1) * NT + J5) - 3
C
         CPSICL(I) = ZF(J5, 1) * ZPCEL(J5, 1) +
     +               ZF(J5, 2) * ZPCEL(J5, 2) +
     +               ZF(J5, 3) * ZPCEL(J5, 3) +
     +               ZF(J5, 4) * ZPCEL(J5, 4) +
     +               ZF(J5, 5) * ZPCEL(J5, 5) +
     +               ZF(J5, 6) * ZPCEL(J5, 6) +
     +               ZF(J5, 7) * ZPCEL(J5, 7) +
     +               ZF(J5, 8) * ZPCEL(J5, 8) +
     +               ZF(J5, 9) * ZPCEL(J5, 9) +
     +               ZF(J5,10) * ZPCEL(J5,10) +
     +               ZF(J5,11) * ZPCEL(J5,11) +
     +               ZF(J5,12) * ZPCEL(J5,12) +
     +               ZF(J5,13) * ZPCEL(J5,13) +
     +               ZF(J5,14) * ZPCEL(J5,14) +
     +               ZF(J5,15) * ZPCEL(J5,15) +
     +               ZF(J5,16) * ZPCEL(J5,16)
C
 5        CONTINUE
 6        CONTINUE
C
C SMOOTH THE NEW SOLUTION WITH BICUBIC SPLINES AND COMPUTE
C DERIVATIVES ON THE (SIGMA; THETA) GRID
C
            IF (NSMOOTH.EQ.1) CALL SMOOTH
C
            DO 7 J7=1,NT
C
            I = 4 * (J7 - 1)
C
            CPSICL(I+2) = 0.
            CPSICL(I+3) = 0.
            CPSICL(I+4) = 0.
C
    7       CONTINUE
C
            BPS( 1) = R0
            BPS(12) = RZ0
C
            IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
            IF (NSURF .EQ. 6) CALL BNDSPL
C
            SPSIM = CPSICL(1)
C
            IF (NSTTP .GE. 2) THEN
C
               CALL RESETI(IC,NISO,1)
               DO 8 JS = 1,NISOO
                  DO 8 JG=1,NISO
                     IF (IC(JG).EQ.0) GOTO 8
                     IS0(JG) = JS-1
                     IF (CSIPR(JG).LE.CSIPRO(JS)) IC(JG)  = 0
 8             CONTINUE
C
               DO 9 J9=1,NISO
C
               ISIPR = IS0(J9)
C
               IF (ISIPR .GE. NISOO - 2) ISIPR = NISOO - 2
               IF (ISIPR .LE. 2)         ISIPR = 2
C
               CID0(J9)   = FCCCC0(CID0O(ISIPR-1),CID0O(ISIPR),
     ,                             CID0O(ISIPR+1),CID0O(ISIPR+2),
     ,                             CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                             CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                             CSIPR(J9))
               D2CID0(J9) = FCCCC0(D2CID0O(ISIPR-1),D2CID0O(ISIPR),
     ,                             D2CID0O(ISIPR+1),D2CID0O(ISIPR+2),
     ,                             CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                             CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                             CSIPR(J9))
               CID2(J9)   = FCCCC0(CID2O(ISIPR-1),CID2O(ISIPR),
     ,                             CID2O(ISIPR+1),CID2O(ISIPR+2),
     ,                             CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                             CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                             CSIPR(J9))
               D2CID2(J9) = FCCCC0(D2CID2O(ISIPR-1),D2CID2O(ISIPR),
     ,                             D2CID2O(ISIPR+1),D2CID2O(ISIPR+2),
     ,                             CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                             CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                             CSIPR(J9))
               TTP(J9)    = FCCCC0(TTPO(ISIPR-1),TTPO(ISIPR),
     ,                             TTPO(ISIPR+1),TTPO(ISIPR+2),
     ,                             CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                             CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                             CSIPR(J9))
               CPPR(J9)   = FCCCC0(CPPRO(ISIPR-1),CPPRO(ISIPR),
     ,                             CPPRO(ISIPR+1),CPPRO(ISIPR+2),
     ,                             CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                             CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                             CSIPR(J9))
C
               IF (NSTTP .EQ. 3) THEN
C
                  TMF(J9)    = FCCCC0(TMFO(ISIPR-1),TMFO(ISIPR),
     ,                                TMFO(ISIPR+1),TMFO(ISIPR+2),
     ,                                CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                                CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                                CSIPR(J9))
                  D2TMF(J9)  = FCCCC0(D2TMFO(ISIPR-1),D2TMFO(ISIPR),
     ,                                D2TMFO(ISIPR+1),D2TMFO(ISIPR+2),
     ,                                CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                                CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                                CSIPR(J9))
C
               ENDIF
C
               IF (NPROFZ .EQ. 1) THEN
C
                  D2CPPR(J9) = FCCCC0(D2CPPRO(ISIPR-1),D2CPPRO(ISIPR),
     ,                                D2CPPRO(ISIPR+1),D2CPPRO(ISIPR+2),
     ,                                CSIPRO(ISIPR-1),CSIPRO(ISIPR),
     ,                                CSIPRO(ISIPR+1),CSIPRO(ISIPR+2),
     ,                                CSIPR(J9))
C
               ENDIF
C
    9          CONTINUE
C
            ENDIF
C
         ELSE IF (KGUESS .EQ. 1) THEN
C
            SPSIM = - 0.1
CYQL2018
            IF (NSTTP.EQ.4) SPSIM = -1.0e-1

            RMAG  = R0
            RZMAG = RZ0
C
            DO 11 J11=1,NS1
C
            DO 10 J10=1,NT
C
            I = (J11 - 1) * NT + J10
C
            CPSICL(4*I-3) = SPSIM * (1 - CSIG(J11)**2)
            CPSICL(4*I-2) = - 2. * SPSIM * CSIG(J11)
            CPSICL(4*I-1) = 0.
            CPSICL(4*I  ) = 0.
C
   10       CONTINUE
   11       CONTINUE
C
            IF (NSTTP .GE. 2) THEN
C
               CALL VZERO(CID0,2*NPISO)
               CALL VZERO(CIDQ,2*NPISO)
               CALL VZERO(CIDR,2*NPISO)
               CALL VZERO(CID2,2*NPISO)
               CALL VZERO(SIGPSI,2*NPMGS*NTP1*NPISO)
               CALL VZERO(TETPSI,2*NPMGS*NTP1*NPISO)
               CALL VZERO(WGTPSI,2*NPMGS*NTP1*NPISO)
               CALL RESETR(TMF,2*NPISO,RC1P)
C
               DO 12 J12=1,NISO
C
               PSIISO(J12)   = SPSIM * (1. - CSIPR(J12)**2)
               TETMAP(1,J12) = 0.
               SIGMAP(1,J12) = CSIPR(J12)
C
   12          CONTINUE
C
               CALL ISOFIND(1,NISO,SIGPSI,TETPSI,WGTPSI,SPSIM,RC0P)
C
               DO 14 J14=1,NISO
C
               CALL CINT(J14,SIGPSI(1,J14),TETPSI(1,J14),WGTPSI(1,J14))
C
   14          CONTINUE
C
               IF (NSTTP .EQ. 3) THEN
C
                  IF (NCSCAL .EQ. 1 .OR. NCSCAL .EQ. 3) THEN
C
                     CALL PRFUNC(1,SPSIM,ZJDB0,0)
C
                     SCALE  = (2. / (QSPEC * ZJDB0))**SCEXP
                     SCALAC = SCALE * SCALAC
C
                  ELSE IF (NCSCAL .EQ. 2) THEN
C
                     PRINT*,'NOREPT : NSTTP=3; NCSCAL=2; NOPT=',NOPT,';'
                     PRINT*,'NBLOPT=',NBLOPT,'; CPRESS=',CPRESS,';'
                     PRINT*,'THIS OPTION IS NOT POSSIBLE'
c                     STOP               
C
                  ENDIF
C
                  IF (NFUNC .EQ. 1) THEN
C
                     CALL SSCAL(NSOUR,SCALE,AT,1)
C
                  ELSE IF (NFUNC .EQ. 2) THEN
C
                     CALL SSCAL(5,SCALE,AT(3),1)
                     CALL SSCAL(5,SCALE,AT2(3),1)
                     CALL SSCAL(5,SCALE,AT3(3),1)
C
                     AT4(3) = SCALE * AT4(3)
C
                  ELSE IF (NFUNC .EQ. 3) THEN
C
                     AT(1) = SCALE * AT(1)
C
                  ELSE IF (NFUNC .EQ. 4) THEN
C
                     CALL SSCAL(NPPF+1,SCALE,RFUN,1)
C
                  ENDIF
C
                  CALL RVAR('SCALE            ',SCALE)
                  CALL RVAR('ACCUMULATED SCALE',SCALAC)
C
               ENDIF
C
               CALL ISOFUN(NISO)
C
            ENDIF
         ENDIF
C
         RETURN
         END
C*DECK C2SC01 
C*CALL PROCESS 
         SUBROUTINE MAGAXE
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SC01 COMPUTE POSITION OF MAGNETIC AXIS. FINDS THE MINIMUM OF PSI  *
*        BY CONJUGATE GRADIENT METHOD                                 *
*        (SEE NUMERICAL RECIPES P.303-306)                            *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION  ZG(2),  ZH(2), ZXI(2)
         parameter (ndivs = 50)
CLiu020922         static zmagold
         REAL zmagold
         zdiv = aspct/float(ndivs)
C
         IF (NMAG .EQ. 0) THEN
C
            ZR0  = R0 + .1 * ASPCT
            ZZ0  = RZ0
            NMAG = 1
C
         ELSE IF (NMAG .EQ. 1) THEN
C
            ZR0 = RMAG
            ZZ0 = RZMAG
C
         ENDIF
C
         IF (RMAG .EQ. R0) ZR0 = R0 + .05 * ASPCT
C
C INITIALIZATION FOR CONJUGATE GRADIENT METHOD
C
         CALL EVLATE(1,ZR0,ZZ0,ZXI(1),ZXI(2),ZPSI1)
         CALL EVLATE(2,ZR0,ZZ0,ZXI(1),ZXI(2),ZPSI0)
         if (zpsi1.lt.0.) goto 100
         write(*,*) ' psivalue ',zpsi1
         stop 'psi>0'
 100     continue
C
         ZG(1)  = - ZXI(1)
         ZG(2)  = - ZXI(2)
         ZH(1)  = ZG(1)
         ZH(2)  = ZG(2)
         ZXI(1) = ZH(1)
         ZXI(2) = ZH(2)
         ZROLD  = ZR0
         ZZOLD  = ZZ0
C
         DO 8 J8=1,50
C
C BRACKET OUT THE MINIMUM IN DIRECTION ZXI
C
         IF (ABS(ZXI(1)) .LT. RC1M14 .AND. ABS(ZXI(2)) .LT. RC1M14) THEN
C
            GOTO 10
C
         ENDIF
C
         ZSCL = SQRT(ZXI(1)**2 + ZXI(2)**2)
         ZR1  = ZR0 + zdiv * ZXI(1) / ZSCL 
         ZZ1  = ZZ0 + zdiv * ZXI(2) / ZSCL
C
         ZDPDXI0 = -ZSCL
C
         DO 2 J2=1,ndivs
C
         CALL EVLATE(2,ZR1,ZZ1,ZDPDR1,ZDPDZ1,ZPSI1)
C
         ZDPDXI1 = (ZDPDR1 * ZXI(1) + ZDPDZ1 * ZXI(2)) / ZSCL
C
         IF (ZDPDXI1 * ZDPDXI0 .GT. 0.) THEN
C
            ZR0     = ZR1
            ZZ0     = ZZ1
            ZDPDXI0 = ZDPDXI1
            ZR1     = ZR0 + zdiv * ZXI(1) / ZSCL 
            ZZ1     = ZZ0 + zdiv * ZXI(2) / ZSCL
C
         ELSE
C
            GOTO 3
C
         ENDIF
C
    2    CONTINUE
C
         PRINT*,'BRACKETING NOT SUCCEEDED'
         STOP
C
    3    CONTINUE
C
C COMPUTE THE MINIMUM IN DIRECTION ZXI BY BISSECTION 
C
         DO 4 J4=1,100
C
         ZRM = .5 * (ZR0 + ZR1)
         ZZM = .5 * (ZZ0 + ZZ1)
C
         CALL EVLATE(2,ZRM,ZZM,ZDPDRM,ZDPDZM,ZPSIM)
C
         ZDPDXIM = (ZDPDRM * ZXI(1) + ZDPDZM * ZXI(2)) / ZSCL
C
         IF ((ABS(ZR1 - ZR0) .LT. RC1M12 .AND.
     +        ABS(ZZ1 - ZZ0) .LT. RC1M12) .OR.
     +        ABS(ZDPDXIM) .LT. RC1M14) GOTO 5
C
         IF (ZDPDXI0 * ZDPDXIM .LT. 0.) THEN
C
            ZR1     = ZRM
            ZZ1     = ZZM
            ZDPDXI1 = ZDPDXIM
C
         ELSE
C
            ZR0     = ZRM
            ZZ0     = ZZM
            ZDPDXI0 = ZDPDXIM
C
         ENDIF

C
    4    CONTINUE
C
         PRINT*,'BISSECTION NOT SUCCEEDED'
         STOP
C
    5    CONTINUE
C
         ZR0 = ZRM
         ZZ0 = ZZM
C
         IF ((ABS(ZROLD - ZR0) .LT. RC1M12 .AND.
     +        ABS(ZZOLD - ZZ0) .LT. RC1M12) .OR.
     +       (ABS(ZDPDRM) .LT. RC1M14 .AND. 
     +        ABS(ZDPDZM) .LT. RC1M14)) GOTO 10
C
C INITIALIZATION FOR NEXT STEP OF CONJUGATE GRADIENT METHOD
C
         ZROLD = ZR0
         ZZOLD = ZZ0
C
         CALL EVLATE(2,ZR0,ZZ0,ZXI(1),ZXI(2),ZPSI0)
C
         ZGG  = ZG(1)**2 + ZG(2)**2
CC
CC FLETCHER-REEVES
CC
C         ZDGG = ZXI(1)**2 + ZXI(2)**2 
C
C POLAK-RIBIERE
C
         ZDGG = ZXI(1) * (ZG(1) + ZXI(1)) + ZXI(2) * (ZG(2) + ZXI(2))
C
         IF (ZGG .EQ. 0.) GOTO 10
C
         ZGAM = ZDGG / ZGG
C
         ZG(1)  = - ZXI(1)
         ZG(2)  = - ZXI(2)
         ZH(1)  = ZG(1) + ZGAM * ZH(1)
         ZH(2)  = ZG(2) + ZGAM * ZH(2)
         ZXI(1) = ZH(1)
         ZXI(2) = ZH(2)
C
    8    CONTINUE
C
         CALL RVAR('ZROLD - ZR0',ZROLD - ZR0)
         CALL RVAR('ZZOLD - ZZ0',ZZOLD - ZZ0)
         CALL RVAR('ZDPDRM',ZDPDRM)
         CALL RVAR('ZDPDZM',ZDPDZM)
         PRINT*,'CONJUGATE GRADIENT NOT CONVERGED'
         write(*,*) ' warning magnetic axis not found'
         goto 20
c        STOP
C
   10    CONTINUE
C
         CALL EVLATE(1,ZR0,ZZ0,ZDPDR0,ZDPDZ0,SPSIM)
C
         RMAG  = ZR0
         RZMAG = ZZ0  
C
CLIU   write(*,*) ' nsym=',nsym
CLIU   write(*,*) ' rzmag =',rzmag
       if (nsym.ne.0) goto 15
       if (abs(rzmag-zmagold) .gt. ezmag) goto 20
       goto 30
cab
cab    sometimes the magnetic axis moves off the midplane after
cab    interpolation from small to large grid.  Then, under-relaxation
cab    usually helps it to go back to z = 0
cab
 15    if (abs(rzmag).lt.ezmag) goto 30
       rzmag = 0.
 20    relax = relax + 0.01
       if (relax.gt.0.80) relax = 0.80
C
 30    zmagold = rzmag
         RETURN
         END
C*DECK C2SC02
C*CALL PROCESS 
         SUBROUTINE EVLATE(KCASE,PR,PZ,PDPDR,PDPDZ,PSI)
C        ##############################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SC02 EVALUATE PSI, D(PSI)/D(R) AND D(PSI)/D(Z) AT (R,Z)           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R      ZBND(5),    ZCPSI(1,16),     ZDFDS(1,16),    ZDFDT(1,16),   
     R      ZF(1,16),   ZTET(5)
C
         ZEPS = 1.E-3
C
         ZRHO    = SQRT((PR - R0)**2 + (PZ - RZ0)**2)
         ZTET(1) = ATAN2(PZ - RZ0,PR - R0)
C
         IF (ZTET(1) .LT. CT(1)) ZTET(1) = ZTET(1) + 2. * CPI
C
         IF (KCASE .EQ. 1) THEN
C
            CALL BOUND(1,ZTET,ZBND)
C
         ELSE IF (KCASE .EQ. 2) THEN
C
            ZTET(2) = ZTET(1) - 2. * ZEPS
            ZTET(3) = ZTET(1) -      ZEPS
            ZTET(4) = ZTET(1) +      ZEPS
            ZTET(5) = ZTET(1) + 2. * ZEPS
C
            CALL BOUND(5,ZTET,ZBND)
C
        ENDIF
C
         ZSIG = ZRHO / ZBND(1)
C
         IS = ISRCHFGE(NS1,CSIG,1,ZSIG)  - 1
         IT = ISRCHFGE(NT1,CT,1,ZTET(1)) - 1
C
         IF (IS .LT. 1)  IS = 1
         IF (IS .GT. NS) IS = NS
         IF (IT .LT. 1)  IT = 1
         IF (IT .GT. NT) IT = NT
C
         ZS1 = CSIG(IS)
         ZS2 = CSIG(IS+1)
         ZT1 = CT(IT)
         ZT2 = CT(IT+1)
C
         CALL PSICEL(IS,IT,1,1,ZCPSI,CPSICL)
C
         IF (KCASE .EQ. 1) THEN
C
            CALL BASIS1(1,1,ZS1,ZS2,ZT1,ZT2,ZSIG,ZTET(1),ZF)
C
            PSI = SDOT(16,ZF,1,ZCPSI,1)
C
         ELSE IF (KCASE .EQ. 2) THEN
C
C EVALUATE FIRST DERIVATIVES WITH RESPECT TO SIGMA AND THETA
C
            CALL BASIS2(1,1,ZS1,ZS2,ZT1,ZT2,ZSIG,ZTET(1),ZDFDS,ZDFDT)
C
            ZDPDS = SDOT(16,ZDFDS,1,ZCPSI,1)
            ZDPDT = SDOT(16,ZDFDT,1,ZCPSI,1)
C
            ZDRSDT = (ZBND(2) + 8*(ZBND(4) - ZBND(3)) - ZBND(5)) / 
     /               (12. * ZEPS)
C
            ZCOST = COS(ZTET(1))
            ZSINT = SIN(ZTET(1))
C
            ZDSDR = (ZDRSDT * ZSINT + ZBND(1) * ZCOST) / ZBND(1)**2
            ZDTDR = - ZSINT / ZRHO
            ZDSDZ = (ZBND(1) * ZSINT - ZDRSDT * ZCOST) / ZBND(1)**2
            ZDTDZ = ZCOST / ZRHO
C
            PDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
            PDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         ENDIF
C
         RETURN
         END
C*DECK C2SD01
C*CALL PROCESS
         SUBROUTINE SETUPA
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SD01 CONSTRUCT MATRIX A. THIS MATRIX IS OBTAINED FROM THE LEFT    *
*        HAND SIDE OF EQ. (27) IN THE PUBLICATION                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R     ZDBDS(NPT,16), ZDBDT(NPT,16), ZS(NPT),       ZS1(NPT),     
     R     ZS2(NPT),      ZT(NPT),       ZT1(NPT),      ZT2(NPT),     
     R     ZV(NPT,16),    ZXA(NPT,16,16)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         CALL VZERO(A,NPBAND*NP4NST)
C
         DO 1 J1=1,NT
C
         ZT1(J1) = CT(J1)
         ZT2(J1) = CT(J1+1)
C
    1    CONTINUE
C
***********************************************************************
*                                                                     *
* 1. SCAN OVER ALL INTERVALS                                          *
*                                                                     *
***********************************************************************
C
         DO 22 J22=1,NS
C
***********************************************************************
*                                                                     *
* 1.1. INITIALIZATION OF LOCAL ARRAYS                                 *
*                                                                     *
***********************************************************************
C
         CALL VZERO(ZXA,256*NPT)
C
***********************************************************************
*                                                                     *
* 1.2 COMPUTE VERTICAL POSITIONS IN A                                 *
*                                                                     *
***********************************************************************
C
         DO 2 J2=1,NT-1
C
         I = (J22 - 1) * NT + J2
C
         NPLAC(J2,1) = NUPDWN(I)
         NPLAC(J2,2) = NUPDWN(I+NT)
         NPLAC(J2,3) = NUPDWN(I+1)
         NPLAC(J2,4) = NUPDWN(I+NT+1)
C
    2    CONTINUE
C
         I = J22 * NT
C
         NPLAC(NT,1) = NUPDWN(I)
         NPLAC(NT,2) = NUPDWN(I+NT)
         NPLAC(NT,3) = NUPDWN(I+1-NT)
         NPLAC(NT,4) = NUPDWN(I+1)
C
         DO 3 J3=1,NT
C
         ZS1(J3) = CSIG(J22)
         ZS2(J3) = CSIG(J22+1)
C
    3    CONTINUE
C
***********************************************************************
*                                                                     *
* 1.3. COMPUTE INTEGRANTS OF VARITIONNAL PROBLEM                      *
*                                                                     *
***********************************************************************
C
         DO 12 J12=1,NWGAUS
C
         DO 4 J4=1,NT
C
         ZS(J4) = RSINT(J22,J12)
         ZT(J4) = RTINT(J4,J12)
C
  4      CONTINUE
C
         CALL BASIS2(NT,NPT,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT)
C
         DO 8 J8=1,NT
C
         ZFRAC = YDRSDT(J8,J12) / YRST(J8,J12)
C
         ZV(J8, 1) = ZDBDT(J8, 1)/RSINT(J22,J12) - ZDBDS(J8, 1)*ZFRAC
         ZV(J8, 2) = ZDBDT(J8, 2)/RSINT(J22,J12) - ZDBDS(J8, 2)*ZFRAC
         ZV(J8, 3) = ZDBDT(J8, 3)/RSINT(J22,J12) - ZDBDS(J8, 3)*ZFRAC
         ZV(J8, 4) = ZDBDT(J8, 4)/RSINT(J22,J12) - ZDBDS(J8, 4)*ZFRAC
         ZV(J8, 5) = ZDBDT(J8, 5)/RSINT(J22,J12) - ZDBDS(J8, 5)*ZFRAC
         ZV(J8, 6) = ZDBDT(J8, 6)/RSINT(J22,J12) - ZDBDS(J8, 6)*ZFRAC
         ZV(J8, 7) = ZDBDT(J8, 7)/RSINT(J22,J12) - ZDBDS(J8, 7)*ZFRAC
         ZV(J8, 8) = ZDBDT(J8, 8)/RSINT(J22,J12) - ZDBDS(J8, 8)*ZFRAC
         ZV(J8, 9) = ZDBDT(J8, 9)/RSINT(J22,J12) - ZDBDS(J8, 9)*ZFRAC
         ZV(J8,10) = ZDBDT(J8,10)/RSINT(J22,J12) - ZDBDS(J8,10)*ZFRAC
         ZV(J8,11) = ZDBDT(J8,11)/RSINT(J22,J12) - ZDBDS(J8,11)*ZFRAC
         ZV(J8,12) = ZDBDT(J8,12)/RSINT(J22,J12) - ZDBDS(J8,12)*ZFRAC
         ZV(J8,13) = ZDBDT(J8,13)/RSINT(J22,J12) - ZDBDS(J8,13)*ZFRAC
         ZV(J8,14) = ZDBDT(J8,14)/RSINT(J22,J12) - ZDBDS(J8,14)*ZFRAC
         ZV(J8,15) = ZDBDT(J8,15)/RSINT(J22,J12) - ZDBDS(J8,15)*ZFRAC
         ZV(J8,16) = ZDBDT(J8,16)/RSINT(J22,J12) - ZDBDS(J8,16)*ZFRAC
C
    8    CONTINUE
C
***********************************************************************
*                                                                     *
* 1.4. PERFORMS DIADIC PRODUCT AND COMPUTE MATRIX CONTRIBUTION        *
*                                                                     *
***********************************************************************
C
         DO 11 J11=1,16
C
         DO 9 J9=1,NT
C
         ZW    = CW(J12) * (CSIG(J22+1) - CSIG(J22)) *
     *                  (CT(J9+1) - CT(J9))
         ZR    = RSINT(J22,J12)*YRST(J9,J12)*COS(RTINT(J9,J12)) + R0
         ZCOEF = ZW * RSINT(J22,J12) / ZR
C
         ZXA(J9, 1,J11) = ZXA(J9, 1,J11) + ZCOEF *
     *      (ZV(J9, 1)*ZV(J9,J11) + ZDBDS(J9, 1)*ZDBDS(J9,J11)) 
         ZXA(J9, 2,J11) = ZXA(J9, 2,J11) + ZCOEF *
     *      (ZV(J9, 2)*ZV(J9,J11) + ZDBDS(J9, 2)*ZDBDS(J9,J11)) 
         ZXA(J9, 3,J11) = ZXA(J9, 3,J11) + ZCOEF *
     *      (ZV(J9, 3)*ZV(J9,J11) + ZDBDS(J9, 3)*ZDBDS(J9,J11)) 
         ZXA(J9, 4,J11) = ZXA(J9, 4,J11) + ZCOEF *
     *      (ZV(J9, 4)*ZV(J9,J11) + ZDBDS(J9, 4)*ZDBDS(J9,J11)) 
         ZXA(J9, 5,J11) = ZXA(J9, 5,J11) + ZCOEF *
     *      (ZV(J9, 5)*ZV(J9,J11) + ZDBDS(J9, 5)*ZDBDS(J9,J11)) 
         ZXA(J9, 6,J11) = ZXA(J9, 6,J11) + ZCOEF *
     *      (ZV(J9, 6)*ZV(J9,J11) + ZDBDS(J9, 6)*ZDBDS(J9,J11)) 
         ZXA(J9, 7,J11) = ZXA(J9, 7,J11) + ZCOEF *
     *      (ZV(J9, 7)*ZV(J9,J11) + ZDBDS(J9, 7)*ZDBDS(J9,J11)) 
         ZXA(J9, 8,J11) = ZXA(J9, 8,J11) + ZCOEF *
     *      (ZV(J9, 8)*ZV(J9,J11) + ZDBDS(J9, 8)*ZDBDS(J9,J11)) 
         ZXA(J9, 9,J11) = ZXA(J9, 9,J11) + ZCOEF *
     *      (ZV(J9, 9)*ZV(J9,J11) + ZDBDS(J9, 9)*ZDBDS(J9,J11)) 
         ZXA(J9,10,J11) = ZXA(J9,10,J11) + ZCOEF *
     *      (ZV(J9,10)*ZV(J9,J11) + ZDBDS(J9,10)*ZDBDS(J9,J11)) 
         ZXA(J9,11,J11) = ZXA(J9,11,J11) + ZCOEF *
     *      (ZV(J9,11)*ZV(J9,J11) + ZDBDS(J9,11)*ZDBDS(J9,J11)) 
         ZXA(J9,12,J11) = ZXA(J9,12,J11) + ZCOEF *
     *      (ZV(J9,12)*ZV(J9,J11) + ZDBDS(J9,12)*ZDBDS(J9,J11)) 
         ZXA(J9,13,J11) = ZXA(J9,13,J11) + ZCOEF *
     *      (ZV(J9,13)*ZV(J9,J11) + ZDBDS(J9,13)*ZDBDS(J9,J11)) 
         ZXA(J9,14,J11) = ZXA(J9,14,J11) + ZCOEF *
     *      (ZV(J9,14)*ZV(J9,J11) + ZDBDS(J9,14)*ZDBDS(J9,J11)) 
         ZXA(J9,15,J11) = ZXA(J9,15,J11) + ZCOEF *
     *      (ZV(J9,15)*ZV(J9,J11) + ZDBDS(J9,15)*ZDBDS(J9,J11)) 
         ZXA(J9,16,J11) = ZXA(J9,16,J11) + ZCOEF *
     *      (ZV(J9,16)*ZV(J9,J11) + ZDBDS(J9,16)*ZDBDS(J9,J11)) 
C
    9    CONTINUE
   11    CONTINUE
   12    CONTINUE
C
***********************************************************************
*                                                                     *
* 1.5. ADD TO MATRIX A                                                *
*                                                                     *
***********************************************************************
*                                                                     *
* 1.5.1. ADD DIAGONAL BLOCS OF ZXA TO A :                             *
*                                                                     *
*                   XA    =====>                       A              *
*                   ==                                 =              *
*                                                                     *
* 4*(J16-1)+1 --> * * * *     4*(NPLAC(J16)-1)+1 --> * * * *          *
*                 ! * * *                            * * *            *
*                 !   * * =====>                     * *              *
*                 !     *                            *                *
*                 !                                                   *
*                -!--------> IXACOL                 -!-----!--> IACOL *
*                 4*(J16-1)+1                        1     4          *
*                                                                     *
***********************************************************************
C
***********************************************************************
*                                                                     *
*    INDEXATION:                                                      *
*    -----------                                                      *
*      J16 ===> VERTICAL POSITION OF BLOC                             *
*      J14 ===> ROW IN BLOC                                           *
*      J13 ===> COLUMN IN BLOC                                        *
*                                                                     *
***********************************************************************
C
C
         DO 16 J16=1,4
C
         DO 15 J15=1,NT
C
         DO 14 J14=1,4
C
         IMAX = 5 - J14
C
         DO 13 J13=1,IMAX
C
         IAROW  = 4 * (NPLAC(J15,J16) - 1) + J14
         IACOL  = J13
         IXAROW = 4 * (J16 - 1) + J14
         IXACOL = IXAROW + J13 - 1
C
         A(IACOL,IAROW) = A(IACOL,IAROW) + ZXA(J15,IXACOL,IXAROW)
C
   13    CONTINUE
   14    CONTINUE
   15    CONTINUE
   16    CONTINUE
C
***********************************************************************
*                                                                     *
*  1.5.2. ADD OUT OF DIAGONAL BLOCS OF ZXA TO A :                     *
*                                                                     *
*          ZXA            =====>              A                       *
*          ===                                =                       *
*                                                                     *
* 4*(J21-1)+1 --> * * * *      4*(NPLAC(J21)-1)+1 --> * * * *         *
*                 * * * *                           * * * *           *
*                 * * * * =====>                  * * * *             *
*                 * * * *                       * * * *               *
*                 !                                   !               *
*               --!-------> IXACOL             -------!--------> IACOL*
*                 4*(J21-1)+1                  4*(KPLAC(J20,J21)-1)+1 *
*                                                                     *
***********************************************************************
C
***********************************************************************
*                                                                     *
*    INDEXATION:                                                      *
*    -----------                                                      *
*      J21 ===> VERTICAL POSITION OF BLOC                             *
*      J20 ===> HORIZONTAL POSITION OF BLOC                           *
*      J18 ===> ROW IN BLOC                                           *
*      J17 ===> COLUMN IN BLOC                                        *
*                                                                     *
***********************************************************************
C
         DO 21 J21=1,4
C
         DO 20 J20=1,4
C
         DO 19 J19=1,NT
C
***********************************************************************
*                                                                     *
* TESTS IF BLOC MUST BE ADDED TO A                                    *
*                                                                     *
***********************************************************************
C
         IF (MPLA(J19,J20,J21) .GT. 1) THEN
C
            DO 18 J18=1,4
C
            DO 17 J17=1,4
C
            IAROW  = 4 * (NPLAC(J19,J21) - 1) + J18
            IACOL  = 4 * (MPLA(J19,J20,J21)-1) + J17 - J18 + 1
            IXAROW = 4 * (J21 - 1) + J18
            IXACOL = 4 * (J20 - 1) + J17
C
            A(IACOL,IAROW) = A(IACOL,IAROW) + ZXA(J19,IXACOL,IXAROW)
C
   17       CONTINUE
   18       CONTINUE
C
         ENDIF
C
   19    CONTINUE
   20    CONTINUE
   21    CONTINUE
   22    CONTINUE
C
***********************************************************************
*                                                                     *
* 2. INTRODUCE LIMIT CONDITIONS                                       *
*                                                                     *
***********************************************************************
C
         CALL LIMITA
C
         RETURN
         END
C*DECK C2SD02
C*CALL PROCESS
         SUBROUTINE SETUPB
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SD02 CONSTRUCT VECTOR B. THIS VECTOR IS OBTAINED FROM THE RIGHT   *
*        HAND SIDE OF EQ. (27) IN THE PUBLICATION                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R   ZJIPHI(NPT),  ZPCEL(NPT,16),    ZPSI(NPT,NPGAUS),   ZR(NPT),
     R   ZXB(NPT,16),  ZCUR1(NPT)
C
C----*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
***********************************************************************
*                                                                     *
* 1. INITIALIZATION                                                   *
*                                                                     *
***********************************************************************
C
         CALL VZERO(B,N4NSNT)
C
         ZCUR = 0.
C
***********************************************************************
*                                                                     *
* 2. SCAN OVER ALL CELLS                                              *
*                                                                     *
***********************************************************************
C
         DO 14 J14=1,NS
C
***********************************************************************
*                                                                     *
* 2.1. INITIALIZATION OF LOCAL ARRAYS                                 *
*                                                                     *
***********************************************************************
C
         DO 2 J2=1,16
C
         DO 1 J1=1,NT
C
         ZXB(J1,J2) = 0.
C
    1    CONTINUE
    2    CONTINUE
C
***********************************************************************
*                                                                     *
* 2.2. COMPUTE ALL QUANTITIES TO DEFINE PSI ON CELL                   *
*                                                                     *
***********************************************************************
C
         DO 4 J4=1,NT
C
         I = (J14 - 1) * NT + J4
C
         ZPCEL(J4,1) = CPSICL(4*I-3)
         ZPCEL(J4,2) = CPSICL(4*I-2)
         ZPCEL(J4,3) = CPSICL(4*I-1)
         ZPCEL(J4,4) = CPSICL(4*I  )
         ZPCEL(J4,5) = CPSICL(4*(I+NT)-3)
         ZPCEL(J4,6) = CPSICL(4*(I+NT)-2)
         ZPCEL(J4,7) = CPSICL(4*(I+NT)-1)
         ZPCEL(J4,8) = CPSICL(4*(I+NT)  )
C
         IF (J4 .EQ. NT) THEN
C
            ZPCEL(J4, 9) = CPSICL(4*(I-NT)+1)
            ZPCEL(J4,10) = CPSICL(4*(I-NT)+2)
            ZPCEL(J4,11) = CPSICL(4*(I-NT)+3)
            ZPCEL(J4,12) = CPSICL(4*(I-NT)+4)
            ZPCEL(J4,13) = CPSICL(4*I+1)
            ZPCEL(J4,14) = CPSICL(4*I+2)
            ZPCEL(J4,15) = CPSICL(4*I+3)
            ZPCEL(J4,16) = CPSICL(4*I+4)
C
         ELSE
C
            ZPCEL(J4, 9) = CPSICL(4*I+1)
            ZPCEL(J4,10) = CPSICL(4*I+2)
            ZPCEL(J4,11) = CPSICL(4*I+3)
            ZPCEL(J4,12) = CPSICL(4*I+4)
            ZPCEL(J4,13) = CPSICL(4*(I+NT)+1)
            ZPCEL(J4,14) = CPSICL(4*(I+NT)+2)
            ZPCEL(J4,15) = CPSICL(4*(I+NT)+3)
            ZPCEL(J4,16) = CPSICL(4*(I+NT)+4)
C
         ENDIF
C
    4    CONTINUE
C
***********************************************************************
*                                                                     *
* 2.3 COMPUTE VERTICAL POSITIONS IN A                                 *
*                                                                     *
***********************************************************************
C
***********************************************************************
*                                                                     *
* 2.4 COMPUTE PSI, R AND SOURCE TERM ON INTEGRATION POINTS            *
*                                                                     *
***********************************************************************
C
         DO 10 J10=1,NWGAUS
C
         DO 5 J5=1,NT
C
         ZPSI(J5,J10) = ZPCEL(J5, 1) * FB(J5,J10, 1,J14) +
     *                  ZPCEL(J5, 2) * FB(J5,J10, 2,J14) +
     *                  ZPCEL(J5, 3) * FB(J5,J10, 3,J14) +
     *                  ZPCEL(J5, 4) * FB(J5,J10, 4,J14) +
     *                  ZPCEL(J5, 5) * FB(J5,J10, 5,J14) +
     *                  ZPCEL(J5, 6) * FB(J5,J10, 6,J14) +
     *                  ZPCEL(J5, 7) * FB(J5,J10, 7,J14) +
     *                  ZPCEL(J5, 8) * FB(J5,J10, 8,J14) +
     *                  ZPCEL(J5, 9) * FB(J5,J10, 9,J14) +
     *                  ZPCEL(J5,10) * FB(J5,J10,10,J14) +
     *                  ZPCEL(J5,11) * FB(J5,J10,11,J14) +
     *                  ZPCEL(J5,12) * FB(J5,J10,12,J14) +
     *                  ZPCEL(J5,13) * FB(J5,J10,13,J14) +
     *                  ZPCEL(J5,14) * FB(J5,J10,14,J14) +
     *                  ZPCEL(J5,15) * FB(J5,J10,15,J14) +
     *                  ZPCEL(J5,16) * FB(J5,J10,16,J14) 
C
         ZR(J5) = RSINT(J14,J10)*YRST(J5,J10)*COS(RTINT(J5,J10)) + R0
C
    5    CONTINUE   
C
         CALL CURENT(NT,ZPSI(1,J10),ZR,ZJIPHI)
C
         DO 7 J7=1,NT
C
***********************************************************************
*                                                                     *
* 2.5. INTEGRATION OF CURRENT                                         *
*                                                                     *
***********************************************************************
C
         ZW     = CW(J10) * (CSIG(J14+1)-CSIG(J14)) * (CT(J7+1)-CT(J7))
         ZINTJ  = - ZJIPHI(J7) * YRST(J7,J10)**2 * RSINT(J14,J10)
C
         ZCUR1(J7) = - ZW * ZINTJ
C
***********************************************************************
*                                                                     *
* 2.7. COMPUTE VECTOR CONTRIBUTIONS                                   *
*                                                                     *
***********************************************************************
C
         ZXB(J7, 1) = ZXB(J7, 1) - ZCUR1(J7) * FB(J7,J10, 1,J14)
         ZXB(J7, 2) = ZXB(J7, 2) - ZCUR1(J7) * FB(J7,J10, 2,J14)
         ZXB(J7, 3) = ZXB(J7, 3) - ZCUR1(J7) * FB(J7,J10, 3,J14)
         ZXB(J7, 4) = ZXB(J7, 4) - ZCUR1(J7) * FB(J7,J10, 4,J14)
         ZXB(J7, 5) = ZXB(J7, 5) - ZCUR1(J7) * FB(J7,J10, 5,J14)
         ZXB(J7, 6) = ZXB(J7, 6) - ZCUR1(J7) * FB(J7,J10, 6,J14)
         ZXB(J7, 7) = ZXB(J7, 7) - ZCUR1(J7) * FB(J7,J10, 7,J14)
         ZXB(J7, 8) = ZXB(J7, 8) - ZCUR1(J7) * FB(J7,J10, 8,J14)
         ZXB(J7, 9) = ZXB(J7, 9) - ZCUR1(J7) * FB(J7,J10, 9,J14)
         ZXB(J7,10) = ZXB(J7,10) - ZCUR1(J7) * FB(J7,J10,10,J14)
         ZXB(J7,11) = ZXB(J7,11) - ZCUR1(J7) * FB(J7,J10,11,J14)
         ZXB(J7,12) = ZXB(J7,12) - ZCUR1(J7) * FB(J7,J10,12,J14)
         ZXB(J7,13) = ZXB(J7,13) - ZCUR1(J7) * FB(J7,J10,13,J14)
         ZXB(J7,14) = ZXB(J7,14) - ZCUR1(J7) * FB(J7,J10,14,J14)
         ZXB(J7,15) = ZXB(J7,15) - ZCUR1(J7) * FB(J7,J10,15,J14)
         ZXB(J7,16) = ZXB(J7,16) - ZCUR1(J7) * FB(J7,J10,16,J14)
C
    7    CONTINUE
         ZCUR = ZCUR + SSUM(NT,ZCUR1,1)
   10    CONTINUE
C
***********************************************************************
*                                                                     *
* 2.8. ADD TO VECTOR B                                                *
*                                                                     *
***********************************************************************
C
CDIR$ IVDEP
         DO 11 J11=1,NT
C
         I = (J14 - 1) * NT + J11
C
         I1 = 4 * (NUPDWN(I)      - 1)
         I2 = 4 * (NUPDWN(I + NT) - 1)
C
         IF (J11 .EQ. NT) THEN
C
            I3 = 4 * (NUPDWN(I-NT+1) - 1)
            I4 = 4 * (NUPDWN(I+1)    - 1)
C
         ELSE
C
            I3 = 4 * (NUPDWN(I+1)    - 1)
            I4 = 4 * (NUPDWN(I+NT+1) - 1)
C
         ENDIF
C
         IBROW1  = I1 + 1
         IBROW2  = I1 + 2
         IBROW3  = I1 + 3
         IBROW4  = I1 + 4
         IBROW5  = I2 + 1
         IBROW6  = I2 + 2
         IBROW7  = I2 + 3
         IBROW8  = I2 + 4
         IBROW9  = I3 + 1
         IBROW10 = I3 + 2
         IBROW11 = I3 + 3
         IBROW12 = I3 + 4
         IBROW13 = I4 + 1
         IBROW14 = I4 + 2
         IBROW15 = I4 + 3
         IBROW16 = I4 + 4
C
         B(IBROW1 ) = B(IBROW1 ) + ZXB(J11, 1)
         B(IBROW2 ) = B(IBROW2 ) + ZXB(J11, 2)
         B(IBROW3 ) = B(IBROW3 ) + ZXB(J11, 3)
         B(IBROW4 ) = B(IBROW4 ) + ZXB(J11, 4)
         B(IBROW5 ) = B(IBROW5 ) + ZXB(J11, 5)
         B(IBROW6 ) = B(IBROW6 ) + ZXB(J11, 6)
         B(IBROW7 ) = B(IBROW7 ) + ZXB(J11, 7)
         B(IBROW8 ) = B(IBROW8 ) + ZXB(J11, 8)
         B(IBROW9 ) = B(IBROW9 ) + ZXB(J11, 9)
         B(IBROW10) = B(IBROW10) + ZXB(J11,10)
         B(IBROW11) = B(IBROW11) + ZXB(J11,11)
         B(IBROW12) = B(IBROW12) + ZXB(J11,12)
         B(IBROW13) = B(IBROW13) + ZXB(J11,13)
         B(IBROW14) = B(IBROW14) + ZXB(J11,14)
         B(IBROW15) = B(IBROW15) + ZXB(J11,15)
         B(IBROW16) = B(IBROW16) + ZXB(J11,16)
C
   11    CONTINUE
   14    CONTINUE
C
***********************************************************************
*                                                                     *
* 3. INTRODUCE LIMIT CONDITIONS                                       *
*                                                                     *
***********************************************************************
C
         CALL LIMITB
C
***********************************************************************
*                                                                     *
* 4. TOTAL TOROIDAL CURRENT                                           *
*                                                                     *
***********************************************************************
C
         CUROLD = ZCUR
C
         RETURN
         END
C*DECK C2SD03
C*CALL PROCESS
         SUBROUTINE LIMITA
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SD03 IMPOSE BOUNDARY CONDITIONS IN A. (SEE EQ. (30) IN PAPER)     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     R   ZC(N4NT,3)
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         CALL CENTER(ZC)
         CALL IDENTA(3,ZC)
C
         DO 1 J1=1,4*NT-3
C
         IDIAG  = J1
         IBAND1 = 1
         IBAND2 = IDIAG + NBAND - 1
C
         CALL AWAY(IDIAG,IBAND1,IBAND2)
C
    1    CONTINUE
C
***********************************************************************
*                                                                     *
* 5. PLASMA OUTSIDE BOUNDARY CONDITION : PSI = ZERO                   *
*                                                                     *
***********************************************************************
C
         DO 2 J2=1,NT
C
         IDIAG  = 4 * (NS * NT + J2 - 1) + 1
         IBAND1 = IDIAG - NBAND + 1
         IBAND2 = IDIAG + NBAND - 1
C
         CALL AWAY(IDIAG,IBAND1,IBAND2)
C
         IDIAG  = 4 * (NS * NT + J2 - 1) + 3
         IBAND1 = IDIAG - NBAND + 1
         IBAND2 = IDIAG + NBAND - 1
C
         CALL AWAY(IDIAG,IBAND1,IBAND2)
C
    2    CONTINUE
C
         RETURN
         END
C*DECK C2SD04
C*CALL PROCESS
         SUBROUTINE LIMITB
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SD04 IMPOSE BOUNDARY CONDITIONS IN B. (SEE EQ. (30) IN PAPER)     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     R   ZC(N4NT,3)
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
***********************************************************************
*                                                                     *
* 1. PSI IDENTIFICATION ON GEOMETRIC CENTER OF MESH :                 *
*                                                                     *
*      ADD COEFFICIENTS OF EQUATIONS NO. 4*(J-1)+1 ,J=1,NT-1 TO       *
*      COEFFICENTS OF EQUATION NO. 4*(NT-1)+1                         *
*                                                                     *
***********************************************************************
C
         CALL CENTER(ZC)
         CALL IDENTB(3,ZC)
C
         DO 1 J1=1,4*NT-3
C
         B(J1) = 0
C
    1    CONTINUE
C
         DO 2 J2=1,NT
C
         B(4*(NS*NT+J2-1)+1) = 0
         B(4*(NS*NT+J2-1)+3) = 0
C
    2    CONTINUE
C
         RETURN
         END
C*DECK C2SD05
C*CALL PROCESS
         SUBROUTINE IDENTA(KVAR,PC)
C        ##########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SD05 PERFORM ROW AND COLUMN OPERATIONS REQUIRED TO IMPOSE         *
*        BOUNDARY CONDITIONS IN A                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     R   PA(N4NT,3),  PAA(N4NT,N4NT),  PC(N4NT,3)
C
         INCLUDE 'BNDIND.inc'
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         CALL VZERO(PAA,N16NT*NPT)
         CALL VZERO(PA,12*NPT)
C
         DO 2 J2=1,4*NT
C
         IDIAG  = J2
         IBAND1 = 4 * NT + 1
         IBAND2 = MIN(IDIAG+NBAND-1, 8*NT)
C
         DO 1 J1=IBAND1,IBAND2
C
         JACOL  = INDCOL(J1,IDIAG)
         JAROW  = INDROW(J1,IDIAG)
         JAACOL = J1 - 4 * NT
C
         PAA(JAACOL,J2) = A(JACOL,JAROW)
C
    1    CONTINUE
    2    CONTINUE
C
         DO 4 J4=1,4*NT
C
         DO 3 J3=1,KVAR
C
         CALL SAXPY(4*NT,PC(J4,J3),PAA(1,J4),1,PA(1,J3),1)
C
    3    CONTINUE
    4    CONTINUE
C
         DO 6 J6=1,KVAR
C
         DO 5 J5=1,4*NT
C
         JCOL = 4 * NT + J5
         JROW = 4 * NT - KVAR + J6
C
         JACOL  = INDCOL(JCOL,JROW)
         JAROW  = INDROW(JCOL,JROW)
C
         A(JACOL,JAROW) = PA(J5,J6)
C
    5    CONTINUE
    6    CONTINUE
C
         CALL VZERO(PAA,N16NT*NPT)
         CALL VZERO(PA,12*NPT)
C
         DO 8 J8=1,4*NT
C
         DO 7 J7=1,4*NT
C
         JCOL = J7
         JROW = J8
C
         JACOL = INDCOL(JCOL,JROW)
         JAROW = INDROW(JCOL,JROW)
C
         PAA(J7,J8) = A(JACOL,JAROW)
C
    7    CONTINUE
    8    CONTINUE
C
         DO 10 J10=1,4*NT
C
         DO 9 J9=1,KVAR
C
         CALL SAXPY(4*NT,PC(J10,J9),PAA(1,J10),1,PA(1,J9),1)
C
    9    CONTINUE
   10    CONTINUE
C
         DO 12 J12=1,4*NT
C
         DO 11 J11=1,KVAR
C
         JPAAR = 4 * NT - KVAR + J11
C
         PAA(J12,JPAAR) = PA(J12,J11)
C
   11    CONTINUE
   12    CONTINUE
C
         CALL VZERO(PA,12*NPT)
C
         JR = 4 * NT - KVAR + 1
C
         DO 14 J14=1,4*NT
C
         DO 13 J13=1,KVAR
C
         CALL SAXPY(KVAR,PC(J14,J13),PAA(J14,JR),N4NT,PA(J13,1),N4NT)
C
   13    CONTINUE
   14    CONTINUE
C
         DO 16 J16=1,KVAR
C
         DO 15 J15=1,KVAR
C
         JCOL = 4 * NT - KVAR + J15
         JROW = 4 * NT - KVAR + J16
C
         JACOL  = INDCOL(JCOL,JROW)
         JAROW  = INDROW(JCOL,JROW)
C
         A(JACOL,JAROW) = PA(J15,J16)
C
   15    CONTINUE
   16    CONTINUE
C
         RETURN
         END
C*DECK C2SD06
C*CALL PROCESS
         SUBROUTINE IDENTB(KVAR,PC)
C        ##########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SD06 PERFORM ROW OPERATIONS REQUIRED TO IMPOSE BOUNDARY           *
*        CONDITIONS IN B                                              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     R   PB(3),   PC(N4NT,3)
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         CALL VZERO(PB,3)
C
         DO 2 J2=1,4*NT
C
         DO 1 J1=1,KVAR
C
         JBROW = J2
C
         PB(J1) = PB(J1) + PC(J2,J1) * B(JBROW)
C
    1    CONTINUE
    2    CONTINUE
C
         DO 3 J3=1,KVAR
C
         JBROW = 4 * NT - KVAR + J3
C
         B(JBROW) = PB(J3)
C
    3    CONTINUE
C
         RETURN
         END
C*DECK C2SD07
C*CALL PROCESS
         SUBROUTINE AWAY(KD,KB1,KB2)
C        ###########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SD07 PERFORMS FOLLOWING LINE-ROW OPERATIONS ON MATRIX A:          *
*                                                                     *
*        1) SET LINE KD TO 0                                          *
*        2) SET ROW KD TO 0                                           *
*        3) SET DIAGONAL TERM A(1,KD) TO 1                            *
*                                                                     *
*        KB1 AND KB2 ARE THE 2 EXTREMAL INDICES OF THE BAND MATRIX.   *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMNUM.inc'
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         INCLUDE 'BNDIND.inc'
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         DO 1 J1=KB1,KB2
C
         IACOL = INDCOL(J1,KD)
         IAROW = INDROW(J1,KD)
C
         A(IACOL,IAROW) = 0
C
    1    CONTINUE
C
         A(1,KD) = 1
C
         RETURN
         END
C*DECK C2SD08
C*CALL PROCESS
         SUBROUTINE CENTER(PC)
C        #####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SD08 EVALUATE THE COFFICIENTS REQUIRED TO IMPOSE BOUNDARY         *
*        CONDITIONS (SEE EQ. (30) IN THE PUBLICATION)                 *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     R   PC(N4NT,3)
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         CALL VZERO(PC,12*NPT)
C
         DO 1 J1=1,NT
C
         ZCOST = COS(CT(J1))
         ZSINT = SIN(CT(J1))
C
         I = 4 * (NUPDWN(J1) - 1)
C
         PC(I+1,1) = 1
         PC(I+2,2) = RHOS(J1) * ZCOST
         PC(I+2,3) = RHOS(J1) * ZSINT
         PC(I+4,2) = DRSDT(J1) * ZCOST - RHOS(J1) * ZSINT
         PC(I+4,3) = DRSDT(J1) * ZSINT + RHOS(J1) * ZCOST
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SE01
C*CALL PROCESS
         SUBROUTINE SOLVIT
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SE01 LEAD GAUSS ELIMINATION FOR THE COMPUTATION OF PSI = A**(-1)*B*
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R   ZB(3),  ZC(N4NT,3)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
***********************************************************************
*                                                                     *
* 1. GAUSS ELIMINATION : SOLVE A * X = B                              *
*                                                                     *
***********************************************************************
C
         CALL DIRECT
         CALL SCOPY(N4NSNT,B,1,CPSILI,1)
c
***********************************************************************
*                                                                     *
* 2. IDENTIFICATIONS AT MESH CENTER                                   *
*                                                                     *
***********************************************************************
C
         DO 1 J1=1,3
C
         ZB(J1)       = B(4*NT-3+J1)
         B(4*NT-3+J1) = 0.
C
    1    CONTINUE
C
         CALL CENTER(ZC)
C
         DO 2 J2=1,NT
C
         I = 4 * (NUPDWN(J2) - 1)
C
         B(I+1) = ZB(1)
         B(I+2) = ZC(I+2,2) * ZB(2) + ZC(I+2,3) * ZB(3)
         B(I+4) = ZC(I+4,2) * ZB(2) + ZC(I+4,3) * ZB(3)
C
    2    CONTINUE
C
***********************************************************************
*                                                                     *
* 3. NEW PSI - SOLUTION                                               *
*                                                                     *
***********************************************************************
C
         CALL SCOPY(N4NSNT,B,1,CPSI,1)
C
***********************************************************************
*                                                                     *
* 4. NEW PSI - SOLUTION IN INVERSE CLOCKWISE NUMEROTATION             *
*                                                                     *
***********************************************************************
C
         DO 3 J3=1,NSTMAX
C
         CPSICL(4*(J3-1)+1) = B(4*(NUPDWN(J3)-1)+1)
         CPSICL(4*(J3-1)+2) = B(4*(NUPDWN(J3)-1)+2)
         CPSICL(4*(J3-1)+3) = B(4*(NUPDWN(J3)-1)+3)
         CPSICL(4*(J3-1)+4) = B(4*(NUPDWN(J3)-1)+4)
C
    3    CONTINUE
C
         RETURN
         END
C*DECK C2SE02
C*CALL PROCESS
         SUBROUTINE DIRECT
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SE02 GAUSS SEIDEL ELIMINATION                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMNUM.inc'
C
***********************************************************************
*                                                                     *
* 1. SOLVE L * Y = B                                                  *
*                                                                     *
***********************************************************************
C
         CALL LYV(A,B,N4NSNT,NP4NST,NBAND,NPBAND)
C
***********************************************************************
*                                                                     *
* 2. SOLVE D * W = Y                                                  *
*                                                                     *
***********************************************************************
C
         CALL DWY(A,B,N4NSNT,NP4NST,NBAND,NPBAND)
C
***********************************************************************
*                                                                     *
* 3. SOLVE LT * X = W                                                 *
*                                                                     *
***********************************************************************
C
         CALL LTXW(A,B,N4NSNT,NP4NST,NBAND,NPBAND)
C
***********************************************************************
*                                                                     *
* 4. RESULT IN B                                                      *
*                                                                     *
***********************************************************************
C
         RETURN
         END
C*DECK C2SE03
C*CALL PROCESS
         SUBROUTINE ERROR1(PSIK,PSIK1)
C        ############################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SE03 COMPUTE ERROR (SEE EQ. (28) IN PUBLICATION)                  *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     I   IS(NPT),       IT(NPT),
     R   PSIK(*),       PSIK1(*),          ZPK(NPT,16),
     R   ZPK1(NPT,16),  ZRES(NPT)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
***********************************************************************
*                                                                     *
* 1. INITIALIZATION                                                   *
*                                                                     *
***********************************************************************
C
         RESIDU = 0.
C
         DO 1 J1=1,NT
C
         IT(J1) = J1
C
    1    CONTINUE
C
***********************************************************************
*                                                                     *
* 2. SCAN OVER ALL CELLS                                              *
*                                                                     *
***********************************************************************
C
         DO 7 J7=1,NS
C
         CALL RESETI(IS,NT,J7)
C
***********************************************************************
*                                                                     *
* 2.1. COMPUTE ALL QUANTITIES TO DEFINE PSI FOR PICARD ITERATION K AND*
*      K + 1                                                          *
*                                                                     *
***********************************************************************
C
         CALL PSICEL(IS,IT,NT,NPT,ZPK,PSIK)
         CALL PSICEL(IS,IT,NT,NPT,ZPK1,PSIK1)
C
***********************************************************************
*                                                                     *
* 2.2 COMPUTE PSI, R AND SOURCE TERM ON INTEGRATION POINTS            *
*                                                                     *
***********************************************************************
C
         DO 4 J4=1,NWGAUS
C
         CALL VZERO(ZRES,NT)
C
         DO 2 J2=1,NT
C
         ZPSIK  = ZPK(J2, 1) * FB(J2,J4, 1,J7) +
     +            ZPK(J2, 2) * FB(J2,J4, 2,J7) +
     +            ZPK(J2, 3) * FB(J2,J4, 3,J7) +
     +            ZPK(J2, 4) * FB(J2,J4, 4,J7) +
     +            ZPK(J2, 5) * FB(J2,J4, 5,J7) +
     +            ZPK(J2, 6) * FB(J2,J4, 6,J7) +
     +            ZPK(J2, 7) * FB(J2,J4, 7,J7) +
     +            ZPK(J2, 8) * FB(J2,J4, 8,J7) +
     +            ZPK(J2, 9) * FB(J2,J4, 9,J7) +
     +            ZPK(J2,10) * FB(J2,J4,10,J7) +
     +            ZPK(J2,11) * FB(J2,J4,11,J7) +
     +            ZPK(J2,12) * FB(J2,J4,12,J7) +
     +            ZPK(J2,13) * FB(J2,J4,13,J7) +
     +            ZPK(J2,14) * FB(J2,J4,14,J7) +
     +            ZPK(J2,15) * FB(J2,J4,15,J7) +
     +            ZPK(J2,16) * FB(J2,J4,16,J7) 
C
         ZPSIK1 = ZPK1(J2, 1) * FB(J2,J4, 1,J7) +
     +            ZPK1(J2, 2) * FB(J2,J4, 2,J7) +
     +            ZPK1(J2, 3) * FB(J2,J4, 3,J7) +
     +            ZPK1(J2, 4) * FB(J2,J4, 4,J7) +
     +            ZPK1(J2, 5) * FB(J2,J4, 5,J7) +
     +            ZPK1(J2, 6) * FB(J2,J4, 6,J7) +
     +            ZPK1(J2, 7) * FB(J2,J4, 7,J7) +
     +            ZPK1(J2, 8) * FB(J2,J4, 8,J7) +
     +            ZPK1(J2, 9) * FB(J2,J4, 9,J7) +
     +            ZPK1(J2,10) * FB(J2,J4,10,J7) +
     +            ZPK1(J2,11) * FB(J2,J4,11,J7) +
     +            ZPK1(J2,12) * FB(J2,J4,12,J7) +
     +            ZPK1(J2,13) * FB(J2,J4,13,J7) +
     +            ZPK1(J2,14) * FB(J2,J4,14,J7) +
     +            ZPK1(J2,15) * FB(J2,J4,15,J7) +
     +            ZPK1(J2,16) * FB(J2,J4,16,J7) 
C
         ZCSURF   = CW(J4) * (CSIG(J7+1)-CSIG(J7)) * (CT(J2+1)-CT(J2))
         ZRES(J2) = ZCSURF * YRST(J2,J4)**2 * RSINT(J7,J4) * 
     *              (ZPSIK1 - ZPSIK)**2
C
    2    CONTINUE
C 
         RESIDU = RESIDU + SSUM(NT,ZRES,1)
C
    4    CONTINUE
    7    CONTINUE
C
         RESIDU = SQRT(RESIDU)
C
         RETURN
         END
C*DECK C2SE04
C*CALL PROCESS
         SUBROUTINE ENERGY
C        #################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SE04  COMPUTE AVERAGED POLOIDAL MAGNETIC FIELD ENERGY OF TORUS    *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     I   IS(NPT),           IT(NPT),
     R   ZPCEL(NPT,16),    
     R   ZDBDS(NPT,16),     ZDBDT(NPT,16),     ZS(NPT),         
     R   ZS1(NPT),          ZS2(NPT),          ZT(NPT),
     R   ZT1(NPT),          ZT2(NPT)        
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
***********************************************************************
*                                                                     *
* 1. INITIALIZATION                                                   *
*                                                                     *
***********************************************************************
C
         ZIBP2 = 0.
C
         DO 1 J1=1,NT
C
         IT(J1)  = J1
         ZT1(J1) = CT(J1)
         ZT2(J1) = CT(J1+1)
C
    1    CONTINUE
C
***********************************************************************
*                                                                     *
* 2. SCAN OVER ALL CELLS                                              *
*                                                                     *
***********************************************************************
C
         DO 8 J8=1,NS
C
         CALL RESETI(IS,NT,J8)
         CALL RESETR(ZS1,NT,CSIG(J8))
         CALL RESETR(ZS2,NT,CSIG(J8+1))
C
***********************************************************************
*                                                                     *
* 2.1. COMPUTE ALL QUANTITIES TO DEFINE PSI ON CELL                   *
*                                                                     *
***********************************************************************
C
         CALL PSICEL(IS,IT,NT,NPT,ZPCEL,CPSICL)
C
***********************************************************************
*                                                                     *
* 2.2 COMPUTE PSI, R AND SOURCE TERM ON INTEGRATION POINTS            *
*                                                                     *
***********************************************************************
C
         DO 7 J7=1,NWGAUS
C
         DO 3 J3=1,NT
C
         ZS(J3) = RSINT(J8,J7)
         ZT(J3) = RTINT(J3,J7)
C
    3    CONTINUE
C
         CALL BASIS2(NT,NPT,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT)
C
         DO 4 J4=1,NT
C
         ZDPDS = ZPCEL(J4, 1) * ZDBDS(J4, 1) +
     +           ZPCEL(J4, 2) * ZDBDS(J4, 2) +
     +           ZPCEL(J4, 3) * ZDBDS(J4, 3) +
     +           ZPCEL(J4, 4) * ZDBDS(J4, 4) +
     +           ZPCEL(J4, 5) * ZDBDS(J4, 5) +
     +           ZPCEL(J4, 6) * ZDBDS(J4, 6) +
     +           ZPCEL(J4, 7) * ZDBDS(J4, 7) +
     +           ZPCEL(J4, 8) * ZDBDS(J4, 8) +
     +           ZPCEL(J4, 9) * ZDBDS(J4, 9) +
     +           ZPCEL(J4,10) * ZDBDS(J4,10) +
     +           ZPCEL(J4,11) * ZDBDS(J4,11) +
     +           ZPCEL(J4,12) * ZDBDS(J4,12) +
     +           ZPCEL(J4,13) * ZDBDS(J4,13) +
     +           ZPCEL(J4,14) * ZDBDS(J4,14) +
     +           ZPCEL(J4,15) * ZDBDS(J4,15) +
     +           ZPCEL(J4,16) * ZDBDS(J4,16)
C
         ZDPDT = ZPCEL(J4, 1) * ZDBDT(J4, 1) +
     +           ZPCEL(J4, 2) * ZDBDT(J4, 2) +
     +           ZPCEL(J4, 3) * ZDBDT(J4, 3) +
     +           ZPCEL(J4, 4) * ZDBDT(J4, 4) +
     +           ZPCEL(J4, 5) * ZDBDT(J4, 5) +
     +           ZPCEL(J4, 6) * ZDBDT(J4, 6) +
     +           ZPCEL(J4, 7) * ZDBDT(J4, 7) +
     +           ZPCEL(J4, 8) * ZDBDT(J4, 8) +
     +           ZPCEL(J4, 9) * ZDBDT(J4, 9) +
     +           ZPCEL(J4,10) * ZDBDT(J4,10) +
     +           ZPCEL(J4,11) * ZDBDT(J4,11) +
     +           ZPCEL(J4,12) * ZDBDT(J4,12) +
     +           ZPCEL(J4,13) * ZDBDT(J4,13) +
     +           ZPCEL(J4,14) * ZDBDT(J4,14) +
     +           ZPCEL(J4,15) * ZDBDT(J4,15) +
     +           ZPCEL(J4,16) * ZDBDT(J4,16) 
C
         ZR     = RSINT(J8,J7) * YRST(J4,J7) * COS(RTINT(J4,J7)) + R0
         ZGRAD2 = (ZDPDS**2 + (ZDPDT / RSINT(J8,J7) -
     -             ZDPDS * YDRSDT(J4,J7) / YRST(J4,J7))**2) /
     /             YRST(J4,J7)**2
         ZW     = CW(J7) * (CSIG(J8+1) - CSIG(J8)) *
     *                     (CT(J4+1) - CT(J4))
         ZJAC   = YRST(J4,J7)**2 * RSINT(J8,J7)
         ZIBP2  = ZIBP2 + ZW * ZJAC * ZGRAD2 / ZR
C
    4    CONTINUE
    7    CONTINUE
    8    CONTINUE
C
         WMAGP = .5  * ZIBP2
C
         RETURN
         END
C*DECK C2SE05
C*CALL PROCESS
         SUBROUTINE SMOOTH
C        #################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SE05  SMMOTH BICUBIC HERMITE SOLUTION WITH BICUBIC SPLINES        *
*         (SEE APPENDIX C OF PUBLICATION)                             *
*                                                                     *
***********************************************************************
C
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'DECLAR.inc'
C
         DIMENSION X1(NPT,NSP1),    Y1(NPT,NSP1),    YD2P1(NPT,NSP1),
     &             X2(2*NSP1,NPT),  Y2(2*NSP1,NPT),  YD2P2(2*NSP1,NPT),
     &             ZA1(NPT,NSP1),   ZB1(NPT,NSP1),   ZWORK(NPT+2*NSP1),
     &             ZA2(2*NSP1,NPT), ZB2(2*NSP1,NPT), ZC2(2*NSP1,NPT),
     &             YP1(NPT),        YPN(NPT)
C
         DO 2 J2=1,NS1
            DO 1 J1=1,NT
C
            X1(J1,J2) = CSIG(J2)
            Y1(J1,J2) = CPSICL(4*((J2-1)*NT+J1)-3)
C
 1          CONTINUE
 2       CONTINUE
C
         DO 3 J3=1,NT
C
            YP1(J3) = CPSICL(4*J3-2)
            YPN(J3) = CPSICL(4*(NS*NT+J3)-2)
C
 3       CONTINUE
C
         CALL MSPLINE(X1,Y1,NS1,NPT,NT,YP1,YPN,YD2P1,ZA1,ZB1,ZWORK)
C
         ZERDPS = 0.
C
         DO 5 J5=2,NS
            DO 4 J4=1,NT
C
               IP = 4 * ((J5 - 1) * NT + J4) - 2
C
               ZDUM      = CPSICL(IP)
               CPSICL(IP) = 
     &            (Y1(J4,J5) - Y1(J4,J5-1)) / 
     &            (X1(J4,J5) - X1(J4,J5-1)) +
     &            YD2P1(J4,J5-1) * (X1(J4,J5) - X1(J4,J5-1)) / 6. +
     &            YD2P1(J4,J5) * (X1(J4,J5) - X1(J4,J5-1)) / 3.
C
               ZERDPS = ZERDPS + (ZDUM - CPSICL(IP))**2
C
 4          CONTINUE
 5       CONTINUE
C
         ZERDPS = SQRT(ZERDPS) / DFLOAT(NS1 * NT)
C
         DO 7 J7=1,NT
            DO 6 J6=1,NS1
C
               X2(J6,J7)     = CT(J7)
               X2(NS1+J6,J7) = CT(J7)
               Y2(J6,J7)     = CPSICL(4*((J6-1)*NT+J7)-3)
               Y2(NS1+J6,J7) = CPSICL(4*((J6-1)*NT+J7)-2)
C
 6          CONTINUE
 7       CONTINUE
C
         CALL MSPLCY(X2,Y2,NT,2*NSP1,2*NS1,RC2PI,YD2P2,
     &               ZA2,ZB2,ZC2,ZWORK)
C
         ZERDPT = 0.
         ZERD2P = 0.
C
         DO 8 J8=2,NS1
C
            IP1   = 4 * (J8 - 1) * NT + 3
            IP2   = 4 * (J8 - 1) * NT + 4
            ZDUM1 = CPSICL(IP1)
            ZDUM2 = CPSICL(IP2)
C
            CPSICL(IP1) = 
     &         (Y2(J8,2) - Y2(J8,1)) / (X2(J8,2) - X2(J8,1)) -
     &         YD2P2(J8,1) * (X2(J8,2) - X2(J8,1)) / 3. -
     &         YD2P2(J8,2) * (X2(J8,2) - X2(J8,1)) / 6.
C
            CPSICL(IP2) = 
     &         (Y2(NS1+J8,2) - Y2(NS1+J8,1)) / 
     &         (X2(NS1+J8,2) - X2(NS1+J8,1)) -
     &         YD2P2(NS1+J8,1) * (X2(NS1+J8,2) - X2(NS1+J8,1)) / 3. -
     &         YD2P2(NS1+J8,2) * (X2(NS1+J8,2) - X2(NS1+J8,1)) / 6.
C
            ZERDPT = ZERDPT + (ZDUM1 - CPSICL(IP1))**2
            ZERD2P = ZERD2P + (ZDUM2 - CPSICL(IP2))**2
C
 8       CONTINUE
C
         DO 10 J10=2,NT
C
            IT = J10 - 1
C
            DO 9 J9=2,NS1
C
               IS    = NS1 + J9
               IP1   = 4 * ((J9 - 1) * NT + J10) - 1
               IP2   = 4 * ((J9 - 1) * NT + J10)
               ZDUM1 = CPSICL(IP1)
               ZDUM2 = CPSICL(IP2)
C
               CPSICL(IP1) = 
     &            (Y2(J9,J10) - Y2(J9,IT)) / (X2(J9,J10) - X2(J9,IT)) +
     &            YD2P2(J9,IT) * (X2(J9,J10) - X2(J9,IT)) / 6. +
     &            YD2P2(J9,J10) * (X2(J9,J10) - X2(J9,IT)) / 3.
C
               CPSICL(IP2) = 
     &            (Y2(IS,J10) - Y2(IS,IT)) / (X2(IS,J10) - X2(IS,IT)) +
     &            YD2P2(IS,IT) * (X2(IS,J10) - X2(IS,IT)) / 6. +
     &            YD2P2(IS,J10) * (X2(IS,J10) - X2(IS,IT)) / 3.
C
            ZERDPT = ZERDPT + (ZDUM1 - CPSICL(IP1))**2
            ZERD2P = ZERD2P + (ZDUM2 - CPSICL(IP2))**2
C
 9          CONTINUE
 10      CONTINUE
C
         ZERDPT = SQRT(ZERDPT) / DFLOAT(NS1 * NT)
         ZERD2P = SQRT(ZERD2P) / DFLOAT(NS1 * NT)
C
         DO 11 J11=1,NS1*NT
C
            CPSI(4*NUPDWN(J11)-3) = CPSICL(4*J11-3)
            CPSI(4*NUPDWN(J11)-2) = CPSICL(4*J11-2)
            CPSI(4*NUPDWN(J11)-1) = CPSICL(4*J11-1)
            CPSI(4*NUPDWN(J11)) = CPSICL(4*J11)
C
 11      CONTINUE
C
         END
C*DECK C2SE06
C*CALL PROCESS
         SUBROUTINE CONVER(K,KCON)
C        #########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SE06 CONVERGENCE TEST : KCON = 0  NO CONVERGENCE                  *
*                           KCON = 1  CONVERGENCE OVER NON-           *
*                                     LINEARITY                       *
*                           KCON = 2  CONVERGENCE OVER CURRENT        *
*                                     PROFILE                         *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMSOL.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         GOTO (1,2) K
C
    1    CONTINUE
C
***********************************************************************
*                                                                     *
* 1. CONVERGENCE TEST ON RESIDU OF ITERATION OVER NONLINEARITY        *
*                                                                     *
***********************************************************************
C
         IF (RESIDU .LT. CEPS) KCON = 1
C
         RETURN
C
    2    CONTINUE
C
***********************************************************************
*                                                                     *
* 2. CONVERGENCE TEST ON RESIDU OF ITERATION OVER MAPPING             *
*                                                                     *
***********************************************************************
C
         IF (RESMAP .LT. 100. * CEPS) THEN
C
            KCON = 2
            RETURN
C
         ENDIF
C
         KCON   = 0
C
         RETURN
         END
C*DECK C2SF01
C*CALL PROCESS
         SUBROUTINE NOREPT(KN,KSHIFT)
C        ############################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SF01 EQUILIBRIUM TRANSFORMATION (SEE SECTIONS 2.3, 5.3 AND 6.4.4  *
*        IN PUBLICATION)                                              *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   ZQ(2*NPISO),   ZS(2*NPISO), zder(2*NPISO)
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         RRAXIS = RMAG
         RZAXIS = RZMAG
C
***********************************************************************
*                                                                     *
* TOKAMAK EQUILIBRIUM                                                 *
*                                                                     *
***********************************************************************
C
         IF (NRFP .EQ. 0) THEN
C
         IF (NCSCAL .EQ. 1) THEN
C
            DO 1 J1=1,KN
C
            ZS(J1) = SQRT(1. - PSIISO(J1) / SPSIM)
            ZQ(J1) = .5 * TMF(J1) * CIDQ(J1) / CPI
C
    1       CONTINUE
cab
      IF (NQMIN .LE. 0) GOTO 103
cab   csspec is redefined only if nqmin > 0
cab
      do 101 j = kn-1,2,-1
      zder(j) = tmf(j)*(tmf(j+1)-tmf(j-1))*cidq(j)**2
     &        + (2.*cpi*qspec)**2*(cidq(j+1)-cidq(j-1))/cidq(j)
      if (zder(j).lt.0.) goto 102
 101  continue
      csspec=0.
      goto 103
 102  csspec=zs(j)-zder(j)*(zs(j+1)-zs(j))/(zder(j+1)-zder(j))
 103  WRITE(*,*) ' NQMIN = ',NQMIN,'   CSSPEC =',CSSPEC
cab
C
            ICS = ISRCHFGE(KN,ZS,1,CSSPEC) - 1
C
            IF (ICS .LT. 2)      ICS = 2
            IF (ICS .GE. KN - 1) ICS = KN - 2
C
            QICS = FCCCC0(ZQ(ICS-1),ZQ(ICS),ZQ(ICS+1),ZQ(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS),ZS(ICS+1),ZS(ICS+2),CSSPEC)
            TICS = FCCCC0(TMF(ICS-1),TMF(ICS),TMF(ICS+1),TMF(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS), ZS(ICS+1),ZS(ICS+2),CSSPEC)
            T0   = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                    ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
C
            ZCSHFT = TICS**2 * ((QSPEC / QICS)**2 - 1.)
C
            CALL TSHIFT(ZCSHFT,KN)
C
            IF (NTMF0 .EQ. 0) THEN
C
               SCALE = 1. / TMF(KN)
C
            ELSE IF (NTMF0 .EQ. 1) THEN
C
               SCALE = 1. / T0
C
            ENDIF
C
            CALL PRNORM(SCALE,KN)
C
            IF (NSTTP .EQ. 3 .AND. (NOPT .EQ. 0 .OR. (NOPT .EQ. 1 .AND.
     &          (NBLOPT .NE. 0 .AND. CPRESS .NE. 1.)))) THEN
C
               IF (NFUNC .EQ. 1) THEN
C
                  CALL SSCAL(NSOUR,SCALE,AT,1)
C
               ELSE IF (NFUNC .EQ. 2) THEN
C
                  CALL SSCAL(5,SCALE,AT(3),1)
                  CALL SSCAL(5,SCALE,AT2(3),1)
                  CALL SSCAL(5,SCALE,AT3(3),1)
C
                  AT4(3) = SCALE * AT4(3)
C
               ELSE IF (NFUNC .EQ. 3) THEN
C
                   scale = scale**scexp
                   AT(1) = SCALE * AT(1)
C
               ELSE IF (NFUNC .EQ. 4) THEN
C
                  CALL SSCAL(NPPF+1,SCALE,RFUN,1)
C
               ENDIF
C
               SCALAC = SCALAC * SCALE
C
               CALL RVAR('SCALE            ',SCALE)
               CALL RVAR('ACCUMULATED SCALE',SCALAC)
C
            ENDIF
C
         ELSE IF (NCSCAL .EQ. 2) THEN
C
            T0 = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                  PSIISO(1),PSIISO(2),PSIISO(3),PSIISO(4),SPSIM)
C
C            IF (NSTTP .LE. 2 .OR. (NSTTP .EQ. 3 .AND. NOPT .EQ. 1 .AND.
C     &          (NBLOPT .EQ. 0 .OR. CPRESS .EQ. 1.))) THEN
C
            IF (NSURF .EQ. 1) THEN
C
               SCALE = 1.
C
            ELSE
C
               SCALE = CURRT / CUROLD
C
            ENDIF
C
            IF (NSTTP.LE.2 .OR. NSTTP.EQ.4) THEN
                CALL PRNORM(SCALE,KN)
C
               IF (NTMF0 .EQ. 0) THEN
C
                  ZCSHFT = 1. - TMF(KN)**2
C
               ELSE IF (NTMF0 .EQ. 1) THEN
C
                  ZCSHFT = 1. - T0**2
C
               ENDIF
C
               CALL TSHIFT(ZCSHFT,KN)
C
            ELSE 
C
C               PRINT*,'NOREPT : NSTTP=3; NCSCAL=2; NOPT=',NOPT,';'
C               PRINT*,'NBLOPT=',NBLOPT,'; CPRESS=',CPRESS,';'
C               PRINT*,'THIS OPTION IS NOT POSSIBLE'
C               STOP               
               scale = scale**scexp
               AT(1) = SCALE  * AT(1)
               SCALAC = SCALAC * SCALE
               CALL RVAR('SCALE            ',SCALE)
               CALL RVAR('ACCUMULATED SCALE',SCALAC)
C
            ENDIF
C
         ELSE IF (NCSCAL .EQ. 3) THEN
C
            ZS(1) = (1. - PSIISO(1) / SPSIM) * ABS(SPSIM) * CIDR(1)
            ZQ(1) = .5 * TMF(1) * CIDQ(1) / CPI
C
            DO 2 J2=2,KN
C
            ZX1 = SQRT(1. - PSIISO(J2-1) / SPSIM)
            ZX2 = SQRT(1. - PSIISO(J2  ) / SPSIM)
C
            ZS(J2) = ZS(J2-1) + ABS(SPSIM) * (ZX2 - ZX1) *
     *               (ZX2 * CIDR(J2) + ZX1 * CIDR(J2-1))
            ZQ(J2) = .5 * TMF(J2) * CIDQ(J2) / CPI
C
    2       CONTINUE
C
            DO 3 J3=1,KN
C
            ZS(J3) = SQRT(ZS(J3) / ZS(KN))
C
    3       CONTINUE
C
            ICS = ISRCHFGE(KN,ZS,1,CSSPEC) - 1
C
            IF (ICS .LT. 2)      ICS = 2
            IF (ICS .GE. KN - 1) ICS = KN - 2
C
            QICS = FCCCC0(ZQ(ICS-1),ZQ(ICS),ZQ(ICS+1),ZQ(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS),ZS(ICS+1),ZS(ICS+2),CSSPEC)
            TICS = FCCCC0(TMF(ICS-1),TMF(ICS),TMF(ICS+1),TMF(ICS+2),
     ,                    ZS(ICS-1),ZS(ICS), ZS(ICS+1),ZS(ICS+2),CSSPEC)
            T0   = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                    ZS(1),ZS(2),ZS(3),ZS(4),0.)
C
            ZCSHFT = TICS**2 * ((QSPEC / QICS)**2 - 1.)
C
            CALL TSHIFT(ZCSHFT,KN)
C
            IF (NTMF0 .EQ. 0) THEN
C
               SCALE = 1. / TMF(KN)
C
            ELSE IF (NTMF0 .EQ. 1) THEN
C
               SCALE = 1. / T0
C
            ENDIF
C
            CALL PRNORM(SCALE,KN)
C
            IF (NSTTP .EQ. 3 .AND. (NOPT .EQ. 0 .OR. (NOPT .EQ. 1 .AND.
     &          (NBLOPT .NE. 0 .AND. CPRESS .NE. 1.)))) THEN
C
               IF (NFUNC .EQ. 1) THEN
C
                  CALL SSCAL(NSOUR,SCALE,AT,1)
C
               ELSE IF (NFUNC .EQ. 2) THEN
C
                  CALL SSCAL(5,SCALE,AT(3),1)
                  CALL SSCAL(5,SCALE,AT2(3),1)
                  CALL SSCAL(5,SCALE,AT3(3),1)
C
                  AT4(3) = SCALE * AT4(3)
C
               ELSE IF (NFUNC .EQ. 3) THEN
C
                  scale = scale**scexp
                  AT(1) = SCALE * AT(1)
C
               ELSE IF (NFUNC .EQ. 4) THEN
C
                  CALL SSCAL(NPPF+1,SCALE,RFUN,1)
C
               ENDIF
C
               SCALAC = SCALAC * SCALE
C
               CALL RVAR('SCALE            ',SCALE)
               CALL RVAR('ACCUMULATED SCALE',SCALAC)
C
            ENDIF
         ENDIF
C
***********************************************************************
*                                                                     *
* REVERSED FIELD PINCH EQUILIBRIUM                                    *
* NOTE THAT FOR RFP, ONE CAN NOT SCALE EQUILIBRIUM BY SHIFTING T**2   *
* BECAUSE ONE CAN NOT SPECIFY T(0) AND CURRT AT THE SAME TIME FOR RFP *
* HENCE, LET'S FIX T(0)=1, ANd ALLOW CURRT TO BE SCALED               *
***********************************************************************
C
         ELSE IF (NRFP .EQ. 1) THEN
C
            IF (NSURF .EQ. 1) THEN
C
               SCALE = 1.
C
            ELSE
C
               SCALE = CURRT / CUROLD
C
            ENDIF
C
            T0 = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                  PSIISO(1),PSIISO(2),PSIISO(3),PSIISO(4),SPSIM)
            SCALE = 1./T0
C
            CALL PRNORM(SCALE,KN)
C
         ENDIF
C
         PSI0   = SPSIM
         CPSRF  = ABS(SPSIM)
C
         IF (NRSCAL .EQ. 1) CALL RSCALE(RRAXIS,KN)
C
         IF (KSHIFT.EQ.0) RETURN
C
         DO 5 J5=1,NSNT
C
         CPSICL(4*J5-3) = CPSICL(4*J5-3) + CPSRF
C
    5    CONTINUE
C
         DO 6 J6=1,KN
C
         PSIISO(J6) = PSIISO(J6) + CPSRF
C
    6    CONTINUE
C
         RETURN
         END
C*DECK C2SF02
C*CALL PROCESS
         SUBROUTINE RSCALE(PR,KN)
C        ########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SF02 SCALE EQUILIBRIUM SO THAT PR = 1. (SEE EQ. (34) IN           *
*        PUBLICATION)                                                 *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         RMAG   = RMAG / PR
         RZMAG  = RZMAG / PR
         R0     = R0 / PR
         RZ0    = RZ0 / PR
         RC     = RC / PR
         RZ0C   = RZ0C / PR
         R0W    = R0W / PR
         RZ0W   = RZ0W / PR
         SPSIM  = SPSIM / PR
         CPSRF  = CPSRF / PR
C
C     ADAPT R0EXP SO THAT R0 STAYS SAME VALUE IN PHYSICAL UNIT
C     AND B0EXP SO THAT R0EXP*B0EXP=CST (SO TOT. CURRENT STAYS SAME
C     AS R SCALE AT CONSTANT CURRENT IS DONE HERE)
C
         R0EXP = R0EXP * PR
         B0EXP = B0EXP / PR
C
         BPS( 1) = BPS( 1) / PR
         BPS( 2) = BPS( 2) / PR
         BPS( 3) = BPS( 3) / PR
         BPS( 6) = BPS( 6) / PR
         BPS(12) = BPS(12) / PR
C
         IF (NSURF .EQ. 6) THEN
C
C NONCONFORMAL WALL
            DO J1=1,NWBPS
               CALL SSCAL(NBPS,1./PR,RRBPS(1,J1),1)
               CALL SSCAL(NBPS,1./PR,RZBPS(1,J1),1)
            ENDDO
            CALL BNDSPL
C
         ENDIF
C
         IF (NSURF .EQ. 7) THEN
            CALL SSCAL(NFOURPB,1./PR,BPSCOS,1)
            CALL SSCAL(NFOURPB,1./PR,BPSSIN,1)
            ALZERO=ALZERO/PR
         ENDIF
C
         CALL SSCAL(N4NSNT,1./PR,CPSI,1)
         CALL SSCAL(N4NSNT,1./PR,CPSICL,1)
         CALL SSCAL(KN,1./PR,PSIISO,1)
         CALL SSCAL(KN,PR,TTP,1)
         CALL SSCAL(KN,PR**2,CPR,1)
         CALL SSCAL(KN,PR**2*PR,CPPR,1)
C
         RETURN
         END
C*DECK C2SF03
C*CALL PROCESS
         SUBROUTINE TSHIFT(PC,KN)
C        ########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SF03 SHIFT T PROFILE. (SEE EQ. (33) IN PUBLICATION)               *
* FOR TOKAMAK ONLY                                                    *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         DO 1 J1=1,KN
C
         TMF(J1) = SQRT(TMF(J1)**2 + PC)
C
    1    CONTINUE
C
         T0 = SQRT(T0**2 + PC)
C
         RETURN
         END
C*DECK C2SF04
C*CALL PROCESS
         SUBROUTINE PRNORM(PC,KN)
C        ########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SF04 SCALE EQUILIBRIUM. (SEE EQ. (32) IN PUBLICATION)             *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL SSCAL(N4NSNT,PC,CPSI,1)
         CALL SSCAL(N4NSNT,PC,CPSICL,1)
         CALL SSCAL(KN,PC,PSIISO,1)
         CALL SSCAL(KN,PC,TMF,1)
         CALL SSCAL(KN,PC,TTP,1)
         CALL SSCAL(KN,PC**2,CPR,1)
         CALL SSCAL(KN,PC,CPPR,1)
C
         SPSIM = PC * SPSIM
         T0    = PC * T0
C
         RETURN
         END
C*DECK C2SG01
C*CALL PROCESS
         SUBROUTINE TEST
C        ###############
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SG01 FOR SOLOVEV EQUILIBRIUM (NSURF = 1), COMPUTE NUMERICAL AND   *
*        ANALYTICAL VALUES OF PSI AT (CSIG(JS),THETA(JT) NODES,       *
*        AS WELL AS ABSOLUTE DIFFERENCES AND RESIDU                   *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMVEV.inc'
C
         DIMENSION
     R   ZND(NPT),  ZNDS(NPT),  ZNDST(NPT),  ZNDT(NPT)
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         INCLUDE 'SOLOV.inc'
C
***********************************************************************
*                                                                     *
* 1. TEST IF NSURF = 1                                                *
*                                                                     *
***********************************************************************
C
         IF (NSURF .NE. 1) RETURN
C
***********************************************************************
*                                                                     *
* 2.1. INITIALIZATION                                                 *
*                                                                     *
***********************************************************************
C
         ZND1   = 0.
         ZNDS1  = 0.
         ZNDT1  = 0.
         ZNDST1 = 0.
C
***********************************************************************
*                                                                     *
* 2.2. SCAN OVER ALL INTERVALS                                        *
*                                                                     *
***********************************************************************
C
         DO 2 J2=1,NS1
C
         DO 1 J1=1,NT
C
         I = (J2 - 1) * NT + J1
C
         ZR = CSIG(J2) * RHOS(J1) * COS(CT(J1)) + R0
         ZZ = CSIG(J2) * RHOS(J1) * SIN(CT(J1)) + RZ0
C
         ZDRDS = RHOS(J1) * COS(CT(J1))
         ZDZDS = RHOS(J1) * SIN(CT(J1))
         ZDRDT = CSIG(J2) * (DRSDT(J1) * COS(CT(J1)) -
     -                       RHOS(J1) * SIN(CT(J1)))
         ZDZDT = CSIG(J2) * (DRSDT(J1) * SIN(CT(J1)) +
     +                       RHOS(J1) * COS(CT(J1)))
         ZD2R  = DRSDT(J1) * COS(CT(J1)) - RHOS(J1) * SIN(CT(J1))
         ZD2Z  = DRSDT(J1) * SIN(CT(J1)) + RHOS(J1) * COS(CT(J1))
C
         CPSI1T(I) = SPS(ZR,ZZ)
         DPDSTH(I) = DSPSDR(ZR,ZZ) * ZDRDS + DSPSDZ(ZR,ZZ) * ZDZDS
         DPDTTH(I) = DSPSDR(ZR,ZZ) * ZDRDT + DSPSDZ(ZR,ZZ) * ZDZDT
         D2PSTT(I) = D2SPSR(ZR,ZZ) * ZDRDT * ZDRDS +
     +               D2SPSZ(ZR,ZZ) * ZDZDT * ZDZDS +
     +               DSPSRZ(ZR,ZZ) * (ZDRDT * ZDZDS +
     +                                ZDRDS * ZDZDT) +
     +               DSPSDR(ZR,ZZ) * ZD2R +
     +               DSPSDZ(ZR,ZZ) * ZD2Z
         DPDSNU(I) = CPSICL(4*(I-1)+2)
         DPDTNU(I) = CPSICL(4*(I-1)+3)
         D2PSTN(I) = CPSICL(4*(I-1)+4)
C
         DIFFP(I)  = CPSI1T(I) - CPSICL(4*(I-1)+1)
         DIFFDS(I) = DPDSTH(I) - DPDSNU(I)
         DIFFDT(I) = DPDTTH(I) - DPDTNU(I)
         DIFFST(I) = D2PSTT(I) - D2PSTN(I)
C
         ZND(J1)   = DIFFP(I)**2
         ZNDS(J1)  = DIFFDS(I)**2
         ZNDT(J1)  = DIFFDT(I)**2
         ZNDST(J1) = DIFFST(I)**2
C
    1    CONTINUE
C
         ZND1   = ZND1   + SSUM(NT,ZND,1)
         ZNDS1  = ZNDS1  + SSUM(NT,ZNDS,1)
         ZNDT1  = ZNDT1  + SSUM(NT,ZNDT,1)
         ZNDST1 = ZNDST1 + SSUM(NT,ZNDST,1)

    2    CONTINUE
C
         RESPSI = SQRT(ZND1)   / NSTMAX
         RESDPS = SQRT(ZNDS1)  / NSTMAX
         RESDPT = SQRT(ZNDT1)  / NSTMAX
         RESDST = SQRT(ZNDST1) / NSTMAX
C
         RETURN
         END
C*DECK C2SG02
C*CALL PROCESS
         SUBROUTINE SOLOVEV
C        ##################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SG02 COMPUTE ANALYTIC PROJECTION ON BICUBIC HERMITE ELEMENTS OF   *
*        SOLOVEV EQUILIBRIUM SOLUTION                                 *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         INCLUDE 'SOLOV.inc'
C
C----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----
C
         DO 2 J2=1,NS1
C
         DO 1 J1=1,NT
C
         I = (J2 - 1) * NT + J1
C
         ZR = CSIG(J2) * RHOS(J1) * COS(CT(J1)) + R0
         ZZ = CSIG(J2) * RHOS(J1) * SIN(CT(J1)) + RZ0
C
         ZDRDS = RHOS(J1) * COS(CT(J1))
         ZDZDS = RHOS(J1) * SIN(CT(J1))
         ZDRDT = CSIG(J2) * (DRSDT(J1) * COS(CT(J1)) -
     -                       RHOS(J1) * SIN(CT(J1)))
         ZDZDT = CSIG(J2) * (DRSDT(J1) * SIN(CT(J1)) +
     +                       RHOS(J1) * COS(CT(J1)))
         ZD2R  = DRSDT(J1) * COS(CT(J1)) - RHOS(J1) * SIN(CT(J1))
         ZD2Z  = DRSDT(J1) * SIN(CT(J1)) + RHOS(J1) * COS(CT(J1))
C
         IF (J2 .NE. NS1) THEN
C
            CPSICL(4*(I-1)+1) = SPS(ZR,ZZ) - SPSI0
            CPSICL(4*(I-1)+3) = DSPSDR(ZR,ZZ) * ZDRDT +
     +                          DSPSDZ(ZR,ZZ) * ZDZDT
C
         ELSE
C
            CPSICL(4*(I-1)+1) = 0.
            CPSICL(4*(I-1)+3) = 0.
C
         ENDIF
C
         CPSICL(4*(I-1)+2) = DSPSDR(ZR,ZZ) * ZDRDS +
     +                       DSPSDZ(ZR,ZZ) * ZDZDS
         CPSICL(4*(I-1)+4) = D2SPSR(ZR,ZZ) * ZDRDT * ZDRDS +
     +                       D2SPSZ(ZR,ZZ) * ZDZDT * ZDZDS +
     +                       DSPSRZ(ZR,ZZ) * (ZDRDT * ZDZDS +
     +                                        ZDRDS * ZDZDT) +
     +                       DSPSDR(ZR,ZZ) * ZD2R +
     +                       DSPSDZ(ZR,ZZ) * ZD2Z
C
    1    CONTINUE
    2    CONTINUE
C
         SPSIM = - SPSI0
         RMAG  = RC
         RZMAG = 0.
C
         CALL OUTPUT(5)
C
         RETURN
         END
C*DECK C2SM01
C*CALL PROCESS
         SUBROUTINE MAPPIN(K)
C        ####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SM01 LEAD MAPPINGS FOR DIFFERENT CODES LINKED TO CHEASE           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   ZS(2*NPISO),    ZTMF(2*NPISO), ZTPR(2*NPISO), ZCPR(2*NPISO),
     R   ZCPPR(2*NPISO), ZQ(2*NPISO),   ZQP(2*NPISO)
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
CMSC
         WRITE(*,'(" MAPPIN  K= ",I4)')K
         GOTO(20,30,50) K
C
   20    CONTINUE
C
         DO 21 J21=1,NPSI1
C
         CALL SURFACE(J21,SIGPSI(1,J21),TETPSI(1,J21),WGTPSI(1,J21),
     ,                CSM(J21))
C
   21    CONTINUE
C
         IF (NSURF .NE. 1) THEN
C
            CP0   = FCCCC0(CPR(1),CPR(2),CPR(3),CPR(4),
     ,                     CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
            DPDP0 = FCCCC0(CPPR(1),CPPR(2),CPPR(3),CPPR(4),
     ,                     CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
            T0    = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                     CSM(1),CSM(2),CSM(3),TMF(4),RC0P)
            DTTP0 = FCCCC0(TTP(1),TTP(2),TTP(3),TTP(4),
     ,                     CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
C
         ELSE
C
            T0    = 1.
            DTTP0 = 0.
            CP0   = CPP * SPSIM
            DPDP0 = CPP
C
         ENDIF
C
         Q0    = FCCCC0(QPSI(1),QPSI(2),QPSI(3),QPSI(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
         RIPR0 = FCCCC0(RIPR(1),RIPR(2),RIPR(3),RIPR(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
         RJDTB0= FCCCC0(RJDOTB(1),RJDOTB(2),RJDOTB(3),RJDOTB(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
         CPSI0 = FCCCC0(CP(1),CP(2),CP(3),CP(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
         DQDP0 = FCCCC0(CDQ(1),CDQ(2),CDQ(3),CDQ(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
         CPDP0 = FCCCC0(CPDP(1),CPDP(2),CPDP(3),CPDP(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
C
         IF (NPROFZ .EQ. 1) THEN
C 
         TEMP0 = FCCCC0(TEMPER(1),TEMPER(2),TEMPER(3),TEMPER(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
         DENS0 = FCCCC0(DENSTY(1),DENSTY(2),DENSTY(3),DENSTY(4),
     ,                  CSM(1),CSM(2),CSM(3),CSM(4),RC0P)
C
         ENDIF
C
         CALL GLOQUA(CSM,CS,NPSI1,1)
C
         RETURN
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
   30    CONTINUE
C
         ASPCTR = 1. / ASPCT
C
         DO 31 J31=1,NPSI
C
         ZS(2*(J31-1)+1) = CSM(J31)
         ZS(2*J31      ) = CS(J31+1)
C
   31    CONTINUE
C
         DO 32 J32=1,2*NPSI
C
         CALL SURFACE(J32,SIGPSI(1,J32),TETPSI(1,J32),WGTPSI(1,J32),
     ,                ZS(J32))
C
   32    CONTINUE
C
         IF (NSURF .NE. 1) THEN
C
            CP0   = FCCCC0(CPR(1),CPR(2),CPR(3),CPR(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
            DPDP0 = FCCCC0(CPPR(1),CPPR(2),CPPR(3),CPPR(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
            T0    = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
            DTTP0 = FCCCC0(TTP(1),TTP(2),TTP(3),TTP(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
C
         ELSE
C
            T0    = 1.
            DTTP0 = 0.
            CP0   = CPP * SPSIM
            DPDP0 = CPP
C
         ENDIF
C
         Q0    = FCCCC0(QPSI(1),QPSI(2),QPSI(3),QPSI(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         RIPR0 = FCCCC0(RIPR(1),RIPR(2),RIPR(3),RIPR(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         RJDTB0= FCCCC0(RJDOTB(1),RJDOTB(2),RJDOTB(3),RJDOTB(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         CPSI0 = FCCCC0(CP(1),CP(2),CP(3),CP(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         DQDP0 = FCCCC0(CDQ(1),CDQ(2),CDQ(3),CDQ(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         CPDP0 = FCCCC0(CPDP(1),CPDP(2),CPDP(3),CPDP(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
C
         IF (NPROFZ .EQ. 1) THEN
C
         TEMP0 = FCCCC0(TEMPER(1),TEMPER(2),TEMPER(3),TEMPER(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         DENS0 = FCCCC0(DENSTY(1),DENSTY(2),DENSTY(3),DENSTY(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
C
         ENDIF
C
         IF (NIDEAL .EQ. 0) THEN
C
C     NUMERICAL INT. TRANSFORM
C
            DO 34 J34=1,2*NPSI
C
            CALL GIJLIN(J34,ZS(J34))
C
            IF (NFFTOPT .EQ. 0) THEN
C
              DO 33 J33=1,MSMAX
C
                CALL FOURIER(J34,SIGPSI(1,J34),TETPSI(1,J34),
     ,                       WGTPSI(1,J34),J33-1)
C
                IF (J34 .EQ. 2*NPSI .AND. REXT .GT. 1.) THEN
C
C NONCONFORMAL WALL: VACUUMNW
c OTHERWISE: USE VACUUM
                  CALL VACUUMNW(J33-1)
C
                ENDIF
C
 33           CONTINUE
C
              CALL GIJREA(J34,ZS(J34))
C
            ELSE
C
C     REAL SPACE OUTPUT
C
              CALL FOURFFT(J34,MSMAX)
              IF (J34 .EQ. 2*NPSI .AND. REXT .GT. 1.)
     +                                     CALL VACUFFT(MSMAX)
C
            ENDIF
C
  34        CONTINUE


C        -----------------------------------------------
C        LIUYQ, 2005-01-29
C        OUTPUT RZ-COORDINATES IN REAL SPACE

         IF (1.EQ.0) THEN
         OPEN(UNIT=26,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='RMZM_R')
         WRITE(26,1021) NMGAUS*NT1,2*NPSI+1
         WRITE(26,1022) 0.0,0.0
         DO I=1,NPSI
            WRITE(26,1022) CSM(I),CSM(I)
            WRITE(26,1022) CS(I+1),CS(I+1)
         ENDDO
         DO J=1,NMGAUS*NT1
            WRITE(26,1022) CHIISO(J),CHIISO(J)
         ENDDO
         DO J=1,NMGAUS*NT1
            WRITE(26,1022) R0,RZ0
         ENDDO
         DO I=1,2*NPSI
            DO J=1,NMGAUS*NT1
               WRITE(26,1022) RRISO(J,I),RZISO(J,I)
            ENDDO
         ENDDO
 1021    FORMAT(I6,1X,I6)
 1022    FORMAT(E15.8,1X,E15.8)
         CLOSE(UNIT=26)
         ENDIF
C        -----------------------------------------------

         IF (NSURF.EQ.6) THEN 
           WRITE(*,*) 'NVACUUMRNW=',NVACUUMRNW
           IF (NVACUUMRNW.EQ.0) CALL VACUUMR
           IF (NVACUUMRNW.EQ.1) CALL VACUUMRNW
           IF (NVACUUMRNW.EQ.2) CALL VACUUMRNW2
         ELSE
           CALL VACUUMR
         ENDIF
C
            DO 36 J36=1,NSMAX
C
            DO 35 J35=1,MSMAX
C
            RM(J35,J36) = DFLOAT(J35 - 1)
            RN(J36)     = DFLOAT(J36 - 1)
C
   35       CONTINUE
   36       CONTINUE
CMSC    output to INP1
             CALL OUTNVW

         ELSE IF (NIDEAL .EQ. 3) THEN
CMSC         moved to nideal=0
C            CALL OUTNVW
CJEM         keep this option...
            CALL OUTNVW
C
         ELSE IF (NIDEAL .EQ. 5) THEN
C
            DO 37 J37=1,2*NPSI
               CALL TPSI(2*NPSI,J37)
               CALL OUTXT(J37,ZS(J37))
 37         CONTINUE
C
         ENDIF
C
         DO 40 J40=1,NPSI1
C
         IF (J40 .EQ. NPSI1) THEN
            I40 = 2*NPSI
         ELSE
            I40 = (J40-1)*2+1
         ENDIF
C
         PSIISO(J40) = PSIISO(I40)
         TMF(J40)    = TMF(I40)
         TTP(J40)    = TTP(I40)
         CPR(J40)    = CPR(I40)
         CPPR(J40)   = CPPR(I40)
         QPSI(J40)   = QPSI(I40)
         CDQ(J40)    = CDQ(I40)
         CP(J40)     = CP(I40)
         CPDP(J40)   = CPDP(I40)
         RIPR(J40)   = RIPR(I40)
         RJDOTB(J40) = RJDOTB(I40)
         CIPR(J40)   = CIPR(I40)
         CID0(J40)   = CID0(I40)
         CID2(J40)   = CID2(I40)
         CIDR(J40)   = CIDR(I40)
         CIDQ(J40)   = CIDQ(I40)
         RIP(J40)    = RIP(I40)
         RIP2(J40)   = RIP2(I40)
         RIB2(J40)   = RIB2(I40)
         RIVOL(J40)  = RIVOL(I40)
         RARE(J40)   = RARE(I40)
         RIIE(J40)   = RIIE(I40)
         RIIR(J40)   = RIIR(I40)
         RLENG(J40)  = RLENG(I40)
         RLENG1(J40) = RLENG1(I40)
         RFCIRC(J40) = RFCIRC(I40)
         RJBSR(J40)  = RJBSR(I40)
         RJBSH(J40)  = RJBSH(I40)
         RJBSOS(J40,1) = RJBSOS(I40,1)
         RJBSOS(J40,2) = RJBSOS(I40,2)
         ARATIO(J40) = ARATIO(I40)
         RJPAR(J40)  = RJPAR(I40)
         RB2AV(J40)  = RB2AV(I40)
         RB2MAX(J40) = RB2MAX(I40)
         RELL(J40)   = RELL(I40)
         RJ1(J40)    = RJ1(I40)
         RJ2(J40)    = RJ2(I40)
         RJ3(J40)    = RJ3(I40)
         RJ4(J40)    = RJ4(I40)
         RJ5(J40)    = RJ5(I40)
         RJ5P(J40)   = RJ5P(I40)
         RJ6(J40)    = RJ6(I40)
C
         IF (NPROFZ .EQ. 1) THEN
C
            DENSTY(J40)   = DENSTY(I40)
            TEMPER(J40)   = TEMPER(I40)
C
         ENDIF
C
         DO 38 J38=1,NMGAUS*NT1
C
         RRISO(J38,J40)  = RRISO(J38,I40)
         RZISO(J38,J40)  = RZISO(J38,I40)
         GPISO(J38,J40)  = GPISO(J38,I40)
         RHOISO(J38,J40) = RHOISO(J38,I40)
         BNDISO(J38,J40) = BNDISO(J38,I40)
         DPSISO(J38,J40) = DPSISO(J38,I40)
         DPTISO(J38,J40) = DPTISO(J38,I40)
         DGNISO(J38,J40) = DGNISO(J38,I40)
         DRNISO(J38,J40) = DRNISO(J38,I40)
         SIGPSI(J38,J40) = SIGPSI(J38,I40)
         TETPSI(J38,J40) = TETPSI(J38,I40)
         WGTPSI(J38,J40) = WGTPSI(J38,I40)
C
 38      CONTINUE
C
         DO 39 J39=1,NT2
C
         SIGMAP(J39,J40) = SIGMAP(J39,I40)
         TETMAP(J39,J40) = TETMAP(J39,I40)
         CHIO(J39,J40)   = CHIO(J39,I40)
         CHIN(J39,J40)   = CHIN(J39,I40)
         BCHIO(J39,J40)  = BCHIO(J39,I40)
         BCHIN(J39,J40)  = BCHIN(J39,I40)
C
 39      CONTINUE
   40    CONTINUE
C
         CALL GLOQUA(CSM,CS,NPSI1,1)
C     
         RETURN
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
   50    CONTINUE
C
         IN  = NPSI * (NMGAUS + 2)
         I1  = 1
         I2  = 1
         I51 = 0
C
         DO 51 J51=1,IN
C
         IF (MOD(J51,NMGAUS+2) .EQ. 3) THEN
            ZS(J51) = CSM(I1)
            I1      = I1 + 1
         ELSE
            ZS(J51) = CSPEN(I2)
            I2      = I2 + 1
         ENDIF
C
         CALL SURFACE(J51,SIGPSI(1,J51),TETPSI(1,J51),WGTPSI(1,J51),
     ,                ZS(J51))
C
         IF (MOD(J51,6) .NE. 3) THEN
C
            I51 = I51 + 1
C
            CALL TPSI(IN,J51)
            CALL OUTPEN(I51)
C
            ZTMF(I51)  = TMF(J51)
            ZTPR(I51)  = TTP(J51) / TMF(J51)
            ZCPR(I51)  = CPR(J51)
            ZCPPR(I51) = CPPR(J51)
            ZQ(I51)    = QPSI(J51)
            ZQP(I51)   = CDQ(J51)
C
         ENDIF
C
   51    CONTINUE
C
         IF (NSURF .NE. 1) THEN
C
            CP0   = FCCCC0(CPR(1),CPR(2),CPR(3),CPR(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
            DPDP0 = FCCCC0(CPPR(1),CPPR(2),CPPR(3),CPPR(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
            T0    = FCCCC0(TMF(1),TMF(2),TMF(3),TMF(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
            DTTP0 = FCCCC0(TTP(1),TTP(2),TTP(3),TTP(4),
     ,                     ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
C
         ELSE
C
            T0    = 1.
            DTTP0 = 0.
            CP0   = CPP * SPSIM
            DPDP0 = CPP
C
         ENDIF
C
         Q0    = FCCCC0(QPSI(1),QPSI(2),QPSI(3),QPSI(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         RIPR0 = FCCCC0(RIPR(1),RIPR(2),RIPR(3),RIPR(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         RJDTB0= FCCCC0(RJDOTB(1),RJDOTB(2),RJDOTB(3),RJDOTB(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         CPSI0 = FCCCC0(CP(1),CP(2),CP(3),CP(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         DQDP0 = FCCCC0(CDQ(1),CDQ(2),CDQ(3),CDQ(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         CPDP0 = FCCCC0(CPDP(1),CPDP(2),CPDP(3),CPDP(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
C
         IF (NPROFZ .EQ. 1) THEN
C
         TEMP0 = FCCCC0(TEMPER(1),TEMPER(2),TEMPER(3),TEMPER(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
         DENS0 = FCCCC0(DENSTY(1),DENSTY(2),DENSTY(3),DENSTY(4),
     ,                  ZS(1),ZS(2),ZS(3),ZS(4),RC0P)
C
         ENDIF
C
         NPESUR = 51
C
         OPEN(UNIT=NPESUR,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     O        FILE='NPESUR')
         REWIND NPESUR
C
         IDCHSE = 222
C
         WRITE(NPESUR) IDCHSE
         WRITE(NPESUR) NPSI, NCHI, NMGAUS
         WRITE(NPESUR) RMAG, RZMAG, CPSRF, T0, Q0
         WRITE(NPESUR) (CSPEN(L),L=1,NPSI*(NMGAUS+1))
         WRITE(NPESUR) (CTPEN(L),L=1,NCHI*(NMGAUS+1))
         WRITE(NPESUR) (ZTMF(L),L=1,NPSI*(NMGAUS+1))
         WRITE(NPESUR) (ZTPR(L),L=1,NPSI*(NMGAUS+1))
         WRITE(NPESUR) (ZCPR(L),L=1,NPSI*(NMGAUS+1))
         WRITE(NPESUR) (ZCPPR(L),L=1,NPSI*(NMGAUS+1))
         WRITE(NPESUR) (ZQ(L),L=1,NPSI*(NMGAUS+1))
         WRITE(NPESUR) (ZQP(L),L=1,NPSI*(NMGAUS+1))
         CLOSE(UNIT=NPESUR,STATUS='KEEP')
C
         DO 56 J56=1,NPSI1
C
         IF (J56 .EQ. NPSI1) THEN
            I56 = IN
         ELSE
            I56 = (J56-1)*(NMGAUS+2)+3
         ENDIF
C
         PSIISO(J56) = PSIISO(I56)
         TMF(J56)    = TMF(I56)
         TTP(J56)    = TTP(I56)
         CPR(J56)    = CPR(I56)
         CPPR(J56)   = CPPR(I56)
         QPSI(J56)   = QPSI(I56)
         CDQ(J56)    = CDQ(I56)
         CP(J56)     = CP(I56)
         CPDP(J56)   = CPDP(I56)
         RIPR(J56)   = RIPR(I56)
         RJDOTB(J56) = RJDOTB(I56)
         CIPR(J56)   = CIPR(I56)
         CID0(J56)   = CID0(I56)
         CID2(J56)   = CID2(I56)
         CIDR(J56)   = CIDR(I56)
         CIDQ(J56)   = CIDQ(I56)
         RIP(J56)    = RIP(I56)
         RIP2(J56)   = RIP2(I56)
         RIB2(J56)   = RIB2(I56)
         RIVOL(J56)  = RIVOL(I56)
         RARE(J56)   = RARE(I56)
         RIIE(J56)   = RIIE(I56)
         RIIR(J56)   = RIIR(I56)
         RLENG(J56)  = RLENG(I56)
         RLENG1(J56) = RLENG1(I56)
         RFCIRC(J56) = RFCIRC(I56)
         RJBSR(J56)  = RJBSR(I56)
         RJBSH(J56)  = RJBSH(I56)
         RJBSOS(J56,1) = RJBSOS(I56,1)
         RJBSOS(J56,2) = RJBSOS(I56,2)
         ARATIO(J56) = ARATIO(I56)
         RJPAR(J56)  = RJPAR(I56)
         RB2AV(J56)  = RB2AV(I56)
         RB2MAX(J56) = RB2MAX(I56)
         RELL(J56)   = RELL(I56)
         RJ1(J56)    = RJ1(I56)
         RJ2(J56)    = RJ2(I56)
         RJ3(J56)    = RJ3(I56)
         RJ4(J56)    = RJ4(I56)
         RJ5(J56)    = RJ5(I56)
         RJ5P(J56)   = RJ5P(I56)
         RJ6(J56)    = RJ6(I56)
C
         IF (NPROFZ .EQ. 1) THEN
C
            DENSTY(J56)   = DENSTY(I56)
            TEMPER(J56)   = TEMPER(I56)
C
         ENDIF
C
         DO 53 J53=1,NMGAUS*NT1
C
         RRISO(J53,J56)  = RRISO(J53,I56)
         RZISO(J53,J56)  = RZISO(J53,I56)
         GPISO(J53,J56)  = GPISO(J53,I56)
         RHOISO(J53,J56) = RHOISO(J53,I56)
         BNDISO(J53,J56) = BNDISO(J53,I56)
         DPSISO(J53,J56) = DPSISO(J53,I56)
         DPTISO(J53,J56) = DPTISO(J53,I56)
         DGNISO(J53,J56) = DGNISO(J53,I56)
         DRNISO(J53,J56) = DRNISO(J53,I56)
         SIGPSI(J53,J56) = SIGPSI(J53,I56)
         TETPSI(J53,J56) = TETPSI(J53,I56)
         WGTPSI(J53,J56) = WGTPSI(J53,I56)
C
   53    CONTINUE
C
         DO 54 J54=1,NT2
C
         SIGMAP(J54,J56) = SIGMAP(J54,I56)
         TETMAP(J54,J56) = TETMAP(J54,I56)
         CHIO(J54,J56)   = CHIO(J54,I56)
         CHIN(J54,J56)   = CHIN(J54,I56)
         BCHIO(J54,J56)  = BCHIO(J54,I56)
         BCHIN(J54,J56)  = BCHIN(J54,I56)
C
   54    CONTINUE
   56    CONTINUE
C
         CALL GLOQUA(CSM,CS,NPSI1,1)
C
         RETURN
         END
C*DECK C2SM02
C*CALL PROCESS
         SUBROUTINE SURFACE(K,PSIGMA,PTETA,PGWGT,PS)
C        ###########################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM02 EVALUATE FLUX SURFACE INTEGRALS ALONG CONSTANT POLOIDAL FLUX *
*        SURFACES                                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         PARAMETER (NPGMAP=NPMGS*NTP1)
C
         DIMENSION
     I   IS0(NPGMAP),       IT0(NPGMAP),       IC(NPGMAP),
     R   PGWGT(*),          PSIGMA(*),         PTETA(*),
     R   ZBND(NPGMAP,5),    ZB2(NPGMAP),       ZDBDS(NPGMAP,16),
     R   ZDBDT(NPGMAP,16),  ZDBDST(NPGMAP,16), ZDCHIN(NPGMAP),
     R   ZDCHIO(NPGMAP),    ZDBCHIN(NPGMAP),   ZDBCHIO(NPGMAP),
     R   ZD2BS2(NPGMAP,16), ZD2BT2(NPGMAP,16), ZPCEL(NPGMAP,16),
     R   ZS1(NPGMAP),       ZS2(NPGMAP),
     R   ZTETA(NPGMAP,5),   ZT1(NPGMAP),
     R   ZT2(0:NPGMAP+1),   ZIVPAR(2*NPISO),   ZB2SHF(0:NPGMAP+1)
C
C     ARRAYS FOR MIN, MAX SEARCH (ADD PERIODIC POINT)
         DIMENSION
     +     ZRISO(0:NPGMAP+1), ZZISO(0:NPGMAP+1), ZGPISO(0:NPGMAP+1)
C
         INCLUDE 'QUAQQQ.inc'
C
         ABSQRT(XX) = SQRT(ABS(XX))
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ZEPS = 1.E-3
C
         ZRIPR1    = 0.
         ZRIPR2    = 0.
         RARE(K)   = 0.
         RIP(K)    = 0.
         RIP2(K)   = 0.
         RIB2(K)   = 0.
         RIVOL(K)  = 0.
         RIIE(K)   = 0.
         RIIR(K)   = 0.
         RLENG(K)  = 0.
         RLENG1(K) = 0.
         ZC2       = 0.
         RJ1(K)    = 0.
         RJ2(K)    = 0.
         RJ3(K)    = 0.
         RJ4(K)    = 0.
         RJ5(K)    = 0.
         RJ5P(K)   = 0.
         RJ6(K)    = 0.
         ZCDQ      = 0.
cab
cab      abs(dpsi/ds)
cab
         ZDPSIS = 2.*PS*CPSRF
cab
C
         IGMAX = NMGAUS * NT1
CC
         CALL RESETI(IC,IGMAX,1)
         DO 1 JT = 1,NT1
            DO 1 JG=1,IGMAX
               IF (IC(JG).EQ.0) GOTO 1
               IT0(JG) = JT-1
               IF (PTETA(JG).LE.CT(JT)) IC(JG)  = 0
 1          CONTINUE
         CALL RESETI(IC,IGMAX,1)
         DO 3 JS = 1,NS1
            DO 3 JG=1,IGMAX
               IF (IC(JG).EQ.0) GOTO 3
               IS0(JG) = JS-1
               IF (PSIGMA(JG).LE.CSIG(JS)) IC(JG)  = 0
 3       CONTINUE
C
         DO 5 J5=1,IGMAX
C
         IF (IS0(J5).LT.  1) IS0(J5) = 1
         IF (IS0(J5).GE.NS1) IS0(J5) = NS
         IF (IT0(J5).LT.  1) IT0(J5) = 1
         IF (IT0(J5).GE.NT1) IT0(J5) = NT
C
         ZTETA(J5,1) = PTETA(J5)
         ZTETA(J5,2) = PTETA(J5) - 2. * ZEPS
         ZTETA(J5,3) = PTETA(J5) -      ZEPS
         ZTETA(J5,4) = PTETA(J5) +      ZEPS
         ZTETA(J5,5) = PTETA(J5) + 2. * ZEPS
C
         ZS1(J5) = CSIG(IS0(J5))
         ZS2(J5) = CSIG(IS0(J5)+1)
         ZT1(J5) = CT(IT0(J5))
 5       ZT2(J5) = CT(IT0(J5)+1)
C
         CALL BOUND(IGMAX,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(IGMAX,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(IGMAX,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(IGMAX,ZTETA(1,4),ZBND(1,4))
         CALL BOUND(IGMAX,ZTETA(1,5),ZBND(1,5))
C
         CALL PSICEL(IS0,IT0,IGMAX,NPGMAP,ZPCEL,CPSICL)
         CALL BASIS3(IGMAX,NPGMAP,ZS1,ZS2,ZT1,ZT2(1),PSIGMA,PTETA,
     &               ZDBDS,ZDBDT,ZDBDST,ZD2BS2,ZD2BT2)
C
         DO 6 J6=1,IGMAX
C
         ZDRSDT = (ZBND(J6,2) + 8*(ZBND(J6,4) - ZBND(J6,3)) -
     -             ZBND(J6,5)) / (12. * ZEPS)
         ZD2RST = (- ZBND(J6,2) + 16. * ZBND(J6,3) -
     -             30. * ZBND(J6,1) + 16. * ZBND(J6,4) -
     -             ZBND(J6,5)) / (12. * ZEPS**2)
C
         ZDPDS = ZDBDS(J6, 1) * ZPCEL(J6, 1) +
     +           ZDBDS(J6, 2) * ZPCEL(J6, 2) +
     +           ZDBDS(J6, 3) * ZPCEL(J6, 3) +
     +           ZDBDS(J6, 4) * ZPCEL(J6, 4) +
     +           ZDBDS(J6, 5) * ZPCEL(J6, 5) +
     +           ZDBDS(J6, 6) * ZPCEL(J6, 6) +
     +           ZDBDS(J6, 7) * ZPCEL(J6, 7) +
     +           ZDBDS(J6, 8) * ZPCEL(J6, 8) +
     +           ZDBDS(J6, 9) * ZPCEL(J6, 9) +
     +           ZDBDS(J6,10) * ZPCEL(J6,10) +
     +           ZDBDS(J6,11) * ZPCEL(J6,11) +
     +           ZDBDS(J6,12) * ZPCEL(J6,12) +
     +           ZDBDS(J6,13) * ZPCEL(J6,13) +
     +           ZDBDS(J6,14) * ZPCEL(J6,14) +
     +           ZDBDS(J6,15) * ZPCEL(J6,15) +
     +           ZDBDS(J6,16) * ZPCEL(J6,16)
C
         ZDPDT = ZDBDT(J6, 1) * ZPCEL(J6, 1) +
     +           ZDBDT(J6, 2) * ZPCEL(J6, 2) +
     +           ZDBDT(J6, 3) * ZPCEL(J6, 3) +
     +           ZDBDT(J6, 4) * ZPCEL(J6, 4) +
     +           ZDBDT(J6, 5) * ZPCEL(J6, 5) +
     +           ZDBDT(J6, 6) * ZPCEL(J6, 6) +
     +           ZDBDT(J6, 7) * ZPCEL(J6, 7) +
     +           ZDBDT(J6, 8) * ZPCEL(J6, 8) +
     +           ZDBDT(J6, 9) * ZPCEL(J6, 9) +
     +           ZDBDT(J6,10) * ZPCEL(J6,10) +
     +           ZDBDT(J6,11) * ZPCEL(J6,11) +
     +           ZDBDT(J6,12) * ZPCEL(J6,12) +
     +           ZDBDT(J6,13) * ZPCEL(J6,13) +
     +           ZDBDT(J6,14) * ZPCEL(J6,14) +
     +           ZDBDT(J6,15) * ZPCEL(J6,15) +
     +           ZDBDT(J6,16) * ZPCEL(J6,16)
C
         ZD2PST = ZDBDST(J6, 1) * ZPCEL(J6, 1) +
     +            ZDBDST(J6, 2) * ZPCEL(J6, 2) +
     +            ZDBDST(J6, 3) * ZPCEL(J6, 3) +
     +            ZDBDST(J6, 4) * ZPCEL(J6, 4) +
     +            ZDBDST(J6, 5) * ZPCEL(J6, 5) +
     +            ZDBDST(J6, 6) * ZPCEL(J6, 6) +
     +            ZDBDST(J6, 7) * ZPCEL(J6, 7) +
     +            ZDBDST(J6, 8) * ZPCEL(J6, 8) +
     +            ZDBDST(J6, 9) * ZPCEL(J6, 9) +
     +            ZDBDST(J6,10) * ZPCEL(J6,10) +
     +            ZDBDST(J6,11) * ZPCEL(J6,11) +
     +            ZDBDST(J6,12) * ZPCEL(J6,12) +
     +            ZDBDST(J6,13) * ZPCEL(J6,13) +
     +            ZDBDST(J6,14) * ZPCEL(J6,14) +
     +            ZDBDST(J6,15) * ZPCEL(J6,15) +
     +            ZDBDST(J6,16) * ZPCEL(J6,16)
C
         ZD2PS2 = ZD2BS2(J6, 1) * ZPCEL(J6, 1) +
     +            ZD2BS2(J6, 2) * ZPCEL(J6, 2) +
     +            ZD2BS2(J6, 3) * ZPCEL(J6, 3) +
     +            ZD2BS2(J6, 4) * ZPCEL(J6, 4) +
     +            ZD2BS2(J6, 5) * ZPCEL(J6, 5) +
     +            ZD2BS2(J6, 6) * ZPCEL(J6, 6) +
     +            ZD2BS2(J6, 7) * ZPCEL(J6, 7) +
     +            ZD2BS2(J6, 8) * ZPCEL(J6, 8) +
     +            ZD2BS2(J6, 9) * ZPCEL(J6, 9) +
     +            ZD2BS2(J6,10) * ZPCEL(J6,10) +
     +            ZD2BS2(J6,11) * ZPCEL(J6,11) +
     +            ZD2BS2(J6,12) * ZPCEL(J6,12) +
     +            ZD2BS2(J6,13) * ZPCEL(J6,13) +
     +            ZD2BS2(J6,14) * ZPCEL(J6,14) +
     +            ZD2BS2(J6,15) * ZPCEL(J6,15) +
     +            ZD2BS2(J6,16) * ZPCEL(J6,16)
C
         ZD2PT2 = ZD2BT2(J6, 1) * ZPCEL(J6, 1) +
     +            ZD2BT2(J6, 2) * ZPCEL(J6, 2) +
     +            ZD2BT2(J6, 3) * ZPCEL(J6, 3) +
     +            ZD2BT2(J6, 4) * ZPCEL(J6, 4) +
     +            ZD2BT2(J6, 5) * ZPCEL(J6, 5) +
     +            ZD2BT2(J6, 6) * ZPCEL(J6, 6) +
     +            ZD2BT2(J6, 7) * ZPCEL(J6, 7) +
     +            ZD2BT2(J6, 8) * ZPCEL(J6, 8) +
     +            ZD2BT2(J6, 9) * ZPCEL(J6, 9) +
     +            ZD2BT2(J6,10) * ZPCEL(J6,10) +
     +            ZD2BT2(J6,11) * ZPCEL(J6,11) +
     +            ZD2BT2(J6,12) * ZPCEL(J6,12) +
     +            ZD2BT2(J6,13) * ZPCEL(J6,13) +
     +            ZD2BT2(J6,14) * ZPCEL(J6,14) +
     +            ZD2BT2(J6,15) * ZPCEL(J6,15) +
     +            ZD2BT2(J6,16) * ZPCEL(J6,16)
C
         ZFP    = (ZDPDS**2 + (ZDPDT / PSIGMA(J6) - ZDPDS * ZDRSDT /
     /            ZBND(J6,1))**2) / ZBND(J6,1)**2
         ZGRADP = SQRT(ZFP)
C
         ZCOST  = COS(PTETA(J6))
         ZSINT  = SIN(PTETA(J6))
C
         ZRHO    = PSIGMA(J6) * ZBND(J6,1)
         ZR      = ZRHO * ZCOST + R0
         ZZ      = ZRHO * ZSINT + RZ0
         ZJPHI   = - ZR * CPPR(K) - TTP(K) / ZR
         ZB2(J6) = (TMF(K)**2 + ZFP) / ZR**2
         ZR2     = ZR**2
C
         ZJAC1 = ZR**(NER-1) * ZGRADP**NEGP * ZDPDS
         ZJAC2 = ZR * ZDPDS
C
         ZINT1 = ZRHO * ZBND(J6,1) / ZJAC1
         ZINT2 = ZRHO * ZBND(J6,1) / ZJAC2
C
         ZDCHIN(J6) = ZINT1 * PGWGT(J6)
         ZDCHIO(J6) = ZINT2 * PGWGT(J6)
C
         ZIPR1 = ZJPHI * ZRHO * ZBND(J6,1) * ZR / ZJAC2
         ZIPR2 = ZRHO * ZBND(J6,1) * ZR / ZJAC2
C
         ZRIPR1 = ZRIPR1 + ZIPR1 * PGWGT(J6)
         ZRIPR2 = ZRIPR2 + ZIPR2 * PGWGT(J6)
C
         ZIDA    = ZRHO * ZBND(J6,1) * ZR     / ZJAC2
         ZINTP   = ZRHO * ZBND(J6,1) * ZR2    / ZJAC2 * CPR(K)
         ZINTP2  = ZRHO * ZBND(J6,1) * ZR2    / ZJAC2 * CPR(K)**2
         ZINTB2  = ZRHO * ZBND(J6,1) * ZR2    / ZJAC2 * ZB2(J6)
         ZINTV   = ZRHO * ZBND(J6,1) * ZR2    / ZJAC2
         ZINTIE  = ZRHO * ZBND(J6,1) * ZR     / ZJAC2 * ZJPHI
         ZINTIR  = ZRHO * ZBND(J6,1) * ZFP    / ZJAC2
         ZINTLE  = ZRHO * ZBND(J6,1) * ZR     / ZJAC2 * ZGRADP
         ZINTLE1 = ZRHO * ZBND(J6,1)          / ZJAC2 * ZGRADP
C
         RARE(K)   = RARE(K)   + ZIDA    * PGWGT(J6)
         RIP(K)    = RIP(K)    + ZINTP   * PGWGT(J6)
         RIP2(K)   = RIP2(K)   + ZINTP2  * PGWGT(J6)
         RIB2(K)   = RIB2(K)   + ZINTB2  * PGWGT(J6)
         RIVOL(K)  = RIVOL(K)  + ZINTV   * PGWGT(J6)
         RIIE(K)   = RIIE(K)   + ZINTIE  * PGWGT(J6)
         RIIR(K)   = RIIR(K)   + ZINTIR  * PGWGT(J6)
         RLENG(K)  = RLENG(K)  + ZINTLE  * PGWGT(J6)
         RLENG1(K) = RLENG1(K) + ZINTLE1 * PGWGT(J6)
         ZC2       = ZC2       + ZINT2   * PGWGT(J6)
C
         ZDSDR = (ZDRSDT * ZSINT + ZBND(J6,1) * ZCOST) / ZBND(J6,1)**2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBND(J6,1) * ZSINT - ZDRSDT * ZCOST) / ZBND(J6,1)**2
         ZDTDZ = ZCOST / ZRHO
C
         ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         Z1     = ZDPDS
         ZDZ1DS = ZD2PS2
         ZDZ1DT = ZD2PST
         Z2     = ZDPDT / PSIGMA(J6) - ZDRSDT * ZDPDS / ZBND(J6,1)
         ZDZ2DS = (ZD2PST - ZDPDT / PSIGMA(J6))  / PSIGMA(J6) -
     -            ZDRSDT * ZD2PS2 / ZBND(J6,1)
         ZDZ2DT = - ZD2PST * ZDRSDT / ZBND(J6,1) + ZD2PT2 / PSIGMA(J6) +
     +            ZDPDS * (ZDRSDT**2-ZD2RST*ZBND(J6,1)) / ZBND(J6,1)**2
C
         ZDFDS = 2 * (Z1 * ZDZ1DS + Z2 * ZDZ2DS) / ZBND(J6,1)**2
         ZDFDT = - 2 * ZFP * ZDRSDT / ZBND(J6,1) +
     +           2 * (Z1 * ZDZ1DT + Z2 * ZDZ2DT) / ZBND(J6,1)**2
C
         ZDFDR = ZDFDS * ZDSDR + ZDFDT * ZDTDR
         ZDFDZ = ZDFDS * ZDSDZ + ZDFDT * ZDTDZ
C
         ZDGDPN = (ZDFDR * ZDPDR + ZDFDZ * ZDPDZ) / ZFP
         ZDRDPN = ZDPDR / ZFP
C
         ZBINT1 = ZRHO * ZBND(J6,1) * ((ZR * ZJPHI -
     /            (1. + .5 * NEGP) * ZDGDPN) / ZFP -
     -            (NER - 2) * ZDRDPN / ZR) / ZJAC1
         ZBINT2 = ZRHO * ZBND(J6,1) * (ZR * ZJPHI - ZDGDPN) /
     /            (ZFP * ZJAC2)
C
         ZDBCHIN(J6) = ZBINT1 * PGWGT(J6)
         ZDBCHIO(J6) = ZBINT2 * PGWGT(J6)
C
C     COMPUTE INTEGRALS FOR LOCAL INTECHANGE FORMULAS (EQ.(22) IN
C     CHEASE PAPER)
C
         ZIJ1  = ZRHO * ZBND(J6,1)          / ZJAC2 / ZFP
         ZIJ2  = ZRHO * ZBND(J6,1) * ZR2    / ZJAC2 / ZFP
         ZIJ3  = ZRHO * ZBND(J6,1) * ZR2**2 / ZJAC2 / ZFP
         ZIJ4  = ZRHO * ZBND(J6,1)          / ZJAC2
         ZIJ5  = ZRHO * ZBND(J6,1) * ZR2    / ZJAC2
         ZIJ6  = ZRHO * ZBND(J6,1)          / ZJAC2 * ZFP
         ZIJ5P = (ZBND(J6,1) / ZDPDS)**2 *
     &           (2. * ZR - R0 - PSIGMA(J6) * ZR * ZD2PS2 / ZDPDS)
         ZIDQ  = (ZBND(J6,1) / ZDPDS)**2 / ZR *
     &           (R0 / ZR - PSIGMA(J6) * ZD2PS2 / ZDPDS)
C
         RJ1(K)  = RJ1(K)  + ZIJ1  * PGWGT(J6)
         RJ2(K)  = RJ2(K)  + ZIJ2  * PGWGT(J6)
         RJ3(K)  = RJ3(K)  + ZIJ3  * PGWGT(J6)
         RJ4(K)  = RJ4(K)  + ZIJ4  * PGWGT(J6)
         RJ5(K)  = RJ5(K)  + ZIJ5  * PGWGT(J6)
         RJ5P(K) = RJ5P(K) + ZIJ5P * PGWGT(J6)
         RJ6(K)  = RJ6(K)  + ZIJ6  * PGWGT(J6)
         ZCDQ    = ZCDQ    + ZIDQ  * PGWGT(J6)
C
         RRISO(J6,K)  = ZR
         RZISO(J6,K)  = ZZ
         GPISO(J6,K)  = ZGRADP
         RHOISO(J6,K) = ZRHO
         BNDISO(J6,K) = ZBND(J6,1)
         DPTISO(J6,K) = ZDPDT
         DPSISO(J6,K) = ZDPDS
         DGNISO(J6,K) = ZDGDPN
         DGRISO(J6,K) = ZDFDR
         DGZISO(J6,K) = ZDFDZ
         DRNISO(J6,K) = ZDRDPN
         DPRISO(J6,K) = ZDPDR
         DPZISO(J6,K) = ZDPDZ
C
 6       CONTINUE
C
         BCHIN(1,K) = 0.
         BCHIO(1,K) = 0.
         CHIN(1,K)  = 0.
         CHIO(1,K)  = 0.
C
         DO 7 J7=1,NT1
C
         I7 = (J7 - 1) * NMGAUS
C
         BCHIN(J7+1,K) = BCHIN(J7,K) + ZDBCHIN(I7+1) + ZDBCHIN(I7+2) +
     &                                 ZDBCHIN(I7+3) + ZDBCHIN(I7+4)
         BCHIO(J7+1,K) = BCHIO(J7,K) + ZDBCHIO(I7+1) + ZDBCHIO(I7+2) +
     &                                 ZDBCHIO(I7+3) + ZDBCHIO(I7+4)
         CHIN(J7+1,K)  = CHIN(J7,K)  + ZDCHIN(I7+1)  + ZDCHIN(I7+2) +
     &                                 ZDCHIN(I7+3)  + ZDCHIN(I7+4)
         CHIO(J7+1,K)  = CHIO(J7,K)  + ZDCHIO(I7+1)  + ZDCHIO(I7+2) +
     &                                 ZDCHIO(I7+3)  + ZDCHIO(I7+4)
C
 7       CONTINUE
C
         ZCP1    = .5 * CHIO(NT2,K) / CPI
         QPSI(K) = ZCP1 * TMF(K)
         CDQ(K)  = .5 * TMF(K) * BCHIO(NT2,K) / CPI +
     +             ZCP1 * TTP(K) / TMF(K)
         CP(K)   = .5 * CHIN(NT2,K) / CPI
         CPDP(K) = .5 * BCHIN(NT2,K) / CPI
C
         CALL SSCAL(NT2,RC1P/ZCP1,CHIO(1,K),1)
         CALL SSCAL(NT2,RC1P/ZCP1,BCHIO(1,K),1)
         CALL SSCAL(NT2,RC1P/CP(K),CHIN(1,K),1)
         CALL SSCAL(NT2,RC1P/CP(K),BCHIN(1,K),1)
C
         CALL SAXPY(NT2,TTP(K)/(TMF(K)*TMF(K))-CDQ(K)/QPSI(K),
     ,              CHIO(1,K),1,BCHIO(1,K),1)
         CALL SAXPY(NT2,-CPDP(K)/CP(K),CHIN(1,K),1,
     ,              BCHIN(1,K),1)
C
         CALL SSCAL(NT2,ZDPSIS,BCHIO(1,K),1)
         CALL SSCAL(NT2,ZDPSIS,BCHIN(1,K),1)
C
         RIPR(K)   = ZRIPR1 / ZRIPR2
         RJDOTB(K) = - CPPR(K) * RIVOL(K) / ZC2 - TTP(K) *
     &                 (1. + RIIR(K) / (TMF(K)**2 * ZC2))
         RB2AV(K)  = RIB2(K) / RIVOL(K)
         RJPAR(K)  = - (TMF(K) * CPPR(K) + TTP(K) * RB2AV(K) / TMF(K))
C
C  NORMALIZE INTEGRALS FOR LOCAL INTERCHANGE STABILITY CRITERIA
C
         RJ1(K)  = RJ1(K) / (2. * CPI)
         RJ2(K)  = RJ2(K) / (2. * CPI)
         RJ3(K)  = RJ3(K) / (2. * CPI)
         RJ4(K)  = RJ4(K) / (2. * CPI)
         RJ5(K)  = RJ5(K) / (2. * CPI)
         RJ5P(K) = RJ5P(K) / (2. * CPI)
         RJ6(K)  = RJ6(K) / (2. * CPI)
C
         ZCDQ = ZCDQ / (2. * CPI)
         CDQ(K) = TMF(K) * ZCDQ + TTP(K) * RJ4(K) / TMF(K)
C
C     FOR FITS OF MIN AND MAX ALONG PTETA, ADD AN EXTRA POINT
C     BEFORE AND AFTER THE JUMP AND SHIFT PTETA TO ZT2 SUCH
C     THAT ZT2(1) JUST AFTER JUMP => ADD ZT2(0) AND ZT2(IGMAX+1)
C
C     FIND JTJUMP
C
         DO JT=1,IGMAX-1
           IF (PTETA(JT+1) .LE. PTETA(JT)) JTJUMP = JT+1
         ENDDO
C     COPY ARRAYS FOR INTERPOLATION WITH JTJUMP SHIFT
C
         DO JT=JTJUMP,IGMAX
           JTEFF = JT-JTJUMP+1
           ZT2(JTEFF)    = PTETA(JT)
           ZRISO(JTEFF)  = RRISO(JT,K)
           ZZISO(JTEFF)  = RZISO(JT,K)
           ZGPISO(JTEFF) = GPISO(JT,K)
           ZB2SHF(JTEFF) = ZB2(JT)
         END DO
         JT2LAST = IGMAX-JTJUMP+1
         DO JT=1,JTJUMP-1
           ZT2(JT2LAST+JT)    = PTETA(JT)
           ZRISO(JT2LAST+JT)  = RRISO(JT,K)
           ZZISO(JT2LAST+JT)  = RZISO(JT,K)
           ZGPISO(JT2LAST+JT) = GPISO(JT,K)
           ZB2SHF(JT2LAST+JT) = ZB2(JT)
         END DO
         ZT2(0) = ZT2(IGMAX) - 2.*CPI
         ZT2(IGMAX+1) = ZT2(1) + 2.*CPI
         ZRISO(0) = ZRISO(IGMAX)
         ZRISO(IGMAX+1) = ZRISO(1)
         ZZISO(0) = ZZISO(IGMAX)
         ZZISO(IGMAX+1) = ZZISO(1)
         ZB2SHF(0) = ZB2SHF(IGMAX)
         ZB2SHF(IGMAX+1) = ZB2SHF(1)
         ZGPISO(0) = ZGPISO(IGMAX)
         ZGPISO(IGMAX+1) = ZGPISO(1)
C
c%OS
         DO J=0,IGMAX
           IF (ZT2(J+1) .LE. ZT2(J)) THEN
             PRINT *,' ERROR IN ZT2 IN SURFACE:'
             PRINT *,' ZT2(J=',J,')= ',ZT2(J),'  >  ','ZT2(J+1)= '
     +         ,ZT2(J+1)
             PRINT *,' IGMAX= ',IGMAX,' K,PS= ',K,PS
             STOP 'ZT2'
           ENDIF
         ENDDO
c%OS
         IRMAX = ISMAX(IGMAX,ZRISO(1),1)
         IRMIN = ISMIN(IGMAX,ZRISO(1),1)
         IZMAX = ISMAX(IGMAX,ZZISO(1),1)
         IZMIN = ISMIN(IGMAX,ZZISO(1),1)
C
         ZTMAX = - .5 * FB1(ZRISO(IRMAX-1),ZRISO(IRMAX),
     ,                      ZRISO(IRMAX+1),ZT2(IRMAX-1),
     ,                      ZT2(IRMAX),ZT2(IRMAX+1)) /
     /                  FB2(ZRISO(IRMAX-1),ZRISO(IRMAX),
     ,                      ZRISO(IRMAX+1),ZT2(IRMAX-1),
     ,                      ZT2(IRMAX),ZT2(IRMAX+1))
C
         ZRMAX = FQQQ0(ZRISO(IRMAX-1),ZRISO(IRMAX),
     ,                 ZRISO(IRMAX+1),ZT2(IRMAX-1),
     ,                 ZT2(IRMAX),ZT2(IRMAX+1),ZTMAX)
C
         ZTMIN = - .5 * FB1(ZRISO(IRMIN-1),ZRISO(IRMIN),
     ,                      ZRISO(IRMIN+1),ZT2(IRMIN-1),
     ,                      ZT2(IRMIN),ZT2(IRMIN+1)) /
     /                  FB2(ZRISO(IRMIN-1),ZRISO(IRMIN),
     ,                      ZRISO(IRMIN+1),ZT2(IRMIN-1),
     ,                      ZT2(IRMIN),ZT2(IRMIN+1))
C
         ZRMIN = FQQQ0(ZRISO(IRMIN-1),ZRISO(IRMIN),
     ,                 ZRISO(IRMIN+1),ZT2(IRMIN-1),
     ,                 ZT2(IRMIN),ZT2(IRMIN+1),ZTMIN)
C
         ZTMAX = - .5 * FB1(ZZISO(IZMAX-1),ZZISO(IZMAX),
     ,                      ZZISO(IZMAX+1),ZT2(IZMAX-1),
     ,                      ZT2(IZMAX),ZT2(IZMAX+1)) /
     /                  FB2(ZZISO(IZMAX-1),ZZISO(IZMAX),
     ,                      ZZISO(IZMAX+1),ZT2(IZMAX-1),
     ,                      ZT2(IZMAX),ZT2(IZMAX+1))
C
         ZZMAX = FQQQ0(ZZISO(IZMAX-1),ZZISO(IZMAX),
     ,                 ZZISO(IZMAX+1),ZT2(IZMAX-1),
     ,                 ZT2(IZMAX),ZT2(IZMAX+1),ZTMAX)
C
         ZTMIN = - .5 * FB1(ZZISO(IZMIN-1),ZZISO(IZMIN),
     ,                      ZZISO(IZMIN+1),ZT2(IZMIN-1),
     ,                      ZT2(IZMIN),ZT2(IZMIN+1)) /
     /                  FB2(ZZISO(IZMIN-1),ZZISO(IZMIN),
     ,                      ZZISO(IZMIN+1),ZT2(IZMIN-1),
     ,                      ZT2(IZMIN),ZT2(IZMIN+1))
C
         ZZMIN = FQQQ0(ZZISO(IZMIN-1),ZZISO(IZMIN),
     ,                 ZZISO(IZMIN+1),ZT2(IZMIN-1),
     ,                 ZT2(IZMIN),ZT2(IZMIN+1),ZTMIN)
C
         ZX1 = ZZMAX - ZZMIN
         ZX2 = ZRMAX - ZRMIN
C
         RELL(K)   = (ZX1 - ZX2) / (ZX1 + ZX2)
         ARATIO(K) = (ZRMAX + ZRMIN) / (ZRMAX - ZRMIN)
C
         IRMAX = ISMAX(IGMAX,ZB2SHF(1),1)
         IRMIN = ISMIN(IGMAX,ZB2SHF(1),1)
C
         ZTMAX = - .5 * FB1(ZB2SHF(IRMAX-1),ZB2SHF(IRMAX),
     ,                      ZB2SHF(IRMAX+1),ZT2(IRMAX-1),
     ,                      ZT2(IRMAX),ZT2(IRMAX+1)) /
     /                  FB2(ZB2SHF(IRMAX-1),ZB2SHF(IRMAX),
     ,                      ZB2SHF(IRMAX+1),ZT2(IRMAX-1),
     ,                      ZT2(IRMAX),ZT2(IRMAX+1))
C
         RB2MAX(K) = FQQQ0(ZB2SHF(IRMAX-1),ZB2SHF(IRMAX),
     ,                     ZB2SHF(IRMAX+1),ZT2(IRMAX-1),
     ,                     ZT2(IRMAX),ZT2(IRMAX+1),ZTMAX)
C
         JB2MIN = IRMIN-1
         ZTB2MIN = - .5 * FB1(ZB2SHF(IRMIN-1),ZB2SHF(IRMIN),
     ,                      ZB2SHF(IRMIN+1),ZT2(IRMIN-1),
     ,                      ZT2(IRMIN),ZT2(IRMIN+1)) /
     /                  FB2(ZB2SHF(IRMIN-1),ZB2SHF(IRMIN),
     ,                      ZB2SHF(IRMIN+1),ZT2(IRMIN-1),
     ,                      ZT2(IRMIN),ZT2(IRMIN+1))
         ZB2MIN = FQQQ0(ZB2SHF(IRMIN-1),ZB2SHF(IRMIN),
     ,                     ZB2SHF(IRMIN+1),ZT2(IRMIN-1),
     ,                     ZT2(IRMIN),ZT2(IRMIN+1),ZTB2MIN)
C
C        FOR MIDPOINT INTEGRATION IN Y = LAMBDA * BMAX OF
C        FCIRC = FRACTION OF CIRCULATING PARTICLES IN HIRSHMAN,
C        PHYS.FLUIDS 31 (10), 1988, P3150, THE CSIPR(J),J=1,NISO-1
C        MESH IS USED FOR Y. BMAX IS ASSUMED TO BE ON THE INSIDE
C        OF THE PLASMA.
C
         DO 9 J9=1,NISO-1
C
         ZIVPAR(J9) = 0.
C
         DO 8 J8=1,IGMAX
C
         ZVPAR  = ZRISO(J8) * PSIGMA(J8) * ZBND(J8,1)**2 *
     &            ABSQRT(1. - CSIPR(J9) * SQRT(ZB2(J8) / RB2MAX(K)))/
     &            DPSISO(J8,K)
C
         ZIVPAR(J9) = ZIVPAR(J9) + ZVPAR * PGWGT(J8)
C
    8    CONTINUE
C
         ZIVPAR(J9) = ZIVPAR(J9) / RIVOL(K)
C
 9       CONTINUE
C
C INTEGRATE FCIRC USING CONSTANT STEPS
C
         RFCIRC(K) = 0.
C
         DO 10 J10=1,NISO-1
C
         RFCIRC(K) = RFCIRC(K) + 0.75 * RB2AV(K) * CSIPR(J10) *
     &               (CSIPRI(J10+1) - CSIPRI(J10)) /
     &               (RB2MAX(K) * ZIVPAR(J10))
C
   10    CONTINUE
C
C  BOOTSTRAP CURRENT DENSITY GIVEN IN M.N.ROSENBLUTH,
C  H.D.HAZELTINE, F.L.HINTON, PHYS. FLUIDS 15, 116, 1972,
C  STORED IN RJBSR(K).
C
         ZSQRTA   = 1 / SQRT(ARATIO(K))
         RJBSR(K) = .5 * ZSQRTA * (.29 * ZSQRTA - .9) * CPPR(K) * R0
C
C  BOOTSTRAP CURRENT DENSITY GIVEN IN HIRSHMAN, PHYS.FLUIDS 31 (10),
C  1988, P3150, STORED IN RJBSH(K). RZION IS THE ION CHARGE,
C  ETAEI = D(LOG(T)) / D(LOG(N)) AND ZX THE FRACTION FTRAPPED / FCIRC.
C
         ZJ0 = TMF(K)
         ZX  = (1. - RFCIRC(K)) / RFCIRC(K)
         ZX1 = (2.40 * ZSQRTA + 1.46) * ZSQRTA / (1 - ZSQRTA)**1.5
         ZX2 = 1.46 * ZSQRTA
C
C  FIX-UP TO GET RID OF INACCURATE BEHAVIOUR NEAR THE ORIGIN
C
         IF (ZX .GT. ZX1) ZX = ZX1
C
         ZDX  = RZION * (1.414 + RZION) + ZX * (0.754 + RZION * (2.657 +
     &          2. * RZION)) + ZX**2 * (0.348 + RZION * (1.243 + RZION))
         ZL31 = ZJ0 * ZX * (0.754 + RZION * (2.21 + RZION) +
     &          ZX * (0.348 + RZION * (1.243 + RZION))) / ZDX
         ZL32 = - ZJ0 * ZX * (0.884 + 2.074 * RZION) / ZDX
         ZA1 = CPPR(K)
C     DP_E/DPSI = RZION/(RZION+1) * CPPR(K)
         ZA2I = ETAEI / (1 + ETAEI) * RZION/(RZION+1) * ZA1
         ZA2E = ZA2I
C
         ZALFA    = - 1.172 / (1. + 0.462 * ZX)
         RJBSH(K) = - ZL31 * (ZA1 + ZALFA * ZA2I / RZION) - ZL32 * ZA2E
C
C  BOOTSTRAP CURRENT DENSITY GIVEN FROM O.SAUTER (INCLUDING COLLISIONALITY)
C  RJBSOS(K). RZION IS THE ION CHARGE, RPOPE=P/P_E, T_I/T_E=CST
C  ZFT THE FRACTION FTRAPPED, NEED ALSO T_E PROFILE => USE ETAEI OR AT4
C  IF AT4(1)=0:    T_E FROM ETAEI AND ASSUME N_E[1E19]=T_E[KEV]
C  IF AT4(1).NE.0: T_E GIVEN AS POLYNOMIALS OF S**2 WITH AT4 IN KEV
C
         ZFT  = 1. - RFCIRC(K)
         ZFT2 = 1.46 * ZSQRTA
C  FIX-UP TO GET RID OF INACCURATE BEHAVIOUR NEAR THE ORIGIN
         IF (ZFT .LE. 0.0) ZFT = ZFT2
C
C     TEMPERATURE OF ELECTRONS
C
         ZMU0 = 4.E-07 * CPI
         ZPMKSA = CPR(K) *B0EXP**2 / ZMU0
         ZPMKSA0 = CPR(1) *B0EXP**2 / ZMU0
         IF (ZPMKSA0 .EQ. 0.) ZPMKSA0 = 1.0E-10
         ZTEMPE0 = SQRT(RPEOP*ABS(ZPMKSA0)/1.E+19/1.602E-16)
         IF (AT4(1) .EQ. 0) THEN
C     T_E IN KEV
           ZTEMPE = ZTEMPE0 * ABS(ZPMKSA/ZPMKSA0)**(1./(1.+ETAEI))
           IF (ZTEMPE.EQ.0 .AND. K.NE.1) THEN
C     USE LAST BUT ONE FLUX SURFACE
             ZPMKSA = CPR(K-1) *B0EXP**2 / ZMU0
             ZTEMPE = ZTEMPE0 * ABS(ZPMKSA/ZPMKSA0)**(1./(1.+ETAEI))
           ENDIF
           IF (ZTEMPE .EQ. 0.0) ZTEMPE = ZTEMPE0
           IF (ZPMKSA .EQ. 0.0) ZPMKSA = ZPMKSA0
C     A2E=PE * 1/TE * DTE/DPSI
           ZA2E = ETAEI / (1 + ETAEI) * RPEOP*CPPR(K)
           ZA2I = ZA2E
         ELSE
           CALL POLYFUN(NSOUR,AT4,1,(1.-PS**2)*SPSIM,ZTEMPE,1,1)
           CALL POLYFUN(NSOUR,AT4,1,(1.-PS**2)*SPSIM,ZTEMPEP,2,1)
           IF (ZTEMPE .EQ. 0.0) ZTEMPE = ZTEMPE0
           ZA2E = RPEOP*CPR(K)*ZTEMPEP/ZTEMPE
           ZA2I = ZA2E
           IF (ZPMKSA.EQ.0.0 .AND. K.NE.1) ZPMKSA=CPR(K-1)*B0EXP**2/ZMU0
           IF (ZPMKSA .EQ. 0.0) ZPMKSA = ZPMKSA0
         ENDIF
C
C     TAUE (FORMULA IN CGS FOR TAUE AND NUESTAR)
C
         ZMASSE = 0.91094E-27
         ZCHARGE = 4.8032E-10
         ZERGKEV = 1000. * 1.6022E-12
         ZLNLAM = 20.39 - LOG(SQRT(ABS(ZPMKSA)*RPEOP/ZTEMPE)/ZTEMPE)
         ZVTE = SQRT(2.*ZTEMPE*ZERGKEV/ZMASSE)
         ZOTAUE = 32.*SQRT(CPI)*ABS(10.*ZPMKSA)*RPEOP*ZCHARGE**4*ZLNLAM
     +     / (3.*ZMASSE**3*ZVTE**5)
C
C     R(TET=ZTB2MIN), B_P(TET=ZTB2MIN), NUESTAR
C
         JTET0 = JB2MIN
         ZRTET0 = FQQQ0(ZRISO(JTET0),ZRISO(JTET0+1),ZRISO(JTET0+2)
     +     ,ZT2(JTET0),ZT2(JTET0+1),ZT2(JTET0+2),ZTB2MIN)
         ZBPOL0 = FQQQ0(ZGPISO(JTET0),ZGPISO(JTET0+1),ZGPISO(JTET0+2)
     +     ,ZT2(JTET0),ZT2(JTET0+1),ZT2(JTET0+2),ZTB2MIN) / ZRTET0
         RNUSTAR(K) = SQRT(2.)*(ZRTET0-RMAG)*R0EXP*100.*SQRT(ZB2MIN)
     +     *ZOTAUE / (ZBPOL0*ZSQRTA**3*ZVTE)
C
C     NON-ZERO NUESTAR CASE
C
         ZFTEF = ZFT
     +     / (1. + SQRT(RZION*RNUSTAR(K)) + 0.25*RNUSTAR(K)/RZION)
         ZFTEF2 = ZFTEF*ZFTEF
         ZFTEF3 = ZFTEF*ZFTEF2
         ZSIG = 1. - (1.+0.36/RZION)*ZFTEF + 0.59/RZION*ZFTEF2
     +     - 0.23/RZION*ZFTEF3
         ZL31 = -(1.+0.62/RZION)*ZFTEF + 0.84/RZION*ZFTEF2
     +     - 0.22/RZION*ZFTEF3
         ZL32 = (0.49+1.04*RZION)/RZION/(1.+0.44*RZION)*(ZFTEF-ZFTEF3)
     +     - (0.82+0.72*RZION)/RZION/(1.+0.245*RZION)*(ZFTEF2-ZFTEF3)
C
C     ZERO NUESTAR CASE
C
         ZFTEF = ZFT
         ZFTEF2 = ZFTEF*ZFTEF
         ZFTEF3 = ZFTEF*ZFTEF2
         ZSIG_0 = 1. - (1.+0.36/RZION)*ZFTEF + 0.59/RZION*ZFTEF2
     +     - 0.23/RZION*ZFTEF3
         ZL31_0 = -(1.+0.62/RZION)*ZFTEF + 0.84/RZION*ZFTEF2
     +     - 0.22/RZION*ZFTEF3
         ZL32_0 = (0.49+1.04*RZION)/RZION/(1.+0.44*RZION)*(ZFTEF-ZFTEF3)
     +     - (0.82+0.72*RZION)/RZION/(1.+0.245*RZION)*(ZFTEF2-ZFTEF3)
C
         ZALFA    = - 1.172 / (1. + 0.462 * ZFT/(1.-ZFT))
         ZA1 = CPPR(K)
         RJBSOS(K,1) = TMF(K)*
     +     (ZL31_0*(ZA1+ZALFA*(1./RPEOP-1.)*ZA2I) + ZL32_0 * ZA2E)
         RJBSOS(K,2) = TMF(K)*
     +     (ZL31*(ZA1+ZALFA*(1./RPEOP-1.)*ZA2I) + ZL32 * ZA2E)
c%OS
C        IF (K.EQ.1) WRITE(6,'("  K  PSI/PSIMAX  PE[CHEASE] ",
C    +     "  NE[E19]     TE[KEV]     L31*A1    L31*ALF*A2I ",
C    +     "  L32*A2E       q         EPS*B/BP")')
C        IF (K.EQ.1 .OR. (K/2)*2.EQ.K) THEN
C          WRITE(*,'(I3,1P10E12.4)') K,PS**2,RPEOP*CPR(K),
C    +       RPEOP*CPR(K)*B0EXP**2/ZMU0/1.E+19/1.602E-16/ZTEMPE,
C    +       ZTEMPE,ZL31_0*ZA1,ZL31_0*ZALFA*(1./RPEOP-1.)*ZA2I,
C    +       ZL32_0*ZA2E,QPSI(K),SQRT(ZB2MIN)/ZBPOL0/ARATIO(K)
C        ENDIF
C%OS
C
C     STOP IF P IS NEGATIF SOMEWHERE
C
         IF (CPR(K) .LT. 0) THEN
           WRITE(6,'(//,58("*"),/,5X,
     +       "ERROR: THE PRESSURE PROFILE IS NEGATIVE => STOP",/,
     +       58("*"))')
           WRITE(*,*) ' K=',K,' S =',PS,' PP=',CPR(K)
           CPR(K) = ABS(CPR(K))
         ENDIF
C
         RETURN
         END
C*DECK C2SM03
C*CALL PROCESS
         SUBROUTINE CHIPSI(NP1,KP)
C        #########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM03 INTERPOLATE THETA(S,CHI) AND SIGMA(S,CHI) WITH CUBIC         *
*        SPLINES, WHERE (S,CHI) IS THE ERATO STABILITY MESH           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   IC(NPCHI),    IC1(NPCHI),   IT0(NPCHI),   IT1(NPCHI),
     R   ZD2BCN(NTP2), ZD2BCO(NTP2), ZD2CHO(NTP2),
     R   ZD2SIG(NTP2), ZD2TET(NTP2), ZTET(NTP2),   
     R   ZA1(NTP2),    ZB1(NTP2),    ZC1(NTP2)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C 
         CHI(1) = CHI(1) + 2. * CPI
C
         CALL SCOPY(NT2,TETMAP(1,KP),1,ZTET,1)
C
         DO 1 J1=2,NT2
C
         IF (ZTET(J1) .LT. ZTET(J1-1)) THEN
C
            ZTET(J1) = ZTET(J1) + 2. * CPI * (1. +
     +                 INT(.5 * ABS(ZTET(J1) - ZTET(J1-1)) / CPI))       
C
         ENDIF
C
 1       CONTINUE
C
         CALL SPLCYP(CHIN(1,KP),CHIO(1,KP),NT1,RC2PI,RC2PI,
     &              ZD2CHO,ZA1,ZB1,ZC1)
         CALL SPLCY(CHIN(1,KP),BCHIN(1,KP),NT1,RC2PI,
     &              ZD2BCN,ZA1,ZB1,ZC1)
         CALL SPLCY(CHIN(1,KP),BCHIO(1,KP),NT1,RC2PI,
     &              ZD2BCO,ZA1,ZB1,ZC1)
         CALL SPLCY(CHIN(1,KP),SIGMAP(1,KP),NT1,RC2PI,
     &              ZD2SIG,ZA1,ZB1,ZC1)
         CALL SPLCYP(CHIN(1,KP),ZTET,NT1,RC2PI,RC2PI,
     &               ZD2TET,ZA1,ZB1,ZC1)
C
         ZD2CHO(NT2) = ZD2CHO(1)
         ZD2BCN(NT2) = ZD2BCN(1) 
         ZD2BCO(NT2) = ZD2BCO(1) 
         ZD2SIG(NT2) = ZD2SIG(1) 
         ZD2TET(NT2) = ZD2TET(1) 
C
         CALL RESETI(IC,NCHI,1)
         CALL RESETI(IC1,NCHI,1)
         DO 3 JG=1,NCHI
            DO 3 JT = 1,NT2
               IF (IC(JG).EQ.0) GOTO 2
               IT0(JG) = JT-1
               IF (CHIN(JT,KP).GE.CHIM(JG)) IC(JG)  = 0
 2             CONTINUE
               IF (KP.EQ.NP1) THEN
                  IF (IC1(JG).EQ.0) GOTO 3
                  IT1(JG) = JT-1
                  IF (CHIN(JT,KP).GE.CHI(JG)) IC1(JG) = 0
               ENDIF
 3       CONTINUE
C
         DO 4 J4=1,NCHI
C
         ICHIM = IT0(J4)
C
         IF (ICHIM .LT. 1)   ICHIM = 1
         IF (ICHIM .GT. NT1) ICHIM = NT1
C
         ZH = CHIN(ICHIM+1,KP) - CHIN(ICHIM,KP)
         ZA = (CHIN(ICHIM+1,KP) - CHIM(J4)) / ZH
         ZB = (CHIM(J4) - CHIN(ICHIM,KP)) / ZH
         ZC = (ZA + 1) * (ZA - 1) * ZH * 
     *        (CHIN(ICHIM+1,KP) - CHIM(J4)) / 6.
         ZD = (ZB + 1) * (ZB - 1) * ZH * 
     *        (CHIM(J4) - CHIN(ICHIM,KP)) / 6.
C 
         EQ13(J4,KP)  = ZA*BCHIN(ICHIM,KP) + ZB*BCHIN(ICHIM+1,KP) +
     +                   ZC*ZD2BCN(ICHIM)   + ZD*ZD2BCN(ICHIM+1)
         EQ22(J4,KP)  = ZA*BCHIO(ICHIM,KP) + ZB*BCHIO(ICHIM+1,KP) +
     +                   ZC*ZD2BCO(ICHIM)   + ZD*ZD2BCO(ICHIM+1)
         EQ24(J4,KP)  = ZA*CHIO(ICHIM,KP)  + ZB*CHIO(ICHIM+1,KP) +
     +                   ZC*ZD2CHO(ICHIM)   + ZD*ZD2CHO(ICHIM+1)
         TETCHI(J4,KP) = ZA*ZTET(ICHIM)     + ZB*ZTET(ICHIM+1) +
     +                   ZC*ZD2TET(ICHIM)   + ZD*ZD2TET(ICHIM+1)
C
         IF (TETCHI(J4,KP) .LT. CT(1))
     =                   TETCHI(J4,KP) = TETCHI(J4,KP) + 2.*CPI
         IF (TETCHI(J4,KP) .GT. CT(NT1)) 
     =                   TETCHI(J4,KP) = TETCHI(J4,KP) - 2.*CPI
C
         IF (KP .EQ. NP1) THEN
C
           SIGCHI(J4,NP1) = 1.
C
         ELSE
C
           SIGCHI(J4,KP) = ZA*SIGMAP(ICHIM,KP) + ZB*SIGMAP(ICHIM+1,KP) +
     +                     ZC*ZD2SIG(ICHIM)    + ZD*ZD2SIG(ICHIM+1)
C
        ENDIF
C
         IF (KP .EQ. NP1) THEN
C
            ICHIO = IT1(J4)
C
            IF (ICHIO .LT. 1)   ICHIO = 1
            IF (ICHIO .GT. NT1) ICHIO = NT1
C
            ZH = CHIN(ICHIO+1,KP) - CHIN(ICHIO,KP)
            ZA = (CHIN(ICHIO+1,KP) - CHI(J4)) / ZH
            ZB = (CHI(J4) - CHIN(ICHIO,KP)) / ZH
            ZC = (ZA + 1) * (ZA - 1) * ZH * 
     *           (CHIN(ICHIO+1,KP) - CHI(J4)) / 6.
            ZD = (ZB + 1) * (ZB - 1) * ZH * 
     *           (CHI(J4) - CHIN(ICHIO,KP)) / 6.
C 
            TETVAC(J4) = ZA*ZTET(ICHIO)    + ZB*ZTET(ICHIO+1) +
     +                   ZC*ZD2TET(ICHIO)  + ZD*ZD2TET(ICHIO+1)
            CHIOLD(J4) = ZA*CHIO(ICHIO,KP) + ZB*CHIO(ICHIO+1,KP) +
     +                   ZC*ZD2CHO(ICHIO)  + ZD*ZD2CHO(ICHIO+1)
C     
         ENDIF
C
 4       CONTINUE
C
         CHI(1)        = CHI(1) - 2.* CPI
         CHIOLD(NCHI1) = CHIOLD(1)
         CHIOLD(1)     = CHIOLD(1) - 2. * CPI
C
         IF (NIDEAL .EQ. 2 .AND. KP .EQ. NP1) THEN
C
C PLASMA BOUNDARY QUANTITIES REQUIRED BY LION ONLY
C
            TETVAC(NCHI1) = TETVAC(1)
C
            CALL SCOPY(NCHI,TETCHI(1,NP1),1,TETVACM,1)
C
            TETVACM(NCHI1) = TETVACM(1)
C
            CALL BOUND(NCHI1,TETVAC ,RHOVAC)
            CALL BOUND(NCHI1,TETVACM,RHOVACM)
C
            DO 15 J15=1,NCHI1
C
            ZR  = RHOVAC(J15) * COS(TETVAC(J15)) + R0
            ZZ  = RHOVAC(J15) * SIN(TETVAC(J15)) + RZ0
            ZRM = RHOVACM(J15) * COS(TETVACM(J15)) + R0
            ZZM = RHOVACM(J15) * SIN(TETVACM(J15)) + RZ0
C
            RHOVAC(J15)  = SQRT((ZR - RMAG)**2 + (ZZ - RZMAG)**2)
            TETVAC(J15)  = ATAN2(ZZ - RZMAG,ZR - RMAG)
            RHOVACM(J15) = SQRT((ZRM - RMAG)**2 + (ZZM - RZMAG)**2)
            TETVACM(J15) = ATAN2(ZZM - RZMAG,ZRM - RMAG)
C
            IF (J15 .NE. 1 .AND. TETVAC(J15) .LT. 0.) THEN
C
               TETVAC(J15) = TETVAC(J15) + 2. * CPI
C
            ENDIF
C
            IF (J15 .NE. 1 .AND. TETVACM(J15) .LT. 0.) THEN
C
               TETVACM(J15) = TETVACM(J15) + 2. * CPI
C
            ENDIF
C
  15        CONTINUE
C
            TETVACM(NCHI1) = TETVACM(1) + 2. * CPI
C
            CALL VLION
C
         ENDIF
C
         RETURN
         END
C*DECK C2SM04
C*CALL PROCESS
         SUBROUTINE ERDATA(KP)
C        #####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM04 COMPUTE EQ'S FOR ERATO (SEE [1], APPENDIX C2, TABLE 2 AND    *
*        SECTION 5.4.1 OF PUBLICATION)
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMPLO.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IS0(NPCHI),       IT0(NPCHI),       IC(NPCHI),
     R   ZBND(NPCHI,5),    ZCURV(NPCHI1),    ZDBDS(NPCHI,16),
     R   ZDBDT(NPCHI,16),  ZDBDST(NPCHI,16), ZD2BS2(NPCHI,16),
     R   ZD2BT2(NPCHI,16), ZPCEL(NPCHI,16),  ZS(NPCHI),
     R   ZS1(NPCHI),       ZS2(NPCHI),       ZTETA(NPCHI,5),
     R   ZT(NPCHI),        ZT1(NPCHI),       ZT2(NPCHI)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ZEPS  = 1.E-3
C
         DO 1 J1=1,NCHI
C
         ZTETA(J1,1) = TETCHI(J1,KP)
         ZTETA(J1,2) = TETCHI(J1,KP) - 2. * ZEPS
         ZTETA(J1,3) = TETCHI(J1,KP) -      ZEPS
         ZTETA(J1,4) = TETCHI(J1,KP) +      ZEPS
         ZTETA(J1,5) = TETCHI(J1,KP) + 2. * ZEPS
C
    1    CONTINUE
C
         CALL BOUND(NCHI,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NCHI,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NCHI,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NCHI,ZTETA(1,4),ZBND(1,4))
         CALL BOUND(NCHI,ZTETA(1,5),ZBND(1,5))
C
         CALL RESETI(IC,NCHI,1)
         DO 2 JT = 1,NT1
            DO 2 JG=1,NCHI
               IF (IC(JG).EQ.0) GOTO 2
               IT0(JG) = JT-1
               IF (TETCHI(JG,KP).LE.CT(JT)) IC(JG)  = 0
 2       CONTINUE
         CALL RESETI(IC,NCHI,1)
         DO 3 JS = 1,NS1
            DO 3 JG=1,NCHI
               IF (IC(JG).EQ.0) GOTO 3
               IS0(JG) = JS-1
               IF (SIGCHI(JG,KP).LE.CSIG(JS)) IC(JG)  = 0
 3       CONTINUE
C
         DO 4 J4=1,NCHI
C
         ZT(J4) = TETCHI(J4,KP)
         ZS(J4) = SIGCHI(J4,KP)
C
         IF (IS0(J4) .GT. NS) IS0(J4) = NS
         IF (IS0(J4) .LT. 1)  IS0(J4) = 1
         IF (IT0(J4) .GT. NT) IT0(J4) = NT
         IF (IT0(J4) .LT. 1)  IT0(J4) = 1
C
         ZS1(J4) = CSIG(IS0(J4))
         ZS2(J4) = CSIG(IS0(J4)+1)
         ZT1(J4) = CT(IT0(J4))
         ZT2(J4) = CT(IT0(J4)+1)
C
 4       CONTINUE
C
         CALL PSICEL(IS0,IT0,NCHI,NPCHI,ZPCEL,CPSICL)
         CALL BASIS3(NCHI,NPCHI,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT,
     &               ZDBDST,ZD2BS2,ZD2BT2)
C
         DO 5 J5=1,NCHI
C
         ZDRSDT = (ZBND(J5,2) + 8*(ZBND(J5,4) - ZBND(J5,3)) -
     -             ZBND(J5,5)) / (12. * ZEPS)
         ZD2RST = (- ZBND(J5,2) + 16. * ZBND(J5,3) -
     -             30. * ZBND(J5,1) + 16. * ZBND(J5,4) -
     -             ZBND(J5,5)) / (12. * ZEPS**2)
C
         ZDPDS = ZDBDS(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDS(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDS(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDS(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDS(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDS(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDS(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDS(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDS(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDS(J5,10) * ZPCEL(J5,10) +
     +           ZDBDS(J5,11) * ZPCEL(J5,11) +
     +           ZDBDS(J5,12) * ZPCEL(J5,12) +
     +           ZDBDS(J5,13) * ZPCEL(J5,13) +
     +           ZDBDS(J5,14) * ZPCEL(J5,14) +
     +           ZDBDS(J5,15) * ZPCEL(J5,15) +
     +           ZDBDS(J5,16) * ZPCEL(J5,16)
C
         ZDPDT = ZDBDT(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDT(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDT(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDT(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDT(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDT(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDT(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDT(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDT(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDT(J5,10) * ZPCEL(J5,10) +
     +           ZDBDT(J5,11) * ZPCEL(J5,11) +
     +           ZDBDT(J5,12) * ZPCEL(J5,12) +
     +           ZDBDT(J5,13) * ZPCEL(J5,13) +
     +           ZDBDT(J5,14) * ZPCEL(J5,14) +
     +           ZDBDT(J5,15) * ZPCEL(J5,15) +
     +           ZDBDT(J5,16) * ZPCEL(J5,16)
C
         ZD2PST = ZDBDST(J5, 1) * ZPCEL(J5, 1) +
     +            ZDBDST(J5, 2) * ZPCEL(J5, 2) +
     +            ZDBDST(J5, 3) * ZPCEL(J5, 3) +
     +            ZDBDST(J5, 4) * ZPCEL(J5, 4) +
     +            ZDBDST(J5, 5) * ZPCEL(J5, 5) +
     +            ZDBDST(J5, 6) * ZPCEL(J5, 6) +
     +            ZDBDST(J5, 7) * ZPCEL(J5, 7) +
     +            ZDBDST(J5, 8) * ZPCEL(J5, 8) +
     +            ZDBDST(J5, 9) * ZPCEL(J5, 9) +
     +            ZDBDST(J5,10) * ZPCEL(J5,10) +
     +            ZDBDST(J5,11) * ZPCEL(J5,11) +
     +            ZDBDST(J5,12) * ZPCEL(J5,12) +
     +            ZDBDST(J5,13) * ZPCEL(J5,13) +
     +            ZDBDST(J5,14) * ZPCEL(J5,14) +
     +            ZDBDST(J5,15) * ZPCEL(J5,15) +
     +            ZDBDST(J5,16) * ZPCEL(J5,16)
C
         ZD2PS2 = ZD2BS2(J5, 1) * ZPCEL(J5, 1) +
     +            ZD2BS2(J5, 2) * ZPCEL(J5, 2) +
     +            ZD2BS2(J5, 3) * ZPCEL(J5, 3) +
     +            ZD2BS2(J5, 4) * ZPCEL(J5, 4) +
     +            ZD2BS2(J5, 5) * ZPCEL(J5, 5) +
     +            ZD2BS2(J5, 6) * ZPCEL(J5, 6) +
     +            ZD2BS2(J5, 7) * ZPCEL(J5, 7) +
     +            ZD2BS2(J5, 8) * ZPCEL(J5, 8) +
     +            ZD2BS2(J5, 9) * ZPCEL(J5, 9) +
     +            ZD2BS2(J5,10) * ZPCEL(J5,10) +
     +            ZD2BS2(J5,11) * ZPCEL(J5,11) +
     +            ZD2BS2(J5,12) * ZPCEL(J5,12) +
     +            ZD2BS2(J5,13) * ZPCEL(J5,13) +
     +            ZD2BS2(J5,14) * ZPCEL(J5,14) +
     +            ZD2BS2(J5,15) * ZPCEL(J5,15) +
     +            ZD2BS2(J5,16) * ZPCEL(J5,16)
C
         ZD2PT2 = ZD2BT2(J5, 1) * ZPCEL(J5, 1) +
     +            ZD2BT2(J5, 2) * ZPCEL(J5, 2) +
     +            ZD2BT2(J5, 3) * ZPCEL(J5, 3) +
     +            ZD2BT2(J5, 4) * ZPCEL(J5, 4) +
     +            ZD2BT2(J5, 5) * ZPCEL(J5, 5) +
     +            ZD2BT2(J5, 6) * ZPCEL(J5, 6) +
     +            ZD2BT2(J5, 7) * ZPCEL(J5, 7) +
     +            ZD2BT2(J5, 8) * ZPCEL(J5, 8) +
     +            ZD2BT2(J5, 9) * ZPCEL(J5, 9) +
     +            ZD2BT2(J5,10) * ZPCEL(J5,10) +
     +            ZD2BT2(J5,11) * ZPCEL(J5,11) +
     +            ZD2BT2(J5,12) * ZPCEL(J5,12) +
     +            ZD2BT2(J5,13) * ZPCEL(J5,13) +
     +            ZD2BT2(J5,14) * ZPCEL(J5,14) +
     +            ZD2BT2(J5,15) * ZPCEL(J5,15) +
     +            ZD2BT2(J5,16) * ZPCEL(J5,16)
C
         ZFP    = (ZDPDS**2 + (ZDPDT / SIGCHI(J5,KP) - ZDPDS * ZDRSDT /
     -             ZBND(J5,1))**2) / ZBND(J5,1)**2
         ZGRADP = SQRT(ZFP)
C
         ZCOST = COS(ZTETA(J5,1))
         ZSINT = SIN(ZTETA(J5,1))
C
         ZRHO   = SIGCHI(J5,KP) * ZBND(J5,1)
         ZR     = ZRHO * ZCOST + R0
         ZZ     = ZRHO * ZSINT + RZ0
         ZJPHI  = - ZR * CPPR(KP) - TTP(KP) / ZR
         ZJAC   = CP(KP) * ZR**NER * ZGRADP**NEGP
C
         ZDSDR = (ZDRSDT * ZSINT + ZBND(J5,1) * ZCOST) / ZBND(J5,1)**2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBND(J5,1) * ZSINT - ZDRSDT * ZCOST) / ZBND(J5,1)**2
         ZDTDZ = ZCOST / ZRHO
C
         ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         Z1     = ZDPDS
         ZDZ1DS = ZD2PS2
         ZDZ1DT = ZD2PST
         Z2     = ZDPDT / SIGCHI(J5,KP) - ZDRSDT * ZDPDS / ZBND(J5,1)
         ZDZ2DS = (ZD2PST - ZDPDT / SIGCHI(J5,KP)) / SIGCHI(J5,KP) -
     -             ZDRSDT * ZD2PS2 / ZBND(J5,1)
         ZDZ2DT = - ZD2PST*ZDRSDT/ZBND(J5,1) + ZD2PT2/SIGCHI(J5,KP) +
     +              ZDPDS * (ZDRSDT**2-ZD2RST*ZBND(J5,1))/ZBND(J5,1)**2
C
         ZDFDS = 2 * (Z1 * ZDZ1DS + Z2 * ZDZ2DS) / ZBND(J5,1)**2
         ZDFDT = - 2 * ZFP * ZDRSDT / ZBND(J5,1) +
     +           2 * (Z1 * ZDZ1DT + Z2 * ZDZ2DT) / ZBND(J5,1)**2
C
         ZDFDR = ZDFDS * ZDSDR + ZDFDT * ZDTDR
         ZDFDZ = ZDFDS * ZDSDZ + ZDFDT * ZDTDZ
C
         ZDGDPN = (ZDFDR * ZDPDR + ZDFDZ * ZDPDZ) / ZFP
         ZDRDPN = ZDPDR / ZFP
         ZDZDPN = ZDPDZ / ZFP
C
         ZDRDCP = - ZJAC * ZDPDZ / ZR
         ZDZDCP =   ZJAC * ZDPDR / ZR
         ZDRDPC = ZDRDPN - .5*EQ13(J5,KP) * ZDRDCP / (CPSRF * CSM(KP))
         ZDZDPC = ZDZDPN - .5*EQ13(J5,KP) * ZDZDCP / (CPSRF * CSM(KP))
         ZDGDCP = ZJAC * (ZDPDR * ZDFDZ - ZDPDZ * ZDFDR) / ZR
         ZDGDPC = ZDGDPN - .5*EQ13(J5,KP) * ZDGDCP / (CPSRF * CSM(KP))
C
         Z3 = (2-NER)*ZDRDPC/ZR - .5*NEGP*ZDGDPC/ZFP - CPDP(KP)/CP(KP)
         Z4 = (2-NER)*ZDRDCP/ZR - .5*NEGP*ZDGDCP/ZFP
C
         IF (NDEQ .LT. 25) GO TO 999
C
         EQ( 1,J5,KP) = CS(KP)
         EQ( 2,J5,KP) = CHI(J5)
C
         IF (KP .LT. NPSI1) THEN
C
            EQ( 3,J5,KP) = CS(KP+1)
C
         ENDIF
C
         EQ( 4,J5,KP) = CHI(J5+1)
         EQ( 5,J5,KP) = CSM(KP)
         EQ( 6,J5,KP) = CHIM(J5)
C
C        EQ( 7,J5,KP) CONTAINS MASS DENSITY
C
         EQ( 7,J5,KP) = 1
         EQ( 8,J5,KP) = GAMMA * CPR(KP) / (Q0 * CPSRF)
         EQ( 9,J5,KP) = TMF(KP)
C
C        EQ(10,J5,KP) CONTAINS HELICAL PARAMETER
C
         EQ(10,J5,KP) = 0
         EQ(11,J5,KP) = QPSI(KP) / Q0
         EQ(12,J5,KP) = PSIISO(KP) * ZR**2 / (Q0 * ZFP)
C
C        EQ13(J5,KP)  IT IS COMPUTED BY SUBROUTINE CHIPSI.
C                     IT CONTAINS NEW BETACHI ON (CHI-MID,S-MID) NODES
C
         EQ(13,J5,KP) = EQ13(J5,KP)
         EQ(14,J5,KP) = ZR**2
         EQ(15,J5,KP) = 4 * CPSRF * CSM(KP) * ZDRDPC / ZR
         EQ(16,J5,KP) = 2 * ZDRDCP / ZR
         EQ(17,J5,KP) = 2 * CSM(KP) * CPSRF * (Z3 - ZR * ZJPHI / ZFP)
         EQ(18,J5,KP) = 2 * PSIISO(KP) * (ZJPHI * ZJPHI / ZFP -
     +                  .5 * ZJPHI * ZDGDPN / (ZR * ZFP) -
     -                  CPPR(KP) * ZDRDPN / ZR) / Q0
         EQ(19,J5,KP) = ZR**2 / ZJAC
         EQ(20,J5,KP) = Z4
         EQ(21,J5,KP) = 2 * CSM(KP) * CPSRF * CDQ(KP) * EQ24(J5,KP) +
     +                  QPSI(KP) * EQ22(J5,KP) -
     -                  TMF(KP) * ZJAC * EQ13(J5,KP) / ZR**2
C
C        EQ22(J5,KP)  IS COMPUTED BY SUBROUTINE CHIPSI
C                     IT CONTAINS OLD BETACHI ON (CHI-MID,S-MID) NODES
C
         EQ(22,J5,KP) = EQ22(J5,KP)
         EQ(23,J5,KP) = ZDGDCP
C
C        EQ24(J5,KP)  IS COMPUTED BY SUBROUTINE CHIPSI
C                     IT CONTAINS OLD CHI ON (CHI-MID,S-MID) NODES
C
         EQ(24,J5,KP) = EQ24(J5,KP)
         CALL ACOPY(NCHI-1,EQ(24,2,KP),NDEQ,EQ(25,1,KP),NDEQ)
         EQ(25,NCHI,KP) = 2. * CPI
C
C  QUANTITIES FOR LION ONLY
C
         IF (NIDEAL .EQ. 2) THEN
C
            EQ(NDEQ-3,J5,KP) = TTP(KP) / TMF(KP)
            EQ(NDEQ-2,J5,KP) = ZJPHI
            EQ(NDEQ-1,J5,KP) = .5 * ZDGDPN / ZFP
            EQ(NDEQ,J5,KP) = Z3
C
         ENDIF
C
 999     CONTINUE
C
         ZTMF2      = TMF(KP)**2
         ZCURV(J5)  = ZFP * (ZR * ZJPHI / (ZFP + ZTMF2) +
     +                ZDRDPN / ZR + .5*ZDGDPN / (ZFP + ZTMF2))
C
         CR(J5,KP)   = ZR
         CZ(J5,KP)   = ZZ
C
         IF (NIDEAL .NE. 0) THEN
C
            CNR1(J5,KP) = ZDPDR
            CNZ1(J5,KP) = ZDPDZ
C
         ELSE
C
            CNR1(J5,KP) = 2 * CSM(KP) * CPSRF * ZDRDPC
            CNZ1(J5,KP) = 2 * CSM(KP) * CPSRF * ZDZDPC
            CNR2(J5,KP) = ZDRDCP
            CNZ2(J5,KP) = ZDZDCP
C
         ENDIF
C
         RSHEAR(J5,KP) = (TTP(KP) / ZR**2 + ZTMF2 * (ZJPHI -
     -                    ZDGDPN / ZR) / (ZFP * ZR)) / TMF(KP)
C
 5       CONTINUE
C
         ZCURV(NCHI1) = ZCURV(1)
C
***********************************************************************
*                                                                     *
*  CORRECT EQ(21,J,KP) WHICH SHOULD BE THE INTEGRAL OF D(JT/R**2)/DPSI*
*                                                                     *
*     IF NCHI IS EVEN                                                 *
*     - THE INTEGRAL IS FROM 0 TO CHI IF (K .LE. NCHI/2+1) AND FROM 0 *
*       TO -CHI IF (K .GT. NCHI/2+1)                                  *
*                                                                     *
*     IF NCHI IS ODD                                                  *
*     - THE INTEGRAL IS FROM 0 TO CHI IF (K .LE. (NCHI+1)/2+1) AND    *
*       FROM 0 TO -CHI IF (K .GT. (NCHI+1)/2+1)                       *
*                                                                     *
***********************************************************************
C
         IF (NDEQ .GE. 25) THEN
C
         ZITJ0 = 4. * CPI * CSM(KP) * CPSRF * CDQ(KP)
C
         JU = NCHI / 2 + 2
C
         IF (MOD(NCHI,2) .EQ. 1) JU = (NCHI + 1) / 2 + 1
C
         DO 6 J6=JU,NCHI
C
         EQ(21,J6,KP) = EQ(21,J6,KP) - ZITJ0
C
    6    CONTINUE
C
         ENDIF
C
C-----------------------------------------------------------------------
C     COMPUTE (R,Z) OF ZERO CURVATURE LINE
C
         DO 7 J7=1,NCHI
C
           IF (NCURV .GE. 4*NPPSI1) GO TO 7
C
         IF (ZCURV(J7) * ZCURV(J7+1) .LE. 0.) THEN
C
            NCURV = NCURV + 1
C
            ZRHO1 = SIGCHI(J7,KP) * ZBND(J7,1)
            ZR1   = ZRHO1 * COS(ZTETA(J7,1)) + R0
            ZZ1   = ZRHO1 * SIN(ZTETA(J7,1)) + RZ0
C
            IF (J7 .NE. NCHI) THEN
C
               ZRHO2 = SIGCHI(J7+1,KP) * ZBND(J7+1,1)
               ZR2   = ZRHO2 * COS(ZTETA(J7+1,1)) + R0
               ZZ2   = ZRHO2 * SIN(ZTETA(J7+1,1)) + RZ0
C
            ELSE
C
               ZRHO2 = SIGCHI(1,KP) * ZBND(1,1)
               ZR2   = ZRHO2 * COS(ZTETA(1,1)) + R0
               ZZ2   = ZRHO2 * SIN(ZTETA(1,1)) + RZ0
C
            ENDIF
C
            RRCURV(NCURV) = (ZCURV(J7) * ZR2 - ZCURV(J7+1) * ZR1) /
     /                      (ZCURV(J7) - ZCURV(J7+1))
            RZCURV(NCURV) = (ZCURV(J7) * ZZ2 - ZCURV(J7+1) * ZZ1) /
     /                      (ZCURV(J7) - ZCURV(J7+1))
C
         ENDIF
C
 7    CONTINUE
C
         RETURN
         END
C*DECK C2SM05
C*CALL PROCESS
         SUBROUTINE CINT(K,PSIGMA,PTETA,PGWGT)
C        #####################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM05 EVALUATE FLUX SURFACE INTEGRALS REQUIRED FOR THE DEFINITION  *
*        OF I* AND I_PARALLEL (SEE EQ. (9) IN PUBLICATION)            *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMSOL.inc'
C
         PARAMETER (NPGMAP = NPMGS * NTP1)
         PARAMETER (ZEPS = 1.E-3)
C
         DIMENSION
     I   IS0(NPGMAP),        IT0(NPGMAP),       IC(NPGMAP),
     R   PSIGMA(*),          PTETA(*),          PGWGT(*),
     R   ZBND(NPGMAP),       ZBND1(NPGMAP,4),   ZDBDS(NPGMAP,16),  
     R   ZDBDT(NPGMAP,16),   ZPCEL(NPGMAP,16), 
     R   ZS1(NPGMAP),        ZS2(NPGMAP),      
     R   ZTETA(NPGMAP,4),    ZT1(NPGMAP),       ZT2(NPGMAP)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CID0(K) = 0.
         CIDR(K) = 0.
         CIDQ(K) = 0.
         CID2(K) = 0.
C
         IGMAP = NMGAUS * NT1
C
         CALL BOUND(NMGAUS*NT1,PTETA,ZBND)
C
         CALL RESETI(IC,IGMAP,1)
         DO 1 JT = 1,NT1
            DO 1 JG=1,IGMAP
               IF (IC(JG).EQ.0) GOTO 1
               IT0(JG) = JT-1
               IF (PTETA(JG).LE.CT(JT)) IC(JG)  = 0
 1       CONTINUE
         CALL RESETI(IC,IGMAP,1)
         DO 3 JS = 1,NS1
            DO 3 JG=1,IGMAP
               IF (IC(JG).EQ.0) GOTO 3
               IS0(JG) = JS-1
               IF (PSIGMA(JG).LE.CSIG(JS)) IC(JG)  = 0
 3       CONTINUE
C     
         DO 5 J5=1,IGMAP
C
         IF (IS0(J5).LT.  1) IS0(J5) = 1
         IF (IS0(J5).GE.NS1) IS0(J5) = NS
         IF (IT0(J5).LT.  1) IT0(J5) = 1
         IF (IT0(J5).GE.NT1) IT0(J5) = NT
C
         ZS1(J5) = CSIG(IS0(J5))
         ZS2(J5) = CSIG(IS0(J5)+1)
         ZT1(J5) = CT(IT0(J5))
 5       ZT2(J5) = CT(IT0(J5)+1)
C
         CALL PSICEL(IS0,IT0,IGMAP,NPGMAP,ZPCEL,CPSICL)
CYQL2018
         IF (NSTTP.LE.2 .OR. NSTTP.EQ.4) THEN
C
         CALL BASIS4(IGMAP,NPGMAP,ZS1,ZS2,ZT1,ZT2,PSIGMA,PTETA,ZDBDS)
C
         DO 6 J6=1,IGMAP
C
         ZDPDS = ZDBDS(J6, 1) * ZPCEL(J6, 1) +
     +           ZDBDS(J6, 2) * ZPCEL(J6, 2) +
     +           ZDBDS(J6, 3) * ZPCEL(J6, 3) +
     +           ZDBDS(J6, 4) * ZPCEL(J6, 4) +
     +           ZDBDS(J6, 5) * ZPCEL(J6, 5) +
     +           ZDBDS(J6, 6) * ZPCEL(J6, 6) +
     +           ZDBDS(J6, 7) * ZPCEL(J6, 7) +
     +           ZDBDS(J6, 8) * ZPCEL(J6, 8) +
     +           ZDBDS(J6, 9) * ZPCEL(J6, 9) +
     +           ZDBDS(J6,10) * ZPCEL(J6,10) +
     +           ZDBDS(J6,11) * ZPCEL(J6,11) +
     +           ZDBDS(J6,12) * ZPCEL(J6,12) +
     +           ZDBDS(J6,13) * ZPCEL(J6,13) +
     +           ZDBDS(J6,14) * ZPCEL(J6,14) +
     +           ZDBDS(J6,15) * ZPCEL(J6,15) +
     +           ZDBDS(J6,16) * ZPCEL(J6,16)
C
         ZCOST = COS(PTETA(J6))
         ZRHO  = PSIGMA(J6) * ZBND(J6)
         ZR    = ZRHO * ZCOST + R0
         ZINT0 = ZRHO * ZBND(J6) / ZDPDS
         ZINT1 = ZR * ZINT0
         ZINT2 = ZINT0 / ZR
c     write(*,*) 'j6 =',j6,' psigma=',psigma(j6),' pteta=',pteta(j6)
c    &       ,' zbnd=',zbnd(j6),' r0=',zr
C
         CID0(K) = CID0(K) + ZINT0 * PGWGT(J6)
         CIDR(K) = CIDR(K) + ZINT1 * PGWGT(J6)
         CIDQ(K) = CIDQ(K) + ZINT2 * PGWGT(J6)
C
 6       CONTINUE
C
         ELSE IF (NSTTP .EQ. 3) THEN
C
         CALL BASIS2(IGMAP,NPGMAP,ZS1,ZS2,ZT1,ZT2,PSIGMA,PTETA,
     +               ZDBDS,ZDBDT)
C
         DO 7 J7=1,IGMAP
C
         ZTETA(J7,1) = PTETA(J7) - 2. * ZEPS
         ZTETA(J7,2) = PTETA(J7) -      ZEPS
         ZTETA(J7,3) = PTETA(J7) +      ZEPS
         ZTETA(J7,4) = PTETA(J7) + 2. * ZEPS
C
 7       CONTINUE
C
         CALL BOUND(IGMAP,ZTETA(1,1),ZBND1(1,1))
         CALL BOUND(IGMAP,ZTETA(1,2),ZBND1(1,2))
         CALL BOUND(IGMAP,ZTETA(1,3),ZBND1(1,3))
         CALL BOUND(IGMAP,ZTETA(1,4),ZBND1(1,4))
C
         DO 8 J8=1,IGMAP
C
         ZDRSDT = (ZBND1(J8,1) + 8*(ZBND1(J8,3) - ZBND1(J8,2)) -
     -             ZBND1(J8,4)) / (12. * ZEPS)
C
         ZDPDS = ZDBDS(J8, 1) * ZPCEL(J8, 1) +
     +           ZDBDS(J8, 2) * ZPCEL(J8, 2) +
     +           ZDBDS(J8, 3) * ZPCEL(J8, 3) +
     +           ZDBDS(J8, 4) * ZPCEL(J8, 4) +
     +           ZDBDS(J8, 5) * ZPCEL(J8, 5) +
     +           ZDBDS(J8, 6) * ZPCEL(J8, 6) +
     +           ZDBDS(J8, 7) * ZPCEL(J8, 7) +
     +           ZDBDS(J8, 8) * ZPCEL(J8, 8) +
     +           ZDBDS(J8, 9) * ZPCEL(J8, 9) +
     +           ZDBDS(J8,10) * ZPCEL(J8,10) +
     +           ZDBDS(J8,11) * ZPCEL(J8,11) +
     +           ZDBDS(J8,12) * ZPCEL(J8,12) +
     +           ZDBDS(J8,13) * ZPCEL(J8,13) +
     +           ZDBDS(J8,14) * ZPCEL(J8,14) +
     +           ZDBDS(J8,15) * ZPCEL(J8,15) +
     +           ZDBDS(J8,16) * ZPCEL(J8,16)
C
         ZDPDT = ZDBDT(J8, 1) * ZPCEL(J8, 1) +
     +           ZDBDT(J8, 2) * ZPCEL(J8, 2) +
     +           ZDBDT(J8, 3) * ZPCEL(J8, 3) +
     +           ZDBDT(J8, 4) * ZPCEL(J8, 4) +
     +           ZDBDT(J8, 5) * ZPCEL(J8, 5) +
     +           ZDBDT(J8, 6) * ZPCEL(J8, 6) +
     +           ZDBDT(J8, 7) * ZPCEL(J8, 7) +
     +           ZDBDT(J8, 8) * ZPCEL(J8, 8) +
     +           ZDBDT(J8, 9) * ZPCEL(J8, 9) +
     +           ZDBDT(J8,10) * ZPCEL(J8,10) +
     +           ZDBDT(J8,11) * ZPCEL(J8,11) +
     +           ZDBDT(J8,12) * ZPCEL(J8,12) +
     +           ZDBDT(J8,13) * ZPCEL(J8,13) +
     +           ZDBDT(J8,14) * ZPCEL(J8,14) +
     +           ZDBDT(J8,15) * ZPCEL(J8,15) +
     +           ZDBDT(J8,16) * ZPCEL(J8,16)
C
         ZCOST = COS(PTETA(J8))
         ZRHO  = PSIGMA(J8) * ZBND(J8)
         ZR    = ZRHO * ZCOST + R0
         ZINT1 = ZR * ZRHO * ZBND(J8) / ZDPDS
         ZINT2 = ZRHO * ZBND(J8) / (ZR * ZDPDS)
         ZINT3 = (PSIGMA(J8) * ZDPDS + (ZDPDT / ZDPDS - PSIGMA(J8) *
     *            ZDRSDT / ZBND(J8)) * (ZDPDT / PSIGMA(J8) - 
     -            ZDPDS * ZDRSDT / ZBND(J8))) / ZR
C
         CIDR(K) = CIDR(K) + ZINT1 * PGWGT(J8)
         CIDQ(K) = CIDQ(K) + ZINT2 * PGWGT(J8)
         CID2(K) = CID2(K) + ZINT3 * PGWGT(J8)
C
 8       CONTINUE
C
         ENDIF
C
         IF (NSTTP .LE. 2) THEN
C
            CID0(K) = CID0(K) / CIDQ(K)
            CID2(K) = CIDR(K) / CIDQ(K)
C
         ELSE IF (NSTTP .EQ. 3) THEN
C
            CID0(K) = CIDR(K) / CIDQ(K)
            CID2(K) = CID2(K) / CIDQ(K)
CYQL2018
         ELSE IF (NSTTP .EQ. 4) THEN
            CID2(K) = 1./CIDQ(K)
 
         ENDIF
C
         RETURN
         END
C*DECK C2SM06
C*CALL PROCESS
         SUBROUTINE PREMAP(K)
C        ####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM06 LEAD SET-UP OF STABILITY S-MESH, TRACING OF CONSTANT FLUX    *
*        SURFACES AND EVALUATION OF PROFILES                          *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IC(2*NPISO),    IS0(2*NPISO),
     R   ZPAR(2*NPISO),  ZS(2*NPISO)
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(CID0,2*NPISO)
         CALL VZERO(CIDQ,2*NPISO)
         CALL VZERO(CIDR,2*NPISO)
         CALL VZERO(CID2,2*NPISO)
         CALL VZERO(SIGPSI,2*NPMGS*NTP1*NPISO)
         CALL VZERO(TETPSI,2*NPMGS*NTP1*NPISO)
         CALL VZERO(WGTPSI,2*NPMGS*NTP1*NPISO)
C
         IF (K .EQ. 1) THEN
C
            IN = NISO
C
            DO 1 J1=1,NISO
C
            PSIISO(J1) = SPSIM * (1. - CSIPR(J1)**2)
C
    1       CONTINUE
C
         ELSE IF (K .EQ. 2) THEN
C
            IN = NPPR+1
C
            DO 2 J2=1,NPPR+1
C
            PSIISO(J2) = SPSIM * (1. - PCSM(J2)**2)
C
    2       CONTINUE
C
         ELSE IF (K .EQ. 3) THEN
C
            IN = NPSI1
C
            DO 3 J3=1,NPSI1
C
            PSIISO(J3) = SPSIM * (1. - PSISCL * CSM(J3)**2)
            ZS(J3)     = CSM(J3)
C
    3       CONTINUE
C
         ELSE IF (K .EQ. 4) THEN
C
            IN = 2 * NPSI
C
            DO 4 J4=1,NPSI
C
            PSIISO(2*(J4-1)+1) = SPSIM * (1. - PSISCL * CSM(J4)**2)
            PSIISO(2*J4      ) = SPSIM * (1. - PSISCL * CS(J4+1)**2)
            ZS(2*(J4-1)+1)     = CSM(J4)
            ZS(2*J4)           = CS(J4+1)
C
    4       CONTINUE
C
         ELSE IF (K .EQ. 5) THEN
C
            IN = NPSI * (NMGAUS + 2)
            I1 = 1
            I2 = 1
C
            DO 5 J5=1,IN
C
            IF (MOD(J5,NMGAUS+2) .EQ. 3) THEN
C
               PSIISO(J5) = SPSIM * (1. - PSISCL * CSM(I1)**2)
               ZS(J5)     = CSM(I1)
               I1         = I1 + 1
C
            ELSE
C
               PSIISO(J5) = SPSIM * (1. - PSISCL * CSPEN(I2)**2)
               ZS(J5)     = CSPEN(I2)
               I2         = I2 + 1
C
            ENDIF
C
    5       CONTINUE
C
         ENDIF
C
         IF (NSTTP .EQ. 3 .AND. K .GT. 2) THEN
C
           CALL RESETI(IC,IN,1)
           DO 6 JS = 1,NISO
             DO 6 JG=1,IN
               IF (IC(JG).EQ.0) GOTO 6
               IS0(JG) = JS-1
               IF (ZS(JG).LE.CSIPR(JS)) IC(JG)  = 0
 6          CONTINUE
C
            DO 7 J7=1,IN
C
            ISIPR = IS0(J7)
C
            IF (ISIPR .GE. NISO - 2) ISIPR = NISO - 2
            IF (ISIPR .LE. 2)        ISIPR = 2
C
            TMF(J7) = FCCCC0(TMFO(ISIPR-1),TMFO(ISIPR),
     ,                       TMFO(ISIPR+1),TMFO(ISIPR+2),
     ,                       CSIPR(ISIPR-1),CSIPR(ISIPR),
     ,                       CSIPR(ISIPR+1),CSIPR(ISIPR+2),
     ,                       ZS(J7))
C
 7          CONTINUE
C
         ENDIF
C
         CALL RMRAD(IN,SPSIM,RC0P,RC0P,ZPAR,SIGMAP,TETMAP,NTP2)
         CALL PROFILE(IN)
C
         RETURN
         END
C*DECK C2SM07
C*CALL PROCESS
         SUBROUTINE GCHI(K)
C        ##################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM07 INTERPOLATE CHI(S,THETA) AND BETA_CHI(S,THETA) AT GAUSSIAN   *
*        QUADRATURE POINTS ALONG CONSTANT POLOIDAL FLUX SURFACES.     *
*        THESE QUANTITIES ARE REQUIRED FOR THE EVALUATION OF THE EQ'S *
*        FOR MARS                                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   ZA1(NTP2),         ZB1(NTP2),    ZC1(NTP2),
     R   ZD2BCN(NTP2),      ZD2CHN(NTP2), ZTET(NTP2), 
     R   ZTET1(NTP2*NPMGS)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C
         CALL SCOPY(NT2,TETMAP(1,K),1,ZTET,1)
         CALL SCOPY(NMGAUS*NT1,TETPSI(1,K),1,ZTET1,1)
C
         DO 1 J1=2,NT2
C
         IF (ZTET(J1) .LT. ZTET(J1-1)) THEN
C
            ZTET(J1) = ZTET(J1) + 2. * CPI * (1. +
     +                 INT(.5 * ABS(ZTET(J1) - ZTET(J1-1)) / CPI))        
C
         ENDIF
C
    1    CONTINUE
C
         DO 3 J3=1,NMGAUS
C
         DO 2 J2=1,NT1
C
         IT = (J2 - 1) * NMGAUS + J3
C
         IF (ZTET1(IT) .LT. ZTET(J2)) THEN
C
            ZTET1(IT) = ZTET1(IT) + 2. * CPI * (1. +
     +                  INT(.5 * ABS(ZTET1(IT) - ZTET(J2)) / CPI))       
C
         ELSE IF (ZTET1(IT) .GT. ZTET(J2+1)) THEN
C
            ZTET1(IT) = ZTET1(IT) - 2. * CPI * (1. +
     +                  INT(.5 * ABS(ZTET1(IT) - ZTET(J2+1)) / CPI))       
C
         ENDIF
C
    2    CONTINUE
    3    CONTINUE
C
         CALL SPLCYP(ZTET,CHIN(1,K),NT1,RC2PI,RC2PI,
     &               ZD2CHN,ZA1,ZB1,ZC1)
         CALL SPLCY(ZTET,BCHIN(1,K),NT1,RC2PI,
     &              ZD2BCN,ZA1,ZB1,ZC1)
C
         ZD2CHN(NT2) = ZD2CHN(1)
         ZD2BCN(NT2) = ZD2BCN(1) 
C
         DO 5 J5=1,NMGAUS
C
         DO 4 J4=1,NT1
C
         IT = (J4 - 1) * NMGAUS + J5
C
         ZH = ZTET(J4+1) - ZTET(J4)
         ZA = (ZTET(J4+1) - ZTET1(IT)) / ZH
         ZB = (ZTET1(IT) - ZTET(J4)) / ZH
         ZC = (ZA + 1) * (ZA - 1) * ZH * 
     *        (ZTET(J4+1) - ZTET1(IT)) / 6.
         ZD = (ZB + 1) * (ZB - 1) * ZH * 
     *        (ZTET1(IT) - ZTET(J4)) / 6.
C 
         CHIISO(IT) = ZA * CHIN(J4,K)  + ZB * CHIN(J4+1,K) +
     +                ZC * ZD2CHN(J4)  + ZD * ZD2CHN(J4+1)
         BCHISO(IT) = ZA * BCHIN(J4,K) + ZB * BCHIN(J4+1,K) +
     +                ZC * ZD2BCN(J4)  + ZD * ZD2BCN(J4+1)
C
    4    CONTINUE
    5    CONTINUE
C
         RETURN
         END
C*DECK C2SM08
C*CALL PROCESS
         SUBROUTINE GIJLIN(KPSI,PS)
C        ##########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM08 COMPUTE EQ'S FOR MARS IN [1], TABLE 3 AT GAUSSIAN QUADRATURE *
*        POINTS ALONG CONSTANT POLOIDAL FLUX SURFACES.                *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C
cab
cab      abs(dpsi/ds)
cab
         ZDPSIS = 2. * PS * CPSRF
         CALL GCHI(KPSI)
cab
         DO 1 J1=1,NMGAUS*NT1
C
         ZR = RRISO(J1,KPSI)
C
C ZJAC = [GRAD-S; GRAD-CHI;GRAD-PHI]
C ZGRADS = |GRAD-S|
C ZBSCHI = BETA_{S,CHI}
C ZGCHI2 = |GRAD-CHI|**2
C
         ZJAC   = ZDPSIS * CP(KPSI) * ZR**NER *
     *            GPISO(J1,KPSI)**NEGP
         ZGRADS = GPISO(J1,KPSI) / ZDPSIS
         ZBSCHI = BCHISO(J1)
         ZGCHI2 = (ZBSCHI * ZGRADS)**2 + (ZR / (ZJAC * ZGRADS))**2
C
C SEE DEFINITIONS OF EQL'S IN H.LUETJENS ET AL., COMPUTER PHYSICS
C COMMUNICATIONS 69, 287 (1992)
C
         EQL(J1, 1) = ZJAC*(ZBSCHI*ZGRADS/ZR)**2 + 1./(ZJAC*ZGRADS**2)
         EQL(J1, 2) = ZJAC * (ZGRADS / ZR)**2
         EQL(J1, 3) = ZR**2 / ZJAC
         EQL(J1, 4) = - ZBSCHI * ZJAC * (ZGRADS / ZR)**2
         EQL(J1, 5) = EQL(J1, 1) * ZJAC**2
         EQL(J1, 6) = EQL(J1, 2) * ZJAC**2
         EQL(J1, 7) = EQL(J1, 3) * ZJAC**2
         EQL(J1, 8) = EQL(J1, 4) * ZJAC**2
         EQL(J1, 9) = ZJAC
         EQL(J1,10) = - ZDPSIS * TTP(KPSI) / TMF(KPSI)
         EQL(J1,11) = -ZJAC * (CPPR(KPSI)+TTP(KPSI)/ZR**2)
         EQL(J1,12) = ZDPSIS
         EQL(J1,13) = ZJAC * TMF(KPSI) / ZR**2
         EQL(J1,14) = CPR(KPSI)
         EQL(J1,15) = ZDPSIS * CPPR(KPSI)
C
C ZDPDR = D(PSI) / D(R) AT Z CONSTANT
C ZDPDZ = D(PSI) / D(Z) AT R CONSTANT
C
         ZDPDR = DPRISO(J1,KPSI)
         ZDPDZ = DPZISO(J1,KPSI)
C
C ZDSDR = D(S)/D(R) AT Z CONSTANT
C ZDSDZ = D(S)/D(Z) AT R CONSTANT
C ZDCDR = D(CHI)/D(R) AT Z CONSTANT
C ZDCDZ = D(CHI)/D(Z) AT R CONSTANT
C
         ZDSDR = ZDPDR/ZDPSIS
         ZDSDZ = ZDPDZ/ZDPSIS
         ZDCDR = ZDSDR*ZBSCHI-ZDSDZ*ZR/(ZJAC*ZGRADS**2)
         ZDCDZ = ZDSDZ*ZBSCHI+ZDSDR*ZR/(ZJAC*ZGRADS**2)
C
C ZB2 = B**2
C ZBR = B_R
C ZBZ = B_Z
C
         ZB2 = (TMF(KPSI)**2+GPISO(J1,KPSI)**2) / ZR**2
         ZBR = -ZDPDZ/ZR
         ZBZ =  ZDPDR/ZR
cab
cab      EQL(*,16) IS D(CHI)/DZ
cab      EQL(*,17) IS DS/DZ 
cab      EQL(*,18) IS BZ
cab      EQL(*,19) IS BR
cab
         EQL(J1,16) = ZDCDZ
         EQL(J1,17) = ZDSDZ
         EQL(J1,18) = ZBZ
         EQL(J1,19) = ZBR
CLIUYQ
         EQL(J1,20) = 1./ZJAC
         EQL(J1,21) = ZR
         EQL(J1,22) = RZISO(J1,KPSI)
C
C QUANTITIES FOR SHEARED FLOWS, INHOMOGEN DENSITY, ETC...
C
C ZDRDCS = D(R)/D(CHI) AT S   CONSTANT
C ZDRDSC = D(R)/D(S)   AT CHI CONSTANT
C ZDGDCS = D(GRAD-PSI**2)/D(CHI) AT S   CONSTANT
C ZDGDSC = D(GRAD-PSI**2)/D(S)   AT CHI CONSTANT
C ZDJACDS = D(JACOBIAN)/DS AT CHI CONSTANT
C
         ZDRDCS = - ZJAC * ZDSDZ / ZR
         ZDRDSC = ZDPSIS*DRNISO(J1,KPSI) -  ZBSCHI * ZDRDCS
         ZDGDCS = ZJAC*(ZDSDR * DGZISO(J1,KPSI) - 
     &                  ZDSDZ * DGRISO(J1,KPSI))/ZR
         ZDGDSC = ZDPSIS*DGNISO(J1,KPSI) - ZBSCHI * ZDGDCS
C
         ZDJACDS = ZJAC*(1./PS+ZDPSIS*CPDP(KPSI)/CP(KPSI)+
     &             NER*ZDRDSC/ZR+.5*NEGP*ZDGDSC/GPISO(J1,KPSI)**2)
C
C OM2RT(KPSI) = OMEGA(S)**2/(K-BOLTZMANN*TEMPERATURE(S))
C DOM2RTDP(KPSI) = D(OM2RT(KPSI))/D(S)
C
         ZGMUNU = 1.
         ZRHO   = 1.
         ZDRHODS= 0.
         ZOMEG  = 0.
         ZDOMDS = 0.
         OM2RT  = 0.
         DOM2RTDS = 0.
C
C FLEXP = EXP(R**2*OMEGA**2/(2.*K-BOLTZMANN*TEMPERATURE(S)))
C DFLEXPDS = D(FLEXP)/D(S)   AT CHI CONSTANT
C DFLEXPDC = D(FLEXP)/D(CHI) AT S   CONSTANT
C
         FLEXP  = EXP(.5*ZR**2*OM2RT)
         DFLEXPDS = FLEXP*ZR*(ZDRDSC*OM2RT+
     &                        .5*ZR*DOM2RTDS)
         DFLEXPDC = FLEXP*ZR*ZDRDCS*OM2RT
C
         EQI(J1, 1) = FLEXP * TMF(KPSI) * ZJAC / ZB2
         EQI(J1, 2) = -FLEXP*ZJAC*EQL(J1,2)*EQL(J1,12)/ZB2
         EQI(J1, 3) = -FLEXP*EQL(J1,4)*EQL(J1,12)**2/ZB2
         EQI(J1, 4) = -FLEXP*EQL(J1,4)*EQL(J1,13)*EQL(J1,12)/ZB2
         EQI(J1, 5) = FLEXP*ZJAC*(ZJAC*EQL(J1,1)-
     &                            (EQL(J1,4)*EQL(J1,12))**2/ZB2)
         EQI(J1, 6) = FLEXP*TMF(KPSI)*ZJAC**2*EQL(J1,4)/ZB2
         EQI(J1, 7) = FLEXP*(ZJAC*ZR)**2*EQL(J1,2)/ZB2
         EQI(J1, 8) = FLEXP*ZJAC*ZB2
         EQI(J1, 9) = EQL(J1,4)*EQL(J1,12)*EQL(J1,15)/ZB2
         EQI(J1,10) = -ZJAC*(TTP(KPSI)/TMF(KPSI)+
     &                       TMF(KPSI)*CPPR(KPSI)/ZB2)
         EQI(J1,11) = - EQI(J1,3)
         EQI(J1,12) = - EQI(J1,4)
         EQI(J1,13) = - EQI(J1,1)
         EQI(J1,14) = - EQI(J1,2)
         EQI(J1,15) = FLEXP*ZDOMDS*TMF(KPSI)*ZJAC*EQL(J1,4)*
     &                EQL(J1,12)/ZB2
         EQI(J1,16) = 2.*FLEXP*ZOMEG*ZJAC**2*ZDCDZ*TMF(KPSI)/ZR**2
         EQI(J1,17) = FLEXP*ZDOMDS*ZR**2*ZJAC*
     &                EQL(J1,2)*EQL(J1,12)/ZB2
         EQI(J1,18) = 2.*FLEXP*ZOMEG*ZJAC**2*ZBZ/ZB2
         EQI(J1,19) = -FLEXP*ZDOMDS*TMF(KPSI)*ZJAC
         EQI(J1,20) = 2*FLEXP*ZOMEG*ZJAC**2*ZDSDZ
         EQI(J1,21) = FLEXP*ZJAC*EQL(J1,4)*EQL(J1,12)
         EQI(J1,22) = FLEXP*ZJAC**2*TMF(KPSI)/EQL(J1,12)
         EQI(J1,23) = -ZGMUNU*ZJAC*(ZJAC*ZDCDZ/ZR-
     &                 ZBR*EQL(J1,12)*EQL(J1,4)/ZB2)**2
         EQI(J1,24) = ZJAC*ZGMUNU*(ZJAC*ZDCDZ/ZR-
     &                ZBR*EQL(J1,12)*EQL(J1,4)/ZB2)*
     &                (ZDSDZ*EQL(J1,13)*ZR/ZB2)
         EQI(J1,25) = EQI(J1,24)
         EQI(J1,26) = -ZGMUNU*ZJAC*(EQL(J1,13)*ZR*ZDSDZ/ZB2)**2
         EQI(J1,27) = ZJAC*EQL(J1,2)*EQL(J1,3)*EQL(J1,13)/ZB2
         EQI(J1,28) = ZJAC*EQL(J1,2)*EQL(J1,3)*EQL(J1,12)/ZB2
         EQI(J1,29) = EQL(J1,3)*EQL(J1,4)*EQL(J1,13)**2/ZB2
         EQI(J1,30) = EQL(J1,3)*EQL(J1,4)*EQL(J1,12)*EQL(J1,13)/ZB2
         EQI(J1,31) = ZJAC*EQL(J1,3)*EQL(J1,4)*EQL(J1,13)/ZB2
         EQI(J1,32) = ZJAC*EQL(J1,1)-(EQL(J1,4)*EQL(J1,12))**2/ZB2
C
         EQ3(J1, 1) = ZRHO
         EQ3(J1, 2) = ZDRHODS
         EQ3(J1, 3) = ZOMEG
         EQ3(J1, 4) = ZDOMDS
         EQ3(J1, 5) = FLEXP
         EQ3(J1, 6) = ZR*ZOMEG**2*ZJAC*(EQL(J1,5)*ZDSDR+
     &                EQL(J1,8)*ZDCDR-ZJAC*ZBR*EQL(J1,12)*EQL(J1,4)/ZB2)
         EQ3(J1, 7) = ZR*ZOMEG**2*TMF(KPSI)*(EQL(J1,6)*ZDCDR+
     &                                       EQL(J1,8)*ZDSDR)/ZB2
         EQ3(J1, 8) = ZR*ZOMEG**2*ZJAC*ZBR
         EQ3(J1, 9) = -EQL(J1,15)
         EQ3(J1,10) = ZJAC*FLEXP
         EQ3(J1,11) = FLEXP*EQL(J1,12)
         EQ3(J1,12) = FLEXP*EQL(J1,13)
         EQ3(J1,13) = ZJAC*DFLEXPDS
         EQ3(J1,14) = DFLEXPDC*EQL(J1,12)**2*EQL(J1,4)/ZB2
         EQ3(J1,15) = ZJAC*DFLEXPDC*TMF(KPSI)/ZB2
         EQ3(J1,16) = DFLEXPDC*EQL(J1,12)
         EQ3(J1,17) = FLEXP*ZDJACDS
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SM09
C*CALL PROCESS
         SUBROUTINE FOURIER(KPSI,PSIGMA,PTETA,PGWGT,KM)
C        ##############################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM09 FOURIER TRANSFORM THE EQ'S COMPUTED IN GIJLIN ACCORDING TO   *
*        [1], EQ. (22)                                                *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
CLIUYQ
         DIMENSION
     R   PSIGMA(*),   PTETA(*),   PGWGT(*),
     R   ZINTI(22),   ZINTR(22),  ZINT2I(32),  ZINT2R(32),
     R   ZINT3I(17),  ZINT3R(17)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(ZINTI,22)
         CALL VZERO(ZINTR,22)
         CALL VZERO(ZINT2I,32)
         CALL VZERO(ZINT2R,32)
         CALL VZERO(ZINT3I,17)
         CALL VZERO(ZINT3R,17)
C
         DO 1 J1=1,NMGAUS*NT1
C
         ZR   = RRISO(J1,KPSI)
         ZGP  = GPISO(J1,KPSI)
         ZRHO = RHOISO(J1,KPSI)
         ZBND = BNDISO(J1,KPSI)
         ZJAC = CP(KPSI) * ZR**(NER-1) * ZGP**NEGP * DPSISO(J1,KPSI)
         ZX   = .5 * PGWGT(J1) * ZRHO * ZBND / (ZJAC * CPI)
         ZARG = KM * CHIISO(J1)

         ZCOEFR =   ZX * COS(ZARG)
         ZCOEFI = - ZX * SIN(ZARG)

CLIUYQ:  compute (RM,RZ) in theta angle, not in chi-angle
C        only for NFTC plot
         ZARGT   = KM * TETPSI(J1,KPSI)
         ZCOEFRT =   .5 * PGWGT(J1) * COS(ZARGT) / CPI
         ZCOEFIT = - .5 * PGWGT(J1) * SIN(ZARGT) / CPI
C
         ZINTR( 1) = ZINTR( 1) + ZCOEFR * EQL(J1, 1) 
         ZINTR( 2) = ZINTR( 2) + ZCOEFR * EQL(J1, 2) 
         ZINTR( 3) = ZINTR( 3) + ZCOEFR * EQL(J1, 3) 
         ZINTR( 4) = ZINTR( 4) + ZCOEFR * EQL(J1, 4) 
         ZINTR( 5) = ZINTR( 5) + ZCOEFR * EQL(J1, 5) 
         ZINTR( 6) = ZINTR( 6) + ZCOEFR * EQL(J1, 6) 
         ZINTR( 7) = ZINTR( 7) + ZCOEFR * EQL(J1, 7) 
         ZINTR( 8) = ZINTR( 8) + ZCOEFR * EQL(J1, 8) 
         ZINTR( 9) = ZINTR( 9) + ZCOEFR * EQL(J1, 9) 
         ZINTR(10) = ZINTR(10) + ZCOEFR * EQL(J1,10) 
         ZINTR(11) = ZINTR(11) + ZCOEFR * EQL(J1,11) 
         ZINTR(12) = ZINTR(12) + ZCOEFR * EQL(J1,12) 
         ZINTR(13) = ZINTR(13) + ZCOEFR * EQL(J1,13) 
         ZINTR(14) = ZINTR(14) + ZCOEFR * EQL(J1,14) 
         ZINTR(15) = ZINTR(15) + ZCOEFR * EQL(J1,15) 
         ZINTR(16) = ZINTR(16) + ZCOEFR * EQL(J1,16) 
         ZINTR(17) = ZINTR(17) + ZCOEFR * EQL(J1,17) 
         ZINTR(18) = ZINTR(18) + ZCOEFR * EQL(J1,18) 
         ZINTR(19) = ZINTR(19) + ZCOEFR * EQL(J1,19) 
         ZINTR(20) = ZINTR(20) + ZCOEFR * EQL(J1,20) 
         ZINTR(21) = ZINTR(21) + ZCOEFR * EQL(J1,21) 
         ZINTR(22) = ZINTR(22) + ZCOEFR * EQL(J1,22) 
C
         ZINT2R( 1) = ZINT2R( 1) + ZCOEFR * EQI(J1, 1) 
         ZINT2R( 2) = ZINT2R( 2) + ZCOEFR * EQI(J1, 2) 
         ZINT2R( 3) = ZINT2R( 3) + ZCOEFR * EQI(J1, 3) 
         ZINT2R( 4) = ZINT2R( 4) + ZCOEFR * EQI(J1, 4) 
         ZINT2R( 5) = ZINT2R( 5) + ZCOEFR * EQI(J1, 5) 
         ZINT2R( 6) = ZINT2R( 6) + ZCOEFR * EQI(J1, 6) 
         ZINT2R( 7) = ZINT2R( 7) + ZCOEFR * EQI(J1, 7) 
         ZINT2R( 8) = ZINT2R( 8) + ZCOEFR * EQI(J1, 8) 
         ZINT2R( 9) = ZINT2R( 9) + ZCOEFR * EQI(J1, 9) 
         ZINT2R(10) = ZINT2R(10) + ZCOEFR * EQI(J1,10) 
         ZINT2R(11) = ZINT2R(11) + ZCOEFR * EQI(J1,11) 
         ZINT2R(12) = ZINT2R(12) + ZCOEFR * EQI(J1,12) 
         ZINT2R(13) = ZINT2R(13) + ZCOEFR * EQI(J1,13) 
         ZINT2R(14) = ZINT2R(14) + ZCOEFR * EQI(J1,14) 
         ZINT2R(15) = ZINT2R(15) + ZCOEFR * EQI(J1,15) 
         ZINT2R(16) = ZINT2R(16) + ZCOEFR * EQI(J1,16) 
         ZINT2R(17) = ZINT2R(17) + ZCOEFR * EQI(J1,17) 
         ZINT2R(18) = ZINT2R(18) + ZCOEFR * EQI(J1,18) 
         ZINT2R(19) = ZINT2R(19) + ZCOEFR * EQI(J1,19) 
         ZINT2R(20) = ZINT2R(20) + ZCOEFR * EQI(J1,20) 
         ZINT2R(21) = ZINT2R(21) + ZCOEFR * EQI(J1,21) 
         ZINT2R(22) = ZINT2R(22) + ZCOEFR * EQI(J1,22) 
         ZINT2R(23) = ZINT2R(23) + ZCOEFR * EQI(J1,23) 
         ZINT2R(24) = ZINT2R(24) + ZCOEFR * EQI(J1,24) 
         ZINT2R(25) = ZINT2R(25) + ZCOEFR * EQI(J1,25) 
         ZINT2R(26) = ZINT2R(26) + ZCOEFR * EQI(J1,26) 
         ZINT2R(27) = ZINT2R(27) + ZCOEFR * EQI(J1,27) 
         ZINT2R(28) = ZINT2R(28) + ZCOEFR * EQI(J1,28) 
         ZINT2R(29) = ZINT2R(29) + ZCOEFR * EQI(J1,29) 
         ZINT2R(30) = ZINT2R(30) + ZCOEFR * EQI(J1,30) 
         ZINT2R(31) = ZINT2R(31) + ZCOEFR * EQI(J1,31) 
         ZINT2R(32) = ZINT2R(32) + ZCOEFR * EQI(J1,32) 
C
         ZINT3R( 1) = ZINT3R( 1) + ZCOEFR * EQ3(J1, 1) 
         ZINT3R( 2) = ZINT3R( 2) + ZCOEFR * EQ3(J1, 2) 
         ZINT3R( 3) = ZINT3R( 3) + ZCOEFR * EQ3(J1, 3) 
         ZINT3R( 4) = ZINT3R( 4) + ZCOEFR * EQ3(J1, 4) 
         ZINT3R( 5) = ZINT3R( 5) + ZCOEFR * EQ3(J1, 5) 
         ZINT3R( 6) = ZINT3R( 6) + ZCOEFR * EQ3(J1, 6) 
         ZINT3R( 7) = ZINT3R( 7) + ZCOEFR * EQ3(J1, 7) 
         ZINT3R( 8) = ZINT3R( 8) + ZCOEFR * EQ3(J1, 8) 
         ZINT3R( 9) = ZINT3R( 9) + ZCOEFR * EQ3(J1, 9) 
         ZINT3R(10) = ZINT3R(10) + ZCOEFR * EQ3(J1,10) 
         ZINT3R(11) = ZINT3R(11) + ZCOEFR * EQ3(J1,11) 
         ZINT3R(12) = ZINT3R(12) + ZCOEFR * EQ3(J1,12) 
         ZINT3R(13) = ZINT3R(13) + ZCOEFR * EQ3(J1,13) 
         ZINT3R(14) = ZINT3R(14) + ZCOEFR * EQ3(J1,14) 
         ZINT3R(15) = ZINT3R(15) + ZCOEFR * EQ3(J1,15) 
         ZINT3R(16) = ZINT3R(16) + ZCOEFR * EQ3(J1,16) 
         ZINT3R(17) = ZINT3R(17) + ZCOEFR * EQ3(J1,17) 
C
         ZINTI( 1) = ZINTI( 1) + ZCOEFI * EQL(J1, 1) 
         ZINTI( 2) = ZINTI( 2) + ZCOEFI * EQL(J1, 2) 
         ZINTI( 3) = ZINTI( 3) + ZCOEFI * EQL(J1, 3) 
         ZINTI( 4) = ZINTI( 4) + ZCOEFI * EQL(J1, 4) 
         ZINTI( 5) = ZINTI( 5) + ZCOEFI * EQL(J1, 5) 
         ZINTI( 6) = ZINTI( 6) + ZCOEFI * EQL(J1, 6) 
         ZINTI( 7) = ZINTI( 7) + ZCOEFI * EQL(J1, 7) 
         ZINTI( 8) = ZINTI( 8) + ZCOEFI * EQL(J1, 8) 
         ZINTI( 9) = ZINTI( 9) + ZCOEFI * EQL(J1, 9) 
         ZINTI(10) = ZINTI(10) + ZCOEFI * EQL(J1,10) 
         ZINTI(11) = ZINTI(11) + ZCOEFI * EQL(J1,11) 
         ZINTI(12) = ZINTI(12) + ZCOEFI * EQL(J1,12) 
         ZINTI(13) = ZINTI(13) + ZCOEFI * EQL(J1,13) 
         ZINTI(14) = ZINTI(14) + ZCOEFI * EQL(J1,14) 
         ZINTI(15) = ZINTI(15) + ZCOEFI * EQL(J1,15) 
         ZINTI(16) = ZINTI(16) + ZCOEFI * EQL(J1,16) 
         ZINTI(17) = ZINTI(17) + ZCOEFI * EQL(J1,17) 
         ZINTI(18) = ZINTI(18) + ZCOEFI * EQL(J1,18) 
         ZINTI(19) = ZINTI(19) + ZCOEFI * EQL(J1,19) 
         ZINTI(20) = ZINTI(20) + ZCOEFI * EQL(J1,20) 
         ZINTI(21) = ZINTI(21) + ZCOEFI * EQL(J1,21) 
         ZINTI(22) = ZINTI(22) + ZCOEFI * EQL(J1,22) 
C
         ZINT2I( 1) = ZINT2I( 1) + ZCOEFI * EQI(J1, 1) 
         ZINT2I( 2) = ZINT2I( 2) + ZCOEFI * EQI(J1, 2) 
         ZINT2I( 3) = ZINT2I( 3) + ZCOEFI * EQI(J1, 3) 
         ZINT2I( 4) = ZINT2I( 4) + ZCOEFI * EQI(J1, 4) 
         ZINT2I( 5) = ZINT2I( 5) + ZCOEFI * EQI(J1, 5) 
         ZINT2I( 6) = ZINT2I( 6) + ZCOEFI * EQI(J1, 6) 
         ZINT2I( 7) = ZINT2I( 7) + ZCOEFI * EQI(J1, 7) 
         ZINT2I( 8) = ZINT2I( 8) + ZCOEFI * EQI(J1, 8) 
         ZINT2I( 9) = ZINT2I( 9) + ZCOEFI * EQI(J1, 9) 
         ZINT2I(10) = ZINT2I(10) + ZCOEFI * EQI(J1,10) 
         ZINT2I(11) = ZINT2I(11) + ZCOEFI * EQI(J1,11) 
         ZINT2I(12) = ZINT2I(12) + ZCOEFI * EQI(J1,12) 
         ZINT2I(13) = ZINT2I(13) + ZCOEFI * EQI(J1,13) 
         ZINT2I(14) = ZINT2I(14) + ZCOEFI * EQI(J1,14) 
         ZINT2I(15) = ZINT2I(15) + ZCOEFI * EQI(J1,15) 
         ZINT2I(16) = ZINT2I(16) + ZCOEFI * EQI(J1,16) 
         ZINT2I(17) = ZINT2I(17) + ZCOEFI * EQI(J1,17) 
         ZINT2I(18) = ZINT2I(18) + ZCOEFI * EQI(J1,18) 
         ZINT2I(19) = ZINT2I(19) + ZCOEFI * EQI(J1,19) 
         ZINT2I(20) = ZINT2I(20) + ZCOEFI * EQI(J1,20) 
         ZINT2I(21) = ZINT2I(21) + ZCOEFI * EQI(J1,21) 
         ZINT2I(22) = ZINT2I(22) + ZCOEFI * EQI(J1,22) 
         ZINT2I(23) = ZINT2I(23) + ZCOEFI * EQI(J1,23) 
         ZINT2I(24) = ZINT2I(24) + ZCOEFI * EQI(J1,24) 
         ZINT2I(25) = ZINT2I(25) + ZCOEFI * EQI(J1,25) 
         ZINT2I(26) = ZINT2I(26) + ZCOEFI * EQI(J1,26) 
         ZINT2I(27) = ZINT2I(27) + ZCOEFI * EQI(J1,27) 
         ZINT2I(28) = ZINT2I(28) + ZCOEFI * EQI(J1,28) 
         ZINT2I(29) = ZINT2I(29) + ZCOEFI * EQI(J1,29) 
         ZINT2I(30) = ZINT2I(30) + ZCOEFI * EQI(J1,30) 
         ZINT2I(31) = ZINT2I(31) + ZCOEFI * EQI(J1,31) 
         ZINT2I(32) = ZINT2I(32) + ZCOEFI * EQI(J1,32) 
C
         ZINT3I( 1) = ZINT3I( 1) + ZCOEFI * EQ3(J1, 1) 
         ZINT3I( 2) = ZINT3I( 2) + ZCOEFI * EQ3(J1, 2) 
         ZINT3I( 3) = ZINT3I( 3) + ZCOEFI * EQ3(J1, 3) 
         ZINT3I( 4) = ZINT3I( 4) + ZCOEFI * EQ3(J1, 4) 
         ZINT3I( 5) = ZINT3I( 5) + ZCOEFI * EQ3(J1, 5) 
         ZINT3I( 6) = ZINT3I( 6) + ZCOEFI * EQ3(J1, 6) 
         ZINT3I( 7) = ZINT3I( 7) + ZCOEFI * EQ3(J1, 7) 
         ZINT3I( 8) = ZINT3I( 8) + ZCOEFI * EQ3(J1, 8) 
         ZINT3I( 9) = ZINT3I( 9) + ZCOEFI * EQ3(J1, 9) 
         ZINT3I(10) = ZINT3I(10) + ZCOEFI * EQ3(J1,10) 
         ZINT3I(11) = ZINT3I(11) + ZCOEFI * EQ3(J1,11) 
         ZINT3I(12) = ZINT3I(12) + ZCOEFI * EQ3(J1,12) 
         ZINT3I(13) = ZINT3I(13) + ZCOEFI * EQ3(J1,13) 
         ZINT3I(14) = ZINT3I(14) + ZCOEFI * EQ3(J1,14) 
         ZINT3I(15) = ZINT3I(15) + ZCOEFI * EQ3(J1,15) 
         ZINT3I(16) = ZINT3I(16) + ZCOEFI * EQ3(J1,16) 
         ZINT3I(17) = ZINT3I(17) + ZCOEFI * EQ3(J1,17) 
C
    1    CONTINUE
C
         IF (MOD(KPSI,2) .EQ. 0) THEN
C
           KP = KPSI / 2 + 1
C
           DG11L(KP,KM+1)   = CMPLX(ZINTR( 1),ZINTI( 1))
           DG22L(KP,KM+1)   = CMPLX(ZINTR( 2),ZINTI( 2))
           DG33L(KP,KM+1)   = CMPLX(ZINTR( 3),ZINTI( 3))
           DG12L(KP,KM+1)   = CMPLX(ZINTR( 4),ZINTI( 4))
           JG11L(KP,KM+1)   = CMPLX(ZINTR( 5),ZINTI( 5))
           JG22L(KP,KM+1)   = CMPLX(ZINTR( 6),ZINTI( 6))
           JG33L(KP,KM+1)   = CMPLX(ZINTR( 7),ZINTI( 7))
           JG12L(KP,KM+1)   = CMPLX(ZINTR( 8),ZINTI( 8))
           JACOBI(KP,KM+1)  = CMPLX(ZINTR( 9),ZINTI( 9))
           J2U(KP,KM+1,1)   = CMPLX(ZINTR(10),ZINTI(10))
           J3U(KP,KM+1,1)   = CMPLX(ZINTR(11),ZINTI(11))
           B2E(KP,KM+1,1)   = CMPLX(ZINTR(12),ZINTI(12))
           B3E(KP,KM+1,1)   = CMPLX(ZINTR(13),ZINTI(13))
           PEQ(KP,KM+1,1)   = CMPLX(ZINTR(14),ZINTI(14))
           DPEDS(KP,KM+1,1) = CMPLX(ZINTR(15),ZINTI(15))
           GCHDZ(KP,KM+1)   = CMPLX(ZINTR(16),ZINTI(16))
           GSDZ (KP,KM+1)   = CMPLX(ZINTR(17),ZINTI(17))
           GBZ(KP,KM+1)     = CMPLX(ZINTR(18),ZINTI(18))
           GBR(KP,KM+1)     = CMPLX(ZINTR(19),ZINTI(19))
           JACOBINV(KP,KM+1)= CMPLX(ZINTR(20),ZINTI(20))
           FRM(KP,KM+1)     = CMPLX(ZINTR(21),ZINTI(21))
           FZM(KP,KM+1)     = CMPLX(ZINTR(22),ZINTI(22))
C
           IDIY2(KP,KM+1)   = CMPLX(ZINT2R( 1),ZINT2I( 1))
           IDIY3(KP,KM+1)   = CMPLX(ZINT2R( 2),ZINT2I( 2))
           IG122(KP,KM+1)   = CMPLX(ZINT2R( 3),ZINT2I( 3))
           IG123(KP,KM+1)   = CMPLX(ZINT2R( 4),ZINT2I( 4))
           INXX(KP,KM+1)    = CMPLX(ZINT2R( 5),ZINT2I( 5))
           INXY(KP,KM+1)    = CMPLX(ZINT2R( 6),ZINT2I( 6))
           INYY(KP,KM+1)    = CMPLX(ZINT2R( 7),ZINT2I( 7))
           INZZ(KP,KM+1)    = CMPLX(ZINT2R( 8),ZINT2I( 8))
           IJ0QX(KP,KM+1)   = CMPLX(ZINT2R( 9),ZINT2I( 9))
           IJ0QY(KP,KM+1)   = CMPLX(ZINT2R(10),ZINT2I(10))
           IGPX2(KP,KM+1)   = CMPLX(ZINT2R(11),ZINT2I(11))
           IGPX3(KP,KM+1)   = CMPLX(ZINT2R(12),ZINT2I(12))
           IGPY2(KP,KM+1)   = CMPLX(ZINT2R(13),ZINT2I(13))
           IGPY3(KP,KM+1)   = CMPLX(ZINT2R(14),ZINT2I(14))
           IDRXX(KP,KM+1)   = CMPLX(ZINT2R(15),ZINT2I(15))
           IRXZ(KP,KM+1)    = CMPLX(ZINT2R(16),ZINT2I(16))
           IDRYX(KP,KM+1)   = CMPLX(ZINT2R(17),ZINT2I(17))
           IRYX(KP,KM+1)    = CMPLX(ZINT2R(18),ZINT2I(18))
           IDRZX(KP,KM+1)   = CMPLX(ZINT2R(19),ZINT2I(19))
           IRZY(KP,KM+1)    = CMPLX(ZINT2R(20),ZINT2I(20))
           VISXZ(KP,KM+1)   = CMPLX(ZINT2R(21),ZINT2I(21))
           VISYZ(KP,KM+1)   = CMPLX(ZINT2R(22),ZINT2I(22))
           IVS11(KP,KM+1)   = CMPLX(ZINT2R(23),ZINT2I(23))
           IVS12(KP,KM+1)   = CMPLX(ZINT2R(24),ZINT2I(24))
           IVS21(KP,KM+1)   = CMPLX(ZINT2R(25),ZINT2I(25))
           IVS22(KP,KM+1)   = CMPLX(ZINT2R(26),ZINT2I(26))
           GSFC(KP,KM+1)    = CMPLX(ZINT2R(27),ZINT2I(27))
           GSCC(KP,KM+1)    = CMPLX(ZINT2R(28),ZINT2I(28))
           GSFS(KP,KM+1)    = CMPLX(ZINT2R(29),ZINT2I(29))
           GSCS(KP,KM+1)    = CMPLX(ZINT2R(30),ZINT2I(30))
           GCFC(KP,KM+1)    = CMPLX(ZINT2R(31),ZINT2I(31))
           GCFS(KP,KM+1)    = CMPLX(ZINT2R(32),ZINT2I(32))
C
           EQRHO(KP,KM+1)   = CMPLX(ZINT3R( 1),ZINT3I( 1))
           DRHOS(KP,KM+1)   = CMPLX(ZINT3R( 2),ZINT3I( 2))
           EQROT(KP,KM+1)   = CMPLX(ZINT3R( 3),ZINT3I( 3))
           DROT(KP,KM+1)    = CMPLX(ZINT3R( 4),ZINT3I( 4))
           FEQ(KP,KM+1)     = CMPLX(ZINT3R( 5),ZINT3I( 5))
           IWSQ1(KP,KM+1)   = CMPLX(ZINT3R( 6),ZINT3I( 6))
           IWSQ2(KP,KM+1)   = CMPLX(ZINT3R( 7),ZINT3I( 7))
           IWSQ3(KP,KM+1)   = CMPLX(ZINT3R( 8),ZINT3I( 8))
           IJ0QZ(KP,KM+1)   = CMPLX(ZINT3R( 9),ZINT3I( 9))
           JACOF(KP,KM+1)   = CMPLX(ZINT3R(10),ZINT3I(10))
           B2F(KP,KM+1)     = CMPLX(ZINT3R(11),ZINT3I(11))
           B3F(KP,KM+1)     = CMPLX(ZINT3R(12),ZINT3I(12))
           JACOS(KP,KM+1)   = CMPLX(ZINT3R(13),ZINT3I(13))
           IGF22(KP,KM+1)   = CMPLX(ZINT3R(14),ZINT3I(14))
           B3FC(KP,KM+1)    = CMPLX(ZINT3R(15),ZINT3I(15))
           B2FC(KP,KM+1)    = CMPLX(ZINT3R(16),ZINT3I(16))
           DJCOF(KP,KM+1)   = CMPLX(ZINT3R(17),ZINT3I(17))
C
         ELSE IF (MOD(KPSI,2) .EQ. 1) THEN
C
           KP = (KPSI + 1) / 2
C
           DG11LM(KP,KM+1)   = CMPLX(ZINTR( 1),ZINTI( 1))
           DG22LM(KP,KM+1)   = CMPLX(ZINTR( 2),ZINTI( 2))
           DG33LM(KP,KM+1)   = CMPLX(ZINTR( 3),ZINTI( 3))
           DG12LM(KP,KM+1)   = CMPLX(ZINTR( 4),ZINTI( 4))
           JG11LM(KP,KM+1)   = CMPLX(ZINTR( 5),ZINTI( 5))
           JG22LM(KP,KM+1)   = CMPLX(ZINTR( 6),ZINTI( 6))
           JG33LM(KP,KM+1)   = CMPLX(ZINTR( 7),ZINTI( 7))
           JG12LM(KP,KM+1)   = CMPLX(ZINTR( 8),ZINTI( 8))
           JACOBM(KP,KM+1)   = CMPLX(ZINTR( 9),ZINTI( 9))
           J2E(KP,KM+1,1)    = CMPLX(ZINTR(10),ZINTI(10))
           J3E(KP,KM+1,1)    = CMPLX(ZINTR(11),ZINTI(11))
           B2U(KP,KM+1,1)    = CMPLX(ZINTR(12),ZINTI(12))
           B3U(KP,KM+1,1)    = CMPLX(ZINTR(13),ZINTI(13))
           PRE(KP,KM+1,1)    = CMPLX(ZINTR(14),ZINTI(14))
           DPEDSM(KP,KM+1,1) = CMPLX(ZINTR(15),ZINTI(15))
           GCHDZM(KP,KM+1)   = CMPLX(ZINTR(16),ZINTI(16))
           GSDZM (KP,KM+1)   = CMPLX(ZINTR(17),ZINTI(17))
           GBZM(KP,KM+1)     = CMPLX(ZINTR(18),ZINTI(18))
           GBRM(KP,KM+1)     = CMPLX(ZINTR(19),ZINTI(19))
C
           IDIY2M(KP,KM+1)   = CMPLX(ZINT2R( 1),ZINT2I( 1))
           IDIY3M(KP,KM+1)   = CMPLX(ZINT2R( 2),ZINT2I( 2))
           IG122M(KP,KM+1)   = CMPLX(ZINT2R( 3),ZINT2I( 3))
           IG123M(KP,KM+1)   = CMPLX(ZINT2R( 4),ZINT2I( 4))
           INXXM(KP,KM+1)    = CMPLX(ZINT2R( 5),ZINT2I( 5))
           INXYM(KP,KM+1)    = CMPLX(ZINT2R( 6),ZINT2I( 6))
           INYYM(KP,KM+1)    = CMPLX(ZINT2R( 7),ZINT2I( 7))
           INZZM(KP,KM+1)    = CMPLX(ZINT2R( 8),ZINT2I( 8))
           IJ0QXM(KP,KM+1)   = CMPLX(ZINT2R( 9),ZINT2I( 9))
           IJ0QYM(KP,KM+1)   = CMPLX(ZINT2R(10),ZINT2I(10))
           IGPX2M(KP,KM+1)   = CMPLX(ZINT2R(11),ZINT2I(11))
           IGPX3M(KP,KM+1)   = CMPLX(ZINT2R(12),ZINT2I(12))
           IGPY2M(KP,KM+1)   = CMPLX(ZINT2R(13),ZINT2I(13))
           IGPY3M(KP,KM+1)   = CMPLX(ZINT2R(14),ZINT2I(14))
           IDRXXM(KP,KM+1)   = CMPLX(ZINT2R(15),ZINT2I(15))
           IRXZM(KP,KM+1)    = CMPLX(ZINT2R(16),ZINT2I(16))
           IDRYXM(KP,KM+1)   = CMPLX(ZINT2R(17),ZINT2I(17))
           IRYXM(KP,KM+1)    = CMPLX(ZINT2R(18),ZINT2I(18))
           IDRZXM(KP,KM+1)   = CMPLX(ZINT2R(19),ZINT2I(19))
           IRZYM(KP,KM+1)    = CMPLX(ZINT2R(20),ZINT2I(20))
           VISXZM(KP,KM+1)   = CMPLX(ZINT2R(21),ZINT2I(21))
           VISYZM(KP,KM+1)   = CMPLX(ZINT2R(22),ZINT2I(22))
           IVS11M(KP,KM+1)   = CMPLX(ZINT2R(23),ZINT2I(23))
           IVS12M(KP,KM+1)   = CMPLX(ZINT2R(24),ZINT2I(24))
           IVS21M(KP,KM+1)   = CMPLX(ZINT2R(25),ZINT2I(25))
           IVS22M(KP,KM+1)   = CMPLX(ZINT2R(26),ZINT2I(26))
           GSFCM(KP,KM+1)    = CMPLX(ZINT2R(27),ZINT2I(27))
           GSCCM(KP,KM+1)    = CMPLX(ZINT2R(28),ZINT2I(28))
           GSFSM(KP,KM+1)    = CMPLX(ZINT2R(29),ZINT2I(29))
           GSCSM(KP,KM+1)    = CMPLX(ZINT2R(30),ZINT2I(30))
           GCFCM(KP,KM+1)    = CMPLX(ZINT2R(31),ZINT2I(31))
           GCFSM(KP,KM+1)    = CMPLX(ZINT2R(32),ZINT2I(32))
C
           EQRHOM(KP,KM+1)   = CMPLX(ZINT3R( 1),ZINT3I( 1))
           DRHOSM(KP,KM+1)   = CMPLX(ZINT3R( 2),ZINT3I( 2))
           EQROTM(KP,KM+1)   = CMPLX(ZINT3R( 3),ZINT3I( 3))
           DROTM(KP,KM+1)    = CMPLX(ZINT3R( 4),ZINT3I( 4))
           FEQM(KP,KM+1)     = CMPLX(ZINT3R( 5),ZINT3I( 5))
           IWSQ1M(KP,KM+1)   = CMPLX(ZINT3R( 6),ZINT3I( 6))
           IWSQ2M(KP,KM+1)   = CMPLX(ZINT3R( 7),ZINT3I( 7))
           IWSQ3M(KP,KM+1)   = CMPLX(ZINT3R( 8),ZINT3I( 8))
           IJ0QZM(KP,KM+1)   = CMPLX(ZINT3R( 9),ZINT3I( 9))
           JACOFM(KP,KM+1)   = CMPLX(ZINT3R(10),ZINT3I(10))
           B2FM(KP,KM+1)     = CMPLX(ZINT3R(11),ZINT3I(11))
           B3FM(KP,KM+1)     = CMPLX(ZINT3R(12),ZINT3I(12))
           JACOSM(KP,KM+1)   = CMPLX(ZINT3R(13),ZINT3I(13))
           IGF22M(KP,KM+1)   = CMPLX(ZINT3R(14),ZINT3I(14))
           B3FCM(KP,KM+1)    = CMPLX(ZINT3R(15),ZINT3I(15))
           B2FCM(KP,KM+1)    = CMPLX(ZINT3R(16),ZINT3I(16))
           DJCOFM(KP,KM+1)   = CMPLX(ZINT3R(17),ZINT3I(17))
C
         ENDIF
C
         RETURN
         END
C*DECK C2SM10
C*CALL PROCESS
         SUBROUTINE PROFILE(KN)
C        ######################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM10 TRACE CONSTANT FLUX SURFACES AND EVALUATE PROFILES           *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMMAP.inc'
C 
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C     ZPSITEST IN BALLIT IS A GUESS FOR PSIISO TO CHECK IF
C     PSIISO(1).LT.CPSICL(1). WRITE THESE VALUES SO THAT ONE CAN COMPARE
C     ON THE OUTPUT WITH ZPSITEST WRITTEN BY BALLIT.
C

         WRITE(6,'(/,"IN PROFILE: PSIISO(1)= ",1PE15.8,"  CPSICL(1)= ",
     +     E15.8)') PSIISO(1), CPSICL(1)
C
         IP = ISRCHFGE(KN,PSIISO,1,CPSICL(1))
C
         IF (IP.LT.1)  IP = 1
         IF (IP.GT.KN) IP = KN
C
C     IF IP.NE.1 THEN SIGPSI,.. NOT DEFINED IN ISOFIND WHICH MAY CAUSE
C     A PROBLEM IN SURFACE FOR EXAMPLE
C
         IF (IP .NE. 1) THEN
           WRITE(6,'(//,"********************  WARNING ***************",
     +       /,5X,"IN PROFILE: IP=",I4,/,5X,
     +       "IP SHOULD BE 1, MAY BE R0 TOO FAR FROM RMAG,",
     +       " DO AN EXTRA ITERATION",/,5x,"USE FILE NOUT WITH NOPT=-2",
     +       /,"********************  WARNING ***************"//)') IP
         ENDIF
C
         CALL ISOFIND(IP,KN,SIGPSI,TETPSI,WGTPSI,SPSIM,RC0P)
C
         DO 1 J1=IP,KN
C
         CALL CINT(J1,SIGPSI(1,J1),TETPSI(1,J1),WGTPSI(1,J1))
C
   1     CONTINUE
C
         IF (IP .GT. 1) THEN
C
            IF (NSTTP.LE.2) THEN
C
               ZCID0 = RMAG
               ZCID2 = RMAG**2
C
            ELSE IF (NSTTP .EQ. 3) THEN
C
               ZCID0 = RMAG**2
               ZCID2 = 0.
C
            ENDIF
C
            DO 2 J2=1,IP-1
C
            CID0(J2) = FCCCC0(ZCID0,CID0(IP),CID0(IP+1),CID0(IP+2),
     ,                        SPSIM,PSIISO(IP),PSIISO(IP+1),
     ,                        PSIISO(IP+2),PSIISO(J2))
            CIDR(J2) = FCCCC0(CIDR(IP),CIDR(IP+1),CIDR(IP+2),CIDR(IP+3),
     ,                        PSIISO(IP),PSIISO(IP+1),PSIISO(IP+2),
     ,                        PSIISO(IP+3),PSIISO(J2))
            CIDQ(J2) = FCCCC0(CIDQ(IP),CIDQ(IP+1),CIDQ(IP+2),CIDQ(IP+3),
     ,                        PSIISO(IP),PSIISO(IP+1),PSIISO(IP+2),
     ,                        PSIISO(IP+3),PSIISO(J2))
            CID2(J2) = FCCCC0(ZCID2,CID2(IP),CID2(IP+1),CID2(IP+2),
     ,                        SPSIM,PSIISO(IP),PSIISO(IP+1),
     ,                        PSIISO(IP+2),PSIISO(J2))
C
    2       CONTINUE
C
         ENDIF
C
C        COMPUTE PROFILES
C
         CALL ISOFUN(KN)
C
         RETURN
         END
C*DECK C2SM11
C*CALL PROCESS
         SUBROUTINE BALOON(KP1,KP2,PSM)
C        ################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM11 EVALUATE BALLOONING STABILITY AND LOCAL INTERCHANGE CRITERIA *
*        REDUCED NUMBER OF PSI SURFACES => KP2-KP1+1.LE.NPPSBAL       *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   ICHI0(NPBLC0),      IS0(NPCHI),         IT0(NPCHI),
     R   IC(NPCHI),          PSM(KP2),
     R   ZBND(NPCHI,5),      ZDBDS(NPCHI,16),    ZDBDT(NPCHI,16),
     R   ZDBDST(NPCHI,16),   ZDCHI0(NPCHI),      ZDPBC(NPPSBAL,NPCHI),
     R   ZDPBP(NPPSBAL,NPCHI),ZD2BS2(NPCHI,16),  ZD2BT2(NPCHI,16),   
     R   ZFP(NPPSBAL,NPCHI), ZJAC(NPPSBAL,NPCHI),ZPCEL(NPCHI,16),    
     R   ZR2(NPPSBAL,NPCHI), ZS(NPCHI),          ZS1(NPCHI),         
     R   ZS2(NPCHI),         ZTETA(NPCHI,5),     ZT(NPCHI),          
     R   ZT1(NPCHI),         ZT2(NPCHI)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IMAX  = 2 * NTURN * NCHI + 1
         IM1   = IMAX + 1
         ZEPS  = 1.E-3
         IPSIBAL = KP2 - KP1 + 1
         IF (IPSIBAL .GT. NPPSBAL) THEN
           PRINT *,' ERROR IN BALOON, KP2-KP1+1= ',IPSIBAL,
     +       ' > NPPSBAL= ',NPPSBAL
           STOP 'BALOON'
         ENDIF
C
         IF (NBLC0 .EQ. 1) THEN
C
            ICHI0(1) = 1
C
         ELSE IF (NBLC0 .GT. 1) THEN
C
            IF (NSYM .EQ. 0) THEN
C
               ZDCHI  = CHIM(NCHI) / (NBLC0 - 1.)
               ICHIMX = NCHI
C
            ELSE IF (NSYM .EQ. 1) THEN
C
               ZDCHI  = CHIM(NCHI/2+1) / (NBLC0 - 1.)
               ICHIMX = NCHI / 2 + 1
C
            ENDIF
C
            DO 2 J2=1,NBLC0
C
            DO 1 J1=1,ICHIMX
C
            ZDCHI0(J1) = ABS((J2-1)*ZDCHI-CHIM(J1))
C
    1       CONTINUE
C
            ICHI0(J2) = ISMIN(ICHIMX,ZDCHI0,1)
C
    2       CONTINUE
C
         ENDIF
C
         DO 8 J8=1,IPSIBAL
C
           J8PSI = J8 + KP1 -1
         DO 3 J3=1,NCHI
C
         ZTETA(J3,1) = TETCHI(J3,J8PSI)
         ZTETA(J3,2) = TETCHI(J3,J8PSI) - 2. * ZEPS
         ZTETA(J3,3) = TETCHI(J3,J8PSI) -      ZEPS
         ZTETA(J3,4) = TETCHI(J3,J8PSI) +      ZEPS
         ZTETA(J3,5) = TETCHI(J3,J8PSI) + 2. * ZEPS
C
    3    CONTINUE
C
         CALL BOUND(NCHI,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NCHI,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NCHI,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NCHI,ZTETA(1,4),ZBND(1,4))
         CALL BOUND(NCHI,ZTETA(1,5),ZBND(1,5))
C
         CALL RESETI(IC,NCHI,1)
         DO 4 JT = 1,NT1
            DO 4 JG=1,NCHI
               IF (IC(JG).EQ.0) GOTO 4
               IT0(JG) = JT-1
               IF (TETCHI(JG,J8PSI).LE.CT(JT)) IC(JG)  = 0
 4       CONTINUE
         CALL RESETI(IC,NCHI,1)
         DO 5 JS = 1,NS1
            DO 5 JG=1,NCHI
               IF (IC(JG).EQ.0) GOTO 5
               IS0(JG) = JS-1
               IF (SIGCHI(JG,J8PSI).LE.CSIG(JS)) IC(JG)  = 0
 5       CONTINUE
C
         DO 6 J6=1,NCHI
C
         ZT(J6) = TETCHI(J6,J8PSI)
         ZS(J6) = SIGCHI(J6,J8PSI)
C
         IF (IS0(J6) .GT. NS) IS0(J6) = NS
         IF (IS0(J6) .LT. 1)  IS0(J6) = 1
         IF (IT0(J6) .GT. NT) IT0(J6) = NT
         IF (IT0(J6) .LT. 1)  IT0(J6) = 1
C
         ZS1(J6) = CSIG(IS0(J6))
         ZS2(J6) = CSIG(IS0(J6)+1)
         ZT1(J6) = CT(IT0(J6))
         ZT2(J6) = CT(IT0(J6)+1)
C
 6       CONTINUE
C
         CALL PSICEL(IS0,IT0,NCHI,NPCHI,ZPCEL,CPSICL)
         CALL BASIS3(NCHI,NPCHI,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT,
     &               ZDBDST,ZD2BS2,ZD2BT2)
C
         ZTMF2  = TMF(J8PSI)**2
C
         DO 7 J7=1,NCHI
C
         ZDRSDT = (ZBND(J7,2) + 8*(ZBND(J7,4) - ZBND(J7,3)) -
     -            ZBND(J7,5)) / (12. * ZEPS)
         ZD2RST = (- ZBND(J7,2) + 16. * ZBND(J7,3) -
     -             30. * ZBND(J7,1) + 16. * ZBND(J7,4) -
     -             ZBND(J7,5)) / (12. * ZEPS**2)
C
         ZDPDS = ZDBDS(J7, 1) * ZPCEL(J7, 1) +
     +           ZDBDS(J7, 2) * ZPCEL(J7, 2) +
     +           ZDBDS(J7, 3) * ZPCEL(J7, 3) +
     +           ZDBDS(J7, 4) * ZPCEL(J7, 4) +
     +           ZDBDS(J7, 5) * ZPCEL(J7, 5) +
     +           ZDBDS(J7, 6) * ZPCEL(J7, 6) +
     +           ZDBDS(J7, 7) * ZPCEL(J7, 7) +
     +           ZDBDS(J7, 8) * ZPCEL(J7, 8) +
     +           ZDBDS(J7, 9) * ZPCEL(J7, 9) +
     +           ZDBDS(J7,10) * ZPCEL(J7,10) +
     +           ZDBDS(J7,11) * ZPCEL(J7,11) +
     +           ZDBDS(J7,12) * ZPCEL(J7,12) +
     +           ZDBDS(J7,13) * ZPCEL(J7,13) +
     +           ZDBDS(J7,14) * ZPCEL(J7,14) +
     +           ZDBDS(J7,15) * ZPCEL(J7,15) +
     +           ZDBDS(J7,16) * ZPCEL(J7,16)
C
         ZDPDT = ZDBDT(J7, 1) * ZPCEL(J7, 1) +
     +           ZDBDT(J7, 2) * ZPCEL(J7, 2) +
     +           ZDBDT(J7, 3) * ZPCEL(J7, 3) +
     +           ZDBDT(J7, 4) * ZPCEL(J7, 4) +
     +           ZDBDT(J7, 5) * ZPCEL(J7, 5) +
     +           ZDBDT(J7, 6) * ZPCEL(J7, 6) +
     +           ZDBDT(J7, 7) * ZPCEL(J7, 7) +
     +           ZDBDT(J7, 8) * ZPCEL(J7, 8) +
     +           ZDBDT(J7, 9) * ZPCEL(J7, 9) +
     +           ZDBDT(J7,10) * ZPCEL(J7,10) +
     +           ZDBDT(J7,11) * ZPCEL(J7,11) +
     +           ZDBDT(J7,12) * ZPCEL(J7,12) +
     +           ZDBDT(J7,13) * ZPCEL(J7,13) +
     +           ZDBDT(J7,14) * ZPCEL(J7,14) +
     +           ZDBDT(J7,15) * ZPCEL(J7,15) +
     +           ZDBDT(J7,16) * ZPCEL(J7,16)
C
         ZD2PST = ZDBDST(J7, 1) * ZPCEL(J7, 1) +
     +            ZDBDST(J7, 2) * ZPCEL(J7, 2) +
     +            ZDBDST(J7, 3) * ZPCEL(J7, 3) +
     +            ZDBDST(J7, 4) * ZPCEL(J7, 4) +
     +            ZDBDST(J7, 5) * ZPCEL(J7, 5) +
     +            ZDBDST(J7, 6) * ZPCEL(J7, 6) +
     +            ZDBDST(J7, 7) * ZPCEL(J7, 7) +
     +            ZDBDST(J7, 8) * ZPCEL(J7, 8) +
     +            ZDBDST(J7, 9) * ZPCEL(J7, 9) +
     +            ZDBDST(J7,10) * ZPCEL(J7,10) +
     +            ZDBDST(J7,11) * ZPCEL(J7,11) +
     +            ZDBDST(J7,12) * ZPCEL(J7,12) +
     +            ZDBDST(J7,13) * ZPCEL(J7,13) +
     +            ZDBDST(J7,14) * ZPCEL(J7,14) +
     +            ZDBDST(J7,15) * ZPCEL(J7,15) +
     +            ZDBDST(J7,16) * ZPCEL(J7,16)
C
         ZD2PS2 = ZD2BS2(J7, 1) * ZPCEL(J7, 1) +
     +            ZD2BS2(J7, 2) * ZPCEL(J7, 2) +
     +            ZD2BS2(J7, 3) * ZPCEL(J7, 3) +
     +            ZD2BS2(J7, 4) * ZPCEL(J7, 4) +
     +            ZD2BS2(J7, 5) * ZPCEL(J7, 5) +
     +            ZD2BS2(J7, 6) * ZPCEL(J7, 6) +
     +            ZD2BS2(J7, 7) * ZPCEL(J7, 7) +
     +            ZD2BS2(J7, 8) * ZPCEL(J7, 8) +
     +            ZD2BS2(J7, 9) * ZPCEL(J7, 9) +
     +            ZD2BS2(J7,10) * ZPCEL(J7,10) +
     +            ZD2BS2(J7,11) * ZPCEL(J7,11) +
     +            ZD2BS2(J7,12) * ZPCEL(J7,12) +
     +            ZD2BS2(J7,13) * ZPCEL(J7,13) +
     +            ZD2BS2(J7,14) * ZPCEL(J7,14) +
     +            ZD2BS2(J7,15) * ZPCEL(J7,15) +
     +            ZD2BS2(J7,16) * ZPCEL(J7,16)
C
         ZD2PT2 = ZD2BT2(J7, 1) * ZPCEL(J7, 1) +
     +            ZD2BT2(J7, 2) * ZPCEL(J7, 2) +
     +            ZD2BT2(J7, 3) * ZPCEL(J7, 3) +
     +            ZD2BT2(J7, 4) * ZPCEL(J7, 4) +
     +            ZD2BT2(J7, 5) * ZPCEL(J7, 5) +
     +            ZD2BT2(J7, 6) * ZPCEL(J7, 6) +
     +            ZD2BT2(J7, 7) * ZPCEL(J7, 7) +
     +            ZD2BT2(J7, 8) * ZPCEL(J7, 8) +
     +            ZD2BT2(J7, 9) * ZPCEL(J7, 9) +
     +            ZD2BT2(J7,10) * ZPCEL(J7,10) +
     +            ZD2BT2(J7,11) * ZPCEL(J7,11) +
     +            ZD2BT2(J7,12) * ZPCEL(J7,12) +
     +            ZD2BT2(J7,13) * ZPCEL(J7,13) +
     +            ZD2BT2(J7,14) * ZPCEL(J7,14) +
     +            ZD2BT2(J7,15) * ZPCEL(J7,15) +
     +            ZD2BT2(J7,16) * ZPCEL(J7,16)
C
         ZFP(J8,J7) = (ZDPDS**2 + (ZDPDT / SIGCHI(J7,J8PSI) -
     -                 ZDPDS * ZDRSDT / ZBND(J7,1))**2) /
     /                ZBND(J7,1)**2
         ZGRADP = SQRT(ZFP(J8,J7))
C
         ZCOST = COS(ZTETA(J7,1))
         ZSINT = SIN(ZTETA(J7,1))
C
         ZRHO        = SIGCHI(J7,J8PSI) * ZBND(J7,1)
         ZR          = ZRHO * ZCOST + R0
         ZR2(J8,J7)  = ZR**2
         ZJAC(J8,J7) = CP(J8PSI) * ZR**NER * ZGRADP**NEGP
C
         ZDSDR = (ZDRSDT * ZSINT + ZBND(J7,1) * ZCOST) / ZBND(J7,1)**2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBND(J7,1) * ZSINT - ZDRSDT * ZCOST) / ZBND(J7,1)**2
         ZDTDZ = ZCOST / ZRHO
C
         ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         Z1     = ZDPDS
         ZDZ1DS = ZD2PS2
         ZDZ1DT = ZD2PST
         Z2     = ZDPDT / SIGCHI(J7,J8PSI) - ZDRSDT * ZDPDS / ZBND(J7,1)
         ZDZ2DS = (ZD2PST-ZDPDT / SIGCHI(J7,J8PSI)) / SIGCHI(J7,J8PSI) -
     -             ZDRSDT * ZD2PS2 / ZBND(J7,1)
         ZDZ2DT = - ZD2PST*ZDRSDT/ZBND(J7,1) + ZD2PT2/SIGCHI(J7,J8PSI) +
     +            ZDPDS * (ZDRSDT**2-ZD2RST*ZBND(J7,1)) / ZBND(J7,1)**2
C
         ZDFDS = 2 * (Z1 * ZDZ1DS + Z2 * ZDZ2DS) / ZBND(J7,1)**2
         ZDFDT = - 2 * ZFP(J8,J7) * ZDRSDT / ZBND(J7,1) +
     +           2 * (Z1 * ZDZ1DT + Z2 * ZDZ2DT) / ZBND(J7,1)**2
C
         ZDFDR = ZDFDS * ZDSDR + ZDFDT * ZDTDR
         ZDFDZ = ZDFDS * ZDSDZ + ZDFDT * ZDTDZ
C
         ZDGDPN = (ZDFDR * ZDPDR + ZDFDZ * ZDPDZ) / ZFP(J8,J7)
         ZDRDPN = ZDPDR / ZFP(J8,J7)
C
         ZDRDCP = - ZJAC(J8,J7) * ZDPDZ / ZR
         ZDGDCP = ZJAC(J8,J7) * (ZDPDR * ZDFDZ - ZDPDZ * ZDFDR) / ZR
         ZDGDPC = ZDGDPN - .5*ZDGDCP *EQ13(J7,J8PSI)/(PSM(J8PSI)*CPSRF)
         ZDRDPC = ZDRDPN - .5*ZDRDCP *EQ13(J7,J8PSI)/(PSM(J8PSI)*CPSRF)
C
         ZDPBP(J8,J7) = CPPR(J8PSI) +(TTP(J8PSI)+.5*ZDGDPN)/ZR2(J8,J7) -
     -                  (ZTMF2+ZFP(J8,J7)) * ZDRDPN / (ZR * ZR2(J8,J7))
         ZDPBC(J8,J7) = .5 * ZDGDCP / ZR2(J8,J7) - (ZTMF2+ZFP(J8,J7)) *
     *                  ZDRDCP / (ZR * ZR2(J8,J7))
C
 7       CONTINUE
 8       CONTINUE
C
         CALL RESETI(NCBAL(KP1),IPSIBAL,0)
         DO J=1,NBLC0
           DO I=KP1,KP2
             NCBLNS(I,J) = 0
           ENDDO
         ENDDO
C
         DO 14 J14=1,NBLC0
C
         CALL VZERO(ABAL,2*NPPSBAL*IM1)
C
         JC = ICHI0(J14)
C
         DO 10 J10=1,IMAX
C
         DO 9 J9=1,IPSIBAL
C
         J9PSI = J9 + KP1 - 1
         ZNU0 = ZJAC(J9,JC) * TMF(J9PSI) / ZR2(J9,JC)
C
         JTOT  = J10 + JC - 2
         JTURN = JTOT / NCHI
         JCHI  = JTOT - JTURN * NCHI + 1
C
C EQ. (18) IN PUBLICATION
C
         ZGBAR = .5 * (ZNU0 * EQ13(JC,J9PSI) + QPSI(J9PSI) * 
     *          (EQ22(JCHI,J9PSI)-EQ22(JC,J9PSI)))/(PSM(J9PSI)*CPSRF)+
     +           CDQ(J9PSI) * (EQ24(JCHI,J9PSI)-EQ24(JC,J9PSI) + 
     +           CPI*(2*(JTURN - NTURN)))
         ZCOEF = ZR2(J9,JCHI) / (TMF(J9PSI)**2 + ZFP(J9,JCHI))
C
C  COEFFICIENT OF D(XSI) / D(CHI) (C_1, EQ. (19) IN PUBLICATION)
C
         ZF = (1. + ZCOEF * (ZGBAR * ZFP(J9,JCHI))**2) / 
     /        (ZJAC(J9,JCHI)**2 * ZFP(J9,JCHI))
C
C  COEFFICIENT OF XSI (C_2, EQ. (17) IN PUBLICATION)
C
         ZG = - 2. * CPPR(J9PSI) * ZCOEF * (ZDPBP(J9,JCHI) - ZCOEF * 
     *          TMF(J9PSI) * ZGBAR * ZDPBC(J9,JCHI) / ZJAC(J9,JCHI))
C
C  COMPUTATTION OF ABAL
C
         ZDY = CHI(JCHI+1) - CHI(JCHI)
C
         ZZA  = ZJAC(J9,JCHI) * ZF
         ZZP  = ZJAC(J9,JCHI) * ZG
         ZX11 = .25 * ZZP + ZZA / (ZDY * ZDY)
         ZX12 = .25 * ZZP - ZZA / (ZDY * ZDY)
         ZX22 = .25 * ZZP + ZZA / (ZDY * ZDY)
C
         ABAL(J9,1,J10)   = ABAL(J9,1,J10)   + ZX11 * ZDY
         ABAL(J9,2,J10)   = ABAL(J9,2,J10)   + ZX12 * ZDY
         ABAL(J9,1,J10+1) = ABAL(J9,1,J10+1) + ZX22 * ZDY
C
 9       CONTINUE
   10    CONTINUE
C
C BOUNDARY CONDITIONS
C
         DO 11 J11=1,IPSIBAL
C
         ABAL(J11,1,1)    = 1.
         ABAL(J11,2,1)    = 0.
         ABAL(J11,2,IMAX) = 0.
         ABAL(J11,1,IM1)  = 1.
         ABAL(J11,2,IM1)  = 0.
C
   11    CONTINUE
C
C FIND THE NUMBER OF NEGATIVE EIGENVALUES
C
         CALL NTRIDG(ABAL,NPPSBAL,1,IPSIBAL,IM1)
C
         DO 13 J13=1,IM1
C
         DO 12 J12=1,IPSIBAL
C
         IF (ABAL(J12,1,J13) .LT. 0.) THEN
C
            J12PSI = J12 + KP1 - 1
            NCBAL(J12PSI)      = NCBAL(J12PSI) + 1
            NCBLNS(J12PSI,J14) = NCBLNS(J12PSI,J14) + 1
C
         ENDIF
C
  12     CONTINUE
  13     CONTINUE
  14     CONTINUE
C
         DO 15 J15=1,NBLC0
C
         CHI0(J15) = CHIM(ICHI0(J15))
C
  15     CONTINUE
C
         RETURN
         END
C*DECK C2SM12
C*CALL PROCESS
         SUBROUTINE GLOQUA(PSM,PS,KP,KCASE)
C        ##################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM12 GLOBAL EQUILIBRIUM QUANTITIES IN TABLE 1 OF PUBLICATION      *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   PS(KP),        PSM(KP),
     R   ZBET(NPISO),   ZJPSI(NPISO),
     R   ZPISOS(NPISO), ZPISOM(NPISO),
     R   ZBND(12*NPT),  ZTET(12*NPT),
     R   ZR(12*NPT)
C     
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         INCLUDE 'CUCCCC.inc'
         INCLUDE 'QUAQQQ.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IBND = 12 * NT
C
         ZDT  = 2.* CPI / DFLOAT(IBND-1)
C
         DO 1 J1=1,IBND
C
         ZTET(J1) = (J1 - 1.) * ZDT
C
   1     CONTINUE
C
         CALL BOUND(IBND,ZTET,ZBND)
C
         DO 2 J2=1,IBND
C
         ZR(J2) = R0  + ZBND(J2) * COS(ZTET(J2))
C
   2     CONTINUE
C
         IMN = ISMIN(IBND,ZR,1)
         IMX = ISMAX(IBND,ZR,1)
C
         ZRMJ = .5 * (ZR(IMN) + ZR(IMX))
C
         DO 3 J3=1,KP
C
         ZPISOS(J3) = CPSRF * PS(J3)**2
         ZPISOM(J3) = CPSRF * PSM(J3)**2
C
    3    CONTINUE
C
         RINDUC(1) = 0.
         AREA      = 0.
         VSURF(1)  = 0.
         ZIP2      = 0.
         ZIP       = 0.
         ZIB2      = 0.
         ZVOL      = 0.
         ZJPSI(1)  = 0.
         RITBS     = 0.
         RITBSC    = 0.
C
         DO 4 J4=2,KP
C
         AREA       = AREA + RARE(J4-1) * (ZPISOS(J4) - ZPISOS(J4-1))
         ZIP2       = ZIP2 + RIP2(J4-1) * (ZPISOS(J4) - ZPISOS(J4-1))
         ZIP        = ZIP  + RIP(J4-1) * (ZPISOS(J4) - ZPISOS(J4-1))
         ZIB2       = ZIB2 + RIB2(J4-1) * (ZPISOS(J4) - ZPISOS(J4-1))
         VSURF(J4)  = VSURF(J4-1)  + RIVOL(J4-1) * 
     &                               (ZPISOS(J4) - ZPISOS(J4-1))
         ZJPSI(J4)  = ZJPSI(J4-1)  + RIIE(J4-1) * 
     &                               (ZPISOS(J4) - ZPISOS(J4-1))
         RINDUC(J4) = RINDUC(J4-1) + RIIR(J4-1) * 
     &                               (ZPISOS(J4) - ZPISOS(J4-1))
C
c%OS         IF (RJPAR(J4-1) .NE. 0.)
c%OS     &   RITBS      = RITBS + RJBSH(J4-1) * RIIE(J4-1) / RJPAR(J4-1) * 
c%OS     &                                     (ZPISOS(J4) - ZPISOS(J4-1))
         IF (RJPAR(J4-1) .NE. 0.)
     &     RITBS  = RITBS + RJBSOS(J4-1,1) * RIIE(J4-1) / RJPAR(J4-1) *
     &                                     (ZPISOS(J4) - ZPISOS(J4-1))
         IF (RJPAR(J4-1) .NE. 0.)
     &     RITBSC = RITBSC + RJBSOS(J4-1,2) * RIIE(J4-1) / RJPAR(J4-1) *
     &                                     (ZPISOS(J4) - ZPISOS(J4-1))
C
    4    CONTINUE
C
         VOLUME  = 2. * CPI * VSURF(KP)
         CPBAR   = ZIP / VSURF(KP)
         RITOT   = ZJPSI(KP)
         RINOR   = RITOT / (ASPCT * TMF(KP))
         BETA    = 2. * ZIP / ZIB2
         BETAP   = 8. * CPI * ZIP / (ZRMJ * RITOT**2)
         BETAS   = 2. * SQRT(ZIP2 * VSURF(KP)) /  ZIB2
         BETAX   = 2. * CPBAR * (ZRMJ / TMF(KP))**2
         CONVF   = .5 * RLENG(KP)**2 * ZRMJ / VOLUME
         RIBSNOR = RITBS / (ASPCT * TMF(KP))
C
         IF (CPBAR .EQ. 0.) THEN
C
            CPPF   = 0.
C
         ELSE
C
            CPPF   = CP0 / CPBAR
C
         ENDIF
C
         DO 10 J10=1,KP
C
C     1) INTERPOLATE J5 TO OBTAIN J5' REQUIRED FOR LOCAL INTERCHANGE
C     STABILITY CRITERIA.     
C
c         I1 = J10
c         IF (j10.eq.1)  I1 = 2
c         IF (j10.eq.KP) I1 = KP-1
C
c         RJ5P(J10) = FQQQ1(RJ5(I1-1),RJ5(I1),RJ5(I1+1),
c     ,                     PSM(I1-1),PSM(I1),PSM(I1+1),
c     ,                     PSM(J10))
c         RJ5P(J10) = RJ5P(J10) / (2.*PSM(J10)*CPSRF)
c
c         CDQ(J10) = FQQQ1(QPSI(I1-1),QPSI(I1),QPSI(I1+1),
c     ,                    PSM(I1-1),PSM(I1),PSM(I1+1),
c     ,                    PSM(J10))
c         CDQ(J10) = CDQ(J10) / (2.*PSM(J10)*CPSRF)
C     
C     2) COMPUTE MERCIER PARAMETER -D_I (EQ.(19) IN CHEASE PAPER)
C     
         SMERCI(J10) = (CPPR(J10)*TMF(J10)*RJ2(J10)/CDQ(J10) - .5)**2 +
     &                  CPPR(J10) * (RJ5P(J10) - CPPR(J10) * RJ3(J10)) * 
     &                 (TMF(J10)**2 * RJ1(J10) + RJ4(J10)) / CDQ(J10)**2
C     
C     3) COMPUTE H OF GREENE, GLASSER, JOHNSON (EQ.(21) IN CHEASE PAPER)
C     
         HMERCR(J10) = TMF(J10) * CPPR(J10) / CDQ(J10) * 
     &                 (RJ2(J10) - RJ5(J10) * (RJ4(J10) + TMF(J10)**2 *
     &                  RJ1(J10)) / (RJ6(J10) + TMF(J10)**2 * RJ4(J10)))
C
C     4) COMPUTE RESISTIVE INTERCHANGE PARAMETER -D_R 
C        (EQ.(20) IN CHEASE PAPER)
C     
         SMERCR(J10) = SMERCI(J10) - (HMERCR(J10) - .5)**2
C
   10    CONTINUE
C
         IF (KCASE .EQ. 2) RETURN
C
         ZBET(1) = 0.
C
         DO 5 J5=2,KP
C
         I1 = J5 - 1
         I2 = J5
         IF (I1.LE.3) I1 = 3
         IF (I2.LE.3) I2 = 3
         IF (I1.GE.KP-1) I1 = KP-1
         IF (I2.GE.KP-1) I2 = KP-1
C
         ZCPPR1 = FCCCC0(CPPR(I1-2),CPPR(I1-1),CPPR(I1),CPPR(I1+1),
     ,                   PSM(I1-2),PSM(I1-1),PSM(I1),PSM(I1+1),
     ,                   PS(J5-1))
         ZCPPR2 = FCCCC0(CPPR(I2-2),CPPR(I2-1),CPPR(I2),CPPR(I2+1),
     ,                   PSM(I2-2),PSM(I2-1),PSM(I2),PSM(I2+1),
     ,                   PS(J5))
C
         ZBET(J5) = ZBET(J5-1) + .5 * (ZPISOS(J5) - ZPISOS(J5-1)) *
     *              (ZCPPR1 * VSURF(J5-1) + ZCPPR2 * VSURF(J5))
C
         BETAB(J5)  = - 8. * CPI * ZBET(J5) / (ZJPSI(J5)**2 * ZRMJ)
         RINDUC(J5) = 4. * CPI * RINDUC(J5) / (ZRMJ * ZJPSI(J5)**2)
         DPRIME(J5) = .5 * RINDUC(J5) + BETAB(J5)
C
    5    CONTINUE
C
         BETAB(1)  = FCCCC0(BETAB(2),BETAB(3),BETAB(4),BETAB(5),
     ,                      PS(2),PS(3),PS(4),PS(5),RC0P)
         RINDUC(1) = FCCCC0(RINDUC(2),RINDUC(3),RINDUC(4),RINDUC(5),
     ,                      PS(2),PS(3),PS(4),PS(5),RC0P)
         DPRIME(1) = .5 * RINDUC(1) + BETAB(1)
C
         DO 6 J6=1,KP
C
         RSURF(J6)  = SQRT(VSURF(J6) / VSURF(KP))
C
   6     CONTINUE
C
         CDRQ(1)  = 0.
         RDEDR(1) = 2. * (RELL(2) - RELL(1)) / (RSURF(3) - RSURF(1))
C
         DO 8 J8=2,KP-1
C
         CDRQ(J8)  = 4 * RSURF(J8) * (QPSI(J8) - QPSI(J8-1)) /
     /               ((QPSI(J8) + QPSI(J8-1)) *
     *               (RSURF(J8+1) - RSURF(J8-1)))
         RDEDR(J8) = 2. * (RELL(J8) - RELL(J8-1)) / 
     /               (RSURF(J8+1) - RSURF(J8-1))
C
    8    CONTINUE
C
         CDRQ(KP)  = 2 * RSURF(KP) * (QPSI(KP) - QPSI(KP - 1)) /
     /                  ((RSURF(KP) - RSURF(KP - 1)) * QPSI(KP))
         RDEDR(KP) = 2. * (RELL(KP) - RELL(KP - 1)) / 
     /                  (RSURF(KP) - RSURF(KP - 1))
C     
         DO 9 J9=2,KP
C
         IF (J9 .EQ. KP) THEN
C
            ZPPR = (CPR(KP)-CPR(KP - 1)) / (RSURF(KP)-RSURF(KP - 1))
C
         ELSE
C
            ZPPR = (CPR(J9+1)-CPR(J9-1)) / (RSURF(J9+1)-RSURF(J9-1))
C
         ENDIF
C
         RDI(J9) = .25 + 2 * RSURF(J9) * ZPPR / (T0*CDRQ(J9))**2 *
     *              (1. - QPSI(J9)**2 * (1. + 1.5*DPRIME(J9)*RDEDR(J9) -
     -               .75 * (2*RELL(J9) + RSURF(J9) * RDEDR(J9))))
         RSY(J9) = .25 + 2 * RSURF(J9) * ZPPR / (T0*CDRQ(J9))**2 *
     *              (1. - QPSI(J9)**2)
C
    9    CONTINUE
C
         RDI(1) = FCCCC0(RDI(2),RDI(3),RDI(4),RDI(5),
     ,                   PS(2),PS(3),PS(4),PS(5),RC0P)
         RSY(1) = FCCCC0(RSY(2),RSY(3),RSY(4),RSY(5),
     ,                   PS(2),PS(3),PS(4),PS(5),RC0P)
C
         IF (NRFP .EQ. 0) RETURN
C
         ZBPHI = 0.
C
         DO 11 J11=2,KP
C
         ZBPHI     = ZBPHI + 4. * CPI * PSM(J11-1) * CPSRF * 
     *               QPSI(J11-1) * (PS(J11) - PS(J11-1))
         RFPBP(J11) = ZJPSI(J11) / RLENG(J11)
         write(*,101) J11,PS(J11),RLENG(J11),ZJPSI(J11)
C
  11     CONTINUE
  101    format(I3,3(E12.4))
C
         ZBPHI = ZBPHI / AREA
C
         RFPF = TMF(KP) * RLENG1(KP) / (RLENG(KP) * ZBPHI)
         RFPT = RFPBP(KP) / ZBPHI
C
         RETURN
         END
C*DECK C2SM13
C*CALL PROCESS
         SUBROUTINE VACUUM(KM)
C        #####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM13 COMPUTE VACUUM EQ'S FOR MARS (SEE SECTION 5.4.2 AND TABLE 2  *
*        IN PUBLICATION)                                              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
C
         DIMENSION
     R   ZBND(NPMGS*NTP1,4), ZDCDT(NPMGS*NTP1),    
     R   ZDR1DC(NPMGS*NTP1), ZDZ1DC(NPMGS*NTP1),    
     R   ZR1(NPMGS*NTP1),    ZZ1(NPMGS*NTP1),
     R   ZS(2*NPV1),         ZTETA(NPMGS*NTP1,4)
C
         COMPLEX     ZCG11,ZCG22,ZCG33,ZCG12,ZCGRM,ZCGZM,ZCDGIJ,ZIW
         COMPLEX     ZIWT
         DIMENSION
     C   ZCG11(2*NPV1),      ZCG22(2*NPV1),     ZCG33(2*NPV1),
     C   ZCG12(2*NPV1),      ZCGRM(2*NPV1),     ZCGZM(2*NPV1),
     C   ZCDGIJ(6)
C        DIMENSION           ZCGJA(2*NPV1)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IEDGE = 2*NPSI
         ZEPS = 1.E-3
C
         DO 1 J1=1,NV
C
         ZS(2*(J1-1)+1) = CSV(J1)
         ZS(2*J1      ) = CSMV(J1)
C
    1    CONTINUE
C
         ZS(2*NV+1) = CSV(NV1)
C
         DO 2 J2=1,NMGAUS*NT1
C
         ZTETA(J2,1) = TETPSI(J2,IEDGE) - 2. * ZEPS
         ZTETA(J2,2) = TETPSI(J2,IEDGE) -      ZEPS
         ZTETA(J2,3) = TETPSI(J2,IEDGE) +      ZEPS
         ZTETA(J2,4) = TETPSI(J2,IEDGE) + 2. * ZEPS
C
    2    CONTINUE
C
         CALL BOUND(NMGAUS*NT1,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,4),ZBND(1,4))
C
         DO 3 J3=1,NMGAUS*NT1
C
         ZCOST = COS(TETPSI(J3,IEDGE))
         ZSINT = SIN(TETPSI(J3,IEDGE))
         ZBNDT = BNDISO(J3,IEDGE)
         ZRHO  = ZBNDT
         ZR    = RRISO(J3,IEDGE)
         ZZ    = RZISO(J3,IEDGE)
         ZGP   = GPISO(J3,IEDGE)
         ZDPDS = DPSISO(J3,IEDGE)
         ZDPDT = 0.
C
         ZDRSDT = (ZBND(J3,1) + 8*(ZBND(J3,3) - ZBND(J3,2)) -
     -             ZBND(J3,4)) / (12. * ZEPS)
         ZJAC   = CP(IEDGE) * ZR**NER * ZGP**NEGP
C
         ZDCDT(J3) = ZRHO * ZBNDT * ZR / (ZJAC * ZDPDS)
         ZR1(J3)   = ZR - R0W
         ZZ1(J3)   = ZZ - RZ0W
C
         ZDSDR = (ZDRSDT * ZSINT + ZBNDT * ZCOST) / ZBNDT**2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBNDT * ZSINT - ZDRSDT * ZCOST) / ZBNDT**2
         ZDTDZ = ZCOST / ZRHO
C
         ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         ZDR1DC(J3) = - ZJAC * ZDPDZ / ZR
         ZDZ1DC(J3) = ZJAC * ZDPDR / ZR
C
    3    CONTINUE
C
         CALL CVZERO(ZCG11,2*NV1)
         CALL CVZERO(ZCG22,2*NV1)
         CALL CVZERO(ZCG33,2*NV1)
         CALL CVZERO(ZCG12,2*NV1)
         CALL CVZERO(ZCGRM,2*NV1)
         CALL CVZERO(ZCGZM,2*NV1)
C        CALL CVZERO(ZCGJA,2*NV1)
C
         DO 9 J9=1,2*NV+1
C
         DO 8 J8=1,NT1
C
         CALL CVZERO(ZCDGIJ,6)
C
         DO 7 J7=1,NMGAUS
C
         IG = (J8-1)*NMGAUS+J7
C
         ZR    = ZR1(IG)
         ZZ    = ZZ1(IG)
         ZS1   = ZS(J9)
         ZDRDC = ZDR1DC(IG)
         ZDZDC = ZDZ1DC(IG)
         ZPWGT = WGTPSI(IG,IEDGE)
         ZARG  = DFLOAT(KM) * CHIISO(IG)
         ZIW   = (0.5/CPI * ZDCDT(IG) * ZPWGT)
     *         *  CMPLX(COS(ZARG),-SIN(ZARG))
         ZARGT  = DFLOAT(KM) * TETPSI(IG,IEDGE)
         ZIWT  = (0.5/CPI * ZPWGT) *  CMPLX(COS(ZARGT),-SIN(ZARGT))
C
         G11L =           ZR**2 + ZZ**2
         G22L = ZS1**2 * (ZDRDC**2 + ZDZDC**2)
         G12L =   ZS1  * (ZR * ZDRDC + ZZ * ZDZDC)
         G33L =          (R0W + ZS1 * ZR)**2
         ZJAC1= ZS1 * (R0W + ZS1 * ZR) * (ZR * ZDZDC - ZZ * ZDRDC)
         RML  = R0W + ZS1*ZR
         ZML  = RZ0W + ZS1*ZZ
C
         ZCDGIJ(1) = ZCDGIJ(1) + ZIW * G11L / ZJAC1
         ZCDGIJ(2) = ZCDGIJ(2) + ZIW * G22L / ZJAC1
         ZCDGIJ(3) = ZCDGIJ(3) + ZIW * G33L / ZJAC1
         ZCDGIJ(4) = ZCDGIJ(4) + ZIW * G12L / ZJAC1
         ZCDGIJ(5) = ZCDGIJ(5) + ZIW * RML
         ZCDGIJ(6) = ZCDGIJ(6) + ZIW * ZML
C
    7    CONTINUE
C
         ZCG11(J9) = ZCG11(J9) + ZCDGIJ(1)
         ZCG22(J9) = ZCG22(J9) + ZCDGIJ(2)
         ZCG33(J9) = ZCG33(J9) + ZCDGIJ(3)
         ZCG12(J9) = ZCG12(J9) + ZCDGIJ(4)
         ZCGRM(J9) = ZCGRM(J9) + ZCDGIJ(5)
         ZCGZM(J9) = ZCGZM(J9) + ZCDGIJ(6)
C
    8    CONTINUE
    9    CONTINUE
C
         CALL CCOPY(NV1,ZCG11(1),2,DG11LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG11(2),2,DG11LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCG22(1),2,DG22LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG22(2),2,DG22LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCG33(1),2,DG33LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG33(2),2,DG33LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCG12(1),2,DG12LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG12(2),2,DG12LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCGRM(1),2,DGRMLV(1,KM+1),1)
         CALL CCOPY(NV1,ZCGZM(1),2,DGZMLV(1,KM+1),1)
C
         RETURN
         END
C*DECK C2SM13
C*CALL PROCESS
         SUBROUTINE VACUUMNW(KM)
C        #####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* NONCONFORMAL WALLS                                                  *
* LIU YQ, DECEMBER 12, 2002                                           *                                                                    
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
         INCLUDE 'COMBND.inc'
C
         DIMENSION
     R   ZBND(NPMGS*NTP1,4), ZDCDT(NPMGS*NTP1,NWBPS0), 
     R   RRBND(NPMGS*NTP1,5,NWBPS0),RZBND(NPMGS*NTP1,5,NWBPS0),
     R   RCBND(NPMGS*NTP1,5,NWBPS0),
     R   ZDR1DT(NPMGS*NTP1,NWBPS0), ZDZ1DT(NPMGS*NTP1,NWBPS0),    
     R   ZS(2*NPV1),         ZTETA(NPMGS*NTP1,5),
     R   ZL(NWBPS0)
C
         COMPLEX     ZCG11,ZCG22,ZCG33,ZCG12,ZCGRM,ZCGZM,ZCDGIJ,ZIW
         COMPLEX     ZIWT
         DIMENSION
     C   ZCG11(2*NPV1),      ZCG22(2*NPV1),     ZCG33(2*NPV1),
     C   ZCG12(2*NPV1),      ZCGRM(2*NPV1),     ZCGZM(2*NPV1),
     C   ZCDGIJ(6)
C        DIMENSION           ZCGJA(2*NPV1)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IEDGE = 2*NPSI
         ZEPS = 1.E-3
C
         DO 1 J1=1,NV
C
         ZS(2*(J1-1)+1) = CSV(J1)
         ZS(2*J1      ) = CSMV(J1)
C
    1    CONTINUE
C
         ZS(2*NV+1) = CSV(NV1)
C
         DO 2 J2=1,NMGAUS*NT1
C
         ZTETA(J2,1) = TETPSI(J2,IEDGE) - 2. * ZEPS
         ZTETA(J2,2) = TETPSI(J2,IEDGE) -      ZEPS
         ZTETA(J2,3) = TETPSI(J2,IEDGE) +      ZEPS
         ZTETA(J2,4) = TETPSI(J2,IEDGE) + 2. * ZEPS
         ZTETA(J2,5) = TETPSI(J2,IEDGE)
C
    2    CONTINUE
C
         CALL BOUND(NMGAUS*NT1,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,4),ZBND(1,4))

         DO JNW = 1,NWBPS
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,1),RRBND(1,1,JNW),
     *                   RZBND(1,1,JNW),RCBND(1,1,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,2),RRBND(1,2,JNW),
     *                   RZBND(1,2,JNW),RCBND(1,2,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,3),RRBND(1,3,JNW),
     *                   RZBND(1,3,JNW),RCBND(1,3,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,4),RRBND(1,4,JNW),
     *                   RZBND(1,4,JNW),RCBND(1,4,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,5),RRBND(1,5,JNW),
     *                   RZBND(1,5,JNW),RCBND(1,5,JNW))
         ENDDO

         J1 = NWBPS-1
         J2 = NWBPS
         ZSA = ZS(2*NW(J1)-1)
         ZSB = ZS(2*NW(J2)-1)
         DO J3=1,NMGAUS*NT1
            DO J4=1,5
               RRBND(J3,J4,J2) = R0W + ZSB/ZSA*(RRBND(J3,J4,J1)-R0W)
               RZBND(J3,J4,J2) = RZ0W + ZSB/ZSA*(RZBND(J3,J4,J1)-RZ0W)
            ENDDO
         ENDDO

C
         CALL VZERO(ZL,NWBPS)

         DO 3 J3=1,NMGAUS*NT1
C
         ZCOST = COS(TETPSI(J3,IEDGE))
         ZSINT = SIN(TETPSI(J3,IEDGE))
         ZBNDT = BNDISO(J3,IEDGE)
         ZRHO  = ZBNDT
         ZR    = RRISO(J3,IEDGE)
         ZZ    = RZISO(J3,IEDGE)
         ZGP   = GPISO(J3,IEDGE)
         ZDPDS = DPSISO(J3,IEDGE)
         ZDPDT = 0.
C
         ZDRSDT = (ZBND(J3,1) + 8*(ZBND(J3,3) - ZBND(J3,2)) -
     -             ZBND(J3,4)) / (12. * ZEPS)
         ZJAC   = CP(IEDGE) * ZR**NER * ZGP**NEGP
C
         ZDCDT(J3,1) = ZRHO * ZBNDT * ZR / (ZJAC * ZDPDS)

         DO JNW = 1,NWBPS
            ZDRDT = (RRBND(J3,1,JNW) + 8*(RRBND(J3,3,JNW) - 
     -               RRBND(J3,2,JNW)) - RRBND(J3,4,JNW)) / (12. * ZEPS)
            ZDZDT = (RZBND(J3,1,JNW) + 8*(RZBND(J3,3,JNW) - 
     -               RZBND(J3,2,JNW)) - RZBND(J3,4,JNW)) / (12. * ZEPS)
            ZDR1DT(J3,JNW) = ZDRDT
            ZDZ1DT(J3,JNW) = ZDZDT

C     REDEFINE ZDCDT FOR EACH SURFACE
            IF (JNW.GE.2) THEN
               ZDCDT(J3,JNW) = SQRT(ZDRDT*ZDRDT + ZDZDT*ZDZDT)
               ZL(JNW) = ZL(JNW) + ZDCDT(J3,JNW)*WGTPSI(J3,IEDGE)
            ENDIF
         ENDDO
    3    CONTINUE

         DO JNW = 2,NWBPS
C            WRITE(*,*) 'ZL(JNW)=',JNW,ZL(JNW)
            DO J3 = 1,NMGAUS*NT1
               ZDCDT(J3,JNW) = ZDCDT(J3,JNW)*2.0*CPI/ZL(JNW)
            ENDDO
         ENDDO

C        IF (KM.EQ.1) THEN
C           WRITE(*,*) 'RRBND(K)  RZBND(K)'
C           DO J3=1,NMGAUS*NT1
C        WRITE(*,1014) (RRBND(J3,5,JNW),RZBND(J3,5,JNW),JNW=1,NWBPS)
C           ENDDO
C
C           WRITE(*,*) 'TETPSI  ZDCDT(K)'
C        DO J3=1,NMGAUS*NT1
C              WRITE(*,1014) ZTETA(J3,5),(ZDCDT(J3,JNW),JNW=1,NWBPS)
C           ENDDO
C        ENDIF
C1014    FORMAT(8(1X,E12.5))

C
         CALL CVZERO(ZCG11,2*NV1)
         CALL CVZERO(ZCG22,2*NV1)
         CALL CVZERO(ZCG33,2*NV1)
         CALL CVZERO(ZCG12,2*NV1)
         CALL CVZERO(ZCGRM,2*NV1)
         CALL CVZERO(ZCGZM,2*NV1)
C        CALL CVZERO(ZCGJA,2*NV1)
C
         DO 10 J1=1,NWBPS-1
         ZSA = ZS(2*NW(J1)-1)
         ZSB = ZS(2*NW(J1+1)-1)
         ZSBA2 = (ZSB-ZSA)**2
         DO 9 J9=2*NW(J1)-1,2*NW(J1+1)-1
         ZS1 = ZS(J9)
         ZL0 = 0.0
         DO J3 = 1,NMGAUS*NT1
            ZDRDT = ZDR1DT(J3,J1) + (ZDR1DT(J3,J1+1)-ZDR1DT(J3,J1))*
     *              (ZS1-ZSA)/(ZSB-ZSA)
            ZDZDT = ZDZ1DT(J3,J1) + (ZDZ1DT(J3,J1+1)-ZDZ1DT(J3,J1))*
     *              (ZS1-ZSA)/(ZSB-ZSA)
            ZL0 = ZL0 + SQRT(ZDRDT*ZDRDT+ZDZDT*ZDZDT)*WGTPSI(J3,IEDGE)
         ENDDO
C
         DO 8 J8=1,NT1
C
         CALL CVZERO(ZCDGIJ,6)
C
         DO 7 J7=1,NMGAUS
C
         IG = (J8-1)*NMGAUS+J7
C
         ZR    = RRBND(IG,5,J1+1) - RRBND(IG,5,J1)
         ZZ    = RZBND(IG,5,J1+1) - RZBND(IG,5,J1)
         ZDRDCA = ZDR1DT(IG,J1)/ZDCDT(IG,1)
         ZDRDCB = ZDR1DT(IG,J1+1)/ZDCDT(IG,1)
         ZDZDCA = ZDZ1DT(IG,J1)/ZDCDT(IG,1)
         ZDZDCB = ZDZ1DT(IG,J1+1)/ZDCDT(IG,1)
         ZPWGT = WGTPSI(IG,IEDGE)
         ZARG  = DFLOAT(KM) * CHIISO(IG)

         ZDRDT = ZDR1DT(IG,J1) + (ZDR1DT(IG,J1+1)-ZDR1DT(IG,J1))*
     *           (ZS1-ZSA)/(ZSB-ZSA)
         ZDZDT = ZDZ1DT(IG,J1) + (ZDZ1DT(IG,J1+1)-ZDZ1DT(IG,J1))*
     *           (ZS1-ZSA)/(ZSB-ZSA)
         ZDCDT0 = SQRT(ZDRDT*ZDRDT + ZDZDT*ZDZDT)*2.0*CPI/ZL0

         ZIW   = (0.5/CPI * ZDCDT(IG,1) * ZPWGT)
     *         *  CMPLX(COS(ZARG),-SIN(ZARG))
         ZARGT  = DFLOAT(KM) * TETPSI(IG,IEDGE)
         ZIWT  = (0.5/CPI * ZPWGT) *  CMPLX(COS(ZARGT),-SIN(ZARGT))
C
         G11L = (ZR**2 + ZZ**2)/ZSBA2

         G22L = (ZDRDCB*(ZS1-ZSA) + ZDRDCA*(ZSB-ZS1))**2 
         G22L = G22L + (ZDZDCB*(ZS1-ZSA) + ZDZDCA*(ZSB-ZS1))**2 
         G22L = G22L/ZSBA2
        
         G12L =  ZR*((ZS1-ZSA)*ZDRDCB+(ZSB-ZS1)*ZDRDCA) + 
     +           ZZ*((ZS1-ZSA)*ZDZDCB+(ZSB-ZS1)*ZDZDCA)
         G12L = G12L/ZSBA2

         G33L =  (RRBND(IG,5,J1+1)*(ZS1-ZSA) + 
     +            RRBND(IG,5,J1)*(ZSB-ZS1))**2/ZSBA2

         RML  = RRBND(IG,5,J1+1)*(ZS1-ZSA) + RRBND(IG,5,J1)*(ZSB-ZS1)
         RML = RML/(ZSB-ZSA)

         ZML  = RZBND(IG,5,J1+1)*(ZS1-ZSA) + RZBND(IG,5,J1)*(ZSB-ZS1)
         ZML = ZML/(ZSB-ZSA)

         ZJAC1 = ZR*((ZS1-ZSA)*ZDZDCB+(ZSB-ZS1)*ZDZDCA) -
     -           ZZ*((ZS1-ZSA)*ZDRDCB+(ZSB-ZS1)*ZDRDCA)
         ZJAC1 = ZJAC1*(RRBND(IG,5,J1+1)*(ZS1-ZSA) + 
     +            RRBND(IG,5,J1)*(ZSB-ZS1))/ZSBA2/(ZSB-ZSA)
C
         ZCDGIJ(1) = ZCDGIJ(1) + ZIW * G11L / ZJAC1
         ZCDGIJ(2) = ZCDGIJ(2) + ZIW * G22L / ZJAC1
         ZCDGIJ(3) = ZCDGIJ(3) + ZIW * G33L / ZJAC1
         ZCDGIJ(4) = ZCDGIJ(4) + ZIW * G12L / ZJAC1
         ZCDGIJ(5) = ZCDGIJ(5) + ZIW * RML
         ZCDGIJ(6) = ZCDGIJ(6) + ZIW * ZML
C
    7    CONTINUE
C
         ZCG11(J9) = ZCG11(J9) + ZCDGIJ(1)
         ZCG22(J9) = ZCG22(J9) + ZCDGIJ(2)
         ZCG33(J9) = ZCG33(J9) + ZCDGIJ(3)
         ZCG12(J9) = ZCG12(J9) + ZCDGIJ(4)
         ZCGRM(J9) = ZCGRM(J9) + ZCDGIJ(5)
         ZCGZM(J9) = ZCGZM(J9) + ZCDGIJ(6)
C
    8    CONTINUE
    9    CONTINUE
   10    CONTINUE
C
C TAKE MAIN VALUE OF THE FUNCTION CROSS THE WALLS
         DO J1=2,NWBPS-1
            J9 = 2*NW(J1)-1
            ZCG11(J9) = ZCG11(J9)*0.5
            ZCG22(J9) = ZCG22(J9)*0.5
            ZCG33(J9) = ZCG33(J9)*0.5
            ZCG12(J9) = ZCG12(J9)*0.5
            ZCGRM(J9) = ZCGRM(J9)*0.5
            ZCGZM(J9) = ZCGZM(J9)*0.5
         ENDDO            

         CALL CCOPY(NV1,ZCG11(1),2,DG11LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG11(2),2,DG11LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCG22(1),2,DG22LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG22(2),2,DG22LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCG33(1),2,DG33LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG33(2),2,DG33LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCG12(1),2,DG12LV(1,KM+1),1)
         CALL CCOPY(NV,ZCG12(2),2,DG12LMV(1,KM+1),1)
         CALL CCOPY(NV1,ZCGRM(1),2,DGRMLV(1,KM+1),1)
         CALL CCOPY(NV1,ZCGZM(1),2,DGZMLV(1,KM+1),1)
C
         RETURN
         END
C
         SUBROUTINE VACUUMRNW
C        ####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* FOR NONCONFORMAL WALLS                                              *
* LIU YQ, AUG. 09, 2003                                               *
*                                                                     *
* NOTE: TETPSI=GAUSSIAN POINTS FOR PHYSICAL ANGLE THETA               *
*       CHIISO=GAUSSIAN POINTS FOR NUMBERICAL ANGLE CHI               * 
*       ZDCDT=DCHIISO/DTETPSI                                         *
*       THEREFORE, DZ/DCHISIO=(DZ/DTETPSI)/(ZDCDT)                    *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMRCH.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
         INCLUDE 'COMBND.inc'
C
      INTEGER NPISOC
      PARAMETER (NPISOC = NPMGS*NTP1)

         DIMENSION
     R   ZBND(NPMGS*NTP1,4), ZDCDT(NPMGS*NTP1,NWBPS0), 
     R   RRBND(NPMGS*NTP1,5,NWBPS0),RZBND(NPMGS*NTP1,5,NWBPS0),
     R   RCBND(NPMGS*NTP1,5,NWBPS0),
     R   ZDR1DT(NPMGS*NTP1,NWBPS0), ZDZ1DT(NPMGS*NTP1,NWBPS0),    
     R   ZS(2*NPV1),         ZTETA(NPMGS*NTP1,5),
     R   ZL(NWBPS0),         CHIISON(NPMGS*NTP1),
     R   ZDCDS(NPMGS*NTP1)
      DIMENSION
     R     ZG11LI(NPISOC),      ZG22LI(NPISOC),
     R     ZG33LI(NPISOC),      ZG12LI(NPISOC),
     R     ZRJAI(NPISOC),       ZCNDI(NPISOC)
         DIMENSION
     R     ICHIISO(2*NPCHI1), 
     R     ZCHI(NP2CHI), ZD2FUN(NPISOC), ZWORK(NPISOC,3),
     R     ZA(NP2CHI), ZB(NP2CHI), ZC(NP2CHI), ZD(NP2CHI),
     R     ZG11L(NP2CHI),      ZG22L(NP2CHI),
     R     ZG33L(NP2CHI),      ZG12L(NP2CHI),
     R     ZRJA(NP2CHI),       ZCND(NP2CHI)
CMSC5/14/04        THE FOLLOWING ADDED TO GET COORDINATES IN VACUUM   
         DIMENSION
     R     ZRANI(NPISOC),       ZZANI(NPISOC),
     R     ZRAN(NP2CHI),        ZZAN(NP2CHI)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
CMSC
         WRITE(*,'(" VACUUMNRW, NWBPS=",I4)') NWBPS
         IEDGE = 2*NPSI
         ZEPS = 1.E-3
C
C
C     PREPARE COEFFICIENTS FOR FFT AND CUBIC SPLINE
C        
         INCHI = 2*NCHI
C
         ZDCHI = 2. * CPI / DFLOAT(INCHI)
         DO I=1,INCHI+1
           ZCHI(I) = DFLOAT(I-1)*ZDCHI
         ENDDO
C
         IGCHISO = NMGAUS*NT1
C     ASSUMES THAT CHIISO HAS ALREADY ARRAY AT BOUNDARY (FROM CALL FOURFFT)
C     (OTHERWISE, ADD CALL GCHI(2*NPSI))
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHI(I)) - 1
C
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
           ZA(I) = (CHIISO(ICHISO+1) - ZCHI(I)) / ZH
           ZB(I) = (ZCHI(I) - CHIISO(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISO(ICHISO+1) - ZCHI(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHI(I) - CHIISO(ICHISO)) / 6.
C
         ENDDO
C
C     PRECOMPUTE SOME ARRAYS
C
         DO 1 J1=1,NV
C
         ZS(2*(J1-1)+1) = CSV(J1)
         ZS(2*J1      ) = CSMV(J1)
C
    1    CONTINUE
C
         ZS(2*NV+1) = CSV(NV1)
C
         DO 2 J2=1,NMGAUS*NT1
C
         ZTETA(J2,1) = TETPSI(J2,IEDGE) - 2. * ZEPS
         ZTETA(J2,2) = TETPSI(J2,IEDGE) -      ZEPS
         ZTETA(J2,3) = TETPSI(J2,IEDGE) +      ZEPS
         ZTETA(J2,4) = TETPSI(J2,IEDGE) + 2. * ZEPS
         ZTETA(J2,5) = TETPSI(J2,IEDGE)
C
    2    CONTINUE
C
         CALL BOUND(NMGAUS*NT1,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,4),ZBND(1,4))

         DO JNW = 1,NWBPS
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,1),RRBND(1,1,JNW),
     *                   RZBND(1,1,JNW),RCBND(1,1,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,2),RRBND(1,2,JNW),
     *                   RZBND(1,2,JNW),RCBND(1,2,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,3),RRBND(1,3,JNW),
     *                   RZBND(1,3,JNW),RCBND(1,3,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,4),RRBND(1,4,JNW),
     *                   RZBND(1,4,JNW),RCBND(1,4,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,5),RRBND(1,5,JNW),
     *                   RZBND(1,5,JNW),RCBND(1,5,JNW))
         ENDDO

         J1 = NWBPS-1
         J2 = NWBPS
         ZSA = ZS(2*NW(J1)-1)
         ZSB = ZS(2*NW(J2)-1)
         DO J3=1,NMGAUS*NT1
            DO J4=1,5
               RRBND(J3,J4,J2) = R0W + ZSB/ZSA*(RRBND(J3,J4,J1)-R0W)
               RZBND(J3,J4,J2) = RZ0W + ZSB/ZSA*(RZBND(J3,J4,J1)-RZ0W)
            ENDDO
         ENDDO

C
         CALL VZERO(ZL,NWBPS)

         DO 3 J3=1,NMGAUS*NT1
C
         ZCOST = COS(TETPSI(J3,IEDGE))
         ZSINT = SIN(TETPSI(J3,IEDGE))
         ZBNDT = BNDISO(J3,IEDGE)
         ZRHO  = ZBNDT
         ZR    = RRISO(J3,IEDGE)
         ZZ    = RZISO(J3,IEDGE)
         ZGP   = GPISO(J3,IEDGE)
         ZDPDS = DPSISO(J3,IEDGE)
         ZDPDT = 0.
C
         ZDRSDT = (ZBND(J3,1) + 8*(ZBND(J3,3) - ZBND(J3,2)) -
     -             ZBND(J3,4)) / (12. * ZEPS)
         ZJAC   = CP(IEDGE) * ZR**NER * ZGP**NEGP
C
         ZDCDT(J3,1) = ZRHO * ZBNDT * ZR / (ZJAC * ZDPDS)

         DO JNW = 1,NWBPS
            ZDRDT = (RRBND(J3,1,JNW) + 8*(RRBND(J3,3,JNW) - 
     -               RRBND(J3,2,JNW)) - RRBND(J3,4,JNW)) / (12. * ZEPS)
            ZDZDT = (RZBND(J3,1,JNW) + 8*(RZBND(J3,3,JNW) - 
     -               RZBND(J3,2,JNW)) - RZBND(J3,4,JNW)) / (12. * ZEPS)
            ZDR1DT(J3,JNW) = ZDRDT
            ZDZ1DT(J3,JNW) = ZDZDT

C     REDEFINE ZDCDT FOR EACH SURFACE
            IF (JNW.GE.2) THEN
               ZDCDT(J3,JNW) = SQRT(ZDRDT*ZDRDT + ZDZDT*ZDZDT)
               ZL(JNW) = ZL(JNW) + ZDCDT(J3,JNW)*WGTPSI(J3,IEDGE)
            ENDIF
         ENDDO
    3    CONTINUE

         DO JNW = 2,NWBPS
            DO J3 = 1,NMGAUS*NT1
               ZDCDT(J3,JNW) = ZDCDT(J3,JNW)*2.0*CPI/ZL(JNW)
            ENDDO
         ENDDO

         DO 10 J1=1,NWBPS-1
         ZSA = ZS(2*NW(J1)-1)
         ZSB = ZS(2*NW(J1+1)-1)
         ZSBA2 = (ZSB-ZSA)**2
         NWUP = 2*NW(J1+1)-2
         IF (J1.EQ.(NWBPS-1)) NWUP = NWUP + 1
         DO 9 J9=2*NW(J1)-1,NWUP
         ZS1 = ZS(J9)
         ZL0 = 0.0
         DO J3 = 1,NMGAUS*NT1
            ZDRDT = ZDR1DT(J3,J1) + (ZDR1DT(J3,J1+1)-ZDR1DT(J3,J1))*
     *              (ZS1-ZSA)/(ZSB-ZSA)
            ZDZDT = ZDZ1DT(J3,J1) + (ZDZ1DT(J3,J1+1)-ZDZ1DT(J3,J1))*
     *              (ZS1-ZSA)/(ZSB-ZSA)
            ZL0 = ZL0 + SQRT(ZDRDT*ZDRDT+ZDZDT*ZDZDT)*WGTPSI(J3,IEDGE)
         ENDDO
C

         J99 = (J9+1)/2
         IF (J99*2-1.NE.J9) J99=0
         JW = 0
         DO L=1,NWBPS
            IF (ABS(NW(L)-J99)<1.0E-1) JW=J99
         ENDDO
         DO L=1,5
            IF (ABS(RCOIL(L)-ZS1)<1.0E-8) JW=J99
         ENDDO

C        MODIFY CHIISO INSIDE VACUUM TO DENSIFY CHI-MESH IN LFS
         CALL BOUNDDCDT(NMGAUS*NT1,ZS1,ZTETA(1,5),CHIISO,CHIISON,
     &                  ZDCDT(1,1),ZDCDS)
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISON,1,ZCHI(I)) - 1
 
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISON(ICHISO+1) - CHIISON(ICHISO)
           ZA(I) = (CHIISON(ICHISO+1) - ZCHI(I)) / ZH
           ZB(I) = (ZCHI(I) - CHIISON(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISON(ICHISO+1) - ZCHI(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHI(I) - CHIISON(ICHISO)) / 6.

         ENDDO

C        TEST
C        IF (J9.EQ.80) THEN
C        WRITE(*,*) 'TEST CHIISON'
C        DO J8=1,NT1*NMGAUS
C           WRITE(*,1001) ZTETA(J8,5),CHIISO(J8),CHIISON(J8),ZDCDS(J8)
C        ENDDO
C1001    FORMAT(4E16.8) 
C        ENDIF
     
         DO 8 J8=1,NT1*NMGAUS
         IG = J8

         ZR    = RRBND(IG,5,J1+1) - RRBND(IG,5,J1)
         ZZ    = RZBND(IG,5,J1+1) - RZBND(IG,5,J1)
         ZRN   = RRBND(IG,5,J1) + ZR*(ZS1-ZSA)/(ZSB-ZSA)
         ZZN   = RZBND(IG,5,J1) + ZZ*(ZS1-ZSA)/(ZSB-ZSA)
         ZDRDT = ZDR1DT(IG,J1) + (ZDR1DT(IG,J1+1)-ZDR1DT(IG,J1))*
     *           (ZS1-ZSA)/(ZSB-ZSA)
         ZDZDT = ZDZ1DT(IG,J1) + (ZDZ1DT(IG,J1+1)-ZDZ1DT(IG,J1))*
     *           (ZS1-ZSA)/(ZSB-ZSA)
         ZDRDC = ZDRDT/ZDCDT(IG,1)         
         ZDZDC = ZDZDT/ZDCDT(IG,1)         
         ZDRDS = ZR/(ZSB-ZSA) - ZDRDC*ZDCDS(IG)
         ZDZDS = ZZ/(ZSB-ZSA) - ZDZDC*ZDCDS(IG)

         ZZG11L = ZDRDS**2 + ZDZDS**2
         ZZG12L = ZDRDS*ZDRDC + ZDZDS*ZDZDC
         ZZG22L = ZDRDC**2 + ZDZDC**2
         ZZG33L = ZRN**2
         ZJAC1  = ZRN*(ZDRDS*ZDZDC - ZDRDC*ZDZDS)

         ZRANI (J8) = ZRN
         ZZANI (J8) = ZZN
         ZRJAI (J8) = ZJAC1
         ZG11LI(J8) = ZZG11L
         ZG22LI(J8) = ZZG22L
         ZG33LI(J8) = ZZG33L
         ZG12LI(J8) = ZZG12L
    8    CONTINUE
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C     3. FOR EACH ARRAY: COMPUTE VALUES ON ZCHI USING A PERIODIC
C     CUBIC SPLINE FIT AND COMPUTE FULL FOURIER TRANSFORM
C
C     EQL
C
CMSC5/14/04  ADDED LINES
         CALL SPLCHI(CHIISON,ZRANI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRAN,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISON,ZZANI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZZAN,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISON,ZRJAI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRJA,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISON,ZG11LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG11L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISON,ZG22LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG22L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISON,ZG12LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG12L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISON,ZG33LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG33L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
C-----------------------------------------------------------------------
C
      IF (J9 .GT. 1) GOTO 20
C
         CLOSE(NRMAR)
         OPEN(NRMAR,FILE='OUTVMAR',STATUS='NEW',FORM='FORMATTED')
         REWIND NRMAR
C
         WRITE(NRMAR,*) NV1, INCHI, NWBPS 
 20   CONTINUE 
C
      WRITE(NRMAR,*) ZS1
      WRITE(NRMAR,*) (ZRJA(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG11L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG22L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG33L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZG12L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZRAN(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZZAN(J),J=1,INCHI)
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
    9    CONTINUE
   10    CONTINUE

C     WRITE POLOIDAL ANGLES
C     WRITE(*,*) ' ZTETA           CHIISO'   
C     DO J=1,NT1*NMGAUS
C        WRITE(*,2004) ZTETA(J,5),CHIISO(J)
C     ENDDO

C     WRITE CONDUCTIVITY COEFFICIENTS FOR ALL RESISTIVE WALLS
      DO J1=2,NWBPS-1
         DO J=1,NT1*NMGAUS
            ZCNDI(J) = RCBND(J,5,J1)
         ENDDO
         CALL SPLCHI(CHIISON,ZCNDI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZCND,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
         WRITE(NRMAR,*) (ZCND(J),J=1,INCHI) 
         DO J=1,NT1*NMGAUS
            ZCNDI(J) = (RCBND(J,1,J1) + 8*(RCBND(J,3,J1) - 
     -               RCBND(J,2,J1)) - RCBND(J,4,J1)) / (12. * ZEPS)
            ZCNDI(J) = ZCNDI(J)/ZDCDT(J,1)
         ENDDO
         CALL SPLCHI(CHIISON,ZCNDI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZCND,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
         WRITE(NRMAR,*) (ZCND(J),J=1,INCHI) 
      ENDDO

      CLOSE(NRMAR)

 2004 FORMAT (3E16.8)

      RETURN
      END

         SUBROUTINE VACUUMRNW2
C        #####################
***********************************************************************
*                                                                     *
* FOR NONCONFORMAL WALLS                                              *
* AND SECOND ORDER VACUUM MESH TO OBTAIN CONTINUOUS METRICS           *
* LIU YQ, FEB. 07, 2012                                               *
*                                                                     *
* NOTE: TETPSI=GAUSSIAN POINTS FOR PHYSICAL ANGLE THETA               *
*       CHIISO=GAUSSIAN POINTS FOR NUMBERICAL ANGLE CHI               * 
*       ZDCDT=DCHIISO/DTETPSI                                         *
*       THEREFORE, DZ/DCHISIO=(DZ/DTETPSI)/(ZDCDT)                    *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMRCH.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
         INCLUDE 'COMBND.inc'
C
      INTEGER NPISOC
      PARAMETER (NPISOC = NPMGS*NTP1)

         DIMENSION
     R   ZDCDT(NPMGS*NTP1,NWBPS0),  CHIISON(NPMGS*NTP1),
     R   RRBND(NPMGS*NTP1,5,NWBPS0),RZBND(NPMGS*NTP1,5,NWBPS0),
     R   RCBND(NPMGS*NTP1,5,NWBPS0),
     R   ZDR1DT(NPMGS*NTP1,NWBPS0), ZDZ1DT(NPMGS*NTP1,NWBPS0),    
     R   ZS(2*NPV1),         ZTETA(NPMGS*NTP1,5),
     R   ZL(NWBPS0),         ZDCDS(NPMGS*NTP1)
         DIMENSION
     R     ZG11LI(NPISOC),      ZG22LI(NPISOC),
     R     ZG33LI(NPISOC),      ZG12LI(NPISOC),
     R     ZRJAI(NPISOC),       ZCNDI(NPISOC)
         DIMENSION
     R     Z2RA(NPISOC), Z2RAS(NPISOC), Z2GR(NPISOC),
     R     Z2ZA(NPISOC), Z2ZAS(NPISOC), Z2GZ(NPISOC),
     R     Z2RAT(NPISOC), Z2RAST(NPISOC), Z2GRT(NPISOC),
     R     Z2ZAT(NPISOC), Z2ZAST(NPISOC), Z2GZT(NPISOC)
         DIMENSION
     R     ICHIISO(2*NPCHI1), 
     R     ZCHI(NP2CHI), ZD2FUN(NPISOC), ZWORK(NPISOC,3),
     R     ZA(NP2CHI), ZB(NP2CHI), ZC(NP2CHI), ZD(NP2CHI),
     R     ZG11L(NP2CHI),      ZG22L(NP2CHI),
     R     ZG33L(NP2CHI),      ZG12L(NP2CHI),
     R     ZRJA(NP2CHI),       ZCND(NP2CHI)
         DIMENSION
     R     ZRANI(NPISOC),       ZZANI(NPISOC),
     R     ZRAN(NP2CHI),        ZZAN(NP2CHI)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
CMSC
         WRITE(*,'(" VACUUMNRW2, NWBPS=",I4)') NWBPS
         IEDGE = 2*NPSI
         ZEPS  = 1.E-3
         Z2GP  = VZ2GP
 
C        PREPARE COEFFICIENTS FOR FFT AND CUBIC SPLINE
         
         INCHI = 2*NCHI
 
         ZDCHI = 2. * CPI / DFLOAT(INCHI)
         DO I=1,INCHI+1
           ZCHI(I) = DFLOAT(I-1)*ZDCHI
         ENDDO
 
         IGCHISO = NMGAUS*NT1
C        ASSUMES THAT CHIISO HAS ALREADY ARRAY AT BOUNDARY (FROM CALL FOURFFT)
C        (OTHERWISE, ADD CALL GCHI(2*NPSI))
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHI(I)) - 1
 
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
           ZA(I) = (CHIISO(ICHISO+1) - ZCHI(I)) / ZH
           ZB(I) = (ZCHI(I) - CHIISO(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISO(ICHISO+1) - ZCHI(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHI(I) - CHIISO(ICHISO)) / 6.

         ENDDO

C        PRECOMPUTE SOME ARRAYS

         DO 1 J1=1,NV

         ZS(2*(J1-1)+1) = CSV(J1)
         ZS(2*J1      ) = CSMV(J1)

    1    CONTINUE

         ZS(2*NV+1) = CSV(NV1)

         DO 2 J2=1,NMGAUS*NT1

         ZTETA(J2,1) = TETPSI(J2,IEDGE) - 2. * ZEPS
         ZTETA(J2,2) = TETPSI(J2,IEDGE) -      ZEPS
         ZTETA(J2,3) = TETPSI(J2,IEDGE) +      ZEPS
         ZTETA(J2,4) = TETPSI(J2,IEDGE) + 2. * ZEPS
         ZTETA(J2,5) = TETPSI(J2,IEDGE)

    2    CONTINUE

         DO JNW = 1,NWBPS
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,1),RRBND(1,1,JNW),
     *                   RZBND(1,1,JNW),RCBND(1,1,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,2),RRBND(1,2,JNW),
     *                   RZBND(1,2,JNW),RCBND(1,2,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,3),RRBND(1,3,JNW),
     *                   RZBND(1,3,JNW),RCBND(1,3,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,4),RRBND(1,4,JNW),
     *                   RZBND(1,4,JNW),RCBND(1,4,JNW))
            CALL BOUNDNW(NMGAUS*NT1,JNW,ZTETA(1,5),RRBND(1,5,JNW),
     *                   RZBND(1,5,JNW),RCBND(1,5,JNW))
         ENDDO

         J1 = NWBPS-1
         J2 = NWBPS
         ZSA = ZS(2*NW(J1)-1)
         ZSB = ZS(2*NW(J2)-1)
         DO J3=1,NMGAUS*NT1
            DO J4=1,5
               RRBND(J3,J4,J2) = R0W + ZSB/ZSA*(RRBND(J3,J4,J1)-R0W)
               RZBND(J3,J4,J2) = RZ0W + ZSB/ZSA*(RZBND(J3,J4,J1)-RZ0W)
            ENDDO
         ENDDO

C
         CALL VZERO(ZL,NWBPS)

         DO 3 J3=1,NMGAUS*NT1
C
         ZBNDT = BNDISO(J3,IEDGE)
         ZRHO  = ZBNDT
         ZR    = RRISO(J3,IEDGE)
         ZZ    = RZISO(J3,IEDGE)
         ZGP   = GPISO(J3,IEDGE)
         ZDPDS = DPSISO(J3,IEDGE)
         ZDPDT = 0.

         ZJAC   = CP(IEDGE) * ZR**NER * ZGP**NEGP
         ZDCDT(J3,1) = ZRHO * ZBNDT * ZR / (ZJAC * ZDPDS)

         DO JNW = 1,NWBPS
            ZDRDT = (RRBND(J3,1,JNW) + 8*(RRBND(J3,3,JNW) - 
     -               RRBND(J3,2,JNW)) - RRBND(J3,4,JNW)) / (12. * ZEPS)
            ZDZDT = (RZBND(J3,1,JNW) + 8*(RZBND(J3,3,JNW) - 
     -               RZBND(J3,2,JNW)) - RZBND(J3,4,JNW)) / (12. * ZEPS)
            ZDR1DT(J3,JNW) = ZDRDT
            ZDZ1DT(J3,JNW) = ZDZDT

C           REDEFINE ZDCDT FOR EACH SURFACE
            IF (JNW.GE.2) THEN
               ZDCDT(J3,JNW) = SQRT(ZDRDT*ZDRDT + ZDZDT*ZDZDT)
               ZL(JNW) = ZL(JNW) + ZDCDT(J3,JNW)*WGTPSI(J3,IEDGE)
            ENDIF
         ENDDO
    3    CONTINUE

         DO JNW = 2,NWBPS
            DO J3 = 1,NMGAUS*NT1
               ZDCDT(J3,JNW) = ZDCDT(J3,JNW)*2.0*CPI/ZL(JNW)
            ENDDO
         ENDDO

C        COMPUTE Z2* QUANTITIES AT PLASMA SURFACE (I.E. J1=1)
         ZSA  = ZS(2*NW(1)-1)
         ZSB  = ZS(2*NW(1+1)-1)
         ZSBA = ZSB-ZSA
         DO IG=1,NT1*NMGAUS
            Z2RA(IG) = RRBND(IG,5,1)
            Z2ZA(IG) = RZBND(IG,5,1)

            ZDPSIS = 2.*CPSRF
            ZR = RRISO(IG,IEDGE)
            ZJAC   = ZDPSIS * CP(IEDGE) * ZR**NER *
     *               GPISO(IG,IEDGE)**NEGP
            ZGRADS = GPISO(IG,IEDGE) / ZDPSIS
            ZBSCHI = BCHISO(IG)
            ZGCHI2 = (ZBSCHI * ZGRADS)**2 + (ZR / (ZJAC * ZGRADS))**2
            ZG11LA = (ZJAC*ZBSCHI*ZGRADS/ZR)**2 + 1./(ZGRADS**2)
            ZG12LA = - ZBSCHI * (ZJAC*ZGRADS / ZR)**2

            ZDPDR = DPRISO(IG,IEDGE)
            ZDPDZ = DPZISO(IG,IEDGE)
            ZDSDR = ZDPDR/ZDPSIS
            ZDSDZ = ZDPDZ/ZDPSIS
            ZDCDZ = ZDSDZ*ZBSCHI+ZDSDR*ZR/(ZJAC*ZGRADS**2)

            Z2RAS(IG) = ZDCDZ*ZJAC/ZR
            Z2ZAS(IG) = ZG11LA*ZDSDZ + ZG12LA*ZDCDZ

            Z2GR(IG) =(RRBND(IG,5,2)-Z2RA(IG)-Z2RAS(IG)*ZSBA)/ZSBA**Z2GP
            Z2GZ(IG) =(RZBND(IG,5,2)-Z2ZA(IG)-Z2ZAS(IG)*ZSBA)/ZSBA**Z2GP
         ENDDO

         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2RA,Z2RAT)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2ZA,Z2ZAT)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2RAS,Z2RAST)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2ZAS,Z2ZAST)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2GR,Z2GRT)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2GZ,Z2GZT)      

C        CHECK 
C        WRITE(*,*) 'VACUUMRNW2 CHECK:',NMGAUS,NT1          
C        DO J=1,NMGAUS*NT1
C          WRITE(*,1001) ZTETA(J,5),RRISO(J,IEDGE),RRBND(J,5,1),Z2RAT(J)
C        ENDDO
C1001    FORMAT(4E16.8)

C        COMPUTE METRICS ELEMENTS INSIDE EACH VACUUM REGION BETWEEN TWO WALLS
         DO 10 J1=1,NWBPS-1
         ZSA = ZS(2*NW(J1)-1)
         ZSB = ZS(2*NW(J1+1)-1)
         ZSBA  = ZSB-ZSA
         NWUP = 2*NW(J1+1)-2
         IF (J1.EQ.(NWBPS-1)) NWUP = NWUP + 1
         DO 9 J9=2*NW(J1)-1,NWUP
         ZS1 = ZS(J9)

         J99 = (J9+1)/2
         IF (J99*2-1.NE.J9) J99=0
         JW = 0
         DO L=1,NWBPS
            IF (ABS(NW(L)-J99)<1.0E-1) JW=J99
         ENDDO
         DO L=1,5
            IF (ABS(RCOIL(L)-ZS1)<1.0E-8) JW=J99
         ENDDO

C        MODIFY CHIISO INSIDE VACUUM TO DENSIFY CHI-MESH IN LFS
         CALL BOUNDDCDT(NMGAUS*NT1,ZS1,ZTETA(1,5),CHIISO,CHIISON,
     &                  ZDCDT(1,1),ZDCDS)
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISON,1,ZCHI(I)) - 1
 
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISON(ICHISO+1) - CHIISON(ICHISO)
           ZA(I) = (CHIISON(ICHISO+1) - ZCHI(I)) / ZH
           ZB(I) = (ZCHI(I) - CHIISON(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISON(ICHISO+1) - ZCHI(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHI(I) - CHIISON(ICHISO)) / 6.

         ENDDO
         
         DO 8 J8=1,NT1*NMGAUS
         IG = J8

         ZDRDT =Z2RAT(IG)+(ZS1-ZSA)*Z2RAST(IG)+(ZS1-ZSA)**Z2GP*Z2GRT(IG)
         ZDZDT =Z2ZAT(IG)+(ZS1-ZSA)*Z2ZAST(IG)+(ZS1-ZSA)**Z2GP*Z2GZT(IG)
         ZDRDC =ZDRDT/ZDCDT(IG,1)
         ZDZDC =ZDZDT/ZDCDT(IG,1)
         ZDRDS =Z2RAS(IG)+Z2GP*(ZS1-ZSA)**(Z2GP-1.)*Z2GR(IG) - 
     &          ZDRDC*ZDCDS(IG)
         ZDZDS =Z2ZAS(IG)+Z2GP*(ZS1-ZSA)**(Z2GP-1.)*Z2GZ(IG) -
     &          ZDZDC*ZDCDS(IG)
         ZR    =Z2RA(IG)+(ZS1-ZSA)*Z2RAS(IG)+(ZS1-ZSA)**Z2GP*Z2GR(IG)
         ZZ    =Z2ZA(IG)+(ZS1-ZSA)*Z2ZAS(IG)+(ZS1-ZSA)**Z2GP*Z2GZ(IG)

         ZZG11L = ZDRDS**2 + ZDZDS**2
         ZZG12L = ZDRDS*ZDRDC + ZDZDS*ZDZDC
         ZZG22L = ZDRDC**2 + ZDZDC**2
         ZZG33L = ZR**2
         ZJAC1  = ZR*(ZDRDS*ZDZDC - ZDRDC*ZDZDS)

         ZG11LI(J8) = ZZG11L
         ZG12LI(J8) = ZZG12L
         ZG22LI(J8) = ZZG22L
         ZG33LI(J8) = ZZG33L
         ZRJAI (J8) = ZJAC1
         ZRANI (J8) = ZR
         ZZANI (J8) = ZZ
    8    CONTINUE
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C     3. FOR EACH ARRAY: COMPUTE VALUES ON ZCHI USING A PERIODIC
C     CUBIC SPLINE FIT AND COMPUTE FULL FOURIER TRANSFORM

         CALL SPLCHI(CHIISON,ZRANI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRAN,INCHI,ICHIISO,ZA,ZB,ZC,ZD)

         CALL SPLCHI(CHIISON,ZZANI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZZAN,INCHI,ICHIISO,ZA,ZB,ZC,ZD)

         CALL SPLCHI(CHIISON,ZRJAI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRJA,INCHI,ICHIISO,ZA,ZB,ZC,ZD)

         CALL SPLCHI(CHIISON,ZG11LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG11L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)

         CALL SPLCHI(CHIISON,ZG22LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG22L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)

         CALL SPLCHI(CHIISON,ZG12LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG12L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)

         CALL SPLCHI(CHIISON,ZG33LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG33L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)

C-----------------------------------------------------------------------

      IF (J9 .GT. 1) GOTO 20

         CLOSE(NRMAR)
         OPEN(NRMAR,FILE='OUTVMAR',STATUS='NEW',FORM='FORMATTED')
         REWIND NRMAR

         WRITE(NRMAR,*) NV1, INCHI, NWBPS 
 20   CONTINUE 

      WRITE(NRMAR,*) ZS1
      WRITE(NRMAR,*) (ZRJA(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG11L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG22L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG33L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZG12L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZRAN(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZZAN(J),J=1,INCHI)
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
    9    CONTINUE
 
C        RECOMPUTE Z2* QUANTITIES AT NEW WALL POSITION ZSA
         IF (J1.LE.NWBPS-2) THEN
         DO IG=1,NT1*NMGAUS
            Z2RA(IG)  = RRBND(IG,5,J1+1)
            Z2ZA(IG)  = RZBND(IG,5,J1+1)
            Z2RAS(IG) = Z2RAS(IG) + Z2GP*ZSBA**(Z2GP-1.)*Z2GR(IG)
            Z2ZAS(IG) = Z2ZAS(IG) + Z2GP*ZSBA**(Z2GP-1.)*Z2GZ(IG)
            ZSBA      = ZS(2*NW(J1+2)-1) - ZS(2*NW(J1+1)-1) 
         Z2GR(IG)=(RRBND(IG,5,J1+2)-Z2RA(IG)-Z2RAS(IG)*ZSBA)/ZSBA**Z2GP
         Z2GZ(IG)=(RZBND(IG,5,J1+2)-Z2ZA(IG)-Z2ZAS(IG)*ZSBA)/ZSBA**Z2GP
         ENDDO

         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2RA,Z2RAT)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2ZA,Z2ZAT)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2RAS,Z2RAST)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2ZAS,Z2ZAST)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2GR,Z2GRT)      
         CALL BOUNDDT2(NMGAUS*NT1,ZTETA(1,5),Z2GZ,Z2GZT)
         ENDIF      
   10    CONTINUE

C     WRITE CONDUCTIVITY COEFFICIENTS FOR ALL RESISTIVE WALLS
      DO J1=2,NWBPS-1
         DO J=1,NT1*NMGAUS
            ZCNDI(J) = RCBND(J,5,J1)
         ENDDO
         CALL SPLCHI(CHIISON,ZCNDI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZCND,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
         WRITE(NRMAR,*) (ZCND(J),J=1,INCHI) 
         DO J=1,NT1*NMGAUS
            ZCNDI(J) = (RCBND(J,1,J1) + 8*(RCBND(J,3,J1) - 
     -               RCBND(J,2,J1)) - RCBND(J,4,J1)) / (12. * ZEPS)
            ZCNDI(J) = ZCNDI(J)/ZDCDT(J,1)
         ENDDO
         CALL SPLCHI(CHIISON,ZCNDI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZCND,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
         WRITE(NRMAR,*) (ZCND(J),J=1,INCHI) 
      ENDDO

      CLOSE(NRMAR)

 2004 FORMAT (3E16.8)

      RETURN
      END

         SUBROUTINE VACUUMR
C        ##################
C
C                                        AUTHOR O. SAUTER, CRPP-EPFL
C                                               A. BONDESON, CHALMERS
***********************************************************************
*                                                                     *
* C2SM23  COMPUTE EQV TERMS FOR MARS USING FAST FFT TRANSFORM ON AN   *
*         2**L EQUIDISTANT CHI-MESH. CUBIC SPLINE INTERPOLATION OF THE*
*        TERMS FROM GAUSSIAN THETA-MESH TO EQUIDISTANT CHI-MESH.      *
*                                                                     *
***********************************************************************
C

         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMRCH.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
C
      INTEGER NPISOC
      PARAMETER (NPISOC = NPMGS*NTP1)
C
         DIMENSION
     R     ZBND(NPMGS*NTP1,4), ZS(2*NPV1), ZTETA(NPMGS*NTP1,4),
     R     ICHIISO(2*NPCHI1), 
     R     ZG11(NPMGS*NTP1), ZG22(NPMGS*NTP1), ZG33(NPMGS*NTP1),
     R     ZG12(NPMGS*NTP1), ZJAC1(NPMGS*NTP1)
      DIMENSION
     R     ZG11LI(NPISOC),      ZG22LI(NPISOC),
     R     ZG33LI(NPISOC),      ZG12LI(NPISOC),
     R     ZRANLI(NPISOC),      ZZANLI(NPISOC),
     R     ZRJAI(NPISOC),
     R     ZRI(NPISOC),         ZZI(NPISOC)
         DIMENSION
     R     ZCHI(NP2CHI), ZD2FUN(NPISOC), ZWORK(NPISOC,3),
     R     ZA(NP2CHI), ZB(NP2CHI), ZC(NP2CHI), ZD(NP2CHI),
     R     ZG11L(NP2CHI),      ZG22L(NP2CHI),
     R     ZG33L(NP2CHI),      ZG12L(NP2CHI),
     R     ZRJA(NP2CHI),
     R     ZRAN(NP2CHI),       ZZAN(NP2CHI)
C
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IEDGE = 2*NPSI
         ZEPS = 1.E-3
C
C     PREPARE COEFFICIENTS FOR FFT AND CUBIC SPLINE
C        
         INCHI = 2*NCHI
C
         ZDCHI = 2. * CPI / DFLOAT(INCHI)
         DO I=1,INCHI+1
           ZCHI(I) = DFLOAT(I-1)*ZDCHI
         ENDDO
C
         IGCHISO = NMGAUS*NT1
C     ASSUMES THAT CHIISO HAS ALREADY ARRAY AT BOUNDARY (FROM CALL FOURFFT)
C     (OTHERWISE, ADD CALL GCHI(2*NPSI))
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHI(I)) - 1
C
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
           ZA(I) = (CHIISO(ICHISO+1) - ZCHI(I)) / ZH
           ZB(I) = (ZCHI(I) - CHIISO(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISO(ICHISO+1) - ZCHI(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHI(I) - CHIISO(ICHISO)) / 6.
C
         ENDDO
C
C     PRECOMPUTE SOME ARRAYS
C
         DO 1 J1=1,NV
C
         ZS(2*(J1-1)+1) = CSV(J1)
         ZS(2*J1      ) = CSMV(J1)
C
    1    CONTINUE
C
         ZS(2*NV+1) = CSV(NV1)
C
         DO 2 J2=1,NMGAUS*NT1
C
         ZTETA(J2,1) = TETPSI(J2,IEDGE) - 2. * ZEPS
         ZTETA(J2,2) = TETPSI(J2,IEDGE) -      ZEPS
         ZTETA(J2,3) = TETPSI(J2,IEDGE) +      ZEPS
         ZTETA(J2,4) = TETPSI(J2,IEDGE) + 2. * ZEPS
C
    2    CONTINUE
C
         CALL BOUND(NMGAUS*NT1,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,4),ZBND(1,4))
C
C     PRECOMPUTE TERMS NON-DEPENDING ON VACUUM S VALUE
C
         DO 3 J3=1,NMGAUS*NT1
C
         ZCOST = COS(TETPSI(J3,IEDGE))
         ZSINT = SIN(TETPSI(J3,IEDGE))
         ZBNDT = BNDISO(J3,IEDGE)
         ZRHO  = ZBNDT
         ZR    = RRISO(J3,IEDGE)
         ZZ    = RZISO(J3,IEDGE)
         ZGP   = GPISO(J3,IEDGE)
         ZDPDS = DPSISO(J3,IEDGE)
         ZDPDT = 0.
C
         ZDRSDT = (ZBND(J3,1) + 8*(ZBND(J3,3) - ZBND(J3,2)) -
     -             ZBND(J3,4)) / (12. * ZEPS)
         ZJAC   = CP(IEDGE) * ZR**NER * ZGP**NEGP
C
         ZDCDT = ZRHO * ZBNDT * ZR / (ZJAC * ZDPDS)
C
         ZDSDR = (ZDRSDT * ZSINT + ZBNDT * ZCOST) / ZBNDT**2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBNDT * ZSINT - ZDRSDT * ZCOST) / ZBNDT**2
         ZDTDZ = ZCOST / ZRHO
C
         ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         ZDRDC = - ZJAC * ZDPDZ / ZR
         ZDZDC = ZJAC * ZDPDR / ZR
C
         ZRW = ZR - R0W
         ZZW = ZZ - RZ0W
         ZG11(J3) = ZRW**2 + ZZW**2
         ZG22(J3) = ZDRDC**2 + ZDZDC**2
         ZG12(J3) = ZRW * ZDRDC + ZZW * ZDZDC
         ZG33(J3) = ZRW
         ZJAC1(J3)= ZRW * ZDZDC - ZZW * ZDRDC

         ZRI(J3) = ZR
         ZZI(J3) = ZZ
C
    3    CONTINUE
C
C     AT EACH S VACUUM LOCATION: COMPUTE TERMS ON GAUSSIAN MESH, INTERPOLATE
C     THEM ON EQUIDISTANT CHI-MESH AND COMPUTE FOURIER TRANSFORM
C
         DO 9 J9=1,2*NV+1
           ZS1   = ZS(J9)
C
C     USE ZBND AND ZTETA MEMORY SPACE
         DO 8 J8=1,NT1*NMGAUS
C
           ZRJAI (J8) = ZS1 * (R0W + ZS1 * ZG33(J8)) * ZJAC1(J8)
           ZG11LI(J8) = ZG11(J8)
           ZG22LI(J8) = ZS1**2 * ZG22(J8)
           ZG33LI(J8) = (R0W + ZS1 * ZG33(J8))**2
           ZG12LI(J8) = ZS1  * ZG12(J8)
           ZRANLI(J8) = R0W + ZS1 * (ZRI(J8)-R0W)
           ZZANLI(J8) = RZ0W + ZS1 * (ZZI(J8)-RZ0W)
C
    8    CONTINUE
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C     3. FOR EACH ARRAY: COMPUTE VALUES ON ZCHI USING A PERIODIC
C     CUBIC SPLINE FIT AND COMPUTE FULL FOURIER TRANSFORM
C
C     EQL
C
         CALL SPLCHI(CHIISO,ZRJAI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRJA,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG11LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG11L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG22LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG22L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG12LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG12L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG33LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG33L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZRANLI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRAN,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZZANLI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZZAN,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
C-----------------------------------------------------------------------
C
      IF (J9 .GT. 1) GOTO 20
C
      NRMAR = 77
         CLOSE(NRMAR)
         OPEN(NRMAR,FILE='OUTVMAR',STATUS='NEW',FORM='FORMATTED')
         REWIND NRMAR
C
         WRITE(NRMAR,*) NV1, INCHI, 2 
C
 20   CONTINUE 
C
      WRITE(NRMAR,*) ZS1
      WRITE(NRMAR,*) (ZRJA(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG11L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG22L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG33L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZG12L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZRAN(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZZAN(J),J=1,INCHI)
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
 9    CONTINUE
      CLOSE(NRMAR)
C
	 RETURN
	 END

	 SUBROUTINE VACUUMRX
C        ##################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM13 COMPUTE VACUUM EQ'S FOR MARS (SEE SECTION 5.4.2 AND TABLE 2  *
*        IN PUBLICATION)                                              *
*                                                                     *
***********************************************************************
C
	 INCLUDE 'DECLAR.inc'
	 INCLUDE 'COMDIM.inc'
	 INCLUDE 'COMCON.inc'
	 INCLUDE 'COMESH.inc'
	 INCLUDE 'COMETA.inc'
         INCLUDE 'COMIOD.inc'
	 INCLUDE 'COMISO.inc'
	 INCLUDE 'COMMAP.inc'
	 INCLUDE 'COMNUM.inc'
	 INCLUDE 'COMPHY.inc'
         INCLUDE 'COMRCH.inc'
	 INCLUDE 'COMSOL.inc'
	 INCLUDE 'COMSUR.inc'
	 INCLUDE 'COMVAC.inc'
C
      INTEGER NPISOC
      PARAMETER (NPISOC = NPMGS*NTP1)
C
	 DIMENSION
     R     ZBND(NPISOC,4),     ZDCDT(NPISOC),    
     R     ZDR1DC(NPISOC),     ZDZ1DC(NPISOC),    
     R     ZR1(NPISOC),        ZZ1(NPISOC),
     R     ZS(2*NPV1),         ZTETA(NPISOC,4)
      DIMENSION
     R     ZG11LI(NPISOC),      ZG22LI(NPISOC),
     R     ZG33LI(NPISOC),      ZG12LI(NPISOC),
     R     ZRJAI(NP2CHI)
         DIMENSION
     R     ICHIISO(NP2CHI),
     R     ZCHI(NP2CHI), ZD2FUN(NPISOC), ZWORK(NPISOC,3),
     R     ZA(NP2CHI), ZB(NP2CHI), ZC(NP2CHI), ZD(NP2CHI),
     R     ZG11L(NP2CHI),      ZG22L(NP2CHI),
     R     ZG33L(NP2CHI),      ZG12L(NP2CHI),
     R     ZRJA(NP2CHI)
C
C        COMPLEX     ZCGJA
C        DIMENSION           ZCGJA(2*NPV1)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
      write(*,*) '********** VACUUMR *********'
	 IEDGE = 2*NPSI
	 ZEPS = 1.E-3
C
	 DO 1 J1=1,NV
C
	 ZS(2*(J1-1)+1) = CSV(J1)
	 ZS(2*J1      ) = CSMV(J1)
C
    1    CONTINUE
C
	 ZS(2*NV+1) = CSV(NV1)
C
	 DO 2 J2=1,NMGAUS*NT1
C
	 ZTETA(J2,1) = TETPSI(J2,IEDGE) - 2. * ZEPS
	 ZTETA(J2,2) = TETPSI(J2,IEDGE) -      ZEPS
	 ZTETA(J2,3) = TETPSI(J2,IEDGE) +      ZEPS
	 ZTETA(J2,4) = TETPSI(J2,IEDGE) + 2. * ZEPS
C
    2    CONTINUE
C
	 CALL BOUND(NMGAUS*NT1,ZTETA(1,1),ZBND(1,1))
	 CALL BOUND(NMGAUS*NT1,ZTETA(1,2),ZBND(1,2))
	 CALL BOUND(NMGAUS*NT1,ZTETA(1,3),ZBND(1,3))
	 CALL BOUND(NMGAUS*NT1,ZTETA(1,4),ZBND(1,4))
C
	 DO 3 J3=1,NMGAUS*NT1
C
	 ZCOST = COS(TETPSI(J3,IEDGE))
	 ZSINT = SIN(TETPSI(J3,IEDGE))
	 ZBNDT = BNDISO(J3,IEDGE)
	 ZRHO  = ZBNDT
	 ZR    = RRISO(J3,IEDGE)
	 ZZ    = RZISO(J3,IEDGE)
	 ZGP   = GPISO(J3,IEDGE)
	 ZDPDS = DPSISO(J3,IEDGE)
	 ZDPDT = 0.
C
	 ZDRSDT = (ZBND(J3,1) + 8*(ZBND(J3,3) - ZBND(J3,2)) -
     -             ZBND(J3,4)) / (12. * ZEPS)
	 ZJAC   = CP(IEDGE) * ZR**NER * ZGP**NEGP
C
	 ZDCDT(J3) = ZRHO * ZBNDT * ZR / (ZJAC * ZDPDS)
	 ZR1(J3)   = ZR - R0W
	 ZZ1(J3)   = ZZ - RZ0W
C
	 ZDSDR = (ZDRSDT * ZSINT + ZBNDT * ZCOST) / ZBNDT**2
	 ZDTDR = - ZSINT / ZRHO
	 ZDSDZ = (ZBNDT * ZSINT - ZDRSDT * ZCOST) / ZBNDT**2
	 ZDTDZ = ZCOST / ZRHO
C
	 ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
	 ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
	 ZDR1DC(J3) = - ZJAC * ZDPDZ / ZR
	 ZDZ1DC(J3) = ZJAC * ZDPDR / ZR
C
    3    CONTINUE
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C     2. COMPUTE EQUIDISTANT CHI MESH, KMMAX INTERVALS
C        
         INCHI = 2*NCHI
C
         ZDCHI = 2. * CPI / DFLOAT(INCHI)
         DO I=1,INCHI+1
           ZCHI(I) = DFLOAT(I-1)*ZDCHI
         ENDDO
C
C     3. PREPARE COEFFICIENTS FOR THE CUBIC SPLINE FIT DEPENDING
C     ONLY ON RELATIVE POSITION OF ZCHI(I) WITH RESPECT TO CHIISO
C
         IGCHISO = NMGAUS*NT1
         KPSI = 2*NPSI
         CALL GCHI(KPSI)
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHI(I)) - 1
C
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
           ZA(I) = (CHIISO(ICHISO+1) - ZCHI(I)) / ZH
           ZB(I) = (ZCHI(I) - CHIISO(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISO(ICHISO+1) - ZCHI(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHI(I) - CHIISO(ICHISO)) / 6.
C
         ENDDO
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
	 DO 9 JV=1,2*NV+1
C
	 ZS1   = ZS(JV)
	 DO 8 IG=1,NT1*NMGAUS
C
	 ZR    = ZR1(IG)
	 ZZ    = ZZ1(IG)
	 ZDRDC = ZDR1DC(IG)
	 ZDZDC = ZDZ1DC(IG)
C
	 ZG11LI(IG) =           ZR**2 + ZZ**2
	 ZG22LI(IG) = ZS1**2 * (ZDRDC**2 + ZDZDC**2)
	 ZG12LI(IG) =   ZS1  * (ZR * ZDRDC + ZZ * ZDZDC)
	 ZG33LI(IG) =          (R0W + ZS1 * ZR)**2
	 ZRJAI(IG) = ZS1 * (R0W + ZS1 * ZR) * (ZR * ZDZDC - ZZ * ZDRDC)
C
    8    CONTINUE
      write(*,*) 'after 8  zrjai(1)=',zrjai(1)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C     3. FOR EACH ARRAY: COMPUTE VALUES ON ZCHI USING A PERIODIC
C     CUBIC SPLINE FIT AND COMPUTE FULL FOURIER TRANSFORM
C
C     EQL
C
         CALL SPLCHI(CHIISO,ZRJAI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRJA,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG11LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG11L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG22LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG22L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG12LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG12L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG33LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG33L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
C-----------------------------------------------------------------------
C
      IF (JV .GT. 1) GOTO 20
C
      NRMAR = 77
         CLOSE(NRMAR)
         OPEN(NRMAR,FILE='OUTVMAR',STATUS='NEW',FORM='UNFORMATTED')
         REWIND NRMAR
C
         WRITE(NRMAR) NV1, INCHI 
C
 20   CONTINUE 
C
      write(*,*) ' after 20, zrja(1)',zrja(1)
      WRITE(NRMAR) ZS1
      WRITE(NRMAR) (ZRJA(J),J=1,INCHI) 
      WRITE(NRMAR) (ZG11L(J),J=1,INCHI) 
      WRITE(NRMAR) (ZG22L(J),J=1,INCHI) 
      WRITE(NRMAR) (ZG33L(J),J=1,INCHI)
      WRITE(NRMAR) (ZG12L(J),J=1,INCHI)
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
    9    CONTINUE
      CLOSE(NRMAR)
C
	 RETURN
	 END
C*DECK C2SM14
C*CALL PROCESS
         SUBROUTINE VLION
C        ################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM14 COMPUTE AUXILIARY QUANTITIES AT PLASMA SURFACE FOR LION      *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   ZT1(NPCHI),     ZT2(NPCHI),
     R   ZBND1(NPCHI),   ZBND2(NPCHI)
C
         ZISQR6 = 1. / SQRT(6.)
         ZEPS   = 1.E-4
C
C         TETVAC(1) = .5 * (TETVAC(1) - TETVAC(2))
C         TETVAC(2) = - TETVAC(1)
C         RHOVAC(1) = .5 * (RHOVAC(1) + RHOVAC(2))
C         RHOVAC(2) = RHOVAC(1)
C         TETVACM(1)        = 0.
C         TETVACM(NCHI1)    = 2. * CPI
C         TETVACM(NCHI/2+1) = CPI
C
C         DO 20 J20=3,NCHI/2+1
C
C         TETVAC(J20) = .5 * (TETVAC(J20) - TETVAC(NCHI-J20+3) + 2.*CPI) 
C         RHOVAC(J20) = .5 * (RHOVAC(J20) + RHOVAC(NCHI-J20+3))
C         TETVAC(NCHI-J20+3) = 2. * CPI - TETVAC(J20)
C         RHOVAC(NCHI-J20+3) = RHOVAC(J20)
C
C  20     CONTINUE
C
C         DO 21 J21=2,NCHI/2
C
C         TETVACM(J21) = .5*(TETVACM(J21) - TETVACM(NCHI-J21+2) + 2.*CPI) 
C         RHOVACM(J21) = .5*(RHOVACM(J21) + RHOVACM(NCHI-J21+2))
C
C  21     CONTINUE
C
         DO 1 J1=1,NCHI
C
         TETVACI(J1,2) = .5 * (TETVAC(J1) + TETVAC(J1+1))
         TETVACI(J1,1) = ZISQR6 * (TETVAC(J1) - TETVAC(J1+1)) + 
     +                   TETVACI(J1,2)
         TETVACI(J1,3) = ZISQR6 * (TETVAC(J1+1) - TETVAC(J1)) + 
     +                   TETVACI(J1,2)
C
    1    CONTINUE
C
         BPS( 1) = RMAG
         BPS(12) = RZMAG
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
         CALL BOUND(NCHI,TETVACI(1,1),RHOVACI(1,1))
         CALL BOUND(NCHI,TETVACI(1,2),RHOVACI(1,2))
         CALL BOUND(NCHI,TETVACI(1,3),RHOVACI(1,3))
C
         DO 4 J4=1,3
C
         DO 2 J2=1,NCHI
C
         ZT1(J2) = TETVACI(J2,J4) - ZEPS
         ZT2(J2) = TETVACI(J2,J4) + ZEPS
C
    2    CONTINUE
C
         CALL BOUND(NCHI,ZT1,ZBND1)
         CALL BOUND(NCHI,ZT2,ZBND2)
C
         DO 3 J3=1,NCHI
C
         DRHOPI(J3,J4) = .5 * (ZBND2(J3) - ZBND1(J3)) / ZEPS
C                
    3    CONTINUE
    4    CONTINUE
C
         BPS( 1) = R0
         BPS(12) = RZ0
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
         RETURN
         END
C*DECK C2SM16
C*CALL PROCESS
         SUBROUTINE STCHPS(KP,N,PCHI,PTMAP,PD2TM,PSMAP,PD2SM,
     ,                     PSIG,PTET,KT0)
C        ####################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM16 EVALUATE (SIGMA,THETA) COODINATES OF (S,CHI) NODES REQUIRED  *
*        BY NOVA-W AND PEST                                           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   KT0(N),
     R   PCHI(N),      PD2SM(NTP2), PD2TM(NTP2), PSIG(N), 
     R   PSMAP(NTP2),  PTET(N),     PTMAP(NTP2)
C
         DO 1 JG=1,N
            IC = 1
            DO 1 JT = 1,NT2
               IF (IC.EQ.0) GOTO 1
               KT0(JG) = JT-1
               IF (CHIN(JT,KP).GE.PCHI(JG)) IC = 0
 1       CONTINUE
C
         DO 2 J2=1,N
C
         ICHI = KT0(J2)
C
         IF (ICHI .LT. 1)   ICHI = 1
         IF (ICHI .GT. NT1) ICHI = NT1
C
         ZH = CHIN(ICHI+1,KP) - CHIN(ICHI,KP)
         ZA = (CHIN(ICHI+1,KP) - PCHI(J2)) / ZH
         ZB = (PCHI(J2) - CHIN(ICHI,KP)) / ZH
         ZC = (ZA + 1) * (ZA - 1) * ZH * 
     *        (CHIN(ICHI+1,KP) - PCHI(J2)) / 6.
         ZD = (ZB + 1) * (ZB - 1) * ZH * 
     *        (PCHI(J2) - CHIN(ICHI,KP)) / 6.
C 
         PTET(J2) = ZA*PTMAP(ICHI) + ZB*PTMAP(ICHI+1) +
     +              ZC*PD2TM(ICHI) + ZD*PD2TM(ICHI+1)
C
         IF (KP .EQ. 2*NPSI) THEN
C
           PSIG(J2) = 1.
C
         ELSE
C
           PSIG(J2) = ZA*PSMAP(ICHI) + ZB*PSMAP(ICHI+1) +
     +                ZC*PD2SM(ICHI) + ZD*PD2SM(ICHI+1)
C
         ENDIF
C
   2     CONTINUE
C
         RETURN
         END
C*DECK C2SM17
C*CALL PROCESS
         SUBROUTINE JNOVAW(KP,N,PTETCP,PSIGCP,PR,PZ,PJAC)
C        ################################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM17 EVALUATE JACOBIAN (EQ. (14) OF PUBLICATION) AT (S,CHI) NODES *
*        REQUIRED BY NOVA-W AND PEST                                  *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IS0(NPCHI+3),       IT0(NPCHI+3),       IC(NPCHI+3),
     R   PSIGCP(N),          PTETCP(N),
     R   PR(N),              
     R   PJAC(N),            PZ(N),
     R   ZBND(NPCHI+3,5),    ZDBDS(NPCHI+3,16),  
     R   ZDBDT(NPCHI+3,16),  ZPCEL(NPCHI+3,16),  
     R   ZS(NPCHI+3),        ZS1(NPCHI+3),       
     R   ZS2(NPCHI+3),       ZTETA(NPCHI+3,5),   
     R   ZT(NPCHI+3),        ZT1(NPCHI+3),       
     R   ZT2(NPCHI+3)
C
         ZEPS  = 1.E-3
C
         DO 1 J1=1,N
C
         ZTETA(J1,1) = PTETCP(J1)
         ZTETA(J1,2) = PTETCP(J1) - 2. * ZEPS
         ZTETA(J1,3) = PTETCP(J1) -      ZEPS
         ZTETA(J1,4) = PTETCP(J1) +      ZEPS
         ZTETA(J1,5) = PTETCP(J1) + 2. * ZEPS
C
    1    CONTINUE
C
         CALL BOUND(N,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(N,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(N,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(N,ZTETA(1,4),ZBND(1,4))
         CALL BOUND(N,ZTETA(1,5),ZBND(1,5))
C
         CALL RESETI(IC,N,1)
         DO 2 JT = 1,NT1
            DO 2 JG=1,N
               IF (IC(JG).EQ.0) GOTO 2
               IT0(JG) = JT-1
               IF (PTETCP(JG).LE.CT(JT)) IC(JG)  = 0
 2       CONTINUE
         CALL RESETI(IC,N,1)
         DO 3 JS = 1,NS1
            DO 3 JG=1,N
               IF (IC(JG).EQ.0) GOTO 3
               IS0(JG) = JS-1
               IF (PSIGCP(JG).LE.CSIG(JS)) IC(JG)  = 0
 3       CONTINUE
C
         DO 4 J4=1,N
C
         ZT(J4) = PTETCP(J4)
         ZS(J4) = PSIGCP(J4)
C
         IF (IS0(J4) .GT. NS) IS0(J4) = NS
         IF (IS0(J4) .LT. 1)  IS0(J4) = 1
         IF (IT0(J4) .GT. NT) IT0(J4) = NT
         IF (IT0(J4) .LT. 1)  IT0(J4) = 1
C
         ZS1(J4) = CSIG(IS0(J4))
         ZS2(J4) = CSIG(IS0(J4)+1)
         ZT1(J4) = CT(IT0(J4))
         ZT2(J4) = CT(IT0(J4)+1)
C
 4       CONTINUE
C
         CALL PSICEL(IS0,IT0,N,NPCHI+3,ZPCEL,CPSICL)
         CALL BASIS2(N,NPCHI+3,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT)
C
         DO 5 J5=1,N
C
         ZDRSDT = (ZBND(J5,2) + 8*(ZBND(J5,4) - ZBND(J5,3)) -
     -             ZBND(J5,5)) / (12. * ZEPS)
C
         ZDPDS = ZDBDS(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDS(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDS(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDS(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDS(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDS(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDS(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDS(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDS(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDS(J5,10) * ZPCEL(J5,10) +
     +           ZDBDS(J5,11) * ZPCEL(J5,11) +
     +           ZDBDS(J5,12) * ZPCEL(J5,12) +
     +           ZDBDS(J5,13) * ZPCEL(J5,13) +
     +           ZDBDS(J5,14) * ZPCEL(J5,14) +
     +           ZDBDS(J5,15) * ZPCEL(J5,15) +
     +           ZDBDS(J5,16) * ZPCEL(J5,16)
C
         ZDPDT = ZDBDT(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDT(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDT(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDT(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDT(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDT(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDT(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDT(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDT(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDT(J5,10) * ZPCEL(J5,10) +
     +           ZDBDT(J5,11) * ZPCEL(J5,11) +
     +           ZDBDT(J5,12) * ZPCEL(J5,12) +
     +           ZDBDT(J5,13) * ZPCEL(J5,13) +
     +           ZDBDT(J5,14) * ZPCEL(J5,14) +
     +           ZDBDT(J5,15) * ZPCEL(J5,15) +
     +           ZDBDT(J5,16) * ZPCEL(J5,16)
C
         ZFP    = (ZDPDS**2 + (ZDPDT / PSIGCP(J5) - ZDPDS * ZDRSDT /
     -             ZBND(J5,1))**2) / ZBND(J5,1)**2
         ZGRADP = SQRT(ZFP)
C
         ZCOST = COS(ZTETA(J5,1))
         ZSINT = SIN(ZTETA(J5,1))
C
         ZRHO   = PSIGCP(J5) * ZBND(J5,1)
         PR(J5) = ZRHO * ZCOST + R0
         PZ(J5) = ZRHO * ZSINT + RZ0
C
         PJAC(J5) = CP(KP) * PR(J5)**NER * ZGRADP**NEGP
C
 5       CONTINUE
C
         RETURN
         END
C*DECK C2SM18
C*CALL PROCESS
         SUBROUTINE TPSI(NP1,KP)
C        #######################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM18 COMPUTE (SIGMA,THETA) COODINATES OF ALL QUADRATURE AND MESH  *
*        POINTS REQUIRED BY PENN (SEE SECTION 5.4.5 OF PUBLICATION)   *
*        OR OF THE MESH POINTS REQUIRED BY XTOR                       *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IC((NPMGS+1)*NPCHI),        IT0((NPMGS+1)*NPCHI),
     R   ZTETA((NPMGS+1)*NPCHI),
     R   ZD2SIG(NTP2), ZD2TET(NTP2), ZTET(NTP2),   ZTETN(NTP2),
     R   ZA1(NTP2),    ZB1(NTP2),    ZC1(NTP2),    ZBND(NTP2)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IF (NIDEAL .EQ. 4) THEN
           IGMAX = NCHI*(NMGAUS+1)
           CALL SCOPY(IGMAX,CTPEN,1,ZTETA,1)
         ELSE IF (NIDEAL .EQ. 5) THEN
           IGMAX = NTNOVA
           CALL SCOPY(IGMAX,CTXT,1,ZTETA,1)
         ELSE
           PRINT*,'TPSI REQUIRES NIDEAL=4 OR 5'
           STOP
         ENDIF
C
         CALL SCOPY(NT2,TETMAP(1,KP),1,ZTET,1)
C
         DO 1 J1=2,NT2
C
         IF (ZTET(J1) .LT. ZTET(J1-1)) THEN
C
            ZTET(J1) = ZTET(J1) + 2. * CPI * (1. +
     +                 INT(.5 * ABS(ZTET(J1) - ZTET(J1-1)) / CPI))       
C
         ENDIF
C
    1    CONTINUE
C
         CALL BOUND(NT2,ZTET,ZBND)
C
         DO 2 J2=1,NT2
C
         ZR        = SIGMAP(J2,KP) * ZBND(J2) * COS(ZTET(J2)) + R0
         ZZ        = SIGMAP(J2,KP) * ZBND(J2) * SIN(ZTET(J2)) + RZ0
         ZTETN(J2) = ATAN2(ZZ - RZMAG,ZR - RMAG)
C
    2    CONTINUE
C
         DO 3 J3=2,NT2
C
         IF (ZTETN(J3) .LT. ZTETN(J3-1)) THEN
C
            ZTETN(J3) = ZTETN(J3) + 2. * CPI * (1. +
     +                 INT(.5 * ABS(ZTETN(J3) - ZTETN(J3-1)) / CPI))       
C
         ENDIF
C
    3    CONTINUE
C
         CALL SPLCY(ZTETN,SIGMAP(1,KP),NT1,RC2PI,ZD2SIG,
     &              ZA1,ZB1,ZC1)
         CALL SPLCYP(ZTETN,ZTET,NT1,RC2PI,RC2PI,ZD2TET,
     &               ZA1,ZB1,ZC1)
C
         ZD2SIG(NT2) = ZD2SIG(1) 
         ZD2TET(NT2) = ZD2TET(1) 
C
         CALL RESETI(IC,IGMAX,1)
         DO 4 JT=1,NT2
           DO 4 JG = 1,IGMAX
               IF (IC(JG).EQ.0) GOTO 4
               IT0(JG) = JT-1
               IF (ZTETA(JG).LE.ZTETN(JT)) IC(JG)  = 0
 4       CONTINUE
C
         DO 5 J5=1,IGMAX
C
         IT = IT0(J5)
C
         IF (IT .LT. 1)   IT = 1
         IF (IT .GT. NT1) IT = NT1
C
         ZH = ZTETN(IT+1) - ZTETN(IT)
         ZA = (ZTETN(IT+1) - ZTETA(J5)) / ZH
         ZB = (ZTETA(J5) - ZTETN(IT)) / ZH
         ZC = (ZA+1) * (ZA-1) * ZH * (ZTETN(IT+1) - ZTETA(J5)) / 6.
         ZD = (ZB+1) * (ZB-1) * ZH * (ZTETA(J5) - ZTETN(IT)) / 6.
C 
         TETPEN(J5) = ZA*ZTET(IT)   + ZB*ZTET(IT+1) +
     +                ZC*ZD2TET(IT) + ZD*ZD2TET(IT+1)
C
         IF (TETPEN(J5) .LT. CT(1))
     =                   TETPEN(J5) = TETPEN(J5) + 2.*CPI
         IF (TETPEN(J5) .GT. CT(NT1)) 
     =                   TETPEN(J5) = TETPEN(J5) - 2.*CPI
C
         IF (KP .EQ. NP1) THEN
C
           SIGPEN(J5) = 1.
C
         ELSE
C
           SIGPEN(J5) = ZA*SIGMAP(IT,KP) + ZB*SIGMAP(IT+1,KP) +
     +                  ZC*ZD2SIG(IT)    + ZD*ZD2SIG(IT+1)
C
         ENDIF
C
 5       CONTINUE
C
         RETURN
         END
C*DECK C2SM19
C*CALL PROCESS
         SUBROUTINE OUTPEN(KSPEN)
C        ########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM19 EVALUATE EQ'S FOR PENN (SEE EQ. (39) OF PUBLICATION)         *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
C
         PARAMETER (NPPOINT = NPCHI * (NPMGS + 1))
C
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMPLO.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IC(NPPOINT),        IS0(NPPOINT),       IT0(NPPOINT),
     R   ZBND(NPPOINT,5),    ZDBDS(NPPOINT,16),  
     R   ZDBDT(NPPOINT,16),  ZDBDST(NPPOINT,16), 
     R   ZDPDR(NPPOINT),     ZDPDZ(NPPOINT),
     R   ZD2PRZ(NPPOINT),    ZD2PR2(NPPOINT),    
     R   ZD2PZ2(NPPOINT),    ZDTTDR(NPPOINT),    
     R   ZDTTDZ(NPPOINT),    ZD2TTRZ(NPPOINT), 
     R   ZD2TTR2(NPPOINT),   ZD2TTZ2(NPPOINT),   
     R   ZD2BS2(NPPOINT,16), ZD2BT2(NPPOINT,16), 
     R   ZPCEL(NPPOINT,16),  ZR(NPPOINT),   
     R   ZS(NPPOINT),        ZS1(NPPOINT),       
     R   ZS2(NPPOINT),       ZTETA(NPPOINT,5),   
     R   ZT(NPPOINT),        ZT1(NPPOINT),       
     R   ZT2(NPPOINT),       ZZ(NPPOINT),
     R   ZSURF(NPPOINT)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ZEPS   = 1.E-3
         NPOINT = NCHI * (NMGAUS + 1)
C
         DO 1 J1=1,NPOINT
           ZTETA(J1,1) = TETPEN(J1)
           ZTETA(J1,2) = TETPEN(J1) - 2. * ZEPS
           ZTETA(J1,3) = TETPEN(J1) -      ZEPS
           ZTETA(J1,4) = TETPEN(J1) +      ZEPS
 1         ZTETA(J1,5) = TETPEN(J1) + 2. * ZEPS
C
         CALL BOUND(NPOINT,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NPOINT,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NPOINT,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NPOINT,ZTETA(1,4),ZBND(1,4))
         CALL BOUND(NPOINT,ZTETA(1,5),ZBND(1,5))
C
         CALL RESETI(IC,NPOINT,1)
         DO 2 JT = 1,NT1
            DO 2 JG=1,NPOINT
               IF (IC(JG).EQ.0) GOTO 2
               IT0(JG) = JT-1
               IF (TETPEN(JG).LE.CT(JT)) IC(JG)  = 0
 2       CONTINUE
         CALL RESETI(IC,NPOINT,1)
         DO 3 JS = 1,NS1
            DO 3 JG=1,NPOINT
               IF (IC(JG).EQ.0) GOTO 3
               IS0(JG) = JS-1
               IF (SIGPEN(JG).LE.CSIG(JS)) IC(JG)  = 0
 3       CONTINUE
C
         DO 4 J4=1,NPOINT
           IF (IS0(J4) .GT. NS) IS0(J4) = NS
           IF (IS0(J4) .LT. 1)  IS0(J4) = 1
           IF (IT0(J4) .GT. NT) IT0(J4) = NT
           IF (IT0(J4) .LT. 1)  IT0(J4) = 1
C
           ZT(J4)  = TETPEN(J4)
           ZS(J4)  = SIGPEN(J4)
           ZS1(J4) = CSIG(IS0(J4))
           ZS2(J4) = CSIG(IS0(J4)+1)
           ZT1(J4) = CT(IT0(J4))
 4         ZT2(J4) = CT(IT0(J4)+1)
C
         CALL PSICEL(IS0,IT0,NPOINT,NPPOINT,ZPCEL,CPSICL)
         CALL BASIS3(NPOINT,NPPOINT,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT,
     &               ZDBDST,ZD2BS2,ZD2BT2)
C
         DO 5 J5=1,NPOINT
C
         ZDRSDT = (ZBND(J5,2) + 8*(ZBND(J5,4) - ZBND(J5,3)) -
     -             ZBND(J5,5)) / (12. * ZEPS)
         ZD2RST = (- ZBND(J5,2) + 16. * ZBND(J5,3) -
     -             30. * ZBND(J5,1) + 16. * ZBND(J5,4) -
     -             ZBND(J5,5)) / (12. * ZEPS**2)
C
         ZDPDS = ZDBDS(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDS(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDS(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDS(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDS(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDS(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDS(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDS(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDS(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDS(J5,10) * ZPCEL(J5,10) +
     +           ZDBDS(J5,11) * ZPCEL(J5,11) +
     +           ZDBDS(J5,12) * ZPCEL(J5,12) +
     +           ZDBDS(J5,13) * ZPCEL(J5,13) +
     +           ZDBDS(J5,14) * ZPCEL(J5,14) +
     +           ZDBDS(J5,15) * ZPCEL(J5,15) +
     +           ZDBDS(J5,16) * ZPCEL(J5,16)
C
         ZDPDT = ZDBDT(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDT(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDT(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDT(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDT(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDT(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDT(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDT(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDT(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDT(J5,10) * ZPCEL(J5,10) +
     +           ZDBDT(J5,11) * ZPCEL(J5,11) +
     +           ZDBDT(J5,12) * ZPCEL(J5,12) +
     +           ZDBDT(J5,13) * ZPCEL(J5,13) +
     +           ZDBDT(J5,14) * ZPCEL(J5,14) +
     +           ZDBDT(J5,15) * ZPCEL(J5,15) +
     +           ZDBDT(J5,16) * ZPCEL(J5,16)
C
         ZD2PST = ZDBDST(J5, 1) * ZPCEL(J5, 1) +
     +            ZDBDST(J5, 2) * ZPCEL(J5, 2) +
     +            ZDBDST(J5, 3) * ZPCEL(J5, 3) +
     +            ZDBDST(J5, 4) * ZPCEL(J5, 4) +
     +            ZDBDST(J5, 5) * ZPCEL(J5, 5) +
     +            ZDBDST(J5, 6) * ZPCEL(J5, 6) +
     +            ZDBDST(J5, 7) * ZPCEL(J5, 7) +
     +            ZDBDST(J5, 8) * ZPCEL(J5, 8) +
     +            ZDBDST(J5, 9) * ZPCEL(J5, 9) +
     +            ZDBDST(J5,10) * ZPCEL(J5,10) +
     +            ZDBDST(J5,11) * ZPCEL(J5,11) +
     +            ZDBDST(J5,12) * ZPCEL(J5,12) +
     +            ZDBDST(J5,13) * ZPCEL(J5,13) +
     +            ZDBDST(J5,14) * ZPCEL(J5,14) +
     +            ZDBDST(J5,15) * ZPCEL(J5,15) +
     +            ZDBDST(J5,16) * ZPCEL(J5,16)
C
         ZD2PS2 = ZD2BS2(J5, 1) * ZPCEL(J5, 1) +
     +            ZD2BS2(J5, 2) * ZPCEL(J5, 2) +
     +            ZD2BS2(J5, 3) * ZPCEL(J5, 3) +
     +            ZD2BS2(J5, 4) * ZPCEL(J5, 4) +
     +            ZD2BS2(J5, 5) * ZPCEL(J5, 5) +
     +            ZD2BS2(J5, 6) * ZPCEL(J5, 6) +
     +            ZD2BS2(J5, 7) * ZPCEL(J5, 7) +
     +            ZD2BS2(J5, 8) * ZPCEL(J5, 8) +
     +            ZD2BS2(J5, 9) * ZPCEL(J5, 9) +
     +            ZD2BS2(J5,10) * ZPCEL(J5,10) +
     +            ZD2BS2(J5,11) * ZPCEL(J5,11) +
     +            ZD2BS2(J5,12) * ZPCEL(J5,12) +
     +            ZD2BS2(J5,13) * ZPCEL(J5,13) +
     +            ZD2BS2(J5,14) * ZPCEL(J5,14) +
     +            ZD2BS2(J5,15) * ZPCEL(J5,15) +
     +            ZD2BS2(J5,16) * ZPCEL(J5,16)
C
         ZD2PT2 = ZD2BT2(J5, 1) * ZPCEL(J5, 1) +
     +            ZD2BT2(J5, 2) * ZPCEL(J5, 2) +
     +            ZD2BT2(J5, 3) * ZPCEL(J5, 3) +
     +            ZD2BT2(J5, 4) * ZPCEL(J5, 4) +
     +            ZD2BT2(J5, 5) * ZPCEL(J5, 5) +
     +            ZD2BT2(J5, 6) * ZPCEL(J5, 6) +
     +            ZD2BT2(J5, 7) * ZPCEL(J5, 7) +
     +            ZD2BT2(J5, 8) * ZPCEL(J5, 8) +
     +            ZD2BT2(J5, 9) * ZPCEL(J5, 9) +
     +            ZD2BT2(J5,10) * ZPCEL(J5,10) +
     +            ZD2BT2(J5,11) * ZPCEL(J5,11) +
     +            ZD2BT2(J5,12) * ZPCEL(J5,12) +
     +            ZD2BT2(J5,13) * ZPCEL(J5,13) +
     +            ZD2BT2(J5,14) * ZPCEL(J5,14) +
     +            ZD2BT2(J5,15) * ZPCEL(J5,15) +
     +            ZD2BT2(J5,16) * ZPCEL(J5,16)
C
         ZCOST  = COS(ZTETA(J5,1))
         ZSINT  = SIN(ZTETA(J5,1))
         ZCOS2T = (ZCOST + ZSINT) * (ZCOST - ZSINT)
         ZSIN2T = 2. * ZSINT * ZCOST
C
         ZRHO   = SIGPEN(J5) * ZBND(J5,1)
         ZRHO2  = ZRHO**2
         ZR(J5) = ZRHO * ZCOST + R0
         ZZ(J5) = ZRHO * ZSINT + RZ0
         ZBND2  = ZBND(J5,1)**2
         ZBND3  = ZBND2 * ZBND(J5,1)
         ZFB    = ZBND2 + 2. * ZDRSDT**2 - ZBND(J5,1) * ZD2RST
C
C  COMPUTE FIRST AND SECOND DERIVATIVES OF THETA-TILD WITH
C  RESPECT TO R AND Z
C
         ZDZ    = ZZ(J5) - RZMAG
         ZDZ2   = ZDZ**2
         ZDR    = ZR(J5) - RMAG
         ZDR2   = ZDR**2
         ZDSUM  = ZDR2 + ZDZ2
         ZDSUM2 = ZDSUM**2
C
         ZDTTDR(J5)  = - ZDZ / ZDSUM
         ZDTTDZ(J5)  =   ZDR / ZDSUM
         ZD2TTRZ(J5) =   (ZDZ2 - ZDR2) / ZDSUM2
         ZD2TTR2(J5) =   2. * ZDR * ZDZ / ZDSUM2
         ZD2TTZ2(J5) = - 2. * ZDR * ZDZ / ZDSUM2
C
C  COMPUTE ALL OTHER DERIVATIVES OF GEOMETRIC QUANTITIES
C  RELATED TO THE EQUILIBRIUM MESH AND INVOLVED IN THE
C  COMPUTATION OF ALL FIRST AND SECOND DERIVATIVES OF PSI
C  WITH RESPECT TO R AND Z.
C
         ZDSDR = (ZDRSDT * ZSINT + ZBND(J5,1) * ZCOST) / ZBND2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBND(J5,1) * ZSINT - ZDRSDT * ZCOST) / ZBND2
         ZDTDZ = ZCOST / ZRHO
C
         Z1 = (ZBND(J5,1) * ZCOS2T + ZDRSDT * ZSIN2T) / (ZRHO * ZBND2)
         Z2 = - ZCOS2T / ZRHO2
         Z3 =   ZSIN2T / ZRHO2
         Z4 = - ZSIN2T / ZRHO2
         Z5 = - ZSINT * ZCOST * ZFB / (ZRHO * ZBND3)
         Z6 =   ZSINT**2 * ZFB / (ZRHO * ZBND3)
         Z7 =   ZCOST**2 * ZFB / (ZRHO * ZBND3)
C
         ZDPDR(J5)  = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ(J5)  = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
         ZD2PRZ(J5) = ZD2PS2 * ZDSDR * ZDSDZ + ZD2PT2 * ZDTDR * ZDTDZ +
     &                ZD2PST * Z1 + ZDPDT * Z2 + ZDPDS * Z5
         ZD2PR2(J5) = ZD2PS2 * ZDSDR**2 + ZD2PT2 * ZDTDR**2 + 
     &                2. * ZD2PST * ZDSDR * ZDTDR + ZDPDT * Z3 +
     &                ZDPDS * Z6
         ZD2PZ2(J5) = ZD2PS2 * ZDSDZ**2 + ZD2PT2 * ZDTDZ**2 + 
     &                2. * ZD2PST * ZDSDZ * ZDTDZ + ZDPDT * Z4 +
     &                ZDPDS * Z7
C
 5       CONTINUE
C
         IMGAUS1 = NMGAUS + 1
C
         DO 6 J6=0,NMGAUS
C
         ZSURF(1) = KSPEN
         ZSURF(2) = J6
C
         WRITE(NPENN) (ZSURF(L),L=1,NCHI)
         WRITE(NPENN) (ZR(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZZ(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZDTTDR(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZDTTDZ(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZD2TTRZ(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZD2TTR2(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZD2TTZ2(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZDPDR(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZDPDZ(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZD2PRZ(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZD2PR2(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
         WRITE(NPENN) (ZD2PZ2(((L-1)*IMGAUS1)+J6+1),L=1,NCHI)
C
 6       CONTINUE
C
         RETURN
         END
C*DECK C2SM20
C*CALL PROCESS
         SUBROUTINE OUTXT(KP,PS)
C        #######################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM20 EVALUATE EQ'S FOR XTOR                                       *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMPLO.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IC(NPCHI1),        IS0(NPCHI1),       IT0(NPCHI1),
     R   ZBND(NPCHI1,5),    ZDBDS(NPCHI1,16),  
     R   ZDBDT(NPCHI1,16),  ZJAC(NPCHI1),      
     R   ZGSS(NPCHI1),      ZGST(NPCHI1),
     R   ZGTT(NPCHI1),      ZGPP(NPCHI1),
     R   ZPCEL(NPCHI1,16),  ZR(NPCHI1),   
     R   ZS(NPCHI1),        ZS1(NPCHI1),       
     R   ZS2(NPCHI1),       ZTETA(NPCHI1,5),   
     R   ZT(NPCHI1),        ZT1(NPCHI1),       
     R   ZT2(NPCHI1),       ZZ(NPCHI1)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ZEPS   = 1.E-3
         ZDPSIS = 2. * PS * CPSRF
C
         DO 1 J1=1,NTNOVA
           ZTETA(J1,1) = TETPEN(J1)
           ZTETA(J1,2) = TETPEN(J1) - 2. * ZEPS
           ZTETA(J1,3) = TETPEN(J1) -      ZEPS
           ZTETA(J1,4) = TETPEN(J1) +      ZEPS
 1         ZTETA(J1,5) = TETPEN(J1) + 2. * ZEPS
C
         CALL BOUND(NTNOVA,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NTNOVA,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NTNOVA,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NTNOVA,ZTETA(1,4),ZBND(1,4))
         CALL BOUND(NTNOVA,ZTETA(1,5),ZBND(1,5))
C
         CALL RESETI(IC,NTNOVA,1)
         DO 2 JT = 1,NT1
            DO 2 JG=1,NTNOVA
               IF (IC(JG).EQ.0) GOTO 2
               IT0(JG) = JT-1
               IF (TETPEN(JG).LE.CT(JT)) IC(JG)  = 0
 2       CONTINUE
         CALL RESETI(IC,NTNOVA,1)
         DO 3 JS = 1,NS1
            DO 3 JG=1,NTNOVA
               IF (IC(JG).EQ.0) GOTO 3
               IS0(JG) = JS-1
               IF (SIGPEN(JG).LE.CSIG(JS)) IC(JG)  = 0
 3       CONTINUE
C
         DO 4 J4=1,NTNOVA
           IF (IS0(J4) .GT. NS) IS0(J4) = NS
           IF (IS0(J4) .LT. 1)  IS0(J4) = 1
           IF (IT0(J4) .GT. NT) IT0(J4) = NT
           IF (IT0(J4) .LT. 1)  IT0(J4) = 1
C
           ZT(J4)  = TETPEN(J4)
           ZS(J4)  = SIGPEN(J4)
           ZS1(J4) = CSIG(IS0(J4))
           ZS2(J4) = CSIG(IS0(J4)+1)
           ZT1(J4) = CT(IT0(J4))
 4         ZT2(J4) = CT(IT0(J4)+1)
C
         CALL PSICEL(IS0,IT0,NTNOVA,NPCHI1,ZPCEL,CPSICL)
         CALL BASIS2(NTNOVA,NPCHI1,ZS1,ZS2,ZT1,ZT2,ZS,ZT,ZDBDS,ZDBDT)
C
         DO 5 J5=1,NTNOVA
C
         ZDRSDT = (ZBND(J5,2) + 8*(ZBND(J5,4) - ZBND(J5,3)) -
     -             ZBND(J5,5)) / (12. * ZEPS)
C
         ZDPDS = ZDBDS(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDS(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDS(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDS(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDS(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDS(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDS(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDS(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDS(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDS(J5,10) * ZPCEL(J5,10) +
     +           ZDBDS(J5,11) * ZPCEL(J5,11) +
     +           ZDBDS(J5,12) * ZPCEL(J5,12) +
     +           ZDBDS(J5,13) * ZPCEL(J5,13) +
     +           ZDBDS(J5,14) * ZPCEL(J5,14) +
     +           ZDBDS(J5,15) * ZPCEL(J5,15) +
     +           ZDBDS(J5,16) * ZPCEL(J5,16)
C
         ZDPDT = ZDBDT(J5, 1) * ZPCEL(J5, 1) +
     +           ZDBDT(J5, 2) * ZPCEL(J5, 2) +
     +           ZDBDT(J5, 3) * ZPCEL(J5, 3) +
     +           ZDBDT(J5, 4) * ZPCEL(J5, 4) +
     +           ZDBDT(J5, 5) * ZPCEL(J5, 5) +
     +           ZDBDT(J5, 6) * ZPCEL(J5, 6) +
     +           ZDBDT(J5, 7) * ZPCEL(J5, 7) +
     +           ZDBDT(J5, 8) * ZPCEL(J5, 8) +
     +           ZDBDT(J5, 9) * ZPCEL(J5, 9) +
     +           ZDBDT(J5,10) * ZPCEL(J5,10) +
     +           ZDBDT(J5,11) * ZPCEL(J5,11) +
     +           ZDBDT(J5,12) * ZPCEL(J5,12) +
     +           ZDBDT(J5,13) * ZPCEL(J5,13) +
     +           ZDBDT(J5,14) * ZPCEL(J5,14) +
     +           ZDBDT(J5,15) * ZPCEL(J5,15) +
     +           ZDBDT(J5,16) * ZPCEL(J5,16)
C
         ZCOST  = COS(ZTETA(J5,1))
         ZSINT  = SIN(ZTETA(J5,1))
C
         ZRHO   = SIGPEN(J5) * ZBND(J5,1)
         ZR(J5) = ZRHO * ZCOST + R0
         ZZ(J5) = ZRHO * ZSINT + RZ0
         ZBND2  = ZBND(J5,1)**2
C
C  COMPUTE FIRST DERIVATIVES OF THETA-TILD WITH
C  RESPECT TO R AND Z
C
         ZDZ    = ZZ(J5) - RZMAG
         ZDR    = ZR(J5) - RMAG
         ZDSUM  = ZDR**2 + ZDZ**2
C
         ZDTTDR  = - ZDZ / ZDSUM
         ZDTTDZ  =   ZDR / ZDSUM
C
C  COMPUTE ALL OTHER DERIVATIVES OF GEOMETRIC QUANTITIES
C  RELATED TO THE EQUILIBRIUM MESH AND INVOLVED IN THE
C  COMPUTATION OF THE FIRST DERIVATIVES OF PSI
C  WITH RESPECT TO R AND Z.
C
         ZDSDR = (ZDRSDT * ZSINT + ZBND(J5,1) * ZCOST) / ZBND2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBND(J5,1) * ZSINT - ZDRSDT * ZCOST) / ZBND2
         ZDTDZ = ZCOST / ZRHO
C
         ZDPDR  = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ  = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         ZJAC(J5)  = (ZDPDR * ZDTTDZ - ZDPDZ * ZDTTDR) / (ZR(J5)*ZDPSIS)
         ZGSS(J5) = (ZDPDR**2 + ZDPDZ**2) / ZDPSIS**2
         ZGTT(J5) = ZDTTDR**2 + ZDTTDZ**2
         ZGPP(J5) = 1./ ZR(J5)**2
         ZGST(J5) = (ZDPDR * ZDTTDR + ZDPDZ * ZDTTDZ) / ZDPSIS

 5       CONTINUE
C
         WRITE(NXTOR) PS,PSIISO(KP),CPR(KP),TMF(KP)
         WRITE(NXTOR) TTP(KP),CPPR(KP),ZDPSIS
         WRITE(NXTOR) (ZR(L),L=1,NTNOVA)
         WRITE(NXTOR) (ZZ(L),L=1,NTNOVA)
         WRITE(NXTOR) (ZJAC(L),L=1,NTNOVA)
         WRITE(NXTOR) (ZGSS(L),L=1,NTNOVA)
         WRITE(NXTOR) (ZGTT(L),L=1,NTNOVA)
         WRITE(NXTOR) (ZGPP(L),L=1,NTNOVA)
         WRITE(NXTOR) (ZGST(L),L=1,NTNOVA)
C
         RETURN
         END
C*DECK C2SM22
C*CALL PROCESS
         SUBROUTINE FOURFFT(KPSI,KMMAX)
C        ##############################
C
C                                        AUTHOR O. SAUTER, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM22 COMPUTE FAST FOURIER TRANSFORM OF THE DIFFERENT EQL, EQI     *
*        AND EQ3 TERMS                                                *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R     ICHIISO(2*NPCHI1),
     R     ZCHIFFT(2*NPCHI1), ZD2FUN(NPMGS*NTP1), ZWORK(NPMGS*NTP1,3),
     R     ZA(2*NPCHI1), ZB(2*NPCHI1), ZC(2*NPCHI1), ZD(2*NPCHI1),
     R     ZFFTEQL(2*NPCHI1,22), ZFFTEQI(2*NPCHI1,32),
     R     ZFFTEQ3(2*NPCHI1,17), ZWORK1(2*NPCHI1)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C     1. COMPUTE EQUIDISTANT CHI MESH, KMMAX INTERVALS
C        
         INCHI = 2*NCHI
C
         ZDCHI = 2. * CPI / DFLOAT(INCHI)
         DO I=1,INCHI+1
           ZCHIFFT(I) = DFLOAT(I-1)*ZDCHI
         ENDDO
C
C     2. PREPARE COEFFICIENTS FOR THE CUBIC SPLINE FIT DEPENDING
C     ONLY ON RELATIVE POSITION OF ZCHIFFT(I) WITH RESPECT TO CHIISO
C
         IGCHISO = NMGAUS*NT1
         CALL GCHI(KPSI)
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHIFFT(I)) - 1
C
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
           ZA(I) = (CHIISO(ICHISO+1) - ZCHIFFT(I)) / ZH
           ZB(I) = (ZCHIFFT(I) - CHIISO(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISO(ICHISO+1) - ZCHIFFT(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHIFFT(I) - CHIISO(ICHISO)) / 6.
C
         ENDDO
C
C     3. FOR EACH ARRAY: COMPUTE VALUES ON ZCHIFFT USING A PERIODIC
C     CUBIC SPLINE FIT AND COMPUTE FULL FOURIER TRANSFORM
C
C     EQL
C
         CALL SPLIFFT(CHIISO,EQL(1,1),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,1),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,2),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,2),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,3),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,3),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,4),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,4),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,5),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,5),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,6),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,6),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,7),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,7),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,8),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,8),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,9),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,9),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,10),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,10),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,11),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,11),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,12),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,12),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,13),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,13),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,14),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,14),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,15),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,15),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,16),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,16),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,17),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,17),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,18),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,18),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,19),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,19),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,20),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,20),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,21),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,21),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQL(1,22),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQL(1,22),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
C
C     EQI
C
         CALL SPLIFFT(CHIISO,EQI(1,1),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,1),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,2),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,2),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,3),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,3),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,4),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,4),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,5),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,5),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,6),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,6),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,7),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,7),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,8),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,8),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,9),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,9),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,10),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,10),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,11),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,11),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,12),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,12),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,13),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,13),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,14),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,14),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,15),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,15),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,16),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,16),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,17),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,17),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,18),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,18),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,19),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,19),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,20),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,20),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,21),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,21),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,22),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,22),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,23),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,23),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,24),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,24),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,25),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,25),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,26),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,26),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,27),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,27),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,28),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,28),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,29),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,29),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,30),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,30),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,31),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,31),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQI(1,32),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQI(1,32),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
C
C     EQ3
C
         CALL SPLIFFT(CHIISO,EQ3(1,1),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,1),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,2),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,2),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,3),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,3),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,4),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,4),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,5),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,5),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,6),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,6),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,7),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,7),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,8),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,8),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,9),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,9),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,10),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,10),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,11),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,11),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,12),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,12),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,13),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,13),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,14),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,14),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,15),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,15),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,16),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,16),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,EQ3(1,17),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTEQ3(1,17),INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
C-----------------------------------------------------------------------
C
C     4. COPY SPECTRUM INTO WANTED ARRAYS
C
C     F_M = CMPLX(FFT(M), -FFT(INCHI+1-M))
C     AS USE EXP(-I*K*X) INSTEAD OF EXP(I*K*X) IN C06FAF, MINUS
C     SIGN FOR IMAGINARY PART
C
         IF (MOD(KPSI,2) .EQ. 0) THEN
           IP = KPSI / 2 + 1
           DG11L(IP,1) =CMPLX(ZFFTEQL(1, 1),0.)
           DG22L(IP,1) =CMPLX(ZFFTEQL(1, 2),0.)
           DG33L(IP,1) =CMPLX(ZFFTEQL(1, 3),0.)
           DG12L(IP,1) =CMPLX(ZFFTEQL(1, 4),0.)
           JG11L(IP,1) =CMPLX(ZFFTEQL(1, 5),0.)
           JG22L(IP,1) =CMPLX(ZFFTEQL(1, 6),0.)
           JG33L(IP,1) =CMPLX(ZFFTEQL(1, 7),0.)
           JG12L(IP,1) =CMPLX(ZFFTEQL(1, 8),0.)
           JACOBI(IP,1)=CMPLX(ZFFTEQL(1, 9),0.)
           J2U(IP,1,1) =CMPLX(ZFFTEQL(1,10),0.)
           J3U(IP,1,1) =CMPLX(ZFFTEQL(1,11),0.)
           B2E(IP,1,1) =CMPLX(ZFFTEQL(1,12),0.)
           B3E(IP,1,1) =CMPLX(ZFFTEQL(1,13),0.)
           PEQ(IP,1,1) =CMPLX(ZFFTEQL(1,14),0.)
           DPEDS(IP,1,1)=CMPLX(ZFFTEQL(1,15),0.)
           GCHDZ(IP,1) =CMPLX(ZFFTEQL(1,16),0.)
           GSDZ (IP,1) =CMPLX(ZFFTEQL(1,17),0.)
           GBZ(IP,1)   =CMPLX(ZFFTEQL(1,18),0.)
           GBR(IP,1)   =CMPLX(ZFFTEQL(1,19),0.)
           JACOBINV(IP,1)=CMPLX(ZFFTEQL(1,20),0.)
           FRM(IP,1)   =CMPLX(ZFFTEQL(1,21),0.)
           FZM(IP,1)   =CMPLX(ZFFTEQL(1,22),0.)
C     EQI
           IDIY2(IP,1) =CMPLX(ZFFTEQI(1, 1),0.)
           IDIY3(IP,1) =CMPLX(ZFFTEQI(1, 2),0.)
           IG122(IP,1) =CMPLX(ZFFTEQI(1, 3),0.)
           IG123(IP,1) =CMPLX(ZFFTEQI(1, 4),0.)
           INXX(IP,1)  =CMPLX(ZFFTEQI(1, 5),0.)
           INXY(IP,1)  =CMPLX(ZFFTEQI(1, 6),0.)
           INYY(IP,1)  =CMPLX(ZFFTEQI(1, 7),0.)
           INZZ(IP,1)  =CMPLX(ZFFTEQI(1, 8),0.)
           IJ0QX(IP,1) =CMPLX(ZFFTEQI(1, 9),0.)
           IJ0QY(IP,1) =CMPLX(ZFFTEQI(1,10),0.)
           IGPX2(IP,1) =CMPLX(ZFFTEQI(1,11),0.)
           IGPX3(IP,1) =CMPLX(ZFFTEQI(1,12),0.)
           IGPY2(IP,1) =CMPLX(ZFFTEQI(1,13),0.)
           IGPY3(IP,1) =CMPLX(ZFFTEQI(1,14),0.)
           IDRXX(IP,1) =CMPLX(ZFFTEQI(1,15),0.)
           IRXZ(IP,1)  =CMPLX(ZFFTEQI(1,16),0.)
           IDRYX(IP,1) =CMPLX(ZFFTEQI(1,17),0.)
           IRYX(IP,1)  =CMPLX(ZFFTEQI(1,18),0.)
           IDRZX(IP,1) =CMPLX(ZFFTEQI(1,19),0.)
           IRZY(IP,1)  =CMPLX(ZFFTEQI(1,20),0.)
           VISXZ(IP,1) =CMPLX(ZFFTEQI(1,21),0.)
           VISYZ(IP,1) =CMPLX(ZFFTEQI(1,22),0.)
           IVS11(IP,1) =CMPLX(ZFFTEQI(1,23),0.)
           IVS12(IP,1) =CMPLX(ZFFTEQI(1,24),0.)
           IVS21(IP,1) =CMPLX(ZFFTEQI(1,25),0.)
           IVS22(IP,1) =CMPLX(ZFFTEQI(1,26),0.)
           GSFC(IP,1)  =CMPLX(ZFFTEQI(1,27),0.)
           GSCC(IP,1)  =CMPLX(ZFFTEQI(1,28),0.)
           GSFS(IP,1)  =CMPLX(ZFFTEQI(1,29),0.)
           GSCS(IP,1)  =CMPLX(ZFFTEQI(1,30),0.)
           GCFC(IP,1)  =CMPLX(ZFFTEQI(1,31),0.)
           GCFS(IP,1)  =CMPLX(ZFFTEQI(1,32),0.)
C     EQ3
           EQRHO(IP,1) =CMPLX(ZFFTEQ3(1, 1),0.)
           DRHOS(IP,1) =CMPLX(ZFFTEQ3(1, 2),0.)
           EQROT(IP,1) =CMPLX(ZFFTEQ3(1, 3),0.)
           DROT(IP,1)  =CMPLX(ZFFTEQ3(1, 4),0.)
           FEQ(IP,1)   =CMPLX(ZFFTEQ3(1, 5),0.)
           IWSQ1(IP,1) =CMPLX(ZFFTEQ3(1, 6),0.)
           IWSQ2(IP,1) =CMPLX(ZFFTEQ3(1, 7),0.)
           IWSQ3(IP,1) =CMPLX(ZFFTEQ3(1, 8),0.)
           IJ0QZ(IP,1) =CMPLX(ZFFTEQ3(1, 9),0.)
           JACOF(IP,1) =CMPLX(ZFFTEQ3(1,10),0.)
           B2F(IP,1)   =CMPLX(ZFFTEQ3(1,11),0.)
           B3F(IP,1)   =CMPLX(ZFFTEQ3(1,12),0.)
           JACOS(IP,1) =CMPLX(ZFFTEQ3(1,13),0.)
           IGF22(IP,1) =CMPLX(ZFFTEQ3(1,14),0.)
           B3FC(IP,1)  =CMPLX(ZFFTEQ3(1,15),0.)
           B2FC(IP,1)  =CMPLX(ZFFTEQ3(1,16),0.)
           DJCOF(IP,1) =CMPLX(ZFFTEQ3(1,17),0.)
C
           IF (KMMAX .EQ. NCHI) THEN
             IMPMAX = KMMAX-1
           ELSE IF (KMMAX .LT. NCHI) THEN
             IMPMAX = KMMAX
           ENDIF
           DO IMP1=2,IMPMAX

             IM2  = IMP1
             IM2P = INCHI+2-IMP1
C     EQL
             DG11L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 1),ZFFTEQL(IM2P, 1))
             DG22L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 2),ZFFTEQL(IM2P, 2))
             DG33L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 3),ZFFTEQL(IM2P, 3))
             DG12L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 4),ZFFTEQL(IM2P, 4))
             JG11L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 5),ZFFTEQL(IM2P, 5))
             JG22L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 6),ZFFTEQL(IM2P, 6))
             JG33L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 7),ZFFTEQL(IM2P, 7))
             JG12L(IP,IMP1) =CMPLX(ZFFTEQL(IM2, 8),ZFFTEQL(IM2P, 8))
             JACOBI(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 9),ZFFTEQL(IM2P, 9))
             J2U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,10),ZFFTEQL(IM2P,10))
             J3U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,11),ZFFTEQL(IM2P,11))
             B2E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,12),ZFFTEQL(IM2P,12))
             B3E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,13),ZFFTEQL(IM2P,13))
             PEQ(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,14),ZFFTEQL(IM2P,14))
             DPEDS(IP,IMP1,1)=CMPLX(ZFFTEQL(IM2,15),ZFFTEQL(IM2P,15))
             GCHDZ(IP,IMP1) =CMPLX(ZFFTEQL(IM2,16),ZFFTEQL(IM2P,16))
             GSDZ (IP,IMP1) =CMPLX(ZFFTEQL(IM2,17),ZFFTEQL(IM2P,17))
             GBZ(IP,IMP1)   =CMPLX(ZFFTEQL(IM2,18),ZFFTEQL(IM2P,18))
             GBR(IP,IMP1)   =CMPLX(ZFFTEQL(IM2,19),ZFFTEQL(IM2P,19))
             JACOBINV(IP,IMP1)=CMPLX(ZFFTEQL(IM2,20),ZFFTEQL(IM2P,20))
             FRM(IP,IMP1)   =CMPLX(ZFFTEQL(IM2,21),ZFFTEQL(IM2P,21))
             FZM(IP,IMP1)   =CMPLX(ZFFTEQL(IM2,22),ZFFTEQL(IM2P,22))
C     EQI
             IDIY2(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 1),ZFFTEQI(IM2P, 1))
             IDIY3(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 2),ZFFTEQI(IM2P, 2))
             IG122(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 3),ZFFTEQI(IM2P, 3))
             IG123(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 4),ZFFTEQI(IM2P, 4))
             INXX(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 5),ZFFTEQI(IM2P, 5))
             INXY(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 6),ZFFTEQI(IM2P, 6))
             INYY(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 7),ZFFTEQI(IM2P, 7))
             INZZ(IP,IMP1)  =CMPLX(ZFFTEQI(IM2, 8),ZFFTEQI(IM2P, 8))
             IJ0QX(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 9),ZFFTEQI(IM2P, 9))
             IJ0QY(IP,IMP1) =CMPLX(ZFFTEQI(IM2,10),ZFFTEQI(IM2P,10))
             IGPX2(IP,IMP1) =CMPLX(ZFFTEQI(IM2,11),ZFFTEQI(IM2P,11))
             IGPX3(IP,IMP1) =CMPLX(ZFFTEQI(IM2,12),ZFFTEQI(IM2P,12))
             IGPY2(IP,IMP1) =CMPLX(ZFFTEQI(IM2,13),ZFFTEQI(IM2P,13))
             IGPY3(IP,IMP1) =CMPLX(ZFFTEQI(IM2,14),ZFFTEQI(IM2P,14))
             IDRXX(IP,IMP1) =CMPLX(ZFFTEQI(IM2,15),ZFFTEQI(IM2P,15))
             IRXZ(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,16),ZFFTEQI(IM2P,16))
             IDRYX(IP,IMP1) =CMPLX(ZFFTEQI(IM2,17),ZFFTEQI(IM2P,17))
             IRYX(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,18),ZFFTEQI(IM2P,18))
             IDRZX(IP,IMP1) =CMPLX(ZFFTEQI(IM2,19),ZFFTEQI(IM2P,19))
             IRZY(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,20),ZFFTEQI(IM2P,20))
             VISXZ(IP,IMP1) =CMPLX(ZFFTEQI(IM2,21),ZFFTEQI(IM2P,21))
             VISYZ(IP,IMP1) =CMPLX(ZFFTEQI(IM2,22),ZFFTEQI(IM2P,22))
             IVS11(IP,IMP1) =CMPLX(ZFFTEQI(IM2,23),ZFFTEQI(IM2P,23))
             IVS12(IP,IMP1) =CMPLX(ZFFTEQI(IM2,24),ZFFTEQI(IM2P,24))
             IVS21(IP,IMP1) =CMPLX(ZFFTEQI(IM2,25),ZFFTEQI(IM2P,25))
             IVS22(IP,IMP1) =CMPLX(ZFFTEQI(IM2,26),ZFFTEQI(IM2P,26))
             GSFC(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,27),ZFFTEQI(IM2P,27))
             GSCC(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,28),ZFFTEQI(IM2P,28))
             GSFS(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,29),ZFFTEQI(IM2P,29))
             GSCS(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,30),ZFFTEQI(IM2P,30))
             GCFC(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,31),ZFFTEQI(IM2P,31))
             GCFS(IP,IMP1)  =CMPLX(ZFFTEQI(IM2,32),ZFFTEQI(IM2P,32))
C     EQ3
             EQRHO(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 1),ZFFTEQ3(IM2P, 1))
             DRHOS(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 2),ZFFTEQ3(IM2P, 2))
             EQROT(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 3),ZFFTEQ3(IM2P, 3))
             DROT(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2, 4),ZFFTEQ3(IM2P, 4))
             FEQ(IP,IMP1)   =CMPLX(ZFFTEQ3(IM2, 5),ZFFTEQ3(IM2P, 5))
             IWSQ1(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 6),ZFFTEQ3(IM2P, 6))
             IWSQ2(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 7),ZFFTEQ3(IM2P, 7))
             IWSQ3(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 8),ZFFTEQ3(IM2P, 8))
             IJ0QZ(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 9),ZFFTEQ3(IM2P, 9))
             JACOF(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,10),ZFFTEQ3(IM2P,10))
             B2F(IP,IMP1)   =CMPLX(ZFFTEQ3(IM2,11),ZFFTEQ3(IM2P,11))
             B3F(IP,IMP1)   =CMPLX(ZFFTEQ3(IM2,12),ZFFTEQ3(IM2P,12))
             JACOS(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,13),ZFFTEQ3(IM2P,13))
             IGF22(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,14),ZFFTEQ3(IM2P,14))
             B3FC(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,15),ZFFTEQ3(IM2P,15))
             B2FC(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,16),ZFFTEQ3(IM2P,16))
             DJCOF(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,17),ZFFTEQ3(IM2P,17))
           ENDDO
           IF (KMMAX.EQ.NCHI) THEN
             IMPMAX1 = IMPMAX+1
             DG11L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 1),0.)
             DG22L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 2),0.)
             DG33L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 3),0.)
             DG12L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 4),0.)
             JG11L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 5),0.)
             JG22L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 6),0.)
             JG33L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 7),0.)
             JG12L(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1, 8),0.)
             JACOBI(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 9),0.)
             J2U(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,10),0.)
             J3U(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,11),0.)
             B2E(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,12),0.)
             B3E(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,13),0.)
             PEQ(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,14),0.)
             DPEDS(IP,1,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1,15),0.)
             GCHDZ(IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,16),0.)
             GSDZ (IP,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,17),0.)
             GBZ(IP,IMPMAX1)   =CMPLX(ZFFTEQL(IMPMAX1,18),0.)
             GBR(IP,IMPMAX1)   =CMPLX(ZFFTEQL(IMPMAX1,19),0.)
             JACOBINV(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1,20),0.)
             FRM(IP,IMPMAX1)   =CMPLX(ZFFTEQL(IMPMAX1,21),0.)
             FZM(IP,IMPMAX1)   =CMPLX(ZFFTEQL(IMPMAX1,22),0.)
C     EQI
             IDIY2(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 1),0.)
             IDIY3(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 2),0.)
             IG122(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 3),0.)
             IG123(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 4),0.)
             INXX(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1, 5),0.)
             INXY(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1, 6),0.)
             INYY(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1, 7),0.)
             INZZ(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1, 8),0.)
             IJ0QX(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 9),0.)
             IJ0QY(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,10),0.)
             IGPX2(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,11),0.)
             IGPX3(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,12),0.)
             IGPY2(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,13),0.)
             IGPY3(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,14),0.)
             IDRXX(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,15),0.)
             IRXZ(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,16),0.)
             IDRYX(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,17),0.)
             IRYX(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,18),0.)
             IDRZX(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,19),0.)
             IRZY(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,20),0.)
             VISXZ(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,21),0.)
             VISYZ(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,22),0.)
             IVS11(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,23),0.)
             IVS12(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,24),0.)
             IVS21(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,25),0.)
             IVS22(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,26),0.)
             GSFC(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,27),0.)
             GSCC(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,28),0.)
             GSFS(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,29),0.)
             GSCS(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,30),0.)
             GCFC(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,31),0.)
             GCFS(IP,IMPMAX1)  =CMPLX(ZFFTEQI(IMPMAX1,32),0.)
C     EQ3
             EQRHO(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 1),0.)
             DRHOS(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 2),0.)
             EQROT(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 3),0.)
             DROT(IP,IMPMAX1)  =CMPLX(ZFFTEQ3(IMPMAX1, 4),0.)
             FEQ(IP,IMPMAX1)   =CMPLX(ZFFTEQ3(IMPMAX1, 5),0.)
             IWSQ1(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 6),0.)
             IWSQ2(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 7),0.)
             IWSQ3(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 8),0.)
             IJ0QZ(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 9),0.)
             JACOF(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1,10),0.)
             B2F(IP,IMPMAX1)   =CMPLX(ZFFTEQ3(IMPMAX1,11),0.)
             B3F(IP,IMPMAX1)   =CMPLX(ZFFTEQ3(IMPMAX1,12),0.)
             JACOS(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1,13),0.)
             IGF22(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1,14),0.)
             B3FC(IP,IMPMAX1)  =CMPLX(ZFFTEQ3(IMPMAX1,15),0.)
             B2FC(IP,IMPMAX1)  =CMPLX(ZFFTEQ3(IMPMAX1,16),0.)
             DJCOF(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1,17),0.)
           ENDIF
         ELSE IF (MOD(KPSI,2) .EQ. 1) THEN
           IP = (KPSI + 1) / 2
           DG11LM(IP,1)=CMPLX(ZFFTEQL(1, 1),0.)
           DG22LM(IP,1)=CMPLX(ZFFTEQL(1, 2),0.)
           DG33LM(IP,1)=CMPLX(ZFFTEQL(1, 3),0.)
           DG12LM(IP,1)=CMPLX(ZFFTEQL(1, 4),0.)
           JG11LM(IP,1)=CMPLX(ZFFTEQL(1, 5),0.)
           JG22LM(IP,1)=CMPLX(ZFFTEQL(1, 6),0.)
           JG33LM(IP,1)=CMPLX(ZFFTEQL(1, 7),0.)
           JG12LM(IP,1)=CMPLX(ZFFTEQL(1, 8),0.)
           JACOBM(IP,1)=CMPLX(ZFFTEQL(1, 9),0.)
           J2E(IP,1,1) =CMPLX(ZFFTEQL(1,10),0.)
           J3E(IP,1,1) =CMPLX(ZFFTEQL(1,11),0.)
           B2U(IP,1,1) =CMPLX(ZFFTEQL(1,12),0.)
           B3U(IP,1,1) =CMPLX(ZFFTEQL(1,13),0.)
           PRE(IP,1,1) =CMPLX(ZFFTEQL(1,14),0.)
           DPEDSM(IP,1,1)=CMPLX(ZFFTEQL(1,15),0.)
           GCHDZM(IP,1)=CMPLX(ZFFTEQL(1,16),0.)
           GSDZM (IP,1)=CMPLX(ZFFTEQL(1,17),0.)
           GBZM(IP,1)  =CMPLX(ZFFTEQL(1,18),0.)
           GBRM(IP,1)  =CMPLX(ZFFTEQL(1,19),0.)
C     EQI
           IDIY2M(IP,1)=CMPLX(ZFFTEQI(1, 1),0.)
           IDIY3M(IP,1)=CMPLX(ZFFTEQI(1, 2),0.)
           IG122M(IP,1)=CMPLX(ZFFTEQI(1, 3),0.)
           IG123M(IP,1)=CMPLX(ZFFTEQI(1, 4),0.)
           INXXM(IP,1) =CMPLX(ZFFTEQI(1, 5),0.)
           INXYM(IP,1) =CMPLX(ZFFTEQI(1, 6),0.)
           INYYM(IP,1) =CMPLX(ZFFTEQI(1, 7),0.)
           INZZM(IP,1) =CMPLX(ZFFTEQI(1, 8),0.)
           IJ0QXM(IP,1)=CMPLX(ZFFTEQI(1, 9),0.)
           IJ0QYM(IP,1)=CMPLX(ZFFTEQI(1,10),0.)
           IGPX2M(IP,1)=CMPLX(ZFFTEQI(1,11),0.)
           IGPX3M(IP,1)=CMPLX(ZFFTEQI(1,12),0.)
           IGPY2M(IP,1)=CMPLX(ZFFTEQI(1,13),0.)
           IGPY3M(IP,1)=CMPLX(ZFFTEQI(1,14),0.)
           IDRXXM(IP,1)=CMPLX(ZFFTEQI(1,15),0.)
           IRXZM(IP,1) =CMPLX(ZFFTEQI(1,16),0.)
           IDRYXM(IP,1)=CMPLX(ZFFTEQI(1,17),0.)
           IRYXM(IP,1) =CMPLX(ZFFTEQI(1,18),0.)
           IDRZXM(IP,1)=CMPLX(ZFFTEQI(1,19),0.)
           IRZYM(IP,1) =CMPLX(ZFFTEQI(1,20),0.)
           VISXZM(IP,1)=CMPLX(ZFFTEQI(1,21),0.)
           VISYZM(IP,1)=CMPLX(ZFFTEQI(1,22),0.)
           IVS11M(IP,1)=CMPLX(ZFFTEQI(1,23),0.)
           IVS12M(IP,1)=CMPLX(ZFFTEQI(1,24),0.)
           IVS21M(IP,1)=CMPLX(ZFFTEQI(1,25),0.)
           IVS22M(IP,1)=CMPLX(ZFFTEQI(1,26),0.)
           GSFCM(IP,1) =CMPLX(ZFFTEQI(1,27),0.)
           GSCCM(IP,1) =CMPLX(ZFFTEQI(1,28),0.)
           GSFSM(IP,1) =CMPLX(ZFFTEQI(1,29),0.)
           GSCSM(IP,1) =CMPLX(ZFFTEQI(1,30),0.)
           GCFCM(IP,1) =CMPLX(ZFFTEQI(1,31),0.)
           GCFSM(IP,1) =CMPLX(ZFFTEQI(1,32),0.)
C     EQ3
           EQRHOM(IP,1)=CMPLX(ZFFTEQ3(1, 1),0.)
           DRHOSM(IP,1)=CMPLX(ZFFTEQ3(1, 2),0.)
           EQROTM(IP,1)=CMPLX(ZFFTEQ3(1, 3),0.)
           DROTM(IP,1) =CMPLX(ZFFTEQ3(1, 4),0.)
           FEQM(IP,1)  =CMPLX(ZFFTEQ3(1, 5),0.)
           IWSQ1M(IP,1)=CMPLX(ZFFTEQ3(1, 6),0.)
           IWSQ2M(IP,1)=CMPLX(ZFFTEQ3(1, 7),0.)
           IWSQ3M(IP,1)=CMPLX(ZFFTEQ3(1, 8),0.)
           IJ0QZM(IP,1)=CMPLX(ZFFTEQ3(1, 9),0.)
           JACOFM(IP,1)=CMPLX(ZFFTEQ3(1,10),0.)
           B2FM(IP,1)  =CMPLX(ZFFTEQ3(1,11),0.)
           B3FM(IP,1)  =CMPLX(ZFFTEQ3(1,12),0.)
           JACOSM(IP,1)=CMPLX(ZFFTEQ3(1,13),0.)
           IGF22M(IP,1)=CMPLX(ZFFTEQ3(1,14),0.)
           B3FCM(IP,1) =CMPLX(ZFFTEQ3(1,15),0.)
           B2FCM(IP,1) =CMPLX(ZFFTEQ3(1,16),0.)
           DJCOFM(IP,1)=CMPLX(ZFFTEQ3(1,17),0.)
C
           IF (KMMAX .EQ. NCHI) THEN
             IMPMAX = KMMAX-1
           ELSE IF (KMMAX .LT. NCHI) THEN
             IMPMAX = KMMAX
           ENDIF
           DO IMP1=2,IMPMAX

             IM2  = IMP1
             IM2P = INCHI+2-IMP1
C     EQL
             DG11LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 1),ZFFTEQL(IM2P, 1))
             DG22LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 2),ZFFTEQL(IM2P, 2))
             DG33LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 3),ZFFTEQL(IM2P, 3))
             DG12LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 4),ZFFTEQL(IM2P, 4))
             JG11LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 5),ZFFTEQL(IM2P, 5))
             JG22LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 6),ZFFTEQL(IM2P, 6))
             JG33LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 7),ZFFTEQL(IM2P, 7))
             JG12LM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 8),ZFFTEQL(IM2P, 8))
             JACOBM(IP,IMP1)=CMPLX(ZFFTEQL(IM2, 9),ZFFTEQL(IM2P, 9))
             J2E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,10),ZFFTEQL(IM2P,10))
             J3E(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,11),ZFFTEQL(IM2P,11))
             B2U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,12),ZFFTEQL(IM2P,12))
             B3U(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,13),ZFFTEQL(IM2P,13))
             PRE(IP,IMP1,1) =CMPLX(ZFFTEQL(IM2,14),ZFFTEQL(IM2P,14))
             DPEDSM(IP,IMP1,1)=CMPLX(ZFFTEQL(IM2,15),ZFFTEQL(IM2P,15))
             GCHDZM(IP,IMP1)=CMPLX(ZFFTEQL(IM2,16),ZFFTEQL(IM2P,16))
             GSDZM (IP,IMP1)=CMPLX(ZFFTEQL(IM2,17),ZFFTEQL(IM2P,17))
             GBZM(IP,IMP1)  =CMPLX(ZFFTEQL(IM2,18),ZFFTEQL(IM2P,18))
             GBRM(IP,IMP1)  =CMPLX(ZFFTEQL(IM2,19),ZFFTEQL(IM2P,19))
C     EQI
             IDIY2M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 1),ZFFTEQI(IM2P, 1))
             IDIY3M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 2),ZFFTEQI(IM2P, 2))
             IG122M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 3),ZFFTEQI(IM2P, 3))
             IG123M(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 4),ZFFTEQI(IM2P, 4))
             INXXM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 5),ZFFTEQI(IM2P, 5))
             INXYM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 6),ZFFTEQI(IM2P, 6))
             INYYM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 7),ZFFTEQI(IM2P, 7))
             INZZM(IP,IMP1) =CMPLX(ZFFTEQI(IM2, 8),ZFFTEQI(IM2P, 8))
             IJ0QXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2, 9),ZFFTEQI(IM2P, 9))
             IJ0QYM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,10),ZFFTEQI(IM2P,10))
             IGPX2M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,11),ZFFTEQI(IM2P,11))
             IGPX3M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,12),ZFFTEQI(IM2P,12))
             IGPY2M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,13),ZFFTEQI(IM2P,13))
             IGPY3M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,14),ZFFTEQI(IM2P,14))
             IDRXXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,15),ZFFTEQI(IM2P,15))
             IRXZM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,16),ZFFTEQI(IM2P,16))
             IDRYXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,17),ZFFTEQI(IM2P,17))
             IRYXM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,18),ZFFTEQI(IM2P,18))
             IDRZXM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,19),ZFFTEQI(IM2P,19))
             IRZYM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,20),ZFFTEQI(IM2P,20))
             VISXZM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,21),ZFFTEQI(IM2P,21))
             VISYZM(IP,IMP1)=CMPLX(ZFFTEQI(IM2,22),ZFFTEQI(IM2P,22))
             IVS11M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,23),ZFFTEQI(IM2P,23))
             IVS12M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,24),ZFFTEQI(IM2P,24))
             IVS21M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,25),ZFFTEQI(IM2P,25))
             IVS22M(IP,IMP1)=CMPLX(ZFFTEQI(IM2,26),ZFFTEQI(IM2P,26))
             GSFCM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,27),ZFFTEQI(IM2P,27))
             GSCCM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,28),ZFFTEQI(IM2P,28))
             GSFSM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,29),ZFFTEQI(IM2P,29))
             GSCSM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,30),ZFFTEQI(IM2P,30))
             GCFCM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,31),ZFFTEQI(IM2P,31))
             GCFSM(IP,IMP1) =CMPLX(ZFFTEQI(IM2,32),ZFFTEQI(IM2P,32))
C     EQ3
             EQRHOM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 1),ZFFTEQ3(IM2P, 1))
             DRHOSM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 2),ZFFTEQ3(IM2P, 2))
             EQROTM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 3),ZFFTEQ3(IM2P, 3))
             DROTM(IP,IMP1) =CMPLX(ZFFTEQ3(IM2, 4),ZFFTEQ3(IM2P, 4))
             FEQM(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2, 5),ZFFTEQ3(IM2P, 5))
             IWSQ1M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 6),ZFFTEQ3(IM2P, 6))
             IWSQ2M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 7),ZFFTEQ3(IM2P, 7))
             IWSQ3M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 8),ZFFTEQ3(IM2P, 8))
             IJ0QZM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2, 9),ZFFTEQ3(IM2P, 9))
             JACOFM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,10),ZFFTEQ3(IM2P,10))
             B2FM(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,11),ZFFTEQ3(IM2P,11))
             B3FM(IP,IMP1)  =CMPLX(ZFFTEQ3(IM2,12),ZFFTEQ3(IM2P,12))
             JACOSM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,13),ZFFTEQ3(IM2P,13))
             IGF22M(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,14),ZFFTEQ3(IM2P,14))
             B3FCM(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,15),ZFFTEQ3(IM2P,15))
             B2FCM(IP,IMP1) =CMPLX(ZFFTEQ3(IM2,16),ZFFTEQ3(IM2P,16))
             DJCOFM(IP,IMP1)=CMPLX(ZFFTEQ3(IM2,17),ZFFTEQ3(IM2P,17))
           ENDDO
           IF (KMMAX.EQ.NCHI) THEN
             IMPMAX1 = IMPMAX+1
             DG11LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 1),0.)
             DG22LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 2),0.)
             DG33LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 3),0.)
             DG12LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 4),0.)
             JG11LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 5),0.)
             JG22LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 6),0.)
             JG33LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 7),0.)
             JG12LM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 8),0.)
             JACOBM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1, 9),0.)
             J2E(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,10),0.)
             J3E(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,11),0.)
             B2U(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,12),0.)
             B3U(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,13),0.)
             PRE(IP,1,IMPMAX1) =CMPLX(ZFFTEQL(IMPMAX1,14),0.)
             DPEDSM(IP,1,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1,15),0.)
             GCHDZM(IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1,16),0.)
             GSDZM (IP,IMPMAX1)=CMPLX(ZFFTEQL(IMPMAX1,17),0.)
             GBZM(IP,IMPMAX1)  =CMPLX(ZFFTEQL(IMPMAX1,18),0.)
             GBRM(IP,IMPMAX1)  =CMPLX(ZFFTEQL(IMPMAX1,19),0.)
C     EQI
             IDIY2M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1, 1),0.)
             IDIY3M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1, 2),0.)
             IG122M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1, 3),0.)
             IG123M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1, 4),0.)
             INXXM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 5),0.)
             INXYM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 6),0.)
             INYYM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 7),0.)
             INZZM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1, 8),0.)
             IJ0QXM(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1, 9),0.)
             IJ0QYM(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,10),0.)
             IGPX2M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,11),0.)
             IGPX3M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,12),0.)
             IGPY2M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,13),0.)
             IGPY3M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,14),0.)
             IDRXXM(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,15),0.)
             IRXZM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,16),0.)
             IDRYXM(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,17),0.)
             IRYXM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,18),0.)
             IDRZXM(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,19),0.)
             IRZYM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,20),0.)
             VISXZM(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,21),0.)
             VISYZM(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,22),0.)
             IVS11M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,23),0.)
             IVS12M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,24),0.)
             IVS21M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,25),0.)
             IVS22M(IP,IMPMAX1)=CMPLX(ZFFTEQI(IMPMAX1,26),0.)
             GSFCM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,27),0.)
             GSCCM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,28),0.)
             GSFSM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,29),0.)
             GSCSM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,30),0.)
             GCFCM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,31),0.)
             GCFSM(IP,IMPMAX1) =CMPLX(ZFFTEQI(IMPMAX1,32),0.)
C     EQ3
             EQRHOM(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1, 1),0.)
             DRHOSM(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1, 2),0.)
             EQROTM(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1, 3),0.)
             DROTM(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1, 4),0.)
             FEQM(IP,IMPMAX1)  =CMPLX(ZFFTEQ3(IMPMAX1, 5),0.)
             IWSQ1M(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1, 6),0.)
             IWSQ2M(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1, 7),0.)
             IWSQ3M(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1, 8),0.)
             IJ0QZM(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1, 9),0.)
             JACOFM(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1,10),0.)
             B2FM(IP,IMPMAX1)  =CMPLX(ZFFTEQ3(IMPMAX1,11),0.)
             B3FM(IP,IMPMAX1)  =CMPLX(ZFFTEQ3(IMPMAX1,12),0.)
             JACOSM(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1,13),0.)
             IGF22M(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1,14),0.)
             B3FCM(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1,15),0.)
             B2FCM(IP,IMPMAX1) =CMPLX(ZFFTEQ3(IMPMAX1,16),0.)
             DJCOFM(IP,IMPMAX1)=CMPLX(ZFFTEQ3(IMPMAX1,17),0.)
           ENDIF
         ENDIF
C
         RETURN
         END
C*DECK C2SM23
C*CALL PROCESS
         SUBROUTINE VACUFFT(KMMAX)
C        #########################
C
C                                        AUTHOR O. SAUTER, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM23  COMPUTE EQV TERMS FOR MARS USING FAST FFT TRANSFORM ON AN   *
*         2**L EQUIDISTANT CHI-MESH. CUBIC SPLINE INTERPOLATION OF THE*
*        TERMS FROM GAUSSIAN THETA-MESH TO EQUIDISTANT CHI-MESH.      *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
C
         DIMENSION
     R     ZBND(NPMGS*NTP1,4), ZS(2*NPV1), ZTETA(NPMGS*NTP1,4),
     R     ICHIISO(2*NPCHI1), ZCHIFFT(2*NPCHI1), ZD2FUN(NPMGS*NTP1),
     R     ZWORK(NPMGS*NTP1,3),ZWORK1(2*NPCHI1),
     R     ZA(2*NPCHI1), ZB(2*NPCHI1), ZC(2*NPCHI1), ZD(2*NPCHI1),
     R     ZFFTG11(2*NPCHI1), ZFFTG22(2*NPCHI1), ZFFTG33(2*NPCHI1),
     R     ZFFTG12(2*NPCHI1), ZFFTJAC(2*NPCHI1),
     R     ZG11(NPMGS*NTP1), ZG22(NPMGS*NTP1), ZG33(NPMGS*NTP1),
     R     ZG12(NPMGS*NTP1), ZJAC1(NPMGS*NTP1)
C
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IEDGE = 2*NPSI
         ZEPS = 1.E-3
C
C     PREPARE COEFFICIENTS FOR FFT AND CUBIC SPLINE
C        
         INCHI = 2*NCHI
C
         ZDCHI = 2. * CPI / DFLOAT(INCHI)
         DO I=1,INCHI+1
           ZCHIFFT(I) = DFLOAT(I-1)*ZDCHI
         ENDDO
C
         IGCHISO = NMGAUS*NT1
C     ASSUMES THAT CHIISO HAS ALREADY ARRAY AT BOUNDARY (FROM CALL FOURFFT)
C     (OTHERWISE, ADD CALL GCHI(2*NPSI))
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHIFFT(I)) - 1
C
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
           ZA(I) = (CHIISO(ICHISO+1) - ZCHIFFT(I)) / ZH
           ZB(I) = (ZCHIFFT(I) - CHIISO(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISO(ICHISO+1) - ZCHIFFT(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHIFFT(I) - CHIISO(ICHISO)) / 6.
C
         ENDDO
C
C     PRECOMPUTE SOME ARRAYS
C
         DO 1 J1=1,NV
C
         ZS(2*(J1-1)+1) = CSV(J1)
         ZS(2*J1      ) = CSMV(J1)
C
    1    CONTINUE
C
         ZS(2*NV+1) = CSV(NV1)
C
         DO 2 J2=1,NMGAUS*NT1
C
         ZTETA(J2,1) = TETPSI(J2,IEDGE) - 2. * ZEPS
         ZTETA(J2,2) = TETPSI(J2,IEDGE) -      ZEPS
         ZTETA(J2,3) = TETPSI(J2,IEDGE) +      ZEPS
         ZTETA(J2,4) = TETPSI(J2,IEDGE) + 2. * ZEPS
C
    2    CONTINUE
C
         CALL BOUND(NMGAUS*NT1,ZTETA(1,1),ZBND(1,1))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,2),ZBND(1,2))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,3),ZBND(1,3))
         CALL BOUND(NMGAUS*NT1,ZTETA(1,4),ZBND(1,4))
C
C     PRECOMPUTE TERMS NON-DEPENDING ON VACUUM S VALUE
C
         DO 3 J3=1,NMGAUS*NT1
C
         ZCOST = COS(TETPSI(J3,IEDGE))
         ZSINT = SIN(TETPSI(J3,IEDGE))
         ZBNDT = BNDISO(J3,IEDGE)
         ZRHO  = ZBNDT
         ZR    = RRISO(J3,IEDGE)
         ZZ    = RZISO(J3,IEDGE)
         ZGP   = GPISO(J3,IEDGE)
         ZDPDS = DPSISO(J3,IEDGE)
         ZDPDT = 0.
C
         ZDRSDT = (ZBND(J3,1) + 8*(ZBND(J3,3) - ZBND(J3,2)) -
     -             ZBND(J3,4)) / (12. * ZEPS)
         ZJAC   = CP(IEDGE) * ZR**NER * ZGP**NEGP
C
         ZDCDT = ZRHO * ZBNDT * ZR / (ZJAC * ZDPDS)
C
         ZDSDR = (ZDRSDT * ZSINT + ZBNDT * ZCOST) / ZBNDT**2
         ZDTDR = - ZSINT / ZRHO
         ZDSDZ = (ZBNDT * ZSINT - ZDRSDT * ZCOST) / ZBNDT**2
         ZDTDZ = ZCOST / ZRHO
C
         ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         ZDRDC = - ZJAC * ZDPDZ / ZR
         ZDZDC = ZJAC * ZDPDR / ZR
C
         ZRW = ZR - R0W
         ZZW = ZZ - RZ0W
         ZG11(J3) = ZRW**2 + ZZW**2
         ZG22(J3) = ZDRDC**2 + ZDZDC**2
         ZG12(J3) = ZRW * ZDRDC + ZZW * ZDZDC
         ZG33(J3) = ZRW
         ZJAC1(J3)= ZRW * ZDZDC - ZZW * ZDRDC
C
    3    CONTINUE
C
C     AT EACH S VACUUM LOCATION: COMPUTE TERMS ON GAUSSIAN MESH, INTERPOLATE
C     THEM ON EQUIDISTANT CHI-MESH AND COMPUTE FOURIER TRANSFORM
C
         DO 9 J9=1,2*NV+1
           ZS1   = ZS(J9)
C
C     USE ZBND AND ZTETA MEMORY SPACE
         DO 8 J8=1,NT1*NMGAUS
C
           ZTETA(J8,1) = ZS1 * (R0W + ZS1 * ZG33(J8)) * ZJAC1(J8)
           ZBND (J8,1) = ZG11(J8) / ZTETA(J8,1)
           ZBND (J8,2) = ZS1**2 * ZG22(J8) / ZTETA(J8,1)
           ZBND (J8,3) = (R0W + ZS1 * ZG33(J8))**2 / ZTETA(J8,1)
           ZBND (J8,4) = ZS1  * ZG12(J8) / ZTETA(J8,1)
C
    8    CONTINUE
C
C     CUBIC SPLINE AND FFT
C
         CALL SPLIFFT(CHIISO,ZBND (1,1),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTG11,INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,ZBND (1,2),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTG22,INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,ZBND (1,3),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTG33,INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
         CALL SPLIFFT(CHIISO,ZBND (1,4),IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZFFTG12,INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
CC         CALL SPLIFFT(CHIISO,ZTETA(1,1),IGCHISO,RC2PI,ZD2FUN,
CC     +     ZWORK,ZFFTJAC,INCHI,ICHIISO,ZA,ZB,ZC,ZD,ZWORK1)
C
C     COPY SPECTRUM FOR M=0,INCHI
C
         IF (MOD(J9,2) .EQ. 1) THEN
           IV = (J9 + 1) / 2
           DG11LV(IV,1) = CMPLX(ZFFTG11(1),0.) 
           DG22LV(IV,1) = CMPLX(ZFFTG22(1),0.) 
           DG33LV(IV,1) = CMPLX(ZFFTG33(1),0.) 
           DG12LV(IV,1) = CMPLX(ZFFTG12(1),0.) 
CC           VJACOB(IV,1) = CMPLX(ZFFTJAC(1),0.) 
C
           IF (KMMAX .EQ. NCHI) THEN
             IMPMAX = KMMAX-1
           ELSE IF (KMMAX .LT. NCHI) THEN
             IMPMAX = KMMAX
           ENDIF
           DO IMP1=2,IMPMAX

             IM2  = IMP1
             IM2P = INCHI+2-IMP1
C
             DG11LV(IV,IMP1) = CMPLX(ZFFTG11(IM2),ZFFTG11(IM2P)) 
             DG22LV(IV,IMP1) = CMPLX(ZFFTG22(IM2),ZFFTG22(IM2P)) 
             DG33LV(IV,IMP1) = CMPLX(ZFFTG33(IM2),ZFFTG33(IM2P)) 
             DG12LV(IV,IMP1) = CMPLX(ZFFTG12(IM2),ZFFTG12(IM2P)) 
CC             VJACOB(IV,IMP1) = CMPLX(ZFFTJAC(IM2),ZFFTJAC(IM2P)) 
           ENDDO
           IF (KMMAX.EQ.NCHI) THEN
             IMPMAX1 = IMPMAX+1
             DG11LV(IV,IMPMAX1) = CMPLX(ZFFTG11(IMPMAX1),0.) 
             DG22LV(IV,IMPMAX1) = CMPLX(ZFFTG22(IMPMAX1),0.) 
             DG33LV(IV,IMPMAX1) = CMPLX(ZFFTG33(IMPMAX1),0.) 
             DG12LV(IV,IMPMAX1) = CMPLX(ZFFTG12(IMPMAX1),0.) 
CC             VJACOB(IV,IMPMAX1) = CMPLX(ZFFTJAC(IMPMAX1),0.) 
           ENDIF
         ELSE IF (MOD(J9,2) .EQ. 0) THEN
           IV = J9 / 2
           DG11LMV(IV,1) = CMPLX(ZFFTG11(1),0.) 
           DG22LMV(IV,1) = CMPLX(ZFFTG22(1),0.) 
           DG33LMV(IV,1) = CMPLX(ZFFTG33(1),0.) 
           DG12LMV(IV,1) = CMPLX(ZFFTG12(1),0.) 
CC           VJACOM(IV,1)  = CMPLX(ZFFTJAC(1),0.) 
C
           IF (KMMAX .EQ. NCHI) THEN
             IMPMAX = KMMAX-1
           ELSE IF (KMMAX .LT. NCHI) THEN
             IMPMAX = KMMAX
           ENDIF
           DO IMP1=2,IMPMAX

             IM2  = IMP1
             IM2P = INCHI+2-IMP1
C
             DG11LMV(IV,IMP1) = CMPLX(ZFFTG11(IM2),ZFFTG11(IM2P)) 
             DG22LMV(IV,IMP1) = CMPLX(ZFFTG22(IM2),ZFFTG22(IM2P)) 
             DG33LMV(IV,IMP1) = CMPLX(ZFFTG33(IM2),ZFFTG33(IM2P)) 
             DG12LMV(IV,IMP1) = CMPLX(ZFFTG12(IM2),ZFFTG12(IM2P)) 
CC             VJACOM(IV,IMP1)  = CMPLX(ZFFTJAC(IM2),ZFFTJAC(IM2P)) 
           ENDDO
           IF (KMMAX.EQ.NCHI) THEN
             IMPMAX1 = IMPMAX+1
             DG11LMV(IV,IMPMAX1) = CMPLX(ZFFTG11(IMPMAX1),0.) 
             DG22LMV(IV,IMPMAX1) = CMPLX(ZFFTG22(IMPMAX1),0.) 
             DG33LMV(IV,IMPMAX1) = CMPLX(ZFFTG33(IMPMAX1),0.) 
             DG12LMV(IV,IMPMAX1) = CMPLX(ZFFTG12(IMPMAX1),0.) 
CC       VJACOM(IV,IMPMAX1)  = CMPLX(ZFFTJAC(IMPMAX1),0.) 
           ENDIF
         ENDIF
C
    9    CONTINUE
C
         RETURN
         END
C*DECK C2SM24
C*CALL PROCESS
         SUBROUTINE SPLIFFT(PCHI,PFUN,KCHI,PPERIOD,PD2FUN,PWORK,PFFTFUN,
     +                      KNFFT,KINDEX,PA,PB,PC,PD,PWORK1)
C        ###############################################################
C
C                                        AUTHOR O. SAUTER, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM24 COMPUTE CUBIC SPLINE WITH PERIODIC B.C. AND FAST FOURIER     *
*        TRANSFORM ASSUMING PCHI REAL                                 *
*                                                                     *
***********************************************************************
C
         DIMENSION PCHI(KCHI), PFUN(KCHI), PD2FUN(KCHI), PWORK(KCHI,3),
     +     PFFTFUN(KNFFT), KINDEX(KNFFT),PWORK1(KNFFT),
     +     PA(KNFFT), PB(KNFFT), PC(KNFFT), PD(KNFFT)
c.......................................................................
C
C     COMPUTE CUBIC SPLINE OF PFUN USING PA, PB, PC, PD PRECOMPUTED
C
         CALL SPLCY(PCHI,PFUN,KCHI,PPERIOD,PD2FUN,PWORK(1,1),PWORK(1,2),
     +     PWORK(1,3))
         DO I=1,KNFFT
           PFFTFUN(I) = PA(I)*PFUN(KINDEX(I)) +
     +       PB(I)*PFUN(KINDEX(I)+1) +
     +       PC(I)*PD2FUN(KINDEX(I)) + PD(I)*PD2FUN(KINDEX(I)+1)
         ENDDO
C
C     COMPUTE FAST FOURIER TRANSFORM WITH NAG
C
         IFAIL = 0
         CALL C06FAE(PFFTFUN,KNFFT,PWORK1,IFAIL)
         IF (IFAIL.NE.0.) THEN
           PRINT*,' IFAIL = ',IFAIL,' IN CO6EAF'
           STOP
         ENDIF
         CALL SSCAL(KNFFT,1./SQRT(DFLOAT(KNFFT)),PFFTFUN,1)
C
         RETURN
         END
C*DECK C2SJ01
C*CALL PROCESS
         SUBROUTINE CURENT(KN,PPSI,PR,PJIPHI)
C        ####################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SI02 COMPUTE CURRENT DENSITY PROFILE FOR PPSI(I), I=1,KN          *
*        (SEE EQ. (3), EQ. (10) AND EQ. (11) OF PUBLICATION)          *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     I   I1(NPT+2*NPISO),    IC(NPT+2*NPISO),
     R   PJIPHI(KN),         PPSI(KN),           PR(KN),
     R   ZCID2(NPT+2*NPISO),
     R   ZCID0(NPT+2*NPISO), ZPPRIM(NPT+2*NPISO),ZFUNC(NPT+2*NPISO), 
     R   ZTMF(NPT+2*NPISO),  ZFUNCD(NPT+2*NPISO) 
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IF (NSURF .EQ. 1) GOTO 300
C
         IF (NPROFZ .EQ. 0) THEN
C
            CALL PPRIME(KN,PPSI,ZPPRIM)
            CALL PRFUNC(KN,PPSI,ZFUNC,0)
CYQL2018
            IF (NSTTP.EQ.4) CALL PRFUNC(KN,PPSI,ZFUNCD,1)
C
         ENDIF         
C 
         IF (NSTTP .EQ. 1 .AND. NPROFZ .EQ. 0) THEN
C
            DO 1 J1=1,KN
C
            PJIPHI(J1) = - PR(J1) * ZPPRIM(J1) - ZFUNC(J1) / PR(J1)
C
    1       CONTINUE
C
         ELSE IF (NSTTP .NE. 1 .OR. NPROFZ .EQ. 1) THEN
C
C   BRACKET OUT [PS(I); PS(I+1)] INTERVAL SUCH THAT
C   PSIISO(I) <= PPSI(J) <= PSIISO(I+1), J=1,...,KN
C   WHEN I* OR <J . B> IS SPECIFIED
C
C
           CALL RESETI(IC,KN,1)
           DO 2 JS = 1,NISO
             DO 2 JG=1,KN
               IF (IC(JG).EQ.0) GOTO 2
               ZS = 1. - PPSI(JG) / SPSIM
               IF (ZS .LT. 0.) ZS = 0.
               I1(JG) = JS-1
               IF (SQRT(ZS).LE.CSIPR(JS)) IC(JG) = 0
 2          CONTINUE
C
C           WRITE(*,*) 'YQL2018-CURENT: ZS JPHI'

            DO 3 J3=1,KN
C
            IF (I1(J3) .LT. 1)      I1(J3) = 1
            IF (I1(J3) .GT. NISO-1) I1(J3) = NISO - 1
C
            ZS = 1. - PPSI(J3) / SPSIM
C
            IF (ZS .LT. 0.) ZS = 0.
C
            ZS = SQRT(ZS)
C
            ZH = CSIPR(I1(J3)+1) - CSIPR(I1(J3))
            ZA = (CSIPR(I1(J3)+1) - ZS) / ZH
            ZB = (ZS - CSIPR(I1(J3))) / ZH
            ZC = (ZA + 1) * (ZA - 1) * ZH * (CSIPR(I1(J3)+1) - ZS) / 6.
            ZD = (ZB + 1) * (ZB - 1) * ZH * (ZS - CSIPR(I1(J3))) / 6.
C
            IF (NPROFZ .EQ. 1) THEN
C
               ZPPRIM(J3) = ZA * CPPR(I1(J3))   + ZB * CPPR(I1(J3)+1) +
     +                      ZC * D2CPPR(I1(J3)) + ZD * D2CPPR(I1(J3)+1)
               ZFUNC(J3)  = ZA * TTP(I1(J3))    + ZB * TTP(I1(J3)+1) +
     +                      ZC * D2TTP(I1(J3))  + ZD * D2TTP(I1(J3)+1)
               PJIPHI(J3) = - PR(J3) * ZPPRIM(J3) - ZFUNC(J3) / PR(J3)
C
            ELSE IF (NPROFZ .EQ. 0) THEN
C 
               ZCID0(J3) = ZA * CID0(I1(J3))   + ZB * CID0(I1(J3)+1) +
     +                     ZC * D2CID0(I1(J3)) + ZD * D2CID0(I1(J3)+1)
               ZCID2(J3) = ZA * CID2(I1(J3))   + ZB * CID2(I1(J3)+1) +
     +                     ZC * D2CID2(I1(J3)) + ZD * D2CID2(I1(J3)+1)
C
               IF (NSTTP .EQ. 2) THEN
C
                  PJIPHI(J3) = ZFUNC(J3) * ZCID0(J3) / PR(J3) +
     +                       ZPPRIM(J3) * (ZCID2(J3) / PR(J3) - PR(J3))
C
               ELSE IF (NSTTP .EQ. 3) THEN
C
                  ZTMF(J3)  = ZA*TMF(I1(J3))   + ZB*TMF(I1(J3)+1) +
     +                          ZC*D2TMF(I1(J3)) + ZD*D2TMF(I1(J3)+1)
C
                  ZX         = 1. + ZCID2(J3) / ZTMF(J3)**2
                  PJIPHI(J3) = ZFUNC(J3)/(ZX*PR(J3))+ZPPRIM(J3) * 
     *                           (ZCID0(J3)/(ZX*PR(J3))-PR(J3))
CYQL2018
               ELSE IF (NSTTP .EQ. 4) THEN
C                 ZX=D(CID2)/D(ZS)
                  ZX      =-1./ZH*CID2(I1(J3))  + 1./ZH*CID2(I1(J3)+1) +
     +                     (1./3.-ZA*ZA)*ZH/2.*D2CID2(I1(J3)) +
     +                     (ZB*ZB-1./3.)*ZH/2.*D2CID2(I1(J3)+1) 
                  PJIPHI(J3) =-ZPPRIM(J3)*PR(J3)-2.*CPI*CPI*ZFUNC(J3)
     &                         /SPSIM*ZCID2(J3)/ZS/PR(J3)*(
     &                         ZFUNCD(J3)*ZCID2(J3) +
     &                         ZFUNC(J3)*ZX )
C                 IF (ZS.LT.0.02) WRITE(*,*)   ZS,PJIPHI(J3)
               ENDIF
            ENDIF
C
    3       CONTINUE

C
         ENDIF
C
         RETURN
C
  300    CONTINUE
C
***********************************************************************
*                                                                     *
*  PROFILES FOR SOLOVEV CASE                                          *
*                                                                     *
***********************************************************************
C
         DO 301 J301=1,KN
C
         PJIPHI(J301) = - PR(J301) * CPP
C
  301    CONTINUE
C
         RETURN
         END
C*DECK C2SJ02
C*CALL PROCESS
         SUBROUTINE ISOFUN(KN)
C        #####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SJ02  COMPUTES T, T-TPRIME, PRESSURE AND P-PRIME PROFILES AT     *
*          CONSTANT FLUX SURFACES                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   ZDNEDP(2*(NPISO+1)), ZDTEDP(2*(NPISO+1)),
     R   ZCPR(2*(NPISO+2)),   ZCPPR(2*(NPISO+2)),
     R   ZD2PPR(2*(NPISO+2)), ZPISO(2*(NPISO+2)),  
     R   ZS(2*(NPISO+2)),     ZWORK(2*(NPISO+2)),
     R   ztemp(2*(npiso+2)),  ZWORK1(2*(NPISO+2))
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IF (NSURF .EQ. 1) GOTO 300
C
         IF (NPROFZ .EQ. 0) THEN
C
            CALL PPRIME(KN,PSIISO,ZCPPR)
C
***********************************************************************
*                                                                     *
*  INTEGRATE P' BY CUBIC SPLINE QUADRATURE                            *
*                                                                     *
***********************************************************************
C
            DO 1 J1=1,KN
C
            ZS(J1)    = SQRT(1 - PSIISO(J1) / SPSIM)
            ZPISO(J1) = PSIISO(J1)
C
    1       CONTINUE
C
            IN = KN
C
            IF (ZS(IN) .NE. 1.) THEN
C
               IN        = KN + 1
               ZS(IN)    = 1.
               ZPISO(IN) = 0.
C
               CALL PPRIME(1,ZPISO(IN),ZCPPR(IN))
C
            ENDIF
C
            CALL SPLINE(ZPISO,ZCPPR,IN,ZD2PPR,ZWORK,ZWORK1)
            ztemp(1) = zd2ppr(1)
            ztemp(in) = zd2ppr(in)
            do 100 j100 = 2,in-1
               ztemp(j100) = zd2ppr(j100)
               if (abs(zd2ppr(j100)).lt.1.e4) goto 100
               if (zd2ppr(j100-1).ge.0..and.zd2ppr(j100).ge.0.
     &             .and.zd2ppr(j100+1).ge.0.) goto 100
               if (zd2ppr(j100-1).le.0..and.zd2ppr(j100).le.0.
     &             .and.zd2ppr(j100+1).le.0.) goto 100
               write(*,*) ' warning inaccurate integration, i =',j100
 100        continue
            if (ztemp(in-1)*ztemp(in).le.0.) ztemp(in) = 0.
            if (ztemp(1)*ztemp(2).le.0.) ztemp(1) = 0.
            do 101 j101 = 1,in
 101           zd2ppr(j101) = ztemp(j101)
C
C
C           STANDARD VERSION
            ZCPR(IN) = PREDGE
C      
C           SPECIAL VERSION TO RECOVER EXPONENTIAL P-PROFILE
C           ZCPR(IN) = PREDGE*(-ZCPPR(1))*(ZPISO(IN)-ZPISO(1))
C
C
            DO 2 J2=IN-1,1,-1
C
            J2P1 = J2 + 1
C
            ZPD =    - (ZPISO(J2P1) - ZPISO(J2)) *
     &              (.5*(ZCPPR(J2) + ZCPPR(J2P1)) - 
     &              (ZD2PPR(J2) + ZD2PPR(J2P1)) *
CMSC     &              (ZPISO(J2P1) - ZPISO(J2))**2/24.)
     &              (ZPISO(J2P1) - ZPISO(J2))**3/48.)
cab
            IF (ZDP.LT.0.) THEN
              WRITE(*,*) 'NEGATIVE PRESSURE INCREMENT, J =',J2
              ZDP = 0.
            ENDIF
            ZCPR(J2) = ZCPR(J2P1) + ZPD
C
    2       CONTINUE
C
            CALL SCOPY(KN,ZCPR,1,CPR,1)
            CALL SCOPY(KN,ZCPPR,1,CPPR,1)
            CALL SPLINE(ZS,CPR,KN,D2CPR,ZWORK,ZWORK1)
            CALL SPLINE(ZS,CPPR,IN,D2CPPR,ZWORK,ZWORK1)
C
         ELSE IF (NPROFZ .EQ. 1) THEN
C
            CALL POLYNM(KN,AP,DENSTY,1)
            CALL POLYNM(KN,AT,TEMPER,1)
            CALL POLYNM(KN,AP,ZDNEDP,2)
            CALL POLYNM(KN,AT,ZDTEDP,2)
            CALL SSCAL(KN,SCALNE,DENSTY,1)
            CALL SSCAL(KN,SCALNE,ZDNEDP,1)
C
            DO 3 J3=1,KN
C
            CPR(J3)  = DENSTY(J3) * TEMPER(J3)
            CPPR(J3) = ZDNEDP(J3) * TEMPER(J3) +
     +                 DENSTY(J3) * ZDTEDP(J3)
C
            IF (CPPR(J3) .GT. 0.) CPPR(J3) = 0.
C
    3       CONTINUE
C
            CALL SPLINE(ZS,CPPR,KN,D2CPPR,ZWORK,ZWORK1)
            CALL SPLINE(ZS,CPR,KN,D2CPR,ZWORK,ZWORK1)
C
         ENDIF
C
         IF (NSTTP .EQ. 1) THEN
C
            CALL PRFUNC(KN,PSIISO,TTP,0)
C
         ELSE IF (NSTTP .EQ. 2) THEN
C
            IF (NPROFZ .EQ. 0) THEN
C
               CALL PRFUNC(KN,PSIISO,CIPR,0)
C
            ELSE IF (NPROFZ .EQ. 1) THEN
C
               DO 4 J4=1,KN
C
               CIPR(J4) = TEMPER(J4) * SQRT(TEMPER(J4))
C
               IF (CIPR(J4) .LT. 0.) CIPR(J4) = 0.
C
    4          CONTINUE
C 
            ENDIF
C
            DO 5 J5=1,KN
C
            TTP(J5) = - CIPR(J5) * CID0(J5) -
     +                  CPPR(J5) * CID2(J5)
            ZS(J5)  = SQRT(1 - PSIISO(J5) / SPSIM)
C
    5       CONTINUE

C           WRITE(*,*) 'YQL2018-2: ZS TTP NSTTP,SPSIM=',NSTTP,SPSIM
C           DO J6=1,KN
C              WRITE(*,*) ZS(J6),TTP(J6)
C           ENDDO 
C
         ELSE IF (NSTTP .EQ. 3) THEN
C
            CALL PRFUNC(KN,PSIISO,CIPR,0)
C
            DO 6 J6=1,KN
C
            TTP(J6) = - (CIPR(J6) + CPPR(J6) * CID0(J6)) / 
     /                  (1. + CID2(J6) / TMF(J6)**2)
            ZS(J6)  = SQRT(1 - PSIISO(J6) / SPSIM)
C
    6       CONTINUE
C
CYQL2018: ADD NSTTP=4 OPTION, WHERE (P',Q) IS SPECIFIED AS INPUT IN EXPEQ
         ELSE IF (NSTTP .EQ. 4) THEN   

            CALL PRFUNC(KN,PSIISO,CIPR,0)
 
            DO J6=1,KN
               ZS(J6)  = SQRT(1 - PSIISO(J6) / SPSIM)
            ENDDO
            CALL SPLINE(ZS,CID2,KN,D2CID2,ZWORK,ZWORK1)
            CALL PRFUNC(KN,PSIISO,ZWORK,1)

            DO J6=1,KN
C              ZX=D(CID2)/D(ZS)
               IF (J6.LT.KN) THEN
                  ZH = ZS(J6+1)-ZS(J6)
                  ZX = (CID2(J6+1)-CID2(J6))/ZH -
     &                 (D2CID2(J6)/3.+D2CID2(J6+1)/6.)*ZH
               ELSE 
                  ZH = ZS(J6)-ZS(J6-1)
                  ZX = (CID2(J6)-CID2(J6-1))/ZH + 
     &                 (D2CID2(J6-1)/6.+D2CID2(J6)/3.)*ZH
               ENDIF
               TTP(J6)=-2.*CPI*CPI*CIPR(J6)/SPSIM*CID2(J6)/ZS(J6)*(
     &                 CIPR(J6)*ZX + CID2(J6)*ZWORK(J6) )
            ENDDO

C           WRITE(*,*) 'YQL2018-4: ZS TTP NSTTP,SPSIM=',NSTTP,SPSIM
C           DO J6=1,KN
C              WRITE(*,*) ZS(J6),TTP(J6)
C           ENDDO 
         ENDIF
C
         CALL SPLINE(ZS,CID0,KN,D2CID0,ZWORK,ZWORK1)
         CALL SPLINE(ZS,CID2,KN,D2CID2,ZWORK,ZWORK1)
C
***********************************************************************
*                                                                     *
*  INTEGRATE T T' BY CUBIC SPLINE QUADRATURE                          *
*                                                                     *
***********************************************************************
C
         CALL SPLINE(PSIISO,TTP,KN,D2TTP,ZWORK,ZWORK1)
C
         IF (NTMF0.EQ.0) THEN
C
            TMF(KN) = 0.5
            IF (NSTTP.EQ.4) TMF(KN)=2.*(CPI*CIPR(KN)/CIDQ(KN))**2
C
            DO 7 J7=KN-1,1,-1
C
            J7P1 = J7 + 1
C
            TMF(J7) = TMF(J7P1) - (PSIISO(J7P1) - PSIISO(J7)) *
     &                            (.5*(TTP(J7) + TTP(J7P1)) - 
     &                            (D2TTP(J7) + D2TTP(J7P1)) *
CMSC     &                            (PSIISO(J7P1) - PSIISO(J7))**2 / 24.)
     &                            (PSIISO(J7P1) - PSIISO(J7))**3 / 48.)
CMSC            WRITE(6,'("J8,PSIISO,TMF,TTP,D2TTP",I5,1p4e12.4)')J8,
CMSC     &           PSIISO(J8),TMF(J8),TTP(J8),D2TTP(J8)
C
    7       CONTINUE
C
         ELSE
C
            TMF(1) = 0.5
CMSC ADDED 11/25/06
            If (NRFP.EQ.1) TMF(1)=AT(3)**2/2.
C
            DO 8 J8=1,KN-1
C
            J8P1 = J8 + 1
C
            TMF(J8P1) = TMF(J8) + (PSIISO(J8P1) - PSIISO(J8)) *
CMSC     &                            (.5*(TTP(J8) + TTP(J8P1)) -
     &                            (.5*(TTP(J8) + TTP(J8P1)) +
     &                             (D2TTP(J8) + D2TTP(J8P1)) *
CMSC     &                             (PSIISO(J8P1) - PSIISO(J8))**2 / 24.)
     &                             (PSIISO(J8P1) - PSIISO(J8))**3 / 48.)
C
    8       CONTINUE
C
         ENDIF
C
C           WRITE(*,1001) (J,ZS(J),TMF(J),CIPR(J),CPPR(J),PSIISO(J),
C    &                     TTP(J),J = 1,KN)
C1001 FORMAT('    J        CS  T**2/2  CIPR  CPPR  PSIISO  TTP',/
C    &      ,(I5,6E12.4))
C

         IF (NRFP.NE.1) THEN
         DO 9 J9=1,KN
C
         IF (TMF(J9).LT.0) THEN
C
            PRINT*,'TMF(',J9,')**2 NEGATIVE. READ MESSAGE PRECEEDING',
     &             'DO LOOP 8 IN SUBROUTINE ISOFUN'
            WRITE(*,1000) (J,ZS(J),TMF(J),TTP(J),D2TTP(J),J = 1,KN)
 1000 FORMAT('    J        CS      T**2          TTP      D2TTP',/
     &      ,(I5,2F10.4,2E12.4))
            WRITE(*,1002) (J,ZS(J),TMF(J),CID0(J),CID2(J),J = 1,KN)
 1002 FORMAT('    J        CS  T**2     CID0       CID2',/
     &      ,(I5,2F10.4,2E12.4))

            STOP
C
         ENDIF
C
         TMF(J9) = SQRT(2. * TMF(J9))
C
    9    CONTINUE
         ELSE
         DO J9=1,KN
CMSC           TMF(J9) = 2.*TMF(J9)
CLIU       IF (TMF(J9).LT.0.) THEN
CLIU          WRITE(6,'("CHANGE INPUT PARAMETER, NO  SOLUTION")')
CLIU          STOP
CLIU       END IF
CLIU       TMF(J9) = SQRT(2.*TMF(J9))
CMSC           IF (TMF(J9).LT.0) THEN
CMSC             TMF(J9) = -SQRT(-TMF(J9))
CLIU       IF (TTP(J9).GE.0.) TMF(J9)=-TMF(J9)
CMSC           ELSE
CMSC             TMF(J9) = SQRT(TMF(J9))
CMSC           ENDIF  
CLIU       NUMERICAL INTEGRATION ABOVE FROM TTP RESULTS IN BAD BEHAVIOR
CLIU       FOR T AROUND FIELD REVERSAL POINT, DUE TO THE SPLINE FOR 
CLIU       T**2/2. INSTEAD WE USE ANALYTICAL FORMULA FOR T HERE
           ZS1 = 1.-PSIISO(J9)/SPSIM
           TMF(J9)=AT(3)+AT(1)*(-ZS1+ZS1**(AT(2)+1)*AT(4)/(AT(2)+1))
         ENDDO  
         ENDIF
C
         CALL SPLINE(ZS,TMF,KN,D2TMF,ZWORK,ZWORK1)
         CALL SPLINE(ZS,TTP,KN,D2TTP,ZWORK,ZWORK1)

CYQL2018 
C        COMPUTE Q
C        WRITE(*,*) 'YQL2018-2: ZS  T Q'
C        DO J9=1,KN
C           WRITE(*,*) ZS(J9), TMF(J9), TMF(J9)*CIDQ(J9)/CPI/2.
C        ENDDO
C
         RETURN
C
  300    CONTINUE
C
***********************************************************************
*                                                                     *
*  PROFILES FOR SOLOVEV CASE                                          *
*                                                                     *
***********************************************************************
C
         DO 301 J301=1,KN
C
         CPR(J301)  = PSIISO(J301) * CPP
         CPPR(J301) = CPP
         TMF(J301)  = 1.
         TTP(J301)  = 0.
C
  301    CONTINUE
C
         RETURN
         END
C*DECK C2SI01
C*CALL PROCESS
         SUBROUTINE PRFUNC(KN,PP,PT,KD)
C        ###########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SI01  EVALUATE FUNCTIONAL FORM OF TT', I* OR I_PARALLEL          *
*          NFUNC = 1 -----> POLYNOMIAL (NIPR; AT(1:NSOUR))            *
*          NFUNC = 2 -----> POL. 3 SECT. (NIPR; AT,AT2,AT3(1:7); AT4(1:3))
*          NFUNC = 3 -----> PRINCETON PROFILE DEFINITION (AT(1:8))    *
*          NFUNC = 4 -----> EXPERIMENTAL DATA (RFUNC(1:NPPF+1))       *
*          NFUNC = 5 -----> COMPLICATED, UNCOMMENTED                  *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R   PP(KN),    PT(KN),   ZS(2*NPISO+NPT)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(PT,KN)
C
         IF (NRFP .EQ. 1) GOTO 30
C
         IF (NFUNC .EQ. 1) THEN
C
***********************************************************************
*                                                                     *
* FUNC GIVEN AS POLYNOMIAL IN PSI/PSIMIN OF DEGREE NSOUR              *
*                                                                     *
***********************************************************************
C
            CALL RESETR(PT,KN,AT(NSOUR))
C
            DO 2 J2 = NSOUR-1,1,-1
C
            DO 1 J1 = 1, KN
C
            IF (NIPR .EQ. 1) THEN
               ZS1 = PP(J1) / SPSIM
            ELSE IF (NIPR .EQ. 2) THEN
               ZS1 = SQRT(PP(J1) / SPSIM)
            ELSE IF (NIPR .EQ. 3) THEN
               ZS1 = (PP(J1) / SPSIM)**(1./3.)
            ELSE IF (NIPR .EQ. 4) THEN
               ZS1 = (PP(J1) / SPSIM)**(1./4.)
            ENDIF
C             
            PT(J1) = PT(J1) * ZS1 + AT(J2)
C
    1       CONTINUE
    2       CONTINUE
C
         ELSE IF (NFUNC .EQ. 2) THEN
C
***********************************************************************
*                                                                     *
*  FUNC IS GIVEN AS A SUM OF NIPR PROFILES                            *
*                                                                     *
*  FOR S=0     TO S=AT(1)  FUNC(PSI) = (ZA0+ZA1*S+ZA2*S**2)           *
*  FOR S=AT(1) TO S=AT(2)  FUNC(PSI) = (ZB0+ZB1*S+ZB2*S**2+ZB3*S**3)  *
*  FOR S=AT(2) TO S=1      FUNC(PSI) = (ZC0+ZC1*S)                    *
*                                                                     *
*         FOR FIRST PROFILE  S = S1 = (1-PSI/PSIM)                    *
*             SECOND PROFILE S = S2 = SQRT(S1)                        *
*             THIRD PROFILE  S = S3 = S1**0.25                        *
*                                                                     *
*  NIPR = 4 : A FOURTH PROFILE IS ADDED                               *
*                                                                     *
*  FUNC(PSI) = AT4(3) * EXP(-((S1- AT4(1)) / AT4(2))**2)              *
*                                                                     *
*  WITH AT4(1)  : CENTRE IN S1                                        *
*  WITH AT4(2)  : WIDTH IN S1                                         *
*  WITH AT4(3)  : HEIGHT                                              *
*                                                                     *
***********************************************************************
*                                                                     *
*  FIRST PROFILE :                                                    *
*                                                                     *
***********************************************************************
C
            DO 3 J3=1,KN
C
            ZS(J3) = 1. - PP(J3) / SPSIM
C
    3       CONTINUE
C
            CALL ATCOEF(KN,ZS,AT,PT,1)
C
            IF (NIPR .EQ. 1) RETURN
C
***********************************************************************
*                                                                     *
*  SECOND PROFILE :                                                   *
*                                                                     *
***********************************************************************
C
            DO 4 J4=1,KN
C
            ZS1 = 1. - PP(J4) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            ZS(J4) = SQRT(ZS1)
C
   4        CONTINUE
C
            CALL ATCOEF(KN,ZS,AT2,PT,1)
C
            IF (NIPR .EQ. 2) RETURN
C
***********************************************************************
*                                                                     *
*  THIRD  PROFILE :                                                   *
*                                                                     *
***********************************************************************
C
            DO 5 J5=1,KN
C
            ZS1 = 1. - PP(J5) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            ZS(J5) = (ZS1)**0.25
C
   5        CONTINUE
C
            CALL ATCOEF(KN,ZS,AT3,PT,1)
C
            IF (NIPR .EQ. 3) RETURN
C
***********************************************************************
*                                                                     *
*  FOURTH PROFILE:                                                    *
*                                                                     *
***********************************************************************
C
            IF (AT4(2) .LE. 0.) RETURN
C
            DO 6 J6=1,KN
C
            ZS(J6) = 1. - PP(J6) / SPSIM
C
            IF (ZS(J6) .LT. 0.) ZS(J6) = 0.
C
    6       CONTINUE
C
            CALL ATCOEF(KN,ZS,AT4,PT,2)
C
         ELSE IF (NFUNC .EQ. 3) THEN
C
***********************************************************************
*                                                                     *
* "PRINCETON" PROFILE WITH TWO EXPONENTS                              *
*  see Manickam et al., Phys.Plasmas 1 p.1601, May 1994               *
*                                                                     *
*  FUNC = AT(1) * (1 - S**AT(3))**AT(2) +                             *
*         AT(1)*AT(4)*AT(6)**2*S*(1-S**AT(8))**AT(7)/                 *
*                               ((S-AT(5))**2+AT(6)**2)               *
*                                                                     *
*         WHERE S = 1 - PSI / SPSIM                                   *
*                                                                     *
***********************************************************************
C
            DO 7 J7=1,KN
C
            ZS1 = 1. - PP(J7) / SPSIM
C
            IF (ZS1 .LT. EPSMCH) ZS1 = EPSMCH
            ZARG1 = 1. - ZS1**at(3)
            if (ZARG1.lt.EPSMCH) ZARG1 = EPSMCH
            ZARG2 = 1. - ZS1**at(8)
            if (ZARG2.lt.EPSMCH) ZARG2 = EPSMCH
            ZARG3 = ZS1 - at(5)
            ZARG4 = 1. - (ZS1/AT(12))**at(11)
            if (ZARG4.lt.EPSMCH) ZARG4 = EPSMCH
C
            PTEMP1 = AT(1) * ZARG1**AT(2)
            PTEMP2 = AT(6)**2*ZS1 * ZARG2**AT(7)
            PTEMP3 = ZARG3 * ZARG3 + AT(6) * AT(6)
            PTEMP2 = AT(4) * AT(1) * PTEMP2 / PTEMP3
C
            PT(J7) = PTEMP1 + PTEMP2
C
    7       CONTINUE
C
         ELSE IF (NFUNC .EQ. 4) THEN
C
***********************************************************************
*                                                                     *
*  INTERPOLATE FUNC WITH CUBIC SPLINES ON RFUN VALUES (IF FUNC IS     *
*  GIVEN BY A SET OF POINTS)                                          *
*                                                                     *
***********************************************************************
C
            CALL PPSPLN(KN,PP,NPPF,FCSM,RFUN,D2RFUN,PT,KD)
C
         ELSE IF (NFUNC .EQ. 5) THEN
C
***********************************************************************
*                                                                     *
*  2 CUBICS                                                           *
*                                                                     *
***********************************************************************
C
         ZA0 = 2.*(1.-AT(2))/AT(1)+AT(3)
         ZA1 =  -3.*(1.-AT(2))/AT(1)-AT(3)
C
         ZB0 = AT(3)-2.*AT(2)/(AT(1)-1.)
         ZB1 = -AT(3)*(AT(1)+2.)+3.*AT(2)*(AT(1)+1.)/(AT(1)-1.)
         ZB2 = AT(3)*(2.*AT(1)+1.)-6.*AT(1)*AT(2)/(AT(1)-1.)
         ZB3 = -AT(1)*AT(3)+AT(2)*(3.*AT(1)-1.)/(AT(1)-1.)
C
         DO 9 J9=1,KN
C
         ZS1 = 1. - PP(J9) / SPSIM
C
         IF (ZS1 .LT. 0.) ZS1 = 0.
C
         IF (NIPR .EQ. 2) THEN
            ZS1 = SQRT(ZS1)
         ELSE IF (NIPR .EQ. 3) THEN
            ZS1 = ZS1**0.25
         ELSE IF (NIPR .EQ. 4) THEN
            ZS1 = ZS1**2
         ENDIF
C
         IF (ZS1 .LE. AT(1)) THEN
            PT(J9) = 1. + ZS1**2 * (ZA1 + ZA0 * ZS1 / AT(1))/AT(1)
         ELSE 
            PT(J9) = (ZB3+ZS1*(ZB2+ZS1*(ZB1+ZS1*ZB0)))/(AT(1)-1.)**2
         ENDIF
C
    9    CONTINUE
C
         ENDIF
C
         RETURN
C
  30     CONTINUE
C
***********************************************************************
*                                                                     *
* REVERSED FIELD PINCH EQUILIBRIUM :                                  *
*                                                                     *
* COMPUTE TTPRIME(PSI) FOR THE RFP ASSUMING THE PROFILES              *
* DETERMINED BY MU(0)=AT(1) , ALFA=AT(2) , Bz(1)=AT(3), g=AT(4)       *
* MU=MU(0)*(1-g*psi**ALFA)=T'                                         *
* ZS1 is normalized psi, the radius is normalized to the major radius *
* The total poloidal flux is abs(psism)                               *
***********************************************************************
C
         ZALF  = AT(2) 
         ZALF1 = ZALF + 1.
C
         DO 31 J31=1,KN
C
         ZS1 = 1. - PP(J31) / SPSIM
C
         IF (ZS1 .LT. 0.) ZS1 = 0.
C
CLiu     ZTT      = AT(3) + AT(1) * (- ZS1 + (1./ZALF1) * ZS1**ZALF1)
CLiu     PT(J31) = ZTT * AT(1) * (1. - ZS1**ZALF) / SPSIM
CMSC     The following modified 11/26/06
CLIU     ZTT     = AT(3)+SPSIM*AT(1)*(ZS1 - AT(4)/ZALF1*ZS1**ZALF1)
CLIU     PT(J31) = -ZTT * AT(1) * (1. - AT(4)*ZS1**ZALF)
CLIU     modified again to fit Roberto's version 22/08/07    
         ZTT     = AT(3)-AT(1)*(ZS1 - AT(4)/ZALF1*ZS1**ZALF1)
         PT(J31) = ZTT * AT(1) * (1. - AT(4)*ZS1**ZALF) /SPSIM
C
  31     CONTINUE
C
         RETURN
         END
C*DECK C2SI02
C*CALL PROCESS
         SUBROUTINE ATCOEF(KN,PX,PAT,PT,K)
C        #################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SI02  EVALUATE TT', I* OR I_PARALLEL IF GIVEN AS POLYNOMIALS IN  *
*         SEVERALS SECTIONS                                           *
*                                                                     *
***********************************************************************
Cab      AT(1) = T1         FIRST INTERVAL BOUNDARY
Cab      AT(2) = T2         SECOND INTERVAL BOUNDARY
Cab      CURRENT PROFILE IS QUADRATIC IN [0,T1], CUBIC IN [T1,T2],
Cab                         LINEAR IN [T2,1]
Cab      AT(3) = F(0)
Cab      AT(4) = F'(0)
Cab      AT(5) = F'(T1)
Cab      AT(6) = F(1)
Cab      AT(7) = F'(T2) = F'(1)
Cab
Cab      Routine has been modified to remove problems when a grid
Cab      point coincides with interval boundaries
C**********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PAT(*),   PT(KN),   PX(KN)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IF (K .EQ. 1) THEN
C
         CALL COPYAT(PAT,ZA0,ZA1,ZA2,ZB0,ZB1,ZB2,ZB3,ZC0,ZC1)
C
         DO 1 J1=1,KN
C
         ZX = PX(J1)
         IF (ZX .LE. PAT(1))
     &   PT(J1) = PT(J1) +  (ZA0 + ZX * (ZA1 + ZX * ZA2))
         IF (PAT(1) .LT. ZX .AND. ZX .LE. PAT(2))
     &   PT(J1) = PT(J1) +  ZB0 + ZX*(ZB1 + ZX*(ZB2 + ZX*ZB3))
         IF (PAT(2) .LT. ZX)
     &   PT(J1) = PT(J1) + ZC0 + ZC1 * ZX
C
    1    CONTINUE
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ELSE IF (K .EQ. 2) THEN
C
         DO 2 J2=1,KN
C
         ZEXP = ((PX(J2) - PAT(1)) / PAT(2))**2
C
         IF (ZEXP .LT. 100.) THEN
C            
            PT(J2) = PT(J2) + PAT(3) * EXP(-ZEXP)
C
         ENDIF
C
    2    CONTINUE
C
         ENDIF
C
         RETURN
         END
C*DECK C2SI03
C*CALL PROCESS
         SUBROUTINE COPYAT(PAT,PA0,PA1,PA2,PB0,PB1,PB2,PB3,PC0,PC1)
C        ##########################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SI03  EVALUATE POLYNOMIAL COEFFICIENTS IF TT', I* OR I_PARALLEL  *
*          ARE PRESCRIBED AS POLYNOMIALS IN SEVERALS SECTIONS         *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION  
     R   PAT(*)
C
         INCLUDE 'CUCDCD.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
***********************************************************************
*                                                                     *
*  FIRST SECTION                                                      *
*                                                                     *
***********************************************************************
         IF (PAT(1) .NE. 0.) THEN
C
           PA0 = PAT(3)
           PA1 = PAT(4)
           PA2 = .5 * (PAT(5) - PAT(4)) / PAT(1)
C
         ELSE IF (PAT(1).EQ.0.0) THEN
C
           PA0 = PAT(3)
           PA1 = 0.
           PA2 = 0.
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  SECOND SECTION                                                     *
*                                                                     *
***********************************************************************
C
         IF (PAT(1) .NE. PAT(2)) THEN
C
           ZF1 = PA0 + PA1 * PAT(1) + PA2 * PAT(1) * PAT(1)
           ZF2 = PAT(6) - PAT(7) + PAT(7) * PAT(2)
           PB3 = FC3(PAT(1),ZF1,PAT(5),PAT(2),ZF2,PAT(7))
           PB2 = FC2(PAT(1),ZF1,PAT(5),PAT(2),ZF2,PAT(7))
           PB1 = FC1(PAT(1),ZF1,PAT(5),PAT(2),ZF2,PAT(7))
           PB0 = FC0(PAT(1),ZF1,PAT(5),PAT(2),ZF2,PAT(7))
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  THIRD SECTION                                                      *
*                                                                     *
***********************************************************************
C
         PC1 = PAT(7)
         PC0 = PAT(6) - PAT(7)
C
         RETURN
         END
C*DECK C2SP01
C*CALL PROCESS
         SUBROUTINE PPRIME(KN,PP,PT)
C        ###########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP01  COMPUTE P-PRIME(PSI)                                       *
*          NPPFUN = 1  -----> POLYNOMIAL (AP(1:NSOUR))                *
*          NPPFUN = 2  -----> POL. IN 3 SECTIONS (AP AND/OR AP2(1:7)) *
*          NPPFUN = 3  -----> PRINCETON DEFINITION (AP OR AP2(1:6))   *
*          NPPFUN = 4  -----> EXPERIMENTAL DATA (RPPF(1:NPPF))        *
*          NPPFUN = 5  -----> GAUSSIAN (AP(1:3))                      *
*          NPPFUN = 6  -----> 6 SECTIONS (AP AND/OR AP2(1:13))        *
*          NPPFUN = 7  -----> VIA PRESSURE (AP OR AP2(1:3))           *
*                                                                     *
***********************************************************************
C
C     THE PROFILES CAN BE RESCALED WITH CPRESS. IN THIS CASE, THE
C     PARAMETERS SHOULD BE MODIFIED IN AUXVAL AND CPRESS SET BACK TO 1.0
C     AS ONE MODIFIES THE INPUT PARAMETERS. CPRESSO KEEPS INPUT CPRESS
C
C-----------------------------------------------------------------------
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   PP(KN),   PT(KN),    ZS(2*NPISO+NPT)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         INCLUDE 'CUCCCC.inc'
         INCLUDE 'CUCDCD.inc'
         INCLUDE 'QUAQDQ.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(PT,KN)
C
         IF (NBLOPT .GE. 2 .OR. NBSOPT .EQ. 2) GOTO 11
         CALL PPRIM0(KN,PP,PT)
         RETURN
C
   11    CONTINUE
         CALL PPSPLN(KN,PP,NPPR,PCSM,RPRM,D2RPRM,PT,0)
         RETURN
         END
C
         SUBROUTINE PPRIM0(KN,PP,PT)
C        ###########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP01  COMPUTE P-PRIME(PSI)                                       *
*          NPPFUN = 1  -----> POLYNOMIAL (AP(1:NSOUR))                *
*          NPPFUN = 2  -----> POL. IN 3 SECTIONS (AP AND/OR AP2(1:7)) *
*          NPPFUN = 3  -----> PRINCETON DEFINITION (AP OR AP2(1:6))   *
*          NPPFUN = 4  -----> EXPERIMENTAL DATA (RPPF(1:NPPF))        *
*          NPPFUN = 5  -----> GAUSSIAN (AP(1:3))                      *
*          NPPFUN = 6  -----> 6 SECTIONS (AP AND/OR AP2(1:13))        *
*          NPPFUN = 7  -----> VIA PRESSURE (AP OR AP2(1:3))           *
*                                                                     *
***********************************************************************
C
C     THE PROFILES CAN BE RESCALED WITH CPRESS. IN THIS CASE, THE
C     PARAMETERS SHOULD BE MODIFIED IN AUXVAL AND CPRESS SET BACK TO 1.0
C     AS ONE MODIFIES THE INPUT PARAMETERS. CPRESSO KEEPS INPUT CPRESS
C
C-----------------------------------------------------------------------
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   PP(KN),   PT(KN),    ZS(2*NPISO+NPT)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         INCLUDE 'CUCCCC.inc'
         INCLUDE 'CUCDCD.inc'
         INCLUDE 'QUAQDQ.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(PT,KN)
C
         IF (NPPFUN .EQ. 1) THEN
C
***********************************************************************
*                                                                     *
* P-PRIME(PSI) DEFINED AS POLYNOMIAL IN PSI/PSIMIN OF DEGREE NSOUR    *
*                                                                     *
***********************************************************************
C
            CALL RESETR(PT,KN,AP(NSOUR))
C
            DO 2 J2=NSOUR-1,1,-1
C
            DO 1 J1=1, KN
C
            ZS1    = PP(J1) / SPSIM
            PT(J1) = PT(J1) * ZS1 + AP(J2)
C
    1       CONTINUE
    2       CONTINUE
C
         ELSE IF (NPPFUN .EQ. 2) THEN
C
***********************************************************************
*                                                                     *
*  NPP =1 OR 2                                                        *
*                                                                     *
*   PPRIME IS GIVEN AS A SUM OF NPP PROFILES                          *
*                                                                     *
*  FOR S=0     TO S=AP(1)  PPRIME(PSI)=-(ZA0+ZA1*S)                   *
*  FOR S=AP(1) TO S=AP(2)  PPRIME(PSI)=-(ZB0+ZB1*S+ZB2*S**2+ZB3*S**3) *
*  FOR S=AP(2) TO S=1      PPRIME(PSI)=-(ZC0+ZC1*S+ZC2*S**2)          *
*                                                                     *
*         FOR FIRST PROFILE  S = S1 = (1-PSI/PSIM)                    *
*             SECOND PROFILE S = S2 = SQRT(S1)                        *
*                                                                     *
***********************************************************************
C
***********************************************************************
*                                                                     *
*  FIRST PROFILE                                                      *
*                                                                     *
***********************************************************************
C
            DO 3 J3=1,KN
C
            ZS(J3) = 1. - PP(J3) / SPSIM
C
            IF (ZS(J3) .LT. 0.) ZS(J3) = 0.
C
   3        CONTINUE
C
            CALL APCOEF(KN,ZS,AP,PT)
C
            IF (NPP .EQ. 1) RETURN
C
***********************************************************************
*                                                                     *
*  SECOND PROFILE                                                     *
*                                                                     *
***********************************************************************
C
            DO 4 J4=1,KN
C
            ZS1 = 1. - PP(J4) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            ZS(J4) = SQRT(ZS1)
C
   4        CONTINUE
C
            CALL APCOEF(KN,ZS,AP2,PT)
C
         ELSE IF (NPPFUN .EQ. 3) THEN
C
***********************************************************************
*                                                                     *
*     USE EITHER A OR AP2 DEPENDING ON NPP=1 OR 2                     *
*                                                                     *
*  P-PRIME(PSI) = AP(1)*(1-S**AP(3))**AP(2) + AP(4)*S**AP(5)          *
*  IF AP(6) .GT. 0. THEN:                                             *
*  P-PRIME(PSI) = AP(1)*(1-S**AP(3))**AP(2) +                         *
*               + AP(4)*S**AP(5)*(1-EXP((S-1)/AP(6)))                 *
*                                                                     *
*  WHERE  S = 1 - PSI / SPSIM                                         *
*                                                                     *
***********************************************************************
C
            ZEPS = 1.E-12
C
            DO 5 J5=1,KN
C
            ZS1 = 1. - PP(J5) / SPSIM
C
            IF (ZS1 .LT. ZEPS) ZS1 = ZEPS
C
            ZPT1 = AP(1) * (1. - ZS1**AP(3))**AP(2)
            ZPT2 = AP(4) * ZS1**AP(5)
            IF (AP(6).GT.0.) ZPT2 = ZPT2 * (1. - EXP((ZS1-1.)/AP(6)))
            PT(J5) = ZPT1 + ZPT2
C
   5        CONTINUE
C
            IF (NPP .EQ. 1) RETURN

            DO J5=1,KN
C
            ZS1 = 1. - PP(J5) / SPSIM
C
            IF (ZS1 .LT. ZEPS) ZS1 = ZEPS
C
            ZPT1 = AP2(1) * (1 - ZS1**AP2(3))**AP2(2)
            ZPT2 = AP2(4) * ZS1**AP2(5)
            IF (AP2(6).GT.0.) ZPT2 = ZPT2 * (1. - EXP((ZS1-1.)/AP2(6)))
            PT(J5) = ZPT1 + ZPT2
C
            ENDDO
C
         ELSE IF (NPPFUN .EQ. 4) THEN
C
***********************************************************************
*                                                                     *
*  INTERPOLATE P' WITH CUBIC SPLINES ON RPPF VALUES (IF P-PRIME IS    *
*  GIVEN BY A SET OF POINTS)                                          *
*                                                                     *
***********************************************************************
C
            CALL PPSPLN(KN,PP,NPPF,FCSM,RPPF,D2RPPF,PT,0)
C
         ELSE IF (NPPFUN .EQ. 5) THEN
C
***********************************************************************
*                                                                     *
*  DEFINE P' WITH A GAUSSIAN CENTERED AT AP(1), OF HEIGHT AP(2) AND   *
*  OF WIDTH AP(3)                                                     *
*                                                                     *
***********************************************************************
C
            ZYS0  = 0.
            ZYS1  = 1.
            ZARG1 = (AP(1) / AP(3))**2
            ZARG2 = ((1. - AP(1)) / AP(3))**2
C
            IF (ZARG1 .LT. 100.) THEN
               ZYS0 = AP(2) * EXP(-ZARG1)
            ENDIF
            IF (ZARG2 .LT. 100.) THEN
               ZYS1 = AP(2) * EXP(-ZARG2)
            ENDIF
C
            ZYSH = MAX(ZYS0,ZYS1)
C
            ZARG1 = (0.1 - AP(1)) / AP(3)
            ZARG2 = (0.9 - AP(1)) / AP(3)
            ZY1   = ZYSH
            ZY2   = ZYSH
C
            IF (ZARG1**2 .LT. 100.) THEN
               ZY1   = - AP(2) * EXP(-ZARG1**2) + ZYSH
            ENDIF
            IF (ZARG2**2 .LT. 100.) THEN
               ZY2   = - AP(2) * EXP(-ZARG2**2) + ZYSH
            ENDIF
C
            ZYP1  = - 2. * ZARG1 * (ZY1 - ZYSH) / AP(3)
            ZYP2  = - 2. * ZARG2 * (ZY2 - ZYSH) / AP(3)
C
            DO 6 J6=1,KN
C
            ZS1 = 1. - PP(J6) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            ZS1 = SQRT(ZS1)
C
            IF (ZS1 .LT. 0.1) THEN
C
               PT(J6) = FCDCD0(0.,0.,0.,0.1,ZY1,ZYP1,ZS1)
C
            ELSE IF (ZS1 .GE. 0.1 .AND. ZS1 .LE. 0.9) THEN
C
               ZEXP = ((ZS1 - AP(1)) / AP(3))**2
C
               IF (ZEXP .LT. 100.) THEN
C            
                  PT(J6) = - AP(2) * EXP(-ZEXP) + ZYSH
C
               ENDIF
C
            ELSE IF (ZS1 .GT. 0.9) THEN
C
               PT(J6) = FCDCD0(0.9,ZY2,ZYP2,1.,0.,0.,ZS1)
C
            ENDIF
C
            IF (PT(J6) .GT. 0.) PT(J6) = 0.
C
   6        CONTINUE
C
         ELSE IF (NPPFUN .EQ. 6) THEN
C
***********************************************************************
*                                                                     *
*  NPP =1 OR 2                                                        *
*                                                                     *
*   PPRIME IS GIVEN AS A SUM OF NPP PROFILES                          *
*                                                                     *
*  FOR S=0     TO S=AP(1)  PPRIME(PSI)=-(ZA0+ZA1*S)                   *
*  FOR S=AP(1) TO S=AP(2)  PPRIME(PSI)=-(ZB0+ZB1*S+ZB2*S**2+ZB3*S**3) *
*  FOR S=AP(2) TO S=AP(3)  PPRIME(PSI)=-(ZC0+ZC1*S+ZC2*S**2)          *
*  FOR S=AP(3) TO S=AP(4)  PPRIME(PSI)=-(ZB0+ZB1*S+ZB2*S**2+ZB3*S**3) *
*  FOR S=AP(4) TO S=AP(5)  PPRIME(PSI)=-(ZA0+ZA1*S)                   *
*  FOR S=AP(5) TO S=1      PPRIME(PSI)=-(ZC0+ZC1*S+ZC2*S**2)          *
*                                                                     *
*     WITH  (AP OR AP2)                                               *
*          F0=AP(6)                                                   *
*                      P1=AP(7)                                       *
*          F2=AP(8)  ; P2=AP(9)                                       *
*          F3=AP(10)                                                  *
*          F4=AP(11) ; P4=AP(12)                                      *
*          F6=AP(13)                                                  *
*                                                                     *
*         FOR FIRST PROFILE  S = S1 = (1-PSI/PSIM)                    *
*             SECOND PROFILE S = S2 = SQRT(S1)                        *
*                                                                     *
***********************************************************************
C
***********************************************************************
*                                                                     *
*  FIRST PROFILE                                                      *
*                                                                     *
***********************************************************************
C
            DO 7 J7=1,KN
C
            ZS(J7) = 1. - PP(J7) / SPSIM
C
            IF (ZS(J7) .LT. 0.) ZS(J7) = 0.
C
   7        CONTINUE
C
            CALL APCOEF2(KN,ZS,AP,PT)
C
            IF (NPP .EQ. 1) RETURN
C
***********************************************************************
*                                                                     *
*  SECOND PROFILE                                                     *
*                                                                     *
***********************************************************************
C
            DO 8 J8=1,KN
C
            ZS1 = 1. - PP(J8) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            ZS(J8) = SQRT(ZS1)
C
   8        CONTINUE
C
            CALL APCOEF2(KN,ZS,AP2,PT)
C
         ELSE IF (NPPFUN .EQ. 7) THEN
C
***********************************************************************
*                                                                     *
*                                                                     *
*  P(PSI) = SPSIM * AP(1) * (1 - S**AP(3))**AP(2)                     *
*  AND P-PRIME(PSI) = D(P(PSI))/D(PSI)                                *
*                                                                     *
*  WHERE  S = 1 - PSI / SPSIM                                         *
*                                                                     *
***********************************************************************
C
            DO 10 J10=1,KN
C
            ZS1 = 1. - PP(J10) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            IF (AP(2).EQ.1..AND.AP(3).EQ.1.) THEN
               PT(J10) = AP(1)
            ELSE IF (AP(2).EQ.1..AND.AP(3).NE.1.) THEN
               PT(J10) = AP(1)*AP(3)**ZS1**(AP(3)-1.)
            ELSE IF (AP(2).NE.1..AND.AP(3).EQ.1.) THEN
               PT(J10) = AP(1)*AP(2)*(1.-ZS1)**(AP(2)-1.)
            ELSE
               PT(J10) = AP(1)*AP(2)*AP(3)*
     &                   (1.-ZS1**AP(3))**(AP(2)-1.)*ZS1**(AP(3)-1.)
            ENDIF
C
  10        CONTINUE
c
            IF (NPP .EQ. 1) RETURN
c
            DO J10=1,KN
C
            ZS1 = 1. - PP(J10) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            IF (AP2(2).EQ.1..AND.AP2(3).EQ.1.) THEN
               PT(J10) = AP2(1)
            ELSE IF (AP2(2).EQ.1..AND.AP2(3).NE.1.) THEN
               PT(J10) = AP2(1)*AP2(3)**ZS1**(AP2(3)-1.)
            ELSE IF (AP2(2).NE.1..AND.AP2(3).EQ.1.) THEN
               PT(J10) = AP2(1)*AP2(2)*(1.-ZS1)**(AP2(2)-1.)
            ELSE
               PT(J10) = AP2(1)*AP2(2)*AP2(3)*
     &                   (1.-ZS1**AP2(3))**(AP2(2)-1.)*ZS1**(AP2(3)-1.)
            ENDIF
C
            ENDDO
C
         ENDIF
C
         RETURN
         END
C*DECK C2SP02
C*CALL PROCESS
         SUBROUTINE BSFUNC(KN,PP,PT)
C        ###########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP02  FRACTION OF BOOTSTRAP CURRENT                              *
*          NBSFUN = 1  -----> POLYNOMIAL                              *
*          NBSFUN = 2  -----> POLYNOMIAL IN 3 SECTIONS                *
*          NBSFUN = 3  -----> PRINCETON DEFINITION                    *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R   PP(KN),   PT(KN),    ZS(2*NPISO+NPT)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         INCLUDE 'CUCCCC.inc'
         INCLUDE 'CUCDCD.inc'
         INCLUDE 'QUAQDQ.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(PT,KN)
C
         IF (NBSFUN .EQ. 1) THEN
C
***********************************************************************
*                                                                     *
* BSFUNC(PSI) DEFINED AS POLYNOMIAL IN PSI/PSIMIN OF DEGREE NSOUR     *
*                                                                     *
***********************************************************************
C
            CALL RESETR(PT,KN,AFBS(NSOUR))
C
            DO 2 J2=NSOUR-1,1,-1
C
            DO 1 J1=1, KN
C
            ZS1    = PP(J1) / SPSIM
            PT(J1) = PT(J1) * ZS1 + AFBS(J2)
C
    1       CONTINUE
    2       CONTINUE
C
         ELSE IF (NBSFUN .EQ. 2) THEN
C
***********************************************************************
*                                                                     *
*  NPP =1 OR 2                                                        *
*                                                                     *
*   BSFUNC IS GIVEN AS A SUM OF NPP PROFILES                          *
*                                                                     *
*  FOR S=0     TO S=AFBS(1)  BSFUNC(PSI)=-(ZA0+ZA1*S)                 *
*  FOR S=AFBS(1) TO S=AFBS(2)  BSFUNC(PSI)=-(ZB0+ZB1*S+ZB2*S**2+      *
*                                           ZB3*S**3)                 *
*  FOR S=AFBS(2) TO S=1      BSFUNC(PSI)=-(ZC0+ZC1*S+ZC2*S**2)        *
*                                                                     *
*         FOR FIRST PROFILE  S = S1 = (1-PSI/PSIM)                    *
*             SECOND PROFILE S = S2 = SQRT(S1)                        *
*                                                                     *
***********************************************************************
C
***********************************************************************
*                                                                     *
*  FIRST PROFILE                                                      *
*                                                                     *
***********************************************************************
C
            DO 3 J3=1,KN
C
            ZS(J3) = 1. - PP(J3) / SPSIM
C
            IF (ZS(J3) .LT. 0.) ZS(J3) = 0.
C
   3        CONTINUE
C
            CALL APCOEF(KN,ZS,AFBS,PT)
C
            IF (NPP .EQ. 1) RETURN
C
***********************************************************************
*                                                                     *
*  SECOND PROFILE                                                     *
*                                                                     *
***********************************************************************
C
            DO 4 J4=1,KN
C
            ZS1 = 1. - PP(J4) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            ZS(J4) = SQRT(ZS1)
C
   4        CONTINUE
C
            CALL APCOEF(KN,ZS,AFBS2,PT)
C
         ELSE IF (NBSFUN .EQ. 3) THEN
C
***********************************************************************
*                                                                     *
*                                                                     *
*   BSFUNC(PSI) = AFBS(1) * (1 - S**AFBS(3))**AFBS(2)                 *
*                                                                     *
*  WHERE  S = 1 - PSI / SPSIM                                         *
*                                                                     *
***********************************************************************
C
            DO 5 J5=1,KN
C
            ZS1 = 1. - PP(J5) / SPSIM
C
            IF (ZS1 .LT. 0.) ZS1 = 0.
C
            zpt1 = AFBS(1) * (1 - ZS1**AFBS(3))**AFBS(2)
            zpt2 = afbs(4) * zs1**afbs(5)
            if (ap(6).gt.0.) zpt2 = zpt2 * (1.-exp((zs1-1.)/afbs(6)))
            pt(j5) = zpt1 + zpt2
C
   5        CONTINUE
C
         ENDIF
C
         RETURN
         END
C*DECK C2SP03
C*CALL PROCESS
         SUBROUTINE PPSPLN(KN,PP,KPP,PS,RPP,D2RPP,PT,KD)
C        ############################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP03  CUBIC SPLINE INTEPOLATION OF EXPERIMENTAL P-PRIME PROFILE  *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   I1(NPT+2*NPISO),  IC(NPT+2*NPISO),
     R   D2RPP(NPISO),   PP(KN),   PS(NPISO),
     R   RPP(NPISO),     PT(KN)
C
C
C   BRACKET OUT [PS(I); PS(I+1)] INTERVAL SUCH THAT
C   PSIISO(I) <= PP(J) <= PSIISO(I+1), J=1,...,KN
C 
         IF (KPP.GT.NPISO) STOP 'KPP>NPISO'

         CALL RESETI(IC,KN,1)
         DO 1 JS = 1,KPP+1
           DO 1 JG=1,KN
             IF (IC(JG).EQ.0) GOTO 1
             ZS1 = 1. - PP(JG) / SPSIM
             IF (ZS1 .LT. 0.) ZS1 = 0.
             I1(JG) = JS-1
             IF (SQRT(ZS1).LE.PS(JS)) IC(JG) = 0
 1       CONTINUE
C
***********************************************************************
*                                                                     *
*  COMPUTE P-PRIME                                                    *
*                                                                     *
***********************************************************************
C
         DO 2 J2=1,KN
C
         IF (I1(J2) .LT. 1)   I1(J2) = 1
         IF (I1(J2) .GT. KPP) I1(J2) = KPP
C
         ZS1 = 1. - PP(J2) / SPSIM
C
         IF (ZS1 .LT. 0.) ZS1 = 0.
C
         ZS1 = SQRT(ZS1)
C
         ZH = PS(I1(J2)+1) - PS(I1(J2))
         ZA = (PS(I1(J2)+1) - ZS1) / ZH
         ZB = (ZS1 - PS(I1(J2))) / ZH
         ZC = (ZA + 1) * (ZA - 1) * ZH * (PS(I1(J2)+1) - ZS1) / 6.
         ZD = (ZB + 1) * (ZB - 1) * ZH * (ZS1 - PS(I1(J2))) / 6.
C 
CYQL2018
         IF (KD.EQ.0) THEN
         PT(J2) = ZA*RPP(I1(J2))   + ZB*RPP(I1(J2)+1) +
     +            ZC*D2RPP(I1(J2)) + ZD*D2RPP(I1(J2)+1)

         ELSE IF (KD.EQ.1) THEN
C        FIRST DERIVATIVE
         PT(J2) =-1./ZH*RPP(I1(J2)) + 1./ZH*RPP(I1(J2)+1) +
     +           (1./3.-ZA*ZA)*ZH/2.*D2RPP(I1(J2)) +
     +           (ZB*ZB-1./3.)*ZH/2.*D2RPP(I1(J2)+1)

         ENDIF

C
    2    CONTINUE
C         
         RETURN
         END
C*DECK C2SP04
C*CALL PROCESS
         SUBROUTINE APCOEF(KN,PX,PAP,PT)
C        ###############################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP04  EVALUATE P' IF GIVEN AS POLYNOMIALS IN SEVERAL SECTIONS    *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PAP(*),  PT(KN),   PX(KN)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL COPYAP(PAP,ZA0,ZA1,ZB0,ZB1,ZB2,ZB3,ZC0,ZC1,ZC2)
C
         DO 1 J1=1,KN
C
         ZSG1 = SIGN(1.,PAP(1) - PX(J1))
         ZSG2 = SIGN(1.,(PX(J1) - PAP(1)) * (PAP(2) - PX(J1)))
         ZSG3 = SIGN(1.,PX(J1) - PAP(2))
C
         ZS1 = MAX(0.,ZSG1)
         ZS2 = MAX(0.,ZSG2)
         ZS3 = MAX(0.,ZSG3)
C
         PT(J1) = PT(J1) - ZS1 * (ZA0 + ZA1 * PX(J1)) -
     -                     ZS2 * (ZB0 + PX(J1) * (ZB1 + PX(J1) * 
     *                                     (ZB2 + PX(J1) * ZB3))) -
     -                     ZS3 * (ZC0 + PX(J1) * (ZC1 + PX(J1) * ZC2))
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SP05
C*CALL PROCESS
         SUBROUTINE COPYAP(PAP,PA0,PA1,PB0,PB1,PB2,PB3,PC0,PC1,PC2)
C        #########################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP05  EVALUATE POLYNOMIAL COEFFICIENTS IF P' IS                  *
*          PRESCRIBED AS POLYNOMIALS IN SEVERALS SECTIONS             *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION  
     R   PAP(*)
C
         INCLUDE 'QUAQDQ.inc'
         INCLUDE 'CUCDCD.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
***********************************************************************
*                                                                     *
*  FIRST SECTION                                                      *
*                                                                     *
***********************************************************************
C
         IF (PAP(1) .NE. 0.) THEN
C
           PA0 = PAP(5)
           PA1 = PAP(4)
C
         ELSE IF (PAP(1) .EQ. 0.) THEN
C
           PA0 = 0.
           PA1 = 0.
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  SECOND SECTION                                                     *
*                                                                     *
***********************************************************************
C
         IF (PAP(1) .NE. PAP(2)) THEN
C
           ZF1 = PA0 + PA1 * PAP(1)
           ZF2 = PAP(3)
           PB3 = FC3(PAP(1),ZF1,PAP(4),PAP(2),ZF2,PAP(7))
           PB2 = FC2(PAP(1),ZF1,PAP(4),PAP(2),ZF2,PAP(7))
           PB1 = FC1(PAP(1),ZF1,PAP(4),PAP(2),ZF2,PAP(7))
           PB0 = FC0(PAP(1),ZF1,PAP(4),PAP(2),ZF2,PAP(7))
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  THIRD SECTION                                                      *
*                                                                     *
***********************************************************************
C
         IF (PAP(2) .NE. 1.) THEN
C
           PC2 = FD2(PAP(2),PAP(3),PAP(7),1.,PAP(6))
           PC1 = FD1(PAP(2),PAP(3),PAP(7),1.,PAP(6))
           PC0 = FD0(PAP(2),PAP(3),PAP(7),1.,PAP(6))
C
         ENDIF
C
         RETURN
         END
C*DECK C2SP06
C*CALL PROCESS
         SUBROUTINE APCOEF2(KN,PX,PAP,PT)
C        ################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP06  EVALUATE P' IF GIVEN AS POLYNOMIALS IN SEVERAL SECTIONS    *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PAP(*),  PT(KN),   PX(KN)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL COPYAPP(PAP,ZA0,ZA1,ZB0,ZB1,ZB2,ZB3,ZC0,ZC1,ZC2,
     &                    ZD0,ZD1,ZD2,ZD3,ZE0,ZE1,ZF0,ZF1,ZF2)
C
         DO 1 J1=1,KN
C
         ZSG1 = SIGN(1.,PAP(1) - PX(J1))
         ZSG2 = SIGN(1.,(PX(J1) - PAP(1)) * (PAP(2) - PX(J1)))
         ZSG3 = SIGN(1.,(PX(J1) - PAP(2)) * (PAP(3) - PX(J1)))
         ZSG4 = SIGN(1.,(PX(J1) - PAP(3)) * (PAP(4) - PX(J1)))
         ZSG5 = SIGN(1.,(PX(J1) - PAP(4)) * (PAP(5) - PX(J1)))
         ZSG6 = SIGN(1.,PX(J1) - PAP(5))
C
         ZS1 = MAX(0.,ZSG1)
         ZS2 = MAX(0.,ZSG2)
         ZS3 = MAX(0.,ZSG3)
         ZS4 = MAX(0.,ZSG4)
         ZS5 = MAX(0.,ZSG5)
         ZS6 = MAX(0.,ZSG6)
C
         PT(J1) = PT(J1) - ZS1 * (ZA0 + ZA1 * PX(J1)) -
     &                     ZS2 * (ZB0 + PX(J1) * (ZB1 + PX(J1) * 
     &                                     (ZB2 + PX(J1) * ZB3))) -
     &                     ZS3 * (ZC0 + PX(J1) * (ZC1 + PX(J1) * ZC2)) -
     &                     ZS4 * (ZD0 + PX(J1) * (ZD1 + PX(J1) * 
     &                                     (ZD2 + PX(J1) * ZD3))) -
     &                     ZS5 * (ZE0 + ZE1 * PX(J1)) -
     &                     ZS6 * (ZF0 + PX(J1) * (ZF1 + PX(J1) * ZF2))
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SP07
C*CALL PROCESS
         SUBROUTINE COPYAPP(PAP,PA0,PA1,PB0,PB1,PB2,PB3,PC0,PC1,PC2,
     &                          PD0,PD1,PD2,PD3,PE0,PE1,PF0,PF1,PF2)
C        ##########################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP07  EVALUATE POLYNOMIAL COEFFICIENTS IF P' IS                  *
*          PRESCRIBED AS POLYNOMIALS IN SEVERALS SECTIONS             *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION  
     R   PAP(*)
C
         INCLUDE 'QUAQDQ.inc'
         INCLUDE 'CUCDCD.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
***********************************************************************
*                                                                     *
*  FIRST SECTION                                                      *
*                                                                     *
***********************************************************************
C
         IF (PAP(1) .NE. 0.) THEN
C
           PA0 = PAP(6)
           PA1 = PAP(7)
C
         ELSE IF (PAP(1) .EQ. 0.) THEN
C
           PA0 = 0.
           PA1 = 0.
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  SECOND SECTION                                                     *
*                                                                     *
***********************************************************************
C
         IF (PAP(1) .NE. PAP(2)) THEN
C
           ZF1 = PA0 + PA1 * PAP(1)
           ZF2 = PAP(8)
           PB3 = FC3(PAP(1),ZF1,PAP(7),PAP(2),ZF2,PAP(9))
           PB2 = FC2(PAP(1),ZF1,PAP(7),PAP(2),ZF2,PAP(9))
           PB1 = FC1(PAP(1),ZF1,PAP(7),PAP(2),ZF2,PAP(9))
           PB0 = FC0(PAP(1),ZF1,PAP(7),PAP(2),ZF2,PAP(9))
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  THIRD SECTION                                                      *
*                                                                     *
***********************************************************************
C
         IF (PAP(2) .NE. PAP(3)) THEN
C
           PC2 = FD2(PAP(2),PAP(8),PAP(9),PAP(3),PAP(10))
           PC1 = FD1(PAP(2),PAP(8),PAP(9),PAP(3),PAP(10))
           PC0 = FD0(PAP(2),PAP(8),PAP(9),PAP(3),PAP(10))
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  FOURTH SECTION                                                     *
*                                                                     *
***********************************************************************
C
         IF (PAP(3) .NE. PAP(4)) THEN
C
           ZP3 = FQDQ1(PAP(2),PAP(8),PAP(9),PAP(3),PAP(10),PAP(3))
           PD3 = FC3(PAP(3),PAP(10),ZP3,PAP(4),PAP(11),PAP(12))
           PD2 = FC2(PAP(3),PAP(10),ZP3,PAP(4),PAP(11),PAP(12))
           PD1 = FC1(PAP(3),PAP(10),ZP3,PAP(4),PAP(11),PAP(12))
           PD0 = FC0(PAP(3),PAP(10),ZP3,PAP(4),PAP(11),PAP(12))
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  FIFTH SECTION                                                      *
*                                                                     *
***********************************************************************
C
         IF (PAP(4) .NE. PAP(5)) THEN
C
           PE0 = PAP(11) - PAP(12) * PAP(4)
           PE1 = PAP(12)
C
         ENDIF
C
***********************************************************************
*                                                                     *
*  SIXTH SECTION                                                      *
*                                                                     *
***********************************************************************
C
         IF (PAP(5) .NE. 1.) THEN
C
           ZF5 = PAP(11) + (PAP(5) - PAP(4)) * PAP(12)
           PF2 = FD2(PAP(5),ZF5,PAP(12),1.,PAP(13))
           PF1 = FD1(PAP(5),ZF5,PAP(12),1.,PAP(13))
           PF0 = FD0(PAP(5),ZF5,PAP(12),1.,PAP(13))
C
         ENDIF
C
         RETURN
         END
C*DECK C2SP08
C*CALL PROCESS
         SUBROUTINE BLTEST
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP08  LEAD EVALUATION OF SURFACE QUANTITIES, BALLOONING STABILITY*
*          AND LOCAL INTERCHANGE CRITERIA AND GLOBAL EQUILIBRIUM      *
*          QUANTITIES DURING BALLOONING OPTIMIZATION                  *
*          LIMIT P' ACCORDING TO EQ. (41) IN PUBLICATION              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IP = ISRCHFGE(NPPR+1,PSIISO,1,CPSICL(1))
C
         IF (IP.LT.1)      IP = 1
         IF (IP.GT.NPPR+1) IP = NPPR+1
C
         DO 1 J1=IP,NPPR+1
C
         CALL SURFACE(J1,SIGPSI(1,J1),TETPSI(1,J1),WGTPSI(1,J1),
     ,                PCSM(J1))
         CALL CHIPSI(NPPR+1,J1)
C
 1       CONTINUE
C
c%OS         IF (NPPR+1-IP+1 .LE. NPPSBAL) THEN
c%OS           CALL BALOON(IP,NPPR+1,PCSM)
c%OS         ELSE
           DO J5=IP,NPPR+1,NPPSBAL
             JEND = MIN(J5+NPPSBAL-1,NPPR+1)
             CALL BALOON(J5,JEND,PCSM)
           ENDDO
c%OS         ENDIF
C
         IF (IP .GT. 1) THEN
C
            DO 2 J2=1,IP-1
C
            CPR(J2)   = CPR(IP)
            CPPR(J2)  = 0.
            NCBAL(J2) = - 1
            SMERCI(J2) = 0.25
C
   2        CONTINUE
C
         ENDIF
C
         CP0   = FCCCC0(CPR(IP),CPR(IP+1),CPR(IP+2),CPR(IP+3),
     ,                  PCSM(IP),PCSM(IP+1),PCSM(IP+2),PCSM(IP+3),0.)
         DPDP0 = FCCCC0(CPPR(IP),CPPR(IP+1),CPPR(IP+2),CPPR(IP+3),
     ,                  PCSM(IP),PCSM(IP+1),PCSM(IP+2),PCSM(IP+3),0.)
C
         CALL GLOQUA(PCSM,PCS,NPPR+1,2)
         CALL RESETI(N2BAL,NPPR+1,0)
cab
      IF (.TRUE.) RETURN
cab
C
cab      DO 3 J3=IP,NPPR+1
C
cab      IF (ABS(CPPR(J3)) .GT. ABS(CFBAL * CDQ(J3)) .OR. 
cab   ,       ABS(CPPR(J3)) .GT. 5.) THEN
C
cab         NCBAL(J3) = 1 
cab         N2BAL(J3) = 1
C
cab       ENDIF
C
    3    CONTINUE
C
         DO 4 J4=NPPR-2,NPPR+1
C
         IF (ABS(CPPR(J4)) .GT. ABS(CPPR(J4-1))) THEN
C
            CPPR(J4) = CPPR(J4-1)
C
            PRINT*,'*****WARNING***** : CPPR(',J4,') = CPPR(',J4-1,')'
            PRINT*,' IMPOSED. NPPR = ',NPPR
C
            NCBAL(J4) = 1 
            N2BAL(J4) = 1
C
         ENDIF
C
    4    CONTINUE
C
         RETURN
         END
C*DECK C2SP09
C*CALL PROCESS
         SUBROUTINE RESPPR
C        #################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP09  INITIALIZE ARRAYS FOR BALLOONING OPTIMIZATION              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE "COMDIM.inc"
         INCLUDE "COMBAL.inc"
         INCLUDE "COMNUM.inc"
         INCLUDE "COMSUR.inc"
C
         NSRCH  = 3
         XPRCNT = 1.
C
         CALL RESETR(XLAMB,NPPR+1,1.)
         CALL VZERO(XPPRMN,NPPR+1)
         CALL VZERO(XPPRMX,NPPR+1)
         CALL VZERO(XPPRDF,NPPR+1)
         CALL VZERO(XP0,NPPR+1)
         CALL VZERO(XP1,NPPR+1)
         CALL VZERO(XP2,NPPR+1)
         CALL VZERO(XP3,NPPR+1)
         CALL VZERO(XP4,NPPR+1)
C
         CALL RESETI(NP0,NPPR+1,-10)
         CALL RESETI(NP1,NPPR+1,0)
         CALL RESETI(NP2,NPPR+1,0)
         CALL RESETI(NP3,NPPR+1,0)
         CALL RESETI(NP4,NPPR+1,0)
C
         RETURN
         END
C*DECK C2SP11
C*CALL PROCESS
         SUBROUTINE PPBSTR(KFIN)
C        #######################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SP10  MODIFY P' WHEN BOOTSTRAP CURRENT DENSITY IS PRESCRIBED     *
*          ACCORDING TO EQ. (42) IN THE PUBLICATION.                  *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION ZBSFR(NPISO),ZDIFF(NPISO),ZWORK(NPISO),ZWORK1(NPISO)
C
         PARAMETER (ZMCOR = 0.50, ZEPS = 1.E-3)
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IP = ISRCHFGE(NPPR+1,PSIISO,1,CPSICL(1))
C
         IF (IP.LT.1)      IP = 1
         IF (IP.GT.NPPR+1) IP = NPPR+1
C
         DO 1 J1=IP,NPPR+1
C
         CALL SURFACE(J1,SIGPSI(1,J1),TETPSI(1,J1),WGTPSI(1,J1),
     ,                PCSM(J1))
         CALL CHIPSI(NPPR+1,J1)
C
   1     CONTINUE
C
         IF (IP .GT. 1) THEN
C
            DO 2 J2=1,IP-1
C
            CPR(J2)   = CPR(IP)
            CPPR(J2)  = 0.
C
   2        CONTINUE
C
         ENDIF
C
         CP0   = FCCCC0(CPR(IP),CPR(IP+1),CPR(IP+2),CPR(IP+3),
     ,                  PCSM(IP),PCSM(IP+1),PCSM(IP+2),PCSM(IP+3),0.)
         DPDP0 = FCCCC0(CPPR(IP),CPPR(IP+1),CPPR(IP+2),CPPR(IP+3),
     ,                  PCSM(IP),PCSM(IP+1),PCSM(IP+2),PCSM(IP+3),0.)
C
         CALL GLOQUA(PCSM,PCS,NPPR+1,2)
         CALL SCOPY(NPPR+1,RPRM,1,XP0,1)
C
         IF (NBSTRP .EQ. 2) THEN
C
            CALL BSFUNC(NPPR+1,PSIISO,ZBSFR)
C
         ENDIF
C
         ZERR    = 0.
C
         DO 3 J3=1,NPPR+1
C
         ZRPRM  = RPRM(J3)
C
         IF (NBSTRP .EQ. 1) THEN
C
c%OS            ZNRPRM = RPRM(J3) * RJPAR(J3) * BSFRAC / (RJBSH(J3)+1.E-6)
            ZNRPRM = RPRM(J3) * RJPAR(J3)*BSFRAC / (RJBSOS(J3,1)+1.E-6)
C
         ELSE IF (NBSTRP .EQ. 2) THEN
C
c%OS            ZNRPRM = RPRM(J3) * RJPAR(J3) * ZBSFR(J3) / 
c%OS     &               (RJBSH(J3)+1.E-6)
            ZNRPRM = RPRM(J3) * RJPAR(J3) * ZBSFR(J3) / 
     &               (RJBSOS(J3,1)+1.E-6)
C
         ENDIF
C
         IF (ZNRPRM .LT. ZRPRM - ZMCOR) ZNRPRM = ZRPRM - ZMCOR
         IF (ZNRPRM .GT. ZRPRM + ZMCOR) ZNRPRM = ZRPRM + ZMCOR
         IF (ZNRPRM .LT. - CFBAL)       ZNRPRM = - CFBAL
C
         RPRM(J3)  = ZNRPRM
C
         IF (RPRM(J3) .GT. 0.) RPRM(J3) = 0.
C
         ZDIFF(J3) = ABS(RPRM(J3) - XP0(J3))
C
         IF (ZDIFF(J3) .GT. ZERR) ZERR = ZDIFF(J3)
C
   3     CONTINUE
C
         WRITE(*,*) 
C
         IF (ZERR .LE. ZEPS) THEN
C
            WRITE(6,110) 'EPSILON > ERROR'
            WRITE(6,*) ZEPS,ZERR
C
            KFIN = 3
C
         ELSE
C
            WRITE(*,*) 'ERROR > EPSILON ',ZERR,ZEPS
C
         ENDIF

         CALL SPLINE(PCSM,RPRM,NPPR+1,D2RPRM,ZWORK,ZWORK1)
C
  110    FORMAT (A,$)
C
         RETURN
         END
C*DECK C2ST01
C*CALL PROCESS
         SUBROUTINE POLYNM(KN,PCOEF,PT,K)
C        ################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2ST01  COMPUTE PLASMA DENSITY OR TEMPERATURE GIVEN AS A POLYNOMIAL*
*          OF DEGREE NSOUR < 11 IN CARTESIAN COORDINATES              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
C
         DIMENSION
     R   PCOEF(*),            PT(KN),             ZP(2*NPISO+NPT),
     R   ZS(2*NPISO+NPT),     ZSIG(2*NPISO+NPT),  ZTET(2*NPISO+NPT),
     R   ZDRODP(2*NPISO+NPT)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         INCLUDE 'CUCCCC.inc'
         INCLUDE 'CUCDCD.inc'
         INCLUDE 'QUAQDQ.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         CALL VZERO(PT,KN)
         CALL RMRAD(KN,SPSIM,0.,PANGLE,ZS,ZSIG,ZTET,1)
C
         IF (K .EQ. 1) THEN
C
            CALL RESETR(ZP,KN,1.)
C
            DO 2 J2=1,NSOUR
C
            DO 1 J1=1,KN
C
            PT(J1) = PT(J1) + PCOEF(J2) * ZP(J1)
            ZP(J1) = ZS(J1) * ZP(J1)
C
    1       CONTINUE
    2       CONTINUE
C
            DO 3 J3=1,KN
C
            IF (PT(J3) .LT. 0.) PT(J3) = 0.
C
    3       CONTINUE
C
         ELSE
C
            CALL RESETR(ZP,KN,1.)
C
            DO 5 J5=2,NSOUR
C
            DO 4 J4=1,KN
C
            PT(J4) = PT(J4) + (J5 - 1) * PCOEF(J5) * ZP(J4)
            ZP(J4) = ZS(J4) * ZP(J4)
C
    4       CONTINUE
    5       CONTINUE
C
            CALL DRHODP(KN,PANGLE,ZS,ZDRODP)
C
            DO 6 J6=1,KN
C
            PT(J6) = PT(J6) * ZDRODP(J6)
C
    6       CONTINUE
C
         ENDIF
C
         RETURN
         END
C*DECK C2ST02
C*CALL PROCESS
         SUBROUTINE DRHODP(KN,PT,PRHO,PDRODP)
C        ####################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2ST02  EVALUATE D(RHO)/D(PSI) FOR A GIVEN SET OF PSI VALUES,      *
*          WHERE RHO IS A NORMALIZED RADIUS ALONG A LINE STARTING     *
*          AT THE MAGNETIC AXIS AND WITH AN ANGLE OF PT WITH          *
*          RESPECT TO THE Z=0 PLANE                                   *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IS0(2*NPISO+NPT),       IT0(2*NPISO+NPT),   
     R   PRHO(KN),               PDRODP(KN),             
     R   ZBND(2*NPISO+NPT),
     R   ZBND1(2*NPISO+NPT),     ZBND2(2*NPISO+NPT), 
     R   ZDFDS(2*NPISO+NPT,16),  ZDFDT(2*NPISO+NPT,16),
     R   ZPCEL(2*NPISO+NPT,16),  ZRHO(2*NPISO+NPT),    
     R   ZS1(2*NPISO+NPT),       ZS2(2*NPISO+NPT),       
     R   ZTETA(2*NPISO+NPT),     ZTETA1(2*NPISO+NPT),
     R   ZTETA2(2*NPISO+NPT),    ZT1(2*NPISO+NPT),       
     R   ZT2(2*NPISO+NPT) ,      ZSIGMA(2*NPISO+NPT)    
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         BPS( 1) = RMAG
         BPS(12) = RZMAG
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
         CALL BOUND(1,PT,ZBNDT0)
C
         BPS( 1) = R0
         BPS(12) = RZ0
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
         ZTEPS = 1.E-4
C
         DO 1 J1=1,KN
C
         ZR         = PRHO(J1) * ZBNDT0 * COS(PT) + RMAG
         ZZ         = PRHO(J1) * ZBNDT0 * SIN(PT) + RZMAG
         ZRHO(J1)   = SQRT((ZR - R0)**2 + (ZZ - RZ0)**2)
         ZTETA(J1)  = ATAN2(ZZ - RZ0,ZR - R0)
C
         IF (ZTETA(J1) .LT. CT(1)) ZTETA(J1) = ZTETA(J1) + 2. * CPI
C
         ZTETA1(J1) = ZTETA(J1) - ZTEPS
         ZTETA2(J1) = ZTETA(J1) + ZTEPS
C
    1    CONTINUE
C
         CALL BOUND(KN,ZTETA,ZBND)
         CALL BOUND(KN,ZTETA1,ZBND1)
         CALL BOUND(KN,ZTETA2,ZBND2)
C
         DO 2 J2=1,KN
C
         ZSIGMA(J2) = ZRHO(J2) / ZBND(J2)
C
         IS0(J2) = ISRCHFGE(NS1,CSIG,1,ZSIGMA(J2)) - 1
C
         IF (IS0(J2) .GT. NS) THEN
C
            IS0(J2) = NS
C
         ELSE IF (IS0(J2) .LT. 1) THEN
C
            IS0(J2) = 1
C
         ENDIF
C
         IT0(J2) = ISRCHFGE(NT1,CT,1,ZTETA(J2)) - 1
C
         IF (IT0(J2) .GT. NT) THEN
C
            IT0(J2) = NT
C
         ELSE IF (IT0(J2) .LT. 1) THEN
C
            IT0(J2) = 1
C
         ENDIF
C
         ZS1(J2) = CSIG(IS0(J2))
         ZS2(J2) = CSIG(IS0(J2)+1)
         ZT1(J2) = CT(IT0(J2))
         ZT2(J2) = CT(IT0(J2)+1)
C
    2    CONTINUE
C
         CALL BASIS2(KN,2*NPISO+NPT,ZS1,ZS2,ZT1,ZT2,ZSIGMA,ZTETA,ZDFDS,
     &               ZDFDT)
         CALL PSICEL(IS0,IT0,KN,2*NPISO+NPT,ZPCEL,CPSICL)
C
         DO 3 J3=1,KN
C
         ZDPDS = ZPCEL(J3, 1) * ZDFDS(J3, 1) +
     +           ZPCEL(J3, 2) * ZDFDS(J3, 2) +
     +           ZPCEL(J3, 3) * ZDFDS(J3, 3) +
     +           ZPCEL(J3, 4) * ZDFDS(J3, 4) +
     +           ZPCEL(J3, 5) * ZDFDS(J3, 5) +
     +           ZPCEL(J3, 6) * ZDFDS(J3, 6) +
     +           ZPCEL(J3, 7) * ZDFDS(J3, 7) +
     +           ZPCEL(J3, 8) * ZDFDS(J3, 8) +
     +           ZPCEL(J3, 9) * ZDFDS(J3, 9) +
     +           ZPCEL(J3,10) * ZDFDS(J3,10) +
     +           ZPCEL(J3,11) * ZDFDS(J3,11) +
     +           ZPCEL(J3,12) * ZDFDS(J3,12) +
     +           ZPCEL(J3,13) * ZDFDS(J3,13) +
     +           ZPCEL(J3,14) * ZDFDS(J3,14) +
     +           ZPCEL(J3,15) * ZDFDS(J3,15) +
     +           ZPCEL(J3,16) * ZDFDS(J3,16) 
C 
        ZDPDT =  ZPCEL(J3, 1) * ZDFDT(J3, 1) +
     +           ZPCEL(J3, 2) * ZDFDT(J3, 2) +
     +           ZPCEL(J3, 3) * ZDFDT(J3, 3) +
     +           ZPCEL(J3, 4) * ZDFDT(J3, 4) +
     +           ZPCEL(J3, 5) * ZDFDT(J3, 5) +
     +           ZPCEL(J3, 6) * ZDFDT(J3, 6) +
     +           ZPCEL(J3, 7) * ZDFDT(J3, 7) +
     +           ZPCEL(J3, 8) * ZDFDT(J3, 8) +
     +           ZPCEL(J3, 9) * ZDFDT(J3, 9) +
     +           ZPCEL(J3,10) * ZDFDT(J3,10) +
     +           ZPCEL(J3,11) * ZDFDT(J3,11) +
     +           ZPCEL(J3,12) * ZDFDT(J3,12) +
     +           ZPCEL(J3,13) * ZDFDT(J3,13) +
     +           ZPCEL(J3,14) * ZDFDT(J3,14) +
     +           ZPCEL(J3,15) * ZDFDT(J3,15) +
     +           ZPCEL(J3,16) * ZDFDT(J3,16) 
C
         ZDRSDT = .5 * (ZBND2(J3) - ZBND1(J3)) / ZTEPS
         ZSIN   = SIN(ZTETA(J3))
         ZCOS   = COS(ZTETA(J3))
C
         ZDSDR = (ZDRSDT * ZSIN + ZBND(J3) * ZCOS) / ZBND(J3)**2
         ZDTDR = - ZSIN / ZRHO(J3) 
         ZDSDZ = (ZBND(J3) * ZSIN - ZDRSDT * ZCOS) / ZBND(J3)**2
         ZDTDZ = ZCOS / ZRHO(J3)  
C
         ZDPDR = ZDPDS * ZDSDR + ZDPDT * ZDTDR
         ZDPDZ = ZDPDS * ZDSDZ + ZDPDT * ZDTDZ
C
         PDRODP(J3) = 1. / (ZBNDT0 * 
     *                 (COS(PT) * ZDPDR + SIN(PT) * ZDPDZ))
C
    3    CONTINUE
C
         RETURN
         END
C*DECK C2ST03
C*CALL PROCESS
         SUBROUTINE POLYFUN(KSOUR,PCOEF,KN,PPSI,PY,KOPT,KOPTX)
C        ###################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2ST03  COMPUTE FUNCTION PY OR ITS NTH DERIVATIVE AS POLYNOMIAL    *
*          OF DEGREE KSOUR IN S,S**2,ETC. DEPENDING ON KOPTX     *
*          KOPT = 1: COMPUTE PY                                       *
*          KOPT = 2: COMPUTE D(PY)/DPSI                               *
*          KOPT = 3: TO BE ADDED WHEN NEEDED                          *
*                                                                     *
*          KOPTX = 1: FUNCTION OF S**2                                *
*          KOPTX = 2: FUNCTION OF PSI/PSIMIN=1.-S**2                  *
*          KOPTX = 3: TO BE ADDED WHEN NEEDED                         *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMPHY.inc'
C
         DIMENSION
     R   PCOEF(KSOUR),            PY(KN), PPSI(KN)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IF (KOPT .EQ. 1) THEN
           J2END = 1
           CALL RESETR(PY,KN,PCOEF(KSOUR))
         ELSE IF (KOPT .EQ. 2) THEN
           J2END = 2
           CALL RESETR(PY,KN,DFLOAT(KSOUR-1)*PCOEF(KSOUR))
         ENDIF
         DO 2 J2 = KSOUR-1,J2END,-1
C
           DO 1 J1 = 1, KN
C
             IF (KOPTX .EQ. 1) THEN
               ZS1 = 1. - PPSI(J1) / SPSIM
               ZDYDPSI = - 1./SPSIM
             ELSE IF (KOPTX .EQ. 2) THEN
               ZS1 = PPSI(J1) / SPSIM
               ZDYDPSI =  1./SPSIM
             ENDIF
C             
             IF (KOPT .EQ. 1) THEN
               PY(J1) = PY(J1) * ZS1 + PCOEF(J2)
             ELSE IF (KOPT .EQ. 2) THEN
               IF (J2 .EQ. KSOUR-1) PY(J1) = PY(J1) * ZDYDPSI
               PY(J1) = PY(J1) * ZS1 + DFLOAT(J2-1)*PCOEF(J2)*ZDYDPSI
             ENDIF
C     
 1         CONTINUE
 2       CONTINUE
C
C
         RETURN
         END
C*DECK C2SU01
C*CALL PROCESS
         SUBROUTINE ISOFIND(K1,K2,PSIGMA,PTETA,PGWGT,PSIAXE,PSIBND)
C        ###########################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SU01  LEAD TRACING OF CONSTANT POLOIDAL FLUX SURFACES            *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     I   IC(2*NPISO),                  IT0(2*NPISO),
     R   PSIGMA(NPMGS*NTP1,2*NPISO),   PTETA(NPMGS*NTP1,2*NPISO), 
     R   PGWGT(NPMGS*NTP1,2*NPISO),    ZSIG(2*NTP1*(NPMGS+1)*NPPSCUB),
     R   ZTET(2*NTP1*(NPMGS+1)*NPPSCUB), ZPISO(2*NPISO),
     R   ZRAC(NPMGS+1),                ZWGT(NPMGS+1)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C
         WRITE(*,*) 'ISOFIND(K1,K2,PSIGMA,PTETA,PGWGT,PSIAXE,PSIBND)'
         WRITE(*,*) K1,K2

         KP = K2 - K1 + 1
C
         CALL GAUSS(NMGAUS,ZRAC,ZWGT)
C     
         CALL RESETI(IC(K1),KP,1)
         DO 1 JT=1,NT1
            DO 1 JP=K1,K2
               IF (IC(JP) .EQ.0) GOTO 1
               IT0(JP) = JT
               IF (TETMAP(1,JP).LT.CT(JT)) IC(JP)  = 0
 1       CONTINUE
C
C     TETMAP STARTS FROM TETMAP(1) AND THEN FOLLOWS THE VALUES
C     OF CT AROUND. IF TETMAP(1) IN [CT(I-1),CT(I)[ THEN (IT0=I):
C     TETMAP(JT=2)     = CT(I)
C     TETMAP   (3)     = CT(I+1)
C     :                :  :
C     TETMAP(NT+1-I+2) = CT(NT+1) = CT(1) - 2.*CPI
C     TETMAP(NT+4-I)   = CT(NT+4-I+I-2 - NT) = CT(2)
C     :                :  :
C     TETMAP(NT+1) = TETMAP(NT+4-I+(I-3)) = CT(2+I-3) = CT(I-1)
C
C     THUS ONE HAS THE VALUES OF CT PLUS ONE EXTRA VALUE AT TETMAP(1)
C     => NT+1 POINTS
C
C     NOTE THAT THE INDEX OF THE JUMP OF 2*PI IS CHANGED BELOW
C
         DO 2 JT=2,NT1
           DO JP=K1,K2
             ZC = 0.
             IT = IT0(JP)+JT-2
             IF (IT.GT.NT1) THEN
               IT = IT-NT
               ZC = 1.
             ENDIF
             TETMAP(JT,JP) = CT(IT) + 2. * ZC * CPI
           ENDDO
 2       CONTINUE
C
C     IF TETMAP(1)=CT(I-1) THEN TETMAP(NT+1)=TETMAP(1), WHICH IS BAD,
C     THUS WE MODIFY TETMAP(NT+1) TO BE IN BETWEEN CT(I-2) AND CT(I-1),
C     I.E IN BETWEEN TETMAP(NT) AND TETMAP(NT+1)
C
         DO JP=K1,K2
           IF (TETMAP(1,JP) .EQ. CT(IT0(JP)-1))
     +       TETMAP(NT+1,JP) = 0.5*(TETMAP(NT,JP) + TETMAP(NT+1,JP))
         ENDDO
C
C     SET TETMAP(NT+2) = TETMAP(1) + 2.*CPI
C     => TETPSI=PTETA IS DEFINED BETWEEN TETMAP(1) AND TETMAP(NT+2)
C     BUT WITH THE JUMP OF 2.*PI ALWAYS AT CT(NT+1)=CT(1)+2*PI
C     THUS AT TETMAP(NT+1-I+2) (NORMALLY, I=1 TO 3)
C
         DO 3 JP=K1,K2
           ZPISO(JP-K1+1) = PSIISO(JP)
 3         TETMAP(NT2,JP) = TETMAP(1,JP) + 2. * CPI  
C           
         DO 4 JG=1,NMGAUS
           DO 4 JT=1,NT1
             DO 4 JP=K1,K2
               ZTDIF = TETMAP(JT+1,JP) - TETMAP(JT,JP)
               ZTADD = TETMAP(JT+1,JP) + TETMAP(JT,JP)
               IG = NMGAUS*(JT-1)+JG
               PTETA(IG,JP) = .5*(ZTADD + ZTDIF*ZRAC(JG))
               IF (PTETA(IG,JP) .GE. CT(NT1)) PTETA(IG,JP) =
     +                             PTETA(IG,JP) - 2. * CPI
               IF (PTETA(IG,JP) .LT. CT(1)) PTETA(IG,JP) =
     +                             PTETA(IG,JP) + 2. * CPI
 4             PGWGT(IG,JP) = ZWGT(JG) * ZTDIF
C
C-----------------------------------------------------------------------
C     COMPUTE ROOTS FOR NPPSCUB ISO-PSI SURFACES AT A TIME (FOR MEMORY)
C
         IKTET = NT1 * (NMGAUS+1)
         DO 100 JPSTART=K1,K2,NPPSCUB
           JK1 = JPSTART
           JK2 = MIN(JPSTART+NPPSCUB-1,K2)
           JPLEN = JK2 - JK1 + 1
C
         DO 7 JT=1,NT1
           DO 5 JP=JK1,JK2
C
C     TETMAP: RESET THE JUMP OF 2.*CPI TO BE AT CT(NT1) => JT+1=NT+1-I+2
C
              IF (TETMAP(JT+1,JP) .GE. CT(NT1)) 
     +                TETMAP(JT+1,JP) = TETMAP(JT+1,JP) - 2.*CPI
 5            ZTET(((JP-JK1)*NT1+JT)*(NMGAUS+1)) = TETMAP(JT+1,JP)
           DO 6 JG=1,NMGAUS
             DO 6 JP=JK1,JK2
 6              ZTET(((JP-JK1)*NT1+JT-1)*(NMGAUS+1)+JG) = 
     +                                  PTETA(NMGAUS*(JT-1)+JG,JP)
 7       CONTINUE
C     
         NROOT = JPLEN * IKTET
CC         ISTART = (JK1-JK1)*IKTET + 1
         CALL CUBRT(NROOT,PSIAXE,PSIBND,ZPISO(JK1-K1+1),ZTET,ZSIG)
C
         DO 10 JT=1,NT1
           DO 8 JP=JK1,JK2
 8           SIGMAP(JT+1,JP) = ZSIG(((JP-JK1)*NT1+JT)*(NMGAUS+1))
           DO 9 JG=1,NMGAUS
             DO 9 JP=JK1,JK2
             PSIGMA(NMGAUS*(JT-1)+JG,JP) = 
     +               ZSIG(((JP-JK1)*NT1+JT-1)*(NMGAUS+1)+JG)
 9        CONTINUE
 10       CONTINUE
C
 100      CONTINUE
C
          DO 11 JP=K1,K2
 11         SIGMAP(NT2,JP) = SIGMAP(1,JP)
C     
           END
C*DECK C2SU02
C*CALL PROCESS
         SUBROUTINE CUBRT(KN,PSIAXE,PSIBND,PISO,PTET,PSIG)
C        #################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SU02  COMPUTE INTERSECTIONS OF CONSTANT FLUX SURFACE WITH        *
*          THETA=PTET(I) LINE.                                        *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IS0(2*NTP1*(NPMGS+1)*NPPSCUB),  IT0(2*NTP1*(NPMGS+1)*NPPSCUB),
     I   ID(2*NTP1*(NPMGS+1)*NPPSCUB),   
     R   PSIG(KN),                     PTET(KN),                
     R   ZBND(2*NTP1*(NPMGS+1)*NPPSCUB), IC(2*NTP1*(NPMGS+1)*NPPSCUB),
     R   ZF(2*NTP1*(NPMGS+1)*NPPSCUB,16),ZFMIN(2*NTP1*(NPMGS+1)*NPPSCUB)
     R  ,ZFMAX(2*NTP1*(NPMGS+1)*NPPSCUB),
     R   ZPCEL(2*NTP1*(NPMGS+1)*NPPSCUB,16),
     R   ZPSI(2*NTP1*(NPMGS+1)*NPPSCUB), PISO(2*NPISO),
     R   ZSIGMN(2*NTP1*(NPMGS+1)*NPPSCUB),
     R   ZSIGMX(2*NTP1*(NPMGS+1)*NPPSCUB), 
     R   ZS1(2*NTP1*(NPMGS+1)*NPPSCUB),  ZS2(2*NTP1*(NPMGS+1)*NPPSCUB),
     R   ZT1(2*NTP1*(NPMGS+1)*NPPSCUB),  ZT2(2*NTP1*(NPMGS+1)*NPPSCUB)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ZEPS = RC1M13
         LOOPMAX = 100
         KTET = NT1 * (NMGAUS+1)
         KP = (KN - 1) / KTET + 1
C
         DO 1 J1=1,KN
           I1 = (J1 - 1) / KTET + 1
           ZSIGMN(J1) = 0.
           ZSIGMX(J1) = 1.
           ZFMIN(J1)  = PSIAXE - PISO(I1)
           ZFMAX(J1)  = PSIBND - PISO(I1)
           IC(J1)     = 1
           PSIG(J1)   = .5
C
           IF (ZFMIN(J1) .EQ. 0.) THEN
             ZSIGMX(J1) = 0.
             PSIG(J1)   = 0.
             ZFMAX(J1)  = ZFMIN(J1)
           ELSE IF (ZFMAX(J1) .EQ. 0.) THEN
             ZSIGMN(J1) = 1.
             PSIG(J1)   = 1.
             ZFMIN(J1)  = ZFMAX(J1)
           ENDIF
 1       CONTINUE
C     
         CALL BOUND(KN,PTET,ZBND)
C
         CALL RESETI(ID,KN,1)
         DO 2 JT = 1,NT1
            DO 2 JG=1,KN
               IF (ID(JG).EQ.0) GOTO 2
               IT0(JG) = JT-1
               IF (PTET(JG).LE.CT(JT)) ID(JG)  = 0
 2       CONTINUE
C
         DO 3 JG=1,KN
            IF (IT0(JG).LT.1) IT0(JG) = 1
            ZT1(JG) = CT(IT0(JG))
 3          ZT2(JG) = CT(IT0(JG)+1)
C     
         DO 10 J10=1,LOOPMAX
C     
         IF (J10.NE.1) THEN
            CALL RESETI(ID,KN,1)
            DO 4 J4=1,KN
            IF (ZSIGMN(J4).GE.ZS1(J4).AND.ZSIGMX(J4).LE.ZS2(J4))
     +              ID(J4) = 0
 4          CONTINUE
            ITEST = ISSUM(KN,ID,1)
            IF (ITEST .EQ. 0) GOTO 7
         ENDIF

         CALL RESETI(ID,KN,1)
         DO 5 JS = 1,NS1
            DO 5 JG=1,KN
               IF (ID(JG).EQ.0) GOTO 5
               IS0(JG) = JS-1
               IF (PSIG(JG).LE.CSIG(JS)) ID(JG)  = 0
 5       CONTINUE
C
         DO 6 J6=1,KN
           IF (IS0(J6) .LT. 1)  IS0(J6) = 1
C
           ZS1(J6) = CSIG(IS0(J6))
           ZS2(J6) = CSIG(IS0(J6)+1)
 6       CONTINUE
C
         CALL PSICEL(IS0,IT0,KN,2*NTP1*(NPMGS+1)*NPPSCUB,ZPCEL,CPSICL)
C
 7       CONTINUE
C
         CALL BASIS1(KN,2*NTP1*(NPMGS+1)*NPPSCUB,ZS1,ZS2,ZT1,ZT2,
     +               PSIG,PTET,ZF)
C
         DO 8 J8=1,KN
           ZPSI(J8) = ZPCEL(J8, 1) * ZF(J8, 1) +
     +                ZPCEL(J8, 2) * ZF(J8, 2) +
     +                ZPCEL(J8, 3) * ZF(J8, 3) +
     +                ZPCEL(J8, 4) * ZF(J8, 4) +
     +                ZPCEL(J8, 5) * ZF(J8, 5) +
     +                ZPCEL(J8, 6) * ZF(J8, 6) +
     +                ZPCEL(J8, 7) * ZF(J8, 7) +
     +                ZPCEL(J8, 8) * ZF(J8, 8) +
     +                ZPCEL(J8, 9) * ZF(J8, 9) +
     +                ZPCEL(J8,10) * ZF(J8,10) +
     +                ZPCEL(J8,11) * ZF(J8,11) +
     +                ZPCEL(J8,12) * ZF(J8,12) +
     +                ZPCEL(J8,13) * ZF(J8,13) +
     +                ZPCEL(J8,14) * ZF(J8,14) +
     +                ZPCEL(J8,15) * ZF(J8,15) +
     +                ZPCEL(J8,16) * ZF(J8,16)
 8      CONTINUE
C
         DO 9 J9=1,KN
           IF (IC(J9) .EQ. 0) GOTO 9
C
           IF (ZFMIN(J9) .NE. 0.) THEN
             I9 = (J9 - 1) / KTET + 1
             ZFMID = ZPSI(J9) - PISO(I9)
          ELSE
             ZFMID = 0.
           ENDIF
C
           ZTEST = .5 * ABS(ZSIGMX(J9) - ZSIGMN(J9))
C
           IF (ZTEST .LE. ZEPS) THEN
              IC(J9)  = 0
              IF (PSIG(J9).GT.1.) PSIG(J9) = 1.
           ELSE IF (ZTEST .GT. ZEPS) THEN
             IF (ZFMIN(J9) * ZFMID .LE. 0.) THEN
               ZSIGMX(J9) = PSIG(J9)
               ZFMAX(J9)  = ZFMID
             ELSE
               ZSIGMN(J9) = PSIG(J9)
               ZFMIN(J9)  = ZFMID
             ENDIF
             PSIG(J9) = .5 * (ZSIGMX(J9) + ZSIGMN(J9))
           ENDIF
C
 9       CONTINUE
C
         ITEST = ISSUM(KN,IC,1)
C
         IF (ITEST .EQ. 0) GOTO 11
 10      CONTINUE
         PRINT*,'CUBRT NOT CONVERGED'
         STOP
 11      CONTINUE
C
         RETURN
         END
C*DECK C2SU03
         SUBROUTINE GAUSS(KPTS,PGAUSX,PGAUSH)
C        ####################################
C
C                                        AUTHORS:
C                                        O.SAUTER,  CRPP-EPFL
C
C  POINTS AND WEIGHTS FOR THE GAUSS INTEGRATION FORMULA FOR KPTS POINTS
C  UP TO 6 POINTS AND 20, I.E. INTEGRATING EXACTLY POLYNOMIAL OF DEGREE
C  UP TO 12 AND 40.
C  TAKEN OUT FROM ABRAMOWITZ
C
C
         INCLUDE 'DECLAR.inc'
         DIMENSION   PGAUSX(KPTS+1),PGAUSH(KPTS+1)
C
         IF (KPTS .EQ. 1) THEN
C
C     1-POINTS FORMULA ==> TRAPEZOIDAL RULE
            PGAUSX(1) = -1.0
            PGAUSX(2) =  1.0
C
            PGAUSH(1) =  1.0
            PGAUSH(2) =  1.0
C
            KPTS = 2
C
         ELSE IF (KPTS .EQ. 2) THEN
C
C     2-POINTS FORMULA
            PGAUSX(1) = -0.577350269189626
            PGAUSX(2) = +0.577350269189626
C
            PGAUSH(1) =  1.0
            PGAUSH(2) =  1.0
C
         ELSE IF (KPTS .EQ. 3) THEN
C
C     3-POINTS FORMULA
            PGAUSX(1) = -0.774596669241483
            PGAUSX(2) =  0.0
            PGAUSX(3) = +0.774596669241483
C
            PGAUSH(1) =  0.555555555555556
            PGAUSH(2) =  0.888888888888889
            PGAUSH(3) =  0.555555555555556
C
         ELSE IF (KPTS .EQ. 4) THEN
C
C     4-POINTS FORMULA
            PGAUSX(1) = -0.861136311594053
            PGAUSX(2) = -0.339981043584856
            PGAUSX(3) =  0.339981043584856
            PGAUSX(4) =  0.861136311594053
C
            PGAUSH(1) = 0.347854845137454
            PGAUSH(2) = 0.652145154862546
            PGAUSH(3) = 0.652145154862546
            PGAUSH(4) = 0.347854845137454
C
         ELSE IF (KPTS .EQ. 5) THEN
C
C     5-POINTS FORMULA
            PGAUSX(1) = -0.906179845938664
            PGAUSX(2) = -0.538469310105683
            PGAUSX(3) =  0.0
            PGAUSX(4) = +0.538469310105683
            PGAUSX(5) = +0.906179845938664
C
            PGAUSH(1) =  0.236926885056189
            PGAUSH(2) =  0.478628670499366
            PGAUSH(3) =  0.568888888888889
            PGAUSH(4) =  0.478628670499366
            PGAUSH(5) =  0.236926885056189
C
         ELSE IF (KPTS .EQ. 6) THEN
C
C     6-POINTS FORMULA
            PGAUSX(1) = -0.932469514203152
            PGAUSX(2) = -0.661209386466265
            PGAUSX(3) = -0.238619186083197
            PGAUSX(4) = +0.238619186083197
            PGAUSX(5) = +0.661209386466265
            PGAUSX(6) = +0.932469514203152
C
            PGAUSH(1) =  0.171324492379170
            PGAUSH(2) =  0.360761573048139
            PGAUSH(3) =  0.467913934572691
            PGAUSH(4) =  0.467913934572691
            PGAUSH(5) =  0.360761573048139
            PGAUSH(6) =  0.171324492379170
C
         ELSE IF (KPTS .EQ. 7) THEN
C
C     7-POINTS FORMULA
            PGAUSX(1) = -0.949107912342759
            PGAUSX(2) = -0.741531185599394
            PGAUSX(3) = -0.405845151377397
            PGAUSX(4) = 0.
            PGAUSX(5) = 0.405845151377397
            PGAUSX(6) = 0.741531185599394
            PGAUSX(7) = 0.949107912342759
C
            PGAUSH(1) =  0.129484966168870
            PGAUSH(2) =  0.279705391489277
            PGAUSH(3) =  0.381830050505119
            PGAUSH(4) =  0.417959183673469
            PGAUSH(5) =  0.381830050505119
            PGAUSH(6) =  0.279705391489277
            PGAUSH(7) =  0.129484966168870
C
         ELSE IF (KPTS .EQ. 8) THEN
C
C     8-POINTS FORMULA
            PGAUSX(1) = -0.960289856497536
            PGAUSX(2) = -0.796666477413627
            PGAUSX(3) = -0.525532409916329
            PGAUSX(4) = -0.183434642495650
            PGAUSX(5) = 0.183434642495650
            PGAUSX(6) = 0.525532409916329
            PGAUSX(7) = 0.796666477413627
            PGAUSX(8) = 0.960289856497536
C
            PGAUSH(1) =  0.101228536290376
            PGAUSH(2) =  0.222381034453374
            PGAUSH(3) =  0.313706645877887
            PGAUSH(4) =  0.362683783378362
            PGAUSH(5) =  0.362683783378362
            PGAUSH(6) =  0.313706645877887
            PGAUSH(7) =  0.222381034453374
            PGAUSH(8) =  0.101228536290376
C
         ELSE IF (KPTS .EQ. 9) THEN
C
C    9-POINTS FORMULA 
            PGAUSX(1) = -0.968160239507626
            PGAUSX(2) = -0.836031107326636
            PGAUSX(3) = -0.613371432700590
            PGAUSX(4) = -0.324253423403809
            PGAUSX(5) =  0.
            PGAUSX(6) = 0.324253423403809
            PGAUSX(7) = 0.613371432700590
            PGAUSX(8) = 0.836031107326636
            PGAUSX(9) = 0.968160239507626
C
            PGAUSH(1) = 0.081274388361574
            PGAUSH(2) = 0.180648160694857 
            PGAUSH(3) = 0.260610696402935
            PGAUSH(4) = 0.312347077040003
            PGAUSH(5) = 0.330239355001260
            PGAUSH(6) = PGAUSH(4)
            PGAUSH(7) = PGAUSH(3)
            PGAUSH(8) = PGAUSH(2)
            PGAUSH(9) = PGAUSH(1)
C
         ELSE IF (KPTS .EQ. 10) THEN
C
C    10-POINTS FORMULA 
            PGAUSX(1)  = -0.973906528517172
            PGAUSX(2)  = -0.865063366688985
            PGAUSX(3)  = -0.679409568299024
            PGAUSX(4)  = -0.433395394129247
            PGAUSX(5)  = -0.148874338981631
            PGAUSX(6)  = 0.148874338981631
            PGAUSX(7)  = 0.433395394129247
            PGAUSX(8)  = 0.679409568299024
            PGAUSX(9)  = 0.865063366688985
            PGAUSX(10) = 0.973906528517172
C
            PGAUSH(1)  = 0.066671344308688
            PGAUSH(2)  = 0.149451349150581
            PGAUSH(3)  = 0.219086362515982
            PGAUSH(4)  = 0.269266719309996
            PGAUSH(5)  = 0.295524224714753
            PGAUSH(6)  = PGAUSH(5)
            PGAUSH(7)  = PGAUSH(4)
            PGAUSH(8)  = PGAUSH(3)
            PGAUSH(9)  = PGAUSH(2)
            PGAUSH(10) = PGAUSH(1)
C
         ELSE IF (KPTS .EQ. 20) THEN
C
C    20-POINTS FORMULA
            PGAUSX(1) = -0.993128599185094
            PGAUSX(2) = -0.963971927277913
            PGAUSX(3) = -0.912234428251325
            PGAUSX(4) = -0.839116971822218
            PGAUSX(5) = -0.746331906460150
            PGAUSX(6) = -0.636053680726515
            PGAUSX(7) = -0.510867001950827
            PGAUSX(8) = -0.373706088715419
            PGAUSX(9) = -0.227785851141645
            PGAUSX(10)= -0.076526521133497
            PGAUSX(11)=  0.076526521133497
            PGAUSX(12)=  0.227785851141645
            PGAUSX(13)=  0.373706088715419
            PGAUSX(14)=  0.510867001950827
            PGAUSX(15)=  0.636053680726515
            PGAUSX(16)=  0.746331906460150
            PGAUSX(17)=  0.839116971822218
            PGAUSX(18)=  0.912234428251325
            PGAUSX(19)=  0.963971927277913
            PGAUSX(20)=  0.993128599185094
C
            PGAUSH(1) =  0.017614007139152
            PGAUSH(2) =  0.040601429800386
            PGAUSH(3) =  0.062672048334109
            PGAUSH(4) =  0.083276741576704
            PGAUSH(5) =  0.101930119817240
            PGAUSH(6) =  0.118194531961518
            PGAUSH(7) =  0.131688638449176
            PGAUSH(8) =  0.142096109318382
            PGAUSH(9) =  0.149172986472603
            PGAUSH(10)=  0.152753387130725
            PGAUSH(11)=  0.152753387130725
            PGAUSH(12)=  0.149172986472603
            PGAUSH(13)=  0.142096109318382
            PGAUSH(14)=  0.131688638449176
            PGAUSH(15)=  0.118194531961518
            PGAUSH(16)=  0.101930119817240
            PGAUSH(17)=  0.083276741576704
            PGAUSH(18)=  0.062672048334109
            PGAUSH(19)=  0.040601429800386
            PGAUSH(20)=  0.017614007139152
C
         ENDIF
C
         DO 1 J1=1,KPTS
C
         PGAUSH(J1) = .5 * PGAUSH(J1)
C
    1    CONTINUE  
C
         RETURN
         END
C*DECK C2SU04
C*CALL PROCESS
         SUBROUTINE RMRAD(KN,PSIAXE,PSIBND,PT,PAR,PSIG,PTET,KINC)
C        ########################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
*  C2SU04  COMPUTE INTERSECTIONS OF CONSTANT FLUX SURFACE WITH        *
*          Z = RZMAG LINE.
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     I   IS0(2*NPISO+NPT),       IT0(2*NPISO+NPT),   
     I   IC(2*NPISO+NPT),        IC1(2*NPISO+NPT),
     R   PAR(KN),                PSIG(*),      
     R   PTET(*),                ZBND(2*NPISO+NPT),       
     R   ZF(2*NPISO+NPT,16),     ZFMIN(2*NPISO+NPT), 
     R   ZFMAX(2*NPISO+NPT),     ZPCEL(2*NPISO+NPT,16),
     R   ZPSI(2*NPISO+NPT),      ZRHO(2*NPISO+NPT),    
     R   ZRTMID(2*NPISO+NPT),    ZRTMIN(2*NPISO+NPT),
     R   ZRTMAX(2*NPISO+NPT),    ZSIGMA(2*NPISO+NPT),    
     R   ZS1(2*NPISO+NPT),       ZS2(2*NPISO+NPT),       
     R   ZTETA(2*NPISO+NPT),     
     R   ZT1(2*NPISO+NPT),       
     R   ZT2(2*NPISO+NPT)     
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ZEPS = RC1M13
C
         BPS( 1) = RMAG
         BPS(12) = RZMAG
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
         CALL BOUND(1,PT,ZBNDT0)
C
         DO 1 J1=1,KN
C
         ZRTMIN(J1) = 0.
         ZRTMAX(J1) = 1.
         ZFMIN(J1)  = PSIAXE - PSIISO(J1)
         ZFMAX(J1)  = PSIBND - PSIISO(J1)
         IC(J1)     = 1
C
         IF (ZFMIN(J1) .EQ. 0.) THEN
C
            ZRTMAX(J1) = 0.
            ZFMAX(J1)  = ZFMIN(J1)
C
         ELSE IF (ZFMAX(J1) .EQ. 0.) THEN
C
            ZRTMIN(J1) = 1.
            ZFMIN(J1)  = ZFMAX(J1)
C
         ENDIF
C
    1    CONTINUE
C
         BPS( 1) = R0
         BPS(12) = RZ0
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
    2    CONTINUE
C
         DO 3 J3=1,KN
C
         ZRTMID(J3) = .5 * (ZRTMAX(J3) + ZRTMIN(J3))
         ZR         = ZRTMID(J3) * ZBNDT0 * COS(PT) + RMAG
         ZZ         = ZRTMID(J3) * ZBNDT0 * SIN(PT) + RZMAG
         ZRHO(J3)   = SQRT((ZR - R0)**2 + (ZZ - RZ0)**2)
         ZTETA(J3)  = ATAN2(ZZ - RZ0,ZR - R0)
C
         IF (ZTETA(J3) .LT. CT(1)) ZTETA(J3) = ZTETA(J3) + 2. * CPI
C
    3    CONTINUE
C
         CALL BOUND(KN,ZTETA,ZBND)
C
         CALL RESETI(IC1,KN,1)
         DO 5 JT = 1,NT1
            DO 5 JG=1,KN
               IF (IC1(JG).EQ.0) GOTO 5
               IT0(JG) = JT-1
               IF (ZTETA(JG).LE.CT(JT)) IC1(JG) = 0
 5       CONTINUE
         CALL RESETI(IC1,KN,1)
         DO 6 JS = 1,NS1
            DO 6 JG=1,KN
               IF (IC1(JG).EQ.0) GOTO 6
               ZSIGMA(JG) = ZRHO(JG) / ZBND(JG)
               IS0(JG) = JS-1
               IF (ZSIGMA(JG).LE.CSIG(JS)) IC1(JG) = 0
 6       CONTINUE
C
         DO 7 J7=1,KN
           IF (IS0(J7) .GT. NS) IS0(J7) = NS
           IF (IS0(J7) .LT. 1)  IS0(J7) = 1
           IF (IT0(J7) .GT. NT) IT0(J7) = NT
           IF (IT0(J7) .LT. 1)  IT0(J7) = 1
C
           ZS1(J7) = CSIG(IS0(J7))
           ZS2(J7) = CSIG(IS0(J7)+1)
           ZT1(J7) = CT(IT0(J7))
 7         ZT2(J7) = CT(IT0(J7)+1)
C
         CALL BASIS1(KN,2*NPISO+NPT,ZS1,ZS2,ZT1,ZT2,ZSIGMA,ZTETA,ZF)
         CALL PSICEL(IS0,IT0,KN,2*NPISO+NPT,ZPCEL,CPSICL)
C
         DO 8 J8=1,KN
C
         ZPSI(J8) = ZPCEL(J8, 1) * ZF(J8, 1) +
     +              ZPCEL(J8, 2) * ZF(J8, 2) +
     +              ZPCEL(J8, 3) * ZF(J8, 3) +
     +              ZPCEL(J8, 4) * ZF(J8, 4) +
     +              ZPCEL(J8, 5) * ZF(J8, 5) +
     +              ZPCEL(J8, 6) * ZF(J8, 6) +
     +              ZPCEL(J8, 7) * ZF(J8, 7) +
     +              ZPCEL(J8, 8) * ZF(J8, 8) +
     +              ZPCEL(J8, 9) * ZF(J8, 9) +
     +              ZPCEL(J8,10) * ZF(J8,10) +
     +              ZPCEL(J8,11) * ZF(J8,11) +
     +              ZPCEL(J8,12) * ZF(J8,12) +
     +              ZPCEL(J8,13) * ZF(J8,13) +
     +              ZPCEL(J8,14) * ZF(J8,14) +
     +              ZPCEL(J8,15) * ZF(J8,15) +
     +              ZPCEL(J8,16) * ZF(J8,16)
C
 8       CONTINUE
C
         DO 9 J9=1,KN
C
         IF (IC(J9) .EQ. 0) GOTO 9
C
         IF (ZFMIN(J9) .NE. 0.) THEN
C
            ZFMID = ZPSI(J9) - PSIISO(J9)
C
         ELSE
C
            ZFMID = 0.
C
         ENDIF
C
         ZTEST = .5 * ABS(ZRTMAX(J9) - ZRTMIN(J9))
C
         IF (ZTEST .LE. ZEPS) THEN
C
            IC(J9)               = 0
            PTET((J9-1)*KINC+1)  = ZTETA(J9)
            PSIG((J9-1)*KINC+1)  = ZSIGMA(J9)
            PAR(J9)              = ZRTMID(J9)
C
         ELSE IF (ZTEST .GT. ZEPS) THEN
C
            IF (ZFMIN(J9) * ZFMID .LE. 0.) THEN
C
               ZRTMAX(J9) = ZRTMID(J9)
               ZFMAX(J9)  = ZFMID
C
            ELSE
C
               ZRTMIN(J9) = ZRTMID(J9)
               ZFMIN(J9)  = ZFMID
C
            ENDIF
         ENDIF
C
 9       CONTINUE
C
         ITEST = ISSUM(KN,IC,1)
C
         IF (ITEST .GE. 1) GOTO 2
C
         RETURN
C
         END
C*DECK C2SX01
C*CALL PROCESS
         SUBROUTINE BOUND(KN,PT,PR)
C        ##########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
C     INTRODUCE JITERMX MAX. NUMBER OF NEWTON ITERATIONS TO PREVENT
C     BEING TRAPPED IN INFINITE LOOPS. APLET 5/1/96
c
***********************************************************************
*                                                                     *
*  C2SX01  DEFINE PLASMA SURFACE                                      *
*             NSURF = 1   SOLOVEV PLASMA SURFACE                      *
*                     2   INTOR-LIKE SURFACE                          *
*                     3   TCV-LIKE SURFACE                            *
*                     4   XPOINT SURFACE                              *
*                     5   HOCTUPOLE SURFACE                           *
*                     6   EXPERIMENTAL POINTS                         *
*                     7   SINE AND COSINE DECOMPOSITION               *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         PARAMETER (NPKNMAX = 12*NPT+2*NPPSCUB*(NPMGS+1)*NTP1)
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R   PT(KN),                 PR(KN),
     R   IC(NPKNMAX),            ZF(NPKNMAX),
     R   ZPT(NPKNMAX),           ZPT0(NPKNMAX),
     R   ZCOSPT(NPKNMAX),        ZSINPT(NPKNMAX),
     R   ZRHOJ(NPKNMAX),         ZRHOJP(NPKNMAX)
C
***********************************************************************
*                                                                     *
* 0.1. FUNCTIONS FOR SOLOVEV SURFACE                                  *
*                                                                     *
***********************************************************************
C
         ZF1(Z)  = BPS(2) * SQRT(1 + 2 * BPS(4) * ZCOS) - BPS(1)
         ZF1P(Z) = - BPS(2) * BPS(4) * ZSIN /
     /             SQRT(1 + 2 * BPS(4) * ZCOS)
         ZF2(Z)  = BPS(5) * BPS(4) * BPS(2) * ZSIN /
     /             SQRT(1 + 2 * BPS(4) * ZCOS)
         ZF2P(Z) = BPS(5) * BPS(4) * BPS(2) * (ZCOS + BPS(4) *
     *             ZSIN**2 / (1 + 2 * BPS(4) * ZCOS)) /
     /             SQRT(1 + 2 * BPS(4) * ZCOS)
C
***********************************************************************
*                                                                     *
* 0.2. FUNCTIONS FOR INTOR-LIKE SURFACE                               *
*                                                                     *
***********************************************************************
C
         ZFI1 (Z) = BPS(3) + BPS(6) * (1. + BPS(8) * ZCOS) * 
     *              COS(Z + BPS(7) * ZSIN - BPS(9) * ZSIN2)
         ZFI1P(Z) = - BPS(6) * (COS(Z+BPS(7)*ZSIN-BPS(9)*ZSIN2) * 
     *                BPS(8) * ZSIN + (1. + BPS(7) * ZCOS - 
     -                2.* BPS(9)*ZCOS2) * SIN(Z + BPS(7) * ZSIN -
     -                BPS(9) * ZSIN2) * (1. + BPS(8) * ZCOS))
         ZFI2 (Z) = BPS(5) * BPS(6) * SIN(Z + BPS(10) * ZSIN2)
         ZFI2P(Z) = BPS(5) * BPS(6) * COS(Z + BPS(10) * ZSIN2) *
     *              (1. + 2. * BPS(10) * ZCOS2)
C
***********************************************************************
*                                                                     *
* 0.3. FUNCTIONS FOR TCV-LIKE SURFACE                                 *
*                                                                     *
***********************************************************************
C
         INTS(Z) = CPI * INT(2. * Z / CPI) 
         ZFV(Z)  = ZEXV + (1. - ZEXV) * ZEXA
         ZDFV(Z) = BPS(10) * (1. - ZEXV) * ZEXA *
     *             SIGN(Z2,ZCOS * ZSIN) / ZSIN**2
C
         T(Z)     = INTPI2 + (Z - INTPI2) * ZFV(Z)
         DT(Z)    = (Z - INTPI2) * ZDFV(Z) +  ZFV(Z)
C
         S(Z)     = COS(ZTPT + BPS(7) * SIN(ZTPT))
         DS(Z)    = DT(Z) * (1. + BPS(7) * COS(ZTPT)) *
     *              SIN(ZTPT + BPS(7) * SIN(ZTPT))
         BEAN(Z)  = 1. + BPS(8) * ZCOS
         TRIP(Z)  = 1. + BPS(11) * COS(2.5 * CPI * ZSIN)
C
         RRVR(Z)  = BPS(6) * BEAN(Z) * TRIP(Z) * S(Z) + BPS(3)
C
         RRVRP(Z) = - BPS(6) * ((BPS(8) * ZSIN * TRIP(Z) + 
     +                           2.5 * CPI * BPS(11) * ZCOS * BEAN(Z) * 
     *                           SIN(2.5 * CPI * ZSIN)) * S(Z) +
     +                          BEAN(Z) * TRIP(Z) * DS(Z))
C
         RZVR(Z)  = BPS(5) * BPS(6) * ZSIN
         RZVRP(Z) = BPS(5) * BPS(6) * ZCOS
C
***********************************************************************
*                                                                     *
* 0.4. FUNCTIONS FOR X-POINT SURFACE                                  *
*                                                                     *
***********************************************************************
C
         ZEX(Z)   = EXP( BPS(7) * LOG(SIN(.5 * (Z - BPS(9)))**2 +
     +                                BPS(8)))
         ZRHO(Z)  = BPS(6) * (1. + BPS(10) * BPS(11) /
     /                                 (ZEX(Z) + BPS(11)))
         ZRHOP(Z) = .5 * BPS(6) * BPS(7) * BPS(10) * BPS(11) * 
     *              ZEX(Z) * SIN(Z - BPS(9)) /
     /              ((SIN(.5 * (Z - BPS(9)))**2 + BPS(8)) *
     *              (ZEX(Z) + BPS(11))**2)
         ZFT(Z)   = (ZRHO(Z) * SIN(Z) * BPS(5) - BPS(12)) * ZCOS -
     -              (ZRHO(Z) * COS(Z + BPS(13) * SIN(Z)) * 
     *              (1. + BPS(14) * COS(Z)) + BPS(3)) * ZSIN
         ZDFDT(Z) = ZRHOP(Z) * SIN(Z) * BPS(5) * ZCOS + 
     +              ZRHO(Z)  * COS(Z) * BPS(5) * ZCOS -
     -              ZRHOP(Z) * COS(Z + BPS(13) * SIN(Z)) * 
     *              (1. + BPS(14) * COS(Z)) * ZSIN + ZRHO(Z) *
     *              (SIN(Z + BPS(13) * SIN(Z)) * 
     *               (1. + BPS(13) * COS(Z)) *
     *               (1. + BPS(14) * COS(Z)) + 
     +               COS(Z + BPS(13) * SIN(Z)) *
     *               BPS(14) * SIN(Z)) * ZSIN
C
***********************************************************************
*                                                                     *
* 0.5. FUNCTIONS FOR HOCTUPOLE TYPE  SURFACE                          *
*                                                                     *
***********************************************************************
C
          ZFIR(Z)= BPS(3) + BPS(6) * ZCOS
     .            /(1.-BPS(5)*(COS(4.*(Z+BPS(8)))-1.))**BPS(7)

          ZFIZ(Z)= BPS(6)*ZSIN
     .            /(1.-BPS(5)*(COS(4.*(Z+BPS(8)))-1.))**BPS(7)

          ZFIRP(Z)=-ZFIZ(Z)
     .           -4.*BPS(5)*BPS(6)*SIN(4.*(Z+BPS(8)))*ZCOS*BPS(7)
     .            /(1.-BPS(5)*(COS(4.*(Z+BPS(8)))-1.))**(BPS(7)+1.)

          ZFIZP(Z)=(ZFIR(Z)-BPS(3))
     .              -4.*BPS(5)*BPS(6)*SIN(4.*(Z+BPS(8)))*ZSIN*BPS(7)
     .            /(1.-BPS(5)*(COS(4.*(Z+BPS(8)))-1.))**(BPS(7)+1.)

C
***********************************************************************
*                                                                     *
* 0.7. FUNCTIONS FOR FOURRIER TRANSFORM                               *
*                                                                     *
***********************************************************************
          RPTET(COSTF) = BPS(3) + RHOPF*COSTF
          ZPTET(SINTF) = RHOPF*SINTF + BPS(6) - BPS(12)
          RPTETP(COSTF,SINTF) = RHOPFP*COSTF - RHOPF*SINTF
          ZPTETP(COSTF,SINTF) = RHOPFP*SINTF + RHOPF*COSTF
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
          ZEPS = RC1M13
c
C MAX NUMBER OF NEWTON ITERATIONS
c
         JITERMX = 40
C
         IF (KN .GT. NPKNMAX) THEN
C
            PRINT*,'DIMENSION OF LOCAL ARRAYS IS TO SMALL: KN= ',KN
            STOP
C
         ENDIF
C
         GO TO (100,200,300,400,500,600,700) NSURF
C
 100     CONTINUE
C
***********************************************************************
*                                                                     *
* 1. SOLOVEV PLASMA SURFACE                                           *
*                                                                     *
***********************************************************************
*                                                                     *
* 1.1. INITIALIZATION                                                 *
*                                                                     *
***********************************************************************
C
C
         DO 110 J = 1,KN
C
         ZPT(J) = PT(J)
         IC(J)  = 1
         ZSIN   = SIN(ZPT(J))
         ZCOS   = COS(ZPT(J))
         ZF(J)  = SIN(PT(J)) * ZF1(ZPT(J)) - COS(PT(J)) * ZF2(ZPT(J))
C
 110     CONTINUE
C
***********************************************************************
*                                                                     *
* 1.2. NEWTON'S RULE                                                  *
*                                                                     *
***********************************************************************
C
         DO 120 JITER=1,JITERMX
C
         DO 121 J = 1,KN
         ZSIN = SIN(ZPT(J))
         ZCOS = COS(ZPT(J))
         ZFP  = SIN(PT(J))*ZF1P(ZPT(J))-COS(PT(J))*ZF2P(ZPT(J))
         ZEP  = ZEPS-ABS(ZF(J))
C
         IF (ZEP .LE. 0.) ZPT(J) = ZPT(J) - ZF(J) / ZFP
C
         ZSIN  = SIN(ZPT(J))
         ZCOS  = COS(ZPT(J))
         ZF(J) = SIN(PT(J))*ZF1(ZPT(J))-COS(PT(J))*ZF2(ZPT(J))
C
         IF (ZEP .GT. 0.) IC(J) = 0
C
 121     CONTINUE
C
         ITEST = ISSUM(KN,IC,1)
C
         IF (ITEST .LE. 0) GO TO 122
C
 120     CONTINUE
C
         WRITE(*,9001)
C
 122     CONTINUE
C
***********************************************************************
*                                                                     *
* 1.3. DEFINE R(ZPT) AND Z(ZPT)                                       *
*                                                                     *
***********************************************************************
C     DEFINE R(ZPT) AND Z(ZPT)
C
         DO 130 J = 1,KN
         ZSIN = SIN(ZPT(J))
         ZCOS = COS(ZPT(J))
         ZR   = ZF1(ZPT(J))
         ZZ   = ZF2(ZPT(J))
C
         PR(J) = SQRT(ZR**2 + ZZ**2)
 130     CONTINUE
C
         RETURN
C
C-----------------------------------------------------------------------
C
 200     CONTINUE
C
***********************************************************************
*                                                                     *
* 2. INTOR-LIKE PLASMA SURFACE                                        *
*                                                                     *
***********************************************************************
*                                                                     *
* 2.1. INITIALIZATION                                                 *
*                                                                     *
***********************************************************************
C
         DO 210 J = 1,KN
         ZPT(J) = PT(J)
         IC(J) = 1
         ZSIN  = SIN(ZPT(J))
         ZCOS  = COS(ZPT(J))
         ZSIN2 = SIN(2.*ZPT(J))
         ZF(J) = SIN(PT(J))*ZFI1(ZPT(J))-COS(PT(J))*ZFI2(ZPT(J))
 210     CONTINUE
C
***********************************************************************
*                                                                     *
* 2.2. NEWTON'S RULE                                                  *
*                                                                     *
***********************************************************************
C
         DO 220 JITER=1,JITERMX
C
         DO 221 J = 1,KN
         ZSIN  = SIN(ZPT(J))
         ZCOS  = COS(ZPT(J))
         ZSIN2 = SIN(2.*ZPT(J))
         ZCOS2 = COS(2.*ZPT(J))
         ZFP   = SIN(PT(J))*ZFI1P(ZPT(J))-COS(PT(J))*ZFI2P(ZPT(J))
         ZEP   = ZEPS-ABS(ZF(J))
C
         IF (ZEP .LE. 0.) ZPT(J) = ZPT(J) - ZF(J) / ZFP
C
         ZSIN  = SIN(ZPT(J))
         ZCOS  = COS(ZPT(J))
         ZSIN2 = SIN(2.*ZPT(J))
         ZF(J) = SIN(PT(J))*ZFI1(ZPT(J))-COS(PT(J))*ZFI2(ZPT(J))
C
         IF (ZEP .GT. 0.) IC(J) = 0
C
 221     CONTINUE
C
         ITEST = ISSUM(KN,IC,1)
C
         IF (ITEST.le.0) GO TO 222
C
 220     CONTINUE
         WRITE(*,9001)
C
 222     CONTINUE
C
         DO 223 J = 1,KN
         ZSIN  = SIN(ZPT(J))
         ZCOS  = COS(ZPT(J))
         ZSIN2 = SIN(2.*ZPT(J))
         ZR    = ZFI1(ZPT(J))
         ZZ    = ZFI2(ZPT(J))
C
         PR(J) = SQRT(ZR*ZR+ZZ*ZZ)
C         WRITE(*,'("J,RR,ZZ,PR",I4,1P3E12.4)')J,ZR,ZZ,PR(J)
 223   CONTINUE

C
         RETURN
C
C-----------------------------------------------------------------------
C
 300     CONTINUE
C
***********************************************************************
*                                                                     *
* 3. TCV-LIKE PLASMA SURFACE                                          *
*                                                                     *
***********************************************************************
*                                                                     *
* 3.1. INITIALIZATION                                                 *
*                                                                     *
***********************************************************************
C
         ZARG = BPS(9) * BPS(10)
C
         IF (ZARG .LE. 50.) THEN
C
            ZEXV = EXP(- ZARG)
C
         ELSE
C
            ZEXV = 0.
C
         ENDIF
C
         DO 311 J311=1,KN
C
         IF (PT(J311) .GE. 0.) THEN
C
            ZPT(J311) = PT(J311) - 2.*CPI * INT(.5*PT(J311)/CPI)
C
         ELSE
C
            ZPT(J311) = PT(J311) + 2.*CPI * (1. - INT(.5*PT(J311)/CPI))
C
         ENDIF
C
 311     CONTINUE  
C
         DO 312 J312=1,KN
C
         IF (ZPT(J312) .GT. CPI) ZPT(J312) = 2. * CPI - ZPT(J312)
C
         ZPT0(J312) = ZPT(J312)
         IC(J312)   = 1
         ZCOS       = COS(ZPT(J312))
         ZSIN       = SIN(ZPT(J312))
         INTPI2     = INTS(ZPT(J312))
C
         IF (ZSIN .NE. 0.) THEN
C
            Z1       = MIN(50.,BPS(10) * ABS(ZCOS / ZSIN))
            ZEXA     = EXP(- Z1)
            ZTPT     = T(ZPT(J312))
            ZF(J312) = ZSIN * RRVR(ZPT(J312)) - 
     -                 ZCOS * RZVR(ZPT(J312))
C
         ELSE
C
            ZF(J312) = 0.
C
         ENDIF
C
 312     CONTINUE
C
***********************************************************************
*                                                                     *
* 3.2. NEWTON'S RULE                                                  *
*                                                                     *
***********************************************************************
C
         DO 320 JITER=1,JITERMX
C
         DO 321 J321=1,KN
C
         ZCOS   = COS(ZPT(J321))
         ZSIN   = SIN(ZPT(J321))
         INTPI2 = INTS(ZPT(J321))
C
         IF (ZSIN .NE. 0.) THEN
C
C
            Z1   = MIN(50.,BPS(10) * ABS(ZCOS / ZSIN))
            Z2   = 1.
            ZEXA = EXP(- Z1)
            ZTPT = T(ZPT(J321))
            ZFP  = SIN(ZPT0(J321)) * RRVRP(ZPT(J321)) -
     -             COS(ZPT0(J321)) * RZVRP(ZPT(J321))
            ZEP  = ZEPS - ABS(ZF(J321))
C
         ELSE
C
C           DUMMY VALUE FOR ZFP TO PERMIT VECTORZATION            
C
            ZFP  = 1.
            ZEP  = ZEPS
C
         ENDIF
C
         IF (ZEP .LE. 0.) ZPT(J321) = ZPT(J321) - ZF(J321) / ZFP
         IF (ZEP .GT. 0.) IC(J321)  = 0
C
 321     CONTINUE
C
         DO 322 J322=1,KN
C
         IF (ZPT(J322) .GT. CPI) ZPT(J322) = 2. * CPI - ZPT(J322)
C
         ZCOS   = COS(ZPT(J322))
         ZSIN   = SIN(ZPT(J322))
         INTPI2 = INTS(ZPT(J322))
C
         IF (ZSIN .NE. 0.) THEN
C
            Z1       = MIN(50.,BPS(10) * ABS(ZCOS / ZSIN))
            ZEXA     = EXP(- Z1)
            ZTPT     = T(ZPT(J322))
            ZF(J322) = SIN(ZPT0(J322)) * RRVR(ZPT(J322)) -
     -                 COS(ZPT0(J322)) * RZVR(ZPT(J322))
C
         ELSE
C
            ZF(J322) = 0.
C
         ENDIF
C
 322     CONTINUE
C
         ITEST = ISSUM(KN,IC,1)
C
         IF (ITEST.le.0) GO TO 323
C
 320     CONTINUE
C
         WRITE(*,9001)
C
 323     CONTINUE
C
         DO 324 J324=1,KN
C
         ZSIN   = SIN(ZPT(J324))
         ZCOS   = COS(ZPT(J324))
         INTPI2 = INTS(ZPT(J324))
C
         IF (ZSIN .NE. 0.) THEN
C
            Z1   = MIN(50.,BPS(10) * ABS(ZCOS / ZSIN))
            ZEXA = EXP(- Z1)
            ZTPT = T(ZPT(J324))
            ZR   = RRVR(ZPT(J324))
            ZZ   = RZVR(ZPT(J324))
C
         ELSE
C
            ZR = BPS(6) * BEAN(ZPT(J324)) * TRIP(ZPT(J324)) * ZCOS + 
     +           BPS(3)
            ZZ = 0.
C
         ENDIF
C
         PR(J324) = SQRT(ZR**2 + ZZ**2)
C
 324     CONTINUE
C
         RETURN
C
C-----------------------------------------------------------------------
C
 400     CONTINUE
C
***********************************************************************
*                                                                     *
* 4. X-POINT SURFACE                                                  *
*                                                                     *
***********************************************************************
*                                                                     *
* 4.1. INITIALIZATION                                                 *
*                                                                     *
***********************************************************************
C
         DO 410 J = 1,KN
         ZPT(J) = PT(J)
         ZSIN  = SIN(PT(J))
         ZCOS  = COS(PT(J))
         ZF(J) = ZFT(PT(J))
         IC(J) = 1
 410     CONTINUE
C
***********************************************************************
*                                                                     *
* 4.2. NEWTON'S RULE                                                  *
*                                                                     *
***********************************************************************
C
         DO 420 JITER=1,JITERMX
C
         DO 421 J = 1,KN
         ZSIN = SIN(PT(J))
         ZCOS = COS(PT(J))
         ZEP  = ZEPS-ABS(ZF(J))
C
         IF (ZEP .LE. 0.) ZPT(J) = ZPT(J) - ZF(J) / ZDFDT(ZPT(J))
C
         ZF(J) = ZFT(ZPT(J))
C
         IF (ZEP .GT. 0.) IC(J) = 0
C
 421     CONTINUE
C
         ITEST = ISSUM(KN,IC,1)
C
         IF (ITEST.le.0) GO TO 422
C
 420     CONTINUE
C
         WRITE(*,9001)
C
 422     CONTINUE
C
         DO 423 J = 1,KN
         ZRO  = ZRHO(ZPT(J))
         ZZZ  = ZPT(J)
         ZR   = BPS(3)+ZRO*COS(ZZZ+BPS(13)*SIN(ZZZ))
     *         *(1.0+BPS(14)*COS(ZZZ))
         ZZ = ZRO*SIN(ZPT(J))*BPS(5)-BPS(12)
         PR(J) = SQRT(ZR**2+ZZ**2)
 423     CONTINUE
C
         RETURN
C
C-----------------------------------------------------------------------
C
 500     CONTINUE
C
***********************************************************************
*                                                                     *
* 5. OCTUPOLE TYPE PLASMA SURFACE                                     *
*                                                                     *
***********************************************************************
*                                                                     *
* 5.1. INITIALIZATION                                                 *
*                                                                     *
***********************************************************************
C
         DO 510 J = 1,KN
         ZPT(J) = PT(J)
         IC(J) = 1
         ZSIN  = SIN(ZPT(J))
         ZCOS  = COS(ZPT(J))
         ZF(J) = SIN(PT(J))*ZFIR(ZPT(J))-COS(PT(J))*ZFIZ(ZPT(J))
 510     CONTINUE
C
***********************************************************************
*                                                                     *
* 5.2. NEWTON'S RULE                                                  *
*                                                                     *
***********************************************************************
C
         DO 520 JITER=1,JITERMX
C
         DO 521 J = 1,KN
         ZSIN = SIN(ZPT(J))
         ZCOS = COS(ZPT(J))
         ZFP  = SIN(PT(J))*ZFIRP(ZPT(J))-COS(PT(J))*ZFIZP(ZPT(J))
         ZEP  = ZEPS-ABS(ZF(J))
C
         IF (ZEP .LE. 0.) ZPT(J) = ZPT(J) - ZF(J) / ZFP
C
         ZSIN  = SIN(ZPT(J))
         ZCOS  = COS(ZPT(J))
         ZF(J) = SIN(PT(J))*ZFIR(ZPT(J))-COS(PT(J))*ZFIZ(ZPT(J))
C
         IF (ZEP .GT. 0.) IC(J) = 0
C
 521     CONTINUE
C
         ITEST = ISSUM(KN,IC,1)
C
         IF (ITEST.le.0) GO TO 522
C
 520     CONTINUE
C
         WRITE(*,9001)
C
 522     CONTINUE
C
         DO 523 J = 1,KN
         ZSIN = SIN(ZPT(J))
         ZCOS = COS(ZPT(J))
         ZR   = ZFIR(ZPT(J))
         ZZ   = ZFIZ(ZPT(J))
C
         PR(J) = SQRT(ZR*ZR+ZZ*ZZ)
 523     CONTINUE
C
         RETURN
C
C-----------------------------------------------------------------------
C
 600     CONTINUE
C
***********************************************************************
*                                                                     *
* 6. PLASMA SURFACE GIVEN BY ARRAY OF R AND Z NODES                   *
*                                                                     *
***********************************************************************
C
         DO 601 J=1,KN
            ZPT(J) = PT(J)
            IF (ZPT(J) .LT. TETBPS(1,1))    ZPT(J) = ZPT(J) + 2. * CPI
            IF (ZPT(J) .GT. TETBPS(NBPS,1)) ZPT(J) = ZPT(J) - 2. * CPI
 601     CONTINUE
C
         DO 602 J=1,KN
         DO 602 JJ=1,NBPS-1
            IF ((ZPT(J)-TETBPS(JJ,1)) * (ZPT(J)-TETBPS(JJ+1,1)).LE.0.)
     &      IC(J) = JJ
 602     CONTINUE
C
         DO 603 J=1,KN
            ZH = TETBPS(IC(J)+1,1) - TETBPS(IC(J),1)
            ZA = (TETBPS(IC(J)+1,1) - ZPT(J)) / ZH
            ZB = (ZPT(J) - TETBPS(IC(J),1)) / ZH
            ZC = (ZA+1)*(ZA-1)*ZH * (TETBPS(IC(J)+1,1) - ZPT(J)) / 6.
            ZD = (ZB+1)*(ZB-1)*ZH * (ZPT(J) - TETBPS(IC(J),1)) / 6.
C
C NONCONFORMAL WALL: ONLY PLASMA SURFACE CONSIDERED
            ZR = ZA * RRBPS(IC(J),1)  + ZB * RRBPS(IC(J)+1,1) +
     +           ZC * D2RBPS(IC(J),1) + ZD * D2RBPS(IC(J)+1,1) - BPS(1)
            ZZ = ZA * RZBPS(IC(J),1)  + ZB * RZBPS(IC(J)+1,1) +
     +           ZC * D2ZBPS(IC(J),1) + ZD * D2ZBPS(IC(J)+1,1) - BPS(12)
C
            PR(J) = SQRT(ZR**2+ZZ**2)
 603     CONTINUE
C
         RETURN
C
C-----------------------------------------------------------------------
C
 700     CONTINUE
C
***********************************************************************
*                                                                     *
* 7. SURFACE CALCULATED BY REAL FOURRIER COEFFICIENTS                 *
*                                                                     *
***********************************************************************
*                                                                     *
* 7.1. INITIALIZATION                                                 *
*                                                                     *
***********************************************************************
C        
         DO 710 J = 1,KN
            ZPT(J) = PT(J)
            IC(J) = 1
            ZCOSPT(J) = COS(PT(J))
            ZSINPT(J) = SIN(PT(J))
 710     CONTINUE
C
***********************************************************************
*                                                                     *
* 7.2. NEWTON'S RULE                                                  *
*      TWO WAYS OF COMPUTING TO SAVE CPU ON VECTORIZED MACHINE        *
*      SCALAR MACHINE SHOULD PROBABLY USE THE 1ST OPTION (K INNER LOOP)
*                                                                     *
***********************************************************************
C
         DO 720 JITER=1,JITERMX
C
           IF (KN .LE. NFOURPB) THEN
C
C     K MOST INNER LOOP
C
             DO 721 J = 1,KN
C
               IF (IC(J) .EQ. 0) GO TO 721
C
               RHOPF  = ALZERO
               RHOPFP = 0.
               DO K=1,NFOURPB
                 RHOPF = RHOPF + BPSCOS(K)*COS(K*ZPT(J))
     +             + BPSSIN(K)*SIN(K*ZPT(J))
                 RHOPFP = RHOPFP - K*BPSCOS(K)*SIN(K*ZPT(J))
     +             + K*BPSSIN(K)*COS(K*ZPT(J))
               ENDDO
C
               ZSIN = SIN(ZPT(J))
               ZCOS = COS(ZPT(J))
               ZRHOJ(J) = RHOPF
C
               ZFP = ZSINPT(J)*RPTETP(ZCOS,ZSIN)
     +           - ZCOSPT(J)*ZPTETP(ZCOS,ZSIN)
               ZFJ = ZSINPT(J)*RPTET(ZCOS)-ZCOSPT(J)*ZPTET(ZSIN)
               ZEP = ZEPS-ABS(ZFJ)
C
               IF (ZEP .LE. 0.) THEN
                 ZPT(J) = ZPT(J) - ZFJ / ZFP
               ELSE
                 IC(J) = 0
               ENDIF
C
 721         CONTINUE
C
           ELSE
C
C     J MOST INNER LOOP
C
             DO J=1,KN
               ZRHOJ (J) = ALZERO
               ZRHOJP(J) = 0.0
             ENDDO
C
             DO 722 K=1,NFOURPB
               ZK = K
               ZAK = BPSCOS(K)
               ZBK = BPSSIN(K)
C
               DO J=1,KN
                 ZRHOJ (J) = ZRHOJ (J) + ZAK*COS(ZK*ZPT(J))
     +             + ZBK*SIN(ZK*ZPT(J))
                 ZRHOJP(J) = ZRHOJP(J) - ZK*ZAK*SIN(ZK*ZPT(J))
     +             + ZK*ZBK*COS(ZK*ZPT(J))
               ENDDO
C
 722         CONTINUE
C
             DO J=1,KN
C
               RHOPF  = ZRHOJ (J)
               RHOPFP = ZRHOJP(J)
               ZSIN = SIN(ZPT(J))
               ZCOS = COS(ZPT(J))
               ZFP = ZSINPT(J)*RPTETP(ZCOS,ZSIN)
     +           - ZCOSPT(J)*ZPTETP(ZCOS,ZSIN)
               ZFJ = ZSINPT(J)*RPTET(ZCOS)-ZCOSPT(J)*ZPTET(ZSIN)
               ZEP = ZEPS - ABS(ZFJ)
C
               IF (ZEP .LE. 0.) THEN
                 ZPT(J) = ZPT(J) - ZFJ / ZFP
               ELSE
                 IC(J) = 0
               ENDIF
C
             ENDDO
C
           ENDIF
C
           ITEST = ISSUM(KN,IC,1)
C        
           IF (ITEST.LE.0) GO TO 723
C
 720     CONTINUE
C
         WRITE(*,9001)
C
 723     CONTINUE
C
         DO J=1,KN
           RHOPF = ZRHOJ(J)
           ZR = RPTET(COS(ZPT(J)))
           ZZ = ZPTET(SIN(ZPT(J)))
           PR(J) = SQRT(ZR*ZR+ZZ*ZZ)
         ENDDO
C
         RETURN
C
C-----------------------------------------------------------------------
C
 9001    FORMAT(/,'WARNING: MAX # OF NEWTON ITERATIONS REACHED IN'
     +     ,' BOUND (INCREASE ZEPS OR JITERMX)')
         END
C*DECK C2SX01
C*CALL PROCESS
         SUBROUTINE BOUNDNW(KN,JNW,PT,PR,PZ,PC)
C        ##########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
*  BOUNDNW: COMPUTE (R,Z) FOR PLASMA SURFACE AND WALL SHAPES AT GIVEN *
*           PT ANGLES                                                 *
*  LIU YQ, DECEMBER 12, 2002                                          *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         PARAMETER (NPKNMAX = 12*NPT+2*NPPSCUB*(NPMGS+1)*NTP1)
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     R   PT(KN),                 PR(KN), PZ(KN), PC(KN),
     R   IC(NPKNMAX),            ZF(NPKNMAX),
     R   ZPT(NPKNMAX),           ZPT0(NPKNMAX),
     R   ZCOSPT(NPKNMAX),        ZSINPT(NPKNMAX),
     R   ZRHOJ(NPKNMAX),         ZRHOJP(NPKNMAX)

         DO 601 J=1,KN
            ZPT(J) = PT(J)
            IF (ZPT(J) .LT. TETBPS(1,JNW))    ZPT(J) = ZPT(J) + 2.*CPI
            IF (ZPT(J) .GT. TETBPS(NBPS,JNW)) ZPT(J) = ZPT(J) - 2.*CPI
 601     CONTINUE
C
         DO 602 J=1,KN
         DO 602 JJ=1,NBPS-1
            IF ((ZPT(J)-TETBPS(JJ,JNW))*(ZPT(J)-TETBPS(JJ+1,JNW)).LE.0.)
     &      IC(J) = JJ
 602     CONTINUE
C
         DO 603 J=1,KN
            ZH = TETBPS(IC(J)+1,JNW) - TETBPS(IC(J),JNW)
            ZA = (TETBPS(IC(J)+1,JNW) - ZPT(J)) / ZH
            ZB = (ZPT(J) - TETBPS(IC(J),JNW)) / ZH
            ZC = (ZA+1)*(ZA-1)*ZH * (TETBPS(IC(J)+1,JNW) - ZPT(J))/6.
            ZD = (ZB+1)*(ZB-1)*ZH * (ZPT(J) - TETBPS(IC(J),JNW)) / 6.
C
            PR(J) = ZA * RRBPS(IC(J),JNW)  + ZB * RRBPS(IC(J)+1,JNW) +
     +       ZC * D2RBPS(IC(J),JNW) + ZD * D2RBPS(IC(J)+1,JNW)
            PZ(J) = ZA * RZBPS(IC(J),JNW)  + ZB * RZBPS(IC(J)+1,JNW) +
     +       ZC * D2ZBPS(IC(J),JNW) + ZD * D2ZBPS(IC(J)+1,JNW)
            PC(J) = ZA * CNDRZ(IC(J),JNW)  + ZB * CNDRZ(IC(J)+1,JNW) +
     +       ZC * D2CBPS(IC(J),JNW) + ZD * D2CBPS(IC(J)+1,JNW)
C
 603     CONTINUE
C
         RETURN
         END

         SUBROUTINE BOUNDDT2(KN,PT,PC,PCT)
C        #################################
***********************************************************************
*  BOUNDDT2: COMPUTE PCT=D(PC)/D(PT)                                  *
*  LIU YQ, FEB. 08, 2012                                              *
***********************************************************************

         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMSOL.inc'
         PARAMETER (NPKNMAX = NPMGS*NTP1)

         DIMENSION
     R   PT(KN),                 PC(KN),             PCT(KN),
     R   ZPTA(NPKNMAX+2),        
     R   ZPTB(NPKNMAX+2),        ZPCB(NPKNMAX+2),
     R   PT2(NPKNMAX,4),         PC2(NPKNMAX,4)

         ZEPS = 1.E-3

C        MAKE SURE THAT PT IS IN INCREASING ORDER 
         DO J=1,KN
            ZPTA(J) = PT(J)
         ENDDO
         
         DO J=2,KN
            IF (ZPTA(J).LT.ZPTA(J-1)) ZPTA(J) = ZPTA(J) + RC2PI
         ENDDO

C        DEFINE ADDITIONAL POINTS PT2
         DO J2=1,KN
            PT2(J2,1) = ZPTA(J2) - 2. * ZEPS
            PT2(J2,2) = ZPTA(J2) -      ZEPS
            PT2(J2,3) = ZPTA(J2) +      ZEPS
            PT2(J2,4) = ZPTA(J2) + 2. * ZEPS
         ENDDO

C        ADD TWO ENDS POINTS TO (PT,PC) ARRAYS
         DO J=1,KN
            ZPTB(J+1) = ZPTA(J)
            ZPCB(J+1) = PC(J)
         ENDDO
         ZDELTA     = CPI - (ZPTA(KN)-ZPTA(1))*.5
         ZPTB(1)    = ZPTA(1)  - ZDELTA
         ZPTB(KN+2) = ZPTA(KN) + ZDELTA
         ZPCB(1)    = (PC(1)+PC(KN))*.5
         ZPCB(KN+2) = ZPCB(1)

C        USE SPLINE TO FIND PC AT OTHER POINTS
         CALL BOUNDNW2(KN,ZPTB,ZPCB,PT2(1,1),PC2(1,1))
         CALL BOUNDNW2(KN,ZPTB,ZPCB,PT2(1,2),PC2(1,2))
         CALL BOUNDNW2(KN,ZPTB,ZPCB,PT2(1,3),PC2(1,3))
         CALL BOUNDNW2(KN,ZPTB,ZPCB,PT2(1,4),PC2(1,4))

C        COMPUTE D(PC)/D(PT)
         DO J3=1,KN
            PCT(J3) = (PC2(J3,1) + 8*(PC2(J3,3) - PC2(J3,2)) -
     &                 PC2(J3,4)) / (12. * ZEPS)
         ENDDO

         RETURN
         END

         SUBROUTINE BOUNDDCDT(KN,ZS,PT,PC,PCN,PCT,ZPCS)
C        #################################
***********************************************************************
*  BOUNDDT2: COMPUTE PCT=D(PC)/D(PT) AT EACH SURFACE IN VACUUM        *
*  LIU YQ, FEB. 27, 2012                                              *
***********************************************************************

         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMESH.inc'
         PARAMETER (NPKNMAX = NPMGS*NTP1)

         DIMENSION
     R   PT(KN),      PC(KN),    PCN(KN),          PCT(KN),    ZPCS(KN),
     R   ZPTA(NPKNMAX+2),        ZPCA(NPKNMAX+2),       
     R   ZPTB(NPKNMAX+2),        ZPCB(NPKNMAX+2),
     R   PT2(NPKNMAX,4),         PC2(NPKNMAX,4)

         DO J=1,KN
            ZPCS(J) = 0.
         ENDDO

         IF (KMETHOD.LE.0) THEN
            DO J=1,KN
               PCN(J) = PC(J)
            ENDDO
            RETURN
         ENDIF    

         ZEPS = 1.E-3

C        MAKE SURE THAT PT IS IN INCREASING ORDER 
         DO J=1,KN
            ZPTA(J) = PT(J)
         ENDDO
         DO J=2,KN
            IF (ZPTA(J).LT.ZPTA(J-1)) ZPTA(J) = ZPTA(J) + RC2PI
         ENDDO

C        MAKE SURE THAT PC IS IN INCREASING ORDER 
         DO J=1,KN
            ZPCA(J) = PC(J)
         ENDDO
         DO J=2,KN
            IF (ZPCA(J).LT.ZPCA(J-1)) ZPCA(J) = ZPCA(J) + RC2PI
         ENDDO

C        DEFINE ADDITIONAL POINTS PT2
         DO J2=1,KN
            PT2(J2,1) = ZPTA(J2) - 2. * ZEPS
            PT2(J2,2) = ZPTA(J2) -      ZEPS
            PT2(J2,3) = ZPTA(J2) +      ZEPS
            PT2(J2,4) = ZPTA(J2) + 2. * ZEPS
         ENDDO

C        ADD TWO END POINTS TO (PT,PC) ARRAYS
         DO J=1,KN
            ZPTB(J+1) = ZPTA(J)
            ZPCB(J+1) = ZPCA(J)
         ENDDO
         ZDELTA     = CPI - (ZPTA(KN)-ZPTA(1))*.5
         ZPTB(1)    = ZPTA(1)  - ZDELTA
         ZPTB(KN+2) = ZPTA(KN) + ZDELTA
         ZDELTA     = CPI - (ZPCA(KN)-ZPCA(1))*.5
         ZPCB(1)    = ZPCA(1)  - ZDELTA
         ZPCB(KN+2) = ZPCA(KN) + ZDELTA

C        REDEFINE PC
         IF (KMETHOD.EQ.1) THEN
         ZALF = 1./(2.-EXP(-(ZS-1.)**2/ZDEL**2))         
         ZDELTA = (ZPCB(1)+ZPCB(KN+2))*.5
         DO J=1,KN+2
            IF (ZPCB(J).LE.ZDELTA) J0=J
         ENDDO
         ZPC1 = ZPCB(1)
         ZPC2 = ZPCB(KN+2)
         ZPC3 = ZPCB(J0)
         DO J=2,J0-1
            ZPCB(J) =ZPC1+((ZPCB(J)-ZPC1)/(ZPC3-ZPC1))**ZALF*(ZPC3-ZPC1)
         ENDDO
         DO J=J0+1,KN+1
            ZPCB(J) =ZPC2-((ZPC2-ZPCB(J))/(ZPC2-ZPC3))**ZALF*(ZPC2-ZPC3)
         ENDDO
         ENDIF

         IF (KMETHOD.EQ.2) THEN
         ZPC1 = ZPCB(1)
         DO J=1,KN+2
            ZPCB(J) = (ZPCB(J) + ZFRC*(ZS-1.)*(ZPTB(J)-ZPTB(1)+ZPC1))/
     &                (1.+ZFRC*(ZS-1.))
         ENDDO
         ENDIF
         
         IF (KMETHOD.EQ.3) THEN
         ZTMP = 0.
         ZPC1 = ZPCB(1)
         DO J=1,NPTS
            ZTMP = ZTMP + 2*ZDEL - ZDEL*(EXP((ZPTB(1)-ZPTS(J))/ZDEL) + 
     &             EXP((ZPTS(J)-ZPTB(KN+2))/ZDEL))
         ENDDO
         ZTMP = RC2PI/ZTMP
         DO J1=2,KN+1
            ZTMP2 = 0.
            DO J=1,NPTS
               IF (ZPTB(J1).LE.ZPTS(J)) THEN
                  ZTMP2 = ZTMP2 + ZDEL*EXP(-ZPTS(J)/ZDEL)*
     &                    (EXP(ZPTB(J1)/ZDEL)-EXP(ZPTB(1)/ZDEL))
               ELSE
                  ZTMP2 = ZTMP2 + 2.*ZDEL - ZDEL*
     &            (EXP((ZPTB(1)-ZPTS(J))/ZDEL) + 
     &             EXP((ZPTS(J)-ZPTB(J1))/ZDEL))
               ENDIF
            ENDDO
            ZTMP3= ZPCB(J1)
            ZPCB(J1) = ZPC1 + ((ZTMP3-ZPC1)+ZFRC*(ZS-1.)*ZTMP*ZTMP2)/
     &                 (1.+ZFRC*(ZS-1.))
            ZPCS(J1-1) = (-(ZTMP3-ZPC1)+ZTMP*ZTMP2)*ZFRC/
     &                 (1.+ZFRC*(ZS-1.))**2
         ENDDO
         ENDIF

C        USE SPLINE TO FIND PC AT OTHER POINTS
         CALL BOUNDNW3(KN,ZPTB,ZPCB,PT2(1,1),PC2(1,1))
         CALL BOUNDNW3(KN,ZPTB,ZPCB,PT2(1,2),PC2(1,2))
         CALL BOUNDNW3(KN,ZPTB,ZPCB,PT2(1,3),PC2(1,3))
         CALL BOUNDNW3(KN,ZPTB,ZPCB,PT2(1,4),PC2(1,4))

C        COMPUTE D(PC)/D(PT)
         DO J3=1,KN
            PCT(J3) = (PC2(J3,1) + 8*(PC2(J3,3) - PC2(J3,2)) -
     &                 PC2(J3,4)) / (12. * ZEPS)
         ENDDO


C        SET NEW VALUES FOR CHIISO         
         DO J=1,KN
            PCN(J) = ZPCB(J+1)
         ENDDO

         RETURN
         END

         SUBROUTINE BOUNDNW2(KN,PT,PC,PT2,PC2)
C        #####################################
***********************************************************************
*  BOUNDNW2: USE SPLINE TO PERFORM INTEPOLATION                       *
*            FOR PERIODIC FUNCTIONS ONLY                              *
*            INPUT  ARRAY: PC(PT)                                     *
*            OUTPUT ARRAY: PC2(PT2)                                   *
*  LIU YQ, FEB. 08, 2012                                              *
***********************************************************************
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMSOL.inc'
         PARAMETER (NPKNMAX = NPMGS*NTP1)

         DIMENSION
     R   PT(KN+2),               PC(KN+2),
     R   PT2(KN),                PC2(KN),
     R   ZPT(NPKNMAX),           ZA1(NPKNMAX+2),
     R   ZB1(NPKNMAX+2),         ZC1(NPKNMAX+2),
     R   D2PC(NPKNMAX+2),          
     I   IC(NPKNMAX),            IT(NPKNMAX)           

         DO J=1,KN
            ZPT(J) = PT2(J)
            IF (ZPT(J) .LT. PT(1))    ZPT(J) = ZPT(J) + 2.*CPI
            IF (ZPT(J) .GT. PT(KN+2)) ZPT(J) = ZPT(J) - 2.*CPI
         ENDDO

         CALL SPLCY(PT,PC,KN+1,RC2PI,D2PC,ZA1,ZB1,ZC1)
         D2PC(KN+2) = D2PC(1)

         CALL RESETI(IC,KN,1)
         DO JG=1,KN
         DO JT=1,KN+2
            IF (IC(JG).EQ.0) GOTO 2
            IT(JG) = JT-1
            IF (PT(JT).GE.ZPT(JG)) IC(JG) = 0
 2          CONTINUE
         ENDDO
         ENDDO

         DO J=1,KN
            ICHI = IT(J)
            IF (ICHI.LT.1) ICHI = 1
            IF (ICHI.GT.KN+1) ICHI = KN+1
            ZH = PT(ICHI+1) - PT(ICHI)
            ZA = (PT(ICHI+1) - ZPT(J)) / ZH
            ZB = (ZPT(J) - PT(ICHI)) / ZH
            ZC = (ZA + 1)*(ZA - 1)*ZH*(PT(ICHI+1)-ZPT(J))/6.
            ZD = (ZB + 1)*(ZB - 1)*ZH*(ZPT(J)-PT(ICHI))/6.
         
            PC2(J) = ZA*PC(ICHI) + ZB*PC(ICHI+1) + 
     &               ZC*D2PC(ICHI) + ZD*D2PC(ICHI+1)
         ENDDO

         RETURN
         END

         SUBROUTINE BOUNDNW3(KN,PT,PC,PT2,PC2)
C        #####################################
***********************************************************************
*  BOUNDNW2: USE SPLINE TO PERFORM INTEPOLATION                       *
*            FOR NON-PERIODIC FUNCTIONS ONLY                          *
*            INPUT  ARRAY: PC(PT)                                     *
*            OUTPUT ARRAY: PC2(PT2)                                   *
*  LIU YQ, FEB. 08, 2012                                              *
***********************************************************************
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMSOL.inc'
         PARAMETER (NPKNMAX = NPMGS*NTP1)

         DIMENSION
     R   PT(KN+2),               PC(KN+2),
     R   PT2(KN),                PC2(KN),
     R   ZPT(NPKNMAX),          
     R   ZWORK(NPKNMAX+2),       ZWORK1(NPKNMAX+2),
     R   D2PC(NPKNMAX+2),          
     I   IC(NPKNMAX),            IT(NPKNMAX)           

         DO J=1,KN
            ZPT(J) = PT2(J)
            IF (ZPT(J) .LT. PT(1))    ZPT(J) = ZPT(J) + 2.*CPI
            IF (ZPT(J) .GT. PT(KN+2)) ZPT(J) = ZPT(J) - 2.*CPI
         ENDDO

         CALL SPLINE(PT,PC,KN+2,D2PC,ZWORK,ZWORK1)
         CALL RESETI(IC,KN,1)
         DO JG=1,KN
         DO JT=1,KN+2
            IF (IC(JG).EQ.0) GOTO 2
            IT(JG) = JT-1
            IF (PT(JT).GE.ZPT(JG)) IC(JG) = 0
 2          CONTINUE
         ENDDO
         ENDDO

         DO J=1,KN
            ICHI = IT(J)
            IF (ICHI.LT.1) ICHI = 1
            IF (ICHI.GT.KN+1) ICHI = KN+1
            ZH = PT(ICHI+1) - PT(ICHI)
            ZA = (PT(ICHI+1) - ZPT(J)) / ZH
            ZB = (ZPT(J) - PT(ICHI)) / ZH
            ZC = (ZA + 1)*(ZA - 1)*ZH*(PT(ICHI+1)-ZPT(J))/6.
            ZD = (ZB + 1)*(ZB - 1)*ZH*(ZPT(J)-PT(ICHI))/6.
            PC2(J) = ZA*PC(ICHI) + ZB*PC(ICHI+1) + 
     &               ZC*D2PC(ICHI) + ZD*D2PC(ICHI+1)
         ENDDO

         RETURN
         END

C*DECK C2SX02
C*CALL PROCESS
         SUBROUTINE BNDINP
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SX02 READ EXPERIMENTAL EQUILIBRIUM BOUNDARY AND/OR EXPERIMENTAL   *
*        PROFILES                                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         CALL IODISK(33)
C
         RETURN
         END
C*DECK C2SX03
C*CALL PROCESS
         SUBROUTINE BNDSPL
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SX03 INITIALIZE QUANTITIES USED FOR SPLINE INTERPOLATION OF       *
*        BOUNDARY GIVEN BY EXPERIMENTAL POINTS                        *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMSOL.inc'
C
         DIMENSION
     &   ZA1(NPBPS), ZB1(NPBPS), ZC1(NPBPS),
     &   ZWKSP(NPBPS)
         DIMENSION
     &   JWKSP(NPBPS)
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C   COMPUTE THETA(R,Z) OF LINES PASSING THROUGH (BPS(1),BPS(12)) AND
C   THE (R,Z) POINTS KNOWN ON THE PLASMA BOUNDARY
C
C NONCONFORMAL WALL: LOOP IN INVERSED ORDER, TO KEEP TETBPS FOR PLASMA 
c SURFACE
         DO 8888 J8=NWBPS,1,-1
         DO 1 J1=1,NBPS-1
            ZDR = RRBPS(J1,J8) - BPS( 1)
            ZDZ = RZBPS(J1,J8) - BPS(12)
C
            IF (ZDR .EQ. 0.) THEN
               IF (ZDZ .GE. 0.) TETBPS(J1,J8) =  .5 * CPI
               IF (ZDZ .LT. 0.) TETBPS(J1,J8) = 1.5 * CPI
            ELSE
               TETBPS(J1,J8) = ATAN2(ZDZ,ZDR)
            ENDIF
C 
            IF (TETBPS(J1,J8).LT.0.) TETBPS(J1,J8)=TETBPS(J1,J8)+2.*CPI
    1    CONTINUE
C
C        SORT TETBPS, RRBPS, RZBPS AND CNDRZ SO THAT ELEMENTS OF TETBPS ARE IN 
C        INCREASING ORDER
C
         CALL SORT(NBPS-1,TETBPS(1,J8),RRBPS(1,J8),
     *             RZBPS(1,J8),CNDRZ(1,J8),ZWKSP,JWKSP)
CSYM
C     SYMMETRIZE BOUNDARY SO THAT LOWER HALF = UPPER HALF
C
         IF (NSYM.EQ.1) THEN
            IBP = 1
            DO 2 J2=1,NBPS-1
               IF (TETBPS(J2,J8).GT.CPI) GOTO 3
               IBP = IBP+1
   2        CONTINUE
   3        CONTINUE
C
            NBPS = 2 * (IBP - 1) + 1
            IF (NBPS.GT.NPBPS) THEN
               PRINT*,'NPBPS SMALLER THAN NBPS'
               STOP
            ENDIF
C
            DO 4 J4=IBP,NBPS-1
               TETBPS(J4,J8) = 2*CPI-TETBPS(NBPS-J4,J8)
               RRBPS(J4,J8)  = RRBPS(NBPS-J4,J8)
               RZBPS(J4,J8)  = -RZBPS(NBPS-J4,J8)
               CNDRZ(J4,J8)  = CNDRZ(NBPS-J4,J8)
   4        CONTINUE
         ENDIF
C
C        COMPUTE SPLINE COEFFICIENTS OF RRBPS AND RZBPS VERSUS TETBPS
C        WITH CYCLIC BOUNDARY CONDITIONS
C
         CALL SPLCY(TETBPS(1,J8),RRBPS(1,J8),NBPS-1,RC2PI,
     *              D2RBPS(1,J8),ZA1,ZB1,ZC1)
         CALL SPLCY(TETBPS(1,J8),RZBPS(1,J8),NBPS-1,RC2PI,
     *              D2ZBPS(1,J8),ZA1,ZB1,ZC1)
         CALL SPLCY(TETBPS(1,J8),CNDRZ(1,J8),NBPS-1,RC2PI,
     *              D2CBPS(1,J8),ZA1,ZB1,ZC1)
C
         TETBPS(NBPS,J8) = TETBPS(1,J8) + 2. * CPI
         RRBPS(NBPS,J8)  = RRBPS(1,J8)
         RZBPS(NBPS,J8)  = RZBPS(1,J8)
         CNDRZ(NBPS,J8)  = CNDRZ(1,J8)
         D2RBPS(NBPS,J8) = D2RBPS(1,J8)
         D2ZBPS(NBPS,J8) = D2ZBPS(1,J8)
         D2CBPS(NBPS,J8) = D2CBPS(1,J8)

 8888    CONTINUE
C
         RETURN
         END
C*DECK C2SX04
C*CALL PROCESS
         SUBROUTINE SUBSZ
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SX04 SHIFT Z POINTS OF BOUNDARY SO THAT SYMMETRY PLANE OF         *
*        EQUILIBRIUM IS Z=0                                           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMESH.inc'
C
C NONCONFORMAL WALL
         DO 1 J2=1,NWBPS
         DO 1 J1=1,NBPS
            RZBPS(J1,J2) = RZBPS(J1,J2) - RZ0
   1     CONTINUE
C
         RETURN
         END
C*DECK C2SX05
C*CALL PROCESS
         SUBROUTINE RZBOUND
C        ##################
C                                        AUTHORS:
C                                        O.SAUTER,  CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SX05  COMPUTE (R,Z) OF PLASMA BOUNDARY ON NBPSOUT POINTS          *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMBND.inc'
C
         DIMENSION
     R      ZTETBPS(NPBPS), ZBND(NPBPS)
C
C-----------------------------------------------------------------------
C
CL       1. SET THETA MESH
C
         ZDTHETA = 2. * CPI / DFLOAT(NBPSOUT-1)
         DO I=1,NBPSOUT
            ZTETBPS(I) = DFLOAT(I-1) * ZDTHETA
         END DO
C
CL       2. COMPUTE R,Z
C
         CALL BOUND(NBPSOUT,ZTETBPS,ZBND)
         DO I=1,NBPSOUT
            RRBPSOU(I) = BPS(1)  + ZBND(I) * COS(ZTETBPS(I))
            RZBPSOU(I) = BPS(12) + ZBND(I) * SIN(ZTETBPS(I))
         END DO
C
         RETURN
         END
C*DECK C2SY01
C*CALL PROCESS
         SUBROUTINE BASIS1(KN,KPN,PS1,PS2,PT1,PT2,PS,PT,PF)
C        ##################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SY01 EVALUATES THE 2-D CUBIC HERMITE BASIS FUNCTIONS AT           *
*        (SIGMA,THETA) = (PS,PT)                                      *
*                                                                     *
* THE 2-D BASIS FUNCTIONS A PRODUCTS OF 1-D FUNCTIONS DEFINED IN      *
* HERMITE.inc. IN CHEASE, THE BASIS FUNCTIONS ARE FUNCTIONS OF        *
* SIGMA AND THETA                                                     *
*                                                                     *
* PS1 AND PS2 ARE THE CELL LIMITS IN THE SIGMA DIRECTION              *
* PT1 AND PT2 ARE THE CELL LIMITS IN THE THETA DIRECTION              *
*                                                                     *
* KN IS THE NUMBER OF POINTS WHERE THE VALUES ARE COMPUTED            *
*                                                                     *
* THE RESULT VECTOR IS CONTAINED BY PF                                *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PS(KPN),    PS1(KPN),    PS2(KPN),    PT(KPN),
     R   PT1(KPN),   PT2(KPN),    PF(KPN,16)
C
         INCLUDE 'HERMIT.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         DO 1 J1=1,KN
C     
         ZFN1S = FN1(PS(J1),PS1(J1),PS2(J1))
         ZFN2S = FN2(PS(J1),PS1(J1),PS2(J1))
         ZFN3S = FN3(PS(J1),PS1(J1),PS2(J1))
         ZFN4S = FN4(PS(J1),PS1(J1),PS2(J1))
C     
         ZFN1T = FN1(PT(J1),PT1(J1),PT2(J1))
         ZFN2T = FN2(PT(J1),PT1(J1),PT2(J1))
         ZFN3T = FN3(PT(J1),PT1(J1),PT2(J1))
         ZFN4T = FN4(PT(J1),PT1(J1),PT2(J1))
         
         PF(J1, 1) = ZFN1S * ZFN1T
         PF(J1, 2) = ZFN2S * ZFN1T
         PF(J1, 3) = ZFN1S * ZFN2T
         PF(J1, 4) = ZFN2S * ZFN2T
C
         PF(J1, 5) = ZFN3S * ZFN1T
         PF(J1, 6) = ZFN4S * ZFN1T
         PF(J1, 7) = ZFN3S * ZFN2T
         PF(J1, 8) = ZFN4S * ZFN2T
C
         PF(J1, 9) = ZFN1S * ZFN3T
         PF(J1,10) = ZFN2S * ZFN3T
         PF(J1,11) = ZFN1S * ZFN4T
         PF(J1,12) = ZFN2S * ZFN4T
C
         PF(J1,13) = ZFN3S * ZFN3T
         PF(J1,14) = ZFN4S * ZFN3T
         PF(J1,15) = ZFN3S * ZFN4T
         PF(J1,16) = ZFN4S * ZFN4T
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SY02
C*CALL PROCESS
         SUBROUTINE BASIS2(KN,KPN,PS1,PS2,PT1,PT2,PS,PT,PDFDS,PDFDT)
C        ###########################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SY02 IS SIMILAR TO C2SY01. IT COMPUTES FIRST DERIVATIVES OF 2-D   *
*        CUBIC HERMITE BASIS FUNCTIONS                                *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PS(KPN),    PS1(KPN),    PS2(KPN),       PT(KPN),
     R   PT1(KPN),   PT2(KPN),    PDFDS(KPN,16),  PDFDT(KPN,16)
C
         INCLUDE 'HERMIT.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         DO 1 J1=1,KN
C
         ZFN1S = FN1(PS(J1),PS1(J1),PS2(J1))
         ZFN2S = FN2(PS(J1),PS1(J1),PS2(J1))
         ZFN3S = FN3(PS(J1),PS1(J1),PS2(J1))
         ZFN4S = FN4(PS(J1),PS1(J1),PS2(J1))
C     
         ZFN1T = FN1(PT(J1),PT1(J1),PT2(J1))
         ZFN2T = FN2(PT(J1),PT1(J1),PT2(J1))
         ZFN3T = FN3(PT(J1),PT1(J1),PT2(J1))
         ZFN4T = FN4(PT(J1),PT1(J1),PT2(J1))
C
         ZDFN1S = DFN1(PS(J1),PS1(J1),PS2(J1))
         ZDFN2S = DFN2(PS(J1),PS1(J1),PS2(J1))
         ZDFN3S = DFN3(PS(J1),PS1(J1),PS2(J1))
         ZDFN4S = DFN4(PS(J1),PS1(J1),PS2(J1))
C     
         ZDFN1T = DFN1(PT(J1),PT1(J1),PT2(J1))
         ZDFN2T = DFN2(PT(J1),PT1(J1),PT2(J1))
         ZDFN3T = DFN3(PT(J1),PT1(J1),PT2(J1))
         ZDFN4T = DFN4(PT(J1),PT1(J1),PT2(J1))
C         
         PDFDS(J1, 1) = ZDFN1S * ZFN1T
         PDFDS(J1, 2) = ZDFN2S * ZFN1T
         PDFDS(J1, 3) = ZDFN1S * ZFN2T
         PDFDS(J1, 4) = ZDFN2S * ZFN2T
C
         PDFDS(J1, 5) = ZDFN3S * ZFN1T
         PDFDS(J1, 6) = ZDFN4S * ZFN1T
         PDFDS(J1, 7) = ZDFN3S * ZFN2T
         PDFDS(J1, 8) = ZDFN4S * ZFN2T
C
         PDFDS(J1, 9) = ZDFN1S * ZFN3T
         PDFDS(J1,10) = ZDFN2S * ZFN3T
         PDFDS(J1,11) = ZDFN1S * ZFN4T
         PDFDS(J1,12) = ZDFN2S * ZFN4T
C
         PDFDS(J1,13) = ZDFN3S * ZFN3T
         PDFDS(J1,14) = ZDFN4S * ZFN3T
         PDFDS(J1,15) = ZDFN3S * ZFN4T
         PDFDS(J1,16) = ZDFN4S * ZFN4T
C         
         PDFDT(J1, 1) = ZFN1S * ZDFN1T
         PDFDT(J1, 2) = ZFN2S * ZDFN1T
         PDFDT(J1, 3) = ZFN1S * ZDFN2T
         PDFDT(J1, 4) = ZFN2S * ZDFN2T
C
         PDFDT(J1, 5) = ZFN3S * ZDFN1T
         PDFDT(J1, 6) = ZFN4S * ZDFN1T
         PDFDT(J1, 7) = ZFN3S * ZDFN2T
         PDFDT(J1, 8) = ZFN4S * ZDFN2T
C
         PDFDT(J1, 9) = ZFN1S * ZDFN3T
         PDFDT(J1,10) = ZFN2S * ZDFN3T
         PDFDT(J1,11) = ZFN1S * ZDFN4T
         PDFDT(J1,12) = ZFN2S * ZDFN4T
C
         PDFDT(J1,13) = ZFN3S * ZDFN3T
         PDFDT(J1,14) = ZFN4S * ZDFN3T
         PDFDT(J1,15) = ZFN3S * ZDFN4T
         PDFDT(J1,16) = ZFN4S * ZDFN4T
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SY03
C*CALL PROCESS
         SUBROUTINE BASIS3(KN,KPN,PS1,PS2,PT1,PT2,PS,PT,PDFDS,PDFDT,
     &                     PDFDST,PD2FS2,PD2FT2)
C        ###########################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SY03 IS SIMILAR TO C2SY01. IT COMPUTES FIRST AND SECOND           *
*        DERIVATIVES OF 2-D CUBIC HERMITE BASIS FUNCTIONS             *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PS(KPN),        PS1(KPN),       PS2(KPN),       PT(KPN),
     R   PT1(KPN),       PT2(KPN),       PDFDS(KPN,16),  PDFDT(KPN,16),
     R   PDFDST(KPN,16), PD2FS2(KPN,16), PD2FT2(KPN,16)
C
         INCLUDE 'HERMIT.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         DO 1 J1=1,KN
         ZFN1S = FN1(PS(J1),PS1(J1),PS2(J1))
         ZFN2S = FN2(PS(J1),PS1(J1),PS2(J1))
         ZFN3S = FN3(PS(J1),PS1(J1),PS2(J1))
         ZFN4S = FN4(PS(J1),PS1(J1),PS2(J1))
C     
         ZFN1T = FN1(PT(J1),PT1(J1),PT2(J1))
         ZFN2T = FN2(PT(J1),PT1(J1),PT2(J1))
         ZFN3T = FN3(PT(J1),PT1(J1),PT2(J1))
         ZFN4T = FN4(PT(J1),PT1(J1),PT2(J1))
C
         ZDFN1S = DFN1(PS(J1),PS1(J1),PS2(J1))
         ZDFN2S = DFN2(PS(J1),PS1(J1),PS2(J1))
         ZDFN3S = DFN3(PS(J1),PS1(J1),PS2(J1))
         ZDFN4S = DFN4(PS(J1),PS1(J1),PS2(J1))
C     
         ZDFN1T = DFN1(PT(J1),PT1(J1),PT2(J1))
         ZDFN2T = DFN2(PT(J1),PT1(J1),PT2(J1))
         ZDFN3T = DFN3(PT(J1),PT1(J1),PT2(J1))
         ZDFN4T = DFN4(PT(J1),PT1(J1),PT2(J1))
C
         ZD2FN1S = D2FN1(PS(J1),PS1(J1),PS2(J1))
         ZD2FN2S = D2FN2(PS(J1),PS1(J1),PS2(J1))
         ZD2FN3S = D2FN3(PS(J1),PS1(J1),PS2(J1))
         ZD2FN4S = D2FN4(PS(J1),PS1(J1),PS2(J1))
C     
         ZD2FN1T = D2FN1(PT(J1),PT1(J1),PT2(J1))
         ZD2FN2T = D2FN2(PT(J1),PT1(J1),PT2(J1))
         ZD2FN3T = D2FN3(PT(J1),PT1(J1),PT2(J1))
         ZD2FN4T = D2FN4(PT(J1),PT1(J1),PT2(J1))
C         
         PDFDS(J1, 1) = ZDFN1S * ZFN1T
         PDFDS(J1, 2) = ZDFN2S * ZFN1T
         PDFDS(J1, 3) = ZDFN1S * ZFN2T
         PDFDS(J1, 4) = ZDFN2S * ZFN2T
C
         PDFDS(J1, 5) = ZDFN3S * ZFN1T
         PDFDS(J1, 6) = ZDFN4S * ZFN1T
         PDFDS(J1, 7) = ZDFN3S * ZFN2T
         PDFDS(J1, 8) = ZDFN4S * ZFN2T
C
         PDFDS(J1, 9) = ZDFN1S * ZFN3T
         PDFDS(J1,10) = ZDFN2S * ZFN3T
         PDFDS(J1,11) = ZDFN1S * ZFN4T
         PDFDS(J1,12) = ZDFN2S * ZFN4T
C
         PDFDS(J1,13) = ZDFN3S * ZFN3T
         PDFDS(J1,14) = ZDFN4S * ZFN3T
         PDFDS(J1,15) = ZDFN3S * ZFN4T
         PDFDS(J1,16) = ZDFN4S * ZFN4T
C         
         PDFDT(J1, 1) = ZFN1S * ZDFN1T
         PDFDT(J1, 2) = ZFN2S * ZDFN1T
         PDFDT(J1, 3) = ZFN1S * ZDFN2T
         PDFDT(J1, 4) = ZFN2S * ZDFN2T
C
         PDFDT(J1, 5) = ZFN3S * ZDFN1T
         PDFDT(J1, 6) = ZFN4S * ZDFN1T
         PDFDT(J1, 7) = ZFN3S * ZDFN2T
         PDFDT(J1, 8) = ZFN4S * ZDFN2T
C
         PDFDT(J1, 9) = ZFN1S * ZDFN3T
         PDFDT(J1,10) = ZFN2S * ZDFN3T
         PDFDT(J1,11) = ZFN1S * ZDFN4T
         PDFDT(J1,12) = ZFN2S * ZDFN4T
C
         PDFDT(J1,13) = ZFN3S * ZDFN3T
         PDFDT(J1,14) = ZFN4S * ZDFN3T
         PDFDT(J1,15) = ZFN3S * ZDFN4T
         PDFDT(J1,16) = ZFN4S * ZDFN4T
C         
         PDFDST(J1, 1) = ZDFN1S * ZDFN1T
         PDFDST(J1, 2) = ZDFN2S * ZDFN1T
         PDFDST(J1, 3) = ZDFN1S * ZDFN2T
         PDFDST(J1, 4) = ZDFN2S * ZDFN2T
C
         PDFDST(J1, 5) = ZDFN3S * ZDFN1T
         PDFDST(J1, 6) = ZDFN4S * ZDFN1T
         PDFDST(J1, 7) = ZDFN3S * ZDFN2T
         PDFDST(J1, 8) = ZDFN4S * ZDFN2T
C
         PDFDST(J1, 9) = ZDFN1S * ZDFN3T
         PDFDST(J1,10) = ZDFN2S * ZDFN3T
         PDFDST(J1,11) = ZDFN1S * ZDFN4T
         PDFDST(J1,12) = ZDFN2S * ZDFN4T
C
         PDFDST(J1,13) = ZDFN3S * ZDFN3T
         PDFDST(J1,14) = ZDFN4S * ZDFN3T
         PDFDST(J1,15) = ZDFN3S * ZDFN4T
         PDFDST(J1,16) = ZDFN4S * ZDFN4T
C
         PD2FS2(J1, 1) = ZD2FN1S * ZFN1T
         PD2FS2(J1, 2) = ZD2FN2S * ZFN1T
         PD2FS2(J1, 3) = ZD2FN1S * ZFN2T
         PD2FS2(J1, 4) = ZD2FN2S * ZFN2T
C
         PD2FS2(J1, 5) = ZD2FN3S * ZFN1T
         PD2FS2(J1, 6) = ZD2FN4S * ZFN1T
         PD2FS2(J1, 7) = ZD2FN3S * ZFN2T
         PD2FS2(J1, 8) = ZD2FN4S * ZFN2T
C
         PD2FS2(J1, 9) = ZD2FN1S * ZFN3T
         PD2FS2(J1,10) = ZD2FN2S * ZFN3T
         PD2FS2(J1,11) = ZD2FN1S * ZFN4T
         PD2FS2(J1,12) = ZD2FN2S * ZFN4T
C
         PD2FS2(J1,13) = ZD2FN3S * ZFN3T
         PD2FS2(J1,14) = ZD2FN4S * ZFN3T
         PD2FS2(J1,15) = ZD2FN3S * ZFN4T
         PD2FS2(J1,16) = ZD2FN4S * ZFN4T
C         
         PD2FT2(J1, 1) = ZFN1S * ZD2FN1T
         PD2FT2(J1, 2) = ZFN2S * ZD2FN1T
         PD2FT2(J1, 3) = ZFN1S * ZD2FN2T
         PD2FT2(J1, 4) = ZFN2S * ZD2FN2T
C
         PD2FT2(J1, 5) = ZFN3S * ZD2FN1T
         PD2FT2(J1, 6) = ZFN4S * ZD2FN1T
         PD2FT2(J1, 7) = ZFN3S * ZD2FN2T
         PD2FT2(J1, 8) = ZFN4S * ZD2FN2T
C
         PD2FT2(J1, 9) = ZFN1S * ZD2FN3T
         PD2FT2(J1,10) = ZFN2S * ZD2FN3T
         PD2FT2(J1,11) = ZFN1S * ZD2FN4T
         PD2FT2(J1,12) = ZFN2S * ZD2FN4T
C
         PD2FT2(J1,13) = ZFN3S * ZD2FN3T
         PD2FT2(J1,14) = ZFN4S * ZD2FN3T
         PD2FT2(J1,15) = ZFN3S * ZD2FN4T
         PD2FT2(J1,16) = ZFN4S * ZD2FN4T
C
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SY04
C*CALL PROCESS
         SUBROUTINE BASIS4(KN,KPN,PS1,PS2,PT1,PT2,PS,PT,PDFDS)
C        #####################################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
C
***********************************************************************
*                                                                     *
* C2SY04 IS SIMILAR TO C2SY01. IT COMPUTES FIRST DERIVATIVES OF 2-D   *
*        CUBIC HERMITE BASIS FUNCTIONS IN THE SIGMA DIRECTION         *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PS(KPN),    PS1(KPN),    PS2(KPN),       PT(KPN),
     R   PT1(KPN),   PT2(KPN),    PDFDS(KPN,16)
C
         INCLUDE 'HERMIT.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         DO 1 J1=1,KN
C     
         ZFN1T = FN1(PT(J1),PT1(J1),PT2(J1))
         ZFN2T = FN2(PT(J1),PT1(J1),PT2(J1))
         ZFN3T = FN3(PT(J1),PT1(J1),PT2(J1))
         ZFN4T = FN4(PT(J1),PT1(J1),PT2(J1))
C
         ZDFN1S = DFN1(PS(J1),PS1(J1),PS2(J1))
         ZDFN2S = DFN2(PS(J1),PS1(J1),PS2(J1))
         ZDFN3S = DFN3(PS(J1),PS1(J1),PS2(J1))
         ZDFN4S = DFN4(PS(J1),PS1(J1),PS2(J1))
C         
         PDFDS(J1, 1) = ZDFN1S * ZFN1T
         PDFDS(J1, 2) = ZDFN2S * ZFN1T
         PDFDS(J1, 3) = ZDFN1S * ZFN2T
         PDFDS(J1, 4) = ZDFN2S * ZFN2T
C
         PDFDS(J1, 5) = ZDFN3S * ZFN1T
         PDFDS(J1, 6) = ZDFN4S * ZFN1T
         PDFDS(J1, 7) = ZDFN3S * ZFN2T
         PDFDS(J1, 8) = ZDFN4S * ZFN2T
C
         PDFDS(J1, 9) = ZDFN1S * ZFN3T
         PDFDS(J1,10) = ZDFN2S * ZFN3T
         PDFDS(J1,11) = ZDFN1S * ZFN4T
         PDFDS(J1,12) = ZDFN2S * ZFN4T
C
         PDFDS(J1,13) = ZDFN3S * ZFN3T
         PDFDS(J1,14) = ZDFN4S * ZFN3T
         PDFDS(J1,15) = ZDFN3S * ZFN4T
         PDFDS(J1,16) = ZDFN4S * ZFN4T
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SY05
C*CALL PROCESS
         SUBROUTINE PSICEL(KS,KT,KN,KPN,PCEL,PSICL)
C        ##########################################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SY05 EXTRACT PSI AND ITS DERIVATIVES AT THE 4 NODES OF A CELL     *
*        FROM THE CPSICL ARRAY, WHICH CONTAINS THE COMPLETE BICUBIC   *
*        HERMITE SOLUTION                                             *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMNUM.inc'
C
         DIMENSION
     I   KS(KPN),         KT(KPN),
     R   PCEL(KPN,16),    PSICL(*)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         DO 1 J1=1,KN
C
         JS = KS(J1)
         JT = KT(J1)
         K  = (JS - 1) * NT + JT
         I4K = 4*K
         I4KPNT = 4*(K+NT)
         I4KMNT = 4*(K-NT)
C
         PCEL(J1,1) = PSICL(I4K-3)
         PCEL(J1,2) = PSICL(I4K-2)
         PCEL(J1,3) = PSICL(I4K-1)
         PCEL(J1,4) = PSICL(I4K  )
         PCEL(J1,5) = PSICL(I4KPNT-3)
         PCEL(J1,6) = PSICL(I4KPNT-2)
         PCEL(J1,7) = PSICL(I4KPNT-1)
         PCEL(J1,8) = PSICL(I4KPNT  )
C
         IF (JT .EQ. NT) THEN
C
            PCEL(J1, 9) = PSICL(I4KMNT+1)
            PCEL(J1,10) = PSICL(I4KMNT+2)
            PCEL(J1,11) = PSICL(I4KMNT+3)
            PCEL(J1,12) = PSICL(I4KMNT+4)
            PCEL(J1,13) = PSICL(I4K+1)
            PCEL(J1,14) = PSICL(I4K+2)
            PCEL(J1,15) = PSICL(I4K+3)
            PCEL(J1,16) = PSICL(I4K+4)
C
         ELSE
C
            PCEL(J1, 9) = PSICL(I4K+1)
            PCEL(J1,10) = PSICL(I4K+2)
            PCEL(J1,11) = PSICL(I4K+3)
            PCEL(J1,12) = PSICL(I4K+4)
            PCEL(J1,13) = PSICL(I4KPNT+1)
            PCEL(J1,14) = PSICL(I4KPNT+2)
            PCEL(J1,15) = PSICL(I4KPNT+3)
            PCEL(J1,16) = PSICL(I4KPNT+4)
C
         ENDIF
C
    1    CONTINUE
C
         RETURN
         END
C*DECK C2SY06
C*CALL PROCESS
         SUBROUTINE PSIBOX(KPSI1)
C       #########################
C
C                                        AUTHOR O. SAUTER, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SY06 EVALUATE PSI ON (R(I),Z(J)), I=1,NRBOX, J=1,NZBOX EQUIDISTANT*
*        MESH SUCH THAT R(1) = RBOXLFT , R(NRBOX) = RBOXLFT + RBOXLEN *
*        Z(1) = -ZBOXLEN/2. , Z(NZBOX) = ZBOXLEN/2.                   *
*                                                                     *
*        ASSUME THAT KPSI1 ISO-SURFACE HAVE BEEN CALCULATED BEFORE    *
*                                                                     *
*        USE SAME ALGORITHM AS IN ROUTINE EVLATE FOR INTERIOR POINTS  *
*        FOR OUTSIDE POINTS, USE CUBIC EXTRAPOLATION WITH PSI AND     *
*        DPSI/DSIGMA AT CSIG(NS) AND CSIG(NS1)                        *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMEQD.inc'
C
         DIMENSION
     R      ZCPSI(1,16),     ZDFDS(1,16),    ZDFDT(1,16),   
     R      ZF(1,16)
C
         INCLUDE 'CUCDCD.inc'
         INCLUDE 'CUCCCC.inc'
         INCLUDE 'QUAQQQ.inc'
C-----------------------------------------------------------------------
C
C        CHECK THAT PLASMA BOUNDARY IS INSIDE BOX
C        NOTE THAT BOX DIMENSIONS ARE IN MKSA
C
         IGMAX = NMGAUS * NT1
         IRMAX = ISMAX(IGMAX,RRISO(1,KPSI1),1)
         IRMIN = ISMIN(IGMAX,RRISO(1,KPSI1),1)
         IZMAX = ISMAX(IGMAX,RZISO(1,KPSI1),1)
         IZMIN = ISMIN(IGMAX,RZISO(1,KPSI1),1)
         ZRMIN = RRISO(IRMIN,KPSI1)
         ZRMAX = RRISO(IRMAX,KPSI1)
         ZZMIN = RZISO(IZMIN,KPSI1)
         ZZMAX = RZISO(IZMAX,KPSI1)
C
         ZRBOXLFT = RBOXLFT / R0EXP
         ZRBOXLEN = RBOXLEN / R0EXP
         ZZBOXLEN = ZBOXLEN / R0EXP
         IFIXBOUN = 0
         IF (ZRBOXLEN.LE.0.0 .OR. ZZBOXLEN.LE.0.0) IFIXBOUN = 1
         IF (ZRMIN.LE.ZRBOXLFT .OR. ZRMAX.GE.ZRBOXLFT+ZRBOXLEN .OR.
     +       ZZMIN.LE.-ZZBOXLEN/2. .OR. ZZMAX.GE.ZZBOXLEN/2.
     +       .OR. IFIXBOUN.EQ.1) THEN
            WRITE(6,'(/,A)') ' ***********************************'
            PRINT *,' BAD VALUES FOR RBOXLFT, RBOXLEN OR ZBOXLEN'
            PRINT *,' THEY HAVE BEEN CHANGED FROM:'
            PRINT *,' RBOXLFT = ',RBOXLFT
            PRINT *,' RBOXLEN = ',RBOXLEN
            PRINT *,' ZBOXLEN = ',ZBOXLEN
            PRINT *,' TO:'
            PRINT *,' '
C
            ZSHIFTR = 1.1 * (ZRMAX-ZRMIN)/DFLOAT(NRBOX)
            ZSHIFTZ = 1.1 * (ZZMAX-ZZMIN)/DFLOAT(NZBOX)
            IF (ZRMIN.LE.ZRBOXLFT .OR. IFIXBOUN.EQ.1)
     +        ZRBOXLFT = ZRMIN - ZSHIFTR
C
            IF (ZRMAX .GE. ZRBOXLFT+ZRBOXLEN .OR. IFIXBOUN.EQ.1)
     +        ZRBOXLEN = ZRMAX - ZRMIN + 2.*ZSHIFTR
C     IF SHIFT ZAXIS
            ZZZLEN = 2.*MAX(ABS(ZZMAX-RZMAG),ABS(ZZMIN-RZMAG))
c%OS            ZZZLEN = 2.*MAX(ABS(ZZMAX),ABS(ZZMIN))
            IF (ZZZLEN .GE. ZZBOXLEN .OR. IFIXBOUN.EQ.1)
     +          ZZBOXLEN = ZZZLEN + 1.2*ZSHIFTZ
C
            RBOXLFT = ZRBOXLFT * R0EXP
            RBOXLEN = ZRBOXLEN * R0EXP
            ZBOXLEN = ZZBOXLEN * R0EXP
C
            PRINT *,' RBOXLFT = ',RBOXLFT
            PRINT *,' RBOXLEN = ',RBOXLEN
            PRINT *,' ZBOXLEN = ',ZBOXLEN
            PRINT *,' PLASMA EDGES:'
            PRINT *,' ZRMIN= ',ZRMIN * R0EXP
            PRINT *,' ZRMAX= ',ZRMAX * R0EXP
            PRINT *,' ZZMIN= ',ZZMIN * R0EXP
            PRINT *,' ZZMAX= ',ZZMAX * R0EXP
C
            PRINT *,' '
         ENDIF

C
C-----------------------------------------------------------------------
CL    1.2 SHIFT Z-AXIS SUCH THAT RZMAG = 0.0 (I.E. RZMAG IS CENTER OF Z-MESH)
C
         ZBOTTOM = RZMAG - 0.5*ZZBOXLEN
c%OS         ZBOTTOM = - 0.5*ZZBOXLEN
C
C-----------------------------------------------------------------------
CL       2. COMPUTE PSI VALUE. NOTE THAT CPSICL WAS SHIFTED BY CPSRF
C        SO SHIFT IT BACK
C
         ZDR = ZRBOXLEN / DFLOAT(NRBOX-1)
         ZDZ = ZZBOXLEN / DFLOAT(NZBOX-1)
         print *,' r0, rz0, rmag, rzmag= ',r0, rz0, rmag, rzmag
         print *,' bps(1),bps(12) = ',bps(1),bps(12),'  r0exp= ',r0exp
         bps(1) = r0
         bps(12) = rz0
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
         DO J=1,NZBOX
            ZZ = ZBOTTOM + DFLOAT(J-1)*ZDZ
C
            DO I=1,NRBOX
               ZR = ZRBOXLFT + DFLOAT(I-1)*ZDR
C
               ZRHO = SQRT((ZR - BPS(1))**2 + (ZZ - BPS(12))**2)
               ZTET = ATAN2(ZZ - BPS(12),ZR - BPS(1))
               IF (ZTET .LT. CT(1)) ZTET = ZTET + 2. * CPI
C
               CALL BOUND(1,ZTET,ZBND)
               IT = ISRCHFGE(NT1,CT,1,ZTET) - 1
               IF (IT .LT. 1)  IT = 1
               IF (IT .GT. NT) IT = NT
               ZT1 = CT(IT)
               ZT2 = CT(IT+1)
C
               ZSIG = ZRHO / ZBND
C
               IF (ZSIG .LE. 1.0) THEN
C        INSIDE POINT
                  IS = ISRCHFGE(NS1,CSIG,1,ZSIG)  - 1
                  IF (IS .LT. 1)  IS = 1
                  IF (IS .GT. NS) IS = NS
C
                  ZS1 = CSIG(IS)
                  ZS2 = CSIG(IS+1)
C
                  CALL PSICEL(IS,IT,1,1,ZCPSI,CPSICL)
                  CALL BASIS1(1,1,ZS1,ZS2,ZT1,ZT2,ZSIG,ZTET,ZF)
C
                  EQDSPSI(I,J) = SDOT(16,ZF,1,ZCPSI,1) - CPSRF
C
               ELSE IF (ZSIG .GT. 1.0) THEN
                  ZSIGM1 = CSIG(NS)
                  ZSIGM3 = ZSIGM1 + 0.33*(1.-ZSIGM1)
                  ZSIGM2 = ZSIGM1 + 0.67*(1.-ZSIGM1)
C
C        EVALUATE PSI AT SIGMA = ZSIGM3
c%OS                  IS = ISRCHFGT(NS1,CSIG,1,ZSIGM3)  - 1
                  IS = NS
                  ZS1 = CSIG(IS)
                  ZS2 = CSIG(IS+1)
                  CALL PSICEL(IS,IT,1,1,ZCPSI,CPSICL)
                  CALL BASIS1(1,1,ZS1,ZS2,ZT1,ZT2,ZSIGM3,ZTET,ZF)
                  ZPSIM3 = SDOT(16,ZF,1,ZCPSI,1)
C
C        EVALUATE PSI AT SIGMA = ZSIGM2
c%OS                  IS = ISRCHFGT(NS1,CSIG,1,ZSIGM2)  - 1
c%OS                  ZS1 = CSIG(IS)
c%OS                  ZS2 = CSIG(IS+1)
                  CALL PSICEL(IS,IT,1,1,ZCPSI,CPSICL)
                  CALL BASIS1(1,1,ZS1,ZS2,ZT1,ZT2,ZSIGM2,ZTET,ZF)
                  ZPSIM2 = SDOT(16,ZF,1,ZCPSI,1)
C
C        EVALUATE PSI AT SIGMA = ZSIGM1
c%OS                  IS = ISRCHFGT(NS1,CSIG,1,ZSIGM1)  - 1
c%OS                  ZS1 = CSIG(IS)
c%OS                  ZS2 = CSIG(IS+1)
                  CALL PSICEL(IS,IT,1,1,ZCPSI,CPSICL)
                  CALL BASIS1(1,1,ZS1,ZS2,ZT1,ZT2,ZSIGM1,ZTET,ZF)
                  ZPSIM1 = SDOT(16,ZF,1,ZCPSI,1)
C        EVALUATE DPSI/DSIGMA AT SIGMA = ZSIGM1
                  CALL BASIS2(1,1,ZS1,ZS2,ZT1,ZT2,ZSIGM1,ZTET,ZDFDS,
     +                ZDFDT)
                  ZDPDSM1 = SDOT(16,ZDFDS,1,ZCPSI,1)
C
C        EVALUATE PSI AT SIGMA = 1.0
                  ZSIG1 = 1.0
                  IS = NS
                  ZS1 = CSIG(IS)
                  ZS2 = CSIG(IS+1)
                  CALL PSICEL(IS,IT,1,1,ZCPSI,CPSICL)
                  CALL BASIS1(1,1,ZS1,ZS2,ZT1,ZT2,ZSIG1,ZTET,ZF)
                  ZPSI1 = SDOT(16,ZF,1,ZCPSI,1)
C
C        EVALUATE DPSI/DSIGMA AT SIGMA = 1.0
                  CALL BASIS2(1,1,ZS1,ZS2,ZT1,ZT2,ZSIG1,ZTET,ZDFDS,
     +                ZDFDT)
                  ZDPDS1 = SDOT(16,ZDFDS,1,ZCPSI,1)
c%OSC
c%OSC        EVALUATE PSI AT ZSIG WITH QUADRATIC
c%OS                  EQDSPSI(I,J) = FQQQ0(ZPSIM3,ZPSIM2,ZPSI1,
c%OS     +              ZSIGM3,ZSIGM2,ZSIG1,ZSIG) - CPSRF
c%OSC
c%OSC        EVALUATE PSI AT ZSIG WITH CUBIC USING DERIVATIVES
c%OS                  EQDSPSI(I,J) = FCDCD0(ZSIGM1,ZPSIM1,ZDPDSM1,
c%OS     +                ZSIG1,ZPSI1,ZDPDS1,ZSIG) - CPSRF
C
C        EVALUATE PSI AT ZSIG WITH CUBIC
                  EQDSPSI(I,J) = FCCCC0(ZPSIM1,ZPSIM3,ZPSIM2,ZPSI1,
     +              ZSIGM1,ZSIGM3,ZSIGM2,ZSIG1,ZSIG) - CPSRF
C
C     KEEP PSI VALUE POSITIVE (>PSILIM)
                  IF (EQDSPSI(I,J) .LE. 0.0) EQDSPSI(I,J) = 1.E-04
C
               ENDIF
            END DO
         END DO
         RETURN
         END
C*DECK C3SA01
C*CALL PROCESS
         SUBROUTINE OUTPUT(K)
C        ####################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SA01  INPUT / OUTPUT                                              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMLAB.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMPLO.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVEV.inc'
         INCLUDE 'COMEQD.inc'
         INCLUDE 'COMDAT.inc'
C
C
         DIMENSION
     R   ZCPSI(NSNT),   ZDPDS(NSNT),   ZDPDT(NSNT),   ZD2PST(NSNT),
     R   ZJDIFF(NPISO)
C-----------------------------------------------------------------------
         INCLUDE 'QUAQQQ.inc'
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         GOTO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,
     ,         160,170,180) K
C
C     PRINT OUT NAMELIST
C
   10    CONTINUE
C
C        CALL WNLLONG(132)
C
         WRITE(6,1001)
         WRITE (6,EQDATA)
C
         RETURN
C
C     READ IN NAMELIST
C
   20    CONTINUE
C
         READ (5,EQDATA)
         IF( (NRFP.EQ.1)) THEN
            NWBPS=2
            NBPS=NCHI
            DTHIN=ACOS(-1.)*2./NCHI*2
            DO IK=1,NBPS
               THIN=DTHIN*(IK*2-1)
               RRBPS(IK,1)=1.+ASPCT*COS(THIN)
               RZBPS(IK,1)=ASPCT*SIN(THIN)
               RRBPS(IK,2)=1.+ASPCT*COS(THIN)*REXT
               RZBPS(IK,2)=ASPCT*SIN(THIN)*REXT
               TETBPS(IK,1)=THIN
               TETBPS(IK,2)=THIN
            END DO
            IF( AT(4).EQ.1.) AT(4)= 0.99999999
         END IF

C        SAVE NAMELIST VARIABLES R0EXP AND B0EXP
C        BECAUSE THESE ARE MODIFIED BY THE CODE
         R0EXP0 = R0EXP
         B0EXP0 = B0EXP
C
         RETURN
C
C     FINAL OUTPUT
C
   30    CONTINUE
C
         IF (NPRPSI .EQ. 1) THEN
C
            CALL SCOPY(NSTMAX,CPSICL(1),4,ZCPSI,1)
            CALL SCOPY(NSTMAX,CPSICL(2),4,ZDPDS,1)
            CALL SCOPY(NSTMAX,CPSICL(3),4,ZDPDT,1)
            CALL SCOPY(NSTMAX,CPSICL(4),4,ZD2PST,1)
C
         ENDIF
C
         WRITE(6,1100)
C
         IF (NSURF .NE. 1) CALL RVAR('SCALE                   ',SCALE)
C
         CALL RVAR('PSI - MIN               ',SPSIM)
         CALL RVAR('R OF MAGAXE             ',RRAXIS)
         CALL RVAR('Z OF MAGAXE             ',RZAXIS)
C
         IF (NOPT .EQ. 0) THEN
C
            CALL RVAR('CHECK OF SOLUTION       ',SCHECK)
            CALL RVAR('POLOIDAL MAGNETIC ENERGY',WMAGP)
C
         ENDIF
C
         IF (NTEST .EQ. 0 .AND. NPRPSI .EQ. 1) THEN
C
            CALL RARRAY('PSI NUMERIQUE',ZCPSI,NSTMAX)
            CALL RARRAY('DPSI/DS',ZDPDS,NSTMAX)
            CALL RARRAY('DPSI/DT',ZDPDT,NSTMAX)
            CALL RARRAY('D2P/DST',ZD2PST,NSTMAX)
C
         ELSE IF (NTEST. EQ. 1) THEN
C
           IF (NPRPSI .EQ. 1) THEN
C
              CALL RARRAY('PSI THEORIQUE',CPSI1T,NSTMAX)
              CALL RARRAY('PSI NUMERIQUE',ZCPSI,NSTMAX)
              CALL RARRAY('ABSOLUT DIFFERENCES OF PSI ON NODES',
     ,                                             DIFFP,NSTMAX)
C
           ENDIF
C
           CALL RVAR('RESIDU OF PSI ON NODES        ',RESPSI)
C
           IF (NPRPSI .EQ. 1) THEN
C
              CALL RARRAY('DPSI/DS THEORIQUE',DPDSTH,NSTMAX)
              CALL RARRAY('DPSI/DS NUMERIQUE',DPDSNU,NSTMAX)
              CALL RARRAY('ABSOLUT DIFFERENCES OF DPSI / DS ON NODES',
     ,                                                   DIFFDS,NSTMAX)
C
           ENDIF
C
           CALL RVAR('RESIDU OF DPSI / DS ON NODES  ',RESDPS)
C
           IF (NPRPSI .EQ. 1) THEN
C
              CALL RARRAY('DPSI/DT THEORIQUE',DPDTTH,NSTMAX)
              CALL RARRAY('DPSI/DT NUMERIQUE',DPDTNU,NSTMAX)
              CALL RARRAY('ABSOLUT DIFFERENCES OF DPSI / DT ON NODES',
     ,                                                 DIFFDT,NSTMAX)
C
           ENDIF
C
           CALL RVAR('RESIDU OF DPSI / DT ON NODES  ',RESDPT)
C
           IF (NPRPSI .EQ. 1) THEN
C
              CALL RARRAY('D2PSI/DST THEORIQUE',D2PSTT,NSTMAX)
              CALL RARRAY('D2PSI/DST NUMERIQUE',D2PSTN,NSTMAX)
              CALL RARRAY('ABSOLUT DIFFERENCES OF D2PSI / DST ON NODES',
     ,                                                 DIFFST,NSTMAX)
C
           ENDIF
C
           CALL RVAR('RESIDU OF D2PSI / DST ON NODES ',RESDST)
C
         ENDIF
C
         RETURN
C
   40    CONTINUE
C
         IF (NSTTP .LT. 2) THEN
C
            IF (NSURF .NE. 1) THEN
C
               WRITE(6,1201) NS,NT,NPSI,NCHI,NINSCA
C
            ELSE
C
               WRITE(6,1202) NS,NT,NPSI,NCHI
C
            ENDIF
C
         ELSE
C
            WRITE(6,1203) NS,NT,NISO,NPSI,NCHI,NINSCA,NINMAP
C
         ENDIF
C
         RETURN
C
   50    CONTINUE
C
         WRITE(6,1250) SPSIM,RMAG,RZMAG
C
         RETURN
C
   60    CONTINUE
C
         WRITE(6,1301) RESIDU,CEPS
C
         RETURN
C
C     PRINT OUT MESH
C
   70    CONTINUE
C
         CALL RARRAY('SIGMA - MESH',CSIG,NS1)
         CALL RARRAY('THETA - MESH',CT,NT1)
         CALL RARRAY('RHOS - MESH',RHOS,NT1)
C
         RETURN
C
   80    CONTINUE
C
         WRITE(6,1302) RESMAP
C
         RETURN
C
   90    CONTINUE
C
         WRITE (6,1450)
C
         RETURN
C
  100    CONTINUE
C
C     SHIFT OF THE MAGNETIC AXIS
C
         WRITE(6,2050) R0,RZ0
         CALL BLINES(1)
C
         RETURN
C
  110    CONTINUE
C
C     FLUX SURFACES AVERAGED QUANTITIES AND BETAS
C
         WRITE(6,2300)
         WRITE(6,2301)
         ZMU0 = 1.256
         ZBPERC = 100. * BETA
         ZBSPER = 100. * BETAS
         ZBXPER = 100. * BETAX
         ZINORM = RINOR/ZMU0
         ZIBSNO = RIBSNOR/ZMU0
         ZGM    = ZBPERC/ZINORM
         ZGMSTA = ZBSPER/ZINORM
         ZGMX   = ZBXPER/ZINORM
         ZBPOL1 = BETAP*CONVF
         ZLI1   = RINDUC(NPSI1)*CONVF
         ZBSF   = RITBS/RITOT
         ZBSFC  = RITBSC/RITOT
         ZCBS1  = 0.
         ZCBS2  = 0.
         IF (BETA.GT.0.) ZCBS1 = ZIBSNO*ZINORM/(ZBXPER*SQRT(ASPCT))
         IF (BETA.GT.0.) ZCBS2 = ZBSF/(SQRT(ASPCT)*ZBPOL1)
         ZBBS  = ZBXPER * ZBSF
C
         CALL RVAR('AVERAGED PRESSURE      ',CPBAR)
         CALL RVAR('TOTAL CURRENT          ',RITOT)
         CALL RVAR('NORMALIZED CURRENT     ',RINOR)
         CALL RVAR('IN    (MA,T,M)         ',ZINORM)
         CALL RVAR('PRESSURE PEAKING FACTOR',CPPF)
         CALL RVAR('CONVERSION FACTOR      ',CONVF)
         CALL RVAR('POLOIDAL BETA          ',BETAP)
         CALL RVAR('POLOIDAL BETA (GA)     ',ZBPOL1)
         CALL RVAR('LI                     ',RINDUC(NPSI1))
         CALL RVAR('LI (GA)                ',ZLI1  )
         CALL RVAR('BETA  [%]              ',ZBPERC)
         CALL RVAR('BETA* [%]              ',ZBSPER)
         CALL RVAR('BETAX [%]              ',ZBXPER)
         CALL RVAR('G     (MA,T,M)         ',ZGM   )
         CALL RVAR('G*    (MA,T,M)         ',ZGMSTA)
         CALL RVAR('GEXP  (MA,T,M)         ',ZGMX  )
         CALL RVAR('F0=IB.S./ITOT (NUE*=0) ',ZBSF  )
         CALL RVAR('IB.S./ITOT (NUE*.NE.0) ',ZBSFC )
         CALL RVAR('CBS1=IB.S.(0)/(G*SQRT(E)) ',ZCBS1 )
         CALL RVAR('CBS2=F0/(BP(1)*SQR(E)) ',ZCBS2 )
         CALL RVAR('BBS(NUE*=0)            ',ZBBS  )
C
         IF (NRFP .EQ. 1) THEN
C
            CALL RVAR('F                      ',RFPF )
            CALL RVAR('THETA                  ',RFPT )
C
         ENDIF
C
         WRITE(6,2302)
C
         CALL RVAR('T(PSI)ON AXIS         ',T0   )
         CALL RVAR('TT-PRIME ON AXIS      ',DTTP0)
         CALL RVAR('Q(PSI) ON AXIS        ',Q0   )
         CALL RVAR('DQ/DPSI ON AXIS       ',DQDP0)
         CALL RVAR('CP(PSI) ON AXIS       ',CPSI0)
         CALL RVAR('DCP/DPSI ON AXIS      ',CPDP0)
         CALL RVAR('PRESSURE ON AXIS      ',CP0  )
         CALL RVAR('DP/DPSI ON AXIS       ',DPDP0)
         CALL RVAR('I - PRIME ON AXIS     ',RIPR0)
         CALL RVAR('<J.B>/<T/R**2> ON AXIS',RJDTB0)
C
         IF (NPROFZ .EQ. 1) THEN
C
            CALL RVAR('DENSITY ON AXIS       ',DENS0)
            CALL RVAR('TEMPERATURE ON AXIS   ',TEMP0)
C
         ENDIF
C
         WRITE(6,2350)
C
         JQMIN = 1
         QMIN = QPSI(1)
         ZS95 = SQRT(0.95)
C
         DO 111 J111=1,NPSI1
C
c%OS         ZJDIFF(J111) = RIPR(J111) - RJBSH(J111)
         ZJDIFF(J111) = RIPR(J111) - RJBSOS(J111,1)
C
         IF ( QPSI(J111) .LT. QMIN ) THEN
           JQMIN = J111
           QMIN = QPSI(JQMIN)
         ENDIF
         IF ( CSM(J111) .LT. ZS95 ) THEN
           JQ95 = J111
         ENDIF
C
 111     CONTINUE
         ZSQMIN=CSM(JQMIN)
C     QUADRATIC INTERPOLATION FOR QMIN
         IF (JQMIN.EQ.1 .OR. JQMIN.EQ.NPSI1) GO TO 1110
         ZS1 = CSM(JQMIN-1)
         ZF1 = QPSI(JQMIN-1)
         ZS2 = CSM(JQMIN)
         ZF2 = QPSI(JQMIN)
         ZS3 = CSM(JQMIN+1)
         ZF3 = QPSI(JQMIN+1)
         ZSQMIN = - FB1(ZF1,ZF2,ZF3,ZS1,ZS2,ZS3)
     +     / FB2(ZF1,ZF2,ZF3,ZS1,ZS2,ZS3) * 0.5
         QMIN = FQQQ0(ZF1,ZF2,ZF3,ZS1,ZS2,ZS3,ZSQMIN)
 1110    CONTINUE
C     QUADRATIC INTERPOLATION FOR Q95
         Q95 = FQQQ0(QPSI(JQ95-1),QPSI(JQ95),QPSI(JQ95+1),
     +     CSM(JQ95-1),CSM(JQ95),CSM(JQ95+1),ZS95)
C
         CALL RVAR('MINIMUM Q VALUE',QMIN)
         CALL RVAR('S VALUE OF QMIN',ZSQMIN)
         CALL RVAR('Q AT 95% FLUX SURFACE',Q95)
         CALL RVAR('Q AT 100% FLUX',QPSI(NPSI1))
         CALL RARRAY('CS  - MESH',CS,NPSI1)
         CALL RARRAY('BETA-POLOIDAL(CS)',BETAB,NPSI1)
         CALL RARRAY('CSM - MESH',CSM,NPSI1)
         CALL RARRAY('PSI(CSM)',PSIISO,NPSI1)
         CALL RARRAY('T (CSM)',TMF,NPSI1)
         CALL RARRAY('T*DT/DPSI(CSM)',TTP,NPSI1)
         CALL RARRAY('P (CSM)',CPR,NPSI1)
         CALL RARRAY('DP/DPSI(CSM)',CPPR,NPSI1)
         CALL RARRAY('Q (CSM)',QPSI,NPSI1)
         CALL RARRAY('DQ/DPSI(CSM)',CDQ,NPSI1)
C        CALL RARRAY('CP(CSM)',CP,NPSI1)
C        CALL RARRAY('D(CP)/DPSI(CSM)',CPDP,NPSI1)
         CALL RARRAY('SHEAR(CS)',CDRQ,NPSI1)
         CALL RARRAY('I-STAR(CSM)',RIPR,NPSI1)
         CALL RARRAY('I//=<J . B> / <T / R**2>',RJDOTB,NPSI1)
         CALL RARRAY('RFCIRC',RFCIRC,NPSI1)
         CALL RARRAY('RB2AV',RB2AV,NPSI1)
C         CALL RARRAY('BS-CURRENT (BOOZER)',RJBSR,NPSI1)
         CALL RARRAY('BS-CURRENT (HIRSHMAN)',RJBSH,NPSI1)
         CALL RARRAY('BS-CURRENT (SAUTER,ZERO COLL.)',RJBSOS(1,1),NPSI1)
         CALL RARRAY('BS-CURRENT (SAUTER)',RJBSOS(1,2),NPSI1)
         CALL RARRAY('NUESTAR (WITHOUT ZEFF)',RNUSTAR,NPSI1)
         CALL RARRAY('J-PARALLEL (<j.B>Eq.43)',RJPAR,NPSI1)
         CALL RARRAY('J-PAR - J-BS(0)',ZJDIFF,NPSI1)
         CALL RARRAY('ELLIPTICITY(CSM)',RELL,NPSI1)
         CALL RARRAY('D(ELL.) / DR (CS)',RDEDR,NPSI1)
         CALL RARRAY('MERCIER BY L.A.TH. WITH E ONLY (CS)',RDI,NPSI1)
         CALL RARRAY('MERCIER SHAFRANOV YOURCHENKO (CS)',RSY,NPSI1)
         CALL RARRAY('ASPECT RATIO(CSM)',ARATIO,NPSI1)
         CALL RARRAY('VSURF=VOLUME/2PI',VSURF,NPSI1)
C
C------- SECTION ADDED TO MATCH WITH GA'S EQUILIBRIUM CODE ---------

         NGA = 3
         OPEN(NGA,FILE='NGA',FORM='FORMATTED')
         REWIND(NGA)
         WRITE(NGA,*) NPSI1
         WRITE(NGA,*) SPSIM
         CALL OARRAY(NGA,'CSM - MESH',CSM,NPSI1)
         CALL OARRAY(NGA,'P (CSM)',CPR,NPSI1)
         CALL OARRAY(NGA,'DP/DPSI(CSM)',CPPR,NPSI1)
         CALL OARRAY(NGA,'Q (CSM)',QPSI,NPSI1)
         CLOSE(NGA)
C
C-------------------------------------------------------------------
C
         IF (NRFP .EQ. 1) THEN
C
            CALL RARRAY('SURFACE AVERAGED POLOIDAL MAGNETIC FIELD',
     ,                  RFPBP,NPSI1)
C
         ENDIF
C
         IF (NPROFZ .EQ. 1) THEN 
C
            CALL RARRAY('PLASMA DENSITY',DENSTY,NPSI1)
            CALL RARRAY('PLASMA TEMPERATURE',TEMPER,NPSI1)
C
         ENDIF
C
         CALL RARRAY('RHO(CS)',RSURF,NPSI1)
C
         RETURN
C
C   COMMON'S SIZES USED IN CLEAR
C
  120    CONTINUE
C
         WRITE(6,2400)
C
         IBLA  = NBLA
         IBAL  = NBAL1 + NBAL2
         IBND  = NBND1 + NBND2
         ICON  = NICON
         IERA  = NERA
         IESH  = NESH1 + NESH2
         IETA  = NETA1 + 2 * (NETA2 + NETA3 + NETA4)
         IINT  = NINT1 + NINT2
         IIOD  = NIOD
         IISO  = NIS2
         IMAP  = NMAP2
         INUM  = NNUM
         IPHY  = NPHY1 + NPHY2
         IPLO  = NPLO1 + NPLO2
         ISOL  = NSOL1 + NSOL2
         ISUR  = NISUR
         IVAC  = NCVAC*2
         IVEV  = NVEV
         IEQD  = NEQD1 + NEQD2
         ITOT = IBLA + IBAL + IBND + ICON + IERA + IESH + IETA + IINT +
     +     IIOD + IISO + IMAP + INUM + IPHY + IPLO + ISOL + ISUR
     +     + IVAC + IVEV + IEQD
C
         CALL IVAR('COMBAL',IBAL)
         CALL IVAR('COMBLA',IBLA)
         CALL IVAR('COMBND',IBND)
         CALL IVAR('COMCON',ICON)
         CALL IVAR('COMERA',IERA)
         CALL IVAR('COMESH',IESH)
         CALL IVAR('COMETA',IETA)
         CALL IVAR('COMINT',IINT)
         CALL IVAR('COMIOD',IIOD)
         CALL IVAR('COMISO',IISO)
         CALL IVAR('COMMAP',IMAP)
         CALL IVAR('COMNUM',INUM)
         CALL IVAR('COMPHY',IPHY)
         CALL IVAR('COMPLO',IPLO)
         CALL IVAR('COMSOL',ISOL)
         CALL IVAR('COMSUR',ISUR)
         CALL IVAR('COMVAC',NCVAC)
         CALL IVAR('COMVEV',IVEV)
         CALL IVAR('TOTAL ',ITOT)
C
         RETURN
C
  130    CONTINUE
C
         WRITE(6,2500)
C
         CALL RARRAY('CHI - VALUES',CHI,NCHI1)
         CALL RARRAY('CHIM - VALUES',CHIM,NCHI)
C
         RETURN
C
  140    CONTINUE
C
         WRITE(6,1500)
C
         RETURN
C
  150    CONTINUE
C
         CALL RARRAY('MERCIER',SMERCI,NPSI1)
         CALL RARRAY('RESISTIVE INTERCHANGE',SMERCR,NPSI1)
         CALL RARRAY('H OF GLASSER, GREENE & JOHNSON',HMERCR,NPSI1)
         CALL IVAR('NTURN',NTURN)
         CALL IARRAY('NCBAL',NCBAL,NPSI1)
         CALL RARRAY('CHI0 VALUES',CHI0,NBLC0)
C
         IF (NBLC0 .LE. 1) RETURN
         IF (.TRUE.) RETURN
C
         DO 151 J151=1,NBLC0
C
         IBL = ISSUM(NPSI1,NCBLNS(1,J151),1)
C
         IF (IBL .NE. 0) THEN
C
            CALL RVAR('CHI0',CHI0(J151))
            CALL IARRAY('NCBAL(CHI0)',NCBLNS(1,J151),NPSI1)
C
         ENDIF
C
  151    CONTINUE
C
         RETURN
C
  160    CONTINUE
C
         CALL PRIQQU
C
         RETURN
C
  170    CONTINUE
C
         WRITE(6,2301)
C
         CALL RVAR('POLOIDAL BETA      ',BETAP)
         CALL RVAR('TOROIDAL BETA      ',BETA )
         CALL RVAR('BETA STAR          ',BETAS)
         CALL RVAR('TOROIDAL BETA EXP. ',BETAX)
C
         WRITE(6,2302)
C
         CALL RVAR('PRESSURE ON AXIS   ',CP0  )
         CALL RVAR('DP/DPSI ON AXIS    ',DPDP0)
         CALL RVAR('IB.S./ITOT (NUE*=0)',RITBS/RITOT)
         CALL RVAR('IB.S./ITOT (NUE*.NE.0)',RITBSC/RITOT)
C
         CALL RARRAY('DP/DPSI(CSM)',CPPR,NPPR+1)
         CALL RARRAY('Q(CSM)',QPSI,NPPR+1)
         CALL RARRAY('DQ/DPSI(CSM)',CDQ,NPPR+1)
C
         IF (NBLOPT .NE. 0) THEN
C
            CALL IARRAY('N2BAL',N2BAL,NPPR+1)
            CALL RARRAY('MERCIER',SMERCI,NPPR+1)
            CALL IARRAY('NCBAL',NCBAL,NPPR+1)
            CALL RARRAY('CHI0 VALUES',CHI0,NBLC0)
C
cab
      IF (.TRUE.) RETURN
cab
            IF (NBLC0 .LE. 1) RETURN
C
            DO 171 J171=1,NBLC0
C
            IBL = ISSUM(NPPR+1,NCBLNS(1,J171),1)
C
            IF (IBL .NE. 0) THEN
C
               CALL RVAR('CHI0',CHI0(J171))
               CALL IARRAY('NCBAL(CHI0)',NCBLNS(1,J171),NPPR+1)
C
            ENDIF
C
  171       CONTINUE
C
         ENDIF
C
         IF (NBSOPT .NE. 0) THEN
C
            DO 172 J172=1,NPPR+1
C
c%OS            ZJDIFF(J172) = RIPR(J172) - RJBSH(J172)
            ZJDIFF(J172) = RIPR(J172) - RJBSOS(J172,1)
C
 172        CONTINUE
C
            CALL RARRAY('RFCIRC',RFCIRC,NPPR+1)
            CALL RARRAY('RB2AV',RB2AV,NPPR+1)
            CALL RARRAY('BS-CURRENT (BOOZER)',RJBSR,NPPR+1)
            CALL RARRAY('BS-CURRENT (HIRSHMAN)',RJBSH,NPPR+1)
            CALL RARRAY('J-PARALLEL (<j.B>Eq.43)',RJPAR,NPPR+1)
            CALL RARRAY('J-PAR - J-BS(0)',ZJDIFF,NPPR+1)
C
         ENDIF
C
         RETURN
C
  180    CONTINUE
C
         WRITE(6,2501) 100. * CPRESS
C
         RETURN         
C
 1001    FORMAT(///,1X,'******************',
     +          //,1X,'NAMELIST VARIABLES',
     +          //,1X,'******************',//)
 1100    FORMAT(///,1X,'******************************',
     ,           //,1X,'FINAL OUTPUT : PSI SOLUTION   ',
     ,           //,1X,'VALUES ON (SIGMA, THETA) NODES',
     ,           //,1X,'******************************')
 1201    FORMAT(//,1X,' NS = ',I4,4X,' NT = ',I4,4X,
     ,                ' NPSI = ',I4,4X,' NCHI = ',I4,4X,
     ,                ' NINSCA = ',I4)
 1202    FORMAT(//,1X,' NS = ',I4,4X,' NT = ',I4,4X,
     ,                ' NPSI = ',I4,4X,' NCHI = ',I4)
 1203    FORMAT(//,1X,' NS = ',I4,4X,' NT = ',I4,4X,' NISO = ',I4,4X,
     ,                ' NPSI = ',I4,4X,' NCHI = ',I4,4X,
     ,                ' NINSCA = ',I4,4X,' NINMAP = ',I4)
 1250    FORMAT(6X,'PSIMIN = ',1PE20.12,3X,'RMAG =   ',1PE20.12,
     ,          3X,'ZMAG =   ',1PE20.12)
 1300    FORMAT(73X,'MITER = ',I3,3X,'RESIDU = ',F12.10,3X,
     +              'EPSLON = ',F12.10)
 1301    FORMAT(73X,'RESIDU = ',E12.6,3X,'EPSLON = ',E12.6,/)
 1302    FORMAT(/,3X,'RESIDU OF ITERATION OVER MAPPING = ',E20.10,/)
 1400    FORMAT(/,1X,'CHECK THE FIRST POINT ON MAGNETIC FLUX SURFACES')
 1450    FORMAT(/,73X,'  *****   ITERATION OVER NON LINEARITY NOT ',
     +          'CONVERGED')
 1500    FORMAT(/,73X,'  *****   ITERATION OVER CURRENT PROFILE NOT ',
     +          'CONVERGED')
 2050    FORMAT(///,1X,'*****************************************',
     ,           //,1X,'POSITION OF THE CALCULATION MESH CENTER :',
     ,              1X,'R0 = ',1PE16.8,3X,'Z0 = ',1PE16.8,
     ,           //,1X,'*****************************************')
 2100    FORMAT(/,1X,'     *****   ITERATION OVER THE SHIFT ',
     +               'NOT CONVERGED')
 2150    FORMAT(/,1X,'COMPARE CURRENTS :',
     +          /,35X,'(1) FROM SOURCE TERM   =  ',1PE13.5,
     +          /,35X,'(2) FROM LINE INTEGRAL =  ',1PE13.5,
     +            10X,'--->  RATIO  (1)/(2)  :  ',1PE13.5)
 2300    FORMAT(///,1X,'*******************************************',
     ,           //,1X,'QUANTITIES COMPUTED IN MAPPIN FROM SOLUTION',
     ,           //,1X,'*******************************************')
 2301    FORMAT(///,1X,'**************************',
     ,           //,1X,'VOLUME AVERAGED QUANTITIES',
     ,           //,1X,'**************************')
 2302    FORMAT(///,1X,'****************************************',
     ,           //,1X,'QUANTITIES EXTRAPOLATED ON MAGNETIC AXIS',
     ,           //,1X,'****************************************')
 2350    FORMAT(///,1X,'*********************************',
     +         //,1X,'FUNCTIONS OF S = SQRT(PSI/PSIMIN) ',
     +         //,1X,'VALUES AT MIDDLE POINTS IN S MESH',
     +         //,1X,'*********************************')
 2400    FORMAT(///,1X,'**************************',
     +         //,1X,'SIZES OF COMMONS FOR CLEAR',
     +         //,1X,'**************************')
 2500    FORMAT(///,1X,'*******************************************',
     ,           //,1X,'QUANTITIES COMPUTED IN MAPPIN TO FEED ERATO',
     ,           //,1X,'*******************************************')
 2501    FORMAT(///,10X,'*******************************************',
     ,           //,10X,'P'' = ',1F10.4,' % OF BALOONING OPTIMIZED P'' '
     ,           //,10X,'*******************************************')
C
         END
C*DECK C3SA02
C*CALL PROCESS
         SUBROUTINE PRIQQU
C        #################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SA02  PRINT EQUILIBRIUM QUANTITIES AT Q-VALUES DEFINED BY QPLACS  *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         DIMENSION
     R     ZBETQ(11),   ZDEDR(11),   ZDI(11),     ZDPRIM(11),  ZELL(11),   
     R     ZHGG(11),    ZLI(11),     ZMERC(11),   ZQR(11),    
     R     ZQS(11),     ZRINT(11),   ZRBAL(11),   ZRMERC(11),  
     R     ZRRINT(11),  ZS(11),      ZSHR(11),    ZVOL(11),    
     R     ZSBAL(11),   ZSMERC(11),  ZSRINT(11),  ZSY(11)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         IF (NMESHA .NE. 2) GOTO 2
C
         WRITE(6,901)
C
         DO 1 J1=1,NPOIDA
C
         IS  = ISRCHFGE(NPSI1,CS ,1,APLACE(J1))
         ISM = ISRCHFGE(NPSI1,CSM,1,APLACE(J1))
C
         IF (IS .LT. 1)    IS = 1
         IF (IS .GT. NPSI) IS = NPSI
C
         IF (ISM .LT. 1)    ISM = 1
         IF (ISM .GT. NPSI) ISM = NPSI
C
         ZBETQ(J1) = BETAB(IS) + (BETAB(IS+1) - BETAB(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
         ZLI(J1)   = RINDUC(IS) + (RINDUC(IS+1) - RINDUC(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
         ZDPRIM(J1)= DPRIME(IS) + (DPRIME(IS+1) - DPRIME(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
         ZQS(J1)   = QPSI(ISM) + (QPSI(ISM+1) - QPSI(ISM)) *
     *               (APLACE(J1) - CSM(ISM)) / (CSM(ISM+1) - CSM(ISM))
         ZSHR(J1)  = CDRQ(IS) + (CDRQ(IS+1) - CDRQ(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
         ZVOL(J1)  = RSURF(IS) + (RSURF(IS+1) - RSURF(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
         ZDI(J1)   = RDI(IS) + (RDI(IS+1) - RDI(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
         ZSY(J1)   = RSY(IS) + (RSY(IS+1) - RSY(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
         ZMERC(J1) = SMERCI(ISM) + (SMERCI(ISM+1) - SMERCI(ISM)) *
     *               (APLACE(J1) - CSM(ISM)) / (CSM(ISM+1) - CSM(ISM))
         ZRINT(J1) = SMERCR(ISM) + (SMERCR(ISM+1) - SMERCR(ISM)) *
     *               (APLACE(J1) - CSM(ISM)) / (CSM(ISM+1) - CSM(ISM))
         ZHGG(J1)  = HMERCR(ISM) + (HMERCR(ISM+1) - HMERCR(ISM)) *
     *               (APLACE(J1) - CSM(ISM)) / (CSM(ISM+1) - CSM(ISM))
         ZELL(J1)  = RELL(ISM) + (RELL(ISM+1) - RELL(ISM)) *
     *               (APLACE(J1) - CSM(ISM)) / (CSM(ISM+1) - CSM(ISM))
         ZDEDR(J1) = RDEDR(IS) + (RDEDR(IS+1) - RDEDR(IS)) *
     *               (APLACE(J1) - CS(IS)) / (CS(IS+1) - CS(IS))
C
    1    CONTINUE
C
         CALL RARRAY(' RATIONAL SURFACE Q',ZQS,NPOIDA)
         CALL RARRAY(' CS(Q)',APLACE,NPOIDA)
         CALL RARRAY(' POLOIDAL BETA(Q)',ZBETQ,NPOIDA)
         CALL RARRAY(' SHEAR(Q)',ZSHR,NPOIDA)
         CALL RARRAY(' SQRT(V/V-TOTAL)(Q)',ZVOL,NPOIDA)
         CALL RARRAY(' ELLIPTICITY',ZELL,NPOIDA)
         CALL RARRAY(' D(E) / DR',ZDEDR,NPOIDA)
         CALL RARRAY(' INDUCTANCE',ZLI,NPOIDA)
         CALL RARRAY(' DELTA PRIME (LOWEST A/R0 ORDER)',ZDPRIM,NPOIDA)
         CALL RARRAY(' MERCIER BY A/R0 EXPANSION WITH E',ZDI,NPOIDA)
         CALL RARRAY(' MERCIER BY SHAFRANOV YOURCHENKO',ZSY,NPOIDA)
C
         IF (NBAL .EQ. 1) THEN
C
            CALL RARRAY(' MERCIER(Q)',ZMERC,NPOIDA)
            CALL RARRAY(' RESISTIVE INTERCHANGE(Q)',ZRINT,NPOIDA)
            CALL RARRAY(' H OF GGG(Q)',ZHGG,NPOIDA)
C
         ENDIF
C
    2    CONTINUE
C
         IF (NBAL .EQ. 1) THEN
C
            WRITE(6,902)
C
            IMERC = 0
            IRINT = 0
            IBAL  = 0
C
            DO 3 J3=1,NPSI
C
            IF (SMERCI(J3) * SMERCI(J3+1) .LT. 0.) THEN
C
               IMERC         = IMERC + 1
               ZSMERC(IMERC) = CSM(J3) + SMERCI(J3) *
     *                         (CSM(J3) - CSM(J3+1)) / 
     /                         (SMERCI(J3+1) - SMERCI(J3))
               ZRMERC(IMERC) = RSURF(J3) + SMERCI(J3) *
     *                         (RSURF(J3) - RSURF(J3+1)) / 
     /                         (SMERCI(J3+1) - SMERCI(J3))
C
            ENDIF
C
            IF (SMERCR(J3) * SMERCR(J3+1) .LT. 0.) THEN
C
               IRINT         = IRINT + 1
               ZSRINT(IRINT) = CSM(J3) + SMERCR(J3) *
     *                         (CSM(J3) - CSM(J3+1)) / 
     /                         (SMERCR(J3+1) - SMERCR(J3))
               ZRRINT(IRINT) = RSURF(J3) + SMERCR(J3) *
     *                         (RSURF(J3) - RSURF(J3+1)) / 
     /                         (SMERCR(J3+1) - SMERCR(J3))
C
            ENDIF
C
            IF (NCBAL(J3) + NCBAL(J3+1) .EQ. 1) THEN
C
               IBAL        = IBAL + 1
               ZSBAL(IBAL) = .5 * (CSM(J3) + CSM(J3+1))
               ZRBAL(IBAL) = .5 * (RSURF(J3) + RSURF(J3+1))
C
            ENDIF
C
    3       CONTINUE
C
            IF (IMERC .GT. 0) THEN
C
               CALL RARRAY('S VALUES FOR WHICH -DI = 0',ZSMERC,IMERC)
               CALL RARRAY('R VALUES FOR WHICH -DI = 0',ZRMERC,IMERC)
C
            ENDIF
C
            IF (IRINT .GT. 0) THEN
C
               CALL RARRAY('S VALUES FOR WHICH -DR = 0',ZSRINT,IRINT)
               CALL RARRAY('R VALUES FOR WHICH -DR = 0',ZRRINT,IRINT)
C
            ENDIF
C
            IF (IBAL .GT. 0) THEN
C
               CALL RARRAY('S VALUES WHERE BAL. CRIT. = 0',ZSBAL,IBAL)
               CALL RARRAY('R VALUES WHERE BAL. CRIT. = 0',ZRBAL,IBAL)
C
            ENDIF
         ENDIF
C
         ZSTEP = DFLOAT(1) / DFLOAT(10)
C
         DO 4 J4=1,11
C
         ZS(J4) = (J4 - 1.) * ZSTEP
C
         IS = ISRCHFGE(NPSI1,CSM,1,ZS(J4))
         IR = ISRCHFGE(NPSI1,RSURF,1,ZS(J4))
C
         IF (IS .LT. 1)    IS = 1
         IF (IS .GT. NPSI) IS = NPSI
         IF (IR .LT. 1)    IR = 1
         IF (IR .GT. NPSI) IR = NPSI
C
         ZQS(J4)   = QPSI(IS) + (QPSI(IS+1) - QPSI(IS)) *
     *               (ZS(J4) - CSM(IS)) / (CSM(IS+1) - CSM(IS))
         ZQR(J4)   = QPSI(IR) + (QPSI(IR+1) - QPSI(IR)) *
     *               (ZS(J4) - RSURF(IR)) / (RSURF(IR+1) - RSURF(IR))
         ZELL(J4)  = RELL(IR) + (RELL(IR+1) - RELL(IR)) *
     *               (ZS(J4) - RSURF(IR)) / (RSURF(IR+1) - RSURF(IR))
         ZDEDR(J4) = RDEDR(IR) + (RDEDR(IR+1) - RDEDR(IR)) *
     *               (ZS(J4) - RSURF(IR)) / (RSURF(IR+1) - RSURF(IR))
C
    4    CONTINUE
C
         WRITE(6,903)
C
         CALL RARRAY('S-VALUES',ZS,11)
         CALL RARRAY('Q-VALUES',ZQS,11)    
C
         WRITE(6,904)
C
         CALL RARRAY('R-VALUES',ZS,11)
         CALL RARRAY('Q-VALUES',ZQR,11)    
         CALL RARRAY('ELLIPTICITY',ZELL,11)    
         CALL RARRAY('D(E) / DR',ZDEDR,11)    
C
  901    FORMAT(///,1X,'******************************',
     ,           //,1X,'OUTPUT ON RATIONAL Q SURFACES ',
     ,           //,1X,'******************************')
  902    FORMAT(///,1X,'****************************************',
     ,           //,1X,'ZEROES OF LOCAL AND BALLOONING CRITERIA ',
     ,           //,1X,'****************************************')
  903    FORMAT(///,1X,'****************************',
     ,           //,1X,'Q VALUES FOR GIVEN S-VALUES ',
     ,           //,1X,'****************************')
  904    FORMAT(///,1X,'****************************',
     ,           //,1X,'Q VALUES FOR GIVEN R-VALUES ',
     ,           //,1X,'****************************')
C
         RETURN
         END
         SUBROUTINE OUTMKSA(KUNIT,KOPT)
C        ##############################
C
C                                        AUTHORS:
C                                        O.SAUTER,  CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SA03 OUTPUT OF SOME EQUILIBRIUM VALUES AND THEIR MKSA VALUES      *
*        USING R0EXP [M] AND B0EXP [T] FOR CONVERSION                 *
*                                                                     *
*        KUNIT: DISK UNIT ON WHICH TO WRITE                           *
*        KOPT = 1: NO HEADER                                          *
*        KOPT = 2: WITH HEADER                                        *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMBND.inc'
C
C----------------------------------------------------------------------
C
CL       1. WRITE SOME VALUES IN MKSA USING R0EXP [M] AND B0EXP [T]
C
         IUNIT = KUNIT
C
         WRITE(IUNIT,*)
         IF (KOPT .EQ. 2) THEN
           WRITE(IUNIT,9100)
         ENDIF
         IF (NRSCAL .EQ. 0) THEN
           WRITE(IUNIT,9101) R0EXP,' R0 [M] USED FOR CONVERTING TO MKSA'
           WRITE(IUNIT,9101) B0EXP,' B0 [T] USED FOR CONVERTING TO MKSA'
         ELSE
           WRITE(IUNIT,9103) R0EXP,' R0 [M] USED FOR CONVERTING TO MKSA'
     +       ,' (CHANGED TO RMAG_EXP)'
           WRITE(IUNIT,9103) B0EXP,' B0 [T] USED FOR CONVERTING TO MKSA'
     +       ,' (CHANGED TO BMAG_EXP)'
         ENDIF
         WRITE(IUNIT,9102) RMAG,' R OF MAGAXE --> [M]   ',RMAG*R0EXP
         WRITE(IUNIT,9102) RZMAG,' Z OF MAGAXE --> [M]   ',RZMAG*R0EXP
         WRITE(IUNIT,9102) SPSIM, ' PSI-AXIS --> [T M**2] ',
     +                    SPSIM*B0EXP*R0EXP**2
         WRITE(IUNIT,9102) 6.28319*SPSIM,' 2*PI*PSI-AXIS -->     ',
     +                    6.28319*SPSIM*B0EXP*R0EXP**2
         WRITE(IUNIT,9102) RITOT, ' TOTAL CURRENT --> [A] ',
     +                    RITOT*R0EXP*B0EXP/4.E-07/3.141593
         WRITE(IUNIT,9101) BETAP,' POLOIDAL BETA'
         WRITE(IUNIT,9101) BETA,' BETA'
         WRITE(IUNIT,9101) BETAS,' BETA* (SQRT(<P**2>))'
         WRITE(IUNIT,9101) BETAX,' BETA_EXP=<P>*2*MU0/B0**2'
         WRITE(IUNIT,9101) RINDUC(NPSI1),' LI'
         WRITE(IUNIT,9101) Q0,' Q_ZERO'
         WRITE(IUNIT,9101) QPSI(NPSI1),' Q_EDGE'
         WRITE(IUNIT,9102) ARATIO(NPSI1),' ASPECT RATIO ; a/R= ',
     +     1./ARATIO(NPSI1)
         WRITE(IUNIT,9101) (1.+RELL(NPSI1))/(1.-RELL(NPSI1)),' b/a'
         WRITE(IUNIT,9102) VOLUME, ' VOLUME -> ',VOLUME*R0EXP**3
         WRITE(IUNIT,9102) AREA,   ' AREA   -> ',AREA*R0EXP**2
         WRITE(IUNIT,9102)RLENG(NPSI1), ' LENGTH -> ',RLENG(NPSI1)*R0EXP
         IRMAX = ISMAX(NBPSOUT,RRBPSOU,1)
         IRMIN = ISMIN(NBPSOUT,RRBPSOU,1)
         IZMAX = ISMAX(NBPSOUT,RZBPSOU,1)
         IZMIN = ISMIN(NBPSOUT,RZBPSOU,1)
         ZRMIN = RRBPSOU(IRMIN)
         ZRMAX = RRBPSOU(IRMAX)
         WRITE(IUNIT,9102) ZRMIN,' RMIN -> RMIN [m] ',ZRMIN*R0EXP
         WRITE(IUNIT,9102) ZRMAX,' RMAX -> RMAX [m] ',ZRMAX*R0EXP
         WRITE(IUNIT,9102) RZBPSOU(IZMIN),' ZMIN -> ZMIN [m] ',
     +       RZBPSOU(IZMIN)*R0EXP
         WRITE(IUNIT,9102) RZBPSOU(IZMAX),' ZMAX -> ZMAX [m] ',
     +       RZBPSOU(IZMAX)*R0EXP
         WRITE(IUNIT,9102) 0.5*(ZRMIN+ZRMAX),' RGEOM -> RGEOM [m] ',
     +       0.5*(ZRMIN+ZRMAX)*R0EXP
         WRITE(IUNIT,9102) 0.5*(ZRMAX-ZRMIN),' MINOR RADIUS -> A [m] ',
     +       0.5*(ZRMAX-ZRMIN)*R0EXP
C
         RETURN
 9100    FORMAT(/,1X,'*************************************',
     ,         //,1X,'SOME QUANTITIES AND THEIR MKSA VALUES',
     ,         //,1X,'*************************************',/)
 9101    FORMAT(1PE18.8,A)
 9102    FORMAT(1PE18.8,A,E18.8)
 9103    FORMAT(1PE18.8,2A)
C
         END
C*DECK C3SB01
C*CALL PROCESS
         SUBROUTINE IODISK(K)
C        ####################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SB01  PERFORM DISK OPERATIONS                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMLAB.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMPLO.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
         INCLUDE 'COMEQD.inc'
         INCLUDE 'COMDAT.inc'
C
         DIMENSION
     R   ZBND(24*NPT),      ZCPR(NPPSI1),      ZCPPR(NPPSI1),
     R   ZPSI(NPPSI1),      ZR(12*NPT+1),
     R   ZT(12*NPT+1),      ZTMF(NPPSI1),      ZTTP(NPPSI1),     
     R   ZZ(12*NPT+1),      ZCID0(NPPSI1),     ZCID2(NPPSI1),
     R   ZWORK(NPISO),      ZSTEMP(NPISO),     ZTEMP(NPISO),
     R   ZD2TMP(NPISO),     ZWORK1(NPISO)
         CHARACTER  ZDATE*8
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     G          170,180,190,200,210,220,230,240,250,260,270,280,290,
     G          300,310,320,330,340,350,360,370,380) K
C
C----------------------------------------------------------------------
C  1. OPEN UNIT NSAVE
C
 10      CONTINUE
C
         OPEN(UNIT=NSAVE,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='NSAVE')
         REWIND NSAVE
C
         RETURN
C
C----------------------------------------------------------------------
C  2.  OPEN UNIT NVAC
C
 20      CONTINUE
C
         OPEN(UNIT=NVAC,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     O        FILE='NVAC')
C
         RETURN
C
C----------------------------------------------------------------------
C  3.   OPEN UNIT NDES
C
 30      CONTINUE
C
         IF (NIDEAL .EQ. 0) THEN
C
            OPEN(UNIT=NDES,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O           FILE='NDES')
C
         ELSE
C
            OPEN(UNIT=NDES,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     O           FILE='NDES')
C
         ENDIF
C
         RETURN
C
C----------------------------------------------------------------------
C  4.   OPEN UNIT MEQ
C
 40      CONTINUE
C
         IF (NIDEAL .EQ. 1) THEN
C
            OPEN(UNIT=MEQ,ACCESS='DIRECT',RECL=8*NDEQ*NCHI,
     O           FORM='UNFORMATTED',FILE='MEQ')
C
         ELSE IF (NIDEAL .EQ. 2) THEN
C
            OPEN(UNIT=MEQ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     O           FILE='MEQ')
C
         ENDIF
C
         RETURN
C
C----------------------------------------------------------------------
C  5.   OPEN UNIT NPRNT
C
 50      CONTINUE
C
         OPEN(UNIT=NPRNT,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='NPRNT')
C
         RETURN
C
C----------------------------------------------------------------------
C  6.   OPEN UNIT NOUT
C
 60      CONTINUE
C
         OPEN(UNIT=NOUT,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     O        FILE='NOUT')
C
         RETURN
C
C----------------------------------------------------------------------
C  7.   OPEN UNIT NIN
C
 70      CONTINUE
C
         OPEN(UNIT=NIN,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     O        STATUS='OLD',FILE='NIN')
         REWIND NIN
C
         RETURN
C
C----------------------------------------------------------------------
C  8.   OPEN UNIT NO AND NOI
C
 80      CONTINUE
C
         OPEN(UNIT=NO,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='EQU01')
         OPEN(UNIT=NOI,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='EQU02')
         OPEN(UNIT=NO3,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='EQU03')
C
         RETURN
C
C----------------------------------------------------------------------
C  9.   OPEN UNIT NUPLO
C
 90      CONTINUE
C
         OPEN(UNIT=NUPLO,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='NUPLO')
C
         RETURN
C
C----------------------------------------------------------------------
C 10.  CLOSE UNIT NSAVE
C
 100     CONTINUE
C
         CLOSE(UNIT=NSAVE,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  11.  CLOSE UNIT NVAC
C
 110     CONTINUE
C
         CLOSE(UNIT=NVAC,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  12.  CLOSE UNIT NDES
C
 120     CONTINUE
C
         CLOSE(UNIT=NDES,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  13.  CLOSE UNIT MEQ
C
 130     CONTINUE
C
         CLOSE(UNIT=MEQ,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  14.  CLOSE UNIT NPRNT
C
 140     CONTINUE
C
         CLOSE(UNIT=NPRNT,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  15.  CLOSE UNIT NOUT
C
 150     CONTINUE
C
         CLOSE(UNIT=NOUT,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  16.  CLOSE UNIT NIN
C
 160     CONTINUE
C
         CLOSE(UNIT=NIN,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  17.  CLOSE UNIT NO AND NOI
C
 170     CONTINUE
C
         CLOSE(UNIT=NO,STATUS='KEEP')
         CLOSE(UNIT=NOI,STATUS='KEEP')
         CLOSE(UNIT=NO3,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  18. CLOSE UNIT NUPLO
C
 180     CONTINUE
C
         CLOSE(UNIT=NUPLO,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  19. WRITE EQUILIBRIUM QUANTITIES EQ
C
 190     CONTINUE
C
         IF (NIDEAL .EQ. 1) THEN
C
            DO 191 J191=1,NPSI
C
            WRITE(UNIT=MEQ,REC=J191) ((EQ(I,J,J191),I=1,NDEQ),J=1,NCHI)
C
  191       CONTINUE
C
         ELSE IF (NIDEAL .EQ. 2) THEN
C
            DO 192 J192=1,NPSI1
C
            WRITE(MEQ) ((EQ(I,J,J192),I=1,NDEQ),J=1,NCHI)
C
  192       CONTINUE
C
         ENDIF
C
         RETURN
C
C----------------------------------------------------------------------
C  20.   STILL OPEN OPTION
C
  200    CONTINUE
C
         RETURN
C
C----------------------------------------------------------------------
C  21. PLOT QUANTITIES FOR ERATO
C
  210    CONTINUE
C
         IF (NIDEAL .EQ. 0) THEN
C
            WRITE(NDES,1015) NPSI
            WRITE(NDES,1015) NCHI
            WRITE(NDES,1016) R0
            WRITE(NDES,1016) RZ0
C
         ENDIF
C
         DO 211 J211=1,NPSI
C
         IF (NIDEAL .NE. 0) THEN
C
            WRITE(UNIT=NDES) (CR(J,J211),J=1,NCHI)
            WRITE(UNIT=NDES) (CZ(J,J211),J=1,NCHI)
            WRITE(UNIT=NDES) (CNR1(J,J211),J=1,NCHI)
            WRITE(UNIT=NDES) (CNZ1(J,J211),J=1,NCHI)
C
         ELSE
C
            WRITE(NDES,1014) (CR(J,J211),J=1,NCHI)
            WRITE(NDES,1014) (CZ(J,J211),J=1,NCHI)
            WRITE(NDES,1014) (CNR1(J,J211),J=1,NCHI)
            WRITE(NDES,1014) (CNZ1(J,J211),J=1,NCHI)
            WRITE(NDES,1014) (CNR2(J,J211),J=1,NCHI)
            WRITE(NDES,1014) (CNZ2(J,J211),J=1,NCHI)
C
         ENDIF
C
  211    CONTINUE
C
         RETURN
C
C----------------------------------------------------------------------
C  22. SURFACE AND LABEL FOR PLOT IN ERATO
C
 220     CONTINUE
C
         INUM = 12 * NT
         ZDT  = 2. * CPI / DFLOAT(INUM)
C
         DO 221 J221=1,INUM
C
         ZT(J221) = .5 * (2*J221 - 3) * ZDT
C
 221     CONTINUE
C
         CALL BOUND(INUM,ZT,ZBND)
C
         DO 222 J222=1,INUM
C
         ZR(J222) = R0  + ZBND(J222) * COS(ZT(J222))
         ZZ(J222) = RZ0 + ZBND(J222) * SIN(ZT(J222))
C
 222     CONTINUE
C
         IF (NIDEAL .NE. 0) THEN
C
            WRITE(UNIT=NDES) INUM
            WRITE(UNIT=NDES) (ZR(J),J=1,INUM),(ZZ(J),J=1,INUM)
            WRITE(UNIT=NDES) LABEL1,LABEL2
            WRITE(UNIT=NDES) LABEL2,LABEL3
            WRITE(UNIT=NDES) LABEL3
            WRITE(UNIT=NDES) LABEL4
C
         ELSE
C
            WRITE(NDES,1015) INUM
            WRITE(NDES,1014) (ZR(J),J=1,INUM),(ZZ(J),J=1,INUM)
C
         ENDIF
C
         RETURN
C
C----------------------------------------------------------------------
C  23. VACUM QUANTITIES FOR ERATO
C
 230     CONTINUE
C
         ZQS = QPSI(NPSI1)
         ZTM = TMF(NPSI1)
C
         WRITE(UNIT=NVAC) (CR(J,NPSI1),J=1,NCHI),
     ,                    (CZ(J,NPSI1),J=1,NCHI),
     ,                    (CHI(J),J=1,NCHI1)
         WRITE(UNIT=NVAC) ZQS,ZTM
         WRITE(UNIT=NVAC) (CNR1(J,NPSI1),J=1,NCHI),
     ,                    (CNZ1(J,NPSI1),J=1,NCHI)
C
C     NOTE: SHOULD BE THE SAME IF STATEMENT AS IN ERDATA, WHICH DECIDES
C     WHETHER OR NOT THE EQ'S ARE DEFINED
         IF (NDEQ .GE. 25) THEN
           WRITE(UNIT=NVAC) (EQ(19,J,NPSI1),J=1,NCHI)
           WRITE(UNIT=NVAC) (EQ(20,J,NPSI1),J=1,NCHI)
         ELSE
           WRITE(6,'(/,"  WARNING IN IODISK: EQ(19) AND EQ(20) NOT",
     +       " WRITTEN ON NVAC AS NDEQ<25",/)')
         ENDIF
         WRITE(UNIT=NVAC) (CHIOLD(J),J=1,NCHI1)
C
         IF (NIDEAL .EQ. 1) RETURN
C
         WRITE(UNIT=NVAC) (TETVAC(J),J=1,NCHI1)
         WRITE(UNIT=NVAC) (RHOVAC(J),J=1,NCHI1)
         WRITE(UNIT=NVAC) (TETVACM(J),J=1,NCHI1)
         WRITE(UNIT=NVAC) (RHOVACM(J),J=1,NCHI1)
         WRITE(UNIT=NVAC) (TETVACI(J,1),J=1,NCHI)
         WRITE(UNIT=NVAC) (TETVACI(J,2),J=1,NCHI)
         WRITE(UNIT=NVAC) (TETVACI(J,3),J=1,NCHI)
         WRITE(UNIT=NVAC) (RHOVACI(J,1),J=1,NCHI)
         WRITE(UNIT=NVAC) (RHOVACI(J,2),J=1,NCHI)
         WRITE(UNIT=NVAC) (RHOVACI(J,3),J=1,NCHI)
         WRITE(UNIT=NVAC) (DRHOPI(J,1),J=1,NCHI)
         WRITE(UNIT=NVAC) (DRHOPI(J,2),J=1,NCHI)
         WRITE(UNIT=NVAC) (DRHOPI(J,3),J=1,NCHI)
C
         RETURN
C
C----------------------------------------------------------------------
C  24. STORE THE EQUILIBRIUM
C
 240     CONTINUE
C
         REWIND NOUT
         WRITE(NOUT) NSURF,THETA0,ASPCT,BEANS,CETA,DELTA,ELONG,RNU,
     W               SGMA,TRIANG,TRIPLT,XI,R0EXP,B0EXP
         WRITE(NOUT) NS,NT,NS1,NT1,NSTMAX,N4NSNT,NSOUR,NPP,NSTTP,NIPR,
     W               NISO,NDIFT,NPPR,NBLOPT,CFBAL,NFUNC,NPPFUN,NBSFUN,
     W               NBSOPT,NBSTRP
         WRITE(NOUT) CUROLD,SCALAC,BSFRAC,RZION,ETAEI,PREDGE
         WRITE(NOUT) (CSIG(J),J=1,NS1)
         WRITE(NOUT) (CT(J),J=1,NT1)
         WRITE(NOUT) (BPS(J),J=1,12)
         WRITE(NOUT) (CSIPR(J),J=1,NISO)
         WRITE(NOUT) (CSIPRI(J),J=1,NISO)
         WRITE(NOUT) (CPSICL(J),J=1,N4NSNT)
         WRITE(NOUT) SPSIM,RMAG,RZMAG,R0,RZ0
         WRITE(NOUT) (CPPR(J),J=1,NISO)
         WRITE(NOUT) (CID0(J),J=1,NISO)
         WRITE(NOUT) (CID2(J),J=1,NISO)
         WRITE(NOUT) (D2CID0(J),J=1,NISO)
         WRITE(NOUT) (D2CID2(J),J=1,NISO)
         WRITE(NOUT) (D2CPPR(J),J=1,NISO)
C
         WRITE(NOUT) (TTP(J),J=1,NISO)
         WRITE(NOUT) (TMF(J),J=1,NISO)
         WRITE(NOUT) (D2TMF(J),J=1,NISO)
C
         WRITE(NOUT) (AT(J),J=1,10)
         WRITE(NOUT) (AT2(J),J=1,10)
         WRITE(NOUT) (AT3(J),J=1,10)
         WRITE(NOUT) (AT4(J),J=1,10)
         WRITE(NOUT) (AP(J),J=1,10)
         WRITE(NOUT) (AP2(J),J=1,10)
         WRITE(NOUT) (AFBS(J),J=1,10)
         WRITE(NOUT) (AFBS2(J),J=1,10)
         WRITE(NOUT) NMESHB,NMESHC,NMESHD,SOLPDB,SOLPDC,SOLPDD,
     W               NPOIDB,NPOIDC,NPOIDD
         WRITE(NOUT) (BPLACE(J),J=1,10)
         WRITE(NOUT) (BWIDTH(J),J=1,10)
         WRITE(NOUT) (CPLACE(J),J=1,10)
         WRITE(NOUT) (CWIDTH(J),J=1,10)
         WRITE(NOUT) (DPLACE(J),J=1,10)
         WRITE(NOUT) (DWIDTH(J),J=1,10)
C
         IF (NBLOPT .EQ. 2 .OR. NBSOPT .EQ. 2)  THEN
C
            WRITE(NOUT) (PCS(J),J=1,NPPR+1)
            WRITE(NOUT) (PCSM(J),J=1,NPPR+1)
            WRITE(NOUT) (RPRM(J),J=1,NPPR+1)
            WRITE(NOUT) (D2RPRM(J),J=1,NPPR+1)
C
         ENDIF    
C
         IF (NSURF .EQ. 6) THEN
C
            WRITE(NOUT) NBPS,NWBPS
            WRITE(NOUT) (RW(J),J=1,NWBPS)
            WRITE(NOUT) ((TETBPS(L,J),L=1,NBPS),J=1,NWBPS)     
            WRITE(NOUT) ((RRBPS(L,J),L=1,NBPS),J=1,NWBPS)
            WRITE(NOUT) ((RZBPS(L,J),L=1,NBPS),J=1,NWBPS)     
            WRITE(NOUT) ((D2RBPS(L,J),L=1,NBPS),J=1,NWBPS)          
            WRITE(NOUT) ((D2ZBPS(L,J),L=1,NBPS),J=1,NWBPS)          
C
         ELSE IF (NSURF .EQ. 7) THEN
            WRITE(NOUT) NFOURPB
            WRITE(NOUT) ALZERO, RC
            WRITE(NOUT) (BPSCOS(L),BPSSIN(L),L=1,NFOURPB)
C
         ENDIF
C
         IF (NFUNC .EQ. 4 .OR. NPPFUN .EQ. 4) THEN
C
            NPPF1 = NPPF + 1
            WRITE(NOUT) NPPF1
            WRITE(NOUT) (FCSM(L),L=1,NPPF1)
            WRITE(NOUT) (RPPF(L),L=1,NPPF1)
            WRITE(NOUT) (RFUN(L),L=1,NPPF1)
            WRITE(NOUT) (D2RPPF(L),L=1,NPPF1)
            WRITE(NOUT) (D2RFUN(L),L=1,NPPF1)
C
         ENDIF
C
         RETURN
C
C----------------------------------------------------------------------
C  25. READ LABELS AND NAMELIST FROM UNIT NSAVE
C
 250     CONTINUE
C
         REWIND NSAVE
         READ (UNIT=NSAVE,FMT=1000) LABEL1
         READ (UNIT=NSAVE,FMT=1000) LABEL2
         READ (UNIT=NSAVE,FMT=1000) LABEL3
         READ (UNIT=NSAVE,FMT=1000) LABEL4
         READ(NSAVE,EQDATA)
C
         RETURN
C
C----------------------------------------------------------------------
C  26. READ THE EQUILIBRIUM
C
 260     CONTINUE
C
         REWIND NIN
C
         IF (NOPT.GE.-2 .AND. NOPT.LE.-1) GO TO 262
C
C-----------------------------------------------------------------------
C     26.1 READ FULL EQUILIBRIUM QUANTITIES
C
         READ(NIN) NSURF,THETA0,ASPCT,BEANS,CETA,DELTA,ELONG,RNU,
     W             SGMA,TRIANG,TRIPLT,XI,ZR0EXP,ZB0EXP
C     USE NAMELIST VALUES IF NEGATIVE
         IF (R0EXP .GE. 0.0) R0EXP = ZR0EXP
         IF (B0EXP .GE. 0.0) B0EXP = ZB0EXP
         R0EXP = ABS(R0EXP)
         B0EXP = ABS(B0EXP)
C
         ISOUROLD = NSOUR
         READ(NIN) NS,NT,NS1,NT1,NSTMAX,N4NSNT,NSOUR,NPP,NSTTP,NIPR,
     W             NISO,NDIFT,NPPR,NBLOPT,CFBAL,NFUNC,NPPFUN,NBSFUN,
     W             NBSOPT,NBSTRP
         READ(NIN) CUROLD,SCALAC,BSFRAC,RZION,ETAEI,PREDGE
         READ(NIN) (CSIG(J),J=1,NS1)
         READ(NIN) (CT(J),J=1,NT1)
         READ(NIN) (BPS(J),J=1,12)
         IF (NSURF .NE. 6) RC = BPS(2)
         IF (NSURF .EQ. 7) RZ0C = BPS(6)
C
         READ(NIN) (CSIPR(J),J=1,NISO)
         READ(NIN) (CSIPRI(J),J=1,NISO)
         READ(NIN) (CPSICL(J),J=1,N4NSNT)
         READ(NIN) SPSIM,RMAG,RZMAG,R0,RZ0
         READ(NIN) (CPPR(J),J=1,NISO)
         READ(NIN) (CID0(J),J=1,NISO)
         READ(NIN) (CID2(J),J=1,NISO)
         READ(NIN) (D2CID0(J),J=1,NISO)
         READ(NIN) (D2CID2(J),J=1,NISO)
         READ(NIN) (D2CPPR(J),J=1,NISO)
C
         READ(NIN) (TTP(J),J=1,NISO)
         READ(NIN) (TMF(J),J=1,NISO)
         CALL SCOPY(NISO,TMF,1,TMFO,1)
         READ(NIN) (D2TMF(J),J=1,NISO)
C
         READ(NIN) (AT(J),J=1,10)
         READ(NIN) (AT2(J),J=1,10)
         READ(NIN) (AT3(J),J=1,10)
         IF (AT4(1).NE.0. .AND. NFUNC.NE.2) THEN
           READ(NIN) (AT3(J),J=1,10)
           NSOUR = ISOUROLD
         ELSE
           READ(NIN) (AT4(J),J=1,10)
         ENDIF
         READ(NIN) (AP(J),J=1,10)
         READ(NIN) (AP2(J),J=1,10)
         READ(NIN) (AFBS(J),J=1,10)
         READ(NIN) (AFBS2(J),J=1,10)
         READ(NIN) NMESHB,NMESHC,NMESHD,SOLPDB,SOLPDC,SOLPDD,
     W             NPOIDB,NPOIDC,NPOIDD
         READ(NIN) (BPLACE(J),J=1,10)
         READ(NIN) (BWIDTH(J),J=1,10)
         READ(NIN) (CPLACE(J),J=1,10)
         READ(NIN) (CWIDTH(J),J=1,10)
         READ(NIN) (DPLACE(J),J=1,10)
         READ(NIN) (DWIDTH(J),J=1,10)
C
         IF (NBLOPT .EQ. 2 .OR. NBSOPT .EQ. 2)  THEN
C
            READ(NIN) (PCS(J),J=1,NPPR+1)
            READ(NIN) (PCSM(J),J=1,NPPR+1)
            READ(NIN) (RPRM(J),J=1,NPPR+1)
            READ(NIN) (D2RPRM(J),J=1,NPPR+1)
C
         ENDIF         
C
         IF (NSURF .EQ. 6) THEN
C
            READ(NIN) NBPS,NWBPS
            READ(NIN) (RW(J),J=1,NWBPS)
            READ(NIN) ((TETBPS(L,J),L=1,NBPS),J=1,NWBPS)     
            READ(NIN) ((RRBPS(L,J),L=1,NBPS),J=1,NWBPS)     
            READ(NIN) ((RZBPS(L,J),L=1,NBPS),J=1,NWBPS)          
            READ(NIN) ((D2RBPS(L,J),L=1,NBPS),J=1,NWBPS)          
            READ(NIN) ((D2ZBPS(L,J),L=1,NBPS),J=1,NWBPS)          
C
         ELSE IF (NSURF .EQ. 7) THEN
            READ(NIN) NFOURPB
            READ(NIN) ALZERO, RC
            READ(NIN) (BPSCOS(L),BPSSIN(L),L=1,NFOURPB)
C
         ENDIF
C
         IF (NFUNC .EQ. 4 .OR. NPPFUN .EQ. 4) THEN
C
            READ(NIN) NPPF1
            NPPF = NPPF1 - 1
            READ(NIN) (FCSM(L),L=1,NPPF1)
            READ(NIN) (RPPF(L),L=1,NPPF1)
            READ(NIN) (RFUN(L),L=1,NPPF1)
            READ(NIN) (D2RPPF(L),L=1,NPPF1)
            READ(NIN) (D2RFUN(L),L=1,NPPF1)
C
         ENDIF
C
         RETURN
C
C-----------------------------------------------------------------------
C     26.2 READ ONLY PSI AND RELATED QUATITIES FOR FIRST GUESS (RESTART)
C
 262     CONTINUE
C
C
         READ(NIN) IDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,
     W             ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM
         READ(NIN) NS,NT,NS1,NT1,IDUM,N4NSNT,IDUM,IDUM,IDUM,IDUM,
     W             NISO,IDUM,IDUM,IDUM,ZDUM,IDUM,IDUM,IDUM,
     W             IDUM,IDUM
         READ(NIN) ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM
         READ(NIN) (CSIG(J),J=1,NS1)
         READ(NIN) (CT(J),J=1,NT1)
C     DO NOT READ BPS, TO ALLOW NAMELIST VALUES TO OVERRIDE
         READ(NIN) (CSIPR(J),J=1,12)
         READ(NIN) (CSIPR(J),J=1,NISO)
         READ(NIN) (CSIPRI(J),J=1,NISO)
         READ(NIN) (CPSICL(J),J=1,N4NSNT)
         READ(NIN) SPSIM,RMAG,RZMAG,R0,RZ0
C
         READ(NIN) (CPPR(J),J=1,NISO)
         READ(NIN) (CID0(J),J=1,NISO)
         READ(NIN) (CID2(J),J=1,NISO)
         READ(NIN) (D2CID0(J),J=1,NISO)
         READ(NIN) (D2CID2(J),J=1,NISO)
         READ(NIN) (D2CPPR(J),J=1,NISO)
C
         READ(NIN) (TTP(J),J=1,NISO)
         READ(NIN) (TMF(J),J=1,NISO)
         CALL SCOPY(NISO,TMF,1,TMFO,1)
         READ(NIN) (D2TMF(J),J=1,NISO)
C
         RETURN
C
C----------------------------------------------------------------------
C  27. WRITE QUANTITIES FOR MARS
C
 270     CONTINUE
C
         REWIND NO
         REWIND NOI
         WRITE(NO,1001) NPSI1,MSMAX,NSMAX,ASPCTR,Q0
         WRITE(NO,1002) (CS(I),I=1,NPSI1)
C
         CALL GENOUT(DG11L (1,1),' DG11L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DG22L (1,1),' DG22L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DG33L (1,1),' DG33L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DG12L (1,1),' DG12L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
C
         CALL GENOUT(DG11LM(1,1),'DG11LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(DG22LM(1,1),'DG22LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(DG33LM(1,1),'DG33LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(DG12LM(1,1),'DG12LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
C
         CALL GENOUT(JG11L (1,1),' JG11L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JG22L (1,1),' JG22L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JG33L (1,1),' JG33L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JG12L (1,1),' JG12L',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
C
         CALL GENOUT(JG11LM(1,1),'JG11LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(JG22LM(1,1),'JG22LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(JG33LM(1,1),'JG33LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(JG12LM(1,1),'JG12LM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
C
         CALL GENOUT(JACOBI(1,1),'JACOBI',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JACOBM(1,1),'JACOBM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
C
         CALL GENOUT( B2U(1,1,1),'   B2U',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT( B3U(1,1,1),'   B3U',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT( J2U(1,1,1),'   J2U',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT( J3U(1,1,1),'   J3U',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT( PRE(1,1,1),'   PRE',NO,NPISO,NPSI,MSMAX,1,RM,RN)
C
         CALL GENOUT( B2E(1,1,1),'   B2E',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT( B3E(1,1,1),'   B3E',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT( J2E(1,1,1),'   J2E',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT( J3E(1,1,1),'   J3E',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT( PEQ(1,1,1),'   PEQ',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
C
         CALL GENOUT(DPEDS(1,1,1),' DPEDS',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DPEDSM(1,1,1),'DPEDSM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
C
cab      NEW QUANTITIES FOR INERTIA IN ROTATING PLASMA
cab
         CALL GENOUT(GCHDZ (1,1),' GCHDZ',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GSDZ  (1,1),'  GSDZ',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GBZ   (1,1),'   GBZ',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GBR   (1,1),'   GBR',NO,NPISO,NPSI1,MSMAX,1,RM,RN)
C
         CALL GENOUT(GCHDZM(1,1),'GCHDZM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GSDZM (1,1),' GSDZM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GBZM  (1,1),'  GBZM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GBRM  (1,1),'  GBRM',NO,NPISO,NPSI,MSMAX,1,RM,RN)
C
       IF (.TRUE.) GOTO 275
         CALL GENOUT(IDIY2(1,1) ,' IDIY2',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IDIY2M(1,1),'IDIY2M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IDIY3(1,1) ,' IDIY3',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IDIY3M(1,1),'IDIY3M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IG122(1,1) ,' IG122',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IG122M(1,1),'IG122M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IG123(1,1) ,' IG123',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IG123M(1,1),'IG123M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(INXX(1,1)  ,'  INXX',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(INXXM(1,1) ,' INXXM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(INXY(1,1)  ,'  INXY',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(INXYM(1,1) ,' INXYM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(INYY(1,1)  ,'  INYY',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(INYYM(1,1) ,' INYYM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(INZZ(1,1)  ,'  INZZ',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(INZZM(1,1) ,' INZZM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IJ0QX(1,1) ,' IJ0QX',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IJ0QXM(1,1),'IJ0QXM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IJ0QY(1,1) ,' IJ0QY',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IJ0QYM(1,1),'IJ0QYM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IGPX2(1,1) ,' IGPX2',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IGPX2M(1,1),'IGPX2M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IGPX3(1,1) ,' IGPX3',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IGPX3M(1,1),'IGPX3M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IGPY2(1,1) ,' IGPY2',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IGPY2M(1,1),'IGPY2M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IGPY3(1,1) ,' IGPY3',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IGPY3M(1,1),'IGPY3M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IDRXX(1,1) ,' IDRXX',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IDRXXM(1,1),'IDRXXM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IRXZ(1,1)  ,'  IRXZ',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IRXZM(1,1) ,' IRXZM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IDRYX(1,1) ,' IDRYX',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IDRYXM(1,1),'IDRYXM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IRYX(1,1)  ,'  IRYX',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IRYXM(1,1) ,' IRYXM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IDRZX(1,1) ,' IDRZX',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IDRZXM(1,1),'IDRZXM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IRZY(1,1)  ,'  IRZY',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IRZYM(1,1) ,' IRZYM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(VISXZ(1,1) ,' VISXZ',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(VISXZM(1,1),'VISXZM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(VISYZ(1,1) ,' VISYZ',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(VISYZM(1,1),'VISYZM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IVS11(1,1) ,' IVS11',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IVS11M(1,1),'IVS11M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IVS12(1,1) ,' IVS12',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IVS12M(1,1),'IVS12M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IVS21(1,1) ,' IVS21',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IVS21M(1,1),'IVS21M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IVS22(1,1) ,' IVS22',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IVS22M(1,1),'IVS22M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GSFC(1,1)  ,'  GSFC',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GSFCM(1,1) ,' GSFCM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GSCC(1,1)  ,'  GSCC',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GSCCM(1,1) ,' GSCCM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GSFS(1,1)  ,'  GSFS',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GSFSM(1,1) ,' GSFSM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GSCS(1,1)  ,'  GSCS',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GSCSM(1,1) ,' GSCSM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GCFC(1,1)  ,'  GCFC',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GCFCM(1,1) ,' GCFCM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(GCFS(1,1)  ,'  GCFS',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(GCFSM(1,1) ,' GCFSM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
C
         CALL GENOUT(EQRHO(1,1) ,'   RHO',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(EQRHOM(1,1),'  RHOM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(DRHOS(1,1) ,' DRHOS',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DRHOSM(1,1),'DRHOSM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(EQROT(1,1) ,'   ROT',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(EQROTM(1,1),'  ROTM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT( DROT(1,1) ,'  DROT',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT( DROTM(1,1),' DROTM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(  FEQ(1,1) ,'   FEQ',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(  FEQM(1,1),'  FEQM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IWSQ1(1,1) ,' IWSQ1',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IWSQ1M(1,1),'IWSQ1M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IWSQ2(1,1) ,' IWSQ2',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IWSQ2M(1,1),'IWSQ2M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IWSQ3(1,1) ,' IWSQ3',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IWSQ3M(1,1),'IWSQ3M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IJ0QZ(1,1) ,' IJ0QZ',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IJ0QZM(1,1),'IJ0QZM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(JACOF(1,1) ,' JACOF',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JACOFM(1,1),'JACOFM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(  B2F(1,1) ,'   B2F',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(  B2FM(1,1),'  B2FM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(  B3F(1,1) ,'   B3F',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(  B3FM(1,1),'  B3FM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(JACOS(1,1) ,' JACOS',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JACOSM(1,1),'JACOSM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(IGF22(1,1) ,' IGF22',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(IGF22M(1,1),'IGF22M',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT( B3FC(1,1) ,'  B3FC',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT( B3FCM(1,1),' B3FCM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT( B2FC(1,1) ,'  B2FC',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT( B2FCM(1,1),' B2FCM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
         CALL GENOUT(DJCOF(1,1) ,' DJCOF',NOI,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DJCOFM(1,1),'DJCOFM',NOI,NPISO,NPSI,MSMAX,1,RM,RN)
C
 275     CONTINUE

C LIUYQ
C WRITE QUANTITIES FOR NFTC

C OBS: AT THIS MOMENT QPSI(I=1,NPSI1) CONTAINS Q-VALUES ON HALF-INTEGER
C      GRID, QPSI(I=NPSI1+1,2*NPSI1) CONTAINS TRASH. THEREFORE, QPSI
C      CAN NOT BE USED FOR EQUILIBRIUM Q-VALUES FOR NFTC.  THE SAME
C      APPLIED TO TMF, CP, AND CPR.
   
         REWIND NO3

         WRITE(NO3,2701) NPSI1,MSMAX,ASPCTR 
         WRITE(NO3,1002) (CS(I),I=1,NPSI1)

C Q-PROFILE
         DO I=2,NPSI1
            ZTEMP(I) = REAL(B3E(I,1,1)/B2E(I,1,1))
         ENDDO
         ZTEMP(1) = QPSI(1)
         WRITE(NO3,2702) (ZTEMP(I),I=1,NPSI1)

C PEQ-PROFILE
         WRITE(NO3,2703) (REAL(PEQ(I,1,1)),I=1,NPSI1)

C A-PROFILE
         DO I=2,NPSI1
            ZTEMP(I) = REAL(1.0/DG33L(I,1))
         ENDDO
         ZTEMP(1) = 0.0
         WRITE(NO3,2702) (ZTEMP(I),I=1,NPSI1)

C F-PROFILE
         DO I=2,NPSI1
            ZTEMP(I) = REAL(B3E(I,1,1)*DG33L(I,1))
         ENDDO
         ZTEMP(1) = TMF(1)
         WRITE(NO3,2702) (ZTEMP(I),I=1,NPSI1)

C PSI'-PROFILE
         WRITE(NO3,2706) (REAL(B2E(I,1,1)),I=1,NPSI1)

         FRM(1,1) = R0
         FZM(1,1) = RZ0

         CALL GENOUT(FRM(1,1),'   FRM',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(FZM(1,1),'   FZM',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)

         CALL GENOUT(JG11L(1,1),' JG11L',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JG12L(1,1),' JG12L',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JG22L(1,1),' JG22L',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(JG33L(1,1),' JG33L',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)

         CALL GENOUT(DG11L(1,1),' DG11L',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DG12L(1,1),' DG12L',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)
         CALL GENOUT(DG22L(1,1),' DG22L',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)

         CALL GENOUT(JACOBINV(1,1),'JACOBV',NO3,NPISO,NPSI1,
     &               MSMAX,1,RM,RN)
         CALL GENOUT(JACOBI(1,1),'JACOBI',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)

         CALL GENOUT( J3U(1,1,1),'   J3U',NO3,NPISO,NPSI1,MSMAX,1,RM,RN)
C VACUUM
         WRITE(NO3,1301) NV1
         WRITE(NO3,1002) (CSV(I),I=1,NV1)
         CALL GENOUT(DG11LV(1,1),'DG11LV',NO3,NPV1,NV1,MSMAX,1,RM,RN)
         CALL GENOUT(DG12LV(1,1),'DG12LV',NO3,NPV1,NV1,MSMAX,1,RM,RN)
         CALL GENOUT(DG22LV(1,1),'DG22LV',NO3,NPV1,NV1,MSMAX,1,RM,RN)
         CALL GENOUT(DGRMLV(1,1),'DGRMLV',NO3,NPV1,NV1,MSMAX,1,RM,RN)
         CALL GENOUT(DGZMLV(1,1),'DGZMLV',NO3,NPV1,NV1,MSMAX,1,RM,RN)


C        -----------------------------------------------
C        LIUYQ, 2005-01-29
C        OUTPUT RZ-COORDINATES IN FOURIER SPACE
C        THESE FOURIER HARMONICS ARE NOT ACCURATE, 
C        BETTER TO RUN SEPERATE CODE FourierRF TO GENERATE RMZM_F

         IF (1.EQ.0) THEN
         OPEN(UNIT=NO4,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='RMZM_F')
         WRITE(NO4,1019) MSMAX,NPSI1,NV1-1,R0EXP
         DO I=1,NPSI1
            WRITE(NO4,1020) CS(I),CS(I),CS(I),B0EXP
         ENDDO
         DO I=2,NV1
            WRITE(NO4,1020) CSV(I),CSV(I),CSV(I),CSV(I)
         ENDDO
         DO J=1,MSMAX
            DO I=1,NPSI1
               ZTEMP(1)=REAL(FRM(I,J))
               ZTEMP(2)=IMAG(FRM(I,J))
               ZTEMP(3)=REAL(FZM(I,J))
               ZTEMP(4)=IMAG(FZM(I,J))
               WRITE(NO4,1020) ZTEMP(1),ZTEMP(2),ZTEMP(3),ZTEMP(4)
            ENDDO
            DO I=2,NV1
               ZTEMP(1)=REAL(DGRMLV(I,J))
               ZTEMP(2)=IMAG(DGRMLV(I,J))
               ZTEMP(3)=REAL(DGZMLV(I,J))
               ZTEMP(4)=IMAG(DGZMLV(I,J))
               WRITE(NO4,1020) ZTEMP(1),ZTEMP(2),ZTEMP(3),ZTEMP(4)
            ENDDO
         ENDDO
 1019    FORMAT(I5,1X,I5,1X,I5,1X,E15.8)
 1020    FORMAT(E15.8,1X,E15.8,1X,E15.8,1X,E15.8)
         CLOSE(UNIT=NO4)
         ENDIF
C        -----------------------------------------------

C        -----------------------------------------------
C        LIUYQ, 2005-01-31
C        OUTPUT JACOBIAN IN FOURIER SPACE

         OPEN(UNIT=NO5,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='JACOB_F')
         WRITE(NO5,1021) MSMAX,NPSI1
         DO J=1,MSMAX
            DO I=1,NPSI1
               ZTEMP(1)=REAL(JACOBI(I,J))
               ZTEMP(2)=IMAG(JACOBI(I,J))
               WRITE(NO5,1022) ZTEMP(1),ZTEMP(2)
            ENDDO
         ENDDO
 1021    FORMAT(I5,1X,I5)
 1022    FORMAT(E15.8,1X,E15.8)
         CLOSE(UNIT=NO5)
C        -----------------------------------------------

         RETURN
C
C----------------------------------------------------------------------
C  28. WRITE QUANTITIES FOR PLOTS
C
 280     CONTINUE
C
         CALL WRTPLOT
C
         RETURN
C
C----------------------------------------------------------------------
C  29.   OPEN UNIT NETVAC
C
 290     CONTINUE
C
         OPEN(UNIT=NETVAC,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='ETAVAC')
         REWIND NETVAC
C
         RETURN
C
C----------------------------------------------------------------------
C  30. WRITE VACUUM QUANTITIES FOR LINEAR CODE
C
 300     CONTINUE
C
         WRITE(NETVAC,1008) NV1,R0W,RZ0W
         WRITE(NETVAC,1002) (CSV(I),I=1,NV1)
C
         CALL GENOUT(DG11LV(1,1),' VG11L',NETVAC,NPV1,NV1,MSMAX,1,RM,RN)
         CALL GENOUT(DG22LV(1,1),' VG22L',NETVAC,NPV1,NV1,MSMAX,1,RM,RN)
         CALL GENOUT(DG33LV(1,1),' VG33L',NETVAC,NPV1,NV1,MSMAX,1,RM,RN)
         CALL GENOUT(DG12LV(1,1),' VG12L',NETVAC,NPV1,NV1,MSMAX,1,RM,RN)
C
         CALL GENOUT(DG11LMV(1,1),'VG11LM',NETVAC,NPV,NV,MSMAX,1,RM,RN)
         CALL GENOUT(DG22LMV(1,1),'VG22LM',NETVAC,NPV,NV,MSMAX,1,RM,RN)
         CALL GENOUT(DG33LMV(1,1),'VG33LM',NETVAC,NPV,NV,MSMAX,1,RM,RN)
         CALL GENOUT(DG12LMV(1,1),'VG12LM',NETVAC,NPV,NV,MSMAX,1,RM,RN)
C       CALL GENOUT(VJACOB (1,1),'VJACOB',NETVAC,NPV1,NV1,MSMAX,1,RM,RN)
C        CALL GENOUT(VJACOM (1,1),'VJACOM',NETVAC,NPV,NV,MSMAX,1,RM,RN)
C
         RETURN
C
C----------------------------------------------------------------------
C  31.  CLOSE UNIT NETVAC
C
 310     CONTINUE
C
         CLOSE(UNIT=NETVAC,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  32.   SAVE QUANTITIES ON "JSOLVER"
C
 320     CONTINUE
C
         JSOLVER = 47
C
         OPEN(UNIT=JSOLVER,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='JSOLVER')
         REWIND JSOLVER
C
         IP = ISRCHFGE(NPSI1,PSIISO,1,CPSICL(1))
C
         IF (IP.LT.1)     IP = 1
         IF (IP.GT.NPSI1) IP = NPSI1
C
         NSTTP = 3
C
         DO 321 J321=IP,NPSI1
C
         CALL CINT(J321,SIGPSI(1,J321),TETPSI(1,J321),WGTPSI(1,J321))
C
 321     CONTINUE
C
         IF (IP .GT. 1) THEN
C
            DO 322 J322=1,IP-1
C
            CID0(J322) = FCCCC0(RMAG**2,CID0(IP),CID0(IP+1),
     ,                          CID0(IP+2),0.,CSM(IP),
     ,                          CSM(IP+1),CSM(IP+2),CSM(J322))
            CID2(J322) = FCCCC0(0.,CID2(IP),CID2(IP+1),
     ,                          CID2(IP+2),0.,CSM(IP),
     ,                          CSM(IP+1),CSM(IP+2),CSM(J322))
C
  322       CONTINUE
C
         ENDIF
C
         DO 323 J323=2,NPSI
C
         IS = ISRCHFGE(NPSI1,CSM,1,CS(J323))
C
         IF (IS .LE. 1)    IS = 2
         IF (IS .GE. NPSI) IS = NPSI - 1
C
         ZPSI(J323) = SPSIM * (1. - CS(J323)**2)
         ZTMF(J323) = FCCCC0(TMF(IS-1),TMF(IS),TMF(IS+1),TMF(IS+2),
     ,                       CSM(IS-1),CSM(IS),CSM(IS+1),CSM(IS+2),
     ,                       CS(J323))
         ZTTP(J323) = FCCCC0(TTP(IS-1),TTP(IS),TTP(IS+1),TTP(IS+2),
     ,                       CSM(IS-1),CSM(IS),CSM(IS+1),CSM(IS+2),
     ,                       CS(J323))
         ZCPR(J323) = FCCCC0(CPR(IS-1),CPR(IS),CPR(IS+1),CPR(IS+2),
     ,                       CSM(IS-1),CSM(IS),CSM(IS+1),CSM(IS+2),
     ,                       CS(J323))
         ZCPPR(J323) = FCCCC0(CPPR(IS-1),CPPR(IS),CPPR(IS+1),CPPR(IS+2),
     ,                        CSM(IS-1),CSM(IS),CSM(IS+1),CSM(IS+2),
     ,                        CS(J323))
         ZCID0(J323) = FCCCC0(CID0(IS-1),CID0(IS),CID0(IS+1),CID0(IS+2),
     ,                        CSM(IS-1),CSM(IS),CSM(IS+1),CSM(IS+2),
     ,                        CS(J323))
         ZCID2(J323) = FCCCC0(CID2(IS-1),CID2(IS),CID2(IS+1),CID2(IS+2),
     ,                        CSM(IS-1),CSM(IS),CSM(IS+1),CSM(IS+2),
     ,                        CS(J323))
C
 323     CONTINUE
C
         ZPSI(1)  = SPSIM
         ZTMF(1)  = T0
         ZTTP(1)  = DTTP0
         ZCPR(1)  = CP0
         ZCPPR(1) = DPDP0
         ZCID0(1) = RMAG**2
         ZCID2(1) = 0.
C
         ZPSI(NPSI1)  = 0.
         ZTMF(NPSI1)  = TMF(NPSI1)
         ZTTP(NPSI1)  = TTP(NPSI1)
         ZCPR(NPSI1)  = CPR(NPSI1)
         ZCPPR(NPSI1) = CPPR(NPSI1)
         ZCID0(NPSI1) = CID0(NPSI1)
         ZCID2(NPSI1) = CID2(NPSI1)
C
         IBND = 12 * NT
C
         IF (NSYM .EQ. 1) THEN
C
            ZDT  = CPI / DFLOAT(IBND-2)
C
            DO 325 J325=1,IBND+1
C
            ZT(J325) = CPI - (J325 - 2) * ZDT
C
 325        CONTINUE
C
         ELSE
C
            ZDT  = 2 * CPI / DFLOAT(IBND)
C
            DO 326 J326=1,IBND+1
C
            ZT(J326) = CPI - (J326 - 2) * ZDT
C
 326        CONTINUE
C
         ENDIF
C
         CALL BOUND(IBND+1,ZT,ZBND)
C
         DO 327 J327=1,IBND+1
C
         ZR(J327) = R0  + ZBND(J327) * COS(ZT(J327))
         ZZ(J327) = RZ0 + ZBND(J327) * SIN(ZT(J327))
C
 327     CONTINUE
C
         IMN = ISMIN(IBND+1,ZR,1)
         IMX = ISMAX(IBND+1,ZR,1)
C
         ZRMJ = .5 * (ZR(IMN) + ZR(IMX))
C
         WRITE(JSOLVER,1200) NPSI1, NSYM, IBND
         WRITE(JSOLVER,1201) RMAG, RZMAG, RITOT
         WRITE(JSOLVER,1201) BETA, BETAP, .5 * RINDUC(NPSI1)
         WRITE(JSOLVER,1201) ZPSI(1), ZPSI(NPSI1), ZRMJ
         WRITE(JSOLVER,1203) (ZPSI(L),L=1,NPSI1)
         WRITE(JSOLVER,1203) (ZCPR(L),L=1,NPSI1)
         WRITE(JSOLVER,1203) (ZCPPR(L),L=1,NPSI1)
         WRITE(JSOLVER,1203) (ZTMF(L),L=1,NPSI1)
         WRITE(JSOLVER,1203) (ZTTP(L),L=1,NPSI1)
         WRITE(JSOLVER,1203) (ZCID0(L),L=1,NPSI1)
         WRITE(JSOLVER,1203) (ZCID2(L),L=1,NPSI1)
         WRITE(JSOLVER,1203) (ZR(L),L=1,IBND+1)
         WRITE(JSOLVER,1203) (ZZ(L),L=1,IBND+1)
C
         CLOSE(UNIT=JSOLVER,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  33.   READ EXPERIMENTAL EQUILIBRIUM ON EXPEQ
C        IF (NSURF=6 .AND. NEQDSK=1 OR 2) EXPEQ IS IN EQDSK FORMAT (MKSA)
C        AND VARIABLES ARE NORMALIZED
C
 330     CONTINUE
C
         NXIN = 48
C
         OPEN(UNIT=NXIN,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='EXPEQ')
         REWIND NXIN
C
         IF (NSURF.EQ.6 .AND. (NEQDSK.EQ.1 .OR. NEQDSK.EQ.2)) GO TO 332
C
C-----------------------------------------------------------------------
CL       33.1 NORMAL EXPEQ FILE
C
         READ(NXIN,*) ASPCT
         READ(NXIN,*) RZ0C
         READ(NXIN,*) PREDGE
         IF (RZ0.EQ.0 .AND. RZ0C.NE.0.0) RZ0 = RZ0C
C
         IF (NSURF .EQ. 6) THEN
C
C NONCONFORMAL WALL
            READ(NXIN,*) NBPS,NWBPS,NDATA
C
            IF (NBPS .GT. NPBPS) THEN
               PRINT*,
     &         'NBPS LARGER THAN NPBPS.RECOMPILE WITH LARGER NPBPS'
               STOP
            ENDIF
C
            NWBPS = NWBPS + 1
            READ(NXIN,*) (RRBPS(L,1),  RZBPS(L,1), L=1,NBPS)
            IF (NDATA.EQ.2) THEN
               DO J=2,NWBPS-1
                  READ(NXIN,*) (RRBPS(L,J),  RZBPS(L,J), L=1,NBPS)
               ENDDO
            ELSE IF (NDATA.EQ.3) THEN
            DO J=2,NWBPS-1
            READ(NXIN,*) (RRBPS(L,J),RZBPS(L,J),CNDRZ(L,J),L=1,NBPS)
            ENDDO
            ENDIF

C NONCONFORMAL WALL: COMPUTE OUTER BOUNDARY (RRBPS(L,NWBPS),RZBPS(L,NWBPS))
C NOTE RESET ASPCT VALUE HERE TO AVOID INCONSISTENCY IN EXPEQ INPUT DATA
            IRMAX     = ISMAX(NBPS,RRBPS(1,1),1)
            IRMIN     = ISMIN(NBPS,RRBPS(1,1),1)
            RPC0      = (RRBPS(IRMIN,1)+RRBPS(IRMAX,1))*.5
            ASPCT     = (RRBPS(IRMAX,1)-RRBPS(IRMIN,1))*.5/RPC0
            TMP_RP    = RRBPS(IRMAX,1) - R0W 
            TMP_HV    = (REXT-1.0)/DFLOAT(NV)
            RW(1)     = 1.0
            RW(NWBPS) = REXT
            DO J=2,NWBPS-1
               IRMAX = ISMAX(NBPS,RRBPS(1,J),1)
               RW(J) = (RRBPS(IRMAX,J)-R0W)/TMP_RP
            ENDDO
            
            DO J=1,NWBPS
               WRITE(*,*) 'J=',J,' RW(J)=',RW(J)
            ENDDO

            DO J=2,NWBPS-1
               IF (RW(J).LE.RW(J-1).OR.RW(J).GE.RW(J+1)) THEN
                  WRITE(*,*) 'ERROR: THE WALL SHAPE IS TOO STRANGE!'
                  STOP 'WALL SHAPE WRONG'
               ENDIF
               DO L=1,NBPS
                  IF (RRBPS(L,J).LE.0.01) THEN
                     WRITE(*,*) 'WALL SHAPE: L,J,R,Z=',L,J,
     &                          RRBPS(L,J),RZBPS(L,J)
C                    STOP 'WALL SHAPE WRONG'
                  ENDIF
               ENDDO
            ENDDO

            TMP_B   = RW(NWBPS-1)
            TMP_REB = REXT/TMP_B
            J       = NWBPS
            DO L=1,NBPS
               RRBPS(L,J) = R0W + TMP_REB*(RRBPS(L,J-1)-R0W)
               RZBPS(L,J) = RZ0W + TMP_REB*(RZBPS(L,J-1)-RZ0W)
               IF (RRBPS(L,J).LE.0.01) THEN
                  WRITE(*,*) 'WALL SHAPE: L,J,R,Z=',L,J,
     &                       RRBPS(L,J),RZBPS(L,J)
C                 STOP 'WALL SHAPE WRONG'
               ENDIF
            ENDDO
         ENDIF
C
         IF (NSURF .EQ. 7) THEN
            READ(NXIN,*) NFOURPB
C     NEEDS VALUES OF RC AND RZ0C TO CORRECTLY RECONSTRUCT PLASMA BOUNDARY
            READ(NXIN,*) ALZERO, RC
            READ(NXIN,*) (BPSCOS(L),BPSSIN(L),L=1,NFOURPB)
         END IF
C
cab      IF (NFUNC.EQ.4 .AND. NPPFUN.EQ.4) THEN
         IF (.TRUE.) THEN
C        
C    BOTH PROFILES GIVEN AS AN ARRAY OF NPPF1 POINTS VS S
C
            READ(NXIN,*) NPPF1,NSTTP
C
            IF (NPPF1 .GT. NPISO) THEN
               PRINT*,
     &         'NPPF1 LARGER THAN NPISO.RECOMPILE WITH LARGER NPISO'
               STOP
            ENDIF
C
            READ(NXIN,*) (FCSM(L), L=1,NPPF1)
            READ(NXIN,*) (RPPF(L), L=1,NPPF1)
            READ(NXIN,*) (RFUN(L), L=1,NPPF1)
C
C============================================================================
C Y.Q. Liu, 2003-11-23
C Simple pressure profile scaling
C Works only if current I^* is given, and ballooning optimization is excluded
C Could be useful for finding beta_N limits
            
            IF (NBLOPT.EQ.0) THEN
               DO L=1,NPPF1
                  RPPF(L) = RPPF(L)*CFBAL
               ENDDO
            ENDIF
C============================================================================ 

            NPPF = NPPF1 - 1
C
          GOTO 333
          ELSE IF (NFUNC.EQ.1 .AND.NPPFUN.EQ.1) THEN
C        
C    BOTH PROFILES GIVEN AS POLYNOMIAL OF PSI/PSIAXIS OF DEGREE NSOUR-1
C
            READ(NXIN,*) NSOUR,NSTTP
C
            IF (NSOUR .GT. 10) THEN
              PRINT *,' NSOUR TOO LARGE IN EXPEQ'
              STOP
            ENDIF
C
            READ(NXIN,*) (AP(L),L=1,NSOUR)
            READ(NXIN,*) (AT(L),L=1,NSOUR)
C        
          ELSE
C        
C        MIXED PROFILES AS POLYNOMIAL AND ARRAY
C
C      P-PRIME PROFILE
            IF (NPPFUN .EQ. 4) THEN
              READ(NXIN,*) NPPF1
              READ(NXIN,*) NSTTP
              IF (NPPF1 .GT. NPISO) THEN
                PRINT*,
     &             'NPPF1 LARGER THAN NPISO.RECOMPILE WITH LARGER NPISO'
                STOP
              ENDIF
              READ(NXIN,*) (FCSM(L), L=1,NPPF1)
              READ(NXIN,*) (RPPF(L), L=1,NPPF1)
              NPPF = NPPF1 - 1
            ELSE IF (NPPFUN .EQ. 1) THEN
              READ(NXIN,*) NSTTP
              READ(NXIN,*) NSOUR
              READ(NXIN,*) (AP(L),L=1,NSOUR)
            ENDIF
C        
C
C      TT-PRIME OR I-PRIME OR .. PROFILE
            IF (NFUNC .EQ. 4) THEN
              READ(NXIN,*) NPPF1
              READ(NXIN,*) NSTTP
              IF (NPPF1 .GT. NPISO) THEN
                PRINT*,
     &             'NPPF1 LARGER THAN NPISO.RECOMPILE WITH LARGER NPISO'
                STOP
              ENDIF
              READ(NXIN,*) (FCSM(L), L=1,NPPF1)
              READ(NXIN,*) (RFUN(L), L=1,NPPF1)
              NPPF = NPPF1 - 1
            ELSE IF (NFUNC .EQ. 1) THEN
              READ(NXIN,*) NSTTP
              READ(NXIN,*) NSOUR
              READ(NXIN,*) (AT(L),L=1,NSOUR)
            ENDIF
         ENDIF
C
         GO TO 333
C
C-----------------------------------------------------------------------
CL       33.2 EQDSK TYPE OF EXPEQ
C
 332     CONTINUE
C
C        READ MKSA VALUES
C
         READ(NXIN,'(48X,3I4)') IDUMMY,INRBOX,INZBOX
         WRITE(*,*) 'INRBOX=',INRBOX,' INZBOX=',INZBOX
         IF (INRBOX .GT. NPISO) THEN
            PRINT *,' NRBOX LARGE THAN NPISO IN EQDSK, SHOULD .LE. ',
     +          NPISO
            STOP
         ENDIF
         IF (INRBOX*INZBOX .GT. NPISO*NPISO) THEN
            PRINT *,' NRBOX*NZBOX .GT. NPISO*NPISO: NRBOX= ',INRBOX,
     +          '  NZBOX= ',INZBOX,'  NPISO= ',NPISO
            STOP
         ENDIF
C
         READ(NXIN,9381) RBOXLEN,ZBOXLEN,R0EXP,RBOXLFT
         READ(NXIN,9381) ZRRAXIS,ZRZAXIS,ZPSIAX,ZDUMMY,B0EXP
         READ(NXIN,9381) CURRT
         READ(NXIN,9381) ZDUMMY
         NPPF1 = INRBOX
         NPPF = NPPF1 -1
C     G
         READ(NXIN,9381) (RFUN(L), L=1,NPPF1)
C     CHECK IF NEED TO RENORMALIZE TO BE CONSISTENT WITH NTMF0
C
         INTMF0 = 0
         IF (ABS(RFUN(1)-R0EXP*B0EXP) .LE. 1.E-04*RFUN(1)) THEN
           INTMF0 = 1
           WRITE(6,'(/,"   EQDSK HAS T(1) = R0*B0",/)')
         ENDIF
         IF (NEQDSK.EQ.2 .AND. INTMF0.NE.NTMF0) THEN
           IF (NTMF0 .EQ. 0) B0EXP = RFUN(NPPF1) / R0EXP
           IF (NTMF0 .EQ. 1) B0EXP = RFUN(1    ) / R0EXP
         ELSE IF (INTMF0 .NE. NTMF0) THEN
           WRITE(6,'(/,"  WARNING: EQDSK HAS NOT SAME NORMALIZATION AS",
     +       " ASK BY NTMF0: => EQUIL. IS DIFFERENT")')
         ENDIF
C
C        EQUIDISTANT MESH IN PSI
C
         ZDPSI = 1. / DFLOAT(NPPF1-1)
         DO I=1,NPPF1
           FCSM(I) = SQRT(DFLOAT(I-1)*ZDPSI)
         END DO
C     P
         READ(NXIN,9381) (RPPF(L), L=1,NPPF1)
         PREDGE = RPPF(NPPF1)
         IF (PREDGE .LT. 0.0) THEN
           WRITE(6,'(//,"   WARNING, PREDGE < 0 IN EXPEQ: PREDGE SET",
     +       " TO 0.0",/)')
           PREDGE = 0.0
         ENDIF
         RPPF = RPPF/(B0EXP**2/CPI/4.0E-7)
         WRITE(*,*) 'S      P'
         DO I=1,NPPF1
            WRITE(*,1302) FCSM(I),RPPF(I)
         ENDDO 
C     GG'
         READ(NXIN,9381) (RFUN(L), L=1,NPPF1)
C     P'
         READ(NXIN,9381) (RPPF(L), L=1,NPPF1)
         READ(NXIN,9381) ((EQDSPSI(I,J),I=1,INRBOX),J=1,INZBOX)
         READ(NXIN,9381) (QPSI(I),I=1,INRBOX)
         READ(NXIN,*) NBPS
         IF (NBPS .GT. NPBPS) THEN
            PRINT *,
     &          'NBPS LARGER THAN NPBPS.RECOMPILE WITH LARGER NPBPS'
            STOP
         ENDIF
C
C NONCONFORMAL WALL: ONLY PLASMA SURFACE CAN BE READ FROM EQDSK
         READ(NXIN,9381) (RRBPS(L,1),  RZBPS(L,1), L=1,NBPS)
         NWBPS = 2
         RW(1)     = 1.0
         RW(NWBPS) = REXT  
         
         DO J=1,NWBPS
            WRITE(*,*) 'J=',J,' RW(J)=',RW(J)
         ENDDO

         TMP_B = RW(NWBPS-1)
         TMP_REB = REXT/TMP_B
         DO J=1,NBPS
            RRBPS(J,NWBPS) = R0EXP + TMP_REB*(RRBPS(J,NWBPS-1)-R0EXP)
            RZBPS(J,NWBPS) = RZ0W + TMP_REB*(RZBPS(J,NWBPS-1)-RZ0W)
         ENDDO
C
C        AUXILIARY PARAMETERS
C        (COMMENT NEXT TWO LINES IF WANT TO CHANGE QSPEC IN NAMELIST)
         QSPEC = QPSI(1)
         CSSPEC = 0.0
C
C        EQUIDISTANT MESH IN PSI
C
         ZDPSI = 1. / DFLOAT(NPPF1-1)
         DO I=1,NPPF1
           FCSM(I) = SQRT(DFLOAT(I-1)*ZDPSI)
         END DO
C
C        IMPOSE B0EXP POSITIVE
C
         B0EXP = ABS(B0EXP)

C        PRINT OUT KEY INFORMATION FROM EQDSK DATA
         I = 1 + FLOOR(0.95*(NPPF1-1))
         QPSI95 = QPSI(I)+(QPSI(I+1)-QPSI(I))*(0.95-FCSM(I)**2)/ZDPSI
         WRITE(*,*) 'KEY EQUILIBRIUM DATA FROM EQDSK:'
         WRITE(*,*) 'EQDSK: R0EXP[m] =',R0EXP
         WRITE(*,*) 'EQDSK: B0EXP[T] =',B0EXP
         WRITE(*,*) 'EQDSK: I0EXP[A] =',CURRT
         WRITE(*,*) 'EQDSK: RAXIS[m] =',ZRRAXIS
         WRITE(*,*) 'EQDSK: ZAXIS[m] =',ZRZAXIS
         WRITE(*,*) 'EQDSK: QZERO    =',QPSI(1) 
         WRITE(*,*) 'EQDSK: QPSI95   =',QPSI95  
         WRITE(*,*) 'EQDSK: QEDGE    =',QPSI(NPPF1)
         WRITE(*,*) 'EQDSK: Q-PROFILE (PSI_P,Q):'
         DO I=1,NPPF1
            WRITE(*,1302) FCSM(I)**2,QPSI(I)
         ENDDO 
C
C        NORMALIZE USING R0EXP AND B0EXP
C
C        NOTE: R0 = R0EXP / R0EXP = 1. MAY NOT BE MIDDLE OF PLASMA
C
         ZR0 = 1.0
         ZCOFR = 1. / R0EXP
         DO J=1,NWBPS
            DO I=1,NBPS
               RRBPS(I,J) = RRBPS(I,J) * ZCOFR
               RZBPS(I,J) = RZBPS(I,J) * ZCOFR
            ENDDO
         ENDDO

C        SPECIAL TREATMENT HERE FOR A SPARC EQUILIBRIUM
         K_SPARC = 0
         IF (K_SPARC.EQ.1) THEN
            R0EXP_SAVE = R0EXP
            B0EXP_SAVE = B0EXP
            R0EXP = 1.0
            B0EXP = 22.49
         ENDIF

         ZMU0 = 4.E-07 * CPI
         ZSIGNCU = SIGN(1.0,CURRT)
         CURRT = ABS(CURRT*ZMU0/R0EXP/B0EXP)
         ZCOFTTP = ZSIGNCU / B0EXP 
         ZCOFPP = ZSIGNCU * ZMU0 * R0EXP * R0EXP / B0EXP
         DO I=1,NPPF1
           RFUN(I) = RFUN(I) * ZCOFTTP
           RPPF(I) = RPPF(I) * ZCOFPP * CFBAL
         END DO
         PREDGE = PREDGE * ZMU0 / B0EXP / B0EXP * CFBAL

C        SPECIAL TREATMENT HERE FOR A SPARC EQUILIBRIUM
         IF (K_SPARC.EQ.1) THEN
            R0EXP = R0EXP_SAVE
            B0EXP = B0EXP_SAVE
         ENDIF
C
C        FIND PLASMA LIMITS
C
         IRMAX = ISMAX(NBPS,RRBPS(1,1),1)
         IRMIN = ISMIN(NBPS,RRBPS(1,1),1)
         IZMAX = ISMAX(NBPS,RZBPS(1,1),1)
         IZMIN = ISMIN(NBPS,RZBPS(1,1),1)
         ZR0 = 0.5 * (RRBPS(IRMIN,1) + RRBPS(IRMAX,1))
         RZ0 = 0.5 * (RZBPS(IZMIN,1) + RZBPS(IZMAX,1))
         ASPCT = (RRBPS(IRMAX,1) - RRBPS(IRMIN,1))
     +     / (RRBPS(IRMAX,1) + RRBPS(IRMIN,1))
         WRITE(6,'(/,3(A,1PE13.4),/)') ' FROM EQDSK: ,ZR0= ',ZR0,
     +     ' RZ0= ',RZ0,' ASPCT= ',ASPCT

C        WRITE EXPEQG.OUT
C        YQL, 2013-03-15
         NXOUT = 49
         OPEN(UNIT=NXOUT,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='EXPEQG.OUT')
         REWIND NXOUT
         WRITE(NXOUT,1303) ASPCT,RZ0,PREDGE
         WRITE(NXOUT,1200) NBPS,1,2
         WRITE(NXOUT,1302) (RRBPS(L,1),RZBPS(L,1), L=1,NBPS)
         WRITE(NXOUT,1301) NPPF1,1
         WRITE(NXOUT,1303) (FCSM(L), L=1,NPPF1)
         WRITE(NXOUT,1303) (RPPF(L),L=1,NPPF1)
         WRITE(NXOUT,1303) (RFUN(L),L=1,NPPF1)
         CLOSE(UNIT=NXOUT)
C
C-----------------------------------------------------------------------
C
 333     CONTINUE
         IF (IDIIICOIL.EQ.1) THEN
            READ (NXIN,'("RMJ,COILWIDTHN=",1p2e12.4)')
     &           RMJ,WIDTHCOILN
            READ (NXIN,'("RCCOILN,ZCCOILN,WCCOILN=",1p3e12.4)')
     &        RCCOILN,ZCCOILN,WCCOILN
            READ (NXIN,'("RICOILN,ZICOILN,WICOILN=",1p3e12.4)')
     &        RICOILN,ZICOILN,WICOILN
         END IF
C
         CLOSE(UNIT=NXIN,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  34.   WRITE EQUIL. IN "EXPERIMENTAL" FORMAT ON EXPEQ.OUT
C
 340     CONTINUE
C
         NXOUT = 50
C
         OPEN(UNIT=NXOUT,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     O        FILE='EXPEQ.OUT')
         REWIND NXOUT
C
         WRITE(NXOUT,1303) ASPCT,RZ0,PREDGE
C
         IF (NSURF .EQ. 7) THEN
            WRITE(NXOUT,1301) NFOURPB
            WRITE(NXOUT,1302) ALZERO, RC
            WRITE(NXOUT,1302) (BPSCOS(L),BPSSIN(L),L=1,NFOURPB)
         ELSE
            WRITE(NXOUT,1200) NBPS,NWBPS-1,NDATA
            WRITE(NXOUT,1302) (RRBPSNEW(L),RZBPSNEW(L), L=1,NBPS)
            IF (NDATA.EQ.2) THEN
               DO J=2,NWBPS-1
                  WRITE(NXOUT,1302) (RRBPS(L,J),RZBPS(L,J),L=1,NBPS)
               ENDDO
            ELSE IF (NDATA.EQ.3) THEN
               DO J=2,NWBPS-1
                  WRITE(NXOUT,1304) (RRBPS(L,J),RZBPS(L,J),
     &                               CNDRZ(L,J),L=1,NBPS)
               ENDDO
            END IF 
         END IF
 
C        WRITE(NXOUT,1301) NPSI1,NPROPT
         WRITE(NXOUT,1301) KSHAVE,NPROPT
         WRITE(NXOUT,1303) (CSM(L)/CSM(KSHAVE), L=1,KSHAVE)
C
         WRITE(NXOUT,1303) (CPPR(L),L=1,KSHAVE)
         IF (NPROPT.NE.2 .AND. NPROPT.NE.3) WRITE(NXOUT,1303) 
     +                                        (TTP(L),L=1,KSHAVE)
         IF (NPROPT .EQ. 2) WRITE(NXOUT,1303) (RIPR(L),L=1,KSHAVE)
         IF (NPROPT .EQ. 3) WRITE(NXOUT,1303) (RJDOTB(L),L=1,KSHAVE)
C
CL       ADD, IN MKSA, EXTRA VALUES OBTAINED BY CHEASE FOR EASIER COMPARISON
C        (NO HEADER)
C        
         CALL OUTMKSA(NXOUT,1)
C
         CLOSE(UNIT=NXOUT,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C  35.   SAVE QUANTITIES FOR PENN CODE ON "PENN"
C
 350     CONTINUE
C
         OPEN(UNIT=NPENN,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     O        FILE='NPENN')
         REWIND NPENN
C
         IDCHSE = 222
         IDATA  = 13
C
         WRITE(NPENN) IDCHSE
         WRITE(NPENN) IDATA
C
         RETURN
C
C----------------------------------------------------------------------
C  36.   CLOSE FILE "PENN"
C
 360     CONTINUE
C
         CLOSE(UNIT=NPENN,STATUS='KEEP')
C
         RETURN
C
C----------------------------------------------------------------------
C
C  37. SAVE QUANTITIES FOR XTOR INTO "OUTXTOR"
C
 370     CONTINUE
C
         OPEN(NXTOR,FILE='OUTXTOR',STATUS='NEW',FORM='UNFORMATTED')
         REWIND NXTOR
C
         WRITE(NXTOR) NPSI1, NTNOVA
         WRITE(NXTOR) RC / ASPCT
C
         RETURN
C
C----------------------------------------------------------------------
C
C        38. SAVE QUANTITIES ON EQDSK IN MKSA USING R0EXP AND B0EXP
C        IF NEQDSK = -2: CHANGE B0EXP SO THAT T(EDGE) = R0EXP*BOEXP
 380     CONTINUE
C
         IEQDSK = 38
         OPEN(IEQDSK,FILE='EQDSK.OUT',FORM='FORMATTED')
C     
         IF (NEQDSK.EQ.-2 .AND. NTMF0.EQ.1) B0EXP = B0EXP / TMF(NPSI1)
C
C        FIRST LINE: 48 CHAR, DUMMY INTEGER, NRBOX, NZBOX
C
         ZDATE='20050328'
C+DATE   CALL DATE(ZDATE)
         IDUMMY = 3
CJEM         WRITE(IEQDSK,9380) 'FROM CHEASE, ALL IN MKSA UNITS     ',
         WRITE(IEQDSK,9380) ' CHEASE    00/00/00      #000000,00000ms',
     +       ZDATE,IDUMMY,NRBOX,NZBOX
C
C        2ND LINE: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZBOXMID
C
         WRITE(IEQDSK,9381) RBOXLEN,ZBOXLEN,R0EXP,RBOXLFT
CJEM     +       ,RC0P
     +     ,RZMAG*R0EXP
C
C        3RD LINE: RMAG, ZMAG, PSIMAG, PSIEDGE, B0
C        ZMAG HAS BEEN SHIFTED TO ZERO IN PSIBOX
CYQL     KEEP ZMAG UNSHIFTED

         ZMAG0 = RZMAG
         ZPSIAX = SPSIM*B0EXP*R0EXP**2
         WRITE(IEQDSK,9381) RRAXIS*R0EXP,ZMAG0*R0EXP,
     +       ZPSIAX,RC0P,B0EXP
C
C        4TH LINE: PLASMA CURRENT, PSIAX1, PSIAX2, RAXIS1, RAXIS2
C
         ZMU0 = 4.E-07 * CPI
         IF (R0EXP.EQ.1 .AND. B0EXP.EQ.1) ZMU0 = 1.0
         WRITE(IEQDSK,9381) RITOT*R0EXP*B0EXP/ZMU0,
     +       ZPSIAX,RC0P,RRAXIS*R0EXP,RC0P
C
C        5TH LINE: ZAXIS1, ZAXIS2, PSI_SEP, R_XPOINT, Z_XPOINT
C
         WRITE(IEQDSK,9381) ZMAG0*R0EXP,RC0P,RC0P,RC0P,RC0P
C
C     EQUISTANT PSI-MESH FOR PROFILES IN S=SQRT(1.-PSI/PSIMIN)
C
         ZDPSI = ABS(SPSIM) / DFLOAT(NRBOX-1)
         DO I=1,NRBOX
           ZSTEMP(I) = SPSIM + DFLOAT(I-1)*ZDPSI
         ENDDO
C
C        6TH ENTRY: T(PSI) (OR G)
C
         CALL SPLINE(CSM,TMF,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,TMF,ZD2TMP,ZTEMP,0)
         ZTEMP(1) = T0
         ZCOF = R0EXP*B0EXP
         WRITE(IEQDSK,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)
C
C        7TH ENTRY: PRESSURE
C
         CALL SPLINE(CSM,CPR,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,CPR,ZD2TMP,ZTEMP,0)
         ZTEMP(1) = CP0
         ZCOF = B0EXP*B0EXP / ZMU0
         WRITE(IEQDSK,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)
C
C        8TH ENTRY: TT' (OR GG')
C
         CALL SPLINE(CSM,TTP,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,TTP,ZD2TMP,ZTEMP,0)
         ZTEMP(1) = DTTP0
         ZCOF = B0EXP
         WRITE(IEQDSK,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)
C
C        9TH ENTRY: P'
C
         CALL SPLINE(CSM,CPPR,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,CPPR,ZD2TMP,ZTEMP,0)
         ZTEMP(1) = DPDP0
         ZCOF = B0EXP / ZMU0 / R0EXP / R0EXP
         WRITE(IEQDSK,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)
C
C        10TH ENTRY: PSI(I,J)
C
         ZCOF = B0EXP * R0EXP**2
         WRITE(IEQDSK,9381) ((EQDSPSI(I,J)*ZCOF,I=1,NRBOX),
     +       J=1,NZBOX)
C
C        11TH ENTRY: Q PROFILE
C
         CALL SPLINE(CSM,QPSI,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,QPSI,ZD2TMP,ZTEMP,0)
         ZTEMP(1) = Q0
         WRITE(IEQDSK,9381) (ZTEMP(I),I=1,NRBOX)
C
C        12TH ENTRY: (R,Z) OF PLASMA BOUNDARY AND DUMMY LIMITER POSITION
C
C     LIMITER ARRAY
CJEM         ILIMITER = 4
         ILIMITER = 5
         ZR(1) = RBOXLFT
CJEM         ZZ(1) = -ZBOXLEN/2.
         ZZ(1) = -ZBOXLEN/2.+ZMAG0*R0EXP
         ZR(2) = RBOXLFT+RBOXLEN
CJEM         ZZ(2) = -ZBOXLEN/2.
         ZZ(2) = -ZBOXLEN/2.+ZMAG0*R0EXP
CJEM         ZR(3) = RBOXLFT
         ZR(3) = RBOXLFT+RBOXLEN
CJEM         ZZ(3) =  ZBOXLEN/2.
         ZZ(3) =  ZBOXLEN/2.+ZMAG0*R0EXP
CJEM         ZR(4) = RBOXLFT+RBOXLEN
         ZR(4) = RBOXLFT
CJEM         ZZ(4) =  ZBOXLEN/2.
         ZZ(4) =  ZBOXLEN/2.+ZMAG0*R0EXP
         ZR(5) = RBOXLFT
         ZZ(5) = -ZBOXLEN/2.+ZMAG0*R0EXP
         WRITE(IEQDSK,1204) NBPSOUT, ILIMITER
         ZCOF = R0EXP
         WRITE(IEQDSK,9381) (RRBPSOU(I)*ZCOF,(RZBPSOU(I))*ZCOF,
     +     I=1,NBPSOUT)
         WRITE(IEQDSK,9381) (ZR(I),ZZ(I),I=1,ILIMITER)
C
C        LAST LINES: SOME EQUILIBRIUM VALUES
C        (NO HEADER)
C
         CALL OUTMKSA(IEQDSK,1)
C
         CLOSE(UNIT=IEQDSK,STATUS='KEEP')
C
         RETURN
C
 1000    FORMAT(A80)
 1001    FORMAT(3I20,2F20.5)
 1002    FORMAT(//,' RADIAL (INTEGER) MESH ',/,(2D30.15))
 1003    FORMAT(1X,13I8)
 1004    FORMAT((1X,40(I3)))
 1005    FORMAT(A132)
 1006    FORMAT((1X,8(1PE15.6)))
 1007    FORMAT(13I8)
 1008    FORMAT(I20,2F20.5)
 1014    FORMAT(5(1X,E12.5))
 1015    FORMAT(I8)
 1016    FORMAT(E12.5)
 1111    FORMAT(A8)
 1113    FORMAT(A24)
 1200    FORMAT(3I5)
 1201    FORMAT(3E15.8)
 1202    FORMAT(2E15.8)
 1203    FORMAT(E15.8)
 1204    FORMAT(5I5)
 1301    FORMAT(I5)
 1302    FORMAT(2E18.8)
 1303    FORMAT(E18.8)
 1304    FORMAT(3E18.8)
 9380    FORMAT(A40,A8,3I4)
 9381    FORMAT(1P,5E16.9)
 2701    FORMAT(2I20,F20.5)
 2702    FORMAT(//,' Q-PROFILE ',/,(2D30.15))
 2703    FORMAT(//,' PRESSURE PROFILE ',/,(2D30.15))
 2704    FORMAT(//,' A-PROFILE ',/,(2D30.15))
 2705    FORMAT(//,' F-PROFILE ',/,(2D30.15))
 2706    FORMAT(//,' PSI_PRIME-PROFILE ',/,(2D30.15))
 2707    FORMAT(//,(1D30.15))
C
         END
C*DECK C3SB02
C*CALL PROCESS
         SUBROUTINE WRTPLOT
C        ##################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SB02  WRITE EQUILIBRIUM PLOT QUANTITIES INTO FILE NUPLO           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBAL.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMERA.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMINT.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMLAB.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMPLO.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMVAC.inc'
         INCLUDE 'COMDAT.inc'
C
         CHARACTER*180 TEXT(1:130)
         CHARACTER*8   CDATE,   CCLOK
         CHARACTER*24  CDATESUN
C
         DIMENSION
     I   IBALL(NPPSI1),     IMERCI(NPPSI1),    IMERCR(NPPSI1)
         DIMENSION
     R   ZABIC(NPCHI1),     ZABIPR(NPISO+1),   ZABIS(NPPSI1+1),       
     R   ZABISG(NSP1),      ZABIT(NTP1),       ZABR (2*NPPSI1+1),
     R   ZABS (NPPSI1),     ZABSM(NPPSI1+1),   ZCHI(NPCHI1),      
     R   ZCSIPR(NPISO+1),   
     R   ZOARS(NPPSI1),     ZOART(NTP1),       ZOBETS(NPPSI1),    
     R   ZOSHR(2*NPPSI1+1), ZOSHS(NPPSI1+1),   ZOJBR(2*NPPSI1+1), 
     R   ZOJBS(NPPSI1+1),   ZOJPR(2*NPPSI1+1), ZOJPS(NPPSI1+1), 
     R   ZOJBSR(2*NPPSI1+1,2),ZOJBSS(NPPSI1+1,2),  ZODIS(NPPSI1+1), 
     R   ZODRS(NPPSI1+1),   ZOBETR(2*NPPSI1+1),ZOTRR(2*NPPSI1+1),
     R   ZOTRS(NPPSI1+1),   ZODQR(2*NPPSI1+1),
     R   ZOFR (2*NPPSI1+1), ZODQS(NPPSI1),     ZOHS (NPPSI1+1),   
     R   ZOIPR(2*NPPSI1+1), ZOIPS(NPPSI1+1),   ZOJR (2*NPPSI1+1), 
     R   ZOPPR(2*NPPSI1+1), ZOPPS(NPPSI1+1),   ZOPR (2*NPPSI1+1), 
     R   ZOPS (NPPSI1+1),   ZOQR (2*NPPSI1+1), ZOQS (NPPSI1+1),   
     R   ZOTR(2*NPPSI1+1),  ZOTS(NPPSI1+1),    ZOTTR(2*NPPSI1+1), 
     R   ZOTTS(NPPSI1+1),   ZPAR(2*NPPSI1),
     R   ZR(12*NPT+1),      ZRCHI(25,NPPSI1),  ZRCURV(4*NPPSI1),  
     R   ZRHOS(NTP1),       
     R   ZRSUR(6*NPT),      ZRTET(NPT),        ZSIG1(NPPSI1),
     R   ZTET(NTP1),        ZTET1(NPPSI1),     ZTSUR(6*NPT),
     R   ZZ(12*NPT+1),      ZZCHI(25,NPPSI1),  ZZCURV(4*NPPSI1),  
     R   ZZTET(NPT), 
     R   ZZSUR(6*NPT)
C
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C
C----------------------------------------------------------------------
C  28. WRITE QUANTITIES FOR PLOTS
C
 280     CONTINUE
C
         REWIND NUPLO
C
         ZMU0 = 1.256
         ZBPERC = 100. * BETA
         ZBSPER = 100. * BETAS
         ZBXPER = 100. * BETAX
         ZINORM = RINOR/ZMU0
         ZIBSNO = RIBSNOR/ZMU0
         ZGM    = ZBPERC/ZINORM
         ZGMSTA = ZBSPER/ZINORM
         ZGMX   = ZBXPER/ZINORM
         ZBPOL1 = BETAP*CONVF
         ZLI1   = RINDUC(NPSI1)*CONVF
         ZBSF   = RITBS/RITOT
         ZBSFC  = RITBSC/RITOT
         ZCBS   = 0.
         IF (BETA.GT.0.) ZCBS1 = ZIBSNO*ZINORM/(ZBXPER*SQRT(ASPCT))
         IF (BETA.GT.0.) ZCBS2 = ZBSF/(SQRT(ASPCT)*ZBPOL1)
         ZBBS  = ZBXPER * ZBSF
C
         TEXT(1) = 'CHEASE - V12.95'
C
         WRITE (TEXT(2),1111) '        '
         WRITE (TEXT(3),1111) '        '
C
C+DATE   IF (COMPTYP .EQ. 'CRAY') THEN
C+DATE     CALL CLOCK(CCLOK)
C+DATE     CALL DATE(CDATE)
C+DATE     WRITE (TEXT(2),1111) CDATE
C+DATE     WRITE (TEXT(3),1111) CCLOK
C
C+DATE   ELSE IF (COMPTYP .EQ. 'SUN' .OR. COMPTYP .EQ. 'SG') THEN
C+DATE     CALL FDATE(CDATE)
C+DATE     WRITE (TEXT(2),1111) CDATESUN
C+DATE     WRITE (TEXT(3),1111) '        '
C+DATE   ENDIF
C
         CALL WHTEXT(LABEL1,TEXT(4))
         CALL WHTEXT(LABEL2,TEXT(5))
         CALL WHTEXT(LABEL3,TEXT(6))
         CALL WHTEXT(LABEL4,TEXT(7))
C
         CALL WHTEXT('PLASMA SURFACE :',TEXT(8))
         CALL WITEXT(NSURF, 'NSURF', TEXT(9),1)
         CALL WRTEXT(ASPCT, 'ASPCT', TEXT(10),1,1)
         CALL WRTEXT(ELONG, 'ELONG', TEXT(11),1,1)
         CALL WRTEXT(TRIANG,'TRIANG',TEXT(12),1,1)
         CALL WRTEXT(DELTA, 'DELTA', TEXT(13),1,1)
         CALL WRTEXT(THETA0,'THETA0',TEXT(14),1,1)
         CALL WRTEXT(BEANS, 'BEANS', TEXT(15),1,1)
         CALL WRTEXT(CETA,  'CETA',  TEXT(16),1,1)
         CALL WRTEXT(SGMA,  'SGMA',  TEXT(17),1,1)
         CALL WRTEXT(TRIPLT,'TRIPLT',TEXT(18),1,1)
         CALL WRTEXT(RNU,   'RNU',   TEXT(19),1,1)
         CALL WRTEXT(XI,    'XI',    TEXT(20),1,1)
         CALL WRTEXT(AREA,  'AREA',  TEXT(21),1,1)
         CALL WHTEXT('EQUILIBRIUM SOLUTION :',TEXT(22))
         CALL WRTEXT(RMAG,       'RMAG',  TEXT(23),1,1)
         CALL WRTEXT(RZMAG,      'ZMAG',  TEXT(24),1,1)
         CALL WRTEXT(PSI0,       'PSIMIN',TEXT(25),1,1)
         CALL WRTEXT(PSISCL,     'PSISCL',TEXT(26),1,1)
         CALL WRTEXT(Q0,         'Q0',    TEXT(27),1,1)
         CALL WRTEXT(QPSI(NPSI1),'QSURF', TEXT(28),1,1)
         CALL WRTEXT(QCYL,       'QCYL',  TEXT(29),1,1)
         CALL WRTEXT(T0,         'T0',    TEXT(30),1,1)
         CALL WRTEXT(TMF(NPSI1), 'TSURF', TEXT(31),1,1)
         CALL WRTEXT(RITOT, 'TOT. CUR.',             TEXT(32),1,1)
         CALL WRTEXT(RINOR, 'NORM. CUR.',            TEXT(33),1,1)
         CALL WRTEXT(ZINORM,'IN (MA,T,M)',           TEXT(34),1,1)
         CALL WRTEXT(ZBSF,  'I-B.S.(0)/I-TOT',       TEXT(35),1,1)
         CALL WRTEXT(ZBSFC, 'I-B.S.(NUE*)/I-TOT',    TEXT(36),1,1)
         CALL WRTEXT(RZION, 'ION CHARGE',            TEXT(37),1,1)
         CALL WRTEXT(ETAEI, 'D(LOG(T))/D(LOG(N))',   TEXT(38),1,1)
         CALL WRTEXT(ZCBS1, 'CBS1=IB.S./(G*SQRT(E))',TEXT(39),1,1)
         CALL WRTEXT(ZCBS2, 'CBS2=F/(BP(1)*SQR(E))', TEXT(40),1,1)
         CALL WRTEXT(ZBBS,  'BBS',                   TEXT(41),1,1)
         CALL WRTEXT(CONVF, 'CONV. FACT.',           TEXT(42),1,1)
         CALL WRTEXT(RINDUC(NPSI1),'LI',              TEXT(43),1,1)
         CALL WRTEXT(ZLI1,         'LI (G.A.)',       TEXT(44),1,1)
         CALL WRTEXT(ZBPERC,       'BETA',            TEXT(45),2,1)
         CALL WRTEXT(ZBSPER,       'BETA*',           TEXT(46),2,1)
         CALL WRTEXT(ZBXPER,       'BETA EXP.',       TEXT(47),2,1)
         CALL WRTEXT(ZGM,          'G (MA,T,M).',     TEXT(48),1,1)
         CALL WRTEXT(ZGMSTA,       'G* (MA,T,M).',    TEXT(49),1,1)
         CALL WRTEXT(ZGMX,         'G EXP. (MA,T,M).',TEXT(50),1,1)
         CALL WRTEXT(BETAP,        'BETA POL.',       TEXT(51),1,1)
         CALL WRTEXT(ZBPOL1,       'BETA POL. (G.A.)',TEXT(52),1,1)
C
         CALL WHTEXT('PROFILES :',TEXT(53))
         CALL WITEXT(NFUNC, 'NFUNC', TEXT(54),1)
         CALL WITEXT(NSTTP, 'NSTTP', TEXT(55),1)
         CALL WITEXT(NIPR,  'NIPR',  TEXT(56),1)
         CALL WITEXT(NSOUR, 'NSOUR', TEXT(57),1)
         CALL WITEXT(NPPFUN,'NPPFUN',TEXT(58),1)
         CALL WITEXT(NPP,   'NPP',   TEXT(59),1)
         CALL WRTEXT(PREDGE,'PREDGE',TEXT(60),1,1)
         CALL WITEXT(NBSOPT,'NBSOPT',TEXT(61),1)
         CALL WITEXT(NBSFUN,'NBSFUN',TEXT(62),1)
         CALL WITEXT(NBSTRP,'NBSTRP',TEXT(63),1)
         CALL WITEXT(NBLOPT,'NBLOPT',TEXT(64),1)
         CALL WITEXT(NBLC0, 'NBLC0', TEXT(65),2)
         CALL WITEXT(NTURN, 'NTURN', TEXT(66),2)
         CALL WITEXT(NPPR  ,'NPPR'  ,TEXT(67),2)
         CALL WRTEXT(CFBAL, 'CFBAL', TEXT(68),1,1)
         CALL WRTEXT(CPRESS,'CPRESS',TEXT(69),1,1)
         CALL WRTEXT(BSFRAC,'BSFRAC',TEXT(70),1,1)
         CALL WRTEXT(AT     ,'AT',   TEXT(71),3,13)
         CALL WRTEXT(AT2    ,'AT2',  TEXT(72),3,13)
         CALL WRTEXT(AT3    ,'AT3',  TEXT(73),3,13)
         CALL WRTEXT(AT4    ,'AT4',  TEXT(74),3,13)
         CALL WRTEXT(AP     ,'AP',   TEXT(75),3,13)
         CALL WRTEXT(AP2    ,'AP2',  TEXT(76),3,13)
         CALL WRTEXT(AFBS   ,'AFBS', TEXT(77),3,13)
         CALL WRTEXT(AFBS2  ,'AFBS2',TEXT(78),3,13)
         CALL WHTEXT('MESHES :',TEXT(79))
         CALL WRTEXT(R0,    'R0',    TEXT(80),1,1)
         CALL WRTEXT(RZ0,   'Z0',    TEXT(81),1,1)
         CALL WRTEXT(EPSLON,'EPSLON',TEXT(82),1,1)
         CALL WITEXT(NS,    'NS',    TEXT(83),2)
         CALL WITEXT(NT,    'NT',    TEXT(84),2)
         CALL WITEXT(NISO,  'NISO',  TEXT(85),2)
         CALL WITEXT(NTNOVA,'NTNOVA',TEXT(86),2)
         CALL WITEXT(NPSI,  'NPSI',  TEXT(87),2)
         CALL WITEXT(NCHI,  'NCHI',  TEXT(88),2)
         CALL WITEXT(NER,   'NER',   TEXT(89),2)
         CALL WITEXT(NEGP,  'NEGP',  TEXT(90),2)
         CALL WHTEXT('S-PACKING :',TEXT(91))
         CALL WITEXT(NMESHA,'NMESHA',TEXT(92),1)
         CALL WRTEXT(SOLPDA,'SOLPDA',TEXT(93),1,1)
         CALL WITEXT(NPOIDA,'NPOIDA',TEXT(94),3)
         CALL WITEXT(NDIFPS,'NDIFPS',TEXT(95),1)
         CALL WHTEXT('I*-PACKING :',TEXT(96))
         CALL WITEXT(NMESHB,'NMESHB',TEXT(97),1)
         CALL WRTEXT(SOLPDB,'SOLPDB',TEXT(98),1,1)
         CALL WITEXT(NPOIDB,'NPOIDB',TEXT(99),3)
         CALL WHTEXT('SIGMA-PACKING :',TEXT(100))
         CALL WITEXT(NMESHC,'NMESHC',TEXT(101),1)
         CALL WRTEXT(SOLPDC,'SOLPDC',TEXT(102),1,1)
         CALL WITEXT(NPOIDC,'NPOIDC',TEXT(103),3)
         CALL WHTEXT('THETA-PACKING :',TEXT(104))
         CALL WITEXT(NMESHD,'NMESHD',TEXT(105),1)
         CALL WRTEXT(SOLPDD,'SOLPDD',TEXT(106),1,1)
         CALL WITEXT(NPOIDD,'NPOIDD',TEXT(107),3)
         CALL WITEXT(NDIFT, 'NDIFT', TEXT(108),1)
         CALL WHTEXT('CHI-PACKING :',TEXT(109))
         CALL WITEXT(NMESHE,'NMESHE',TEXT(110),1)
         CALL WRTEXT(SOLPDE,'SOLPDE',TEXT(111),1,1)
         CALL WITEXT(NPOIDE,'NPOIDE',TEXT(112),3)
         CALL WHTEXT('NORMALIZATION :',TEXT(113))
         CALL WITEXT(NCSCAL,'NCSCAL',TEXT(114),1)
         CALL WITEXT(NTMF0, 'NTMF0', TEXT(115),1)
         CALL WITEXT(NRSCAL,'NRSCAL',TEXT(116),1)
         CALL WRTEXT(SCALE, 'SCALE', TEXT(117),1,1)
         CALL WRTEXT(CSSPEC,'CSSPEC',TEXT(118),1,1)
         CALL WRTEXT(QSPEC, 'QSPEC', TEXT(119),1,1)
         CALL WRTEXT(CURRT, 'CURRT', TEXT(120),1,1)
C
C  COMPUTE MAXIMUM AND MINIMUM OF R AND Z
C  REFERENCE IS THE MAGNETIC AXIS
C
         ZOART(1) = 0.
         ZABIT(1) = 0.
C
         CALL BOUND(NT1,CT,ZRHOS)
C
         DO 281 J281=1,NT
C
         ZR(J281)      = ZRHOS(J281) * COS(CT(J281))
         ZZ(J281)      = ZRHOS(J281) * SIN(CT(J281))
         ZRTET(J281)   = ZR(J281)
         ZZTET(J281)   = ZZ(J281)
         ZTET(J281)    = CT(J281) - CT(1)
         ZOART(J281+1) = ZOART(J281) + .25 * (CT(J281+1) - CT(J281)) *
     *                   (ZRHOS(J281)**2 + ZRHOS(J281+1)**2)
         ZABIT(J281+1) = DFLOAT(J281) / DFLOAT(NT)
C
  281    CONTINUE
C
         ZTET(NT1) = 2. * CPI
C
         ZRMAX = ZR(ISMAX(NT,ZR,1))
         ZRMIN = ZR(ISMIN(NT,ZR,1))
         ZZMAX = ZZ(ISMAX(NT,ZZ,1))
         ZZMIN = ZZ(ISMIN(NT,ZZ,1))
C
C  COMPUTE THE SURFACE
C
         BPS( 1) = RMAG
         BPS(12) = RZMAG
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
         ZDT = 2. * CPI / (6 * NT - 1)
C
         DO 282 J282=1,6*NT
C
         ZTSUR(J282) = (J282 - 1) * ZDT
C
  282    CONTINUE
C
         CALL BOUND(6*NT,ZTSUR,ZRSUR)
         CALL BOUND(1,PANGLE    ,ZBND1)
         CALL BOUND(1,PANGLE+CPI,ZBND2)
C
         BPS( 1) = R0
         BPS(12) = RZ0
C
         IF (NSURF .NE. 1) BPS(3) = BPS(2) - BPS(1)
         IF (NSURF .EQ. 6) CALL BNDSPL
C
         DO 283 J283=1,6*NT
C
         ZZSUR(J283) = ZRSUR(J283) * SIN(ZTSUR(J283))
         ZRSUR(J283) = ZRSUR(J283) * COS(ZTSUR(J283))
C
  283    CONTINUE
C
C  COMPUTE THE VALUES OF R ON RADIUS USED FOR PROFILE DEFINITION
C
         CALL RMRAD(NPSI1,0.,CPSRF,PANGLE,ZPAR(1),ZSIG1,ZTET1,1)
         CALL RMRAD(NPSI1,0.,CPSRF,PANGLE+CPI,ZPAR(NPSI1+1),ZSIG1,
     ,              ZTET1,1)
C
         DO 284 J284=1,NPSI1
C
         I1 = NPSI1 + J284 + 1
         I2 = NPSI1 - J284 + 1
C
         ZABR(I1)  =   ZPAR(J284)       * ZBND1
         ZABR(I2)  = - ZPAR(NPSI1+J284) * ZBND2
         ZOQR(I2)  = QPSI(J284)
         ZOQR(I1)  = QPSI(J284)
         ZOPPR(I2) = CPPR(J284)
         ZOPPR(I1) = CPPR(J284)
         ZOPR(I2)  = CPR(J284)
         ZOPR(I1)  = CPR(J284)
         ZOTTR(I2) = TTP(J284)
         ZOTTR(I1) = TTP(J284)
         ZOTR(I2)  = TMF(J284)
         ZOTR(I1)  = TMF(J284)
         ZOIPR(I2) = RIPR(J284)
         ZOIPR(I1) = RIPR(J284)
         ZOJBR(I2) = RJDOTB(J284)
         ZOJBR(I1) = RJDOTB(J284)
         ZODQR(I2) = CDQ(J284)
         ZODQR(I1) = CDQ(J284)
         ZOSHR(I2) = CDRQ(J284)
         ZOSHR(I1) = CDRQ(J284)
         ZOFR(I2)  = CSM(J284)**2 * CPSRF
         ZOFR(I1)  = CSM(J284)**2 * CPSRF
         ZOJR(I1)  = - (ZABR(I1) + RMAG) * ZOPPR(I1) - 
     -                  ZOTTR(I1) / (ZABR(I1) + RMAG)
         ZOJR(I2)  = - (ZABR(I2) + RMAG) * ZOPPR(I2) - 
     -                  ZOTTR(I2) / (ZABR(I2) + RMAG)
         ZOBETR(I2) = BETAB(J284)
         ZOBETR(I1) = BETAB(J284)
         ZOJPR(I2)  = RJPAR(J284)
         ZOJPR(I1)  = RJPAR(J284)
         ZOJBSR(I2,1) = RJBSOS(J284,1)
         ZOJBSR(I1,1) = RJBSOS(J284,1)
         ZOJBSR(I2,2) = RJBSOS(J284,2)
         ZOJBSR(I1,2) = RJBSOS(J284,2)
         ZOTRR(I2)  = 1. - RFCIRC(J284)
         ZOTRR(I1)  = 1. - RFCIRC(J284)
C
  284    CONTINUE
C
         ZABR(NPSI1+1)  = 0.
         ZOQR(NPSI1+1)  = Q0
         ZOJBR(NPSI1+1) = RJDTB0
         ZOPPR(NPSI1+1) = DPDP0
         ZOPR(NPSI1+1)  = CP0
         ZOTTR(NPSI1+1) = DTTP0
         ZOTR(NPSI1+1)  = T0
         ZOIPR(NPSI1+1) = RIPR0
         ZOFR(NPSI1+1)   = 0.
         ZOJR(NPSI1+1)   = - (RMAG * DPDP0 + DTTP0 / RMAG)
         ZODQR(NPSI1+1)  = DQDP0
         ZOSHR(NPSI1+1)  = 0.
         ZOBETR(NPSI1+1) = BETAB(1)
         ZOTRR(NPSI1+1)  = 0.
         ZOJPR(NPSI1+1)  = RJPAR(1)
         ZOJBSR(NPSI1+1,1) = 0.
         ZOJBSR(NPSI1+1,2) = 0.
C
C  SET S-VALUES IN ZABS
C
         ZABSM(1) = 0.
         ZABIS(1) = 0.
         ZOQS(1)  = Q0
         ZODQS(1) = DQDP0
         ZOJBS(1) = RJDTB0
         ZOTRS(1) = 0.
         ZOHS(1)  = FCCCC0(HMERCR(1),HMERCR(2),HMERCR(3),HMERCR(4),
     ,                     CSM(1),CSM(2),CSM(3),CSM(4),0.)
         ZODIS(1) = FCCCC0(SMERCI(1),SMERCI(2),SMERCI(3),SMERCI(4),
     ,                     CSM(1),CSM(2),CSM(3),CSM(4),0.)
         ZODRS(1) = FCCCC0(SMERCR(1),SMERCR(2),SMERCR(3),SMERCR(4),
     ,                     CSM(1),CSM(2),CSM(3),CSM(4),0.)
         ZOPPS(1) = DPDP0
         ZOPS(1)  = CP0
         ZOTTS(1) = DTTP0
         ZOTS(1)  = T0
         ZOIPS(1) = RIPR0
         ZOIPS(1) = RIPR0
         ZOJPS(1) = RJPAR(1)
         ZOJBSS(1,1) = 0.
         ZOJBSS(1,2) = 0.
C
         DO 285 J285=1,NPSI1
C
         ZABSM(J285+1) = CSM(J285)
         ZABS(J285)    = CS(J285)
         ZABIS(J285+1) = (J285 - .5) / DFLOAT(NPSI)
         ZOHS(J285+1)  = HMERCR(J285)
         ZODIS(J285+1) = SMERCI(J285)
         ZODRS(J285+1) = SMERCR(J285)
         ZOBETS(J285)  = BETAB(J285)
         ZOQS(J285+1)  = QPSI(J285)
         ZODQS(J285+1) = CDQ(J285)
         ZOSHS(J285)   = CDRQ(J285)
         ZOPPS(J285+1) = CPPR(J285)
         ZOPS(J285+1)  = CPR(J285)
         ZOTTS(J285+1) = TTP(J285)
         ZOTS(J285+1)  = TMF(J285)
         ZOIPS(J285+1) = RIPR(J285)
         ZOJBS(J285+1) = RJDOTB(J285)
         ZOARS(J285)   = RSURF(J285)
         ZOJPS(J285+1) = RJPAR(J285)
         ZOJBSS(J285+1,1)= RJBSOS(J285,1)
         ZOJBSS(J285+1,2)= RJBSOS(J285,2)
         ZOTRS(J285+1) = 1. - RFCIRC(J285)
C
 285     CONTINUE
C
         ZABIS(NPSI1+1) = 1.
C
C  CONSTANT CHI LINES
C
         JSCHI = NCHI / 25 + 1
C
         DO 288 J288=1,NPSI1
C
         JNB = 0
C
         DO 287 J287=1,NCHI,JSCHI
C
         JNB             = JNB + 1
         ZRCHI(JNB,J288) = CR(J287,J288) - R0
         ZZCHI(JNB,J288) = CZ(J287,J288) - RZ0
C
 287     CONTINUE
 288     CONTINUE
C
         INBCHI = JNB
C
C  BALLOONING AND MERCIER
C
         WRITE(*,*) 'CS   DMERCI  DMERCR   IBALL'
         DO 289 J289=1,NPSI1
C
         IMERCI(J289) = 0
         IMERCR(J289) = 0
         IBALL(J289) = 0
C
         IF (SMERCI(J289) .LT. 0.) IMERCI(J289) = 1
         IF (SMERCR(J289) .LT. 0.) IMERCR(J289) = 1
         IF (NCBAL(J289) .NE. 0)  IBALL(J289) = 1
C
         WRITE(*,7979) CS(J289),SMERCI(J289),SMERCR(J289),IBALL(J289)
 289     CONTINUE
 7979    FORMAT(3E13.5,I2)
C
         DO 293 J293=1,NCURV
C
         ZRCURV(J293) = RRCURV(J293) - R0
         ZZCURV(J293) = RZCURV(J293) - RZ0
C
  293    CONTINUE
C
         ZCSIPR(1)      = 0.
         ZABIPR(NISO+1) = 1.
C
         DO 294 J294=1,NISO
C
         ZABIPR(J294)   = DFLOAT(J294-1) / DFLOAT(NISO)
         ZCSIPR(J294+1) = CSIPR(J294)
C
  294    CONTINUE
C
         DO 295 J295=1,NCHI1
C
         ZABIC(J295) = DFLOAT(J295-1) / DFLOAT(NCHI)
         ZCHI(J295)  = CHI(J295) - CHI(1)
C
  295    CONTINUE
C
         DO 296 J296=1,NS1
C
         ZABISG(J296) = DFLOAT(J296-1) / DFLOAT(NS)
C
  296    CONTINUE
C
         INSUR  = 6 * NT
         INS    = NPSI1 + 1
         INR    = 2 * NPSI1 + 1
         INTEXT = 120
C
C  WRITE DATA NECESSARY FOR PLOT ON THE FILE NUPLO
C
C
         WRITE(NUPLO,1003) INSUR,NCHI,NCHI1,NPSI,NPSI1,NS,NS1,NT,NT1,
     W                     INS,INR,INBCHI,INTEXT,NCURV,NMESHA,NMESHB,
     W                     NMESHC,NMESHD,NMESHE,NPOIDA,NPOIDB,
     W                     NPOIDC,NPOIDD,NPOIDE,NISO,NMGAUS,NPROFZ,
     W                     NRFP
         DO 297 J297=1,INTEXT
C
         WRITE(NUPLO,1005) TEXT(J297)
C
 297     CONTINUE
C
         WRITE(NUPLO,1004)(IBALL(J),J=1,NPSI1)
         WRITE(NUPLO,1004)(IMERCI(J),J=1,NPSI1)
         WRITE(NUPLO,1004)(IMERCR(J),J=1,NPSI1)
C
         WRITE(NUPLO,1006) SOLPDA,SOLPDB,SOLPDC,SOLPDD,SOLPDE,
     W                     ZRMAX,ZRMIN,ZZMAX,ZZMIN,PANGLE
         WRITE(NUPLO,1006)(APLACE(J),J=1,10)
         WRITE(NUPLO,1006)(AWIDTH(J),J=1,10)
         WRITE(NUPLO,1006)(BPLACE(J),J=1,10)
         WRITE(NUPLO,1006)(BWIDTH(J),J=1,10)
         WRITE(NUPLO,1006)(CPLACE(J),J=1,10)
         WRITE(NUPLO,1006)(CWIDTH(J),J=1,10)
         WRITE(NUPLO,1006)(DPLACE(J),J=1,10)
         WRITE(NUPLO,1006)(DWIDTH(J),J=1,10)
         WRITE(NUPLO,1006)(EPLACE(J),J=1,10)
         WRITE(NUPLO,1006)(EWIDTH(J),J=1,10)
         WRITE(NUPLO,1006) (ZTET(J),J=1,NT1)
         WRITE(NUPLO,1006) (CSIG(J),J=1,NS1)
         WRITE(NUPLO,1006) (CS(J),J=1,NPSI1)
         WRITE(NUPLO,1006) (ZCHI(J),J=1,NCHI1)
         WRITE(NUPLO,1006) (ZCSIPR(J),J=1,NISO+1)
         WRITE(NUPLO,1006)(ZRTET(J),J=1,NT)
         WRITE(NUPLO,1006)(ZZTET(J),J=1,NT)
         WRITE(NUPLO,1006)(ZRSUR(J),J=1,INSUR)
         WRITE(NUPLO,1006)(ZTSUR(J),J=1,INSUR)
         WRITE(NUPLO,1006)(ZZSUR(J),J=1,INSUR)
         WRITE(NUPLO,1006)(ZABIS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZABIT(J),J=1,NT1)
         WRITE(NUPLO,1006)(ZABIC(J),J=1,NCHI1)
         WRITE(NUPLO,1006)(ZOART(J),J=1,NT1)
         WRITE(NUPLO,1006)(ZABIPR(J),J=1,NISO+1)
         WRITE(NUPLO,1006)(ZABISG(J),J=1,NS1)
         WRITE(NUPLO,1006)(ZABSM(J),J=1,INS)
         WRITE(NUPLO,1006)(ZABR (J),J=1,INR)
         WRITE(NUPLO,1006)(ZOQS (J),J=1,INS)
         WRITE(NUPLO,1006)(ZOQR (J),J=1,INR)
         WRITE(NUPLO,1006)(ZODQS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZODQR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOSHS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOSHR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOJBS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOJBR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOJBSS(J,1),J=1,INS)
         WRITE(NUPLO,1006)(ZOJBSS(J,2),J=1,INS)
c%OS
         DO J=1,NPSI1
           ZOJBSS(J+1,2) = RJBSH(J)
         ENDDO
         WRITE(NUPLO,1006)(ZOJBSS(J,2),J=1,INS)
c%OS
         WRITE(NUPLO,1006)(ZOJBSR(J,1),J=1,INR)
         WRITE(NUPLO,1006)(ZOJBSR(J,2),J=1,INR)
         WRITE(NUPLO,1006)(ZOJPS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOJPR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOTRS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOTRR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOHS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZODIS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZODRS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOPPS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOPPR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOPS (J),J=1,INS)
         WRITE(NUPLO,1006)(ZOPR (J),J=1,INR)
         WRITE(NUPLO,1006)(ZOTTS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOTTR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOTS(J) ,J=1,INS)
         WRITE(NUPLO,1006)(ZOTR (J),J=1,INR)
         WRITE(NUPLO,1006)(ZOIPS(J),J=1,INS)
         WRITE(NUPLO,1006)(ZOIPR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOBETR(J),J=1,INR)
         WRITE(NUPLO,1006)(ZOBETS(J),J=1,NPSI1)
         WRITE(NUPLO,1006)(ZOFR (J),J=1,INR)
         WRITE(NUPLO,1006)(ZOARS(J),J=1,NPSI1)
         WRITE(NUPLO,1006)(ZOJR (J),J=1,INR)
         WRITE(NUPLO,1006)(ZABS(J),J=1,NPSI1)
C
         DO J=1,NPSI1
           DO I=1,NMGAUS*NT1
             RRISO(I,J) = RRISO(I,J) - R0
             RZISO(I,J) = RZISO(I,J) - RZ0
           ENDDO
           DO I=1,NCHI
             CR(I,J) = CR(I,J) - R0
             CZ(I,J) = CZ(I,J) - RZ0
           ENDDO
         ENDDO
C
         DO 301 J301=1,NPSI1
C
         WRITE(NUPLO,1006)(RRISO(I,J301),I=1,NMGAUS*NT1)
         WRITE(NUPLO,1006)(RZISO(I,J301),I=1,NMGAUS*NT1)
C
  301    CONTINUE
C
         WRITE(NUPLO,1006)(ZRCURV(J),J=1,NCURV)
         WRITE(NUPLO,1006)(ZZCURV(J),J=1,NCURV)
         WRITE(NUPLO,1006)((ZRCHI(J,I),J=1,INBCHI),I=1,NPSI1)
         WRITE(NUPLO,1006)((ZZCHI(J,I),J=1,INBCHI),I=1,NPSI1)
         WRITE(NUPLO,1006)((RSHEAR(J,I),J=1,NCHI),I=1,NPSI1)
         WRITE(NUPLO,1006)((CR(J,I),J=1,NCHI),I=1,NPSI1)
         WRITE(NUPLO,1006)((CZ(J,I),J=1,NCHI),I=1,NPSI1)
C
         DO J=1,NPSI1
           DO I=1,NMGAUS*NT1
             RRISO(I,J) = RRISO(I,J) + R0
             RZISO(I,J) = RZISO(I,J) + RZ0
           ENDDO
           DO I=1,NCHI
             CR(I,J) = CR(I,J) + R0
             CZ(I,J) = CZ(I,J) + RZ0
           ENDDO
         ENDDO
C
         RETURN
C
 1003    FORMAT(1X,13I8)
 1004    FORMAT((1X,40(I3)))
 1005    FORMAT(A132)
 1006    FORMAT((1X,8(1PE15.6)))
 1007    FORMAT(13I8)
 1111    FORMAT(A)
C        
         END
C*DECK C3SB03
C*CALL PROCESS
         SUBROUTINE SHAVE
C        ################
C
C                                        AUTHORS:
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SB03  SHAVE AWAY ALL POLOIDAL FLUX SURFACES WHICH DO NOT SATISFY  *
*         Q<QSHAVE OR PSI<PSISHAVE, WHICHEVER COMES FIRST.            *
*         THE SURFACE Q=QSHAVE OR PSI=PSISHAVE IS SAVED IN DISK FILE  *
*         NRZPEL, AND AN BE RESUSED AS INPUT WITH NSURF=6 FOR A       *
*         SUBSEQUENT RUN                                              *
*                                                                     *
***********************************************************************
C
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
      DO 10 J = NPSI1,1,-1
      ZQ = .5 * TMF(J) * CIDQ(J) / CPI
      WRITE(*,*) 'SHAVING Q: Q=',ZQ
      IF (ZQ.LE.QSHAVE.AND.CS(J)*CS(J).LE.PSISHAVE) GOTO 20
 10   CONTINUE
      WRITE(6,*) ' Q =',ZQ,' < QSHAVE =',QSHAVE
      WRITE(6,*) ' PSI =',CS(J)*CS(J),' < PSISHAVE =',PSISHAVE
      RETURN
C
 20   CONTINUE
      KSHAVE = J
      WRITE(*,*) 'NPSI1=',NPSI1,' KSHAVE=',KSHAVE
      WRITE(6,1000) 100.*CS(J)*CS(J),ZQ
      CALL SURFRZ(J,SIGPSI(1,J),TETPSI(1,J))
      RETURN
 1000 FORMAT(' R-Z COORDINATES OF ',F5.1,
     &       ' % FLUX SURFACE OUTPUT ON FILE RZPEEL',/,
     &       ' SURFACE Q =',F7.2)
      END
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C*DECK C3SB04
C*CALL PROCESS
         SUBROUTINE SURFRZ(K,PSIGMA,PTETA)
C        #################################
C                                        AUTHORS:
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SB04  COMPUTE AND SAVE R,Z COORDINATES OF FLUX SURFACE Q = QSHAVE *
*         OR PSI=PSISHAVE                                             *
*                                                                     *
***********************************************************************
C
C     WRITE R,Z COORDINATES OF K:TH FLUX SURFACE ON FILE RZPEEL
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMBND.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         PARAMETER (NPGMAP=NPMGS*NTP1)
C
         DIMENSION
     R   PSIGMA(*),         PTETA(*),
     R   ZBND(NPGMAP),      ZTETA(NPGMAP),
     R   D2RBP(NPGMAP),     D2ZBP(NPGMAP),
     R   ZA1(NPGMAP),       ZB1(NPGMAP),
     R   ZC1(NPGMAP),       CHINEW(NPBPS),
     I   IC(NPBPS),         IT(NPBPS)
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
         ZEPS = 1.E-3
C
         IGMAX = NMGAUS * NT1
C
         DO 1 J1=1,IGMAX
C
         ZTETA(J1) = PTETA(J1)
C
    1    CONTINUE
C
         CALL BOUND(IGMAX,ZTETA,ZBND)
C
         DO 4 J4=1,IGMAX
C 
         ZCOST  = COS(PTETA(J4))
         ZSINT  = SIN(PTETA(J4))
C
         ZRHO    = PSIGMA(J4) * ZBND(J4)
         ZR      = ZRHO * ZCOST + R0
         ZZ      = ZRHO * ZSINT + RZ0
C
         RRISO(J4,K)  = ZR
         RZISO(J4,K)  = ZZ
C
    4    CONTINUE
C
CYQL 2005-03-22
CSAVE SHAVED SURFACE TO EXPEQ.OUT VIA (RRBPSNEW,RZBPSNEW)
      DO J=1,IGMAX
         IF (ZTETA(J).LT.0.0) ZTETA(J)=ZTETA(J)+RC2PI
      ENDDO

C     CHIH = (ZTETA(IGMAX)-ZTETA(1))/(NBPS-1)
C     CHINEW(1) = ZTETA(1)
      CHIH = RC2PI/(NBPS-1)
      CHINEW(1) = 0.0
      DO J=1,NBPS-1
         CHINEW(J+1) = CHINEW(J) + CHIH
      ENDDO

      CALL SPLCY(ZTETA,RRISO(1,K),IGMAX-1,RC2PI,D2RBP,ZA1,ZB1,ZC1)
      CALL SPLCY(ZTETA,RZISO(1,K),IGMAX-1,RC2PI,D2ZBP,ZA1,ZB1,ZC1)
      D2RBP(IGMAX) = D2RBP(1)
      D2ZBP(IGMAX) = D2ZBP(1)
      
      CALL RESETI(IC,NBPS,1)
      DO JG=1,NBPS
         DO JT=1,IGMAX
            IF (IC(JG).EQ.0) GOTO 2
            IT(JG) = JT-1
            IF (ZTETA(JT).GE.CHINEW(JG)) IC(JG) = 0
 2          CONTINUE
         ENDDO
      ENDDO

      DO J=1,NBPS
         ICHI = IT(J)
         IF (ICHI.LT.1) ICHI = 1
         IF (ICHI.GT.IGMAX-1) ICHI = IGMAX-1
         ZH = ZTETA(ICHI+1) - ZTETA(ICHI)
         ZA = (ZTETA(ICHI+1) - CHINEW(J)) / ZH
         ZB = (CHINEW(J) - ZTETA(ICHI)) / ZH
         ZC = (ZA + 1)*(ZA - 1)*ZH*(ZTETA(ICHI+1)-CHINEW(J))/6.
         ZD = (ZB + 1)*(ZB - 1)*ZH*(CHINEW(J)-ZTETA(ICHI))/6.
         
         RRBPSNEW(J) = ZA*RRISO(ICHI,K) + ZB*RRISO(ICHI+1,K) + 
     &                 ZC*D2RBP(ICHI) + ZD*D2RBP(ICHI+1)
         RZBPSNEW(J) = ZA*RZISO(ICHI,K) + ZB*RZISO(ICHI+1,K) + 
     &                 ZC*D2ZBP(ICHI) + ZD*D2ZBP(ICHI+1)
      ENDDO

      OPEN(NRZPEL,FILE='RZPEEL',FORM='FORMATTED')
      WRITE(NRZPEL,1000) IGMAX
      WRITE(NRZPEL,1010) (RRISO(J,K),RZISO(J,K),J=1,IGMAX)
 1000 FORMAT(I5)
 1010 FORMAT(2E18.8)
      CLOSE(NRZPEL)
      RETURN
      END
C*DECK C3SB05
C*CALL PROCESS
         SUBROUTINE NERAT
C        ################
C                                        AUTHORS:
C                                        A. ROY,  CRPP-EPFL
***********************************************************************
*                                                                     *
* C3SB05  PRODUCE ERATO NAMELIST                                      *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMBLA.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMESH.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
C
         INCLUDE 'NEWRUN.inc'
         INCLUDE 'CUCCCC.inc'
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
***********************************************************************
*                                                                     *
*  . COMPUTE RINOR AND QCYL FOR ERATO DIAGNOSTICS                     *
*                                                                     *
***********************************************************************
C
         RINOR = RITOT / (ASPCT * TMF(NPSI1))
         QCYL  = 2 * AREA / (RITOT * RMAG**2)
C
         IF (NIDEAL.NE.1 .AND. NIDEAL.NE.2) RETURN
C
         CALL VZERO(ANGLE,16)
         CALL VZERO(QTILDA,NPPSI1)
         CALL VZERO(QS,NPPSI1)
C
         ALARG     = 24.
         AL0       = -1.5
         ANGLE(1)  = 0.
         ANGLE(2)  = 90.
C
         CALL VZERO(ANGLE(3),13)
C
         ARROW     = 0.05
         CST       = 1. / (Q0 * CPSRF)
         B2R2      = 0.
         EPSCON    = 1.E-04
         EPSMAC    = RC1M12
         REXT      = 1.E+05
         P0        = 0.
         QIAXE     = 1. / Q0
         WNL       = 0.
         WNTORE    = 2.
         WASPCT    = 1.1
         WCURV     = 0.
         WVERT     = 11.
         ITEST     = 1
         NAL0AUTO  = 0
         NPLEQ     = 0
         SMAX      = 2.
         NFIG      = 2
         NITMAX    = 10
         NV        = 7
         NVIT      = 1
C
C  DISK UNITS
C
         NDES      = 16
         MEQ       = 4
         NSAVE     = 8
         NVAC      = 17
         NUA1      = 31
         NUA2      = 32
         NUB1      = 33
         NUB2      = 34
         NUX       = 35
C
         CALL VZERO(WALL,10)
C
         NLDIAG(1) = .TRUE.
         NLDIAG(2) = .TRUE.
C
         NUPL      = 21
         NUSG      = 22
         NUWA      = 66
         NWALL     = 0
C
         NLEINQ    = .FALSE.
         NLGREN    = .TRUE.
C
         IF (NSYM .EQ. 0) THEN
C
            NLSYM = .FALSE.
C
         ELSE
C
            NLSYM = .TRUE.
C
         ENDIF
C
         QS(1)     = Q0
         QS(2)     = FCCCC0(QPSI(1),QPSI(2),QPSI(3),QPSI(4),
     ,                      CSM(1),CSM(2),CSM(3),CSM(4),CS(2))
         QTILDA(1) = QPSI(1)
         QTILDA(2) = QPSI(2)
C
         DO 3 J3=3,NPSI-1
C
         Z1        = QPSI(J3-2)
         Z2        = QPSI(J3-1)
         Z3        = QPSI(J3  )
         Z4        = QPSI(J3+1)
C
         QS(J3)     = FCCCC0(Z1,Z2,Z3,Z4,CSM(J3-2),
     ,                       CSM(J3-1),CSM(J3),CSM(J3+1),CS(J3))
         QTILDA(J3) = QPSI(J3)
C
    3    CONTINUE
C
         Z1        = QPSI(NPSI-2)
         Z2        = QPSI(NPSI-1)
         Z3        = QPSI(NPSI  )
         Z4        = QPSI(NPSI+1)
C
         QS(NPSI)  = FCCCC0(Z1,Z2,Z3,Z4,CSM(NPSI-2),
     ,                      CSM(NPSI-1),CSM(NPSI),CS(NPSI+1),CS(NPSI))
         QTILDA(NPSI)  = QPSI(NPSI)
         QS(NPSI1)     = QPSI(NPSI1)
         QTILDA(NPSI1) = QPSI(NPSI1)
         QSURF         = QPSI(NPSI1)
         TSURF         = TMF(NPSI1)
C
         DO 4 J4=3,20
C
         NLDIAG(J4) = .FALSE.
C
    4    CONTINUE
C
         READ(5,NEWRUN)
C
         IF (NLGREN .OR. (REXT .LE. 1.)) NV = 0
         IF ((.NOT. NLGREN) .AND. (NV.EQ.0) .AND. (REXT.GT.1.)) THEN
            STOP 'NLGREN?'
         ENDIF
C
         WRITE(NSAVE,NEWRUN)
C
         WRITE(NSAVE,98) NPSI1
         WRITE(NSAVE,99) (QS(J),J=1,NPSI1)
         WRITE(NSAVE,98) NPSI1
         WRITE(NSAVE,99) (QTILDA(J),J=1,NPSI1)
C
         WRITE(NPRNT,*)
         WRITE(NPRNT,*)
         WRITE(NPRNT,*)
         WRITE(NPRNT,*)
         WRITE(NPRNT,NEWRUN)
C
         CALL RARRAY('QS',QS,NPSI1)
         CALL RARRAY('QTILDA',QTILDA,NPSI1)
C
   98    FORMAT(I5)
   99    FORMAT(4E30.13)
C
         RETURN
         END
C*DECK C3SB06
C*CALL PROCESS
      SUBROUTINE GENOUT(VECTOR,ALPHA,NUNIT,NCRAY,NR,MSMAX,NSMAX,RM,RN)
C     ================================================================
C
C                                        AUTHORS:
C                                        G. VLAD,  ENEA FRASCATI
C
C.....GENERAL  INPUT ROUTINE:...........................................
C.....                       VECTOR        :  VECTOR....................
C.....                       ALPHA         :  ALPHANUMERIC STRING 6 CH..
C.....                       NUNIT         :  INPUT  UNIT NUMBER........
C.....                       NCRAY         :  RADIAL VECTOR DIMENSION...
C.....                       NR,MSMAX,NSMAX:  CURRENT VECTOR DIM........
C.....                       RM,RN         :  POLOIDAL AND TOROIDAL # ..
C
C...!!IF NSMAX .GT. 1, MSMAX MUST BE THE TRUE DIMENSION OF THE VECTORS..
C.....TO PRESERVE THE CORRECT ORDERING OF THE VECTOR ELEMENTS...........
C
C
         INCLUDE 'DECLAR.inc'
         COMPLEX    VECTOR
         DIMENSION
     C   VECTOR(*)
C
         DIMENSION             
     R   RM(*),   RN(*)
C
         CHARACTER*6      ALPHA
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
      DO 100 NS=1,NSMAX
        DO 100 MS=1,MSMAX
          WRITE(NUNIT,1010) ALPHA,RM(MS+(NS-1)*MSMAX),RN(NS)
          WRITE(NUNIT,1020)
     &         (VECTOR(I+(MS-1)*NCRAY+(NS-1)*MSMAX*NCRAY),I=1,NR)
  100 CONTINUE
C
 1010 FORMAT(//,1X,A,2F20.0)
 1020 FORMAT(2D30.15)
C
      RETURN
      END
C*DECK MR01
C*CALL PROCESS
         SUBROUTINE ALDLT(A,EPS,N,M,MP,NSING)
C        ------------------------------------
C
C     DECOMPOSE A=L*D*LT                                             
C                                                                    
C     VERSION 1C           13.9.74     RALF GRUBER    CRPP LAUSANNE  
C                                                                    
C     A IS A BAND MATRIX WITH HALF WIDTH M AND LENGTH N              
C     L CONTAINS 1 IN THE DIAGONAL                                   
C     AS OUTPUT D REPLACES THE DIAGONAL OF A AND                     
C     LT WITHOUT ITS DIAGONAL THE REST OF A                          
C     ALL CALCULATIONS ARE PERFORMED IN A                            
C     NSING = -1 WHEN A IS SINGULAR                                  
C
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   A(N*MP)
C
C     INITIALIZE
C
         M1  = MP - 1
         IKD = 0
         AD  = ABS(A(1)) * EPS
C
C     SCAN OVER THE WHOLE LENGTH OF A
C
         DO  4  JIB=2,N
            DIAG=A(IKD+1)
C
C     TEST FOR ZERO PIVOT
C
            IF (ABS(DIAG) .LT. AD) THEN
               NSING = -1
               RETURN
            ENDIF
C
C     RESTRICTION OF LOOP FOR NOT EXCEEDING BAND MATRIX
C
            LOPBND = M
            I1     = N - JIB + 2
C
            IF (I1 .LT. M) LOPBND = I1
C
C     DIAGONAL ELEMENT BEFORE GAUSS ELIMINATION
C
            IJ = IKD + MP + 1
            AD = ABS(A(IJ)) * EPS
C
C     SETS THE ROW OF THE TRANSPOSED LEFT HAND SIDE MATRIX LT
C
            DO 3 JJB=2,LOPBND
               ITOP    = IKD + JJB
               TOP     = A(ITOP)
               A(ITOP) = A(ITOP) / DIAG
C
C     GAUSS RECTANGULAR RULE GOING DOWNWARDS
C
               CALL SAXPY(JJB-1,-TOP,A(IKD+2),1,A(ITOP+M1),M1)
   3        CONTINUE
         IKD = IKD + MP
   4     CONTINUE
C
C     LAST DIAGONAL ELEMENT
C
         IKD = (N - 1) * MP + 1
         IJ  = IKD - MP
C
         IF (ABS(A(IKD)) .LT. ABS(A(IJ))*EPS) THEN
            NSING = -1
         ENDIF
C
         RETURN
         END
C*DECK MR02
C*CALL PROCESS
         SUBROUTINE LYV(A,U,N,NP,M,MP)
C        -----------------------------
C
C     SOLVE L*Y=V                                                    
C                                                                     
C     VERSION 1C           13.9.74     RALF GRUBER    CRPP LAUSANNE   
C                                                                     
C     L IS STORED AS LT IN THE OFF DIAGONAL PART OF THE BAND MATRIX A 
C     WITH HALF BAND WIDTH M AND LENGTH N .                           
C     U ARE THE KV VECTORS OF N COMPONENTS , REPRESENTING AT INPUT    
C     THE RIGHT SIDE PART OF THE SYSTEM OF LINEAR EQUATIONS .         
C     ALL CALCULATIONS ARE PERFORMED IN U .                           
C
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   A(N*MP),   U(NP)
C
         M1 = MP - 1
C
C     FIRST COMPONENT UNCHANGED (L CONTAINS 1 ON DIAGONAL)
C
         IKD = 1
C
C     SCAN OVER ALL COMPONENTS
C
         DO 2 J2=2,N
            LOPBND = M
C
            IF (J2 .LT. M) LOPBND = J2
C
            IKD = IKD + MP
            IJV = J2 + 1
C
c%OS  needs the lowest index as start-up, to keep positive dimension in SDOT
            UJ = U(J2)-SDOT(LOPBND-1,U(IJV-LOPBND),-1,
     +                               A(IKD-(LOPBND-1)*M1),-M1)
C
            U(IJV-1) = UJ
   2     CONTINUE
C
         RETURN
         END
C*DECK MR03
C*CALL PROCESS
         SUBROUTINE DWY(A,U,N,NP,M,MP)
C        ---------------------------------
C
C     SOLVE D*W=Y                                                    
C                                                                    
C     VERSION 1C           13.9.74     RALF GRUBER    CRPP LAUSANNE  
C                                                                    
C  U               KV SETS OF N-VECTORS, (BOTH X AND Y)              
C  A               DIAGONAL MATRIX D                                 
C                    (DIAGONAL IS FIRST COMPONENT OF BAND OF LENGTH M)
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   A(N*MP),   U(NP)
C
C---------------------------------------------------------------------
C
         IKD = 1
C
C      STORE VALUE
C
         DO 1 J1=1,N
            U(J1) = U(J1) / A(IKD)
            IKD = IKD + MP
   1     CONTINUE
C
         RETURN
         END
C*DECK MR04
C*CALL PROCESS
         SUBROUTINE LTXW(A,U,N,NP,M,MP)
C        ------------------------------
C
C     SOLVE LT*X=W                                                   
C                                                                    
C     VERSION 1C           13.9.74     RALF GRUBER    CRPP LAUSANNE  
C                                                                    
C     LT IS STORED IN THE OFF DIAGONAL PART OF THE BAND MATRIX A     
C     WITH HALF BAND WIDTH M AND LENGTH N .                          
C     U ARE THE KV VECTORS OF N COMPONENTS , REPRESENTING AT INPUT   
C     THE RIGHT SIDE PART OF THE SYSTEM OF LINEAR EQUATIONS .        
C     ALL CALCULATIONS ARE PERFORMED IN U .                          
C                                                                    
C     FIRST COMPONENT UNCHANGED ( LT CONTAINS 1 ON DIAGONAL)         
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   A(N*MP),   U(NP)
C
C     SCAN OVER ALL COMPONENTS
C
         IKD = (N - 1) * MP
C
         DO 2 J2=2,N
            LOPBND = M
C
            IF (J2 .LT. M) LOPBND = J2
C
            IJV = N - J2
            IKD = IKD - MP
C
            UJ = U(IJV+1)-SDOT(LOPBND-1,U(IJV+2),1,A(IKD+2),1)
C
            U(IJV+1) = UJ
  2      CONTINUE
C
         RETURN
         END
C*DECK MSP01
C*CALL PROCESS
         SUBROUTINE SPLINE(X,Y,N,YP2,A,B)
C        ################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MSP03  CUBIC SPLINE INTERPOLATION (CUBIC INTERP. FOR YP1, YPN)      *
*        X = DISCRETIZATION GRID                                      *
*        Y = FUNCTION TO BE INTERPOLATED                              *
*        N = NUMBER OF X GRID POINTS                                  *
*        YP2 = SECOND DERIVATIVE OF Y WITH RESPECT TO X               *
*        A,B = WORK ARRAYS                                            *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION X(N),Y(N),YP2(N),A(N),B(N)
C
         INCLUDE 'CUCCCC.inc'
C
         YP1 = FCCCC1(Y(1),Y(2),Y(3),Y(4),X(1),X(2),X(3),X(4),X(1))
         YPN = FCCCC1(Y(N-3),Y(N-2),Y(N-1),Y(N),
     ,                X(N-3),X(N-2),X(N-1),X(N),X(N))
C
         DO 2 J2=2,N-1
C
            A(J2)   = (X(J2+1) - X(J2-1)) / 3.
            B(J2)   = (X(J2+1) - X(J2  )) / 6.
            YP2(J2) = (Y(J2+1)-Y(J2)) / (X(J2+1)-X(J2))-
     &                (Y(J2)-Y(J2-1)) / (X(J2)-X(J2-1))
C
 2       CONTINUE
C
         A(1)   = (X(2) - X(1)) / 3.
         B(1)   = (X(2) - X(1)) / 6.
         YP2(1) = (Y(2)-Y(1)) / (X(2)-X(1)) - YP1
         A(N)   = (X(N) - X(N-1)) / 3.
         YP2(N) = YPN - (Y(N)-Y(N-1)) / (X(N)-X(N-1))
C
         CALL TRIDAG(A,B,YP2,N,RC1M14)
C
         RETURN
         END
C*DECK MSP02
C*CALL PROCESS
         SUBROUTINE MSPLINE(X,Y,N,MD,M,YP1,YPN,YP2,A,B,WORK)
C        ###################################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MSP02  SOLVE M CUBIC SPLINE INTERPOLATIONS IN PARALLEL.             *
*        X = DISCRETIZATION GRID                                      *
*        Y = FUNCTION TO BE INTERPOLATED                              *
*        N = NUMBER OF X GRID POINTS                                  *
*        MD = LEADING DIMENSION OF ARRAYS                             *
*        M = NUMBER OF SPLINES TO BE SOLVED                           *
*        YP1 = ARRAY OF FIRST DERIVATIVES AT LEFT OF X GRID           *
*        YPN = ARRAY OF FIRST DERIVATIVES AT RIGHT OF X GRID          *
*        YP2 = SECOND DERIVATIVE OF Y WITH RESPECT TO X               *
*        A,B,WORK = WORK ARRAYS                                       *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(MD,N), B(MD,N), WORK(MD),
     &             X(MD,N), Y(MD,N), YP2(MD,N),
     &             YP1(MD), YPN(MD)

         DO 2 J2=2,N-1
            DO 1 J1=1,M
C
               A(J1,J2)   = (X(J1,J2+1) - X(J1,J2-1)) / 3.
               B(J1,J2)   = (X(J1,J2+1) - X(J1,J2  )) / 6.
               YP2(J1,J2) = (Y(J1,J2+1)-Y(J1,J2)) / 
     &                      (X(J1,J2+1)-X(J1,J2)) -
     &                      (Y(J1,J2)-Y(J1,J2-1)) / 
     &                      (X(J1,J2)-X(J1,J2-1))
C
 1          CONTINUE
 2       CONTINUE
C
         DO 3 J3=1,M
C
            A(J3,1)   = (X(J3,2) - X(J3,1)) / 3.
            B(J3,1)   = (X(J3,2) - X(J3,1)) / 6.
            YP2(J3,1) = (Y(J3,2)-Y(J3,1)) / 
     &                  (X(J3,2)-X(J3,1)) - YP1(J3)
            A(J3,N)   = (X(J3,N) - X(J3,N-1)) / 3.
            YP2(J3,N) = YPN(J3) - (Y(J3,N)-Y(J3,N-1)) / 
     &                  (X(J3,N)-X(J3,N-1))
C
 3       CONTINUE
C
         CALL TRIDAGM(A,B,YP2,WORK,N,MD,M,RC1M14)
C
         RETURN
         END
C*DECK MSP03
C*CALL PROCESS
         SUBROUTINE SPLCY(X,Y,N,PERIOD,YP2,A,B,C)
C        ########################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MSP03  CUBIC SPLINE INTERPOLATION, Y PERIODIC.                      *
*        X = DISCRETIZATION GRID                                      *
*        Y = FUNCTION TO BE INTERPOLATED                              *
*        N = NUMBER OF X GRID POINTS                                  *
*        PERIOD = Y(X+PERIOD) = Y(X)                                  *
*        YP2 = SECOND DERIVATIVE OF Y WITH RESPECT TO X               *
*        A,B,C = WORK ARRAYS                                          *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(N), B(N), C(N),
     &             X(N), Y(N), YP2(N)
C
         DO 2 J2=2,N-1
C
            A(J2)   = (X(J2+1) - X(J2-1)) / 3.
            B(J2)   = (X(J2+1) - X(J2  )) / 6.
            C(J2)   = 0.
            YP2(J2) = (Y(J2+1)-Y(J2)) / (X(J2+1)-X(J2))-
     &                (Y(J2)-Y(J2-1)) / (X(J2)-X(J2-1))
C
 2       CONTINUE
C
         A(1)   = (PERIOD + X(2) - X(N)) / 3.
         B(1)   = (X(2) - X(1)) / 6.
         C(1)   = (PERIOD + X(1) - X(N)) / 6.
         YP2(1) = (Y(2) - Y(1)) / (X(2) - X(1)) - 
     &            (Y(1) - Y(N)) / (PERIOD + X(1) - X(N))
         C(N-1) = B(N-1)
         C(N)   = (PERIOD + X(1) - X(N-1)) / 3.
         YP2(N) = (Y(1) - Y(N)) / (PERIOD + X(1) - X(N))- 
     &            (Y(N) - Y(N-1)) / (X(N) - X(N-1))
C
 3       CONTINUE
C
         CALL TRICYC(A,B,C,YP2,N,RC1M14)
C
         RETURN
         END
C*DECK MSP04
C*CALL PROCESS
         SUBROUTINE MSPLCY(X,Y,N,MD,M,PERIOD,YP2,A,B,C,WORK)
C        ###################################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MSP04  SOLVE M CUBIC SPLINE INTERPOLATIONS IN PARALLEL OF PERIODIC  *
*        FUNCTINS.                                                    *
*        X = DISCRETIZATION GRID                                      *
*        Y = FUNCTION TO BE INTERPOLATED                              *
*        N = NUMBER OF X GRID POINTS                                  *
*        MD = LEADING DIMENSION OF ARRAYS                             *
*        M = NUMBER OF SPLINES TO BE SOLVED                           *
*        PERIOD = Y(M,X+PERIOD) = Y(M,X)                              *
*        YP2 = SECOND DERIVATIVE OF Y WITH RESPECT TO X               *
*        A,B,C,WORK = WORK ARRAYS                                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(MD,N), B(MD,N), C(MD,N), WORK(MD),
     &             X(MD,N), Y(MD,N), YP2(MD,N)
C
         DO 2 J2=2,N-1
            DO 1 J1=1,M
C
               A(J1,J2)   = (X(J1,J2+1) - X(J1,J2-1)) / 3.
               B(J1,J2)   = (X(J1,J2+1) - X(J1,J2  )) / 6.
               C(J1,J2)   = 0.
               YP2(J1,J2) = (Y(J1,J2+1)-Y(J1,J2)) / 
     &                      (X(J1,J2+1)-X(J1,J2)) -
     &                      (Y(J1,J2)-Y(J1,J2-1)) / 
     &                      (X(J1,J2)-X(J1,J2-1))
C
 1          CONTINUE
 2       CONTINUE
C
         DO 3 J3=1,M
C
            A(J3,1)   = (PERIOD + X(J3,2) - X(J3,N)) / 3.
            B(J3,1)   = (X(J3,2) - X(J3,1)) / 6.
            C(J3,1)   = (PERIOD + X(J3,1) - X(J3,N)) / 6.
            YP2(J3,1) = (Y(J3,2)-Y(J3,1)) / (X(J3,2)-X(J3,1)) - 
     &                  (Y(J3,1)-Y(J3,N)) / (PERIOD+X(J3,1)-X(J3,N))
            C(J3,N-1) = B(J3,N-1)
            C(J3,N)   = (PERIOD + X(J3,1) - X(J3,N-1)) / 3.
            YP2(J3,N) = (Y(J3,1)-Y(J3,N)) / (PERIOD+X(J3,1)-X(J3,N))- 
     &                  (Y(J3,N)-Y(J3,N-1)) / (X(J3,N)-X(J3,N-1))
C
 3       CONTINUE
C
         CALL TRICYCM(A,B,C,YP2,WORK,N,MD,M,RC1M14)
C
         RETURN
         END
C*DECK MSP05
C*CALL PROCESS
         SUBROUTINE SPLCYP(X,Y,N,XPERIOD,YPERIOD,YP2,A,B,C)
C        ##################################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MSP05  CUBIC SPLINE INTERPOLATION, X AND Y PERIODIC.                *
*        X = DISCRETIZATION GRID                                      *
*        Y = FUNCTION TO BE INTERPOLATED                              *
*        N = NUMBER OF X GRID POINTS                                  *
*        XPERIOD = X+XPERIOD = X                                      *
*        YPERIOD = Y(X+YPERIOD) = Y(X)                                *
*        YP2 = SECOND DERIVATIVE OF Y WITH RESPECT TO X               *
*        A,B,C = WORK ARRAYS                                          *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(N), B(N), C(N),
     &             X(N), Y(N), YP2(N)
C
         DO 2 J2=2,N-1
C
            A(J2)   = (X(J2+1) - X(J2-1)) / 3.
            B(J2)   = (X(J2+1) - X(J2  )) / 6.
            C(J2)   = 0.
            YP2(J2) = (Y(J2+1)-Y(J2)) / (X(J2+1)-X(J2))-
     &                (Y(J2)-Y(J2-1)) / (X(J2)-X(J2-1))
C
 2       CONTINUE
C
         A(1)   = (XPERIOD + X(2) - X(N)) / 3.
         B(1)   = (X(2) - X(1)) / 6.
         C(1)   = (XPERIOD + X(1) - X(N)) / 6.
         YP2(1) = (Y(2) - Y(1)) / (X(2) - X(1)) - 
     &            (YPERIOD + Y(1) - Y(N)) / (XPERIOD + X(1) - X(N))
         C(N-1) = B(N-1)
         C(N)   = (XPERIOD + X(1) - X(N-1)) / 3.
         YP2(N) = (YPERIOD + Y(1) - Y(N)) / (XPERIOD + X(1) - X(N)) - 
     &            (Y(N) - Y(N-1)) / (X(N) - X(N-1))
C
 3       CONTINUE
C
         CALL TRICYC(A,B,C,YP2,N,RC1M14)
C
         RETURN
         END
C*DECK MRD01
C*CALL PROCESS
         SUBROUTINE NTRIDG(A,KDIMA,ND1,ND2,N)
C        ####################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MRD01 PERFORM L*U DECOMPOSITION FOR ND2 - ND1 TRIDIAGONAL SYMETRIC  *
*       REAL MATRIXES IN PARALLEL.                                    *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
C
         DIMENSION
     D      DIAG(NPISO),    A(KDIMA,2,N)
C
         DO 4 J4=2,N
C
         ITOP = J4 - 1
C
         DO 2 J2=ND1,ND2
C
         DIAG(J2) = A(J2,1,ITOP)
C
    2    CONTINUE
C
         DO 3 J3=ND1,ND2
C
         TOP    = A(J3,2,ITOP)
C
         A(J3,2,ITOP) = A(J3,2,ITOP) / DIAG(J3)
         A(J3,1,J4)   = A(J3,1,J4) - TOP * A(J3,2,ITOP)
C
    3    CONTINUE
    4    CONTINUE
C
         RETURN
         END
C*DECK MRD02
C*CALL PROCESS
         SUBROUTINE TRIDAGM(A,B,R,DIAG,N,MD,M,EPS)
C        #########################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MRD02 PERFORM L*U DECOMPOSITION AND BACK-SUBSTITUTION OF M          *
*       TRIDIAGONAL SYSTEMS                                           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(MD,N),B(MD,N),R(MD,N),DIAG(MD)
C
C DECOMPOSE AND FORWARD SUBSTITUTION
C
         DO 4 J4=2,N
C
            ITOP = J4 - 1
C
            DO 2 J2=1,M
C
               DIAG(J2) = A(J2,ITOP)
C
 2          CONTINUE
C
            IMN = ISAMIN(M,DIAG,1)
C
            IF (DIAG(IMN) .LT. EPS) THEN
C
               WRITE(*,*) ' ZERO PIVOT I = ',J4,', M = ',IMN
               STOP
C
            ENDIF
C
            DO 3 J3=1,M
C
               TOP  = B(J3,ITOP)
C
               B(J3,ITOP) = B(J3,ITOP) / DIAG(J3)
               R(J3,ITOP) = R(J3,ITOP) / DIAG(J3)
C
               A(J3,J4) = A(J3,J4) - TOP * B(J3,ITOP)
               R(J3,J4) = R(J3,J4) - TOP * R(J3,ITOP)
C 
 3          CONTINUE
 4       CONTINUE
C
C CHECK LAST PIVOT
C
         IMN = ISAMIN(M,A(1,N),1)
C
         IF (A(IMN,N) .LT. EPS) THEN
C
            WRITE(*,*) ' ZERO PIVOT I = ',N,', M = ',IMN
            STOP
C
         ENDIF
C
C BACKSUBSTITUTION
C
         DO 5 J5=1,M
C
            R(J5,N) = R(J5,N) / A(J5,N)
C
 5       CONTINUE
c
         DO 10 J10=N-1,1,-1
            DO 6 J6=1,M
C
               R(J6,J10) = R(J6,J10) - B(J6,J10) * R(J6,J10+1)
C
  6         CONTINUE
 10      CONTINUE
C
         RETURN
         END
C*DECK MRD03
C*CALL PROCESS
         SUBROUTINE TRICYC(A,B,C,R,N,EPS)
C        ################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MRD03 PERFORM L*U DECOMPOSITION AND BACK-SUBSTITUTION OF 1          *
*       TRIDIAGONAL SYSTEMS WITH PERIODIC BOUNDARY CONDITIONS         *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(N),B(N),C(N),R(N)
C
C DECOMPOSE AND FORWARD SUBSTITUTION
C
         DO 4 J4=2,N-1
C
            ITOP = J4 - 1
            DIAG = A(ITOP)
C
            IF (DIAG .LT. EPS) THEN
C
               WRITE(*,*) ' ZERO PIVOT I = ',J4,' A(I-1)=',DIAG
               STOP
C
            ENDIF
C
            BTOP = B(ITOP)
            CTOP = C(ITOP)
C
            B(ITOP) = B(ITOP) / DIAG
            C(ITOP) = C(ITOP) / DIAG
            R(ITOP) = R(ITOP) / DIAG
C
            A(J4) = A(J4) - BTOP * B(ITOP)
            C(J4) = C(J4) - BTOP * C(ITOP)
            R(J4) = R(J4) - BTOP * R(ITOP)
C
            C(N) = C(N) - CTOP * C(ITOP)
            R(N) = R(N) - CTOP * R(ITOP)
C
 4       CONTINUE
C
         ITOP = N - 1
         DIAG = A(ITOP)
C
         IF (DIAG .LT. EPS) THEN
C
            WRITE(*,*) ' ZERO PIVOT I = ',N-1
            STOP
C
         ENDIF
C
         TOP  = C(ITOP)
C
         C(ITOP) = C(ITOP) / DIAG
         R(ITOP) = R(ITOP) / DIAG
         C(N)    = C(N) - TOP * C(ITOP)
         R(N)    = R(N) - TOP * R(ITOP)
C
C CHECK LAST PIVOT
C
         IF (C(N) .LT. EPS) THEN
C
            WRITE(*,*) ' ZERO PIVOT I = ',N
            STOP
C
         ENDIF
C
C BACKSUBSTITUTION
C
         R(N) = R(N) / C(N)
C
         DO 9 J9=1,N-1
C
            R(J9) = R(J9) - C(J9) * R(N)
C
 9       CONTINUE
C
         DO 11 J11=N-2,1,-1
C
            R(J11) = R(J11) - B(J11) * R(J11+1)
C
 11      CONTINUE
C
         RETURN
         END
C*DECK MRD04
C*CALL PROCESS
         SUBROUTINE TRICYCM(A,B,C,R,DIAG,N,MD,M,EPS)
C        ###########################################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MRD02 PERFORM L*U DECOMPOSITION AND BACK-SUBSTITUTION OF M          *
*       TRIDIAGONAL SYSTEMS                                           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(MD,N),B(MD,N),C(MD,N),R(MD,N),DIAG(MD)
C
C DECOMPOSE AND FORWARD SUBSTITUTION
C
         DO 4 J4=2,N-1
C
            ITOP = J4 - 1
C
            DO 2 J2=1,M
C
               DIAG(J2) = A(J2,ITOP)
C
 2          CONTINUE
C
            IMN = ISAMIN(M,DIAG,1)
C
            IF (DIAG(IMN) .LT. EPS) THEN

C
               WRITE(*,*) ' ZERO PIVOT I = ',J4,', M = ',IMN
               STOP
C
            ENDIF

C
            DO 3 J3=1,M
C
               BTOP = B(J3,ITOP)
               CTOP = C(J3,ITOP)
C
               B(J3,ITOP) = B(J3,ITOP) / DIAG(J3)
               C(J3,ITOP) = C(J3,ITOP) / DIAG(J3)
               R(J3,ITOP) = R(J3,ITOP) / DIAG(J3)
C
               A(J3,J4) = A(J3,J4) - BTOP * B(J3,ITOP)
               C(J3,J4) = C(J3,J4) - BTOP * C(J3,ITOP)
               R(J3,J4) = R(J3,J4) - BTOP * R(J3,ITOP)
C
               C(J3,N) = C(J3,N) - CTOP * C(J3,ITOP)
               R(J3,N) = R(J3,N) - CTOP * R(J3,ITOP)

C
 3          CONTINUE
 4       CONTINUE
C
         ITOP = N - 1
C
         DO 5 J5=1,M
C
            DIAG(J5) = A(J5,ITOP)
C
 5       CONTINUE
C
         IMN = ISAMIN(M,DIAG,1)
C
         IF (DIAG(IMN) .LT. EPS) THEN
C
            WRITE(*,*) ' ZERO PIVOT I = ',N-1,', M = ',IMN
            STOP
C
         ENDIF
C
         DO 6 J6=1,M
C
            TOP  = C(J6,ITOP)
C
            C(J6,ITOP) = C(J6,ITOP) / DIAG(J6)
            R(J6,ITOP) = R(J6,ITOP) / DIAG(J6)
            C(J6,N)    = C(J6,N) - TOP * C(J6,ITOP)
            R(J6,N)    = R(J6,N) - TOP * R(J6,ITOP)
C
 6       CONTINUE
C
C CHECK LAST PIVOT
C
         IMN = ISAMIN(M,C(1,N),1)
C
         IF (C(IMN,N) .LT. EPS) THEN
C
            WRITE(*,*) ' ZERO PIVOT I = ',N,', M = ',IMN
            STOP
C
         ENDIF
C
C BACKSUBSTITUTION
C
         DO 7 J7=1,M
C
            R(J7,N) = R(J7,N) / C(J7,N)

C
 7       CONTINUE
C
         DO 9 J9=1,N-1
            DO 8 J8=1,M
C
               R(J8,J9) = R(J8,J9) - C(J8,J9) * R(J8,N)
C
 8          CONTINUE
 9       CONTINUE
C
         DO 11 J11=N-2,1,-1
            DO 10 J10=1,M
C
               R(J10,J11) = R(J10,J11) - B(J10,J11) * R(J10,J11+1)
C
 10         CONTINUE
 11      CONTINUE
C
         RETURN
         END
C*DECK MRD05
C*CALL PROCESS
         SUBROUTINE TRIDAG(A,B,R,N,EPS)
C        ##############################
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* MRD02 PERFORM L*U DECOMPOSITION AND BACK-SUBSTITUTION OF A          *
*       TRIDIAGONAL SYSTEMS                                           *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         DIMENSION A(N),B(N),R(N)
C
C DECOMPOSE AND FORWARD SUBSTITUTION
C
         DO 4 J4=2,N
C
            ITOP = J4 - 1
C
            DIAG = A(ITOP)
C
            IF (DIAG .LT. EPS) THEN
C
               WRITE(*,*) ' ZERO PIVOT I = ',J4
               STOP
C
            ENDIF
C
            TOP  = B(ITOP)
C
            B(ITOP) = B(ITOP) / DIAG
            R(ITOP) = R(ITOP) / DIAG
C
            A(J4) = A(J4) - TOP * B(ITOP)
            R(J4) = R(J4) - TOP * R(ITOP)
C 
 4       CONTINUE
C
C CHECK LAST PIVOT
C
         IF (A(N) .LT. EPS) THEN
C
            WRITE(*,*) ' ZERO PIVOT I = ',N
            STOP
C
         ENDIF
C
C BACKSUBSTITUTION
C
         R(N) = R(N) / A(N)
c
         DO 10 J10=N-1,1,-1
C
            R(J10) = R(J10) - B(J10) * R(J10+1)
C
 10      CONTINUE
C
         RETURN
         END
C*DECK MRD06
C*CALL PROCESS
        SUBROUTINE SORT(N,A,B,C,D,WORK,IWORK)
C       -----------------------------------
C BRUTE FORCE SORTING OF A, B, C, D SO THAT A IS IN INCREASING
C ORDER. COULD DONE BE MORE ELEGANTLY WITH A HEAPSORT ALGORITHM...
C IS ONLY USED TO SORT EXPERIMENTAL BOUNDARY POINTS ONCE PER RUN.
C
        INCLUDE 'DECLAR.inc'
        DIMENSION A(N),B(N),C(N),D(N),WORK(N),IWORK(N)
C
        IMX = ISMAX(N,A,1)
        EPS = ABS(A(IMX))/10.
        CALL SCOPY(N,A,1,WORK,1)
C
        DO 1 J1=1,N
          IWORK(J1) = ISMIN(N,WORK,1)
 1        WORK(IWORK(J1)) = WORK(IMX)+EPS
C
        CALL SCOPY(N,A,1,WORK,1)
C
        DO 2 J2=1,N
 2        A(J2)=WORK(IWORK(J2))
C
        CALL SCOPY(N,B,1,WORK,1)
C
        DO 3 J3=1,N
 3        B(J3)=WORK(IWORK(J3))
C
        CALL SCOPY(N,C,1,WORK,1)
C
        DO 4 J4=1,N
 4        C(J4)=WORK(IWORK(J4))
C
        CALL SCOPY(N,D,1,WORK,1)
C
        DO 5 J5=1,N
 5        D(J5)=WORK(IWORK(J5))
C
        RETURN
        END
C*DECK U1
C*CALL PROCESS
         SUBROUTINE PAGE
C        ---------------
C
C  FETCH A NEW PAGE
C
         INCLUDE 'DECLAR.inc'
         WRITE (6,9900)
C
         RETURN
 9900    FORMAT('1')
         END
C*DECK U2
C*CALL PROCESS
         SUBROUTINE BLINES(K)
C        --------------------
C
C  INSERT K BLANK LINES ON OUTPUT CHANEL UNIT 6
C
         INCLUDE 'DECLAR.inc'
         DO 1 J=1,K
  1         WRITE (6,9900)
C
         RETURN
 9900    FORMAT(' ')
         END
C*DECK U10
C*CALL PROCESS
         SUBROUTINE MESAGE(KMESS)
C        ------------------------
C
C  PRINT A MESSAGE ON OUPUT CHANNEL UNIT 6 (JUMP A BLANK LINE FIRST)
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KMESS
C
         WRITE (6,9900) KMESS
C
         RETURN
 9900    FORMAT(/,1X,A)
         END
C*DECK U20
C*CALL PROCESS
         SUBROUTINE RVAR(KNAME,PVALUE)
C        -----------------------------
***********************************************************************
*                                                                     *
* U.20   PRINT NAME AND VALUE OF REAL VARIABLE                        *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME
C
         WRITE (6,9900) KNAME, PVALUE
C
         RETURN
 9900    FORMAT(/,1X,A,' =',1PE17.8)
         END
C*DECK U21
C*CALL PROCESS
         SUBROUTINE RVAR2(KN1,PV1,KN2,PV2)
C        #################################
***********************************************************************
*                                                                     *
* U.21   PRINT NAME AND VALUE OF TWO REAL VARIABLES                   *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KN1, KN2
C
         WRITE (6,9900) KN1,PV1,KN2,PV2
C
         RETURN
 9900    FORMAT(/,1X,A,' =',1PE17.8,T36,A,' =',1PE17.8)
         END
C*DECK U22
C*CALL PROCESS
         SUBROUTINE RVAR3(KN1,PV1,KN2,PV2,KN3,PV3)
C        #########################################
***********************************************************************
*                                                                     *
* U.22   PRINT NAME AND VALUE OF THREE REAL VARIABLES                 *
*                                                                     *
***********************************************************************
C
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KN1, KN2, KN3
C
         WRITE (6,9900) KN1,PV1,KN2,PV2,KN3,PV3
C
         RETURN
 9900    FORMAT(/,1X,A,' =',1PE17.8,T36,A,' =',
     +                                     1PE17.8,T71,A,' =',1PE17.8)
         END
C*DECK U23
C*CALL PROCESS
         SUBROUTINE IVAR(KNAME,KVALUE)
C        #############################
***********************************************************************
*                                                                     *
* U.23   PRINT NAME AND VALUE OF INTEGER VARIABLE                     *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME
C
         WRITE (6,9900) KNAME, KVALUE
C
         RETURN
 9900    FORMAT(/,1X,A,' =',I12)
         END
C*DECK U24
C*CALL PROCESS
         SUBROUTINE IVAR2(KN1,KV1,KN2,KV2)
C        ---------- -----
C
C U.24   PRINT NAME AND VALUE OF TWO INTEGER VARIABLES
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KN1, KN2
C
         WRITE (6,9900) KN1,KV1,KN2,KV2
C
         RETURN
 9900    FORMAT(/1X,A,' =',I12,T36,A,' =',I12)
         END
C*DECK U25
C*CALL PROCESS
         SUBROUTINE IVAR3(KN1,KV1,KN2,KV2,KN3,KV3)
C        ---------- -----
C
C U.25   PRINT NAME AND VALUE OF THREE INTEGER VARIABLES
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KN1, KN2, KN3
C
         WRITE (6,9900) KN1,KV1,KN2,KV2,KN3,KV3
C
         RETURN
 9900    FORMAT(/,1X,A,' =',I12,T36,A,' =',I12,T71,A,' =',I12)
         END
C*DECK U26
C*CALL PROCESS
         SUBROUTINE HVAR(KNAME,KVALUE)
C        ---------- ----
C
C U.26   PRINT NAME AND VALUE OF CHARACTER-TYPE VARIABLE
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME, KVALUE
C
C
         WRITE (6,9900) KNAME, KVALUE
C
         RETURN
 9900    FORMAT(/,1X,A,' =',A)
         END
C*DECK U27
C*CALL PROCESS
         SUBROUTINE LVAR(KNAME,KLVAL)
C        ---------- ----
C
C U.27   PRINT NAME AND VALUE OF LOGICAL VARIABLE
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME
         LOGICAL KLVAL
C
         IF (KLVAL)       WRITE (6,9901) KNAME
         IF (.NOT. KLVAL) WRITE (6,9902) KNAME
C
         RETURN
 9901    FORMAT(/,1X,A,' =      .TRUE.')
 9902    FORMAT(/,1X,A,' =      .FALSE.')
         END
C*DECK U30 
C*CALL PROCESS
         SUBROUTINE RARRAY(KNAME,PA,KDIM)
C        ---------- ------
C
C U.30   PRINT NAME AND VALUES OF REAL ARRAY
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME
         DIMENSION
     R   PA(KDIM)
C
         CALL BLINES(1)
         WRITE (6,9900) KNAME
         WRITE (6,9901) (PA(J),J=1,KDIM)
C
         RETURN
 9900    FORMAT(1X,A)
 9901    FORMAT((1X,8(1PE13.4)))
         END
C*DECK U31
C*CALL PROCESS
         SUBROUTINE IARRAY(KNAME,KA,KDIM)
C        ---------- ------
C
C U.31   PRINT NAME AND VALUES OF INTEGER ARRAY
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME
         DIMENSION KA(KDIM)
C
         CALL BLINES(1)
         WRITE (6,9900) KNAME
         WRITE (6,9901) (KA(J),J=1,KDIM)
C
         RETURN
 9900    FORMAT(1X,A)
 9901    FORMAT((1X,8(I13)))
         END
C*DECK U32
C*CALL PROCESS
         SUBROUTINE LARRAY(KNAME,KLA,KDIM)
C        ---------- ------
C
C U.32   PRINT NAME AND VALUES OF LOGICAL ARRAY
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         LOGICAL KLA
         CHARACTER*(*) KNAME
         DIMENSION KLA(KDIM)
C
         CALL BLINES(1)
         WRITE (6,9900) KNAME
         WRITE (6,9901) (KLA(J),J=1,KDIM)
C
         RETURN
 9900    FORMAT(1X,A)
 9901    FORMAT((1X,10(L12)))
         END
C*DECK U33
C*CALL PROCESS
         SUBROUTINE HARRAY(KNAME,KA,KDIM)
C        ---------- ------
C
C U.33   PRINT NAME AND VALUES OF CHARACTER ARRAY
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME, KA
         DIMENSION KA(KDIM)
C
         CALL BLINES(1)
         WRITE (6,9900) KNAME
         WRITE (6,9901) (KA(J),J=1,KDIM)
C
         RETURN
 9900    FORMAT(1X,A)
 9901    FORMAT((1X,10(2X,A10)))
         END
C*DECK U34
C*CALL PROCESS
         SUBROUTINE SARRAY(KNAME,PSCALE,PA,KDIM)
C        ---------- ------
C
C U.34   PRINT NAME AND SCALED VALUES OF REAL ARRAY
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME
         DIMENSION
     R   PA(KDIM)
C
         CALL BLINES(1)
         WRITE (6,9900) KNAME
         WRITE (6,9901) (PSCALE*PA(J),J=1,KDIM)
C
         RETURN
 9900    FORMAT(1X,A)
 9901    FORMAT((1X,10(1PE13.4)))
         END
C*DECK U35
C*CALL PROCESS
         SUBROUTINE OARRAY(NCHAN,KNAME,PA,KDIM)
C        ---------- ------
C
C U.35   SAME AS RARRAY BUT TO FILE NUMBER NCHAN
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) KNAME
         DIMENSION
     R   PA(KDIM)
C
         WRITE (NCHAN,9900) KNAME
         WRITE (NCHAN,9901) (PA(J),J=1,KDIM)
C
         RETURN
 9900    FORMAT(1X,A)
 9901    FORMAT((1X,10(1PE13.4)))
         END
C*DECK U36
C*CALL PROCESS
         SUBROUTINE WRTEXT(VAR,CTEXT,TEXT,K,KN)
C        --------------------------------------
C
         INCLUDE 'DECLAR.inc'
         DIMENSION VAR(*)
         CHARACTER*(*) TEXT,CTEXT
C
         IF (K .EQ. 1) THEN
            WRITE(TEXT,9900) CTEXT,VAR(1)
         ELSE IF (K .EQ. 2) THEN
            WRITE(TEXT,9901) CTEXT,VAR(1)
         ELSE IF (K .EQ. 3) THEN
            WRITE(TEXT,9902) CTEXT,(VAR(L),L=1,KN)
         ENDIF
C
         RETURN
 9900    FORMAT(A,' = ',1PG13.5)
 9901    FORMAT(A,' = ',1PG13.5,'%')
 9902    FORMAT(A,' = ',13(1PG12.4))
         END
C*DECK U37
C*CALL PROCESS
         SUBROUTINE WITEXT(IVAR,CTEXT,TEXT,K)
C        ------------------------------------
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) TEXT,CTEXT
C
         IF (K .EQ. 1) THEN
            WRITE(TEXT,9900) CTEXT,IVAR
         ELSE IF (K .EQ. 2) THEN
            WRITE(TEXT,9901) CTEXT,IVAR
         ELSE IF (K .EQ. 3) THEN
            WRITE(TEXT,9902) CTEXT,IVAR
         ENDIF
C
         RETURN
 9900    FORMAT(A,' = ',I1)
 9901    FORMAT(A,' = ',I4)
 9902    FORMAT(A,' = ',I2)
         END
C*DECK U38
C*CALL PROCESS
         SUBROUTINE WHTEXT(CTEXT,TEXT)
C        -----------------------------
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*) TEXT,CTEXT
C
         WRITE(TEXT,9900) CTEXT
C
         RETURN
 9900    FORMAT(A)
         END
C*DECK U40
C*CALL PROCESS
         SUBROUTINE RESETR(PA,KDIM,PVALUE)
C        #################################
C
C U.40   RESET REAL ARRAY TO SPECIFIED VALUE
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PA(KDIM)
C
         DO 1 J1=1,KDIM
C
         PA(J1) = PVALUE
C
    1    CONTINUE
C
         RETURN
         END
C*DECK U41
C*CALL PROCESS
         SUBROUTINE RESETI(KA,KDIM,KVALUE)
C        #################################
C
C U.41   RESET INTEGER ARRAY TO SPECIFIED VALUE
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     I   KA(KDIM)
C
         DO 1 J1=1,KDIM
C
         KA(J1) = KVALUE
C
    1    CONTINUE
C
         RETURN
         END
C*DECK U42
C*CALL PROCESS
         SUBROUTINE RESETH(KA,KDIM,KVALUE)
C        #################################
C
C U.42   RESET CHARACTER-TYPE ARRAY TO SPECIFIED VALUE
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         CHARACTER*(*)
     H   KA(KDIM), KVALUE
C
         DO 1 J1=1,KDIM
C
         KA(J1) = KVALUE
C
    1    CONTINUE
C
         RETURN
         END
C*DECK U43
C*CALL PROCESS
         SUBROUTINE RESETL(KLA,KDIM,KLVAL)
C        #################################
C
C U.43   RESET LOGICAL ARRAY TO SPECIFIED VALUE
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
         INCLUDE 'DECLAR.inc'
         LOGICAL
     L     KLA(KDIM),
     L     KLVAL
C
         DO 1 J1=1,KDIM
C
         KLA(J1) = KLVAL
C
    1    CONTINUE
C
         RETURN
         END
C*DECK U44
C*CALL PROCESS
         SUBROUTINE RESETC(PC,N,PCX)
C        #################################
C
C U.44   RESET COMPLEX ARRAY TO SPECIFIED VALUE
C
C---*----*----*----*----*----*----*----*----*----*----*----*----*----*
C
C RESET ELEMENTS OF COMPLEX ARRAY PC TO PCX.
C
         INCLUDE 'DECLAR.inc'
         COMPLEX*16    PC, PCX
         DIMENSION
     C   PC(N)
C
         DO 1 J1=1,N
            PC(J1) = PCX
    1    CONTINUE
C
         RETURN
         END
C*DECK U49
C*CALL PROCESS
         SUBROUTINE SCOPYR(RF, N, X, NX, Y, NY) 
C        --------------------------------------
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C
C Y = Y + RF * (X - Y)
C X IS INCREMENTED BY NX
C Y IS INCREMENTED BY NY
C
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   X(N*NX), Y(N*NY)
C
         IF (N .LE. 0) RETURN
C
         Y(1) = Y(1) + RF*(X(1)-Y(1))
C
         IF (N .EQ. 1) RETURN
C
         IF (NX.LT.0 .OR. NY.LT.0) THEN
            PRINT*,'NEGATIVE INCREMENTS FORBIDDEN'
            STOP
         ENDIF
C
         NM1 = N - 1
C
         DO 1 J=1,NM1 
            Y(J*NY+1) = Y(J*NY+1) + RF*(X(J*NX+1)-Y(J*NY+1))
   1     CONTINUE
C
         RETURN
         END
C*DECK U50
C*CALL PROCESS
         SUBROUTINE RUNTIM
C        -----------------
C
C  UPDATE CPU TIME (SECS) AND PRINT IT
C
         INCLUDE 'DECLAR.inc'
         COMMON /COMTIM/ STIME
         DIMENSION TARRAY(2)
         REAL*4 TARRAY
C
         CPTIME = 0.
C
C+DATE   IF (COMPTYP .EQ. 'CRAY') THEN
C+DATE     CALL SECOND(CPTIME)
C+DATE     CPTIME = CPTIME - STIME
C+DATE     WRITE (6,9900) CPTIME
C+DATE   ELSE IF (COMPTYP .EQ. 'SUN' .OR. COMPTYP .EQ. 'SG') THEN
C+DATE     ZETIME = ETIME(TARRAY) 
C+DATE     CPTIME = TARRAY(1)
C+DATE     WRITE (6,9900) CPTIME
C+DATE   ENDIF
C
         RETURN
C
 9900    FORMAT(/,1X,'CPU TIME USED SO FAR =',1PE14.6,' SECS')
C
         END

C*DECK MAT7
C*CALL PROCESS
         SUBROUTINE VZERO(PV,N)
C        ----------------------
C
C  SET REAL ARRAY PV TO ZERO.
C
         INCLUDE 'DECLAR.inc'
         DIMENSION
     R   PV(N)
C
         DO 1 J=1,N
            PV(J) = 0.
    1    CONTINUE
C
         RETURN
         END
C*DECK MAT8
C*CALL PROCESS
         SUBROUTINE CVZERO(CV,N)
C        -------------------------
C
C SET COMPLEX ARRAY CV TO 0.
C
         INCLUDE 'DECLAR.inc'
         COMPLEX    CV
         DIMENSION
     C   CV(N)
C
         DO 1 J=1,N
            CV(J) = (0.,0.)
    1    CONTINUE
C
         RETURN
         END
C*DECK MAT9
C*CALL PROCESS
         SUBROUTINE ACOPY(N,X,NX,Y,NY)
C        -----------------------------
C
C COPIES REAL ARRAY X INTO REAL ARRAY Y WITH INTERMEDIATE
C STORAGE IN ZX. X IS INCREMENTED BY NX AND Y BY NY.
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
C
         DIMENSION
     R   X(*),   Y(*),   ZX(NPMGS*NTP2+2*NPISO)
C
         IF (N .GT. NPMGS*NTP2+2*NPISO) THEN
C
            PRINT*,'N TOO LARGE FOR SUBROUTINE ACOPY'
            STOP
C
         ENDIF
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
            ZX(J) = X(J*NX+1) 
   1     CONTINUE
C
         DO 2 J=1,NM1
            Y(J*NY+1) = ZX(J)
   2     CONTINUE
C
         RETURN
         END
C*DECK MAT10
C*CALL PROCESS
         SUBROUTINE ICOPY(N,IX,NX,IY,NY)
C        -------------------------------
C
C COPIES INTEGER ARRAY IX INTO INTEGER ARRAY IY WITH INTERMEDIATE
C STORAGE IN JX. IX IS INCREMENTED BY NX AND IY BY NY.
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
C
         DIMENSION
     I   JX(NPMGS*NTP2+2*NPISO),   IX(*),   IY(*)
C
         IF (N .GT. NPMGS*NTP2+2*NPISO) THEN
C
            PRINT*,'N TOO LARGE FOR SUBROUTINE ICOPY'
            STOP
C
         ENDIF
C
         IF (N .LE. 0) RETURN
C
         IY(1) = IX(1)
C
         IF (N .EQ. 1) RETURN
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            JX(J) = IX(J*NX+1)
   1     CONTINUE
C
         DO 2 J=1,NM1
            IY(J*NY+1) = JX(J)
   2     CONTINUE
C
         RETURN
         END
C*DECK MAT11
C*CALL PROCESS
         SUBROUTINE CCOPY(N,CX,NX,CY,NY)
C        -------------------------------
C
C COPIES COMPLEX ARRAY CX INTO COMPLEX ARRAY CY.
C CX IS INCREMENTED BY NX AND CY BY NY.
C
         INCLUDE 'DECLAR.inc'
         COMPLEX    CX, CY
         DIMENSION
     C   CX(N*NX),  CY(N*NY)
C
         IF (N .LE. 0) RETURN
C
         CY(1) = CX(1)
C
         IF (N .EQ. 1) RETURN
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            CY(J*NY+1) = CX(J*NX+1)
    1    CONTINUE
C
         RETURN
         END
C*DECK CRAY01
C*CALL PROCESS
         FUNCTION ISMAX(N,PV,NX)
C        -----------------------
C
C  FIND ELEMENT WITH MINIMUM VALUE IN REAL ARRAY PV
C  PV IS INCREMENTED BY NX
C
         INCLUDE 'DECLAR.inc'
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
            IF (PV(I) .GT. PV(ISM)) ISM = I         
    1    CONTINUE   
    2    CONTINUE   
C
         ISMAX = ISM
C
         RETURN
         END
C*DECK CRAY02
C*CALL PROCESS
         FUNCTION ISMIN(N,PV,NX)
C        -----------------------
C
C  FIND ELEMENT WITH MAXIMUM VALUE IN REAL ARRAY PV.
C  PV IS INCREMENTED BY NX
C
         INCLUDE 'DECLAR.inc'
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
            IF (PV(I) .LT. PV(ISM)) ISM = I         
    1    CONTINUE   
    2    CONTINUE   
C
         ISMIN = ISM
C
         RETURN
         END
C*DECK CRAY03
C*CALL PROCESS
         FUNCTION ISAMIN(N,PV,NX)
C        ------------------------
C
C  FIND ELEMENT WITH MINIMUM ABSOLUTE VALUE IN REAL ARRAY PV.
C  PV IS INCREMENTED BY NX
C
         INCLUDE 'DECLAR.inc'
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
            IF (ABS(PV(I)) .LT. ABS(PV(ISM))) ISM = I         
    1    CONTINUE   
    2    CONTINUE   
C
         ISAMIN = ISM
C
         RETURN
         END
C*DECK CRAY04
C*CALL PROCESS
         FUNCTION ISRCHFGT(N,PV,NX,TARGET)
C        ---------------------------------
C
C  FIND FIRST ELEMENT IN REAL ARRAY PV WHICH IS GREATER THAN TARGET
C  PV IS INCREMENTED BY NX
C
         INCLUDE 'DECLAR.inc'
         DIMENSION 
     R   PV(N*NX)
C
         I = 1
         IF (NX .LT. 0) STOP 'NX<0'
C
         DO 1 J=1,N
            IF (PV(I) .GT. TARGET) GOTO 2
            I = I + NX
    1    CONTINUE
    2    ISRCHFGT = I
C
         RETURN
         END
C*DECK CRAY05
C*CALL PROCESS
         FUNCTION ISRCHFGE(N,PV,NX,TARGET)
C        ---------------------------------
C
C  FIND FIRST ELEMENT IN REAL ARRAY PV WHICH IS GREATER OR EQUAL
C  THAN TARGET. PV IS INCREMENTED BY NX
C
         INCLUDE 'DECLAR.inc'
         DIMENSION 
     R   PV(N*NX)
C
         I = 1
         IF (NX .LT. 0) STOP 'NX<0'
C
         DO 1 J=1,N
            IF (PV(I) .GE. TARGET) GOTO 2
            I = I + NX
    1    CONTINUE
    2    ISRCHFGE = I
C
         RETURN
         END
C*DECK CRAY06
C*CALL PROCESS
         SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C        --------------------------------------
C
         INCLUDE 'DECLAR.inc'
         REAL SX(*), SY(*)
C
         IF (N.LE.0 .OR. SA.EQ.0.) RETURN
         IX = 1
         IY = 1
         IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
         IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
         DO 1 I = 1,N
           SY(IY) = SY(IY) + SA*SX(IX)
           IX = IX + INCX
           IY = IY + INCY
   1     CONTINUE
         RETURN
         END
C*DECK CRAY07
C*CALL PROCESS
         SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
C        -----------------------------------
C
         INCLUDE 'DECLAR.inc'
         DIMENSION SX(*),SY(*)
C
         IF (N .LE. 0) RETURN
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX+1
         IF (INCY.LT.0) IY = (-N+1)*INCY+1
         DO 1 I=1,N
            SY(IY) = SX(IX)
            IX = IX+INCX
            IY = IY+INCY
   1     CONTINUE
C
         RETURN
         END
C*DECK CRAY08
C*CALL PROCESS
         FUNCTION SDOT(N,SX,INCX,SY,INCY)
C        --------------------------------
C
         INCLUDE 'DECLAR.inc'
         DIMENSION SX(*), SY(*)
C
         SDOT = 0.
         IF (N .LE. 0) RETURN
         IX = 1
         IY = 1
         IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
         IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
         DO 1 I = 1,N
           SDOT = SDOT + SX(IX)*SY(IY)
           IX = IX + INCX
           IY = IY + INCY
   1     CONTINUE
         RETURN
         END
C*DECK CRAY09
C*CALL PROCESS
         SUBROUTINE SSCAL(N,SA,SX,INCX)
C        ------------------------------
C
         INCLUDE 'DECLAR.inc'
         DIMENSION SX(*)
C
         IF (N.LE.0) RETURN
         IX = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX+1
         DO 1 I=1,N
            SX(IX) = SA*SX(IX)
            IX = IX + INCX
 1       CONTINUE
         RETURN
         END
C*DECK CRAY10
C*CALL PROCESS
         FUNCTION SSUM(N,PV,NX)
C        ----------------------
C
C  SUMS ALL ELEMENTS OF REAL ARRAY PV
C
         INCLUDE 'DECLAR.inc'
         DIMENSION 
     R   PV(N*NX)
C
         SS = 0.
C
         IF (N .LE. 0) GOTO 2
C
         I = 1
         IF (NX .LT. 0) STOP 'NX<0'
         SS = PV(I)
C
         IF (N .EQ. 1) GOTO 2
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            I = I + NX
            SS = SS + PV(I)
    1    CONTINUE
    2    CONTINUE
C   
         SSUM = SS
C
         RETURN
         END
C*DECK CRAY11
C*CALL PROCESS
         FUNCTION ISSUM(N,IV,NX)
C        -----------------------
C
C  SUMS ALL ELEMENTS OF INTEGER ARRAY IC
C
         INCLUDE 'DECLAR.inc'
         DIMENSION 
     R   IV(N*NX)
C
         IS = 0
C
         IF (N .LE. 0) GOTO 2
C
         I = 1
         IF (NX .LT. 0) STOP 'NX<0'
         IS = IV(I)
C
         IF (N .EQ. 1) GOTO 2
C
         NM1 = N - 1
C
         DO 1 J=1,NM1
            I = I + NX
            IS = IS + IV(I)
    1    CONTINUE
    2    CONTINUE
C
         ISSUM = IS
C
         RETURN
         END
         SUBROUTINE C06FAE(X,N,WORK,IFAIL)
C        ---------------------------------
         DIMENSION X(N), WORK(N)
C
         PRINT*,'YOU TRIED TO RUN FFT''S FOR MARS EQ''S WITHOUT'
         PRINT*,'LOADING THE NAG ROUTINE C06FAE. LOAD THE NAG'
         PRINT*,'LIBRARY AND REMOVE DUMMY C06FAE ROUTINE IN'
         PRINT*,'THE CHEASE SOURCE'
         STOP
         RETURN
         END
C---------------------------------------------------------------------
         SUBROUTINE GIJREA(KPSI,PS)
C        ##########################
C
C                                        AUTHORS:
C                                        H. LUTJENS,  CRPP-EPFL
C                                        A. BONDESON, CRPP-EPFL
***********************************************************************
*                                                                     *
* COMPUTE EQ'S FOR MARS IN REAL SPACE [1], TABLE 3                    *
* AT GAUSSIAN QUADRATURE POINTS ALONG CONSTANT POLOIDAL FLUX SURFACES *
*                                                                     *
***********************************************************************
C
         INCLUDE 'DECLAR.inc'
         INCLUDE 'COMDIM.inc'
         INCLUDE 'COMCON.inc'
         INCLUDE 'COMETA.inc'
         INCLUDE 'COMIOD.inc'
         INCLUDE 'COMISO.inc'
         INCLUDE 'COMMAP.inc'
         INCLUDE 'COMNUM.inc'
         INCLUDE 'COMPHY.inc'
         INCLUDE 'COMRCH.inc'
         INCLUDE 'COMSOL.inc'
         INCLUDE 'COMSUR.inc'
         INCLUDE 'COMESH.inc'
C
         INTEGER NPISOC
         PARAMETER (NPISOC = NPMGS*NTP1)
         DIMENSION
     R     ICHIISO(NP2CHI),
     R     ZCHI(NP2CHI), ZD2FUN(NPISOC), ZWORK(NPISOC,3),
     R     ZA(NP2CHI), ZB(NP2CHI), ZC(NP2CHI), ZD(NP2CHI),
     R     ZRJA(NP2CHI),ZG11L(NP2CHI),ZG22L(NP2CHI),
     R     ZG12L(NP2CHI),ZG33L(NP2CHI),ZRR(NP2CHI),ZRZ(NP2CHI),
     R     ZRJAI(NPISOC),ZG11LI(NPISOC),ZG22LI(NPISOC),
     R     ZG12LI(NPISOC),ZG33LI(NPISOC),ZRRI(NPISOC),ZRZI(NPISOC),
     R     ZRDCDZ(NP2CHI),ZRDSDZ(NP2CHI),ZRBZ(NP2CHI),
     R     ZRDCDZI(NPISOC),ZRDSDZI(NPISOC),ZRBZI(NPISOC)
C
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C
C
C        abs(dpsi/ds)
C  
         ZDPSIS = 2. * PS * CPSRF
         CALL GCHI(KPSI)
C 
C  1. QUANTITIES ON GAUSSIAN INTEGRATION POINTS
C
         DO 1 J1=1,NMGAUS*NT1
C
         ZR = RRISO(J1,KPSI)
C
C ZJAC = 1/[GRAD-S * GRAD-CHI x GRAD-PHI]
C ZGRADS = |GRAD-S|
C ZBSCHI = BETA_{S,CHI}
C ZGCHI2 = |GRAD-CHI|**2
C
         ZJAC   = ZDPSIS * CP(KPSI) * ZR**NER *
     *            GPISO(J1,KPSI)**NEGP
         ZGRADS = GPISO(J1,KPSI) / ZDPSIS
         ZBSCHI = BCHISO(J1)
         ZGCHI2 = (ZBSCHI * ZGRADS)**2 + (ZR / (ZJAC * ZGRADS))**2
C
C SEE DEFINITIONS OF EQL'S IN H.LUETJENS ET AL., COMPUTER PHYSICS
C COMMUNICATIONS 69, 287 (1992)
C
        ZRRI  (J1) = ZR
        ZRZI  (J1) = RZISO(J1,KPSI)
        ZRJAI (J1) = ZJAC
        ZG11LI(J1) = (ZJAC*ZBSCHI*ZGRADS/ZR)**2 + 1./(ZGRADS**2)
        ZG22LI(J1) =  (ZJAC*ZGRADS/ZR)**2
        ZG33LI(J1) = ZR**2
        ZG12LI(J1) = - ZBSCHI * (ZJAC*ZGRADS / ZR)**2

        ZDPDR = DPRISO(J1,KPSI)
        ZDPDZ = DPZISO(J1,KPSI)
        ZDSDR = ZDPDR/ZDPSIS
        ZDSDZ = ZDPDZ/ZDPSIS
        ZDCDZ = ZDSDZ*ZBSCHI+ZDSDR*ZR/(ZJAC*ZGRADS**2)
        ZBZ =  ZDPDR/ZR

        ZRDCDZI(J1) = ZDCDZ
        ZRDSDZI(J1) = ZDSDZ
        ZRBZI(J1) = ZBZ
 1    CONTINUE
C
C----*----*----*---*----*----*----*----*----*----*----*----*----*----*-
C
C     2. COMPUTE EQUIDISTANT CHI MESH, KMMAX INTERVALS
C        
         INCHI = 2*NCHI
C
         ZDCHI = 2. * CPI / DFLOAT(INCHI)
         DO I=1,INCHI+1
           ZCHI(I) = DFLOAT(I-1)*ZDCHI
         ENDDO
C
C     3. PREPARE COEFFICIENTS FOR THE CUBIC SPLINE FIT DEPENDING
C     ONLY ON RELATIVE POSITION OF ZCHI(I) WITH RESPECT TO CHIISO
C
         IGCHISO = NMGAUS*NT1
         CALL GCHI(KPSI)
         DO I=1,INCHI
           ICHISO = ISRCHFGE(IGCHISO,CHIISO,1,ZCHI(I)) - 1
C
           IF (ICHISO .LT. 1) THEN
             ICHISO = 1
           ELSE IF (ICHISO .GT. IGCHISO) THEN
             ICHISO = IGCHISO
           ENDIF
           ICHIISO(I) = ICHISO
C
           ZH = CHIISO(ICHISO+1) - CHIISO(ICHISO)
           ZA(I) = (CHIISO(ICHISO+1) - ZCHI(I)) / ZH
           ZB(I) = (ZCHI(I) - CHIISO(ICHISO)) / ZH
           ZC(I) = (ZA(I) + 1) * (ZA(I) - 1) * ZH * 
     *       (CHIISO(ICHISO+1) - ZCHI(I)) / 6.
           ZD(I) = (ZB(I) + 1) * (ZB(I) - 1) * ZH * 
     *       (ZCHI(I) - CHIISO(ICHISO)) / 6.
C
         ENDDO
C
C     3. FOR EACH ARRAY: COMPUTE VALUES ON ZCHI USING A PERIODIC
C     CUBIC SPLINE FIT AND COMPUTE FULL FOURIER TRANSFORM
C
C     EQL
C
         CALL SPLCHI(CHIISO,ZRJAI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRJA,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZRRI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRR,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZRZI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRZ,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG11LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG11L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG22LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG22L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG12LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG12L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZG33LI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZG33L,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZRDCDZI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRDCDZ,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZRDSDZI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRDSDZ,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
         CALL SPLCHI(CHIISO,ZRBZI,IGCHISO,RC2PI,ZD2FUN,
     +     ZWORK,ZRBZ,INCHI,ICHIISO,ZA,ZB,ZC,ZD)
C
C-----------------------------------------------------------------------
C
      IF (KPSI .GT. 1) GOTO 20
C
C        OPEN(NRMAR,FILE='OUTRMAR',STATUS='NEW',FORM='UNFORMATTED')
          OPEN(NRMAR,FILE='OUTRMAR',STATUS='NEW',FORM='FORMATTED')
         REWIND NRMAR
C
C        WRITE(NRMAR) NPSI1, INCHI 
C        WRITE(NRMAR) 1. / ASPCT
         WRITE(NRMAR,*) NPSI1, INCHI, MSMAX 
         WRITE(NRMAR,*) 1./ASPCT,R0EXP0,B0EXP0
 20   CONTINUE 
C
      WRITE(NRMAR,*) PS,PSIISO(KPSI),CPR(KPSI),TMF(KPSI)
      WRITE(NRMAR,*) TTP(KPSI),CPPR(KPSI),ZDPSIS
      WRITE(NRMAR,*) (ZRR(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZRZ(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZRJA(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG11L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG22L(J),J=1,INCHI) 
      WRITE(NRMAR,*) (ZG33L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZG12L(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZRDCDZ(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZRDSDZ(J),J=1,INCHI)
      WRITE(NRMAR,*) (ZRBZ(J),J=1,INCHI)

      RETURN
      END
C

         SUBROUTINE SPLCHI(PCHI,PFUN,KCHI,PPERI,PD2FUN,PWORK,PEQFUN,
     +                      KNFFT,KINDEX,PA,PB,PC,PD)
C        ##############################################################
C
C                                        AUTHOR O. SAUTER, CRPP-EPFL
***********************************************************************
*                                                                     *
* C2SM24 COMPUTE CUBIC SPLINE WITH PERIODIC B.C.
*                                                                     *
***********************************************************************
C
      DIMENSION PCHI(KCHI), PFUN(KCHI), PD2FUN(KCHI), PWORK(KCHI,3),
     +   PEQFUN(KNFFT), KINDEX(KNFFT),
     +   PA(KNFFT), PB(KNFFT), PC(KNFFT), PD(KNFFT)
c.......................................................................
C
C     COMPUTE CUBIC SPLINE OF PFUN USING PA, PB, PC, PD PRECOMPUTED
C
         CALL SPLCY(PCHI,PFUN,KCHI,PPERI,PD2FUN,PWORK(1,1),PWORK(1,2),
     +     PWORK(1,3))
         DO 10 I=1,KNFFT
           K = KINDEX(I)
           PEQFUN(I) = PA(I)*PFUN(K)   + PB(I)*PFUN(K+1) +
     +                  PC(I)*PD2FUN(K) + PD(I)*PD2FUN(K+1)
 10   CONTINUE
C
      RETURN
      END

!*DECK C2SM15
!*CALL PROCESS
      SUBROUTINE OUTNVW
!        #################
!
!                                        AUTHORS:
!                                        H. LUTJENS,  CRPP-EPFL
!                                        A. BONDESON, CRPP-EPFL
!**********************************************************************
!                                                                     *
! C2SM15 COMPUTE EQ'S FOR NOVA-W AND PEST (SEE SECTION 5.4.3 OF       *
!        PUBLICATION                                                  *
!     *
!**********************************************************************
!     
      INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(10)
C     USE globals
C     USE interpol
C     IMPLICIT NONE
!     
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMSOL.inc'
      INCLUDE 'COMNUM.inc'
      INCLUDE 'COMMAP.inc'
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMSUR.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'CUCCCC.inc'
      INCLUDE 'COMESH.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMIOD.inc'
      INTEGER       	::  	M ! <outnvw.f90>
      INTEGER       	::  	L ! <outnvw.f90>
      REAL(RKIND)   	::  	ZRMJ ! <outnvw.f90>
      INTEGER       	::  	ISMAX ! <outnvw.f90>
      INTEGER       	::  	IMX ! <outnvw.f90>
      INTEGER       	::  	ISMIN ! <outnvw.f90>
      INTEGER       	::  	IMN ! <outnvw.f90>
      REAL(RKIND)   	::  	ZZ ! <outnvw.f90>
      REAL(RKIND)   	::  	ZR ! <outnvw.f90>
      INTEGER       	::  	J7 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZBND ! <outnvw.f90>
      INTEGER       	::  	J6 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZDT ! <outnvw.f90>
      INTEGER       	::  	IBND ! <outnvw.f90>
      REAL(RKIND)   	::  	ZTMF ! <outnvw.f90>
      REAL(RKIND)   	::  	ZDQ ! <outnvw.f90>
      REAL(RKIND)   	::  	ZQ ! <outnvw.f90>
      REAL(RKIND)   	::  	ZCPPR ! <outnvw.f90>
      REAL(RKIND)   	::  	ZCPR ! <outnvw.f90>
      REAL(RKIND)   	::  	AXX ! <outnvw.f90>
      INTEGER*4       	::  	NXX ! <outnvw.f90>
      REAL(RKIND)   	::  	SSUM ! <outnvw.f90>
      REAL(RKIND)   	::  	ZJMAG ! <outnvw.f90>
      REAL(RKIND)   	::  	ZJACM0 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZJAC0 ! <outnvw.f90>
      INTEGER       	::  	J5 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZFBP ! <outnvw.f90>
      REAL(RKIND)   	::  	ZFB ! <outnvw.f90>
      REAL(RKIND)   	::  	ZTP ! <outnvw.f90>
      REAL(RKIND)   	::  	ZPSIM ! <outnvw.f90>
      REAL(RKIND)   	::  	ZPSI ! <outnvw.f90>
      REAL(RKIND)   	::  	ZJACM ! <outnvw.f90>
      REAL(RKIND)   	::  	ZZCPM ! <outnvw.f90>
      REAL(RKIND)   	::  	ZRCPM ! <outnvw.f90>
      REAL(RKIND)   	::  	ZJAC ! <outnvw.f90>
      REAL(RKIND)   	::  	ZZCP ! <outnvw.f90>
      REAL(RKIND)   	::  	ZRCP ! <outnvw.f90>
      INTEGER       	::  	J3 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZTETCPM ! <outnvw.f90>
      REAL(RKIND)   	::  	ZSIGCPM ! <outnvw.f90>
      INTEGER       	::  	IT0 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZTETCP ! <outnvw.f90>
      REAL(RKIND)   	::  	ZSIGCP ! <outnvw.f90>
      REAL(RKIND)   	::  	ZD2TET ! <outnvw.f90>
      REAL(RKIND)   	::  	ZC1 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZB1 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZA1 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZD2SIG ! <outnvw.f90>
      INTEGER       	::  	J2 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZTET ! <outnvw.f90>
      INTEGER       	::  	I3 ! <outnvw.f90>
      INTEGER       	::  	J4 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZCHIM ! <outnvw.f90>
      REAL(RKIND)   	::  	ZCHI ! <outnvw.f90>
      INTEGER       	::  	J1 ! <outnvw.f90>
      REAL(RKIND)   	::  	ZDCHI ! <outnvw.f90>
      DIMENSION 
     &     NXX(3),                IT0(NPCHI+3),
     &     ZA1(NTP2),             ZBND(12*NPT),
     &     ZB1(NTP2),             ZC1(NTP2),
     &     AXX(5),                ZCHI(NPCHI+3),
     &     ZCHIM(NPCHI+3),        ZCPR(NPPSNVW),
     &     ZCPPR(NPPSNVW),        ZDQ(NPPSNVW),
     &     ZD2SIG(NTP2),          ZD2TET(NTP2),
     &     ZFB(NPPSNVW),          ZFBP(NPPSNVW),
     &     ZJAC(NPCHI+3,NPPSNVW), ZJAC0(NPCHI+3),
     &     ZJACM(NPCHI+3,NPPSNVW),ZJACM0(NPCHI+3),
     &     ZPSI(NPPSNVW),         ZPSIM(NPPSNVW),
     &     ZQ(NPPSNVW),           ZR(12*NPT),
     &     ZRCP(NPCHI+3,NPPSNVW), ZRCPM(NPCHI+3),
     &     ZSIGCP(NPCHI+3),       ZSIGCPM(NPCHI+3),
     &     ZTET(NTP2+12*NPT),     ZTETCP(NPCHI+3),
     &     ZTETCPM(NPCHI+3),      ZTMF(NPPSNVW),
     &     ZTP(NPPSNVW),          ZZ(12*NPT),
     &     ZZCP(NPCHI+3,NPPSNVW), ZZCPM(NPCHI+3)
!     
!---- *----*----*---*----*----*----*----*----*----*----*----*----*----*-
!     
      IF (NSYM .EQ. 1) THEN
!     
         ZDCHI = CPI / REAL(NTNOVA-1)
!     
      ELSE
!     
         ZDCHI = 2._RKIND * CPI / REAL(NTNOVA)
!     
      ENDIF
!     
      DO 1 J1=1,NTNOVA+3
!     
         ZCHI(J1)  = (J1 - 3._RKIND)  * ZDCHI
         ZCHIM(J1) = (J1 - 3.5_RKIND) * ZDCHI
!     
 1    CONTINUE
!     
      ZCHI(1)  = ZCHI(1)  + 2._RKIND * CPI
      ZCHI(2)  = ZCHI(2)  + 2._RKIND * CPI
      ZCHIM(1) = ZCHIM(1) + 2._RKIND * CPI
      ZCHIM(2) = ZCHIM(2) + 2._RKIND * CPI
      ZCHIM(3) = ZCHIM(3) + 2._RKIND * CPI
!     
      DO 4 J4=2,2*NPSI,2
!     
         I3 = J4 / 2  + 1
!     
         CALL SCOPY(NT2,TETMAP(1,J4),1,ZTET,1)
!     
         DO 2 J2=2,NT2
!     
            IF (ZTET(J2) .LT. ZTET(J2-1)) THEN
!     
               ZTET(J2) = ZTET(J2) + 2._RKIND * CPI * (1._RKIND + 
     &              INT(.5_RKIND * ABS(ZTET(J2) - ZTET(J2-1)) / CPI))       
!     
            ENDIF
!     
    2    CONTINUE
!     
         CALL SPLCY(CHIN(1,J4),SIGMAP(1,J4),NT1,RC2PI,
     &        ZD2SIG,ZA1,ZB1,ZC1)
         CALL SPLCYP(CHIN(1,J4),ZTET,NT1,RC2PI,RC2PI,
     &        ZD2TET,ZA1,ZB1,ZC1)
!     
         ZD2SIG(NT2) = ZD2SIG(1) 
         ZD2TET(NT2) = ZD2TET(1) 
!     
         CALL STCHPS(J4,NTNOVA+3,ZCHI ,ZTET,ZD2TET,SIGMAP(1,J4),
     &        ZD2SIG, ZSIGCP, ZTETCP, IT0)
         CALL STCHPS(J4,NTNOVA+3,ZCHIM,ZTET,ZD2TET,SIGMAP(1,J4),
     &        ZD2SIG,ZSIGCPM,ZTETCPM, IT0)
!     
         DO 3 J3=1,NTNOVA+3
!     
            IF (ZTETCP(J3) .GT.CT(NT1)) ZTETCP(J3)  = ZTETCP(J3)  
     &           - 2._RKIND*CPI
            IF (ZTETCPM(J3).GT.CT(NT1)) ZTETCPM(J3) = ZTETCPM(J3)
     &           - 2._RKIND*CPI
!     
    3    CONTINUE
!     
         CALL JNOVAW(J4,NTNOVA+3,ZTETCP, ZSIGCP,ZRCP(1,I3),ZZCP(1,I3),
     &        ZJAC(1,I3))
         CALL JNOVAW(J4,NTNOVA+3,ZTETCPM,ZSIGCPM,ZRCPM,ZZCPM,
     &        ZJACM(1,I3))
!     
         ZPSI(I3)    = PSIISO(J4)   - CPSRF
         ZPSIM(I3-1) = PSIISO(J4-1) - CPSRF
         ZTP(I3)     = TTP(J4) / TMF(J4)
         ZFB(I3)     = TMF(J4) / QPSI(J4)
         ZFBP(I3)    = (TTP(J4) / (TMF(J4) * QPSI(J4)) - 
     &        TMF(J4) * CDQ(J4) / QPSI(J4)**2)
!     
 4    CONTINUE
!     
      CALL RESETR(ZRCP(1,1),NTNOVA+3,RMAG)
      CALL RESETR(ZZCP(1,1),NTNOVA+3,RZMAG)
!     
      DO 5 J5=1,NTNOVA+3
!     
         ZJAC0(J5)  = FCCCC0(ZJAC(J5,2),ZJAC(J5,3),
     &        ZJAC(J5,4),ZJAC(J5,5),
     &        CS(2),CS(3),CS(4),CS(5),RC0P)
!     
         ZJACM0(J5) = FCCCC0(ZJACM(J5,2),ZJACM(J5,3),
     &        ZJACM(J5,4),ZJACM(J5,5),
     &        CS(2),CS(3),CS(4),CS(5),RC0P)
!     
 5    CONTINUE
!     
      ZJMAG = .5_RKIND*(SSUM(NTNOVA+3,ZJAC0,1) + 
     &     SSUM(NTNOVA+3,ZJACM0,1)) / 
     &     REAL(NTNOVA+3)
!     
      CALL RESETR(ZJAC(1,1), NTNOVA+3,ZJMAG)
      CALL RESETR(ZJACM(1,1),NTNOVA+3,ZJMAG)
!     
      NXX(1) = NTNOVA
      NXX(2) = NPSI1
      NXX(3) = NSYM
      AXX(1) = ZDCHI
      AXX(2) = SQRT(2._RKIND * CPI) / REAL(NPSI,RKIND)
      AXX(3) = RMAG 
      AXX(4) = RZMAG
      AXX(5) = TMF(2*NPSI)
!     
      CALL SCOPY(NPSI,CPR(1) ,2,ZCPR,1)
      CALL SCOPY(NPSI,CPPR(2),2,ZCPPR(2),1)
      CALL SCOPY(NPSI,QPSI(2),2,ZQ(2),1)
      CALL SCOPY(NPSI,CDQ(2),2,ZDQ(2),1)
      CALL SCOPY(NPSI,TMF(1),2,ZTMF(1),1)
!     
      ZCPPR(1) = DPDP0
      ZQ(1)    = Q0
      ZDQ(1)   = DQDP0
      ZTP(1)   = DTTP0 / T0
      ZFB(1)   = FCCCC0(ZFB(2),ZFB(3),ZFB(4),ZFB(5),
     &     CS(2),CS(3),CS(4),CS(5),RC0P)
      ZFBP(1)  = FCCCC0(ZFBP(2),ZFBP(3),ZFBP(4),ZFBP(5),
     &     CS(2),CS(3),CS(4),CS(5),RC0P)
      ZPSI(1)  = - CPSRF
!     
      ZTMF(NPSI1) = TMF(2*NPSI)
!     
      IBND = 12 * NT
!     
      ZDT  = 2._RKIND* CPI / REAL(IBND-1,RKIND)
!     
      DO 6 J6=1,IBND
!     
         ZTET(J6) = (J6 - 1._RKIND) * ZDT
!     
 6    CONTINUE
!     
      CALL BOUND(IBND,ZTET,ZBND)
!     
      DO 7 J7=1,IBND
!     
         ZR(J7) = R0  + ZBND(J7) * COS(ZTET(J7))
         ZZ(J7) = RZ0 + ZBND(J7) * SIN(ZTET(J7))
!     
 7    CONTINUE
!     
      IMN = ISMIN(IBND,ZR,1)
      IMX = ISMAX(IBND,ZR,1)
!     
      ZRMJ = .5_RKIND * (ZR(IMN) + ZR(IMX))
!     
      CALL SSCAL(NPSI1,ZRMJ,ZFB,1)
      CALL SSCAL(NPSI1,ZRMJ,ZFBP,1)
!     
!     -------------------------------------------------------------------
!     
!     THE MODIFIED OPEN STATEMENT AND SOLUTION RENORMALIZATION WERE
!     ADDED 8/2003 BY J. MENARD, PPPL TO MAKE THE OUTPUT HAVE THE 
!     CORRECT PHYSICAL UNITS AND MAKE THE FILE COMPATIBLE WITH DCON
!
      OPEN(INP1,FILE='INP1',FORM='UNFORMATTED')
!MSC   OPEN(INP1,FILE='INP1',FORM='UNFORMATTED',CONVERT="BIG_ENDIAN")
!     
      AXX(3) = AXX(3) * R0EXP
      AXX(4) = AXX(4) * R0EXP
      AXX(5) = AXX(5) * R0EXP
      
      ZCPR   = ZCPR  * R0EXP**0 * B0EXP**2
      ZCPPR  = ZCPPR / R0EXP**2 * B0EXP**1
      ZQ     = ZQ    * R0EXP**0 * B0EXP**0
      ZDQ    = ZDQ   / R0EXP**2 / B0EXP**1
      ZTMF   = ZTMF  * R0EXP**1 * B0EXP**1
      ZTP    = ZTP   / R0EXP**1 * B0EXP**0
      ZFB    = ZFB   * R0EXP**1 * B0EXP**1
      ZFBP   = ZFBP  / R0EXP**1 * B0EXP**0
      ZPSI   = ZPSI  * R0EXP**2 * B0EXP**1
      ZPSIM  = ZPSIM * R0EXP**2 * B0EXP**1
      ZRCP   = ZRCP  * R0EXP**1 * B0EXP**0
      ZZCP   = ZZCP  * R0EXP**1 * B0EXP**0
      ZJACM  = ZJACM * R0EXP**1 / B0EXP**1
      ZJAC   = ZJAC  * R0EXP**1 / B0EXP**1
!     
!     -------------------------------------------------------------------
!     
      WRITE(INP1) (NXX(L),L=1,3)
      WRITE(INP1) (AXX(L),L=1,5)
      WRITE(INP1) (ZCPR(L),L=1,NPSI)
      WRITE(INP1) (ZCPPR(L),L=1,NPSI1)
      WRITE(INP1) (ZQ(L),L=1,NPSI1)
      WRITE(INP1) (ZDQ(L),L=1,NPSI1)
      WRITE(INP1) (ZTMF(L),L=1,NPSI1)
      WRITE(INP1) (ZTP(L),L=1,NPSI1)
      WRITE(INP1) (ZFB(L),L=1,NPSI1)
      WRITE(INP1) (ZFBP(L),L=1,NPSI1)
      WRITE(INP1) (ZPSI(L),L=1,NPSI1)
      WRITE(INP1) (ZPSIM(L),L=1,NPSI)
      WRITE(INP1) ((ZRCP(L,M),L=1,NTNOVA+3),M=1,NPSI1)
      WRITE(INP1) ((ZZCP(L,M),L=1,NTNOVA+3),M=1,NPSI1)
      WRITE(INP1) ((ZJACM(L,M),L=1,NTNOVA+3),M=1,NPSI1)
      WRITE(INP1) ((ZJAC(L,M),L=1,NTNOVA+3),M=1,NPSI1)
!
      CLOSE(INP1,STATUS='KEEP')
!
!CJEM WRITE FORMATTED ASCII FILE FOR OS COMPATIBILITY
!
      OPEN(INP1,FILE='INP1_FORMATTED',FORM='FORMATTED')
      WRITE(INP1,1000) (NXX(L),L=1,3)
      WRITE(INP1,1001) (AXX(L),L=1,5)
      WRITE(INP1,1001) (ZCPR(L),L=1,NPSI)
      WRITE(INP1,1001) (ZCPPR(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZQ(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZDQ(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZTMF(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZTP(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZFB(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZFBP(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZPSI(L),L=1,NPSI1)
      WRITE(INP1,1001) (ZPSIM(L),L=1,NPSI)
      WRITE(INP1,1001) ((ZRCP(L,M),L=1,NTNOVA+3),M=1,NPSI1)
      WRITE(INP1,1001) ((ZZCP(L,M),L=1,NTNOVA+3),M=1,NPSI1)
      WRITE(INP1,1001) ((ZJACM(L,M),L=1,NTNOVA+3),M=1,NPSI1)
      WRITE(INP1,1001) ((ZJAC(L,M),L=1,NTNOVA+3),M=1,NPSI1)
 1000 FORMAT(3I5)
 1001 FORMAT(5E22.15)
      CLOSE(INP1,STATUS='KEEP')
!     
      RETURN
      END
