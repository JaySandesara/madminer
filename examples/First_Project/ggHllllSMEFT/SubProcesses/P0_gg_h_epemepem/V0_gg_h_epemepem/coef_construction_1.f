      SUBROUTINE COEF_CONSTRUCTION_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=64)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=2)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=80, NLOOPGROUPS=9, NCTAMPS=46)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=126)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=47,NLOOPWAVEFUNCS=98)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/I_SO/I_SO
      INTEGER I_LIB
      COMMON/I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/WL/WL,PL

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 3
      CALL VVV3L2P0_1(PL(0,0),W(1,1),GC_10,ZERO,ZERO,PL(0,1),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL VVV3L2P0_1(PL(0,1),W(1,2),GC_10,ZERO,ZERO,PL(0,2),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL VVS4L2P0_1(PL(0,2),W(1,12),GC_426,ZERO,ZERO,PL(0,3),COEFS)
      CALL UPDATE_WL_2_2(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,3))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,3),4,4,1,1,1,47,H)
C     Coefficient construction for loop diagram with ID 4
      CALL VVS4L2P0_1(PL(0,2),W(1,13),GC_426,ZERO,ZERO,PL(0,4),COEFS)
      CALL UPDATE_WL_2_2(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,4))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,4),4,4,1,1,1,48,H)
C     Coefficient construction for loop diagram with ID 5
      CALL VVVS9L2P0_1(PL(0,1),W(1,2),W(1,12),GC_457,ZERO,ZERO,PL(0,5)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,5))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,5),2,4,2,2,1,49,H)
C     Coefficient construction for loop diagram with ID 6
      CALL VVVS9L2P0_1(PL(0,1),W(1,2),W(1,13),GC_457,ZERO,ZERO,PL(0,6)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,6))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,6),2,4,2,2,1,50,H)
C     Coefficient construction for loop diagram with ID 7
      CALL VVV3L2P0_1(PL(0,0),W(1,2),GC_10,ZERO,ZERO,PL(0,7),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,7))
      CALL VVVS9L2P0_1(PL(0,7),W(1,1),W(1,12),GC_457,ZERO,ZERO,PL(0,8)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,8))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,8),2,4,3,2,1,51,H)
C     Coefficient construction for loop diagram with ID 8
      CALL VVVS9L2P0_1(PL(0,7),W(1,1),W(1,13),GC_457,ZERO,ZERO,PL(0,9)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,9))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,9),2,4,3,2,1,52,H)
C     Coefficient construction for loop diagram with ID 9
      CALL VVVV4L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,10)
     $ ,COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,10))
      CALL VVVV8L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,11)
     $ ,COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,11))
      CALL VVVV9L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,12)
     $ ,COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,12))
      CALL VVS4L2P0_1(PL(0,10),W(1,12),GC_426,ZERO,ZERO,PL(0,13),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,10),4,COEFS,4,4,WL(1,0,1,13))
      CALL VVS4L2P0_1(PL(0,11),W(1,12),GC_426,ZERO,ZERO,PL(0,14),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,14))
      CALL VVS4L2P0_1(PL(0,12),W(1,12),GC_426,ZERO,ZERO,PL(0,15),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,12),4,COEFS,4,4,WL(1,0,1,15))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,13),2,4,4,2,1,53,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,14),2,4,4,2,1,54,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,15),2,4,4,2,1,55,H)
C     Coefficient construction for loop diagram with ID 10
      CALL VVS4L2P0_1(PL(0,10),W(1,13),GC_426,ZERO,ZERO,PL(0,16),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,10),4,COEFS,4,4,WL(1,0,1,16))
      CALL VVS4L2P0_1(PL(0,11),W(1,13),GC_426,ZERO,ZERO,PL(0,17),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,17))
      CALL VVS4L2P0_1(PL(0,12),W(1,13),GC_426,ZERO,ZERO,PL(0,18),COEFS)
      CALL UPDATE_WL_0_2(WL(1,0,1,12),4,COEFS,4,4,WL(1,0,1,18))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,16),2,4,4,2,1,56,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,17),2,4,4,2,1,57,H)
      CALL CREATE_LOOP_COEFS(WL(1,0,1,18),2,4,4,2,1,58,H)
C     Coefficient construction for loop diagram with ID 11
      CALL FFV10L1_2(PL(0,0),W(1,1),GC_11,MDL_MT,MDL_WT,PL(0,19),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,19))
      CALL FFV10L1_2(PL(0,19),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,20)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,19),4,COEFS,4,4,WL(1,0,1,20))
      CALL FFS3L1_2(PL(0,20),W(1,14),GC_684,MDL_MT,MDL_WT,PL(0,21)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,21))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,21),3,4,5,1,1,59,H)
C     Coefficient construction for loop diagram with ID 12
      CALL FFV26L1_2(PL(0,0),W(1,1),GC_458,MDL_MT,MDL_WT,PL(0,22)
     $ ,COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,22))
      CALL FFV10L1_2(PL(0,22),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,23)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,22),4,COEFS,4,4,WL(1,0,1,23))
      CALL FFS3L1_2(PL(0,23),W(1,12),GC_684,MDL_MT,MDL_WT,PL(0,24)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,23),4,COEFS,4,4,WL(1,0,1,24))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,24),3,4,5,1,1,60,H)
C     Coefficient construction for loop diagram with ID 13
      CALL FFV26L1_2(PL(0,19),W(1,2),GC_458,MDL_MT,MDL_WT,PL(0,25)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,19),4,COEFS,4,4,WL(1,0,1,25))
      CALL FFS3L1_2(PL(0,25),W(1,12),GC_684,MDL_MT,MDL_WT,PL(0,26)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,25),4,COEFS,4,4,WL(1,0,1,26))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,26),3,4,5,1,1,61,H)
C     Coefficient construction for loop diagram with ID 14
      CALL FFS3L1_2(PL(0,20),W(1,17),GC_684,MDL_MT,MDL_WT,PL(0,27)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,27))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,27),3,4,5,1,1,62,H)
C     Coefficient construction for loop diagram with ID 15
      CALL FFS3L1_2(PL(0,20),W(1,18),GC_684,MDL_MT,MDL_WT,PL(0,28)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,28))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,28),3,4,5,1,1,63,H)
C     Coefficient construction for loop diagram with ID 16
      CALL FFS3L1_2(PL(0,20),W(1,19),GC_684,MDL_MT,MDL_WT,PL(0,29)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,29))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,29),3,4,5,1,1,64,H)
C     Coefficient construction for loop diagram with ID 17
      CALL FFS3L1_2(PL(0,20),W(1,12),GC_691,MDL_MT,MDL_WT,PL(0,30)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,30))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,30),3,4,5,1,1,65,H)
C     Coefficient construction for loop diagram with ID 18
      CALL FFS3L1_2(PL(0,20),W(1,20),GC_684,MDL_MT,MDL_WT,PL(0,31)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,31))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,31),3,4,5,1,1,66,H)
C     Coefficient construction for loop diagram with ID 19
      CALL FFS3L1_2(PL(0,20),W(1,12),GC_684,MDL_MT,MDL_WT,PL(0,32)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,32))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,32),3,4,5,1,1,67,H)
C     Coefficient construction for loop diagram with ID 20
      CALL FFS3L1_2(PL(0,20),W(1,21),GC_684,MDL_MT,MDL_WT,PL(0,33)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,33))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,33),3,4,5,1,1,68,H)
C     Coefficient construction for loop diagram with ID 21
      CALL FFS3L1_2(PL(0,20),W(1,23),GC_684,MDL_MT,MDL_WT,PL(0,34)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,34))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,34),3,4,5,1,1,69,H)
C     Coefficient construction for loop diagram with ID 22
      CALL FFS3L1_2(PL(0,20),W(1,25),GC_684,MDL_MT,MDL_WT,PL(0,35)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,35))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,35),3,4,5,1,1,70,H)
C     Coefficient construction for loop diagram with ID 23
      CALL FFS3L1_2(PL(0,20),W(1,27),GC_684,MDL_MT,MDL_WT,PL(0,36)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,36))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,36),3,4,5,1,1,71,H)
C     Coefficient construction for loop diagram with ID 24
      CALL FFS3L1_2(PL(0,20),W(1,29),GC_684,MDL_MT,MDL_WT,PL(0,37)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,37))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,37),3,4,5,1,1,72,H)
C     Coefficient construction for loop diagram with ID 25
      CALL FFS3L1_2(PL(0,20),W(1,30),GC_684,MDL_MT,MDL_WT,PL(0,38)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,38))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,38),3,4,5,1,1,73,H)
C     Coefficient construction for loop diagram with ID 26
      CALL FFS3L1_2(PL(0,23),W(1,13),GC_684,MDL_MT,MDL_WT,PL(0,39)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,23),4,COEFS,4,4,WL(1,0,1,39))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,39),3,4,5,1,1,74,H)
C     Coefficient construction for loop diagram with ID 27
      CALL FFS3L1_2(PL(0,25),W(1,13),GC_684,MDL_MT,MDL_WT,PL(0,40)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,25),4,COEFS,4,4,WL(1,0,1,40))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,40),3,4,5,1,1,75,H)
C     Coefficient construction for loop diagram with ID 28
      CALL FFS3L1_2(PL(0,20),W(1,33),GC_684,MDL_MT,MDL_WT,PL(0,41)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,41))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,41),3,4,5,1,1,76,H)
C     Coefficient construction for loop diagram with ID 29
      CALL FFS3L1_2(PL(0,20),W(1,34),GC_684,MDL_MT,MDL_WT,PL(0,42)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,42))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,42),3,4,5,1,1,77,H)
C     Coefficient construction for loop diagram with ID 30
      CALL FFS3L1_2(PL(0,20),W(1,35),GC_684,MDL_MT,MDL_WT,PL(0,43)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,43))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,43),3,4,5,1,1,78,H)
C     Coefficient construction for loop diagram with ID 31
      CALL FFS3L1_2(PL(0,20),W(1,13),GC_691,MDL_MT,MDL_WT,PL(0,44)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,44))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,44),3,4,5,1,1,79,H)
C     Coefficient construction for loop diagram with ID 32
      CALL FFS3L1_2(PL(0,20),W(1,36),GC_684,MDL_MT,MDL_WT,PL(0,45)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,45))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,45),3,4,5,1,1,80,H)
C     Coefficient construction for loop diagram with ID 33
      CALL FFS3L1_2(PL(0,20),W(1,13),GC_684,MDL_MT,MDL_WT,PL(0,46)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,46))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,46),3,4,5,1,1,81,H)
C     Coefficient construction for loop diagram with ID 34
      CALL FFS3L1_2(PL(0,20),W(1,37),GC_684,MDL_MT,MDL_WT,PL(0,47)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,47))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,47),3,4,5,1,1,82,H)
C     Coefficient construction for loop diagram with ID 35
      CALL FFS3L1_2(PL(0,20),W(1,39),GC_684,MDL_MT,MDL_WT,PL(0,48)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,48))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,48),3,4,5,1,1,83,H)
C     Coefficient construction for loop diagram with ID 36
      CALL FFS3L1_2(PL(0,20),W(1,41),GC_684,MDL_MT,MDL_WT,PL(0,49)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,49))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,49),3,4,5,1,1,84,H)
C     Coefficient construction for loop diagram with ID 37
      CALL FFS3L1_2(PL(0,20),W(1,43),GC_684,MDL_MT,MDL_WT,PL(0,50)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,50))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,50),3,4,5,1,1,85,H)
C     Coefficient construction for loop diagram with ID 38
      CALL FFS3L1_2(PL(0,20),W(1,45),GC_684,MDL_MT,MDL_WT,PL(0,51)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,51))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,51),3,4,5,1,1,86,H)
C     Coefficient construction for loop diagram with ID 39
      CALL FFS3L1_2(PL(0,20),W(1,46),GC_684,MDL_MT,MDL_WT,PL(0,52)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,52))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,52),3,4,5,1,1,87,H)
C     Coefficient construction for loop diagram with ID 40
      CALL FFS3L1_2(PL(0,20),W(1,47),GC_684,MDL_MT,MDL_WT,PL(0,53)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,53))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,53),3,4,5,1,1,88,H)
C     Coefficient construction for loop diagram with ID 41
      CALL FFV10L2_1(PL(0,0),W(1,1),GC_11,MDL_MT,MDL_WT,PL(0,54),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,54))
      CALL FFV10L2_1(PL(0,54),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,55)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,54),4,COEFS,4,4,WL(1,0,1,55))
      CALL FFS3L2_1(PL(0,55),W(1,14),GC_684,MDL_MT,MDL_WT,PL(0,56)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,56))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,56),3,4,5,1,1,89,H)
C     Coefficient construction for loop diagram with ID 42
      CALL FFV26L2_1(PL(0,0),W(1,1),GC_458,MDL_MT,MDL_WT,PL(0,57)
     $ ,COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,57))
      CALL FFV10L2_1(PL(0,57),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,58)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,57),4,COEFS,4,4,WL(1,0,1,58))
      CALL FFS3L2_1(PL(0,58),W(1,12),GC_684,MDL_MT,MDL_WT,PL(0,59)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,58),4,COEFS,4,4,WL(1,0,1,59))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,59),3,4,5,1,1,90,H)
C     Coefficient construction for loop diagram with ID 43
      CALL FFS3L2_1(PL(0,55),W(1,17),GC_684,MDL_MT,MDL_WT,PL(0,60)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,60))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,60),3,4,5,1,1,91,H)
C     Coefficient construction for loop diagram with ID 44
      CALL FFS3L2_1(PL(0,55),W(1,18),GC_684,MDL_MT,MDL_WT,PL(0,61)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,61))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,61),3,4,5,1,1,92,H)
C     Coefficient construction for loop diagram with ID 45
      CALL FFS3L2_1(PL(0,55),W(1,19),GC_684,MDL_MT,MDL_WT,PL(0,62)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,62))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,62),3,4,5,1,1,93,H)
C     Coefficient construction for loop diagram with ID 46
      CALL FFV26L2_1(PL(0,54),W(1,2),GC_458,MDL_MT,MDL_WT,PL(0,63)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,54),4,COEFS,4,4,WL(1,0,1,63))
      CALL FFS3L2_1(PL(0,63),W(1,12),GC_684,MDL_MT,MDL_WT,PL(0,64)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,63),4,COEFS,4,4,WL(1,0,1,64))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,64),3,4,5,1,1,94,H)
C     Coefficient construction for loop diagram with ID 47
      CALL FFS3L2_1(PL(0,55),W(1,20),GC_684,MDL_MT,MDL_WT,PL(0,65)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,65))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,65),3,4,5,1,1,95,H)
C     Coefficient construction for loop diagram with ID 48
      CALL FFS3L2_1(PL(0,55),W(1,12),GC_691,MDL_MT,MDL_WT,PL(0,66)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,66))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,66),3,4,5,1,1,96,H)
C     Coefficient construction for loop diagram with ID 49
      CALL FFS3L2_1(PL(0,55),W(1,12),GC_684,MDL_MT,MDL_WT,PL(0,67)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,67))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,67),3,4,5,1,1,97,H)
C     Coefficient construction for loop diagram with ID 50
      CALL FFS3L2_1(PL(0,55),W(1,21),GC_684,MDL_MT,MDL_WT,PL(0,68)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,68))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,68),3,4,5,1,1,98,H)
C     Coefficient construction for loop diagram with ID 51
      CALL FFVS19L2_1(PL(0,54),W(1,2),W(1,12),GC_197,MDL_MT,MDL_WT
     $ ,PL(0,69),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,54),4,COEFS,4,4,WL(1,0,1,69))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,69),2,4,6,1,1,99,H)
C     Coefficient construction for loop diagram with ID 52
      CALL FFS3L2_1(PL(0,55),W(1,23),GC_684,MDL_MT,MDL_WT,PL(0,70)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,70))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,70),3,4,5,1,1,100,H)
C     Coefficient construction for loop diagram with ID 53
      CALL FFS3L2_1(PL(0,55),W(1,25),GC_684,MDL_MT,MDL_WT,PL(0,71)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,71))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,71),3,4,5,1,1,101,H)
C     Coefficient construction for loop diagram with ID 54
      CALL FFS3L2_1(PL(0,55),W(1,27),GC_684,MDL_MT,MDL_WT,PL(0,72)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,72))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,72),3,4,5,1,1,102,H)
C     Coefficient construction for loop diagram with ID 55
      CALL FFS3L2_1(PL(0,55),W(1,29),GC_684,MDL_MT,MDL_WT,PL(0,73)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,73))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,73),3,4,5,1,1,103,H)
C     Coefficient construction for loop diagram with ID 56
      CALL FFS3L2_1(PL(0,55),W(1,30),GC_684,MDL_MT,MDL_WT,PL(0,74)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,74))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,74),3,4,5,1,1,104,H)
C     Coefficient construction for loop diagram with ID 57
      CALL FFS3L2_1(PL(0,58),W(1,13),GC_684,MDL_MT,MDL_WT,PL(0,75)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,58),4,COEFS,4,4,WL(1,0,1,75))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,75),3,4,5,1,1,105,H)
C     Coefficient construction for loop diagram with ID 58
      CALL FFS3L2_1(PL(0,55),W(1,33),GC_684,MDL_MT,MDL_WT,PL(0,76)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,76))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,76),3,4,5,1,1,106,H)
C     Coefficient construction for loop diagram with ID 59
      CALL FFS3L2_1(PL(0,55),W(1,34),GC_684,MDL_MT,MDL_WT,PL(0,77)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,77))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,77),3,4,5,1,1,107,H)
C     Coefficient construction for loop diagram with ID 60
      CALL FFS3L2_1(PL(0,55),W(1,35),GC_684,MDL_MT,MDL_WT,PL(0,78)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,78))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,78),3,4,5,1,1,108,H)
C     Coefficient construction for loop diagram with ID 61
      CALL FFS3L2_1(PL(0,63),W(1,13),GC_684,MDL_MT,MDL_WT,PL(0,79)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,63),4,COEFS,4,4,WL(1,0,1,79))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,79),3,4,5,1,1,109,H)
C     Coefficient construction for loop diagram with ID 62
      CALL FFS3L2_1(PL(0,55),W(1,36),GC_684,MDL_MT,MDL_WT,PL(0,80)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,80))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,80),3,4,5,1,1,110,H)
C     Coefficient construction for loop diagram with ID 63
      CALL FFS3L2_1(PL(0,55),W(1,13),GC_691,MDL_MT,MDL_WT,PL(0,81)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,81))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,81),3,4,5,1,1,111,H)
C     Coefficient construction for loop diagram with ID 64
      CALL FFS3L2_1(PL(0,55),W(1,13),GC_684,MDL_MT,MDL_WT,PL(0,82)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,82))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,82),3,4,5,1,1,112,H)
C     Coefficient construction for loop diagram with ID 65
      CALL FFS3L2_1(PL(0,55),W(1,37),GC_684,MDL_MT,MDL_WT,PL(0,83)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,83))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,83),3,4,5,1,1,113,H)
C     Coefficient construction for loop diagram with ID 66
      CALL FFVS19L2_1(PL(0,54),W(1,2),W(1,13),GC_197,MDL_MT,MDL_WT
     $ ,PL(0,84),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,54),4,COEFS,4,4,WL(1,0,1,84))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,84),2,4,6,1,1,114,H)
C     Coefficient construction for loop diagram with ID 67
      CALL FFS3L2_1(PL(0,55),W(1,39),GC_684,MDL_MT,MDL_WT,PL(0,85)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,85))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,85),3,4,5,1,1,115,H)
C     Coefficient construction for loop diagram with ID 68
      CALL FFS3L2_1(PL(0,55),W(1,41),GC_684,MDL_MT,MDL_WT,PL(0,86)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,86))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,86),3,4,5,1,1,116,H)
C     Coefficient construction for loop diagram with ID 69
      CALL FFS3L2_1(PL(0,55),W(1,43),GC_684,MDL_MT,MDL_WT,PL(0,87)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,87))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,87),3,4,5,1,1,117,H)
C     Coefficient construction for loop diagram with ID 70
      CALL FFS3L2_1(PL(0,55),W(1,45),GC_684,MDL_MT,MDL_WT,PL(0,88)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,88))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,88),3,4,5,1,1,118,H)
C     Coefficient construction for loop diagram with ID 71
      CALL FFS3L2_1(PL(0,55),W(1,46),GC_684,MDL_MT,MDL_WT,PL(0,89)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,89))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,89),3,4,5,1,1,119,H)
C     Coefficient construction for loop diagram with ID 72
      CALL FFS3L2_1(PL(0,55),W(1,47),GC_684,MDL_MT,MDL_WT,PL(0,90)
     $ ,COEFS)
      CALL UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,90))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,90),3,4,5,1,1,120,H)
C     Coefficient construction for loop diagram with ID 73
      CALL FFV10L2_1(PL(0,0),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,91),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,91))
      CALL FFVS19L2_1(PL(0,91),W(1,1),W(1,12),GC_197,MDL_MT,MDL_WT
     $ ,PL(0,92),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,91),4,COEFS,4,4,WL(1,0,1,92))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,92),2,4,7,1,1,121,H)
C     Coefficient construction for loop diagram with ID 74
      CALL FFVS19L2_1(PL(0,91),W(1,1),W(1,13),GC_197,MDL_MT,MDL_WT
     $ ,PL(0,93),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,91),4,COEFS,4,4,WL(1,0,1,93))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,93),2,4,7,1,1,122,H)
C     Coefficient construction for loop diagram with ID 75
      CALL FFVV32L2_1(PL(0,0),W(1,1),W(1,2),GC_460,MDL_MT,MDL_WT,PL(0
     $ ,94),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,94))
      CALL FFS3L2_1(PL(0,94),W(1,12),GC_684,MDL_MT,MDL_WT,PL(0,95)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,94),4,COEFS,4,4,WL(1,0,1,95))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,95),2,4,8,1,1,123,H)
C     Coefficient construction for loop diagram with ID 76
      CALL FFS3L2_1(PL(0,94),W(1,13),GC_684,MDL_MT,MDL_WT,PL(0,96)
     $ ,COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,94),4,COEFS,4,4,WL(1,0,1,96))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,96),2,4,8,1,1,124,H)
C     Coefficient construction for loop diagram with ID 77
      CALL FFVVS63L2_1(PL(0,0),W(1,1),W(1,2),W(1,12),GC_203,MDL_MT
     $ ,MDL_WT,PL(0,97),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,97))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,97),1,4,9,1,1,125,H)
C     Coefficient construction for loop diagram with ID 78
      CALL FFVVS63L2_1(PL(0,0),W(1,1),W(1,2),W(1,13),GC_203,MDL_MT
     $ ,MDL_WT,PL(0,98),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,98))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,98),1,4,9,1,1,126,H)
C     At this point, all loop coefficients needed for (QCD=2), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

