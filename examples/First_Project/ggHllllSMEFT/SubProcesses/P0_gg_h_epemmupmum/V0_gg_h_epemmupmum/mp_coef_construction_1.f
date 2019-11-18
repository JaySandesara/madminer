      SUBROUTINE MP_COEF_CONSTRUCTION_1(P,NHEL,H,IC)
C     
      USE POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=64)

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=1)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=40, NLOOPGROUPS=9, NCTAMPS=23)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=63)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=27,NLOOPWAVEFUNCS=58)
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 2
      CALL MP_VVV3L2P0_1(PL(0,0),W(1,1),GC_10,ZERO,ZERO,PL(0,1),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL MP_VVV3L2P0_1(PL(0,1),W(1,2),GC_10,ZERO,ZERO,PL(0,2),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL MP_VVS4L2P0_1(PL(0,2),W(1,10),GC_426,ZERO,ZERO,PL(0,3)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_2(WL(1,0,1,2),4,COEFS,4,4,WL(1,0,1,3))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,3),4,4,1,1,1,24,H)
C     Coefficient construction for loop diagram with ID 3
      CALL MP_VVVS9L2P0_1(PL(0,1),W(1,2),W(1,10),GC_457,ZERO,ZERO,PL(0
     $ ,4),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,4))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,4),2,4,2,2,1,25,H)
C     Coefficient construction for loop diagram with ID 4
      CALL MP_VVV3L2P0_1(PL(0,0),W(1,2),GC_10,ZERO,ZERO,PL(0,5),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,5))
      CALL MP_VVVS9L2P0_1(PL(0,5),W(1,1),W(1,10),GC_457,ZERO,ZERO,PL(0
     $ ,6),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,6))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,6),2,4,3,2,1,26,H)
C     Coefficient construction for loop diagram with ID 5
      CALL MP_VVVV4L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,7)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,7))
      CALL MP_VVVV8L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,8)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,8))
      CALL MP_VVVV9L2P0_1(PL(0,0),W(1,1),W(1,2),GC_12,ZERO,ZERO,PL(0,9)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,9))
      CALL MP_VVS4L2P0_1(PL(0,7),W(1,10),GC_426,ZERO,ZERO,PL(0,10)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_2(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,10))
      CALL MP_VVS4L2P0_1(PL(0,8),W(1,10),GC_426,ZERO,ZERO,PL(0,11)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_2(WL(1,0,1,8),4,COEFS,4,4,WL(1,0,1,11))
      CALL MP_VVS4L2P0_1(PL(0,9),W(1,10),GC_426,ZERO,ZERO,PL(0,12)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_2(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,12))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,10),2,4,4,2,1,27,H)
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,11),2,4,4,2,1,28,H)
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,12),2,4,4,2,1,29,H)
C     Coefficient construction for loop diagram with ID 6
      CALL MP_FFV10L1_2(PL(0,0),W(1,1),GC_11,MDL_MT,MDL_WT,PL(0,13)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,13))
      CALL MP_FFV10L1_2(PL(0,13),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,14)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,14))
      CALL MP_FFS3L1_2(PL(0,14),W(1,11),GC_684,MDL_MT,MDL_WT,PL(0,15)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,15))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,15),3,4,5,1,1,30,H)
C     Coefficient construction for loop diagram with ID 7
      CALL MP_FFV26L1_2(PL(0,0),W(1,1),GC_458,MDL_MT,MDL_WT,PL(0,16)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,16))
      CALL MP_FFV10L1_2(PL(0,16),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,17)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,16),4,COEFS,4,4,WL(1,0,1,17))
      CALL MP_FFS3L1_2(PL(0,17),W(1,10),GC_684,MDL_MT,MDL_WT,PL(0,18)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,17),4,COEFS,4,4,WL(1,0,1,18))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,18),3,4,5,1,1,31,H)
C     Coefficient construction for loop diagram with ID 8
      CALL MP_FFV26L1_2(PL(0,13),W(1,2),GC_458,MDL_MT,MDL_WT,PL(0,19)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,19))
      CALL MP_FFS3L1_2(PL(0,19),W(1,10),GC_684,MDL_MT,MDL_WT,PL(0,20)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,19),4,COEFS,4,4,WL(1,0,1,20))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,20),3,4,5,1,1,32,H)
C     Coefficient construction for loop diagram with ID 9
      CALL MP_FFS3L1_2(PL(0,14),W(1,14),GC_684,MDL_MT,MDL_WT,PL(0,21)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,21))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,21),3,4,5,1,1,33,H)
C     Coefficient construction for loop diagram with ID 10
      CALL MP_FFS3L1_2(PL(0,14),W(1,15),GC_684,MDL_MT,MDL_WT,PL(0,22)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,22))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,22),3,4,5,1,1,34,H)
C     Coefficient construction for loop diagram with ID 11
      CALL MP_FFS3L1_2(PL(0,14),W(1,16),GC_684,MDL_MT,MDL_WT,PL(0,23)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,23))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,23),3,4,5,1,1,35,H)
C     Coefficient construction for loop diagram with ID 12
      CALL MP_FFS3L1_2(PL(0,14),W(1,10),GC_691,MDL_MT,MDL_WT,PL(0,24)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,24))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,24),3,4,5,1,1,36,H)
C     Coefficient construction for loop diagram with ID 13
      CALL MP_FFS3L1_2(PL(0,14),W(1,17),GC_684,MDL_MT,MDL_WT,PL(0,25)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,25))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,25),3,4,5,1,1,37,H)
C     Coefficient construction for loop diagram with ID 14
      CALL MP_FFS3L1_2(PL(0,14),W(1,10),GC_684,MDL_MT,MDL_WT,PL(0,26)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,26))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,26),3,4,5,1,1,38,H)
C     Coefficient construction for loop diagram with ID 15
      CALL MP_FFS3L1_2(PL(0,14),W(1,18),GC_684,MDL_MT,MDL_WT,PL(0,27)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,27))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,27),3,4,5,1,1,39,H)
C     Coefficient construction for loop diagram with ID 16
      CALL MP_FFS3L1_2(PL(0,14),W(1,20),GC_684,MDL_MT,MDL_WT,PL(0,28)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,28))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,28),3,4,5,1,1,40,H)
C     Coefficient construction for loop diagram with ID 17
      CALL MP_FFS3L1_2(PL(0,14),W(1,22),GC_684,MDL_MT,MDL_WT,PL(0,29)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,29))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,29),3,4,5,1,1,41,H)
C     Coefficient construction for loop diagram with ID 18
      CALL MP_FFS3L1_2(PL(0,14),W(1,24),GC_684,MDL_MT,MDL_WT,PL(0,30)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,30))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,30),3,4,5,1,1,42,H)
C     Coefficient construction for loop diagram with ID 19
      CALL MP_FFS3L1_2(PL(0,14),W(1,26),GC_684,MDL_MT,MDL_WT,PL(0,31)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,31))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,31),3,4,5,1,1,43,H)
C     Coefficient construction for loop diagram with ID 20
      CALL MP_FFS3L1_2(PL(0,14),W(1,27),GC_684,MDL_MT,MDL_WT,PL(0,32)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,32))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,32),3,4,5,1,1,44,H)
C     Coefficient construction for loop diagram with ID 21
      CALL MP_FFV10L2_1(PL(0,0),W(1,1),GC_11,MDL_MT,MDL_WT,PL(0,33)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,33))
      CALL MP_FFV10L2_1(PL(0,33),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,34)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,33),4,COEFS,4,4,WL(1,0,1,34))
      CALL MP_FFS3L2_1(PL(0,34),W(1,11),GC_684,MDL_MT,MDL_WT,PL(0,35)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,35))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,35),3,4,5,1,1,45,H)
C     Coefficient construction for loop diagram with ID 22
      CALL MP_FFV26L2_1(PL(0,0),W(1,1),GC_458,MDL_MT,MDL_WT,PL(0,36)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,36))
      CALL MP_FFV10L2_1(PL(0,36),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,37)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,36),4,COEFS,4,4,WL(1,0,1,37))
      CALL MP_FFS3L2_1(PL(0,37),W(1,10),GC_684,MDL_MT,MDL_WT,PL(0,38)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,37),4,COEFS,4,4,WL(1,0,1,38))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,38),3,4,5,1,1,46,H)
C     Coefficient construction for loop diagram with ID 23
      CALL MP_FFS3L2_1(PL(0,34),W(1,14),GC_684,MDL_MT,MDL_WT,PL(0,39)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,39))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,39),3,4,5,1,1,47,H)
C     Coefficient construction for loop diagram with ID 24
      CALL MP_FFS3L2_1(PL(0,34),W(1,15),GC_684,MDL_MT,MDL_WT,PL(0,40)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,40))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,40),3,4,5,1,1,48,H)
C     Coefficient construction for loop diagram with ID 25
      CALL MP_FFS3L2_1(PL(0,34),W(1,16),GC_684,MDL_MT,MDL_WT,PL(0,41)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,41))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,41),3,4,5,1,1,49,H)
C     Coefficient construction for loop diagram with ID 26
      CALL MP_FFV26L2_1(PL(0,33),W(1,2),GC_458,MDL_MT,MDL_WT,PL(0,42)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,33),4,COEFS,4,4,WL(1,0,1,42))
      CALL MP_FFS3L2_1(PL(0,42),W(1,10),GC_684,MDL_MT,MDL_WT,PL(0,43)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,42),4,COEFS,4,4,WL(1,0,1,43))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,43),3,4,5,1,1,50,H)
C     Coefficient construction for loop diagram with ID 27
      CALL MP_FFS3L2_1(PL(0,34),W(1,17),GC_684,MDL_MT,MDL_WT,PL(0,44)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,44))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,44),3,4,5,1,1,51,H)
C     Coefficient construction for loop diagram with ID 28
      CALL MP_FFS3L2_1(PL(0,34),W(1,10),GC_691,MDL_MT,MDL_WT,PL(0,45)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,45))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,45),3,4,5,1,1,52,H)
C     Coefficient construction for loop diagram with ID 29
      CALL MP_FFS3L2_1(PL(0,34),W(1,10),GC_684,MDL_MT,MDL_WT,PL(0,46)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,46))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,46),3,4,5,1,1,53,H)
C     Coefficient construction for loop diagram with ID 30
      CALL MP_FFS3L2_1(PL(0,34),W(1,18),GC_684,MDL_MT,MDL_WT,PL(0,47)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,47))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,47),3,4,5,1,1,54,H)
C     Coefficient construction for loop diagram with ID 31
      CALL MP_FFVS19L2_1(PL(0,33),W(1,2),W(1,10),GC_197,MDL_MT,MDL_WT
     $ ,PL(0,48),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,33),4,COEFS,4,4,WL(1,0,1,48))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,48),2,4,6,1,1,55,H)
C     Coefficient construction for loop diagram with ID 32
      CALL MP_FFS3L2_1(PL(0,34),W(1,20),GC_684,MDL_MT,MDL_WT,PL(0,49)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,49))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,49),3,4,5,1,1,56,H)
C     Coefficient construction for loop diagram with ID 33
      CALL MP_FFS3L2_1(PL(0,34),W(1,22),GC_684,MDL_MT,MDL_WT,PL(0,50)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,50))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,50),3,4,5,1,1,57,H)
C     Coefficient construction for loop diagram with ID 34
      CALL MP_FFS3L2_1(PL(0,34),W(1,24),GC_684,MDL_MT,MDL_WT,PL(0,51)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,51))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,51),3,4,5,1,1,58,H)
C     Coefficient construction for loop diagram with ID 35
      CALL MP_FFS3L2_1(PL(0,34),W(1,26),GC_684,MDL_MT,MDL_WT,PL(0,52)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,52))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,52),3,4,5,1,1,59,H)
C     Coefficient construction for loop diagram with ID 36
      CALL MP_FFS3L2_1(PL(0,34),W(1,27),GC_684,MDL_MT,MDL_WT,PL(0,53)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,34),4,COEFS,4,4,WL(1,0,1,53))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,53),3,4,5,1,1,60,H)
C     Coefficient construction for loop diagram with ID 37
      CALL MP_FFV10L2_1(PL(0,0),W(1,2),GC_11,MDL_MT,MDL_WT,PL(0,54)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,54))
      CALL MP_FFVS19L2_1(PL(0,54),W(1,1),W(1,10),GC_197,MDL_MT,MDL_WT
     $ ,PL(0,55),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,54),4,COEFS,4,4,WL(1,0,1,55))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,55),2,4,7,1,1,61,H)
C     Coefficient construction for loop diagram with ID 38
      CALL MP_FFVV32L2_1(PL(0,0),W(1,1),W(1,2),GC_460,MDL_MT,MDL_WT
     $ ,PL(0,56),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,56))
      CALL MP_FFS3L2_1(PL(0,56),W(1,10),GC_684,MDL_MT,MDL_WT,PL(0,57)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,56),4,COEFS,4,4,WL(1,0,1,57))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,57),2,4,8,1,1,62,H)
C     Coefficient construction for loop diagram with ID 39
      CALL MP_FFVVS63L2_1(PL(0,0),W(1,1),W(1,2),W(1,10),GC_203,MDL_MT
     $ ,MDL_WT,PL(0,58),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,58))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,58),1,4,9,1,1,63,H)
C     At this point, all loop coefficients needed for (QCD=2), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      MP_LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

