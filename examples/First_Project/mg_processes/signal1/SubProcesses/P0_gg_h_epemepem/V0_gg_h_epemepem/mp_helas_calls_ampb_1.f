      SUBROUTINE MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
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
      PARAMETER (NBORNAMPS=2)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=80, NLOOPGROUPS=9, NCTAMPS=46)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=126)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=47,NLOOPWAVEFUNCS=98)
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
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_IXXXXX(P(0,3),ZERO,NHEL(3),-1*IC(3),W(1,3))
      CALL MP_OXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL MP_IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))
      CALL MP_OXXXXX(P(0,6),ZERO,NHEL(6),+1*IC(6),W(1,6))
      CALL MP_VVS4_3(W(1,1),W(1,2),GC_426,MDL_MH,MDL_WH,W(1,7))
      CALL MP_FFV12_2_3(W(1,3),W(1,4),GC_317,GC_243,MDL_MZ,MDL_WZ,W(1
     $ ,8))
      CALL MP_FFV12_2_3(W(1,5),W(1,6),GC_317,GC_243,MDL_MZ,MDL_WZ,W(1
     $ ,9))
C     Amplitude(s) for born diagram with ID 1
      CALL MP_VVS3_0(W(1,8),W(1,9),W(1,7),GC_588,AMP(1))
      CALL MP_FFV12_2_3(W(1,3),W(1,6),GC_317,GC_243,MDL_MZ,MDL_WZ,W(1
     $ ,10))
      CALL MP_FFV12_2_3(W(1,5),W(1,4),GC_317,GC_243,MDL_MZ,MDL_WZ,W(1
     $ ,11))
C     Amplitude(s) for born diagram with ID 2
      CALL MP_VVS3_0(W(1,10),W(1,11),W(1,7),GC_588,AMP(2))
      CALL MP_VVS3_3(W(1,9),W(1,8),GC_588,MDL_MH,MDL_WH,W(1,12))
C     Counter-term amplitude(s) for loop diagram number 3
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,12),UVGC_1217_229_1EPS,AMPL(2,1)
     $ )
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,12),UVGC_1217_229_1EPS,AMPL(2,2)
     $ )
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,12),UVGC_1217_229_1EPS,AMPL(2,3)
     $ )
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,12),UVGC_1217_229_1EPS,AMPL(2,4)
     $ )
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,12),UVGC_1217_229_1EPS,AMPL(2,5)
     $ )
      CALL MP_VVS7_0(W(1,1),W(1,2),W(1,12),UVGC_991_428_1EPS,AMPL(2,6))
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,12),UVGC_1217_230_1EPS,AMPL(2,7)
     $ )
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,12),UVGC_1217_231,AMPL(1,8))
      CALL MP_VVS5_0(W(1,1),W(1,2),W(1,12),UVGC_992_429_1EPS,AMPL(2,9))
      CALL MP_VVS8_0(W(1,1),W(1,2),W(1,12),R2GC_694_147,AMPL(1,10))
      CALL MP_VVS3_3(W(1,11),W(1,10),GC_588,MDL_MH,MDL_WH,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 4
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,13),UVGC_1217_229_1EPS,AMPL(2
     $ ,11))
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,13),UVGC_1217_229_1EPS,AMPL(2
     $ ,12))
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,13),UVGC_1217_229_1EPS,AMPL(2
     $ ,13))
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,13),UVGC_1217_229_1EPS,AMPL(2
     $ ,14))
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,13),UVGC_1217_229_1EPS,AMPL(2
     $ ,15))
      CALL MP_VVS7_0(W(1,1),W(1,2),W(1,13),UVGC_991_428_1EPS,AMPL(2,16)
     $ )
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,13),UVGC_1217_230_1EPS,AMPL(2
     $ ,17))
      CALL MP_VVS4_0(W(1,1),W(1,2),W(1,13),UVGC_1217_231,AMPL(1,18))
      CALL MP_VVS5_0(W(1,1),W(1,2),W(1,13),UVGC_992_429_1EPS,AMPL(2,19)
     $ )
      CALL MP_VVS8_0(W(1,1),W(1,2),W(1,13),R2GC_694_147,AMPL(1,20))
      CALL MP_FFVS15_3_4(W(1,5),W(1,6),W(1,8),GC_571,GC_573,MDL_MH
     $ ,MDL_WH,W(1,14))
C     Counter-term amplitude(s) for loop diagram number 11
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,14),R2GC_709_162,AMPL(1,21))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL MP_VVS3_6_0(W(1,1),W(1,2),W(1,12),R2GC_900_292,R2GC_701_154
     $ ,AMPL(1,22))
      CALL MP_FFV10P0_3(W(1,5),W(1,6),GC_3,ZERO,ZERO,W(1,15))
      CALL MP_FFV10P0_3(W(1,3),W(1,4),GC_3,ZERO,ZERO,W(1,16))
      CALL MP_VVS4_3(W(1,15),W(1,16),GC_590,MDL_MH,MDL_WH,W(1,17))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,17),R2GC_709_162,AMPL(1,23))
      CALL MP_VVS3_4_3(W(1,16),W(1,9),GC_681,GC_591,MDL_MH,MDL_WH,W(1
     $ ,18))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,18),R2GC_709_162,AMPL(1,24))
      CALL MP_VVS3_4_3(W(1,15),W(1,8),GC_681,GC_591,MDL_MH,MDL_WH,W(1
     $ ,19))
C     Counter-term amplitude(s) for loop diagram number 16
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,19),R2GC_709_162,AMPL(1,25))
      CALL MP_VVS4_3(W(1,9),W(1,8),GC_589,MDL_MH,MDL_WH,W(1,20))
C     Counter-term amplitude(s) for loop diagram number 18
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,20),R2GC_709_162,AMPL(1,26))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,12),R2GC_709_162,AMPL(1,27))
      CALL MP_VVS3_3(W(1,9),W(1,8),GC_683,MDL_MH,MDL_WH,W(1,21))
C     Counter-term amplitude(s) for loop diagram number 20
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,21),R2GC_709_162,AMPL(1,28))
      CALL MP_FFV10_6_9_3(W(1,5),W(1,6),GC_617,GC_644,GC_535,MDL_MZ
     $ ,MDL_WZ,W(1,22))
      CALL MP_VVS3_3(W(1,22),W(1,8),GC_588,MDL_MH,MDL_WH,W(1,23))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,23),R2GC_709_162,AMPL(1,29))
      CALL MP_FFV12_2_3(W(1,5),W(1,6),GC_638,GC_645,MDL_MZ,MDL_WZ,W(1
     $ ,24))
      CALL MP_VVS3_3(W(1,24),W(1,8),GC_588,MDL_MH,MDL_WH,W(1,25))
C     Counter-term amplitude(s) for loop diagram number 22
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,25),R2GC_709_162,AMPL(1,30))
      CALL MP_FFV10_6_9_3(W(1,3),W(1,4),GC_617,GC_644,GC_535,MDL_MZ
     $ ,MDL_WZ,W(1,26))
      CALL MP_VVS3_3(W(1,9),W(1,26),GC_588,MDL_MH,MDL_WH,W(1,27))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,27),R2GC_709_162,AMPL(1,31))
      CALL MP_FFV12_2_3(W(1,3),W(1,4),GC_638,GC_645,MDL_MZ,MDL_WZ,W(1
     $ ,28))
      CALL MP_VVS3_3(W(1,9),W(1,28),GC_588,MDL_MH,MDL_WH,W(1,29))
C     Counter-term amplitude(s) for loop diagram number 24
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,29),R2GC_709_162,AMPL(1,32))
      CALL MP_FFVS15_3_4(W(1,5),W(1,4),W(1,10),GC_571,GC_573,MDL_MH
     $ ,MDL_WH,W(1,30))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,30),R2GC_709_162,AMPL(1,33))
C     Counter-term amplitude(s) for loop diagram number 26
      CALL MP_VVS3_6_0(W(1,1),W(1,2),W(1,13),R2GC_900_292,R2GC_701_154
     $ ,AMPL(1,34))
      CALL MP_FFV10P0_3(W(1,5),W(1,4),GC_3,ZERO,ZERO,W(1,31))
      CALL MP_FFV10P0_3(W(1,3),W(1,6),GC_3,ZERO,ZERO,W(1,32))
      CALL MP_VVS4_3(W(1,31),W(1,32),GC_590,MDL_MH,MDL_WH,W(1,33))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,33),R2GC_709_162,AMPL(1,35))
      CALL MP_VVS3_4_3(W(1,32),W(1,11),GC_681,GC_591,MDL_MH,MDL_WH,W(1
     $ ,34))
C     Counter-term amplitude(s) for loop diagram number 29
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,34),R2GC_709_162,AMPL(1,36))
      CALL MP_VVS3_4_3(W(1,31),W(1,10),GC_681,GC_591,MDL_MH,MDL_WH,W(1
     $ ,35))
C     Counter-term amplitude(s) for loop diagram number 30
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,35),R2GC_709_162,AMPL(1,37))
      CALL MP_VVS4_3(W(1,11),W(1,10),GC_589,MDL_MH,MDL_WH,W(1,36))
C     Counter-term amplitude(s) for loop diagram number 32
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,36),R2GC_709_162,AMPL(1,38))
C     Counter-term amplitude(s) for loop diagram number 33
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,13),R2GC_709_162,AMPL(1,39))
      CALL MP_VVS3_3(W(1,11),W(1,10),GC_683,MDL_MH,MDL_WH,W(1,37))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,37),R2GC_709_162,AMPL(1,40))
      CALL MP_FFV10_6_9_3(W(1,5),W(1,4),GC_617,GC_644,GC_535,MDL_MZ
     $ ,MDL_WZ,W(1,38))
      CALL MP_VVS3_3(W(1,38),W(1,10),GC_588,MDL_MH,MDL_WH,W(1,39))
C     Counter-term amplitude(s) for loop diagram number 35
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,39),R2GC_709_162,AMPL(1,41))
      CALL MP_FFV12_2_3(W(1,5),W(1,4),GC_638,GC_645,MDL_MZ,MDL_WZ,W(1
     $ ,40))
      CALL MP_VVS3_3(W(1,40),W(1,10),GC_588,MDL_MH,MDL_WH,W(1,41))
C     Counter-term amplitude(s) for loop diagram number 36
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,41),R2GC_709_162,AMPL(1,42))
      CALL MP_FFV10_6_9_3(W(1,3),W(1,6),GC_617,GC_644,GC_535,MDL_MZ
     $ ,MDL_WZ,W(1,42))
      CALL MP_VVS3_3(W(1,11),W(1,42),GC_588,MDL_MH,MDL_WH,W(1,43))
C     Counter-term amplitude(s) for loop diagram number 37
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,43),R2GC_709_162,AMPL(1,43))
      CALL MP_FFV12_2_3(W(1,3),W(1,6),GC_638,GC_645,MDL_MZ,MDL_WZ,W(1
     $ ,44))
      CALL MP_VVS3_3(W(1,11),W(1,44),GC_588,MDL_MH,MDL_WH,W(1,45))
C     Counter-term amplitude(s) for loop diagram number 38
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,45),R2GC_709_162,AMPL(1,44))
      CALL MP_FFVS15_3_4(W(1,3),W(1,6),W(1,11),GC_571,GC_573,MDL_MH
     $ ,MDL_WH,W(1,46))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,46),R2GC_709_162,AMPL(1,45))
      CALL MP_FFVS15_3_4(W(1,3),W(1,4),W(1,9),GC_571,GC_573,MDL_MH
     $ ,MDL_WH,W(1,47))
C     Counter-term amplitude(s) for loop diagram number 40
      CALL MP_VVS3_0(W(1,1),W(1,2),W(1,47),R2GC_709_162,AMPL(1,46))
C     At this point, all CT amps needed for (QCD=2), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

