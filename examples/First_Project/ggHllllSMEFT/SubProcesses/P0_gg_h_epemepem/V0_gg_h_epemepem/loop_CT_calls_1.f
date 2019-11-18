      SUBROUTINE LOOP_CT_CALLS_1(P,NHEL,H,IC)
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
      IF (FILTER_SO.AND.CTCALL_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     CutTools call for loop numbers 1,2
      CALL LOOP_3(1,2,12,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),4,I_SO
     $ ,1)
C     CutTools call for loop numbers 3,4
      CALL LOOP_2_3(1,2,1,2,12,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,2)
C     CutTools call for loop numbers 5,6
      CALL LOOP_2_3(1,2,2,1,12,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,3)
C     CutTools call for loop numbers 7,8,9,10,11,12
      CALL LOOP_2_3(2,1,1,2,12,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,4)
C     CutTools call for loop numbers 13,14,15,16,17,18,19,20,21,22,23,2
C     4,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,
C     46,47,48,49,50,51,52,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69
C     ,70,71,72,73,74
      CALL LOOP_3(1,2,14,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,5)
C     CutTools call for loop numbers 53,68
      CALL LOOP_2_3(1,2,1,2,12,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,6)
C     CutTools call for loop numbers 75,76
      CALL LOOP_2_3(1,2,2,1,12,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,7)
C     CutTools call for loop numbers 77,78
      CALL LOOP_2_3(2,1,1,2,12,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,8)
C     CutTools call for loop numbers 79,80
      CALL LOOP_1_3(3,1,2,12,DCMPLX(MDL_MT),1,I_SO,9)
C     At this point, all reductions needed for (QCD=2), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 5000

      GOTO 1001
 5000 CONTINUE
      CTCALL_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END
