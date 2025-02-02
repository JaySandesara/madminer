      SUBROUTINE SMATRIX_2(P,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.6.7, 2019-10-16
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: b~ g > h > e+ e- mu+ mu- b~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: c~ g > h > e+ e- mu+ mu- c~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: d~ g > h > e+ e- mu+ mu- d~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: s~ g > h > e+ e- mu+ mu- s~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: u~ g > h > e+ e- mu+ mu- u~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INCLUDE 'nexternal.inc'
      INTEGER     NCOMB
      PARAMETER ( NCOMB=128)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
      DOUBLE PRECISION       WGT_ME_BORN,WGT_ME_REAL
      COMMON /C_WGT_ME_TREE/ WGT_ME_BORN,WGT_ME_REAL
C     
C     LOCAL VARIABLES 
C     
      INTEGER IHEL,IDEN,I,T_IDENT(NCOMB)
      REAL*8 MATRIX_2
      REAL*8 T,T_SAVE(NCOMB)
      SAVE T_SAVE,T_IDENT
      INTEGER NHEL(NEXTERNAL,NCOMB)
      DATA (NHEL(I,   1),I=1,7) /-1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,   2),I=1,7) /-1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,   3),I=1,7) /-1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,   4),I=1,7) /-1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,   5),I=1,7) /-1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,   6),I=1,7) /-1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,   7),I=1,7) /-1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,   8),I=1,7) /-1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,   9),I=1,7) /-1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  10),I=1,7) /-1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  11),I=1,7) /-1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  12),I=1,7) /-1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  13),I=1,7) /-1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  14),I=1,7) /-1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  15),I=1,7) /-1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  16),I=1,7) /-1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  17),I=1,7) /-1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  18),I=1,7) /-1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  19),I=1,7) /-1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  20),I=1,7) /-1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  21),I=1,7) /-1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  22),I=1,7) /-1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  23),I=1,7) /-1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  24),I=1,7) /-1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  25),I=1,7) /-1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  26),I=1,7) /-1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  27),I=1,7) /-1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  28),I=1,7) /-1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  29),I=1,7) /-1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  30),I=1,7) /-1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  31),I=1,7) /-1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  32),I=1,7) /-1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  33),I=1,7) /-1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  34),I=1,7) /-1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  35),I=1,7) /-1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  36),I=1,7) /-1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  37),I=1,7) /-1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  38),I=1,7) /-1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  39),I=1,7) /-1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  40),I=1,7) /-1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  41),I=1,7) /-1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  42),I=1,7) /-1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  43),I=1,7) /-1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  44),I=1,7) /-1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  45),I=1,7) /-1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  46),I=1,7) /-1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  47),I=1,7) /-1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  48),I=1,7) /-1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  49),I=1,7) /-1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  50),I=1,7) /-1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  51),I=1,7) /-1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  52),I=1,7) /-1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  53),I=1,7) /-1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  54),I=1,7) /-1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  55),I=1,7) /-1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  56),I=1,7) /-1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  57),I=1,7) /-1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  58),I=1,7) /-1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  59),I=1,7) /-1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  60),I=1,7) /-1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  61),I=1,7) /-1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  62),I=1,7) /-1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  63),I=1,7) /-1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  64),I=1,7) /-1, 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  65),I=1,7) / 1,-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  66),I=1,7) / 1,-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  67),I=1,7) / 1,-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  68),I=1,7) / 1,-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  69),I=1,7) / 1,-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  70),I=1,7) / 1,-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  71),I=1,7) / 1,-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  72),I=1,7) / 1,-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  73),I=1,7) / 1,-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  74),I=1,7) / 1,-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  75),I=1,7) / 1,-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  76),I=1,7) / 1,-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  77),I=1,7) / 1,-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  78),I=1,7) / 1,-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  79),I=1,7) / 1,-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  80),I=1,7) / 1,-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  81),I=1,7) / 1,-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  82),I=1,7) / 1,-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  83),I=1,7) / 1,-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  84),I=1,7) / 1,-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  85),I=1,7) / 1,-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  86),I=1,7) / 1,-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  87),I=1,7) / 1,-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  88),I=1,7) / 1,-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  89),I=1,7) / 1,-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  90),I=1,7) / 1,-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  91),I=1,7) / 1,-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  92),I=1,7) / 1,-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  93),I=1,7) / 1,-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  94),I=1,7) / 1,-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  95),I=1,7) / 1,-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  96),I=1,7) / 1,-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  97),I=1,7) / 1, 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  98),I=1,7) / 1, 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  99),I=1,7) / 1, 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I, 100),I=1,7) / 1, 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I, 101),I=1,7) / 1, 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I, 102),I=1,7) / 1, 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I, 103),I=1,7) / 1, 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I, 104),I=1,7) / 1, 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I, 105),I=1,7) / 1, 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I, 106),I=1,7) / 1, 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I, 107),I=1,7) / 1, 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I, 108),I=1,7) / 1, 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I, 109),I=1,7) / 1, 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I, 110),I=1,7) / 1, 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I, 111),I=1,7) / 1, 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I, 112),I=1,7) / 1, 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I, 113),I=1,7) / 1, 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I, 114),I=1,7) / 1, 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I, 115),I=1,7) / 1, 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I, 116),I=1,7) / 1, 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I, 117),I=1,7) / 1, 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I, 118),I=1,7) / 1, 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I, 119),I=1,7) / 1, 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I, 120),I=1,7) / 1, 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I, 121),I=1,7) / 1, 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I, 122),I=1,7) / 1, 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I, 123),I=1,7) / 1, 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I, 124),I=1,7) / 1, 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I, 125),I=1,7) / 1, 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I, 126),I=1,7) / 1, 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I, 127),I=1,7) / 1, 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I, 128),I=1,7) / 1, 1,-1, 1,-1, 1,-1/
      LOGICAL GOODHEL(NCOMB)
      DATA GOODHEL/NCOMB*.FALSE./
      INTEGER NTRY
      DATA NTRY/0/
      DATA IDEN/96/
C     ----------
C     BEGIN CODE
C     ----------
      NTRY=NTRY+1
      ANS = 0D0
      DO IHEL=1,NCOMB
        IF (GOODHEL(IHEL) .OR. NTRY .LT. 2) THEN
          IF (NTRY.LT.2) THEN
C           for the first ps-point, check for helicities that give
C           identical matrix elements
            T=MATRIX_2(P ,NHEL(1,IHEL))
            T_SAVE(IHEL)=T
            T_IDENT(IHEL)=-1
            DO I=1,IHEL-1
              IF (T.EQ.0D0) EXIT
              IF (T_SAVE(I).EQ.0D0) CYCLE
              IF (ABS(T/T_SAVE(I)-1D0) .LT. 1D-12) THEN
C               WRITE (*,*) 'FOUND IDENTICAL',T,IHEL,T_SAVE(I),I
                T_IDENT(IHEL) = I
              ENDIF
            ENDDO
          ELSE
            IF (T_IDENT(IHEL).GT.0) THEN
C             if two helicity states are identical, dont recompute
              T=T_SAVE(T_IDENT(IHEL))
              T_SAVE(IHEL)=T
            ELSE
              T=MATRIX_2(P ,NHEL(1,IHEL))
              T_SAVE(IHEL)=T
            ENDIF
          ENDIF
C         add to the sum of helicities
          ANS=ANS+T
          IF (T .NE. 0D0 .AND. .NOT. GOODHEL(IHEL)) THEN
            GOODHEL(IHEL)=.TRUE.
          ENDIF
        ENDIF
      ENDDO
      ANS=ANS/DBLE(IDEN)
      WGT_ME_REAL=ANS
      END


      REAL*8 FUNCTION MATRIX_2(P,NHEL)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.6.7, 2019-10-16
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: b~ g > h > e+ e- mu+ mu- b~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: c~ g > h > e+ e- mu+ mu- c~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: d~ g > h > e+ e- mu+ mu- d~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: s~ g > h > e+ e- mu+ mu- s~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     Process: u~ g > h > e+ e- mu+ mu- u~ NP<=2 QCD<=5 QED<=4
C      WEIGHTED<=28 [ all = QCD ]
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=1)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=8, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      INTEGER IC(NEXTERNAL)
      DATA IC /NEXTERNAL*1/
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 ZTEMP, AMP(NGRAPHS), JAMP(NCOLOR), W(8,NWAVEFUNCS)
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    4/
C     1 T(2,1,7)
C     ----------
C     BEGIN CODE
C     ----------
      CALL OXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL IXXXXX(P(0,3),ZERO,NHEL(3),-1*IC(3),W(1,3))
      CALL OXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX(P(0,6),ZERO,NHEL(6),+1*IC(6),W(1,6))
      CALL IXXXXX(P(0,7),ZERO,NHEL(7),-1*IC(7),W(1,7))
      CALL FFV1P0_3(W(1,7),W(1,1),GC_11,ZERO,ZERO,W(1,8))
      CALL FFV12_2_3(W(1,3),W(1,4),GC_317,GC_243,MDL_MZ,MDL_WZ,W(1,7))
      CALL FFV12_2_3(W(1,5),W(1,6),GC_317,GC_243,MDL_MZ,MDL_WZ,W(1,4))
      CALL VVS4_3(W(1,8),W(1,2),GC_426,MDL_MH,MDL_WH,W(1,6))
C     Amplitude(s) for diagram number 1
      CALL VVS3_0(W(1,7),W(1,4),W(1,6),GC_588,AMP(1))
      JAMP(1)=-AMP(1)
      MATRIX_2 = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
        ENDDO
        MATRIX_2 = MATRIX_2+ZTEMP*DCONJG(JAMP(I))/DENOM(I)
      ENDDO
      END

