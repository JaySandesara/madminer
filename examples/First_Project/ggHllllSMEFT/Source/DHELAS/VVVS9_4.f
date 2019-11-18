C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) +
C      P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)
C     
      SUBROUTINE VVVS9_4(V1, V2, V3, COUP, M4, W4,S4)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP12
      COMPLEX*16 V3(*)
      COMPLEX*16 TMP11
      COMPLEX*16 TMP0
      COMPLEX*16 P2(0:3)
      COMPLEX*16 TMP10
      REAL*8 W4
      COMPLEX*16 P3(0:3)
      COMPLEX*16 S4(5)
      COMPLEX*16 DENOM
      COMPLEX*16 P4(0:3)
      COMPLEX*16 P1(0:3)
      REAL*8 M4
      COMPLEX*16 COUP
      COMPLEX*16 TMP1
      COMPLEX*16 TMP18
      COMPLEX*16 TMP9
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP13
      COMPLEX*16 TMP8
      P1(0) = V1(1)
      P1(1) = V1(2)
      P1(2) = V1(3)
      P1(3) = V1(4)
      P2(0) = V2(1)
      P2(1) = V2(2)
      P2(2) = V2(3)
      P2(3) = V2(4)
      P3(0) = V3(1)
      P3(1) = V3(2)
      P3(2) = V3(3)
      P3(3) = V3(4)
      S4(1) = +V1(1)+V2(1)+V3(1)
      S4(2) = +V1(2)+V2(2)+V3(2)
      S4(3) = +V1(3)+V2(3)+V3(3)
      S4(4) = +V1(4)+V2(4)+V3(4)
      P4(0) = -S4(1)
      P4(1) = -S4(2)
      P4(2) = -S4(3)
      P4(3) = -S4(4)
      TMP11 = (P2(0)*V1(5)-P2(1)*V1(6)-P2(2)*V1(7)-P2(3)*V1(8))
      TMP0 = (V3(5)*P1(0)-V3(6)*P1(1)-V3(7)*P1(2)-V3(8)*P1(3))
      TMP9 = (V2(5)*P3(0)-V2(6)*P3(1)-V2(7)*P3(2)-V2(8)*P3(3))
      TMP8 = (V2(5)*P1(0)-V2(6)*P1(1)-V2(7)*P1(2)-V2(8)*P1(3))
      TMP18 = (V3(5)*V1(5)-V3(6)*V1(6)-V3(7)*V1(7)-V3(8)*V1(8))
      TMP1 = (V3(5)*P2(0)-V3(6)*P2(1)-V3(7)*P2(2)-V3(8)*P2(3))
      TMP10 = (V2(5)*V3(5)-V2(6)*V3(6)-V2(7)*V3(7)-V2(8)*V3(8))
      TMP13 = (V2(5)*V1(5)-V2(6)*V1(6)-V2(7)*V1(7)-V2(8)*V1(8))
      TMP12 = (P3(0)*V1(5)-P3(1)*V1(6)-P3(2)*V1(7)-P3(3)*V1(8))
      DENOM = COUP/(P4(0)**2-P4(1)**2-P4(2)**2-P4(3)**2 - M4 * (M4 -CI
     $ * W4))
      S4(5)= DENOM*(TMP10*(-CI*(TMP12)+CI*(TMP11))+(TMP13*(-CI*(TMP1)
     $ +CI*(TMP0))+TMP18*(-CI*(TMP8)+CI*(TMP9))))
      END


