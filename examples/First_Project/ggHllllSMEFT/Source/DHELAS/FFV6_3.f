C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-1)*ProjP(-1,1)
C     
      SUBROUTINE FFV6_3(F1, F2, COUP, M3, W3,V3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP5
      COMPLEX*16 V3(8)
      REAL*8 W3
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 P3(0:3)
      COMPLEX*16 F2(*)
      REAL*8 OM3
      COMPLEX*16 DENOM
      COMPLEX*16 COUP
      OM3 = 0D0
      IF (M3.NE.0D0) OM3=1D0/M3**2
      V3(1) = +F1(1)+F2(1)
      V3(2) = +F1(2)+F2(2)
      V3(3) = +F1(3)+F2(3)
      V3(4) = +F1(4)+F2(4)
      P3(0) = -V3(1)
      P3(1) = -V3(2)
      P3(2) = -V3(3)
      P3(3) = -V3(4)
      TMP5 = (F1(7)*(F2(5)*(P3(0)-P3(3))-F2(6)*(P3(1)+CI*(P3(2))))
     $ +F1(8)*(F2(5)*(+CI*(P3(2))-P3(1))+F2(6)*(P3(0)+P3(3))))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI
     $ * W3))
      V3(5)= DENOM*(-CI)*(F1(7)*F2(5)+F1(8)*F2(6)-P3(0)*OM3*TMP5)
      V3(6)= DENOM*(-CI)*(F1(7)*F2(6)+F1(8)*F2(5)-P3(1)*OM3*TMP5)
      V3(7)= DENOM*(-CI)*(-CI*(F1(8)*F2(5))+CI*(F1(7)*F2(6))-P3(2)*OM3
     $ *TMP5)
      V3(8)= DENOM*(-CI)*(F1(7)*F2(5)-F1(8)*F2(6)-P3(3)*OM3*TMP5)
      END


