C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Gamma(3,2,-1)*ProjP(-1,1)) + Coup(2) * (Gamma(3,2,-1)*ProjM(-1,1))
C     
      SUBROUTINE FFVS15_3_4(F1, F2, V3, COUP1, COUP2, M4, W4,S4)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 DENOM
      COMPLEX*16 V3(*)
      COMPLEX*16 P4(0:3)
      COMPLEX*16 COUP2
      COMPLEX*16 TMP20
      COMPLEX*16 F1(*)
      COMPLEX*16 COUP1
      COMPLEX*16 F2(*)
      REAL*8 M4
      REAL*8 W4
      COMPLEX*16 TMP19
      COMPLEX*16 S4(5)
      S4(1) = +F1(1)+F2(1)+V3(1)
      S4(2) = +F1(2)+F2(2)+V3(2)
      S4(3) = +F1(3)+F2(3)+V3(3)
      S4(4) = +F1(4)+F2(4)+V3(4)
      P4(0) = -S4(1)
      P4(1) = -S4(2)
      P4(2) = -S4(3)
      P4(3) = -S4(4)
      TMP19 = (F1(7)*(F2(5)*(V3(5)-V3(8))-F2(6)*(V3(6)+CI*(V3(7))))
     $ +F1(8)*(F2(5)*(+CI*(V3(7))-V3(6))+F2(6)*(V3(5)+V3(8))))
      TMP20 = (F1(5)*(F2(7)*(V3(5)+V3(8))+F2(8)*(V3(6)+CI*(V3(7))))
     $ +F1(6)*(F2(7)*(V3(6)-CI*(V3(7)))+F2(8)*(V3(5)-V3(8))))
      DENOM = 1D0/(P4(0)**2-P4(1)**2-P4(2)**2-P4(3)**2 - M4 * (M4 -CI*
     $  W4))
      S4(5)= DENOM*(+CI*(COUP1*TMP19+COUP2*TMP20))
      END


