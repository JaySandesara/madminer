C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Gamma(3,2,-1)*ProjM(-1,1) + Gamma(3,2,-1)*ProjP(-1,1)) + Coup(2) * (Gamma(3,2,-1)*ProjP(-1,1)) + Coup(3) * (Gamma(3,2,-1)*ProjM(-1,1) + (2*Gamma(3,2,-1)*ProjP(-1,1))/3.)
C     
      SUBROUTINE MP_FFV10_6_9P0_3(F1, F2, COUP1, COUP2, COUP3, M3, W3
     $ ,V3)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 F2(*)
      COMPLEX*32 V3(8)
      REAL*16 W3
      COMPLEX*32 COUP2
      REAL*16 M3
      COMPLEX*32 P3(0:3)
      COMPLEX*32 F1(*)
      COMPLEX*32 DENOM
      COMPLEX*32 COUP1
      COMPLEX*32 COUP3
      V3(1) = +F1(1)+F2(1)
      V3(2) = +F1(2)+F2(2)
      V3(3) = +F1(3)+F2(3)
      V3(4) = +F1(4)+F2(4)
      P3(0) = -V3(1)
      P3(1) = -V3(2)
      P3(2) = -V3(3)
      P3(3) = -V3(4)
      DENOM = 1Q0/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 -CI*
     $  W3))
      V3(5)= DENOM*(-2Q0/3Q0 * CI)*(COUP1*3Q0/2Q0*(F1(5)*F2(7)+F1(6)
     $ *F2(8)+F1(7)*F2(5)+F1(8)*F2(6))+(+3Q0/2Q0*(COUP2*(F1(7)*F2(5)
     $ +F1(8)*F2(6)))+COUP3*3Q0/2Q0*(F1(5)*F2(7)+F1(6)*F2(8)+2Q0/3Q0
     $ *(F1(7)*F2(5)+F1(8)*F2(6)))))
      V3(6)= DENOM*(-2Q0/3Q0 * CI)*(COUP1*3Q0/2Q0*(F1(7)*F2(6)+F1(8)
     $ *F2(5)-F1(5)*F2(8)-F1(6)*F2(7))+(+3Q0/2Q0*(COUP2*(F1(7)*F2(6)
     $ +F1(8)*F2(5)))+COUP3*3Q0/2Q0*(+2Q0/3Q0*(F1(7)*F2(6)+F1(8)*F2(5))
     $ -F1(5)*F2(8)-F1(6)*F2(7))))
      V3(7)= DENOM*(-CI)*(COUP1*(-CI*(F1(5)*F2(8)+F1(8)*F2(5))+CI
     $ *(F1(6)*F2(7)+F1(7)*F2(6)))+(COUP3*(-2Q0/3Q0 * CI*(F1(8)*F2(5))
     $ -CI*(F1(5)*F2(8))+CI*(F1(6)*F2(7))+2Q0/3Q0 * CI*(F1(7)*F2(6)))
     $ +COUP2*(-CI*(F1(8)*F2(5))+CI*(F1(7)*F2(6)))))
      V3(8)= DENOM*(-2Q0/3Q0 * CI)*(COUP1*3Q0/2Q0*(F1(6)*F2(8)+F1(7)
     $ *F2(5)-F1(5)*F2(7)-F1(8)*F2(6))+(+3Q0/2Q0*(COUP2*(F1(7)*F2(5)
     $ -F1(8)*F2(6)))+COUP3*3Q0/2Q0*(F1(6)*F2(8)-2Q0/3Q0*(F1(8)*F2(6))
     $ +2Q0/3Q0*(F1(7)*F2(5))-F1(5)*F2(7))))
      END

