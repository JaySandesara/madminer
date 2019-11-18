C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Metric(1,2)) + Coup(2) * (P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2))
C     
      SUBROUTINE VVS3_4_0(V1, V2, S3, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 S3(*)
      COMPLEX*16 TMP11
      COMPLEX*16 P2(0:3)
      COMPLEX*16 COUP1
      COMPLEX*16 COUP2
      COMPLEX*16 P1(0:3)
      COMPLEX*16 VERTEX
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
      TMP11 = (P2(0)*V1(5)-P2(1)*V1(6)-P2(2)*V1(7)-P2(3)*V1(8))
      TMP8 = (V2(5)*P1(0)-V2(6)*P1(1)-V2(7)*P1(2)-V2(8)*P1(3))
      TMP13 = (V2(5)*V1(5)-V2(6)*V1(6)-V2(7)*V1(7)-V2(8)*V1(8))
      TMP2 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      VERTEX = -S3(5)*(TMP13*(-CI*(TMP2*COUP2)+CI*(COUP1))+CI*(TMP8
     $ *TMP11*COUP2))
      END

