C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Metric(1,2)) + Coup(2) * (P(1,2)*P(2,1) - (P(1,3)*P(2,1))/2. - (P(1,2)*P(2,3))/2. - (P(-1,1)*P(-1,2)*Metric(1,2))/2. + (P(-1,2)**2*Metric(1,2))/2. + (P(-1,1)*P(-1,3)*Metric(1,2))/2. + P(-1,2)*P(-1,3)*Metric(1,2))
C     
      SUBROUTINE VVS3_6_0(V1, V2, S3, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 S3(*)
      COMPLEX*16 TMP11
      COMPLEX*16 COUP2
      COMPLEX*16 TMP12
      COMPLEX*16 COUP1
      COMPLEX*16 P2(0:3)
      COMPLEX*16 P3(0:3)
      COMPLEX*16 TMP16
      COMPLEX*16 TMP15
      COMPLEX*16 TMP14
      COMPLEX*16 P1(0:3)
      COMPLEX*16 VERTEX
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
      P3(0) = S3(1)
      P3(1) = S3(2)
      P3(2) = S3(3)
      P3(3) = S3(4)
      TMP9 = (V2(5)*P3(0)-V2(6)*P3(1)-V2(7)*P3(2)-V2(8)*P3(3))
      TMP8 = (V2(5)*P1(0)-V2(6)*P1(1)-V2(7)*P1(2)-V2(8)*P1(3))
      TMP2 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP15 = (P1(0)*P3(0)-P1(1)*P3(1)-P1(2)*P3(2)-P1(3)*P3(3))
      TMP14 = (P2(0)*P2(0)-P2(1)*P2(1)-P2(2)*P2(2)-P2(3)*P2(3))
      TMP16 = (P2(0)*P3(0)-P2(1)*P3(1)-P2(2)*P3(2)-P2(3)*P3(3))
      TMP11 = (P2(0)*V1(5)-P2(1)*V1(6)-P2(2)*V1(7)-P2(3)*V1(8))
      TMP13 = (V2(5)*V1(5)-V2(6)*V1(6)-V2(7)*V1(7)-V2(8)*V1(8))
      TMP12 = (P3(0)*V1(5)-P3(1)*V1(6)-P3(2)*V1(7)-P3(3)*V1(8))
      VERTEX = -1D0/2D0 * S3(5)*(COUP2*(TMP13*(-CI*(TMP2)+CI*(TMP14
     $ +TMP15)+2D0 * CI*(TMP16))+(TMP11*(-CI*(TMP9)+2D0 * CI*(TMP8))
     $ -CI*(TMP8*TMP12)))+2D0 * CI*(TMP13*COUP1))
      END


