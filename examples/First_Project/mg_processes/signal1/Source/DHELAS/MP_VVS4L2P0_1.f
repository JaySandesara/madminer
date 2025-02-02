C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)
C     
      SUBROUTINE MP_VVS4L2P0_1(P2, S3, COUP, M1, W1, P1, COEFF)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 TMP2
      COMPLEX*32 S3(*)
      REAL*16 M1
      INCLUDE 'coef_specs.inc'
      COMPLEX*32 COEFF(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 P2(0:3)
      REAL*16 W1
      COMPLEX*32 P1(0:3)
      COMPLEX*32 COUP
      P1(0) = +P2(0)+S3(1)
      P1(1) = +P2(1)+S3(2)
      P1(2) = +P2(2)+S3(3)
      P1(3) = +P2(3)+S3(4)
      TMP2 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      COEFF(1,0,1)= COUP*S3(5)*(-CI*(TMP2)+CI*(P1(0)*P2(0)))
      COEFF(2,0,1)= COUP*CI * P1(0)*P2(1)*S3(5)
      COEFF(3,0,1)= COUP*CI * P1(0)*P2(2)*S3(5)
      COEFF(4,0,1)= COUP*CI * P1(0)*P2(3)*S3(5)
      COEFF(1,1,1)= 0Q0
      COEFF(2,1,1)= COUP*CI * P2(1)*S3(5)
      COEFF(3,1,1)= COUP*CI * P2(2)*S3(5)
      COEFF(4,1,1)= COUP*CI * P2(3)*S3(5)
      COEFF(1,2,1)= COUP*S3(5)*(+CI*(P2(1)+P1(1)))
      COEFF(2,2,1)= COUP*CI * P1(0)*S3(5)
      COEFF(3,2,1)= 0Q0
      COEFF(4,2,1)= 0Q0
      COEFF(1,3,1)= COUP*S3(5)*(+CI*(P2(2)+P1(2)))
      COEFF(2,3,1)= 0Q0
      COEFF(3,3,1)= COUP*CI * P1(0)*S3(5)
      COEFF(4,3,1)= 0Q0
      COEFF(1,4,1)= COUP*S3(5)*(+CI*(P2(3)+P1(3)))
      COEFF(2,4,1)= 0Q0
      COEFF(3,4,1)= 0Q0
      COEFF(4,4,1)= COUP*CI * P1(0)*S3(5)
      COEFF(1,5,1)= 0Q0
      COEFF(2,5,1)= 0Q0
      COEFF(3,5,1)= 0Q0
      COEFF(4,5,1)= 0Q0
      COEFF(1,6,1)= 0Q0
      COEFF(2,6,1)= COUP*CI * S3(5)
      COEFF(3,6,1)= 0Q0
      COEFF(4,6,1)= 0Q0
      COEFF(1,7,1)= COUP*CI * S3(5)
      COEFF(2,7,1)= 0Q0
      COEFF(3,7,1)= 0Q0
      COEFF(4,7,1)= 0Q0
      COEFF(1,8,1)= 0Q0
      COEFF(2,8,1)= 0Q0
      COEFF(3,8,1)= COUP*CI * S3(5)
      COEFF(4,8,1)= 0Q0
      COEFF(1,9,1)= 0Q0
      COEFF(2,9,1)= 0Q0
      COEFF(3,9,1)= 0Q0
      COEFF(4,9,1)= 0Q0
      COEFF(1,10,1)= COUP*CI * S3(5)
      COEFF(2,10,1)= 0Q0
      COEFF(3,10,1)= 0Q0
      COEFF(4,10,1)= 0Q0
      COEFF(1,11,1)= 0Q0
      COEFF(2,11,1)= 0Q0
      COEFF(3,11,1)= 0Q0
      COEFF(4,11,1)= COUP*CI * S3(5)
      COEFF(1,12,1)= 0Q0
      COEFF(2,12,1)= 0Q0
      COEFF(3,12,1)= 0Q0
      COEFF(4,12,1)= 0Q0
      COEFF(1,13,1)= 0Q0
      COEFF(2,13,1)= 0Q0
      COEFF(3,13,1)= 0Q0
      COEFF(4,13,1)= 0Q0
      COEFF(1,14,1)= COUP*CI * S3(5)
      COEFF(2,14,1)= 0Q0
      COEFF(3,14,1)= 0Q0
      COEFF(4,14,1)= 0Q0
      COEFF(1,0,2)= COUP*-CI * P1(1)*P2(0)*S3(5)
      COEFF(2,0,2)= COUP*-S3(5)*(+CI*(P1(1)*P2(1)+TMP2))
      COEFF(3,0,2)= COUP*-CI * P1(1)*P2(2)*S3(5)
      COEFF(4,0,2)= COUP*-CI * P1(1)*P2(3)*S3(5)
      COEFF(1,1,2)= COUP*-CI * P1(1)*S3(5)
      COEFF(2,1,2)= COUP*-S3(5)*(+CI*(P2(0)+P1(0)))
      COEFF(3,1,2)= 0Q0
      COEFF(4,1,2)= 0Q0
      COEFF(1,2,2)= COUP*-CI * P2(0)*S3(5)
      COEFF(2,2,2)= 0Q0
      COEFF(3,2,2)= COUP*-CI * P2(2)*S3(5)
      COEFF(4,2,2)= COUP*-CI * P2(3)*S3(5)
      COEFF(1,3,2)= 0Q0
      COEFF(2,3,2)= COUP*S3(5)*(+CI*(P2(2)+P1(2)))
      COEFF(3,3,2)= COUP*-CI * P1(1)*S3(5)
      COEFF(4,3,2)= 0Q0
      COEFF(1,4,2)= 0Q0
      COEFF(2,4,2)= COUP*S3(5)*(+CI*(P2(3)+P1(3)))
      COEFF(3,4,2)= 0Q0
      COEFF(4,4,2)= COUP*-CI * P1(1)*S3(5)
      COEFF(1,5,2)= 0Q0
      COEFF(2,5,2)= COUP*-CI * S3(5)
      COEFF(3,5,2)= 0Q0
      COEFF(4,5,2)= 0Q0
      COEFF(1,6,2)= COUP*-CI * S3(5)
      COEFF(2,6,2)= 0Q0
      COEFF(3,6,2)= 0Q0
      COEFF(4,6,2)= 0Q0
      COEFF(1,7,2)= 0Q0
      COEFF(2,7,2)= 0Q0
      COEFF(3,7,2)= 0Q0
      COEFF(4,7,2)= 0Q0
      COEFF(1,8,2)= 0Q0
      COEFF(2,8,2)= 0Q0
      COEFF(3,8,2)= 0Q0
      COEFF(4,8,2)= 0Q0
      COEFF(1,9,2)= 0Q0
      COEFF(2,9,2)= 0Q0
      COEFF(3,9,2)= COUP*-CI * S3(5)
      COEFF(4,9,2)= 0Q0
      COEFF(1,10,2)= 0Q0
      COEFF(2,10,2)= COUP*CI * S3(5)
      COEFF(3,10,2)= 0Q0
      COEFF(4,10,2)= 0Q0
      COEFF(1,11,2)= 0Q0
      COEFF(2,11,2)= 0Q0
      COEFF(3,11,2)= 0Q0
      COEFF(4,11,2)= 0Q0
      COEFF(1,12,2)= 0Q0
      COEFF(2,12,2)= 0Q0
      COEFF(3,12,2)= 0Q0
      COEFF(4,12,2)= COUP*-CI * S3(5)
      COEFF(1,13,2)= 0Q0
      COEFF(2,13,2)= 0Q0
      COEFF(3,13,2)= 0Q0
      COEFF(4,13,2)= 0Q0
      COEFF(1,14,2)= 0Q0
      COEFF(2,14,2)= COUP*CI * S3(5)
      COEFF(3,14,2)= 0Q0
      COEFF(4,14,2)= 0Q0
      COEFF(1,0,3)= COUP*-CI * P1(2)*P2(0)*S3(5)
      COEFF(2,0,3)= COUP*-CI * P1(2)*P2(1)*S3(5)
      COEFF(3,0,3)= COUP*-S3(5)*(+CI*(P1(2)*P2(2)+TMP2))
      COEFF(4,0,3)= COUP*-CI * P1(2)*P2(3)*S3(5)
      COEFF(1,1,3)= COUP*-CI * P1(2)*S3(5)
      COEFF(2,1,3)= 0Q0
      COEFF(3,1,3)= COUP*-S3(5)*(+CI*(P2(0)+P1(0)))
      COEFF(4,1,3)= 0Q0
      COEFF(1,2,3)= 0Q0
      COEFF(2,2,3)= COUP*-CI * P1(2)*S3(5)
      COEFF(3,2,3)= COUP*S3(5)*(+CI*(P2(1)+P1(1)))
      COEFF(4,2,3)= 0Q0
      COEFF(1,3,3)= COUP*-CI * P2(0)*S3(5)
      COEFF(2,3,3)= COUP*-CI * P2(1)*S3(5)
      COEFF(3,3,3)= 0Q0
      COEFF(4,3,3)= COUP*-CI * P2(3)*S3(5)
      COEFF(1,4,3)= 0Q0
      COEFF(2,4,3)= 0Q0
      COEFF(3,4,3)= COUP*S3(5)*(+CI*(P2(3)+P1(3)))
      COEFF(4,4,3)= COUP*-CI * P1(2)*S3(5)
      COEFF(1,5,3)= 0Q0
      COEFF(2,5,3)= 0Q0
      COEFF(3,5,3)= COUP*-CI * S3(5)
      COEFF(4,5,3)= 0Q0
      COEFF(1,6,3)= 0Q0
      COEFF(2,6,3)= 0Q0
      COEFF(3,6,3)= 0Q0
      COEFF(4,6,3)= 0Q0
      COEFF(1,7,3)= 0Q0
      COEFF(2,7,3)= 0Q0
      COEFF(3,7,3)= COUP*CI * S3(5)
      COEFF(4,7,3)= 0Q0
      COEFF(1,8,3)= COUP*-CI * S3(5)
      COEFF(2,8,3)= 0Q0
      COEFF(3,8,3)= 0Q0
      COEFF(4,8,3)= 0Q0
      COEFF(1,9,3)= 0Q0
      COEFF(2,9,3)= COUP*-CI * S3(5)
      COEFF(3,9,3)= 0Q0
      COEFF(4,9,3)= 0Q0
      COEFF(1,10,3)= 0Q0
      COEFF(2,10,3)= 0Q0
      COEFF(3,10,3)= 0Q0
      COEFF(4,10,3)= 0Q0
      COEFF(1,11,3)= 0Q0
      COEFF(2,11,3)= 0Q0
      COEFF(3,11,3)= 0Q0
      COEFF(4,11,3)= 0Q0
      COEFF(1,12,3)= 0Q0
      COEFF(2,12,3)= 0Q0
      COEFF(3,12,3)= 0Q0
      COEFF(4,12,3)= 0Q0
      COEFF(1,13,3)= 0Q0
      COEFF(2,13,3)= 0Q0
      COEFF(3,13,3)= 0Q0
      COEFF(4,13,3)= COUP*-CI * S3(5)
      COEFF(1,14,3)= 0Q0
      COEFF(2,14,3)= 0Q0
      COEFF(3,14,3)= COUP*CI * S3(5)
      COEFF(4,14,3)= 0Q0
      COEFF(1,0,4)= COUP*-CI * P1(3)*P2(0)*S3(5)
      COEFF(2,0,4)= COUP*-CI * P1(3)*P2(1)*S3(5)
      COEFF(3,0,4)= COUP*-CI * P1(3)*P2(2)*S3(5)
      COEFF(4,0,4)= COUP*-S3(5)*(+CI*(P1(3)*P2(3)+TMP2))
      COEFF(1,1,4)= COUP*-CI * P1(3)*S3(5)
      COEFF(2,1,4)= 0Q0
      COEFF(3,1,4)= 0Q0
      COEFF(4,1,4)= COUP*-S3(5)*(+CI*(P2(0)+P1(0)))
      COEFF(1,2,4)= 0Q0
      COEFF(2,2,4)= COUP*-CI * P1(3)*S3(5)
      COEFF(3,2,4)= 0Q0
      COEFF(4,2,4)= COUP*S3(5)*(+CI*(P2(1)+P1(1)))
      COEFF(1,3,4)= 0Q0
      COEFF(2,3,4)= 0Q0
      COEFF(3,3,4)= COUP*-CI * P1(3)*S3(5)
      COEFF(4,3,4)= COUP*S3(5)*(+CI*(P2(2)+P1(2)))
      COEFF(1,4,4)= COUP*-CI * P2(0)*S3(5)
      COEFF(2,4,4)= COUP*-CI * P2(1)*S3(5)
      COEFF(3,4,4)= COUP*-CI * P2(2)*S3(5)
      COEFF(4,4,4)= 0Q0
      COEFF(1,5,4)= 0Q0
      COEFF(2,5,4)= 0Q0
      COEFF(3,5,4)= 0Q0
      COEFF(4,5,4)= COUP*-CI * S3(5)
      COEFF(1,6,4)= 0Q0
      COEFF(2,6,4)= 0Q0
      COEFF(3,6,4)= 0Q0
      COEFF(4,6,4)= 0Q0
      COEFF(1,7,4)= 0Q0
      COEFF(2,7,4)= 0Q0
      COEFF(3,7,4)= 0Q0
      COEFF(4,7,4)= COUP*CI * S3(5)
      COEFF(1,8,4)= 0Q0
      COEFF(2,8,4)= 0Q0
      COEFF(3,8,4)= 0Q0
      COEFF(4,8,4)= 0Q0
      COEFF(1,9,4)= 0Q0
      COEFF(2,9,4)= 0Q0
      COEFF(3,9,4)= 0Q0
      COEFF(4,9,4)= 0Q0
      COEFF(1,10,4)= 0Q0
      COEFF(2,10,4)= 0Q0
      COEFF(3,10,4)= 0Q0
      COEFF(4,10,4)= COUP*CI * S3(5)
      COEFF(1,11,4)= COUP*-CI * S3(5)
      COEFF(2,11,4)= 0Q0
      COEFF(3,11,4)= 0Q0
      COEFF(4,11,4)= 0Q0
      COEFF(1,12,4)= 0Q0
      COEFF(2,12,4)= COUP*-CI * S3(5)
      COEFF(3,12,4)= 0Q0
      COEFF(4,12,4)= 0Q0
      COEFF(1,13,4)= 0Q0
      COEFF(2,13,4)= 0Q0
      COEFF(3,13,4)= COUP*-CI * S3(5)
      COEFF(4,13,4)= 0Q0
      COEFF(1,14,4)= 0Q0
      COEFF(2,14,4)= 0Q0
      COEFF(3,14,4)= 0Q0
      COEFF(4,14,4)= 0Q0
      END


