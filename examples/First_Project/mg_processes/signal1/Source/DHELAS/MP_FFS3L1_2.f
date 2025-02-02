C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Identity(2,1)
C     
      SUBROUTINE MP_FFS3L1_2(P1, S3, COUP, M2, W2, P2, COEFF)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 S3(*)
      INCLUDE 'coef_specs.inc'
      COMPLEX*32 COEFF(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      REAL*16 W2
      COMPLEX*32 P2(0:3)
      REAL*16 M2
      COMPLEX*32 P1(0:3)
      COMPLEX*32 COUP
      P2(0) = +P1(0)+S3(1)
      P2(1) = +P1(1)+S3(2)
      P2(2) = +P1(2)+S3(3)
      P2(3) = +P1(3)+S3(4)
      COEFF(1,0,1)= COUP*CI * S3(5)*M2
      COEFF(2,0,1)= 0Q0
      COEFF(3,0,1)= COUP*-S3(5)*(+CI*(P2(0)+P2(3)))
      COEFF(4,0,1)= COUP*-S3(5)*(+CI*(P2(1))-P2(2))
      COEFF(1,1,1)= 0Q0
      COEFF(2,1,1)= 0Q0
      COEFF(3,1,1)= COUP*-CI * S3(5)
      COEFF(4,1,1)= 0Q0
      COEFF(1,2,1)= 0Q0
      COEFF(2,2,1)= 0Q0
      COEFF(3,2,1)= 0Q0
      COEFF(4,2,1)= COUP*-CI * S3(5)
      COEFF(1,3,1)= 0Q0
      COEFF(2,3,1)= 0Q0
      COEFF(3,3,1)= 0Q0
      COEFF(4,3,1)= COUP*S3(5)
      COEFF(1,4,1)= 0Q0
      COEFF(2,4,1)= 0Q0
      COEFF(3,4,1)= COUP*-CI * S3(5)
      COEFF(4,4,1)= 0Q0
      COEFF(1,0,2)= 0Q0
      COEFF(2,0,2)= COUP*CI * S3(5)*M2
      COEFF(3,0,2)= COUP*-S3(5)*(P2(2)+CI*(P2(1)))
      COEFF(4,0,2)= COUP*S3(5)*(-CI*(P2(0))+CI*(P2(3)))
      COEFF(1,1,2)= 0Q0
      COEFF(2,1,2)= 0Q0
      COEFF(3,1,2)= 0Q0
      COEFF(4,1,2)= COUP*-CI * S3(5)
      COEFF(1,2,2)= 0Q0
      COEFF(2,2,2)= 0Q0
      COEFF(3,2,2)= COUP*-CI * S3(5)
      COEFF(4,2,2)= 0Q0
      COEFF(1,3,2)= 0Q0
      COEFF(2,3,2)= 0Q0
      COEFF(3,3,2)= COUP*-S3(5)
      COEFF(4,3,2)= 0Q0
      COEFF(1,4,2)= 0Q0
      COEFF(2,4,2)= 0Q0
      COEFF(3,4,2)= 0Q0
      COEFF(4,4,2)= COUP*CI * S3(5)
      COEFF(1,0,3)= COUP*S3(5)*(-CI*(P2(0))+CI*(P2(3)))
      COEFF(2,0,3)= COUP*S3(5)*(+CI*(P2(1))-P2(2))
      COEFF(3,0,3)= COUP*CI * S3(5)*M2
      COEFF(4,0,3)= 0Q0
      COEFF(1,1,3)= COUP*-CI * S3(5)
      COEFF(2,1,3)= 0Q0
      COEFF(3,1,3)= 0Q0
      COEFF(4,1,3)= 0Q0
      COEFF(1,2,3)= 0Q0
      COEFF(2,2,3)= COUP*CI * S3(5)
      COEFF(3,2,3)= 0Q0
      COEFF(4,2,3)= 0Q0
      COEFF(1,3,3)= 0Q0
      COEFF(2,3,3)= COUP*-S3(5)
      COEFF(3,3,3)= 0Q0
      COEFF(4,3,3)= 0Q0
      COEFF(1,4,3)= COUP*CI * S3(5)
      COEFF(2,4,3)= 0Q0
      COEFF(3,4,3)= 0Q0
      COEFF(4,4,3)= 0Q0
      COEFF(1,0,4)= COUP*S3(5)*(P2(2)+CI*(P2(1)))
      COEFF(2,0,4)= COUP*-S3(5)*(+CI*(P2(0)+P2(3)))
      COEFF(3,0,4)= 0Q0
      COEFF(4,0,4)= COUP*CI * S3(5)*M2
      COEFF(1,1,4)= 0Q0
      COEFF(2,1,4)= COUP*-CI * S3(5)
      COEFF(3,1,4)= 0Q0
      COEFF(4,1,4)= 0Q0
      COEFF(1,2,4)= COUP*CI * S3(5)
      COEFF(2,2,4)= 0Q0
      COEFF(3,2,4)= 0Q0
      COEFF(4,2,4)= 0Q0
      COEFF(1,3,4)= COUP*S3(5)
      COEFF(2,3,4)= 0Q0
      COEFF(3,3,4)= 0Q0
      COEFF(4,3,4)= 0Q0
      COEFF(1,4,4)= 0Q0
      COEFF(2,4,4)= COUP*-CI * S3(5)
      COEFF(3,4,4)= 0Q0
      COEFF(4,4,4)= 0Q0
      END


