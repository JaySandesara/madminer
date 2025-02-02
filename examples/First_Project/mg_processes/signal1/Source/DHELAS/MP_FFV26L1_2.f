C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     -(P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1)) +
C      P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) -
C      P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) +
C      P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)
C     
      SUBROUTINE MP_FFV26L1_2(P1, V3, COUP, M2, W2, P2, COEFF)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 V3(*)
      INCLUDE 'coef_specs.inc'
      COMPLEX*32 COEFF(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      REAL*16 W2
      COMPLEX*32 P2(0:3)
      COMPLEX*32 P3(0:3)
      REAL*16 M2
      COMPLEX*32 P1(0:3)
      COMPLEX*32 COUP
      P3(0) = V3(1)
      P3(1) = V3(2)
      P3(2) = V3(3)
      P3(3) = V3(4)
      P2(0) = +P1(0)+V3(1)
      P2(1) = +P1(1)+V3(2)
      P2(2) = +P1(2)+V3(3)
      P2(3) = +P1(3)+V3(4)
      COEFF(1,0,1)= COUP*2Q0 * M2*(V3(6)*P3(2)-CI*(V3(8)*P3(0))+CI
     $ *(V3(5)*P3(3))-V3(7)*P3(1))
      COEFF(2,0,1)= COUP*2Q0 * M2*(P3(0)*(V3(7)-CI*(V3(6)))+(P3(1)*(
     $ +CI*(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)*(V3(7)-CI
     $ *(V3(6))))))
      COEFF(3,0,1)= COUP*2Q0*(P2(1)*(P3(0)*(+CI*(V3(6))-V3(7))+(P3(1)
     $ *(-1Q0)*(+CI*(V3(5)+V3(8)))+(P3(2)*(V3(5)+V3(8))+P3(3)*(+CI
     $ *(V3(6))-V3(7)))))+(P2(2)*(P3(0)*(V3(6)+CI*(V3(7)))+(P3(1)*(
     $ -1Q0)*(V3(5)+V3(8))+(P3(2)*(-1Q0)*(+CI*(V3(5)+V3(8)))+P3(3)
     $ *(V3(6)+CI*(V3(7))))))+(P2(0)*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI
     $ *(V3(8)*P3(0))-V3(6)*P3(2))+P2(3)*(V3(7)*P3(1)-CI*(V3(5)*P3(3))
     $ +CI*(V3(8)*P3(0))-V3(6)*P3(2)))))
      COEFF(4,0,1)= COUP*2Q0*(P2(0)*(P3(0)*(+CI*(V3(6))-V3(7))+(P3(1)
     $ *(-1Q0)*(+CI*(V3(5)+V3(8)))+(P3(2)*(V3(5)+V3(8))+P3(3)*(+CI
     $ *(V3(6))-V3(7)))))+(P2(3)*(P3(0)*(V3(7)-CI*(V3(6)))+(P3(1)*(+CI
     $ *(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)*(V3(7)-CI
     $ *(V3(6))))))+(P2(1)*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(6)*P3(2))+P2(2)*(V3(5)*P3(3)-CI*(V3(6)*P3(2))+CI
     $ *(V3(7)*P3(1))-V3(8)*P3(0)))))
      COEFF(1,1,1)= 0Q0
      COEFF(2,1,1)= 0Q0
      COEFF(3,1,1)= COUP*2Q0*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(6)*P3(2))
      COEFF(4,1,1)= COUP*2Q0*(P3(0)*(+CI*(V3(6))-V3(7))+(P3(1)*(-1Q0)
     $ *(+CI*(V3(5)+V3(8)))+(P3(2)*(V3(5)+V3(8))+P3(3)*(+CI*(V3(6))
     $ -V3(7)))))
      COEFF(1,2,1)= 0Q0
      COEFF(2,2,1)= 0Q0
      COEFF(3,2,1)= COUP*2Q0*(P3(0)*(+CI*(V3(6))-V3(7))+(P3(1)*(-1Q0)
     $ *(+CI*(V3(5)+V3(8)))+(P3(2)*(V3(5)+V3(8))+P3(3)*(+CI*(V3(6))
     $ -V3(7)))))
      COEFF(4,2,1)= COUP*2Q0*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(6)*P3(2))
      COEFF(1,3,1)= 0Q0
      COEFF(2,3,1)= 0Q0
      COEFF(3,3,1)= COUP*2Q0*(P3(0)*(V3(6)+CI*(V3(7)))+(P3(1)*(-1Q0)
     $ *(V3(5)+V3(8))+(P3(2)*(-1Q0)*(+CI*(V3(5)+V3(8)))+P3(3)*(V3(6)
     $ +CI*(V3(7))))))
      COEFF(4,3,1)= COUP*2Q0*(V3(5)*P3(3)-CI*(V3(6)*P3(2))+CI*(V3(7)
     $ *P3(1))-V3(8)*P3(0))
      COEFF(1,4,1)= 0Q0
      COEFF(2,4,1)= 0Q0
      COEFF(3,4,1)= COUP*2Q0*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(6)*P3(2))
      COEFF(4,4,1)= COUP*2Q0*(P3(0)*(V3(7)-CI*(V3(6)))+(P3(1)*(+CI
     $ *(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)*(V3(7)-CI
     $ *(V3(6))))))
      COEFF(1,0,2)= COUP*2Q0 * M2*(P3(0)*(-1Q0)*(V3(7)+CI*(V3(6)))
     $ +(P3(1)*(-CI*(V3(8))+CI*(V3(5)))+(P3(2)*(V3(5)-V3(8))+P3(3)
     $ *(V3(7)+CI*(V3(6))))))
      COEFF(2,0,2)= COUP*2Q0 * M2*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI
     $ *(V3(8)*P3(0))-V3(6)*P3(2))
      COEFF(3,0,2)= COUP*2Q0*(P2(0)*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(
     $ -CI*(V3(5))+CI*(V3(8)))+(P3(2)*(V3(8)-V3(5))-P3(3)*(V3(7)+CI
     $ *(V3(6))))))+(P2(3)*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(-CI*(V3(5)
     $ )+CI*(V3(8)))+(P3(2)*(V3(8)-V3(5))-P3(3)*(V3(7)+CI*(V3(6))))))
     $ +(P2(1)*(V3(6)*P3(2)-CI*(V3(8)*P3(0))+CI*(V3(5)*P3(3))-V3(7)
     $ *P3(1))+P2(2)*(V3(5)*P3(3)-CI*(V3(6)*P3(2))+CI*(V3(7)*P3(1))
     $ -V3(8)*P3(0)))))
      COEFF(4,0,2)= COUP*2Q0*(P2(1)*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(
     $ -CI*(V3(5))+CI*(V3(8)))+(P3(2)*(V3(8)-V3(5))-P3(3)*(V3(7)+CI
     $ *(V3(6))))))+(P2(2)*(P3(0)*(+CI*(V3(7))-V3(6))+(P3(1)*(V3(5)
     $ -V3(8))+(P3(2)*(-CI*(V3(5))+CI*(V3(8)))+P3(3)*(V3(6)-CI*(V3(7)))
     $ )))+(P2(0)*(V3(6)*P3(2)-CI*(V3(8)*P3(0))+CI*(V3(5)*P3(3))-V3(7)
     $ *P3(1))+P2(3)*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI*(V3(8)*P3(0))
     $ -V3(6)*P3(2)))))
      COEFF(1,1,2)= 0Q0
      COEFF(2,1,2)= 0Q0
      COEFF(3,1,2)= COUP*2Q0*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(-CI
     $ *(V3(5))+CI*(V3(8)))+(P3(2)*(V3(8)-V3(5))-P3(3)*(V3(7)+CI*(V3(6)
     $ )))))
      COEFF(4,1,2)= COUP*2Q0*(V3(6)*P3(2)-CI*(V3(8)*P3(0))+CI*(V3(5)
     $ *P3(3))-V3(7)*P3(1))
      COEFF(1,2,2)= 0Q0
      COEFF(2,2,2)= 0Q0
      COEFF(3,2,2)= COUP*2Q0*(V3(6)*P3(2)-CI*(V3(8)*P3(0))+CI*(V3(5)
     $ *P3(3))-V3(7)*P3(1))
      COEFF(4,2,2)= COUP*2Q0*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(-CI
     $ *(V3(5))+CI*(V3(8)))+(P3(2)*(V3(8)-V3(5))-P3(3)*(V3(7)+CI*(V3(6)
     $ )))))
      COEFF(1,3,2)= 0Q0
      COEFF(2,3,2)= 0Q0
      COEFF(3,3,2)= COUP*2Q0*(V3(5)*P3(3)-CI*(V3(6)*P3(2))+CI*(V3(7)
     $ *P3(1))-V3(8)*P3(0))
      COEFF(4,3,2)= COUP*2Q0*(P3(0)*(+CI*(V3(7))-V3(6))+(P3(1)*(V3(5)
     $ -V3(8))+(P3(2)*(-CI*(V3(5))+CI*(V3(8)))+P3(3)*(V3(6)-CI*(V3(7)))
     $ )))
      COEFF(1,4,2)= 0Q0
      COEFF(2,4,2)= 0Q0
      COEFF(3,4,2)= COUP*2Q0*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(-CI
     $ *(V3(5))+CI*(V3(8)))+(P3(2)*(V3(8)-V3(5))-P3(3)*(V3(7)+CI*(V3(6)
     $ )))))
      COEFF(4,4,2)= COUP*2Q0*(V3(7)*P3(1)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(6)*P3(2))
      COEFF(1,0,3)= COUP*2Q0*(P2(1)*(P3(0)*(+CI*(V3(6))-V3(7))+(P3(1)
     $ *(-CI*(V3(5))+CI*(V3(8)))+(P3(2)*(V3(5)-V3(8))+P3(3)*(V3(7)-CI
     $ *(V3(6))))))+(P2(2)*(P3(0)*(V3(6)+CI*(V3(7)))+(P3(1)*(V3(8)
     $ -V3(5))+(P3(2)*(-CI*(V3(5))+CI*(V3(8)))-P3(3)*(V3(6)+CI*(V3(7)))
     $ )))+(P2(0)*(V3(7)*P3(1)-CI*(V3(8)*P3(0))+CI*(V3(5)*P3(3))-V3(6)
     $ *P3(2))+P2(3)*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI*(V3(8)*P3(0))
     $ -V3(7)*P3(1)))))
      COEFF(2,0,3)= COUP*2Q0*(P2(0)*(P3(0)*(V3(7)-CI*(V3(6)))+(P3(1)*(
     $ -CI*(V3(8))+CI*(V3(5)))+(P3(2)*(V3(8)-V3(5))+P3(3)*(+CI*(V3(6))
     $ -V3(7)))))+(P2(3)*(P3(0)*(V3(7)-CI*(V3(6)))+(P3(1)*(-CI*(V3(8))
     $ +CI*(V3(5)))+(P3(2)*(V3(8)-V3(5))+P3(3)*(+CI*(V3(6))-V3(7)))))
     $ +(P2(1)*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI*(V3(8)*P3(0))-V3(7)
     $ *P3(1))+P2(2)*(V3(5)*P3(3)-CI*(V3(7)*P3(1))+CI*(V3(6)*P3(2))
     $ -V3(8)*P3(0)))))
      COEFF(3,0,3)= COUP*2Q0 * M2*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI
     $ *(V3(8)*P3(0))-V3(7)*P3(1))
      COEFF(4,0,3)= COUP*2Q0 * M2*(P3(0)*(+CI*(V3(6))-V3(7))+(P3(1)*(
     $ -CI*(V3(5))+CI*(V3(8)))+(P3(2)*(V3(5)-V3(8))+P3(3)*(V3(7)-CI
     $ *(V3(6))))))
      COEFF(1,1,3)= COUP*2Q0*(V3(7)*P3(1)-CI*(V3(8)*P3(0))+CI*(V3(5)
     $ *P3(3))-V3(6)*P3(2))
      COEFF(2,1,3)= COUP*2Q0*(P3(0)*(V3(7)-CI*(V3(6)))+(P3(1)*(-CI
     $ *(V3(8))+CI*(V3(5)))+(P3(2)*(V3(8)-V3(5))+P3(3)*(+CI*(V3(6))
     $ -V3(7)))))
      COEFF(3,1,3)= 0Q0
      COEFF(4,1,3)= 0Q0
      COEFF(1,2,3)= COUP*2Q0*(P3(0)*(+CI*(V3(6))-V3(7))+(P3(1)*(-CI
     $ *(V3(5))+CI*(V3(8)))+(P3(2)*(V3(5)-V3(8))+P3(3)*(V3(7)-CI*(V3(6)
     $ )))))
      COEFF(2,2,3)= COUP*2Q0*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(7)*P3(1))
      COEFF(3,2,3)= 0Q0
      COEFF(4,2,3)= 0Q0
      COEFF(1,3,3)= COUP*2Q0*(P3(0)*(V3(6)+CI*(V3(7)))+(P3(1)*(V3(8)
     $ -V3(5))+(P3(2)*(-CI*(V3(5))+CI*(V3(8)))-P3(3)*(V3(6)+CI*(V3(7)))
     $ )))
      COEFF(2,3,3)= COUP*2Q0*(V3(5)*P3(3)-CI*(V3(7)*P3(1))+CI*(V3(6)
     $ *P3(2))-V3(8)*P3(0))
      COEFF(3,3,3)= 0Q0
      COEFF(4,3,3)= 0Q0
      COEFF(1,4,3)= COUP*2Q0*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(7)*P3(1))
      COEFF(2,4,3)= COUP*2Q0*(P3(0)*(V3(7)-CI*(V3(6)))+(P3(1)*(-CI
     $ *(V3(8))+CI*(V3(5)))+(P3(2)*(V3(8)-V3(5))+P3(3)*(+CI*(V3(6))
     $ -V3(7)))))
      COEFF(3,4,3)= 0Q0
      COEFF(4,4,3)= 0Q0
      COEFF(1,0,4)= COUP*2Q0*(P2(0)*(P3(0)*(-1Q0)*(V3(7)+CI*(V3(6)))
     $ +(P3(1)*(+CI*(V3(5)+V3(8)))+(P3(2)*(V3(5)+V3(8))-P3(3)*(V3(7)
     $ +CI*(V3(6))))))+(P2(3)*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(-1Q0)
     $ *(+CI*(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)*(V3(7)
     $ +CI*(V3(6))))))+(P2(1)*(V3(7)*P3(1)-CI*(V3(8)*P3(0))+CI*(V3(5)
     $ *P3(3))-V3(6)*P3(2))+P2(2)*(V3(5)*P3(3)-CI*(V3(7)*P3(1))+CI
     $ *(V3(6)*P3(2))-V3(8)*P3(0)))))
      COEFF(2,0,4)= COUP*2Q0*(P2(1)*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(
     $ -1Q0)*(+CI*(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)
     $ *(V3(7)+CI*(V3(6))))))+(P2(2)*(P3(0)*(+CI*(V3(7))-V3(6))+(P3(1)
     $ *(V3(5)+V3(8))+(P3(2)*(-1Q0)*(+CI*(V3(5)+V3(8)))+P3(3)*(+CI
     $ *(V3(7))-V3(6)))))+(P2(0)*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI
     $ *(V3(8)*P3(0))-V3(7)*P3(1))+P2(3)*(V3(6)*P3(2)-CI*(V3(5)*P3(3))
     $ +CI*(V3(8)*P3(0))-V3(7)*P3(1)))))
      COEFF(3,0,4)= COUP*2Q0 * M2*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(
     $ -1Q0)*(+CI*(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)
     $ *(V3(7)+CI*(V3(6))))))
      COEFF(4,0,4)= COUP*2Q0 * M2*(V3(7)*P3(1)-CI*(V3(8)*P3(0))+CI
     $ *(V3(5)*P3(3))-V3(6)*P3(2))
      COEFF(1,1,4)= COUP*2Q0*(P3(0)*(-1Q0)*(V3(7)+CI*(V3(6)))+(P3(1)*(
     $ +CI*(V3(5)+V3(8)))+(P3(2)*(V3(5)+V3(8))-P3(3)*(V3(7)+CI*(V3(6)))
     $ )))
      COEFF(2,1,4)= COUP*2Q0*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(7)*P3(1))
      COEFF(3,1,4)= 0Q0
      COEFF(4,1,4)= 0Q0
      COEFF(1,2,4)= COUP*2Q0*(V3(7)*P3(1)-CI*(V3(8)*P3(0))+CI*(V3(5)
     $ *P3(3))-V3(6)*P3(2))
      COEFF(2,2,4)= COUP*2Q0*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(-1Q0)*(
     $ +CI*(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)*(V3(7)+CI
     $ *(V3(6))))))
      COEFF(3,2,4)= 0Q0
      COEFF(4,2,4)= 0Q0
      COEFF(1,3,4)= COUP*2Q0*(V3(5)*P3(3)-CI*(V3(7)*P3(1))+CI*(V3(6)
     $ *P3(2))-V3(8)*P3(0))
      COEFF(2,3,4)= COUP*2Q0*(P3(0)*(+CI*(V3(7))-V3(6))+(P3(1)*(V3(5)
     $ +V3(8))+(P3(2)*(-1Q0)*(+CI*(V3(5)+V3(8)))+P3(3)*(+CI*(V3(7))
     $ -V3(6)))))
      COEFF(3,3,4)= 0Q0
      COEFF(4,3,4)= 0Q0
      COEFF(1,4,4)= COUP*2Q0*(P3(0)*(V3(7)+CI*(V3(6)))+(P3(1)*(-1Q0)*(
     $ +CI*(V3(5)+V3(8)))+(P3(2)*(-1Q0)*(V3(5)+V3(8))+P3(3)*(V3(7)+CI
     $ *(V3(6))))))
      COEFF(2,4,4)= COUP*2Q0*(V3(6)*P3(2)-CI*(V3(5)*P3(3))+CI*(V3(8)
     $ *P3(0))-V3(7)*P3(1))
      COEFF(3,4,4)= 0Q0
      COEFF(4,4,4)= 0Q0
      END


