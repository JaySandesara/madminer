C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Gamma(3,2,-1)*ProjM(-1,1) + Gamma(3,2,-1)*ProjP(-1,1)) + Coup(2) * (Gamma(3,2,-1)*ProjP(-1,1)) + Coup(3) * (Gamma(3,2,-1)*ProjM(-1,1) + (2*Gamma(3,2,-1)*ProjP(-1,1))/3.)
C     
      SUBROUTINE MP_FFV10_6_9L1_2(P1, V3, COUP1, COUP2, COUP3, M2, W2,
     $  P2, COEFF)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 V3(*)
      INCLUDE 'coef_specs.inc'
      COMPLEX*32 COEFF(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      REAL*16 W2
      COMPLEX*32 P2(0:3)
      COMPLEX*32 COUP1
      REAL*16 M2
      COMPLEX*32 COUP2
      COMPLEX*32 COUP3
      COMPLEX*32 P1(0:3)
      P2(0) = +P1(0)+V3(1)
      P2(1) = +P1(1)+V3(2)
      P2(2) = +P1(2)+V3(3)
      P2(3) = +P1(3)+V3(4)
      COEFF(1,0,1)= (COUP1*(P2(0)*(-1Q0)*(+CI*(V3(5)+V3(8)))+(P2(1)*(
     $ +CI*(V3(6))-V3(7))+(P2(2)*(V3(6)+CI*(V3(7)))+P2(3)*(+CI*(V3(5)
     $ +V3(8))))))+COUP3*(P2(0)*(-1Q0)*(+CI*(V3(5)+V3(8)))+(P2(1)*(+CI
     $ *(V3(6))-V3(7))+(P2(2)*(V3(6)+CI*(V3(7)))+P2(3)*(+CI*(V3(5)
     $ +V3(8)))))))
      COEFF(2,0,1)= (COUP1*(P2(0)*(V3(7)-CI*(V3(6)))+(P2(1)*(+CI*(V3(5)
     $ +V3(8)))+(P2(2)*(-1Q0)*(V3(5)+V3(8))+P2(3)*(V3(7)-CI*(V3(6))))))
     $ +COUP3*(P2(0)*(V3(7)-CI*(V3(6)))+(P2(1)*(+CI*(V3(5)+V3(8)))
     $ +(P2(2)*(-1Q0)*(V3(5)+V3(8))+P2(3)*(V3(7)-CI*(V3(6)))))))
      COEFF(3,0,1)= M2*(COUP1*(+CI*(V3(5)+V3(8)))+COUP3*(+CI*(V3(5)
     $ +V3(8))))
      COEFF(4,0,1)= M2*(COUP1*(+CI*(V3(6))-V3(7))+COUP3*(+CI*(V3(6))
     $ -V3(7)))
      COEFF(1,1,1)= (-1Q0)*(COUP1*(+CI*(V3(5)+V3(8)))+COUP3*(+CI*(V3(5)
     $ +V3(8))))
      COEFF(2,1,1)= (-1Q0)*(COUP1*(+CI*(V3(6))-V3(7))+COUP3*(+CI*(V3(6)
     $ )-V3(7)))
      COEFF(3,1,1)= 0Q0
      COEFF(4,1,1)= 0Q0
      COEFF(1,2,1)= (COUP1*(+CI*(V3(6))-V3(7))+COUP3*(+CI*(V3(6))-V3(7)
     $ ))
      COEFF(2,2,1)= (COUP1*(+CI*(V3(5)+V3(8)))+COUP3*(+CI*(V3(5)+V3(8))
     $ ))
      COEFF(3,2,1)= 0Q0
      COEFF(4,2,1)= 0Q0
      COEFF(1,3,1)= (COUP1*(V3(6)+CI*(V3(7)))+COUP3*(V3(6)+CI*(V3(7))))
      COEFF(2,3,1)= (COUP1*(-1Q0)*(V3(5)+V3(8))-COUP3*(V3(5)+V3(8)))
      COEFF(3,3,1)= 0Q0
      COEFF(4,3,1)= 0Q0
      COEFF(1,4,1)= (COUP1*(+CI*(V3(5)+V3(8)))+COUP3*(+CI*(V3(5)+V3(8))
     $ ))
      COEFF(2,4,1)= (-1Q0)*(COUP1*(+CI*(V3(6))-V3(7))+COUP3*(+CI*(V3(6)
     $ )-V3(7)))
      COEFF(3,4,1)= 0Q0
      COEFF(4,4,1)= 0Q0
      COEFF(1,0,2)= (COUP1*(P2(0)*(-1Q0)*(V3(7)+CI*(V3(6)))+(P2(1)*(
     $ -CI*(V3(8))+CI*(V3(5)))+(P2(2)*(V3(5)-V3(8))+P2(3)*(V3(7)+CI
     $ *(V3(6))))))+COUP3*(P2(0)*(-1Q0)*(V3(7)+CI*(V3(6)))+(P2(1)*(-CI
     $ *(V3(8))+CI*(V3(5)))+(P2(2)*(V3(5)-V3(8))+P2(3)*(V3(7)+CI*(V3(6)
     $ ))))))
      COEFF(2,0,2)= (COUP1*(P2(0)*(-CI*(V3(5))+CI*(V3(8)))+(P2(1)
     $ *(V3(7)+CI*(V3(6)))+(P2(2)*(+CI*(V3(7))-V3(6))+P2(3)*(-CI*(V3(5)
     $ )+CI*(V3(8))))))+COUP3*(P2(0)*(-CI*(V3(5))+CI*(V3(8)))+(P2(1)
     $ *(V3(7)+CI*(V3(6)))+(P2(2)*(+CI*(V3(7))-V3(6))+P2(3)*(-CI*(V3(5)
     $ )+CI*(V3(8)))))))
      COEFF(3,0,2)= M2*(COUP1*(V3(7)+CI*(V3(6)))+COUP3*(V3(7)+CI*(V3(6)
     $ )))
      COEFF(4,0,2)= M2*(COUP1*(-CI*(V3(8))+CI*(V3(5)))+COUP3*(-CI
     $ *(V3(8))+CI*(V3(5))))
      COEFF(1,1,2)= (-1Q0)*(COUP1*(V3(7)+CI*(V3(6)))+COUP3*(V3(7)+CI
     $ *(V3(6))))
      COEFF(2,1,2)= (COUP1*(-CI*(V3(5))+CI*(V3(8)))+COUP3*(-CI*(V3(5))
     $ +CI*(V3(8))))
      COEFF(3,1,2)= 0Q0
      COEFF(4,1,2)= 0Q0
      COEFF(1,2,2)= (COUP1*(-CI*(V3(8))+CI*(V3(5)))+COUP3*(-CI*(V3(8))
     $ +CI*(V3(5))))
      COEFF(2,2,2)= (COUP1*(V3(7)+CI*(V3(6)))+COUP3*(V3(7)+CI*(V3(6))))
      COEFF(3,2,2)= 0Q0
      COEFF(4,2,2)= 0Q0
      COEFF(1,3,2)= (COUP1*(V3(5)-V3(8))+COUP3*(V3(5)-V3(8)))
      COEFF(2,3,2)= (COUP1*(+CI*(V3(7))-V3(6))+COUP3*(+CI*(V3(7))-V3(6)
     $ ))
      COEFF(3,3,2)= 0Q0
      COEFF(4,3,2)= 0Q0
      COEFF(1,4,2)= (COUP1*(V3(7)+CI*(V3(6)))+COUP3*(V3(7)+CI*(V3(6))))
      COEFF(2,4,2)= (COUP1*(-CI*(V3(5))+CI*(V3(8)))+COUP3*(-CI*(V3(5))
     $ +CI*(V3(8))))
      COEFF(3,4,2)= 0Q0
      COEFF(4,4,2)= 0Q0
      COEFF(1,0,3)= M2*(V3(5)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ -V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(2,0,3)= -M2*(V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ -V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(3,0,3)= (COUP1*(P2(0)*(-CI*(V3(5))+CI*(V3(8)))+(P2(1)*(+CI
     $ *(V3(6))-V3(7))+(P2(2)*(V3(6)+CI*(V3(7)))+P2(3)*(-CI*(V3(5))+CI
     $ *(V3(8))))))+(COUP2*(P2(0)*(-CI*(V3(5))+CI*(V3(8)))+(P2(1)*(+CI
     $ *(V3(6))-V3(7))+(P2(2)*(V3(6)+CI*(V3(7)))+P2(3)*(-CI*(V3(5))+CI
     $ *(V3(8))))))+COUP3*(P2(0)*2Q0/3Q0*(-CI*(V3(5))+CI*(V3(8)))
     $ +(P2(1)*2Q0/3Q0*(+CI*(V3(6))-V3(7))+(P2(2)*2Q0/3Q0*(V3(6)+CI
     $ *(V3(7)))+2Q0/3Q0*(P2(3)*(-CI*(V3(5))+CI*(V3(8)))))))))
      COEFF(4,0,3)= (COUP1*(P2(0)*(+CI*(V3(6))-V3(7))+(P2(1)*(-CI
     $ *(V3(5))+CI*(V3(8)))+(P2(2)*(V3(5)-V3(8))+P2(3)*(V3(7)-CI*(V3(6)
     $ )))))+(COUP2*(P2(0)*(+CI*(V3(6))-V3(7))+(P2(1)*(-CI*(V3(5))+CI
     $ *(V3(8)))+(P2(2)*(V3(5)-V3(8))+P2(3)*(V3(7)-CI*(V3(6))))))
     $ +COUP3*(P2(0)*2Q0/3Q0*(+CI*(V3(6))-V3(7))+(P2(1)*2Q0/3Q0*(-CI
     $ *(V3(5))+CI*(V3(8)))+(P2(2)*2Q0/3Q0*(V3(5)-V3(8))+2Q0/3Q0*(P2(3)
     $ *(V3(7)-CI*(V3(6)))))))))
      COEFF(1,1,3)= 0Q0
      COEFF(2,1,3)= 0Q0
      COEFF(3,1,3)= (V3(5)*(-1Q0)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI
     $ *(COUP3))+V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(4,1,3)= (V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ -V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(1,2,3)= 0Q0
      COEFF(2,2,3)= 0Q0
      COEFF(3,2,3)= (V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ -V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(4,2,3)= (V3(5)*(-1Q0)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI
     $ *(COUP3))+V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(1,3,3)= 0Q0
      COEFF(2,3,3)= 0Q0
      COEFF(3,3,3)= (V3(6)*(COUP1+COUP2+2Q0/3Q0*(COUP3))+V3(7)*(+CI
     $ *(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(4,3,3)= 2Q0/3Q0*(-3Q0/2Q0*(V3(8)*(COUP1+COUP2+2Q0/3Q0
     $ *(COUP3)))+V3(5)*3Q0/2Q0*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(1,4,3)= 0Q0
      COEFF(2,4,3)= 0Q0
      COEFF(3,4,3)= (V3(5)*(-1Q0)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI
     $ *(COUP3))+V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(4,4,3)= (-1Q0)*(V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI
     $ *(COUP3))-V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(1,0,4)= -M2*(V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ +V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(2,0,4)= M2*(V3(5)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ +V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(3,0,4)= (COUP1*(P2(0)*(V3(7)+CI*(V3(6)))+(P2(1)*(-1Q0)*(
     $ +CI*(V3(5)+V3(8)))+(P2(2)*(-1Q0)*(V3(5)+V3(8))+P2(3)*(V3(7)+CI
     $ *(V3(6))))))+(COUP2*(P2(0)*(V3(7)+CI*(V3(6)))+(P2(1)*(-1Q0)*(
     $ +CI*(V3(5)+V3(8)))+(P2(2)*(-1Q0)*(V3(5)+V3(8))+P2(3)*(V3(7)+CI
     $ *(V3(6))))))+COUP3*(P2(0)*2Q0/3Q0*(V3(7)+CI*(V3(6)))+(P2(1)*(
     $ -2Q0/3Q0)*(+CI*(V3(5)+V3(8)))+(P2(2)*(-2Q0/3Q0)*(V3(5)+V3(8))
     $ +2Q0/3Q0*(P2(3)*(V3(7)+CI*(V3(6)))))))))
      COEFF(4,0,4)= (COUP1*(P2(0)*(-1Q0)*(+CI*(V3(5)+V3(8)))+(P2(1)
     $ *(V3(7)+CI*(V3(6)))+(P2(2)*(+CI*(V3(7))-V3(6))+P2(3)*(+CI*(V3(5)
     $ +V3(8))))))+(COUP2*(P2(0)*(-1Q0)*(+CI*(V3(5)+V3(8)))+(P2(1)
     $ *(V3(7)+CI*(V3(6)))+(P2(2)*(+CI*(V3(7))-V3(6))+P2(3)*(+CI*(V3(5)
     $ +V3(8))))))+COUP3*(P2(0)*(-2Q0/3Q0)*(+CI*(V3(5)+V3(8)))+(P2(1)
     $ *2Q0/3Q0*(V3(7)+CI*(V3(6)))+(P2(2)*2Q0/3Q0*(+CI*(V3(7))-V3(6))
     $ +2Q0/3Q0*(P2(3)*(+CI*(V3(5)+V3(8)))))))))
      COEFF(1,1,4)= 0Q0
      COEFF(2,1,4)= 0Q0
      COEFF(3,1,4)= (V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ +V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(4,1,4)= (-1Q0)*(V3(5)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI
     $ *(COUP3))+V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(1,2,4)= 0Q0
      COEFF(2,2,4)= 0Q0
      COEFF(3,2,4)= (-1Q0)*(V3(5)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI
     $ *(COUP3))+V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(4,2,4)= (V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ +V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(1,3,4)= 0Q0
      COEFF(2,3,4)= 0Q0
      COEFF(3,3,4)= (-2Q0/3Q0)*(+3Q0/2Q0*(V3(8)*(COUP1+COUP2+2Q0/3Q0
     $ *(COUP3)))+V3(5)*3Q0/2Q0*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(4,3,4)= (V3(6)*(-1Q0)*(COUP1+COUP2+2Q0/3Q0*(COUP3))+V3(7)
     $ *(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      COEFF(1,4,4)= 0Q0
      COEFF(2,4,4)= 0Q0
      COEFF(3,4,4)= (V3(6)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ +V3(7)*(COUP1+COUP2+2Q0/3Q0*(COUP3)))
      COEFF(4,4,4)= (V3(5)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3))
     $ +V3(8)*(+CI*(COUP1+COUP2)+2Q0/3Q0 * CI*(COUP3)))
      END


