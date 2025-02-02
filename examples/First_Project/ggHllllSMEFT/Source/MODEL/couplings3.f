ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP3()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_457 = (4.000000D+00*MDL_CPG*G*MDL_VEV0)/MDL_LAMBDA__EXP__2
      GC_458 = (MDL_CTG*MDL_COMPLEXI*G*MDL_VEV0)/(MDL_LAMBDA__EXP__2
     $ *MDL_SQRT__2)
      GC_460 = -((MDL_CTG*MDL_G__EXP__2*MDL_VEV0)/(MDL_LAMBDA__EXP__2
     $ *MDL_SQRT__2))
      R2GC_694_147 = (-3.100000D+01*MDL_CPG*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_VEV0)/(1.600000D+01*PI**2*MDL_LAMBDA__EXP__2)
      R2GC_709_162 = -(MDL_COMPLEXI*MDL_G__EXP__2*MDL_MT__EXP__2)
     $ /(8.000000D+00*PI**2*MDL_VEV0)
      R2GC_900_292 = (MDL_C3PL1*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_MT__EXP__2*MDL_VEV0)/(1.600000D+01*PI**2
     $ *MDL_LAMBDA__EXP__2)+(MDL_C3PL2*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_MT__EXP__2*MDL_VEV0)/(1.600000D+01*PI**2
     $ *MDL_LAMBDA__EXP__2)-(MDL_CDP*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_MT__EXP__2*MDL_VEV0)/(8.000000D+00*PI**2
     $ *MDL_LAMBDA__EXP__2)-(MDL_CLL1221*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_MT__EXP__2*MDL_VEV0)/(1.600000D+01*PI**2
     $ *MDL_LAMBDA__EXP__2)+(MDL_CPDC*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_MT__EXP__2*MDL_VEV0)/(3.200000D+01*PI**2
     $ *MDL_LAMBDA__EXP__2)+(MDL_CTP*MDL_COMPLEXI*MDL_G__EXP__2*MDL_MT
     $ *MDL_VEV0__EXP__2)/(8.000000D+00*PI**2*MDL_LAMBDA__EXP__2
     $ *MDL_SQRT__2)
      UVGC_1217_229_1EPS = (MDL_CPG*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_VEV0)/(6.000000D+00*PI**2*MDL_LAMBDA__EXP__2)
      UVGC_1217_230_1EPS = -(MDL_CPG*MDL_COMPLEXI*MDL_G__EXP__2
     $ *MDL_VEV0)/(1.600000D+01*PI**2*MDL_LAMBDA__EXP__2)
      UVGC_991_428_1EPS = (-2.500000D+01*MDL_CPG*MDL_COMPLEXI
     $ *MDL_G__EXP__2*MDL_VEV0)/(1.600000D+01*PI**2*MDL_LAMBDA__EXP__2)
      UVGC_1217_231 = (MDL_CPG*MDL_COMPLEXI*MDL_G__EXP__2*MDL_VEV0
     $ *REGLOG(DCMPLX(MDL_MT__EXP__2/MDL_MU_R__EXP__2)))/(6.000000D+00
     $ *PI**2*MDL_LAMBDA__EXP__2)
      GC_10 = -G
      GC_11 = MDL_COMPLEXI*G
      GC_12 = MDL_COMPLEXI*MDL_G__EXP__2
      GC_197 = (MDL_CTG*MDL_COMPLEXI*G)/(MDL_LAMBDA__EXP__2
     $ *MDL_SQRT__2)
      GC_203 = -((MDL_CTG*MDL_G__EXP__2)/(MDL_LAMBDA__EXP__2
     $ *MDL_SQRT__2))
      R2GC_701_154 = (MDL_CTG*MDL_COMPLEXI*MDL_G__EXP__2*MDL_MT)
     $ /(3.000000D+00*PI**2*MDL_LAMBDA__EXP__2*MDL_SQRT__2)
      UVGC_992_429_1EPS = (MDL_CTG*MDL_COMPLEXI*MDL_G__EXP__2*MDL_MT)
     $ /(PI**2*MDL_LAMBDA__EXP__2*MDL_SQRT__2)
      END
