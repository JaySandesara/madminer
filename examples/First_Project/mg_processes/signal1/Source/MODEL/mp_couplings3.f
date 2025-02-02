ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MP_COUP3()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      MP__GC_457 = (4.000000E+00_16*MP__MDL_CPG*MP__G*MP__MDL_VEV0)
     $ /MP__MDL_LAMBDA__EXP__2
      MP__GC_458 = (MP__MDL_CTG*MP__MDL_COMPLEXI*MP__G*MP__MDL_VEV0)
     $ /(MP__MDL_LAMBDA__EXP__2*MP__MDL_SQRT__2)
      MP__GC_460 = -((MP__MDL_CTG*MP__MDL_G__EXP__2*MP__MDL_VEV0)
     $ /(MP__MDL_LAMBDA__EXP__2*MP__MDL_SQRT__2))
      MP__R2GC_694_147 = (-3.100000E+01_16*MP__MDL_CPG
     $ *MP__MDL_COMPLEXI*MP__MDL_G__EXP__2*MP__MDL_VEV0)/(1.600000E
     $ +01_16*MP__PI**2*MP__MDL_LAMBDA__EXP__2)
      MP__R2GC_709_162 = -(MP__MDL_COMPLEXI*MP__MDL_G__EXP__2
     $ *MP__MDL_MT__EXP__2)/(8.000000E+00_16*MP__PI**2*MP__MDL_VEV0)
      MP__R2GC_900_292 = (MP__MDL_C3PL1*MP__MDL_COMPLEXI
     $ *MP__MDL_G__EXP__2*MP__MDL_MT__EXP__2*MP__MDL_VEV0)/(1.600000E
     $ +01_16*MP__PI**2*MP__MDL_LAMBDA__EXP__2)+(MP__MDL_C3PL2
     $ *MP__MDL_COMPLEXI*MP__MDL_G__EXP__2*MP__MDL_MT__EXP__2
     $ *MP__MDL_VEV0)/(1.600000E+01_16*MP__PI**2*MP__MDL_LAMBDA__EXP__2)
     $ -(MP__MDL_CDP*MP__MDL_COMPLEXI*MP__MDL_G__EXP__2
     $ *MP__MDL_MT__EXP__2*MP__MDL_VEV0)/(8.000000E+00_16*MP__PI**2
     $ *MP__MDL_LAMBDA__EXP__2)-(MP__MDL_CLL1221*MP__MDL_COMPLEXI
     $ *MP__MDL_G__EXP__2*MP__MDL_MT__EXP__2*MP__MDL_VEV0)/(1.600000E
     $ +01_16*MP__PI**2*MP__MDL_LAMBDA__EXP__2)+(MP__MDL_CPDC
     $ *MP__MDL_COMPLEXI*MP__MDL_G__EXP__2*MP__MDL_MT__EXP__2
     $ *MP__MDL_VEV0)/(3.200000E+01_16*MP__PI**2*MP__MDL_LAMBDA__EXP__2)
     $ +(MP__MDL_CTP*MP__MDL_COMPLEXI*MP__MDL_G__EXP__2*MP__MDL_MT
     $ *MP__MDL_VEV0__EXP__2)/(8.000000E+00_16*MP__PI**2
     $ *MP__MDL_LAMBDA__EXP__2*MP__MDL_SQRT__2)
      MP__UVGC_1217_229_1EPS = (MP__MDL_CPG*MP__MDL_COMPLEXI
     $ *MP__MDL_G__EXP__2*MP__MDL_VEV0)/(6.000000E+00_16*MP__PI**2
     $ *MP__MDL_LAMBDA__EXP__2)
      MP__UVGC_1217_230_1EPS = -(MP__MDL_CPG*MP__MDL_COMPLEXI
     $ *MP__MDL_G__EXP__2*MP__MDL_VEV0)/(1.600000E+01_16*MP__PI**2
     $ *MP__MDL_LAMBDA__EXP__2)
      MP__UVGC_991_428_1EPS = (-2.500000E+01_16*MP__MDL_CPG
     $ *MP__MDL_COMPLEXI*MP__MDL_G__EXP__2*MP__MDL_VEV0)/(1.600000E
     $ +01_16*MP__PI**2*MP__MDL_LAMBDA__EXP__2)
      MP__UVGC_1217_231 = (MP__MDL_CPG*MP__MDL_COMPLEXI
     $ *MP__MDL_G__EXP__2*MP__MDL_VEV0*MP_REGLOG(CMPLX((MP__MDL_MT__EXP
     $__2/MP__MDL_MU_R__EXP__2),KIND=16)))/(6.000000E+00_16*MP__PI**2
     $ *MP__MDL_LAMBDA__EXP__2)
      MP__GC_10 = -MP__G
      MP__GC_11 = MP__MDL_COMPLEXI*MP__G
      MP__GC_12 = MP__MDL_COMPLEXI*MP__MDL_G__EXP__2
      MP__GC_197 = (MP__MDL_CTG*MP__MDL_COMPLEXI*MP__G)
     $ /(MP__MDL_LAMBDA__EXP__2*MP__MDL_SQRT__2)
      MP__GC_203 = -((MP__MDL_CTG*MP__MDL_G__EXP__2)
     $ /(MP__MDL_LAMBDA__EXP__2*MP__MDL_SQRT__2))
      MP__R2GC_701_154 = (MP__MDL_CTG*MP__MDL_COMPLEXI
     $ *MP__MDL_G__EXP__2*MP__MDL_MT)/(3.000000E+00_16*MP__PI**2
     $ *MP__MDL_LAMBDA__EXP__2*MP__MDL_SQRT__2)
      MP__UVGC_992_429_1EPS = (MP__MDL_CTG*MP__MDL_COMPLEXI
     $ *MP__MDL_G__EXP__2*MP__MDL_MT)/(MP__PI**2*MP__MDL_LAMBDA__EXP__2
     $ *MP__MDL_SQRT__2)
      END
