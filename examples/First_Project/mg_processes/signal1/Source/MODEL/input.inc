ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DOUBLE PRECISION MDL_SQRT__AS,MDL_G__EXP__2,MDL_G__EXP__3
     $ ,MDL_G__EXP__4,MDL_MU_R__EXP__2,MDL_CQE,MDL_CQMU
     $ ,MDL_CQQ11INTERNAL,MDL_CQQ13INTERNAL,MDL_CQQ81INTERNAL
     $ ,MDL_CQQ83INTERNAL,MDL_CQTA,MDL_CTE,MDL_CTMU,MDL_CTTA
     $ ,MDL_CTU1INTERNAL,MDL_CTU8INTERNAL,MDL_CW0,MDL_SQRT__2,MDL_MUH0
     $ ,MDL_MW__EXP__2,MDL_MZ__EXP__2,MDL_SW0,MDL_C3PQ3INTERNAL
     $ ,MDL_C3PQIINTERNAL,MDL_CF3PL1,MDL_CF3PL2,MDL_CFLL
     $ ,MDL_CPQ3INTERNAL,MDL_CPQIINTERNAL,MDL_NB__2__EXP__0_25
     $ ,MDL_LAMBDA__EXP__2,MDL_CTB,MDL_CW0__EXP__2,MDL_SW0__EXP__2
     $ ,MDL_MH__EXP__2,MDL_SW0__EXP__3,MDL_CW0__EXP__3,MDL_MT__EXP__2
     $ ,MDL_VEV0,MDL_VEV0__EXP__2,MDL_DT,MDL_DV,MDL_DWB,MDL_DLAM
     $ ,MDL_EE0,MDL_CW,MDL_MUH,MDL_SW,MDL_G1,MDL_GW,MDL_LAM,MDL_VEV
     $ ,MDL_EE,MDL_EE__EXP__2,MDL_AEW,AEWM1,MDL_EE0__EXP__2
     $ ,MDL_EE0__EXP__3,MDL_VEV0__EXP__3,MDL_LAMBDA,MDL_CPDC,MDL_CPWB
     $ ,MDL_CDP,MDL_CP,MDL_CWWW,MDL_CG,MDL_CPG,MDL_CPW,MDL_CPBB
     $ ,MDL_CPL1,MDL_CPL2,MDL_CPL3,MDL_C3PL1,MDL_C3PL2,MDL_C3PL3
     $ ,MDL_CPE,MDL_CPMU,MDL_CPTA,MDL_CPQMI,MDL_CPQ3I,MDL_CPQ3
     $ ,MDL_CPQM,MDL_CPU,MDL_CPT,MDL_CPD,MDL_CTP,MDL_CTZ,MDL_CTW
     $ ,MDL_CTG,MDL_CQQ83,MDL_CQQ81,MDL_CQU8,MDL_CTQ8,MDL_CQD8
     $ ,MDL_CTU8,MDL_CTD8,MDL_CQQ13,MDL_CQQ11,MDL_CQU1,MDL_CTQ1
     $ ,MDL_CQD1,MDL_CTU1,MDL_CTD1,MDL_CQQ8,MDL_CQQ1,MDL_CQT1,MDL_CTT1
     $ ,MDL_CQT8,MDL_CQLM1,MDL_CQLM2,MDL_CQL31,MDL_CQL32,MDL_CQE1
     $ ,MDL_CQE2,MDL_CTL1,MDL_CTL2,MDL_CTE1,MDL_CTE2,MDL_CQLM3
     $ ,MDL_CQL33,MDL_CQE3,MDL_CTL3,MDL_CTE3,MDL_CTLS3,MDL_CTLT3
     $ ,MDL_CBLS3,MDL_CLL1111,MDL_CLL2222,MDL_CLL3333,MDL_CLL1122
     $ ,MDL_CLL1133,MDL_CLL2233,MDL_CLL1221,MDL_CLL1331,MDL_CLL2332
     $ ,MDL_MRENO,MDL_GF,AS,MDL_YMT

      COMMON/PARAMS_R/ MDL_SQRT__AS,MDL_G__EXP__2,MDL_G__EXP__3
     $ ,MDL_G__EXP__4,MDL_MU_R__EXP__2,MDL_CQE,MDL_CQMU
     $ ,MDL_CQQ11INTERNAL,MDL_CQQ13INTERNAL,MDL_CQQ81INTERNAL
     $ ,MDL_CQQ83INTERNAL,MDL_CQTA,MDL_CTE,MDL_CTMU,MDL_CTTA
     $ ,MDL_CTU1INTERNAL,MDL_CTU8INTERNAL,MDL_CW0,MDL_SQRT__2,MDL_MUH0
     $ ,MDL_MW__EXP__2,MDL_MZ__EXP__2,MDL_SW0,MDL_C3PQ3INTERNAL
     $ ,MDL_C3PQIINTERNAL,MDL_CF3PL1,MDL_CF3PL2,MDL_CFLL
     $ ,MDL_CPQ3INTERNAL,MDL_CPQIINTERNAL,MDL_NB__2__EXP__0_25
     $ ,MDL_LAMBDA__EXP__2,MDL_CTB,MDL_CW0__EXP__2,MDL_SW0__EXP__2
     $ ,MDL_MH__EXP__2,MDL_SW0__EXP__3,MDL_CW0__EXP__3,MDL_MT__EXP__2
     $ ,MDL_VEV0,MDL_VEV0__EXP__2,MDL_DT,MDL_DV,MDL_DWB,MDL_DLAM
     $ ,MDL_EE0,MDL_CW,MDL_MUH,MDL_SW,MDL_G1,MDL_GW,MDL_LAM,MDL_VEV
     $ ,MDL_EE,MDL_EE__EXP__2,MDL_AEW,AEWM1,MDL_EE0__EXP__2
     $ ,MDL_EE0__EXP__3,MDL_VEV0__EXP__3,MDL_LAMBDA,MDL_CPDC,MDL_CPWB
     $ ,MDL_CDP,MDL_CP,MDL_CWWW,MDL_CG,MDL_CPG,MDL_CPW,MDL_CPBB
     $ ,MDL_CPL1,MDL_CPL2,MDL_CPL3,MDL_C3PL1,MDL_C3PL2,MDL_C3PL3
     $ ,MDL_CPE,MDL_CPMU,MDL_CPTA,MDL_CPQMI,MDL_CPQ3I,MDL_CPQ3
     $ ,MDL_CPQM,MDL_CPU,MDL_CPT,MDL_CPD,MDL_CTP,MDL_CTZ,MDL_CTW
     $ ,MDL_CTG,MDL_CQQ83,MDL_CQQ81,MDL_CQU8,MDL_CTQ8,MDL_CQD8
     $ ,MDL_CTU8,MDL_CTD8,MDL_CQQ13,MDL_CQQ11,MDL_CQU1,MDL_CTQ1
     $ ,MDL_CQD1,MDL_CTU1,MDL_CTD1,MDL_CQQ8,MDL_CQQ1,MDL_CQT1,MDL_CTT1
     $ ,MDL_CQT8,MDL_CQLM1,MDL_CQLM2,MDL_CQL31,MDL_CQL32,MDL_CQE1
     $ ,MDL_CQE2,MDL_CTL1,MDL_CTL2,MDL_CTE1,MDL_CTE2,MDL_CQLM3
     $ ,MDL_CQL33,MDL_CQE3,MDL_CTL3,MDL_CTE3,MDL_CTLS3,MDL_CTLT3
     $ ,MDL_CBLS3,MDL_CLL1111,MDL_CLL2222,MDL_CLL3333,MDL_CLL1122
     $ ,MDL_CLL1133,MDL_CLL2233,MDL_CLL1221,MDL_CLL1331,MDL_CLL2332
     $ ,MDL_MRENO,MDL_GF,AS,MDL_YMT


      DOUBLE COMPLEX MDL_COMPLEXI

      COMMON/PARAMS_C/ MDL_COMPLEXI


