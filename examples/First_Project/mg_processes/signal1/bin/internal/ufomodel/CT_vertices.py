# This file was automatically created by FeynRules 2.4.72
# Mathematica version: 11.3.0 for Mac OS X x86 (64-bit) (March 7, 2018)
# Date: Thu 8 Aug 2019 15:29:47


from object_library import all_vertices, all_CTvertices, Vertex, CTVertex
import particles as P
import CT_couplings as C
import lorentz as L


V_1 = CTVertex(name = 'V_1',
               type = 'R2',
               particles = [ P.g, P.g, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.VVSS1, L.VVSS10, L.VVSS4, L.VVSS6 ],
               loop_particles = [ [ [P.g] ], [ [P.t] ] ],
               couplings = {(0,1,0):C.R2GC_693_146,(0,2,1):[ C.R2GC_898_290, C.R2GC_705_158 ],(0,3,1):C.R2GC_706_159,(0,0,1):C.R2GC_708_161})

V_2 = CTVertex(name = 'V_2',
               type = 'R2',
               particles = [ P.g, P.g, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.VVSS1, L.VVSS10, L.VVSS2, L.VVSS4, L.VVSS7 ],
               loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.g] ], [ [P.t] ] ],
               couplings = {(0,1,4):C.R2GC_693_146,(0,2,0):C.R2GC_905_302,(0,2,2):C.R2GC_905_303,(0,2,3):C.R2GC_905_304,(0,2,5):C.R2GC_905_305,(0,3,5):C.R2GC_956_412,(0,3,1):[ C.R2GC_956_413, C.R2GC_705_158 ],(0,4,1):C.R2GC_782_195,(0,0,1):C.R2GC_708_161})

V_3 = CTVertex(name = 'V_3',
               type = 'R2',
               particles = [ P.g, P.g, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.VVSS1, L.VVSS10, L.VVSS4, L.VVSS6 ],
               loop_particles = [ [ [P.g] ], [ [P.t] ] ],
               couplings = {(0,1,0):C.R2GC_693_146,(0,2,1):[ C.R2GC_899_291, C.R2GC_705_158 ],(0,3,1):C.R2GC_706_159,(0,0,1):C.R2GC_708_161})

V_4 = CTVertex(name = 'V_4',
               type = 'R2',
               particles = [ P.g, P.g, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.VVS3, L.VVS6, L.VVS8 ],
               loop_particles = [ [ [P.g] ], [ [P.t] ] ],
               couplings = {(0,2,0):C.R2GC_694_147,(0,0,1):[ C.R2GC_709_162, C.R2GC_900_292 ],(0,1,1):C.R2GC_701_154})

V_5 = CTVertex(name = 'V_5',
               type = 'R2',
               particles = [ P.g, P.g, P.g ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVV3 ],
               loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ], [ [P.t] ] ],
               couplings = {(0,0,0):C.R2GC_1211_65,(0,0,1):C.R2GC_1211_66,(0,0,2):C.R2GC_1224_83})

V_6 = CTVertex(name = 'V_6',
               type = 'R2',
               particles = [ P.g, P.g, P.g, P.H ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVVS9 ],
               loop_particles = [ [ [P.g] ], [ [P.t] ] ],
               couplings = {(0,0,0):C.R2GC_1218_75,(0,0,1):C.R2GC_1218_76})

V_7 = CTVertex(name = 'V_7',
               type = 'R2',
               particles = [ P.g, P.g, P.g, P.g ],
               color = [ 'd(-1,1,3)*d(-1,2,4)', 'd(-1,1,3)*f(-1,2,4)', 'd(-1,1,4)*d(-1,2,3)', 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,3)*f(-1,1,4)', 'd(-1,2,4)*f(-1,1,3)', 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)', 'Identity(1,2)*Identity(3,4)', 'Identity(1,3)*Identity(2,4)', 'Identity(1,4)*Identity(2,3)' ],
               lorentz = [ L.VVVV11, L.VVVV2, L.VVVV3, L.VVVV5 ],
               loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ], [ [P.t] ] ],
               couplings = {(2,1,0):C.R2GC_1111_5,(2,1,1):C.R2GC_1111_6,(0,1,0):C.R2GC_1111_5,(0,1,1):C.R2GC_1111_6,(4,1,0):C.R2GC_1109_1,(4,1,1):C.R2GC_1109_2,(3,1,0):C.R2GC_1109_1,(3,1,1):C.R2GC_1109_2,(8,1,0):C.R2GC_1110_3,(8,1,1):C.R2GC_1110_4,(6,1,0):C.R2GC_1215_73,(6,1,1):C.R2GC_1215_74,(6,1,2):C.R2GC_1226_85,(7,1,0):C.R2GC_1214_71,(7,1,1):C.R2GC_1214_72,(7,1,2):C.R2GC_1226_85,(5,1,0):C.R2GC_1109_1,(5,1,1):C.R2GC_1109_2,(1,1,0):C.R2GC_1109_1,(1,1,1):C.R2GC_1109_2,(11,0,0):C.R2GC_1113_8,(11,0,1):C.R2GC_1113_9,(10,0,0):C.R2GC_1113_8,(10,0,1):C.R2GC_1113_9,(9,0,1):C.R2GC_1112_7,(2,2,0):C.R2GC_1111_5,(2,2,1):C.R2GC_1111_6,(0,2,0):C.R2GC_1111_5,(0,2,1):C.R2GC_1111_6,(4,2,0):C.R2GC_1109_1,(4,2,1):C.R2GC_1109_2,(3,2,0):C.R2GC_1109_1,(3,2,1):C.R2GC_1109_2,(8,2,0):C.R2GC_1214_71,(8,2,1):C.R2GC_1214_72,(8,2,2):C.R2GC_1226_85,(6,2,0):C.R2GC_1213_69,(6,2,1):C.R2GC_1213_70,(6,2,2):C.R2GC_1225_84,(5,2,0):C.R2GC_1109_1,(5,2,1):C.R2GC_1109_2,(1,2,0):C.R2GC_1109_1,(1,2,1):C.R2GC_1109_2,(7,2,0):C.R2GC_1110_3,(7,2,1):C.R2GC_1110_4,(0,3,0):C.R2GC_1111_5,(0,3,1):C.R2GC_1111_6,(2,3,0):C.R2GC_1111_5,(2,3,1):C.R2GC_1111_6,(5,3,0):C.R2GC_1109_1,(5,3,1):C.R2GC_1109_2,(1,3,0):C.R2GC_1109_1,(1,3,1):C.R2GC_1109_2,(7,3,0):C.R2GC_1212_67,(7,3,1):C.R2GC_1212_68,(7,3,2):C.R2GC_1225_84,(4,3,0):C.R2GC_1109_1,(4,3,1):C.R2GC_1109_2,(3,3,0):C.R2GC_1109_1,(3,3,1):C.R2GC_1109_2,(8,3,0):C.R2GC_1212_67,(8,3,1):C.R2GC_1212_68,(8,3,2):C.R2GC_1225_84})

V_8 = CTVertex(name = 'V_8',
               type = 'R2',
               particles = [ P.g, P.g, P.g, P.g, P.H ],
               color = [ 'd(-1,1,3)*d(-1,2,4)', 'd(-1,1,3)*f(-1,2,4)', 'd(-1,1,4)*d(-1,2,3)', 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,3)*f(-1,1,4)', 'd(-1,2,4)*f(-1,1,3)', 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)', 'Identity(1,2)*Identity(3,4)', 'Identity(1,3)*Identity(2,4)', 'Identity(1,4)*Identity(2,3)' ],
               lorentz = [ L.VVVVS10, L.VVVVS2, L.VVVVS3, L.VVVVS5 ],
               loop_particles = [ [ [P.g] ], [ [P.t] ] ],
               couplings = {(0,1,0):C.R2GC_695_148,(2,1,0):C.R2GC_695_148,(6,1,0):C.R2GC_1220_79,(6,1,1):C.R2GC_1219_78,(7,1,0):C.R2GC_1219_77,(7,1,1):C.R2GC_1219_78,(5,1,0):C.R2GC_699_152,(1,1,0):C.R2GC_699_152,(4,1,0):C.R2GC_699_152,(3,1,0):C.R2GC_699_152,(8,1,0):C.R2GC_696_149,(11,0,0):C.R2GC_698_151,(10,0,0):C.R2GC_698_151,(9,0,0):C.R2GC_697_150,(0,2,0):C.R2GC_695_148,(2,2,0):C.R2GC_695_148,(6,2,0):C.R2GC_1222_82,(6,2,1):C.R2GC_1221_81,(8,2,0):C.R2GC_1219_77,(8,2,1):C.R2GC_1219_78,(5,2,0):C.R2GC_699_152,(1,2,0):C.R2GC_699_152,(7,2,0):C.R2GC_696_149,(4,2,0):C.R2GC_699_152,(3,2,0):C.R2GC_699_152,(0,3,0):C.R2GC_695_148,(2,3,0):C.R2GC_695_148,(7,3,0):C.R2GC_1221_80,(7,3,1):C.R2GC_1221_81,(8,3,0):C.R2GC_1221_80,(8,3,1):C.R2GC_1221_81,(5,3,0):C.R2GC_699_152,(1,3,0):C.R2GC_699_152,(4,3,0):C.R2GC_699_152,(3,3,0):C.R2GC_699_152})

V_9 = CTVertex(name = 'V_9',
               type = 'R2',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS16, L.FFS21, L.FFS25, L.FFS3 ],
               loop_particles = [ [ [P.g, P.t] ] ],
               couplings = {(0,3,0):[ C.R2GC_1278_107, C.R2GC_1310_137 ],(0,0,0):C.R2GC_1150_31,(0,1,0):C.R2GC_1145_26,(0,2,0):C.R2GC_787_198})

V_10 = CTVertex(name = 'V_10',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS14, L.FFS15, L.FFS24, L.FFS6, L.FFS7, L.FFS8 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,2):C.R2GC_1291_119,(0,0,1):[ C.R2GC_1291_120, C.R2GC_1279_108 ],(0,3,2):C.R2GC_1196_59,(0,5,1):C.R2GC_856_248,(0,2,0):C.R2GC_712_165,(0,1,2):C.R2GC_1132_17,(0,4,1):C.R2GC_1281_110})

V_11 = CTVertex(name = 'V_11',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1, L.FFS13, L.FFS18, L.FFS19, L.FFS20, L.FFS7 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,1,0):C.R2GC_1149_30,(0,5,0):C.R2GC_1282_111,(0,3,0):C.R2GC_1283_112,(0,0,0):[ C.R2GC_1277_106, C.R2GC_1293_123 ],(0,2,0):C.R2GC_786_197,(0,4,0):C.R2GC_1146_27})

V_12 = CTVertex(name = 'V_12',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS10, L.FFS12, L.FFS4, L.FFS5, L.FFS7, L.FFS9 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,2,2):C.R2GC_1292_121,(0,2,1):[ C.R2GC_1292_122, C.R2GC_1280_109 ],(0,5,2):C.R2GC_1197_60,(0,0,1):C.R2GC_857_249,(0,1,2):C.R2GC_785_196,(0,3,0):C.R2GC_1116_11,(0,4,1):C.R2GC_1281_110})

V_13 = CTVertex(name = 'V_13',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):[ C.R2GC_711_164, C.R2GC_939_403 ]})

V_14 = CTVertex(name = 'V_14',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS19, L.FFS7 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,1,0):C.R2GC_936_400,(0,0,0):C.R2GC_721_174})

V_15 = CTVertex(name = 'V_15',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.a, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS27, L.FFVS3, L.FFVS43 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,1,2):C.R2GC_1285_113,(0,1,1):C.R2GC_1285_114,(0,2,0):C.R2GC_713_166,(0,0,2):C.R2GC_789_200})

V_16 = CTVertex(name = 'V_16',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):[ C.R2GC_1228_86, C.R2GC_1316_143 ]})

V_17 = CTVertex(name = 'V_17',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS19, L.FFS7 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,1,0):C.R2GC_941_405,(0,0,0):C.R2GC_736_184})

V_18 = CTVertex(name = 'V_18',
                type = 'R2',
                particles = [ P.s__tilde__, P.c, P.a, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1190_53})

V_19 = CTVertex(name = 'V_19',
                type = 'R2',
                particles = [ P.s__tilde__, P.c, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS7 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1189_52})

V_20 = CTVertex(name = 'V_20',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):[ C.R2GC_711_164, C.R2GC_939_403 ]})

V_21 = CTVertex(name = 'V_21',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS19, L.FFS7 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,1,0):C.R2GC_946_409,(0,0,0):C.R2GC_721_174})

V_22 = CTVertex(name = 'V_22',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.a, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1191_54})

V_23 = CTVertex(name = 'V_23',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS7 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1189_52})

V_24 = CTVertex(name = 'V_24',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.a, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1191_54})

V_25 = CTVertex(name = 'V_25',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS7 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1189_52})

V_26 = CTVertex(name = 'V_26',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):[ C.R2GC_711_164, C.R2GC_939_403 ]})

V_27 = CTVertex(name = 'V_27',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS19, L.FFS7 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,1,0):C.R2GC_946_409,(0,0,0):C.R2GC_721_174})

V_28 = CTVertex(name = 'V_28',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.a, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS13, L.FFVS3, L.FFVS8 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,1,2):C.R2GC_1286_115,(0,1,1):C.R2GC_1286_116,(0,2,0):C.R2GC_714_167,(0,0,2):C.R2GC_788_199})

V_29 = CTVertex(name = 'V_29',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV16, L.FFV27 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,1,0):C.R2GC_1166_38,(0,0,0):[ C.R2GC_1228_86, C.R2GC_1316_143 ],(0,3,0):C.R2GC_822_230,(0,2,0):C.R2GC_823_231})

V_30 = CTVertex(name = 'V_30',
                type = 'R2',
                particles = [ P.d__tilde__, P.u, P.a, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1190_53})

V_31 = CTVertex(name = 'V_31',
                type = 'R2',
                particles = [ P.d__tilde__, P.u, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS7 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1189_52})

V_32 = CTVertex(name = 'V_32',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):[ C.R2GC_1228_86, C.R2GC_1316_143 ]})

V_33 = CTVertex(name = 'V_33',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS19, L.FFS7 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,1,0):C.R2GC_941_405,(0,0,0):C.R2GC_736_184})

V_34 = CTVertex(name = 'V_34',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_1115_10})

V_35 = CTVertex(name = 'V_35',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_1115_10})

V_36 = CTVertex(name = 'V_36',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_1115_10})

V_37 = CTVertex(name = 'V_37',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1115_10})

V_38 = CTVertex(name = 'V_38',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV19 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,1,0):C.R2GC_1294_124,(0,0,0):C.R2GC_1115_10,(0,2,0):C.R2GC_825_233})

V_39 = CTVertex(name = 'V_39',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1115_10})

V_40 = CTVertex(name = 'V_40',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.W__minus__, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ] ],
                couplings = {(0,1,0):C.R2GC_1297_127,(0,1,1):C.R2GC_885_277,(0,0,0):C.R2GC_725_178})

V_41 = CTVertex(name = 'V_41',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.W__minus__, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,1,0):C.R2GC_739_187,(0,0,0):C.R2GC_741_189})

V_42 = CTVertex(name = 'V_42',
                type = 'R2',
                particles = [ P.s__tilde__, P.c, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):[ C.R2GC_1188_51, C.R2GC_1195_58 ]})

V_43 = CTVertex(name = 'V_43',
                type = 'R2',
                particles = [ P.s__tilde__, P.c, P.W__minus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1192_55})

V_44 = CTVertex(name = 'V_44',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.W__minus__, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,1,0):C.R2GC_739_187,(0,0,0):C.R2GC_725_178})

V_45 = CTVertex(name = 'V_45',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.W__minus__, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,1,0):C.R2GC_739_187,(0,0,0):C.R2GC_725_178})

V_46 = CTVertex(name = 'V_46',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):[ C.R2GC_1188_51, C.R2GC_1312_139 ],(0,1,0):C.R2GC_892_284})

V_47 = CTVertex(name = 'V_47',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.W__minus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3, L.FFVS8 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_1296_126,(0,1,0):C.R2GC_867_259})

V_48 = CTVertex(name = 'V_48',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.W__minus__, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS16, L.FFVS3, L.FFVS43 ],
                loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,2,1):C.R2GC_1297_127,(0,0,1):C.R2GC_1299_129,(0,1,0):C.R2GC_885_277,(0,3,0):C.R2GC_871_263})

V_49 = CTVertex(name = 'V_49',
                type = 'R2',
                particles = [ P.d__tilde__, P.u, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):[ C.R2GC_1188_51, C.R2GC_1195_58 ]})

V_50 = CTVertex(name = 'V_50',
                type = 'R2',
                particles = [ P.d__tilde__, P.u, P.W__minus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1192_55})

V_51 = CTVertex(name = 'V_51',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.W__minus__, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,1,0):C.R2GC_739_187,(0,0,0):C.R2GC_741_189})

V_52 = CTVertex(name = 'V_52',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ] ],
                couplings = {(0,1,0):C.R2GC_1298_128,(0,1,1):C.R2GC_886_278,(0,0,0):C.R2GC_724_177})

V_53 = CTVertex(name = 'V_53',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV20 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):[ C.R2GC_1188_51, C.R2GC_1312_139 ],(0,1,0):C.R2GC_893_285})

V_54 = CTVertex(name = 'V_54',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.W__plus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3, L.FFVS43 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_1296_126,(0,1,0):C.R2GC_868_260})

V_55 = CTVertex(name = 'V_55',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,1,0):C.R2GC_738_186,(0,0,0):C.R2GC_740_188})

V_56 = CTVertex(name = 'V_56',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,1,0):C.R2GC_738_186,(0,0,0):C.R2GC_724_177})

V_57 = CTVertex(name = 'V_57',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):[ C.R2GC_1188_51, C.R2GC_1195_58 ]})

V_58 = CTVertex(name = 'V_58',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.W__plus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1192_55})

V_59 = CTVertex(name = 'V_59',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):[ C.R2GC_1188_51, C.R2GC_1195_58 ]})

V_60 = CTVertex(name = 'V_60',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.W__plus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1192_55})

V_61 = CTVertex(name = 'V_61',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,1,0):C.R2GC_738_186,(0,0,0):C.R2GC_724_177})

V_62 = CTVertex(name = 'V_62',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS16, L.FFVS3, L.FFVS4, L.FFVS7 ],
                loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,2,1):C.R2GC_1298_128,(0,0,1):C.R2GC_1300_130,(0,1,0):C.R2GC_886_278,(0,4,0):C.R2GC_870_262,(0,3,0):C.R2GC_872_264})

V_63 = CTVertex(name = 'V_63',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS3 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,1,0):C.R2GC_738_186,(0,0,0):C.R2GC_740_188})

V_64 = CTVertex(name = 'V_64',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV13, L.FFV2, L.FFV6 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,3,0):[ C.R2GC_719_172, C.R2GC_938_402 ],(0,1,0):[ C.R2GC_720_173, C.R2GC_940_404 ],(0,4,0):C.R2GC_729_182,(0,0,0):C.R2GC_728_181,(0,2,0):C.R2GC_730_183})

V_65 = CTVertex(name = 'V_65',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.Z, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS18, L.FFVS3 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,2,0):C.R2GC_937_401,(0,0,0):C.R2GC_726_179,(0,1,0):C.R2GC_727_180})

V_66 = CTVertex(name = 'V_66',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.Z, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS27, L.FFVS3, L.FFVS43 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,1,2):C.R2GC_1308_133,(0,1,1):C.R2GC_1308_134,(0,2,0):C.R2GC_932_396,(0,0,2):C.R2GC_812_221})

V_67 = CTVertex(name = 'V_67',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV11, L.FFV2, L.FFV6 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,3,0):[ C.R2GC_1259_104, C.R2GC_943_407 ],(0,1,0):[ C.R2GC_1268_105, C.R2GC_945_408 ],(0,4,0):C.R2GC_746_192,(0,0,0):C.R2GC_1311_138,(0,2,0):C.R2GC_1313_140})

V_68 = CTVertex(name = 'V_68',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.Z, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS18, L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,2,0):C.R2GC_942_406,(0,0,0):C.R2GC_742_190,(0,1,0):C.R2GC_743_191})

V_69 = CTVertex(name = 'V_69',
                type = 'R2',
                particles = [ P.s__tilde__, P.c, P.Z, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1193_56})

V_70 = CTVertex(name = 'V_70',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV13, L.FFV2, L.FFV6 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,3,0):[ C.R2GC_719_172, C.R2GC_948_411 ],(0,1,0):[ C.R2GC_720_173, C.R2GC_940_404 ],(0,4,0):C.R2GC_729_182,(0,0,0):C.R2GC_728_181,(0,2,0):C.R2GC_730_183})

V_71 = CTVertex(name = 'V_71',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.Z, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS18, L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,2,0):C.R2GC_947_410,(0,0,0):C.R2GC_726_179,(0,1,0):C.R2GC_727_180})

V_72 = CTVertex(name = 'V_72',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.Z, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1194_57})

V_73 = CTVertex(name = 'V_73',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.Z, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_1194_57})

V_74 = CTVertex(name = 'V_74',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV13, L.FFV2, L.FFV6 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,3,0):[ C.R2GC_719_172, C.R2GC_948_411 ],(0,1,0):[ C.R2GC_720_173, C.R2GC_940_404 ],(0,4,0):C.R2GC_729_182,(0,0,0):C.R2GC_728_181,(0,2,0):C.R2GC_730_183})

V_75 = CTVertex(name = 'V_75',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.Z, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS18, L.FFVS3 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,2,0):C.R2GC_947_410,(0,0,0):C.R2GC_726_179,(0,1,0):C.R2GC_727_180})

V_76 = CTVertex(name = 'V_76',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.Z, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS13, L.FFVS3, L.FFVS8 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,1,2):C.R2GC_1309_135,(0,1,1):C.R2GC_1309_136,(0,2,0):C.R2GC_933_397,(0,0,2):C.R2GC_811_220})

V_77 = CTVertex(name = 'V_77',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV11, L.FFV16, L.FFV2, L.FFV21, L.FFV27, L.FFV6, L.FFV8 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,4,0):[ C.R2GC_1259_104, C.R2GC_1314_141 ],(0,1,0):[ C.R2GC_1268_105, C.R2GC_1317_144 ],(0,7,0):C.R2GC_1315_142,(0,0,0):C.R2GC_1311_138,(0,2,0):C.R2GC_1313_140,(0,8,0):C.R2GC_830_238,(0,5,0):C.R2GC_828_236,(0,6,0):C.R2GC_831_239,(0,3,0):C.R2GC_832_240})

V_78 = CTVertex(name = 'V_78',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.Z, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS17, L.FFVS18, L.FFVS21, L.FFVS3, L.FFVS46, L.FFVS47 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,4,0):C.R2GC_1302_131,(0,0,0):C.R2GC_1303_132,(0,2,0):C.R2GC_1155_32,(0,1,0):C.R2GC_810_219,(0,5,0):C.R2GC_806_215,(0,6,0):C.R2GC_814_223,(0,3,0):C.R2GC_816_225})

V_79 = CTVertex(name = 'V_79',
                type = 'R2',
                particles = [ P.d__tilde__, P.u, P.Z, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS3 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_1193_56})

V_80 = CTVertex(name = 'V_80',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV1, L.FFV10, L.FFV11, L.FFV2, L.FFV6 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,3,0):[ C.R2GC_1259_104, C.R2GC_943_407 ],(0,1,0):[ C.R2GC_1268_105, C.R2GC_945_408 ],(0,4,0):C.R2GC_746_192,(0,0,0):C.R2GC_1311_138,(0,2,0):C.R2GC_1313_140})

V_81 = CTVertex(name = 'V_81',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.Z, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS15, L.FFVS18, L.FFVS3 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,2,0):C.R2GC_942_406,(0,0,0):C.R2GC_742_190,(0,1,0):C.R2GC_743_191})

V_82 = CTVertex(name = 'V_82',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.W__minus__, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS8 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_866_258})

V_83 = CTVertex(name = 'V_83',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.a, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS20, L.FFVS32 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,1,0):C.R2GC_790_201,(0,0,0):C.R2GC_792_203})

V_84 = CTVertex(name = 'V_84',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.a, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS1, L.FFVS18, L.FFVS21, L.FFVS47 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_737_185,(0,1,0):C.R2GC_1147_28,(0,3,0):C.R2GC_791_202,(0,2,0):C.R2GC_793_204})

V_85 = CTVertex(name = 'V_85',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.Z, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVS20, L.FFVS32, L.FFVS33 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,2,0):C.R2GC_807_216,(0,1,0):C.R2GC_813_222,(0,0,0):C.R2GC_815_224})

V_86 = CTVertex(name = 'V_86',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.g, P.G__minus__ ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFVS11, L.FFVS3, L.FFVS6 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,1,2):C.R2GC_1198_61,(0,1,1):C.R2GC_1198_62,(0,2,0):C.R2GC_715_168,(0,0,2):C.R2GC_715_168})

V_87 = CTVertex(name = 'V_87',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g, P.G0 ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFVS36 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_799_208})

V_88 = CTVertex(name = 'V_88',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g, P.H ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFVS1, L.FFVS18, L.FFVS37 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_1121_12,(0,2,0):C.R2GC_798_207,(0,1,0):C.R2GC_1148_29})

V_89 = CTVertex(name = 'V_89',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.g, P.g, P.G__minus__ ],
                color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                lorentz = [ L.FFVVS19, L.FFVVS20, L.FFVVS21, L.FFVVS32 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g], [P.g, P.t] ], [ [P.g] ], [ [P.g, P.t] ] ],
                couplings = {(0,3,1):C.R2GC_962_429,(2,2,2):C.R2GC_1243_87,(2,2,0):C.R2GC_1245_91,(2,2,3):C.R2GC_1245_90,(1,2,2):C.R2GC_1247_92,(1,2,0):C.R2GC_1249_95,(1,2,3):C.R2GC_1249_96,(2,1,2):C.R2GC_1247_92,(2,1,0):C.R2GC_1249_96,(2,1,3):C.R2GC_1249_95,(1,1,2):C.R2GC_1243_87,(1,1,0):C.R2GC_1245_90,(1,1,3):C.R2GC_1245_91,(2,0,0):C.R2GC_1135_18,(2,0,3):C.R2GC_1135_19,(1,0,0):C.R2GC_1135_19,(1,0,3):C.R2GC_1135_18})

V_90 = CTVertex(name = 'V_90',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g, P.g, P.G0 ],
                color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                lorentz = [ L.FFVVS12, L.FFVVS2, L.FFVVS5, L.FFVVS8 ],
                loop_particles = [ [ [P.g] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,1):C.R2GC_804_213,(2,2,0):C.R2GC_1255_101,(2,2,1):C.R2GC_1255_102,(1,2,0):C.R2GC_1251_97,(1,2,1):C.R2GC_1251_98,(2,1,0):C.R2GC_1251_97,(2,1,1):C.R2GC_1251_98,(1,1,0):C.R2GC_1255_101,(1,1,1):C.R2GC_1255_102,(2,3,1):C.R2GC_1144_25,(1,3,1):C.R2GC_1144_25})

V_91 = CTVertex(name = 'V_91',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g, P.g, P.H ],
                color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                lorentz = [ L.FFVVS1, L.FFVVS16, L.FFVVS20, L.FFVVS21, L.FFVVS3, L.FFVVS41, L.FFVVS42, L.FFVVS9 ],
                loop_particles = [ [ [P.g] ], [ [P.g, P.t] ] ],
                couplings = {(0,1,1):C.R2GC_805_214,(2,4,0):C.R2GC_1254_100,(1,4,0):C.R2GC_1253_99,(2,0,1):C.R2GC_1142_23,(1,0,1):C.R2GC_1143_24,(2,7,1):C.R2GC_1141_22,(1,7,1):C.R2GC_1141_22,(2,3,1):C.R2GC_1143_24,(1,3,1):C.R2GC_1142_23,(2,2,0):C.R2GC_1253_99,(1,2,0):C.R2GC_1254_100,(2,6,1):C.R2GC_1143_24,(1,6,1):C.R2GC_1142_23,(2,5,0):C.R2GC_1253_99,(1,5,0):C.R2GC_1254_100})

V_92 = CTVertex(name = 'V_92',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g, P.g ],
                color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                lorentz = [ L.FFVV1, L.FFVV10, L.FFVV11, L.FFVV2, L.FFVV21, L.FFVV22, L.FFVV4, L.FFVV5 ],
                loop_particles = [ [ [P.g] ], [ [P.g, P.t] ] ],
                couplings = {(0,7,1):C.R2GC_827_235,(2,3,0):C.R2GC_1290_118,(1,3,0):C.R2GC_1289_117,(2,0,1):C.R2GC_1163_35,(1,0,1):C.R2GC_1164_36,(2,6,1):C.R2GC_1162_34,(1,6,1):C.R2GC_1162_34,(2,2,1):C.R2GC_1164_36,(1,2,1):C.R2GC_1163_35,(2,1,0):C.R2GC_1289_117,(1,1,0):C.R2GC_1290_118,(2,5,1):C.R2GC_1164_36,(1,5,1):C.R2GC_1163_35,(2,4,0):C.R2GC_1289_117,(1,4,0):C.R2GC_1290_118})

V_93 = CTVertex(name = 'V_93',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVVS28 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_876_268})

V_94 = CTVertex(name = 'V_94',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVVS28 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_875_267})

V_95 = CTVertex(name = 'V_95',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.a, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVV15 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_895_287})

V_96 = CTVertex(name = 'V_96',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.a, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVVS35 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_878_270})

V_97 = CTVertex(name = 'V_97',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.W__minus__, P.W__plus__, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVVS33 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_859_251})

V_98 = CTVertex(name = 'V_98',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVVS55 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_861_253})

V_99 = CTVertex(name = 'V_99',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__, P.H ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFVVS78 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_860_252})

V_100 = CTVertex(name = 'V_100',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV40 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_889_281})

V_101 = CTVertex(name = 'V_101',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS35, L.FFVVS38 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_853_245,(0,1,0):C.R2GC_862_254})

V_102 = CTVertex(name = 'V_102',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS35, L.FFVVS38 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_852_244,(0,1,0):C.R2GC_864_256})

V_103 = CTVertex(name = 'V_103',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV17, L.FFVV20 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_888_280,(0,1,0):C.R2GC_891_283})

V_104 = CTVertex(name = 'V_104',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.W__plus__, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS24, L.FFVVS28 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,1,0):C.R2GC_855_247,(0,0,0):C.R2GC_865_257})

V_105 = CTVertex(name = 'V_105',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS43 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_869_261})

V_106 = CTVertex(name = 'V_106',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.g, P.G__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS24, L.FFVS3, L.FFVS38 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,2):C.R2GC_1199_63,(0,1,1):C.R2GC_1199_64,(0,2,0):C.R2GC_716_169,(0,0,2):C.R2GC_716_169})

V_107 = CTVertex(name = 'V_107',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.g, P.g, P.G__plus__ ],
                 color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVS40, L.FFVVS41, L.FFVVS42, L.FFVVS72 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g], [P.g, P.t] ], [ [P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,3,1):C.R2GC_961_428,(2,2,2):C.R2GC_1247_92,(2,2,0):C.R2GC_1247_94,(2,2,3):C.R2GC_1247_93,(1,2,2):C.R2GC_1243_87,(1,2,0):C.R2GC_1243_88,(1,2,3):C.R2GC_1243_89,(2,1,2):C.R2GC_1243_87,(2,1,0):C.R2GC_1243_89,(2,1,3):C.R2GC_1243_88,(1,1,2):C.R2GC_1247_92,(1,1,0):C.R2GC_1247_93,(1,1,3):C.R2GC_1247_94,(2,0,0):C.R2GC_1139_20,(2,0,3):C.R2GC_1139_21,(1,0,0):C.R2GC_1139_21,(1,0,3):C.R2GC_1139_20})

V_108 = CTVertex(name = 'V_108',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.W__minus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS60 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_877_269})

V_109 = CTVertex(name = 'V_109',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS45 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_873_265})

V_110 = CTVertex(name = 'V_110',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS45 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_874_266})

V_111 = CTVertex(name = 'V_111',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.a, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV24 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_894_286})

V_112 = CTVertex(name = 'V_112',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.W__minus__, P.W__plus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS73 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_858_250})

V_113 = CTVertex(name = 'V_113',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.W__minus__, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS45, L.FFVVS67 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,1,0):C.R2GC_865_257,(0,0,0):C.R2GC_854_246})

V_114 = CTVertex(name = 'V_114',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS60, L.FFVVS76 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_850_242,(0,1,0):C.R2GC_862_254})

V_115 = CTVertex(name = 'V_115',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS60, L.FFVVS76 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_851_243,(0,1,0):C.R2GC_863_255})

V_116 = CTVertex(name = 'V_116',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV29, L.FFVV39 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_887_279,(0,1,0):C.R2GC_890_282})

V_117 = CTVertex(name = 'V_117',
                 type = 'R2',
                 particles = [ P.u__tilde__, P.u ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_710_163})

V_118 = CTVertex(name = 'V_118',
                 type = 'R2',
                 particles = [ P.c__tilde__, P.c ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_710_163})

V_119 = CTVertex(name = 'V_119',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1, L.FF2, L.FF3 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):[ C.R2GC_1257_103, C.R2GC_1295_125 ],(0,0,0):[ C.R2GC_710_163, C.R2GC_1165_37 ],(0,2,0):C.R2GC_1160_33})

V_120 = CTVertex(name = 'V_120',
                 type = 'R2',
                 particles = [ P.d__tilde__, P.d ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_710_163})

V_121 = CTVertex(name = 'V_121',
                 type = 'R2',
                 particles = [ P.s__tilde__, P.s ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.R2GC_710_163})

V_122 = CTVertex(name = 'V_122',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_710_163})

V_123 = CTVertex(name = 'V_123',
                 type = 'R2',
                 particles = [ P.g, P.g ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VV2, L.VV3, L.VV4 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ], [ [P.t] ] ],
                 couplings = {(0,2,1):C.R2GC_692_145,(0,0,2):C.R2GC_703_156,(0,1,0):C.R2GC_901_293})

V_124 = CTVertex(name = 'V_124',
                 type = 'R2',
                 particles = [ P.u__tilde__, P.u, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_722_175})

V_125 = CTVertex(name = 'V_125',
                 type = 'R2',
                 particles = [ P.c__tilde__, P.c, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_722_175})

V_126 = CTVertex(name = 'V_126',
                 type = 'R2',
                 particles = [ P.d__tilde__, P.d, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_722_175})

V_127 = CTVertex(name = 'V_127',
                 type = 'R2',
                 particles = [ P.s__tilde__, P.s, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.R2GC_722_175})

V_128 = CTVertex(name = 'V_128',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS2 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_722_175})

V_129 = CTVertex(name = 'V_129',
                 type = 'R2',
                 particles = [ P.g, P.g, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVV1 ],
                 loop_particles = [ [ [P.b] ], [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,1):C.R2GC_906_306,(0,0,2):C.R2GC_906_307,(0,0,0):C.R2GC_925_372,(0,0,3):C.R2GC_925_373,(0,0,4):C.R2GC_925_374,(0,0,5):C.R2GC_925_375})

V_130 = CTVertex(name = 'V_130',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVS1, L.VVS2 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_911_316,(0,0,1):C.R2GC_911_317,(0,0,2):C.R2GC_911_318,(0,0,3):C.R2GC_911_319,(0,1,3):C.R2GC_700_153})

V_131 = CTVertex(name = 'V_131',
                 type = 'R2',
                 particles = [ P.a, P.a, P.g, P.g ],
                 color = [ 'Identity(3,4)' ],
                 lorentz = [ L.VVVV11 ],
                 loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                 couplings = {(0,0,0):[ C.R2GC_902_294, C.R2GC_926_376 ],(0,0,1):[ C.R2GC_902_295, C.R2GC_926_377 ]})

V_132 = CTVertex(name = 'V_132',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.Z ],
                 color = [ 'Identity(2,3)' ],
                 lorentz = [ L.VVVV11 ],
                 loop_particles = [ [ [P.b] ], [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,1):C.R2GC_908_310,(0,0,2):C.R2GC_908_311,(0,0,0):C.R2GC_929_384,(0,0,3):C.R2GC_929_385,(0,0,4):C.R2GC_929_386,(0,0,5):C.R2GC_929_387})

V_133 = CTVertex(name = 'V_133',
                 type = 'R2',
                 particles = [ P.g, P.g, P.Z, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVV11 ],
                 loop_particles = [ [ [P.b] ], [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,1):C.R2GC_910_314,(0,0,2):C.R2GC_910_315,(0,0,0):C.R2GC_931_392,(0,0,3):C.R2GC_931_393,(0,0,4):C.R2GC_931_394,(0,0,5):C.R2GC_931_395})

V_134 = CTVertex(name = 'V_134',
                 type = 'R2',
                 particles = [ P.g, P.g, P.W__minus__, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVV11 ],
                 loop_particles = [ [ [P.b, P.t] ], [ [P.b, P.t], [P.c, P.s], [P.d, P.u] ], [ [P.c, P.s], [P.d, P.u] ] ],
                 couplings = {(0,0,1):C.R2GC_973_444,(0,0,0):C.R2GC_980_473,(0,0,2):C.R2GC_980_474})

V_135 = CTVertex(name = 'V_135',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.g ],
                 color = [ 'd(2,3,4)' ],
                 lorentz = [ L.VVVV11 ],
                 loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                 couplings = {(0,0,0):[ C.R2GC_903_296, C.R2GC_928_382 ],(0,0,1):[ C.R2GC_903_297, C.R2GC_928_383 ]})

V_136 = CTVertex(name = 'V_136',
                 type = 'R2',
                 particles = [ P.g, P.g, P.g, P.Z ],
                 color = [ 'd(1,2,3)', 'f(1,2,3)' ],
                 lorentz = [ L.VVVV1, L.VVVV11 ],
                 loop_particles = [ [ [P.b] ], [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(1,0,1):C.R2GC_907_308,(1,0,2):C.R2GC_907_309,(1,0,0):C.R2GC_927_378,(1,0,3):C.R2GC_927_379,(1,0,4):C.R2GC_927_380,(1,0,5):C.R2GC_927_381,(0,1,1):C.R2GC_909_312,(0,1,2):C.R2GC_909_313,(0,1,0):C.R2GC_930_388,(0,1,3):C.R2GC_930_389,(0,1,4):C.R2GC_930_390,(0,1,5):C.R2GC_930_391})

V_137 = CTVertex(name = 'V_137',
                 type = 'R2',
                 particles = [ P.g, P.g, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVS1, L.VVVS2 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,1,0):C.R2GC_918_344,(0,1,1):C.R2GC_918_345,(0,1,2):C.R2GC_918_346,(0,1,3):C.R2GC_918_347,(0,0,3):C.R2GC_704_157})

V_138 = CTVertex(name = 'V_138',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.G0 ],
                 color = [ 'Identity(2,3)' ],
                 lorentz = [ L.VVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_912_320,(0,0,1):C.R2GC_912_321,(0,0,2):C.R2GC_912_322,(0,0,3):C.R2GC_912_323})

V_139 = CTVertex(name = 'V_139',
                 type = 'R2',
                 particles = [ P.g, P.g, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_924_368,(0,0,1):C.R2GC_924_369,(0,0,2):C.R2GC_924_370,(0,0,3):C.R2GC_924_371})

V_140 = CTVertex(name = 'V_140',
                 type = 'R2',
                 particles = [ P.g, P.g, P.g, P.G0 ],
                 color = [ 'd(1,2,3)', 'f(1,2,3)' ],
                 lorentz = [ L.VVVS10, L.VVVS3, L.VVVS4 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(1,1,0):C.R2GC_914_328,(1,1,1):C.R2GC_914_329,(1,1,2):C.R2GC_914_330,(1,1,3):C.R2GC_914_331,(0,0,0):C.R2GC_913_324,(0,0,1):C.R2GC_913_325,(0,0,2):C.R2GC_913_326,(0,0,3):C.R2GC_913_327,(1,2,3):C.R2GC_702_155})

V_141 = CTVertex(name = 'V_141',
                 type = 'R2',
                 particles = [ P.g, P.g, P.W__minus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVS1, L.VVVS10, L.VVVS2 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c, P.s], [P.d, P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,2,0):C.R2GC_916_336,(0,2,3):C.R2GC_916_337,(0,2,4):C.R2GC_916_338,(0,2,5):C.R2GC_916_339,(0,1,1):C.R2GC_977_459,(0,1,2):C.R2GC_977_460,(0,0,1):C.R2GC_779_193})

V_142 = CTVertex(name = 'V_142',
                 type = 'R2',
                 particles = [ P.g, P.g, P.W__plus__, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVS1, L.VVVS10, L.VVVS2 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c, P.s], [P.d, P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,2,0):C.R2GC_917_340,(0,2,3):C.R2GC_917_341,(0,2,4):C.R2GC_917_342,(0,2,5):C.R2GC_917_343,(0,1,1):C.R2GC_977_459,(0,1,2):C.R2GC_977_460,(0,0,1):C.R2GC_780_194})

V_143 = CTVertex(name = 'V_143',
                 type = 'R2',
                 particles = [ P.g, P.g, P.G0, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS2, L.VVSS3 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_904_298,(0,0,1):C.R2GC_904_299,(0,0,2):C.R2GC_904_300,(0,0,3):C.R2GC_904_301,(0,1,3):C.R2GC_707_160})

V_144 = CTVertex(name = 'V_144',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.a ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV7 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_824_232})

V_145 = CTVertex(name = 'V_145',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV7 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_971_442})

V_146 = CTVertex(name = 'V_146',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.Z, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV7 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_972_443})

V_147 = CTVertex(name = 'V_147',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.g ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVV8 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_826_234})

V_148 = CTVertex(name = 'V_148',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.g, P.Z ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVV27, L.FFVV37 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_829_237,(0,1,0):C.R2GC_833_241})

V_149 = CTVertex(name = 'V_149',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.g, P.W__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVV36 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_896_288})

V_150 = CTVertex(name = 'V_150',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.g, P.W__minus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVV16 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_897_289})

V_151 = CTVertex(name = 'V_151',
                 type = 'R2',
                 particles = [ P.u__tilde__, P.u, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_737_185})

V_152 = CTVertex(name = 'V_152',
                 type = 'R2',
                 particles = [ P.c__tilde__, P.c, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_737_185})

V_153 = CTVertex(name = 'V_153',
                 type = 'R2',
                 particles = [ P.d__tilde__, P.d, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_723_176})

V_154 = CTVertex(name = 'V_154',
                 type = 'R2',
                 particles = [ P.s__tilde__, P.s, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.R2GC_723_176})

V_155 = CTVertex(name = 'V_155',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_723_176})

V_156 = CTVertex(name = 'V_156',
                 type = 'R2',
                 particles = [ P.u__tilde__, P.u, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_1121_12})

V_157 = CTVertex(name = 'V_157',
                 type = 'R2',
                 particles = [ P.c__tilde__, P.c, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_1121_12})

V_158 = CTVertex(name = 'V_158',
                 type = 'R2',
                 particles = [ P.d__tilde__, P.d, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_1121_12})

V_159 = CTVertex(name = 'V_159',
                 type = 'R2',
                 particles = [ P.s__tilde__, P.s, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.R2GC_1121_12})

V_160 = CTVertex(name = 'V_160',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.R2GC_1121_12})

V_161 = CTVertex(name = 'V_161',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS16 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_795_206})

V_162 = CTVertex(name = 'V_162',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS16 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_966_435})

V_163 = CTVertex(name = 'V_163',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.Z, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS16 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_970_441})

V_164 = CTVertex(name = 'V_164',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.g, P.H ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS18 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_803_212})

V_165 = CTVertex(name = 'V_165',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.g, P.Z, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS46, L.FFVVS48 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_809_218,(0,1,0):C.R2GC_819_228})

V_166 = CTVertex(name = 'V_166',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.g, P.W__minus__, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS34 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_880_272})

V_167 = CTVertex(name = 'V_167',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.g, P.W__plus__, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS69 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_881_273})

V_168 = CTVertex(name = 'V_168',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.a, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS12 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_794_205})

V_169 = CTVertex(name = 'V_169',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS12 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_965_434})

V_170 = CTVertex(name = 'V_170',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.Z, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS12 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_969_440})

V_171 = CTVertex(name = 'V_171',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.g, P.G0 ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS14 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_802_211})

V_172 = CTVertex(name = 'V_172',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.g, P.Z, P.G0 ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS47, L.FFVVS49 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_808_217,(0,1,0):C.R2GC_820_229})

V_173 = CTVertex(name = 'V_173',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.g, P.W__minus__, P.G0 ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS34 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_882_274})

V_174 = CTVertex(name = 'V_174',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.g, P.W__plus__, P.G0 ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS69 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_879_271})

V_175 = CTVertex(name = 'V_175',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.a, P.a, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS72 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_959_424,(0,0,1):C.R2GC_959_425})

V_176 = CTVertex(name = 'V_176',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.a, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS72 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_963_430,(0,0,1):C.R2GC_963_431})

V_177 = CTVertex(name = 'V_177',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.Z, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS72 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_967_436,(0,0,1):C.R2GC_967_437})

V_178 = CTVertex(name = 'V_178',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.a, P.g, P.G__plus__ ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS69, L.FFVVS74 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.R2GC_718_171,(0,0,1):C.R2GC_800_209})

V_179 = CTVertex(name = 'V_179',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.b, P.g, P.Z, P.G__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS69, L.FFVVS74 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_935_399,(0,1,1):C.R2GC_817_226})

V_180 = CTVertex(name = 'V_180',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.g, P.W__minus__, P.G__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS69 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_883_275})

V_181 = CTVertex(name = 'V_181',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.a, P.a, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS32 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_960_426,(0,0,1):C.R2GC_960_427})

V_182 = CTVertex(name = 'V_182',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.a, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS32 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_964_432,(0,0,1):C.R2GC_964_433})

V_183 = CTVertex(name = 'V_183',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.Z, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS32 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_968_438,(0,0,1):C.R2GC_968_439})

V_184 = CTVertex(name = 'V_184',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.a, P.g, P.G__minus__ ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS29, L.FFVVS34 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_717_170,(0,1,1):C.R2GC_801_210})

V_185 = CTVertex(name = 'V_185',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.t, P.g, P.Z, P.G__minus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS29, L.FFVVS34 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.R2GC_934_398,(0,0,1):C.R2GC_818_227})

V_186 = CTVertex(name = 'V_186',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.g, P.W__plus__, P.G__minus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS34 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_884_276})

V_187 = CTVertex(name = 'V_187',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.g, P.g, P.g ],
                 color = [ 'T(3,-1,1)*T(4,-2,-1)*T(5,2,-2)', 'T(3,-1,1)*T(4,2,-2)*T(5,-2,-1)', 'T(3,-2,-1)*T(4,-1,1)*T(5,2,-2)', 'T(3,-2,-1)*T(4,2,-2)*T(5,-1,1)', 'T(3,2,-2)*T(4,-1,1)*T(5,-2,-1)', 'T(3,2,-2)*T(4,-2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV1, L.FFVVV13, L.FFVVV14, L.FFVVV15, L.FFVVV16, L.FFVVV17, L.FFVVV18, L.FFVVV2, L.FFVVV3, L.FFVVV4, L.FFVVV5, L.FFVVV6, L.FFVVV8 ],
                 loop_particles = [ [ [P.b, P.G__plus__] ], [ [P.G0, P.t], [P.H, P.t] ] ],
                 couplings = {(5,9,1):C.R2GC_1180_44,(4,9,1):C.R2GC_1179_43,(3,9,1):C.R2GC_1179_43,(5,11,1):C.R2GC_1179_43,(4,11,1):C.R2GC_1180_44,(2,11,1):C.R2GC_1179_43,(5,8,1):C.R2GC_1179_43,(3,8,1):C.R2GC_1180_44,(1,8,1):C.R2GC_1179_43,(3,10,1):C.R2GC_1179_43,(1,10,1):C.R2GC_1180_44,(0,10,1):C.R2GC_1179_43,(4,0,1):C.R2GC_1179_43,(2,0,1):C.R2GC_1180_44,(0,0,1):C.R2GC_1179_43,(1,7,1):C.R2GC_1179_43,(2,7,1):C.R2GC_1179_43,(0,12,1):C.R2GC_1180_44,(5,4,0):C.R2GC_1125_16,(4,4,0):C.R2GC_1124_15,(3,4,0):C.R2GC_1124_15,(5,6,0):C.R2GC_1124_15,(4,6,0):C.R2GC_1125_16,(2,6,0):C.R2GC_1124_15,(5,3,0):C.R2GC_1124_15,(3,3,0):C.R2GC_1125_16,(1,3,0):C.R2GC_1124_15,(3,5,0):C.R2GC_1124_15,(1,5,0):C.R2GC_1125_16,(0,5,0):C.R2GC_1124_15,(4,1,0):C.R2GC_1124_15,(2,1,0):C.R2GC_1125_16,(0,1,0):C.R2GC_1124_15,(1,2,0):C.R2GC_1124_15,(2,2,0):C.R2GC_1124_15,(0,2,0):C.R2GC_1125_16,(0,2,1):C.R2GC_1180_44})

V_188 = CTVertex(name = 'V_188',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b, P.g, P.g, P.g ],
                 color = [ 'T(3,-1,1)*T(4,-2,-1)*T(5,2,-2)', 'T(3,-1,1)*T(4,2,-2)*T(5,-2,-1)', 'T(3,-2,-1)*T(4,-1,1)*T(5,2,-2)', 'T(3,-2,-1)*T(4,2,-2)*T(5,-1,1)', 'T(3,2,-2)*T(4,-1,1)*T(5,-2,-1)', 'T(3,2,-2)*T(4,-2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV10, L.FFVVV11, L.FFVVV12, L.FFVVV7, L.FFVVV8, L.FFVVV9 ],
                 loop_particles = [ [ [P.G__plus__, P.t] ] ],
                 couplings = {(5,0,0):C.R2GC_1125_16,(4,0,0):C.R2GC_1124_15,(3,0,0):C.R2GC_1124_15,(5,2,0):C.R2GC_1124_15,(4,2,0):C.R2GC_1125_16,(2,2,0):C.R2GC_1124_15,(5,5,0):C.R2GC_1124_15,(3,5,0):C.R2GC_1125_16,(1,5,0):C.R2GC_1124_15,(3,1,0):C.R2GC_1124_15,(1,1,0):C.R2GC_1125_16,(0,1,0):C.R2GC_1124_15,(4,3,0):C.R2GC_1124_15,(2,3,0):C.R2GC_1125_16,(0,3,0):C.R2GC_1124_15,(1,4,0):C.R2GC_1124_15,(2,4,0):C.R2GC_1124_15,(0,4,0):C.R2GC_1125_16})

V_189 = CTVertex(name = 'V_189',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.a, P.g, P.g ],
                 color = [ 'T(4,-1,1)*T(5,2,-1)', 'T(4,2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV12, L.FFVVV14, L.FFVVV16, L.FFVVV17, L.FFVVV18, L.FFVVV4, L.FFVVV5, L.FFVVV8 ],
                 loop_particles = [ [ [P.b, P.G__plus__] ], [ [P.G0, P.t], [P.H, P.t] ] ],
                 couplings = {(1,5,1):C.R2GC_1123_14,(0,5,1):C.R2GC_1122_13,(1,6,1):C.R2GC_1123_14,(0,6,1):C.R2GC_1122_13,(1,0,1):C.R2GC_1122_13,(0,0,1):C.R2GC_1123_14,(1,7,1):C.R2GC_1122_13,(0,7,1):C.R2GC_1123_14,(1,2,0):C.R2GC_1122_13,(0,2,0):C.R2GC_1123_14,(1,4,0):C.R2GC_1123_14,(1,4,1):C.R2GC_1122_13,(0,4,0):C.R2GC_1122_13,(0,4,1):C.R2GC_1123_14,(1,3,0):C.R2GC_1122_13,(0,3,0):C.R2GC_1123_14,(1,1,0):C.R2GC_1123_14,(1,1,1):C.R2GC_1122_13,(0,1,0):C.R2GC_1122_13,(0,1,1):C.R2GC_1123_14})

V_190 = CTVertex(name = 'V_190',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b, P.a, P.g, P.g ],
                 color = [ 'T(4,-1,1)*T(5,2,-1)', 'T(4,2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV10, L.FFVVV11, L.FFVVV12, L.FFVVV8 ],
                 loop_particles = [ [ [P.G__plus__, P.t] ] ],
                 couplings = {(1,0,0):C.R2GC_1170_40,(0,0,0):C.R2GC_1169_39,(1,2,0):C.R2GC_1169_39,(0,2,0):C.R2GC_1170_40,(1,1,0):C.R2GC_1170_40,(0,1,0):C.R2GC_1169_39,(1,3,0):C.R2GC_1169_39,(0,3,0):C.R2GC_1170_40})

V_191 = CTVertex(name = 'V_191',
                 type = 'R2',
                 particles = [ P.t__tilde__, P.t, P.g, P.g, P.Z ],
                 color = [ 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVV1, L.FFVVV13, L.FFVVV14, L.FFVVV15, L.FFVVV16, L.FFVVV2, L.FFVVV3, L.FFVVV4 ],
                 loop_particles = [ [ [P.b, P.G__plus__] ], [ [P.G0, P.t], [P.H, P.t] ] ],
                 couplings = {(1,7,1):C.R2GC_1185_50,(0,7,1):C.R2GC_1184_49,(1,6,1):C.R2GC_1184_49,(0,6,1):C.R2GC_1185_50,(1,0,1):C.R2GC_1185_50,(0,0,1):C.R2GC_1184_49,(1,5,1):C.R2GC_1184_49,(0,5,1):C.R2GC_1185_50,(1,4,0):C.R2GC_1182_45,(1,4,1):C.R2GC_1182_46,(0,4,0):C.R2GC_1183_47,(0,4,1):C.R2GC_1183_48,(1,3,0):C.R2GC_1183_47,(1,3,1):C.R2GC_1183_48,(0,3,0):C.R2GC_1182_45,(0,3,1):C.R2GC_1182_46,(1,1,0):C.R2GC_1182_45,(1,1,1):C.R2GC_1182_46,(0,1,0):C.R2GC_1183_47,(0,1,1):C.R2GC_1183_48,(1,2,0):C.R2GC_1183_47,(1,2,1):C.R2GC_1183_48,(0,2,0):C.R2GC_1182_45,(0,2,1):C.R2GC_1182_46})

V_192 = CTVertex(name = 'V_192',
                 type = 'R2',
                 particles = [ P.b__tilde__, P.b, P.g, P.g, P.Z ],
                 color = [ 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVV10, L.FFVVV7, L.FFVVV8, L.FFVVV9 ],
                 loop_particles = [ [ [P.G__plus__, P.t] ] ],
                 couplings = {(1,0,0):C.R2GC_1174_42,(0,0,0):C.R2GC_1173_41,(1,3,0):C.R2GC_1173_41,(0,3,0):C.R2GC_1174_42,(1,1,0):C.R2GC_1174_42,(0,1,0):C.R2GC_1173_41,(1,2,0):C.R2GC_1173_41,(0,2,0):C.R2GC_1174_42})

V_193 = CTVertex(name = 'V_193',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.Z, P.H ],
                 color = [ 'Identity(2,3)' ],
                 lorentz = [ L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_919_348,(0,0,1):C.R2GC_919_349,(0,0,2):C.R2GC_919_350,(0,0,3):C.R2GC_919_351})

V_194 = CTVertex(name = 'V_194',
                 type = 'R2',
                 particles = [ P.g, P.g, P.Z, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_915_332,(0,0,1):C.R2GC_915_333,(0,0,2):C.R2GC_915_334,(0,0,3):C.R2GC_915_335})

V_195 = CTVertex(name = 'V_195',
                 type = 'R2',
                 particles = [ P.g, P.g, P.W__minus__, P.W__plus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVVS10 ],
                 loop_particles = [ [ [P.b, P.t] ], [ [P.c, P.s], [P.d, P.u] ] ],
                 couplings = {(0,0,0):C.R2GC_976_457,(0,0,1):C.R2GC_976_458})

V_196 = CTVertex(name = 'V_196',
                 type = 'R2',
                 particles = [ P.g, P.g, P.g, P.Z, P.H ],
                 color = [ 'd(1,2,3)', 'f(1,2,3)' ],
                 lorentz = [ L.VVVVS1, L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(1,0,0):C.R2GC_923_364,(1,0,1):C.R2GC_923_365,(1,0,2):C.R2GC_923_366,(1,0,3):C.R2GC_923_367,(0,1,0):C.R2GC_922_360,(0,1,1):C.R2GC_922_361,(0,1,2):C.R2GC_922_362,(0,1,3):C.R2GC_922_363})

V_197 = CTVertex(name = 'V_197',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.W__minus__, P.G__plus__ ],
                 color = [ 'Identity(2,3)' ],
                 lorentz = [ L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c, P.s], [P.d, P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_978_461,(0,0,3):C.R2GC_978_462,(0,0,4):C.R2GC_978_463,(0,0,5):C.R2GC_978_464,(0,0,1):C.R2GC_978_465,(0,0,2):C.R2GC_978_466})

V_198 = CTVertex(name = 'V_198',
                 type = 'R2',
                 particles = [ P.g, P.g, P.W__minus__, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c, P.s], [P.d, P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_974_445,(0,0,3):C.R2GC_974_446,(0,0,4):C.R2GC_974_447,(0,0,5):C.R2GC_974_448,(0,0,1):C.R2GC_974_449,(0,0,2):C.R2GC_974_450})

V_199 = CTVertex(name = 'V_199',
                 type = 'R2',
                 particles = [ P.g, P.g, P.g, P.W__minus__, P.G__plus__ ],
                 color = [ 'd(1,2,3)', 'f(1,2,3)' ],
                 lorentz = [ L.VVVVS1, L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(1,0,0):C.R2GC_957_414,(1,0,2):C.R2GC_957_415,(1,0,3):C.R2GC_957_416,(1,0,4):C.R2GC_957_417,(1,0,1):C.R2GC_957_418,(0,1,0):C.R2GC_920_352,(0,1,2):C.R2GC_920_353,(0,1,3):C.R2GC_920_354,(0,1,4):C.R2GC_920_355})

V_200 = CTVertex(name = 'V_200',
                 type = 'R2',
                 particles = [ P.a, P.g, P.g, P.W__plus__, P.G__minus__ ],
                 color = [ 'Identity(2,3)' ],
                 lorentz = [ L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c, P.s], [P.d, P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_979_467,(0,0,3):C.R2GC_979_468,(0,0,4):C.R2GC_979_469,(0,0,5):C.R2GC_979_470,(0,0,1):C.R2GC_979_471,(0,0,2):C.R2GC_979_472})

V_201 = CTVertex(name = 'V_201',
                 type = 'R2',
                 particles = [ P.g, P.g, P.W__plus__, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c, P.s], [P.d, P.u] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.R2GC_975_451,(0,0,3):C.R2GC_975_452,(0,0,4):C.R2GC_975_453,(0,0,5):C.R2GC_975_454,(0,0,1):C.R2GC_975_455,(0,0,2):C.R2GC_975_456})

V_202 = CTVertex(name = 'V_202',
                 type = 'R2',
                 particles = [ P.g, P.g, P.g, P.W__plus__, P.G__minus__ ],
                 color = [ 'd(1,2,3)', 'f(1,2,3)' ],
                 lorentz = [ L.VVVVS1, L.VVVVS10 ],
                 loop_particles = [ [ [P.b] ], [ [P.b, P.t] ], [ [P.c], [P.u] ], [ [P.d], [P.s] ], [ [P.t] ] ],
                 couplings = {(1,0,0):C.R2GC_958_419,(1,0,2):C.R2GC_958_420,(1,0,3):C.R2GC_958_421,(1,0,4):C.R2GC_958_422,(1,0,1):C.R2GC_958_423,(0,1,0):C.R2GC_921_356,(0,1,2):C.R2GC_921_357,(0,1,3):C.R2GC_921_358,(0,1,4):C.R2GC_921_359})

V_203 = CTVertex(name = 'V_203',
                 type = 'UV',
                 particles = [ P.g, P.g, P.G0, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS5, L.VVSS8, L.VVSS9 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1216_226,(0,0,2):C.UVGC_1216_227,(0,0,3):C.UVGC_1216_228,(0,2,1):C.UVGC_988_425,(0,1,3):C.UVGC_1005_4})

V_204 = CTVertex(name = 'V_204',
                 type = 'UV',
                 particles = [ P.g, P.g, P.G__minus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS5, L.VVSS8, L.VVSS9 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.b, P.t] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1216_226,(0,0,3):C.UVGC_1216_227,(0,0,4):C.UVGC_1216_228,(0,2,2):C.UVGC_988_425,(0,1,1):C.UVGC_1005_4})

V_205 = CTVertex(name = 'V_205',
                 type = 'UV',
                 particles = [ P.g, P.g, P.H, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVSS5, L.VVSS8, L.VVSS9 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1216_226,(0,0,2):C.UVGC_1216_227,(0,0,3):C.UVGC_1216_228,(0,2,1):C.UVGC_988_425,(0,1,3):C.UVGC_1005_4})

V_206 = CTVertex(name = 'V_206',
                 type = 'UV',
                 particles = [ P.g, P.g, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VVS4, L.VVS5, L.VVS7 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1217_229,(0,0,2):C.UVGC_1217_230,(0,0,3):C.UVGC_1217_231,(0,2,1):C.UVGC_991_428,(0,1,3):C.UVGC_992_429})

V_207 = CTVertex(name = 'V_207',
                 type = 'UV',
                 particles = [ P.g, P.g, P.g ],
                 color = [ 'f(1,2,3)' ],
                 lorentz = [ L.VVV3 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1211_210,(0,0,1):C.UVGC_1211_211,(0,0,2):C.UVGC_1211_212,(0,0,3):[ C.UVGC_1211_213, C.UVGC_1224_245 ]})

V_208 = CTVertex(name = 'V_208',
                 type = 'UV',
                 particles = [ P.g, P.g, P.g, P.H ],
                 color = [ 'f(1,2,3)' ],
                 lorentz = [ L.VVVS9 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1218_232,(0,0,1):C.UVGC_1218_233,(0,0,2):C.UVGC_1218_234,(0,0,3):C.UVGC_1218_235})

V_209 = CTVertex(name = 'V_209',
                 type = 'UV',
                 particles = [ P.g, P.g, P.g, P.g ],
                 color = [ 'd(-1,1,3)*d(-1,2,4)', 'd(-1,1,3)*f(-1,2,4)', 'd(-1,1,4)*d(-1,2,3)', 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,3)*f(-1,1,4)', 'd(-1,2,4)*f(-1,1,3)', 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)', 'Identity(1,2)*Identity(3,4)', 'Identity(1,3)*Identity(2,4)', 'Identity(1,4)*Identity(2,3)' ],
                 lorentz = [ L.VVVV11, L.VVVV2, L.VVVV3, L.VVVV5 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(2,1,1):C.UVGC_1110_111,(2,1,2):C.UVGC_1110_110,(0,1,1):C.UVGC_1110_111,(0,1,2):C.UVGC_1110_110,(4,1,1):C.UVGC_1109_108,(4,1,2):C.UVGC_1109_109,(3,1,1):C.UVGC_1109_108,(3,1,2):C.UVGC_1109_109,(8,1,1):C.UVGC_1110_110,(8,1,2):C.UVGC_1110_111,(6,1,0):C.UVGC_1214_220,(6,1,1):C.UVGC_1215_224,(6,1,2):C.UVGC_1215_225,(6,1,3):[ C.UVGC_1214_223, C.UVGC_1226_247 ],(7,1,0):C.UVGC_1214_220,(7,1,1):C.UVGC_1214_221,(7,1,2):C.UVGC_1214_222,(7,1,3):[ C.UVGC_1214_223, C.UVGC_1226_247 ],(5,1,1):C.UVGC_1109_108,(5,1,2):C.UVGC_1109_109,(1,1,1):C.UVGC_1109_108,(1,1,2):C.UVGC_1109_109,(11,0,1):C.UVGC_1113_114,(11,0,2):C.UVGC_1113_115,(10,0,1):C.UVGC_1113_114,(10,0,2):C.UVGC_1113_115,(9,0,1):C.UVGC_1112_112,(9,0,2):C.UVGC_1112_113,(2,2,1):C.UVGC_1110_111,(2,2,2):C.UVGC_1110_110,(0,2,1):C.UVGC_1110_111,(0,2,2):C.UVGC_1110_110,(4,2,1):C.UVGC_1109_108,(4,2,2):C.UVGC_1109_109,(3,2,1):C.UVGC_1109_108,(3,2,2):C.UVGC_1109_109,(8,2,0):C.UVGC_1214_220,(8,2,1):C.UVGC_1214_221,(8,2,2):C.UVGC_1214_222,(8,2,3):[ C.UVGC_1214_223, C.UVGC_1226_247 ],(6,2,0):C.UVGC_1212_214,(6,2,1):C.UVGC_1213_218,(6,2,2):C.UVGC_1213_219,(6,2,3):[ C.UVGC_1212_217, C.UVGC_1225_246 ],(5,2,1):C.UVGC_1109_108,(5,2,2):C.UVGC_1109_109,(1,2,1):C.UVGC_1109_108,(1,2,2):C.UVGC_1109_109,(7,2,1):C.UVGC_1110_110,(7,2,2):C.UVGC_1110_111,(0,3,1):C.UVGC_1110_111,(0,3,2):C.UVGC_1110_110,(2,3,1):C.UVGC_1110_111,(2,3,2):C.UVGC_1110_110,(5,3,1):C.UVGC_1109_108,(5,3,2):C.UVGC_1109_109,(1,3,1):C.UVGC_1109_108,(1,3,2):C.UVGC_1109_109,(7,3,0):C.UVGC_1212_214,(7,3,1):C.UVGC_1212_215,(7,3,2):C.UVGC_1212_216,(7,3,3):[ C.UVGC_1212_217, C.UVGC_1225_246 ],(4,3,1):C.UVGC_1109_108,(4,3,2):C.UVGC_1109_109,(3,3,1):C.UVGC_1109_108,(3,3,2):C.UVGC_1109_109,(8,3,0):C.UVGC_1212_214,(8,3,1):C.UVGC_1212_215,(8,3,2):C.UVGC_1212_216,(8,3,3):[ C.UVGC_1212_217, C.UVGC_1225_246 ]})

V_210 = CTVertex(name = 'V_210',
                 type = 'UV',
                 particles = [ P.g, P.g, P.g, P.g, P.H ],
                 color = [ 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)' ],
                 lorentz = [ L.VVVVS2, L.VVVVS3, L.VVVVS5 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1219_236,(0,0,1):C.UVGC_1219_237,(0,0,2):C.UVGC_1219_238,(0,0,3):C.UVGC_1219_239,(1,0,0):C.UVGC_1219_236,(1,0,1):C.UVGC_1219_237,(1,0,2):C.UVGC_1219_238,(1,0,3):C.UVGC_1219_239,(0,1,0):C.UVGC_1221_240,(0,1,1):C.UVGC_1221_241,(0,1,2):C.UVGC_1221_242,(0,1,3):C.UVGC_1221_243,(2,1,0):C.UVGC_1219_236,(2,1,1):C.UVGC_1219_237,(2,1,2):C.UVGC_1219_238,(2,1,3):C.UVGC_1219_239,(1,2,0):C.UVGC_1221_240,(1,2,1):C.UVGC_1221_241,(1,2,2):C.UVGC_1221_242,(1,2,3):C.UVGC_1221_243,(2,2,0):C.UVGC_1221_240,(2,2,1):C.UVGC_1221_241,(2,2,2):C.UVGC_1221_242,(2,2,3):C.UVGC_1221_243})

V_211 = CTVertex(name = 'V_211',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS16, L.FFS21, L.FFS23, L.FFS3 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,3,0):[ C.UVGC_1278_350, C.UVGC_1310_415 ],(0,0,0):C.UVGC_1150_153,(0,1,0):C.UVGC_1145_148,(0,2,0):C.UVGC_1009_8})

V_212 = CTVertex(name = 'V_212',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS14, L.FFS15, L.FFS22, L.FFS6, L.FFS7 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ], [ [P.g, P.t], [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):[ C.UVGC_1291_382, C.UVGC_1279_351 ],(0,0,2):[ C.UVGC_1291_383, C.UVGC_1279_352 ],(0,0,1):[ C.UVGC_1291_384, C.UVGC_1279_353 ],(0,3,3):C.UVGC_1196_203,(0,1,2):C.UVGC_1132_133,(0,2,0):C.UVGC_995_430,(0,4,0):C.UVGC_1281_357,(0,4,2):C.UVGC_1281_358,(0,4,1):C.UVGC_1281_359})

V_213 = CTVertex(name = 'V_213',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS1, L.FFS13, L.FFS17, L.FFS19, L.FFS20, L.FFS7 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1149_152,(0,5,0):C.UVGC_1282_360,(0,3,0):C.UVGC_1283_361,(0,0,0):[ C.UVGC_1277_349, C.UVGC_1293_388 ],(0,4,0):C.UVGC_1146_149,(0,2,0):C.UVGC_1008_7})

V_214 = CTVertex(name = 'V_214',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS11, L.FFS4, L.FFS5, L.FFS7, L.FFS9 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ], [ [P.g, P.t], [P.b, P.g, P.t] ] ],
                 couplings = {(0,1,0):[ C.UVGC_1292_385, C.UVGC_1280_354 ],(0,1,2):[ C.UVGC_1292_386, C.UVGC_1280_355 ],(0,1,1):[ C.UVGC_1292_387, C.UVGC_1280_356 ],(0,4,3):C.UVGC_1197_204,(0,2,0):C.UVGC_1116_123,(0,0,2):C.UVGC_1007_6,(0,3,0):C.UVGC_1281_357,(0,3,2):C.UVGC_1281_358,(0,3,1):C.UVGC_1281_359})

V_215 = CTVertex(name = 'V_215',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.a, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS28, L.FFVS3, L.FFVS30, L.FFVS44 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,2,0):C.UVGC_1273_343,(0,2,2):C.UVGC_1273_344,(0,1,0):C.UVGC_1285_364,(0,1,2):C.UVGC_1285_365,(0,1,1):C.UVGC_1285_366,(0,3,0):C.UVGC_996_431,(0,0,2):C.UVGC_1011_10})

V_216 = CTVertex(name = 'V_216',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.c, P.a, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1190_191,(0,0,1):C.UVGC_1190_192})

V_217 = CTVertex(name = 'V_217',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.c, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS7 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1189_189,(0,0,1):C.UVGC_1189_190})

V_218 = CTVertex(name = 'V_218',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.d, P.a, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1191_193,(0,0,1):C.UVGC_1191_194})

V_219 = CTVertex(name = 'V_219',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.d, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS7 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1189_189,(0,0,1):C.UVGC_1189_190})

V_220 = CTVertex(name = 'V_220',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.s, P.a, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1191_193,(0,0,1):C.UVGC_1191_194})

V_221 = CTVertex(name = 'V_221',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.s, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS7 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1189_189,(0,0,1):C.UVGC_1189_190})

V_222 = CTVertex(name = 'V_222',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.a, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS10, L.FFVS14, L.FFVS3, L.FFVS9 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1274_345,(0,0,2):C.UVGC_1274_346,(0,2,0):C.UVGC_1286_367,(0,2,2):C.UVGC_1286_368,(0,2,1):C.UVGC_1286_369,(0,3,0):C.UVGC_997_432,(0,1,2):C.UVGC_1010_9})

V_223 = CTVertex(name = 'V_223',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV1, L.FFV10, L.FFV14, L.FFV17, L.FFV24 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1166_172,(0,0,0):[ C.UVGC_1228_249, C.UVGC_1316_423 ],(0,3,0):C.UVGC_1050_49,(0,2,0):C.UVGC_1307_408,(0,4,0):C.UVGC_1051_50})

V_224 = CTVertex(name = 'V_224',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.u, P.a, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1190_191,(0,0,1):C.UVGC_1190_192})

V_225 = CTVertex(name = 'V_225',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.u, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFS7 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1189_189,(0,0,1):C.UVGC_1189_190})

V_226 = CTVertex(name = 'V_226',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.b, P.g] ], [ [P.g] ], [ [P.ghG] ] ],
                 couplings = {(0,0,0):C.UVGC_1115_119,(0,0,2):C.UVGC_1115_120,(0,0,3):C.UVGC_1115_121,(0,0,1):C.UVGC_1115_122})

V_227 = CTVertex(name = 'V_227',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.c, P.g] ], [ [P.g] ], [ [P.ghG] ] ],
                 couplings = {(0,0,0):C.UVGC_1115_119,(0,0,2):C.UVGC_1115_120,(0,0,3):C.UVGC_1115_121,(0,0,1):C.UVGC_1115_122})

V_228 = CTVertex(name = 'V_228',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.d, P.g] ], [ [P.g] ], [ [P.ghG] ] ],
                 couplings = {(0,0,0):C.UVGC_1115_119,(0,0,2):C.UVGC_1115_120,(0,0,3):C.UVGC_1115_121,(0,0,1):C.UVGC_1115_122})

V_229 = CTVertex(name = 'V_229',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1115_119,(0,0,1):C.UVGC_1115_120,(0,0,2):C.UVGC_1115_121,(0,0,3):C.UVGC_1115_122})

V_230 = CTVertex(name = 'V_230',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1, L.FFV10, L.FFV14, L.FFV15, L.FFV18, L.FFV26 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,3):C.UVGC_1294_389,(0,0,0):C.UVGC_1115_119,(0,0,1):C.UVGC_1115_120,(0,0,2):C.UVGC_1115_121,(0,0,3):C.UVGC_1229_250,(0,5,0):C.UVGC_1114_116,(0,5,1):C.UVGC_1114_117,(0,5,2):C.UVGC_1114_118,(0,4,3):C.UVGC_1056_55,(0,3,3):C.UVGC_1055_54,(0,2,3):C.UVGC_1208_207})

V_231 = CTVertex(name = 'V_231',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u, P.g ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFV1 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1115_119,(0,0,1):C.UVGC_1115_120,(0,0,2):C.UVGC_1115_121,(0,0,3):C.UVGC_1115_122})

V_232 = CTVertex(name = 'V_232',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.c, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):[ C.UVGC_1188_187, C.UVGC_1195_201 ],(0,0,1):[ C.UVGC_1188_188, C.UVGC_1195_202 ]})

V_233 = CTVertex(name = 'V_233',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.c, P.W__minus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1192_195,(0,0,1):C.UVGC_1192_196})

V_234 = CTVertex(name = 'V_234',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV4, L.FFV5 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,2,0):C.UVGC_1284_362,(0,2,2):C.UVGC_1284_363,(0,0,0):[ C.UVGC_1188_187, C.UVGC_1312_417 ],(0,0,2):[ C.UVGC_1258_321, C.UVGC_1312_418 ],(0,0,1):[ C.UVGC_1188_188, C.UVGC_1312_419 ],(0,1,1):C.UVGC_1103_101})

V_235 = CTVertex(name = 'V_235',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS10, L.FFVS3, L.FFVS9 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1231_253,(0,0,2):C.UVGC_1231_254,(0,1,0):C.UVGC_1296_391,(0,1,2):C.UVGC_1296_392,(0,1,1):C.UVGC_1296_393,(0,2,1):C.UVGC_1081_79})

V_236 = CTVertex(name = 'V_236',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.W__minus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS15, L.FFVS3, L.FFVS30, L.FFVS44 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,1):C.UVGC_1297_394,(0,0,1):C.UVGC_1299_396,(0,0,0):C.UVGC_1299_397,(0,3,0):C.UVGC_1084_82,(0,2,1):C.UVGC_1234_258})

V_237 = CTVertex(name = 'V_237',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.u, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):[ C.UVGC_1188_187, C.UVGC_1195_201 ],(0,0,1):[ C.UVGC_1188_188, C.UVGC_1195_202 ]})

V_238 = CTVertex(name = 'V_238',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.u, P.W__minus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1192_195,(0,0,1):C.UVGC_1192_196})

V_239 = CTVertex(name = 'V_239',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2, L.FFV23, L.FFV25 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,2,0):C.UVGC_1284_362,(0,2,2):C.UVGC_1284_363,(0,0,0):[ C.UVGC_1188_187, C.UVGC_1312_417 ],(0,0,2):[ C.UVGC_1258_321, C.UVGC_1312_418 ],(0,0,1):[ C.UVGC_1188_188, C.UVGC_1312_419 ],(0,1,1):C.UVGC_1103_101})

V_240 = CTVertex(name = 'V_240',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3, L.FFVS30, L.FFVS44 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1231_253,(0,1,2):C.UVGC_1231_254,(0,0,0):C.UVGC_1296_391,(0,0,2):C.UVGC_1296_392,(0,0,1):C.UVGC_1296_393,(0,2,1):C.UVGC_1081_79})

V_241 = CTVertex(name = 'V_241',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.d, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):[ C.UVGC_1188_187, C.UVGC_1195_201 ],(0,0,1):[ C.UVGC_1188_188, C.UVGC_1195_202 ]})

V_242 = CTVertex(name = 'V_242',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.d, P.W__plus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1192_195,(0,0,1):C.UVGC_1192_196})

V_243 = CTVertex(name = 'V_243',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.s, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV2 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):[ C.UVGC_1188_187, C.UVGC_1195_201 ],(0,0,1):[ C.UVGC_1188_188, C.UVGC_1195_202 ]})

V_244 = CTVertex(name = 'V_244',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.s, P.W__plus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1192_195,(0,0,1):C.UVGC_1192_196})

V_245 = CTVertex(name = 'V_245',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.W__plus__, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS10, L.FFVS15, L.FFVS3, L.FFVS9 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,2,1):C.UVGC_1298_395,(0,1,1):C.UVGC_1300_398,(0,1,0):C.UVGC_1300_399,(0,3,0):C.UVGC_1083_81,(0,0,1):C.UVGC_1233_257})

V_246 = CTVertex(name = 'V_246',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS28, L.FFVS3, L.FFVS30, L.FFVS44 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,2,0):C.UVGC_1270_339,(0,2,2):C.UVGC_1270_340,(0,1,0):C.UVGC_1308_409,(0,1,2):C.UVGC_1308_410,(0,1,1):C.UVGC_1308_411,(0,3,0):C.UVGC_1117_124,(0,0,2):C.UVGC_1041_40})

V_247 = CTVertex(name = 'V_247',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.c, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1193_197,(0,0,1):C.UVGC_1193_198})

V_248 = CTVertex(name = 'V_248',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.d, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1194_199,(0,0,1):C.UVGC_1194_200})

V_249 = CTVertex(name = 'V_249',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.s, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1194_199,(0,0,1):C.UVGC_1194_200})

V_250 = CTVertex(name = 'V_250',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS10, L.FFVS14, L.FFVS3, L.FFVS9 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1269_337,(0,0,2):C.UVGC_1269_338,(0,2,0):C.UVGC_1309_412,(0,2,2):C.UVGC_1309_413,(0,2,1):C.UVGC_1309_414,(0,3,0):C.UVGC_1118_125,(0,1,2):C.UVGC_1040_39})

V_251 = CTVertex(name = 'V_251',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFV1, L.FFV10, L.FFV11, L.FFV14, L.FFV17, L.FFV2, L.FFV22, L.FFV24, L.FFV6 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,5,0):[ C.UVGC_1259_322, C.UVGC_1314_421 ],(0,1,0):[ C.UVGC_1268_336, C.UVGC_1317_424 ],(0,8,0):C.UVGC_1315_422,(0,0,0):C.UVGC_1311_416,(0,2,0):C.UVGC_1313_420,(0,6,0):C.UVGC_1059_58,(0,4,0):C.UVGC_1061_60,(0,3,0):C.UVGC_1306_407,(0,7,0):C.UVGC_1062_61})

V_252 = CTVertex(name = 'V_252',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS15, L.FFVS18, L.FFVS19, L.FFVS29, L.FFVS3, L.FFVS42, L.FFVS45 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,4,0):C.UVGC_1302_401,(0,0,0):C.UVGC_1303_402,(0,1,0):C.UVGC_1155_160,(0,6,0):C.UVGC_1036_35,(0,5,0):C.UVGC_1043_42,(0,2,0):C.UVGC_1271_341,(0,3,0):C.UVGC_1045_44})

V_253 = CTVertex(name = 'V_253',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.u, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS3 ],
                 loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1193_197,(0,0,1):C.UVGC_1193_198})

V_254 = CTVertex(name = 'V_254',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS10, L.FFVS9 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1230_251,(0,0,2):C.UVGC_1230_252,(0,1,1):C.UVGC_1080_78})

V_255 = CTVertex(name = 'V_255',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS2, L.FFVS25, L.FFVS35 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1276_348,(0,2,0):C.UVGC_1012_11,(0,1,0):C.UVGC_1014_13})

V_256 = CTVertex(name = 'V_256',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS18, L.FFVS19, L.FFVS29, L.FFVS42 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1275_347,(0,0,0):C.UVGC_1147_150,(0,3,0):C.UVGC_1013_12,(0,2,0):C.UVGC_1015_14})

V_257 = CTVertex(name = 'V_257',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS2, L.FFVS25, L.FFVS34, L.FFVS35 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1272_342,(0,2,0):C.UVGC_1037_36,(0,3,0):C.UVGC_1042_41,(0,1,0):C.UVGC_1044_43})

V_258 = CTVertex(name = 'V_258',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.g, P.G__minus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS10, L.FFVS12, L.FFVS3, L.FFVS5 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.b, P.g] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ], [ [P.g, P.t], [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1239_266,(0,0,2):C.UVGC_1239_267,(0,0,3):C.UVGC_1239_268,(0,0,4):C.UVGC_1239_269,(0,2,5):C.UVGC_1198_205,(0,3,1):C.UVGC_1022_21,(0,1,4):C.UVGC_1022_21})

V_259 = CTVertex(name = 'V_259',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.G0 ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS2, L.FFVS26, L.FFVS39 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1242_277,(0,0,2):C.UVGC_1242_278,(0,0,3):C.UVGC_1242_279,(0,2,3):C.UVGC_1025_24,(0,1,1):C.UVGC_989_426})

V_260 = CTVertex(name = 'V_260',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1, L.FFVS18, L.FFVS19, L.FFVS23, L.FFVS40 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,3):C.UVGC_1121_128,(0,2,0):C.UVGC_1241_274,(0,2,2):C.UVGC_1241_275,(0,2,3):C.UVGC_1241_276,(0,1,3):C.UVGC_1148_151,(0,4,3):C.UVGC_1024_23,(0,3,1):C.UVGC_990_427})

V_261 = CTVertex(name = 'V_261',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.g, P.g, P.G__minus__ ],
                 color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVS19, L.FFVVS20, L.FFVVS21, L.FFVVS30 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.b, P.g] ], [ [P.b, P.g], [P.g, P.t] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(2,2,0):C.UVGC_1243_280,(2,2,3):C.UVGC_1243_281,(2,2,4):C.UVGC_1243_282,(2,2,1):C.UVGC_1246_289,(2,2,5):C.UVGC_1246_290,(1,2,0):C.UVGC_1247_291,(1,2,3):C.UVGC_1247_292,(1,2,4):C.UVGC_1247_293,(1,2,1):C.UVGC_1249_298,(1,2,5):C.UVGC_1249_299,(2,1,0):C.UVGC_1247_291,(2,1,3):C.UVGC_1247_292,(2,1,4):C.UVGC_1247_293,(2,1,1):C.UVGC_1250_300,(2,1,5):C.UVGC_1250_301,(1,1,0):C.UVGC_1243_280,(1,1,3):C.UVGC_1243_281,(1,1,4):C.UVGC_1243_282,(1,1,1):C.UVGC_1245_287,(1,1,5):C.UVGC_1245_288,(2,0,1):C.UVGC_1135_138,(2,0,5):C.UVGC_1135_139,(1,0,1):C.UVGC_1135_139,(1,0,5):C.UVGC_1135_138,(0,3,2):C.UVGC_1138_141})

V_262 = CTVertex(name = 'V_262',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.g, P.G0 ],
                 color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVS10, L.FFVVS2, L.FFVVS21, L.FFVVS42, L.FFVVS5, L.FFVVS8 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(2,4,0):C.UVGC_1255_315,(2,4,1):C.UVGC_1255_316,(2,4,2):C.UVGC_1255_317,(2,4,3):C.UVGC_1256_319,(1,4,0):C.UVGC_1251_302,(1,4,1):C.UVGC_1251_303,(1,4,2):C.UVGC_1251_304,(1,4,3):C.UVGC_1251_305,(2,1,0):C.UVGC_1251_302,(2,1,1):C.UVGC_1251_303,(2,1,2):C.UVGC_1251_304,(2,1,3):C.UVGC_1252_306,(1,1,0):C.UVGC_1255_315,(1,1,1):C.UVGC_1255_316,(1,1,2):C.UVGC_1255_317,(1,1,3):C.UVGC_1255_318,(2,5,3):C.UVGC_1144_147,(1,5,3):C.UVGC_1144_147,(0,0,3):C.UVGC_1032_31,(2,2,3):C.UVGC_1035_34,(1,2,3):C.UVGC_1031_30,(2,3,3):C.UVGC_1030_29,(1,3,3):C.UVGC_1034_33})

V_263 = CTVertex(name = 'V_263',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.g, P.H ],
                 color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVS1, L.FFVVS20, L.FFVVS21, L.FFVVS3, L.FFVVS41, L.FFVVS42, L.FFVVS58, L.FFVVS9 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(2,3,0):C.UVGC_1254_311,(2,3,1):C.UVGC_1254_312,(2,3,2):C.UVGC_1254_313,(2,3,3):C.UVGC_1254_314,(1,3,0):C.UVGC_1253_307,(1,3,1):C.UVGC_1253_308,(1,3,2):C.UVGC_1253_309,(1,3,3):C.UVGC_1253_310,(2,0,3):C.UVGC_1142_145,(1,0,3):C.UVGC_1143_146,(2,7,3):C.UVGC_1141_144,(1,7,3):C.UVGC_1141_144,(2,2,3):C.UVGC_1143_146,(1,2,3):C.UVGC_1142_145,(2,1,0):C.UVGC_1253_307,(2,1,1):C.UVGC_1253_308,(2,1,2):C.UVGC_1253_309,(2,1,3):C.UVGC_1253_310,(1,1,0):C.UVGC_1254_311,(1,1,1):C.UVGC_1254_312,(1,1,2):C.UVGC_1254_313,(1,1,3):C.UVGC_1254_314,(2,5,3):C.UVGC_1143_146,(1,5,3):C.UVGC_1142_145,(2,4,0):C.UVGC_1253_307,(2,4,1):C.UVGC_1253_308,(2,4,2):C.UVGC_1253_309,(2,4,3):C.UVGC_1253_310,(1,4,0):C.UVGC_1254_311,(1,4,1):C.UVGC_1254_312,(1,4,2):C.UVGC_1254_313,(1,4,3):C.UVGC_1254_314,(0,6,3):C.UVGC_1033_32})

V_264 = CTVertex(name = 'V_264',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.g ],
                 color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVV1, L.FFVV10, L.FFVV11, L.FFVV2, L.FFVV21, L.FFVV22, L.FFVV26, L.FFVV4 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(2,3,0):C.UVGC_1290_378,(2,3,1):C.UVGC_1290_379,(2,3,2):C.UVGC_1290_380,(2,3,3):C.UVGC_1290_381,(1,3,0):C.UVGC_1289_374,(1,3,1):C.UVGC_1289_375,(1,3,2):C.UVGC_1289_376,(1,3,3):C.UVGC_1289_377,(2,0,3):C.UVGC_1163_169,(1,0,3):C.UVGC_1164_170,(2,7,3):C.UVGC_1162_168,(1,7,3):C.UVGC_1162_168,(2,2,3):C.UVGC_1164_170,(1,2,3):C.UVGC_1163_169,(2,1,0):C.UVGC_1289_374,(2,1,1):C.UVGC_1289_375,(2,1,2):C.UVGC_1289_376,(2,1,3):C.UVGC_1289_377,(1,1,0):C.UVGC_1290_378,(1,1,1):C.UVGC_1290_379,(1,1,2):C.UVGC_1290_380,(1,1,3):C.UVGC_1290_381,(2,5,3):C.UVGC_1164_170,(1,5,3):C.UVGC_1163_169,(2,4,0):C.UVGC_1289_374,(2,4,1):C.UVGC_1289_375,(2,4,2):C.UVGC_1289_376,(2,4,3):C.UVGC_1289_377,(1,4,0):C.UVGC_1290_378,(1,4,1):C.UVGC_1290_379,(1,4,2):C.UVGC_1290_380,(1,4,3):C.UVGC_1290_381,(0,6,3):C.UVGC_1058_57})

V_265 = CTVertex(name = 'V_265',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS22, L.FFVVS36 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1235_259,(0,0,2):C.UVGC_1235_260,(0,1,1):C.UVGC_1088_86})

V_266 = CTVertex(name = 'V_266',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS22, L.FFVVS36 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1237_263,(0,0,2):C.UVGC_1237_264,(0,1,1):C.UVGC_1087_85})

V_267 = CTVertex(name = 'V_267',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.a, P.W__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV12, L.FFVV18 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1288_372,(0,0,2):C.UVGC_1288_373,(0,1,1):C.UVGC_1105_103})

V_268 = CTVertex(name = 'V_268',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.W__plus__, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS22, L.FFVVS26 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,1):C.UVGC_1238_265,(0,1,0):C.UVGC_1090_88})

V_269 = CTVertex(name = 'V_269',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.W__plus__, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS22, L.FFVVS27 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1261_325,(0,0,2):C.UVGC_1261_326,(0,1,1):C.UVGC_1073_71})

V_270 = CTVertex(name = 'V_270',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS50, L.FFVVS6 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,1):C.UVGC_1263_328,(0,0,0):C.UVGC_1075_73})

V_271 = CTVertex(name = 'V_271',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS56, L.FFVVS79 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,1):C.UVGC_1262_327,(0,1,0):C.UVGC_1074_72})

V_272 = CTVertex(name = 'V_272',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV25, L.FFVV41 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,1):C.UVGC_1301_400,(0,1,0):C.UVGC_1100_98})

V_273 = CTVertex(name = 'V_273',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS22, L.FFVVS23, L.FFVVS26 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1264_329,(0,0,2):C.UVGC_1264_330,(0,2,1):C.UVGC_1068_66,(0,1,1):C.UVGC_1076_74})

V_274 = CTVertex(name = 'V_274',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS22, L.FFVVS23, L.FFVVS26 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1266_333,(0,0,2):C.UVGC_1266_334,(0,2,1):C.UVGC_1066_64,(0,1,1):C.UVGC_1078_76})

V_275 = CTVertex(name = 'V_275',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV12, L.FFVV13, L.FFVV14 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1305_405,(0,0,2):C.UVGC_1305_406,(0,2,1):C.UVGC_1099_97,(0,1,1):C.UVGC_1102_100})

V_276 = CTVertex(name = 'V_276',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.W__plus__, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS22, L.FFVVS36, L.FFVVS39 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,1):C.UVGC_1267_335,(0,1,0):C.UVGC_1070_68,(0,2,0):C.UVGC_1079_77})

V_277 = CTVertex(name = 'V_277',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVS30, L.FFVS44 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1232_255,(0,0,2):C.UVGC_1232_256,(0,1,1):C.UVGC_1082_80})

V_278 = CTVertex(name = 'V_278',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.g, P.G__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS22, L.FFVS3, L.FFVS30, L.FFVS41 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.b, P.g] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ], [ [P.g, P.t], [P.b, P.g, P.t] ] ],
                 couplings = {(0,2,0):C.UVGC_1240_270,(0,2,2):C.UVGC_1240_271,(0,2,3):C.UVGC_1240_272,(0,2,4):C.UVGC_1240_273,(0,1,5):C.UVGC_1199_206,(0,3,1):C.UVGC_1023_22,(0,0,4):C.UVGC_1023_22})

V_279 = CTVertex(name = 'V_279',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.g, P.g, P.G__plus__ ],
                 color = [ 'Identity(1,2)*Identity(3,4)', 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVS40, L.FFVVS41, L.FFVVS42, L.FFVVS70 ],
                 loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.u] ], [ [P.b, P.g] ], [ [P.b, P.g], [P.g, P.t] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ] ],
                 couplings = {(2,2,0):C.UVGC_1247_291,(2,2,3):C.UVGC_1247_292,(2,2,4):C.UVGC_1247_293,(2,2,1):C.UVGC_1248_296,(2,2,5):C.UVGC_1248_297,(1,2,0):C.UVGC_1243_280,(1,2,3):C.UVGC_1243_281,(1,2,4):C.UVGC_1243_282,(1,2,1):C.UVGC_1243_283,(1,2,5):C.UVGC_1243_284,(2,1,0):C.UVGC_1243_280,(2,1,3):C.UVGC_1243_281,(2,1,4):C.UVGC_1243_282,(2,1,1):C.UVGC_1244_285,(2,1,5):C.UVGC_1244_286,(1,1,0):C.UVGC_1247_291,(1,1,3):C.UVGC_1247_292,(1,1,4):C.UVGC_1247_293,(1,1,1):C.UVGC_1247_294,(1,1,5):C.UVGC_1247_295,(2,0,1):C.UVGC_1139_142,(2,0,5):C.UVGC_1139_143,(1,0,1):C.UVGC_1139_143,(1,0,5):C.UVGC_1139_142,(0,3,2):C.UVGC_1137_140})

V_280 = CTVertex(name = 'V_280',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.W__minus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS41, L.FFVVS59, L.FFVVS62 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,2,1):C.UVGC_1238_265,(0,0,0):C.UVGC_1089_87,(0,1,0):C.UVGC_1091_89})

V_281 = CTVertex(name = 'V_281',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS61, L.FFVVS62 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1235_259,(0,1,2):C.UVGC_1235_260,(0,0,1):C.UVGC_1085_83})

V_282 = CTVertex(name = 'V_282',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS61, L.FFVVS62 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1236_261,(0,1,2):C.UVGC_1236_262,(0,0,1):C.UVGC_1086_84})

V_283 = CTVertex(name = 'V_283',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.a, P.W__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV30, L.FFVV31 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1287_370,(0,1,2):C.UVGC_1287_371,(0,0,1):C.UVGC_1104_102})

V_284 = CTVertex(name = 'V_284',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__minus__, P.W__plus__, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS62, L.FFVVS65 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1260_323,(0,0,2):C.UVGC_1260_324,(0,1,1):C.UVGC_1072_70})

V_285 = CTVertex(name = 'V_285',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.W__minus__, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS42, L.FFVVS44, L.FFVVS62, L.FFVVS77 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,2,1):C.UVGC_1267_335,(0,0,0):C.UVGC_1069_67,(0,3,0):C.UVGC_1079_77,(0,1,0):C.UVGC_1071_69})

V_286 = CTVertex(name = 'V_286',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS62, L.FFVVS64, L.FFVVS66 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1264_329,(0,0,2):C.UVGC_1264_330,(0,1,1):C.UVGC_1067_65,(0,2,1):C.UVGC_1076_74})

V_287 = CTVertex(name = 'V_287',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS62, L.FFVVS64, L.FFVVS66 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1265_331,(0,0,2):C.UVGC_1265_332,(0,1,1):C.UVGC_1065_63,(0,2,1):C.UVGC_1077_75})

V_288 = CTVertex(name = 'V_288',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV31, L.FFVV34, L.FFVV35 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1304_403,(0,0,2):C.UVGC_1304_404,(0,1,1):C.UVGC_1098_96,(0,2,1):C.UVGC_1101_99})

V_289 = CTVertex(name = 'V_289',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF5 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1002_3})

V_290 = CTVertex(name = 'V_290',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF5 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_1002_3})

V_291 = CTVertex(name = 'V_291',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF1, L.FF2, L.FF3, L.FF4, L.FF6 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):[ C.UVGC_1257_320, C.UVGC_1295_390 ],(0,0,0):C.UVGC_1165_171,(0,3,0):[ C.UVGC_1227_248, C.UVGC_1209_208 ],(0,4,0):C.UVGC_1006_5,(0,2,0):C.UVGC_1160_167})

V_292 = CTVertex(name = 'V_292',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF5 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_1002_3})

V_293 = CTVertex(name = 'V_293',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF5 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1002_3})

V_294 = CTVertex(name = 'V_294',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FF5 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_1002_3})

V_295 = CTVertex(name = 'V_295',
                 type = 'UV',
                 particles = [ P.g, P.g ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VV1, L.VV5 ],
                 loop_particles = [ [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1108_106,(0,0,1):C.UVGC_1108_107,(0,1,2):[ C.UVGC_1210_209, C.UVGC_1223_244 ]})

V_296 = CTVertex(name = 'V_296',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.a ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV6 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1052_51})

V_297 = CTVertex(name = 'V_297',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV6 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1167_173})

V_298 = CTVertex(name = 'V_298',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.Z, P.Z ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVV3, L.FFVV4, L.FFVV6 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1053_52,(0,2,0):C.UVGC_1168_174,(0,1,0):C.UVGC_1054_53})

V_299 = CTVertex(name = 'V_299',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.g ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVV9 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1057_56})

V_300 = CTVertex(name = 'V_300',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.Z ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVV28, L.FFVV38 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1060_59,(0,1,0):C.UVGC_1063_62})

V_301 = CTVertex(name = 'V_301',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.g, P.W__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVV33 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1106_104})

V_302 = CTVertex(name = 'V_302',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.g, P.W__minus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVV19 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1107_105})

V_303 = CTVertex(name = 'V_303',
                 type = 'UV',
                 particles = [ P.u__tilde__, P.u, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.g, P.u] ] ],
                 couplings = {(0,0,0):C.UVGC_1121_128})

V_304 = CTVertex(name = 'V_304',
                 type = 'UV',
                 particles = [ P.c__tilde__, P.c, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.c, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_1121_128})

V_305 = CTVertex(name = 'V_305',
                 type = 'UV',
                 particles = [ P.d__tilde__, P.d, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.d, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_1121_128})

V_306 = CTVertex(name = 'V_306',
                 type = 'UV',
                 particles = [ P.s__tilde__, P.s, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.g, P.s] ] ],
                 couplings = {(0,0,0):C.UVGC_1121_128})

V_307 = CTVertex(name = 'V_307',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.g, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVS1 ],
                 loop_particles = [ [ [P.b, P.g] ] ],
                 couplings = {(0,0,0):C.UVGC_1121_128})

V_308 = CTVertex(name = 'V_308',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.a, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS15 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1017_16})

V_309 = CTVertex(name = 'V_309',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS15 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1154_159})

V_310 = CTVertex(name = 'V_310',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.Z, P.Z, P.H ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS15, L.FFVVS4, L.FFVVS9 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1018_17,(0,0,0):C.UVGC_1159_166,(0,2,0):C.UVGC_1021_20})

V_311 = CTVertex(name = 'V_311',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.g, P.H ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS17 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1029_28})

V_312 = CTVertex(name = 'V_312',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.Z, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS51, L.FFVVS53 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1039_38,(0,1,0):C.UVGC_1048_47})

V_313 = CTVertex(name = 'V_313',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.g, P.W__minus__, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS37 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1093_91})

V_314 = CTVertex(name = 'V_314',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.g, P.W__plus__, P.H ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS68 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1094_92})

V_315 = CTVertex(name = 'V_315',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.a, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS11 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1016_15})

V_316 = CTVertex(name = 'V_316',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS11 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1153_158})

V_317 = CTVertex(name = 'V_317',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.Z, P.Z, P.G0 ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS11, L.FFVVS7, L.FFVVS8 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1019_18,(0,0,0):C.UVGC_1158_165,(0,2,0):C.UVGC_1020_19})

V_318 = CTVertex(name = 'V_318',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.g, P.G0 ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS13 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1028_27})

V_319 = CTVertex(name = 'V_319',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.Z, P.G0 ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS52, L.FFVVS54 ],
                 loop_particles = [ [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1038_37,(0,1,0):C.UVGC_1049_48})

V_320 = CTVertex(name = 'V_320',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.g, P.W__minus__, P.G0 ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS37 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1095_93})

V_321 = CTVertex(name = 'V_321',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.g, P.W__plus__, P.G0 ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS68 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1092_90})

V_322 = CTVertex(name = 'V_322',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.a, P.a, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS71 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1133_134,(0,0,1):C.UVGC_1133_135})

V_323 = CTVertex(name = 'V_323',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.a, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS71 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1151_154,(0,0,1):C.UVGC_1151_155})

V_324 = CTVertex(name = 'V_324',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.Z, P.Z, P.G__plus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS71 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1156_161,(0,0,1):C.UVGC_1156_162})

V_325 = CTVertex(name = 'V_325',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.a, P.g, P.G__plus__ ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS68, L.FFVVS75 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1001_2,(0,0,1):C.UVGC_1026_25})

V_326 = CTVertex(name = 'V_326',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.b, P.g, P.Z, P.G__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS68, L.FFVVS75 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1120_127,(0,1,1):C.UVGC_1046_45})

V_327 = CTVertex(name = 'V_327',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.W__minus__, P.G__plus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS68 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1096_94})

V_328 = CTVertex(name = 'V_328',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.a, P.a, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS31 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1134_136,(0,0,1):C.UVGC_1134_137})

V_329 = CTVertex(name = 'V_329',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.a, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS31 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1152_156,(0,0,1):C.UVGC_1152_157})

V_330 = CTVertex(name = 'V_330',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.Z, P.Z, P.G__minus__ ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.FFVVS31 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1157_163,(0,0,1):C.UVGC_1157_164})

V_331 = CTVertex(name = 'V_331',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.a, P.g, P.G__minus__ ],
                 color = [ 'T(4,2,1)' ],
                 lorentz = [ L.FFVVS25, L.FFVVS37 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1000_1,(0,1,1):C.UVGC_1027_26})

V_332 = CTVertex(name = 'V_332',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.t, P.g, P.Z, P.G__minus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS25, L.FFVVS37 ],
                 loop_particles = [ [ [P.b, P.g] ], [ [P.g, P.t] ] ],
                 couplings = {(0,1,0):C.UVGC_1119_126,(0,0,1):C.UVGC_1047_46})

V_333 = CTVertex(name = 'V_333',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.W__plus__, P.G__minus__ ],
                 color = [ 'T(3,2,1)' ],
                 lorentz = [ L.FFVVS37 ],
                 loop_particles = [ [ [P.b, P.g, P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_1097_95})

V_334 = CTVertex(name = 'V_334',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.g, P.g ],
                 color = [ 'T(3,-1,1)*T(4,-2,-1)*T(5,2,-2)', 'T(3,-1,1)*T(4,2,-2)*T(5,-2,-1)', 'T(3,-2,-1)*T(4,-1,1)*T(5,2,-2)', 'T(3,-2,-1)*T(4,2,-2)*T(5,-1,1)', 'T(3,2,-2)*T(4,-1,1)*T(5,-2,-1)', 'T(3,2,-2)*T(4,-2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV1, L.FFVVV13, L.FFVVV14, L.FFVVV15, L.FFVVV16, L.FFVVV17, L.FFVVV18, L.FFVVV2, L.FFVVV3, L.FFVVV4, L.FFVVV5, L.FFVVV6, L.FFVVV8 ],
                 loop_particles = [ [ [P.b, P.G__plus__] ], [ [P.G0, P.t], [P.H, P.t] ] ],
                 couplings = {(5,9,1):C.UVGC_1180_180,(4,9,1):C.UVGC_1179_179,(3,9,1):C.UVGC_1179_179,(5,11,1):C.UVGC_1179_179,(4,11,1):C.UVGC_1180_180,(2,11,1):C.UVGC_1179_179,(5,8,1):C.UVGC_1179_179,(3,8,1):C.UVGC_1180_180,(1,8,1):C.UVGC_1179_179,(3,10,1):C.UVGC_1179_179,(1,10,1):C.UVGC_1180_180,(0,10,1):C.UVGC_1179_179,(4,0,1):C.UVGC_1179_179,(2,0,1):C.UVGC_1180_180,(0,0,1):C.UVGC_1179_179,(1,7,1):C.UVGC_1179_179,(2,7,1):C.UVGC_1179_179,(0,12,1):C.UVGC_1180_180,(5,4,0):C.UVGC_1125_132,(4,4,0):C.UVGC_1124_131,(3,4,0):C.UVGC_1124_131,(5,6,0):C.UVGC_1124_131,(4,6,0):C.UVGC_1125_132,(2,6,0):C.UVGC_1124_131,(5,3,0):C.UVGC_1124_131,(3,3,0):C.UVGC_1125_132,(1,3,0):C.UVGC_1124_131,(3,5,0):C.UVGC_1124_131,(1,5,0):C.UVGC_1125_132,(0,5,0):C.UVGC_1124_131,(4,1,0):C.UVGC_1124_131,(2,1,0):C.UVGC_1125_132,(0,1,0):C.UVGC_1124_131,(1,2,0):C.UVGC_1124_131,(2,2,0):C.UVGC_1124_131,(0,2,0):C.UVGC_1125_132,(0,2,1):C.UVGC_1180_180})

V_335 = CTVertex(name = 'V_335',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.g, P.g, P.g ],
                 color = [ 'T(3,-1,1)*T(4,-2,-1)*T(5,2,-2)', 'T(3,-1,1)*T(4,2,-2)*T(5,-2,-1)', 'T(3,-2,-1)*T(4,-1,1)*T(5,2,-2)', 'T(3,-2,-1)*T(4,2,-2)*T(5,-1,1)', 'T(3,2,-2)*T(4,-1,1)*T(5,-2,-1)', 'T(3,2,-2)*T(4,-2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV10, L.FFVVV11, L.FFVVV12, L.FFVVV7, L.FFVVV8, L.FFVVV9 ],
                 loop_particles = [ [ [P.G__plus__, P.t] ] ],
                 couplings = {(5,0,0):C.UVGC_1125_132,(4,0,0):C.UVGC_1124_131,(3,0,0):C.UVGC_1124_131,(5,2,0):C.UVGC_1124_131,(4,2,0):C.UVGC_1125_132,(2,2,0):C.UVGC_1124_131,(5,5,0):C.UVGC_1124_131,(3,5,0):C.UVGC_1125_132,(1,5,0):C.UVGC_1124_131,(3,1,0):C.UVGC_1124_131,(1,1,0):C.UVGC_1125_132,(0,1,0):C.UVGC_1124_131,(4,3,0):C.UVGC_1124_131,(2,3,0):C.UVGC_1125_132,(0,3,0):C.UVGC_1124_131,(1,4,0):C.UVGC_1124_131,(2,4,0):C.UVGC_1124_131,(0,4,0):C.UVGC_1125_132})

V_336 = CTVertex(name = 'V_336',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.a, P.g, P.g ],
                 color = [ 'T(4,-1,1)*T(5,2,-1)', 'T(4,2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV12, L.FFVVV14, L.FFVVV16, L.FFVVV17, L.FFVVV18, L.FFVVV4, L.FFVVV5, L.FFVVV8 ],
                 loop_particles = [ [ [P.b, P.G__plus__] ], [ [P.G0, P.t], [P.H, P.t] ] ],
                 couplings = {(1,5,1):C.UVGC_1123_130,(0,5,1):C.UVGC_1122_129,(1,6,1):C.UVGC_1123_130,(0,6,1):C.UVGC_1122_129,(1,0,1):C.UVGC_1122_129,(0,0,1):C.UVGC_1123_130,(1,7,1):C.UVGC_1122_129,(0,7,1):C.UVGC_1123_130,(1,2,0):C.UVGC_1122_129,(0,2,0):C.UVGC_1123_130,(1,4,0):C.UVGC_1123_130,(1,4,1):C.UVGC_1122_129,(0,4,0):C.UVGC_1122_129,(0,4,1):C.UVGC_1123_130,(1,3,0):C.UVGC_1122_129,(0,3,0):C.UVGC_1123_130,(1,1,0):C.UVGC_1123_130,(1,1,1):C.UVGC_1122_129,(0,1,0):C.UVGC_1122_129,(0,1,1):C.UVGC_1123_130})

V_337 = CTVertex(name = 'V_337',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.a, P.g, P.g ],
                 color = [ 'T(4,-1,1)*T(5,2,-1)', 'T(4,2,-1)*T(5,-1,1)' ],
                 lorentz = [ L.FFVVV10, L.FFVVV11, L.FFVVV12, L.FFVVV8 ],
                 loop_particles = [ [ [P.G__plus__, P.t] ] ],
                 couplings = {(1,0,0):C.UVGC_1170_176,(0,0,0):C.UVGC_1169_175,(1,2,0):C.UVGC_1169_175,(0,2,0):C.UVGC_1170_176,(1,1,0):C.UVGC_1170_176,(0,1,0):C.UVGC_1169_175,(1,3,0):C.UVGC_1169_175,(0,3,0):C.UVGC_1170_176})

V_338 = CTVertex(name = 'V_338',
                 type = 'UV',
                 particles = [ P.t__tilde__, P.t, P.g, P.g, P.Z ],
                 color = [ 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVV1, L.FFVVV13, L.FFVVV14, L.FFVVV15, L.FFVVV16, L.FFVVV2, L.FFVVV3, L.FFVVV4 ],
                 loop_particles = [ [ [P.b, P.G__plus__] ], [ [P.G0, P.t], [P.H, P.t] ] ],
                 couplings = {(1,7,1):C.UVGC_1185_186,(0,7,1):C.UVGC_1184_185,(1,6,1):C.UVGC_1184_185,(0,6,1):C.UVGC_1185_186,(1,0,1):C.UVGC_1185_186,(0,0,1):C.UVGC_1184_185,(1,5,1):C.UVGC_1184_185,(0,5,1):C.UVGC_1185_186,(1,4,0):C.UVGC_1182_181,(1,4,1):C.UVGC_1182_182,(0,4,0):C.UVGC_1183_183,(0,4,1):C.UVGC_1183_184,(1,3,0):C.UVGC_1183_183,(1,3,1):C.UVGC_1183_184,(0,3,0):C.UVGC_1182_181,(0,3,1):C.UVGC_1182_182,(1,1,0):C.UVGC_1182_181,(1,1,1):C.UVGC_1182_182,(0,1,0):C.UVGC_1183_183,(0,1,1):C.UVGC_1183_184,(1,2,0):C.UVGC_1183_183,(1,2,1):C.UVGC_1183_184,(0,2,0):C.UVGC_1182_181,(0,2,1):C.UVGC_1182_182})

V_339 = CTVertex(name = 'V_339',
                 type = 'UV',
                 particles = [ P.b__tilde__, P.b, P.g, P.g, P.Z ],
                 color = [ 'T(3,-1,1)*T(4,2,-1)', 'T(3,2,-1)*T(4,-1,1)' ],
                 lorentz = [ L.FFVVV10, L.FFVVV7, L.FFVVV8, L.FFVVV9 ],
                 loop_particles = [ [ [P.G__plus__, P.t] ] ],
                 couplings = {(1,0,0):C.UVGC_1174_178,(0,0,0):C.UVGC_1173_177,(1,3,0):C.UVGC_1173_177,(0,3,0):C.UVGC_1174_178,(1,1,0):C.UVGC_1174_178,(0,1,0):C.UVGC_1173_177,(1,2,0):C.UVGC_1173_177,(0,2,0):C.UVGC_1174_178})

