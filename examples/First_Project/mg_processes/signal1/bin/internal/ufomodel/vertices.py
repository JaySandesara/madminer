# This file was automatically created by FeynRules 2.4.72
# Mathematica version: 11.3.0 for Mac OS X x86 (64-bit) (March 7, 2018)
# Date: Thu 8 Aug 2019 15:24:53


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.G0, P.G0, P.G0, P.G0, P.G0, P.G0 ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_85})

V_2 = Vertex(name = 'V_2',
             particles = [ P.G0, P.G0, P.G0, P.G0, P.G__minus__, P.G__plus__ ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_83})

V_3 = Vertex(name = 'V_3',
             particles = [ P.G0, P.G0, P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_82})

V_4 = Vertex(name = 'V_4',
             particles = [ P.G__minus__, P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__, P.G__plus__ ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_84})

V_5 = Vertex(name = 'V_5',
             particles = [ P.G0, P.G0, P.G0, P.G0, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_83})

V_6 = Vertex(name = 'V_6',
             particles = [ P.G0, P.G0, P.G__minus__, P.G__plus__, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_81})

V_7 = Vertex(name = 'V_7',
             particles = [ P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_82})

V_8 = Vertex(name = 'V_8',
             particles = [ P.G0, P.G0, P.H, P.H, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_83})

V_9 = Vertex(name = 'V_9',
             particles = [ P.G__minus__, P.G__plus__, P.H, P.H, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSSSS1 ],
             couplings = {(0,0):C.GC_83})

V_10 = Vertex(name = 'V_10',
              particles = [ P.H, P.H, P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSSS1 ],
              couplings = {(0,0):C.GC_85})

V_11 = Vertex(name = 'V_11',
              particles = [ P.G0, P.G0, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.SSSS1, L.SSSS6 ],
              couplings = {(0,0):C.GC_214,(0,1):C.GC_14})

V_12 = Vertex(name = 'V_12',
              particles = [ P.G0, P.G0, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_407})

V_13 = Vertex(name = 'V_13',
              particles = [ P.G0, P.G0, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.SSSS1, L.SSSS2 ],
              couplings = {(0,0):C.GC_210,(0,1):C.GC_13})

V_14 = Vertex(name = 'V_14',
              particles = [ P.G0, P.G0, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_405})

V_15 = Vertex(name = 'V_15',
              particles = [ P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.SSSS1, L.SSSS2, L.SSSS5 ],
              couplings = {(0,0):C.GC_212,(0,2):C.GC_71,(0,1):C.GC_89})

V_16 = Vertex(name = 'V_16',
              particles = [ P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_406})

V_17 = Vertex(name = 'V_17',
              particles = [ P.G0, P.G0, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1, L.SSSS2, L.SSSS3 ],
              couplings = {(0,0):C.GC_211,(0,2):C.GC_89,(0,1):C.GC_70})

V_18 = Vertex(name = 'V_18',
              particles = [ P.G0, P.G0, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_405})

V_19 = Vertex(name = 'V_19',
              particles = [ P.G0, P.G0, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_526})

V_20 = Vertex(name = 'V_20',
              particles = [ P.G__minus__, P.G__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1, L.SSSS2 ],
              couplings = {(0,0):C.GC_209,(0,1):C.GC_13})

V_21 = Vertex(name = 'V_21',
              particles = [ P.G__minus__, P.G__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_405})

V_22 = Vertex(name = 'V_22',
              particles = [ P.G__minus__, P.G__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_526})

V_23 = Vertex(name = 'V_23',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1, L.SSSS6 ],
              couplings = {(0,0):C.GC_213,(0,1):C.GC_14})

V_24 = Vertex(name = 'V_24',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_407})

V_25 = Vertex(name = 'V_25',
              particles = [ P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS1 ],
              couplings = {(0,0):C.GC_527})

V_26 = Vertex(name = 'V_26',
              particles = [ P.G0, P.G0, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1, L.SSS2, L.SSS4 ],
              couplings = {(0,0):C.GC_408,(0,1):C.GC_424,(0,2):C.GC_417})

V_27 = Vertex(name = 'V_27',
              particles = [ P.G0, P.G0, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_551})

V_28 = Vertex(name = 'V_28',
              particles = [ P.G__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1, L.SSS4 ],
              couplings = {(0,0):C.GC_408,(0,1):C.GC_538})

V_29 = Vertex(name = 'V_29',
              particles = [ P.G__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_550})

V_30 = Vertex(name = 'V_30',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1, L.SSS5 ],
              couplings = {(0,0):C.GC_409,(0,1):C.GC_539})

V_31 = Vertex(name = 'V_31',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_552})

V_32 = Vertex(name = 'V_32',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS1 ],
              couplings = {(0,0):C.GC_537})

V_33 = Vertex(name = 'V_33',
              particles = [ P.G0, P.G0, P.G0, P.G0, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSS1 ],
              couplings = {(0,0):C.GC_420})

V_34 = Vertex(name = 'V_34',
              particles = [ P.G0, P.G0, P.G__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSS1 ],
              couplings = {(0,0):C.GC_418})

V_35 = Vertex(name = 'V_35',
              particles = [ P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSS1 ],
              couplings = {(0,0):C.GC_419})

V_36 = Vertex(name = 'V_36',
              particles = [ P.G0, P.G0, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSS1 ],
              couplings = {(0,0):C.GC_420})

V_37 = Vertex(name = 'V_37',
              particles = [ P.G__minus__, P.G__plus__, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSS1 ],
              couplings = {(0,0):C.GC_420})

V_38 = Vertex(name = 'V_38',
              particles = [ P.H, P.H, P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSSS1 ],
              couplings = {(0,0):C.GC_421})

V_39 = Vertex(name = 'V_39',
              particles = [ P.a, P.a, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_398,(0,0):C.GC_6})

V_40 = Vertex(name = 'V_40',
              particles = [ P.a, P.a, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_635})

V_41 = Vertex(name = 'V_41',
              particles = [ P.a, P.a, P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSSSS1 ],
              couplings = {(0,0):C.GC_190})

V_42 = Vertex(name = 'V_42',
              particles = [ P.a, P.G0, P.G__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VSSSS7 ],
              couplings = {(0,0):C.GC_173})

V_43 = Vertex(name = 'V_43',
              particles = [ P.a, P.G0, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VSSS1 ],
              couplings = {(0,0):C.GC_450})

V_44 = Vertex(name = 'V_44',
              particles = [ P.a, P.G0, P.H ],
              color = [ '1' ],
              lorentz = [ L.VSS2 ],
              couplings = {(0,0):C.GC_600})

V_45 = Vertex(name = 'V_45',
              particles = [ P.a, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VSS2 ],
              couplings = {(0,0):C.GC_3})

V_46 = Vertex(name = 'V_46',
              particles = [ P.a, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VSS2 ],
              couplings = {(0,0):C.GC_624})

V_47 = Vertex(name = 'V_47',
              particles = [ P.a, P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VSSSS8 ],
              couplings = {(0,0):C.GC_172})

V_48 = Vertex(name = 'V_48',
              particles = [ P.G0, P.G__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSSS4 ],
              couplings = {(0,0):C.GC_88})

V_49 = Vertex(name = 'V_49',
              particles = [ P.G0, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.SSS3 ],
              couplings = {(0,0):C.GC_423})

V_50 = Vertex(name = 'V_50',
              particles = [ P.a, P.a, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.VVSS5 ],
              couplings = {(0,0):C.GC_397})

V_51 = Vertex(name = 'V_51',
              particles = [ P.a, P.a, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS5 ],
              couplings = {(0,0):C.GC_397})

V_52 = Vertex(name = 'V_52',
              particles = [ P.a, P.a, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS4 ],
              couplings = {(0,0):C.GC_590})

V_53 = Vertex(name = 'V_53',
              particles = [ P.g, P.g, P.G0, P.G0 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.VVSS5 ],
              couplings = {(0,0):C.GC_92})

V_54 = Vertex(name = 'V_54',
              particles = [ P.g, P.g, P.G__minus__, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.VVSS5 ],
              couplings = {(0,0):C.GC_92})

V_55 = Vertex(name = 'V_55',
              particles = [ P.g, P.g, P.H, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.VVSS5 ],
              couplings = {(0,0):C.GC_92})

V_56 = Vertex(name = 'V_56',
              particles = [ P.g, P.g, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.VVS4 ],
              couplings = {(0,0):C.GC_426})

V_57 = Vertex(name = 'V_57',
              particles = [ P.a, P.W__minus__, P.G0, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_153,(0,0):C.GC_247})

V_58 = Vertex(name = 'V_58',
              particles = [ P.a, P.W__minus__, P.G0, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_627})

V_59 = Vertex(name = 'V_59',
              particles = [ P.a, P.W__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_152,(0,0):C.GC_246})

V_60 = Vertex(name = 'V_60',
              particles = [ P.a, P.W__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_626})

V_61 = Vertex(name = 'V_61',
              particles = [ P.a, P.W__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3, L.VVS4 ],
              couplings = {(0,1):C.GC_438,(0,0):C.GC_471})

V_62 = Vertex(name = 'V_62',
              particles = [ P.a, P.W__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_668})

V_63 = Vertex(name = 'V_63',
              particles = [ P.a, P.W__plus__, P.G0, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_153,(0,0):C.GC_247})

V_64 = Vertex(name = 'V_64',
              particles = [ P.a, P.W__plus__, P.G0, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_627})

V_65 = Vertex(name = 'V_65',
              particles = [ P.a, P.W__plus__, P.G__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_154,(0,0):C.GC_248})

V_66 = Vertex(name = 'V_66',
              particles = [ P.a, P.W__plus__, P.G__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_628})

V_67 = Vertex(name = 'V_67',
              particles = [ P.a, P.W__plus__, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3, L.VVS4 ],
              couplings = {(0,1):C.GC_439,(0,0):C.GC_472})

V_68 = Vertex(name = 'V_68',
              particles = [ P.a, P.W__plus__, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_673})

V_69 = Vertex(name = 'V_69',
              particles = [ P.W__minus__, P.W__plus__, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_101,(0,0):C.GC_225})

V_70 = Vertex(name = 'V_70',
              particles = [ P.W__minus__, P.W__plus__, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_602})

V_71 = Vertex(name = 'V_71',
              particles = [ P.W__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_101,(0,0):C.GC_225})

V_72 = Vertex(name = 'V_72',
              particles = [ P.W__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_604})

V_73 = Vertex(name = 'V_73',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_101,(0,0):C.GC_225})

V_74 = Vertex(name = 'V_74',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_603})

V_75 = Vertex(name = 'V_75',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3, L.VVS4 ],
              couplings = {(0,1):C.GC_431,(0,0):C.GC_463})

V_76 = Vertex(name = 'V_76',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_666})

V_77 = Vertex(name = 'V_77',
              particles = [ P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVV3, L.VVV4, L.VVV5 ],
              couplings = {(0,2):C.GC_321,(0,0):C.GC_4,(0,1):C.GC_534})

V_78 = Vertex(name = 'V_78',
              particles = [ P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVV3 ],
              couplings = {(0,0):C.GC_598})

V_79 = Vertex(name = 'V_79',
              particles = [ P.a, P.Z, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_399,(0,0):C.GC_659})

V_80 = Vertex(name = 'V_80',
              particles = [ P.a, P.Z, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_400,(0,0):C.GC_345})

V_81 = Vertex(name = 'V_81',
              particles = [ P.a, P.Z, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_660})

V_82 = Vertex(name = 'V_82',
              particles = [ P.a, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_399,(0,0):C.GC_659})

V_83 = Vertex(name = 'V_83',
              particles = [ P.a, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS3, L.VVS4 ],
              couplings = {(0,1):C.GC_591,(0,0):C.GC_681})

V_84 = Vertex(name = 'V_84',
              particles = [ P.W__minus__, P.Z, P.G0, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_319,(0,0):C.GC_8})

V_85 = Vertex(name = 'V_85',
              particles = [ P.W__minus__, P.Z, P.G0, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_630})

V_86 = Vertex(name = 'V_86',
              particles = [ P.W__minus__, P.Z, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_318,(0,0):C.GC_7})

V_87 = Vertex(name = 'V_87',
              particles = [ P.W__minus__, P.Z, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_631})

V_88 = Vertex(name = 'V_88',
              particles = [ P.W__minus__, P.Z, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3, L.VVS4 ],
              couplings = {(0,1):C.GC_513,(0,0):C.GC_410})

V_89 = Vertex(name = 'V_89',
              particles = [ P.W__minus__, P.Z, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_679})

V_90 = Vertex(name = 'V_90',
              particles = [ P.W__plus__, P.Z, P.G0, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_319,(0,0):C.GC_8})

V_91 = Vertex(name = 'V_91',
              particles = [ P.W__plus__, P.Z, P.G0, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_630})

V_92 = Vertex(name = 'V_92',
              particles = [ P.W__plus__, P.Z, P.G__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_320,(0,0):C.GC_9})

V_93 = Vertex(name = 'V_93',
              particles = [ P.W__plus__, P.Z, P.G__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_629})

V_94 = Vertex(name = 'V_94',
              particles = [ P.W__plus__, P.Z, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3, L.VVS4 ],
              couplings = {(0,1):C.GC_514,(0,0):C.GC_411})

V_95 = Vertex(name = 'V_95',
              particles = [ P.W__plus__, P.Z, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS3 ],
              couplings = {(0,0):C.GC_674})

V_96 = Vertex(name = 'V_96',
              particles = [ P.Z, P.Z, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_396,(0,0):C.GC_394})

V_97 = Vertex(name = 'V_97',
              particles = [ P.Z, P.Z, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_661})

V_98 = Vertex(name = 'V_98',
              particles = [ P.Z, P.Z, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4, L.VVSS5 ],
              couplings = {(0,1):C.GC_395,(0,0):C.GC_393})

V_99 = Vertex(name = 'V_99',
              particles = [ P.Z, P.Z, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS4 ],
              couplings = {(0,0):C.GC_662})

V_100 = Vertex(name = 'V_100',
               particles = [ P.Z, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS4, L.VVSS5 ],
               couplings = {(0,1):C.GC_396,(0,0):C.GC_394})

V_101 = Vertex(name = 'V_101',
               particles = [ P.Z, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS4 ],
               couplings = {(0,0):C.GC_663})

V_102 = Vertex(name = 'V_102',
               particles = [ P.Z, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS3, L.VVS4 ],
               couplings = {(0,1):C.GC_589,(0,0):C.GC_588})

V_103 = Vertex(name = 'V_103',
               particles = [ P.Z, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS3 ],
               couplings = {(0,0):C.GC_683})

V_104 = Vertex(name = 'V_104',
               particles = [ P.W__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVV2, L.VVV3, L.VVV5 ],
               couplings = {(0,2):C.GC_155,(0,1):C.GC_244,(0,0):C.GC_528})

V_105 = Vertex(name = 'V_105',
               particles = [ P.W__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVV3 ],
               couplings = {(0,0):C.GC_615})

V_106 = Vertex(name = 'V_106',
               particles = [ P.ghA, P.ghWm__tilde__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_669})

V_107 = Vertex(name = 'V_107',
               particles = [ P.ghA, P.ghWm__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_3})

V_108 = Vertex(name = 'V_108',
               particles = [ P.ghA, P.ghWm__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_599})

V_109 = Vertex(name = 'V_109',
               particles = [ P.ghA, P.ghWp__tilde__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_671})

V_110 = Vertex(name = 'V_110',
               particles = [ P.ghA, P.ghWp__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_4})

V_111 = Vertex(name = 'V_111',
               particles = [ P.ghA, P.ghWp__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_598})

V_112 = Vertex(name = 'V_112',
               particles = [ P.ghA, P.ghZ__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_680})

V_113 = Vertex(name = 'V_113',
               particles = [ P.ghWm, P.ghA__tilde__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_472})

V_114 = Vertex(name = 'V_114',
               particles = [ P.ghWm, P.ghA__tilde__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_672})

V_115 = Vertex(name = 'V_115',
               particles = [ P.ghWm, P.ghA__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_3})

V_116 = Vertex(name = 'V_116',
               particles = [ P.ghWm, P.ghA__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_599})

V_117 = Vertex(name = 'V_117',
               particles = [ P.ghWm, P.ghWm__tilde__, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_461})

V_118 = Vertex(name = 'V_118',
               particles = [ P.ghWm, P.ghWm__tilde__, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_667})

V_119 = Vertex(name = 'V_119',
               particles = [ P.ghWm, P.ghWm__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_462})

V_120 = Vertex(name = 'V_120',
               particles = [ P.ghWm, P.ghWm__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_665})

V_121 = Vertex(name = 'V_121',
               particles = [ P.ghWm, P.ghWm__tilde__, P.a ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_4})

V_122 = Vertex(name = 'V_122',
               particles = [ P.ghWm, P.ghWm__tilde__, P.a ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_622})

V_123 = Vertex(name = 'V_123',
               particles = [ P.ghWm, P.ghWm__tilde__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_244})

V_124 = Vertex(name = 'V_124',
               particles = [ P.ghWm, P.ghWm__tilde__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_616})

V_125 = Vertex(name = 'V_125',
               particles = [ P.ghWm, P.ghZ__tilde__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_555})

V_126 = Vertex(name = 'V_126',
               particles = [ P.ghWm, P.ghZ__tilde__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_678})

V_127 = Vertex(name = 'V_127',
               particles = [ P.ghWm, P.ghZ__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_244})

V_128 = Vertex(name = 'V_128',
               particles = [ P.ghWm, P.ghZ__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_615})

V_129 = Vertex(name = 'V_129',
               particles = [ P.ghWp, P.ghA__tilde__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_471})

V_130 = Vertex(name = 'V_130',
               particles = [ P.ghWp, P.ghA__tilde__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_670})

V_131 = Vertex(name = 'V_131',
               particles = [ P.ghWp, P.ghA__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_4})

V_132 = Vertex(name = 'V_132',
               particles = [ P.ghWp, P.ghA__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_598})

V_133 = Vertex(name = 'V_133',
               particles = [ P.ghWp, P.ghWp__tilde__, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_464})

V_134 = Vertex(name = 'V_134',
               particles = [ P.ghWp, P.ghWp__tilde__, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_664})

V_135 = Vertex(name = 'V_135',
               particles = [ P.ghWp, P.ghWp__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_462})

V_136 = Vertex(name = 'V_136',
               particles = [ P.ghWp, P.ghWp__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_665})

V_137 = Vertex(name = 'V_137',
               particles = [ P.ghWp, P.ghWp__tilde__, P.a ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_3})

V_138 = Vertex(name = 'V_138',
               particles = [ P.ghWp, P.ghWp__tilde__, P.a ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_625})

V_139 = Vertex(name = 'V_139',
               particles = [ P.ghWp, P.ghWp__tilde__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_245})

V_140 = Vertex(name = 'V_140',
               particles = [ P.ghWp, P.ghWp__tilde__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_619})

V_141 = Vertex(name = 'V_141',
               particles = [ P.ghWp, P.ghZ__tilde__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_554})

V_142 = Vertex(name = 'V_142',
               particles = [ P.ghWp, P.ghZ__tilde__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_676})

V_143 = Vertex(name = 'V_143',
               particles = [ P.ghWp, P.ghZ__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_245})

V_144 = Vertex(name = 'V_144',
               particles = [ P.ghWp, P.ghZ__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_618})

V_145 = Vertex(name = 'V_145',
               particles = [ P.ghZ, P.ghA__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_680})

V_146 = Vertex(name = 'V_146',
               particles = [ P.ghZ, P.ghWm__tilde__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_556})

V_147 = Vertex(name = 'V_147',
               particles = [ P.ghZ, P.ghWm__tilde__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_675})

V_148 = Vertex(name = 'V_148',
               particles = [ P.ghZ, P.ghWm__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_244})

V_149 = Vertex(name = 'V_149',
               particles = [ P.ghZ, P.ghWm__tilde__, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_615})

V_150 = Vertex(name = 'V_150',
               particles = [ P.ghZ, P.ghWp__tilde__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_553})

V_151 = Vertex(name = 'V_151',
               particles = [ P.ghZ, P.ghWp__tilde__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_677})

V_152 = Vertex(name = 'V_152',
               particles = [ P.ghZ, P.ghWp__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_245})

V_153 = Vertex(name = 'V_153',
               particles = [ P.ghZ, P.ghWp__tilde__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_618})

V_154 = Vertex(name = 'V_154',
               particles = [ P.ghZ, P.ghZ__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_587})

V_155 = Vertex(name = 'V_155',
               particles = [ P.ghZ, P.ghZ__tilde__, P.H ],
               color = [ '1' ],
               lorentz = [ L.UUS1 ],
               couplings = {(0,0):C.GC_682})

V_156 = Vertex(name = 'V_156',
               particles = [ P.g, P.g, P.g ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVV3, L.VVV5 ],
               couplings = {(0,1):C.GC_193,(0,0):C.GC_10})

V_157 = Vertex(name = 'V_157',
               particles = [ P.ghG, P.ghG__tilde__, P.g ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.UUV1 ],
               couplings = {(0,0):C.GC_10})

V_158 = Vertex(name = 'V_158',
               particles = [ P.g, P.g, P.g, P.G0, P.G0 ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVVSS5 ],
               couplings = {(0,0):C.GC_194})

V_159 = Vertex(name = 'V_159',
               particles = [ P.g, P.g, P.g, P.G__minus__, P.G__plus__ ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVVSS5 ],
               couplings = {(0,0):C.GC_194})

V_160 = Vertex(name = 'V_160',
               particles = [ P.g, P.g, P.g, P.H, P.H ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVVSS5 ],
               couplings = {(0,0):C.GC_194})

V_161 = Vertex(name = 'V_161',
               particles = [ P.g, P.g, P.g, P.H ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVVS9 ],
               couplings = {(0,0):C.GC_457})

V_162 = Vertex(name = 'V_162',
               particles = [ P.g, P.g, P.g, P.g ],
               color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
               lorentz = [ L.VVVV12, L.VVVV13, L.VVVV4, L.VVVV6, L.VVVV8, L.VVVV9 ],
               couplings = {(0,3):C.GC_199,(1,1):C.GC_199,(2,0):C.GC_199,(1,4):C.GC_12,(0,2):C.GC_12,(2,5):C.GC_12})

V_163 = Vertex(name = 'V_163',
               particles = [ P.g, P.g, P.g, P.g, P.G0, P.G0 ],
               color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
               lorentz = [ L.VVVVSS1, L.VVVVSS3, L.VVVVSS4 ],
               couplings = {(1,1):C.GC_200,(0,0):C.GC_200,(2,2):C.GC_200})

V_164 = Vertex(name = 'V_164',
               particles = [ P.g, P.g, P.g, P.g, P.G__minus__, P.G__plus__ ],
               color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
               lorentz = [ L.VVVVSS1, L.VVVVSS3, L.VVVVSS4 ],
               couplings = {(1,1):C.GC_200,(0,0):C.GC_200,(2,2):C.GC_200})

V_165 = Vertex(name = 'V_165',
               particles = [ P.g, P.g, P.g, P.g, P.H, P.H ],
               color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
               lorentz = [ L.VVVVSS1, L.VVVVSS3, L.VVVVSS4 ],
               couplings = {(1,1):C.GC_200,(0,0):C.GC_200,(2,2):C.GC_200})

V_166 = Vertex(name = 'V_166',
               particles = [ P.g, P.g, P.g, P.g, P.H ],
               color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
               lorentz = [ L.VVVVS4, L.VVVVS7, L.VVVVS8 ],
               couplings = {(1,1):C.GC_459,(0,0):C.GC_459,(2,2):C.GC_459})

V_167 = Vertex(name = 'V_167',
               particles = [ P.g, P.g, P.g, P.g, P.g ],
               color = [ 'f(-2,1,2)*f(-1,-2,3)*f(4,5,-1)', 'f(-2,1,2)*f(-1,-2,4)*f(3,5,-1)', 'f(-2,1,2)*f(-1,-2,5)*f(3,4,-1)', 'f(-2,1,3)*f(-1,-2,2)*f(4,5,-1)', 'f(-2,1,3)*f(-1,-2,4)*f(2,5,-1)', 'f(-2,1,3)*f(-1,-2,5)*f(2,4,-1)', 'f(-2,1,4)*f(-1,-2,2)*f(3,5,-1)', 'f(-2,1,4)*f(-1,-2,3)*f(2,5,-1)', 'f(-2,1,4)*f(-1,-2,5)*f(2,3,-1)', 'f(-2,1,5)*f(-1,-2,2)*f(3,4,-1)', 'f(-2,1,5)*f(-1,-2,3)*f(2,4,-1)', 'f(-2,1,5)*f(-1,-2,4)*f(2,3,-1)', 'f(-2,2,3)*f(-1,-2,1)*f(4,5,-1)', 'f(-2,2,3)*f(-1,-2,4)*f(1,5,-1)', 'f(-2,2,3)*f(-1,-2,5)*f(1,4,-1)', 'f(-2,2,4)*f(-1,-2,1)*f(3,5,-1)', 'f(-2,2,4)*f(-1,-2,3)*f(1,5,-1)', 'f(-2,2,4)*f(-1,-2,5)*f(1,3,-1)', 'f(-2,2,5)*f(-1,-2,1)*f(3,4,-1)', 'f(-2,2,5)*f(-1,-2,3)*f(1,4,-1)', 'f(-2,2,5)*f(-1,-2,4)*f(1,3,-1)', 'f(-2,3,4)*f(-1,-2,1)*f(2,5,-1)', 'f(-2,3,4)*f(-1,-2,2)*f(1,5,-1)', 'f(-2,3,4)*f(-1,-2,5)*f(1,2,-1)', 'f(-2,3,5)*f(-1,-2,1)*f(2,4,-1)', 'f(-2,3,5)*f(-1,-2,2)*f(1,4,-1)', 'f(-2,3,5)*f(-1,-2,4)*f(1,2,-1)', 'f(-2,4,5)*f(-1,-2,1)*f(2,3,-1)', 'f(-2,4,5)*f(-1,-2,2)*f(1,3,-1)', 'f(-2,4,5)*f(-1,-2,3)*f(1,2,-1)' ],
               lorentz = [ L.VVVVV1, L.VVVVV10, L.VVVVV11, L.VVVVV12, L.VVVVV13, L.VVVVV14, L.VVVVV15, L.VVVVV17, L.VVVVV18, L.VVVVV2, L.VVVVV4, L.VVVVV5, L.VVVVV6, L.VVVVV7, L.VVVVV8 ],
               couplings = {(24,11):C.GC_206,(21,12):C.GC_205,(18,12):C.GC_206,(15,11):C.GC_205,(28,9):C.GC_206,(22,2):C.GC_206,(9,2):C.GC_205,(3,9):C.GC_205,(29,10):C.GC_206,(16,3):C.GC_206,(10,3):C.GC_205,(0,10):C.GC_205,(26,6):C.GC_205,(20,5):C.GC_205,(4,5):C.GC_206,(1,6):C.GC_206,(25,1):C.GC_206,(6,1):C.GC_205,(19,4):C.GC_206,(7,4):C.GC_205,(23,8):C.GC_205,(17,7):C.GC_205,(5,7):C.GC_206,(2,8):C.GC_206,(27,0):C.GC_206,(12,0):C.GC_205,(13,13):C.GC_206,(11,13):C.GC_205,(14,14):C.GC_205,(8,14):C.GC_206})

V_168 = Vertex(name = 'V_168',
               particles = [ P.g, P.g, P.g, P.g, P.g, P.g ],
               color = [ 'f(-3,1,2)*f(-2,3,4)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,1,2)*f(-2,3,5)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,1,2)*f(-2,3,6)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,1,2)*f(-2,4,5)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,1,2)*f(-2,4,6)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,1,2)*f(-2,5,6)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,1,3)*f(-2,2,4)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,1,3)*f(-2,2,5)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,1,3)*f(-2,2,6)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,1,3)*f(-2,4,5)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,1,3)*f(-2,4,6)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,1,3)*f(-2,5,6)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,1,4)*f(-2,2,3)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,1,4)*f(-2,2,5)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,1,4)*f(-2,2,6)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,1,4)*f(-2,3,5)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,1,4)*f(-2,3,6)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,1,4)*f(-2,5,6)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,1,5)*f(-2,2,3)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,1,5)*f(-2,2,4)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,1,5)*f(-2,2,6)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,1,5)*f(-2,3,4)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,1,5)*f(-2,3,6)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,1,5)*f(-2,4,6)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,1,6)*f(-2,2,3)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,1,6)*f(-2,2,4)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,1,6)*f(-2,2,5)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,1,6)*f(-2,3,4)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,1,6)*f(-2,3,5)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,1,6)*f(-2,4,5)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,2,3)*f(-2,1,4)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,2,3)*f(-2,1,5)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,2,3)*f(-2,1,6)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,2,3)*f(-2,4,5)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,2,3)*f(-2,4,6)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,2,3)*f(-2,5,6)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,2,4)*f(-2,1,3)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,2,4)*f(-2,1,5)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,2,4)*f(-2,1,6)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,2,4)*f(-2,3,5)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,2,4)*f(-2,3,6)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,2,4)*f(-2,5,6)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,2,5)*f(-2,1,3)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,2,5)*f(-2,1,4)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,2,5)*f(-2,1,6)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,2,5)*f(-2,3,4)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,2,5)*f(-2,3,6)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,2,5)*f(-2,4,6)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,2,6)*f(-2,1,3)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,2,6)*f(-2,1,4)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,2,6)*f(-2,1,5)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,2,6)*f(-2,3,4)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,2,6)*f(-2,3,5)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,2,6)*f(-2,4,5)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,3,4)*f(-2,1,2)*f(-1,-2,-3)*f(5,6,-1)', 'f(-3,3,4)*f(-2,1,5)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,3,4)*f(-2,1,6)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,3,4)*f(-2,2,5)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,3,4)*f(-2,2,6)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,3,4)*f(-2,5,6)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,3,5)*f(-2,1,2)*f(-1,-2,-3)*f(4,6,-1)', 'f(-3,3,5)*f(-2,1,4)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,3,5)*f(-2,1,6)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,3,5)*f(-2,2,4)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,3,5)*f(-2,2,6)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,3,5)*f(-2,4,6)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,3,6)*f(-2,1,2)*f(-1,-2,-3)*f(4,5,-1)', 'f(-3,3,6)*f(-2,1,4)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,3,6)*f(-2,1,5)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,3,6)*f(-2,2,4)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,3,6)*f(-2,2,5)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,3,6)*f(-2,4,5)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,4,5)*f(-2,1,2)*f(-1,-2,-3)*f(3,6,-1)', 'f(-3,4,5)*f(-2,1,3)*f(-1,-2,-3)*f(2,6,-1)', 'f(-3,4,5)*f(-2,1,6)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,4,5)*f(-2,2,3)*f(-1,-2,-3)*f(1,6,-1)', 'f(-3,4,5)*f(-2,2,6)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,4,5)*f(-2,3,6)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,4,6)*f(-2,1,2)*f(-1,-2,-3)*f(3,5,-1)', 'f(-3,4,6)*f(-2,1,3)*f(-1,-2,-3)*f(2,5,-1)', 'f(-3,4,6)*f(-2,1,5)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,4,6)*f(-2,2,3)*f(-1,-2,-3)*f(1,5,-1)', 'f(-3,4,6)*f(-2,2,5)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,4,6)*f(-2,3,5)*f(-1,-2,-3)*f(1,2,-1)', 'f(-3,5,6)*f(-2,1,2)*f(-1,-2,-3)*f(3,4,-1)', 'f(-3,5,6)*f(-2,1,3)*f(-1,-2,-3)*f(2,4,-1)', 'f(-3,5,6)*f(-2,1,4)*f(-1,-2,-3)*f(2,3,-1)', 'f(-3,5,6)*f(-2,2,3)*f(-1,-2,-3)*f(1,4,-1)', 'f(-3,5,6)*f(-2,2,4)*f(-1,-2,-3)*f(1,3,-1)', 'f(-3,5,6)*f(-2,3,4)*f(-1,-2,-3)*f(1,2,-1)' ],
               lorentz = [ L.VVVVVV1, L.VVVVVV10, L.VVVVVV11, L.VVVVVV12, L.VVVVVV13, L.VVVVVV14, L.VVVVVV15, L.VVVVVV16, L.VVVVVV2, L.VVVVVV3, L.VVVVVV4, L.VVVVVV5, L.VVVVVV6, L.VVVVVV7, L.VVVVVV9 ],
               couplings = {(65,10):C.GC_208,(71,12):C.GC_207,(77,12):C.GC_208,(83,10):C.GC_207,(41,8):C.GC_208,(53,2):C.GC_208,(76,2):C.GC_207,(88,8):C.GC_207,(35,9):C.GC_208,(52,5):C.GC_208,(64,5):C.GC_207,(87,9):C.GC_207,(34,4):C.GC_207,(40,3):C.GC_207,(69,3):C.GC_208,(81,4):C.GC_208,(17,9):C.GC_207,(23,4):C.GC_208,(80,4):C.GC_207,(86,9):C.GC_208,(11,8):C.GC_207,(22,3):C.GC_208,(68,3):C.GC_207,(85,8):C.GC_208,(9,2):C.GC_207,(15,5):C.GC_207,(61,5):C.GC_208,(73,2):C.GC_208,(4,10):C.GC_207,(14,5):C.GC_208,(49,5):C.GC_207,(78,10):C.GC_208,(3,12):C.GC_208,(19,3):C.GC_207,(37,3):C.GC_208,(72,12):C.GC_207,(2,12):C.GC_207,(8,2):C.GC_208,(48,2):C.GC_207,(66,12):C.GC_208,(1,10):C.GC_208,(18,4):C.GC_207,(31,4):C.GC_208,(60,10):C.GC_207,(6,8):C.GC_208,(12,9):C.GC_208,(30,9):C.GC_207,(36,8):C.GC_207,(47,14):C.GC_208,(82,14):C.GC_207,(46,6):C.GC_208,(70,6):C.GC_207,(33,7):C.GC_207,(39,1):C.GC_207,(63,1):C.GC_208,(75,7):C.GC_208,(29,7):C.GC_208,(74,7):C.GC_207,(28,1):C.GC_208,(62,1):C.GC_207,(10,14):C.GC_207,(16,6):C.GC_207,(67,6):C.GC_208,(79,14):C.GC_208,(25,1):C.GC_207,(38,1):C.GC_208,(13,6):C.GC_208,(43,6):C.GC_207,(24,7):C.GC_207,(32,7):C.GC_208,(7,14):C.GC_208,(42,14):C.GC_207,(59,0):C.GC_208,(89,0):C.GC_207,(51,11):C.GC_208,(58,11):C.GC_207,(21,11):C.GC_207,(55,11):C.GC_208,(5,0):C.GC_207,(20,11):C.GC_208,(50,11):C.GC_207,(84,0):C.GC_208,(0,0):C.GC_208,(54,0):C.GC_207,(45,13):C.GC_207,(57,13):C.GC_208,(27,13):C.GC_208,(56,13):C.GC_207,(26,13):C.GC_207,(44,13):C.GC_208})

V_169 = Vertex(name = 'V_169',
               particles = [ P.b__tilde__, P.t, P.G0, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS3 ],
               couplings = {(0,0):C.GC_134})

V_170 = Vertex(name = 'V_170',
               particles = [ P.b__tilde__, P.t, P.G__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS3 ],
               couplings = {(0,0):C.GC_133})

V_171 = Vertex(name = 'V_171',
               particles = [ P.b__tilde__, P.t, P.G__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS3 ],
               couplings = {(0,0):C.GC_134})

V_172 = Vertex(name = 'V_172',
               particles = [ P.b__tilde__, P.t, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS2, L.FFSS3 ],
               couplings = {(0,0):C.GC_432,(0,1):C.GC_63})

V_173 = Vertex(name = 'V_173',
               particles = [ P.t__tilde__, P.t, P.G0, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1 ],
               couplings = {(0,0):C.GC_140})

V_174 = Vertex(name = 'V_174',
               particles = [ P.t__tilde__, P.t, P.G0, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1 ],
               couplings = {(0,0):C.GC_139})

V_175 = Vertex(name = 'V_175',
               particles = [ P.t__tilde__, P.t, P.G0, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS2 ],
               couplings = {(0,0):C.GC_137})

V_176 = Vertex(name = 'V_176',
               particles = [ P.t__tilde__, P.t, P.G__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS2 ],
               couplings = {(0,0):C.GC_137})

V_177 = Vertex(name = 'V_177',
               particles = [ P.t__tilde__, P.t, P.G0, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS1 ],
               couplings = {(0,0):C.GC_139})

V_178 = Vertex(name = 'V_178',
               particles = [ P.t__tilde__, P.t, P.H, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS2 ],
               couplings = {(0,0):C.GC_138})

V_179 = Vertex(name = 'V_179',
               particles = [ P.t__tilde__, P.t, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_435})

V_180 = Vertex(name = 'V_180',
               particles = [ P.t__tilde__, P.t, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS6, L.FFSS7 ],
               couplings = {(0,0):C.GC_27,(0,2):C.GC_95,(0,1):C.GC_435})

V_181 = Vertex(name = 'V_181',
               particles = [ P.t__tilde__, P.t, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS4, L.FFSS7 ],
               couplings = {(0,0):C.GC_29,(0,2):C.GC_96,(0,1):C.GC_434})

V_182 = Vertex(name = 'V_182',
               particles = [ P.t__tilde__, P.t, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_436})

V_183 = Vertex(name = 'V_183',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3 ],
               couplings = {(0,0):C.GC_691})

V_184 = Vertex(name = 'V_184',
               particles = [ P.t__tilde__, P.t, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3 ],
               couplings = {(0,0):C.GC_684})

V_185 = Vertex(name = 'V_185',
               particles = [ P.vt__tilde__, P.ta__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF1, L.FFFF3, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF8 ],
               couplings = {(0,1):C.GC_107,(0,4):C.GC_130,(0,2):C.GC_129,(0,3):C.GC_129,(0,0):C.GC_127,(0,5):C.GC_69})

V_186 = Vertex(name = 'V_186',
               particles = [ P.ta__plus__, P.ta__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF1, L.FFFF10, L.FFFF14, L.FFFF15, L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF4, L.FFFF5, L.FFFF6, L.FFFF9 ],
               couplings = {(0,7):C.GC_110,(0,10):C.GC_131,(0,8):C.GC_128,(0,9):C.GC_128,(0,0):C.GC_126,(0,1):C.GC_117,(0,2):C.GC_125,(0,3):C.GC_144,(0,6):C.GC_131,(0,4):C.GC_128,(0,5):C.GC_128,(0,11):C.GC_126})

V_187 = Vertex(name = 'V_187',
               particles = [ P.t__tilde__, P.b, P.G0, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS4 ],
               couplings = {(0,0):C.GC_135})

V_188 = Vertex(name = 'V_188',
               particles = [ P.t__tilde__, P.b, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS4 ],
               couplings = {(0,0):C.GC_136})

V_189 = Vertex(name = 'V_189',
               particles = [ P.t__tilde__, P.b, P.G__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSSS4 ],
               couplings = {(0,0):C.GC_135})

V_190 = Vertex(name = 'V_190',
               particles = [ P.t__tilde__, P.b, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS5 ],
               couplings = {(0,1):C.GC_433,(0,0):C.GC_63})

V_191 = Vertex(name = 'V_191',
               particles = [ P.ta__plus__, P.ta__minus__, P.b__tilde__, P.b ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF10, L.FFFF3, L.FFFF7, L.FFFF8 ],
               couplings = {(0,1):C.GC_37,(0,0):C.GC_117,(0,2):C.GC_69,(0,3):C.GC_69})

V_192 = Vertex(name = 'V_192',
               particles = [ P.t__tilde__, P.b, P.ta__plus__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF16, L.FFFF17, L.FFFF18, L.FFFF3, L.FFFF8, L.FFFF9 ],
               couplings = {(0,3):C.GC_107,(0,2):C.GC_130,(0,0):C.GC_129,(0,1):C.GC_129,(0,4):C.GC_69,(0,5):C.GC_127})

V_193 = Vertex(name = 'V_193',
               particles = [ P.t__tilde__, P.b, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS14, L.FFS7 ],
               couplings = {(0,0):C.GC_686,(0,1):C.GC_415})

V_194 = Vertex(name = 'V_194',
               particles = [ P.t__tilde__, P.b, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS14 ],
               couplings = {(0,0):C.GC_688})

V_195 = Vertex(name = 'V_195',
               particles = [ P.t__tilde__, P.t, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS13, L.FFS19, L.FFS7 ],
               couplings = {(0,0):C.GC_685,(0,2):C.GC_546,(0,1):C.GC_428})

V_196 = Vertex(name = 'V_196',
               particles = [ P.t__tilde__, P.t, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS13 ],
               couplings = {(0,0):C.GC_690})

V_197 = Vertex(name = 'V_197',
               particles = [ P.b__tilde__, P.t, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS4, L.FFS7 ],
               couplings = {(0,0):C.GC_687,(0,1):C.GC_415})

V_198 = Vertex(name = 'V_198',
               particles = [ P.b__tilde__, P.t, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS4 ],
               couplings = {(0,0):C.GC_689})

V_199 = Vertex(name = 'V_199',
               particles = [ P.a, P.W__minus__, P.G0, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_310})

V_200 = Vertex(name = 'V_200',
               particles = [ P.a, P.W__minus__, P.G__minus__, P.G__plus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_309})

V_201 = Vertex(name = 'V_201',
               particles = [ P.a, P.W__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_510})

V_202 = Vertex(name = 'V_202',
               particles = [ P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSS2, L.VSS3 ],
               couplings = {(0,0):C.GC_240,(0,1):C.GC_533})

V_203 = Vertex(name = 'V_203',
               particles = [ P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_606})

V_204 = Vertex(name = 'V_204',
               particles = [ P.W__minus__, P.G0, P.G0, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS5 ],
               couplings = {(0,0):C.GC_258})

V_205 = Vertex(name = 'V_205',
               particles = [ P.W__minus__, P.G0, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS10 ],
               couplings = {(0,0):C.GC_261})

V_206 = Vertex(name = 'V_206',
               particles = [ P.W__minus__, P.G0, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSS3 ],
               couplings = {(0,0):C.GC_480})

V_207 = Vertex(name = 'V_207',
               particles = [ P.W__minus__, P.G0, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSS5 ],
               couplings = {(0,0):C.GC_482})

V_208 = Vertex(name = 'V_208',
               particles = [ P.W__minus__, P.G0, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSSS11 ],
               couplings = {(0,0):C.GC_261})

V_209 = Vertex(name = 'V_209',
               particles = [ P.W__minus__, P.G__minus__, P.G__plus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS3 ],
               couplings = {(0,0):C.GC_262})

V_210 = Vertex(name = 'V_210',
               particles = [ P.W__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSS6 ],
               couplings = {(0,0):C.GC_483})

V_211 = Vertex(name = 'V_211',
               particles = [ P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_238})

V_212 = Vertex(name = 'V_212',
               particles = [ P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_612})

V_213 = Vertex(name = 'V_213',
               particles = [ P.W__minus__, P.W__minus__, P.G0, P.G0, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_231})

V_214 = Vertex(name = 'V_214',
               particles = [ P.W__minus__, P.W__minus__, P.G0, P.G__plus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_232})

V_215 = Vertex(name = 'V_215',
               particles = [ P.W__minus__, P.W__minus__, P.G__plus__, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_230})

V_216 = Vertex(name = 'V_216',
               particles = [ P.W__minus__, P.W__minus__, P.G0, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_468})

V_217 = Vertex(name = 'V_217',
               particles = [ P.W__minus__, P.W__minus__, P.G__plus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_467})

V_218 = Vertex(name = 'V_218',
               particles = [ P.W__minus__, P.W__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSS4 ],
               couplings = {(0,0):C.GC_531})

V_219 = Vertex(name = 'V_219',
               particles = [ P.a, P.a, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS2 ],
               couplings = {(0,0):C.GC_187})

V_220 = Vertex(name = 'V_220',
               particles = [ P.a, P.a, P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS2 ],
               couplings = {(0,0):C.GC_188})

V_221 = Vertex(name = 'V_221',
               particles = [ P.a, P.a, P.W__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVS6 ],
               couplings = {(0,0):C.GC_455})

V_222 = Vertex(name = 'V_222',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVVSS1, L.VVVSS5 ],
               couplings = {(0,1):C.GC_179,(0,0):C.GC_299})

V_223 = Vertex(name = 'V_223',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS1, L.VVVSS5 ],
               couplings = {(0,1):C.GC_179,(0,0):C.GC_298})

V_224 = Vertex(name = 'V_224',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS1, L.VVVSS5 ],
               couplings = {(0,1):C.GC_179,(0,0):C.GC_299})

V_225 = Vertex(name = 'V_225',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVS5, L.VVVS9 ],
               couplings = {(0,1):C.GC_451,(0,0):C.GC_506})

V_226 = Vertex(name = 'V_226',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV15, L.VVVV7 ],
               couplings = {(0,0):C.GC_340,(0,1):C.GC_5})

V_227 = Vertex(name = 'V_227',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV7 ],
               couplings = {(0,0):C.GC_633})

V_228 = Vertex(name = 'V_228',
               particles = [ P.a, P.W__minus__, P.Z, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS1, L.VVVSS3 ],
               couplings = {(0,0):C.GC_305,(0,1):C.GC_323})

V_229 = Vertex(name = 'V_229',
               particles = [ P.a, P.W__minus__, P.Z, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS1, L.VVVSS3 ],
               couplings = {(0,0):C.GC_306,(0,1):C.GC_324})

V_230 = Vertex(name = 'V_230',
               particles = [ P.a, P.W__minus__, P.Z, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVS5, L.VVVS7 ],
               couplings = {(0,0):C.GC_509,(0,1):C.GC_515})

V_231 = Vertex(name = 'V_231',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV10, L.VVVV14 ],
               couplings = {(0,1):C.GC_189,(0,0):C.GC_249})

V_232 = Vertex(name = 'V_232',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV10 ],
               couplings = {(0,0):C.GC_632})

V_233 = Vertex(name = 'V_233',
               particles = [ P.a, P.W__plus__, P.G0, P.G__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_310})

V_234 = Vertex(name = 'V_234',
               particles = [ P.a, P.W__plus__, P.G__minus__, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_311})

V_235 = Vertex(name = 'V_235',
               particles = [ P.a, P.W__plus__, P.G__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_511})

V_236 = Vertex(name = 'V_236',
               particles = [ P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VSS2, L.VSS3 ],
               couplings = {(0,0):C.GC_239,(0,1):C.GC_532})

V_237 = Vertex(name = 'V_237',
               particles = [ P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_605})

V_238 = Vertex(name = 'V_238',
               particles = [ P.W__plus__, P.G0, P.G0, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS5 ],
               couplings = {(0,0):C.GC_258})

V_239 = Vertex(name = 'V_239',
               particles = [ P.W__plus__, P.G0, P.G__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS10 ],
               couplings = {(0,0):C.GC_260})

V_240 = Vertex(name = 'V_240',
               particles = [ P.W__plus__, P.G0, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSS3 ],
               couplings = {(0,0):C.GC_480})

V_241 = Vertex(name = 'V_241',
               particles = [ P.W__plus__, P.G0, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSS5 ],
               couplings = {(0,0):C.GC_481})

V_242 = Vertex(name = 'V_242',
               particles = [ P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_238})

V_243 = Vertex(name = 'V_243',
               particles = [ P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_612})

V_244 = Vertex(name = 'V_244',
               particles = [ P.W__plus__, P.G0, P.G__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSSS6 ],
               couplings = {(0,0):C.GC_259})

V_245 = Vertex(name = 'V_245',
               particles = [ P.W__plus__, P.G__minus__, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS2 ],
               couplings = {(0,0):C.GC_258})

V_246 = Vertex(name = 'V_246',
               particles = [ P.W__plus__, P.G__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSS4 ],
               couplings = {(0,0):C.GC_480})

V_247 = Vertex(name = 'V_247',
               particles = [ P.W__minus__, P.W__plus__, P.G0, P.G0, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_229})

V_248 = Vertex(name = 'V_248',
               particles = [ P.W__minus__, P.W__plus__, P.G__minus__, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_229})

V_249 = Vertex(name = 'V_249',
               particles = [ P.W__minus__, P.W__plus__, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_466})

V_250 = Vertex(name = 'V_250',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_191})

V_251 = Vertex(name = 'V_251',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_191})

V_252 = Vertex(name = 'V_252',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_191})

V_253 = Vertex(name = 'V_253',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS6 ],
               couplings = {(0,0):C.GC_456})

V_254 = Vertex(name = 'V_254',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVVSS3, L.VVVSS5 ],
               couplings = {(0,1):C.GC_297,(0,0):C.GC_181})

V_255 = Vertex(name = 'V_255',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS3, L.VVVSS5 ],
               couplings = {(0,1):C.GC_297,(0,0):C.GC_180})

V_256 = Vertex(name = 'V_256',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS3, L.VVVSS5 ],
               couplings = {(0,1):C.GC_297,(0,0):C.GC_181})

V_257 = Vertex(name = 'V_257',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVS7, L.VVVS9 ],
               couplings = {(0,1):C.GC_505,(0,0):C.GC_452})

V_258 = Vertex(name = 'V_258',
               particles = [ P.a, P.a, P.a, P.W__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVVV3 ],
               couplings = {(0,0):C.GC_341})

V_259 = Vertex(name = 'V_259',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV15, L.VVVV7 ],
               couplings = {(0,0):C.GC_307,(0,1):C.GC_226})

V_260 = Vertex(name = 'V_260',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVV7 ],
               couplings = {(0,0):C.GC_601})

V_261 = Vertex(name = 'V_261',
               particles = [ P.a, P.a, P.W__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVVV19 ],
               couplings = {(0,0):C.GC_192})

V_262 = Vertex(name = 'V_262',
               particles = [ P.W__plus__, P.W__plus__, P.G0, P.G0, P.G__minus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_231})

V_263 = Vertex(name = 'V_263',
               particles = [ P.W__plus__, P.W__plus__, P.G0, P.G__minus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_228})

V_264 = Vertex(name = 'V_264',
               particles = [ P.W__plus__, P.W__plus__, P.G__minus__, P.G__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_230})

V_265 = Vertex(name = 'V_265',
               particles = [ P.W__plus__, P.W__plus__, P.G0, P.G__minus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_465})

V_266 = Vertex(name = 'V_266',
               particles = [ P.W__plus__, P.W__plus__, P.G__minus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_467})

V_267 = Vertex(name = 'V_267',
               particles = [ P.W__plus__, P.W__plus__, P.G__minus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSS4 ],
               couplings = {(0,0):C.GC_531})

V_268 = Vertex(name = 'V_268',
               particles = [ P.a, P.a, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS2 ],
               couplings = {(0,0):C.GC_186})

V_269 = Vertex(name = 'V_269',
               particles = [ P.a, P.a, P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS2 ],
               couplings = {(0,0):C.GC_188})

V_270 = Vertex(name = 'V_270',
               particles = [ P.a, P.a, P.W__plus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVS6 ],
               couplings = {(0,0):C.GC_455})

V_271 = Vertex(name = 'V_271',
               particles = [ P.a, P.W__plus__, P.Z, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS1, L.VVVSS3 ],
               couplings = {(0,0):C.GC_304,(0,1):C.GC_322})

V_272 = Vertex(name = 'V_272',
               particles = [ P.a, P.W__plus__, P.Z, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS1, L.VVVSS3 ],
               couplings = {(0,0):C.GC_306,(0,1):C.GC_324})

V_273 = Vertex(name = 'V_273',
               particles = [ P.a, P.W__plus__, P.Z, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVS5, L.VVVS7 ],
               couplings = {(0,0):C.GC_509,(0,1):C.GC_515})

V_274 = Vertex(name = 'V_274',
               particles = [ P.a, P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVVV21 ],
               couplings = {(0,0):C.GC_313})

V_275 = Vertex(name = 'V_275',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_233})

V_276 = Vertex(name = 'V_276',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_233})

V_277 = Vertex(name = 'V_277',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_233})

V_278 = Vertex(name = 'V_278',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS6 ],
               couplings = {(0,0):C.GC_469})

V_279 = Vertex(name = 'V_279',
               particles = [ P.a, P.a, P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVVVV8 ],
               couplings = {(0,0):C.GC_315})

V_280 = Vertex(name = 'V_280',
               particles = [ P.a, P.Z, P.G0, P.G0, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_391})

V_281 = Vertex(name = 'V_281',
               particles = [ P.a, P.Z, P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_392})

V_282 = Vertex(name = 'V_282',
               particles = [ P.a, P.Z, P.G__minus__, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_391})

V_283 = Vertex(name = 'V_283',
               particles = [ P.a, P.Z, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_586})

V_284 = Vertex(name = 'V_284',
               particles = [ P.Z, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSS1, L.VSS2 ],
               couplings = {(0,1):C.GC_344,(0,0):C.GC_643})

V_285 = Vertex(name = 'V_285',
               particles = [ P.Z, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_640})

V_286 = Vertex(name = 'V_286',
               particles = [ P.Z, P.G0, P.G0, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS4 ],
               couplings = {(0,0):C.GC_359})

V_287 = Vertex(name = 'V_287',
               particles = [ P.Z, P.G0, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS7 ],
               couplings = {(0,0):C.GC_358})

V_288 = Vertex(name = 'V_288',
               particles = [ P.Z, P.G0, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS12 ],
               couplings = {(0,0):C.GC_360})

V_289 = Vertex(name = 'V_289',
               particles = [ P.Z, P.G0, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VSSS8 ],
               couplings = {(0,0):C.GC_569})

V_290 = Vertex(name = 'V_290',
               particles = [ P.Z, P.G0, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSS1 ],
               couplings = {(0,0):C.GC_568})

V_291 = Vertex(name = 'V_291',
               particles = [ P.Z, P.G0, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSS7 ],
               couplings = {(0,0):C.GC_570})

V_292 = Vertex(name = 'V_292',
               particles = [ P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_343})

V_293 = Vertex(name = 'V_293',
               particles = [ P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSS2 ],
               couplings = {(0,0):C.GC_642})

V_294 = Vertex(name = 'V_294',
               particles = [ P.Z, P.G0, P.G0, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSSS9 ],
               couplings = {(0,0):C.GC_356})

V_295 = Vertex(name = 'V_295',
               particles = [ P.Z, P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSSSS8 ],
               couplings = {(0,0):C.GC_357})

V_296 = Vertex(name = 'V_296',
               particles = [ P.Z, P.G__minus__, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSSS1 ],
               couplings = {(0,0):C.GC_356})

V_297 = Vertex(name = 'V_297',
               particles = [ P.Z, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VSSS2 ],
               couplings = {(0,0):C.GC_567})

V_298 = Vertex(name = 'V_298',
               particles = [ P.W__minus__, P.Z, P.G0, P.G0, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_220})

V_299 = Vertex(name = 'V_299',
               particles = [ P.W__minus__, P.Z, P.G0, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_219})

V_300 = Vertex(name = 'V_300',
               particles = [ P.W__minus__, P.Z, P.G0, P.G0, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_217})

V_301 = Vertex(name = 'V_301',
               particles = [ P.W__minus__, P.Z, P.G__minus__, P.G__plus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_222})

V_302 = Vertex(name = 'V_302',
               particles = [ P.W__minus__, P.Z, P.G0, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_218})

V_303 = Vertex(name = 'V_303',
               particles = [ P.W__minus__, P.Z, P.G__plus__, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_215})

V_304 = Vertex(name = 'V_304',
               particles = [ P.W__minus__, P.Z, P.G0, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_559})

V_305 = Vertex(name = 'V_305',
               particles = [ P.W__minus__, P.Z, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_562})

V_306 = Vertex(name = 'V_306',
               particles = [ P.W__minus__, P.Z, P.G0, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_560})

V_307 = Vertex(name = 'V_307',
               particles = [ P.W__minus__, P.Z, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_557})

V_308 = Vertex(name = 'V_308',
               particles = [ P.W__minus__, P.Z, P.Z, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS4 ],
               couplings = {(0,0):C.GC_187})

V_309 = Vertex(name = 'V_309',
               particles = [ P.W__minus__, P.Z, P.Z, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS4 ],
               couplings = {(0,0):C.GC_188})

V_310 = Vertex(name = 'V_310',
               particles = [ P.W__minus__, P.Z, P.Z, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVS8 ],
               couplings = {(0,0):C.GC_455})

V_311 = Vertex(name = 'V_311',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV15, L.VVVV7 ],
               couplings = {(0,0):C.GC_308,(0,1):C.GC_227})

V_312 = Vertex(name = 'V_312',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV7 ],
               couplings = {(0,0):C.GC_634})

V_313 = Vertex(name = 'V_313',
               particles = [ P.W__plus__, P.Z, P.G0, P.G0, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_220})

V_314 = Vertex(name = 'V_314',
               particles = [ P.W__plus__, P.Z, P.G0, P.G__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_219})

V_315 = Vertex(name = 'V_315',
               particles = [ P.W__plus__, P.Z, P.G0, P.G0, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_221})

V_316 = Vertex(name = 'V_316',
               particles = [ P.W__plus__, P.Z, P.G__minus__, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_216})

V_317 = Vertex(name = 'V_317',
               particles = [ P.W__plus__, P.Z, P.G0, P.G__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_218})

V_318 = Vertex(name = 'V_318',
               particles = [ P.W__plus__, P.Z, P.G__minus__, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_223})

V_319 = Vertex(name = 'V_319',
               particles = [ P.W__plus__, P.Z, P.G0, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_561})

V_320 = Vertex(name = 'V_320',
               particles = [ P.W__plus__, P.Z, P.G__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_558})

V_321 = Vertex(name = 'V_321',
               particles = [ P.W__plus__, P.Z, P.G0, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_560})

V_322 = Vertex(name = 'V_322',
               particles = [ P.W__plus__, P.Z, P.G__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_563})

V_323 = Vertex(name = 'V_323',
               particles = [ P.W__plus__, P.Z, P.Z, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVSS4 ],
               couplings = {(0,0):C.GC_186})

V_324 = Vertex(name = 'V_324',
               particles = [ P.W__plus__, P.Z, P.Z, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVSS4 ],
               couplings = {(0,0):C.GC_188})

V_325 = Vertex(name = 'V_325',
               particles = [ P.W__plus__, P.Z, P.Z, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVS8 ],
               couplings = {(0,0):C.GC_455})

V_326 = Vertex(name = 'V_326',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS5 ],
               couplings = {(0,0):C.GC_312})

V_327 = Vertex(name = 'V_327',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS5 ],
               couplings = {(0,0):C.GC_312})

V_328 = Vertex(name = 'V_328',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS5 ],
               couplings = {(0,0):C.GC_312})

V_329 = Vertex(name = 'V_329',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS9 ],
               couplings = {(0,0):C.GC_512})

V_330 = Vertex(name = 'V_330',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVVV9 ],
               couplings = {(0,0):C.GC_314})

V_331 = Vertex(name = 'V_331',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVVV16 ],
               couplings = {(0,0):C.GC_235})

V_332 = Vertex(name = 'V_332',
               particles = [ P.a, P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVVVV17 ],
               couplings = {(0,0):C.GC_237})

V_333 = Vertex(name = 'V_333',
               particles = [ P.Z, P.Z, P.G0, P.G0, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_404})

V_334 = Vertex(name = 'V_334',
               particles = [ P.Z, P.Z, P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_403})

V_335 = Vertex(name = 'V_335',
               particles = [ P.Z, P.Z, P.G0, P.G0, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_402})

V_336 = Vertex(name = 'V_336',
               particles = [ P.Z, P.Z, P.H, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_404})

V_337 = Vertex(name = 'V_337',
               particles = [ P.Z, P.Z, P.G0, P.G0, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_401})

V_338 = Vertex(name = 'V_338',
               particles = [ P.Z, P.Z, P.G__minus__, P.G__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSSS1 ],
               couplings = {(0,0):C.GC_401})

V_339 = Vertex(name = 'V_339',
               particles = [ P.Z, P.Z, P.G0, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_593})

V_340 = Vertex(name = 'V_340',
               particles = [ P.Z, P.Z, P.H, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_594})

V_341 = Vertex(name = 'V_341',
               particles = [ P.Z, P.Z, P.G__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSSS1 ],
               couplings = {(0,0):C.GC_592})

V_342 = Vertex(name = 'V_342',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_234})

V_343 = Vertex(name = 'V_343',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_234})

V_344 = Vertex(name = 'V_344',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVSS2 ],
               couplings = {(0,0):C.GC_234})

V_345 = Vertex(name = 'V_345',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVVVS6 ],
               couplings = {(0,0):C.GC_470})

V_346 = Vertex(name = 'V_346',
               particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVVVV8 ],
               couplings = {(0,0):C.GC_224})

V_347 = Vertex(name = 'V_347',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVVV20 ],
               couplings = {(0,0):C.GC_236})

V_348 = Vertex(name = 'V_348',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV2, L.FFV7 ],
               couplings = {(0,0):C.GC_1,(0,2):C.GC_595,(0,1):C.GC_530})

V_349 = Vertex(name = 'V_349',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_620})

V_350 = Vertex(name = 'V_350',
               particles = [ P.b__tilde__, P.b, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_50,(0,1):C.GC_171})

V_351 = Vertex(name = 'V_351',
               particles = [ P.b__tilde__, P.b, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_30,(0,1):C.GC_87})

V_352 = Vertex(name = 'V_352',
               particles = [ P.b__tilde__, P.b, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_547,(0,0):C.GC_422})

V_353 = Vertex(name = 'V_353',
               particles = [ P.b__tilde__, P.b, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_28,(0,1):C.GC_86})

V_354 = Vertex(name = 'V_354',
               particles = [ P.t__tilde__, P.b, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_166})

V_355 = Vertex(name = 'V_355',
               particles = [ P.t__tilde__, P.b, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_165})

V_356 = Vertex(name = 'V_356',
               particles = [ P.t__tilde__, P.b, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS30 ],
               couplings = {(0,1):C.GC_350,(0,0):C.GC_446})

V_357 = Vertex(name = 'V_357',
               particles = [ P.t__tilde__, P.b, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_65})

V_358 = Vertex(name = 'V_358',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV13, L.FFV2 ],
               couplings = {(0,0):C.GC_2,(0,1):C.GC_595,(0,2):C.GC_529})

V_359 = Vertex(name = 'V_359',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_621})

V_360 = Vertex(name = 'V_360',
               particles = [ P.c__tilde__, P.c, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_53,(0,1):C.GC_178})

V_361 = Vertex(name = 'V_361',
               particles = [ P.c__tilde__, P.c, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_33,(0,1):C.GC_100})

V_362 = Vertex(name = 'V_362',
               particles = [ P.c__tilde__, P.c, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_548,(0,0):C.GC_430})

V_363 = Vertex(name = 'V_363',
               particles = [ P.c__tilde__, P.c, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_31,(0,1):C.GC_99})

V_364 = Vertex(name = 'V_364',
               particles = [ P.s__tilde__, P.c, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_169})

V_365 = Vertex(name = 'V_365',
               particles = [ P.s__tilde__, P.c, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_170})

V_366 = Vertex(name = 'V_366',
               particles = [ P.s__tilde__, P.c, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_449})

V_367 = Vertex(name = 'V_367',
               particles = [ P.s__tilde__, P.c, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_67})

V_368 = Vertex(name = 'V_368',
               particles = [ P.s__tilde__, P.c, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_66})

V_369 = Vertex(name = 'V_369',
               particles = [ P.s__tilde__, P.c, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_416})

V_370 = Vertex(name = 'V_370',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV2, L.FFV7 ],
               couplings = {(0,0):C.GC_1,(0,2):C.GC_595,(0,1):C.GC_530})

V_371 = Vertex(name = 'V_371',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_620})

V_372 = Vertex(name = 'V_372',
               particles = [ P.d__tilde__, P.d, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_52,(0,1):C.GC_171})

V_373 = Vertex(name = 'V_373',
               particles = [ P.d__tilde__, P.d, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_34,(0,1):C.GC_87})

V_374 = Vertex(name = 'V_374',
               particles = [ P.d__tilde__, P.d, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_549,(0,0):C.GC_422})

V_375 = Vertex(name = 'V_375',
               particles = [ P.d__tilde__, P.d, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_32,(0,1):C.GC_86})

V_376 = Vertex(name = 'V_376',
               particles = [ P.u__tilde__, P.d, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_169})

V_377 = Vertex(name = 'V_377',
               particles = [ P.u__tilde__, P.d, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_168})

V_378 = Vertex(name = 'V_378',
               particles = [ P.u__tilde__, P.d, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_448})

V_379 = Vertex(name = 'V_379',
               particles = [ P.u__tilde__, P.d, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_68})

V_380 = Vertex(name = 'V_380',
               particles = [ P.u__tilde__, P.d, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_66})

V_381 = Vertex(name = 'V_381',
               particles = [ P.u__tilde__, P.d, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_416})

V_382 = Vertex(name = 'V_382',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV12, L.FFV2 ],
               couplings = {(0,0):C.GC_3,(0,1):C.GC_596,(0,2):C.GC_530})

V_383 = Vertex(name = 'V_383',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_623})

V_384 = Vertex(name = 'V_384',
               particles = [ P.e__plus__, P.e__minus__, P.a, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_44,(0,1):C.GC_174})

V_385 = Vertex(name = 'V_385',
               particles = [ P.e__plus__, P.e__minus__, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_18,(0,1):C.GC_91})

V_386 = Vertex(name = 'V_386',
               particles = [ P.e__plus__, P.e__minus__, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_541,(0,0):C.GC_425})

V_387 = Vertex(name = 'V_387',
               particles = [ P.e__plus__, P.e__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_16,(0,1):C.GC_90})

V_388 = Vertex(name = 'V_388',
               particles = [ P.ve__tilde__, P.e__minus__, P.a, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_157})

V_389 = Vertex(name = 'V_389',
               particles = [ P.ve__tilde__, P.e__minus__, P.a, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_156})

V_390 = Vertex(name = 'V_390',
               particles = [ P.ve__tilde__, P.e__minus__, P.a, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_440})

V_391 = Vertex(name = 'V_391',
               particles = [ P.ve__tilde__, P.e__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_56})

V_392 = Vertex(name = 'V_392',
               particles = [ P.ve__tilde__, P.e__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_54})

V_393 = Vertex(name = 'V_393',
               particles = [ P.ve__tilde__, P.e__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_412})

V_394 = Vertex(name = 'V_394',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV12, L.FFV2 ],
               couplings = {(0,0):C.GC_3,(0,1):C.GC_596,(0,2):C.GC_530})

V_395 = Vertex(name = 'V_395',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_623})

V_396 = Vertex(name = 'V_396',
               particles = [ P.mu__plus__, P.mu__minus__, P.a, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_46,(0,1):C.GC_175})

V_397 = Vertex(name = 'V_397',
               particles = [ P.mu__plus__, P.mu__minus__, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_22,(0,1):C.GC_94})

V_398 = Vertex(name = 'V_398',
               particles = [ P.mu__plus__, P.mu__minus__, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_543,(0,0):C.GC_427})

V_399 = Vertex(name = 'V_399',
               particles = [ P.mu__plus__, P.mu__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_20,(0,1):C.GC_93})

V_400 = Vertex(name = 'V_400',
               particles = [ P.vm__tilde__, P.mu__minus__, P.a, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_160})

V_401 = Vertex(name = 'V_401',
               particles = [ P.vm__tilde__, P.mu__minus__, P.a, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_159})

V_402 = Vertex(name = 'V_402',
               particles = [ P.vm__tilde__, P.mu__minus__, P.a, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_442})

V_403 = Vertex(name = 'V_403',
               particles = [ P.vm__tilde__, P.mu__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_59})

V_404 = Vertex(name = 'V_404',
               particles = [ P.vm__tilde__, P.mu__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_57})

V_405 = Vertex(name = 'V_405',
               particles = [ P.vm__tilde__, P.mu__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_413})

V_406 = Vertex(name = 'V_406',
               particles = [ P.c__tilde__, P.s, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_169})

V_407 = Vertex(name = 'V_407',
               particles = [ P.c__tilde__, P.s, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_168})

V_408 = Vertex(name = 'V_408',
               particles = [ P.c__tilde__, P.s, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_448})

V_409 = Vertex(name = 'V_409',
               particles = [ P.c__tilde__, P.s, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_68})

V_410 = Vertex(name = 'V_410',
               particles = [ P.c__tilde__, P.s, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_66})

V_411 = Vertex(name = 'V_411',
               particles = [ P.c__tilde__, P.s, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_416})

V_412 = Vertex(name = 'V_412',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV2, L.FFV7 ],
               couplings = {(0,0):C.GC_1,(0,2):C.GC_595,(0,1):C.GC_530})

V_413 = Vertex(name = 'V_413',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_620})

V_414 = Vertex(name = 'V_414',
               particles = [ P.s__tilde__, P.s, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_52,(0,1):C.GC_171})

V_415 = Vertex(name = 'V_415',
               particles = [ P.s__tilde__, P.s, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_34,(0,1):C.GC_87})

V_416 = Vertex(name = 'V_416',
               particles = [ P.s__tilde__, P.s, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_549,(0,0):C.GC_422})

V_417 = Vertex(name = 'V_417',
               particles = [ P.s__tilde__, P.s, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_32,(0,1):C.GC_86})

V_418 = Vertex(name = 'V_418',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV12, L.FFV2 ],
               couplings = {(0,0):C.GC_3,(0,1):C.GC_596,(0,2):C.GC_530})

V_419 = Vertex(name = 'V_419',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_623})

V_420 = Vertex(name = 'V_420',
               particles = [ P.ta__plus__, P.ta__minus__, P.a, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_48,(0,1):C.GC_177})

V_421 = Vertex(name = 'V_421',
               particles = [ P.ta__plus__, P.ta__minus__, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_26,(0,1):C.GC_98})

V_422 = Vertex(name = 'V_422',
               particles = [ P.ta__plus__, P.ta__minus__, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_545,(0,0):C.GC_429})

V_423 = Vertex(name = 'V_423',
               particles = [ P.ta__plus__, P.ta__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_24,(0,1):C.GC_97})

V_424 = Vertex(name = 'V_424',
               particles = [ P.vt__tilde__, P.ta__minus__, P.a, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_163})

V_425 = Vertex(name = 'V_425',
               particles = [ P.vt__tilde__, P.ta__minus__, P.a, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_162})

V_426 = Vertex(name = 'V_426',
               particles = [ P.vt__tilde__, P.ta__minus__, P.a, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_444})

V_427 = Vertex(name = 'V_427',
               particles = [ P.vt__tilde__, P.ta__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_62})

V_428 = Vertex(name = 'V_428',
               particles = [ P.vt__tilde__, P.ta__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_60})

V_429 = Vertex(name = 'V_429',
               particles = [ P.vt__tilde__, P.ta__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_414})

V_430 = Vertex(name = 'V_430',
               particles = [ P.b__tilde__, P.t, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_166})

V_431 = Vertex(name = 'V_431',
               particles = [ P.b__tilde__, P.t, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_167})

V_432 = Vertex(name = 'V_432',
               particles = [ P.b__tilde__, P.t, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS10, L.FFVS3 ],
               couplings = {(0,0):C.GC_351,(0,1):C.GC_447})

V_433 = Vertex(name = 'V_433',
               particles = [ P.b__tilde__, P.t, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_64})

V_434 = Vertex(name = 'V_434',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV13, L.FFV2, L.FFV26 ],
               couplings = {(0,0):C.GC_2,(0,1):C.GC_595,(0,2):C.GC_529,(0,3):C.GC_565})

V_435 = Vertex(name = 'V_435',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_621})

V_436 = Vertex(name = 'V_436',
               particles = [ P.t__tilde__, P.t, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_51,(0,1):C.GC_176})

V_437 = Vertex(name = 'V_437',
               particles = [ P.d__tilde__, P.u, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_169})

V_438 = Vertex(name = 'V_438',
               particles = [ P.d__tilde__, P.u, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_170})

V_439 = Vertex(name = 'V_439',
               particles = [ P.d__tilde__, P.u, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_449})

V_440 = Vertex(name = 'V_440',
               particles = [ P.d__tilde__, P.u, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_67})

V_441 = Vertex(name = 'V_441',
               particles = [ P.d__tilde__, P.u, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_66})

V_442 = Vertex(name = 'V_442',
               particles = [ P.d__tilde__, P.u, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_416})

V_443 = Vertex(name = 'V_443',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV13, L.FFV2 ],
               couplings = {(0,0):C.GC_2,(0,1):C.GC_595,(0,2):C.GC_529})

V_444 = Vertex(name = 'V_444',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10 ],
               couplings = {(0,0):C.GC_621})

V_445 = Vertex(name = 'V_445',
               particles = [ P.u__tilde__, P.u, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_53,(0,1):C.GC_178})

V_446 = Vertex(name = 'V_446',
               particles = [ P.u__tilde__, P.u, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_33,(0,1):C.GC_100})

V_447 = Vertex(name = 'V_447',
               particles = [ P.u__tilde__, P.u, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS19, L.FFS7 ],
               couplings = {(0,1):C.GC_548,(0,0):C.GC_430})

V_448 = Vertex(name = 'V_448',
               particles = [ P.u__tilde__, P.u, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS3, L.FFSS7 ],
               couplings = {(0,0):C.GC_31,(0,1):C.GC_99})

V_449 = Vertex(name = 'V_449',
               particles = [ P.e__plus__, P.ve, P.a, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_157})

V_450 = Vertex(name = 'V_450',
               particles = [ P.e__plus__, P.ve, P.a, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_158})

V_451 = Vertex(name = 'V_451',
               particles = [ P.e__plus__, P.ve, P.a, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_441})

V_452 = Vertex(name = 'V_452',
               particles = [ P.e__plus__, P.ve, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_55})

V_453 = Vertex(name = 'V_453',
               particles = [ P.e__plus__, P.ve, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_54})

V_454 = Vertex(name = 'V_454',
               particles = [ P.e__plus__, P.ve, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_412})

V_455 = Vertex(name = 'V_455',
               particles = [ P.ve__tilde__, P.ve, P.a, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_45})

V_456 = Vertex(name = 'V_456',
               particles = [ P.ve__tilde__, P.ve, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_597})

V_457 = Vertex(name = 'V_457',
               particles = [ P.ve__tilde__, P.ve, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_17})

V_458 = Vertex(name = 'V_458',
               particles = [ P.ve__tilde__, P.ve, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_540})

V_459 = Vertex(name = 'V_459',
               particles = [ P.ve__tilde__, P.ve, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_15})

V_460 = Vertex(name = 'V_460',
               particles = [ P.mu__plus__, P.vm, P.a, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_160})

V_461 = Vertex(name = 'V_461',
               particles = [ P.mu__plus__, P.vm, P.a, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_161})

V_462 = Vertex(name = 'V_462',
               particles = [ P.mu__plus__, P.vm, P.a, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_443})

V_463 = Vertex(name = 'V_463',
               particles = [ P.mu__plus__, P.vm, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_58})

V_464 = Vertex(name = 'V_464',
               particles = [ P.mu__plus__, P.vm, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_57})

V_465 = Vertex(name = 'V_465',
               particles = [ P.mu__plus__, P.vm, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_413})

V_466 = Vertex(name = 'V_466',
               particles = [ P.vm__tilde__, P.vm, P.a, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_47})

V_467 = Vertex(name = 'V_467',
               particles = [ P.vm__tilde__, P.vm, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_597})

V_468 = Vertex(name = 'V_468',
               particles = [ P.vm__tilde__, P.vm, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_21})

V_469 = Vertex(name = 'V_469',
               particles = [ P.vm__tilde__, P.vm, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_542})

V_470 = Vertex(name = 'V_470',
               particles = [ P.vm__tilde__, P.vm, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_19})

V_471 = Vertex(name = 'V_471',
               particles = [ P.ta__plus__, P.vt, P.a, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_163})

V_472 = Vertex(name = 'V_472',
               particles = [ P.ta__plus__, P.vt, P.a, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_164})

V_473 = Vertex(name = 'V_473',
               particles = [ P.ta__plus__, P.vt, P.a, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_445})

V_474 = Vertex(name = 'V_474',
               particles = [ P.ta__plus__, P.vt, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_61})

V_475 = Vertex(name = 'V_475',
               particles = [ P.ta__plus__, P.vt, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_60})

V_476 = Vertex(name = 'V_476',
               particles = [ P.ta__plus__, P.vt, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_414})

V_477 = Vertex(name = 'V_477',
               particles = [ P.vt__tilde__, P.vt, P.a, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_49})

V_478 = Vertex(name = 'V_478',
               particles = [ P.vt__tilde__, P.vt, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_597})

V_479 = Vertex(name = 'V_479',
               particles = [ P.vt__tilde__, P.vt, P.G0, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_25})

V_480 = Vertex(name = 'V_480',
               particles = [ P.vt__tilde__, P.vt, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFS7 ],
               couplings = {(0,0):C.GC_544})

V_481 = Vertex(name = 'V_481',
               particles = [ P.vt__tilde__, P.vt, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFSS3 ],
               couplings = {(0,0):C.GC_23})

V_482 = Vertex(name = 'V_482',
               particles = [ P.b__tilde__, P.b, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_483 = Vertex(name = 'V_483',
               particles = [ P.c__tilde__, P.c, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_484 = Vertex(name = 'V_484',
               particles = [ P.d__tilde__, P.d, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_485 = Vertex(name = 'V_485',
               particles = [ P.s__tilde__, P.s, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_486 = Vertex(name = 'V_486',
               particles = [ P.t__tilde__, P.t, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV10, L.FFV26 ],
               couplings = {(0,0):C.GC_11,(0,1):C.GC_458})

V_487 = Vertex(name = 'V_487',
               particles = [ P.u__tilde__, P.u, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_488 = Vertex(name = 'V_488',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_279,(0,1):C.GC_256})

V_489 = Vertex(name = 'V_489',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_278,(0,1):C.GC_255})

V_490 = Vertex(name = 'V_490',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_494,(0,0):C.GC_478})

V_491 = Vertex(name = 'V_491',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_291})

V_492 = Vertex(name = 'V_492',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_281,(0,1):C.GC_290})

V_493 = Vertex(name = 'V_493',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_496,(0,0):C.GC_502})

V_494 = Vertex(name = 'V_494',
               particles = [ P.s__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_495 = Vertex(name = 'V_495',
               particles = [ P.s__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_611})

V_496 = Vertex(name = 'V_496',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_497 = Vertex(name = 'V_497',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_498 = Vertex(name = 'V_498',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_499 = Vertex(name = 'V_499',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_477})

V_500 = Vertex(name = 'V_500',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_256})

V_501 = Vertex(name = 'V_501',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_281,(0,1):C.GC_255})

V_502 = Vertex(name = 'V_502',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_496,(0,0):C.GC_478})

V_503 = Vertex(name = 'V_503',
               particles = [ P.e__plus__, P.e__minus__, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_267,(0,1):C.GC_264})

V_504 = Vertex(name = 'V_504',
               particles = [ P.e__plus__, P.e__minus__, P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_266,(0,1):C.GC_263})

V_505 = Vertex(name = 'V_505',
               particles = [ P.e__plus__, P.e__minus__, P.W__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_486,(0,0):C.GC_484})

V_506 = Vertex(name = 'V_506',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_270,(0,1):C.GC_276})

V_507 = Vertex(name = 'V_507',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_269,(0,1):C.GC_275})

V_508 = Vertex(name = 'V_508',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_488,(0,0):C.GC_492})

V_509 = Vertex(name = 'V_509',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_256})

V_510 = Vertex(name = 'V_510',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_281,(0,1):C.GC_255})

V_511 = Vertex(name = 'V_511',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_496,(0,0):C.GC_478})

V_512 = Vertex(name = 'V_512',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_273,(0,1):C.GC_288})

V_513 = Vertex(name = 'V_513',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_272,(0,1):C.GC_287})

V_514 = Vertex(name = 'V_514',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_490,(0,0):C.GC_500})

V_515 = Vertex(name = 'V_515',
               particles = [ P.b__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV5 ],
               couplings = {(0,1):C.GC_437,(0,0):C.GC_241})

V_516 = Vertex(name = 'V_516',
               particles = [ P.b__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_610})

V_517 = Vertex(name = 'V_517',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_253})

V_518 = Vertex(name = 'V_518',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_253})

V_519 = Vertex(name = 'V_519',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_253})

V_520 = Vertex(name = 'V_520',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS10, L.FFVS3 ],
               couplings = {(0,0):C.GC_148,(0,1):C.GC_476})

V_521 = Vertex(name = 'V_521',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_279,(0,1):C.GC_285})

V_522 = Vertex(name = 'V_522',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_278,(0,1):C.GC_284})

V_523 = Vertex(name = 'V_523',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3, L.FFVS30 ],
               couplings = {(0,1):C.GC_494,(0,0):C.GC_498,(0,2):C.GC_151})

V_524 = Vertex(name = 'V_524',
               particles = [ P.d__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_525 = Vertex(name = 'V_525',
               particles = [ P.d__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_611})

V_526 = Vertex(name = 'V_526',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_527 = Vertex(name = 'V_527',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_528 = Vertex(name = 'V_528',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_529 = Vertex(name = 'V_529',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_477})

V_530 = Vertex(name = 'V_530',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_291})

V_531 = Vertex(name = 'V_531',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_281,(0,1):C.GC_290})

V_532 = Vertex(name = 'V_532',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_496,(0,0):C.GC_502})

V_533 = Vertex(name = 'V_533',
               particles = [ P.e__plus__, P.ve, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_534 = Vertex(name = 'V_534',
               particles = [ P.e__plus__, P.ve, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_607})

V_535 = Vertex(name = 'V_535',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_250})

V_536 = Vertex(name = 'V_536',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_250})

V_537 = Vertex(name = 'V_537',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_250})

V_538 = Vertex(name = 'V_538',
               particles = [ P.e__plus__, P.ve, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_473})

V_539 = Vertex(name = 'V_539',
               particles = [ P.ve__tilde__, P.ve, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_267})

V_540 = Vertex(name = 'V_540',
               particles = [ P.ve__tilde__, P.ve, P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_266})

V_541 = Vertex(name = 'V_541',
               particles = [ P.ve__tilde__, P.ve, P.W__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_486})

V_542 = Vertex(name = 'V_542',
               particles = [ P.mu__plus__, P.vm, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_543 = Vertex(name = 'V_543',
               particles = [ P.mu__plus__, P.vm, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_608})

V_544 = Vertex(name = 'V_544',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_251})

V_545 = Vertex(name = 'V_545',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_251})

V_546 = Vertex(name = 'V_546',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_251})

V_547 = Vertex(name = 'V_547',
               particles = [ P.mu__plus__, P.vm, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_474})

V_548 = Vertex(name = 'V_548',
               particles = [ P.vm__tilde__, P.vm, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_270})

V_549 = Vertex(name = 'V_549',
               particles = [ P.vm__tilde__, P.vm, P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_269})

V_550 = Vertex(name = 'V_550',
               particles = [ P.vm__tilde__, P.vm, P.W__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_488})

V_551 = Vertex(name = 'V_551',
               particles = [ P.ta__plus__, P.vt, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_552 = Vertex(name = 'V_552',
               particles = [ P.ta__plus__, P.vt, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_609})

V_553 = Vertex(name = 'V_553',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_252})

V_554 = Vertex(name = 'V_554',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_252})

V_555 = Vertex(name = 'V_555',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_252})

V_556 = Vertex(name = 'V_556',
               particles = [ P.ta__plus__, P.vt, P.W__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_475})

V_557 = Vertex(name = 'V_557',
               particles = [ P.vt__tilde__, P.vt, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_273})

V_558 = Vertex(name = 'V_558',
               particles = [ P.vt__tilde__, P.vt, P.W__minus__, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_272})

V_559 = Vertex(name = 'V_559',
               particles = [ P.vt__tilde__, P.vt, P.W__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_490})

V_560 = Vertex(name = 'V_560',
               particles = [ P.b__tilde__, P.b, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_279,(0,1):C.GC_256})

V_561 = Vertex(name = 'V_561',
               particles = [ P.b__tilde__, P.b, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_280,(0,1):C.GC_257})

V_562 = Vertex(name = 'V_562',
               particles = [ P.b__tilde__, P.b, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_495,(0,0):C.GC_479})

V_563 = Vertex(name = 'V_563',
               particles = [ P.t__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV25 ],
               couplings = {(0,1):C.GC_437,(0,0):C.GC_241})

V_564 = Vertex(name = 'V_564',
               particles = [ P.t__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_610})

V_565 = Vertex(name = 'V_565',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_253})

V_566 = Vertex(name = 'V_566',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_253})

V_567 = Vertex(name = 'V_567',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_253})

V_568 = Vertex(name = 'V_568',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS30 ],
               couplings = {(0,1):C.GC_148,(0,0):C.GC_476})

V_569 = Vertex(name = 'V_569',
               particles = [ P.c__tilde__, P.c, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_291})

V_570 = Vertex(name = 'V_570',
               particles = [ P.c__tilde__, P.c, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_283,(0,1):C.GC_292})

V_571 = Vertex(name = 'V_571',
               particles = [ P.c__tilde__, P.c, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_497,(0,0):C.GC_503})

V_572 = Vertex(name = 'V_572',
               particles = [ P.d__tilde__, P.d, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_256})

V_573 = Vertex(name = 'V_573',
               particles = [ P.d__tilde__, P.d, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_283,(0,1):C.GC_257})

V_574 = Vertex(name = 'V_574',
               particles = [ P.d__tilde__, P.d, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_497,(0,0):C.GC_479})

V_575 = Vertex(name = 'V_575',
               particles = [ P.u__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_576 = Vertex(name = 'V_576',
               particles = [ P.u__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_611})

V_577 = Vertex(name = 'V_577',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_578 = Vertex(name = 'V_578',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_579 = Vertex(name = 'V_579',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_580 = Vertex(name = 'V_580',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_477})

V_581 = Vertex(name = 'V_581',
               particles = [ P.e__plus__, P.e__minus__, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_267,(0,1):C.GC_264})

V_582 = Vertex(name = 'V_582',
               particles = [ P.e__plus__, P.e__minus__, P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_268,(0,1):C.GC_265})

V_583 = Vertex(name = 'V_583',
               particles = [ P.e__plus__, P.e__minus__, P.W__plus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_487,(0,0):C.GC_485})

V_584 = Vertex(name = 'V_584',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_585 = Vertex(name = 'V_585',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_607})

V_586 = Vertex(name = 'V_586',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_250})

V_587 = Vertex(name = 'V_587',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_250})

V_588 = Vertex(name = 'V_588',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_250})

V_589 = Vertex(name = 'V_589',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_473})

V_590 = Vertex(name = 'V_590',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_270,(0,1):C.GC_276})

V_591 = Vertex(name = 'V_591',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_271,(0,1):C.GC_277})

V_592 = Vertex(name = 'V_592',
               particles = [ P.mu__plus__, P.mu__minus__, P.W__plus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_489,(0,0):C.GC_493})

V_593 = Vertex(name = 'V_593',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_594 = Vertex(name = 'V_594',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_608})

V_595 = Vertex(name = 'V_595',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_251})

V_596 = Vertex(name = 'V_596',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_251})

V_597 = Vertex(name = 'V_597',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_251})

V_598 = Vertex(name = 'V_598',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_474})

V_599 = Vertex(name = 'V_599',
               particles = [ P.c__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_600 = Vertex(name = 'V_600',
               particles = [ P.c__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_611})

V_601 = Vertex(name = 'V_601',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_602 = Vertex(name = 'V_602',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_603 = Vertex(name = 'V_603',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_254})

V_604 = Vertex(name = 'V_604',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_477})

V_605 = Vertex(name = 'V_605',
               particles = [ P.s__tilde__, P.s, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_256})

V_606 = Vertex(name = 'V_606',
               particles = [ P.s__tilde__, P.s, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_283,(0,1):C.GC_257})

V_607 = Vertex(name = 'V_607',
               particles = [ P.s__tilde__, P.s, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_497,(0,0):C.GC_479})

V_608 = Vertex(name = 'V_608',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_273,(0,1):C.GC_288})

V_609 = Vertex(name = 'V_609',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_274,(0,1):C.GC_289})

V_610 = Vertex(name = 'V_610',
               particles = [ P.ta__plus__, P.ta__minus__, P.W__plus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_491,(0,0):C.GC_501})

V_611 = Vertex(name = 'V_611',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_241})

V_612 = Vertex(name = 'V_612',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_609})

V_613 = Vertex(name = 'V_613',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_252})

V_614 = Vertex(name = 'V_614',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_252})

V_615 = Vertex(name = 'V_615',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_252})

V_616 = Vertex(name = 'V_616',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_475})

V_617 = Vertex(name = 'V_617',
               particles = [ P.t__tilde__, P.t, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_279,(0,1):C.GC_285})

V_618 = Vertex(name = 'V_618',
               particles = [ P.t__tilde__, P.t, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_280,(0,1):C.GC_286})

V_619 = Vertex(name = 'V_619',
               particles = [ P.t__tilde__, P.t, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS10, L.FFVS15, L.FFVS3 ],
               couplings = {(0,2):C.GC_495,(0,1):C.GC_499,(0,0):C.GC_150})

V_620 = Vertex(name = 'V_620',
               particles = [ P.u__tilde__, P.u, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_282,(0,1):C.GC_291})

V_621 = Vertex(name = 'V_621',
               particles = [ P.u__tilde__, P.u, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_283,(0,1):C.GC_292})

V_622 = Vertex(name = 'V_622',
               particles = [ P.u__tilde__, P.u, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_497,(0,0):C.GC_503})

V_623 = Vertex(name = 'V_623',
               particles = [ P.ve__tilde__, P.ve, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_267})

V_624 = Vertex(name = 'V_624',
               particles = [ P.ve__tilde__, P.ve, P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_268})

V_625 = Vertex(name = 'V_625',
               particles = [ P.ve__tilde__, P.ve, P.W__plus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_487})

V_626 = Vertex(name = 'V_626',
               particles = [ P.vm__tilde__, P.vm, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_270})

V_627 = Vertex(name = 'V_627',
               particles = [ P.vm__tilde__, P.vm, P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_271})

V_628 = Vertex(name = 'V_628',
               particles = [ P.vm__tilde__, P.vm, P.W__plus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_489})

V_629 = Vertex(name = 'V_629',
               particles = [ P.vt__tilde__, P.vt, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_273})

V_630 = Vertex(name = 'V_630',
               particles = [ P.vt__tilde__, P.vt, P.W__plus__, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_274})

V_631 = Vertex(name = 'V_631',
               particles = [ P.vt__tilde__, P.vt, P.W__plus__, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_491})

V_632 = Vertex(name = 'V_632',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV2, L.FFV6, L.FFV7 ],
               couplings = {(0,1):C.GC_243,(0,3):C.GC_316,(0,0):C.GC_613,(0,2):C.GC_641})

V_633 = Vertex(name = 'V_633',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV7 ],
               couplings = {(0,0):C.GC_653,(0,1):C.GC_636})

V_634 = Vertex(name = 'V_634',
               particles = [ P.b__tilde__, P.b, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_380,(0,1):C.GC_355})

V_635 = Vertex(name = 'V_635',
               particles = [ P.b__tilde__, P.b, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_377,(0,1):C.GC_354})

V_636 = Vertex(name = 'V_636',
               particles = [ P.b__tilde__, P.b, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_380,(0,1):C.GC_355})

V_637 = Vertex(name = 'V_637',
               particles = [ P.b__tilde__, P.b, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_580,(0,0):C.GC_566})

V_638 = Vertex(name = 'V_638',
               particles = [ P.t__tilde__, P.b, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_335})

V_639 = Vertex(name = 'V_639',
               particles = [ P.t__tilde__, P.b, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_334})

V_640 = Vertex(name = 'V_640',
               particles = [ P.t__tilde__, P.b, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3, L.FFVS30 ],
               couplings = {(0,1):C.GC_347,(0,0):C.GC_522})

V_641 = Vertex(name = 'V_641',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV13, L.FFV2, L.FFV6 ],
               couplings = {(0,2):C.GC_242,(0,1):C.GC_316,(0,0):C.GC_614,(0,3):C.GC_658})

V_642 = Vertex(name = 'V_642',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV13, L.FFV2 ],
               couplings = {(0,1):C.GC_654,(0,0):C.GC_636})

V_643 = Vertex(name = 'V_643',
               particles = [ P.c__tilde__, P.c, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_382,(0,1):C.GC_390})

V_644 = Vertex(name = 'V_644',
               particles = [ P.c__tilde__, P.c, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_383,(0,1):C.GC_389})

V_645 = Vertex(name = 'V_645',
               particles = [ P.c__tilde__, P.c, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_382,(0,1):C.GC_390})

V_646 = Vertex(name = 'V_646',
               particles = [ P.c__tilde__, P.c, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_581,(0,0):C.GC_585})

V_647 = Vertex(name = 'V_647',
               particles = [ P.s__tilde__, P.c, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_338})

V_648 = Vertex(name = 'V_648',
               particles = [ P.s__tilde__, P.c, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_339})

V_649 = Vertex(name = 'V_649',
               particles = [ P.s__tilde__, P.c, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_525})

V_650 = Vertex(name = 'V_650',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV2, L.FFV6, L.FFV7 ],
               couplings = {(0,1):C.GC_243,(0,3):C.GC_316,(0,0):C.GC_613,(0,2):C.GC_641})

V_651 = Vertex(name = 'V_651',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV7 ],
               couplings = {(0,0):C.GC_655,(0,1):C.GC_636})

V_652 = Vertex(name = 'V_652',
               particles = [ P.d__tilde__, P.d, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_384,(0,1):C.GC_355})

V_653 = Vertex(name = 'V_653',
               particles = [ P.d__tilde__, P.d, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_381,(0,1):C.GC_354})

V_654 = Vertex(name = 'V_654',
               particles = [ P.d__tilde__, P.d, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_384,(0,1):C.GC_355})

V_655 = Vertex(name = 'V_655',
               particles = [ P.d__tilde__, P.d, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_582,(0,0):C.GC_566})

V_656 = Vertex(name = 'V_656',
               particles = [ P.u__tilde__, P.d, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_338})

V_657 = Vertex(name = 'V_657',
               particles = [ P.u__tilde__, P.d, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_337})

V_658 = Vertex(name = 'V_658',
               particles = [ P.u__tilde__, P.d, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_524})

V_659 = Vertex(name = 'V_659',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV12, L.FFV2, L.FFV6, L.FFV9 ],
               couplings = {(0,2):C.GC_243,(0,1):C.GC_317,(0,0):C.GC_617,(0,4):C.GC_535,(0,3):C.GC_644})

V_660 = Vertex(name = 'V_660',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV12, L.FFV2 ],
               couplings = {(0,1):C.GC_645,(0,0):C.GC_638})

V_661 = Vertex(name = 'V_661',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_366,(0,1):C.GC_362})

V_662 = Vertex(name = 'V_662',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_363,(0,1):C.GC_361})

V_663 = Vertex(name = 'V_663',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_366,(0,1):C.GC_362})

V_664 = Vertex(name = 'V_664',
               particles = [ P.e__plus__, P.e__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_573,(0,0):C.GC_571})

V_665 = Vertex(name = 'V_665',
               particles = [ P.ve__tilde__, P.e__minus__, P.Z, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_326})

V_666 = Vertex(name = 'V_666',
               particles = [ P.ve__tilde__, P.e__minus__, P.Z, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_325})

V_667 = Vertex(name = 'V_667',
               particles = [ P.ve__tilde__, P.e__minus__, P.Z, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_516})

V_668 = Vertex(name = 'V_668',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV12, L.FFV2, L.FFV6, L.FFV9 ],
               couplings = {(0,2):C.GC_243,(0,1):C.GC_317,(0,0):C.GC_617,(0,4):C.GC_536,(0,3):C.GC_651})

V_669 = Vertex(name = 'V_669',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV12, L.FFV2 ],
               couplings = {(0,1):C.GC_647,(0,0):C.GC_637})

V_670 = Vertex(name = 'V_670',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_370,(0,1):C.GC_376})

V_671 = Vertex(name = 'V_671',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_367,(0,1):C.GC_375})

V_672 = Vertex(name = 'V_672',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_370,(0,1):C.GC_376})

V_673 = Vertex(name = 'V_673',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_575,(0,0):C.GC_578})

V_674 = Vertex(name = 'V_674',
               particles = [ P.vm__tilde__, P.mu__minus__, P.Z, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_329})

V_675 = Vertex(name = 'V_675',
               particles = [ P.vm__tilde__, P.mu__minus__, P.Z, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_328})

V_676 = Vertex(name = 'V_676',
               particles = [ P.vm__tilde__, P.mu__minus__, P.Z, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_518})

V_677 = Vertex(name = 'V_677',
               particles = [ P.c__tilde__, P.s, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_338})

V_678 = Vertex(name = 'V_678',
               particles = [ P.c__tilde__, P.s, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_337})

V_679 = Vertex(name = 'V_679',
               particles = [ P.c__tilde__, P.s, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_524})

V_680 = Vertex(name = 'V_680',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV2, L.FFV6, L.FFV7 ],
               couplings = {(0,1):C.GC_243,(0,3):C.GC_316,(0,0):C.GC_613,(0,2):C.GC_641})

V_681 = Vertex(name = 'V_681',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV7 ],
               couplings = {(0,0):C.GC_655,(0,1):C.GC_636})

V_682 = Vertex(name = 'V_682',
               particles = [ P.s__tilde__, P.s, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_384,(0,1):C.GC_355})

V_683 = Vertex(name = 'V_683',
               particles = [ P.s__tilde__, P.s, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_381,(0,1):C.GC_354})

V_684 = Vertex(name = 'V_684',
               particles = [ P.s__tilde__, P.s, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_384,(0,1):C.GC_355})

V_685 = Vertex(name = 'V_685',
               particles = [ P.s__tilde__, P.s, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_582,(0,0):C.GC_566})

V_686 = Vertex(name = 'V_686',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV10, L.FFV12, L.FFV2, L.FFV6 ],
               couplings = {(0,2):C.GC_243,(0,1):C.GC_317,(0,0):C.GC_617,(0,3):C.GC_657})

V_687 = Vertex(name = 'V_687',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV12, L.FFV2 ],
               couplings = {(0,1):C.GC_649,(0,0):C.GC_639})

V_688 = Vertex(name = 'V_688',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_374,(0,1):C.GC_388})

V_689 = Vertex(name = 'V_689',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_371,(0,1):C.GC_387})

V_690 = Vertex(name = 'V_690',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_374,(0,1):C.GC_388})

V_691 = Vertex(name = 'V_691',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_577,(0,0):C.GC_584})

V_692 = Vertex(name = 'V_692',
               particles = [ P.vt__tilde__, P.ta__minus__, P.Z, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_332})

V_693 = Vertex(name = 'V_693',
               particles = [ P.vt__tilde__, P.ta__minus__, P.Z, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_331})

V_694 = Vertex(name = 'V_694',
               particles = [ P.vt__tilde__, P.ta__minus__, P.Z, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_520})

V_695 = Vertex(name = 'V_695',
               particles = [ P.b__tilde__, P.t, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_335})

V_696 = Vertex(name = 'V_696',
               particles = [ P.b__tilde__, P.t, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_336})

V_697 = Vertex(name = 'V_697',
               particles = [ P.b__tilde__, P.t, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS10, L.FFVS3 ],
               couplings = {(0,0):C.GC_346,(0,1):C.GC_523})

V_698 = Vertex(name = 'V_698',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV13, L.FFV2, L.FFV26, L.FFV6 ],
               couplings = {(0,2):C.GC_242,(0,1):C.GC_316,(0,0):C.GC_614,(0,4):C.GC_656,(0,3):C.GC_564})

V_699 = Vertex(name = 'V_699',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV13, L.FFV2 ],
               couplings = {(0,1):C.GC_652,(0,0):C.GC_636})

V_700 = Vertex(name = 'V_700',
               particles = [ P.t__tilde__, P.t, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_378,(0,1):C.GC_386})

V_701 = Vertex(name = 'V_701',
               particles = [ P.t__tilde__, P.t, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_379,(0,1):C.GC_385})

V_702 = Vertex(name = 'V_702',
               particles = [ P.t__tilde__, P.t, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_378,(0,1):C.GC_386})

V_703 = Vertex(name = 'V_703',
               particles = [ P.t__tilde__, P.t, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3, L.FFVS31 ],
               couplings = {(0,1):C.GC_579,(0,0):C.GC_583,(0,2):C.GC_348})

V_704 = Vertex(name = 'V_704',
               particles = [ P.d__tilde__, P.u, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_338})

V_705 = Vertex(name = 'V_705',
               particles = [ P.d__tilde__, P.u, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_339})

V_706 = Vertex(name = 'V_706',
               particles = [ P.d__tilde__, P.u, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_525})

V_707 = Vertex(name = 'V_707',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV10, L.FFV13, L.FFV2, L.FFV6 ],
               couplings = {(0,2):C.GC_242,(0,1):C.GC_316,(0,0):C.GC_614,(0,3):C.GC_658})

V_708 = Vertex(name = 'V_708',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV13, L.FFV2 ],
               couplings = {(0,1):C.GC_654,(0,0):C.GC_636})

V_709 = Vertex(name = 'V_709',
               particles = [ P.u__tilde__, P.u, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_382,(0,1):C.GC_390})

V_710 = Vertex(name = 'V_710',
               particles = [ P.u__tilde__, P.u, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_383,(0,1):C.GC_389})

V_711 = Vertex(name = 'V_711',
               particles = [ P.u__tilde__, P.u, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_382,(0,1):C.GC_390})

V_712 = Vertex(name = 'V_712',
               particles = [ P.u__tilde__, P.u, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS15, L.FFVS3 ],
               couplings = {(0,1):C.GC_581,(0,0):C.GC_585})

V_713 = Vertex(name = 'V_713',
               particles = [ P.e__plus__, P.ve, P.Z, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_326})

V_714 = Vertex(name = 'V_714',
               particles = [ P.e__plus__, P.ve, P.Z, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_327})

V_715 = Vertex(name = 'V_715',
               particles = [ P.e__plus__, P.ve, P.Z, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_517})

V_716 = Vertex(name = 'V_716',
               particles = [ P.ve__tilde__, P.ve, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_342})

V_717 = Vertex(name = 'V_717',
               particles = [ P.ve__tilde__, P.ve, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_646})

V_718 = Vertex(name = 'V_718',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_364})

V_719 = Vertex(name = 'V_719',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_365})

V_720 = Vertex(name = 'V_720',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_364})

V_721 = Vertex(name = 'V_721',
               particles = [ P.ve__tilde__, P.ve, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_572})

V_722 = Vertex(name = 'V_722',
               particles = [ P.mu__plus__, P.vm, P.Z, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_329})

V_723 = Vertex(name = 'V_723',
               particles = [ P.mu__plus__, P.vm, P.Z, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_330})

V_724 = Vertex(name = 'V_724',
               particles = [ P.mu__plus__, P.vm, P.Z, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_519})

V_725 = Vertex(name = 'V_725',
               particles = [ P.vm__tilde__, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_342})

V_726 = Vertex(name = 'V_726',
               particles = [ P.vm__tilde__, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_648})

V_727 = Vertex(name = 'V_727',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_368})

V_728 = Vertex(name = 'V_728',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_369})

V_729 = Vertex(name = 'V_729',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_368})

V_730 = Vertex(name = 'V_730',
               particles = [ P.vm__tilde__, P.vm, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_574})

V_731 = Vertex(name = 'V_731',
               particles = [ P.ta__plus__, P.vt, P.Z, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_332})

V_732 = Vertex(name = 'V_732',
               particles = [ P.ta__plus__, P.vt, P.Z, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_333})

V_733 = Vertex(name = 'V_733',
               particles = [ P.ta__plus__, P.vt, P.Z, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_521})

V_734 = Vertex(name = 'V_734',
               particles = [ P.vt__tilde__, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_342})

V_735 = Vertex(name = 'V_735',
               particles = [ P.vt__tilde__, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_650})

V_736 = Vertex(name = 'V_736',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_372})

V_737 = Vertex(name = 'V_737',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_373})

V_738 = Vertex(name = 'V_738',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_372})

V_739 = Vertex(name = 'V_739',
               particles = [ P.vt__tilde__, P.vt, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.FFVS3 ],
               couplings = {(0,0):C.GC_576})

V_740 = Vertex(name = 'V_740',
               particles = [ P.b__tilde__, P.b, P.b__tilde__, P.b ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF10, L.FFFF11, L.FFFF13, L.FFFF14, L.FFFF2, L.FFFF3 ],
               couplings = {(1,4):C.GC_41,(0,5):C.GC_41,(0,0):C.GC_102,(2,0):C.GC_103,(1,2):C.GC_102,(3,2):C.GC_103,(1,1):C.GC_102,(3,1):C.GC_103,(0,3):C.GC_102,(2,3):C.GC_103})

V_741 = Vertex(name = 'V_741',
               particles = [ P.b__tilde__, P.b, P.c__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF14, L.FFFF3 ],
               couplings = {(0,1):C.GC_38,(1,1):C.GC_42,(0,0):C.GC_118,(1,0):C.GC_119})

V_742 = Vertex(name = 'V_742',
               particles = [ P.t__tilde__, P.b, P.s__tilde__, P.c ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_112,(1,0):C.GC_114})

V_743 = Vertex(name = 'V_743',
               particles = [ P.b__tilde__, P.b, P.d__tilde__, P.d ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF14, L.FFFF3 ],
               couplings = {(0,1):C.GC_39,(1,1):C.GC_43,(0,0):C.GC_102,(1,0):C.GC_103})

V_744 = Vertex(name = 'V_744',
               particles = [ P.b__tilde__, P.b, P.e__plus__, P.e__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF14, L.FFFF3 ],
               couplings = {(0,1):C.GC_35,(0,0):C.GC_104})

V_745 = Vertex(name = 'V_745',
               particles = [ P.e__plus__, P.e__minus__, P.e__plus__, P.e__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_72,(0,1):C.GC_72})

V_746 = Vertex(name = 'V_746',
               particles = [ P.mu__plus__, P.e__minus__, P.e__plus__, P.mu__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_73,(0,1):C.GC_75})

V_747 = Vertex(name = 'V_747',
               particles = [ P.b__tilde__, P.b, P.mu__plus__, P.mu__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF14, L.FFFF3 ],
               couplings = {(0,1):C.GC_36,(0,0):C.GC_111})

V_748 = Vertex(name = 'V_748',
               particles = [ P.mu__plus__, P.mu__minus__, P.mu__plus__, P.mu__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_77,(0,1):C.GC_77})

V_749 = Vertex(name = 'V_749',
               particles = [ P.b__tilde__, P.b, P.s__tilde__, P.s ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF14, L.FFFF3 ],
               couplings = {(0,1):C.GC_39,(1,1):C.GC_43,(0,0):C.GC_102,(1,0):C.GC_103})

V_750 = Vertex(name = 'V_750',
               particles = [ P.ta__plus__, P.e__minus__, P.e__plus__, P.ta__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_74,(0,1):C.GC_76})

V_751 = Vertex(name = 'V_751',
               particles = [ P.ta__plus__, P.mu__minus__, P.mu__plus__, P.ta__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_78,(0,1):C.GC_79})

V_752 = Vertex(name = 'V_752',
               particles = [ P.ta__plus__, P.ta__minus__, P.ta__plus__, P.ta__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_80,(0,1):C.GC_80})

V_753 = Vertex(name = 'V_753',
               particles = [ P.ve__tilde__, P.e__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_105})

V_754 = Vertex(name = 'V_754',
               particles = [ P.vm__tilde__, P.mu__minus__, P.b__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_106})

V_755 = Vertex(name = 'V_755',
               particles = [ P.e__plus__, P.e__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF10, L.FFFF14, L.FFFF15, L.FFFF3 ],
               couplings = {(0,3):C.GC_108,(0,0):C.GC_104,(0,1):C.GC_123,(0,2):C.GC_122})

V_756 = Vertex(name = 'V_756',
               particles = [ P.mu__plus__, P.mu__minus__, P.t__tilde__, P.t ],
               color = [ 'Identity(3,4)' ],
               lorentz = [ L.FFFF10, L.FFFF14, L.FFFF15, L.FFFF3 ],
               couplings = {(0,3):C.GC_109,(0,0):C.GC_111,(0,1):C.GC_124,(0,2):C.GC_132})

V_757 = Vertex(name = 'V_757',
               particles = [ P.t__tilde__, P.b, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF11, L.FFFF12, L.FFFF13, L.FFFF2, L.FFFF3 ],
               couplings = {(1,3):C.GC_40,(0,4):C.GC_113,(1,2):C.GC_102,(2,2):C.GC_103,(1,0):C.GC_115,(2,0):C.GC_116,(1,1):C.GC_120,(2,1):C.GC_121})

V_758 = Vertex(name = 'V_758',
               particles = [ P.u__tilde__, P.d, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_112,(1,0):C.GC_114})

V_759 = Vertex(name = 'V_759',
               particles = [ P.c__tilde__, P.s, P.b__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_112,(1,0):C.GC_114})

V_760 = Vertex(name = 'V_760',
               particles = [ P.c__tilde__, P.c, P.t__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF10, L.FFFF14, L.FFFF15, L.FFFF3 ],
               couplings = {(0,3):C.GC_39,(1,3):C.GC_43,(0,0):C.GC_118,(1,0):C.GC_119,(0,1):C.GC_141,(1,1):C.GC_142,(0,2):C.GC_145,(1,2):C.GC_146})

V_761 = Vertex(name = 'V_761',
               particles = [ P.d__tilde__, P.d, P.t__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF10, L.FFFF14, L.FFFF15, L.FFFF3 ],
               couplings = {(0,3):C.GC_38,(1,3):C.GC_42,(0,0):C.GC_102,(1,0):C.GC_103,(0,1):C.GC_141,(1,1):C.GC_142,(0,2):C.GC_120,(1,2):C.GC_121})

V_762 = Vertex(name = 'V_762',
               particles = [ P.s__tilde__, P.s, P.t__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF10, L.FFFF14, L.FFFF15, L.FFFF3 ],
               couplings = {(0,3):C.GC_38,(1,3):C.GC_42,(0,0):C.GC_102,(1,0):C.GC_103,(0,1):C.GC_141,(1,1):C.GC_142,(0,2):C.GC_120,(1,2):C.GC_121})

V_763 = Vertex(name = 'V_763',
               particles = [ P.t__tilde__, P.t, P.t__tilde__, P.t ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'Identity(1,4)*Identity(2,3)', 'T(-1,2,1)*T(-1,4,3)', 'T(-1,2,3)*T(-1,4,1)' ],
               lorentz = [ L.FFFF10, L.FFFF11, L.FFFF12, L.FFFF13, L.FFFF14, L.FFFF15, L.FFFF2, L.FFFF3 ],
               couplings = {(1,6):C.GC_41,(0,7):C.GC_41,(0,0):C.GC_115,(2,0):C.GC_116,(1,3):C.GC_115,(3,3):C.GC_116,(1,1):C.GC_115,(3,1):C.GC_116,(1,2):C.GC_143,(0,4):C.GC_115,(2,4):C.GC_116,(0,5):C.GC_143})

V_764 = Vertex(name = 'V_764',
               particles = [ P.t__tilde__, P.b, P.d__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_112,(1,0):C.GC_114})

V_765 = Vertex(name = 'V_765',
               particles = [ P.b__tilde__, P.b, P.u__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF14, L.FFFF3 ],
               couplings = {(0,1):C.GC_38,(1,1):C.GC_42,(0,0):C.GC_118,(1,0):C.GC_119})

V_766 = Vertex(name = 'V_766',
               particles = [ P.t__tilde__, P.t, P.u__tilde__, P.u ],
               color = [ 'Identity(1,2)*Identity(3,4)', 'T(-1,2,1)*T(-1,4,3)' ],
               lorentz = [ L.FFFF10, L.FFFF14, L.FFFF15, L.FFFF3 ],
               couplings = {(0,3):C.GC_39,(1,3):C.GC_43,(0,0):C.GC_141,(1,0):C.GC_142,(0,1):C.GC_118,(1,1):C.GC_119,(0,2):C.GC_145,(1,2):C.GC_146})

V_767 = Vertex(name = 'V_767',
               particles = [ P.t__tilde__, P.b, P.e__plus__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_105})

V_768 = Vertex(name = 'V_768',
               particles = [ P.b__tilde__, P.b, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_108})

V_769 = Vertex(name = 'V_769',
               particles = [ P.e__plus__, P.e__minus__, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_72})

V_770 = Vertex(name = 'V_770',
               particles = [ P.mu__plus__, P.mu__minus__, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_73})

V_771 = Vertex(name = 'V_771',
               particles = [ P.ta__plus__, P.ta__minus__, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_74})

V_772 = Vertex(name = 'V_772',
               particles = [ P.t__tilde__, P.t, P.ve__tilde__, P.ve ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF10, L.FFFF3 ],
               couplings = {(0,1):C.GC_35,(0,0):C.GC_123})

V_773 = Vertex(name = 'V_773',
               particles = [ P.ve__tilde__, P.ve, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_72,(0,1):C.GC_72})

V_774 = Vertex(name = 'V_774',
               particles = [ P.e__plus__, P.mu__minus__, P.vm__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_75})

V_775 = Vertex(name = 'V_775',
               particles = [ P.e__plus__, P.ta__minus__, P.vt__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_76})

V_776 = Vertex(name = 'V_776',
               particles = [ P.t__tilde__, P.b, P.mu__plus__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_106})

V_777 = Vertex(name = 'V_777',
               particles = [ P.mu__plus__, P.e__minus__, P.ve__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_75})

V_778 = Vertex(name = 'V_778',
               particles = [ P.vm__tilde__, P.ve, P.ve__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_73,(0,1):C.GC_75})

V_779 = Vertex(name = 'V_779',
               particles = [ P.b__tilde__, P.b, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_109})

V_780 = Vertex(name = 'V_780',
               particles = [ P.e__plus__, P.e__minus__, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_73})

V_781 = Vertex(name = 'V_781',
               particles = [ P.mu__plus__, P.mu__minus__, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_77})

V_782 = Vertex(name = 'V_782',
               particles = [ P.ta__plus__, P.ta__minus__, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_78})

V_783 = Vertex(name = 'V_783',
               particles = [ P.t__tilde__, P.t, P.vm__tilde__, P.vm ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF10, L.FFFF3 ],
               couplings = {(0,1):C.GC_36,(0,0):C.GC_124})

V_784 = Vertex(name = 'V_784',
               particles = [ P.vm__tilde__, P.vm, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_77,(0,1):C.GC_77})

V_785 = Vertex(name = 'V_785',
               particles = [ P.mu__plus__, P.ta__minus__, P.vt__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_79})

V_786 = Vertex(name = 'V_786',
               particles = [ P.ta__plus__, P.e__minus__, P.ve__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_76})

V_787 = Vertex(name = 'V_787',
               particles = [ P.vt__tilde__, P.ve, P.ve__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_74,(0,1):C.GC_76})

V_788 = Vertex(name = 'V_788',
               particles = [ P.ta__plus__, P.mu__minus__, P.vm__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_79})

V_789 = Vertex(name = 'V_789',
               particles = [ P.vt__tilde__, P.vm, P.vm__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_78,(0,1):C.GC_79})

V_790 = Vertex(name = 'V_790',
               particles = [ P.b__tilde__, P.b, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_110})

V_791 = Vertex(name = 'V_791',
               particles = [ P.e__plus__, P.e__minus__, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_74})

V_792 = Vertex(name = 'V_792',
               particles = [ P.mu__plus__, P.mu__minus__, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_78})

V_793 = Vertex(name = 'V_793',
               particles = [ P.ta__plus__, P.ta__minus__, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF3 ],
               couplings = {(0,0):C.GC_80})

V_794 = Vertex(name = 'V_794',
               particles = [ P.t__tilde__, P.t, P.vt__tilde__, P.vt ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFFF10, L.FFFF3 ],
               couplings = {(0,1):C.GC_37,(0,0):C.GC_125})

V_795 = Vertex(name = 'V_795',
               particles = [ P.vt__tilde__, P.vt, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.FFFF2, L.FFFF3 ],
               couplings = {(0,0):C.GC_80,(0,1):C.GC_80})

V_796 = Vertex(name = 'V_796',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS10 ],
               couplings = {(0,0):C.GC_147})

V_797 = Vertex(name = 'V_797',
               particles = [ P.t__tilde__, P.t, P.a, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_353})

V_798 = Vertex(name = 'V_798',
               particles = [ P.t__tilde__, P.t, P.a, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS19 ],
               couplings = {(0,0):C.GC_352})

V_799 = Vertex(name = 'V_799',
               particles = [ P.t__tilde__, P.t, P.Z, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_349})

V_800 = Vertex(name = 'V_800',
               particles = [ P.b__tilde__, P.t, P.g, P.G__minus__ ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS10 ],
               couplings = {(0,0):C.GC_195})

V_801 = Vertex(name = 'V_801',
               particles = [ P.t__tilde__, P.t, P.g, P.G0 ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_198})

V_802 = Vertex(name = 'V_802',
               particles = [ P.t__tilde__, P.t, P.g, P.H ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS19 ],
               couplings = {(0,0):C.GC_197})

V_803 = Vertex(name = 'V_803',
               particles = [ P.b__tilde__, P.t, P.g, P.g, P.G__minus__ ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_201})

V_804 = Vertex(name = 'V_804',
               particles = [ P.t__tilde__, P.t, P.g, P.g, P.G0 ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS57 ],
               couplings = {(0,0):C.GC_204})

V_805 = Vertex(name = 'V_805',
               particles = [ P.t__tilde__, P.t, P.g, P.g, P.H ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS63 ],
               couplings = {(0,0):C.GC_203})

V_806 = Vertex(name = 'V_806',
               particles = [ P.t__tilde__, P.t, P.g, P.g ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVV32 ],
               couplings = {(0,0):C.GC_460})

V_807 = Vertex(name = 'V_807',
               particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_182})

V_808 = Vertex(name = 'V_808',
               particles = [ P.b__tilde__, P.t, P.a, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_184})

V_809 = Vertex(name = 'V_809',
               particles = [ P.b__tilde__, P.t, P.a, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV12 ],
               couplings = {(0,0):C.GC_454})

V_810 = Vertex(name = 'V_810',
               particles = [ P.t__tilde__, P.t, P.a, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_185})

V_811 = Vertex(name = 'V_811',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_294})

V_812 = Vertex(name = 'V_812',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS6 ],
               couplings = {(0,0):C.GC_296})

V_813 = Vertex(name = 'V_813',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS43 ],
               couplings = {(0,0):C.GC_295})

V_814 = Vertex(name = 'V_814',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV23 ],
               couplings = {(0,0):C.GC_504})

V_815 = Vertex(name = 'V_815',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_300})

V_816 = Vertex(name = 'V_816',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_302})

V_817 = Vertex(name = 'V_817',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV12 ],
               couplings = {(0,0):C.GC_508})

V_818 = Vertex(name = 'V_818',
               particles = [ P.t__tilde__, P.t, P.W__plus__, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS22 ],
               couplings = {(0,0):C.GC_303})

V_819 = Vertex(name = 'V_819',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS30 ],
               couplings = {(0,0):C.GC_149})

V_820 = Vertex(name = 'V_820',
               particles = [ P.t__tilde__, P.b, P.g, P.G__plus__ ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFVS30 ],
               couplings = {(0,0):C.GC_196})

V_821 = Vertex(name = 'V_821',
               particles = [ P.t__tilde__, P.b, P.g, P.g, P.G__plus__ ],
               color = [ 'f(-1,3,4)*T(-1,2,1)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_202})

V_822 = Vertex(name = 'V_822',
               particles = [ P.t__tilde__, P.t, P.a, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_185})

V_823 = Vertex(name = 'V_823',
               particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_182})

V_824 = Vertex(name = 'V_824',
               particles = [ P.t__tilde__, P.b, P.a, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_183})

V_825 = Vertex(name = 'V_825',
               particles = [ P.t__tilde__, P.b, P.a, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV31 ],
               couplings = {(0,0):C.GC_453})

V_826 = Vertex(name = 'V_826',
               particles = [ P.t__tilde__, P.b, P.W__minus__, P.W__plus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_293})

V_827 = Vertex(name = 'V_827',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_303})

V_828 = Vertex(name = 'V_828',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_300})

V_829 = Vertex(name = 'V_829',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVVS62 ],
               couplings = {(0,0):C.GC_301})

V_830 = Vertex(name = 'V_830',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVV31 ],
               couplings = {(0,0):C.GC_507})

