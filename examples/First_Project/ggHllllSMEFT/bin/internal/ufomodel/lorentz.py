# This file was automatically created by FeynRules 2.4.72
# Mathematica version: 11.3.0 for Mac OS X x86 (64-bit) (March 7, 2018)
# Date: Thu 8 Aug 2019 15:25:02


from object_library import all_lorentz, Lorentz

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot
try:
   import form_factors as ForFac 
except ImportError:
   pass


FF1 = Lorentz(name = 'FF1',
              spins = [ 2, 2 ],
              structure = 'P(-1,1)*Gamma(-1,2,1)')

FF2 = Lorentz(name = 'FF2',
              spins = [ 2, 2 ],
              structure = 'ProjM(2,1) + ProjP(2,1)')

FF3 = Lorentz(name = 'FF3',
              spins = [ 2, 2 ],
              structure = 'P(-1,1)**2*ProjM(2,1) + P(-1,1)**2*ProjP(2,1)')

FF4 = Lorentz(name = 'FF4',
              spins = [ 2, 2 ],
              structure = 'P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,2,-2)*ProjP(-2,1)')

FF5 = Lorentz(name = 'FF5',
              spins = [ 2, 2 ],
              structure = '-(P(-1,1)*Gamma(-1,2,1)) + P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,2,-2)*ProjP(-2,1)')

FF6 = Lorentz(name = 'FF6',
              spins = [ 2, 2 ],
              structure = '(P(-1,1)*Gamma(-1,2,1))/2. + P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,2,-2)*ProjP(-2,1)')

VV1 = Lorentz(name = 'VV1',
              spins = [ 3, 3 ],
              structure = 'P(1,2)*P(2,2)')

VV2 = Lorentz(name = 'VV2',
              spins = [ 3, 3 ],
              structure = 'Metric(1,2)')

VV3 = Lorentz(name = 'VV3',
              spins = [ 3, 3 ],
              structure = 'P(-1,2)**2*Metric(1,2)')

VV4 = Lorentz(name = 'VV4',
              spins = [ 3, 3 ],
              structure = 'P(1,2)*P(2,2) - (3*P(-1,2)**2*Metric(1,2))/2.')

VV5 = Lorentz(name = 'VV5',
              spins = [ 3, 3 ],
              structure = 'P(1,2)*P(2,2) - P(-1,2)**2*Metric(1,2)')

UUS1 = Lorentz(name = 'UUS1',
               spins = [ -1, -1, 1 ],
               structure = '1')

UUV1 = Lorentz(name = 'UUV1',
               spins = [ -1, -1, 3 ],
               structure = 'P(3,2) + P(3,3)')

SSS1 = Lorentz(name = 'SSS1',
               spins = [ 1, 1, 1 ],
               structure = '1')

SSS2 = Lorentz(name = 'SSS2',
               spins = [ 1, 1, 1 ],
               structure = 'P(-1,1)*P(-1,2)')

SSS3 = Lorentz(name = 'SSS3',
               spins = [ 1, 1, 1 ],
               structure = 'P(-1,1)*P(-1,2) - P(-1,1)*P(-1,3)')

SSS4 = Lorentz(name = 'SSS4',
               spins = [ 1, 1, 1 ],
               structure = 'P(-1,1)*P(-1,3) + P(-1,2)*P(-1,3)')

SSS5 = Lorentz(name = 'SSS5',
               spins = [ 1, 1, 1 ],
               structure = 'P(-1,1)*P(-1,2) + P(-1,1)*P(-1,3) + P(-1,2)*P(-1,3)')

FFS1 = Lorentz(name = 'FFS1',
               spins = [ 2, 2, 1 ],
               structure = 'Gamma5(2,1)')

FFS2 = Lorentz(name = 'FFS2',
               spins = [ 2, 2, 1 ],
               structure = 'P(-1,1)*Gamma(-1,2,1) + (P(-1,3)*Gamma(-1,2,1))/2.')

FFS3 = Lorentz(name = 'FFS3',
               spins = [ 2, 2, 1 ],
               structure = 'Identity(2,1)')

FFS4 = Lorentz(name = 'FFS4',
               spins = [ 2, 2, 1 ],
               structure = 'ProjM(2,1)')

FFS5 = Lorentz(name = 'FFS5',
               spins = [ 2, 2, 1 ],
               structure = 'P(-1,1)**2*ProjM(2,1)')

FFS6 = Lorentz(name = 'FFS6',
               spins = [ 2, 2, 1 ],
               structure = 'P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1)')

FFS7 = Lorentz(name = 'FFS7',
               spins = [ 2, 2, 1 ],
               structure = 'P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1)')

FFS8 = Lorentz(name = 'FFS8',
               spins = [ 2, 2, 1 ],
               structure = 'P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1) - (P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1))/2.')

FFS9 = Lorentz(name = 'FFS9',
               spins = [ 2, 2, 1 ],
               structure = 'P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1)')

FFS10 = Lorentz(name = 'FFS10',
                spins = [ 2, 2, 1 ],
                structure = 'P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1) + (3*P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1))/2.')

FFS11 = Lorentz(name = 'FFS11',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjM(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjM(-3,1) + 3*P(-1,1)**2*ProjM(2,1) + 4*P(-1,1)*P(-1,3)*ProjM(2,1) + 3*P(-1,3)**2*ProjM(2,1)')

FFS12 = Lorentz(name = 'FFS12',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjM(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjM(-3,1) + (9*P(-1,1)**2*ProjM(2,1))/2. + 7*P(-1,1)*P(-1,3)*ProjM(2,1) + (9*P(-1,3)**2*ProjM(2,1))/2.')

FFS13 = Lorentz(name = 'FFS13',
                spins = [ 2, 2, 1 ],
                structure = 'ProjM(2,1) - ProjP(2,1)')

FFS14 = Lorentz(name = 'FFS14',
                spins = [ 2, 2, 1 ],
                structure = 'ProjP(2,1)')

FFS15 = Lorentz(name = 'FFS15',
                spins = [ 2, 2, 1 ],
                structure = 'P(-1,1)**2*ProjP(2,1)')

FFS16 = Lorentz(name = 'FFS16',
                spins = [ 2, 2, 1 ],
                structure = 'ProjM(2,1) + ProjP(2,1)')

FFS17 = Lorentz(name = 'FFS17',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma5(-3,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3) + P(-2,3)*P(-1,1)*Gamma5(-3,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4) - 6*P(-1,1)**2*ProjM(2,1) - 4*P(-1,1)*P(-1,3)*ProjM(2,1) - 3*P(-1,3)**2*ProjM(2,1) + 6*P(-1,1)**2*ProjP(2,1) + 4*P(-1,1)*P(-1,3)*ProjP(2,1) + 3*P(-1,3)**2*ProjP(2,1)')

FFS18 = Lorentz(name = 'FFS18',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma5(-3,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3) + P(-2,3)*P(-1,1)*Gamma5(-3,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4) - 9*P(-1,1)**2*ProjM(2,1) - 7*P(-1,1)*P(-1,3)*ProjM(2,1) - (9*P(-1,3)**2*ProjM(2,1))/2. + 9*P(-1,1)**2*ProjP(2,1) + 7*P(-1,1)*P(-1,3)*ProjP(2,1) + (9*P(-1,3)**2*ProjP(2,1))/2.')

FFS19 = Lorentz(name = 'FFS19',
                spins = [ 2, 2, 1 ],
                structure = 'P(-1,3)*Gamma(-1,2,-2)*ProjP(-2,1)')

FFS20 = Lorentz(name = 'FFS20',
                spins = [ 2, 2, 1 ],
                structure = 'P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,3)*Gamma(-1,2,-2)*ProjP(-2,1)')

FFS21 = Lorentz(name = 'FFS21',
                spins = [ 2, 2, 1 ],
                structure = 'P(-1,1)*Gamma(-1,2,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1))/2. + P(-1,1)*Gamma(-1,2,-2)*ProjP(-2,1) + (P(-1,3)*Gamma(-1,2,-2)*ProjP(-2,1))/2.')

FFS22 = Lorentz(name = 'FFS22',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjP(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjP(-3,1) + 3*P(-1,1)**2*ProjP(2,1) + 4*P(-1,1)*P(-1,3)*ProjP(2,1) + 3*P(-1,3)**2*ProjP(2,1)')

FFS23 = Lorentz(name = 'FFS23',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjM(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjM(-3,1) + 6*P(-1,1)**2*ProjM(2,1) + 4*P(-1,1)*P(-1,3)*ProjM(2,1) + 3*P(-1,3)**2*ProjM(2,1) + P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjP(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjP(-3,1) + 6*P(-1,1)**2*ProjP(2,1) + 4*P(-1,1)*P(-1,3)*ProjP(2,1) + 3*P(-1,3)**2*ProjP(2,1)')

FFS24 = Lorentz(name = 'FFS24',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjP(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjP(-3,1) + (9*P(-1,1)**2*ProjP(2,1))/2. + 7*P(-1,1)*P(-1,3)*ProjP(2,1) + (9*P(-1,3)**2*ProjP(2,1))/2.')

FFS25 = Lorentz(name = 'FFS25',
                spins = [ 2, 2, 1 ],
                structure = 'P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjM(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjM(-3,1) + 9*P(-1,1)**2*ProjM(2,1) + 7*P(-1,1)*P(-1,3)*ProjM(2,1) + (9*P(-1,3)**2*ProjM(2,1))/2. + P(-2,3)*P(-1,1)*Gamma(-2,2,-4)*Gamma(-1,-4,-3)*ProjP(-3,1) + P(-2,3)*P(-1,1)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*ProjP(-3,1) + 9*P(-1,1)**2*ProjP(2,1) + 7*P(-1,1)*P(-1,3)*ProjP(2,1) + (9*P(-1,3)**2*ProjP(2,1))/2.')

FFV1 = Lorentz(name = 'FFV1',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,1)')

FFV2 = Lorentz(name = 'FFV2',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1)')

FFV3 = Lorentz(name = 'FFV3',
               spins = [ 2, 2, 3 ],
               structure = '(13*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + (8*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (4*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/11.')

FFV4 = Lorentz(name = 'FFV4',
               spins = [ 2, 2, 3 ],
               structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3.')

FFV5 = Lorentz(name = 'FFV5',
               spins = [ 2, 2, 3 ],
               structure = '-(P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1)) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1)')

FFV6 = Lorentz(name = 'FFV6',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjP(-1,1)')

FFV7 = Lorentz(name = 'FFV7',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) - 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV8 = Lorentz(name = 'FFV8',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) - (5*Gamma(3,2,-1)*ProjP(-1,1))/3.')

FFV9 = Lorentz(name = 'FFV9',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + (2*Gamma(3,2,-1)*ProjP(-1,1))/3.')

FFV10 = Lorentz(name = 'FFV10',
                spins = [ 2, 2, 3 ],
                structure = 'Gamma(3,2,-1)*ProjM(-1,1) + Gamma(3,2,-1)*ProjP(-1,1)')

FFV11 = Lorentz(name = 'FFV11',
                spins = [ 2, 2, 3 ],
                structure = 'Gamma(3,2,-1)*ProjM(-1,1) + (8*Gamma(3,2,-1)*ProjP(-1,1))/5.')

FFV12 = Lorentz(name = 'FFV12',
                spins = [ 2, 2, 3 ],
                structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV13 = Lorentz(name = 'FFV13',
                spins = [ 2, 2, 3 ],
                structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 4*Gamma(3,2,-1)*ProjP(-1,1)')

FFV14 = Lorentz(name = 'FFV14',
                spins = [ 2, 2, 3 ],
                structure = 'P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2) - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1)')

FFV15 = Lorentz(name = 'FFV15',
                spins = [ 2, 2, 3 ],
                structure = 'P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2) - P(-1,3)*Gamma(-1,2,-2)*Gamma(3,-2,-3)*ProjM(-3,1) - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1)')

FFV16 = Lorentz(name = 'FFV16',
                spins = [ 2, 2, 3 ],
                structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFV17 = Lorentz(name = 'FFV17',
                spins = [ 2, 2, 3 ],
                structure = '(P(-1,3)*Gamma(-1,2,-2)*Gamma(3,-2,1))/3. + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFV18 = Lorentz(name = 'FFV18',
                spins = [ 2, 2, 3 ],
                structure = 'P(3,1)*Identity(2,1) - P(3,2)*Identity(2,1) - (P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/24. + (163*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/72. - (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/24. - (83*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/36. - (P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/24. + (163*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/72. - (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/24. - (83*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/36.')

FFV19 = Lorentz(name = 'FFV19',
                spins = [ 2, 2, 3 ],
                structure = 'P(3,1)*Identity(2,1) - P(3,2)*Identity(2,1) - (5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/42. + (101*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/84. - (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/42. - (37*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/28. - (5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/42. + (101*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/84. - (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/42. - (37*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/28.')

FFV20 = Lorentz(name = 'FFV20',
                spins = [ 2, 2, 3 ],
                structure = '(11*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/13. + (7*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/13. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (5*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/13.')

FFV21 = Lorentz(name = 'FFV21',
                spins = [ 2, 2, 3 ],
                structure = '(13*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + (8*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (4*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (7*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/11. + (13*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11. + (5*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11.')

FFV22 = Lorentz(name = 'FFV22',
                spins = [ 2, 2, 3 ],
                structure = '(P(-1,3)*Gamma(-1,2,-2)*Gamma(3,-2,1))/3. + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3. + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3.')

FFV23 = Lorentz(name = 'FFV23',
                spins = [ 2, 2, 3 ],
                structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3.')

FFV24 = Lorentz(name = 'FFV24',
                spins = [ 2, 2, 3 ],
                structure = 'P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFV25 = Lorentz(name = 'FFV25',
                spins = [ 2, 2, 3 ],
                structure = '-(P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1)) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFV26 = Lorentz(name = 'FFV26',
                spins = [ 2, 2, 3 ],
                structure = '-(P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1)) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFV27 = Lorentz(name = 'FFV27',
                spins = [ 2, 2, 3 ],
                structure = '(5*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (5*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

VSS1 = Lorentz(name = 'VSS1',
               spins = [ 3, 1, 1 ],
               structure = 'P(1,2)')

VSS2 = Lorentz(name = 'VSS2',
               spins = [ 3, 1, 1 ],
               structure = 'P(1,2) - P(1,3)')

VSS3 = Lorentz(name = 'VSS3',
               spins = [ 3, 1, 1 ],
               structure = 'P(1,2) + P(1,3)/3.')

VVS1 = Lorentz(name = 'VVS1',
               spins = [ 3, 3, 1 ],
               structure = 'Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,2)')

VVS2 = Lorentz(name = 'VVS2',
               spins = [ 3, 3, 1 ],
               structure = 'Epsilon(1,2,-1,-2)*P(-2,2)*P(-1,1) - Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,1) - Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,2)')

VVS3 = Lorentz(name = 'VVS3',
               spins = [ 3, 3, 1 ],
               structure = 'Metric(1,2)')

VVS4 = Lorentz(name = 'VVS4',
               spins = [ 3, 3, 1 ],
               structure = 'P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)')

VVS5 = Lorentz(name = 'VVS5',
               spins = [ 3, 3, 1 ],
               structure = 'P(1,2)*P(2,1) - (P(-1,1)*P(-1,2)*Metric(1,2))/2. + (P(-1,2)**2*Metric(1,2))/2. + (P(-1,2)*P(-1,3)*Metric(1,2))/2.')

VVS6 = Lorentz(name = 'VVS6',
               spins = [ 3, 3, 1 ],
               structure = 'P(1,2)*P(2,1) - (P(1,3)*P(2,1))/2. - (P(1,2)*P(2,3))/2. - (P(-1,1)*P(-1,2)*Metric(1,2))/2. + (P(-1,2)**2*Metric(1,2))/2. + (P(-1,1)*P(-1,3)*Metric(1,2))/2. + P(-1,2)*P(-1,3)*Metric(1,2)')

VVS7 = Lorentz(name = 'VVS7',
               spins = [ 3, 3, 1 ],
               structure = 'P(1,2)*P(2,1) - (9*P(1,3)*P(2,1))/25. - (9*P(1,2)*P(2,3))/25. - (16*P(-1,1)*P(-1,2)*Metric(1,2))/25. + (9*P(-1,2)**2*Metric(1,2))/25. + (3*P(-1,1)*P(-1,3)*Metric(1,2))/25. + (12*P(-1,2)*P(-1,3)*Metric(1,2))/25. - (6*P(-1,3)**2*Metric(1,2))/25.')

VVS8 = Lorentz(name = 'VVS8',
               spins = [ 3, 3, 1 ],
               structure = 'P(1,2)*P(2,1) - (13*P(1,3)*P(2,1))/62. - (13*P(1,2)*P(2,3))/62. + (P(1,3)*P(2,3))/62. - (20*P(-1,1)*P(-1,2)*Metric(1,2))/31. + (19*P(-1,2)**2*Metric(1,2))/62. + (2*P(-1,1)*P(-1,3)*Metric(1,2))/31. + (23*P(-1,2)*P(-1,3)*Metric(1,2))/62. - (13*P(-1,3)**2*Metric(1,2))/62.')

VVV1 = Lorentz(name = 'VVV1',
               spins = [ 3, 3, 3 ],
               structure = '-(Epsilon(1,2,3,-1)*P(-1,1)) + Epsilon(1,2,3,-1)*P(-1,2)')

VVV2 = Lorentz(name = 'VVV2',
               spins = [ 3, 3, 3 ],
               structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(1,2)*Metric(2,3)')

VVV3 = Lorentz(name = 'VVV3',
               spins = [ 3, 3, 3 ],
               structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

VVV4 = Lorentz(name = 'VVV4',
               spins = [ 3, 3, 3 ],
               structure = 'P(3,2)*Metric(1,2) - P(2,3)*Metric(1,3) - P(1,2)*Metric(2,3) + P(1,3)*Metric(2,3)')

VVV5 = Lorentz(name = 'VVV5',
               spins = [ 3, 3, 3 ],
               structure = '-(P(1,2)*P(2,3)*P(3,1)) + P(1,3)*P(2,1)*P(3,2) + P(-1,2)*P(-1,3)*P(3,1)*Metric(1,2) - P(-1,1)*P(-1,3)*P(3,2)*Metric(1,2) - P(-1,2)*P(-1,3)*P(2,1)*Metric(1,3) + P(-1,1)*P(-1,2)*P(2,3)*Metric(1,3) + P(-1,1)*P(-1,3)*P(1,2)*Metric(2,3) - P(-1,1)*P(-1,2)*P(1,3)*Metric(2,3)')

SSSS1 = Lorentz(name = 'SSSS1',
                spins = [ 1, 1, 1, 1 ],
                structure = '1')

SSSS2 = Lorentz(name = 'SSSS2',
                spins = [ 1, 1, 1, 1 ],
                structure = 'P(-1,1)*P(-1,3) + P(-1,2)*P(-1,3) + P(-1,1)*P(-1,4) + P(-1,2)*P(-1,4)')

SSSS3 = Lorentz(name = 'SSSS3',
                spins = [ 1, 1, 1, 1 ],
                structure = 'P(-1,1)*P(-1,2) + P(-1,3)*P(-1,4)')

SSSS4 = Lorentz(name = 'SSSS4',
                spins = [ 1, 1, 1, 1 ],
                structure = 'P(-1,1)*P(-1,2) - P(-1,1)*P(-1,3) - P(-1,2)*P(-1,4) + P(-1,3)*P(-1,4)')

SSSS5 = Lorentz(name = 'SSSS5',
                spins = [ 1, 1, 1, 1 ],
                structure = 'P(-1,1)*P(-1,2) + (P(-1,1)*P(-1,3))/2. + (P(-1,2)*P(-1,3))/2. + (P(-1,1)*P(-1,4))/2. + (P(-1,2)*P(-1,4))/2. + P(-1,3)*P(-1,4)')

SSSS6 = Lorentz(name = 'SSSS6',
                spins = [ 1, 1, 1, 1 ],
                structure = 'P(-1,1)*P(-1,2) + P(-1,1)*P(-1,3) + P(-1,2)*P(-1,3) + P(-1,1)*P(-1,4) + P(-1,2)*P(-1,4) + P(-1,3)*P(-1,4)')

FFSS1 = Lorentz(name = 'FFSS1',
                spins = [ 2, 2, 1, 1 ],
                structure = 'Identity(2,1)')

FFSS2 = Lorentz(name = 'FFSS2',
                spins = [ 2, 2, 1, 1 ],
                structure = 'ProjM(2,1)')

FFSS3 = Lorentz(name = 'FFSS3',
                spins = [ 2, 2, 1, 1 ],
                structure = 'P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,4)*Gamma(-1,2,-2)*ProjM(-2,1)')

FFSS4 = Lorentz(name = 'FFSS4',
                spins = [ 2, 2, 1, 1 ],
                structure = 'ProjM(2,1) - ProjP(2,1)')

FFSS5 = Lorentz(name = 'FFSS5',
                spins = [ 2, 2, 1, 1 ],
                structure = 'ProjP(2,1)')

FFSS6 = Lorentz(name = 'FFSS6',
                spins = [ 2, 2, 1, 1 ],
                structure = 'ProjM(2,1) + ProjP(2,1)')

FFSS7 = Lorentz(name = 'FFSS7',
                spins = [ 2, 2, 1, 1 ],
                structure = 'P(-1,3)*Gamma(-1,2,-2)*ProjP(-2,1) - P(-1,4)*Gamma(-1,2,-2)*ProjP(-2,1)')

FFFF1 = Lorentz(name = 'FFFF1',
                spins = [ 2, 2, 2, 2 ],
                structure = 'ProjM(2,1)*ProjM(4,3)')

FFFF2 = Lorentz(name = 'FFFF2',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-3,1)*ProjM(-2,3)')

FFFF3 = Lorentz(name = 'FFFF3',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-3,3)*ProjM(-2,1)')

FFFF4 = Lorentz(name = 'FFFF4',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-2,-4,-3)*Gamma(-2,2,-6)*Gamma(-1,-6,-5)*Gamma(-1,4,-4)*ProjM(-5,1)*ProjM(-3,3)')

FFFF5 = Lorentz(name = 'FFFF5',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-2,-4,-3)*Gamma(-2,4,-6)*Gamma(-1,-6,-5)*Gamma(-1,2,-4)*ProjM(-5,3)*ProjM(-3,1)')

FFFF6 = Lorentz(name = 'FFFF6',
                spins = [ 2, 2, 2, 2 ],
                structure = 'Gamma(-2,-6,-5)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*Gamma(-1,4,-6)*ProjM(-5,3)*ProjM(-3,1)')

FFFF7 = Lorentz(name = 'FFFF7',
                spins = [ 2, 2, 2, 2 ],
                structure = 'ProjM(4,3)*ProjP(2,1)')

FFFF8 = Lorentz(name = 'FFFF8',
                spins = [ 2, 2, 2, 2 ],
                structure = 'ProjM(2,1)*ProjP(4,3)')

FFFF9 = Lorentz(name = 'FFFF9',
                spins = [ 2, 2, 2, 2 ],
                structure = 'ProjP(2,1)*ProjP(4,3)')

FFFF10 = Lorentz(name = 'FFFF10',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-1,2,-3)*Gamma(-1,4,-2)*ProjM(-2,3)*ProjP(-3,1)')

FFFF11 = Lorentz(name = 'FFFF11',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-2,3)*ProjP(-3,1)')

FFFF12 = Lorentz(name = 'FFFF12',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjP(-3,1)*ProjP(-2,3)')

FFFF13 = Lorentz(name = 'FFFF13',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-1,2,-3)*Gamma(-1,4,-2)*ProjM(-2,1)*ProjP(-3,3)')

FFFF14 = Lorentz(name = 'FFFF14',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-2,1)*ProjP(-3,3)')

FFFF15 = Lorentz(name = 'FFFF15',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjP(-3,3)*ProjP(-2,1)')

FFFF16 = Lorentz(name = 'FFFF16',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-2,-4,-3)*Gamma(-2,2,-6)*Gamma(-1,-6,-5)*Gamma(-1,4,-4)*ProjP(-5,1)*ProjP(-3,3)')

FFFF17 = Lorentz(name = 'FFFF17',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-2,-4,-3)*Gamma(-2,4,-6)*Gamma(-1,-6,-5)*Gamma(-1,2,-4)*ProjP(-5,3)*ProjP(-3,1)')

FFFF18 = Lorentz(name = 'FFFF18',
                 spins = [ 2, 2, 2, 2 ],
                 structure = 'Gamma(-2,-6,-5)*Gamma(-2,-4,-3)*Gamma(-1,2,-4)*Gamma(-1,4,-6)*ProjP(-5,3)*ProjP(-3,1)')

FFVS1 = Lorentz(name = 'FFVS1',
                spins = [ 2, 2, 3, 1 ],
                structure = 'Gamma(3,2,1)')

FFVS2 = Lorentz(name = 'FFVS2',
                spins = [ 2, 2, 3, 1 ],
                structure = '-(P(-1,3)*Gamma5(-2,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)) + P(-1,3)*Gamma5(-2,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)')

FFVS3 = Lorentz(name = 'FFVS3',
                spins = [ 2, 2, 3, 1 ],
                structure = 'Gamma(3,2,-1)*ProjM(-1,1)')

FFVS4 = Lorentz(name = 'FFVS4',
                spins = [ 2, 2, 3, 1 ],
                structure = 'P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1)')

FFVS5 = Lorentz(name = 'FFVS5',
                spins = [ 2, 2, 3, 1 ],
                structure = '(P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + (71*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/36. - (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/12. - (59*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/36. + P(3,1)*ProjM(2,1) - P(3,2)*ProjM(2,1) - P(3,4)*ProjM(2,1)')

FFVS6 = Lorentz(name = 'FFVS6',
                spins = [ 2, 2, 3, 1 ],
                structure = '(19*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/126. + (13*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/9. - (7*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/18. - (85*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/63. + P(3,1)*ProjM(2,1) - P(3,2)*ProjM(2,1) - P(3,4)*ProjM(2,1)')

FFVS7 = Lorentz(name = 'FFVS7',
                spins = [ 2, 2, 3, 1 ],
                structure = '(13*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (4*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/11.')

FFVS8 = Lorentz(name = 'FFVS8',
                spins = [ 2, 2, 3, 1 ],
                structure = '(13*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + (8*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (4*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/11.')

FFVS9 = Lorentz(name = 'FFVS9',
                spins = [ 2, 2, 3, 1 ],
                structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3.')

FFVS10 = Lorentz(name = 'FFVS10',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '-(P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1)) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1)')

FFVS11 = Lorentz(name = 'FFVS11',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(-7*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/18. + (121*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/126. - (7*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/18. + (19*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/126. - (163*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/126. + (19*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/126. + P(3,1)*ProjM(2,1) - P(3,2)*ProjM(2,1) + P(3,4)*ProjM(2,1)')

FFVS12 = Lorentz(name = 'FFVS12',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(-5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/12. + (17*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/9. - (5*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/12. + (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3. - (83*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/36. + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3. + P(3,1)*ProjM(2,1) - P(3,2)*ProjM(2,1) + P(3,4)*ProjM(2,1)')

FFVS13 = Lorentz(name = 'FFVS13',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(11*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/13. + (7*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/13. + (11*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/13. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (5*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/13. + P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1)')

FFVS14 = Lorentz(name = 'FFVS14',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3. + P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1)')

FFVS15 = Lorentz(name = 'FFVS15',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'Gamma(3,2,-1)*ProjP(-1,1)')

FFVS16 = Lorentz(name = 'FFVS16',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'Gamma(3,2,-1)*ProjM(-1,1) - (11*Gamma(3,2,-1)*ProjP(-1,1))/3.')

FFVS17 = Lorentz(name = 'FFVS17',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'Gamma(3,2,-1)*ProjM(-1,1) - (11*Gamma(3,2,-1)*ProjP(-1,1))/6.')

FFVS18 = Lorentz(name = 'FFVS18',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'Gamma(3,2,-1)*ProjM(-1,1) + Gamma(3,2,-1)*ProjP(-1,1)')

FFVS19 = Lorentz(name = 'FFVS19',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2) - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1)')

FFVS20 = Lorentz(name = 'FFVS20',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) - P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) - P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS21 = Lorentz(name = 'FFVS21',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS22 = Lorentz(name = 'FFVS22',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + (95*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/36. - (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/12. - (83*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/36. + P(3,1)*ProjP(2,1) - P(3,2)*ProjP(2,1) - P(3,4)*ProjP(2,1)')

FFVS23 = Lorentz(name = 'FFVS23',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2) + (17*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/19. - (36*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/19. + (17*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/19. - (36*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/19.')

FFVS24 = Lorentz(name = 'FFVS24',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(19*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/126. + (13*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/9. - (7*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/18. - (85*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/63. + P(3,1)*ProjP(2,1) - P(3,2)*ProjP(2,1) - P(3,4)*ProjP(2,1)')

FFVS25 = Lorentz(name = 'FFVS25',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) - P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS26 = Lorentz(name = 'FFVS26',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(-19*P(-1,3)*Gamma5(-2,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2))/36. + (19*P(-1,3)*Gamma5(-2,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3))/36. - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) - P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS27 = Lorentz(name = 'FFVS27',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(13*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/11. + (8*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/11. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (4*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11.')

FFVS28 = Lorentz(name = 'FFVS28',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3.')

FFVS29 = Lorentz(name = 'FFVS29',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS30 = Lorentz(name = 'FFVS30',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '-(P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1)) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS31 = Lorentz(name = 'FFVS31',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '-(P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1)) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) - P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS32 = Lorentz(name = 'FFVS32',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(5*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + (11*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/9. + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (13*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/9. - (5*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. - (11*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/9. - P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) - (13*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/9.')

FFVS33 = Lorentz(name = 'FFVS33',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(13*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + (8*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (4*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/11. - P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) - (7*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/11. - P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) - (13*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11. - (5*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11. - (13*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11.')

FFVS34 = Lorentz(name = 'FFVS34',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3. - P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) - (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. - P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) - P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) - (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3. - P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS35 = Lorentz(name = 'FFVS35',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + (P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/2. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/2. - P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) - (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. - (P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/2. - P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) - (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/2.')

FFVS36 = Lorentz(name = 'FFVS36',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(3,1)*Gamma5(2,1) - P(3,2)*Gamma5(2,1) + (5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/42. - (101*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/84. + (7*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/36. + (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/42. + (37*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/28. - (19*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/252. - (5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/42. + (101*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/84. - (7*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/36. - (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/42. - (37*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/28. + (19*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/252.')

FFVS37 = Lorentz(name = 'FFVS37',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(3,1)*Identity(2,1) - P(3,2)*Identity(2,1) - (5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/42. + (101*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/84. - (7*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/36. - (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/42. - (37*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/28. + (19*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/252. - (5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/42. + (101*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/84. - (7*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/36. - (5*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/42. - (37*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/28. + (19*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/252.')

FFVS38 = Lorentz(name = 'FFVS38',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(-7*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/18. + (121*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/126. - (7*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/18. + (19*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/126. - (163*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/126. + (19*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/126. + P(3,1)*ProjP(2,1) - P(3,2)*ProjP(2,1) + P(3,4)*ProjP(2,1)')

FFVS39 = Lorentz(name = 'FFVS39',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(3,1)*Gamma5(2,1) - P(3,2)*Gamma5(2,1) + (4*P(-1,3)*Gamma5(-2,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2))/9. - (4*P(-1,3)*Gamma5(-2,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3))/9. + (P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/24. - (131*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/72. + (5*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/24. + (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/24. + (67*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/36. - (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/6. - (P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/24. + (131*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/72. - (5*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/24. - (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/24. - (67*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/36. + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/6.')

FFVS40 = Lorentz(name = 'FFVS40',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(-4*P(-1,3)*Gamma(-1,-2,1)*Gamma(3,2,-2))/9. + P(3,1)*Identity(2,1) - P(3,2)*Identity(2,1) - (P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/24. + (163*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/72. - (5*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/24. - (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/24. - (67*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/36. + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/6. - (P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/24. + (163*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/72. - (5*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/24. - (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/24. - (67*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/36. + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/6.')

FFVS41 = Lorentz(name = 'FFVS41',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(-5*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/12. + (11*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/9. - (5*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/12. + (P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3. - (59*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/36. + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3. + P(3,1)*ProjP(2,1) - P(3,2)*ProjP(2,1) + P(3,4)*ProjP(2,1)')

FFVS42 = Lorentz(name = 'FFVS42',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + (P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/2. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/2. + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + (P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/2. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/2.')

FFVS43 = Lorentz(name = 'FFVS43',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(11*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/13. + (7*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/13. + (11*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/13. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (5*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/13. + P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS44 = Lorentz(name = 'FFVS44',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3. + P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS45 = Lorentz(name = 'FFVS45',
                 spins = [ 2, 2, 3, 1 ],
                 structure = 'P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/3. + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (2*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/3. + P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1)')

FFVS46 = Lorentz(name = 'FFVS46',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(13*P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + (8*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (4*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/11. + P(-1,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (7*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/11. + P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1) + (13*P(-1,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11. + (5*P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11. + (13*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/11.')

FFVS47 = Lorentz(name = 'FFVS47',
                 spins = [ 2, 2, 3, 1 ],
                 structure = '(5*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/3. + (11*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjM(-2,1))/9. + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1) + (13*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjM(-2,1))/9. + (5*P(-1,3)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/3. + (11*P(-1,4)*Gamma(-1,2,-3)*Gamma(3,-3,-2)*ProjP(-2,1))/9. + P(-1,3)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1) + (13*P(-1,4)*Gamma(-1,-3,-2)*Gamma(3,2,-3)*ProjP(-2,1))/9.')

FFVV1 = Lorentz(name = 'FFVV1',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,-1,1)*Gamma(4,2,-1)')

FFVV2 = Lorentz(name = 'FFVV2',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,2,-1)*Gamma(4,-1,1)')

FFVV3 = Lorentz(name = 'FFVV3',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1)')

FFVV4 = Lorentz(name = 'FFVV4',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Identity(2,1)*Metric(3,4)')

FFVV5 = Lorentz(name = 'FFVV5',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) - (22*Identity(2,1)*Metric(3,4))/17.')

FFVV6 = Lorentz(name = 'FFVV6',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) - (4*Identity(2,1)*Metric(3,4))/11.')

FFVV7 = Lorentz(name = 'FFVV7',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) - (2*Identity(2,1)*Metric(3,4))/7.')

FFVV8 = Lorentz(name = 'FFVV8',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) + (56*Identity(2,1)*Metric(3,4))/65.')

FFVV9 = Lorentz(name = 'FFVV9',
                spins = [ 2, 2, 3, 3 ],
                structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) + (160*Identity(2,1)*Metric(3,4))/127.')

FFVV10 = Lorentz(name = 'FFVV10',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1)')

FFVV11 = Lorentz(name = 'FFVV11',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1)')

FFVV12 = Lorentz(name = 'FFVV12',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1)')

FFVV13 = Lorentz(name = 'FFVV13',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/3. - (4*Metric(3,4)*ProjM(2,1))/3.')

FFVV14 = Lorentz(name = 'FFVV14',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (10*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/19. - (20*Metric(3,4)*ProjM(2,1))/19.')

FFVV15 = Lorentz(name = 'FFVV15',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (5*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/7. - (6*Metric(3,4)*ProjM(2,1))/7.')

FFVV16 = Lorentz(name = 'FFVV16',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (11*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/8. - (Metric(3,4)*ProjM(2,1))/20.')

FFVV17 = Lorentz(name = 'FFVV17',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/5. - (6*Metric(3,4)*ProjM(2,1))/5.')

FFVV18 = Lorentz(name = 'FFVV18',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (19*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/10. - 2*Metric(3,4)*ProjM(2,1)')

FFVV19 = Lorentz(name = 'FFVV19',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (158*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/77. - (4*Metric(3,4)*ProjM(2,1))/11.')

FFVV20 = Lorentz(name = 'FFVV20',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/3. - (10*Metric(3,4)*ProjM(2,1))/3.')

FFVV21 = Lorentz(name = 'FFVV21',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVV22 = Lorentz(name = 'FFVV22',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1)')

FFVV23 = Lorentz(name = 'FFVV23',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVV24 = Lorentz(name = 'FFVV24',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) + 2*Metric(3,4)*ProjP(2,1)')

FFVV25 = Lorentz(name = 'FFVV25',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,-1,1)*Gamma(4,2,-1) - Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1)')

FFVV26 = Lorentz(name = 'FFVV26',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Identity(2,1)*Metric(3,4) - (Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/2. - (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/2. - (Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/2. - (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/2.')

FFVV27 = Lorentz(name = 'FFVV27',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (11*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/8. - (Metric(3,4)*ProjM(2,1))/20. + (5*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/8. + (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/4. + (29*Metric(3,4)*ProjP(2,1))/20.')

FFVV28 = Lorentz(name = 'FFVV28',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (158*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/77. - (4*Metric(3,4)*ProjM(2,1))/11. + (50*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/77. - (31*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/77. + (188*Metric(3,4)*ProjP(2,1))/77.')

FFVV29 = Lorentz(name = 'FFVV29',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + 2*Metric(3,4)*ProjP(2,1)')

FFVV30 = Lorentz(name = 'FFVV30',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - 8*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) + 16*Metric(3,4)*ProjP(2,1)')

FFVV31 = Lorentz(name = 'FFVV31',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVV32 = Lorentz(name = 'FFVV32',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) + Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVV33 = Lorentz(name = 'FFVV33',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - (31*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/50. + (94*Metric(3,4)*ProjP(2,1))/25.')

FFVV34 = Lorentz(name = 'FFVV34',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/8. - 2*Metric(3,4)*ProjP(2,1)')

FFVV35 = Lorentz(name = 'FFVV35',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/3. - (4*Metric(3,4)*ProjP(2,1))/3.')

FFVV36 = Lorentz(name = 'FFVV36',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (2*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/5. + (58*Metric(3,4)*ProjP(2,1))/25.')

FFVV37 = Lorentz(name = 'FFVV37',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) + (56*Metric(3,4)*ProjM(2,1))/65. + Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) + (56*Metric(3,4)*ProjP(2,1))/65.')

FFVV38 = Lorentz(name = 'FFVV38',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) + (160*Metric(3,4)*ProjM(2,1))/127. + Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) + (160*Metric(3,4)*ProjP(2,1))/127.')

FFVV39 = Lorentz(name = 'FFVV39',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/3. - (10*Metric(3,4)*ProjP(2,1))/3.')

FFVV40 = Lorentz(name = 'FFVV40',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/9. + (2*Metric(3,4)*ProjM(2,1))/9. + (14*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/9. + (4*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/3. - (8*Metric(3,4)*ProjP(2,1))/9.')

FFVV41 = Lorentz(name = 'FFVV41',
                 spins = [ 2, 2, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + 4*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) + 4*Metric(3,4)*ProjM(2,1) + 7*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + 10*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) - 8*Metric(3,4)*ProjP(2,1)')

VSSS1 = Lorentz(name = 'VSSS1',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2)')

VSSS2 = Lorentz(name = 'VSSS2',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2) - P(1,3)')

VSSS3 = Lorentz(name = 'VSSS3',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2) + P(1,3)')

VSSS4 = Lorentz(name = 'VSSS4',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2) + P(1,3) - 2*P(1,4)')

VSSS5 = Lorentz(name = 'VSSS5',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2) - P(1,4)/2.')

VSSS6 = Lorentz(name = 'VSSS6',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2) - P(1,3)/2. - P(1,4)/2.')

VSSS7 = Lorentz(name = 'VSSS7',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2) - P(1,3)/3. - P(1,4)/3.')

VSSS8 = Lorentz(name = 'VSSS8',
                spins = [ 3, 1, 1, 1 ],
                structure = 'P(1,2) + P(1,3) + P(1,4)')

VVSS1 = Lorentz(name = 'VVSS1',
                spins = [ 3, 3, 1, 1 ],
                structure = 'P(1,2)*P(2,1)')

VVSS2 = Lorentz(name = 'VVSS2',
                spins = [ 3, 3, 1, 1 ],
                structure = 'Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,2) - Epsilon(1,2,-1,-2)*P(-2,4)*P(-1,2) - 2*Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,4)')

VVSS3 = Lorentz(name = 'VVSS3',
                spins = [ 3, 3, 1, 1 ],
                structure = 'Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,1) - Epsilon(1,2,-1,-2)*P(-2,4)*P(-1,1) - Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,2) + Epsilon(1,2,-1,-2)*P(-2,4)*P(-1,2)')

VVSS4 = Lorentz(name = 'VVSS4',
                spins = [ 3, 3, 1, 1 ],
                structure = 'Metric(1,2)')

VVSS5 = Lorentz(name = 'VVSS5',
                spins = [ 3, 3, 1, 1 ],
                structure = 'P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)')

VVSS6 = Lorentz(name = 'VVSS6',
                spins = [ 3, 3, 1, 1 ],
                structure = 'P(1,3)*P(2,1) + P(1,4)*P(2,1) + P(1,2)*P(2,3) + P(1,2)*P(2,4) + 2*P(-1,1)*P(-1,2)*Metric(1,2) - 2*P(-1,2)**2*Metric(1,2) - P(-1,1)*P(-1,3)*Metric(1,2) - 3*P(-1,2)*P(-1,3)*Metric(1,2) - P(-1,1)*P(-1,4)*Metric(1,2) - 3*P(-1,2)*P(-1,4)*Metric(1,2)')

VVSS7 = Lorentz(name = 'VVSS7',
                spins = [ 3, 3, 1, 1 ],
                structure = 'Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,1) - Epsilon(1,2,-1,-2)*P(-2,4)*P(-1,1) - Epsilon(1,2,-1,-2)*P(-2,3)*P(-1,2) + Epsilon(1,2,-1,-2)*P(-2,4)*P(-1,2) + (complex(0,1)*P(1,3)*P(2,1))/2. + (complex(0,1)*P(1,4)*P(2,1))/2. + (complex(0,1)*P(1,2)*P(2,3))/2. + (complex(0,1)*P(1,2)*P(2,4))/2. + complex(0,1)*P(-1,1)*P(-1,2)*Metric(1,2) - complex(0,1)*P(-1,2)**2*Metric(1,2) - (complex(0,1)*P(-1,1)*P(-1,3)*Metric(1,2))/2. - (3*complex(0,1)*P(-1,2)*P(-1,3)*Metric(1,2))/2. - (complex(0,1)*P(-1,1)*P(-1,4)*Metric(1,2))/2. - (3*complex(0,1)*P(-1,2)*P(-1,4)*Metric(1,2))/2.')

VVSS8 = Lorentz(name = 'VVSS8',
                spins = [ 3, 3, 1, 1 ],
                structure = 'P(1,2)*P(2,1) - (P(-1,1)*P(-1,2)*Metric(1,2))/2. + (P(-1,2)**2*Metric(1,2))/2. + (P(-1,2)*P(-1,3)*Metric(1,2))/2. + (P(-1,2)*P(-1,4)*Metric(1,2))/2.')

VVSS9 = Lorentz(name = 'VVSS9',
                spins = [ 3, 3, 1, 1 ],
                structure = 'P(1,2)*P(2,1) - (9*P(1,3)*P(2,1))/25. - (9*P(1,4)*P(2,1))/25. - (9*P(1,2)*P(2,3))/25. - (9*P(1,2)*P(2,4))/25. - (16*P(-1,1)*P(-1,2)*Metric(1,2))/25. + (9*P(-1,2)**2*Metric(1,2))/25. + (3*P(-1,1)*P(-1,3)*Metric(1,2))/25. + (12*P(-1,2)*P(-1,3)*Metric(1,2))/25. - (6*P(-1,3)**2*Metric(1,2))/25. + (3*P(-1,1)*P(-1,4)*Metric(1,2))/25. + (12*P(-1,2)*P(-1,4)*Metric(1,2))/25. - (12*P(-1,3)*P(-1,4)*Metric(1,2))/25. - (6*P(-1,4)**2*Metric(1,2))/25.')

VVSS10 = Lorentz(name = 'VVSS10',
                 spins = [ 3, 3, 1, 1 ],
                 structure = 'P(1,2)*P(2,1) - (13*P(1,3)*P(2,1))/62. - (13*P(1,4)*P(2,1))/62. - (13*P(1,2)*P(2,3))/62. + (P(1,3)*P(2,3))/62. + (P(1,4)*P(2,3))/62. - (13*P(1,2)*P(2,4))/62. + (P(1,3)*P(2,4))/62. + (P(1,4)*P(2,4))/62. - (20*P(-1,1)*P(-1,2)*Metric(1,2))/31. + (19*P(-1,2)**2*Metric(1,2))/62. + (2*P(-1,1)*P(-1,3)*Metric(1,2))/31. + (23*P(-1,2)*P(-1,3)*Metric(1,2))/62. - (13*P(-1,3)**2*Metric(1,2))/62. + (2*P(-1,1)*P(-1,4)*Metric(1,2))/31. + (23*P(-1,2)*P(-1,4)*Metric(1,2))/62. - (13*P(-1,3)*P(-1,4)*Metric(1,2))/31. - (13*P(-1,4)**2*Metric(1,2))/62.')

VVVS1 = Lorentz(name = 'VVVS1',
                spins = [ 3, 3, 3, 1 ],
                structure = '-(Epsilon(1,2,3,-1)*P(-1,1)) + Epsilon(1,2,3,-1)*P(-1,2)')

VVVS2 = Lorentz(name = 'VVVS2',
                spins = [ 3, 3, 3, 1 ],
                structure = '-(Epsilon(1,2,3,-1)*P(-1,1)) - 2*Epsilon(1,2,3,-1)*P(-1,2)')

VVVS3 = Lorentz(name = 'VVVS3',
                spins = [ 3, 3, 3, 1 ],
                structure = '-(Epsilon(1,2,3,-1)*P(-1,4))')

VVVS4 = Lorentz(name = 'VVVS4',
                spins = [ 3, 3, 3, 1 ],
                structure = '-(Epsilon(1,2,3,-1)*P(-1,1)) - Epsilon(1,2,3,-1)*P(-1,2) - Epsilon(1,2,3,-1)*P(-1,3) - 10*Epsilon(1,2,3,-1)*P(-1,4)')

VVVS5 = Lorentz(name = 'VVVS5',
                spins = [ 3, 3, 3, 1 ],
                structure = 'P(3,1)*Metric(1,2) - P(2,1)*Metric(1,3)')

VVVS6 = Lorentz(name = 'VVVS6',
                spins = [ 3, 3, 3, 1 ],
                structure = 'P(3,1)*Metric(1,2) + P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) - P(1,2)*Metric(2,3)')

VVVS7 = Lorentz(name = 'VVVS7',
                spins = [ 3, 3, 3, 1 ],
                structure = 'P(2,3)*Metric(1,3) - P(1,3)*Metric(2,3)')

VVVS8 = Lorentz(name = 'VVVS8',
                spins = [ 3, 3, 3, 1 ],
                structure = 'P(3,2)*Metric(1,2) + P(2,3)*Metric(1,3) - P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

VVVS9 = Lorentz(name = 'VVVS9',
                spins = [ 3, 3, 3, 1 ],
                structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

VVVS10 = Lorentz(name = 'VVVS10',
                 spins = [ 3, 3, 3, 1 ],
                 structure = 'P(3,4)*Metric(1,2) + P(2,4)*Metric(1,3) + P(1,4)*Metric(2,3)')

VVVV1 = Lorentz(name = 'VVVV1',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Epsilon(1,2,3,4)')

VVVV2 = Lorentz(name = 'VVVV2',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3)')

VVVV3 = Lorentz(name = 'VVVV3',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,3)*Metric(2,4)')

VVVV4 = Lorentz(name = 'VVVV4',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)')

VVVV5 = Lorentz(name = 'VVVV5',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,2)*Metric(3,4)')

VVVV6 = Lorentz(name = 'VVVV6',
                spins = [ 3, 3, 3, 3 ],
                structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,1)*P(4,2)*Metric(1,2) + P(2,1)*P(4,2)*Metric(1,3) + P(2,4)*P(4,3)*Metric(1,3) - P(2,1)*P(3,2)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,2)*P(3,1)*Metric(2,4) + P(1,3)*P(3,4)*Metric(2,4) - P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) - P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,3)*Metric(3,4) - P(1,3)*P(2,4)*Metric(3,4)')

VVVV7 = Lorentz(name = 'VVVV7',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)')

VVVV8 = Lorentz(name = 'VVVV8',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)')

VVVV9 = Lorentz(name = 'VVVV9',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)')

VVVV10 = Lorentz(name = 'VVVV10',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'Metric(1,4)*Metric(2,3) - (Metric(1,3)*Metric(2,4))/2. - (Metric(1,2)*Metric(3,4))/2.')

VVVV11 = Lorentz(name = 'VVVV11',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) + Metric(1,2)*Metric(3,4)')

VVVV12 = Lorentz(name = 'VVVV12',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(3,2)*P(4,3)*Metric(1,2) - P(2,4)*P(4,1)*Metric(1,3) - P(2,3)*P(4,2)*Metric(1,3) + P(2,4)*P(3,1)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) + P(1,3)*P(4,2)*Metric(2,3) - P(1,2)*P(4,3)*Metric(2,3) - P(1,4)*P(3,1)*Metric(2,4) - P(1,3)*P(3,2)*Metric(2,4) + P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,4)*P(2,1)*Metric(3,4) + P(1,2)*P(2,3)*Metric(3,4) - P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV13 = Lorentz(name = 'VVVV13',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,2)*Metric(1,2) + P(3,1)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,1)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) - P(2,4)*P(3,2)*Metric(1,4) - P(1,3)*P(4,1)*Metric(2,3) - P(1,4)*P(4,2)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) + P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,4)*P(3,2)*Metric(2,4) - P(1,2)*P(3,4)*Metric(2,4) + P(1,3)*P(2,1)*Metric(3,4) + P(1,2)*P(2,4)*Metric(3,4) - P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV14 = Lorentz(name = 'VVVV14',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,2)*P(4,1)*Metric(1,2) - P(3,1)*P(4,2)*Metric(1,2) + P(3,4)*P(4,2)*Metric(1,2) + P(3,1)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) + P(2,1)*P(4,2)*Metric(1,3) - P(2,1)*P(4,3)*Metric(1,3) + P(2,4)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) - P(2,1)*P(3,2)*Metric(1,4) - P(2,4)*P(3,2)*Metric(1,4) - P(2,3)*P(3,4)*Metric(1,4) - P(1,2)*P(4,1)*Metric(2,3) - P(1,3)*P(4,1)*Metric(2,3) - P(1,4)*P(4,2)*Metric(2,3) - P(1,4)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,2)*Metric(1,4)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) + P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(-1,3)*P(-1,4)*Metric(1,4)*Metric(2,3) + P(1,2)*P(3,1)*Metric(2,4) + P(1,4)*P(3,2)*Metric(2,4) - P(1,2)*P(3,4)*Metric(2,4) + P(1,3)*P(3,4)*Metric(2,4) - P(-1,1)*P(-1,2)*Metric(1,3)*Metric(2,4) - P(-1,3)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,3)*P(2,1)*Metric(3,4) + P(1,4)*P(2,3)*Metric(3,4) + P(1,2)*P(2,4)*Metric(3,4) - P(1,3)*P(2,4)*Metric(3,4) - P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

VVVV15 = Lorentz(name = 'VVVV15',
                 spins = [ 3, 3, 3, 3 ],
                 structure = 'P(3,4)*P(4,1)*Metric(1,2) + P(3,4)*P(4,2)*Metric(1,2) + P(3,1)*P(4,3)*Metric(1,2) + P(3,2)*P(4,3)*Metric(1,2) + P(2,3)*P(4,1)*Metric(1,3) - P(2,4)*P(4,1)*Metric(1,3) - P(2,3)*P(4,2)*Metric(1,3) - P(2,1)*P(4,3)*Metric(1,3) - P(2,3)*P(3,1)*Metric(1,4) + P(2,4)*P(3,1)*Metric(1,4) - P(2,4)*P(3,2)*Metric(1,4) - P(2,1)*P(3,4)*Metric(1,4) - P(1,3)*P(4,1)*Metric(2,3) + P(1,3)*P(4,2)*Metric(2,3) - P(1,4)*P(4,2)*Metric(2,3) - P(1,2)*P(4,3)*Metric(2,3) + P(-1,1)*P(-1,3)*Metric(1,4)*Metric(2,3) + P(-1,2)*P(-1,4)*Metric(1,4)*Metric(2,3) - P(1,4)*P(3,1)*Metric(2,4) - P(1,3)*P(3,2)*Metric(2,4) + P(1,4)*P(3,2)*Metric(2,4) - P(1,2)*P(3,4)*Metric(2,4) + P(-1,2)*P(-1,3)*Metric(1,3)*Metric(2,4) + P(-1,1)*P(-1,4)*Metric(1,3)*Metric(2,4) + P(1,3)*P(2,1)*Metric(3,4) + P(1,4)*P(2,1)*Metric(3,4) + P(1,2)*P(2,3)*Metric(3,4) + P(1,2)*P(2,4)*Metric(3,4) - P(-1,1)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,3)*Metric(1,2)*Metric(3,4) - P(-1,1)*P(-1,4)*Metric(1,2)*Metric(3,4) - P(-1,2)*P(-1,4)*Metric(1,2)*Metric(3,4)')

SSSSS1 = Lorentz(name = 'SSSSS1',
                 spins = [ 1, 1, 1, 1, 1 ],
                 structure = '1')

FFSSS1 = Lorentz(name = 'FFSSS1',
                 spins = [ 2, 2, 1, 1, 1 ],
                 structure = 'Gamma5(2,1)')

FFSSS2 = Lorentz(name = 'FFSSS2',
                 spins = [ 2, 2, 1, 1, 1 ],
                 structure = 'Identity(2,1)')

FFSSS3 = Lorentz(name = 'FFSSS3',
                 spins = [ 2, 2, 1, 1, 1 ],
                 structure = 'ProjM(2,1)')

FFSSS4 = Lorentz(name = 'FFSSS4',
                 spins = [ 2, 2, 1, 1, 1 ],
                 structure = 'ProjP(2,1)')

FFVSS1 = Lorentz(name = 'FFVSS1',
                 spins = [ 2, 2, 3, 1, 1 ],
                 structure = 'Gamma(3,2,-1)*ProjM(-1,1)')

FFVSS2 = Lorentz(name = 'FFVSS2',
                 spins = [ 2, 2, 3, 1, 1 ],
                 structure = 'Gamma(3,2,-1)*ProjP(-1,1)')

FFVVS1 = Lorentz(name = 'FFVVS1',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma(3,-1,1)*Gamma(4,2,-1)')

FFVVS2 = Lorentz(name = 'FFVVS2',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2)')

FFVVS3 = Lorentz(name = 'FFVVS3',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma(3,2,-1)*Gamma(4,-1,1)')

FFVVS4 = Lorentz(name = 'FFVVS4',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1)')

FFVVS5 = Lorentz(name = 'FFVVS5',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1)')

FFVVS6 = Lorentz(name = 'FFVVS6',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1) - Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2)')

FFVVS7 = Lorentz(name = 'FFVVS7',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1) + Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2)')

FFVVS8 = Lorentz(name = 'FFVVS8',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Gamma5(2,1)*Metric(3,4)')

FFVVS9 = Lorentz(name = 'FFVVS9',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'Identity(2,1)*Metric(3,4)')

FFVVS10 = Lorentz(name = 'FFVVS10',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1) + Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2) - 2*Gamma5(2,1)*Metric(3,4)')

FFVVS11 = Lorentz(name = 'FFVVS11',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1) + Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2) - (4*Gamma5(2,1)*Metric(3,4))/11.')

FFVVS12 = Lorentz(name = 'FFVVS12',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1) + Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2) - (2*Gamma5(2,1)*Metric(3,4))/7.')

FFVVS13 = Lorentz(name = 'FFVVS13',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1) + Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2) - (56*Gamma5(2,1)*Metric(3,4))/235.')

FFVVS14 = Lorentz(name = 'FFVVS14',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(-1,1)*Gamma(3,2,-2)*Gamma(4,-2,-1) + Gamma5(-1,1)*Gamma(3,-2,-1)*Gamma(4,2,-2) - (4*Gamma5(2,1)*Metric(3,4))/95.')

FFVVS15 = Lorentz(name = 'FFVVS15',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) - (4*Identity(2,1)*Metric(3,4))/11.')

FFVVS16 = Lorentz(name = 'FFVVS16',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) - (2*Identity(2,1)*Metric(3,4))/7.')

FFVVS17 = Lorentz(name = 'FFVVS17',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) - (56*Identity(2,1)*Metric(3,4))/235.')

FFVVS18 = Lorentz(name = 'FFVVS18',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) + Gamma(3,-1,1)*Gamma(4,2,-1) - (4*Identity(2,1)*Metric(3,4))/95.')

FFVVS19 = Lorentz(name = 'FFVVS19',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Metric(3,4)*ProjM(2,1)')

FFVVS20 = Lorentz(name = 'FFVVS20',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1)')

FFVVS21 = Lorentz(name = 'FFVVS21',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1)')

FFVVS22 = Lorentz(name = 'FFVVS22',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1)')

FFVVS23 = Lorentz(name = 'FFVVS23',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/3. - (4*Metric(3,4)*ProjM(2,1))/3.')

FFVVS24 = Lorentz(name = 'FFVVS24',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (3*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/7. - (10*Metric(3,4)*ProjM(2,1))/7.')

FFVVS25 = Lorentz(name = 'FFVVS25',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (77*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/158. - (14*Metric(3,4)*ProjM(2,1))/79.')

FFVVS26 = Lorentz(name = 'FFVVS26',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (10*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/19. - (20*Metric(3,4)*ProjM(2,1))/19.')

FFVVS27 = Lorentz(name = 'FFVVS27',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/10. - (4*Metric(3,4)*ProjM(2,1))/5.')

FFVVS28 = Lorentz(name = 'FFVVS28',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (5*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/7. - (6*Metric(3,4)*ProjM(2,1))/7.')

FFVVS29 = Lorentz(name = 'FFVVS29',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (8*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/11. - (2*Metric(3,4)*ProjM(2,1))/55.')

FFVVS30 = Lorentz(name = 'FFVVS30',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) - 2*Metric(3,4)*ProjM(2,1)')

FFVVS31 = Lorentz(name = 'FFVVS31',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) - (4*Metric(3,4)*ProjM(2,1))/11.')

FFVVS32 = Lorentz(name = 'FFVVS32',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) - (2*Metric(3,4)*ProjM(2,1))/7.')

FFVVS33 = Lorentz(name = 'FFVVS33',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/6. - (2*Metric(3,4)*ProjM(2,1))/3.')

FFVVS34 = Lorentz(name = 'FFVVS34',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (11*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/8. - (Metric(3,4)*ProjM(2,1))/20.')

FFVVS35 = Lorentz(name = 'FFVVS35',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/5. - (6*Metric(3,4)*ProjM(2,1))/5.')

FFVVS36 = Lorentz(name = 'FFVVS36',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (19*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/10. - 2*Metric(3,4)*ProjM(2,1)')

FFVVS37 = Lorentz(name = 'FFVVS37',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (158*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/77. - (4*Metric(3,4)*ProjM(2,1))/11.')

FFVVS38 = Lorentz(name = 'FFVVS38',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/3. - (10*Metric(3,4)*ProjM(2,1))/3.')

FFVVS39 = Lorentz(name = 'FFVVS39',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + 3*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) - 4*Metric(3,4)*ProjM(2,1)')

FFVVS40 = Lorentz(name = 'FFVVS40',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Metric(3,4)*ProjP(2,1)')

FFVVS41 = Lorentz(name = 'FFVVS41',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVVS42 = Lorentz(name = 'FFVVS42',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1)')

FFVVS43 = Lorentz(name = 'FFVVS43',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-1)*Gamma(4,-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVVS44 = Lorentz(name = 'FFVVS44',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) - 2*Metric(3,4)*ProjP(2,1)')

FFVVS45 = Lorentz(name = 'FFVVS45',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) + 2*Metric(3,4)*ProjP(2,1)')

FFVVS46 = Lorentz(name = 'FFVVS46',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Identity(2,1)*Metric(3,4) - 20*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - (55*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/2. - (55*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/2. - 20*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVVS47 = Lorentz(name = 'FFVVS47',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(2,1)*Metric(3,4) + 20*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (55*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/2. - (55*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/2. - 20*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVVS48 = Lorentz(name = 'FFVVS48',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Identity(2,1)*Metric(3,4) - (95*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/4. - (95*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/4. - (95*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/4. - (95*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/4.')

FFVVS49 = Lorentz(name = 'FFVVS49',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(2,1)*Metric(3,4) + (95*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/4. + (95*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/4. - (95*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/4. - (95*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/4.')

FFVVS50 = Lorentz(name = 'FFVVS50',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + 4*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) + 4*Metric(3,4)*ProjM(2,1) - 7*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - 10*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) + 8*Metric(3,4)*ProjP(2,1)')

FFVVS51 = Lorentz(name = 'FFVVS51',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Identity(2,1)*Metric(3,4) - (11*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/4. - (79*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/14. - (79*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/14. - (11*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/4.')

FFVVS52 = Lorentz(name = 'FFVVS52',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(2,1)*Metric(3,4) + (11*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/4. + (79*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/14. - (79*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/14. - (11*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/4.')

FFVVS53 = Lorentz(name = 'FFVVS53',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Identity(2,1)*Metric(3,4) - (235*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/56. - (235*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/56. - (235*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/56. - (235*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/56.')

FFVVS54 = Lorentz(name = 'FFVVS54',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma5(2,1)*Metric(3,4) + (235*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/56. + (235*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/56. - (235*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/56. - (235*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/56.')

FFVVS55 = Lorentz(name = 'FFVVS55',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/9. + (2*Metric(3,4)*ProjM(2,1))/9. - (14*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/9. - (4*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/3. + (8*Metric(3,4)*ProjP(2,1))/9.')

FFVVS56 = Lorentz(name = 'FFVVS56',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,-1,1)*Gamma(4,2,-1) - Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1)')

FFVVS57 = Lorentz(name = 'FFVVS57',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) - Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVVS58 = Lorentz(name = 'FFVVS58',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Identity(2,1)*Metric(3,4) - (Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1))/2. - (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/2. - (Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/2. - (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/2.')

FFVVS59 = Lorentz(name = 'FFVVS59',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - 2*Metric(3,4)*ProjP(2,1)')

FFVVS60 = Lorentz(name = 'FFVVS60',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + 2*Metric(3,4)*ProjP(2,1)')

FFVVS61 = Lorentz(name = 'FFVVS61',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - 8*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) + 16*Metric(3,4)*ProjP(2,1)')

FFVVS62 = Lorentz(name = 'FFVVS62',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVVS63 = Lorentz(name = 'FFVVS63',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) + Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)')

FFVVS64 = Lorentz(name = 'FFVVS64',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) - (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/8. - 2*Metric(3,4)*ProjP(2,1)')

FFVVS65 = Lorentz(name = 'FFVVS65',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/4. + Metric(3,4)*ProjP(2,1)')

FFVVS66 = Lorentz(name = 'FFVVS66',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/3. - (4*Metric(3,4)*ProjP(2,1))/3.')

FFVVS67 = Lorentz(name = 'FFVVS67',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (3*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/7. - (10*Metric(3,4)*ProjP(2,1))/7.')

FFVVS68 = Lorentz(name = 'FFVVS68',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (77*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/158. - (14*Metric(3,4)*ProjP(2,1))/79.')

FFVVS69 = Lorentz(name = 'FFVVS69',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (8*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/11. - (2*Metric(3,4)*ProjP(2,1))/55.')

FFVVS70 = Lorentz(name = 'FFVVS70',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) - 2*Metric(3,4)*ProjP(2,1)')

FFVVS71 = Lorentz(name = 'FFVVS71',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) - (4*Metric(3,4)*ProjP(2,1))/11.')

FFVVS72 = Lorentz(name = 'FFVVS72',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) - (2*Metric(3,4)*ProjP(2,1))/7.')

FFVVS73 = Lorentz(name = 'FFVVS73',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (9*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/7. + (2*Metric(3,4)*ProjP(2,1))/7.')

FFVVS74 = Lorentz(name = 'FFVVS74',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (11*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/8. - (Metric(3,4)*ProjP(2,1))/20.')

FFVVS75 = Lorentz(name = 'FFVVS75',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (158*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/77. - (4*Metric(3,4)*ProjP(2,1))/11.')

FFVVS76 = Lorentz(name = 'FFVVS76',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/3. - (10*Metric(3,4)*ProjP(2,1))/3.')

FFVVS77 = Lorentz(name = 'FFVVS77',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + 3*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) - 4*Metric(3,4)*ProjP(2,1)')

FFVVS78 = Lorentz(name = 'FFVVS78',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + (7*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1))/9. + (2*Metric(3,4)*ProjM(2,1))/9. + (14*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1))/9. + (4*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1))/3. - (8*Metric(3,4)*ProjP(2,1))/9.')

FFVVS79 = Lorentz(name = 'FFVVS79',
                  spins = [ 2, 2, 3, 3, 1 ],
                  structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) + 4*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) + 4*Metric(3,4)*ProjM(2,1) + 7*Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) + 10*Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1) - 8*Metric(3,4)*ProjP(2,1)')

FFVVV1 = Lorentz(name = 'FFVVV1',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,-2,-1)*Gamma(4,-1,1)*Gamma(5,2,-2)')

FFVVV2 = Lorentz(name = 'FFVVV2',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,-1,1)*Gamma(4,-2,-1)*Gamma(5,2,-2)')

FFVVV3 = Lorentz(name = 'FFVVV3',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,-2,-1)*Gamma(4,2,-2)*Gamma(5,-1,1)')

FFVVV4 = Lorentz(name = 'FFVVV4',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-2,-1)*Gamma(5,-1,1)')

FFVVV5 = Lorentz(name = 'FFVVV5',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,-1,1)*Gamma(4,2,-2)*Gamma(5,-2,-1)')

FFVVV6 = Lorentz(name = 'FFVVV6',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,2,-2)*Gamma(4,-1,1)*Gamma(5,-2,-1)')

FFVVV7 = Lorentz(name = 'FFVVV7',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,-3,-2)*Gamma(4,-2,-1)*Gamma(5,2,-3)*ProjM(-1,1)')

FFVVV8 = Lorentz(name = 'FFVVV8',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,-2,-1)*Gamma(4,-3,-2)*Gamma(5,2,-3)*ProjM(-1,1)')

FFVVV9 = Lorentz(name = 'FFVVV9',
                 spins = [ 2, 2, 3, 3, 3 ],
                 structure = 'Gamma(3,-3,-2)*Gamma(4,2,-3)*Gamma(5,-2,-1)*ProjM(-1,1)')

FFVVV10 = Lorentz(name = 'FFVVV10',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,2,-3)*Gamma(4,-3,-2)*Gamma(5,-2,-1)*ProjM(-1,1)')

FFVVV11 = Lorentz(name = 'FFVVV11',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,-2,-1)*Gamma(4,2,-3)*Gamma(5,-3,-2)*ProjM(-1,1)')

FFVVV12 = Lorentz(name = 'FFVVV12',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,2,-3)*Gamma(4,-2,-1)*Gamma(5,-3,-2)*ProjM(-1,1)')

FFVVV13 = Lorentz(name = 'FFVVV13',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,-3,-2)*Gamma(4,-2,-1)*Gamma(5,2,-3)*ProjP(-1,1)')

FFVVV14 = Lorentz(name = 'FFVVV14',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,-2,-1)*Gamma(4,-3,-2)*Gamma(5,2,-3)*ProjP(-1,1)')

FFVVV15 = Lorentz(name = 'FFVVV15',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,-3,-2)*Gamma(4,2,-3)*Gamma(5,-2,-1)*ProjP(-1,1)')

FFVVV16 = Lorentz(name = 'FFVVV16',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,2,-3)*Gamma(4,-3,-2)*Gamma(5,-2,-1)*ProjP(-1,1)')

FFVVV17 = Lorentz(name = 'FFVVV17',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,-2,-1)*Gamma(4,2,-3)*Gamma(5,-3,-2)*ProjP(-1,1)')

FFVVV18 = Lorentz(name = 'FFVVV18',
                  spins = [ 2, 2, 3, 3, 3 ],
                  structure = 'Gamma(3,2,-3)*Gamma(4,-2,-1)*Gamma(5,-3,-2)*ProjP(-1,1)')

VSSSS1 = Lorentz(name = 'VSSSS1',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,2) - P(1,3)')

VSSSS2 = Lorentz(name = 'VSSSS2',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,2) + P(1,3) - 2*P(1,4)')

VSSSS3 = Lorentz(name = 'VSSSS3',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,2) - P(1,3)/2. - P(1,4)/2.')

VSSSS4 = Lorentz(name = 'VSSSS4',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,2) + P(1,3) + P(1,4) - 3*P(1,5)')

VSSSS5 = Lorentz(name = 'VSSSS5',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,2) + P(1,3) - 2*P(1,5)')

VSSSS6 = Lorentz(name = 'VSSSS6',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,3) + P(1,4) - 2*P(1,5)')

VSSSS7 = Lorentz(name = 'VSSSS7',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,2) - P(1,5)')

VSSSS8 = Lorentz(name = 'VSSSS8',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,2) + P(1,3) - P(1,4) - P(1,5)')

VSSSS9 = Lorentz(name = 'VSSSS9',
                 spins = [ 3, 1, 1, 1, 1 ],
                 structure = 'P(1,4) - P(1,5)')

VSSSS10 = Lorentz(name = 'VSSSS10',
                  spins = [ 3, 1, 1, 1, 1 ],
                  structure = 'P(1,2) - P(1,4)/2. - P(1,5)/2.')

VSSSS11 = Lorentz(name = 'VSSSS11',
                  spins = [ 3, 1, 1, 1, 1 ],
                  structure = 'P(1,3) - P(1,4)/2. - P(1,5)/2.')

VSSSS12 = Lorentz(name = 'VSSSS12',
                  spins = [ 3, 1, 1, 1, 1 ],
                  structure = 'P(1,2) - P(1,3)/3. - P(1,4)/3. - P(1,5)/3.')

VVSSS1 = Lorentz(name = 'VVSSS1',
                 spins = [ 3, 3, 1, 1, 1 ],
                 structure = 'Metric(1,2)')

VVVSS1 = Lorentz(name = 'VVVSS1',
                 spins = [ 3, 3, 3, 1, 1 ],
                 structure = 'P(3,1)*Metric(1,2) - P(2,1)*Metric(1,3)')

VVVSS2 = Lorentz(name = 'VVVSS2',
                 spins = [ 3, 3, 3, 1, 1 ],
                 structure = 'P(3,1)*Metric(1,2) + P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) - P(1,2)*Metric(2,3)')

VVVSS3 = Lorentz(name = 'VVVSS3',
                 spins = [ 3, 3, 3, 1, 1 ],
                 structure = 'P(2,3)*Metric(1,3) - P(1,3)*Metric(2,3)')

VVVSS4 = Lorentz(name = 'VVVSS4',
                 spins = [ 3, 3, 3, 1, 1 ],
                 structure = 'P(3,2)*Metric(1,2) + P(2,3)*Metric(1,3) - P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

VVVSS5 = Lorentz(name = 'VVVSS5',
                 spins = [ 3, 3, 3, 1, 1 ],
                 structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

VVVVS1 = Lorentz(name = 'VVVVS1',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Epsilon(1,2,3,4)')

VVVVS2 = Lorentz(name = 'VVVVS2',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,4)*Metric(2,3)')

VVVVS3 = Lorentz(name = 'VVVVS3',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,3)*Metric(2,4)')

VVVVS4 = Lorentz(name = 'VVVVS4',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)')

VVVVS5 = Lorentz(name = 'VVVVS5',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,2)*Metric(3,4)')

VVVVS6 = Lorentz(name = 'VVVVS6',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)')

VVVVS7 = Lorentz(name = 'VVVVS7',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)')

VVVVS8 = Lorentz(name = 'VVVVS8',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)')

VVVVS9 = Lorentz(name = 'VVVVS9',
                 spins = [ 3, 3, 3, 3, 1 ],
                 structure = 'Metric(1,4)*Metric(2,3) - (Metric(1,3)*Metric(2,4))/2. - (Metric(1,2)*Metric(3,4))/2.')

VVVVS10 = Lorentz(name = 'VVVVS10',
                  spins = [ 3, 3, 3, 3, 1 ],
                  structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) + Metric(1,2)*Metric(3,4)')

VVVVV1 = Lorentz(name = 'VVVVV1',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) - P(4,1)*Metric(1,3)*Metric(2,5) + P(3,1)*Metric(1,4)*Metric(2,5) - P(5,1)*Metric(1,2)*Metric(3,4) + P(2,1)*Metric(1,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5)')

VVVVV2 = Lorentz(name = 'VVVVV2',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,2)*Metric(1,4)*Metric(2,3) - P(4,2)*Metric(1,5)*Metric(2,3) + P(3,2)*Metric(1,5)*Metric(2,4) - P(3,2)*Metric(1,4)*Metric(2,5) - P(5,2)*Metric(1,2)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) + P(4,2)*Metric(1,2)*Metric(3,5) - P(1,2)*Metric(2,4)*Metric(3,5)')

VVVVV3 = Lorentz(name = 'VVVVV3',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,4)*Metric(2,3) - (P(5,2)*Metric(1,4)*Metric(2,3))/2. - (P(5,3)*Metric(1,4)*Metric(2,3))/2. - P(4,1)*Metric(1,5)*Metric(2,3) + (P(4,2)*Metric(1,5)*Metric(2,3))/2. + (P(4,3)*Metric(1,5)*Metric(2,3))/2. - (P(5,1)*Metric(1,3)*Metric(2,4))/2. + P(5,2)*Metric(1,3)*Metric(2,4) - (P(5,3)*Metric(1,3)*Metric(2,4))/2. + (P(3,1)*Metric(1,5)*Metric(2,4))/2. - (P(3,2)*Metric(1,5)*Metric(2,4))/2. + (P(4,1)*Metric(1,3)*Metric(2,5))/2. - P(4,2)*Metric(1,3)*Metric(2,5) + (P(4,3)*Metric(1,3)*Metric(2,5))/2. - (P(3,1)*Metric(1,4)*Metric(2,5))/2. + (P(3,2)*Metric(1,4)*Metric(2,5))/2. - (P(5,1)*Metric(1,2)*Metric(3,4))/2. - (P(5,2)*Metric(1,2)*Metric(3,4))/2. + P(5,3)*Metric(1,2)*Metric(3,4) + (P(2,1)*Metric(1,5)*Metric(3,4))/2. - (P(2,3)*Metric(1,5)*Metric(3,4))/2. + (P(1,2)*Metric(2,5)*Metric(3,4))/2. - (P(1,3)*Metric(2,5)*Metric(3,4))/2. + (P(4,1)*Metric(1,2)*Metric(3,5))/2. + (P(4,2)*Metric(1,2)*Metric(3,5))/2. - P(4,3)*Metric(1,2)*Metric(3,5) - (P(2,1)*Metric(1,4)*Metric(3,5))/2. + (P(2,3)*Metric(1,4)*Metric(3,5))/2. - (P(1,2)*Metric(2,4)*Metric(3,5))/2. + (P(1,3)*Metric(2,4)*Metric(3,5))/2.')

VVVVV4 = Lorentz(name = 'VVVVV4',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,3)*Metric(1,4)*Metric(2,3) - P(4,3)*Metric(1,5)*Metric(2,3) - P(5,3)*Metric(1,3)*Metric(2,4) + P(4,3)*Metric(1,3)*Metric(2,5) + P(2,3)*Metric(1,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) - P(2,3)*Metric(1,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5)')

VVVVV5 = Lorentz(name = 'VVVVV5',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) + P(4,1)*Metric(1,3)*Metric(2,5) - P(3,1)*Metric(1,4)*Metric(2,5) - P(5,1)*Metric(1,2)*Metric(3,4) + P(2,1)*Metric(1,5)*Metric(3,4) + P(3,1)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5)')

VVVVV6 = Lorentz(name = 'VVVVV6',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) - P(5,1)*Metric(1,3)*Metric(2,4) + P(3,1)*Metric(1,5)*Metric(2,4) + P(4,1)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) - P(3,1)*Metric(1,2)*Metric(4,5) + P(2,1)*Metric(1,3)*Metric(4,5)')

VVVVV7 = Lorentz(name = 'VVVVV7',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,4)*Metric(1,3)*Metric(2,4) - P(3,4)*Metric(1,4)*Metric(2,5) - P(5,4)*Metric(1,2)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) + P(2,4)*Metric(1,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) + P(3,4)*Metric(1,2)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5)')

VVVVV8 = Lorentz(name = 'VVVVV8',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(3,5)*Metric(1,5)*Metric(2,4) - P(4,5)*Metric(1,3)*Metric(2,5) - P(2,5)*Metric(1,5)*Metric(3,4) + P(1,5)*Metric(2,5)*Metric(3,4) + P(4,5)*Metric(1,2)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5) - P(3,5)*Metric(1,2)*Metric(4,5) + P(2,5)*Metric(1,3)*Metric(4,5)')

VVVVV9 = Lorentz(name = 'VVVVV9',
                 spins = [ 3, 3, 3, 3, 3 ],
                 structure = 'P(5,1)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + 2*P(3,4)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) + P(4,1)*Metric(1,3)*Metric(2,5) - P(4,5)*Metric(1,3)*Metric(2,5) - P(3,1)*Metric(1,4)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + 2*P(3,5)*Metric(1,4)*Metric(2,5) - P(5,1)*Metric(1,2)*Metric(3,4) + P(5,4)*Metric(1,2)*Metric(3,4) + P(2,1)*Metric(1,5)*Metric(3,4) - 2*P(2,4)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) - P(4,1)*Metric(1,2)*Metric(3,5) + P(4,5)*Metric(1,2)*Metric(3,5) + P(2,1)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) - 2*P(2,5)*Metric(1,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) + 2*P(3,1)*Metric(1,2)*Metric(4,5) - P(3,4)*Metric(1,2)*Metric(4,5) - P(3,5)*Metric(1,2)*Metric(4,5) - 2*P(2,1)*Metric(1,3)*Metric(4,5) + P(2,4)*Metric(1,3)*Metric(4,5) + P(2,5)*Metric(1,3)*Metric(4,5)')

VVVVV10 = Lorentz(name = 'VVVVV10',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(4,2)*Metric(1,5)*Metric(2,3) + P(5,2)*Metric(1,3)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) - P(5,2)*Metric(1,2)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) + P(3,2)*Metric(1,2)*Metric(4,5) - P(1,2)*Metric(2,3)*Metric(4,5)')

VVVVV11 = Lorentz(name = 'VVVVV11',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,2)*Metric(1,4)*Metric(2,3) - P(5,2)*Metric(1,3)*Metric(2,4) + P(4,2)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) - P(4,2)*Metric(1,2)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(3,2)*Metric(1,2)*Metric(4,5) - P(1,2)*Metric(2,3)*Metric(4,5)')

VVVVV12 = Lorentz(name = 'VVVVV12',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,3)*Metric(1,4)*Metric(2,3) - P(4,3)*Metric(1,3)*Metric(2,5) - P(5,3)*Metric(1,2)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(4,3)*Metric(1,2)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVV13 = Lorentz(name = 'VVVVV13',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(4,3)*Metric(1,5)*Metric(2,3) - P(5,3)*Metric(1,3)*Metric(2,4) + P(5,3)*Metric(1,2)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(4,3)*Metric(1,2)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) + P(2,3)*Metric(1,3)*Metric(4,5) - P(1,3)*Metric(2,3)*Metric(4,5)')

VVVVV14 = Lorentz(name = 'VVVVV14',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,4)*Metric(1,4)*Metric(2,3) - P(3,4)*Metric(1,5)*Metric(2,4) - P(5,4)*Metric(1,2)*Metric(3,4) + P(2,4)*Metric(1,5)*Metric(3,4) - P(2,4)*Metric(1,4)*Metric(3,5) + P(1,4)*Metric(2,4)*Metric(3,5) + P(3,4)*Metric(1,2)*Metric(4,5) - P(1,4)*Metric(2,3)*Metric(4,5)')

VVVVV15 = Lorentz(name = 'VVVVV15',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,4)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,3)*Metric(2,4) + P(3,4)*Metric(1,5)*Metric(2,4) - P(3,4)*Metric(1,4)*Metric(2,5) - P(2,4)*Metric(1,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) + P(2,4)*Metric(1,3)*Metric(4,5) - P(1,4)*Metric(2,3)*Metric(4,5)')

VVVVV16 = Lorentz(name = 'VVVVV16',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,1)*Metric(1,4)*Metric(2,3) + P(5,2)*Metric(1,4)*Metric(2,3) - P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) + P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) - P(5,3)*Metric(1,3)*Metric(2,4) - P(5,4)*Metric(1,3)*Metric(2,4) - P(3,1)*Metric(1,5)*Metric(2,4) + P(3,4)*Metric(1,5)*Metric(2,4) - P(4,2)*Metric(1,3)*Metric(2,5) + P(4,3)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(3,4)*Metric(1,4)*Metric(2,5) - 2*P(5,1)*Metric(1,2)*Metric(3,4) - 2*P(5,2)*Metric(1,2)*Metric(3,4) + 2*P(5,3)*Metric(1,2)*Metric(3,4) + 2*P(5,4)*Metric(1,2)*Metric(3,4) + 2*P(2,1)*Metric(1,5)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,4)*Metric(1,5)*Metric(3,4) + 2*P(1,2)*Metric(2,5)*Metric(3,4) - P(1,3)*Metric(2,5)*Metric(3,4) - P(1,4)*Metric(2,5)*Metric(3,4) + P(4,1)*Metric(1,2)*Metric(3,5) + P(4,2)*Metric(1,2)*Metric(3,5) - 2*P(4,3)*Metric(1,2)*Metric(3,5) - P(2,1)*Metric(1,4)*Metric(3,5) + P(2,3)*Metric(1,4)*Metric(3,5) - P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) + P(3,1)*Metric(1,2)*Metric(4,5) + P(3,2)*Metric(1,2)*Metric(4,5) - 2*P(3,4)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,4)*Metric(1,3)*Metric(4,5) - P(1,2)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5)')

VVVVV17 = Lorentz(name = 'VVVVV17',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(4,5)*Metric(1,5)*Metric(2,3) - P(3,5)*Metric(1,4)*Metric(2,5) - P(2,5)*Metric(1,5)*Metric(3,4) + P(1,5)*Metric(2,5)*Metric(3,4) - P(4,5)*Metric(1,2)*Metric(3,5) + P(2,5)*Metric(1,4)*Metric(3,5) + P(3,5)*Metric(1,2)*Metric(4,5) - P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVV18 = Lorentz(name = 'VVVVV18',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(4,5)*Metric(1,5)*Metric(2,3) - P(3,5)*Metric(1,5)*Metric(2,4) - P(4,5)*Metric(1,3)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) - P(2,5)*Metric(1,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) + P(2,5)*Metric(1,3)*Metric(4,5) - P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVV19 = Lorentz(name = 'VVVVV19',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,1)*Metric(1,4)*Metric(2,3) - P(5,2)*Metric(1,4)*Metric(2,3) - P(4,1)*Metric(1,5)*Metric(2,3) + 2*P(4,2)*Metric(1,5)*Metric(2,3) - P(4,5)*Metric(1,5)*Metric(2,3) - P(5,1)*Metric(1,3)*Metric(2,4) + P(5,2)*Metric(1,3)*Metric(2,4) + P(3,1)*Metric(1,5)*Metric(2,4) - 2*P(3,2)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) + 2*P(4,1)*Metric(1,3)*Metric(2,5) - P(4,2)*Metric(1,3)*Metric(2,5) - P(4,5)*Metric(1,3)*Metric(2,5) - 2*P(3,1)*Metric(1,4)*Metric(2,5) + P(3,2)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) - P(4,1)*Metric(1,2)*Metric(3,5) - P(4,2)*Metric(1,2)*Metric(3,5) + 2*P(4,5)*Metric(1,2)*Metric(3,5) + P(2,1)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5) + P(3,1)*Metric(1,2)*Metric(4,5) + P(3,2)*Metric(1,2)*Metric(4,5) - 2*P(3,5)*Metric(1,2)*Metric(4,5) - P(2,1)*Metric(1,3)*Metric(4,5) + P(2,5)*Metric(1,3)*Metric(4,5) - P(1,2)*Metric(2,3)*Metric(4,5) + P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVV20 = Lorentz(name = 'VVVVV20',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,3)*Metric(1,4)*Metric(2,3) - P(5,4)*Metric(1,4)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - P(4,5)*Metric(1,5)*Metric(2,3) - P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) + P(3,4)*Metric(1,5)*Metric(2,4) - P(3,5)*Metric(1,5)*Metric(2,4) - P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,4)*Metric(1,4)*Metric(2,5) + P(3,5)*Metric(1,4)*Metric(2,5) - P(2,3)*Metric(1,5)*Metric(3,4) - P(2,4)*Metric(1,5)*Metric(3,4) + 2*P(2,5)*Metric(1,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) + P(1,4)*Metric(2,5)*Metric(3,4) - 2*P(1,5)*Metric(2,5)*Metric(3,4) - P(2,3)*Metric(1,4)*Metric(3,5) + 2*P(2,4)*Metric(1,4)*Metric(3,5) - P(2,5)*Metric(1,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - 2*P(1,4)*Metric(2,4)*Metric(3,5) + P(1,5)*Metric(2,4)*Metric(3,5) + 2*P(2,3)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,3)*Metric(2,3)*Metric(4,5) + P(1,4)*Metric(2,3)*Metric(4,5) + P(1,5)*Metric(2,3)*Metric(4,5)')

VVVVV21 = Lorentz(name = 'VVVVV21',
                  spins = [ 3, 3, 3, 3, 3 ],
                  structure = 'P(5,2)*Metric(1,4)*Metric(2,3) + P(5,3)*Metric(1,4)*Metric(2,3) - 2*P(5,4)*Metric(1,4)*Metric(2,3) + P(4,2)*Metric(1,5)*Metric(2,3) + P(4,3)*Metric(1,5)*Metric(2,3) - 2*P(4,5)*Metric(1,5)*Metric(2,3) - P(5,3)*Metric(1,3)*Metric(2,4) + P(5,4)*Metric(1,3)*Metric(2,4) - P(3,2)*Metric(1,5)*Metric(2,4) + P(3,5)*Metric(1,5)*Metric(2,4) - P(4,3)*Metric(1,3)*Metric(2,5) + P(4,5)*Metric(1,3)*Metric(2,5) - P(3,2)*Metric(1,4)*Metric(2,5) + P(3,4)*Metric(1,4)*Metric(2,5) - P(5,2)*Metric(1,2)*Metric(3,4) + P(5,4)*Metric(1,2)*Metric(3,4) - P(2,3)*Metric(1,5)*Metric(3,4) + P(2,5)*Metric(1,5)*Metric(3,4) + P(1,2)*Metric(2,5)*Metric(3,4) + P(1,3)*Metric(2,5)*Metric(3,4) - P(1,4)*Metric(2,5)*Metric(3,4) - P(1,5)*Metric(2,5)*Metric(3,4) - P(4,2)*Metric(1,2)*Metric(3,5) + P(4,5)*Metric(1,2)*Metric(3,5) - P(2,3)*Metric(1,4)*Metric(3,5) + P(2,4)*Metric(1,4)*Metric(3,5) + P(1,2)*Metric(2,4)*Metric(3,5) + P(1,3)*Metric(2,4)*Metric(3,5) - P(1,4)*Metric(2,4)*Metric(3,5) - P(1,5)*Metric(2,4)*Metric(3,5) + 2*P(3,2)*Metric(1,2)*Metric(4,5) - P(3,4)*Metric(1,2)*Metric(4,5) - P(3,5)*Metric(1,2)*Metric(4,5) + 2*P(2,3)*Metric(1,3)*Metric(4,5) - P(2,4)*Metric(1,3)*Metric(4,5) - P(2,5)*Metric(1,3)*Metric(4,5) - 2*P(1,2)*Metric(2,3)*Metric(4,5) - 2*P(1,3)*Metric(2,3)*Metric(4,5) + 2*P(1,4)*Metric(2,3)*Metric(4,5) + 2*P(1,5)*Metric(2,3)*Metric(4,5)')

SSSSSS1 = Lorentz(name = 'SSSSSS1',
                  spins = [ 1, 1, 1, 1, 1, 1 ],
                  structure = '1')

VVSSSS1 = Lorentz(name = 'VVSSSS1',
                  spins = [ 3, 3, 1, 1, 1, 1 ],
                  structure = 'Metric(1,2)')

VVVVSS1 = Lorentz(name = 'VVVVSS1',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)')

VVVVSS2 = Lorentz(name = 'VVVVSS2',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)')

VVVVSS3 = Lorentz(name = 'VVVVSS3',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)')

VVVVSS4 = Lorentz(name = 'VVVVSS4',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)')

VVVVSS5 = Lorentz(name = 'VVVVSS5',
                  spins = [ 3, 3, 3, 3, 1, 1 ],
                  structure = 'Metric(1,4)*Metric(2,3) - (Metric(1,3)*Metric(2,4))/2. - (Metric(1,2)*Metric(3,4))/2.')

VVVVVV1 = Lorentz(name = 'VVVVVV1',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) - Metric(1,4)*Metric(2,6)*Metric(3,5) - Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,3)*Metric(2,5)*Metric(4,6)')

VVVVVV2 = Lorentz(name = 'VVVVVV2',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,4)*Metric(2,6)*Metric(3,5) - Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVV3 = Lorentz(name = 'VVVVVV3',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,3)*Metric(2,6)*Metric(4,5) - Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,3)*Metric(2,5)*Metric(4,6) + Metric(1,2)*Metric(3,5)*Metric(4,6)')

VVVVVV4 = Lorentz(name = 'VVVVVV4',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,5)*Metric(2,4)*Metric(3,6) - Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,4)*Metric(2,3)*Metric(5,6) - Metric(1,3)*Metric(2,4)*Metric(5,6)')

VVVVVV5 = Lorentz(name = 'VVVVVV5',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,4)*Metric(3,5) - Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,3)*Metric(2,5)*Metric(4,6) - Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) - Metric(1,3)*Metric(2,4)*Metric(5,6)')

VVVVVV6 = Lorentz(name = 'VVVVVV6',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,3)*Metric(2,5)*Metric(4,6) - Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6)')

VVVVVV7 = Lorentz(name = 'VVVVVV7',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,4)*Metric(2,6)*Metric(3,5) - Metric(1,5)*Metric(2,4)*Metric(3,6) - Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,2)*Metric(3,5)*Metric(4,6) - Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6)')

VVVVVV8 = Lorentz(name = 'VVVVVV8',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) + Metric(1,5)*Metric(2,6)*Metric(3,4) - (Metric(1,6)*Metric(2,4)*Metric(3,5))/2. - (Metric(1,4)*Metric(2,6)*Metric(3,5))/2. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/2. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/2. - (Metric(1,6)*Metric(2,3)*Metric(4,5))/2. - (Metric(1,3)*Metric(2,6)*Metric(4,5))/2. + Metric(1,2)*Metric(3,6)*Metric(4,5) - (Metric(1,5)*Metric(2,3)*Metric(4,6))/2. - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. + Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - 2*Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV9 = Lorentz(name = 'VVVVVV9',
                  spins = [ 3, 3, 3, 3, 3, 3 ],
                  structure = 'Metric(1,5)*Metric(2,6)*Metric(3,4) + Metric(1,6)*Metric(2,4)*Metric(3,5) - Metric(1,4)*Metric(2,6)*Metric(3,5) - Metric(1,5)*Metric(2,4)*Metric(3,6) - Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,4)*Metric(2,3)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV10 = Lorentz(name = 'VVVVVV10',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) - Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,3)*Metric(2,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV11 = Lorentz(name = 'VVVVVV11',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,5)*Metric(2,4)*Metric(3,6) - Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV12 = Lorentz(name = 'VVVVVV12',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,4)*Metric(2,6)*Metric(3,5) - Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) - Metric(1,3)*Metric(2,5)*Metric(4,6) + Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,4)*Metric(2,3)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV13 = Lorentz(name = 'VVVVVV13',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,4)*Metric(2,6)*Metric(3,5) - Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,3)*Metric(2,4)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV14 = Lorentz(name = 'VVVVVV14',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - Metric(1,5)*Metric(2,4)*Metric(3,6) - Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) - Metric(1,3)*Metric(2,5)*Metric(4,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV15 = Lorentz(name = 'VVVVVV15',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,6)*Metric(2,4)*Metric(3,5) + Metric(1,6)*Metric(2,3)*Metric(4,5) - Metric(1,3)*Metric(2,6)*Metric(4,5) - Metric(1,5)*Metric(2,3)*Metric(4,6) + Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV16 = Lorentz(name = 'VVVVVV16',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,5)*Metric(2,6)*Metric(3,4) - Metric(1,4)*Metric(2,6)*Metric(3,5) - Metric(1,5)*Metric(2,4)*Metric(3,6) + Metric(1,4)*Metric(2,5)*Metric(3,6) - Metric(1,3)*Metric(2,5)*Metric(4,6) + Metric(1,2)*Metric(3,5)*Metric(4,6) + Metric(1,3)*Metric(2,4)*Metric(5,6) - Metric(1,2)*Metric(3,4)*Metric(5,6)')

VVVVVV17 = Lorentz(name = 'VVVVVV17',
                   spins = [ 3, 3, 3, 3, 3, 3 ],
                   structure = 'Metric(1,6)*Metric(2,5)*Metric(3,4) - (Metric(1,5)*Metric(2,6)*Metric(3,4))/2. + Metric(1,6)*Metric(2,4)*Metric(3,5) - (Metric(1,4)*Metric(2,6)*Metric(3,5))/2. - (Metric(1,5)*Metric(2,4)*Metric(3,6))/2. - (Metric(1,4)*Metric(2,5)*Metric(3,6))/2. - 2*Metric(1,6)*Metric(2,3)*Metric(4,5) + Metric(1,3)*Metric(2,6)*Metric(4,5) + Metric(1,2)*Metric(3,6)*Metric(4,5) + Metric(1,5)*Metric(2,3)*Metric(4,6) - (Metric(1,3)*Metric(2,5)*Metric(4,6))/2. - (Metric(1,2)*Metric(3,5)*Metric(4,6))/2. + Metric(1,4)*Metric(2,3)*Metric(5,6) - (Metric(1,3)*Metric(2,4)*Metric(5,6))/2. - (Metric(1,2)*Metric(3,4)*Metric(5,6))/2.')

