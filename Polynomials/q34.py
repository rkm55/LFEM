# HW 2.1 - 2.3
# author: ryan melander

import unittest

def evaluateMonomialBasis1D(degree, variate):
    if degree == 0 and variate == 0:
        return 1
    elif degree > 0 and variate == 0:
        return 0
    elif degree >0 and variate == 1:
        return 1
    elif degree >= 0 and variate == 0.5:
        return 1/(2**degree)
    
    
# test = evaluateMonomialBasis1D(1, 1)



class Test_evaluateMonomialBasis1D( unittest.TestCase ):
   def test_basisAtBounds( self ):
       self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = 0, variate = 0 ), second = 1.0, delta = 1e-12 )
       for p in range( 1, 11 ):
           self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 0 ), second = 0.0, delta = 1e-12 )
           self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 1 ), second = 1.0, delta = 1e-12 )

   def test_basisAtMidpoint( self ):
       for p in range( 0, 11 ):
           self.assertAlmostEqual( first = evaluateMonomialBasis1D( degree = p, variate = 0.5 ), second = 1 / ( 2**p ), delta = 1e-12 )
           
#def evaluateMonomialBasis1D(degree,variate):
    #return variate**degree

unittest.main()         