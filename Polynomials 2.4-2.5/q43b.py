# main HW 2.4 - 2.5
# author: ryan melander

import unittest
import numpy
import math


def computeRiemannQuadrature(fun, num_points):
    x, w = getRiemannQuadrature(num_points)
    y = numpy.zeros(len(x))
    for i in range (0, len(y)):
        y[i] = fun(x[i])
    s = numpy.sum(numpy.multiply(y,w))
    return s


# function from q43a
def getRiemannQuadrature(num_points):
    if num_points <= 0:
        raise ValueError('num_points_MUST_BE_INTEGER_GEQ_1')
    else:
        p = 1/num_points
        x = numpy.arange(p-1, 1, 2*p)
        w = numpy.full(len(x), 2*p)
        return x, w



class Test_computeRiemannQuadrature( unittest.TestCase ):
    def test_integrate_constant_one( self ):
        constant_one = lambda x : 1
        for num_points in range( 1, 100 ):
            self.assertAlmostEqual( first = computeRiemannQuadrature( fun = constant_one, num_points = num_points ), second = 2.0, delta = 1e-12 )

    def test_integrate_linear( self ):
        linear = lambda x : x
        for num_points in range( 1, 100 ):
            self.assertAlmostEqual( first = computeRiemannQuadrature( fun = linear, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_quadratic( self ):
        linear = lambda x : x**2
        error = []
        for num_points in range( 1, 100 ):
            error.append( abs( (2.0 / 3.0) - computeRiemannQuadrature( fun = linear, num_points = num_points ) ) )
        self.assertTrue( numpy.all( numpy.diff( error ) <= 0.0 ) )

    def test_integrate_sin( self ):
        sin = lambda x : math.sin(x)
        error = []
        for num_points in range( 1, 100 ):
            self.assertAlmostEqual( first = computeRiemannQuadrature( fun = sin, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_cos( self ):
        cos = lambda x : math.cos(x)
        error = []
        for num_points in range( 1, 100 ):
            error.append( abs( (2.0 / 3.0) - computeRiemannQuadrature( fun = cos, num_points = num_points ) ) )
        self.assertTrue( numpy.all( numpy.diff( error ) <= 0.0 ) )
        
unittest.main()