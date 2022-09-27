# main HW 2.4 - 2.5
# author: ryan melander

import unittest
import math
import numpy


def computeNewtonCotesQuadrature(fun, num_points):
    x, w = getNewtonCotesQuadrature(num_points)
    f = numpy.zeros(len(x))
    for i in range (0, len(f)):
        f[i] = fun(x[i])
    s = numpy.sum(numpy.multiply(f,w))
    return s


def getNewtonCotesQuadrature(num_points):
    if num_points <= 0 or num_points >= 7:
        raise ValueError('num_points_MUST_BE_INTEGER_IN_[1,6]')
    elif num_points == 1:
        x = numpy.array([0])
        w = numpy.array([2])
    else:
        a = -1
        b = 1
        n = num_points - 1
        h = (b-a)/n
        x = numpy.arange(-1, 1+h, h)
        if num_points == 2:
            w = [h/2, h/2]
        elif num_points == 3:
            w = [h/3, 4*h/3, h/3]
        elif num_points == 4:
            w = [3*h/8, 9*h/8, 9*h/8, 3*h/8]
        elif num_points == 5:
            w = [14*h/45, 64*h/45, 24*h/45, 64*h/45, 14*h/45]
        elif num_points == 6:
            w = [95*h/288, 375*h/288, 250*h/288, 250*h/288, 375*h/288, 95*h/288]
        w = numpy.array(w)
    return x, w



class Test_computeNewtonCotesQuadrature( unittest.TestCase ):
    def test_integrate_constant_one( self ):
        constant_one = lambda x : 1 * x**0
        for degree in range( 1, 6 ):
            num_points = degree + 1
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = constant_one, num_points = num_points ), second = 2.0, delta = 1e-12 )

    def test_exact_poly_int( self ):
        for degree in range( 1, 6 ):
            num_points = degree + 1
            poly_fun = lambda x : ( x + 1.0 ) ** degree
            indef_int = lambda x : ( ( x + 1 ) ** ( degree + 1) ) / ( degree + 1 )
            def_int = indef_int(1.0) - indef_int(-1.0)
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = poly_fun, num_points = num_points ), second = def_int, delta = 1e-12 )

    def test_integrate_sin( self ):
        sin = lambda x : math.sin(x)
        for num_points in range( 1, 7 ):
            self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = sin, num_points = num_points ), second = 0.0, delta = 1e-12 )

    def test_integrate_cos( self ):
        cos = lambda x : math.cos(x)
        self.assertAlmostEqual( first = computeNewtonCotesQuadrature( fun = cos, num_points = 6 ), second = 2*math.sin(1), delta = 1e-4 )
        
unittest.main()