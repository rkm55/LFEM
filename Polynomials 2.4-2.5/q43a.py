# main HW 2.4 - 2.5
# author: ryan melander

import unittest
import numpy


def getRiemannQuadrature(num_points):
    if num_points <= 0:
        raise ValueError('num_points_MUST_BE_INTEGER_GEQ_1')
    else:
        p = 1/num_points
        x = numpy.arange(p-1, 1, 2*p)
        w = numpy.full(len(x), 2*p)
        return x, w



class Test_getRiemannQuadrature( unittest.TestCase ):
    def test_zero_points( self ):
        with self.assertRaises( Exception ) as context:
            getRiemannQuadrature( num_points = 0 )
        self.assertEqual( "num_points_MUST_BE_INTEGER_GEQ_1", str( context.exception ) )

    def test_one_point( self ):
        x, w = getRiemannQuadrature( num_points = 1 )
        self.assertAlmostEqual( first = x, second = 0.0 )
        self.assertAlmostEqual( first = w, second = 2.0 )
        self.assertIsInstance( obj = x, cls = numpy.ndarray )
        self.assertIsInstance( obj = w, cls = numpy.ndarray )

    def test_two_point( self ):
        x, w = getRiemannQuadrature( num_points = 2 )
        self.assertTrue( numpy.allclose( x, [ -0.50, 0.50 ] ) )
        self.assertTrue( numpy.allclose( w, [ 1.0, 1.0 ] ) )
        self.assertIsInstance( obj = x, cls = numpy.ndarray )
        self.assertIsInstance( obj = w, cls = numpy.ndarray )

    def test_three_point( self ):
        x, w = getRiemannQuadrature( num_points = 3 )
        self.assertTrue( numpy.allclose( x, [ -2.0/3.0, 0.0, 2.0/3.0 ] ) )
        self.assertTrue( numpy.allclose( w, [ 2.0/3.0, 2.0/3.0, 2.0/3.0 ] ) )
        self.assertIsInstance( obj = x, cls = numpy.ndarray )
        self.assertIsInstance( obj = w, cls = numpy.ndarray )

    def test_many_points( self ):
        for num_points in range( 1, 100 ):
            x, w = getRiemannQuadrature( num_points = num_points )
            self.assertTrue( len( x ) == num_points )
            self.assertTrue( len( w ) == num_points )
            self.assertIsInstance( obj = x, cls = numpy.ndarray )
            self.assertIsInstance( obj = w, cls = numpy.ndarray )
            
unittest.main()