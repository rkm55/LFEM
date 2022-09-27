# main HW 2.4 - 2.5
# author: ryan melander

import unittest
import numpy

# def getNewtonCotesQuadrature(num_points):
#     if num_points <= 0 or num_points >= 7:
#         raise ValueError('num_points_MUST_BE_INTEGER_IN_[1,6]')
#     else:
#         p = 1/num_points
#         x = numpy.arange(p-1, 1, 2*p)
#         w = numpy.full(len(x), 2*p)
#         return x, w

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



class Test_getNewtonCotesQuadrature( unittest.TestCase ):
    def test_incorrect_num_points( self ):
        with self.assertRaises( Exception ) as context:
            getNewtonCotesQuadrature( num_points = 0 )
        self.assertEqual( "num_points_MUST_BE_INTEGER_IN_[1,6]", str( context.exception ) )
        with self.assertRaises( Exception ) as context:
            getNewtonCotesQuadrature( num_points = 7 )
        self.assertEqual( "num_points_MUST_BE_INTEGER_IN_[1,6]", str( context.exception ) )

    def test_return_types( self ):
        for num_points in range( 1, 7 ):
            x, w = getNewtonCotesQuadrature( num_points = num_points )
            self.assertIsInstance( obj = x, cls = numpy.ndarray )
            self.assertIsInstance( obj = w, cls = numpy.ndarray )
            self.assertTrue( len( x ) == num_points )
            self.assertTrue( len( w ) == num_points )
            
unittest.main()