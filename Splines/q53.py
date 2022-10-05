# main splines
# author: ryan melander

# Write a function that computes a piecewise polynomial approximation, 
# using an interpolatory basis, of a univariate functions

import unittest
import numpy
import math
from q52 import generateMesh1D

def computeSolution(target_fun, domain, num_elems, degree):
    node_coords, ien_array = generateMesh1D(domain[0], domain[1], num_elems, degree)
    
    return 0

# def evaluateBernsteinBasis1D( variate, degree, basis_idx):
#     p = degree
#     i = basis_idx
#     v = (variate + 1)/2
#     t1 = math.comb(p,i)
#     t2 = v**i
#     t3 = (1-v)**(p-i)
#     val = t1*t2*t3
#     return val

class Test_computeSolution( unittest.TestCase ):
    def test_single_linear_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x, domain = [-1.0, 1.0 ], num_elems = 1, degree = 1 )
        gold_solution = numpy.array( [ -1.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
    
    def test_single_quad_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 1, degree = 2 )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
    
    def test_two_linear_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 2, degree = 1 )
        gold_solution = numpy.array( [ 1.0, 0.0, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
    
    def test_four_quad_element_poly( self ):
        test_solution, node_coords, ien_array = computeSolution( target_fun = lambda x : x**2, domain = [-1.0, 1.0 ], num_elems = 4, degree = 1 )
        gold_solution = numpy.array( [ 1.0, 0.25, 0.0, 0.25, 1.0 ] )
        self.assertTrue( numpy.allclose( test_solution, gold_solution ) )
        
unittest.main()