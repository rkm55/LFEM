# HW 2.1 - 2.3
# author: ryan melander

import unittest
import numpy as np

def evaluateLagrangeBasis1D(variate,degree,basis_idx):
    step = 2/degree
    xj = np.arange(-1,2,step)
    val = 1 
    for j in range(0,degree+1):
        if j == basis_idx:
            val = val
        else:
            val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
    return val

# def evaluateLagrangeBasis1D(variate,degree,basis_idx):
#     step = 2/degree
#     xj = np.arange(-1,2,step)
#     if degree == 1:
#         xj  = [-1, 1]
#     elif degree == 2:
#         xj  = [-1, 0, 1]
#     val = 1 
#     for j in range(0,degree+1):
#         if j == basis_idx:
#             val = val
#         else:
#             val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
#     return val

# def evaluateLagrangeBasis1D(variate,degree,basis_idx):
#     if degree == 1:
#         xj = [-1,1]
#         if basis_idx == 0:
#             val = (variate - xj[1]) / (xj[basis_idx] - xj[1])
#         else:
#             val = (variate - xj[0]) / (xj[basis_idx] - xj[0])
#     elif degree == 2:
#         xj = [-1,0,1]
#         if basis_idx == 0:
#             val = (variate - xj[1]) / (xj[basis_idx] - xj[1])*(variate - xj[2])/(xj[basis_idx] - xj[2])
#         elif basis_idx == 1:
#             val = (variate - xj[0]) / (xj[basis_idx] - xj[0])*(variate - xj[2])/(xj[basis_idx] - xj[2])
#         else:
#             val = (variate - xj[0]) / (xj[basis_idx] - xj[0])*(variate - xj[1])/(xj[basis_idx] - xj[1])
#     return val


class Test_evaluateLagrangeBasis1D( unittest.TestCase ):
    def test_linearLagrange( self ):
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 1, basis_idx = 0 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 1, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 1, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 1, basis_idx = 1 ), second = 1.0, delta = 1e-12 )

    def test_quadraticLagrange( self ):
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 2, basis_idx = 0 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 2, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = -1, degree = 2, basis_idx = 2 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate =  0, degree = 2, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate =  0, degree = 2, basis_idx = 1 ), second = 1.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate =  0, degree = 2, basis_idx = 2 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 2, basis_idx = 0 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 2, basis_idx = 1 ), second = 0.0, delta = 1e-12 )
        self.assertAlmostEqual( first = evaluateLagrangeBasis1D( variate = +1, degree = 2, basis_idx = 2 ), second = 1.0, delta = 1e-12 )
        
unittest.main()