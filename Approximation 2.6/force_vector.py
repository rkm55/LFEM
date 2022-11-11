# approximation 2.6 force vector
# author: ryan melander

import numpy
import unittest
import basis 
import quadrature


def assembleForceVector(target_fun,domain,degree,solution_basis):
    nodes = degree + 1
    num_basis_vec = degree + 1
    F = numpy.zeros(num_basis_vec)
    qp, w = quadrature.computeGaussLegendreQuadrature(2*nodes)
    qp_fun = ((domain[-1] - domain[0])/2)*qp + (domain[0] + domain[-1])/2
    derivative = (domain[-1] - domain[0]) / 2
    for A in range(0,num_basis_vec):
        for k in range(0,len(qp)):
            F[A] += solution_basis(qp[k],degree,A) * target_fun(qp_fun[k]) * w[k] * derivative
    return F

# test = assembleForceVector(target_fun = lambda x: 2*x + numpy.pi, domain = [0,0.5], degree = 3, solution_basis = basis.evalLagrangeBasis1D)

if __name__ == '__main__':
    class Test_assembleForceVector( unittest.TestCase ):
        def test_legendre_const_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: numpy.pi, domain = [0, 1], degree = 1, solution_basis = basis.evalLegendreBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi, 0.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_legendre_linear_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: 2*x + numpy.pi, domain = [0, 1], degree = 1, solution_basis = basis.evalLegendreBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi + 1.0, 1.0/3.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_legendre_quadratic_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: x**2.0, domain = [0, 1], degree = 1, solution_basis = basis.evalLegendreBasis1D )
            gold_force_vector = numpy.array( [ 1.0/3.0, 1.0/6.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
            test_force_vector = assembleForceVector( target_fun = lambda x: x**2.0, domain = [0, 1], degree = 2, solution_basis = basis.evalLegendreBasis1D )
            gold_force_vector = numpy.array( [ 1.0/3.0, 1.0/6.0, 1.0/30.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
    
        def test_lagrange_const_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: numpy.pi, domain = [0, 1], degree = 1, solution_basis = basis.evalLagrangeBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi / 2.0, numpy.pi / 2.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_lagrange_linear_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: 2*x + numpy.pi, domain = [0, 1], degree = 1, solution_basis = basis.evalLagrangeBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi/2.0 + 1.0/3.0, numpy.pi/2.0 + 2.0/3.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_lagrange_quadratic_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: x**2.0, domain = [0, 1], degree = 1, solution_basis = basis.evalLagrangeBasis1D )
            gold_force_vector = numpy.array( [ 1.0/12.0, 1.0/4.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
            test_force_vector = assembleForceVector( target_fun = lambda x: x**2.0, domain = [0, 1], degree = 2, solution_basis = basis.evalLagrangeBasis1D )
            gold_force_vector = numpy.array( [ -1.0/60.0, 1.0/5.0, 3.0/20.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
    
        def test_bernstein_const_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: numpy.pi, domain = [0, 1], degree = 1, solution_basis = basis.evalBernsteinBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi / 2.0, numpy.pi / 2.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_bernstein_linear_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: 2*x + numpy.pi, domain = [0, 1], degree = 1, solution_basis = basis.evalBernsteinBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi/2.0 + 1.0/3.0, numpy.pi/2.0 + 2.0/3.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_bernstein_quadratic_force_fun( self ):
            test_force_vector = assembleForceVector( target_fun = lambda x: x**2.0, domain = [0, 1], degree = 1, solution_basis = basis.evalBernsteinBasis1D )
            gold_force_vector = numpy.array( [ 1.0/12.0, 1.0/4.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
            test_force_vector = assembleForceVector( target_fun = lambda x: x**2.0, domain = [0, 1], degree = 2, solution_basis = basis.evalBernsteinBasis1D )
            gold_force_vector = numpy.array( [ 1.0/30.0, 1.0/10.0, 1.0/5.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
            
    unittest.main()