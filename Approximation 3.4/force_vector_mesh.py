# approximation 3.4 force vector
# author: ryan melander

import numpy
import unittest
import basis 
import quadrature
import mesh


def assembleForceVector(target_fun, node_coords, ien_array, solution_basis):
    F = numpy.zeros(len(node_coords))
    for i in range(len(ien_array)):
        p = len(ien_array[i]) - 1
        node_idx = ien_array[i][0]
        domain = [node_coords[ien_array[i][0]], node_coords[ien_array[i][-1]]]
        qp, w = quadrature.computeGaussLegendreQuadrature(p+1)
        qp_fun = ((domain[-1] - domain[0])/2)*qp + (domain[0] + domain[-1])/2
        derivative = (domain[-1] - domain[0]) / 2
        for A in range(0, p + 1):
                for k in range(0,len(qp)):
                    F[node_idx] += solution_basis(qp[k], p, A, domain) * target_fun(qp_fun[k]) * w[k] * derivative
                node_idx += 1
    return F


# Testing
# target_fun = lambda x: x**2.0
# solution_basis = basis.evalLagrangeBasis1D
# domain = [ 0, 1 ]
# degree = [ 3, 3 ]
# node_coords, ien_array = mesh.generateMesh( domain[0], domain[1], degree )
# F = numpy.zeros(len(node_coords))
# for i in range(len(ien_array)):
#     p = len(ien_array[i]) - 1
#     node_idx = ien_array[i][0]
#     domain = [node_coords[ien_array[i][0]], node_coords[ien_array[i][-1]]]
#     qp, w = quadrature.computeGaussLegendreQuadrature(p+1)
#     qp_fun = ((domain[-1] - domain[0])/2)*qp + (domain[0] + domain[-1])/2
#     derivative = (domain[-1] - domain[0]) / 2
#     for A in range(0, p + 1):
#             for k in range(0,len(qp)):
#                 F[node_idx] += solution_basis(qp[k], p, A, domain) * target_fun(qp_fun[k]) * w[k] * derivative
#             node_idx += 1
    


# def assembleForceVector(target_fun,domain,degree,solution_basis):
#     nodes = degree + 1
#     num_basis_vec = degree + 1
#     F = numpy.zeros(num_basis_vec)
#     qp, w = quadrature.computeGaussLegendreQuadrature(nodes)
#     qp_fun = (qp+1)/2
#     derivative = (domain[-1] - domain[0]) / 2
#     for A in range(0,num_basis_vec):
#         for k in range(0,len(qp)):
#             F[A] += solution_basis(qp[k],degree,A) * target_fun(qp_fun[k]) * w[k] * derivative
#     return F



if __name__ == '__main__':
    class Test_assembleForceVector( unittest.TestCase ):
        def test_lagrange_const_force_fun( self ):
            domain = [ 0, 1 ]
            degree = [ 3, 3 ]
            target_fun = lambda x: numpy.pi
            node_coords, ien_array = mesh.generateMesh( domain[0], domain[1], degree )
            test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi / 16.0, 3.0 * numpy.pi / 16.0, 3.0 * numpy.pi / 16.0, numpy.pi / 8.0, 3.0 * numpy.pi / 16.0, 3.0 * numpy.pi / 16.0, numpy.pi / 16.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_lagrange_linear_force_fun( self ):
            domain = [ 0, 1 ]
            degree = [ 3, 3 ]
            target_fun = lambda x: 2*x + numpy.pi
            node_coords, ien_array = mesh.generateMesh( domain[0], domain[1], degree )
            test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
            gold_force_vector = numpy.array( [ 0.20468287, 0.62654862, 0.73904862, 0.51769908, 0.81404862, 0.92654862, 0.31301621 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_lagrange_quadratic_force_fun( self ):
            domain = [ 0, 1 ]
            degree = [ 3, 3 ]
            target_fun = lambda x: x**2.0
            node_coords, ien_array = mesh.generateMesh( domain[0], domain[1], degree )
            test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalLagrangeBasis1D )
            gold_force_vector = numpy.array( [ 1.04166667e-03, 0, 2.81250000e-02, 3.33333333e-02, 6.56250000e-02, 1.50000000e-01, 5.52083333e-02 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
    
        def test_bernstein_const_force_fun( self ):
            domain = [ 0, 1 ]
            degree = [ 3, 3 ]
            target_fun = lambda x: numpy.pi
            node_coords, ien_array = mesh.generateMesh( domain[0], domain[1], degree )
            test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
            gold_force_vector = numpy.array( [ numpy.pi / 8.0, numpy.pi / 8.0, numpy.pi / 8.0, numpy.pi / 4.0, numpy.pi / 8.0, numpy.pi / 8.0, numpy.pi / 8.0 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_bernstein_linear_force_fun( self ):
            domain = [ 0, 1 ]
            degree = [ 3, 3 ]
            target_fun = lambda x: 2*x + numpy.pi
            node_coords, ien_array = mesh.generateMesh( domain[0], domain[1], degree )
            test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
            gold_force_vector = numpy.array( [ 0.41769908, 0.44269908, 0.46769908, 1.03539816, 0.56769908, 0.59269908, 0.61769908 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
        
        def test_bernstein_quadratic_force_fun( self ):
            domain = [ 0, 1 ]
            degree = [ 3, 3 ]
            target_fun = lambda x: x**2.0
            node_coords, ien_array = mesh.generateMesh( domain[0], domain[1], degree )
            test_force_vector = assembleForceVector( target_fun = target_fun, node_coords = node_coords, ien_array = ien_array, solution_basis = basis.evalBernsteinBasis1D )
            gold_force_vector = numpy.array( [ 1/480, 1/160, 1/80, 1/15, 1/16, 13/160, 49/480 ] )
            self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
            
    unittest.main()