# quadrature functions
# author: ryan melander

import scipy
import numpy


def computeGaussLegendreQuadrature( n ):
    M = numpy.zeros( 2*n, dtype = "double" )
    M[0] = 2.0
    x0 = numpy.linspace( -1, 1, n )
    sol = scipy.optimize.least_squares( lambda x : objFun( M, x ), x0, bounds = (-1, 1), ftol = 1e-14, xtol = 1e-14, gtol = 1e-14 )
    qp = sol.x
    w = solveLinearMomentFit( M, qp )
    return qp, w

def assembleLinearMomentFitSystem( degree, pts ):
    A = numpy.zeros( shape = ( degree + 1, len( pts ) ), dtype = "double" )
    for m in range(0,degree + 1):
        for n in range(0,len(pts)):
            A[m,n] = evalLegendreBasis1D(m,pts[n])
    return A

def solveLinearMomentFit( M, pts ):
    degree = len( M ) - 1
    A = assembleLinearMomentFitSystem( degree, pts )
    sol = scipy.optimize.lsq_linear( A, M )
    w = sol.x
    return w

def objFun( M, pts ):
    degree = len( M ) - 1
    A = assembleLinearMomentFitSystem( degree, pts )
    w = solveLinearMomentFit( M, pts )
    obj_val = M - numpy.matmul(A,w)
    return obj_val

def evalLegendreBasis1D(degree,variate):
    if degree == 0:
        val = 1
    elif degree == 1:
        val = variate
    else:
        i = degree - 1
        val = ((i + 1)**(-1)) * ((2*i+1)*variate*evalLegendreBasis1D(i,variate) - i*evalLegendreBasis1D(i-1,variate))
    return val