# main HW 2.4 - 2.5
# author: ryan melander
# testing file

import scipy
# from scipy import optimize
import numpy
# import unittest
# from src import basis

def computeGaussLegendreQuadrature( n ):
    M = numpy.zeros( 2*n, dtype = "double" ) # changed to 2*n-1
    M[0] = 2.0
    x0 = numpy.linspace( -1, 1, n )
    sol = scipy.optimize.least_squares( lambda x : objFun( M, x ), x0, bounds = (-1, 1), ftol = 1e-14, xtol = 1e-14, gtol = 1e-14 )
    qp = sol.x
    w = solveLinearMomentFit( M, qp )
    return qp, w

def assembleLinearMomentFitSystem( degree, pts ):
    A = numpy.zeros( shape = ( degree + 1, len( pts ) ), dtype = "double" ) 
    # use evaluate legendre function code from last assignment
    # nested loop to fill matrix
    for m in range(0, degree + 1):
        for n in range(0, len(pts)):
            A[m,n] = evalLegendreBasis1D(m, pts[n])
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
    ## YOUR CODE GOES HERE
    # M - Aw = 0 solve for left side as obj_val
    return obj_val

def evalLegendreBasis1D( degree, variate):
    if degree == 0:
        val = 1.0
    elif degree == 1:
        val = variate
    else:
        i = degree - 1
        term_1 = i * evalLegendreBasis1D(degree = i-1, variate = variate)
        term_2 = (2*i + 1) * variate * evalLegendreBasis1D( degree = i, variate = variate)
        val = (term_2 - term_1)/(i + 1)
    return val


########################## Testing

n = 3
qp, w = computeGaussLegendreQuadrature(n)
print(n, ' qp = ', qp, ' w = ', w)























