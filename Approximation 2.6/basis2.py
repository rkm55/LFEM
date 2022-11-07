# Basis 2 file taken from class github

# import unittest
import math
import numpy
# import sympy
# import joblib

def affine_mapping_1D( domain, target_domain, x ):
    A = numpy.array( [ [ 1.0, domain[0] ], [ 1.0, domain[1] ] ] )
    b = numpy.array( [target_domain[0], target_domain[1] ] )
    c = numpy.linalg.solve( A, b )
    fx = c[0] + c[1] * x
    return fx

def changeOfBasis( b1, b2, x1 ):
    T = numpy.linalg.solve( b2, b1 )
    x2 = numpy.dot( T, x1 )
    return x2, T

def evalMonomialBasis1D( degree, basis_idx, domain, variate ):
    return variate ** basis_idx

def evalLagrangeBasis1D( degree, basis_idx, domain, variate ):
    variate = affine_mapping_1D( domain, [-1, 1], variate )
    nodes = numpy.linspace( -1.0, 1.0, degree + 1 )
    basis_val = 1
    for i in range( 0, degree + 1 ):
        if ( i != basis_idx ):
            basis_val *= ( variate - nodes[i] ) / ( nodes[basis_idx] - nodes[i] )
    return basis_val

def evalBernsteinBasis1D( degree, basis_idx, domain, variate ):
    variate = affine_mapping_1D( domain, [0, 1], variate )
    if ( variate < 0.0 ) or ( variate > +1.0 ):
        raise Exception( "NOT_IN_DOMAIN" )
    term_1 = math.comb( degree, basis_idx )
    term_2 = variate ** basis_idx
    term_3 = ( 1.0 - variate ) ** ( degree - basis_idx )
    basis_val = term_1 * term_2 * term_3
    return basis_val

def evalLegendreBasis1D( degree, basis_idx, domain, variate ):
    variate = affine_mapping_1D( domain, [-1, 1], variate )
    if ( basis_idx == 0 ):
        basis_val = 1.0
    elif ( basis_idx == 1 ):
        basis_val = variate
    else:
        i = basis_idx - 1
        term_1 = ( ( 2 * i ) + 1 ) * variate * evalLegendreBasis1D( degree = degree, basis_idx = i, domain = [-1, 1], variate = variate )
        term_2 = i * evalLegendreBasis1D( degree = degree, basis_idx = i - 1, domain = [-1, 1], variate = variate )
        basis_val = ( term_1 - term_2 ) / ( i + 1 )
    return basis_val 

