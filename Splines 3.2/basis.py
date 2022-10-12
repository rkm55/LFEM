# basis functions
# author: ryan melander

import numpy
import math


def evalMonomialBasis1D(degree,variate):
    return variate**degree

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

def evalLagrangeBasis1D(variate,degree,basis_idx):
    step = 2/degree
    xj = numpy.arange(-1,2,step)
    val = 1 
    for j in range(0,degree+1):
        if j == basis_idx:
            val = val
        else:
            val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
    return val

def evalBernsteinBasis1D( variate, degree, basis_idx):
    p = degree
    i = basis_idx
    v = (variate + 1)/2
    t1 = math.comb(p,i)
    t2 = v**i
    t3 = (1-v)**(p-i)
    val = t1*t2*t3
    return val