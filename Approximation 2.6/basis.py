# basis functions
# author: ryan melander

import numpy
import math
import sympy


def evalMonomialBasis1D(variate, degree):
    return variate**degree

def evalLegendreBasis1D(variate,degree,basis_idx):
    if basis_idx == 0:
        val = 1
    elif basis_idx == 1:
        val = variate
    else:
        i = basis_idx - 1
        val = ((i + 1)**(-1)) * ((2*i+1)*variate*evalLegendreBasis1D(variate,degree,i) - i*evalLegendreBasis1D(variate,degree,i-1))
    return val

def evalLegendreBasis_idx(variate, degree, basis_idx):
    i = basis_idx
    term1 = 1/((2**i)*math.factorial(i))
    x = sympy.Symbol('x')
    fun = (x**2 - 1)**i
    term2 = sympy.diff(fun, x, i)
    val = term1*term2.subs(x, variate)
    return val

def evalLagrangeBasis1D(variate, degree, basis_idx):
    step = 2/degree
    xj = numpy.arange(-1,2,step)
    val = 1 
    for j in range(0,degree+1):
        if j == basis_idx:
            val = val
        else:
            val = val*(variate - xj[j])/(xj[basis_idx] - xj[j])
    return val

def evalBernsteinBasis1D(variate, degree, basis_idx):
    p = degree
    i = basis_idx
    v = (variate + 1)/2
    t1 = math.comb(p,i)
    t2 = v**i
    t3 = (1-v)**(p-i)
    val = t1*t2*t3
    return val