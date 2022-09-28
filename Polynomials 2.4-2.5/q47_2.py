# main HW 2.4 - 2.5
# author: ryan melander

# second question 47 (written)
# calculates gauss legendre quadrature for monomials
# compare to definite integrals of the functions f

import math

for p in range(0,7):
    f = lambda x : x**p
    val = (8*f(0)/9) + (5/9)*(f(-1*math.sqrt(3/5))+f(math.sqrt(3/5)))
    print('degree = ', p, ' value = ', val)