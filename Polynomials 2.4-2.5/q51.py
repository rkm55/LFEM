# main HW 2.4 - 2.5
# author: ryan melander

# import sympy
import math
import numpy as np
import matplotlib.pyplot as plt

# x = sympy.Symbol('x')
# fun1 = sympy.sin(x)
# fun2 = x - ((x**3)/6) - ((x**5)/120)
# fun3 = math.sqrt(sympy.integrate((fun1 - fun2)**2, (x, 0, 1)))

# analytical error
fun3 = lambda x : (x**11)/158400 + (x**9)/3240 + (x**7)/252 + (1/60)*(x**5)*math.cos(x) - (1/12)*(x**4)*math.sin(x) + x/2 - (1/4)*math.sin(2*x)
xx = np.arange(0,1.01,0.01)
er_a = np.zeros(len(xx), dtype = "double")
for i in range(0,len(xx)):
    er_a[i] = fun3(xx[i])


plt.figure(0)
plt.plot(xx,er_a, '-r', label = 'Analytical Error')
plt.legend(loc="upper left")
plt.xlim(0, 1)
plt.ylim(0, 0.25)
plt.xlabel('Range')
plt.ylabel('Error')
plt.savefig('plot_q51.png', dpi=350)