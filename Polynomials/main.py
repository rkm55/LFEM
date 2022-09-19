# main HW 2.1 - 2.3
# author: ryan melander

import numpy as np
import sympy 
import matplotlib.pyplot as plt

def taylorExpansion( fun, a, order ):
    x = list( fun.atoms( sympy.Symbol ) )[0]
    t = 0
    for i in range( 0, order + 1 ):
       df = sympy.diff( fun, x, i )
       term = ( df.subs( x, a ) / sympy.factorial( i ) ) * ( x - a )**i
       t += term
    return t
   
xx = np.arange(-1,1.1,0.01)  
xxx = sympy.Symbol('xxx')  
plt.plot(xx,np.sin(np.pi*xx), "-k", label = "sin(x)")
Order0 = taylorExpansion(sympy.sin(np.pi*xxx),0,0)
Orders0 = np.zeros(len(xx))
for j in range(0,len(xx)):
    Orders0[j] = float(Order0.subs(xxx,xx[j]))
plt.plot(xx, Orders0, "m", label = "Order 0")
Orders = np.zeros((len(xx),4))
colors = ["b","r","g","c"]
for i in range(1,5):
    k = 2*i - 1
    Order = taylorExpansion(sympy.sin(np.pi*xxx),0,k)
    for j in range(0,len(xx)):
        Orders[j,i-1] = float(Order.subs(xxx,xx[j]))
    plt.plot(xx, Orders[:,i-1], colors[i-1], label = f"Order {str(k)}")
plt.legend(loc="upper left")
plt.xlim(-1, 1)
plt.ylim(-3, 3)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.savefig('plot_q32.png', dpi=600)