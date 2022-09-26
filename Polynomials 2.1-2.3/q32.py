# main HW 2.1 - 2.3
# author: ryan melander

import numpy as np
#import math
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


# taylor expansion
xxx = sympy.Symbol('xxx')  
fun1 = sympy.sin(np.pi*xxx)
fun2 = sympy.exp(xxx)
fun3 = sympy.erfc(xxx)
xx = np.arange(-2,2.01,0.01)
yy = np.zeros(len(xx))
for i in range (0,len(xx)):
    yy[i] = math.erfc(xx[i])
plt.figure(0)
plt.plot(xx,yy, "-k", label = "sin(x)")
Order0 = taylorExpansion(fun3,0,0)
Orders0 = np.zeros(len(xx))
for j in range(0,len(xx)):
    Orders0[j] = float(Order0.subs(xxx,xx[j]))
plt.figure(0)
plt.plot(xx, Orders0, "m", label = "Order 0")
Orders = np.zeros((len(xx),4))
colors = ["b","r","g","c"]
for i in range(1,5):
    k = 2*i - 1
    Order = taylorExpansion(fun3,0,k)
    for j in range(0,len(xx)):
        Orders[j,i-1] = float(Order.subs(xxx,xx[j]))
    plt.figure(0)
    plt.plot(xx, Orders[:,i-1], colors[i-1], label = f"Order {str(k)}")
plt.legend(loc="upper left")
plt.xlim(-2, 2)
plt.ylim(-1, 3)
plt.xlabel('x')
plt.ylabel('f(x)')
#plt.savefig('plot_q32_3.png', dpi=600)
#plt.figure(0) to name figures


# taylor expansion error
x = sympy.Symbol('x')
f1 = sympy.sin(np.pi*x)
f2 = sympy.exp(x)
f3 = sympy.erfc(x)


# f1
order = np.arange(0,11,1)
e1 = np.zeros(len(order))
for i in range(0,11):
    t = taylorExpansion(f1,0,order[i])
    e1[i] = float(sympy.integrate(abs(f1 - t),(x,-1,1)))
plt.figure(1)
plt.plot(order,e1, '-ok', label = r'$sin(\pi*x)$')
plt.legend(loc="upper right")
plt.yscale('log')
plt.grid(True)
plt.xlim(0, 10)
plt.ylim(0.001, 10)
plt.xlabel('Order')
plt.ylabel('| Error |')
plt.savefig('plot_q32_4.png', dpi=600)

# f2
order = np.arange(0,11,1)
e2 = np.zeros(len(order))
for i in range(0,11):
    t = taylorExpansion(f2,0,order[i])
    e2[i] = float(sympy.integrate(abs(f2 - t),(x,-1,1)))
plt.figure(2)
plt.plot(order,e2, '-ok', label = r'$e^{x}$')
plt.legend(loc="upper right")
plt.yscale('log')
plt.grid(True)
plt.xlim(0, 10)
plt.ylim(0.00001, 10)
plt.xlabel('Order')
plt.ylabel('| Error |')
plt.savefig('plot_q32_5.png', dpi=600)

# f3
order = np.arange(0,11,1)
e3 = np.zeros(len(order))
for i in range(0,11):
    t = taylorExpansion(f3,0,order[i])
    e3[i] = float(sympy.integrate(abs(f3 - t),(x,-2,2)))
plt.figure(3)
plt.plot(order,e3, '-ok', label = r'$erfc(x)$')
plt.legend(loc="upper right")
plt.yscale('log')
plt.grid(True)
plt.xlim(0, 10)
plt.ylim(0.000000001, 10)
plt.xlabel('Order')
plt.ylabel('| Error |')
plt.savefig('plot_q32_6.png', dpi=600)
















