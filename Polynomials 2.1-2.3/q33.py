# HW 2.1 - 2.3
# author: ryan melander

import numpy as np
#import math
#import sympy 
import matplotlib.pyplot as plt

x = np.arange(0,1.01,0.01)
y0 = np.ones(101)
y1 = x
y2 = np.square(x)
y3 = np.power(x,3)
y4 = np.power(x,4)
y5 = np.power(x,5)
y6 = np.power(x,6)
y7 = np.power(x,7)
y8 = np.power(x,8)
y9 = np.power(x,9)
y10 = np.power(x,10)

plt.figure(0)
plt.plot(x, y0, '-k', label = 'p = 0')
plt.plot(x, y1, '-b', label = 'p = 1')
plt.plot(x, y2, '-g', label = 'p = 2')
plt.plot(x, y3, '-r', label = 'p = 3')
plt.plot(x, y4, '-c', label = 'p = 4')
plt.plot(x, y5, '-m', label = 'p = 5')
plt.plot(x, y6, '-y', label = 'p = 6')
plt.plot(x, y7, '--k', label = 'p = 7')
plt.plot(x, y8, '--b', label = 'p = 8')
plt.plot(x, y9, '--g', label = 'p = 9')
plt.plot(x, y10, '--r', label = 'p = 10')
plt.legend(loc="upper left")
plt.xlim(0, 1)
plt.ylim(0, 2)
#plt.xlabel('P')
#plt.ylabel('| Error |')
plt.savefig('plot_q33.png', dpi=400)
