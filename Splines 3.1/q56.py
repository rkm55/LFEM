# main splines 3.1
# author: ryan melander

# import scipy
import numpy
import math
from q54 import evaluateSolutionAt
import mesh
import basis
import matplotlib.pyplot as plt

# quadratic, sin(pi*x), [-1,1]
target_fun = lambda x : math.sin(math.pi*x)
degree = 2
domain = [-1,1]
eval_basis = basis.evalLagrangeBasis1D
num_elems = 3
x = numpy.linspace(-1,1,7)
solns = []
for i in range(0, len(x)):
    xx = x[i]
    coeff, node_coords, ien_array = mesh.computeSolution(target_fun, domain, num_elems, degree)
    sol_at_point = evaluateSolutionAt(xx, coeff, node_coords, ien_array, eval_basis)
    solns.append(sol_at_point)
    
plt.figure(0)
plt.plot(x, solns, '--ok', label = 'sin(pi*x)')
x1 = numpy.linspace(-1,1,100)
y1 = numpy.sin(numpy.pi*x1)
plt.plot(x1, y1, '-r')
for i in range(0, len(x)-1, 2):
    plt.axvline(x[i], color = 'k', linestyle = '--')
# matplotlib.pyplot.axvline(x, color, xmin, xmax, linestyle)
plt.plot()
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.legend(loc="upper left")
plt.xlabel('Number of Elements')
plt.ylabel('|Error|')
plt.savefig('fig_q56', dpi=325)