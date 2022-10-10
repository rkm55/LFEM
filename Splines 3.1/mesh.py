# mesh function(s)
# author: ryan melander

import numpy

def generateMesh(xmin, xmax, num_elems, degree):
    if degree <= 0 or degree >= 3:
        raise ValueError('degree_MUST_BE_1=LINEAR_OR_2=QUADRATIC')
    elif degree == 1:
        node_coords = numpy.linspace(xmin, xmax, num_elems + 1)
        n = len(node_coords) - 1
        ien_array = []
        for i in range(0, n):
            ien = numpy.array([i, i + 1])
            ien_array.append(ien)
        ien_array = numpy.asarray(ien_array)
    elif degree == 2:
        node_coords = numpy.linspace(xmin, xmax, 2*num_elems + 1)
        n = len(node_coords) - 1
        ien_array = []
        for i in range(0, n, 2):
            ien = numpy.array([i, i + 1, i + 2])
            ien_array.append(ien)
        ien_array = numpy.asarray(ien_array)
    return node_coords, ien_array