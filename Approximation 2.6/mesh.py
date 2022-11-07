# mesh function(s)
# author: ryan melander

import numpy

def generateMesh(xmin,xmax,num_elems,degree):
    if degree <= 0:
        raise ValueError('degree_MUST_BE_>=_1')
    ien_array = []
    node_coords = numpy.linspace(xmin,xmax,(degree*num_elems)+1)
    for i in range(0,len(node_coords)-1,degree):
        ien = numpy.linspace(i,i+degree,degree+1)
        ien = ien.astype(int)
        ien = list(ien)
        ien_array.append(ien)
    ien_array = numpy.asarray(ien_array)
    return node_coords, ien_array


def generateMeshElemWise(xmin,xmax,degree):
    ien_array = []
    key = []
    last_node = 0
    num_elems = len(degree)
    node_coords = []
    last_coord = xmin
    domain = xmax - xmin
    for i in range(0, num_elems):
        ien = numpy.linspace(last_node, degree[i] + last_node, degree[i] + 1)
        ien = list(ien.astype(int))
        last_node = ien[-1]
        key.append(i)
        ien_array.append(ien)
        elem_coords = list(numpy.linspace(last_coord, (domain/num_elems) + last_coord, degree[i] + 1))
        last_coord = elem_coords[-1]
        node_coords.extend(elem_coords)
    ien_array = dict(zip(key, ien_array))
    node_coords = numpy.unique(numpy.asarray(node_coords))
    return node_coords, ien_array


def computeSolution(target_fun, domain, num_elems, degree):
    node_coords, ien_array = generateMesh(domain[0], domain[1], num_elems, degree)
    test_solution = []
    for i in range(0,len(node_coords)):
        y = target_fun(node_coords[i])
        test_solution.append(y)
    test_solution = numpy.asarray(test_solution)
    return test_solution, node_coords, ien_array