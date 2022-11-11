import numpy
import matplotlib
import bext
import basis

def evaluateElementBernsteinBasisAtParamCoord( uspline, elem_id, param_coord ):
    elem_degree = bext.getElementDegree(uspline, elem_id)
    elem_bernstein_basis = numpy.zeros( elem_degree + 1 )
    elem_domain = bext.getElementDomain(uspline, elem_id)
    for n in range( 0, elem_degree + 1 ):
        elem_bernstein_basis[n] = basis.evalBernsteinBasis1D(param_coord, elem_degree, n, elem_domain)
    return elem_bernstein_basis

def evaluateElementSplineBasisAtParamCoord( uspline, elem_id, param_coord ):
    elem_ext_operator = bext.getElementExtractionOperator(uspline, elem_id)
    elem_bernstein_basis = evaluateElementBernsteinBasisAtParamCoord( uspline, elem_id, param_coord )
    elem_spline_basis = sum(numpy.transpose(elem_ext_operator*elem_bernstein_basis))
    return elem_spline_basis 

def plotUsplineBasis( uspline, color_by ):
    num_pts = 100
    xi = numpy.linspace( 0, 1, num_pts )
    num_elems = bext.getNumElems(uspline)
    for elem_idx in range( 0, num_elems ):
        elem_id = bext.elemIdFromElemIdx(uspline, elem_idx)
        elem_domain = bext.getElementDomain(uspline, elem_id)
        elem_node_ids = bext.getElementNodeIds(uspline, elem_id)
        elem_degree = bext.getElementDegree(uspline, elem_id)
        x = numpy.linspace( elem_domain[0], elem_domain[1], num_pts )
        y = numpy.zeros( shape = ( elem_degree + 1, num_pts ) )
        for i in range( 0, num_pts ):
            y[:,i] = evaluateElementSplineBasisAtParamCoord( uspline, elem_id, x )
        # Do plotting for the current element
        for n in range( 0, elem_degree + 1 ):
            if color_by == "element":
                color = getLineColor( elem_idx )
            elif color_by == "node":
                color = getLineColor( elem_node_ids[n] )
            matplotlib.pyplot.plot( x, y[n,:], color = getLineColor( color ) )
    # plt.show()

def getLineColor( idx ):
    colors = list( matplotlib.colors.TABLEAU_COLORS.keys() )
    num_colors = len( colors )
    mod_idx = idx % num_colors
    return matplotlib.colors.TABLEAU_COLORS[ colors[ mod_idx ] ]

uspline = bext.readBEXT( "quadratic_bspline_match.json" )
plotUsplineBasis( uspline, "element" )
plotUsplineBasis( uspline, "node" )