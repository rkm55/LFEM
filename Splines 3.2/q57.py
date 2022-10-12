# main splines 3.1
# author: ryan melander

import unittest
import numpy


# generate mesh with support for element wise specification of polynomial degree
# degree input as array where length = num_elems and each entry is the degree of
# that element

def generateMesh(xmin,xmax,degree):
    ien_array = []
    key = []
    last_node = 0
    num_elems = len(degree)
    node_coords = []
    last_coord = xmin
    domain = xmax - xmin
    for i in range(0, num_elems):
        # ien_array
        ien = numpy.linspace(last_node, degree[i] + last_node, degree[i] + 1)
        ien = list(ien.astype(int))
        last_node = ien[-1]
        key.append(i)
        ien_array.append(ien)
        # node_coords
        elem_coords = list(numpy.linspace(last_coord, (domain/num_elems) + last_coord, degree[i] + 1))
        last_coord = elem_coords[-1]
        node_coords.extend(elem_coords)
    ien_array = dict(zip(key, ien_array))
    node_coords = numpy.unique(numpy.asarray(node_coords))
    return node_coords, ien_array


# testing
node_coords, ien_array = generateMesh(0, 4, [1,2,3,4])



class Test_generateMesh( unittest.TestCase ):
    def test_make_1_linear_elem( self ):
        gold_node_coords = numpy.array( [ 0.0, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1 ] }
        node_coords, ien_array = generateMesh( xmin = 0.0, xmax = 1.0, degree = [ 1 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_1_quadratic_elem( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.5, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1, 2 ] }
        node_coords, ien_array = generateMesh( xmin = 0.0, xmax = 1.0, degree = [ 2 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_2_linear_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.5, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1 ], 1: [ 1, 2 ] }
        node_coords, ien_array = generateMesh( xmin = 0.0, xmax = 1.0, degree = [ 1, 1 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_2_quadratic_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.25, 0.5, 0.75, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1, 2 ], 1: [ 2, 3, 4 ] }
        node_coords, ien_array = generateMesh( xmin = 0.0, xmax = 1.0, degree = [ 2, 2 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_4_linear_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.25, 0.5, 0.75, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1 ], 1: [ 1, 2 ], 2: [ 2, 3 ], 3: [ 3, 4 ] }
        node_coords, ien_array = generateMesh( xmin = 0.0, xmax = 1.0, degree = [ 1, 1, 1, 1 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_4_quadratic_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0 ] )
        gold_ien_array = { 0: [ 0, 1, 2 ], 1: [ 2, 3, 4 ], 2: [ 4, 5, 6 ], 3: [ 6, 7, 8 ] }
        node_coords, ien_array = generateMesh( xmin = 0.0, xmax = 1.0, degree = [ 2, 2, 2, 2 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )
    
    def test_make_4_p_refine_elems( self ):
        gold_node_coords = numpy.array( [ 0.0, 1.0, 1.5, 2.0, (2.0 + 1.0/3.0), (2.0 + 2.0/3.0), 3.0, 3.25, 3.5, 3.75, 4.0 ] )
        gold_ien_array = { 0: [ 0, 1 ], 1: [ 1, 2, 3 ], 2: [ 3, 4, 5, 6 ], 3: [ 6, 7, 8, 9, 10 ] }
        node_coords, ien_array = generateMesh( xmin = 0.0, xmax = 4.0, degree = [ 1, 2, 3, 4 ] )
        self.assertIsInstance( obj = node_coords, cls = numpy.ndarray )
        self.assertIsInstance( obj = ien_array, cls = dict )
        self.assertTrue( numpy.allclose( node_coords, gold_node_coords ) )
        self.assertDictEqual( d1 = gold_ien_array, d2 = ien_array )

unittest.main()