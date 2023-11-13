# Final Code Review

import unittest
import numpy
import uspline
import bext

def assembleForceVector (target_fun, uspline_bext):
    #stuff here
    return



class Test_assembleForceVector( unittest.TestCase ):
    def test_const_force_fun_two_element_linear_bspline( self ):
        target_fun = lambda x: numpy.pi
        spline_space = { "domain": [-1, 1], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ numpy.pi / 2.0, numpy.pi, numpy.pi / 2.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_linear_force_fun_two_element_linear_bspline( self ):
        target_fun = lambda x: x
        spline_space = { "domain": [-1, 1], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ -1.0/3.0, 0.0, 1.0/3.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_quadratic_force_fun_two_element_linear_bspline( self ):
        target_fun = lambda x: x**2
        spline_space = { "domain": [-1, 1], "degree": [ 1, 1 ], "continuity": [ -1, 0, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ 1.0/4.0, 1.0/6.0, 1.0/4.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_const_force_fun_two_element_quadratic_bspline( self ):
        target_fun = lambda x: numpy.pi
        spline_space = { "domain": [-1, 1], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ numpy.pi/3.0, 2.0*numpy.pi/3.0, 2.0*numpy.pi/3.0, numpy.pi/3.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_linear_force_fun_two_element_quadratic_bspline( self ):
        target_fun = lambda x: x
        spline_space = { "domain": [-1, 1], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ -1.0/4.0, -1.0/6.0, 1.0/6.0, 1.0/4.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )

    def test_quadratic_force_fun_two_element_quadratic_bspline( self ):
        target_fun = lambda x: x**2
        spline_space = { "domain": [-1, 1], "degree": [ 2, 2 ], "continuity": [ -1, 1, -1 ] }
        uspline.make_uspline_mesh( spline_space, "temp_uspline" )
        uspline_bext = bext.readBEXT( "temp_uspline.json" )
        test_force_vector = assembleForceVector( target_fun = target_fun, uspline_bext = uspline_bext )
        gold_force_vector = numpy.array( [ 2.0/10.0, 2.0/15.0, 2.0/15.0, 2.0/10.0 ] )
        self.assertTrue( numpy.allclose( test_force_vector, gold_force_vector ) )
