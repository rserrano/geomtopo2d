"""
Copyright 2018 Geomodelr, Inc. 
rserrano at geomodelr.com

This file is part of Geomtopo2d. Geomtopo2d is free software: 
you can redistribute it and/or modify it under the terms of 
the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or
(at your option) any later version.

Geomtopo2d is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License
along with Geomtopo2d.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function, division

import unittest
import math
import numpy as np
from numpy import linalg as la
import numpy.random as nprnd
import cProfile, pstats, StringIO
import sys

import graphs
import polygons
import geometry
import grid 
from shapely.geometry import Polygon
from itertools import product

PROFILE = False

class Test(unittest.TestCase):
    
    def setUp(self):
        if PROFILE:
            # Profile GeoModelR
            self.pr = cProfile.Profile()
            self.pr.enable()
    
    def tearDown(self):
        if PROFILE:
            # Print profiling of GeoModelR.
            self.pr.disable()
            s = StringIO.StringIO()
            sortby = 'cumulative'
            ps = pstats.Stats(self.pr, stream=s).sort_stats(sortby)
            ps.print_stats()
            print >> sys.stderr, s.getvalue()
    
    def test_vector_angle(self):
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [ 0,-1] ), 0.0)
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [ 1,-1] ), 5.49778714378)
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [ 1, 0] ), 4.71238898038)
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [ 1, 1] ), 3.92699081699)
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [ 0, 1] ), 3.14159265359)
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [-1, 1] ), 2.35619449019)
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [-1, 0] ), 1.57079632679)
        self.assertAlmostEqual(geometry.vector_angle( [0, 1], [-1,-1] ), 0.785398163397)
    
    def test_graphs( self ):
        points = nprnd.uniform(0.0, 512.0, (1000,2))
        edges = graphs.relative_neighborhood_graph(points)
        holed, graph, points = polygons.obtain_polygons(edges, points)
        # polygon outupt: the polygons with holes, the areas of the polygons, the graph of neighbour polygons, 
        self.assertEqual(max(map( len, holed )), 1, "The relative neighborhood graph does not have holes")
        self.assertEqual(type(holed[0][0][0]), type(int()))
        self.assertEqual(type(holed[-1][-1][-1]), type(int()))
        self.assertEqual(len(graph), len(holed))
        for i in xrange(len(graph)):
            if len(graph[i]):
                self.assertIn(i, graph[graph[i][0]])
            else:
                print i, "g", graph[i], "p", holed[i]
    
    def test_grid( self ):
        ex = ["ABAA",
              "ABBA",
              "ABAA",
              "ACCC"]
        # Get the boundaries of the interpolation.
        edgese, pointse = grid.boundary_edges( ex )
        self.assertEqual( pointse, [(0.0, 0.0), (0.0, 3.0), (3.0, 3.0), (3.0, 0.0), (0.5, 3.0), (3.0, 2.5), (1.5, 0.0), (0.5, 0.0)] )
        self.assertEqual( edgese, [(0, 1), (1, 4), (4, 2), (2, 5), (5, 3), (3, 6), (6, 7), (7, 0)] )
        
        # Get all the points where polygons are different.
        pointsb = grid.point_breaks( ex )
        self.assertEqual(pointsb, [(0.5, 1.0), (2.0, 0.5), (2.5, 1.0), (2.0, 1.5), (0.5, 2.0), (1.5, 2.0), (1.0, 2.5), (2.0, 2.5)])
        
        # Get the mid points.
        edgesm, pointsm = grid.create_mid_points_and_edges( ex, pointse+pointsb )
        ot = set(pointsm)-set(pointse+pointsb)
        ot = sorted(list(ot))
        self.assertEqual(ot, [(0.5, 2.5), (1.5, 2.5)])
        self.assertEqual(len(set(pointse+pointsb)-set(pointsm)), 0)
        self.assertEqual(edgesm, [(8, 7), (9, 6), (10, 9), (12, 8), (13, 11), (10, 11), (4, 16), (14, 16), (12, 16), (15, 17), (13, 17), (14, 17), (5, 15)])
        
        polygons, cls, points = grid.polygons_from_grid( ex )
        self.assertEqual(cls, ['A', 'B', 'A', 'C'])
        self.assertEqual(polygons, [[[9, 6, 3, 5, 15, 17, 13, 11, 10]], [[10, 11, 13, 17, 14, 16, 12, 8, 7, 6, 9]], [[16, 4, 1, 0, 7, 8, 12]], [[14, 17, 15, 5, 2, 4, 16]]])
        ex = ["ABA",
              "ABB",
              "ABA",
              "ACC"]
        polygons, cls, points = grid.polygons_from_grid( ex )
        self.assertEqual(points, [(0.0, 0.0), (0.0, 3.0), (2.0, 3.0), (2.0, 0.0), (0.5, 3.0), (2.0, 2.5), (2.0, 1.5), (2.0, 0.5), 
                                  (1.5, 0.0), (0.5, 0.0), (0.5, 1.0), (0.5, 2.0), (1.5, 2.0), (1.0, 2.5), (0.5, 2.5), (1.5, 2.5)])
        self.assertEqual(polygons, [[[11, 10, 9, 8, 7, 6, 12, 15, 13, 14]], [[13, 15, 5, 2, 4, 14]], [[8, 3, 7]], [[14, 4, 1, 0, 9, 10, 11]], [[6, 5, 15, 12]]])
        self.assertEqual(cls, ['B', 'C', 'A', 'A', 'A'])
        ex = ["ABAA",
              "ABBA",
              "ABAA"]
        
        polygons, cls, points = grid.polygons_from_grid( ex )
        self.assertEqual(polygons, [[[6, 3, 2, 5, 11, 10, 9]], [[5, 4, 8, 7, 6, 9, 10, 11]], [[4, 1, 0, 7, 8]]])
        self.assertEqual(cls, ['A', 'B', 'A'])
        self.assertEqual(points, [(0.0, 0.0), (0.0, 2.0), (3.0, 2.0), (3.0, 0.0), (0.5, 2.0), (1.5, 2.0), (1.5, 0.0), (0.5, 0.0), (0.5, 1.0), (2.0, 0.5), (2.5, 1.0), (2.0, 1.5)])
        ex = ["AAAA",
              "ABCA",
              "AAAA",
              "CCAA",
              "CCCC"]
        polygons, cls, points = grid.polygons_from_grid( ex )
        self.assertEqual(polygons, [[[17, 9, 6, 7, 16, 8]], [[4, 0, 3, 5, 15, 14, 13], [9, 17, 12, 11, 10, 16, 7, 6]], [[16, 10, 11, 12, 17, 8]], [[4, 13, 14, 15, 5, 2, 1]]])
        self.assertEqual(points, [(0.0, 0.0), (0.0, 4.0), (3.0, 4.0), (3.0, 0.0), (0.0, 2.5), (3.0, 3.5), (0.5, 1.0), (1.0, 0.5), (1.5, 1.0), 
                                  (1.0, 1.5), (2.0, 0.5), (2.5, 1.0), (2.0, 1.5), (1.0, 2.5), (1.5, 3.0), (2.0, 3.5), (1.5, 0.5), (1.5, 1.5)])
        self.assertEqual(cls, ['B', 'A', 'C', 'C'])
        print grid.visi
def main(args=None):
    unittest.main()

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == "-p":
        del sys.argv[1]
        PROFILE = True
    main()
