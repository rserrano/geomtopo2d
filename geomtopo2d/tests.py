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

import matplotlib.pyplot as plt
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
            self.assertIn(i, graph[graph[i][0]])

#     def test_pseudo_circle_from_points( self ):
#         v = np.array([0.5, 0.3])
#         v /= la.norm(v)
#         pts = [ t*v for t in xrange(0, 20) ]
#         a, b, c = geometry.line_fit2(np.array(pts))
#         self.assertAlmostEqual(a[0],v[0])
#         self.assertAlmostEqual(a[1],v[1])
#         self.assertAlmostEqual(b[0], 0)
#         self.assertAlmostEqual(b[1], 0)
#         
#         pts.append([10, 10])
#         a, b, c = geometry.line_fit2(np.array(pts))
#         self.assertAlmostEqual(a[0],v[0])
#         self.assertAlmostEqual(a[1],v[1])
#         self.assertAlmostEqual(b[0], 0)
#         self.assertAlmostEqual(b[1], 0)
#         return
#         r = 3
#         points = []
#         for t in xrange( 10 ):
#             x = round(r*math.cos(2*math.pi*(t+5)/20))
#             y = round(r*math.sin(2*math.pi*(t+5)/20))
#             points.append([x,y])
#         for t in xrange( 1, 11 ):
#             points.append( [-3, t] )
#         
#         ln = np.array(points)[10:,:]
#         self.assertAlmostEqual(geometry.line_fit2(ln)[2],0.0)
#         ln = ln[:,[1,0]]
#         self.assertAlmostEqual(geometry.line_fit2(ln)[2],0.0)
#         
#         r = 3
#         points = []
#         for t in xrange( 10 ):
#             x = round(r*math.cos(2*math.pi*t/10))
#             y = round(r*math.sin(2*math.pi*t/10))
#             points.append([x,y])
#         
#         r, c, s= geometry.pseudo_circle_from_points2(np.array(points))
#         
#         self.assertAlmostEqual(r, 2.99628191397)
#         self.assertAlmostEqual(s, 0.047217255364194841)
#         
#         r = 3
#         points = []
#         for t in xrange( 5 ):
#             x = round(r*math.cos(2*math.pi*(t+5)/20))
#             y = round(r*math.sin(2*math.pi*(t+5)/20))
#             points.append([x,y])
#         r, c, s= geometry.pseudo_circle_from_points2(np.array(points))
#         
#         self.assertAlmostEqual(r, 2.9962819139658272)
#         self.assertAlmostEqual(s, 0.066775282914078068)
#         
#         r = 3
#         points = []
#         for t in xrange( 10 ):
#             x = round(r*math.cos(2*math.pi*(t+5)/20))
#             y = round(r*math.sin(2*math.pi*(t+5)/20))
#             points.append([x,y])
#         for t in xrange( 1, 11 ):
#             points.append( [t, -3] )
#         
#         points = np.array(points)
#         r, c, s= geometry.pseudo_circle_from_points2(points)
#         self.assertAlmostEqual(r, 4.8351351883071061)
#         self.assertAlmostEqual(s, 0.3184479828758664)
#         
#         r, c, s = geometry.pseudo_circle_from_points2(points[:10,:])
#         self.assertAlmostEqual(r, 2.9962819139658281)
#         self.assertAlmostEqual(s, 0.047217255364194841)
#         
#         w = geometry.what_is_fit ( points[10:,:] )
#         self.assertEqual(w, ("l", 0.0))
#         w = geometry.what_is_fit ( points[:10,:] )
#         self.assertEqual(w[0], "c")
#         
#         for t in xrange( 11, 31 ):
#             points = np.append(points, np.array([[t, -3]]), axis=0 )
#         
#         w = geometry.what_is_fit ( points[10:,:] )
#         self.assertEqual(w, ("l", 0.0))
#         w = geometry.what_is_fit( points )
#         self.assertAlmostEqual(w[1], 0.23176365560706205)
#         self.assertEqual(w[0], 'l')
# 
#         r = 10
#         points = []
#         for t in xrange( 20 ):
#             x = round(r*math.cos(2*math.pi*(t+10)/40))
#             y = round(r*math.sin(2*math.pi*(t+10)/40))
#             points.append([x,y])
#         for t in xrange( 1, 21 ):
#             points.append( [t, -10] )
#         
#         points = np.array(points)
#         w = geometry.what_is_fit( points[:20,:] )
#         self.assertEqual(w[0], 'c')
#         w = geometry.what_is_fit( points[20:,:] )
#         self.assertEqual(w[0], 'l')
#         w = geometry.what_is_fit( points )
#         self.assertEqual(w[0], 'c')
#         
#         r = 30
#         points = []
#         for t in xrange( 40 ):
#             x = round(r*math.cos(2*math.pi*(t+20)/80))
#             y = round(r*math.sin(2*math.pi*(t+20)/80))
#             points.append([x,y])
#         
#         for t in xrange( 1, 21 ):
#             points.append( [3*t, -30] )
#         
#         xl, yl = points[-1]
#         for t in xrange( 41 ):
#             x = round(xl+r*math.cos(2*math.pi*(t+60)/80))
#             y = round(r*math.sin(2*math.pi*(t+60)/80))
#             points.append([x,y])
#         
#         xl, yl = points[-1]
#         for t in xrange( 1, 21 ):
#             points.append( [xl-3*t, yl] )
#         
# 
#         line = [ i for i in xrange(len(points)) ]
#         line.append(0)
#         
#         r = 30
#         points = []
#         
#         for t in xrange( 30 ):
#             points.append( [t, round(2.0*t/3.0)] )
#         xl, yl = points[-1]
#         for t in xrange( 30 ):
#             points.append( [xl+t, yl-round(2.0*t/3.0)] )
#         xl, yl = points[-1]
#         points = list(reversed(points))
#         
#         r = 20
#         for t in xrange( 40 ):
#             x = round(r*math.cos(2*math.pi*(t+40)/80))
#             y = round(r*math.sin(2*math.pi*(t+40)/80))
#             points.append([x,y])
#         line = [ i for i in xrange(len(points)) ]
#         line.append(0)
        
def main(args=None):
    unittest.main()

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == "-p":
        del sys.argv[1]
        PROFILE = True
    main()
