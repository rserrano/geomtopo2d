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

from scipy.spatial import Delaunay, distance, cKDTree
from scipy.stats import norm as normal
import itertools
import math
import numpy as np
from numpy import linalg as la

def srtedg( e ):
    """
    Small routine to order an edge from its lower to its upper index.
    """
    if e[0] > e[1]:
        return (e[1], e[0])
    return (e[0], e[1])

def delaunay_graph( points ):
    """
    Calculates the delaunay triangulation of a set of points, and returns
    its edges.
    """
    # Calculate the Delaunay triangulation.
    triangulation = Delaunay(points)
    # Get all its edges
    edgs = set()
    for tri in triangulation.simplices:
        triedgs = []
        for edg in itertools.combinations(tri,2):
            edgs.add(srtedg(edg))
    
    return list(edgs)

def gabriel_graph( points, matrix=None, tree_c=None, dt_c=None ):
    """
    Calculates the gabriel graph of a set of points, and then returns
    its edges.
    The gabriel graph is the pair of points, that given a circle (sphere, hypersphere),
    with its center in the middle of the points that passes through those points,
    the circle does not contain any other point of the sample.
    """
    if matrix is None:
        matrix = {}
    initial = delaunay_graph( points )
    if dt_c is not None:
        dt_c['dt'] = initial
    tree = cKDTree( points )
    if tree_c is not None:
        tree_c['tree'] = tree
    gabriel = []
    for e0, e1 in initial:
        p0 = np.array(points[e0])
        p1 = np.array(points[e1])
        mid = (p0+p1)/2.0
        dis = la.norm(p1-p0)
        rad = dis/2.0
        cls = tree.query_ball_point( mid, rad )
        matrix[(e0, e1)] = dis
        for i in cls:
            if i != e0 and i != e1:
                break
        else:
            gabriel.append( ( e0, e1 ) )
    return gabriel

def relative_neighborhood_graph( points, matrix=None, dt_c=None ):
    """
    Calculates the relative neighborhood graph.
    The relative neighborhood graph is the graph that
    given two points, it has an edge if there's no point
    that's closer to both of them.
    """
    if matrix is None:
        matrix = {}
    tree_c = {}
    initial = gabriel_graph( points, matrix, tree_c, dt_c )
    tree = tree_c['tree']
    rng = []
    sin60 = math.sin(math.radians(60))
    for e0, e1 in initial:
        p0 = points[e0]
        p1 = points[e1]
        p0 = np.array(p0)
        p1 = np.array(p1)
        mid = (p0+p1) / 2.0
        dis = matrix[(e0, e1)]
        rad = dis*sin60
        cls = tree.query_ball_point( mid, rad )
        for i in cls:
            if i != e0 and i != e1:
                ei0 = srtedg( (i, e0) )
                if ei0 in matrix:
                    d0 = matrix[ei0]
                else:
                    d0 = la.norm( np.array(points[i])-p0 )
                    matrix[ei0] = d0
                ei1 = srtedg( (i, e1) )
                if ei1 in matrix:
                    d1 = matrix[ei1]
                else:
                    d1 = la.norm( np.array(points[i])-p1 )
                    matrix[ei1] = d1
                
                if d0 < dis and d1 < dis:
                    break
        else:
            rng.append( ( e0, e1 ) )
    return rng

