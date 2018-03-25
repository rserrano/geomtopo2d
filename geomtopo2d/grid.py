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

from polygons import obtain_polygons

def boundary_edges( grid ):
    # These points will always be.
    points = [[-0.5, -0.5], [-0.5, len(grid)+0.5], [len(grid)+0.5, len(grid)+0.5], [len(grid)+0.5, -0.5]]
    edges = []
    ## First obtain the surrounding edges.
    # Obtain the left edge.
    last = 0
    for x in xrange(len(grid)-1):
        if grid[x][0] != grid[x+1][0]:
            edges.append( [last, len(points)] )
            last = len(points)
            points.append( [-0.5, x+0.5] )
    edges.append( [ last, 1 ] )
    # Obtain the bottom edge.
    last = 1
    for x in xrange(len(grid[0])-1):
        if grid[-1][x] != grid[-1][x+1]:
            edges.append( [last, len(points)] )
            last = len(points)
            points.append( [x+0.5, len(points)+0.5 ] )
    edges.append( [ last, 2 ] )
    # Obtain the right edge.
    last = 2
    for x in xrange(len(grid)-1):
        i = len(grid)-x-1
        if grid[i][0] != grid[i-1][0]:
            edges.append( [last, len(points)] )
            last = len(points)
            points.append( [len(points[0])+0.5, i-0.5] )
    edges.append( [ last, 3 ] )
    # Obtain the top edge.
    last = 3
    for x in xrange(len(grid[0])-1):
        i = len(grid)-x-1
        if grid[0][i] != grid[0][i-1]:
            edges.append( [last, len(points)] )
            last = len(points)
            points.append( [i-0.5, len(points)+0.5] )
    edges.append( [ last, 0 ] )
    return (points, edges)

def polygons_from_grid( grid ):
    """
    Given a grid in the form:
    [[A, B, A, A],
     [A, B, B, A],
     [A, B, A, A],
     [A, C, C, C]]
    Returns the polygons that surround the different values. In this, the algorithm sees the grid as:
    ·········
    ·A·B·A·A·
    ·········
    ·A·B·B·A·
    ·········
    ·A·B·A·A·
    ·········
    ·A·C·C·C·
    ·········
    Where the left, top is [-0.5, -0.5], the right, bottom is [4.5, 4.5].
    It then returns the polygons that surround the points.
    """
    
    
    ## Then obtain the points where the areas change.
    ## Find the points where there are more than one possibility to join, and add a point in the middle.
    ## Join the edges.
    ## Pass the edges to obtain polygons and return.
