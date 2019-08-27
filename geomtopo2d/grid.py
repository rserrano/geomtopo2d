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

from .polygons import obtain_polygons
import rtree
import math
import numpy as np

def boundary_edges( grid ):
    """
    Given a grid, obtains the surrounding edges of the polygon.
    """
    # These points will always be.
    points = [(0.0, 0.0), (0.0, len(grid)-1.0), (len(grid[0])-1.0, len(grid)-1.0), (len(grid[0])-1.0, 0.0)]
    edges = []
    # Obtain the left edge.
    last = 0
    for x in range(len(grid)-1):
        if grid[x][0] != grid[x+1][0]:
            edges.append( (last, len(points)) )
            last = len(points)
            points.append( (0.0, x+0.5) )
    edges.append( ( last, 1 ) )
    # Obtain the bottom edge.
    last = 1
    for x in range(len(grid[0])-1):
        if grid[-1][x] != grid[-1][x+1]:
            edges.append( (last, len(points)) )
            last = len(points)
            points.append( (x+0.5, len(grid)-1.0 ) )
    edges.append( ( last, 2 ) )
    # Obtain the right edge.
    last = 2
    for x in range(len(grid)-1):
        i = len(grid)-x-1
        if grid[i][-1] != grid[i-1][-1]:
            edges.append( (last, len(points)) )
            last = len(points)
            points.append( (len(grid[0])-1.0, i-0.5) )
    edges.append( ( last, 3 ) )
    # Obtain the top edge.
    last = 3
    for x in range(len(grid[0])-1):
        i = len(grid[0])-x-1
        if grid[0][i] != grid[0][i-1]:
            edges.append( (last, len(points)) )
            last = len(points)
            points.append( (i-0.5, 0.0) )
    edges.append( ( last, 0 ) )
    return (edges, points)

def point_breaks( grid ):
    points = []
    for y in range( 1, len(grid)-1 ):
        for x in range( 1, len(grid[0])-1 ):
            # Add breaks in the left line.
            if x == 1:
                if grid[y][x] != grid[y][x-1]:
                    points.append( (x-0.5, float(y)) )
            # Add breaks in the top line.
            if y == 1:
                if grid[y][x] != grid[y-1][x]:
                    points.append( (float(x), y-0.5) )
            # Add all the rest of breaks.
            if grid[y][x] != grid[y][x+1]:
                points.append( (x+0.5, float(y)) )
            if grid[y][x] != grid[y+1][x]:
                points.append( (float(x), y+0.5) )
    return points

def create_mid_points_and_edges( grid, points ):
    spoints = {}
    for i, k in enumerate(points):
        spoints[k] = i
    edges = []
    # Add all midpoints and append the given edges.
    for y in range( len(grid)-1 ):
        for x in range( len(grid[0])-1 ):
            midx = x+0.5
            midy = y+0.5
            cnt = 0
            es = []
            for d in [ [0, 0.5], [0.5, 0], [0, -0.5], [-0.5, 0] ]:
                ex = midx + d[0]
                ey = midy + d[1]
                if ( ex, ey ) in spoints:
                    cnt += 1
                    es.append( spoints[(ex, ey)] )
            if cnt == 0:
                continue
            elif cnt == 1:
                raise Exception("error tying polygons")
            elif cnt == 2:
                edges.append( (es[0], es[1]) )
            elif cnt > 2:
                l = len(spoints)
                spoints[ ( midx, midy ) ] = l
                for e in es:
                    edges.append( (e, l) )

    ret_points = [ -1 for i in range( len( spoints ) ) ]
    for k, v in spoints.items():
        ret_points[v] = k
    return edges, ret_points

def classify_polygon( grid, polygon, points ):
    for r in polygon:
        for idx, n in enumerate(r):
            nn = r[idx+1]
            pn = points[n]
            pnn = points[nn]
            # Grid or middle node.
            loc = (pn[0]-math.floor(pn[0]), pn[1]-math.floor(pn[1]))

            # If it's an grid node, then classify grid value.
            if loc[0] == 0.0 and loc[1] == 0.0:
                return grid[int(pn[1])][int(pn[0])]
            # Leaving middle nodes, classify at next.
            if loc[0] == 0.5 and loc[1] == 0.5:
                continue
            # Check if it's diagonal or not.
            v = ( pnn[0]-pn[0], pnn[1]-pn[1] )
            if v[0] == 0.0:
                if loc[0] == 0.0:
                    if v[1] < 0.0:
                        return grid[int(pn[1]-0.5)][int(pn[0])]
                    else:
                        return grid[int(pn[1]+0.5)][int(pn[0])]
                if v[1] < 0.0:
                    return grid[int(pn[1])][int(pn[0]+0.5)]
                elif v[1] > 0.0:
                    return grid[int(pn[1])][int(pn[0]-0.5)]
                else:
                    raise Exception("Error classifying polygons.")
            elif v[1] == 0.0:
                if loc[1] == 0.0:
                    if v[0] < 0.0:
                        return grid[int(pn[1])][int(pn[0]-0.5)]
                    else:
                        return grid[int(pn[1])][int(pn[0]+0.5)]
                if v[0] < 0.0:
                    return grid[int(pn[1]-0.5)][int(pn[0])]
                elif v[0] > 0.0:
                    return grid[int(pn[1]+0.5)][int(pn[0])]
                else:
                    raise Exception("Error classifying polygons.")
            elif v[0] < 0.0:
                if v[1] < 0.0:
                    return grid[int(pn[1]-0.5)][int(pn[0])]
                elif v[1] > 0.0:
                    return grid[int(pn[1])][int(pn[0]-0.5)]
                else:
                    raise Exception("Error classifying polygons.")
            elif v[0] > 0.0:
                if v[1] < 0.0:
                    return grid[int(pn[1])][int(pn[0]+0.5)]
                elif v[1] > 0.0:
                    return grid[int(pn[1]+0.5)][int(pn[0])]
                else:
                    raise Exception("Error classifying polygons.")


def polygons_from_grid( grid ):
    """
    Given a grid in the form:
    [[A, B, A, A],
     [A, B, B, A],
     [A, B, A, A],
     [A, C, C, C]]
    Returns the polygons that surround the different values. In this, the algorithm sees the grid as:
    *********
    *A*B*A*A*
    *********
    *A*B*B*A*
    *********
    *A*B*A*A*
    *********
    *A*C*C*C*
    *********
    Where the left, top is [-0.5, -0.5], the right, bottom is [3.5, 3.5].
    It then returns the polygons that surround the points.
    """
    # First obtain the surrounding edges.
    edges, points = boundary_edges( grid )
    
    # Then obtain the points where the areas change.
    points += point_breaks( grid )
    
    # Find the points where there are more than one possibility to join, and add a point in the middle. 
    # Also create the edges in the interior.
    edgesm, points = create_mid_points_and_edges( grid, points )
    edges += edgesm
    # Pass the edges to obtain polygons and return.
    polygons, graph, points = obtain_polygons( edges, points )
    classification = list(map( lambda p: classify_polygon(grid, p, points), polygons ))
    return polygons, classification, points

