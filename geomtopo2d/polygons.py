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

import math
import numpy as np
import rtree
from numpy import linalg as la
from shapely.geometry import Polygon
from .geometry import vector_angle

def separate_lines(edgs):
    """
    Given a general graph, as represented by edges,
    it splits the graph at the points where three or
    more edges converge, and returns
    the edges as lines.
    Parameters
    ----------
    edgs : dict
        A set of edges to join together.
    Results
    -------
    dict : the edges joined at the corners.
    """
    # Find the points with respective edges.
    pnh = {}
    for edg in edgs:
        for p in edg:
            if not p in pnh:
                pnh[p] = [edg]
            else:
                pnh[p].append(edg)
    
    rem_edges = set(edgs)
    cor_edges = []
    end_edges = []

    for edg in rem_edges:
        n0 = edg[0]
        n1 = edg[1]
        c0 = len(pnh[n0])
        c1 = len(pnh[n1])
        
        # Return immediately if it's the end of a line.
        if c0 == 1:
            end_edges.append( (n0, n1) )
            continue
        elif c1 == 1:
            end_edges.append( (n1, n0) )
            continue
        # Save to return if it's in a 3 point.
        if c0 > 2:
            cor_edges.append( (n0, n1) )
            continue
        if c1 > 2:
            cor_edges.append( (n1, n0) )
    
    def new_line_edge():
        while ( len( end_edges ) ):
            e = end_edges.pop()
            if e in rem_edges:
                return e
        
        while ( len( cor_edges ) ):
            e = cor_edges.pop()
            if e in rem_edges:
                return e
        
        for e in rem_edges:
            return e
        return None
    
    def remove_edge(e):
        """
        Removes an edge from the set of edges, even if it's in the oposite direction.
        """
        if e in rem_edges:
            rem_edges.remove(e)
        else:
            if (e[1], e[0]) in rem_edges:
                rem_edges.remove((e[1], e[0]))
            else:
                print( e )
                assert False, "An edge not in to remove" 
    
    equal_edge = lambda e1, e2: (e1[0] == e2[0] and e1[1] == e2[1]) or \
                                (e1[1] == e2[0] and e1[0] == e2[1])
    
    def next_edge(e):
        """
        Picks the next edge from this node or None if it's in 
        the end of a line or if the corner is divided in 3.
        """
        # Get possible edges from the list.
        pose = pnh[e[1]]
        
        # If the possible edges are more than 2 or less than two, return None.
        if len(pose) != 2:
            return None
        # Now, for all the possible e, 
        for pe in pose:
            if equal_edge(pe, e) or not pe in rem_edges:
                continue
            if pe[1] == e[1]:
                return (pe[1], pe[0])
            return (pe[0], pe[1])
        return None
    
    # Goes through every edge, getting the next.
    e = new_line_edge()
    lines = []
    while ( e is not None ):
        lines.append(list(e))
        remove_edge(e)
        ne = next_edge(e)
        while( ne is not None ):
            lines[-1].append(ne[1])
            remove_edge(ne)
            ne = next_edge(ne)
        e = new_line_edge()
    return lines

def signed_polygon_area( polygon, points ):
    area = 0.0
    for i in range( len( polygon ) ):
        ni = (i+1)%len(polygon)
        area += ( points[polygon[i]][0]*points[polygon[ni]][1] - points[polygon[i]][1]*points[polygon[ni]][0] )
    return area/2.0

def tie_polygons( lines, points ):
    """
    It creates a set of polygons from a set of lines, separating them and ordering using
    a right hand rule.
    Parameters
    ----------
    lines : list
        The set of lines returned by separate lines.
    Results
    -------
    polygons: list
        The set of polygons as point, list.
    graph_conn: list
        The lines surounding the polygons
    graph_dual: list
        The polygons surounding the lines.
    """
    # Separate loops in two, let's try to avoid a bug.
    nlines = len(lines)
    for i in range( nlines ):
        if lines[i][0] == lines[i][-1]:
            leni = len(lines[i])
            # if leni <= 2:
            #     print lines[i]
            #     assert False, "A loop is only two edges or less"
            lines.append(lines[i][(leni//2):])
            lines[i] = lines[i][:(leni//2)+1]
    
    ends = {}
    for i, l in enumerate(lines):
        
        if l[-1] in ends:
            ends[l[-1]].append(( i, True))
        else:
            ends[l[-1]] = [(i, True)]
        
        if l[0] in ends:
            ends[l[0]].append(( i, False))
        else:
            ends[l[0]] = [(i, False)]
    
    # Set of unclassified lines.
    unclass = set([(i, True) for i in range(len(lines))] + [(i, False) for i in range(len(lines))])
   
    # Eliminate all lines that don't suround anything
    to_remove = set()
    for e, ls in ends.items():
        if len(ls) == 1:
            l = ls[0]
            ln = lines[l[0]]
            to_remove.add(( ln[-1], l[0]))
            to_remove.add(( ln[0], l[0]))
            
            if (l[0], True) in unclass:
                unclass.remove((l[0], True))
                unclass.remove((l[0], False))
    
    while ( len(to_remove) ):
        e, l = to_remove.pop() 
        if not e in ends:
            continue
        
        lres = []
        for li in ends[e]:
            if li[0] != l:
                lres.append(li)
    
        if len(lres) > 1:
            ends[e] = lres
        else:
            ends.pop(e)
        if len(lres) == 1:
            l = lres[0]
            ln = lines[l[0]]
            to_remove.add(( ln[-1], l[0]))
            to_remove.add(( ln[0], l[0]))
            if (l[0], True) in unclass:
                unclass.remove((l[0], True))
                unclass.remove((l[0], False))
    
    def next_lines( nd, lidx ):
        vs = []
        
        for l, end in lidx:
            if end:
                v = np.array(points[nd])-np.array(points[lines[l][-2]])
            else:
                v = np.array(points[nd])-np.array(points[lines[l][1]])
            vs.append( ( l, v, end ) )
        
        vsang = list(map( lambda v: vector_angle(vs[0][1], -v[1]), vs[1:] ))
        vsord = [0] + sorted( [i for i in range(1, len(vs))], key=lambda i: vsang[i-1] )
        
        nxt = []
        for i in range(len(lidx)):
            nxt.append( lidx[vsord[i]] )
        return nxt
    
    # Generate the nexts dictionary.
    nexts = {}
    for e, ls in ends.items():
        lo = list(map( lambda l: l[0], ls ))
        # if len(lo) != len(set(lo)):
        #     print e
        #     print lo
        #     print map( lambda l: lines[l], lo )
        #     assert False
        # if len(ls) <= 1:
        #     print ls
        #     assert False
        nexts[e] = next_lines( e, ls )

    def next_line( nend, cl ):
        for i, n in enumerate(nend):
            if n == cl:
                for j in range( 1, len(nend) ):
                    nx = nend[(i+j)%len(nend)]
                    cl = ( nx[0], not nx[1] )
                    return cl
                assert False
        assert False
    # Tie all the polygons.
    all_polygons = []
    graph_dual = [[] for i in range(len(lines))]
    graph_conn = []
    cpol = 0
    while ( len( unclass ) ):
        cl = unclass.pop() # The current line.
        poly_conn = [cl[0]]
        graph_dual[cl[0]].append(cpol)
        stcl = cl
        if cl[1]:
            polygon = lines[cl[0]][:]
        else:
            polygon = lines[cl[0]][::-1]
        
        end = polygon[-1]
        while ( True ):
            pr = cl
            cl = next_line( nexts[end], cl )
            graph_dual[cl[0]].append(cpol)
            poly_conn.append(cl[0])
            if cl == stcl:
                break
            if cl[1]:
                polygon += lines[cl[0]][1:]
            else:
                polygon += (lines[cl[0]][::-1])[1:]
            try:
                unclass.remove(cl)
            except Exception as e:
                print( pr, nexts[end], cl, lines[cl[0]], polygon )
                raise e
            end = polygon[-1]
        
        all_polygons.append( polygon )
        graph_conn.append(poly_conn)
        cpol += 1
 
    return ( all_polygons, graph_conn, graph_dual )

def containments_from_to( polygons, contain, contained, points ):
    shcontain = []
    tree = rtree.index.Index()
    for i, n in enumerate(contain):
        ptpol = [points[p] for p in polygons[n]]
        shpol = Polygon(ptpol)
        shcontain.append( shpol )
        tree.insert( i, shpol.bounds )
    
    containments = [ [] for p in contain ]
    for i, n in enumerate(contained):
        pol = Polygon([ points[p] for p in polygons[n] ])
        for j in tree.intersection( pol.bounds ):
            if shcontain[j].contains(pol):
                containments[j].append(i)
                break
    return containments

def containments_all( polygons, to_search, points ):
    """
    Given a set of polygons with points, it finds which polygons have which holes, and substracts
    them from them, to return holed polygons.
    """
    shpolygons = []
    tree = rtree.index.Index()
    for i in to_search:
        ptpol = [points[n] for n in polygons[i]]
        shpol = Polygon(ptpol)
        
        # assert (not len(shpolygons)) or shpol.area <= shpolygons[-1].area
        shpolygons.append( shpol )
        
    containments = [[] for p in to_search]
    
    for i in range(len(shpolygons)):
        pos = []
        for j in tree.intersection( shpolygons[i].bounds ):
            if shpolygons[j].contains(shpolygons[i]):
                containments[j].append(i)
        if len( pos ):
            containments[max(pos)].append(i)
        tree.insert(i, shpolygons[i].bounds)
    
    return containments

def topology_relations( polygons, graph_conn, graph_dual, points ):
    """
    From the polygons, and its connnectivity, it creates a set of
    polygons with holes, with areas, and with a graph of connectivity.
    Parameters
    ----------
    polygons:
        The polygons returned by tie_polygons.
    graph_conn:
        The polygon to edge connectivity.
    graph_dual:
        Edge to polygon connectivity.
    Results ( tuple )
    -----------------
    polygons:
        The set of polygons with holes.
    graph:
        The graph filtered.
    parent:
        The body that's a parent of this one.
    parent_info:
        The basic information of the given parent.
    """
    
    # Create a direct graph where the polygon is connected to other polygons.
    graph = []
    for i, conn in enumerate(graph_conn):
        graph.append([])
        for l in conn:
            for k in graph_dual[l]:
                if k != i:
                    graph[-1].append(k)
    
    del graph_conn
    del graph_dual
    pos_polygons = []
    neg_polygons = []
    
    areas = []
    # Polygons with positive are what remains,
    # Polygons with negative area are either 
    # the whole covering or the polygon holes.
    for i, polygon in enumerate(polygons):
        area = signed_polygon_area( polygon, points )
        if area >= 0.0:
            pos_polygons.append(i)
            areas.append(area)
        else:
            neg_polygons.append(i)
            areas.append(-area)
    
    # By sorting the negative area polygons we are able to 
    # find which is the covering that contains the rest.
    neg_polygons.sort(key=lambda p: areas[p], reverse=True)
    
    # Find which polygon is a hole to another polygon.
    parent = [-1 for i in range(len(polygons))]
    coverings = [[] for i in neg_polygons]
    
    # Make parent array of all the polygons that can be reached by neg_polygons, (holes or covering).
    for i, polygon in enumerate(neg_polygons):
        seeds = [polygon]
        while len(seeds):
            cur = seeds.pop()
            if not parent[cur] == -1:
                continue
            parent[cur] = i
            coverings[i].append(cur)
            for p in graph[cur]:
                seeds.append(p)
    
    # Add the holes to its respective polygon.
    cover_contains = containments_all( polygons, neg_polygons, points )    
    holed_polygons = [[polygon] for polygon in polygons]
    
    for ri, contain in enumerate(reversed(cover_contains)):
        i = len(cover_contains)-(1+ri)
        inside = coverings[i][1:]
        spec_conts = containments_from_to( polygons, inside, [ neg_polygons[x] for x in contain ], points )
        for j, conts in enumerate(spec_conts):
            outpol = inside[j]
            for hole in conts:
        
                # Organize graph.
                inpol = neg_polygons[contain[hole]]
                parent[inpol] = i
                for k in graph[inpol]:
                    parent[k] = i
                    if graph[k] == inpol:
                        graph[k] = outpol
                graph[outpol] += graph[inpol]
                graph[inpol] = []
                # assert polygons[inpol] is not None
                holed_polygons[outpol].append(polygons[inpol])
                holed_polygons[inpol] = None
    
    
    return ( holed_polygons, areas, graph, parent, neg_polygons ) 
    
def reduce_everything( holed, areas, graph, parent, all_parents, points ):
    len_start = len(holed)
    len_parents = len(all_parents)
    rem_parents = sorted(list(set(parent))) # Fuck the police.
    
    parent_info = [ (i, areas[all_parents[i]]) for i in rem_parents ]
    
    trans = {}
    for i, p in enumerate(parent_info):
        trans[p[0]] = i
    parent_info = [ p[1] for p in parent_info ]
    for i in range(len(parent)):
        parent[i] = trans[parent[i]]
    
    # Reduce polygons, removing neg polygons.
    trans = [ -1 for i in range(len(holed))]
    to_rem = set( all_parents )
    cnt = 0
    for i in range(len( holed )):
        if not i in to_rem:
            trans[i] = cnt
            cnt += 1
    
    for i in range(len(graph)):
        if i in to_rem:
            graph[i] = None
            continue
        rem = []
        for j in range(len(graph[i])):
            if not graph[i][j] in to_rem:
                rem.append(trans[graph[i][j]])
        graph[i] = rem
    
    rever = [ i for i, c in enumerate(trans) if c >= 0 ]
    holed  = [ holed[i] for  i in rever ]
    areas  = [ areas[i] for  i in rever ]
    graph  = [ graph[i] for  i in rever ]
    parent = [ parent[i] for i in rever ]
    
    # Reduce points from all polygons.
    rem_points = set()
    for pol in holed:
        # assert pol is not None
        for ring in pol:
            for p in ring:
                rem_points.add(p)
    
    rem_points = sorted(list(rem_points))
    
    trans = [ -1 for i in range(len(points)) ]
    for i, n in enumerate(rem_points):
        trans[n] = i
    
    for pol in holed:
        for ring in pol:
            for i in range(len(ring)):
                # assert trans[ring[i]] >= 0
                ring[i] = trans[ring[i]]
    # Reduce the graphs to a single polygon (edges don't exist anymore, so order is unnecessary and can't really know anything about it).
    graph = [ list(set(g)) for g in graph ]
    # assert (len(holed) + len_parents) == len_start
    return ( holed, areas, graph, parent, parent_info, [points[i] for i in rem_points] )

def obtain_polygons( edges, points ):
    """
    Given a graph in 2D space, with their points attached,
    it obtains a set of polygons in the given graph.
    It also discards the edges that don't suround anything.
    """
    if type(points) != np.array:
        points = np.array(points)
    lines = separate_lines( edges )
    polygons, conn, dual = tie_polygons( lines, points )
    holed, areas, graph, parent, all_parents = topology_relations( polygons, conn, dual, points )
    holed, areas, graph, parent, parent_info, points = reduce_everything( holed, areas, graph, parent, all_parents, points )
    holed = list(map( lambda p: list(map( lambda r: r[:-1], p )), holed ))
    return ( holed, graph, list(map( tuple, points )) )

