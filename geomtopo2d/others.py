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

import numpy as np
from numpy import linalg as la
import math

# THIS LIBRARY HAS UNCLASSIFIED FUNCTIONS.

def circle_from_points(p0, p1, p2, thresh_angle):
    """
    Finds a circle that passes through 3 points.
    """
    # Find the vectors.
    v1 = p1-p0
    v2 = p2-p0
    # p0 is in the middle, p1 is behind, p2 is afterwards.
    ang = math.acos( np.dot(v1/la.norm(v1), v2/la.norm(v2)) )
    if ang < thresh_angle:
        # This is a sharp edge.
        return ( -1.0, None )
    
    n = np.cross(v1, v2)
    if la.norm(n) < TOL:
        return ( float('inf'), None )
    
    # First parameter of the equation.
    y0 = np.dot(n, p0)

    # Find the other two parameters.
    m0 = np.dot(p0, p0)
    m1 = np.dot(p1, p1)
    m2 = np.dot(p2, p2)
    y1 = (m1 - m0) / 2.0
    y2 = (m2 - m0) / 2.0
    
    c = la.solve( np.array([n, v1, v2]), np.array([y0, y1, y2]) )
    return ( la.norm( c-p0 ) , c )

def pseudo_circle_from_points3(pts, nrm):
    """
    Finds an approximate circle that passes through n points.
    """
    
    # Find the vectors.
    vcts = np.zeros((pts.shape[0]-1, pts.shape[1]))
    for i in range(pts.shape[0]-1):
        vcts[i] = pts[i+1,:]-pts[0,:]
    
    # First parameter of the equation.
    y0 = np.dot(nrm, pts[0,:])
    
    # Find the other parameters.
    m = np.array([ np.dot(pts[i,:], pts[i,:]) for i in range(pts.shape[0]) ])
    y = (m[1:]-m[0])/2.0
    
    y = np.append(y, y0)
    A = np.append(vcts, [nrm], axis=0)
    try:
        x = la.lstsq( A, y )[0]
    except:
        return float('inf')
    
    return ( la.norm( x-pts[0,:] ) , x )

def pseudo_circle_from_points2(pts):
    """
    Finds an approximate circle that passes through n points.
    """
    # Find the vectors.
    A = np.zeros((pts.shape[0]-1, pts.shape[1]))
    for i in range(pts.shape[0]-1):
        A[i] = pts[i+1,:]-pts[0,:]
    
    # Find the other parameters.
    m = np.array([ np.dot(pts[i,:], pts[i,:]) for i in range(pts.shape[0]) ])
    y = (m[1:]-m[0])/2.0
    
    x = la.lstsq( A, y )[0]
    sr = 0.0
    for i in xrange(pts.shape[0]):
        sr += la.norm(pts[i,:]-x)
    sr /= pts.shape[0]
    dif = sr-la.norm( pts-x, axis=1 )
    sm = la.norm(dif)
    sm /= pts.shape[0]
    
    return (sr, x, sm )

def line_fit2( pts ):
    """
    Bit hacky linear squares fit for a line.
    """
    vct = pts[-1,:]-pts[0,:]
    nrm = la.norm(vct)
    # This means I pushed some "circular" thing that's flat A.F.
    if nrm < TOL:
        vct = pts[pts.shape[0]/2,:]-pts[0,:]
        vct /= la.norm(vct)
    else:
        vct /= la.norm(vct)
    ort = np.array([-vct[1], vct[0]])
    # Find the vectors.
    y = np.zeros((pts.shape[0]-2,))
    x = np.zeros((pts.shape[0]-2,))
    for i in range(0, pts.shape[0]-2):
        v = pts[i+1,:]-pts[0,:]
        x[i] = np.dot(v, vct)
        y[i] = np.dot(v, ort)
    A = np.column_stack([x, np.ones(x.shape)])
    
    m, c = la.lstsq(A, y)[0]
    
    cr = ort*c + pts[0,:]
    mr = ort*m + vct
    mr /= la.norm(mr)
    ortr = np.array([-mr[1], mr[0]])
    
    sm = np.zeros((pts.shape[0],))
    mnt = float('inf')
    mxt = -float('inf')
    for i in xrange(pts.shape[0]):
        p = pts[i,:]-cr
        sm[i] = np.dot( ortr, p )
        t = np.dot( mr, p )
        mnt = min(mnt, t)
        mxt = max(mxt, t)
    
    sm  = la.norm(sm)
    sm /= pts.shape[0]
    
    return (mr, cr, (mnt, mxt), sm)

def what_is_fit(pts):
    r, c, sc = pseudo_circle_from_points2(pts)
    ml, cl, tr, sl = line_fit2( pts )
    if sl < sc:
        return ("l", sl, ml, cl)
    return ("c", sc, r, c)

def triangles_angle( t1, t2 ):
    n1 = np.cross( t1[2]-t1[0], t1[1]-t1[0] )
    n2 = np.cross( t2[2]-t2[0], t2[1]-t2[0] )
    if la.norm(n1) < TOL or la.norm(n2) < TOL:
        return math.pi
    d1 = la.norm(n1)
    d2 = la.norm(n2)
    if d1 < TOL or d2 < TOL:
        return -1.0
    n1 /= d1
    n2 /= d2
    dp = np.dot(n1, n2)
    return math.pi - math.acos( min(dp, 1.0) )

def triangle_normal( t ):
    n = np.cross( t[2]-t[0], t[1]-t[0] )
    d = la.norm(n)
    if d < TOL:
        return None
    n /= d
    return n

def find_pseudo_normal( pts ):
    """
    Normal of 4 or more points.
    """
    vcts = np.zeros((pts.shape[0]-1, 3))
    for i in range(pts.shape[0]-1):
        vcts[i] = pts[i+1,:]-pts[0,:]
    csum = np.zeros(3)
    n = 0
    for i in range(vcts.shape[0]):
        c = np.cross(vcts[i], vcts[(i+1)%vcts.shape[0]])
        d = la.norm(c)
        if d < TOL:
            continue
        c /= d
        if la.norm(csum) > la.norm(csum+c):
            csum += -c
        else:
            csum += c
        n += 1
    
    if n == 0:
        return None
    
    csum /= n
    return csum

def find_perpendicular_basis( p0 ):
    mn = float('inf')
    mni = None
    for i in range(3):
        if abs(p0[i]) < mn:
            mn = abs(p0[i])
            mni = i
    
    o1 = (mni+1)%3
    o2 = (mni+2)%3
    
    p1 = np.zeros(3)
    p1[o1] = p0[o2]
    p1[o2] = -p0[o1]
    
    p1 /= la.norm(p1)
    
    p2 = np.cross( p0, p1 )
    p2 /= la.norm(p2)
    
    return p1, p2

def find_curvature( pts, n ):
    normal = find_pseudo_normal(pts)
    if normal is None:
        return ( float('inf'), None, None )
    
    p0 = normal - (np.dot(n, normal)*n)
    d0 = la.norm(p0)
    if d0 < TOL:
        return ( float('inf'), None, None )
    
    p0 /= d0
    
    p1, p2 = find_perpendicular_basis(p0)
    projpts = np.zeros(pts.shape)
    for i in xrange(pts.shape[0]):
        u = np.dot(p1, pts[i,:])
        v = np.dot(p2, pts[i,:])
        projpts[i,:] = u*p1 + v*p2
    
    r, c = pseudo_circle_from_points3(projpts, p0)
    
    return (r, c, p0)

def pseudo_3d_projection( points3d, points2d ):
    """
    Find the matrix that transforms the points in the cad model to the 
    [ a b c d ] [xm]   [xp]
    [ e f g h ] [ym] = [yp]
    [ 0 0 0 1 ] [zm]   [ 1]
                [ 1]
    points in the picture.
    """
    xp = points2d[:,0]
    yp = points2d[:,1]

    A = np.append( points3d, np.ones((points3d.shape[0],1)), axis=1 )
    abcd = la.lstsq( points3d, xp )[0]
    efgh = la.lstsq( points3d, yp )[0]
    matrix = np.array([abcd, efgh, np.array([0, 0, 0, 1], dtype=float)])
    return matrix


