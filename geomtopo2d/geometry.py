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

import numpy as np
from numpy import linalg as la
import math

TOL = np.finfo(np.float32).eps

def vector_angle( v1, v2 ):
    v1 = np.array(v1, dtype=float)
    v1 /= la.norm(v1)
    v2 = np.array(v2, dtype=float)
    v2 /= la.norm(v2)
    c = np.cross(v1, v2)
    d = np.dot(v1, v2)
    ret = (math.pi - math.atan2(c, d)) % (2*math.pi)
    return ret

def vector_angle3( v1, v2 ):
    v1 = np.array(v1, dtype=float)
    nv1 = la.norm( v1 )
    if nv1 < TOL:
        return math.pi
    v1 /= nv1
    
    v2 = np.array(v2, dtype=float)
    nv2 = la.norm( v2 )
    if nv2 < TOL:
        return math.pi
    v2 /= nv2
    
    c = la.norm(np.cross(v1, v2))
    d = np.dot(v1, v2)
    ret = (math.pi - math.atan2(c, d)) % (2*math.pi)
    return ret
    
