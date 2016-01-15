#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module, based on Shapely, enable to have access to contours
   and Patch for astrobjects"""

import numpy as np
try:
    from shapely.geometry import MultiPoint, polygon
    HAS_SHAPELY = True
except ImportError:
    _ERRORMESSAGE = "you do not have shapely installed."
    print _ERRORMESSAGE + " Some function will not be accessible."
    print "--> pip install Shapely"
    HAS_SHAPELY = False


def get_contour_polygon(x,y):
    """
    """
    if not HAS_SHAPELY:
        raise ImportError(_ERRORMESSAGE)
    
    points = MultiPoint(np.asarray([x,y]).T)
    return points.convex_hull

def polygon_to_vertices(polygon_):
    """
    """
    # ----------------
    # - Input Test
    if not HAS_SHAPELY:
        raise ImportError(_ERRORMESSAGE)
    if type(polygon_) != polygon.Polygon:
        raise TypeError("the given polygon must be a shapely polygon")
    
    # ----------------
    # - Convertion
    return np.concatenate([np.asarray(polygon_.exterior)]
                          + [np.asarray(r) for r in polygon_.interiors])

def polygon_to_patch(polygon_,**kwargs):
    """
    **kwargs is the properties of the returned matplotlib polygon patch.
    """
    from matplotlib.patches import Polygon
    return Polygon(polygon_to_vertices(polygon_), **kwargs)

