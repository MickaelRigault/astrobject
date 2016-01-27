#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module, based on Shapely, enable to have access to contours
   and Patch for astrobjects"""

import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.patches import Polygon
try:
    from shapely.geometry import MultiPoint,Point, polygon, multipolygon
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

def point_in_contours(x,y,contours):
    """
    """
    if not HAS_SHAPELY:
        raise ImportError(_ERRORMESSAGE)
    if "__iter__" not in dir(x):
        return contours.contains(Point([x,y]))
    return contours.contains(MultiPoint(np.asarray([x,y]).T))
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
    return Polygon(polygon_to_vertices(polygon_), **kwargs)


def draw_polygon(ax,polygon_,**kwargs):
    """
    """
    from .tools import kwargs_update

    defaultprop = dict(fc=mpl.cm.Blues(0.6,0.1),
                       ec=mpl.cm.binary(0.8,0.9),lw=2)
    
    prop = kwargs_update(defaultprop,**kwargs)
    
    if type(polygon_) is polygon.Polygon:
        patch = polygon_to_patch(polygon_,**prop)
        ax.add_patch(patch)
    elif type(polygon_) is multipolygon.MultiPolygon:
        for p_ in polygon_:
            draw_polygon(ax,p_,**prop)
    else:
        raise TypeError("Only Polygon implemented")
