#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module, based on Shapely, enable to have access to contours
   and Patch for astrobjects"""
import warnings
import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.patches import Polygon

# -- Local module
from .decorators import make_method

# -- Shapely
try:
    from shapely.geometry import MultiPoint,Point, polygon, multipolygon, MultiPolygon
    HAS_SHAPELY = True
    
except ImportError:
    _ERRORMESSAGE = "you do not have shapely installed."
    warnings.warn(_ERRORMESSAGE + " Some function will not be accessible."+"\n"+"--> pip install Shapely")
    HAS_SHAPELY = False


def get_contour_polygon(x,y):
    """
    """
    if not HAS_SHAPELY:
        raise ImportError(_ERRORMESSAGE)
    
    points = MultiPoint(np.asarray([x,y]).T)
    return points.convex_hull


def point_in_contours(x,y,contours, all=True):
    """
    This methods check if the given coordinate (x,y)
    is within the contours, which could be a:
      - matplotlib.patches.Polygon ; or
      - Shapely.geometry.Polygon

    x and y could be array of coordinates

    if x and y are list of coordinates, if 'all' is True
    an single boolean value is returned, answering the qestion:
    'are all the given points in the contours?' Otherwise
    a list of array is returned for each individual points
    Return
    ------
    bool (bool-array)
    """
    if not HAS_SHAPELY:
        raise ImportError(_ERRORMESSAGE)

    # ----------------------
    # - Matplotlib Polygon
    if type(contours) is Polygon:
        if "__iter__" not in dir(x):
            return contours.contains_point([x,y])
        
        return np.asarray([contours.contains_point([x_,y_])
                           for x_,y_ in zip(x,y)],dtype=bool)
    
    # ----------------------
    # - Shapely Polygon
    if "__iter__" not in dir(x):
        return contours.contains(Point(x,y))
    if all:
        return contours.contains(MultiPoint(np.asarray([x,y]).T))
    return [contours.contains(Point(x_,y_))
            for x_,y_ in zip(x,y) ]

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
    return np.concatenate([np.asarray(polygon_.exterior)] + [np.asarray(r) for r in polygon_.interiors])

def polygon_to_patch(polygon_,**kwargs):
    """
    **kwargs is the properties of the returned matplotlib polygon patch.
    """
    return Polygon(polygon_to_vertices(polygon_), **kwargs)

def patch_to_polygon(mpl_patch):
    """ converts a matplotlib patch into a shapely polygon. If mpl_patch is a list
     of patches, the returned polygon will be a 'union cascade' MultiPolygon  """
    if not hasattr(mpl_patch,"__iter__"):
        return get_contour_polygon(*mpl_patch.get_verts().T)
    
    from shapely.ops import cascaded_union
    return cascaded_union([get_contour_polygon(*p_.get_verts().T) for p_ in mpl_patch])
    

@make_method(mpl.Axes)
def draw_polygon(ax,polygon_,**kwargs):
    """
    draw the polygon within the current axis
    """
    from .tools import kwargs_update
    # ----------------
    # - Styling
    defaultprop = dict(fc=mpl.cm.Blues(0.6,0.1),
                       ec=mpl.cm.binary(0.8,0.9),lw=2)
    
    prop = kwargs_update(defaultprop,**kwargs)
    if prop["ec"] is None:
        prop["ec"] = prop["fc"]

    # ----------------
    # - Patching
    if type(polygon_) is polygon.Polygon:
        patch = polygon_to_patch(polygon_,**prop)
        return ax.add_patch(patch)
    elif type(polygon_) is multipolygon.MultiPolygon:
        return [draw_polygon(ax,p_,**prop) for p_ in polygon_]
    else:
        raise TypeError("Only Polygon implemented")



def get_voronoy_multipolygon(x, y ,  edges=None):
    """
    Parameters
    ----------
    rect: [xmin, xmax, ymin, ymax] -optional-
        Edges of the voronoy map
    """
    from scipy.spatial import Voronoi
    from itertools import product
    
    # - Define the Grid
    # --------------------
    flagok = np.isnan(x) * np.isnan(y)
    x,y    = x[~flagok], y[~flagok]
    # define extremum to avoid non finished edges
    ext_x,ext_y = np.asarray(list (product([x.min()-x.max()*10,x.mean(), x.max()+x.max()*10],
                                    [y.min()-y.max()*10,y.mean(), y.max()+y.max()*10]))).T

    xy = np.asarray([np.concatenate([x,ext_x]),
                     np.concatenate([y,ext_y])]).T
    
    npoint = np.shape(xy)[0]
    vor = Voronoi(xy)
    xy_poly = [[vor.vertices[x] for x in vor.regions[vor.point_region[i]]
                        if x>=0  # this others could be saved
                        ] for i in range(npoint)]

    
    #edges = polygon.Polygon([[0,0],[0,4100],[2100,4100],[2100,0]])
    polygons = []
    for xy_ in xy_poly:
        try:
            p_ = edges.intersection(polygon.Polygon(xy_))
        except:
            continue
        if p_.area>0:
            polygons.append(p_)
        
    return MultiPolygon(polygons)
    
