#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" This module gather the add-on methods for matplotlib axes and figure"""

import numpy as np
import matplotlib.pyplot as mpl

from .tools import kwargs_update
from .decorators import make_method
from .skyplot import convert_radec_azel


__all__ = ["specplot","skyplot","figout"]

# ========================== #
# =  Axes Add-on           = #
# ========================== #

# --------------------- #
# - Spectrum          - #
# --------------------- #

@make_method(mpl.Axes)
def specplot(ax,x,y,var=None,
             color=None,bandprop={},**kwargs):
    """This function in a build-in axes method that enable to quickly and
    easily plot a spectrum.
    """
    # -----------------------
    # - Properties of plot
    default_kwargs = dict(
        color=mpl.cm.Blues(0.8),
        ls="-",lw=1,marker=None,zorder=6,
        )
    if color is not None:
        default_kwargs["color"] = color
    propplot = kwargs_update(default_kwargs,**kwargs)
    # -- Plot 
    pl = ax.plot(x,y,**propplot)
    
    # -----------------------
    # - Properties of band
    if var is not None:
        default_band   = dict(
            color=propplot["color"],alpha=0.3,
            zorder=3,label="_no_legend_"
            )
        bandprop = kwargs_update(default_band,**bandprop)
        # -- Band
        fill = ax.fill_between(x,y+np.sqrt(var),y-np.sqrt(var),
                        **bandprop)
    else:
        fill = None
        
    return pl,fill

# --------------------- #
# - Skyplot           - #
# --------------------- #

@make_method(mpl.Axes)
def skyplot(ax, ra, dec, var=None,
            color=None, bandprop={}, **kwargs):
    """This function in a build-in axes method that allows easily plotting points 
    on the sky.
    """
    # -----------------------
    # - Properties of plot
    default_kwargs = dict(marker='o', markersize=5, linestyle='none')
    if color is not None:
        default_kwargs["color"] = color
    propplot = kwargs_update(default_kwargs,**kwargs)

    az, el = convert_radec_azel(ra, dec)
    # -- Plot 
    pl = ax.plot(az, el, **propplot)
    
    return pl


@make_method(mpl.Axes)
def skyscatter(ax, ra, dec, var=None,
               bandprop={}, **kwargs):
    """This function in a build-in axes method that allows simple scatter plots
    easily plot a spectrum.
    """
    # -----------------------
    # - Properties of plot
    default_kwargs = dict(marker='o', s=30, #edgecolor='none',
                          cmap=mpl.cm.RdYlBu_r)

    propplot = kwargs_update(default_kwargs,**kwargs)

    az, el = convert_radec_azel(ra, dec)
    # -- Plot 
    sc = ax.scatter(az, el, **propplot)
    
    return sc


# --------------------- #
# - Map Plots         - #
# --------------------- #
@make_method(mpl.Axes)
def voronoi_patchs(ax, xy, c=None,vmax=None,vmin=None,
                   cmap=mpl.cm.jet,**kwargs):
    """
    """
    from scipy.spatial import Voronoi
    from matplotlib.patches import Polygon
    from matplotlib.collections import PolyCollection
    
    npoint = np.shape(xy)[0]
    vor = Voronoi(xy)

    def _test_regions_match_(r1,r2):
        return (r1,r2) in vor.ridge_dict.keys() or (r2,r1) in vor.ridge_dict.keys()

    
    xy_poly = [[vor.vertices[x] for x in vor.regions[vor.point_region[i]]
                        if x>=0  # this others could be saved
                        ] for i in range(npoint)]
    if c is not None:
        vmin = c.min() if vmin is None else vmin
        vmax = c.max() if vmax is None else vmax
        color = cmap((c-vmin)/(vmax-vmin))
    else:
        color = "None"
        
    collec = PolyCollection(xy_poly,facecolors=color,**kwargs)
    ax.add_collection(collec)    

    
def voronoi_grid(ax,xy):
    """
    TO BE REMOVED
    """
    def _adjust_bounds(ax, points):
        ptp_bound = points.ptp(axis=0)
        ax.set_xlim(points[:,0].min() - 0.1*ptp_bound[0],
                    points[:,0].max() + 0.1*ptp_bound[0])
        ax.set_ylim(points[:,1].min() - 0.1*ptp_bound[1],
                    points[:,1].max() + 0.1*ptp_bound[1])

    from scipy.spatial import Voronoi
    from matplotlib.patches import Polygon
    vor = Voronoi(xy)

    
    if vor.points.shape[1] != 2:
        raise ValueError("Voronoi diagram is not 2-D")

    ax.plot(vor.points[:,0], vor.points[:,1], '.')

    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            ax.plot(vor.vertices[simplex,0], vor.vertices[simplex,1], 'k-')
                        
    ptp_bound = vor.points.ptp(axis=0)

    center = vor.points.mean(axis=0)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.any(simplex < 0):
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex
            print i,simplex, simplex[simplex >= 0][0]
            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            ax.plot([vor.vertices[i,0], far_point[0]],
                    [vor.vertices[i,1], far_point[1]], 'k--')

    _adjust_bounds(ax, vor.points)
        

# ========================== #
# =  Figure Add-on         = #
# ========================== #
@make_method(mpl.Figure)
def figout(fig,savefile=None,show=True,add_thumbnails=False,
           dpi=200):
    """This methods parse the show/savefile to know if the figure
    shall the shown or saved."""
    
    if savefile in ["dont_show","_dont_show_","_do_not_show_"]:
        show = False
        savefile = None

    if savefile is not None:
        fig.savefig(savefile+'.png',dpi=dpi)
        fig.savefig(savefile+'.pdf')
        if add_thumbnails:
            fig.savefig(savefile+"_thumb"+'.png',dpi=dpi/10.)
            
    elif show:
        fig.canvas.draw()
        fig.show()
        
