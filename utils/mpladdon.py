#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" This module gather the add-on methods for matplotlib axes and figure"""

import numpy as np
import matplotlib.pyplot as mpl

from .tools import kwargs_update
from .decorators import make_method
# remark: do no import .plot.*** at this level
#         otherwise import loop with astrobject
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
def skyplot(ax, ra, dec, color=None, **kwargs):
    """
    Build-in axes method for easy plotting of points on the sky.

    ax: [matplotlib.axes]          where the data will be displayed
                                   should be using Mollweide or Hammer
                                   projection (not strictly enforced)

    ra, dec: [N array, N array]    arrays of sky coordinates (RA, Dec)
                                   in degrees

    - options -

    color: [string]                color to be used for the plot

    - kwargs goes to ax.plot -

    Return
    ------
    pl (output of ax.plot)
    """
    from .plot.skyplot import convert_radec_azel
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
def skyscatter(ax, ra, dec, **kwargs):
    """
    Build-in axes method for scatter plots of points on the sky.

    ax: [matplotlib.axes]          where the data will be displayed
                                   should be using Mollweide or Hammer
                                   projection (not strictly enforced)

    ra, dec: [N array, N array]    arrays of sky coordinates (RA, Dec)
                                   in degrees

    - kwargs goes to ax.scatter -

    Return
    ------
    sc (output of ax.scatter)
    """
    from .plot.skyplot import convert_radec_azel
    # -----------------------
    # - Properties of plot
    default_kwargs = dict(marker='o', s=30,
                          cmap=mpl.cm.RdYlBu_r)

    propplot = kwargs_update(default_kwargs,**kwargs)

    az, el = convert_radec_azel(ra, dec)
    # -- Plot 
    sc = ax.scatter(az, el, **propplot)
    
    return sc

@make_method(mpl.Axes)
def skyhist(ax, ra=None, dec=None, values=None, bins=None, steps=None, max_stepsize=5, edge=1e-6,
            vmin=None, vmax=None, cmap=mpl.cm.Blues, cblabel=None, **kwargs):
    """This function in a build-in axes method that makes a sky histogram of the 
    coordinates.
    
    """
    """
    Build-in axes method for 2d-histograms of points on the sky.

    ax: [matplotlib.axes]          where the data will be displayed
                                   should be using Mollweide or Hammer
                                   projection (not strictly enforced)

    - options -

    ra, dec: [N array, N array]    arrays of sky coordinates (RA, Dec)
                                   in degrees

    values: [N array]              array of values for each bin (can be used 
                                   instead of ra and dec; must number of bins as
                                   length)
    
    bins [Bins object]             object that bins the coordinates and
                                   provides the boundaries for drawing
                                   polygons 
                                   (see the documentation of 
                                   utils/plot/skybins.py)

    steps [int]                    number of steps between bin verices
                                   (if None, it will be determined based on
                                   max_stepsize)
    
    max_stepsize [float]           maximal stepsize used to determined
                                   steps if None
    
    edges [float]                  edge to be left near RA = 180 deg

    vmin,vmax: [float,float]       minimal and maximal values for the colormapping.

    cmap: [mpl.cm]                 a colormap
    
    clabel: [string]               colorbar label

    - kwargs goes to matplotlib.collections.PolyCollection -

    Return
    ------
    collection (output of matplotlib.collections.PolyCollection)
    cbar       (output of ax.colorbar)
    """
    # -----------------------
    # - Properties of plot
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection, PolyCollection

    from .plot.skybins import SkyBins
    from .plot.skyplot import convert_radec_azel
    
    if bins is None:
        bins = SkyBins()

    if clabel is None:
        clabel = ''

    if values is None:
        values = bins.hist(ra, dec)
    else:
        values = np.array(values)

    patches = []
    p_idx = []
    for k in xrange(len(values)):
        radec_bd = bins.boundary(k, steps=steps, max_stepsize=max_stepsize,
                                 edge=edge)
        for r, d in radec_bd:
            coord_bd = np.asarray(convert_radec_azel(r, d, edge=edge)).T
            patches.append(coord_bd)
            p_idx.append(k)

    c = np.asarray(values[np.asarray(p_idx)])
    vmin = c.min() if vmin is None else vmin
    vmax = c.max() if vmax is None else vmax
    color = cmap((c-vmin)/float(vmax-vmin))

    collec = PolyCollection(patches, facecolors=color, **kwargs)
    collec.set_edgecolor('face')
    
    # -- Plot 
    ax.add_collection(collec)
    
    axcar = ax.insert_ax(fraction=0.9, space=.08, pad=0.05, location='bottom')
    cbar = axcar.colorbar(cmap, vmin=vmin, vmax=vmax, label=cblabel)

    return collec, cbar


# --------------------- #
# - WCS Plots         - #
# --------------------- #
@make_method(mpl.Axes)
def wcsplot(ax, wcs, exp_order=10,
            fc=mpl.cm.Blues(0.6,0.1),ec=mpl.cm.binary(0.8,1),
            draw_corner=False,
            **kwargs):
    """
    Comments TO BE DONE
    """
    from matplotlib.patches import Polygon
    # -----------------
    # - verticles
    if "has_contours" not in dir(wcs) or not wcs.has_contours():
        npoints = 2+exp_order
        width = np.linspace(0,wcs._naxis1,npoints)
        heigh = np.linspace(0,wcs._naxis2,npoints)
        v1 = np.asarray([np.ones(npoints-1)*0, width[:-1]]).T
        v2 = np.asarray([heigh[:-1], np.ones(npoints-1)*wcs._naxis1]).T
        v3 = np.asarray([np.ones(npoints-1)*wcs._naxis2, width[::-1][:-1]]).T
        v4 = np.asarray([heigh[::-1][:-1], np.ones(npoints-1)*0]).T
        v = np.asarray([wcs.pix2wcs(i,j)
                        for i,j in np.concatenate([v1,v2,v3,v4],axis=0)])
    else:
        from .shape import polygon_to_vertices
        v = polygon_to_vertices(wcs.contours)
        
    poly = Polygon(v,fc=fc,ec=ec,lw=1,**kwargs)
    # ------------------
    # - Draw
    # The point used
    pl = ax.plot(v.T[0],v.T[1],ls="None",marker="o",mfc=fc,mec=ec,
            visible=draw_corner)
    # The actual Patch
    ax.add_patch(poly)
    # ------------------
    # - Returns
    return pl, poly
    
# --------------------- #
# - Patchs Plots      - #
# --------------------- #
@make_method(mpl.Axes)
def voronoi_patchs(ax, xy, c=None, vmax=None, vmin=None,
                   cmap=mpl.cm.jet,
                   cbar=True,cblabel="",cbarprop={},
                   **kwargs):
    """
    ax: [matplotlib.axes]          where the data will be displaid

    xy: [N-2D array]               the combination upon which the voronoi tesselation
                                   will be built. (scipy.spatial.Voronoi)
                                   
    - options -
    
    c: [array]                     an array of value for the color of the patchs. You
                                   can also provide 1 float between 0 and 1 (see cmap)

    vmin,vmax: [float,float]       minimal and maximal values for the colormapping.
                                   If None the c-array's minimal/maximal value will be
                                   set.
    - other options -
    
    cmap: [mpl.cm]                 a colormap
    
    cbar: [bool or ax]             provide here an ax where the colorbar should be
                                   drawn. You can also set True to have a default one
                                   or set False to avoid having a colorbar.

    
    - kwargs goes to matplotlib.collections.PolyCollection -                        

    Return
    ------
    collection (output of matplotlib.collections.PolyCollection)
    """
    from scipy.spatial import Voronoi
    from matplotlib.collections import PolyCollection
    
    # --------------------
    # - Define the Cells
    npoint = np.shape(xy)[0]
    vor = Voronoi(xy)
    xy_poly = [[vor.vertices[x] for x in vor.regions[vor.point_region[i]]
                        if x>=0  # this others could be saved
                        ] for i in range(npoint)]

    # --------------------
    # - Color of the Cells
    if c is not None:
        if "__iter__" not in dir(c):
            c = [c]*npoint
        else:
            c = np.asarray(c)
            vmin = c.min() if vmin is None else vmin
            vmax = c.max() if vmax is None else vmax
            color = cmap((c-vmin)/(vmax-vmin))
        edgecolors = kwargs.pop("edgecolors","0.5")
    else:
        color = "None"
        edgecolors = kwargs.pop("edgecolors","k")
        

    collec = PolyCollection(xy_poly,facecolors=color,edgecolors=edgecolors,
                            **kwargs)
    ax.add_collection(collec)

    # ----------------- #
    # - color bar     - #
    # ----------------- #
    _display_colorbar_ = not (color is "None" or not cbar)
    if _display_colorbar_:
        # - this means it is not an ax
        if "imshow" not in dir(cbar):
            axcar = ax.insert_ax(space=.05,pad=0.03,location="right")
        else:
            axcar = cbar
        calpha = cbarprop.pop('alpha',kwargs.pop("alpha",None))
        return collec, axcar.colorbar(cmap,vmin=vmin,vmax=vmax,label=cblabel,
                                      alpha=calpha,**cbarprop)
    # - no cbar
    return collec
    
# --------------------------- #
# -    Color Bar            - #
# --------------------------- #
@make_method(mpl.Axes)
def colorbar(ax,cmap,vmin=0,vmax=1,label="",
             fontsize="large",npoint=256,
            **kwargs):
    """
    """
    import matplotlib
    if "orientation" not in kwargs.keys():
        bbox = ax.get_position()
        orientiation = "vertical" if bbox.xmax - bbox.xmin < bbox.ymax - bbox.ymin \
          else "horizontal"
        kwargs["orientation"] = orientiation

    norm    = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    c_bar   = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                              norm=norm,**kwargs)
    
    c_bar.set_label(label,fontsize=fontsize)
    if "ticks" in kwargs.keys() and "ticklabels" not in kwargs.keys():
        c_bar.ax.set_xticklabels([r"%s"%v for v in kwargs["ticks"]])
        
    ax.set_xticks([]) if kwargs["orientation"] == "vertical" \
      else ax.set_yticks([])
    return c_bar

# ========================== #
# =  axes manipulation     = #
# ========================== #
@make_method(mpl.Axes)
def insert_ax(ax,space,pad,fraction=None,
              location="right"):
    """
    Split the given axis to insert a new one at this location.
    This is typically useful to add colorbar for instance.

    Parameters
    ----------
    
    ax: [mpl.axes]                  The axis that will be split.

    space: [float]                  Extra space between the *ax* and the *newax* that
                                    will be created from the freed space.

    pad: [float]                    Extra space at the right/top of the *newax*

    fraction: [float /None]         The location where the axes will be
                                    split, given in faction of *ax* size.
                                    e.g., 0.5 means in the middle, while
                                    0.8 is at the far right / top (see *location*)
                                    If None: the new axis + given axis will fill exactly
                                    the same space as the given ax.
    - options -
    
    location: [string]              where you whish to add the *newax*. Only
                                    right and top are defined.

    Return
    ------
    axes (the *newax*, the input *ax* is reshape)
    """
    known_location = ["right","top","left","bottom"]
    if location not in known_location:
        raise ValueError("'%s' is not a correct location"%location," These are: "+", ".join(known_location))

    if fraction is None:
        fraction = 1- (space + pad)
        
    if location == "right":
        bboxax,_,bboxnew,_ = ax.get_position().splitx(fraction,fraction+space,
                                                      fraction+space+pad)
    elif location == "left":
        bboxnew,_,bboxax,_ = ax.get_position().splitx(pad,pad+space,
                                                      fraction+space+pad)
    elif location == "top":
        bboxax,_,bboxnew,_ = ax.get_position().splity(fraction,fraction+space,
                                                      fraction+space+pad)
    else: # bottom
        bboxnew,_,bboxax,_ = ax.get_position().splity(pad,pad+space,
                                                      fraction+space+pad)

    ax.set_position(bboxax)
    return ax.figure.add_axes(bboxnew)

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
