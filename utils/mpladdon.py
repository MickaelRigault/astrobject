#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" This module gather the add-on methods for matplotlib axes and figure"""

import numpy as np
import matplotlib.pyplot as mpl

from .tools import kwargs_update
from .decorators import make_method


__all__ = ["specplot","figout"]

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
        
