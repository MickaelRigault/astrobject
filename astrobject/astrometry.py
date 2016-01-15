#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module manage the Astrometry issue"""

import numpy as np
from astLib      import astWCS


def get_wcs(wcsinput,verbose=True):
    """
    """
    if wcsinput is None:
        return None
    
    # ----------------------
    # - astLib.astWCS Format
    if wcsinput.__module__ == "astLib.astWCS":
        return convert_to_wcs(wcsinput)
    
    # --------------
    # - Good format
    if "__nature__" in dir(wcsinput) and wcsinput.__nature__ == "WCS":
        return wcsinput
    if verbose:
        print "WARNING: only astLib.astWCS / astrobject.astrometry.WCS wcs solution are implemented"
        print " ----> No wcs solution"
        
    return None

def astlibwcs_to_wcs(astlibwcs):
    """
    """
    if "headerSource" not in dir(otherformat):
        raise TypeError("The given 'astlibwcs' is not an astLib.WCS")
    
    return WCS(otherformat.headerSource)

##########################################
#                                        #
# Advance WCS Astrobject                 #
#                                        #
##########################################
class WCS(astWCS.WCS):
    """
    """
    __nature__ = "WCS"
    
    
    @property
    def contours(self,exp_order=5):
        # -- Load it the fist time you use it
        if "_contour" not in dir(self):
            from ..utils import shape
            if not shape.HAS_SHAPELY:
                self._contour = None
            else:
                # This could be improved
                npoints = 2+exp_order
                # -- This defines the contours
                width = np.linspace(0,self.header["NAXIS1"],npoints)
                heigh = np.linspace(0,self.header["NAXIS2"],npoints)
                v1 = np.asarray([np.ones(npoints-1)*0, width[:-1]]).T
                v2 = np.asarray([heigh[:-1], np.ones(npoints-1)*self.header["NAXIS1"]]).T
                v3 = np.asarray([np.ones(npoints-1)*self.header["NAXIS2"], width[::-1][:-1]]).T
                v4 = np.asarray([heigh[::-1][:-1], np.ones(npoints-1)*0]).T
                ra,dec = np.asarray([self.pix2wcs(i,j)
                                for i,j in np.concatenate([v1,v2,v3,v4],axis=0)]).T
                
                self._contour = shape.get_contour_polygon(ra,dec)
            
        return self._contour

    def has_contours(self):
        """Test if this have contours defined"""
        return self.contours is not None
