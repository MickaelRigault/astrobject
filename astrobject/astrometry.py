#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module manage the Astrometry issue"""

import numpy as np
from ..utils import shape
from astLib      import astWCS
from astropy.wcs import WCS as pyWCS
from astropy.io  import fits
from astropy.coordinates.angle_utilities import angular_separation

def wcs(filename=None,header=None,extension=0):
    """
    """
    if filename is not None:
        header = fits.getheader(filename,ext=extension)
        return WCS(header)
    if header is None:
        raise ValueError("'filename' or 'header' must be given")
    return WCS(header)

# =========================== #
# =
# =========================== #
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
    
    return wcs(otherformat.headerSource)

##########################################
#                                        #
# Advance WCS Astrobject                 #
#                                        #
##########################################
class WCS(pyWCS):
    """
    """
    __nature__ = "WCS"
    
    # -------------------- #
    # - Mimic astLib     - #
    # -------------------- #
    def pix2wcs(self,x,y):
        if "__iter__" not in dir(x):
            return self.wcs_pix2world([[x,y]],1)[0]
        return self.wcs_pix2world([x,y],1)
    
    def wcs2pix(self,ra,dec):
        if "__iter__" not in dir(ra):
            return self.wcs_world2pix([[ra,dec]],0)[0]
        return self.wcs_world2pix([ra,dec],0)

    def coordsAreInImage(self,ra,dec):
        """
        """
        return shape.point_in_contours(ra,dec,self.contours)
        
    
    @property
    def central_coords(self):
        return self.wcs_pix2world([[self._naxis1/2.,self._naxis2/2.]],1)[0]

    @property
    def edge_size(self):
        """return the p1,p2 and p2,p3 size (The order is counter-clockwise starting with the bottom left corner. p1,p2,p3,p4)"""
        p1,p2,p3,p4 = self.calc_footprint()
        return angular_separation(*np.concatenate([p1,p2])),angular_separation(*np.concatenate([p2,p3]))
    
    @property
    def diag_size(self):
        """return the largest diagonal size"""
        p1,p2,p3,p4 = self.calc_footprint()
        return np.max([angular_separation(*np.concatenate([p1,p3])),angular_separation(*np.concatenate([p2,p4]))])
    
    @property
    def contours(self,exp_order=5):
        # -- Load it the fist time you use it
        if "_contour" not in dir(self):
            
            if not shape.HAS_SHAPELY:
                self._contour = None
            else:
                # This could be improved
                npoints = 2+exp_order
                # -- This defines the contours
                width = np.linspace(0,self._naxis1,npoints)
                heigh = np.linspace(0,self._naxis2,npoints)
                v1 = np.asarray([np.ones(npoints-1)*0, width[:-1]]).T
                v2 = np.asarray([heigh[:-1], np.ones(npoints-1)*self._naxis1]).T
                v3 = np.asarray([np.ones(npoints-1)*self._naxis2, width[::-1][:-1]]).T
                v4 = np.asarray([heigh[::-1][:-1], np.ones(npoints-1)*0]).T
                ra,dec = np.asarray([self.pix2wcs(i,j)
                                for i,j in np.concatenate([v1,v2,v3,v4],axis=0)]).T
                
                self._contour = shape.get_contour_polygon(ra,dec)
            
        return self._contour

    def has_contours(self):
        """Test if this have contours defined"""
        return self.contours is not None

"""
class WCS(astWCS.WCS):

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
        return self.contours is not None
"""
