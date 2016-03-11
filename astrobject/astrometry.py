#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module manage the Astrometry issue"""

import numpy as np
from ..utils import shape
import warnings
from astLib      import astWCS
from astropy.wcs import WCS as pyWCS
from astropy.io  import fits
from astropy.coordinates.angle_utilities import angular_separation
from astropy import units

def wcs(filename=None,header=None,extension=0):
    """
    """
    if filename is not None:
        header = fits.getheader(filename,ext=extension)
    if header is None:
        raise ValueError("'filename' or 'header' must be given")
    
    if "PV1_5" in header.keys():
        warnings.warn("WARNING potential issue with wcs SCAM parameters")
        warnings.warn("WARNING astLib based wcs class used")
        return _WCSbackup(header,mode = "pyfits")
    
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
    def pix2world(self,x,y):
        if "__iter__" not in dir(x):
            return self.wcs_pix2world([[x,y]],1)[0]
        
        return self.wcs_pix2world(np.asarray([x,y]).T.tolist(),1)
    
    def world2pix(self,ra,dec):
        if "__iter__" not in dir(ra):
            return self.wcs_world2pix([[ra,dec]],0)[0]
         
        return self.wcs_world2pix(np.asarray([ra,dec]).T.tolist(),0)

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
    def pix_indeg(self):
        """Size in degree. Returns astropy Quantity in degree"""
        [cd1_1,cd1_2],[cd2_1,cd2_2] = self.wcs.cd
        pxl = np.sqrt(cd1_1**2+cd2_1**2),np.sqrt(cd1_2**2+cd2_2**2)
    
        if (pxl[0]-pxl[1])/pxl[0] < 1e-2 : # equal to 1%
            return  pxl[0] * units.Unit(self.wcs.cunit[0]).in_units("degree") * units.degree
        
        return  pxl[0] * units.Unit(self.wcs.cunit[0]).in_units("degree")* units.degree,\
          pxl[1] * units.Unit(self.wcs.cunit[1]).in_units("degree")* units.degree
          
    
    @property
    def contours(self,exp_order=5):
        # -- Load it the fist time you use it
        if "_contour" not in dir(self):
            from ..utils import shape
            if not shape.HAS_SHAPELY:
                self._contour = None
            else:
                a2,a1 = self._naxis1,self._naxis2
                # This could be improved
                npoints = 2+exp_order
                # -- This defines the contours
                width = np.linspace(0,a1,npoints)
                heigh = np.linspace(0,a2,npoints)
                v1 = np.asarray([np.ones(npoints-1)*0, width[:-1]]).T
                v2 = np.asarray([heigh[:-1], np.ones(npoints-1)*a1]).T
                v3 = np.asarray([np.ones(npoints-1)*a2, width[::-1][:-1]]).T
                v4 = np.asarray([heigh[::-1][:-1], np.ones(npoints-1)*0]).T
                
                ra,dec = np.asarray([self.pix2world(i,j)
                                for i,j in np.concatenate([v1,v2,v3,v4],axis=0)]).T
                self._contour = shape.get_contour_polygon(ra,dec)
            
        return self._contour

    @property
    def contours_pxl(self,**kwargs):
        """Based on the contours (in wcs) and wcs2pxl, this draws the pixels contours"""
        x,y = np.asarray([self.world2pix(ra_,dec_) for ra_,dec_ in
                          np.asarray(self.contours.exterior.xy).T]).T # switch ra and dec ;  checked
        return shape.get_contour_polygon(x,y)
    
    def has_contours(self):
        """Test if this have contours defined"""
        return self.contours is not None


class _WCSbackup(astWCS.WCS ):

    __nature__ = "WCS"

    # -------------------- #
    # - Mimic astLib     - #
    # -------------------- #
    def pix2world(self,x,y):
        """
        """
        if "__iter__" in dir(x):
            x,y = np.asarray(x),np.asarray(y)
        return self.pix2wcs(x,y)
        
    def world2pix(self,ra,dec):
        """
        """
        if "__iter__" in dir(ra):
            ra,dec = np.asarray(ra),np.asarray(dec)
        return self.wcs2pix(ra,dec)
        

    def coordsAreInImage(self,ra,dec):
        """
        """
        return shape.point_in_contours(ra,dec,self.contours)
        
    @property
    def central_coords(self):
        return self.pix2world(self._naxis1/2.,self._naxis2/2.)

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

    # ===================== #
    # = 
    # ===================== #
    def calc_footprint(self,center=True):
        naxis1,naxis2 = self._naxis1,self._naxis2
        if center == True:
            corners = np.array([[1, 1],
                                [1, naxis2],
                                [naxis1, naxis2],
                                [naxis1, 1]], dtype = np.float64)
        else:
            corners = np.array([[0.5, 0.5],
                                [0.5, naxis2 + 0.5],
                                [naxis1 + 0.5, naxis2 + 0.5],
                                [naxis1 + 0.5, 0.5]], dtype = np.float64)
        return np.asarray([self.pix2world(*p_) for p_ in corners])

        
    @property
    def _naxis1(self):
        return self.header["NAXIS1"]
    
    @property
    def _naxis2(self):
        return self.header["NAXIS2"]
    
    @property
    def pix_indeg(self):
        return  self.getPixelSizeDeg() * units.degree
        
          
    # ===================== #
    # = 
    # ===================== #
    
    @property
    def contours(self,exp_order=5):
        # -- Load it the fist time you use it
        if "_contour" not in dir(self):
            from ..utils import shape
            if not shape.HAS_SHAPELY:
                self._contour = None
            else:
                a2,a1 = self._naxis1,self._naxis2
                # This could be improved
                npoints = 2+exp_order
                # -- This defines the contours
                width = np.linspace(0,a1,npoints)
                heigh = np.linspace(0,a2,npoints)
                v1 = np.asarray([np.ones(npoints-1)*0, width[:-1]]).T
                v2 = np.asarray([heigh[:-1], np.ones(npoints-1)*a1]).T
                v3 = np.asarray([np.ones(npoints-1)*a2, width[::-1][:-1]]).T
                v4 = np.asarray([heigh[::-1][:-1], np.ones(npoints-1)*0]).T
                
                ra,dec = np.asarray([self.pix2world(i,j)
                                for i,j in np.concatenate([v1,v2,v3,v4],axis=0)]).T
                
                self._contour = shape.get_contour_polygon(ra,dec)
            
        return self._contour

    @property
    def contours_pxl(self,**kwargs):
        """Based on the contours (in wcs) and wcs2pxl, this draws the pixels contours"""
        x,y = np.asarray([self.world2pix(ra_,dec_) for ra_,dec_ in
                          np.asarray(self.contours.exterior.xy).T]).T # switch ra and dec ;  checked
        return shape.get_contour_polygon(x,y)
    
    def has_contours(self):
        return self.contours is not None
