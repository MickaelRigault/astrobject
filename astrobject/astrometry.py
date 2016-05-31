#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module manage the Astrometry issue"""

import numpy as np
import warnings
# - local
from .utils import shape

try:
    from astLib      import astWCS
    _HAS_ASTLIB_ = True
except ImportError:
    warnings.warn("ImportError - astLib can not be imported ; PTF-like images won't load.")
    _HAS_ASTLIB_ = False
    
# - astropy
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
        if not _HAS_ASTLIB_:
            raise ImportError(" astLib could not be imported. Astrometry inaccessible for this file."+\
                              " Such file (PTF-like) might be accessible with astLib installed.")
                              
        warnings.warn("WARNING potential issue with wcs SCAM parameters")
        warnings.warn("WARNING astLib based wcs class used")
        try:
            return _WCSbackup(header,mode = "pyfits")
        except:
            from pyfits import getheader
            print "I had to use pyfits for %s"%filename
            return _WCSbackup(getheader(filename,ext=extension),mode = "pyfits")
        
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
        return astlibwcs_to_wcs(wcsinput)
    
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
class _MotherWCS_( object ):
    """ This shoud be merged with WCS once the aslLib dependency is cleaned """
    __nature__ = "WCS"

    # ---------------- #
    # - Converstion  - #
    # ---------------- #
    def units_to_pixels(self,units_, target=None):
        """units should be a parsable string or a astropy units"""
        
        if type(units_) == str:
            # - astropy units
            if units_.lower() in ["pixel","pixels","pxl","pix","pxls","pixs"]:
                return units.Quantity(1.)
            elif units_.lower() in ["kpc"]:
                if target is None:
                    raise AttributeError("You need a target to convert kpc->arcsec")
                units_ = target.arcsec_per_kpc * units.arcsec
            elif units_ in dir(units):
                units_ = 1*units.Unit(units_)
            else:
                TypeError("unparsable units %s: astropy.units or pixels/pxl or kpc could be provided."%units_)
        # ----------------
        # - astropy Units
        if type(units_) is units.core.Unit:
            units_ = 1*units_

        if type(units_) is units.quantity.Quantity:
            return units_.to("degree") / self.pix_indeg
        
        raise TypeError("unparsable units: astropy.units or pixels/pxl or kpc could be provided.")

    # ---------------- #
    # - Offset Stuff - #
    # ---------------- #
    def set_offset(self,offset0,offset1, width=None, height=None):
        self._image_width  = width
        self._image_height = height
        self._reload_contours = True
        self._image_offset = -np.asarray([offset0,offset1])


    @property
    def image_width(self):
        return self._image_width if "_image_width" in dir(self) and \
          self._image_width is not None else self._naxis2
          
    @property
    def image_height(self):
        return self._image_height if "_image_height" in dir(self) and \
          self._image_height is not None else self._naxis1
                    
    @property
    def image_offset(self):
        if "_image_offset" not in dir(self) or self._image_offset is None:
            self.set_offset(0,0,None,None)
        return self._image_offset
    
    @property
    def central_coords(self):
        return self.pix2world(self.image_height/2.,self.image_width/2.)

    @property
    def _central_coords_nooffset(self):
        return self.pix2world(self._naxis1/2.,self._naxis2/2., withoffset=False)
    # ----------------
    # - Contours
    # ----------------
    @property
    def contours(self,exp_order=5):
        # -- Load it the fist time you use it
        if "_contour" not in dir(self) or\
           ("_reload_contours" in dir(self) and self._reload_contours):
            self._reload_contours = False
            from .utils import shape
            if not shape.HAS_SHAPELY:
                self._contour = None
            else:
                a2,a1 = self.image_height,self.image_width
                # This could be improved
                npoints = 2+exp_order
                # -- This defines the contours
                width = np.linspace(0,a1,npoints)
                height = np.linspace(0,a2,npoints)
                v1 = np.asarray([np.ones(npoints-1)*0, width[:-1]]).T
                v2 = np.asarray([height[:-1], np.ones(npoints-1)*a1]).T
                v3 = np.asarray([np.ones(npoints-1)*a2, width[::-1][:-1]]).T
                v4 = np.asarray([height[::-1][:-1], np.ones(npoints-1)*0]).T
                
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

    
class WCS(pyWCS,_MotherWCS_):
    """
    """
    # -------------------- #
    # - Mimic astLib     - #
    # -------------------- #
    def pix2world(self,x,y, withoffset=True):
        """ Convert pixel to world coordinate """
        offset = self.image_offset[::-1] if withoffset else np.asarray([0,0])
        
        if "__iter__" not in dir(x):
            xoffset = x-offset[0]+1
            yoffset = y-offset[1]+1
            return self.wcs_pix2world([[xoffset,yoffset]],
                                        1)[0]
        
        xyoffset = np.asarray([x,y]).T-offset+1    
        return self.wcs_pix2world(xyoffset.tolist(),1)
    
    def world2pix(self,ra,dec,withoffset=True):
        """
        """
        offset = self.image_offset[::-1] if withoffset else np.asarray([0,0])
        if "__iter__" not in dir(ra):
            x,y = self.wcs_world2pix([[ra,dec]],0)[0]
            return x+offset[0],y+offset[1]
            
        x,y = np.asarray(self.wcs_world2pix(np.asarray([ra,dec]).T.tolist(),0)).T
        return np.asarray([np.asarray(x)+offset[0],
                            np.asarray(y)+offset[1]]).T

    def coordsAreInImage(self,ra,dec):
        """
        """
        return shape.point_in_contours(ra,dec,self.contours)
        
    
    
    @property
    def edge_size(self):
        """return the p1,p2 and p2,p3 size (The order is counter-clockwise starting
        with the bottom left corner. p1,p2,p3,p4)"""
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
          
if _HAS_ASTLIB_:
    class _WCSbackup(astWCS.WCS,_MotherWCS_):

        # -------------------- #
        # - Mimic astLib     - #
        # -------------------- #
        def pix2world(self,x,y, withoffset=True):
            """
            """
            offset = self.image_offset if withoffset else np.asarray([0,0])
            x,y = (np.asarray([x,y]).T-offset[::-1]-0.5).T
            return np.asarray(self.pix2wcs(x,y))
    
        def world2pix(self,ra,dec, withoffset=True):
            """
            """
            offset = self.image_offset if withoffset else np.asarray([0,0])
            if "__iter__" not in dir(ra):
                ra,dec = [ra],[dec]
            ra,dec = np.asarray(ra),np.asarray(dec)
            
            x,y = np.asarray(self.wcs2pix(ra,dec)).T

            return zip(np.asarray(x)+offset[1]+0.5,np.asarray(y)+offset[0]+0.5)
        

        def coordsAreInImage(self,ra,dec):
            """
            """
            return shape.point_in_contours(ra,dec,self.contours)
        
        @property
        def edge_size(self):
            """return the p1,p2 and p2,p3 size (The order is counter-clockwise starting with
            the bottom left corner. p1,p2,p3,p4)"""
            p1,p2,p3,p4 = self.calc_footprint()
            return angular_separation(*np.concatenate([p1,p2])),angular_separation(*np.concatenate([p2,p3]))
    
        @property
        def diag_size(self):
            """return the largest diagonal size"""
            p1,p2,p3,p4 = self.calc_footprint()
            return np.max([angular_separation(*np.concatenate([p1,p3])),
                           angular_separation(*np.concatenate([p2,p4]))])

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
        
