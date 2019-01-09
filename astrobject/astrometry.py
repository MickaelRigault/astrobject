#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module manage the Astrometry issue"""

# Moved to ccdphot #

import numpy as np
import warnings

# - astropy
from astropy     import units
from astropy.wcs import WCS as astropyWCS
from astropy.io  import fits
from astropy.coordinates.angle_utilities import angular_separation
    
# - local
from .utils import shape
from .utils.tools import is_arraylike

def get_wcs(wcsinput, verbose=True, **kwargs):
    """
    """
    if type(wcsinput)==str:
        return WCS.load_from_fitsfile(wcsinput, **kwargs)
    elif fits.Header in wcsinput.__class__.__mro__:
        return WCS.load_from_header(wcsinput, **kwargs)
    elif WCS in wcsinput.__class__.__mro__:
        return wcsinput

    raise TypeError("Unrecognized WCS input format (given type: %s) fitsfile / fits.header or WCS implemented"%type(wcsinput))

# = Older
def wcs(filename=None, header=None, extension=0):
    """ loads the WCS solution for the given data. """
    print("DECREPATED, used get_wcs")
    return get_wcs(filename) if filename is not None else get_wcs(header, extension=extension)

##########################################
#                                        #
# Advance WCS Astrobject                 #
#                                        #
##########################################
class WCS( astropyWCS ):
    """ More user friendly methods than the default astropy onces (inherit astropy.wcs.WCS) """
    __nature__ = "WCS"

    # ---------------- #
    #  Statics method  #
    # ---------------- #
    @staticmethod
    def load_from_header(header):
        """ Load the wcs solution from a given header. This is similar to WCS(header)"""
        return WCS(header=header)
    
    @staticmethod
    def load_from_fitsfile(fitsfile, extension=0):
        """ Load the header from the given file and return a WCS object"""
        return WCS.load_from_header( fits.getheader(fitsfile,ext=extension) )
    
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
        self._reload_contours_pxl = True
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
        if not hasattr(self,"_contour") or\
           ( hasattr(self,"_reload_contours") and self._reload_contours ):
            self._reload_contours = False
            from .utils import shape
            if not shape.HAS_SHAPELY:
                self._contour = None
            else:
                a2,a1 = self.image_height,self.image_width
                if a2==0 or a1 == 0:
                    self._contour = None
                else:
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
        if not hasattr(self,"_contour_pxl") or\
           ( hasattr(self,"_reload_contours_pxl") and self._reload_contours_pxl ):
            x,y = np.asarray([self.world2pix(ra_,dec_) for ra_,dec_ in
                            np.asarray(self.contours.exterior.xy).T]).T # switch ra and dec ;  checked
            self._contour_pxl = shape.get_contour_polygon(x,y)
        return self._contour_pxl
    
    def has_contours(self):
        return self.contours is not None
    
    def pix2world(self,x,y, withoffset=True):
        """ Convert pixel to world coordinate """
        offset = self.image_offset[::-1] if withoffset else np.asarray([0,0])
        
        if not is_arraylike(x):
            xoffset = x-offset[0]+1
            yoffset = y-offset[1]+1
            return self.all_pix2world([[xoffset,yoffset]],
                                        0)[0]
        
        xyoffset = np.asarray([x,y]).T-offset+1    
        return self.all_pix2world(xyoffset.tolist(),0)
    
    def world2pix(self,ra,dec,withoffset=True):
        """
        """
        offset = self.image_offset[::-1] if withoffset else np.asarray([0,0])
        if not is_arraylike(ra):
            x,y = self.all_world2pix([[ra,dec]],0)[0]
            return x+offset[0],y+offset[1]
            
        x,y = np.asarray(self.all_world2pix(np.asarray([ra,dec]).T.tolist(),0)).T
        return np.asarray([np.asarray(x)+offset[0],
                            np.asarray(y)+offset[1]]).T

    def coordsAreInImage(self,ra,dec):
        """
        """
        return shape.point_in_contours(ra,dec,self.contours)

    def get_rotations(self, as_rotationmatrix=False, skew_limit=0.1):
        """ returns x and y rotation (in radian)
        if as_rotationmatrix, return
          [[cos(rot), -sin(rot)],[sin(rot), cos(rot)]]
        """
        sign = np.sign(np.linalg.det(self.wcs.cd))
        xrot = np.arctan2(sign * self.wcs.cd[0,1], sign * self.wcs.cd[0,0])
        yrot = np.arctan2(-self.wcs.cd[1,0], self.wcs.cd[1,1])
        if as_rotationmatrix:
            if np.abs(xrot-yrot)>skew_limit/180*np.pi:
                raise NotImplementedError("no rotation matrix implemented for skewed rotation")
            return np.asarray([[np.cos(yrot), -np.sin(yrot)],[np.sin(yrot), np.cos(yrot)]])
        
        return xrot,yrot

    @property
    def rotation_indeg(self, skew_limit=0.1):
        """ rotation of the north with respect of the column axis
        1 value is returned is skewness lower then 0.1 deg (i.e. xrot~yrot)
        else xrot, and yrot are returned
        """
        xrot,yrot = self.get_rotations()*180 / np.pi
        if np.abs(xrot-yrot)>skew_limit:
            return xrot,yrot
        return yrot

    @property
    def rotmatrix(self, skew_limit=0.1):
        """ rotation matrix given the wcs.cd ( see self.get_rotations() ) """
        return self.get_rotations(True, skew_limit=skew_limit)
        
    
    @property
    def edge_size(self):
        """return the p1,p2 and p2,p3 size (The order is counter-clockwise starting
        with the bottom left corner. p1,p2,p3,p4)"""
        lower_left, lower_right, upper_right, upper_left = self.calc_footprint()
        return angular_separation(*np.concatenate([lower_left,lower_right])),angular_separation(*np.concatenate([lower_left,upper_left]))
    
    @property
    def diag_size(self):
        """return the largest diagonal size"""
        lower_left, lower_right, upper_right, upper_left = self.calc_footprint()
        return np.max([angular_separation(*np.concatenate([lower_left,upper_right])),angular_separation(*np.concatenate([lower_right,upper_left]))])

    @property
    def pix_indeg(self):
        """Size in degree. Returns astropy Quantity in degree"""
        [cd1_1,cd1_2],[cd2_1,cd2_2] = self.wcs.cd if hasattr(self.wcs,"cd") else self.pixel_scale_matrix
        pxl = np.sqrt(cd1_1**2+cd2_1**2),np.sqrt(cd1_2**2+cd2_2**2)
    
        if (pxl[0]-pxl[1])/pxl[0] < 1e-2 : # equal to 1%
            return  pxl[0] * units.Unit(self.wcs.cunit[0]).in_units("degree") * units.degree
        
        return  pxl[0] * units.Unit(self.wcs.cunit[0]).in_units("degree")* units.degree,\
          pxl[1] * units.Unit(self.wcs.cunit[1]).in_units("degree")* units.degree
