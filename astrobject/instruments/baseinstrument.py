#! /usr/bin/env python
# -*- coding: utf-8 -*-

from ..utils.decorators import _autogen_docstring_inheritance
from ..astrobject.photometry import photopoint,Image

__all__ = ["Instrument"]



class Instrument( Image ):
    """
    """
    def __build__(self,data_index=0):
        """This is a slightly advanced Image object"""
        super(Instrument,self).__build__(data_index=data_index)

    # ----------- #
    #  PhotoPoint #
    # ----------- #
    @_autogen_docstring_inheritance(Image.get_aperture,"Image.get_aperture")
    def get_photopoint(self,x,y,r_pixels=None,
                       aptype="circle",
                       **kwargs):
        #
        # Be returns a PhotoPoint
        #
        count,err,flag  = self.get_aperture(x,y,r_pixels=r_pixels,aptype=aptype,
                                           **kwargs)
        flux = self.count_to_flux(count)
        var  = self.count_to_flux(err)**2
        return photopoint(self.lbda,flux,var,source="image",
                          instrument_name=self.instrument_name)
    
    @_autogen_docstring_inheritance(Image.get_target_aperture,"Image.get_target_aperture")
    def get_target_photopoint(self,r_pixels=None,
                              aptype="circle",**kwargs):
        #
        # Be returns a PhotoPoint
        #
        if not self.has_target():
            return AttributeError("No 'target' loaded")
        
        xpix,ypix = self.coords_to_pixel(self.target.ra,self.target.dec)
        return self.get_photopoint(xpix,ypix,r_pixels=r_pixels,
                                   aptype="circle",**kwargs)
    
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def count_to_flux(self,counts):
        return counts* 10**(-(2.406+self.mab0) / 2.5 ) / (self.lbda**2)
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def band(self):
        raise NotImplementedError("'band' must be implemented")

    @property
    def band_info(self):
        raise NotImplementedError("'band_info' must be implemented")

    # -- Derived values
    @property
    def lbda(self):
        return self.band_info["lbda"]

    @property
    def mjd_obstime(self):
        raise NotImplementedError("'obstime' must be implemented")
    
    @property
    def mab0(self):
        return self.band_info["ABmag0"]

    @property
    def _gain(self):
        raise NotImplementedError("'_gain' must be implemented (even for None)")
        
    
    
    


