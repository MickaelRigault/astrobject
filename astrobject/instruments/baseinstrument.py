#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np

# - astropy
from astropy.io import fits as pf
from astropy.io import ascii

# - local
from ..photometry import Image, get_photopoint
from ..baseobject import BaseObject, WCSHandler
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools import kwargs_update, mag_to_flux, load_pkl, dump_pkl, is_arraylike

__all__ = ["Instrument"]

class Instrument( Image ):
    """
    """
    PROPERTIES         = ["bandname"]
    DERIVED_PROPERTIES = ["bandpass"]
    
    instrument_name = "TO_BE_DEFINED"
    INFO            = {} # set one
    
    def __build__(self,data_index=0):
        """This is a slightly advanced Image object"""
        super(Instrument,self).__build__(data_index=data_index)

    # ---------------- #
    #  PhotoPoints     #
    # ---------------- #
    
    def _aperture_to_photopoint_(self,count,err,flag):
        """ convert the aperture output to a photopoints """
        
        flux = self.count_to_flux(count)
        var  = self.count_to_flux(err)**2
        
        # ------------------
        # - One Photopoint
        if not is_arraylike(flux):
            return get_photopoint(lbda=self.lbda, flux=flux, var=var,
                            source="image",mjd=self.mjd,
                            zp=self.mab0,bandname=self.bandname,
                            instrument_name=self.instrument_name)
        # -----------------------
        # - Several Photopoints
        return [get_photopoint(lbda=self.lbda,flux=flux_,var=var_,
                            source="image",mjd=self.mjd,
                            zp=self.mab0,bandname=self.bandpass.name,
                            instrument_name=self.instrument_name)
                            for flux_,var_ in zip(flux,var)]
    
        
    @_autogen_docstring_inheritance(Image.get_aperture,"Image.get_aperture")
    def get_photopoint(self,x, y, radius=None, runits="pixels",
                       aptype="circle", wcs_coords=False, getlist=False,
                       **kwargs):
        #
        # Returns a PhotoPoint
        #
        pp = self._aperture_to_photopoint_(*self.get_aperture(x,y,radius=radius,runits=runits,
                                                            aptype=aptype,wcs_coords=wcs_coords,
                                                            **kwargs))
        # ------------------
        # - One Photopoint
        if not is_arraylike(pp) or len(pp)==1 or getlist:
            return pp
        
        # -----------------------
        # - Several Photopoints\
        from ..collections import get_photomap
        return get_photomap(pp,np.asarray([x,y]).T,
                        wcs=self.wcs,
                        catalogue=self.catalogue.get_subcatalogue(fovmask=True, catmag_range=[1,30]) \
                          if self.has_catalogue() else None,
                        wcs_coords=wcs_coords)
    
    @_autogen_docstring_inheritance(Image.get_target_aperture,
                                    "Image.get_target_aperture")
    def get_target_photopoint(self,radius=None, runits="pixels",
                              aptype="circle", **kwargs):
        #
        # Returns a PhotoPoint
        #
        if not self.has_target():
            return AttributeError("No 'target' loaded")
        
        xpix,ypix = self.coords_to_pixel(self.target.ra,self.target.dec)
        pp = self.get_photopoint(xpix,ypix,radius=radius,runits=runits,
                                   aptype="circle",**kwargs)
        pp.set_target(self.target)
        return pp
    
    @_autogen_docstring_inheritance(Image.get_host_aperture,
                                    "Image.get_host_aperture")
    def get_host_photopoint(self,scaleup=2.5,**kwargs):
        #
        # Be returns a PhotoPoint
        #
        if not self.has_target():
            return AttributeError("No 'target' loaded")
        ap_output = self.get_host_aperture(scaleup=scaleup,**kwargs)
        if ap_output is None:
            return None
        
        pp = self._aperture_to_photopoint_(*ap_output)
        pp.set_target(self.target)
        return pp

    
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def count_to_flux(self,counts):
        """ converts counts into flux """
        return counts * 10**(-(2.406+self.mab0) / 2.5 ) / (self.lbda**2)

    # =========================== #
    # = Internal Catalogue      = #
    # =========================== #
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    
    # ------------------
    # - Band Information
    @property
    def bandname(self):
        if self._properties['bandname'] is None:
            raise NotImplementedError("'band' must be implemented")
        return self._properties['bandname']
    
    def set_bandname(self, value, reset_bandpass=True):
        """ Change the name of the bandname.
        
        """
        if value is not None:
            if type(value) != str and type(value) != np.string_:
                raise TypeError("The bandname must be a string", type(value))
        
        self._properties["bandname"] = value
        if reset_bandpass:
            self._derived_properties["bandpass"] = None
        
    @property
    def bandpass(self):
        """ Object containing the basic information associated to the bandname """
        # - No bandname ?
        if self._derived_properties["bandpass"] is not None:
            return self._derived_properties["bandpass"]
        
        if self.bandname is None:
            raise AttributeError("No bandname given")
        
        # - Should this use sncosmo
        try:
            from sncosmo import get_bandpass
            has_sncosmo = True
        except ImportError:
            warnings.warn("sncosmo is not installed. Could not access the bandpass")
            has_sncosmo = False
            
        if has_sncosmo:
            try:
                self._derived_properties["bandpass"] = get_bandpass(self.bandname)
                use_default_bandpass = False
            except:
                use_default_bandpass = True
        else:
            use_default_bandpass = True
            
        if use_default_bandpass:
            wave_eff = self.INFO[self.bandname]["lbda"] \
              if self.bandname in self.INFO else np.NaN
            self._derived_properties["bandpass"] = _DefaultBandpass_(self.bandname, wave_eff)
            
        return self._derived_properties["bandpass"]
    
    @property
    def lbda(self):
        """ effective wavelength """
        return self.bandpass.wave_eff

    # -- Derived values
    @property
    def mab0(self):
        raise NotImplementedError("'mab0' must be implemented")

    @property
    def _gain(self):
        raise NotImplementedError("'_gain' must be implemented (even for None)")

    @property
    def _dataunits_to_electron(self):
        """The gain converts ADU->electron. The Data shouls be in ADU/s"""
        return self._gain * self.exposuretime
    
    @property
    def mjd(self):
        raise NotImplementedError("'mjd' (Modified Julian Date) must be implemented")

############################################
#                                          #
# Base Instrument: CATALOGUE               #
#                                          #
############################################
#from ..catalogue.basecatalogue import Catalogue

    
############################################
#                                          #
# Backup Plan if no SNCOSMO                #
#                                          #
############################################
    
class _DefaultBandpass_( BaseObject ):
    """ """
    PROPERTIES = ["bandname", "wave_eff"]
    def __init__(self, bandname, wave_eff):
        """ """
        self._properties["bandname"] = bandname
        self._properties["wave_eff"] = wave_eff

    @property
    def name(self):
        return self._properties["bandname"]
    @property
    def wave_eff(self):
        return self._properties["wave_eff"]
