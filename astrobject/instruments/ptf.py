#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from astropy import units
from .baseinstrument import Instrument
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools import kwargs_update

__all__ = ["ptf","PTF_INFO"]

PTF_INFO = {"telescope":{
                 "lon": None,
                 "lat": None}
            }

DATAINDEX = 0


# -------------------- #
# - Instrument Info  - #
# -------------------- #


def ptf(*args,**kwargs):
    return PTF(*args,**kwargs)

def is_ptf_file(filename):
    """This test if the input file is a ptf one"""
    # TO  BE CHANGED
    return True if pf.getheader(filename).get("ORIGIN") == "Palomar Transient Factory" \
      else False

def which_band_is_file(filename):
    """This resuts the band of the given file if it is a
    ptf one"""
    if not is_ptf_file(filename):
        return None
    return pf.getheader(filename).get("FILTER")



#######################################
#                                     #
#   PTF Image Object                  #
#                                     #
#######################################
class PTF( Instrument ):
    """
    This is the image object custom for PTF data
    """
    instrument_name = "PTF"

    PROPERTIES         = ["mab0"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    def __build__(self,**kwargs):
        """
        """
        # -- Load the basic builds
        super(PTF,self).__build__()
        # -- How to read the image
        self._build_properties = kwargs_update(self._build_properties,
                                               **dict(
                    data_index = DATAINDEX,
                    header_exptime = "EXPTIME"
                    ))

    @_autogen_docstring_inheritance(Instrument.set_catalogue,"Instrument.set_catalogue")
    def set_catalogue(self,catalogue,force_it=True,**kwargs):
        #
        # - Add the bandname key_mag setting
        #
        if catalogue.source_name =="SDSS":
            key_mag = "%smag"%self.bandname[0].lower()
            key_magerr = "e_%smag"%self.bandname[0].lower()
            if key_mag not in catalogue.data.keys():
                warnings.warn("WARNING No %s in the catalogue data. Cannot assign a key_mag"%key_mag)
            catalogue.set_mag_keys(key_mag,key_magerr)

        super(PTF,self).set_catalogue(catalogue,force_it=force_it,**kwargs)
    

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    # --------------
    # - Image Data
    @property
    def data(self):
        # raw data in PTF are in count
        return self._derived_properties["data"] #/ self.exposuretime
    
    @property
    def bandname(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return "sdss%s"%(self.header["FILTER"][0])

    # --------------
    # - ZERO POINTS
    @property
    def mab0(self):
        if self._properties["mab0"] is None:
            self.mab0 = 20#self.header["ZP"]
            
        return self._properties["mab0"]

    @mab0.setter
    def mab0(self,value):
        if value<10:
            raise ValueError("given mab0 lozer then 10 looks unlikely")
        self._properties["mab0"] = value

    # --------------
    # - FWHM
    @property
    def fwhm(self,fromheader=True):
        if not self.has_fwhm():
            if fromheader:
                self.set_fwhm(self.header["SEEING"] * units.arcsec)
            elif not self.has_sepobjects():
                raise AttributeError("'fwhm' is not defined and no sepobjects loaded."+\
                                     " You could use the header value (set fromheader to True)")
            fwhm_pxl = self.sepobjects.get_fwhm_pxl(isolated_only=True,stars_only=True)
            self.set_fwhm(fwhm_pxl/self.units_to_pixels("arcsec").value*units.arcsec)
            
        return self._derived_properties["fwhm"]
    # --------------
    # - Times
    @property
    def mjd(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.header["OBSMJD"]

    # --------------
    # - Low Level
    @property
    def ccd_id(self):
        return self.header["CCDID"]
    @property
    def _gain(self):
        " In e-/adu"
        return float(self.header["GAIN"])
    @property
    def _readout_noise(self):
        """read noise e-"""
        return float(self.header["READNOI"])
        
    # =========================== #
    # = Internal Tools          = #
    # =========================== #
