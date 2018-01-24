#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from astropy import time
from .baseinstrument import Instrument
from ..utils.decorators import _autogen_docstring_inheritance

__all__ = ["snifs","SNIFS_INFO"]

SDSS_INFO = {
             "g":{"lbda":4686,"ABmag0":22.50}, # assuming pure SDSS-SNIFS equality
             "i":{"lbda":7480,"ABmag0":22.50},
             "bands":["g","i"],
             "telescope":{
                 "lon": 19.823,
                 "lat":-155.469}
            }

DATAINDEX = 0


# -------------------- #
# - Instrument Info  - #
# -------------------- #
    
def snifs(*args,**kwargs):
    return SNIFS_P(*args,**kwargs)

def is_snifs_file(filename):
    """This test if the input file is a SNIFS one"""
    # not great but this is the structure of MJC images
    return pf.getheader(filename).get("ORIGIN") == "hyades2.lbl.gov"  

def which_band_is_file(filename):
    """This resuts the band of the given file if it is a
    sdss one"""
    if not is_snifs_file(filename):
        return None
    return filename.split("snifs_")[-1].split(".")[0]

def which_obs_mjd(filename):
    """ read the sdss-filename and return the
    modified julian date """
    if not is_snifs_file(filename):
        return None
    return time.Time(pf.getheader(filename)["DATE"]).mjd



#######################################
#                                     #
#   SNIFS Image Object                #
#                                     #
#######################################
class SNIFS_P( Instrument ):
    """
    Class Build to use SNIFS images processed by M. J. Childress.
    See http://cdsads.u-strasbg.fr/abs/2013ApJ...770..107C.
    """
    instrument_name = "SNIFS_P"

    # ================== #
    #   Properties       #
    # ================== #
    def set_var(self, variance_image):
        """ """
        if np.shape(variance_image) != np.shape(self.data):
            raise ValueError("the given variance image do not have the same shape as the data")
        self._properties["var"] = np.asarray(variance_image)

        
    @property
    def bandname(self):
        """ band of the image. SNIFS bands are assimilated to SDSS' ones"""
        f_ = self.header["FILTER"] if "FILTER" in self.header.keys() \
          else self.filename.split("snifs_")[-1].split(".")[0]
        return "sdss"+f_.lower()

    # - ZeroPoint
    @property
    def mab0(self):
        """The ABmag zero point of SDSS data"""
        if "ZP" in self.header:
            return self.header["ZP"]
        print("WARNING default zp=24 used")
        return 24

    def set_mab0(self, zp, zperr, source="unknonw"):
        """ """
        self.header["ZP"]       = zp
        self.header["ZPERR"]    = zperr
        self.header["ZPSOURCE"] = source

    
    # - Low Level image information
    @property
    def _dataunits_to_electron(self):
        """ Should not be used in the MJC formated SNIFS Photometry """
        return None
    
    @property
    def _gain(self):
        raise self.header["GAIN"]

    
    @property
    def mjd(self):
        """ Time ob the observation in Modified Julian Date """
        return time.Time(self.header["DATE"]).mjd

