#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np

# - astropy
from astropy.io      import fits as pf
from astropy         import time

# - local dependencies
from .baseinstrument    import Instrument
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools      import kwargs_update

GALEX_INFO= {"fuv":{"lbda":1516,"ABmag0":18.82},
             "nuv":{"lbda":2267,"ABmag0":20.08},
             "bands":["fuv","nuv"],
             "telescope":{
                 "lon": None,
                 "lat": None}
            }
    
DATAINDEX = 0


# -------------------- #
# - Instrument Info  - #
# -------------------- #
    
def galex(*args,**kwargs):
    return GALEX(*args,**kwargs)

def is_galex_file(filename):
    """This test if the input file is a GALEX one. Test if 'MPSTYPE' is in the header """
    # not great but this is the structure of MJC images
    return "MPSTYPE" in pf.getheader(filename).keys()

def which_band_is_file(filename):
    """This resuts the band of the given file if it is a
    sdss one"""
    if not is_galex_file(filename):
        return None
    return "fuv" if "-fd-" in filename else "nuv" if "-nd-" in filename \
      else "unknown"

def which_obs_mjd(filename):
    """ read the galex-filename and return the
    modified julian date """
    if not is_galex_file(filename):
        return None
    h_ = pf.getheader(filename)
    return time.Time(h_["OBS-DATE"]+"T"+h_["TIME-OBS"]).mjd


#######################################
#                                     #
#   GALEX Image Object                #
#                                     #
#######################################
class GALEX( Instrument ):
    """
    Class Build to use GALEX 
    """
    PROPERTIES = ["sky"]
    
    instrument_name = "GALEX"

    def set_sky(self, filename=None, skydata=None,
                set_background=True):
        """ Provide the Sky file image or directly the skydata (a Galex Object)
        (in galex filename format they have the skybg label).
        This will set the sky that you can access using self.sky.
        The Variance will be updated except if derive_variance is False
        
        Returns
        -------
        Void
        """
        # - set the sky
        if filename is not None:
            self._properties["sky"] = GALEX(filename, background=0,
                                            dataslice0=self._build_properties["dataslice0"],
                                            dataslice1=self._build_properties["dataslice1"])
            
        elif skydata is not None and GALEX not in skydata.__class__.__mro__:
            raise TypeError("Skydata must be a GALEX instrument file")
        else:
            self._properties["sky"] = skydata

        if self.has_target() and self.has_sky():
            self.sky.set_target(self.target)
            
        # - skybg are is the background of the image:
        if set_background:
            self.set_background(self.sky.rawdata, force_it=True)

    def _derive_variance_(self):
        """ Build the variance image based on the sky+data (=raw int data) assuming pure photon noise. """
        # Pure Photon Noise
        self._properties["var"] = np.sqrt(self.rawdata*self.exposuretime) / self.exposuretime
                
    def _get_default_background_(self,*args,**kwargs):
        return np.zeros(np.shape(self.rawdata))
    
    # ================== #
    #   Properties       #
    # ================== #        
    @property
    def bandname(self):
        """ band of the image. GALEX bands are assimilated to SDSS' ones"""
        return "fuv" if "-fd-" in self.filename else "nuv" if "-nd-" in self.filename \
            else "unknown"

    @property
    def mjd(self):
        """ Time ob the observation in Modified Julian Date """
        return time.Time(self.header["OBS-DATE"]+"T"+self.header["TIME-OBS"]).mjd
    
    @property
    def mab0(self):
        """The ABmag zero point of SDSS data"""
        return GALEX_INFO[self.bandname]["ABmag0"]
    
    # - Low Level image information
    @property
    def _dataunits_to_electron(self):
        """  """
        return None

    # ----------------
    #  GALEX Specific
    @property
    def var(self):
        """ variance image. Poisson noise only in Galex. Based on rawdata of 'int' """
        if self._properties["var"] is None:
            self._derive_variance_()
        return self._properties["var"]
    
    @property
    def sky(self):
        """ """
        return self._properties["sky"]
    
    def has_sky(self):
        return self.sky is not None
    
    @property
    def survey_type(self):
        return self.header["MPSTYPE"]
