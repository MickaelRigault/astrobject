#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from .baseinstrument import Instrument
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools import kwargs_update

__all__  = ["wise", "WISE_INFO"]

WISE_INFO= {"wisew1":{"lbda":33680.0,"ABmag0":22.5},
            "wisew2":{"lbda":46180.0,"ABmag0":22.5},
            "wisew3":{"lbda":120820.0,"ABmag0":22.5},
            "wisew4":{"lbda":221940.0,"ABmag0":22.5},
            "bands":["wisew1","wisew2","wisew3","wisew4"],
            "telescope":{
                 "lon": None,
                 "lat": None}
            }

DATAINDEX = 0

# -------------------- #
# - Instrument Info  - #
# -------------------- #
    
def wise(*args, **kwargs):
    return WISE(*args, **kwargs)

def is_wise_file(filename):
    """This test if the input file is an unWISE file"""
    # No useful info in the header! So must rely on filename itself
    return "unwise" in filename

def which_band_is_file(filename):
    """This returns the band of the given file if it is WISE"""
    if not is_wise_file(filename):
        return None
    if   "-w1-" in filename: return "wisew1"
    elif "-w2-" in filename: return "wisew2"
    elif "-w3-" in filename: return "wisew3"
    elif "-w4-" in filename: return "wisew4"
    else: return "unknown"

def which_obs_mjd(filename):
    """unWISE FITS files are coadds so no MJD info"""
    return None

#######################################
#                                     #
#   WISE Image Object                 #
#                                     #
#######################################
class WISE( Instrument ):
    """
    Class Build to use WISE
    """
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    instrument_name = "WISE"
    INFO            = WISE_INFO
    # ================== #
    #   Properties       #
    # ================== #        
    @property
    def bandname(self):
        """ band of the instrument. Change it using set_bandname() """
        if self._properties['bandname'] is None:
            self._properties['bandname'] = "wisew1" if "-w1-" in self.filename \
                            else "wisew2" if "-w2-" in self.filename \
                            else "wisew3" if "-w3-" in self.filename \
                            else "wisew4" if "-w4-" in self.filename \
                            else "unknown"
        return self._properties['bandname']

    @property
    def mjd(self):
        """ Time of the observation in Modified Julian Date """
        # unWISE FITS files are coadds so no MJD info
        return None
    
    @property
    def mab0(self):
        """ The ABmag zero point of WISE data """
        return GALEX_INFO[self.bandname]["ABmag0"]
    
    # - Low Level image information
    @property
    def _dataunits_to_electron(self):
        """  """
        return None

    @property
    def _gain(self):
        """ The gain of the instrument """
        return None

    def set_invvar(self, filename):
	""" Get the inverse-variance map from separate file 
        which should have the same file name as the image file 
        except with '-img-' replace with '-invar-'. """
	self._properties["var"] = 1.0/pf.getdata(filename)

    
