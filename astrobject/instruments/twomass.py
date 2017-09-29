#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from .baseinstrument import Instrument
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools import kwargs_update

__all__  = ["twomass", "TWOMASS_INFO"]

# filter transmission from FSPS allfilters.dat
# effective wavelength from sncosmo
TWOMASS_INFO= {"2massJ":{"lbda":12408.375977},
               "2massH":{"lbda":16513.664599},
               "2massKs":{"lbda":21655.392611},
               "bands":["2massJ","2massH","2massKs"],
               "telescope":{
                   "lon": None,
                   "lat": None}
            }

DATAINDEX = 0

# -------------------- #
# - Instrument Info  - #
# -------------------- #
    
def twomass(*args, **kwargs):
    return TWOMASS(*args, **kwargs)

def is_twomass_file(filename):
    """This test if the input file is a TWOMASS file"""
    return pf.getheader(filename).get("ORIGIN") == "2MASS"

def which_band_is_file(filename):
    """This returns the band of the given file if it is TWOMASS"""
    if not is_twomass_file(filename):
        return None
    return pf.getheader(filename).get("FILTER")

def which_obs_mjd(filename):
    """ read the 2MASS filename and return the
    modified julian date """
    if not is_twomass_file(filename):
        return None
    return get_mjd(pf.getheader(filename))

# -------------------- #
# - Inside tools     - #
# -------------------- #
def get_mjd(tmheader):
    """ read the header of the 2MASS file and return the
    correct modified julian date"""
    from astropy import time
    YYMMDD, UTIME = tmheader["UT_DATE"], tmheader["UT"]
    YY, MM, DD = YYMMDD[0:2], YYMMDD[2:4], YYMMDD[4:6]
    YY = "19"+YY if int(YY)>1 else "20"+YY
    return time.Time('-'.join([YY, MM, DD])+' '+UTIME).mjd
    
#######################################
#                                     #
#   TWOMASS Image Object              #
#                                     #
#######################################
class TWOMASS( Instrument ):
    """
    Class Build to use TWOMASS
    """
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    instrument_name = "TWOMASS"
    INFO            = TWOMASS_INFO
    # ================== #
    #   Properties       #
    # ================== #        
    @property
    def bandname(self):
        """ band of the instrument. Change it using set_bandname() """
        if self._properties['bandname'] is None:
            if self.header is None:
                raise AttributeError("no header loaded")
            self._properties['bandname'] = "2massJ" if self.header['FILTER']=='j' \
                            else "2massH" if self.header['FILTER']=='h' \
                            else "2massKs" 
        return self._properties['bandname']

    @property
    def mjd(self):
        """ Time of the observation in Modified Julian Date """
        if self.header is None:
            raise AttributeError("no header loaded ")
        return get_mjd(self.header)
    
    @property
    def mab0(self):
        """ 
        The ABmag zero point of TWOMASS data. Each image has ZP, 
        but must also convert Vega -> AB (from Brown et al. 2014)
        [0.89, 1.37, 1.84], # Vega->AB, Brown 2014
        """
        imageZP = self.header['MAGZP']
        Vega2AB = 0.89 if self.header['FILTER']=='j' \
                else 1.37 if self.header['FILTER']=='h' \
                else 1.84 
        return imageZP + Vega2AB
    
    # - Low Level image information
    @property
    def _dataunits_to_electron(self):
        """  """
        return None

    @property
    def _gain(self):
        """ The gain of the instrument """
        return None
    
    # =========================== #
    # = Internal Tools          = #
    # =========================== #
    # -----------------
    # - Background hacking

    def _get_default_background_(self,*args,**kwargs):
        return np.zeros(np.shape(self.rawdata))
