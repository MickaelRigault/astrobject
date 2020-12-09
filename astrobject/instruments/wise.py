#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from .baseinstrument import Instrument
from ..utils.tools import kwargs_update

__all__  = ["wise", "WISE_INFO"]

""" 
WISE images have Vega ZP=22.5. Need to add AB offsets from
http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab
Filter transmissions from FSPS allfilters.dat
Effective wavelength from sncosmo
"""
WISE_INFO= {"wisew1":{"lbda":33791.912981,"ABmag0":22.5+2.699},
            "wisew2":{"lbda":46292.961143,"ABmag0":22.5+3.339},
            "wisew3":{"lbda":123321.610431,"ABmag0":22.5+5.174},
            "wisew4":{"lbda":222532.734640,"ABmag0":22.5+6.620},
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
    """This test if the input file is a WISE file"""
    return "WISE" in pf.getheader(filename).get("ORIGIN")

def which_band_is_file(filename):
    """This returns the band of the given file if it is WISE"""
    if not is_wise_file(filename):
        return None
    return f"w{pf.getheader(filename).get('BAND')}"

def which_obs_mjd(filename):
    """WISE FITS files are coadds so no MJD info"""
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
    
    def __init__(self, filename=None, uncfilename=None, 
                 astrotarget=None,data_index=0,
                 dataslice0=None,dataslice1=None,
                 empty=False, **kwargs):
        """
        Initalize the image by giving its filelocation (*filename*). This
        will load it using the load() method.

        Parameters
        ----------
        filename: [string.fits]    fits file from where the image will be loaded
                                   - Trick - Set None and no image will be loaded
        
        uncfilename: [string.fits] fits file from where the error map will be loaded
                                   - Trick - Set None and no image will be loaded

        astrotarget: [AstroTarget] An AstroTarget object you which to associate
                                   to this image. 
                                   
        empty: [bool]              Does not do anything, just loads an empty object.
                                   (Careful with that)

        Returns
        -------
        Void
        """
        super(WISE,self).__init__(filename=filename, astrotarget=astrotarget, data_index=data_index, 
                                  dataslice0=dataslice0, dataslice1=dataslice1, empty=empty, **kwargs)
        if uncfilename:
            self.set_var(uncfilename)

    def set_var(self, filename):
        """ Get the error map from separate file 
        which should have the same file name as the image file 
        except with 'int' replaced with 'unc'. """
        self._properties["var"] = pf.getdata(filename).byteswap().newbyteorder()
        
    
    
    # ================== #
    #   Properties       #
    # ================== #        
    @property
    def bandname(self):
        """ band of the instrument. Change it using set_bandname() """
        if self._properties['bandname'] is None:
            if self.header is None:
                raise AttributeError("no header loaded")
            self._properties['bandname'] = f"wisew{self.header['BAND']}"
        return self._properties['bandname']

    @property
    def mjd(self):
        """ Time of the observation in Modified Julian Date """
        # unWISE FITS files are coadds so no MJD info
        return None
    
    @property
    def mab0(self):
        """ The ABmag zero point of WISE data """
        return WISE_INFO[self.bandname]["ABmag0"]
    
    # - Low Level image information
    @property
    def _dataunits_to_electron(self):
        """  """
        return None

    @property
    def _gain(self):
        """ The gain of the instrument """
        return None

    
