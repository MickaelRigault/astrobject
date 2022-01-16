#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from .baseinstrument import Instrument
from ..utils.tools import kwargs_update
from ..photometry import get_photopoint

__all__  = ["wise", "WISE_INFO"]

""" 
WISE images have Vega ZP=22.5. Need to add AB offsets from
http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab
Filter transmissions from FSPS allfilters.dat
Effective wavelength from FITS header
Counts to flux conversion parameters from 
https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_3f.html
"""
WISE_INFO= {"wisew1":{"lbda":33680,"ABmag0":22.5+2.699, "f0":306.682, "s_f0":4.600, "c":1.9350E-06},
            "wisew2":{"lbda":46180,"ABmag0":22.5+3.339, "f0":170.663, "s_f0":2.600, "c":2.7048E-06},
            "wisew3":{"lbda":120820,"ABmag0":22.5+5.174, "f0":29.0448, "s_f0":0.436, "c":1.8326e-06},
            "wisew4":{"lbda":221940,"ABmag0":22.5+6.620, "f0":8.2839, "s_f0":0.124, "c":5.2269E-05},
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
        _var = pf.getdata(filename).byteswap().newbyteorder()
        self._properties["var"] = _var**2
    
    def _aperture_to_photopoint_(self,count,err,flag):
        """ convert the aperture output to a photopoints """
        
        flux, unc = self.count_to_flux(count, err)
        flux, var = np.atleast_1d(flux), np.atleast_1d(unc**2)
        is_single = len(flux)
        
        # ------------------
        # - One Photopoint
        if is_single:
            return get_photopoint(lbda=self.lbda, flux=flux[0], var=var[0],
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
    
    def count_to_flux(self, counts, counts_err=None):
        """ converts counts into flux """
        from ..utils.tools import convert_flux_unit
        _c = WISE_INFO[self.bandname]["c"]
        _f0 = WISE_INFO[self.bandname]["f0"]
        _s_f0 = WISE_INFO[self.bandname]["s_f0"]
        
        _flux = convert_flux_unit(_c * counts, "Jy", "AA", self.lbda)
        if counts_err is not None:
            _flux_unc = _c * np.sqrt(counts**2 * ((_s_f0/_f0)**2 + 0.8483*self.magzp_unc**2) + counts_err**2)
            _flux_unc = convert_flux_unit(_flux_unc, "Jy", "AA", self.lbda)
        return _flux, _flux_unc if counts_err is not None else None

    def flux_to_count(self, flux, flux_err=None):
        """ converts flux into counts """
        from ..utils.tools import convert_flux_unit
        _c = WISE_INFO[self.bandname]["c"]
        _f0 = WISE_INFO[self.bandname]["f0"]
        _s_f0 = WISE_INFO[self.bandname]["s_f0"]
        
        _counts = convert_flux_unit(flux, "AA", "Jy", self.lbda) / _c
        # TO BE DONE
        #if flux_err is not None:
        #    _counts_err = _c * np.sqrt(counts_err**2 * ((_s_f0/_f0)**2 + 0.8483*self.magzp_unc**2) + counts**2)
        #    _flux_unc = convert_flux_unit(_flux_unc, "Jy", "AA", self.lbda)
        return _flux, None#_flux_unc if counts_err is not None else None
    
    
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
    
    @property
    def magzp(self):
        """ Relative photometric zero point """
        return self.header['MAGZP']
    
    @property
    def magzp_unc(self):
        """ 1-sigma uncertainty in zero point """
        return self.header['MAGZPUNC']
    
    
    # - Low Level image information
    @property
    def _dataunits_to_electron(self):
        """  """
        return None

    @property
    def _gain(self):
        """ The gain of the instrument """
        return None

    
