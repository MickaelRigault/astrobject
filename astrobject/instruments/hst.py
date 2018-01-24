#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from astropy.table import Table
from .baseinstrument import Instrument
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools import kwargs_update

__all__ = ["hst"]


_FWHM_WFC3 = Table(data=[[2000,3000,4000,5000,6000,7000,8000,9000,10000,11000],
                         [0.083,0.075,0.070,0.067,0.067,0.070,0.074,0.078,0.084,0.089]],
                   names=["wavelength","FWHM"])

# F225W filter from David Rubin: Averaged chip 1 and chip 2
# transmission from SYNPHOT (STSci software)
HST_INFO = {
    "telescope":{
                 "lon":np.NaN,
                 "lat":np.NaN},
    # Windhorst et al 2010
    "f225w":{"lbda":2372.069991, "fwhm":0.092,
             "ABmag0":24.06, # abmag @ e-/s
             "skybkgd":25.46 # mag per arcsec 2
             },
    "f275w":{"lbda":2708.529114, "fwhm":0.087,
             "ABmag0":24.14,
             "skybkgd":25.64 # mag per arcsec 2
            }
        
             
    }

DATAINDEX = 1
_NAN_OUTSIDE_DATA = True
# -------------------- #
# - Instrument Info  - #
# -------------------- #  
def hst(*args,**kwargs):
    return HST(*args,**kwargs)

def is_hst_file(filename):
    """This test if the given file is an HST one"""
    return pf.getheader(filename).get("TELESCOP") == "HST" 

def which_band_is_file(filename):
    """
    """
    if not is_hst_file(filename):
        return None
    return pf.getheader(filename).get("FILTER")


def get_psf(wavelength_angstrom,show=True):
    """This function enables to get the fwhm in arcsec
    for the given wavelength. This works for the WFC3.
    http://documents.stsci.edu/hst/wfc3/documents/handbooks/cycle21/c06_uvis07.html
    """
    from scipy.interpolate import UnivariateSpline
    interpolate = UnivariateSpline(_FWHM_WFC3["wavelength"],_FWHM_WFC3["FWHM"],k=4,s=0)
    returned = interpolate(wavelength_angstrom)
    # --------- #
    # - Show  - #
    # --------- #
    if show:
        _lbda = np.linspace(_FWHM_WFC3["wavelength"][0],_FWHM_WFC3["wavelength"][-1],100)
        import matplotlib.pyplot as mpl
        # ------------
        # - Ax Plot
        fig = mpl.figure(figsize=[10,8])
        ax  = fig.add_axes([0.12,0.12,0.8,0.8])
        ax.set_xlabel(r"$\mathrm{Wavelenth\ in\ \AA}$",fontsize="x-large")
        ax.set_ylabel(r"$\mathrm{FWHM\ in\ arcsec}$",fontsize="x-large")
        # ------------
        # - Da Plot
        ax.plot(_FWHM_WFC3["wavelength"],_FWHM_WFC3["FWHM"],ls="None",marker="o",
                ms=10, label=r"$\mathrm{STScI\ data}$",
                color=mpl.cm.Blues(0.8),mew=0,)
        color_interpol = mpl.cm.Greens(0.8)
        # - Interpolate
        ax.plot(_lbda,interpolate(_lbda),ls="-",color=color_interpol,
                label=r"$\mathrm{Interpolated\ values}$")
        # - Estimated Value
        ax.plot(wavelength_angstrom,returned,marker="s",ms=10,
                mfc=color_interpol,mew=1,alpha=0.8,ls="None",
                label=r"$\mathrm{Estimated\ value}$")
        ax.axvline(wavelength_angstrom,ls="--",
                   color=color_interpol,alpha=0.5)
        ax.axhline(returned,ls="--",
                   color=color_interpol,alpha=0.5)
        # -- Legend
        ax.legend(loc="best", frameon=False, fontsize="large",numpoints=1)
        fig.show()
        
    return returned
        
########################################
#                                      #
# HST INSTRUMENT CLASS                 #
#                                      #
########################################    
class HST( Instrument ):
    """This is the umage object custom for HST data"""
    
    instrument_name = "HST"
    INFO            = HST_INFO
    
    PROPERTIES         = ["used_amplifier"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    def __build__(self,**kargs):
        """  """
        super(HST,self).__build__()
        # -- How to read the image
        self._build_properties = kwargs_update(self._build_properties,
                                               **dict(
                    data_index = DATAINDEX,
                    error_index = DATAINDEX,
                    header_exptime = "EXPTIME"
                    ))

    # ------------------------ #
    # - Speciality           - #
    # ------------------------ #
    @_autogen_docstring_inheritance(Instrument.set_catalogue,"Instrument.set_catalogue")
    def set_catalogue(self,catalogue,force_it=True,**kwargs):
        #
        # - Add the bandname key_mag setting
        #
        if catalogue.source_name =="SDSS":
            key_mag = "%smag"%"r"
            key_magerr = "e_%smag"%"r"
            if key_mag not in catalogue.data.keys():
                warnings.warn("WARNING No %s in the catalogue data. Cannot assign a key_mag"%key_mag)
            catalogue.set_mag_keys(key_mag,key_magerr)
            
        super(HST,self).set_catalogue(catalogue,force_it=force_it,**kwargs)
        
    def get_contours(self,pixel=True):
        """ """
        from ..utils import shape
        if pixel:
            x = [2060,9,-1,self.width,2060,19,9,2060]
            y = [1154,1027,-1,114,self.height,2056,1038,1166]
            #x,y= [-1,19,self.width,2060],[-1,2056,114,self.height]
            return shape.polygon.Polygon(zip(x,y))

        
        x,y = np.asarray([self.pixel_to_coords(x_,y_) for x_,y_ in
                          np.asarray(self.get_contours(pixel=True).exterior.xy).T]).T # switch ra and dec ;  checked
        return shape.get_contour_polygon(x,y)

        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #    
    # --------------
    # - Image Data
    @property
    def exposuretime(self):
        if self._side_properties['exptime'] is None:
            # -- It has not be set manually, maybe check the header
            self._side_properties['exptime'] = \
              np.float(self.fits[0].header[self._build_properties["header_exptime"]])
              
        # -- You have it ? This will stay None if not
        return self._side_properties['exptime']
    
    @property
    def bandname(self):
        """ Bandname of the instrument. Change it using set_bandname() """
        if self._properties['bandname'] is None:
            if self.header is None:
                raise AttributeError("no header loaded ")
            self._properties['bandname'] = self.fits[0].header["FILTER"]
        return self._properties['bandname']

    # --------------------
    # - Band Information
    @property
    def mjd(self):
        """This is the Modify Julien Date at the start of the Exposure"""
        return self.fits[0].header["EXPSTART"]

    @property
    def mab0(self):
        """zeropoint in ABmag
        http://www.stsci.edu/hst/wfc3/phot_zp_lbn
        ABMAG_ZEROPOINT = -2.5 Log (PHOTFLAM) - 21.10 - 5 Log (PHOTPLAM) + 18.6921 
        """
        return -2.5*np.log10(self.header["PHOTFLAM"]) - 21.10 \
          - 5*np.log10(self.header["PHOTPLAM"]) + 18.6921
          
    # -------------------
    # - Instrument Info
    @property
    def fwhm(self):
        """
        The FWHM in arcsec is based on WFC3 handbook.
        See the module's get_psf() 
        """
        if self._derived_properties["fwhm"] is None:
            self.set_fwhm(get_psf(self.bandpass.wave_eff, show=False) )
            
        return self._derived_properties["fwhm"]

    @property
    def used_amplifier(self):
        """This is the amplifier used for the analysis"""
        if self._properties["used_amplifier"] is None:
            warnings.warn("Default 'C' amplifier set for the hst image")
            self._properties["used_amplifier"] = "C"
            
        return self._properties["used_amplifier"]

    @used_amplifier.setter
    def used_amplifier(self,value):
        """Change the value of the amplifier [A,B,C, or D]"""
        if value.upper() not in ["A","B","C","D"]:
            raise ValueError("the amplifier must be A,B,C, or D")
        
        self._properties["used_amplifier"] = value.upper()

    @property
    def _flagdata(self):
        """Based on PIL ; this is the mask following the pixel contours.
        Most likely, this is by far the fastest technique to do so
        """
        from PIL import Image, ImageDraw
        back = Image.new('RGBA', (self.width, self.height), (100,0,0,0))
        mask = Image.new('RGBA', (self.width, self.height))
        # PIL *needs* (!!) [(),()] format [[],[]] won
        tuplemask = [(x_[0],x_[1]) for x_ in np.asarray(self.get_contours(pixel=True).exterior.xy).T]
        ImageDraw.Draw(mask).polygon( tuplemask,fill=(255,255,255,127),outline=(255,255,255,255))
        back.paste(mask,mask=mask)
        invalue = 494 # because of the color chosen above
        return (np.sum(np.array(back),axis=2)==invalue)
    
    @property
    def datamask(self):
        """ """
        if self._side_properties["datamask"] is None and _NAN_OUTSIDE_DATA:
            self._side_properties["datamask"] = ~self._flagdata
        return self._side_properties["datamask"]
    
    # ------------------------
    # - Image Data Reduction
    @property
    def _gain(self):
        """This is the calibrated gain"""
        return self.fits[0].header["ATODGN%s"%self.used_amplifier]
    
    @property
    def _readnoise(self):
        """This is the calibrated read noise"""
        return self.fits[0].header["READNSE%s"%self.used_amplifier]
    
    @property
    def _biaslevel(self):
        return self.fits[0].header["BIASLEV%s"%self.used_amplifier]

