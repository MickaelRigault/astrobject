#! /usr/bin/env python
# -*- coding: utf-8 -*-
#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np

# - astropy
from astropy.io      import fits as pf
from astropy         import time

# - local dependencies
from .baseinstrument    import Instrument
from ..photometry       import get_photopoint
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools      import kwargs_update

# filter transmissions from Tonry et al. (2012):
# http://iopscience.iop.org/0004-637X/750/2/99/suppdata/apj425122t3_mrt.txt
# effective wavelengths from sncosmo
PANSTARRS_INFO= {"ps1.g":{"lbda":4866.457871,"ABmag0":25.0},
                 "ps1.r":{"lbda":6214.623038,"ABmag0":25.0},
                 "ps1.g":{"lbda":7544.570357,"ABmag0":25.0},
                 "ps1.g":{"lbda":8679.482571,"ABmag0":25.0},
                 "ps1.g":{"lbda":9633.284241,"ABmag0":25.0},
                 "bands":["ps1.g","ps1.r","ps1.i","ps1.z","ps1.y"],
                 "telescope":{
                     "lon": -156.2571, 
                     "lat": 20.7083}}
"""
Info from Magnier 2016

Seeing:
-------
the median image quality for the 3π survey is FWHM = (1.31, 1.19, 1.11, 1.07, 1.02) 
arcseconds for (gP1,rP1,iP1,zP1,yP1), with a floor of ∼ 0.7 arcseconds


First, 3 PSF-convolved galaxy models (Sersic, DeVaucouleurs, Exponential)
"""
                 
DATAINDEX = 0


# -------------------- #
# - Instrument Info  - #
# -------------------- #
    
def panstarrs(warps_image, weight_image=None, exptime_image=None, **kwargs):
    """ """
    return PanSTARRS(warps_image,weightfilename=weight_image, exptimefilename=exptime_image,
                        **kwargs)

def is_panstarrs_file(filename):
    """This test if the input file is a GALEX one. Test if 'MPSTYPE' is in the header """
    # not great but this is the structure of MJC images
    return "PSCAMERA" in pf.getheader(filename).keys()

def which_band_is_file(filename):
    """This resuts the band of the given file if it is a
    sdss one"""
    if not is_panstarrs_file(filename):
        return None
        
    h_ = pf.getheader(filename)
    return h_.header.get("HIERARCH FPA.FILTER",None).split(".")[0]
    
def which_obs_mjd(filename):
    """ read the galex-filename and return the
    modified julian date """
    if not is_panstarrs_file(filename):
        return None
    h_ = pf.getheader(filename)
    return h_.get("MJD-OBS")



class PanSTARRS( Instrument ):
    """ """
    PROPERTIES = ["weightmap","exptimemap"]
    
    instrument_name = "PanSTARRS"
    INFO            = PANSTARRS_INFO
    
    def __init__(self,filename=None, weightfilename=None,
                maskfilename=None, exptimefilename=None,
                background=None, astrotarget=None,
                data_index=0, dataslice0=None, dataslice1=None,
                 empty=False, **kwargs):
        """
        Initalize the image by giving its filelocation (*filename*). This
        will load it using the load() method.

        Parameters
        ----------
        filename: [string.fits]    
            fits file from where the image will be loaded
            Set None and no image will be loaded

        weightfilename: [string.fits]
            fits file containing the weight map.

        astrotarget: [AstroTarget] An AstroTarget object you which to associate
                                   to this image. 
                                   
        empty: [bool]              Does not do anything, just loads an empty object.
                                   (Careful with that)
                                   
        """
        self.__build__(data_index=data_index)

        if empty:
            return

        if filename is not None:
            force_it = kwargs.pop("force_it",True)
            self.load(filename,force_it=force_it,
                      dataslice0=dataslice0,
                      dataslice1=dataslice1,
                      background=0,
                      **kwargs)

        if weightfilename is not None:
            self.set_weightimage(weightfilename)
            
        if exptimefilename is not None:
            self.set_exptimeimage(exptimefilename)
            
        if maskfilename is not None:
            self.set_datamask(maskfilename)
            
        # - Set the target if any
        if astrotarget is not None:
            self.set_target(astrotarget)
            if self.has_weightimage():
                self.weightimage.set_target(astrotarget)

        # - Background
        self.set_background(background, force_it=True)


    def __build__(self, *args, **kwargs):
        """ """
        super(PanSTARRS, self).__build__(*args, **kwargs)
        self._build_properties["bkgdbox"] = {"bh":100,"bw":100,"fh":3,"fw":3}

                
    def set_weightimage(self, weightmap):
        """ Attach the Weight Maps of the Warps images. 
        (Weight images are variance maps see
        https://confluence.stsci.edu/display/PANSTARRS/PS1+Weight+image)
        """
        self._properties["weightmap"] = PanSTARRS(weightmap, astrotarget=self.target, background=0)

    def set_exptimeimage(self, exptimemap):
        """ Attach the Exposure Time Maps of the Warps images. 
        (Times images are necessary to get the accurate 'counts per second'
        https://confluence.stsci.edu/display/PANSTARRS/PS1+Weight+image)
        """
        self._properties["exptimemap"] = PanSTARRS(exptimemap, astrotarget=self.target, background=0)._sourcedata
        
    def set_datamask(self, maskmap):
        """ The filename of the data containing the mask or the mask itself"""
        if type(maskmap) == str:
            mask = PanSTARRS(maskmap,  background=0)
            super(PanSTARRS,self).set_datamask(mask._sourcedata)
        else:
            super(PanSTARRS,self).set_datamask(maskmap)
        
    # ---------------------
    # - PS structure with wt

    @property
    def weightimage(self):
        """ The weight map. (see set_weightimage())"""
        return self._properties["weightmap"]
    
    def has_weightimage(self):
        """ Has the weight image been set? True means yes """
        return self.weightimage is not None

    @property
    def exptimeimage(self):
        """ The weight map. (see set_weightimage())"""
        return self._properties["exptimemap"]
    
    def has_exptimeimage(self):
        """ Has the weight image been set? True means yes """
        return self.exptimeimage is not None

    @property
    def exposuretime(self):
        """ Effective exposure time. If you loaded a exptime image, this is its data. """
        if not self.has_exptimeimage():
            return super(PanSTARRS, self).exposuretime
        
        return self.exptimeimage

    @property
    def _sourcedata(self):
        """ The data in the PanStarrs images are not per second. """
        return super(PanSTARRS,self).rawdata
    
    @property
    def rawdata(self):
        """ The data in the PanStarrs images are not per second but this is. """
        # see https://confluence.stsci.edu/display/PANSTARRS/PS1+Image+Cutout+Service#PS1ImageCutoutService-Scriptedimagedownloadsandimagecutoutextractions
        return self._sourcedata / self.exposuretime
        #v_per_sec = self._bzero + self._bscale * self._sourcedata / self.exposuretime
        #return self._boffset + self._bsoften * (10**(0.4*v) - 10**(-0.4*v))
    
    @property
    def var(self):
        """ The weight maps since:
           '''Weight images are variance maps'''
           - https://confluence.stsci.edu/display/PANSTARRS/PS1+Weight+image
        """
        if not self.has_weightimage():
            return None
        # Since rawdata have to be devided by exptime -> variance byt exptime**2
        return self.weightimage.rawdata  / self.exposuretime 

    # ------------------
    # Unusual PanSTARRS Stuffs
    @property
    def _bzero(self):
        """ """
        self.header["BZERO"]
        
    @property
    def _bscale(self):
        """ """
        self.header["BSCALE"]
        
    @property
    def _bsoften(self):
        """ """
        self.header["BSOFTEN"]
        
    @property
    def _boffset(self):
        """ """
        self.header["BOFFSET"]

        
    # ---------------------
    # - Generic Properties
    @property
    def mab0(self):
        return self.header["HIERARCH FPA.ZP"]
    
    @property
    def bandname(self):
        """ band of the image. """
        if self._properties['bandname'] is None:
            self._properties['bandname'] = "ps1.%s"%self.header.get("HIERARCH FPA.FILTER",None).split(".")[0]
        return self._properties['bandname']

    @property
    def mjd(self):
        return self.header["MJD-OBS"]
