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

# filter transmissions from FSPS allfilters.dat
# effective wavelength from sncosmo
GALEX_INFO= {"fuv":{"lbda":1538.620702,"ABmag0":18.82},
             "nuv":{"lbda":2315.663104,"ABmag0":20.08},
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
    PROPERTIES         = ["sky"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    instrument_name = "GALEX"
    INFO            = GALEX_INFO
    
    def set_sky(self, filename=None,
                    skydata=None,
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


    def _aperture_to_photopoint_(self, *args, **kwargs):
        """ convert the aperture output to a photopoints. 
        Galex specialty with SkyBackground """

        if not self.has_sky():
            return super(GALEX,self)._aperture_to_photopoint_( *args, **kwargs )

        datacps, bkgdcps = args
        datacounts = datacps[0]*self.exposuretime
        bkgdcounts = bkgdcps[0]*self.exposuretime
        
        # ------------------
        # - One Photopoint
        if "__iter__" not in dir(datacounts):
            return get_photopoint(lbda=self.lbda,
                                  datacounts=datacounts,
                                  bkgdcounts=bkgdcounts,
                                  exptime=self.exposuretime,
                                  source="image",mjd=self.mjd,
                                  zp=self.mab0,
                                  bandname=self.bandname,
                                  instrument_name=self.instrument_name)
        # -----------------------
        # - Several Photopoints
        return [get_photopoint(lbda=self.lbda,
                                datacounts=datacounts_,
                                bkgdcounts=bkgdcounts_,
                                exptime=self.exposuretime,
                                source="image",mjd=self.mjd,
                                zp=self.mab0,bandname=self.bandname,
                                   instrument_name=self.instrument_name)
                    for datacounts_,bkgdcounts_ in zip(datacounts,bkgdcounts)]
        
    def get_aperture(self, x, y, radius=None,
                     runits="pixels",wcs_coords=False,
                     aptype="circle",subpix=5,
                     ellipse_args={"a":None,"b":None,"theta":None},
                     annular_args={"rin":None,"rout":None},
                     on="data",
                     **kwargs):
        """
        This method uses K. Barary's Sextractor python module SEP
        to measure the aperture photometry. See details here:
        sep.readthedocs.org/en/v0.4.x/apertures.html
        *Remark* you can use the radec option to give x,y as radec coordinates.
        This will only works if you have a wcs solution loaded.
        
        Parameters
        ----------

        x: [float]                  
            Pixel coordinate of the second ("fast") axis
            (so x in the in imshow). In pixels
                                   
        y: [float]                  
            Pixel coordinate of the first ("slow") axis
            (so y in the in imshow). In pixels
                                   
        (other aperture arguments are required be depend on the type of
        aperture photometry *aptype* is choosen)
                                   
        - options - 
        
        on: [string] -optional-
            On which variable should the aperture be made? 
            By default `self.data`. If you are not sure, do not change this.


        aptype: [string]           
            Type of aperture photometry used.
            -circle  => set radius
            -ellipse => set all ellipse_args entries
            -circan  => set all annulus_args entries
            -ellipan => set all ellipse_args entries
            (no other type allowed)
                                               
        
        radius: [float]            
            Size of the circle radius.
            (This is used only if aptype is circle)
                                   
        runits: [str/astropy.units] 
            The unit of the radius (used to convert radius in pixels)

        wcs_coords: [bool]         
            Set True if x,y are ra,dec coordinates
        
        subpix: [int]              
            Division of the real pixel to perform the
            circle to square overlap.

        ellipse_args: [dict]       
            The ellipse parameter that must be filled if
            atype is ellipse of ellipan.

        annular_args: [dict]       
            The annular parameter that must be filled if
            atype is circan of ellipan.

        
        - other options ; not exhautive ; goes to sep.sum_*aptype* - 

        gain: [float]              
            You can manually set the image gain. Otherwise
            this method will look for the key 'self._gain'
            and set None if it does not find it.

        var: [2d-array/float]      
            You can manually set the variance of the image.
            Otherwise, this will look for self.var. If this
            is None and and sepbackground has been created
            the sepbackground.rms()**2 will be used as a
            variance proxy. If not, None is set to var.
                                   
        See other sep options like: 'mask=None, maskthresh=0.0, bkgann=None...'
        
        Return
        ------
        if not sky signal set:
            sum, sumerr, flags (0 if no flag given)
        if a sky signal exists:
            sum, sumerr, flags (0 if no flag given) on raw data +
            sum, sumerr, flags (0 if no flag given) on sky raw data
        """
        propfit = kwargs_update( dict(radius=radius,
                       runits=runits,wcs_coords=wcs_coords,
                       aptype=aptype,subpix=subpix,
                       ellipse_args=ellipse_args,
                       annular_args=annular_args,on=on), **kwargs)
        # = Returns the regular aperture
        if not self.has_sky():
            warnings.warn("No Sky Background set. Only regular aperture photometry returned")
            return super(GALEX, self).get_aperture( x, y, **propfit)
        # = Returns both data and background
        _ = propfit.pop("on",None)
        var = propfit.pop("var",None)
        return [super(GALEX, self).get_aperture( x, y, on="rawdata",var=var, **propfit),
                super(GALEX, self).get_aperture( x, y, on="sky.rawdata",var=var, **propfit)]
    
    # ================== #
    #   Internal         #
    # ================== #    
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
        """ band of the instrument. Change it using set_bandname() """
        if self._properties['bandname'] is None:
            self._properties['bandname'] = "fuv" if "-fd-" in self.filename else "nuv" if "-nd-" in self.filename \
              else "unknown"
        return self._properties['bandname']

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
