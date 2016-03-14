#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pyfits as pf
from .baseinstrument import Instrument,get_bandpass
from ...utils.decorators import _autogen_docstring_inheritance

__all__ = ["sdss","SDSS_INFO"]

SDSS_INFO = {"u":{"lbda":3551,"ABmag0":22.46},
             "g":{"lbda":4686,"ABmag0":22.50},
             "r":{"lbda":6166,"ABmag0":22.50},
             "i":{"lbda":7480,"ABmag0":22.50},
             "z":{"lbda":8932,"ABmag0":22.52},
             "bands":["u","g","r","i","z"],
             "telescope":{
                 "lon": 32.780,
                 "lat":-105.82}
            }

DATAINDEX = 0

# -------------------- #
# - Instrument Info  - #
# -------------------- #
    
def sdss(*args,**kwargs):
    return SDSS(*args,**kwargs)

def is_sdss_file(filename):
    """This test if the input file is a SDSS one"""
    return True if pf.getheader(filename).get("ORIGIN") == "SDSS" \
      else False

def which_band_is_file(filename):
    """This resuts the band of the given file if it is a
    sdss one"""
    if not is_sdss_file(filename):
        return None
    return pf.getheader(filename).get("FILTER")

    
# -------------------- #
# - Inside tools     - #
# -------------------- #
def get_darkvariance(camcol,band,run=None):
    """
    http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    """
    DARK_VAR_CCD = {
            0:{"u":9.61,   "g":15.6025,"r":1.8225,
               "i":7.84,     "z":0.81},
            1:{"u":12.6025,"g":1.44,   "r":1.00,
               "i":[5.76,6.25],"z":1.0},
            2:{"u":8.7025, "g":1.3225, "r":1.3225,
               "i":4.6225,   "z":1.0},
            3:{"u":12.6025,"g":1.96,   "r":1.3225,
               "i":[6.25,7.5625],"z":[9.61,12.6025]},
            4:{"u":9.3025, "g":1.1025, "r":0.81,
               "i":7.84,     "z":[1.8225,2.1025]},
            5:{"u":7.0225, "g":1.8225, "r":0.9025,
               "i":5.0625,   "z":1.21}
            }
         
    dark = DARK_VAR_CCD[camcol-1][band]
    # ----------
    # - output
    if type(dark) == float:
        return dark
    if run is None:
        raise ValueError("there is two dark-variance possibilites for "+\
                         " *camcol* %d, *band* %s "%(
                            camcol-1,band) + "Please, provide a *run*")
    
    return dark[1] if run>1500 else dark[0]
    
def get_gain(camcol,band,run=None):
    """
    http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    """
    GAIN_CCD = {
            0:{"u":1.62, "g":3.32, "r":4.71,
               "i":5.165,"z":4.745},
            1:{"u":[1.595,1.825],"g":3.855,"r":4.6,
               "i":6.565,"z":5.155},
            2:{"u":1.59,"g":3.845,"r":4.72,
               "i":4.86,"z":4.885},
            3:{"u":1.6,"g":3.995,"r":4.76,
               "i":4.885,"z":4.775},
            4:{"u":1.47,"g":4.05,"r":4.725,
               "i":4.64,"z":3.48},
            5:{"u":2.17,"g":4.035,"r":4.895,
               "i":4.76,"z":4.69}}
    
    gain = GAIN_CCD[camcol-1][band]
    # ----------
    # - output
    if type(gain) == float:
        return gain
    if run is None:
        raise ValueError("there is two gain possibilites for *camcol* %d, *band* %s "%(
            camcol-1,band) + "Please, provide a *run*")
    
    return gain[1] if run>1100 else gain[0]

#######################################
#                                     #
#   SDSS Image Object                 #
#                                     #
#######################################
class SDSS( Instrument ):
    """
    This is the image object custom for SDSS data
    This method is based on
    http://www.sdss.org/dr12/algorithms/magnitudes/
    """
    instrument_name = "SDSS"
    
    def __build__(self,data_index=DATAINDEX):
        """
        """
        # -- Load the basic builds
        self._derived_properties_keys.append("error")
        self._derived_properties_keys.append("sky")
        self._derived_properties_keys.append("skyparam")
        super(SDSS,self).__build__(data_index=data_index)

    @_autogen_docstring_inheritance(Instrument.set_catalogue,"Instrument.set_catalogue")
    def set_catalogue(self,catalogue,force_it=True,**kwargs):
        #
        # - Add the bandname key_mag setting
        #
        if catalogue.source_name =="SDSS":
            key_mag = "%smag"%self.bandname[-1]
            key_magerr = "e_%smag"%self.bandname[-1]
            if key_mag not in catalogue.data.keys():
                print "WARNING No %s in the catalogue data. Cannot assign a key_mag"%key_mag
            catalogue.set_mag_keys(key_mag,key_magerr)
            
        super(SDSS,self).set_catalogue(catalogue,force_it=force_it,**kwargs)
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    # --------------
    # - Sky Data 
    @property
    def sky(self):
        if self._derived_properties["sky"] is None:
            self._define_sky_(force_it=True)
        return self._derived_properties["sky"]

    @property
    def _rawsky_prop(self):
        if self.fits is None:
            raise AttributeError("no 'fits' loaded")
        # -- if first time you call it, load it.
        if self._derived_properties["skyparam"] is None:
            self._derived_properties["skyparam"] = self.fits[2].data[0]
        # -- return it
        return self._derived_properties["skyparam"]

    # --------------------
    # - Band Information
    @property
    def bandname(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return "sdss"+self.header["FILTER"]

    @property
    def mjd_obstime(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        from astropy import time
        dateobs = self.header["DATE-OBS"] if "/" not in self.header["DATE-OBS"]\
          else "19"+"-".join(self.header["DATE-OBS"].split("/")[::-1])
            
        return time.Time("%sT%s"%(dateobs,self.header["TAIHMS"])).mjd
        

    # -- Derived values
    @property
    def mab0(self):
        """The ABmag zero point of SDSS data"""
        if self.bandname  == "sdssu":
            return 22.46
        if self.bandname  == "sdssz":
            return 22.52
        return 22.5
    
    # ------------
    # - Internal
    @property
    def _cimg(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.header["NMGY"]
    
    @property
    def _gain(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return get_gain(self.header["CAMCOL"],self.bandname[-1],
                        self.header["RUN"])

    @property
    def _darkvariance(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return get_darkvariance(self.header["CAMCOL"],self.bandname[-1],
                                self.header["RUN"])
    # =========================== #
    # = Internal Tools          = #
    # =========================== #
    # -----------------
    # - Sky Background
    def _get_brightsource_mask_(self,stellar_uppermag=15,gal_uppermag=15,**kwargs):
        """
        """
        if not self.has_catalogue():
            self.download_catalogue(column_filters=\
                                    {"%smag"%self.sdss.bandname[-1]:"12..20"},
                                    )
        brightstars_coords = [self.catalogue.ra[(self.catalogue.mag>stellar_uppermag)\
                                                & self.catalogue.starmask],
                             self.catalogue.dec[(self.catalogue.mag>stellar_uppermag)\
                                                 & self.catalogue.starmask],
                             self.catalogue.mag[(self.catalogue.mag>stellar_uppermag)\
                                                 & self.catalogue.starmask],
                                                 ]
        brightgal_coords  = [self.catalogue.ra[(self.catalogue.mag>gal_uppermag)\
                                                & ~self.catalogue.starmask],
                              self.catalogue.dec[(self.catalogue.mag>gal_uppermag)\
                                                 & ~self.catalogue.starmask],
                            self.catalogue.mag[(self.catalogue.mag>gal_uppermag)\
                                                 & ~self.catalogue.starmask],
                                                 ]
        def _mag_to_pixel_mask_(self,mag):
            """
            based on the 15mag->32pixels 12->1600 pixel from
            http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1105.1960
            """
            return 0 if mag>15 else  int(1600 - (mag - 12)*522)
            
        return [[ra,dec,self._mag_to_pixel_mask_(mag)]
                for ra,dec,mag in np.asarray(brightstars_coords).T if mag<15]
    
    def _define_sky_(self,force_it=False):
        """This methods convert the sky background registered in the fits file
        such that it maps the data values"""
        from scipy import interpolate
        
        if force_it is False and self.sky is not None:
            return

        rawsky,yinterpol,xinterpol = self._rawsky_prop
        xshape,yshape = np.shape(rawsky)
        XI, YI = np.meshgrid( range(xshape), range(yshape))
        t = interpolate.bisplrep(XI,YI,rawsky.T)
        self._derived_properties["sky"] = interpolate.bisplev(xinterpol,yinterpol,t)
        return t

