#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from .baseinstrument import Instrument
from ..utils.decorators import _autogen_docstring_inheritance

__all__ = ["sdss","SDSS_INFO"]

# filter transmission & effective wavelength from sncosmo
SDSS_INFO = {"u":{"lbda":3594.325367,"ABmag0":22.46},
             "g":{"lbda":4717.599777,"ABmag0":22.50},
             "r":{"lbda":6186.799970,"ABmag0":22.50},
             "i":{"lbda":7506.238080,"ABmag0":22.50},
             "z":{"lbda":8918.301484,"ABmag0":22.52},
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
    return pf.getheader(filename).get("ORIGIN") == "SDSS"

def which_band_is_file(filename):
    """This resuts the band of the given file if it is a
    sdss one"""
    if not is_sdss_file(filename):
        return None
    return pf.getheader(filename).get("FILTER")

def which_obs_mjd(filename):
    """ read the sdss-filename and return the
    modified julian date """
    if not is_sdss_file(filename):
        return None
    return get_mjd(pf.getheader(filename))

# -------------------- #
# - Inside tools     - #
# -------------------- #
def get_mjd(sdssheader):
    """ read the header of the sdss file and return the
    correct modified julian date"""
    from astropy import time
    dateobs = sdssheader["DATE-OBS"] if "/" not in sdssheader["DATE-OBS"]\
      else "19"+"-".join(sdssheader["DATE-OBS"].split("/")[::-1])
    return time.Time("%sT%s"%(dateobs,sdssheader["TAIHMS"])).mjd
    
    
def get_darkvariance(camcol,band,run=None):
    """
    data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
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
    data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
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


def hdu2_to_sky(hdu2, kind="linear"):
    """ This converts the hdu2 entry of the sdss frame (sky) to a good sky array """
    from scipy import interpolate
    sky, x, y = hdu2.data[0]
    xold = np.arange(np.shape(sky)[1])
    yold = np.arange(np.shape(sky)[0])
    tck  = interpolate.interp2d(xold,yold,sky,kind=kind)
    return tck(x,y)

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

    fits composition:
    -----------------
    (data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html)

    
    ..hdu-0: the corrected frame, what is normally in the "fpC" files
        The "image", a 2048x1489 array of floating point values, the calibrated and
        sky-subtracted version of the fpC "corrected frame" files produced by photo.
        Units are in nanomaggies.
        The header has additional WCS information.
        These have been altered from the fpC versions to correct for any offsets
        from the astrom solutions found in the asTrans files.
    
    ..hdu-1: the flat-field and calibration vector
        The "calibvec", a 2048-element array of "float32" values, encompassing
        the flat-field correction to apply, multiplied by the calibration.
        Translates the counts in the original image into nanomaggies.
        The calibrations have ALREADY been applied to the HDU0, so this calibration
        vector is only to be used to decalibrate an image back into counts.
     
    ..hdu2: the sky image
        The "sky", an approximately 256x192 array of "float32" values (there are some
        variations around the y-size of the array), with information about how to
        interpolate it to the full image size.
        These sky values are determined from the global sky fits across the run
        (not photo). There are cases at the end of runs where the y-size of ALLSKY
        is less than 186, which requires one to extrapolate off the end of the ALLSKY
        image in the y-direction when using XINTERP and YINTERP.
        This extrapolation should be a constant extrapolation (not linear).
        The sky values have ALREADY been subtracted from HDU0, so this sky estimate is
        only to be used to return an image to its (near) original set of values.

    ..hdu3:  the asTrans structure
        Detailed astrometric information for the field. Basically the asTrans
        structure as found in the asTrans files.
        (data.sdss3.org/datamodel/files/PHOTO_REDUX/RERUN/RUN/astrom/asTrans.html)
    
    """
    instrument_name = "SDSS"
    INFO            = SDSS_INFO
    
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["sky"]

    def __build__(self,data_index=DATAINDEX):
        """
        """
        # -- Load the basic builds
        super(SDSS,self).__build__(data_index=data_index)
        self._build_properties["bkgdbox"]['bh'] = 200
        self._build_properties["bkgdbox"]['bw'] = 200
        
    @_autogen_docstring_inheritance(Instrument.set_catalogue,"Image.set_catalogue")
    def set_catalogue(self,catalogue, **kwargs):
        #
        # - Add the bandname key_mag setting
        #
        if catalogue.source_name == "SDSS":
            key_mag = "%smag"%self.bandname[-1]
            key_magerr = "e_%smag"%self.bandname[-1]
            if key_mag not in catalogue.data.keys():
                warnings.warn("WARNING No %s in the catalogue data. Cannot assign a key_mag"%key_mag)
            catalogue.set_mag_keys(key_mag,key_magerr)
            
        super(SDSS,self).set_catalogue(catalogue,**kwargs)
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    # --------------
    # - Sky Data 
    @property
    def sky(self):
        """ The sky background of the image. This is in counts, shaped to follow the data slicing """
        if self._derived_properties["sky"] is None:
            self._update_sky_()
        return self._derived_properties["sky"]
            
    @property
    def var(self):
        """ Variance estimated from the SDSS rawdata, skynoise and Dark Variance.
        see http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html """
        if self._derived_properties["var"] is None:
            self._update_var_()
        return self._derived_properties["var"]
        
    @property
    def _dataunits_to_electron(self):
        """ No _dataunits_to_electron for sep sum_circle ; the accurate variance is provided """
        return None
    
    # --------------------
    # - Band Information
    @property
    def bandname(self):
        if self._properties['bandname'] is None:
            if self.header is None:
                raise AttributeError("no header loaded ")
            self._properties['bandname'] = "sdss"+self.header["FILTER"]
        return self._properties['bandname']

    @property
    def mjd(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return get_mjd(self.header)
        

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
        if self.fits is None:
            raise AttributeError("no fits loaded ")
        return np.asarray([self.fits[1].data[self._build_properties['dataslice1'][0]:\
                                             self._build_properties['dataslice1'][1]].tolist()]\
          * self._dataslicing[2])
    
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
    # - Update hacking
    def _update_(self,*args,**kwargs):
        if self.fits is not None:
            self._update_sky_()
        super(SDSS,self)._update_(*args,**kwargs)
    
    def _update_data_(self,*args,**kwargs):
        super(SDSS,self)._update_data_(*args,**kwargs)
        self._update_var_()

    def _update_sky_(self):
        self._derived_properties["sky"] = \
          hdu2_to_sky(self.fits[2])[self._build_properties['dataslice0'][0]:self._build_properties['dataslice0'][1],
                                    self._build_properties['dataslice1'][0]:self._build_properties['dataslice1'][1]]

    def _update_var_(self):
        """
        self.rawdata / self._cimg -> data in counts
        self.rawdata / self._cimg + self.sky -> original counts
        (self.rawdata / self._cimg + self.sky) / gain -> original photoelectron
        """
        self._derived_properties["var"] = \
          dn_err = ((self.rawdata / self._cimg + self.sky) / self._gain  + self._darkvariance) * self._cimg **2

    # -------------------
    # - Background hacking
    def _get_default_background_(self,*args,**kwargs):
        return np.zeros(np.shape(self.rawdata))
    
