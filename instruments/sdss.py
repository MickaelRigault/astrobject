#! /usr/bin/env python
# -*- coding: utf-8 -*-

from ..astrobject.photometry import Image


__all__ = ["sdss"]

SDSS_INFO = {
     "u":{"lbda":3551,"ABmag0":22.46},
     "g":{"lbda":4686,"ABmag0":22.50},
     "r":{"lbda":6166,"ABmag0":22.50},
     "i":{"lbda":7480,"ABmag0":22.50},
     "z":{"lbda":8932,"ABmag0":22.52},
    }
# -------------------- #
# - Instrument Info  - #
# -------------------- #
    
def sdss(*args,**kwargs):
    return SDSS(*args,**kwargs)

# -------------------- #
# - Inside tools     - #
# -------------------- #
def get_darkvariance(camcol,band,run=None):
    """
    http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    """
    DARK_VAR_CCD = {
            0:{"u":9.61,   "g":15.6025,"r":1.8225,
               "i":7.84,"z":0.81},
            1:{"u":12.6025,"g":1.44,   "r":1.00,
               "i":[5.76,6.25],"z":1.0},
            2:{"u":8.7025, "g":1.3225, "r":1.3225,
               "i":4.6225,"z":1.0},
            3:{"u":12.6025,"g":1.96,   "r":1.3225,
               "i":[6.25,7.5625],"z":[9.61,12.6025]},
            4:{"u":9.3025, "g":1.1025, "r":0.81,"i":7.84,
               "z":[1.8225,2.1025]},
            5:{"u":7.0225, "g":1.8225, "r":0.9025,
               "i":5.0625,"z":1.21}
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
class SDSS( Image ):
    """This is the image object custom for SDSS data"""
    def __build__(self):
        """
        """
        print "NOT READY YET"
        # -- Load the basic builds
        self._derived_properties_keys.append("error")
        super(SDSS,self).__build__()


    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def band(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.header["FILTER"]
    
    @property
    def _cimg(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.header["NMGY"]
    
    @property
    def _gain(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return get_gain(self.header["CAMCOL"],self.band,
                        self.header["RUN"])

    @property
    def _darkvariance(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return get_darkvariance(self.header["CAMCOL"],self.band,
                                self.header["RUN"])
    # =========================== #
    # = Internal Tools          = #
    # =========================== #
    def _update_error_(self):
        """
        This method follows the indication provided by SDSS' guide:
        data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
        """
        nelect = self._gain*(self.rawdata / self._cimg + 0) # should be skyimg)
        self._derived_properties["errors"] = None
        
