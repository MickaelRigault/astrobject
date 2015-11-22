#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from ..astrobject.photometry import pf
from .baseinstrument import Instrument

ZTF_INFO = {"u":{"lbda":3551,"ABmag0":22.46},
            "g":{"lbda":4686,"ABmag0":22.50},
            "r":{"lbda":6166,"ABmag0":22.50},
            "i":{"lbda":7480,"ABmag0":22.50},
            "z":{"lbda":8932,"ABmag0":22.52},
            }

class ZTF( Instrument ):
    """This class is a temporary class containing the
    ZTF observations"""
    
    instrument_name = "ZTF"
    
    def __build__(self,data_index=0):
        """
        """
        # Nothing more here but ready to add info
        # -- Load the basic builds
        super(ZTF,self).__build__()
        self._build_properties = dict(
            data_index = data_index,
            header_exptime = "EXPTIME"
            )
    # --------------
    # - Image Data
    # --------------  
    @property
    def band(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.header["FILTER"]

    @property
    def band_info(self):
        return ZTF_INFO[self.band]

    @property
    def _gain(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.header["GAIN"]


###############################
#                             #
# Simulated object: Fake ZTF  #
#                             #
###############################
class SimulatedZTF( ZTF ):
    """This module enables to return create a fake ZTF Instrument"""
    def __init__(self,band_info=None,
                 exptime=None,gain=None,
                 empty=False):
        """
        
        ..example:
              ztf.SimulatedZTF({"lbda":4000,"ABmag0":22,"name":"ztfg"},
                               exptime=60,gain=1)
        
        """
        self.__build__()
        if empty:
            return
        
        self.set_inst_header(band_info,
                             exptime,gain)
        
    def __build__(self,**kwargs):
        """
        """
        self._side_properties_keys.append("band_info")
        super(SimulatedZTF,self).__build__(**kwargs)
        
    # =============================== #
    # = Module to fake the header   = #
    # =============================== #
    def set_inst_header(self,band_info,exptime,gain):
        """
        """
        if type(band_info) is not dict or \
          "ABmag0" not in band_info.keys() or "lbda" not in band_info.keys() \
          or "name" not in band_info.keys():
            raise TypeError("'band_info' must be dict having at least ABmag0, lbda, and name")

        self._side_properties['band_info'] = band_info
        if self.header is None:
            self._properties['header'] = pf.Header()

        self.header['GAIN'] = (gain, "Faked Gain")
        self.header['FILTER'] = (band_info["name"], "Faked band")
        self.header['EXPTIME'] = (exptime, "Faked band")
        self.header['INSTRUM'] = ("SimulatedZTF", "Object created for test")
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def band_info(self):
        return self._side_properties['band_info']
    
    # - band in haeder['Filter'] is ZTF default
    @property
    def lbda(self):
        return self.band_info['lbda']
        
    @lbda.setter
    def lbda(self,value):
        self._side_properties['band_info']['lbda'] = value

    @property
    def mab0(self):
        return self._side_properties['band_info']['ABmag0']
        
    @mab0.setter
    def mab0(self,value):
        self._side_properties['band_info']['ABmag0'] = value

    @property
    def exposuretime(self):
        if self.header is None:
            raise AttributeError("no 'header' loaded")        
        return self.header[self._build_properties['header_exptime']]
    
    @exposuretime.setter
    def exposuretime(self,value):
        self.header.update(self._build_properties['header_exptime'],value)
    
    
