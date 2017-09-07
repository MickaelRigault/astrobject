#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits as pf
from .baseinstrument import Instrument

ZTF_INFO = {"g":{"lbda":4686,"ABmag0":22.50},
            "r":{"lbda":6166,"ABmag0":22.50},
            "i":{"lbda":7480,"ABmag0":22.50},
            "telescope":{
                 "lon": 32.780,
                 "lat":-105.82
                }
            }

class ZTF( Instrument ):
    """This class is a temporary class containing the
    ZTF observations"""
    
    instrument_name = "ZTF"
    INFO            = ZTF_INFO
    
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

