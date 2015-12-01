#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pyfits as pf
from .baseinstrument import Instrument

__all__ = ["hst"]

HST_INFO = {
    "telescope":{
                 "lon":np.NaN,
                 "lat":np.NaN}
    }

def hst(*args,**kwargs):
    return HST(*args,**kwargs)

def is_hst_file(filename):
    """This test if the given file is an HST one"""
    return True if pf.getheader(filename).get("TELESCOP") == "HST" \
      else False

class HST( Instrument ):
    """This is the umage object custom for HST data"""

    instrument_name = "HST"
    
    def __build__(self,**kargs):
        """
        """
        super(HST,self).__build__()
        # -- How to read the image
        self._build_properties = dict(
                data_index = 1,
                error_index = 1,
                header_exptime = "EXPTIME"
                )
        
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
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.fits[0].header["FILTER"]

    # --------------------
    # - Band Information
    @property
    def mjd_obstime(self):
        """This is the Modify Julien Date at the start of the Exposure"""
        return self.fits[0].header["EXPSTART"]
