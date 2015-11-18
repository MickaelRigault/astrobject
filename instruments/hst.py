#! /usr/bin/env python
# -*- coding: utf-8 -*-

from .baseinstrument import Instrument

__all__ = ["hst"]

HST_INFO = {
    "f225w":{"lbda":2359,"ABmag0":None,"instrument":"UVIS"}
    }

def hst(*args,**kwargs):
    return HST(*args,**kwargs)

def is_hst_file(filename):
    """This test if the given file is an HST one"""
    raise NotImplementedError("To Be Done")

class HST( Instrument ):
    """This is the umage object custom for HST data"""

    instrument_name = "HST"
    
    def __build__(self):
        """
        """
        super(HST,self).__build__()
        # -- How to read the image
        self._build_properties = dict(
                data_index = 1,
                error_index = 1,
                #header_exptime = "EXPTIME"
                )
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    # --------------
    # - Image Data    
    @property
    def band(self):
        if self.header is None:
            raise AttributeError("no header loaded ")
        return self.header["FILTER"]

    @property
    def band_info(self):
        return HST_INFO[self.band]

    @property
    def lbda(self):
        return self.band_info["lbda"]
