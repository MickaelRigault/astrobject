#! /usr/bin/env python
# -*- coding: utf-8 -*-

from ..astrobject.photometry import Image


__all__ = ["hst"]

SDSS_INFO {
    "f225w":{"lbda":2359,"ABmag0":None,"instrument":"UVIS"},
    }

def hst(*args,**kwargs):
    return HST(*args,**kwargs)


class HST( Image ):
    """This is the umage object custom for HST data"""
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
        
