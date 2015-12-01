#! /usr/bin/env python
# -*- coding: utf-8 -*-

# -- here load all the object that could be parsed
from ...utils.decorators import _autogen_docstring_inheritance

from .sdss   import sdss,SDSS_INFO,is_sdss_file
from .hst    import hst,HST_INFO,is_hst_file
from .stella import stella,STELLA_INFO,is_stella_file

__all__ = ["instrument","KNOWN_INSTRUMENTS"]

KNOWN_INSTRUMENTS = ["sdss","hst","stella"]


def instrument(filename,astrotarget=None,**kwargs):
    """This method parse the input file to know to which
    instrument this filename is associated to.
    This methods will return the corresponding object"""
    
    # - SDSS 
    if is_sdss_file(filename):
        return sdss(filename,astrotarget=astrotarget,
                    **kwargs)
    # - HST     
    if is_hst_file(filename):
        return hst(filename,astrotarget=astrotarget,
                    **kwargs)
    # - STELLA
    if is_stella_file(filename):
        return stella(filename,astrotarget=astrotarget,
                    **kwargs)
    
    # - Nothing else...
    raise ValueError("'filename' does not belongs to a known instrument "+"\n"+\
                     "these are:"+", ".join(KNOWN_INSTRUMENTS))
                     
