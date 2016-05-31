#! /usr/bin/env python
# -*- coding: utf-8 -*-

# -- here load all the object that could be parsed
from ..utils.decorators import _autogen_docstring_inheritance

import catalogues
import sdss
import hst
import stella
import ptf

__all__ = ["get_instrument","get_catalogue","fetch_catalogue"]+["instrument","catalogue"] # To be removed

KNOWN_INSTRUMENTS = ["sdss","hst","stella","ptf"]

def catalogue(source,radec,radius,extracolums=[],column_filters={"rmag":"5..25"},**kwargs):
    print "DECREPATED: catalogue->fetch_catalogue"
    return get_catalogue(source,radec,radius,extracolums=extracolums,
                         column_filters=column_filters,**kwargs)

def fetch_catalogue(source,radec,radius,extracolums=[],column_filters={"rmag":"5..25"},**kwargs):
    """
    """
    if source.lower() == "sdss":
        return catalogues.fetch_sdss_catalogue(radec,radius,
                                    extracolums=extracolums,
                                    column_filters=column_filters,**kwargs)
    if source.lower() in ["2mass","mass"]:
        return catalogues.fetch_2mass_catalogue(radec,radius,
                                    extracolums=extracolums,
                                    column_filters=column_filters,**kwargs)
    if source.lower() == "wise":
        return catalogues.fetch_wise_catalogue(radec,radius,
                                    extracolums=extracolums,
                                    column_filters=column_filters,**kwargs)
    
    raise NotImplementedError("Only the SDSS, 2MASS source catalogues implemented")

def get_catalogue(filename, source,**kwargs):
    """
    """
    if source == "sdss":
        return catalogues.SDSSCatalogue(filename,**kwargs)
    
    if source == "2mass":
        return catalogues.MASSCatalogue(filename,**kwargs)
    
    if source.lower() == "wise":
        return catalogues.WISECatalogue(filename,**kwargs)
    
    raise NotImplementedError("Only the SDSS, 2MASS source catalogues implemented")

def instrument(filename,astrotarget=None,**kwargs):
    print "DECREPATED: instrument->get_instrument"
    return get_instrument(filename,astrotarget=astrotarget,
                          **kwargs)

def get_instrument(filename,astrotarget=None,**kwargs):
    """This method parse the input file to know to which
    instrument this filename is associated to.
    This methods will return the corresponding object"""
    
    # - SDSS 
    if sdss.is_sdss_file(filename):
        return sdss.sdss(filename,astrotarget=astrotarget,
                    **kwargs)
    # - HST     
    if hst.is_hst_file(filename):
        return hst.hst(filename,astrotarget=astrotarget,
                    **kwargs)
    # - STELLA
    if stella.is_stella_file(filename):
        return stella.stella(filename,astrotarget=astrotarget,
                    **kwargs)
    # - PTF
    if ptf.is_ptf_file(filename):
        return ptf.ptf(filename,astrotarget=astrotarget,
                    **kwargs)
    
    # - Nothing else...
    raise ValueError("'filename' does not belong to a known instrument "+"\n"+\
                     "these are:"+", ".join(KNOWN_INSTRUMENTS))

def which_band_is_file(filename):
    """
    """
    # - SDSS 
    if sdss.is_sdss_file(filename):
        return sdss.which_band_is_file(filename)
    
    # - HST     
    if hst.is_hst_file(filename):
        return hst.which_band_is_file(filename)
    
    # - STELLA
    if stella.is_stella_file(filename):
        return stella.which_band_is_file(filename)

    # - PTF
    if ptf.is_ptf_file(filename):
        return ptf.which_band_is_file(filename)
    
    # - Nothing else...
    raise ValueError("'filename' does not belong to a known instrument "+"\n"+\
                     "these are:"+", ".join(KNOWN_INSTRUMENTS))
                     
def is_known_instrument_file(filename):
    """
    This function test if the given filename is a known object that
    could be loaded by the function *instrument*.
    This try to return an empty instrument to do so, this will loop
    through the known *is_`inst`_file*

    Return
    ------
    bool
    """
    try:
        inst = get_instrument(filename,empty=True)
    except:
        return False
    return True

def get_instrument_wcs(filename):
    """
    This function enables to load the wcs solution only, not opening the
    full image. This might be useful to avoid opening large images
    """
    from ..  import astrometry
    # - SDSS 
    if sdss.is_sdss_file(filename):
        index = sdss.DATAINDEX
    
    # - HST     
    elif hst.is_hst_file(filename):
        index = hst.DATAINDEX
    # - STELLA
    elif stella.is_stella_file(filename):
        index = stella.DATAINDEX
    # - PTF
    elif ptf.is_ptf_file(filename):
        index = ptf.DATAINDEX
        
    else:
        # - Nothing else...
        raise ValueError("'filename' does not belong to a known instrument "+"\n"+\
                        "these are:"+", ".join(KNOWN_INSTRUMENTS))
    # -- good to go
    return astrometry.wcs(filename,extension=index)
