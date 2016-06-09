#! /usr/bin/env python
# -*- coding: utf-8 -*-

# -- here load all the object that could be parsed
from ..utils.decorators import _autogen_docstring_inheritance

import catalogues
import sdss
import hst
import stella
import ptf

__all__ = ["get_instrument","get_catalogue","fetch_catalogue"]

KNOWN_INSTRUMENTS = ["sdss","hst","stella","ptf"]

def catalogue(source,radec,radius,extracolums=[],column_filters={"rmag":"5..25"},**kwargs):
    print "DECREPATED: catalogue->fetch_catalogue"
    return get_catalogue(source,radec,radius,extracolums=extracolums,
                         column_filters=column_filters,**kwargs)

def fetch_catalogue(source,radec,radius,extracolums=[],column_filters={"rmag":"5..25"},**kwargs):
    """ Download a catalogue from internet (Vizier)
    (Module based on astroquery.)

    Parameters
    ----------

    source: [string]
        Source of the catalogue. Could be: sdss/2mass/wise

    radec: [string / 2D array]
        Coordinates [deg] of the central point of the catalogue.
        Format - string ('%f %f'%(ra,dec)) or 2D array [ra,dec]

    radius: [string]
        Radius around the central location (radec) of the catalogue
        Format: string 'ValueUnits' e.g.: 0.5degree radius would be '0.5d'

    extracolumns: [array of string] -optional-
        Additional column to download for the catalogue. Entries depend on the
        requested catalogue.

    column_filters: [dict] -optional-
        Criteria the requested catalogue entries should have to be downloaded.
        This should be a dictionary with keys being column names (see Vizier)
        Column name depend on the catalogues.
        Example: If you only want stars for the sdss catalogue add  'cl':6

    **kwargs goes to vizier.Vizier

    Return
    ------
    Catalogue
    """
    if type(radec) is not str:
        if len(radec) != 2: raise TypeError("radec must be a string ('ra dec') or a 2D array ([ra,dec])")
        radec = "%f %f"(radec[0],radec[1])
    
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
    
    raise NotImplementedError("Only the SDSS, 2MASS and WISE source catalogues implemented")

def get_catalogue(filename, source,**kwargs):
    """ Reads the given catalogue file and open its corresponding Catalogue

    Parameters:
    -----------
    filename: [string]
        location of the catalogue file to open.
        
    soource: [string]
        Name of the catalogue associated to the data (sdss/2mass/wise/etc)

    **kwargs goes to the Catalogue __init__

    Return
    ------
    Catalogue (corresponding Child's class)
    """
    if source == "sdss":
        return catalogues.SDSSCatalogue(filename,**kwargs)
    
    if source == "2mass":
        return catalogues.MASSCatalogue(filename,**kwargs)
    
    if source.lower() == "wise":
        return catalogues.WISECatalogue(filename,**kwargs)
    
    raise NotImplementedError("Only the SDSS, 2MASS source catalogues implemented")

def get_instrument(filename,astrotarget=None,**kwargs):
    """ Reads the given file and open its corresponding Instrument object.
    Known instruments are (might not be exhaustive): SDSS / HST / PTF  

    Parameters
    ----------
    filename: [string]
        location of the data file to open (a fits file)

    astrotarget: [AstroTarget] -optional-
        Target associated to the image. The target should be within the
        image's boundaries.

    Return
    ------
    Instrument (the corresponding Child's object)
    """
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
    """ Read to filename and return the name of the photometric
    band associated to the given data file.
    Known instruments are (might not be exhaustive): SDSS / HST / PTF
    
    Parameters
    ----------
    filename: [string]
        location of the data file to open (a fits file)

    Return
    ------
    String (name of the photometric band)
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
