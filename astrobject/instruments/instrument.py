#! /usr/bin/env python
# -*- coding: utf-8 -*-

# -- here load all the object that could be parsed
from ..utils.decorators import _autogen_docstring_inheritance

import catalogues
import sdss
import hst
import stella
import ptf
import snifs

__all__ = ["get_instrument","get_catalogue","fetch_catalogue"]

KNOWN_INSTRUMENTS = ["sdss","snifs","hst","stella","ptf"]


def fetch_catalogue(source,radec,radius,extracolumns=[],column_filters={"rmag":"5..25"},**kwargs):
    """ Download a catalogue from internet (Vizier)
    (Module based on astroquery.)

    NB: The gaia catalogue DR1 only have G-band magnitude (in vega magnitude) and has no error.
        There is FG and e_FG entries that are the fluxes (and errors) in e-/s in the gaia G-band
    
    Parameters
    ----------
    source: [string]
        Source of the catalogue. Could be: sdss/2mass/wise/gaia

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

    Returns
    -------
    Catalogue
    """
    
    if type(radec) is not str:
        if len(radec) != 2: raise TypeError("radec must be a string ('ra dec') or a 2D array ([ra,dec])")
        radec = "%f %f"(radec[0],radec[1])

    # Test if catalog exist
    if not hasattr(catalogues, "fetch_%s_catalogue"%source.lower()):
        raise ValueError("Unknown catalog %s"%source.lower())
    
    return eval("catalogues.fetch_%s_catalogue(radec,radius,"%source.lower()+\
                                    "extracolumns=extracolumns,"+\
                                    "column_filters=column_filters,**kwargs)")

def get_catalogue(datafile, source, **kwargs):
    """ Reads the given catalogue file and open its corresponding Catalogue

    Parameters:
    -----------
    datafile: [string/ data]
        location of the catalogue file to open or data already opened
        
    soource: [string]
        Name of the catalogue associated to the data (sdss/2mass/wise/etc)
        If the file is already open, you can set None but better set the good
        source if you know it
        
    **kwargs goes to the Catalogue __init__

    Return
    ------
    Catalogue (corresponding Child's class)
    """
    if type(datafile) == str:
        # Test if catalog exist
        if not hasattr(catalogues, "%sCatalogue"%source.upper()):
            raise ValueError("Unknown catalog %s"%source.upper())
        
        return eval("catalogues.%sCatalogue(datafile,**kwargs)"%source.upper())

    return create_catalogue(datafile, source=source)
        
def create_catalogue(datacatalogue, source=None, **kwargs):
    """ load a catalogue given the data. You can set the source if you know it, otherwise
    the function will try to find out """
    from astropy import table
    if table.table.Table not in datacatalogue.__class__.__mro__:
        try:
            datacatalogue = table.Table(datacatalogue)
        except:
            raise TypeError("Cannot convert the input datacatalogue into astropy Table")
    # - SDSS
    if "objID" in datacatalogue.colnames    or (source is not None and source.lower()=="sdss"):
        cat = catalogues.SDSSCatalogue()
    # - Gaia
    elif "Source" in datacatalogue.colnames or (source is not None and source.lower()=="gaia"):
        cat = catalogues.GAIACatalogue()
    elif "Jmag" in datacatalogue.colnames   or (source is not None and source.lower() in ["mass","2mass"]):
        cat = catalogues.MASSCatalogue()
    else:
        raise TypeError("Unknwon catalogue source")
    
    cat.create(datacatalogue, None)
    return cat

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
    # - SNIFS
    if snifs.is_snifs_file(filename):
        return snifs.snifs(filename,astrotarget=astrotarget,
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

    # - SDSS 
    if snifs.is_snifs_file(filename):
        return snifs.which_band_is_file(filename)

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

def which_obs_mjd(filename):
    """ Read to filename and return the observation time of the given data file.
    
    Parameters
    ----------
    filename: [string]
        location of the data file to open (a fits file)

    Return
    ------
    float (MJD)
    """
    # - SDSS 
    if sdss.is_sdss_file(filename):
        return sdss.which_obs_mjd(filename)

    # - SDSS 
    if snifs.is_snifs_file(filename):
        return snifs.which_obs_mjd(filename)

    # - HST     
    if hst.is_hst_file(filename):
        return hst.which_obs_mjd(filename)
    
    # - STELLA
    if stella.is_stella_file(filename):
        return stella.which_obs_mjd(filename)

    # - PTF
    if ptf.is_ptf_file(filename):
        return ptf.which_obs_mjd(filename)
    
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

    # - SNIFS
    elif snifs.is_snifs_file(filename):
        index = snifs.DATAINDEX
    
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
