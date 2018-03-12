#! /usr/bin/env python
# -*- coding: utf-8 -*-

# -- here load all the object that could be parsed
from ..utils.decorators import _autogen_docstring_inheritance


__all__ = ["get_instrument","get_catalogue","fetch_catalogue"]


from . import catalogues
from . import sdss
from . import hst
from . import ptf
from . import snifs
from . import galex
from . import panstarrs
from . import wise
from . import twomass


KNOWN_INSTRUMENTS = ["sdss","galex","hst","panstarrs","snifs","ptf","stella","wise","twomass"]


def fetch_catalogue(source, radec, radius, extracolumns=[], column_filters={"rmag":"5..25"},**kwargs):
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

def get_instrument(filename, astrotarget=None,**kwargs):
    """ Reads the given file and open its corresponding Instrument object.
    Known instruments are (might not be exhaustive): SDSS / HST / PTF  / GALEX / SNIFS

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
    for instrument in KNOWN_INSTRUMENTS:
        if not eval("%s.is_%s_file(filename)"%(instrument,instrument)):
            continue
        return eval("%s.%s(filename,astrotarget=astrotarget,**kwargs)"%(instrument,instrument))

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
    for instrument in KNOWN_INSTRUMENTS:
        if not eval("%s.is_%s_file(filename)"%(instrument,instrument)):
            continue
        return eval("%s.which_band_is_file(filename)"%instrument)

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
    for instrument in KNOWN_INSTRUMENTS:
        if not eval("%s.is_%s_file(filename)"%(instrument,instrument)):
            continue
        return eval("%s.which_obs_mjd(filename)"%instrument)

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
    from .. import astrometry
    
    index="notset"
    for instrument in KNOWN_INSTRUMENTS:
        if not eval("%s.is_%s_file(filename)"%(instrument,instrument)):
            continue
        index = eval("%s.DATAINDEX"%instrument)
        break
    
    if index == "notset":
        raise ValueError("'filename' does not belong to a known instrument "+"\n"+\
                        "these are:"+", ".join(KNOWN_INSTRUMENTS))
    
    # -- good to go
    return astrometry.wcs(filename,extension=index)
