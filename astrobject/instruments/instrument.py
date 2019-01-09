#! /usr/bin/env python
# -*- coding: utf-8 -*-

# -- here load all the object that could be parsed
from ..utils.decorators import _autogen_docstring_inheritance


__all__ = ["get_instrument"]

from . import sdss
from . import hst
from . import ptf
from . import snifs
from . import galex
from . import panstarrs
from . import wise
from . import twomass


KNOWN_INSTRUMENTS = ["sdss","galex","hst","panstarrs","snifs","ptf","stella","wise","twomass"]


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
