#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module has the basic functionalities for the few file io of astrobject"""
# Note to dev. This shall eventually be removed ; all linked to sncosmo and astropy

import os
import numpy as np
import warnings
from glob import glob
            
DEFAULT_SF98_DUSTDIR = os.path.dirname(__file__)+"/data/dustmap/"
DEFAULT_SFD_FILES = ["SFD_dust_4096_ngp.fits","SFD_mask_4096_ngp.fits",
                     "SFD_dust_4096_sgp.fits","SFD_mask_4096_sgp.fits"]

############################################
#                                          #
#  Dust maps IO                            #
#                                          #
############################################
def get_default_sfd98_dir(download_if_needed=True):
    """
    """
    if not test_dustdir(DEFAULT_SF98_DUSTDIR):
        if dowload_if_needed:
            download_sf98_dust_map()
        else:
            warnings.warn("No dust map downloaded yet. Do so to have access to SFD98 data")
            return None
        
    return DEFAULT_SF98_DUSTDIR
        

def download_sf98_dust_map(verbose=True):
    """
    This downloads the dustmap from sncosmo
    """
    try:
        import wget
    except:
        raise ImportError("no wget installed. pip install wget")

    if verbose:
        print "Dust map downloading starting. 4 files, 167 Mo in total"

    try:
        os.makedirs(DEFAULT_SF98_DUSTDIR)
        if verbose: print "Dust files saved here: %s"%DEFAULT_SF98_DUSTDIR
    except:
        if verbose: print "%s path already exists. Dust map downloading their"%DEFAULT_SF98_DUSTDIR
    for f in ["SFD_dust_4096_ngp.fits","SFD_mask_4096_ngp.fits",
              "SFD_dust_4096_sgp.fits","SFD_mask_4096_sgp.fits"]:
        wget.download("http://sncosmo.github.io/data/dust/"+f,out=DEFAULT_SF98_DUSTDIR)

    if verbose:
        print "SF98 dust map downloaded."


def test_dustdir(dustdir=DEFAULT_SF98_DUSTDIR,
                 hard_test=True):
    """
    """

    filesin = glob(dustdir+"*.fits")
    if len(filesin)==0:
        return False
    if hard_test and np.asarray([dustdir+f not in filesin
                                 for f in DEFAULT_SFD_FILES]).any():
        return False
    return True


