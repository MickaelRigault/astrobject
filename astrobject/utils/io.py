#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module has the basic functionalities for the few file io of astrobject"""
# Note to dev. This shall eventually be removed ; all linked to sncosmo and astropy

import os
import numpy as np
import warnings
from glob import glob
            
DEFAULT_SFD_FILES = ["SFD_dust_4096_ngp.fits","SFD_dust_4096_sgp.fits"]

############################################
#                                          #
#  Dust maps IO                            #
#                                          #
############################################
def get_default_sfd98_dir(download_if_needed=True):
    """
    """
    # -------------------------
    # - Reads astropy/sncosmo
    a = [l for l in open(os.path.expanduser('~')+"/.astropy/config/sncosmo.cfg").read().splitlines()
         if l.startswith("sfd98_dir")]
    if len(a) == 0:
        raise IOError("You need to define *sfd98_dir* in the sncosmo/astropy config file: see ~/.astropy/config/sncosmo.cfg")
    
    sfd98_dir = a[-1].split('=')[-1].replace(" ","")

    # -------------------------
    # - Do you have the files
    if not test_dustdir(sfd98_dir):
        if download_if_needed:
            download_sf98_dust_map(sfd98_dir)
        else:
            warnings.warn("No dust map downloaded yet. Do so to have access to SFD98 data")
            return None
        
    return sfd98_dir


def download_sf98_dust_map(outdir,verbose=True):
    """
    This downloads the dustmap from sncosmo
    """
    try:
        import wget
    except:
        raise ImportError("no wget installed. pip install wget")

    if verbose:
        print("Dust map downloading starting. 4 files, 167 Mo in total")

    try:
        os.makedirs(outdir)
        if verbose: print("Dust files saved here: %s"%outdir)
    except:
        if verbose: print("%s path already exists. Dust map downloading their"%outdir)
            
    for f in ["SFD_dust_4096_ngp.fits","SFD_mask_4096_ngp.fits",
              "SFD_dust_4096_sgp.fits","SFD_mask_4096_sgp.fits"]:
        wget.download("http://sncosmo.github.io/data/dust/"+f,out=outdir)

    if verbose:
        print("SF98 dust map downloaded.")


def test_dustdir(dustdir,hard_test=True):
    """
    """

    filesin = glob(dustdir+"/*.fits")
    if len(filesin)==0:
        print("no data in %s"%dustdir)
        return False
    if hard_test and np.asarray([dustdir+"/"+f not in filesin
                                 for f in DEFAULT_SFD_FILES]).any():
        return False
    return True


