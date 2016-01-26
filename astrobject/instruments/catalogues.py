#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from sncosmo import get_bandpass

from .baseinstrument import Catalogue
# -- here load all the object that could be parsed
from ...utils.decorators import _autogen_docstring_inheritance


#################################
#                               #
# BASIC SDSS: Catalogue         #
#                               #
#################################
def fetch_sdss_catalogue(center,radius,extracolums=[],column_filters={"rmag":"13..22"}):
    """
    """
    from .sdss import SDSS_INFO
    try:
        from astroquery import vizier
    except:
        raise ImportError("install astroquery. (pip install astroquery)")
    
    # -----------
    # - DL info
    columns = ["cl","objID",
               "RAJ2000","e_RAJ2000","DEJ2000","e_DEJ2000",
               #"ObsDate","Q"#"mode","SDSS9",
               ]
    for band in SDSS_INFO["bands"]:
        columns.append("%smag"%band)
        columns.append("e_%smag"%band)
    
    columns = columns+extracolums
    # - WARNING if discovered that some of the bandmag were missing if too many colums requested
    c = vizier.Vizier(catalog="V/139", columns=columns, column_filters=column_filters)
    c.ROW_LIMIT = 100000
    try:
        t = c.query_region(center,radius=radius).values()[0]
    except:
        raise IOError("Error while querying the given coords. You might not have an internet connection")
    
    cat = SDSSCatalogue(empty=True)
    cat.create(t.columns,None,
               key_class="cl",value_star=6,
               key_ra="RAJ2000",key_dec="DEJ2000")
    return cat

# ------------------- #
# - SDSS CATALOGUE  - #
# ------------------- #
class SDSSCatalogue( Catalogue ):
    """
    """
    source_name = "SDSS"
    
    def __init__(self, catalogue_file=None,empty=False,
                 key_mag=None,key_magerr=None,key_ra=None,key_dec=None):
        """
        """
        self.__build__(data_index=2,key_mag=key_mag,
                       key_magerr=key_magerr,
                       key_ra=key_ra,key_dec=key_dec)
        if empty:
            return
        
        self.load(catalogue_file)        
    
    @_autogen_docstring_inheritance(Catalogue.set_mag_keys,"Catalogue.set_mag_keys")
    def set_mag_keys(self,key_mag,key_magerr):
        #
        # add lbda def
        #
        super(SDSSCatalogue,self).set_mag_keys(key_mag,key_magerr)
        if key_mag is not None:
            bandpass = get_bandpass("sdss%s"%key_mag[0])
            self.lbda = bandpass.wave_eff
    
#################################
#                               #
# BASIC 2MASS: Catalogue        #
#                               #
#################################
def fetch_2mass_catalogue(center,radius,extracolums=[],column_filters={"Jmag":"10..30"}):
    """
    """
    try:
        from astroquery import vizier
    except:
        raise ImportError("install astroquery. (pip install astroquery)")
    
    # -----------
    # - DL info
    columns = ["2MASS",
               "RAJ2000","DEJ2000",
               ]
        
    for band in ["J","H","K"]:
        columns.append("%smag"%band)
        columns.append("e_%smag"%band)
    
    columns = columns+extracolums
    # - WARNING if discovered that some of the bandmag were missing if too many colums requested
    c = vizier.Vizier(catalog="II/246", columns=columns, column_filters=column_filters)
    c.ROW_LIMIT = 100000
    try:
        t = c.query_region(center,radius=radius).values()[0]
    except:
        raise IOError("Error while querying the given coords. You might not have an internet connection")
    
    cat = MASSCatalogue(empty=True)
    cat.create(t.columns,None,
               key_class="PointSource",value_star=None,
               key_ra="RAJ2000",key_dec="DEJ2000")
    return cat

# ------------------- #
# - 2MASS CATALOGUE - #
# ------------------- #
class MASSCatalogue( Catalogue ):
    """
    """
    source_name = "2MASS"
    
    def __init__(self, catalogue_file=None,empty=False,
                 key_mag=None,key_magerr=None,key_ra=None,key_dec=None):
        """
        """
        self.__build__(data_index=2,key_mag=key_mag,
                       key_magerr=key_magerr,
                       key_ra=key_ra,key_dec=key_dec)
        if empty:
            return
        
        self.load(catalogue_file)        
    
    @_autogen_docstring_inheritance(Catalogue.set_mag_keys,"Catalogue.set_mag_keys")
    def set_mag_keys(self,key_mag,key_magerr):
        #
        # add lbda def
        #
        super(MASSCatalogue,self).set_mag_keys(key_mag,key_magerr)
        if key_mag is not None:
            self.lbda = "TO BE DEFINED"

    # ----------------------- #
    # -  CATALOGUE HACK     - #
    # ----------------------- #
    @property
    def mag(self):
        if not self._is_keymag_set_(verbose=False):
            print "No 'key_mag' defined. J band used by default. -> To change: set_mag_keys() "
            self.set_mag_keys("Jmag","e_Jmag")
            
        return super(MASSCatalogue,self).mag

    # ------------------------------
    # - All points are Point Sources
    @property
    def _objecttype(self):
        print "All Loaded data are %s"%self._build_properties["key_class"]
        return np.ones(self.nobjects)

    @property
    def starmask(self):
        """ This will tell which of the datapoints is a star
        Remark, you need to have defined key_class and value_star
        in the __build_properties to be able to have access to this mask
        ==> In 2MASS PointSource catalogue, all data are stars
        """
        return np.ones(self.nobjects_in_fov,dtype="bool")  #not self.fovmask already in objecttype


#################################
#                               #
# BASIC WISE: Catalogue         #
#                               #
#################################
def fetch_wise_catalogue(center,radius,extracolums=[],column_filters={"Jmag":"10..30"}):
    """
    """
    try:
        from astroquery import vizier
    except:
        raise ImportError("install astroquery. (pip install astroquery)")
    
    # -----------
    # - DL info
    columns = ["AllWISE","ID",
               "RAJ2000","DEJ2000",
               ]
        
    for band in ["J","H","K","W1","W2","W3","W4"]:
        columns.append("%smag"%band)
        columns.append("e_%smag"%band)
    
    columns = columns+extracolums
    # - WARNING if discovered that some of the bandmag were missing if too many colums requested
    c = vizier.Vizier(catalog="II/328", columns=columns, column_filters=column_filters)
    c.ROW_LIMIT = 100000
    try:
        t = c.query_region(center,radius=radius).values()[0]
    except:
        raise IOError("Error while querying the given coords. You might not have an internet connection")
    
    cat = WISECatalogue(empty=True)
    cat.create(t.columns,None,
               key_class="ToBeDone",value_star=6,
               key_ra="RAJ2000",key_dec="DEJ2000")
    return cat

# ------------------- #
# - WISE CATALOGUE  - #
# ------------------- #
class WISECatalogue( Catalogue ):
    """
    """
    source_name = "WISE"
    
    def __init__(self, catalogue_file=None,empty=False,
                 key_mag=None,key_magerr=None,key_ra=None,key_dec=None):
        """
        """
        
        print "STAR vs. GALAXY PARSING NOT READY YET"
        
        self.__build__(data_index=2,key_mag=key_mag,
                       key_magerr=key_magerr,
                       key_ra=key_ra,key_dec=key_dec)
        if empty:
            return
        
        self.load(catalogue_file)        
    
    @_autogen_docstring_inheritance(Catalogue.set_mag_keys,"Catalogue.set_mag_keys")
    def set_mag_keys(self,key_mag,key_magerr):
        #
        # add lbda def
        #
        super(WISECatalogue,self).set_mag_keys(key_mag,key_magerr)
        if key_mag is not None:
            self.lbda = "TO BE DEFINED"
            
    @property
    def mag(self):
        if not self._is_keymag_set_(verbose=False):
            print "No 'key_mag' defined. W1 band used by default. -> To change: set_mag_keys() "
            self.set_mag_keys("W1mag","e_W1mag")
            
        return super(WISECatalogue,self).mag

    