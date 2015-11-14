#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module defines the basic low-level objects"""


import numpy as np
import copy

from ..utils.tools import load_pkl,dump_pkl

from astropy import cosmology
from astropy import units
COSMO_DEFAULT = cosmology.Planck13


__version__ = 0.1
__all__     = ["astrobject"]


def astrobject(name,zcmb,ra,dec,
               type_=None,forced_mwebmv = None,**kwargs):
    """
    = Initialize the AstroObject Function =

    Parameters
    ----------
        
    name: [string]             name of the astro-object

    zcmb: [float]              redshift in the cmb frame. This will be used
                               to derive distmeter etc. depending on the cosmology

    ra: [float]                right-ascention of the object. (degree favored).
                               *ra* and *dec* must have the same unit.
                                   
    dec: [float]               declination of the object. (degree favored)
                               *ra* and *dec* must have the same unit.
                                   
        - options -
        
    type_:[string]             type of the astro-object (galaxy, sn, Ia, Ibc...)
                               (no predifined list type so far, but this could append)

    forced_mwebmv: [float]     Force the Milky way extinction for this object.
                               Otherwise this
                               extinction depend on the object *radec*.
                               ->Use this use Caution<-
                                   
        - kwargs options, potentially non-exhaustive -
        
    cosmo:[astropy.cosmology]  the cosmology used to derive the distances and so on.

    load_from: [dict]          a dictionary having an entry for each of the fundamental
                               parameter (name, ra, dec, zcmb...)

    empty: [bool]              Does not do anything, just loads an empty object.
                              (Careful with that)
                                   
    Returns
    -------
    AstroObject
    """
    return AstroObject(name=name,zcmb=zcmb,
                       ra=ra, dec=dec,type_=type_,
                       forced_mwebmv=forced_mwebmv,
                       **kwargs).copy() # dont forget the copy
    


#######################################
#                                     #
# Base Object Classes                 #
#                                     #
#######################################

class BaseObject( object ):
    """The zero Level Structure that has the basic tools"""
    _properties_keys         = []
    _side_properties_keys    = []
    _derived_properties_keys = []

    def __init__(self):
        self.__build__()
        
    def __build__(self):
        """Create the properties dictionary"""
        self._properties = {}
        for k in self._properties_keys:
            self._properties[k] = None
            
        self._side_properties = {}
        for k in self._side_properties_keys:
            self._side_properties[k] = None

        self._derived_properties = {}
        for k in self._derived_properties_keys:
            self._derived_properties[k] = None
            
    
    def copy(self):
        """returns an independent copy of the current object"""
        newobject = copy.deepcopy(self)
        # -- This to be 100% sure
        newobject._properties         = copy.deepcopy(self._properties)
        newobject._side_properties    = copy.deepcopy(self._side_properties)
        newobject._derived_properties = copy.deepcopy(self._derived_properties)
        newobject._update_()
        
        return newobject
        
                
    def _update_(self):
        """Adapte the derived properties as a function of the main ones"""
        return

    # ================ #
    # = Properties   = #
    # ================ #
    @property
    def _fundamental_parameters(self):
        return self._properties.keys()
    
    
class AstroObject( BaseObject ):
    """
    This is the default astrophysical object that has basic
    information attached to it, like zcmb, ra, dec, mwebmv.

    You will also have access to derived parameter like the distances
    (meter and mpc) or the angular/physical size ratio based on the
    given cosmology and redshift (zcmb).
    """
    # -- This is set to ease inheritance and tests -- #
    # If you change that, function that needs astrobject
    # could crash (i.e. astroimages.Aperture)
     
    __nature__ = "AstroObject"
    
    # -------------------- #
    # Internal Properties  #
    # -------------------- #
    _properties_keys         = ["zcmb","ra","dec","name"]
    _side_properties_keys    = ["cosmology","litrature_name",
                                "type","mwebmv"]
    _derived_properties_keys = ["distmeter","distmpc",
                                "arc_per_kpc"]
    

    # =========================== #
    # = Initialization          = #
    # =========================== #
    def __init__(self,name=None,zcmb=None,
                 ra=None,dec=None,
                 type_=None,cosmo=COSMO_DEFAULT,
                 load_from=None,empty=False,**kwargs
                 ):
        """
         = Initialize the AstroObject Function =
        
        Parameters
        ----------
        
        name: [string]             name of the astro-object

        zcmb: [float]              redshift in the cmb frame. This will be used
                                   to derive distmeter etc. depending on the cosmology

        ra: [float]                right-ascention of the object. (degree favored).
                                   *ra* and *dec* must have the same unit.
                                   
        dec: [float]               declination of the object. (degree favored)
                                   *ra* and *dec* must have the same unit.

        type_:[string]             type of the astro-object (galaxy, sn, Ia, Ibc...)
                                   (no predifined list type so far. It could append)

        cosmo:[astropy.cosmology]  the cosmology used to derive the distances etc.

        load_from: [dict]          a dictionary having an entry for each of the
                                   fundamental
                                   parameter (name, ra, dec, zcmb...)

        empty: [bool]              Does not do anything, just loads an empty object.
                                   (Careful with that)
        Returns
        -------
        Void
        """
        self.__build__()
        
        if empty:
            return
        
        if load_from is not None:
            self.load(load_from,**kwargs)
            return

        self.define(name,zcmb,ra,dec,
                    cosmo=cosmo,type_=type_,
                    **kwargs)
        return
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def define(self,name,zcmb,ra,dec,
               cosmo=COSMO_DEFAULT,type_=None,
               forced_mwebmv=None):
        """
        This function enables you to define the fundamental object parameter
        upon which the entire object will be defined.
        """
        self.name = name
        self.zcmb = zcmb
        self.radec = ra,dec
        self.type  = type_
        self.set_cosmo(COSMO_DEFAULT)
        self.set_mwebmv(forced_mwebmv,force_it=True)
        self._update_()

    def define_from_properties(self,properties,
                               side_properties={}):
        """
        """
        pass
        
    def writeto(self,output_file,**kwargs):
        """
        This function save the basic information of the object
        in the given output_file. options goes to dump_pkl
        """
        dump_pkl(self._properties,output_file,**kwargs)
        return

    def load(self,properties,side_properties=None):
        """
        """
        if type(properties) == str:
            candidatepkl = load_pkl(properties)
        elif type(properties) == dict:
            candidatepkl = properties
            
        for k in self._fundamental_parameters:
            if k not in candidatepkl.keys():
                raise AttributeError("'%s' is a requested fundamental "+\
                                     "parameter and is not in the input dictionnary"%k)

        self.define(candidatepkl["name"],candidatepkl["zcmb"],
                    candidatepkl["ra"],candidatepkl["dec"],
                    cosmo = candidatepkl.pop("cosmo",COSMO_DEFAULT),
                    type_ = candidatepkl.pop("type",None),
                    forced_mwebmv = candidatepkl.pop("mwebmv",None))
        
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    # ------------------ #
    # - Object Name    - #
    # ------------------ #
    @property
    def name(self):
        return self._properties["name"]
    
    @property
    def _litrature_name(self):
        return self._side_properties["litrature_name"]
    
    @name.setter
    def name(self,value):
        self._properties["name"] = value
        self._check_litrature_name_()
                          
    # ------------------ #
    # - Redshift       - #
    # ------------------ #
    @property
    def zcmb(self):
        return self._properties["zcmb"]
    
    @zcmb.setter
    def zcmb(self,value):
        self._properties["zcmb"] = value
        self._update_distance_()
    
    # ------------------ #
    # - Coordinate     - #
    # ------------------ #
    @property
    def ra(self):
        return self._properties["ra"]
    @property
    def dec(self):
        return self._properties["dec"]
    @property
    def radec(self):
        return self.ra,self.dec
    
    @radec.setter
    def radec(self,value):
        if np.shape(value) != (2,):
            raise SyntaxError("radec must have two input parameters: ra,dec")
        
        self._properties["ra"],self._properties["dec"] = value
        self._update_mwebmv_()

    # ------------------ #
    # - MW Extinction  - #
    # ------------------ #
    @property
    def mwebmv(self):
        return self._side_properties["mwebmv"]
    
    def set_mwebmv(self,value,force_it=False):
        if force_it is False:
            raise AssertionError("You should not manually change mwebmv, "+\
                                 "it is bound to the object coordinate."+\
                    "Set force_it to True if you really know what you are doing.")
        self._side_properties["mwebmv"] = value
    
    # ------------------ #
    # - SN Type        - #
    # ------------------ #
    @property
    def type(self):
        return self._side_properties["type"]

    @type.setter
    def type(self,value):
        self._side_properties["type"] = value

        
    # ========================= #
    # = Derived Values        = #
    # ========================= #
    # ------------------ #
    # - Distances      - #
    # ------------------ #
    # If you wnat to change these, change
    # the cosmology or the redshift   
    @property
    def distmeter(self):
        return self._derived_properties["distmeter"]
    @property
    def distmpc(self):
        return self._derived_properties["distmpc"]
    @property
    def arcsec_per_kpc(self):
        return self._derived_properties["arcsec_per_kpc"]
    
    # -------------------- #
    # - COSMOLOGY        - #
    # -------------------- #
    @property
    def cosmo(self):
        return self._side_properties["cosmology"]
    
    def set_cosmo(self,astropycosmo):
        """change the object according to the given cosmology.
        this cosmology must be an astropy one"""
        if "astropy" not in astropycosmo.__module__:
            raise ValueError("'astropycosmo' must be an astropy cosmology object")
        
        self._side_properties["cosmology"] = astropycosmo
        self._update_distance_()
        
    # ========================= #
    # = Internal Tools        = #
    # ========================= #
    def _update_(self):
        self._update_mwebmv_()
        self._update_distance_()
        self._check_litrature_name_()
        
    def _update_mwebmv_(self,verbose=False):
        if verbose:print "_update_mwebmv_ to be done"
            
    def _check_litrature_name_(self,verbose=False):
        if verbose:print "_check_litrature_name_ to be done"
        if self.name is None:
            self._side_properties["litrature_name"] = None
        
    def _update_distance_(self):
        if self.cosmo is None:
            return
        if self.zcmb is None:
            return

        self._derived_properties["distmpc"]   = \
          self.cosmo.luminosity_distance(self.zcmb).value
          
        self._derived_properties["distmeter"] = \
          self._derived_properties["distmpc"] * units.mpc.in_units("m")
          
        self._derived_properties["arcsec_per_kpc"] = \
          self.cosmo.arcsec_per_kpc_proper(self.zcmb).value
          