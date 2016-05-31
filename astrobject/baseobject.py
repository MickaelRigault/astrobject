#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module defines the basic low-level objects"""


import numpy as np
import copy
import warnings

from astropy import units
from .utils.tools import load_pkl,dump_pkl



__version__ = 0.1
__all__     = ["BaseObject","get_target"]

def get_target(name=None,zcmb=None, ra=None, dec=None,
               type_=None,mwebmv=None,zcmberr=None,**kwargs):
    """
    = Initialize the AstroTarget Function =

    Parameters
    ----------
        
    name: [string]             name of the astro-object

    zcmb,zcmberr: [float]      redshift/error in the cmb frame. This will be used
                               to derive distmeter etc. depending on the cosmology

    ra: [float]                right-ascention of the object. (degree favored).
                               *ra* and *dec* must have the same unit.
                                   
    dec: [float]               declination of the object. (degree favored)
                               *ra* and *dec* must have the same unit.
                               
    - options -
        
    type_:[string]             type of the astro-object (galaxy, sn, Ia, Ibc...)
                               (no predifined list type so far, but this could append)

    mwebmv (or forced_mwebmv): [float]
                               Force the Milky way extinction for this object.
                               Otherwise this
                               extinction depend on the object *radec*.
                               ->Use this use Caution<-

    - kwargs options, potentially non-exhaustive ; kwargs goes to AstroTarget -
        
    cosmo:[astropy.cosmology]  the cosmology used to derive the distances and so on.

    load_from: [dict]          a dictionary having an entry for each of the fundamental
                               parameter (name, ra, dec, zcmb...)

    empty: [bool]              Does not do anything, just loads an empty object.
                              (Careful with that)
                              
    Returns
    -------
    AstroTarget
    """
    if name is None and "object" in kwargs.keys():
        name = kwargs.pop("object")
        
    forced_mwebmv = kwargs.pop("forced_mwebmv",kwargs.pop("MWebmv",
                                             kwargs.pop("MWebv",None)))
    
    dec = kwargs.pop("Dec",dec)
    ra = kwargs.pop("Ra",ra)
    empty = kwargs.pop("empty",False)
    
    return AstroTarget(name=name,zcmb=zcmb,ra=ra,dec=dec,
                       type_=type_,
                       forced_mwebmv=forced_mwebmv,
                       empty=empty).copy()    





#######################################
#                                     #
# Base Object Classes                 #
#                                     #
#######################################    
class BaseObject( object ):
    """The zero level structure to handle the
    _proprerties, _side_properties, _derived_properties
    tricks.

    The Class inheriting from BaseObject could have the following
    global-class variables:

    ```
    PROPERTIES         = ["a","b"]
    SIDE_PROPERTIES    = ["t"]
    DERIVED_PROPERTIES = ["whynot"]
    ```
    if so, object created from this class or any inheriting class
    could have access to the self._properties["a"], self._properties["b"]
    self._side_properties["t"] or self._derived_properties["whynot"]
    parameters.

    BaseObject also have a basic .copy() method.
    """

    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    def __new__(cls,*arg,**kwargs):
        """ Upgrade of the New function to enable the
        _properties,_side_properties and _derived_properties tricks
        """
        obj = super(BaseObject,cls).__new__(cls)
        obj._properties_keys         = copy.copy([])
        obj._side_properties_keys    = copy.copy([])
        obj._derived_properties_keys = copy.copy([])
        # ---------------------------------------------
        # - load also the properties of all the Parents
        for c in obj.__class__.__mro__:
            if "PROPERTIES" in dir(c):
                obj._properties_keys += c.PROPERTIES
            if "SIDE_PROPERTIES" in dir(c):
                obj._side_properties_keys += c.SIDE_PROPERTIES
            if "DERIVED_PROPERTIES" in dir(c):
                obj._derived_properties_keys += c.DERIVED_PROPERTIES

        # -----------------
        # - avoid doublon
        obj._properties_keys         = np.unique(obj._properties_keys).tolist()
        obj._side_properties_keys    = np.unique(obj._side_properties_keys).tolist()
        obj._derived_properties_keys = np.unique(obj._derived_properties_keys).tolist()
        
        # -- keys
        if "_properties" not in dir(obj):
            obj._properties = {}
            obj._side_properties = {}
            obj._derived_properties = {}

        # -- fill empty
        for k in obj._properties_keys:
            if k in obj._properties.keys():
                warnings.warn("%s is already build.  Conflit => key ignored"%k)
                continue
            obj._properties[k] = None
            
        for k in obj._side_properties_keys:
            if k in obj._side_properties.keys():
                warnings.warn("%s is already build.  Conflit => key ignored"%k)
                continue
            obj._side_properties[k] = None
            
        for k in obj._derived_properties_keys:
            if k in obj._derived_properties.keys():
                warnings.warn("%s is already build.  Conflit => key ignored"%k)
                continue
            obj._derived_properties[k] = None

            
        return obj
    
    def __init__(self):
        self.__build__()
        
    def __build__(self):
        """Create the properties dictionary"""
        pass

    def copy(self, empty=False):
        """returns an independent copy of the current object"""
        
        # Create an empty object
        newobject = self.__class__(empty=True)
        if empty:
            return
        # And fill it !
        for prop in ["_properties","_side_properties",
                     "_derived_properties","_build_properties"
                     ]:
            if prop not in dir(self):
                continue
            try: # Try to deep copy because but some time it does not work (e.g. wcs) 
                newobject.__dict__[prop] = copy.deepcopy(self.__dict__[prop])
            except:
                newobject.__dict__[prop] = copy.copy(self.__dict__[prop])
        # This be sure things are correct
        newobject._update_()
        # and return it
        return newobject
        
                
    def _update_(self):
        """Adapte the derived properties as a function of the main ones"""
        pass

    # ================ #
    # = Properties   = #
    # ================ #
    @property
    def _fundamental_parameters(self):
        return self._properties.keys()
    
class AstroTarget( BaseObject ):
    """
    This is the default astrophysical object that has basic
    information attached to it, like zcmb, ra, dec, mwebmv.

    You will also have access to derived parameter like the distances
    (meter and mpc) or the angular/physical size ratio based on the
    given cosmology and redshift (zcmb).
    """
    # -- This is set to ease inheritance and tests -- #
    # If you change that, function that needs astrotarget
    # could crash (i.e. astroimages.Aperture)
     
    __nature__ = "AstroTarget"

    PROPERTIES         = ["zcmb","ra","dec","name","zcmb.err"]
    SIDE_PROPERTIES    = ["cosmology","literature_name","type","mwebmv","sfd98_dir"]
    DERIVED_PROPERTIES = ["distmeter","distmpc","arc_per_kpc"]

    # -------------------- #
    # Internal Properties  #
    # -------------------- #
    # =========================== #
    # = Initialization          = #
    # =========================== #
    def __init__(self,name=None,zcmb=None,zcmberr=None,
                 ra=None,dec=None,
                 type_=None,cosmo=None,
                 load_from=None,empty=False,sfd98_dir=None,
                 **kwargs):
        """
         = Initialize the AstroTarget Function =
        
        Parameters
        ----------
        
        name: [string]             name of the astro-object

        zcmb,zcmberr: [float]      redshift/error in the cmb frame. This will be used
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
        
        sfd98_dir: [string]        directory of SFD98 dust map files if not set
                                   in .astropy/config/sncosmo.cfg
                                   files SFD_[dust,mask]_4096_[n,s]gp.fits are
                                   available at http://sncosmo.github.io/data/dust/[file]

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
                    cosmo=cosmo,type_=type_,sfd98_dir=sfd98_dir,
                    **kwargs)
        return
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def define(self,name,zcmb,ra,dec,
               cosmo=None,type_=None,sfd98_dir=None,
               forced_mwebmv=None,zcmberr=None):
        """
        This function enables you to define the fundamental object parameter
        upon which the entire object will be defined.
        """
        self.name = name
        self.set_zcmb(zcmb,zcmberr)
        self.radec = ra,dec
        self.type  = type_
        if cosmo is None:
            from astropy.cosmology import Planck15
            cosmo = Planck15
            warnings.warn("Planck 2015 cosmology used by default")
            
        self.set_cosmo(cosmo)
        self._update_()
        if forced_mwebmv is not None:
            self.set_mwebmv(forced_mwebmv,force_it=True)
        
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

        cosmo = candidatepkl.pop("cosmo")
        self.define(candidatepkl["name"],candidatepkl["zcmb"],
                    candidatepkl["ra"],candidatepkl["dec"],
                    cosmo = cosmo,
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
    def _literature_name(self):
        return self._side_properties["literature_name"]
    
    @name.setter
    def name(self,value):
        self._properties["name"] = value
        self._check_literature_name_()
                          
    # ------------------ #
    # - Redshift       - #
    # ------------------ #
    @property
    def zcmb(self):
        return self._properties["zcmb"]
    
    def set_zcmb(self,value, error=None):
        self._properties["zcmb"] = value
        self._properties["zcmb.err"] = error
        self._update_distance_()

    @property
    def zcmb_err(self):
        return self._properties["zcmb.err"]
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
        self._side_properties["mwebmv"] = None

    # ------------------ #
    # - MW Extinction  - #
    # ------------------ #
    @property
    def mwebmv(self):
        if self._side_properties["mwebmv"] is None:
            self._update_mwebmv_()
        return self._side_properties["mwebmv"]
    
    def set_mwebmv(self,value,force_it=False):
        if force_it is False:
            raise AssertionError("You should not manually change mwebmv, "+\
                                 "it is bound to the object coordinate."+\
                    "Set force_it to True if you really know what you are doing.")
        self._side_properties["mwebmv"] = value
    
    @property
    def _sfd98_dir(self):
        """Director where the maps are. Default option set"""
        if self._side_properties["sfd98_dir"] is None:
            from .utils.io import get_default_sfd98_dir
            self._side_properties["sfd98_dir"] = get_default_sfd98_dir(download_if_needed=True)
            
        return self._side_properties["sfd98_dir"]

    def set_sfd98_dir(self, value):
        self._side_properties['sfd98_dir'] = value
        self._update_mwebmv_()

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
        """ converts the distmpc in meter """
        return self.distmpc * units.Mpc.in_units("m")

    @property
    def distmeter_err(self):
        """ converts the distmpc_err in meter """
        return self.distmpc_err * units.Mpc.in_units("m")
    
    @property
    def distmpc(self):
        if self.cosmo is None or self.zcmb is None:
            raise AttributeError("'cosmo' and 'redshift' required.")
        return self.cosmo.luminosity_distance(self.zcmb).value

    @property
    def distmpc_err(self):
        """ assymetric error to account for the non linearity of
        the distance measurement
        Returns:
        --------
        [-lower, +upper] (both in absolute value).
        (np.NaN,np.NaN are returned if no error on the redshift
        """
        if self.zcmb is None or self.zcmb_err is None:
            return np.NaN,np.NaN
        
        return np.asarray([
          self.cosmo.luminosity_distance(self.zcmb-self.zcmb_err).value-self.distmpc,
          self.distmpc-self.cosmo.luminosity_distance(self.zcmb+self.zcmb_err).value])
                
    
    
    @property
    def arcsec_per_kpc(self):
        if self.cosmo is None or self.zcmb is None:
            raise AttributeError("'cosmo' and 'redshift' requiered.")
        return self.cosmo.arcsec_per_kpc_proper(self.zcmb).value
    
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
        """
        Does not call _update_mewbmv_ because extinction should 
        only be fetched on demand.
        """
        self._update_distance_()
        self._check_literature_name_()
        
    def _update_mwebmv_(self):
        try:
            from sncosmo import SFD98Map
            m = SFD98Map(mapdir=self._sfd98_dir)
            self.set_mwebmv(m.get_ebv((self.ra, self.dec)), force_it=True)
        except IOError:
            warnings.warn("MW E(B-V) could not be fetched. Please set sfd98_dir to the map driectory.")
                
            
    def _check_literature_name_(self,verbose=False):
        if verbose:print "_check_literature_name_ to be done"
        if self.name is None:
            self._side_properties["literature_name"] = None
        
    def _update_distance_(self):
        """Do nothing so far since the distance are derived on the fly"""
        pass
          
