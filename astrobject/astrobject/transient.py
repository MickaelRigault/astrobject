#! /usr/bin/env python
# -*- coding: utf-8 -*-

from .baseobject import AstroTarget

__all__     = ["transient"]

def snIa(*args,**kwargs):
    print "not ready yet"
    return Supernova(*args,**kwargs)

def transient(lightcurve=None,name=None,
              zcmb=None, ra=None, dec=None,
              type_=None,forced_mwebmv=None,
              mjd=None,
              datasource={},**kwargs):
    """
    = Initialize the AstroTarget Function =

    Parameters
    ----------
       ** You have 2 way to initialize the object **
       1) By giving the following parameters
       2) By setting datasource as dict having these parameters

    ..Method 1:

    lightcurve: [LightCurve]   The lightcurve associated to the transient
    
    name: [string]             name of the astro-object

    zcmb: [float]              redshift in the cmb frame. This will be used
                               to derive distmeter etc. depending on the cosmology

    ra: [float]                right-ascention of the object. (degree favored).
                               *ra* and *dec* must have the same unit.
                                   
    dec: [float]               declination of the object. (degree favored)
                               *ra* and *dec* must have the same unit.

    mjd: [float/array]         typical time(s) for the transient (could be a list)
                               If None and a lightcurve is set. The first time of the
                               LightCurve point will be used.
    - options -
        
    type_:[string]             type of the astro-object (galaxy, sn, Ia, Ibc...)
                               (no predifined list type so far, but this could append)

    forced_mwebmv: [float]     Force the Milky way extinction for this object.
                               Otherwise this
                               extinction depend on the object *radec*.
                               ->Use this use Caution<-

    ..Method 2:

    datasource: [dict]         Dictionnary having {name,ra,dec,...} values in the
                               above parameter will be used if key not found.

    - kwargs options, potentially non-exhaustive ; kwargs goes to AstroTarget -
        
    cosmo:[astropy.cosmology]  the cosmology used to derive the distances and so on.

    load_from: [dict]          a dictionary having an entry for each of the fundamental
                               parameter (name, ra, dec, zcmb...)

    empty: [bool]              Does not do anything, just loads an empty object.
                              (Careful with that)
                              
    Code
    ----
    return AstroTarget(name=datasource.pop("name",name),
                       zcmb=datasource.pop("zcmb",zcmb),
                       ra=datasource.pop("ra",ra),
                       dec=datasource.pop("dec",dec),
                       type_=datasource.pop("type",type_),
                       forced_mwebmv=datasource.pop("forced_mwebmv",forced_mwebmv),
                       **kwargs).copy()

                       
    Returns
    -------
    AstroTarget
    """
    if "name" in datasource.keys() and "object" not in datasource.keys():
        datasource["object"] = datasource["name"]

    # datasource is a bit overkill since **datasource would have done the same job
    return Transient( lightcurve=datasource.pop("lightcurve",lightcurve),
                      name=datasource.pop("object",name),
                      zcmb=datasource.pop("zcmb",zcmb),
                      ra=datasource.pop("ra",ra),
                      dec=datasource.pop("dec",dec),
                      type_=datasource.pop("type",type_),
                      forced_mwebmv=datasource.pop("forced_mwebmv",forced_mwebmv),
                      mjd=datasource.pop("mjd",mjd),
                      **kwargs).copy() # dont forget the copy
#######################################
#                                     #
# Base Object Classes: Transient      #
#                                     #
#######################################
class Transient( AstroTarget ):
    """This instance enables time information in the astroTarget
    it has a lightcurve"""

    PROPERTIES         = ['lightcurve']
    SIDE_PROPERTIES    = ['mjd']
    DERIVED_PROPERTIES = ["fits","lbda","raw_lbda"]
    
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,lightcurve=None,
                 name=None,zcmb=None,
                 ra=None,dec=None,type_=None,
                 mjd=None,empty=False,**kwargs):
        """ """
        super(Transient,self).__init__(name=name,zcmb=zcmb,
                            ra=ra,dec=dec,type_=type_,
                            empty=empty,**kwargs)
        
        self.set_lightcurve(lightcurve)
        self.mjd = mjd
        


    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def set_lightcurve(self,lightcurve,force_it=False):
        """
        """
        if self.has_lightcurve() and force_it is False:
            raise AssertionError("A lightcurve is already defined. "+\
                                 " Set force_it to true to overwrite it.")
        if lightcurve is None:
            self._properties["lightcurve"] = None
            
        elif "__nature__" not in dir(lightcurve) or \
          lightcurve.__nature__ != "LightCurve":
          raise TypeError("'lightcurve' must be a astrobject LightCurve, ",lightcurve)
      
        

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def lightcurve(self):
        return self._properties['lightcurve']
    
    def has_lightcurve(self):
        return False if self.lightcurve is None \
          else True
          
    @property
    def mjd(self):
        """typical time associated to the transient"""
        if self.has_lightcurve() and self._side_properties["mjd"] is None:
            return self.lightcurve.times[0]
        
        return self._side_properties["mjd"]
    
    @mjd.setter
    def mjd(self,value):
        self._side_properties["mjd"] = value
#######################################
#                                     #
# Base Object Classes: Supernova      #
#                                     #
#######################################  
class Supernova( Transient ):
    """This object is a supernova, it has an explosion date and a
    lightcurve"""
