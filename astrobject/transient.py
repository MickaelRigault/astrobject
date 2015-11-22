#! /usr/bin/env python
# -*- coding: utf-8 -*-

from .baseobject import AstroTarget

__all__     = ["snIa"]

#######################################
#                                     #
# Base Object Classes: Transient      #
#                                     #
#######################################
class Transient( AstroTarget ):
    """This instance enables time information in the astroTarget
    it has a lightcurve"""

    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,lightcurve=None,
                 name=None,zcmb=None,
                 ra=None,dec=None,type_=None,
                 empty=False,**kwargs):
        
        super(Transient,self).__init__(name=name,zcmb=zcmb,
                            ra=ra,dec=dec,type_=type_,
                            empty=empty,**kwargs)
        
        self.set_lightcurve(lightcurve)

        
    def __build__(self):
        self._properties_keys.append('lightcurve')
        super(Transient,self).__build__()

    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def set_lightcurve(self,lightcurve,force_it=False):
        """
        """
        if "__nature__" not in dir(lightcurve) or \
          lightcurve.__nature__ != "LightCurve":
          raise TypeError("'lightcurve' must be a astrobject LightCurve")
      
        if self.has_lightcurve() and force_it is False:
            raise AssertionError("A lightcurve is already defined. "+\
                                 " Set force_it to true to overwrite it.")

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def lightcurve(self):
        return self._properties['lightcurve']
    
    def has_lightcurve(self):
        return False if self.lightcurve is None \
          else True

#######################################
#                                     #
# Base Object Classes: Supernova      #
#                                     #
#######################################  
class Supernova( Transient )
    """This object is a supernova, it has an explosion date and a
    lightcurve"""
