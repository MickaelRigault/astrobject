#! /usr/bin/env python
# -*- coding: utf-8 -*-

from ..astrobject.baseobject import BaseObject

class Survey( BaseObject ):
    """
    Basic survey object
    (far from finished)
    """
    _properties_keys         = ["lon","lat","exptime"]
    _side_properties_keys    = []
    _derived_properties_keys = ["maglimit"]
    
    def __init__(self, filters=None,empty=None):
        """
        """
        self.__build__()
        self.lon = 30
        if empty:
            return 
        #self.filters = filters

    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def set_targets(self, astrotargets):
        """
        astrotarget have ra,dec, zcmb ... distance
        supernova in additation has a abs-mag, light-curve
                                 => obs mag
        """
        
    def recover_targets(self):
        """
        bunch threshold...
        """

    def recover_lightcurves(self):
        """
        """
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def lon(self):
        return self._properties["lon"]

    @lon.setter
    def lon(self,value):
        if value < 20:
            raise ValueError("should be greater than 20 (%f)"%value)
        self._properties["lon"] = value
            
