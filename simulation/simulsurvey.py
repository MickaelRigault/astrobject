#! /usr/bin/env python
# -*- coding: utf-8 -*-

from ..astrobject.baseobject import BaseObject
from .simultarget import transient_generator


class SimulSurvey( BaseObject ):
    """
    Basic survey object
    (far from finished)
    """
    _properties_keys         = ["targetlist","instprop"]
    _side_properties_keys    = []
    _derived_properties_keys = []
    
    def __init__(self):
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

    @property
    def lat(self):
        return self._properties["lat"]
            
