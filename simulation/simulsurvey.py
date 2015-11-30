#! /usr/bin/env python
# -*- coding: utf-8 -*-

from ..astrobject.baseobject import BaseObject
from .simultarget import transient_generator


#######################################
#                                     #
# Survey: Simulation Base             #
#                                     #
#######################################
class SimulSurvey( BaseObject ):
    """
    Basic survey object
    (far from finished)
    """
    _properties_keys         = ["targetlist","instprop"]
    _side_properties_keys    = []
    _derived_properties_keys = []
    
    def __init__(self,targetlist=None,
                 instprop=None,
                 empty=False):
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
    def set_targets(self, targetlist):
        """
        astrotarget have ra,dec, zcmb ... distance
        supernova in additation has a abs-mag, light-curve
                                 => obs mag
        """
        

    def set_instrument(self,bias=None,bands=None):
        """
        """
        if self.instprop is None:
            self.instprop = {}
            
        if bias is not None:
            
            self.instprop['bias'] = bias
        
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
            
    @property
    def bands(self):
        if self.instprop is None:
            raise AttributeError("no 'instprop' defined")

    @property
    def instprop(self):
        """This dictionary contains the basic instrument properties"""
        return self._properties['instprop']

    def is_instrument_set(self):
        """This function test if the given instrument is ok to be used"""
        if self.instprop is None or self.instprop.keys() == 0:
            raise AttributeError("'instprop' is None or empty")
        
        
