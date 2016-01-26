#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains basic classes related to the surveys functions"""

import numpy as np
from .baseobject import BaseObject


class Pointing( BaseObject ):
    """
    This simple instance contains the basic information about the pointing
    """
    _propoerties_keys = ["fov","ra","dec","rotation"]

    def __init__(self,fov=None,
                 ra=None,dec=None,rotation=None,
                 empty = False,**kwargs):
        """
        Parameters
        ----------

        fov: [float]               
        """

    def contains_target(self,ra,dec):
        """
        """

    def move(self,ra,dec,rotation=None):
        """
        This function enables to move the coordinate of the center of the fov
        DEV REMARK: The rotation should derived from the coords....
        """

    def display(self,ax,**kwargs):
        """Draw on the given axes the object's fov"""
        
        

    
class Cadence( BaseObject ):
    """This instance manage the instruments bands, pointing and timing"""

    __nature__ = "Cadence"
    
    _properties_keys = ["npointings"]

    def __init__(self,npointings=1000,empty=False):
        """
        """
        print "TEMPORARY TEST CADENCE SOFTWARE"
        self.__build__()
        if empty:
            return
        
        self.npointings = npointings

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def npointings(self):
        return self._properties["npointings"]
    
    @npointings.setter
    def npointings(self, value):
        if value <= 0:
            raise ValueError("npointings cannot be lower or equal to zero")
        self._properties["npointings"] = value
        
    # -----------------
    # - TMP Skynoises
    @property
    def skynoises(self):
        """
        """
        sigma_PSF = 0.5 #is the standard deviation of the PSF in pixels 
        sigma_pixel = 100  #is the background noise in a single pixel in counts
        # see http://sncosmo.readthedocs.org/en/v1.1.x/_modules/sncosmo/simulation.html
        return np.ones(self.npointings)* (4 * pi * sigma_PSF * sigma_pixel)

    @property
    def mjds(self):
        """
        """
        return np.linspace(57700.0,58900.0,self.npointings) + np.random.rand(self.npointings)*50

    @property
    def bands(self):
        """
        """
        default_instruments = ["desg","desr"]
        return [default_instruments[int(i)] for i in np.round(np.random.rand(1000))]
    
