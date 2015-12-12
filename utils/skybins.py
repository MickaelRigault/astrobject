#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes for binning sky coordinates
"""

import numpy as np

from ..astrobject.baseobject import BaseObject,astrotarget

_d2r = np.pi / 180 

__all__ = ["SkyBins","HealpixBins"]

class SkyBins( BaseObject ):
    """
    """
    _properties_keys = ["ra_min", "ra_max", "dec_min", "dec_max"]

    _derived_keys = ["nbins"]

    def __init__(self, ra_range=(-180, 180), dec_range=(-90, 90), 
                 ra_nbins=18, dec_nbins=10, dec_sin=True, empty=False):
        """
        """
        self.__build__()
        if empty:
            return
        self.create(ra_range=ra_range, dec_range=dec_range, 
                    ra_nbins=ra_nbins, dec_nbins=dec_nbins,dec_sin=dec_sin)
        

    def create(self, ra_range=(-180, 180), dec_range=(-90, 90), 
               ra_nbins=18, dec_nbins=10, dec_sin=True):
        """
        """
        ra_base = np.linspace(ra_range[0], ra_range[1], ra_nbins+1)
        if dec_sin:
            dec_base = np.arcsin(np.linspace(np.sin(dec_range[0] * _d2r), 
                                             np.sin(dec_range[1] * _d2r), 
                                             dec_nbins+1)) / _d2r
        else:
            dec_base = np.linspace(dec_range[0], dec_range[1], dec_nbins+1)

        ra_min, dec_min = np.meshgrid(ra_base[:-1], dec_base[:-1])
        ra_max, dec_max = np.meshgrid(ra_base[1:], dec_base[1:])

        self._properties['ra_min'] = ra_min.flatten() 
        self._properties['ra_max'] = ra_min.flatten() 
        self._properties['dec_min'] = dec_min.flatten() 
        self._properties['dec_max'] = dec_max.flatten() 
        
        #self.__update__()

    # =========== #
    # = Binning = #
    # =========== #
    def coord2bin(self, ra, dec):
        """
        """
        k = np.where((ra > self.ra_min) &
                     (ra <= self.ra_max) &
                     (dec > self.dec_min) &
                     (dec <= self.dec_max))[0]

        if len(k) == 0:
            return np.asarray([np.NaN])
        else:
            return k

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def ra_min(self):
        """list of minimum ra of each bin"""
        return self._properties["ra_min"]

    @property
    def ra_max(self):
        """list of maximum ra of each bin"""
        return self._properties["ra_max"]

    @property
    def dec_min(self):
        """list of minimum dec of each bin"""
        return self._properties["dec_min"]

    @property
    def dec_max(self):
        """list of maximum dec of each bin"""
        return self._properties["dec_max"]

    # --------------------
    # - Derived Properties
    @property
    def nbins(self):
        """number of bins"""
        return len(self._properties["ra_min"])
    


class HealpixBins( BaseObject ):
    """
    """
    _properties_keys = ["ra_min", "ra_max", "dec_min", "dec_max"]

    _derived_keys = ["nbins"]
