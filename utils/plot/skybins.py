#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes for binning sky coordinates
"""

import warnings
import numpy as np

from ...astrobject.baseobject import BaseObject
#import astrobject.astrobject import baseobject as b

try:
    import healpy as hp
    HEALPY_IMPORTED = True
except ImportError:
    HEALPY_IMPORTED = False

_d2r = np.pi / 180 

__all__ = ["SkyBins","HealpixBins"]

class SkyBins( BaseObject ):
    """
    Object to collect all information for rectangular binning.
    
    # Currently requires that bins do not cross the RA = 180 deg line
    # to work properly. (No warnings issued). Will fix this later.
    """
    _properties_keys = ["ra_min", "ra_max", "dec_min", "dec_max",
                        "max_stepsize"]

    _derived_keys = ["nbins"]

    def __init__(self, ra_range=(-180, 180), dec_range=(-90, 90), 
                 ra_nbins=12, dec_nbins=6, dec_sin=True, empty=False,
                 max_stepsize=5):
        """
        """
        self.__build__()
        if empty:
            return
        self.create(ra_range=ra_range, dec_range=dec_range, 
                    ra_nbins=ra_nbins, dec_nbins=dec_nbins,dec_sin=dec_sin,
                    max_stepsize=max_stepsize)
        

    def create(self, ra_range=(-180, 180), dec_range=(-90, 90), 
               ra_nbins=18, dec_nbins=10, dec_sin=True, max_stepsize=5):
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

        self._properties["ra_min"] = ra_min.flatten() 
        self._properties["ra_max"] = ra_max.flatten() 
        self._properties["dec_min"] = dec_min.flatten() 
        self._properties["dec_max"] = dec_max.flatten() 
                         
        self._properties["max_stepsize"] = max_stepsize
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

    def hist(self, ra, dec, verbose=True):
        """
        Return the counts per bin for a list a coordinates.
        If verbose, it will issue warning the number of coordinates outside the bins
        and notify the user if coordinates where in multiple bins.
        """
        outside = 0
        binned = np.zeros(self.nbins)
        
        for r, d in zip(ra, dec):
            k = self.coord2bin(r, d)
            if np.isnan(k[0]):
                outside += 1
                continue
            elif len(k) > 1 and verbose:
                warnings.warn("Some bins appear to be overlapping.")
            binned[k] += 1

        if outside > 0 and verbose:
            warnings.warn("%i points lay outside the binned area.")

        return binned
        
    # ========================== #
    # = Utilities for plotting = #
    # ========================== #
    def boundary(self, k, steps=None, max_stepsize=None):
        """
        Return boundary of a bin; used for drawing polygons.
        If steps is None, max_stepsize is used to automatically determine 
        the appropriate step size.
        """
        if max_stepsize is None:
            max_stepsize = self.max_stepsize
        lkwargs = dict(steps=steps, max_stepsize=max_stepsize)
        
        # Start in top left and work clockwise 
        ra1, dec1 = self._draw_line(self.ra_min[k], self.dec_max[k], 
                                    self.ra_max[k], self.dec_max[k], **lkwargs)
        ra2, dec2 = self._draw_line(self.ra_max[k], self.dec_max[k], 
                                    self.ra_max[k], self.dec_min[k], **lkwargs)
        ra3, dec3 = self._draw_line(self.ra_max[k], self.dec_min[k], 
                                    self.ra_min[k], self.dec_min[k], **lkwargs)
        ra4, dec4 = self._draw_line(self.ra_min[k], self.dec_min[k], 
                                    self.ra_min[k], self.dec_max[k], **lkwargs)
        
        ra = np.concatenate([ra1, ra2[1:], ra3[1:], ra4[1:-1]])
        dec = np.concatenate([dec1, dec2[1:], dec3[1:], dec4[1:-1]])
        
        return ra, dec
        
    def _draw_line(self, ra1, dec1, ra2, dec2, steps=None, 
                   max_stepsize=None):
        """
        Return 'line' between (ra1, dec1) and (ra2, dec2).
        If steps is None, max_stepsize is used to automatically determine 
        the appropriate step size.
        Note: This treats coordinates as Euclidean in (ra, dec) space.
        Therefore it only works for lines of constant ra or dec. 
        """
        if max_stepsize is None:
            max_stepsize = self.max_stepsize
        
        if ra1 == ra2:
            if steps is None:
                steps = self._determine_steps(dec1, dec2, 
                                              max_stepsize=max_stepsize)
            ra = ra1 * np.ones(steps)
            dec = np.linspace(dec1, dec2, steps)
        elif dec1 == dec2:
            if steps is None:
                steps = self._determine_steps(ra1, ra2, 
                                              max_stepsize=max_stepsize)
            ra = np.linspace(ra1, ra2, steps)
            dec = dec1 * np.ones(steps)
        else:
            raise ValueError("Either ra1 and ra2 or dec1 and dec2 must be the same.")
        
        return ra, dec

    def _determine_steps(self, x1, x2, max_stepsize=None):
        """
        """
        if max_stepsize is None:
            max_stepsize = self.max_stepsize

        return np.ceil(np.abs(float(x2 - x1)) / max_stepsize) + 1

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

    @property
    def max_stepsize(self):
        """maximum stepsize for boundary in degrees"""
        return self._properties["max_stepsize"]


    # --------------------
    # - Derived Properties
    @property
    def nbins(self):
        """number of bins"""
        return len(self._properties["ra_min"])
    


class HealpixBins( BaseObject ):
    """
    HEALPix Binner
    """
    _properties_keys = ["nside", "nest", "max_stepsize"]

    _derived_keys = ["nbins"] 

    def __init__(self, nside=8, nest=True, max_stepsize=5, empty=False):
        """
        """
        self.__build__()
        if empty:
            return
        self.create(nside=nside, nest=nest, max_stepsize=max_stepsize)

    def create(self, nside=8, nest=True, max_stepsize=5):
        """
        """
        if not HEALPY_IMPORTED:
            raise ValueError("HealpixBins requires the healpy package.")

        self._properties["nside"] = nside
        self._properties["nest"] = nest
                         
        self._properties["max_stepsize"] = max_stepsize
        #self.__update__()

    # =========== #
    # = Binning = #
    # =========== #
    def hist(self, ra, dec):
        """
        Return the counts per bin for a list a coordinates.
        """
        pixels = hp.ang2pix(self.nside,(90 - dec) * _d2r, ra * _d2r, 
                            nest=self.nest)
        binned = np.histogram(pixels, bins=range(hp.nside2npix(self.nside) + 1))
        
        return binned[0]

    # ========================== #
    # = Utilities for plotting = #
    # ========================== #
    def boundary(self, k, steps=None, max_stepsize=None):
        """
        Return boundary of a bin; used for drawing polygons.
        If steps is None, max_stepsize is used to automatically determine
        the appropriate step size.
        """
        if max_stepsize is None:
            max_stepsize = self.max_stepsize

        if steps is None:
            steps = self._determine_steps(max_stepsize=max_stepsize)

        bd = hp.boundaries(self.nside, k, step=steps, nest=self.nest)
        dec_raw, ra_raw = hp.vec2ang(np.transpose(bd))

        ra = ((ra_raw / _d2r + 180) % 360) - 180
        dec = 90 - dec_raw / _d2r

        return ra, dec

    def _determine_steps(self, max_stepsize=None):
        """
        """
        if max_stepsize is None:
            max_stepsize = self.max_stepsize

        return np.ceil(hp.nside2resol(self.nside) / _d2r / max_stepsize) + 1

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def nside(self):
        """HEALPix NSIDE parameter"""
        return self._properties["nside"]

    @property
    def nest(self):
        """HEALPix NESTED parameter"""
        return self._properties["nest"]

    @property
    def max_stepsize(self):
        """maximum stepsize for boundary in degrees"""
        return self._properties["max_stepsize"]

    # --------------------
    # - Derived Properties
    @property
    def nbins(self):
        """number of bins"""
        return hp.nside(self._properties["nside"])
