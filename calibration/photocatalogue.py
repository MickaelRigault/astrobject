#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module use the instrument objects to perform photometry"""
import numpy as np
from scipy.stats import sigmaclip
from astropy import units
from ..astrobject.baseobject import BaseObject
from ..astrobject.instruments.instrument import instrument 

class PhotoCatalogue( BaseObject ):
    """
    """
    __nature__ = "PhotoCatalogue"

    _properties_keys = ["instrument"]
    _side_properties_keys = ["catmag_range"]
    _derived_properties_keys = []

    
    def __init__(self,instrument=None,empty=True,**kwargs):
        """
        Parameters
        ----------

        instrument: [astrobject's instrument or filename]

        Return
        ------
        
        """
        self.__build__()

        if instrument is not None:
            if type(instrument) == str:
                self.load_instrument(instrument,**kwargs)
            else:
                self.create(instrument,**kwargs)

    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    # ---------------- #
    # - Set Function - #
    # ---------------- #
    # --------------- #
    # - ZeroPoints  - #
    # --------------- #
    def derive_zeropoint(self,catmag_range=None, clipping=3):
        """
        Parameters
        ----------
        catmag_range: [2D-array / None]
                                   The magnitude range (of the catalogue) used to
                                   derive the zero point.
                                   If None the current one (self.catmag_range) will
                                   be used. otherwise this will overwirte it.
                                   
        """
        self._test_instrument_()
        if "__iter__" not in dir(clipping):
            clipping = [clipping,clipping]
        # -- apply the masking 
        # ------------------
        # - Instrument
        radius = self.fetch_optimal_radius()
        pmap = self.instrument.get_stars_photomap(radius,catmag_range=catmag_range)

        delta_mag = pmap.refmap.mag- pmap.mag
        print delta_mag
        flaggood = [delta_mag==delta_mag]
        clipped = sigmaclip(delta_mag[flaggood] , *clipping)
        return clipped[0].mean(),clipped[0].std(),len(clipped[0])
    

    def estimate_fwhm(self,catmag_range=None,
                     isolated_only=True,in_units="pixel"):
        """
        Parameters
        ----------
        catmag_range: [2D-array / None]
                                   The magnitude range (of the catalogue) used to
                                   derive the zero point.
                                   If None the current one (self.catmag_range) will
                                   be used. otherwise this will overwirte it.

        isolated_only: [bool]      Use only the isolated stars (you should)

        units: [pixel or astropy units]

        Return
        ------
        3-float array (mean clipped (a,b [in 'units'] and theta [in rad])           
        """

        # ---------------- #
        # - output       - #
        # ---------------- #
        if type(in_units) == str:
            if in_units.lower() in ["pixels","pixel"]:
                return psf_a,psf_b,psf_t
            raise ValueError("the given 'units' must be 'pixels' or an astropy units")
        
        if "in_units" not in dir(in_units):
            raise TypeError("the given 'units' must be 'pixels' or an astropy units")
        
        print "TO BE DONE"

    def fetch_optimal_radius(self,catmag_range=None,isolated_only=True,**kwargs):
        """
        """
        print "This wont always work"
        return self.header["FWHM"]/2 / 0.32247221105013674 *3
        #apdict = cal.get_aperture_photometries(catmag_range=[14,18],isolated_only=isolated_only,
        #                              **kwargs)
        
        
    def get_aperture_photometries(self,pixel_range=[2,10],bins=10,
                                  catmag_range=None,isolated_only=True,
                                  **kwargs
                                  ):
        """This method loop over the given range of pixel radius
        and return the associated counts

        **kwargs goes to instrument.get_aperture ; so K Barbary's sep
        """
        self._test_instrument_()
            
        # -- apply the masking 
        if catmag_range is not None:
            self.catmag_range = catmag_range
        
        kwardsmask = dict(stars_only=True,isolated_only=isolated_only,
                          catmag_range=self.catmag_range)
        mask = self.sepobjects.get_indexes(**kwardsmask)
        # ---------------
        # - get the coordinate
        x,y = self.sepobjects.data["x"][mask],self.sepobjects.data["y"][mask]
        
        
        radius = np.linspace(pixel_range[0],pixel_range[1],bins)
        ap = {}
        for r_pixels in radius:
            counts,err,flag = self.instrument.get_aperture(x,y,r_pixels=r_pixels,**kwargs)
            ap[r_pixels] = {"counts":counts,
                            "errors":err,
                            "flag":flag}
        return ap

        
        
    # --------------- #
    # - IO Methods  - #
    # --------------- #
    def load_instrument(self,filename,astrotarget=None,
                        **kwargs):
        """
        """
        if type(filename) != str:
            raise TypeError("the 'filename' must be a string (location of the file)")

        inst = instrument(filename,astrotarget=astrotarget)
        self.create(inst,**kwargs)
        
        
    def create(self,instrument, force_it=False,
               load_catalogue=True,do_extraction=True,
               verbose=True):
        """
        """
        if self.has_instrument() and force_it is False:
            raise AttributeError("'instrument' is already defined."+\
                    " Set force_it to True if you really known what you are doing")
                    
        if "__nature__" not in dir(instrument) or instrument.__nature__ not in ["Image","BaseInstrument"]:
            raise TypeError("The given 'instrument' is not an astrobject's Image or Instrument")

        # -- Things looks good
        if load_catalogue and not instrument.has_catalogue():
            if verbose: print "Downloading the associated catalogue"
            instrument.download_catalogue()
        if do_extraction and not instrument.has_sepobjects():
            if verbose: print "Extracting the Image's sources (sep.extract)"
            instrument.sep_extract()

        self._properties["instrument"] = instrument
            
    
    # --------------- #
    # - Get Methods - #
    # --------------- #
    
    
    # --------------- #
    # - Set Methods - #
    # --------------- #
    
    # =========================== #
    # = Plotting Methods        = #
    # =========================== #
    def show(self):
        
        # -- Setting -- #
        from ..utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        if ax is None:
            fig = mpl.figure(figsize=[8,6])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        elif "hist" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure

        # -- Properties -- #
        print "TO BE DONE"
    
    
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _test_instrument_(self):
        if not self.has_instrument():
            raise AttributeError("no 'instrument' loaded")

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def instrument(self):
        return self._properties["instrument"]

    def has_instrument(self):
        """ Test if an instrument is loaded """
        return self.instrument is not None
        
    # -------------------
    # - Derived Properties
    @property
    def catalogue(self):
        self._test_instrument_()
        return self.instrument.catalogue
    
    @property
    def sepobjects(self):
        self._test_instrument_()
        return self.instrument.sepobjects

    @property
    def header(self):
        self._test_instrument_()
        return self.instrument.header
    # -------------------
    # - Mask
    @property
    def catmag_range(self):
        if self._side_properties["catmag_range"] is None:
            self._side_properties["catmag_range"] = [None,None]
        return self._side_properties["catmag_range"]

    @catmag_range.setter
    def catmag_range(self,value):
        if np.shape(value) !=  (2,):
            raise TypeError("'catmag_range' must be a 2d-array (minmag,maxmag)")
        self._side_properties["catmag_range"] = value
    
    @property
    def magmask(self):
        return self.sepobjects.get_refmag_mask(*self.catmag_range)
    
