#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module use the instrument objects to perform photometry"""
import numpy as np
from scipy.stats import sigmaclip
from astropy import units
from ..astrobject.baseobject import BaseObject
from ..utils.tools import kwargs_update
from ..utils.decorators import _autogen_docstring_inheritance
from ..astrobject.instruments.instrument import instrument 

class PhotoCatalogue( BaseObject ):
    """
    """
    __nature__ = "PhotoCatalogue"

    _properties_keys = ["instrument"]
    _side_properties_keys = ["catmag_range"]
    _derived_properties_keys = ["fwhm","zp_calibration",
                                "photo_apertures"]

    
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
    # = Main Methods            = #
    # =========================== #
    def calibrate(self, catmag_range=[None,None],
                  isolated_only=True,clipping=3,
                  redo_aperture_photos=False):
        """
        """
        maskdict = dict(catmag_range=catmag_range,
                        isolated_only=isolated_only)
        if not self.instrument.has_apertures_photos() or redo_aperture_photos or \
          (catmat_range is not None and catmag_range != self.aperture_photos["properties"]["catmag_range"]):
            self.instrument.set_apertures_photos(**maskdict)
        if not self.instrument.has_fwhm():
            self.instrument.derive_fwhm(set_it=True,force_it=False)
            
        self.derive_zeropoint(clipping=clipping,**maskdict)

        
    # --------------- #
    # - ZeroPoints  - #
    # --------------- #
    def derive_zeropoint(self,catmag_range=None,
                         isolated_only=True,
                         fwhm_scale=3,clipping=3,
                         set_it=True,force_it=False):
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
        # ------------------
        # - Input
        if "__iter__" not in dir(clipping):
            clipping = [clipping,clipping]

        # ------------------
        # - Instrument
        radius_pixel=self.fwhm.value /2/self.instrument.pixel_size_arcsec * fwhm_scale
        pmap = self.instrument.get_stars_photomap(radius_pixel,
                                                  catmag_range=catmag_range,
                                                  isolated_only=isolated_only)

        # -----------------
        # - Cleaning 
        delta_mag = pmap.refmap.mag- pmap.mag
        flaggood = [delta_mag==delta_mag]
        clipped = sigmaclip(delta_mag[flaggood] , *clipping)

        # -----------------
        # - Output
        calibration = {"zp":self.instrument.mab0 + clipped[0].mean(),
                       "zp_std":clipped[0].std(),
                       "npoint_used":len(clipped[0]),
                       "properties": {"catmag_range":catmag_range,
                                        "isolated_only":isolated_only,
                                        "radius_pixel":radius_pixel,
                                        "clipping":clipping,
                                        "clipping_bounds":np.asarray(clipped[1:]) + self.instrument.mab0 ,
                                        }
                        }
        if set_it:
            self.set_zpcalibration(calibration)
        else:
            return calibration
        

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
               verbose=True,catmag_range=[None,None],**kwargs):
        """
        """
        if self.has_instrument() and not force_it :
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
    def set_zpcalibration(self,zpdict,force_it=False):
        """
        """
        if self.has_zpcalibration() and not force_it:
            raise AttributeError("the 'zeropoint calibration' is already set."+\
                    " Set force_it to True if you really known what you are doing")
                    
        cal = {"zp":None,"zp_err":None,"derivation_prop":{}}
        self._derived_properties['zp_calibration'] = kwargs_update(cal,**zpdict)
        
        if self.zpcalibration["zp"] is not None:
            try:
                self.instrument.mab0 = self.zpcalibration["zp"]
            except:
                print "WARNING no '@setter.mab0' defined for the given instrument. Can't update it."
                
            
    
        
        
                
    # =========================== #
    # = Plotting Methods        = #
    # =========================== #
    def show_aperture_photos(self,savefile=None,ax=None,show=True,
                                   cmap=None,cbar=True,set_axes_labels=True,
                                   show_fwhm=True,
                                   **kwargs):
        return self.instrument.show_aperture_photos(savefile=savefile,ax=ax,show=show,
                                    cmap=cmap,cbar=cbar,set_axes_labels=set_axes_labels,
                                   show_fwhm=show_fwhm,
                                   **kwargs)
    
        
    def show_zpcalibration(self,savefile=None,ax=None,show=True,
                           colorin=None,colorout=None,
                           catmag_range=None,
                           percentshade=[1,2,5,10],
                           **kwargs):
        """
        """
        # ---------------
        # -- Setting 
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

        show_catmag_range = self.zpcalibration["properties"]["catmag_range"] \
          if catmag_range is None else catmag_range
          
                                                  
        pmap = self.instrument.get_stars_photomap(self.zpcalibration["properties"]["radius_pixel"],
                                                  catmag_range=show_catmag_range,
                                                  isolated_only=self.zpcalibration["properties"]["isolated_only"],
                                                  )
        detlamag =  pmap.refmap.mag - pmap.mag
        clipping = self.zpcalibration["properties"]["clipping_bounds"] - self.zpcalibration["zp"]
        maskout = (detlamag<clipping[0]) + (detlamag>clipping[1])
        maskmagout = (pmap.refmap.mag<self.zpcalibration["properties"]["catmag_range"][0]+0.01) +\
           (pmap.refmap.mag>self.zpcalibration["properties"]["catmag_range"][1]-0.01)
        
        # -------------------------
        # - Properties
        if colorin is None:
            colorin = mpl.cm.Blues(0.7)
        if colorout is None:
            colorout = mpl.cm.Reds(0.7)

        default_prop = {"marker":"o", "ms":10,"ls":"None","zorder":5}
        prop = kwargs_update(default_prop,**kwargs)
        # -------------------------
        # - Da plot itself
        # -- The Good
        pl = ax.plot(pmap.refmap.mag[~maskout & ~maskmagout],detlamag[~maskout & ~maskmagout],
                     color=colorin,**prop)
        # -- The Bad
        pl.append(ax.plot(pmap.refmap.mag[maskout& ~maskmagout],detlamag[maskout& ~maskmagout],
                          color=colorout,**prop))
        # -- And the Ugly
        if maskmagout.any():
            propout_mag = kwargs_update(default_prop,**{"mec":colorin,"mfc":"None","mew":2})
            pl.append(ax.plot(pmap.refmap.mag[maskmagout],detlamag[maskmagout],
                              **propout_mag))
        
        pl.append(ax.errorbar(pmap.refmap.mag,detlamag,
                              yerr=pmap.magerr,ecolor="0.7",zorder=3,
                              label="_no_legend_",ms=0,ls="None"))
        fancy = []
        fancy.append(ax.axhline(0,ls="-",color="0.7",zorder=2))
                
        # -------------------
        # - Percent Curves
        
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["pl"]     = pl

        fig.figout(savefile=savefile,show=show)
        
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
    # - Derived Values
    @property
    def apertures_curve(self):
        return self.instrument.apertures_curve[1]

    @property
    def apertures_radii(self):
        return self.instrument.apertures_curve[0]

    @property
    def fwhm(self):
        return self.instrument.fwhm
    
    # ZERO POINTS
    @property
    def zpcalibration(self):
        """
        """
        if self._derived_properties['zp_calibration'] is None:
            self.set_zpcalibration({}) 
        return self._derived_properties['zp_calibration']
    
    def has_zpcalibration(self):
        return not (self._derived_properties['zp_calibration'] is None \
                    or self.zpcalibration['zp'] is None)

    # -------------------
    # - Mask
        
