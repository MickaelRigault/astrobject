#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module use the instrument objects to perform photometry"""
import numpy as np
from scipy.stats import sigmaclip
from astropy import units

from ..astrobject.baseobject import BaseObject
from ..utils.tools import kwargs_update
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

        # -----------------
        # - TO BE REMOVED
        self.derive_fwhm()
        
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def calibrate(self, catmag_range=[None,None],
                  isolated_only=True,clipping=3):
        """
        """
        maskdict = dict(catmag_range=catmag_range,
                        isolated_only=isolated_only)
        
        self.derive_fwhm(set_it=True,force_it=False,**maskdict)
        delta_zp,delta_zp_std,npoint_used = self.derive_zeropoint(clipping=clipping,**maskdict)
        
        self.instrument.mab0 += delta_zp
        
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
        radius_pixel=self.fwhm/2*fwhm_scale/self.instrument.pixel_size_arcsec
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
        
    def derive_fwhm(self,catmag_range=None,
                      isolated_only=True,
                      set_it=True,force_it=False):
        """
        Parameters
        ----------
        catmag_range: [2D-array / None]
                                   The magnitude range (of the catalogue) used to
                                   derive the zero point.
                                   If None the current one (self.catmag_range) will
                                   be used. otherwise this will overwirte it.

        isolated_only: [bool]      Use only the isolated stars (you should)

        Return
        ------
        float (if set_it is False)
        """
        print "*WARNING* STELLA BASE FOR KNOW"
        fwhm = self.header["FWHM"]
        
        if set_it:
            self.set_fwhm(fwhm,force_it=force_it)
        else:
            return fwhm
        
        
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
    def get_aperture_photometries(self,rpixel_range=[0.5,20],bins=18,
                                  catmag_range=None,isolated_only=True,
                                  **kwargs
                                  ):
        """This method loop over the given range of pixel radius
        and return the associated counts

        **kwargs goes to instrument.get_aperture ; so K Barbary's sep
        """
        # ------------------
        # - Input
        self._test_instrument_()
        if catmag_range is not None:
            self.catmag_range = catmag_range

        # ----------------
        # - Aperture Info
        radius = np.linspace(rpixel_range[0],rpixel_range[1],bins)
        
        kwardsmask = dict(stars_only=True,isolated_only=isolated_only,
                          catmag_range=self.catmag_range)
        mask = self.sepobjects.get_indexes(**kwardsmask)
        catidx = self.sepobjects.get_indexes(cat_indexes=True,**kwardsmask)
        # ---------------------
        # - get the coordinate
        x,y = self.sepobjects.data["x"][mask],self.sepobjects.data["y"][mask]

        # ---------------
        # - SHOUD BE FASTER
        ap = {}
        for r_pixels in radius:
            counts,err,flag = self.instrument.get_aperture(x,y,r_pixels=r_pixels,**kwargs)
            ap[r_pixels] = {"counts":counts,
                            "errors":err,
                            "flag":flag}
        return {"radius":radius,
                "counts":np.asarray([ap[d]["counts"] for d in radius]),
                "errors":np.asarray([ap[d]["errors"] for d in radius]),
                "idx_catalogue":catidx,
                "idx":mask,
                "property":kwargs_update( dict(catmag_range=catmag_range,
                                        isolated_only=isolated_only),
                                        **kwargs)
                }
    
    # --------------- #
    # - Set Methods - #
    # --------------- #
    def set_fwhm(self,value,force_it=True):
        """
        """
        if self.has_fwhm() and not force_it:
            raise AttributeError("'fwhm' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        if value<0:
            raise ValueError("the 'fwhm' must be positive")
        self._derived_properties['fwhm'] = value

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
                raise AttributeError("Not '@setter.mab0' defined for the given instrument.")
            
    def set_photo_apertures(self,rpixel_range=[0.5,20],bins=18,
                            catmag_range=None,isolated_only=True,
                            **kwargs):
        """
        """

        self._derived_properties["photo_apertures"] = \
          self.get_aperture_photometries(rpixel_range=rpixel_range,bins=bins,
                                         catmag_range=catmag_range,
                                        isolated_only=isolated_only,
                                         **kwargs)
        
         
        
        
                
    # =========================== #
    # = Plotting Methods        = #
    # =========================== #

        
    def show_aperture_photometries(self,savefile=None,ax=None,show=True,
                                   cmap=None,cbar=True,**kwargs
                                   ):
        """
        """
        if not self.has_photo_apertures():
            raise AttributeError("no 'photo_apertures' defined")
        
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
        
        # -- Properties -- #

        # -----------------
        # - Data To show
        counts = self.photo_apertures['counts']
        # ------------------
        # - Colors  
        catmag = self.catalogue.mag[self.photo_apertures['idx_catalogue']]
        catmag_range = self.photo_apertures['property']['catmag_range']
        if catmag_range is None:
            catmag_range = [catmag.min(),catmag.max()]
        elif catmag_range[0] is None:
            catmag_range[0] = catmag.min()
        elif catmag_range[1] is None:
            catmag_range[1] = catmag.max()
            
        cmap = mpl.cm.jet if cmap is None else cmap
        colors = cmap((catmag-catmag_range[0])/(catmag_range[1]-catmag_range[0]))
        
        pl = [ax.plot(self.photo_apertures['radius'],c_ / np.mean(c_), color=color,**kwargs)
              for c_,color in zip(counts.T,colors)]

        # ----------------------
        # - colorbar
        # - this means it is not an ax
        if cbar:
            from ..utils.mpladdon import colorbar,insert_ax
            if "imshow" not in dir(cbar):
                axcar = ax.insert_ax(space=.05,pad=0.03,location="right")
            else:
                axcar = cbar
            calpha = kwargs.pop("alpha",None)
            cbar = axcar.colorbar(cmap,vmin=catmag_range[0],vmax=catmag_range[1],
                                label="catalogue magnitude",alpha=calpha)
            
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["pl"]     = pl
    
        fig.figout(savefile=savefile,show=show)
        
        return self._plot
    
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
    # FWHM
    @property
    def fwhm(self):
        if not self.has_fwhm():
            raise AttributeError("'fwhm' is not defined")
        return self._derived_properties["fwhm"]

    def has_fwhm(self):
        return not self._derived_properties["fwhm"] is None

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

    # APERTURE PHOTOMETRIES
    @property
    def photo_apertures(self):
        """
        """
        if self._derived_properties['photo_apertures'] is None:
            self._derived_properties['photo_apertures'] = {}
            
        return self._derived_properties['photo_apertures']


    def has_photo_apertures(self):
        return not len(self.photo_apertures.keys())==0

    
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
        
