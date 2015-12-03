#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sncosmo import get_bandpass
import pyfits as pf
import numpy as np
from ...utils.decorators import _autogen_docstring_inheritance
from ...utils.tools import kwargs_update

from ...astrobject.photometry import photopoint,Image
from ..baseobject import BaseObject
__all__ = ["Instrument"]



class Instrument( Image ):
    """
    """
    def __build__(self,data_index=0):
        """This is a slightly advanced Image object"""
        super(Instrument,self).__build__(data_index=data_index)

    # ----------- #
    #  PhotoPoint #
    # ----------- #
    @_autogen_docstring_inheritance(Image.get_aperture,"Image.get_aperture")
    def get_photopoint(self,x,y,r_pixels=None,
                       aptype="circle",
                       **kwargs):
        #
        # Be returns a PhotoPoint
        #
        count,err,flag  = self.get_aperture(x,y,r_pixels=r_pixels,aptype=aptype,
                                           **kwargs)
        flux = self.count_to_flux(count)
        var  = self.count_to_flux(err)**2
        return photopoint(self.lbda,flux,var,source="image",
                          instrument_name=self.instrument_name)
    
    @_autogen_docstring_inheritance(Image.get_target_aperture,"Image.get_target_aperture")
    def get_target_photopoint(self,r_pixels=None,
                              aptype="circle",**kwargs):
        #
        # Be returns a PhotoPoint
        #
        if not self.has_target():
            return AttributeError("No 'target' loaded")
        
        xpix,ypix = self.coords_to_pixel(self.target.ra,self.target.dec)
        return self.get_photopoint(xpix,ypix,r_pixels=r_pixels,
                                   aptype="circle",**kwargs)
    
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def count_to_flux(self,counts):
        return counts* 10**(-(2.406+self.mab0) / 2.5 ) / (self.lbda**2)
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    # ------------------
    # - Band Information
    @property
    def bandname(self):
        raise NotImplementedError("'band' must be implemented")
    
    @property
    def bandpass(self):
        return get_bandpass(self.bandname)
    
    @property
    def lbda(self):
        return self.bandpass.wave_eff

    # -- Derived values
    @property
    def mab0(self):
        raise NotImplementedError("'mab0' must be implemented")

    @property
    def _gain(self):
        raise NotImplementedError("'_gain' must be implemented (even for None)")
        
    @property
    def mjd_obstime(self):
        raise NotImplementedError("'obstime' must be implemented")
    

class Catalogue( BaseObject ):
    """
    """
    __nature__ = "Catalogue"
    source_name = "_not_defined_"
    
    _properties_keys = ["filename","data","header"]
    _side_properties_keys = ["wcs","fovmask"]
    _derived_properties_keys = ["fits","wcsx","wcsy"]

    
    def __init__(self, catalogue_file=None,empty=False,data_index=0):
        """
        """
        self.__build__(data_index=data_index)
        if empty:
            return
        self.load(catalogue_file)

    def __build__(self,data_index=0):
        """
        """
        super(Catalogue,self).__build__()
        self._build_properties = dict(
            data_index = data_index,
            header_ra = "X_WORLD",
            header_dec = "Y_WORLD",
            header_mag = "MAG",
            header_magerr = "MAGERR"
            )
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def set_wcs(self,wcs,force_it=False,update_fovmask=True):
        """
        """
        if self.has_wcs() and force_it is False:
            raise AttributeError("'wcs' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        if wcs is not None and "coordsAreInImage" not in dir(wcs):
            raise TypeError("'wcs' solution not recognize, should be an astLib.astWCS")

        self._side_properties["wcs"] = wcs
        
        if update_fovmask:
            if self.has_wcs():
                self.set_fovmask(wcs=self.wcs,update=False)
            else:
                print "None wcs"
                self._load_default_fovmask_()

            
    def set_fovmask(self, wcs=None,
                    ra_range=None, dec_range=None,
                    update=True):
        """
        This methods enable to define the mask of catalgue objects within the
        given field of view.
        
        Parameters
        ----------
        - options -
        update: [bool]             True to have a consistent object. Set False
                                   only if you know what you are doing
        Return
        ------
        Void
        """
        if wcs is not None:
            if "coordsAreInImage" not in dir(wcs):
                raise TypeError("'wcs' solution not recognize, should be an astLib.astWCS")
            
            self.fovmask = np.asarray([wcs.coordsAreInImage(ra,dec)
                                       for ra,dec in zip(self.data[self._build_properties["header_ra"]],
                                                         self.data[self._build_properties["header_dec"]])])
        elif ra_range is None or dec_range is None:
            raise AttributeError("please provide either 'wcs' and ra_range *and* dec_range")
        else:
            self.fovmask = (self.ra>ra_range[0]) & (self.ra<ra_range[0]) \
              (self.dec>dec_range[0]) & (self.dec<dec_range[0]) \

    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def load(self,catalogue_file,**kwargs):
        """
        """
        fits   = pf.open(catalogue_file)
        header = fits[self._build_properties["data_index"]].header
        data   = fits[self._build_properties["data_index"]].data

        self.create(data,header,**kwargs)
        self._properties["filename"] = catalogue_file
        self._derived_properties["fits"] = fits
        
    def create(self,data,header,force_it=True):
        """
        """
        if self.has_data() and force_it is False:
            raise AttributeError("'data' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        self._properties["data"] = data
        self._properties["header"] = header if header is not None \
          else pf.Header()

    # --------------------- #
    # PLOT METHODS          #
    # --------------------- #
    def display(self,ax,wcs_coords=True,draw=True,
                maskout=None,**kwargs):
        """
        This methods enable to show all the known sources in the
        image's field of view.
        This will only works if a catalogue has been set

        Parameters
        ----------

        ax: [matplotlib.axes]      the axes where the catalogue should be
                                   displaid

        maskout: [array]
        
        Return
        ------
        None (if no data) / ax.plot returns
        """
        if not self.has_data():
            print "Catalogue has no 'data' to display"
            return
        
        if wcs_coords:
            x,y = self.ra,self.dec
        else:
            x,y = self.wcs_xy

        if maskout is not None:
            x,y = x[~maskout],y[~maskout]
            
        default_prop = dict(
            ls="None",marker="o",mfc="b",mec="k",alpha=0.5,
            label="%s-catalgue"%self.source_name,
            scalex=False,scaley=False
            )
        prop = kwargs_update(default_prop,**kwargs)
        pl = ax.plot(x,y,**prop)
        
        if draw:
            ax.figure.canvas.draw()
        return pl
    
        
    # =========================== #
    # Properties and Settings     #
    # =========================== #    
    @property
    def data(self):
        return self._properties["data"]
    
    def has_data(self):
        return False if self.data is None\
          else True

    @property
    def header(self):
        return self._properties["header"]
    
    @property
    def wcs(self):
        return self._side_properties["wcs"]
        
    def has_wcs(self):
        return False if self.wcs is None\
          else True
        
    
    @property
    def fovmask(self):
        if self._side_properties["fovmask"] is None:
            self._load_default_fovmask_()
            
        return self._side_properties["fovmask"]
    
    def _load_default_fovmask_(self):
        self._side_properties["fovmask"] = np.ones(self.header["NAXIS2"],dtype=bool)
        
    @fovmask.setter
    def fovmask(self,newmask):
        if len(newmask) != self.header["NAXIS2"]:
            raise ValueError("the given 'mask' must have the size of 'ra'")
        self._side_properties["fovmask"] = newmask
        
    # ------------
    # - on flight
    # - coords
    @property
    def ra(self):
        """Barycenter position along world x axis"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["header_ra"]][self.fovmask]

    @property
    def dec(self):
        """arycenter position along world y axis"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["header_dec"]][self.fovmask]
    
    # - mag
    @property
    def mag(self):
        """Generic magnitude"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["header_mag"]][self.fovmask]

    @property
    def mag_err(self):
        """Generic magnitude RMS error"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["header_magerr"]][self.fovmask]

    
    # ------------
    # - Derived
    @property
    def fits(self):
        return self._derived_properties["fits"]

    @property
    def wcs_xy(self):
        if self.has_wcs():
            return np.asarray([self.wcs.wcs2pix(ra_,dec_) for ra_,dec_ in zip(self.ra,self.dec)]).T
        raise AttributeError("no 'wcs' solution loaded")
    
