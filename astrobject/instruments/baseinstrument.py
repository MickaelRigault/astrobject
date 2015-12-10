#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sncosmo import get_bandpass
from astropy import coordinates
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
    _side_properties_keys = ["wcs","fovmask","matchedmask"]
    _derived_properties_keys = ["fits","wcsx","wcsy"]

    
    def __init__(self, catalogue_file=None,empty=False,
                 data_index=0,key_mag=None,key_magerr=None):
        """
        """
        self.__build__(data_index=data_index,
                       key_mag=key_mag,key_magerr=key_magerr)
        if empty:
            return
        if catalogue_file is not None:
            self.load(catalogue_file)

    def __build__(self,data_index=0,
                  key_mag=None,key_magerr=None):
        """
        """
        super(Catalogue,self).__build__()
        self._build_properties = dict(
            data_index = data_index,
            key_ra = "X_WORLD",
            key_dec = "Y_WORLD"
            )
        self.set_mag_keys(key_mag,key_magerr)
        
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
                                       for ra,dec in zip(self.data[self._build_properties["key_ra"]],
                                                         self.data[self._build_properties["key_dec"]])])
        elif ra_range is None or dec_range is None:
            raise AttributeError("please provide either 'wcs' and ra_range *and* dec_range")
        else:
            self.fovmask = (self.ra>ra_range[0]) & (self.ra<ra_range[0]) \
              (self.dec>dec_range[0]) & (self.dec<dec_range[0]) \

    def set_matchedmask(self,matchedmask):
        """This methods enable to set to matchedmask, this mask is an addon
        mask that indicate which point from the catalogue (after the fov cut)
        has been matched by for instance a sextractor/sep extraction
        """
        if matchedmask is None or len(matchedmask) == 0:
            return

        if type(matchedmask[0]) is bool:
            # - it already is a mask, good
            self._side_properties["matchedmask"] = np.asarray(matchedmask,dtype=bool)
            return

        if type(matchedmask[0]) in [int,np.int32,np.int64]:
            # - it must be a list of matched index (from SkyCoord matching fuction e.g.)
            self._side_properties["matchedmask"] = np.asarray([i in matchedmask for i in range(len(self.ra))],
                                                              dtype=bool)
            return
        raise TypeError("the format of the given 'matchedmask' is not recongnized. "+\
                        "You could give a booleen mask array or a list of accepted index")
                        
        
    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def load(self,catalogue_file,**kwargs):
        """
        """
        fits   = pf.open(catalogue_file)
        header = fits[self._build_properties["data_index"]].header
        data   = fits[self._build_properties["data_index"]].data
        if type(data) == pf.fitsrec.FITS_rec:
            from astrobject.utils.tools import fitsrec_to_dict
            from astropy.table import TableColumns
            data = TableColumns(fitsrec_to_dict(data))
            
        self.create(data,header,**kwargs)
        self._properties["filename"] = catalogue_file
        self._derived_properties["fits"] = fits

        
    def create(self,data,header,force_it=True,**build):
        """
        """
        if self.has_data() and force_it is False:
            raise AttributeError("'data' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        self._properties["data"] = data
        self._properties["header"] = header if header is not None \
          else pf.Header()
        self._build_properties = kwargs_update(self._build_properties,**build)
      
    def writeto(self,savefile,force_it=True):
        """
        """
        raise NotImplementedError("to be done")
    # --------------------- #
    # Set Methods           #
    # --------------------- #
    def set_mag_keys(self,key_mag,key_magerr):
        """
        """
        self._build_properties["key_mag"] = key_mag
        self._build_properties["key_magerr"] = key_magerr
        
    # --------------------- #
    # PLOT METHODS          #
    # --------------------- #
    def display(self,ax,wcs_coords=True,draw=True,
                apply_machedmask=True,
                show_nonmatched=True,propout={},
                **kwargs):
        """
        This methods enable to show all the known sources in the
        image's field of view.
        This will only works if a catalogue has been set

        Parameters
        ----------

        ax: [matplotlib.axes]      the axes where the catalogue should be
                                   displaid

        
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
        # -------------- #
        # - mask
        mask = self.matchedmask if self.has_matchedmask() \
          else np.ones(len(self.ra),dtype=bool)

        starmask = self.starmask if self.has_starmask() \
          else np.ones(len(self.ra),dtype=bool)
          
        # -- in / out star / notstar
        x_starin,y_starin = x[mask & starmask], y[mask & starmask]
        x_nostarin,y_nostarin = x[mask & ~starmask], y[mask & ~starmask]
        
        x_starout,y_starout = x[~mask & starmask],y[~mask & starmask]
        x_nostarout,y_nostarout = x[~mask & ~starmask], y[~mask & ~starmask]
        # -- Properties
        colorin = "b"
        colorout = "r"
        default_prop = dict(
            ls="None",marker="o",mfc="b",alpha=0.7,ms=6,mew=0,
            label="%s-catalgue"%self.source_name,
            )
        prop = kwargs_update(default_prop,**kwargs)

        # -- plot loop
        axout = []
        for x_,y_,show_,propextra in [[x_starin,  y_starin,True,
                                       {}],
                                      [x_nostarin,y_nostarin,True,
                                        dict(mfc="None",mec=colorin,mew=2,alpha=0.6)],
                                      [x_starout, y_starout, show_nonmatched,
                                        dict(mfc=colorout,ms=4,alpha=0.5)],
                                      [x_nostarout,y_nostarout,show_nonmatched,
                                        dict(mfc="None",mec=colorout,mew=1,ms=4,alpha=0.5)],
                                ]:
            prop_ = kwargs_update(prop,**propextra)
            if len(x_)>0 and show_:
                axout.append(ax.plot(x_,y_,**prop_))
        
        if draw:
            ax.figure.canvas.draw()
            
        return axout
    
        
    # =========================== #
    # Properties and Settings     #
    # =========================== #
    # -------
    # - data
    @property
    def data(self):
        return self._properties["data"]
    
    def has_data(self):
        return False if self.data is None\
          else True
          
    @property
    def ndata(self):
        if self.data is None:
            return 0
        
        if self.header is not None and "NAXIS2" in self.header:
            return self.header["NAXIS2"]
        
        return len(self.data.values()[0])
    
    # ---------------
    # - header / wcs 
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
        self._side_properties["fovmask"] = np.ones(self.ndata,dtype=bool)
        
    @fovmask.setter
    def fovmask(self,newmask):
        if len(newmask) != self.ndata:
            raise ValueError("the given 'mask' must have the size of 'ra'")
        self._side_properties["fovmask"] = newmask


    @property
    def matchedmask(self):
        return self._side_properties["matchedmask"]

    def has_matchedmask(self):
        return False if self.matchedmask is None\
          else True
          
    # ------------
    # - on flight
    # - coords
    @property
    def ra(self):
        """Barycenter position along world x axis"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["key_ra"]][self.fovmask]

    @property
    def dec(self):
        """arycenter position along world y axis"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["key_dec"]][self.fovmask]

    @property
    def sky_radec(self):
        """This is an advanced radec methods tight to astropy SkyCoords"""
        return coordinates.SkyCoord(ra=self.ra,dec=self.dec, unit="deg")
    
    # - mag
    @property
    def mag(self):
        """Generic magnitude"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        if not self._is_keymag_set_():
            raise AttributeError("no 'key_mag' defined. see self.set_mag_keys ")
        
        return self.data[self._build_properties["key_mag"]][self.fovmask]

    @property
    def mag_err(self):
        """Generic magnitude RMS error"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        if not self._is_keymag_set_():
            raise AttributeError("no 'key_magerr' defined. see self.set_mag_keys ")
        
        return self.data[self._build_properties["key_magerr"]][self.fovmask]

    def _is_keymag_set_(self,verbose=True):
        """this method test if the keymag has been set"""
        if self._build_properties["key_mag"] is None or \
          self._build_properties["key_magerr"] is None:
            if verbose:
                print "No 'key_mag'/'key_magerr' set ; call 'set_mag_keys'. List of potential keys: "\
                  +", ".join([k for k in self.data.keys() if "mag" in k or "MAG" in k])
            return False
        return True
    
    @property
    def objectclass(self):
        if "key_class" not in self._build_properties.keys():
            raise AttributeError("no 'key_class' provided in the _build_properties.")
        
        return self.data[self._build_properties["key_class"]][self.fovmask]
        
    @property
    def starmask(self):
        """ This will tell which of the datapoints is a star
        Remark, you need to have defined key_class and value_star
        in the __build_properties to be able to have access to this mask
        """
        if "value_star" not in self._build_properties.keys():
            return None

        flag = self.objectclass == self._build_properties["value_star"]
        return flag.data  #not self.fovmask already in objectclass

    def has_starmask(self):
        return False if self.starmask is None \
          else True
          
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
    
