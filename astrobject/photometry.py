#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module defines the photometric objects"""


import numpy  as np
import pyfits as pf

from astropy     import units
from astLib      import astWCS

from .baseobject   import BaseObject 
from ..utils.tools import kwargs_update


__all__ = ["image"]


def image(filename=None,astrotarget=None,**kwargs):
    """
    Initalize the image by giving its filelocation (*filename*). This
    will load it using the load() method.

    Parameters
    ----------

    - options -

    filename: [string.fits]    fits file from where the image will be loaded

    astrotarget: [AstroTarget] An AstroTarget object you which to associate
                               to this image.
                                                                  
    - kwargs options, potentially non-exhaustive -
    
    empty: [bool]              Set True to load an empty object.
                               
    Return
    ------
    Image
    """
    return Image(filename,astrotarget=astrotarget,
                 **kwargs)#.copy()


#######################################
#                                     #
# Base Object Classes                 #
#                                     #
#######################################
class Image( BaseObject ):
    """
    """
    # -------------------- #
    # Internal Properties  #
    # -------------------- #
    _properties_keys         = ["filename","rawdata","header",
                                "var","background"]
    _side_properties_keys    = ["wcs","target",
                                "exptime"]
    _derived_properties_keys = ["fits","data"]
    
    # Where in the fitsfile the data are
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,filename=None,
                 astrotarget=None,
                 empty=False,**kwargs):
        """
        Initalize the image by giving its filelocation (*filename*). This
        will load it using the load() method.

        Parameters
        ----------
        filename: [string.fits]    fits file from where the image will be loaded
                                   - Trick - Set None and no image will be loaded

        astrotarget: [AstroTarget] An AstroTarget object you which to associate
                                   to this image. 
                                   
        empty: [bool]              Does not do anything, just loads an empty object.
                                   (Careful with that)
                                   
        """
        self.__build__()
        if empty:
            return
        
        if filename is not None:
            force_it = kwargs.pop("force_it",True)
            self.load(filename,force_it=force_it,
                      **kwargs)
            
    def __build__(self):
        #
        # Improvement of BaseObject
        # including the _object_properties
        #
        super(Image,self).__build__()
        # -- How to read the image
        self._build_properties = dict(
                data_index = 0,
                header_exptime = "EXPTIME"
                )

    # =========================== #
    # = Main Methods            = #
    # =========================== #
    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def load(self,filename,index=None,
             force_it=False):
        """
        This enables to load a fitsfile image and will create
        the basic data and wcs solution if possible.
        *var* (error) and *background* has to be defined
        separately has this strongly depend on the instrument

        Parameters
        ----------
        filename: [string.fits]    The file containing the fits data

        - options -
        
        index: [int]               The fits entry that contains the data
                                   If None, this will fetch it in the build_properties

        force_it: [bool]           If the data already exist, this method
                                   will raise an execption except if you set
                                   *force_it* to True. Be Careful with this.
        Return
        ------
        Void
        """
        # -- Check if you will not overwrite anything
        if self.data is not None and force_it is False:
            raise AttributeError("'data' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        index = self._build_properties["data_index"] if index is None \
          else index
          
        # -------------------------- #
        #  Check The input           #
        # -------------------------- #
        try:
            fits = pf.open(filename)
        except:
            raise TypeError("'filename' cannot be open by pyfits.")

        data = fits[index].data
        if data is None:
            raise ValueError("no 'data' in the given fits file for the given index")
        
        # -------------------------- #
        #  Everythin looks good !    #
        # -------------------------- #
        try:
            wcs_ = astWCS.WCS(filename,
                              extensionName=index)
        except:
            wcs_ = None
            
        self.create(data,None,wcs_,
                    background= None,
                    header=fits[index].header,
                    filename=filename,fits=fits,
                    force_it=True)
        
    def create(self,rawdata,variance,wcs,
               background=None,header=None,exptime=None,
               filename=None,fits=None,force_it=False):
        """
        Create the image-object using by filling its different component.
        Each of them can be None, but in that case, their corresponding method
        and derived values won't be avialable.

        Parameters
        ----------
        force_it: [bool]           If the data already exist, this method
                                   will raise an execption except if you set
                                   *force_it* to True. Be Careful with this.

        Return
        ------
        Void
        """
        # -- Check if you will not overwrite anything
        if self.data is not None and force_it is False:
            raise AttributeError("'data' is already defined."+\
                    " Set force_it to True if you really known what you are doing")
        # ********************* #
        # * Create the object * #
        # ********************* #
        self._properties["filename"]     = filename
        self._derived_properties["fits"] = fits
        # basics data and wcs
        self._properties["rawdata"]      = np.asarray(rawdata,dtype="float")
        self._properties["header"]       = pf.Header if header is None else header
        self._properties["var"]          = variance
        self.set_background(background)
        # -- Side, exposure time
        self._side_properties["exptime"] = \
          self.header[self._build_properties["header_exptime"]] if exptime is None \
          else exptime
          
        if wcs is not None and wcs.__module__ != "astLib.astWCS":
            print "WARNING: only astLib.astWCS wcs solution is implemented"
            print " ----> No wcs solution loaded"
            self._side_properties["wcs"]     = None
        else:
            self._side_properties["wcs"]     = wcs

        self._update_()

    # ------------------- #
    # - Target          - #
    # ------------------- #
    def set_target(self,newtarget,test_inclusion=True):
        """
        Change (or create) an object associated to the given image.
        This function will test if the object is withing the image
        boundaries (expect if *test_inclusion* is set to False).
        Set newtarget to None to remove the association between this
        object and a target
        """
        if newtarget is None:
            self._side_properties['target'] = None
            return
        
        # -- Input Test -- #
        if newtarget.__nature__ != "AstroTarget":
            raise TypeError("'newtarget' should be (or inherite) an AstroTarget")
        
        if test_inclusion:
            if self.has_wcs() is False:
                print "WARNING: because there is no wcs solution, "+\
                  "I can't test the inclusion of the new astrotarget"
            else:
                if not self.wcs.coordsAreInImage(*newtarget.radec):
                    raise ValueError("The new 'target' is not inside the image "+\
                                      " boundaries"+ "\n"+\
                                     "--> object radec: %.3f,%.4f"%(newtarget.ra,
                                                                    newtarget.dec))
        # -- Seems Ok -- #
        self._side_properties["target"] = newtarget.copy()

    def set_background(self,background,
                       force_it=False,check=True):
        """
        # Might not make sense in image #
        This is a method that might strongly depend on the instrument.
        As a default (background = None) this uses Sextractor background
        estimation from 'get_sep_background'. Give background or overwrite
        this method for your specific instrument.

        Parameters
        ----------
        background: [2d-array]     This is the float-array containing the background
                                   level of the data.
                                   *Important* if None is set, this method will call
                                   Sextractor's sep-Background module to estimate it.

        - option -
        force_it: [bool]           If the data already exist, this method
                                   will raise an execption except if you set
                                   *force_it* to True. Be Careful with this.
                                   
        check: [bool]              If True, this will check tha the given background 
                                   is consistant with the existing data (rawdata)
                                   if any. Set that to False to avoid this check.

        Return
        ------
        Void
        """
        # -- Check if you will not overwrite anything
        if self.background is not None and force_it is False:
            raise AttributeError("'background' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        background  = self.get_sep_background() if background is None else background
        
        # Shape test
        if self.rawdata is not None and np.shape(background) != self.shape:
            raise ValueError("The given background must have rawdata's shape")
        # -- Looks good
        self._properties['background'] = np.asarray(background)
        self._update_data_()
    # ------------------- #
    # - get Methods     - #
    # ------------------- #
    def get_pixels_around(self,pixel_x,pixel_y,radius_in_pixel):
        """
        """
        # - This should be moved in a c-lib
        raise NotImplementedError("this will be kyle's SEP. Soon.")

    def pixel_to_coords(self,pixel_x,pixel_y):
        """get the coordinate (ra,dec; degree) associated to the given pixel (x,y)"""
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded.")
        
        return self.wcs.pix2wcs(pixel_x,pixel_y)

    def coords_to_pixel(self,ra,dec):
        """Return the pixel (x,y) associated to the given ra,dec (degree) coordinate"""
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded.")
        
        return self.wcs.wcs2pix(ra,dec)

    # ------------------- #
    # - SEP Tools       - #
    # ------------------- #
    def get_sep_background(self):
        """
        This module is based on K. Barbary's python module of Sextractor: sep.
        
        """
        if "_background" not in dir(self):
            if self.rawdata is None:
                raise ValueError("no 'rawdata' loaded. Cannot get a background")
            from sep import Background
            self._background = Background(self.rawdata)

        return self._background.back()
        
    def sep_extract(self):
        """
        This module is based on K. Barbary's python module of Sextractor SEP.
        """
        print "to be done"
        
    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def writeto(self):
        """
        """
        print "to be done"

    def show(self,toshow="data",savefile=None,logscale=True,
             ax=None,show=True,wcs_coords=False,proptarget={},
             **kwargs):
        """
        """
        # ----------------
        # - Input test
        if toshow not in dir(self):
            raise ValueError("'%s' is not a known image parameter"%toshow)
        valuetoshow = eval("self.%s"%toshow)
        if valuetoshow is None:
            raise AttributeError("no '%s' to show (=None)"%toshow)
        
        # -- Setting -- #
        from ..utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
        # ----------- #
        # -  What
        x = np.log10(valuetoshow) if logscale else valuetoshow
        
        # ----------- #
        # - How
        default_prop = {
            "interpolation":"nearest",
            "origin":"lower"
            }
        if self.has_wcs() and wcs_coords:
            default_prop["extent"] = self.worldcoords_boundaries
            ax.set_xlabel(r"$\mathrm{Ra\ [deg]}$",fontsize = "x-large")
            ax.set_ylabel(r"$\mathrm{Dec\ [deg]}$",fontsize = "x-large")
        else:
            wcs_coords = False
            #default_prop["extent"] = [0,self.width,0,self.height]
            ax.set_xlabel("x",fontsize = "x-large")
            ax.set_ylabel("y",fontsize = "x-large")
            
        prop = kwargs_update(default_prop,**kwargs)

        # ----------- #
        # - Do It
        im = ax.imshow(x,**prop)
        # - add target
        pl_tgt = None if self.has_target() is False \
          else self.display_target(ax,wcs_coords=wcs_coords,
                                   **proptarget)
        # ----------- #
        # - Recordit
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["imshow"] = im
        self._plot["target_plot"] = pl_tgt
        self._plot["prop"]   = prop
        self._plot["wcs_coords"] = wcs_coords
        
        fig.figout(savefile=savefile,show=show)
        
        return self._plot

    def display_target(self,ax,wcs_coords=True,**kwargs):
        """If a target is loaded, use this to display the target on the
        given ax"""
        if self.has_target() is False:
            print "No target to display"
            return
        # --------------------
        # - Fancy
        default_markerprop = {
            "marker":"s",
            "mfc":"w",
            "mec":"k","mew":2,
            "zorder":12,
            "scalex":False,"scaley":False
            }
        prop = kwargs_update(default_markerprop,**kwargs)
    
        if wcs_coords:
            pl = ax.plot(self.target.ra,self.target.dec,**prop)
        else:
            radec_pixel = self.coords_to_pixel(*self.target.radec)
            pl = ax.plot(radec_pixel[0],radec_pixel[1],**prop)
                
        return pl

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def filename(self):
        return self._properties["filename"]

    @property
    def fits(self):
        return self._derived_properties["fits"]
    
    @property
    def data(self):
        return self._derived_properties["data"]
    # -------
    # -- This under defines the data
    # Raw Data
    @property
    def rawdata(self):
        return self._properties["rawdata"]

    @rawdata.setter
    def rawdata(self):
        if self.rawdata is not None and np.shape(value) != self.shape:
            raise ValueError("'rawdata' cannot change shape using this setter. ")
        
        self._properties['rawdata'] = np.asarray(value)
        self._update_data_()
        
    # Background
    @property
    def background(self):
        return self._properties["background"]
    
    @background.setter
    def background(self,value):
        if self.background is not None and np.shape(value) != self.shape:
            raise ValueError("'background' cannot change shape using this setter. ")
        
        self._properties['background'] = np.asarray(value)
        self._update_data_()
        
    # -- Header stuff
    @property
    def header(self):
        return self._properties["header"]
    
    @property
    def exposuretime(self):
        return float(self._side_properties['exptime']) \
          if self._side_properties['exptime'] is not None else None
    
    # -- derived values
    @property # based on rawdata for pratical reason (like background check)
    def shape(self):
        if self.rawdata is None:
            raise AttributeError("No rawdata loaded.")
        return np.shape(self.rawdata)
    
    @property
    def width(self):
        return self.shape[1]
    @property
    def height(self):
        return self.shape[0]

    
    @property
    def var(self):
        return self._properties["var"]
    
    def has_var(self):
        return False if self.var is None \
          else True

    # ------------      
    # -- wcs tools
    @property
    def wcs(self):
        return self._side_properties["wcs"]
    
    def has_wcs(self):
        return False if self.wcs is None \
          else True
          
    @property
    def pixel_size_deg(self):
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded")
        return self.wcs.getPixelSizeDeg()
    
    @property
    def pixel_size_arcsec(self):
        return self.pixel_size_deg*units.degree.in_units("arcsec")

    @property
    def worldcoords_boundaries(self):
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded")
        return self.wcs.getImageMinMaxWCSCoords()
    
    # ----------------      
    # -- target tools
    @property
    def target(self):
        return self._side_properties["target"]
    
    def has_target(self):
        return False if self.target is None \
          else True
    
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _update_(self):
        """The module derives the 'derived_properties' based on the
        fundamental once
        """
        # -- Make sure the fundamental update (if any) are made
        super(Image,self)._update_()
        # - Data
        self._update_data_()

    def _update_data_(self):
        """
        """
        self._derived_properties["data"] = self.rawdata - self.background
        
