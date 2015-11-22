#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module defines the photometric objects"""


import numpy  as np
import pyfits as pf

from astropy     import units
from astLib      import astWCS

from .baseobject   import BaseObject 
from ..utils.tools import kwargs_update,flux_2_mag


__all__ = ["image","photopoint","lightcurve"]


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

def photopoint(lbda,flux,var,source=None,
               instrument_name=None,**kwargs):
    """
    Get the basic object containing  photometric point information
    
    Parameters
    ----------

    lbda: [float]              The central wavelength associated to the photometric
                               points.

    flux: [float]              The flux (*not magnitude*) of the photometric point.

    var: [float]               The variance associated to the point's flux.

    - option -

    empty: [bool]              Set True to return an empty object.

    - other options ; not exhautive ; goes to 'create' -

    source: [string]           Staten the origin of the point (e.g. image, ...)

    instrument_name:[string]   Give a name of the intrument that enable to take the
                               photometric point.

    Return
    ------
    PhotoPoint
    """
    return PhotoPoint(lbda,flux,var,source=source,
                      instrument_name=instrument_name,
                      **kwargs)


def lightcurve(datapoints,times,**kwargs):
    """
    This functions enables to create a lightcurve object
    from the given datapoints and times. The datapoints can
    either be a list of 'PhotoPoint' or a dictionnary having
    {flux, variance, lbda_s}

    Parameters
    ----------
    datapoints: [list/dict]        Two formes fo data can be given to this function.
                                   1) datapoints = list-of-PhotoPoint
                                   This is a simple list (or array) of astrobject
                                   PhotoPoint. Its size must corresponds to 'times' one
                                   2) a dictionnary having 3 parameters: lbda, fluxes,
                                   and variances. Each have to have the size of 'times'
                                   except lbda that could be simple float since all
                                   must share the same lbda.

    times: [array]                 This is the list/array of times corresponding to
                                   each data points.

    - options -

    Return
    ------
    LightCurve
    
    """
    if type(datapoints) is dict:
        # this must be a dictsource:
        datapoints = dictsource_2_photopoints(datapoints)
        
    if type(datapoints) is list or type(datapoints) is np.ndarray:
        # - This should be a list of PhotoPoints, otherwise create will break
        if len(times) != len(datapoints):
            raise ValueError("'times and datapoints must have the same size")
        return LightCurve(datapoints, times)
    
    raise TypeError("'datapoints' must be a list of photopoints or a dictsource")
        

def dictsource_2_photopoints(dictsource,**kwargs):
    """This fuctions enable to convert a dictionnary
    into a list of photopoints.
    This uses 'photopoint' to load the list.

    Parameters
    ----------
    dictsource: [dictionnary]      This dictionnary must have the following entries:
                                   {
                                   'flux': [list of fluxes],
                                   'var': [list of associated variances],
                                   'lbda': float-if-unique wavelength/[list of lbda]
                                   - options -
                                   'source': string/[list of sources],
                                   'instrument_name': string/[list of instruments],
                                   }
                                   The array must all have the same size.

                                   
    - kwargs options goes to all the looped PhotoPoint __init__-
    
    Return
    ------
    list-of-PhotoPoints
    """
    if type(dictsource) is not dict:
        raise TypeError("'dictsource' must be a dictionnary")
    
    # - Does is have the requiered basic information:
    for musthave in ["fluxes","variances","lbda"]:
        if musthave not in dictsource.keys():
            raise TypeError("'dictsource' must have the following entries: "+ \
                            "fluxes","variances","lbda")
    # - Good.
    # -------
    # - array-data
    fluxes,variances = dictsource["flux"],dictsource["var"]
    # - potentially float or None
    lbda = dictsource["lbda"]
    source = dictsource["source"] if "source" in dictsource.keys() else None
    instrument_name = dictsource["instrument_name"] \
      if "instrument_name" in dictsource.keys() else None
    
    lbda = [lbda]*len(fluxes) if hasattr(lbda,'__iter__') is False else lbda
    source = [source]*len(fluxes) if hasattr(source,'__iter__') is False else source
    instrument_name = [instrument_name]*len(fluxes)\
       if hasattr(instrument_name,'__iter__') is False else instrument_name

    # -- Let's create the list of PhotoPoints
    return [photopoint(lbda,flux,var,source=source,
                      instrument_name=instrument_name,
                      **kwargs)
            for flux,var,lbda,source,ins in zip(fluxes,variances,lbda,
                                                source,instrument_name)]

    
    
    
#######################################
#                                     #
# Base Object Classes: Image          #
#                                     #
#######################################
class Image( BaseObject ):
    """
    """
    __nature__ = "Image"
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
                 astrotarget=None,data_index=0,
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
        self.__build__(data_index=data_index)

        if empty:
            return
        
        if filename is not None:
            force_it = kwargs.pop("force_it",True)
            self.load(filename,force_it=force_it,
                      **kwargs)
        # - Set the target if any
        if astrotarget is not None:
            self.set_target(astrotarget)
            
    def __build__(self,data_index=0):
        #
        # Improvement of BaseObject
        # including the _object_properties
        #
        super(Image,self).__build__()
        # -- How to read the image
        self._build_properties = dict(
                data_index = data_index,
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
        self._properties["header"]       = pf.Header() if header is None else header
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

        background  = self._get_default_background_() if background is None \
           else background
        
        # Shape test
        if self.rawdata is not None and np.shape(background) != self.shape:
            raise ValueError("The given background must have rawdata's shape")
        # -- Looks good
        self._properties['background'] = np.asarray(background)
        self._update_data_()

    # ------------------- #
    # - get Methods     - #
    # ------------------- #
    # ----------- #
    #  Aperture   #
    # ----------- #
    def get_target_aperture(self,r_pixels,aptype="circle",subpix=5,**kwargs):
        """If a target is loaded, this will get the target coords and run
        'get_aperture'.
        """
        if self.target is None:
            raise AttributeError("No 'target' loaded")
        xpix,ypix = self.coords_to_pixel(self.target.ra,self.target.dec)
        return self.get_aperture(xpix,ypix,
                                 r_pixels=r_pixels,aptype=aptype,
                                 subpix=subpix,**kwargs)
    
    def get_aperture(self,x,y,r_pixels=None,aptype="circle",subpix=5,
                     ellipse_args={"a":None,"b":None,"theta":None},
                     annular_args={"rin":None,"rout":None},
                     **kwargs):
        """
        This method uses K. Barary's Sextractor python module SEP
        to measure the aperture photometry. See details here:
        sep.readthedocs.org/en/v0.4.x/apertures.html

        
        Parameters
        ----------

        x:[float]                  Pixel coordinate of the second ("fast") axis
                                   (so x in the in imshow). In pixels
                                   
        y:[float]                  Pixel coordinate of the first ("slow") axis
                                   (so y in the in imshow). In pixels
                                   
        (other aperture arguments are required be depend on the type of
        aperture photometry *aptype* is choosen)
                                   
        - options - 

        aptype: [string]           Type of Aperture photometry used.
                                   -circle  => set r_pixels
                                   -ellipse => set all ellipse_args entries
                                   -circan  => set all annulus_args entries
                                   -ellipan => set all ellipse_args entries
                                              set all annulus_args entries
                                    (no other type allowed)
                                               
        
        r_pixels: [float]          Size of the circle radius. In pixels.
                                   (This is used only if aptype is circle)

        subpix: [int]              Division of the real pixel to perform the
                                   circle to square overlap.

        ellipse_args: [dict]       The ellipse parameter that must be filled if
                                   atype is ellipse of ellipan.

        annular_args: [dict]       The annular parameter that must be filled if
                                   atype is circan of ellipan.

        
        - other options ; not exhautive ; goes to sep.sum_*aptype* - 

        gain: [float]              You can manually set the image gain. Otherwise
                                   this method will look for the key 'self._gain'
                                   and set None if it does not find it.

        var: [2d-array/float]      You can manually set the variance of the image.
                                   Otherwise, this will look for self.var. If this
                                   is None and and sepbackground has been created
                                   the sepbackground.rms()**2 will be used as a
                                   variance proxy. If not, None is set to var.
                                   
        See other sep options like: 'mask=None, maskthresh=0.0, bkgann=None...'
        
        
        Return
        ------
        sum, sumerr, flags (0 if no flag given)
        """
        import sep
        # - This should be moved in a c-lib
        if aptype not in ["circle","circann","ellipse","ellipan"]:
            raise ValueError("the given aptype (%s) is not a "+\
                             "known/implemeted sep aperture"%aptype)
                             
        # -------------
        # - SEP Input 
        gain = None if "_gain" not in dir(self) else self._gain
        gain = kwargs.pop("gain",gain)
        # -----------------
        # - Variance Trick
        var = self._sepbackground.rms()**2 if self.var is None \
          and "_sepbackground" in dir(self) else None
        var = kwargs.pop("var",var)

        # ------------------------
        # - Do the Aperture Photo
        
        # - Circle
        if aptype == "circle":
            return sep.sum_circle(self.data,x,y,r_pixels,subpix=subpix,
                            var=var,gain=gain,**kwargs)
        # - Annulus
        if aptype == "circann":
            if np.asarray([k is None for k in annular_args.values()]).any():
                raise ValueError("You must set the annular arguments 'annular_arg'")
            
            rin,rout = [annular_args[k] for k in ["rin","rout"]]
            return sep.sum_circann(self.data,x,y,rin,rout,
                                  r_pixels,subpix=subpix,
                                  var=var,gain=gain,**kwargs)
        # - Ellipse
        if aptype == "ellipse":
            if np.asarray([k is None for k in ellipse_args.values()]).any():
                raise ValueError("You must set the ellipse arguments 'ellipse_arg'")
            
            a,b,theta = [annular_args[k] for k in ["a","b","theta"]]
            return sep.sum_ellipse(self.data,x,y,a,b,theta,
                                    r_pixels,subpix=subpix,
                                    var=var,gain=gain,**kwargs)
        # - Elliptical Annulus
        if aptype == "ellipan":
            if np.asarray([k is None for k in ellipse_args.values()]).any():
                raise ValueError("You must set the ellipse arguments 'ellipse_arg'")
            if np.asarray([k is None for k in annular_args.values()]).any():
                raise ValueError("You must set the annular arguments 'annular_arg'")
            
            rin,rout = [annular_args[k] for k in ["rin","rout"]]
            a,b,theta = [ellipse_args[k] for k in ["a","b","theta"]]
            return sep.sum_ellipan(self.data,x,y,a,b,theta,rin,rout,
                                    r_pixels,subpix=subpix,
                                    var=var,gain=gain,**kwargs)
            

    # ------------------- #
    # - WCS Tools       - #
    # ------------------- #
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
        if "_sepbackground" not in dir(self):
            if self.rawdata is None:
                raise ValueError("no 'rawdata' loaded. Cannot get a background")
            from sep import Background
            self._sepbackground = Background(self.rawdata)

        return self._sepbackground.back()
        
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
    def _get_default_background_(self,*args,**kwargs):
        return self.get_sep_background(*args,**kwargs)
        
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
        
#######################################
#                                     #
# Base Object Classes: PhotoPoint     #
#                                     #
#######################################
class PhotoPoint( BaseObject ):
    """This Class hold the basic information associated to
    a photometric point"""

    __nature__ = "PhotoPoint"

    _properties_keys = ["lbda","flux","var"]
    _side_properties_keys = ["source","intrument_name"]
    _derived_properties_keys = []
    
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,lbda,flux,var,
                 empty=False,**kwargs):
        """
        Initialize the PhotoPoint object

        Parameters
        ----------

        lbda: [float]              The central wavelength associated to the photometric
                                   points.

        flux: [float]              The flux (*not magnitude*) of the photometric point.

        var: [float]               The variance associated to the point's flux.
        
        - option -

        empty: [bool]              Set True to return an empty object.
        
        - other options ; not exhautive ; goes to 'create' -

        source: [string]           Staten the origin of the point (e.g. image, ...)

        instrument_name:[string]   Give a name of the intrument that enable to take the
                                   photometric point.
        
        Return
        ------
        Void
        """
        self.__build__()
        if empty:
            return

        self.create(lbda,flux,var,**kwargs)    
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def create(self,lbda,flux,var,
               source=None,instrument_name=None,
               force_it=False):
        """This method create the object"""
        if self.flux is not None and not force_it:
            raise AttributeError("object is already defined."+\
                    " Set force_it to True if you really known what you are doing")
        # ****************** #
        # * Creation       * #
        # ****************** #
        self._properties["lbda"] = np.float(lbda)
        self._properties["flux"] = np.float(flux)
        self._properties["var"]  = np.float(var)

        self._side_properties["source"] = source
        self._side_properties["instrument_name"] = instrument_name
        self._update_()

    def display(self,ax,toshow="flux",**kwargs):
        """This method enable to display the current point
        in the given matplotlib axes"""
        # - Test if there is data
        if not self.has_data():
            raise AttributeError("no data to display")
        # -----------
        # - Input
        if toshow == "flux":
            y = self.flux
            dy= np.sqrt(self.var) if self.var is not None else None
        elif toshow == "mag":
            y = self.mag
            dy= np.sqrt(self.magvar) if self.magvar is not None else None
        else:
            raise ValueError("%s is not a known parameter"%toshow)
        # -----------
        # - Fancy
        default_prop = dict(marker="o",ecolor="0.7",
                        zorder=3)
        prop = kwargs_update(default_prop,**kwargs)
        # -----------
        # - Input
        pl = ax.errorbar(self.lbda,y,yerr=dy,**prop)
        self._plot = {}
        self._plot["ax"] = ax
        self._plot["marker"] = pl
        self._plot["prop"] = prop
        return self._plot
    
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def lbda(self):
        return self._properties["lbda"]

    @property
    def flux(self):
        return self._properties["flux"]

    @property
    def var(self):
        return self._properties["var"]

    def has_data(self):
        return not self.flux is None
    # ------------
    # - Side
    @property
    def source(self):
        return self._side_properties["source"]
    @source.setter
    def source(self,value):
        self._side_properties["source"] = value

    @property
    def instrument_name(self):
        return self._side_properties["instrument_name"]
    @instrument_name.setter
    def instrument_name(self,value):
        self._side_properties["instrument_name"] = value
        
    # ------------
    # - Derived 
    @property
    def mag(self):
        return flux_2_mag(self.flux,self.dflux,self.lbda)[0]
    @property
    def magvar(self):
        return flux_2_mag(self.flux,self.dflux,self.lbda)[1] ** 2



#######################################
#                                     #
# Base Object Classes: LightCurve     #
#                                     #
#######################################
def LightCurve( BaseObject ):
    """This object gather several PhotoPoint and there corresponding time
    to forme a lightcurve"""
    
    __nature__ = "LightCurve"

    _properties_keys = ["photopoints","times"]
    _side_properties_keys = []
    _derived_properties_keys = ["lbda"]
    
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,
                 photopoints=None,
                 times=None,
                 empty=False):
        """
        """
        self.__build__()
        if empty:
            return
        
        if photopoints is not None or times is not None:
            self.create(photopoints, times)
            
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def create(self, photopoints, times):
        """
        """
        if len(photopoints) != len(times):
            raise ValueError("photopoints and times must have the same length")
        
        # ********************* #
        # * Create the Object * #
        # ********************* #
        for p,t in zip(photopoints,times):
            self.add_photopoint(p,t)
            

    # ------------------ #
    # - I/O Photopoint - #
    # ------------------ #
    def create_new_photopoint(self,flux,var,lbda,time,
                              **kwargs):
        """
        """
        new_photopoint = photopoint(lbda,flux,var,**kwargs)
        self.add_photopoint(new_photopoint,time)
        
        
    def add_photopoint(self,photopoint,time):
        """
        """
        # ----------------------- #
        # - Test the bands      - #
        # ----------------------- #
        if "__nature__" not in photopoint or \
          photopoint.__nature__ != "PhotoPoint":
            raise TypeError("'photopoint' must be a list of astrobject PhotoPoint")
        
        self._test_lbda_(photopoint)
        
        if self.lbda is None:
            self._derived_properties['lbda'] = photopoint.lbda
            
        if self._properties['photopoints'] is None:
            self._properties['photopoints'] = []
            self._properties['times'] = []
            
            
        self._properties['photopoints'].append(photopoint)
        self._properties['times'].append(time)


    # ------------------ #
    # - Show           - #
    # ------------------ #
    def show(self,savefile=None,ax=None,
             inmag=False,show=True,**kwargs):
        """
        """
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

        pl = self.display(ax,inmag=inmag,**kwargs)
        self._plot['ax'] = ax
        self._plot['fig'] = fig
        self._plot['plot'] = pl
        self._plot['prop'] = kwargs
        fig.figout(savefile=savefile,show=show)
        
    def display(self,ax,inmag=False,**kwargs):
        """
        """
        if inmag:
            y,dy = self.mags,self.magsvar
        else:
            y,dy = self.fluxes,self.fluxesvar

        default_prop = dict(marker="o",ms=15,mec="k",mfc="b",
                            alpha=0.8,ecolor="0.7")
        prop = kwargs_update(default_prop,**kwargs)
        
        pl = ax.errorbar(self.times,y,yerr=dy,**prop)
        return pl
    
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _test_lbda_(self,photopoint):
        """This module check that the lbda in the given
        photopoint is the same as the other ones.
        Remark that is None won't be penalized. 
        """
        if self.lbda is None or photopoint.lbda is None:
            return
        if self.lbda != photopoint.lbda:
            raise ValueError("the given 'photopoint' does not share the other's lbda")
        
        
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def photopoints(self):
        return self._properties['photopoints']
    
    @property
    def times(self):
        return self._properties['times']
    
    @property
    def fluxes(self):
        return [p.flux for p in self.photopoints]
    
    @property
    def fluxesvar(self):
        return [p.var for p in self.photopoints]

    @property
    def mags(self):
        return [p.mag for p in self.photopoints]
    
    @property
    def magsvar(self):
        return [p.magvar for p in self.photopoints]
    
    @property
    def lbda(self):
        return self._derived_properties['lbda']
    
