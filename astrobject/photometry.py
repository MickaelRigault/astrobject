#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module defines the photometric objects"""


import numpy  as np
import pyfits as pf
from scipy.stats import sigmaclip

from astropy     import units,coordinates
from astLib      import astWCS

from .baseobject   import BaseObject 
from ..utils.tools import kwargs_update,flux_to_mag


__all__ = ["image","sexobjects",
           "photopoint","lightcurve","photomap"]


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

def photomap(fluxes=None,variances=None,ra=None,dec=None,lbda=None,
             photopoints=None,**kwargs):
    """
    Get the object gathering photometric and coordinate informations
    
    Parameters
    ----------

    flux: [array]              The fluxes (*not magnitude*) of the photometric points.

    variances: [array]         The variance associated to the point's fluxes.

    ra,dec: [array,array]      Coordinates *in degree* of the given points

    
    - option -

    lbda: [float]              The central wavelength associated to the photometric
                               points.

    photopoints: [array]       list of photopoints. Set that instead of fluxes / variances
                               and lbda. (ra and dec style requiered)
    
    empty: [bool]              Set True to return an empty object.

    
    - other options ; not exhautive ; goes to 'create' -


    Return
    ------
    PhotoMap
    """
    return PhotoMap(fluxes=fluxes,variances=variances,
                    ra=ra,dec=dec,lbda=lbda,
                    photopoints=photopoints,**kwargs)

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
        

def sexobjects(sexoutput):
    """
    """
    return SexObjects(sexoutput)
    
# ========================== #
#  Internal Tool             #
# ========================== #
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
    _side_properties_keys    = ["wcs",
                                "target","catalogue",
                                "exptime"] # maybe Exptime should just be on flight
        
    _derived_properties_keys = ["fits","data","sepobjects"]
    
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
        if self.rawdata is not None and force_it is False:
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
        if self.rawdata is not None and force_it is False:
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
        self._side_properties["exptime"] = exptime
          
        if wcs is not None and wcs.__module__ != "astLib.astWCS":
            print "WARNING: only astLib.astWCS wcs solution is implemented"
            print " ----> No wcs solution loaded"
            self._side_properties["wcs"]     = None
        else:
            self._side_properties["wcs"]     = wcs

        self._update_()


    # ------------------- #
    # - Set Methods     - #
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


    def set_catalogue(self,catalogue,force_it=False):
        """
        This Methods enables to load a Catalogue object
        """
        if self.has_catalogue() and force_it is False:
            raise AttributeError("'catalogue' already defined"+\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(catalogue) or catalogue.__nature__ != "Catalogue":
            raise TypeError("the input 'catalogue' must be an astrobject catalogue")
        # -------------------------
        # - Add the world_2_pixel
        if not self.has_wcs():
            print "WARNING: no wcs solution. catalogue recorded, but most likely is useless"
        else:
            # -- Lets only consider the objects in the FoV
            catalogue.set_wcs(self.wcs,force_it=True)
            # -- Lets save the pixel values
            #catalogue.x,catalogue.y = np.asarray([self.coords_to_pixel(ra_,dec_)
            #                                      for ra_,dec_ in zip(catalogue.ra,catalogue.dec)]).T
            if self.has_sepobjects():
                self.sepobjects.set_catalogue(catalogue,force_it=True)
                self.sepobjects.match_catalogue()
            
                
        # --------
        # - set it
        self._side_properties["catalogue"] = catalogue
        

    # ------------------- #
    # - download Methods- #
    # ------------------- #
    def download_catalogue(self,source="sdss",
                           set_it=True,force_it=False,**kwargs):
        """
        """
        from .instruments import instrument
        radec = "%s %s"%(self.wcs.getCentreWCSCoords()[0],self.wcs.getCentreWCSCoords()[1])
        radius = np.max(self.wcs.getHalfSizeDeg())*np.sqrt(2)
        cat = instrument.catalogue(source=source,radec=radec,radius="%sd"%radius,
                                    **kwargs)
        if not set_it:
            return cat
        
        self.set_catalogue(cat,force_it=force_it)
                    
    # ------------------- #
    # - get Methods     - #
    # ------------------- #
    # ----------- #
    #  Aperture   #
    # ----------- #
    def get_stars_aperture(self, r_pixels,aptype="circle",
                           isolated_only=True, catmag_range=[None,None],
                           **kwargs):
        """
        This methods fetch in the catalogue and sep_extracted objects (sepobjects)
        the matched points.
        It then perform a large aperture

        Return
        ------
        sep_idx, catalogue_idx, get_aperture's output (sum,sumerr,flag)
        """
        # -------------- #
        # - Input Test - #
        # -------------- #
        if not self.has_catalogue():
            raise AttributeError("no 'catalogue' set. This is required")
        
        if not self.has_sepobjects():
            raise AttributeError("no 'sepobjects'.  This is required")

        if isolated_only and not self.catalogue.around_defined():
            raise AttributeError("no 'around' not defined in the catalogue."+ \
                                 " Run catalogue.defined_around or set isolated_only to False")
            
        # -------------------- #
        # - Pairing          - #
        # -------------------- #
        # note that might be a bit slow if called for the first time
        kwardsmask = dict(stars_only=True,isolated_only=isolated_only,
                          catmag_range=catmag_range)
        idx = self.sepobjects.get_indexes(**kwardsmask)
        cat_idx = self.sepobjects.get_indexes(cat_indexes=True,**kwardsmask)
        
        x,y = self.catalogue.wcs_xy.T[cat_idx].T
        
        # -------------------- #
        # - Get it           - #
        # -------------------- #
        return idx, cat_idx, self.get_aperture(x,y,r_pixels=r_pixels,aptype=aptype,**kwargs)
    
    
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
        var = kwargs.pop("var",self.var)

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
        
    def sep_extract(self,thresh=None,returnobjects=False,
                    set_catalogue=True,match_catalogue=True,**kwargs):
        """
        This module is based on K. Barbary's python module of Sextractor SEP.

        Parameters
        ----------

        - options -

        thresh: [float]            Threshold pixel value for detection.
                                   If None is set, the globalrms background from sep
                                   will be used.
                                   Additional information from sep.extract:
                                   "If an err array is not given, this is interpreted
                                    as an absolute threshold. If err is given,
                                    this is interpreted as a relative threshold:
                                    the absolute threshold at pixel (j, i) will be
                                    thresh * err[j, i]."

        set_catalogue: [bool]      If the current instance has a catalogue, it will be transfered
                                   to the SexOutput object created. Set False to avoid that.
                                   
                                   
        returnobjects: [bool]      Change the output of this function. if True the
                                   extracted objects are recorded and returned (self.sepobjects)
                                   if not they are just recorded.

        
        - others options -

        kwargs                     goes to set.extract
                                   (sep.readthedocs.org/en/v0.5.x/api/sep.extract.html)
        
        Return
        -------
        Void [or ndarray(sep.extract output) is returnobjects set to True]
        """
        from sep import extract
        if thresh is None:
            thresh = self._get_sep_extract_threshold_()
            
        o = extract(self.data,thresh,**kwargs)
        self._derived_properties["sepobjects"] = sexobjects(o)

        if self.has_wcs():
            self.sepobjects.set_wcs(self.wcs)

        if set_catalogue and self.has_catalogue():
            self.sepobjects.set_catalogue(self.catalogue,force_it=True)
            if match_catalogue:
                self.sepobjects.match_catalogue()
                    
        if returnobjects:
            return self.sepobjects

    def _get_sep_extract_threshold_(self):
        """this will be used as a default threshold for sep_extract"""
        if "_sepbackground" not in dir(self):
                _ = self.get_sep_background()
        return self._sepbackground.globalrms*1.5 

    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def writeto(self):
        """
        """
        print "to be done"


    # ------------------- #
    # - Plot Methods    - #
    # ------------------- #        
    def show(self,toshow="data",savefile=None,logscale=True,
             ax=None,show=True,wcs_coords=False,
             zoomon=None,zoompxl=200,
             show_sepobjects=True,propsep={},
             show_catalogue=True,
             proptarget={},
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
        # - Do It     #
        im = ax.imshow(x,**prop)
        # - add target
        pl_tgt = None if self.has_target() is False \
          else self.display_target(ax,wcs_coords=wcs_coords,
                                   **proptarget)

        if show_sepobjects and self.has_sepobjects():
            self.sepobjects.display(ax,world_coords=wcs_coords,
                                    **propsep)
        if show_catalogue and self.has_catalogue():
            self.display_catalogue(ax,wcs_coords=wcs_coords)
        # ----------- #
        # - Zoom      #
        if zoomon is not None:
            # -- Zoom on target
            if type(zoomon) is str and zoomon=="target" and self.has_target():
                    coords_zoom = self.coords_to_pixel(*self.target.radec) if not wcs_coords \
                    else self.target.radec
            elif np.shape(zoomon) == (2,):
                coords_zoom = zoomon
            else:
                print "WARNING can not parse to zoom on input"
                coords_zoom = None
        else:
            coords_zoom = None

        if coords_zoom is not None:
            width = zoompxl if not wcs_coords else zoompxl*self.pixel_size_deg
            ax.set_xlim(coords_zoom[0]-width,coords_zoom[0]+width)
            ax.set_ylim(coords_zoom[1]-width,coords_zoom[1]+width)
            
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

    def display_target(self,ax,wcs_coords=True,draw=True,**kwargs):
        """If a target is loaded, use this to display the target on the
        given ax"""
        if self.has_target() is False:
            print "No 'target' to display"
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

        if draw:
            ax.figure.canvas.draw()        
        return pl

    def display_sepobjects(self,ax=None,wcs_coords=True,draw=True,**kwargs):
        """If sep_extract has been ran, you have an sepobjects entry.
        This entry will be red and parsed here.
        """
        # --------------
        # - ax parsing
        if ax is None and ("_plot" not in dir(self) or "ax" not in self._plot.keys()):
            raise ValueError('no ax defined')
        ax = self._plot['ax'] if ax is None else ax
        
        self.sepobjects.display(ax,world_coords=wcs_coords,draw=draw,
                                **kwargs)

    def display_catalogue(self,ax=None,wcs_coords=True,draw=True,**kwargs):
        """This methods enable to show all the known sources in the
        image's field of view.
        This will only works if a catalogue has been set"""
        if not self.has_catalogue():
            print "No 'catalogue' to display"
            return

        # --------------
        # - ax parsing
        if ax is None and ("_plot" not in dir(self) or "ax" not in self._plot.keys()):
            raise ValueError('no ax defined')
        ax = self._plot['ax'] if ax is None else ax
        
        # --------------------------------
        # - Draw the matched-points
        return self.catalogue.display(ax,wcs_coords=wcs_coords,draw=draw,
                                      scalex=False,scaley=False,
                                      **kwargs)

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
        if self._side_properties['exptime'] is None:
            # -- It has not be set manually, maybe check the header
            if self._build_properties["header_exptime"] in self.header:
                self._side_properties['exptime'] = \
                  np.float(self.header[self._build_properties["header_exptime"]])
        # -- You have it ? This will stay None if not
        return self._side_properties['exptime']
    # ----------------------
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
        if self._properties["var"] is None:
           self._properties["var"] = self._get_default_variance_()
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

          
    # ----------------      
    # -- SEP OUTPUT
    @property
    def sepobjects(self):
        return self._derived_properties["sepobjects"]

    def has_sepobjects(self):
        return True if self.sepobjects is not None and self.sepobjects.has_data() \
          else False
    
    # ----------------      
    # -- CATALOGUE
    @property
    def catalogue(self):
        return self._side_properties["catalogue"]
    
    def has_catalogue(self):
        return False if self.catalogue is None\
          else True
        
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _get_default_variance_(self):
        """
        """
        if "_sepbackground" in dir(self):
            return self._sepbackground.rms()**2
        return None
        
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
    _side_properties_keys = ["source","intrument_name","target"]
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

    def set_target(self,newtarget):
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
        
        # -- Seems Ok -- #
        self._side_properties["target"] = newtarget.copy()
        
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
    
    # -- source
    @property
    def source(self):
        return self._side_properties["source"]
    @source.setter
    def source(self,value):
        self._side_properties["source"] = value

    # -- instrument
    @property
    def instrument_name(self):
        return self._side_properties["instrument_name"]
    
    @instrument_name.setter
    def instrument_name(self,value):
        self._side_properties["instrument_name"] = value
        
    # -- Target
    @property
    def target(self):
        return self._side_properties['target']

    def has_target(self):
        return False if self.target is None \
          else True
          
    # ------------
    # - Derived 
    @property
    def mag(self):
        return flux_to_mag(self.flux,np.sqrt(self.var),self.lbda)[0]
    @property
    def magvar(self):
        return flux_to_mag(self.flux,np.sqrt(self.var),self.lbda)[1] ** 2



#######################################
#                                     #
# Base Object Classes: LightCurve     #
#                                     #
#######################################
class LightCurve( BaseObject ):
    """This object gather several PhotoPoint and there corresponding time
    to forme a lightcurve"""
    
    __nature__ = "LightCurve"

    _properties_keys = ["photopoints","times"]
    _side_properties_keys = []
    _derived_properties_keys = ["lbda"]
    
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,photopoints=None,
                 times=None,empty=False):
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
        if "__nature__" not in dir(photopoint) or \
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



#######################################
#                                     #
# Base Object Classes: PhotoMap       #
#                                     #
#######################################
class PhotoMap( BaseObject ):
    """
    """
    __nature__ = "PhotoMap"
    
    _properties_keys = ["fluxes","variances","ra","dec","lbda"]
    _side_properties_keys = ["refmap", "wcs", "sep_params"]
    _derived_properties_keys = [""]
    
    def __init__(self,fluxes=None,variances=None,
                 ra=None,dec=None,lbda=None,
                 photopoints=None,
                 empty=False):
        """
        """
        self.__build__()
        if empty:
            return

        if photopoints is not None:
            self.create_from_photopoints(photopoints,ra,dec)
            return
        
        if fluxes is not None and variances is not None:
            self.create(fluxes,variances,ra,dec,lbda)
            
            
        
    # ========================== #
    # = Init Methods           = #
    # ========================== #
    def create(self, fluxes, variances, ra, dec,lbda):
        """
        This method is the core function that create the object.
        All the other initializing method rely on this one.
        
        
        Parameters
        ----------

        fluxes: [array]            The fluxes of each points

        variances: [array]         Associated variances

        ra,dec: [array,array]      Poisiton in the sky *in degree*

        lbda: [None, float]        wavelength associated to the points

        Return
        ------
        void
        """
        
        if len(fluxes) != len(variances):
            raise ValueError("fluxes and variances must have the same size")

        if len(fluxes) != len(ra):
            raise ValueError("fluxes and ra must have the same size")

        if len(fluxes) != len(dec):
            raise ValueError("ra and dec must have the same size")

        if "__iter__" in dir(lbda):
            raise TypeError("lbda must be a float or an int. not an array")

        # ********************* #
        # * Creation          * #
        # ********************* #
        self._properties["lbda"] = np.float(lbda) if lbda is not None else None
        self._properties["fluxes"] = np.asarray(fluxes,dtype= float)
        self._properties["variances"] = np.asarray(variances,dtype= float)
        self._properties["ra"] = np.asarray(ra,dtype= float)
        self._properties["dec"] = np.asarray(dec,dtype= float)
        
        
    def create_from_photopoints(self,photopoints,
                                ra,dec):
        """
        This methods enable to extract, from a list of photpoints
        the requiered input for this instance. Remark that if
        the photpoints have a target set and if this target has
        ra and dec defined, no need to inpout radec. Otherwisem you do.

        Parameters
        ----------

        photopoints: [array of PhotoPoints]
                                   a list of photopoints that contains
                                   the flux/error/lbda [target] parameters

        Return
        ------
        Void
        """
        fluxes, variances, lbda = np.asarray([ [p.flux, p.var, p.lbda]
                                                for p in photopoints ]).T
        # -- ra,dec
        ra,dec = radec if radec is not None else None,None
        if ra is not None and len(ra) != len(flux):
            raise ValueError("radec size must be equal to flux's one")
        if ra is None:
            try:
                ra,dec = np.asarray([ [p.target.ra,p.target.dec] 
                                 for p in photopoints]).T
            except:
                raise AttributeError("no radec provided and no photopoint.target to find it.")
        
        # ------------ #
        # call create
        # ------------ #
        self.create(fluxes, variances, ra, dec, lbda)

        
    def writeto(self,savefile):
        """
        """
        print "to be done"
                
    def load(self,filename):
        """
        """
        print "to be done"
    # ========================== #
    # = Get Methods            = #
    # ========================== #
    def get_photopoint(self,index):
        """
        """
        print "to be done"
        
    def get_index_around(self,ra, dec, radius):
        """
        Parameters:
        ----------
        ra, dec : [float]          Position in the sky *in degree*

        radius: [float with unit]  distance of search, using astropy.units

        Return
        ------
        (SkyCoords.search_around_sky output)
        """
        sky = coordinates.SkyCoord(ra=ra*units.degree, dec=dec*units.degree)
        return self.sky_radec.search_around_sky(sky, radius)
    
    # ========================== #
    # = Set Methods            = #
    # ========================== #
    def set_wcs(self,wcs,force_it=False):
        """
        """
        if self.has_wcs() and force_it is False:
            raise AttributeError("'wcs' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        if wcs is not None and "coordsAreInImage" not in dir(wcs):
            raise TypeError("'wcs' solution not recognize, should be an astLib.astWCS")

        self._side_properties["wcs"] = wcs

    def set_refmap(self, photomap,force_it=False):
        """
        Set here the reference map associated to the current one.
        The ordering of the point must be consistant with that of the current
        photomap. 
        """
        if self.has_refmap() and force_it is False:
            raise AttributeError("'refmap' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(photomap) or photomap.__nature__ != "PhotoMap":
            raise TypeError("The given reference map must be a astrobject PhotoMap")

        self._side_properties["refmap"] = photomap


    def set_sep_params(self, sep_params,force_it=True):
        """
        sep_param is a dictionnary containing some SEP info. Setting it here
        will enable to use the associated function (e.g. some 'get' keyword)
        """
        if self.has_sep_params() and force_it is False:
            raise AttributeError("'sep_params' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        self._side_properties["sep_params"] = sep_params
        
    def get(self,key):
        """
        """
        if key is None:
            return None
        _sep_keys_ = self.sep_params.keys()
        _derived_keys_ = ["flux_ratio","scaled_flux_ratio"]
        help_text = " Known keys are: "+", ".join(_sep_keys_+_derived_keys_)
        # ---------------
        # - key parsing
        
        if key in ["help","keys","keylist"]:
            print help_text
            return
        
        if key in _sep_keys_:
            return self.sep_params[key]
        
        if key in _derived_keys_:
            if key == "flux_ratio":
                return self.fluxes / self.refmap.fluxes
            if key == "scaled_flux_ratio":
                ratio = self.get("flux_ratio")
                return ratio - np.median(ratio)
            
        raise NotImplementedError("'%s' access is not implemented"+help_text)
    
    # ========================== #
    # = Plot Methods           = #
    # ========================== #
    def display_voronoi(self,ax,toshow="scaled_flux_ratio",wcs_coords=False,**kwargs):
        """
        This methods enables to display in the given axes patchs following the voronoi
        tesselation based on the photomap radec/xy (see wcs_coords) parameters.
        This methods won't scale the ax according to the voronoi tesselation.

        Parameters
        ----------

        ax: [matplotlib.Axes]      where the voronoi patch collection should be drawn
        
        - options -

        toshow: [string]           keyword set to self.get that parses it. Use 'help'
                                   for additional information.
                                   The value of this key will define the color of the
                                   patchs.
        
        - kwargs goes to mpladdon.voronoi_patchs -

        non-exhaustive list of key arguments: *vmin, vmax, cmap,
        cbar [bool or ax], cblabel, any matplotlib.PolyCollection keys*
        
        Return
        ------
        PatchCollection (see mpladdon.voronoi_patchs)
        """
        from ..utils.mpladdon import voronoi_patchs
        
        # -----------------
        # - What to show
        if wcs_coords:
            x,y = self.ra,self.dec
        else:
            xy = np.asarray(self.wcs_xy)
            
        # -----------------
        # - The plot itself
        out = ax.voronoi_patchs(xy,self.get(toshow),**kwargs)
        ax.figure.canvas.draw()
        return out
        
    
    # =========================== #
    # Properties and Settings     #
    # =========================== #
    @property
    def npoints(self):
        if self.ra is None:
            return 0
        return len(self.ra)
    
    # ----------------
    # - PhotoPoints
    @property
    def lbda(self):
        return self._properties["lbda"]
    
    @property
    def fluxes(self):
        return self._properties["fluxes"]

    @property
    def variances(self):
        return self._properties["variances"]

    @property
    def _mag_magerr(self):
        return flux_to_mag(self.fluxes,np.sqrt(self.variances), self.lbda)
    
    @property
    def mag(self):
        return self._mag_magerr[0]
    @property
    def magerr(self):
        return self._mag_magerr[1]
    # --------------
    # - radec
    @property
    def ra(self):
        return self._properties["ra"]
    
    @property
    def dec(self):
        return self._properties["dec"]

    @property
    def sky_radec(self):
        """This is an advanced radec methods tight to astropy SkyCoords"""
        ra,dec = self.ra,self.dec
        return coordinates.SkyCoord(ra=ra*units.degree,dec=dec*units.degree)
    
    # --------------- #
    # - Side Prop   - #
    # --------------- #
    
    # ----------------      
    # -- WCS
    @property
    def wcs(self):
        return self._side_properties['wcs']

    def has_wcs(self):
        return False if self.wcs is None \
          else True
          
    @property
    def wcs_xy(self):
        if not self.has_wcs():
            raise AttributeError("no 'wcs' solution defined")
        
        return self.wcs.wcs2pix(self.ra,self.dec)
        
    # ----------------      
    # -- References
    @property
    def refmap(self):
        return self._side_properties["refmap"]
    
    
    def has_refmap(self):
        return False if self.refmap is None\
          else True

    # ----------------      
    # -- sep Parameters
    @property
    def sep_params(self):
        if self._side_properties["sep_params"] is None:
            self._side_properties["sep_params"] = {}
        return self._side_properties["sep_params"]
    
    def has_sep_params(self):
        return False if len(self.sep_params.keys())==0 \
          else True




#######################################
#                                     #
# Base Object Classes: SEPObjects     #
#                                     #
#######################################
class SexObjects( BaseObject ):
    """This instance parse the ourput from Sextractor/SEP and have
    convinient associated functions"""

    _properties_keys = ["data"]
    _side_properties_keys = ["catalogue","wcs"]
    _derived_properties_keys = ["catmatch","starmask","pairdict"]
    
    
    def __init__(self,sexoutput=None,wcs=None,empty=False):
        """
        = Load the SexObject =

        Parameters
        ----------

        sexoutput: [ndarray/file]  The ndarray as output by sextractor/sep or the
                                   associated file

        - options -

        empty: [bool]              return an empy object
        """
        self.__build__()
        if empty:
            return

        self.create(sexoutput)
        if wcs is not None:
            self.set_wcs(wcs)
    # ====================== #
    # Main Methods           #
    # ====================== #
    def create(self,sexoutput,force_it=False):
        """
        """
        if self.has_data() and force_it is False:
            raise AttributeError("'data' is already defined."+\
                    " Set force_it to True if you really known what you are doing")
                    
        sexdata = self._read_sexoutput_input_(sexoutput)
        if sexdata is None:
            print "WARNING empty imput data. Empty object loaded"
            return
        
        # ****************** #
        # * Creation       * #
        # ****************** #
        self._properties["data"] = sexdata

    # ------------------ #
    # - SET Methods    - #
    # ------------------ #
    def set_wcs(self,wcs,force_it=False):
        """
        """
        if self.has_wcs() and force_it is False:
            raise AttributeError("'wcs' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        if wcs is not None and "coordsAreInImage" not in dir(wcs):
            raise TypeError("'wcs' solution not recognize, should be an astLib.astWCS")

        self._side_properties["wcs"] = wcs
        

    def set_catalogue(self,catalogue,force_it=True,
                      default_isolation_def = 5*units.arcsec):
        """
        Parameters
        ---------

        default_isolation_def: [ang dist]
                                   If 'define_around' has not be ran yet,
                                   this scale will be used.
                                   This is important to know which object is
                                   isolated.
                                   
                                   
        Return
        ------
        
        """
        if self.has_catalogue() and force_it is False:
            raise AttributeError("'catalogue' already defined"+\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(catalogue) or catalogue.__nature__ != "Catalogue":
            raise TypeError("the input 'catalogue' must be an astrobject catalogue")
        # -------------------------
        # - Add the world_2_pixel
        if not catalogue.has_wcs():
            print "WARNING the given 'catalogue' has no pixel coordinates. Cannot load it"
            return
        
        if not catalogue.around_defined():
            catalogue.define_around(default_isolation_def)
            
        self._side_properties["catalogue"] = catalogue

    def match_catalogue(self,catalogue=None,force_it=False,arcsec_size=2):
        """This methods enable to attached a given sexobject entry
        to a catalog value.
        You can set a catalogue.
        """
        # --------------
        # - input 
        if catalogue is not None:
            self.set_catalogue(catalogue,force_it=force_it)

        if not self.has_catalogue():
            raise AttributeError("No 'catalogue' defined or given.")
        
        if self.has_wcs():
            # -- matching are made in degree space
            idxcatalogue, idxsexobjects, d2d, d3d = self.sky_radec.search_around_sky(self.catalogue.sky_radec, arcsec_size*units.arcsec)
        else:
            raise NotImplementedError("You currently need a wcs solution in the SexObjects to match a catalogue")
        # --------------------
        # - Save the results
        self._derived_properties["catmatch"] = {
            "idx_catalogue":idxcatalogue,
            "idx":idxsexobjects,
            "angsep":d2d
            }
        
        self.catalogue.set_matchedmask(idxcatalogue)
        self._derived_properties["starmask"] = np.asarray([i in self.catmask[self.catstarmask]
                                                for i in range(len(self.data["x"]))],dtype=bool)
        
    # ------------------ #
    # - get Methods    - #
    # ------------------ #
    def get(self,key):
        """This function enable to get from the data the values of the given keys
        or derived values, like ellipticity. Set 'help' for help."""
        if not self.has_data():
            raise AttributeError("no 'data' defined")

        # -- These are the default key values
        _data_keys_ = self.data.dtype.fields.keys()
        _matching_keys_ = ["angsep"]
        _derived_keys_ = ["elongation","ellipticity"]
        help_text = " Known keys are: "+", ".join(_data_keys_+_matching_keys_+_derived_keys_)
        if key in ["help","keys","keylist"]:
            print help_text
            return
        # -- These are from the data
        if key in _data_keys_:
            return self.data[key]
        
        # -- These are the catalogue values
        if key in _matching_keys_:
            if not self.has_catmatch():
                raise AttributeError("no 'catmatch' defined. The matching has not been ran")
            return self.catmatch[key]

        # -- These are derived values
        if key in _derived_keys_:
            if key == "elongation":
                return self.get("a") / self.get("b")
            
            if key == "ellipticity":
                return 1. - 1. / self.get("elongation")

        raise ValueError("Cannot parse '%s'."%key +\
                          help_text)

    def get_median_ellipse(self,apply_catmask=True,
                            stars_only=True, isolated_only=True,
                            catmag_range=[None,None],clipping=[3,3]):
        
        """This methods look for the stars and return the mean ellipse parameters"""
        if apply_catmask:
            idx = self.get_indexes(isolated_only=isolated_only,stars_only=stars_only,
                                    catmag_range=catmag_range)
            
        # -- add input test
         # -- apply the masking
        a_clipped,_alow,_ahigh = sigmaclip(self.get("a")[idx],*clipping)
        b_clipped,_blow,_bhigh = sigmaclip(self.get("b")[idx],*clipping)
        t_clipped,_tlow,_thigh = sigmaclip(self.get("theta")[idx],*clipping)
        # - so        
        psf_a,psf_b,psf_t = a_clipped.mean(),b_clipped.mean(),t_clipped.mean()
        m = np.sqrt(len(a_clipped)-1)
        
        return [psf_a,np.std(a_clipped)/m],[psf_b,np.std(t_clipped)/m],\
        [psf_t,np.std(t_clipped)/m]
        
        
    def get_pairing(self, stars_only=True, isolated_only=False,
                    catmag_range=[None,None]):
        """
        """
        pair = {}
        
        for i,d in self.pairdict.items():
            if stars_only and not d['is_star']:
                continue
            if isolated_only and not d['is_isolated']:
                continue
            if catmag_range[0] is not None and d["cat_mag"]<catmag_range[0]:
                continue
            if catmag_range[1] is not None and d["cat_mag"]>catmag_range[1]:
                continue
            
            pair[i] = d
            
        return pair

    def idx_to_mask(self,idx):
        mask = np.zeros(self.nobjects,dtype=bool)
        for i in idx:
            mask[i] = True
        return mask

    # ---------------- #
    # - get Mask     - #
    # ---------------- #
    def get_mask(self,isolated_only,stars_only,
                 catmag_range=[None,None]):
        """
        """
        mask = np.asarray(self.get_catmag_mask(*catmag_range))
        if isolated_only:
            mask = mask & self.catisolatedmask
        if stars_only:
            mask = mask & self.catstarmask
        
        return mask

    def get_indexes(self,isolated_only=False,stars_only=False,
                    catmag_range=[None,None], cat_indexes=False):
        """
        """
        id = "idx_catalogue" if cat_indexes else "idx"
        return self.catmatch[id][self.get_mask(isolated_only=isolated_only,
                                                  stars_only=stars_only,
                                                  catmag_range=catmag_range
                                                  )]
    
                 
    def get_catmag_mask(self,magmin,magmax):
        """
        return the boolen mask of which matched point
        belong to the given magnitude range.
        Set None for no limit
        Return
        ------
        bool mak array
        """
        if not self.has_catalogue():
            raise AttributeError("no 'catalogue' loaded")
        if not self.has_catmatch():
            raise AttributeError("catalogue has not been matched to the data")
    
        mags = self.catalogue.mag[self.catmatch["idx_catalogue"]]
        magmin = np.min(mags) if magmin is None else magmin
        magmax = np.max(mags) if magmax is None else magmax
        return (mags>=magmin) & (mags<=magmax)
        

    def get_photomap(self, matched_only=True,
                     stars_only=False, isolated_only=False):
        """
        """
        print "to be done"
    # ------------------- #
    # - PLOT Methods    - #
    # ------------------- #
    # - Display
    def display(self,ax,world_coords=True,draw=True,
                apply_catmask=True,
                stars_only=False, isolated_only=False,
                catmag_range=[None,None],
                **kwargs):
        """
        """
        if not self.has_data():
            print "WARNING [Sexobjects] No data to display"
            return
        
        from matplotlib.patches import Ellipse
        
        # -- maskout non matched one if requested
        if apply_catmask:
            kwargsmask = dict(isolated_only=isolated_only,stars_only=stars_only,
                              catmag_range=catmag_range)
            mask = self.get_indexes(**kwargsmask)
        
        x,y = self.data["x"][mask],self.data["y"][mask]
        # -------------
        # - Properties
        default_prop = dict(ls="None",marker="o",mec="k",mfc="0.7",alpha=0.5,
                            scalex=False,scaley=False,
                            label="_no_legend_")
        prop = kwargs_update(default_prop,**kwargs)

        
        
        ells = [Ellipse([x,y],a*5,b*5,t*units.radian.in_units("degree"))
                for x,y,a,b,t in zip(self.data["x"][mask],self.data["y"][mask],
                                     self.data["a"][mask],self.data["b"][mask],
                                     self.data["theta"][mask])]
        
        #pl = ax.plot(x,y,**prop)
        for ell in ells:
            ell.set_clip_box(ax.bbox)
            ell.set_facecolor("None")
            ell.set_edgecolor("k")
            ax.add_patch(ell)
        if draw:
            ax.figure.canvas.draw()
        #return pl

    # - Histograms
    def show_hist(self,toshow="a",ax=None,
                  savefile=None,show=True,
                  apply_catmask=True,
                  stars_only=False, isolated_only=False,
                  catmag_range=[None,None],
                  **kwargs):
        """This methods enable to show the histogram of any given
        key."""
        # -- Properties -- #
        if apply_catmask:
            kwargsmask = dict(isolated_only=isolated_only,stars_only=stars_only,
                              catmag_range=catmag_range)
            mask = self.get_indexes(**kwargsmask)
        
        v = self.get(toshow)[mask]

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
        default_prop = dict(histtype="step",fill=True,
                            fc=mpl.cm.Blues(0.7,0.5),ec="k",lw=2)
        prop = kwargs_update(default_prop,**kwargs)

        out = ax.hist(v,**prop)

        self._plot['ax'] = ax
        self._plot['fig'] = fig
        self._plot['hist'] = out
        self._plot['prop'] = kwargs
        fig.figout(savefile=savefile,show=show)

    # - Ellipse
    def show_ellipses(self,ax=None,
                      savefile=None,show=True,
                      apply_catmask=True,
                      stars_only=False, isolated_only=False,
                      catmag_range=[None,None],
                      **kwargs):
        """
        """
        if not self.has_data():
            print "WARNING [Sexobjects] No data to display"
            return
        
        from matplotlib.patches import Ellipse,Polygon
        from ..utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        # ------------------- #
        # - axes            - #
        # ------------------- #
        if ax is None:
            fig = mpl.figure(figsize=[6,6])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        elif "hist" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure

        # ------------------- #
        # - axes            - #
        # ------------------- #
        if apply_catmask:
            kwargsmask = dict(isolated_only=isolated_only,stars_only=stars_only,
                              catmag_range=catmag_range)
            mask = self.get_indexes(**kwargsmask)
        # -------------
        # - Properties
        
        
        ells = [Ellipse([0,0],2.,2*b/a,t*units.radian.in_units("degree"))
                for a,b,t in zip(self.data["a"][mask],self.data["b"][mask],
                                 self.data["theta"][mask])]
        # -- Show the typical angle
        psf_a,psf_b,psf_theta = self.get_median_ellipse(apply_catmask=apply_catmask,
                                                       **kwargsmask)
        ellipticity = 1- psf_b[0]/psf_a[0]
        # - cos/ sin what angle in radian
        
        ax.plot([0,np.cos(psf_theta[0])*ellipticity],[0,np.sin(psf_theta[0])*ellipticity],ls="-",lw=2,
                 color=mpl.cm.Blues(0.8),zorder=8)

        Cone_error = Polygon( [ [0,0],[np.cos(psf_theta[0]-psf_theta[1])*ellipticity,
                                       np.sin(psf_theta[0]-psf_theta[1])*ellipticity],
                                [np.cos(psf_theta[0]+psf_theta[1])*ellipticity,
                                 np.sin(psf_theta[0]+psf_theta[1])*ellipticity],
                                 [0,0]]
                                )
        Cone_error.set_facecolor(mpl.cm.Blues(0.8))
        Cone_error.set_edgecolor(mpl.cm.Blues(0.8))
        Cone_error.set_alpha(0.3)
        ax.add_patch(Cone_error)

        # -- Show the Ellipses
        for ell in ells:
            ell.set_clip_box(ax.bbox)
            ell.set_facecolor("None")
            ell.set_edgecolor("k")
            ell.set_alpha(0.1)
            ax.add_patch(ell)

        ax.set_ylim(-1.1,1.1)
        ax.set_xlim(-1.1,1.1)
        # -- show the center
        ax.axvline(0,ls="--",color="k",alpha=0.2)
        ax.axhline(0,ls="--",color="k",alpha=0.2)
        
        self._plot['ax'] = ax
        self._plot['fig'] = fig
        self._plot['prop'] = kwargs
        fig.figout(savefile=savefile,show=show)
        
    # =========================== #
    # Internal Methods            #
    # =========================== #
    def _read_sexoutput_input_(self,sexoutput):
        """This method enable to parse the input and therefore
        allow flexible inputs"""

        # ---------------
        # -- no input
        if sexoutput is None:
            return
        
        # ---------------
        # -- data input
        if type(sexoutput) == np.ndarray:
            # let's check
            try:
                _ = sexoutput["cxx"]
            except :
                raise TypeError("the given 'sexoutput' ndarray has no 'cxx' key."+"\n"+\
                                " It most likely is not a sextrator/sep output ")
            return sexoutput
        
        # ---------------
        # -- file input
        if type(sexoutput) == str:
            if sexoutput.endswith(".pkl"):
                from ..utils.tools import load_pkl
                self._read_sexoutput_input_(load_pkl(sexoutput))
                
            # -- So this is a sextrator output ?
            return

    def _get_pairing_(self):
        """
        This is slow, should be optimized
        """
        if not self.has_catmatch():
            raise AttributeError("no matching performed")
        
        idxself,idxcat = self.catmatch["idx"],self.catmatch["idx_catalogue"]
        pair = {}
        for i,icat in zip(idxself,idxcat):
            if i in pair.keys():
                pair[i]["other_cat_match"].append(icat)
                continue
            
            pair[i] = {"cat_idx":icat,
                       "cat_mag":self.catalogue.mag[icat],
                       "cat_magerr":self.catalogue.mag_err[icat],
                       #"cat_ra":self.catalogue.ra[icat],
                       #"cat_dec":self.catalogue.dec[icat],
                       "is_isolated":self.catalogue.isolatedmask[icat] if self.catalogue.around_defined() else None,
                       "is_star":self.catalogue.starmask[icat],
                       #"sep_a":self.get("a")[i],
                       #"sep_b":self.get("b")[i],
                       "sep_x":self.get("x")[i],
                       "sep_y":self.get("y")[i],
                       #"sep_theta":self.get("theta")[i],
                       #"sep_flux":self.get("flux")[i],
                       "other_cat_match":[]
                       }
            
        return pair
    # =========================== #
    # Properties and Settings     #
    # =========================== #
    # ----------------      
    # -- Data
    @property
    def data(self):
        return self._properties["data"]
    
    def has_data(self):
        return False if self.data is None\
          else True
          
    @property
    def nobjects(self):
        if not self.has_data():
            return None
        return len(self.data)

    
    # ----------------
    # - WCS
    @property
    def wcs(self):
        return self._side_properties['wcs']

    def has_wcs(self):
        return False if self.wcs is None else True
    
    @property
    def radec(self):
        if not self.has_wcs():
            raise AttributeError("no 'wcs' solution avialable. Cannot have radec")
        return np.asarray([self.wcs.pix2wcs(x_,y_) for x_,y_ in zip(self.data["x"],self.data["y"])]).T

    @property
    def xy(self):
        return np.asarray([self.data["x"],self.data["y"]])

    @property
    def sky_radec(self):
        """This is an advanced radec methods tight to astropy SkyCoords"""
        ra,dec = self.radec
        return coordinates.SkyCoord(ra=ra*units.degree,dec=dec*units.degree)
    
    # ----------------      
    # -- CATALOGUE
    @property
    def catalogue(self):
        return self._side_properties["catalogue"]
    
    def has_catalogue(self):
        return False if self.catalogue is None\
          else True

    @property
    def pairdict(self):
        if self._derived_properties["pairdict"] is None:
            self._derived_properties["pairdict"] = self._get_pairing_()
            
        return self._derived_properties["pairdict"]
        
    @property
    def catmask(self):
        """this array tells you if sextractor objects have a catalogue match"""
        if not self.has_catmatch() and self.has_data():
            return np.ones(self.nobjects,dtype=bool)
        return self._derived_properties["catmatch"]["idx"]

    @property
    def starmask(self):
        """ This is which are stars amoung *ALL* the detected sources, even those without matcatching"""
        return self._derived_properties["starmask"]
    

    @property
    def isolatedmask(self):
        return self.catalogue.isolatedmask[self.catmcatch["idx_catalogue"]]
    
    @property
    def catstarmask(self):
        """return a bool mask of the stars. If it can not do it, not data are masked (full True)"""
        if not self.has_catmatch() and self.has_data():
            return np.ones(self.nobjects,dtype=bool)
        
        if not self.catalogue.has_starmask():
            return np.ones(self.nobjects,dtype=bool)
        
        return self.catalogue.starmask[ self.catmatch["idx_catalogue"]]

    @property
    def catisolatedmask(self):
        """return a bool mask of the stars. If it can not do it, not data are masked (full True)"""
        if not self.has_catmatch() and self.has_data():
            return np.ones(self.nobjects,dtype=bool)
        
        if not self.catalogue.has_starmask():
            return np.ones(self.nobjects,dtype=bool)
        
        return self.catalogue.isolatedmask[ self.catmatch["idx_catalogue"]]
    
    
    # ----------------------
    # - Catalogue Matching
    @property
    def catmatch(self):
        """This is the match dictionnary"""
        if self._derived_properties["catmatch"] is None:
            self._derived_properties["catmatch"] = {}
            
        return self._derived_properties["catmatch"]

    
    def has_catmatch(self):
        return False if self.catmatch is None or len(self.catmatch.keys())==0 \
          else True
    
