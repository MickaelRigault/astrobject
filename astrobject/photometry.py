#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module defines the photometric objects"""

import warnings
import numpy  as np

from astropy.io import fits as pf

from scipy.stats import sigmaclip

from astropy     import units,coordinates
from astropy.table import Table

from . import astrometry
from .baseobject   import BaseObject 
from ..utils.tools import kwargs_update,flux_to_mag

__all__ = ["image","photopoint"]

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

def photopoint(lbda=None,flux=None,var=None,
               zp=None,bandname=None,
               mjd=None,source=None,
               instrument_name=None,**kwargs):
    """
    Get the basic object containing  photometric point information
    
    Parameters
    ----------

    lbda: [float]              The central wavelength associated to the photometric
                               points.

    flux: [float]              The flux (*not magnitude*) of the photometric point.

    var: [float]               The variance associated to the point's flux.

    zp: [float]                Zeropoint of the instrument's image

    bandname: [string]         Name of the Bandpass though which the observation is made
        
    - option -
    
    empty: [bool]              Set True to return an empty object.

    - other options ; not exhautive ; goes to 'create' -

    mjd: [float]               Modified Julian Data of the observation
    
    source: [string]           Staten the origin of the point (e.g. image, ...)

    instrument_name:[string]   Give a name of the intrument that enable to take the
                               photometric point.

    Return
    ------
    PhotoPoint
    """
    # -------------
    # - Parser
    if "variance" in kwargs.keys() and var is None:
        var = kwargs.pop("variance")
    if "wavelength" in kwargs.keys() and lbda is None:
        lbda = kwargs.pop("wavelength")
    if "zpsystem" in kwargs.keys():
        kwargs["zpsys"] = kwargs.pop("zpsystem")
        
    return PhotoPoint(lbda,flux,var,source=source,
                      instrument_name=instrument_name,
                      mjd=mjd,zp=zp,bandname=bandname,
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
    _side_properties_keys    = ["wcs","datamask",
                                "target","catalogue",
                                "exptime"] # maybe Exptime should just be on flight
        
    _derived_properties_keys = ["fits","data","sepobjects","backgroundmask",
                                "apertures_photos","fwhm"]
    
    # Where in the fitsfile the data are
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,filename=None,
                 astrotarget=None,data_index=0,
                 dataslice0=[0,-1],dataslice1=[0,-1],
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
                      dataslice0=dataslice0,dataslice1=dataslice1,
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
                header_exptime = "EXPTIME",
                dataslice0="undefined",
                dataslice1="undefined"
                )

    # =========================== #
    # = Main Methods            = #
    # =========================== #
    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def load(self,filename,index=None,
             force_it=False,
             dataslice0=[0,-1],dataslice1=[0,-1]):
        """
        This enables to load a fitsfile image and will create
        the basic data and wcs solution if possible.
        *var* (error) and *background* has to be defined
        separately has this strongly depend on the instrument

        Parameters
        ----------
        filename: [string.fits]    The file containing the fits data

        - options -

        dataslice0/1 [2D-array]    load only the data within the given boundaries.
                                   The 0-offset will be accessible in self._dataoffset
                                   and will be passed to the wcs solution.
        
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
        #  fits file and wcs         #
        # -------------------------- #
        try:
            fits = pf.open(filename,memmap=True)
        except:
            raise TypeError("'filename' cannot be open by pyfits.")

        try:
            wcs_ = astrometry.wcs(filename,extension=index)
        except:
            wcs_ = None
            
        # ---------- #
        # - Data   - #
        # ---------- #
        data = fits[index].data[dataslice0[0]:dataslice0[1],
                                dataslice1[0]:dataslice1[1]]
        self._build_properties["dataslice0"] = dataslice0
        self._build_properties["dataslice1"] = dataslice1            
        # -------------------------- #
        #  Everythin looks good !    #
        # -------------------------- #
        self.create(data,None,wcs_,
                    background= None,
                    header=fits[index].header,
                    filename=filename,fits=fits,
                    force_it=True)
        
            
    def create(self,rawdata,variance,wcs,mask=None,
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
        # --------------
        # - Set instance
        self._properties["filename"]      = filename
        self._derived_properties["fits"]  = fits
        self._properties["header"]        = pf.Header() if header is None else header
        self._side_properties["exptime"]  = exptime
        # --------------
        # - Read data
        self._set_data_(rawdata,mask,
                        variance=variance,
                        background=background)
        # - WCS solution
        self.set_wcs(wcs)
        
    def reload_data(self, dataslice0, dataslice1,
                    variance=None,background=None,mask=None):
        """ Change the slicing of the data.
        This will update the background, variance and the wcs solution's offset.
        If you had an sepobjects loaded, this will update it."""
        
        rawdata = self.fits[self._build_properties["data_index"]].data[
            dataslice0[0]:dataslice0[1],
            dataslice1[0]:dataslice1[1]]
        # -- no more sepobject attached
        reload_sep = self.has_sepobjects()

        self._derived_properties["sepobjects"] = None
        self._set_data_(rawdata,mask=mask,
                        variance=variance,
                        background=background)
        
        self._build_properties["dataslice0"] = dataslice0
        self._build_properties["dataslice1"] = dataslice1
        if self.has_wcs():
            self.wcs.set_offset(*self._dataslicing)
            if self.has_catalogue():
                self.catalogue.set_fovmask(wcs=self.wcs)
        
        if reload_sep: self.sep_extract()
        
    # ------------------- #
    # - Set Methods     - #
    # ------------------- #
    def set_target(self,newtarget,test_inclusion=True):
        """
        Change (or create) an object associated to the given image.

        This function will test if the object (the target) is withing
        the image boundaries (expect if 'test_inclusion' is set to False).
        Set 'newtarget' to None to remove the association between this
        object and a target.

        Return
        ------
        Void
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
        
    def set_wcs(self,wcs,force_it=False):
        """
        """
        if self.has_wcs() and not force_it:
            raise AttributeError("A wcs solution is already loaded."\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(wcs) or wcs.__nature__ != "WCS":
            raise TypeError("The given wcs solution is not a astrobject WCS instance")
        self._side_properties["wcs"] = astrometry.get_wcs(wcs)
        self.wcs.set_offset(*self._dataslicing)
        
    def set_catalogue(self,catalogue,force_it=False):
        """
        A Catalogue object will be loaded to this instance, you could then access it
        as 'self.catalogue' ('self.has_catalogue()' will be True).

        The wcs solution of this
        instance will be passed to the 'calague' one, as well as the sepobject
        if you have one. This way, the detected object will directly be match to the
        'catalogue'

        Catalogue: must be an astrobject Catalogue
        
        Return
        ------
        Void
        """
        if self.has_catalogue() and force_it is False:
            raise AttributeError("'catalogue' already defined"+\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(catalogue) or catalogue.__nature__ != "Catalogue":
            raise TypeError("the input 'catalogue' must be an astrobject catalogue")
        # -------------------------
        # - Add the world_2_pixel
        
        if not self.has_wcs():
            print "WARNING: no wcs solution."+"\n"+\
              " catalogue recorded, but most likely is useless"
            
        else:
            # -- Lets only consider the objects in the FoV
            catalogue.set_wcs(self.wcs,force_it=True)
            if catalogue.nobjects_in_fov < 1:
                warnings.warn("WARNING No object in the field of view,"+"\n"+\
                              "  -> catalogue not loaded")
                return
            
            # -- Lets save the pixel values
            if self.has_sepobjects():
                self.sepobjects.set_catalogue(catalogue,force_it=True)
                self.sepobjects.match_catalogue()
            
                
        # --------
        # - set it
        self._side_properties["catalogue"] = catalogue
        
    def set_background(self,background,
                       force_it=False,check=True):
        """
        This is a method that might strongly depend on the instrument.
        As a default (background = None) this uses Sextractor background
        estimation from 'get_sep_background'.

        Give background or overwrite this method for your specific instrument.

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

        if background is None:            
            background  = self._get_default_background_()
            self._uses_default_background = True
        else:
            self._uses_default_background = False
            
        # Shape test
        if self.rawdata is not None and np.shape(background) != self.shape:
            raise ValueError("The given background must have rawdata's shape")
        # -- Looks good
        self._properties['background'] = np.asarray(background)
        self._update_data_(update_background=False)
        
    def set_fwhm(self,value,force_it=True):
        """
        value is the value of the fwhm. If no units is provided (astropy units)
        arcsec will be assumed.
        """
        if self.has_fwhm() and not force_it:
            raise AttributeError("'fwhm' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        if value<0:
            raise ValueError("the 'fwhm' must be positive")
        if type(value) is not units.quantity.Quantity:
            value = value*units.arcsec
            
        self._derived_properties['fwhm'] = value
        
    # --------------------- #
    # - download Methods  - #
    # --------------------- #
    def download_catalogue(self,source="sdss",
                           set_it=True,force_it=False,
                           **kwargs):
        """
        Downloads a catalogue of the given 'source'. This methods requires an
        internet connection. This downloaded catalogue is then loaded in this
        instance using the 'set_catalogue' method. If 'set_it' is set to False,
        the catalogue is returned instead of being loaded.
        
        kwargs goes to instrument.catalogue
        example: column_filters={"gmag":"13..22"}

        Return:
        ------
        Void (or the Catalogue if set_it is False)
        """
        from .instruments import instrument
        radec = "%s %s"%(self.wcs._central_coords_nooffset[0],
                         self.wcs._central_coords_nooffset[1])
        radius = self.wcs.diag_size/1.8 # not 2 to have some room around
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
    def get_stars_aperture(self, radius,runits="pixels",aptype="circle",
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

        if isolated_only and not self.catalogue._is_around_defined():
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
        return idx, cat_idx, self.get_aperture(x,y,radius=radius,runits=runits,
                                               aptype=aptype,**kwargs)
    
    def get_target_aperture(self,radius,runits="pixels",aptype="circle",subpix=5,
                            **kwargs):
        """If a target is loaded, this will get the target coords and run
        'get_aperture'.
        """
        if self.target is None:
            raise AttributeError("No 'target' loaded")
        xpix,ypix = self.coords_to_pixel(self.target.ra,self.target.dec)
        return self.get_aperture(xpix,ypix,
                                 radius=radius,runits=runits,
                                 aptype=aptype,
                                 subpix=subpix,**kwargs)

    def get_aperture(self,x,y,radius=None,
                     runits="pixels",wcs_coords=False,
                     aptype="circle",subpix=5,
                     ellipse_args={"a":None,"b":None,"theta":None},
                     annular_args={"rin":None,"rout":None},
                     **kwargs):
        """
        This method uses K. Barary's Sextractor python module SEP
        to measure the aperture photometry. See details here:
        sep.readthedocs.org/en/v0.4.x/apertures.html
        *Remark* you can use the radec option to give x,y as radec coordinates.
        This will only works if you have a wcs solution loaded.
        
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
                                   -circle  => set radius
                                   -ellipse => set all ellipse_args entries
                                   -circan  => set all annulus_args entries
                                   -ellipan => set all ellipse_args entries
                                              set all annulus_args entries
                                    (no other type allowed)
                                               
        
        radius: [float]            Size of the circle radius. In pixels.
                                   (This is used only if aptype is circle)
                                   
        runits: [str/astropy.units] The unit of the radius

        wcs_coords: [bool]         Set True if x,y are ra,dec coordinates
        
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

        # ---------------
        # Input Parsing
        # - radius
        r_pixels = radius*self.units_to_pixels(runits)
        # - position
        if wcs_coords and not self.has_wcs():
            raise AttributeError("you cannot provide ra,dec coordinate without a wcs solution. cannot convert them into pixel coords")
        if wcs_coords:
            x,y = self.coords_to_pixel(x,y).T 
            
            
        # -------------
        # - SEP Input 
        gain = None if "_dataunits_to_electron" not in dir(self) else \
          self._dataunits_to_electron
        gain = kwargs.pop("gain",gain)
        # -----------------
        # - Variance Trick
        var = kwargs.pop("var",self.var)

        # ------------------------
        # - Do the Aperture Photo
        # ----------------------
        # - GAIN MEANS Conversion factor between data array units and poisson counts,
        # gain = ADU per Electron, data = ADU/s
        # - Circle
        if aptype == "circle":
            return sep.sum_circle(self.data,x,y,r_pixels,subpix=subpix,
                            var=var,gain=gain,**kwargs)

        # - Annulus
        if aptype == "circann":
            
            if np.asarray([k is None for k in annular_args.values()]).any():
                raise ValueError("You must set the annular arguments 'annular_arg'")
            
            rin,rout = [annular_args[k] for k in ["rin","rout"]]
            sepout= sep.sum_circann(self.data,x,y,rin,rout,
                                  r_pixels,subpix=subpix,
                                  var=var,gain=gain,**kwargs)
        # - Ellipse
        elif aptype == "ellipse":
            if np.asarray([k is None for k in ellipse_args.values()]).any():
                raise ValueError("You must set the ellipse arguments 'ellipse_arg'")
            
            a,b,theta = [annular_args[k] for k in ["a","b","theta"]]
            sepout= sep.sum_ellipse(self.data,x,y,a,b,theta,
                                    r_pixels,subpix=subpix,
                                    var=var,gain=gain,**kwargs)
        # - Elliptical Annulus
        elif aptype == "ellipan":
            if np.asarray([k is None for k in ellipse_args.values()]).any():
                raise ValueError("You must set the ellipse arguments 'ellipse_arg'")
            if np.asarray([k is None for k in annular_args.values()]).any():
                raise ValueError("You must set the annular arguments 'annular_arg'")
            
            rin,rout = [annular_args[k] for k in ["rin","rout"]]
            a,b,theta = [ellipse_args[k] for k in ["a","b","theta"]]
            sepout= sep.sum_ellipan(self.data,x,y,a,b,theta,rin,rout,
                                    r_pixels,subpix=subpix,
                                    var=var,gain=gain,**kwargs)
        # ===================================
        # = Return of the aperture extraction
        return sepout
    
    # ------------------- #
    # - WCS Tools       - #
    # ------------------- #
    def get_contours(self,pixel=True):
        """Contours (Shapely) of the image. This refere's to the wcs solution"""
        if not self.has_wcs():
            raise NotImplementedError("Only WCS contours implemented and this has no wcs solution")
        
        if pixel:
            return self.wcs.contours_pxl
        return self.wcs.contours

    # - Conversion tools
    def units_to_pixels(self,units_):
        """units should be a parsable string or a astropy units. Uses the wcs method units_to_pixels()
        In addition, you have access to the 'psf [==fwhm]' unit"""
        if type(units_) == str and units_.lower() in ["fwhm","psf"]:
            units_ = self.fwhm
            
        return self.wcs.units_to_pixels(units_,target=self.target)
    
            
    def pixel_to_coords(self,pixel_x,pixel_y):
        """get the coordinate (ra,dec; degree) associated to the given pixel (x,y)"""
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded.")
        pixelsoffset = self._dataoffset
        return self.wcs.pix2world(pixel_x,pixel_y)

    def coords_to_pixel(self,ra,dec):
        """Return the pixel (x,y) associated to the given ra,dec (degree) coordinate"""
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded.")
        # Remark the np.asarray are only required by the astLib wcs solution
        return np.asarray(self.wcs.world2pix(ra,dec))

    # ------------------- #
    # - SEP Tools       - #
    # ------------------- #
    def get_sep_background(self,update_background=True,cleaning_sep=True,**kwargs):
        """
        This module is based on K. Barbary's python module of Sextractor: sep.
        
        """
        self._measure_sep_background_(**kwargs)
        if self.has_sepobjects():
            # -----------------
            # - First loop get the first exposure
            if "_rmsep" in dir(self) and self._rmsep:
                self._derived_properties["sepobjects"] = None
            # -- No need to conserve that
            del self._rmsep
            
            return self._sepbackground.back()
        
        self.set_background(self._sepbackground.back(),force_it=True)
        self.sep_extract(match_catalogue= not cleaning_sep)
        self._rmsep=cleaning_sep
        return self.get_sep_background(update_background=True)

    
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

        set_catalogue: [bool]      If the current instance has a catalogue, it will be
                                   transfered to the SexOutput object created.
                                   Set False to avoid that.
                                   
                                   
        returnobjects: [bool]      Change the output of this function. if True the
                                   extracted objects are recorded and returned
                                   (self.sepobjects) if not they are just recorded.

        
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
             ax=None,show=True,zoomon=None,zoom=200,zunits="pixels",
             show_sepobjects=False,propsep={},
             show_catalogue=False,proptarget={},
             **kwargs):
        """
        Display the 2D-image. The displayed information can be "data", "background",
        "rawdata", "variance" or anything known by self.`toshow` that is a 2D image.

        Parameters
        -----------

        - options -

        // Data/display options
        
        toshow: [string/array]     key of the element you want to show. Could be
                                   data, rawdata, background, variance.

        logscale: [bool]           The values of the 2D array is scaled in log10
                                   before being displayed.
                                   
        // Additional Information

        show_sepobjects: [bool]   Add in the plot the contours about the detected
                                  sources. extract_sep() must have been ran. If a
                                  catalogue is loaded, only the catalogue matched
                                  sources will be displayed.
                                  
        propsep: [dict]           If the sepobjects are display, this set the
                                  properties of the ellipse (mpl.Patches)

        show_catalogue: [bool]   Add in the location of the stars (full markers) and
                                 galaxies (open marker) of the catalogue objects.
                                 If no catalogue set to the instance this is set
                                 to False.

        proptarget: [dict]       Change the properties of the marker showing the target
                                 location. If no target set to the instance, this won't
                                 affect anything.
                                 
        // Zoom
        
        zoomon: [string/2D/None]   Zoom on the image centered in the given 'zoomin'
                                   [x,y] pixel coordinates. You can also set 'target'
                                   and if a target is loaded, this will parse 'target'
                                   to the target's [x,y] coordinates.
                                   Set None not to zoom.

        zoom: [float]              If you zoomin (zoonin!=None) this set the  box
                                   of the zoom (width and height of 2*zoom).
                                   
        zunits: [string]          units of the zoom value (pixels, fhhm, kpc, arcsec..)
        
        // Data/display options
        
        savefile: [None/string]    Give the name of the file where the plot should be
                                   saved (no extention, .pdf and .png will be made)
                                   If None, nothing will be saved, the plot will be
                                   shown instead (but see *show*)

        show: [bool]               If the image is not recorded (savefile=None) the
                                   plot will be shown except if this is set to False
                                   
        ax: [None/mpl.Axes]        Where the imshow will be displayed. If None, this
                                   will create a figure and a unique axis to set 'ax'

        -- Additional options --
        
        **kwargs goes to matplotib's imshow (vmin, vmax ...).
          *Important* use string values for 'vmin' and 'vmax' to considered them are
          percentile of the data values.

        Return
        ------
        dictionnary {'figure','ax','plot','prop'}

        """
        # ----------------
        # - Input test
        if type(toshow) is np.ndarray:
            if np.shape(toshow) != (self.height, self.width):
                raise ValueError("the given ndarray 'toshow' must have the size of the image (%d,%d)"%(self.height, self.width))
            valuetoshow = toshow
        elif type(toshow) is not str:
            raise TypeError("'toshow' must be a string (self.'toshow') or a np.ndarray)"%toshow)
        elif toshow not in dir(self):
            raise ValueError("'%s' is not a known image parameter"%toshow)
        else:
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
            ax.set_xlabel("x",fontsize = "x-large")
            ax.set_ylabel("y",fontsize = "x-large")
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
            
        prop = kwargs_update(default_prop,**kwargs)

        # ----------- #
        # - Do It     #
        
        if "vmin" in prop.keys() and type(prop["vmin"]) is str:
            prop["vmin"] = np.percentile(x[x==x],float(prop["vmin"]))
        if "vmax" in prop.keys() and type(prop["vmax"]) is str:
            prop["vmax"] = np.percentile(x[x==x],float(prop["vmax"]))
        im = ax.imshow(x,**prop)
        
        # - add target
        pl_tgt = None if self.has_target() is False \
          else self.display_target(ax,wcs_coords=False,
                                   **proptarget)

        if show_sepobjects and self.has_sepobjects():
            self.sepobjects.display(ax,
                                    **propsep)
        if show_catalogue and self.has_catalogue():
            self.display_catalogue(ax,wcs_coords=False)
        # ----------- #
        # - Zoom      #
        if zoomon is not None:
            # -- Zoom on target
            if type(zoomon) is str and zoomon=="target" and self.has_target():
                    coords_zoom = self.coords_to_pixel(*self.target.radec) 
            elif np.shape(zoomon) == (2,):
                coords_zoom = zoomon
            else:
                warnings.warn("(Image.show) can not parse to zoom on input")
                coords_zoom = None
        else:
            coords_zoom = None

        if coords_zoom is not None:
            width = zoom * self.units_to_pixels(zunits)
            ax.set_xlim(coords_zoom[0]-width,coords_zoom[0]+width)
            ax.set_ylim(coords_zoom[1]-width,coords_zoom[1]+width)
            
        # ----------- #
        # - Recordit
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["imshow"] = im
        self._plot["prop"]   = prop
        self._plot["target_plot"] = pl_tgt
        
        fig.figout(savefile=savefile,show=show)
        
        return self._plot

    
    def show_hist(self,toshow="data",savefile=None,logscale=True,
                ax=None,show=True,proptarget={},
                **kwargs):
        """
        This methods enable to show histogram distribution of the concatenate data.
        
        """
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
        elif "hist" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
        # ----------- #
        # -  What
        _x = np.log10(valuetoshow) if logscale else valuetoshow
        x = np.concatenate(_x)
        x = x[x==x]
        # ----------- #
        # - How
        default_prop = {
            "histtype":"step",
            "fill":"True",
            "lw":2,"ec":mpl.cm.Blues(0.8,0.8),
            "fc":mpl.cm.Blues(0.6,0.3),
            "bins":100,
            "range":np.percentile(x,[1,99])
            }
            
        prop = kwargs_update(default_prop,**kwargs)

        # ----------- #
        # - Do It     #
        pl = ax.hist(x,**prop)

         # ----------- #
        # - Recordit
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["hist"] = pl
        self._plot["prop"]   = prop
        
        fig.figout(savefile=savefile,show=show)
        
        return self._plot


    def show_background(self,savefile=None,show=True,
                        **kwargs):
        """
        """
        # -- Setting -- #
        from ..utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        # ----------- #
        # - Where     #
        fig  = mpl.figure(figsize=[16,8])
        axB  = fig.add_axes([0.1, 0.1,0.4,0.8])
        axM  = fig.add_axes([0.52,0.1,0.4,0.8])
        
        # ----------- #
        # - What      #
        backgroudsource = self.rawdata.copy()
        backgroudsource[self.backgroundmask] = np.NaN
        # ----------- #
        # - How
        # -- labels
        prop = kwargs_update({"vmin":"3","vmax":"97"},**kwargs)
        axB.set_xlabel("x",fontsize = "x-large")
        axM.set_xlabel("x",fontsize = "x-large")
        axB.set_ylabel("y",fontsize = "x-large")
        # -- titles
        axB.set_title(r"$\mathrm{background}$",fontsize = "xx-large")
        axM.set_title(r"$\mathrm{masked\ rawdata\ (log)\ used\ to\ create\ the\ background}$",fontsize = "xx-large")
        
        if "logscale" in kwargs.keys():
            print "No logscale option available for show_background."
            _ = kwargs.pop("logscale")
        # ----------- #
        # - Do It     #
        _plotb  = self.show("background",ax=axB,logscale=False,
                        show=False,savefile=None,**prop)
        _plotbm = self.show(backgroudsource,ax=axM,logscale=True,
                        show=False,savefile=None,**prop)
        
        # ----------- #
        # - Done     #
        del backgroudsource
        # ----------- #
        # - Recordit
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["axes"]   = [axB,axM]
        self._plot["imshows"] = [_plotb["imshow"],_plotbm["imshow"]]
        
        fig.figout(savefile=savefile,show=show)
        
        return self._plot
    
    # ---------------------- #
    # - Plot-Displays      - #
    # ---------------------- #
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

    def display_sepobjects(self,ax=None,draw=True,**kwargs):
        """If sep_extract has been ran, you have an sepobjects entry.
        This entry will be red and parsed here.
        """
        # --------------
        # - ax parsing
        if ax is None and ("_plot" not in dir(self) or "ax" not in self._plot.keys()):
            raise ValueError('no ax defined')
        ax = self._plot['ax'] if ax is None else ax
        
        self.sepobjects.display(ax,draw=draw,
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

        
    @property
    def _dataslicing(self):
        """ return offset0, offset1, height, width """
        if self._build_properties["dataslice0"] == "undefined":
            warnings.warn("No _build_properties['dataslice0'] defined. 0 assumed")
            offset0 = [0,-1]
        else:
            offset0=self._build_properties["dataslice0"]
        if self._build_properties["dataslice1"] == "undefined":
            warnings.warn("No _build_properties['dataslice1'] defined. 0 assumed")
            offset1 = [0,-1]
        else:
            offset1=self._build_properties["dataslice1"]

        if offset0[1] <0: offset0[1] = self.height+offset0[1]+2
        if offset1[1] <0: offset1[1] = self.width +offset1[1]+2
        return offset0[0],offset1[0],offset0[1]-offset0[0],offset1[1]-offset1[0]
    
    @property
    def _dataoffset(self):
        return self._dataslicing[:2]
            
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

    @property
    def backgroundmask(self):
        return self._derived_properties["backgroundmask"]
    
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
    def datamask(self):
        return self._side_properties["datamask"]
    
    def has_datamask(self):
        """ Test you if defined a datamask key. This mask is True
        for data you wish to remove"""
        return self.datamask is not None
        
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
        return self.wcs.pix_indeg
          
    @property
    def pixel_size_arcsec(self):
        """Pixel size in arcsec. Based on wcs solution"""
        if type(self.pixel_size_deg) is not units.quantity.Quantity:
            return [ps.to("arcsec") for ps in self.pixel_size_deg]
        return self.pixel_size_deg.to("arcsec")

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
          
    @property
    def sepmask(self,r=10):
        if not self.has_sepobjects():
            raise AttributeError("No sepobjects loaded. Run sep_extract")
        return self.sepobjects.get_ellipse_mask(self.width,self.height,r=r)
            
    # ----------------      
    # -- CATALOGUE
    @property
    def catalogue(self):
        return self._side_properties["catalogue"]
    
    def has_catalogue(self):
        return False if self.catalogue is None\
          else True

    
    # FWHM
    @property
    def fwhm(self):
        """ full width half maximum in arcsec.
        If you did not set it manually,
        this will use sepobjects' get_fwhm_pxl() method to set it.
        """
        if not self.has_fwhm():
            if self.has_sepobjects():
                fwhm_pxl = self.sepobjects.get_fwhm_pxl(isolated_only=True,
                                                        stars_only=True)
                self.set_fwhm(fwhm_pxl/self.units_to_pixels("arcsec").value*\
                              units.arcsec)
            else:
                raise AttributeError("'fwhm' is not defined and no sepobjects loaded.")
        return self._derived_properties["fwhm"]

    def has_fwhm(self):               
        return not self._derived_properties["fwhm"] is None
    
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _set_data_(self,rawdata,mask=None,
                   variance=None,
                   background=None):
        """ Change the instance fondamental: the rawdata """
        
        self._properties["rawdata"]       = self._read_rawdata_(rawdata)
        self._side_properties["datamask"] = mask
        if self.has_datamask():
            self._properties["rawdata"][self.datamask] = np.NaN

        # VARIANCE
        self._properties["var"]          = variance
        if self.has_datamask() and self.has_var():
            self._properties["var"][self.datamask] = np.NaN
        # BACKGROUND   
        self.set_background(background, force_it=True)
        # --> update
        self._update_(update_background=False)
    
    def _read_rawdata_(self,rawdata):
        return np.asarray(rawdata,dtype="float")
    
    def _get_default_variance_(self):
        """
        """
        if "_sepbackground" in dir(self):
            return self._sepbackground.rms()**2
        return None
        
    def _get_default_background_(self,*args,**kwargs):
        
        return self.get_sep_background(*args,**kwargs)

    
    def _measure_sep_background_(self,**kwargs):
        """
        """
        from sep import Background
        
        if self.rawdata is None:
            raise ValueError("no 'rawdata' loaded. Cannot get a background")
        
        # --------------
        # -- Mask issue
        if self.has_datamask():
            maskdata = self.datamask
        else:
            maskdata = None
            
        if self.has_sepobjects():
            masksep = self.sepmask
            mask = masksep+maskdata if maskdata is not None else\
              masksep
        else:
            mask = maskdata
        # ---------------
        # - masking sign tested
        self._derived_properties["backgroundmask"] = mask
        self._sepbackground_prop = kwargs_update({"mask":self.backgroundmask ,
                                                  "bw":100,"bh":100},
                                                  **kwargs)
        
        self._sepbackground = Background(self.rawdata,
                                         **self._sepbackground_prop)
        
    def _update_(self,update_background=True):
        """The module derives the 'derived_properties' based on the
        fundamental once
        """
        # -- Make sure the fundamental update (if any) are made
        super(Image,self)._update_()
        # - Data
        self._update_data_(update_background=update_background)

    def _update_data_(self,update_background=True):
        """
        """
        if update_background:
            if "_sepbackground" in dir(self):
                self._measure_sep_background_(**self._sepbackground_prop)
                if self._uses_default_background:
                    self.set_background(self._get_default_background_(),force_it=True)
                
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

    _properties_keys = ["lbda","flux","var","mjd","bandname","zp"]
    _side_properties_keys = ["source","intrument_name","target","zpsys",
                             "meta"]
    _derived_properties_keys = []
    
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,lbda=None,flux=None,var=None,mjd=None,
                 bandname=None,zp=None,zpsys="ab",
                 empty=False,**kwargs):
        """
        Initialize the PhotoPoint object

        Parameters
        ----------

        lbda: [float]              The central wavelength associated to the photometric
                                   points.

        flux: [float]              The flux (*not magnitude*) of the photometric point.

        var: [float]               The variance associated to the point's flux.

        mjd: [float]               Modify Julian Date of the Observation

        bandpass: [string]         Name of the Bandpass though which the observation is made

        zp: [float]                Zeropoint of the instrument's image
        
        - option -

        empty: [bool]              Set True to return an empty object.
        
        - other options ; not exhautive ; goes to 'create' -

        source: [string]           Staten the origin of the point (e.g. image, ...)

        instrument_name:[string]   Give a name of the intrument that enable to take the
                                   photometric point.

        *META*  any additional key will be stored as a dictionary accessible with the meta entry
        and more generally with the get() method.
        
        Return
        ------
        Void
        """
        self.__build__()
        if empty:
            return
        prop = kwargs_update(dict(mjd=mjd,zp=zp,bandname=bandname),
                             **kwargs)
        self.create(lbda,flux,var,**prop)    
        
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def create(self,lbda,flux,var,
               source=None,instrument_name=None,
               mjd=None,zp=None,bandname=None,zpsys="ab",
               force_it=False,**meta):
        """
        This method creates the object by setting the fundamental parameters.

        
        Parameters:
        -----------

        lbda,flux,var: [floats]    'lbda' is the effective wavelength (in Angstrom) of the
                                   bandpass through which the 'flux' and the 'var' (variance) have
                                   been measured.

        - options -

        mjd: [float]               Modified Julian date of the observation.

        zp: [float]                the zeropoint (in ABmag) of the instrument used to
                                   derive the photopoint. This is convinient when several photopoints
                                   are combined.
                                   
        source: [string]           is method used the measure this photometric point (e.g. image, spectrum...)
        

        instrument_name: [string]  set here the name of the instrument used to derive this photopoint
                                   (e.g. sdss)

        bandname: [string]         name of the band used to derive the flux. If you set this you can
                                   access the 'sncosmo bandpass'.

        Returns:
        --------
        Void
        
        """
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
        self._side_properties["meta"] = meta
        # -- Interactive ones
        self._side_properties["zpsys"] = zpsys
        self.mjd = mjd
        self.zp = zp
        self.bandname = bandname
        
        self._update_()

    def display(self,ax,toshow="flux",function_of_time=False,**kwargs):
        """This method enable to display the current point
        in the given matplotlib axes"""
        # - Test if there is data
        if not self.has_data():
            raise AttributeError("no data to display")
        # -----------
        # - Input
        y = self.get(toshow)
        if toshow == "flux":
            dy= np.sqrt(self.var) if self.var is not None else None
        elif toshow == "mag":
            dy= np.sqrt(self.magvar) if self.magvar is not None else None
            

        # -----------
        # - Fancy
        default_prop = dict(marker="o",ecolor="0.7",
                        zorder=3)
        prop = kwargs_update(default_prop,**kwargs)
        # -----------
        # - Input
        x_ = self.lbda if not function_of_time else self.mjd
        pl = ax.errorbar(x_,y,yerr=dy,**prop)
        self._plot = {}
        self._plot["ax"] = ax
        self._plot["plot"] = pl
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


    def get(self, key, safeexit=False):
        """ Generic method to access information of the instance.
        Taken either from the instance itself self.`key` or from the meta parameters.

        *Remark* 'key' could be an list of keys.

        Returns:
        --------
        Value (or list)
        """
        ## Tested, try except faster than if key in dir(self) and enable things like key="meta.keys()"

        if "__iter__" in dir(key):
            return [self.get(key_,safeexit=safeexit) for key_ in key]
        
        try: # if key in dir(self):
            return eval("self.%s"%key)
        except:
            if key in self.meta.keys():
                return self.meta[key]
            elif not safeexit:
                raise ValueError("No instance or meta key %s "%key)
            warnings.warn("No instance or meta key %s. NaN returned"%key)
            return np.NaN
        
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

    @property
    def mjd(self):
        return self._properties["mjd"]
    
    def has_mjd(self):
        return not self.mjd is None
    
    @mjd.setter
    def mjd(self,value):
        if value < 48987:
            print "Information mjd prior to 1993..."
        elif value > 58849.0:
            print "Information mjd posterior to 2020..."
        
        self._properties["mjd"] = value

    @property
    def zp(self):
        return self._properties["zp"]
    
    @zp.setter
    def zp(self,value):
        if value <0:
            raise ValueError("a zp cannot be negative")
        if value >35:
            warnings.warn("(PhotoPoint instance) the given zp is above 35...")
        
        self._properties["zp"] = value
        
    @property
    def bandname(self):
        return self._properties["bandname"]
    
    @bandname.setter
    def bandname(self,value):
        if type(value) != str and type(value) != np.string_:
            raise TypeError("The bandname must be a string", type(value))
        self._properties["bandname"] = value

    @property
    def bandpass(self):
        """
        """
        if self.bandname is None:
            raise AttributeError("No bandname given")
        
        from sncosmo import get_bandpass
        return get_bandpass(self.bandname)
        
    # ------------
    # - Side
    @property
    def zpsys(self):
        return self._side_properties["zpsys"]
    
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

    # -- Meta
    @property
    def meta(self):
        if self._side_properties["meta"] is None:
            self._side_properties["meta"]= {}
        return self._side_properties["meta"]
    
    # ------------
    # - Derived
    # magnitudes
    @property
    def mag(self):
        return flux_to_mag(self.flux,np.sqrt(self.var),self.lbda)[0]

    @property
    def magvar(self):
        return flux_to_mag(self.flux,np.sqrt(self.var),self.lbda)[1] ** 2    

    @property
    def magabs(self):
        if not self.has_target():
            raise AttributeError("No target defined, I can't get the distance")
        return self.mag - 5*(np.log10(self.target.distmpc*1.e6) - 1)



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
        print "Decrepated... To be moved to PhotoPointCollection'"
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
    # - Get Method     - #
    # ------------------ #
    def get_data(self):
        """
        This will be an SN cosmo like data. This is based on astropy table
        """
        from astropy.table.table import Table
        return Table(data=[self.times, self.bandnames,
                           self.fluxes,np.sqrt(self.fluxesvar),
                           self.zps,[p.zpsys for p in self.photopoints]
                           ],
                    names=["time","band","flux","fluxerr","zp","zpsys"])
    
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
    def npoints(self):
        return len(self.photopoints)
    @property
    def photopoints(self):
        return self._properties['photopoints']
    
    @property
    def times(self):
        return self._properties['times']
    
    @property
    def lbda(self):
        return self._derived_properties['lbda']

    # -------------
    # - On the flight
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
    def bandnames(self):
        # - a unique lbda is tested while loading so this should be always the same
        return [p.bandname for p in self.photopoints]

    @property
    def zps(self):
        # - a unique lbda is tested while loading so this should be always the same
        return [p.zp for p in self.photopoints]

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
    # = Main Methods       = #
    # ====================== #
    def create(self,sexoutput,force_it=False):
        """
        """
        if self.has_data() and force_it is False:
            raise AttributeError("'data' is already defined."+\
                    " Set force_it to True if you really known what you are doing")
                    
        sexdata = self._read_sexoutput_input_(sexoutput)
        if sexdata is None:
            warning.warn("(SexObjects.create) empty input data. *Empty SexObjects loaded*")
            return
        
        # ****************** #
        # * Creation       * #
        # ****************** #
        self._properties["data"] = sexdata

        
    def match_catalogue(self,catalogue=None,force_it=False,arcsec_size=2):
        """
        This methods enable to attached a given sexobject entry
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

    def idx_to_mask(self,idx):
        """
        Change the list of index to a True/False mask: index in the list
        are the True values

        Return:
        -------
        mask (boolean array)
        """
        mask = np.zeros(self.nobjects,dtype=bool)
        for i in idx:
            mask[i] = True
        return np.asarray(mask,dtype=bool)
        
    # ------------------ #
    # - SETTER         - #
    # ------------------ #
    def set_wcs(self,wcs,force_it=False):
        """
        """
        if self.has_wcs() and force_it is False:
            raise AttributeError("'wcs' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        self._side_properties["wcs"] = astrometry.get_wcs(wcs)
        

    def set_catalogue(self,catalogue,force_it=True,
                      default_isolation_def = 10*units.arcsec):
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
        Voids
        """
        if self.has_catalogue() and force_it is False:
            raise AttributeError("'catalogue' already defined"+\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(catalogue) or catalogue.__nature__ != "Catalogue":
            raise TypeError("the input 'catalogue' must be an astrobject catalogue")

        # -------------------------
        # - Add the world_2_pixel
        if not catalogue.has_wcs():
            warnings.warn("WARNING the given 'catalogue' has no pixel coordinates. Cannot load it")
            return
        
        if catalogue.nobjects_in_fov < 1:
            warnings.warn("WARNING No object in the field of view, catalogue not loaded")
            return

        if not catalogue._is_around_defined():
            catalogue.define_around(default_isolation_def)
        self._side_properties["catalogue"] = catalogue
                
    # ------------------ #
    # - GETTER         - #
    # ------------------ #
    def get(self,key,mask=None):
        """
        This function enable to get from the data the values of the given keys
        or derived values, like ellipticity. Set 'help' for help.

        *Remark* 'key' could be an list of keys.

        Returns:
        --------
        array (or list of)
        """
        if not self.has_data():
            raise AttributeError("no 'data' defined")

        if "__iter__" in dir(key):
            return [self.get(key_,mask=mask) for key_ in key]
        
        # -- These are the default key values
        _data_keys_ = self.data.keys()
        _matching_keys_ = ["angsep"]
        _derived_keys_ = ["elongation","ellipticity"]
        help_text = " Known keys are: "+", ".join(_data_keys_+_matching_keys_+_derived_keys_)

        if key in ["help","keys","keylist"]:
            print help_text
            return
        
        # -- These are from the data
        if key in _data_keys_:
            val_ = self.data[key]
        # -- These are the catalogue values
        elif key in _matching_keys_:
            if not self.has_catmatch():
                raise AttributeError("no 'catmatch' defined. The matching has not been ran")
            val_ = self.catmatch[key]
        # -- These are derived values
        elif key in _derived_keys_:
            if key == "elongation":
                val_ = self.get("a") / self.get("b")
            elif key == "ellipticity":
                val_ = 1. - 1. / self.get("elongation")
        else:
            raise ValueError("Cannot parse '%s'."%key +\
                          help_text)
                          
        return val_ if mask is None else val_[mask]

    def get_fwhm_pxl(self,stars_only=True,isolated_only=True,
                    catmag_range=[None,None]):
        """
        This is defined as 2 * sqrt(ln(2) * (a^2 + b^2))
        """
        mask = self.get_indexes(isolated_only=isolated_only,
                                stars_only=stars_only,
                                catmag_range=catmag_range)
        return np.median(2 * np.sqrt( np.log(2) * (self.get('a',mask)**2 + self.get('b',mask)**2)))
    
    # -----------------
    # - get ellipse    
    def get_ellipse_mask(self,width,height, r=3, apply_catmask=False):
        """
        This method returns a boolean mask of the detected ellipses
        on the given width x height pixels image

        (this method is based on the sep mask_ellipse function)
        
        Parameters:
        -----------
        r: [float]                 The scale of the ellipse (r=1 is a typical contours included the
                                   object ; 2 enables to get the tail of most of the bright sources

        apply_catmask: [bool]      Only mask the detected object associated with the current catalogue.
                                   If no catalogue loaded, this will be set to False in any case.

        Returns:
        -------
        2D-bool array (height x width)
        """
        from sep import mask_ellipse
        ellipsemask = np.asarray(np.zeros((height,width)),dtype="bool")
        mask = None if not apply_catmask else self.catmask
        # -- Apply the mask to falsemask
        mask_ellipse(ellipsemask,
                     self.get('x',mask=mask),self.get('y',mask=mask),
                     self.get('a',mask=mask),self.get('b',mask=mask),
                     self.get('theta',mask=mask),
                     r=r)
        
        return ellipsemask

    def get_detected_ellipses(self,scaleup=5,apply_catmask=True,
                stars_only=False, isolated_only=False,
                catmag_range=[None,None]):
        """
        """
        if not self.has_data():
            print "WARNING [Sexobjects] No data to display"
            return
        
        from matplotlib.patches import Ellipse
        
        # -- maskout non matched one if requested
        if not self.has_catalogue():
            apply_catmask = False
            
        mask = None if not apply_catmask else\
          self.get_indexes(isolated_only=isolated_only,stars_only=stars_only,
                        catmag_range=catmag_range)
            
        # -------------
        # - Properties
        return [Ellipse([x,y],a*scaleup,b*scaleup,
                        t*units.radian.in_units("degree"))
                for x,y,a,b,t in zip(self.get("x",mask=mask),self.get("y",mask=mask),
                                     self.get("a",mask=mask),self.get("b",mask=mask),
                                     self.get("theta",mask=mask) )]
    
    def get_median_ellipse(self,mask=None,clipping=[3,3]):
        
        """This methods look for the stars and return the mean ellipse parameters"""
        if not self.has_catalogue():
            apply_catmask = False
          
        # -- apply the masking
        a_clipped,_alow,_ahigh = sigmaclip(self.get("a",mask=mask),*clipping)
        b_clipped,_blow,_bhigh = sigmaclip(self.get("b",mask=mask),*clipping)
        t_clipped,_tlow,_thigh = sigmaclip(self.get("theta",mask=mask),*clipping)
        # - so        
        psf_a,psf_b,psf_t = a_clipped.mean(),b_clipped.mean(),t_clipped.mean()
        m = np.sqrt(len(a_clipped)-1)
        
        return [psf_a,np.std(a_clipped)/m],[psf_b,np.std(t_clipped)/m],\
        [psf_t,np.std(t_clipped)/m]

        
    # ---------------- #
    # - get Mask     - #
    # ---------------- #
    def get_mask(self,isolated_only,stars_only,
                 catmag_range=[None,None]):
        """
        This main masking method builds a boolean array following
        the requested cuts. Remark that this implies doing a catalogue cuts
        since stars and magnitude information arises from the catalogue.

        Returns
        -------
        array (dtype=bool)
        """
        if catmag_range[0] is None:
            warnings.warn("No lower catalogue magnitude limit given. 13mag set.")
            catmag_range[0] = 13
            
        mask = np.asarray(self.get_catmag_mask(*catmag_range))
        if isolated_only:
            mask = mask & self.catisolatedmask
        if stars_only:
            mask = mask & self.catstarmask
        
        return mask

    def get_catmag_mask(self,magmin,magmax):
        """
        return the boolen mask of which matched point
        belong to the given magnitude range.
        Set None for no limit

        (see also get_mask)
        
        Returns
        -------
        array (dtype=bool)
        """
        if not self.has_catalogue():
            raise AttributeError("no 'catalogue' loaded")
        if not self.has_catmatch():
            raise AttributeError("catalogue has not been matched to the data")
    
        mags = self.catalogue.mag[self.catmatch["idx_catalogue"]]
        magmin = np.min(mags) if magmin is None else magmin
        magmax = np.max(mags) if magmax is None else magmax
        return (mags>=magmin) & (mags<=magmax)

    def get_indexes(self,isolated_only=False,stars_only=False,
                    catmag_range=[None,None], cat_indexes=False):
        """
        Converts mask into index. Particularly useful to access the catalogue
        values (set cat_indexes to True)
        
        Returns
        -------
        array of indexes
        """
        id = "idx_catalogue" if cat_indexes else "idx"
        return np.unique(self.catmatch[id][self.get_mask(isolated_only=isolated_only,
                                                  stars_only=stars_only,
                                                  catmag_range=catmag_range
                                                  )])

    def get_photomap(self, matched_only=True,
                     stars_only=False, isolated_only=False):
        """
        """
        print "to be done"

    # ------------------- #
    # - PLOT Methods    - #
    # ------------------- #    
    # - Display
    def display(self,ax,draw=True,
                apply_catmask=True,
                stars_only=False, isolated_only=False,
                catmag_range=[None,None]):
        """
        """
        ells = self.get_detected_ellipses(scaleup=5,apply_catmask=apply_catmask,
                                          stars_only=stars_only,
                                          isolated_only=isolated_only,
                                          catmag_range=catmag_range)
        
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
        if not self.has_catalogue():
            apply_catmask = False
            
        mask = None if not apply_catmask else\
          self.get_indexes(isolated_only=isolated_only,stars_only=stars_only,
                        catmag_range=catmag_range)
        
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
                      apply_catmask=True,stars_only=False,
                      isolated_only=False,catmag_range=[None,None],
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
        if not self.has_catalogue():
            apply_catmask = False
            
        mask = None if not apply_catmask else\
          self.get_indexes(isolated_only=isolated_only,stars_only=stars_only,
                        catmag_range=catmag_range)
        # -------------
        # - Properties
        ells = [Ellipse([0,0],2.,2*b/a,t*units.radian.in_units("degree"))
                for a,b,t in zip(self.get("a",mask=mask),self.get("b",mask=mask),
                                 self.get("theta",mask))]
        # -- Show the typical angle
        psf_a,psf_b,psf_theta = self.get_median_ellipse(mask=mask)
        ellipticity = 1- psf_b[0]/psf_a[0]
        # - cos/ sin what angle in radian
        
        ax.plot([0,np.cos(psf_theta[0])*ellipticity],
                [0,np.sin(psf_theta[0])*ellipticity],ls="-",lw=2,
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
            return Table(sexoutput)
        
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
                       "is_isolated":self.catalogue.isolatedmask[icat] \
                       if self.catalogue._is_around_defined() else None,
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
        return np.asarray(self.wcs.pix2world(self.data["x"],
                                           self.data["y"])).T

    @property
    def xy(self):
        return np.asarray([self.get("x"),self.get("y")])

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
    
