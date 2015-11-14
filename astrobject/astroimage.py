#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy  as np
import pyfits as pf

from .baseobject import BaseObject 
from astropy     import units
from astLib      import astWCS

from ..utils.tools import kwargs_update
from ..utils.decorators import _autogen_docstring_inheritance


__all__ = ["image","aperture"]


def image(filename,**kwargs):
    """
    Initalize the image by giving its filelocation (*filename*). This
    will load it using the load() method.

    Parameters
    ----------

    filename: [string.fits]    fits file from where the image will be loaded

    - options -

    - kwargs options, potentially non-exhaustive -
    
    empty: [bool]              Set True to load an empty object.
                               
    Return
    ------
    Image
    """
    return Image(filename,**kwargs).copy()
    

def aperture(astrobject_,filename,**kwargs):
    """
    Initalize the image by giving its filelocation (*filename*). This
    will load it using the load() method.

    Parameters
    ----------

    astrobject_: [AstroObject] An astrobject associated to the given image.
    
    filename: [string.fits]    fits file from where the image will be loaded

    - options -

    - kwargs options, potentially non-exhaustive -
    
    empty: [bool]              Set True to load an empty object.
                               
    Return
    ------
    Image
    """
    return Aperture(astrobject_,filename,**kwargs).copy()



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
    _properties_keys         = ["filename","data","header",
                                "variance"]
    _side_properties_keys    = ["wcs"]
    _derived_properties_keys = ["fits"]
    
    # Where in the fitsfile the data are
    _image_property = dict(
        data_index = 0
        )
    
    # =========================== #
    # = Initialization          = #
    # =========================== #
    def __init__(self,filename,empty=False,**kwargs):
        """
        Initalize the image by giving its filelocation (*filename*). This
        will load it using the load() method.

        Parameters
        ----------
        filename: [string.fits]    fits file from where the image will be loaded
                                   - Trick - Set None and no image will be loaded
        
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
            
        
    def load(self,filename,index=_image_property["data_index"],
             force_it=False):
        """
        This enables to load a fitsfile image and will create
        the basic data and wcs solution if possible.
        *Variance* (error) and *background* has to be defined
        separately has this strongly depend on the instrument

        Parameters
        ----------
        filename: [string.fits]    The file containing the fits data

        - options -
        
        index: [int]               The fits entry that contains the data

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
            
        self.define(data,None,fits[index].header,wcs_,
                    filename,fits)
        
    def define(self,data,variance,header,wcs,
               filename,fits):
        """
        Define image using by filling its different component. Each of them
        can be None, but in that case, their corresponding method and derived
        values won't be avialable.
        """
        self._properties["filename"]     = filename
        self._derived_properties["fits"] = fits
        # basics data and wcs
        self._properties["data"]         = np.asarray(data,dtype="float")
        self._properties["header"]       = header
        
        if wcs is not None and wcs.__module__ != "astLib.astWCS":
            print "WARNING: only astLib.astWCS wcs solution is implemented"
            print " ----> No wcs solution loaded"
            self._side_properties["wcs"]     = None
        else:
            self._side_properties["wcs"]     = wcs

        self._update_()
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    # ------------------- #
    # - get Methods     - #
    # ------------------- #
    def get_pixels_around(self,pixel_x,pixel_y,radius_in_pixel):
        """
        """
        # - This should be moved in a c-lib
        weight = np.zeros(self.shape)
        
        
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
    # - I/O Methods     - #
    # ------------------- #
    def writeto(self):
        """
        """
        print "to be done"

    def show(self,logscale=True,ax=None,show=True,
             wcs_coords=False,
             **kwargs):
        """
        """
        if self.data is None:
            raise AttributeError("no 'data' to show")
        
        # -- Setting -- #
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax = fig.add_axes([0.1,0.1,0.8,0.8])
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")

        # ----------- #
        # -  What
        if logscale:
            x = np.log10(self.data)
        else:
            x = self.data
            
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

        # ----------- #
        # - Recordit
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["imshow"] = im
        self._plot["prop"]   = prop
        self._plot["wcs_coords"] = wcs_coords
        if show:
            fig.show()
        
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
        return self._properties["data"]

    @property
    def header(self):
        return self._properties["header"]
    
    # -- derived values
    @property
    def shape(self):
        if self.data is None:
            raise AttributeError("No data loaded.")
        return np.shape(self.data)
    
    @property
    def width(self):
        return self.shape[1]
    @property
    def height(self):
        return self.shape[0]

    
    @property
    def variance(self):
        return self._properties["variance"]
    def has_variance(self):
        return False if self.variance is None \
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

    @property
    def _worldcoords_bbox(self):
        from matplotlib.transforms  import Bbox
        xmin,xmax,ymin,ymax = self.worldcoords_boundaries
        bbox= Bbox([[xmin,ymin],[xmax,ymax]])
        return bbox.rotated(self.wcs.getRotationDeg()*units.degree.in_units("radian"))
        
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    



    
class Aperture( Image ):
    """
    An aperture is an Image connected to an AstroObject that has
    a position in a sky and a redshift. This combination enables
    to derive aperture photometry around this AstroObject. 
    """
    
    @_autogen_docstring_inheritance(Image.__init__,"Image.__init__")
    def __init__(self,astrobject,*args,**kwargs):
        #
        # - Add the astrobject associated tools
        #
        super(Aperture,self).__init__(*args,**kwargs)
        self._properties_keys.append("object")
        self.change_object(astrobject)


    # =========================== #
    # = Image Hack              = #
    # =========================== #
    @_autogen_docstring_inheritance(Image.show,"Image.show")
    def show(self,markerprop={},**kwargs):
        #
        # Add object indication
        #
        show = kwargs.pop("show",True)
        super(Aperture,self).show(show=False,**kwargs)
        # --- Can you show the object location ?
        if self.has_wcs() is False:
            print "WARNING: no wcs solution so I cannot display the object location"
            if show:
                self._plot["figure"].show()
            return
        # ------------
        # - Ok... then
        # ---------------------
        # - get back the values
        fig,ax = self._plot["figure"],self._plot["ax"]

        # ----------- #
        # - Fancy
        default_markerprop = {
            "marker":"s",
            "mfc":"w",
            "mec":"k","mew":2,
            "zorder":12,
            "scalex":False,"scaley":False
            }
        prop = kwargs_update(default_markerprop,**markerprop)
    
        if self._plot["wcs_coords"]:
            self._plot["object_marker"]      = \
                ax.plot(self.object.ra,self.object.dec,**prop)
        else:
            radec_pixel = self.coords_to_pixel(*self.object.radec)
            self._plot["object_marker"]      = \
                ax.plot(radec_pixel[0],radec_pixel[1],**prop)
                
        self._plot["object_marker_prop"] = prop
        
        if show:
            fig.show()
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def object(self):
        return self._properties["object"]
    
    def change_object(self,newastrobject,test_inclusion=True):
        """
        Change the object associated to the given image. This function
        will test if the object is withing the image boundaries (expect if
        *test_inclusion* is set to False).
        """
        # -- Input Test -- #
        # - Is that even an AstroObject
        if newastrobject.__nature__ != "AstroObject":
            raise TypeError("'newastrobject' should be (or inherite) an AstroObject")
        
        # - Is that a relevant one ?
        if test_inclusion:
            if self.has_wcs() is False:
                print "WARNING: because there is no wcs solution, "+\
                  "I can't test the inclusion of the new astrobject"
            else:
                if not self.wcs.coordsAreInImage(*newastrobject.radec):
                    raise ValueError("The new 'astrobject' is not inside the image "+\
                                      " boundaries"+ "\n"+\
                                     "--> object radec: %.3f,%.4f"%(newastrobject.ra,
                                                                    newastrobject.dec))
        # -- Seems Ok -- #
        self._properties["object"] = newastrobject.copy()
        
