#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module gather the spectral-objects and there associated tools"""
import numpy  as np
import pyfits as pf
from .baseobject import BaseObject


__all__ = ["spectrum","cube"]


def spectrum(filename,**kwargs):
    return Spectrum(filename,**kwargs).copy()

def cube():
    raise NotImplementedError("Not ready yet")


# -------------------- #
# - Inside tools     - #
# -------------------- #
def lbda2headerparameters(lbda):
    """
    Converts the input wavelength-array (*lbda*) into the three fundamental
    header parameters.

    Return
    ------
    npix (int), step (float), start (float)
    """
    npix,start = len(lbda),lbda[0]
    step       = (lbda[-1] - start) / (npix-1.)
    return npix, step, start

def headerparameter2lbda(npix,step,start):
    """
    Creates the wavelength array from its given fundamenetal parameters.
    Remark that *start* and *step* (together) could be given in log scale
    (velocity step).

    Parameters
    ----------
    start: [float]             The first wavelength of the array.
                               (could in log ; in that case *step* must be too)
        
    step: [float]              The step between two wavelength bin.
                               (could in log ; in that case *start* must be too)
                                   
    npix: [int]                The number of wavelength bin.

    Return
    ------
    npix-length 1d-array (wavelength)
    """
    return np.arange(npix) * step + start

#######################################
#                                     #
# Spectro Object Classes              #
#                                     #
#######################################
class Spectrum( BaseObject ):
    """
    """
    _properties_keys = ["y","var","header",
                        "name","filename",
                        "npix","step","start"]
        
    _derived_properties_keys = ["fits","lbda",
                               "raw_lbda"]
        
    _side_properties_keys = ["astrobject"]

    
   

    def __init__(self,filename=None,empty=False,**kwargs):
        """
        """
        self.__build__()
        
        if empty:
            return
        
        if filename is not None:
            self.load(filename,**kwargs)

    def __build__(self):
        """
        """
        super(Spectrum,self).__build__()
         # -- How to read the image
        self._spectro_properties = dict(
            header_step   = "CDELT1",
            header_start  = "CRVAL1",
            header_npix   = "NAXIS1",
            data_index    = 0,
            variance_index= 1, # if any
            )

    # =========================== #
    # = Main Methods            = #
    # =========================== #
    # ----------------------- #
    # - Affect the Spectrum - #
    # ----------------------- #
    def reshape(self,x_model,k=4):
        """
        Change the structure of the object to match the input *x_model*,
        this will notably affect lbda, y, v and the header.
        This function calls set_lbda and the y and v setter.

        Parameters
        ----------

        x_model:[array]            The wavelength array to which the object should be
                                   projected. Important, the step of the array must be
                                   constant. You can give the array in log scale
                                   (velocity_step)

        k: [int]                   Interpolation parameter of scipy's UnivariateSpline
        
        Return
        ------
        Voids
        """
        from ..utils.tools import shape_ajustment
        newy = shape_ajustment(self.lbda,self.y,x_model,k=k)
        newv = np.abs(shape_ajustment(self.lbda,self.v,x_model,k=k)) \
           if self.has_var() else None
           
        self.set_lbda(x_model,check=False)
        self.y = newy
        self.v = newv
        
        
    def shift(self,velocity_km_s,
              keep_original_lbda=False):
        """
        
              
        Parameter
        ---------

        velocity_km_s: [float]

        - option -

        keep_original_lbda: [bool] 
        
        Return
        ------
        Void
        """

    def truncate(self,min_lbda=None,max_lbda=None):
        """The trunctates the wavelength of the spectrum to match the given
        value (*not in log scale*) and will reshape the object consequently"""
        # ---------------
        # - Input test
        if self.lbda is None:
            raise AttributeError("no wavelength (lbda) loaded. Can't truncate it.")
        if min_lbda is None and max_lbda is None:
            print "WARNING both lbda boundaries are None. Nothing to be done"
            return
        # ----------------------
        # - Defines the new lbda
        min_lbda = self.lbda[0] if min_lbda is None else min_lbda
        max_lbda = self.lbda[1] if max_lbda is None else max_lbda
        flaggood = (self.lbda <= max_lbda) & (self.lbda >= min_lbda)
        new_lbda = self._derived_properties["rawlbda"][flaggood]
        # ---------------
        # - And Do it
        self.reshape(new_lbda)
        
    # ----------------------- #
    # - Return a Spectrum   - #
    # ----------------------- #
    def get_reshuffled(self):
        """The will redraw the current spectrum based on the flux and variance
        of the current object. It won't affect this object but will return an
        equivalent one statistically speaking.

        Return
        ------
        Spectrum
        """
        # ---------------
        # - Input test
        if self.y is None:
            raise AttributeError("no flux loaded in the object.")
        if self.has_var() is False:
            raise AttributeError("no variance in the object, so nothing to do.")
        newspec = self.copy()
        newspec.y  = self.y + np.random.normal(loc=self.y, scale=np.sqrt(self.v))
        return newspec
    
    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def load(self,filename,
             index=None,variance_index=None,
             force_it=False):
        """
        This function reads the fitsfile (*filename*) and parse the data to
        then call the 'create' function. 

        Parameters
        ----------
        
        filename: [string.fits]    fitsfile containing the spectral data
             
        - option -
        
        index: [int]               The fits entry where the data are registered.

        variance_index: [int]      The fits entry where the variance or error are
                                   registered.
                                   This will only work if the hdu table VARIANCE (VAR)
                                   or ERROR (ERR) exists.
        
        force_it: [bool]           If the object already is loaded, it could be
                                   hazardous to reload it. Hence this function will
                                   raise an exception in that case except if
                                   *force_it* is set to True.
        Return
        ------
        Void
        """
        # ---------------
        # - Input Test 
        if type(filename) != str:
            raise TypeError("'filename' must be a string")
        if self.y is not None and force_it is False:
            raise AttributeError("The object already is defined."+\
                    " Set force_it to True if you really known what you are doing")
                    
        # -------------------
        # - Parse the input
        index = self._spectro_properties["data_index"] if index is None \
          else index
          
        variance_index = self._spectro_properties["variance_index"] \
          if variance_index is None \
          else variance_index
          
        # ---------------
        # - Data Test                    
        fits   = pf.open(filename)
        data   = fits[index].data
        header = fits[index].header
        if data is None:
            raise ValueError("no 'data' in the given fits file for"+\
                             " the given index (%d)"%index)
        # ----------------- #
        # - Looks good    - #
        # ----------------- #
        # -------------
        # - wavelength
        tmp_lbda = headerparameter2lbda(
            header.get(self._spectro_properties["header_npix"] ),
            header.get(self._spectro_properties["header_step"] ),
            header.get(self._spectro_properties["header_start"]),
            ).copy()

        # -------------
        # - variance
        if len(fits) > 1:
            if fits[variance_index].name in ['VARIANCE','VAR','VARIANCES']:
                variance = fits[variance_index].data
            elif fits[variance_index].name in ['ERROR','ERR','ERRORS']:
                variance = fits[variance_index].data**2
            
        # ------------------------------
        # - Actually Create the Object
        self._derived_properties["fits"]  = fits
        self._properties['filename']      = filename
        self.create(tmp_lbda,data,
                    header=header,variance=variance,
                    name = header["OBJECT"],
                    force_it=True)
        # - Done ! 
        # ----------
        
    # - Main Methods that create the object. The others (load...) end up call it.
    def create(self,lbda,flux,header=None,
               variance=None,name=None,
               force_it=False):
        """
        This creates the object based on the fundamental parameter

        Parameters:
        -----------

        lbda: [array]              The wavelength-array of the object. It could be
                                   given in log scale (velocity_step).

        flux: [array]              The flux of the given object. It will define *y*
        
        - options -

        variance: [array]          The variance associated to the object's flux (*y*).
                                   It will set *v*

        name: [str]                The object's name. If None, this will fetch for the
                                   keywork 'OBJECT' in the header and set name to it if
                                   it finds it.
        
        force_it: [bool]           If the object already is loaded, it could be
                                   hazardous to reload it. Hence this function will
                                   raise an exception in that case except if
                                   *force_it* is set to True.

                                   
        Return
        ------
        Void
        """
        # ----------------
        # - Test the input
        if len(lbda) != len(flux):
            raise ValueError("'flux' and 'lbda' must have the same dimension")

        if variance is not None and len(lbda) != len(variance):
            raise ValueError("If 'variance' is given 'variance' and 'lbda' "+\
                             "must have the same dimension")
        if self.y is not None and force_it is False:
            raise AttributeError("The object already is defined."+\
                    " Set force_it to True if you really known what you are doing")
                    
        # ************************ #
        #  Creation of the Object  #
        # ************************ #
        self._properties['header']  = header
        self._properties['y']       = np.asarray(flux).copy()
        
        self.set_lbda(lbda)
        self._properties['name']    = name
        if variance is not None:
            self._properties['var'] = np.asarray(variance).copy()
        else:
            self._properties['var'] = None
            
        # ------------
        # - Checkout 
        if name is None:
            if "OBJECT" in header.keys():
                self.name = header["OBJECT"]
        
    def writeto(self):
        """
        """
        print "to be done"

    def show(self,savefile=None,ax=None,show=True,
             add_thumbnails=False,**kwargs):
        """
        """
        from ..utils.mpladdon import specplot, figout
        from matplotlib.pyplot import figure
        if self.y is None:
            raise AttributeError("no data to show")

        # --------------
        # - Set the Figure
        self._plot = {}
        if ax is None:
            
            fig = figure(figsize=[8,5])
            ax = fig.add_axes([0.1,0.1,0.8,0.8])
        elif "plot" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. ")
        else:
            fig = ax.figure
            
        # --------------
        # - Do the plot    
        pl,fill = ax.specplot(self.lbda,self.y,var=self.v,**kwargs)
        # --------------
        # - output
        self._plot['ax']     = ax
        self._plot['figure'] = fig
        self._plot['plot']   = pl
        self._plot['band']   = fill
        self._plot['prop']   = kwargs
        
        fig.figout(savefile=savefile,show=show,add_thumbnails=add_thumbnails)
        
        return self._plot
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    # --------------------
    # - Files In
    @property
    def filename(self):
        return self._properties["filename"]

    @property
    def fits(self):
        return self._derived_properties["fits"]

    # --------------------
    # - wavelength
    @property
    def lbda(self):
        return self._derived_properties["lbda"]
    
    def set_lbda(self, x, check=True):
        """
        Change the current wavelength (lbda) to *x*.
        If this object has a flux (y), it will check that the new lbda
        if consistant with the flux (y) struture except if *check* is
        set to False. Most likely you do not want to do that.
        """
        if self.y is not None and check and len(x) != len(self.y):
            raise ValueError("'x' must have the length of self.y (%d)"%len(self.y))

        npix, step, start = lbda2headerparameters(np.asarray(x))
        # -- use the real function.
        self.create_lbda(npix,step,start,force_it=True,check=False)

        
    def create_lbda(self,npix,step,start,
                    force_it=False,check=True):
        """
        Define the basic component of the wavelength. They will be
        dumped in the header and then lbda will be loaded.
        Remark that you can give start and step in log scale (velocity step).
        
        Parameters
        ----------
        start: [float]             The first wavelength of the array.
                                   (could in log ; in that case *step* must be too)
        
        step: [float]              The step between two wavelength bin.
                                   (could in log ; in that case *start* must be too)
                                   
        npix: [int]                The number of wavelength bin.

        - option -

        force_it: [bool]           If the object already is loaded, these parameters
                                   already exist and it could be dangerous to change
                                   them. In that case, this function will raise an
                                   exception except if *force_it* is set to True.
                                   
        check: [bool]              If the object already has a flux (*y*), npix must
                                   corresponds to the size of *y*. If not, this will
                                   raise an exception except if check is set to False.
                                   -> You should not do that.
        
        Return
        ------
        Void
        """
        # ----------------
        # - Input Tests
        if self.lbda is not None and force_it is False:
            raise AttributeError("'lbda' is already defined."+\
                  " Set force_it to True if you really known what you are doing")
            
        if self.y is not None and check and npix != len(self.y):
            raise ValueError("'npix' must have the length of self.y (%d)"%len(self.y))

        # ----------------
        # - Input Tests
        self.header.update(self._spectro_properties['header_npix'],  npix, "")
        self.header.update(self._spectro_properties['header_step'],  step, "")
        self.header.update(self._spectro_properties['header_start'], start,"")
        self._load_lbda_()

        
    # --------------------
    # - flux and variance        
    @property
    def y(self):
        return self._properties['y']

    @y.setter
    def y(self,value):
        if self.lbda is not None and len(value) != self.npix:
            raise ValueError("'y' must have the same length as 'lbda'")
        
        self._properties['y'] = np.asarray(value)
        
    @property
    def v(self):
        return self._properties['var']
    
    @v.setter
    def v(self,value):
        if value is None:
            self._properties['var'] = None
            return
        if self.lbda is not None and len(value) != self.npix:
            raise ValueError("'v' must have the same length as 'lbda'")
        self._properties['var'] = np.asarray(value)
        
    def has_var(self):
        return True if self.v is not None \
          else False
          
    # --------------------
    # - Header
    @property
    def header(self):
        return self._properties['header']

    @property
    def name(self):
        return self._properties['name']
    
    @name.setter
    def name(self,newname):
        if type(newname) is not str:
            raise TypeError("the given 'name' must be a string")
        self._properties['name'] = newname
        
    # --------------------
    # - Inner properties
    @property
    def npix(self):
        return self._properties["npix"]
    @property
    def step(self):
        return self._properties["step"]
    @property
    def start(self):
        return self._properties["start"]

    def has_velocity_step(self,min_linear_start=100):
        if self.start is None:
            raise AttributeError("no *start* defined (so no lbda)")
        
        return True if self.start <= min_linear_start \
          else False
          
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _update_(self):
        """Loads the derived properties based on the fundamental ones"""
        if self.header is not None:
            self._load_lbda_()

        
    def _load_lbda_(self):
        """This function reads the header to create the lbda-array"""
        if self.header is None:
            raise AttributeError("no header loaded, can't load the wavelength")
        
        self._properties["npix"] = \
          self.header.get(self._spectro_properties['header_npix'])
        self._properties["step"] = \
          self.header.get(self._spectro_properties['header_step'])
        self._properties["start"] = \
          self.header.get(self._spectro_properties['header_start'])

        # -- This is the wavelength given by the header.
        #    Could be in velocity step.
        self._derived_properties['rawlbda'] = \
          np.arange(self.npix) * self.step + self.start

        if self.has_velocity_step():
            self._derived_properties['lbda'] = \
              np.exp(self._derived_properties['rawlbda'])
        else:
            self._derived_properties['lbda'] = \
              self._derived_properties['rawlbda']

    
class Cube( BaseObject ):
    """
    """
    _properties = dict(
        lbda= None,
        data= None,
        var = None,
        header= None,
        name= None,
        )
    
    _side_properties = dict(
        astroobject= None # an astrobject
        )
    
    _derived_properties = dict(
        )

    
    def __init__(self):
        """
        """
        raise NotImplementedError("Ongoing")
