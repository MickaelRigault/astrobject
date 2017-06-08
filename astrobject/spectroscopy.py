#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module gather the spectral-objects and there associated tools"""

import warnings
import numpy  as np
import matplotlib.pyplot as mpl

from astropy.io import fits as pf

from .baseobject import TargetHandler
from .utils.tools import shape_ajustment

__all__ = ["get_spectrum"]


def get_spectrum(filename, **kwargs):
    """ create a spectrum having wavelength, flux, variance
    and associated methods

    Parameters
    ----------
    filename: [string]
        location of the spectrum fits-file to open.

    Return
    ------
    Spectrum
    """
    return Spectrum(filename,**kwargs).copy()    

def get_cube():
    raise NotImplementedError("Not ready yet")


# -------------------- #
# - Useful tools     - #
# -------------------- #

    
def merge_spectra(spectrum1,spectrum2):
    """This will merge the two given spectra. They must be Spectrum-object"""
    #
    # Tested with variance
    #
    new_step  = np.min([spectrum1.step,  spectrum2.step ])
    new_start = np.min([spectrum1.start, spectrum2.start ])
    new_end   = np.max([np.max(spectrum1.lbda), np.max(spectrum2.lbda) ])
    # - The new lbda
    new_lbda  = np.arange(new_start,new_end,new_step)

    # - Copy the spectrum not to affect them
    copyspec1 = spectrum1.get_reshaped(new_lbda)
    copyspec2 = spectrum2.get_reshaped(new_lbda)

    # - weight merging
    if copyspec1.has_var() is False or copyspec2.has_var() is False:
        merging_w1,merging_w2 = np.ones(copyspec1.npix), np.ones(copyspec2.npix)
        both_has_var = False
    else:
        merging_w1,merging_w2 = 1./copyspec1.v, 1./copyspec2.v
        both_has_var = True
    # - clean the nans
    # nanflags
    flagout1  = (copyspec1.lbda > spectrum1.lbda[-1]) +\
       (copyspec1.lbda < spectrum1.lbda[0])
    flagout2  = (copyspec2.lbda > spectrum2.lbda[-1]) + \
      (copyspec2.lbda < spectrum2.lbda[0])
    flagout1 = (copyspec1.y != copyspec1.y) + flagout1
    flagout2 = (copyspec2.y != copyspec2.y) + flagout2
    # technic, data = 1 ; weight = 0
    copyspec1.y[flagout1] = 1 ; merging_w1[flagout1] = 0
    copyspec2.y[flagout2] = 1 ; merging_w2[flagout2] = 0
    # clean total weight
    total_w = merging_w1+merging_w2
    flagzero = (total_w<=0)
    total_w[flagzero] = np.NaN
    # -----------------
    # - actual merging
    merged_data = (copyspec1.y*merging_w1 + copyspec2.y*merging_w2) / total_w
    merged_var = 1/np.sum([merging_w1,merging_w2], axis=0) if both_has_var \
        else None
    merged_spec = Spectrum(empty=True)
    merged_spec.create(new_lbda,merged_data,
                        header=copyspec1.header,
                        variance=merged_var,
                        force_it=True)
    return merged_spec
        
    
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
class BaseSpectrum( TargetHandler ):
    PROPERTIES = ["wave","flux","errors"]

    # =================== #
    #   Main Methods      #
    # =================== #
    def __init__(self, wave=None, flux=None, errors=None, astrotarget=None):
        """ """
        self.__build__()
        if wave is not None:
            self.set_wave(wave)
        if flux is not None:
            self.set_flux(flux, errors)
        if astrotarget is not None:
            self.set_target(astrotarget)
        
    # --------- #
    # SETTER    #
    # --------- #
    def set_wave(self, wave, force_it=False):
        """ Sets the wavelength. 
        If a flux is already set, this is test if the size of the 
        given `wave` matches that of the flux. 
        Set force_it to True to ignore that.
        
        """
        if self.has_flux() and len(wave) != len(self.flux) and not force_it:
            raise ValueError("the size given wavelength does not match the current flux. Set force_it to True to ignore this test")
        self._properties['wave'] = np.asarray(wave)
          
    def set_flux(self, flux, errors=None, force_it=True):
        """ Sets the flux of the spectrum """
        if self.wave is not None and len(flux) != len(self.wave) and not force_it:
            raise ValueError("the size given flux does not match the current wavelength. Set force_it to True to ignore this test")
        if errors is not None and len(flux) != len(errors):
            raise ValueError("flux and errors must have the same dimension")
        
        self._properties['flux'] = np.asarray(flux)
        self._properties['errors'] = np.asarray(errors) \
          if errors is not None else None
    
    # --------- #
    # PLOTTER   #
    # --------- #
    def show(self, ax=None, err_onzero=False, bandprop={},
                 savefile=None, show=True, **kwargs):
        """ """
        from astrobject.utils.mpladdon import specplot, figout
        from matplotlib.pyplot import figure
        self._plot = {}
        if ax is None:
            fig = mpl.figure(figsize=[8,5])
            ax  = fig.add_axes([0.12,0.12,0.8,0.8])
        else:
            fig = ax.figure

        
        pl = ax.specplot(self.wave if self.wave is not None else np.arange(len(self.flux)),
                             self.flux, var=self.errors**2 if self.has_errors() else None,
                        err_onzero=err_onzero, bandprop=bandprop,
                        **kwargs)
        self._plot['ax']     = ax
        self._plot['figure'] = fig
        self._plot['plot']   = pl
        self._plot['band']   = bandprop
        self._plot['prop']   = kwargs
        
        fig.figout(savefile=savefile, show=show)
        
        return self._plot
    
    # =================== #
    #   Properties        #
    # =================== #
    @property
    def flux(self):
        """ Observed flux """
        return self._properties['flux']
    
    @property
    def errors(self):
        """ Error on the observed flux """
        return self._properties['errors']

    def has_flux(self):
        """ Has the flux value been set? """
        return self.flux is not None
    
    def has_errors(self):
        """ Have the errors been set? """
        return self.errors is not None
    
    @property
    def wave(self):
        """ wavelength solution """
        return self._properties['wave']

    
class Spectrum( BaseSpectrum ):
    """
    """
    PROPERTIES         = ["header","name","filename","npix","step","start"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["fits","lbda","raw_lbda"]

    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,filename=None,empty=False,**kwargs):
        """
        """
        self.__build__()
        
        if empty:
            return
        
        if filename is not None:
            self.load(filename,**kwargs)

    def __build__(self):
        """ build the object's structure """
        """
        self._properties_keys += ["y","var","header",
                                "name","filename",
                                "npix","step","start"]
        self._side_properties_keys += ["target"]        
        self._derived_properties_keys += ["fits","lbda",
                               "raw_lbda"]
        """
        super(Spectrum,self).__build__()
        
         # -- How to read the image
        self._build_properties = dict(
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
        newy = shape_ajustment(self.lbda,self.flux,x_model,k=k)
        newe = np.abs(shape_ajustment(self.lbda,self.errors,x_model,k=k)) \
           if self.has_errors() else None
           
        self.set_lbda(x_model,check=False)
        self.set_flux(newy, newe)
        
        
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
        print "To Be Done"
        
    def truncate(self,min_lbda=None,max_lbda=None):
        """
        The trunctates the wavelength of the spectrum to match the given
        value (*not in log scale*) and will reshape the object consequently
        """
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
    def get_reshaped(self,x_model,**kwargs):
        """This will copy the current object and reshape it following the
        given x_model. This calls 'reshape' to do so.

        Parameter
        ---------

        x_model:[array]            The wavelength array to which the object should be
                                   projected. Important, the step of the array must be
                                   constant. You can give the array in log scale
                                   (velocity_step)

        - kwargs options ; potentially non-exhaustive ; goes to 'reshape' -

        k: [int]                   Interpolation parameter of scipy's UnivariateSpline

        Return
        ------
        Spectrum
        """
        
        newspec = self.copy()
        newspec.reshape(x_model,**kwargs)
        return newspec
        
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
        if self.has_flux():
            raise AttributeError("no flux loaded in the object.")
        if self.has_errors() is False:
            raise AttributeError("no variance in the object, so nothing to do.")
        newspec = self.copy()
        newspec.set_flux(self.flux + np.random.normal(loc=self.flux, scale=self.errors))
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
        if self.has_flux() and force_it is False:
            raise AttributeError("The object already is defined."+\
                    " Set force_it to True if you really known what you are doing")
                    
        # -------------------
        # - Parse the input
        index = self._build_properties["data_index"] if index is None \
          else index
          
        variance_index = self._build_properties["variance_index"] \
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
            header.get(self._build_properties["header_npix"] ),
            header.get(self._build_properties["header_step"] ),
            header.get(self._build_properties["header_start"]),
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
                    name = header.get("OBJECT","unknown") if header is not None\
                     else "unknown",
                    force_it=True)
        # - Done ! 
        # ----------
        
    # - Main Methods that create the object. The others (load...) end up call it.
    def create(self, lbda, flux, header=None,
               variance=None, x_model=None,
               name=None, force_it=False, 
               fixformat=True, k=4,**kwargs):
        """
        This creates the object based on the fundamental parameter

        Parameters:
        -----------

        lbda: [array]
            The wavelength-array of the object. It could be
            given in log scale (velocity_step).
            The step should be constant. If not, it will be reformated
            to be so (warning raised).
            
        flux: [array]
            The flux of the given object. It will define *y*


        header: [pyfits.Header] -optional-
            Give the header of the object. This will be recorded
            when you save the object ('writeto'). If None, this
            will be converted in an empty header.
        
        variance: [array] -optional-
            The variance associated to the object's flux (*y*).
            It will set *v*

        x_model: [array] -optional-
            Force the shape of the wavelength to be the given x_model.

        name: [str] -optional-
            The object's name. If None, this will fetch for the
            keywork 'OBJECT' in the header and set name to it if
            it finds it.
        
        force_it: [bool] -optional-
            If the object already is loaded, it could be hazardous to reload it.
            Hence this function will raise an exception in that case except if
            *force_it* is set to True.

        fixformat: [bool] -optional-
            If the given lbda does not have a constant step, this method will
            create a new lbda with a constant step (and correct the flux and variance
            correspondingly) except if fixformat is set to False. In that case, a
            TypeError will be raised.

        **kwargs goes to the shape_adjustment function (so to scipy's UniverseSpline:
        k and s parameter e.g.)
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
        if self.has_flux() and force_it is False:
            raise AttributeError("The object already is defined."+\
                    " Set force_it to True if you really known what you are doing")

        # =============== #
        # =  Modeling   = #
        # =============== #
        if x_model is not None:
            if np.std(x_model[:-1]-x_model[1:])>1e-1:
                if fixformat:
                    warnings.warn("Non constant step on the given spectrum. Create one as such.")
                    x_model = np.linspace(x_model[0],x_model[-1], len(x_model))
                else:
                    raise TypeError("The given wavelength array does not have a constant step.")
                
                
        elif np.std(lbda[:-1]-lbda[1:])>1e-1:
            if fixformat:
                warnings.warn("Non constant step on the given spectrum. Create one as such.")
                x_model = np.linspace(lbda[0],lbda[-1], len(lbda))
            else:
                raise TypeError("The given wavelength array does not have a constant step.")
        
        if x_model is not None: # test if similar step    
            flux = shape_ajustment(lbda, flux, x_model, **kwargs)
            variance = np.abs(shape_ajustment(lbda, variance, x_model, **kwargs)) \
              if variance is not None else None
            lbda = x_model
            
        # ************************ #
        #  Creation of the Object  #
        # ************************ #
        # -- Header 
        self._properties['header'] = pf.Header() if header is None else header
        # -- Wavelength
        self.set_lbda(lbda)
        # -- Data
        self.set_flux(flux, np.sqrt(variance) if variance is not None else None)
        # -- Name
        self.name = self.header["OBJECT"] if name is None and "OBJECT" in self.header.keys() \
          else str(name) if name is not None else "unknown"

                
    def writeto(self,savefile,saveerror=False,
                overwrite=False):
        """
        save the object in a fitsfile (*savefile*).

        Parameters
        ----------
        savefile: [str]            The name (including its Path) of
                                   the file where you wish to register
                                   the spectrum (.fits)                                   
        - options -
                                       
        saveerror:  [bool]         Set this to True if you wish to record the error
                                   and not the variance in you first hdu-table.
                                   if False, the table will be called VARIANCE and
                                   have self.errors**2; if True, the table will be called
                                   ERROR and have self.errors
                                   
        overwrite: [bool]          If the file already exist, saving it in the same
                                   file will overwrite it. Set True to allow that.


        Return
        ------
        Void 
        """
        hdu     = self._get_hdu_array_(use_error=saveerror)
        hdulist = pf.HDUList(hdu)
        hdulist.writeto(savefile,clobber=overwrite)
        
        
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
        If this object has a flux (y), this method will check that
        the new lbda (*x*) if consistant with the struture of the
        flux (*y*)  except if *check* is set to False.
        (Most likely you do not want to do that.)

        Parameter
        ---------

        x: [array]                 The new wavelength of the object. (1d array)
                                   A log of wavelength array can be set (velocity_step)

        - options -

        check: [bool]              If True, this will check tha the new wavelength *x*
                                   is consistant with the existing flux (*y*) if any.
                                   Set that to False to avoid this check.
        Return
        ------
        Void
        """
        if self.has_flux() and check and len(x) != len(self.flux):
            raise ValueError("'x' must have the length of self.flux (%d)"%len(self.flux))

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
            
        if self.has_flux() and check and npix != len(self.flux):
            raise ValueError("'npix' must have the length of self.flux (%d)"%len(self.flux))

        # ----------------
        # - Input Tests
        self.header.set(self._build_properties['header_npix'],  npix, "")
        self.header.set(self._build_properties['header_step'],  step, "")
        self.header.set(self._build_properties['header_start'], start,"")
        self._load_lbda_()

    @property
    def wave(self):
        """ pointer towards self.lbda """
        return self.lbda
    
    def set_wave(self,*args, **kwargs):
        raise NotImplementedError("use set_lbda")
    
    # --------------------
    # - flux and variance        
    @property
    def y(self):
        print "DECREPATED, use self.flux"
        return self.flux
        
    @property
    def v(self):
        print "DECREPATED, use self.errors**2"
        return self.errors**2 if self.has_errors() else None
        
    def has_var(self):
        print "DECREPATED, use self.has_errors()"
        return self.has_errors()
          
    # --------------------
    # - Header
    @property
    def header(self):
        if self._properties["header"] is None:
            self._properties["header"] = pf.Header()
        return self._properties['header']
        
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

    def _get_hdu_array_(self,use_error):
        """
        This method create the current object's hdu containing
        the primary-hdu (self.flux) and the error/variance if it exists.
        Do pyfits.HDUList(return_of_this) to have an hdulist.
        
        Parameter
        ---------

        use_error: [bool]          Shall this register the error (set True) or
                                   the variance (set False) in the *hdu*.

        Return
        ------
        pyfits hdu 2darray (primary, [extension])
        """
        self._hdu = [pf.PrimaryHDU(self.flux,self.header)]
        # - If there is error, saveit
        if self.has_var():
            if use_error:
                extention = pf.ImageHDU(self.errors, name='ERROR')
            else:
                extention = pf.ImageHDU(self.errors**2, name='VARIANCE')
            self._hdu.append(extention)
            
        # -- Return the current object's hdu
        return self._hdu
            
        
    def _load_lbda_(self):
        """This function reads the header to create the lbda-array"""
        if self.header is None:
            raise AttributeError("no header loaded, can't load the wavelength")
        
        self._properties["npix"] = \
          self.header.get(self._build_properties['header_npix'])
        self._properties["step"] = \
          self.header.get(self._build_properties['header_step'])
        self._properties["start"] = \
          self.header.get(self._build_properties['header_start'])

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

