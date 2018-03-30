#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np


# - astropy
from astropy import coordinates,units
from astropy.table.table import Table,TableColumns, Column
from astropy.io import fits as pf
from astropy.io import ascii

# - local
from .. import astrometry
from ..photometry import Image, get_photopoint
from ..baseobject import BaseObject, WCSHandler
from ..utils.decorators import _autogen_docstring_inheritance
from ..utils.tools import kwargs_update, mag_to_flux, load_pkl, dump_pkl, is_arraylike
from ..utils       import shape

__all__ = ["Instrument"]

class Instrument( Image ):
    """
    """
    PROPERTIES         = ["bandname"]
    DERIVED_PROPERTIES = ["bandpass"]
    
    instrument_name = "TO_BE_DEFINED"
    INFO            = {} # set one
    
    def __build__(self,data_index=0):
        """This is a slightly advanced Image object"""
        super(Instrument,self).__build__(data_index=data_index)

    # ---------------- #
    #  PhotoPoints     #
    # ---------------- #
    
    def _aperture_to_photopoint_(self,count,err,flag):
        """ convert the aperture output to a photopoints """
        
        flux = self.count_to_flux(count)
        var  = self.count_to_flux(err)**2
        
        # ------------------
        # - One Photopoint
        if not is_arraylike(flux):
            return get_photopoint(lbda=self.lbda, flux=flux, var=var,
                            source="image",mjd=self.mjd,
                            zp=self.mab0,bandname=self.bandname,
                            instrument_name=self.instrument_name)
        # -----------------------
        # - Several Photopoints
        return [get_photopoint(lbda=self.lbda,flux=flux_,var=var_,
                            source="image",mjd=self.mjd,
                            zp=self.mab0,bandname=self.bandpass.name,
                            instrument_name=self.instrument_name)
                            for flux_,var_ in zip(flux,var)]
    
        
    @_autogen_docstring_inheritance(Image.get_aperture,"Image.get_aperture")
    def get_photopoint(self,x, y, radius=None, runits="pixels",
                       aptype="circle", wcs_coords=False, getlist=False,
                       **kwargs):
        #
        # Returns a PhotoPoint
        #
        pp = self._aperture_to_photopoint_(*self.get_aperture(x,y,radius=radius,runits=runits,
                                                            aptype=aptype,wcs_coords=wcs_coords,
                                                            **kwargs))
        # ------------------
        # - One Photopoint
        if not is_arraylike(pp) or len(pp)==1 or getlist:
            return pp
        
        # -----------------------
        # - Several Photopoints\
        from ..collections import get_photomap
        return get_photomap(pp,np.asarray([x,y]).T,
                        wcs=self.wcs,
                        catalogue=self.catalogue.get_subcatalogue(fovmask=True, catmag_range=[1,30]) \
                          if self.has_catalogue() else None,
                        wcs_coords=wcs_coords)
    
    @_autogen_docstring_inheritance(Image.get_target_aperture,
                                    "Image.get_target_aperture")
    def get_target_photopoint(self,radius=None, runits="pixels",
                              aptype="circle", **kwargs):
        #
        # Returns a PhotoPoint
        #
        if not self.has_target():
            return AttributeError("No 'target' loaded")
        
        xpix,ypix = self.coords_to_pixel(self.target.ra,self.target.dec)
        pp = self.get_photopoint(xpix,ypix,radius=radius,runits=runits,
                                   aptype="circle",**kwargs)
        pp.set_target(self.target)
        return pp
    
    @_autogen_docstring_inheritance(Image.get_host_aperture,
                                    "Image.get_host_aperture")
    def get_host_photopoint(self,scaleup=2.5,**kwargs):
        #
        # Be returns a PhotoPoint
        #
        if not self.has_target():
            return AttributeError("No 'target' loaded")
        ap_output = self.get_host_aperture(scaleup=scaleup,**kwargs)
        if ap_output is None:
            return None
        
        pp = self._aperture_to_photopoint_(*ap_output)
        pp.set_target(self.target)
        return pp

    
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def count_to_flux(self,counts):
        """ converts counts into flux """
        return counts * 10**(-(2.406+self.mab0) / 2.5 ) / (self.lbda**2)

    # =========================== #
    # = Internal Catalogue      = #
    # =========================== #
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    
    # ------------------
    # - Band Information
    @property
    def bandname(self):
        if self._properties['bandname'] is None:
            raise NotImplementedError("'band' must be implemented")
        return self._properties['bandname']
    
    def set_bandname(self, value, reset_bandpass=True):
        """ Change the name of the bandname.
        
        """
        if value is not None:
            if type(value) != str and type(value) != np.string_:
                raise TypeError("The bandname must be a string", type(value))
        
        self._properties["bandname"] = value
        if reset_bandpass:
            self._derived_properties["bandpass"] = None
        
    @property
    def bandpass(self):
        """ Object containing the basic information associated to the bandname """
        # - No bandname ?
        if self._derived_properties["bandpass"] is not None:
            return self._derived_properties["bandpass"]
        
        if self.bandname is None:
            raise AttributeError("No bandname given")
        
        # - Should this use sncosmo
        try:
            from sncosmo import get_bandpass
            has_sncosmo = True
        except ImportError:
            warnings.warn("sncosmo is not installed. Could not access the bandpass")
            has_sncosmo = False
            
        if has_sncosmo:
            try:
                self._derived_properties["bandpass"] = get_bandpass(self.bandname)
                use_default_bandpass = False
            except:
                use_default_bandpass = True
        else:
            use_default_bandpass = True
            
        if use_default_bandpass:
            wave_eff = self.INFO[self.bandname]["lbda"] \
              if self.bandname in self.INFO else np.NaN
            self._derived_properties["bandpass"] = _DefaultBandpass_(self.bandname, wave_eff)
            
        return self._derived_properties["bandpass"]
    
    @property
    def lbda(self):
        """ effective wavelength """
        return self.bandpass.wave_eff

    # -- Derived values
    @property
    def mab0(self):
        raise NotImplementedError("'mab0' must be implemented")

    @property
    def _gain(self):
        raise NotImplementedError("'_gain' must be implemented (even for None)")

    @property
    def _dataunits_to_electron(self):
        """The gain converts ADU->electron. The Data shouls be in ADU/s"""
        return self._gain * self.exposuretime
    
    @property
    def mjd(self):
        raise NotImplementedError("'mjd' (Modified Julian Date) must be implemented")

############################################
#                                          #
# Base Instrument: CATALOGUE               #
#                                          #
############################################

class Catalogue( WCSHandler ):
    """
    """
    __nature__ = "Catalogue"
    source_name = "_not_defined_"

    PROPERTIES         = ["filename","data","header"]
    SIDE_PROPERTIES    = ["fovcontours","fovmask","matchedmask",
                          "lbda","excluded_list"]
    DERIVED_PROPERTIES = ["fits","naround","naround_nofovcut","contours"]


    def __init__(self, catalogue_file=None,
                 data_index=0,
                 key_mag=None,key_magerr=None,
                 key_ra=None,key_dec=None, empty=False):
        """ initialize the calalogue object

        Parameters
        ----------
        
        catalogue_file: [string] -optional-
            location of a catalogue file.

        empty: [bool] -optional-
            get an empty object.
            
        // settings

        data_index: [int] -optional-
            the fits file entry where the header containing the data
            is recorded.

        key_mag, key_magerr: [string, string] -optional-
            provide the keys associated to the magnitude of the object
            in the given catalogue.
            This is necessary to access 'mag'

        key_ra, key_dec: [string, string] -optional-
            proivde the keys associated to the location (ra and dec)
            of the entries of the catalogue.
            There is a default search (for DE/DEC etc.) that should catch
            any typical case. If not, you should provide it here.

        
        """
        self.__build__(data_index=data_index,
                       key_mag=key_mag,key_magerr=key_magerr,
                       key_ra=key_ra,key_dec=key_dec)
        if empty:
            return
        
        if catalogue_file is not None:
            self.load(catalogue_file)

    def __build__(self,data_index=0,
                  key_mag=None,key_magerr=None,
                  key_ra=None,key_dec=None, key_id=None):
        """
        """
        super(Catalogue,self).__build__()
        self._build_properties = dict(
            data_index = data_index,
            key_id=key_id
            )
        self.set_mag_keys(key_mag,key_magerr)
        self.set_coord_keys(key_ra,key_dec)
            
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    
    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def load(self,catalogue_file,**kwargs):
        """ load the given file and create the object.
        
        kwargs can have any build option like key_ra, key_mag etc.
        """
        # ---------------------
        # - Parsing the input
        if catalogue_file.endswith(".fits"):
            # loading from fits file
            fits   = pf.open(catalogue_file)
            header = fits[self._build_properties["data_index"]].header
            data   = fits[self._build_properties["data_index"]].data
            if type(data) == pf.fitsrec.FITS_rec:
                from astrobject.utils.tools import fitsrec_to_dict
                data = TableColumns(fitsrec_to_dict(data))
                
        elif catalogue_file.endswith(".pkl"):
            # loading from pkl
            fits = None
            header = None
            data = load_pkl(catalogue_file)
            if not type(data) is Table:
                try:
                    data = Table(data)
                except:
                    warnings.warn("Convertion of 'data' into astropy Table failed")
        else:
            fits   = None
            header = None
            format_ = kwargs.pop("format","ascii")
            data   = Table.read(catalogue_file,format=format_,**kwargs)
            
        # ---------------------
        # - Calling Creates
        self.create(data, header, **kwargs)
        self._properties["filename"] = catalogue_file
        self._derived_properties["fits"] = fits

        
    def create(self,data,header,force_it=True,**build):
        """ builds the catalogue

        Parameters
        ----------

        data: [astropy.TableColumns, dictionnary or numpy.ndarray (withkey)]
            the data associated to the catalogue. It must have
            ra, dec, and magnitude entries. The keys associated
            to these are set in _build_properties (key_ra,
            key_dec, key_mag ...). See also set_mag_keys

        header: [pyfits.Header / None]
            Header containing the data for the catalogue (if any)


        force_it: [bool] -optional-
            if data already exists, set force_it to true to overwrite it.

        **build goes to the build dictorty (key_mag, data_slice etc.)
        
        Returns
        -------
        Void
        """
        if self.has_data() and force_it is False:
            raise AttributeError("'data' is already defined."+\
                    " Set force_it to True if you really known what you are doing")
    
        self._properties["data"] = Table(data)
        self._properties["header"] = header if header is not None \
          else pf.Header()
        self.set_starsid(build.pop("key_class",None),build.pop("value_star",None))
        self._build_properties = kwargs_update(self._build_properties,**build)
        # -------------------------------
        # - Try to get the fundamentals
        if self._build_properties['key_ra'] is None:
            self._automatic_key_match_("ra")
            
        if self._build_properties['key_dec'] is None:
            self._automatic_key_match_("dec")
    
        self._update_contours_()

    def set_starsid(self,key, value, testkey=True):
        """ provide the information on how to identify a star.
        A star will then be any entry for which its `key` equals to `value`
        Set testkey to False to avoid checking if the given key exist in the catalogue
        """
        if self.has_data() and key is not None and key not in self.data.keys() and testkey:
            raise ValueError("%s is not a known data entry. Set testkey to avoid the test"%key)
        if key is not None and value is None:
            warnings.warn("No star value set (star_value = None)")
            
        self._build_properties['key_class'] = key
        self._build_properties['value_star'] = value
        
    def extract(self,contours):
        """  get a subpart of the existing catalogue based on the given 'contours'.
        'contours' is a shapely.Polygon or a matplotlib.patches.Polygon
        (see shape.point_in_contours)
        """
        mask = shape.point_in_contours(self._ra,self._dec,contours) # WRONG
        copy_ = self.copy()
        copy_.create(self.data[mask],None,force_it=True)
        return copy_

    
    def join(self,datatable):
        """ add data to the current catalogue.
        This is based on astropy.Table join:
          "
            The join() method allows one to merge these two tables into a single table
            based on matching values in the “key columns”.
          "
        We use the join_type='outer'
        (http://docs.astropy.org/en/stable/table/operations.html)
        """
        # ---------------------
        # - Input Test
        if type(datatable) is not Table:
            try:
                datatable = Table(datatable)
            except:
                raise TypeError("the given datatable is not an astropy. Table and cannot be converted into.")

        from astropy.table import join
        
        self._properties["data"] = join(self.data,datatable,join_type='outer')
        self._update_fovmask_()

    def merge(self,catalogue_):
        """ Combing a given catalogue to the current instance.
        This makes use of the join() method.
        """
        if "__nature__" not in dir(catalogue_) or catalogue_.__nature__ != "Catalogue":
            raise TypeError("the input 'catalogue' must be an astrobject catalogue")

        self.join(catalogue_.data)
        self._derived_properties["contours"] = self.contours.union(catalogue_.contours)
        
            
    def writeto(self,savefile,format="ascii",force_it=True,
                fill_values=[(ascii.masked, "nan"), ("--","nan"),("","nan")],
                **kwargs):
        """ save the catalogue as pkl or fits files.
        The fits wil be used if this has a header, the pkl otherwise
        """
        # -- First file
        if not self.header is None and len(self.header.keys())>0 and format=="fits":
            self._writeto_fits_(savefile,force_it=force_it,**kwargs)
        # -- pkl file
        elif savefile.endswith(".pkl") and format=="pkl":
            self._writeto_pkl_(savefile,force_it=force_it,**kwargs)
        else:
            self.data.write(savefile,format=format,**kwargs)

    # --------------------- #
    # Unset Methods         #
    # --------------------- #
    def reset(self):
        """This method removed masking on the catalogue, i.e. no more FoV defined, no more matching mask"""
        # ---------------------
        # - No more matching
        self.set_matchedmask(None)
        # ---------------------
        # - No more fovmask 
        self.remove_fovmask()

    def remove_fovmask(self):
        """
        """
        self._load_default_fovmask_()
        self._side_properties["fovcontours"] = None
        
    # --------------------- #
    # Set Methods           #
    # --------------------- #
    def set_mag_keys(self,key_mag,key_magerr):
        """ provide the catalogue entry  associated with magnitude """
        self._build_properties["key_mag"] = key_mag
        self._build_properties["key_magerr"] = key_magerr

    def set_coord_keys(self,key_ra,key_dec):
        """ provide the catalogue entry  associated with coordinates (Ra and Dec) """
        self._build_properties["key_ra"] = key_ra
        self._build_properties["key_dec"] = key_dec
        
    def set_wcs(self,wcs,force_it=False,update_fovmask=True):
        """
        """
        super(Catalogue, self).set_wcs(wcs, force_it=force_it)
        
        if update_fovmask:
            if self.has_wcs() and shape.HAS_SHAPELY and self.wcs.has_contours():
                self.set_fovmask(wcs=self.wcs,update=False)
            else:
                warnings.warn("loading default fovmask since no wcs solution or no Shapely or wcssolution without image")
                self._load_default_fovmask_()

    def set_fovmask(self, wcs=None, fovcontours=None,
                    update=True):
        """
        This methods enable to define the mask of catalgue objects within the
        given field of view.
        
        Parameters
        ----------
        fovcontours: [shapely polygon]
        
        
        update: [bool] -optional-
            True to have a consistent object. Set False
            only if you know what you are doing
            
        Returns
        -------
        Void
        """
        if wcs is not None:
            fovcontours = wcs.contours
            
        elif fovcontours is None:
            raise ValueError("Either wcs or fovcontours must be provided")
        
        self.fovmask = self.get_contour_mask(fovcontours, infov=False)
        self._side_properties["fovcontours"] = fovcontours
        
    def set_matchedmask(self,matchedmask):
        """
        This methods enable to set to matchedmask, this mask is an addon
        mask that indicate which point from the catalogue (after the fov cut)
        has been matched by for instance a sextractor/sep extraction.
        
        Set None to remove the matching association
        """
        if matchedmask is None:
            self._side_properties["matchedmask"] = None
            return
        if len(matchedmask) == 0:
            self._side_properties["matchedmask"] = np.zeros(self.nobjects_in_fov)
            return
        if type(matchedmask[0]) is bool:
            # - it already is a mask, good
            self._side_properties["matchedmask"] = np.asarray(matchedmask,dtype=bool)
            return

        if type(matchedmask[0]) in [int,np.int32,np.int64]:
            # it must be a list of matched index (from SkyCoord matching fuction e.g.)
            self._side_properties["matchedmask"] = \
              np.asarray([i in matchedmask for i in range(len(self.ra))],dtype=bool)
            return
        
        raise TypeError("the format of the given 'matchedmask' is not recongnized. "+\
                "You could give a booleen mask array or a list of accepted index")
                     
    def set_ingalaxymask(self, galaxycontours, reset=False):
        """
        The will update the current ingalaxy mask with the given contours.
        If a galaxy mask already exist, this will only check the None value,
        use reset=True
        to restart the process from stratch (slower).

        This will update the fundamental self.data. This way, saving it will save
        this information.
        """
        
        # -- Sarting points
        galmask = np.asarray([np.NaN]*self.nobjects) \
          if not self.has_ingalaxymask() or reset else self._ingalaxymask
        
        # -- ID to work with, i.e. are not None and are in the FoV
        idx = np.asarray([i for i in np.arange(self.nobjects)
               if galmask[i]!=galmask[i] and self.fovmask[i]])
        
        if len(idx) == 0:
            warnings.warn("No new coordinates needs a 'ingalaxy' to check")
            return
         
        # -- ID that are galaxies
        idxgal    = np.asarray(idx[~self._starmask[idx]])
        idxnotgal = np.asarray(idx[ self._starmask[idx]])
        
        # -- ID that are not galaxies but are in
        from ...utils.shape import Point
        listingal = idxnotgal[\
            np.asarray([galaxycontours.contains(p_) for p_ in \
            [Point(x_,y_) for x_,y_ in \
             self.wcs.world2pix(self._ra[idxnotgal],self._dec[idxnotgal])]
             ],dtype=bool)]
        
        # -- And let's set what should be
        gal = idxgal.tolist() + listingal.tolist()
        galmask[idx] = False # Not in galaxy
        galmask[np.asarray(gal)] = True# Except if they are

        self.data.add_column(Column(galmask,name="ingalaxy"),
                             rename_duplicate="ingalaxy" in self.data.keys())
        
    # --------------------- #
    #  convertion methods   #
    # --------------------- #
    def id_to_idx(self, id_, mask=None, infov=True):
        """ uses get_value_idx to return the idx of the given id.

        development note: id entry depends on the catalogue source.

        Returns
        -------
        [int]
        """
        if "key_id" not in self._build_properties.keys() or self._build_properties["key_id"] is None:
            raise AttributeError("no 'key_id' set in the _build_properties for this catalogue")
        
        return self.get_value_idx(self._build_properties["key_id"], id_,
                                  mask=mask, infov=infov)

    def idx_to_id(self, idx, **kwargs):
        """ """
        if "key_id" not in self._build_properties.keys() or self._build_properties["key_id"] is None:
            raise AttributeError("no 'key_id' set in the _build_properties for this catalogue")
        if not is_arraylike(idx):
            idx = [idx]
            
        return self.get(self._build_properties["key_id"], idx, **kwargs)

    # --------------------- #
    # Get Methods           #
    # --------------------- #
    def get(self, key, mask=None, infov=True):
        """
        get any 'key' known by the instance (self.`key`) or more generally
        in the data. You can mask the returned data (*CAUTION* some values
        have default fovmasking (ra,dec... use _`key` like _ra to have the
        none fov mask values.
        
        *Remark* 'key' could be an list of keys.

        Returns:
        --------
        array (or list of)
        """
        if is_arraylike(key):
            return [self.get(key_,mask=mask) for key_ in key]
        
        if key in dir(self):
            val_ = eval("self.%s"%key) if infov else eval("self._%s"%key)
        elif key in self.data.keys():
            val_ = self.data[key][self.fovmask] if infov else self.data[key]
        else:
            raise ValueError("Unknown key %s"%key)
        
        return val_ if mask is None else val_[mask]
    
    def get_value_idx(self, key, value, mask=None, infov=True):
        """ get the index of the entry for which `key`'s value is `value`
        If you provide a mask, this will be the index of the mask entry that
        are True.

        Parameters
        ----------
        key: [string]
            Entry to search
        value: [any]
            value that will be tested.
        Returns
        -------
        int / None
        """
        values = self.get(key, mask=mask, infov=infov)
        id_ = np.where(values==value)
        
        return None if len(id_)==0 else id_[0] if len(id_)==1 else id_
    
    def get_subcatalogue(self, contours=None, catmag_range=[None,None],
                         stars_only=False, isolated_only=False,
                         fovmask=False,
                         **kwargs):
        """ Returns a value of a sub fov of the catalogue. Only objects within the contours' fov
        will be contained in the returned catalogue """

        mask = self.get_mask(contours=contours, catmag_range=catmag_range,
                            stars_only=stars_only, isolated_only=isolated_only,
                            fovmask=False)
        if fovmask:
            mask = mask * self.fovmask
        
        subcat = self.__class__(empty=True)
        subcat.create(self.data[mask], None, force_it=True,**self._build_properties)
        return subcat

    def get_mask(self,catmag_range=[None,None],stars_only=False,
                 isolated_only=False, nonstars_only=False,
                 contours=None, notingalaxy=False, matched=False,
                 fovmask=True):
        """ This returns a bolean mask following the argument cuts. """
        
        mask = np.ones(self.nobjects_in_fov, dtype="bool") if fovmask else\
          np.ones(self.nobjects, dtype="bool")

        # - stars etc.        
        if stars_only:
            mask *= self.starmask if fovmask else self._starmask
        if nonstars_only:
            if stars_only:
                warnings.warn("WARNING you ask for both stars_only and nonstars_only !!!")
            mask *= ~self.starmask if fovmask else ~self._starmask
        # - isolation
        if isolated_only:
            mask *= self.isolatedmask if fovmask else self._isolatedmask

        if matched:
            if not self.has_matchedmask():
                raise AttributeError("No matching set. See set_matchedmask()")
            if not fovmask:
                raise AttributeError("Matching association only made in combination with fovmask=True")
            
            mask *= self.matchedmask
            
        # - not in galaxy
        if notingalaxy:
            if not self.has_ingalaxymask():
                raise AttributeError("No 'ingalaxymask' set. See set_ingalaxymask()")
            if not fovmask:
                if [None] in self._ingalaxymask:
                    raise ValueError("Some of the requested index have no 'ingalaxy'"+\
                    " information (None value). Please run set_ingalaxymask for them"+\
                    " (Info - the not in fovmask option might be a problem )")
                mask *= ~np.asarray(self._ingalaxymask,dtype="bool")
            else:
                if [None] in self.ingalaxymask:
                    raise ValueError("Some of the requested index have no 'ingalaxy'"+\
                    " information (None value). Please run set_ingalaxymask()")
                mask *= ~np.asarray(self.ingalaxymask,dtype="bool")
                
        # - contours
        if contours is not None:
            mask *= self.get_contour_mask(contours,infov=fovmask)
        # - magcut
        if catmag_range[0] is not None or catmag_range[1] is not None:
            if catmag_range[0] is None:
                catmag_range[0] = np.nanmin(self.mag if fovmask else self._mag)
            if catmag_range[1] is None:
                catmag_range[1] = np.nanmax(self.mag if fovmask else self._mag)
                
            magmask = (self.mag>=catmag_range[0]) & (self.mag<=catmag_range[1]) \
              if fovmask \
              else (self._mag>=catmag_range[0]) & (self._mag<=catmag_range[1])
            mask *= magmask
         
        return mask

    def get_mask_around(self,ra,dec,radius,runits="arcsec",wcs_coords=True,
                        stars_only=False, isolated_only=False, infov=True,**kwargs):
        """
        return the mask of the ra dec 
        -- kwargs goes to get_mask() --
        """
        maskbase = self.get_mask(stars_only=stars_only, isolated_only=isolated_only,
                                 fovmask=infov,**kwargs)
        maskaround = self.idx_to_mask(self.get_idx_around(ra,dec,radius, runits=runits,
                                                          wcs_coords=wcs_coords,infov=infov)[0],
                                                          infov=infov)
        return maskaround * maskbase
        

    def get_nearest_idx(self, ra, dec, wcs_coords=True, mask=None, infov=True):
        """ get the index of the nearest (masked-)catalogue entry
        
        Parameters
        ----------
        ra, dec : [float (or array-of), float (or array-of)]
            Coordinates that should be matched

        wcs_coords: [bool] -optional-
            True if the ra and dec are given in degree. Set to False
            if you provided the pixel coordinates.

        mask: [bool-array] -optional-
            Mask the catalogue to only given the nearest entry after of the given
            selection.
            ===
            **CAUTION** The idx will then be that of the **mask-catalogue**
            ===
            (use np.argwhere(mask)[get_nearest_idx(..)[0]] to get the index
            of the unmasked catalogue)
            
        Returns
        -------
        idx, sep2d
        """
        # --------------
        # - Input 
        if not wcs_coords and not self.has_wcs():
            raise AttributeError("Needs a wcs solution to get pixel coordinates")
        if not wcs_coords:
            ra,dec = np.asarray(self.wcs.pix2world(ra,dec)).T
        if not is_arraylike(ra):
            ra = [ra]
            dec = [dec]
        skytarget = coordinates.SkyCoord(ra*units.degree,dec*units.degree)
        
        # -------------
        # - Cat matching
        catsky = self.sky_radec if infov else self._sky_radec
        if mask is not None:
            catsky = catsky[mask]
            
        return coordinates.match_coordinates_sky(skytarget, catsky)[:2]
        
    def get_idx_around(self,ra,dec,radius,runits="arcsec",wcs_coords=True,
                       infov=True):
        """
        Returns the catalogue indexes of the elements within `radius` `runits`around
        the `ra` `dec` location.
        Returns:
        --------
        2xN index array (idx, angular sep. N is the number of matching elements.)
        """
        # --------------
        # - Input 
        if not wcs_coords and not self.has_wcs():
            raise AttributeError("Needs a wcs solution to get pixel coordinates")
        if not wcs_coords:
            ra,dec = np.asarray(self.wcs.pix2world(ra,dec)).T
        if not is_arraylike(ra):
            ra = [ra]
            dec = [dec]
            
        skytarget = coordinates.SkyCoord(ra*units.degree,dec*units.degree)
        catsky = self.sky_radec if infov else self._sky_radec
        return catsky.search_around_sky(skytarget,
                                        radius*units.Unit(runits))[1:3]
        
    def get_contour_mask(self, contours, infov=True):
        """  returns a boolean array for the given contours """
        if contours is None:
            return np.asarray([True for ra in self._ra])
        
        if type(contours) != shape.polygon.Polygon and\
           type(contours) != shape.multipolygon.MultiPolygon:
            raise TypeError("contours must be a shapely Polygon or MultiPolygon")
        _ra = self.ra if infov else self._ra
        _dec = self.dec if infov else self._dec
        return np.asarray([shape.point_in_contours(ra,dec, contours)
                           for ra,dec in zip(_ra,_dec)])

    # --------------------- #
    # Convertors            #
    # --------------------- #

    
    # --------------------- #
    # Convertors            #
    # --------------------- #
    def idx_to_mask(self,idx, infov=False):
        mask = np.zeros(self.nobjects,dtype=bool) if not infov else np.zeros(self.nobjects_in_fov,dtype=bool)
        for i in idx:
            mask[i] = True
        return mask

    # --------------------- #
    # Exclusion             #
    # --------------------- #
    def exclude_source(self,key,value):
        """ Exclude the cases where data[key] == value. Then use the exclusionmask """
        if key not in self.data.keys():
            raise ValueError("the given key (%s) is not known by self.data"%key)
        ids_to_exclude = np.argwhere(self.data[key]==value)
        if len(ids_to_exclude)==0:
            warnings.warn("WARNING: No value excluded. not match")
            return
        self._side_properties["excluded_list"] = self.excluded_list.tolist()+ids_to_exclude.tolist()
        
    def clear_excluded_list(self):
        """ empty the exclusion list """
        self._side_properties["excluded_list"] = None
        
    # --------------------- #
    # PLOT METHODS          #
    # --------------------- #
    def display(self,ax,wcs_coords=True,draw=True,
                apply_machedmask=True,draw_contours=True,
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
            print("Catalogue has no 'data' to display")
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
        for x_,y_,show_,propextra in [[x_starin,  y_starin,True,{}],
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

        if self.contours is not None and wcs_coords and draw_contours:
            shape.draw_polygon(ax,self.contours,ec=None)
            
        if draw:
            ax.figure.canvas.draw()
            
        return axout
    
    # =========================== #
    # Internal Methods            #
    # =========================== #
    # ------------------
    # --  Save files
    def _writeto_pkl_(self,savefile,force_it=True):
        """This the current catalogue has a pkl file"""
        if not savefile.endswith(".pkl"):
            savefile +=".pkl"
            
        dump_pkl(self.data,savefile)
    
    def _writeto_fits_(self,savefile,force_it=True):
        raise NotImplementedError("to be done")
    # ------------------
    # --  Update
    def _update_fovmask_(self):
        """
        """
        if self.has_wcs():
            self.set_fovmask(wcs=self.wcs)
        return
    
    def _update_contours_(self):
        """
        """
        if not shape.HAS_SHAPELY:
            self._derived_properties["contours"] = None
        else:
            # -- This roughly take 0.2s for 1e4 objects
            self._derived_properties["contours"] = \
                shape.get_contour_polygon(np.asarray(self._ra),
                                        np.asarray(self._dec))
        
    # ------------------
    # --  Key match
    def _automatic_key_match_(self, key,build_key=None):
        """
        """
        try:
            knownkeys = list(self.data.keys())
        except:
            warnings.warn("WARNING no automatic key available (data.keys() failed)")
            return
        vkey = [k for k in knownkeys if key in k.lower() if not k.startswith("e_")]
        if len(vkey) >1:
            warnings.warn("WARNING ambiguous %s key. Use set_coord_keys. "%key+", ".join(vkey))
            return
        if len(vkey) ==0:
            if key.lower() == "dec":
                return self._automatic_key_match_("de",build_key="dec")
            print("WARNING no match found for the %s key. Use set_coord_keys. "%key)
            print("       known keys: "+", ".join(knownkeys))
            return

        if build_key is None:
            build_key = key
        self._build_properties['key_%s'%build_key] = vkey[0]


    def _automatic_key_ra_match_(self):
        """
        """
        self._automatic_key_match_("ra")

    def _automatic_key_dec_match_(self):
        """
        """
        self._automatic_key_match_("dec")

    
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
    def nobjects(self):
        if self.data is None:
            return 0
        
        if self.header is not None and "NAXIS2" in self.header:
            return self.header["NAXIS2"]

        if type(self.data) is dict:
            return len(self.data.values()[0])
        
        return len(self.data)

    @property
    def nobjects_in_fov(self):
        return len(self.ra)
    
    # ---------------
    # - header / wcs 
    @property
    def header(self):
        return self._properties["header"]
    
    @property
    def fovmask(self):
        if self._side_properties["fovmask"] is None:
            self._load_default_fovmask_()
        return self._side_properties["fovmask"]
    
    def _load_default_fovmask_(self):
        self._side_properties["fovmask"] = np.ones(self.nobjects,dtype=bool)
        
    @fovmask.setter
    def fovmask(self,newmask):
        if len(newmask) != self.nobjects:
            raise ValueError("the given 'mask' must have the size of 'ra'")
        self._side_properties["fovmask"] = newmask

    # -- Exclusion
    @property
    def excluded_list(self):
        if self._side_properties["excluded_list"] is None:
            self._side_properties["excluded_list"] = []
        return np.asarray(self._side_properties["excluded_list"])
    
    def has_excluded_cases(self):
        return len(self.excluded_list)>0

    @property
    def excludedmask(self):
        return self._excludedmask[self.fovmask]
    
    @property
    def _excludedmask(self):
        return self.idx_to_mask(self.excluded_list,infov=False)
        
    # -- Matching
    @property
    def matchedmask(self):
        return self._side_properties["matchedmask"]

    def has_matchedmask(self):
        return self.matchedmask is not None

    # -- InGalaxy Mask
    @property
    def ingalaxymask(self):
        return self._ingalaxymask[self.fovmask] \
          if self.has_ingalaxymask() else None
    
    @property
    def _ingalaxymask(self):
        return None if "ingalaxy" not in self.data.keys() else\
          self.data["ingalaxy"]
    
    def has_ingalaxymask(self):
        return self._ingalaxymask is not None
    
    # ------------
    # - on flight
    # - coords
    @property
    def ra(self):
        """Barycenter position along world x axis"""
        return self._ra[self.fovmask]
    
    @property
    def _ra(self):
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["key_ra"]]
    
    @property
    def dec(self):
        """arycenter position along world y axis"""
        return self._dec[self.fovmask]
    
    @property
    def _dec(self):
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        return self.data[self._build_properties["key_dec"]]
        
    @property
    def sky_radec(self):
        """This is an advanced radec methods tight to astropy SkyCoords"""
        return coordinates.SkyCoord(ra=self.ra,dec=self.dec, unit="deg")

    @property
    def _sky_radec(self):
        """This is an advanced radec methods tight to astropy SkyCoords"""
        return coordinates.SkyCoord(ra=self._ra,dec=self._dec, unit="deg")
    
    # - mag
    @property
    def mag(self):
        """Generic magnitude"""
        return self._mag[self.fovmask]

    @property
    def _mag(self):
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        if not self._is_keymag_set_():
            raise AttributeError("no 'key_mag' defined. see self.set_mag_keys ")
        
        return self.data[self._build_properties["key_mag"]]
        
    @property
    def mag_err(self):
        """Generic magnitude RMS error"""
        return self._mag_err[self.fovmask]

    @property
    def _mag_err(self):
        """Generic magnitude RMS error"""
        if not self.has_data():
            raise AttributeError("no 'data' loaded")
        if not self._is_keymag_set_():
            raise AttributeError("no 'key_magerr' defined. see self.set_mag_keys ")
        
        return self.data[self._build_properties["key_magerr"]]
    
    # ----------------
    # - Fluxes
    @property
    def lbda(self):
        return self._side_properties["lbda"]

    @lbda.setter
    def lbda(self,value):
        self._side_properties["lbda"] = value
        
    @property
    def _flux_fluxerr(self):
        if self.lbda is None:
            raise AttributeError("set 'lbda' first (self.lbda = ...)")
        
        return mag_to_flux(self.mag,self.mag_err,self.lbda)
    
    @property
    def flux(self):
        return self._flux_fluxerr[0]
    
    @property
    def flux_err(self):
        return self._flux_fluxerr[1]
    
    def _is_keymag_set_(self,verbose=True):
        """this method test if the keymag has been set"""
        if self._build_properties["key_mag"] is None or \
          self._build_properties["key_magerr"] is None:
            if verbose:
                print("No 'key_mag'/'key_magerr' set ; call 'set_mag_keys'." +\
                  " List of potential keys: "\
                  +", ".join([k for k in self.data.keys() if "mag" in k or "MAG" in k]))
            return False
        return True
    
    @property
    def objecttype(self):
        return self._objecttype[self.fovmask]

    @property
    def _objecttype(self):
        if "key_class" not in self._build_properties.keys():
            raise AttributeError("no 'key_class' provided in the _build_properties.")
        
        return np.asarray(self.data[self._build_properties["key_class"]])
    
    @property
    def starmask(self, infov=True):
        """ This will tell which of the datapoints is a star
        Remark, you need to have defined key_class and value_star
        in the __build_properties to be able to have access to this mask
        """
        if "value_star" not in self._build_properties.keys() or \
          self._build_properties["value_star"] is None:
            return None
        return (self.objecttype == self._build_properties["value_star"])

    @property
    def _starmask(self):
        """ This will tell which of the datapoints is a star
        Remark, you need to have defined key_class and value_star
        in the __build_properties to be able to have access to this mask
        """
        if "value_star" not in self._build_properties.keys():
            return None
        return (self._objecttype == self._build_properties["value_star"])
        
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
            return np.asarray(self.wcs.world2pix(self.ra,self.dec)).T
        raise AttributeError("no 'wcs' solution loaded")

    # ----------------------
    # - Alone Object
    def define_around(self,ang_distance):
        """
        """
        # -- FoV cut
        idxcatalogue = self.sky_radec.search_around_sky(self.sky_radec,
                                                        ang_distance)[0]
        self._derived_properties["naround"] = np.bincount(idxcatalogue)
        # -- no FoV cut
        _idxcatalogue = self._sky_radec.search_around_sky(self._sky_radec,
                                                        ang_distance)[0]
        self._derived_properties["naround_nofovcut"] = np.bincount(_idxcatalogue)
        
    def _is_around_defined(self):
        return self._derived_properties["naround"] is not None
        
    @property
    def nobjects_around(self):
        """ Number of object around an object. """
        if not self._is_around_defined():
            print("INFORMATION: run 'define_around' to set nobject_around")
        return self._derived_properties["naround"]

    @property
    def _nobjects_around(self):
        """ Number of object around an object. No FoV Cut"""
        if not self._is_around_defined():
            print("INFORMATION: run 'define_around' to set nobject_around")
        return self._derived_properties["naround_nofovcut"]
    
    @property
    def isolatedmask(self):
        """
        """
        if not self._is_around_defined():
            raise AttributeError("no 'nobjects_around' parameter derived."+\
                                 " Run 'define_around'")
        return (self.nobjects_around == 1)

    @property
    def _isolatedmask(self):
        """
        """
        if not self._is_around_defined():
            raise AttributeError("no 'nobjects_around' parameter derived."+\
                                 " Run 'define_around'")
        return (self._nobjects_around == 1)
    

    # ----------------------
    # - Shapely
    @property
    def contours(self):
        return self._derived_properties["contours"]

    @property
    def contours_pxl(self,**kwargs):
        """Based on contours (in wcs) and wcs2pxl, this draws the pixels contours"""
        if not self.has_wcs():
            raise AttributeError("no wcs solution loaded. You need one.")
        x,y = np.asarray([self.wcs.world2pix(ra_,dec_) for ra_,dec_ in
                          np.asarray(self.contours.exterior.xy).T]).T
        # switch ra and dec ;  checked
        return shape.get_contour_polygon(x,y)
    
    @property
    def fovcontours(self):
        return self._side_properties["fovcontours"]

    @property
    def fovcontours_pxl(self,**kwargs):
        """Based on contours (in wcs) and wcs2pxl, this draws the pixels contours"""
        if not self.has_wcs():
            raise AttributeError("no wcs solution loaded. You need one.")
        x,y = np.asarray([self.wcs.world2pix(ra_,dec_) for ra_,dec_ in
                          np.asarray(self.fovcontours.exterior.xy).T]).T
        # switch ra and dec ;  checked
        return shape.get_contour_polygon(x,y)

    
############################################
#                                          #
# Backup Plan if no SNCOSMO                #
#                                          #
############################################
    
class _DefaultBandpass_( BaseObject ):
    """ """
    PROPERTIES = ["bandname", "wave_eff"]
    def __init__(self, bandname, wave_eff):
        """ """
        self._properties["bandname"] = bandname
        self._properties["wave_eff"] = wave_eff

    @property
    def name(self):
        return self._properties["bandname"]
    @property
    def wave_eff(self):
        return self._properties["wave_eff"]
