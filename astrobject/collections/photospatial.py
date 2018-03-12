#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np

from astropy            import coordinates, units, table
from ..baseobject       import WCSHandler, CatalogueHandler
from ..photometry       import get_photopoint
from ..collection       import BaseCollection, PhotoPointCollection
from ..utils.tools      import kwargs_update, dump_pkl, load_pkl, is_arraylike
from ..utils.decorators import _autogen_docstring_inheritance


__all__ = ["get_photomap","get_sepobject"]


def get_photomap(photopoints=None,coords=None,wcs_coords=True, **kwargs):
    """ Create a collection of photopoints

    Parameters
    ----------
    photopoints: [array of photopoint]
        Array of the photopoint you want to unify within a unique Collection
        (a PhotoMap)

    coords: [2D-array]
        List of coordinates [[x0,y0],[x1,y1] etc.] associated to the given
        photopoint array.

    wcs_coords: [bool]
        True if the coordinate have been given in Ra,Dec, False if they are
        given in pixel coordinates.

    **kwargs goes to PhotoMap init

    Returns
    -------
    PhotoMap (Collection of PhotoPoints)
    """
    return PhotoMap(photopoints=photopoints,
                    coords=coords,
                    wcs_coords=wcs_coords,
                    **kwargs)

def get_sepobject(sepoutput, ppointkwargs={},
                  use_peakposition=False, wcs_extension=None,
                   **kwargs):
    """ Create a Spectial collection of photopoints dedicated
    to the SEP's output.

    Parameters
    ----------
    sepoutput: [array]
        Whatever is returned by sep.extract()

    ppointkwargs: [dict]
        Potential meta value associated to the individual photopoints
        of the PhotoMap. Should have the format {'k':v} where v is an array
        of the same size as the number of photopoint or being a unique
        value.
        See set_meta() of PhotoPoint

    wcs_extension: [int/None] -optional-
        If the data (a fits file) has a wcs recorded. Provide here the extension.
        If None, no wcs solution will be loaded.
        
    **kwargs goes to SepObject init (Child of PhotoMap)

    Returns
    -------
    SepObject (Child of PhotoMap)
    """
    xkey,ykey = ["x","y"] if not use_peakposition else ["xpeak","ypeak"]
    
    # - Test aninput parsing made at this level
    inputphotomap = parse_sepoutput(sepoutput,**ppointkwargs) if is_sepoutput(sepoutput) else \
      parse_sepsex_fitsdata(sepoutput)

    # - good to go.
    pmap = SepObject(photopoints=inputphotomap["ppoints"],
                     coords=inputphotomap["coords"],
                     wcs_coords=False, **kwargs)

    ids_sorting = [pmap.coords_to_id(x,y)
                       for x,y in zip(inputphotomap["meta"][xkey].data,
                                      inputphotomap["meta"][ykey].data)]
        
    [pmap.set_meta(k,inputphotomap["meta"][k].data, ids=ids_sorting)
         for k in inputphotomap["meta"].keys()]
    
    # WCS Solution
    # --------------
    if wcs_extension is not None:
        from ..astrometry import wcs
        try:
            wcssolution = wcs(sepoutput, extension=wcs_extension)
        except:
            warnings.warn("Failing loading the wcs solution")
            wcssolution = None
            
        if wcssolution is not None:
            pmap.set_wcs(wcssolution)
        
    return pmap



def get_photomapcollection( photomaps, wcs_extension=None, **kwargs ):
    """ Get a PhotoMap collection

    Parameters
    ----------
    photomaps: [list of photomaps or filenames of]
        The photomaps of roughly a same area that you want to combine into a collection
        if a list of filename is given, get_sepobject() will be used to read them.
        This could be e.g. sep or sextractor outputs

    wcs_extension: [int/None] -optional-
        If the given data are fitfiles, you can specify the extension of the header
        containing the wcs solution if any. Leave this to None if no wcs solution have
        to be loaded.

    **kwargs goes to the PhotoMapCollection __init__ (catalogue etc.)
    
    Returns
    -------
    PhotoMapCollection
    """
    if not is_arraylike(photomaps):
        raise TypeError("photomaps must be a list/array")

    # -- loading data
    if type(photomaps[0]) == str:
        photomaps = [get_sepobject(filename, wcs_extension=wcs_extension, **kwargs) for filename in photomaps]
        if wcs_extension is not None:
            [p._derive_radec_() for p in photomaps if p.has_wcs()]

    return PhotoMapCollection(photomaps, **kwargs)
# ======================= #
#                         #
# Internal Functions      #
#                         #
# ======================= #
# ------------ #
#  From SEP    #
# ------------ #
def is_sepoutput(array):
    """
    """

    if type(array) != np.ndarray:
        return False
    tsep = table.Table(array)
    return "cxx" in tsep.keys()

def parse_sepoutput(sepoutput,lbda=None,**kwargs):
    """ convert entries from sep or formatted sextractor output.

    Parameters
    ----------
    sepoutput: [ndarray or dict]
        output of sep.extract or sextractor output formatted by `parse_sepsex_fitsdata()`

    lbda: [float/None] -optional-
        provide the wavelength [in angstrom] of the photometric measurements.
        
    Returns
    -------
    dictionary (formatted to by the input of SepObject)

    
    """
    if table.table.Table not in sepoutput.__class__.__mro__:
        try:
            sepoutput = table.Table(sepoutput)
        except:
            raise TypeError("the given 'sexoutput' cannot be converted into astropy's Table")

    ppoints = [get_photopoint(lbda=lbda,flux=t_["flux"],var=None,
                          source="sepextract",**kwargs)
                          for t_ in sepoutput]
    coords = np.asarray([sepoutput["x"],sepoutput["y"]]).T
    return {"ppoints":ppoints,"coords":coords,"meta":sepoutput[[t_ for t_ in sepoutput.keys() if t_ not in ["flux"]]]}
    
# ------------------- #
#  From Sextractor    #
# ------------------- #
def parse_sepsex_fitsdata( filename, lbda=None, fluxkey="flux_auto", **kwargs):
    """ This checks is the fits file is from sextractor """
    from astropy.io import fits
    # ----------
    # - Input
    if type(filename) == str:
        try:
            header = fits.getheader(filename)
            data   = fits.getdata(filename)
        except:
            raise TypeError("cannot get the data and header from the given file %s. Not a fits file?"%filename)
    else:
        data = filename
    
    # ----------
    #  formating
    if type(data) in [fits.fitsrec.FITS_rec]:
        data = {col.name.lower():data[col.name] for col in data.columns}
    elif hasattr(data,"columns"):
        data = {col.lower():data[col] for col in data.columns}

    # ----------
    #  formating    
    if "x_world" in data.keys():
        #  SEXTRACTOR
        data["ra"]        = data["x_world"]
        data["dec"]       = data["y_world"]
        data["x"]         = data["x_image"]
        data["y"]         = data["y_image"]
        data["a"]         = data["a_image"]
        data["b"]         = data["b_image"]
        data["a.err"]     = data["erra_image"]
        data["b.err"]     = data["errb_image"]
        data["theta"]     = data["theta_image"]
        
        data["theta.err"] = data["errtheta_image"]
        
        if fluxkey.lower() not in data.keys():
            raise ValueError("Unknwon key entry (%s)"%fluxkey +" known flux entries: "+", ".join([k for k in data.keys() if "flux" in k]))
        data["flux"] = data[fluxkey.lower()]
        data["flux.err"] = data[fluxkey.lower().replace("_","err_")] if fluxkey.lower().replace("_","err_") in data.keys() else None

    return parse_sepoutput(data,lbda, **kwargs)

######################################
#                                    #
# PhotoMetric Mapping                #
#                                    #
######################################

class PhotoMap( PhotoPointCollection, WCSHandler, CatalogueHandler ):
    """
    """

    PROPERTIES         = []
    SIDE_PROPERTIES    = ["refmap"]
    DERIVED_PROPERTIES = ["catmatch"]

    
    def __init__(self, photopoints=None, coords=None,
                 filein=None, empty=False, wcs=None,
                 catalogue=None,
                 **kwargs):
        """

        Parameters
        ----------
        photopoints: [list/array of photopoints] -optional-
            List of astrobject's PhotoPoints that will set this Collection.

        coords: [list of coordinates] -optional-
            List of Nx2 array [[x,y],[x,y]...] location of the PhotoPoints
            
        filein: [string] -optional-
            File containing PhotoMap data.

        wcs: [astrobject's WCS] -optional-
            Attach a wcs solution to the instance.
            (ignored if photopoints are not given)
            
        catalogue: [astrobject's Catalogue] -optional-
            Attach a catalogue to the instance.
            (ignored if photopoints are not given)
            
        empty: [bool] -optional-
            return an empty instance of this object

        **kwargs goes to load / create
        
        Returns
        -------
        Void (defines the object)
        """
        self.__build__()
        if empty:
            return
        if filein is not None:
            self.load(filein,**kwargs)
            
        if photopoints is not None:
            self.create(photopoints,coords,wcs=wcs,catalogue=catalogue,**kwargs)

    # Super Loading including catmatch information
    def _load_pkl_(self,filename=None, **kwargs):
        # Super
        datain = super(PhotoMap, self)._load_pkl_(filename,get_loaded_data=True,**kwargs)
        # Do you also have a catalog?
        if "catmatch" in datain.keys():
            if "idx_catalogue" in datain["catmatch"].keys() and \
               "idx" in datain["catmatch"].keys():
                self._derived_properties["catmatch"] = datain["catmatch"]
        
    def _read_table_comment_(self, comments):
        """
        """
        for c in comments:
            key,value = c.split()
            if key == "wcsid":
                self._build_properties["wcsid"]=bool(value)
            else:
                print("comment not saved: ", c)

    def create_from_table(self, table, idkey=None):
        """ """
        super(PhotoMap, self).create_from_table(table, idkey="id")
        if "x" in self._table.keys() and "y" in self._table.keys():
            self._build_properties["wcsid"] = False
            
        # - Set the ID for list_id
        if idkey is None or idkey not in table.colnames:
            idkey = "id"
            if idkey in table.colnames:
                idkey +="1"
            if not self._wcsid:
                table.add_column(table.Column(self.coords_to_id(self.get("x"),self.get("y")),idkey))
            else:
                table.add_column(table.Column(self.coords_to_id(self.get("ra"),self.get("dec")),idkey))
                
        self._build_properties["idkey"] = idkey
            
    # =============================== #
    # = Building Methods            = #
    # =============================== #
    # ---------- #
    #  GETTER    #
    # ---------- #
    def getcat(self, param, catindex):
        """ get method based on catindex entry(ies).
        This converts the catindex into index and return the
        values for the requested parameter(s) (param). If the catindex
        has no associated index entry, the corresponding returned value
        will be nan.
        
        Returns
        -------
        array
        """
        if self.catmatch is None or len(self.catmatch) == 0:
            raise AttributeError("No catalogue matched")
        
        # - output array  
        values_out = np.ones([len(catindex), len(param)])*np.NaN if is_arraylike(param) else\
          np.ones(len(catindex))*np.NaN
        # - flag has data
        flag_hasidx = self.is_catindex_detected(catindex)
        # - value within the maskin are filled
        if np.any(flag_hasidx):
            values_out[flag_hasidx] = self.get(param, mask=self.catindex_to_index(np.asarray(catindex)[flag_hasidx]))
            
        return values_out

    def is_catindex_detected(self, catindex):
        """ check if the given catalog index has been matched """
        return np.in1d(catindex, self.catmatch["idx_catalogue"])
    
    def get_skycoords(self):
        """ """
        ra,dec = self.radec.T
        return coordinates.SkyCoord(ra=ra*units.degree,dec=dec*units.degree)
    
    # ---------------- #
    # - SUPER THEM   - #
    # ---------------- #
    def create(self, photopoints, coords, wcs_coords=True,
               wcs=None,refmaps=None,catalogue=None,**kwargs):
        """ wcs_coords means that the coordinate are given in ra,dec """
        
        super(PhotoMap,self).create(photopoints,ids=coords,**kwargs)
        self._build_properties["wcsid"] = wcs_coords
        
        if wcs is not None:
            self.set_wcs(wcs)
        if refmaps is not None:
            self.set_refmaps(wcs)
        if catalogue is not None:
            self.set_catalogue(catalogue)
        
    def add_photopoint(self,photopoint,coords,refphotopoint=None):
        """
        Add a new photopoint to the photomap.
        """
        if coords is None:
            raise ValueError("The coordinates of the photopoint are necessary")
        
        if type(coords) is str or type(coords) is np.string_:
            # -- If this failse, the coordinate format is not correct
            coords = self.id_to_coords(coords)
        elif np.shape(coords) != (2,):
            raise TypeError("the given coordinate must be a 2d array (x,y)/(ra,dec)")
        
        super(PhotoMap,self).add_photopoint(photopoint,self.coords_to_id(*coords))
        # -----------------------
        # - Complet the refmap
        if self.has_refmap():
            if refphotopoint is None:
                refphotopoint = get_photopoint(empty=True)
            self.remap.add_photopoint(refphotopoint,coords)
        
    # =============================== #
    # = Main Methods                = #
    # =============================== #
    def measure_local_photometry(self, ra, dec, radius,runits=False,
                                 catindex_exclusion=None):
        """
        Calibrate the ppoint photometry using the current photomap
        """
        # ------------------- #
        # - Get the Indexes - #
        # ------------------- #
        if not self.has_catalogue():
            raise AttributeError("No Catalogue loaded. Set one.")
        
        # -- all existing catalogue sources
        cat_indexes, angsep = self.get_idx_around(ra,dec,radius,runits=False,
                                    catalogue_idx=True)
        if len(cat_indexes)<3:
            raise ValueError("Not enough points in the given location")
        
        # -- known sources in the photomap
        catindex = cat_indexes[self.is_catindex_detected(catindex)]
        index = self.catindex_to_index(catindex)

        # ------------------- #
        # - Get Photometry  - #
        # ------------------- #        
        return catindex,index
    
    # -- Id <-> Coords
    def coords_to_id(self,x,y):
        """ this enables to have consistent coordinate/id mapping """
        if not is_arraylike(x):
            return "%.8f,%.8f"%(x,y)
        return ["%.8f,%.8f"%(x_,y_) for x_,y_ in zip(x,y)]
    
    def id_to_coords(self,id):
        """ convert photopoint's id into coordinates. Only reading here, see radec, xy """
        return np.asarray(id.split(","),dtype="float")

    # -- Catalogue <-> PhotoMap
    def match_catalogue(self,catalogue=None,
                        force_it=False, deltadist=None,
                        **kwargs):
        """
        This methods enable to attached a given sepobject entry
        to a catalog value.
        You can set a catalogue.

        **kwargs goes to the get_skycoords() method.
        
        Returns
        -------
        Void
        """
        # --------------
        # - input 
        if catalogue is not None:
            self.set_catalogue(catalogue, force_it=force_it)

        if not self.has_catalogue():
            raise AttributeError("No 'catalogue' defined or given.")
        
        if deltadist is None:
            deltadist = 3*units.arcsec
        elif type(deltadist) is not units.quantity.Quantity:
            deltadist = deltadist*units.arcsec
            
        # -- matching are made in degree space
        skyradec = self.get_skycoords(**kwargs)
        idxcatalogue, idxsepobjects, d2d, d3d = skyradec.search_around_sky(self.catalogue.sky_radec, deltadist)
        
        # --------------------
        # - Save the results
        self._derived_properties["catmatch"] = {
            "idx_catalogue":idxcatalogue,
            "idx":idxsepobjects,
            "angsep_arcsec":d2d.to("arcsec").value,
            "catquery":self._catquery
            }
        
        self.catalogue.set_matchedmask(idxcatalogue)

    def index_to_catindex(self,index, cleanindex=False):
        """ give the catalogue index corresponding to the given photomap index.
        The cleanindex option enable to automatically remove the unknown index values.
        If you do not. this will raise a ValueError indicating the unknown index
        """
        if not is_arraylike(index):
            catindex = [index]
            
        index = np.asarray(index)
        
        bool_ = np.in1d( index, self.catmatch["idx"] )
        if not np.all(bool_) and  not cleanindex:
            raise ValueError("Unknown photomap index(es):"+", ".join(index[bool_]))

        catindexmask = np.argwhere(np.in1d(self.catmatch["idx"], index[bool_]))
        # all the index[bool_] entries have a catidx associated. But some may have several
        # If so len(catindexmask) would be greater. Hence we run a slower, but accurate technique
        if len(catindexmask)>len(index[bool_]):
            warnings.warn("At least one index has several catindex matched. Closest used")
            return np.asarray([self.catmatch["idx_catalogue"][np.argwhere(self.catmatch["idx"]==i)[np.argmin(self.catmatch["angsep_arcsec"][self.catmatch["idx"]==i])]][0]
                for i in index[bool_]])
                    
        return np.concatenate(self.catmatch["idx_catalogue"][catindexmask])
                        
        
    def catindex_to_index(self,catindex, cleanindex=False):
        """ give the photomap index corresponding to the given catalogue index.
        The cleanindex option enable to automatically remove the unknown index values.
        If you do not. this will raise a ValueError indicating the unknown index
        """
        if not is_arraylike(catindex):
            catindex = [catindex]
            
        catindex = np.asarray(catindex)
        
        bool_ = np.in1d( catindex, self.catmatch["idx_catalogue"] )
        if not np.all(bool_) and not cleanindex:
            raise ValueError("Unknown catalogue index(es):"+", ".join(catindex[bool_]))

        
        indexmask =  np.argwhere(np.in1d(self.catmatch["idx_catalogue"], catindex[bool_]))
        # all the catindex[bool_] entries have a idx associated. But some may have several
        # If so len(indexmask) would be greater. Hence we run a slower, but accurate technique
        if len(indexmask)>len(catindex[bool_]):
            warnings.warn("At least one catalogue entry has several index matched. First used  ")
            return np.asarray([self.catmatch["idx"][np.argwhere(self.catmatch["idx_catalogue"]==i)[np.argmin(self.catmatch["angsep_arcsec"][self.catmatch["idx_catalogue"]==i])]][0]
                for i in catindex[bool_]])
                    
        return np.concatenate(self.catmatch["idx"][indexmask])

    # ------------- #
    # - SETTER    - #
    # ------------- #
    def set_refmap(self, photomap,force_it=False):
        """
        Set here the reference map associated to the current one.
        The ordering of the point must be consistant with that of the current
        photomap. 
        """
        if self.has_refmap() and force_it is False:
            raise AttributeError("'refmap' is already defined."+\
                    " Set force_it to True if you really known what you are doing")

        self._side_properties["refmap"] = photomap

        
    #  WCS Solution
    # -----------------
    @_autogen_docstring_inheritance(WCSHandler.set_wcs,"WCSHandler.set_wcs")
    def set_wcs(self, wcs, force_it=False,**kwargs):
        #
        #
        #
        super(PhotoMap, self).set_wcs(wcs, force_it=force_it,**kwargs)
        if not self.has_wcs():
            return
        if not self._wcsid:
            self._derive_radec_()
        else:
            self._derive_xy_()
            
    #  Super Catalogue
    # -----------------
    @_autogen_docstring_inheritance(CatalogueHandler.set_catalogue,"CatalogueHandler.get_indexes")
    def set_catalogue(self,catalogue, reset=True,
                      match_catalogue=True,match_angsep=3*units.arcsec, **kwargs):
        #
        # + reset and matching
        #
        if reset:
            catalogue.reset()

        _ = kwargs.pop("fast_setup", None)
        super(PhotoMap, self).set_catalogue(catalogue, fast_setup=True, **kwargs)
        
        if self.has_catalogue() and match_catalogue:
            self.match_catalogue(deltadist=match_angsep)

            
    @_autogen_docstring_inheritance(CatalogueHandler.download_catalogue,"CatalogueHandler.download_catalogue")
    def download_catalogue(self, source="sdss",
                           set_it=True,force_it=False,
                           radec=None, radius=None,
                           **kwargs):
        #
        # default definition of radec and radius
        #
        if radec is None:
            radec = np.mean(self.radec, axis=0)
            radec = "%s %s"%(radec[0],radec[1])
        if radius is None:
            radius = np.max(np.std(self.radec,axis=0)*5)
            
        return super(PhotoMap, self).download_catalogue(source=source,
                                                        set_it=set_it, force_it=force_it,
                                                        radec=radec, radius=radius,
                                                        **kwargs)
    # ------------- #
    # - GETTER    - #
    # ------------- #
    # -- Access Index
    def get_host_idx(self, target=None, coords=None,  catid=None,
                     radius=30, runits="kpc", scaleup=3, max_galdist=3,
                     verbose=False):

        """
        It first searches galaxies in a given radius and then returns
        the index of the one minimizing the elliptical radius (not necesseraly
        the nearest in angular distance).

        Parameters
        ----------

        // location
        
        target or coords: [astrobject.Target or [ra,dec]] -one required-
             provide the target location

        // host setting

        catid: [float/string/None] -optional-
            You can force the id for the host by providing it here.
            If None catid will be ignore and this will run the host-search
            
        // search options
        
        radius: [float/None] -optional-
            Distance used for the first galaxy search.
            
        runits: [string, astropy.units] -optional-
            unit of the radius.

        max_galdist: [float/None] -optional-
            Size (in ellipse radius) above which the target is
            assumed too far away for the nearest host and is consequently
            assumed hostless.
            Set this to None to ignore this test.

        // general options

        verbose: [bool]
            Print additional information
        
        Returns
        -------
        [int] / None # None if no host detected
        """
        if catid is not None:
            return self.catindex_to_index(self.catalogue.id_to_idx(catid))
        
        # --- input test --- #
        if target is None:
            if coords is None:
                raise ValueError("Both target and coords are None. Give one.")
            ra,dec = coords
        else:
            ra,dec = target.ra,target.dec
        if ra is None or dec is None:
            raise ValueError("Ra and/or Dec are/is None")

        # -- distances -- #
        from scipy.spatial import distance
        idxaround,idxaround_dist = self.get_idx_around(ra, dec,
                                                       radius*self.units_to_pixels(runits, target=target)*self.wcs.pix_indeg.value,
                                                       runits="degree")
        if verbose:
            print(idxaround,idxaround_dist)
            
        # --- No host around --- #
        if len(idxaround) == 0:
            print("No detected host within the given search limits: ", radius, runits)
            return None
            
        x,y,a,b,theta = self.get_ellipse_values(mask=idxaround)
        pix_x, pix_y  = self.coords_to_pixel(ra,dec)
        if len(idxaround) == 1:
            x,y,a,b,theta = [x],[y],[a],[b],[theta]
            
        dist = np.asarray([[idx_,distance.pdist([[0,0],
                                                 np.dot(np.asarray([[np.cos(theta[i]), -np.sin(theta[i])],[np.sin(theta[i]), np.cos(theta[i])]]).T,
                                                        [pix_x-x[i],pix_y-y[i]])], # coords_aligned = rotation matrix . coord_xy
                                                        "wminkowski", w=[1./(a[i]*scaleup),1./(b[i]*scaleup)])]
                                    for i,idx_ in enumerate(idxaround)
                                        if not self.has_catalogue() or \
                                        (idx_ in self.catmatch["idx"] and not \
                                         (self.catalogue.has_starmask() and np.all(self.catalogue.starmask[self.index_to_catindex([idx_])])))
                                    ]).T
        if len(dist) == 0:
            print("No detected host within the given search limits: ", radius, runits)
            return None
        
        # this is the composition of dist
        idx, dist_in_radius = dist
        
        i_nearest = np.argmin(dist_in_radius)
        if max_galdist is not None and dist_in_radius[i_nearest]>max_galdist:
            print("No nearby host identified"+(" for %s"%(target.name) if target is not None else ""))
            print(" Maximum allowed galaxy radius %.1f ; nearest galaxy [in its radius] %.1f"%(max_galdist, dist_in_radius[i_nearest]))
            return None
        return idx[i_nearest]

            
    def get_nearest_idx(self, ra, dec, wcs_coords=True,
                        catmatch=True, catalogue_idx=False,
                        **kwargs):
        """ get the index of the nearest object

        **kwargs goes to the catalogue's get_mask() method:
            Aviablable masking with there default values:
                - catmag_range=[None, None],
                - stars_only=False,
                - isolated_only=False,
                - nonstars_only=False,
                - contours=None,
                - notingalaxy=False
        """
        if not catmatch:
            raise NotImplementedError("only catalogue based associated so far")

        elif not self.has_catalogue():
            raise AttributeError("No catalogue set.")
        else:
            masking = kwargs_update({"matched":True},**kwargs)
            mask = self.catalogue.get_mask(**masking)
            id_, dist_  = self.catalogue.get_nearest_idx(ra, dec, wcs_coords=wcs_coords, mask=mask)
            if mask is not None:
                id_ = np.argwhere(mask)[id_[0]]
            
            return id_ if catalogue_idx else self.catindex_to_index(id_), dist_
        
    def get_idx_around(self, ra, dec, radius,runits="arcmin", catalogue_idx=False):
        """
        Parameters:
        ----------
        ra, dec : [float]          Position in the sky *in degree*

        radius:  [float]           Distance of search (see runits for unit)

        runits: [string]           Unit of the distance. (an astropy known units)
        
        Return
        ------
        2d-array (index, angular separation)
        """
        radius = radius*units.Unit(runits)
        
        if not is_arraylike(ra):
            ra = [ra]
            dec = [dec]
            
        sky = coordinates.SkyCoord(ra=ra*units.degree, dec=dec*units.degree)
        if not catalogue_idx:
            ra_,dec_ = self.radec.T
            return coordinates.SkyCoord(ra=ra_*units.degree,
                        dec=dec_*units.degree).search_around_sky(sky, radius)[1:3]
        
        if not self.has_catalogue():
            raise ValueError("no catalogue loaded, do not set catalogue_idx to True.")
            
        return self.catalogue.get_idx_around(ra, dec, radius.value,
                                             runits=radius.unit.name,
                                             wcs_coords=True)
        
    # -- Masking Tool
    def get_indexes(self,isolated_only=False,
                    stars_only=False,nonstars_only=False,
                    catmag_range=[None,None],
                    contours=None, cat_indexes=False):
        """
        Based on catalogue's 'get_mask()' method, this returns the list of indexes
        matching the input's criteria. Set cat_indexes to True to have the catalogue
        indexes references instead of the instance's one. 
        
        Returns
        -------
        array of indexes
        """
        if not self.has_catmatch():
            raise AttributeError("No catalogue matching. Run match_catalogue()")
        # -----------------------
        # - Catalogue Based cuts
        id_ = "idx_catalogue" if cat_indexes else "idx"
        return self.catmatch[id_][\
                        self.catalogue.get_mask(catmag_range=catmag_range,
                                                stars_only=stars_only,
                                                nonstars_only=nonstars_only,
                                                isolated_only=isolated_only,
                                                contours=contours, fovmask=True)[\
                                            self.catmatch["idx_catalogue"]]]
                                            
    
    
    # ------------- #
    # - PLOTTER   - #
    # ------------- #
    def display_voronoy(self,ax=None,toshow="flux",wcs_coords=False,
                        mask=None,
                        show_nods=False,**kwargs):
        """
        Show a Voronoy cell map colored as a function of the given 'toshow'.
        You can see the nods used to define the cells setting show_nods to True.    
        """
        import matplotlib.pyplot as mpl
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            ax.set_xlabel("x" if not wcs_coords else "Ra",fontsize = "x-large")
            ax.set_ylabel("y" if not wcs_coords else "Dec",fontsize = "x-large")
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure

        from ..utils.mpladdon import voronoi_patchs
        # -----------------
        # - The plot itself
        xy = self.radec if wcs_coords else self.xy
        out = ax.voronoi_patchs(xy if mask is None else xy[mask],
                                self.get(toshow,mask=mask),**kwargs)
        if show_nods:
            x,y = xy.T if mask is None else xy[mask].T
            ax.plot(x,y,marker=".",ls="None",color="k",
                    scalex=False,scaley=False,
                    alpha=kwargs.pop("alpha",0.8),zorder=kwargs.pop("zorder",3)+1)
        ax.figure.canvas.draw()
        return out

    # =============================== #
    #   Internal Tools                #
    # =============================== #
    def _derive_radec_(self):
        """ derives the ra dec values from the wcs solution """
        if not self.has_wcs():
            raise AttributeError("You need a wcs solution to convert pixels to radec. None set.")

        ra,dec = self.wcs.pix2world(*self.xy.T).T
        self.set_meta("ra",ra)
        self.set_meta("dec",dec)

    def _derive_xy_():
        """ derives the x y values from the wcs solution """
        if not self.has_wcs():
            raise AttributeError("You need a wcs solution to convert radec to pixels. None set.")
        x,y = self.wcs.world2pix(*self.radec.T).T
        self.set_meta("x",x)
        self.set_meta("y",y)
    # =============================== #
    #   properties                    #
    # =============================== #
    @property
    def radec(self):
        """ world coordinates in degree """
        if self._wcsid and not self.fromtable:
            return self._coords
        if "ra" not in self.getkeys or "dec" not in self.getkeys:
            self._derive_radec_()
            
        return self.get(["ra","dec"])
    
    @property
    def xy(self):
        """ image coordinates in pixels """
        if not self._wcsid and not self.fromtable:
            return self._coords
        
        if "x" not in self.getkeys or "y" not in self.getkeys:
            self._derive_xy_()
            
        return self.get(["x","y"])
    
    @property
    def _coords(self):
        """ given id that are coordinates for PhotoMaps. see self.radec and self.xy """
        return np.asarray([self.id_to_coords(id_) for id_ in self.list_id])
    
    @property
    def list_id(self):
        """ """
        if not self.fromtable:
            list_ = super(PhotoMap, self).list_id
            return np.asarray(list_)[np.argsort([float(l.split(",")[1]) for l in list_])] # y aligned
        
        return self.get(self._build_properties["idkey"])
    
    @property
    def _wcsid(self):
        if "wcsid" not in self._build_properties.keys() or \
          self._build_properties["wcsid"] is None:
            warnings.warn("No information about the nature of the coordinate ids (wcs/pixel?). **WCS system assumed**")
            self._build_properties["wcsid"] = True
        return self._build_properties["wcsid"]
    
    
    # -------------
    # - RefMap
    @property
    def refmap(self):
        return self._side_properties['refmap']
    
    def has_refmap(self):
        return self.refmap is not None
    
    # ------------- #
    # - Derived   - #
    # ------------- #
    @property
    def contours(self):
        from ...utils import shape
        return shape.get_contour_polygon(*self.radec.T)
    
    @property
    def contours_pxl(self):
        from ...utils import shape
        return shape.get_contour_polygon(*self.xy.T)
    
    # -------------
    # - Catalogue
    @property
    def catmatch(self):
        """This is the match dictionnary"""
        if self._derived_properties["catmatch"] is None:
            self._derived_properties["catmatch"] = {}
            
        return self._derived_properties["catmatch"]

    def has_catmatch(self):
        return (self.catmatch is not None and len(self.catmatch.keys())>0)
    
    @property
    def _fulldata(self):
        """ data in the saved format """
        fulldata = super(PhotoMap,self)._fulldata
        if self.has_catmatch():
            fulldata["catmatch"] = self.catmatch
        return fulldata
        
######################################
#                                    #
# Sextractor Output                  #
#                                    #
######################################
class SepObject( PhotoMap ):
    """ Child of PhotoMap that make uses of all the meta keys that are sextractor
    output """

    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["galaxy_contours","ingalaxy_mask"]
        
    # -------------------------- #
    # -  Plotting tools        - #
    # -------------------------- #
    def display(self,ax, scaleup=2.5,draw=True,
                apply_catmask=True,
                stars_only=False, isolated_only=False,
                catmag_range=[None,None], **kwargs):
        """ Show the ellipses of the sep extracted sources """
        
        mask = None if not apply_catmask or not self.has_catalogue() else\
          self.get_indexes(isolated_only=isolated_only,stars_only=stars_only,
                            catmag_range=catmag_range, cat_indexes=False)
                
        self.display_ellipses(ax, mask,scaleup=scaleup,
                              draw=draw, **kwargs )

    # - Ellipse
    def display_ellipses(self, ax, idx=None, scaleup=2.5,
                         draw=True, fc="None", ec="k"):
        """ draw the requested ellipse on the axis """
        ells = self.get_detected_ellipses(scaleup=scaleup,
                                          mask=idx if is_arraylike(idx) or idx is None\
                                           else [idx],
                                          contours=False)
        for ell in ells:
            ax.add_artist(ell)
            ell.set_clip_box(ax.bbox)
            ell.set_facecolor(fc)
            ell.set_edgecolor(ec)
        if draw:
            ax.figure.canvas.draw()

    def show_ellipses(self,ax=None,
                      savefile=None,show=True,
                      apply_catmask=True,stars_only=False,nonstars_only=False,
                      isolated_only=False,catmag_range=[None,None],
                      **kwargs):
        """ Display ellipses of the extracted sources (see masking options) """
        if not self.has_data():
            warnings.warn("WARNING [Sepobjects] No data to display")
            return
        
        from matplotlib.patches import Ellipse, Polygon
        import matplotlib.pyplot    as mpl
        from ..utils.mpladdon   import figout
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
          self.get_indexes(isolated_only=isolated_only,
                           stars_only=stars_only,nonstars_only=nonstars_only,
                        catmag_range=catmag_range)
        # -------------
        # - Properties
        ells = [Ellipse([0,0],2.,2*b/a,t*units.radian.in_units("degree"))
                for a,b,t in self.get(["a","b","theta"],mask=mask)]
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


    # ----------- #
    #   GETTER    #
    # ----------- #
    def get_skycoords(self, position="average", default_radec=True):
        """ astropy SkyCoord object corresponding to ra and dec positions

        Parameters
        ----------
        default_radec: [bool] -optional-
            If there are 'ra' and 'dec' entries on the meta, should this use them?
            Remark: if no wcs solution loaded, this is forced to True
            
        position: [bool] -optional-
            If there is a wcs solution and if this has to build the ra and dec values
            how this should be made?
            could be:
                - center:  This uses the center of the ellipse 
                - peak:    This uses the center of mass of the ellipses
                - average: The average of 'center' and 'peak'
        Returns
        -------
        astropy's SkyCoord
        """
        if not default_radec and not self.has_wcs():
            warnings.warn("Without wcs solution, you need to have ra and dec in the meta entries")
            default_radec = True
            
        if "ra" in self.getkeys and "dec" in self.getkeys and default_radec:
            ra,dec = self.get("ra"), self.get("dec")
        elif not self.has_wcs():
            raise AttributeError("no wcs solution loaded and no 'ra', 'dec' entries in the meta")
        else:
            x,y,xpeak,ypeak = self.get(["x", "y", "xpeak", "ypeak"]).T
            if position.lower() == "center":
                ra,dec = self.pixel_to_coords(x,y).T
            elif position.lower() == "peak":
                ra,dec =self.pixel_to_coords(xpeak,ypeak).T 
            elif position.lower() in ["mean", "average"]:
                ra,dec = self.pixel_to_coords(np.mean([x, xpeak], axis=0),np.mean([y, ypeak], axis=0)).T
            else:
                raise ValueError("Unknown matching position %s"%position)
        
        return coordinates.SkyCoord(ra=ra*units.degree,dec=dec*units.degree)

    # -------------------------- #
    # Super It                 - #
    # -------------------------- #
    @_autogen_docstring_inheritance(PhotoMap.get_indexes,"PhotoMap.get_indexes")
    def get_indexes(self,isolated_only=False,stars_only=False,nonstars_only=False,
                    catmag_range=[None,None],contours=None,
                    notwithin_galaxies=False,
                    cat_indexes=False):
        #
        # Add the notwithin_galaxies masking
        #
        catmasking = {"isolated_only":isolated_only,
                      "stars_only":stars_only,"nonstars_only":nonstars_only,
                       "catmag_range":catmag_range,"contours":contours}
        # ---------------------------
        # - Catalogue based masking
        mask = super(SepObject,self).get_indexes(cat_indexes=cat_indexes,
                                                 **catmasking)
        # ---------------------------
        # - SEP based masking
        if notwithin_galaxies:
            from ...utils.shape import point_in_contours
            masking = mask if not cat_indexes else \
              self.get_indexes(notwithin_galaxies=False,cat_indexes=False,
                               **catmasking)
            mask = mask[~self._ingalaxy_mask[masking]]
            
        # -----------
        # - Returns
        return mask
        
    # -------------------------- #
    # Other gets               - #
    # -------------------------- #
    def get_fwhm_pxl(self,stars_only=True,isolated_only=True,
                    catmag_range=[None,None], default_around=10*units.arcsec):
        """
        FWHM estimated defined as the median of 2 * sqrt(ln(2) * (a^2 + b^2))
        in units of a/b ; must be pixels
        """
        if isolated_only and self.has_catalogue() and not self.catalogue._is_around_defined():
            warnings.warn("No isolation defined. default value set: %s"%default_around)
            self.catalogue.define_around(default_around)
            
        mask = self.get_indexes(isolated_only=isolated_only,
                                stars_only=stars_only,
                                catmag_range=catmag_range)
        return np.median(2 * np.sqrt( np.log(2) * (self.get('a',mask)**2 + self.get('b',mask)**2)))

    
    def get_median_ellipse(self,mask=None,clipping=[3,3]):
        
        """ This methods look for the stars and return the mean ellipse parameters 

        Returns
        -------
        [a, std_a], [b, std_b], [theta, std_theta]
        """
        from scipy.stats import sigmaclip
        if not self.has_catalogue():
            warnings.warn("No catalogue loaded, not catalogue masking avialable")
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
        
    def get_detected_ellipses(self, scaleup=2.5, mask=None, contours=False):
        """
        Get the matplotlib Patches (Ellipse) defining the detected object. You can
        select the returned ellipses using the apply_catmask, stars_only,
        isolated_only and catmag_range cuts.

        'scaleup' scale the radius used of the ellipses.
        2.5 means that most of the visible light will be within the inside the returned ellipses.

        'contours' means that the returned value is a shapely MultiPolgon not a list of patches
        Return
        ------
        list of patches
        """
        from matplotlib.patches import Ellipse
        
        # -------------
        # - Properties
        x_,y_,a_,b_,t_ = self.get_ellipse_values(mask=mask)
        ells = [Ellipse([x,y],a*scaleup*2,b*scaleup*2,
                            t*units.radian.in_units("degree"))
                    for x,y,a,b,t in zip(x_,y_,a_,b_,t_)]
                
        if not contours:
            return ells
        from ...utils.shape import patch_to_polygon
        return patch_to_polygon(ells)

    def get_ellipse_values(self, mask=None):
        """ returns the x, y, a, b, theta for the sep objects """
        if mask is not None and len(mask)==1:
            return self.get(["x","y","a","b","theta"], mask=mask)[0]
        return self.get(["x","y","a","b","theta"], mask=mask).T

    def project_ellipse(self, idx, target_wcs):
        """ grab the ellipse coordinate of the given idx entry and project them
        in the new wcs system using self.project_ellipse_to_wcs() """
        x,y,a,b,theta = self.get_ellipse_values(idx if is_arraylike(idx) else [idx])
        return self.project_ellipse_to_wcs(x,y,a,b,theta, target_wcs)
    
    def project_ellipse_to_wcs(self, x,y,a,b,theta, target_wcs):
        """ Get the x,y, a, b, theta values in the new wcs system. """
        xa,ya = np.cos(theta)*a + x, np.sin(theta)*a + y
        xb,yb = np.cos(theta-3*np.pi/2)*b + x, np.sin(theta-3*np.pi/2)*b + y
        # - Move to Ra Dec where things are the same for everyone
        centroid, vertex, covertex = self.wcs.pix2world([x,xa,xb],[y,ya,yb])
        # - Move back this points using the new wcs
        x_new, y_new  = target_wcs.world2pix(*centroid)
        xa_new,ya_new = target_wcs.world2pix(*vertex)
        xb_new,yb_new = target_wcs.world2pix(*covertex)
        # - convert back to x,y,a,b,theta
        a_new = np.sqrt((xa_new-x_new)**2+(ya_new-y_new)**2)
        b_new = np.sqrt((xb_new-x_new)**2+(yb_new-y_new)**2)
        theta_new = np.arctan((ya_new-y_new)/(xa_new-x_new))
        return x_new, y_new, a_new, b_new, theta_new
        
    def get_ellipse_mask(self,width, height, r=3, apply_catmask=False):
        """
        This method returns a boolean mask of the detected ellipses
        on the given width x height pixels image

        (this method is based on the sep mask_ellipse function)
        
        Parameters:
        -----------
        r: [float]                 The scale of the ellipse (r=1 is a typical
                                   contours included the object ; 2 enables to
                                   get the tail of most of the bright sources

        apply_catmask: [bool]      Only mask the detected object associated with
                                   the current catalogue. If no catalogue loaded,
                                   this will be set to False in any case.

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
    
    # =============================== #
    # = properties                  = #
    # =============================== #
    @property
    def galaxy_contours(self, scaleup=5):
        """ """
        if self._derived_properties["galaxy_contours"] is None:
            self._derived_properties["galaxy_contours"] = \
              self.get_detected_ellipses(contours=True,
                                         scaleup=scaleup,
                                         mask = self.get_indexes(nonstars_only=True,
                                                                 stars_only=False,
                                                                 cat_indexes=False))
        return self._derived_properties["galaxy_contours"]
        
    @property
    def _ingalaxy_mask(self):
        """ the array masking the detected sources with or without of galaxies"""
        if not self.has_catalogue():
            raise AttributeError("No Catalogue loaded. Requested for the galaxy mask")
        
        if self._derived_properties["ingalaxy_mask"] is None:
            from ...utils.shape import Point
            # -- This won't do anything if this is already loaded
            #    otherwise, loads the catalogue with the ingalaxy info
            self.catalogue.set_ingalaxymask(self.galaxy_contours)
            # -- Ok let's build the mask, first, from the Catalogue values
            galaxymask = np.zeros(self.nsources)
            catgalmask = np.asarray(self.catalogue.ingalaxymask,dtype="bool")
            idxIn  = self.catmatch["idx"][ catgalmask[self.catmatch["idx_catalogue"]]]
            idxOut = self.catmatch["idx"][~catgalmask[self.catmatch["idx_catalogue"]]]
            idxtbd = np.asarray([i for i in range(self.nsources) if i not in
                      idxIn.tolist()+idxOut.tolist()])
        
            listingal = idxtbd[\
                np.asarray([self.galaxy_contours.contains(p_) for p_ in \
                [Point(x_,y_) for x_,y_ in self.get(["x","y"],mask=idxtbd)]
                ],dtype=bool)]
            
            galaxymask[idxIn] = True
            galaxymask[idxOut] = False
            galaxymask[idxtbd] = False
            galaxymask[listingal] = True # this is a subpart of idxtbd
            self._derived_properties["ingalaxy_mask"] = np.asarray(galaxymask,
                                                                   dtype="bool")
              
        return self._derived_properties["ingalaxy_mask"]




######################################
#                                    #
#   Collection of PhotoMaps          #
#                                    #
######################################

class PhotoMapCollection( BaseCollection, CatalogueHandler ):
    """ Collection of PhotoMaps """

    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    def __init__(self, photomaps=None, empty=False, catalogue=None, **kwargs):
        """ """
        self.__build__()
        if empty:
            return

        # - Load the Catalogue if given
        if catalogue is not None:
            self.set_catalogue(catalogue, **kwargs)
            
        # - Load any given photomap
        if photomaps is not None:
            if is_arraylike(photomaps):
                [self.add_photomap(photomap, **kwargs) for photomap in photomaps]
            else:
                self.add_photomap(photomaps, **kwargs)


    def writeto(self, filename, include_catatogue=True):
        """ Dumpt the data (could including the catalog) in the given file """
        if include_catatogue:
            dump_pkl(self._fulldata, filename)
        else:
            dump_pkl({"photomaps" : {id_:self.photomaps[id_]._fulldata for id_ in self.list_id},
                      "catalogue" : None}, filename)

    def load(self, filename, match_catalogue=False, **kwargs):
        """ load the data from the given file
         **kwargs goes to the match_catalogue()
         """
        
        datain = load_pkl(filename)
        if "catalogue" in datain.keys() and datain["catalogue"] is not None:
            from ..instruments.instrument import get_catalogue
            self.set_catalogue( get_catalogue( table.Table(datain["catalogue"]), None ) )
            
        if "photomaps" in datain.keys():
            for k,v in datain["photomaps"].items():
                pmap = PhotoMap() if "a" not in v["data"] or "b" not in v["data"] else SepObject()
                pmap._load_pkl_(fulldata=v)
                self.add_photomap(pmap, id_=k, match_catalogue=False)

        if match_catalogue:
            self.match_catalogue( **kwargs)
            
    # =================== #
    #   Main Methods      #
    # =================== #
    def add_photomap(self, photomap, id_=None,
                     match_catalogue=True, **kwargs):
        """ 
        This method enables to load the given photomap in the
        self.photomaps container (dict).
        
        Parameters
        ----------
        new_image: [string or astrobject's Image (or children)]

        - option -
        id_: [any]
            key used to access the photomap from the _handler.
            id_ is set the photopoint's bandname if None

        **kwargs goes to catalogue matching.
        Return
        ------
        Void
        """
        # --------------
        # - Define ID
        if photomap is None:
            return

        if PhotoMap not in photomap.__class__.__mro__:
            raise TypeError("The givem photomap must be (or inherate of) an astrobject's PhotoMap")
        
        if self.has_sources():
            if id_ is None or id_ in self.list_id:
                id_= ""
                i = 1
                while id_+"%d"%i in self.list_id:
                    i+=1
                id_ = id_+"%d"%i
        else:
            id_ = "1"
            
        # -------------------
        # - Match the catalogue
        if self.has_catalogue() and match_catalogue:
            photomap.match_catalogue(force_it=True, **kwargs)
        # -------------------
        # - Record the map            
        self.photomaps[id_]  = photomap

    # -------------- #
    #    GETTER      #
    # -------------- #
    def getcat(self, param, catindex):
        """ get method based on catindex entry(ies). """
        return np.asarray([self.photomaps[id_].getcat(param, catindex)
                for id_ in self.list_id])
    
    # -------------- #
    #  Catalogue     #
    # -------------- #
    def match_catalogue(self, deltadist=3*units.arcsec, **kwargs):
        """ Match the catalogue to all the PhotoMaps """
        if not self.has_catalogue():
            raise AttributeError("No catalogue to match")
        
        [self.photomaps[id_].match_catalogue(deltadist=3*units.arcsec,
                                              **kwargs)
         for id_ in self.list_id]

    def get_catindexes(self, inclusion=0,
                       isolated_only=False, stars_only=False,
                       nonstars_only=False, catmag_range=[None, None],
                       contours=None, **kwargs):
        """ Get the index of the catalogue that has been mathced.

        Parameters
        ----------
        inclusion: [0,1 or >1] -optional-
            - inclusion [0,1]:
                Minimal fraction of photomaps that should share a catalogue entry.
                If 0 this means that any catalogue entry detected at least once
                will be returned. If 1, only catalogue entry shared by all Photomaps
                will be returned.
            - inclusion >1:
                Minimal number of time a given catalogue entry is matched within a photomaps
                For instance if include is 3, only catalogue entries detected at least 3 times
                will be used.

        Returns
        -------
        list (indexes of the catalogue)
        """
        catindexes = [self.photomaps[id_].get_indexes(isolated_only=isolated_only, stars_only=stars_only,
                                                      nonstars_only=nonstars_only, catmag_range=catmag_range,
                                                      contours=contours,
                                                      cat_indexes=True, **kwargs)
                      for id_ in self.list_id]
            
        # For speed, let's do the obvious cases first
        if inclusion==1:
            return list(frozenset(catindexes[0]).intersection(*catindexes[1:]))
        if len(catindexes) == 1:
            return catindexes[0]

        detected_once = np.unique(np.concatenate(catindexes))
        if inclusion==0:
            return detected_once
        
        # Ok then let's measure it.
        fractionin = np.sum([self.photomaps[i].is_catindex_detected(detected_once) for i in self.list_id], axis=0)/float(self.nsources) if inclusion<1\
          else np.sum([self.photomaps[i].is_catindex_detected(detected_once) for i in self.list_id], axis=0)

        return detected_once[fractionin>=inclusion]        

    def catindex_to_index(self, catindex, cleanindex=False):
        """ Get the index of the individual photomaps associated to the given catindex
        
        Parameters
        ----------
        catindex: [int of list of]
            catalogue entry that should be matched with the photomap indexes

            
        Returns
        -------
        list [[idx],[idx],...]
        """
        return [self.photomaps[id_].catindex_to_index(catindex,cleanindex=cleanindex)
                for id_ in self.list_id]

    
    # - Super Mother CatalogueHandler
    @_autogen_docstring_inheritance(CatalogueHandler.download_catalogue,"CatalogueHandler.download_catalogue")
    def download_catalogue(self, source="sdss",
                           set_it=True,force_it=False,
                           radec=None, radius=None,
                           match_catalogue=True,match_angsep=2*units.arcsec,
                           **kwargs):
        #
        # default definition of radec and radius
        #
        if radec is None or radius is None:
            ra,dec = np.concatenate(self.get(["ra","dec"]), axis=0).T
            
        if radec is None:
            radec = np.mean([ra,dec],axis=1)
            radec = "%s %s"%(radec[0],radec[1])
        if radius is None:
            radius = np.max(np.std([ra,dec],axis=1)*5)
            
        out = super(PhotoMapCollection, self).download_catalogue(source=source,
                                                        set_it=set_it, force_it=force_it,
                                                        radec=radec, radius=radius,
                                                        **kwargs)
        if self.has_catalogue() and match_catalogue:
            self.match_catalogue(deltadist=match_angsep)

        return out
    
    @_autogen_docstring_inheritance(CatalogueHandler.set_catalogue,"CatalogueHandler.set_catalogue")
    def set_catalogue(self, catalogue, force_it=False, **kwargs):
        # 
        super(PhotoMapCollection, self).set_catalogue( catalogue, force_it=force_it, **kwargs)
        [self.photomaps[id_].set_catalogue(self.catalogue, fast_setup=True, **kwargs)
         for id_ in self.list_id]

                
    # =================== #
    #   Properties        #
    # =================== #
    @property
    def photomaps(self):
        """ """
        return self._handler

    @property
    def _fulldata(self):
        return {"photomaps" : {id_:self.photomaps[id_]._fulldata for id_ in self.list_id},
                "catalogue" : np.asarray(self.catalogue.data)}
