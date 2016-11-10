#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np

from astropy        import coordinates, units

from ..baseobject   import WCSHandler, CatalogueHandler
from ..photometry   import get_photopoint
from ..collection   import PhotoPointCollection
from ..utils.tools import kwargs_update
from ..utils.decorators import _autogen_docstring_inheritance
__all__ = ["get_photomap","get_sepobject"]



def get_photomap(photopoints=None,coords=None,wcs_coords=True,**kwargs):
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

    Return
    ------
    PhotoMap (Collection of PhotoPoints)
    """
    return PhotoMap(photopoints=photopoints,
                    coords=coords,
                    wcs_coords=wcs_coords,
                    **kwargs)


def get_sepobject(sepoutput, ppointkwargs={},
                  use_peakposition=False,
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

    **kwargs goes to SepObject init (Child of PhotoMap)

    Return
    ------
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
        
    return pmap
    
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
    from astropy import table
    if type(array) != np.ndarray:
        return False
    tsep = table.Table(array)
    return "cxx" in tsep.keys()

def parse_sepoutput(sepoutput,lbda=None,**kwargs):
    """
    """
    from astropy import table
    if type(sepoutput) != np.ndarray and type(sepoutput) != dict:
        raise TypeError("the given 'sexoutput' is not an ndarray ; This is not a sepoutput")

    # - so that any entry now has the same shape
    tsep = table.Table(sepoutput)
    ppoints = [get_photopoint(lbda,t_["flux"],None,
                          source="sepextract",**kwargs)
                          for t_ in tsep]
        
    coords = np.asarray([tsep["x"],tsep["y"]]).T
    return {"ppoints":ppoints,"coords":coords,"meta":tsep[[t_ for t_ in tsep.keys() if t_ not in ["flux"]]]}
    
# ------------------- #
#  From Sextractor    #
# ------------------- #
def parse_sepsex_fitsdata( filename, lbda=None, fluxkey="flux_auto", **kwargs):
    """ This checks is the fits file is from sextractor
    Parameters
    ----------
    
    dataindex: [int/None] -optional-
        
    """
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

    
    def __init__(self, photopoints=None,coords=None,
                 filein=None,empty=False,wcs=None,catalogue=None,
                 **kwargs):
        """
        """
        self.__build__()
        if empty:
            return
        if filein is not None:
            self.load(filein,**kwargs)
            
        if photopoints is not None:
            self.create(photopoints,coords,wcs=wcs,catalogue=catalogue,**kwargs)

    def _read_table_comment_(self, comments):
        """
        """
        for c in comments:
            key,value = c.split()
            if key == "wcsid":
                self._build_properties["wcsid"]=bool(value)
            else:
                print "comment not saved: ", c
                
    # =============================== #
    # = Building Methods            = #
    # =============================== #
    # ---------- #
    #  GETTER    #
    # ---------- #
    def get_skycoords(self):
        """ """
        ra,dec = self.radec.T
        return coordinates.SkyCoord(ra=ra*units.degree,dec=dec*units.degree)
    
    # ---------------- #
    # - SUPER THEM   - #
    # ---------------- #
    def create(self,photopoints,coords, wcs_coords=True,
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
            
    def writeto(self,filename,format="ascii.commented_header",**kwargs):
        comments = "".join(["# %s %s \n "%(key, value) for key,value in [["wcsid",self._wcsid]]])
        super(PhotoMap,self).writeto(filename,format,comment=comments,**kwargs)
        
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
        catindex = cat_indexes[np.in1d( cat_indexes, self.catmatch["idx_catalogue"] )]
        index = self.catindex_to_index(catindex)

        # ------------------- #
        # - Get Photometry  - #
        # ------------------- #        
        return catindex,index
    
    # -- Id <-> Coords
    def coords_to_id(self,x,y):
        """ this enables to have consistent coordinate/id mapping """
        return "%.8f,%.8f"%(x,y)
    
    def id_to_coords(self,id):
        """ convert photopoint's id into coordinates. Only reading here, see radec, xy """
        return np.asarray(id.split(","),dtype="float")

    # -- Catalogue <-> PhotoMap
    def match_catalogue(self,catalogue=None,
                        force_it=False, deltadist=5*units.arcsec,
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
        
        # -- matching are made in degree space
        skyradec = self.get_skycoords(**kwargs)
        idxcatalogue, idxsepobjects, d2d, d3d = skyradec.search_around_sky(self.catalogue.sky_radec, deltadist)
        
        # --------------------
        # - Save the results
        self._derived_properties["catmatch"] = {
            "idx_catalogue":idxcatalogue,
            "idx":idxsepobjects,
            "angsep":d2d
            }
        
        self.catalogue.set_matchedmask(idxcatalogue)

    def index_to_catindex(self,index, cleanindex=False):
        """ give the catalogue index corresponding to the given photomap index.
        The cleanindex option enable to automatically remove the unknown index values.
        If you do not. this will raise a ValueError indicating the unknown index
        """
        if "__iter__" not in dir(index):
            catindex = [index]
        index = np.asarray(index)
        
        bool_ = np.in1d( index, self.catmatch["idx"] )
        if not np.all(bool_) and  not cleanindex:
            raise ValueError("Unknown photomap index(es):"+", ".join(index[bool_]))
            
        return np.concatenate(self.catmatch["idx_catalogue"][\
                        np.argwhere(np.in1d(self.catmatch["idx"], index[bool_]))])
                        
        
    def catindex_to_index(self,catindex, cleanindex=False):
        """ give the photomap index corresponding to the given catalogue index.
        The cleanindex option enable to automatically remove the unknown index values.
        If you do not. this will raise a ValueError indicating the unknown index
        """
        if "__iter__" not in dir(catindex):
            catindex = [catindex]
        catindex = np.asarray(catindex)
        
        bool_ = np.in1d( catindex, self.catmatch["idx_catalogue"] )
        if not np.all(bool_) and not cleanindex:
            raise ValueError("Unknown catalogue index(es):"+", ".join(catindex[bool_]))
            
        return np.concatenate(self.catmatch["idx"][\
                        np.argwhere(np.in1d(self.catmatch["idx_catalogue"], catindex[bool_]))])

    
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

    #  Super Catalogue
    # -----------------
    @_autogen_docstring_inheritance(CatalogueHandler.set_catalogue,"CatalogueHandler.get_indexes")
    def set_catalogue(self,catalogue,force_it=True,
                      reset=True,
                      match_catalogue=True):
        #
        # + reset and matching
        #
        if reset:
            catalogue.reset()
            
        super(PhotoMap, self).set_catalogue(catalogue,force_it=True)
        
        if self.has_catalogue() and match_catalogue:
            self.match_catalogue()

            
    @_autogen_docstring_inheritance(CatalogueHandler.download_catalogue,"CatalogueHandler.download_catalogue")
    def download_catalogue(self, source="sdss",
                           set_it=True,force_it=False,
                           radec=None, radius_degree=None,
                           **kwargs):
        #
        # default definition of radec and radius_degree
        #
        if radec is None:
            radec = np.mean(self.radec, axis=0)
            radec = "%s %s"%(radec[0],radec[1])
        if radius_degree is None:
            radius_degree = np.max(np.std(self.radec,axis=0)*5)
            
        return super(PhotoMap, self).download_catalogue(source="sdss",
                                                        set_it=set_it, force_it=force_it,
                                                        radec=radec, radius_degree=radius_degree,
                                                        **kwargs)
    # ------------- #
    # - GETTER    - #
    # ------------- #
    # -- Access Index
    def get_host_idx(self, target=None, coords=None,  catid=None,
                     radius=30, runits="kpc", scaleup=3, max_galdist=3):

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
        # --- No host around --- #
        if len(idxaround) == 0:
            print("No detected host within the given search limits: ", radius, runits)
            return None
        
        x,y,a,b,theta = self.get(["x","y","a","b","theta"], mask = idxaround).T
        pix_x, pix_y = self.coords_to_pixel(ra,dec)

        dist = np.asarray([[idx_,distance.pdist([[0,0],
                                                 np.dot(np.asarray([[np.cos(theta[i]), -np.sin(theta[i])],[np.sin(theta[i]), np.cos(theta[i])]]).T,
                                                        [pix_x-x[i],pix_y-y[i]])], # coords_aligned = rotation matrix . coord_xy
                                                        "wminkowski", w=[1/(a[i]*scaleup),1/(b[i]*scaleup)])]
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
            print "No nearby host identified"+(" for %s"%(target.name) if target is not None else "")
            print " Maximum allowed galaxy radius %.1f ; nearest galaxy [in its radius] %.1f"%(max_galdist, dist_in_radius[i_nearest])
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
        
        if "__iter__" not in dir(ra):
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

    # -- Ellipse Shape
    def get_idx_ellipse(self, idx):
        """ get the ellipse parameters ("x","y","a","b","theta") for the given idx.
        
        idx could be a single idx, a list of index or a bolean mask.
        """
        return self.get(["x","y","a","b","theta"],
                        mask=idx if hasattr(idx, "__iter__") or idx is None else [idx])[0]
        
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

        from ...utils.mpladdon import voronoi_patchs
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
    # = properties                  = #
    # =============================== #
    @property
    def radec(self):
        """ """
        if self._wcsid:
            return self._coords
        
        if not self.has_wcs():
            if "ra" in self.metakeys and "dec" in self.metakeys:
                return self.get(["ra","dec"])
            else:
                raise AttributeError("You need a wcs solution to convert pixels to radec or to have 'ra' and 'dec' in the metakeys. None set.")
            
        return self.wcs.pix2world(*self._coords.T)
    
    @property
    def xy(self):
        """ """
        if not self._wcsid:
            return self._coords
        if not self.has_wcs():
            raise AttributeError("You need a wcs solution to convert pixels to radec. None set.")
        return self.wcs.world2pix(*self._coords.T)
    
    @property
    def _coords(self):
        """ given id that are coordinates for PhotoMaps. see self.radec and self.xy """
        return np.asarray([self.id_to_coords(id_) for id_ in self.list_id])
    
    @property
    def _wcsid(self):
        if self._build_properties["wcsid"] is None:
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
                         draw=True,
                         fc="None", ec="k"):
        """ draw the requested ellipse on the axis """
        ells = self.get_detected_ellipses(scaleup=scaleup,
                                          mask=idx if hasattr(idx, "__iter__") or idx is None\
                                           else [idx],
                                          contours=False)
        for ell in ells:
            ell.set_clip_box(ax.bbox)
            ell.set_facecolor(fc)
            ell.set_edgecolor(ec)
            ax.add_patch(ell)
        if draw:
            ax.figure.canvas.draw()
        
    def show_ellipses(self,ax=None,
                      savefile=None,show=True,
                      apply_catmask=True,stars_only=False,nonstars_only=True,
                      isolated_only=False,catmag_range=[None,None],
                      **kwargs):
        """ Display ellipses of the extracted sources (see masking options) """
        if not self.has_data():
            print "WARNING [Sepobjects] No data to display"
            return
        
        from matplotlib.patches import Ellipse,Polygon
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
    def get_skycoords(self, position="average", meta_radec=True):
        """ astropy SkyCoord object corresponding to ra and dec positions

        Parameters
        ----------
        meta_radec: [bool] -optional-
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
        if not meta_radec and not self.has_wcs():
            warnings.warn("Without wcs solution, you need to have ra and dec in the meta entries")
            meta_radec = True
            
        if "ra" in self.metakeys and "dec" in self.metakeys and meta_radec:
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
        
        """ This methods look for the stars and return the mean ellipse parameters """
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
        
    def get_detected_ellipses(self,scaleup=2.5, mask=None, contours=False):
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
        ells = [Ellipse([x,y],a*scaleup*2,b*scaleup*2,
                        t*units.radian.in_units("degree"))
                for x,y,a,b,t in self.get(["x","y","a","b","theta"],mask=mask)]
        if not contours:
            return ells
        from ...utils.shape import patch_to_polygon
        return patch_to_polygon(ells)
    
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
