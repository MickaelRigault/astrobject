#! /usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np

from astropy        import coordinates, units

from ..photometry   import get_photopoint
from ..collection   import PhotoPointCollection
from ..utils.tools import kwargs_update
from ..utils.decorators import _autogen_docstring_inheritance
__all__ = ["get_photomap","get_sepobject"]



def get_photomap(photopoints=None,coords=None,wcs_coords=True,**kwargs):
    """
    """
    return PhotoMap(photopoints=photopoints,
                    coords=coords,
                    wcs_coords=wcs_coords,
                    **kwargs)


def get_sepobject(sepoutput, wcs_coords=False,ppointkwargs={},
                   **kwargs):
    """
    """
    if is_sepoutput(sepoutput):
        inputphotomap = parse_sepoutput(sepoutput,**ppointkwargs)
        pmap = SepObject(photopoints=inputphotomap["ppoints"],
                         coords=inputphotomap["coords"],
                         wcs_coords=wcs_coords, **kwargs)
        ids_sorting = [pmap.coords_to_id(x,y)
                       for x,y in zip(inputphotomap["meta"]["x"].data,
                                      inputphotomap["meta"]["y"].data)]
        [pmap.set_meta(k,inputphotomap["meta"][k].data, ids=ids_sorting)
         for k in inputphotomap["meta"].keys()]
        
        return pmap
        
# ======================= #
#                         #
# Internal Functions      #
#                         #
# ======================= #
def is_sepoutput(array):
    """
    """
    from astropy import table
    if type(array) != np.ndarray:
        return False
    tsep = table.Table(array)
    if "cxx" not in tsep.keys():
        return False
    return True

def parse_sepoutput(sepoutput,lbda=None,**kwargs):
    """
    """
    from astropy import table
    if type(sepoutput) != np.ndarray:
        raise TypeError("the given 'sexoutput' is not an ndarray ; This is not a sepoutput")

    tsep = table.Table(sepoutput)
    if "cxx" not in tsep.keys():
        raise TypeError("the given 'sexoutput' ndarray has no 'cxx' key."+\
                            "\n"+" It most likely is not a sextrator/sep output ")
    ppoints = [get_photopoint(lbda,t_["flux"],None,
                          source="sepextract",**kwargs)
                          for t_ in tsep]
    coords = np.asarray([tsep["x"],tsep["y"]]).T
    return {"ppoints":ppoints,"coords":coords,"meta":tsep[[t_ for t_ in tsep.keys() if t_ not in ["flux"]]]}
    
    
    


######################################
#                                    #
# PhotoMetric Mapping                #
#                                    #
######################################

class PhotoMap( PhotoPointCollection ):
    """
    """
    #__nature__ = "PhotoMap"

    PROPERTIES         = []
    SIDE_PROPERTIES    = ["refmap","wcs","catalogue"]
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
                        force_it=False, deltadist=1*units.arcsec):
        """
        This methods enable to attached a given sepobject entry
        to a catalog value.
        You can set a catalogue.
        """
        # --------------
        # - input 
        if catalogue is not None:
            self.set_catalogue(catalogue, force_it=force_it)

        if not self.has_catalogue():
            raise AttributeError("No 'catalogue' defined or given.")
        
        if self.has_wcs():
            # -- matching are made in degree space
            ra,dec = self.radec.T
            skyradec = coordinates.SkyCoord(ra=ra*units.degree,dec=dec*units.degree)
            idxcatalogue, idxsepobjects, d2d, d3d = skyradec.search_around_sky(self.catalogue.sky_radec, deltadist)
        else:
            raise NotImplementedError("You currently need a wcs solution in the SepObject to match a catalogue")
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
            raise ValueError("Unknown photomap index(es):"+","%join(index[bool_]))
            
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
            raise ValueError("Unknown catalogue index(es):"+","%join(catindex[bool_]))
            
        return np.concatenate(self.catmatch["idx"][\
                        np.argwhere(np.in1d(self.catmatch["idx_catalogue"], catindex[bool_]))])

    
    # ------------- #
    # - SETTER    - #
    # ------------- #
    def set_wcs(self,wcs,force_it=False):
        """
        """
        if self.has_wcs() and force_it is False:
            raise AttributeError("'wcs' is already defined."+\
                    " Set force_it to True if you really known what you are doing")
                    
        from ..astrometry import get_wcs
        self._side_properties["wcs"] = get_wcs(wcs,verbose=True)

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

    def set_catalogue(self,catalogue,force_it=True,
                      reset=True,
                      match_catalogue=True):
        """
        Attach to this instance a catalogue.
        
        Parameters
        ---------
        Return
        ------
        Voids
        """
        if self.has_catalogue() and force_it is False:
            raise AttributeError("'catalogue' already defined"+\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(catalogue) or catalogue.__nature__ != "Catalogue":
            raise TypeError("the input 'catalogue' must be an astrobject catalogue")

        # ----------------------
        # - Clean the catalogue 
        if reset: catalogue.reset()
        
        # -------------------------
        # - Add the world_2_pixel            
        if self.has_wcs():
            if catalogue.has_wcs():
                warnings.warn("WARNING the given 'catalogue' already has a wcs solution."+\
                              " This did not overwrite it.")
            else:
                catalogue.set_wcs(self.wcs)
                
        # ---------------------
        # - Test consistancy
        if catalogue.nobjects_in_fov < 1:
            warnings.warn("WARNING No object in the FoV, catalogue not loaded")
            return
        
        # ---------------------
        # - Good To Go
        self._side_properties["catalogue"] = catalogue
        if match_catalogue:
            self.match_catalogue()

    # ------------- #
    # - GETTER    - #
    # ------------- #
    # -- Access Index
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
            raise AttributeError("You need a wcs solution to convert pixels to radec. None set.")
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
    # - WCS
    @property
    def wcs(self):
        return self._side_properties['wcs']
    
    def has_wcs(self):
        return self.wcs is not None
    
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
    def catalogue(self):
        return self._side_properties["catalogue"]
    
    def has_catalogue(self):
        return self.catalogue is not None

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
    def display(self,ax,draw=True,apply_catmask=True,
                stars_only=False, isolated_only=False,
                catmag_range=[None,None]):
        """ Show the ellipses of the sep extracted sources """
        mask = None if not apply_catmask else\
          self.get_indexes(isolated_only=isolated_only,stars_only=stars_only,
                            catmag_range=catmag_range, cat_indexes=False)
                
        ells = self.get_detected_ellipses(scaleup=5, mask=mask,
                                          contours=False)
        for ell in ells:
            ell.set_clip_box(ax.bbox)
            ell.set_facecolor("None")
            ell.set_edgecolor("k")
            ax.add_patch(ell)
        if draw:
            ax.figure.canvas.draw()

    # - Ellipse
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
        from ...utils.mpladdon import figout
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
        
    def get_detected_ellipses(self,scaleup=5,mask=None, contours=False):
        """
        Get the matplotlib Patches (Ellipse) defining the detected object. You can
        select the returned ellipses using the apply_catmask, stars_only,
        isolated_only and catmag_range cuts.

        'scaleup' scale the radius used of the ellipses. 5 means that most of the
        visible light will be within the inside the returned ellipses.

        'contours' means that the returned value is a shapely MultiPolgon not a list of patches
        Return
        ------
        list of patches
        """
        from matplotlib.patches import Ellipse
        
        # -- maskout non matched one if requested
        if not self.has_catalogue():
            print "no catalogue mask applied"
            apply_catmask = False
            
        # -------------
        # - Properties
        ells = [Ellipse([x,y],a*scaleup,b*scaleup,
                        t*units.radian.in_units("degree"))
                for x,y,a,b,t in self.get(["x","y","a","b","theta"],mask=mask)]
        if not contours:
            return ells
        from ...utils.shape import patch_to_polygon
        return patch_to_polygon(ells)
    
    def get_ellipse_mask(self,width,height, r=3, apply_catmask=False):
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
    def galaxy_contours(self, scaleup=10):
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
