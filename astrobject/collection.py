#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module contain the collection of astrobject"""
import numpy as np
import warnings
from astropy import coordinates, table, units
import matplotlib.pyplot as mpl

from .baseobject import BaseObject, TargetHandler, CatalogueHandler
from .photometry import get_photopoint, Image
from .instruments import instrument as inst 
from .utils.tools import kwargs_update, load_pkl, dump_pkl, is_arraylike
from .utils.shape import draw_polygon, HAS_SHAPELY

__all__ = ["get_imagecollection"]



def get_imagecollection(images,**kwargs):
    """ Create an Image collection from the given list of images """
    return ImageCollection(images=images,**kwargs)


#######################################
#                                     #
# Base Collection                     #
#                                     #
#######################################
class BaseCollection(BaseObject):
    """
    """
    __nature__ = "Collection"

    PROPERTIES         = ["handler"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    def __init__(self):
        """ """
        self.__build__()
         
        
    # ========================== #
    # = Generic Methods        = #
    # ========================== #
    def remove(self,id):
        """This methods enable to remove the given data from the current instance"""
        self._test_source_(id)
        to_rm = self._handler.pop(id)
        del to_rm

    def rename(self,id,newid,force_it=False):
        """This methods enable to change the id with which the data is accessed.
        Remark: If the 'newid' already exists, an exception is raised except
        if force_it is set to true. In that a case, this will overwrite
        the 'newid' entry.
        """
        # -- Test if it is ok
        self._test_source_(id)
        if newid in self.list_id and not force_it:
            raise ValueError("the newid '%s' already exists. Set force_it to true to overwrite it")
        
        # -- Do the switch
        self._handler[newid] = self._handler.pop(id)

    # ----------------- #
    # - Getter        - #
    # ----------------- #
    def get(self, param, mask=None, softout=False, **kwargs):
        """ get data from the object this contains (get of collected object called) """

        ids = self.list_id if mask is None else np.asarray(self.list_id)[mask]

        def has_no_get(id_, softout=softout):
            if softout:
                warnings.warn("The object labeled %s has no get() method"%id_)
                return []
            
            raise AttributeError("The object labeled %s has no get() method"%id_)
            
        return np.asarray([self._handler[id_].get(param, **kwargs)
                if hasattr(self._handler[id_],"get") else has_no_get(id_)
                for id_ in ids])
            
    
    def get_id(self, containing):
        """ returns all the id containing the 'containing'
        part within their text """
        return [id_ for id_ in self.list_sources if containing in id_]
    # ----------------- #
    # - Setter        - #
    # ----------------- #

    # ========================== #
    # = Internal Tools         = #
    # ========================== #
    def _test_source_(self,source):
        """This methods tests if the given 'id' exists"""
        if source not in self.list_sources:
            raise ValueError("'%s' is not a known source "%source +"\n"\
                             "These are: "+", ".join(self.list_sources))
        return True
    
    # ========================== #
    # = Properties             = #
    # ========================== #
    # ---------------------
    # - Data Handlers
    @property
    def _handler(self):
        """
        This is where the data are recorded. Any data.
        """
        # - Handler are dictionnaries
        if self._properties["handler"] is None:
            self._properties["handler"] = {}
            
        return self._properties["handler"]

    @property
    def nsources(self):
        return len(self._handler.keys())
        
    def has_sources(self):
        return len(self._handler.keys()) > 0
    
    @property
    def list_sources(self):
        """This is the list of the known data id"""
        if not self.has_sources():
            raise AttributeError("no data loaded")
        
        return np.sort( list(self._handler.keys()) ).tolist()

    @property
    def list_id(self):
        """This is the list of the known data id (equivalent to list_sources). Empty list return if not data"""
        return [] if not self.has_sources() else \
          self.list_sources

    def _test_id_(self,id_):
        """ Is the given id known? """
        return self._test_source_(id_)
    
# =============================== #
#                                 #
#  Astro-oriented BaseCollection  #
#                                 #
# =============================== #
class Collection( TargetHandler, BaseCollection ):
    """
    """
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    
    def __build__(self,*args,**kwargs):
        """
        """
        super(Collection,self).__build__(*args,**kwargs)
        self._build_properties = {}


    # ========================== #
    # = Properties             = #
    # ========================== #
    @property
    def ndata(self):
        return self.nsources
        
    def has_data(self):
        return self.has_sources()
    


#######################################
#                                     #
# Image Collection                    #
#                                     #
#######################################
class ImageCollection( Collection, CatalogueHandler ):
    """
    """
    PROPERTIES         = []
    SIDE_PROPERTIES    = ["hostcollection"]
    DERIVED_PROPERTIES = []
    
    def __init__(self, images=None, empty=False, catalogue=None, **kwargs):
        """
        """
        self.__build__()
        if empty:
            return
        
        self.create(images,catalogue=catalogue,**kwargs)

        
    def create(self,images, catalogue=None, imageid=None,**kwargs):
        """
        Loop over the images end add them to the instance. See add_image
        """
        if not is_arraylike(images):
            images = [images]
        # ---------------------------
        # - dealing with catalogues
        if catalogue is None:
            catalogue = True
        elif type(catalogue) is not bool: 
            self.set_catalogue(catalogue)

        if imageid is None:
            imageid = [None]*len(images)
            
        [self.add_image(i_,load_catalogue=catalogue,imageid=id_, **kwargs) for i_,id_ in zip(images, imageid) 
         if i_ is not None]

    def set_target(self,newtarget, set_to_images=True):
        """
        Change (or create) an object associated to the given image.
        This function will test if the object is withing the image
        boundaries (expect if *test_inclusion* is set to False).
        Set newtarget to None to remove the association between this
        object and a target

        set_to_images [bool]       If any images is loaded and this is within the image FoV
                                   shall we set the target to the image ?
        """
        super(ImageCollection,self).set_target(newtarget)
        if set_to_images and self.has_data():
            for id_ in self.list_id:
                if self.images[id_]["image"] is not None:
                    try:
                        self.images[id_]["image"].set_target(self.target)
                    except ValueError:
                        warnings.warn("the new target is not in %s's FoV "%id_)
        
                    
        
    # ========================= #
    # = Images IO             = #
    # ========================= #
    def add_image(self,new_image,load_catalogue=True, imageid = None, **kwargs):
        """
        This method enables to save the given image in the self.images container (dict).
        To save speed and memory the given is saved as it is i.e., if its a filename
        the filename is recorded, if its an astrobject Image (or children) this Image
        is recorded. If its a file, it will be loaded only when needed.
        * Important * if a filename is given in input, this methods *TEST if the file
        corresponds to a known instrument*
        ( see astrobject's instrument.is_known_instrument_file )

        Parameters
        ----------
        new_image: [string or astrobject's Image (or children)]

        - option -

        load_catalogue: [bool]     run the 'download_catalogue' method to
                                   download the corresponding catalogue if
                                   it does not exist yet.

        imageid: [string]
            Provide the name of the image. If None, a default one will be used.

        Return
        ------
        Void
        """
        # ------------------------ #
        # What should be loaded  - #
        # ------------------------ #
        if type(new_image) == str:
            idloaded = self._add_image_file_(new_image,imageid=imageid)
        
        elif Image in new_image.__class__.__mro__:
            idloaded = self._add_astroimage_(new_image,imageid=imageid)
        else:
            # -- Issue
            raise TypeError("the given new_image must be a image-file or an astrobject imnage")
        if load_catalogue:
            self.download_catalogue(id_=idloaded)
            
    def set_catalogue(self, catalogue, force_it=False,
                      fast_setup=False):
        """ attach a catalogue to the current instance.
        you can then access it through 'self.catalogue'.

        The current instance's wcs solution is passed to the calague.

        Parameters
        ----------
        fast_setup: [bool] -optional-
            No Test to automatic wcs association, pure setting
        Returns
        -------
        Void
        """
        super(ImageCollection, self).set_catalogue(catalogue, force_it=force_it,
                                                   fast_setup=fast_setup)
        if self.has_host():
            self.host.set_catalogue(self.catalogue)
            
        #for id_ in self.list_id:
        #    if self.images[id_]["image"] is not None:
        #        self.images[id_]["image"].set_catalogue(self.catalogue)
                    
    def load_images(self,verbose=True, **kwargs):
        """
        """
        [self._load_image_(id_, **kwargs) for id_ in self.list_id
         if self.images[id_]["image"] is None]
        # - done
            
    # ========================== #
    # = Get                    = #
    # ========================== #
    def get_image(self,id,load_catalogue=True,
                  dataslice0=None,dataslice1=None,
                  reload=True,
                  **kwargs):
        """

        **kwargs goes to instrument's init if the image is loaded for the first time
        and to reload_data() if not and dataslicing differ from what have been loaded before.
        
        """
        self._test_id_(id)            
        # -----------------------
        # - image not loaded yet
        if self.images[id]["image"] is None:
            self._load_image_(id,dataslice0=dataslice0,
                                 dataslice1=dataslice1,
                              **kwargs)

        im = self.images[id]["image"]
        # ---------------
        # Target
        if not im.has_target() and self.has_target():
            im.set_target(self.target, test_inclusion=HAS_SHAPELY)
            
        # ---------------
        # Catalogue
        if load_catalogue and self.has_catalogue() and not im.has_catalogue():
            # --- Check if the current catalogue is big enough
            im.set_catalogue(self.catalogue.copy())

        return im

    def get_subcatalogue(self, kind="all", isolated_only=False,
                         stars_only=False,catmag_range=[None,None]):
        """
        Get a copy of a catalogue. If could be the entire catalogue (kind=None or kind="all")
        or the shared fov catalogue (kind="shared" or "union") or the combined fov catalogue
        (kind="combined" or "intersection")
        """
        if not self.has_catalogue():
            raise AttributeError("No Catalogue loaded")
        # ---------------
        # - Which mask
        if kind is None or kind.lower() in ["all"]:
            contours = None
        elif kind.lower() in ["shared","union"]:
            contours = self.contours_shared
        elif kind.lower() in ["combined","intersection"]:
            contours = self.contours_combined
        else:
            raise ValueError("I could not parse the given kind %s (must be 'all', 'shared' or 'combined')"%kind)
        # ---------------
        # - get the new catalogue
        return self.catalogue.get_subcatalogue(contours=contours,
                                               isolated_only=isolated_only,
                                               stars_only=stars_only,
                                               catmag_range=catmag_range)

    def get_pixels_slicing(self,central_coords, radius, runits="arcmin", id_=None,
                           wcs_coords=True, clean=True):
        """
        The dataslice0,dataslice1 for the given 'central_coords' with a +/- 'radius' 'runits' value.
        This slicing will be cleaned (pixels within the image FoV) except if clean is False
        """
        # ---------------
        # - one id_
        if id_ is not None:
            self._test_id_(id_)
            radec = np.asarray(central_coords) if not wcs_coords else\
              np.asarray(self.images[id_]["wcs"].world2pix(*central_coords))
            pixels   = radius*self.images[id_]["wcs"].units_to_pixels(runits).value
            dataslice0,dataslice1 = [[radec_-pixels,radec_+pixels] for radec_ in radec[::-1]]
            if clean:
                # -- cleaning
                max0 = self.images[id_]["wcs"]._naxis1
                max1 = self.images[id_]["wcs"]._naxis2
                if dataslice0[0] <0    : dataslice0[0] = 0
                if dataslice1[0] <0    : dataslice1[0] = 0
                if dataslice0[1] >max0 : dataslice0[1] = max0
                if dataslice1[1] >max1 : dataslice1[1] = max1
            return {"dataslice0":np.asarray(dataslice0,dtype=int).tolist(),
                    "dataslice1":np.asarray(dataslice1,dtype=int).tolist()}
        
        # ---------------
        # - get all id_
        if id_ is None:
            out = {}
            for id_ in self.list_id:
                out[id_] = self.get_pixels_slicing(central_coords, radius, runits="arcmin", id_=id_)
            return out
        
    # --------------------- #
    # - Get Extraction    - #
    # --------------------- #
            
          
    # --------------------- #
    # - Get Target        - #
    # --------------------- #
    def get_target_ids(self,target=None):
        """
        This methods enable to access the list of 'data' ids that contains the given target.

        Parameters:
        ----------
        - options -
        target: [None/ AstroTarget]   If no target is given and a target has been set
                                      (set_target) the current target will be used.
                                      If you give here a target, it will use it *without*
                                      changing the current target.
                                      -> Error
                                      This method raise a ValueError is no target accessible.

        Return:
        -------
        list (ids)
        """
        # -----------------
        # - Target
        if target is None:
            if not self.has_target():
                raise ValueError("no target set to the instance and no input target")
            target = self.target
        elif "__nature__" not in dir(target) or target.__nature__ != "AstroTarget":
            raise TypeError("The given 'target' must be an astrobject AstroTarget")
        # -----------------
        # - ID information
        return [id_ for id_ in self.list_id
            if self.images[id_]["wcs"].coordsAreInImage(target.ra,target.dec)]
            
    def get_target_collection(self,target=None):
        """
        This modules enables to get a new ImageCollection containing only
        the images corresponding to the given target (target in the image FoV).

        Parameters:
        -----------

        - options -

        target: [None/ AstroTarget]   If no target is given and a target has been set
                                      (set_target) the current target will be used.
                                      If you give here a target, it will use it *without*
                                      changing the current target.
                                      -> Error
                                      This method raise a ValueError is no target accessible.
                                      (see self.get_target_ids)

        Return:
        -------
        ImageCollection
        """
        im = ImageCollection([self.images[id_]["file"]
                              for id_ in self.get_target_ids(target)])
        
        if self.has_catalogue():
            im.set_catalogue(self.catalogue)
            
        im.set_target(target)
        return im

    def get_target_photopoints(self,radius,runits="kpc",target=None,
                               onflight=False,**kwargs):
        """
        This modules enables to get a new ImageCollection containing only
        the images corresponding to the given target (target in the image FoV).

        Parameters:
        -----------

        radius, runits: [float, string]           
            The aperture photometry radius and its associated unit.

        target: [None/ AstroTarget] -optional-
            If no target is given and a target has been set
            (set_target) the current target will be used.
            If you give here a target, it will use it *without*
            changing the current target.
            -> Error
            This method raise a ValueError is no target accessible.
            (see self.get_target_ids)

        onflight: [bool] -optional-
            Use this to avoid to load the images in the current instance
            If True, this will instead load a copy of the current instance
            using 'get_target_collection' which will then be delaited.
            -> Use this to save cach-memory once the method is finished.
        
        ** kwargs goes to get_target_photopoints (-> get_aperture), including notably:
               aptype, subpix, ellipse_args, annular_args

           
        Return:
        -------
        PhotoPointCollection
        """
        images = self.get_target_images(target=target, onflight=onflight, reload=False)
        
        photopoints = PhotoPointCollection(photopoints=[image_.get_target_photopoint(radius=radius,runits=runits,
                                            **kwargs) for image_ in images], ids=self.list_id)
        photopoints.set_target(self.target)
        # ------------------------------
        # - Shall we delete the files ?
        if onflight:
            del images
            del target_imcoll
            
        return photopoints

    def get_target_images(self, target=None, onflight=False, **kwargs):
        """

        Parameters
        ----------
        target: [None/ AstroTarget]   If no target is given and a target has been set
                                      (set_target) the current target will be used.
                                      If you give here a target, it will use it *without*
                                      changing the current target.
                                      -> Error
                                      This method raise a ValueError is no target accessible.
                                      (see self.get_target_ids)

        onflight: [bool]              Use this to avoid to load the images in the current instance
                                      If True, this will instead load a copy of the current instance
                                      using 'get_target_collection' which will then be delaited.
                                      -> Use this to save cach-memory once the method is finished.

        Returns
        -------
        list of images
        """
        if not onflight:
            ids = self.get_target_ids(target=target) if HAS_SHAPELY else self.list_id
            images = [self.get_image(_id,load_catalogue=False, **kwargs)
                      for _id in ids]
        else:
            target_imcoll = self.get_target_collection(target=target)
            images = [target_imcoll.get_image(_id,load_catalogue=False)
                      for _id in target_imcoll.list_id]
            del target_imcoll
        return images

    def get_host_photopoints(self, scaleup=3, reference_id=None, catid=None,
                             target=None, prop_hostsearch={},
                             verbose=False, force_ref_ellipse=False,
                             **kwargs):
        """
        This modules enables to get a new ImageCollection containing only
        the images corresponding to the given target (target in the image FoV).

        Parameters:
        -----------

        target: [None/ AstroTarget]   If no target is given and a target has been set
                                      (set_target) the current target will be used.
                                      If you give here a target, it will use it *without*
                                      changing the current target.
                                      -> Error
                                      This method raise a ValueError is no target accessible.
                                      (see self.get_target_ids)

        catid: [float/string/None] -optional-
            You can force the id for the host by providing it here.
            If None catid will be ignore and this will run the host-search


        idx_hostprop: [idx / None] -optional-
        
              TO BE DONE
        
        prop_hostsearch: [dict] -optional-

              TO BE DONE

                      
        -- other options --
        
        kwargs goes to get_thost_photopoints (-> get_aperture), including notably:
               aptype, subpix, ellipse_args, annular_args

           
        Return:
        -------
        PhotoPointCollection
        """
        if target is not None:
            self.set_target(target)

        if catid is not None:
            self.host.set_catid(catid)
            
        if reference_id is not None:
            self.host.set_reference_image(reference_id,
                                          fetch_hostid= self.host.host_catid is None and self.has_catalogue(),
                                          **prop_hostsearch)


        return self.host.get_photopoints(scaleup=scaleup, verbose=verbose,
                                         force_ref_ellipse=force_ref_ellipse,
                                          **kwargs)
    
    # ========================== #
    # = Set                    = #
    # ========================== #
    # --------------------- #
    # - Catalogue Methods - #
    # --------------------- #
    def download_catalogue(self,id_=None,
                           source="sdss",
                           **kwargs):
        """
        """
        loop_id = self.list_id if id_ is None else [id_] if not is_arraylike(id_) else id_
        for id_ in loop_id:
            # -- Loop over the images, check what is needed
            if not self.has_catalogue() or \
              (HAS_SHAPELY and not self.catalogue.contours.contains(self.images[id_]["wcs"].contours)):
                print("Fetching a new catalogue")
                new_cat = self._get_id_catalogue_(id_,source=source,**kwargs)
                
                if not self.has_catalogue():
                    self.set_catalogue(new_cat)
                else:
                    self.catalogue.merge(new_cat)
            if self.images[id_]["image"] is not None and not self.images[id_]["image"].has_catalogue():
                self.images[id_]["image"].set_catalogue(self.catalogue)
                
    # ========================== #
    # = Show                   = #
    # ========================== #
    def show_stamps(self, savefile=None, show=True,
                    bands=None, maxcol=3, zoom=50, zunit="arcsec",
                    localcircle=[None, "kpc"],
                    add_refellipse=True, scaleup=3,
                    ellprop={}, **kwargs):
        """ """
        from matplotlib.patches import Ellipse
        from .utils.mpladdon    import figout
        
        if bands is None:
            bands = self.list_id
            
        if len(bands) < maxcol:
            maxcol == len(bands)

        if add_refellipse and not self.host.has_refid():
            warnings.warn("No Reference ID defined. No Reference ellipse available")
            add_refellipse = False
            
        nrow = ( len(bands)-1 ) // maxcol + 1
        fig = mpl.figure(figsize=[10,3*nrow])

        prop = kwargs_update( dict(facecolor="None", edgecolor="k", lw=2),
                                     **ellprop)
        for i,b_ in enumerate(bands):
            ax = fig.add_subplot(nrow, maxcol, i+1)
            im = self.get_image(b_)
            im.show(ax=ax, zoomon="target", zoom=zoom*im.units_to_pixels(zunit),
                    show=False, localcircle=localcircle, proptarget={"circleprop":ellprop},
                    **kwargs)
            
            if add_refellipse:
                x,y,a,b,theta = self.host.get_ref_ellipse_in(b_)
                # 2* since Ellipse in Matplotlib want diameter
                ell = Ellipse([x,y], 2*a*scaleup, 2*b*scaleup,
                                theta*units.radian.in_units("degree"),
                                **prop)
                ell.set_clip_box(ax.bbox)
                ax.add_patch(ell)
            
            # - Highlight reference
            if self.has_host() and b_ == self.host.refid:
                [s_.set_linewidth(2) for s_ in ax.spines.values()]
                    
            # - Cleaning
            ax.text(0.95,0.95, b_, va="top", ha="right", transform=ax.transAxes,
                            bbox=dict(facecolor=mpl.cm.binary(0.1,0.8), edgecolor="k"))
            ax.set_xticks([])
            ax.set_yticks([])
        
        fig.figout(savefile=savefile, show=show)
        
    def show_skypatches(self,ax=None,
                        savefile=None,show=True,
                        fc="None",ec="k",
                        targetprop={},
                        show_catalogue=True,show_target=True,
                        **kwargs):
        """
        Plot On the sky the 
        """
        from .utils.mpladdon import wcsplot,figout,mpl
        
        # ----------------------
        # -   Axis Init
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure

        # ----------------------
        # -  loop over the images
        _prop = kwargs_update({"zorder":5},**kwargs)
        pl = [ax.wcsplot(self.images[id_]["wcs"],fc=fc,ec=ec,**_prop)
              for id_ in self.list_id if self.images[id_]["wcs"] is not None]

        # ----------------------
        # - Catalogue
        if self.has_catalogue() and show_catalogue:
            self.catalogue.display(ax,zorder=3)
            draw_polygon(ax,self.contours_shared,zorder=4,fc=mpl.cm.Greens(0.4,0.3))
            
        # ----------------------
        # - Target
        if self.has_target() and show_target:
            default_markerprop = {
                "marker":"o","mfc":"w",
                "mec":"k","mew":1,"zorder":12,
                "scalex":False,"scaley":False
                }
            prop = kwargs_update(default_markerprop,**targetprop)
            pl = ax.plot(self.target.ra,self.target.dec,**prop)
            
        # ----------------------
        # -  output
        fig.figout(savefile=savefile,show=show)
    
    # ========================== #
    # = INTERNAL               = #
    # ========================== #
    def _get_id_catalogue_(self,id_,source="sdss",radius_degree=None,**kwargs):
        """
        """
        print("downloading catalogue for %s"%id_)
        radec = "%s %s"%(self.images[id_]['wcs'].central_coords[0],
                         self.images[id_]['wcs'].central_coords[1])
        img_radius= self.images[id_]['wcs'].diag_size/ 1.5
        radius = img_radius if radius_degree is None or radius_degree<img_radius \
          else radius_degree
          
        return inst.fetch_catalogue(source=source,radec=radec,radius="%sd"%radius,
                             **kwargs)
    # =============== #
    # = IO images   = #
    # =============== #
    def _add_image_file_(self,imagefile, imageid=None):
        """
        """
        # -- Load a file
        if not inst.is_known_instrument_file(imagefile):
            raise TypeError("the given new_image file is not an image of a known instrument")
        
        ID = imagefile.split("/")[-1] if imageid is None else imageid
            
        self.images[ID] = {
            "file":imagefile,
            "image":None
            }
        try:
            self.images[ID]["wcs"] = inst.get_instrument_wcs(imagefile)
        except:
            warnings.warn("WARNING Failure to load the wcs solution")
            self.images[ID]["wcs"] = None
            
        return ID
    
    def _add_astroimage_(self,astroimage, imageid=None):
        """
        """
        # -- Load an image
        if imageid is not None:
            ID = imageid
        elif astroimage.filename is not None:
            ID = astroimage.filename.split("/")[-1]
        else:
            ID = "image%d"%len(self.images.keys())

        self.images[ID] = {
            "file":astroimage.filename,
            "image":astroimage,
            "wcs":astroimage.wcs if astroimage.has_wcs() else None
            }
        return ID
        
    def _load_image_(self,id,set_target=True, set_catalogue=True, **kwargs):
        """
        """
        self.images[id]["image"] = inst.get_instrument(self.images[id]["file"],**kwargs)
        if self.has_catalogue():
            self.images[id]["image"].set_catalogue(self.catalogue.copy())
        if self.has_target():
            self.images[id]["image"].set_target(self.target)
            
    # ========================== #
    # = Properties             = #
    # ========================== #    
    @property
    def images(self):
        """
        """
        return self._handler
    
    # ------------------ #
    # - Host Tricks    - #
    # ------------------ #
    def has_host(self):
        """ Tests if a host has been set"""
        return self._side_properties["hostcollection"] is not None
    
    @property
    def host(self):
        """ HostImageCollection linked to the current instance.
        It requires that this istance has a target set."""
        if self._side_properties["hostcollection"] is None:
            if not self.has_target():
                warnings.warn(" You need to set a target to have access to the host collection")
                return None
            host = HostImageCollection(empty=True)
            host.link_to(self)
            self._side_properties["hostcollection"] = host
            
        return self._side_properties["hostcollection"]
    # ------------------ #
    # - On flight prop - #
    # ------------------ #
    @property
    def bandnames(self):
        """This returns an array of the recorded images' bandname"""
        # -- This might be included in the add/remove image tools
        return [self.images[id_]["image"].bandname
                  if self.images[id_]["image"] is not None \
                  else inst.which_band_is_file(self.images[id_]["file"])
                for id_ in self.list_id]

    # ------------------
    # - Shapely Based
    @property
    def contours_shared(self):
        """ returns the contour (shapely polygon) of the shared covered sky area"""
        for i,d in enumerate(self.images.values()):
            fov = d["wcs"].contours if i == 0 else\
              fov.intersection(d["wcs"].contours)
        return fov
    
    @property
    def contours_combined(self):
        """ returns the contour (shapely polygon) of the combined covered sky area"""
        for i,d in enumerate(self.images.values()):
            fov = d["wcs"].contours if i == 0 else\
              fov.union(d["wcs"].contours)
        return fov

#######################################
#                                     #
# Host Associated Image collection    #
#                                     #
#######################################
class HostImageCollection( ImageCollection ):
    """ """
    PROPERTIES         = ["refid"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["ref_ellipsedata","hostcatid","host_catindex"]

    # ================== #
    # = Main Methods   = #
    # ================== #
    def link_to(self, imagecollection):
        """ anchor an Image collection to this object """
        if ImageCollection not in imagecollection.__class__.__mro__:
            raise ValueError("The given imagecollection must be (or inherite from) an astrobject's ImageCollection")
        
        self._properties["handler"]     = imagecollection._handler
        if imagecollection.target is not None:
            self.set_target(imagecollection.target)
            
        if imagecollection.catalogue is not None:
            self.set_catalogue(imagecollection.catalogue)
        
    # ---------- #
    # - SETTER - #
    # ---------- #
    def set_reference_image(self, id_, fetch_hostid=False, extractprop={}, **kwargs):
        """
        **kwargs goes to the fetch_catalogue_id() method
        """
        if id_ not in self.list_id:
            raise ValueError("%s is not a known id. These are: "+", ".join(self.list_id))
        
        self._properties["refid"] = id_
        if not self.refimage.has_sepobjects():
            self.refimage.sep_extract(**extractprop)

        if not self.refimage.has_catalogue() and self.has_catalogue():
            self.refimage.set_catalogue(self.catalogue)
            
        if fetch_hostid:
            self.fetch_catalogue_id(**kwargs)

    def fetch_catalogue_id(self, refid=None,scaleup=3,
                          radius=30, runits="kpc",
                          max_galdist=3, verbose=True):
        """ look inside the reference image for the host ID.
        Once identified, access it throught the property 'host_catid'

        Parameters
        ----------
        refid: [None or any image id] -optional-
            Directly change/set the reference id.
            this will only be made if this is not None

        // host search optional
        
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
        Void (sets self.catalogue)
        """
        if refid is not None:
            self.set_reference_image(refid)
        if not self.has_catalogue():
            raise AttributeError("No Catalogue set. Cannot Fetch a catalogue id")
            
        # refimage will load automatically
        # raise AttributeError if id not set
        if not self.refimage.has_sepobjects():
            self.refimage.sep_extract()
            
        host_idx = self.refimage.sepobjects.get_host_idx(self.target,
                                                         radius=radius, runits=runits,
                                                         scaleup=scaleup,
                                                         max_galdist=max_galdist)
        # = No Host
        if host_idx is None:
            message = "No host identified for %s. Used reference ID %s"%(self.target.name, self.refid)
            if verbose: print(message)
            warnings.warn(message)
            return
        
        # = Fetch catalog index
        catidx  = self.refimage.sepobjects.index_to_catindex([host_idx])
        if len(catidx)>0:
            catidx = [[id_] for id_ in catidx if not self.refimage.catalogue.starmask[id_]][0]
        if len(catidx)==0:
            # - this should not happen
            raise ValueError("No galaxy associated to the target in image band %s"%band_hostprop)

        self._derived_properties["host_catindex"] = catidx
        self._derived_properties["hostcatid"]     = self.refimage.catalogue.idx_to_id(catidx)

    def set_catid(self, value):
        """ provide the catalogue ID of the host """
        if value not in self.catalogue.data[self.catalogue._build_properties["key_id"]]:
            raise ValueError("%s is not a known catalog ID entry"%value)
        
        self._derived_properties["host_catindex"] = self.catalogue.id_to_idx(value)
        self._derived_properties["hostcatid"]     = value if is_arraylike(value) else [value]
        
    # ---------- #
    # - GETTER - #
    # ---------- #
    
    def get_sepid(self, image_id, catid=None, **kwargs):
        """ get the sepobject entry of the host for the given
        image_id.
        You can select a given catid, otherwise,
        this will use the current catid if it exists
        Returns
        -------
        int (sepobject id)
        """
        self._test_id_(image_id)

        # - Basic image loading
        if self.images[image_id]["image"] is None:
            self._load_image_(image_id)
            
        if not self.images[image_id]["image"].has_sepobjects():
            self.images[image_id]["image"].sep_extract()
            
        # - Cat ID provided
        if catid is not None or self.host_catid is not None:
            catidx = self.images[image_id]["image"].catalogue.id_to_idx(self.host_catid
                                                            if self.host_catid is not None else catid)
            return self.images[image_id]["image"].sepobjects.catindex_to_index(catidx) if \
               catidx in self.images[image_id]["image"].sepobjects.catmatch["idx_catalogue"] else\
               None
        else:
            return self.images[image_id]["image"].sepobjects.get_host_idx(self.target, **kwargs)

    # ------------------ #
    #  Global Ellipse    #
    # ------------------ #
    def get_ref_ellipse(self):
        """ Returns the x, y, a, b, theta parameters 
        if sep object associated to the catalogue id. """
        if not self.has_refid():
            raise AttributeError("Reference id not set. run self.set_reference_image()")
        
        return self.get_ellipse(self.refid)[0]

    def get_ref_ellipse_in(self,id_):
        """ get the reference ellipse and projects it from the reference image
        into the given `id_` image
        Returns
        -------
        x, y, a, b, theta
        """
        x, y, a, b, theta = self.get_ref_ellipse()
        return self.refimage.sepobjects.project_ellipse_to_wcs(x, y, a, b, theta,
                                                        self.images[id_]["wcs"])

    def get_ellipse(self, id_, **kwargs):
        """ get the host ellipse contours for the given image id.
        
        Returns
        -------
        [x,y,a,b,theta], bool (False means used default a,b,theta)
        """
        self._test_id_(id_)
        if not self.images[id_]["image"].has_sepobjects():
            self.images[id_]["image"].sep_extract()
            
        idx = self.get_sepid(id_, catid = self.host_catid, **kwargs)
        
        # = has the data
        if idx is None:
            if not self.has_refid():
                raise AttributeError("No reference image has been set. No sep id detected for the given host.")
            catidx = self.images[id_]["image"].catalogue.id_to_idx(self.host_catid)
            x,y = self.images[id_]["image"].coords_to_pixel(*np.asarray(self.images[id_]["image"].catalogue.get(["ra","dec"],
                                                                            mask= catidx)).T[0])
            return np.concatenate([[x],[y],self.ref_host_ellipse[2:]]), False
        
        # = has the data
        return self.images[id_]["image"].sepobjects.get_ellipse_values(idx if is_arraylike(idx) else [idx]).T, True

    def get_petrorad(self, refid=None):
        """ Measures the Petrosian Radius as defined by SDSS 
        based on data from the given refid.
        SDSS uses twice this radius to get the petroMag
        """
        import sep
        from scipy import optimize
        if refid is None:
            data_     = self.refimage.data
            x,y,a,b,t = self.get_ref_ellipse()
        else:
            data_ = self.get_image(refid).data
            x,y,a,b,t = self.get_ellipse(refid)[0]

        def rp_lim_res(rp):
            return np.abs(0.2- sep.sum_ellipann(data_, x,y,a,b,t, rp*0.8,rp*1.2)[0] / sep.sum_ellipse(data_, x,y,a,b,t, rp)[0] )
        
        return optimize.fmin(rp_lim_res,
                             2., disp=0)[0]

    
    def get_photopoints(self, scaleup=3, refid=None,
                        verbose=False, use_cat_mags=False,
                        force_ref_ellipse=False,
                        **kwargs):
        """
        This modules enables to get a new ImageCollection containing only
        the images corresponding to the given target (target in the image FoV).

        Parameters:
        -----------

        use_cat_mags: [bool] -optional-
            Force the use of the calague magnitudes instead of the extracted ones.
            [In dev]
        
                      
        -- other options --
        
        kwargs goes to get_thost_photopoints (-> get_aperture), including notably:
               aptype, subpix, ellipse_args, annular_args

           
        Return:
        -------
        PhotoPointCollection
        """
        # ---------------------
        # - Is instance ready?

        # Reference ID
        if refid is not None:
            self.set_reference_image(refid)
            
        if use_cat_mags and not self.has_catalogue():
            warnings.warn("No catalogue loaded, cannot the catalogue magnitude. use_cat_mags set to False")
            use_cat_mags = False
            
        # ---------------------
        # - Get the PhotoPoints
        # No Host identified
        if self.host_catid is None and self.has_catalogue():
            print("not catid set.")
            photopoints = PhotoPointCollection(empty=True)
            
        # Host identified
        else:
            pps= []
            for id_ in self.list_id:
                #  Our data
                img = self.get_image(id_)
                if not use_cat_mags:
                    x,y,a,b,theta = self.get_ref_ellipse_in(id_) if force_ref_ellipse else self.get_ellipse(id_)[0]
                    pps.append(img.get_photopoint(x, y, radius=scaleup, runits="pixels",
                                                  ellipse_args=dict(a=a, b=b, theta=theta),
                                                  aptype="ellipse", getlist=True,**kwargs))
                #  Catalogue Data
                else:
                    if "sdss" not in img.bandname and self.catalogue.source_name.lower() =="sdss":
                        raise NotImplementedError("Only the catalogue magnitude measurement made for SDSS entries")
                    from astrobject.utils import tools
                    mag_cat = img.bandname[-1]+"mag"
                    mag, mag_err = self.catalogue.data[ self.catalogue.id_to_idx( self.host_catid)][mag_cat,"e_"+mag_cat][0].data
                    flux, err = tools.mag_to_flux(mag, mag_err,img.lbda)
                    pps.append(get_photopoint(lbda=img.lbda, flux=flux, var=err**2,
                                        source="catalogue",mjd=img.mjd,
                                        zp=img.mab0, bandname=img.bandpass.name,
                                        instrument_name=self.catalogue.source_name))
                    
            photopoints = PhotoPointCollection(photopoints=pps, ids=self.list_id)
        
        # - Set the target
        photopoints.set_target(self.target)
        
        return photopoints

    # ------------ #
    # - PLOTTER  - #
    # ------------ #
    def display_ref_ellipse(self, ax, id_, scaleup, **prop):
        """ """
        from matplotlib.patches import Ellipse
        x,y,a,b,theta = self.get_ref_ellipse_in(id_)
        prop = kwargs_update({"facecolor":"None", "edgecolor":"k"},**kwargs_)
        ell= Ellipse([x,y],a*scaleup,b*scaleup,
                    theta*units.radian.in_units("degree"),
                    **prop)

        ell.set_clip_box(ax.bbox)
        ax.add_patch(ell)
        if draw:
            ax.figure.canvas.draw()
        return ell
    
    def show(self, fig=None, savefile=None, show=True,
             scaleup=3, ellipse_prop={}, id_to_show=None,
              **kwargs):
        """
        """
        from matplotlib.patches import Ellipse
        from astrobject.utils.mpladdon import figout
        
        def add_ellipse(ax, host_ellipse, scaleup=scaleup, draw=True,**kwargs_):
            x,y,a,b,theta = host_ellipse
            prop = kwargs_update({"facecolor":"None", "edgecolor":"k"},**kwargs_)
            ell= Ellipse([x,y],a*scaleup*2,b*scaleup*2,
                    theta*units.radian.in_units("degree"),
                    **prop)
            ell.set_clip_box(ax.bbox)
            ax.add_patch(ell)
            if draw:
                ax.figure.canvas.draw()
            return ell

        # -------------
        # - Input
        if id_to_show is None:
            id_to_show = self.list_id
        ndata = len(id_to_show)
        # -------------
        # - Fig & Axes
        if fig is None:
            fig = mpl.figure(figsize=[12,6])

        axes = [fig.add_subplot(1,ndata,i+1) for i in range(ndata)]

        # -------------
        # - Properties
        show_prop = kwargs_update(dict(show_catalogue=True,show_sepobjects=False,\
                                       zoomon="target", zoom=30, zunits="kpc"),
                                       **kwargs)
        
        for id_,ax in zip(id_to_show,axes):
            if self.images[id_]["image"] is None: self._load_image_(id_)
                
            self.images[id_]["image"].show(ax =ax, show=False, **show_prop)
            try:
                ellip_data, notforced = self.get_ellipse(id_)
                add_ellipse(ax, ellip_data, ls= "-" if notforced else "--",
                            **kwargs_update({"lw":2},**ellipse_prop))
            except:
                    warnings.warn("No Ellipse draw possible for %s"%id_)

            ax.text(0.05,0.95, "band: %s"% self.images[id_]["image"].bandname, 
                    fontsize="large",va="top",ha="left",
                    transform=ax.transAxes, bbox=dict(facecolor='w', alpha=0.9))
            
        # = Titles
        if self.has_catalogue():
            fig.suptitle("%s ID of %s-host: %s"%(self.catalogue.source_name, self.target.name, self.host_catid[0] if self.host_catid is not None else "no"), 
                        y=0.90, va="top",fontsize="x-large")
        else:
            fig.suptitle("%s-host"%self.target.name,  y=0.90, va="top",fontsize="x-large")
        
        fig.figout(savefile=savefile, show=show)
        
    # ================== #
    # =  Properties    = #
    # ================== #

    # = Reference
    @property
    def host(self):
        return None
    
    @property
    def refid(self):
        """ The id defined as reference image. """
        return self._properties["refid"]

    def has_refid(self):
        """ test if a reference id has been set."""
        return self.refid is not None
    
    @property
    def ref_ellipse(self):
        """ """
        return self.get_ellipse(self.refid)[0]
    
    @property
    def refimage(self):
        """ """
        if not self.has_refid():
            raise AttributeError("You did not set the refid")
        if self.images[self.refid]["image"] is None:
            self._load_image_(self.refid)
            
        return self.images[self.refid]["image"]

    # ----------
    # - IDs
    # ----------
    # Catalogue
    @property
    def host_catid(self):
        """ """
        return self._derived_properties["hostcatid"]
    
    @property
    def _ref_ellipsedata(self):
        """ """
        if self._derived_properties["ref_ellipsedata"] is None:
            self._derived_properties["ref_ellipsedata"] = {}
                
        return self._derived_properties["ref_ellipsedata"]
    

    

    
#######################################
#                                     #
# Image Collection -> SED SOURCE      #
#                                     #
#######################################
class PhotoPointCollection( Collection ):
    """
    """
    PROPERTIES         = []
    SIDE_PROPERTIES    = ["table"]
    DERIVED_PROPERTIES = ["getkeys","fromtable"]
    
    def __init__(self, photopoints=None,filein=None,empty=False,**kwargs):
        """
        """
        self.__build__()
        if empty:
            return
        
        if filein is not None:
            self.load(filein,**kwargs)
            
        if photopoints is not None:
            self.create(photopoints,**kwargs)
            
    # =============================== #
    # = Main Methods                = #
    # =============================== #

    # ------------------- #
    # - I/O PhotoPoint  - #
    # ------------------- #
    def create(self,photopoints,ids=None):
        """
        """
        if not is_arraylike(photopoints):
            photopoints = [photopoints]
        if not is_arraylike(ids):
            ids = [ids]*len(photopoints)
            
        if len(ids) != len(photopoints):
            print("ids: (%d)"%len(ids))
            print("photopoints: (%d)"%len(photopoints))
            raise ValueError("photopoints and ids must have the same size")
        
        [self.add_photopoint(p_,id_) for p_,id_ in zip(photopoints,ids)]

    def add_photopoint(self,photopoint,id_=None):
        """
        This method enables to register the given photpoint in the self.photopoints
        container (dict).
        
        Parameters
        ----------
        new_image: [string or astrobject's Image (or children)]

        - option -
        id_: [any]                 key used to access the photopoint from the _handler.
                                   id_ is set the photopoint's bandname if None
        
        Return
        ------
        Void
        """
        # --------------
        # - Define ID
        if photopoint is None:
            return
        
        ID = id_ if id_ is not None else photopoint.bandname
        if self.has_data():
            if ID in self.list_id:
                i = 1
                while ID+"-%d"%i in self.list_id:
                    i+=1
                ID = ID+"-%d"%i
        # -------------------
        # - Record the point        
        self.photopoints[ID]                = photopoint
        self._derived_properties["getkeys"] = None
        
    def create_from_table(self, table, idkey=None):
        """ provide a table (tested for astropy table; should work with pandas).
        The columns are the keys value that you can 'get()' ; the rows are the
        individual photopoint data.
        """
        
        if not hasattr(table, "columns"):
            raise TypeError("table is not a Table (no columns attribute)")
        
        self._side_properties["table"]        = table
        # the associated entries
        self._derived_properties["fromtable"] = True
        self._build_properties["idkey"]       = "number"


    def get_photopoint(self, id):
        """ returns a copy of the photometric point"""
        if not self.fromtable:
            return self.photopoints[id].copy()
        
        return get_photopoint(**{k:self._table[self._build_properties["idkey"]==id][k] for k in self._table.keys()})
    
    # ------ #
    #  IO    #
    # ------ #
    def writeto(self,filename,astable=False,**kwargs):
        """
        Save the PhotoPointCollection following the astropy's Table 'write()' method.
        Such data will be accessible using the 'load()' method of the class.
        """
        if astable:
            warnings.warn("Only the data table can be saved in ascii format")
            self.data.write(filename,format=format,**kwargs)
            return
        
        dump_pkl(self._fulldata, filename)
            
    def load(self,filename,**kwargs):
        """ e.g. kwargs you can specify the format (format='...') """
        if filename.endswith("pkl"):
            self._load_pkl_(filename)
        else:
            self._load_table_(filename,**kwargs)

    def _load_pkl_(self, filename=None, get_loaded_data=False, fulldata=None ):
        """ """
        if fulldata is None:
            fulldata = load_pkl(filename)
            
        if "data" not in fulldata.keys():
            raise TypeError("Wrong pkl data format {data:data, ...}")

        if "build" in fulldata.keys():
            self._build_properties = kwargs_update(self._build_properties,**fulldata["build"])
            
        data = self.create_from_table(table.Table(fulldata["data"]),
                                      idkey = None if "idkey" not in self._build_properties.keys() else\
                                              self._build_properties["idkey"])
        
        if "target" in fulldata.keys() and fulldata["target"] is not None:
            from .baseobject import get_target
            self.set_target(get_target(**fulldata["target"]))
            
        
            
        if get_loaded_data:
            return fulldata

        
    # - load for ascii format
    def _load_table_(self, filename, format="ascii",**kwargs):
        data = table.Table.read(filename,format=format)
        self.create_from_table(data)
        if "comments" in data.meta.keys():
            self._read_table_comment_(data.meta["comments"])
            
    def _read_table_comment_(self, comment):
        """ """
        warnings.warn("The following comments on the table are not saved:"+"\n".join(comment))
        pass
        
    # ------------------- #
    # - Getter          - #
    # ------------------- #
    def get_data(self, which="complet", format="dict"):
        """ return a copy of the requested data in the requested format.
        For individual parameter information see get()

        Parameter
        ---------
        which: [string] complet/data
            Which can of data to you want to get.
            Default 'complet' is the combination of data and meta

        format: [string] dict/table
            The default format is the python dict(). You can also access
            the astropy table format.

        Returns
        -------
        dict/table (see format) or None if no data available
        """
        if not self.fromtable:
            data = self.data_complet.copy() if which is "complet" else \
            self.data.copy() 
        else:
            data = self._table.copy()
            
        if data is None:
            return None

        if format == "table":
            return data
        if format == "dict":
            dt = {}
            dt_ = dict(data)
            for k in dt_.keys():
                dt[k] = dt_[k].data
            return dt
        raise ValueError("unknown format.")
      
    def get(self,param, mask=None, safeexit=True):
        """Loop over the photopoints (following the list_id sorting) to get the
         parameters using the photopoint.get() method.

        If safeexit is True NaN will be set to cases where the get(key) does not exist.
        Otherwise this raises a ValueError.
        
        Returns
        -------
        list
        """
        if not self.fromtable:
            return super(PhotoPointCollection, self).get(param, mask=mask, safeexit=safeexit)
        # not np.asarray(self._table[mask, param]) to have pure float
        if mask is not None:
            if is_arraylike(param):
                return np.asarray([self._table[mask][p] for p in param]).T
            return np.asarray(self._table[mask][param])
        else:
            if is_arraylike(param):
                return np.asarray([self._table[p] for p in param]).T
            return np.asarray(self._table[param])

    # ------------------- #
    # - Setter          - #
    # ------------------- #
    def set_meta(self, key, values, ids=None):
        """
        Assign to all photopoints of the instance a 'value' registered
        as 'key' in their meta dict. The 'values' could either be a
        constant (all photopoint will then have the same) or an N-array
        where N is the number of photopoints.
        The set meta 'values' would then be accessible using the 'get()' method.

        Returns:
        --------
        Void
        """
        #if ids is not None:
            # -- This below is really time consuming. don't !
            #[self._test_id_(id_) for id_ in ids]
        id_to_loop = self.list_id if ids is None else ids
        nids = len(id_to_loop)
        if not is_arraylike(values):
            values = [values]*nids
        elif len(values)==1:
            values = [values[0]]*nids
        elif len(values) !=nids:
            raise ValueError("the given 'values', must provide a value for each photopoints (%d)"\
                             +" or set ids to the correct id list"\
                             +" or a unique value  %d given"%(self.nsources,len(values)))
        if self.fromtable:
            if key in self._table.colnames:
                self._table.remove_column(key)
            self._table.add_column(table.Column(values,key))
        else:
            for id_,v in zip(id_to_loop,values):
                self.photopoints[id_].meta[key] = v 
            
            
    # ------------------- #
    # - Plot Stuff      - #
    # ------------------- #
    def show(self,mode="sed",savefile=None,toshow="flux",ax=None,
             cmap="jet",show=True,**kwargs):
        """
        """
        if mode.lower()=="sed":
            kwargs["function_of_time"] = False
            return self._show_(savefile=savefile,toshow=toshow,ax=ax,
                        cmap=cmap,show=show,**kwargs)
        if mode.lower()=="lightcurve":
            kwargs["function_of_time"] = True
            return self._show_(savefile=savefile,toshow=toshow,ax=ax,
                        cmap=cmap,show=show,**kwargs)

    def show_hist(self,toshow="a",ax=None,mask=None,
                  savefile=None,show=True,
                  **kwargs):
        """This methods enable to show the histogram of any given
        key."""
        # -- Properties -- #
        
        v = self.get(toshow,mask=mask)

        # -- Setting -- #
        from .utils.mpladdon import figout
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


    # ========================== #
    # = Internal Stufff        = #
    # ========================== #
    def _show_(self,savefile=None,toshow="flux",ax=None,
                cmap="jet",show=True,**kwargs):
        """
        """
        from .utils.mpladdon import figout
        self._plot = {}
        # --------------------
        # --- INPUTS
        # - axes
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            ax.set_xlabel(r"Wavelength [$\AA$]" if not kwargs["function_of_time"] else "Time [MJD]",
                              fontsize="large")
            ax.set_ylabel(r"Flux [per $\lambda$]", fontsize="large")
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
            
        # - Colors
        if type(cmap) is str:
            cmap = eval("mpl.cm.%s"%cmap)
        colors = kwargs.pop("color",[cmap((i-np.min(self.lbdas))/(np.max(self.lbdas)-np.min(self.lbdas))) for i in self.lbdas])
        prop   = kwargs_update({"ms":10},**kwargs)
        # -------------
        # - Da Plot
        pl = [self.photopoints[id_].display(ax,toshow=toshow,color=c_,**prop)["plot"]
              for id_,c_ in zip(self.list_id,colors)]

        
        # -------------
        # - Out
        self._plot['ax'] = ax
        self._plot['fig'] = fig
        self._plot['plot'] = pl
        self._plot['prop'] = kwargs
        fig.figout(savefile=savefile,show=show)

    
    # ========================== #
    # = Properties             = #
    # ========================== #    
    @property
    def photopoints(self):
        """
        """
        return self._handler

    @property
    def fromtable(self):
        """ were the data loaded from a table """
        if self._derived_properties["fromtable"] is None:
            self._derived_properties["fromtable"] = False
        return self._derived_properties["fromtable"]

    @property
    def _table(self):
        """ The table containing the data, if fromtable"""
        return self._side_properties["table"]

    # ------------------ #
    # - On flight prop - #
    # ------------------ #
    # -- Bands
    @property
    def bandnames(self):
        """This returns an array of the recorded images' bandname (self.bandnames<==>self.get('bandname'))"""
        # -- This might be included in the add/remove image tools
        return self.get("bandname")

    @property
    def nbands(self):
        """The amount of diffenrent recorded bands"""
        return np.unique(self.bandnames)

    # -- Flux Variances etc.
    @property
    def lbdas(self):
        """The wavelength of every photopoints (self.lbdas<==>self.get(lbda'))"""
        return self.get("lbda")

    @property
    def fluxes(self):
        """The fluxes of every photopoints (self.fluxs<==>self.get('flux'))"""
        return self.get("flux")
    
    @property
    def fluxvars(self):
        """The flux variances of every photopoints (self.fluxvars<==>self.get('var'))"""
        return self.get("var")
    
    @property
    def mags(self):
        """The (AB) mag of every photopoints (self.mags<==>self.get('mag'))"""
        return self.get("mag")

    @property
    def magvars(self):
        """The (AB) mag variances of every photopoints (self.magvars<==>self.get('magvar'))"""
        return self.get("magvar")
    
    # -- Times
    @property
    def mjds(self):
        """The modified Julian Dates of every photopoints (self.mdjs<==>self.get('mjd'))"""
        return self.get("mjd")

    
    # --------------
    # - Derived
    @property
    def data(self):
        """ builds an astropy table containing the fundatemental parameters of the collection """
        if not self.fromtable:
            # - Main
            maindata = [self.list_id,self.fluxes,self.fluxvars,self.lbdas,self.mjds,self.bandnames,
                        self.get("zp"),self.get("zpsys")]
            mainnames= ["id","flux","var","lbda","mjd","bandname","zp","zpsys"]
            # - Meta
            return table.Table(data=maindata,
                               names=mainnames)
        return self._table
    
    @property
    def meta(self):
        """  builds an astropy table containing the meta information about the 
        sources of the collection """
        raise NotImplementedError("self.meta does not exist in PhotoPointCollections anymore")
    @property
    def getkeys(self):
        """ List of keys known in data and meta you can call by get().
        Remark that these are not the only thing you can call with get.
        """
        if self._derived_properties["getkeys"] is None:
            if self.fromtable:
                self._derived_properties["getkeys"] = list(self._table.keys())
            else:
                self._derived_properties["getkeys"] = np.sort(np.unique( list(self.data.keys())))
                
        return self._derived_properties["getkeys"]
    
    @property
    def _fulldata(self):
        """ Combination of the data and the metadata """
        
        if self.fromtable:
            data = self._table
        else:
            data = self.data
        return {"data": np.asarray(data),
                "target": self.target.data if self.has_target() else None,
                "build": self._build_properties}
        



class TargetPhotoPointCollection( PhotoPointCollection ):
    """
    PhotoPoint Collection associated to a single target.
    (like lightcurves or Photo SED)
    """
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["mw_corrected"]

    # --------------------- #
    # Correct Extinction    #
    # --------------------- #
    def correct_extinction(self, embv, r_v=3.1, law="fitzpatrick99"):
        """ """
        [self.photopoints[b].apply_extinction(embv,r_v=r_v,law=law) for b in self.list_id]
        
    def correct_mw_extinction(self, r_v=3.1, law="fitzpatrick99"):
        """ Correct for the MW extinction"""
        if self.mw_corrected:
            warnings.warn("MW correction already applied. Skiped")
            return
        
        if not self.has_target():
            raise AttributeError("No target attached to this TargetPhotoPointCollection")

        mwebmv = self.target.mwebmv
        if mwebmv is None:
            raise ValueError("No MW extinction for the given target (self.target.mwebmv is None)")
        
        self.correct_extinction(mwebmv,r_v=r_v,law=law)
        self._derived_properties["mw_corrected"] = True

    # ================ #
    #  Properties      #
    # ================ #
    @property
    def mw_corrected(self):
        """ Has the Milky Way extinction correction ran already ?
        CAUTION it does not mean that the MW extinction is not already corrected.
        For Instance, the provided photopoints could already be MW corrected
        """
        if self._derived_properties["mw_corrected"] is None:
            self._derived_properties["mw_corrected"] = False
            
        return self._derived_properties["mw_corrected"]
