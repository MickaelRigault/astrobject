#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module contain the collection of astrobject"""
import numpy as np
import warnings
from astropy import coordinates, table, units
import matplotlib.pyplot as mpl

from .baseobject import BaseObject, TargetHandler
from .photometry import get_photopoint
from .instruments import instrument as inst
from .utils.tools import kwargs_update
from .utils.shape import draw_polygon, HAS_SHAPELY

__all__ = ["ImageCollection"]



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
        
        return np.sort(self._handler.keys()).tolist()
    

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

    def _test_id_(self,id_):
        """
        """
        return self._test_source_(id_)
    # ========================== #
    # = Properties             = #
    # ========================== #
    @property
    def ndata(self):
        return self.nsources
        
    def has_data(self):
        return self.has_sources()
    
    @property
    def list_id(self):
        """This is the list of the known data id (equivalent to list_sources). Empty list return if not data"""
        return [] if not self.has_sources() else \
          self.list_sources 

#######################################
#                                     #
# Image Collection                    #
#                                     #
#######################################
class ImageCollection( Collection ):
    """
    """
    PROPERTIES         = []
    SIDE_PROPERTIES    = ["catalogue", "hostcollection"]
    DERIVED_PROPERTIES = []
    
    def __init__(self,images=None,empty=False,catalogue=None):
        """
        """
        self.__build__()
        if empty:
            return
        
        self.create(images,catalogue=catalogue)

        
    def create(self,images,catalogue=None):
        """
        Loop over the images end add them to the instance. See add_image
        """
        if "__iter__" not in dir(images):
            images = [images]
        # ---------------------------
        # - dealing with catalogues
        if catalogue is None:
            catalogue = True
        elif type(catalogue) is not bool: 
            self.set_catalogue(catalogue)

        [self.add_image(i_,load_catalogue=catalogue) for i_ in images
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
                        warnings.wanr("the new target is not in %s's FoV "%id_)
        
                    
        
    # ========================= #
    # = Images IO             = #
    # ========================= #
    def add_image(self,new_image,load_catalogue=True):
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
        Return
        ------
        Void
        """
        # ------------------------ #
        # What should be loaded  - #
        # ------------------------ #
        if type(new_image) == str:
            idloaded = self._add_image_file_(new_image)
        
        elif "__nature__" in dir(new_image) and new_image.__nature__ != "Image":
            idloaded = self._add_astroimage_(new_image)
        else:
            # -- Issue
            raise TypeError("the given new_image must be a image-file or an astrobject imnage")
        if load_catalogue:
            self.download_catalogue(id_=idloaded)
            

    def load_images(self,verbose=True):
        """
        """
        [self._load_image_(id_) for id_ in self.list_id
         if self.images[id_]["image"] is None]
        # - done
            
    # ========================== #
    # = Get                    = #
    # ========================== #
    def get_image(self,id,load_catalogue=True,
                  dataslice0=[0,-1],dataslice1=[0,-1],
                  **kwargs):
        """

        **kwargs goes to instrument's init if the image is loaded for the first time
        and to reload_data() if not and dataslicing differ from what have been loaded before.
        
        """
        self._test_id_(id)
        # --------------------
        # - image already exist
        if self.images[id]["image"] is None:
            self._load_image_(id,dataslice0=dataslice0,
                                 dataslice1=dataslice1,
                              **kwargs)
            return self.get_image(id,dataslice0=dataslice0,
                                    dataslice1=dataslice1)

        im = self.images[id]["image"]
        # ------------------
        # - Check if slicing ok
        if dataslice0 != im._build_properties["dataslice0"] or \
          dataslice1 != im._build_properties["dataslice1"]:
            # -- Remark: This update the self.images[id]["image"]
            im.reload_data(dataslice0, dataslice1,**kwargs)
            
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

        radius_kpc: [float]           The aperture photometry radius in kpc.
                                      (kpc will be converted in pixels for every images)
                                      
                                      
        - options -

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
        
        -- other options --
        
        kwargs goes to get_target_photopoints (-> get_aperture), including notably:
               aptype, subpix, ellipse_args, annular_args

           
        Return:
        -------
        PhotoPointCollection
        """
        images = self.get_target_images(target=target, onflight=onflight)
        
        photopoints = PhotoPointCollection(photopoints=[image_.get_target_photopoint(radius=radius,runits=runits,
                                            **kwargs) for image_ in images])
        photopoints.set_target(self.target)
        # ------------------------------
        # - Shall we delete the files ?
        if onflight:
            del images
            del target_imcoll
            
        return photopoints

    def get_target_images(self, target=None, onflight=False):
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
            images = [self.get_image(_id,load_catalogue=False)
                      for _id in ids]
        else:
            target_imcoll = self.get_target_collection(target=target)
            images = [target_imcoll.get_image(_id,load_catalogue=False)
                      for _id in target_imcoll.list_id]
            del target_imcoll
        return images

    def get_host_photopoints(self, scaleup=2.5, reference_id=None,
                             target=None, prop_hostsearch={},
                             verbose=False,**kwargs):
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
            
        if reference_id is not None:
            self.host.set_reference_image(reference_id, fetch_hostid=True, **prop_hostsearch)
            
        return self.host.get_photopoints(scaleup=scaleup, verbose=verbose, **kwargs)
    
    # ========================== #
    # = Set                    = #
    # ========================== #
    # --------------------- #
    # - Catalogue Methods - #
    # --------------------- #
    def set_catalogue(self,catalogue,force_it=False):
        """
        """
        if self.has_catalogue() and force_it is False:
            raise AttributeError("'catalogue' already defined"+\
                    " Set force_it to True if you really known what you are doing")

        if "__nature__" not in dir(catalogue) or catalogue.__nature__ != "Catalogue":
            raise TypeError("the input 'catalogue' must be an astrobject catalogue")
        # Important:
        # No WCS solution loaded here, this will be made image by image
        # --------
        # - set it
        self._side_properties["catalogue"] = catalogue
        

    def download_catalogue(self,id_=None,
                           source="sdss",
                           **kwargs):
        """
        """
        loop_id = self.list_id if id_ is None else [id_] if "__iter__" not in dir(id_) else id_
        for id_ in loop_id:
            # -- Loop over the images, check what is needed
            if not self.has_catalogue() or \
              not self.catalogue.contours.contains(self.images[id_]["wcs"].contours):
                print "Fetching a new catalogue"
                new_cat = self._get_id_catalogue_(id_,source=source,**kwargs)
                
                if not self.has_catalogue():
                    self.set_catalogue(new_cat)
                else:
                    self.catalogue.merge(new_cat)
    
    # ========================== #
    # = Show                   = #
    # ========================== #
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
        print "downloading catalogue for %s"%id_
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
    def _add_image_file_(self,imagefile):
        """
        """
        # -- Load a file
        if not inst.is_known_instrument_file(imagefile):
            raise TypeError("the given new_image file is not an image of a known instrument")

        ID = imagefile.split("/")[-1]
        self.images[ID] = {
            "file":imagefile,
            "image":None
            }
        try:
            self.images[ID]["wcs"] = inst.get_instrument_wcs(imagefile)
        except:
            print "WARNING Failure to load the wcs solution"
            self.images[ID]["wcs"] = None
            
        return ID
    
    def _add_astroimage_(self,astroimage):
        """
        """
        # -- Load an image
        if new_image.filename is not None:
            ID = new_image.filename.split("/")[-1]
        else:
            ID = "image%d"%len(self.images.keys())

        self.images[ID] = {
            "file":new_image.filename,
            "image":new_image,
            "wcs":new_image.wcs if new_image.has_wcs() else None
            }
        return ID
        
    def _load_image_(self,id,set_target=True, set_catalogue=True, **kwargs):
        """
        """
        self.images[id]["image"] = inst.get_instrument(self.images[id]["file"],**kwargs)
        if self.has_target():
            self.images[id]["image"].set_target(self.target)
        if self.has_catalogue():
            self.images[id]["image"].set_catalogue(self.catalogue.copy())
            
    # ========================== #
    # = Properties             = #
    # ========================== #    
    @property
    def images(self):
        """
        """
        return self._handler

    @property
    def _imageids(self):
        print "_imageids to be changed to self.list_id"
        return self.list_id
    
    # -- Catalogue
    @property
    def catalogue(self):
        return self._side_properties['catalogue']

    def has_catalogue(self):
        return not self.catalogue is None 

    # ------------------ #
    # - Host Tricks    - #
    # ------------------ #
    @property
    def host(self):
        """ HostImageCollection linked to the current instance.
        It requires that this istance has a target set."""
        if self._side_properties["hostcollection"] is None:
            if not self.has_target():
                raise AttributeError(" You need to set a target to have access to the host collection")
            if not self.has_catalogue():
                self.download_catalogue()
                
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
    def set_reference_image(self, id_, fetch_hostid=False, **kwargs):
        """
        **kwargs goes to the fetch_catalogue_id() method
        """
        if id_ not in self.list_id:
            raise ValueError("%s is not a known id. These are: "+", ".join(self.list_id))
        
        self._properties["refid"] = id_
        if fetch_hostid:
            self.fetch_catalogue_id(**kwargs)

    def fetch_catalogue_id(self, refid=None,
                          radius=30, runits="kpc",
                          max_galdist=2.5, verbose=True):
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
        """
        if refid is not None:
            self.set_reference_image(refid)
        
        # refimage will load automatically
        # raise AttributeError if id not set
        if not self.refimage.has_sepobjects():
            self.refimage.sep_extract()
            
        host_idx = self.refimage.sepobjects.get_host_idx(self.target,
                                                         radius=radius, runits=runits,
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
        
    # ---------- #
    # - GETTER - #
    # ---------- #
    def get_host_ellipse(self, id_):
        """
        Returns
        -------
        [x,y,a,b,theta], bool (False means used default a,b,theta)
        """
        self._test_id_(id_)

        if not self.has_refid():
            raise AttributeError("No reference image has been set. Do so.")
        if self.host_catid is None:
            raise AttributeError("No host has been found in the reference image.")
        
        
        if self.images[id_]["image"] is None:
            self._load_image_(id_)
        if not self.images[id_]["image"].has_sepobjects():
            self.images[id_]["image"].sep_extract()
            
        
        catidx = self.images[id_]["image"].catalogue.id_to_idx(self.host_catid)
        # = has the data
        if catidx not in self.images[id_]["image"].sepobjects.catmatch["idx_catalogue"]:
            x,y = self.images[id_]["image"].coords_to_pixel(*np.asarray(self.images[id_]["image"].catalogue.get(["ra","dec"],
                                                                            mask= catidx)).T[0]
                                            )
            return np.concatenate([[x],[y],self.ref_host_ellipse[2:]]), False
        
        # = has the data
        idx = self.images[id_]["image"].sepobjects.catindex_to_index(catidx)
        return self.images[id_]["image"].sepobjects.get(["x","y","a","b","theta"], mask=idx)[0], True
    
    def get_photopoints(self, scaleup=2.5,refid=None,
                        verbose=False,**kwargs):
        """
        This modules enables to get a new ImageCollection containing only
        the images corresponding to the given target (target in the image FoV).

        Parameters:
        -----------


        onflight: [bool]              Use this to avoid to load the images in the current instance
                                      If True, this will instead load a copy of the current instance
                                      using 'get_target_collection' which will then be delaited.
                                      -> Use this to save cach-memory once the method is finished.
                      
        -- other options --
        
        kwargs goes to get_thost_photopoints (-> get_aperture), including notably:
               aptype, subpix, ellipse_args, annular_args

           
        Return:
        -------
        PhotoPointCollection
        """
        # ---------------------
        # - Is instance ready?
        # Catalogue
        if not self.has_catalogue():
            self.download_catalogue()
        # Reference ID
        if refid is not None:
            self.set_reference_image(refid)
        if not self.has_refid():
            raise AttributeError("You need to set (or provide) a reference image")

        
        # ---------------------
        # - Get the PhotoPoints
        # No Host identified
        if self.host_catid is None:
            photopoints = PhotoPointCollection(empty=True)
        # Host identified
        else:
            pps= []
            for id_ in self.list_id:
                [x,y,a,b,theta], notforced = self.get_host_ellipse(id_)
                if verbose:
                    print "x,y position", x,y
                    print "a,b,theta, ",a,b,theta
                    print "band, ", self.images[id_]["image"].bandname
                    print "ellipse param from reference image", not notforced
    
                pps.append(self.images[id_]["image"].get_photopoint(x, y, radius=scaleup, runits="pixels",
                                                     ellipse_args=dict(a=a, b=b, theta=theta),
                                                     aptype="ellipse", getlist=True,**kwargs))
                
            photopoints = PhotoPointCollection(photopoints=pps)
        
        # - Set the target
        photopoints.set_target(self.target)
        
        return photopoints

    # ------------ #
    # - PLOTTER  - #
    # ------------ #
    def show(self, fig=None, savefile=None, show=True,
             scaleup=2.5, ellipse_prop={}, id_to_show=None,
              **kwargs):
        """
        """
        from matplotlib.patches import Ellipse
        from astrobject.utils.mpladdon import figout
        
        def add_ellipse(ax, host_ellipse, scaleup=2.5, draw=True,**kwargs_):
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
            fig = mpl.figure(figsize=[18,10])

        axes = [fig.add_subplot(1,ndata,i+1) for i in range(ndata)]

        # -------------
        # - Properties
        show_prop = kwargs_update(dict(show_catalogue=True,show_sepobjects=False,\
                                       zoomon="target", zoom=30, zunits="kpc"),
                                       **kwargs)
        
        for id_,ax in zip(id_to_show,axes):
            if self.images[id_]["image"] is None: self._load_image_(id_)
                
            self.images[id_]["image"].show(ax =ax, show=False, **show_prop)
            ellip_data, notforced = self.get_host_ellipse(id_)
            add_ellipse(ax, ellip_data, ls= "-" if notforced else "--",
                        **kwargs_update({"lw":2},**ellipse_prop))

            ax.text(0.05,0.95, "band: %s"% self.images[id_]["image"].bandname, 
                    fontsize="large",va="top",ha="left",
                    transform=ax.transAxes, bbox=dict(facecolor='w', alpha=0.9))
            
        # = Titles
        fig.suptitle("%s ID of %s-host: %s"%(self.catalogue.source_name, self.target.name, self.host_catid[0]), 
                    y=0.90, va="top",fontsize="x-large")
        
        fig.figout(savefile=savefile, show=show)
        
    # ================== #
    # =  Properties    = #
    # ================== #

    # = Reference
    
    @property
    def refid(self):
        """ The id defined as reference image. """
        return self._properties["refid"]

    def has_refid(self):
        """ test if a reference id has been set."""
        return self.refid is not None

    @property
    def ref_host_ellipse(self):
        """ """
        return self.get_host_ellipse(self.refid)[0]
    
    @property
    def refimage(self):
        """ """
        if not self.has_refid():
            raise AttributeError("You did not set the refid")
        if self.images[self.refid]["image"] is None:
            self._load_image_(self.refid)
            
        return self.images[self.refid]["image"]

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
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
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
        if "__iter__" not in dir(photopoints):
            photopoints = [photopoints]
        if "__iter__" not in dir(ids):
            ids = [ids]*len(photopoints)
            
        if len(ids) != len(photopoints):
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
        self.photopoints[ID] = photopoint

    def writeto(self,filename,format="ascii",**kwargs):
        """
        Save the PhotoPointCollection following the astropy's Table 'write()' method.
        Such data will be accessible using the 'load()' method of the class.
        """
        complet = table.join(self.data,self.meta,join_type='left',keys='id') if self.meta is not None\
          else self.data
        complet.write(filename,format=format,**kwargs)
        
    def load(self,filename,format="ascii",**kwargs):
        """
        """
        data = table.Table.read(filename,format=format)
        ids,pps = [],[]
        # - This way to be able to call create that might be more
        # generic than add_photopoints
        for pp in data:
            d = {}
            for i,k in enumerate(data.keys()):
                d[k] = pp[i]
            ids.append(d.pop("id"))
            pps.append(get_photopoint(**d))
            
        self.create(pps,ids,**kwargs)
        if "comments" in data.meta.keys():
            self._read_table_comment_(data.meta["comments"])

    def _read_table_comment_(self, comment):
        """ """
        warnings.warn("The following comments on the table are not saved:"+"\n".join(comment))
        pass
        
    # ------------------- #
    # - Getter          - #
    # ------------------- #
    def get(self,param, mask=None, safeexit=True):
        """Loop over the photopoints (following the list_id sorting) to get the
         parameters using the photopoint.get() method.

        If safeexit is True NaN will be set to cases where the get(key) does not exist.
        Otherwise this raises a ValueError.
        
        Returns
        -------
        list
        """
        ids = self.list_id if mask is None else np.asarray(self.list_id)[mask]
        return np.asarray([self.photopoints[id_].get(param,safeexit=safeexit)
                            for id_ in ids])

    # ------------------- #
    # - Setter          - #
    # ------------------- #
    def set_meta(self,key,values, ids=None):
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
        if "__iter__" not in dir(values):
            values = [values]*nids
        elif len(values)==1:
            values = [values[0]]*nids
        elif len(values) !=nids:
            raise ValueError("the given 'values', must provide a value for each photopoints (%d)"\
                             +" or set ids to the correct id list"\
                             +" or a unique value  %d given"%(self.nsources,len(values)))
        # -- latest check
            
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
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
            
        # - Colors
        if type(cmap) is str:
            cmap = eval("mpl.cm.%s"%cmap)
        colors = kwargs.pop("color",[cmap((i-np.min(self.lbdas))/(np.max(self.lbdas)-np.min(self.lbdas))) for i in self.lbdas])
        prop = kwargs_update({"ms":10},**kwargs)
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

    @property
    def metakeys(self):
        """ This is the list of all the meta keys known by at least one photopoints
        You can access these get from self.get(key)
        """
        return np.unique(np.concatenate(self.get("meta.keys()"))).tolist()
    
    # --------------
    # - Derived
    @property
    def data(self):
        """ builds an astropy table containing the fundatemental parameters of the collection """
        # - Main
        maindata = [self.list_id,self.fluxes,self.fluxvars,self.lbdas,self.mjds,self.bandnames,
                    self.get("zp"),self.get("zpsys")]
        mainnames= ["id","flux","var","lbda","mjd","bandname","zp","zpsys"]
        # - Meta
        #metaname = self.metakeys
        #metadata = [self.get(metak) for metak in metaname]
        return table.Table(data=maindata,
                     names=mainnames)
    @property
    def meta(self):
        """  builds an astropy table containing the meta information about the sources of the collection """
        if len(self.metakeys) == 0:
            return None
        metanames = self.metakeys
        metadata = [self.get(metak) for metak in metanames]
        return table.Table(data=[self.list_id]+metadata,
                            names=["id"]+metanames)
    
