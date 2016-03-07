#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module contain the collection of astrobject"""
import numpy as np

from .baseobject import BaseObject
from .instruments import instrument as inst
from ..utils.tools import kwargs_update
from ..utils.shape import draw_polygon

__all__ = ["ImageCollection"]

#######################################
#                                     #
# Base Collection                     #
#                                     #
#######################################
class BaseCollection(BaseObject):
    """
    """
    __nature__ = "Collection"

    _properties_keys = ["handler"]
    
    def __init__(self):
        """
        """
        self.__build__()
        return 

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
class Collection( BaseCollection ):
    """
    """
    def __build__(self,*args,**kwargs):
        """
        """
        self._side_properties_keys.append("target")
        super(Collection,self).__build__(*args,**kwargs)
        
    # ----------------- #
    # - Setter        - #
    # ----------------- #
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
        """This is the list of the known data id (equivalent to list_sources)"""
        return self.list_sources
    # ---------------------
    # - Target
    @property
    def target(self):
        """The associated target if loaded"""
        return self._side_properties['target']
    
    def has_target(self):
        return self.target is not None

#######################################
#                                     #
# Image Collection                    #
#                                     #
#######################################
class ImageCollection( Collection ):
    """
    """
    _side_properties_keys = ["catalogue"]
    
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
    def get_image(self,id,load_catalogue=True):
        """
        """
        self._test_id_(id)
        # --------------------
        # - image already exist
        if self.images[id]["image"] is None:
            self._load_image_(id)
            return self.get_image(id)
        
        im = self.images[id]["image"]
        
        # ---------------
        # Target
        if not im.has_target() and self.has_target():
            im.set_target(self.target)
        # ---------------
        # Catalogue
        if load_catalogue and self.has_catalogue() and not im.has_catalogue():
            # --- Check if the current catalogue is big enough
            im.set_catalogue(self.catalogue.copy())

        return im

    def get_catalogue(self,kind="all", isolated_only=False, stars_only=False):
        """
        Get a copy of a catalogue. If could be the entire catalogue (kind=None or kind="all")
        or the shared fov catalogue (kind="shared" or "union") or the combined fov catalogue
        (kind="combined" or "intersection")
        """
        if not self.has_catalogue():
            raise AttributeError("No Catalogue loaded")
        
        if kind is None or kind.lower() in ["all"]:
            contours = None
        elif kind.lower() in ["shared","union"]:
            contours = self.contours_shared
        elif kind.lower() in ["combined","intersection"]:
            contours = self.contours_combined
        else:
            raise ValueError("I could not parse the given kind %s (must be 'all', 'shared' or 'combined')"%kind)
        
        return self.catalogue.get_subcatalogue(contours=contours,isolated_only=isolated_only, stars_only=stars_only)
        
        
    
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
        im.set_target(target)
        return im

    def get_target_photopoints(self,radius_kpc,target=None,
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
        if not onflight:
            ids = self.get_target_ids(target=target)
            images = [self.get_image(_id,load_catalogue=False)
                      for _id in ids]
        else:
            target_imcoll = self.get_target_collection(target=target)
            images = [target_imcoll(_id,load_catalogue=False)
                      for _id in target_imcoll.list_id]

        
        photopoints = PhotoPointCollection(photopoints=[image_.get_target_photopoint(radius_kpc/image_.target.arcsec_per_kpc/image_.pixel_size_arcsec.value,
                                            **kwargs) for image_ in images])
        
        # ------------------------------
        # - Shall we delete the files ?
        if onflight:
            del images
            del target_imcoll
            
        return photopoints
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
        from ..utils.mpladdon import wcsplot,figout,mpl
        
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
          
        return inst.catalogue(source=source,radec=radec,radius="%sd"%radius,
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
        
    def _load_image_(self,id):
        """
        """
        self.images[id]["image"] = inst.instrument(self.images[id]["file"])
        
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
# Image Collection -> SED SOURCE      #
#                                     #
#######################################
class PhotoPointCollection( Collection ):
    """
    """
    
    def __init__(self, photopoints=None,empty=False):
        """
        """
        self.__build__()
        if empty:
            return
        if photopoints is not None:
            self.create(photopoints)
    # =============================== #
    # = Main Methods                = #
    # =============================== #

    # ------------------- #
    # - I/O PhotoPoint  - #
    # ------------------- #
    def create(self,photopoints):
        """
        """
        if "__iter__" not in dir(photopoints):
            photopoints = [photopoints]
        [self.add_photopoint(p_) for p_ in photopoints]
        
    def add_photopoint(self,photopoint):
        """
        This method enables to register the given photpoint in the self.photopoints
        container (dict).
        
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
        # --------------
        # - Define ID
        ID = photopoint.bandname
        if self.has_data():
            if ID in self.list_id:
                i = 1
                while ID+"-%d"%i in self.list_id:
                    i+=1
                ID = ID+"-%d"%i
        # -------------------
        # - Record the point        
        self.photopoints[ID] = photopoint

    # ------------------- #
    # - Getter          - #
    # ------------------- #


        
    # ------------------- #
    # - Plot Stuff      - #
    # ------------------- #
    def show_sed(self,savefile=None,toshow="flux",ax=None,
                cmap="jet",show=True,**kwargs):
        """
        """
        kwargs["function_of_time"] = False
        return self._show_(savefile=savefile,toshow=toshow,ax=ax,
                    cmap=cmap,show=show,**kwargs)

    def show_lightcurve(self,savefile=None,toshow="flux",ax=None,
                        cmap="jet",show=True,**kwargs):
        """
        """
        kwargs["function_of_time"] = True
        return self._show_(savefile=savefile,toshow=toshow,ax=ax,
                    cmap=cmap,show=show,**kwargs)

        
    # ========================== #
    # = Internal Stufff        = #
    # ========================== #
    def _show_(self,savefile=None,toshow="flux",ax=None,
                cmap="jet",show=True,**kwargs):
        """
        """
        import matplotlib.pyplot as mpl
        from ..utils.mpladdon import figout
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

        
    def _get_list_prop_(self,param):
        """
        """
        return [eval("self.photopoints[id_].%s"%param)
                for id_ in self.list_id]
    
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
        """This returns an array of the recorded images' bandname"""
        # -- This might be included in the add/remove image tools
        return self._get_list_prop_("bandname")

    @property
    def nbands(self):
        """The amount of diffenrent recorded bands"""
        return np.unique(self.bandnames)

    # -- Flux Variances etc.
    @property
    def lbdas(self):
        """The wavelength of every photopoints"""
        return self._get_list_prop_("lbda")

    @property
    def fluxes(self):
        """The fluxes of every photopoints"""
        return self._get_list_prop_("flux")
    
    @property
    def fluxvars(self):
        """The flux variances of every photopoints"""
        return self._get_list_prop_("var")
    
    @property
    def mags(self):
        """The (AB) mag of every photopoints"""
        return self._get_list_prop_("mag")

    @property
    def magvars(self):
        """The (AB) mag variances of every photopoints"""
        return self._get_list_prop_("magvar")
    
    # -- Times
    @property
    def mjds(self):
        """The modified Julian Dates of every photopoints"""
        return self._get_list_prop_("mjd")
