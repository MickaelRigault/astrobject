#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module contain the collection of astrobject"""
import numpy as np

from .baseobject import BaseObject
from .instruments import instrument as inst
from ..utils.tools import kwargs_update


__all__ = ["ImageCollection"]

#######################################
#                                     #
# Base Collection                     #
#                                     #
#######################################
class Collection( BaseObject ):
    """
    """
    __nature__ = "Collection"

    _properties_keys = ["handler"]
    _side_properties_keys = ["target"]
    
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
        self._test_id_(id)
        to_rm = self._handler.pop(id)
        del to_rm

    def rename(self,id,newid,force_it=False):
        """This methods enable to change the id with which the data is accessed.
        Remark: If the 'newid' already exists, an exception is raised except
        if force_it is set to true. In that a case, this will overwrite
        the 'newid' entry.
        """
        # -- Test if it is ok
        self._test_id_(id)
        if newid in self.list_id and not force_it:
            raise ValueError("the newid '%s' already exists. Set force_it to true to overwrite it")
        
        # -- Do the switch
        self._handler[newid] = self._handler.pop(id)
        
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

    # ========================== #
    # = Internal Tools         = #
    # ========================== #
    
    def _test_id_(self,id):
        """This methods tests if the given 'id' exists"""
        if id not in self.list_id:
            raise ValueError("%s is not a known image ID"%id +"\n"\
                             "These are: "+", ".join(self.list_id))
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
    def ndata(self):
        return len(self._handler.keys())
        
    def has_data(self):
        return len(self._handler.keys()) > 0
    
    @property
    def list_id(self):
        """This is the list of the known data id"""
        if not self.has_data():
            raise AttributeError("no data loaded")
        
        return np.sort(self._handler.keys()).tolist()
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
            self._add_image_file_(new_image)
        
        elif "__nature__" in dir(new_image) and new_image.__nature__ != "Image":
            self._add_astroimage_(new_image)
        else:
            # -- Issue
            raise TypeError("the given new_image must be a image-file or an astrobject imnage")
        if load_catalogue:
            self.download_catalogue()
            

    def load_images(self,verbose=True):
        """
        """
        [self._load_image_(id_) for id_ in self.list_id
         if self.images[id_]["image"] is None]
        # - done
            
    # ========================== #
    # = Get                    = #
    # ========================== #
    def get_target_collection(self,target=None):
        """
        return a collection of images associated to the target.
        if no target is given and a target has been set (set_target) the current
        target will be used. If you give here a target, it will use it *without*
        changing the current target. This method raise a ValueError is no target
        accessible.
        """
        if target is None:
            if not self.has_target():
                raise ValueError("no target set to the instance and no input target")
            target = self.target
        elif "__nature__" not in dir(target) or target.__nature__ != "AstroTarget":
            raise TypeError("The given 'target' must be an astrobject AstroTarget")
        
        good_images = [self.images[id_]["file"] for id_ in self.list_id
                       if self.images[id_]["wcs"].coordsAreInImage(target.ra,target.dec)]
            
        im = ImageCollection(good_images)
        
        im.set_target(target)
        return im
            
    def get_image(self,id):
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
        if not im.has_catalogue() and self.has_catalogue():
            # --- Check if the current catalogue is big enough
            im.set_catalogue(self.catalogue.copy())

        return im

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


    def download_catalogue(self,source="sdss",
                           **kwargs):
        """
        """
        for id_ in self.list_id:
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
        pl = [ax.wcsplot(self.images[id_]["wcs"],fc=fc,ec=ec,**kwargs)
              for id_ in self.list_id if self.images[id_]["wcs"] is not None]

        # ----------------------
        # - Catalogue
        if self.has_catalogue() and show_catalogue:
            self.catalogue.display(ax)
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
        print ID
        self.images[ID] = {
            "file":imagefile,
            "image":None

            }
        try:
            self.images[ID]["wcs"] = inst.get_instrument_wcs(imagefile)
        except:
            print "WARNING Failure to load the wcs solution"
            self.images[ID]["wcs"] = None
        
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
    
    
#######################################
#                                     #
# Image Collection -> SED SOURCE      #
#                                     #
#######################################
class PhotoCollection( Collection ):
    """
    """
    
    def __init__(self, photopoints=None,empty=True):
        """
        """
        self.__build__()
        
        
    
