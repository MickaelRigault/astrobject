#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""This module contain the collection of astrobject"""
import numpy as np

from .baseobject import BaseObject
from .instruments import instrument as inst


#######################################
#                                     #
# Image Collection                    #
#                                     #
#######################################
class ImageCollection( BaseObject ):
    """
    """
    _properties_keys = ["images"]
    _side_properties_keys = ["catalogue","target"]
    
    def __init__(self,images=None,empty=False):
        """
        """
        self.__build__()
        if empty:
            return
        
        self.create(images)

        
    def create(self,images):
        """
        Loop over the images end add them to the instance. See add_image
        """
        if "__iter__" not in dir(images):
            images = [images]
        
        [self.add_image(i_) for i_ in images
         if i_ is not None]
        
    # ========================= #
    # = Images IO             = #
    # ========================= #
    def add_image(self,new_image):
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

    def remove_image(self,id):
        """
        """
        self._test_id_(id)
        
        to_rm = self.images.pop(id)
        del to_rm

    def load_images(self,verbose=True):
        """
        """
        if verbose:
            print "you have %s images to load. This might take some time"%(self.nimages)
        [self._load_image_(id_) for id_ in self._imageids
         if self.images[id_]["image"] is None]
        # - done
            
    # ========================== #
    # = Get                    = #
    # ========================== #
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
            if np.max(im.wcs.getHalfSizeDeg())*2 > self._radiuscatalogue:
                self.download_catalogue(radius=np.max(im.wcs.getHalfSizeDeg())*2)
                
            im.set_catalogue(self.catalogue.copy())

        # --------------
        # - TEMPORARY
        return im

    # ========================== #
    # = Set                    = #
    # ========================== #
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
                           force_it=False,set_it=True,
                           radec=None,radius=2,
                           **kwargs):
        """
        remark, a catalogue size of 1 or 2 degree have a comparable speed (~1.5s)
        """
        if not self.has_target():
            if radec is None or radius is None:
                raise ValueError("No Target loaded, no ra,dec,radius given... no catalogue to load")
            ra,dec = radec
        else:
            ra,dec = self.target.ra,self.target.dec

        self._radiuscatalogue = radius
        
        radec = "%s %s"%(ra,dec)
        cat = inst.catalogue(source=source,radec=radec,radius="%sd"%self._radiuscatalogue,
                             **kwargs)
        if not set_it:
            return cat
        
        self.set_catalogue(cat,force_it=force_it)

    def show_skypatches(self,ax=None,
                        savefile=None,show=True,
                        fc="None",ec="k",**kwargs):
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
              for id_ in self._imageids if self.images[id_]["wcs"] is not None]

        # ----------------------
        # -  output
        fig.figout(savefile=savefile,show=show)
    # ========================== #
    # = INTERNAL               = #
    # ========================== #
    def _test_id_(self,id):
        """
        """
        if not self.has_images:
            raise AttributeError("no images loaded")
        
        if id not in self.images.keys():
            raise ValueError("%s is not a known image ID"%id +"\n"\
                             "These are: "+", ".join(self.images.keys()))
        return True

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
    # -- Images
    @property
    def nimages(self):
        return len(self.images.keys())
    
    @property
    def images(self):
        """
        """
        if self._properties["images"] is None:
            self._properties["images"] = {}
        return self._properties["images"]

    def has_images(self):
        return len(self.images.keys()) > 0
    
    @property
    def _imageids(self):
        return np.sort(self.images.keys())
    
    # -- Target
    @property
    def target(self):
        return self._side_properties['target']

    def has_target(self):
        return not self.target is None
    
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
                for id_ in self._imageids]
    

    
