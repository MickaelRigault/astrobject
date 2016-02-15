#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from .baseobject  import BaseObject
from .instruments import instrument
from ..utils.decorators import _autogen_docstring_inheritance

__all__ = ["multibands","multimages"]

def multimages(target,images,**kwargs):
    """
    """
    astroimages = []
    for img in images:
        if type(img) == str:
            img = instrument.instrument(img)
        elif "__nature__" not in dir(img) or img.__nature__ != "Image":
            continue
        astroimages.append(img)
        
    return MultiImages(astrotarget=target,
                       images=astroimages,
                       **kwargs)

def multibands(photopoints,labels,target=None,**kwargs):
    """
    """
    return MultiBands(photopoints,labels,
                      astrotarget=target,**kwargs)

#######################################
#                                     #
# SED-Data: BASE : MultiBands         #
#                                     #
#######################################

class MultiBands( BaseObject ):
    """
    This is the object gathering several observations
    from potentially several instruments for a given target
    """
    __nature__ = "MultiBand"
    _properties_keys = ["photopoints"]
    _side_properties_keys = ["target"]
    _derived_properties_keys = []

    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,photopoints=None,labels=None,
                 astrotarget=None,
                 empty=False):
        """
        """
        self.__build__()
        if empty:
            return
        self.set_target(astrotarget)
        
        if photopoints is not None:
            if labels is None or len(labels) != len(photopoints):
                raise ValueError("Please provide a label for each photopoint")
            
            for photopoint,label in zip(photopoints,labels):
                self.add_photopoint(photopoint,label,update=False)
                
            self._update_()
            
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    # ------------------- #
    # - Get Methods     - #
    # ------------------- #
    def get_photopoint(self,label):
        """This methods return the image associated to the
        given label"""
        self._test_label_(label)
        return self.photopoints[label]

    # --------------------- #
    # add/mv/rm photopoint  #
    # --------------------- #
    def add_photopoint(self,photopoint,label,
                  force_it=False,check=True,
                  update=True):
        """
        This module enables to add a new photopoint to this object .
        This photopoint will be added if it does not already exist
        (set force_it to true to overwrite the previous one.).
        """
        if photopoint.__nature__ != "PhotoPoint":
            raise TypeError("the given 'photopoint' is not an astrobject 'PhotoPoint'")
        # -- Convert it in dictionnary
        if self.photopoints is None:
            self._properties['photopoints'] = {}
            
        self._properties['photopoints'][label] = photopoint
        self._update_()

    def rename_label(self,label,newlabel,force_it=False):
        """This methods allows you to change an photopoint label name
        *label* to an other *newlabel*. 
        """
        if label not in self.labellist:
            raise ValueError('%s is not a known label'%label)
        if newlabel in self.labellist and force_it is False:
            raise ValueError("%s already exists, find another 'newlabel'"%newlabel)
    
        self._properties['photopoints'][newlabel] = self._properties['photopoints'].pop(label)
        self._update_()
        
    def remove_label(self,label):
        """remove the given photopoint-label from the source"""
        self._test_label_(label)
        _ = self.photopoints.pop(label)
        self._update_()

    # ------------------- #
    # - Target          - #
    # ------------------- #
    def set_target(self,newtarget):
        """
        Change (or create) an object associated to the given image.
        Set newtarget to None to remove the association between this
        object and a target
        """
        if newtarget is None:
            self._side_properties['target'] = None
            return
        # -- Input Test -- #
        if newtarget.__nature__ != "AstroTarget":
            raise TypeError("'newtarget' should be (or inherite) an AstroTarget")
        
        self._side_properties["target"] = newtarget.copy()

        
    # ------------------- #
    # - I/O Methods     - #
    # ------------------- #
    def show(self,savefile=None,ax=None,
             toshow="flux",**kwargs):
        """
        """
        # ----------------
        # - The Axis 
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure

        [self.display_photopoint(ax,label,toshow=toshow,**kwargs)
         for label in self.labellist]
        fig.show()         
        
    def display_photopoint(self,ax,label,**kwargs):
        """This method enable to add the band associated
        to the given *label* to the given *ax*"""
        # - Is that an existing label ?
        self._test_label_(label)
        photo = self.get_photopoint(label)
        photo.display(ax,**kwargs)
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def photopoints(self):
        return self._properties['photopoints']
    
    @property
    def labellist(self):
        if self.has_images is False:
            return None
        return np.sort(self.images.keys())
    
    def has_photopoints(self):
        return False if self.photopoints is None \
          or len(self.labellist) == 0 else True
          
    def has_label(self,label):
        """
        """
        return not (self.has_images() is False or
                    label not in self.images.keys())
        
    @property
    def target(self):
        return self._side_properties["target"]
    
    def has_target(self):
        return False if self.target is None \
          else True

    # ------------
    # - Derived
    # from Target
    @property
    def arcsec_per_kpc(self):
        if not self.has_target():
            return AttributeError("no 'target' loaded")
        return self.target.arcsec_per_kpc

    # from Images
    @property
    def instruments(self):
        return [self.get_photopoint(label).instrument_name
                for label in self.labellist]
    
    @property
    def lbdas(self):
        return [self.get_photopoint(label).lbda
                for label in self.labellist]
    # =========================== #
    # = Internal Tools          = #
    # =========================== #
    def _set_target_to_image_(self,label):
        """This methods enable to set the object's astrotarget
        to the given labeled image"""
        ### DONE ###
        self._test_label_(label)
        self.images[label].set_target(self.target)
        
    # -------------
    # - Image Label
    def _test_label_(self,label):
        """This internal methods raise an ValueError exception
        if the label in not known"""
        if not self.has_label(label):
            raise ValueError("%s is not a known 'label'"%label)


#######################################
#                                     #
# SED-Data: Advanced : MultiImages    #
#                                     #
#######################################

class MultiImages( MultiBands ):
    """This object parses the given images to
    build access the MultiBand tools"""
    
    
    # =========================== #
    # = Constructor             = #
    # =========================== #
    def __init__(self,astrotarget=None,
                 images=None,labels=None,
                 empty=False):
        """
        """
        self.__build__()
        if empty:
            return
        self.set_target(astrotarget)
        
        if images is not None:
            labels = [None]*len(images) if labels is None or\
              len(labels) != len(images) else labels
            
            for image,label in zip(images,labels):
                self.add_image(image,label,update=False)
                
            self._update_()
            
    def __build__(self):
        """
        """
        self._properties_keys.append("images")
        super(MultiImages,self).__build__()
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    @_autogen_docstring_inheritance(MultiBands.set_target,"MultiBands.set_target ; HERE TARGET CANNOT BE NONE")
    def set_target(self,newtarget,**kwargs):
        """
        """
        if newtarget is None:
            raise ValueError("The 'MultiImages' object requires to have a target")
        
        super(MultiImages,self).set_target(newtarget,**kwargs)
        
    def get_image(self,label):
        """This methods return the image associated to the
        given label"""
        self._test_label_(label)
        return self.images[label]

    def measure_target_photopoint(self,label,r_pixels,**kwargs):
        """
        """
        img = self.get_image(label)
        return img.get_target_photopoint(r_pixels=r_pixels,**kwargs)
    
    # ----------------------- #
    # - Add / Mv / Rv Image - #
    # ----------------------- #
    def add_image(self,image,label=None,
                  force_it=False,check=True,
                  update=True,r_pixels=5,**kwargs):
        """
        This module enables to add to the object a new image.
        This image will be added if it does not already exist
        (see force_it) and if the current object is contained with
        the new image field of view.
        """
        if image.__nature__ != "Image":
            raise TypeError("the given 'image' is not an astrobject 'Image'")
        # -- Convert it in dictionnary
        if self.images is None:
            self._properties['images'] = {}

        label = image.bandname if label is None else label
        # -- Set the image
        self._properties['images'][label] = image
        # -- Insert the target inthere
        self._set_target_to_image_(label)
        # -- Get the corresponding photopoint
        self.add_photopoint(image.get_target_photopoint(r_pixels,**kwargs),
                            label)
        # - Done
        self._update_()

    @_autogen_docstring_inheritance(MultiBands.rename_label,"MultiBands.rename_label")
    def rename_label(self,label,newlabel,force_it=False):
        #
        # - Add the images label update
        #
        super(MultiImages,self).rename_label(label,newlabel,force_it=force_it)
        self._properties['images'][newlabel] = self._properties['images'].pop(label)
        self._update_()
        
    @_autogen_docstring_inheritance(MultiBands.remove_label,"MultiBands.remove_label") 
    def remove_label(self,label):
        """remove the given image-label from the images source"""
        super(MultiImages,self).remove_label(label)
        _ = self.images.pop(label)
        self._update_()
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def images(self):
        return self._properties["images"]
        
    def has_images(self):
        return False if self.images is None \
          or len(self.images.keys()) == 0 else True

    @property
    def bands(self):
        return [self.get_image(label).bandname
                for label in self.labellist]
    # =========================== #
    # = Internal Tools          = #
    # =========================== #
    def _set_target_to_image_(self,label):
        """This methods enable to set the object's astrotarget
        to the given labeled image"""
        self._test_label_(label)
        self.images[label].set_target(self.target)
