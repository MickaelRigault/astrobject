#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This pupose of this method is the have a generator of fake astrotarget, like SN"""

from ..astrobject.baseobject import BaseObject,astrotarget
from ..astrobject.transient import transient
from ..utils import random, skyplot

__all__ = ["transient_generator","generate_transients"]

def transient_generator(npoints, zrange,
                        zpdf="flat",
                        ra_range=[-180,180], dec_range=[-90,90],
                        **kwargs):
    """
    This model return the object that enable to create and change
    the kind of transient you which to set in the sky.

    # - HERE COPY PASTE THE TransientGenerator INIT - #
    # - TO BE DONE
    
    """
    return TransientGenerator(ntransients=npoints, zrange=zrange,
                              zpdf=zpdf,ra_range=ra_range,dec_range=dec_range,
                              **kwargs)
    
def generate_transients(npoints, zrange,**kwargs):
    """This module call transient_generator to create the
    TransientGenerator object and then returns the associated
    TransientGenerator.transients

    # - HERE COPY PASTE the transient_generator docstring
    """
    return transient_generator(npoints, zrange,**kwargs).transients


#######################################
#                                     #
# Generator: Any Transient            #
#                                     #
#######################################
class TransientGenerator( BaseObject ):
    """
    """
    _properties_keys = ["ntransients","redshift_coverage",
                        "skycoverage"]
        
    _derived_keys    = ["transients","transientsources"]

    def __init__(self,ntransients=10,
                 zrange=[0.0,0.2], zpdf="flat", # How deep 
                 ra_range=(-180,180),dec_range=(-90,90), # Where, see also kwargs
                 empty=False,**kwargs):
        """
        """
        self.__build__()
        if empty:
            return
        self.create(ntransients,zrange,
                    zpdf=zpdf,
                    ra_range=ra_range,dec_range=dec_range,
                    **kwargs)

    def create(self,ntransients,zrange,zpdf="flat",
               ra_range=(-180,180),dec_range=(-90,90),mw_exclusion=10):
        """
        """
        # == Add the Input Test == #
        #   TO BE DONE
        
        # *************** #
        # * Create      * #
        # *************** #
        self._properties["ntransient"] = ntransients
        self._properties["redshift_coverage"] = {"zrange":zrange,
                                                 "pdfkind":zpdf}

        # -- This will be directly used as random.radec inputs
        self._properties["sky_coverage"] = {"ra_range":ra_range,
                                           "dec_range":dec_range,
                                           "mw_exclusion":mw_exclusion,
                                           }
        
        self._update_()

    # --------------------------- #
    # - Plots Methods           - #
    # --------------------------- #
    def show_skycoverage(self):
        """This function enable to draw on the sky the position of the
        transients"""
        
        
    def display_target_onsky(self,ax):
        """This function enable to draw on the given ax the sky coverage of
        the given transients"""

    
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _update_(self):
        """This module create the derived values based on the
        fundamental ones"""
        if self.ntransient is None:
            self._derived_properties["transientsources"] = None
            return
        if self.zrange is None:
            self._derived_properties["transientsources"] = None
            return
        
        simul_z = random.redshift(self.ntransient,self.zrange,
                                  pdfkind=self.redshift_coverage["pdfkind"])
        simul_ra, simul_dec = random.radec(self.ntransient,**self.sky_coverage)

        # The sources are created and not the transient for save space
        self._derived_properties["transientsources"] = \
                                 [ dict(lightcurve=None,name="simul%d"%i,
                                        zcmb=z_,ra=ra_,dec=dec_)
                                    for i,z_,ra_,dec_ in zip(range(self.ntransient),
                                                            simul_z,simul_ra, simul_dec)]
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def zrange(self):
        """redshift range used to draw transient"""
        if self.redshift_coverage is None:
            raise AttributeError("'redshift_coverage' has not been defined")
        return self.redshift_coverage["zrange"]

    # -- The fundamental Parameters
    @property
    def ntransient(self):
        """number of transient requested"""
        return self._properties["ntransient"]

    # --------------------
    # - Target Coverage
    @property
    def redshift_coverage(self):
        return self._properties["redshift_coverage"]
        
    @property
    def sky_coverage(self):
        """where the transient could be drawn"""
        return self._properties["sky_coverage"]

    # --------------------
    # - Derived Properties
    @property
    def transients(self):
        """loops over the transientsources to load the transients objects"""
        if self.transientsources is None:
            raise AttributeError("No 'transientsources' defined.")
        return [transient(**s) for s in self.transientsources]
    
    @property
    def transientsources(self):
        """dictionary containing the fundamental parameters that enable to load the transient objects"""
        return self._derived_properties["transientsources"]
