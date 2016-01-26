#! /usr/bin/env python
# -*- coding: utf-8 -*-


import warnings
import numpy as np

import sncosmo
from astropy.table import Table

from ..astrobject.baseobject import BaseObject
from ..utils.tools import kwargs_update
from .simultarget import sn_generator,transient_generator


__all__ = ["SimulSurvey"] # to be changed

#######################################
#                                     #
# Survey: Simulation Base             #
#                                     #
#######################################
class SimulSurvey( BaseObject ):
    """
    Basic survey object
    (far from finished)
    """
    _properties_keys         = ["generator","instruments","cadence"]
    _side_properties_keys    = []
    _derived_properties_keys = ["observations"]
    
    def __init__(self,generator=None,
                 instprop=None,
                 empty=False):
        """
        Parameters:
        ----------
        generator: [simultarget.transient_generator or derived child like sn_generator]

        
        """
        self.__build__()
        if empty:
            return
    
    # =========================== #
    # = Main Methods            = #
    # =========================== #

    # ---------------------- #
    # - Get Methods        - #
    # ---------------------- #
    def get_lightcurves(self):
        """
        """
        if not self.is_set():
            raise AttributeError("cadence, genetor or instrument not set")
        self.observations.sort("time")
        params = [{"z":self.generator.zcmb[i],
                    "c":self.generator.color[i],
                    "x1":self.generator.x1[i],
                    "x0":self.generator.lightcurve["x0"][i],
                    "t0":self.generator.mjd[i]}
                    for i in range(self.generator.ntransient)
                    ]
        return sncosmo.realize_lcs(self.observations,
                        sncosmo.Model(source=self.generator.lightcurve_model),
                        params)
        
    

        
    # ---------------------- #
    # - Setter Methods     - #
    # ---------------------- #
    
    # -------------
    # - Targets
    def set_target_generator(self, generator):
        """
        """
        if "__nature__" not in dir(generator) or\
          generator.__nature__ != "TransientGenerator":
            raise TypeError("generator must be an astrobject TransientGenerator")

        if not generator.has_lightcurves():
            warnings.warn("No lightcurves set in the given transient generator")

        self._properties["generator"] = generator

    # -------------
    # - Cadences
    def set_cadence(self,cadence):
        """
        """
        # ----------------------
        # - Load cadence here
        if "__nature__" not in dir(cadence) or \
          cadence.__nature__ != "Cadence":
            raise TypeError("the input 'cadence' must be an astrobject Cadence")
        self._properties["cadence"] = cadence
        
        # ----------------------------
        # - Set back the observations
        self._reset_observations_()

    # -------------
    # - Instruments
    def set_instruments(self,properties):
        """
        properties must be a dictionary containing the
        instruments' information (bandname,gain,zp,zpsys) related
        to each bands
        
        example..
        ---------
        properties = {"desg":{"gain":1,"zp":30,"zpsys":'ab'},
                      "desr":{"gain":1,"zp":30,"zpsys":'ab'}}
        """
        for band,d in properties.items():
            gain,zp,zpsys=d.pop("gain"),d.pop("zp"),d.pop("zpsys","ab")
            if gain is None or zp is None:
                raise ValueError('gain or zp is None or not defined for %s'%band)
            self.add_instrument(band,gain,zp,zpsys,
                                update=False,**d)
        
        self._reset_observations_()

    # ---------------------- #
    # - Add Stuffs         - #
    # ---------------------- #
    def add_instrument(self,bandname,gain,zp,zpsys="ab",
                       force_it=True,update=True,**kwargs):
        """
        kwargs could be any properties you wish to save with the instrument
        """
        if self.instruments is None:
            self._properties["instruments"] = {}
            
        if bandname in self.instruments.keys() and not force_it:
            raise AttributeError("%s is already defined."+\
                                 " Set force_it to True to overwrite it. ")
                                 
        instprop = {"gain":gain,"zp":zp,"zpsys":zpsys}
        self.instruments[bandname] = kwargs_update(instprop,**kwargs)
        
        if update:
            self._reset_observations_()
            
    # ---------------------- #
    # - Recover Methods    - #
    # ---------------------- #
    #def recover_targets(self):
    #    """
    #    bunch threshold...
    #    """
    #
    #def recover_lightcurves(self):
    #    """
    #    """

    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _update_lc_(self):
        """
        """
        # -----------------------------
        # -- Do you have all you need ?
        if not self.is_set():
            return

    def _reset_observations_(self):
        """
        """
        self._derived_properties["observations"]
        
    def _load_observations_(self):
        """
        """
        # -------------
        # - Input test
        if self.cadence is None or self.instruments is None:
            raise AttributeError("Cadence or Instruments is not set.")
        
        # -----------------------
        # - Check if instruments exists
        all_instruments = np.unique(self.cadence.bands)
        if not np.all([i in self.instruments.keys() for i in all_instruments]):
            raise ValueError("Some of the instrument in cadence have not been defined."+"\n"+
                             "given instruments :"+", ".join(all_instruments.tolist())+"\n"+
                             "known instruments :"+", ".join(self.instruments.keys()))
        
        # -----------------------
        # - Lets build the table
        self._derived_properties["observations"] = Table(
            {"time":self.cadence.mjds,
             "band":self.cadence.bands,
             "skynoise":self.cadence.skynoises,
             "gain":[self.instruments[b]["gain"] for b in self.cadence.bands],
             "zp":[self.instruments[b]["zp"] for b in self.cadence.bands],
             "zpsys":[self.instruments[b]["zpsys"] for b in self.cadence.bands]
            })
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def instruments(self):
        """The basic information relative to the instrument used for the survey"""
        return self._properties["instruments"]
    
    @property
    def generator(self):
        """The instance that enable to create fake targets"""
        return self._properties["generator"]

    @property
    def cadence(self):
        """This is a table containing where the telescope is pointed with which band"""
        return self._properties["cadence"]

    def is_set(self):
        """This parameter is True if this has cadence, instruments and genetor set"""
        return not (self.instruments is None or \
                    self.generator is None or \
                    self.cadence is None)
                    
    # ------------------
    # - Derived values
    @property
    def observations(self):
        """Observations derived from cadence and instrument properties.
        Remark that the first time this is called, observations will be recorded"""
        
        if self._derived_properties["observations"] is None:
            self._load_observations_()
            
        return self._derived_properties["observations"]
    
