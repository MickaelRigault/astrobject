#! /usr/bin/env python
# -*- coding: utf-8 -*-


import warnings
import numpy as np

import sncosmo
from astropy.table import Table, vstack

from ..astrobject.baseobject import BaseObject
from ..utils.tools import kwargs_update
from ..utils.plot.skybins import SurveyField, SurveyFieldBins 
from .simultarget import sn_generator,transient_generator


__all__ = ["SimulSurvey", "SurveyPlan"] # to be changed

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
            raise AttributeError("cadence, generator or instrument not set")
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
    
#######################################
#                                     #
# Survey: Plan object                 #
#                                     #
#######################################
class SurveyPlan( BaseObject ):
    """
    Survey Plan
    contains the list of observation times, bands and pointings and
    can return that times and bands, which a transient is observed at/with.
    A list of fields can be given to simplify adding observations and avoid 
    lookups whether an object is in a certain field.
    Currently assumes a single instrument, especially for FoV width and height.
    [This may be useful for the cadence property of SimulSurvey]
    """
    _properties_keys         = ["cadence", "width", "height"]
    _side_properties_keys    = ["fields"]
    _derived_properties_keys = []
    
    def __init__(self, mjd=None, ra=None, dec=None, band=None, obs_field=None,
                 width=7., height=7., fields=None, empty=False):
        """
        Parameters:
        ----------
        TBA
        
        """
        self.__build__()
        if empty:
            return
    
        self.create(mjd=mjd,ra=ra,dec=dec,band=band,obs_field=obs_field,
                    fields=fields)

    def create(self, mjd=None, ra=None, dec=None, band=None, obs_field=None,
               width=7., height=7., fields=None):
        """
        """
        self._properties["width"] = float(width)
        self._properties["height"] = float(height)
        self.set_fields(**fields)

        self.add_observation(mjd,band,ra=ra,dec=dec,field=obs_field)

    # =========================== #
    # = Main Methods            = #
    # =========================== #

    # ---------------------- #
    # - Get Methods        - #
    # ---------------------- #
    
    # ---------------------- #
    # - Setter Methods     - #
    # ---------------------- #
    def set_fields(self, ra=None, dec=None, **kwargs):
        """
        """
        kwargs["width"] = kwargs.get("width", self.width)
        kwargs["height"] = kwargs.get("height", self.height)

        self._side_properties["fields"] = SurveyFieldBins(ra, dec, **kwargs)

        if self.cadence is not None and np.any(np.isnan(self.cadence['field'])):
            warnings.warning("cadence was already set, field pointing will be updated")
            self._update_field_radec()

    def add_observation(self, mjd, band, ra=None, dec=None, field=None):
        """
        """
        if ra is None and dec is None and field is None:
            raise ValueError("Either field or ra and dec must to specified.")
        elif ra is None and dec is None:
            if self.fields is None:
                raise ValueError("Survey fields not defined.")
            else:
                ra = self.fields.ra[field]
                dec = self.fields.dec[field]
        elif field is None:
            field = np.array([np.nan for r in ra])

        new_obs = Table({"MJD": mjd,
                         "Band": band,
                         "RA": ra,
                         "Dec": dec,
                         "Field": field})

        if self._properties["cadence"] is None:
            self._properties["cadence"] = new_obs
        else:
            self._properties["cadence"] = vstack((self._properties["cadence"], 
                                                  new_obs))

    # ================================== #
    # = Observation time determination = #
    # ================================== #
    def observed_on(self, ra, dec, mjd_range=None):
        """
        mjd_range must be (2,N)-array 
        where N is the length of ra and dec
        """
        single_coord = None

        # first get the observation times and bands for pointings without a
        # field number use this to determine whether ra and dec were arrays or
        # floats (since this is done in SurveyField.coord_in_field there is no
        # need to redo this)
        for k, obs in enumerate(self.cadence[np.isnan(self.cadence["Field"])]):
            tmp_f = SurveyField(obs["RA"], obs["Dec"], 
                                self.width, self.height)
            b = tmp_f.coord_in_field(ra, dec)
            
            # Setup output as dictionaries that can be converted to Tables and
            # sorted later
            if k == 0:
                if type(b) is np.bool_:
                    single_coord = True
                    out = {'MJD': [], 'Band': []}
                else:
                    single_coord = False
                    out = [{'MJD': [], 'Band': []} for r in ra]

            if single_coord:
                if b:
                    out['MJD'].extend(obs['MJD'].quantity.values)
                    out['Band'].extend(obs['Band'].quantity.values)
            else:
                for l in np.where(b)[0]:
                    out[l]['MJD'].extend(obs['MJD'].quantity.values)
                    out[l]['Band'].extend(obs['Band'].quantity.values)

        # Now get the other observations (those with a field number)
        if (self.fields is not None and 
            not np.all(np.isnan(self.cadence["Field"]))):
            b = self.fields.coord2field(ra, dec)
            
            # if all pointings were in fields create new dicts, otherwise append
            if single_coord is None:
                if type(b) is not list:
                    single_coord = True
                    out = {'MJD': [], 'Band': []}
                else:
                    single_coord = False
                    out = [{'MJD': [], 'Band': []} for r in ra]
            
            if single_coord:
                for l in b:
                    mask = (self.cadence['Field'] == l)
                    out['MJD'].extend(self.cadence['MJD'][mask].quantity.value)
                    out['Band'].extend(self.cadence['Band'][mask])
            else:
                for k, idx in enumerate(b):
                    for l in idx:
                        mask = (self.cadence['Field'] == l)
                        out[k]['MJD'].extend(self.cadence['MJD'][mask].quantity.value)
                        out[k]['Band'].extend(self.cadence['Band'][mask])

        # Make Tables and sort by MJD
        if single_coord:
            table = Table(out, meta={'RA': ra, 'Dec': dec})
            idx = np.argsort(table['MJD'])
            if mjd_range is None:
                return table[idx]
            else:
                t = table[idx]
                return t[(t['MJD'] >= mjd_range[0]) &
                         (t['MJD'] <= mjd_range[1])]
        else:
            tables = [Table(a, meta={'RA': r, 'Dec': d}) for a, r, d 
                      in zip(out, ra, dec)]
            idx = [np.argsort(t['MJD']) for t in tables]
            if mjd_range is None:
                return [t[i] for t, i in zip(tables, idx)]
            else:
                ts = [t[i] for t, i in zip(tables, idx)]
                return [t[(t['MJD'] >= mjd_range[0][k]) &
                          (t['MJD'] <= mjd_range[1][k])] 
                        for k, t in enumerate(ts)]

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def cadence(self):
        """Table of observations"""
        return self._properties["cadence"]

    @property
    def width(self):
        """field width"""
        return self._properties["width"]

    @property
    def height(self):
        """field height"""
        return self._properties["height"]

    @property
    def fields(self):
        """Observation fields"""
        return self._side_properties["fields"]
