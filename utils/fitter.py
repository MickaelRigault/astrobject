#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes for lightcurve fits and their statistics

[Currently not using the sncosmo model covariance option when fitting.
 This will lead to underestimated fit uncertainties.]
"""

import numpy as np
import sncosmo
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, Column

from ..astrobject.baseobject import BaseObject

__all__ = ["Fitter"]

class Fitter( BaseObject ):
    """
    Basic lightcurve fitter class
    """

    _properties_keys = ['model', 'cosmo', 'lightcurves']

    _side_properties = ['source', 'fit_param']

    _derived_properties_keys = ['fit', 'raw_fit', 'true']

    def __init__(self, lightcurves, source='salt2', model=None, 
                 cosmo=FlatLambdaCDM(H0=70, Om0=0.3),
                 fit_param=['t0', 'x0', 'x1', 'c']):
        """
        Arguments:
        ----------
        lightcurves -- list of astropy Tables containing the lightcurves
                       The meta dictionary of each table must contain the fixed
                       fit parameters, e.g. redshift

        Parameters:
        -----------
        model -- sncosmo model to be fit; if it has no MW dust correction,
                 the effect will be added

        source -- sncosmo source for the fit, ignored if model is not None; 
                  can be string for built-in model

        fit_param -- list of parameters to be fit
        """
        if model is None:
            self.set_source(source)
        else:
            self.set_model(model)

        self.set_cosmo(cosmo, update=False)

        self._properties['lightcurves'] = lightcurves
        self._side_properties['fit_param'] = fit_param

        self._update_()

    def _update_(self):
        """
        """
        self.fit_all()
        self._update_HR()

    def _update_HR_(self):
        """
        Not implemented yet, requires calculation of peak magnitudes
        """
        pass

    def fit_all(self):
        """
        """
        if self.raw_fit is not None:
            self._derived_properties["raw_fit"] = []
            
        for lc in self.lightcurves:
            for pname in [name for name in self.model.param_names
                          if name not in self.fit_param]:
                if pname not in lc.meta.keys():
                    raise KeyError("Parameter '%s' not in lightcurve meta dict"%pname)
                self.model.set(**{pname: lc.meta[pname]})

            res, _ = sncosmo.fit_lc(lc, self.model, self.fit_param)

            self._derived_keys['raw_fit'].append(res)
                
    # ------------
    # - Properties
    @property
    def model(self):
        """Light curve model (derived from source if not set)"""
        if self._properties["model"] is None:
            self.set_model(sncosmo.Model(source=self.source))
        
        return self._properties["model"]

    def set_model(self, model):
        """
        Set the transient model.
        If it does not have MW dust effect, the effect is added.
        """
        if model.__class__ is not sncosmo.models.Model:
            raise TypeError("model must be sncosmo.model.Model")

        if "mwebv" not in model.param_names:
            model.add_effect(sncosmo.CCM89Dust(), 'mw', 'obs')

        self._side_properties["model"] = model

    @property
    def cosmo(self):
        """Astropy cosmology object used to derive Hubble Residuals"""
        return self._properties["cosmo"]    

    def set_cosmo(self, cosmo, update=True):
        """
        """
        # TODO: Add check that cosmo is an astropy cosmology object
        self._properties["cosmo"] = cosmo
        
        if update:
            self._update_HR_()

    @property
    def lightcurves(self):
        """List of Astropy tables containing the lightcurves"""
        return self._properties["lightcurves"]    

    # -----------------
    # - Side properties
    @property
    def source(self):
        """Lightcurve model source (only needed if model not provided)"""        
        return self._side_properties["source"]

    def set_source(self, source):
        """Set Lightcurve model source"""        
        self._side_properties["source"] = source

    @property
    def fit_param(self):
        """Lightcurve fit parameter names"""        
        return self._side_properties["fit_param"]
    # -----------------
    # - Derived properties
    @property
    def raw_fit(self):
        """Raw fit ouput of sncosmo.fit_lc()"""        
        return self._derived_properties["raw_fit"]

    @property
    def fit(self):
        """Summary of fit results (not implement yet)"""        
        pass

     @property
     def true(self):
        """True values for simulated lcs taken from meta dict of lightcurves"""
        return self._derived_properties["raw_fit"]

 
