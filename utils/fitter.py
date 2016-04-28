#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes for lightcurve fits and their statistics

[Currently not using the sncosmo model covariance option when fitting.
 This will lead to underestimated fit uncertainties.]
"""

import numpy as np
import re
from itertools import combinations

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

    _side_properties_keys = ['source', 'fit_param']

    _derived_properties_keys = ['fit', 'raw_fit', 'true']

    def __init__(self, lightcurves, source='salt2', model=None, 
                 cosmo=FlatLambdaCDM(H0=70, Om0=0.3), empty=False,
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
        self.__build__()
        if empty:
            return

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
        self._update_covnames_()
        self._update_HR_()

    def _update_HR_(self):
        """
        Not implemented yet, requires calculation of peak magnitudes
        """
        pass

    def fit_all(self):
        """
        """
        if self.raw_fit is None:
            self._derived_properties["raw_fit"] = []
            
        for lc in self.lightcurves:
            for pname in [name for name in self.model.param_names
                          if name not in self.fit_param and name != 'mwr_v']:
                if pname not in lc.meta.keys():
                    raise KeyError("Parameter '%s' not in lightcurve meta dict"%pname)
                self.model.set(**{pname: lc.meta[pname]})

            res, _ = sncosmo.fit_lc(lc, self.model, self.fit_param)

            self._derived_properties['raw_fit'].append(res)

    def get_param(self, param):
        """
        Return array of values for parameter 'param' 
        """
        if self._derived_properties["fit"] is not None:
            if param not in self.fit.colnames:
                if param.startswith('cov_'):
                    param = self._get_cov_reverse_(param)
                else:
                    raise ValueError('Unknown parameter name')                
            return self._derived_properties["fit"][param]

        if param in self.model.param_names:
            k = self.model.param_names.index(param)
            return np.array([res['parameters'][k] for res in self.raw_fit])
        elif param.startswith('mag_'):
            raise ImplementationError('Peak magnitudes not implemented yet')
        elif param.startswith('err_'):
            if not param[4:].startswith('mag'):
                return np.array([res['errors'][param[4:]] 
                                 for res in self.raw_fit])
            else:
                raise ImplementationError('Peak magnitudes not implemented yet')
        elif param.startswith('cov_'):
            return self._get_cov_(param)
        else:
            raise ValueError('Unknown parameter name: %s'%param)

    def _make_param_table_(self):
        """
        """
        self._derived_properties["fit"] = None
                
        new = Table()
        for param in [name for name in self.model.param_names
                      if name not in self.fit_param and name != 'mwr_v']:
            new.add_column(name=param, data=self.get_param(param))
        for param in self.fit_param:
            new.add_column(name=param, data=self.get_param(param))
            new.add_column(name='err_%s'%param, 
                           data=self.get_param('err_%s'%param))
        for n1, n2 in combinations(self.fit_param, 2):
            new.add_column(name='cov_%s_%s'%(n1,n2), 
                           data=self.get_param('cov_%s_%s'%(n1,n2)))

    def _get_cov_idx_(self, param, return_names=False):
        """
        Find indices for covariance name
        Can also return the parameter names instead

        TODO: Figure out how this should handle peak magnitudes 
        """
        pattern = '(%s)'%('|'.join(self.fit_param))
        matches = re.findall(patter, param[4:])
        
        if len(matches) != 2:
            raise ValueError('Unknown covariance name: %s'%param)
        
        if return_names:
            return tuple(matches)

        return (self.fit_param.index(matches[0]),
                self.fit_param.index(matches[1]))

    def _get_cov_reverse_(self, param):
        """
        Reverse order of parameter names in covariance name
        """
        n1, n2 = self._get_cov_idx_(param, return_names=True)
        
        return 'cov_%s_%s'%(n1, n2)

    def _get_cov_(self, param):
        """
        Get array of covariances
        TODO: Extend to calculating the 
        """
        k0, k1 = self._get_cov_idx_(param)
        
        return np.array([res['covariance'][k0,k1] for res in self.raw_fit])

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

        self._properties["model"] = model

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
        if self._derived_properties["raw_fit"] is None:
            raise ValueError("Lightcurves not fit yet. Run Fitter.fit_all().")

        return self._derived_properties["raw_fit"]

    @property
    def fit(self):
        """Summary of fit results (not implement yet)"""        
        if self._derived_properties["fit"] is None:
            self._process_raw_fits_()

        return self._derived_properties["fit"]

    @property
    def true(self):
        """True values for simulated lcs taken from meta dict of lightcurves"""
        return self._derived_properties["raw_fit"]

 
