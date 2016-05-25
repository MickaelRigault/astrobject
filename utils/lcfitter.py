#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes for lightcurve fits and their statistics

[Currently not using the sncosmo model covariance option when fitting.
 This will lead to underestimated fit uncertainties.]

TODO:
- Calculate band mags and unertainties
- Calculate HRs
- Make SNFitter that can also correct for stretch and colour (given alpha/beta)
"""

import numpy as np
import cPickle
import glob
import re
import os
from itertools import combinations
from collections import OrderedDict as odict

import sncosmo
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, Column

from .. import BaseObject

__all__ = ["LCFitter"]

class LCFitter( BaseObject ):
    """ Basic lightcurve fitter class """
        
    PROPERTIES         = ['model', 'cosmo', 'lightcurves']
    SIDE_PROPERTIES    = ['source', 'fit_param', 'bands']
    DERIVED_PROPERTIES = ['fit', 'raw_fit', 'idx_good', 'good_lightcurves', 
                          'true']

    def __init__(self, lightcurves=None, source='salt2', model=None, 
                 cosmo=FlatLambdaCDM(H0=70, Om0=0.3), empty=False,
                 fit_param=['t0', 'x0', 'x1', 'c'], load_from=None, 
                 modelcov=False):
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

        load_from -- file with previously saved lcs and fits

        modelcov -- Use model covariance (currently not properly implemented; 
                    using hacked sncosmo.fitting) [Not implemented yet]
        """
        self.__build__()
        if empty:
            return

        if model is None:
            self.set_source(source)
        else:
            self.set_model(model)

        self.set_cosmo(cosmo, update=False)

        if load_from is None:
            if lightcurves is None:
                raise ValueError('Please provide lightcurves or savefile')
            self._properties['lightcurves'] = lightcurves
            self._side_properties['fit_param'] = fit_param
            self._update_(modelcov=modelcov)
        else:
            self.load(load_from)
        
    # --------------------------- #
    # - Update and fit Methods  - #
    # --------------------------- #
    def _update_(self, modelcov=False):
        """
        """
        self.fit_all(modelcov=modelcov)

    def fit_all(self, modelcov=False):
        """
        """
        if modelcov:
            # TODO: Implement fitting with modelcov
            # from hacked_sncosmo_fitting import fit_lc_hacked
            raise NotImplementedError("Model covariance coming soon.")

        param0 = self._get_current_model_param_()
        self._derived_properties["raw_fit"] = []
        self._derived_properties["idx_good"] = []

        for k, lc in enumerate(self.lightcurves):
            for pname in [name for name in self.model.param_names
                          if name not in self.fit_param and name != 'mwr_v']:
                if pname not in lc.meta.keys():
                    raise KeyError("Parameter '%s' not in lightcurve meta dict"%pname)
                self.model.set(**{pname: lc.meta[pname]})

            try: 
                if modelcov:
                    res, _ = fit_lc_hacked(lc, self.model, self.fit_param)
                else:
                    res, _ = sncosmo.fit_lc(lc, self.model, self.fit_param)
                
                if res['covariance'] is not None:
                    self._derived_properties['raw_fit'].append(res)
                    self._derived_properties['idx_good'].append(k)
                else:
                    print "Light curve fit #%i failed"%k
            except sncosmo.fitting.DataQualityError:
                print "Light curve fit #%i failed because of data quality"%k
            except RuntimeError:
                print "Light curve fit #%i failed to converge"%k

        self.model.set(**param0)

    def _make_param_table_(self):
        """
        """
        self._derived_properties["fit"] = None
                
        new = Table()
        new.add_column(Column(name='n_points', data=self.get_param('n_points')))
        new.add_column(Column(name='chisq', data=self.get_param('chisq')))
        for param in [name for name in self.model.param_names
                      if name not in self.fit_param and name != 'mwr_v']:
            new.add_column(Column(name=param, data=self.get_param(param)))
        
        full_list = self.fit_param+self._get_mag_names_()
        for param in full_list:
            new.add_column(Column(name=param, data=self.get_param(param)))
            new.add_column(Column(name='err_%s'%param, 
                           data=self.get_param('err_%s'%param)))

        for n1, n2 in combinations(full_list, 2):
            new.add_column(Column(name='cov_%s_%s'%(n1,n2), 
                           data=self.get_param('cov_%s_%s'%(n1,n2))))

        self._derived_properties["fit"] = new

    # --------------------------- #
    # - Save/load Methods       - #
    # --------------------------- #

    def save(self, savefile):
        """
        Save lightcurves and fits to a Pickle file.
        """
        out = {'lightcurves': self.lightcurves,
               'fit_param': self.fit_param,
               'raw_fit': self.raw_fit,
               'idx_good': self.idx_good}
            
        cPickle.dump(out, open(savefile, 'w'))

    def load(self, savefile):
        """
        """
        data = cPickle.load(open(savefile))

        self._properties['lightcurves'] = data['lightcurves']
        self._side_properties['fit_param'] = data['fit_param']
        self._derived_properties['raw_fit'] = data['raw_fit']
        self._derived_properties['idx_good'] = data['idx_good'] 
        
        if self._derived_properties['fit'] is not None:
            self._derived_properties['fit'] = None
            
    # --------------------------- #
    # - Get Methods             - #
    # --------------------------- #
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
            return np.array(self._derived_properties["fit"][param])
        
        if param == 'n_points':
            return np.array([len(self.lightcurves[k]) for k in self.idx_good])
        elif param == 'chisq':
            return np.array([res['chisq'] for res in self.raw_fit])
        elif param in self.model.param_names:
            k = self.model.param_names.index(param)
            return np.array([res['parameters'][k] for res in self.raw_fit])
        elif param.startswith('mag_'):
            if param[4:] not in self.bands.keys():
                raise ValueError('Unknown band: %s'%param[4:])
            return self.get_bandmag(**self.bands[param[4:]])
        elif param.startswith('err_'):
            if not param[4:].startswith('mag'):
                return np.array([res['errors'][param[4:]] 
                                 for res in self.raw_fit])
            else:
                return self.get_bandmag_err(**self.bands[param[8:]])
        elif param.startswith('cov_'):
            return self._get_cov_(param)
        else:
            raise ValueError('Unknown parameter name: %s'%param)

    def get_param_true(self, param):
        """
        Get true (i.e. input) values of the parameter param
        """
        if param.startswith('mag'):
            return self.get_bandmag(self.bands[param[4:]], true_param=True)
        else:
            return np.array([lc.meta[param] for lc in self.good_lightcurves])

    def get_hr(self, band, absmag=-19.3, corr=None):
        """
        band must be registered band 
        corr is dictionary of lc parameter names for correction and
        nuisance parameter values, e.g. {'x1': 0.13, 'c': -3}
        (Note the sign change for beta. This is because the correction is always
        performed by adding the values)
        """
        if corr is None:
            corr = {}

        mag = self.get_param('mag_%s'%band)
        distmod = self.cosmo.distmod(self.get_param('z')).value
        hr = mag - absmag - distmod

        for name, value in corr.items():
            hr += value * self.get_param(name)

        return hr

    def get_hr_err(self, band, corr=None):
        """
        band must be registered band 
        corr is dictionary of lc parameter names for correction and
        nuisance parameter values, e.g. {'x1': 0.13, 'c': -3}
        (Note the sign change for beta. This is because the correction is always
        performed by adding the values)
        """
        if corr is None:
            corr = {}

        var = self.get_param('err_mag_%s'%band)**2

        for name, value in corr.items():
            var += value**2 * self.get_param('err_%s'%name)**2
            var += 2*value * self.get_param('cov_mag_%s_%s'%(band,name))

        for n0, n1 in combinations(corr.keys(), 2):
            var +- 2*corr[n0]*corr[n1] * self.get_param('cov_%s_%s'%(n0,n1)) 

        return np.sqrt(var)

    def _get_cov_idx_(self, param, return_names=False):
        """
        Find indices for covariance name
        Can also return the parameter names instead

        TODO: Figure out how this should handle peak magnitudes 
        """
        full_list = self.fit_param+self._get_mag_names_()
        pattern = '(%s)'%('|'.join(full_list))
        matches = re.findall(pattern, param[4:])
        
        if len(matches) != 2:
            raise ValueError('Unknown covariance name: %s'%param)
        elif matches[0] == matches[1]:
            raise ValueError('Use "err_%s" to get variance'%matches[0])

        if return_names:
            return tuple(matches)

        return (full_list.index(matches[0]), full_list.index(matches[1]))

    def _get_cov_reverse_(self, param):
        """
        Reverse order of parameter names in covariance name
        """
        n1, n2 = self._get_cov_idx_(param, return_names=True)
        
        return 'cov_%s_%s'%(n2, n1)

    def _get_cov_(self, param):
        """
        Get array of covariances
        TODO: Extend to calculating the band covariances as well
        """
        k0, k1 = self._get_cov_idx_(param)

        if k0 < len(self.fit_param) and k1 < len(self.fit_param):
            return np.array([(res['covariance'][k0,k1] 
                              if res['covariance'] is not None
                              else np.nan) 
                             for res in self.raw_fit])
        else:
            grad0 = self._get_param_grad_(k0)
            grad1 = self._get_param_grad_(k1)
            return np.array([g0.dot(res['covariance']).dot(g1) 
                             for g0, g1, res 
                             in zip(grad0, grad1, self.raw_fit)])

    def _get_param_grad_(self, idx):
        """
        """
        if idx < len(self.fit_param):
            out = np.zeros(len(self.fit_param))
            out[idx] = 1.
            return [out for k in range(len(self.raw_fit))]
        else:
            idx -= len(self.fit_param) 
            return self._get_bandmag_gradient_(**self.bands.values()[idx])

    def _get_param_dicts_(self, true_param=False):
        """
        Returns odered dictionaries of fixed parameters, fit parameters and
        their uncertainties (all that is needed to get bandmags + uncertainties)
        """
        out = []
        for lc, res in zip(self.good_lightcurves, self.raw_fit):
            fixed = {p: lc.meta[p] for p in [n for n in self.model.param_names
                                             if n not in self.fit_param 
                                             and n != 'mwr_v']}
            fit = odict()
            for p in self.fit_param:
                if true_param:
                    fit[p] = lc.meta[p]
                else:
                    fit[p] = res['parameters'][res['param_names'].index(p)]

            out.append((fixed, fit, res['errors']))

        return out

    def get_bandmag(self, band='bessellb', magsys='vega', t=0, restframe=True,
                    remove_mw_dust=True, true_param=False):
        """
        Returns the magnitudes of transient according to lightcurve parameters
        [Modified from simultarget, should find a common place to store this]
        """
        # Save old params, so you can restore them 
        param0 = self._get_current_model_param_()
 
        out = []
        for fixed, fit, err  in self._get_param_dicts_(true_param=true_param):
            self.model.set(**fixed)
            self.model.set(**fit)

            if restframe:
                self.model.set(z=0)
            if remove_mw_dust:
                self.model.set(mwebv=0)

            out.append(self.model.bandmag(band, magsys, fit['t0'] + t))

        self.model.set(**param0)

        return np.array(out)

    def get_bandmag_err(self, **kwargs):
        """
        kwargs passed to _get_bandmag_gradient_
        """
        grad = self._get_bandmag_gradient_(**kwargs)

        return np.sqrt(np.array([g.dot(res['covariance']).dot(g) 
                                 for g, res in zip(grad, self.raw_fit)]))

    def _get_bandmag_gradient_(self, band='bessellb', magsys='vega', t=0, 
                               restframe=True, remove_mw_dust=True,
                               true_param=False):
        """
        Return gradient of _get_bandmag as function of param
        param, sig must be dictionaries of means and uncertainties
        Best use odicts to make sure that the order of the components is correct
        """
        # Save old params, so you can restore them 
        param0 = self._get_current_model_param_()

        out = []

        for fixed, fit, err  in self._get_param_dicts_(true_param=true_param):
            grad = []

            self.model.set(**fixed)
            self.model.set(**fit)

            if restframe:
                self.model.set(z=0)
            if remove_mw_dust:
                self.model.set(mwebv=0)

            for key,val in fit.items():
                if key != 't0':
                    self.model.set(**fit)
                    h = err[key] / 100.

                    self.model.set(**{key: val - h})
                    m0 = self.model.bandmag(band, magsys, fit['t0'] + t)
                    
                    self.model.set(**{key: val + h})
                    m1 = self.model.bandmag(band, magsys, fit['t0'] + t)

                    grad.append((m1 - m0) / (2. * h))
                else:
                    grad.append(0.)

            out.append(np.array(grad))

        self.model.set(**param0)

        return out

    def _get_current_model_param_(self):
        """
        Get current model parameter, so it can be reset to them
        """
        return {name: value for name, value 
                in zip(self.model.param_names,
                       self.model.parameters)}

    def _get_mag_names_(self):
        """
        """
        return ['mag_%s'%name for name in self.bands.keys()]

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
        
    @property
    def bands(self):
        """
        Bandpasses for which to calculcate peak magnitudes
        Ordered Dictionary containing the sncosmo bandpass names or objects 
        and the magnitude system 
        Keys are short names for the bands, e.g {'B': {'band': 'bessellb', 
                                                       'magsys': 'vega']}
        """        
        if self._side_properties["bands"] is None:
            return {}
        return self._side_properties["bands"]

    def add_band(self, name, band, magsys='vega', update=True):
        """
        Add a bandpass to Fitter.band

        name -- short name for the band (used for Fitter.get_param etc.)
        band -- sncosmo bandpass object or name of registered band

        Options:
        update -- update Fitter.fit table if it exists (This will recreate the 
                  *whole* table, so it could be slow. But it is necessary for 
                  subsequent calls of Fitter.get_param)
        """
        if name in self.bands.keys():
            raise ValueError('Name alreay in use. Please select other name '
                             + 'or delete the old band first.')
        
        if self._side_properties["bands"] is None:
            self._side_properties["bands"] = odict()

        self._side_properties["bands"][name] = {'band': band, 
                                                'magsys': magsys}

        if update and self._derived_properties["fit"] is not None:
            self._derived_properties["fit"] = None
            self._make_param_table_()

    def remove_band(self, name, update=True):
        """
        Add a bandpass to Fitter.band

        name -- short name for the band (used for Fitter.get_param etc.)

        Options:
        update -- update Fitter.fit table if it exists (This will recreate the 
                  *whole* table, so it could be slow.)
        """
        if name not in self.bands.keys():
            raise ValueError('Unknown name.')
        
        del self._side_properties["bands"][name]

        if update and self._derived_properties["fit"] is not None:
            self._derived_properties["fit"] = None
            self._make_param_table_()

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
            self._make_param_table_()

        return self._derived_properties["fit"]

    @property
    def idx_good(self):
        """Indices of good fits"""
        return self._derived_properties["idx_good"]

    @property
    def good_lightcurves(self):
        """List of lightcurves that could be fit"""
        return [self.lightcurves[k] for k in self.idx_good]
