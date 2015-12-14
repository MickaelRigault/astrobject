#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This pupose of this method is the have a generator of fake astrotarget, like SN"""

import numpy as np
from numpy.random import uniform, normal
import sncosmo
import warnings

from ..astrobject.baseobject import BaseObject,astrotarget
from ..astrobject.transient import transient

from ..utils.tools import kwargs_extract,kwargs_update
from ..utils import random

_d2r = np.pi / 180

__all__ = ["transient_generator","generate_transients"]

def transient_generator(zrange,ratekind="basic",
                        ra_range=[-180,180], dec_range=[-90,90],
                        ntransients=None,
                        **kwargs):
    """
    This model return the object that enable to create and change
    the kind of transient you which to set in the sky.

    # - HERE COPY PASTE THE TransientGenerator INIT - #
    # - TO BE DONE
    
    """
    return TransientGenerator(ratekind=ratekind,ntransients=ntransients,
                               zrange=zrange,
                              ra_range=ra_range,dec_range=dec_range,
                              **kwargs)
    
def generate_transients(zrange,**kwargs):
    """This module call transient_generator to create the
    TransientGenerator object and then returns the associated
    TransientGenerator.transients

    # - HERE COPY PASTE the transient_generator docstring
    """
    return transient_generator(zrange,**kwargs).transients


#######################################
#                                     #
# Generator: Any Transient            #
#                                     #
#######################################
class TransientGenerator( BaseObject ):
    """
    """
    __nature__ = "TransientGenerator"
    
    _properties_keys = ["transient_coverage",
                        "event_coverage"]
        
    _derived_properties_keys = ["ratefunc","simul_parameters"]

    

    def __init__(self,zrange=[0.0,0.2], ratekind="basic", # How deep
                 mdj_range=[57754.0,58849.0],
                 ra_range=(-180,180),dec_range=(-90,90), # Where, see also kwargs
                 ntransients=None, empty=False,**kwargs):
        """
        """
        self.__build__()
        if empty:
            return
        
        self.create(zrange,
                    ratekind=ratekind, ntransients=ntransients,
                    ra_range=ra_range, dec_range=dec_range,
                    mdj_range=mdj_range,
                    **kwargs)

    def create(self,zrange,ratekind="basic",ntransients=None,
               mdj_range=[57754.0,58849.0],
               ra_range=(-180,180),dec_range=(-90,90),
               mw_exclusion=0,transientpror={}):
        """
        """
        # == Add the Input Test == #
        #   TO BE DONE
        
        # *************** #
        # * Create      * #
        # *************** #
        # -- This will be directly used as random.radec inputs
        self.set_event_parameters(update=False,
                                  **{"ra_range":ra_range,"dec_range":dec_range,
                                   "zcmb_range":zrange,"mdj_range":mdj_range,
                                   "mw_exclusion":mw_exclusion})
        
        self.set_transient_parameters(ratekind=ratekind,
                                      update=False,**transientpror)
        
        self._update_()
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    # --------------------------- #
    # - Set Methods             - #
    # --------------------------- #
    def set_event_parameters(self,update=True,**kwargs):
        """
        Change the properties associated to the transient events.
        Known properties: "ra_range","dec_range","zcmb_range","mdj_range",
                          "mw_exclusion"

        Set update to True to update the derived properties
        """
        known_event_prop = ["ra_range","dec_range","zcmb_range",
                            "mw_exclusion","mdj_range"]
            
        for k in kwargs.keys():
            if k not in known_event_prop:
                raise ValueError("'%s' is not an 'event' parameter."%k +\
                                 " These are: "+", ".join(known_event_prop))
                                 
        self._properties["event_coverage"] = kwargs_update(self.event_coverage,
                                                           **kwargs)
        if update:
            self._update_()
                            
    def set_transient_parameters(self,ratekind="basic",
                                 update=True,type_=None):
        """
        This method will define the transient properties.
        """
        if self._properties["transient_coverage"] is None:
            self._properties["transient_coverage"] = {}
            
        # -- if this works, you are good to go
        f = RateGenerator().get_ratefunc(transient=type_,ratekind=ratekind)
        
        # - you are good to fill it
        self._properties["transient_coverage"]["transienttype"] = type_
        self._properties["transient_coverage"]["ratekind"] = ratekind
        self._derived_properties['ratefunction'] = f
        
        if update:
            self._update_()
            
    # --------------------------- #
    # - Get Methods             - #
    # --------------------------- #
    def get_transients(self,index=None):
        """loops over the transientsources to load the transients objects.
        This method could be a bit slow..."""
        return [transient(**s) for s in self.get_transientsource(index)]
    
    def get_transientsource(self,index=None):
        """dictionary containing the fundamental parameters that enable to
        load the transient objects"""
        if index is not None and "__iter__" not in dir(index):
            index = [index]
        
        return [dict(name="simul%d"%i,ra=self.ra[i],dec=self.dec[i], zcmb=self.zcmb[i],
                mjd=self.mjd[i],
                lightcurve=self.lightcurve[i] if self.has_lightcurves() else None)
                for i in range(self.ntransient) if i is None or i in index]
            
    # --------------------------- #
    # - Plots Methods           - #
    # --------------------------- #
    def show_skycoverage(self, ax=None, savefile=None, show=True, cscale=None, 
                         cblabel=None, **kwargs):
        """This function enable to draw on the sky the position of the
        transients"""
        import matplotlib.pyplot as mpl
        from ..utils.mpladdon import figout, skyplot
        from ..utils.plot.skyplot import ax_skyplot
        self._plot = {}

        # ------------------
        # - Color Scale 
        if cscale == 'zcmb':
            c = np.asarray([t['zcmb'] for t in self.transientsources])
            if cblabel is None:
                cblabel = r"$\mathrm{Redshift}$"
        elif type(cscale) != str and hasattr(cscale, '__iter__'):
            c = cscale
        elif cscale is not None:
            raise ValueError('cscale must be array or predefined string, '+\
                             ' e.g. "redshift"')
        
        # ------------------
        # - Axis definition
        if ax is None:
            ax_default = dict(fig=None, figsize=(12, 6), 
                                   rect=[0.1, 0.1, 0.8, 0.8], 
                                   projection='mollweide')
            if cscale is not None:
                ax_default['figsize'] = (12,8)
                
            ax_kw, kwargs = kwargs_extract(ax_default, **kwargs)
            fig, ax = ax_skyplot(**ax_kw)
        elif ("MollweideTransform" not in dir(ax) and
              "HammerTransform" not in dir(ax)):
            raise TypeError("The given 'ax' most likely is not a matplotlib axis "+\
                            "with Mollweide or Hammer projection. Transform "+\
                            "function not found.")
        else:
            fig = ax.fig

        # ------------------
        # - Actual plotting
        if cscale is None:
            pl = ax.skyplot(self.ra, self.dec, **kwargs)
            cb = None
        else:
            pl = ax.skyscatter(self.ra, self.dec, c=c, **kwargs)
            cb = fig.colorbar(pl, orientation='horizontal', shrink=0.85, pad=0.08)
            if cblabel is not None:
                cb.set_label(cblabel, fontsize="x-large") 

        # ------------------- #
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["plot"]   = pl
        if cb is not None:
            self._plot["cbar"] = cb

        fig.figout(savefile=savefile,show=show)
        
        return self._plot        

    def hist_skycoverage(self, ax=None, savefile=None, show=True, 
                         cblabel=r"$N_{SNe}$", **kwargs):
        """This function draws a sky histogram of the transient coverage"""
        import matplotlib.pyplot as mpl
        from ..utils.mpladdon import figout, skyhist
        from ..utils.plot.skyplot import ax_skyplot
        self._plot = {}

        if ax is None:
            ax_default = dict(fig=None, figsize=(12, 8),
                              rect=[0.1, 0.1, 0.8, 0.8],
                              projection='mollweide')
            ax_kw, kwargs = kwargs_extract(ax_default, **kwargs)
            fig, ax = ax_skyplot(**ax_kw)
        elif ("MollweideTransform" not in dir(ax) and
              "HammerTransform" not in dir(ax)):
            raise TypeError("The given 'ax' most likely is not a matplotlib axis "+\
                        "with Mollweide or Hammer projection. Transform "+\
                        "function not found.")
        else:
            fig = ax.fig

        # maybe these arrays can be integrate into the generator
        p = ax.skyhist(self.ra, self.dec, **kwargs)
        cb = fig.colorbar(p, orientation='horizontal', shrink=0.85, pad=0.08)
        if cblabel:
            cb.set_label(cblabel, fontsize="x-large") 

        # ------------------- #
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["patches"] = p
        if cb is not None:
            self._plot["cbar"] = cb

        fig.figout(savefile=savefile,show=show)
        return self._plot 
    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _update_simulation_(self):
        """
        """
        # -----------------------
        # - Redshift from Rate
        self.simul_parameters["zcmb"] = \
          list(sncosmo.zdist(self.zcmb_range[0], self.zcmb_range[1],
                             time=self.timescale, area=self.coveredarea,
                             ratefunc=self.ratefunc))

        self.simul_parameters["mjd"] = self._simulate_mjd_()
        self.simul_parameters["ra"], self.simul_parameters["dec"] = \
          random.radec(self.ntransient,
                       ra_range=self.ra_range,
                       dec_range=self.dec_range,
                       mw_exclusion=self._get_event_property_("mw_exclusion"))
        
    def _update_(self):
        """This module create the derived values based on the
        fundamental ones"""
        # --------------
        # - update the actual simulation
        self._update_simulation_()        

        
    def _simulate_mjd_(self):
        """
        Be default, this is a random flat time distribution returning a float
        per transient.
        Simple overwrite this function in a child-class to have more advanced
        properties.
        The Simulated mjd will be stored in transientsources.
        """
        return np.random.rand(self.ntransient)*self.timescale + self.mjd_range[0]

    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    def _get_event_property_(self,key):
        """
        """
        if self.transient_coverage is None:
            raise AttributeError("'transient_coverage' has not been defined")
        
        return self.event_coverage[key]
    # -----------------
    # - Ranges
    @property
    def zcmb_range(self):
        """zcmb range used to draw transient"""
        return self._get_event_property_("zcmb_range")
    
    @property
    def ra_range(self):
        return self._get_event_property_("ra_range")
    
    @property
    def dec_range(self):
        return self._get_event_property_("dec_range")
    
    @property
    def mjd_range(self):
        """zcmb range used to draw transient"""
        return self._get_event_property_("mdj_range")
    # -----------------
    # - Rates
    @property
    def ratefunc(self):
        return self._derived_properties["ratefunction"]

    # -------------------------------
    # - Derived Transient Properties
    @property
    def timescale(self):
        return self.mjd_range[1] - self.mjd_range[0]
       
    @property
    def ntransient(self):
        """number of transient requested"""
        return len(self.zcmb)

    @property
    def coveredarea(self):
        """Covered area in degree squared"""
        mw_exclusion = self._get_event_property_("mw_exclusion")
        
        # Area in steradians without accounting for MW exclusion
        area_sr = ((np.sin(self.dec_range[1] * _d2r) 
                    - np.sin(self.dec_range[0] * _d2r)) 
                   * (self.ra_range[1] - self.ra_range[0]) * _d2r)

        if mw_exclusion > 0:
            if self.ra_range == [-180, 180] and self.dec_range == [-90, 90]:
                area_sr -= 4 * np.pi * np.sin(mw_exclusion * _d2r) 
            else:
                # Make sure the warning is issued every time
                warnings.simplefilter('always', UserWarning)
                warnings.warn("MW exclusion was ignored when calculating covered area.")

        return area_sr / _d2r ** 2
        
    # --------------------
    # - Target Coverage
    @property
    def transient_coverage(self):
        if self._properties["transient_coverage"] is None:
            self._properties["transient_coverage"] = {}
        return self._properties["transient_coverage"]
        
    @property
    def event_coverage(self):
        """where and when the transient could be drawn"""
        if self._properties["event_coverage"] is None:
            self._properties["event_coverage"] = {}
            
        return self._properties["event_coverage"]

    # --------------------
    # - Derived Properties
    @property
    def simul_parameters(self):
        if self._derived_properties["simul_parameters"] is None:
            self._derived_properties["simul_parameters"] = {}
        return self._derived_properties["simul_parameters"]

    # -----------------
    # - transient info
    @property
    def zcmb(self):
        """Simulated zcmb based on the given rate function"""
        return self.simul_parameters["zcmb"]

    @property
    def mjd(self):
        """Loop over the transient sources to get the mjd"""
        return np.asarray(self.simul_parameters['mjd'])
    
    @property
    def ra(self):
        """Loop over the transient sources to get the ra"""
        return np.asarray(self.simul_parameters['ra'])
    
    @property
    def dec(self):
        """Loop over the transient sources to get the ra"""        
        return np.asarray(self.simul_parameters['dec'])
    
    @property
    def lightcurve(self):
        """Check if there is a lightcurve simul_parameter defined"""
        return None if "lightcurve" not in self.simul_parameters.keys()\
          else np.asarray(self.simul_parameters['lightcurve'])
          
    def has_lightcurves(self):
        return False if self.lightcurve is None \
          else True
            
    
#######################################
#                                     #
# Generator: SN Ia                    #
#                                     #
#######################################
class SNIaGenerator( TransientGenerator ):
    """
    This child Class enable to add the SN properties to the
    transient generator.
    """
    
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    # -------------------- #
    # - Hacked Methods   - #
    # -------------------- #
    
    def set_transient_parameters(self,ratekind="basic",
                                 lcsimulation="basic",
                                 lcmodel="salt2",type_="Ia",
                                 update=True,lcsimul_prop={}):
        """
        Add to the TransientGenerator the SN Ia properties
        """
        # - If this works, you are good to go
        f = LightCurveGenerator().get_lightcurve_func(transient=type_,simulation=lcsimulation)

        super(SNIaGenerator,self).set_transient_parameters(type_=type_,ratekind=ratekind,update=False)
        print "loading lightcurve parameters..."
        self._properties["transient_coverage"] = {
            "lcsimulation":lcsimulation,
            "lcmodel":lcmodel,
            "lcprop_function":f
            }
        
        if update:
            self._update_()

    # ----------------- #
    # - Set Methods   - #
    # ----------------- #
        
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
          
#######################################
#                                     #
# Generator: SN Rate                  #
#                                     #
#######################################
class _PropertyGenerator_(BaseObject):
    """
    """
    __nature__ = "PropertyGenerator"

    # ========================== #
    # = Internal               = #
    # ========================== #
    def _parse_rate_(self,key="rate"):
        return [m.split(key+"_")[-1] for m in dir(self)
                if m.startswith(key)]
    
class RateGenerator( _PropertyGenerator_ ):
    """
    This follows the SN cosmo ratefunc
    
    ratefunc : callable
        A callable that accepts a single float (redshift) and returns the
        comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.
        The default is a function that returns ``1.e-4``.
    
    """
    __nature__ = "RateGenerator"
    #
    # User: Add any rate_TTT_BLA(self,z) to define a new rate for the TTT transient.
    #       where this return the comoving volumetric rate at each redshift in
    #        units of yr^-1 Mpc^-3.
    # 
    #       This rate will be accessible from
    #       get_ratefunc(transient="TTT",ratekind="BLA")
    #

    def get_ratefunc(self,transient="Ia",ratekind="basic"):
        """
        Parameters
        ----------

        Return
        ------
        function(z)
        """
        # ---------------
        # - Transients
        if transient is None or transient == "":
            avialable_rates = self.known_rates
            transient = None
        elif transient == "Ia":
            avialable_rates = self.known_Ia_rates
        else:
            raise ValueError("'%s' is not a known transient"%transient)
    
        # ---------------
        # - Rate Kinds
        if ratekind not in avialable_rates:
            raise ValueError("not '%s' rate kind for '%s'"%(ratekind,transient)+\
                             "These are: "+",".join(avialable_rates))
        # -- the output
        if transient is None:
            return eval("self.rate_%s"%(ratekind))
        return eval("self.rate_%s_%s"%(transient,ratekind))
    
    # ========================== #
    # = Rates                  = #
    # ========================== #
    def rate_basic(self,z):
        """
        Basic default rate function in sncosmo: returns ``1.e-4``.
        (comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.)
        """
        return 1.e-4
    # ----------------- #
    # - Ia rate       - #
    # ----------------- #
    def rate_Ia_basic(self,z):
        """
        Basic default rate function in sncosmo: returns ``1.e-4``.
        (comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.)
        """
        return self.rate_basic(z)
    
    def rate_Ia_ptflike(self,z):
        """
        """
        return self.rate_basic(z)
    # ========================== #
    # = *** Rates              = #
    # ========================== #
    

    
    # ========================== #
    # = Properties             = #
    # ========================== #
    @property
    def known_rates(self):
        return self._parse_rate_("rate")
    
    @property
    def known_Ia_rates(self):
        return self._parse_rate_("rate_Ia")

class LightCurveGenerator( _PropertyGenerator_ ):
    """
    """

    # ========================== #
    # = Method                 = #
    # ========================== #
    def get_lightcurve_func(self,transient="Ia",simulation="basic"):
        """
        Parameters
        ----------

        Return
        ------
        a lightcurve generator function
        """
        # ---------------
        # - Transients
        if transient is None or transient == "":
            avialable_lc = self.known_lightcurve_simulation
            if len(avialable_lc) == 0:
                raise NotImplementedError("no lightcurve simulation implemented")
            
            transient = None
            
        elif transient == "Ia":
            avialable_lc = self.known_Ia_lightcurve_simulation
        else:
            raise ValueError("'%s' is not a known transient"%transient)
    
        # ---------------
        # - Rate Kinds
        if simulation not in avialable_lc:
            raise ValueError("not '%s' rate kind for '%s'"%(simulation,transient)+\
                             "These are: "+",".join(avialable_lc))
        # -- the output
        if transient is None:
            return eval("self.lightcurve_%s"%(simulation))
        return eval("self.lightcurve_%s_%s"%(transient,simulation))
    
    def set_model(self,model):
        """
        """
        self._side_properties["model"] = model
        
    # ========================== #
    # = LC Kind                = #
    # ========================== #
    
    # ----------------- #
    # - Ia LC         - #
    # ----------------- #
    def lightcurve_Ia_basic(self,redshifts,
                            color_mean=0,color_sigma=0.1,
                            stretch_mean=0,stretch_sigma=1,
                            x0_mean=1e-5,x0_sigma=0.1,
                            model="salt2",
                            ):
        """
        """
        # ----------------
        # - Models        
        self.set_model(sncosmo.Model(source=model))
        x0 = []
        for z in redshifts:
            self.model.set(z=z)
            mabs = normal(-19.3, 0.3)
            self.model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
            x0.append(self.model.get('x0'))
            
        ntransient = len(redshifts)
        
        return {'x0':np.array(x0),
                'x1':normal(stretch_mean, stretch_sigma, ntransient),
                'c':normal(color_mean, color_sigma, ntransient)}
    def lightcurve_Ia_realistic():
        raise NotImplementedError("To be done")
    
    def lightcurve_Ia_hostdependent():
        raise NotImplementedError("To be done")
    # ========================== #
    # = Properties             = #
    # ========================== #
    @property
    def model(self):
        return self._side_properties['model']
    
    @property
    def known_lightcurve_simulation(self):
        return self._parse_rate_("lightcurve")
    
    @property
    def known_Ia_lightcurve_simulation(self):
        return self._parse_rate_("lightcurve_Ia")
    
    
