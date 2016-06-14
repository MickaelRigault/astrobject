#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Module containing the Basic Photometric Diagnostic largely used in Astro/Cosmo """

import warnings
import numpy   as np
from scipy import stats

from ..collection   import PhotoPointCollection
from ..utils.tools import kwargs_update


__all__ = ["get_massestimator"]


def get_massestimator(photopoints=None,samplers=None,dist_param=None,
                      nsamples=10000,
                      **kwargs):
    """ Get the photopoint collection instance that enables
    to derive the stellar mass of a galaxy.

    The natural non-parametrized mass estimation is based on
    rest-frame 'i' and 'g' band measurements following
    eq. 8 of Taylor et al 2011 (2011MNRAS.418.1587T)
    ```
      log M∗/[M⊙] = 1.15 + 0.70(g − i) − 0.4Mi
    ```
    
    You can otherwize provide 'samplers' of the mass estimate
    or a parametrized version of the 'samplers' using 'dist_param'.
    It has been tested, the natural method returns a galaxy mass
    pdf well described by a loggamma distribution (3-parameters).
    By providing a 'dist_param', which parametrizes the scipy's loggamma
    object, this will draw a fake sampler to set the MassEstimate.
    
    Parameters
    ----------
    = Setter: one must be used. If several are given, the first is used. =
    
    photopoints: [2 PhotoPoints: photopoint_g,photopoint_i]
        Astrobject's PhotoPoint of the g and i band measurements of the galaxy
        Magnitude should be rest-frame magnitude in the sdss system.
        
    samplers: [N-float array]
        Array of potential mass estimate. From this sampling of masses is drawn
        the pdf of the galaxy's stellar mass estimate.
        
    dist_param: [3-floats]
        Parameters defining scipy's loggamma functions. The mass estimate
        will be base on the random drawing of samplers from this distribution.

    = other options =
    
    nsamples: [int] -optional-
        Number of sampler used to estimae the galaxy's stellar mass.
        If you set the MassEstimator with the 'samplers' input, nsamples
        is ignored.


    **kwargs goes to MassEstimate's __init__ method
    
    Return
    ------
    MassEstimate (PhotoPointCollection's child)
    """
    # ------------- #
    # - Setting   - #
    # ------------- #
    if photopoints is not None:
        photopoint_g, photopoint_i = photopoints
        if "i" not in photopoint_i.bandname:
            warnings.warn("The PhotoPoint with the bandname %s is used as a 'i' band photopoint for the MassEstimator"%photopoint_i.bandname)
        if "g" not in photopoint_g.bandname:
            warnings.warn("The PhotoPoint with the bandname %s is used as a 'g' band photopoint for the MassEstimator"%photopoint_g.bandname)
            
        return MassEstimate(photopoint_g, photopoint_i,**kwargs)

    if samplers is not None:
        m = MassEstimate(**kwargs)
        m.set_samplers(samplers)
        return m
    
    if dist_param is not None:
        samplers = stats.loggamma(*dist_param).rvs(nsamples)
        return get_massestimator(samplers=samplers, **kwargs)
    
    if "empty" in kwargs.keys():
        return MassEstimate(**kwargs)
    
    raise ValueError("You need to provide a setter (photopoints, samplers or dist_param) or to set empty=True")


##########################################
#                                        #
# Mass Measurements                      #
#                                        #
##########################################
def taylor_mass_relation(mag_i,mag_g, distmpc):
    """ Estimate the stellar-mass of a galaxy based on its i and g bands.
    
    This fuction uses the eq. 8 of Taylor et al 2011
    http://adsabs.harvard.edu/abs/2011MNRAS.418.1587T

    log M∗/[M⊙] = 1.15 + 0.70(g − i) − 0.4Mi
    
    Error floor of 0.1dex assumed to account for the spread of the
    (g-i) vs. log M∗/[M⊙] relation. See Taylor et al. 2011 (5.3).
    This error is added in quadrature with the errors of g and i.

    Parameters
    ----------
    
    mag_i, mag_g: [(array of) float,(array of) float]
        Restframe magnitude in "i" and "g" band respectively (calibrated for sdss)

    distmpc: [float]
        Distance in Mega parsec of the galaxy

    Return
    ------
    float (the mass of the target galaxy)
    """
    # Rewrite of the Taylor Equation
    # Mi = i - 5*(np.log10(ppoint_i.target.distmpc*1.e6) - 1)
    #      1.15 + 0.70*(g - i) - 0.4*i + (0.4*5*(np.log10(distmpc*1.e6) - 1))
    #      1.15 + 0.70*g - 1.1*i + (0.4*5*(np.log10(distmpc*1.e6) - 1))
    return 1.15 + 0.70*mag_g - 1.1*mag_i + (0.4*5*(np.log10(distmpc*1.e6) - 1))
        
# ==================================== #
#                                      #
#      Da Class for Mass Estimate      #
#                                      #
# ==================================== #
class MassEstimate( PhotoPointCollection ):
    """
    """
    PROPERTIES = ["nsamplers"]
    SIDE_PROPERTIES = ["relation_magdisp"]
    DERIVED_PROPERTIES = ["loggamma_param","samplers"]
    
    def __init__(self, ppoint_g=None, ppoint_i=None,
                 nsamplers=10000, relation_dispersion=0.1, empty=False):
        """ """
        self.__build__()
        if empty:
            return
        if ppoint_g is not None or ppoint_i is not None:
            self.create([ppoint_i, ppoint_g], ids=["i","g"])
        
        self.taylordispersion = relation_dispersion
        self.nsamplers = nsamplers
    # =================== #
    # = Main Methods    = #
    # =================== #
    # --------- #
    #  SETTER   #
    # --------- #
    def set_samplers(self, samplers):
        """ Manually set the mass-samplers from which the
        galaxy's stellar mass will be estimated.
        """
        self._derived_properties['samplers'] = samplers
        self._derived_properties['loggamma_param'] = None
        
    # --------- #
    #  GETTER   #
    # --------- #
    def draw_samplers(self):
        """ Return a monte-carlo distribution of the masses
        based on Taylor et al. 2011 relation and the given
        photopoints magnitudes
        """
        if "i" not in self.photopoints.keys() or "g" not in self.photopoints.keys():
            raise AttributeError("You need 'i' and 'g' photopoints to be able to draw the samplers")

        nodisp_sampler = \
          taylor_mass_relation(self.photopoints["i"].magdist.rvs(self.nsamplers),
                               self.photopoints["g"].magdist.rvs(self.nsamplers),
                               self.target.distmpc)
        addon_dispersion = np.random.normal(loc=0, scale=self.taylordispersion,
                                            size=self.nsamplers)
        
        self.set_samplers( nodisp_sampler + addon_dispersion)
    
    def get_estimate(self):
        """
        """
        if not self.has_samplers():
            self.draw_samplers()
            
        v = np.percentile(self.samplers, [16, 50, 84])
        return v[1], v[2]-v[1], v[1]-v[0]
    
    # --------- #
    #  PLOT     #
    # --------- #
    def show(self, savefile=None, show=True, ax=None,
             propmodel={},**kwargs):
        """
        """
        from astrobject.utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        # -------------
        # - Input
        if ax is None:
            fig = mpl.figure(figsize=[8,6])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            ax.set_xlabel(r"$\mathrm{\log(M/M\odot)}$",fontsize = "x-large")
            ax.set_ylabel(r"$\mathrm{frequency}$",fontsize = "x-large")
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
        # -------------
        # - Prop
        prop = kwargs_update(dict(histtype="step", bins=50, normed=True,
                    lw="2",fill=True, fc=mpl.cm.Blues(0.5,0.4),
                    ec=mpl.cm.Blues(1.,1.)),
                    **kwargs)
        propmodel_ = kwargs_update(dict(scalex=False, color="g",lw=2),
                    **propmodel)
        # -------------
        # - Samplers
        if not self.has_samplers():
            warnings.warn("Samplers created for the plot method")
            self.draw_samplers()
            
        h = ax.hist(self.samplers,**prop)

        med,highmed,lowmed = self.get_estimate()
        # - show estimate
        x = np.linspace(med-lowmed*10,med+highmed*10,10000)
        pl = ax.plot(x,self.massdist.pdf(x), **propmodel_)
        ax.axvline(med, color="0.5", zorder=2)
        ax.axvline(med-lowmed, color="0.6", ls="--", zorder=2)
        ax.axvline(med+highmed, color="0.6", ls="--", zorder=2)
        
        self._plot["figure"] = fig
        self._plot["ax"] = ax
        self._plot["plot"] = [h,pl]
        
        fig.figout(savefile=savefile,show=show)
        
        return self._plot
    
    # =================== #
    # = Properties      = #
    # =================== #
    @property
    def target(self):
        if self.has_data():
            return self.photopoints["i"].target if self.photopoints["i"].has_target() \
            else self.photopoints["g"].target
        return None
    
    def has_target(self):
        """ test if the instance has a target (via the photopoints) """
        return self.target is not None

    @property
    def massdist(self):
        """
        distribution (loggamma) of the mass.
        
        Return
        ------
        scipy's stats.loggamma initialized

        Example
        -------
        get the pdf of the mass distribution: self.massdist.pdf(x)
        """
        return stats.loggamma(*self._mass_loggamma_param)
    
    @property
    def _mass_loggamma_param(self):
        """ """
        if self._derived_properties['loggamma_param'] is None:
            self._derived_properties['loggamma_param'] = \
              stats.loggamma.fit(self.samplers, loc=np.median(self.samplers))
        return self._derived_properties['loggamma_param']
                
    # --------------
    # - Taylor Stuff
    @property
    def samplers(self):
        return self._derived_properties['samplers']
    
    def has_samplers(self):
        return self.samplers is not None
    
    @property
    def nsamplers(self, default=10000):
        if self._properties["nsamplers"] is None:
            self._properties["nsamplers"] = default
        return self._properties["nsamplers"]
    
    @nsamplers.setter
    def nsamplers(self, value):
        """ Set the number of random drawing used to build the mass pdf """
        if value<10:
            raise ValueError("set a value greater than 10!")
        self._properties["nsamplers"] = value

    @property
    def taylordispersion(self, default=0.1):
        """ Natural dispersion of the Taylor et al. 2011
        relation. They quote 0.1mag
        http://adsabs.harvard.edu/abs/2011MNRAS.418.1587T """
        if self._side_properties["relation_magdisp"] is None:
            self._side_properties["relation_magdisp"] = default
        return self._side_properties["relation_magdisp"]

    @taylordispersion.setter
    def taylordispersion(self, value):
        """ Set the intrinsic dispersion of the Taylor Relation """
        if value < 0:
            raise ValueError("The intrinsic dispersion must be positive")
        self._side_properties["relation_magdisp"] = value
    
