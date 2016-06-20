#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Module containing the Basic Photometric Diagnostic largely used in Astro/Cosmo """

import warnings
import numpy   as np
from scipy import stats
from ..baseobject  import Samplers
from ..collection  import PhotoPointCollection
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
class MassEstimate( Samplers, PhotoPointCollection ):
    """
    """
    PROPERTIES = []
    SIDE_PROPERTIES = ["relation_magdisp"]
    DERIVED_PROPERTIES = []
    
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

    # --------- #
    #  GETTER   #
    # --------- #
    def draw_samplers(self, nsamplers=None):
        """ Return a monte-carlo distribution of the masses
        based on Taylor et al. 2011 relation and the given
        photopoints magnitudes
        
        You can change the number of samplers used by setting nsamplers.
        """
        if "i" not in self.photopoints.keys() or "g" not in self.photopoints.keys():
            raise AttributeError("You need 'i' and 'g' photopoints to be able to draw the samplers")
        if nsamplers is not None:
            self.nsamplers = nsamplers
        nodisp_sampler = \
          taylor_mass_relation(self.photopoints["i"].magdist.rvs(self.nsamplers),
                               self.photopoints["g"].magdist.rvs(self.nsamplers),
                               self.target.distmpc)
        addon_dispersion = np.random.normal(loc=0, scale=self.taylordispersion,
                                            size=self.nsamplers)
        
        self.set_samplers( nodisp_sampler + addon_dispersion)
    
    # =================== #
    #  Internal           #
    # =================== #
    # this is currently default in astrobject, but it might change
    # Let's set it here
    def _set_rvdist_(self):
        """ set the rvdistribution.
        This method defines which kind of rv_continuous distribution you use
        """
        return stats.loggamma(*stats.loggamma.fit(self.samplers,
                                        loc=np.median(self.samplers)))
    @property
    def rvdist_info(self):
        """ information about the rvdistribution """
        return r"$\mathrm{loggamma(%.1e, %.1e, %.1e)}$"%(self.rvdist.args)
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

    # --------------
    # - Taylor Stuff
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
    
