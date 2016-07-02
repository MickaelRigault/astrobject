#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Module containing the Basic Photometric Diagnostic largely used in Astro/Cosmo """

import warnings
import numpy   as np
from scipy import stats
from ..baseobject  import Samplers, TargetHandler
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
def taylor_mass_relation(mag_i, gi_color, distmpc):
    """ Estimate the stellar-mass of a galaxy based on its i and g bands.
    
    This fuction uses the eq. 8 of Taylor et al 2011
    http://adsabs.harvard.edu/abs/2011MNRAS.418.1587T

    log M∗/[M⊙] = 1.15 + 0.70(g − i) − 0.4Mi
    
    Error floor of 0.1dex assumed to account for the spread of the
    (g-i) vs. log M∗/[M⊙] relation. See Taylor et al. 2011 (5.3).
    This error is added in quadrature with the errors of g and i.

    Parameters
    ----------
    
    mag_i, gi_color: [(array of) float,(array of) float]
        Restframe magnitude in "i" band and "g-i" color respectively (calibrated for sdss)

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

    return 1.15 + 0.70*(gi_color) - 0.4*mag_i + (0.4*5*(np.log10(distmpc*1.e6) - 1))

######################
#                    #
#   G-I PRIOR        #
#                    #
######################
def g_i_prior(x, image_norm=False):
    """ prior pdf distribution following The g-i distribution
    given by Rebecca Lange et al 2014 (GAMA).

    See details in the PriorGmI class.
    
    Return
    ------
    array
    """
    p = PriorGmI()
    coef = 1. if not image_norm else p._Lange_plotamp
    return p.pdf(x) * coef

LANGE_ETAL_IMAGE_LOC = '/Users/mrigault/Desktop/these/code_source/Data/Literature_Images/Lange14_g_i_distrib_019_124.png'
class PriorGmI( Samplers ):
    """ Samplers object as Prior.
    From Lange et al 2014"""

    # =============== #
    #  Main Methods   #
    # =============== #
    def pdf(self, x):
        """ the probability distrubution function for the binormal prior model. """
        return self.coef_blue * stats.norm.pdf(x, **self.prop_blue ) +\
          (1-self.coef_blue) * stats.norm.pdf(x, **self.prop_red)

    def show(self, ax=None, savefile=None, **kwargs):
        """ """
        from ..utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        # - Setting - #
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            ax.set_xlabel(r"$\mathrm{g-i\ [mag]}$",fontsize = "x-large")
            ax.set_ylabel(r"$\mathrm{frequency}$",fontsize = "x-large")
            ax.set_title(r"$\mathrm{Image\ from\ Lange\,et\,al\ 2014}$",fontsize = "x-large")
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
        
        prop = kwargs_update(dict(ls="-", color="k", lw=2, label=r"$\mathrm{prior\ reconstruction}$"),**kwargs)
        # ---
        ax.imshow(self.image, extent=(0.162, 1.14, -45, 800), aspect='auto')
        x = np.linspace(0,1.2, 1000)
        
        ax.plot(x,self._Lange_plotamp * self.pdf(x), **kwargs )
        fig.figout(savefile=savefile)
        
    def draw_samplers(self, nsamplers=None):
        """ draw sampler and set then to self.samplers based on the binormal
        distributions. """
        if nsamplers is not None:
            self.nsamplers = nsamplers
            
        s =  np.concatenate([stats.norm.rvs(size=self.nsamplers * self.blue_to_red_ratio,**self.prop_blue),
                            stats.norm.rvs(size=self.nsamplers ,**self.prop_red)])
        np.random.shuffle(s)
        self.set_samplers(s[:self.nsamplers])



    # ================= #
    #   Properties      #
    # ================= #
    # ================= #
    #   Properties      #
    # ================= #
    @property
    def image(self):
        """ """
        import matplotlib.image as mpimg
        return mpimg.imread(LANGE_ETAL_IMAGE_LOC)
    
    @property
    def prop_blue(self):
        """ """
        return dict(loc=0.49, scale=0.11)
    
    @property
    def prop_red(self):
        """ """
        return dict(loc=0.85, scale=0.15)

    @property
    def coef_blue(self):
        return 1-1./(self.blue_to_red_ratio + 1)

    
    @property
    def blue_to_red_ratio(self):
        """ """
        return 2.5
    @property
    def _Lange_plotamp(self):
        return 50 * (self.blue_to_red_ratio + 1)
    
    
    

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
    DERIVED_PROPERTIES = ["gi_samplers"]
    
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

        self.photopoints["i"].magsamplers.draw_samplers(nsamplers=self.nsamplers*2) # assumes less than half will be nan
        self.photopoints["g"].magsamplers.draw_samplers(nsamplers=self.nsamplers*2)
        neff = np.min([self.photopoints["g"].magsamplers.nsamplers,self.photopoints["i"].magsamplers.nsamplers, self.nsamplers])
        if neff != self.nsamplers:
            warnings.warn("Reduced effective number of sampler (%d -> %d)for the mass because of too many nans"%(self.nsamplers, neff))
            self.nsamplers = neff

        # -- g and i samplers
        g_ = self.photopoints["g"].magsamplers.samplers[:self.nsamplers]
        i_ = self.photopoints["i"].magsamplers.samplers[:self.nsamplers]
        
        # -- Set the gi_sampler consequently
        self.gi_samplers.set_samplers(g_ - i_)
        
        # -- Convert that to masses    
        nodisp_sampler   = taylor_mass_relation(i_,
                                                self.gi_samplers.resample(self.nsamplers,prior=g_i_prior,xrange=[-0.5,2]),
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
        return r"$\mathrm{loggamma(%.1e,\, %.1e,\, %.1e)}$"%(self.rvdist.args)
    # =================== #
    # = Properties      = #
    # =================== #
    @property
    def data(self):
        """ PhotoPoint Collections and samplers data """
        return {
            "estimate" : self.get_estimate(),
            "samplers" : self.samplers,
            "g": self.photopoints["g"].data,
            "i": self.photopoints["i"].data}
    
    @property
    def target(self):
        
        if self._side_properties["target"] is None and self.has_data():
            if self.photopoints["i"].has_target():
                self.set_target(self.photopoints["i"].target)
            elif self.photopoints["g"].has_target():
                self.set_target(self.photopoints["g"].target)
                
        return self._side_properties["target"]


    @property
    def gi_samplers(self):
        """ Samplers of the g-i color """
        if self._derived_properties["gi_samplers"] is None:
            self._derived_properties["gi_samplers"] = Samplers()
        
        return self._derived_properties["gi_samplers"]
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
    
