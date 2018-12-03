#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Module containing the Basic Photometric Diagnostic largely used in Astro/Cosmo """

import warnings
import numpy   as np
from scipy import stats
from ..photometry  import get_photopoint
from ..baseobject  import Samplers, TargetHandler
from ..collection  import TargetPhotoPointCollection
from ..utils.tools import kwargs_update, load_pkl


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
    = One of these must be set. If several are given, the first is used. =
    
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
        if photopoint_i is not None and "i" not in photopoint_i.bandname:
            warnings.warn("The PhotoPoint with the bandname %s is used as a 'i' band photopoint for the MassEstimator"%photopoint_i.bandname)
        if photopoint_g is not None and "g" not in photopoint_g.bandname:
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
def gi_prior_snfhost(x):
    """ """
    return g_i_prior(x, which="snhost", image_norm=False)

def gi_prior_lange(x):
    """ """
    return g_i_prior(x, which="lange14", image_norm=False)

def gi_prior_flat(x, inrange=[-0.5,2]):
    """ """
    flagin = (x>inrange[0]) & (x<inrange[1])
    p_ = np.zeros(len(x))
    p_[flagin] = 1
    return p_

def gi_prior_normal(x, loc=0.7, scale=0.5):
    """ """
    return stats.norm.pdf(x, loc=loc, scale=scale)



def g_i_prior(x, which="snhost", image_norm=False):
    """ prior pdf distribution following The g-i distribution
    given by Rebecca Lange et al 2014 (GAMA).

    See details in the PriorGmI class.
    
    Return
    ------
    array
    """
    p = PriorGmI()
    p.set_prop_gaussians(**p.get_typical_prior(which))
    coef = 1. if not image_norm else p._Lange_plotamp
    return p.pdf(x) * coef



LANGE_ETAL_IMAGE_LOC = '/Users/mrigault/Desktop/these/code_source/Data/Literature_Images/Lange14_g_i_distrib_019_124.png'
class PriorGmI( Samplers ):
    """ Samplers object as Prior.
    From Lange et al 2014"""

    PROPERTIES = ["mub", "mur", "sigmab", "sigmar", "b_coef"]
    
    # =============== #
    #  Main Methods   #
    # =============== #
    def pdf(self, x):
        """ the probability distrubution function for the binormal prior model. """
        return self.coef_blue * stats.norm.pdf(x, **self.prop_blue ) +\
          (1-self.coef_blue) * stats.norm.pdf(x, **self.prop_red)

    def show(self, ax=None, savefile=None, show_legend=False, **kwargs):
        """ """
        from ..utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        # - Setting - #
        if ax is None:
            fig = mpl.figure(figsize=[8,5])
            ax  = fig.add_axes([0.1,0.13,0.8,0.8])
            ax.set_xlabel(r"$\mathrm{g-i\ [mag]}$",fontsize = "x-large")
            ax.set_ylabel(r"$\mathrm{frequency}$",fontsize = "x-large")
            ax.set_title(r"$\mathrm{Image\ from\ Lange\,et\,al\ 2014}$",fontsize = "x-large")
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
        
        prop = kwargs_update(dict(ls="-", color="k", lw=2, label=r"$\mathrm{Prior}(g-i)$"),**kwargs)
        # ---
        ax.imshow(self.image, extent=(0.162, 1.14, -45, 800), aspect='auto')
        x = np.linspace(0,1.2, 1000)
        
        ax.plot(x,self._Lange_plotamp * self.pdf(x), **kwargs )
        # -- 0 line
        ax.axhline(0, ls="--", color="k", alpha=0.5, zorder=6 )
        if show_legend:
            ax.legend(loc="upper right", fontsize="large")
        fig.figout(savefile=savefile)
        
    #def draw_samplers(self, nsamplers=None):
    #    """ draw sampler and set then to self.samplers based on the binormal
    #    distributions. """
    #    if nsamplers is not None:
    #        self.nsamplers = nsamplers
    #        
    #    s =  np.concatenate([stats.norm.rvs(size=self.nsamplers * self.blue_to_red_ratio,**self.prop_blue),
    #                        stats.norm.rvs(size=self.nsamplers ,**self.prop_red)])
    #    np.random.shuffle(s)
    #    self.set_samplers(s[:self.nsamplers])

    def set_prop_gaussians(self, mub, mur,
                               sigmab, sigmar,
                               b_coef):
        """ Change the properties of the two gaussians parameters.
        Parameters
        ----------
        mub, mur: [float, flat]
            Central values of the blue and red gaussians

        sigmab, sigmar: [float, float]
            Dispersion of the ble and red gaussians

        b_coef: [0<float<1] 
            Amplitude of the blue gaussian with relative to the
            red one.

        Returns
        -------
        Void
        """
        self._properties["mub"], self._properties["sigmab"] = mub, sigmab
        self._properties["mur"], self._properties["sigmar"] = mur, sigmar
        self._properties["b_coef"] = b_coef


    def get_typical_prior(self, which="snhost"):
        """ Select pre made structure:
        - snhost
        - lange14
        """
        if which == "snhost":
            return dict(mur=1.25, mub=0.85, sigmar=0.1, sigmab=0.3, b_coef=0.9)
        if which == "lange14":
            return dict(mur=0.85, mub=0.49, sigmar=0.15, sigmab=0.11, b_coef=0.7)
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
        # Lange Default: dict(loc=0.49, scale=0.11)
        if self._properties["mub"] is None:
            self._properties["mub"] = 0.85
        if self._properties["sigmab"] is None:
            self._properties["sigmab"] = 0.3
            
        return dict(loc   = self._properties["mub"],
                    scale = self._properties["sigmab"]) 
    
    @property
    def prop_red(self):
        """ """
        # Lange Default: dict(loc=0.85, scale=0.15)
        if self._properties["mur"] is None:
            self._properties["mur"] = 1.25
        if self._properties["sigmar"] is None:
            self._properties["sigmar"] = 0.1
            
        return dict(loc   = self._properties["mur"],
                    scale = self._properties["sigmar"]) 

    @property
    def coef_blue(self):
        if self._properties["b_coef"] is None:
            self._properties["b_coef"] = 0.9 # works well for local
        return self._properties["b_coef"]
        #return 1-1./(self.blue_to_red_ratio + 1)

    """
    @property
    def blue_to_red_ratio(self):
        return 2.5
    
    @property
    def _Lange_plotamp(self):
        return 50 * (self.blue_to_red_ratio + 1)
    """
    
    

# ==================================== #
#                                      #
#      Da Class for Mass Estimate      #
#                                      #
# ==================================== #
class GmISamplers( Samplers ):
    """ """
    @property
    def _default_sampling_xrange(self):
        """ """
        noprior_xrange = super(GmISamplers,self)._default_sampling_xrange
        return [np.nanmax([noprior_xrange[0], -2]), np.nanmin([noprior_xrange[1], 3.])]
    
class MassEstimate( Samplers, TargetPhotoPointCollection ):
    """
    """
    PROPERTIES = []
    SIDE_PROPERTIES = ["relation_magdisp","gi_prior"]
    DERIVED_PROPERTIES = ["gi_samplers", "gi_priored_sample", "samplers_noprior"]
    
    def __init__(self, ppoint_g=None, ppoint_i=None,
                 nsamplers=10000, relation_dispersion=0.1,
                 target=None,empty=False):
        """ """
        self.__build__()
        if empty:
            return
        if ppoint_g is not None or ppoint_i is not None:
            self.create([ppoint_i, ppoint_g], ids=["i","g"])
        
        self.taylordispersion = relation_dispersion
        self.nsamplers = nsamplers

        if target is not None:
            self.set_target(target)
            
    # =================== #
    # = Main Methods    = #
    # =================== #
    # --------- #
    #  SETTER   #
    # --------- #

    # --------- #
    #  GETTER   #
    # --------- #
    def get_estimate(self, noprior=False):
        """ Estimation of the mass based on the current samplers.
        You can get the mass without the g-i prior by setting noprior to True
        """
        if not noprior:
            return super(MassEstimate, self).get_estimate()
        
        v = np.percentile(self.samplers_noprior, [16, 50, 84])
        return v[1], v[2]-v[1], v[1]-v[0]
    
    # --------- #
    #  OTHERS   #
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

        self._set_color_samplers_()
        i_ = self.photopoints["i"].photosamplers.mag[:self.nsamplers]
        # -- Convert that to masses    
        nodisp_sampler   = taylor_mass_relation(i_,self.gi_priored_sample,
                                                self.target.distmpc)
        
        nodisp_samplernp = taylor_mass_relation(i_,self.gi_samplers.samplers,
                                                self.target.distmpc)
        
        addon_dispersion = np.random.normal(loc=0, scale=self.taylordispersion,
                                            size=self.nsamplers)
        
        self.set_samplers( nodisp_sampler + addon_dispersion)
        
        self._derived_properties["samplers_noprior"] = nodisp_samplernp + addon_dispersion
    
    # ----------- #
    #   PLOTTER   #
    # ----------- #
    def show_details(self, savefile=None, show=True, kind="Local",
                        figsize=[10,5], **kwargs):
        """ Display and advanced figure showing the details on
        how the mass estimate is drawn
        """
        import matplotlib.pyplot as mpl
        from astrobject.utils.mpladdon import figout

        # -- Axes settings 
        fig     = mpl.figure(figsize=figsize)
        axcolor = fig.add_axes([0.05,0.12,0.56,0.78])
        ax      = fig.add_axes([0.07,0.65,0.2,0.18],zorder=9)
        axcolor_prior = axcolor.twinx()
        axmass  = fig.add_axes([0.65,0.12,0.3,0.78])
        # -----------
        # - Property
        prop = kwargs_update(dict(histtype="step", bins=20, normed=True, 
                            lw="2",fill=False,
                            ec=mpl.cm.binary(0.8,1.), zorder=6), **kwargs)
        
        # - Plots    - #
        _ = self.photopoints["i"].photosamplers.show(mag=True,
                                            ax=ax, fancy_xticklabel=False,
                                            show_legend=False, show_model=False, show=False,
                                            ec=mpl.cm.copper(0.5),fc=mpl.cm.copper(0.5,0.2),
                                            bins=prop["bins"], show_estimate=False,xscale=False)
        _ = self.photopoints["g"].photosamplers.show(mag=True,
                                            ax=ax, fancy_xticklabel=False,
                                            show_legend=False,show_model=False,show=False,
                                            ec=mpl.cm.Greens(0.8),fill=False,
                                            bins=prop["bins"], show_estimate=False,xscale=False)
        # - Main g-i
        x = np.linspace(-1,3,1e3)

        # - likelihood        
        axcolor.hist(self.gi_samplers.samplers, label=r"$\mathrm{Likelihood}$",**prop)
        
        # - prior
        axcolor_prior.plot(x, self.gi_prior(x)/self.gi_prior(x).max(),
                           "k--", alpha=0.7, lw=2, scalex=False, zorder=3,
                           label=r"$\mathrm{Prior}$")
        # - Posterior
        kde = stats.gaussian_kde(self.gi_priored_sample)
        axcolor_prior.fill_between(x, kde.pdf(x)/kde.pdf(x).max(), lw=2, 
                                   edgecolor=mpl.cm.binary(0.9,0.7),
                                   facecolor=mpl.cm.Blues(0.8, 0.5),
                                   label=r"$\mathrm{Posterior}$")

        
        
        # - Derived Mass
        self.show(ax=axmass, fancy_xticklabel=False,
                show_legend=False, show_model=False,show=False,
                show_estimate=True, kde=True,
                 edgecolor=mpl.cm.binary(0.9,1),
                  facecolor=mpl.cm.binary(0.5, 0.2))
        
        v_ = self.get_estimate()
        axmass.text(0.05,0.95,
                    r"$\log(\mathrm{M_*/M_{\odot}}) = %.2f ^{+%.2f}_{-%.2f}$"%(v_[0],v_[1],v_[2]),
                    transform=axmass.transAxes, va="top",ha="left",
                    color="k", bbox=dict(facecolor='w', alpha=0.5, edgecolor="None"),
                   fontsize="x-large")
                        
        # -- infor and fancy
        
        ax.set_xlabel(r"$\mathrm{AB\ magnitude}$", fontsize="xx-large")
        axcolor.set_xlabel(r"$g-i\ \mathrm{color\ [mag]}$", fontsize="xx-large")
        axmass.set_xlabel(r"$\log(\mathrm{M_*/M_{\odot}})$", fontsize="xx-large")
        axmass.set_ylim(0, axmass.get_ylim()[-1])


        axcolor.set_xlim(-0.5,2)

        axcolor_prior.legend(loc="upper right", fontsize="x-large", frameon=False)
        axcolor_prior.set_ylim(0,1.0*1.2)
        axcolor.set_ylim(0,axcolor.get_ylim()[1]*1.2)

        _ = [ax_.set_yticks([]) for ax_ in fig.axes]
        fig.suptitle(r"$\mathrm{Derivation\ of\ the\ %s\ Stellar\ Mass:\ %s}$"%(kind,self.target.name),
                    fontsize="xx-large")
        
        fig.figout(savefile=savefile, show=show)

    # --------- #
    #  I/O      #
    # --------- #
    def load(self, data, set_samplers=True, **kwargs):
        """ load and set the current Sampler.

        Parameters
        ----------
        data: [dict or string]
            Could be:
            - a dictionary [dict] as formated by the self.data
            - a filename [string] which contains the dictionary

        set_samplers: [bool]
            Set to the current instance the sampler saved in `data` (if any)
            
        Return
        ------
        Void
        """
        if type(data) is str:
            data = load_pkl(data)
        if type(data) is not dict:
            raise ValueError("The given data is not a dictionary")

        if "g" not in data or "i" not in data:
            raise TypeError("The data does not contain 'g' and/or 'i' band data")
        # -- loading g and i bands

        if data["g"] is not None:
            g = get_photopoint(**data["g"])
        else:
            g = None
        if data["i"] is not None:
            i = get_photopoint(**data["i"])
        else:
            i = None
        self.__init__(g,i, **kwargs)
        
        if set_samplers and "samplers" in data and data["samplers"] is not None:
            #self._set_color_samplers_()
            self.set_samplers(data["samplers"])

        if data is not None and "meta" in data.keys():
            for k,v in self.data['meta']:
                self.set_meta(k,v)
            
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

    def _set_color_samplers_(self):
        """ """
        self.photopoints["i"].draw_photosamplers(nsamplers=self.nsamplers*2) # assumes less than half will be nan
        self.photopoints["g"].draw_photosamplers(nsamplers=self.nsamplers*2)
        
        neff = np.min([self.photopoints["g"].photosamplers.nsamplers,
                       self.photopoints["i"].photosamplers.nsamplers,self.nsamplers])
        if neff != self.nsamplers:
            warnings.warn("Reduced effective number of sampler (%d -> %d)for the mass because of too many nans"%(self.nsamplers, neff))
            self.nsamplers = neff

        # -- g and i samplers
        g_ = self.photopoints["g"].photosamplers._magsamples.get_random_samplers(self.nsamplers)
        i_ = self.photopoints["i"].photosamplers._magsamples.get_random_samplers(self.nsamplers)
        
        # -- Set the gi_sampler consequently
        self.gi_samplers.set_samplers(g_ - i_)
        # -- and draw the priored sampling (could be done automatically be here it is explicit)
        self._derived_properties["gi_priored_sample"] = \
          self.gi_samplers.resample(self.nsamplers, prior=self.gi_prior,
                                    xrange=self.gi_samplers._default_sampling_xrange, rand_nsample=1e3)

    def change_gi_prior(self, prior, redraw=True):
        """ Set a new prior function. 
        The first argument of this function must be the g-i color"""
        self._side_properties["gi_prior"] = prior
        if redraw:
            self.draw_samplers()
            
    # =================== #
    # = Properties      = #
    # =================== #
    @property
    def data(self):
        """ PhotoPoint Collections and samplers data """

        if "i" not in self.photopoints.keys() or "g" not in self.photopoints.keys():
            return {"estimate": [np.NaN,np.NaN,np.NaN],
                    "samplers":None,
                    "g": None if "g" not in self.photopoints.keys() else
                    self.photopoints["g"].data,
                    "i":None if "i" not in self.photopoints.keys() else
                    self.photopoints["i"].data} # to be improved
        
        return {
            "estimate" : self.get_estimate(),
            "samplers" : self.samplers,
            "g": self.photopoints["g"].data if self.photopoints["g"].meta is None
                 else kwargs_update(self.photopoints["g"].meta, **self.photopoints["g"].data),
            "i": self.photopoints["i"].data if self.photopoints["i"].meta is None
                 else kwargs_update(self.photopoints["i"].meta, **self.photopoints["i"].data)}

    @property
    def samplers_noprior(self):
        """ Samplers without prior applied on the g-i color"""
        return self._derived_properties["samplers_noprior"]

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
            self._derived_properties["gi_samplers"] = GmISamplers()
        
        return self._derived_properties["gi_samplers"]

    @property
    def gi_priored_sample(self):
        """ sample dranw fro the gi_samplers assuming the a g-i prior """
        if self._derived_properties["gi_priored_sample"] is None:
            if not self.gi_samplers.has_samplers():
                self._set_color_samplers_()
            
            self._derived_properties["gi_priored_sample"] = \
              self.gi_samplers.resample(self.nsamplers,prior=self.gi_prior,
                                        xrange=self.gi_samplers._default_sampling_xrange)
              
        return self._derived_properties["gi_priored_sample"]

    @property
    def gi_prior(self):
        """ Prior Function """
        if self._side_properties["gi_prior"] is None:
            self._side_properties["gi_prior"] = gi_prior_snfhost
        return self._side_properties["gi_prior"]
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
    
