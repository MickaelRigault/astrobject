#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module defines the basic low-level objects"""

import copy
import warnings

import numpy   as np

# - PropObject
from propobject import BaseObject

# - Local Tools
from .utils.tools import load_pkl, dump_pkl, kwargs_update

__all__     = ["get_target"]


#######################################
#                                     #
#  Function to get the Objects        #
#                                     #
#######################################

def get_target(name=None,zcmb=None, ra=None, dec=None, mjd=None,
               type_=None, mwebmv=None,zcmb_err=None, 
               **kwargs):
    """ Create an AstroTarget

    Parameters:
    -----------
    (all are optional, you may create an empty Target)
    
    name: [string]
        name of the astro-object

    zcmb,zcmb_err: [float,float]
        redshift and its error in the cmb frame. This will be used
        to derive distances etc. given the used cosmology

    ra, dec: [float, float]
        right-ascention / declination of the object. (degree favored).
                                   
    mjd: [float]:
        Modified Julian Date associated to the target
        
    type_:[string]
        type of the astro-object (galaxy, sn, Ia, Ibc...)
        (no predifined list type so far, but this could append)

    mwebmv: [float]
        Force the Milky way extinction for this object.
        Otherwise this extinction depend on the object coordinate.
        ->Use this use Caution<-

    **kwargs  goes to AstroTarget's init could e.g. be:
        cosmo:[astropy.cosmology]
            Cosmology used to derive the distances and so on. Default is Planck15

        empty: [bool]
            To force having an empty Target no matter what you gave before
                              
    Return
    ------
    AstroTarget
    """
    if name is None and "object" in kwargs.keys():
        name = kwargs.pop("object")
    forced_mwebmv = mwebmv if mwebmv is not None else kwargs.pop("forced_mwebmv",kwargs.pop("MWebmv", kwargs.pop("MWebv",None)))
    
    dec = kwargs.pop("Dec",dec)
    ra = kwargs.pop("Ra",ra)
    empty = kwargs.pop("empty",False)
    
    return AstroTarget(name=name,zcmb=zcmb,
                       ra=ra,dec=dec,mjd=mjd,
                       type_=type_,
                       forced_mwebmv=forced_mwebmv,
                       empty=empty,**kwargs).copy()    


#######################################
#                                     #
#    Low Level Object Classes         #
#                                     #
#######################################

class Samplers( BaseObject ):
    """ Class to Handle the Samplers for scipy.stats rv_continuous distribution """
    PROPERTIES = ["nsamplers"]
    DERIVED_PROPERTIES = ["samplers","rvdist"]

    def __init__(self, sampler=None, empty=False):
        """ """
        self.__build__()
        if empty:
            return
        if sampler is not None:
            self.set_samplers(sampler)

    def __call__(self):
        return self.samplers
    
    # =================== #
    #   Main Methods      #
    # =================== #
    def pdf(self, x, **kwargs):
        """ Probability distribution function based on the estimated rvdist (gaussian_kde by default) """
        from scipy.stats import rv_continuous
        return self.rvdist.pdf(x, **kwargs) if rv_continuous in self.rvdist.__class__.__mro__ \
          else self.rvdist.evaluate(x, **kwargs)
    
    def resample(self, size, prior=None,
                 xrange=None, rand_nsample=1e3, **kwargs):
        """ draw and return a new set on samplers given the estimated rvdist.

        Prior: You can add a prior distribution. The resampling will then
               combine the Samplers.pdf and the prior function to randomly
               draw samplers

        Parameters:
        -----------

        size: [int]
            Size of the sample you want to create from the rvdist (uses rvdist.rvs)
            (of resample for gaussian_kde)

        // Prior functionality -optional-
        
        prior: [function] -optional-
            Probability Distribution Function (pdf) for the prior. Its first argument
            should be 'x', the absysce where the pdf is estimated [pdf(x)].

        xrange: [2-array (min, max)] -optional-
            If you set a prior, the random sampling will be estimated within this interval
            If nothing provided, the sampling will be that of the current sampler

        rand_nsample: [int] -optional-
            number of value used for the random sampling. This should typically
            an order of magnitude larger than the typical plot bining you are using
            (1000 by default)
        
        **kwargs are any option to the prior function: this uses prior(x, **kwargs)
        """
        if prior is None:
            from scipy.stats import rv_continuous
            return self.rvdist.rvs(size) if rv_continuous in self.rvdist.__class__.__mro__ \
              else self.rvdist.rvs(size, xrange=xrange, nsample=rand_nsample)

        # -----------
        # - Inputs
        if xrange is None or len(xrange) != 2:
            # -- it is not a kde
            xrange = self._default_sampling_xrange
            
        # - The random sampler seed
        x = np.linspace(xrange[0],xrange[1],rand_nsample)
        # the new pdf
        
        pdf = self.pdf(x) * prior(x, **kwargs)
        return np.random.choice(x, p= pdf / pdf.sum(), size=size)

    @property
    def _default_sampling_xrange(self):
        """ usefule tools for plotting """
        from scipy.stats import rv_continuous
        dataset = self.rvdist.rvs(1000) if rv_continuous in self.rvdist.__class__.__mro__ \
              else self.rvdist.dataset
        scale = np.nanmax(dataset) - np.nanmin(dataset)
        return [np.nanmin(dataset) - scale*0.05, np.nanmax(dataset) + scale*0.05]
    
    # --------- #
    #  SETTER   #
    # --------- #
    def set_samplers(self, samplers):
        """ Set the samplers to use """
        self._derived_properties["samplers"] = samplers
        self.nsamplers = len(self.samplers)
        
    # --------- #
    #  GETTER   #
    # --------- #
    def draw_samplers(self):
        """ To Be Implemented in the child class"""
        raise NotImplementedError(" The draw_samplers() method has not been implemented ")

    def get_samplers_cdf(self, value):
        """returns the fraction of samplers below or equal to the given value """
        return len(self.samplers[self.samplers<=value]) / np.float(len(self.samplers))

    def get_random_samplers(self, n):
        """ return n ramdomly selected samplers """
        if not self.has_samplers():
            self.draw_samplers()
        
        def get_shuffledcopy(x):
            x_ = x.copy()
            np.random.shuffle(x_)
            return x_
        
        return get_shuffledcopy(self.samplers)[:n]
        
    def get_estimate(self):
        """ Estimation of the True parameters based on the current samplers.
        Return
        ------
        value [+sigma, -sigma] (such that value = 50% , value+sigma=84%, value-sigma=16%)
        """
        if not self.has_samplers():
            self.draw_samplers()
            
        v = np.percentile(self.samplers, [16, 50, 84])
        return v[1], v[2]-v[1], v[1]-v[0]
    
    # --------- #
    #  PLOT     #
    # --------- #
    def show(self, savefile=None, show=True, ax=None,
             show_model=True, propmodel={}, xlabel="",
             fancy_xticklabel=False, kde=False,
             show_legend=True, logscale=False,
             show_estimate= True,xscale=True,yscale=True,
             **kwargs):
        """ Show the samplers and the derived rv_distribution (scipy.stats)

        Parameters
        ----------

        kde: [bool] -optional-
            Show a kde shape instead of an histogram.
            
        Return
        ------
        dict (plot information)
        """
        
        from .utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        # -------------
        # - Input
        if ax is None:
            fig = mpl.figure(figsize=[8,6])
            ax  = fig.add_axes([0.16,0.16,0.73,0.73])
            ax.set_xlabel(xlabel ,fontsize = "x-large")
            ax.set_ylabel(r"$\mathrm{frequency}$",fontsize = "x-large")
        elif "imshow" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
        # -------------
        # - Prop
        if not kde:
            prop = kwargs_update(dict(histtype="step", bins=40, normed=True,
                        lw="2",fill=True, fc=mpl.cm.Blues(0.5,0.4),
                        ec=mpl.cm.Blues(1.,1.)),
                        **kwargs)
        else:
            prop = kwargs_update(dict(lw="2",edgecolor=mpl.cm.binary(0.9,1),
                                      facecolor=mpl.cm.binary(0.5, 0.2)),
                                **kwargs)
            
        propmodel_ = kwargs_update(dict(scalex=False, color="g",lw=2,
                                        label=self.rvdist_info),
                    **propmodel)
        # -------------
        # - Samplers
        med,highmed,lowmed = self.get_estimate()
        if not logscale:
            xrange = self.samplers.min()-lowmed,self.samplers.max()+highmed
        else:
            xrange = np.log10(med)-np.abs((np.log10(med)-np.log10(med-lowmed)))*7,\
              np.log10(med)+np.abs((np.log10(med)-np.log10(med+highmed)))*7
            
        x = np.linspace(xrange[0],xrange[1],1e3)
        
        if not self.has_samplers():
            warnings.warn("Samplers created for the plot method")
            self.draw_samplers()

        if not kde:
            h = ax.hist(self.samplers if not logscale else np.log10(self.samplers), **prop)
        else:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(self.samplers if not logscale else np.log10(self.samplers))
            
            normed = 1. if "normed" in prop.keys() and not prop["normed"] else kde.pdf(x).max()
            for name in prop.keys():
                if name in ["bins", "prop","histtype"]:
                    prop.pop(name)
                    
            h = ax.fill_between(x, kde.pdf(x)/normed, **prop)
        
            
        # - show estimate
        pl = ax.plot(x,self.rvdist.pdf(x), **propmodel_) if show_model and self.rvdist is not None and not logscale else None
        if show_estimate:
            if not logscale:
                ax.axvline(med, color="0.5", zorder=2)
                ax.axvline(med-lowmed, color="0.6", ls="--", zorder=2)
                ax.axvline(med+highmed, color="0.6", ls="--", zorder=2)
            else:
                ax.axvline(np.log10(med), color="0.5", zorder=2)
                ax.axvline(np.log10(med-lowmed), color="0.6", ls="--", zorder=2)
                ax.axvline(np.log10(med+highmed), color="0.6", ls="--", zorder=2)
        
        # - Legend
        if show_legend:
            legend_ = ax.legend(loc="upper left", frameon=False, fontsize="large")

        # -- out
        if yscale:
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.2) # for the legend
        if xscale:
            ax.set_xlim(xrange[0],xrange[1])
        # - Fancy
        
        if fancy_xticklabel:
            ax.set_xticklabels(["{:.2e}".format(float(a_)) if a_ is not None and a_ !="" else ""
                                for a_ in ax.get_xticks()],
                                rotation=30, ha="right")
            
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["plot"]   = [h,pl]
        self._plot["prop"]   = prop
        if show_legend:
            self._plot["legend"] = legend_
        fig.figout(savefile=savefile,show=show)
        
        return self._plot

    # --------- #
    #  I/O      #
    # --------- #
    def writeto(self, fileout):
        """ dump the data property inside the given file """
        
        dump_pkl(self.data, fileout)

        
    # =================== #
    #   Properties        #
    # =================== #
    @property
    def data(self):
        """ dictionary containing the basic information about the instance """
        return {"estimate":self.get_estimate(),
                "samplers":self.samplers}

    @property
    def rvdist(self):
        """
        distribution of the samplers
        
        Return
        ------
        scipy's stats.rv_continuous initialized (frozen)

        Example
        -------
        get the pdf of the distribution: self.rvdist.pdf(x)
        """
        if self._derived_properties['rvdist'] is None or \
          np.any(self._samplers10_rvdist != self.samplers[:10]):
            # Test if the samplers did not changed, based on the first 10
            self._derived_properties['rvdist'] = self._set_rvdist_()
            self._samplers10_rvdist = self.samplers[:10].copy()
            
        return self._derived_properties['rvdist']

    def _set_rvdist_(self):
        """ set the rvdistribution.
        This method defines which kind of rv_continuous distribution you use
        """
        from astrobject.utils.decorators import make_method
        from scipy import stats
        kde = stats.gaussian_kde(self.samplers)
        
        @make_method(stats.kde.gaussian_kde)
        def cdf(kde, value):
            return kde.integrate_box(-np.inf, value)
        
        @make_method(stats.kde.gaussian_kde)
        def rvs(kde, size, xrange=None, nsample=1e3):
            """ random distribution following the estimated pdf.
            random values drawn from a finite sampling of `nsample` points
            """
            if xrange is None:
                scale = np.nanmax(kde.dataset) - np.nanmin(kde.dataset)
                xrange = np.nanmin(kde.dataset) - scale*0.1, np.nanmax(kde.dataset) + scale*0.1
            
            x = np.linspace(xrange[0], xrange[1], nsample)
            return np.random.choice(x, p= kde.pdf(x) / kde.pdf(x).sum(), size=size)
            
        return kde
    
    @property
    def rvdist_info(self):
        """ information about the rvdistribution """
        return ""
    
    # -------------- #
    #   Samplers   - #
    # -------------- #    
    @property
    def samplers(self):
        return self._derived_properties["samplers"]
    
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
        if value is None:
            raise ValueError("value given for nsamplers cannot be 'None' ")
        if value<10:
            raise ValueError("set a value greater than 10! %s"%value)
        self._properties["nsamplers"] = value

###########################
#                         #
#   AstroTarget           #
#                         #
###########################
class AstroTarget( BaseObject ):
    """
    This is the default astrophysical object that has basic
    information attached to it, like zcmb, ra, dec, mwebmv.

    You will also have access to derived parameter like the distances
    (meter and mpc) or the angular/physical size ratio based on the
    given cosmology and redshift (zcmb).
    """
    # -- This is set to ease inheritance and tests -- #
    # If you change that, function that needs astrotarget
    # could crash (i.e. astroimages.Aperture)
     
    __nature__ = "AstroTarget"

    PROPERTIES         = ["zcmb","ra","dec","name","zcmb.err"]
    SIDE_PROPERTIES    = ["cosmology","literature_name",
                          "type","mwebmv","sfd98_dir","mjd"]
    DERIVED_PROPERTIES = ["distmeter","distmpc","arc_per_kpc"]

    # -------------------- #
    # Internal Properties  #
    # -------------------- #
    # =========================== #
    # = Initialization          = #
    # =========================== #
    def __init__(self,name=None,zcmb=None,zcmb_err=None,
                 ra=None,dec=None, mjd=None,
                 type_=None,cosmo=None,
                 load_from=None,empty=False,sfd98_dir=None,
                 **kwargs):
        """
         = Initialize the AstroTarget Function =
        
        Parameters
        ----------
        
        name: [string]             name of the astro-object

        zcmb,zcmb_err: [float]      redshift/error in the cmb frame. This will be used
                                   to derive distmeter etc. depending on the cosmology

        ra: [float]                right-ascention of the object. (degree favored).
                                   *ra* and *dec* must have the same unit.
                                   
        dec: [float]               declination of the object. (degree favored)
                                   *ra* and *dec* must have the same unit.
                                   
        mjd: [float]               modified Julian Date associated to the target
        
        type_:[string]             type of the astro-object (galaxy, sn, Ia, Ibc...)
                                   (no predifined list type so far. It could append)

        cosmo:[astropy.cosmology]  the cosmology used to derive the distances etc.

        load_from: [dict]          a dictionary having an entry for each of the
                                   fundamental
                                   parameter (name, ra, dec, zcmb...)

        empty: [bool]              Does not do anything, just loads an empty object.
                                   (Careful with that)
        
        sfd98_dir: [string]        directory of SFD98 dust map files if not set
                                   in .astropy/config/sncosmo.cfg
                                   files SFD_[dust,mask]_4096_[n,s]gp.fits are
                                   available at http://sncosmo.github.io/data/dust/[file]

        Returns
        -------
        Void
        """
        self.__build__()
        
        if empty:
            return
        
        if load_from is not None:
            self.load(load_from,**kwargs)
            return

        self.define(name,zcmb,ra,dec,mjd=mjd,
                    cosmo=cosmo,type_=type_,sfd98_dir=sfd98_dir,
                    **kwargs)
        return
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def define(self,name,zcmb,ra,dec,mjd=None,
               cosmo=None,type_=None,sfd98_dir=None,
               forced_mwebmv=None,zcmb_err=None):
        """
        This function enables you to define the fundamental object parameter
        upon which the entire object will be defined.
        """
        self.name = name
        self.set_zcmb(zcmb,zcmb_err)
        self.radec = ra,dec

        self._side_properties["mjd"] = mjd
        self.type  = type_
        if cosmo is None:
            from astropy.cosmology import Planck15
            cosmo = Planck15
            print("Planck 2015 cosmology used by default")
            warnings.warn("Planck 2015 cosmology used by default")
            
        self.set_cosmo(cosmo)
        self._update_()
        if forced_mwebmv is not None:
            self.set_mwebmv(forced_mwebmv,force_it=True)
        
    def writeto(self,output_file,**kwargs):
        """
        This function save the basic information of the object
        in the given output_file. options goes to dump_pkl
        """
        dump_pkl(self._properties,output_file,**kwargs)
        return

    def load(self,properties,side_properties=None):
        """
        """
        if type(properties) == str:
            candidatepkl = load_pkl(properties)
        elif type(properties) == dict:
            candidatepkl = properties
            
        for k in self._fundamental_parameters:
            if k not in candidatepkl.keys():
                raise AttributeError("'%s' is a requested fundamental "+\
                                     "parameter and is not in the input dictionnary"%k)

        cosmo = candidatepkl.pop("cosmo")
        self.define(candidatepkl["name"],candidatepkl["zcmb"],
                    candidatepkl["ra"],candidatepkl["dec"],
                    mjd=candidatepkl.pop("mjd",None),
                    cosmo = cosmo,
                    type_ = candidatepkl.pop("type",None),
                    forced_mwebmv = candidatepkl.pop("mwebmv",None))
        
        
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def data(self):
        """ dictionary containing the target information """
        return  dict(name=self.name, ra=self.ra, dec=self.dec, mjd=self.mjd,
                     type_=self.type, cosmo= (self.cosmo.name if self.cosmo is not None else None),
                     zcmb = self.zcmb, zcmb_err=self.zcmb_err)
    
    # ------------------ #
    # - Object Name    - #
    # ------------------ #
    @property
    def name(self):
        return self._properties["name"]
    
    @property
    def _literature_name(self):
        return self._side_properties["literature_name"]
    
    @name.setter
    def name(self,value):
        self._properties["name"] = value
        self._check_literature_name_()
                          
    # ------------------ #
    # - Redshift       - #
    # ------------------ #
    @property
    def zcmb(self):
        return self._properties["zcmb"]
    
    def set_zcmb(self,value, error=None):
        self._properties["zcmb"] = value
        self._properties["zcmb.err"] = error
        self._update_distance_()

    @property
    def zcmb_err(self):
        return self._properties["zcmb.err"]
    
    # ------------------ #
    # - Coordinate     - #
    # ------------------ #
    @property
    def ra(self):
        return self._properties["ra"]
    @property
    def dec(self):
        return self._properties["dec"]
    @property
    def radec(self):
        return self.ra,self.dec
    
    @radec.setter
    def radec(self,value):
        if np.shape(value) != (2,):
            raise SyntaxError("radec must have two input parameters: ra,dec")
        
        self._properties["ra"],self._properties["dec"] = value
        self._side_properties["mwebmv"] = None

    @property
    def mjd(self):
        """ Modified Julian date associated to the target """
        return self._side_properties["mjd"]
    # ------------------ #
    # - MW Extinction  - #
    # ------------------ #
    @property
    def mwebmv(self):
        if self._side_properties["mwebmv"] is None:
            self._update_mwebmv_()
        return self._side_properties["mwebmv"]
    
    def set_mwebmv(self,value,force_it=False):
        if force_it is False:
            raise AssertionError("You should not manually change mwebmv, "+\
                                 "it is bound to the object coordinate."+\
                    "Set force_it to True if you really know what you are doing.")
        self._side_properties["mwebmv"] = value
    
    @property
    def _sfd98_dir(self):
        """Director where the maps are. Default option set"""
        if self._side_properties["sfd98_dir"] is None:
            from .utils.io import get_default_sfd98_dir
            self._side_properties["sfd98_dir"] = get_default_sfd98_dir(download_if_needed=True)
            
        return self._side_properties["sfd98_dir"]

    def set_sfd98_dir(self, value):
        self._side_properties['sfd98_dir'] = value
        self._update_mwebmv_()

    # ------------------ #
    # - SN Type        - #
    # ------------------ #
    @property
    def type(self):
        return self._side_properties["type"]

    @type.setter
    def type(self,value):
        self._side_properties["type"] = value

        
    # ========================= #
    # = Derived Values        = #
    # ========================= #
    # ------------------ #
    # - Distances      - #
    # ------------------ #
    # If you wnat to change these, change
    # the cosmology or the redshift   
    @property
    def distmeter(self):
        """ converts the distmpc in meter """
        return self._distance.to("m").value

    @property
    def distmeter_err(self):
        """ converts the distmpc_err in meter """
        from astropy    import units
        return self.distmpc_err * units.Mpc.in_units("m")
    
    @property
    def distmpc(self):
        """ distance to the target in megaparsec """
        return self._distance.to("Mpc").value

    @property
    def distmpc_err(self):
        """ assymetric error to account for the non linearity of
        the distance measurement
        Returns:
        --------
        [-lower, +upper] (both in absolute value).
        (np.NaN,np.NaN are returned if no error on the redshift
        """
        if self.zcmb is None or self.zcmb_err is None:
            return np.NaN,np.NaN
        
        return np.asarray([
          self.cosmo.luminosity_distance(self.zcmb-self.zcmb_err).value-self.distmpc,
          self.distmpc-self.cosmo.luminosity_distance(self.zcmb+self.zcmb_err).value])
                
    @property
    def _distance(self):
        """ astropy based distance """
        if self.cosmo is None or self.zcmb is None:
            raise AttributeError("'cosmo' and 'redshift' required.")
        return self.cosmo.luminosity_distance(self.zcmb)

    @property
    def arcsec_per_kpc(self):
        if self.cosmo is None or self.zcmb is None:
            raise AttributeError("'cosmo' and 'redshift' requiered.")
        return self.cosmo.arcsec_per_kpc_proper(self.zcmb).value
    
    # -------------------- #
    # - COSMOLOGY        - #
    # -------------------- #
    @property
    def cosmo(self):
        return self._side_properties["cosmology"]
    
    def set_cosmo(self,astropycosmo):
        """change the object according to the given cosmology.
        this cosmology must be an astropy one"""
        if "astropy" not in astropycosmo.__module__:
            raise ValueError("'astropycosmo' must be an astropy cosmology object")
        
        self._side_properties["cosmology"] = astropycosmo
        self._update_distance_()

    # ========================= #
    # = Internal Tools        = #
    # ========================= #
    def _update_(self):
        """
        Does not call _update_mewbmv_ because extinction should 
        only be fetched on demand.
        """
        self._update_distance_()
        self._check_literature_name_()
        
    def _update_mwebmv_(self):
        try:
            from sncosmo import SFD98Map
            m = SFD98Map(mapdir=self._sfd98_dir)
            self.set_mwebmv(m.get_ebv((self.ra, self.dec)), force_it=True)
        except:
            warnings.warn("MW E(B-V) could not be fetched. Please set sfd98_dir to the map directory (and check that sncosmo is installed).")
            
                
            
    def _check_literature_name_(self,verbose=False):
        if verbose:print("_check_literature_name_ to be done")
        if self.name is None:
            self._side_properties["literature_name"] = None
        
    def _update_distance_(self):
        """Do nothing so far since the distance are derived on the fly"""
        pass
          


###########################
#                         #
#   Handlers              #
#                         #
###########################
class WCSHandler( BaseObject ):
    """  Virtual Class to deal with wcs solutions """
    PROPERTIES         = []
    SIDE_PROPERTIES    = ["wcs"]
    DERIVED_PROPERTIES = []

    # ================ #
    #   Main Class     #
    # ================ #
    # ---------- #
    # CONVERTION #
    # ---------- #
    def units_to_pixels(self,units_, target=None):
        """units should be a parsable string or a astropy units. Uses the wcs method units_to_pixels()
        In addition, you have access to the 'psf [==fwhm]' unit"""
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded.")
            
        return self.wcs.units_to_pixels(units_,target=target)
    
    def pixel_to_coords(self,pixel_x,pixel_y):
        """get the coordinate (ra,dec; degree) associated to the given pixel (x,y)"""
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded.")
        #pixelsoffset = self._dataoffset
        return self.wcs.pix2world(pixel_x,pixel_y)

    def coords_to_pixel(self,ra,dec):
        """Return the pixel (x,y) associated to the given ra,dec (degree) coordinate"""
        if self.has_wcs() is False:
            raise AttributeError("no wcs solution loaded.")
        # Remark the np.asarray are only required by the astLib wcs solution
        return np.asarray(self.wcs.world2pix(ra,dec))
    
    # ---------- #
    #  SETTER    #
    # ---------- #
    def set_wcs(self, wcs, force_it=False):
        """ Attach a wcs solution to the current object """
        if self.has_wcs() and not force_it:
            raise AttributeError("A wcs solution is already loaded."\
                    " Set force_it to True if you really known what you are doing")

        from .astrometry import _MotherWCS_,get_wcs
        if wcs is not None and _MotherWCS_ not in wcs.__class__.__mro__:
            raise TypeError("The given wcs solution is not a astrobject WCS instance")
        self._side_properties["wcs"] = get_wcs(wcs) if wcs is not None else None
        
    # ================ #
    #    Properties    #
    # ================ #
    @property
    def wcs(self):
        """ wcs class """
        return self._side_properties["wcs"]
    
    def has_wcs(self):
        """ Test if the current object has a wcs solution loaded """
        return self.wcs is not None

class TargetHandler( BaseObject ):
    """ Virtual Class to deal with targets """
    
    PROPERTIES         = []
    SIDE_PROPERTIES    = ["target"]
    DERIVED_PROPERTIES = []

    def set_target(self, newtarget):
        """ set a target to the object. This target will be accessible as
        self.target. (test if a target is already set using self.has_target() )

        Parameters
        ----------
        newtarget: [AstroTarget]:
            The target object you want to attach to the object.

        Return
        ------
        Void
        """
        if newtarget is None:
            self._side_properties['target'] = None
            return
        
        # -- Input Test -- #
        if AstroTarget not in newtarget.__class__.__mro__:
            raise TypeError("'newtarget' should be (or inherite of) an AstroTarget")
        
        # -- Seems Ok -- #
        self._side_properties["target"] = newtarget.copy()

    # ================= #
    # = Property      = #
    # ================= #
    @property
    def target(self):
        """ target containing the basic object information such as:
        coordinates, redshift etc """
        return self._side_properties["target"]
    
    def has_target(self):
        """ Test if a target has been set. True means yes """
        return self.target is not None


class CatalogueHandler( BaseObject ):
    """ """
    PROPERTIES         = []
    SIDE_PROPERTIES    = ["catalogue"]
    DERIVED_PROPERTIES = ["catquery"]


    # ================= #
    #   Main            #
    # ================= #
    def download_catalogue(self,source="sdss",
                           set_it=True,force_it=False,
                           radec=None, radius=None, r_unit="degree",
                           **kwargs):
        """
        Downloads a catalogue of the given 'source'. This methods requires an
        internet connection. This downloaded catalogue is then loaded in this
        instance using the 'set_catalogue' method. If 'set_it' is set to False,
        the catalogue is returned instead of being loaded.
        
        kwargs goes to instrument.catalogue
        example: column_filters={"gmag":"13..22"}

        Return:
        ------
        Void (or the Catalogue if set_it is False)
        """
        from .instruments.instrument import fetch_catalogue
        # - RaDec
        if radec is None:
            if hasattr(self,"wcs") and self.has_wcs():
                radec = self.wcs.central_coords
            else:
                raise ValueError("Without wcs solution you need to provide 'radec'")

        
        # - Radius
        if radius is None:
            if hasattr(self,"wcs") and self.has_wcs():
                radius = self.wcs.diag_size/1.8 # not 2 to have some room around
            else:
                raise ValueError("Without wcs solution you need to provide 'radius'")
            
        elif type(radius) == str:
            radius = float(radius)
        
        # - Catalogue object found.
        query = kwargs_update(dict(source=source,radec=radec,radius=radius, r_unit=r_unit),
                               **kwargs)
        cat = fetch_catalogue(**query)
        
        # - set it of return it
        if not set_it:
            return cat

        self._derived_properties["catquery"] = query
        self.set_catalogue(cat,force_it=force_it)
        
    def set_catalogue(self, catalogue, force_it=False,
                      fast_setup=False):
        """ attach a catalogue to the current instance.
        you can then access it through 'self.catalogue'.

        The current instance's wcs solution is passed to the calague.

        Parameters
        ----------
        fast_setup: [bool] -optional-
            No Test to automatic wcs association, pure setting
        Returns
        -------
        Void
        """
        from .instruments.baseinstrument import Catalogue

        if not fast_setup:
            if self.has_catalogue() and force_it is False:
                raise AttributeError("'catalogue' already defined"+\
                        " Set force_it to True if you really known what you are doing")
    
            if Catalogue not in catalogue.__class__.__mro__:
                raise TypeError("the input 'catalogue' must be an astrobject Catalogue")
            
            if hasattr(self,"wcs") and self.has_wcs():
                catalogue.set_wcs(self.wcs, force_it=True)
                if catalogue.nobjects_in_fov < 1:
                    warnings.warn("WARNING No object in the field of view,"+"\n"+\
                                "  -> catalogue not loaded")
                    return
                
        # --------
        # - set it
        self._side_properties["catalogue"] = catalogue

    
    # ================= #
    #   Property        #
    # ================= #
    @property
    def catalogue(self):
        """ Catalogue object """
        return self._side_properties["catalogue"]
    
    def has_catalogue(self):
        """ Test if a catalogue has been assigned. True means yes """
        return self.catalogue is not None

    @property
    def _catquery(self):
        """ """
        if self._derived_properties["catquery"] is None:
            self._derived_properties["catquery"] = {}
        return self._derived_properties["catquery"]
    
