#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This modules has some basics numpy based statisitical tools """
import numpy         as np
from scipy       import stats, special
from .tools      import is_arraylike

__all__ = ["kfold_it",
           "pearson_coef","spearman_rank_coef","ks_test","aicc"]

def kfold_it(a,foldpc=10,nsample=1000):
    """This methods enables to get a kflod sample of 100-foldpc% randomly
    picked among the given data. `nsample` times realization are returned
    a could be a [NxM] array. The kfolding is performed on the 0th axis
    """
    
    a=np.asarray(a).T
    if foldpc<0 or foldpc>100:
        raise ValueError("foldpc (the percentage of removed data) must be between 0 and 100")
    size = int(len(a)*(1-foldpc/100.))
    
    indexes = np.arange(len(a))
    def get_shuffled(x):
        np.random.shuffle(x)
        return x

    return np.asarray([a[get_shuffled(indexes)[:size]].T for i in range(nsample)])

def pearson_coef(a,b,verbose=False):
    """
    This return rho,p_val,signifiance
    """
    a,b = np.asarray(a),np.asarray(b)
    if len(a) != len(a):
        raise ValueError("A and B must have the same size")

    rho,P_val = stats.pearsonr(a,b)
    significance =  rho * np.sqrt( (len(a)-2) /(1-rho**2) )
    if verbose:
        print("rho,pvalues,significance: ",rho,P_val,significance)
        
    return rho,P_val,significance

def spearman_rank_coef(a,b,axis=None,verbose=True):
    """
    based of scipy.stats.spearmanr but add the error of rho
    see http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient

    if all return are false, return_all is set True
    ----
    
    """
    a,b = np.asarray(a),np.asarray(b)
    if len(a) != len(a):
        raise ValueError("a and b must have the same size")
    
    rho,P_val = stats.spearmanr(a,b,axis=axis)
    significance =  rho*np.sqrt( (len(a)-2) /(1-rho**2) )
    return rho,P_val,significance

def get_kfolded_significance(stat,argstat,
                             foldpc=10,nsample=100):
    """This functions draws `nsample` folded samples from `argstat`
    and which should returns an array of associated statistics significance"""
    if stat=="spearman_rank_coef":
        return [spearman_rank_coef(*arg)[-1] for arg in
                kfold_it(argstat,foldpc=foldpc,nsample=nsample)]
    elif stat=="pearson_coef":
        return [pearson_coef(*arg)[-1] for arg in
                kfold_it(argstat,foldpc=foldpc,nsample=nsample)]
    else:
        raise NotImplementedError("Only spearman_rank_coef and pearson_coef stat implemented")
    
def ks_test(a,b):
    """
    Compare if 2 sample can be identical or not within a confiance value.
    The nulle hypothese is given by the p-value.
    ** This function return  *TRUE* if A and B are *DIFFERENT*
    at the `limit_test` % confident level **
    
    limit_test : if p-value < `limit_test` [in %] return True
    ----
    return Booleen
    """
    a,b = np.asarray(a),np.asarray(b)
    if len(a) != len(a):
        raise ValueError("a and b must have the same size")
    
    return stats.ks_2samp(a,b)
    

def aicc(k,L,n,logL_given=False):
    """
    ***
    Akaike_information_criterion
    ***
    
    * Only differences in AIC are meaningful *
    * AIC = \chi^2 + 2k for model comparisons. *
    
    This function has to be minimized, see
    http://en.wikipedia.org/wiki/Akaike_information_criterion
    
    where k is the number of parameters,
    L is the maximized value of the likelihood function
    n the number of element fitted (sample size)
    """
    if logL_given:
        AIC = 2*k - 2*L
    else:
        AIC = 2*k - 2*np.log(L)
        
    return AIC + (2*k*(k+1))/(n-k-1)


# ========================== #
#                            #
#   Sampler Tools            #
#                            #
# ========================== #
def gaussian_kde(*arg,**kwargs):
    """ scipy.stats' gaussian kde with 2 additional methods:
    rvs and cdf
    """
    from .decorators import make_method
    kde = stats.gaussian_kde(*arg,**kwargs)
    
    @make_method(stats.kde.gaussian_kde)
    def cdf(kde, value):
        return kde.integrate_box(-np.inf, value)
        
    @make_method(stats.kde.gaussian_kde)
    def rvs(kde, size, xrange=None, nsample=1e3):
        """ random distribution following the estimated pdf.
        random values drawn from a finite sampling of `nsample` points
        """
        # faster than resample
        if xrange is None:
            scale = np.nanmax(kde.dataset) - np.nanmin(kde.dataset)
            xrange = np.nanmin(kde.dataset) - scale*0.1, np.nanmax(kde.dataset) + scale*0.1
            
        x = np.linspace(xrange[0], xrange[1], nsample)
        return np.random.choice(x, p= kde.pdf(x) / kde.pdf(x).sum(), size=size)
            
    return kde

def continuous_poisson(mu, *args, **kwargs):
    """ Similar to a poisson distribution from scipy.stats except that it accept non integer values.
    
    Please see caveats on the input parameters descriptions
    Implemented and tested continuous variables
    - pdf
    - cdf
    - rvs
    
    Parameters
    ----------
    mu: [float]
        characteristic parameter of the poisson distribution
        Caveats: for consistency mu<0.25 is set to 0 i.e. pdf(x) = exp(-x)
        
    Returns
    -------
    array (pdf with the size of x). 
    """
    p_ = poissoncont_gen(name="poissoncont", longname='Continuous Poisson')
    return p_(mu)

class poissoncont_gen(stats._discrete_distns.poisson_gen):
    """
    Child of poisson distribution. 
    pdf, rvs and cdf have been implemented to support continuous input
    """
    def _xrange_nonzero_(self, mu):
            """ extreme x values are fixed to 0:
            - x<=0 
            - x value much greater than mu [x> max(30,mu+5*sqrt(mu))] 
            (for computing efficiency reasons.)
            """
            return np.max([0,mu-np.sqrt(mu)*6]), np.max([50,mu+np.sqrt(mu)*6])
        
    def pdf(self, x, mu, **kwargs):
            """ Hand made continuous version of the poisson distribution 
        
            Parameters
            ----------
            x: [array]
                value where the pdf will be estimated pdf(x)
                Caveats: extreme x values are fixed to 0:
                - x<=0 or below mu-5*sqrt(mu)
                - x value much greater than mu [x> max(50,mu+6*sqrt(mu))] 
                (for computing efficiency reasons.)

            Returns
            -------
            array (size of x)
            """
            xrange = self._xrange_nonzero_(mu)
            
            if is_arraylike(x):
                x = np.asarray(x) 
                Poisson_continuous = np.zeros(len(x))    
                flag_measure_it    = (x>=xrange[0]) & (x<xrange[1])
                Poisson_continuous[flag_measure_it] =  np.exp(-x[flag_measure_it]) if mu<=0.25 else \
                  mu**x[flag_measure_it] * np.exp(-mu) / special.gamma(1+x[flag_measure_it])
                return Poisson_continuous      

            return 0 if (x>=xrange[0]) or (x<xrange[1]) else np.exp(-x) if mu<=0.25 else \
                  mu**x * np.exp(-mu) / special.gamma(1+x)
            
    def rvs(self, mu,size=None, nsample=1e3, **kwargs):
            """ random distribution following the estimated pdf.
            random values drawn from a finite sampling of `nsample` points
            
            Caveats: For computing efficiency reasons the random points are drawn
                     from central values i.e. these are impossible:
                     - x<=0 or below mu-5*sqrt(mu)
                     - x value much greater than mu [x> max(50,mu+6*sqrt(mu))] 
                 
            Returns
            -------
            array (`size` random points)
            """
            # faster than resample
            xrange = self._xrange_nonzero_(mu)
            x = np.linspace(xrange[0], xrange[1], nsample)
            return np.random.choice(x, p= self.pdf(x, mu) / self.pdf(x, mu).sum(), size=size)

    def cdf(self, x, mu, **kwargs):
        """ Cumulative distribution function for continuous poisson distribution """
        if not hasattr(self, "_cdfsample"):
            # - just get it once
            self._cdfsample = self.rvs(mu, size=1e4, nsample=1e4)
        if is_arraylike(x):
            return np.asarray([float(len(self._cdfsample[self._cdfsample<x_]))/ 1e4
                        for x_ in x])
        return float(len(self._cdfsample[self._cdfsample<x]))/ 1e4
# ========================== #
#                            #
#   Outlier Rejection        #
#                            #
# ========================== #
def grubbs_criterion(npoints, alpha=0.05, twosided=True):
    """ The value (in sigma) above which a point is an outlier
    Caution the grubbs criterion is design to reject 1 outlier maxinum
    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm
    and see http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm

    This is really close to Chauvenet's criterion with p=0.1
    
    Parameters:
    -----------
    npoints: [int]
        The size of your dataset
        
    alpha: [float] -optional-
        Rejection criterion.

    twosided: [bool] -optional-
        Set this to False if you do a one-sided test.
        
    Return:
    -------
    float (numbers of sigma)
    """
    from scipy import stats
    significance_level = alpha / (2.*npoints) if twosided else alpha / np.float(npoints)
    t = stats.t.isf(significance_level, npoints-2)
    return ((npoints-1) / np.sqrt(npoints)) * (np.sqrt(t**2 / (npoints-2 + t**2)))

def chauvenet_criterion(npoints, p=0.5):
    """ Equivalent, in sigma, of the Chauvenet's criterion.
    
    (https://en.wikipedia.org/wiki/Chauvenet%27s_criterion)
    
    By default the Chauvenet's criterium reject points you have that have
    less than 50% chance to exist (p=0.5).
    You can change this value by setting p to an other quantity.

    This criterium is a 2-tailed rejection, multiply p by 2 to have a 1-tailed value.

    Parameters
    ----------
    npoints: [int]
        Number of points in you sample

    p: [float] -optional-
        value that defines the sigma cliping. (see above)
    
    Return
    ------
    float (number of sigma)
    """
    
    
    return np.abs(stats.norm.ppf(p/(2.*npoints), loc=0., scale=1.))

def chauvenet_rejection(mu_disp, sigma_disp, data, errors, npoints, p=0.5):
    """ flag data following the Chauvenet's criterium.
    
    (https://en.wikipedia.org/wiki/Chauvenet%27s_criterion)
    
    By default the Chauvenet's criterium reject points you have that have
    less than 50% chance to exist (p=0.5).
    You can change this value by setting p to an other quantity.

    This criterium is a 2-tailed rejection, multiply p by 2 to have a 1-tailed value.

    Parameters
    ----------
    mu_disp, simga_disp: [float, float]
        central value and dispersion of the model-gaussian that the data should follow

    data, errors: [float, float] (or array of)
        data and errors (without intrinsic dispersion) of the datapoint(s)

    npoints: [int]
        number of datapoint in your sample

    p: [float] -optional-
        value that defines the sigma cliping. (see above)
    
    Return
    ------
    bool (array of bool)
    NB: False means 'not to be rejected'
    """
    return np.abs((data - mu_disp)/np.sqrt(sigma_disp**2+errors**2)) > chauvenet_criterion(npoints,p=p)    
