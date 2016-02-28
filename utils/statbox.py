#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This modules has some basics numpy based statisitical tools """
import numpy         as np
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
    from scipy.stats import pearsonr
    a,b = np.asarray(a),np.asarray(b)
    if len(a) != len(a):
        raise ValueError("A and B must have the same size")

    rho,P_val = pearsonr(a,b)
    significance =  rho * np.sqrt( (len(a)-2) /(1-rho**2) )
    if verbose:
        print "rho,pvalues,significance: ",rho,P_val,significance
        
    return rho,P_val,significance

def spearman_rank_coef(a,b,axis=None,verbose=True):
    """
    based of scipy.stats.spearmanr but add the error of rho
    see http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient

    if all return are false, return_all is set True
    ----
    
    """
    from scipy.stats import spearmanr
    a,b = np.asarray(a),np.asarray(b)
    if len(a) != len(a):
        raise ValueError("a and b must have the same size")
    
    rho,P_val = spearmanr(a,b,axis=axis)
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
    from scipy.stats import ks_2samp
    a,b = np.asarray(a),np.asarray(b)
    if len(a) != len(a):
        raise ValueError("a and b must have the same size")
    
    return ks_2samp(a,b)
    

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
